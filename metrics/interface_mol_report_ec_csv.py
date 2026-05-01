#!/usr/bin/env python3
"""
Extract interface rows and/or summary totals from FoldKit *EC* interface analyser text output.

Inputs:
  - metrics/interface_analyser_asu_ec.py output (ASU, EC)
  - metrics/interface_analyser_lattice_ec.py output (lattice, EC)

Modes:
  - details (default): per-interface rows (EC metrics) + optional summary rows
  - summary: summary totals only (one row per structure per report)
"""

from __future__ import annotations

import argparse
import csv
import fnmatch
import os
import re
import sys
from collections import defaultdict
from typing import Any


# --- context lines (multi-structure / multi-set output) ---
_RE_SET = re.compile(r"^Set '([^']*)' \(patterns:")
_RE_SET_DQ = re.compile(r'^Set "([^"]*)" \(patterns:')
_RE_PROGRESS = re.compile(r"^\[(\d+)/(\d+)\]\s+(.+?)\s*$")
_RE_ANALYSING = re.compile(r"^(?:Analysing|Analyzing) interfaces in (.+?)\.\.\.\s*$")

# --- summary header ---
_RE_RESULTS_HEADER = re.compile(r"^INTERFACE ANALYSIS RESULTS:\s*$")

# --- summary scalar lines (common) ---
_RE_SUM_TOTAL_INTERFACES = re.compile(r"^Total interfaces:\s*(\d+)\s*$")
_RE_SUM_TOTAL_BSA = re.compile(r"^Total buried surface area:\s*([\d.]+)\s*Å²\s*$")
_RE_SUM_AVG_BSA = re.compile(r"^Average buried area per interface:\s*([\d.]+)\s*Å²\s*$")

# When a reference chain is used, newer reports can include an explicit all-chains
# summary section with indented lines. These are optional; parsers should ignore
# them if absent.
_RE_SUM_ALL_HEADER = re.compile(r"^All-chains interface summary\b")
_RE_SUM_REF_HEADER = re.compile(r"^Reference-chain interface summary\b")
_RE_SUM_ALL_TOTAL_INTERFACES = re.compile(r"^\s{2}Total interfaces:\s*(\d+)\s*$")
_RE_SUM_ALL_TOTAL_BSA = re.compile(r"^\s{2}Total buried surface area:\s*([\d.]+)\s*Å²\s*$")
_RE_SUM_ALL_AVG_BSA = re.compile(r"^\s{2}Average buried area per interface:\s*([\d.]+)\s*Å²\s*$")
_RE_SUM_ALL_ISO_CHAIN = re.compile(r"^\s{4}Chain\s+(.+?):\s*([\d.]+)\s*Å²\s*$")
_RE_SUM_ALL_ISO_SUM = re.compile(r"^\s{4}Sum isolated SASA:\s*([\d.]+)\s*Å²\s*$")

# --- isolated SASA block ---
_RE_SUM_ISO_CHAIN = re.compile(r"^\s+Chain\s+(.+?):\s*([\d.]+)\s*Å²\s*$")
_RE_SUM_ISO_SUM = re.compile(
    r"^\s+Sum of isolated SASA .*:\s*([\d.]+)\s*Å²\s*$"
)

# --- lattice reference chain block (may appear in EC lattice reports) ---
_RE_SUM_LAT_REF = re.compile(r"^Lattice reference chain\s+(.+?)\s+\(multi-copy / cluster SASA\):\s*$")
_RE_SUM_LAT_SASA_ISO = re.compile(r"^\s+SASA isolated \(chain alone\):\s*([\d.]+)\s*Å²\s*$")
_RE_SUM_LAT_SASA_CLUSTER = re.compile(
    r"^\s+SASA in full model .*:\s*([\d.]+)\s*Å²\s*$"
)
_RE_SUM_LAT_REF_BSA = re.compile(
    r"^\s+Reference-chain BSA \(SASA_iso.*\):\s*([\d.]+)\s*Å²\s*$"
)
_RE_SUM_LAT_NORM_REF_DIV = re.compile(
    r"^\s+Normalisation divisors - reference chain: residues=(\d+) atoms=(\d+) mass=([\d.]+) Da \(([\d.]+) kDa\)\s*$"
)
_RE_SUM_LAT_ISO_PR_REF = re.compile(
    r"^\s+SASA isolated per residue \(/ reference chain\):\s*([\d.eE+-]+)\s*Å²\s*$"
)
_RE_SUM_LAT_ISO_PA_REF = re.compile(
    r"^\s+SASA isolated per atom \(/ reference chain\):\s*([\d.eE+-]+)\s*Å²\s*$"
)
_RE_SUM_LAT_ISO_PK_REF = re.compile(
    r"^\s+SASA isolated per kDa protein \(/ reference chain\):\s*([\d.eE+-]+)\s*Å²/kDa\s*$"
)
_RE_SUM_LAT_CL_PR_REF = re.compile(
    r"^\s+SASA in cluster per residue \(/ reference chain\):\s*([\d.eE+-]+)\s*Å²\s*$"
)
_RE_SUM_LAT_CL_PA_REF = re.compile(
    r"^\s+SASA in cluster per atom \(/ reference chain\):\s*([\d.eE+-]+)\s*Å²\s*$"
)
_RE_SUM_LAT_CL_PK_REF = re.compile(
    r"^\s+SASA in cluster per kDa protein \(/ reference chain\):\s*([\d.eE+-]+)\s*Å²/kDa\s*$"
)
_RE_SUM_LAT_REF_BSA_PR_REF = re.compile(
    r"^\s+Reference-chain BSA per residue \(/ reference chain\):\s*([\d.eE+-]+)\s*Å²\s*$"
)
_RE_SUM_LAT_REF_BSA_PA_REF = re.compile(
    r"^\s+Reference-chain BSA per atom \(/ reference chain\):\s*([\d.eE+-]+)\s*Å²\s*$"
)
_RE_SUM_LAT_REF_BSA_PK_REF = re.compile(
    r"^\s+Reference-chain BSA per kDa protein \(/ reference chain\):\s*([\d.eE+-]+)\s*Å²/kDa\s*$"
)
_RE_SUM_LAT_BURIAL_FRAC = re.compile(
    r"^\s+Lattice burial fraction .*:\s*([\d.]+)\s*\(\s*([\d.]+)\s*%\s*\)\s*$"
)
_RE_SUM_LAT_CONTACT_RES_FRAC = re.compile(
    r"^\s+Fraction of residues with cross-chain neighbour within .*:\s*([\d.]+)\s*\(\s*([\d.]+)\s*%\s*\)\s*$"
)

# --- lattice EC summary lines ---
_RE_SUM_LAT_EC_R_BSA = re.compile(r"^\s+Lattice EC \(r, BSA-weighted Fisher-z\):\s*([\d.]+)\s*$")
_RE_SUM_LAT_EC_R_NP = re.compile(
    r"^\s+Lattice EC \(r, n_pairs-weighted Fisher-z\):\s*([\d.]+)\s+\(total_n_pairs=(\d+)\)\s*$"
)
_RE_SUM_LAT_EC_D_BSA = re.compile(
    r"^\s+Lattice EC density \(BSA-weighted\):\s*([\d.eE+-]+)\s+r/Å²\s+\(per\s+(.+?)\)\s*$"
)
_RE_SUM_LAT_EC_D_NP = re.compile(
    r"^\s+Lattice EC density \(n_pairs-weighted\):\s*([\d.eE+-]+)\s+r/Å²\s+\(per\s+(.+?)\)\s*$"
)
_RE_SUM_LAT_EC_PARTNER = re.compile(
    r"^\s+Partner\s+(.+?):\s+r=([-\d.]+)\s+\(n_pairs=(\d+)\)\s*$"
)

# --- interface block ---
_RE_INTERFACE_START = re.compile(r"^Interface\s+(\d+):\s*(.+?)\s*$")

# --- metrics under each interface (EC analyser interface blocks) ---
_RE_CONTACT_COUNT = re.compile(r"^\s+Contact count \(within limits\):\s*(\d+)\s*$")
_RE_CONTACT_TYPES = re.compile(
    r"H-bonds[^:]*:\s*(\d+)\s+electrostatic[^:]*:\s*(\d+)\s+hydrophobic[^:]*:\s*(\d+)\s+van der Waals[^:]*:\s*(\d+)"
)
_RE_BURIED = re.compile(
    r"Buried surface area:\s*([\d.]+)\s*Å²\s+Contact area:\s*([\d.]+)\s*Å²"
)
_RE_SHAPE = re.compile(r"Complementarity \(shape\):\s*([\d.]+)")
_RE_EC_R = re.compile(r"EC \(r\):\s*([-\d.]+)\s+\(n_pairs=(\d+)\)")
_RE_EC_DENS = re.compile(r"EC density:\s*([-\d.eE+]+)\s+r/Å²\s+\(per\s+(.+?)\)")
_RE_RMSD = re.compile(r"Interface RMSD \(CA\):\s*([\d.]+)\s*Å")
_RE_RMSD_NA = re.compile(r"Interface RMSD \(CA\):\s*N/A")
_RE_DIST = re.compile(
    r"Distance \(Å\): min=([\d.]+)\s+max=([\d.]+)\s+avg=([\d.]+)"
)
_RE_POLARITY = re.compile(
    r"^\s+Polarity\s+(.+?):\s+charged=(\d+)\s+polar=(\d+)\s+apolar=(\d+)\s+other=(\d+)\s*$"
)
_RE_ACCESS = re.compile(
    r"^\s+Accessibility\s+(.+?):\s+avg SASA=([\d.]+)\s*Å²\s+accessible_fraction=([\d.]+)\s*$"
)

_RE_BLOCK_BREAK = re.compile(
    r"^(={10,}|\[\d+/\d+\]|Analysing interfaces in |Analyzing interfaces in |INTERFACE ANALYSIS RESULTS:|Set ['\"]|Done\.|Error:)"
)


def _split_chain_pair(pair: str) -> tuple[str, str] | None:
    pair = pair.strip()
    if "-" not in pair:
        return None
    a, b = pair.split("-", 1)
    return a.strip(), b.strip()


def _normalise_chain_id(s: str) -> str:
    return s.strip()


def structure_matches_report_name(basename: str, patterns: list[str]) -> bool:
    if not patterns:
        return True
    b = (basename or "").strip()
    stem = os.path.splitext(b)[0]
    for raw in patterns:
        pat = raw.strip()
        if not pat:
            continue
        if pat in ("*", "all"):
            return True
        if any(c in pat for c in "*?["):
            if fnmatch.fnmatchcase(b, pat) or fnmatch.fnmatchcase(stem, pat):
                return True
            continue
        if b == pat or stem == pat:
            return True
    return False


def _parse_interface_block(lines: list[str]) -> dict[str, Any]:
    out: dict[str, Any] = {}
    pol_idx = 0
    acc_idx = 0
    for line in lines:
        m = _RE_CONTACT_COUNT.match(line)
        if m:
            out["contact_count"] = int(m.group(1))
            continue
        m = _RE_CONTACT_TYPES.search(line)
        if m:
            out["hbonds"] = int(m.group(1))
            out["electrostatic"] = int(m.group(2))
            out["hydrophobic"] = int(m.group(3))
            out["van_der_waals"] = int(m.group(4))
            continue
        m = _RE_BURIED.search(line)
        if m:
            out["buried_surface_area_A2"] = float(m.group(1))
            out["contact_area_A2"] = float(m.group(2))
            continue
        m = _RE_SHAPE.search(line)
        if m:
            out["complementarity_shape"] = float(m.group(1))
            continue
        m = _RE_EC_R.search(line)
        if m:
            out["ec_r"] = float(m.group(1))
            out["ec_n_pairs"] = int(m.group(2))
            continue
        m = _RE_EC_DENS.search(line)
        if m:
            out["ec_density"] = float(m.group(1))
            out["ec_density_denominator"] = str(m.group(2)).strip()
            out["ec_density_percent"] = float(m.group(1)) * 100.0
            continue
        if _RE_RMSD.search(line):
            m = _RE_RMSD.search(line)
            assert m
            out["interface_rmsd_ca"] = float(m.group(1))
            continue
        if _RE_RMSD_NA.search(line):
            out["interface_rmsd_ca"] = ""
            continue
        m = _RE_DIST.search(line)
        if m:
            out["distance_min_A"] = float(m.group(1))
            out["distance_max_A"] = float(m.group(2))
            out["distance_avg_A"] = float(m.group(3))
            continue
        m = _RE_POLARITY.match(line)
        if m:
            pol_idx += 1
            prefix = f"polarity_{pol_idx}"
            out[f"{prefix}_chain"] = _normalise_chain_id(m.group(1))
            out[f"{prefix}_charged"] = int(m.group(2))
            out[f"{prefix}_polar"] = int(m.group(3))
            out[f"{prefix}_apolar"] = int(m.group(4))
            out[f"{prefix}_other"] = int(m.group(5))
            continue
        m = _RE_ACCESS.match(line)
        if m:
            acc_idx += 1
            prefix = f"accessibility_{acc_idx}"
            out[f"{prefix}_chain"] = _normalise_chain_id(m.group(1))
            out[f"{prefix}_avg_sasa_A2"] = float(m.group(2))
            out[f"{prefix}_accessible_fraction"] = float(m.group(3))
            continue
    return out


def parse_ec_report_text(text: str) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    lines = text.splitlines()
    ifaces: list[dict[str, Any]] = []
    summaries: list[dict[str, Any]] = []

    set_label: str | None = None
    structure: str | None = None

    in_summary = False
    current_sum: dict[str, Any] | None = None
    iso_by_chain: dict[str, float] = {}
    iso_by_chain_all: dict[str, float] = {}
    ec_by_partner: dict[str, dict[str, Any]] = {}
    sum_section = "ref"

    def flush_summary() -> None:
        nonlocal in_summary, current_sum, iso_by_chain, iso_by_chain_all, ec_by_partner, sum_section
        if not in_summary or current_sum is None:
            return
        if ec_by_partner:
            current_sum["lattice_ec_by_partner_chain"] = ";".join(
                f"{k}:r={ec_by_partner[k]['ec_r']:.3f},n_pairs={int(ec_by_partner[k]['ec_n_pairs'])}"
                for k in sorted(ec_by_partner.keys())
            )
        else:
            current_sum["lattice_ec_by_partner_chain"] = ""
        summaries.append(current_sum)
        in_summary = False
        current_sum = None
        iso_by_chain = {}
        iso_by_chain_all = {}
        ec_by_partner = {}
        sum_section = "ref"

    i = 0
    while i < len(lines):
        line = lines[i]

        ms = _RE_SET.match(line) or _RE_SET_DQ.match(line)
        if ms:
            flush_summary()
            set_label = ms.group(1)
            i += 1
            continue

        mp = _RE_PROGRESS.match(line)
        if mp:
            flush_summary()
            structure = mp.group(3).strip()
            i += 1
            continue

        ma = _RE_ANALYSING.match(line)
        if ma:
            flush_summary()
            structure = os.path.basename(ma.group(1).strip())
            i += 1
            continue

        if _RE_RESULTS_HEADER.match(line):
            in_summary = True
            current_sum = {
                "set_label": set_label or "",
                "structure_basename": structure or "",
            }
            iso_by_chain = {}
            iso_by_chain_all = {}
            ec_by_partner = {}
            sum_section = "ref"
            i += 1
            continue

        if in_summary and current_sum is not None:
            if _RE_SUM_ALL_HEADER.match(line):
                sum_section = "all"
                i += 1
                continue
            if _RE_SUM_REF_HEADER.match(line):
                sum_section = "ref"
                i += 1
                continue

            if sum_section == "all":
                m = _RE_SUM_ALL_TOTAL_INTERFACES.match(line)
                if m:
                    current_sum["all_total_interfaces"] = int(m.group(1))
                    i += 1
                    continue
                m = _RE_SUM_ALL_TOTAL_BSA.match(line)
                if m:
                    current_sum["all_total_buried_surface_area_A2"] = float(m.group(1))
                    i += 1
                    continue
                m = _RE_SUM_ALL_AVG_BSA.match(line)
                if m:
                    current_sum["all_average_buried_area_per_interface_A2"] = float(m.group(1))
                    i += 1
                    continue
                m = _RE_SUM_ALL_ISO_CHAIN.match(line)
                if m:
                    cid = _normalise_chain_id(m.group(1))
                    iso_by_chain_all[cid] = float(m.group(2))
                    i += 1
                    continue
                m = _RE_SUM_ALL_ISO_SUM.match(line)
                if m:
                    current_sum["all_sasa_isolated_sum_A2"] = float(m.group(1))
                    current_sum["all_sasa_isolated_by_chain"] = ";".join(
                        f"{k}={iso_by_chain_all[k]:.1f}" for k in sorted(iso_by_chain_all.keys())
                    ) if iso_by_chain_all else ""
                    i += 1
                    continue

            m = _RE_SUM_TOTAL_INTERFACES.match(line)
            if m:
                current_sum["total_interfaces"] = int(m.group(1))
                i += 1
                continue
            m = _RE_SUM_TOTAL_BSA.match(line)
            if m:
                current_sum["total_buried_surface_area_A2"] = float(m.group(1))
                i += 1
                continue
            m = _RE_SUM_AVG_BSA.match(line)
            if m:
                current_sum["average_buried_area_per_interface_A2"] = float(m.group(1))
                i += 1
                continue

            m = _RE_SUM_ISO_CHAIN.match(line)
            if m:
                cid = _normalise_chain_id(m.group(1))
                iso_by_chain[cid] = float(m.group(2))
                i += 1
                continue
            m = _RE_SUM_ISO_SUM.match(line)
            if m:
                current_sum["sasa_isolated_sum_A2"] = float(m.group(1))
                current_sum["sasa_isolated_by_chain"] = ";".join(
                    f"{k}={iso_by_chain[k]:.1f}" for k in sorted(iso_by_chain.keys())
                ) if iso_by_chain else ""
                i += 1
                continue

            # lattice reference info (if present)
            m = _RE_SUM_LAT_REF.match(line)
            if m:
                current_sum["lattice_reference_chain"] = _normalise_chain_id(m.group(1))
                i += 1
                continue
            m = _RE_SUM_LAT_SASA_ISO.match(line)
            if m:
                current_sum["sasa_reference_isolated_A2"] = float(m.group(1))
                i += 1
                continue
            m = _RE_SUM_LAT_SASA_CLUSTER.match(line)
            if m:
                current_sum["sasa_reference_in_cluster_A2"] = float(m.group(1))
                i += 1
                continue
            m = _RE_SUM_LAT_REF_BSA.match(line)
            if m:
                current_sum["reference_chain_BSA_A2"] = float(m.group(1))
                i += 1
                continue
            m = _RE_SUM_LAT_NORM_REF_DIV.match(line)
            if m:
                current_sum["reference_chain_residue_count"] = int(m.group(1))
                current_sum["reference_chain_atom_count"] = int(m.group(2))
                current_sum["reference_chain_mass_Da"] = float(m.group(3))
                current_sum["reference_chain_mass_kDa"] = float(m.group(4))
                i += 1
                continue
            m = _RE_SUM_LAT_ISO_PR_REF.match(line)
            if m:
                current_sum["sasa_reference_isolated_per_residue_reference_chain_A2"] = float(m.group(1))
                i += 1
                continue
            m = _RE_SUM_LAT_ISO_PA_REF.match(line)
            if m:
                current_sum["sasa_reference_isolated_per_atom_reference_chain_A2"] = float(m.group(1))
                i += 1
                continue
            m = _RE_SUM_LAT_ISO_PK_REF.match(line)
            if m:
                current_sum["sasa_reference_isolated_per_kDa_reference_chain_A2"] = float(m.group(1))
                i += 1
                continue
            m = _RE_SUM_LAT_CL_PR_REF.match(line)
            if m:
                current_sum["sasa_reference_in_cluster_per_residue_reference_chain_A2"] = float(m.group(1))
                i += 1
                continue
            m = _RE_SUM_LAT_CL_PA_REF.match(line)
            if m:
                current_sum["sasa_reference_in_cluster_per_atom_reference_chain_A2"] = float(m.group(1))
                i += 1
                continue
            m = _RE_SUM_LAT_CL_PK_REF.match(line)
            if m:
                current_sum["sasa_reference_in_cluster_per_kDa_reference_chain_A2"] = float(m.group(1))
                i += 1
                continue
            m = _RE_SUM_LAT_REF_BSA_PR_REF.match(line)
            if m:
                current_sum["reference_chain_BSA_per_residue_reference_chain_A2"] = float(m.group(1))
                i += 1
                continue
            m = _RE_SUM_LAT_REF_BSA_PA_REF.match(line)
            if m:
                current_sum["reference_chain_BSA_per_atom_reference_chain_A2"] = float(m.group(1))
                i += 1
                continue
            m = _RE_SUM_LAT_REF_BSA_PK_REF.match(line)
            if m:
                current_sum["reference_chain_BSA_per_kDa_reference_chain_A2"] = float(m.group(1))
                i += 1
                continue
            m = _RE_SUM_LAT_BURIAL_FRAC.match(line)
            if m:
                current_sum["lattice_burial_fraction"] = float(m.group(1))
                current_sum["lattice_burial_fraction_percent"] = float(m.group(2))
                i += 1
                continue
            m = _RE_SUM_LAT_CONTACT_RES_FRAC.match(line)
            if m:
                current_sum["lattice_contact_residue_fraction"] = float(m.group(1))
                current_sum["lattice_contact_residue_fraction_percent"] = float(m.group(2))
                i += 1
                continue

            m = _RE_SUM_LAT_EC_R_BSA.match(line)
            if m:
                current_sum["lattice_ec_r_bsa_weighted"] = float(m.group(1))
                i += 1
                continue
            m = _RE_SUM_LAT_EC_R_NP.match(line)
            if m:
                current_sum["lattice_ec_r_npairs_weighted"] = float(m.group(1))
                current_sum["lattice_ec_total_npairs"] = int(m.group(2))
                i += 1
                continue
            m = _RE_SUM_LAT_EC_D_BSA.match(line)
            if m:
                current_sum["lattice_ec_density_bsa_weighted"] = float(m.group(1))
                current_sum["lattice_ec_density_denominator"] = str(m.group(2)).strip()
                current_sum["lattice_ec_density_bsa_weighted_percent"] = float(m.group(1)) * 100.0
                i += 1
                continue
            m = _RE_SUM_LAT_EC_D_NP.match(line)
            if m:
                current_sum["lattice_ec_density_npairs_weighted"] = float(m.group(1))
                current_sum["lattice_ec_density_denominator"] = str(m.group(2)).strip()
                current_sum["lattice_ec_density_npairs_weighted_percent"] = float(m.group(1)) * 100.0
                i += 1
                continue
            m = _RE_SUM_LAT_EC_PARTNER.match(line)
            if m:
                cid = _normalise_chain_id(m.group(1))
                ec_by_partner[cid] = {
                    "partner_chain": cid,
                    "ec_r": float(m.group(2)),
                    "ec_n_pairs": int(m.group(3)),
                }
                i += 1
                continue

        mi = _RE_INTERFACE_START.match(line)
        if mi:
            iface_num = int(mi.group(1))
            pair_raw = mi.group(2).strip()
            i += 1
            body: list[str] = []
            while i < len(lines):
                nxt = lines[i]
                if _RE_INTERFACE_START.match(nxt):
                    break
                if _RE_BLOCK_BREAK.match(nxt):
                    break
                body.append(nxt)
                i += 1

            pair = _split_chain_pair(pair_raw)
            c1, c2 = (pair if pair else ("", ""))
            metrics = _parse_interface_block(body)
            rec: dict[str, Any] = {
                "set_label": set_label or "",
                "structure_basename": structure or "",
                "interface_number": iface_num,
                "chain_pair": pair_raw,
                "chain1_id": c1,
                "chain2_id": c2,
            }
            rec.update(metrics)
            ifaces.append(rec)
            continue

        i += 1

    flush_summary()

    return ifaces, summaries


def filter_by_structures(records: list[dict[str, Any]], patterns: list[str]) -> list[dict[str, Any]]:
    if not patterns:
        return list(records)
    return [
        r
        for r in records
        if structure_matches_report_name(str(r.get("structure_basename", "")), patterns)
    ]


def filter_by_molecules(records: list[dict[str, Any]], molecules: set[str]) -> list[dict[str, Any]]:
    if not molecules:
        return [{**dict(r), "matched_molecules": "", "focus_chain": ""} for r in records]
    out: list[dict[str, Any]] = []
    for r in records:
        c1 = _normalise_chain_id(str(r.get("chain1_id", "")))
        c2 = _normalise_chain_id(str(r.get("chain2_id", "")))
        matched = sorted({x for x in (c1, c2) if x in molecules})
        if not matched:
            continue
        row = dict(r)
        row["matched_molecules"] = ";".join(matched)
        row["focus_chain"] = ""
        out.append(row)
    return out


def _collect_chains_ordered(args: argparse.Namespace) -> tuple[list[str], set[str]]:
    ordered: list[str] = []
    seen: set[str] = set()

    def add(m: str) -> None:
        m = _normalise_chain_id(m)
        if not m or m in seen:
            return
        seen.add(m)
        ordered.append(m)

    for m in args.molecule:
        add(m)
    if args.chains:
        for part in args.chains.split(","):
            add(part)
    return ordered, seen


def apply_group_by_chain_layout(
    rows: list[dict[str, Any]],
    ordered_chains: list[str],
    enabled: bool,
) -> list[dict[str, Any]]:
    if not rows:
        return rows
    if not enabled or len(ordered_chains) <= 1:
        out: list[dict[str, Any]] = []
        for r in rows:
            row = dict(r)
            row["focus_chain"] = ordered_chains[0] if len(ordered_chains) == 1 else ""
            out.append(row)
        return out

    order_idx = {c: i for i, c in enumerate(ordered_chains)}
    out = []
    for r in rows:
        c1 = _normalise_chain_id(str(r.get("chain1_id", "")))
        c2 = _normalise_chain_id(str(r.get("chain2_id", "")))
        in_iface = {c1, c2}
        for foc in ordered_chains:
            if foc not in in_iface:
                continue
            row = dict(r)
            row["focus_chain"] = foc
            out.append(row)
    out.sort(
        key=lambda x: (
            order_idx.get(str(x.get("focus_chain")), 99),
            str(x.get("structure_basename", "")),
            int(x.get("interface_number") or 0),
        )
    )
    return out


def group_by_structure(rows: list[dict[str, Any]]) -> dict[str, list[dict[str, Any]]]:
    g: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for r in rows:
        g[str(r.get("structure_basename") or "")].append(r)
    return dict(g)


def _sort_rows(rows: list[dict[str, Any]]) -> None:
    rows.sort(
        key=lambda r: (
            str(r.get("structure_basename", "")),
            0 if r.get("interface_number") in ("", None) else 1,
            int(r.get("interface_number") or 0),
            str(r.get("focus_chain", "")),
        )
    )


_CSV_SUMMARY_FIELDS = [
    "report_txt",
    "structure_basename",
    "all_total_interfaces",
    "all_total_buried_surface_area_A2",
    "all_average_buried_area_per_interface_A2",
    "all_sasa_isolated_sum_A2",
    "all_sasa_isolated_by_chain",
    "total_interfaces",
    "total_buried_surface_area_A2",
    "average_buried_area_per_interface_A2",
    "sasa_isolated_sum_A2",
    "sasa_isolated_by_chain",
    "lattice_reference_chain",
    "sasa_reference_isolated_A2",
    "sasa_reference_in_cluster_A2",
    "reference_chain_BSA_A2",
    "reference_chain_BSA_per_residue_reference_chain_A2",
    "reference_chain_BSA_per_atom_reference_chain_A2",
    "reference_chain_BSA_per_kDa_reference_chain_A2",
    "reference_chain_residue_count",
    "reference_chain_atom_count",
    "reference_chain_mass_Da",
    "reference_chain_mass_kDa",
    "sasa_reference_isolated_per_residue_reference_chain_A2",
    "sasa_reference_isolated_per_atom_reference_chain_A2",
    "sasa_reference_isolated_per_kDa_reference_chain_A2",
    "sasa_reference_in_cluster_per_residue_reference_chain_A2",
    "sasa_reference_in_cluster_per_atom_reference_chain_A2",
    "sasa_reference_in_cluster_per_kDa_reference_chain_A2",
    "lattice_burial_fraction",
    "lattice_burial_fraction_percent",
    "lattice_contact_residue_fraction",
    "lattice_contact_residue_fraction_percent",
    "lattice_ec_r_bsa_weighted",
    "lattice_ec_r_npairs_weighted",
    "lattice_ec_total_npairs",
    "lattice_ec_density_bsa_weighted",
    "lattice_ec_density_npairs_weighted",
    "lattice_ec_density_bsa_weighted_percent",
    "lattice_ec_density_npairs_weighted_percent",
    "lattice_ec_by_partner_chain",
]

_CSV_DETAILS_FIELDS = [
    # Summary-mode columns first (same order as --mode summary)
    *_CSV_SUMMARY_FIELDS,
    # Then per-interface columns
    "interface_number",
    "chain_pair",
    "chain1_id",
    "chain2_id",
    "contact_count",
    "hbonds",
    "electrostatic",
    "hydrophobic",
    "van_der_waals",
    "buried_surface_area_A2",
    "contact_area_A2",
    "complementarity_shape",
    "ec_r",
    "ec_n_pairs",
    "ec_density",
    "ec_density_percent",
    "distance_min_A",
    "distance_max_A",
    "distance_avg_A",
    "polarity_1_chain",
    "polarity_1_charged",
    "polarity_1_polar",
    "polarity_1_apolar",
    "polarity_1_other",
    "polarity_2_chain",
    "polarity_2_charged",
    "polarity_2_polar",
    "polarity_2_apolar",
    "polarity_2_other",
    "accessibility_1_chain",
    "accessibility_1_avg_sasa_A2",
    "accessibility_1_accessible_fraction",
    "accessibility_2_chain",
    "accessibility_2_avg_sasa_A2",
    "accessibility_2_accessible_fraction",
]


def write_csv(rows: list[dict[str, Any]], path: str | None, fieldnames: list[str]) -> None:
    use_stdout = path is None
    f = sys.stdout if use_stdout else open(path, "w", newline="", encoding="utf-8")
    try:
        w = csv.DictWriter(f, fieldnames=fieldnames, extrasaction="ignore")
        w.writeheader()
        for r in rows:
            w.writerow({k: r.get(k, "") for k in fieldnames})
    finally:
        if not use_stdout:
            f.close()


def _collect_pdb_patterns(args: argparse.Namespace) -> list[str]:
    pats: list[str] = []
    for p in args.pdb:
        p = p.strip()
        if p:
            pats.append(p)
    if args.pdbs_csv:
        for part in args.pdbs_csv.split(","):
            part = part.strip()
            if part:
                pats.append(part)
    return pats


def main() -> None:
    ap = argparse.ArgumentParser(
        description=(
            "Extract interfaces and/or totals from FoldKit *EC* interface analyser text output, "
            "filter by structure and/or chain ID(s), and write CSV."
        )
    )
    ap.add_argument(
        "report",
        nargs="+",
        help="One or more analyser text report path(s). Use - to read one report from stdin.",
    )
    ap.add_argument("--pdb", action="append", default=[], metavar="PATTERN", help="Only include this structure (basename/stem/glob). Repeatable.")
    ap.add_argument("--pdbs", dest="pdbs_csv", metavar="LIST", help="Comma-separated structure names/patterns (same as --pdb).")
    ap.add_argument("-m", "--molecule", action="append", default=[], metavar="CHAIN", help="Chain ID to keep (repeatable).")
    ap.add_argument("--chains", metavar="LIST", help="Comma-separated chain IDs (alternative to repeating -m).")
    ap.add_argument("--group-by-chain", action="store_true", help="With 2+ chains: duplicate A–B rows per focus_chain to group output by chain.")
    ap.add_argument(
        "--mode",
        choices=["details", "summary"],
        default="details",
        help="details = per-interface rows (+ optional summary rows); summary = summary totals only.",
    )
    ap.add_argument(
        "--include-summary",
        action="store_true",
        default=True,
        help="When --mode details, include parsed per-structure totals as extra CSV row(s). (default: on)",
    )
    ap.add_argument(
        "--no-summary",
        action="store_false",
        dest="include_summary",
        help="When --mode details, do not include totals rows.",
    )
    ap.add_argument("-o", "--output", metavar="PATH", help="CSV path, '-' for stdout. If omitted, writes to stdout.")
    args = ap.parse_args()

    pdb_patterns = _collect_pdb_patterns(args)
    ordered_chains, chains = _collect_chains_ordered(args)
    if not pdb_patterns and not chains:
        ap.error("Specify at least one structure (--pdb/--pdbs) and/or chain (-m/--chains).")

    reports: list[tuple[str, str]] = []
    for rp in args.report:
        rp = str(rp)
        if rp == "-":
            reports.append(("<stdin>", sys.stdin.read()))
            continue
        if not os.path.isfile(rp):
            ap.error(f"File not found: {rp}")
        with open(rp, encoding="utf-8", errors="replace") as f:
            reports.append((rp, f.read()))

    all_rows: list[dict[str, Any]] = []
    for report_path, text in reports:
        iface_rows, sum_rows = parse_ec_report_text(text)
        sum_rows = filter_by_structures(sum_rows, pdb_patterns)
        if args.mode == "summary":
            for s in sum_rows:
                row = dict(s)
                row["report_txt"] = os.path.basename(str(report_path))
                all_rows.append(row)
            continue

        iface_rows = filter_by_structures(iface_rows, pdb_patterns)
        iface_rows = filter_by_molecules(iface_rows, chains)
        iface_rows = apply_group_by_chain_layout(iface_rows, ordered_chains, args.group_by_chain)
        # In details mode, include summary as extra row(s) (one per structure) when enabled.
        if args.include_summary:
            for s in sum_rows:
                row = dict(s)
                row["report_txt"] = os.path.basename(str(report_path))
                all_rows.append(row)

        for r in iface_rows:
            row = dict(r)
            row["report_txt"] = os.path.basename(str(report_path))
            all_rows.append(row)

    if not all_rows:
        print("No matching rows.", file=sys.stderr)
        return

    _sort_rows(all_rows)
    out_path = None if not args.output or args.output == "-" else args.output
    fieldnames = _CSV_SUMMARY_FIELDS if args.mode == "summary" else _CSV_DETAILS_FIELDS
    write_csv(all_rows, out_path, fieldnames=fieldnames)


if __name__ == "__main__":
    main()

