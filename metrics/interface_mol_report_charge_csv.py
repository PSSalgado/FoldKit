#!/usr/bin/env python3
"""
Extract interface rows + summary rows from FoldKit *charge* interface analyser text output.

Inputs:
  - metrics/interface_analyser_asu_charge.py output (ASU, charge)
  - metrics/interface_analyser_lattice_charge.py output (lattice, charge)

This script understands both ASU and lattice report text and can optionally
include summary totals as extra CSV row(s) per structure.
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

# --- isolated SASA block ---
_RE_SUM_ISO_CHAIN = re.compile(r"^\s+Chain\s+(.+?):\s*([\d.]+)\s*Å²\s*$")
_RE_SUM_ISO_SUM = re.compile(
    r"^\s+Sum of isolated SASA .*:\s*([\d.]+)\s*Å²\s*$"
)

# When `--reference-chain` is used, newer lattice reports include an explicit
# all-chains summary section with indented lines (optional).
_RE_SUM_ALL_HEADER = re.compile(r"^All-chains interface summary\b")
_RE_SUM_REF_HEADER = re.compile(r"^Reference-chain interface summary\b")
_RE_SUM_ALL_TOTAL_INTERFACES = re.compile(r"^\s{2}Total interfaces:\s*(\d+)\s*$")
_RE_SUM_ALL_TOTAL_BSA = re.compile(r"^\s{2}Total buried surface area:\s*([\d.]+)\s*Å²\s*$")
_RE_SUM_ALL_AVG_BSA = re.compile(r"^\s{2}Average buried area per interface:\s*([\d.]+)\s*Å²\s*$")
_RE_SUM_ALL_ISO_CHAIN = re.compile(r"^\s{4}Chain\s+(.+?):\s*([\d.]+)\s*Å²\s*$")
_RE_SUM_ALL_ISO_SUM = re.compile(r"^\s{4}Sum isolated SASA:\s*([\d.]+)\s*Å²\s*$")

# --- lattice reference chain block ---
_RE_SUM_LAT_REF = re.compile(r"^Lattice reference chain\s+(.+?)\s+\(multi-copy / cluster SASA\):\s*$")
_RE_SUM_LAT_SASA_ISO = re.compile(r"^\s+SASA isolated \(chain alone\):\s*([\d.]+)\s*Å²\s*$")
_RE_SUM_LAT_SASA_CLUSTER = re.compile(
    r"^\s+SASA in full model .*:\s*([\d.]+)\s*Å²\s*$"
)
_RE_SUM_LAT_REF_BSA = re.compile(
    r"^\s+Reference-chain BSA \(SASA_iso.*\):\s*([\d.]+)\s*Å²\s*$"
)
_RE_SUM_LAT_NORM_REF_DIV = re.compile(
    r"^\s+Normalization divisors - reference chain: residues=(\d+) atoms=(\d+) mass=([\d.]+) Da \(([\d.]+) kDa\)\s*$"
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

# lattice charge complementarity (atom-contact) + density
_RE_SUM_LAT_CC = re.compile(
    r"^\s+Lattice charge complementarity:\s*([\d.]+)\s*\(\s*([\d.]+)\s*%\s*\)\s*(?:\((?:.*)\))?\s*$"
)
_RE_SUM_LAT_CC_COUNTS = re.compile(
    r"^\s+Lattice charge complementarity:\s*([\d.]+).*charged–charged contacts:\s*(\d+)\s+opposite,\s*(\d+)\s+same\)\s*$"
)
_RE_SUM_LAT_CCD = re.compile(
    r"^\s+Lattice charge complementarity density:\s*([\d.]+)\s+opposite-sign contacts/Å²\s+\(per\s+(.+?)\)\s*$"
)

# lattice extended metrics (residue-pair)
_RE_SUM_LAT_RP = re.compile(
    r"^\s+Lattice charge complementarity \(residue-pair\):\s*([\d.]+)\s*\(\s*([\d.]+)\s*%\s*\)\s*\(pairs:\s*(\d+)\s+opposite,\s*(\d+)\s+same\)\s*$"
)
_RE_SUM_LAT_RPW = re.compile(
    r"^\s+Lattice charge complementarity \(residue-pair, 1/d²-weighted\):\s*([\d.]+)\s*\(\s*([\d.]+)\s*%\s*\)\s*$"
)


# --- interface block ---
_RE_INTERFACE_START = re.compile(r"^Interface\s+(\d+):\s*(.+?)\s*$")

# --- metrics under each interface (charge-specific) ---
_RE_CONTACT_COUNT = re.compile(r"^\s+Contact count \(within limits\):\s*(\d+)\s*$")
_RE_BURIED = re.compile(
    r"Buried surface area:\s*([\d.]+)\s*Å²\s+Contact area:\s*([\d.]+)\s*Å²"
)
_RE_CHARGE_OPP = re.compile(
    r"Complementarity \(charge\):\s*([\d.]+)\s+\(charged.charged contacts:\s*(\d+)\s+opposite,\s*(\d+)\s+same\)"
)
_RE_CHARGE_NONE = re.compile(
    r"Complementarity \(charge\):\s*([\d.]+)\s+\(no charged"
)
_RE_CHARGE_DENS = re.compile(
    r"Complementarity \(charge\) density:\s*([\d.]+)\s+opposite-sign contacts"
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


def _normalize_chain_id(s: str) -> str:
    return s.strip()


def structure_stem_from_basename(basename: str) -> str:
    name = (basename or "").strip() or "structure"
    stem = os.path.splitext(name)[0]
    safe = "".join(ch if ch.isalnum() or ch in "-_." else "_" for ch in stem)
    return safe.strip("._") or "structure"


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
    for line in lines:
        m = _RE_CONTACT_COUNT.match(line)
        if m:
            out["contact_count"] = int(m.group(1))
            continue
        m = _RE_BURIED.search(line)
        if m:
            out["buried_surface_area_A2"] = float(m.group(1))
            out["contact_area_A2"] = float(m.group(2))
            continue
        m = _RE_CHARGE_OPP.search(line)
        if m:
            out["complementarity_charge"] = float(m.group(1))
            out["charged_contacts_opposite"] = int(m.group(2))
            out["charged_contacts_same"] = int(m.group(3))
            continue
        m = _RE_CHARGE_NONE.search(line)
        if m:
            out["complementarity_charge"] = float(m.group(1))
            out["charged_contacts_opposite"] = ""
            out["charged_contacts_same"] = ""
            continue
        m = _RE_CHARGE_DENS.search(line)
        if m:
            out["charge_complementarity_density"] = float(m.group(1))
            continue
    return out


def parse_charge_report_text(text: str) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    """
    Return (interface_records, summary_records).

    - interface_records: one row per Interface N block.
    - summary_records: one row per structure (parsed from the summary section).
    """
    lines = text.splitlines()
    ifaces: list[dict[str, Any]] = []
    summaries: list[dict[str, Any]] = []

    set_label: str | None = None
    structure: str | None = None

    # summary state
    in_summary = False
    current_sum: dict[str, Any] | None = None
    iso_by_chain: dict[str, float] = {}
    iso_by_chain_all: dict[str, float] = {}
    lat_ref_chain: str | None = None
    sum_section = "ref"

    def flush_summary() -> None:
        nonlocal in_summary, current_sum, iso_by_chain, iso_by_chain_all, lat_ref_chain, sum_section
        if not in_summary or current_sum is None:
            return
        summaries.append(current_sum)
        in_summary = False
        current_sum = None
        iso_by_chain = {}
        iso_by_chain_all = {}
        lat_ref_chain = None
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
            lat_ref_chain = None
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
                    cid = _normalize_chain_id(m.group(1))
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
                cid = _normalize_chain_id(m.group(1))
                iso_by_chain[cid] = float(m.group(2))
                i += 1
                continue
            m = _RE_SUM_ISO_SUM.match(line)
            if m:
                current_sum["sasa_isolated_sum_A2"] = float(m.group(1))
                # store packed dict as "A=123.4;B=..."
                if iso_by_chain:
                    current_sum["sasa_isolated_by_chain"] = ";".join(
                        f"{k}={iso_by_chain[k]:.1f}" for k in sorted(iso_by_chain.keys())
                    )
                else:
                    current_sum["sasa_isolated_by_chain"] = ""
                i += 1
                continue

            m = _RE_SUM_LAT_REF.match(line)
            if m:
                lat_ref_chain = _normalize_chain_id(m.group(1))
                current_sum["lattice_reference_chain"] = lat_ref_chain
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
                i += 1
                continue
            m = _RE_SUM_LAT_CONTACT_RES_FRAC.match(line)
            if m:
                current_sum["lattice_contact_residue_fraction"] = float(m.group(1))
                i += 1
                continue

            m = _RE_SUM_LAT_CC_COUNTS.match(line)
            if m:
                current_sum["lattice_charge_complementarity"] = float(m.group(1))
                current_sum["lattice_charge_contacts_opposite"] = int(m.group(2))
                current_sum["lattice_charge_contacts_same"] = int(m.group(3))
                i += 1
                continue
            m = _RE_SUM_LAT_CC.match(line)
            if m:
                current_sum["lattice_charge_complementarity"] = float(m.group(1))
                i += 1
                continue
            m = _RE_SUM_LAT_CCD.match(line)
            if m:
                current_sum["lattice_charge_complementarity_density"] = float(m.group(1))
                current_sum["lattice_charge_complementarity_density_denominator"] = str(m.group(2)).strip()
                i += 1
                continue

            m = _RE_SUM_LAT_RP.match(line)
            if m:
                current_sum["lattice_charge_residue_pair_score"] = float(m.group(1))
                current_sum["lattice_charge_residue_pair_opposite"] = int(m.group(3))
                current_sum["lattice_charge_residue_pair_same"] = int(m.group(4))
                i += 1
                continue
            m = _RE_SUM_LAT_RPW.match(line)
            if m:
                current_sum["lattice_charge_residue_pair_weighted_score"] = float(m.group(1))
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

        # end-of-structure / next header: if we were in summary and are leaving it, store current summary
        if in_summary and current_sum is not None and (
            _RE_INTERFACE_START.match(line)
            or _RE_PROGRESS.match(line)
            or _RE_SET.match(line)
            or _RE_SET_DQ.match(line)
        ):
            summaries.append(current_sum)
            in_summary = False
            current_sum = None
            iso_by_chain = {}
            lat_ref_chain = None

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
        c1 = _normalize_chain_id(str(r.get("chain1_id", "")))
        c2 = _normalize_chain_id(str(r.get("chain2_id", "")))
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
        m = _normalize_chain_id(m)
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
        c1 = _normalize_chain_id(str(r.get("chain1_id", "")))
        c2 = _normalize_chain_id(str(r.get("chain2_id", "")))
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
            0 if str(r.get("row_type", "")) == "summary" else 1,
            int(r.get("interface_number") or 0),
            str(r.get("focus_chain", "")),
        )
    )


_CSV_FIELDNAMES = [
    "report_path",
    "row_type",  # interface | summary
    "set_label",
    "structure_basename",
    "interface_number",
    "chain_pair",
    "chain1_id",
    "chain2_id",
    "matched_molecules",
    "focus_chain",
    # interface metrics (charge-only)
    "contact_count",
    "buried_surface_area_A2",
    "contact_area_A2",
    "complementarity_charge",
    "charged_contacts_opposite",
    "charged_contacts_same",
    "charge_complementarity_density",
    # summary metrics (structure-level totals; lattice fields empty for ASU)
    # When `--reference-chain` is set: optional all-chains totals before reference-scoped totals.
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
    "lattice_contact_residue_fraction",
    "lattice_charge_complementarity",
    "lattice_charge_contacts_opposite",
    "lattice_charge_contacts_same",
    "lattice_charge_complementarity_density",
    "lattice_charge_complementarity_density_denominator",
    "lattice_charge_residue_pair_score",
    "lattice_charge_residue_pair_opposite",
    "lattice_charge_residue_pair_same",
    "lattice_charge_residue_pair_weighted_score",
]


def write_csv(rows: list[dict[str, Any]], path: str | None) -> None:
    use_stdout = path is None
    f = sys.stdout if use_stdout else open(path, "w", newline="", encoding="utf-8")
    try:
        w = csv.DictWriter(f, fieldnames=_CSV_FIELDNAMES, extrasaction="ignore")
        w.writeheader()
        for r in rows:
            w.writerow({k: r.get(k, "") for k in _CSV_FIELDNAMES})
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
            "Extract interfaces and totals from FoldKit *charge* interface analyser text output, "
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
        "--include-summary",
        action="store_true",
        default=True,
        help="Include parsed per-structure totals as extra CSV row(s). (default: on)",
    )
    ap.add_argument(
        "--no-summary",
        action="store_false",
        dest="include_summary",
        help="Do not include totals rows; output interfaces only.",
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
        iface_rows, sum_rows = parse_charge_report_text(text)
        iface_rows = filter_by_structures(iface_rows, pdb_patterns)
        iface_rows = filter_by_molecules(iface_rows, chains)
        iface_rows = apply_group_by_chain_layout(iface_rows, ordered_chains, args.group_by_chain)
        for r in iface_rows:
            row = dict(r)
            row["report_path"] = report_path
            row["row_type"] = "interface"
            all_rows.append(row)

        if args.include_summary:
            sum_rows = filter_by_structures(sum_rows, pdb_patterns)
            for s in sum_rows:
                row = dict(s)
                row["report_path"] = report_path
                row["row_type"] = "summary"
                row.setdefault("matched_molecules", "")
                row.setdefault("focus_chain", "")
                all_rows.append(row)

    if not all_rows:
        print("No matching rows.", file=sys.stderr)
        return

    _sort_rows(all_rows)
    out_path = None if not args.output or args.output == "-" else args.output
    write_csv(all_rows, out_path)


if __name__ == "__main__":
    main()

