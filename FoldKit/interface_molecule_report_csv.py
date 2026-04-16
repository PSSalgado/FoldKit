#!/usr/bin/env python3
"""
Extract interface rows from FoldKit interface analyser text output, optionally
restricting to specific structure files (PDB basenames) and/or chain IDs, and
write CSV table(s).

The analyser labels each interface with a chain pair like ``A-B``. With a chain
filter, rows are kept where either side matches a requested chain (case-sensitive).
With a structure filter, only rows from matching PDB sections are kept. You can
use **chains only**, **PDBs only**, or **both**.

Examples::

  python FoldKit/interface_molecule_report_csv.py results.txt -m A -m B -o interfaces_AB.csv
  python FoldKit/interface_molecule_report_csv.py results.txt --pdb model_01.pdb --chains A,B --output-dir ./out
  python FoldKit/interface_molecule_report_csv.py results.txt --pdbs model_01.pdb,model_02.pdb -o 'out/{}.csv'
  python FoldKit/interface_molecule_report_csv.py results.txt --chains A,B --combine-regex '^(model\\d+[^_]*)_' --output-dir ./out
  python FoldKit/interface_molecule_report_csv.py results.txt --chains A,B --combine-regex '^(model\\d+(?:a|del)?)_' --output-dir ./out
  python FoldKit/interface_molecule_report_csv.py results.txt --chains A,B --combine-regex '^(model\\d+)_' --output-dir ./out
  python FoldKit/interface_molecule_report_csv.py results.txt --chains A,B --combine-regex '^(model1a|model1del|model1)_' --output-dir ./out
  python FoldKit/interface_molecule_report_csv.py results.txt --chains A,B --combine-glob 'model1*' --combine-glob 'model2*' --output-dir ./out
  python FoldKit/interface_molecule_report_csv.py results.txt --chains A,B --combine-glob 'model1a*' --combine-glob 'model1del*' --combine-glob 'model1_*' --output-dir ./out
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

from cli_log import add_log_args, setup_log_from_args

# --- context lines (multi-structure / multi-set output) ---
_RE_SET = re.compile(r"^Set '([^']*)' \(patterns:")
_RE_SET_DQ = re.compile(r'^Set "([^"]*)" \(patterns:')
_RE_PROGRESS = re.compile(r"^\[(\d+)/(\d+)\]\s+(.+?)\s*$")
_RE_ANALYSING = re.compile(r"^(?:Analysing|Analyzing) interfaces in (.+?)\.\.\.\s*$")

# --- interface block ---
_RE_INTERFACE_START = re.compile(r"^Interface\s+(\d+):\s*(.+?)\s*$")

# --- metrics (indented lines under each interface) ---
_RE_CONTACT_COUNT = re.compile(r"^\s+Contact count \(within limits\):\s*(\d+)\s*$")
_RE_CONTACT_TYPES = re.compile(
    r"H-bonds[^:]*:\s*(\d+)\s+electrostatic[^:]*:\s*(\d+)\s+hydrophobic[^:]*:\s*(\d+)\s+van der Waals[^:]*:\s*(\d+)"
)
_RE_BURIED = re.compile(
    r"Buried surface area:\s*([\d.]+)\s*Å²\s+Contact area:\s*([\d.]+)\s*Å²"
)
_RE_SHAPE = re.compile(r"Complementarity \(shape\):\s*([\d.]+)")
_RE_CHARGE_OPP = re.compile(
    r"Complementarity \(charge\):\s*([\d.]+)\s+\(charged.charged contacts:\s*(\d+)\s+opposite,\s*(\d+)\s+same\)"
)
_RE_CHARGE_NONE = re.compile(
    r"Complementarity \(charge\):\s*([\d.]+)\s+\(no charged"
)
_RE_CHARGE_DENS = re.compile(
    r"Complementarity \(charge\) density:\s*([\d.]+)\s+opposite-sign contacts"
)
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


def _normalize_chain_id(s: str) -> str:
    return s.strip()


def structure_stem_from_basename(basename: str) -> str:
    """Basename like ``model_01.pdb`` → safe filename stem ``model_01``."""
    name = (basename or "").strip() or "structure"
    stem = os.path.splitext(name)[0]
    safe = "".join(ch if ch.isalnum() or ch in "-_." else "_" for ch in stem)
    return safe.strip("._") or "structure"


def structure_matches_report_name(basename: str, patterns: list[str]) -> bool:
    """
    True if ``basename`` (as in the report, e.g. ``model_01.pdb``) matches
    any pattern: exact basename or stem, or :func:`fnmatch.fnmatchcase` if the
    pattern contains ``*``, ``?``, or ``[``.
    """
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


def filter_by_structures(
    records: list[dict[str, Any]], patterns: list[str]
) -> list[dict[str, Any]]:
    """Keep only records whose ``structure_basename`` matches ``patterns``."""
    if not patterns:
        return list(records)
    return [
        r
        for r in records
        if structure_matches_report_name(str(r.get("structure_basename", "")), patterns)
    ]


def _parse_interface_block(lines: list[str]) -> dict[str, Any]:
    """Parse body lines (excluding the 'Interface N: pair' header) into flat fields."""
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
            out[f"{prefix}_chain"] = _normalize_chain_id(m.group(1))
            out[f"{prefix}_charged"] = int(m.group(2))
            out[f"{prefix}_polar"] = int(m.group(3))
            out[f"{prefix}_apolar"] = int(m.group(4))
            out[f"{prefix}_other"] = int(m.group(5))
            continue
        m = _RE_ACCESS.match(line)
        if m:
            acc_idx += 1
            prefix = f"accessibility_{acc_idx}"
            out[f"{prefix}_chain"] = _normalize_chain_id(m.group(1))
            out[f"{prefix}_avg_sasa_A2"] = float(m.group(2))
            out[f"{prefix}_accessible_fraction"] = float(m.group(3))
            continue
    return out


def parse_interface_analyser_text(text: str) -> list[dict[str, Any]]:
    """
    Parse full report text into a list of interface records (before molecule filter).
    Each record includes context: set_label, structure_basename, interface_number,
    chain_pair, chain1_id, chain2_id, and all parsed metrics.
    """
    lines = text.splitlines()
    records: list[dict[str, Any]] = []

    set_label: str | None = None
    structure: str | None = None

    i = 0
    while i < len(lines):
        line = lines[i]

        ms = _RE_SET.match(line)
        if ms:
            set_label = ms.group(1)
            i += 1
            continue
        ms = _RE_SET_DQ.match(line)
        if ms:
            set_label = ms.group(1)
            i += 1
            continue

        mp = _RE_PROGRESS.match(line)
        if mp:
            structure = mp.group(3).strip()
            i += 1
            continue

        ma = _RE_ANALYSING.match(line)
        if ma:
            structure = os.path.basename(ma.group(1).strip())
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
            records.append(rec)
            continue

        i += 1

    return records


def filter_by_molecules(
    records: list[dict[str, Any]], molecules: set[str]
) -> list[dict[str, Any]]:
    """Keep interfaces where chain1_id or chain2_id is in ``molecules``."""
    if not molecules:
        return [{**dict(r), "focus_chain": ""} for r in records]
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


_CSV_FIELDNAMES = [
    "set_label",
    "structure_basename",
    "interface_number",
    "chain_pair",
    "chain1_id",
    "chain2_id",
    "matched_molecules",
    "focus_chain",
    "contact_count",
    "hbonds",
    "electrostatic",
    "hydrophobic",
    "van_der_waals",
    "buried_surface_area_A2",
    "contact_area_A2",
    "complementarity_shape",
    "complementarity_charge",
    "charged_contacts_opposite",
    "charged_contacts_same",
    "charge_complementarity_density",
    "interface_rmsd_ca",
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


def write_csv(rows: list[dict[str, Any]], path: str | None) -> None:
    """Write CSV to ``path``; ``path is None`` writes to stdout."""
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


def group_by_structure(
    rows: list[dict[str, Any]],
) -> dict[str, list[dict[str, Any]]]:
    g: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for r in rows:
        key = str(r.get("structure_basename") or "")
        g[key].append(r)
    return dict(g)


def combine_group_key(stem: str, rx: Any) -> str:
    """
    Map a structure filename stem to a merge key: first capturing group of ``rx``,
    or ``stem`` if there is no match (each structure stays in its own group).
    """
    m = rx.match(stem)
    if m and m.lastindex and m.lastindex >= 1:
        return m.group(1)
    return stem


def _glob_pattern_to_group_key(pat: str) -> str:
    """
    Filename stem for merged CSV: text before the first ``*`` or ``?``, with a
    trailing ``_`` removed (so ``model1_*`` → ``model1``, not ``model1_``).
    """
    if "*" in pat:
        key = pat.split("*", 1)[0] or pat
        return key.rstrip("_") or key
    if "?" in pat:
        key = pat.split("?", 1)[0] or pat
        return key.rstrip("_") or key
    return pat


def combine_group_key_glob(stem: str, patterns: list[str]) -> str:
    """
    First ``patterns`` entry (in order) that matches ``stem`` wins; group key from
    that pattern via :func:`_glob_pattern_to_group_key`. No match → ``stem``.
    Put **more specific** patterns first (e.g. ``model1a*`` before ``model1_*``),
    because a broad ``model1*`` also matches ``model1del_*`` and ``model1a_*``.
    """
    for pat in patterns:
        if fnmatch.fnmatchcase(stem, pat):
            return _glob_pattern_to_group_key(pat)
    return stem


def _sort_merged_rows(rows: list[dict[str, Any]]) -> None:
    rows.sort(
        key=lambda r: (
            str(r.get("structure_basename", "")),
            int(r.get("interface_number") or 0),
            str(r.get("focus_chain", "")),
        )
    )


def build_output_batches(
    by_st: dict[str, list[dict[str, Any]]],
    combine_rx: Any | None,
    combine_globs: list[str] | None,
) -> list[tuple[str, list[dict[str, Any]]]]:
    """
    Return (filename_stem, rows) for each CSV. With no combiner, one file per
    structure. With ``combine_rx`` or ``combine_globs``, merge rows that map to
    the same group key.
    """
    if not combine_rx and not combine_globs:
        return [
            (structure_stem_from_basename(bn), by_st[bn])
            for bn in sorted(by_st.keys())
        ]

    merged: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for basename in sorted(by_st.keys()):
        stem = os.path.splitext(os.path.basename(basename.strip()))[0]
        if combine_rx:
            gkey = combine_group_key(stem, combine_rx)
        else:
            assert combine_globs is not None
            gkey = combine_group_key_glob(stem, combine_globs)
        merged[gkey].extend(by_st[basename])

    for rows in merged.values():
        _sort_merged_rows(rows)

    out: list[tuple[str, list[dict[str, Any]]]] = []
    for gkey in sorted(merged.keys()):
        stem_safe = structure_stem_from_basename(f"{gkey}.pdb")
        out.append((stem_safe, merged[gkey]))
    return out


def _collect_chains_ordered(
    args: argparse.Namespace,
) -> tuple[list[str], set[str]]:
    """Preserve first-seen order for --group-by-chain; set is used for filtering."""
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
    """
    With ``enabled`` and two or more chains: one CSV row per (interface, focus_chain)
    for each requested chain that participates, sorted by chain order (A rows, then B, …).
    An A–B interface yields two rows (focus_chain A, focus_chain B). Otherwise one row
    per interface; ``focus_chain`` is set for a single requested chain, else empty.
    """
    if not rows:
        return rows
    if not enabled or len(ordered_chains) <= 1:
        out: list[dict[str, Any]] = []
        for r in rows:
            row = dict(r)
            if len(ordered_chains) == 1:
                row["focus_chain"] = ordered_chains[0]
            else:
                row["focus_chain"] = ""
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
            "Filter FoldKit interface analyser text output by optional structure (PDB) "
            "name and/or chain ID(s), and write CSV of matching interfaces."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Give at least one of: structure selection (--pdb / --pdbs) and/or chain selection
(-m / --molecule and/or --chains). Chain IDs are case-sensitive. Structure patterns match the basename
as printed in the report (e.g. model_01.pdb or stem model_01); use * or ? for
fnmatch-style globs.

Examples:
  %(prog)s results.txt -m A -o interfaces_chain_A.csv
  %(prog)s results.txt --pdb model_01.pdb --chains A,B --output-dir ./out
  %(prog)s results.txt --pdbs model_01.pdb,model_02.pdb -m A -m B -o 'out/{}.csv'
  %(prog)s results.txt --pdb '*.pdb' --output-dir ./out
  %(prog)s results.txt --chains A,B --group-by-chain --output-dir ./out
  %(prog)s results.txt --chains A,B --combine-regex '^(model\\d+[^_]*)_' --output-dir ./out
  %(prog)s results.txt --chains A,B --combine-regex '^(model\\d+(?:a|del)?)_' --output-dir ./out
  %(prog)s results.txt --chains A,B --combine-regex '^(model\\d+)_' --output-dir ./out
  %(prog)s results.txt --chains A,B --combine-regex '^(model1a|model1del|model1)_' --output-dir ./out
  %(prog)s results.txt --chains A,B --combine-glob 'model1*' --combine-glob 'model2*' --output-dir ./out
  %(prog)s results.txt --chains A,B --combine-glob 'model1a*' --combine-glob 'model1del*' --combine-glob 'model1_*' --output-dir ./out
""",
    )
    ap.add_argument(
        "report",
        help="Text report path (typically results.txt from interface_analyser_asu_charge.py or interface_analyser_asu_ec.py, or lattice variants); use - for stdin.",
    )
    ap.add_argument(
        "--pdb",
        "--structure",
        action="append",
        default=[],
        metavar="PATTERN",
        help=(
            "Only include interfaces from this structure: report basename or stem "
            "(e.g. model_01.pdb or model_01). Repeat for multiple structures. "
            "Patterns with wildcards (*, ?, [) are matched as globs."
        ),
    )
    ap.add_argument(
        "--pdbs",
        dest="pdbs_csv",
        metavar="LIST",
        help="Comma-separated structure names/patterns (same as --pdb).",
    )
    ap.add_argument(
        "-m",
        "--molecule",
        action="append",
        default=[],
        metavar="CHAIN",
        help="Chain ID: keep interfaces involving this chain (repeat for more chains).",
    )
    ap.add_argument(
        "--chains",
        metavar="LIST",
        help="Comma-separated chain IDs (alternative to repeating -m).",
    )
    ap.add_argument(
        "-o",
        "--output",
        metavar="PATH",
        help=(
            "Single CSV path, '-' for stdout, or a template with '{}' for the structure "
            "stem (e.g. 'out/{}.csv'). One file per PDB when the report has multiple "
            "structures: use this pattern or --output-dir."
        ),
    )
    ap.add_argument(
        "-d",
        "--output-dir",
        metavar="DIR",
        help=(
            "Write one CSV per structure: DIR/<pdb_stem>.csv (pdb_stem from the report "
            "basename, e.g. model_01.pdb → model_01.csv). For multi-PDB reports."
        ),
    )
    ap.add_argument(
        "--group-by-chain",
        action="store_true",
        help=(
            "With two or more chains: one row per (interface, focus_chain) so each PDB "
            "file lists all interfaces involving chain A, then all involving B (column "
            "focus_chain). A–B interfaces appear twice. Sort order follows the chain list "
            "you gave (e.g. --chains A,B)."
        ),
    )
    ap.add_argument(
        "--combine-regex",
        metavar="REGEX",
        help=(
            "Merge structures whose stem matches the same first capturing group (output "
            "CSV stem). Generic suffix before the next underscore: '^(model\\\\d+[^_]*)_' "
            "(model1_2m → model1, model1a_2m → model1a); narrower: '^(model\\\\d+)_' or "
            "'^(model\\\\d+(?:a|del)?)_'. Shell-style: --combine-glob. Requires "
            "--output-dir or -o '.../{}.csv' unless one merge group only."
        ),
    )
    ap.add_argument(
        "--combine-glob",
        action="append",
        default=[],
        dest="combine_globs",
        metavar="PATTERN",
        help=(
            "Like --combine-regex but uses fnmatch on the stem. Output stem is the "
            "text before the first * (trailing _ removed, so model1_* → model1.csv). "
            "Repeat for multiple groups; list more specific patterns first "
            "(e.g. model1a* then model1del* then model1_*). Mutually exclusive with "
            "--combine-regex."
        ),
    )
    add_log_args(ap)
    args = ap.parse_args()
    setup_log_from_args(args, script_path=__file__, inputs=[getattr(args, "report", "")], pattern=None)

    pdb_patterns = _collect_pdb_patterns(args)
    ordered_chains, chains = _collect_chains_ordered(args)
    if not pdb_patterns and not chains:
        ap.error(
            "Specify at least one structure (--pdb / --pdbs) and/or chain "
            "(-m / --molecule / --chains)."
        )

    if args.report == "-":
        text = sys.stdin.read()
    else:
        if not os.path.isfile(args.report):
            ap.error(f"File not found: {args.report}")
        with open(args.report, encoding="utf-8", errors="replace") as f:
            text = f.read()

    if args.output_dir and args.output:
        ap.error("Use either --output-dir or -o, not both.")

    globs = [x.strip() for x in (args.combine_globs or []) if x.strip()]

    combine_rx: Any | None = None
    if args.combine_regex:
        if globs:
            ap.error("Use either --combine-regex or --combine-glob, not both.")
        try:
            combine_rx = re.compile(args.combine_regex)
        except re.error as e:
            ap.error(f"Invalid --combine-regex: {e}")
        if combine_rx.groups < 1:
            ap.error(
                "--combine-regex must include at least one capturing group, e.g. "
                r"'^(model\d+[^_]*)_', '^(model\d+(?:a|del)?)_', or '^(model\d+)_'."
            )

    all_recs = parse_interface_analyser_text(text)
    filtered = filter_by_structures(all_recs, pdb_patterns)
    filtered = filter_by_molecules(filtered, chains)
    filtered = apply_group_by_chain_layout(filtered, ordered_chains, args.group_by_chain)
    by_st = group_by_structure(filtered)
    n_struct = len(by_st)
    n_rows = len(filtered)

    if n_rows == 0:
        print(
            f"No matching interfaces (parsed {len(all_recs)} interface block(s)).",
            file=sys.stderr,
        )
        return

    batches = build_output_batches(by_st, combine_rx, globs if globs else None)

    out_arg = args.output
    if (combine_rx or globs) and len(batches) > 1:
        if not args.output_dir and not (
            out_arg and out_arg != "-" and "{}" in out_arg
        ):
            ap.error(
                "Multiple merge groups; use --output-dir DIR or -o 'path/{}.csv' "
                "({} = group stem), or omit -o for stdout."
            )

    # --- Per-structure or combined-group files ---
    if args.output_dir:
        os.makedirs(args.output_dir, exist_ok=True)
        written: list[tuple[str, int]] = []
        for file_stem, rows in batches:
            path = os.path.join(args.output_dir, f"{file_stem}.csv")
            write_csv(rows, path)
            written.append((path, len(rows)))
        parts = ", ".join(f"{p} ({n} row(s))" for p, n in written)
        tag = (
            " (merged by --combine-regex)"
            if combine_rx
            else (" (merged by --combine-glob)" if globs else "")
        )
        print(
            f"Wrote {len(written)} file(s){tag}: {parts} "
            f"[{n_rows} row(s) from {len(all_recs)} interface block(s) parsed].",
            file=sys.stderr,
        )
        return

    if out_arg and out_arg != "-" and "{}" in out_arg:
        written = []
        for file_stem, rows in batches:
            path = out_arg.replace("{}", file_stem)
            parent = os.path.dirname(path)
            if parent:
                os.makedirs(parent, exist_ok=True)
            write_csv(rows, path)
            written.append((path, len(rows)))
        parts = ", ".join(f"{p} ({n} row(s))" for p, n in written)
        tag = (
            " (merged by --combine-regex)"
            if combine_rx
            else (" (merged by --combine-glob)" if globs else "")
        )
        print(
            f"Wrote {len(written)} file(s){tag}: {parts} "
            f"[{n_rows} row(s) from {len(all_recs)} interface block(s) parsed].",
            file=sys.stderr,
        )
        return

    # --- Single destination (stdout or one path) ---
    out_path: str | None = None if not out_arg or out_arg == "-" else out_arg
    if len(batches) > 1 and out_path and "{}" not in out_arg:
        ap.error(
            "Multiple output groups; use --output-dir DIR or -o 'path/{}.csv' "
            "(with '{}'), or omit -o to merge all rows to stdout."
        )

    rows_out = batches[0][1] if batches else filtered
    write_csv(rows_out, out_path)
    dest = "stdout" if out_path is None else out_path
    print(
        f"Wrote {len(rows_out)} row(s) to {dest} "
        f"({n_struct} structure(s); {len(all_recs)} interface block(s) parsed).",
        file=sys.stderr,
    )


if __name__ == "__main__":
    main()
