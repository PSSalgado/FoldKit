#!/usr/bin/env python3
"""
Extract atom–atom contact rows from ``contact_analyser.py`` text output (and
optional ``*_asu_contacts.txt`` sidecar files), filter by structure basename
and/or chain ID, and write CSV table(s).

Columns: set_label, structure_basename, chain1, chain2, res1, atom1, res2,
atom2, distance_A, contact_type.

Examples::

  python metrics/contact_molecule_report_csv.py contact_results.txt -m A -m B -o contacts_AB.csv
  python metrics/contact_molecule_report_csv.py contact_results.txt --pdb model_01.pdb --chains A,B --output-dir ./out
  python metrics/contact_molecule_report_csv.py contact_results_model_01_asu_contacts.txt --structure-basename model_01.pdb -m A -o model_01_contacts.csv
"""

from __future__ import annotations

import argparse
import csv
import os
import re
import sys
from typing import Any

_REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if _REPO_ROOT not in sys.path:
    sys.path.insert(0, _REPO_ROOT)

from utils.cli_log import add_log_args, setup_log_from_args

# Sibling script (shared structure matching and CSV batching)
_SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
if _SCRIPT_DIR not in sys.path:
    sys.path.insert(0, _SCRIPT_DIR)

from interface_molecule_report_csv import (  # noqa: E402
    build_output_batches,
    filter_by_structures,
    group_by_structure,
)

_RE_SET = re.compile(r"^Set '([^']*)' \(patterns:")
_RE_SET_DQ = re.compile(r'^Set "([^"]*)" \(patterns:')
_RE_PROGRESS = re.compile(r"^\[(\d+)/(\d+)\]\s+(.+?)\s*$")
_RE_ANALYSING = re.compile(r"^(?:Analysing|Analyzing) contacts in (.+?)\.\.\.\s*$")

_RE_CONTACT_HEADER = re.compile(
    r"^Contact details \(chain1 chain2 res1 atom1 res2 atom2 distance type\)", re.I
)
_RE_ASU_SIDECAR_HDR = re.compile(r"^#\s*All\s+\d+\s+ASU\s+contacts", re.I)

_RE_BLOCK_BREAK = re.compile(
    r"^(={10,}|\[\d+/\d+\]|Analysing contacts in |Analyzing contacts in |CONTACT ANALYSIS RESULTS:|"
    r"Set ['\"]|Done\.|Error:|===\s*$)"
)


def _parse_contact_data_line(line: str) -> dict[str, Any] | None:
    """
    One contact line: chain1 chain2 res1 atom1 res2 atom2 <dist> [Å] <type>.
    Inline report uses two decimals and a literal Å before type; sidecar omits Å.
    """
    s = line.strip()
    if not s or s.startswith("#"):
        return None
    if _RE_CONTACT_HEADER.match(s):
        return None
    if s.startswith("Full list") or s.startswith("Contact details"):
        return None

    tok = s.split()
    if len(tok) < 8:
        return None
    if tok[-2] in ("Å", "A"):
        try:
            dist = float(tok[-3])
        except ValueError:
            return None
        ctype = tok[-1]
        core = tok[:-3]
    else:
        try:
            dist = float(tok[-2])
        except ValueError:
            return None
        ctype = tok[-1]
        core = tok[:-2]
    if len(core) != 6:
        return None
    return {
        "chain1": core[0],
        "chain2": core[1],
        "res1": core[2],
        "atom1": core[3],
        "res2": core[4],
        "atom2": core[5],
        "distance_A": dist,
        "contact_type": ctype,
    }


def parse_contact_analyser_text(
    text: str, *, default_structure: str = ""
) -> list[dict[str, Any]]:
    """
    Parse contact_analyser output into contact records with set_label and
    structure_basename context.
    """
    lines = text.splitlines()
    records: list[dict[str, Any]] = []

    set_label: str | None = None
    structure: str | None = None
    in_details = False

    def flush_context_break(line: str) -> bool:
        nonlocal in_details
        if not in_details:
            return False
        if _RE_BLOCK_BREAK.match(line.strip()):
            in_details = False
            return True
        return False

    i = 0
    while i < len(lines):
        line = lines[i]

        if flush_context_break(line):
            i += 1
            continue

        ms = _RE_SET.match(line)
        if ms:
            set_label = ms.group(1)
            in_details = False
            i += 1
            continue
        ms = _RE_SET_DQ.match(line)
        if ms:
            set_label = ms.group(1)
            in_details = False
            i += 1
            continue

        mp = _RE_PROGRESS.match(line)
        if mp:
            structure = os.path.basename(mp.group(3).strip())
            in_details = False
            i += 1
            continue

        ma = _RE_ANALYSING.match(line)
        if ma:
            structure = os.path.basename(ma.group(1).strip())
            in_details = False
            i += 1
            continue

        if _RE_CONTACT_HEADER.match(line.strip()):
            in_details = True
            i += 1
            continue

        if _RE_ASU_SIDECAR_HDR.match(line.strip()):
            in_details = True
            i += 1
            continue

        if in_details:
            parsed = _parse_contact_data_line(line)
            if parsed:
                rec = {
                    "set_label": set_label or "",
                    "structure_basename": structure or "",
                }
                rec.update(parsed)
                records.append(rec)
            i += 1
            continue

        i += 1

    fb = (default_structure or "").strip()
    if fb:
        for r in records:
            if not str(r.get("structure_basename", "")).strip():
                r["structure_basename"] = fb

    if not records and fb:
        for line in lines:
            parsed = _parse_contact_data_line(line)
            if parsed:
                records.append(
                    {
                        "set_label": "",
                        "structure_basename": fb,
                        **parsed,
                    }
                )

    return records


def filter_by_molecules(
    records: list[dict[str, Any]], molecules: set[str]
) -> list[dict[str, Any]]:
    if not molecules:
        return list(records)
    out = []
    for r in records:
        c1 = str(r.get("chain1", "")).strip()
        c2 = str(r.get("chain2", "")).strip()
        if c1 in molecules or c2 in molecules:
            out.append(r)
    return out


_CSV_FIELDNAMES = [
    "set_label",
    "structure_basename",
    "chain1",
    "chain2",
    "res1",
    "atom1",
    "res2",
    "atom2",
    "distance_A",
    "contact_type",
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


def _collect_chains_and_patterns(args: argparse.Namespace) -> tuple[set[str], list[str]]:
    chains: set[str] = set()
    for m in args.molecule:
        m = m.strip()
        if m:
            chains.add(m)
    if args.chains:
        for part in args.chains.split(","):
            part = part.strip()
            if part:
                chains.add(part)

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
    return chains, pats


def main() -> None:
    ap = argparse.ArgumentParser(
        description=(
            "Filter contact_analyser.py text output by optional structure name "
            "and/or chain ID(s), and write CSV of matching atom–atom contacts."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Chain IDs are case-sensitive. Structure patterns match the basename as printed
in the report (e.g. model_01.pdb); use * or ? for fnmatch-style globs.

Use --structure-basename when the input is only an *_asu_contacts.txt sidecar
(with no progress lines) so rows are tagged with that PDB basename for filtering.
""",
    )
    ap.add_argument(
        "report",
        help="Text report path (contact_analyser.py -o); use - for stdin.",
    )
    ap.add_argument(
        "--pdb",
        "--structure",
        action="append",
        default=[],
        metavar="PATTERN",
        help="Only contacts from this structure basename/stem (repeat or use --pdbs).",
    )
    ap.add_argument(
        "--pdbs",
        dest="pdbs_csv",
        metavar="LIST",
        help="Comma-separated structure names/patterns.",
    )
    ap.add_argument(
        "-m",
        "--molecule",
        action="append",
        default=[],
        metavar="CHAIN",
        help="Chain ID: keep contacts touching this chain (repeat for more).",
    )
    ap.add_argument(
        "--chains",
        metavar="LIST",
        help="Comma-separated chain IDs (alternative to repeating -m).",
    )
    ap.add_argument(
        "--structure-basename",
        metavar="NAME",
        help=(
            "Assign this basename (e.g. model_01.pdb) to rows with no structure "
            "context (sidecar-only files)."
        ),
    )
    ap.add_argument(
        "-o",
        "--output",
        metavar="PATH",
        help="CSV path, '-' for stdout, or template 'out/{}.csv' per structure stem.",
    )
    ap.add_argument(
        "-d",
        "--output-dir",
        metavar="DIR",
        help="One CSV per structure: DIR/<stem>.csv.",
    )
    ap.add_argument(
        "--combine-regex",
        metavar="REGEX",
        help=(
            "Merge structures whose stem shares the same first capturing group "
            "(same grouping rule as the interface report CSV extractors)."
        ),
    )
    ap.add_argument(
        "--combine-glob",
        action="append",
        default=[],
        dest="combine_globs",
        metavar="PATTERN",
        help="Shell-style combine groups (repeat); mutually exclusive with --combine-regex.",
    )
    add_log_args(ap)
    args = ap.parse_args()
    setup_log_from_args(
        args,
        script_path=__file__,
        inputs=[getattr(args, "report", "")],
        pattern=None,
    )

    chains, pdb_patterns = _collect_chains_and_patterns(args)
    if not pdb_patterns and not chains:
        ap.error(
            "Specify at least one of: structure (--pdb / --pdbs) and/or chain "
            "(-m / --chains)."
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
            ap.error("--combine-regex must include at least one capturing group.")

    default_st = (args.structure_basename or "").strip()
    if not default_st and args.report != "-":
        base = os.path.basename(args.report)
        if "_asu_contacts.txt" in base:
            # e.g. contact_results_model_01_asu_contacts.txt → stem model_01.pdb
            if base.endswith("_asu_contacts.txt"):
                stem = base[: -len("_asu_contacts.txt")]
                if "_" in stem:
                    stem = stem.rsplit("_", 1)[-1]
                default_st = f"{stem}.pdb"

    all_recs = parse_contact_analyser_text(text, default_structure=default_st)
    filtered = filter_by_structures(all_recs, pdb_patterns)
    filtered = filter_by_molecules(filtered, chains)
    by_st = group_by_structure(filtered)
    n_rows = len(filtered)

    if n_rows == 0:
        print(
            f"No matching contacts (parsed {len(all_recs)} contact line(s) total).",
            file=sys.stderr,
        )
        return

    batches = build_output_batches(by_st, combine_rx, globs if globs else None)

    out_arg = args.output
    if (combine_rx or globs) and len(batches) > 1:
        if not args.output_dir and not (out_arg and out_arg != "-" and "{}" in out_arg):
            ap.error(
                "Multiple merge groups; use --output-dir DIR or -o 'path/{}.csv', "
                "or omit -o for stdout."
            )

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
            f"[{n_rows} row(s) from {len(all_recs)} contact line(s) parsed].",
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
            f"[{n_rows} row(s) from {len(all_recs)} contact line(s) parsed].",
            file=sys.stderr,
        )
        return

    out_path: str | None = None if not out_arg or out_arg == "-" else out_arg
    if len(batches) > 1 and out_path and "{}" not in (out_arg or ""):
        ap.error(
            "Multiple output groups; use --output-dir DIR or -o 'path/{}.csv', "
            "or omit -o to merge all rows to stdout."
        )

    rows_out = batches[0][1] if batches else filtered
    write_csv(rows_out, out_path)
    dest = "stdout" if out_path is None else out_path
    print(
        f"Wrote {len(rows_out)} row(s) to {dest} "
        f"({len(by_st)} structure(s); {len(all_recs)} contact line(s) parsed).",
        file=sys.stderr,
    )


if __name__ == "__main__":
    main()
