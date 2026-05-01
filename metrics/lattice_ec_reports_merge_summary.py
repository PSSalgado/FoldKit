#!/usr/bin/env python3
"""
Merge several `interface_analyser_lattice_ec.py` text reports into one summary CSV.

Reuses the same summary-line parser as `interface_mol_report_ec_csv.py`, so columns
stay aligned with that tool. Adds one extra column:

  reference_chain_BSA_A2 =
      sasa_reference_isolated_A2 - sasa_reference_in_cluster_A2

(when both values are present).

Example::

  PYTHONPATH=/path/to/FoldKit python3 metrics/lattice_ec_reports_merge_summary.py \\
    /path/to/lattice_ec_*_results.txt \\
    -o /path/to/combined_ec_summary.csv

Paths may include shell globs (expand quotes so the shell expands the glob).
"""

from __future__ import annotations

import argparse
import csv
import glob
import os
import sys

_REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if _REPO_ROOT not in sys.path:
    sys.path.insert(0, _REPO_ROOT)

from metrics.interface_mol_report_ec_csv import _CSV_SUMMARY_FIELDS, parse_ec_report_text  # noqa: E402


def _expand_paths(inputs: list[str]) -> list[str]:
    out: list[str] = []
    for raw in inputs:
        raw = os.path.expanduser(raw.strip())
        if not raw:
            continue
        if any(ch in raw for ch in "*?["):
            matches = sorted(glob.glob(raw))
            if matches:
                out.extend(matches)
            else:
                print(f"Warning: skipping unmatched glob {raw!r}", file=sys.stderr)
        elif os.path.isfile(raw):
            out.append(raw)
        else:
            print(f"Warning: skipping missing path {raw!r}", file=sys.stderr)
    # stable unique order
    seen: set[str] = set()
    uniq: list[str] = []
    for p in out:
        ap = os.path.abspath(p)
        if ap not in seen:
            seen.add(ap)
            uniq.append(ap)
    return uniq


def _maybe_drop_redundant_all_columns(
    fieldnames: list[str],
    rows: list[dict[str, object]],
) -> list[str]:
    """
    In `--chains A`-style reports, the "all-chains interface summary" is already
    scoped by the focus-chain filter, so the `all_*` summary scalars become
    identical to the reference/top-level scalars. In that common case, keep the
    combined CSV compact by dropping those redundant columns entirely.
    """
    pairs = [
        ("all_total_interfaces", "total_interfaces", "int"),
        ("all_total_buried_surface_area_A2", "total_buried_surface_area_A2", "float"),
        ("all_average_buried_area_per_interface_A2", "average_buried_area_per_interface_A2", "float"),
        ("all_sasa_isolated_sum_A2", "sasa_isolated_sum_A2", "float"),
        ("all_sasa_isolated_by_chain", "sasa_isolated_by_chain", "str"),
    ]

    def _empty(v: object) -> bool:
        return v is None or v == ""

    def _equal(v_all: object, v_ref: object, kind: str) -> bool:
        if kind == "int":
            try:
                return int(v_all) == int(v_ref)
            except Exception:
                return False
        if kind == "float":
            try:
                return abs(float(v_all) - float(v_ref)) <= 1e-6
            except Exception:
                return False
        # str: exact compare after stringify (covers "A=123.4;B=..." ordering already normalised upstream)
        return str(v_all) == str(v_ref)

    # Only drop if for every pair, every row is either missing one side, or equal.
    for k_all, k_ref, kind in pairs:
        for r in rows:
            v_all = r.get(k_all)
            v_ref = r.get(k_ref)
            if _empty(v_all) or _empty(v_ref):
                continue
            if not _equal(v_all, v_ref, kind):
                return fieldnames

    drop = {k_all for (k_all, _k_ref, _kind) in pairs}
    return [f for f in fieldnames if f not in drop]


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(
        description="Merge lattice EC analyser reports into one summary CSV with reference-chain BSA column.",
    )
    ap.add_argument(
        "reports",
        nargs="+",
        metavar="PATH",
        help="One or more lattice_ec_*_results.txt paths (shell globs ok).",
    )
    ap.add_argument(
        "-o",
        "--output",
        required=True,
        metavar="CSV",
        help="Output CSV path.",
    )
    args = ap.parse_args(argv)

    paths = _expand_paths(list(args.reports))
    if not paths:
        ap.error("No report files found.")

    fieldnames = list(_CSV_SUMMARY_FIELDS)
    if "reference_chain_BSA_A2" not in fieldnames:
        insert_at = fieldnames.index("sasa_reference_in_cluster_A2") + 1
        fieldnames = fieldnames[:insert_at] + ["reference_chain_BSA_A2"] + fieldnames[insert_at:]
    rows: list[dict[str, object]] = []

    for path in paths:
        with open(path, encoding="utf-8", errors="replace") as f:
            text = f.read()
        _ifaces, summaries = parse_ec_report_text(text)
        if not summaries:
            print(f"Warning: no summary section parsed in {path}", file=sys.stderr)
            continue
        # One structure per report file → take first summary row.
        row = dict(summaries[0])
        row["report_txt"] = os.path.basename(str(path))
        iso = row.get("sasa_reference_isolated_A2")
        clu = row.get("sasa_reference_in_cluster_A2")
        try:
            if iso not in ("", None) and clu not in ("", None):
                row["reference_chain_BSA_A2"] = float(iso) - float(clu)
        except (TypeError, ValueError):
            row["reference_chain_BSA_A2"] = ""

        rows.append(row)

    rows.sort(key=lambda r: (str(r.get("structure_basename", "")), str(r.get("report_txt", ""))))

    fieldnames = _maybe_drop_redundant_all_columns(fieldnames, rows)

    os.makedirs(os.path.dirname(os.path.abspath(args.output)) or ".", exist_ok=True)
    with open(args.output, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames, extrasaction="ignore")
        w.writeheader()
        for r in rows:
            w.writerow({k: r.get(k, "") for k in fieldnames})

    print(f"Wrote {len(rows)} row(s) to {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
