#!/usr/bin/env python3
"""
Post-process FoldKit lattice packing outputs into a single TXT and/or CSV.

This is a companion to `metrics/lattice_packing_analyser.py`.

It can:
- Concatenate multiple per-structure packing TXT reports into one combined TXT.
- Convert per-structure packing JSON outputs (or a combined JSON bundle) into a CSV table.

Examples (from repository root, in the directory containing outputs):

  # If lattice_packing_analyser wrote per-structure files:
  python metrics/lattice_packing_report_csv.py packing_out/ -o lattice_packing.csv --output-txt lattice_packing.txt

  # Explicit globs:
  python metrics/lattice_packing_report_csv.py --json-glob 'packing_out/packing_*.json' -o packing.csv
  python metrics/lattice_packing_report_csv.py --txt-glob 'packing_out/packing_*.txt' --output-txt packing.txt

  # If lattice_packing_analyser wrote a combined JSON bundle:
  python metrics/lattice_packing_report_csv.py lattice_packing_combined.json -o lattice_packing.csv
"""

from __future__ import annotations

import argparse
import csv
import glob
import json
import os
import sys
from typing import Any, Iterable

_REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if _REPO_ROOT not in sys.path:
    sys.path.insert(0, _REPO_ROOT)

from utils.cli_log import add_log_args, setup_log_from_args


_PRESET_COLUMNS: dict[str, list[str]] = {
    # Matches the key summary metrics users typically want.
    "summary": [
        "input",
        "chain_count",
        "residue_count",
        "total_atoms",
        "volume_source",
        "volume_volume_a3",
        "lattice_atom_density_atoms_per_a3",
        "lattice_mass_density_da_per_a3",
        "lattice_packing_density_fraction",
        "lattice_packing_density_percent",
        "lattice_matthews_a3_per_da",
        "estimated_solvent_content_fraction",
        "estimated_solvent_content_percent",
        "estimated_protein_volume_a3",
        "estimated_solvent_volume_a3",
        "total_atomic_volume_a3",
        "total_mass_da",
        "volume_max",
        "volume_min",
        "volume_pad_a",
        "volume_span",
        "error",
    ],
}


def _iter_paths(inputs: list[str]) -> list[str]:
    out: list[str] = []
    for raw in inputs:
        if not raw:
            continue
        # Expand globs explicitly so quoting is optional in shells.
        matches = glob.glob(raw)
        if matches:
            out.extend(matches)
        else:
            out.append(raw)
    # Normalise, drop duplicates while preserving order.
    seen: set[str] = set()
    normed: list[str] = []
    for p in out:
        ap = os.path.abspath(p)
        if ap in seen:
            continue
        seen.add(ap)
        normed.append(ap)
    return normed


def _discover_from_dirs(
    paths: list[str],
    json_patterns: list[str],
    txt_patterns: list[str],
) -> tuple[list[str], list[str]]:
    json_files: list[str] = []
    txt_files: list[str] = []
    for p in paths:
        if os.path.isdir(p):
            for pat in json_patterns:
                json_files.extend(glob.glob(os.path.join(p, pat)))
            for pat in txt_patterns:
                txt_files.extend(glob.glob(os.path.join(p, pat)))
        elif os.path.isfile(p):
            if p.lower().endswith(".json"):
                json_files.append(p)
            elif p.lower().endswith(".txt"):
                txt_files.append(p)
        else:
            # Non-existent: ignore (matches interface_molecule_report_csv behaviour: hard error is noisy here)
            continue
    json_files = sorted({os.path.abspath(x) for x in json_files})
    txt_files = sorted({os.path.abspath(x) for x in txt_files})
    return json_files, txt_files


def _flatten(obj: Any, prefix: str = "") -> dict[str, Any]:
    """
    Flatten nested dicts to a 1-level mapping with underscore keys.
    Lists are JSON-encoded to keep CSV cells scalar.
    """
    out: dict[str, Any] = {}
    if isinstance(obj, dict):
        for k, v in obj.items():
            key = f"{prefix}{k}" if not prefix else f"{prefix}_{k}"
            if isinstance(v, dict):
                out.update(_flatten(v, key))
            elif isinstance(v, list):
                out[key] = json.dumps(v, ensure_ascii=False)
            else:
                out[key] = v
        return out
    out[prefix or "value"] = obj
    return out


def _load_json_records(json_paths: Iterable[str]) -> list[dict[str, Any]]:
    """
    Return a list of flat per-structure rows.

    Supports:
    - Per-structure JSON output from lattice_packing_analyser: a single dict of results.
    - Combined bundle JSON output from lattice_packing_analyser: {n_structures, structures:[{set_label, structure_basename, structure_stem, result:{...}}]}
    """
    rows: list[dict[str, Any]] = []
    for p in json_paths:
        try:
            with open(p, "r", encoding="utf-8", errors="replace") as f:
                data = json.load(f)
        except Exception as e:
            rows.append({"source_json": os.path.basename(p), "error": f"Failed to read JSON: {e}"})
            continue

        if isinstance(data, dict) and "structures" in data and isinstance(data.get("structures"), list):
            for rec in data.get("structures") or []:
                if not isinstance(rec, dict):
                    continue
                base: dict[str, Any] = {
                    "set_label": rec.get("set_label", ""),
                    "structure_basename": rec.get("structure_basename", ""),
                    "structure_stem": rec.get("structure_stem", ""),
                    "source_json": os.path.basename(p),
                }
                res = rec.get("result")
                flat = _flatten(res if isinstance(res, dict) else {"value": res})
                rows.append({**base, **flat})
            continue

        # Per-structure dict: add a few helpful context fields derived from the result.
        if isinstance(data, dict):
            input_path = str(data.get("input", "") or "")
            base = {
                "set_label": "",
                "structure_basename": os.path.basename(input_path) if input_path else "",
                "structure_stem": os.path.splitext(os.path.basename(input_path))[0] if input_path else "",
                "source_json": os.path.basename(p),
            }
            rows.append({**base, **_flatten(data)})
            continue

        rows.append({"source_json": os.path.basename(p), "error": "Unrecognised JSON format"})
    return rows


def _to_float_or_none(x: Any) -> float | None:
    try:
        return float(x)
    except Exception:
        return None


def _ensure_percent_fields(rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    """
    Guarantee percent fields exist when fraction fields are present.
    This keeps downstream CSV consistent even if upstream JSON changes.
    """
    out: list[dict[str, Any]] = []
    for r in rows:
        rr = dict(r)

        if rr.get("lattice_packing_density_percent", None) in (None, ""):
            frac = _to_float_or_none(rr.get("lattice_packing_density_fraction", None))
            if frac is not None:
                rr["lattice_packing_density_percent"] = frac * 100.0

        if rr.get("estimated_solvent_content_percent", None) in (None, ""):
            frac = _to_float_or_none(rr.get("estimated_solvent_content_fraction", None))
            if frac is not None:
                rr["estimated_solvent_content_percent"] = frac * 100.0

        out.append(rr)
    return out


def _write_csv(
    rows: list[dict[str, Any]],
    out_csv_path: str | None,
    *,
    columns: list[str] | None = None,
) -> None:
    if not rows:
        raise SystemExit("No rows to write.")

    all_cols = sorted({k for r in rows for k in r.keys()})
    if columns is not None:
        cols = [c for c in columns if c in all_cols]
        # Keep explicitly requested columns even if currently absent (write blank cells).
        missing = [c for c in columns if c not in cols]
        cols = cols + missing
    else:
        # Stable, readable column order.
        preferred = [
            "set_label",
            "structure_basename",
            "structure_stem",
            "source_json",
            "input",
            "chain_count",
            "residue_count",
            "total_atoms",
            "volume_source",
            "volume_volume_a3",
            "lattice_atom_density_atoms_per_a3",
            "lattice_mass_density_da_per_a3",
            "lattice_packing_density_fraction",
            "lattice_packing_density_percent",
            "lattice_matthews_a3_per_da",
            "estimated_solvent_content_fraction",
            "estimated_solvent_content_percent",
            "warning",
            "warning_volume_source",
            "error",
        ]
        cols = [c for c in preferred if c in all_cols] + [c for c in all_cols if c not in preferred]

    out_f = sys.stdout if (out_csv_path is None or out_csv_path == "-") else open(out_csv_path, "w", newline="", encoding="utf-8")
    try:
        w = csv.DictWriter(out_f, fieldnames=cols, extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)
    finally:
        if out_f is not sys.stdout:
            out_f.close()


def _concat_txt(txt_paths: list[str], out_txt_path: str) -> None:
    os.makedirs(os.path.dirname(os.path.abspath(out_txt_path)) or ".", exist_ok=True)
    with open(out_txt_path, "w", encoding="utf-8", newline="\n") as out:
        out.write("FOLDKIT LATTICE PACKING (CONCATENATED REPORT)\n\n")
        for p in sorted(txt_paths):
            out.write("=" * 72 + "\n")
            out.write(f"{os.path.basename(p)}\n")
            out.write("=" * 72 + "\n")
            try:
                with open(p, "r", encoding="utf-8", errors="replace") as f:
                    out.write(f.read().rstrip() + "\n\n")
            except Exception as e:
                out.write(f"LATTICE PACKING ANALYSIS\nError: Failed to read TXT: {e}\n\n")


def main() -> None:
    ap = argparse.ArgumentParser(
        description="Convert FoldKit lattice packing outputs into a combined TXT and/or CSV.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Notes:
- For CSV generation, JSON is preferred (more robust than parsing TXT).
- When a combined TXT/JSON was written directly from lattice_packing_analyser.py
  (single paths to --output-txt / --output-json), pass those paths here instead of raw analyser inputs.
""",
    )
    ap.add_argument(
        "inputs",
        nargs="*",
        help="Paths/globs/directories containing packing JSON/TXT files (default: current directory).",
    )
    ap.add_argument(
        "--json-glob",
        action="append",
        dest="json_globs",
        metavar="GLOB",
        help="Glob(s) used when scanning directories (repeatable). Default includes '*_lattice_packing.json' and 'packing_*.json'.",
    )
    ap.add_argument(
        "--txt-glob",
        action="append",
        dest="txt_globs",
        metavar="GLOB",
        help="Glob(s) used when scanning directories (repeatable). Default includes '*_lattice_packing.txt' and 'packing_*.txt'.",
    )
    ap.add_argument(
        "--output-txt",
        metavar="FILE",
        default=None,
        help="Write concatenated TXT report to this path.",
    )
    ap.add_argument(
        "--output-csv",
        "-o",
        metavar="FILE",
        default=None,
        help="Write CSV table to this path (or '-' for stdout).",
    )
    ap.add_argument(
        "--preset",
        choices=sorted(_PRESET_COLUMNS.keys()),
        default=None,
        help="CSV column preset. Currently supported: 'summary'.",
    )
    ap.add_argument(
        "--columns",
        metavar="COLS",
        default=None,
        help="Comma-separated list of CSV columns to write (overrides --preset).",
    )
    add_log_args(ap)
    args = ap.parse_args()
    setup_log_from_args(args, script_path=__file__, inputs=list(getattr(args, "inputs", []) or []), pattern=None)

    inputs = list(getattr(args, "inputs", []) or [])
    if not inputs:
        inputs = ["."]

    json_patterns = list(getattr(args, "json_globs", None) or []) or ["*_lattice_packing.json", "packing_*.json", "*.json"]
    txt_patterns = list(getattr(args, "txt_globs", None) or []) or ["*_lattice_packing.txt", "packing_*.txt", "*.txt"]

    expanded = _iter_paths(inputs)
    json_files, txt_files = _discover_from_dirs(expanded, json_patterns, txt_patterns)

    if not getattr(args, "output_txt", None) and not getattr(args, "output_csv", None):
        raise SystemExit("Specify at least one of --output-txt and/or --output-csv (-o).")

    if getattr(args, "output_txt", None):
        if not txt_files:
            raise SystemExit("No TXT files found to concatenate.")
        _concat_txt(txt_files, str(args.output_txt))

    if getattr(args, "output_csv", None):
        if not json_files:
            raise SystemExit("No JSON files found to convert to CSV.")
        rows = _ensure_percent_fields(_load_json_records(json_files))
        cols: list[str] | None = None
        if getattr(args, "columns", None):
            cols = [c.strip() for c in str(args.columns).split(",") if c.strip()]
        elif getattr(args, "preset", None):
            cols = list(_PRESET_COLUMNS.get(str(args.preset), []))
        _write_csv(rows, str(args.output_csv) if args.output_csv else None, columns=cols)


if __name__ == "__main__":
    main()

