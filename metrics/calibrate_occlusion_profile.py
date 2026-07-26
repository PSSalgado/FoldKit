#!/usr/bin/env python3
"""
Derive fixed 0-10 scoring bands for the three lattice *occlusion* metrics from
FoldKit-computed multi-copy lattices, and write them as a versioned profile JSON.

The occlusion metrics are BSAmolA/kDa, LOImolA (lattice burial fraction), and the
Mol. A lattice contact-residue fraction. Unlike the volume axes, these have no
directly transferable literature band: published crystal-packing statistics are
either per-interface (Janin & Rodier 1995; Luo et al. 2015) or describe 3-D packing
of monomeric proteins (Carugo & Argos 1997, 15-49% surface in contacts). FoldKit
instead measures one focal chain against *all* neighbours in a multi-copy 2-D
lattice (for example SlpA), where much of the molecule still faces solvent on the
environment and cell-wall sides, so absolute occlusion is far lower. The bands are
therefore calibrated empirically on FoldKit-computed multi-copy lattices.

Reads one or more metrics tables in either layout accepted by
``lattice_metrics_radar.load_rows`` (``combined_lattice_vs_ec.csv`` with structures
as rows, or a transposed summary table), pools the structures, and records for each
occlusion metric the observed span, the rounded scoring band, and the provenance.

Bands are rounded *outward*. If the 5-point grid would turn a positive empirical
LOI/contact minimum into zero, that lower edge uses a 1%-point grid instead,
preserving a nonzero baseline for an assembled lattice. Upper edges use the
coarser 5-point grid. BSA/kDa uses 5 units at both edges.

Examples
--------
Rebuild the default expanded-lattice profile from an experimental cohort::

    PYTHONPATH=. python metrics/calibrate_occlusion_profile.py \\
      --input ./expanded_lattices/combined_lattice_vs_ec.csv \\
      --profile bbox --version 1 \\
      --population "Expanded-lattice / multi-copy occlusion (calibrated on FoldKit multi-copy lattices)"

Restrict to a compact subset (tighter band for SlpA lattices)::

    PYTHONPATH=. python metrics/calibrate_occlusion_profile.py \\
      --input ./compact_slpa_15m/combined_lattice_vs_ec.csv \\
      --profile slayer-compact --version 1 \\
      --stems s2_opt2427_15mA,s7_cd630_15mA,s7b_r7404_15mA,s11_ox247o_15mA,s4del_r2d2_15mA
"""

from __future__ import annotations

import argparse
import json
import math
import os
import sys

_REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if _REPO_ROOT not in sys.path:
    sys.path.insert(0, _REPO_ROOT)

from metrics.lattice_metrics_radar import (  # noqa: E402
    METRICS,
    OCCLUSION_COLUMNS,
    PROFILE_DIR,
    _round_outward_asymmetric,
    _to_float,
    load_rows,
    occlusion_lower_round_step,
    occlusion_round_step,
    profile_filename,
)

_BAND_METHODS = ("minmax", "p5p95")


def _percentile(sorted_vals: list[float], q: float) -> float:
    """Linear-interpolation percentile (q in [0, 100]) on an already-sorted list."""
    if not sorted_vals:
        return math.nan
    if len(sorted_vals) == 1:
        return sorted_vals[0]
    pos = (q / 100.0) * (len(sorted_vals) - 1)
    lo = math.floor(pos)
    hi = math.ceil(pos)
    if lo == hi:
        return sorted_vals[int(pos)]
    frac = pos - lo
    return sorted_vals[lo] * (1.0 - frac) + sorted_vals[hi] * frac


def _display_for(column: str) -> str:
    for m in METRICS:
        if m.column == column:
            return m.display
    return column


def collect_values(
    inputs: list[str],
    *,
    stems: list[str] | None = None,
) -> tuple[dict[str, list[float]], list[dict[str, object]]]:
    """Pool occlusion values across input tables. Returns (values_by_column, per_structure)."""
    keep = set(stems or [])
    values: dict[str, list[float]] = {col: [] for col in OCCLUSION_COLUMNS}
    per_structure: list[dict[str, object]] = []
    seen: set[str] = set()

    for path in inputs:
        rows = load_rows(path)
        source = os.path.basename(path)
        for row in rows:
            stem = str(row.get("structure_stem", "")).strip()
            if not stem or (keep and stem not in keep):
                continue
            if stem in seen:
                print(f"Warning: duplicate structure {stem!r} ignored ({path}).", file=sys.stderr)
                continue
            seen.add(stem)
            record: dict[str, object] = {"structure_stem": stem, "source": source}
            usable = False
            for col in OCCLUSION_COLUMNS:
                val = _to_float(row.get(col))
                if math.isnan(val):
                    record[col] = None
                    continue
                values[col].append(val)
                record[col] = round(val, 4)
                usable = True
            if usable:
                per_structure.append(record)
            else:
                seen.discard(stem)

    if keep:
        missing = sorted(keep - seen)
        if missing:
            raise ValueError(f"--stems not found in the input table(s): {missing}")
    return values, per_structure


def build_profile(
    inputs: list[str],
    *,
    profile: str,
    version: int,
    band_method: str = "minmax",
    stems: list[str] | None = None,
    population: str = "",
    notes: list[str] | None = None,
) -> dict[str, object]:
    """Compute observed spans and rounded bands, returning the profile document."""
    if band_method not in _BAND_METHODS:
        raise ValueError(f"Unknown --band {band_method!r}; choose from {_BAND_METHODS}.")

    values, per_structure = collect_values(inputs, stems=stems)
    if not per_structure:
        raise ValueError("No structures with usable occlusion metrics in the input table(s).")

    metric_limits: dict[str, object] = {}
    for col in OCCLUSION_COLUMNS:
        vals = sorted(values[col])
        if not vals:
            raise ValueError(f"No usable values for {_display_for(col)} ({col}).")
        obs_lo, obs_hi = vals[0], vals[-1]
        if band_method == "p5p95":
            band_lo, band_hi = _percentile(vals, 5.0), _percentile(vals, 95.0)
        else:
            band_lo, band_hi = obs_lo, obs_hi
        lower_step = occlusion_lower_round_step(col, band_lo)
        upper_step = occlusion_round_step(col)
        lo, hi = _round_outward_asymmetric(
            band_lo,
            band_hi,
            lower_step=lower_step,
            upper_step=upper_step,
        )
        metric_limits[col] = {
            "display": _display_for(col),
            "n": len(vals),
            "observed_min": round(obs_lo, 4),
            "observed_max": round(obs_hi, 4),
            "band_min_unrounded": round(band_lo, 4),
            "band_max_unrounded": round(band_hi, 4),
            "lower_round_step": lower_step,
            "upper_round_step": upper_step,
            "limits": [lo, hi],
        }

    return {
        "profile": profile,
        "version": int(version),
        "band_method": band_method,
        "population": population,
        "rounding": (
            "Bands are rounded outward. When a 5-point grid would reduce a positive "
            "LOI/contact minimum to zero, the lower limit uses a 1-point grid instead; "
            "upper limits use a 5-point grid. BSA/kDa uses 5 units at both edges."
        ),
        "sources": [os.path.abspath(p) for p in inputs],
        "n_structures": len(per_structure),
        "metric_limits": metric_limits,
        "structures": per_structure,
        "notes": list(notes or []),
    }


def write_profile(doc: dict[str, object], *, output_dir: str | None = None) -> str:
    """Write the profile JSON next to the shipped profiles (or a chosen directory)."""
    target_dir = output_dir or PROFILE_DIR
    os.makedirs(target_dir, exist_ok=True)
    path = os.path.join(
        target_dir, profile_filename(str(doc["profile"]), int(doc["version"]))
    )
    with open(path, "w", encoding="utf-8") as f:
        json.dump(doc, f, indent=2)
        f.write("\n")
    return path


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(
        description=(
            "Derive fixed 0-10 scoring bands for BSAmolA/kDa, LOImolA, and lattice "
            "contact-residue % from FoldKit-computed multi-copy lattices, and write "
            "them as a versioned occlusion profile JSON."
        )
    )
    ap.add_argument(
        "--input",
        required=True,
        nargs="+",
        metavar="CSV",
        help="One or more metrics tables (combined_lattice_vs_ec.csv or a transposed summary table).",
    )
    ap.add_argument("--profile", required=True, metavar="NAME", help="Profile name, for example 'bbox'.")
    ap.add_argument("--version", type=int, default=1, help="Profile version number (default: 1).")
    ap.add_argument(
        "--band",
        choices=_BAND_METHODS,
        default="minmax",
        help=(
            "Band before rounding: 'minmax' uses the full observed span (default; "
            "appropriate for small curated cohorts), 'p5p95' trims to the 5th-95th percentile."
        ),
    )
    ap.add_argument(
        "--stems",
        default=None,
        metavar="STEM,STEM,...",
        help="Comma-separated subset of structure stems to calibrate on (default: all).",
    )
    ap.add_argument(
        "--population",
        default="",
        help="Short description of the reference population recorded in the profile.",
    )
    ap.add_argument(
        "--note",
        action="append",
        default=None,
        metavar="TEXT",
        help="Provenance/literature note to record (repeatable).",
    )
    ap.add_argument(
        "--output-dir",
        default=None,
        metavar="DIR",
        help=f"Where to write the profile (default: {PROFILE_DIR}).",
    )
    args = ap.parse_args(argv)

    stems = [s.strip() for s in str(args.stems).split(",") if s.strip()] if args.stems else None
    try:
        doc = build_profile(
            [os.path.abspath(os.path.expanduser(p)) for p in args.input],
            profile=args.profile,
            version=int(args.version),
            band_method=args.band,
            stems=stems,
            population=args.population,
            notes=args.note,
        )
        path = write_profile(doc, output_dir=args.output_dir)
    except ValueError as exc:
        ap.error(str(exc))
        return 2  # unreachable; ap.error raises SystemExit

    limits = doc["metric_limits"]
    assert isinstance(limits, dict)
    print(f"Wrote {path}")
    print(f"Profile {doc['profile']} v{doc['version']} from {doc['n_structures']} structures:")
    for col, info in limits.items():
        assert isinstance(info, dict)
        lo, hi = info["limits"]
        print(
            f"  {info['display']}: observed {info['observed_min']:g}-{info['observed_max']:g}"
            f" -> band {lo:g}-{hi:g} (lower step {info['lower_round_step']:g}, "
            f"upper step {info['upper_round_step']:g}, n={info['n']})"
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
