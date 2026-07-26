#!/usr/bin/env python3
"""
Radar (spider) charts of six key lattice packing/occlusion metrics.

Two input layouts are auto-detected:
  - structures as rows: ``combined_lattice_vs_ec.csv`` from
    ``metrics/lattice_compare_batch.py`` (``structure_stem`` + metric columns);
  - metrics as rows: a transposed summary table whose first column holds metric
    labels (``Packing density (%)``, ``BSAmolA/kDa (Å²)``, ...) and whose header
    row lists the structure names.

Labels are matched ignoring case, punctuation, and unit glyphs, so snake_case
column names and human-readable display labels both work.

Each metric is mapped to a configurable score range (default 0-10) against a
fixed raw-value band, so scores stay comparable between runs that use the same
score range. Volume and occlusion axes share the same scale
vocabulary via ``--scale-mode`` / ``--occlusion-scale-mode``:
``cryst1``, ``bbox``, ``empirical``, ``slayer-compact``, ``user``
(default ``bbox``). See ``metrics/radar_profiles/``. Values outside a band are
clipped to the score-range endpoints and recorded in ``radar_scale.json``. Two metrics are inverted
(Matthews, estimated solvent) so that a *larger polygon always means a more
occluded / tighter lattice*.

The scale-mode tag is appended to ``--output-dir`` and carried in every
filename (for example ``<output-dir>_bbox/radar_grid_bbox.png``), so runs with
different scaling never overwrite one another without adding a nested folder.
When the volume and
occlusion families use different modes the tag records both
(``vol-<mode>_occ-<mode>``, for example ``vol-slayer-compact_occ-empirical``). Pass
``--no-scale-tag`` to write flat, untagged names directly in ``--output-dir``.
Files (``<tag>`` shown, omitted with ``--no-scale-tag``):
  - radar_grid_<tag>.<fmt>            small-multiples radar, one panel per
                                      structure (paginated as
                                      radar_grid_<tag>_01.<fmt>, ... when
                                      ``--max-per-page`` splits the cohort)
  - radar_scores_heatmap_<tag>.<fmt>  structures x six metrics (optional;
                                      omit with ``--no-heatmap``)
  - radar_scores_<tag>.csv            per-structure plot scores + raw values
  - radar_scale_<tag>.json            resolved bands, clipping, and scale tag

Comparison across proteins is by design without overlaying polygons: every panel
shares identical geometry and score rings, carries a faint dashed ghost silhouette
(leave-one-out median by default, or ``--reference-stem`` / ``--ghost``), and an
optional companion score heatmap gives the quantitative cross-structure view. Use
``--display deviation`` to plot signed cohort deviations instead of band scores
when within-cohort differences are a small fraction of the physical band.

Requires PYTHONPATH pointing at the FoldKit repo root (same as other metrics
scripts) and matplotlib for figure output.
"""

from __future__ import annotations

import argparse
import csv
import math
import os
import re
import sys
from typing import NamedTuple


class MetricSpec(NamedTuple):
    column: str
    display: str
    short: str
    invert: bool
    aliases: tuple[str, ...] = ()


# Default axis order (clockwise from LOI at 12 o'clock): LOI -> Matthews ->
# contact -> packing -> BSA/kDa -> solvent. ``aliases`` let transposed summary
# tables use human-readable row labels, and also accept the upstream analyser
# field names. Min-max scaling is invariant to a constant factor, so
# fraction-vs-percent aliases score identically.
METRICS: tuple[MetricSpec, ...] = (
    MetricSpec(
        "lattice_burial_fraction_percent",
        "LOImolA (%)",
        "LOI (%)",
        False,
        (
            "Lattice occlusion index (LOImolA) (%)",
            "Lattice occlusion index (%)",
            "lattice_burial_fraction",
        ),
    ),
    MetricSpec(
        "matthews_a3_per_Da",
        "Matthews-like (A^3/Da)",
        "Matthews\n(A3/Da, inv)",
        True,
        (
            "Matthews-like coefficient (Å3/Da)",
            "Matthews coefficient (Å3/Da)",
            "lattice_matthews_a3_per_da",
            "matthews_coefficient",
        ),
    ),
    MetricSpec(
        "lattice_contact_residue_fraction_percent",
        "Mol. A lattice contact residues (%)",
        "Contact\nres. (%)",
        False,
        ("Lattice contact residues (%)", "lattice_contact_residue_fraction"),
    ),
    MetricSpec(
        "packing_density_percent",
        "Packing density (%)",
        "Packing\ndensity (%)",
        False,
        ("lattice_packing_density_percent", "packing_efficiency_percent"),
    ),
    MetricSpec(
        "reference_chain_BSA_per_kDa_reference_chain_A2",
        "BSAmolA/kDa (A^2)",
        "BSA/kDa\n(A2)",
        False,
        ("BSA/kDa (Å²)", "reference_chain_BSA_per_kDa_A2"),
    ),
    MetricSpec(
        "estimated_solvent_percent",
        "Estimated solvent (%)",
        "Solvent\n(%, inv)",
        True,
        ("estimated_solvent_content_percent", "solvent_content_percent"),
    ),
)

# Short tokens accepted by --metric-order (plus full column names).
_METRIC_ORDER_ALIASES: dict[str, str] = {
    "packing": "packing_density_percent",
    "packing_density": "packing_density_percent",
    "packing_density_percent": "packing_density_percent",
    "matthews": "matthews_a3_per_Da",
    "matthews_a3_per_da": "matthews_a3_per_Da",
    "solvent": "estimated_solvent_percent",
    "estimated_solvent": "estimated_solvent_percent",
    "estimated_solvent_percent": "estimated_solvent_percent",
    "bsa": "reference_chain_BSA_per_kDa_reference_chain_A2",
    "bsa_kda": "reference_chain_BSA_per_kDa_reference_chain_A2",
    "bsa/kda": "reference_chain_BSA_per_kDa_reference_chain_A2",
    "reference_chain_bsa_per_kda_reference_chain_a2": (
        "reference_chain_BSA_per_kDa_reference_chain_A2"
    ),
    "loi": "lattice_burial_fraction_percent",
    "loimola": "lattice_burial_fraction_percent",
    "lattice_burial_fraction_percent": "lattice_burial_fraction_percent",
    "contact": "lattice_contact_residue_fraction_percent",
    "contacts": "lattice_contact_residue_fraction_percent",
    "lattice_contact_residue_fraction_percent": (
        "lattice_contact_residue_fraction_percent"
    ),
}

# Preset spoke orders. ``interleaved`` alternates volume and occlusion with
# biophysically matched neighbours (packing↔contact, Matthews↔LOI, solvent↔BSA)
# rather than zigzagging the default block order.
METRIC_ORDER_PRESETS: dict[str, tuple[str, ...]] = {
    "default": tuple(m.column for m in METRICS),
    "interleaved": (
        "packing_density_percent",
        "lattice_contact_residue_fraction_percent",
        "matthews_a3_per_Da",
        "lattice_burial_fraction_percent",
        "estimated_solvent_percent",
        "reference_chain_BSA_per_kDa_reference_chain_A2",
    ),
}
_METRICS_BY_COLUMN: dict[str, MetricSpec] = {m.column: m for m in METRICS}

_VALID_FORMATS = ("png", "svg", "pdf", "tiff")
DEFAULT_SCORE_RANGE = (0.0, 10.0)

# Radial encodings. 'score' maps raw values onto the fixed band; 'deviation'
# standardises each axis about the cohort centre, which is the only way to
# resolve cohorts whose spread is a small fraction of a physically-motivated band.
DISPLAY_MODES = ("score", "deviation")
DEVIATION_SCALES = ("sd", "mad", "range")
# 1 / Phi^-1(3/4): scales the MAD to a standard-deviation equivalent for a normal sample.
_MAD_TO_SD = 1.4826
DEFAULT_DEVIATION_LIMIT = 0.0  # 0 = auto-fit the ring to the observed spread

# Ghost silhouette references. The leave-one-out variants exclude the panel's own
# structure, so with an odd cohort the ghost can no longer land exactly on the
# polygon it is meant to be compared against.
GHOST_MODES = ("loo-median", "loo-mean", "cohort-median", "cohort-mean")

# Crystallographic literature bands (unit-cell / CRYST1; no expansion).
# Matthews & solvent: Matthews 1968; Kantardjieff & Rupp 2003 (typical practical VM band).
# Packing: crystallographic packing coefficient ≈ 1 − solvent (Andersson & Hovmöller 2000).
CRYSTAL_MATTHEWS_A3_PER_DA = (1.5, 4.0)
CRYSTAL_SOLVENT_PERCENT = (25.0, 80.0)
CRYSTAL_PACKING_PERCENT = (20.0, 75.0)
# FoldKit heuristic protein density (Da/Å³) used in lattice_packing_analyser solvent %.
_PROTEIN_DENSITY_DA_PER_A3 = 0.81

_VOLUME_COLUMNS = (
    "packing_density_percent",
    "matthews_a3_per_Da",
    "estimated_solvent_percent",
)

# Shared vocabulary for --scale-mode and --occlusion-scale-mode.
SCALE_MODES = ("cryst1", "bbox", "empirical", "slayer-compact", "user")
# Occlusion keeps deprecated 'cohort' (raw within-run min-max, no padding).
OCCLUSION_SCALE_MODES = (*SCALE_MODES, "cohort")
_SCALE_MODE_ALIASES = {"crystal": "cryst1"}

# Occlusion axes: whole-neighbourhood burial of the reference chain in a multi-copy
# lattice. These have no transferable published band (see metrics_details.md), so the
# defaults are empirical profiles calibrated on FoldKit-computed multi-copy lattices
# (including SlpA; see metrics_details.md).
OCCLUSION_COLUMNS = (
    "reference_chain_BSA_per_kDa_reference_chain_A2",
    "lattice_burial_fraction_percent",
    "lattice_contact_residue_fraction_percent",
)
# Columns whose fixed bands are expressed in percent, so a fraction-valued input column
# (0-1 rather than 0-100) would be scored against the wrong scale.
_PERCENT_COLUMNS = (
    "packing_density_percent",
    "estimated_solvent_percent",
    "lattice_burial_fraction_percent",
    "lattice_contact_residue_fraction_percent",
)
PROFILE_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "radar_profiles")
# Shipped profile versions: occlusion uses bare name (bbox_v1.json);
# volume uses "{name}_volume_vN.json".
_OCCLUSION_PROFILE_VERSION: dict[str, int] = {"bbox": 1, "slayer-compact": 1}
# Occlusion mode → profile name.
_OCCLUSION_PROFILE_FOR_MODE: dict[str, str] = {
    "bbox": "bbox",
    "slayer-compact": "slayer-compact",
}
_VOLUME_PROFILE_VERSION: dict[str, int] = {"bbox": 1, "slayer-compact": 1}
# Within-run empirical bands: span at or below this uses ±1 integer padding;
# larger spans round outward onto the column reporting grid.
EMPIRICAL_CLOSE_SPAN = 5.0

# Backward-compatible alias used by older call sites / docs snippets.
_SCALE_MODES = SCALE_MODES


def _normalize_scale_mode(mode: str | None, *, default: str) -> str:
    raw = (mode or default).strip().lower()
    return _SCALE_MODE_ALIASES.get(raw, raw)


def _mode_slug(mode: str) -> str:
    """Filesystem-safe slug for a scale-mode name (for example ``slayer-compact``)."""
    return re.sub(r"[^a-z0-9]+", "-", str(mode).strip().lower()).strip("-")


def resolve_metric_order(order: str | list[str] | tuple[str, ...] | None) -> tuple[MetricSpec, ...]:
    """
    Resolve a spoke order from a preset name or an explicit token / column list.

    Presets: ``default`` (LOI, Matthews, contact, packing, BSA/kDa, solvent),
    ``interleaved`` (packing, contact, Matthews, LOI, solvent, BSA/kDa —
    matched neighbours). A custom order is a comma-separated string or sequence
    of short tokens / full column names that must list each of the six metrics
    exactly once.
    """
    if order is None or order == "":
        return METRICS
    if isinstance(order, (list, tuple)):
        tokens = [str(t).strip() for t in order if str(t).strip()]
        preset_key = None
    else:
        raw = str(order).strip()
        preset_key = raw.lower()
        if preset_key in METRIC_ORDER_PRESETS:
            tokens = list(METRIC_ORDER_PRESETS[preset_key])
        else:
            tokens = [tok.strip() for tok in raw.split(",") if tok.strip()]

    if not tokens:
        raise ValueError(
            "--metric-order is empty; use a preset "
            f"({', '.join(METRIC_ORDER_PRESETS)}) or six comma-separated tokens."
        )

    columns: list[str] = []
    for tok in tokens:
        key = tok.strip().lower().replace(" ", "_")
        col = _METRIC_ORDER_ALIASES.get(key)
        if col is None:
            # Accept exact column name with original case.
            if tok in _METRICS_BY_COLUMN:
                col = tok
            else:
                known = ", ".join(METRIC_ORDER_PRESETS) + "; tokens: packing, matthews, solvent, bsa, loi, contact"
                raise ValueError(f"Unknown metric-order token {tok!r}. Known: {known}.")
        columns.append(col)

    if len(columns) != len(METRICS):
        raise ValueError(
            f"--metric-order must list all {len(METRICS)} metrics exactly once "
            f"(got {len(columns)}: {columns})."
        )
    if len(set(columns)) != len(columns):
        raise ValueError(f"--metric-order has duplicate metrics: {columns}.")
    missing = [m.column for m in METRICS if m.column not in set(columns)]
    if missing:
        raise ValueError(f"--metric-order is missing: {missing}.")
    return tuple(_METRICS_BY_COLUMN[c] for c in columns)


def metric_order_tag(metrics: tuple[MetricSpec, ...]) -> str:
    """Short tag for non-default spoke orders (empty string for the default)."""
    cols = tuple(m.column for m in metrics)
    if cols == METRIC_ORDER_PRESETS["default"]:
        return ""
    for name, preset in METRIC_ORDER_PRESETS.items():
        if cols == preset:
            return f"order-{_mode_slug(name)}"
    return "order-custom"


# Column placed at 12 o'clock on the radar (angle 0 with theta_offset = π/2).
DEFAULT_SPOKE_TOP_COLUMN = "lattice_burial_fraction_percent"


def spoke_angles(
    metrics: tuple[MetricSpec, ...],
    *,
    top_column: str = DEFAULT_SPOKE_TOP_COLUMN,
) -> list[float]:
    """
    Polar angles for ``metrics`` in clockwise order, with ``top_column`` at the top.

    Matplotlib polar uses ``theta_offset=π/2`` and ``theta_direction=-1`` in the
    panel drawer, so angle 0 is 12 o'clock. Metric / CSV / heatmap order is
    unchanged; only the drawn spoke positions rotate.
    """
    n = len(metrics)
    if n == 0:
        return []
    top_i = 0
    for i, m in enumerate(metrics):
        if m.column == top_column:
            top_i = i
            break
    step = 2.0 * math.pi / float(n)
    return [((i - top_i) % n) * step for i in range(n)]


def scale_mode_tag(volume_mode: str, occlusion_mode: str) -> str:
    """Directory / filename tag describing the scale modes in use.

    When both metric families share a mode the tag is that mode name
    (for example ``slayer-compact``). When they differ the tag carries both, so the
    output directory records the mix (for example ``vol-slayer-compact_occ-empirical``).
    """
    v = _mode_slug(volume_mode)
    o = _mode_slug(occlusion_mode)
    if v == o:
        return v
    return f"vol-{v}_occ-{o}"


def output_run_tag(
    volume_mode: str,
    occlusion_mode: str,
    *,
    display_mode: str = "score",
    deviation_scale: str = "sd",
    ghost_mode: str = "loo-median",
    reference_stem: str | None = None,
    metrics: tuple[MetricSpec, ...] = METRICS,
    edge_delta: bool = False,
) -> str:
    """Compose the full output tag, including display / ghost / order choices when non-default."""
    parts = [scale_mode_tag(volume_mode, occlusion_mode)]
    order_tag = metric_order_tag(metrics)
    if order_tag:
        parts.append(order_tag)
    if display_mode == "deviation":
        parts.append(f"dev-{_mode_slug(deviation_scale)}")
    if edge_delta:
        parts.append("edges-delta")
    if reference_stem:
        parts.append(f"ref-{_mode_slug(reference_stem)}")
    elif ghost_mode != "loo-median":
        parts.append(f"ghost-{_mode_slug(ghost_mode)}")
    return "_".join(p for p in parts if p)


def _tagged_name(stem: str, ext: str, tag: str, *, page: int | None = None) -> str:
    """Compose an output filename, inserting ``tag`` before the extension."""
    suffix = f"_{tag}" if tag else ""
    page_suffix = f"_{page:02d}" if page is not None else ""
    return f"{stem}{suffix}{page_suffix}.{ext}"


def _round_outward(lo: float, hi: float, *, step: float) -> tuple[float, float]:
    """Expand [lo, hi] outward to the nearest multiples of ``step``."""
    if step <= 0 or not math.isfinite(step):
        raise ValueError(f"Invalid rounding step: {step}")
    if not (math.isfinite(lo) and math.isfinite(hi)):
        raise ValueError(f"Non-finite limits for rounding: {lo}, {hi}")
    if hi < lo:
        lo, hi = hi, lo
    lo_r = step * math.floor(lo / step)
    hi_r = step * math.ceil(hi / step)
    if hi_r <= lo_r:
        hi_r = lo_r + step
    return lo_r, hi_r


def _round_outward_asymmetric(
    lo: float,
    hi: float,
    *,
    lower_step: float,
    upper_step: float,
) -> tuple[float, float]:
    """Expand a band using independent lower- and upper-edge reporting grids."""
    if lower_step <= 0 or upper_step <= 0:
        raise ValueError(
            f"Invalid asymmetric rounding steps: lower={lower_step}, upper={upper_step}"
        )
    if not (math.isfinite(lo) and math.isfinite(hi)):
        raise ValueError(f"Non-finite limits for rounding: {lo}, {hi}")
    if hi < lo:
        lo, hi = hi, lo
    lo_r = lower_step * math.floor(lo / lower_step)
    hi_r = upper_step * math.ceil(hi / upper_step)
    if hi_r <= lo_r:
        hi_r = lo_r + upper_step
    return lo_r, hi_r


def _round_outward_to_0_or_5(lo: float, hi: float) -> tuple[float, float]:
    """Expand [lo, hi] to the next outward multiples of 5 (…, 0, 5, 10, 15, …)."""
    return _round_outward(lo, hi, step=5.0)


def _round_outward_to_half(lo: float, hi: float) -> tuple[float, float]:
    """Expand [lo, hi] to the next outward multiples of 0.5 (Matthews grid)."""
    return _round_outward(lo, hi, step=0.5)


# Reporting grid for each axis: limits are rounded to the nearest grid value below the
# band minimum and above the band maximum, whatever the scale mode.
_ROUND_STEP_BY_COLUMN: dict[str, float] = {
    "matthews_a3_per_Da": 0.5,
    "estimated_solvent_percent": 5.0,
    "packing_density_percent": 5.0,
    "reference_chain_BSA_per_kDa_reference_chain_A2": 5.0,
    "lattice_burial_fraction_percent": 5.0,
    "lattice_contact_residue_fraction_percent": 5.0,
}


def _round_limits_to_grid(
    limits: dict[str, tuple[float, float]],
    *,
    skip: set[str] | None = None,
) -> dict[str, tuple[float, float]]:
    """Round each band outward onto its reporting grid, leaving ``skip`` columns exact."""
    skipped = skip or set()
    out: dict[str, tuple[float, float]] = {}
    for col, (lo, hi) in limits.items():
        if col in skipped:
            out[col] = (float(lo), float(hi))
            continue
        step = _ROUND_STEP_BY_COLUMN.get(col)
        out[col] = _round_outward(lo, hi, step=step) if step else (float(lo), float(hi))
    return out


def _solvent_from_matthews(vm: float) -> float:
    """FoldKit heuristic: solvent% = 100 * (1 - 1/(0.81 * VM))."""
    if not math.isfinite(vm) or vm <= 0:
        return math.nan
    return 100.0 * (1.0 - 1.0 / (_PROTEIN_DENSITY_DA_PER_A3 * vm))


def resolve_volume_limits(
    *,
    scale_mode: str,
    matthews_values: list[float],
    packing_values: list[float] | None = None,
    solvent_values: list[float] | None = None,
    matthews_min: float | None = None,
    matthews_max: float | None = None,
    solvent_min: float | None = None,
    solvent_max: float | None = None,
    packing_min: float | None = None,
    packing_max: float | None = None,
    profile_dir: str | None = None,
) -> tuple[dict[str, tuple[float, float]], dict[str, object]]:
    """
    Resolve (min, max) raw-value ranges for the three volume metrics.

    scale_mode (shared vocabulary with occlusion):
      - cryst1: literature unit-cell bands (alias: crystal)
      - bbox: fixed expanded-lattice volume profile
      - slayer-compact: tighter SlpA-core volume profile
      - empirical: within-run span with ±1 integer pad (close) or reporting-grid
        outward rounding (wide)
      - user: requires --matthews-min/--matthews-max; solvent/packing optional

    Explicitly supplied --*-min/--*-max values win and stay exact.
    """
    mode = _normalize_scale_mode(scale_mode, default="bbox")
    if mode not in SCALE_MODES:
        raise ValueError(f"Unknown scale mode {scale_mode!r}; choose from {SCALE_MODES}.")

    meta: dict[str, object] = {"scale_mode": mode}
    limits: dict[str, tuple[float, float]] = {}
    skip_round: set[str] = set()

    if mode == "cryst1":
        limits = {
            "matthews_a3_per_Da": CRYSTAL_MATTHEWS_A3_PER_DA,
            "estimated_solvent_percent": CRYSTAL_SOLVENT_PERCENT,
            "packing_density_percent": CRYSTAL_PACKING_PERCENT,
        }
    elif mode in ("bbox", "slayer-compact"):
        doc = load_volume_profile(mode, profile_dir=profile_dir)
        metric_limits = doc["metric_limits"]
        assert isinstance(metric_limits, dict)
        for col in _VOLUME_COLUMNS:
            lo, hi = metric_limits[col]["limits"]
            limits[col] = (float(lo), float(hi))
            skip_round.add(col)  # profiles already store rounded bands
        meta["volume_profile"] = {
            "profile": doc.get("profile"),
            "version": doc.get("version"),
            "path": doc.get("_path"),
            "population": doc.get("population", ""),
            "n_structures": doc.get("n_structures"),
            "notes": doc.get("notes", []),
        }
    elif mode == "empirical":
        value_map = {
            "matthews_a3_per_Da": matthews_values,
            "estimated_solvent_percent": solvent_values or [],
            "packing_density_percent": packing_values or [],
        }
        # Solvent can be derived from Matthews when the input column is absent.
        if not any(math.isfinite(v) for v in value_map["estimated_solvent_percent"]):
            value_map["estimated_solvent_percent"] = [
                _solvent_from_matthews(v)
                for v in matthews_values
                if math.isfinite(v)
            ]
        per_metric: dict[str, object] = {}
        for col in _VOLUME_COLUMNS:
            vals = [v for v in value_map[col] if math.isfinite(v)]
            if not vals:
                raise ValueError(f"--scale-mode empirical needs usable values for {col}.")
            wide_step = _ROUND_STEP_BY_COLUMN.get(col, 5.0)
            obs_lo = min(vals)
            lower_step = volume_lower_round_step(col, obs_lo)
            lo, hi, rule = empirical_band(
                vals, wide_step=wide_step, lower_step=lower_step
            )
            limits[col] = (lo, hi)
            skip_round.add(col)
            per_metric[col] = {
                "observed_min": round(min(vals), 4),
                "observed_max": round(max(vals), 4),
                "span": round(max(vals) - min(vals), 4),
                "rule": rule,
                "limits": [lo, hi],
                "n": len(vals),
            }
        meta["empirical_bands"] = per_metric
        meta["empirical_close_span"] = EMPIRICAL_CLOSE_SPAN
    else:  # user
        if matthews_min is None or matthews_max is None:
            raise ValueError(
                "--scale-mode user requires --matthews-min and --matthews-max "
                "(optional --solvent-* / --packing-* otherwise derived from Matthews)."
            )
        vm_lo, vm_hi = float(matthews_min), float(matthews_max)
        if vm_hi <= vm_lo:
            raise ValueError(f"Matthews max must exceed min (got {vm_lo}, {vm_hi}).")
        if solvent_min is None or solvent_max is None:
            s_at_lo = _solvent_from_matthews(vm_lo)
            s_at_hi = _solvent_from_matthews(vm_hi)
            sol_lo, sol_hi = sorted((s_at_lo, s_at_hi))
        else:
            sol_lo, sol_hi = float(solvent_min), float(solvent_max)
        if packing_min is None or packing_max is None:
            pack_lo, pack_hi = sorted((100.0 - sol_hi, 100.0 - sol_lo))
        else:
            pack_lo, pack_hi = float(packing_min), float(packing_max)
        limits = {
            "matthews_a3_per_Da": (vm_lo, vm_hi),
            "estimated_solvent_percent": (sol_lo, sol_hi),
            "packing_density_percent": (pack_lo, pack_hi),
        }

    # Explicitly supplied bands are used verbatim; everything else snaps to the grid
    # unless already finalised by a profile / empirical rule.
    explicit: set[str] = set()
    if matthews_min is not None and matthews_max is not None:
        explicit.add("matthews_a3_per_Da")
        limits["matthews_a3_per_Da"] = (float(matthews_min), float(matthews_max))
    if solvent_min is not None and solvent_max is not None:
        explicit.add("estimated_solvent_percent")
        limits["estimated_solvent_percent"] = (float(solvent_min), float(solvent_max))
    if packing_min is not None and packing_max is not None:
        explicit.add("packing_density_percent")
        limits["packing_density_percent"] = (float(packing_min), float(packing_max))

    meta["limits_before_round"] = {
        k: (round(v[0], 3), round(v[1], 3)) for k, v in limits.items()
    }
    meta["explicit_limits"] = sorted(explicit)
    limits = _round_limits_to_grid(limits, skip=explicit | skip_round)

    for col, (lo, hi) in limits.items():
        if not (hi > lo):
            raise ValueError(f"Invalid limits for {col}: min={lo}, max={hi} (need max > min).")

    meta["limits"] = {k: (round(v[0], 4), round(v[1], 4)) for k, v in limits.items()}
    return limits, meta


def occlusion_round_step(column: str) -> float:
    """Upper-edge reporting-grid step used for an occlusion band."""
    return _ROUND_STEP_BY_COLUMN.get(column, 5.0)


def occlusion_lower_round_step(column: str, band_min: float | None = None) -> float:
    """
    Lower-edge reporting-grid step used for an occlusion band.

    LOI and contact-residue fractions use a 1%-point lower grid only when the
    coarser 5-point grid would round a positive empirical minimum down to zero.
    Otherwise they retain the 5-point grid. Their upper edge, and both BSA/kDa
    edges, always use the coarser 5-unit grid.
    """
    if column in (
        "lattice_burial_fraction_percent",
        "lattice_contact_residue_fraction_percent",
    ):
        if band_min is None or 0.0 < band_min < 5.0:
            return 1.0
    return occlusion_round_step(column)


def volume_lower_round_step(column: str, band_min: float | None = None) -> float:
    """
    Lower-edge reporting-grid step for a volume band.

    Packing density uses a 1%-point lower grid when a 5-point grid would collapse
    a positive assembled-lattice minimum to zero (assembled lattices are never empty).
    Matthews and solvent keep their usual steps at both edges.
    """
    if column == "packing_density_percent":
        if band_min is None or 0.0 < band_min < 5.0:
            return 1.0
    return _ROUND_STEP_BY_COLUMN.get(column, 5.0)


def profile_filename(profile: str, version: int) -> str:
    """File name for a named profile, for example slayer-compact v1 -> slayer_compact_v1.json."""
    slug = re.sub(r"[^a-z0-9]+", "_", str(profile).strip().lower()).strip("_")
    return f"{slug}_v{int(version)}.json"


def load_occlusion_profile(
    profile: str,
    *,
    version: int | None = None,
    profile_dir: str | None = None,
) -> dict[str, object]:
    """Load a shipped occlusion profile JSON written by ``calibrate_occlusion_profile.py``."""
    import json

    name = _OCCLUSION_PROFILE_FOR_MODE.get(profile, profile)
    ver = version if version is not None else _OCCLUSION_PROFILE_VERSION.get(name)
    if ver is None:
        raise ValueError(f"No shipped version registered for occlusion profile {profile!r}.")
    directory = profile_dir or PROFILE_DIR
    path = os.path.join(directory, profile_filename(name, ver))
    if not os.path.isfile(path):
        raise ValueError(
            f"Occlusion profile not found: {path}. Regenerate it with "
            "metrics/calibrate_occlusion_profile.py."
        )
    with open(path, encoding="utf-8") as f:
        doc = json.load(f)
    metric_limits = doc.get("metric_limits")
    if not isinstance(metric_limits, dict):
        raise ValueError(f"Malformed occlusion profile (no metric_limits): {path}")
    for col in OCCLUSION_COLUMNS:
        entry = metric_limits.get(col)
        if not isinstance(entry, dict) or "limits" not in entry:
            raise ValueError(f"Occlusion profile {path} is missing limits for {col}.")
    doc["_path"] = path
    return doc


def load_volume_profile(
    profile: str,
    *,
    version: int | None = None,
    profile_dir: str | None = None,
) -> dict[str, object]:
    """Load a shipped volume profile JSON (``{name}_volume_vN.json``)."""
    import json

    ver = version if version is not None else _VOLUME_PROFILE_VERSION.get(profile)
    if ver is None:
        raise ValueError(f"No shipped version registered for volume profile {profile!r}.")
    directory = profile_dir or PROFILE_DIR
    path = os.path.join(directory, profile_filename(f"{profile}_volume", ver))
    if not os.path.isfile(path):
        raise ValueError(f"Volume profile not found: {path}.")
    with open(path, encoding="utf-8") as f:
        doc = json.load(f)
    metric_limits = doc.get("metric_limits")
    if not isinstance(metric_limits, dict):
        raise ValueError(f"Malformed volume profile (no metric_limits): {path}")
    for col in _VOLUME_COLUMNS:
        entry = metric_limits.get(col)
        if not isinstance(entry, dict) or "limits" not in entry:
            raise ValueError(f"Volume profile {path} is missing limits for {col}.")
    doc["_path"] = path
    return doc


def empirical_band(
    values: list[float],
    *,
    close_span: float = EMPIRICAL_CLOSE_SPAN,
    wide_step: float = 5.0,
    lower_step: float | None = None,
    clamp_nonnegative: bool = True,
) -> tuple[float, float, str]:
    """
    Derive a within-run scoring band from observed values.

    Close cohorts (max − min ≤ ``close_span``): expand to the previous integer
    below the minimum and the next integer above the maximum
    (``floor(min) − 1`` … ``ceil(max) + 1``).

    Larger spreads: round outward onto ``wide_step`` (0/5 for %, 0.5 for Matthews).
    When ``lower_step`` differs from ``wide_step``, the lower edge uses that finer
    grid so a positive minimum is not collapsed to zero.

    Returns ``(lo, hi, rule)`` where ``rule`` is ``"close_integer_pad"`` or
    ``"wide_grid"``.
    """
    vals = [v for v in values if math.isfinite(v)]
    if not vals:
        raise ValueError("Cannot derive an empirical band from empty values.")
    obs_lo, obs_hi = min(vals), max(vals)
    span = obs_hi - obs_lo
    if span <= close_span:
        lo = float(math.floor(obs_lo) - 1)
        hi = float(math.ceil(obs_hi) + 1)
        rule = "close_integer_pad"
        pad = 1.0
    else:
        lo_step = lower_step if lower_step is not None else wide_step
        if lo_step != wide_step:
            lo, hi = _round_outward_asymmetric(
                obs_lo, obs_hi, lower_step=lo_step, upper_step=wide_step
            )
        else:
            lo, hi = _round_outward(obs_lo, obs_hi, step=wide_step)
        rule = "wide_grid"
        pad = wide_step
    if clamp_nonnegative and lo < 0.0:
        lo = 0.0
    if hi <= lo:
        hi = lo + pad
    return lo, hi, rule


def empirical_occlusion_band(
    values: list[float],
    *,
    column: str | None = None,
    close_span: float = EMPIRICAL_CLOSE_SPAN,
    clamp_nonnegative: bool = True,
) -> tuple[float, float, str]:
    """
    Within-run occlusion band.

    Close cohorts use ±1 integer padding (same as ``empirical_band``). Wide
    cohorts round outward onto the occlusion reporting grid; LOI/contact use an
    asymmetric lower step so a positive observed minimum is not collapsed to 0.
    """
    vals = [v for v in values if math.isfinite(v)]
    if not vals:
        raise ValueError("Cannot derive an empirical occlusion band from empty values.")
    obs_lo, obs_hi = min(vals), max(vals)
    span = obs_hi - obs_lo
    if span <= close_span:
        lo = float(math.floor(obs_lo) - 1)
        hi = float(math.ceil(obs_hi) + 1)
        rule = "close_integer_pad"
        pad = 1.0
    else:
        upper_step = 5.0
        lower_step = (
            occlusion_lower_round_step(column, obs_lo) if column else upper_step
        )
        lo, hi = _round_outward_asymmetric(
            obs_lo, obs_hi, lower_step=lower_step, upper_step=upper_step
        )
        rule = "wide_0_or_5"
        pad = upper_step
    if clamp_nonnegative and lo < 0.0:
        lo = 0.0
    if hi <= lo:
        hi = lo + pad
    return lo, hi, rule


def resolve_occlusion_limits(
    *,
    occlusion_scale_mode: str,
    occlusion_values: dict[str, list[float]] | None = None,
    loi_min: float | None = None,
    loi_max: float | None = None,
    contact_min: float | None = None,
    contact_max: float | None = None,
    bsa_kda_min: float | None = None,
    bsa_kda_max: float | None = None,
    profile_dir: str | None = None,
) -> tuple[dict[str, tuple[float, float]], dict[str, object]]:
    """
    Resolve (min, max) raw-value ranges for the three occlusion metrics.

    occlusion_scale_mode (shared vocabulary with volume):
      - cryst1: no transferable crystal LOI band → within-run empirical fallback
      - bbox: fixed multi-copy expanded-lattice occlusion profile
      - slayer-compact: tighter SlpA-core occlusion profile
      - empirical: within-run ±1 pad or 0/5 grid
      - cohort: legacy raw within-run min-max (no padding)
      - user: requires all three explicit --loi/--contact/--bsa-kda pairs
    """
    mode = _normalize_scale_mode(occlusion_scale_mode, default="bbox")
    if mode not in OCCLUSION_SCALE_MODES:
        raise ValueError(
            f"Unknown occlusion scale mode {occlusion_scale_mode!r}; "
            f"choose from {OCCLUSION_SCALE_MODES}."
        )

    overrides: dict[str, tuple[float, float]] = {}
    for col, lo, hi, flag in (
        ("lattice_burial_fraction_percent", loi_min, loi_max, "--loi"),
        ("lattice_contact_residue_fraction_percent", contact_min, contact_max, "--contact"),
        ("reference_chain_BSA_per_kDa_reference_chain_A2", bsa_kda_min, bsa_kda_max, "--bsa-kda"),
    ):
        if (lo is None) != (hi is None):
            raise ValueError(f"{flag}-min and {flag}-max must be given together.")
        if lo is not None and hi is not None:
            if float(hi) <= float(lo):
                raise ValueError(f"{flag}-max must exceed {flag}-min (got {lo}, {hi}).")
            overrides[col] = (float(lo), float(hi))

    meta: dict[str, object] = {"occlusion_scale_mode": mode}
    limits: dict[str, tuple[float, float]] = {}

    if mode == "user":
        missing = [c for c in OCCLUSION_COLUMNS if c not in overrides]
        if missing:
            raise ValueError(
                "--occlusion-scale-mode user requires --loi-min/--loi-max, "
                "--contact-min/--contact-max, and --bsa-kda-min/--bsa-kda-max."
            )
    elif mode in ("empirical", "cryst1"):
        if mode == "cryst1":
            meta["occlusion_note"] = (
                "cryst1 has no transferable LOI/contact/BSA band; "
                "using within-run empirical fallback"
            )
        values = occlusion_values or {}
        per_metric: dict[str, object] = {}
        for col in OCCLUSION_COLUMNS:
            vals = [v for v in values.get(col, []) if math.isfinite(v)]
            if not vals:
                raise ValueError(
                    f"--occlusion-scale-mode {mode} needs usable values for {col}."
                )
            lo, hi, rule = empirical_occlusion_band(vals, column=col)
            limits[col] = (lo, hi)
            per_metric[col] = {
                "observed_min": round(min(vals), 4),
                "observed_max": round(max(vals), 4),
                "span": round(max(vals) - min(vals), 4),
                "rule": rule,
                "limits": [lo, hi],
                "n": len(vals),
            }
        meta["empirical_bands"] = per_metric
        meta["empirical_close_span"] = EMPIRICAL_CLOSE_SPAN
    elif mode == "cohort":
        meta["occlusion_note"] = (
            "within-run cohort min-max (legacy; prefer --occlusion-scale-mode empirical)"
        )
    else:
        # bbox / slayer-compact → shipped occlusion profile
        profile_name = _OCCLUSION_PROFILE_FOR_MODE.get(mode, mode)
        doc = load_occlusion_profile(profile_name, profile_dir=profile_dir)
        metric_limits = doc["metric_limits"]
        assert isinstance(metric_limits, dict)
        for col in OCCLUSION_COLUMNS:
            lo, hi = metric_limits[col]["limits"]
            limits[col] = (float(lo), float(hi))
        meta["occlusion_profile"] = {
            "profile": doc.get("profile"),
            "requested_mode": mode,
            "version": doc.get("version"),
            "path": doc.get("_path"),
            "population": doc.get("population", ""),
            "band_method": doc.get("band_method", ""),
            "n_structures": doc.get("n_structures"),
            "sources": doc.get("sources", []),
            "observed": {
                col: [
                    metric_limits[col].get("observed_min"),
                    metric_limits[col].get("observed_max"),
                ]
                for col in OCCLUSION_COLUMNS
            },
            "notes": doc.get("notes", []),
        }

    limits.update(overrides)
    meta["explicit_occlusion_limits"] = sorted(overrides)
    meta["occlusion_limits"] = {k: (round(v[0], 4), round(v[1], 4)) for k, v in limits.items()}
    if mode == "cohort" and not limits:
        meta["occlusion_note"] = (
            "within-run cohort min-max (legacy; prefer --occlusion-scale-mode empirical)"
        )
    return limits, meta


def _norm_label(value: object) -> str:
    """Fold a column header / row label to a comparison key (unit glyphs and case ignored)."""
    s = str(value or "").strip().lower()
    for src, dst in (("å", "a"), ("³", "3"), ("²", "2"), ("µ", "u"), ("μ", "u")):
        s = s.replace(src, dst)
    return re.sub(r"[^a-z0-9]+", "", s)


def _build_label_index() -> dict[str, int]:
    index: dict[str, int] = {}
    for i, m in enumerate(METRICS):
        for label in (m.column, m.display, m.short, *m.aliases):
            key = _norm_label(label)
            if key:
                index.setdefault(key, i)
    return index


_LABEL_TO_METRIC: dict[str, int] = _build_label_index()


def _to_float(value: object) -> float:
    if value is None:
        return math.nan
    s = str(value).strip().replace("%", "").replace("\u2009", "").replace("\u00a0", "")
    if not s:
        return math.nan
    try:
        return float(s)
    except ValueError:
        return math.nan


def _parse_csv_list(raw: str | None) -> list[str]:
    if not raw:
        return []
    return [tok.strip() for tok in str(raw).split(",") if tok.strip()]


def _median(values: list[float]) -> float:
    vals = sorted(v for v in values if not math.isnan(v))
    if not vals:
        return math.nan
    n = len(vals)
    mid = n // 2
    if n % 2:
        return vals[mid]
    return 0.5 * (vals[mid - 1] + vals[mid])


def _mean(values: list[float]) -> float:
    vals = [v for v in values if not math.isnan(v)]
    return sum(vals) / len(vals) if vals else math.nan


def _stdev(values: list[float]) -> float:
    """Sample standard deviation (n-1); NaN when fewer than two finite values."""
    vals = [v for v in values if not math.isnan(v)]
    if len(vals) < 2:
        return math.nan
    mu = sum(vals) / len(vals)
    return math.sqrt(sum((v - mu) ** 2 for v in vals) / (len(vals) - 1))


def _mad(values: list[float]) -> float:
    """Median absolute deviation about the median, rescaled to an SD equivalent."""
    vals = [v for v in values if not math.isnan(v)]
    if not vals:
        return math.nan
    med = _median(vals)
    return _MAD_TO_SD * _median([abs(v - med) for v in vals])


def _missing_metrics_error(found: set[int], path: str, *, orientation: str) -> None:
    missing = [m.display for i, m in enumerate(METRICS) if i not in found]
    if missing:
        raise ValueError(
            f"Input CSV ({orientation}) missing required metric(s): "
            + ", ".join(missing)
            + f" ({path})"
        )


def _rows_from_wide(raw: list[list[str]], path: str) -> list[dict[str, str]]:
    """One row per structure: 'structure_stem' plus a column per metric."""
    header = raw[0]
    by_key = {_norm_label(h): pos for pos, h in enumerate(header)}
    stem_pos = by_key.get(_norm_label("structure_stem"))
    if stem_pos is None:
        raise ValueError(f"Input CSV missing required column 'structure_stem': {path}")

    positions: dict[int, int] = {}
    for pos, h in enumerate(header):
        idx = _LABEL_TO_METRIC.get(_norm_label(h))
        if idx is not None and idx not in positions:
            positions[idx] = pos
    _missing_metrics_error(set(positions), path, orientation="structures as rows")

    out: list[dict[str, str]] = []
    for row in raw[1:]:
        if not any(str(c).strip() for c in row):
            continue
        stem = str(row[stem_pos]).strip() if stem_pos < len(row) else ""
        if not stem:
            continue
        rec: dict[str, str] = {"structure_stem": stem}
        for idx, pos in positions.items():
            rec[METRICS[idx].column] = str(row[pos]).strip() if pos < len(row) else ""
        out.append(rec)
    return out


def _rows_from_transposed(raw: list[list[str]], path: str) -> list[dict[str, str]]:
    """Metrics as rows, structures as columns (first column holds the metric label)."""
    header = raw[0]
    stems = [str(c).strip() for c in header[1:]]

    values_by_metric: dict[int, list[str]] = {}
    for row in raw[1:]:
        if not row:
            continue
        idx = _LABEL_TO_METRIC.get(_norm_label(row[0]))
        if idx is None or idx in values_by_metric:
            continue
        values_by_metric[idx] = [str(c).strip() for c in row[1:]]
    _missing_metrics_error(set(values_by_metric), path, orientation="metrics as rows")

    out: list[dict[str, str]] = []
    for col, stem in enumerate(stems):
        if not stem:
            continue
        rec: dict[str, str] = {"structure_stem": stem}
        for idx, vals in values_by_metric.items():
            rec[METRICS[idx].column] = vals[col] if col < len(vals) else ""
        out.append(rec)
    return out


def detect_orientation(raw: list[list[str]]) -> str:
    """Return 'wide' (structures as rows) or 'transposed' (metrics as rows)."""
    header = raw[0] if raw else []
    header_keys = {_norm_label(h) for h in header}
    if _norm_label("structure_stem") in header_keys:
        return "wide"
    if any(_norm_label(h) in _LABEL_TO_METRIC for h in header[1:]):
        return "wide"
    first_col = [row[0] for row in raw[1:] if row]
    if any(_norm_label(c) in _LABEL_TO_METRIC for c in first_col):
        return "transposed"
    return "unknown"


def load_rows(path: str) -> list[dict[str, str]]:
    """
    Read a metrics table in either orientation and return one dict per structure.

    Accepted layouts (auto-detected):
      - structures as rows: ``combined_lattice_vs_ec.csv`` from ``lattice_compare_batch.py``
        (needs ``structure_stem`` plus the six metric columns);
      - metrics as rows: a transposed summary table whose first column holds metric
        labels (for example ``Packing density (%)``) and whose header lists structure names.

    Metric labels are matched case-insensitively, ignoring unit glyphs and punctuation,
    so both snake_case column names and display labels work.
    """
    if not os.path.isfile(path):
        raise ValueError(f"Input CSV not found: {path}")
    with open(path, newline="", encoding="utf-8") as f:
        raw = [list(row) for row in csv.reader(f)]
    raw = [row for row in raw if any(str(c).strip() for c in row)]
    if not raw:
        raise ValueError(f"Input CSV is empty: {path}")

    orientation = detect_orientation(raw)
    if orientation == "wide":
        return _rows_from_wide(raw, path)
    if orientation == "transposed":
        return _rows_from_transposed(raw, path)
    raise ValueError(
        "Could not recognise the metrics table layout in "
        f"{path}: expected 'structure_stem' plus metric columns (structures as rows), "
        "or metric labels in the first column (metrics as rows)."
    )


def _order_stems(
    stems: list[str],
    *,
    subset: list[str],
    order: list[str],
    order_file_stems: list[str],
) -> list[str]:
    present = list(dict.fromkeys(stems))
    if subset:
        keep = set(subset)
        unknown = [s for s in subset if s not in present]
        if unknown:
            raise ValueError(f"--stems not found in input: {unknown}")
        present = [s for s in present if s in keep]
    wanted = list(dict.fromkeys([*order, *order_file_stems]))
    if not wanted:
        return present
    unknown = [s for s in wanted if s not in present]
    if unknown:
        raise ValueError(f"Structure order entries not found in input: {unknown}")
    lead = [s for s in wanted if s in present]
    rest = [s for s in present if s not in set(lead)]
    return [*lead, *rest]


def _load_order_file(path: str | None) -> list[str]:
    if not path:
        return []
    p = os.path.expanduser(str(path).strip())
    if not os.path.isfile(p):
        raise ValueError(f"Structure order file not found: {p}")
    out: list[str] = []
    with open(p, encoding="utf-8", errors="replace") as f:
        for line in f:
            raw = line.strip()
            if not raw or raw.startswith("#"):
                continue
            out.append(raw.split(",")[0].strip())
    return out


def _validate_score_range(score_range: tuple[float, float]) -> tuple[float, float]:
    """Return a finite, increasing score range or raise a CLI-friendly error."""
    try:
        score_min, score_max = (float(score_range[0]), float(score_range[1]))
    except (TypeError, ValueError, IndexError) as exc:
        raise ValueError("--score-range requires exactly two numeric values: MIN MAX.") from exc
    if not (math.isfinite(score_min) and math.isfinite(score_max)):
        raise ValueError("--score-range values must be finite.")
    if score_max <= score_min:
        raise ValueError(
            f"--score-range MAX must exceed MIN (got {score_min:g} {score_max:g})."
        )
    return score_min, score_max


def score_metric(
    x: float,
    xmin: float,
    xmax: float,
    invert: bool,
    *,
    score_range: tuple[float, float] = DEFAULT_SCORE_RANGE,
) -> float:
    """
    Map a raw value onto ``score_range``; equal cohort -> midpoint; NaN -> NaN.

    Fixed bands can be exceeded by an individual structure, so the score is clamped
    to the plotted ring. ``clipped_values`` reports which values were outside a band.
    """
    score_min, score_max = _validate_score_range(score_range)
    if math.isnan(x):
        return math.nan
    if not (xmax > xmin):
        return (score_min + score_max) / 2.0
    span = score_max - score_min
    raw = score_min + span * (x - xmin) / (xmax - xmin)
    score = score_max - (raw - score_min) if invert else raw
    return min(max(score, score_min), score_max)


def clipped_values(
    kept: list[str],
    raw_by_stem: dict[str, list[float]],
    used_limits: dict[str, tuple[float, float]],
    *,
    metrics: tuple[MetricSpec, ...] = METRICS,
) -> list[dict[str, object]]:
    """List raw values outside their raw-value band (score saturates at a range endpoint)."""
    out: list[dict[str, object]] = []
    for i, m in enumerate(metrics):
        lo, hi = used_limits.get(m.column, (math.nan, math.nan))
        if not (math.isfinite(lo) and math.isfinite(hi) and hi > lo):
            continue
        for stem in kept:
            val = raw_by_stem[stem][i]
            if math.isnan(val) or lo <= val <= hi:
                continue
            out.append(
                {
                    "structure_stem": stem,
                    "metric": m.column,
                    "value": round(val, 4),
                    "limits": [lo, hi],
                    "side": "below" if val < lo else "above",
                }
            )
    return out


def compute_scores(
    rows: list[dict[str, str]],
    stems: list[str],
    *,
    volume_limits: dict[str, tuple[float, float]] | None = None,
    occlusion_limits: dict[str, tuple[float, float]] | None = None,
    score_range: tuple[float, float] = DEFAULT_SCORE_RANGE,
    metrics: tuple[MetricSpec, ...] = METRICS,
) -> tuple[dict[str, list[float]], dict[str, list[float]], list[str], dict[str, tuple[float, float]]]:
    """
    Return (scores_by_stem, raw_by_stem, kept_stems, used_limits).

    Volume metrics (packing / Matthews / solvent) use ``volume_limits`` and the
    occlusion metrics (BSA/kDa, LOI, contact residue %) use ``occlusion_limits``
    when supplied; any axis without a fixed band falls back to cohort min–max.
    Spoke order follows ``metrics``.
    """
    score_range = _validate_score_range(score_range)
    by_stem = {str(r.get("structure_stem", "")).strip(): r for r in rows}
    raw_by_stem: dict[str, list[float]] = {}
    for stem in stems:
        row = by_stem.get(stem, {})
        raw_by_stem[stem] = [_to_float(row.get(m.column)) for m in metrics]

    kept = [s for s in stems if not all(math.isnan(v) for v in raw_by_stem[s])]
    dropped = [s for s in stems if s not in set(kept)]
    if dropped:
        print(f"Warning: skipping structures with no usable metrics: {dropped}", file=sys.stderr)

    fixed_limits: dict[str, tuple[float, float]] = {}
    fixed_limits.update(volume_limits or {})
    fixed_limits.update(occlusion_limits or {})

    used_limits: dict[str, tuple[float, float]] = {}
    mins: list[float] = []
    maxs: list[float] = []
    for i, m in enumerate(metrics):
        col_vals = [raw_by_stem[s][i] for s in kept if not math.isnan(raw_by_stem[s][i])]
        if m.column in fixed_limits:
            lo, hi = fixed_limits[m.column]
            used_limits[m.column] = (float(lo), float(hi))
            mins.append(float(lo))
            maxs.append(float(hi))
        elif col_vals:
            lo, hi = min(col_vals), max(col_vals)
            used_limits[m.column] = (lo, hi)
            mins.append(lo)
            maxs.append(hi)
        else:
            used_limits[m.column] = (math.nan, math.nan)
            mins.append(math.nan)
            maxs.append(math.nan)

    # Cohort min-max is invariant to a constant factor, but a fixed percent band is not:
    # warn rather than silently score a fraction-valued column against a 0-100 band.
    for i, m in enumerate(metrics):
        if m.column not in fixed_limits or m.column not in _PERCENT_COLUMNS:
            continue
        vals = [raw_by_stem[s][i] for s in kept if not math.isnan(raw_by_stem[s][i])]
        if vals and 0.0 < max(vals) <= 1.0 and fixed_limits[m.column][1] > 1.0:
            print(
                f"Warning: {m.display} values are all <= 1 but the scoring band is "
                f"{fixed_limits[m.column][0]:g}-{fixed_limits[m.column][1]:g}; the input "
                "column may be a fraction where a percentage is expected.",
                file=sys.stderr,
            )

    scores_by_stem: dict[str, list[float]] = {}
    for stem in kept:
        scores_by_stem[stem] = [
            score_metric(
                raw_by_stem[stem][i],
                mins[i],
                maxs[i],
                metrics[i].invert,
                score_range=score_range,
            )
            for i in range(len(metrics))
        ]
    return scores_by_stem, raw_by_stem, kept, used_limits


def _mean_score(scores: list[float]) -> float:
    vals = [v for v in scores if not math.isnan(v)]
    return sum(vals) / len(vals) if vals else math.nan


def compute_deviations(
    kept: list[str],
    raw_by_stem: dict[str, list[float]],
    *,
    deviation_scale: str = "sd",
    metrics: tuple[MetricSpec, ...] = METRICS,
) -> tuple[dict[str, list[float]], dict[str, object]]:
    """
    Standardise every axis about the cohort centre, returning (dev_by_stem, meta).

    Deviations are computed from the *raw* values rather than from band scores, so
    a cohort occupying a narrow slice of a wide physical band keeps its full
    resolution instead of collapsing (or clipping) against the band edges. The
    sign convention matches score mode: inverted metrics are negated so that a
    larger radius still means a tighter / more occluded lattice.

    ``deviation_scale`` selects the spread estimator:
      - ``sd``    : (x - mean) / sample SD, the classic z-score.
      - ``mad``   : (x - median) / (1.4826 * MAD), robust to a single outlier.
      - ``range`` : (x - midrange) / half-range, mapping the cohort exactly onto
                    [-1, +1] with no distributional assumption.
    """
    if deviation_scale not in DEVIATION_SCALES:
        raise ValueError(
            f"Unsupported --deviation-scale {deviation_scale!r}; choose from {DEVIATION_SCALES}."
        )

    per_axis: list[dict[str, object]] = []
    dev_by_stem: dict[str, list[float]] = {s: [] for s in kept}
    degenerate: list[str] = []

    for i, m in enumerate(metrics):
        vals = [raw_by_stem[s][i] for s in kept if not math.isnan(raw_by_stem[s][i])]
        if deviation_scale == "sd":
            centre, spread = _mean(vals), _stdev(vals)
        elif deviation_scale == "mad":
            centre, spread = _median(vals), _mad(vals)
            # A tied-heavy small cohort can give MAD == 0 while the values still
            # differ; fall back to the SD so the axis does not silently flatten.
            if vals and (math.isnan(spread) or spread <= 0.0):
                spread = _stdev(vals)
        else:
            centre = 0.5 * (min(vals) + max(vals)) if vals else math.nan
            spread = 0.5 * (max(vals) - min(vals)) if vals else math.nan

        usable = bool(vals) and math.isfinite(spread) and spread > 0.0
        if vals and not usable:
            degenerate.append(m.column)
        sign = -1.0 if m.invert else 1.0
        for stem in kept:
            x = raw_by_stem[stem][i]
            if math.isnan(x) or not usable:
                dev_by_stem[stem].append(math.nan if math.isnan(x) else 0.0)
            else:
                dev_by_stem[stem].append(sign * (x - centre) / spread)
        per_axis.append(
            {
                "metric": m.column,
                "centre": None if math.isnan(centre) else round(centre, 4),
                "spread": None if (math.isnan(spread) or not usable) else round(spread, 4),
                "inverted": m.invert,
            }
        )

    meta: dict[str, object] = {
        "display_mode": "deviation",
        "deviation_scale": deviation_scale,
        "deviation_axes": per_axis,
    }
    if degenerate:
        meta["deviation_degenerate_axes"] = degenerate
        print(
            "Warning: zero spread on "
            + ", ".join(degenerate)
            + "; those axes are drawn at the cohort centre.",
            file=sys.stderr,
        )
    return dev_by_stem, meta


def deviation_limit(
    kept: list[str],
    dev_by_stem: dict[str, list[float]],
    *,
    requested: float = DEFAULT_DEVIATION_LIMIT,
) -> float:
    """Symmetric radial limit for deviation mode; 0 auto-fits to the observed spread."""
    if requested and requested > 0:
        return float(requested)
    finite = [
        abs(v)
        for s in kept
        for v in dev_by_stem[s]
        if not math.isnan(v)
    ]
    peak = max(finite) if finite else 0.0
    # Round up to a half-unit so the rings land on readable values, and keep a
    # little headroom so the outermost vertex is not drawn on the frame itself.
    return max(1.0, math.ceil((peak * 1.08) * 2.0) / 2.0)


def _ghost_values(
    values_by_stem: dict[str, list[float]],
    kept: list[str],
    reference_stem: str | None,
    *,
    ghost_mode: str = "loo-median",
    metrics: tuple[MetricSpec, ...] = METRICS,
) -> tuple[dict[str, list[float]], str]:
    """
    Per-panel ghost silhouettes, as (ghost_by_stem, label).

    The leave-one-out modes build each panel's reference from the other n-1
    structures. That keeps the ghost off the polygon it is being compared with
    (an odd cohort's median is always one of its own members) and turns the gap
    between the two outlines into a genuine "this structure vs. its peers" read.
    """
    n_axes = len(metrics)
    if reference_stem:
        if reference_stem not in values_by_stem:
            raise ValueError(f"--reference-stem not found among plotted structures: {reference_stem}")
        ref = values_by_stem[reference_stem]
        return {s: ref for s in kept}, f"reference: {reference_stem}"
    if ghost_mode not in GHOST_MODES:
        raise ValueError(f"Unsupported --ghost {ghost_mode!r}; choose from {GHOST_MODES}.")

    centre = _mean if ghost_mode.endswith("mean") else _median
    stat_label = "mean" if ghost_mode.endswith("mean") else "median"

    if ghost_mode.startswith("loo"):
        # A single-structure cohort has no "others" to average, so fall back to
        # the cohort statistic (which for n=1 is just that structure).
        if len(kept) < 2:
            shared = [centre([values_by_stem[s][i] for s in kept]) for i in range(n_axes)]
            return {s: shared for s in kept}, f"cohort {stat_label}"
        ghost_by_stem = {
            stem: [
                centre([values_by_stem[s][i] for s in kept if s != stem])
                for i in range(n_axes)
            ]
            for stem in kept
        }
        return ghost_by_stem, f"leave-one-out {stat_label} (n={len(kept) - 1} others)"

    shared = [centre([values_by_stem[s][i] for s in kept]) for i in range(n_axes)]
    return {s: shared for s in kept}, f"cohort {stat_label}"


def _write_scores_csv(
    path: str,
    kept: list[str],
    values_by_stem: dict[str, list[float]],
    raw_by_stem: dict[str, list[float]],
    *,
    value_suffix: str = "score",
    metrics: tuple[MetricSpec, ...] = METRICS,
) -> None:
    header = ["structure_stem"]
    for m in metrics:
        header.append(f"{m.column}__{value_suffix}")
    for m in metrics:
        header.append(f"{m.column}__raw")
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(header)
        for stem in kept:
            values = values_by_stem[stem]
            raws = raw_by_stem[stem]
            row: list[object] = [stem]
            row.extend("" if math.isnan(v) else round(v, 3) for v in values)
            row.extend("" if math.isnan(v) else round(v, 4) for v in raws)
            w.writerow(row)


def _closed(seq: list[float]) -> list[float]:
    return [*seq, seq[0]]


def _format_score(value: float) -> str:
    """Compact score label that preserves useful non-integer tick values."""
    return f"{value:.2f}".rstrip("0").rstrip(".")


def _score_ticks(score_range: tuple[float, float]) -> list[float]:
    """Four evenly spaced rings (25%, 50%, 75%, and the upper endpoint)."""
    score_min, score_max = score_range
    span = score_max - score_min
    return [score_min + span * fraction for fraction in (0.25, 0.5, 0.75, 1.0)]


# Polygon edge encoding for adjacent spoke disagreement (shared across the run).
# Width/alpha map |Δplot| / (smax−smin) onto a fixed interval so panels are comparable.
# A sub-linear power expands mid-range gaps (typical |Δ| ≪ full span) so zig-zags read clearly.
EDGE_DELTA_LW_MIN = 0.7
EDGE_DELTA_LW_MAX = 6.5
EDGE_DELTA_ALPHA_MIN = 0.18
EDGE_DELTA_ALPHA_MAX = 1.00
EDGE_DELTA_POWER = 0.55  # visual = (|Δ|/span) ** power; <1 boosts moderate disagreements
DEFAULT_POLYGON_LW = 1.8


def edge_delta_styles(
    values: list[float],
    score_range: tuple[float, float],
    *,
    lw_min: float = EDGE_DELTA_LW_MIN,
    lw_max: float = EDGE_DELTA_LW_MAX,
    alpha_min: float = EDGE_DELTA_ALPHA_MIN,
    alpha_max: float = EDGE_DELTA_ALPHA_MAX,
    power: float = EDGE_DELTA_POWER,
) -> list[tuple[float, float]]:
    """
    Per-edge (linewidth, alpha) from adjacent plot-value gaps.

    Edge *i* joins spokes *i* and *i+1* (wrapping). Fraction is
    ``(|v_i − v_{i+1}| / (smax − smin)) ** power``, clamped to [0, 1], using the
    same score/deviation range as the radar rings so the mapping is run-wide and
    comparable across panels. ``power < 1`` expands mid-range gaps so typical
    zig-zags (well below the full span) remain visible. Returns one
    ``(lw, alpha)`` per spoke (the edge leaving that spoke clockwise).
    """
    score_min, score_max = score_range
    span = score_max - score_min
    n = len(values)
    if n == 0:
        return []
    if not (math.isfinite(span) and span > 0):
        return [(lw_min, alpha_min)] * n
    if not (math.isfinite(power) and power > 0):
        power = 1.0
    out: list[tuple[float, float]] = []
    for i in range(n):
        a = values[i]
        b = values[(i + 1) % n]
        if math.isnan(a) or math.isnan(b):
            out.append((lw_min, alpha_min))
            continue
        frac = min(1.0, abs(a - b) / span)
        t = frac**power
        lw = lw_min + (lw_max - lw_min) * t
        alpha = alpha_min + (alpha_max - alpha_min) * t
        out.append((lw, alpha))
    return out


def _draw_panel(
    ax,
    angles,
    values,
    raw,
    ghost,
    color,
    title,
    annotate_values,
    score_range: tuple[float, float],
    *,
    display_mode: str = "score",
    metrics: tuple[MetricSpec, ...] = METRICS,
    edge_delta: bool = False,
) -> None:
    score_min, score_max = score_range
    score_span = score_max - score_min
    # Matplotlib polar axes cannot host a negative radius, so signed displays
    # (deviation mode) are shifted onto [0, span] while the tick labels keep
    # the original signed values.
    offset = -score_min if score_min < 0.0 else 0.0
    r_min, r_max = score_min + offset, score_max + offset

    def _to_radius(v: float) -> float:
        if math.isnan(v):
            return r_min
        return min(max(v + offset, r_min), r_max)

    plot_values = [_to_radius(v) for v in values]
    ax.set_theta_offset(math.pi / 2)
    ax.set_theta_direction(-1)
    ax.set_xticks(angles)
    ax.set_xticklabels([m.short for m in metrics], fontsize=7)
    ax.set_ylim(r_min, r_max)
    ticks = _score_ticks(score_range)
    ax.set_yticks([t + offset for t in ticks])
    ax.set_yticklabels([_format_score(v) for v in ticks], fontsize=6, color="0.4")
    ax.set_rlabel_position(180.0 / len(metrics))
    ax.grid(True, color="0.85", linewidth=0.6)
    for spine in ax.spines.values():
        spine.set_edgecolor("0.7")

    # In deviation mode the cohort centre is the mid-ring; draw it more clearly so
    # "above / below peers" is readable at a glance.
    if display_mode == "deviation":
        mid = 0.5 * (r_min + r_max)
        ax.plot(
            _closed(list(angles)),
            [mid] * (len(metrics) + 1),
            color="0.75",
            linewidth=0.8,
            linestyle=":",
            zorder=1,
        )

    ang_c = _closed(list(angles))
    ghost_c = _closed([_to_radius(v) for v in ghost])
    ax.plot(ang_c, ghost_c, color="0.55", linewidth=1.0, linestyle="--", zorder=2)

    val_c = _closed(plot_values)
    ax.fill(ang_c, val_c, color=color, alpha=0.14 if edge_delta else 0.22, zorder=2)
    if edge_delta:
        # Variable edge weight from adjacent plot-value gaps; shared |Δ|/span scale.
        styles = edge_delta_styles(values, score_range)
        n = len(metrics)
        for i in range(n):
            j = (i + 1) % n
            lw, alpha = styles[i]
            ax.plot(
                [angles[i], angles[j]],
                [plot_values[i], plot_values[j]],
                color=color,
                linewidth=lw,
                alpha=alpha,
                solid_capstyle="round",
                solid_joinstyle="round",
                zorder=3,
            )
    else:
        ax.plot(ang_c, val_c, color=color, linewidth=DEFAULT_POLYGON_LW, zorder=3)

    if annotate_values:
        for ang, sc, rv in zip(angles, plot_values, raw):
            label = "n/a" if math.isnan(rv) else f"{rv:.1f}"
            # Push labels outward from the vertex; clamp to a minimum radius so
            # near-centre vertices on opposite spokes do not overlap.
            r_text = min(
                max(sc + 0.10 * score_span, r_min + 0.16 * score_span),
                r_max + 0.04 * score_span,
            )
            ax.annotate(
                label,
                xy=(ang, r_text),
                fontsize=5.5,
                color="0.12",
                ha="center",
                va="center",
                zorder=4,
            )
    ax.set_title(title, fontsize=9.0, pad=10)


def _render_grid_pages(
    kept: list[str],
    values_by_stem: dict[str, list[float]],
    raw_by_stem: dict[str, list[float]],
    ghost_by_stem: dict[str, list[float]],
    ghost_label: str,
    *,
    output_dir: str,
    out_format: str,
    dpi: int,
    title: str | None,
    annotate_values: bool,
    max_per_page: int,
    score_range: tuple[float, float],
    scale_note: str = "",
    tag: str = "",
    display_mode: str = "score",
    metrics: tuple[MetricSpec, ...] = METRICS,
    edge_delta: bool = False,
) -> list[str]:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    cmap = plt.get_cmap("tab10")
    angles = spoke_angles(metrics)
    if display_mode == "deviation":
        base_title = title or "Lattice packing / occlusion (cohort deviations)"
        footnote = (
            "Radius = signed deviation from cohort centre (larger = tighter / more occluded). "
            "Dotted mid-ring = cohort centre; dashed ghost = " + ghost_label + ". "
            "Spoke labels are raw metric values."
        )
    else:
        base_title = title or "Lattice packing / occlusion radar"
        footnote = (
            "Larger polygon = tighter / more occluded. Inverted axes: Matthews, Solvent. "
            "Dashed ghost = " + ghost_label + "."
        )
    if edge_delta:
        footnote += (
            " Edge width/opacity ∝ (|adjacent plot-value difference| / score-span)"
            f"^{EDGE_DELTA_POWER:g} (shared scale; emphasises zig-zags)."
        )
    if scale_note:
        footnote = scale_note + " " + footnote

    if max_per_page and max_per_page > 0:
        pages = [kept[i : i + max_per_page] for i in range(0, len(kept), max_per_page)]
    else:
        pages = [kept]

    written: list[str] = []
    multi = len(pages) > 1
    for page_idx, page_stems in enumerate(pages, start=1):
        n = len(page_stems)
        ncols = max(1, int(math.ceil(math.sqrt(n))))
        nrows = int(math.ceil(n / ncols))
        fig, _axes = plt.subplots(
            nrows,
            ncols,
            subplot_kw={"polar": True},
            figsize=(3.6 * ncols, 4.1 * nrows),
        )
        axes = list(_axes.flat) if hasattr(_axes, "flat") else [_axes]
        for idx, stem in enumerate(page_stems):
            color = cmap((kept.index(stem)) % 10)
            _draw_panel(
                axes[idx],
                angles,
                values_by_stem[stem],
                raw_by_stem[stem],
                ghost_by_stem[stem],
                color,
                stem,
                annotate_values,
                score_range,
                display_mode=display_mode,
                metrics=metrics,
                edge_delta=edge_delta,
            )
        for ax in axes[n:]:
            ax.set_visible(False)

        page_title = base_title if not multi else f"{base_title} (page {page_idx}/{len(pages)})"
        # Reserve room under the suptitle for the top row's panel titles, and keep
        # enough inter-row space that a panel title never touches the row above.
        if nrows == 1:
            top, bottom, hspace = 0.86, 0.12, 0.30
        else:
            top, bottom, hspace = 0.90, 0.06, 0.48
        fig.subplots_adjust(top=top, bottom=bottom, left=0.06, right=0.94, hspace=hspace, wspace=0.45)
        fig.suptitle(page_title, fontsize=13, y=0.985)
        fig.text(0.5, 0.015, footnote, ha="center", va="bottom", fontsize=7.0, color="0.35")

        if multi:
            out_path = os.path.join(
                output_dir, _tagged_name("radar_grid", out_format, tag, page=page_idx)
            )
        else:
            out_path = os.path.join(output_dir, _tagged_name("radar_grid", out_format, tag))
        fig.savefig(out_path, dpi=dpi, bbox_inches="tight")
        plt.close(fig)
        written.append(out_path)
    return written


def _render_score_heatmap(
    kept: list[str],
    values_by_stem: dict[str, list[float]],
    *,
    output_dir: str,
    out_format: str,
    dpi: int,
    title: str | None,
    score_range: tuple[float, float],
    tag: str = "",
    display_mode: str = "score",
    deviation_scale: str = "sd",
    metrics: tuple[MetricSpec, ...] = METRICS,
) -> str:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np

    data = np.array(
        [[math.nan if math.isnan(v) else v for v in values_by_stem[s]] for s in kept],
        dtype=float,
    )
    fig_h = max(2.2, 0.42 * len(kept) + 1.6)
    fig_w = max(5.0, 0.95 * len(metrics) + 3.0)
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))
    score_min, score_max = score_range
    if display_mode == "deviation":
        cmap = plt.get_cmap("coolwarm").copy()
        cmap.set_bad("0.85")
        default_title = f"Lattice metric deviations ({deviation_scale}; ±{_format_score(score_max)})"
        cbar_label = (
            f"Signed deviation ({deviation_scale}; "
            f"0 = cohort centre, + = tighter / more occluded)"
        )
    else:
        cmap = plt.get_cmap("viridis").copy()
        cmap.set_bad("0.85")
        range_label = f"{_format_score(score_min)}–{_format_score(score_max)}"
        default_title = f"Lattice metric scores ({range_label})"
        cbar_label = (
            f"Score ({_format_score(score_min)} = loosest, "
            f"{_format_score(score_max)} = tightest)"
        )
    im = ax.imshow(
        np.ma.masked_invalid(data),
        cmap=cmap,
        vmin=score_min,
        vmax=score_max,
        aspect="auto",
    )
    ax.set_xticks(range(len(metrics)))
    ax.set_xticklabels([m.short.replace("\n", " ") for m in metrics], fontsize=7, rotation=30, ha="right")
    ax.set_yticks(range(len(kept)))
    ax.set_yticklabels(kept, fontsize=7)
    ax.set_title(title or default_title, fontsize=11, pad=10)

    for i in range(len(kept)):
        for j in range(len(metrics)):
            v = data[i, j]
            if math.isnan(v):
                continue
            relative = (v - score_min) / (score_max - score_min)
            if display_mode == "deviation":
                txt_color = "black" if 0.25 < relative < 0.75 else "white"
            else:
                txt_color = "white" if relative < 0.6 else "black"
            ax.text(j, i, f"{v:.1f}", ha="center", va="center", fontsize=6, color=txt_color)

    cbar = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    cbar.set_label(cbar_label, fontsize=8)
    cbar.ax.tick_params(labelsize=7)
    fig.tight_layout()
    out_path = os.path.join(output_dir, _tagged_name("radar_scores_heatmap", out_format, tag))
    fig.savefig(out_path, dpi=dpi, bbox_inches="tight")
    plt.close(fig)
    return out_path


def _occlusion_scale_note(
    occlusion_meta: dict[str, object],
    used_limits: dict[str, tuple[float, float]],
    *,
    n_clipped: int = 0,
    score_range: tuple[float, float] = DEFAULT_SCORE_RANGE,
) -> str:
    """One-line footnote describing how the occlusion axes were scaled."""
    mode = str(occlusion_meta.get("occlusion_scale_mode", "bbox"))
    overrides = occlusion_meta.get("explicit_occlusion_limits") or []
    bsa = used_limits.get("reference_chain_BSA_per_kDa_reference_chain_A2", (math.nan, math.nan))
    loi = used_limits.get("lattice_burial_fraction_percent", (math.nan, math.nan))
    con = used_limits.get("lattice_contact_residue_fraction_percent", (math.nan, math.nan))
    bands = (
        f"BSA/kDa {bsa[0]:g}–{bsa[1]:g} Å²; "
        f"LOI {loi[0]:g}–{loi[1]:g}%; "
        f"contact {con[0]:g}–{con[1]:g}%"
    )

    if mode == "cohort" and not overrides:
        return "Occlusion axes: cohort min–max (legacy within-run)."
    if mode == "user":
        note = f"Occlusion axes: user-defined ({bands})."
    elif mode == "cohort":
        note = f"Occlusion axes: cohort min–max with explicit overrides ({bands})."
    elif mode in ("empirical", "cryst1"):
        prefix = (
            "Occlusion axes: cryst1→empirical fallback "
            if mode == "cryst1"
            else "Occlusion axes: within-run empirical "
        )
        note = (
            f"{prefix}"
            f"(close ≤{occlusion_meta.get('empirical_close_span', EMPIRICAL_CLOSE_SPAN):g} "
            f"→ floor(min)−1…ceil(max)+1; wider → 0/5 grid; {bands})."
        )
        if overrides:
            note += " Overridden: " + ", ".join(str(c) for c in overrides) + "."
    else:
        profile = occlusion_meta.get("occlusion_profile") or {}
        assert isinstance(profile, dict)
        version = profile.get("version", "")
        n = profile.get("n_structures", "")
        note = (
            f"Occlusion axes: {mode} profile v{version} "
            f"(n={n} FoldKit multi-copy lattices; {bands})."
        )
        if overrides:
            note += " Overridden: " + ", ".join(str(c) for c in overrides) + "."
    if n_clipped:
        score_min, score_max = score_range
        note += (
            f" {n_clipped} value(s) clipped to "
            f"{_format_score(score_min)}/{_format_score(score_max)} "
            "(see radar_scale JSON)."
        )
    return note


def generate_radar(
    *,
    input_csv: str,
    output_dir: str,
    stems: list[str] | None = None,
    structure_order: list[str] | None = None,
    structure_order_file: str | None = None,
    reference_stem: str | None = None,
    ghost_mode: str = "loo-median",
    display_mode: str = "score",
    deviation_scale: str = "sd",
    deviation_limit_value: float = DEFAULT_DEVIATION_LIMIT,
    sort_by: str = "none",
    annotate_values: bool = True,
    out_format: str = "png",
    dpi: int = 200,
    title: str | None = None,
    max_per_page: int = 0,
    score_range: tuple[float, float] = DEFAULT_SCORE_RANGE,
    tag_outputs: bool = True,
    scale_mode: str = "bbox",
    matthews_min: float | None = None,
    matthews_max: float | None = None,
    solvent_min: float | None = None,
    solvent_max: float | None = None,
    packing_min: float | None = None,
    packing_max: float | None = None,
    occlusion_scale_mode: str = "bbox",
    loi_min: float | None = None,
    loi_max: float | None = None,
    contact_min: float | None = None,
    contact_max: float | None = None,
    bsa_kda_min: float | None = None,
    bsa_kda_max: float | None = None,
    metric_order: str | list[str] | tuple[str, ...] | None = None,
    write_heatmap: bool = True,
    edge_delta: bool = False,
) -> list[str]:
    """Build radar grid, optional score heatmap, and radar_scores.csv. Returns written paths."""
    if out_format not in _VALID_FORMATS:
        raise ValueError(f"Unsupported --format {out_format!r}; choose from {_VALID_FORMATS}.")
    if display_mode not in DISPLAY_MODES:
        raise ValueError(f"Unsupported --display {display_mode!r}; choose from {DISPLAY_MODES}.")
    if ghost_mode not in GHOST_MODES:
        raise ValueError(f"Unsupported --ghost {ghost_mode!r}; choose from {GHOST_MODES}.")
    score_range = _validate_score_range(score_range)
    metrics = resolve_metric_order(metric_order)

    rows = load_rows(input_csv)
    all_stems = [str(r.get("structure_stem", "")).strip() for r in rows]
    all_stems = [s for s in all_stems if s]
    if not all_stems:
        raise ValueError(f"No structures found in {input_csv}.")

    ordered = _order_stems(
        all_stems,
        subset=stems or [],
        order=structure_order or [],
        order_file_stems=_load_order_file(structure_order_file),
    )

    # Peek volume values for empirical mode / provenance before scoring.
    ordered_set = set(ordered)
    matthews_vals: list[float] = []
    packing_vals: list[float] = []
    solvent_vals: list[float] = []
    occlusion_values: dict[str, list[float]] = {col: [] for col in OCCLUSION_COLUMNS}
    for r in rows:
        if str(r.get("structure_stem", "")).strip() not in ordered_set:
            continue
        matthews_vals.append(_to_float(r.get("matthews_a3_per_Da")))
        packing_vals.append(_to_float(r.get("packing_density_percent")))
        solvent_vals.append(_to_float(r.get("estimated_solvent_percent")))
        for col in OCCLUSION_COLUMNS:
            val = _to_float(r.get(col))
            if math.isfinite(val):
                occlusion_values[col].append(val)
    volume_limits, scale_meta = resolve_volume_limits(
        scale_mode=scale_mode,
        matthews_values=matthews_vals,
        packing_values=packing_vals,
        solvent_values=solvent_vals,
        matthews_min=matthews_min,
        matthews_max=matthews_max,
        solvent_min=solvent_min,
        solvent_max=solvent_max,
        packing_min=packing_min,
        packing_max=packing_max,
    )
    occlusion_limits, occlusion_meta = resolve_occlusion_limits(
        occlusion_scale_mode=occlusion_scale_mode,
        occlusion_values=occlusion_values,
        loi_min=loi_min,
        loi_max=loi_max,
        contact_min=contact_min,
        contact_max=contact_max,
        bsa_kda_min=bsa_kda_min,
        bsa_kda_max=bsa_kda_max,
    )

    scores_by_stem, raw_by_stem, kept, used_limits = compute_scores(
        rows,
        ordered,
        volume_limits=volume_limits,
        occlusion_limits=occlusion_limits,
        score_range=score_range,
        metrics=metrics,
    )
    if not kept:
        raise ValueError("No structures with usable metrics to plot.")

    display_meta: dict[str, object] = {"display_mode": display_mode}
    if display_mode == "deviation":
        values_by_stem, display_meta = compute_deviations(
            kept,
            raw_by_stem,
            deviation_scale=deviation_scale,
            metrics=metrics,
        )
        limit = deviation_limit(
            kept,
            values_by_stem,
            requested=deviation_limit_value,
        )
        plot_range = (-limit, limit)
        value_suffix = "deviation"
        display_meta["deviation_limit"] = limit
        display_meta["requested_deviation_limit"] = deviation_limit_value
    else:
        values_by_stem = scores_by_stem
        plot_range = score_range
        value_suffix = "score"

    if sort_by == "mean":
        kept = sorted(kept, key=lambda s: (-_safe_sort(_mean_score(values_by_stem[s])), s))

    ghost_by_stem, ghost_label = _ghost_values(
        values_by_stem,
        kept,
        reference_stem,
        ghost_mode=ghost_mode,
        metrics=metrics,
    )

    mode = str(scale_meta.get("scale_mode", scale_mode))
    vm_lim = used_limits.get("matthews_a3_per_Da", (math.nan, math.nan))
    sol_lim = used_limits.get("estimated_solvent_percent", (math.nan, math.nan))
    pack_lim = used_limits.get("packing_density_percent", (math.nan, math.nan))
    bands = (
        f"Matthews {vm_lim[0]:g}–{vm_lim[1]:g}; "
        f"solvent {sol_lim[0]:g}–{sol_lim[1]:g}%; "
        f"packing {pack_lim[0]:g}–{pack_lim[1]:g}%"
    )
    if mode == "cryst1":
        volume_note = f"Volume axes: cryst1/unit-cell literature ({bands})."
    elif mode == "empirical":
        volume_note = (
            f"Volume axes: within-run empirical "
            f"(close ≤{scale_meta.get('empirical_close_span', EMPIRICAL_CLOSE_SPAN):g} "
            f"→ floor(min)−1…ceil(max)+1; wider → reporting-grid outward; {bands})."
        )
    elif mode == "user":
        volume_note = f"Volume axes: user-defined ({bands})."
    elif mode in ("bbox", "slayer-compact"):
        profile = scale_meta.get("volume_profile") or {}
        assert isinstance(profile, dict)
        version = profile.get("version", "")
        volume_note = (
            f"Volume axes: {mode} v{version} expanded-lattice profile ({bands})."
        )
    else:
        volume_note = f"Volume axes: {mode} ({bands})."

    clipped = clipped_values(kept, raw_by_stem, used_limits, metrics=metrics)
    # Clipping is a score-mode concern (fixed bands); in deviation mode the ring
    # is fitted to the cohort itself so band clipping does not apply to the plot.
    clip_count = len(clipped) if display_mode == "score" else 0
    occlusion_note = _occlusion_scale_note(
        occlusion_meta,
        used_limits,
        n_clipped=clip_count,
        score_range=score_range,
    )
    order_note = ""
    order_tag = metric_order_tag(metrics)
    if order_tag:
        order_note = (
            " Spoke order: "
            + ", ".join(m.short for m in metrics)
            + "."
        )
    if display_mode == "deviation":
        scale_note = (
            f"Display: cohort deviation ({deviation_scale}; "
            f"ring ±{_format_score(plot_range[1])}). "
            + volume_note
            + " "
            + occlusion_note
            + order_note
            + " (raw bands recorded for provenance; not used for the radial mapping)."
        )
    else:
        scale_note = (
            volume_note
            + " "
            + occlusion_note
            + order_note
            + f" Score range: {_format_score(score_range[0])}–{_format_score(score_range[1])}."
        )

    # Fold the resolved scale modes into the output directory and filenames so
    # runs with different scaling never overwrite one another. A single tag is
    # used when both families share a mode; a mixed run records both. Append the
    # tag to the requested directory rather than creating another nested level.
    occl_mode = str(occlusion_meta.get("occlusion_scale_mode", occlusion_scale_mode))
    tag = (
        output_run_tag(
            mode,
            occl_mode,
            display_mode=display_mode,
            deviation_scale=deviation_scale,
            ghost_mode=ghost_mode,
            reference_stem=reference_stem,
            metrics=metrics,
            edge_delta=edge_delta,
        )
        if tag_outputs
        else ""
    )
    if tag:
        output_dir = f"{output_dir}_{tag}"

    os.makedirs(output_dir, exist_ok=True)
    written: list[str] = []

    scores_csv = os.path.join(output_dir, _tagged_name("radar_scores", "csv", tag))
    _write_scores_csv(
        scores_csv,
        kept,
        values_by_stem,
        raw_by_stem,
        value_suffix=value_suffix,
        metrics=metrics,
    )
    written.append(scores_csv)

    scale_json = os.path.join(output_dir, _tagged_name("radar_scale", "json", tag))
    with open(scale_json, "w", encoding="utf-8") as f_out:
        import json

        json.dump(
            {
                **scale_meta,
                **occlusion_meta,
                **display_meta,
                "scale_tag": tag,
                "metric_order": [m.column for m in metrics],
                "metric_order_tag": order_tag or "default",
                "edge_delta": edge_delta,
                "edge_delta_scale": (
                    {
                        "mode": "plot_span",
                        "lw_min": EDGE_DELTA_LW_MIN,
                        "lw_max": EDGE_DELTA_LW_MAX,
                        "alpha_min": EDGE_DELTA_ALPHA_MIN,
                        "alpha_max": EDGE_DELTA_ALPHA_MAX,
                        "power": EDGE_DELTA_POWER,
                    }
                    if edge_delta
                    else None
                ),
                "ghost_mode": ghost_mode if not reference_stem else "reference",
                "ghost_label": ghost_label,
                "reference_stem": reference_stem,
                "score_range": list(score_range),
                "plot_range": list(plot_range),
                "used_limits": {k: list(v) for k, v in used_limits.items()},
                "clipped": clipped if display_mode == "score" else [],
                "n_clipped": clip_count,
                "structures": kept,
            },
            f_out,
            indent=2,
        )
    written.append(scale_json)

    written.extend(
        _render_grid_pages(
            kept,
            values_by_stem,
            raw_by_stem,
            ghost_by_stem,
            ghost_label,
            output_dir=output_dir,
            out_format=out_format,
            dpi=dpi,
            title=title,
            annotate_values=annotate_values,
            max_per_page=max_per_page,
            score_range=plot_range,
            scale_note=scale_note,
            tag=tag,
            display_mode=display_mode,
            metrics=metrics,
            edge_delta=edge_delta,
        )
    )
    if write_heatmap:
        written.append(
            _render_score_heatmap(
                kept,
                values_by_stem,
                output_dir=output_dir,
                out_format=out_format,
                dpi=dpi,
                title=title,
                score_range=plot_range,
                tag=tag,
                display_mode=display_mode,
                deviation_scale=deviation_scale,
                metrics=metrics,
            )
        )
    print(scale_note, file=sys.stderr)
    return written


def _safe_sort(value: float) -> float:
    return -math.inf if math.isnan(value) else value


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(
        description=(
            "Radar charts of six lattice packing/occlusion metrics. Accepts either "
            "combined_lattice_vs_ec.csv (structures as rows) or a transposed summary table "
            "(metrics as rows). Volume and occlusion axes share the same scale vocabulary "
            "(cryst1, bbox, empirical, slayer-compact, user). Matthews and solvent "
            "are inverted so a larger polygon means tighter / more occluded."
        )
    )
    ap.add_argument(
        "--input",
        required=True,
        metavar="PATH",
        help="Metrics table: combined_lattice_vs_ec.csv, or a transposed summary table "
        "with metric labels in the first column.",
    )
    ap.add_argument("-d", "--output-dir", required=True, metavar="DIR", help="Output directory.")
    ap.add_argument(
        "--stems",
        metavar="STEM,STEM,...",
        default=None,
        help="Comma-separated subset of structure stems to include (default: all).",
    )
    ap.add_argument(
        "--structure-order",
        metavar="STEM,STEM,...",
        default=None,
        help="Panel / row order (listed stems first, remainder follows input order).",
    )
    ap.add_argument(
        "--structure-order-file",
        metavar="PATH",
        default=None,
        help="One structure stem per line (applied after --structure-order).",
    )
    ap.add_argument(
        "--reference-stem",
        metavar="STEM",
        default=None,
        help="Use this structure's silhouette as the ghost outline "
        "(overrides --ghost; default ghost is leave-one-out median).",
    )
    ap.add_argument(
        "--ghost",
        choices=GHOST_MODES,
        default="loo-median",
        dest="ghost_mode",
        help=(
            "Ghost silhouette when --reference-stem is not set. "
            "'loo-median' (default) / 'loo-mean' = leave-one-out of the other n-1 "
            "structures (never coincides with the panel); "
            "'cohort-median' / 'cohort-mean' = shared silhouette including the panel."
        ),
    )
    ap.add_argument(
        "--display",
        choices=DISPLAY_MODES,
        default="score",
        dest="display_mode",
        help=(
            "Radial encoding. 'score' (default) maps each raw value onto the selected "
            "band and --score-range. 'deviation' standardises each axis about the "
            "cohort centre so small within-cohort differences use the full ring "
            "(ignores --score-range for the plot; raw bands are still recorded)."
        ),
    )
    ap.add_argument(
        "--deviation-scale",
        choices=DEVIATION_SCALES,
        default="sd",
        help=(
            "Spread estimator for --display deviation: "
            "'sd' = (x-mean)/sample SD (classic z-score); "
            "'mad' = (x-median)/(1.4826*MAD), robust to a single outlier; "
            "'range' = (x-midrange)/half-range, mapping the cohort onto [-1, +1]."
        ),
    )
    ap.add_argument(
        "--deviation-limit",
        type=float,
        default=DEFAULT_DEVIATION_LIMIT,
        metavar="Z",
        help=(
            "Symmetric ring radius for --display deviation (default: 0 = auto-fit "
            "to the observed spread, rounded up to a half-unit)."
        ),
    )
    ap.add_argument(
        "--sort-by",
        choices=("none", "mean"),
        default="none",
        help="Panel / heatmap ordering: 'mean' sorts by mean score (tighter first).",
    )
    ap.add_argument(
        "--metric-order",
        default="default",
        metavar="PRESET|TOKENS",
        help=(
            "Spoke order. Presets: 'default' (LOI, Matthews, contact, packing, "
            "BSA/kDa, solvent; LOI at 12 o'clock); "
            "'interleaved' (packing, contact, Matthews, LOI, solvent, BSA/kDa — "
            "matched packing/occlusion neighbours). Or a comma-separated list of the six "
            "tokens packing,matthews,solvent,bsa,loi,contact (each exactly once)."
        ),
    )
    ap.add_argument(
        "--scale-mode",
        choices=[*SCALE_MODES, "crystal"],
        default="bbox",
        help=(
            "Volume-axis scoring band (shared vocabulary with --occlusion-scale-mode). "
            "'bbox' = fixed expanded-lattice profile (default); "
            "'slayer-compact' = tighter SlpA-core profile; "
            "'cryst1' = unit-cell literature (alias: crystal); "
            "'empirical' = within-run ±1 pad when close, else reporting-grid outward; "
            "'user' = require --matthews-min/--matthews-max. "
            "Explicit --*-min/--*-max overrides stay exact."
        ),
    )
    ap.add_argument(
        "--matthews-min",
        type=float,
        default=None,
        metavar="X",
        help="Matthews lower limit (Å³/Da). Required for --scale-mode user; optional override otherwise.",
    )
    ap.add_argument(
        "--matthews-max",
        type=float,
        default=None,
        metavar="X",
        help="Matthews upper limit (Å³/Da). Required for --scale-mode user; optional override otherwise.",
    )
    ap.add_argument(
        "--solvent-min",
        type=float,
        default=None,
        metavar="PCT",
        help="Estimated solvent lower limit (%%). Optional; derived in user mode if omitted.",
    )
    ap.add_argument(
        "--solvent-max",
        type=float,
        default=None,
        metavar="PCT",
        help="Estimated solvent upper limit (%%). Optional; derived in user mode if omitted.",
    )
    ap.add_argument(
        "--packing-min",
        type=float,
        default=None,
        metavar="PCT",
        help="Packing density lower limit (%%). Optional; derived in user mode if omitted.",
    )
    ap.add_argument(
        "--packing-max",
        type=float,
        default=None,
        metavar="PCT",
        help="Packing density upper limit (%%). Optional; derived in user mode if omitted.",
    )
    ap.add_argument(
        "--occlusion-scale-mode",
        choices=OCCLUSION_SCALE_MODES,
        default="bbox",
        help=(
            "Occlusion-axis scoring band (shared vocabulary with --scale-mode). "
            "'bbox' = fixed expanded-lattice profile (default); "
            "'slayer-compact' = tighter SlpA-core profile; "
            "'cryst1' = falls back to within-run empirical (no crystal LOI band); "
            "'empirical' = within-run ±1 pad when close, else 0/5 grid; "
            "'cohort' = raw within-run min–max (legacy); "
            "'user' = require all three --loi/--contact/--bsa-kda pairs."
        ),
    )
    ap.add_argument(
        "--loi-min",
        type=float,
        default=None,
        metavar="PCT",
        help="LOImolA lower limit (%%). Overrides the profile band; exact (not rounded).",
    )
    ap.add_argument(
        "--loi-max",
        type=float,
        default=None,
        metavar="PCT",
        help="LOImolA upper limit (%%). Overrides the profile band; exact (not rounded).",
    )
    ap.add_argument(
        "--contact-min",
        type=float,
        default=None,
        metavar="PCT",
        help="Lattice contact-residue lower limit (%%). Overrides the profile band; exact.",
    )
    ap.add_argument(
        "--contact-max",
        type=float,
        default=None,
        metavar="PCT",
        help="Lattice contact-residue upper limit (%%). Overrides the profile band; exact.",
    )
    ap.add_argument(
        "--bsa-kda-min",
        type=float,
        default=None,
        metavar="X",
        help="BSAmolA/kDa lower limit (Å²). Overrides the profile band; exact.",
    )
    ap.add_argument(
        "--bsa-kda-max",
        type=float,
        default=None,
        metavar="X",
        help="BSAmolA/kDa upper limit (Å²). Overrides the profile band; exact.",
    )
    ap.add_argument(
        "--no-annotate-values",
        action="store_true",
        help="Do not print raw metric values at each spoke (annotations on by default).",
    )
    ap.add_argument(
        "--score-range",
        nargs=2,
        type=float,
        default=DEFAULT_SCORE_RANGE,
        metavar=("MIN", "MAX"),
        help="Score/radar range as MIN MAX (default: 0 10). Raw metric bands are "
        "unchanged; this only remaps their scores, rings, heatmap, and CSV values.",
    )
    ap.add_argument(
        "--format",
        choices=_VALID_FORMATS,
        default="png",
        help="Output image format (default: png).",
    )
    ap.add_argument("--dpi", type=int, default=200, help="Raster output DPI (default: 200).")
    ap.add_argument("--title", default=None, help="Optional figure title.")
    ap.add_argument(
        "--max-per-page",
        type=int,
        default=0,
        help="Max radar panels per page (0 = single page). Splits into "
        "radar_grid_<tag>_01.<fmt>, ...",
    )
    ap.add_argument(
        "--edge-delta",
        action="store_true",
        help=(
            "Encode adjacent spoke disagreement on polygon edges: linewidth and opacity "
            "scale with |Δplot| / (smax−smin) using a shared mapping across all panels "
            "(highlights zig-zags; ghost stays uniform)."
        ),
    )
    ap.add_argument(
        "--no-heatmap",
        action="store_true",
        help="Skip the companion score/deviation heatmap (radar grid, scores CSV, "
        "and scale JSON are still written).",
    )
    ap.add_argument(
        "--no-scale-tag",
        action="store_true",
        help="Write outputs directly in --output-dir without the scale-mode "
        "directory suffix/filename tag (default: tag by scale mode, for example "
        "<output-dir>_bbox/radar_grid_bbox.png; mixed runs use "
        "vol-<mode>_occ-<mode>).",
    )
    args = ap.parse_args(argv)

    try:
        written = generate_radar(
            input_csv=os.path.abspath(args.input),
            output_dir=os.path.abspath(args.output_dir),
            stems=_parse_csv_list(args.stems),
            structure_order=_parse_csv_list(args.structure_order),
            structure_order_file=args.structure_order_file,
            reference_stem=(args.reference_stem or None),
            ghost_mode=args.ghost_mode,
            display_mode=args.display_mode,
            deviation_scale=args.deviation_scale,
            deviation_limit_value=float(args.deviation_limit),
            sort_by=args.sort_by,
            annotate_values=not args.no_annotate_values,
            out_format=args.format,
            dpi=int(args.dpi),
            title=args.title,
            max_per_page=int(args.max_per_page),
            score_range=(float(args.score_range[0]), float(args.score_range[1])),
            tag_outputs=not args.no_scale_tag,
            scale_mode=args.scale_mode,
            matthews_min=args.matthews_min,
            matthews_max=args.matthews_max,
            solvent_min=args.solvent_min,
            solvent_max=args.solvent_max,
            packing_min=args.packing_min,
            packing_max=args.packing_max,
            occlusion_scale_mode=args.occlusion_scale_mode,
            loi_min=args.loi_min,
            loi_max=args.loi_max,
            contact_min=args.contact_min,
            contact_max=args.contact_max,
            bsa_kda_min=args.bsa_kda_min,
            bsa_kda_max=args.bsa_kda_max,
            metric_order=args.metric_order,
            write_heatmap=not args.no_heatmap,
            edge_delta=bool(args.edge_delta),
        )
    except ValueError as exc:
        ap.error(str(exc))
        return 2  # unreachable; ap.error raises SystemExit
    except ImportError as exc:
        print(f"Error: matplotlib/numpy required for radar output: {exc}", file=sys.stderr)
        return 3

    print("Wrote:\n  " + "\n  ".join(written))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
