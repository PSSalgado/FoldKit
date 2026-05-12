#!/usr/bin/env python3
"""
Shared argparse helpers for rectangular interface-matrix heatmaps
(interface_analysis_matrix.py, lattice_compare_batch.py).

Per-metric colour/scale overrides and generic heatmap flags are defined in
``utils/foldkit_heatmap.py`` (:func:`~utils.foldkit_heatmap.add_generic_heatmap_args`,
:func:`~utils.foldkit_heatmap.add_per_metric_heatmap_override_args`). This module adds
interface-only options (axis transpose, BSA robust split, cell annotations, …).
"""

from __future__ import annotations

import argparse
from collections.abc import Sequence, Set
METRIC_CHOICES: tuple[str, ...] = (
    "bsa",
    "ec",
    "ec_density",
    "contacts",
    "charged_opposite",
    "charged_same",
)

DEFAULT_METRICS: tuple[str, ...] = ("bsa", "ec", "ec_density", "contacts")


def heatmap_cell_text_map_from_cli(
    metric_choices: Sequence[str],
    items: list[str],
    *,
    cell_text_allowed: Set[str],
) -> dict[str, str]:
    """Parse METRIC=FIELD tokens. FIELD must be 'same' or allowed per-interface column name."""
    mc = tuple(metric_choices)
    mc_set = frozenset(mc)
    out: dict[str, str] = {}
    for raw in items:
        item = raw.strip()
        if "=" not in item:
            raise ValueError(f"Invalid --heatmap-cell-text {raw!r}; expected METRIC=FIELD.")
        mk, fk = item.split("=", 1)
        mk, fk = mk.strip(), fk.strip()
        if mk not in mc_set:
            raise ValueError(
                f"Unknown heatmap metric {mk!r} in --heatmap-cell-text (choices: {', '.join(mc)})."
            )
        if fk not in cell_text_allowed:
            raise ValueError(
                f"Unknown annotation field {fk!r} in --heatmap-cell-text "
                f"(use 'same' or a per-interface CSV column name)."
            )
        out[mk] = fk
    return out


def add_interface_matrix_heatmap_argument_group(
    ap: argparse.ArgumentParser,
    metric_choices: Sequence[str],
    *,
    default_structures_x_axis: bool = False,
    default_colorbar_orientation: str = "vertical",
    default_short_labels: bool = False,
) -> argparse._ArgumentGroup:
    hg = ap.add_argument_group(
        "heatmap",
        "Figure options (same spirit as utils/foldkit_heatmap.py and ranking/rmsd_to_csv.py).",
    )
    try:
        from utils.foldkit_heatmap import (  # noqa: PLC0415
            add_generic_heatmap_args,
            add_per_metric_heatmap_override_args,
        )
    except ImportError:
        add_generic_heatmap_args = None
        add_per_metric_heatmap_override_args = None
    hg.add_argument(
        "--heatmap-structures-x-axis",
        action=argparse.BooleanOptionalAction,
        default=bool(default_structures_x_axis),
        help=(
            "When true: structures on the x-axis and interfaces on the y-axis (transpose matrices and heatmaps). "
            "Use --no-heatmap-structures-x-axis for the legacy layout (interfaces on x). "
            "lattice_compare_batch defaults to true."
        ),
    )
    hg.add_argument(
        "--heatmap-x-order",
        default=None,
        metavar="A,B,C,…",
        help=(
            "Comma-separated order for the x-axis labels (partial lists allowed). "
            "Applies to structures when using --heatmap-structures-x-axis; otherwise applies to "
            "interface canonical pairs."
        ),
    )
    if add_per_metric_heatmap_override_args is not None:
        add_per_metric_heatmap_override_args(hg, metric_choices, prefix="heatmap-")
    if add_generic_heatmap_args is not None:
        add_generic_heatmap_args(
            hg,
            prefix="heatmap-",
            default_colorbar_orientation=default_colorbar_orientation,
            default_short_labels=default_short_labels,
        )
    hg.add_argument(
        "--heatmap-annotate-metrics",
        nargs="*",
        choices=list(metric_choices),
        default=[],
        metavar="NAME",
        help=(
            "Draw one numeric label per cell on these heatmaps (repeatable metric names). "
            "Defaults to the colour value; use --heatmap-cell-text to show a different CSV column."
        ),
    )
    hg.add_argument(
        "--heatmap-cell-text",
        action="append",
        default=None,
        metavar="METRIC=FIELD",
        help=(
            "For heatmap METRIC, print FIELD in each cell (colour still follows that heatmap's metric). "
            "FIELD is a per-interface CSV column name, or 'same' for the colour values (default when omitted). "
            "Example: --heatmap-cell-text ec=ec_density_per_1000_A2 colours by EC (r) with EC per 1000 Å² labels. "
            "Repeat the flag for multiple metrics."
        ),
    )
    hg.add_argument(
        "--heatmap-annotate-fmt",
        default="{:.2f}",
        metavar="FMT",
        help="Format string for single-metric cell annotations (default: '{:.2f}').",
    )
    hg.add_argument(
        "--heatmap-annotate-fontsize",
        type=float,
        default=7.0,
        metavar="PT",
        help="Font size for annotated heatmap cells (default: 7).",
    )
    hg.add_argument(
        "--heatmap-bsa-positive-range",
        action="store_true",
        help=(
            "For BSA heatmaps only: set default vmin/vmax from positive finite cells "
            "(RMSD-style autoscale)."
        ),
    )
    hg.add_argument(
        "--heatmap-bsa-annotate-contacts",
        action="store_true",
        help=(
            "Convenience: for the BSA heatmap, annotate cells with contact_count values "
            "(colour still uses BSA). Equivalent to: --heatmap-annotate-metrics bsa "
            "--heatmap-cell-text bsa=contact_count."
        ),
    )
    hg.add_argument(
        "--heatmap-bsa-robust",
        action="store_true",
        help=(
            "BSA only: split the colour bar between the bulk and any dominant interface(s). A "
            "canonical chain pair is auto-flagged as an outlier when its minimum BSA across "
            "structures exceeds --heatmap-bsa-outlier-factor (default 3) times the median of "
            "the remaining cells; the cap is the rounded max of non-outlier cells. Cells below "
            "the cap use the BSA colormap (--heatmap-cmap-by-metric bsa=… if set, else "
            "--heatmap-cmap); above use --heatmap-bsa-above-cmap (default Reds). "
            "Above-cap pairs are listed in a figure note."
        ),
    )
    hg.add_argument(
        "--heatmap-bsa-outlier-factor",
        type=float,
        default=3.0,
        metavar="K",
        help=(
            "Multiplier (> 1) for the auto outlier rule: a canonical chain pair is an outlier "
            "when its min positive BSA across structures > K x median(remaining cells). "
            "Default: 3.0."
        ),
    )
    hg.add_argument(
        "--heatmap-bsa-split-at",
        type=float,
        default=None,
        metavar="VALUE",
        help=(
            "BSA only: explicit colour-bar split value (Å²). Wins over auto detection. Use "
            "when the auto rule does not isolate the desired interface(s)."
        ),
    )
    hg.add_argument(
        "--heatmap-bsa-above-cmap",
        default="Reds",
        metavar="NAME",
        help="Matplotlib colour map for cells above the BSA cap (default: Reds).",
    )
    hg.add_argument(
        "--heatmap-bsa-no-outlier-note",
        action="store_true",
        help="Disable the figure-level note listing above-cap BSA pairs.",
    )
    return hg


def finalize_interface_matrix_heatmap_args(
    args: argparse.Namespace,
    ap: argparse.ArgumentParser,
    *,
    metric_choices: Sequence[str],
    cell_text_allowed: Set[str],
) -> None:
    """Populate derived Namespace attributes; call ``ap.error`` on invalid CLI."""
    try:
        args.heatmap_cell_text_map = heatmap_cell_text_map_from_cli(
            metric_choices,
            list(args.heatmap_cell_text or []),
            cell_text_allowed=cell_text_allowed,
        )
    except ValueError as exc:
        ap.error(str(exc))
    try:
        from utils.foldkit_heatmap import (  # noqa: PLC0415
            metric_edges_map_from_cli,
            metric_float_map_from_cli,
            metric_str_map_from_cli,
        )

        args.heatmap_vmin_by_metric = metric_float_map_from_cli(
            metric_choices,
            list(args.heatmap_vmin_by_metric or []),
        )
        args.heatmap_vmax_by_metric = metric_float_map_from_cli(
            metric_choices,
            list(args.heatmap_vmax_by_metric or []),
        )
        args.heatmap_scale_by_metric = metric_str_map_from_cli(
            metric_choices,
            list(args.heatmap_scale_by_metric or []),
            allowed=frozenset({"linear", "log10", "log1p", "clip_p95", "clip_p98"}),
        )
        args.heatmap_boundaries_by_metric = metric_edges_map_from_cli(
            metric_choices,
            list(args.heatmap_boundaries_by_metric or []),
        )
        args.heatmap_cmap_by_metric = metric_str_map_from_cli(
            metric_choices,
            list(args.heatmap_cmap_by_metric or []),
            allowed=None,
        )
        args.heatmap_diverging_center_by_metric = metric_str_map_from_cli(
            metric_choices,
            list(getattr(args, "heatmap_diverging_center_by_metric", None) or []),
            allowed=frozenset({"none", "median"}),
            normalize_value=str.lower,
        )
    except ValueError as exc:
        ap.error(str(exc))

    from utils.foldkit_heatmap import parse_boundaries_csv as parse_boundaries_csv_main  # noqa: PLC0415

    gb_edges = parse_boundaries_csv_main(getattr(args, "heatmap_boundaries", None))
    if len(gb_edges) == 1:
        ap.error("--heatmap-boundaries needs at least two distinct comma-separated edge values.")

    if bool(getattr(args, "heatmap_bsa_annotate_contacts", False)):
        ann = list(getattr(args, "heatmap_annotate_metrics", None) or [])
        if "bsa" not in ann:
            ann.append("bsa")
        args.heatmap_annotate_metrics = ann
        m = dict(getattr(args, "heatmap_cell_text_map", None) or {})
        m.setdefault("bsa", "contact_count")
        args.heatmap_cell_text_map = m


def apply_default_lattice_compare_heatmap_annotations(args: argparse.Namespace) -> None:
    """Default BSA + EC annotations for lattice_compare_batch (contact count / EC per 1000 Å²)."""
    ann = list(getattr(args, "heatmap_annotate_metrics", None) or [])
    if not ann:
        ann = ["bsa", "ec"]
        args.heatmap_annotate_metrics = ann
    m = dict(getattr(args, "heatmap_cell_text_map", None) or {})
    if "bsa" in ann:
        m.setdefault("bsa", "contact_count")
    if "ec" in ann:
        m.setdefault("ec", "ec_density_per_1000_A2")
    args.heatmap_cell_text_map = m
