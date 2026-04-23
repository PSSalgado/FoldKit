#!/usr/bin/env python3
"""
Matplotlib square-matrix heatmaps for FoldKit (RMSD, Dali Z, …).

Used by ``rmsd_to_csv`` and ``dalilite_pairs`` without those entry points
importing each other. Documented in the main ``README.md`` (``rmsd_to_csv.py`` and
``dalilite_pairs.py`` sections) and in ``metrics/metrics_details.md`` §6.3.1.
CLI details: ``python ranking/rmsd_to_csv.py --help`` and
``python ranking/dalilite_pairs.py --help``.
"""
from __future__ import annotations

import argparse
import os
import re
import sys
import tempfile

from structure_phylogeny import _natural_sort_key


def short_heatmap_label(name: str) -> str:
    """
    Shorten model names for heatmap axis labels only (CSV data unchanged).

    Example: S1refAF3_CWB2s.pdb -> S1_CWB2s — leading letter+digits from the first underscore
    segment, plus the segment after the last underscore. If the pattern does not apply, returns
    the basename without extension.
    """
    s = (name or "").strip()
    if not s:
        return s
    base = os.path.splitext(s)[0] if "." in s else s
    parts = base.split("_")
    if len(parts) < 2:
        return base
    m = re.match(r"^([A-Za-z]\d+)", parts[0])
    if not m:
        return base
    return f"{m.group(1)}_{parts[-1]}"


def median_rmsd_from_square_array(arr) -> float | None:
    """Median of off-diagonal finite RMSDs (prefer positive), before diagonal masking."""
    try:
        import numpy as np
    except ImportError:
        return None
    n = int(arr.shape[0])
    if n < 2:
        return None
    mask = ~np.eye(n, dtype=bool)
    vals = arr[mask]
    vals = vals[np.isfinite(vals)]
    pos = vals[vals > 0]
    use = pos if len(pos) > 0 else vals
    if len(use) == 0:
        return None
    return float(np.median(use))


def median_finite_offdiag_square_array(arr) -> float | None:
    """Median of all finite off-diagonal entries (for Z-score or signed matrices)."""
    try:
        import numpy as np
    except ImportError:
        return None
    n = int(arr.shape[0])
    if n < 2:
        return None
    mask = ~np.eye(n, dtype=bool)
    vals = arr[mask]
    vals = vals[np.isfinite(vals)]
    if len(vals) == 0:
        return None
    return float(np.median(vals))


def make_heatmap_norm(
    disp_vmin: float,
    disp_vmax: float,
    diverging_center: str | None,
    vcenter: float | None,
):
    """
    Linear ``Normalize``, or ``TwoSlopeNorm`` with ``vcenter`` at the median for diverging colour maps.
    If vcenter is outside (vmin, vmax) or degenerate, falls back to linear scale.
    """
    import matplotlib.colors as mcolors

    if diverging_center is None or diverging_center == "none" or vcenter is None:
        return mcolors.Normalize(vmin=disp_vmin, vmax=disp_vmax)
    vc = float(vcenter)
    span = disp_vmax - disp_vmin
    eps = max(1e-12 * max(abs(disp_vmax), abs(disp_vmin), 1.0), 1e-15 * max(span, 1.0))
    lo, hi = disp_vmin + eps, disp_vmax - eps
    if lo >= hi:
        return mcolors.Normalize(vmin=disp_vmin, vmax=disp_vmax)
    vc = min(max(vc, lo), hi)
    return mcolors.TwoSlopeNorm(vmin=disp_vmin, vcenter=vc, vmax=disp_vmax)


def _nice_tick_values_linear_segment(vmin: float, vmax: float, nbins: int = 6):
    """Return sorted unique tick values in [vmin, vmax] for linear scales (Å)."""
    import numpy as np
    import matplotlib.ticker as mticker

    if vmax < vmin:
        vmin, vmax = vmax, vmin
    if not (np.isfinite(vmin) and np.isfinite(vmax)):
        return np.array([vmin, vmax], dtype=float)
    if vmax - vmin < 1e-30:
        return np.array([vmin, vmax], dtype=float)
    loc = mticker.MaxNLocator(nbins=nbins + 1, min_n_ticks=2, symmetric=False)
    ticks = loc.tick_values(vmin, vmax)
    ticks = np.asarray(ticks, dtype=float)
    ticks = ticks[np.isfinite(ticks)]
    ticks = ticks[(ticks >= vmin - 1e-12) & (ticks <= vmax + 1e-12)]
    if ticks.size == 0:
        return np.array([vmin, vmax], dtype=float)
    return np.unique(ticks)


def _set_heatmap_colorbar_ticks(cbar, norm) -> None:
    """
    Set colour-bar ticks in data units. TwoSlopeNorm (median centre) is non-linear in display
    space, so we pick nice ticks separately below and above vcenter instead of a uniform locator.
    """
    import matplotlib.colors as mcolors
    import matplotlib.ticker as mticker
    import numpy as np

    vmin, vmax = float(norm.vmin), float(norm.vmax)
    if isinstance(norm, mcolors.TwoSlopeNorm):
        vc = float(norm.vcenter)
        span = max(vmax - vmin, 1e-30)
        if abs(vc - vmin) < 1e-9 * span or abs(vc - vmax) < 1e-9 * span:
            ticks = _nice_tick_values_linear_segment(vmin, vmax)
        else:
            t_lo = _nice_tick_values_linear_segment(vmin, vc, nbins=5)
            t_hi = _nice_tick_values_linear_segment(vc, vmax, nbins=5)
            ticks = np.unique(np.concatenate([t_lo, t_hi]))
    else:
        ticks = _nice_tick_values_linear_segment(vmin, vmax, nbins=8)

    cbar.set_ticks(ticks)
    fmt = mticker.ScalarFormatter()
    fmt.set_scientific(False)
    fmt.set_useOffset(False)
    orient = getattr(cbar, "orientation", "vertical")
    if orient == "vertical":
        cbar.ax.yaxis.set_major_formatter(fmt)
    else:
        cbar.ax.xaxis.set_major_formatter(fmt)


def _apply_heatmap_y_axis_right(ax) -> None:
    """Put row (y) tick labels on the right; hide y ticks on the left."""
    ax.yaxis.set_ticks_position("right")
    ax.yaxis.set_label_position("right")
    ax.tick_params(axis="y", which="major", left=False, right=True, labelleft=False, labelright=True)


def _apply_label_order(
    args: argparse.Namespace,
    ids: list[str],
    matrix: list[list[float]],
    *,
    context: str,
) -> tuple[list[str], list[list[float]]]:
    """Apply --order / --order-dir to ids and matrix; context is a label for error messages."""
    if args.order_dir:
        glob_pat = args.order_glob or "*.pdb"
        root_abs = os.path.abspath(args.order_dir)
        if not os.path.isdir(root_abs):
            print("Error: --order-dir is not a directory:", root_abs, file=sys.stderr)
            sys.exit(1)
        order_list = order_labels_from_scan(root_abs, glob_pat)
        if order_list:
            ids, matrix = reorder_matrix(ids, matrix, order_list)
            print(
                f"[{context}] Applied --order-dir {root_abs} glob {glob_pat!r}: {len(order_list)} labels",
                file=sys.stderr,
            )
        else:
            print(
                f"Warning [{context}]: no files matched under --order-dir; using default matrix order.",
                file=sys.stderr,
            )
    elif args.order:
        order_list = parse_order(args.order)
        if order_list:
            ids, matrix = reorder_matrix(ids, matrix, order_list)
            print(f"[{context}] Applied --order: {len(order_list)} labels", file=sys.stderr)
        else:
            print(
                f"Warning [{context}]: --order empty or file not found; using default order.",
                file=sys.stderr,
            )
    elif args.order_glob:
        print("Error: --order-glob requires --order-dir.", file=sys.stderr)
        sys.exit(1)
    return ids, matrix


def _image_magic_ok(path: str, fmt: str) -> bool:
    """Return True if the file looks like the claimed image format (magic bytes)."""
    if not os.path.exists(path) or os.path.getsize(path) < 8:
        return False
    with open(path, "rb") as f:
        head = f.read(12)
    if fmt == "png" and head[:8] == b"\x89PNG\r\n\x1a\n":
        return True
    if fmt == "pdf" and head[:4] == b"%PDF":
        return True
    if fmt == "svg" and (b"<?xml" in head or b"<svg" in head[:100]):
        return True
    return False


def _heatmap_vector_format(fmt: str) -> bool:
    """Use per-cell vector graphics for these backends (imshow embeds a bitmap in SVG/PDF)."""
    return fmt in ("svg", "pdf")


# Second-channel hatch: cycle none → vertical → horizontal → diagonal (thin black lines).
# Matplotlib tightens hatch spacing by repeating | - / ; too few repeats look empty, too many
# merge into solid black (especially on small cells). Tune these if your matrix size changes a lot.
_HATCH_LINEWIDTH_PT = 0.15
_HATCH_LINE_COLOR = "grey"
_HATCH_REPEAT_CELLS = 4
_HATCH_REPEAT_LEGEND = 5

# Cell border: draw a thin black grid around each cell.
_CELL_BORDER_LINEWIDTH_PT = 0.3
_CELL_BORDER_COLOR = "black"


def _hatch_pattern_for_bin_index(
    bin_index: int,
    *,
    legend: bool = False,
    repeat_cells: int = _HATCH_REPEAT_CELLS,
    repeat_legend: int = _HATCH_REPEAT_LEGEND,
) -> str | None:
    """Bin 0: no hatch; then dense |, -, /; repeat for higher bin indices."""
    k = int(bin_index) % 4
    if k == 0:
        return None
    n = int(repeat_legend) if legend else int(repeat_cells)
    n = max(1, n)
    if k == 1:
        return "|" * n
    if k == 2:
        return "-" * n
    return "/" * n


def _hatch_bin_edges_from_values(hatch_val_arr, n_equal_bins: int) -> list[float] | None:
    """Build bin boundaries [lo, …, hi] for n_equal_bins intervals from off-diagonal hatch values."""
    import numpy as np

    m = ~np.eye(hatch_val_arr.shape[0], dtype=bool)
    vals = hatch_val_arr[m]
    vals = vals[np.isfinite(vals)]
    if len(vals) == 0:
        return None
    lo, hi = float(np.min(vals)), float(np.max(vals))
    if hi <= lo + 1e-9:
        # One effective bin; use a small synthetic span for searchsorted
        return [lo, lo + 1.0]
    n_edges = n_equal_bins + 1
    return [float(x) for x in np.linspace(lo, hi, n_edges)]


def _hatch_value_bin_index(v: float, full_edges: list[float]) -> int:
    """Index 0..len(full_edges)-2 for v in the corresponding half-open interval; clipped."""
    import numpy as np

    if not full_edges or len(full_edges) < 2:
        return 0
    n_bins = len(full_edges) - 1
    idx = int(np.searchsorted(np.asarray(full_edges, dtype=float), v, side="right")) - 1
    return int(np.clip(idx, 0, n_bins - 1))


def _hatch_bin_legend_label(
    b: int, full_edges: list[float], value_name: str
) -> str:
    a, c = full_edges[b], full_edges[b + 1]
    return f"{value_name} {a:.0f}–{c:.0f}"


def _relative_luminance_srgb(rgba: tuple[float, float, float, float]) -> float:
    """
    Approximate perceived brightness from an sRGB RGBA tuple in [0,1].
    Used to choose a contrasting hatch colour on dark cells.
    """
    r, g, b, _a = [float(x) for x in rgba]
    # Simple Rec.709 luminance in display space (good enough for hatch contrast).
    return 0.2126 * r + 0.7152 * g + 0.0722 * b


def _format_cell_value(v: float, fmt: str) -> str:
    try:
        return fmt.format(v)
    except Exception:
        # Fallback: tolerate printf-style formats if given accidentally.
        try:
            return fmt % v
        except Exception:
            return str(v)


def _draw_heatmap_vector_cells(
    ax,
    arr,
    cmap_obj,
    norm,
    *,
    hatch_value_arr: object | None = None,
    hatch_bin_edges: list[float] | None = None,
    hatch_channel_name: str = "Hatch values",
    hatch_repeat_cells: int = _HATCH_REPEAT_CELLS,
    hatch_repeat_legend: int = _HATCH_REPEAT_LEGEND,
    hatch_linewidth_pt: float = _HATCH_LINEWIDTH_PT,
    cell_border_linewidth_pt: float = _CELL_BORDER_LINEWIDTH_PT,
    cell_border_color: str = _CELL_BORDER_COLOR,
    hatch_color_light: str = "black",
    hatch_color_dark: str = "white",
    hatch_color_luma_threshold: float = 0.45,
    annotate_value_arr: object | None = None,
    annotate_fmt: str = "{:.0f}",
    annotate_fontsize: float = 6.0,
    annotate_color_light: str = "black",
    annotate_color_dark: str = "white",
    annotate_color_luma_threshold: float = 0.45,
    annotate_hatch_value_arr: object | None = None,
    annotate_hatch_fmt: str = "{:.0f}",
    annotate_hatch_fontsize: float | None = None,
    triangle: str | None = None,
) -> None:
    """
    Draw one matplotlib Rectangle per matrix cell, matching imshow(origin='lower') layout.
    Row 0 is at the bottom, column 0 at the left, so the diagonal runs bottom-left → top-right.
    Each patch has a stable SVG id (gid) for editing. No bitmaps.
    ``norm`` is a matplotlib Normalize or TwoSlopeNorm (vmin/vmax for under-colour).
    If ``hatch_value_arr`` and ``hatch_bin_edges`` are set, hatches mark binned second-channel
    values; the main matrix still sets face colour (e.g. RMSD, Dali Z).
    """
    import matplotlib.colors as mcolors
    import numpy as np
    from matplotlib.patches import Rectangle

    nrow, ncol = int(arr.shape[0]), int(arr.shape[1])
    vmin_plot = float(norm.vmin)
    under_rgba = mcolors.to_rgba(cmap_obj.get_under())
    use_hatch = (
        hatch_value_arr is not None
        and hatch_bin_edges is not None
        and len(hatch_bin_edges) >= 2
    )
    use_annot = annotate_value_arr is not None
    dual_annot = annotate_hatch_value_arr is not None
    n_bins = (len(hatch_bin_edges) - 1) if use_hatch else 0
    afs2 = float(annotate_hatch_fontsize) if annotate_hatch_fontsize is not None else max(
        3.25, float(annotate_fontsize) * 0.88
    )
    for i in range(nrow):
        for j in range(ncol):
            # Note: we plot with imshow(origin='lower') semantics (row 0 at bottom).
            # The diagonal runs bottom-left → top-right. The region "under" that diagonal
            # (bottom-right half, including diagonal) corresponds to i <= j.
            if triangle == "lower" and i > j:
                continue
            v = arr[i, j]
            if np.isnan(v):
                fc = (1.0, 1.0, 1.0, 1.0)
            elif v < vmin_plot:
                fc = under_rgba
            else:
                fc = cmap_obj(norm(v))
            y0 = i - 0.5
            x0 = j - 0.5
            hatch: str | None = None
            edgecolor = cell_border_color
            linewidth = float(cell_border_linewidth_pt)
            if use_hatch and i != j:
                hv = hatch_value_arr[i, j]
                if np.isfinite(hv):
                    bi = _hatch_value_bin_index(float(hv), hatch_bin_edges)
                    hatch = _hatch_pattern_for_bin_index(
                        bi,
                        legend=False,
                        repeat_cells=hatch_repeat_cells,
                        repeat_legend=hatch_repeat_legend,
                    )
            rect = Rectangle(
                (x0, y0),
                1.0,
                1.0,
                linewidth=linewidth,
                edgecolor=edgecolor,
                facecolor=fc,
                hatch=None,
            )
            rect.set_gid(f"rmsd_cell_r{i}_c{j}")
            ax.add_patch(rect)
            if hatch is not None:
                # Matplotlib uses the patch edge colour for hatch stroke colour.
                # To keep a stable black grid border while adapting hatch contrast,
                # draw the hatch as a separate transparent overlay patch.
                luma = _relative_luminance_srgb(fc)
                hatch_ec = hatch_color_dark if luma <= float(hatch_color_luma_threshold) else hatch_color_light
                hrect = Rectangle(
                    (x0, y0),
                    1.0,
                    1.0,
                    linewidth=0.0,
                    edgecolor=hatch_ec,
                    facecolor="none",
                    hatch=hatch,
                )
                hrect.set_gid(f"rmsd_hatch_r{i}_c{j}")
                ax.add_patch(hrect)
            if use_annot and (triangle != "lower" or i <= j) and i != j:
                av = annotate_value_arr[i, j]
                luma = _relative_luminance_srgb(fc)
                tc = (
                    annotate_color_dark
                    if luma <= float(annotate_color_luma_threshold)
                    else annotate_color_light
                )
                if dual_annot:
                    aha = annotate_hatch_value_arr[i, j] if annotate_hatch_value_arr is not None else np.nan
                    if not (np.isfinite(av) or np.isfinite(aha)):
                        continue
                    if np.isfinite(av) and np.isfinite(aha):
                        text = (
                            f"{_format_cell_value(float(av), annotate_fmt)}\n"
                            f"{_format_cell_value(float(aha), annotate_hatch_fmt)}"
                        )
                    elif np.isfinite(av):
                        text = _format_cell_value(float(av), annotate_fmt)
                    else:
                        text = _format_cell_value(float(aha), annotate_hatch_fmt)
                    ax.text(
                        j,
                        i,
                        text,
                        ha="center",
                        va="center",
                        fontsize=afs2,
                        linespacing=0.92,
                        color=tc,
                        clip_on=True,
                    )
                elif np.isfinite(av):
                    ax.text(
                        j,
                        i,
                        _format_cell_value(float(av), annotate_fmt),
                        ha="center",
                        va="center",
                        fontsize=float(annotate_fontsize),
                        color=tc,
                        clip_on=True,
                    )
    ax.set_xlim(-0.5, ncol - 0.5)
    ax.set_ylim(-0.5, nrow - 0.5)
    # Keep cells square (same scale on x and y); avoids stretch when a horizontal colour bar reshapes the axes box.
    ax.set_aspect("equal", adjustable="box")
    if use_hatch and n_bins > 0:
        from matplotlib.patches import Patch

        name = (hatch_channel_name or "Hatch values").strip() or "Hatch values"
        handles: list[Patch] = []
        for b in range(n_bins):
            hp = _hatch_pattern_for_bin_index(
                b,
                legend=True,
                repeat_cells=hatch_repeat_cells,
                repeat_legend=hatch_repeat_legend,
            )
            handles.append(
                Patch(
                    facecolor="0.85",
                    edgecolor=cell_border_color if hp is not None else "none",
                    linewidth=float(cell_border_linewidth_pt) if hp is not None else 0.0,
                    hatch=hp,
                    label=_hatch_bin_legend_label(b, hatch_bin_edges, name),
                )
            )
        leg = ax.legend(
            handles=handles,
            loc="lower left",
            bbox_to_anchor=(0.01, 0.99),
            borderaxespad=0.0,
            frameon=True,
            fancybox=True,
            framealpha=0.92,
            title=f"Hatch: {name}",
            fontsize=6,
            title_fontsize=7,
        )
        leg.set_in_layout(False)


def _axes_lower_edge_figure(fig, ax) -> float:
    """Lowest y in figure coordinates (0–1, origin bottom) of the axes tight bbox (includes tick labels)."""
    import numpy as np

    fig.canvas.draw()
    r = fig.canvas.get_renderer()
    bb = ax.get_tightbbox(r)
    if bb is None:
        pos = ax.get_position()
        return float(pos.y0)
    corners = np.array([[bb.x0, bb.y0], [bb.x1, bb.y0]], dtype=float)
    fig_y = fig.transFigure.inverted().transform(corners)[:, 1]
    return float(np.min(fig_y))


def _heatmap_data_bbox_in_figure(ax, nrows: int, ncols: int) -> tuple[float, float, float, float]:
    """
    Bounding box in figure coordinates (0–1) for the heatmap cell region
    x in [-0.5, ncols-0.5], y in [-0.5, nrows-0.5] — matches the square grid, not the full axes box.
    """
    import numpy as np

    fig = ax.figure
    fig.canvas.draw()
    pts = np.array(
        [
            [-0.5, -0.5],
            [ncols - 0.5, -0.5],
            [ncols - 0.5, nrows - 0.5],
            [-0.5, nrows - 0.5],
        ],
        dtype=float,
    )
    disp = ax.transData.transform(pts)
    f = fig.transFigure.inverted().transform(disp)
    xmin, xmax = float(np.min(f[:, 0])), float(np.max(f[:, 0]))
    ymin, ymax = float(np.min(f[:, 1])), float(np.max(f[:, 1]))
    return xmin, ymin, xmax - xmin, ymax - ymin


def _add_raster_colorbar(
    ax,
    cmap_obj,
    norm,
    *,
    nrows: int,
    ncols: int,
    label: str = "RMSD (Å)",
    orientation: str = "vertical",
    mappable=None,
    vertical_colorbar_on_left: bool = False,
    ticks: "list[float] | None" = None,
):
    """
    Colour bar as a raster, sized to match the heatmap grid extent (vertical: same height as the
    square; horizontal: same width), not the full axes.
    Call after tight_layout() on the heatmap axes.
    If vertical_colorbar_on_left is True, place the vertical bar to the left of the heatmap grid
    (e.g. when the y-axis labels are on the right).
    """
    import numpy as np
    from matplotlib.cm import ScalarMappable

    fig = ax.figure
    orient = orientation.lower() if orientation else "vertical"
    if orient not in ("vertical", "horizontal"):
        orient = "vertical"

    if mappable is None:
        sm = ScalarMappable(cmap=cmap_obj, norm=norm)
        sm.set_array(np.linspace(float(norm.vmin), float(norm.vmax), 256))
        mappable = sm

    fig.canvas.draw()
    xmin, ymin, bw, bh = _heatmap_data_bbox_in_figure(ax, nrows, ncols)
    if bw < 1e-9 or bh < 1e-9:
        pos = ax.get_position()
        xmin, ymin, bw, bh = pos.x0, pos.y0, pos.width, min(pos.width, pos.height)

    pad = 0.012
    if orient == "vertical":
        cw = max(0.007, min(0.024, 0.016 * ncols / max(nrows, ncols, 8)))
        if vertical_colorbar_on_left:
            right = xmin - pad
            left = right - cw
            if left < 0.03:
                cur = float(fig.subplotpars.left)
                fig.subplots_adjust(left=min(0.22, cur + (0.04 - left)))
                fig.canvas.draw()
                xmin, ymin, bw, bh = _heatmap_data_bbox_in_figure(ax, nrows, ncols)
                right = xmin - pad
                left = right - cw
            cax = fig.add_axes([max(0.015, left), ymin, cw, bh])
        else:
            left = xmin + bw + pad
            if left + cw > 0.995:
                fig.subplots_adjust(right=0.99)
                fig.canvas.draw()
                xmin, ymin, bw, bh = _heatmap_data_bbox_in_figure(ax, nrows, ncols)
                left = xmin + bw + pad
            cax = fig.add_axes([left, ymin, cw, bh])
    else:
        ch = max(0.007, min(0.032, 0.022))
        hpad = 0.022
        bottom = 0.0
        for _ in range(10):
            lo = _axes_lower_edge_figure(fig, ax)
            bottom = lo - hpad - ch
            if bottom >= 0.02:
                break
            cur = float(fig.subplotpars.bottom)
            fig.subplots_adjust(bottom=min(0.45, cur + 0.07))
            fig.canvas.draw()
            xmin, ymin, bw, bh = _heatmap_data_bbox_in_figure(ax, nrows, ncols)
        if bottom < 0.01:
            bottom = 0.01
        cax = fig.add_axes([xmin, bottom, bw, ch])

    cbar = fig.colorbar(mappable, cax=cax, orientation=orient, extend="neither", label=label)
    if ticks is not None and len(ticks) > 0:
        _ticks = [float(x) for x in ticks]
        cbar.set_ticks(_ticks)
    else:
        _set_heatmap_colorbar_ticks(cbar, norm)
    if hasattr(cbar, "solids") and cbar.solids is not None:
        cbar.solids.set_rasterized(True)
    for p in cbar.ax.patches:
        if hasattr(p, "set_rasterized"):
            p.set_rasterized(True)
    return cbar


def _savefig_vector_heatmap(fig, tmp_path: str, fmt: str, dpi: int = 150) -> None:
    """SVG/PDF: editable text; keep heatmap cells vector; leave colour bar raster (single image)."""
    import matplotlib as mpl

    if fig.axes:
        main_ax = fig.axes[0]
        for artist in main_ax.get_children():
            if hasattr(artist, "set_rasterized"):
                try:
                    artist.set_rasterized(False)
                except (AttributeError, TypeError):
                    pass
    rc_extra: dict = {}
    if fmt == "svg":
        rc_extra["svg.fonttype"] = "none"
    with mpl.rc_context(rc_extra):
        fig.savefig(tmp_path, format=fmt, dpi=dpi, bbox_inches="tight")


def plot_heatmap(
    ids: list[str],
    matrix: list[list[float]],
    out_path: str,
    title: str = "RMSD heatmap",
    cmap: str = "viridis_r",
    vmin: float | None = None,
    vmax: float | None = None,
    short_labels: bool = False,
    diverging_center: str | None = None,
    vcenter: float | None = None,
    colorbar_orientation: str = "vertical",
    y_axis_right: bool = False,
    *,
    cbar_label: str = "RMSD (Å)",
    cbar_ticks: list[float] | None = None,
    cbar_tick_step: float | None = None,
    autoscale_positive_offdiag_only: bool = True,
    hatch_value_matrix: list[list[float | int | None]] | None = None,
    hatch_n_equal_bins: int = 4,
    hatch_bin_edges: list[float] | None = None,
    hatch_channel_name: str = "Hatch values",
    hatch_repeat_cells: int = _HATCH_REPEAT_CELLS,
    hatch_repeat_legend: int = _HATCH_REPEAT_LEGEND,
    hatch_linewidth_pt: float = _HATCH_LINEWIDTH_PT,
    hatch_color_light: str = "black",
    hatch_color_dark: str = "white",
    hatch_color_luma_threshold: float = 0.45,
    hatch_draw_patterns: bool = True,
    annotate: str | None = None,
    annotate_fmt: str = "{:.0f}",
    annotate_fontsize: float = 6.0,
    annotate_color_light: str = "black",
    annotate_color_dark: str = "white",
    annotate_color_luma_threshold: float = 0.45,
    annotate_hatch_fmt: str = "{:.0f}",
    annotate_hatch_fontsize: float | None = None,
    cell_border_linewidth_pt: float = _CELL_BORDER_LINEWIDTH_PT,
    cell_border_color: str = _CELL_BORDER_COLOR,
    triangle: str | None = None,
) -> None:
    """
    Draw heatmap with matplotlib and save to out_path as an image (PNG/SVG/PDF).

    RMSD tables: leave cbar_label and autoscale_positive_offdiag_only at defaults.
    Other square matrices (e.g. Dali Z): set cbar_label (e.g. \"Dali Z\") and
    autoscale_positive_offdiag_only=False so autoscale and median use all finite off-diagonals.

    Optional **second channel** (hatch patterns): if ``hatch_value_matrix`` is set (same shape
    as ``matrix``; diagonal ignored) and ``hatch_draw_patterns`` is True, each off-diagonal cell
    uses a hatch pattern for binned values from that matrix (e.g. Dali ``n_core``). Colour still
    encodes the main ``matrix``. Set ``hatch_draw_patterns=False`` to omit patterns while still
    passing ``hatch_value_matrix`` for ``annotate='main+hatch'`` (numeric n_core in cells only).

    ``annotate`` may be ``\"main+hatch\"`` (aliases: ``main_hatch``, ``both``): two-line cell labels
    from the main matrix and ``hatch_value_matrix`` (e.g. Dali Z + n_core), using
    ``annotate_hatch_fmt`` / ``annotate_hatch_fontsize`` for the second line. Requires
    ``hatch_value_matrix``.
    """
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        import numpy as np
    except ImportError as e:
        print("Matplotlib is required for heatmap output. Install with: pip install matplotlib", file=sys.stderr)
        raise SystemExit(1) from e

    out_path = os.path.abspath(out_path)
    plot_dir = os.path.dirname(out_path)
    if plot_dir and not os.path.isdir(plot_dir):
        os.makedirs(plot_dir, exist_ok=True)

    n = len(ids)
    if n < 2:
        print("Need at least 2 models to plot a heatmap.", file=sys.stderr)
        return

    # Build float array from matrix (handle None, nan, or mixed types)
    arr = np.zeros((n, n), dtype=np.float64)
    for i in range(n):
        for j in range(n):
            v = matrix[i][j]
            if v is None or (isinstance(v, float) and (v != v)):
                arr[i, j] = np.nan
            else:
                arr[i, j] = float(v)

    # Default colour scale from finite off-diagonal values (RMSD: prefer positive only).
    mask = ~np.eye(n, dtype=bool)
    off = arr[mask]
    off_finite = off[np.isfinite(off)]
    if autoscale_positive_offdiag_only:
        pos = off_finite[off_finite > 0]
        if len(pos) > 0:
            data_min = float(np.min(pos))
            data_max = float(np.max(pos))
        elif len(off_finite) > 0:
            data_min = float(np.min(off_finite))
            data_max = float(np.max(off_finite))
        else:
            print("No finite off-diagonal values to plot.", file=sys.stderr)
            return
    elif len(off_finite) > 0:
        data_min = float(np.min(off_finite))
        data_max = float(np.max(off_finite))
    else:
        print("No finite off-diagonal values to plot.", file=sys.stderr)
        return
    disp_vmin = float(vmin) if vmin is not None else data_min
    disp_vmax = float(vmax) if vmax is not None else data_max
    if disp_vmin > disp_vmax:
        print(
            "Error: heatmap scale requires vmin < vmax (got vmin=%s vmax=%s)."
            % (disp_vmin, disp_vmax),
            file=sys.stderr,
        )
        raise SystemExit(1)
    if disp_vmax - disp_vmin <= 1e-15:
        mid = 0.5 * (disp_vmin + disp_vmax)
        eps = 1e-6 * max(abs(mid), 1.0)
        disp_vmin, disp_vmax = mid - eps, mid + eps

    dc = diverging_center if diverging_center else "none"
    vc_for_norm: float | None = None
    if dc == "median":
        if vcenter is not None:
            vc_for_norm = vcenter
        elif autoscale_positive_offdiag_only:
            vc_for_norm = median_rmsd_from_square_array(arr)
        else:
            vc_for_norm = median_finite_offdiag_square_array(arr)
        if vc_for_norm is None:
            print("Warning: could not compute median; using linear colour scale.", file=sys.stderr)
            dc = "none"
        elif vcenter is None:
            print("Heatmap: median centre = %.4g" % vc_for_norm, file=sys.stderr)
    norm = make_heatmap_norm(disp_vmin, disp_vmax, dc, vc_for_norm if dc == "median" else None)

    user_ticks: list[float] | None = None
    if cbar_ticks is not None and len(cbar_ticks) > 0:
        user_ticks = [float(x) for x in cbar_ticks]
    elif cbar_tick_step is not None:
        step = float(cbar_tick_step)
        if step > 0:
            import numpy as np

            lo = float(norm.vmin)
            hi = float(norm.vmax)
            start = np.floor(lo / step) * step
            vals = np.arange(start, hi + 0.5 * step, step, dtype=float)
            vals = vals[(vals >= lo - 1e-9) & (vals <= hi + 1e-9)]
            user_ticks = [float(x) for x in vals]
        else:
            print(f"Warning: cbar_tick_step must be > 0 (got {cbar_tick_step!r}); ignoring.", file=sys.stderr)

    # Diagonal: below scale so it maps to the colour map "under" colour (white)
    np.fill_diagonal(arr, disp_vmin - 1.0)

    try:
        cmap_obj = plt.get_cmap(cmap).copy()
    except ValueError:
        print("Unknown colour map '%s'. Common: viridis_r, plasma_r, RdYlBu_r, coolwarm, YlOrRd." % cmap, file=sys.stderr)
        cmap_obj = plt.get_cmap("viridis_r").copy()
    cmap_obj.set_under(color="white")

    ext = os.path.splitext(out_path)[1].lower()
    fmt = "png" if ext == ".png" else "pdf" if ext == ".pdf" else "svg" if ext == ".svg" else "png"

    hatch_arr: np.ndarray | None = None
    hatch_full_edges: list[float] | None = None
    hatch_draw_patterns = bool(hatch_draw_patterns)
    if hatch_value_matrix is not None:
        if len(hatch_value_matrix) != n or any(len(row) != n for row in hatch_value_matrix):
            print("Error: hatch_value_matrix must match matrix shape (%d x %d)." % (n, n), file=sys.stderr)
            raise SystemExit(1)
        hatch_arr = np.zeros((n, n), dtype=np.float64)
        for i in range(n):
            for j in range(n):
                v = hatch_value_matrix[i][j]
                if v is None or i == j:
                    hatch_arr[i, j] = np.nan
                else:
                    hatch_arr[i, j] = float(v)

    ann = (annotate or "").strip().lower()
    if ann in ("", "none", "off", "false", "no"):
        ann = ""
    elif ann in ("main", "matrix", "values"):
        ann = "main"
    elif ann in ("main+hatch", "main_hatch", "both", "z_ncore", "main_and_hatch"):
        ann = "main+hatch"
    elif ann in ("hatch", "second", "second_channel", "second-channel"):
        ann = "hatch"
    else:
        print(f"Warning: unknown annotate mode {annotate!r}; using none.", file=sys.stderr)
        ann = ""

    annotate_arr = None
    annotate_hatch_stored = None
    if ann == "main":
        annotate_arr = arr.copy()
    elif ann == "main+hatch":
        if hatch_arr is None:
            print(
                "Error: annotate 'main+hatch' requires hatch_value_matrix (e.g. n_core alongside Dali Z).",
                file=sys.stderr,
            )
            raise SystemExit(1)
        annotate_arr = arr.copy()
        annotate_hatch_stored = hatch_arr.copy()
    elif ann == "hatch":
        annotate_arr = hatch_arr.copy() if hatch_arr is not None else None
        if annotate_arr is None:
            print("Warning: annotate='hatch' requested but no hatch_value_matrix provided; skipping annotations.", file=sys.stderr)
            ann = ""

    # Hatch patterns: optional (hatch_draw_patterns); never together with annotate='hatch' (values only).
    want_hatch_patterns = (
        hatch_draw_patterns and hatch_arr is not None and ann != "hatch"
    )
    if want_hatch_patterns:
        if hatch_bin_edges is not None:
            fe = sorted({float(x) for x in hatch_bin_edges})
            if len(fe) < 2:
                print("Error: hatch_bin_edges needs at least two values (bin boundaries).", file=sys.stderr)
                raise SystemExit(1)
            m2 = ~np.eye(hatch_arr.shape[0], dtype=bool)
            v2 = hatch_arr[m2]
            v2 = v2[np.isfinite(v2)]
            dlo, dhi = float(np.min(v2)), float(np.max(v2))
            if fe[0] > dlo:
                fe = [dlo] + fe
            if fe[-1] < dhi:
                fe = fe + [dhi]
            hatch_full_edges = sorted(fe)
        else:
            hatch_full_edges = _hatch_bin_edges_from_values(
                hatch_arr, max(1, int(hatch_n_equal_bins))
            )
        if hatch_full_edges is None:
            print("Warning: could not bin hatch values; hatch patterns omitted.", file=sys.stderr)
            want_hatch_patterns = False

    hatch_patterns_arr = hatch_arr if want_hatch_patterns else None
    hatch_patterns_edges = hatch_full_edges if want_hatch_patterns else None
    if ann == "hatch" and hatch_arr is not None:
        print("Heatmap: annotate hatch values (patterns disabled).", file=sys.stderr)

    tri = (triangle or "").strip().lower()
    if tri in ("lower", "lower_triangle", "lower-triangle", "lowertri", "lowertriangular"):
        tri = "lower"
    elif tri in ("", "none", "full"):
        tri = ""
    else:
        print(f"Warning: unknown triangle mode {triangle!r}; using full matrix.", file=sys.stderr)
        tri = ""

    if tri == "lower":
        # Keep diagonal and the region under it (bottom-right, i <= j); hide i > j.
        il = np.tril_indices(n, k=-1)
        arr[il] = np.nan
        if hatch_arr is not None:
            hatch_arr[il] = np.nan
        if hatch_patterns_arr is not None:
            hatch_patterns_arr[il] = np.nan
        if annotate_arr is not None:
            annotate_arr[il] = np.nan
        if annotate_hatch_stored is not None:
            annotate_hatch_stored[il] = np.nan

    use_patches = _heatmap_vector_format(fmt) or (hatch_patterns_arr is not None)

    cb_orient = colorbar_orientation.lower() if colorbar_orientation else "vertical"
    if cb_orient not in ("vertical", "horizontal"):
        cb_orient = "vertical"
    fig_h = max(5, n * 0.35)
    if cb_orient == "horizontal":
        fig_h += 1.25
    if hatch_patterns_arr is not None:
        fig_h = max(fig_h, n * 0.35 + 0.5)
    fig_w = max(6, n * 0.4)
    if hatch_patterns_arr is not None:
        fig_w += 1.8
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))
    im = None
    import matplotlib as mpl
    rc_hatch = {"hatch.linewidth": float(hatch_linewidth_pt)}
    with mpl.rc_context(rc_hatch):
        if use_patches:
            if hatch_patterns_arr is not None and hatch_patterns_edges is not None:
                _draw_heatmap_vector_cells(
                    ax,
                    arr,
                    cmap_obj,
                    norm,
                    hatch_value_arr=hatch_patterns_arr,
                    hatch_bin_edges=hatch_patterns_edges,
                    hatch_channel_name=hatch_channel_name,
                    hatch_repeat_cells=hatch_repeat_cells,
                    hatch_repeat_legend=hatch_repeat_legend,
                    hatch_linewidth_pt=hatch_linewidth_pt,
                    hatch_color_light=hatch_color_light,
                    hatch_color_dark=hatch_color_dark,
                    hatch_color_luma_threshold=hatch_color_luma_threshold,
                    annotate_value_arr=annotate_arr,
                    annotate_fmt=annotate_fmt,
                    annotate_fontsize=annotate_fontsize,
                    annotate_color_light=annotate_color_light,
                    annotate_color_dark=annotate_color_dark,
                    annotate_color_luma_threshold=annotate_color_luma_threshold,
                    annotate_hatch_value_arr=annotate_hatch_stored,
                    annotate_hatch_fmt=annotate_hatch_fmt,
                    annotate_hatch_fontsize=annotate_hatch_fontsize,
                    cell_border_linewidth_pt=cell_border_linewidth_pt,
                    cell_border_color=cell_border_color,
                    triangle=tri or None,
                )
            else:
                _draw_heatmap_vector_cells(
                    ax,
                    arr,
                    cmap_obj,
                    norm,
                    annotate_value_arr=annotate_arr,
                    annotate_fmt=annotate_fmt,
                    annotate_fontsize=annotate_fontsize,
                    annotate_color_light=annotate_color_light,
                    annotate_color_dark=annotate_color_dark,
                    annotate_color_luma_threshold=annotate_color_luma_threshold,
                    annotate_hatch_value_arr=annotate_hatch_stored,
                    annotate_hatch_fmt=annotate_hatch_fmt,
                    annotate_hatch_fontsize=annotate_hatch_fontsize,
                    cell_border_linewidth_pt=cell_border_linewidth_pt,
                    cell_border_color=cell_border_color,
                    triangle=tri or None,
                )
        else:
            im = ax.imshow(
                arr,
                aspect="equal",
                origin="lower",
                cmap=cmap_obj,
                interpolation="nearest",
                norm=norm,
            )
            # Add per-cell borders for raster heatmaps.
            ax.set_xticks(np.arange(-0.5, n, 1), minor=True)
            ax.set_yticks(np.arange(-0.5, n, 1), minor=True)
            ax.grid(
                which="minor",
                color=cell_border_color,
                linestyle="-",
                linewidth=float(cell_border_linewidth_pt),
            )
            ax.tick_params(which="minor", bottom=False, left=False)
            if annotate_arr is not None:
                rafs2 = (
                    float(annotate_hatch_fontsize)
                    if annotate_hatch_fontsize is not None
                    else max(3.25, float(annotate_fontsize) * 0.88)
                )
                # For raster heatmaps, add text labels after imshow.
                for i in range(n):
                    for j in range(n):
                        if tri == "lower" and i > j:
                            continue
                        if i == j:
                            continue
                        av = annotate_arr[i, j]
                        if annotate_hatch_stored is not None:
                            aha = annotate_hatch_stored[i, j]
                            if not (np.isfinite(av) or np.isfinite(aha)):
                                continue
                        elif not np.isfinite(av):
                            continue
                        # Determine face colour for contrast.
                        v = arr[i, j]
                        if np.isnan(v):
                            continue
                        if v < float(norm.vmin):
                            fc = cmap_obj.get_under()
                            rgba = (1.0, 1.0, 1.0, 1.0) if fc is None else plt.matplotlib.colors.to_rgba(fc)
                        else:
                            rgba = cmap_obj(norm(v))
                        luma = _relative_luminance_srgb(rgba)
                        tc = annotate_color_dark if luma <= float(annotate_color_luma_threshold) else annotate_color_light
                        if annotate_hatch_stored is not None:
                            aha = annotate_hatch_stored[i, j]
                            if np.isfinite(av) and np.isfinite(aha):
                                ttxt = (
                                    f"{_format_cell_value(float(av), annotate_fmt)}\n"
                                    f"{_format_cell_value(float(aha), annotate_hatch_fmt)}"
                                )
                            elif np.isfinite(av):
                                ttxt = _format_cell_value(float(av), annotate_fmt)
                            else:
                                ttxt = _format_cell_value(float(aha), annotate_hatch_fmt)
                            ax.text(
                                j,
                                i,
                                ttxt,
                                ha="center",
                                va="center",
                                fontsize=rafs2,
                                linespacing=0.92,
                                color=tc,
                                clip_on=True,
                            )
                        else:
                            ax.text(
                                j,
                                i,
                                _format_cell_value(float(av), annotate_fmt),
                                ha="center",
                                va="center",
                                fontsize=float(annotate_fontsize),
                                color=tc,
                                clip_on=True,
                            )
    ax.set_xticks(range(n))
    ax.set_yticks(range(n))
    tick_labels = [short_heatmap_label(x) if short_labels else x for x in ids]
    ax.set_xticklabels(tick_labels, rotation=45, ha="right")
    ax.set_yticklabels(tick_labels)
    ax.set_title(title)
    if y_axis_right:
        _apply_heatmap_y_axis_right(ax)
    cbar_left = y_axis_right and cb_orient == "vertical"
    if cb_orient == "horizontal":
        r_edge = 0.86 if y_axis_right else 0.94
        if hatch_patterns_arr is not None:
            r_edge = max(0.70, r_edge - 0.08)
        plt.tight_layout(rect=(0.06, 0.22, r_edge, 0.94))
    elif cbar_left:
        plt.tight_layout(rect=(0.16, 0.06, 0.90, 0.92))
    else:
        r_layout = 0.78 if hatch_patterns_arr is not None else 0.88
        plt.tight_layout(rect=(0.08, 0.06, r_layout, 0.92))
    if _heatmap_vector_format(fmt):
        _add_raster_colorbar(
            ax,
            cmap_obj,
            norm,
            nrows=n,
            ncols=n,
            label=cbar_label,
            orientation=cb_orient,
            mappable=None,
            vertical_colorbar_on_left=cbar_left,
            ticks=user_ticks,
        )
    else:
        _add_raster_colorbar(
            ax,
            cmap_obj,
            norm,
            nrows=n,
            ncols=n,
            label=cbar_label,
            orientation=cb_orient,
            mappable=im,
            vertical_colorbar_on_left=cbar_left,
            ticks=user_ticks,
        )

    # Save to a temp file in the same directory, then move to target (avoids partial/CSV overwrite)
    fd, tmp_path = tempfile.mkstemp(suffix=ext, prefix="rmsd_heatmap_", dir=plot_dir or ".")
    try:
        os.close(fd)
        try:
            if _heatmap_vector_format(fmt):
                _savefig_vector_heatmap(fig, tmp_path, fmt, dpi=150)
            else:
                fig.savefig(tmp_path, format=fmt, dpi=150, bbox_inches="tight")
        except Exception as e:
            plt.close(fig)
            if os.path.exists(tmp_path):
                try:
                    os.remove(tmp_path)
                except OSError:
                    pass
            print("Error saving figure:", e, file=sys.stderr)
            if out_path.lower().endswith((".pdf", ".svg")):
                print("Try PNG instead: --plot", os.path.splitext(out_path)[0] + ".png", file=sys.stderr)
            raise SystemExit(1) from e
        plt.close(fig)

        size = os.path.getsize(tmp_path) if os.path.exists(tmp_path) else 0
        min_ok = 3000 if fmt in ("svg", "pdf") else 4000
        if size < min_ok or not _image_magic_ok(tmp_path, fmt):
            with open(tmp_path, "rb") as f:
                preview = f.read(200)
            try:
                os.remove(tmp_path)
            except OSError:
                pass
            print("Heatmap save produced a too-small or invalid file (%d bytes). First bytes: %r" % (size, preview[:80]), file=sys.stderr)
            if b"Model" in preview or (b"," in preview and size < 5000):
                print("The file looks like CSV. Do not use the same path for -o and --plot (e.g. use -o out.csv --plot heatmap.png).", file=sys.stderr)
            else:
                print("Try saving with a simple test: python -c \"import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt; import numpy as np; plt.imshow(np.random.rand(10,10)); plt.savefig('test_heatmap.png'); print(open('test_heatmap.png','rb').read()[:8])\"", file=sys.stderr)
            raise SystemExit(1)

        os.replace(tmp_path, out_path)
        print("Wrote", out_path)
    except Exception:
        if os.path.exists(tmp_path):
            try:
                os.remove(tmp_path)
            except OSError:
                pass
        raise


def _load_square_matrix_csv(csv_path: str) -> tuple[list[str], list[list[float | None]]]:
    """
    Load a square matrix CSV with header: Model,<id1>,<id2>,... and one row per id.
    Cells may be '-', '', 'NA', 'NaN' for missing.
    Returns (ids, matrix) where matrix is list[list[float|None]].
    """
    import csv

    path = os.path.abspath(csv_path)
    with open(path, newline="") as f:
        r = csv.reader(f)
        header = next(r, None)
        if not header or len(header) < 2:
            raise ValueError(f"matrix CSV has no header or too few columns: {path}")
        ids = [x.strip() for x in header[1:] if str(x).strip()]
        if not ids:
            raise ValueError(f"matrix CSV has empty label header: {path}")
        rows: dict[str, list[str]] = {}
        for row in r:
            if not row:
                continue
            key = (row[0] or "").strip()
            if not key:
                continue
            rows[key] = [str(x) for x in row[1:]]
    missing = [x for x in ids if x not in rows]
    if missing:
        raise ValueError(f"matrix CSV missing rows for {len(missing)} label(s) (e.g. {missing[:3]}): {path}")
    n = len(ids)
    mat: list[list[float | None]] = []
    for la in ids:
        vals = rows[la]
        if len(vals) < n:
            vals = vals + ["-"] * (n - len(vals))
        use: list[float | None] = []
        for s in vals[:n]:
            ss = (s or "").strip()
            if ss in ("", "-", "NA", "NaN", "nan", "N/A"):
                use.append(None)
            else:
                use.append(float(ss))
        mat.append(use)
    return ids, mat


def _load_pair_table_hatch_matrix(
    pairs_csv: str,
    ids: list[str],
    *,
    a_col: str = "label_a",
    b_col: str = "label_b",
    value_col: str = "n_core",
) -> list[list[float | int | None]]:
    """Build an NxN hatch matrix from a pairwise CSV keyed by (label_a,label_b)."""
    import csv

    path = os.path.abspath(pairs_csv)
    ncore: dict[tuple[str, str], float | int] = {}
    with open(path, newline="") as f:
        r = csv.DictReader(f)
        for row in r:
            a = (row.get(a_col) or "").strip()
            b = (row.get(b_col) or "").strip()
            v = row.get(value_col)
            if not a or not b or v is None:
                continue
            try:
                ncore[(a, b)] = float(v)
            except ValueError:
                continue
    mat: list[list[float | int | None]] = []
    for i, la in enumerate(ids):
        row: list[float | int | None] = []
        for j, lb in enumerate(ids):
            if i == j:
                row.append(None)
            else:
                val = ncore.get((la, lb), ncore.get((lb, la)))
                row.append(val)
        mat.append(row)
    return mat


def _load_pair_table_square_matrix(
    pairs_csv: str,
    *,
    a_col: str = "label_a",
    b_col: str = "label_b",
    value_col: str = "pct_id",
    ids: list[str] | None = None,
) -> tuple[list[str], list[list[float | None]]]:
    """
    Build a square matrix from a pairwise CSV keyed by (label_a,label_b).
    If ids is None, labels are derived from the union of a_col/b_col and natural-sorted.
    Diagonal entries are set to None.
    """
    import csv

    path = os.path.abspath(pairs_csv)
    vals: dict[tuple[str, str], float] = {}
    labels: set[str] = set()
    with open(path, newline="") as f:
        r = csv.DictReader(f)
        for row in r:
            a = (row.get(a_col) or "").strip()
            b = (row.get(b_col) or "").strip()
            v = row.get(value_col)
            if not a or not b:
                continue
            labels.add(a)
            labels.add(b)
            if v is None:
                continue
            try:
                vals[(a, b)] = float(v)
            except ValueError:
                continue

    use_ids = ids if ids is not None else sorted(labels, key=_natural_sort_key)
    n = len(use_ids)
    mat: list[list[float | None]] = []
    for i, la in enumerate(use_ids):
        row: list[float | None] = []
        for j, lb in enumerate(use_ids):
            if i == j:
                row.append(None)
            else:
                row.append(vals.get((la, lb), vals.get((lb, la))))
        mat.append(row)
    return use_ids, mat


def _parse_float_list_csv(s: str) -> list[float] | None:
    if s is None:
        return None
    parts = [p.strip() for p in str(s).split(",") if p.strip()]
    if len(parts) < 2:
        return None
    return [float(p) for p in parts]


def main(argv: list[str] | None = None) -> None:
    ap = argparse.ArgumentParser(
        description="Plot a square matrix CSV as a heatmap (optional hatch second channel).",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    mg = ap.add_mutually_exclusive_group(required=True)
    mg.add_argument(
        "--matrix",
        metavar="CSV",
        help="Square matrix CSV with header 'Model,<id1>,<id2>,…' and one row per id.",
    )
    mg.add_argument(
        "--matrix-pairs",
        metavar="CSV",
        help="Pairwise CSV table to turn into a square matrix (uses --matrix-pairs-*-col).",
    )
    ap.add_argument(
        "--out",
        required=True,
        metavar="PATH",
        help="Output image path (.png, .svg, or .pdf).",
    )
    ap.add_argument("--title", default="Heatmap", help="Figure title.")
    ap.add_argument("--cmap", default="viridis_r", help="Matplotlib colormap name.")
    ap.add_argument("--vmin", type=float, default=None, help="Colour scale floor (data units).")
    ap.add_argument("--vmax", type=float, default=None, help="Colour scale ceiling (data units).")
    ap.add_argument(
        "--diverging-center",
        choices=("none", "median"),
        default="none",
        help="Use TwoSlopeNorm centred at the median (good for diverging colormaps).",
    )
    ap.add_argument(
        "--short-labels",
        action="store_true",
        help="Shorten axis labels (FoldKit convention) for the figure only.",
    )
    ap.add_argument(
        "--colorbar-orientation",
        choices=("vertical", "horizontal"),
        default="vertical",
        help="Colour bar layout.",
    )
    ap.add_argument(
        "--y-axis-right",
        action="store_true",
        help="Draw row tick labels on the right (useful for long labels).",
    )
    ap.add_argument(
        "--cbar-label",
        default="Value",
        help="Colour bar label.",
    )
    ap.add_argument(
        "--cbar-ticks",
        default=None,
        metavar="T0,T1,…",
        help="Optional comma-separated explicit colour bar ticks in data units.",
    )
    ap.add_argument(
        "--cbar-tick-step",
        type=float,
        default=None,
        metavar="STEP",
        help="Optional colour bar tick spacing in data units (e.g. 5 for 50,55,...).",
    )
    ap.add_argument(
        "--matrix-pairs-a-col",
        default="label_a",
        help="Column name for the first label in --matrix-pairs.",
    )
    ap.add_argument(
        "--matrix-pairs-b-col",
        default="label_b",
        help="Column name for the second label in --matrix-pairs.",
    )
    ap.add_argument(
        "--matrix-pairs-value-col",
        default="pct_id",
        help="Column name for the matrix value in --matrix-pairs (e.g. pct_id, z_score, n_core).",
    )
    ap.add_argument(
        "--matrix-pairs-order-from",
        default=None,
        metavar="CSV",
        help="Optional: square matrix CSV whose label order should be reused when building from --matrix-pairs.",
    )
    ap.add_argument(
        "--triangle",
        choices=("full", "lower"),
        default="full",
        help="Plot only part of the matrix: full, or lower (under diagonal, including diagonal).",
    )
    ap.add_argument(
        "--cell-border-linewidth-pt",
        type=float,
        default=_CELL_BORDER_LINEWIDTH_PT,
        help="Per-cell border line width (points).",
    )
    ap.add_argument(
        "--cell-border-color",
        default=_CELL_BORDER_COLOR,
        help="Per-cell border colour.",
    )
    ap.add_argument(
        "--autoscale-positive-offdiag-only",
        action="store_true",
        help="Autoscale from positive off-diagonal values only (RMSD-style). Default is to use all finite off-diagonal values.",
    )
    ap.add_argument(
        "--annotate",
        choices=("none", "main", "hatch", "main+hatch"),
        default="none",
        help="Cell labels: main (colour matrix), hatch (second channel), or main+hatch (two lines; requires --hatch-*).",
    )
    ap.add_argument(
        "--annotate-fmt",
        default="{:.0f}",
        help="Python format string for annotations, e.g. '{:.0f}', '{:.1f}', '{:.2g}'.",
    )
    ap.add_argument(
        "--annotate-fontsize",
        type=float,
        default=6.0,
        help="Annotation font size (points).",
    )
    ap.add_argument(
        "--annotate-color-light",
        default="black",
        help="Annotation colour on light cells.",
    )
    ap.add_argument(
        "--annotate-color-dark",
        default="white",
        help="Annotation colour on dark cells (for contrast).",
    )
    ap.add_argument(
        "--annotate-color-threshold",
        type=float,
        default=0.45,
        metavar="LUMA",
        help="Relative luminance threshold in [0,1]; cells darker than this use --annotate-color-dark.",
    )
    ap.add_argument(
        "--annotate-hatch-fmt",
        default="{:.0f}",
        help="With --annotate main+hatch: format for the second line (hatch / n_core).",
    )

    hatch = ap.add_argument_group("Hatch (second channel)")
    hatch.add_argument(
        "--hatch-matrix",
        default=None,
        metavar="CSV",
        help="Optional: square matrix CSV (same labels) for hatch values.",
    )
    hatch.add_argument(
        "--hatch-pairs",
        default=None,
        metavar="CSV",
        help="Optional: pairwise table CSV (e.g. dali_pairs_*.csv) with label_a/label_b and a value column.",
    )
    hatch.add_argument(
        "--hatch-pairs-a-col",
        default="label_a",
        help="Column name for the first label in --hatch-pairs.",
    )
    hatch.add_argument(
        "--hatch-pairs-b-col",
        default="label_b",
        help="Column name for the second label in --hatch-pairs.",
    )
    hatch.add_argument(
        "--hatch-pairs-value-col",
        default="n_core",
        help="Column name for the hatch value in --hatch-pairs (e.g. n_core).",
    )
    hatch.add_argument("--hatch-name", default="Hatch values", help="Legend title for hatch values.")
    hatch.add_argument("--hatch-bins", type=int, default=4, help="Number of equal-width bins for hatch values.")
    hatch.add_argument(
        "--hatch-edges",
        default=None,
        metavar="E0,E1,…",
        help="Optional comma-separated bin boundaries for hatch values; extended to data min/max if needed.",
    )
    hatch.add_argument(
        "--hatch-repeat-cells",
        type=int,
        default=_HATCH_REPEAT_CELLS,
        help="Hatch density in cells (smaller = more spaced lines).",
    )
    hatch.add_argument(
        "--hatch-repeat-legend",
        type=int,
        default=_HATCH_REPEAT_LEGEND,
        help="Hatch density in legend patches.",
    )
    hatch.add_argument(
        "--hatch-linewidth-pt",
        type=float,
        default=_HATCH_LINEWIDTH_PT,
        help="Hatch line thickness (points).",
    )
    hatch.add_argument(
        "--hatch-color-light",
        default="black",
        help="Hatch colour on light cells.",
    )
    hatch.add_argument(
        "--hatch-color-dark",
        default="white",
        help="Hatch colour on dark cells (for contrast).",
    )
    hatch.add_argument(
        "--hatch-color-threshold",
        type=float,
        default=0.45,
        metavar="LUMA",
        help="Relative luminance threshold in [0,1]; cells darker than this use --hatch-color-dark.",
    )

    args = ap.parse_args(argv)

    if args.matrix:
        ids, matrix = _load_square_matrix_csv(args.matrix)
    else:
        ids_order: list[str] | None = None
        if args.matrix_pairs_order_from:
            ids_order, _m = _load_square_matrix_csv(args.matrix_pairs_order_from)
        ids, matrix = _load_pair_table_square_matrix(
            args.matrix_pairs,
            a_col=args.matrix_pairs_a_col,
            b_col=args.matrix_pairs_b_col,
            value_col=args.matrix_pairs_value_col,
            ids=ids_order,
        )

    hatch_value_matrix = None
    if args.hatch_matrix and args.hatch_pairs:
        ap.error("Use only one of --hatch-matrix or --hatch-pairs.")
    if args.hatch_matrix:
        h_ids, hatch_value_matrix = _load_square_matrix_csv(args.hatch_matrix)
        if h_ids != ids:
            raise SystemExit("Error: --hatch-matrix labels do not match --matrix labels (same order required).")
    elif args.hatch_pairs:
        hatch_value_matrix = _load_pair_table_hatch_matrix(
            args.hatch_pairs,
            ids,
            a_col=args.hatch_pairs_a_col,
            b_col=args.hatch_pairs_b_col,
            value_col=args.hatch_pairs_value_col,
        )

    hatch_edges = _parse_float_list_csv(args.hatch_edges) if args.hatch_edges else None

    cbar_ticks = _parse_float_list_csv(args.cbar_ticks) if args.cbar_ticks else None
    plot_heatmap(
        ids,
        matrix,
        args.out,
        title=args.title,
        cmap=args.cmap,
        vmin=args.vmin,
        vmax=args.vmax,
        short_labels=bool(args.short_labels),
        diverging_center=args.diverging_center if args.diverging_center != "none" else "none",
        vcenter=None,
        colorbar_orientation=args.colorbar_orientation,
        y_axis_right=bool(args.y_axis_right),
        cbar_label=args.cbar_label,
        cbar_ticks=cbar_ticks,
        cbar_tick_step=args.cbar_tick_step,
        triangle=args.triangle,
        cell_border_linewidth_pt=float(args.cell_border_linewidth_pt),
        cell_border_color=str(args.cell_border_color),
        autoscale_positive_offdiag_only=bool(args.autoscale_positive_offdiag_only),
        annotate=args.annotate,
        annotate_fmt=str(args.annotate_fmt),
        annotate_fontsize=float(args.annotate_fontsize),
        annotate_color_light=str(args.annotate_color_light),
        annotate_color_dark=str(args.annotate_color_dark),
        annotate_color_luma_threshold=float(args.annotate_color_threshold),
        annotate_hatch_fmt=str(args.annotate_hatch_fmt),
        hatch_value_matrix=hatch_value_matrix,
        hatch_n_equal_bins=max(1, int(args.hatch_bins)),
        hatch_bin_edges=hatch_edges,
        hatch_channel_name=args.hatch_name,
        hatch_repeat_cells=max(1, int(args.hatch_repeat_cells)),
        hatch_repeat_legend=max(1, int(args.hatch_repeat_legend)),
        hatch_linewidth_pt=float(args.hatch_linewidth_pt),
        hatch_color_light=str(args.hatch_color_light),
        hatch_color_dark=str(args.hatch_color_dark),
        hatch_color_luma_threshold=float(args.hatch_color_threshold),
    )


if __name__ == "__main__":
    main()
