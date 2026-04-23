#!/usr/bin/env python3
"""
Convert pairwise RMSD logs (or an existing RMSD CSV) to square CSV table(s), optional heatmaps.

Heatmap rendering lives in ``foldkit_heatmap`` (also used by ``dalilite_pairs``).

Single-file mode (pass one INPUT path; omit --scan-dir):
  Reads one LSQ text, SSM text, or square rmsd_table_*.csv (same format detection as
  structure_phylogeny). Builds a symmetric pairwise RMSD matrix with natural-sorted
  labels, writes one CSV (default: beside INPUT as rmsd_table_<suffix>.csv), and
  optionally --plot for one matplotlib heatmap. Use --order / --order-dir to reorder
  rows and columns to match PDB stems or a list.

Batch mode (--scan-dir):
  Recursively finds logs matching --scan-glob (default rmsd_values_*.txt); for each,
  writes rmsd_table_<suffix>.csv next to the log and stacks all blocks into
  combined_rmsd_table.csv with a Subdomain column (--combined-subdomain-order for block order).
  Optional --heatmap-dir writes figures
  per file plus combined_rmsd_heatmap.<fmt> (--heatmap-format png|pdf|svg, default svg;
  SVG/PDF: vector patches per matrix cell; colour bar is one embedded raster for a smooth scale (text stays editable).
  Per-subdomain heatmaps use the same colour scale as the merged CSV by default (--no-shared-heatmap-scale to use local auto scale).
  --combined-heatmap PATH for a custom combined path (extension overrides --heatmap-format);
  --no-combined-heatmap to skip the combined figure.

Examples:
  python ranking/rmsd_to_csv.py /path/to/project/rmsd_SSM_values.txt -o /path/to/project/rmsd_table.csv
  python ranking/rmsd_to_csv.py /path/to/project/rmsd_SSM_values.txt --plot /path/to/project/heatmap.pdf
  python ranking/rmsd_to_csv.py /path/to/project/rmsd_values.txt --order-dir /path/to/project/models --order-glob 'model_*.pdb'
  python ranking/rmsd_to_csv.py --scan-dir /path/to/project/run1/ --scan-glob 'rmsd_values_*.txt'
  python ranking/rmsd_to_csv.py --scan-dir /path/to/project/run1/ --heatmap-dir /path/to/project/figs

Outputs feed ranking/structure_phylogeny.py and ranking/create_rmsd_heatmap.R.
"""

from __future__ import annotations

import argparse
import csv
import glob
import os
import re
import sys
import tempfile
from pathlib import Path

_REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if _REPO_ROOT not in sys.path:
    sys.path.insert(0, _REPO_ROOT)

from cli_log import add_log_args, setup_log_from_args

from foldkit_heatmap import (
    _add_raster_colorbar,
    _apply_heatmap_y_axis_right,
    _draw_heatmap_vector_cells,
    _heatmap_vector_format,
    _image_magic_ok,
    _savefig_vector_heatmap,
    make_heatmap_norm,
    plot_heatmap,
    short_heatmap_label,
)


def _resolved_path(path: str) -> str:
    """Expand ~ and resolve to absolute path (Python does not expand ~ in abspath)."""
    return os.path.abspath(os.path.expanduser(path))


# Reuse parsing and matrix logic from structure_phylogeny
try:
    from structure_phylogeny import (
        _natural_sort_key,
        alignments_pair_counts_to_matrix,
        alignments_to_matrix,
        detect_format,
        parse_lsq_rmsd_txt,
        parse_lsq_rmsd_txt_with_pair_counts,
        parse_rmsd_csv,
        parse_ssm_rmsd_txt,
        parse_ssm_rmsd_txt_with_pair_counts,
    )
except ImportError:
    # Run from repository root (python ranking/rmsd_to_csv.py ...)
    _dir = os.path.dirname(os.path.abspath(__file__))
    if _dir not in sys.path:
        sys.path.insert(0, _dir)
    from structure_phylogeny import (
        _natural_sort_key,
        alignments_pair_counts_to_matrix,
        alignments_to_matrix,
        detect_format,
        parse_lsq_rmsd_txt,
        parse_lsq_rmsd_txt_with_pair_counts,
        parse_rmsd_csv,
        parse_ssm_rmsd_txt,
        parse_ssm_rmsd_txt_with_pair_counts,
    )


# Colour = RMSD; hatch = binned counts from Coot log (LSQ/SSM only).
_HATCH_CHANNEL_LABELS_RMSD = {
    "lsq_txt": "Matched atoms (Coot LSQ)",
    "ssm_txt": "Aligned residues (Coot SSM)",
}


def rmsd_file_to_matrix(path: str):
    """Load RMSD file and return (ids, matrix) as in structure_phylogeny."""
    fmt = detect_format(path)
    if fmt == "lsq_txt":
        alignments = parse_lsq_rmsd_txt(path)
    elif fmt == "ssm_txt":
        alignments = parse_ssm_rmsd_txt(path)
    else:
        alignments = parse_rmsd_csv(path)
    if not alignments:
        return None, None
    return alignments_to_matrix(alignments)


def rmsd_file_to_matrix_with_hatch(
    path: str,
) -> tuple[
    list[str] | None,
    list[list[float]] | None,
    str,
    list[list[float | None]] | None,
    str | None,
]:
    """
    Load RMSD and optional per-pair counts for heatmap hatching: LSQ → matched atoms;
    SSM → aligned residues. CSV input has no second channel (returns hatch None).
    """
    fmt = detect_format(path)
    counts: dict = {}
    if fmt == "lsq_txt":
        alignments, counts = parse_lsq_rmsd_txt_with_pair_counts(path)
    elif fmt == "ssm_txt":
        alignments, counts = parse_ssm_rmsd_txt_with_pair_counts(path)
    else:
        alignments = parse_rmsd_csv(path)
    if not alignments:
        return None, None, fmt, None, None
    ids, matrix = alignments_to_matrix(alignments)
    if not counts or fmt not in _HATCH_CHANNEL_LABELS_RMSD:
        return ids, matrix, fmt, None, None
    hmat = alignments_pair_counts_to_matrix(ids, counts)
    has_any = any(
        hmat[i][j] is not None
        for i in range(len(ids))
        for j in range(len(ids))
        if i != j
    )
    if not has_any:
        return ids, matrix, fmt, None, None
    return (
        ids,
        matrix,
        fmt,
        hmat,
        _HATCH_CHANNEL_LABELS_RMSD[fmt],
    )


def reorder_matrix(ids: list[str], matrix: list[list[float]], order: list[str]) -> tuple[list[str], list[list[float]]]:
    """
    Return (ordered_ids, reordered_matrix).

    order lists desired labels; any missing from matrix are appended.

    Note: LSQ/SSM logs sometimes embed extra suffixes in labels (e.g. "_rechain") while
    an --order-dir scan typically yields raw structure stems. To make ordering robust,
    we also try a conservative normalised match when an exact label is absent.
    """

    def _norm(s: str) -> str:
        s = (s or "").strip()
        # Common FoldKit pipeline suffixes that are not part of the "structure identity"
        for suf in ("_rechain", "_renamed", "_trimmed", "_trim"):
            if s.endswith(suf):
                s = s[: -len(suf)]
        return s

    id_to_idx = {x: i for i, x in enumerate(ids)}
    # Normalised id -> matching original ids (handle potential collisions)
    norm_to_ids: dict[str, list[str]] = {}
    for x in ids:
        nx = _norm(x)
        norm_to_ids.setdefault(nx, []).append(x)

    ordered_ids: list[str] = []
    for want in order:
        if want in id_to_idx:
            ordered_ids.append(want)
            continue
        nw = _norm(want)
        cands = norm_to_ids.get(nw) or []
        if len(cands) == 1:
            ordered_ids.append(cands[0])
    missing = [x for x in ids if x not in ordered_ids]
    ordered_ids = ordered_ids + missing
    n = len(ordered_ids)
    new_matrix = [[0.0] * n for _ in range(n)]
    for i, a in enumerate(ordered_ids):
        for j, b in enumerate(ordered_ids):
            new_matrix[i][j] = matrix[id_to_idx[a]][id_to_idx[b]]
    return ordered_ids, new_matrix


def reorder_hatch_matrix(
    ids: list[str],
    hatch: list[list[float | None]],
    order: list[str],
) -> tuple[list[str], list[list[float | None]]]:
    """Same index permutation as ``reorder_matrix`` for the optional hatch matrix."""

    def _norm(s: str) -> str:
        s = (s or "").strip()
        for suf in ("_rechain", "_renamed", "_trimmed", "_trim"):
            if s.endswith(suf):
                s = s[: -len(suf)]
        return s

    id_to_idx = {x: i for i, x in enumerate(ids)}
    norm_to_ids: dict[str, list[str]] = {}
    for x in ids:
        nx = _norm(x)
        norm_to_ids.setdefault(nx, []).append(x)

    ordered_ids: list[str] = []
    for want in order:
        if want in id_to_idx:
            ordered_ids.append(want)
            continue
        nw = _norm(want)
        cands = norm_to_ids.get(nw) or []
        if len(cands) == 1:
            ordered_ids.append(cands[0])
    missing = [x for x in ids if x not in ordered_ids]
    ordered_ids = ordered_ids + missing
    n = len(ordered_ids)
    new_h = [[None] * n for _ in range(n)]
    for i, a in enumerate(ordered_ids):
        for j, b in enumerate(ordered_ids):
            new_h[i][j] = hatch[id_to_idx[a]][id_to_idx[b]]
    return ordered_ids, new_h


def parse_order(order_arg: str) -> list[str]:
    """Parse --order: path to file (one label per line) or comma-separated list."""
    order_arg = order_arg.strip()
    if not order_arg:
        return []
    if os.path.isfile(order_arg):
        with open(order_arg, "r") as f:
            return [line.strip() for line in f if line.strip()]
    return [x.strip() for x in order_arg.split(",") if x.strip()]


def order_labels_from_scan(root: str, glob_pat: str) -> list[str]:
    """Recursive rglob under root; sort paths naturally; labels are file stems (unique, first path wins)."""
    root_p = Path(root).resolve()
    if not root_p.is_dir():
        return []
    paths = [p for p in root_p.rglob(glob_pat) if p.is_file()]
    paths.sort(key=lambda p: _natural_sort_key(str(p.relative_to(root_p))))
    out: list[str] = []
    seen: set[str] = set()
    for p in paths:
        stem = p.stem
        if stem not in seen:
            out.append(stem)
            seen.add(stem)
    return out


def write_rmsd_csv(ids: list[str], matrix: list[list[float]], out_path: str) -> None:
    """Write square RMSD table: Model column + one column per model, diagonal '-'."""
    out_dir = os.path.dirname(out_path)
    if out_dir and not os.path.isdir(out_dir):
        os.makedirs(out_dir, exist_ok=True)
    with open(out_path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["Model"] + ids)
        for i, name in enumerate(ids):
            row = [name]
            for j in range(len(ids)):
                if i == j:
                    row.append("-")
                else:
                    v = matrix[i][j]
                    if v is None or (isinstance(v, float) and (v != v or v < 0)):  # nan or None
                        row.append("")
                    else:
                        row.append(f"{v:.4f}")
            w.writerow(row)


def subdomain_label_from_rmsd_path(path: str) -> str:
    """Short label for rmsd_table_<label>.csv / combined Subdomain column (legacy create_rmsd_table names)."""
    base = os.path.basename(path)
    if base.startswith("rmsd_values_") and base.endswith(".txt"):
        return base[12:-4]
    if base.startswith("rmsd_SSM_values") and base.endswith(".txt"):
        rest = base[len("rmsd_SSM_values") : -4].strip("_")
        return rest if rest else "ssm"
    return os.path.splitext(base)[0]


_KNOWN_HEATMAP_EXT = frozenset({".png", ".pdf", ".svg"})


def ensure_batch_heatmap_path(path: str, heatmap_format: str) -> str:
    """If path ends with .png/.pdf/.svg, return as-is; else append .{heatmap_format}."""
    _, ext = os.path.splitext(path)
    if ext.lower() in _KNOWN_HEATMAP_EXT:
        return path
    fmt = heatmap_format.lower().lstrip(".")
    if fmt not in ("png", "pdf", "svg"):
        fmt = "png"
    return f"{path}.{fmt}"


def discover_rmsd_files(scan_root: str, glob_pattern: str) -> list[str]:
    """Recursive glob under scan_root; return sorted existing file paths."""
    root = os.path.abspath(scan_root)
    hits = glob.glob(os.path.join(root, "**", glob_pattern), recursive=True)
    return sorted(p for p in hits if os.path.isfile(p))


def _combined_subdomain_matches_token(sub: str, tok: str) -> bool:
    """
    True if subdomain label should be grouped under token tok (from --combined-subdomain-order).

    Uses a boundary after '_' so D1 does not match …_D10s. Token may be e.g. D1, CWB2, ID.
    """
    t = (tok or "").strip()
    if not t:
        return False
    sl = sub.lower()
    tl = t.lower()
    if sl == tl:
        return True
    pat = r"(^|_)" + re.escape(tl) + r"(s|_|$)"
    return re.search(pat, sl) is not None


def _combined_subdomain_block_index(sub: str, order: list[str]) -> int:
    """First matching token index, or len(order) for natural-sort tail."""
    for i, tok in enumerate(order):
        if _combined_subdomain_matches_token(sub, tok):
            return i
    return len(order)


def _combined_sorted_row_specs(
    tables: list,
    subdomain_order: list[str] | None = None,
) -> list[tuple[str, str, list[str], list[list[float]], int]]:
    """(Subdomain, model, ids, matrix, row_index) for every row, sorted by subdomain then model."""
    specs: list[tuple[str, str, list[str], list[list[float]], int]] = []
    for t in tables:
        sub, ids, matrix = t[0], t[1], t[2]
        for i, name in enumerate(ids):
            specs.append((sub, name, ids, matrix, i))
    if subdomain_order:
        specs.sort(
            key=lambda t: (
                _combined_subdomain_block_index(t[0], subdomain_order),
                _natural_sort_key(t[0]),
                _natural_sort_key(t[1]),
            )
        )
    else:
        specs.sort(key=lambda t: (_natural_sort_key(t[0]), _natural_sort_key(t[1])))
    return specs


def _combined_structure_column_order(
    specs: list[tuple[str, str, list[str], list[list[float]], int]],
    all_models: set[str],
) -> list[str]:
    """Structure columns follow row order: first time each model appears as a row (sorted rows), then any stragglers."""
    order: list[str] = []
    seen: set[str] = set()
    for _sub, name, _ids, _matrix, _i in specs:
        if name not in seen:
            seen.add(name)
            order.append(name)
    rest = all_models - seen
    if rest:
        order.extend(sorted(rest, key=_natural_sort_key))
    return order


def write_combined_rmsd_csv(
    tables: list[tuple[str, list[str], list[list[float]]]],
    out_path: str,
    subdomain_order: list[str] | None = None,
) -> None:
    """Wide CSV: rows sorted by (Subdomain, Model); structure columns match row order (first-seen model order)."""
    if not tables:
        return
    all_models = {x for t in tables for x in t[1]}
    specs = _combined_sorted_row_specs(tables, subdomain_order=subdomain_order)
    col_order = _combined_structure_column_order(specs, all_models)
    out_dir = os.path.dirname(out_path)
    if out_dir and not os.path.isdir(out_dir):
        os.makedirs(out_dir, exist_ok=True)
    with open(out_path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["Model"] + col_order + ["Subdomain"])
        for sub, name, ids, matrix, i in specs:
            id_to_idx = {lab: j for j, lab in enumerate(ids)}
            row = [name]
            for col in col_order:
                if col == name:
                    row.append("-")
                elif col in id_to_idx:
                    j = id_to_idx[col]
                    v = matrix[i][j]
                    if v is None or (isinstance(v, float) and (v != v or v < 0)):
                        row.append("")
                    else:
                        row.append(f"{v:.4f}")
                else:
                    row.append("")
            row.append(sub)
            w.writerow(row)


def disp_limits_from_rmsd_values_array(
    arr,
    vmin: float | None,
    vmax: float | None,
) -> tuple[float, float] | None:
    """
    vmin/vmax for heatmaps from a float array of RMSD values (NaN allowed).
    Prefers positive finite values for auto limits; honors vmin/vmax when set.
    """
    try:
        import numpy as np
    except ImportError:
        return None
    finite = arr[np.isfinite(arr)]
    pos = finite[finite > 0]
    if len(pos) > 0:
        data_min = float(np.min(pos))
        data_max = float(np.max(pos))
    elif len(finite) > 0:
        data_min = float(np.min(finite))
        data_max = float(np.max(finite))
    else:
        return None

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
    return disp_vmin, disp_vmax


def disp_vmin_vmax_for_combined_tables(
    tables: list[tuple[str, list[str], list[list[float]]]],
    vmin: float | None,
    vmax: float | None,
    subdomain_order: list[str] | None = None,
) -> tuple[float, float] | None:
    """Same limits as the merged table / combined CSV (in-memory construction)."""
    try:
        import numpy as np
    except ImportError:
        return None
    if not tables:
        return None
    all_models = {x for t in tables for x in t[1]}
    specs = _combined_sorted_row_specs(tables, subdomain_order=subdomain_order)
    col_order = _combined_structure_column_order(specs, all_models)
    n_r, n_c = len(specs), len(col_order)
    arr = np.full((n_r, n_c), np.nan, dtype=np.float64)
    for i, (_sub, name, ids, matrix, ri) in enumerate(specs):
        id_to_idx = {lab: j for j, lab in enumerate(ids)}
        for j, cid in enumerate(col_order):
            if cid == name:
                continue
            if cid not in id_to_idx:
                continue
            v = matrix[ri][id_to_idx[cid]]
            if v is None or (isinstance(v, float) and (v != v or v < 0)):
                continue
            arr[i, j] = float(v)

    return disp_limits_from_rmsd_values_array(arr, vmin, vmax)


def load_combined_rmsd_csv_matrix(csv_path: str):
    """
    Read wide combined_rmsd_table.csv into a numpy array and metadata.
    Returns (arr, col_ids, model_names, rows) or None if missing/invalid/empty.
    """
    try:
        import numpy as np
    except ImportError:
        return None
    csv_path = os.path.abspath(csv_path)
    if not os.path.isfile(csv_path):
        return None
    with open(csv_path, newline="") as f:
        reader = csv.DictReader(f)
        if not reader.fieldnames or "Model" not in reader.fieldnames or "Subdomain" not in reader.fieldnames:
            return None
        col_ids = [c for c in reader.fieldnames if c not in ("Model", "Subdomain")]
        rows = list(reader)
    if not rows or not col_ids:
        return None
    n_r, n_c = len(rows), len(col_ids)
    arr = np.full((n_r, n_c), np.nan, dtype=np.float64)
    model_names = [r.get("Model", "").strip() for r in rows]
    for i, r in enumerate(rows):
        m = model_names[i]
        for j, cid in enumerate(col_ids):
            if cid == m:
                continue
            v = _parse_rmsd_table_cell(r.get(cid, ""))
            if v is not None:
                arr[i, j] = v
    return arr, col_ids, model_names, rows


def vmin_vmax_from_combined_csv_path(
    csv_path: str,
    vmin: float | None,
    vmax: float | None,
) -> tuple[float, float] | None:
    """Colour limits from the merged CSV on disk (same rules as the combined heatmap)."""
    data = load_combined_rmsd_csv_matrix(csv_path)
    if data is None:
        return None
    arr, _col_ids, _model_names, _rows = data
    return disp_limits_from_rmsd_values_array(arr, vmin, vmax)


def median_rmsd_from_combined_csv_path(csv_path: str) -> float | None:
    """Median of all finite RMSD entries in the merged CSV (same cells as the combined heatmap)."""
    try:
        import numpy as np
    except ImportError:
        return None
    data = load_combined_rmsd_csv_matrix(csv_path)
    if data is None:
        return None
    arr, _, _, _ = data
    finite = arr[np.isfinite(arr)]
    if len(finite) == 0:
        return None
    return float(np.median(finite))


def median_rmsd_from_combined_array(arr) -> float | None:
    """Median of all finite values in the combined RMSD array."""
    try:
        import numpy as np
    except ImportError:
        return None
    finite = arr[np.isfinite(arr)]
    if len(finite) == 0:
        return None
    return float(np.median(finite))


def _apply_label_order(
    args: argparse.Namespace,
    ids: list[str],
    matrix: list[list[float]],
    *,
    context: str,
    hatch_matrix: list[list[float | None]] | None = None,
) -> tuple[list[str], list[list[float]], list[list[float | None]] | None]:
    """Apply --order / --order-dir to ids and matrix; optionally the same permutation for hatch."""
    if args.order_dir:
        glob_pat = args.order_glob or "*.pdb"
        root_abs = os.path.abspath(args.order_dir)
        if not os.path.isdir(root_abs):
            print("Error: --order-dir is not a directory:", root_abs, file=sys.stderr)
            sys.exit(1)
        order_list = order_labels_from_scan(root_abs, glob_pat)
        if order_list:
            old_ids = ids
            ids, matrix = reorder_matrix(old_ids, matrix, order_list)
            if hatch_matrix is not None:
                ids, hatch_matrix = reorder_hatch_matrix(
                    old_ids, hatch_matrix, order_list
                )
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
            old_ids = ids
            ids, matrix = reorder_matrix(old_ids, matrix, order_list)
            if hatch_matrix is not None:
                ids, hatch_matrix = reorder_hatch_matrix(
                    old_ids, hatch_matrix, order_list
                )
            print(f"[{context}] Applied --order: {len(order_list)} labels", file=sys.stderr)
        else:
            print(
                f"Warning [{context}]: --order empty or file not found; using default order.",
                file=sys.stderr,
            )
    elif args.order_glob:
        print("Error: --order-glob requires --order-dir.", file=sys.stderr)
        sys.exit(1)
    return ids, matrix, hatch_matrix


def _parse_rmsd_table_cell(s: str) -> float | None:
    t = (s or "").strip()
    if not t or t == "-":
        return None
    try:
        return float(t)
    except ValueError:
        return None


def plot_combined_rmsd_heatmap_from_csv(
    csv_path: str,
    out_path: str,
    title: str = "Combined RMSD (all subdomains)",
    cmap: str = "viridis_r",
    vmin: float | None = None,
    vmax: float | None = None,
    short_labels: bool = False,
    diverging_center: str | None = None,
    vcenter: float | None = None,
    colorbar_orientation: str = "vertical",
    y_axis_right: bool = False,
) -> None:
    """Heatmap from wide combined_rmsd_table.csv (Model + structure columns + Subdomain)."""
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        import numpy as np
    except ImportError as e:
        print("Matplotlib required for combined heatmap. Install with: pip install matplotlib", file=sys.stderr)
        raise SystemExit(1) from e

    data = load_combined_rmsd_csv_matrix(csv_path)
    if data is None:
        print("Error: Combined CSV not found or invalid:", csv_path, file=sys.stderr)
        raise SystemExit(1)
    arr, col_ids, model_names, rows = data

    limits = disp_limits_from_rmsd_values_array(arr, vmin, vmax)
    if limits is None:
        print("No finite RMSD values in combined CSV to plot.", csv_path, file=sys.stderr)
        return
    disp_vmin, disp_vmax = limits

    dc = diverging_center if diverging_center else "none"
    vc_for_norm: float | None = None
    if dc == "median":
        vc_for_norm = vcenter if vcenter is not None else median_rmsd_from_combined_array(arr)
        if vc_for_norm is None:
            print("Warning: could not compute median for combined heatmap; using linear colour scale.", file=sys.stderr)
            dc = "none"
        elif vcenter is None:
            print("Combined heatmap: median centre = %.4g Å" % vc_for_norm, file=sys.stderr)
    norm = make_heatmap_norm(disp_vmin, disp_vmax, dc, vc_for_norm if dc == "median" else None)

    n_r, n_c = arr.shape

    arr_disp = arr.copy()
    for i, m in enumerate(model_names):
        for j, cid in enumerate(col_ids):
            if cid == m or not np.isfinite(arr[i, j]):
                arr_disp[i, j] = disp_vmin - 1.0

    if short_labels:
        row_labels = [short_heatmap_label(r.get("Model", "").strip()) for r in rows]
        col_tick_labels = [short_heatmap_label(c) for c in col_ids]
    else:
        row_labels = [f"{r.get('Subdomain', '')} | {r.get('Model', '')}" for r in rows]
        col_tick_labels = list(col_ids)

    try:
        cmap_obj = plt.get_cmap(cmap).copy()
    except ValueError:
        print("Unknown colour map '%s'. Common: viridis_r, plasma_r, RdYlBu_r, coolwarm, YlOrRd." % cmap, file=sys.stderr)
        cmap_obj = plt.get_cmap("viridis_r").copy()
    cmap_obj.set_under(color="white")

    out_path = os.path.abspath(out_path)
    plot_dir = os.path.dirname(out_path)
    if plot_dir and not os.path.isdir(plot_dir):
        os.makedirs(plot_dir, exist_ok=True)

    y_fs = max(4, min(9, 140 // max(n_r, 1)))
    x_fs = max(5, min(10, 200 // max(n_c, 1)))

    ext = os.path.splitext(out_path)[1].lower()
    fmt = "png" if ext == ".png" else "pdf" if ext == ".pdf" else "svg" if ext == ".svg" else "png"

    cb_orient = colorbar_orientation.lower() if colorbar_orientation else "vertical"
    if cb_orient not in ("vertical", "horizontal"):
        cb_orient = "vertical"
    fig_h = max(6, n_r * 0.22)
    if cb_orient == "horizontal":
        fig_h += 1.45
    fig, ax = plt.subplots(figsize=(max(8, n_c * 0.35), fig_h))
    im = None
    if _heatmap_vector_format(fmt):
        _draw_heatmap_vector_cells(ax, arr_disp, cmap_obj, norm)
    else:
        im = ax.imshow(
            arr_disp,
            aspect="equal",
            origin="lower",
            cmap=cmap_obj,
            interpolation="nearest",
            norm=norm,
        )
    ax.set_xticks(range(n_c))
    ax.set_yticks(range(n_r))
    ax.set_xticklabels(col_tick_labels, rotation=45, ha="right", fontsize=x_fs)
    ax.set_yticklabels(row_labels, fontsize=y_fs)
    ax.set_title(title)
    if y_axis_right:
        _apply_heatmap_y_axis_right(ax)
    cbar_left = y_axis_right and cb_orient == "vertical"
    if cb_orient == "horizontal":
        r_edge = 0.86 if y_axis_right else 0.94
        plt.tight_layout(rect=(0.06, 0.22, r_edge, 0.94))
    elif cbar_left:
        plt.tight_layout(rect=(0.16, 0.06, 0.90, 0.92))
    else:
        plt.tight_layout(rect=(0.08, 0.06, 0.88, 0.92))
    if _heatmap_vector_format(fmt):
        _add_raster_colorbar(
            ax,
            cmap_obj,
            norm,
            nrows=n_r,
            ncols=n_c,
            label="RMSD (Å)",
            orientation=cb_orient,
            mappable=None,
            vertical_colorbar_on_left=cbar_left,
        )
    else:
        _add_raster_colorbar(
            ax,
            cmap_obj,
            norm,
            nrows=n_r,
            ncols=n_c,
            label="RMSD (Å)",
            orientation=cb_orient,
            mappable=im,
            vertical_colorbar_on_left=cbar_left,
        )

    fd, tmp_path = tempfile.mkstemp(suffix=ext, prefix="rmsd_combined_heatmap_", dir=plot_dir or ".")
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
            print("Error saving combined heatmap:", e, file=sys.stderr)
            raise SystemExit(1) from e
        plt.close(fig)

        size = os.path.getsize(tmp_path) if os.path.exists(tmp_path) else 0
        min_ok = 3000 if fmt in ("svg", "pdf") else 4000
        if size < min_ok or not _image_magic_ok(tmp_path, fmt):
            try:
                os.remove(tmp_path)
            except OSError:
                pass
            print("Combined heatmap save produced an invalid file (%d bytes)." % size, file=sys.stderr)
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


def run_single_file(args: argparse.Namespace) -> None:
    path = os.path.abspath(args.input)
    if not os.path.isfile(path):
        print("Error: Input file not found:", path, file=sys.stderr)
        sys.exit(1)

    if args.plot_hatch:
        pack = rmsd_file_to_matrix_with_hatch(path)
        if not pack[0] or not pack[1]:
            print("Error: No RMSD pairs found in", path, file=sys.stderr)
            sys.exit(1)
        ids, matrix, _fmt, hatch_mat, hatch_legend = pack
    else:
        ids, matrix = rmsd_file_to_matrix(path)
        hatch_mat, hatch_legend = None, None
    if not ids or not matrix:
        print("Error: No RMSD pairs found in", path, file=sys.stderr)
        sys.exit(1)

    print("Loaded", len(ids), "models, matrix", len(ids), "x", len(ids))

    ids, matrix, hatch_mat = _apply_label_order(
        args, ids, matrix, context=os.path.basename(path), hatch_matrix=hatch_mat
    )

    if args.output:
        out_csv = args.output
        ext = os.path.splitext(out_csv)[1].lower()
        if ext in (".png", ".svg", ".pdf"):
            print("Note: -o/--output is for the CSV table. You passed an image extension (%s)." % ext, file=sys.stderr)
            print("To save the heatmap image, use: --plot %s" % out_csv, file=sys.stderr)
            print("Writing CSV to same path with .csv extension.", file=sys.stderr)
            out_csv = os.path.splitext(out_csv)[0] + ".csv"
    else:
        base = os.path.basename(path)
        if base.startswith("rmsd_SSM_values") and base.endswith(".txt"):
            suffix = base[15:-4] if len(base) > 19 else "ssm"
        elif base.startswith("rmsd_values_") and base.endswith(".txt"):
            suffix = base[12:-4] if len(base) > 16 else "lsq"
        else:
            suffix = "rmsd"
        out_csv = os.path.join(os.path.dirname(path), f"rmsd_table_{suffix}.csv")

    write_rmsd_csv(ids, matrix, out_csv)
    print("Wrote", out_csv)

    if args.plot:
        ph: dict = {}
        if args.plot_hatch and hatch_mat is not None and hatch_legend:
            ph = {
                "hatch_value_matrix": hatch_mat,
                "hatch_channel_name": hatch_legend,
                "hatch_n_equal_bins": max(1, int(args.heatmap_hatch_bins)),
            }
            if args.heatmap_hatch_edges:
                parts = [p.strip() for p in str(args.heatmap_hatch_edges).split(",") if p.strip()]
                if len(parts) < 2:
                    print(
                        "Error: --heatmap-hatch-edges needs at least two comma-separated values.",
                        file=sys.stderr,
                    )
                    sys.exit(1)
                ph["hatch_bin_edges"] = [float(p) for p in parts]
        elif args.plot_hatch:
            print(
                "Warning: --plot-hatch ignored (no pair counts in this input, or not LSQ/SSM .txt).",
                file=sys.stderr,
            )
        plot_heatmap(
            ids,
            matrix,
            args.plot,
            title=args.title,
            cmap=args.cmap,
            vmin=args.vmin,
            vmax=args.vmax,
            short_labels=args.short_heatmap_labels,
            diverging_center=args.heatmap_diverging_center,
            vcenter=None,
            colorbar_orientation=args.heatmap_colorbar_orientation,
            y_axis_right=args.heatmap_y_axis_right,
            **ph,
        )


def run_scan_dir(args: argparse.Namespace) -> None:
    root = _resolved_path(args.scan_dir)
    if not os.path.isdir(root):
        print("Error: --scan-dir is not a directory:", root, file=sys.stderr)
        sys.exit(1)

    sources = discover_rmsd_files(root, args.scan_glob)
    if not sources:
        print(
            "Error: No files matched",
            repr(args.scan_glob),
            "under",
            root,
            file=sys.stderr,
        )
        sys.exit(1)

    print(f"Scan {root!r}: {len(sources)} file(s) matching {args.scan_glob!r}", file=sys.stderr)

    if args.heatmap_dir:
        hdir = _resolved_path(args.heatmap_dir)
        os.makedirs(hdir, exist_ok=True)

    tables: list[
        tuple[str, list[str], list[list[float]], list[list[float | None]] | None, str | None]
    ] = []
    for path in sources:
        label = subdomain_label_from_rmsd_path(path)
        if args.plot_hatch and args.heatmap_dir:
            pack = rmsd_file_to_matrix_with_hatch(path)
            if not pack[0] or not pack[1]:
                print(f"Warning: skip (no pairs): {path}", file=sys.stderr)
                continue
            ids, matrix, _fmt, hatch_mat, hatch_leg = pack
        else:
            ids, matrix = rmsd_file_to_matrix(path)
            hatch_mat, hatch_leg = None, None
        if not ids or not matrix:
            print(f"Warning: skip (no pairs): {path}", file=sys.stderr)
            continue
        ids, matrix, hatch_mat = _apply_label_order(
            args, ids, matrix, context=label, hatch_matrix=hatch_mat
        )
        out_csv = os.path.join(os.path.dirname(path), f"rmsd_table_{label}.csv")
        write_rmsd_csv(ids, matrix, out_csv)
        print("Wrote", out_csv)
        tables.append((label, ids, matrix, hatch_mat, hatch_leg))

    if not tables:
        print("Error: No valid RMSD tables produced.", file=sys.stderr)
        sys.exit(1)

    combined_path: str | None = None
    combined_sub_order: list[str] | None = None
    if args.combined_subdomain_order:
        combined_sub_order = parse_order(args.combined_subdomain_order)
        if not combined_sub_order:
            print("Warning: --combined-subdomain-order produced no tokens; using default subdomain sort.", file=sys.stderr)

    if not args.no_combined:
        combined_path = _resolved_path(args.combined_csv or os.path.join(root, "combined_rmsd_table.csv"))
        write_combined_rmsd_csv(
            tables,
            combined_path,
            subdomain_order=combined_sub_order if combined_sub_order else None,
        )
        print("Wrote", combined_path)

    use_shared_scale = bool(args.heatmap_dir) and not args.no_shared_heatmap_scale
    shared_vmin: float | None = None
    shared_vmax: float | None = None
    if use_shared_scale:
        sv = None
        if combined_path is not None:
            sv = vmin_vmax_from_combined_csv_path(combined_path, args.vmin, args.vmax)
            if sv is not None:
                shared_vmin, shared_vmax = sv
                print(
                    "Using merged CSV colour scale for per-subdomain heatmaps: vmin=%s vmax=%s"
                    % (shared_vmin, shared_vmax),
                    file=sys.stderr,
                )
        if sv is None:
            sv = disp_vmin_vmax_for_combined_tables(
                tables,
                args.vmin,
                args.vmax,
                subdomain_order=combined_sub_order if combined_sub_order else None,
            )
            if sv is not None:
                shared_vmin, shared_vmax = sv
                note = (
                    " (--no-combined)"
                    if combined_path is None
                    else " (fallback: same limits as merged table)"
                )
                print(
                    "Using merged-table scale from memory%s: vmin=%s vmax=%s"
                    % (note, shared_vmin, shared_vmax),
                    file=sys.stderr,
                )
        if sv is None:
            print(
                "Warning: could not derive merged-table colour scale; per-subdomain heatmaps use local auto scale.",
                file=sys.stderr,
            )

    global_median: float | None = None
    if (
        args.heatmap_diverging_center == "median"
        and combined_path is not None
        and os.path.isfile(combined_path)
    ):
        global_median = median_rmsd_from_combined_csv_path(combined_path)
        if global_median is not None:
            print(
                "Diverging centre (median from merged CSV): %.4g Å" % global_median,
                file=sys.stderr,
            )

    if args.heatmap_dir:
        hf = args.heatmap_format
        hdir_abs = _resolved_path(args.heatmap_dir)
        for label, ids, matrix, hatch_mat, hatch_leg in tables:
            hp = os.path.join(hdir_abs, f"rmsd_heatmap_{label}.{hf}")
            if use_shared_scale and shared_vmin is not None and shared_vmax is not None:
                pv, px = shared_vmin, shared_vmax
            else:
                pv, px = args.vmin, args.vmax
            ph: dict = {}
            if args.plot_hatch and hatch_mat is not None and hatch_leg:
                ph = {
                    "hatch_value_matrix": hatch_mat,
                    "hatch_channel_name": hatch_leg,
                    "hatch_n_equal_bins": max(1, int(args.heatmap_hatch_bins)),
                }
                if args.heatmap_hatch_edges:
                    parts = [p.strip() for p in str(args.heatmap_hatch_edges).split(",") if p.strip()]
                    if len(parts) < 2:
                        print(
                            "Error: --heatmap-hatch-edges needs at least two comma-separated values.",
                            file=sys.stderr,
                        )
                        sys.exit(1)
                    ph["hatch_bin_edges"] = [float(p) for p in parts]
            elif args.plot_hatch:
                print(
                    f"Warning: --plot-hatch: no pair counts for subdomain {label!r} (use LSQ/SSM .txt).",
                    file=sys.stderr,
                )
            plot_heatmap(
                ids,
                matrix,
                hp,
                title=f"{args.title} ({label})" if args.title else f"RMSD heatmap ({label})",
                cmap=args.cmap,
                vmin=pv,
                vmax=px,
                short_labels=args.short_heatmap_labels,
                diverging_center=args.heatmap_diverging_center,
                vcenter=global_median,
                colorbar_orientation=args.heatmap_colorbar_orientation,
                y_axis_right=args.heatmap_y_axis_right,
                **ph,
            )

    if not args.no_combined and combined_path is not None:
        want_combined_hm = (args.heatmap_dir or args.combined_heatmap) and not args.no_combined_heatmap
        if want_combined_hm:
            hf = args.heatmap_format
            if args.combined_heatmap:
                ch_out = ensure_batch_heatmap_path(
                    _resolved_path(args.combined_heatmap),
                    hf,
                )
            else:
                ch_out = os.path.join(
                    _resolved_path(args.heatmap_dir),
                    f"combined_rmsd_heatmap.{hf}",
                )
            comb_title = f"{args.title} (combined)" if args.title else "Combined RMSD (all subdomains)"
            plot_combined_rmsd_heatmap_from_csv(
                combined_path,
                ch_out,
                title=comb_title,
                cmap=args.cmap,
                vmin=args.vmin,
                vmax=args.vmax,
                short_labels=args.short_heatmap_labels,
                diverging_center=args.heatmap_diverging_center,
                vcenter=global_median,
                colorbar_orientation=args.heatmap_colorbar_orientation,
                y_axis_right=args.heatmap_y_axis_right,
            )


def main() -> None:
    class _HelpFmt(
        argparse.RawDescriptionHelpFormatter,
        argparse.ArgumentDefaultsHelpFormatter,
    ):
        pass

    epilog = """
SINGLE-FILE MODE (provide INPUT; do not use --scan-dir):
  Read one pairwise RMSD log (Coot LSQ or SSM style, as written by extract_rmsd / superposition
  workflows) or an existing square CSV with a Model column and one column per structure.
  Formats match structure_phylogeny.py (LSQ: Alignment + core rmsd; SSM: Superposing/Aligning
  + core rmsd; CSV: re-read and rewrite).
  All structure names from the file are placed in a symmetric matrix (natural-sorted row/column
  order unless --order or --order-dir is set). Write one square CSV: column Model, then one column
  per structure, diagonal '-', off-diagonal RMSD (Å). Default output path is next to INPUT:
  rmsd_table_<suffix>.csv where <suffix> comes from the input basename (e.g. rmsd_values_run1.txt
  -> run1). Use -o to choose another path. Optional --plot saves one heatmap image for that matrix.

BATCH MODE (--scan-dir DIR):
  Recursively find every file under DIR matching --scan-glob (default: rmsd_values_*.txt).
  For each match: parse like single-file mode, apply --order / --order-dir / --order-glob if set
  (same rules as single-file, repeated for every matched log), then write rmsd_table_<suffix>.csv
  beside that file (<suffix> from the log filename). Then stack every subdomain into one wide
  combined CSV (default DIR/combined_rmsd_table.csv): rows sorted by Subdomain then Model (override block
  order with --combined-subdomain-order, e.g. D1,D2,ID,CWB2); structure
  columns follow that row order (first-seen model per sorted rows). Columns Model + those names +
  Subdomain. Empty cells where a subdomain has no pair. Use --no-combined to skip the merged file. For figures in batch mode use --heatmap-dir
  (one SVG per log plus combined_rmsd_heatmap.svg from the merged CSV by default), or set --combined-heatmap
  PATH to write only the combined figure (file extension sets format if present; else --heatmap-format).
  Use --heatmap-format png|pdf|svg for auto-named batch outputs (default svg). Use --no-combined-heatmap
  to skip the combined figure when using --heatmap-dir. Single-file mode uses --plot (extension sets format).
  Per-subdomain heatmaps use the merged CSV colour scale by default; use --no-shared-heatmap-scale for local auto scale.
  --vmin/--vmax apply to the merged scale and the combined heatmap when set.

Examples (from repository root):
  python ranking/rmsd_to_csv.py /path/to/project/rmsd_SSM_values.txt -o /path/to/project/rmsd_table.csv
  python ranking/rmsd_to_csv.py --scan-dir /path/to/project/run1/ --scan-glob 'rmsd_SSM_values*.txt'
"""
    ap = argparse.ArgumentParser(
        description="Pairwise RMSD logs or RMSD CSV -> square CSV table(s); optional matplotlib heatmaps.",
        epilog=epilog,
        formatter_class=_HelpFmt,
    )
    ap.add_argument(
        "input",
        nargs="?",
        default=None,
        help="Single-file mode only: path to rmsd_values_*.txt, rmsd_SSM_values*.txt, or rmsd_table_*.csv. "
        "Omit when using --scan-dir.",
    )
    ap.add_argument(
        "--scan-dir",
        metavar="DIR",
        default=None,
        help="Batch mode: recursively scan DIR for --scan-glob matches; write per-file CSVs plus "
        "combined_rmsd_table.csv (unless --no-combined / --combined-csv). Incompatible with INPUT.",
    )
    ap.add_argument(
        "--scan-glob",
        metavar="PATTERN",
        default="rmsd_values_*.txt",
        help="Batch only: glob relative to each subdirectory under --scan-dir (default: rmsd_values_*.txt). "
        "Example: rmsd_SSM_values*.txt",
    )
    ap.add_argument(
        "--combined-csv",
        metavar="PATH",
        default=None,
        help="Batch only: output path for the wide merged table (default: <scan-dir>/combined_rmsd_table.csv).",
    )
    ap.add_argument(
        "--combined-subdomain-order",
        metavar="LIST_OR_FILE",
        default=None,
        help="Batch only: comma-separated tokens or path (one token per line) giving subdomain block order in "
        "combined_rmsd_table.csv and the combined heatmap. Each token matches one subdomain label with an "
        "underscore boundary (e.g. D1 matches …_D1s but not …_D10s). Unmatched subdomains follow in natural sort.",
    )
    ap.add_argument(
        "--no-combined",
        action="store_true",
        help="Batch only: write per-file rmsd_table_*.csv but not combined_rmsd_table.csv.",
    )
    ap.add_argument(
        "--heatmap-dir",
        metavar="DIR",
        default=None,
        help="Batch only: directory for rmsd_heatmap_<suffix>.<fmt> per matched log, and by default "
        "combined_rmsd_heatmap.<fmt> from the merged CSV (unless --no-combined or --no-combined-heatmap). "
        "<fmt> is set by --heatmap-format (default svg).",
    )
    ap.add_argument(
        "--heatmap-format",
        choices=("png", "pdf", "svg"),
        default="svg",
        help="Batch only: image format extension for auto-named heatmaps under --heatmap-dir and for the "
        "default combined path. Ignored when --combined-heatmap ends with .png/.pdf/.svg. "
        "SVG/PDF: vector matrix cells; colour bar is one embedded image (PNG: raster matrix).",
    )
    ap.add_argument(
        "--combined-heatmap",
        metavar="PATH",
        default=None,
        help="Batch only: write merged-table heatmap to PATH. Extension .png/.pdf/.svg selects format; "
        "if omitted, append .<heatmap-format>. Default with --heatmap-dir only: "
        "<heatmap-dir>/combined_rmsd_heatmap.<heatmap-format>. Requires merged CSV (not --no-combined).",
    )
    ap.add_argument(
        "--no-combined-heatmap",
        action="store_true",
        help="Batch only: with --heatmap-dir, do not write combined_rmsd_heatmap.<fmt>.",
    )
    ap.add_argument(
        "--no-shared-heatmap-scale",
        action="store_true",
        help="Batch only: with --heatmap-dir, use auto vmin/vmax per subdomain matrix instead of the "
        "merged CSV / combined heatmap scale.",
    )
    ap.add_argument(
        "-o", "--output",
        default=None,
        help="Single-file only: output CSV path. Default: directory of INPUT, name rmsd_table_<suffix>.csv.",
    )
    ap.add_argument(
        "--plot",
        metavar="PATH",
        default=None,
        help="Single-file only: save one RMSD heatmap image to PATH (PNG/SVG/PDF; PNG recommended headless). "
        "SVG/PDF: vector matrix cells; colour bar is one embedded image.",
    )
    ap.add_argument(
        "--title",
        default="RMSD heatmap",
        help="Heatmap title; batch mode appends (suffix) for each figure from --heatmap-dir.",
    )
    order_group = ap.add_mutually_exclusive_group()
    order_group.add_argument(
        "--order",
        metavar="FILE_OR_LIST",
        default=None,
        help="Reorder matrix rows/columns: comma-separated labels, or path to a file with one label per line. "
        "Labels must match names in the RMSD file. Single-file: applies once; batch: applied to each matched log.",
    )
    order_group.add_argument(
        "--order-dir",
        metavar="DIR",
        default=None,
        help="Reorder matrix from a recursive file listing under DIR: each label is basename without extension, "
        "order follows relative paths (natural numeric sort). Use --order-glob (default *.pdb). Mutually exclusive "
        "with --order.",
    )
    ap.add_argument(
        "--order-glob",
        metavar="PATTERN",
        default=None,
        help="Glob for --order-dir only (default if --order-dir set: *.pdb). Example: model_*.pdb.",
    )
    ap.add_argument(
        "--cmap",
        default="viridis_r",
        help="Matplotlib colour map (`cmap`) for any heatmap (single-file --plot or batch --heatmap-dir).",
    )
    ap.add_argument(
        "--vmin",
        type=float,
        default=None,
        metavar="VAL",
        help="Heatmap colour scale floor (Å); default: data min in the matrix or merged CSV (unless "
        "--no-shared-heatmap-scale for per-subdomain plots).",
    )
    ap.add_argument(
        "--vmax",
        type=float,
        default=None,
        metavar="VAL",
        help="Heatmap colour scale ceiling (Å); default: data max in the matrix or merged CSV (unless "
        "--no-shared-heatmap-scale for per-subdomain plots).",
    )
    ap.add_argument(
        "--short-heatmap-labels",
        action="store_true",
        help="Heatmap figures only: shorten axis labels, e.g. S1refAF3_CWB2s.pdb -> S1_CWB2s (leading "
        "letter+digits from first segment + last underscore segment). CSV tables unchanged.",
    )
    ap.add_argument(
        "--heatmap-diverging-center",
        choices=("none", "median"),
        default="none",
        help="Use matplotlib TwoSlopeNorm with centre at the median RMSD so diverging colormaps "
        "(e.g. PuOr_r, RdBu_r) map low vs high around the median. Batch: median from merged CSV for "
        "all figures when combined_rmsd_table.csv is written; with --no-combined, per-matrix median. "
        "Single-file: median of that matrix. Default: none (linear scale).",
    )
    ap.add_argument(
        "--heatmap-colorbar-orientation",
        choices=("vertical", "horizontal"),
        default="vertical",
        help="Colour bar layout for heatmap figures (single-file --plot or batch --heatmap-dir). "
        "horizontal places the bar below the matrix.",
    )
    ap.add_argument(
        "--heatmap-y-axis-right",
        action="store_true",
        help="Heatmap figures: draw row (y) tick labels on the right. With a vertical colour bar, "
        "the bar is placed on the left so it does not overlap the labels.",
    )
    ap.add_argument(
        "--plot-hatch",
        action="store_true",
        help="With --plot or --heatmap-dir: add a second channel (hatch pattern) binned from Coot "
        "per-pair counts in LSQ/SSM .txt logs — LSQ: matched atoms; SSM: aligned residues. "
        "Legend titles state the source. Ignored for CSV-only input (no counts in file).",
    )
    ap.add_argument(
        "--heatmap-hatch-bins",
        type=int,
        default=4,
        metavar="N",
        help="With --plot-hatch: number of equal-width bins for the count axis (default: 4).",
    )
    ap.add_argument(
        "--heatmap-hatch-edges",
        metavar="E0,E1,…",
        default=None,
        help="With --plot-hatch: optional comma-separated bin boundaries; extended to data min/max if needed.",
    )
    add_log_args(ap)
    args = ap.parse_args()
    summary_log = setup_log_from_args(
        args,
        script_path=__file__,
        inputs=[x for x in [getattr(args, "input", None), getattr(args, "scan_dir", None)] if x],
        pattern=getattr(args, "scan_glob", None),
    )
    if summary_log is not None:
        if args.scan_dir:
            summary_log.task("RMSD to CSV (batch: --scan-dir)")
            summary_log.kv("scan_dir", args.scan_dir)
            summary_log.kv("scan_glob", args.scan_glob)
            summary_log.kv("combined_csv", args.combined_csv or "(default)")
            summary_log.kv("no_combined", bool(args.no_combined))
            if args.heatmap_dir:
                summary_log.kv("heatmap_dir", args.heatmap_dir)
        else:
            summary_log.task("RMSD to CSV (single file)")
            summary_log.kv("input", args.input)
            summary_log.kv("output", args.output or "(default)")
            if args.plot:
                summary_log.kv("plot", args.plot)

    if args.scan_dir and args.input:
        stray = args.input.strip()
        if stray.lower() in ("png", "pdf", "svg") and not os.path.isfile(stray):
            ap.error(
                "Do not pass a bare format name as INPUT with --scan-dir (got %r). "
                "Put a space before --heatmap-format, e.g. .../figs --heatmap-format %s"
                % (stray, stray.lower())
            )
        else:
            ap.error(
                "Do not pass INPUT when using --scan-dir (remove the extra argument %r)" % stray
            )
    if not args.scan_dir and not args.input:
        ap.error("Provide INPUT, or use --scan-dir DIR for batch mode")
    if args.scan_dir and args.plot:
        ap.error("In --scan-dir mode use --heatmap-dir for figures, not --plot (which is for a single output file)")
    if args.heatmap_dir and not args.scan_dir:
        ap.error("--heatmap-dir requires --scan-dir")
    if args.combined_csv and not args.scan_dir:
        ap.error("--combined-csv is only used with --scan-dir")
    if args.no_combined and not args.scan_dir:
        ap.error("--no-combined is only used with --scan-dir")
    if args.combined_heatmap and not args.scan_dir:
        ap.error("--combined-heatmap is only used with --scan-dir")
    if args.no_combined_heatmap and not args.scan_dir:
        ap.error("--no-combined-heatmap is only used with --scan-dir")
    if args.no_combined_heatmap and args.combined_heatmap:
        ap.error("Do not combine --combined-heatmap with --no-combined-heatmap")
    if args.combined_heatmap and args.no_combined:
        ap.error("--combined-heatmap requires a merged CSV; omit --no-combined")
    if args.no_shared_heatmap_scale and not args.scan_dir:
        ap.error("--no-shared-heatmap-scale is only used with --scan-dir")

    if args.scan_dir:
        run_scan_dir(args)
    else:
        run_single_file(args)


if __name__ == "__main__":
    main()
