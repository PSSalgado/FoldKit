#!/usr/bin/env python3
"""
Convert pairwise RMSD logs (or an existing RMSD CSV) to square CSV table(s), optional heatmaps.

Single-file mode (pass one INPUT path; omit --scan-dir):
  Reads one LSQ text, SSM text, or square rmsd_table_*.csv (same format detection as
  structure_phylogeny). Builds a symmetric pairwise RMSD matrix with natural-sorted
  labels, writes one CSV (default: beside INPUT as rmsd_table_<suffix>.csv), and
  optionally --plot for one matplotlib heatmap. Use --order / --order-dir to reorder
  rows and columns to match PDB stems or a list.

Batch mode (--scan-dir):
  Recursively finds logs matching --scan-glob (default rmsd_values_*.txt); for each,
  writes rmsd_table_<suffix>.csv next to the log and stacks all blocks into
  combined_rmsd_table.csv with a Subdomain column. Optional --heatmap-dir writes figures
  per file plus combined_rmsd_heatmap.<fmt> (--heatmap-format png|pdf|svg, default png;
  SVG and PDF draw each matrix cell as a vector square; the colour bar may still be a bitmap strip);
  --combined-heatmap PATH for a custom combined path (extension overrides --heatmap-format);
  --no-combined-heatmap to skip the combined figure.

Examples:
  python rmsd_to_csv.py path/to/rmsd_SSM_values.txt -o path/to/rmsd_table.csv
  python rmsd_to_csv.py path/to/rmsd_SSM_values.txt --plot heatmap.pdf
  python rmsd_to_csv.py path/to/rmsd_values.txt --order-dir /path/to/root --order-glob 'model_*.pdb'
  python rmsd_to_csv.py --scan-dir /path/to/base --scan-glob 'rmsd_values_*.txt'
  python rmsd_to_csv.py --scan-dir /path/to/base --heatmap-dir /path/to/figs

Outputs feed structure_phylogeny.py and create_rmsd_heatmap.R.
"""

from __future__ import annotations

import argparse
import csv
import glob
import os
import sys
import tempfile
from pathlib import Path

# Reuse parsing and matrix logic from structure_phylogeny
try:
    from structure_phylogeny import (
        _natural_sort_key,
        alignments_to_matrix,
        detect_format,
        parse_lsq_rmsd_txt,
        parse_rmsd_csv,
        parse_ssm_rmsd_txt,
    )
except ImportError:
    # Run from repo root or FoldKit/
    _dir = os.path.dirname(os.path.abspath(__file__))
    if _dir not in sys.path:
        sys.path.insert(0, _dir)
    from structure_phylogeny import (
        _natural_sort_key,
        alignments_to_matrix,
        detect_format,
        parse_lsq_rmsd_txt,
        parse_rmsd_csv,
        parse_ssm_rmsd_txt,
    )


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


def reorder_matrix(ids: list[str], matrix: list[list[float]], order: list[str]) -> tuple[list[str], list[list[float]]]:
    """Return (ordered_ids, reordered_matrix). order lists desired labels; any missing from matrix are appended."""
    id_to_idx = {x: i for i, x in enumerate(ids)}
    ordered_ids = [x for x in order if x in id_to_idx]
    missing = [x for x in ids if x not in ordered_ids]
    ordered_ids = ordered_ids + missing
    n = len(ordered_ids)
    new_matrix = [[0.0] * n for _ in range(n)]
    for i, a in enumerate(ordered_ids):
        for j, b in enumerate(ordered_ids):
            new_matrix[i][j] = matrix[id_to_idx[a]][id_to_idx[b]]
    return ordered_ids, new_matrix


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


def _combined_sorted_row_specs(
    tables: list[tuple[str, list[str], list[list[float]]]],
) -> list[tuple[str, str, list[str], list[list[float]], int]]:
    """(Subdomain, model, ids, matrix, row_index) for every row, sorted by subdomain then model."""
    specs: list[tuple[str, str, list[str], list[list[float]], int]] = []
    for sub, ids, matrix in tables:
        for i, name in enumerate(ids):
            specs.append((sub, name, ids, matrix, i))
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
) -> None:
    """Wide CSV: rows sorted by (Subdomain, Model); structure columns match row order (first-seen model order)."""
    if not tables:
        return
    all_models = {x for _, ids, _ in tables for x in ids}
    specs = _combined_sorted_row_specs(tables)
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


def disp_vmin_vmax_for_combined_tables(
    tables: list[tuple[str, list[str], list[list[float]]]],
    vmin: float | None,
    vmax: float | None,
) -> tuple[float, float] | None:
    """Color limits matching plot_combined_rmsd_heatmap_from_csv (merged table, same vmin/vmax rules)."""
    try:
        import numpy as np
    except ImportError:
        return None
    if not tables:
        return None
    all_models = {x for _, ids, _ in tables for x in ids}
    specs = _combined_sorted_row_specs(tables)
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


def _draw_heatmap_vector_cells(ax, arr, cmap_obj, vmin: float, vmax: float):
    """
    Draw one matplotlib Rectangle per matrix cell, matching imshow(origin='upper') layout.
    Returns a ScalarMappable suitable for ``plt.colorbar(..., extend='min')``.
    """
    import matplotlib.colors as mcolors
    import numpy as np
    from matplotlib.cm import ScalarMappable
    from matplotlib.patches import Rectangle

    nrow, ncol = int(arr.shape[0]), int(arr.shape[1])
    norm = mcolors.Normalize(vmin=vmin, vmax=vmax)
    under_rgba = mcolors.to_rgba(cmap_obj.get_under())
    for i in range(nrow):
        for j in range(ncol):
            v = arr[i, j]
            if np.isnan(v):
                fc = (1.0, 1.0, 1.0, 1.0)
            elif v < vmin:
                fc = under_rgba
            else:
                fc = cmap_obj(norm(v))
            y0 = (nrow - 1 - i) - 0.5
            x0 = j - 0.5
            ax.add_patch(
                Rectangle(
                    (x0, y0),
                    1.0,
                    1.0,
                    linewidth=0,
                    edgecolor="none",
                    facecolor=fc,
                )
            )
    ax.set_xlim(-0.5, ncol - 0.5)
    ax.set_ylim(nrow - 0.5, -0.5)
    ax.set_aspect("auto")
    sm = ScalarMappable(cmap=cmap_obj, norm=norm)
    sm.set_array(np.linspace(vmin, vmax, 256))
    return sm


def plot_heatmap(
    ids: list[str],
    matrix: list[list[float]],
    out_path: str,
    title: str = "RMSD heatmap",
    cmap: str = "viridis_r",
    vmin: float | None = None,
    vmax: float | None = None,
) -> None:
    """Draw heatmap with matplotlib and save to out_path as an image (PNG/SVG/PDF)."""
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        import numpy as np
    except ImportError as e:
        print("Matplotlib required for --plot. Install with: pip install matplotlib", file=sys.stderr)
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

    # Default colour scale: off-diagonal pairwise RMSDs only, and > 0 (exclude diagonal zeros).
    mask = ~np.eye(n, dtype=bool)
    off = arr[mask]
    off_finite = off[np.isfinite(off)]
    pos = off_finite[off_finite > 0]
    if len(pos) > 0:
        data_min = float(np.min(pos))
        data_max = float(np.max(pos))
    elif len(off_finite) > 0:
        data_min = float(np.min(off_finite))
        data_max = float(np.max(off_finite))
    else:
        print("No finite off-diagonal RMSD values to plot.", file=sys.stderr)
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

    fig, ax = plt.subplots(figsize=(max(6, n * 0.4), max(5, n * 0.35)))
    if _heatmap_vector_format(fmt):
        sm = _draw_heatmap_vector_cells(ax, arr, cmap_obj, disp_vmin, disp_vmax)
        plt.colorbar(sm, ax=ax, label="RMSD (Å)", extend="min")
    else:
        im = ax.imshow(
            arr,
            aspect="auto",
            cmap=cmap_obj,
            interpolation="nearest",
            vmin=disp_vmin,
            vmax=disp_vmax,
        )
        plt.colorbar(im, ax=ax, label="RMSD (Å)", extend="min")
    ax.set_xticks(range(n))
    ax.set_yticks(range(n))
    ax.set_xticklabels(ids, rotation=45, ha="right")
    ax.set_yticklabels(ids)
    ax.set_title(title)
    plt.tight_layout()

    # Save to a temp file in the same directory, then move to target (avoids partial/CSV overwrite)
    fd, tmp_path = tempfile.mkstemp(suffix=ext, prefix="rmsd_heatmap_", dir=plot_dir or ".")
    try:
        os.close(fd)
        try:
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

    csv_path = os.path.abspath(csv_path)
    if not os.path.isfile(csv_path):
        print("Error: Combined CSV not found:", csv_path, file=sys.stderr)
        raise SystemExit(1)

    with open(csv_path, newline="") as f:
        reader = csv.DictReader(f)
        if not reader.fieldnames or "Model" not in reader.fieldnames or "Subdomain" not in reader.fieldnames:
            print("Error: Combined CSV must have Model and Subdomain columns:", csv_path, file=sys.stderr)
            raise SystemExit(1)
        col_ids = [c for c in reader.fieldnames if c not in ("Model", "Subdomain")]
        rows = list(reader)

    if not rows or not col_ids:
        print("Error: Combined CSV is empty or has no structure columns:", csv_path, file=sys.stderr)
        raise SystemExit(1)

    # Preserve CSV row and column order (columns follow row order as written by write_combined_rmsd_csv).
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

    finite = arr[np.isfinite(arr)]
    pos = finite[finite > 0]
    if len(pos) > 0:
        data_min = float(np.min(pos))
        data_max = float(np.max(pos))
    elif len(finite) > 0:
        data_min = float(np.min(finite))
        data_max = float(np.max(finite))
    else:
        print("No finite RMSD values in combined CSV to plot.", file=sys.stderr)
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

    arr_disp = arr.copy()
    for i, m in enumerate(model_names):
        for j, cid in enumerate(col_ids):
            if cid == m or not np.isfinite(arr[i, j]):
                arr_disp[i, j] = disp_vmin - 1.0

    row_labels = [f"{r.get('Subdomain', '')} | {r.get('Model', '')}" for r in rows]

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

    fig, ax = plt.subplots(figsize=(max(8, n_c * 0.35), max(6, n_r * 0.22)))
    if _heatmap_vector_format(fmt):
        sm = _draw_heatmap_vector_cells(ax, arr_disp, cmap_obj, disp_vmin, disp_vmax)
        plt.colorbar(sm, ax=ax, label="RMSD (Å)", extend="min")
    else:
        im = ax.imshow(
            arr_disp,
            aspect="auto",
            cmap=cmap_obj,
            interpolation="nearest",
            vmin=disp_vmin,
            vmax=disp_vmax,
        )
        plt.colorbar(im, ax=ax, label="RMSD (Å)", extend="min")
    ax.set_xticks(range(n_c))
    ax.set_yticks(range(n_r))
    ax.set_xticklabels(col_ids, rotation=45, ha="right", fontsize=x_fs)
    ax.set_yticklabels(row_labels, fontsize=y_fs)
    ax.set_title(title)
    plt.tight_layout()

    fd, tmp_path = tempfile.mkstemp(suffix=ext, prefix="rmsd_combined_heatmap_", dir=plot_dir or ".")
    try:
        os.close(fd)
        try:
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

    ids, matrix = rmsd_file_to_matrix(path)
    if not ids or not matrix:
        print("Error: No RMSD pairs found in", path, file=sys.stderr)
        sys.exit(1)

    print("Loaded", len(ids), "models, matrix", len(ids), "x", len(ids))

    ids, matrix = _apply_label_order(args, ids, matrix, context=os.path.basename(path))

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
        plot_heatmap(
            ids,
            matrix,
            args.plot,
            title=args.title,
            cmap=args.cmap,
            vmin=args.vmin,
            vmax=args.vmax,
        )


def run_scan_dir(args: argparse.Namespace) -> None:
    root = os.path.abspath(args.scan_dir)
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
        hdir = os.path.abspath(args.heatmap_dir)
        os.makedirs(hdir, exist_ok=True)

    tables: list[tuple[str, list[str], list[list[float]]]] = []
    for path in sources:
        label = subdomain_label_from_rmsd_path(path)
        ids, matrix = rmsd_file_to_matrix(path)
        if not ids or not matrix:
            print(f"Warning: skip (no pairs): {path}", file=sys.stderr)
            continue
        ids, matrix = _apply_label_order(args, ids, matrix, context=label)
        out_csv = os.path.join(os.path.dirname(path), f"rmsd_table_{label}.csv")
        write_rmsd_csv(ids, matrix, out_csv)
        print("Wrote", out_csv)
        tables.append((label, ids, matrix))

    if not tables:
        print("Error: No valid RMSD tables produced.", file=sys.stderr)
        sys.exit(1)

    shared_vmin: float | None = None
    shared_vmax: float | None = None
    if args.heatmap_dir and args.shared_heatmap_scale:
        sv = disp_vmin_vmax_for_combined_tables(tables, args.vmin, args.vmax)
        if sv is not None:
            shared_vmin, shared_vmax = sv
            print(
                "Using combined-table colour scale for per-log heatmaps: vmin=%s vmax=%s"
                % (shared_vmin, shared_vmax),
                file=sys.stderr,
            )
        else:
            print(
                "Warning: could not derive combined scale; per-log heatmaps use local auto scale.",
                file=sys.stderr,
            )

    if args.heatmap_dir:
        hf = args.heatmap_format
        hdir_abs = os.path.abspath(args.heatmap_dir)
        for label, ids, matrix in tables:
            hp = os.path.join(hdir_abs, f"rmsd_heatmap_{label}.{hf}")
            pv = shared_vmin if shared_vmin is not None else args.vmin
            px = shared_vmax if shared_vmax is not None else args.vmax
            plot_heatmap(
                ids,
                matrix,
                hp,
                title=f"{args.title} ({label})" if args.title else f"RMSD heatmap ({label})",
                cmap=args.cmap,
                vmin=pv,
                vmax=px,
            )

    if not args.no_combined:
        combined_path = args.combined_csv or os.path.join(root, "combined_rmsd_table.csv")
        combined_path = os.path.abspath(combined_path)
        write_combined_rmsd_csv(tables, combined_path)
        print("Wrote", combined_path)

        want_combined_hm = (args.heatmap_dir or args.combined_heatmap) and not args.no_combined_heatmap
        if want_combined_hm:
            hf = args.heatmap_format
            if args.combined_heatmap:
                ch_out = ensure_batch_heatmap_path(
                    os.path.abspath(args.combined_heatmap),
                    hf,
                )
            else:
                ch_out = os.path.join(
                    os.path.abspath(args.heatmap_dir),
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
  combined CSV (default DIR/combined_rmsd_table.csv): rows sorted by Subdomain then Model; structure
  columns follow that row order (first-seen model per sorted rows). Columns Model + those names +
  Subdomain. Empty cells where a subdomain has no pair. Use --no-combined to skip the merged file. For figures in batch mode use --heatmap-dir
  (one PNG per log plus combined_rmsd_heatmap.png from the merged CSV), or set --combined-heatmap
  PATH to write only the combined figure (file extension sets format if present; else --heatmap-format).
  Use --heatmap-format png|pdf|svg for auto-named batch outputs (default png). Use --no-combined-heatmap
  to skip the combined figure when using --heatmap-dir. Single-file mode uses --plot (extension sets format).
  Use --shared-heatmap-scale so per-log heatmaps use the same vmin/vmax as the combined heatmap (requires
  --heatmap-dir; honors --vmin/--vmax when set).
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
        "<fmt> is set by --heatmap-format (default png).",
    )
    ap.add_argument(
        "--heatmap-format",
        choices=("png", "pdf", "svg"),
        default="png",
        help="Batch only: image format extension for auto-named heatmaps under --heatmap-dir and for the "
        "default combined path. Ignored when --combined-heatmap ends with .png/.pdf/.svg. "
        "SVG and PDF draw each matrix cell as a vector patch (PNG uses a raster grid); "
        "the colour bar may still use a raster strip.",
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
        "--shared-heatmap-scale",
        action="store_true",
        help="Batch only: with --heatmap-dir, use the same vmin/vmax as the combined heatmap "
        "(from the merged table; honors --vmin/--vmax when set) for every per-log heatmap.",
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
        "SVG/PDF use a vector cell per matrix entry; the colour bar may still use a raster strip.",
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
        help="Heatmap colour scale floor (Å); default: min of positive off-diagonal RMSDs in that matrix "
        "(or in the merged table when using --shared-heatmap-scale).",
    )
    ap.add_argument(
        "--vmax",
        type=float,
        default=None,
        metavar="VAL",
        help="Heatmap colour scale ceiling (Å); default: max of positive off-diagonal RMSDs in that matrix "
        "(or in the merged table when using --shared-heatmap-scale).",
    )
    args = ap.parse_args()

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
    if args.shared_heatmap_scale and not args.scan_dir:
        ap.error("--shared-heatmap-scale is only used with --scan-dir")
    if args.shared_heatmap_scale and not args.heatmap_dir:
        ap.error("--shared-heatmap-scale requires --heatmap-dir")

    if args.scan_dir:
        run_scan_dir(args)
    else:
        run_single_file(args)


if __name__ == "__main__":
    main()
