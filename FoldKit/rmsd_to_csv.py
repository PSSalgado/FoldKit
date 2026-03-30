#!/usr/bin/env python3
"""
Convert pairwise RMSD file (LSQ or SSM, e.g. from all-vs-all) to a square CSV table
suitable for structure_phylogeny.py and create_rmsd_heatmap.R.

Usage:
  python rmsd_to_csv.py path/to/rmsd_SSM_values.txt -o path/to/rmsd_table.csv
  python rmsd_to_csv.py path/to/rmsd_SSM_values.txt --plot heatmap.pdf
  python rmsd_to_csv.py path/to/rmsd_values.txt --order-dir /path/to/root --order-glob 'model_*.pdb'
"""

from __future__ import annotations

import argparse
import csv
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

    # Default color scale: off-diagonal pairwise RMSDs only, and > 0 (exclude diagonal zeros).
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
    # Diagonal: below scale so it maps to the colormap "under" color (white)
    np.fill_diagonal(arr, disp_vmin - 1.0)

    try:
        cmap_obj = plt.get_cmap(cmap).copy()
    except ValueError:
        print("Unknown colormap '%s'. Common: viridis_r, plasma_r, RdYlBu_r, coolwarm, YlOrRd." % cmap, file=sys.stderr)
        cmap_obj = plt.get_cmap("viridis_r").copy()
    cmap_obj.set_under(color="white")

    fig, ax = plt.subplots(figsize=(max(6, n * 0.4), max(5, n * 0.35)))
    im = ax.imshow(arr, aspect="auto", cmap=cmap_obj, interpolation="nearest", vmin=disp_vmin, vmax=disp_vmax)
    ax.set_xticks(range(n))
    ax.set_yticks(range(n))
    ax.set_xticklabels(ids, rotation=45, ha="right")
    ax.set_yticklabels(ids)
    ax.set_title(title)
    plt.colorbar(im, ax=ax, label="RMSD (Å)", extend="min")
    plt.tight_layout()

    ext = os.path.splitext(out_path)[1].lower()
    fmt = "png" if ext == ".png" else "pdf" if ext == ".pdf" else "svg" if ext == ".svg" else "png"

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


def main() -> None:
    ap = argparse.ArgumentParser(
        description="Convert pairwise RMSD file (LSQ/SSM) to CSV table and optionally plot heatmap."
    )
    ap.add_argument(
        "input",
        help="RMSD file: rmsd_values_*.txt, rmsd_SSM_values*.txt, or rmsd_table_*.csv",
    )
    ap.add_argument(
        "-o", "--output",
        default=None,
        help="Output CSV path (default: same dir as input, rmsd_table_<suffix>.csv)",
    )
    ap.add_argument(
        "--plot",
        metavar="PATH",
        default=None,
        help="Save heatmap figure to PATH (e.g. heatmap.png, heatmap.svg, or heatmap.pdf; PNG/SVG are more reliable on headless systems)",
    )
    ap.add_argument(
        "--title",
        default="RMSD heatmap",
        help="Title for heatmap (default: RMSD heatmap)",
    )
    order_group = ap.add_mutually_exclusive_group()
    order_group.add_argument(
        "--order",
        metavar="FILE_OR_LIST",
        default=None,
        help="Order of labels: path to file (one label per line) or comma-separated list. Applies to CSV and heatmap.",
    )
    order_group.add_argument(
        "--order-dir",
        metavar="DIR",
        default=None,
        help="Build label order from files under DIR (recursive). Labels are basenames without extension; "
        "sorted by relative path with natural numeric order. Use with --order-glob (default *.pdb). "
        "Incompatible with --order.",
    )
    ap.add_argument(
        "--order-glob",
        metavar="PATTERN",
        default=None,
        help="Glob for --order-dir (default when --order-dir is set: *.pdb). Example: model_*.pdb",
    )
    ap.add_argument(
        "--cmap",
        default="viridis_r",
        help="Matplotlib colormap for heatmap (default: viridis_r). Examples: viridis_r, plasma_r, RdYlBu_r, coolwarm, YlOrRd.",
    )
    ap.add_argument(
        "--vmin",
        type=float,
        default=None,
        metavar="VAL",
        help="Heatmap color scale minimum (Å); default: min of positive off-diagonal values.",
    )
    ap.add_argument(
        "--vmax",
        type=float,
        default=None,
        metavar="VAL",
        help="Heatmap color scale maximum (Å); default: max of positive off-diagonal values.",
    )
    args = ap.parse_args()

    path = os.path.abspath(args.input)
    if not os.path.isfile(path):
        print("Error: Input file not found:", path, file=sys.stderr)
        sys.exit(1)

    ids, matrix = rmsd_file_to_matrix(path)
    if not ids or not matrix:
        print("Error: No RMSD pairs found in", path, file=sys.stderr)
        sys.exit(1)

    print("Loaded", len(ids), "models, matrix", len(ids), "x", len(ids))

    if args.order_dir:
        glob_pat = args.order_glob or "*.pdb"
        root_abs = os.path.abspath(args.order_dir)
        if not os.path.isdir(root_abs):
            print("Error: --order-dir is not a directory:", root_abs, file=sys.stderr)
            sys.exit(1)
        order_list = order_labels_from_scan(root_abs, glob_pat)
        if order_list:
            ids, matrix = reorder_matrix(ids, matrix, order_list)
            print("Applied --order-dir", root_abs, "glob", repr(glob_pat), ":", len(order_list), "labels")
        else:
            print(
                "Warning: no files matched under --order-dir; using default matrix order.",
                file=sys.stderr,
            )
    elif args.order:
        order_list = parse_order(args.order)
        if order_list:
            ids, matrix = reorder_matrix(ids, matrix, order_list)
            print("Applied --order:", len(order_list), "labels")
        else:
            print("Warning: --order empty or file not found; using default order.", file=sys.stderr)
    elif args.order_glob:
        print("Error: --order-glob requires --order-dir.", file=sys.stderr)
        sys.exit(1)

    # CSV output (-o is for the CSV table; --plot is for the heatmap image)
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


if __name__ == "__main__":
    main()
