#!/usr/bin/env python3
"""
DALI-style all-vs-all structural comparison, ranking, and phylogeny.

This script mirrors the typical DALI workflow at the level of downstream analysis:
1) Use an all-vs-all similarity measure (DALI Z-scores) to rank structures.
2) Convert similarities into a distance matrix.
3) Build a dendrogram / phylogeny from that distance matrix (prefer NJ).

It does NOT call the public DALI server. Instead, it consumes pairwise Z-scores
that you already computed with DALI (e.g. DaliLite or other DALI workflows).

Input:
  - Pairwise Z-score table (TSV/CSV/space-delimited): model_01, model_02, zscore
    Example lines:
      model_01   model_02   12.3
      model_01,model_03,8.7
    Comment lines beginning with '#' are ignored.

Outputs:
  - A ranking table (CSV) summarising average/max Z per structure.
  - A Newick tree inferred from a Z-score-derived distance matrix.
  - Optional plot if ete3 or biopython+matplotlib are installed.

Distance transforms (Z -> distance):
  - inv:        d = 1 / (1 + Z)
  - maxminus:   d = (Zmax - Z) / max(1e-9, Zmax)   (normalised to [0,1] roughly)
  - exp:        d = exp(-Z / scale)                (scale configurable)

Notes:
  - DALI Z-scores are similarities (higher = more similar). Trees require distances.
  - Neighbor-joining (NJ) is used when scikit-bio or Biopython are available.
    If not, a pure-Python UPGMA fallback is used.
"""

from __future__ import annotations

import argparse
import csv
import math
import os
import re
import sys
from pathlib import Path

_DALI_PKG = os.path.dirname(os.path.abspath(__file__))
if _DALI_PKG not in sys.path:
    sys.path.insert(0, _DALI_PKG)

_REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if _REPO_ROOT not in sys.path:
    sys.path.insert(0, _REPO_ROOT)

from cli_log import add_log_args, setup_log_from_args


def read_pairwise_zscores(path: str) -> dict[tuple[str, str], float]:
    """
    Read pairwise Z-scores from a delimited text file.
    Accepts comma, tab, or whitespace separators. Expects 3 fields: a, b, z.
    """
    pairs: dict[tuple[str, str], float] = {}
    with open(path, "r", newline="") as f:
        for raw in f:
            line = raw.strip()
            if not line or line.startswith("#"):
                continue
            # Split on comma, tab, or whitespace (robust to mixed separators).
            fields = [x for x in re.split(r"[,\t ]+", line) if x]
            if len(fields) < 3:
                continue
            a, b = fields[0].strip(), fields[1].strip()
            try:
                z = float(fields[2])
            except ValueError:
                continue
            if a == b:
                continue
            pairs[(a, b)] = z
            pairs[(b, a)] = z
    return pairs


def zscores_to_distance_matrix(
    zscores: dict[tuple[str, str], float],
    transform: str,
    exp_scale: float = 10.0,
) -> tuple[list[str], list[list[float]], float]:
    """
    Build a symmetric distance matrix from pairwise Z-scores.
    Returns (labels, matrix, zmax).
    Missing pairs are filled with the minimum observed Z (weak similarity).
    """
    labels = sorted(set(x for ab in zscores.keys() for x in ab))
    n = len(labels)
    if n == 0:
        return [], [], 0.0

    observed = [v for (a, b), v in zscores.items() if a < b]
    zmax = max(observed) if observed else 0.0
    zmin = min(observed) if observed else 0.0

    def to_dist(z: float) -> float:
        if transform == "inv":
            return 1.0 / (1.0 + max(0.0, z))
        if transform == "maxminus":
            if zmax <= 0:
                return 1.0
            return max(0.0, (zmax - z) / max(1e-9, zmax))
        if transform == "exp":
            s = max(1e-9, float(exp_scale))
            return float(math.exp(-max(0.0, z) / s))
        raise ValueError(f"Unknown transform: {transform}")

    matrix = [[0.0] * n for _ in range(n)]
    for i, a in enumerate(labels):
        for j, b in enumerate(labels):
            if i == j:
                continue
            z = zscores.get((a, b))
            if z is None:
                z = zmin
            matrix[i][j] = matrix[j][i] = to_dist(float(z))
    return labels, matrix, zmax


def rank_structures(zscores: dict[tuple[str, str], float]) -> list[dict[str, object]]:
    """
    Rank structures by average Z-score to all others (descending).
    Also reports max Z and count of comparisons.
    """
    labels = sorted(set(x for ab in zscores.keys() for x in ab))
    out = []
    for a in labels:
        vals = []
        for b in labels:
            if a == b:
                continue
            z = zscores.get((a, b))
            if z is not None:
                vals.append(float(z))
        if vals:
            out.append(
                {
                    "label": a,
                    "n_pairs": len(vals),
                    "avg_z": sum(vals) / len(vals),
                    "max_z": max(vals),
                }
            )
        else:
            out.append({"label": a, "n_pairs": 0, "avg_z": 0.0, "max_z": 0.0})
    out.sort(key=lambda r: (float(r["avg_z"]), float(r["max_z"])), reverse=True)
    return out


def _ensure_numeric_matrix(matrix: list[list[float]]) -> list[list[float]]:
    flat = [x for row in matrix for x in row if not (isinstance(x, float) and math.isnan(x))]
    fill = (max(flat) * 1.5) if flat else 1.0
    out = []
    for i, row in enumerate(matrix):
        out.append(
            [
                0.0 if i == j else (fill if (isinstance(x, float) and math.isnan(x)) else float(x))
                for j, x in enumerate(row)
            ]
        )
    return out


def _upgma_newick(labels: list[str], matrix: list[list[float]]) -> str:
    """Pure-Python UPGMA fallback; returns Newick."""
    n = len(labels)
    if n == 0:
        return "();"
    if n == 1:
        return f"({labels[0]});"
    if n == 2:
        d = matrix[0][1]
        return f"({labels[0]}:{d/2:.6f},{labels[1]}:{d/2:.6f});"

    cluster_size = [1] * n + [0] * (n - 1)
    height = [0.0] * (2 * n - 1)
    parent = [None] * (2 * n - 1)
    dist: dict[tuple[int, int], float] = {}
    for i in range(n):
        for j in range(i + 1, n):
            dist[(i, j)] = float(matrix[i][j])

    next_node = n
    active = set(range(n))
    for _ in range(n - 1):
        best = float("inf")
        best_pair = None
        for i in active:
            for j in active:
                if i >= j:
                    continue
                d = dist.get((i, j) if i < j else (j, i), float("inf"))
                if d < best:
                    best = d
                    best_pair = (i, j)
        if best_pair is None:
            break
        i, j = best_pair
        ni, nj = cluster_size[i], cluster_size[j]
        h = best / 2.0
        height[next_node] = h
        cluster_size[next_node] = ni + nj
        parent[i] = parent[j] = next_node
        active.discard(i)
        active.discard(j)
        active.add(next_node)
        for k in active:
            if k == next_node:
                continue
            d_ik = dist.get((i, k) if i < k else (k, i), 0.0)
            d_jk = dist.get((j, k) if j < k else (k, j), 0.0)
            new_d = (d_ik * ni + d_jk * nj) / (ni + nj)
            dist[(k, next_node) if k < next_node else (next_node, k)] = new_d
        next_node += 1

    def newick_node(node: int) -> str:
        if node < n:
            return labels[node]
        children = [c for c in range(2 * n - 1) if parent[c] == node]
        if len(children) != 2:
            return "();"
        L, R = children[0], children[1]
        h = height[node]
        bl_L = h - (height[L] if L >= n else 0.0)
        bl_R = h - (height[R] if R >= n else 0.0)
        return f"({newick_node(L)}:{bl_L:.6f},{newick_node(R)}:{bl_R:.6f})"

    root_node = next_node - 1
    return newick_node(root_node) + ";"


def build_tree_newick(labels: list[str], matrix: list[list[float]], root: str | None, midpoint_root: bool) -> str:
    """
    Prefer neighbour-joining if available; otherwise UPGMA fallback.
    Rooting:
      - scikit-bio: supports midpoint and root_at(outgroup set)
      - biopython: supports outgroup rooting; midpoint rooting not implemented here
      - fallback: produces rooted UPGMA by construction; --root is ignored
    """
    matrix_clean = _ensure_numeric_matrix(matrix)

    try:
        import skbio

        dm = skbio.DistanceMatrix(matrix_clean, labels)
        tree = skbio.tree.nj(dm)
        if root and root in labels:
            tree = tree.root_at({root})
        elif midpoint_root:
            tree = tree.root_at_midpoint()
        return str(tree)
    except ImportError:
        pass

    try:
        from Bio.Phylo.TreeConstruction import DistanceMatrix, DistanceTreeConstructor
        from io import StringIO
        from Bio import Phylo

        dm = DistanceMatrix(labels, matrix_clean)
        tree = DistanceTreeConstructor().nj(dm)
        if root and root in labels:
            for c in tree.get_terminals():
                if c.name == root:
                    tree.root_with_outgroup(c)
                    break
        s = StringIO()
        Phylo.write(tree, s, "newick")
        return s.getvalue()
    except ImportError:
        pass

    return _upgma_newick(labels, matrix_clean)


def write_ranking_csv(rows: list[dict[str, object]], out_csv: str) -> None:
    Path(os.path.dirname(out_csv) or ".").mkdir(parents=True, exist_ok=True)
    with open(out_csv, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["label", "n_pairs", "avg_z", "max_z"])
        w.writeheader()
        for r in rows:
            w.writerow(r)


def plot_tree(newick: str, out_path: str, title: str) -> None:
    try:
        from ete3 import Tree, TreeStyle, TextFace

        t = Tree(newick)
        ts = TreeStyle()
        ts.title.add_face(TextFace(title, fsize=12), column=0)
        t.render(out_path, tree_style=ts)
        return
    except ImportError:
        pass

    try:
        from Bio import Phylo
        from io import StringIO
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        tree = Phylo.read(StringIO(newick), "newick")
        fig, ax = plt.subplots(1, 1, figsize=(10, 8))
        Phylo.draw(tree, axes=ax, do_show=False)
        ax.set_title(title)
        plt.savefig(out_path, dpi=150, bbox_inches="tight")
        plt.close()
    except ImportError:
        raise SystemExit("Install ete3 or biopython+matplotlib to generate --plot output.")


class _ArgHelp(argparse.ArgumentDefaultsHelpFormatter, argparse.RawDescriptionHelpFormatter):
    """Defaults in --help plus preserved epilog newlines."""


def main() -> None:
    ap = argparse.ArgumentParser(
        description="DALI-style analysis from all-vs-all Z-scores: ranking + tree (Newick).",
        formatter_class=_ArgHelp,
        epilog="""
Examples (from repository root):
  python ranking/dali_phylogeny.py pairwise_zscores.tsv -o dali_tree.nwk
  python ranking/dali_phylogeny.py pairwise_zscores.tsv --output-ranking dali_ranking.csv --plot tree.pdf
""",
    )
    ap.add_argument(
        "zscores",
        help="Pairwise Z-score table: label_a label_b zscore (TSV/CSV/space-delimited).",
    )
    ap.add_argument(
        "-o",
        "--output-tree",
        default="dali_tree.nwk",
        help="Output Newick tree path.",
    )
    ap.add_argument(
        "--output-ranking",
        default="dali_ranking.csv",
        help="Output ranking CSV path.",
    )
    ap.add_argument(
        "--transform",
        choices=["inv", "maxminus", "exp"],
        default="inv",
        help="Z-score to distance transform.",
    )
    ap.add_argument(
        "--exp-scale",
        type=float,
        default=10.0,
        help="Scale for exp transform: d=exp(-Z/scale).",
    )
    ap.add_argument(
        "--root",
        default=None,
        help="Outgroup label to root on (if supported by tree builder).",
    )
    ap.add_argument(
        "--no-midpoint-root",
        action="store_true",
        help="Disable midpoint rooting (NJ only).",
    )
    ap.add_argument(
        "--plot",
        default=None,
        metavar="PATH",
        help="Write a tree plot to PATH.",
    )
    add_log_args(ap)
    args = ap.parse_args()
    summary_log = setup_log_from_args(args, script_path=__file__, inputs=[args.zscores], pattern=None)
    if summary_log is not None:
        summary_log.task("DALI phylogeny from pairwise Z-scores (dali_phylogeny.py)")

    if not os.path.isfile(args.zscores):
        if summary_log is not None:
            summary_log.error(f"File not found: {args.zscores}")
        ap.error(f"File not found: {args.zscores}")
    _main_after_input(args, summary_log=summary_log)


def _main_after_input(args: argparse.Namespace, summary_log=None) -> None:
    z = read_pairwise_zscores(args.zscores)
    if not z:
        if summary_log is not None:
            summary_log.error("No pairwise Z-scores parsed from input.")
        raise SystemExit("No pairwise Z-scores parsed from input.")

    ranking = rank_structures(z)
    write_ranking_csv(ranking, args.output_ranking)

    labels, matrix, zmax = zscores_to_distance_matrix(z, args.transform, args.exp_scale)
    newick = build_tree_newick(labels, matrix, args.root, midpoint_root=not args.no_midpoint_root)
    Path(os.path.dirname(args.output_tree) or ".").mkdir(parents=True, exist_ok=True)
    with open(args.output_tree, "w") as f:
        f.write(newick)

    print(f"Parsed {len(labels)} structures; max Z observed: {zmax:.3f}")
    if summary_log is not None:
        summary_log.task(f"Parsed {len(labels)} structures; max Z observed: {zmax:.3f}")
    print(f"Wrote ranking: {args.output_ranking}")
    print(f"Wrote tree: {args.output_tree}")
    if args.plot:
        plot_tree(newick, args.plot, title="DALI-style tree (from Z-scores)")
        print(f"Wrote plot: {args.plot}")


if __name__ == "__main__":
    main()

