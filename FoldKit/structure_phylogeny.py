#!/usr/bin/env python3
"""
Structure-based phylogenetic tree from RMSD distances.

Reads pairwise RMSD values (from Coot LSQ/SSM superpositions or from rmsd_table_*.csv),
builds a symmetric distance matrix, and constructs a phylogenetic tree (rooted or unrooted) using
neighbour-joining (NJ). Optionally computes pairwise distances from PDB files using
TM-align (if installed).

Input formats supported:
  - rmsd_values_*.txt (LSQ): "Alignment: Aligning A to B" then "core rmsd achieved: X.XX"
  - rmsd_SSM_values*.txt (SSM): "Superposing … onto …" then "INFO: core rmsd" block with number
  - Other *.txt names: format is detected from file content (SSM vs LSQ) when the basename is ambiguous
  - rmsd_table_*.csv: CSV with Model column and one column per model (from rmsd_to_csv.py)
"""

from __future__ import annotations

import argparse
import os
import re
import subprocess
import sys
from collections import defaultdict
from pathlib import Path
import math
import csv

from cli_log import add_log_args, setup_log_from_args


def _parse_rmsd_number(line: str) -> float | None:
    m = re.search(r"(?:INFO:\s*)?core\s+rmsd(?:\s+achieved)?\s*:\s*([\d.]+)", line, re.I)
    if m:
        return float(m.group(1))
    m = re.search(r"rmsd\s*[=:]?\s*([\d.]+)", line, re.I)
    return float(m.group(1)) if m else None


def parse_lsq_rmsd_txt(path: str) -> dict[tuple[str, str], float]:
    alignments = {}
    current_pair = None
    with open(path, "r") as f:
        for line in f:
            line = line.strip()
            if line.startswith("Alignment:"):
                text = line.replace("Alignment:", "").strip()
                match = re.match(r"Aligning\s+(.+?)\s+to\s+(.+?)$", text)
                if match:
                    a, b = match.group(1).strip(), match.group(2).strip()
                    current_pair = (a, b)
            elif current_pair and "core rmsd" in line:
                rmsd = _parse_rmsd_number(line)
                if rmsd is not None:
                    a, b = current_pair
                    alignments[(a, b)] = alignments[(b, a)] = rmsd
                current_pair = None
    return alignments


def parse_ssm_rmsd_txt(path: str) -> dict[tuple[str, str], float]:
    alignments = {}
    current_pair = None
    with open(path, "r") as f:
        for line in f:
            s = line.strip()
            if "Superposing" in s and " onto " in s:
                # "Superposing model (chain X) onto reference (chain Y)" or
                # "  Superposing model onto ref_name" (all-vs-all)
                m = re.match(r"Superposing\s+(.+?)\s+\([^)]*\)\s+onto\s+(.+?)(?:\s+\(|$)", s)
                if m:
                    current_pair = (m.group(1).strip(), m.group(2).strip())
                else:
                    m = re.match(r"Superposing\s+(.+?)\s+onto\s+(.+)", s)
                    current_pair = (m.group(1).strip(), m.group(2).strip()) if m else None
            elif "Aligning" in s and " to " in s:
                m = re.match(r"Aligning\s+(.+?)\s+to\s+(.+?)$", s)
                current_pair = (m.group(1).strip(), m.group(2).strip()) if m else None
            elif current_pair and ("INFO: core rmsd" in s or "core rmsd" in s):
                rmsd = _parse_rmsd_number(s)
                if rmsd is not None:
                    a, b = current_pair
                    alignments[(a, b)] = alignments[(b, a)] = rmsd
                current_pair = None
    return alignments


def parse_rmsd_csv(path: str) -> dict[tuple[str, str], float]:
    import csv as csv_module

    alignments = {}
    with open(path, "r", newline="") as f:
        rows = list(csv_module.DictReader(f))
    if not rows or "Model" not in rows[0]:
        return alignments
    models = [r["Model"] for r in rows]
    for i, row in enumerate(rows):
        for j, m in enumerate(models):
            if i == j:
                continue
            val = row.get(m, "").strip()
            if val and val != "-":
                try:
                    rmsd = float(val)
                    alignments[(models[i], m)] = alignments[(m, models[i])] = rmsd
                except ValueError:
                    pass
    return alignments


def _natural_sort_key(s: str):
    """Sort key so 'model_02' comes before 'model_10' (natural/numeric order)."""
    import re
    return [int(x) if x.isdigit() else x.lower() for x in re.split(r"(\d+)", s)]


def alignments_to_matrix(
    alignments: dict[tuple[str, str], float],
) -> tuple[list[str], list[list[float]]]:
    ids = sorted(set(a for pair in alignments for a in pair), key=_natural_sort_key)
    n = len(ids)
    matrix = [[0.0] * n for _ in range(n)]
    for i, a in enumerate(ids):
        for j, b in enumerate(ids):
            if i == j:
                continue
            v = alignments.get((a, b)) or alignments.get((b, a))
            if v is not None:
                matrix[i][j] = matrix[j][i] = v
            else:
                vals = [
                    alignments.get((a, x)) or alignments.get((x, a))
                    for x in ids
                    if x != a
                    and (
                        alignments.get((a, x)) is not None
                        or alignments.get((x, a)) is not None
                    )
                ]
                matrix[i][j] = matrix[j][i] = (
                    (sum(vals) / len(vals)) if vals else float("nan")
                )
    return ids, matrix


def _rmsd_to_similarity(rmsd: float, mode: str, tau: float) -> float:
    """
    Convert RMSD distance to a similarity.

    - neg_rmsd: s = -rmsd
    - exp:      s = exp(-rmsd/tau)
    """
    d = float(rmsd)
    if mode == "neg_rmsd":
        return -d
    if mode == "exp":
        t = max(1e-9, float(tau))
        return float(math.exp(-d / t))
    raise ValueError(f"Unknown similarity mode: {mode}")


def _matrix_to_similarity(ids: list[str], rmsd_matrix: list[list[float]], mode: str, tau: float) -> list[list[float]]:
    n = len(ids)
    sim = [[0.0] * n for _ in range(n)]
    for i in range(n):
        for j in range(n):
            if i == j:
                sim[i][j] = 0.0
            else:
                sim[i][j] = _rmsd_to_similarity(rmsd_matrix[i][j], mode, tau)
    return sim


def _per_query_zscores(ids: list[str], sim: list[list[float]]) -> list[list[float]]:
    """
    Compute per-query z-scores from a similarity matrix.
    For each i, z_ij = (s_ij - mean_i) / std_i over j != i.
    """
    n = len(ids)
    z = [[0.0] * n for _ in range(n)]
    for i in range(n):
        vals = [sim[i][j] for j in range(n) if j != i]
        if not vals:
            continue
        mu = sum(vals) / len(vals)
        var = sum((v - mu) ** 2 for v in vals) / len(vals)
        sd = math.sqrt(var)
        if sd <= 1e-12:
            continue
        for j in range(n):
            if i == j:
                continue
            z[i][j] = (sim[i][j] - mu) / sd
    return z


def _symmetrize_matrix(m: list[list[float]]) -> list[list[float]]:
    n = len(m)
    out = [[0.0] * n for _ in range(n)]
    for i in range(n):
        for j in range(n):
            if i == j:
                continue
            out[i][j] = out[j][i] = (float(m[i][j]) + float(m[j][i])) / 2.0
    return out


def _z_to_distance(z: float, transform: str, scale: float, zmax: float) -> float:
    """
    Convert a (possibly negative) z-score into a distance.
    Negative z indicates below-average similarity; we clamp at 0 for transforms that expect similarity.
    """
    zp = max(0.0, float(z))
    if transform == "inv":
        return 1.0 / (1.0 + zp)
    if transform == "exp":
        s = max(1e-9, float(scale))
        return float(math.exp(-zp / s))
    if transform == "maxminus":
        if zmax <= 1e-12:
            return 1.0
        return max(0.0, (zmax - zp) / max(1e-9, zmax))
    raise ValueError(f"Unknown z-distance transform: {transform}")


def zscore_distance_matrix(
    ids: list[str],
    rmsd_matrix: list[list[float]],
    z_mode: str,
    similarity: str,
    tau: float,
    z_transform: str,
    z_scale: float,
) -> tuple[list[list[float]], list[dict[str, object]]]:
    """
    Build a distance matrix for tree construction.

    If z_mode == 'off': returns the original RMSD matrix and empty ranking.
    If z_mode == 'per_query': RMSD -> similarity -> per-query z -> symmetrise -> z->distance.

    Returns: (distance_matrix, ranking_rows)
    """
    if z_mode == "off":
        return rmsd_matrix, []

    if z_mode != "per_query":
        raise ValueError(f"Unsupported z_mode: {z_mode}")

    sim = _matrix_to_similarity(ids, rmsd_matrix, similarity, tau)
    z_asym = _per_query_zscores(ids, sim)
    z_sym = _symmetrize_matrix(z_asym)

    # zmax over off-diagonal (after clamp) for maxminus normalisation
    zmax = 0.0
    n = len(ids)
    for i in range(n):
        for j in range(i + 1, n):
            zmax = max(zmax, max(0.0, float(z_sym[i][j])))

    dist = [[0.0] * n for _ in range(n)]
    for i in range(n):
        for j in range(n):
            if i == j:
                continue
            dist[i][j] = dist[j][i] = _z_to_distance(z_sym[i][j], z_transform, z_scale, zmax)

    # Ranking by average per-query z (asymmetric, like DALI query-centric ranking)
    ranking = []
    for i, label in enumerate(ids):
        vals = [float(z_asym[i][j]) for j in range(n) if j != i]
        if vals:
            ranking.append(
                {
                    "label": label,
                    "n_pairs": len(vals),
                    "avg_z": sum(vals) / len(vals),
                    "max_z": max(vals),
                }
            )
        else:
            ranking.append({"label": label, "n_pairs": 0, "avg_z": 0.0, "max_z": 0.0})
    ranking.sort(key=lambda r: (float(r["avg_z"]), float(r["max_z"])), reverse=True)
    return dist, ranking


def write_ranking_csv(rows: list[dict[str, object]], out_csv: str) -> None:
    if not rows:
        return
    Path(os.path.dirname(out_csv) or ".").mkdir(parents=True, exist_ok=True)
    with open(out_csv, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["label", "n_pairs", "avg_z", "max_z"])
        w.writeheader()
        for r in rows:
            w.writerow(r)


def compute_rmsd_matrix_from_pdbs(
    pdb_paths: list[str],
) -> tuple[list[str], list[list[float]]]:
    labels = [os.path.splitext(os.path.basename(p))[0] for p in pdb_paths]
    n = len(labels)
    matrix = [[0.0] * n for _ in range(n)]
    try:
        subprocess.run(["TMalign", "-h"], capture_output=True, check=False)
    except (FileNotFoundError, OSError):
        raise SystemExit(
            "TM-align not found. Install it or use precomputed RMSD input."
        )
    for i in range(n):
        for j in range(i + 1, n):
            try:
                out = subprocess.run(
                    ["TMalign", pdb_paths[i], pdb_paths[j]],
                    capture_output=True,
                    text=True,
                    timeout=120,
                )
                for line in out.stdout.splitlines():
                    if "TM-score=" in line or "TM-score =" in line:
                        m = re.search(r"TM-score\s*=\s*([\d.]+)", line)
                        if m:
                            tm = float(m.group(1))
                            matrix[i][j] = matrix[j][i] = max(0.0, 1.0 - tm)
                            break
                else:
                    matrix[i][j] = matrix[j][i] = 1.0
            except (subprocess.TimeoutExpired, ValueError):
                matrix[i][j] = matrix[j][i] = 1.0
    return labels, matrix


def _ensure_numeric_matrix(
    ids: list[str], matrix: list[list[float]]
) -> list[list[float]]:
    import math

    flat = [
        x
        for row in matrix
        for x in row
        if not (isinstance(x, float) and math.isnan(x))
    ]
    fill = (max(flat) * 1.5) if flat else 1.0
    return [
        [
            0.0
            if i == j
            else (fill if isinstance(x, float) and math.isnan(x) else x)
            for j, x in enumerate(row)
        ]
        for i, row in enumerate(matrix)
    ]


def _biopython_path_to_leaf(clade, leaf_name, path_so_far):
    """Path from root to leaf as list of (clade, branch_length)."""
    bl = clade.branch_length if clade.branch_length is not None else 0.0
    here = path_so_far + [(clade, bl)]
    if clade.name == leaf_name:
        return here
    for child in clade.clades:
        r = _biopython_path_to_leaf(child, leaf_name, here)
        if r is not None:
            return r
    return None


def _biopython_path_between_leaves(tree, name_a, name_b):
    """Path from leaf name_a to leaf name_b as list of (clade, branch_length)."""
    path_a = _biopython_path_to_leaf(tree.root, name_a, [])
    path_b = _biopython_path_to_leaf(tree.root, name_b, [])
    if path_a is None or path_b is None:
        return None
    # Find LCA index (last common clade)
    k = 0
    while k < min(len(path_a), len(path_b)) and path_a[k][0] is path_b[k][0]:
        k += 1
    # Path A to B: reverse path_a from A up to LCA, include LCA, then path_b from LCA to B.
    # path_a[k-1] is (LCA, bl_lca); path_a[k:] is LCA->...->A; path_b[k:] is LCA->...->B.
    # Include LCA so the edge (LCA's parent -> LCA) is counted; otherwise total path length is wrong.
    rev_up = list(reversed(path_a[k:]))  # (A, bl), ..., (child_toward_A, bl)
    lca_entry = path_a[k - 1]  # (LCA, bl_lca)
    down = path_b[k:]  # (child_toward_B, bl), ..., (B, bl)
    path_ab = [(rev_up[i][0], rev_up[i][1]) for i in range(len(rev_up))]
    path_ab.append(lca_entry)
    path_ab.extend((down[i][0], down[i][1]) for i in range(len(down)))
    return path_ab


def _biopython_midpoint_root(tree, ids, matrix_clean):
    """Root the Bio.Phylo tree at the midpoint of the two farthest leaves."""
    n = len(ids)
    if n < 3:
        return
    # Pair with max distance in matrix
    best = -1.0
    idx_a, idx_b = 0, 1
    for i in range(n):
        for j in range(i + 1, n):
            d = matrix_clean[i][j]
            if d > best:
                best = d
                idx_a, idx_b = i, j
    name_a, name_b = ids[idx_a], ids[idx_b]
    path_ab = _biopython_path_between_leaves(tree, name_a, name_b)
    if not path_ab:
        return
    total_len = sum(bl for _, bl in path_ab)
    if total_len <= 0:
        return
    mid = total_len / 2.0
    acc = 0.0
    for i, (clade, bl) in enumerate(path_ab):
        acc += bl
        if acc >= mid:
            # Midpoint is on the edge we just traversed; outgroup is the clade on the A side
            outgroup_clade = path_ab[i - 1][0] if i > 0 else path_ab[0][0]
            tree.root_with_outgroup(outgroup_clade)
            return
    # Fallback: root with last (B)
    tree.root_with_outgroup(path_ab[-1][0])


def build_nj_tree(
    ids: list[str],
    matrix: list[list[float]],
    root: str | None = None,
    midpoint_root: bool = True,
) -> str:
    matrix_clean = _ensure_numeric_matrix(ids, matrix)

    try:
        import skbio

        dm = skbio.DistanceMatrix(matrix_clean, ids)
        tree = skbio.tree.nj(dm)
        if root and root in ids:
            tree = tree.root_at({root})
        elif midpoint_root:
            tree = tree.root_at_midpoint()
        return str(tree)
    except ImportError:
        pass

    try:
        from Bio.Phylo.TreeConstruction import DistanceMatrix, DistanceTreeConstructor
        from Bio.Phylo import BaseTree
        from io import StringIO
        from Bio import Phylo

        # Biopython DistanceMatrix expects lower triangle format: row i has i+1 elements [matrix[i][0], ..., matrix[i][i]]
        n_bio = len(ids)
        matrix_lower = [
            [matrix_clean[i][j] for j in range(i + 1)]
            for i in range(n_bio)
        ]
        dm = DistanceMatrix(ids, matrix_lower)
        constructor = DistanceTreeConstructor()
        tree = constructor.nj(dm)
        if root and root in ids:
            for c in tree.get_terminals():
                if c.name == root:
                    tree.root_with_outgroup(c)
                    break
        elif midpoint_root:
            _biopython_midpoint_root(tree, ids, matrix_clean)
        s = StringIO()
        Phylo.write(tree, s, "newick")
        return s.getvalue()
    except ImportError:
        pass

    try:
        import scipy.cluster.hierarchy as sch
        from scipy.spatial.distance import squareform

        n = len(ids)
        cond = squareform(
            [matrix_clean[i][j] for i in range(n) for j in range(i + 1, n)]
        )
        Z = sch.linkage(cond, method="average")

        def to_newick(i):
            if i < n:
                return ids[i]
            L, R = int(Z[i - n, 0]), int(Z[i - n, 1])
            return "(" + to_newick(L) + "," + to_newick(R) + ")"

        return "(" + to_newick(2 * n - 2) + ");"
    except ImportError:
        pass

    # Pure-Python UPGMA fallback (average-linkage) when no tree library available
    return _upgma_newick(ids, matrix_clean)


def _upgma_newick(ids: list[str], matrix: list[list[float]]) -> str:
    """Build Newick string via UPGMA (average-linkage) in pure Python."""
    n = len(ids)
    if n == 0:
        return "();"
    if n == 1:
        return f"({ids[0]});"
    if n == 2:
        d = matrix[0][1]
        return f"({ids[0]}:{d/2:.6f},{ids[1]}:{d/2:.6f});"

    # Cluster indices: 0..n-1 = leaves, n..2n-2 = internal
    cluster_size = [1] * n + [0] * (n - 1)
    height = [0.0] * (2 * n - 1)
    parent = [None] * (2 * n - 1)
    dist = {}
    for i in range(n):
        for j in range(i + 1, n):
            dist[(i, j)] = matrix[i][j]

    next_node = n
    active = set(range(n))
    for _ in range(n - 1):
        best = float("inf")
        best_pair = None
        for i in active:
            for j in active:
                if i >= j:
                    continue
                key = (i, j) if i < j else (j, i)
                d = dist.get(key, float("inf"))
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
            key_ik = (i, k) if i < k else (k, i)
            key_jk = (j, k) if j < k else (k, j)
            d_ik = dist.get(key_ik, 0.0)
            d_jk = dist.get(key_jk, 0.0)
            new_d = (d_ik * ni + d_jk * nj) / (ni + nj)
            key_new = (k, next_node) if k < next_node else (next_node, k)
            dist[key_new] = new_d
        next_node += 1

    def newick_node(node):
        if node < n:
            return ids[node]
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


def plot_tree(
    newick: str, out_path: str | None, title: str = "Structure-based phylogeny"
) -> None:
    try:
        from ete3 import Tree, TreeStyle, TextFace

        t = Tree(newick)
        ts = TreeStyle()
        ts.title.add_face(TextFace(title, fsize=12), column=0)
        if out_path:
            t.render(out_path, tree_style=ts)
    except ImportError:
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
            if out_path:
                plt.savefig(out_path, dpi=150, bbox_inches="tight")
            plt.close()
        except ImportError:
            if out_path:
                raise SystemExit(
                    "Install ete3 or biopython+matplotlib to generate tree plot."
                )


def _sniff_rmsd_text_sample(sample: str) -> str | None:
    """Return 'ssm_txt' or 'lsq_txt' from the start of a file, or None."""
    if "Superposing" in sample and " onto " in sample:
        return "ssm_txt"
    if "Alignment:" in sample and "Aligning" in sample:
        return "lsq_txt"
    return None


def detect_format(path: str) -> str:
    b = os.path.basename(path).lower()
    if path.lower().endswith(".csv"):
        return "csv"
    if "rmsd_ssm" in b or "ssm_values" in b:
        return "ssm_txt"
    # Peek generic .txt (custom names like rmsd_1m*.txt, coot_log_*.txt, etc.)
    if path.lower().endswith(".txt"):
        try:
            with open(path, "r") as f:
                head = f.read(8192)
            sniffed = _sniff_rmsd_text_sample(head)
            if sniffed:
                return sniffed
        except OSError:
            pass
    return "lsq_txt"


def main() -> None:
    ap = argparse.ArgumentParser(
        description="Structure-based phylogenetic tree from RMSD (Coot LSQ/SSM or CSV)."
    )
    ap.add_argument(
        "input",
        nargs="?",
        help="RMSD file: rmsd_values_*.txt, rmsd_SSM_values*.txt, or rmsd_table_*.csv.",
    )
    ap.add_argument(
        "-o",
        "--output",
        default="structure_tree.nwk",
        help="Output Newick path.",
    )
    ap.add_argument("--plot", metavar="PATH", help="Save tree figure.")
    ap.add_argument(
        "--root", metavar="LABEL", help="Root on this tip (outgroup)."
    )
    ap.add_argument(
        "--unrooted",
        action="store_true",
        help="Build an unrooted tree (no midpoint or outgroup rooting).",
    )
    ap.add_argument(
        "--no-midpoint-root",
        action="store_true",
        help="Do not midpoint-root (same as --unrooted).",
    )
    ap.add_argument(
        "--from-pdb",
        metavar="DIR",
        nargs="?",
        const=".",
        default=None,
        help="Compute distances from PDBs in DIR with TM-align.",
    )
    ap.add_argument(
        "--format",
        choices=["auto", "lsq_txt", "ssm_txt", "csv"],
        default="auto",
        help="RMSD input format (default: auto-detect from filename and content).",
    )
    ap.add_argument(
        "--pseudo-z",
        choices=["off", "per_query"],
        default="off",
        help="Compute per-query pseudo Z-scores from RMSD: RMSD->similarity->z->distance (DALI-like ranking mode).",
    )
    ap.add_argument(
        "--similarity",
        choices=["exp", "neg_rmsd"],
        default="exp",
        help="Similarity transform used before z-scoring. exp: exp(-RMSD/tau); neg_rmsd: -RMSD.",
    )
    ap.add_argument(
        "--tau",
        type=float,
        default=2.0,
        help="Tau for exp similarity: s=exp(-RMSD/tau).",
    )
    ap.add_argument(
        "--zdist",
        choices=["inv", "exp", "maxminus"],
        default="inv",
        help="Transform z-score -> distance for tree building (uses z+ = max(0,z)).",
    )
    ap.add_argument(
        "--zscale",
        type=float,
        default=1.0,
        help="Scale for --zdist exp: d=exp(-z+/zscale).",
    )
    ap.add_argument(
        "--ranking-csv",
        default=None,
        help="Optional: write per-query pseudo-z ranking table to this CSV path.",
    )
    add_log_args(ap)
    args = ap.parse_args()
    setup_log_from_args(args, script_path=__file__, inputs=[getattr(args, "input", "")], pattern=None)

    if args.from_pdb is not None:
        pdb_dir = os.path.abspath(args.from_pdb)
        paths = [
            str(p)
            for p in sorted(Path(pdb_dir).glob("*.pdb"))
            + sorted(Path(pdb_dir).glob("*.cif"))
        ]
        if not paths:
            print("No PDB/CIF in", pdb_dir, file=sys.stderr)
            sys.exit(1)
        print("Computing distances with TM-align...")
        ids, matrix = compute_rmsd_matrix_from_pdbs(paths)
    else:
        if not args.input or not os.path.isfile(args.input):
            ap.error("Input RMSD file required (or use --from-pdb DIR).")
        fmt = (
            args.format
            if args.format != "auto"
            else detect_format(args.input)
        )
        if fmt == "lsq_txt":
            alignments = parse_lsq_rmsd_txt(args.input)
        elif fmt == "ssm_txt":
            alignments = parse_ssm_rmsd_txt(args.input)
        else:
            alignments = parse_rmsd_csv(args.input)
        if not alignments:
            print("No RMSD pairs in", args.input, file=sys.stderr)
            sys.exit(1)
        ids, matrix = alignments_to_matrix(alignments)
        print("Loaded", len(ids), "labels,", len(alignments) // 2, "pairs.")

    # Optional: DALI-like per-query pseudo-z workflow (RMSD -> similarity -> z -> distance)
    try:
        matrix_for_tree, ranking = zscore_distance_matrix(
            ids,
            matrix,
            z_mode=args.pseudo_z,
            similarity=args.similarity,
            tau=args.tau,
            z_transform=args.zdist,
            z_scale=args.zscale,
        )
    except ValueError as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(2)
    if args.ranking_csv:
        write_ranking_csv(ranking, args.ranking_csv)

    use_unrooted = args.unrooted or args.no_midpoint_root
    if use_unrooted and args.root:
        print("Warning: --root is ignored when using --unrooted or --no-midpoint-root.", file=sys.stderr)
    newick = build_nj_tree(
        ids,
        matrix_for_tree,
        root=None if use_unrooted else args.root,
        midpoint_root=not use_unrooted,
    )
    with open(args.output, "w") as f:
        f.write(newick)
    print("Wrote", args.output)
    if args.plot:
        plot_tree(newick, args.plot, title="Structure-based phylogeny (RMSD)")


if __name__ == "__main__":
    main()
