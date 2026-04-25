#!/usr/bin/env python3
"""
Analyse a Caver 3.0 PyMOL plugin tunnel output folder.

Inputs expected (relative to caver output dir):
- analysis/tunnel_profiles.csv
- analysis/residues.txt
- analysis/bottlenecks.csv (optional; used for reporting bottleneck radius/xyz)

Computes:
- Tunnel electrostatic complementarity (EC) using a tunnel analogue of the McCoy et al. facing-point idea.
- Optional rolling-window "local" values along the tunnel for colouring (applies to EC and hydropathy).
- Hydropathy (Kyte–Doolittle) for residues lining each tunnel point.

Examples (from repository root):
  # Per-point CSV + PNG plot, coloured by EC.
  python metrics/caver_tunnel_analysis.py /path/to/caver_output_dir \
    --tunnel 26 --protein-pdb /path/to/protein.pdb --shell-a 3.0 --local-window 11 \
    --color-by ec --output-csv tunnel_26_points.csv --plot-out tunnel_26_ec.png

  # Per-point CSV + PNG plot, coloured by hydropathy
  # (orange = hydrophobic, white = neutral, green = hydrophilic).
  python metrics/caver_tunnel_analysis.py /path/to/caver_output_dir \
    --tunnel 26 --protein-pdb /path/to/protein.pdb --shell-a 3.0 \
    --color-by hydropathy --output-csv tunnel_26_points.csv --plot-out tunnel_26_hydro.png

  # Several Caver run directories: default CSV/PNG under each directory (tunnel_<N>_points.csv, etc.).
  python metrics/caver_tunnel_analysis.py run_a/ run_b/ run_c/ \
    --tunnel 26 --protein-pdb /path/to/protein.pdb --shell-a 3.0

  # Same, but one shared protein per run and explicit paths get a per-run suffix (stem_<run_basename>.ext).
  python metrics/caver_tunnel_analysis.py run_a/ run_b/ \
    --tunnel 26 --protein-pdb /path/to/protein.pdb -o /path/to/out/points.csv --plot-out /path/to/out/fig.png

  # Finer plot mesh (smoother fill): higher --plot-upsample (default 6), optional --plot-upsample-max 800.
  # Report/plot local opening as diameter: add --diameter (2× Caver radius in CSV/plot/JSON).
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import os
import re
import sys
from dataclasses import dataclass
from typing import Any, Iterable

try:
    from Bio.PDB import PDBParser  # type: ignore
    from Bio.PDB.NeighborSearch import NeighborSearch  # type: ignore
    from Bio.PDB.PDBExceptions import PDBConstructionWarning  # type: ignore
    import warnings

    warnings.simplefilter("ignore", PDBConstructionWarning)
    _BIOPYTHON_OK = True
except Exception:
    PDBParser = None  # type: ignore
    NeighborSearch = None  # type: ignore
    _BIOPYTHON_OK = False


# Per-residue EC sign weights for Caver residues.txt provenance (not the Kyte–Doolittle table).
_EC_RESIDUE_VALUE = {
    "ASP": -1.0,
    "GLU": -1.0,
    "LYS": 1.0,
    "ARG": 1.0,
}

_KYTE_DOOLITTLE = {
    "ILE": 4.5,
    "VAL": 4.2,
    "LEU": 3.8,
    "PHE": 2.8,
    "CYS": 2.5,
    "MET": 1.9,
    "ALA": 1.8,
    "GLY": -0.4,
    "THR": -0.7,
    "SER": -0.8,
    "TRP": -0.9,
    "TYR": -1.3,
    "PRO": -1.6,
    "HIS": -3.2,
    "GLU": -3.5,
    "GLN": -3.5,
    "ASP": -3.5,
    "ASN": -3.5,
    "LYS": -3.9,
    "ARG": -4.5,
}


def kyte_doolittle_hydropathy(resname3: str) -> float | None:
    r = (resname3 or "").upper().strip()
    return _KYTE_DOOLITTLE.get(r)

def _heavy_atoms_biopython(residue) -> list:
    out = []
    for a in residue.get_atoms():
        el = (getattr(a, "element", "") or "").strip().upper()
        if el == "H":
            continue
        out.append(a)
    return out


def _residue_sign_ec(resname3: str, *, his_is_positive: bool = True) -> int:
    """
    FoldKit EC-style sign:
    - ASP/GLU = -1
    - ARG/LYS/HIS = +1 (if his_is_positive)
    """
    r = (resname3 or "").upper().strip()
    if r in ("ASP", "GLU"):
        return -1
    if r in ("ARG", "LYS"):
        return 1
    if r == "HIS":
        return 1 if his_is_positive else 0
    return 0


def coulomb_potential(
    points_xyz: list[tuple[float, float, float]] | list[list[float]],
    ec_coords: list[list[float]],
    ec_values: list[float],
    *,
    epsilon_r: float = 80.0,
    r_min: float = 1.0,
) -> list[float]:
    """
    Lightweight electrostatic potential estimate, consistent with the EC module:
    phi(p) = Σ_i w_i / (epsilon_r * max(|r_i - p|, r_min))

    Units are arbitrary (relative scale only) because EC weights are heuristic.
    """
    import numpy as np

    pts = np.asarray(points_xyz, dtype=float)
    coords = np.asarray(ec_coords, dtype=float)
    w = np.asarray(ec_values, dtype=float)
    if pts.size == 0:
        return []
    if coords.size == 0 or w.size == 0:
        return [0.0 for _ in range(int(pts.shape[0]))]

    # Broadcast distances: [P, A, 3] -> [P, A]
    d = np.linalg.norm(pts[:, None, :] - coords[None, :, :], axis=2)
    d = np.maximum(d, float(r_min))
    phi = (w[None, :] / (float(epsilon_r) * d)).sum(axis=1)
    return [float(x) for x in phi]


def pearson_r(x: list[float], y: list[float]) -> float | None:
    """
    Pearson correlation (like electrostatic_complementarity._pearson_r):
    returns 0.0 when variance is zero, and None when lengths mismatch / empty.
    """
    import numpy as np

    xx = np.asarray(x, dtype=float)
    yy = np.asarray(y, dtype=float)
    if xx.size == 0 or yy.size == 0 or xx.size != yy.size:
        return None
    x0 = xx - float(np.mean(xx))
    y0 = yy - float(np.mean(yy))
    den = float(np.sqrt(float(np.sum(x0 * x0) * np.sum(y0 * y0))))
    if not math.isfinite(den) or den <= 0.0:
        return 0.0
    r = float(np.sum(x0 * y0) / den)
    return max(-1.0, min(1.0, r))


def _ec_plot_cmap_rwb():
    """
    Diverging map: red (low / negative) -> white (0) -> blue (high / positive),
    as used for EC r in many viewers.
    """
    from matplotlib.colors import LinearSegmentedColormap  # type: ignore

    return LinearSegmentedColormap.from_list("ec_rwb", ["#b2182b", "#f7f7f7", "#2166ac"])


def _ec_diverging_norm_from_values(color_values: list[float]) -> tuple[Any, float, float]:
    """
    TwoSlopeNorm centred at 0. Uses finite data min/max when they straddle 0; otherwise
    a symmetric range around 0 so white still marks zero (e.g. all-positive r).
    """
    import numpy as np
    from matplotlib.colors import TwoSlopeNorm  # type: ignore

    a = np.asarray(color_values, dtype=float)
    a = a[np.isfinite(a)]
    if a.size == 0:
        lo, hi = -1.0, 1.0
    else:
        lo, hi = float(np.min(a)), float(np.max(a))
        if lo == hi:
            s = max(abs(lo), 1e-6)
            lo, hi = -s, s
        elif lo < 0.0 < hi:
            pass
        else:
            s = max(abs(lo), abs(hi), 1e-6)
            lo, hi = -s, s
    if not (lo < 0.0 < hi):
        s = max(abs(lo), abs(hi), 1e-6)
        lo, hi = -s, s
    norm = TwoSlopeNorm(vmin=lo, vcenter=0.0, vmax=hi)
    return norm, lo, hi


def _diverging_norm_fixed(vmin: float, vmax: float) -> tuple[Any, float, float]:
    """TwoSlopeNorm centred at 0 with explicit limits."""
    from matplotlib.colors import TwoSlopeNorm  # type: ignore

    lo, hi = float(vmin), float(vmax)
    if not (lo < 0.0 < hi):
        s = max(abs(lo), abs(hi), 1e-6)
        lo, hi = -s, s
    return TwoSlopeNorm(vmin=lo, vcenter=0.0, vmax=hi), lo, hi


def _diverging_norm_from_series(color_values: list[float], *, fallback: tuple[float, float]) -> tuple[Any, float, float]:
    """
    Diverging scaling centred at 0.

    Uses finite data limits when they straddle 0; otherwise expands to a symmetric range about 0
    so the neutral point remains white. Falls back to the provided limits when there is no data.
    """
    import numpy as np

    a = np.asarray(color_values, dtype=float)
    a = a[np.isfinite(a)]
    if a.size == 0:
        return _diverging_norm_fixed(fallback[0], fallback[1])
    return _ec_diverging_norm_from_values([float(x) for x in a.tolist()])


def _rolling_mean_ignore_nan(values: list[float], *, window: int) -> list[float]:
    """
    Rolling mean with an odd window, ignoring non-finite values.

    Intended for plot-only smoothing of per-point colour series.
    """
    win = int(window or 0)
    if win < 3:
        return list(values)
    if win % 2 == 0:
        win += 1
    half = win // 2
    out: list[float] = []
    n = len(values)
    for i in range(n):
        lo = max(0, i - half)
        hi = min(n, i + half + 1)
        seg = values[lo:hi]
        seg_ok = [float(x) for x in seg if isinstance(x, (int, float)) and math.isfinite(float(x))]
        out.append(float(sum(seg_ok) / len(seg_ok)) if seg_ok else float("nan"))
    return out


def _compress_unique_d(
    d: Any, *ys: Any
) -> tuple[Any, ...]:
    """Sort by d and average r, v (etc.) at duplicate distance values for stable interpolation."""
    import numpy as np

    d = np.asarray(d, dtype=float).ravel()
    ys = [np.asarray(y, dtype=float).ravel() for y in ys]
    for y in ys:
        if y.shape[0] != d.shape[0]:
            raise ValueError("length mismatch in _compress_unique_d")
    o = np.argsort(d)
    d = d[o]
    ys = [y[o] for y in ys]
    if d.size == 0:
        return (d,) + tuple(ys)
    du: list[float] = []
    out: list[list[float]] = [[] for _ in ys]
    i = 0
    n = d.size
    while i < n:
        j = i
        v = float(d[i])
        while j < n and float(d[j]) == v:
            j += 1
        du.append(v)
        for k, y in enumerate(ys):
            out[k].append(float(np.nanmean(y[i:j])))
        i = j
    return (np.array(du, dtype=float),) + tuple(np.array(t, dtype=float) for t in out)


def _upsample_along_distance(
    d: Any,
    r: Any,
    vals: Any,
    *,
    factor: float,
    max_points: int,
) -> tuple[Any, Any, Any]:
    """
    Insert additional distance samples so the fill uses smaller quads and appears smoother.
    """
    import numpy as np

    d = np.asarray(d, dtype=float).ravel()
    r = np.asarray(r, dtype=float).ravel()
    vals = np.asarray(vals, dtype=float).ravel()
    n = d.size
    if n < 2 or float(factor) <= 1.0:
        return d, r, vals
    du, r_u, v_u = _compress_unique_d(d, r, vals)
    if du.size < 2:
        return d, r, vals
    n_out = int(min(int(max_points), max(du.size, int((du.size - 1) * float(factor)) + 1)))
    if n_out <= du.size:
        return d, r, vals
    d_lo, d_hi = float(du[0]), float(du[-1])
    if not (d_hi > d_lo):
        return d, r, vals
    d_new = np.linspace(d_lo, d_hi, n_out)
    v_fill = v_u.copy()
    x_idx = np.arange(du.size, dtype=float)
    m = np.isfinite(v_fill)
    if not m.all():
        if m.sum() == 0:
            pass
        elif m.sum() == 1:
            v_fill[:] = v_fill[m][0]
        else:
            v_fill[~m] = np.interp(
                x_idx[~m],
                x_idx[m],
                v_fill[m],
            )
    v_new = np.interp(d_new, du, v_fill)
    try:
        from scipy.interpolate import PchipInterpolator  # type: ignore

        pr = PchipInterpolator(du, r_u, extrapolate=True)
        r_new = np.clip(np.asarray(pr(d_new), dtype=float), 0.0, None)
    except Exception:
        r_new = np.maximum(0.0, np.interp(d_new, du, r_u))
    return d_new, r_new, v_new


@dataclass(frozen=True)
class LiningResidue:
    chain_id: str
    resseq: int
    resname: str  # 3-letter
    n_snapshots: int | None = None
    n_side_snapshots: int | None = None


@dataclass(frozen=True)
class TunnelProfile:
    snapshot: str
    tunnel_cluster: int
    tunnel: int
    length: float | None
    bottleneck_radius: float | None
    distance: list[float]
    x: list[float]
    y: list[float]
    z: list[float]
    r: list[float]


def _read_text(path: str) -> str:
    with open(path, "r", encoding="utf-8", errors="replace") as f:
        return f.read()


def _find_analysis_path(base_dir: str, rel: str) -> str:
    p = os.path.join(base_dir, rel)
    if not os.path.isfile(p):
        raise FileNotFoundError(p)
    return p


def parse_residues_txt(path: str, tunnel_cluster: int) -> list[LiningResidue]:
    """
    Parse Caver `analysis/residues.txt` and return lining residues for one tunnel cluster.

    Expected section format:
      == Tunnel cluster 26 ==
      # C    res  AA       N   sideN atoms
        B     43 GLU       1       1 ...
    """
    lines = _read_text(path).splitlines()
    header = f"== Tunnel cluster {int(tunnel_cluster)} =="
    start = None
    for i, line in enumerate(lines):
        if line.strip() == header:
            start = i + 1
            break
    if start is None:
        raise ValueError(f"Tunnel cluster {tunnel_cluster} not found in residues.txt")

    out: list[LiningResidue] = []
    i = start
    while i < len(lines):
        line = lines[i].rstrip("\n")
        s = line.strip()
        if not s:
            i += 1
            continue
        if s.startswith("== Tunnel cluster "):
            break
        if s.startswith("#") or s.startswith("==="):
            i += 1
            continue

        parts = s.split()
        # Minimal expected columns: chain, res, AA, N, sideN
        if len(parts) < 5:
            i += 1
            continue
        chain_id = parts[0]
        try:
            resseq = int(parts[1])
        except Exception:
            i += 1
            continue
        resname = parts[2].upper()
        try:
            n = int(parts[3])
        except Exception:
            n = None
        try:
            side_n = int(parts[4])
        except Exception:
            side_n = None
        out.append(LiningResidue(chain_id=chain_id, resseq=resseq, resname=resname, n_snapshots=n, n_side_snapshots=side_n))
        i += 1

    return out


def _float_or_none(x: str) -> float | None:
    s = (x or "").strip()
    if not s or s == "-" or s.lower() == "nan":
        return None
    try:
        v = float(s)
    except Exception:
        return None
    if not math.isfinite(v):
        return None
    return v


def parse_tunnel_profiles_csv(path: str, tunnel_cluster: int, *, tunnel: int | None = None, snapshot: str | None = None) -> TunnelProfile:
    """
    Parse `analysis/tunnel_profiles.csv` for one tunnel cluster (+ optional tunnel id + optional snapshot).

    Caver writes multiple rows per tunnel: one for each axis label (X/Y/Z/distance/length/R).
    """
    tunnel_cluster = int(tunnel_cluster)
    want_tunnel = tunnel_cluster if tunnel is None else int(tunnel)

    rows = []
    with open(path, "r", encoding="utf-8", errors="replace", newline="") as f:
        r = csv.reader(f, skipinitialspace=True)
        header = next(r, None)
        if not header:
            raise ValueError("Empty tunnel_profiles.csv")
        for rec in r:
            if len(rec) < 13:
                continue
            snap = rec[0].strip()
            tc = int(rec[1])
            t = int(rec[2])
            if tc != tunnel_cluster or t != want_tunnel:
                continue
            if snapshot and snap != snapshot:
                continue
            axis = (rec[12] or "").strip()
            vals = [v.strip() for v in rec[13:] if v.strip() != ""]
            rows.append((snap, tc, t, axis, rec, vals))

    if not rows:
        raise ValueError(f"No rows found for tunnel_cluster={tunnel_cluster} tunnel={want_tunnel} in tunnel_profiles.csv")

    # Use the first matching snapshot unless user constrained it.
    snap0 = rows[0][0]
    rows = [x for x in rows if x[0] == snap0]

    # Common scalar metadata live in the first few columns.
    rec0 = rows[0][4]
    length = _float_or_none(rec0[10])
    bottleneck_radius = _float_or_none(rec0[5])

    by_axis: dict[str, list[float]] = {}
    for (_snap, _tc, _t, axis, _rec, vals) in rows:
        axis_norm = axis.strip()
        if not axis_norm:
            continue
        fvals = []
        for v in vals:
            fv = _float_or_none(v)
            if fv is None:
                continue
            fvals.append(float(fv))
        by_axis[axis_norm] = fvals

    required = ["distance", "X", "Y", "Z", "R"]
    missing = [k for k in required if k not in by_axis or not by_axis[k]]
    if missing:
        raise ValueError(f"Missing required axis series in tunnel_profiles.csv for tunnel {tunnel_cluster}: {missing}")

    n = min(len(by_axis["distance"]), len(by_axis["X"]), len(by_axis["Y"]), len(by_axis["Z"]), len(by_axis["R"]))
    return TunnelProfile(
        snapshot=snap0,
        tunnel_cluster=tunnel_cluster,
        tunnel=want_tunnel,
        length=length,
        bottleneck_radius=bottleneck_radius,
        distance=by_axis["distance"][:n],
        x=by_axis["X"][:n],
        y=by_axis["Y"][:n],
        z=by_axis["Z"][:n],
        r=by_axis["R"][:n],
    )


def parse_bottlenecks_csv(path: str, tunnel_cluster: int, *, tunnel: int | None = None, snapshot: str | None = None) -> dict[str, Any] | None:
    tunnel_cluster = int(tunnel_cluster)
    want_tunnel = tunnel_cluster if tunnel is None else int(tunnel)
    with open(path, "r", encoding="utf-8", errors="replace", newline="") as f:
        r = csv.reader(f, skipinitialspace=True)
        header = next(r, None)
        if not header:
            return None
        # Snapshot,Tunnel cluster,Tunnel,Throughput,Cost,Bottleneck X,Bottleneck Y,Bottleneck Z,Bottleneck R,Bottleneck residues
        for rec in r:
            if len(rec) < 10:
                continue
            snap = rec[0].strip()
            tc = int(rec[1])
            t = int(rec[2])
            if tc != tunnel_cluster or t != want_tunnel:
                continue
            if snapshot and snap != snapshot:
                continue
            return {
                "snapshot": snap,
                "tunnel_cluster": tc,
                "tunnel": t,
                "throughput": _float_or_none(rec[3]),
                "cost": _float_or_none(rec[4]),
                "bottleneck_x": _float_or_none(rec[5]),
                "bottleneck_y": _float_or_none(rec[6]),
                "bottleneck_z": _float_or_none(rec[7]),
                "bottleneck_r": _float_or_none(rec[8]),
                "bottleneck_residues_raw": (rec[9] or "").strip(),
            }
    return None


def _parse_centerline_pdb(path: str) -> tuple[list[float], list[float], list[float], list[float]]:
    """
    Parse a Caver tunnel PDB like `*_026.pdb` which contains pseudo-atoms for the
    centreline and stores the local radius in the last numeric column.

    Returns x,y,z,r arrays in file order.
    """
    xs: list[float] = []
    ys: list[float] = []
    zs: list[float] = []
    rs: list[float] = []
    for line in _read_text(path).splitlines():
        if not line.startswith("ATOM"):
            continue
        # Caver tunnel PDBs are not strict PDB format; rely on whitespace splitting.
        parts = line.split()
        if len(parts) < 9:
            continue
        try:
            x = float(parts[5])
            y = float(parts[6])
            z = float(parts[7])
            r = float(parts[8])
        except Exception:
            continue
        xs.append(x)
        ys.append(y)
        zs.append(z)
        rs.append(r)
    if not xs:
        raise ValueError(f"No ATOM records parsed from centreline PDB: {path}")
    return xs, ys, zs, rs


def _dist_along_polyline(xs: list[float], ys: list[float], zs: list[float]) -> list[float]:
    d: list[float] = [0.0]
    for i in range(1, len(xs)):
        dx = xs[i] - xs[i - 1]
        dy = ys[i] - ys[i - 1]
        dz = zs[i] - zs[i - 1]
        d.append(d[-1] + float(math.sqrt(dx * dx + dy * dy + dz * dz)))
    return d


def residue_ec(resname3: str, *, his_ec: float = 0.1) -> float:
    """EC-style signed value for one residue (HIS uses `his_ec` as partial +1)."""
    r = (resname3 or "").upper().strip()
    if r == "HIS":
        return float(his_ec)
    return float(_EC_RESIDUE_VALUE.get(r, 0.0))


def summarise_lining_residues(residues: list[LiningResidue], *, his_ec: float = 0.1) -> dict[str, Any]:
    """
    Summary of unique lining residues from Caver (EC sign model + Kyte–Doolittle where defined).
    """
    uniq = {(r.chain_id, r.resseq, r.resname) for r in residues}
    uniq_list = [LiningResidue(chain_id=c, resseq=i, resname=aa) for (c, i, aa) in sorted(uniq)]
    n = len(uniq_list)
    ec_values = [residue_ec(r.resname, his_ec=his_ec) for r in uniq_list]
    overall = float(sum(ec_values))
    hydros = [kyte_doolittle_hydropathy(r.resname) for r in uniq_list]
    n_def = sum(1 for h in hydros if h is not None)
    n_miss = n - n_def
    mean_h: float | None
    if n_def:
        mean_h = float(sum(h for h in hydros if h is not None) / n_def)
    else:
        mean_h = None
    return {
        "n_lining_residues": n,
        "overall_ec": overall,
        "mean_ec_per_residue": (overall / n) if n else 0.0,
        "mean_hydropathy_kyte_doolittle": mean_h,
        "n_hydropathy_defined": n_def,
        "n_hydropathy_missing": n_miss,
        "lining_residues": [{"chain": r.chain_id, "res": r.resseq, "aa": r.resname} for r in uniq_list],
    }


def _load_protein_neighbor_search(pdb_path: str):
    if not _BIOPYTHON_OK or PDBParser is None or NeighborSearch is None:
        raise SystemExit(
            "BioPython is required for protein mapping. Install it (e.g. `pip install biopython`) "
            "or run without --protein-pdb."
        )
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("protein", pdb_path)
    atoms = [a for a in structure.get_atoms()]
    return structure, NeighborSearch(atoms)


def _residue_id(res) -> tuple[str, int, str]:
    """
    Return a stable (chain, resseq, resname) identifier for a Biopython residue.
    """
    try:
        chain_id = res.get_parent().id
    except Exception:
        chain_id = ""
    hetflag, resseq, icode = getattr(res, "id", (" ", None, " "))
    resname = getattr(res, "resname", "").strip().upper()
    return str(chain_id), int(resseq) if resseq is not None else -1, f"{resname}{str(icode).strip() or ''}"


def map_profile_to_local_properties(
    prof: TunnelProfile,
    protein_pdb: str,
    *,
    shell_a: float = 3.0,
    ec_epsilon_r: float = 4.0,
    ec_r_min: float = 1.0,
    local_window: int = 0,
) -> dict[str, Any]:
    """
    For each centreline node, collect protein residues with any atom within (radius + shell_a)
    and compute local EC wall potentials and mean Kyte–Doolittle hydropathy.
    """
    structure, ns = _load_protein_neighbor_search(protein_pdb)

    local_phi_wall_pos: list[float] = []
    local_phi_wall_neg: list[float] = []
    local_ec_r_window: list[float] = []
    local_hydro: list[float] = []
    local_hydro_window: list[float] = []
    local_n_res: list[int] = []
    local_res_ids: list[list[tuple[str, int, str]]] = []

    union: dict[tuple[str, int, str], dict[str, Any]] = {}

    # Per-residue heavy-atom sites + per-atom EC weights (signed, spread over heavy atoms).
    # Keyed by the stable residue id produced by `_residue_id`.
    residue_ec_atoms: dict[tuple[str, int, str], tuple[list[list[float]], list[float]]] = {}
    for res in structure.get_residues():
        rid = _residue_id(res)
        aa3 = rid[2][:3]
        if len(aa3) != 3 or not aa3.isalpha():
            continue
        sign = _residue_sign_ec(aa3, his_is_positive=True)
        if sign == 0:
            continue
        atoms = _heavy_atoms_biopython(res)
        if not atoms:
            continue
        per_atom = float(sign) / float(len(atoms))
        coords = [[float(x) for x in a.coord] for a in atoms]
        ec_w = [per_atom for _ in atoms]
        residue_ec_atoms[rid] = (coords, ec_w)

    # Build a stable "wall direction" per node from the centreline tangent.
    # Facing points are taken as p ± r*u, where u is a unit vector perpendicular to the tangent.
    # This is an approximation of "facing points across a tunnel".
    import numpy as np

    pts = np.column_stack([np.asarray(prof.x, dtype=float), np.asarray(prof.y, dtype=float), np.asarray(prof.z, dtype=float)])
    rr = np.asarray(prof.r, dtype=float)
    if pts.shape[0] >= 2:
        tang = np.zeros_like(pts)
        tang[1:-1] = pts[2:] - pts[:-2]
        tang[0] = pts[1] - pts[0]
        tang[-1] = pts[-1] - pts[-2]
        # Normalise tangents.
        tnorm = np.linalg.norm(tang, axis=1)
        tnorm = np.maximum(tnorm, 1e-9)
        tang = tang / tnorm[:, None]
    else:
        tang = np.array([[0.0, 0.0, 1.0]], dtype=float)

    ref = np.array([0.0, 0.0, 1.0], dtype=float)
    alt = np.array([0.0, 1.0, 0.0], dtype=float)

    for idx, (x, y, z, r0) in enumerate(zip(prof.x, prof.y, prof.z, prof.r, strict=False)):
        radius = float(max(0.0, r0)) + float(max(0.0, shell_a))
        # NeighborSearch expects a list-like coordinate.
        near_atoms = ns.search([float(x), float(y), float(z)], radius, level="A")
        residues = {}
        for a in near_atoms:
            res = a.get_parent()
            if res is None:
                continue
            rid = _residue_id(res)
            # Filter out non-standard residues/water by requiring a 3-letter AA name.
            aa3 = rid[2][:3]
            if len(aa3) != 3 or not aa3.isalpha():
                continue
            residues[rid] = aa3

        ids = sorted(residues.keys())
        local_res_ids.append(ids)
        local_n_res.append(len(ids))

        # Coulomb-style potential at this point from EC weights on the selected residues.
        phi_ec_coords: list[list[float]] = []
        phi_ec_values: list[float] = []
        hydros: list[float] = []
        for rid in ids:
            aa3 = residues[rid]
            h = kyte_doolittle_hydropathy(aa3)
            if h is not None:
                hydros.append(float(h))
            atom_pack = residue_ec_atoms.get(rid)
            if atom_pack is not None:
                c0, w0 = atom_pack
                phi_ec_coords.extend(c0)
                phi_ec_values.extend(w0)

            if rid not in union:
                union[rid] = {
                    "chain": rid[0],
                    "res": rid[1],
                    "aa": aa3,
                }

        # EC-style facing-wall potentials.
        t = tang[idx] if idx < tang.shape[0] else np.array([0.0, 0.0, 1.0], dtype=float)
        # Choose a reference vector not parallel to t.
        ref_use = ref if abs(float(np.dot(t, ref))) < 0.95 else alt
        u = np.cross(t, ref_use)
        un = float(np.linalg.norm(u))
        if not math.isfinite(un) or un <= 1e-9:
            u = np.array([1.0, 0.0, 0.0], dtype=float)
            un = 1.0
        u = u / un
        r_wall = float(max(0.0, r0))
        p = np.array([float(x), float(y), float(z)], dtype=float)
        p_pos = p + u * r_wall
        p_neg = p - u * r_wall

        phi_pos = coulomb_potential(
            [p_pos.tolist()], phi_ec_coords, phi_ec_values, epsilon_r=float(ec_epsilon_r), r_min=float(ec_r_min)
        )[0]
        phi_neg = coulomb_potential(
            [p_neg.tolist()], phi_ec_coords, phi_ec_values, epsilon_r=float(ec_epsilon_r), r_min=float(ec_r_min)
        )[0]
        local_phi_wall_pos.append(float(phi_pos))
        local_phi_wall_neg.append(float(phi_neg))

        local_hydro.append(float(sum(hydros) / len(hydros)) if hydros else float("nan"))

    # Global tunnel-EC: correlation across nodes of phi_pos vs -phi_neg.
    tunnel_ec_r = pearson_r(local_phi_wall_pos, [-x for x in local_phi_wall_neg])
    n_pairs = int(len(local_phi_wall_pos))

    # Optional rolling window along the tunnel (used for colouring).
    # Historically this was called `ec_window` because it only applied to EC.
    win = int(local_window or 0)
    if win and win >= 3:
        if win % 2 == 0:
            win += 1
        half = win // 2
        for i in range(n_pairs):
            lo = max(0, i - half)
            hi = min(n_pairs, i + half + 1)
            rr0 = pearson_r(local_phi_wall_pos[lo:hi], [-x for x in local_phi_wall_neg[lo:hi]])
            local_ec_r_window.append(float(rr0) if rr0 is not None else float("nan"))
    else:
        local_ec_r_window = [float("nan") for _ in range(n_pairs)]

    # Rolling window for hydropathy (mean of the per-node lining-residue mean).
    if win and win >= 3:
        if win % 2 == 0:
            win += 1
        half = win // 2
        for i in range(n_pairs):
            lo = max(0, i - half)
            hi = min(n_pairs, i + half + 1)
            seg = local_hydro[lo:hi]
            seg_ok = [float(x) for x in seg if x is not None and isinstance(x, (int, float)) and math.isfinite(float(x))]
            local_hydro_window.append(float(sum(seg_ok) / len(seg_ok)) if seg_ok else float("nan"))
    else:
        local_hydro_window = list(local_hydro)

    union_list = [union[k] for k in sorted(union.keys())]

    return {
        "protein_pdb": os.path.abspath(protein_pdb),
        "shell_a": float(shell_a),
        "overall_from_mapped_residues": {"n_residues": int(len(union_list))},
        "tunnel_ec": {
            "r": tunnel_ec_r,
            "n_pairs": n_pairs,
            "epsilon_r": float(ec_epsilon_r),
            "r_min": float(ec_r_min),
            "window": int(win) if win and win >= 3 else 0,
        },
        "local": {
            "phi_wall_pos": local_phi_wall_pos,
            "phi_wall_neg": local_phi_wall_neg,
            "ec_r_window": local_ec_r_window,
            "hydropathy_mean": local_hydro,
            "hydropathy_mean_window": local_hydro_window,
            "n_residues": local_n_res,
            "residue_ids": local_res_ids,
        },
        "mapped_residues": union_list,
    }


def _maybe_import_matplotlib():
    try:
        import matplotlib  # noqa: F401
        import matplotlib.pyplot as plt  # noqa: F401
        from matplotlib.collections import LineCollection  # noqa: F401
    except Exception as e:
        raise SystemExit(
            "matplotlib is required for plotting. Install it (e.g. `pip install matplotlib`) "
            f"or rerun with --no-plot. Import error: {e}"
        )


def _axes_lower_edge_figure(fig, ax) -> float:
    """Lowest y in figure coordinates (0–1) of the axes tight bbox (includes tick labels)."""
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


def _add_raster_colorbar_tunnel_horizontal(
    fig,
    ax,
    mappable,
    *,
    label: str,
) -> object:
    """
    Horizontal colour bar below the main axes, same width as the tunnel plot.
    Gradient is rasterised; ticks and label stay vector.
    """
    ch = max(0.007, min(0.032, 0.022))
    hpad = 0.018
    bottom = 0.0
    for _ in range(10):
        fig.canvas.draw()
        lo = _axes_lower_edge_figure(fig, ax)
        bottom = lo - hpad - ch
        if bottom >= 0.02:
            break
        cur = float(fig.subplotpars.bottom)
        fig.subplots_adjust(bottom=min(0.45, cur + 0.05))
    if bottom < 0.01:
        bottom = 0.01
    fig.canvas.draw()
    pos = ax.get_position()
    xmin, bw = float(pos.x0), float(pos.width)
    cax = fig.add_axes([xmin, bottom, bw, ch])
    cbar = fig.colorbar(
        mappable,
        cax=cax,
        orientation="horizontal",
        extend="neither",
        label=label,
    )
    if hasattr(cbar, "solids") and cbar.solids is not None:
        cbar.solids.set_rasterized(True)
    for p in cbar.ax.patches:
        if hasattr(p, "set_rasterized"):
            p.set_rasterized(True)
    return cbar


def _savefig_tunnel_vector(
    fig,
    out_path: str,
    *,
    dpi: int,
    rasterise_artists: list[object] | None = None,
) -> None:
    """
    SVG/PDF: keep text/axes editable; optionally rasterise selected artists (e.g. the tunnel fill),
    and always rasterise the colour bar gradient.
    """
    import os

    import matplotlib as mpl

    ext = os.path.splitext(out_path)[1].lower()
    if fig.axes:
        main_ax = fig.axes[0]
        for artist in main_ax.get_children():
            if hasattr(artist, "set_rasterized"):
                try:
                    artist.set_rasterized(False)
                except (AttributeError, TypeError):
                    pass
    if rasterise_artists:
        for a in rasterise_artists:
            if hasattr(a, "set_rasterized"):
                try:
                    a.set_rasterized(True)
                except (AttributeError, TypeError):
                    pass
    rc_extra: dict = {}
    if ext == ".svg":
        rc_extra["svg.fonttype"] = "none"
    fmt = ext[1:] if ext.startswith(".") else ext
    with mpl.rc_context(rc_extra):
        fig.savefig(out_path, format=fmt, dpi=dpi, bbox_inches="tight")


def plot_profile(
    prof: TunnelProfile,
    *,
    out_path: str,
    cbar_label: str,
    color_values: list[float],
    title: str,
    cmap: Any,
    vmin: float | None = None,
    vmax: float | None = None,
    norm: Any | None = None,
    annotate: dict[str, Any] | None = None,
    upsample_factor: float = 6.0,
    upsample_max_points: int = 800,
    as_diameter: bool = False,
    rasterise_fill: bool = False,
    rasterise_fill_dpi: int = 300,
    xlim: tuple[float, float] | None = None,
    invert_y: bool = True,
) -> None:
    _maybe_import_matplotlib()
    import matplotlib.pyplot as plt
    from matplotlib.collections import LineCollection
    from matplotlib.collections import PolyCollection
    import numpy as np

    # Vertical projection: distance is the vertical axis; horizontal is Caver radius, or 2× for diameter.
    d_anno = np.asarray(prof.distance, dtype=float)
    r_caver = np.asarray(prof.r, dtype=float)
    wfac = 2.0 if as_diameter else 1.0

    vals = np.asarray(color_values, dtype=float)
    if vals.shape[0] != d_anno.shape[0]:
        raise ValueError("color_values length must match profile length")
    d, r, vals = _upsample_along_distance(
        d_anno, r_caver, vals, factor=float(upsample_factor), max_points=int(upsample_max_points)
    )
    vals_mid = (vals[:-1] + vals[1:]) / 2.0

    def _segments(x: np.ndarray, y: np.ndarray) -> np.ndarray:
        pts = np.column_stack([x, y])
        return np.stack([pts[:-1], pts[1:]], axis=1)

    x_rt = wfac * r
    x_lf = -wfac * r
    seg_rt = _segments(x_rt, d)
    seg_lf = _segments(x_lf, d)

    # Portrait figure: tunnel length (y) usually dominates; extra height gives the profile more
    # vertical room. Width matches a single-column size; fixed margins avoid tight_layout squashing.
    fig, ax = plt.subplots(figsize=(7.0, 20.0), dpi=160)
    # Filled tunnel interior, coloured by per-segment value.
    # Build quads between consecutive points (x,y):
    #   (+r[i], d[i]) -> (+r[i+1], d[i+1]) -> (-r[i+1], d[i+1]) -> (-r[i], d[i])
    quads = np.zeros((len(d) - 1, 4, 2), dtype=float)
    quads[:, 0, 0] = x_rt[:-1]
    quads[:, 0, 1] = d[:-1]
    quads[:, 1, 0] = x_rt[1:]
    quads[:, 1, 1] = d[1:]
    quads[:, 2, 0] = x_lf[1:]
    quads[:, 2, 1] = d[1:]
    quads[:, 3, 0] = x_lf[:-1]
    quads[:, 3, 1] = d[:-1]

    pc = PolyCollection(quads, cmap=cmap, edgecolors="none", alpha=0.85)
    pc.set_array(vals_mid)
    if norm is not None:
        pc.set_norm(norm)
    elif vmin is not None or vmax is not None:
        pc.set_clim(vmin=vmin, vmax=vmax)
    ax.add_collection(pc)

    # Monochrome outline for legibility.
    lc_rt = LineCollection(seg_rt, colors=["#333333"], linewidths=1.2)
    lc_lf = LineCollection(seg_lf, colors=["#333333"], linewidths=1.2)
    ax.add_collection(lc_rt)
    ax.add_collection(lc_lf)

    ax.axvline(0.0, color="#999999", linewidth=0.8)
    ax.set_xlabel("Diameter (Å)" if as_diameter else "Radius (Å)")
    ax.set_ylabel("Distance along tunnel (Å)")
    ax.set_title(title)
    r_span = float(wfac * max(np.max(r_caver), float(np.max(r))))
    if xlim is not None:
        ax.set_xlim(float(xlim[0]), float(xlim[1]))
    else:
        ax.set_xlim(-r_span * 1.12, r_span * 1.12)
    ax.set_ylim(float(np.min(d_anno)), float(np.max(d_anno)))
    if invert_y:
        ax.invert_yaxis()
    ax.set_aspect("auto")
    # Fixed margins: tight_layout() shrinks the main axes; we set margins once, then add margin labels.
    fig.subplots_adjust(left=0.1, right=0.9, top=0.9, bottom=0.1)

    if annotate and isinstance(annotate, dict):
        pts = annotate.get("points") or []
        if pts:
            from matplotlib.transforms import blended_transform_factory

            # x in axes (>1) to the right of the spines, y in data: text in the margin; larger x = further right.
            trans_text = blended_transform_factory(ax.transAxes, ax.transData)
            ax_text_x = 1.16
            for p in pts:
                if not isinstance(p, dict):
                    continue
                idx = int(p.get("idx", -1))
                if idx < 0 or idx >= len(d_anno):
                    continue
                lab = str(p.get("label", "") or "")
                residues = str(p.get("residues", "") or "")
                y0 = float(d_anno[idx])
                r0 = float(r_caver[idx])
                x0 = wfac * r0

                ax.plot([0.0], [y0], marker="o", markersize=3.0, color="#111111", zorder=6)
                ax.plot([0.0, x0], [y0, y0], color="#111111", linewidth=0.8, alpha=0.7, zorder=6)

                if as_diameter:
                    txt = f"{lab}  d={y0:.2f} Å  D={x0:.2f} Å"
                else:
                    txt = f"{lab}  d={y0:.2f} Å  r={r0:.2f} Å"
                if residues:
                    txt += f"\n{residues}"
                ax.annotate(
                    txt,
                    xy=(x0, y0),
                    xycoords="data",
                    xytext=(ax_text_x, y0),
                    textcoords=trans_text,
                    va="center",
                    ha="left",
                    fontsize=7.0,
                    clip_on=False,
                    zorder=7,
                    arrowprops={"arrowstyle": "-", "color": "#333333", "lw": 0.55},
                )

    os.makedirs(os.path.dirname(os.path.abspath(out_path)) or ".", exist_ok=True)
    _add_raster_colorbar_tunnel_horizontal(fig, ax, pc, label=cbar_label)
    # Colour bar may call subplots_adjust(bottom=...); keep left/right/top margins so the tunnel panel width is stable.
    b = float(fig.subplotpars.bottom)
    fig.subplots_adjust(left=0.1, right=0.9, top=0.9, bottom=b)
    out_abs = os.path.abspath(out_path)
    ext = os.path.splitext(out_abs)[1].lower()
    dpi = 160
    if ext in (".svg", ".pdf"):
        _savefig_tunnel_vector(
            fig,
            out_abs,
            dpi=int(rasterise_fill_dpi) if rasterise_fill else dpi,
            rasterise_artists=[pc] if rasterise_fill else None,
        )
    else:
        fig.savefig(out_abs, bbox_inches="tight", dpi=dpi)
    plt.close(fig)


def _write_points_csv(
    out_csv_path: str,
    *,
    prof: TunnelProfile,
    caver_output_dir: str,
    bottleneck: dict[str, Any] | None,
    lining_summary: dict[str, Any],
    mapped: dict[str, Any] | None,
    residue_list_max: int,
    as_diameter: bool = False,
) -> None:
    """
    Write one CSV row per profile point.
    Includes a residue list (semicolon-separated) when mapping is enabled.
    """
    os.makedirs(os.path.dirname(os.path.abspath(out_csv_path)) or ".", exist_ok=True)
    wfac = 2.0 if as_diameter else 1.0
    col_width = "diameter_A" if as_diameter else "radius_A"
    col_bneck_prof = "bottleneck_diameter_A_profiles" if as_diameter else "bottleneck_radius_A_profiles"
    col_bneck_csv = "bottleneck_diameter_A_bottlenecks_csv" if as_diameter else "bottleneck_radius_A_bottlenecks_csv"

    if mapped is not None:
        loc = mapped["local"]
        phi_wall_pos = list(loc.get("phi_wall_pos") or [])
        phi_wall_neg = list(loc.get("phi_wall_neg") or [])
        ec_r_window = list(loc.get("ec_r_window") or [])
        n_series = list(loc["n_residues"])
        ids_series = list(loc["residue_ids"])
        tunnel_ec = mapped.get("tunnel_ec") or {}
        tunnel_ec_r = tunnel_ec.get("r")
        tunnel_ec_n = tunnel_ec.get("n_pairs")
        protein_pdb = mapped.get("protein_pdb")
        shell_a = mapped.get("shell_a")
    else:
        phi_wall_pos = [float("nan")] * len(prof.distance)
        phi_wall_neg = [float("nan")] * len(prof.distance)
        ec_r_window = [float("nan")] * len(prof.distance)
        n_series = [0] * len(prof.distance)
        ids_series = [[] for _ in prof.distance]
        tunnel_ec_r = None
        tunnel_ec_n = None
        protein_pdb = None
        shell_a = None

    bneck_r = None if bottleneck is None else bottleneck.get("bottleneck_r")
    bneck_x = None if bottleneck is None else bottleneck.get("bottleneck_x")
    bneck_y = None if bottleneck is None else bottleneck.get("bottleneck_y")
    bneck_z = None if bottleneck is None else bottleneck.get("bottleneck_z")

    fieldnames = [
        "caver_output_dir",
        "tunnel_cluster",
        "snapshot",
        "protein_pdb",
        "shell_a",
        "idx",
        "distance_A",
        "x",
        "y",
        "z",
        col_width,
        "phi_wall_pos",
        "phi_wall_neg",
        "local_ec_r_window",
        "local_hydropathy_mean",
        "local_hydropathy_mean_window",
        "n_residues",
        "residues",
        "tunnel_length_A",
        col_bneck_prof,
        col_bneck_csv,
        "bottleneck_x",
        "bottleneck_y",
        "bottleneck_z",
        "lining_residue_ec_overall",
        "tunnel_ec_r",
        "tunnel_ec_n_pairs",
    ]

    with open(out_csv_path, "w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()
        for i, d0 in enumerate(prof.distance):
            ids = ids_series[i] if i < len(ids_series) else []
            labels_all = [f"{c}:{r}:{name[:3]}" for (c, r, name) in ids]
            if residue_list_max > 0 and len(labels_all) > residue_list_max:
                labels = labels_all[:residue_list_max] + [f"...(+{len(labels_all) - residue_list_max})"]
            else:
                labels = labels_all

            row: dict[str, Any] = {
                "caver_output_dir": caver_output_dir,
                "tunnel_cluster": prof.tunnel_cluster,
                "snapshot": prof.snapshot,
                "protein_pdb": protein_pdb,
                "shell_a": shell_a,
                "idx": i,
                "distance_A": d0,
                "x": prof.x[i],
                "y": prof.y[i],
                "z": prof.z[i],
                col_width: wfac * (prof.r[i] if i < len(prof.r) else float("nan")),
                "phi_wall_pos": phi_wall_pos[i] if i < len(phi_wall_pos) else "",
                "phi_wall_neg": phi_wall_neg[i] if i < len(phi_wall_neg) else "",
                "local_ec_r_window": ec_r_window[i] if i < len(ec_r_window) else "",
                "local_hydropathy_mean": (loc.get("hydropathy_mean") or [])[i]
                if (mapped is not None and i < len(loc.get("hydropathy_mean") or []))
                else "",
                "local_hydropathy_mean_window": (loc.get("hydropathy_mean_window") or [])[i]
                if (mapped is not None and i < len(loc.get("hydropathy_mean_window") or []))
                else "",
                "n_residues": n_series[i] if i < len(n_series) else "",
                "residues": ";".join(labels),
                "tunnel_length_A": prof.length,
                col_bneck_prof: (wfac * prof.bottleneck_radius) if prof.bottleneck_radius is not None else "",
                col_bneck_csv: (wfac * bneck_r) if bneck_r is not None else "",
                "bottleneck_x": bneck_x,
                "bottleneck_y": bneck_y,
                "bottleneck_z": bneck_z,
                "lining_residue_ec_overall": lining_summary["overall_ec"],
                "tunnel_ec_r": tunnel_ec_r,
                "tunnel_ec_n_pairs": tunnel_ec_n,
            }
            w.writerow(row)


def _write_residues_csv(out_csv_path: str, *, mapped: dict[str, Any]) -> None:
    os.makedirs(os.path.dirname(os.path.abspath(out_csv_path)) or ".", exist_ok=True)
    rows = list(mapped.get("mapped_residues") or [])
    fieldnames = ["chain", "res", "aa"]
    with open(out_csv_path, "w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()
        for r in rows:
            w.writerow({k: r.get(k, "") for k in fieldnames})


def _slug_from_caver_base(caver_base: str) -> str:
    """Basename of a Caver output directory, safe for use as a filename suffix."""
    stem = os.path.basename(os.path.abspath(caver_base).rstrip(os.sep)) or "caver_out"
    stem = re.sub(r"[^\w.-]+", "_", stem).strip("._")
    return (stem or "caver_out")[:120]


def _resolve_batch_output_path(
    explicit: str | None,
    *,
    default_path: str,
    n_inputs: int,
    caver_base: str,
) -> str:
    """Single run: use explicit or default. Multiple runs: explicit path gets _<caver_basename> before extension."""
    if explicit is None:
        return default_path
    p = os.path.abspath(os.path.expanduser(str(explicit)))
    if n_inputs <= 1:
        return p
    stem, ext = os.path.splitext(p)
    slug = _slug_from_caver_base(caver_base)
    return f"{stem}_{slug}{ext}" if ext else f"{stem}_{slug}"


def _broadcast_list(values: list[str], n: int, *, label: str) -> list[str | None]:
    if n <= 0:
        return []
    if len(values) == n:
        return list(values)
    if len(values) == 1:
        return [values[0]] * n
    if len(values) == 0:
        return [None] * n
    raise SystemExit(
        f"Error: {label}: expected 0, 1, or {n} value(s) (one per Caver directory), got {len(values)}."
    )


def _run_single_caver(
    args: argparse.Namespace,
    *,
    base: str,
    protein_pdb: str,
    centerline_pdb: str | None,
    n_inputs: int,
) -> None:
    tunnel_cluster = int(args.tunnel)
    his_ec = float(args.his_ec)

    residues_path = _find_analysis_path(base, os.path.join("analysis", "residues.txt"))
    profiles_path = _find_analysis_path(base, os.path.join("analysis", "tunnel_profiles.csv"))
    bottlenecks_path = os.path.join(base, "analysis", "bottlenecks.csv")

    residues = parse_residues_txt(residues_path, tunnel_cluster)
    lining = summarise_lining_residues(residues, his_ec=his_ec)

    # Build the centreline profile from tunnel_profiles.csv (preferred) or a tunnel PDB.
    if os.path.isfile(profiles_path):
        prof = parse_tunnel_profiles_csv(profiles_path, tunnel_cluster)
    elif centerline_pdb:
        xs, ys, zs, rs = _parse_centerline_pdb(os.path.abspath(os.path.expanduser(str(centerline_pdb))))
        dist = _dist_along_polyline(xs, ys, zs)
        prof = TunnelProfile(
            snapshot="(unknown)",
            tunnel_cluster=tunnel_cluster,
            tunnel=tunnel_cluster,
            length=dist[-1] if dist else None,
            bottleneck_radius=min(rs) if rs else None,
            distance=dist,
            x=xs,
            y=ys,
            z=zs,
            r=rs,
        )
    else:
        raise SystemExit("Missing analysis/tunnel_profiles.csv. Pass --centerline-pdb to analyse from a centreline tunnel PDB.")
    bottleneck = None
    if os.path.isfile(bottlenecks_path):
        try:
            bottleneck = parse_bottlenecks_csv(bottlenecks_path, tunnel_cluster, snapshot=prof.snapshot)
        except Exception:
            bottleneck = None

    rmin = min(prof.r) if prof.r else None
    rmax = max(prof.r) if prof.r else None
    bneck_prof = prof.bottleneck_radius
    if bool(getattr(args, "diameter", False)):
        width_summary: dict[str, Any] = {
            "output_tunnel_width": "diameter",
            "profile_diameter_min_A": (2.0 * rmin) if rmin is not None else None,
            "profile_diameter_max_A": (2.0 * rmax) if rmax is not None else None,
            "bottleneck_diameter_A_from_profiles": (2.0 * bneck_prof) if bneck_prof is not None else None,
        }
    else:
        width_summary = {
            "output_tunnel_width": "radius",
            "profile_radius_min_A": rmin,
            "profile_radius_max_A": rmax,
            "bottleneck_from_profiles_A": bneck_prof,
        }

    summary: dict[str, Any] = {
        "caver_output_dir": base,
        "tunnel_cluster": tunnel_cluster,
        "snapshot": prof.snapshot,
        "tunnel_length_A": prof.length,
        "profile_n_points": len(prof.distance),
        "profile_distance_range_A": [min(prof.distance), max(prof.distance)] if prof.distance else None,
        **width_summary,
        "bottleneck": bottleneck,
        "lining": {
            "n_lining_residues": lining["n_lining_residues"],
            "overall_ec": lining["overall_ec"],
            "mean_ec_per_residue": lining["mean_ec_per_residue"],
            "mean_hydropathy_kyte_doolittle": lining["mean_hydropathy_kyte_doolittle"],
            "n_hydropathy_defined": lining["n_hydropathy_defined"],
            "n_hydropathy_missing": lining["n_hydropathy_missing"],
        },
    }

    mapped = map_profile_to_local_properties(
        prof,
        os.path.abspath(os.path.expanduser(str(protein_pdb))),
        shell_a=float(args.shell_a),
        ec_epsilon_r=float(args.ec_epsilon_r),
        ec_r_min=float(args.ec_r_min),
        local_window=int(getattr(args, "local_window", 0) or 0),
    )
    summary["mapped"] = {
        "protein_pdb": mapped["protein_pdb"],
        "shell_a": mapped["shell_a"],
        "overall_from_mapped_residues": mapped["overall_from_mapped_residues"],
        "tunnel_ec": mapped.get("tunnel_ec"),
    }

    norm_for_plot: Any = None
    cbar_label = "EC r"
    smooth_win = int(getattr(args, "colour_smooth_window", 0) or 0)
    if str(args.color_by) == "hydropathy":
        win = int(getattr(args, "local_window", 0) or 0)
        if win >= 3:
            color_values = list(mapped["local"].get("hydropathy_mean_window") or [])
            title = f"Tunnel {tunnel_cluster} hydropathy (local, window={win})"
            cbar_label = "Kyte–Doolittle hydropathy (local, rolling window)"
        else:
            color_values = list(mapped["local"].get("hydropathy_mean") or [])
            title = f"Tunnel {tunnel_cluster} hydropathy (Kyte–Doolittle mean)"
            cbar_label = "Kyte–Doolittle hydropathy (local mean)"
        if smooth_win >= 3:
            color_values = _rolling_mean_ignore_nan([float(x) for x in color_values], window=smooth_win)
        # Orange (hydrophobic) -> white (0) -> green (hydrophilic)
        try:
            from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm  # type: ignore

            cmap = LinearSegmentedColormap.from_list("hydro_owg", ["#2ca25f", "#ffffff", "#f28e2b"])
            if str(getattr(args, "color_scale", "fixed")) == "data":
                norm_for_plot, vmin, vmax = _diverging_norm_from_series(
                    [float(x) for x in color_values], fallback=(-4.5, 4.5)
                )
            else:
                norm_for_plot, vmin, vmax = _diverging_norm_fixed(-4.5, 4.5)
        except Exception:
            cmap = "PiYG"
            vmin, vmax = -4.5, 4.5
            norm_for_plot = None
    else:
        win = int(getattr(args, "local_window", 0) or 0)
        if win >= 3:
            color_values = list(mapped["local"].get("ec_r_window") or [])
            title = f"Tunnel {tunnel_cluster} EC (local, window={win})"
            cbar_label = "EC r (local, rolling window)"
        else:
            r0 = (mapped.get("tunnel_ec") or {}).get("r")
            color_values = [float(r0) if r0 is not None else float("nan")] * len(prof.distance)
            title = f"Tunnel {tunnel_cluster} EC (global r={r0})"
            cbar_label = "EC r (global)"
        if smooth_win >= 3:
            color_values = _rolling_mean_ignore_nan([float(x) for x in color_values], window=smooth_win)
        # Red (negative) -> white (0) -> blue (positive); scale from segment means (same as draw order)
        _cv = [float(x) for x in color_values]
        if len(_cv) >= 2:
            v_plot = [(_cv[i] + _cv[i + 1]) / 2.0 for i in range(len(_cv) - 1)]
        else:
            v_plot = _cv
        if str(getattr(args, "color_scale", "fixed")) == "fixed":
            norm_for_plot, vmin, vmax = _diverging_norm_fixed(-1.0, 1.0)
        else:
            norm_for_plot, vmin, vmax = (
                _ec_diverging_norm_from_values(v_plot) if v_plot else _ec_diverging_norm_from_values(_cv)
            )
        cmap = _ec_plot_cmap_rwb()

    default_plot = os.path.join(base, f"tunnel_{tunnel_cluster}_ec.png")
    plot_out = _resolve_batch_output_path(
        args.plot_out, default_path=default_plot, n_inputs=n_inputs, caver_base=base
    )
    default_csv = os.path.join(base, f"tunnel_{tunnel_cluster}_points.csv")
    output_csv = _resolve_batch_output_path(
        args.output_csv, default_path=default_csv, n_inputs=n_inputs, caver_base=base
    )

    if not args.no_plot:
        import numpy as np

        def _fmt_res_list(idx: int) -> str:
            if mapped is None:
                return ""
            ids = mapped["local"]["residue_ids"][idx] if idx < len(mapped["local"]["residue_ids"]) else []
            labels_all = [f"{c}:{r}:{name[:3]}" for (c, r, name) in ids]
            max_list = int(getattr(args, "residue_list_max", 30))
            if max_list > 0 and len(labels_all) > max_list:
                labels = labels_all[:max_list] + [f"...(+{len(labels_all) - max_list})"]
            else:
                labels = labels_all
            return ";".join(labels)

        idx_top = 0
        idx_bottom = max(0, len(prof.distance) - 1)
        idx_bneck = int(np.argmin(np.asarray(prof.r, dtype=float))) if prof.r else 0
        ann = {
            "points": [
                {"label": "top", "idx": idx_top, "residues": _fmt_res_list(idx_top)},
                {"label": "bottleneck", "idx": idx_bneck, "residues": _fmt_res_list(idx_bneck)},
                {"label": "bottom", "idx": idx_bottom, "residues": _fmt_res_list(idx_bottom)},
            ]
        }
        plot_profile(
            prof,
            out_path=str(plot_out),
            cbar_label=cbar_label,
            color_values=color_values,
            title=title,
            cmap=cmap,
            vmin=vmin,
            vmax=vmax,
            norm=norm_for_plot,
            annotate=ann,
            upsample_factor=float(getattr(args, "plot_upsample", 6.0)),
            upsample_max_points=int(getattr(args, "plot_upsample_max", 800)),
            as_diameter=bool(getattr(args, "diameter", False)),
            rasterise_fill=bool(getattr(args, "rasterise_fill", False)),
            rasterise_fill_dpi=int(getattr(args, "rasterise_fill_dpi", 300)),
            xlim=tuple(getattr(args, "xlim")) if getattr(args, "xlim", None) is not None else None,
            invert_y=bool(getattr(args, "invert_y", True)),
        )
        summary["plot"] = {
            "path": os.path.abspath(str(plot_out)),
            "color_by": str(args.color_by),
            "cbar_label": cbar_label,
            "plot_upsample": float(getattr(args, "plot_upsample", 6.0)),
            "plot_upsample_max": int(getattr(args, "plot_upsample_max", 800)),
            "as_diameter": bool(getattr(args, "diameter", False)),
        }
        if str(args.color_by) == "ec":
            summary["plot"]["ec_color_vmin"] = vmin
            summary["plot"]["ec_color_vmax"] = vmax

    _write_points_csv(
        str(output_csv),
        prof=prof,
        caver_output_dir=base,
        bottleneck=bottleneck,
        lining_summary={"overall_ec": lining.get("overall_ec", 0.0)},
        mapped=mapped,
        residue_list_max=int(getattr(args, "residue_list_max", 30)),
        as_diameter=bool(getattr(args, "diameter", False)),
    )
    summary["points_csv"] = os.path.abspath(str(output_csv))

    if args.residues_csv:
        residues_path_out = _resolve_batch_output_path(
            args.residues_csv,
            default_path=os.path.join(base, "mapped_residues.csv"),
            n_inputs=n_inputs,
            caver_base=base,
        )
        _write_residues_csv(str(residues_path_out), mapped=mapped)
        summary["mapped_residues_csv"] = os.path.abspath(str(residues_path_out))

    if args.json_out:
        json_path_out = _resolve_batch_output_path(
            args.json_out,
            default_path=os.path.join(base, f"tunnel_{tunnel_cluster}_summary.json"),
            n_inputs=n_inputs,
            caver_base=base,
        )
        os.makedirs(os.path.dirname(os.path.abspath(str(json_path_out))) or ".", exist_ok=True)
        with open(str(json_path_out), "w", encoding="utf-8", newline="\n") as f:
            json.dump({**summary, "lining_residues": lining["lining_residues"]}, f, indent=2, sort_keys=True)
            f.write("\n")

    print(f"[caver_tunnel_analysis] wrote {os.path.abspath(str(output_csv))}")
    if not args.no_plot:
        print(f"[caver_tunnel_analysis] wrote {os.path.abspath(str(plot_out))}")


def main() -> None:
    ap = argparse.ArgumentParser(
        description="Compute tunnel electrostatic complementarity (EC) and generate a 2D tunnel profile figure.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    ap.add_argument(
        "caver_output_dirs",
        nargs="+",
        metavar="CAVER_DIR",
        help="One or more Caver output directories (each contains analysis/).",
    )
    ap.add_argument("--tunnel", type=int, required=True, help="Tunnel cluster number (e.g. 26).")
    ap.add_argument(
        "--color-by",
        choices=("ec", "hydropathy"),
        default="ec",
        help="Colour the plot by EC or Kyte–Doolittle hydropathy.",
    )
    ap.add_argument(
        "--color-scale",
        choices=("fixed", "data"),
        default="fixed",
        help=(
            "Colour scale limits. fixed: use reference limits (EC: [-1, +1]; hydropathy: [-4.5, +4.5]). "
            "data: set limits from the plotted values (diverging scales remain centred at 0)."
        ),
    )
    ap.add_argument(
        "--colour-scale",
        choices=("fixed", "data"),
        default=None,
        dest="color_scale",
        help="Alias for --color-scale.",
    )
    ap.add_argument(
        "--his-ec",
        type=float,
        default=0.1,
        dest="his_ec",
        help=(
            "HIS partial EC value in the residues.txt provenance sum (HIS=+1 uses this weight). "
            "EC wall potentials in mapping use full EC sign rules (HIS=+1 to match interface EC)."
        ),
    )
    ap.add_argument(
        "--local-window",
        type=int,
        default=0,
        metavar="N",
        help=(
            "Window size (odd integer) for rolling 'local' values along the tunnel when colouring. "
            "Applies to EC correlation and hydropathy. 0 disables."
        ),
    )
    ap.add_argument(
        "--colour-smooth-window",
        type=int,
        default=0,
        metavar="N",
        dest="colour_smooth_window",
        help=(
            "Optional additional rolling-window smoothing (odd integer) applied to the per-point colour series "
            "for plotting only (does not change CSV/JSON). 0 disables."
        ),
    )
    ap.add_argument(
        "--color-smooth-window",
        type=int,
        default=None,
        metavar="N",
        dest="colour_smooth_window",
        help="Alias for --colour-smooth-window.",
    )
    ap.add_argument(
        "--rasterise-fill",
        action="store_true",
        help=(
            "For SVG/PDF output, rasterise the tunnel interior colour fill (high-DPI bitmap) while keeping "
            "axes, text, and outlines as vector. This can substantially reduce file size for highly upsampled plots. "
            "This option does not change the sampling used to construct the plotted tunnel."
        ),
    )
    ap.add_argument(
        "--rasterise-fill-dpi",
        type=int,
        default=300,
        metavar="DPI",
        help="DPI used for the rasterised tunnel fill in SVG/PDF when --rasterise-fill is enabled.",
    )
    ap.add_argument(
        "--xlim",
        nargs=2,
        type=float,
        metavar=("XMIN", "XMAX"),
        default=None,
        help=(
            "Optional fixed x-axis limits for the profile plot (e.g. --xlim -21 21). "
            "When set, this overrides the automatic radius/diameter-based range."
        ),
    )
    ap.add_argument(
        "--invert-y",
        action="store_true",
        default=True,
        dest="invert_y",
        help="Invert the profile y-axis so the tunnel top is at the top of the figure (default).",
    )
    ap.add_argument(
        "--no-invert-y",
        action="store_false",
        dest="invert_y",
        help="Do not invert the profile y-axis.",
    )
    ap.add_argument(
        "--ec-epsilon-r",
        type=float,
        default=4.0,
        help="Relative dielectric used for tunnel EC wall potentials (matches EC default: 4.0).",
    )
    ap.add_argument(
        "--ec-r-min",
        type=float,
        default=1.0,
        help="Minimum distance clamp (Å) for tunnel EC wall potentials (matches EC default: 1.0).",
    )
    ap.add_argument(
        "--protein-pdb",
        nargs="+",
        required=True,
        metavar="PDB",
        help=(
            "Protein structure PDB path(s) for EC mapping: one file (used for every Caver directory) "
            "or one file per directory in the same order as CAVER_DIR arguments."
        ),
    )
    ap.add_argument(
        "--shell-a",
        type=float,
        default=3.0,
        help="Extra shell added to tunnel radius when selecting nearby residues (Å).",
    )
    ap.add_argument(
        "--residue-list-max",
        type=int,
        default=30,
        help="Max residues to write per profile point (0 = no limit).",
    )
    ap.add_argument("--no-plot", action="store_true", help="Skip generating plots.")
    ap.add_argument(
        "--diameter",
        action="store_true",
        help=(
            "Use tunnel opening diameter (2× Caver radius) for the profile plot and for width columns in CSV/JSON "
            "(diameter_A, bottleneck_diameter_*, etc.)."
        ),
    )
    ap.add_argument(
        "--plot-upsample",
        type=float,
        default=6.0,
        metavar="FACTOR",
        dest="plot_upsample",
        help=(
            "Resample the tunnel along distance for the plot: effective point count is multiplied by this "
            "(1 = use Caver points as-is; larger = smoother fill, capped by --plot-upsample-max). "
            "Radius uses a monotonic spline when SciPy is available, else linear. "
            "This option affects the plotted geometry; it is independent of --rasterise-fill."
        ),
    )
    ap.add_argument(
        "--plot-upsample-max",
        type=int,
        default=800,
        metavar="N",
        dest="plot_upsample_max",
        help="Maximum number of distance samples in the profile plot (limits cost for long tunnels).",
    )
    ap.add_argument(
        "--plot-out",
        default=None,
        help=(
            "Output plot path. Extension selects format (e.g. .png, .pdf, .svg). "
            "Default per directory: <dir>/tunnel_<N>_ec.png. With multiple CAVER_DIR, a shared path "
            "gets _<basename> inserted before the extension."
        ),
    )
    ap.add_argument(
        "--output-csv",
        "-o",
        default=None,
        metavar="FILE",
        help=(
            "Per-point CSV path. Default per directory: <dir>/tunnel_<N>_points.csv. "
            "With multiple CAVER_DIR, a shared path gets _<basename> before the extension."
        ),
    )
    ap.add_argument(
        "--residues-csv",
        default=None,
        metavar="FILE",
        help="Optional mapped-residue CSV; with multiple CAVER_DIR, path is suffixed like --output-csv.",
    )
    ap.add_argument(
        "--json-out",
        default=None,
        help="Optional JSON summary path; with multiple CAVER_DIR, path is suffixed like --output-csv.",
    )
    ap.add_argument(
        "--centerline-pdb",
        action="append",
        default=None,
        metavar="PDB",
        help=(
            "Optional centreline tunnel PDB if analysis/tunnel_profiles.csv is missing. "
            "Repeat once per Caver directory (same order) or once for all runs."
        ),
    )
    args = ap.parse_args()

    caver_bases = [os.path.abspath(os.path.expanduser(p)) for p in args.caver_output_dirs]
    for b in caver_bases:
        if not os.path.isdir(b):
            raise SystemExit(f"Error: not a directory: {b}")

    n = len(caver_bases)
    proteins = [os.path.abspath(os.path.expanduser(p)) for p in (args.protein_pdb or [])]
    if len(proteins) == 1:
        proteins = proteins * n
    elif len(proteins) != n:
        raise SystemExit(f"Error: --protein-pdb: expected 1 or {n} path(s), got {len(proteins)}.")

    cl_list_in = [str(x) for x in (args.centerline_pdb or [])]
    centerlines: list[str | None] = _broadcast_list(cl_list_in, n, label="--centerline-pdb")

    for base, pdb, cl in zip(caver_bases, proteins, centerlines):
        _run_single_caver(args, base=base, protein_pdb=pdb, centerline_pdb=cl, n_inputs=n)


if __name__ == "__main__":
    main()

