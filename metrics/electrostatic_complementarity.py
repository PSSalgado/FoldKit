#!/usr/bin/env python3
"""
Electrostatic complementarity (EC)
=================================

FoldKit uses an electrostatic complementarity (EC) score: a Pearson correlation (r)
between electrostatic potentials evaluated at *facing* points across an
interface, with opposite-sign potentials indicating complementarity.

Method reference (definition of EC as a correlation across facing surface points):
- McCoy, A. J.; Epa, V. C.; Colman, P. M. (1997). *Electrostatic complementarity at protein/protein interfaces.*
  Journal of Molecular Biology 268, 570–584. DOI: 10.1006/jmbi.1997.0987.

This file now provides a lightweight, self-contained EC implementation that is:
- dependency-light (uses Biopython if available, like the rest of FoldKit),
- deterministic,
- fast enough for typical ASU/lattice interface sizes.

Important limitation: PDB/CIF inputs generally lack partial charges, so this
implementation approximates charges using residue-level sign
(ASP/GLU=-1, ARG/LYS/HIS=+1) distributed evenly across a residue's heavy atoms.
"""

from __future__ import annotations

import argparse
import math
import os
import sys
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np

_REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if _REPO_ROOT not in sys.path:
    sys.path.insert(0, _REPO_ROOT)

from cli_log import add_log_args, setup_log_from_args

try:
    from Bio.PDB.PDBParser import PDBParser
    from Bio.PDB.NeighborSearch import NeighborSearch

    BIOPYTHON_AVAILABLE = True
except Exception:
    PDBParser = None  # type: ignore[assignment]
    NeighborSearch = None  # type: ignore[assignment]
    BIOPYTHON_AVAILABLE = False


def _residue_charge_sign(resname: str) -> int:
    name = (resname or "").strip().upper()
    if name in ("ASP", "GLU"):
        return -1
    if name in ("ARG", "LYS", "HIS"):
        return 1
    return 0


def _heavy_atoms(residue) -> List:
    out = []
    for a in residue.get_atoms():
        el = (getattr(a, "element", "") or "").strip().upper()
        if el == "H":
            continue
        out.append(a)
    return out


def _charges_for_chain_atoms(chain) -> Tuple[np.ndarray, np.ndarray]:
    """
    Return (coords[N,3], charges[N]) for heavy atoms in `chain`.
    Charge model: residue sign distributed uniformly across heavy atoms.
    """
    coords = []
    charges = []
    for res in chain.get_residues():
        q = _residue_charge_sign(getattr(res, "resname", ""))
        if q == 0:
            continue
        atoms = _heavy_atoms(res)
        if not atoms:
            continue
        per_atom = float(q) / float(len(atoms))
        for a in atoms:
            coords.append(np.asarray(a.coord, dtype=float))
            charges.append(per_atom)
    if not coords:
        return np.zeros((0, 3), dtype=float), np.zeros((0,), dtype=float)
    return np.vstack(coords), np.asarray(charges, dtype=float)


def _contact_atoms(chain_a, chain_b, contact_distance: float) -> Tuple[List, List]:
    """
    Return contact heavy atoms on A and B within `contact_distance` (Å), based on
    any atom–atom proximity. Uses NeighborSearch for speed when available.
    """
    atoms_a = [a for a in chain_a.get_atoms() if (getattr(a, "element", "") or "").strip().upper() != "H"]
    atoms_b = [a for a in chain_b.get_atoms() if (getattr(a, "element", "") or "").strip().upper() != "H"]
    if not atoms_a or not atoms_b:
        return [], []

    if NeighborSearch is None:
        a_hits = set()
        b_hits = set()
        for a in atoms_a:
            pa = np.asarray(a.coord, dtype=float)
            for b in atoms_b:
                d = float(np.linalg.norm(pa - np.asarray(b.coord, dtype=float)))
                if d <= float(contact_distance):
                    a_hits.add(a)
                    b_hits.add(b)
        return list(a_hits), list(b_hits)

    ns_b = NeighborSearch(atoms_b)
    a_hits = []
    b_hits = set()
    for a in atoms_a:
        nearby = ns_b.search(a.coord, float(contact_distance), level="A")
        if nearby:
            a_hits.append(a)
            for b in nearby:
                b_hits.add(b)
    return a_hits, list(b_hits)


def _pair_facing_points(points_a: np.ndarray, points_b: np.ndarray) -> Optional[np.ndarray]:
    """
    For each point in A, find the index of nearest neighbour in B.
    Returns idx array of shape (n_a,), or None if pairing isn't possible.
    """
    if points_a.size == 0 or points_b.size == 0:
        return None
    idx = np.empty((points_a.shape[0],), dtype=int)
    for i in range(points_a.shape[0]):
        d2 = np.sum((points_b - points_a[i]) ** 2, axis=1)
        idx[i] = int(np.argmin(d2))
    return idx


def _coulomb_potential(
    sample_points: np.ndarray,
    charge_coords: np.ndarray,
    charge_values: np.ndarray,
    *,
    epsilon_r: float = 4.0,
    r_min: float = 1.0,
) -> np.ndarray:
    """
    Simple Coulomb sum with softening:
      phi(p) = sum_k q_k / (epsilon_r * max(r_min, ||p - r_k||))
    """
    if sample_points.size == 0 or charge_coords.size == 0 or charge_values.size == 0:
        return np.zeros((sample_points.shape[0],), dtype=float)
    eps = float(epsilon_r) if float(epsilon_r) > 0 else 4.0
    rmin = float(r_min) if float(r_min) > 0 else 1.0
    out = np.zeros((sample_points.shape[0],), dtype=float)
    for i, p in enumerate(sample_points):
        r = np.linalg.norm(charge_coords - p[None, :], axis=1)
        r = np.maximum(r, rmin)
        out[i] = float(np.sum(charge_values / (eps * r)))
    return out


def _pearson_r(x: np.ndarray, y: np.ndarray) -> Optional[float]:
    if x.size == 0 or y.size == 0 or x.size != y.size:
        return None
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    x0 = x - float(np.mean(x))
    y0 = y - float(np.mean(y))
    den = float(np.sqrt(np.sum(x0 * x0) * np.sum(y0 * y0)))
    # If there is no variance in one or both vectors, the correlation is undefined.
    # For FoldKit reporting we return 0.0 (no linear relationship) instead of None,
    # so callers don't drop EC fields and fall back to charge-complementarity.
    if not math.isfinite(den) or den <= 0.0:
        return 0.0
    r = float(np.sum(x0 * y0) / den)
    return max(-1.0, min(1.0, r))


def compute_ec_complementarity_detailed(
    chain_a,
    chain_b,
    *,
    contact_distance: float = 5.0,
    epsilon_r: float = 4.0,
    r_min: float = 1.0,
) -> Dict[str, object]:
    """
    Compute an EC correlation score for a chain pair.

    This follows the facing-point correlation definition described by
    McCoy, Epa & Colman (1997) (J Mol Biol 268:570–584; DOI: 10.1006/jmbi.1997.0987).

    Returns dict with keys:
    - r: correlation coefficient (float) or None if not computable
    - n_pairs: number of facing point pairs used (int)
    - details: small diagnostics dict
    """
    a_atoms, b_atoms = _contact_atoms(chain_a, chain_b, float(contact_distance))
    if not a_atoms or not b_atoms:
        return {"r": None, "n_pairs": 0, "details": {"reason": "no_contact_atoms"}}

    pts_a = np.vstack([np.asarray(a.coord, dtype=float) for a in a_atoms])
    pts_b = np.vstack([np.asarray(b.coord, dtype=float) for b in b_atoms])
    nn = _pair_facing_points(pts_a, pts_b)
    if nn is None:
        return {"r": None, "n_pairs": 0, "details": {"reason": "pairing_failed"}}

    coords_a, q_a = _charges_for_chain_atoms(chain_a)
    coords_b, q_b = _charges_for_chain_atoms(chain_b)

    phi_a = _coulomb_potential(pts_a, coords_b, q_b, epsilon_r=epsilon_r, r_min=r_min)
    phi_b = _coulomb_potential(pts_b[nn], coords_a, q_a, epsilon_r=epsilon_r, r_min=r_min)
    r = _pearson_r(phi_a, -phi_b)
    return {
        "r": r,
        "n_pairs": int(phi_a.size),
        "details": {
            "n_contact_atoms_a": int(len(a_atoms)),
            "n_contact_atoms_b": int(len(b_atoms)),
            "epsilon_r": float(epsilon_r),
            "r_min": float(r_min),
            "charge_model": "residue_sign_uniform_over_heavy_atoms",
        },
    }


def compute_ec_complementarity(chain_a, chain_b, *, contact_distance: float = 5.0) -> Optional[float]:
    """Convenience wrapper returning just `r` (or None)."""
    return compute_ec_complementarity_detailed(chain_a, chain_b, contact_distance=contact_distance).get("r")  # type: ignore[return-value]


def analyze_structure_ec(
    structure_path: str,
    *,
    focus_chains: Optional[Sequence[str]] = None,
    mode: str = "full",
    contact_distance: float = 5.0,
) -> Dict[str, object]:
    """
    Structure-level driver used by the CLI.
    Returns a dict with `interfaces`: list of chain-pair results.
    """
    if not BIOPYTHON_AVAILABLE or PDBParser is None:
        return {"error": "BioPython not available", "interfaces": []}
    if str(mode).strip().lower() != "full":
        return {"error": "Only mode='full' is supported", "interfaces": []}

    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("_ec", structure_path)
    chains = list(structure.get_chains())
    if focus_chains:
        fs = {str(c).strip() for c in focus_chains if str(c).strip()}
        chains = [c for c in chains if c.id in fs]

    interfaces = []
    for i in range(len(chains)):
        for j in range(i + 1, len(chains)):
            ca = chains[i]
            cb = chains[j]
            d = compute_ec_complementarity_detailed(ca, cb, contact_distance=float(contact_distance))
            interfaces.append(
                {
                    "chain_a": ca.id,
                    "chain_b": cb.id,
                    "ec_r": d.get("r"),
                    "ec_n_pairs": d.get("n_pairs", 0),
                }
            )
    return {"interfaces": interfaces}


__all__ = [
    "compute_ec_complementarity",
    "compute_ec_complementarity_detailed",
    "analyze_structure_ec",
]


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Compute electrostatic complementarity (EC) for interfaces in a structure file.",
    )
    parser.add_argument(
        "input",
        help="Input PDB/CIF file.",
    )
    parser.add_argument(
        "--chains",
        metavar="IDS",
        help=(
            "Optional chain IDs to include (comma-separated). If omitted, all chains are considered."
        ),
    )
    parser.add_argument(
        "--mode",
        choices=("full", "fast"),
        default="full",
        help=(
            "Backend mode. 'full' computes EC as a facing-point Pearson correlation as in "
            "McCoy, Epa & Colman (1997) J Mol Biol 268:570–584 (DOI: 10.1006/jmbi.1997.0987)."
        ),
    )
    parser.add_argument(
        "--output",
        "-o",
        metavar="FILE",
        help="Output file. If omitted, write to stdout.",
    )
    add_log_args(parser)
    args = parser.parse_args()
    summary_log = setup_log_from_args(
        args,
        script_path=__file__,
        inputs=[getattr(args, "input", "")],
        pattern=None,
    )

    focus_chains = None
    if getattr(args, "chains", None):
        focus_chains = [c.strip() for c in str(args.chains).split(",") if c.strip()]

    out = sys.stdout
    if args.output:
        out = open(args.output, "w")

    try:
        if summary_log is not None:
            summary_log.task("Electrostatic complementarity (EC) analysis")
            summary_log.task(f"Mode: {str(args.mode).strip().lower()}")
            summary_log.task(f"Focus chains: {','.join(focus_chains) if focus_chains else 'ALL'}")
        results = analyze_structure_ec(
            args.input,
            focus_chains=focus_chains,
            mode=str(args.mode).strip().lower(),
        )

        # Minimal, stable text output: one line per interface.
        out.write("ELECTROSTATIC COMPLEMENTARITY (EC)\n")
        out.write(f"Input: {os.path.basename(args.input)}\n\n")

        interfaces = results.get("interfaces", []) if isinstance(results, dict) else []
        if not interfaces:
            out.write("No interfaces reported.\n")
            return

        out.write("Interface\tEC_r\tn_pairs\n")
        for row in interfaces:
            a = row.get("chain_a", "?")
            b = row.get("chain_b", "?")
            r = row.get("ec_r", None)
            n = row.get("ec_n_pairs", None)
            r_str = "NA" if r is None else f"{float(r):.3f}"
            n_str = "NA" if n is None else str(int(n))
            out.write(f"{a}-{b}\t{r_str}\t{n_str}\n")
    finally:
        if out is not sys.stdout:
            out.close()


if __name__ == "__main__":
    main()

