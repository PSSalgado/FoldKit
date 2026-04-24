#!/usr/bin/env python3
"""
Lattice packing analyser
========================

Compute packing-style metrics for a symmetry-expanded coordinate model (a supercell
or lattice fragment) where many copies of a molecule are present in one Cartesian
frame (typically P1).

This script is intended for inputs like `supercell.pdb` with a user-chosen number
of chain copies (for example 12 chains arranged in a 3×4 grid). It does not apply
space-group operators; it analyses the coordinates exactly as supplied.

Volume definition
-----------------
Packing-like metrics require a volume. Two options are supported:

- **CRYST1 volume**: use the unit cell from the PDB `CRYST1` record (if present).
- **Bounding-box volume**: compute the axis-aligned bounding box of all atoms and
  use its volume (optionally padded by a constant margin, Å).

The default behaviour is:
- use CRYST1 volume if present;
- otherwise fall back to the bounding box.
"""

from __future__ import annotations

import argparse
import json
import math
import os
import sys
from datetime import datetime

_REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if _REPO_ROOT not in sys.path:
    sys.path.insert(0, _REPO_ROOT)

import numpy as np

try:
    from Bio.PDB.PDBParser import PDBParser
    from Bio.PDB.PDBExceptions import PDBConstructionWarning
    import warnings

    warnings.simplefilter("ignore", PDBConstructionWarning)
    BIOPYTHON_AVAILABLE = True
except Exception:
    PDBParser = None  # type: ignore[assignment]
    BIOPYTHON_AVAILABLE = False

from cli_log import add_log_args, setup_log_from_args
from metrics.interface_analyser_base import collect_structure_paths, filter_paths_by_patterns


_ATOMIC_VOLUMES_A3 = {
    "C": 16.44,
    "N": 11.51,
    "O": 9.13,
    "S": 19.86,
    "P": 17.5,
    "H": 5.15,
    "CA": 25.99,
    "MG": 22.57,
    "ZN": 14.6,
    "FE": 11.8,
    "MN": 15.6,
    "CU": 11.8,
    "NA": 23.7,
    "K": 43.2,
    "CL": 22.7,
    "BR": 26.5,
}

_ATOMIC_WEIGHTS_DA = {
    "C": 12.011,
    "N": 14.007,
    "O": 15.999,
    "S": 32.06,
    "P": 30.974,
    "H": 1.008,
    "CA": 40.078,
    "MG": 24.305,
    "ZN": 65.38,
    "FE": 55.845,
    "MN": 54.938,
    "CU": 63.546,
    "NA": 22.99,
    "K": 39.098,
    "CL": 35.45,
    "BR": 79.904,
}


def _read_cryst1(pdb_path: str) -> dict | None:
    """
    Parse PDB CRYST1 fields (a,b,c,alpha,beta,gamma,space_group).
    Returns None if no CRYST1 record exists or parsing fails.
    """
    try:
        with open(pdb_path, "r", encoding="utf-8", errors="replace") as f:
            for line in f:
                if not line.startswith("CRYST1"):
                    continue
                a = float(line[6:15].strip())
                b = float(line[15:24].strip())
                c = float(line[24:33].strip())
                alpha = float(line[33:40].strip())
                beta = float(line[40:47].strip())
                gamma = float(line[47:54].strip())
                sg = (line[55:66].strip() or "P1").strip()
                return {"a": a, "b": b, "c": c, "alpha": alpha, "beta": beta, "gamma": gamma, "space_group": sg}
    except Exception:
        return None
    return None


def _unit_cell_volume_a3(cryst1: dict) -> float:
    a = float(cryst1["a"])
    b = float(cryst1["b"])
    c = float(cryst1["c"])
    alpha, beta, gamma = np.radians([float(cryst1["alpha"]), float(cryst1["beta"]), float(cryst1["gamma"])])
    v = a * b * c * float(
        np.sqrt(
            1.0
            + 2.0 * np.cos(alpha) * np.cos(beta) * np.cos(gamma)
            - np.cos(alpha) ** 2
            - np.cos(beta) ** 2
            - np.cos(gamma) ** 2
        )
    )
    return float(v)


def _bbox_volume_a3(coords: np.ndarray, pad_a: float = 0.0) -> dict:
    if coords.size == 0:
        return {"volume_a3": 0.0, "min": None, "max": None, "pad_a": float(pad_a)}
    pad = float(max(0.0, pad_a))
    cmin = np.min(coords, axis=0) - pad
    cmax = np.max(coords, axis=0) + pad
    span = np.maximum(cmax - cmin, 0.0)
    vol = float(span[0] * span[1] * span[2])
    return {
        "volume_a3": vol,
        "min": [float(x) for x in cmin],
        "max": [float(x) for x in cmax],
        "span": [float(x) for x in span],
        "pad_a": pad,
    }


def _element_symbol(atom) -> str:
    el = (getattr(atom, "element", "") or "").strip()
    if el:
        return el.upper()
    name = (getattr(atom, "name", "") or "").strip()
    return (name[:1] or "C").upper()


def analyse_lattice_packing(
    pdb_path: str,
    *,
    expected_chains: int | None = None,
    require_cryst1: bool = False,
    volume_source: str = "auto",
    bbox_pad_a: float = 0.0,
    allow_cryst1_mismatch: bool = False,
) -> dict:
    if not BIOPYTHON_AVAILABLE or PDBParser is None:
        return {"error": "BioPython not available for structure parsing"}

    cryst1 = _read_cryst1(pdb_path)
    if require_cryst1 and not cryst1:
        return {"error": "CRYST1 record required but not found"}

    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("lattice", pdb_path)

    chains = list(structure.get_chains())
    chain_count = len(chains)
    if expected_chains is not None and int(expected_chains) != chain_count:
        return {
            "error": f"Expected {int(expected_chains)} chain(s) but found {chain_count}",
            "chain_count": chain_count,
        }

    coords = []
    total_atoms = 0
    residue_count = 0
    element_counts: dict[str, int] = {}
    total_mass_da = 0.0
    total_atomic_volume_a3 = 0.0

    for chain in chains:
        for residue in chain:
            residue_count += 1
            for atom in residue:
                total_atoms += 1
                el = _element_symbol(atom)
                element_counts[el] = element_counts.get(el, 0) + 1
                total_mass_da += _ATOMIC_WEIGHTS_DA.get(el, _ATOMIC_WEIGHTS_DA["C"])
                total_atomic_volume_a3 += _ATOMIC_VOLUMES_A3.get(el, _ATOMIC_VOLUMES_A3["C"])
                try:
                    coords.append(np.asarray(atom.coord, dtype=float))
                except Exception:
                    pass

    coords_arr = np.vstack(coords) if coords else np.zeros((0, 3), dtype=float)

    # Choose volume.
    #
    # NOTE: `vol_meta` is used to carry warnings (and defaults) through the
    # volume selection logic. Be careful with dict merge order: defaults must
    # never overwrite computed values like bbox/cryst1 volume.
    vol_meta = {"source": None, "volume_a3": None}
    source = (volume_source or "auto").strip().lower()
    if source not in ("auto", "cryst1", "bbox"):
        return {"error": "volume_source must be one of: auto, cryst1, bbox"}

    # Guard against a common mismatch: multi-copy supercell coordinates but CRYST1
    # describes the original (smaller) crystallographic unit cell.
    if cryst1 and (expected_chains is not None) and int(expected_chains) > 1 and not allow_cryst1_mismatch:
        if source == "cryst1":
            return {
                "error": (
                    "volume_source=cryst1 with expected_chains>1 is likely inconsistent: the CRYST1 record "
                    "often describes the original unit cell, while this file contains a multi-copy supercell. "
                    "Use --volume-source bbox (recommended) or pass --allow-cryst1-mismatch to override."
                )
            }
        if source == "auto":
            # Fall back to bbox for consistency, but record a warning.
            source = "bbox"
            vol_meta["warning"] = (
                "CRYST1 record present but expected_chains>1. Assuming CRYST1 describes the original unit cell; "
                "using bounding-box volume instead. To force CRYST1 volume, pass --volume-source cryst1 "
                "and --allow-cryst1-mismatch."
            )

    if source in ("auto", "cryst1") and cryst1:
        v = _unit_cell_volume_a3(cryst1)
        # Merge defaults first so they can't overwrite computed volume fields.
        vol_meta = {**vol_meta, "source": "cryst1", "volume_a3": float(v), "cryst1": cryst1}
    elif source == "cryst1" and not cryst1:
        return {"error": "volume_source=cryst1 but CRYST1 record not found"}
    else:
        bbox = _bbox_volume_a3(coords_arr, pad_a=bbox_pad_a)
        # Merge defaults first so they can't overwrite computed bbox volume.
        vol_meta = {**vol_meta, "source": "bbox", **bbox}

    volume_a3 = float(vol_meta.get("volume_a3") or 0.0)
    if not math.isfinite(volume_a3) or volume_a3 <= 0.0:
        return {"error": "Computed volume is not positive", "volume": vol_meta}

    # Density-style metrics.
    atom_density = float(total_atoms) / volume_a3
    mass_density_da_per_a3 = float(total_mass_da) / volume_a3
    packing_density = float(total_atomic_volume_a3) / volume_a3

    # Matthews-like ratio and a solvent estimate (heuristic, same constant as packing_metrics.py).
    matthews_a3_per_da = volume_a3 / float(total_mass_da) if total_mass_da > 0 else None
    protein_volume_a3 = float(total_mass_da) / 0.81 if total_mass_da > 0 else None
    solvent_volume_a3 = (volume_a3 - protein_volume_a3) if protein_volume_a3 is not None else None
    solvent_fraction = (float(solvent_volume_a3) / float(volume_a3)) if solvent_volume_a3 is not None else None
    solvent_percent = (
        float(max(0.0, min(100.0, float(solvent_fraction) * 100.0))) if solvent_fraction is not None else None
    )

    metric_descriptions = {
        "volume.volume_a3": "Volume used for lattice-wide metrics (Å³); from CRYST1 unit cell or a coordinate bounding box.",
        "lattice_atom_density_atoms_per_a3": "Atom count divided by the chosen volume (atoms/Å³).",
        "lattice_mass_density_da_per_a3": "Total mass (sum of atomic weights, Da) divided by the chosen volume (Da/Å³).",
        "lattice_packing_density_fraction": "Packing density as a fraction: Σ(atomic volumes)/volume (dimensionless).",
        "lattice_packing_density_percent": "Packing density expressed as a percentage (100×fraction).",
        "lattice_matthews_a3_per_da": "Matthews-like ratio: volume/mass (Å³/Da) for the whole supercell (heuristic).",
        "estimated_solvent_content_fraction": "Heuristic solvent fraction: (volume - mass/0.81)/volume (dimensionless).",
        "estimated_solvent_content_percent": "Heuristic solvent fraction expressed as a percentage (100×fraction).",
    }

    out = {
        "input": os.path.basename(pdb_path),
        "timestamp": datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
        "metric_descriptions": metric_descriptions,
        "chain_count": chain_count,
        "residue_count": residue_count,
        "total_atoms": total_atoms,
        "element_composition": dict(sorted(element_counts.items())),
        "total_mass_da": float(total_mass_da),
        "total_atomic_volume_a3": float(total_atomic_volume_a3),
        "volume": vol_meta,
        "lattice_atom_density_atoms_per_a3": atom_density,
        "lattice_mass_density_da_per_a3": mass_density_da_per_a3,
        "lattice_packing_density_fraction": packing_density,
        "lattice_packing_density_percent": packing_density * 100.0,
        "lattice_matthews_a3_per_da": matthews_a3_per_da,
        "estimated_protein_volume_a3": protein_volume_a3,
        "estimated_solvent_volume_a3": solvent_volume_a3,
        "estimated_solvent_content_fraction": solvent_fraction,
        "estimated_solvent_content_percent": solvent_percent,
    }

    # Small sanity note for space group.
    if cryst1 and str(cryst1.get("space_group", "")).strip() and str(cryst1.get("space_group", "")).strip() != "P1":
        out["warning"] = f"CRYST1 space_group is {cryst1.get('space_group')!r} (expected 'P1' for a supercell lattice fragment)"
    # Carry any volume-selection warning through to the top level for visibility.
    if isinstance(vol_meta, dict) and vol_meta.get("warning"):
        out["warning_volume_source"] = vol_meta.get("warning")

    return out


def main() -> None:
    ap = argparse.ArgumentParser(
        description="Analyse packing-style metrics for a symmetry-expanded lattice/supercell coordinate model.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples (from repository root):
  python metrics/lattice_packing_analyser.py supercell.pdb --volume-source bbox --bbox-pad 2.0
  python metrics/lattice_packing_analyser.py supercell.pdb --expected-chains 12 --volume-source bbox
  python metrics/lattice_packing_analyser.py supercells/ --sets set_a set_b --output-dir "lattice_pack_{}" --output-txt "lattice_{}.txt"
""",
    )
    ap.add_argument(
        "input",
        nargs="+",
        help="Input PDB file(s), directory, or glob pattern (e.g. supercell.pdb or '*.pdb').",
    )
    ap.add_argument(
        "--expected-chains",
        type=int,
        default=None,
        help="If set, require exactly this many chains in the structure (e.g. 12).",
    )
    ap.add_argument(
        "--require-cryst1",
        action="store_true",
        help="Fail if the PDB has no CRYST1 record.",
    )
    ap.add_argument(
        "--volume-source",
        choices=("auto", "cryst1", "bbox"),
        default="auto",
        help="Volume definition for lattice-wide density/packing metrics.",
    )
    ap.add_argument(
        "--bbox-pad",
        type=float,
        default=0.0,
        metavar="A",
        help="Padding (Å) applied to the coordinate bounding box when --volume-source=bbox.",
    )
    ap.add_argument(
        "--allow-cryst1-mismatch",
        action="store_true",
        help=(
            "Allow volume_source=cryst1 when expected_chains>1. Use only if CRYST1 describes the supercell "
            "volume (not the original unit cell)."
        ),
    )
    ap.add_argument(
        "--set",
        "-s",
        action="append",
        dest="sets",
        metavar="PATTERNS",
        help=(
            "Comma-separated patterns; file is included only if its basename contains ALL patterns. "
            "Repeat for multiple sets (e.g. -s set_a,set_b -s set_c,set_d)."
        ),
    )
    ap.add_argument(
        "--sets",
        dest="sets_multi",
        nargs="+",
        metavar="SET",
        help=(
            "Multiple set names in one go (one pattern per set). E.g. --sets set_a set_b set_c is "
            "equivalent to --set set_a --set set_b --set set_c."
        ),
    )
    ap.add_argument(
        "--output-dir",
        metavar="DIR",
        default=None,
        help=(
            "Output directory for per-structure JSON/TXT files. If it contains '{}', it is replaced "
            "by the set label. Required when analysing multiple structures unless --output-json uses '{}'."
        ),
    )
    ap.add_argument(
        "--output-json",
        "-o",
        metavar="FILE",
        help=(
            "Write JSON output. Single structure: file path (default: stdout). "
            "Multiple structures: use '{}' in the filename for per-structure output (replaced by file stem)."
        ),
    )
    ap.add_argument(
        "--output-txt",
        metavar="FILE",
        help=(
            "Write a human-readable text report. Single structure: file path. "
            "Multiple structures: use '{}' in the filename for per-structure output (replaced by file stem)."
        ),
    )
    ap.add_argument(
        "--dry-run",
        action="store_true",
        help="Print which files would be processed (per set) and exit without running.",
    )
    ap.add_argument(
        "--verbose",
        action="store_true",
        help="Write additional progress details to stderr (captured in --log).",
    )
    add_log_args(ap)
    args = ap.parse_args()
    log_setup = setup_log_from_args(
        args,
        script_path=__file__,
        inputs=list(getattr(args, "input", []) or []),
        pattern=None,
    )
    if log_setup is not None:
        log_setup.task("Lattice packing analysis (lattice_packing_analyser.py)")
        log_setup.kv("output_dir", getattr(args, "output_dir", None) or "(default)")
        log_setup.kv("output_json", getattr(args, "output_json", None) or "")
        log_setup.kv("output_txt", getattr(args, "output_txt", None) or "")
        log_setup.kv("dry_run", bool(getattr(args, "dry_run", False)))

    def _log(msg: str) -> None:
        if getattr(args, "verbose", False) or getattr(args, "log", None):
            print(msg, file=sys.stderr, flush=True)

    # Expand inputs to structure paths.
    paths = collect_structure_paths(list(getattr(args, "input", []) or []))
    if not paths:
        _log("[FoldKit] lattice_packing_analyser: no structure files found")
        if log_setup is not None:
            log_setup.error("No structure files found.")
        raise SystemExit(1)
    if log_setup is not None:
        log_setup.kv("n_inputs", len(paths))

    # Build list of (label, patterns, filtered_paths). If --set/--sets not used, one set with all paths.
    set_list: list[tuple[str, list[str], list[str]]] = []
    if getattr(args, "sets_multi", None):
        for name in args.sets_multi:
            patterns = [str(name).strip()] if str(name).strip() else []
            if patterns:
                label = patterns[0]
                filtered = filter_paths_by_patterns(paths, patterns)
                set_list.append((label, patterns, filtered))
    if getattr(args, "sets", None):
        for s in args.sets:
            patterns = [p.strip() for p in str(s).split(",") if p.strip()]
            if patterns:
                label = "_".join(patterns)
                filtered = filter_paths_by_patterns(paths, patterns)
                set_list.append((label, patterns, filtered))
    if not set_list:
        set_list = [("all", [], list(paths))]

    multi_any = any(len(filtered) > 1 for (_lab, _pat, filtered) in set_list) or len(paths) > 1

    # Validate output routing for multi-structure runs.
    out_dir_template = getattr(args, "output_dir", None)
    out_json_template = getattr(args, "output_json", None)
    out_txt_template = getattr(args, "output_txt", None)
    if multi_any:
        per_structure_json = bool(out_json_template and "{}" in str(out_json_template))
        per_structure_txt = bool(out_txt_template and "{}" in str(out_txt_template))
        combined_json = bool(out_json_template and "{}" not in str(out_json_template))
        combined_txt = bool(out_txt_template and "{}" not in str(out_txt_template))

        if not out_dir_template and not per_structure_json and not combined_json:
            _log(
                "[FoldKit] lattice_packing_analyser: "
                "multiple structures detected; set --output-dir or use --output-json with '{}' (per-structure) "
                "or a single JSON path (combined output)"
            )
            raise SystemExit(2)
        if out_txt_template and not per_structure_txt and not combined_txt and not out_dir_template:
            _log(
                "[FoldKit] lattice_packing_analyser: "
                "multiple structures detected; use --output-txt with '{}' (per-structure) or a single TXT path "
                "(combined output), or provide --output-dir"
            )
            raise SystemExit(2)

    if getattr(args, "dry_run", False):
        for label, patterns, filtered in set_list:
            if not filtered:
                _log(f"[FoldKit] set {label!r} (patterns={patterns}): 0 file(s)")
                continue
            _log(f"[FoldKit] set {label!r} (patterns={patterns}): {len(filtered)} file(s)")
            for p in filtered:
                _log(f"  - {os.path.basename(p)}")
        return

    def _fmt(x, ndp: int = 3) -> str:
        if x is None:
            return "N/A"
        try:
            xf = float(x)
        except Exception:
            return str(x)
        if not math.isfinite(xf):
            return "N/A"
        return f"{xf:.{ndp}f}"

    def _write_txt_report_to_stream(res0: dict, out_txt) -> None:
        if "error" in res0:
            out_txt.write("LATTICE PACKING ANALYSIS\n")
            out_txt.write(f"Error: {res0.get('error')}\n")
            return
        vol = (res0 or {}).get("volume", {}) if isinstance(res0, dict) else {}
        out_txt.write("LATTICE PACKING ANALYSIS\n")
        out_txt.write(f"Input: {res0.get('input')}\n")
        out_txt.write(
            f"Chains: {res0.get('chain_count')}  Residues: {res0.get('residue_count')}  Atoms: {res0.get('total_atoms')}\n"
        )
        out_txt.write(f"Volume source: {vol.get('source')}  Volume (Å³): {_fmt(vol.get('volume_a3'), 1)}\n")
        if vol.get("source") == "cryst1":
            cryst1 = vol.get("cryst1") or {}
            out_txt.write(
                "CRYST1: "
                f"a={_fmt(cryst1.get('a'), 3)} b={_fmt(cryst1.get('b'), 3)} c={_fmt(cryst1.get('c'), 3)}  "
                f"alpha={_fmt(cryst1.get('alpha'), 2)} beta={_fmt(cryst1.get('beta'), 2)} gamma={_fmt(cryst1.get('gamma'), 2)}  "
                f"space_group={cryst1.get('space_group', 'P1')}\n"
            )
        if vol.get("source") == "bbox":
            out_txt.write(f"BBox span (Å): {vol.get('span')}  pad (Å): {_fmt(vol.get('pad_a'), 2)}\n")
        out_txt.write("\n")
        out_txt.write("DENSITIES / PACKING\n")
        out_txt.write(f"Atom density (atoms/Å³): {_fmt(res0.get('lattice_atom_density_atoms_per_a3'), 6)}\n")
        out_txt.write(f"Mass density (Da/Å³): {_fmt(res0.get('lattice_mass_density_da_per_a3'), 6)}\n")
        out_txt.write(
            "Packing density: "
            f"{_fmt(res0.get('lattice_packing_density_fraction'), 6)} "
            f"({_fmt(res0.get('lattice_packing_density_percent'), 2)}%)\n"
        )
        out_txt.write("\n")
        out_txt.write("MATTHEWS-LIKE (heuristic)\n")
        out_txt.write(f"Matthews (Å³/Da): {_fmt(res0.get('lattice_matthews_a3_per_da'), 3)}\n")
        out_txt.write(f"Estimated protein volume (Å³): {_fmt(res0.get('estimated_protein_volume_a3'), 1)}\n")
        out_txt.write(f"Estimated solvent volume (Å³): {_fmt(res0.get('estimated_solvent_volume_a3'), 1)}\n")
        out_txt.write(
            "Estimated solvent content: "
            f"{_fmt(res0.get('estimated_solvent_content_fraction'), 4)} "
            f"({_fmt(res0.get('estimated_solvent_content_percent'), 1)}%)\n"
        )
        warn = res0.get("warning")
        if warn:
            out_txt.write(f"\nWarning: {warn}\n")

    def _write_txt_report(res0: dict, path_txt: str) -> None:
        with open(path_txt, "w", encoding="utf-8", newline="\n") as out_txt:
            _write_txt_report_to_stream(res0, out_txt)

    # Optional combined outputs (when templates are single paths without '{}').
    combined_json_path = str(out_json_template) if (multi_any and out_json_template and "{}" not in str(out_json_template)) else None
    combined_txt_path = str(out_txt_template) if (multi_any and out_txt_template and "{}" not in str(out_txt_template)) else None

    combined_records = []
    combined_txt_f = None
    if combined_txt_path:
        combined_txt_f = open(combined_txt_path, "w", encoding="utf-8", newline="\n")
        combined_txt_f.write("FOLDKIT LATTICE PACKING (COMBINED REPORT)\n\n")

    # Run sets / structures.
    had_error = False
    for label, patterns, filtered in set_list:
        if not filtered:
            _log(f"[FoldKit] set {label!r}: no files match (patterns={patterns})")
            continue

        out_dir = None
        if out_dir_template:
            out_dir = str(out_dir_template)
            if "{}" in out_dir:
                out_dir = out_dir.replace("{}", label)
            os.makedirs(out_dir, exist_ok=True)

        for p in filtered:
            stem = os.path.splitext(os.path.basename(p))[0]
            _log(f"[FoldKit] analysing {os.path.basename(p)} (set={label!r})")
            res0 = analyse_lattice_packing(
                p,
                expected_chains=getattr(args, "expected_chains", None),
                require_cryst1=bool(getattr(args, "require_cryst1", False)),
                volume_source=str(getattr(args, "volume_source", "auto")),
                bbox_pad_a=float(getattr(args, "bbox_pad", 0.0)),
                allow_cryst1_mismatch=bool(getattr(args, "allow_cryst1_mismatch", False)),
            )
            if isinstance(res0, dict) and "error" in res0:
                had_error = True
                _log(f"[FoldKit] error: {res0.get('error')}")

            # JSON output routing.
            if combined_json_path is not None:
                combined_records.append(
                    {
                        "set_label": label,
                        "structure_basename": os.path.basename(p),
                        "structure_stem": stem,
                        "result": res0,
                    }
                )
            elif out_json_template and "{}" in str(out_json_template):
                out_json_path = str(out_json_template).replace("{}", stem)
                with open(out_json_path, "w", encoding="utf-8", newline="\n") as f:
                    json.dump(res0, f, indent=2, sort_keys=True)
                    f.write("\n")
            elif out_dir:
                out_json_path = os.path.join(out_dir, f"{stem}_lattice_packing.json")
                with open(out_json_path, "w", encoding="utf-8", newline="\n") as f:
                    json.dump(res0, f, indent=2, sort_keys=True)
                    f.write("\n")
            else:
                # Single-structure: stdout by default, but avoid duplicating into the log.
                # SummaryLog (--log) does not tee stdout; write JSON to stdout unless --output-json is set.
                out_json = sys.stdout
                if out_json_template:
                    out_json = open(str(out_json_template), "w", encoding="utf-8", newline="\n")
                try:
                    json.dump(res0, out_json, indent=2, sort_keys=True)
                    out_json.write("\n")
                finally:
                    if out_json_template and out_json is not sys.stdout:
                        out_json.close()

            # TXT output routing (optional).
            if combined_txt_f is not None:
                combined_txt_f.write("=" * 72 + "\n")
                combined_txt_f.write(f"{os.path.basename(p)}\n")
                combined_txt_f.write("=" * 72 + "\n")
                _write_txt_report_to_stream(res0, combined_txt_f)
                combined_txt_f.write("\n\n")
            elif out_txt_template and "{}" in str(out_txt_template):
                out_txt_path = str(out_txt_template).replace("{}", stem)
                _write_txt_report(res0, out_txt_path)
            elif out_dir and out_txt_template:
                out_txt_path = os.path.join(out_dir, f"{stem}_lattice_packing.txt")
                _write_txt_report(res0, out_txt_path)
            elif out_txt_template and not multi_any:
                _write_txt_report(res0, str(out_txt_template))

    if combined_txt_f is not None:
        combined_txt_f.close()

    if combined_json_path is not None:
        combined = {
            "n_structures": len(combined_records),
            "structures": combined_records,
        }
        with open(combined_json_path, "w", encoding="utf-8", newline="\n") as f:
            json.dump(combined, f, indent=2, sort_keys=True)
            f.write("\n")

    if had_error:
        raise SystemExit(1)


if __name__ == "__main__":
    main()

