"""
Extract core residues from PDB/mmCIF by ResidueKey and write renumbered PDBs.
"""

from __future__ import annotations

import copy
import os
import warnings
from typing import Iterable

from utils.ssm_log_parse import ResidueKey

try:
    from Bio.PDB.MMCIFParser import MMCIFParser
    from Bio.PDB.PDBIO import PDBIO
    from Bio.PDB.PDBParser import PDBParser
    from Bio.PDB.PDBExceptions import PDBConstructionWarning
    from Bio.PDB.Structure import Structure as StructureClass

    warnings.simplefilter("ignore", PDBConstructionWarning)
    BIOPYTHON_AVAILABLE = True
except ImportError:
    BIOPYTHON_AVAILABLE = False
    StructureClass = None  # type: ignore[misc, assignment]


def _normalize_icode(icode: str) -> str:
    s = (icode or "").strip()
    return s[0] if s else " "


def residue_key_from_bio(residue) -> ResidueKey | None:
    if residue.id[0] != " ":
        return None
    icode = residue.id[2] if len(residue.id) > 2 else " "
    return ResidueKey(
        chain=residue.get_parent().id,
        resnum=int(residue.id[1]),
        icode=_normalize_icode(icode),
    )


def _key_tuple(key: ResidueKey) -> tuple[str, int, str]:
    ic = key.icode if key.icode else " "
    return (key.chain, key.resnum, ic)


def load_structure(path: str):
    if not BIOPYTHON_AVAILABLE:
        raise RuntimeError("BioPython is required (pip install biopython)")
    path = os.path.abspath(os.path.expanduser(path))
    label = os.path.splitext(os.path.basename(path))[0]
    if path.lower().endswith(".cif"):
        parser = MMCIFParser(QUIET=True)
    else:
        parser = PDBParser(QUIET=True)
    return parser.get_structure(label, path)


def find_structure_file(directory: str, stem: str) -> str | None:
    """Return path to ``stem.pdb`` or ``stem.cif`` under directory."""
    directory = os.path.abspath(os.path.expanduser(directory))
    for ext in (".pdb", ".cif", ".mmcif"):
        candidate = os.path.join(directory, stem + ext)
        if os.path.isfile(candidate):
            return candidate
    return None


def find_native_structure(
    models_dir: str,
    label: str,
    *,
    name_map: dict[str, str] | None = None,
) -> str | None:
    """Resolve native structure path, optionally via a short-stem name map."""
    if name_map and label in name_map:
        original = name_map[label].strip()
        if original:
            stem = os.path.splitext(os.path.basename(original))[0]
            found = find_structure_file(models_dir, stem)
            if found:
                return found
    return find_structure_file(models_dir, label)


def count_core_ca_present(
    structure_path: str,
    core_keys: Iterable[ResidueKey],
) -> int:
    """Return how many core ResidueKeys have a Cα atom in the structure."""
    if not os.path.isfile(structure_path):
        return 0
    structure = load_structure(structure_path)
    model = structure[0]
    wanted = {_key_tuple(k) for k in core_keys}
    found = 0
    for chain in model:
        for residue in chain:
            key = residue_key_from_bio(residue)
            if key is None:
                continue
            if _key_tuple(key) not in wanted:
                continue
            if "CA" in residue:
                found += 1
    return found


def extract_core_structure(
    structure_path: str,
    core_keys: list[ResidueKey],
    *,
    out_chain: str = "A",
    renumber_from: int = 1,
):
    """
    Build a one-chain structure containing core residues in ``core_keys`` order.

    Residues are renumbered sequentially from ``renumber_from``.
    """
    structure = load_structure(structure_path)
    model = structure[0]

    key_to_residue = {}
    for chain in model:
        for residue in chain:
            key = residue_key_from_bio(residue)
            if key is None:
                continue
            key_to_residue[_key_tuple(key)] = residue

    from Bio.PDB.Model import Model
    from Bio.PDB.Chain import Chain

    out_struct = StructureClass("core")
    out_model = Model(0)
    out_struct.add(out_model)
    out_chain_obj = Chain(out_chain)
    out_model.add(out_chain_obj)

    resnum = renumber_from
    for core_key in core_keys:
        kt = _key_tuple(core_key)
        if kt not in key_to_residue:
            raise KeyError(
                "Residue {} not found in {}".format(core_key, structure_path)
            )
        src = key_to_residue[kt]
        new_res = src.copy()
        new_res.id = (" ", resnum, " ")
        out_chain_obj.add(new_res)
        resnum += 1

    return out_struct


def extract_ca_coords(structure) -> list[tuple[float, float, float]]:
    """Return Cα coordinates in chain/residue order (one chain per trimmed model)."""
    model = structure[0]
    coords: list[tuple[float, float, float]] = []
    for chain in model:
        for residue in chain:
            if residue.id[0] != " ":
                continue
            if "CA" not in residue:
                continue
            atom = residue["CA"]
            vec = atom.get_vector().get_array()
            coords.append((float(vec[0]), float(vec[1]), float(vec[2])))
    return coords


def _ca_rmsd_angstrom(
    coords_a: list[tuple[float, float, float]],
    coords_b: list[tuple[float, float, float]],
) -> float:
    """RMSD between paired Cα coordinates (same length, no superposition)."""
    if len(coords_a) != len(coords_b):
        raise ValueError(
            "Cα count mismatch: {} vs {}".format(len(coords_a), len(coords_b))
        )
    n = len(coords_a)
    if n == 0:
        return float("nan")
    sse = 0.0
    for (x1, y1, z1), (x2, y2, z2) in zip(coords_a, coords_b):
        dx = x1 - x2
        dy = y1 - y2
        dz = z1 - z2
        sse += dx * dx + dy * dy + dz * dz
    return (sse / n) ** 0.5


def pairwise_core_ca_rmsd(
    labelled_structures: list[tuple[str, object]],
) -> dict[tuple[str, str], float]:
    """
    Pairwise Cα RMSD on renumbered core residues.

    Structures must already share the anchor coordinate frame (same core column
    order and count per model).
    """
    from utils.ssm_log_parse import _natural_sort_key

    labels = sorted((lab for lab, _ in labelled_structures), key=_natural_sort_key)
    struct_map = {lab: struct for lab, struct in labelled_structures}
    coords_map = {lab: extract_ca_coords(struct_map[lab]) for lab in labels}
    n_core = len(next(iter(coords_map.values()), []))

    results: dict[tuple[str, str], float] = {}
    for i, la in enumerate(labels):
        for lb in labels[i:]:
            rmsd = _ca_rmsd_angstrom(coords_map[la], coords_map[lb])
            results[(la, lb)] = round(rmsd, 4)
    if labels and n_core:
        for lab in labels:
            results[(lab, lab)] = 0.0
    return results


def write_structure_pdb(structure, path: str) -> None:
    path = os.path.abspath(os.path.expanduser(path))
    os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
    io = PDBIO()
    io.set_structure(structure)
    io.save(path)


def merge_structures_to_pdb(
    labelled_structures: list[tuple[str, object]],
    path: str,
) -> None:
    """Write one PDB with each input structure on its own chain (labels → chain IDs)."""
    if not BIOPYTHON_AVAILABLE:
        raise RuntimeError("BioPython is required (pip install biopython)")

    from Bio.PDB.Model import Model
    from Bio.PDB.Chain import Chain

    combined = StructureClass("core_superimposed")
    model = Model(0)
    combined.add(model)

    for idx, (_label, struct) in enumerate(labelled_structures):
        chain_id = _chain_id_for_index(idx)
        src_model = struct[0]
        out_chain = Chain(chain_id)
        model.add(out_chain)
        for src_chain in src_model:
            for residue in src_chain:
                out_chain.add(residue.copy())

    write_structure_pdb(combined, path)


def _chain_id_for_index(index: int) -> str:
    """Map 0-based index to PDB chain ID (A–Z, then AA, AB, …)."""
    alphabet = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
    if index < 26:
        return alphabet[index]
    first = alphabet[(index // 26) - 1]
    second = alphabet[index % 26]
    return first + second


MIN_EQUIVALENCE_CA_PAIRS = 3
TRIM_RENUMBER_FROM = 1


def find_residue_by_key(struct, key: ResidueKey):
    """Return protein residue matching ``key`` in a BioPython structure."""
    ic = key.icode if key.icode else " "
    model = struct[0]
    for chain in model:
        if chain.id != key.chain:
            continue
        for residue in chain:
            if residue.id[0] != " ":
                continue
            ric = residue.id[2] if len(residue.id) > 2 else " "
            if residue.id[1] == key.resnum and ric == ic:
                return residue
    return None


def get_trim_residue_by_index(
    trim_struct,
    index: int,
    *,
    renumber_from: int = TRIM_RENUMBER_FROM,
):
    """Residue at ``index`` in a single-chain renumbered trim structure."""
    chain = next(trim_struct[0].get_chains())
    resnum = renumber_from + index
    rid = (" ", resnum, " ")
    if rid not in chain:
        raise KeyError(
            "trim residue index {} (seq {}) missing".format(index, resnum)
        )
    return chain[rid]


def trim_key_index(keys: list[ResidueKey], key: ResidueKey) -> int:
    kt = _key_tuple(key)
    for i, k in enumerate(keys):
        if _key_tuple(k) == kt:
            return i
    raise KeyError("{} not among {} trim keys".format(key, len(keys)))


def equivalence_ca_atom_pairs_from_trims(
    columns: list[dict[str, ResidueKey]],
    ref_label: str,
    mov_label: str,
    ref_trim_struct,
    mov_trim_struct,
    ref_keys: list[ResidueKey],
    mov_keys: list[ResidueKey],
) -> tuple[list, list]:
    """Paired Cα atoms for equivalence columns between two trim models."""
    ref_atoms: list = []
    mov_atoms: list = []
    for col in columns:
        if ref_label not in col or mov_label not in col:
            continue
        try:
            ri = trim_key_index(ref_keys, col[ref_label])
            mi = trim_key_index(mov_keys, col[mov_label])
        except KeyError:
            continue
        ref_res = get_trim_residue_by_index(ref_trim_struct, ri)
        mov_res = get_trim_residue_by_index(mov_trim_struct, mi)
        if "CA" not in ref_res or "CA" not in mov_res:
            continue
        ref_atoms.append(ref_res["CA"])
        mov_atoms.append(mov_res["CA"])
    return ref_atoms, mov_atoms


def equivalence_ca_atom_pairs_native_to_trim(
    trim_keys: list[ResidueKey],
    ref_trim_struct,
    mov_native_struct,
) -> tuple[list, list]:
    """Paired Cα atoms between a trim reference and a full native structure."""
    ref_atoms: list = []
    mov_atoms: list = []
    for i, key in enumerate(trim_keys):
        ref_res = get_trim_residue_by_index(ref_trim_struct, i)
        mov_res = find_residue_by_key(mov_native_struct, key)
        if mov_res is None or "CA" not in ref_res or "CA" not in mov_res:
            continue
        ref_atoms.append(ref_res["CA"])
        mov_atoms.append(mov_res["CA"])
    return ref_atoms, mov_atoms


def superimpose_moving_onto_reference(
    ref_atoms: list,
    mov_atoms: list,
    moving_struct,
) -> float:
    """Apply a Kabsch superposition; return post-fit Cα RMSD (Å)."""
    if not BIOPYTHON_AVAILABLE:
        raise RuntimeError("BioPython is required (pip install biopython)")
    from Bio.PDB.Superimposer import Superimposer

    if len(ref_atoms) < MIN_EQUIVALENCE_CA_PAIRS:
        raise ValueError(
            "need at least {} Cα pairs, got {}".format(
                MIN_EQUIVALENCE_CA_PAIRS, len(ref_atoms)
            )
        )
    sup = Superimposer()
    sup.set_atoms(ref_atoms, mov_atoms)
    sup.apply(list(moving_struct.get_atoms()))
    return float(sup.rms)


def equivalence_identity_stats(
    columns: list[dict[str, ResidueKey]],
    ref_label: str,
    mov_label: str,
    ref_native_path: str,
    mov_native_path: str,
) -> tuple[int, float | None]:
    """Return (n_shared_columns, pct_identity) from native residue names."""
    ref_native = load_structure(ref_native_path)
    mov_native = load_structure(mov_native_path)
    n = 0
    ident = 0
    for col in columns:
        if ref_label not in col or mov_label not in col:
            continue
        ref_res = find_residue_by_key(ref_native, col[ref_label])
        mov_res = find_residue_by_key(mov_native, col[mov_label])
        if ref_res is None or mov_res is None:
            continue
        n += 1
        if ref_res.get_resname().strip() == mov_res.get_resname().strip():
            ident += 1
    pct = round(100.0 * ident / n, 2) if n else None
    return n, pct


def ssm_block_ca_atom_pairs(
    block,
    ref_native_struct,
    mov_native_struct,
) -> tuple[list, list]:
    """Paired Cα atoms from one Coot SSM superposition block (native structures)."""
    ref_atoms: list = []
    mov_atoms: list = []
    for mov_key, ref_key, _dist in block.equivalences:
        ref_res = find_residue_by_key(ref_native_struct, ref_key)
        mov_res = find_residue_by_key(mov_native_struct, mov_key)
        if ref_res is None or mov_res is None:
            continue
        if "CA" not in ref_res or "CA" not in mov_res:
            continue
        ref_atoms.append(ref_res["CA"])
        mov_atoms.append(mov_res["CA"])
    return ref_atoms, mov_atoms


def ssm_block_identity_stats(
    block,
    ref_native_struct,
    mov_native_struct,
) -> tuple[int, float | None]:
    """Return (n_pairs, pct_identity) from native residue names in an SSM block."""
    n = 0
    ident = 0
    for mov_key, ref_key, _dist in block.equivalences:
        ref_res = find_residue_by_key(ref_native_struct, ref_key)
        mov_res = find_residue_by_key(mov_native_struct, mov_key)
        if ref_res is None or mov_res is None:
            continue
        n += 1
        if ref_res.get_resname().strip() == mov_res.get_resname().strip():
            ident += 1
    pct = round(100.0 * ident / n, 2) if n else None
    return n, pct


def ordered_protein_residue_keys(
    structure_path: str,
    *,
    chain: str | None = None,
) -> list[ResidueKey]:
    """Residue keys in PDB chain order (standard amino acids only)."""
    structure = load_structure(structure_path)
    keys: list[ResidueKey] = []
    for ch in structure[0]:
        if chain is not None and ch.id != chain:
            continue
        for residue in ch:
            key = residue_key_from_bio(residue)
            if key is not None:
                keys.append(key)
    return keys


def _residue_keys_consecutive_in_chain(a: ResidueKey, b: ResidueKey) -> bool:
    if a.chain != b.chain:
        return False
    ia = a.icode if a.icode else " "
    ib = b.icode if b.icode else " "
    if ia != ib and ia != " " and ib != " ":
        return False
    return b.resnum == a.resnum + 1


def longest_contiguous_residue_keys(
    structure_path: str,
    keys: list[ResidueKey],
) -> list[ResidueKey]:
    """
    Longest contiguous subsequence of ``keys`` in native PDB chain order.

    When SSM/MSA matched residues are sparse or gapped in sequence, this keeps
    only the longest run of consecutive residue numbers on the original chain.
    """
    if not keys:
        return []
    selected = {_key_tuple(k) for k in keys}
    chains = sorted({k.chain for k in keys})
    in_order: list[ResidueKey] = []
    for chain_id in chains:
        for key in ordered_protein_residue_keys(structure_path, chain=chain_id):
            if _key_tuple(key) in selected:
                in_order.append(key)
    if not in_order:
        return []

    best_start = 0
    best_len = 1
    cur_start = 0
    cur_len = 1
    for i in range(1, len(in_order)):
        if _residue_keys_consecutive_in_chain(in_order[i - 1], in_order[i]):
            cur_len += 1
        else:
            if cur_len > best_len:
                best_len = cur_len
                best_start = cur_start
            cur_start = i
            cur_len = 1
    if cur_len > best_len:
        best_len = cur_len
        best_start = cur_start
    return in_order[best_start : best_start + best_len]
