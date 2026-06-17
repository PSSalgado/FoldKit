"""
Coot SSM script templates (Python strings passed to ``coot --script``).

Used by superimpose_coot_SSM.py and SSM_aligned_core.py.
"""

from __future__ import annotations

from superimposition.aligned_pdb_names import COOT_ALIGNED_PDB_HELPERS_PY

_POST_BATCH = """
print("FOLDKIT_ALIGNMENTS_DONE: All superpositions written to disk. Exiting Coot.")
print("Aligned structures saved to: {output_dir}")
coot_real_exit(0)
"""

# Shared helpers: chain detection and SSM with CA-only fallback (matches Step 1 script).
COOT_SSM_SUPERPOSE_PY = """
def find_available_chains(mol_id):
    chains = []
    try:
        chain_list = chain_ids(mol_id)
        if chain_list:
            chains = chain_list
        else:
            for chain_id in ["A", "B", "C", "D", "E", "F", "G", "H", "I", "J"]:
                try:
                    if n_atoms_in_chain(mol_id, chain_id) > 0:
                        chains.append(chain_id)
                except:
                    continue
    except Exception as e:
        print("Error getting chains for molecule " + str(mol_id) + ": " + str(e))
        chains = ["A"]
    return chains

def validate_chain_exists(mol_id, chain_id):
    try:
        available_chains = chain_ids(mol_id)
        if available_chains and chain_id in available_chains:
            return True
        try:
            return chain_length(mol_id, chain_id) > 0
        except:
            pass
        return False
    except:
        return False

def resolve_chain(mol_id, preferred):
    if validate_chain_exists(mol_id, preferred):
        return preferred
    chains = find_available_chains(mol_id)
    if chains:
        print("Using chain '" + chains[0] + "' instead of '" + preferred + "'")
        return chains[0]
    print("No chains found - using 'A' as fallback")
    return "A"

def ensure_secondary_structure(mol_id):
    \"\"\"SSM graph building needs secondary structure; assign from geometry if missing.\"\"\"
    try:
        ss = get_header_secondary_structure_info(mol_id)
        has_ss = isinstance(ss, dict) and (ss.get('helices') or ss.get('strands'))
        if not has_ss:
            add_header_secondary_structure_info(mol_id)
    except Exception as e:
        print("Warning: could not assign secondary structure for molecule " + str(mol_id) + ": " + str(e))

def foldkit_ssm_has_alignment(reference_mol, model_mol, ref_chain_id, model_chain_id):
    \"\"\"Return True when post-SSM CA pairs exist; defer to host log audit if API missing.\"\"\"
    ref_ca = "//" + ref_chain_id + "//CA"
    model_ca = "//" + model_chain_id + "//CA"
    try:
        pairs = match_atoms(reference_mol, model_mol, ref_ca, model_ca, 5.0)
        if pairs is None:
            return False
        n_pairs = len(pairs) if hasattr(pairs, "__len__") else int(pairs)
        if n_pairs <= 0:
            print("SSM alignment check: 0 matched CA pairs after superposition")
            return False
        return True
    except Exception as e:
        print("SSM alignment check skipped (API unavailable): " + str(e))
        return True

def foldkit_ssm_superpose(reference_mol, model_mol, ref_chain_id, model_chain_id, model_name):
    ref_chain_id = resolve_chain(reference_mol, ref_chain_id)
    model_chain_id = resolve_chain(model_mol, model_chain_id)
    ensure_secondary_structure(reference_mol)
    ensure_secondary_structure(model_mol)
    ref_selection = "//" + ref_chain_id + "//"
    model_selection = "//" + model_chain_id + "//"
    print("Attempting SSM: " + model_name + " (chain " + model_chain_id + ") onto reference (chain " + ref_chain_id + ")")
    print("Reference selection: " + ref_selection)
    print("Model selection: " + model_selection)
    try:
        superpose_with_atom_selection(reference_mol, model_mol, ref_selection, model_selection, 0)
        if not foldkit_ssm_has_alignment(reference_mol, model_mol, ref_chain_id, model_chain_id):
            print("SSM superposition produced no verifiable alignment for " + model_name)
            return False
        print("SSM superposition successful for " + model_name)
        return True
    except Exception as e:
        print("Error during SSM superposition of " + model_name + ": " + str(e))
    try:
        print("Trying SSM with CA atoms only for " + model_name + "...")
        ref_ca = "//" + ref_chain_id + "//CA"
        model_ca = "//" + model_chain_id + "//CA"
        superpose_with_atom_selection(reference_mol, model_mol, ref_ca, model_ca, 0)
        if not foldkit_ssm_has_alignment(reference_mol, model_mol, ref_chain_id, model_chain_id):
            print("CA-only SSM produced no verifiable alignment for " + model_name)
            return False
        print("CA-only SSM superposition successful for " + model_name)
        return True
    except Exception as e2:
        print("SSM superposition failed for " + model_name + ": " + str(e2))
        return False
"""


def create_one_to_many_ssm_script(
    reference_file: str,
    model_files: list[str],
    output_dir: str,
    *,
    keep_coot_open: bool = True,
    aligned_tag: str = "_SSMaligned2_",
) -> str:
    """One-to-many SSM: each model superimposed onto a fixed reference."""
    script_content = (
        COOT_ALIGNED_PDB_HELPERS_PY
        + COOT_SSM_SUPERPOSE_PY
        + """
import os
import sys

set_nomenclature_errors_on_read("ignore")
set_show_symmetry_master(0)

reference_mol = read_pdb("{reference_file}")
ref_name = os.path.splitext(os.path.basename("{reference_file}"))[0]
reference_chain = resolve_chain(reference_mol, chain_ids(reference_mol)[0])
graphics_to_ca_representation(reference_mol)

if not os.path.exists("{output_dir}"):
    os.makedirs("{output_dir}")

for model_path in {model_files}:
    model_mol = read_pdb(model_path)
    model_name = os.path.splitext(os.path.basename(model_path))[0]
    model_chain = resolve_chain(model_mol, chain_ids(model_mol)[0])
    print("Superposing " + model_name + " onto " + ref_name)
    graphics_to_ca_representation(model_mol)

    if foldkit_ssm_superpose(reference_mol, model_mol, reference_chain, model_chain, model_name):
        output_name = os.path.join("{output_dir}", foldkit_aligned_pdb_basename(model_name, ref_name, "__ALIGNED_TAG__"))
        write_pdb_file(model_mol, output_name)
        print("Saved aligned model: " + output_name)
    else:
        print("Skipping write for failed superposition: " + model_name)
    close_molecule(model_mol)

close_molecule(reference_mol)

__POST_CLOSE__
"""
    )
    post_open = """
print("FOLDKIT_ALIGNMENTS_DONE: All superpositions written to disk. Reloading structures in Coot for inspection...")
handle_read_draw_molecule_with_recentre("{reference_file}", 0)
ref_mol = graphics_n_molecules() - 1
set_molecule_bonds_colour_map_rotation(ref_mol, 0)
graphics_to_ca_representation(int(ref_mol))

for model_path in {model_files}:
    model_name = os.path.splitext(os.path.basename(model_path))[0]
    aligned_path = os.path.join("{output_dir}", foldkit_aligned_pdb_basename(model_name, ref_name, "__ALIGNED_TAG__"))
    if os.path.exists(aligned_path):
        handle_read_draw_molecule_with_recentre(aligned_path, 0)
        mol = graphics_n_molecules() - 1
        set_molecule_bonds_colour_map_rotation(mol, 21 * (mol - ref_mol))
        graphics_to_ca_representation(int(mol))

x, y, z = molecule_centre(ref_mol)
set_rotation_centre(x, y, z)
"""
    post_batch = _POST_BATCH.format(output_dir=output_dir)
    post_close = post_open if keep_coot_open else post_batch
    tpl = script_content.replace("__POST_CLOSE__", post_close)
    tpl = tpl.replace("__ALIGNED_TAG__", aligned_tag)
    return tpl.format(
        reference_file=reference_file,
        output_dir=output_dir,
        model_files=model_files,
    )


# Alias used by superimpose_coot_SSM.py
create_coot_script = create_one_to_many_ssm_script


def create_all_vs_all_ssm_script(
    model_files: list[str],
    output_dir: str,
    *,
    keep_coot_open: bool = True,
) -> str:
    """All-vs-all SSM superposition."""
    script_content = (
        COOT_ALIGNED_PDB_HELPERS_PY
        + COOT_SSM_SUPERPOSE_PY
        + """
import os
import sys

set_nomenclature_errors_on_read("ignore")
set_show_symmetry_master(0)

if not os.path.exists("{output_dir}"):
    os.makedirs("{output_dir}")

model_files = {model_files}

for ref_path in model_files:
    reference_mol = read_pdb(ref_path)
    reference_chain = resolve_chain(reference_mol, chain_ids(reference_mol)[0])
    ref_name = os.path.splitext(os.path.basename(ref_path))[0]
    graphics_to_ca_representation(reference_mol)
    print("Superposing onto reference: " + ref_name)

    for model_path in model_files:
        if model_path == ref_path:
            continue
        model_name = os.path.splitext(os.path.basename(model_path))[0]
        print("  Superposing " + model_name + " onto " + ref_name)
        model_mol = read_pdb(model_path)
        model_chain = resolve_chain(model_mol, chain_ids(model_mol)[0])
        graphics_to_ca_representation(model_mol)
        if foldkit_ssm_superpose(reference_mol, model_mol, reference_chain, model_chain, model_name):
            output_name = os.path.join("{output_dir}", foldkit_aligned_pdb_basename(model_name, ref_name, "_SSMaligned2_"))
            write_pdb_file(model_mol, output_name)
        else:
            print("Skipping write for failed superposition: " + model_name + " onto " + ref_name)
        close_molecule(model_mol)
    close_molecule(reference_mol)

__ALL_VS_ALL_TAIL__
"""
    )
    tail_open = """
print("FOLDKIT_ALIGNMENTS_DONE: All pairwise superpositions written to disk. Reloading structures in Coot for inspection...")
unique_files = list(set(model_files))
for i, file_path in enumerate(unique_files):
    mol = read_pdb(file_path)
    graphics_to_ca_representation(mol)
    set_molecule_bonds_colour_map_rotation(mol, 20 * i)
print("Aligned structures saved to: {output_dir}")
"""
    tail_batch = """
print("FOLDKIT_ALIGNMENTS_DONE: All pairwise superpositions written to disk. Exiting Coot.")
for i in range(graphics_n_molecules()):
    close_molecule(i)
print("All-vs-all SSM superposition complete (batch). Aligned structures saved to: {output_dir}")
coot_real_exit(0)
"""
    tail = tail_open if keep_coot_open else tail_batch
    script_content = script_content.replace("__ALL_VS_ALL_TAIL__", tail)
    return script_content.format(output_dir=output_dir, model_files=model_files)


def create_full_to_core_ssm_script(
    pairs: list[tuple[str, str, str]],
    output_dir: str,
    *,
    keep_coot_open: bool = False,
    aligned_tag: str = "_SSMfull2core_",
) -> str:
    """
    Superimpose full native structures onto their aligned core references.

    Each entry in ``pairs`` is (reference_core_path, moving_full_path, label_stem).
    """
    script_content = (
        COOT_ALIGNED_PDB_HELPERS_PY
        + COOT_SSM_SUPERPOSE_PY
        + """
import os
import sys

set_nomenclature_errors_on_read("ignore")
set_show_symmetry_master(0)

if not os.path.exists("{output_dir}"):
    os.makedirs("{output_dir}")

pairs = {pairs}
aligned_tag = "{aligned_tag}"

for ref_path, mov_path, label in pairs:
    reference_mol = read_pdb(ref_path)
    reference_chain = resolve_chain(reference_mol, chain_ids(reference_mol)[0])
    ref_name = os.path.splitext(os.path.basename(ref_path))[0]
    model_name = label
    print("Superposing full " + model_name + " onto core " + ref_name)

    model_mol = read_pdb(mov_path)
    model_chain = resolve_chain(model_mol, chain_ids(model_mol)[0])
    graphics_to_ca_representation(reference_mol)
    graphics_to_ca_representation(model_mol)

    if foldkit_ssm_superpose(reference_mol, model_mol, reference_chain, model_chain, model_name):
        output_name = os.path.join("{output_dir}", foldkit_aligned_pdb_basename(model_name, ref_name, aligned_tag))
        write_pdb_file(model_mol, output_name)
        print("Saved aligned model: " + output_name)
    else:
        print("Skipping write for failed superposition: full " + model_name)
    close_molecule(reference_mol)
    close_molecule(model_mol)

__POST_CLOSE__
"""
    )
    post_batch = _POST_BATCH.format(output_dir=output_dir)
    post_close = post_batch if not keep_coot_open else """
print("FOLDKIT_ALIGNMENTS_DONE: Full-to-core superpositions written to disk.")
"""
    tpl = script_content.replace("__POST_CLOSE__", post_close)
    return tpl.format(
        output_dir=output_dir,
        pairs=pairs,
        aligned_tag=aligned_tag,
    )
