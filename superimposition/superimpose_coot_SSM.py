import hashlib
import os
import sys
from subprocess import Popen, PIPE, STDOUT
import glob
import re
import fnmatch

_REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if _REPO_ROOT not in sys.path:
    sys.path.insert(0, _REPO_ROOT)

from superimposition.superimpose_pattern_match import find_ref_model_matches

from cli_log import setup_log_from_argv

_argv_no_log, _ = setup_log_from_argv(
    script_path=__file__,
    argv=sys.argv[1:],
    inputs=[sys.argv[1]] if len(sys.argv) > 1 else [],
    pattern=None,
)
sys.argv = [sys.argv[0]] + _argv_no_log

if _ is not None:
    _.task("Coot SSM superposition (superimpose_coot_SSM.py)")
    try:
        _.kv("argv", _argv_no_log)
    except Exception:
        pass

def sanitize_pattern_for_filename(pattern):
    """Make a filter pattern safe for use in log/rmsd filenames."""
    if not pattern:
        return ""
    s = re.sub(r'[\*?\ /\\]', '_', pattern)
    return s[:64]  # limit length

def expand_output_dir_pattern(template, ref_name=None, filter_pattern=None):
    """
    Expand --output-dir pattern: [reference_name] -> ref stem, [pattern] or [filter] -> sanitised filter.
    Also replaces literal *filter* with the sanitised filter (convenience placeholder).
    """
    if not template:
        return template
    s = template
    ref = ref_name if ref_name is not None else ""
    pat = sanitize_pattern_for_filename(filter_pattern) if filter_pattern else ""
    s = s.replace("[reference_name]", ref).replace("[pattern]", pat).replace("[filter]", pat)
    if pat:
        s = s.replace("*filter*", pat)
    return s


def _announce_superposition_finished(
    log_file, exit_code, summary="", rmsd_format=None
):
    """
    Announce completion on stderr and append the same text to the Coot log.

    While Coot runs, its stdout is only copied to the log file (not echoed), so the
    terminal can stay quiet. Using stderr here avoids losing the message when Python
    stdout is fully buffered or redirected (non-TTY, some IDEs, wrappers).
    """
    lines = [
        "",
        "=" * 72,
    ]
    if exit_code == 0:
        lines.append(
            "Finished: all superpositions for this run are complete (Coot exit code 0)."
        )
    else:
        lines.append(
            "Coot process ended (exit code {}). Check the log if outputs look wrong.".format(
                exit_code
            )
        )
    if summary:
        lines.append(summary)
    lines.append("Log: {}".format(log_file))
    lines.append("=" * 72)
    if rmsd_format:
        lines.append("")
        lines.append("To extract RMSD values, run:")
        lines.append(
            "python ranking/extract_rmsd.py --format {} {}".format(rmsd_format, log_file)
        )
    text = "\n".join(lines) + "\n"
    print(text, file=sys.stderr, end="", flush=True)
    try:
        with open(log_file, "a", encoding="utf-8", errors="replace") as logf:
            logf.write(text)
    except OSError:
        pass


# Coot prints a line containing this marker when alignments are saved (before reload or exit).
# SSM: Coot User Manual "Secondary Structure Matching (SSM)" — superpose_with_atom_selection
# (last argument is move_copy_flag: 0 = transform moving molecule in place before write_pdb_file).
COOT_STDOUT_ALIGNMENTS_DONE_MARKER = "FOLDKIT_ALIGNMENTS_DONE"


def _pipe_coot_stdout_line(line, log_f, echo_all_stdout=False):
    """Append Coot stdout to the log; echo marker lines to stderr; optionally echo every line to stdout."""
    log_f.write(line)
    log_f.flush()
    stripped = line.rstrip("\n\r")
    if COOT_STDOUT_ALIGNMENTS_DONE_MARKER in line:
        print(stripped, file=sys.stderr, flush=True)
    if echo_all_stdout:
        print(stripped)


# Function to create Coot script for SSM superposition
def create_coot_script(
    reference_file,
    model_files,
    output_dir,
    keep_coot_open=True,
    aligned_tag="_SSMaligned2_",
):
    """One-to-many SSM. If keep_coot_open is False, skip reload-for-display and exit Coot.

    aligned_tag: substring between model basename and reference basename in output PDB names
    (default _SSMaligned2_; use _SSMaligned_ for --pattern mode, matching LSQ pattern naming).
    """
    script_content = """
import os
import sys

# Global setting for nomenclature errors
set_nomenclature_errors_on_read("ignore")

# Turn off symmetry display
set_show_symmetry_master(0)

# Load reference structure
reference_mol = read_pdb("{0}")
reference_chain = chain_ids(reference_mol)[0]
ref_name = os.path.splitext(os.path.basename("{0}"))[0]
graphics_to_ca_representation(reference_mol)

# Create output directory if it doesn't exist
if not os.path.exists("{1}"):
    os.makedirs("{1}")

# Process each model
for model_path in {2}:
    # Load model
    model_mol = read_pdb(model_path)
    model_chain = chain_ids(model_mol)[0]
    model_name = os.path.splitext(os.path.basename(model_path))[0]
    print("Superposing " + model_name + " onto " + ref_name)
    graphics_to_ca_representation(model_mol)

    try:
        # Coot: superpose_with_atom_selection is SSM (see manual). Last arg is move_copy_flag:
        # 0 = transform moving mol in place (required before write_pdb_file on this handle).
        superpose_with_atom_selection(reference_mol, model_mol, 
                                    "//" + reference_chain + "//", 
                                    "//" + model_chain + "//",
                                    0)
            
    except Exception as e:
        print("Error during superposition:", e)
        continue
    
    # Create output name with new format
    output_name = os.path.join("{1}", model_name + "__ALIGNED_TAG__" + ref_name + ".pdb")
    
    # Save aligned structure
    write_pdb_file(model_mol, output_name)

# Close all existing molecules
for i in range(graphics_n_molecules()):
    close_molecule(i)

__POST_CLOSE__
"""
    post_open = """
print("FOLDKIT_ALIGNMENTS_DONE: All superpositions written to disk. Reloading structures in Coot for inspection...")
# Start new Coot window with reference structure
handle_read_draw_molecule_with_recentre("{0}", 0)
ref_mol = graphics_n_molecules() - 1
set_molecule_bonds_colour_map_rotation(ref_mol, 0)
graphics_to_ca_representation(int(ref_mol))

# Load and display all aligned structures
for model_path in {2}:
    model_name = os.path.splitext(os.path.basename(model_path))[0]
    aligned_path = os.path.join("{1}", model_name + "__ALIGNED_TAG__" + ref_name + ".pdb")
    handle_read_draw_molecule_with_recentre(aligned_path, 0)
    mol = graphics_n_molecules() - 1
    set_molecule_bonds_colour_map_rotation(mol, 21 * (mol - ref_mol))
    graphics_to_ca_representation(int(mol))

# Get centre coordinates of reference molecule
x, y, z = molecule_centre(ref_mol)
set_rotation_centre(x, y, z)

# coot_real_exit(0)  # Comment out to keep Coot window open
"""
    post_batch = """
print("FOLDKIT_ALIGNMENTS_DONE: All superpositions written to disk. Exiting Coot.")
print("Aligned structures saved to: {1}")
coot_real_exit(0)
"""
    post_close = post_open if keep_coot_open else post_batch
    tpl = script_content.replace("__POST_CLOSE__", post_close)
    tpl = tpl.replace("__ALIGNED_TAG__", aligned_tag)
    return tpl.format(
        reference_file,  # {0}
        output_dir,     # {1}
        model_files,    # {2}
    )

def create_all_vs_all_ssm_script(model_files, output_dir, keep_coot_open=True):
    """Create Coot script for all-vs-all SSM superposition."""
    script_content = """
import os
import sys

# Global setting for nomenclature errors
set_nomenclature_errors_on_read("ignore")

# Turn off symmetry display
set_show_symmetry_master(0)

# Create output directory if it doesn't exist
if not os.path.exists("{0}"):
    os.makedirs("{0}")

# Get all model files
model_files = {1}

# Keep track of loaded molecules for final display
final_molecules = []

# Perform all-vs-all superposition
for ref_path in model_files:
    # Load reference structure
    reference_mol = read_pdb(ref_path)
    reference_chain = chain_ids(reference_mol)[0]
    ref_name = os.path.splitext(os.path.basename(ref_path))[0]

    graphics_to_ca_representation(reference_mol)

    print("Superposing onto reference: " + ref_name)

    # Process each model against this reference
    for model_path in model_files:
        # Skip if model is the same as reference
        if model_path == ref_path:
            continue

        model_name = os.path.splitext(os.path.basename(model_path))[0]
        print("  Superposing " + model_name + " onto " + ref_name)

        # Load model
        model_mol = read_pdb(model_path)
        model_chain = chain_ids(model_mol)[0]

        graphics_to_ca_representation(model_mol)

        try:
            superpose_with_atom_selection(reference_mol, model_mol,
                                        "//" + reference_chain + "//",
                                        "//" + model_chain + "//",
                                        0)

        except Exception as e:
            print("Error during SSM superposition of " + model_name + " onto " + ref_name + ": " + str(e))
            close_molecule(model_mol)
            continue

        # Create output name with new format
        output_name = os.path.join("{0}", model_name + "_SSMaligned2_" + ref_name + ".pdb")

        # Save aligned structure
        write_pdb_file(model_mol, output_name)

        if model_path not in [mol[1] for mol in final_molecules]:
            final_molecules.append((model_mol, model_path))
        else:
            close_molecule(model_mol)

    if ref_path not in [mol[1] for mol in final_molecules]:
        final_molecules.append((reference_mol, ref_path))
    else:
        close_molecule(reference_mol)

__ALL_VS_ALL_TAIL__
"""
    tail_open = """
print("FOLDKIT_ALIGNMENTS_DONE: All pairwise superpositions written to disk. Reloading structures in Coot for inspection...")
# Reload final set for visual inspection
print("\\nReloading structures for visual inspection...")
unique_files = list(set(model_files))

for i, file_path in enumerate(unique_files):
    mol = read_pdb(file_path)
    graphics_to_ca_representation(mol)
    set_molecule_bonds_colour_map_rotation(mol, 20 * i)
    print("Loaded for display: " + os.path.basename(file_path))

print("\\nAll-vs-all SSM superposition complete. All structures loaded in Coot for visual inspection.")
print("Aligned structures saved to: {0}")

# coot_real_exit(0)  # Comment out to keep Coot window open
"""
    tail_batch = """
print("FOLDKIT_ALIGNMENTS_DONE: All pairwise superpositions written to disk. Exiting Coot.")
for i in range(graphics_n_molecules()):
    close_molecule(i)
print("\\nAll-vs-all SSM superposition complete (batch). Aligned structures saved to: {0}")
coot_real_exit(0)
"""
    tail = tail_open if keep_coot_open else tail_batch
    script_content = script_content.replace("__ALL_VS_ALL_TAIL__", tail)
    return script_content.format(output_dir, model_files)


def create_coot_script_explicit_chains(
    reference_file, model_files, output_dir, ref_chain, model_chain, keep_coot_open=True
):
    """One-to-many SSM with explicit reference/model chains (mmdb selections; in-place superposition)."""
    script_content = """
import os
import sys

set_nomenclature_errors_on_read("ignore")
set_show_symmetry_master(0)

def get_chain_id(file_path, specified_chain):
    \"\"\"
    Get appropriate chain ID based on file format and specified chain.
    For PDB files: use specified chain as-is (typically single character)
    For CIF files: may need to adapt chain naming (can be longer strings)
    \"\"\"
    file_ext = os.path.splitext(file_path.lower())[1]

    if file_ext == '.pdb':
        return specified_chain
    elif file_ext == '.cif':
        return specified_chain
    else:
        return specified_chain

def find_available_chains(mol_id):
    \"\"\"Find available chains in a molecule using Coot functions\"\"\"
    chains = []
    try:
        chain_list = chain_ids(mol_id)
        if chain_list:
            chains = chain_list
        else:
            common_chains = ["A", "B", "C", "D", "E", "F", "G", "H", "I", "J"]
            for chain_id in common_chains:
                try:
                    n_atoms = n_atoms_in_chain(mol_id, chain_id)
                    if n_atoms > 0:
                        chains.append(chain_id)
                except:
                    continue
    except Exception as e:
        print("Error getting chains for molecule " + str(mol_id) + ": " + str(e))
        chains = ["A"]

    return chains

def validate_chain_exists(mol_id, chain_id):
    \"\"\"Validate that a specific chain exists in the molecule\"\"\"
    try:
        available_chains = chain_ids(mol_id)
        if available_chains and chain_id in available_chains:
            return True
        try:
            length = chain_length(mol_id, chain_id)
            return length > 0
        except:
            pass
        return False
    except:
        return False

reference_mol = read_pdb(r"{reference_file}")
ref_name = os.path.splitext(os.path.basename(r"{reference_file}"))[0]

if not os.path.exists(r"{output_dir}"):
    os.makedirs(r"{output_dir}")

ref_chain_id = get_chain_id(r"{reference_file}", "{ref_chain}")
model_chain_id_template = "{model_chain}"

print("Reference chain: " + ref_chain_id)
print("Model chain template: " + model_chain_id_template)

ref_chains = find_available_chains(reference_mol)
print("Available chains in reference: " + str(ref_chains))

if not validate_chain_exists(reference_mol, ref_chain_id):
    print("Warning: Reference chain '" + ref_chain_id + "' not found in reference structure")
    if ref_chains:
        print("Using first available chain '" + ref_chains[0] + "' instead")
        ref_chain_id = ref_chains[0]
    else:
        print("No chains found in reference structure - using 'A' as fallback")
        ref_chain_id = "A"

for model_path in {model_files}:
    model_mol = read_pdb(model_path)
    model_name = os.path.splitext(os.path.basename(model_path))[0]

    model_chain_id = get_chain_id(model_path, model_chain_id_template)

    print("Superposing " + model_name + " (chain " + model_chain_id + ") onto reference (chain " + ref_chain_id + ")")

    model_chains = find_available_chains(model_mol)
    print("Available chains in " + model_name + ": " + str(model_chains))

    if not validate_chain_exists(model_mol, model_chain_id):
        print("Warning: Model chain '" + model_chain_id + "' not found in " + model_name)
        if model_chains:
            print("Using first available chain '" + model_chains[0] + "' instead")
            model_chain_id = model_chains[0]
        else:
            print("No chains found in " + model_name + " - using 'A' as fallback")
            model_chain_id = "A"

    try:
        ref_selection = "/1/" + ref_chain_id + "/*"
        model_selection = "/1/" + model_chain_id + "/*"

        print("Attempting SSM superposition with chain " + model_chain_id + " -> chain " + ref_chain_id)
        print("Reference selection: " + ref_selection)
        print("Model selection: " + model_selection)

        superpose_with_atom_selection(reference_mol, model_mol, ref_selection, model_selection, 0)
        print("SSM superposition successful for " + model_name)
    except Exception as e:
        print("Error during SSM superposition of " + model_name + ": " + str(e))
        try:
            print("Trying SSM with CA atoms only...")
            ref_ca_selection = "/1/" + ref_chain_id + "/*/CA"
            model_ca_selection = "/1/" + model_chain_id + "/*/CA"
            superpose_with_atom_selection(reference_mol, model_mol, ref_ca_selection, model_ca_selection, 0)
            print("CA-only SSM superposition successful for " + model_name)
        except Exception as e2:
            print("SSM superposition failed for " + model_name + ": " + str(e2))
            continue

    output_name = os.path.join(r"{output_dir}", model_name + "_SSMaligned2_" + ref_name + ".pdb")
    write_pdb_file(model_mol, output_name)
    print("Saved aligned model: " + output_name)

for i in range(graphics_n_molecules()):
    close_molecule(i)

__POST_CLOSE__
"""
    post_open = """
print("FOLDKIT_ALIGNMENTS_DONE: All superpositions written to disk. Reloading structures in Coot for inspection...")
handle_read_draw_molecule_with_recentre(r"{reference_file}", 0)
ref_mol = graphics_n_molecules() - 1
set_molecule_bonds_colour_map_rotation(ref_mol, 0)
graphics_to_ca_representation(int(ref_mol))

for model_path in {model_files}:
    model_name = os.path.splitext(os.path.basename(model_path))[0]
    aligned_path = os.path.join(r"{output_dir}", model_name + "_SSMaligned2_" + ref_name + ".pdb")
    if os.path.exists(aligned_path):
        handle_read_draw_molecule_with_recentre(aligned_path, 0)
        mol = graphics_n_molecules() - 1
        set_molecule_bonds_colour_map_rotation(mol, 21 * (mol - ref_mol))
        graphics_to_ca_representation(int(mol))

x, y, z = molecule_centre(ref_mol)
set_rotation_centre(x, y, z)
# coot_real_exit(0)
"""
    post_batch = """
print("FOLDKIT_ALIGNMENTS_DONE: All superpositions written to disk. Exiting Coot.")
print("Aligned structures saved to: " + r"{output_dir}")
coot_real_exit(0)
"""
    post_close = post_open if keep_coot_open else post_batch
    script_content = script_content.replace("__POST_CLOSE__", post_close)
    return script_content.format(
        reference_file=reference_file,
        model_files=model_files,
        output_dir=output_dir,
        ref_chain=ref_chain,
        model_chain=model_chain,
    )


def create_axb_ssm_script(
    ref_files,
    model_files_B,
    output_dirs_per_ref,
    keep_coot_open=False,
):
    """
    AxB SSM: for each reference in ref_files, superpose each model in model_files_B
    (skipping same path). Uses first chain per structure; Coot SSM via superpose_with_atom_selection.
    output_dirs_per_ref[i] is the output directory for ref_files[i].
    keep_coot_open False (default): coot_real_exit(0) after all alignments (batch).
    keep_coot_open True (--interactive): reload inputs and leave Coot open.
    """
    exit_line = "# coot_real_exit(0)  # AxB --interactive: keep Coot open after reload"
    if not keep_coot_open:
        exit_line = "coot_real_exit(0)  # AxB default: non-interactive batch (exit when done)"

    script_content = """
import os
import sys

set_nomenclature_errors_on_read("ignore")
set_show_symmetry_master(0)

ref_configs = {ref_configs}
model_paths = {model_paths}
keep_coot_open = {keep_coot_open_py}

for ref_path, out_dir in ref_configs:
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    reference_mol = read_pdb(ref_path)
    reference_chain = chain_ids(reference_mol)[0]
    ref_name = os.path.splitext(os.path.basename(ref_path))[0]
    print("SSM AxB: reference " + ref_name)
    for model_path in model_paths:
        if model_path == ref_path:
            continue
        model_mol = read_pdb(model_path)
        model_chain = chain_ids(model_mol)[0]
        model_name = os.path.splitext(os.path.basename(model_path))[0]
        print("  Superposing " + model_name + " onto " + ref_name)
        try:
            superpose_with_atom_selection(reference_mol, model_mol,
                "//" + reference_chain + "//",
                "//" + model_chain + "//",
                0)
        except Exception as e:
            print("Error SSM " + model_name + " onto " + ref_name + ": " + str(e))
            close_molecule(model_mol)
            continue
        output_name = os.path.join(out_dir, model_name + "_SSMaligned2_" + ref_name + ".pdb")
        write_pdb_file(model_mol, output_name)
        close_molecule(model_mol)
    close_molecule(reference_mol)

for i in range(graphics_n_molecules()):
    close_molecule(i)

if keep_coot_open:
    print("FOLDKIT_ALIGNMENTS_DONE: All AxB superpositions written to disk. Reloading structures in Coot for inspection...")
    print("\\nReloading structures for visual inspection (AxB)...")
    unique_paths = []
    seen = set()
    for p, _ in ref_configs:
        ap = os.path.abspath(p)
        if ap not in seen:
            seen.add(ap)
            unique_paths.append(p)
    for p in model_paths:
        ap = os.path.abspath(p)
        if ap not in seen:
            seen.add(ap)
            unique_paths.append(p)
    for i, file_path in enumerate(unique_paths):
        mol = read_pdb(file_path)
        graphics_to_ca_representation(mol)
        set_molecule_bonds_colour_map_rotation(mol, 20 * i)
        print("Loaded for display: " + os.path.basename(file_path))
    print("SSM AxB complete. Aligned PDBs under each reference output directory.")
else:
    print("FOLDKIT_ALIGNMENTS_DONE: All AxB superpositions written to disk. Exiting Coot.")

__EXIT_LINE__
""".replace("__EXIT_LINE__", exit_line)

    ref_configs = list(zip(ref_files, output_dirs_per_ref))
    return script_content.format(
        ref_configs=repr(ref_configs),
        model_paths=repr(model_files_B),
        keep_coot_open_py=repr(bool(keep_coot_open)),
    )


def create_axb_ssm_explicit_chains_script(
    ref_files,
    model_files_B,
    output_dirs_per_ref,
    ref_chain,
    model_chain,
    keep_coot_open=False,
):
    """AxB SSM with explicit chains (mmdb selections; superpose_with_atom_selection, move_copy_flag=0)."""
    exit_line = "# coot_real_exit(0)  # AxB --interactive: keep Coot open after reload"
    if not keep_coot_open:
        exit_line = "coot_real_exit(0)  # AxB default: non-interactive batch (exit when done)"

    script_content = """
import os
import sys

set_nomenclature_errors_on_read("ignore")
set_show_symmetry_master(0)

def get_chain_id(file_path, specified_chain):
    return specified_chain

def find_available_chains(mol_id):
    chains = []
    try:
        chain_list = chain_ids(mol_id)
        if chain_list:
            chains = chain_list
        else:
            common_chains = ["A", "B", "C", "D", "E", "F", "G", "H", "I", "J"]
            for chain_id in common_chains:
                try:
                    n_atoms = n_atoms_in_chain(mol_id, chain_id)
                    if n_atoms > 0:
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
            length = chain_length(mol_id, chain_id)
            return length > 0
        except:
            pass
        return False
    except:
        return False

ref_configs = {ref_configs}
model_paths = {model_paths}
ref_chain_spec = "{ref_chain}"
model_chain_spec = "{model_chain}"
keep_coot_open = {keep_coot_open_py}

for ref_path, out_dir in ref_configs:
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    reference_mol = read_pdb(ref_path)
    ref_name = os.path.splitext(os.path.basename(ref_path))[0]
    ref_chain_id = get_chain_id(ref_path, ref_chain_spec)
    if not validate_chain_exists(reference_mol, ref_chain_id):
        ref_chains = find_available_chains(reference_mol)
        ref_chain_id = ref_chains[0] if ref_chains else "A"
    print("SSM AxB (explicit chains): reference " + ref_name + " chain " + ref_chain_id)

    for model_path in model_paths:
        if model_path == ref_path:
            continue
        model_mol = read_pdb(model_path)
        model_name = os.path.splitext(os.path.basename(model_path))[0]
        model_chain_id = get_chain_id(model_path, model_chain_spec)
        if not validate_chain_exists(model_mol, model_chain_id):
            model_chains = find_available_chains(model_mol)
            model_chain_id = model_chains[0] if model_chains else "A"
        try:
            ref_selection = "/1/" + ref_chain_id + "/*"
            model_selection = "/1/" + model_chain_id + "/*"
            superpose_with_atom_selection(reference_mol, model_mol, ref_selection, model_selection, 0)
        except Exception as e:
            try:
                ref_ca = "/1/" + ref_chain_id + "/*/CA"
                model_ca = "/1/" + model_chain_id + "/*/CA"
                superpose_with_atom_selection(reference_mol, model_mol, ref_ca, model_ca, 0)
            except Exception as e2:
                print("Error SSM " + model_name + " onto " + ref_name + ": " + str(e2))
                close_molecule(model_mol)
                continue
        output_name = os.path.join(out_dir, model_name + "_SSMaligned2_" + ref_name + ".pdb")
        write_pdb_file(model_mol, output_name)
        close_molecule(model_mol)
    close_molecule(reference_mol)

for i in range(graphics_n_molecules()):
    close_molecule(i)

if keep_coot_open:
    print("FOLDKIT_ALIGNMENTS_DONE: All AxB superpositions written to disk. Reloading structures in Coot for inspection...")
    print("\\nReloading structures for visual inspection (AxB explicit)...")
    unique_paths = []
    seen = set()
    for p, _ in ref_configs:
        ap = os.path.abspath(p)
        if ap not in seen:
            seen.add(ap)
            unique_paths.append(p)
    for p in model_paths:
        ap = os.path.abspath(p)
        if ap not in seen:
            seen.add(ap)
            unique_paths.append(p)
    for i, file_path in enumerate(unique_paths):
        mol = read_pdb(file_path)
        graphics_to_ca_representation(mol)
        set_molecule_bonds_colour_map_rotation(mol, 20 * i)
    print("SSM AxB (explicit chains) complete.")
else:
    print("FOLDKIT_ALIGNMENTS_DONE: All AxB superpositions written to disk. Exiting Coot.")

__EXIT_LINE__
""".replace("__EXIT_LINE__", exit_line)

    ref_configs = list(zip(ref_files, output_dirs_per_ref))
    return script_content.format(
        ref_configs=repr(ref_configs),
        model_paths=repr(model_files_B),
        ref_chain=ref_chain,
        model_chain=model_chain,
        keep_coot_open_py=repr(bool(keep_coot_open)),
    )


def create_all_vs_all_ssm_explicit_chains(
    model_files, output_dir, ref_chain, model_chain, keep_coot_open=True
):
    """All-vs-all SSM with explicit chains per structure."""
    script_content = """
import os
import sys

set_nomenclature_errors_on_read("ignore")
set_show_symmetry_master(0)

def get_chain_id(file_path, specified_chain):
    file_ext = os.path.splitext(file_path.lower())[1]
    if file_ext == '.pdb':
        return specified_chain
    elif file_ext == '.cif':
        return specified_chain
    return specified_chain

def find_available_chains(mol_id):
    chains = []
    try:
        chain_list = chain_ids(mol_id)
        if chain_list:
            chains = chain_list
        else:
            common_chains = ["A", "B", "C", "D", "E", "F", "G", "H", "I", "J"]
            for chain_id in common_chains:
                try:
                    n_atoms = n_atoms_in_chain(mol_id, chain_id)
                    if n_atoms > 0:
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
            length = chain_length(mol_id, chain_id)
            return length > 0
        except:
            pass
        return False
    except:
        return False

if not os.path.exists(r"{output_dir}"):
    os.makedirs(r"{output_dir}")

model_files = {model_files}
ref_chain_spec = "{ref_chain}"
model_chain_spec = "{model_chain}"

for ref_path in model_files:
    reference_mol = read_pdb(ref_path)
    ref_name = os.path.splitext(os.path.basename(ref_path))[0]
    ref_chain_id = get_chain_id(ref_path, ref_chain_spec)
    if not validate_chain_exists(reference_mol, ref_chain_id):
        ref_chains = find_available_chains(reference_mol)
        ref_chain_id = ref_chains[0] if ref_chains else "A"
    graphics_to_ca_representation(reference_mol)
    print("Superposing onto reference: " + ref_name + " (chain " + ref_chain_id + ")")

    for model_path in model_files:
        if model_path == ref_path:
            continue
        model_mol = read_pdb(model_path)
        model_name = os.path.splitext(os.path.basename(model_path))[0]
        model_chain_id = get_chain_id(model_path, model_chain_spec)
        if not validate_chain_exists(model_mol, model_chain_id):
            model_chains = find_available_chains(model_mol)
            model_chain_id = model_chains[0] if model_chains else "A"
        graphics_to_ca_representation(model_mol)
        print("  Superposing " + model_name + " onto " + ref_name)

        try:
            ref_selection = "/1/" + ref_chain_id + "/*"
            model_selection = "/1/" + model_chain_id + "/*"
            superpose_with_atom_selection(reference_mol, model_mol, ref_selection, model_selection, 0)
        except Exception as e:
            try:
                ref_ca_selection = "/1/" + ref_chain_id + "/*/CA"
                model_ca_selection = "/1/" + model_chain_id + "/*/CA"
                superpose_with_atom_selection(reference_mol, model_mol, ref_ca_selection, model_ca_selection, 0)
            except Exception as e2:
                print("Error SSM " + model_name + " onto " + ref_name + ": " + str(e2))
                close_molecule(model_mol)
                continue

        output_name = os.path.join(r"{output_dir}", model_name + "_SSMaligned2_" + ref_name + ".pdb")
        write_pdb_file(model_mol, output_name)
        close_molecule(model_mol)

    close_molecule(reference_mol)

__ALL_VS_ALL_EXPLICIT_TAIL__
"""
    tail_open = """
print("FOLDKIT_ALIGNMENTS_DONE: All pairwise superpositions written to disk. Reloading structures in Coot for inspection...")
print("\\nReloading structures for visual inspection...")
unique_files = list(set(model_files))
for i, file_path in enumerate(unique_files):
    mol = read_pdb(file_path)
    graphics_to_ca_representation(mol)
    set_molecule_bonds_colour_map_rotation(mol, 20 * i)
print("\\nAll-vs-all SSM (explicit chains) complete. Aligned structures saved to: " + r"{output_dir}")
"""
    tail_batch = """
print("FOLDKIT_ALIGNMENTS_DONE: All pairwise superpositions written to disk. Exiting Coot.")
for i in range(graphics_n_molecules()):
    close_molecule(i)
print("\\nAll-vs-all SSM (explicit chains) complete (batch). Aligned structures saved to: " + r"{output_dir}")
coot_real_exit(0)
"""
    tail = tail_open if keep_coot_open else tail_batch
    script_content = script_content.replace("__ALL_VS_ALL_EXPLICIT_TAIL__", tail)
    return script_content.format(
        model_files=model_files,
        output_dir=output_dir,
        ref_chain=ref_chain,
        model_chain=model_chain,
    )


def _pattern_has_glob_chars(pattern: str) -> bool:
    return any(ch in pattern for ch in ["*", "?", "["])

def _matches_filter(basename: str, pattern: str) -> bool:
    """
    Match `pattern` against a file basename.
    - If pattern contains glob chars (* ? [), use shell-style wildcard matching.
    - Otherwise, preserve legacy behaviour: substring match.
    Also matches against the filename stem (without extension) for convenience.
    """
    stem, _ = os.path.splitext(basename)
    if _pattern_has_glob_chars(pattern):
        return (
            fnmatch.fnmatchcase(basename, pattern)
            or fnmatch.fnmatchcase(stem, pattern)
        )
    return (pattern in basename) or (pattern in stem)


def run_pattern_ssm_superposition(
    ref_dir,
    model_dir,
    ref_pattern,
    model_pattern,
    target_pattern=None,
    divider=None,
    strict_position=False,
    ref_file_pattern="*.pdb",
    model_file_pattern="*.cif",
    output_suffix="_SSMaligned_",
    keep_coot_open=True,
):
    """Pair reference and model files by patterns, then SSM-superpose each match in Coot (same pairing as LSQ --pattern)."""
    matches = find_ref_model_matches(
        ref_dir,
        model_dir,
        ref_pattern,
        model_pattern,
        target_pattern,
        divider,
        strict_position,
        ref_file_pattern,
        model_file_pattern,
    )
    if not matches:
        print("No matching files found to process.")
        return

    for ref_file, subdir, model_files in matches:
        subdir_name = os.path.basename(subdir)
        ref_name = os.path.splitext(os.path.basename(ref_file))[0]
        output_dir = "{}{}{}".format(subdir_name, output_suffix, ref_name)

        if os.path.exists(output_dir):
            print("Directory {} already exists. Skipping this set.".format(output_dir))
            continue

        os.makedirs(output_dir)
        print("\nProcessing {}...".format(subdir_name))
        print("Reference: {}".format(os.path.basename(ref_file)))
        print("Models to align: {}".format(len(model_files)))

        script_file = "temp_coot_script.py"
        with open(script_file, "w") as f:
            f.write(
                create_coot_script(
                    ref_file,
                    model_files,
                    output_dir,
                    keep_coot_open=keep_coot_open,
                    aligned_tag="_SSMaligned_",
                )
            )

        try:
            log_file = os.path.join(output_dir, "coot_log.txt")
            with open(log_file, "w") as log:
                log.write("# SSM alignment log\n")
                log.write(
                    "# Mode: pattern (pair reference and model files by name patterns; same pairing rules as LSQ --pattern)\n"
                )
                log.write("# Reference: {}\n".format(ref_file))
                log.write("# Models directory: {}\n".format(subdir))
                log.write("# Number of models: {}\n\n".format(len(model_files)))

            process = Popen(
                ["coot", "--script", script_file],
                stdout=PIPE,
                stderr=STDOUT,
                universal_newlines=True,
                bufsize=1,
            )
            with open(log_file, "a") as log:
                if process.stdout is not None:
                    while True:
                        output = process.stdout.readline()
                        if output == "" and process.poll() is not None:
                            break
                        if output:
                            _pipe_coot_stdout_line(output, log)
                else:
                    print("Warning: process.stdout is None, cannot capture output.")
            rc = process.wait()
            _announce_superposition_finished(
                log_file,
                rc,
                "SSM pattern mode — subdirectory: {}.".format(subdir_name),
                rmsd_format="ssm",
            )
        except Exception as e:
            print("Error during superposition: {}".format(e))
        finally:
            if os.path.exists(script_file):
                os.remove(script_file)


def main_pattern_mode():
    """CLI: --pattern [options] reference_dir model_dir ref_pattern model_pattern [target_pattern]"""
    args = sys.argv[2:]
    strict_position = False
    divider = None
    target_pattern = None
    ref_file_pattern = "*.pdb"
    model_file_pattern = "*.cif"
    output_suffix = "_SSMaligned_"
    keep_coot_open = True

    i = 0
    while i < len(args):
        if args[i] == "--strict-position":
            strict_position = True
            args.pop(i)
        elif args[i] == "--not-interactive":
            keep_coot_open = False
            args.pop(i)
        elif args[i].startswith("--divider="):
            divider = args[i].split("=", 1)[1]
            args.pop(i)
        elif args[i].startswith("--ref-file-pattern="):
            ref_file_pattern = args[i].split("=", 1)[1]
            args.pop(i)
        elif args[i].startswith("--model-file-pattern="):
            model_file_pattern = args[i].split("=", 1)[1]
            args.pop(i)
        elif args[i].startswith("--output-suffix="):
            output_suffix = args[i].split("=", 1)[1]
            args.pop(i)
        else:
            i += 1

    if len(args) < 4:
        print("Pattern mode: pair reference and model files by name patterns (superimpose_pattern_match).")
        print("This mode cannot be combined with --all-vs-all, --reference, --filter, AxB filters, etc.")
        print("")
        print(
            "Usage: python superimposition/superimpose_coot_SSM.py --pattern [options] reference_dir model_dir ref_pattern model_pattern [target_pattern]"
        )
        print("\nOptions:")
        print("  --strict-position        ref_pattern before divider and target_pattern after")
        print("  --divider=STRING")
        print("  --ref-file-pattern=GLOB  (default: \"*.pdb\")")
        print("  --model-file-pattern=GLOB (default: \"*.cif\")")
        print("  --output-suffix=STRING   output dir suffix (default: \"_SSMaligned_\")")
        print("  --not-interactive        exit Coot after each alignment set (default: keep open)")
        print("\nExamples:")
        print("  python superimposition/superimpose_coot_SSM.py --pattern /path/to/refs /path/to/models ref_id model_id")
        print("  python superimposition/superimpose_coot_SSM.py --pattern --not-interactive /path/to/refs /path/to/models ref_id model_id")
        sys.exit(1)

    ref_dir = os.path.abspath(args[0])
    model_dir = os.path.abspath(args[1])
    ref_pattern = args[2]
    model_pattern = args[3]
    if len(args) >= 5:
        target_pattern = args[4]

    if not os.path.exists(ref_dir):
        print("Error: Reference directory '{}' does not exist.".format(ref_dir))
        sys.exit(1)
    if not os.path.exists(model_dir):
        print("Error: Model directory '{}' does not exist.".format(model_dir))
        sys.exit(1)

    run_pattern_ssm_superposition(
        ref_dir,
        model_dir,
        ref_pattern,
        model_pattern,
        target_pattern,
        divider,
        strict_position,
        ref_file_pattern,
        model_file_pattern,
        output_suffix,
        keep_coot_open=keep_coot_open,
    )


def main():
    if len(sys.argv) > 1 and sys.argv[1] == "--pattern":
        main_pattern_mode()
        return

    # Parse arguments (same style as superimpose_coot_LSQ.py: options may appear anywhere)
    filter_pattern = None
    ref_filter = None
    model_filter = None
    output_dir_arg = None
    directories = []
    reference_file = None
    reference_from_flag = None
    all_vs_all = False
    ref_chain = None
    model_chain = None
    axb_keep_coot_open = False
    legacy_keep_coot_open = True
    positionals = []

    if len(sys.argv) < 2:
        print("Usage: python superimposition/superimpose_coot_SSM.py [options] reference_file dir1 [dir2 ...]")
        print("       python superimposition/superimpose_coot_SSM.py [options] --reference=REF_FILE --filter=TAG dir1 [dir2 ...]")
        print("       python superimposition/superimpose_coot_SSM.py [options] --all-vs-all dir1 [dir2 ...]")
        print(
            "       python superimposition/superimpose_coot_SSM.py --pattern [options] reference_dir model_dir ref_pattern model_pattern [target_pattern]"
        )
        print("Options:")
        print("  --filter=TAG         Substring or glob (* ? [) on basename / stem (single set)")
        print("  --ref-filter=TAG     Like --filter, but selects reference set (A) for AxB mode")
        print("  --model-filter=TAG   Like --filter, but selects model set (B) for AxB mode")
        print("  --reference=FILE     Reference structure (one-to-many); alternative to leading positional")
        print("  --ref=FILE           Same as --reference=")
        print("  --output-dir=DIR     Output dir; may contain [reference_name], [filter], *filter*")
        print("  --out-dir=DIR        Alias for --output-dir=")
        print("  --ref-chain=CHAIN    Use this chain on the reference (implies explicit-chain SSM)")
        print("  --model-chain=CHAIN  Use this chain on each model (implies explicit-chain SSM)")
        print("  --all-vs-all         Each structure as reference in turn; with --ref-filter/--model-filter: AxB mode")
        print("  --interactive        AxB only: after all alignments, reload structures and keep Coot open")
        print("  --not-interactive    One-to-many and single-set all-vs-all: exit Coot when done (default: keep open)")
        print("Without --ref-chain/--model-chain: first chain per structure (Coot SSM superposition).")
        print("Examples:")
        print("  python superimposition/superimpose_coot_SSM.py /path/to/reference.pdb dir1 dir2 dir3")
        print("  python superimposition/superimpose_coot_SSM.py --reference=/path/to/reference.pdb dir1 dir2")
        print("  python superimposition/superimpose_coot_SSM.py --ref-chain=B --model-chain=B reference.pdb models/")
        print("  python superimposition/superimpose_coot_SSM.py --all-vs-all --filter=set_a /path/to/models/")
        print("  python superimposition/superimpose_coot_SSM.py --pattern /path/to/refs /path/to/models ref_id model_id")
        print("")
        print("Pattern mode (--pattern) is a separate entry point, not combinable with the lines above:")
        print("  --pattern must be the FIRST argument after the script name.")
        print("  Do not use --pattern together with --all-vs-all, --reference/--ref, --filter,")
        print("  --ref-filter, --model-filter, or a leading reference positional in the same command.")
        sys.exit(1)

    for arg in sys.argv[1:]:
        if arg == "--all-vs-all":
            all_vs_all = True
        elif arg == "--interactive":
            axb_keep_coot_open = True
        elif arg == "--not-interactive":
            legacy_keep_coot_open = False
        elif arg.startswith("--filter="):
            filter_pattern = arg.split("=", 1)[1]
        elif arg.startswith("--ref-filter="):
            ref_filter = arg.split("=", 1)[1]
        elif arg.startswith("--model-filter="):
            model_filter = arg.split("=", 1)[1]
        elif arg.startswith("--output-dir=") or arg.startswith("--out-dir="):
            output_dir_arg = arg.split("=", 1)[1]
        elif arg.startswith("--reference="):
            reference_from_flag = os.path.abspath(arg.split("=", 1)[1])
        elif arg.startswith("--ref="):
            reference_from_flag = os.path.abspath(arg.split("=", 1)[1])
        elif arg.startswith("--ref-chain="):
            ref_chain = arg.split("=", 1)[1]
        elif arg.startswith("--model-chain="):
            model_chain = arg.split("=", 1)[1]
        elif arg.startswith("--"):
            print("Error: Unknown option '{}'.".format(arg))
            sys.exit(1)
        else:
            positionals.append(os.path.abspath(arg))

    explicit_chains = ref_chain is not None or model_chain is not None
    if explicit_chains:
        if ref_chain is None:
            ref_chain = "A"
        if model_chain is None:
            model_chain = "A"

    if all_vs_all and reference_from_flag is not None:
        print(
            "Error: Cannot use --reference/--ref together with --all-vs-all.\n"
            "  One-to-many: use --reference=FILE (or a leading positional) without --all-vs-all.\n"
            "  AxB (two sets): use --all-vs-all --ref-filter=... --model-filter=... dir_A dir_B; "
            "each structure in set A is a reference in turn (no single --reference)."
        )
        sys.exit(1)

    if all_vs_all:
        if not positionals:
            print("Error: With --all-vs-all, at least one directory is required.")
            sys.exit(1)
        directories = positionals
        reference_file = None
    else:
        if reference_from_flag is not None:
            reference_file = reference_from_flag
            if not os.path.isfile(reference_file):
                print("Error: Reference file '{}' does not exist or is not a file.".format(reference_file))
                sys.exit(1)
            if not positionals:
                print("Error: With --reference/--ref, at least one model directory is required.")
                sys.exit(1)
            directories = positionals
        else:
            if len(positionals) < 2:
                print(
                    "Error: Reference file and at least one directory required, "
                    "or use --reference=FILE with directories, "
                    "or use --all-vs-all with directory/ies."
                )
                sys.exit(1)
            reference_file = positionals[0]
            if not os.path.isfile(reference_file):
                print("Error: Reference '{}' is not a file.".format(reference_file))
                sys.exit(1)
            if explicit_chains:
                ref_ext = os.path.splitext(reference_file.lower())[1]
                if ref_ext not in (".pdb", ".cif"):
                    print(
                        f"Warning: Reference '{reference_file}' may not be .pdb or .cif; "
                        "explicit-chain mode expects PDB or mmCIF."
                    )
            directories = positionals[1:]

    # Collect model files from all specified directories
    model_files = []
    for directory in directories:
        if not os.path.exists(directory):
            print(f"Warning: Directory '{directory}' does not exist. Skipping.")
            continue
        
        # Get all PDB and CIF files from the directory and subdirectories
        dir_models = []
        # Search directly in the directory
        dir_models.extend(glob.glob(os.path.join(directory, "*.pdb")))
        dir_models.extend(glob.glob(os.path.join(directory, "*.cif")))
        
        # Search recursively in subdirectories
        dir_models.extend(glob.glob(os.path.join(directory, "**", "*.pdb"), recursive=True))
        dir_models.extend(glob.glob(os.path.join(directory, "**", "*.cif"), recursive=True))
        
        # Remove duplicates that might occur from searching both directly and recursively
        dir_models = list(set(dir_models))
        
        if not dir_models:
            print(f"Warning: No PDB or CIF files found in '{directory}' or its subdirectories. Skipping.")
            continue
        model_files.extend(dir_models)
    
    # Remove reference file from the list if it's in the input directories (one-to-many mode only)
    if reference_file:
        reference_name = os.path.basename(reference_file)
        model_files = [f for f in model_files if os.path.basename(f) != reference_name]

    # Skip prior alignment outputs when collecting inputs (avoids re-superposing outputs)
    model_files = [f for f in model_files if "SSMaligned" not in f and "LSQaligned" not in f]

    if not model_files:
        print("Error: No PDB or CIF files found in the specified directories.")
        sys.exit(1)

    # Two-set AxB mode: --ref-filter / --model-filter (requires --all-vs-all)
    use_two_sets = bool(ref_filter or model_filter)
    if use_two_sets and not all_vs_all:
        print("Error: --ref-filter/--model-filter currently require --all-vs-all (AxB mode).")
        sys.exit(1)

    if use_two_sets and filter_pattern:
        print("Warning: --filter is ignored when --ref-filter/--model-filter are used.")
        filter_pattern = None

    ref_files = None
    model_files_B = None

    if use_two_sets:
        # Build reference set A
        if ref_filter:
            ref_files = [
                f for f in model_files
                if _matches_filter(os.path.basename(f), ref_filter)
            ]
        else:
            ref_files = list(model_files)

        # Build model set B
        if model_filter:
            model_files_B = [
                f for f in model_files
                if _matches_filter(os.path.basename(f), model_filter)
            ]
        else:
            model_files_B = list(model_files)

        if not ref_files:
            if ref_filter is not None:
                print(f"Error: No structures match --ref-filter='{ref_filter}'.")
            else:
                print(
                    "Error: Reference set (A) is empty (no structure files available; "
                    "check input directories)."
                )
            sys.exit(1)
        if not model_files_B:
            if model_filter is not None:
                print(f"Error: No structures match --model-filter='{model_filter}'.")
            else:
                print(
                    "Error: Model set (B) is empty (no structure files available; "
                    "check input directories)."
                )
            sys.exit(1)

        print(f"AxB mode: {len(ref_files)} references (A), {len(model_files_B)} models (B).")

    # Single-set filter (legacy --filter)
    if not use_two_sets and filter_pattern:
        filtered_models = []
        for model in model_files:
            if _matches_filter(os.path.basename(model), filter_pattern):
                filtered_models.append(model)
        if not filtered_models:
            print(f"Error: No models match the filter pattern '{filter_pattern}'.")
            print("  (Matched against filename and filename without extension.)")
            print(f"  Sample of {min(10, len(model_files))} candidate basenames:")
            for m in model_files[:10]:
                print(f"    {os.path.basename(m)}")
            if len(model_files) > 10:
                print(f"    ... and {len(model_files) - 10} more")
            sys.exit(1)
        model_files = filtered_models
        print(f"Applied filter '{filter_pattern}', found {len(model_files)} matching models.")

    if all_vs_all and not use_two_sets and len(model_files) < 2:
        print("Error: All-vs-all requires at least 2 structure files.")
        sys.exit(1)

    if explicit_chains:
        if reference_file:
            print(f"Reference file: {reference_file}")
        print(f"Explicit chains — reference: {ref_chain}, model: {model_chain}")
        print(f"Directories: {directories}")
        if filter_pattern:
            print(f"Filter pattern: {filter_pattern}")
        if all_vs_all:
            print("Mode: all-vs-all (explicit chains)")
        print()

    # Convenience printer for structure lists (for user feedback)
    def _print_structures(label, files):
        print(f"{label} ({len(files)} structures):")
        for f in files:
            print(f"  - {os.path.basename(f)}")

    if use_two_sets:
        # AxB mode: one Coot script for all ref x model pairs (batch exits Coot unless --interactive)
        _print_structures("Reference set (A)", ref_files)
        _print_structures("Model set (B)", model_files_B)

        ref_stems = [os.path.splitext(os.path.basename(p))[0] for p in ref_files]
        stem_counts = {}
        for s in ref_stems:
            stem_counts[s] = stem_counts.get(s, 0) + 1

        output_dirs_per_ref = []
        for ref_path, ref_stem in zip(ref_files, ref_stems):
            if stem_counts[ref_stem] > 1:
                ref_name = "{}__{}".format(
                    ref_stem,
                    hashlib.sha256(os.path.abspath(ref_path).encode("utf-8")).hexdigest()[:8],
                )
            else:
                ref_name = ref_stem
            if output_dir_arg:
                od = expand_output_dir_pattern(
                    output_dir_arg, ref_name=ref_name, filter_pattern=None
                )
            else:
                od = f"SSMaligned2_{ref_name}"
            output_dirs_per_ref.append(od)

        if any(stem_counts[s] > 1 for s in stem_counts):
            print(
                "Note: Duplicate reference basenames: each output dir uses the basename plus "
                "__<8 hex chars> from the reference file path so directories stay distinct."
            )

        # One prompt if any output dir already exists
        existing = [d for d in output_dirs_per_ref if os.path.exists(d)]
        if existing:
            response = input(
                f"{len(existing)} output director(y/ies) already exist. "
                "Files may be overwritten. Continue? (y/n): "
            )
            if response.lower() != 'y':
                print("Operation cancelled.")
                sys.exit(0)

        for od in output_dirs_per_ref:
            if not os.path.exists(od):
                os.makedirs(od)

        tag_parts = ["AxB"]
        if ref_filter:
            tag_parts.append(f"ref_{sanitize_pattern_for_filename(ref_filter)}")
        if model_filter:
            tag_parts.append(f"model_{sanitize_pattern_for_filename(model_filter)}")
        log_suffix = "_".join(tag_parts)
        log_basename = f"coot_log_{log_suffix}.txt"
        log_file = log_basename
        print(
            "Mode: SSM AxB (non-interactive batch, Coot exits when done)"
            if not axb_keep_coot_open
            else "Mode: SSM AxB (interactive: Coot stays open after reload)"
        )
        print(f"Creating log file at: {log_file}")

        script_file = "temp_coot_script.py"
        if os.path.exists(script_file):
            print(f"Warning: Temporary file '{script_file}' exists and will be overwritten")

        if explicit_chains:
            script_body = create_axb_ssm_explicit_chains_script(
                ref_files,
                model_files_B,
                output_dirs_per_ref,
                ref_chain,
                model_chain,
                keep_coot_open=axb_keep_coot_open,
            )
        else:
            script_body = create_axb_ssm_script(
                ref_files,
                model_files_B,
                output_dirs_per_ref,
                keep_coot_open=axb_keep_coot_open,
            )

        with open(script_file, "w") as f:
            f.write(script_body)

        try:
            print("Starting Coot process (SSM AxB)...")
            with open(log_file, "w") as log:
                log.write("# SSM AxB alignment log\n")
                log.write("# Directories: {}\n".format(", ".join(directories)))
                if explicit_chains:
                    log.write(
                        "# Ref chain: {}  Model chain: {}\n".format(ref_chain, model_chain)
                    )
                if ref_filter:
                    log.write(f"# Ref filter: {ref_filter}\n")
                if model_filter:
                    log.write(f"# Model filter: {model_filter}\n")
                log.write(
                    f"# AxB keep_coot_open ( --interactive ): {axb_keep_coot_open}\n"
                )
                log.write(f"# References (A): {len(ref_files)}\n")
                log.write(f"# Models (B): {len(model_files_B)}\n\n")

            process = Popen(
                ["coot", "--script", script_file],
                stdout=PIPE,
                stderr=STDOUT,
                universal_newlines=True,
                bufsize=1,
            )
            with open(log_file, "a") as log:
                if process.stdout is not None:
                    while True:
                        output = process.stdout.readline()
                        if output == "" and process.poll() is not None:
                            break
                        if output:
                            _pipe_coot_stdout_line(output, log)
                else:
                    print("Warning: process.stdout is None, cannot capture output.")
            rc = process.wait()
            _announce_superposition_finished(
                log_file,
                rc,
                "SSM AxB: {} reference(s), {} model(s) in set B.".format(
                    len(ref_files), len(model_files_B)
                ),
                rmsd_format="ssm",
            )
        finally:
            if os.path.exists(script_file):
                os.remove(script_file)

        return

    # Legacy single-set modes (one-to-many and all-vs-all)
    print(f"Found {len(model_files)} structure files to align:")
    for f in model_files:
        print(f"  - {os.path.basename(f)}")

    # Output directory: user-specified or default
    if all_vs_all:
        if output_dir_arg:
            output_dir = expand_output_dir_pattern(output_dir_arg, ref_name=None, filter_pattern=filter_pattern)
        elif filter_pattern:
            output_dir = "SSMaligned_all_vs_all_{}".format(sanitize_pattern_for_filename(filter_pattern) or filter_pattern)
        else:
            output_dir = "SSMaligned_all_vs_all"
    else:
        ref_name = os.path.splitext(os.path.basename(reference_file))[0]
        if output_dir_arg:
            output_dir = expand_output_dir_pattern(output_dir_arg, ref_name=ref_name, filter_pattern=filter_pattern)
        else:
            output_dir = f"SSMaligned2_{ref_name}"

    # Log filename: include pattern when specified
    log_suffix = sanitize_pattern_for_filename(filter_pattern) if filter_pattern else ""
    log_basename = f"coot_log_{log_suffix}.txt" if log_suffix else "coot_log.txt"
    log_file = os.path.join(output_dir, log_basename)

    # Create output directory with safety check
    if os.path.exists(output_dir):
        response = input(f"Directory '{output_dir}' already exists. Files may be overwritten. Continue? (y/n): ")
        if response.lower() != 'y':
            print("Operation cancelled.")
            sys.exit(0)

    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    print(f"Creating log file at: {log_file}")

    # Create temporary script file
    script_file = "temp_coot_script.py"
    if os.path.exists(script_file):
        print(f"Warning: Temporary file '{script_file}' exists and will be overwritten")

    if all_vs_all:
        with open(script_file, "w") as f:
            if explicit_chains:
                f.write(
                    create_all_vs_all_ssm_explicit_chains(
                        model_files,
                        output_dir,
                        ref_chain,
                        model_chain,
                        keep_coot_open=legacy_keep_coot_open,
                    )
                )
            else:
                f.write(
                    create_all_vs_all_ssm_script(
                        model_files, output_dir, keep_coot_open=legacy_keep_coot_open
                    )
                )
    else:
        with open(script_file, "w") as f:
            if explicit_chains:
                f.write(
                    create_coot_script_explicit_chains(
                        reference_file,
                        model_files,
                        output_dir,
                        ref_chain,
                        model_chain,
                        keep_coot_open=legacy_keep_coot_open,
                    )
                )
            else:
                f.write(
                    create_coot_script(
                        reference_file,
                        model_files,
                        output_dir,
                        keep_coot_open=legacy_keep_coot_open,
                    )
                )

    # Run Coot with the script and capture output
    try:
        mode_note = " (all-vs-all)" if all_vs_all else ""
        if explicit_chains:
            mode_note += ", explicit chains"
        print("Starting Coot process..." + mode_note)
        with open(log_file, 'w') as log:
            if all_vs_all:
                log.write("# SSM all-vs-all alignment log\n")
                log.write("# Directories: {}\n".format(", ".join(directories)))
                if explicit_chains:
                    log.write(
                        "# Ref chain: {}  Model chain: {}\n".format(ref_chain, model_chain)
                    )
                if filter_pattern:
                    log.write("# Filter pattern: {}\n".format(filter_pattern))
                log.write(
                    "# Keep Coot open after run: {}\n".format(legacy_keep_coot_open)
                )
                log.write("# Number of models: {}\n\n".format(len(model_files)))
            else:
                log.write("# SSM alignment log\n")
                if explicit_chains:
                    log.write("# Reference: {} (chain {})\n".format(reference_file, ref_chain))
                    log.write("# Model chain: {}\n".format(model_chain))
                else:
                    log.write("# Reference: {}\n".format(reference_file))
                log.write("# Directories: {}\n".format(", ".join(directories)))
                log.write("# Number of models: {}\n".format(len(model_files)))
                if filter_pattern:
                    log.write("# Filter pattern: {}\n".format(filter_pattern))
                log.write(
                    "# Keep Coot open after run: {}\n".format(legacy_keep_coot_open)
                )
                log.write("\n")
        process = Popen(["coot", "--script", script_file], 
                       stdout=PIPE, stderr=STDOUT,
                       universal_newlines=True, bufsize=1)
        with open(log_file, 'a') as log:
            if process.stdout is not None:
                while True:
                    output = process.stdout.readline()
                    if output == '' and process.poll() is not None:
                        break
                    if output:
                        _pipe_coot_stdout_line(output, log)
            else:
                print("Warning: process.stdout is None, cannot capture output.")
        rc = process.wait()
        if all_vs_all:
            summ = "SSM all-vs-all: {} structures (pairwise alignments in log).".format(
                len(model_files)
            )
        else:
            summ = "SSM one-to-many: {} model(s) superposed to reference.".format(
                len(model_files)
            )
        _announce_superposition_finished(log_file, rc, summ, rmsd_format="ssm")
    finally:
        if os.path.exists(script_file):
            os.remove(script_file)

if __name__ == "__main__":
    main()
