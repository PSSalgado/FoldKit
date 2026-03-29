import os
import sys
from subprocess import Popen, PIPE, STDOUT
import glob
import re
import fnmatch

def sanitize_pattern_for_filename(pattern):
    """Make a filter pattern safe for use in log/rmsd filenames."""
    if not pattern:
        return ""
    s = re.sub(r'[\*?\ /\\]', '_', pattern)
    return s[:64]  # limit length

def expand_output_dir_pattern(template, ref_name=None, filter_pattern=None):
    """
    Expand --output-dir pattern: [reference_name] -> ref stem, [pattern] or [filter] -> sanitized filter.
    Also replaces literal *filter* with the sanitized filter (convenience placeholder).
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

# Function to create Coot script for SSM superposition
def create_coot_script(reference_file, model_files, output_dir):
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

# Create output directory if it doesn't exist
if not os.path.exists("{1}"):
    os.makedirs("{1}")

# Process each model
for model_path in {2}:
    # Load model
    model_mol = read_pdb(model_path)
    model_chain = chain_ids(model_mol)[0]
    
    try:
        # Use SSM superposition (1 for SSM mode)
        superpose_with_atom_selection(reference_mol, model_mol, 
                                    "//" + reference_chain + "//", 
                                    "//" + model_chain + "//",
                                    1)
            
    except Exception as e:
        print("Error during superposition:", e)
        continue
    
    # Create output name with new format
    model_name = os.path.splitext(os.path.basename(model_path))[0]
    output_name = os.path.join("{1}", model_name + "_SSMaligned2_" + ref_name + ".pdb")
    
    # Save aligned structure
    write_pdb_file(model_mol, output_name)

# Close all existing molecules
for i in range(graphics_n_molecules()):
    close_molecule(i)

# Start new Coot window with reference structure
handle_read_draw_molecule_with_recentre("{0}", 0)
ref_mol = graphics_n_molecules() - 1
set_molecule_bonds_colour_map_rotation(ref_mol, 0)
graphics_to_ca_representation(int(ref_mol))

# Load and display all aligned structures
for model_path in {2}:
    model_name = os.path.splitext(os.path.basename(model_path))[0]
    aligned_path = os.path.join("{1}", model_name + "_SSMaligned2_" + ref_name + ".pdb")
    handle_read_draw_molecule_with_recentre(aligned_path, 0)
    mol = graphics_n_molecules() - 1
    set_molecule_bonds_colour_map_rotation(mol, 21 * (mol - ref_mol))
    graphics_to_ca_representation(int(mol))

# Get center coordinates of reference molecule
x, y, z = molecule_centre(ref_mol)
set_rotation_centre(x, y, z)

# Exit Coot
# coot_real_exit(0)  # Comment out to keep Coot window open
""".format(
        reference_file,  # {0}
        output_dir,     # {1}
        model_files,    # {2}
    )
    
    return script_content

def create_all_vs_all_ssm_script(model_files, output_dir):
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
            # Use SSM superposition (1 for SSM mode)
            superpose_with_atom_selection(reference_mol, model_mol,
                                        "//" + reference_chain + "//",
                                        "//" + model_chain + "//",
                                        1)

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
""".format(output_dir, model_files)
    return script_content


def create_coot_script_explicit_chains(reference_file, model_files, output_dir, ref_chain, model_chain):
    """One-to-many SSM with explicit reference/model chains (mmdb selections, move_copy_flag 0)."""
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
""".format(
        reference_file=reference_file,
        model_files=model_files,
        output_dir=output_dir,
        ref_chain=ref_chain,
        model_chain=model_chain,
    )
    return script_content


def create_all_vs_all_ssm_explicit_chains(model_files, output_dir, ref_chain, model_chain):
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

print("\\nReloading structures for visual inspection...")
unique_files = list(set(model_files))
for i, file_path in enumerate(unique_files):
    mol = read_pdb(file_path)
    graphics_to_ca_representation(mol)
    set_molecule_bonds_colour_map_rotation(mol, 20 * i)
print("\\nAll-vs-all SSM (explicit chains) complete. Aligned structures saved to: " + r"{output_dir}")
""".format(
        model_files=model_files,
        output_dir=output_dir,
        ref_chain=ref_chain,
        model_chain=model_chain,
    )
    return script_content


def _pattern_has_glob_chars(pattern: str) -> bool:
    return any(ch in pattern for ch in ["*", "?", "["])

def _matches_filter(basename: str, pattern: str) -> bool:
    """
    Match `pattern` against a file basename.
    - If pattern contains glob chars (* ? [), use shell-style wildcard matching.
    - Otherwise, preserve legacy behavior: substring match.
    Also matches against the filename stem (without extension) for convenience.
    """
    stem, _ = os.path.splitext(basename)
    if _pattern_has_glob_chars(pattern):
        return (
            fnmatch.fnmatchcase(basename, pattern)
            or fnmatch.fnmatchcase(stem, pattern)
        )
    return (pattern in basename) or (pattern in stem)

def main():
    # Parse arguments
    filter_pattern = None
    output_dir_arg = None
    directories = []
    reference_file = None
    all_vs_all = False
    ref_chain = None
    model_chain = None
    args = sys.argv[1:]

    if not args:
        print("Usage: python superimpose_coot_SSM.py [options] reference_file dir1 [dir2 ...]")
        print("       python superimpose_coot_SSM.py [options] --all-vs-all dir1 [dir2 ...]")
        print("Options:")
        print("  --filter=TAG         Substring or glob (* ? [) on basename / stem")
        print("  --output-dir=DIR     Output dir; may contain [reference_name], [filter], *filter*")
        print("  --ref-chain=CHAIN    Use this chain on the reference (implies explicit-chain SSM)")
        print("  --model-chain=CHAIN  Use this chain on each model (implies explicit-chain SSM)")
        print("  --all-vs-all         Each structure as reference in turn")
        print("Without --ref-chain/--model-chain: first chain per structure, SSM mode 1 (legacy).")
        print("Examples:")
        print("  python superimpose_coot_SSM.py /path/to/reference.pdb dir1 dir2 dir3")
        print("  python superimpose_coot_SSM.py --ref-chain=B --model-chain=B ref.pdb models/")
        print("  python superimpose_coot_SSM.py --all-vs-all --filter=tag_a /path/to/models/")
        sys.exit(1)

    # Allow --output-dir anywhere
    for a in args:
        if a.startswith("--output-dir=") or a.startswith("--out-dir="):
            output_dir_arg = a.split("=", 1)[1]

    # Parse flags, then positionals
    i = 0
    while i < len(args):
        if args[i] == "--all-vs-all":
            all_vs_all = True
            i += 1
        elif args[i].startswith("--filter="):
            filter_pattern = args[i].split("=", 1)[1]
            i += 1
        elif args[i].startswith("--output-dir=") or args[i].startswith("--out-dir="):
            i += 1
        elif args[i].startswith("--ref-chain="):
            ref_chain = args[i].split("=", 1)[1]
            i += 1
        elif args[i].startswith("--model-chain="):
            model_chain = args[i].split("=", 1)[1]
            i += 1
        else:
            break

    explicit_chains = ref_chain is not None or model_chain is not None
    if explicit_chains:
        if ref_chain is None:
            ref_chain = "A"
        if model_chain is None:
            model_chain = "A"

    # Positionals: either (reference, dir1, ...) or (dir1, ...) when --all-vs-all
    positionals = [os.path.abspath(arg) for arg in args[i:] if not arg.startswith("--")]

    if all_vs_all:
        if not positionals:
            print("Error: With --all-vs-all, at least one directory is required.")
            sys.exit(1)
        directories = positionals
        reference_file = None
    else:
        if len(positionals) < 2:
            print("Error: Reference file and at least one directory required (or use --all-vs-all with directory/ies).")
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

    # Apply filter if specified
    if filter_pattern:
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

    if all_vs_all and len(model_files) < 2:
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
                        model_files, output_dir, ref_chain, model_chain
                    )
                )
            else:
                f.write(create_all_vs_all_ssm_script(model_files, output_dir))
    else:
        with open(script_file, "w") as f:
            if explicit_chains:
                f.write(
                    create_coot_script_explicit_chains(
                        reference_file, model_files, output_dir, ref_chain, model_chain
                    )
                )
            else:
                f.write(create_coot_script(reference_file, model_files, output_dir))

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
                log.write("# Number of models: {}\n\n".format(len(model_files)))
            else:
                log.write("# SSM alignment log\n")
                if explicit_chains:
                    log.write("# Reference: {} (chain {})\n".format(reference_file, ref_chain))
                    log.write("# Model chain: {}\n".format(model_chain))
                else:
                    log.write("# Reference: {}\n".format(reference_file))
                log.write("# Directories: {}\n".format(", ".join(directories)))
                if filter_pattern:
                    log.write("# Filter pattern: {}\n".format(filter_pattern))
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
                        log.write(output)
                        log.flush()
            else:
                print("Warning: process.stdout is None, cannot capture output.")
        process.wait()
        print("Coot process completed")
    finally:
        if os.path.exists(script_file):
            os.remove(script_file)
        print(f"Log file written to: {log_file}")

if __name__ == "__main__":
    main()
