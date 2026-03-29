import os
import sys
from subprocess import Popen, PIPE, STDOUT
import glob
import itertools
import re
import fnmatch

def sanitize_pattern_for_filename(pattern):
    """Make a filter pattern safe for use in log/rmsd filenames."""
    if not pattern:
        return ""
    s = re.sub(r'[\*?\ /\\]', '_', pattern)
    return s[:64]

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

# Function to create Coot script for superposition
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
for model_path in {3}:
    # Load model
    model_mol = read_pdb(model_path)
    model_chain = chain_ids(model_mol)[0]
    
    try:
        # Use LSQ superposition (0 for LSQ mode)
        superpose_with_atom_selection(reference_mol, model_mol, 
                                    "//" + reference_chain + "//", 
                                    "//" + model_chain + "//",
                                    0)
            
    except Exception as e:
        print("Error during superposition:", e)
        continue
    
    # Create output name with new format
    model_name = os.path.splitext(os.path.basename(model_path))[0]
    output_name = os.path.join("{1}", model_name + "_LSQaligned2_" + ref_name + ".pdb")
    
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
for model_path in {3}:
    model_name = os.path.splitext(os.path.basename(model_path))[0]
    aligned_path = os.path.join("{1}", model_name + "_LSQaligned2_" + ref_name + ".pdb")
    handle_read_draw_molecule_with_recentre(aligned_path, 0)
    mol = graphics_n_molecules() - 1
    set_molecule_bonds_colour_map_rotation(mol, 21 * (mol - ref_mol))
    graphics_to_ca_representation(int(mol))

# Get center coordinates of reference molecule
x, y, z = molecule_centre(ref_mol)
set_rotation_centre(x, y, z)

# Exit Coot
# coot_real_exit(0)  # Comment out to keep Coot window open
"""
    return script_content.format(
        reference_file,  # {0}
        output_dir,     # {1}
        output_dir,     # {2} - for display purposes
        model_files     # {3}
    )

def create_lsq_script(reference_file, model_files, output_dir, aligned_tag="_LSQaligned2_"):
    """Create Coot script for LSQ superposition.

    aligned_tag: substring between model basename and reference basename in output PDB names
    (default _LSQaligned2_; use _LSQaligned_ for --pattern mode compatibility).
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

# Set reference to C-alpha representation
graphics_to_ca_representation(reference_mol)

# Create output directory if it doesn't exist
if not os.path.exists("{1}"):
    os.makedirs("{1}")

# Process each model
loaded_molecules = [reference_mol]  # Keep track of loaded molecules for display

for model_path in {2}:
    # Load model
    model_mol = read_pdb(model_path)
    model_chain = chain_ids(model_mol)[0]
    
    # Set model to C-alpha representation
    graphics_to_ca_representation(model_mol)
    
    # Add to loaded molecules list
    loaded_molecules.append(model_mol)
    
    try:
        # Use LSQ superposition (0 for LSQ mode)
        superpose_with_atom_selection(reference_mol, model_mol, 
                                    "//" + reference_chain + "//", 
                                    "//" + model_chain + "//",
                                    0)
            
    except Exception as e:
        print("Error during superposition:", e)
        continue
    
    # Create output name with new format
    model_name = os.path.splitext(os.path.basename(model_path))[0]
    output_name = os.path.join("{1}", model_name + "__ALIGNED_TAG__" + ref_name + ".pdb")
    
    # Save aligned structure
    write_pdb_file(model_mol, output_name)

# Set color rotation to make models visually distinguishable
for i, mol in enumerate(loaded_molecules):
    set_molecule_bonds_colour_map_rotation(mol, 20 * i)

print("Superposition complete. All structures loaded in Coot for visual inspection.")
print("Reference: " + ref_name)
print("Aligned structures saved to: {1}")

# Keep Coot open for manual visualization
# coot_real_exit(0)  # Commented out to keep Coot window open
""".replace("__ALIGNED_TAG__", aligned_tag).format(reference_file, output_dir, model_files)
    return script_content

def create_all_vs_all_lsq_script(model_files, output_dir):
    """Create Coot script for all-vs-all LSQ superposition."""
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
    
    # Set reference to C-alpha representation
    graphics_to_ca_representation(reference_mol)
    
    print("Using reference: " + ref_name)
    
    # Process each model against this reference
    for model_path in model_files:
        # Skip if model is the same as reference
        if model_path == ref_path:
            continue
            
        model_name = os.path.splitext(os.path.basename(model_path))[0]
        print("  Aligning " + model_name + " to " + ref_name)
        
        # Load model
        model_mol = read_pdb(model_path)
        model_chain = chain_ids(model_mol)[0]
        
        # Set model to C-alpha representation
        graphics_to_ca_representation(model_mol)
        
        try:
            # Use LSQ superposition (0 for LSQ mode)
            superpose_with_atom_selection(reference_mol, model_mol, 
                                        "//" + reference_chain + "//", 
                                        "//" + model_chain + "//",
                                        0)
                
        except Exception as e:
            print("Error during superposition of " + model_name + " to " + ref_name + ": " + str(e))
            close_molecule(model_mol)
            continue
        
        # Create output name with new format
        output_name = os.path.join("{0}", model_name + "_LSQaligned2_" + ref_name + ".pdb")
        
        # Save aligned structure
        write_pdb_file(model_mol, output_name)
        
        # Close model molecule to free memory (but keep one copy for final display)
        if model_path not in [mol[1] for mol in final_molecules]:
            final_molecules.append((model_mol, model_path))
        else:
            close_molecule(model_mol)
    
    # Close reference molecule to free memory (but keep one copy for final display)
    if ref_path not in [mol[1] for mol in final_molecules]:
        final_molecules.append((reference_mol, ref_path))
    else:
        close_molecule(reference_mol)

# Reload final set for visual inspection (load original files)
print("\\nReloading structures for visual inspection...")
loaded_for_display = []
unique_files = list(set(model_files))

for i, file_path in enumerate(unique_files):
    mol = read_pdb(file_path)
    graphics_to_ca_representation(mol)
    set_molecule_bonds_colour_map_rotation(mol, 20 * i)
    loaded_for_display.append(mol)
    print("Loaded for display: " + os.path.basename(file_path))

print("\\nAll-vs-all superposition complete. All structures loaded in Coot for visual inspection.")
print("Aligned structures saved to: {0}")

# Keep Coot open for manual visualization
# coot_real_exit(0)  # Commented out to keep Coot window open
""".format(output_dir, model_files)
    return script_content


def _normalize_ref_model_pattern(pattern):
    """Normalize a pattern by removing special characters, converting to lowercase."""
    if pattern is None:
        return None
    return re.sub(r"[_\-\s]", "", pattern.lower())


def find_ref_model_matches(
    ref_dir,
    model_dir,
    ref_pattern,
    model_pattern,
    target_pattern=None,
    divider=None,
    strict_position=False,
    ref_file_pattern="*.pdb",
    model_file_pattern="*.cif",
):
    """Find matching reference and model files based on filename / subdirectory patterns."""
    print(f"\nLooking for reference files in: {ref_dir}")
    print(f"Using reference file pattern: {ref_file_pattern}")
    print(f"Using reference search pattern: '{ref_pattern}'")
    if target_pattern:
        print(f"Using target search pattern: '{target_pattern}'")
    if divider:
        print(f"Using divider: '{divider}'")

    norm_ref_pattern = _normalize_ref_model_pattern(ref_pattern)
    norm_model_pattern = _normalize_ref_model_pattern(model_pattern)

    ref_file_glob = os.path.join(ref_dir, ref_file_pattern)
    all_ref_files = glob.glob(ref_file_glob)
    print(f"Found {len(all_ref_files)} total files matching '{ref_file_pattern}' in reference directory")

    ref_files = []
    for ref_file in all_ref_files:
        filename = os.path.basename(ref_file)
        norm_filename = _normalize_ref_model_pattern(filename)

        if divider is None:
            if norm_ref_pattern in norm_filename:
                if target_pattern is None or _normalize_ref_model_pattern(target_pattern) in norm_filename:
                    ref_files.append(ref_file)
                    print(f"Matched reference: {filename}")
            continue

        if divider in filename:
            norm_divider = _normalize_ref_model_pattern(divider)
            if strict_position and target_pattern is not None:
                parts = filename.split(divider)
                if len(parts) != 2:
                    continue
                first_part = _normalize_ref_model_pattern(parts[0])
                second_part = _normalize_ref_model_pattern(parts[1])
                if norm_ref_pattern in first_part and _normalize_ref_model_pattern(target_pattern) in second_part:
                    ref_files.append(ref_file)
                    print(f"Matched reference: {filename}")
            else:
                if norm_divider in norm_filename and norm_ref_pattern in norm_filename:
                    if target_pattern is None or _normalize_ref_model_pattern(target_pattern) in norm_filename:
                        ref_files.append(ref_file)
                        print(f"Matched reference: {filename}")

    if not ref_files:
        if divider is not None and strict_position and target_pattern is not None:
            print(
                f"No reference files found matching pattern '{ref_pattern}' before '{divider}' and '{target_pattern}' after"
            )
        elif divider is not None and target_pattern is not None:
            print(
                f"No reference files found containing both patterns '{ref_pattern}' and '{target_pattern}' with divider '{divider}'"
            )
        else:
            print(f"No reference files found containing pattern '{ref_pattern}'")
            if target_pattern is not None:
                print(f"and '{target_pattern}'")
        return []

    print(f"Found {len(ref_files)} matching reference files")

    print(f"\nLooking for model files in: {model_dir}")
    print(f"Using model file pattern: {model_file_pattern}")
    print(f"Using model search pattern: '{model_pattern}'")

    if not os.path.exists(model_dir):
        print(f"Error: Model directory '{model_dir}' does not exist")
        return []

    print("Contents of model directory:")
    try:
        for item in os.listdir(model_dir):
            item_path = os.path.join(model_dir, item)
            type_str = "DIR" if os.path.isdir(item_path) else "FILE"
            print(f"  {type_str} - {item}")
    except Exception as e:
        print(f"Error listing model directory: {e}")

    subdirs = [
        item
        for item in os.listdir(model_dir)
        if os.path.isdir(os.path.join(model_dir, item))
    ]
    model_dir_has_subdirs = len(subdirs) > 0
    print(f"Model directory has {len(subdirs)} subdirectories")

    if model_dir_has_subdirs:
        matching_subdirs = []
        for subdir in subdirs:
            subdir_path = os.path.join(model_dir, subdir)
            if norm_model_pattern in _normalize_ref_model_pattern(subdir):
                matching_subdirs.append(subdir_path)
                print(f"Matched subdirectory by name: {subdir}")
                continue
            has_matching_files = False
            for file_path in glob.glob(os.path.join(subdir_path, model_file_pattern)):
                if norm_model_pattern in _normalize_ref_model_pattern(os.path.basename(file_path)):
                    has_matching_files = True
                    break
            if has_matching_files:
                matching_subdirs.append(subdir_path)
                print(f"Matched subdirectory by file content: {subdir}")

        if not matching_subdirs:
            print(f"No subdirectories found matching pattern: {model_pattern}")
            return []

        print(f"Found {len(matching_subdirs)} matching subdirectories")
        model_files_by_subdir = {}
        for subdir in matching_subdirs:
            mf = glob.glob(os.path.join(subdir, model_file_pattern))
            if mf:
                model_files_by_subdir[subdir] = mf
                print(f"Found {len(mf)} model files in {os.path.basename(subdir)}")
        if not model_files_by_subdir:
            print("No model files found in matching subdirectories")
            return []

        matches = []
        for ref_file in ref_files:
            for subdir, model_files in model_files_by_subdir.items():
                matches.append((ref_file, subdir, model_files))
        return matches

    model_file_glob = os.path.join(model_dir, model_file_pattern)
    flat_models = glob.glob(model_file_glob)
    print(f"Found {len(flat_models)} total files matching '{model_file_pattern}' in model directory")

    matching_models = []
    for model_file in flat_models:
        filename = os.path.basename(model_file)
        if norm_model_pattern in _normalize_ref_model_pattern(filename):
            matching_models.append(model_file)
            print(f"Matched model file: {filename}")

    if not matching_models:
        print(f"No model files found matching pattern: {model_pattern}")
        return []

    print(f"Found {len(matching_models)} matching model files")
    matches = []
    for ref_file in ref_files:
        matches.append((ref_file, model_dir, matching_models))
    return matches


def run_pattern_lsq_superposition(
    ref_dir,
    model_dir,
    ref_pattern,
    model_pattern,
    target_pattern=None,
    divider=None,
    strict_position=False,
    ref_file_pattern="*.pdb",
    model_file_pattern="*.cif",
    output_suffix="_LSQaligned_",
):
    """Pair reference and model files by patterns, then LSQ-superpose each match in Coot (legacy trim_superimposeLSQ_pattern behavior)."""
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
        output_dir = f"{subdir_name}{output_suffix}{ref_name}"

        if os.path.exists(output_dir):
            print(f"Directory {output_dir} already exists. Skipping this set.")
            continue

        os.makedirs(output_dir)
        print(f"\nProcessing {subdir_name}...")
        print(f"Reference: {os.path.basename(ref_file)}")
        print(f"Models to align: {len(model_files)}")

        script_file = "temp_coot_script.py"
        with open(script_file, "w") as f:
            f.write(
                create_lsq_script(
                    ref_file,
                    model_files,
                    output_dir,
                    aligned_tag="_LSQaligned_",
                )
            )

        try:
            log_file = os.path.join(output_dir, "coot_log.txt")
            with open(log_file, "w") as log:
                log.write("# LSQ alignment log\n")
                log.write(f"# Reference: {ref_file}\n")
                log.write(f"# Models directory: {subdir}\n")
                log.write(f"# Number of models: {len(model_files)}\n\n")

            process = Popen(
                ["coot", "--script", script_file],
                stdout=PIPE,
                stderr=STDOUT,
                universal_newlines=True,
                bufsize=1,
            )
            with open(log_file, "a") as log:
                while True:
                    output = process.stdout.readline()
                    if output == "" and process.poll() is not None:
                        break
                    if output:
                        log.write(output)
                        log.flush()
            process.wait()
            print(f"Superposition completed for {subdir_name}")
        except Exception as e:
            print(f"Error during superposition: {e}")
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
    output_suffix = "_LSQaligned_"

    i = 0
    while i < len(args):
        if args[i] == "--strict-position":
            strict_position = True
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
        print(
            "Usage: python superimpose_coot_LSQ.py --pattern [options] reference_dir model_dir ref_pattern model_pattern [target_pattern]"
        )
        print("\nOptions:")
        print("  --strict-position        ref_pattern before divider and target_pattern after")
        print("  --divider=STRING")
        print("  --ref-file-pattern=GLOB  (default: \"*.pdb\")")
        print("  --model-file-pattern=GLOB (default: \"*.cif\")")
        print("  --output-suffix=STRING   output dir suffix (default: \"_LSQaligned_\")")
        print("\nExample:")
        print("  python superimpose_coot_LSQ.py --pattern /path/to/refs /path/to/models ref_tag model_tag")
        sys.exit(1)

    ref_dir = os.path.abspath(args[0])
    model_dir = os.path.abspath(args[1])
    ref_pattern = args[2]
    model_pattern = args[3]
    if len(args) >= 5:
        target_pattern = args[4]

    if not os.path.exists(ref_dir):
        print(f"Error: Reference directory '{ref_dir}' does not exist.")
        sys.exit(1)
    if not os.path.exists(model_dir):
        print(f"Error: Model directory '{model_dir}' does not exist.")
        sys.exit(1)

    run_pattern_lsq_superposition(
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
    )


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
    if len(sys.argv) > 1 and sys.argv[1] == "--pattern":
        main_pattern_mode()
        return

    # Check for minimum arguments
    if len(sys.argv) < 3:
        print("Usage: python superimpose_coot_LSQ.py [--output-dir=DIR] --all-vs-all --filter=TAG dir1 [dir2 ...]")
        print("       python superimpose_coot_LSQ.py [--output-dir=DIR] --reference=REF_FILE --filter=TAG dir1 [dir2 ...]")
        print("       python superimpose_coot_LSQ.py [--output-dir=DIR] dir1 [dir2 ...]")
        print(
            "       python superimpose_coot_LSQ.py --pattern [options] reference_dir model_dir ref_pattern model_pattern [target_pattern]"
        )
        print("Examples:")
        print("  python superimpose_coot_LSQ.py --reference=/path/to/ref.cif --filter=tag_a /path/to/models/")
        print("  python superimpose_coot_LSQ.py --all-vs-all --filter=tag_a /path/to/models/")
        print("  python superimpose_coot_LSQ.py --pattern /path/to/refs /path/to/models ref_tag model_tag")
        sys.exit(1)
    
    # Parse arguments
    all_vs_all = False
    filter_pattern = None
    reference_file = None
    output_dir_arg = None
    directories = []
    
    # Allow --output-dir anywhere (e.g. after positionals); only this pass sets it
    for a in sys.argv[1:]:
        if a.startswith("--output-dir=") or a.startswith("--out-dir="):
            output_dir_arg = a.split("=", 1)[1]
    
    for arg in sys.argv[1:]:
        if arg == "--all-vs-all":
            all_vs_all = True
        elif arg.startswith("--filter="):
            filter_pattern = arg.split("=")[1]
        elif arg.startswith("--reference=") or arg.startswith("--ref="):
            reference_file = arg.split("=")[1]
            reference_file = os.path.abspath(reference_file)
            if not os.path.exists(reference_file):
                print("Error: Reference file '{}' does not exist.".format(reference_file))
                sys.exit(1)
        elif arg.startswith("--output-dir=") or arg.startswith("--out-dir="):
            pass  # already collected above
        elif not arg.startswith("--"):
            directories.append(os.path.abspath(arg))
    
    if not directories:
        print("Error: No directories specified.")
        sys.exit(1)
    
    # Check for conflicting options
    if all_vs_all and reference_file:
        print("Error: Cannot use both --all-vs-all and --reference options together.")
        sys.exit(1)
    
    # Collect model files from all specified directories
    model_files = []
    for directory in directories:
        if not os.path.exists(directory):
            print("Warning: Directory '{}' does not exist. Skipping.".format(directory))
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
            print("Warning: No PDB or CIF files found in '{}' or its subdirectories. Skipping.".format(directory))
            continue
        
        model_files.extend(dir_models)
    
    if not model_files:
        print("Error: No PDB or CIF files found in any of the specified directories.")
        sys.exit(1)
    
    # Apply filter if specified
    if filter_pattern:
        filtered_models = []
        for model in model_files:
            if _matches_filter(os.path.basename(model), filter_pattern):
                filtered_models.append(model)
        
        if not filtered_models:
            print("Error: No models match the filter pattern '{}'.".format(filter_pattern))
            print("  (Matched against filename and filename without extension.)")
            print("  Sample of {} candidate basenames:".format(min(10, len(model_files))))
            for m in model_files[:10]:
                print("    {}".format(os.path.basename(m)))
            if len(model_files) > 10:
                print("    ... and {} more".format(len(model_files) - 10))
            sys.exit(1)
        
        model_files = filtered_models
    
    print("Found {} structure files matching criteria:".format(len(model_files)))
    for f in model_files:
        print("  - {}".format(os.path.basename(f)))
    
    if all_vs_all:
        # All-vs-all superposition
        print("\nPerforming all-vs-all superposition...")
        
        # Output directory: user-specified (patterns [reference_name], [pattern]/[filter]) or default
        if output_dir_arg:
            output_dir = expand_output_dir_pattern(output_dir_arg, ref_name=None, filter_pattern=filter_pattern)
        elif filter_pattern:
            output_dir = "LSQaligned_all_vs_all_{}".format(sanitize_pattern_for_filename(filter_pattern) or filter_pattern)
        else:
            output_dir = "LSQaligned_all_vs_all"
        
        if os.path.exists(output_dir):
            response = input("Directory '{}' already exists. Files may be overwritten. Continue? (y/n): ".format(output_dir))
            if response.lower() != 'y':
                print("Operation cancelled.")
                sys.exit(0)
        
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)
        
        log_suffix = sanitize_pattern_for_filename(filter_pattern) if filter_pattern else ""
        log_basename = "coot_log_{}.txt".format(log_suffix) if log_suffix else "coot_log.txt"
        log_file = os.path.join(output_dir, log_basename)
        print("Creating log file at: {}".format(log_file))
        
        # Create temporary script file
        script_file = "temp_coot_script.py"
        if os.path.exists(script_file):
            print("Warning: Temporary file '{}' exists and will be overwritten".format(script_file))
        
        with open(script_file, "w") as f:
            f.write(create_all_vs_all_lsq_script(model_files, output_dir))
        
        # Run Coot with the script and capture output
        try:
            print("Starting Coot process...")
            
            # Initialize main log file with header
            with open(log_file, 'w') as log:
                log.write("# LSQ all-vs-all alignment log\n")
                log.write("# Directories: {}\n".format(", ".join(directories)))
                if filter_pattern:
                    log.write("# Filter pattern: {}\n".format(filter_pattern))
                log.write("# Number of models: {}\n\n".format(len(model_files)))
            
            # Run process and capture output
            process = Popen(["coot", "--script", script_file], 
                           stdout=PIPE, stderr=STDOUT,
                           universal_newlines=True, bufsize=1)
            
            # Read and write output in real-time
            with open(log_file, 'a') as log:
                while True:
                    output = process.stdout.readline()
                    if output == '' and process.poll() is not None:
                        break
                    if output:
                        log.write(output)
                        log.flush()
            
            # Wait for process to complete
            process.wait()
            print("Coot process completed")
                
        finally:
            # Clean up temporary script
            if os.path.exists(script_file):
                os.remove(script_file)
            print("Log file written to: {}".format(log_file))
        
        # Print summary
        num_alignments = len(model_files) * (len(model_files) - 1)
        print("\nCompleted {} alignments.".format(num_alignments))
        print("Output files saved to: {}".format(output_dir))
        
    else:
        # One-to-many superposition
        if reference_file:
            # Reference file specified from command line
            print("Using reference file: {}".format(os.path.basename(reference_file)))
            
            # Remove reference file from the list of models to align if it's in there
            model_files = [f for f in model_files if os.path.abspath(f) != reference_file]
            
            if not model_files:
                print("Error: No model files remain after removing reference file.")
                sys.exit(1)
        else:
            # Original behavior - ask for reference file
            print("\nPlease select a reference file from the list above (enter the number):")
            for i, f in enumerate(model_files):
                print("  {}. {}".format(i+1, os.path.basename(f)))
            
            try:
                choice = int(input("Enter number: "))
                if choice < 1 or choice > len(model_files):
                    print("Invalid selection.")
                    sys.exit(1)
                
                reference_file = model_files[choice-1]
                print("Selected reference: {}".format(os.path.basename(reference_file)))
                
                # Remove reference file from the list of models to align
                model_files = [f for f in model_files if f != reference_file]
            except (ValueError, KeyboardInterrupt):
                print("Invalid input or operation cancelled.")
                sys.exit(1)
        
        # Output directory: user-specified (patterns [reference_name], [pattern]/[filter]) or default
        ref_name = os.path.splitext(os.path.basename(reference_file))[0]
        if output_dir_arg:
            output_dir = expand_output_dir_pattern(output_dir_arg, ref_name=ref_name, filter_pattern=filter_pattern)
        else:
            output_dir = "LSQaligned2_{}".format(ref_name)
        
        if os.path.exists(output_dir):
            response = input("Directory '{}' already exists. Files may be overwritten. Continue? (y/n): ".format(output_dir))
            if response.lower() != 'y':
                print("Operation cancelled.")
                sys.exit(0)
        
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)
        
        log_suffix = sanitize_pattern_for_filename(filter_pattern) if filter_pattern else ""
        log_basename = "coot_log_{}.txt".format(log_suffix) if log_suffix else "coot_log.txt"
        log_file = os.path.join(output_dir, log_basename)
        print("Creating log file at: {}".format(log_file))
        
        # Create temporary script file
        script_file = "temp_coot_script.py"
        if os.path.exists(script_file):
            print("Warning: Temporary file '{}' exists and will be overwritten".format(script_file))
        
        with open(script_file, "w") as f:
            f.write(create_lsq_script(reference_file, model_files, output_dir))
        
        # Run Coot with the script and capture output
        try:
            print("Starting Coot process...")
            
            # Initialize main log file with header
            with open(log_file, 'w') as log:
                log.write("# LSQ alignment log\n")
                log.write("# Reference: {}\n".format(reference_file))
                log.write("# Directories: {}\n".format(", ".join(directories)))
                if filter_pattern:
                    log.write("# Filter pattern: {}\n".format(filter_pattern))
            
            # Run process and capture output
            process = Popen(["coot", "--script", script_file], 
                           stdout=PIPE, stderr=STDOUT,
                           universal_newlines=True, bufsize=1)
            
            # Read and write output in real-time
            with open(log_file, 'a') as log:
                while True:
                    output = process.stdout.readline()
                    if output == '' and process.poll() is not None:
                        break
                    if output:
                        print(output.strip())
                        log.write(output)
                        log.flush()
                        
                # Wait for process to complete
                process.wait()
                
                print("Coot process completed.")
                print("Alignment complete. Output saved to: {}".format(output_dir))
                print("Log file available at: {}".format(log_file))
                print("\nTo extract RMSD values, run:")
                print("python extract_rmsd.py --format lsq {}".format(log_file))
        
        except Exception as e:
            print("Error running Coot: {}".format(e))
            
        except KeyboardInterrupt:
            print("\nProcess interrupted.")
            
        finally:
            # Clean up temporary script
            if os.path.exists(script_file):
                os.remove(script_file)

if __name__ == "__main__":
    main() 