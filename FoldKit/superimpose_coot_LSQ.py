import hashlib
import os
import sys
from subprocess import Popen, PIPE, STDOUT
import glob
import itertools
import re
import fnmatch

from superimpose_pattern_match import find_ref_model_matches

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


def _announce_superposition_finished(log_file, exit_code, summary=""):
    """
    Print a clear completion line to the terminal. While Coot runs, stdout/stderr
    are only copied to the log file, so users otherwise see no progress until the end.
    """
    print(flush=True)
    print("=" * 72, flush=True)
    if exit_code == 0:
        print(
            "Finished: all superpositions for this run are complete (Coot exit code 0).",
            flush=True,
        )
    else:
        print(
            "Coot process ended (exit code {}). Check the log if outputs look wrong.".format(
                exit_code
            ),
            flush=True,
        )
    if summary:
        print(summary, flush=True)
    print("Log: {}".format(log_file), flush=True)
    print("=" * 72, flush=True)


# Function to create Coot script for superposition
def create_coot_script(reference_file, model_files, output_dir, keep_coot_open=True):
    """Legacy one-to-many template. If keep_coot_open is False, skip reload-for-display and exit Coot."""
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

__POST_CLOSE__
"""
    post_open = """
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

# coot_real_exit(0)  # Comment out to keep Coot window open
"""
    post_batch = """
print("Aligned structures saved to: {1}")
coot_real_exit(0)
"""
    post_close = post_open if keep_coot_open else post_batch
    tpl = script_content.replace("__POST_CLOSE__", post_close)
    return tpl.format(
        reference_file,  # {0}
        output_dir,     # {1}
        output_dir,     # {2} - for display purposes
        model_files,    # {3}
    )

def create_lsq_script(
    reference_file, model_files, output_dir, aligned_tag="_LSQaligned2_", keep_coot_open=True
):
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

__TAIL__
"""
    tail_open = """
# Set color rotation to make models visually distinguishable
for i, mol in enumerate(loaded_molecules):
    set_molecule_bonds_colour_map_rotation(mol, 20 * i)

print("Superposition complete. All structures loaded in Coot for visual inspection.")
print("Reference: " + ref_name)
print("Aligned structures saved to: {1}")

# Keep Coot open for manual visualization
# coot_real_exit(0)  # Commented out to keep Coot window open
"""
    tail_batch = """
for i in range(graphics_n_molecules()):
    close_molecule(i)
print("Superposition complete (batch). Aligned structures saved to: {1}")
coot_real_exit(0)
"""
    tail = tail_open if keep_coot_open else tail_batch
    script_content = script_content.replace("__TAIL__", tail)
    script_content = script_content.replace("__ALIGNED_TAG__", aligned_tag)
    return script_content.format(reference_file, output_dir, model_files)

def create_all_vs_all_lsq_script(model_files, output_dir, keep_coot_open=True):
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

__ALL_VS_ALL_TAIL__
"""
    tail_open = """
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
"""
    tail_batch = """
for i in range(graphics_n_molecules()):
    close_molecule(i)
print("\\nAll-vs-all superposition complete (batch). Aligned structures saved to: {0}")
coot_real_exit(0)
"""
    tail = tail_open if keep_coot_open else tail_batch
    script_content = script_content.replace("__ALL_VS_ALL_TAIL__", tail)
    return script_content.format(output_dir, model_files)


def create_axb_lsq_script(
    ref_files, model_files_B, output_dirs_per_ref, keep_coot_open=False
):
    """
    AxB LSQ: for each reference in ref_files, least-squares align each model in model_files_B
    (skipping same path). First chain per structure.
    keep_coot_open False (default): call coot_real_exit(0) after all alignments (batch).
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
    graphics_to_ca_representation(reference_mol)
    print("LSQ AxB: reference " + ref_name)
    for model_path in model_paths:
        if model_path == ref_path:
            continue
        model_mol = read_pdb(model_path)
        model_chain = chain_ids(model_mol)[0]
        model_name = os.path.splitext(os.path.basename(model_path))[0]
        print("  Aligning " + model_name + " to " + ref_name)
        graphics_to_ca_representation(model_mol)
        try:
            superpose_with_atom_selection(reference_mol, model_mol,
                "//" + reference_chain + "//",
                "//" + model_chain + "//",
                0)
        except Exception as e:
            print("Error LSQ " + model_name + " to " + ref_name + ": " + str(e))
            close_molecule(model_mol)
            continue
        output_name = os.path.join(out_dir, model_name + "_LSQaligned2_" + ref_name + ".pdb")
        write_pdb_file(model_mol, output_name)
        close_molecule(model_mol)
    close_molecule(reference_mol)

for i in range(graphics_n_molecules()):
    close_molecule(i)

if keep_coot_open:
    print("\\nReloading structures for visual inspection (LSQ AxB)...")
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
    print("LSQ AxB complete.")

__EXIT_LINE__
""".replace("__EXIT_LINE__", exit_line)

    ref_configs = list(zip(ref_files, output_dirs_per_ref))
    return script_content.format(
        ref_configs=repr(ref_configs),
        model_paths=repr(model_files_B),
        keep_coot_open_py=repr(bool(keep_coot_open)),
    )


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
    keep_coot_open=True,
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
                    keep_coot_open=keep_coot_open,
                )
            )

        try:
            log_file = os.path.join(output_dir, "coot_log.txt")
            with open(log_file, "w") as log:
                log.write("# LSQ alignment log\n")
                log.write(
                    "# Mode: pattern (pair reference and model files by name patterns; same pairing rules as SSM --pattern)\n"
                )
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
                if process.stdout is not None:
                    while True:
                        output = process.stdout.readline()
                        if output == "" and process.poll() is not None:
                            break
                        if output:
                            log.write(output)
                            log.flush()
                else:
                    print("Warning: process.stdout is None, cannot capture output.")
            rc = process.wait()
            _announce_superposition_finished(
                log_file,
                rc,
                "LSQ pattern mode — subdirectory: {}.".format(subdir_name),
            )
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
            "Usage: python superimpose_coot_LSQ.py --pattern [options] reference_dir model_dir ref_pattern model_pattern [target_pattern]"
        )
        print("\nOptions:")
        print("  --strict-position        ref_pattern before divider and target_pattern after")
        print("  --divider=STRING")
        print("  --ref-file-pattern=GLOB  (default: \"*.pdb\")")
        print("  --model-file-pattern=GLOB (default: \"*.cif\")")
        print("  --output-suffix=STRING   output dir suffix (default: \"_LSQaligned_\")")
        print("  --not-interactive        exit Coot after each alignment set (default: keep open)")
        print("\nExamples:")
        print("  python superimpose_coot_LSQ.py --pattern /path/to/refs /path/to/models ref_id model_id")
        print(
            "  python superimpose_coot_LSQ.py --pattern --ref-file-pattern=\"*final.pdb\" --model-file-pattern=\"*.pdb\" "
            "/path/to/refs /path/to/models ref_id model_id"
        )
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
        keep_coot_open=keep_coot_open,
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
        print("\nAdditional options:")
        print("  --filter=TAG         Substring or glob on basename / stem (single-set all-vs-all)")
        print("  --ref-filter=TAG     Like --filter, but selects reference set (A) for AxB mode")
        print("  --model-filter=TAG   Like --filter, but selects model set (B) for AxB mode")
        print("  --reference=FILE     Reference structure (one-to-many); alternative to interactive pick")
        print("  --ref=FILE           Same as --reference=")
        print("  --output-dir=DIR     Output dir; placeholders [reference_name], [filter], *filter*")
        print("  --out-dir=DIR        Alias for --output-dir=")
        print("  --ref-chain=CHAIN    Explicit reference chain (one-to-many; AxB LSQ still uses first chain)")
        print("  --model-chain=CHAIN  Explicit model chain (same)")
        print("  --interactive        AxB only: after all alignments, reload structures and keep Coot open")
        print("  --not-interactive    One-to-many and single-set all-vs-all: exit Coot when done (default: keep open)")
        print("\nExamples:")
        print("  python superimpose_coot_LSQ.py --reference=/path/to/ref.cif --filter=set_a /path/to/models/")
        print("  python superimpose_coot_LSQ.py --all-vs-all --filter=set_a /path/to/models/")
        print("  python superimpose_coot_LSQ.py --all-vs-all --ref-filter=set_a --model-filter=set_b models_dir1 models_dir2")
        print("  python superimpose_coot_LSQ.py --pattern /path/to/refs /path/to/models ref_id model_id")
        print("")
        print("Pattern mode (--pattern) is a separate entry point, not combinable with the lines above:")
        print("  --pattern must be the FIRST argument after the script name.")
        print("  Do not use --pattern together with --all-vs-all, --reference/--ref, --filter,")
        print("  --ref-filter, --model-filter, or interactive directory-only pick in the same command.")
        sys.exit(1)
    
    # Parse arguments
    all_vs_all = False
    filter_pattern = None
    ref_filter = None
    model_filter = None
    reference_file = None
    output_dir_arg = None
    directories = []
    ref_chain = None
    model_chain = None
    axb_keep_coot_open = False
    legacy_keep_coot_open = True

    # Allow --output-dir anywhere (e.g. after positionals); only this pass sets it
    for a in sys.argv[1:]:
        if a.startswith("--output-dir=") or a.startswith("--out-dir="):
            output_dir_arg = a.split("=", 1)[1]
    
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
        elif arg.startswith("--reference=") or arg.startswith("--ref="):
            reference_file = arg.split("=", 1)[1]
            reference_file = os.path.abspath(reference_file)
            if not os.path.exists(reference_file):
                print("Error: Reference file '{}' does not exist.".format(reference_file))
                sys.exit(1)
        elif arg.startswith("--ref-chain="):
            ref_chain = arg.split("=", 1)[1]
        elif arg.startswith("--model-chain="):
            model_chain = arg.split("=", 1)[1]
        elif arg.startswith("--output-dir=") or arg.startswith("--out-dir="):
            pass  # already collected above
        elif not arg.startswith("--"):
            directories.append(os.path.abspath(arg))
    
    if not directories:
        print("Error: No directories specified.")
        sys.exit(1)
    
    # Check for conflicting options
    if all_vs_all and reference_file:
        print(
            "Error: Cannot use both --all-vs-all and --reference/--ref together.\n"
            "  One-to-many: use --reference=FILE without --all-vs-all.\n"
            "  AxB (two sets): use --all-vs-all --ref-filter=... --model-filter=... dir_A dir_B; "
            "each structure in set A is a reference in turn (no single --reference)."
        )
        sys.exit(1)

    explicit_chains = ref_chain is not None or model_chain is not None
    if explicit_chains:
        if ref_chain is None:
            ref_chain = "A"
        if model_chain is None:
            model_chain = "A"
    
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
            print("Error: No models match the filter pattern '{}'.".format(filter_pattern))
            print("  (Matched against filename and filename without extension.)")
            print("  Sample of {} candidate basenames:".format(min(10, len(model_files))))
            for m in model_files[:10]:
                print("    {}".format(os.path.basename(m)))
            if len(model_files) > 10:
                print("    ... and {} more".format(len(model_files) - 10))
            sys.exit(1)
        
        model_files = filtered_models

    # Convenience printer for structure lists (for user feedback)
    def _print_structures(label, files):
        print(f"{label} ({len(files)} structures):")
        for f in files:
            print(f"  - {os.path.basename(f)}")

    if use_two_sets:
        # AxB mode: one Coot script for all ref x model pairs (batch exits Coot unless --interactive)
        _print_structures("Reference set (A)", ref_files)
        _print_structures("Model set (B)", model_files_B)

        if explicit_chains:
            print(
                "Note: LSQ AxB uses the first chain in each structure; "
                "--ref-chain/--model-chain are not applied in this mode yet."
            )

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
                od = f"LSQaligned2_{ref_name}"
            output_dirs_per_ref.append(od)

        if any(stem_counts[s] > 1 for s in stem_counts):
            print(
                "Note: Duplicate reference basenames: each output dir uses the basename plus "
                "__<8 hex chars> from the reference file path so directories stay distinct."
            )

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
            "Mode: LSQ AxB (non-interactive batch, Coot exits when done)"
            if not axb_keep_coot_open
            else "Mode: LSQ AxB (interactive: Coot stays open after reload)"
        )
        print(f"Creating log file at: {log_file}")

        script_file = "temp_coot_script.py"
        if os.path.exists(script_file):
            print("Warning: Temporary file '{}' exists and will be overwritten".format(script_file))

        script_body = create_axb_lsq_script(
            ref_files,
            model_files_B,
            output_dirs_per_ref,
            keep_coot_open=axb_keep_coot_open,
        )

        with open(script_file, "w") as f:
            f.write(script_body)

        try:
            print("Starting Coot process (LSQ AxB)...")
            with open(log_file, "w") as log:
                log.write("# LSQ AxB alignment log\n")
                log.write("# Directories: {}\n".format(", ".join(directories)))
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
                            log.write(output)
                            log.flush()
                else:
                    print("Warning: process.stdout is None, cannot capture output.")
            rc = process.wait()
            _announce_superposition_finished(
                log_file,
                rc,
                "LSQ AxB: {} reference(s), {} model(s) in set B.".format(
                    len(ref_files), len(model_files_B)
                ),
            )
        finally:
            if os.path.exists(script_file):
                os.remove(script_file)

        return
    
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
            f.write(
                create_all_vs_all_lsq_script(
                    model_files, output_dir, keep_coot_open=legacy_keep_coot_open
                )
            )
        
        # Run Coot with the script and capture output
        try:
            print("Starting Coot process...")
            
            # Initialize main log file with header
            with open(log_file, 'w') as log:
                log.write("# LSQ all-vs-all alignment log\n")
                log.write("# Directories: {}\n".format(", ".join(directories)))
                if filter_pattern:
                    log.write("# Filter pattern: {}\n".format(filter_pattern))
                log.write(
                    "# Keep Coot open after run: {}\n".format(legacy_keep_coot_open)
                )
                log.write("# Number of models: {}\n\n".format(len(model_files)))
            
            # Run process and capture output
            process = Popen(["coot", "--script", script_file], 
                           stdout=PIPE, stderr=STDOUT,
                           universal_newlines=True, bufsize=1)
            
            # Read and write output in real-time
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
            
            rc = process.wait()
            num_alignments = len(model_files) * (len(model_files) - 1)
            _announce_superposition_finished(
                log_file,
                rc,
                "LSQ all-vs-all: {} structures, {} directed pairwise alignments. Output directory: {}.".format(
                    len(model_files),
                    num_alignments,
                    output_dir,
                ),
            )
                
        finally:
            # Clean up temporary script
            if os.path.exists(script_file):
                os.remove(script_file)
        
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
            f.write(
                create_lsq_script(
                    reference_file,
                    model_files,
                    output_dir,
                    keep_coot_open=legacy_keep_coot_open,
                )
            )
        
        # Run Coot with the script and capture output
        try:
            print("Starting Coot process...")
            
            # Initialize main log file with header
            with open(log_file, 'w') as log:
                log.write("# LSQ alignment log\n")
                log.write("# Reference: {}\n".format(reference_file))
                log.write("# Directories: {}\n".format(", ".join(directories)))
                log.write("# Number of models: {}\n".format(len(model_files)))
                if filter_pattern:
                    log.write("# Filter pattern: {}\n".format(filter_pattern))
                log.write(
                    "# Keep Coot open after run: {}\n".format(legacy_keep_coot_open)
                )

            # Run process and capture output
            process = Popen(["coot", "--script", script_file], 
                           stdout=PIPE, stderr=STDOUT,
                           universal_newlines=True, bufsize=1)
            
            # Read and write output in real-time
            with open(log_file, 'a') as log:
                if process.stdout is not None:
                    while True:
                        output = process.stdout.readline()
                        if output == '' and process.poll() is not None:
                            break
                        if output:
                            print(output.strip())
                            log.write(output)
                            log.flush()
                else:
                    print("Warning: process.stdout is None, cannot capture output.")

            rc = process.wait()
            _announce_superposition_finished(
                log_file,
                rc,
                "LSQ one-to-many: output directory {}.".format(output_dir),
            )
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