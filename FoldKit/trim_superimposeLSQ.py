import os
import sys
from subprocess import Popen, PIPE, STDOUT

from trim_models import collect_models_from_directories, matches_filter_ci, run_trim_job

def create_coot_script(reference_file, model_files, output_dir):
    """Create Coot script for superposition with visualization."""
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
""".format(
        reference_file,  # {0}
        output_dir,     # {1}
        output_dir,     # {2} (unused)
        model_files,    # {3}
    )
    
    return script_content

def create_lsq_script(reference_file, model_files, output_dir):
    """Create Coot script for LSQ superposition."""
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

# Exit Coot
coot_real_exit(0)
"""
    return script_content.format(reference_file, output_dir, model_files)

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

# Perform all-vs-all superposition
for ref_path in model_files:
    # Load reference structure
    reference_mol = read_pdb(ref_path)
    reference_chain = chain_ids(reference_mol)[0]
    ref_name = os.path.splitext(os.path.basename(ref_path))[0]
    
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
        
        try:
            # Use LSQ superposition (0 for LSQ mode)
            superpose_with_atom_selection(reference_mol, model_mol, 
                                        "//" + reference_chain + "//", 
                                        "//" + model_chain + "//",
                                        0)
                
        except Exception as e:
            print("Error during superposition of " + model_name + " to " + ref_name + ": " + str(e))
            continue
        
        # Create output name with new format
        output_name = os.path.join("{0}", model_name + "_LSQaligned2_" + ref_name + ".pdb")
        
        # Save aligned structure
        write_pdb_file(model_mol, output_name)
    
    # Close reference molecule to free memory
    close_molecule(reference_mol)

# Exit Coot
coot_real_exit(0)
"""
    return script_content.format(output_dir, model_files)

def main():
    # Check for minimum arguments
    if len(sys.argv) < 2:
        print("Usage: python trim_superimposeLSQ.py [--trim] [--trim-only] [--all-vs-all] [--filter=tag_a,tag_b,...] directory1 [directory2 ...]")
        print("       python trim_models.py [--filter=tag_a,tag_b,...] directory1 [directory2 ...]   # standalone trim (no Coot)")
        print("Options:")
        print("  --trim           Trim models to match the shortest model's residue range, then superpose")
        print("  --trim-only      Only trim (no Coot); implies --trim. Same trimming logic as trim_models.py.")
        print("  --all-vs-all     Perform all-vs-all superposition")
        print("  --filter=TAG,... Comma-separated substrings; model basename must contain one of them")
        sys.exit(1)

    # Parse arguments
    trim_models = False
    trim_only = False
    all_vs_all = False
    filter_patterns = []
    directories = []
    
    for arg in sys.argv[1:]:
        if arg == "--trim":
            trim_models = True
        elif arg == "--trim-only":
            trim_only = True
            trim_models = True
        elif arg == "--all-vs-all":
            all_vs_all = True
        elif arg.startswith("--filter="):
            filter_patterns = arg.split("=")[1].split(",")
        elif not arg.startswith("--"):
            directories.append(os.path.abspath(arg))

    if trim_only and all_vs_all:
        print("Error: --trim-only cannot be used with --all-vs-all.")
        sys.exit(1)
    
    if not directories:
        print("Error: No input directories specified.")
        sys.exit(1)

    all_model_files = collect_models_from_directories(directories)
    if not all_model_files:
        print("Error: No PDB or CIF files found in any of the specified directories.")
        sys.exit(1)

    trimmed_files = []
    if trim_models:
        trimmed_files = run_trim_job(all_model_files, filter_patterns)
    elif filter_patterns:
        for pattern in filter_patterns:
            filtered_models = [
                (model, src_dir)
                for model, src_dir in all_model_files
                if matches_filter_ci(os.path.basename(model), pattern)
            ]
            if not filtered_models:
                print(f"Warning: No models match the filter pattern '{pattern}'. Skipping.")
                continue
            print(f"\nProcessing filter pattern: '{pattern}'")
            print(f"Found {len(filtered_models)} matching models.")
            trimmed_files.extend([model for model, _ in filtered_models])
    else:
        trimmed_files = [model for model, _ in all_model_files]
    
    # Now perform superposition on all trimmed files
    if not trimmed_files:
        print("Error: No models available after trimming.")
        sys.exit(1)

    if trim_only:
        print(f"\nTrim-only mode: wrote {len(trimmed_files)} structure(s); skipping Coot superposition.")
        sys.exit(0)
    
    print(f"\nFound {len(trimmed_files)} models for superposition.")
    
    # Create output directory name
    if all_vs_all:
        if filter_patterns:
            if len(filter_patterns) == 1:
                output_dir = f"LSQ_{filter_patterns[0]}_all_vs_all"
            else:
                output_dir = f"LSQ_{'_'.join(filter_patterns)}_all_vs_all"
        else:
            output_dir = "LSQ_all_vs_all"
    else:
        output_dir = "LSQaligned"
    
    # Create output directory with safety check
    if os.path.exists(output_dir):
        response = input(f"Directory '{output_dir}' already exists. Files may be overwritten. Continue? (y/n): ")
        if response.lower() != 'y':
            print("Operation cancelled.")
            sys.exit(0)
    
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    # Set up logging to file
    log_file = os.path.join(output_dir, "coot_log.txt")
    print(f"Creating log file at: {log_file}")
    
    # Create temporary script file
    script_file = "temp_coot_script.py"
    if os.path.exists(script_file):
        print(f"Warning: Temporary file '{script_file}' exists and will be overwritten")
    
    if all_vs_all:
        # All-vs-all superposition
        with open(script_file, "w") as f:
            f.write(create_all_vs_all_lsq_script(trimmed_files, output_dir))
        
        # Run Coot with the script and capture output
        try:
            print("Starting Coot process for all-vs-all superposition...")
            
            # Initialize main log file with header
            with open(log_file, 'w') as log:
                log.write("# LSQ all-vs-all alignment log\n")
                log.write(f"# Directories: {', '.join(directories)}\n")
                if filter_patterns:
                    log.write(f"# Filter patterns: {', '.join(filter_patterns)}\n")
                if trim_models:
                    log.write("# Models were trimmed to match shortest model in each pattern group\n")
                log.write(f"# Number of models: {len(trimmed_files)}\n\n")
            
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
            print(f"Log file written to: {log_file}")
        
        # Print summary
        num_alignments = len(trimmed_files) * (len(trimmed_files) - 1)
        print(f"\nCompleted {num_alignments} alignments.")
        print(f"Output files saved to: {output_dir}")
        
    else:
        # Original behavior - ask for reference file
        print("\nPlease select a reference file from the list above (enter the number):")
        for i, f in enumerate(trimmed_files):
            print(f"  {i+1}. {os.path.basename(f)}")
        
        try:
            choice = int(input("Enter number: "))
            if choice < 1 or choice > len(trimmed_files):
                print("Invalid selection.")
                sys.exit(1)
            
            reference_file = trimmed_files[choice-1]
            print(f"Selected reference: {os.path.basename(reference_file)}")
            
            # Remove reference file from the list of models to align
            model_files = [f for f in trimmed_files if f != reference_file]
            
            # Create output directory name based on reference model
            ref_name = os.path.splitext(os.path.basename(reference_file))[0]
            output_dir = f"LSQaligned2_{ref_name}"
            
            # Create output directory with safety check
            if os.path.exists(output_dir):
                response = input(f"Directory '{output_dir}' already exists. Files may be overwritten. Continue? (y/n): ")
                if response.lower() != 'y':
                    print("Operation cancelled.")
                    sys.exit(0)
            
            # Create output directory if it doesn't exist
            if not os.path.exists(output_dir):
                os.makedirs(output_dir)
            
            # Set up logging to file
            log_file = os.path.join(output_dir, "coot_log.txt")
            print(f"Creating log file at: {log_file}")
            
            # Create temporary script file
            script_file = "temp_coot_script.py"
            if os.path.exists(script_file):
                print(f"Warning: Temporary file '{script_file}' exists and will be overwritten")
            
            with open(script_file, "w") as f:
                f.write(create_lsq_script(reference_file, model_files, output_dir))
            
            # Run Coot with the script and capture output
            try:
                print("Starting Coot process...")
                
                # Initialize main log file with header
                with open(log_file, 'w') as log:
                    log.write("# LSQ alignment log\n")
                    log.write(f"# Reference: {reference_file}\n")
                    log.write(f"# Directories: {', '.join(directories)}\n")
                    if filter_patterns:
                        log.write(f"# Filter patterns: {', '.join(filter_patterns)}\n")
                    if trim_models:
                        log.write("# Models were trimmed to match shortest model in each pattern group\n")
                
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
                print(f"Log file written to: {log_file}")
            
        except ValueError:
            print("Invalid input. Please enter a number.")
            sys.exit(1)

if __name__ == "__main__":
    main() 