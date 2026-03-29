#!/usr/bin/env python3

import os
import sys
import argparse
import re
import glob

def extract_filename(path):
    """Extract just the filename from a full path."""
    return os.path.basename(path)

def extract_rmsd_values(log_file, output_file=None, debug=False):
    """
    Extract RMSD and related values from a Coot log file,
    extracting the exact model/reference pairs from the log content.
    
    Parameters:
    - log_file: Path to the Coot log file
    - output_file: Path to save the extracted information (optional)
    - debug: Whether to create a debug file with all decisions
    """
    if not os.path.exists(log_file):
        print(f"Error: Log file '{log_file}' not found.")
        return
    
    # Get directory name containing the log file
    directory = os.path.dirname(log_file)
    
    # Extract subdomain name from directory path
    # Look for pattern LSQ_[subdomain]_all_vs_all
    subdomain_match = re.search(r'LSQ_([^_]+)_all_vs_all', directory)
    subdomain = subdomain_match.group(1) if subdomain_match else "unknown"
    
    # Get output filename if not specified - in same directory as log file
    if not output_file:
        output_file = os.path.join(directory, f"rmsd_values_{subdomain}.txt")
    
    debug_file = os.path.join(directory, "rmsd_debug.txt") if debug else None
    
    print(f"Extracting RMSD values from {log_file}")
    print(f"Writing to {output_file}")
    
    # Read the entire log file
    with open(log_file, 'r') as log:
        log_content = log.readlines()
    
    # Output file
    with open(output_file, 'w') as out:
        reference_file = None
        current_model = None
        rmsd_block = []
        in_rmsd_block = False
        alignment_direction = None
        
        for i, line in enumerate(log_content):
            line = line.strip()
            
            # Look for reference file declaration
            if "Using reference:" in line:
                ref_name = line.split("Using reference:")[1].strip()
                reference_file = ref_name + ".pdb"
                if debug:
                    print(f"Found reference file declaration: {reference_file}")
            
            # Look for alignment direction
            elif "Aligning" in line and "to" in line:
                alignment_direction = line
                if debug:
                    print(f"Found alignment direction: {alignment_direction}")
            
            # Look for model file
            elif "Reading coordinate file:" in line and reference_file is not None:
                path = line.split("Reading coordinate file:")[1].strip()
                current_model = extract_filename(path)
                if debug:
                    print(f"Found model file: {current_model}")
            
            # Look for superposition start
            elif "superposing..." in line and current_model and reference_file:
                in_rmsd_block = True
                rmsd_block = []
                
                # Write the alignment header with direction
                out.write(f"Alignment: {alignment_direction}\n")
                if debug:
                    print(f"Writing alignment: {alignment_direction}")
            
            # Collect RMSD metrics
            elif in_rmsd_block:
                if 'core rmsd' in line or any(x in line for x in [
                    'number of residues', 'number of aligned', 'number of gaps', 
                    'number of misdirections', 'number of SSE', 'sequence identity'
                ]):
                    rmsd_block.append(line)
                    # Write the RMSD line immediately
                    out.write(f"{line}\n")
                
                # Check if we've reached the end of an RMSD block
                elif len(rmsd_block) > 0 and "core rmsd" in rmsd_block[0]:
                    # We've collected at least some RMSD info and have moved past it
                    in_rmsd_block = False
                    
                    # Add a blank line after each block for readability
                    out.write("\n")
                    
                    # Reset alignment direction for next block
                    alignment_direction = None
            
            # Check for next model (next molecule number)
            elif "Molecule " in line and " read successfully" in line:
                # This indicates a new molecule was loaded, potentially a new model
                # We don't reset reference_file because it should be constant
                in_rmsd_block = False  # End any previous RMSD block
    
    print(f"Successfully extracted RMSD values from {log_file}")
    print(f"Output written to {output_file}")
    
    return output_file

def find_and_extract_all(base_dir, output_dir=None, debug=False):
    """Find all coot_log.txt files and extract RMSD values from them."""
    print(f"Searching for log files in {base_dir}")
    
    # Find all coot_log.txt files
    log_files = glob.glob(os.path.join(base_dir, "**", "coot_log.txt"), recursive=True)
    print(f"Found {len(log_files)} log files")
    
    output_files = []
    
    for log_file in log_files:
        # If output_dir is specified, create a corresponding output file path
        # Otherwise, use the default (same directory as the log file)
        if output_dir:
            # Get relative path from base_dir to ensure we maintain directory structure
            rel_path = os.path.relpath(os.path.dirname(log_file), base_dir)
            
            # Create output directory
            out_subdir = os.path.join(output_dir, rel_path)
            os.makedirs(out_subdir, exist_ok=True)
            
            # Create output file in the new location
            out_file = os.path.join(out_subdir, "rmsd_values.txt")
        else:
            # Default: same directory as log file
            out_file = None
        
        # Extract RMSD values
        result = extract_rmsd_values(log_file, out_file, debug)
        if result:
            output_files.append(result)
    
    print(f"Successfully processed {len(output_files)} of {len(log_files)} log files")
    return output_files

def main():
    # Set up argument parser
    parser = argparse.ArgumentParser(description="Extract RMSD values from Coot log files")
    
    # Define main operation modes as subparsers
    subparsers = parser.add_subparsers(dest="command", help="Command")
    
    # Single file processing
    file_parser = subparsers.add_parser("file", help="Process a single log file")
    file_parser.add_argument("log_file", help="Path to the Coot log file")
    file_parser.add_argument("--output", "-o", help="Output file path")
    file_parser.add_argument("--debug", "-d", action="store_true", help="Create debug file")
    
    # Directory processing
    dir_parser = subparsers.add_parser("dir", help="Process all log files in a directory")
    dir_parser.add_argument("base_dir", help="Base directory to search for log files")
    dir_parser.add_argument("--output-dir", "-o", help="Output directory for RMSD values")
    dir_parser.add_argument("--debug", "-d", action="store_true", help="Create debug files")
    
    args = parser.parse_args()
    
    # Process based on command
    if args.command == "file":
        extract_rmsd_values(args.log_file, args.output, args.debug)
    elif args.command == "dir":
        find_and_extract_all(args.base_dir, args.output_dir, args.debug)
    else:
        parser.print_help()

if __name__ == "__main__":
    main()