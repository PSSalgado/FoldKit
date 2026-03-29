#!/usr/bin/env python3

import os
import glob
import re
import csv
import pandas as pd
from collections import defaultdict

def extract_model_name(filename):
    """Extract model name from PDB filename."""
    return os.path.splitext(filename)[0]

def process_rmsd_file(rmsd_file):
    """Process a single RMSD values file and return a dictionary of alignments."""
    alignments = defaultdict(dict)
    
    with open(rmsd_file, 'r') as f:
        current_alignment = None
        current_rmsd = None
        
        for line in f:
            line = line.strip()
            
            # Look for alignment line
            if line.startswith('Alignment:'):
                alignment_text = line.replace('Alignment:', '').strip()
                # Extract model names from "Aligning X to Y"
                match = re.match(r'Aligning (.*?) to (.*?)$', alignment_text)
                if match:
                    model1, model2 = match.groups()
                    current_alignment = (model1, model2)
            
            # Look for RMSD value
            elif 'core rmsd achieved:' in line:
                rmsd_match = re.search(r'core rmsd achieved:\s*([\d.]+)', line)
                if rmsd_match:
                    current_rmsd = float(rmsd_match.group(1))
                    if current_alignment:
                        model1, model2 = current_alignment
                        alignments[model1][model2] = current_rmsd
                        # Also store reverse alignment
                        alignments[model2][model1] = current_rmsd
    
    return alignments

def create_rmsd_table(base_dir):
    """Create CSV tables for each subdomain's RMSD values."""
    # Find all rmsd_values_*.txt files
    rmsd_files = glob.glob(os.path.join(base_dir, "**", "rmsd_values_*.txt"), recursive=True)
    
    # Create list to store dataframes for combined table
    all_dataframes = []
    
    for rmsd_file in rmsd_files:
        print(f"Processing {rmsd_file}")
        
        # Get subdomain from filename
        subdomain = os.path.basename(rmsd_file).replace('rmsd_values_', '').replace('.txt', '')
        
        # Process the RMSD file
        alignments = process_rmsd_file(rmsd_file)
        
        # Get all unique model names
        all_models = sorted(set(alignments.keys()))
        
        # Create output CSV filename
        output_file = os.path.join(os.path.dirname(rmsd_file), f"rmsd_table_{subdomain}.csv")
        
        # Create a dataframe for this subdomain
        df_data = []
        for model1 in all_models:
            row = {'Model': model1}
            for model2 in all_models:
                if model1 == model2:
                    row[model2] = '-'
                else:
                    rmsd = alignments.get(model1, {}).get(model2, '')
                    row[model2] = str(rmsd) if rmsd != '' else ''
            df_data.append(row)
            
        df = pd.DataFrame(df_data)
        
        # Write CSV file for individual subdomain
        df.to_csv(output_file, index=False)
        print(f"Created {output_file}")
        
        # Add subdomain column for combined table
        df['Subdomain'] = subdomain
        all_dataframes.append(df)
    
    # Create combined CSV if we have data
    if all_dataframes:
        # Combine all dataframes
        combined_df = pd.concat(all_dataframes, ignore_index=True)
        
        # Create output CSV filename for combined data
        combined_output_file = os.path.join(base_dir, "combined_rmsd_table.csv")
        
        # Write combined CSV file
        combined_df.to_csv(combined_output_file, index=False)
        print(f"Created combined table: {combined_output_file}")

def main():
    import argparse
    parser = argparse.ArgumentParser(description="Create CSV tables from RMSD values files")
    parser.add_argument("base_dir", help="Base directory containing RMSD value files")
    args = parser.parse_args()
    
    create_rmsd_table(args.base_dir)

if __name__ == "__main__":
    main() 