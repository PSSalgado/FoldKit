#!/usr/bin/env python3

import os
import sys
import fnmatch

_REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if _REPO_ROOT not in sys.path:
    sys.path.insert(0, _REPO_ROOT)

from cli_log import setup_log_from_argv

_argv_no_log, _ = setup_log_from_argv(
    script_path=__file__,
    argv=sys.argv[1:],
    inputs=[sys.argv[1]] if len(sys.argv) > 1 else [],
    pattern=None,
)
sys.argv = [sys.argv[0]] + _argv_no_log

try:
    from Bio import PDB
    from Bio.PDB.Polypeptide import protein_letters_3to1
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
    from Bio import SeqIO
    BIOPYTHON_AVAILABLE = True
except ImportError:
    print("Error: BioPython is required for this script.")
    print("Install it with: pip install biopython")
    sys.exit(1)

def extract_sequence_from_pdb(pdb_file):
    """
    Extract amino acid sequence from a PDB file.
    Returns a tuple of (sequence, protein_id) or (None, None) if extraction fails.
    """
    try:
        # Get protein ID from filename (without extension)
        protein_id = os.path.splitext(os.path.basename(pdb_file))[0]
        
        # Parse PDB file
        parser = PDB.PDBParser()
        structure = parser.get_structure(protein_id, pdb_file)
        
        # Extract sequence from first model
        model = structure[0]
        sequence = ""
        
        # Iterate through all chains
        for chain in model:
            for residue in chain:
                # Check if residue is an amino acid
                if residue.id[0] == ' ':
                    try:
                        # Convert 3-letter code to 1-letter code
                        aa = protein_letters_3to1[residue.resname]
                        sequence += aa
                    except KeyError:
                        # Skip non-standard amino acids
                        continue
        
        return sequence, protein_id
    
    except Exception as e:
        print(f"Error processing {pdb_file}: {str(e)}")
        return None, None

def find_pdb_files(base_dir, pattern):
    """
    Recursively find all PDB files matching the pattern in the given directory.
    Uses fnmatch for proper wildcard pattern matching.
    """
    pdb_files = []
    for root, dirs, files in os.walk(base_dir):
        for file in files:
            if fnmatch.fnmatch(file, pattern):
                pdb_files.append(os.path.join(root, file))
    return pdb_files

def main():
    if len(sys.argv) < 3 or len(sys.argv) > 4:
        print("Usage: python extract_sequences.py <base_directory> <output_fasta_file> [file_pattern]")
        print("Example: python extract_sequences.py /path/to/pdbs output.fasta *_pattern.pdb")
        sys.exit(1)
    
    base_dir = sys.argv[1]
    output_file = sys.argv[2]
    
    # Set default pattern or use provided pattern
    file_pattern = "*_D1.pdb"  # Default pattern
    if len(sys.argv) == 4:
        file_pattern = sys.argv[3]
    
    print(f"Searching for PDB files matching pattern: {file_pattern}")
    
    # Find all PDB files matching the pattern
    pdb_files = find_pdb_files(base_dir, file_pattern)
    
    if not pdb_files:
        print(f"No files matching '{file_pattern}' found in {base_dir}")
        sys.exit(1)
    
    print(f"Found {len(pdb_files)} PDB files")
    
    # Extract sequences and create SeqRecord objects
    records = []
    for pdb_file in pdb_files:
        sequence, protein_id = extract_sequence_from_pdb(pdb_file)
        if sequence:
            record = SeqRecord(
                Seq(sequence),
                id=protein_id,
                description=f"Extracted from {os.path.basename(pdb_file)}"
            )
            records.append(record)
            print(f"Extracted sequence from {protein_id}")
    
    # Write sequences to FASTA file
    if records:
        SeqIO.write(records, output_file, "fasta")
        print(f"\nSuccessfully wrote {len(records)} sequences to {output_file}")
    else:
        print("No sequences were successfully extracted")

if __name__ == "__main__":
    main() 