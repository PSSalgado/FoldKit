#!/usr/bin/env python3
"""
Crystal Packing Analysis Pipeline
=================================

A comprehensive tool for analyzing crystal lattice packing differences
with focus on 12-molecule arrangements (4x3 grids).

Author: Crystal Analysis Pipeline
Date: 2025
"""

import glob
import os
import sys
import argparse
import numpy as np
import pandas as pd
from pathlib import Path
import warnings
warnings.filterwarnings('ignore')


def collect_structure_paths(inputs):
    """Expand inputs (files, directories, glob patterns) to a list of structure file paths."""
    paths = []
    for arg in inputs:
        arg = os.path.abspath(os.path.expanduser(arg))
        if os.path.isdir(arg):
            for ext in ('*.pdb', '*.cif', '*.ent'):
                paths.extend(glob.glob(os.path.join(arg, ext)))
        elif '*' in arg or '?' in arg:
            paths.extend(glob.glob(arg))
        elif os.path.isfile(arg):
            paths.append(arg)
    return sorted(set(paths))


def filter_paths_by_patterns(paths, patterns):
    """Return paths whose basename contains every pattern in `patterns`."""
    if not patterns:
        return list(paths)
    return [p for p in paths if all(pat in os.path.basename(p) for pat in patterns)]

# Import analysis modules - these will be created separately
try:
    from packing_metrics import PackingMetricsCalculator
except ImportError:
    print("Warning: packing_metrics module not available")
    PackingMetricsCalculator = None

try:
    from interface_analyzer import InterfaceAnalyzer
except ImportError:
    print("Warning: interface_analyzer module not available")
    InterfaceAnalyzer = None

try:
    from contact_analyzer import ContactAnalyzer
except ImportError:
    print("Warning: contact_analyzer module not available")
    ContactAnalyzer = None

class CrystalPackingAnalyzer:
    """Main class for crystal packing analysis pipeline."""
    
    def __init__(self, output_dir="crystal_analysis_output"):
        """Initialize the analyzer with output directory."""
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(exist_ok=True)
        
        # Initialize analyzers (handle cases where modules are not available)
        self.packing_calc = PackingMetricsCalculator() if PackingMetricsCalculator else None
        self.interface_analyzer = InterfaceAnalyzer() if InterfaceAnalyzer else None
        self.contact_analyzer = ContactAnalyzer() if ContactAnalyzer else None
        
        self.results = {}
        
        # Log which analyzers are available
        available = []
        if self.packing_calc: available.append("packing_metrics")
        if self.interface_analyzer: available.append("interface_analyzer")
        if self.contact_analyzer: available.append("contact_analyzer")
        
        print(f"Available analyzers: {', '.join(available) if available else 'None'}")
    
    def analyze_single_structure(self, pdb_file, structure_id=None):
        """
        Perform comprehensive analysis on a single crystal structure.
        
        Parameters:
        -----------
        pdb_file : str
            Path to PDB file
        structure_id : str, optional
            Identifier for the structure
            
        Returns:
        --------
        dict : Analysis results
        """
        if structure_id is None:
            structure_id = Path(pdb_file).stem
            
        print(f"Analyzing structure: {structure_id}")
        
        results = {
            'structure_id': structure_id,
            'pdb_file': pdb_file
        }
        
        try:
            # Phase 1: Basic packing metrics
            if self.packing_calc:
                print("  - Calculating packing metrics...")
                packing_metrics = self.packing_calc.calculate_metrics(pdb_file)
                results['packing_metrics'] = packing_metrics
            else:
                print("  - Skipping packing metrics (module not available)")
                results['packing_metrics'] = {'error': 'packing_metrics module not available'}
            
            # Phase 2: Interface analysis
            if self.interface_analyzer:
                print("  - Analyzing interfaces...")
                interface_data = self.interface_analyzer.analyze_interfaces(pdb_file)
                results['interface_analysis'] = interface_data
            else:
                print("  - Skipping interface analysis (module not available)")
                results['interface_analysis'] = {'error': 'interface_analyzer module not available'}
            
            # Phase 3: Contact analysis
            if self.contact_analyzer:
                print("  - Analyzing crystal contacts...")
                contact_data = self.contact_analyzer.analyze_contacts(pdb_file)
                results['contact_analysis'] = contact_data
            else:
                print("  - Skipping contact analysis (module not available)")
                results['contact_analysis'] = {'error': 'contact_analyzer module not available'}
            
            # Save individual results
            self._save_single_results(results, structure_id)
            
        except Exception as e:
            print(f"  ERROR: {str(e)}")
            results['error'] = str(e)
            
        return results
    
    def compare_structures(self, structure_results):
        """
        Write all per-structure results into one JSON file (batch export).

        Parameters
        ----------
        structure_results : list
            List of dicts from ``analyze_single_structure``.

        Returns
        -------
        dict
            ``{'n_structures': n, 'structures': structure_results}`` (JSON-serializable).
        """
        import json

        combined = {
            'n_structures': len(structure_results),
            'structures': structure_results,
        }
        out_path = self.output_dir / 'batch_analysis_results.json'
        with open(out_path, 'w') as f:
            json.dump(self._make_serializable(combined), f, indent=2, default=str)
        print(f"Wrote combined results for {len(structure_results)} structure(s) to {out_path}")
        return combined
    
    def _save_single_results(self, results, structure_id):
        """Save results for a single structure."""
        output_file = self.output_dir / f"{structure_id}_analysis.json"
        
        # Convert numpy arrays to lists for JSON serialization
        serializable_results = self._make_serializable(results)
        
        import json
        with open(output_file, 'w') as f:
            json.dump(serializable_results, f, indent=2)
    
    def _make_serializable(self, obj):
        """Convert numpy arrays and other non-serializable objects to JSON-safe format."""
        if isinstance(obj, dict):
            return {k: self._make_serializable(v) for k, v in obj.items()}
        elif isinstance(obj, (list, tuple)):
            return [self._make_serializable(v) for v in obj]
        elif hasattr(obj, 'tolist') and callable(getattr(obj, 'tolist')):  # numpy array
            return obj.tolist()
        elif hasattr(obj, 'item') and callable(getattr(obj, 'item')):  # numpy scalar
            return obj.item()
        else:
            return obj
    
def main():
    """Main function for command-line interface."""
    parser = argparse.ArgumentParser(
        description="Crystal Packing Analysis Pipeline",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Analyze single structure
  python crystal_packing_analyzer.py --input model_01.pdb

  # Analyze all PDBs in a directory
  python crystal_packing_analyzer.py --input models/

  # Combine per-structure JSON into one file (batch export)
  python crystal_packing_analyzer.py --input *.pdb --compare

  # Filter by set(s): one output dir per set (use {} in --output)
  python crystal_packing_analyzer.py --input *.pdb --sets set_a set_b --output "crystal_analysis_{}"
  python crystal_packing_analyzer.py --input *.pdb --sets set_a set_b --output analysis_output --dry-run
        """
    )

    parser.add_argument(
        '--input', '-i',
        nargs='+',
        required=True,
        help='Input PDB/CIF file(s), directory (all *.pdb/*.cif inside), or glob (e.g. *.pdb)'
    )

    parser.add_argument(
        '--output', '-o',
        default='crystal_analysis_output',
        help='Output directory. Use "{}" for one directory per set (e.g. crystal_analysis_{}).'
    )

    parser.add_argument(
        '--set', '-s',
        action='append',
        dest='sets',
        metavar='PATTERNS',
        help='Comma-separated patterns; file included only if basename contains ALL. Repeat for multiple sets.',
    )
    parser.add_argument(
        '--sets',
        dest='sets_multi',
        nargs='+',
        metavar='SET',
        help='Multiple set names in one go (one pattern per set). E.g. --sets set_a set_b set_c.',
    )

    parser.add_argument(
        '--compare', '-c',
        action='store_true',
        help='After processing, write batch_analysis_results.json (all structures in one file)'
    )

    parser.add_argument(
        '--verbose', '-v',
        action='store_true',
        help='Verbose output'
    )

    parser.add_argument(
        '--dry-run',
        action='store_true',
        help='Print which files would be processed per set and output dir(s), then exit.'
    )

    args = parser.parse_args()

    # Expand inputs: directories and globs to file list
    paths = collect_structure_paths(args.input)
    if not paths:
        print("No structure files found.", file=sys.stderr)
        sys.exit(1)

    set_list = []
    if getattr(args, 'sets_multi', None):
        for name in args.sets_multi:
            patterns = [name.strip()] if name.strip() else []
            if patterns:
                label = patterns[0]
                filtered = filter_paths_by_patterns(paths, patterns)
                set_list.append((label, patterns, filtered))
    if args.sets:
        for s in args.sets:
            patterns = [p.strip() for p in s.split(',') if p.strip()]
            if patterns:
                label = '_'.join(patterns)
                filtered = filter_paths_by_patterns(paths, patterns)
                set_list.append((label, patterns, filtered))
    if not set_list:
        set_list = [("all", [], list(paths))]

    if args.dry_run:
        print("Dry run: no analysis will be performed.\n", file=sys.stderr)
        for label, patterns, filtered in set_list:
            if '{}' in args.output:
                dest = args.output.replace('{}', label)
            elif len(set_list) > 1 and args.output:
                dest = args.output.rstrip(os.sep) + '_' + label
            else:
                dest = args.output
            print(f"Set {label!r} (patterns: {patterns}): {len(filtered)} file(s) -> {dest}", file=sys.stderr)
            for p in filtered:
                print(f"  {os.path.basename(p)}", file=sys.stderr)
            print(file=sys.stderr)
        return

    for label, patterns, filtered in set_list:
        if not filtered:
            print(f"No files match set {label!r} (patterns: {patterns}), skipping.", file=sys.stderr)
            continue

        if '{}' in args.output:
            out_dir = args.output.replace('{}', label)
        elif len(set_list) > 1 and args.output:
            out_dir = args.output.rstrip(os.sep) + '_' + label
        else:
            out_dir = args.output
        if args.verbose:
            print(f"Processing set {label!r}: {len(filtered)} structure(s) -> {out_dir}")
        analyzer = CrystalPackingAnalyzer(out_dir)
        all_results = []
        for pdb_file in filtered:
            if not os.path.exists(pdb_file):
                print(f"Warning: File {pdb_file} not found, skipping...", file=sys.stderr)
                continue
            result = analyzer.analyze_single_structure(pdb_file)
            all_results.append(result)

        if args.compare and all_results:
            analyzer.compare_structures(all_results)

        print(f"\nSet {label!r} complete. Results saved to: {out_dir}")
        print(f"Number of structures analyzed: {len(all_results)}")

if __name__ == "__main__":
    main() 