#!/usr/bin/env python3
"""
Crystal Packing Analysis Pipeline
=================================

A comprehensive tool for analysing crystal lattice packing differences
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
    from interface_analyser_base import InterfaceAnalyser, InterfaceAnalyserEC
except ImportError:
    print("Warning: interface_analyser_base module not available")
    InterfaceAnalyser = None
    InterfaceAnalyserEC = None

try:
    from contact_analyser import ContactAnalyser
except ImportError:
    print("Warning: contact_analyser module not available")
    ContactAnalyser = None

class CrystalPackingAnalyser:
    """Main class for crystal packing analysis pipeline."""
    
    def __init__(
        self,
        output_dir="crystal_analysis_output",
        *,
        interface_context: str = "asu",
        interface_metrics: str = "charge",
    ):
        """Initialise the analyser with output directory."""
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(exist_ok=True)
        self.interface_context = (interface_context or "asu").strip().lower()
        
        # Initialise analysers (handle cases where modules are not available)
        self.packing_calc = PackingMetricsCalculator() if PackingMetricsCalculator else None
        mode = (interface_metrics or "charge").strip().lower()
        self.interface_metrics = mode
        if mode == "ec":
            self.interface_analyser = InterfaceAnalyserEC(ec_mode="full") if InterfaceAnalyserEC else None
        else:
            self.interface_analyser = InterfaceAnalyser() if InterfaceAnalyser else None
        self.contact_analyser = ContactAnalyser() if ContactAnalyser else None
        
        self.results = {}
        
        # Log which analysers are available
        available = []
        if self.packing_calc: available.append("packing_metrics")
        if self.interface_analyser:
            available.append(
                "interface_analyser_asu_ec / interface_analyser_lattice_ec"
                if self.interface_metrics == "ec"
                else "interface_analyser_asu_charge / interface_analyser_lattice_charge"
            )
        if self.contact_analyser: available.append("contact_analyser")
        
        print(f"Available analysers: {', '.join(available) if available else 'None'}")
    
    def analyse_single_structure(self, pdb_file, structure_id=None, focus_chains=None, reference_chain_id=None):
        """
        Perform comprehensive analysis on a single crystal structure.
        
        Parameters:
        -----------
        pdb_file : str
            Path to PDB file
        structure_id : str, optional
            Identifier for the structure
        focus_chains : None | list[str] | set[str]
            If provided, chain-focused analysis is applied to stages that support it
            (currently interface analysis and contact analysis). Only chain pairs where
            at least one chain is in this set are analysed.
        reference_chain_id : None | str
            Passed to interface analysis: focal chain for multi-copy lattice metrics
            (isolated vs embedded SASA, burial fraction, cross-chain contact residue fraction).
            
        Returns:
        --------
        dict : Analysis results
        """
        if structure_id is None:
            structure_id = Path(pdb_file).stem
            
        print(f"Analysing structure: {structure_id}")
        
        results = {
            'structure_id': structure_id,
            'pdb_file': pdb_file,
            # Explicit mode metadata to simplify downstream aggregation.
            'interface_context': str(self.interface_context),
            'interface_metrics': str(getattr(self, "interface_metrics", "charge")),
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
            if self.interface_analyser:
                print("  - Analysing interfaces...")
                fn = getattr(self.interface_analyser, "analyse_interfaces", None) or getattr(
                    self.interface_analyser, "analyze_interfaces", None
                )
                if fn is None:
                    interface_data = {"error": "interface analyser has no analyse_interfaces/analyze_interfaces method"}
                else:
                    ref = reference_chain_id if self.interface_context == "lattice" else None
                    interface_data = fn(
                        pdb_file,
                        focus_chains=focus_chains,
                        reference_chain_id=ref,
                    )
                # Normalise metadata inside interface_analysis.summary for consumers
                if isinstance(interface_data, dict):
                    summary = interface_data.get("summary")
                    if not isinstance(summary, dict):
                        summary = {}
                        interface_data["summary"] = summary
                    summary["interface_context"] = str(self.interface_context)
                    summary["interface_metrics"] = str(getattr(self, "interface_metrics", "charge"))
                results['interface_analysis'] = interface_data
            else:
                print("  - Skipping interface analysis (module not available)")
                results['interface_analysis'] = {'error': 'interface_analyser module not available'}
            
            # Phase 3: Contact analysis
            if self.contact_analyser:
                print("  - Analysing crystal contacts...")
                contact_data = self.contact_analyser.analyse_contacts(pdb_file, focus_chains=focus_chains)
                results['contact_analysis'] = contact_data
            else:
                print("  - Skipping contact analysis (module not available)")
                results['contact_analysis'] = {'error': 'contact_analyser module not available'}
            
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
            List of dicts from ``analyse_single_structure``.

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
        description="Crystal packing analysis pipeline (packing metrics, interfaces, contacts).",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Analyse single structure
  python crystal_packing_analyser.py --input model_01.pdb

  # Analyse all PDBs in a directory
  python crystal_packing_analyser.py --input models/

  # Combine per-structure JSON into one file (batch export)
  python crystal_packing_analyser.py --input *.pdb --compare

  # Filter by set(s): one output dir per set (use {} in --output)
  python crystal_packing_analyser.py --input *.pdb --sets set_a set_b --output "crystal_analysis_{}"
  python crystal_packing_analyser.py --input *.pdb --sets set_a set_b --output analysis_output --dry-run
        """
    )

    parser.add_argument(
        '--input', '-i',
        nargs='+',
        required=True,
        help='Input PDB/CIF file(s), directory (all *.pdb/*.cif inside), or glob (e.g. *.pdb).'
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
        help='After processing, write batch_analysis_results.json (all structures in one file).'
    )

    parser.add_argument(
        '--verbose', '-v',
        action='store_true',
        help='Verbose output.'
    )

    parser.add_argument(
        '--dry-run',
        action='store_true',
        help='Print which files would be processed per set and target output directories, then exit.'
    )
    parser.add_argument(
        '--chains',
        metavar='IDS',
        help="Focus analysis on specific chain IDs (comma-separated). Passed through to interface/contact stages that support it. Example: --chains A or --chains A,B",
    )
    parser.add_argument(
        '--reference-chain',
        metavar='ID',
        dest='reference_chain_id',
        help='Focal chain ID for multi-copy assemblies; passed to interface analysis (lattice SASA / burial / contact-residue fraction).',
    )
    parser.add_argument(
        "--interface-metrics",
        choices=("charge", "ec"),
        default="charge",
        help=(
            "Interface metric family. 'charge' reports charge-tag complementarity; "
            "'ec' reports electrostatic complementarity (EC; McCoy method)."
        ),
    )
    parser.add_argument(
        '--interface-context',
        choices=('asu', 'lattice'),
        default='asu',
        help=(
            "Interface context. 'asu' analyses pairwise interfaces in the asymmetric unit; "
            "'lattice' computes lattice-style reference-chain metrics (requires --reference-chain)."
        ),
    )

    args = parser.parse_args()

    focus_chains = None
    if getattr(args, 'chains', None):
        focus_chains = [c.strip() for c in str(args.chains).split(',') if c.strip()]
    reference_chain_id = getattr(args, 'reference_chain_id', None)
    if reference_chain_id:
        reference_chain_id = str(reference_chain_id).strip() or None
    if str(getattr(args, "interface_context", "asu")).strip().lower() == "lattice" and not reference_chain_id:
        print("Error: --interface-context lattice requires --reference-chain.", file=sys.stderr)
        sys.exit(1)

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
        analyser = CrystalPackingAnalyser(
            out_dir,
            interface_context=getattr(args, 'interface_context', 'asu'),
            interface_metrics=getattr(args, "interface_metrics", "charge"),
        )
        all_results = []
        for pdb_file in filtered:
            if not os.path.exists(pdb_file):
                print(f"Warning: File {pdb_file} not found, skipping...", file=sys.stderr)
                continue
            result = analyser.analyse_single_structure(
                pdb_file,
                focus_chains=focus_chains,
                reference_chain_id=reference_chain_id,
            )
            all_results.append(result)

        if args.compare and all_results:
            analyser.compare_structures(all_results)

        print(f"\nSet {label!r} complete. Results saved to: {out_dir}")
        print(f"Number of structures analysed: {len(all_results)}")

if __name__ == "__main__":
    main() 