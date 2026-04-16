#!/usr/bin/env python3
"""
Packing Metrics Calculator
=========================

Calculate basic crystal packing metrics including Matthews coefficient,
solvent content, packing density, and related parameters.
"""

import numpy as np
import warnings
try:
    from Bio.PDB.PDBParser import PDBParser
    from Bio.PDB.PDBIO import PDBIO, Select
    from Bio.PDB.PDBExceptions import PDBConstructionWarning
    warnings.simplefilter('ignore', PDBConstructionWarning)
    BIOPYTHON_AVAILABLE = True
except ImportError:
    PDBParser = None  # type: ignore[misc, assignment]
    PDBIO = None  # type: ignore[misc, assignment]
    Select = None  # type: ignore[misc, assignment]
    PDBConstructionWarning = None  # type: ignore[misc, assignment]
    print("Warning: BioPython not available. Some functionality may be limited.")
    BIOPYTHON_AVAILABLE = False

class PackingMetricsCalculator:
    """Calculator for basic crystal packing metrics."""
    
    def __init__(self):
        """Initialise the calculator."""
        if BIOPYTHON_AVAILABLE and PDBParser is not None:
            self.parser = PDBParser(QUIET=True)
        else:
            self.parser = None
        
        # Atomic volumes for common elements (Å³)
        self.atomic_volumes = {
            'C': 16.44, 'N': 11.51, 'O': 9.13, 'S': 19.86,
            'P': 17.5, 'H': 5.15, 'Ca': 25.99, 'Mg': 22.57,
            'Zn': 14.6, 'Fe': 11.8, 'Mn': 15.6, 'Cu': 11.8,
            'Na': 23.7, 'K': 43.2, 'Cl': 22.7, 'Br': 26.5
        }
        
        # Average atomic weights (Da)
        self.atomic_weights = {
            'C': 12.011, 'N': 14.007, 'O': 15.999, 'S': 32.06,
            'P': 30.974, 'H': 1.008, 'Ca': 40.078, 'Mg': 24.305,
            'Zn': 65.38, 'Fe': 55.845, 'Mn': 54.938, 'Cu': 63.546,
            'Na': 22.990, 'K': 39.098, 'Cl': 35.45, 'Br': 79.904
        }
    
    def calculate_metrics(self, pdb_file):
        """
        Calculate comprehensive packing metrics for a crystal structure.
        
        Parameters:
        -----------
        pdb_file : str
            Path to PDB file
            
        Returns:
        --------
        dict : Dictionary containing all calculated metrics
        """
        if not BIOPYTHON_AVAILABLE or self.parser is None:
            return {'error': 'BioPython not available for structure parsing'}
            
        try:
            structure = self.parser.get_structure('crystal', pdb_file)
            
            # Get unit cell parameters from PDB header
            unit_cell_params = self._extract_unit_cell_params(pdb_file)
            
            results = {}
            
            # Basic unit cell information
            results['unit_cell'] = unit_cell_params
            results['unit_cell_volume'] = self._calculate_unit_cell_volume(unit_cell_params)
            
            # Molecular weight and composition
            molecular_data = self._analyse_molecular_content(structure)
            results.update(molecular_data)
            
            # Matthews coefficient and solvent content
            matthews_data = self._calculate_matthews_metrics(
                results['molecular_weight'], 
                results['unit_cell_volume']
            )
            results.update(matthews_data)
            
            # Packing density metrics
            packing_data = self._calculate_packing_density(structure, results['unit_cell_volume'])
            results.update(packing_data)
            
            # Void volume analysis
            void_data = self._analyse_void_spaces(structure, results['unit_cell_volume'])
            results.update(void_data)
            
            return results
            
        except Exception as e:
            return {'error': f"Failed to calculate packing metrics: {str(e)}"}
    
    def _extract_unit_cell_params(self, pdb_file):
        """Extract unit cell parameters from PDB file header."""
        unit_cell = {'a': 1.0, 'b': 1.0, 'c': 1.0, 'alpha': 90.0, 'beta': 90.0, 'gamma': 90.0, 'space_group': 'P1'}
        
        try:
            with open(pdb_file, 'r') as f:
                for line in f:
                    if line.startswith('CRYST1'):
                        unit_cell['a'] = float(line[6:15].strip())
                        unit_cell['b'] = float(line[15:24].strip())
                        unit_cell['c'] = float(line[24:33].strip())
                        unit_cell['alpha'] = float(line[33:40].strip())
                        unit_cell['beta'] = float(line[40:47].strip())
                        unit_cell['gamma'] = float(line[47:54].strip())
                        unit_cell['space_group'] = line[55:66].strip()
                        break
        except Exception:
            pass
            
        return unit_cell
    
    def _calculate_unit_cell_volume(self, unit_cell):
        """Calculate unit cell volume from parameters."""
        a, b, c = unit_cell['a'], unit_cell['b'], unit_cell['c']
        alpha, beta, gamma = np.radians([unit_cell['alpha'], unit_cell['beta'], unit_cell['gamma']])
        
        volume = a * b * c * np.sqrt(
            1 + 2*np.cos(alpha)*np.cos(beta)*np.cos(gamma) 
            - np.cos(alpha)**2 - np.cos(beta)**2 - np.cos(gamma)**2
        )
        
        return volume
    
    def _analyse_molecular_content(self, structure):
        """Analyse molecular weight and atomic composition."""
        results = {}
        
        total_atoms = 0
        total_weight = 0.0
        element_counts = {}
        chain_count = 0
        residue_count = 0
        
        for model in structure:
            for chain in model:
                chain_count += 1
                for residue in chain:
                    residue_count += 1
                    for atom in residue:
                        element = atom.element.strip().upper()
                        if element == '':
                            element = atom.name[0].upper()
                        
                        total_atoms += 1
                        element_counts[element] = element_counts.get(element, 0) + 1
                        
                        # Add atomic weight
                        if element in self.atomic_weights:
                            total_weight += self.atomic_weights[element]
                        else:
                            # Default to carbon weight for unknown elements
                            total_weight += self.atomic_weights['C']
        
        results['total_atoms'] = total_atoms
        results['molecular_weight'] = total_weight
        results['element_composition'] = element_counts
        results['chain_count'] = chain_count
        results['residue_count'] = residue_count
        
        return results
    
    def _calculate_matthews_metrics(self, molecular_weight, unit_cell_volume):
        """Calculate Matthews coefficient and solvent content."""
        results = {}
        
        # Matthews coefficient (Å³/Da)
        matthews_coeff = unit_cell_volume / molecular_weight if molecular_weight > 0 else 0
        results['matthews_coefficient'] = matthews_coeff
        
        # Solvent content estimation using Matthews relationship
        # Typical protein density ~1.35 g/cm³ = 0.81 Da/Å³
        protein_volume = molecular_weight / 0.81  # Da / (Da/Å³) = Å³
        solvent_volume = unit_cell_volume - protein_volume
        solvent_content = solvent_volume / unit_cell_volume * 100 if unit_cell_volume > 0 else 0
        
        results['estimated_protein_volume'] = protein_volume
        results['estimated_solvent_volume'] = solvent_volume
        results['solvent_content_percent'] = max(0, min(100, solvent_content))
        
        return results
    
    def _calculate_packing_density(self, structure, unit_cell_volume):
        """Calculate packing density based on atomic volumes."""
        results = {}
        
        total_atomic_volume = 0.0
        
        for model in structure:
            for chain in model:
                for residue in chain:
                    for atom in residue:
                        element = atom.element.strip().upper()
                        if element == '':
                            element = atom.name[0].upper()
                        
                        if element in self.atomic_volumes:
                            total_atomic_volume += self.atomic_volumes[element]
                        else:
                            # Default to carbon volume
                            total_atomic_volume += self.atomic_volumes['C']
        
        packing_density = total_atomic_volume / unit_cell_volume if unit_cell_volume > 0 else 0
        packing_efficiency = packing_density * 100
        
        results['total_atomic_volume'] = total_atomic_volume
        results['packing_density'] = packing_density
        results['packing_efficiency_percent'] = packing_efficiency
        results['void_fraction'] = 1 - packing_density
        
        return results
    
    def _analyse_void_spaces(self, structure, unit_cell_volume):
        """Analyse void spaces and their distribution."""
        results = {}
        
        # Get all atom coordinates
        coords = []
        radii = []
        
        for model in structure:
            for chain in model:
                for residue in chain:
                    for atom in residue:
                        coords.append(atom.coord)
                        
                        # Estimate atomic radius from volume
                        element = atom.element.strip().upper()
                        if element == '':
                            element = atom.name[0].upper()
                        
                        if element in self.atomic_volumes:
                            volume = self.atomic_volumes[element]
                        else:
                            volume = self.atomic_volumes['C']
                        
                        # Radius from volume (assuming sphere: V = 4/3 * π * r³)
                        radius = (3 * volume / (4 * np.pi)) ** (1/3)
                        radii.append(radius)
        
        coords = np.array(coords)
        radii = np.array(radii)
        
        # Basic void analysis
        if len(coords) > 0:
            # Calculate average nearest-neighbour distance
            try:
                from scipy.spatial.distance import pdist
                distances = pdist(coords)
                avg_nn_distance = np.mean(distances) if len(distances) > 0 else 0
                min_distance = np.min(distances) if len(distances) > 0 else 0
            except ImportError:
                # Fallback if scipy not available
                if len(coords) > 1:
                    # Simple pairwise distance calculation for small sets
                    distances = []
                    for i in range(len(coords)):
                        for j in range(i+1, len(coords)):
                            dist = np.linalg.norm(coords[i] - coords[j])
                            distances.append(dist)
                    avg_nn_distance = np.mean(distances) if distances else 0
                    min_distance = np.min(distances) if distances else 0
                else:
                    avg_nn_distance = 0
                    min_distance = 0
            
            results['average_interatomic_distance'] = avg_nn_distance
            results['minimum_interatomic_distance'] = min_distance
            results['atom_density'] = len(coords) / unit_cell_volume if unit_cell_volume > 0 else 0
        else:
            results['average_interatomic_distance'] = 0
            results['minimum_interatomic_distance'] = 0
            results['atom_density'] = 0
        
        return results


def filter_paths_by_patterns(paths, patterns):
    """
    Return paths whose basename contains every pattern in `patterns`.
    Patterns are matched as substrings (e.g. 'set_a' matches 'prefix_set_a.pdb').
    """
    import os
    if not patterns:
        return list(paths)
    result = []
    for p in paths:
        name = os.path.basename(p)
        if all(pat in name for pat in patterns):
            result.append(p)
    return result


def collect_structure_paths(inputs):
    """
    Expand inputs (files, directories, glob patterns) to a list of structure file paths.
    Accepts: single files, multiple files, a directory (all *.pdb, *.cif inside), or globs (e.g. *.pdb).
    """
    import glob
    import os
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


def main():
    """CLI: single file, multiple files, directory, or glob patterns; optional multi-set filtering."""
    import argparse
    import os
    import sys

    parser = argparse.ArgumentParser(
        description="Calculate crystal packing metrics (Matthews coefficient, solvent content, etc.).",
        epilog="""Examples:
  python packing_metrics.py model_01.pdb
  python packing_metrics.py dir/
  python packing_metrics.py *.pdb -o results.txt
  python packing_metrics.py *.pdb --sets set_a set_b set_c -o "analysis_{}.txt"
  python packing_metrics.py *.pdb --per-structure -o "{}_metrics.txt"
  python packing_metrics.py *.pdb --per-structure --sets set_a set_b -o "{}_metrics.txt"
""",
    )
    parser.add_argument(
        'input',
        nargs='+',
        help='PDB/CIF file(s), directory (all *.pdb/*.cif inside), or glob pattern (e.g. *.pdb).',
    )
    parser.add_argument(
        '--quiet', '-q',
        action='store_true',
        help='One-line summary per file only.',
    )
    parser.add_argument(
        '--output', '-o',
        metavar='FILE',
        help='Output file. Single file for all sets, or use "{}" for one file per set. Per-structure: use --per-structure with "{}" for file stem.',
    )
    parser.add_argument(
        '--per-structure', '-p',
        action='store_true',
        help='Write one output file per structure; -o must contain "{}" (replaced by file stem). Can combine with --set/--sets to filter which files are processed.',
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
    from cli_log import add_log_args, setup_log_from_args
    add_log_args(parser)
    args = parser.parse_args()
    setup_log_from_args(args, script_path=__file__, inputs=list(getattr(args, "input", []) or []), pattern=None)

    paths = collect_structure_paths(args.input)
    if not paths:
        print("No structure files found.", file=sys.stderr)
        sys.exit(1)

    # Per-structure mode: one output file per input file (optionally filtered by --set/--sets)
    if args.per_structure:
        if not args.output or '{}' not in args.output:
            print("Error: --per-structure requires -o with '{}' in the path (e.g. -o '{}_metrics.txt').", file=sys.stderr)
            sys.exit(1)
        paths_to_process = list(paths)
        if getattr(args, 'sets_multi', None) or args.sets:
            filtered_set = set()
            if getattr(args, 'sets_multi', None):
                for name in args.sets_multi:
                    patterns = [name.strip()] if name.strip() else []
                    if patterns:
                        filtered_set.update(filter_paths_by_patterns(paths, patterns))
            if args.sets:
                for s in args.sets:
                    patterns = [p.strip() for p in s.split(',') if p.strip()]
                    if patterns:
                        filtered_set.update(filter_paths_by_patterns(paths, patterns))
            paths_to_process = sorted(filtered_set)
        if not paths_to_process:
            print("No structure files match the filter.", file=sys.stderr)
            sys.exit(1)
        for p in paths_to_process:
            stem = os.path.splitext(os.path.basename(p))[0]
            out_path = args.output.replace('{}', stem)
            print(f"Writing {os.path.basename(p)} -> {out_path}", file=sys.stderr)
            with open(out_path, 'w') as out_stream:
                _run_analysis([p], args.quiet, out_stream)
        return

    # Build list of (label, filtered_paths). If --set/--sets not used, one set with all paths.
    set_list = []
    if getattr(args, 'sets_multi', None):
        for name in args.sets_multi:
            patterns = [name.strip()] if name.strip() else []
            if patterns:
                label = patterns[0]
                filtered = filter_paths_by_patterns(paths, patterns)
                set_list.append((label, filtered))
    if args.sets:
        for s in args.sets:
            patterns = [p.strip() for p in s.split(',') if p.strip()]
            if patterns:
                label = '_'.join(patterns)
                filtered = filter_paths_by_patterns(paths, patterns)
                set_list.append((label, filtered))
    if not set_list:
        set_list = [("all", list(paths))]

    single_output_file = args.output and '{}' not in args.output
    if single_output_file and args.output:
        out_stream = open(args.output, 'w')
        print(f"Writing all sets to {args.output}", file=sys.stderr)
    else:
        out_stream = None

    try:
        for label, filtered in set_list:
            if not filtered:
                print(f"No files match set {label!r}, skipping.", file=sys.stderr)
                continue
            if args.output and not single_output_file:
                out_path = args.output.replace('{}', label)
                out_stream = open(out_path, 'w')
                print(f"Writing set {label!r} ({len(filtered)} file(s)) to {out_path}", file=sys.stderr)
            elif not single_output_file:
                out_stream = sys.stdout
                if len(set_list) > 1:
                    print(f"\n=== Set {label!r} ({len(filtered)} file(s)) ===\n", file=sys.stderr)
            if single_output_file and len(set_list) > 1:
                print(f"\n{'='*50}\nSet {label!r}\n{'='*50}", file=out_stream)

            try:
                _run_analysis(filtered, args.quiet, out_stream)
            finally:
                if args.output and not single_output_file:
                    out_stream.close()

        if single_output_file and out_stream:
            out_stream.close()
    except Exception:
        if single_output_file and out_stream:
            out_stream.close()
        raise


def _run_analysis(paths, quiet, out_stream):
    """Run packing metrics for paths and print to out_stream (file or stdout)."""
    import os
    import sys
    calculator = PackingMetricsCalculator()
    for i, pdb_file in enumerate(paths):
        if not quiet and len(paths) > 1:
            print(f"\n{'='*50}\n[{i+1}/{len(paths)}] {os.path.basename(pdb_file)}\n{'='*50}", file=out_stream)
        elif not quiet:
            print(f"Analysing {pdb_file}...", file=out_stream)

        results = calculator.calculate_metrics(pdb_file)

        if 'error' in results:
            print(f"Error: {results['error']}", file=out_stream)
            continue

        if quiet:
            vm = results.get('matthews_coefficient', 'N/A')
            solv = results.get('solvent_content_percent', 'N/A')
            print(f"{os.path.basename(pdb_file)}\tVm={vm}\tsolvent%={solv}", file=out_stream)
        else:
            print("\nPACKING METRICS RESULTS:", file=out_stream)
            print("=" * 30, file=out_stream)
            for key, value in results.items():
                if isinstance(value, dict):
                    print(f"{key}:", file=out_stream)
                    for subkey, subvalue in value.items():
                        print(f"  {subkey}: {subvalue}", file=out_stream)
                else:
                    print(f"{key}: {value}", file=out_stream)

    if len(paths) > 1 and not quiet:
        print(f"\nDone. Processed {len(paths)} file(s).", file=out_stream)


if __name__ == "__main__":
    main() 