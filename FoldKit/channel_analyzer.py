#!/usr/bin/env python3
"""
Channel Analyzer
===============

Analyze solvent channels and void spaces in crystal structures.
"""

import glob
import os
import sys
import numpy as np
from pathlib import Path

try:
    from Bio.PDB.PDBParser import PDBParser
    BIOPYTHON_AVAILABLE = True
except ImportError:
    PDBParser = None  # type: ignore[misc, assignment]
    BIOPYTHON_AVAILABLE = False

# Minimum channel size (Å): bottlenecks with radius less than this are discarded.
# In _identify_bottlenecks(), radius = nearest_distance/2 with nearest_distance < 4 Å,
# so radii are always < 2 Å; use a threshold ≤ 2 so the filter can retain some bottlenecks.
MIN_CHANNEL_RADIUS = 1.0


class ChannelAnalyzer:
    """Analyzer for solvent channels and void spaces."""
    
    def __init__(self, probe_radius=1.4):
        """
        Initialize the channel analyzer.
        
        Parameters:
        -----------
        probe_radius : float
            Probe radius for channel analysis (Å)
        """
        if BIOPYTHON_AVAILABLE and PDBParser is not None:
            self.parser = PDBParser(QUIET=True)
        else:
            self.parser = None
            
        self.probe_radius = probe_radius
    
    def analyze_channels(self, pdb_file):
        """
        Analyze solvent channels in a structure.
        
        Parameters:
        -----------
        pdb_file : str
            Path to PDB file
            
        Returns:
        --------
        dict : Channel analysis results
        """
        if not BIOPYTHON_AVAILABLE or self.parser is None:
            return {'error': 'BioPython not available for channel analysis'}
        
        try:
            structure = self.parser.get_structure('crystal', pdb_file)
            
            results = {
                'channel_analysis': {},
                'void_spaces': {},
                'bottleneck_analysis': {}
            }
            
            # Basic void analysis
            void_data = self._analyze_void_spaces(structure)
            results['void_spaces'] = void_data
            
            # Channel connectivity (simplified)
            channel_data = self._analyze_channel_connectivity(structure)
            results['channel_analysis'] = channel_data
            
            # Bottleneck analysis (simplified)
            bottleneck_data = self._analyze_bottlenecks(structure)
            results['bottleneck_analysis'] = bottleneck_data
            
            return results
            
        except Exception as e:
            return {'error': f"Failed to analyze channels: {str(e)}"}
    
    def _analyze_void_spaces(self, structure):
        """Analyze void spaces in the structure."""
        atoms = list(structure.get_atoms())
        
        if not atoms:
            return {'void_fraction': 0, 'accessible_volume': 0}
        
        # Get coordinates and radii
        coords = np.array([atom.coord for atom in atoms])
        
        # Estimate void spaces using grid-based approach
        grid_spacing = 1.0  # Å
        
        # Calculate bounding box
        min_coords = np.min(coords, axis=0) - 5.0
        max_coords = np.max(coords, axis=0) + 5.0
        
        # Create grid
        x_range = np.arange(min_coords[0], max_coords[0], grid_spacing)
        y_range = np.arange(min_coords[1], max_coords[1], grid_spacing)
        z_range = np.arange(min_coords[2], max_coords[2], grid_spacing)
        
        total_grid_points = len(x_range) * len(y_range) * len(z_range)
        accessible_points = 0
        
        # Sample subset of grid points for efficiency
        sample_size = min(10000, total_grid_points)
        sampled_points = self._generate_sample_points(
            min_coords, max_coords, sample_size
        )
        
        # Check accessibility of sampled points
        for point in sampled_points:
            if self._is_point_accessible(point, coords):
                accessible_points += 1
        
        # Estimate void fraction
        void_fraction = accessible_points / len(sampled_points) if len(sampled_points) > 0 else 0
        
        # Estimate accessible volume
        total_volume = np.prod(max_coords - min_coords)
        accessible_volume = total_volume * void_fraction
        
        return {
            'void_fraction': void_fraction,
            'accessible_volume': accessible_volume,
            'total_volume_estimate': total_volume,
            'sampled_points': len(sampled_points),
            'accessible_points': accessible_points
        }
    
    def _generate_sample_points(self, min_coords, max_coords, sample_size):
        """Generate random sample points in the bounding box."""
        points = []
        for _ in range(sample_size):
            x = np.random.uniform(min_coords[0], max_coords[0])
            y = np.random.uniform(min_coords[1], max_coords[1])
            z = np.random.uniform(min_coords[2], max_coords[2])
            points.append([x, y, z])
        return np.array(points)
    
    def _is_point_accessible(self, point, atom_coords, atom_radius=1.7):
        """Check if a point is accessible (not inside any atom)."""
        distances = np.linalg.norm(atom_coords - point, axis=1)
        min_distance = np.min(distances)
        
        # Point is accessible if it's at least probe_radius + atom_radius away
        return min_distance > (self.probe_radius + atom_radius)
    
    def _analyze_channel_connectivity(self, structure):
        """Analyze channel connectivity (simplified approach)."""
        atoms = list(structure.get_atoms())
        
        if len(atoms) < 10:
            return {'connectivity_score': 0, 'channel_count': 0}
        
        # Simple connectivity analysis based on void distribution
        coords = np.array([atom.coord for atom in atoms])
        
        # Calculate local density variations
        density_variations = self._calculate_density_variations(coords)
        
        # Estimate connectivity based on density patterns
        connectivity_score = self._estimate_connectivity(density_variations)
        
        return {
            'connectivity_score': connectivity_score,
            'channel_count': max(1, int(connectivity_score * 5)),  # Rough estimate
            'density_variation_score': np.std(density_variations)
        }
    
    def _calculate_density_variations(self, coords):
        """Calculate local density variations."""
        n_samples = min(100, len(coords))
        sample_indices = np.random.choice(len(coords), n_samples, replace=False)
        
        densities = []
        for idx in sample_indices:
            center = coords[idx]
            
            # Count neighbors within 8Å
            distances = np.linalg.norm(coords - center, axis=1)
            neighbor_count = np.sum(distances <= 8.0)
            
            densities.append(neighbor_count)
        
        return np.array(densities)
    
    def _estimate_connectivity(self, density_variations):
        """Estimate connectivity based on density patterns."""
        if len(density_variations) == 0:
            return 0
        
        # Higher variation suggests more channels/voids
        cv = np.std(density_variations) / np.mean(density_variations) if np.mean(density_variations) > 0 else 0
        
        # Normalize to 0-1 range
        connectivity = min(1.0, cv / 2.0)
        
        return connectivity
    
    def _analyze_bottlenecks(self, structure):
        """Analyze channel bottlenecks (simplified approach)."""
        atoms = list(structure.get_atoms())
        
        if not atoms:
            return {'bottleneck_radius': 0, 'bottleneck_count': 0}
        
        coords = np.array([atom.coord for atom in atoms])
        
        # Find regions with restricted access; keep only channels with radius >= MIN_CHANNEL_RADIUS
        bottlenecks = [
            b for b in self._identify_bottlenecks(coords)
            if b['radius'] >= MIN_CHANNEL_RADIUS
        ]

        if not bottlenecks:
            return {'bottleneck_radius': 0, 'bottleneck_count': 0}
        
        # Calculate bottleneck statistics
        bottleneck_radii = [b['radius'] for b in bottlenecks]
        
        return {
            'bottleneck_count': len(bottlenecks),
            'average_bottleneck_radius': np.mean(bottleneck_radii),
            'minimum_bottleneck_radius': np.min(bottleneck_radii),
            'bottleneck_positions': [b['position'].tolist() for b in bottlenecks]
        }
    
    def _identify_bottlenecks(self, coords):
        """Identify bottleneck positions (simplified)."""
        bottlenecks = []
        if len(coords) < 2:
            return bottlenecks

        # Sample points and find local minima in accessible radius
        n_samples = min(50, len(coords))
        sample_indices = np.random.choice(len(coords), n_samples, replace=False)

        for idx in sample_indices:
            center = coords[idx]

            # Find nearest neighbors (excluding self)
            distances = np.linalg.norm(coords - center, axis=1)
            sorted_d = np.sort(distances)
            nearest_distance = sorted_d[1]  # index 0 is self
            
            # If this point has a very close neighbor, it might be a bottleneck
            if nearest_distance < 4.0:  # Arbitrary threshold
                bottlenecks.append({
                    'position': center,
                    'radius': nearest_distance / 2.0  # Rough estimate
                })
        
        return bottlenecks

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


LIMIT_INLINE = 50


def _run_analysis(analyzer, paths, out_stream, output_file_path=None):
    """Run channel analysis on paths and write results to out_stream."""
    for i, pdb_file in enumerate(paths):
        stem = os.path.splitext(os.path.basename(pdb_file))[0]
        base = (os.path.splitext(output_file_path)[0] + "_" + stem) if output_file_path else stem
        if len(paths) > 1:
            print(f"\n{'='*50}\n[{i+1}/{len(paths)}] {os.path.basename(pdb_file)}\n{'='*50}", file=out_stream)
        else:
            print(f"Analyzing channels in {pdb_file}...", file=out_stream)

        results = analyzer.analyze_channels(pdb_file)

        if 'error' in results:
            print(f"Error: {results['error']}", file=out_stream)
            continue

        print("\nCHANNEL ANALYSIS RESULTS:", file=out_stream)
        print("=" * 40, file=out_stream)
        void_data = results.get('void_spaces') or {}
        if isinstance(void_data, dict):
            print("\nVoid spaces:", file=out_stream)
            print(f"  Void fraction: {void_data.get('void_fraction', 0):.3f}", file=out_stream)
            print(f"  Accessible volume: {void_data.get('accessible_volume', 0):.1f} Å³", file=out_stream)
            print(f"  Total volume (bounding box): {void_data.get('total_volume_estimate', 0):.1f} Å³", file=out_stream)
            print(f"  Sampled points: {void_data.get('sampled_points', 0)} (accessible: {void_data.get('accessible_points', 0)})", file=out_stream)
        channel_data = results.get('channel_analysis') or {}
        if isinstance(channel_data, dict):
            print("\nChannel connectivity:", file=out_stream)
            print(f"  Connectivity score: {channel_data.get('connectivity_score', 0):.3f}", file=out_stream)
            print(f"  Estimated channel count: {channel_data.get('channel_count', 0)}", file=out_stream)
            print(f"  Density variation score: {channel_data.get('density_variation_score', 0):.3f}", file=out_stream)
        bottleneck_data = results.get('bottleneck_analysis') or {}
        if isinstance(bottleneck_data, dict):
            nb = bottleneck_data.get('bottleneck_count', 0)
            print("\nBottlenecks:", file=out_stream)
            print(f"  Count: {nb}", file=out_stream)
            if nb > 0:
                print(f"  Average radius: {bottleneck_data.get('average_bottleneck_radius', 0):.2f} Å", file=out_stream)
                print(f"  Minimum radius: {bottleneck_data.get('minimum_bottleneck_radius', 0):.2f} Å", file=out_stream)
                positions = bottleneck_data.get('bottleneck_positions', [])
                if isinstance(positions, list) and positions:
                    print(f"  Positions (x, y, z Å):", file=out_stream)
                    if len(positions) <= LIMIT_INLINE:
                        for idx, pos in enumerate(positions):
                            if isinstance(pos, (list, tuple)) and len(pos) >= 3:
                                print(f"    {idx+1}: ({pos[0]:.1f}, {pos[1]:.1f}, {pos[2]:.1f})", file=out_stream)
                    else:
                        for idx, pos in enumerate(positions[:LIMIT_INLINE]):
                            if isinstance(pos, (list, tuple)) and len(pos) >= 3:
                                print(f"    {idx+1}: ({pos[0]:.1f}, {pos[1]:.1f}, {pos[2]:.1f})", file=out_stream)
                        extra_path = f"{base}_bottleneck_positions.txt"
                        with open(extra_path, 'w') as ef:
                            ef.write(f"# All {len(positions)} bottleneck positions (x, y, z Å)\n")
                            for idx, pos in enumerate(positions):
                                if isinstance(pos, (list, tuple)) and len(pos) >= 3:
                                    ef.write(f"{idx+1} {pos[0]:.1f} {pos[1]:.1f} {pos[2]:.1f}\n")
                        print(f"  Full list ({len(positions)} positions) written to {extra_path}", file=out_stream)

    if len(paths) > 1:
        print(f"\nDone. Processed {len(paths)} file(s).", file=out_stream)


def main():
    """CLI: single file, multiple files, directory, or glob patterns; optional multi-set filtering and dry-run."""
    import argparse

    parser = argparse.ArgumentParser(
        description="Analyze solvent channels and void spaces in crystal structures.",
        epilog="""Examples:
  python channel_analyzer.py model_01.pdb
  python channel_analyzer.py dir/
  python channel_analyzer.py *.pdb -o results.txt
  python channel_analyzer.py *.pdb --set set_a,set_b -o by_set.txt
  python channel_analyzer.py *.pdb --sets set_a set_b set_c -o "output_{}.txt"
  python channel_analyzer.py *.pdb --per-structure -o "{}_channel.txt"
  python channel_analyzer.py *.pdb --sets set_a set_b -o results.txt --dry-run
""",
    )
    parser.add_argument('input', nargs='+', help='PDB/CIF file(s), directory, or glob pattern (e.g. *.pdb)')
    parser.add_argument(
        '--set', '-s', action='append', dest='sets', metavar='PATTERNS',
        help='Comma-separated patterns; file included only if basename contains ALL. Repeat for multiple sets.',
    )
    parser.add_argument(
        '--sets', dest='sets_multi', nargs='+', metavar='SET',
        help='Multiple set names in one go (one pattern per set). E.g. --sets set_a set_b set_c.',
    )
    parser.add_argument(
        '--output', '-o', metavar='FILE',
        help='Output file. Single file for all sets, or use "{}" for one file per set. Per-structure: use --per-structure with "{}" for file stem.',
    )
    parser.add_argument(
        '--per-structure', '-p', action='store_true',
        help='Write one output file per structure; -o must contain "{}" (replaced by file stem). Can combine with --set/--sets to filter which files are processed.',
    )
    parser.add_argument('--dry-run', action='store_true', help='Print which files would be processed per set and exit.')
    args = parser.parse_args()

    paths = collect_structure_paths(args.input)
    if not paths:
        print("No structure files found.", file=sys.stderr)
        sys.exit(1)

    # Per-structure mode: one output file per input file (optionally filtered by --set/--sets)
    if args.per_structure:
        if not args.output or '{}' not in args.output:
            print("Error: --per-structure requires -o with '{}' in the path (e.g. -o '{}_channel.txt').", file=sys.stderr)
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
        if args.dry_run:
            print("Dry run (per-structure mode): no analysis will be performed.\n", file=sys.stderr)
            for p in paths_to_process:
                stem = os.path.splitext(os.path.basename(p))[0]
                dest = args.output.replace('{}', stem)
                print(f"  {os.path.basename(p)} -> {dest}", file=sys.stderr)
            return
        analyzer = ChannelAnalyzer()
        for p in paths_to_process:
            stem = os.path.splitext(os.path.basename(p))[0]
            out_path = args.output.replace('{}', stem)
            print(f"Writing {os.path.basename(p)} -> {out_path}", file=sys.stderr)
            with open(out_path, 'w') as out:
                _run_analysis(analyzer, [p], out, output_file_path=None)
        return

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
        single_file = args.output and '{}' not in args.output
        for label, patterns, filtered in set_list:
            dest = args.output if single_file else (args.output.replace('{}', label) if args.output and '{}' in args.output else None)
            dest = dest or "stdout"
            print(f"Set {label!r} (patterns: {patterns}): {len(filtered)} file(s) -> {dest}", file=sys.stderr)
            for p in filtered:
                print(f"  {os.path.basename(p)}", file=sys.stderr)
            print(file=sys.stderr)
        return

    single_output_file = args.output and '{}' not in args.output
    out = sys.stdout
    if single_output_file:
        out = open(args.output, 'w')
        print(f"Writing all sets to {args.output}", file=sys.stderr)

    analyzer = ChannelAnalyzer()
    for idx, (label, patterns, filtered) in enumerate(set_list):
        if not filtered:
            print(f"No files match set {label!r} (patterns: {patterns}), skipping.", file=sys.stderr)
            continue

        out_path = None
        current_out = sys.stdout
        if single_output_file:
            current_out = out
        elif args.output:
            out_path = args.output.replace('{}', label)
            current_out = open(out_path, 'w')
            print(f"Writing set {label!r} ({len(filtered)} file(s)) to {out_path}", file=sys.stderr)
        elif len(set_list) > 1:
            print(f"\n=== Set {label!r} ({len(filtered)} file(s)) ===\n", file=sys.stderr)

        try:
            if len(set_list) > 1 and (single_output_file or not args.output):
                print(f"\n{'='*50}\nSet {label!r} (patterns: {patterns})\n{'='*50}", file=current_out)
            _run_analysis(analyzer, filtered, current_out, out_path if out_path is not None else (args.output if single_output_file else None))
        finally:
            if args.output and not single_output_file:
                current_out.close()

    if single_output_file:
        out.close()


if __name__ == "__main__":
    main() 