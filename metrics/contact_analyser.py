#!/usr/bin/env python3
"""
Contact Analyser
===============

Analyse crystal contacts and symmetry-related interactions
in protein crystal structures.
"""

import glob
import os
import sys
import numpy as np
import warnings
from pathlib import Path

_REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if _REPO_ROOT not in sys.path:
    sys.path.insert(0, _REPO_ROOT)

from cli_log import add_log_args, setup_log_from_args

# Distance limits (Å) for contact types - contacts outside these are ignored
H_BOND_MAX = 3.5
ELECTROSTATIC_MAX = 5.0
HYDROPHOBIC_MIN, HYDROPHOBIC_MAX = 3.5, 4.5
VDW_MAX = 5.0
# Initial distance filter: keep pairs within this so type-specific limits can be applied
MAX_CONTACT_DISTANCE = max(ELECTROSTATIC_MAX, VDW_MAX, HYDROPHOBIC_MAX)


def _contact_within_limit(contact_type: str, distance: float) -> bool:
    """True if distance is within the allowed limit for this contact type."""
    if contact_type == 'hydrogen_bond':
        return distance <= H_BOND_MAX
    if contact_type == 'electrostatic':
        return distance <= ELECTROSTATIC_MAX
    if contact_type == 'hydrophobic':
        return HYDROPHOBIC_MIN <= distance <= HYDROPHOBIC_MAX
    if contact_type == 'van_der_waals':
        return distance <= VDW_MAX
    return False


try:
    from Bio.PDB.PDBParser import PDBParser
    from Bio.PDB.PDBExceptions import PDBConstructionWarning
    warnings.simplefilter('ignore', PDBConstructionWarning)
    BIOPYTHON_AVAILABLE = True
except ImportError:
    PDBParser = None  # type: ignore[misc, assignment]
    PDBConstructionWarning = None  # type: ignore[misc, assignment]
    BIOPYTHON_AVAILABLE = False

class ContactAnalyser:
    """Analyser for crystal contacts and symmetry interactions."""
    
    def __init__(self, contact_distance=4.5):
        """
        Initialise the contact analyser.
        
        Parameters:
        -----------
        contact_distance : float
            Maximum distance for considering atoms in contact (Å)
        """
        if BIOPYTHON_AVAILABLE and PDBParser is not None:
            self.parser = PDBParser(QUIET=True)
        else:
            self.parser = None
            
        self.contact_distance = contact_distance
    
    def analyse_contacts(self, pdb_file, focus_chains=None):
        """
        Analyse crystal contacts in a structure.
        
        Parameters:
        -----------
        pdb_file : str
            Path to PDB file
            
        Returns:
        --------
        dict : Contact analysis results
        """
        if not BIOPYTHON_AVAILABLE or self.parser is None:
            return {'error': 'BioPython not available for contact analysis'}
        
        try:
            structure = self.parser.get_structure('crystal', pdb_file)
            
            # Extract unit cell and space group information
            unit_cell_info = self._extract_crystal_info(pdb_file)
            
            results = {
                'unit_cell_info': unit_cell_info,
                'crystal_contacts': [],
                'contact_summary': {}
            }
            
            # Analyse contacts within the asymmetric unit
            asu_contacts = self._analyse_asu_contacts(structure, focus_chains=focus_chains)
            results['asu_contacts'] = asu_contacts
            
            # Analyse potential crystal contacts (simplified)
            # Note: Full crystal contact analysis requires symmetry operations
            crystal_contacts = self._estimate_crystal_contacts(structure, unit_cell_info)
            results['crystal_contacts'] = crystal_contacts
            
            # Calculate summary statistics
            summary = self._calculate_contact_summary(asu_contacts, crystal_contacts)
            results['contact_summary'] = summary
            
            return results
            
        except Exception as e:
            return {'error': f"Failed to analyse contacts: {str(e)}"}
    
    def _extract_crystal_info(self, pdb_file):
        """Extract crystallographic information from PDB file."""
        crystal_info = {
            'space_group': 'P1',
            'unit_cell': {'a': 1.0, 'b': 1.0, 'c': 1.0, 'alpha': 90.0, 'beta': 90.0, 'gamma': 90.0}
        }
        
        try:
            with open(pdb_file, 'r') as f:
                for line in f:
                    if line.startswith('CRYST1'):
                        crystal_info['unit_cell'] = {
                            'a': float(line[6:15].strip()),
                            'b': float(line[15:24].strip()),
                            'c': float(line[24:33].strip()),
                            'alpha': float(line[33:40].strip()),
                            'beta': float(line[40:47].strip()),
                            'gamma': float(line[47:54].strip())
                        }
                        crystal_info['space_group'] = line[55:66].strip()
                        break
        except Exception:
            pass
            
        return crystal_info
    
    def _analyse_asu_contacts(self, structure, focus_chains=None):
        """Analyse contacts within the asymmetric unit."""
        contacts = []
        
        # Get all chains in the structure
        chains = list(structure.get_chains())

        focus_set = None
        if focus_chains:
            focus_set = {str(c).strip() for c in focus_chains if str(c).strip()}
        
        if len(chains) < 2:
            return {
                'intra_asu_contacts': [],
                'contact_count': 0,
                'contact_density': 0
            }
        
        # Analyse contacts between chains
        for i, chain1 in enumerate(chains):
            for j, chain2 in enumerate(chains):
                if i >= j:  # Avoid duplicates and self-comparison
                    continue
                if focus_set is not None and (chain1.id not in focus_set and chain2.id not in focus_set):
                    continue
                
                chain_contacts = self._find_chain_contacts(chain1, chain2)
                contacts.extend(chain_contacts)

        # Denominator must cover every chain that contributes atoms to the counted
        # contacts. With --chains, pairs are (focus, partner); using only focus
        # chains in the denominator inflates density by omitting partner atoms.
        chain_ids_in_contacts: set[str] = set()
        for c in contacts:
            c1, c2 = c.get('chain1'), c.get('chain2')
            if c1 is not None:
                chain_ids_in_contacts.add(str(c1))
            if c2 is not None:
                chain_ids_in_contacts.add(str(c2))
        if focus_set is not None and chain_ids_in_contacts:
            chain_by_id = {c.id: c for c in chains}
            chains_for_density = [
                chain_by_id[cid]
                for cid in chain_ids_in_contacts
                if cid in chain_by_id
            ]
            # Stable order: follow original structure chain order where possible
            order = {c.id: i for i, c in enumerate(chains)}
            chains_for_density.sort(key=lambda c: order.get(c.id, 9999))
        else:
            chains_for_density = chains

        return {
            'intra_asu_contacts': contacts,
            'contact_count': len(contacts),
            'contact_density': self._calculate_contact_density(contacts, chains_for_density)
        }
    
    def _find_chain_contacts(self, chain1, chain2):
        """Find contacts between two chains."""
        contacts = []
        
        atoms1 = list(chain1.get_atoms())
        atoms2 = list(chain2.get_atoms())
        
        for atom1 in atoms1:
            for atom2 in atoms2:
                distance = np.linalg.norm(atom1.coord - atom2.coord)
                if distance > self.contact_distance:
                    continue
                contact_type = self._classify_contact_type(atom1, atom2)
                if not _contact_within_limit(contact_type, float(distance)):
                    continue
                contact = {
                    'chain1': chain1.id,
                    'chain2': chain2.id,
                    'residue1': f"{atom1.parent.resname}{atom1.parent.id[1]}",
                    'residue2': f"{atom2.parent.resname}{atom2.parent.id[1]}",
                    'atom1': atom1.name,
                    'atom2': atom2.name,
                    'distance': distance,
                    'contact_type': contact_type
                }
                contacts.append(contact)

        return contacts
    
    def _classify_contact_type(self, atom1, atom2):
        """Classify the type of contact based on atom types."""
        element1 = atom1.element.strip().upper() or atom1.name[0].upper()
        element2 = atom2.element.strip().upper() or atom2.name[0].upper()
        
        # Hydrogen bond criteria
        h_bond_elements = {'N', 'O', 'S'}
        if (element1 in h_bond_elements and element2 in h_bond_elements):
            return 'hydrogen_bond'
        
        # Hydrophobic contact
        hydrophobic_elements = {'C'}
        if element1 in hydrophobic_elements and element2 in hydrophobic_elements:
            return 'hydrophobic'
        
        # Electrostatic (simplified)
        charged_elements = {'N', 'O', 'S', 'P'}
        if element1 in charged_elements or element2 in charged_elements:
            return 'electrostatic'
        
        return 'van_der_waals'
    
    def _calculate_contact_density(self, contacts, chains):
        """Calculate contact density."""
        if not contacts or not chains:
            return 0
        
        # Calculate total atoms in all chains
        total_atoms = sum(len(list(chain.get_atoms())) for chain in chains)
        
        return len(contacts) / total_atoms if total_atoms > 0 else 0
    
    def _estimate_crystal_contacts(self, structure, unit_cell_info):
        """
        Estimate crystal contacts (simplified approach).
        
        Note: This is a simplified estimation. Proper crystal contact
        analysis requires applying crystallographic symmetry operations.
        """
        # For now, return a placeholder that would be expanded with
        # proper symmetry operation implementation
        
        chains = list(structure.get_chains())
        surface_atoms = self._identify_surface_atoms(structure)
        
        # Estimate potential crystal contacts based on surface atoms
        potential_contacts = len(surface_atoms)
        
        # Rough estimation of crystal contacts
        unit_cell = unit_cell_info['unit_cell']
        cell_volume = self._calculate_cell_volume(unit_cell)
        
        return {
            'potential_crystal_contacts': potential_contacts,
            'surface_atoms': len(surface_atoms),
            'estimated_contact_area': potential_contacts * 4.0,  # ~4 Å² per contact
            'packing_efficiency': self._estimate_packing_efficiency(structure, cell_volume)
        }
    
    def _identify_surface_atoms(self, structure):
        """Identify surface atoms (simplified approach)."""
        surface_atoms = []
        
        # Get all atoms
        all_atoms = list(structure.get_atoms())
        
        # Simple surface identification: atoms with fewer neighbours
        for atom in all_atoms:
            neighbour_count = 0
            for other_atom in all_atoms:
                if atom != other_atom:
                    distance = np.linalg.norm(atom.coord - other_atom.coord)
                    if distance <= 5.0:  # Neighbour cutoff
                        neighbour_count += 1
            
            # Consider atoms with < 12 neighbours as surface atoms
            if neighbour_count < 12:
                surface_atoms.append(atom)
        
        return surface_atoms
    
    def _calculate_cell_volume(self, unit_cell):
        """Calculate unit cell volume."""
        a, b, c = unit_cell['a'], unit_cell['b'], unit_cell['c']
        alpha = np.radians(unit_cell['alpha'])
        beta = np.radians(unit_cell['beta'])
        gamma = np.radians(unit_cell['gamma'])
        
        volume = a * b * c * np.sqrt(
            1 + 2*np.cos(alpha)*np.cos(beta)*np.cos(gamma) 
            - np.cos(alpha)**2 - np.cos(beta)**2 - np.cos(gamma)**2
        )
        
        return volume
    
    def _estimate_packing_efficiency(self, structure, cell_volume):
        """Estimate crystal packing efficiency."""
        # Calculate molecular volume (simplified)
        atom_count = len(list(structure.get_atoms()))
        
        # Rough atomic volume estimation
        avg_atomic_volume = 20.0  # Å³ per atom (rough average)
        molecular_volume = atom_count * avg_atomic_volume
        
        return molecular_volume / cell_volume if cell_volume > 0 else 0
    
    def _calculate_contact_summary(self, asu_contacts, crystal_contacts):
        """Calculate summary statistics for contacts."""
        summary = {}
        
        # ASU contact summary
        asu_contact_list = asu_contacts.get('intra_asu_contacts', [])
        
        if asu_contact_list:
            distances = [contact['distance'] for contact in asu_contact_list]
            contact_types = [contact['contact_type'] for contact in asu_contact_list]
            
            summary['asu_stats'] = {
                'total_contacts': len(asu_contact_list),
                'average_distance': np.mean(distances),
                'min_distance': np.min(distances),
                'max_distance': np.max(distances),
                'contact_type_distribution': self._count_contact_types(contact_types)
            }
        else:
            summary['asu_stats'] = {
                'total_contacts': 0,
                'average_distance': 0,
                'min_distance': 0,
                'max_distance': 0,
                'contact_type_distribution': {}
            }
        
        # Crystal contact summary
        summary['crystal_stats'] = {
            'potential_contacts': crystal_contacts.get('potential_crystal_contacts', 0),
            'surface_atoms': crystal_contacts.get('surface_atoms', 0),
            'estimated_contact_area': crystal_contacts.get('estimated_contact_area', 0),
            'packing_efficiency': crystal_contacts.get('packing_efficiency', 0)
        }
        
        return summary
    
    def _count_contact_types(self, contact_types):
        """Count occurrences of each contact type."""
        type_counts = {}
        for contact_type in contact_types:
            type_counts[contact_type] = type_counts.get(contact_type, 0) + 1
        return type_counts

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


LIMIT_INLINE = 200


def _run_analysis(analyser, paths, out_stream, output_file_path=None, focus_chains=None):
    """Run contact analysis on paths and write results to out_stream."""
    for i, pdb_file in enumerate(paths):
        stem = os.path.splitext(os.path.basename(pdb_file))[0]
        base = (os.path.splitext(output_file_path)[0] + "_" + stem) if output_file_path else stem
        if len(paths) > 1:
            print(f"\n{'='*50}\n[{i+1}/{len(paths)}] {os.path.basename(pdb_file)}\n{'='*50}", file=out_stream)
        else:
            print(f"Analysing contacts in {pdb_file}...", file=out_stream)

        results = analyser.analyse_contacts(pdb_file, focus_chains=focus_chains)

        if 'error' in results:
            print(f"Error: {results['error']}", file=out_stream)
            continue

        print("\nCONTACT ANALYSIS RESULTS:", file=out_stream)
        print("=" * 40, file=out_stream)
        crystal_info = results.get('unit_cell_info')
        if isinstance(crystal_info, dict):
            print(f"Space group: {crystal_info.get('space_group', '')}", file=out_stream)
            unit_cell = crystal_info.get('unit_cell')
            if isinstance(unit_cell, dict):
                print(f"Unit cell: a={unit_cell.get('a', 0):.2f}, b={unit_cell.get('b', 0):.2f}, c={unit_cell.get('c', 0):.2f}", file=out_stream)
        summary = results.get('contact_summary')
        if isinstance(summary, dict):
            asu_stats = summary.get('asu_stats', {}) or {}
            crystal_stats = summary.get('crystal_stats', {}) or {}
            if isinstance(asu_stats, dict):
                print(f"\nASU contacts: {asu_stats.get('total_contacts', 0)}", file=out_stream)
                print(f"  Distance (Å): min={asu_stats.get('min_distance', 0):.2f} max={asu_stats.get('max_distance', 0):.2f} avg={asu_stats.get('average_distance', 0):.2f}", file=out_stream)
                type_dist = asu_stats.get('contact_type_distribution') or {}
                if isinstance(type_dist, dict) and type_dist:
                    print(f"  Contact types: {type_dist}", file=out_stream)
            if isinstance(crystal_stats, dict):
                print(f"Surface atoms: {crystal_stats.get('surface_atoms', 0)}", file=out_stream)
                print(f"Estimated contact area: {crystal_stats.get('estimated_contact_area', 0):.1f} Å²", file=out_stream)
                print(f"Packing efficiency: {crystal_stats.get('packing_efficiency', 0):.3f}", file=out_stream)
        asu_contacts = results.get('asu_contacts')
        if isinstance(asu_contacts, dict):
            intra = asu_contacts.get('intra_asu_contacts', [])
            if isinstance(intra, list) and intra:
                print(f"\nContact details (chain1 chain2 res1 atom1 res2 atom2 distance type):", file=out_stream)
                if len(intra) <= LIMIT_INLINE:
                    for c in intra:
                        if isinstance(c, dict):
                            print(f"  {c.get('chain1','')} {c.get('chain2','')} {c.get('residue1','')} {c.get('atom1','')} {c.get('residue2','')} {c.get('atom2','')} {c.get('distance',0):.2f} Å {c.get('contact_type','')}", file=out_stream)
                else:
                    for c in intra[:LIMIT_INLINE]:
                        if isinstance(c, dict):
                            print(f"  {c.get('chain1','')} {c.get('chain2','')} {c.get('residue1','')} {c.get('atom1','')} {c.get('residue2','')} {c.get('atom2','')} {c.get('distance',0):.2f} Å {c.get('contact_type','')}", file=out_stream)
                    extra_path = f"{base}_asu_contacts.txt"
                    with open(extra_path, 'w') as ef:
                        ef.write(f"# All {len(intra)} ASU contacts: chain1 chain2 res1 atom1 res2 atom2 distance type\n")
                        for c in intra:
                            if isinstance(c, dict):
                                ef.write(f"{c.get('chain1','')} {c.get('chain2','')} {c.get('residue1','')} {c.get('atom1','')} {c.get('residue2','')} {c.get('atom2','')} {c.get('distance',0):.2f} {c.get('contact_type','')}\n")
                    print(f"  Full list ({len(intra)} contacts) written to {extra_path}", file=out_stream)

    if len(paths) > 1:
        print(f"\nDone. Processed {len(paths)} file(s).", file=out_stream)


def main():
    """CLI: single file, multiple files, directory, or glob patterns; optional multi-set filtering and dry-run."""
    import argparse

    parser = argparse.ArgumentParser(
        description="Analyse crystal contacts and symmetry-related interactions.",
        epilog="""Examples:
  python contact_analyser.py model_01.pdb
  python contact_analyser.py dir/
  python contact_analyser.py *.pdb -o contact_results.txt
  python contact_analyser.py *.pdb --set set_a,set_b -o by_set.txt
  python contact_analyser.py *.pdb --sets set_a set_b set_c -o "contact_{}.txt"
  python contact_analyser.py *.pdb --per-structure -o "{}_contact.txt"
  python contact_analyser.py *.pdb --per-structure --sets set_a set_b -o "{}_contact.txt"
  python contact_analyser.py *.pdb --sets set_a set_b -o contact_results.txt --dry-run

Filter text output to CSV by PDB and/or chain: contact_molecule_report_csv.py
  python contact_molecule_report_csv.py contact_results.txt -m A -m B -o contacts_AB.csv
  python contact_molecule_report_csv.py contact_results.txt --pdbs model_01.pdb --chains A,B --output-dir ./out
  python contact_molecule_report_csv.py contact_results_model_01_asu_contacts.txt --structure-basename model_01.pdb -m A -o model_01_contacts.csv
""",
    )
    parser.add_argument('input', nargs='+', help='PDB/CIF file(s), directory, or glob pattern (e.g. *.pdb).')
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
    parser.add_argument(
        '--chains',
        metavar='IDS',
        help="Focus analysis on specific chain IDs (comma-separated). Only chain pairs where at least one chain is in this list are analysed. Example: --chains A or --chains A,B",
    )
    parser.add_argument('--dry-run', action='store_true', help='Print which files would be processed per set, then exit.')
    add_log_args(parser)
    args = parser.parse_args()
    summary_log = setup_log_from_args(
        args,
        script_path=__file__,
        inputs=list(getattr(args, "input", []) or []),
        pattern=None,
    )
    if summary_log is not None:
        summary_log.task("Contact analysis (contact_analyser.py)")
        summary_log.kv("per_structure", bool(getattr(args, "per_structure", False)))
        summary_log.kv("output", getattr(args, "output", None) or "stdout")
        summary_log.kv("dry_run", bool(getattr(args, "dry_run", False)))

    focus_chains = None
    if getattr(args, 'chains', None):
        focus_chains = [c.strip() for c in str(args.chains).split(',') if c.strip()]

    paths = collect_structure_paths(args.input)
    if not paths:
        print("No structure files found.", file=sys.stderr)
        if summary_log is not None:
            summary_log.error("No structure files found.")
        sys.exit(1)
    if summary_log is not None:
        summary_log.kv("n_inputs", len(paths))

    # Per-structure mode: one output file per input file (optionally filtered by --set/--sets)
    if args.per_structure:
        if not args.output or '{}' not in args.output:
            print("Error: --per-structure requires -o with '{}' in the path (e.g. -o '{}_contact.txt').", file=sys.stderr)
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
            if summary_log is not None:
                summary_log.error("No structure files match the filter.")
            sys.exit(1)
        if args.dry_run:
            print("Dry run (per-structure mode): no analysis will be performed.\n", file=sys.stderr)
            for p in paths_to_process:
                stem = os.path.splitext(os.path.basename(p))[0]
                dest = args.output.replace('{}', stem)
                print(f"  {os.path.basename(p)} -> {dest}", file=sys.stderr)
            return
        analyser = ContactAnalyser()
        for p in paths_to_process:
            stem = os.path.splitext(os.path.basename(p))[0]
            out_path = args.output.replace('{}', stem)
            print(f"Writing {os.path.basename(p)} -> {out_path}", file=sys.stderr)
            if summary_log is not None:
                summary_log.task(f"Analyse contacts: {os.path.basename(p)}")
            with open(out_path, 'w') as out:
                _run_analysis(analyser, [p], out, output_file_path=None, focus_chains=focus_chains)
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

    analyser = ContactAnalyser()
    for idx, (label, patterns, filtered) in enumerate(set_list):
        if not filtered:
            print(f"No files match set {label!r} (patterns: {patterns}), skipping.", file=sys.stderr)
            continue
        if summary_log is not None:
            summary_log.task(f"Set {label!r}: {len(filtered)} structure(s)")

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
            _run_analysis(analyser, filtered, current_out, out_path if out_path is not None else (args.output if single_output_file else None), focus_chains=focus_chains)
        finally:
            if args.output and not single_output_file:
                current_out.close()

    if single_output_file:
        out.close()


if __name__ == "__main__":
    main() 