#!/usr/bin/env python3
"""
Interface Analyzer
=================

Analyze interfaces between molecules in crystal structures using the approach
described in Biopython's Interface Analysis wiki:
https://biopython.org/wiki/Interface_Analysis

- **Contact residues** (contact_residues_chain1/2): First identified by NeighborSearch (residue pairs within
  distance threshold, different chains). Then filtered by applying the same criteria as the contact list:
  keep only residues that have at least one atom pair (with the other chain) satisfying both the distance
  cutoff and the type-specific distance limits (H-bond, hydrophobic, etc.).
- **Buried surface area:** Shrake-Rupley SASA when available:
  BSA = SASA(chain1 alone) + SASA(chain2 alone) - SASA(complex).
  Falls back to a per-atom estimate if SASA is not available.
"""

import copy
import glob
import os
import sys
import numpy as np
import warnings
from pathlib import Path

try:
    from Bio.PDB.PDBParser import PDBParser
    from Bio.PDB.NeighborSearch import NeighborSearch
    from Bio.PDB.PDBExceptions import PDBConstructionWarning
    from Bio.PDB import Structure, Model, Superimposer
    warnings.simplefilter('ignore', PDBConstructionWarning)
    BIOPYTHON_AVAILABLE = True
except ImportError:
    PDBParser = None  # type: ignore[misc, assignment]
    NeighborSearch = None  # type: ignore[misc, assignment]
    PDBConstructionWarning = None  # type: ignore[misc, assignment]
    Structure = None  # type: ignore[misc, assignment]
    Model = None  # type: ignore[misc, assignment]
    Superimposer = None  # type: ignore[misc, assignment]
    BIOPYTHON_AVAILABLE = False

try:
    from Bio.PDB.SASA import ShrakeRupley
    SASA_AVAILABLE = True
except ImportError:
    ShrakeRupley = None  # type: ignore[misc, assignment]
    SASA_AVAILABLE = False

# Distance limits (Å) for contact types - contacts outside these are ignored
H_BOND_MAX = 3.5
ELECTROSTATIC_MAX = 5.0
HYDROPHOBIC_MIN, HYDROPHOBIC_MAX = 3.5, 4.5
VDW_MAX = 5.0


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


class InterfaceAnalyzer:
    """Analyzer for protein-protein interfaces in crystal structures."""

    def __init__(self, contact_distance=5.0):
        """
        Initialize the interface analyzer.

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

        # Atomic radii for surface area calculations (Å)
        self.atomic_radii = {
            'C': 1.70, 'N': 1.55, 'O': 1.52, 'S': 1.80,
            'P': 1.80, 'H': 1.20, 'Ca': 2.31, 'Mg': 1.73,
            'Zn': 1.39, 'Fe': 1.32, 'Mn': 1.35, 'Cu': 1.32,
            'Na': 2.27, 'K': 2.75, 'Cl': 1.75, 'Br': 1.85
        }
    
    def analyze_interfaces(self, pdb_file):
        """
        Analyze all interfaces in a crystal structure.
        
        Parameters:
        -----------
        pdb_file : str
            Path to PDB file
            
        Returns:
        --------
        dict : Interface analysis results
        """
        if not BIOPYTHON_AVAILABLE or self.parser is None:
            return {'error': 'BioPython not available for interface analysis'}
        
        try:
            structure = self.parser.get_structure('crystal', pdb_file)
            if structure is None:
                return {'error': 'Failed to parse structure'}
            results = {
                'interfaces': [],
                'summary': {}
            }
            # Get all chains
            chains = list(structure.get_chains())
            
            if len(chains) < 2:
                return {'error': 'Need at least 2 chains for interface analysis'}
            
            # Analyze all pairwise interfaces
            total_buried_area = 0
            total_contact_area = 0
            total_interface_rmsd = 0.0
            interface_count = 0
            rmsd_interface_count = 0
            
            for i, chain1 in enumerate(chains):
                for j, chain2 in enumerate(chains):
                    if i >= j:  # Avoid duplicates and self-comparison
                        continue
                    # Use the model that contains both chains (required for multi-model e.g. NMR)
                    if chain1.parent is not chain2.parent:
                        continue
                    model = chain1.parent
                    interface_data = self._analyze_pairwise_interface(chain1, chain2, model)
                    
                    if interface_data['contact_count'] > 0 or interface_data.get('buried_surface_area', 0) > 0:
                        interface_data['chain_pair'] = f"{chain1.id}-{chain2.id}"
                        results['interfaces'].append(interface_data)
                        
                        total_buried_area += interface_data.get('buried_surface_area', 0)
                        total_contact_area += interface_data.get('contact_area', 0)
                        rmsd_val = interface_data.get('interface_rmsd_ca')
                        if rmsd_val is not None and isinstance(rmsd_val, (int, float)):
                            total_interface_rmsd += float(rmsd_val)
                            rmsd_interface_count += 1
                        interface_count += 1
            
            # Calculate summary statistics
            results['summary'] = {
                'total_interfaces': interface_count,
                'total_buried_surface_area': total_buried_area,
                'total_contact_area': total_contact_area,
                'average_buried_area_per_interface': total_buried_area / interface_count if interface_count > 0 else 0,
                'average_contact_area_per_interface': total_contact_area / interface_count if interface_count > 0 else 0,
                'average_interface_rmsd_ca': total_interface_rmsd / rmsd_interface_count if rmsd_interface_count > 0 else 0,
            }
            
            return results
            
        except Exception as e:
            return {'error': f"Failed to analyze interfaces: {str(e)}"}
    
    def _analyze_pairwise_interface(self, chain1, chain2, model):
        """
        Analyze interface between two chains using Biopython wiki approach:
        - Interface residues from NeighborSearch (residue pairs within threshold, different chains).
        - BSA from ShrakeRupley SASA when available: BSA = SASA(c1) + SASA(c2) - SASA(complex).
        - Polarity and accessibility summaries for contact residues.
        - Interface RMSD for matched contact residues (CA atoms) when possible.
        """
        results = {
            'chain1_id': chain1.id,
            'chain2_id': chain2.id,
        }

        # Get atoms from both chains
        atoms1 = list(chain1.get_atoms())
        atoms2 = list(chain2.get_atoms())
        atom_list = atoms1 + atoms2

        # Interface residues from NeighborSearch (Biopython wiki): residue pairs within threshold, different chains
        interface_residues1 = set()
        interface_residues2 = set()
        if NeighborSearch is not None and atom_list:
            ns = NeighborSearch(atom_list)
            for res_a, res_b in ns.search_all(self.contact_distance, 'R'):
                if res_a.parent.id == res_b.parent.id:
                    continue
                if res_a.parent.id == chain1.id and res_b.parent.id == chain2.id:
                    interface_residues1.add(res_a)
                    interface_residues2.add(res_b)
                elif res_a.parent.id == chain2.id and res_b.parent.id == chain1.id:
                    interface_residues1.add(res_b)
                    interface_residues2.add(res_a)

        # Find raw contacts, then filter by type-specific distance limits (same criteria for contact list)
        raw_contacts = self._find_contacts(atoms1, atoms2)
        contacts, type_counts = self._filter_contacts_by_limits(raw_contacts)
        results['contact_count'] = len(contacts)
        results['contact_type_counts'] = type_counts
        results['hydrogen_bonds'] = type_counts.get('hydrogen_bond', 0)

        # Contact residues: NeighborSearch interface residues filtered by the same criteria (distance + type limits)
        contact_residues1 = {r for r in interface_residues1 if self._residue_passes_contact_criteria(r, atoms2)}
        contact_residues2 = {r for r in interface_residues2 if self._residue_passes_contact_criteria(r, atoms1)}

        # If NeighborSearch was unavailable or returned nothing, derive contact residues from the filtered contacts
        if (not contact_residues1 and not contact_residues2) and contacts:
            for c in contacts:
                contact_residues1.add(c['residue1'])
                contact_residues2.add(c['residue2'])

        # Buried surface area (Biopython wiki): always compute for this chain pair
        buried_area = self._calculate_buried_surface_area(chain1, chain2, model)
        results['buried_surface_area'] = buried_area

        # Optional per-residue SASA for accessibility; computed once per chain pair
        residue_sasa = self._compute_residue_sasa_for_chains(model, (chain1.id, chain2.id))

        # Interface RMSD (CA atoms) for matched contact residues, when possible
        interface_rmsd = self._calculate_interface_rmsd(chain1, chain2, contact_residues1, contact_residues2)
        results['interface_rmsd_ca'] = interface_rmsd

        if len(contacts) == 0:
            results.update({
                'contact_area': 0,
                'interface_complementarity': 0,
                'charge_complementarity': 0,
                'charge_complementarity_opposite': 0,
                'charge_complementarity_same': 0,
                'charge_complementarity_total': 0,
                'charge_complementarity_density': None,
                'charge_complementarity_density_denominator': None,
                'contact_residues_chain1': [f"{r.resname}{r.id[1]}" for r in sorted(contact_residues1, key=lambda r: r.id[1])],
                'contact_residues_chain2': [f"{r.resname}{r.id[1]}" for r in sorted(contact_residues2, key=lambda r: r.id[1])],
                'num_contact_residues_chain1': len(contact_residues1),
                'num_contact_residues_chain2': len(contact_residues2),
                'average_contact_distance': 0,
                'min_contact_distance': 0,
                'max_contact_distance': 0,
                'contact_distance_std': 0,
            })
            # Even if there are no atom-level contacts after filtering, report polarity/accessibility
            results.update(self._summarize_polarity_and_accessibility(contact_residues1, contact_residues2, residue_sasa, chain1.id, chain2.id))
            return results

        # Analyze contact details (on filtered contacts)
        contact_analysis = self._analyze_contacts(contacts)
        results.update(contact_analysis)

        # Contact area (unchanged)
        contact_area = self._calculate_contact_area(contacts)
        results['contact_area'] = contact_area

        # Contact residues: NeighborSearch interface residues filtered by same criteria (distance + type) as contact list
        results['contact_residues_chain1'] = [f"{r.resname}{r.id[1]}" for r in sorted(contact_residues1, key=lambda r: r.id[1])]
        results['contact_residues_chain2'] = [f"{r.resname}{r.id[1]}" for r in sorted(contact_residues2, key=lambda r: r.id[1])]
        results['num_contact_residues_chain1'] = len(contact_residues1)
        results['num_contact_residues_chain2'] = len(contact_residues2)

        complementarity = self._calculate_interface_complementarity(contacts)
        results['interface_complementarity'] = complementarity

        charge_comp_result = self._calculate_charge_complementarity(contacts)
        results['charge_complementarity'] = charge_comp_result['score']
        results['charge_complementarity_opposite'] = charge_comp_result['opposite']
        results['charge_complementarity_same'] = charge_comp_result['same']
        results['charge_complementarity_total'] = charge_comp_result['total']

        # Charge complementarity density: opposite-sign charged contacts per unit interface area
        # (uses contact area when > 0, else BSA, so the metric is comparable across interfaces)
        area_for_density = contact_area if contact_area >= 0.01 else (results.get('buried_surface_area') or 0)
        if area_for_density > 0:
            results['charge_complementarity_density'] = charge_comp_result['opposite'] / float(area_for_density)
            results['charge_complementarity_density_denominator'] = 'contact_area' if contact_area >= 0.01 else 'buried_surface_area'
        else:
            results['charge_complementarity_density'] = None
            results['charge_complementarity_density_denominator'] = None

        # Polarity and accessibility summaries for contact residues
        results.update(self._summarize_polarity_and_accessibility(contact_residues1, contact_residues2, residue_sasa, chain1.id, chain2.id))

        return results
    
    def _residue_passes_contact_criteria(self, residue, other_chain_atoms):
        """True if residue has at least one atom that, with some atom in other_chain_atoms, satisfies the same distance + type-specific criteria as the contact list."""
        for atom in residue.get_atoms():
            for other in other_chain_atoms:
                d = float(np.linalg.norm(atom.coord - other.coord))
                if d <= self.contact_distance:
                    contact_type = self._classify_contact_type(atom, other)
                    if _contact_within_limit(contact_type, d):
                        return True
        return False
    
    def _find_contacts(self, atoms1, atoms2):
        """Find all contacts between two sets of atoms."""
        contacts = []
        
        # Simple O(n²) approach for smaller sets
        for atom1 in atoms1:
            for atom2 in atoms2:
                distance = np.linalg.norm(atom1.coord - atom2.coord)
                if distance <= self.contact_distance:
                    contacts.append({
                        'atom1': atom1,
                        'atom2': atom2,
                        'distance': distance,
                        'residue1': atom1.parent,
                        'residue2': atom2.parent
                    })
        
        return contacts

    def _classify_contact_type(self, atom1, atom2):
        """Classify contact type from atom elements (same scheme as contact_analyzer)."""
        element1 = atom1.element.strip().upper() or atom1.name[0].upper()
        element2 = atom2.element.strip().upper() or atom2.name[0].upper()
        h_bond_elements = {'N', 'O', 'S'}
        if element1 in h_bond_elements and element2 in h_bond_elements:
            return 'hydrogen_bond'
        if element1 == 'C' and element2 == 'C':
            return 'hydrophobic'
        charged_elements = {'N', 'O', 'S', 'P'}
        if element1 in charged_elements or element2 in charged_elements:
            return 'electrostatic'
        return 'van_der_waals'

    def _filter_contacts_by_limits(self, raw_contacts):
        """Keep only contacts within type-specific distance limits; return (filtered list, type_counts)."""
        type_counts = {'hydrogen_bond': 0, 'electrostatic': 0, 'hydrophobic': 0, 'van_der_waals': 0}
        filtered = []
        for c in raw_contacts:
            contact_type = self._classify_contact_type(c['atom1'], c['atom2'])
            d = float(c['distance'])
            if _contact_within_limit(contact_type, d):
                filtered.append(c)
                type_counts[contact_type] = type_counts.get(contact_type, 0) + 1
        return filtered, type_counts

    def _analyze_contacts(self, contacts):
        """Analyze contact statistics."""
        if not contacts:
            return {'average_contact_distance': 0}
        
        distances = [contact['distance'] for contact in contacts]
        
        return {
            'average_contact_distance': np.mean(distances),
            'min_contact_distance': np.min(distances),
            'max_contact_distance': np.max(distances),
            'contact_distance_std': np.std(distances)
        }
    
    def _calculate_buried_surface_area(self, chain1, chain2, model):
        """
        Buried surface area using the Biopython Interface Analysis approach:
        BSA = SASA(chain1 alone) + SASA(chain2 alone) - SASA(complex).
        Uses Bio.PDB.SASA.ShrakeRupley when available; otherwise falls back
        to a per-contact-atom estimate (~15 Å² per atom).
        """
        if not SASA_AVAILABLE or ShrakeRupley is None or Structure is None or Model is None:
            return self._estimate_buried_surface_area_fallback(chain1, chain2)

        try:
            sr = ShrakeRupley(probe_radius=1.4, n_points=100)
            # SASA of complex (both chains in one model)
            sr.compute(model, level='S')
            sasa_complex = getattr(model, 'sasa', None)
            if sasa_complex is None:
                return self._estimate_buried_surface_area_fallback(chain1, chain2)

            # Single-chain structures for SASA in isolation (use copies to avoid mutating model).
            # Use a fresh ShrakeRupley per chain so internal state from one compute() does not affect another.
            def _single_chain_sasa(chain):
                try:
                    c_copy = copy.deepcopy(chain)
                except (RecursionError, TypeError):
                    return None
                st = Structure("single")
                m = Model(0)
                m.add(c_copy)
                st.add(m)
                sr_chain = ShrakeRupley(probe_radius=1.4, n_points=100)
                sr_chain.compute(st, level='S')
                return getattr(st, 'sasa', None)

            sasa1 = _single_chain_sasa(chain1)
            sasa2 = _single_chain_sasa(chain2)
            if sasa1 is None or sasa2 is None:
                return self._estimate_buried_surface_area_fallback(chain1, chain2)

            bsa = float(sasa1) + float(sasa2) - float(sasa_complex)
            return max(0.0, bsa)
        except Exception:
            return self._estimate_buried_surface_area_fallback(chain1, chain2)

    def _estimate_buried_surface_area_fallback(self, chain1, chain2):
        """Fallback BSA estimate when SASA is unavailable (~15 Å² per contact atom)."""
        atoms1 = set(chain1.get_atoms())
        atoms2 = set(chain2.get_atoms())
        # Approximate contact atoms as those within contact_distance of the other chain
        contact_atoms = set()
        for a1 in atoms1:
            for a2 in atoms2:
                if np.linalg.norm(a1.coord - a2.coord) <= self.contact_distance:
                    contact_atoms.add(a1)
                    contact_atoms.add(a2)
        return len(contact_atoms) * 15.0 if contact_atoms else 0.0
    
    def _calculate_contact_area(self, contacts):
        """Calculate total contact area."""
        # Simplified: estimate contact area based on contact distance
        total_area = 0
        
        for contact in contacts:
            # Estimate contact patch area based on atomic radii and distance
            atom1 = contact['atom1']
            atom2 = contact['atom2']
            
            element1 = atom1.element.strip().upper() or atom1.name[0].upper()
            element2 = atom2.element.strip().upper() or atom2.name[0].upper()
            
            radius1 = self.atomic_radii.get(element1, 1.7)
            radius2 = self.atomic_radii.get(element2, 1.7)
            
            # Simple contact area estimation
            distance = contact['distance']
            max_radius = max(radius1, radius2)
            
            if distance < (radius1 + radius2):
                # Overlapping spheres - calculate contact area
                area = np.pi * min(radius1, radius2)**2 * (1 - distance/(radius1 + radius2))
                total_area += max(0.0, float(area))
        
        return total_area
    
    def _calculate_interface_complementarity(self, contacts):
        """Calculate interface shape complementarity (simplified)."""
        if not contacts:
            return 0
        
        # Simplified complementarity based on contact distance distribution
        distances = [contact['distance'] for contact in contacts]
        
        # Good complementarity = tight distribution around optimal distance
        optimal_distance = 3.5  # Typical van der Waals contact distance
        deviations = [abs(d - optimal_distance) for d in distances]
        
        # Convert to complementarity score (0-1, higher is better)
        avg_deviation = np.mean(deviations)
        complementarity = max(0.0, float(1.0 - avg_deviation / 2.0))  # Normalize
        
        return complementarity

    def _residue_charge(self, resname: str) -> int:
        """
        Return residue charge sign for complementarity calculation.
        +1 for positively charged (ARG, LYS, HIS),
        -1 for negatively charged (ASP, GLU),
        0 for all others.
        """
        name = resname.upper()
        if name in ('ASP', 'GLU'):
            return -1
        if name in ('ARG', 'LYS', 'HIS'):
            return 1
        return 0

    def _calculate_charge_complementarity(self, contacts):
        """
        Calculate charge complementarity score (0-1).

        Defined as the fraction of charged–charged contacts that are
        oppositely charged:

            score = opposite / (opposite + same)

        Counts are over all atom–atom contacts: each contact where both
        residues are charged (ASP/GLU/ARG/LYS/HIS) is classified as
        opposite-sign (attractive) or same-sign (repulsive). Contacts
        with a neutral residue are ignored. Returns dict with 'score',
        'opposite', 'same', 'total' for diagnostics.
        """
        if not contacts:
            return {'score': 0.0, 'opposite': 0, 'same': 0, 'total': 0}

        opposite = 0
        same = 0

        for c in contacts:
            r1 = c.get('residue1')
            r2 = c.get('residue2')
            if r1 is None or r2 is None:
                continue
            q1 = self._residue_charge(getattr(r1, 'resname', ''))
            q2 = self._residue_charge(getattr(r2, 'resname', ''))
            if q1 == 0 or q2 == 0:
                continue
            prod = q1 * q2
            if prod < 0:
                opposite += 1
            elif prod > 0:
                same += 1

        total = opposite + same
        score = opposite / float(total) if total else 0.0
        return {'score': score, 'opposite': opposite, 'same': same, 'total': total}

    # --- Polarity and accessibility helpers ---------------------------------

    def _calculate_interface_rmsd(self, chain1, chain2, contact_residues1, contact_residues2):
        """
        Calculate RMSD between matching contact residues of two chains using CA atoms.

        This is mainly meaningful for homodimeric interfaces where residue IDs / names correspond.
        Returns RMSD (Å) on success (including 0.0 for perfect superposition), or None if
        not computable (unavailable Superimposer, no/few matched residues, etc.).
        """
        # Superimposer is only available when BIOPYTHON_AVAILABLE was True
        if not BIOPYTHON_AVAILABLE or Superimposer is None:
            return None

        if not contact_residues1 or not contact_residues2:
            return None

        # Build quick lookup for residues in chain2 by (resname, residue number, insertion code)
        res_index_chain2 = {}
        for r in contact_residues2:
            key = (r.resname, r.id[1], r.id[2] if len(r.id) > 2 else '')
            res_index_chain2[key] = r

        fixed_atoms = []
        moving_atoms = []

        for r1 in contact_residues1:
            key = (r1.resname, r1.id[1], r1.id[2] if len(r1.id) > 2 else '')
            r2 = res_index_chain2.get(key)
            if r2 is None:
                continue
            # Use CA atoms for RMSD; skip if missing
            try:
                ca1 = r1['CA']
                ca2 = r2['CA']
            except KeyError:
                continue
            fixed_atoms.append(ca1)
            moving_atoms.append(ca2)

        if len(fixed_atoms) < 3:
            # Need at least 3 points for a stable superposition
            return None

        try:
            sup = Superimposer()
            sup.set_atoms(fixed_atoms, moving_atoms)
            rms = getattr(sup, 'rms', None)
            if rms is None:
                return None
            return float(rms)
        except Exception:
            return None

    def _classify_residue_polarity(self, residue):
        """
        Classify residue into simple polarity categories:
        - 'charged': ARG, LYS, ASP, GLU, HIS
        - 'polar': SER, THR, ASN, GLN, TYR, CYS
        - 'apolar': ALA, VAL, LEU, ILE, MET, PHE, TRP, PRO, GLY
        - 'other': anything else (e.g. ligands)
        """
        name = (residue.resname or "").strip().upper()
        charged = {'ARG', 'LYS', 'ASP', 'GLU', 'HIS'}
        polar = {'SER', 'THR', 'ASN', 'GLN', 'TYR', 'CYS'}
        apolar = {'ALA', 'VAL', 'LEU', 'ILE', 'MET', 'PHE', 'TRP', 'PRO', 'GLY'}
        if name in charged:
            return 'charged'
        if name in polar:
            return 'polar'
        if name in apolar:
            return 'apolar'
        return 'other'

    def _compute_residue_sasa_for_chains(self, model, chain_ids):
        """
        Compute per-residue SASA for a subset of chains in the given model.
        Returns a dict keyed by (chain_id, residue.id) → SASA (Å²).

        Uses ShrakeRupley(level='R') when available; otherwise returns an empty dict.
        """
        if not SASA_AVAILABLE or ShrakeRupley is None:
            return {}
        try:
            # Work on a deepcopy of the model to avoid mutating the original (which may be
            # reused for structure-level SASA in _calculate_buried_surface_area for other pairs).
            try:
                model_copy = copy.deepcopy(model)
            except (RecursionError, TypeError):
                # Do not run sr.compute on the original model; it would overwrite structure-level
                # SASA attributes and corrupt BSA calculations for other chain pairs.
                return {}

            sr = ShrakeRupley(probe_radius=1.4, n_points=100)
            sr.compute(model_copy, level='R')

            residue_sasa = {}
            for chain in model_copy.get_chains():
                if chain.id not in chain_ids:
                    continue
                for residue in chain.get_residues():
                    sasa = getattr(residue, 'sasa', None)
                    if sasa is not None:
                        residue_sasa[(chain.id, residue.id)] = float(sasa)
            return residue_sasa
        except Exception:
            return {}

    def _summarize_polarity_and_accessibility(self, contact_residues1, contact_residues2, residue_sasa, chain1_id, chain2_id):
        """
        Summarize polarity and accessibility for contact residues on both chains.

        Accessibility is reported as:
        - average_sasa: mean SASA over contact residues with SASA information
        - accessible_fraction: fraction of contact residues with SASA > 0
        """
        def summarize_chain(residues, chain_id):
            polarity_counts = {'charged': 0, 'polar': 0, 'apolar': 0, 'other': 0}
            for r in residues:
                cat = self._classify_residue_polarity(r)
                polarity_counts[cat] = polarity_counts.get(cat, 0) + 1

            total_res = len(residues)
            polarity_fractions = {}
            if total_res > 0:
                for k, v in polarity_counts.items():
                    polarity_fractions[k] = v / total_res

            # Accessibility from per-residue SASA
            sasa_values = []
            accessible = 0
            for r in residues:
                sasa = residue_sasa.get((chain_id, r.id))
                if sasa is None:
                    continue
                sasa_values.append(sasa)
                if sasa > 0.0:
                    accessible += 1

            if sasa_values:
                avg_sasa = float(np.mean(sasa_values))
                accessible_fraction = accessible / len(sasa_values)
            else:
                avg_sasa = 0.0
                accessible_fraction = 0.0

            return {
                'polarity_counts': polarity_counts,
                'polarity_fractions': polarity_fractions,
                'average_sasa': avg_sasa,
                'accessible_fraction': accessible_fraction,
                'residues_with_sasa': len(sasa_values),
            }

        chain1_summary = summarize_chain(contact_residues1, chain1_id)
        chain2_summary = summarize_chain(contact_residues2, chain2_id)

        return {
            'polarity_chain1': chain1_summary['polarity_counts'],
            'polarity_chain2': chain2_summary['polarity_counts'],
            'polarity_fractions_chain1': chain1_summary['polarity_fractions'],
            'polarity_fractions_chain2': chain2_summary['polarity_fractions'],
            'accessibility_chain1': {
                'average_sasa': chain1_summary['average_sasa'],
                'accessible_fraction': chain1_summary['accessible_fraction'],
                'residues_with_sasa': chain1_summary['residues_with_sasa'],
            },
            'accessibility_chain2': {
                'average_sasa': chain2_summary['average_sasa'],
                'accessible_fraction': chain2_summary['accessible_fraction'],
                'residues_with_sasa': chain2_summary['residues_with_sasa'],
            },
        }

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
    """
    Return paths whose basename contains every pattern in `patterns`.
    Patterns are matched as substrings (e.g. 'tag1' and 'tag2' match 'prefix_tag1_tag2.pdb').
    """
    if not patterns:
        return list(paths)
    result = []
    for p in paths:
        name = os.path.basename(p)
        if all(pat in name for pat in patterns):
            result.append(p)
    return result


def _run_analysis(analyzer, paths, out_stream):
    """Run interface analysis on paths and write results to out_stream."""
    for i, pdb_file in enumerate(paths):
        if len(paths) > 1:
            print(f"\n{'='*50}\n[{i+1}/{len(paths)}] {os.path.basename(pdb_file)}\n{'='*50}", file=out_stream)
        else:
            print(f"Analyzing interfaces in {pdb_file}...", file=out_stream)

        results = analyzer.analyze_interfaces(pdb_file)

        if 'error' in results:
            print(f"Error: {results['error']}", file=out_stream)
            continue

        print("\nINTERFACE ANALYSIS RESULTS:", file=out_stream)
        print("=" * 40, file=out_stream)
        summary = results.get('summary')
        if isinstance(summary, dict):
            print(f"Total interfaces: {summary.get('total_interfaces', 0)}", file=out_stream)
            print(f"Total buried surface area: {summary.get('total_buried_surface_area', 0):.1f} Å²", file=out_stream)
            print(f"Average buried area per interface: {summary.get('average_buried_area_per_interface', 0):.1f} Å²", file=out_stream)
        interfaces = results.get('interfaces', [])
        if isinstance(interfaces, list):
            for j, interface in enumerate(interfaces):
                if isinstance(interface, dict):
                    print(f"\nInterface {j+1}: {interface.get('chain_pair', '')}", file=out_stream)
                    print(f"  Contact count (within limits): {interface.get('contact_count', 0)}", file=out_stream)
                    tc = interface.get('contact_type_counts') or {}
                    if isinstance(tc, dict):
                        print(f"  H-bonds (≤3.5 Å): {tc.get('hydrogen_bond', 0)}  electrostatic (≤5 Å): {tc.get('electrostatic', 0)}  hydrophobic (3.5–4.5 Å): {tc.get('hydrophobic', 0)}  van der Waals (≤5 Å): {tc.get('van_der_waals', 0)}", file=out_stream)
                    print(f"  Buried surface area: {interface.get('buried_surface_area', 0):.1f} Å²  Contact area: {interface.get('contact_area', 0):.2f} Å²", file=out_stream)
                    print(f"  Complementarity (shape): {interface.get('interface_complementarity', 0):.3f}", file=out_stream)
                    cc_total = interface.get('charge_complementarity_total', 0)
                    if cc_total > 0:
                        cc_opp = interface.get('charge_complementarity_opposite', 0)
                        cc_same = interface.get('charge_complementarity_same', 0)
                        print(f"  Complementarity (charge): {interface.get('charge_complementarity', 0):.3f}  (charged–charged contacts: {cc_opp} opposite, {cc_same} same)", file=out_stream)
                    else:
                        print(f"  Complementarity (charge): {interface.get('charge_complementarity', 0):.3f}  (no charged–charged contacts)", file=out_stream)
                    ccd = interface.get('charge_complementarity_density')
                    if ccd is not None:
                        denom = interface.get('charge_complementarity_density_denominator') or 'area'
                        print(f"  Complementarity (charge) density: {ccd:.4f} opposite-sign contacts/Å²  (per {denom})", file=out_stream)
                    rmsd_val = interface.get('interface_rmsd_ca')
                    if rmsd_val is not None:
                        print(f"  Interface RMSD (CA): {rmsd_val:.3f} Å", file=out_stream)
                    else:
                        print("  Interface RMSD (CA): N/A", file=out_stream)
                    print(f"  Distance (Å): min={interface.get('min_contact_distance', 0):.2f} max={interface.get('max_contact_distance', 0):.2f} avg={interface.get('average_contact_distance', 0):.2f}", file=out_stream)
                    c1 = interface.get('chain1_id', '')
                    c2 = interface.get('chain2_id', '')
                    pol1 = interface.get('polarity_chain1') or {}
                    pol2 = interface.get('polarity_chain2') or {}
                    print(f"  Polarity {c1}: charged={pol1.get('charged', 0)} polar={pol1.get('polar', 0)} apolar={pol1.get('apolar', 0)} other={pol1.get('other', 0)}", file=out_stream)
                    print(f"  Polarity {c2}: charged={pol2.get('charged', 0)} polar={pol2.get('polar', 0)} apolar={pol2.get('apolar', 0)} other={pol2.get('other', 0)}", file=out_stream)
                    acc1 = interface.get('accessibility_chain1') or {}
                    acc2 = interface.get('accessibility_chain2') or {}
                    print(f"  Accessibility {c1}: avg SASA={acc1.get('average_sasa', 0.0):.1f} Å²  accessible_fraction={acc1.get('accessible_fraction', 0.0):.2f}", file=out_stream)
                    print(f"  Accessibility {c2}: avg SASA={acc2.get('average_sasa', 0.0):.1f} Å²  accessible_fraction={acc2.get('accessible_fraction', 0.0):.2f}", file=out_stream)

    if len(paths) > 1:
        print(f"\nDone. Processed {len(paths)} file(s).", file=out_stream)


def main():
    """CLI: single file, multiple files, directory, or glob patterns; optional multi-set filtering and dry-run."""
    import argparse

    parser = argparse.ArgumentParser(
        description="Analyze interfaces between molecules in crystal structures.",
        epilog="""Examples:
  python interface_analyzer.py structure.pdb
  python interface_analyzer.py dir/
  python interface_analyzer.py *.pdb -o results.txt
  python interface_analyzer.py *.pdb --set tag1,tag2 -o output_set.txt
  python interface_analyzer.py *.pdb --sets set1 set2 set3 -o "output_{}.txt"
  python interface_analyzer.py *.pdb --per-structure -o "{}_interface.txt"
  python interface_analyzer.py *.pdb --per-structure --sets tag_a tag_b -o "{}_interface.txt"
  python interface_analyzer.py *.pdb --sets set1 set2 set3 set4 -o "output_{}.txt" --dry-run
""",
    )
    parser.add_argument(
        'input',
        nargs='+',
        help='PDB/CIF file(s), directory, or glob pattern (e.g. *.pdb)',
    )
    parser.add_argument(
        '--set', '-s',
        action='append',
        dest='sets',
        metavar='PATTERNS',
        help='Comma-separated patterns; file is included only if basename contains ALL patterns. Repeat for multiple sets (e.g. -s tag1,tag2 -s set1).',
    )
    parser.add_argument(
        '--sets',
        dest='sets_multi',
        nargs='+',
        metavar='SET',
        help='Multiple set names in one go (one pattern per set). E.g. --sets set1 set2 set3 is equivalent to --set set1 --set set2 --set set3.',
    )
    parser.add_argument(
        '--output', '-o',
        metavar='FILE',
        help='Output file. Single file: -o results.txt (all sets concatenated). Per-set: use "{}" for set label. Per-structure: use --per-structure with "{}" for file stem.',
    )
    parser.add_argument(
        '--per-structure', '-p',
        action='store_true',
        help='Write one output file per structure; -o must contain "{}" (replaced by file stem). Can combine with --set/--sets to filter which files are processed.',
    )
    parser.add_argument(
        '--dry-run',
        action='store_true',
        help='Print which files would be processed for each set and exit without running. Writes summary to interface_analyzer_dryrun.txt in the current directory.',
    )
    args = parser.parse_args()

    paths = collect_structure_paths(args.input)

    # Per-structure mode: one output file per input file (optionally filtered by --set/--sets)
    if args.per_structure:
        if not args.output or '{}' not in args.output:
            print("Error: --per-structure requires -o with '{}' in the path (e.g. -o '{}_interface.txt').", file=sys.stderr)
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
        structure_list = [(os.path.splitext(os.path.basename(p))[0], [p]) for p in paths_to_process]
        if args.dry_run:
            lines = ["Dry run (per-structure mode): no analysis will be performed.", ""]
            for stem, [path] in structure_list:
                dest = args.output.replace('{}', stem)
                lines.append(f"  {os.path.basename(path)} -> {dest}")
            lines.append("")
            out_text = "\n".join(lines)
            dryrun_file = os.path.join(os.getcwd(), "interface_analyzer_dryrun.txt")
            try:
                with open(dryrun_file, "w") as f:
                    f.write(out_text)
                print(out_text, file=sys.stderr, flush=True)
                print(f"Dry-run summary written to {dryrun_file}", file=sys.stderr, flush=True)
            except OSError as e:
                print(out_text, file=sys.stderr, flush=True)
                print(f"Could not write {dryrun_file}: {e}", file=sys.stderr, flush=True)
            return
        analyzer = InterfaceAnalyzer()
        for stem, [path] in structure_list:
            out_path = args.output.replace('{}', stem)
            print(f"Writing {os.path.basename(path)} -> {out_path}", file=sys.stderr)
            with open(out_path, 'w') as out:
                _run_analysis(analyzer, [path], out)
        return

    # Build list of (label, filtered_paths). If --set/--sets not used, one set with all paths.
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
        lines = ["Dry run: no analysis will be performed.", ""]
        if not paths:
            lines.extend([
                "No structure files found.",
                "The shell expands '*.pdb' in the current working directory.",
                f"Current directory: {os.getcwd()}",
                f"Input received: {args.input}",
                "Run from the directory that contains your PDB files, or pass explicit paths/directories.",
                "",
            ])
        single_file = args.output and '{}' not in args.output
        for label, patterns, filtered in set_list:
            if single_file:
                dest = args.output
            else:
                dest = (args.output.replace('{}', label) if args.output and '{}' in args.output
                        else (args.output if len(set_list) == 1 and args.output else None))
            dest = dest or "stdout"
            lines.append(f"Set {label!r} (patterns: {patterns}): {len(filtered)} file(s) -> {dest}")
            for p in filtered:
                lines.append(f"  {os.path.basename(p)}")
            lines.append("")
        out_text = "\n".join(lines)
        # Always write to a file so output is visible even if stderr is not
        dryrun_file = os.path.join(os.getcwd(), "interface_analyzer_dryrun.txt")
        try:
            with open(dryrun_file, "w") as f:
                f.write(out_text)
            # Also print to stderr
            print(out_text, file=sys.stderr, flush=True)
            print(f"Dry-run summary written to {dryrun_file}", file=sys.stderr, flush=True)
        except OSError as e:
            print(out_text, file=sys.stderr, flush=True)
            print(f"Could not write {dryrun_file}: {e}", file=sys.stderr, flush=True)
        return

    if not paths:
        print("No structure files found.", file=sys.stderr)
        print("Run from the directory that contains your PDB files, or pass explicit paths/directories.", file=sys.stderr)
        sys.exit(1)

    # Single output file for all sets, or one file per set when -o contains '{}'
    single_output_file = args.output and '{}' not in args.output
    if single_output_file:
        out = open(args.output, 'w')
        print(f"Writing all sets to {args.output}", file=sys.stderr)

    analyzer = InterfaceAnalyzer()
    for idx, (label, patterns, filtered) in enumerate(set_list):
        if not filtered:
            print(f"No files match set {label!r} (patterns: {patterns}), skipping.", file=sys.stderr)
            continue

        if args.output and not single_output_file:
            out_path = args.output.replace('{}', label)
            out = open(out_path, 'w')
            print(f"Writing analysis for set {label!r} ({len(filtered)} file(s)) to {out_path}", file=sys.stderr)
        elif not single_output_file:
            out = sys.stdout
            if len(set_list) > 1:
                print(f"\n=== Set {label!r} ({len(filtered)} file(s)) ===\n", file=sys.stderr)

        try:
            if len(set_list) > 1 and (single_output_file or not args.output):
                print(f"\n{'='*50}\nSet {label!r} (patterns: {patterns})\n{'='*50}", file=out)
            _run_analysis(analyzer, filtered, out)
        finally:
            if args.output and not single_output_file:
                out.close()

    if single_output_file:
        out.close()


if __name__ == "__main__":
    try:
        main()
    except Exception as e:
        print(f"Error: {e}", file=sys.stderr, flush=True)
        raise 