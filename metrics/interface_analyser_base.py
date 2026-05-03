#!/usr/bin/env python3
"""
Interface Analyser
=================

Analyse interfaces between molecules in crystal structures using the approach
described in Biopython's Interface Analysis wiki:
https://biopython.org/wiki/Interface_Analysis

- **Contact residues** (contact_residues_chain1/2): First identified by NeighborSearch (residue pairs within
  distance threshold, different chains). Then filtered by applying the same criteria as the contact list:
  keep only residues that have at least one atom pair (with the other chain) satisfying both the distance
  cutoff and the type-specific distance limits (H-bond, hydrophobic, etc.).
- **Buried surface area:** Shrake-Rupley SASA when available:
  BSA = SASA(chain1 alone) + SASA(chain2 alone) − SASA(two-chain complex).
  The complex uses only those two chains (deep-copied), so pairwise BSA is
  well-defined in multi-chain assemblies (e.g. A–B vs A–C vs A–D).
  Falls back to a per-atom estimate if SASA is not available.
"""

import copy
import glob
import math
import os
import sys
import threading
import warnings
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from typing import Any

import numpy as np

try:
    from Bio.PDB.PDBParser import PDBParser
    from Bio.PDB.NeighborSearch import NeighborSearch
    from Bio.PDB.PDBExceptions import PDBConstructionWarning
    from Bio.PDB.Superimposer import Superimposer
    # Bio.PDB.Structure / Bio.PDB.Model are modules; import the actual classes.
    from Bio.PDB.Structure import Structure as StructureClass
    from Bio.PDB.Model import Model as ModelClass
    warnings.simplefilter('ignore', PDBConstructionWarning)
    BIOPYTHON_AVAILABLE = True
except ImportError:
    PDBParser = None  # type: ignore[misc, assignment]
    NeighborSearch = None  # type: ignore[misc, assignment]
    PDBConstructionWarning = None  # type: ignore[misc, assignment]
    Superimposer = None  # type: ignore[misc, assignment]
    StructureClass = None  # type: ignore[misc, assignment]
    ModelClass = None  # type: ignore[misc, assignment]
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


def _foldkit_interface_workers() -> int:
    """
    Parallel worker count for pairwise interface analysis and EC loops.

    Environment variable ``FOLDKIT_INTERFACE_WORKERS``:
    - unset or invalid: 1 (sequential, backward compatible)
    - integer >= 2: cap at 64
    """
    raw = os.environ.get("FOLDKIT_INTERFACE_WORKERS", "").strip()
    if not raw:
        return 1
    try:
        n = int(raw)
    except ValueError:
        return 1
    if n <= 1:
        return 1
    return min(n, 64)


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


# Match `metrics/lattice_packing_analyser.py` for mass normalisation (Da).
_ATOMIC_WEIGHTS_DA: dict[str, float] = {
    "C": 12.011,
    "N": 14.007,
    "O": 15.999,
    "S": 32.06,
    "P": 30.974,
    "H": 1.008,
    "CA": 40.078,
    "MG": 24.305,
    "ZN": 65.38,
    "FE": 55.845,
    "MN": 54.938,
    "CU": 63.546,
    "NA": 22.99,
    "K": 39.098,
    "CL": 35.45,
    "BR": 79.904,
}


def _atom_element_upper(atom: Any) -> str:
    el = (getattr(atom, "element", "") or "").strip().upper()
    if el:
        return el
    name = (getattr(atom, "name", "") or "").strip()
    return (name[:1] or "C").upper()


def _atom_mass_da(atom: Any) -> float:
    return float(_ATOMIC_WEIGHTS_DA.get(_atom_element_upper(atom), _ATOMIC_WEIGHTS_DA["C"]))


def _polymer_chain_res_atom_mass(chain: Any) -> tuple[int, int, float]:
    """
    Polymer residues in a chain: hetero flag ' ' (standard PDB) and not water.
    Returns (n_residues, n_atoms, mass_da).
    """
    skip_res = {"HOH", "WAT", "SOL"}
    n_res = 0
    n_atoms = 0
    mass_da = 0.0
    for res in chain:
        try:
            rn = (getattr(res, "resname", "") or "").strip().upper()
        except Exception:
            continue
        if rn in skip_res:
            continue
        try:
            het = res.id[0]
        except Exception:
            het = " "
        if het != " ":
            continue
        n_res += 1
        for atom in res:
            n_atoms += 1
            mass_da += _atom_mass_da(atom)
    return n_res, n_atoms, mass_da


def _model_polymer_res_atom_mass(model: Any) -> tuple[int, int, float]:
    tr, ta, tm = 0, 0, 0.0
    for ch in model:
        r, a, m = _polymer_chain_res_atom_mass(ch)
        tr += r
        ta += a
        tm += m
    return tr, ta, tm


class InterfaceAnalyser:
    """Analyser for protein-protein interfaces in crystal structures."""

    def __init__(self, contact_distance=5.0, *, skip_accessibility_sasa: bool = False):
        """
        Initialise the interface analyser.

        Parameters:
        -----------
        contact_distance : float
            Maximum distance for considering atoms in contact (Å)
        skip_accessibility_sasa : bool
            When True, skip per-residue SASA in the full model used only for contact-residue
            accessibility summaries (faster on large assemblies; those fields reflect missing SASA).
        """
        if BIOPYTHON_AVAILABLE and PDBParser is not None:
            self.parser = PDBParser(QUIET=True)
        else:
            self.parser = None

        self.contact_distance = contact_distance
        self.skip_accessibility_sasa = bool(skip_accessibility_sasa)
        self._sasa_cache_lock = threading.Lock()

        # Atomic radii for surface area calculations (Å)
        self.atomic_radii = {
            'C': 1.70, 'N': 1.55, 'O': 1.52, 'S': 1.80,
            'P': 1.80, 'H': 1.20, 'Ca': 2.31, 'Mg': 1.73,
            'Zn': 1.39, 'Fe': 1.32, 'Mn': 1.35, 'Cu': 1.32,
            'Na': 2.27, 'K': 2.75, 'Cl': 1.75, 'Br': 1.85
        }
    
    def analyse_interfaces(self, pdb_file, focus_chains=None, reference_chain_id=None):
        """
        Analyse all interfaces in a crystal structure.
        
        Parameters:
        -----------
        pdb_file : str
            Path to PDB file
        focus_chains : None | list[str] | set[str]
            If provided, only analyse chain pairs where at least one chain ID
            is in this set (e.g. focus_chains={'A'} analyses all interfaces that
            involve chain A).
        reference_chain_id : None | str
            If set (e.g. central chain in a symmetry-expanded multi-copy PDB),
            compute lattice-style metrics: SASA of that chain alone vs sum of its
            per-residue SASA in the full model (all chains occluding), and
            fraction of reference residues with a cross-chain neighbour within
            contact_distance. Normally pairwise totals are computed twice so
            ``summary_all_chains`` / ``all_*`` aggregate **every** chain pair (subject to
            ``focus_chains``) while top-level ``summary`` stays reference-chain-filtered.
            When ``focus_chains`` is exactly ``{reference_chain_id}``, both summaries use the
            same pair list and only one pairwise pass is run.
            
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

            focus_set = None
            if focus_chains:
                focus_set = {str(c).strip() for c in focus_chains if str(c).strip()}

            ref_id = str(reference_chain_id).strip() if reference_chain_id else ""
            
            if len(chains) < 2:
                return {'error': 'Need at least 2 chains for interface analysis'}

            # Cleared each structure load so cached SASA maps never leak stale geometry between runs.
            self._model_residue_sasa_cache = {}
            self._chain_isolated_sasa_cache: dict[int, float | None] = {}

            def _analyse_pairs(*, reference_filter: str | None) -> tuple[list[dict[str, Any]], dict[str, Any]]:
                """
                Analyse interfaces for either:
                - reference_filter=None: all chain–chain pairs (subject to focus_set)
                - reference_filter='A': only interfaces involving that chain (subject to focus_set)
                """
                out_intfs: list[dict[str, Any]] = []
                total_buried_area = 0.0
                total_contact_area = 0.0
                total_interface_rmsd = 0.0
                interface_count = 0
                rmsd_interface_count = 0

                pair_specs: list[tuple[Any, Any, Any]] = []
                for i, chain1 in enumerate(chains):
                    for j, chain2 in enumerate(chains):
                        if i >= j:  # Avoid duplicates and self-comparison
                            continue
                        if focus_set is not None and (chain1.id not in focus_set and chain2.id not in focus_set):
                            continue
                        if reference_filter and chain1.id != reference_filter and chain2.id != reference_filter:
                            continue
                        # Use the model that contains both chains (required for multi-model e.g. NMR)
                        if chain1.parent is not chain2.parent:
                            continue
                        model = chain1.parent
                        pair_specs.append((chain1, chain2, model))

                workers_iface = _foldkit_interface_workers()
                if workers_iface <= 1 or len(pair_specs) <= 1:
                    interface_rows = [
                        self._analyse_pairwise_interface(c1, c2, m) for (c1, c2, m) in pair_specs
                    ]
                else:
                    def _iface_job(spec: tuple[Any, Any, Any]) -> dict[str, Any]:
                        c1, c2, m = spec
                        return self._analyse_pairwise_interface(c1, c2, m)

                    with ThreadPoolExecutor(max_workers=workers_iface) as pool:
                        interface_rows = list(pool.map(_iface_job, pair_specs))

                for interface_data, (chain1, chain2, _m) in zip(interface_rows, pair_specs):
                    if interface_data['contact_count'] > 0 or interface_data.get('buried_surface_area', 0) > 0:
                        interface_data['chain_pair'] = f"{chain1.id}-{chain2.id}"
                        out_intfs.append(interface_data)

                        total_buried_area += float(interface_data.get('buried_surface_area', 0) or 0.0)
                        total_contact_area += float(interface_data.get('contact_area', 0) or 0.0)
                        rmsd_val = interface_data.get('interface_rmsd_ca')
                        if rmsd_val is not None and isinstance(rmsd_val, (int, float)):
                            total_interface_rmsd += float(rmsd_val)
                            rmsd_interface_count += 1
                        interface_count += 1

                summary0: dict[str, Any] = {
                    'total_interfaces': int(interface_count),
                    'total_buried_surface_area': float(total_buried_area),
                    'total_contact_area': float(total_contact_area),
                    'average_buried_area_per_interface': float(total_buried_area) / float(interface_count) if interface_count > 0 else 0.0,
                    'average_contact_area_per_interface': float(total_contact_area) / float(interface_count) if interface_count > 0 else 0.0,
                    'average_interface_rmsd_ca': float(total_interface_rmsd) / float(rmsd_interface_count) if rmsd_interface_count > 0 else 0.0,
                }
                return out_intfs, summary0

            # All-chains summary (subject to focus_set) plus reference-chain subset when given.
            # When focus is exactly {reference chain}, enumerating pairs with vs without reference_filter
            # yields the same set — skip duplicate expensive pairwise passes.
            single_focus_matches_ref = bool(
                ref_id and focus_set is not None and len(focus_set) == 1 and ref_id in focus_set
            )
            if single_focus_matches_ref:
                interfaces_ref, summary_ref = _analyse_pairs(reference_filter=ref_id)
                interfaces_all, summary_all = interfaces_ref, summary_ref
            else:
                interfaces_all, summary_all = _analyse_pairs(reference_filter=None)
                interfaces_ref, summary_ref = _analyse_pairs(reference_filter=ref_id or None)

            # Backwards compatibility: keep `interfaces` as the reference-filtered list when
            # reference_chain_id is set; otherwise include all interfaces.
            results['interfaces'] = interfaces_ref if ref_id else interfaces_all

            # Backwards compatibility: keep top-level summary fields as the reference-filtered
            # values (current behaviour), and add all-chains fields alongside.
            results['summary'] = dict(summary_ref)
            results['summary']['summary_reference_chain'] = dict(summary_ref)
            results['summary']['summary_all_chains'] = dict(summary_all)
            for k, v in summary_all.items():
                results['summary'][f"all_{k}"] = v

            # Optional unit-suffixed aliases for downstream tabular consumers.
            # JSON keys in this module are historically unitless (e.g. *_buried_surface_area),
            # while some CSV/flattening tools use explicit unit suffixes (e.g. *_A2).
            results['summary']['total_buried_surface_area_A2'] = results['summary'].get('total_buried_surface_area')
            results['summary']['average_buried_area_per_interface_A2'] = results['summary'].get('average_buried_area_per_interface')
            results['summary']['total_contact_area_A2'] = results['summary'].get('total_contact_area')
            results['summary']['average_contact_area_per_interface_A2'] = results['summary'].get('average_contact_area_per_interface')
            results['summary']['all_total_buried_surface_area_A2'] = results['summary'].get('all_total_buried_surface_area')
            results['summary']['all_average_buried_area_per_interface_A2'] = results['summary'].get('all_average_buried_area_per_interface')

            # Total solvent-accessible surface (SASA) of each chain in isolation, for all chains
            # that appear in any reported interface (e.g. focus A with A–B and A–C → A, B, C).
            involved_ids_ref = set()
            for intf in interfaces_ref:
                involved_ids_ref.add(intf.get('chain1_id'))
                involved_ids_ref.add(intf.get('chain2_id'))
            involved_ids_ref.discard(None)

            involved_ids_all = set()
            for intf in interfaces_all:
                involved_ids_all.add(intf.get('chain1_id'))
                involved_ids_all.add(intf.get('chain2_id'))
            involved_ids_all.discard(None)

            chain_by_id = {}
            for ch in chains:
                if (ch.id in involved_ids_ref or ch.id in involved_ids_all) and ch.id not in chain_by_id:
                    chain_by_id[ch.id] = ch
            for cid in sorted(involved_ids_ref | involved_ids_all):
                if cid not in chain_by_id:
                    for model in structure:
                        for ch in model:
                            if ch.id == cid:
                                chain_by_id[cid] = ch
                                break
                        if cid in chain_by_id:
                            break

            def _sasa_dict(ids: set[str]) -> dict[str, float]:
                out0: dict[str, float] = {}
                for cid in sorted(ids):
                    ch = chain_by_id.get(cid)
                    if ch is None:
                        continue
                    sasa_val = self._compute_sasa_chain_isolated(ch)
                    if sasa_val is not None:
                        out0[cid] = float(sasa_val)
                return out0

            sasa_isolated_by_chain_ref = _sasa_dict({str(x) for x in involved_ids_ref if x})
            sasa_isolated_by_chain_all = _sasa_dict({str(x) for x in involved_ids_all if x})

            # Existing fields: keep them reference-filtered when ref_id is set; otherwise all.
            sasa_isolated_by_chain = sasa_isolated_by_chain_ref if ref_id else sasa_isolated_by_chain_all

            if sasa_isolated_by_chain:
                results['summary']['sasa_isolated_by_chain'] = sasa_isolated_by_chain
                results['summary']['sasa_isolated_sum'] = sum(sasa_isolated_by_chain.values())
            else:
                results['summary']['sasa_isolated_by_chain'] = {}
                results['summary']['sasa_isolated_sum'] = None

            # Unit-suffixed aliases for SASA sums (Å²) to match CSV column conventions.
            results['summary']['sasa_isolated_sum_A2'] = results['summary'].get('sasa_isolated_sum')

            # Additional explicit summaries (always present).
            results['summary']['all_sasa_isolated_by_chain'] = sasa_isolated_by_chain_all
            results['summary']['all_sasa_isolated_sum'] = (
                sum(sasa_isolated_by_chain_all.values()) if sasa_isolated_by_chain_all else None
            )
            results['summary']['reference_sasa_isolated_by_chain'] = sasa_isolated_by_chain_ref
            results['summary']['reference_sasa_isolated_sum'] = (
                sum(sasa_isolated_by_chain_ref.values()) if sasa_isolated_by_chain_ref else None
            )

            results['summary']['all_sasa_isolated_sum_A2'] = results['summary'].get('all_sasa_isolated_sum')
            results['summary']['reference_sasa_isolated_sum_A2'] = results['summary'].get('reference_sasa_isolated_sum')

            if ref_id:
                lattice = self._compute_lattice_reference_metrics(structure, ref_id)
                if lattice:
                    results['summary'].update(lattice)

            return results
            
        except Exception as e:
            return {'error': f"Failed to analyse interfaces: {str(e)}"}


def _fisher_z(r: float) -> float:
    r = max(min(float(r), 0.999999), -0.999999)
    return 0.5 * math.log((1.0 + r) / (1.0 - r))


def _inv_fisher_z(z: float) -> float:
    ez2 = math.exp(2.0 * float(z))
    return (ez2 - 1.0) / (ez2 + 1.0)


def _weighted_fisher_mean(pairs):
    num = 0.0
    den = 0.0
    for r, w in pairs:
        if r is None:
            continue
        w = float(w or 0.0)
        if w <= 0.0:
            continue
        num += w * _fisher_z(float(r))
        den += w
    if den <= 0.0:
        return None, 0.0
    return _inv_fisher_z(num / den), den


class InterfaceAnalyserEC(InterfaceAnalyser):
    """
    Electrostatic complementarity (EC) extension.

    This analyser disables charge-sign complementarity metrics and augments
    per-interface and lattice-wide results with EC fields.

    EC is computed as a facing-point Pearson correlation of electrostatic
    potentials across the interface, following:
    McCoy, Epa & Colman (1997) J Mol Biol 268:570–584. DOI: 10.1006/jmbi.1997.0987.
    """

    def __init__(
        self,
        contact_distance=5.0,
        ec_mode: str = "full",
        *,
        skip_accessibility_sasa: bool = False,
        ec_max_contact_points: int | None = None,
    ):
        super().__init__(
            contact_distance=contact_distance,
            skip_accessibility_sasa=skip_accessibility_sasa,
        )
        mode = (ec_mode or "full").strip().lower()
        self.ec_mode = mode
        self._ec_detail_fn = None
        self.ec_max_contact_points = (
            int(ec_max_contact_points) if ec_max_contact_points is not None else None
        )
        if self.ec_max_contact_points is not None and self.ec_max_contact_points < 4:
            self.ec_max_contact_points = 4

        if mode != "full":
            raise ValueError("ec_mode must be 'full'")

        try:
            from metrics.electrostatic_complementarity import compute_ec_complementarity_detailed

            self._ec_detail_fn = compute_ec_complementarity_detailed
        except Exception:
            try:
                from electrostatic_complementarity import compute_ec_complementarity_detailed  # type: ignore

                self._ec_detail_fn = compute_ec_complementarity_detailed
            except Exception:
                self._ec_detail_fn = None

    def analyse_interfaces(self, pdb_file, focus_chains=None, reference_chain_id=None):
        return self.analyze_interfaces(
            pdb_file, focus_chains=focus_chains, reference_chain_id=reference_chain_id
        )

    def _calculate_charge_complementarity(self, contacts):
        return {"score": 0.0, "opposite": 0, "same": 0, "total": 0}

    def _compute_lattice_charge_metrics(self, *args, **kwargs):
        return None

    def analyze_interfaces_sasa_only(self, pdb_file, focus_chains=None, reference_chain_id=None):
        """SASA/BSA + lattice summaries only (no EC); for two-phase runs with JSON sidecar."""
        return super().analyse_interfaces(
            pdb_file, focus_chains=focus_chains, reference_chain_id=reference_chain_id
        )

    def apply_ec_phase(self, pdb_file, results, reference_chain_id=None):
        """
        Given a results dict from the SASA/BSA phase, parse ``pdb_file`` once for chains
        and fill per-interface EC fields plus lattice EC summary (mutates ``results``).
        """
        if "error" in results or self._ec_detail_fn is None or self.parser is None:
            return results

        try:
            structure = self.parser.get_structure("_ec_iface", pdb_file)
        except Exception:
            return results

        chain_by_id = {}
        for ch in structure.get_chains():
            if ch.id not in chain_by_id:
                chain_by_id[ch.id] = ch

        interfaces = results.get("interfaces", []) or []
        ec_kw: dict[str, Any] = {"contact_distance": self.contact_distance}
        if getattr(self, "ec_max_contact_points", None) is not None:
            ec_kw["ec_max_contact_points"] = self.ec_max_contact_points

        ec_by_pair: dict[tuple[str, str], dict] = {}

        def _ec_job(idx_iface: tuple[int, dict[str, Any]]) -> tuple[int, dict[str, Any] | None]:
            _idx, iface0 = idx_iface
            c1 = chain_by_id.get(iface0.get("chain1_id"))
            c2 = chain_by_id.get(iface0.get("chain2_id"))
            if c1 is None or c2 is None:
                return _idx, None
            d0 = self._ec_detail_fn(c1, c2, **ec_kw)
            if not isinstance(d0, dict):
                return _idx, None
            return _idx, d0

        workers_ec = _foldkit_interface_workers()
        ec_tasks = list(enumerate(interfaces))
        if workers_ec <= 1 or len(ec_tasks) <= 1:
            ec_outputs = [_ec_job(t) for t in ec_tasks]
        else:
            with ThreadPoolExecutor(max_workers=workers_ec) as pool:
                ec_outputs = list(pool.map(_ec_job, ec_tasks))

        for idx, d in ec_outputs:
            iface = interfaces[idx]
            if not isinstance(d, dict):
                continue

            r_val = d.get("r")
            n_pairs = int(d.get("n_pairs", 0) or 0)
            iface["ec_r"] = r_val
            iface["ec_n_pairs"] = n_pairs

            id1, id2 = iface.get("chain1_id"), iface.get("chain2_id")
            if id1 is not None and id2 is not None:
                ec_by_pair[tuple(sorted((str(id1), str(id2))))] = d

            bsa = float(iface.get("buried_surface_area", 0) or 0.0)
            if r_val is not None and bsa > 0:
                iface["ec_density"] = float(r_val) / bsa
                iface["ec_density_denominator"] = "buried_surface_area"
            else:
                iface["ec_density"] = None
                iface["ec_density_denominator"] = None

        ref_id = str(reference_chain_id).strip() if reference_chain_id else ""
        if not ref_id:
            return results

        ref_chain = chain_by_id.get(ref_id)
        if ref_chain is None:
            return results

        partner_ids = set()
        for iface in interfaces:
            a = iface.get("chain1_id")
            b = iface.get("chain2_id")
            if a == ref_id and b:
                partner_ids.add(b)
            elif b == ref_id and a:
                partner_ids.add(a)

        per_partner = []
        by_partner = {}
        total_npairs = 0
        # Lattice EC density is normalised by the reference buried area:
        #   reference_buried_area = SASA_iso(reference) - SASA_cluster(reference)
        # The lattice SASA fields live under results['summary'].
        summary0 = results.get("summary")
        density_den = 0.0
        if isinstance(summary0, dict):
            sri = summary0.get("sasa_reference_isolated")
            src = summary0.get("sasa_reference_in_cluster")
            if isinstance(sri, (int, float)) and isinstance(src, (int, float)):
                density_den = float(sri) - float(src)
        if density_den < 0.0:
            density_den = 0.0

        for pid in sorted(partner_ids):
            partner = chain_by_id.get(pid)
            if partner is None:
                continue
            d = ec_by_pair.get(tuple(sorted((ref_id, str(pid)))))
            if not isinstance(d, dict):
                d = self._ec_detail_fn(ref_chain, partner, **ec_kw)
            if not isinstance(d, dict):
                continue
            r_val = d.get("r")
            n_pairs = int(d.get("n_pairs", 0) or 0)
            total_npairs += n_pairs

            bsa = 0.0
            for iface in interfaces:
                a = iface.get("chain1_id")
                b = iface.get("chain2_id")
                if (a == ref_id and b == pid) or (a == pid and b == ref_id):
                    bsa = float(iface.get("buried_surface_area", 0) or 0.0)
                    break

            by_partner[str(pid)] = {
                "partner_chain_id": pid,
                "ec_r": r_val,
                "ec_n_pairs": n_pairs,
                "buried_surface_area": bsa,
            }
            per_partner.append((r_val, bsa, n_pairs))

        bsa_pairs = [(r, bsa) for (r, bsa, _n) in per_partner]
        np_pairs = [(r, n) for (r, _bsa, n) in per_partner]
        r_bsa, bsa_sum = _weighted_fisher_mean(bsa_pairs)
        r_np, np_sum = _weighted_fisher_mean(np_pairs)

        # Store lattice EC fields under `summary` because the reporting code reads from there.
        summary = results.get("summary")
        if not isinstance(summary, dict):
            summary = {}
            results["summary"] = summary

        # Make denominator visible/inspectable in outputs.
        summary["reference_buried_area"] = float(density_den) if density_den > 0 else None

        summary["lattice_ec_r_bsa_weighted"] = r_bsa
        summary["lattice_ec_r_npairs_weighted"] = r_np
        summary["lattice_ec_total_npairs"] = int(total_npairs)
        summary["lattice_ec_density_denominator"] = (
            "reference_buried_area" if density_den > 0 else None
        )
        summary["lattice_ec_density_bsa_weighted"] = (
            (float(r_bsa) / density_den) if (r_bsa is not None and density_den > 0) else None
        )
        summary["lattice_ec_density_npairs_weighted"] = (
            (float(r_np) / density_den) if (r_np is not None and density_den > 0) else None
        )
        summary["lattice_ec_by_partner_chain"] = by_partner
        summary["lattice_ec_density_denominator_weight_sums"] = {
            "bsa_weight_sum": float(bsa_sum),
            "npairs_weight_sum": float(np_sum),
        }

        return results

    def analyze_interfaces(self, pdb_file, focus_chains=None, reference_chain_id=None):
        results = self.analyze_interfaces_sasa_only(
            pdb_file, focus_chains=focus_chains, reference_chain_id=reference_chain_id
        )
        return self.apply_ec_phase(pdb_file, results, reference_chain_id)

    def _analyse_pairwise_interface(self, chain1, chain2, model):
        """
        Analyse interface between two chains using Biopython wiki approach:
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

        ns_atoms2 = NeighborSearch(atoms2) if (NeighborSearch is not None and atoms2) else None
        ns_atoms1 = NeighborSearch(atoms1) if (NeighborSearch is not None and atoms1) else None
        # Contact residues: NeighborSearch interface residues filtered by the same criteria (distance + type limits)
        contact_residues1 = {
            r for r in interface_residues1 if self._residue_passes_contact_criteria(r, atoms2, ns_atoms2)
        }
        contact_residues2 = {
            r for r in interface_residues2 if self._residue_passes_contact_criteria(r, atoms1, ns_atoms1)
        }

        # If NeighborSearch was unavailable or returned nothing, derive contact residues from the filtered contacts
        if (not contact_residues1 and not contact_residues2) and contacts:
            for c in contacts:
                contact_residues1.add(c['residue1'])
                contact_residues2.add(c['residue2'])

        # Buried surface area (Biopython wiki): always compute for this chain pair
        buried_area = self._calculate_buried_surface_area(chain1, chain2, model)
        results['buried_surface_area'] = buried_area

        # Optional per-residue SASA for accessibility; cached per model across chain pairs
        if self.skip_accessibility_sasa:
            residue_sasa = {}
        else:
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

        # Analyse contact details (on filtered contacts)
        contact_analysis = self._analyse_contacts(contacts)
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
    
    def _residue_passes_contact_criteria(self, residue, other_chain_atoms, ns_other=None):
        """True if residue has at least one atom that, with some atom in other_chain_atoms, satisfies the same distance + type-specific criteria as the contact list."""
        cutoff = float(self.contact_distance)
        if ns_other is not None:
            for atom in residue.get_atoms():
                nearby = ns_other.search(atom.coord, cutoff, 'A')
                if not nearby:
                    continue
                for other in nearby:
                    d = float(np.linalg.norm(atom.coord - other.coord))
                    if d <= cutoff:
                        contact_type = self._classify_contact_type(atom, other)
                        if _contact_within_limit(contact_type, d):
                            return True
            return False
        for atom in residue.get_atoms():
            for other in other_chain_atoms:
                d = float(np.linalg.norm(atom.coord - other.coord))
                if d <= cutoff:
                    contact_type = self._classify_contact_type(atom, other)
                    if _contact_within_limit(contact_type, d):
                        return True
        return False
    
    def _find_contacts(self, atoms1, atoms2):
        """Find all contacts between two sets of atoms (within ``contact_distance``)."""
        contacts = []
        if not atoms1 or not atoms2:
            return contacts
        cutoff = float(self.contact_distance)
        if NeighborSearch is not None:
            ns = NeighborSearch(atoms2)
            for atom1 in atoms1:
                nearby = ns.search(atom1.coord, cutoff, 'A')
                if not nearby:
                    continue
                for atom2 in nearby:
                    distance = float(np.linalg.norm(atom1.coord - atom2.coord))
                    if distance <= cutoff:
                        contacts.append({
                            'atom1': atom1,
                            'atom2': atom2,
                            'distance': distance,
                            'residue1': atom1.parent,
                            'residue2': atom2.parent
                        })
            return contacts
        for atom1 in atoms1:
            for atom2 in atoms2:
                distance = np.linalg.norm(atom1.coord - atom2.coord)
                if distance <= cutoff:
                    contacts.append({
                        'atom1': atom1,
                        'atom2': atom2,
                        'distance': distance,
                        'residue1': atom1.parent,
                        'residue2': atom2.parent
                    })
        return contacts

    def _classify_contact_type(self, atom1, atom2):
        """Classify contact type from atom elements (same scheme as contact_analyser)."""
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

    def _analyse_contacts(self, contacts):
        """Analyse contact statistics."""
        if not contacts:
            return {'average_contact_distance': 0}
        
        distances = [contact['distance'] for contact in contacts]
        
        return {
            'average_contact_distance': np.mean(distances),
            'min_contact_distance': np.min(distances),
            'max_contact_distance': np.max(distances),
            'contact_distance_std': np.std(distances)
        }
    
    def _find_chain_and_model(self, structure, chain_id):
        """Return (chain, model) for the first chain with given id, or (None, None)."""
        for model in structure:
            for chain in model:
                if chain.id == chain_id:
                    return chain, model
        return None, None

    def _sum_residue_sasa_for_chain_in_full_model(self, model, chain_id):
        """
        Sum per-residue SASA (Å²) for chain_id after Shrake–Rupley on a deepcopy
        of the full model (all chains present → mates occlude the probe).

        Reuses ``_model_residue_sasa_cache`` when pairwise analysis already computed
        full-model residue SASA (same geometry), avoiding a second assembly-wide SR pass.
        """
        if not SASA_AVAILABLE or ShrakeRupley is None:
            return None
        # Reuse full-model residue SASA whenever present; independent of
        # ``skip_accessibility_sasa`` (that flag only skips *computing* maps for accessibility).
        mid = id(model)
        with self._sasa_cache_lock:
            cache = getattr(self, "_model_residue_sasa_cache", None)
            if isinstance(cache, dict) and mid in cache:
                full_map = cache[mid]
                if any(cid == chain_id for cid, _ in full_map):
                    return sum(
                        float(v)
                        for (cid, _rid), v in full_map.items()
                        if cid == chain_id
                    )
        try:
            model_copy = copy.deepcopy(model)
        except (RecursionError, TypeError):
            return None
        try:
            sr = ShrakeRupley(probe_radius=1.4, n_points=100)
            sr.compute(model_copy, level='R')
            total = 0.0
            for chain in model_copy.get_chains():
                if chain.id != chain_id:
                    continue
                for residue in chain:
                    s = getattr(residue, 'sasa', None)
                    if s is not None:
                        total += float(s)
            return total
        except Exception:
            return None

    def _lattice_contact_residue_fraction(self, reference_chain, model):
        """
        Fraction of residues in reference_chain that have at least one atom
        within contact_distance of an atom on a different chain (same model).
        """
        if NeighborSearch is None:
            return None
        atoms = list(model.get_atoms())
        if not atoms:
            return None
        ref_cid = reference_chain.id
        try:
            ns = NeighborSearch(atoms)
            contacted = set()
            for atom in reference_chain.get_atoms():
                for nb in ns.search(atom.coord, self.contact_distance, 'A'):
                    if nb.parent.parent.id != ref_cid:
                        contacted.add(atom.parent)
                        break
            res_list = list(reference_chain.get_residues())
            if not res_list:
                return None
            return len(contacted) / len(res_list)
        except Exception:
            return None

    def _lattice_charge_complementarity(self, reference_chain, model):
        """
        Charge complementarity between a reference chain and *all other chains*
        in the same model (multi-copy lattice context).

        This uses the same atom–atom contact definition as pairwise interfaces:
        - atom pairs within contact_distance
        - filtered by type-specific distance limits (H-bond, hydrophobic, etc.)

        Returns a dict with score/opposite/same/total, or None if not computable.
        """
        if NeighborSearch is None:
            return None
        ref_atoms = list(reference_chain.get_atoms())
        if not ref_atoms:
            return None

        ref_cid = reference_chain.id
        other_atoms = []
        for ch in model.get_chains():
            if ch.id == ref_cid:
                continue
            other_atoms.extend(list(ch.get_atoms()))
        if not other_atoms:
            return None

        try:
            ns = NeighborSearch(other_atoms)
            raw_contacts = []
            for a1 in ref_atoms:
                for a2 in ns.search(a1.coord, self.contact_distance, 'A'):
                    d = float(np.linalg.norm(a1.coord - a2.coord))
                    if d <= self.contact_distance:
                        raw_contacts.append({
                            'atom1': a1,
                            'atom2': a2,
                            'distance': d,
                            'residue1': a1.parent,
                            'residue2': a2.parent,
                        })
            if not raw_contacts:
                return {'score': 0.0, 'opposite': 0, 'same': 0, 'total': 0}
            filtered, _type_counts = self._filter_contacts_by_limits(raw_contacts)
            return self._calculate_charge_complementarity(filtered)
        except Exception:
            return None

    def _compute_lattice_charge_metrics(self, reference_chain, model):
        """
        Lattice-wide charge metrics for the reference chain vs all other chains.

        Computes:
        - atom-contact based charge complementarity (existing definition)
        - residue-pair based charge complementarity (unique charged residue pairs)
        - distance-weighted charge complementarity (weights ~ 1/d^2 per unique residue pair, d = min atom distance)
        - partner-resolved versions of the above (per neighbour chain ID)

        Returns JSON-serializable dict, or None if not computable.
        """
        if NeighborSearch is None:
            return None
        ref_atoms = list(reference_chain.get_atoms())
        if not ref_atoms:
            return None

        ref_cid = reference_chain.id
        other_atoms = []
        for ch in model.get_chains():
            if ch.id == ref_cid:
                continue
            other_atoms.extend(list(ch.get_atoms()))
        if not other_atoms:
            return None

        try:
            ns = NeighborSearch(other_atoms)
            raw_contacts = []
            for a1 in ref_atoms:
                for a2 in ns.search(a1.coord, self.contact_distance, 'A'):
                    d = float(np.linalg.norm(a1.coord - a2.coord))
                    if d <= self.contact_distance:
                        raw_contacts.append({
                            'atom1': a1,
                            'atom2': a2,
                            'distance': d,
                            'residue1': a1.parent,
                            'residue2': a2.parent,
                        })
            if not raw_contacts:
                return {
                    'atom': {'score': 0.0, 'opposite': 0, 'same': 0, 'total': 0},
                    'residue_pair': {'score': 0.0, 'opposite': 0, 'same': 0, 'total': 0},
                    'residue_pair_weighted': {'score': 0.0, 'opposite_weight': 0.0, 'same_weight': 0.0, 'total_weight': 0.0},
                    'by_partner_chain': {},
                }

            filtered, _type_counts = self._filter_contacts_by_limits(raw_contacts)

            # --- Atom-contact based score (existing)
            atom_cc = self._calculate_charge_complementarity(filtered)

            # --- Aggregate contacts by charged residue pair and partner chain
            # Key: (ref_res_id, partner_chain_id, partner_res_id)
            pair_min_dist = {}
            pair_charge_prod = {}
            by_partner_pairs = {}

            for c in filtered:
                r1 = c.get('residue1')
                r2 = c.get('residue2')
                if r1 is None or r2 is None:
                    continue
                q1 = self._residue_charge(getattr(r1, 'resname', ''))
                q2 = self._residue_charge(getattr(r2, 'resname', ''))
                if q1 == 0 or q2 == 0:
                    continue
                # partner chain ID for residue2
                try:
                    partner_cid = r2.parent.id
                except Exception:
                    continue
                key = (r1.id, partner_cid, r2.id)
                d = float(c.get('distance', 0.0))
                if key not in pair_min_dist or d < pair_min_dist[key]:
                    pair_min_dist[key] = d
                pair_charge_prod[key] = q1 * q2
                by_partner_pairs.setdefault(partner_cid, set()).add(key)

            # Residue-pair based counts and weighted sums
            rp_opp = 0
            rp_same = 0
            w_opp = 0.0
            w_same = 0.0

            for key, prod in pair_charge_prod.items():
                if prod < 0:
                    rp_opp += 1
                elif prod > 0:
                    rp_same += 1
                d = pair_min_dist.get(key)
                if d is None or d <= 0:
                    continue
                w = 1.0 / (d * d)
                if prod < 0:
                    w_opp += w
                elif prod > 0:
                    w_same += w

            rp_total = rp_opp + rp_same
            rp_score = (rp_opp / float(rp_total)) if rp_total else 0.0
            w_total = w_opp + w_same
            w_score = (w_opp / float(w_total)) if w_total else 0.0

            # Partner-resolved metrics (residue-pair based + weighted)
            by_partner = {}
            for partner_cid, keys in by_partner_pairs.items():
                p_opp = 0
                p_same = 0
                pw_opp = 0.0
                pw_same = 0.0
                for key in keys:
                    prod = pair_charge_prod.get(key, 0)
                    if prod < 0:
                        p_opp += 1
                    elif prod > 0:
                        p_same += 1
                    d = pair_min_dist.get(key)
                    if d is None or d <= 0:
                        continue
                    w = 1.0 / (d * d)
                    if prod < 0:
                        pw_opp += w
                    elif prod > 0:
                        pw_same += w
                p_total = p_opp + p_same
                pw_total = pw_opp + pw_same
                by_partner[str(partner_cid)] = {
                    'residue_pair': {
                        'score': (p_opp / float(p_total)) if p_total else 0.0,
                        'opposite': int(p_opp),
                        'same': int(p_same),
                        'total': int(p_total),
                    },
                    'residue_pair_weighted': {
                        'score': (pw_opp / float(pw_total)) if pw_total else 0.0,
                        'opposite_weight': float(pw_opp),
                        'same_weight': float(pw_same),
                        'total_weight': float(pw_total),
                    },
                }

            return {
                'atom': {
                    'score': float(atom_cc.get('score', 0.0)),
                    'opposite': int(atom_cc.get('opposite', 0)),
                    'same': int(atom_cc.get('same', 0)),
                    'total': int(atom_cc.get('total', 0)),
                },
                'residue_pair': {
                    'score': float(rp_score),
                    'opposite': int(rp_opp),
                    'same': int(rp_same),
                    'total': int(rp_total),
                },
                'residue_pair_weighted': {
                    'score': float(w_score),
                    'opposite_weight': float(w_opp),
                    'same_weight': float(w_same),
                    'total_weight': float(w_total),
                },
                'by_partner_chain': by_partner,
            }
        except Exception:
            return None

    def _compute_lattice_reference_metrics(self, structure, reference_chain_id):
        """
        Empirical lattice-compactness style metrics for one reference chain in a
        multi-copy model (e.g. symmetry-expanded assembly).

        Returns a dict for summary (JSON-serializable floats) or empty dict.
        """
        ref_chain, ref_model = self._find_chain_and_model(structure, reference_chain_id)
        if ref_chain is None or ref_model is None:
            return {
                'lattice_reference_chain': reference_chain_id,
                'lattice_metrics_error': 'reference chain not found in structure',
            }

        out = {'lattice_reference_chain': reference_chain_id}

        ref_res, ref_atoms, ref_mass_da = _polymer_chain_res_atom_mass(ref_chain)
        out['reference_chain_residue_count'] = int(ref_res)
        out['reference_chain_atom_count'] = int(ref_atoms)
        out['reference_chain_mass_Da'] = float(ref_mass_da)
        out['reference_chain_mass_kDa'] = float(ref_mass_da) / 1000.0 if ref_mass_da > 0 else 0.0

        sasa_iso = self._compute_sasa_chain_isolated(ref_chain)
        sasa_cluster = self._sum_residue_sasa_for_chain_in_full_model(ref_model, reference_chain_id)
        out['sasa_reference_isolated'] = sasa_iso
        out['sasa_reference_in_cluster'] = sasa_cluster
        if sasa_iso is not None and sasa_cluster is not None:
            out["reference_chain_BSA"] = float(sasa_iso) - float(sasa_cluster)
        else:
            out["reference_chain_BSA"] = None

        def _norm_sasa_triplet(
            sasa_val: float | None,
            n_res: int,
            n_atoms: int,
            mass_da: float,
            *,
            prefix: str,
        ) -> None:
            if sasa_val is None:
                return
            v = float(sasa_val)
            if n_res > 0:
                out[f"{prefix}_per_residue_reference_chain_A2"] = v / float(n_res)
            if n_atoms > 0:
                out[f"{prefix}_per_atom_reference_chain_A2"] = v / float(n_atoms)
            if mass_da > 0:
                out[f"{prefix}_per_kDa_reference_chain_A2"] = v / (mass_da / 1000.0)

        _norm_sasa_triplet(sasa_iso, ref_res, ref_atoms, ref_mass_da, prefix="sasa_reference_isolated")
        _norm_sasa_triplet(sasa_cluster, ref_res, ref_atoms, ref_mass_da, prefix="sasa_reference_in_cluster")
        _norm_sasa_triplet(out.get("reference_chain_BSA"), ref_res, ref_atoms, ref_mass_da, prefix="reference_chain_BSA")

        if sasa_iso is not None and sasa_iso > 0 and sasa_cluster is not None:
            raw = 1.0 - (float(sasa_cluster) / float(sasa_iso))
            out['lattice_burial_fraction'] = float(max(0.0, min(1.0, raw)))
        else:
            out['lattice_burial_fraction'] = None

        crf = self._lattice_contact_residue_fraction(ref_chain, ref_model)
        out['lattice_contact_residue_fraction'] = crf

        cc = self._compute_lattice_charge_metrics(ref_chain, ref_model)
        out['lattice_charge_metrics'] = cc
        # Backwards-compatible scalar fields (atom-contact based)
        if cc is None:
            out['lattice_charge_complementarity'] = None
            out['lattice_charge_complementarity_opposite'] = None
            out['lattice_charge_complementarity_same'] = None
            out['lattice_charge_complementarity_total'] = None
            out['lattice_charge_complementarity_density'] = None
            out['lattice_charge_complementarity_density_denominator'] = None
        else:
            atom_cc = (cc.get('atom') or {}) if isinstance(cc, dict) else {}
            out['lattice_charge_complementarity'] = float(atom_cc.get('score', 0.0))
            out['lattice_charge_complementarity_opposite'] = int(atom_cc.get('opposite', 0))
            out['lattice_charge_complementarity_same'] = int(atom_cc.get('same', 0))
            out['lattice_charge_complementarity_total'] = int(atom_cc.get('total', 0))

            if sasa_iso is not None and sasa_cluster is not None:
                ref_buried_area = float(sasa_iso) - float(sasa_cluster)
            else:
                ref_buried_area = 0.0
            if ref_buried_area > 0.0:
                out['lattice_charge_complementarity_density'] = float(atom_cc.get('opposite', 0)) / ref_buried_area
                out['lattice_charge_complementarity_density_denominator'] = 'reference_buried_area'
            else:
                out['lattice_charge_complementarity_density'] = None
                out['lattice_charge_complementarity_density_denominator'] = None

        return out

    def _compute_sasa_chain_isolated(self, chain):
        """
        Shrake–Rupley SASA (Å²) for one chain alone (probe 1.4 Å, n_points=100).
        Returns None if SASA is unavailable or computation fails.

        Memoised per chain object during ``analyse_interfaces`` so pairwise BSA does not
        repeat isolated SASA for the same chain across many lattice partners.
        """
        if not SASA_AVAILABLE or ShrakeRupley is None or StructureClass is None or ModelClass is None:
            return None
        oid = id(chain)
        icache_fast = getattr(self, "_chain_isolated_sasa_cache", None)
        if isinstance(icache_fast, dict) and oid in icache_fast:
            return icache_fast[oid]
        with self._sasa_cache_lock:
            icache = getattr(self, "_chain_isolated_sasa_cache", None)
            if not isinstance(icache, dict):
                icache = {}
                self._chain_isolated_sasa_cache = icache
            if oid in icache:
                return icache[oid]
            try:
                c_copy = copy.deepcopy(chain)
            except (RecursionError, TypeError):
                icache[oid] = None
                return None
            try:
                st = StructureClass("single")
                m = ModelClass(0)
                m.add(c_copy)
                st.add(m)
                sr_chain = ShrakeRupley(probe_radius=1.4, n_points=100)
                sr_chain.compute(st, level='S')
                val = getattr(st, 'sasa', None)
                out = float(val) if val is not None else None
                icache[oid] = out
                return out
            except Exception:
                icache[oid] = None
                return None

    def _calculate_buried_surface_area(self, chain1, chain2, _full_model):
        """
        Buried surface area using the Biopython Interface Analysis approach:
        BSA = SASA(chain1 alone) + SASA(chain2 alone) − SASA(complex).

        The complex SASA is computed on a model that contains **only** these
        two chains (deep copies at their current coordinates). That matches
        the standard pairwise definition for each chain pair in a multi-chain
        file (e.g. A–B, A–C, A–D) without folding in burial by other molecules.

        ``_full_model`` is ignored; the Shrake–Rupley path builds a two-chain
        complex only. The parameter is kept so callers can pass the parent model
        unchanged.

        Uses Bio.PDB.SASA.ShrakeRupley when available; otherwise falls back
        to a per-contact-atom estimate (~15 Å² per atom).
        """
        if not SASA_AVAILABLE or ShrakeRupley is None or ModelClass is None:
            return self._estimate_buried_surface_area_fallback(chain1, chain2)

        try:
            c1 = copy.deepcopy(chain1)
            c2 = copy.deepcopy(chain2)
            m_pair = ModelClass(0)
            m_pair.add(c1)
            m_pair.add(c2)
            sr = ShrakeRupley(probe_radius=1.4, n_points=100)
            sr.compute(m_pair, level='S')
            sasa_complex = getattr(m_pair, 'sasa', None)
            if sasa_complex is None:
                return self._estimate_buried_surface_area_fallback(chain1, chain2)

            sasa1 = self._compute_sasa_chain_isolated(chain1)
            sasa2 = self._compute_sasa_chain_isolated(chain2)
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
        complementarity = max(0.0, float(1.0 - avg_deviation / 2.0))  # Normalise
        
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
        One full-model Shrake–Rupley pass per model is cached on ``self._model_residue_sasa_cache``.
        """
        if self.skip_accessibility_sasa:
            return {}
        if not SASA_AVAILABLE or ShrakeRupley is None:
            return {}
        wanted = frozenset(chain_ids)
        try:
            mid = id(model)
            # All cache reads/writes for this model go through one lock so misses cannot race
            # unlocked checks against fills from other threads.
            with self._sasa_cache_lock:
                cache = getattr(self, "_model_residue_sasa_cache", None)
                if cache is None:
                    cache = {}
                    self._model_residue_sasa_cache = cache

                if mid in cache:
                    full_map = cache[mid]
                    return {k: v for k, v in full_map.items() if k[0] in wanted}

                try:
                    model_copy = copy.deepcopy(model)
                except (RecursionError, TypeError):
                    # Memoise failure so other threads / pair jobs skip repeat deepcopy work,
                    # and post-lock reads never assume cache[mid] exists without assignment.
                    cache[mid] = {}
                    return {}

                sr = ShrakeRupley(probe_radius=1.4, n_points=100)
                sr.compute(model_copy, level='R')

                full_map = {}
                for chain in model_copy.get_chains():
                    for residue in chain.get_residues():
                        sasa = getattr(residue, "sasa", None)
                        if sasa is not None:
                            full_map[(chain.id, residue.id)] = float(sasa)
                cache[mid] = full_map
                return {k: v for k, v in full_map.items() if k[0] in wanted}
        except Exception:
            return {}

    def _summarize_polarity_and_accessibility(self, contact_residues1, contact_residues2, residue_sasa, chain1_id, chain2_id):
        """
        Summarise polarity and accessibility for contact residues on both chains.

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
    Patterns are matched as substrings (e.g. 'set_a' and 'set_b' match 'prefix_set_a_set_b.pdb').
    """
    if not patterns:
        return list(paths)
    result = []
    for p in paths:
        name = os.path.basename(p)
        if all(pat in name for pat in patterns):
            result.append(p)
    return result


def _run_analysis(
    analyser,
    paths,
    out_stream,
    focus_chains=None,
    reference_chain_id=None,
    summary_log=None,
    *,
    analyse_hook=None,
    after_results=None,
):
    """Run interface analysis on paths and write results to out_stream.

    ``analyse_hook(analyser, pdb_file) -> results`` overrides the default
    ``analyse_interfaces`` call (used for phased lattice EC runs).

    ``after_results(pdb_file, results)`` runs after a successful report for each file.
    """
    for i, pdb_file in enumerate(paths):
        if len(paths) > 1:
            print(f"\n{'='*50}\n[{i+1}/{len(paths)}] {os.path.basename(pdb_file)}\n{'='*50}", file=out_stream)
        else:
            print(f"Analysing interfaces in {pdb_file}...", file=out_stream)

        if analyse_hook is not None:
            results = analyse_hook(analyser, pdb_file)
        else:
            results = analyser.analyse_interfaces(
                pdb_file, focus_chains=focus_chains, reference_chain_id=reference_chain_id
            )

        if 'error' in results:
            print(f"Error: {results['error']}", file=out_stream)
            try:
                if summary_log is not None:
                    summary_log.error(f"{os.path.basename(pdb_file)}: {results['error']}")
            except Exception:
                pass
            continue

        print("\nINTERFACE ANALYSIS RESULTS:", file=out_stream)
        print("=" * 40, file=out_stream)
        summary = results.get('summary')
        if isinstance(summary, dict):
            # Always print a header for the totals section. In reference-chain mode, print a
            # two-part header (all-chains totals, then reference-chain-scoped totals).
            if reference_chain_id:
                if 'all_total_interfaces' in summary:
                    print("\nAll-chains interface summary (all chain–chain pairs in file):", file=out_stream)
                    print(f"  Total interfaces: {summary.get('all_total_interfaces', 0)}", file=out_stream)
                    print(
                        f"  Total buried surface area: {summary.get('all_total_buried_surface_area', 0):.1f} Å²",
                        file=out_stream,
                    )
                    print(
                        f"  Average buried area per interface: {summary.get('all_average_buried_area_per_interface', 0):.1f} Å²",
                        file=out_stream,
                    )
                    sasa_iso_all = summary.get('all_sasa_isolated_by_chain') or {}
                    sasa_sum_all = summary.get('all_sasa_isolated_sum')
                    if sasa_iso_all and sasa_sum_all is not None:
                        print("  Isolated SASA (chains involved in any all-chains interface):", file=out_stream)
                        for cid in sorted(sasa_iso_all.keys()):
                            print(f"    Chain {cid}: {sasa_iso_all[cid]:.1f} Å²", file=out_stream)
                        print(f"    Sum isolated SASA: {sasa_sum_all:.1f} Å²", file=out_stream)

                print("\nReference-chain interface summary (only interfaces involving reference):", file=out_stream)
            else:
                print("\nInterface summary:", file=out_stream)
            print(f"Total interfaces: {summary.get('total_interfaces', 0)}", file=out_stream)
            print(f"Total buried surface area: {summary.get('total_buried_surface_area', 0):.1f} Å²", file=out_stream)
            print(f"Average buried area per interface: {summary.get('average_buried_area_per_interface', 0):.1f} Å²", file=out_stream)
            sasa_iso = summary.get('sasa_isolated_by_chain') or {}
            sasa_sum = summary.get('sasa_isolated_sum')
            if sasa_iso and sasa_sum is not None:
                print("Isolated SASA (Shrake–Rupley, probe 1.4 Å; each chain alone in solvent):", file=out_stream)
                for cid in sorted(sasa_iso.keys()):
                    print(f"  Chain {cid}: {sasa_iso[cid]:.1f} Å²", file=out_stream)
                print(f"  Sum of isolated SASA (all chains in any reported interface): {sasa_sum:.1f} Å²", file=out_stream)
            elif summary.get('total_interfaces', 0) > 0:
                print("Isolated SASA: N/A (SASA unavailable or computation failed)", file=out_stream)
            err_lat = summary.get('lattice_metrics_error')
            if err_lat:
                print(f"Lattice reference metrics: {err_lat}", file=out_stream)
            elif summary.get('lattice_reference_chain'):
                lc = summary['lattice_reference_chain']
                print(f"Lattice reference chain {lc} (multi-copy / cluster SASA):", file=out_stream)
                sri = summary.get('sasa_reference_isolated')
                src = summary.get('sasa_reference_in_cluster')
                if sri is not None:
                    print(f"  SASA isolated (chain alone): {sri:.1f} Å²", file=out_stream)
                else:
                    print("  SASA isolated: N/A", file=out_stream)
                if src is not None:
                    print(f"  SASA in full model (sum per-residue, chain {lc}): {src:.1f} Å²", file=out_stream)
                else:
                    print("  SASA in full model: N/A", file=out_stream)
                if sri is not None and src is not None:
                    print(f"  Reference-chain BSA (SASA_iso − SASA_cluster): {float(sri) - float(src):.1f} Å²", file=out_stream)

                rr = summary.get("reference_chain_residue_count")
                ra = summary.get("reference_chain_atom_count")
                rm = summary.get("reference_chain_mass_Da")
                if rr is not None and ra is not None and rm is not None:
                    rmk = float(rm) / 1000.0 if float(rm) > 0 else 0.0
                    print(
                        f"  Normalisation divisors - reference chain: residues={int(rr)} atoms={int(ra)} "
                        f"mass={float(rm):.1f} Da ({rmk:.3f} kDa)",
                        file=out_stream,
                    )

                def _fmt_opt(label: str, key: str, unit: str = "Å²") -> None:
                    v = summary.get(key)
                    if v is None:
                        return
                    print(f"  {label}: {float(v):.4g} {unit}", file=out_stream)

                _fmt_opt(
                    "SASA isolated per residue (/ reference chain)",
                    "sasa_reference_isolated_per_residue_reference_chain_A2",
                )
                _fmt_opt(
                    "SASA isolated per atom (/ reference chain)",
                    "sasa_reference_isolated_per_atom_reference_chain_A2",
                )
                _fmt_opt(
                    "SASA isolated per kDa protein (/ reference chain)",
                    "sasa_reference_isolated_per_kDa_reference_chain_A2",
                    unit="Å²/kDa",
                )
                _fmt_opt(
                    "SASA in cluster per residue (/ reference chain)",
                    "sasa_reference_in_cluster_per_residue_reference_chain_A2",
                )
                _fmt_opt(
                    "SASA in cluster per atom (/ reference chain)",
                    "sasa_reference_in_cluster_per_atom_reference_chain_A2",
                )
                _fmt_opt(
                    "SASA in cluster per kDa protein (/ reference chain)",
                    "sasa_reference_in_cluster_per_kDa_reference_chain_A2",
                    unit="Å²/kDa",
                )
                _fmt_opt(
                    "Reference-chain BSA per residue (/ reference chain)",
                    "reference_chain_BSA_per_residue_reference_chain_A2",
                )
                _fmt_opt(
                    "Reference-chain BSA per atom (/ reference chain)",
                    "reference_chain_BSA_per_atom_reference_chain_A2",
                )
                _fmt_opt(
                    "Reference-chain BSA per kDa protein (/ reference chain)",
                    "reference_chain_BSA_per_kDa_reference_chain_A2",
                    unit="Å²/kDa",
                )

                lbf = summary.get('lattice_burial_fraction')
                if lbf is not None:
                    print(f"  Lattice burial fraction (1 − SASA_cluster/SASA_iso): {lbf:.3f} ({100.0*float(lbf):.1f}%)", file=out_stream)
                else:
                    print("  Lattice burial fraction: N/A", file=out_stream)
                crf = summary.get('lattice_contact_residue_fraction')
                if crf is not None:
                    cd = analyser.contact_distance
                    print(f"  Fraction of residues with cross-chain neighbour within {cd:.1f} Å: {crf:.3f} ({100.0*float(crf):.1f}%)", file=out_stream)
                else:
                    print("  Cross-chain contact residue fraction: N/A", file=out_stream)
                # If lattice EC (weighted summaries) is present, omit all charge-based lattice metrics.
                if summary.get('lattice_ec_r_bsa_weighted') is None and summary.get('lattice_ec_r_npairs_weighted') is None:
                    lcc_total = summary.get('lattice_charge_complementarity_total')
                    lcc = summary.get('lattice_charge_complementarity')
                    if lcc is None:
                        print("  Lattice charge complementarity: N/A", file=out_stream)
                    else:
                        if isinstance(lcc_total, (int, float)) and lcc_total > 0:
                            lcc_opp = summary.get('lattice_charge_complementarity_opposite', 0)
                            lcc_same = summary.get('lattice_charge_complementarity_same', 0)
                            print(f"  Lattice charge complementarity: {lcc:.3f} ({100.0*float(lcc):.1f}%)  (charged–charged contacts: {lcc_opp} opposite, {lcc_same} same)", file=out_stream)
                        else:
                            print(f"  Lattice charge complementarity: {lcc:.3f} ({100.0*float(lcc):.1f}%)  (no charged–charged contacts)", file=out_stream)
                        lccd = summary.get('lattice_charge_complementarity_density')
                        if lccd is not None:
                            denom = summary.get('lattice_charge_complementarity_density_denominator') or 'area'
                            print(f"  Lattice charge complementarity density: {lccd:.4f} opposite-sign contacts/Å²  (per {denom})", file=out_stream)
                        # Extended lattice charge metrics (if present)
                        lcm = summary.get('lattice_charge_metrics')
                        if isinstance(lcm, dict):
                            rp = lcm.get('residue_pair') or {}
                            rpw = lcm.get('residue_pair_weighted') or {}
                            if isinstance(rp, dict) and 'score' in rp:
                                rp_score = float(rp.get('score', 0.0))
                                print(f"  Lattice charge complementarity (residue-pair): {rp_score:.3f} ({100.0*rp_score:.1f}%)  (pairs: {int(rp.get('opposite', 0))} opposite, {int(rp.get('same', 0))} same)", file=out_stream)
                            if isinstance(rpw, dict) and 'score' in rpw:
                                rpw_score = float(rpw.get('score', 0.0))
                                print(f"  Lattice charge complementarity (residue-pair, 1/d²-weighted): {rpw_score:.3f} ({100.0*rpw_score:.1f}%)", file=out_stream)
                            by_partner = lcm.get('by_partner_chain') or {}
                            if isinstance(by_partner, dict) and by_partner:
                                print("  Lattice charge complementarity by partner chain (residue-pair):", file=out_stream)
                                for partner_cid in sorted(by_partner.keys()):
                                    rec = by_partner.get(partner_cid) or {}
                                    prp = rec.get('residue_pair') or {}
                                    prpw = rec.get('residue_pair_weighted') or {}
                                    score = float(prp.get('score', 0.0)) if isinstance(prp, dict) else 0.0
                                    opp = int(prp.get('opposite', 0)) if isinstance(prp, dict) else 0
                                    same = int(prp.get('same', 0)) if isinstance(prp, dict) else 0
                                    wscore = float(prpw.get('score', 0.0)) if isinstance(prpw, dict) else 0.0
                                    print(f"    Partner {partner_cid}: score={score:.3f} ({100.0*score:.1f}%) (pairs: {opp} opposite, {same} same)  weighted_score={wscore:.3f} ({100.0*wscore:.1f}%)", file=out_stream)
                # Lattice EC (weighted summaries, if present)
                r_bsa = summary.get('lattice_ec_r_bsa_weighted')
                r_np = summary.get('lattice_ec_r_npairs_weighted')
                if r_bsa is not None or r_np is not None:
                    if r_bsa is not None:
                        print(f"  Lattice EC (r, BSA-weighted Fisher-z): {float(r_bsa):.3f}", file=out_stream)
                    else:
                        print("  Lattice EC (r, BSA-weighted Fisher-z): N/A", file=out_stream)
                    if r_np is not None:
                        tot_np = summary.get('lattice_ec_total_npairs', 0)
                        print(f"  Lattice EC (r, n_pairs-weighted Fisher-z): {float(r_np):.3f}  (total_n_pairs={int(tot_np) if isinstance(tot_np, (int, float)) else 0})", file=out_stream)
                    else:
                        print("  Lattice EC (r, n_pairs-weighted Fisher-z): N/A", file=out_stream)

                    denom = summary.get('lattice_ec_density_denominator') or 'reference_buried_area'
                    d_bsa = summary.get('lattice_ec_density_bsa_weighted')
                    d_np = summary.get('lattice_ec_density_npairs_weighted')
                    if d_bsa is not None:
                        print(f"  Lattice EC density (BSA-weighted): {float(d_bsa):.6f} r/Å²  (per {denom})", file=out_stream)
                    if d_np is not None:
                        print(f"  Lattice EC density (n_pairs-weighted): {float(d_np):.6f} r/Å²  (per {denom})", file=out_stream)

                    byp = summary.get('lattice_ec_by_partner_chain')
                    if byp:
                        print("  Lattice EC by partner chain:", file=out_stream)
                        if isinstance(byp, dict):
                            for partner_cid in sorted(byp.keys()):
                                rec = byp.get(partner_cid) or {}
                                r = rec.get('ec_r')
                                n = rec.get('ec_n_pairs', 0)
                                if r is None:
                                    print(f"    Partner {partner_cid}: r=N/A (n_pairs={int(n) if isinstance(n, (int, float)) else 0})", file=out_stream)
                                else:
                                    print(f"    Partner {partner_cid}: r={float(r):.3f} (n_pairs={int(n) if isinstance(n, (int, float)) else 0})", file=out_stream)
                        elif isinstance(byp, list):
                            for rec in byp:
                                if not isinstance(rec, dict):
                                    continue
                                partner_cid = rec.get('partner_chain_id', '?')
                                r = rec.get('ec_r')
                                n = rec.get('ec_n_pairs', 0)
                                if r is None:
                                    print(f"    Partner {partner_cid}: r=N/A (n_pairs={int(n) if isinstance(n, (int, float)) else 0})", file=out_stream)
                                else:
                                    print(f"    Partner {partner_cid}: r={float(r):.3f} (n_pairs={int(n) if isinstance(n, (int, float)) else 0})", file=out_stream)
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
                    # If EC is present, omit all charge-based interface metrics.
                    if interface.get('ec_r') is None:
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
                    # EC per interface (if present)
                    e_r = interface.get('ec_r')
                    e_n = interface.get('ec_n_pairs')
                    if e_r is not None:
                        print(f"  EC (r): {float(e_r):.3f}  (n_pairs={int(e_n) if isinstance(e_n, (int, float)) else 0})", file=out_stream)
                        e_d = interface.get('ec_density')
                        if e_d is not None:
                            denom = interface.get('ec_density_denominator') or 'buried_surface_area'
                            print(f"  EC density: {float(e_d):.6f} r/Å²  (per {denom})", file=out_stream)
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

        if after_results is not None:
            after_results(pdb_file, results)

    if len(paths) > 1:
        print(f"\nDone. Processed {len(paths)} file(s).", file=out_stream)


def main():
    """CLI: single file, multiple files, directory, or glob patterns; optional multi-set filtering and dry-run."""
    import argparse

    parser = argparse.ArgumentParser(
        description="Analyse interfaces between molecules in crystal structures.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""Examples (repository root; replace paths, filenames, and chain IDs):
  # ASU interface metrics (charge-tag; single structure)
  python metrics/interface_analyser_asu_charge.py /path/to/project/model_01.pdb

  # ASU interface metrics (charge-tag; batch; write merged text report)
  python metrics/interface_analyser_asu_charge.py /path/to/project/models/*.pdb -o /path/to/project/results_charge.txt

  # ASU interface metrics (charge-tag; batch; set filtering)
  python metrics/interface_analyser_asu_charge.py /path/to/project/models/*.pdb --set set_a,set_b -o /path/to/project/results_charge_set_a_set_b.txt

  # ASU interface metrics (charge-tag; per-set output)
  python metrics/interface_analyser_asu_charge.py /path/to/project/models/*.pdb --sets set_a set_b set_c -o "/path/to/project/results_charge_{}.txt"

  # ASU interface metrics (charge-tag; per-structure output)
  python metrics/interface_analyser_asu_charge.py /path/to/project/models/*.pdb --per-structure -o "/path/to/project/out/{}_interfaces_charge.txt"

  # Dry run (no analysis; list resolved inputs)
  python metrics/interface_analyser_asu_charge.py /path/to/project/models/*.pdb --sets set_a set_b -o "/path/to/project/results_charge_{}.txt" --dry-run

  # Lattice interface metrics (multi-copy model; reference-chain report)
  # Note: when --reference-chain is set, the report header totals (Total interfaces, total buried surface area,
  # isolated SASA list and sum) are scoped to the reference chain and its partner chains.
  python metrics/interface_analyser_lattice_charge.py /path/to/project/expanded_assembly.pdb --reference-chain A -o /path/to/project/lattice_charge.txt
  python metrics/interface_analyser_lattice_ec.py /path/to/project/expanded_assembly.pdb --reference-chain A -o /path/to/project/lattice_ec.txt

  # ASU interface metrics (electrostatic complementarity, EC)
  python metrics/interface_analyser_asu_ec.py /path/to/project/model_01.pdb -o /path/to/project/results_ec.txt

Interface report CSV extractors (charge vs EC):
  # Charge reports (ASU or lattice) → per-interface rows; optional parsed summary rows
  python metrics/interface_mol_report_charge_csv.py /path/to/project/results_charge*.txt -m A -o /path/to/project/charge_interfaces_A.csv

  # EC reports (ASU or lattice)
  # - Summary only (one row per structure)
  python metrics/interface_mol_report_ec_csv.py /path/to/project/results_ec*.txt -m A --mode summary -o /path/to/project/ec_summary.csv

  # - Summary columns followed by per-interface columns (interfaces involving chain A)
  python metrics/interface_mol_report_ec_csv.py /path/to/project/results_ec*.txt -m A -o /path/to/project/ec_summary_and_interfaces_A.csv

Legacy (interfaces only): python metrics/interface_molecule_report_csv.py results.txt -m A -o interfaces_A.csv
""",
    )
    parser.add_argument(
        'input',
        nargs='+',
        help='PDB/CIF file(s), directory, or glob pattern (e.g. *.pdb).',
    )
    parser.add_argument(
        '--set', '-s',
        action='append',
        dest='sets',
        metavar='PATTERNS',
        help='Comma-separated patterns; file is included only if basename contains ALL patterns. Repeat for multiple sets (e.g. -s set_a,set_b -s set_c,set_d).',
    )
    parser.add_argument(
        '--sets',
        dest='sets_multi',
        nargs='+',
        metavar='SET',
        help='Multiple set names in one go (one pattern per set). E.g. --sets set_a set_b set_c is equivalent to --set set_a --set set_b --set set_c.',
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
        help='Print which files would be processed for each set and exit without running. Writes summary to interface_analyser_dryrun.txt in the current directory.',
    )
    parser.add_argument(
        '--chains',
        metavar='IDS',
        help="Focus analysis on specific chain IDs (comma-separated). Only chain pairs where at least one chain is in this list are analysed. Example: --chains A or --chains A,B",
    )
    parser.add_argument(
        '--reference-chain',
        metavar='ID',
        dest='reference_chain_id',
        help='For symmetry-expanded / multi-copy PDBs: chain ID of the focal molecule (e.g. central copy). Reports isolated SASA vs sum of per-residue SASA in the full model, lattice_burial_fraction, and fraction of residues with a cross-chain neighbour within contact_distance.',
    )
    args = parser.parse_args()

    focus_chains = None
    if getattr(args, 'chains', None):
        focus_chains = [c.strip() for c in str(args.chains).split(',') if c.strip()]
    reference_chain_id = getattr(args, 'reference_chain_id', None)
    if reference_chain_id:
        reference_chain_id = str(reference_chain_id).strip() or None

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
            dryrun_file = os.path.join(os.getcwd(), "interface_analyser_dryrun.txt")
            try:
                with open(dryrun_file, "w") as f:
                    f.write(out_text)
                print(out_text, file=sys.stderr, flush=True)
                print(f"Dry-run summary written to {dryrun_file}", file=sys.stderr, flush=True)
            except OSError as e:
                print(out_text, file=sys.stderr, flush=True)
                print(f"Could not write {dryrun_file}: {e}", file=sys.stderr, flush=True)
            return
        analyser = InterfaceAnalyser()
        for stem, [path] in structure_list:
            out_path = args.output.replace('{}', stem)
            print(f"Writing {os.path.basename(path)} -> {out_path}", file=sys.stderr)
            with open(out_path, 'w') as out:
                _run_analysis(
                    analyser, [path], out,
                    focus_chains=focus_chains,
                    reference_chain_id=reference_chain_id,
                )
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
                "Run from a directory that contains the structure files, or pass explicit paths/directories.",
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
        dryrun_file = os.path.join(os.getcwd(), "interface_analyser_dryrun.txt")
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
        print("Run from a directory that contains the structure files, or pass explicit paths/directories.", file=sys.stderr)
        sys.exit(1)

    # Single output file for all sets, or one file per set when -o contains '{}'
    single_output_file = args.output and '{}' not in args.output
    if single_output_file:
        out = open(args.output, 'w')
        print(f"Writing all sets to {args.output}", file=sys.stderr)

    analyser = InterfaceAnalyser()
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
            _run_analysis(
                analyser, filtered, out,
                focus_chains=focus_chains,
                reference_chain_id=reference_chain_id,
            )
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