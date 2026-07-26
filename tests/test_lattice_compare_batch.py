#!/usr/bin/env python3
"""Tests for metrics/lattice_compare_batch.py."""

from __future__ import annotations

import contextlib
import csv
import io
import json
import os
import sys
import tempfile
import unittest

_REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def _ec_txt_common(*, basename_stem: str) -> str:
    """Minimal lattice EC report text parseable by parse_ec_report_text."""
    pdb_name = f"{basename_stem}.pdb"
    return f"""Analysing interfaces in /fake/path/{pdb_name}...

INTERFACE ANALYSIS RESULTS:
========================================

Reference-chain interface summary (only interfaces involving reference):
Total interfaces: 1
Total buried surface area: 60.0 Å²
Average buried area per interface: 60.0 Å²

Lattice reference chain A (multi-copy / cluster SASA):
  SASA isolated (chain alone): 5000.0 Å²
  SASA in full model (sum per-residue, chain A): 4000.0 Å²
  Reference-chain BSA (SASA_iso − SASA_cluster): 1000.0 Å²
  Normalisation divisors - reference chain: residues=10 atoms=100 mass=5000.0 Da (5.000 kDa)
  Reference-chain BSA per residue (/ reference chain): 100.0000 Å²
  Reference-chain BSA per atom (/ reference chain): 10.0000 Å²
  Reference-chain BSA per kDa protein (/ reference chain): 200.0000 Å²/kDa
  Lattice burial fraction (1 − SASA_cluster/SASA_iso): 0.200 (20.0%)
  Fraction of residues with cross-chain neighbour within 5.0 Å: 0.300 (30.0%)
  Lattice EC (r, BSA-weighted Fisher-z): 0.850
  Lattice EC (r, n_pairs-weighted Fisher-z): 0.900 (total_n_pairs=100)
  Lattice EC density (BSA-weighted): 0.001 r/Å² (per reference_buried_area)

Interface 1: A-B
    Contact count (within limits): 42
    Buried surface area: 60.0 Å² Contact area: 30.0 Å²
    EC (r): 0.70 (n_pairs=50)
    EC density: 0.012 r/Å² (per pairwise_BSA)
"""


class LatticeCompareBatchTests(unittest.TestCase):
    def test_join_matrices_and_combined_scaling(self) -> None:
        sys.path.insert(0, _REPO)
        from metrics import lattice_compare_batch as lcb  # noqa: PLC0415

        with tempfile.TemporaryDirectory() as td:
            for stem in ("alpha", "beta"):
                pj = os.path.join(td, f"{stem}_lattice_packing.json")
                with open(pj, "w", encoding="utf-8") as f:
                    json.dump(
                        {
                            "input": f"/data/{stem}.pdb",
                            "lattice_atom_density_atoms_per_a3": 0.002,
                            "lattice_mass_density_da_per_a3": 0.003,
                            "lattice_packing_density_percent": 35.0,
                            "lattice_matthews_a3_per_da": 2.5,
                            "estimated_solvent_content_percent": 42.0,
                        },
                        f,
                    )
                rt = os.path.join(td, f"{stem}_ec.txt")
                with open(rt, "w", encoding="utf-8") as f:
                    f.write(_ec_txt_common(basename_stem=stem))

            rc = lcb.main(
                [
                    "--packing",
                    td,
                    "--ec-reports",
                    os.path.join(td, "*.txt"),
                    "--reference-chain",
                    "A",
                    "--output-dir",
                    td,
                    "--skip-heatmaps",
                    "--no-default-heatmap-annotations",
                ]
            )
            self.assertEqual(rc, 0)

            combined = os.path.join(td, "combined_lattice_vs_ec.csv")
            self.assertTrue(os.path.isfile(combined))
            with open(combined, newline="", encoding="utf-8") as f:
                rdr = csv.DictReader(f)
                header = rdr.fieldnames or []
                rows = list(rdr)
            banned = {"structure_basename", "set_label", "source_json", "ec_report_txt"}
            self.assertTrue(all(x not in header for x in banned))
            self.assertEqual(header[0], "structure_stem")
            self.assertEqual(header[1], "total_interfaces")
            tail = [
                "lattice_reference_chain",
                "atom_density_per_1000_A3",
                "mass_density_Da_per_1000_A3",
                "packing_density_percent",
                "matthews_a3_per_Da",
                "estimated_solvent_percent",
            ]
            self.assertEqual(header[-len(tail) :], tail)
            self.assertEqual(len(rows), 2)
            by_stem = {r["structure_stem"]: r for r in rows}
            self.assertAlmostEqual(float(by_stem["alpha"]["atom_density_per_1000_A3"]), 2.0)
            self.assertAlmostEqual(float(by_stem["alpha"]["mass_density_Da_per_1000_A3"]), 3.0)
            self.assertAlmostEqual(float(by_stem["alpha"]["lattice_ec_density_bsa_weighted_per_1000_A2"]), 1.0)
            self.assertEqual(int(by_stem["beta"]["total_interfaces"]), 1)

            mat_bsa = os.path.join(td, "lattice_compare_matrix_bsa.csv")
            self.assertTrue(os.path.isfile(mat_bsa))
            with open(mat_bsa, newline="", encoding="utf-8") as f:
                bsa_rows = list(csv.DictReader(f))
            # Default: structures on x-axis (transpose); CSV rows are interfaces.
            self.assertEqual(len(bsa_rows), 1)
            self.assertEqual(bsa_rows[0]["chain_pair_canonical"], "A-B")
            self.assertAlmostEqual(float(bsa_rows[0]["alpha.pdb"]), 60.0)
            self.assertAlmostEqual(float(bsa_rows[0]["beta.pdb"]), 60.0)

            afull = os.path.join(td, "alpha_full.csv")
            self.assertTrue(os.path.isfile(afull))
            self.assertFalse(os.path.isfile(os.path.join(td, "alpha_packing_full.csv")))
            self.assertFalse(os.path.isfile(os.path.join(td, "alpha_lattice_ec_full.csv")))
            with open(afull, newline="", encoding="utf-8") as f:
                mr = csv.DictReader(f)
                mh = list(mr.fieldnames or [])
                mrows = list(mr)
            self.assertEqual(mh[0], "section")
            self.assertIn("structure_stem", mh)
            self.assertNotIn("ec_density_percent", mh)
            self.assertNotIn("distance_min_A", mh)
            sections = [r["section"] for r in mrows]
            self.assertIn("packing", sections)
            self.assertIn("ec_summary", sections)
            self.assertIn("ec_interface", sections)

    def test_strict_match_errors_on_mismatch(self) -> None:
        sys.path.insert(0, _REPO)
        from metrics import lattice_compare_batch as lcb  # noqa: PLC0415

        with tempfile.TemporaryDirectory() as td:
            for stem in ("alpha", "beta"):
                pj = os.path.join(td, f"{stem}_lattice_packing.json")
                with open(pj, "w", encoding="utf-8") as f:
                    json.dump({"input": f"/data/{stem}.pdb"}, f)
            with open(os.path.join(td, "ec_alpha.txt"), "w", encoding="utf-8") as f:
                f.write(_ec_txt_common(basename_stem="alpha"))

            with contextlib.redirect_stderr(io.StringIO()):
                with self.assertRaises(SystemExit):
                    lcb.main(
                        [
                            "--packing",
                            td,
                            "--ec-reports",
                            os.path.join(td, "ec_alpha.txt"),
                            "--reference-chain",
                            "A",
                            "--output-dir",
                            td,
                            "--strict-match",
                            "--skip-heatmaps",
                            "--no-default-heatmap-annotations",
                        ]
                    )


if __name__ == "__main__":
    unittest.main()
