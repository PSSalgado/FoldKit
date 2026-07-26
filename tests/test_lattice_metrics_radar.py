#!/usr/bin/env python3
"""Tests for metrics/lattice_metrics_radar.py."""

from __future__ import annotations

import contextlib
import csv
import io
import math
import os
import sys
import tempfile
import unittest

_REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if _REPO not in sys.path:
    sys.path.insert(0, _REPO)


_METRIC_COLUMNS = (
    "packing_density_percent",
    "matthews_a3_per_Da",
    "estimated_solvent_percent",
    "reference_chain_BSA_per_kDa_reference_chain_A2",
    "lattice_burial_fraction_percent",
    "lattice_contact_residue_fraction_percent",
)


def _write_combined_csv(path: str, rows: list[dict[str, object]], *, drop_columns=()) -> None:
    fieldnames = ["structure_stem", *[c for c in _METRIC_COLUMNS if c not in set(drop_columns)]]
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames, extrasaction="ignore")
        w.writeheader()
        for r in rows:
            w.writerow(r)


def _sample_rows() -> list[dict[str, object]]:
    return [
        {
            "structure_stem": "alpha",
            "packing_density_percent": 30.0,
            "matthews_a3_per_Da": 2.0,
            "estimated_solvent_percent": 40.0,
            "reference_chain_BSA_per_kDa_reference_chain_A2": 150.0,
            "lattice_burial_fraction_percent": 20.0,
            "lattice_contact_residue_fraction_percent": 25.0,
        },
        {
            "structure_stem": "beta",
            "packing_density_percent": 40.0,
            "matthews_a3_per_Da": 3.0,
            "estimated_solvent_percent": 55.0,
            "reference_chain_BSA_per_kDa_reference_chain_A2": 250.0,
            "lattice_burial_fraction_percent": 35.0,
            "lattice_contact_residue_fraction_percent": 45.0,
        },
    ]


class ScoreLogicTests(unittest.TestCase):
    def test_min_max_map_to_0_and_10(self) -> None:
        from metrics import lattice_metrics_radar as lmr  # noqa: PLC0415

        self.assertAlmostEqual(lmr.score_metric(30.0, 30.0, 40.0, invert=False), 0.0)
        self.assertAlmostEqual(lmr.score_metric(40.0, 30.0, 40.0, invert=False), 10.0)
        self.assertAlmostEqual(lmr.score_metric(35.0, 30.0, 40.0, invert=False), 5.0)

    def test_custom_score_range_and_inversion(self) -> None:
        from metrics import lattice_metrics_radar as lmr  # noqa: PLC0415

        score_range = (1.0, 5.0)
        self.assertAlmostEqual(
            lmr.score_metric(30.0, 30.0, 40.0, invert=False, score_range=score_range),
            1.0,
        )
        self.assertAlmostEqual(
            lmr.score_metric(40.0, 30.0, 40.0, invert=False, score_range=score_range),
            5.0,
        )
        self.assertAlmostEqual(
            lmr.score_metric(35.0, 30.0, 40.0, invert=False, score_range=score_range),
            3.0,
        )
        self.assertAlmostEqual(
            lmr.score_metric(30.0, 30.0, 40.0, invert=True, score_range=score_range),
            5.0,
        )

    def test_invalid_score_range_errors(self) -> None:
        from metrics import lattice_metrics_radar as lmr  # noqa: PLC0415

        with self.assertRaises(ValueError):
            lmr.score_metric(5.0, 0.0, 10.0, invert=False, score_range=(5.0, 5.0))
        with self.assertRaises(ValueError):
            lmr.compute_scores(_sample_rows(), ["alpha", "beta"], score_range=(10.0, 0.0))

    def test_inversion_flips_direction(self) -> None:
        from metrics import lattice_metrics_radar as lmr  # noqa: PLC0415

        # Higher Matthews (looser) must yield a lower plot score.
        low = lmr.score_metric(2.0, 2.0, 3.0, invert=True)
        high = lmr.score_metric(3.0, 2.0, 3.0, invert=True)
        self.assertAlmostEqual(low, 10.0)
        self.assertAlmostEqual(high, 0.0)
        self.assertGreater(low, high)

    def test_equal_cohort_maps_to_mid(self) -> None:
        from metrics import lattice_metrics_radar as lmr  # noqa: PLC0415

        self.assertAlmostEqual(lmr.score_metric(5.0, 5.0, 5.0, invert=False), 5.0)
        self.assertAlmostEqual(lmr.score_metric(5.0, 5.0, 5.0, invert=True), 5.0)

    def test_nan_passthrough(self) -> None:
        from metrics import lattice_metrics_radar as lmr  # noqa: PLC0415

        self.assertTrue(math.isnan(lmr.score_metric(math.nan, 0.0, 10.0, invert=False)))

    def test_compute_scores_matthews_inverted_across_cohort(self) -> None:
        from metrics import lattice_metrics_radar as lmr  # noqa: PLC0415

        rows = _sample_rows()
        scores_by_stem, raw_by_stem, kept, _limits = lmr.compute_scores(rows, ["alpha", "beta"])
        self.assertEqual(kept, ["alpha", "beta"])
        matthews_idx = [m.column for m in lmr.METRICS].index("matthews_a3_per_Da")
        # alpha has lower Matthews (tighter) -> higher inverted score than beta.
        self.assertAlmostEqual(scores_by_stem["alpha"][matthews_idx], 10.0)
        self.assertAlmostEqual(scores_by_stem["beta"][matthews_idx], 0.0)
        packing_idx = [m.column for m in lmr.METRICS].index("packing_density_percent")
        self.assertAlmostEqual(scores_by_stem["alpha"][packing_idx], 0.0)
        self.assertAlmostEqual(scores_by_stem["beta"][packing_idx], 10.0)
        self.assertAlmostEqual(raw_by_stem["alpha"][matthews_idx], 2.0)

    def test_compute_scores_uses_custom_range(self) -> None:
        from metrics import lattice_metrics_radar as lmr  # noqa: PLC0415

        scores, _raw, kept, _limits = lmr.compute_scores(
            _sample_rows(), ["alpha", "beta"], score_range=(0.0, 100.0)
        )
        packing_idx = [m.column for m in lmr.METRICS].index("packing_density_percent")
        self.assertEqual(kept, ["alpha", "beta"])
        self.assertAlmostEqual(scores["alpha"][packing_idx], 0.0)
        self.assertAlmostEqual(scores["beta"][packing_idx], 100.0)


class VolumeScaleTests(unittest.TestCase):
    def test_round_outward_steps(self) -> None:
        from metrics import lattice_metrics_radar as lmr  # noqa: PLC0415

        self.assertEqual(lmr._round_outward_to_0_or_5(5.1, 13.6), (5.0, 15.0))
        self.assertEqual(lmr._round_outward_to_0_or_5(78.2, 94.1), (75.0, 95.0))
        self.assertEqual(lmr._round_outward_to_0_or_5(6.0, 22.0), (5.0, 25.0))
        # Extended Matthews uses 0.5 steps.
        self.assertEqual(lmr._round_outward_to_half(5.1, 13.6), (5.0, 14.0))
        self.assertEqual(lmr._round_outward_to_half(5.4, 14.4), (5.0, 14.5))

    def test_cryst1_limits_and_crystal_alias(self) -> None:
        from metrics import lattice_metrics_radar as lmr  # noqa: PLC0415

        limits, meta = lmr.resolve_volume_limits(scale_mode="cryst1", matthews_values=[2.5, 2.8])
        self.assertEqual(limits["matthews_a3_per_Da"], (1.5, 4.0))
        self.assertEqual(limits["estimated_solvent_percent"], (25.0, 80.0))
        self.assertEqual(limits["packing_density_percent"], (20.0, 75.0))
        self.assertEqual(meta["scale_mode"], "cryst1")
        _, meta2 = lmr.resolve_volume_limits(scale_mode="crystal", matthews_values=[2.5])
        self.assertEqual(meta2["scale_mode"], "cryst1")

    def test_bbox_volume_profile(self) -> None:
        from metrics import lattice_metrics_radar as lmr  # noqa: PLC0415

        limits, meta = lmr.resolve_volume_limits(scale_mode="bbox", matthews_values=[9.0])
        self.assertEqual(limits["matthews_a3_per_Da"], (8.0, 45.0))
        self.assertEqual(limits["estimated_solvent_percent"], (80.0, 100.0))
        self.assertEqual(limits["packing_density_percent"], (2.0, 15.0))
        self.assertEqual(meta["scale_mode"], "bbox")
        self.assertEqual(meta["volume_profile"]["version"], 1)

    def test_retired_slayer_mode_is_rejected(self) -> None:
        from metrics import lattice_metrics_radar as lmr  # noqa: PLC0415

        self.assertNotIn("slayer", lmr.SCALE_MODES)
        with self.assertRaises(ValueError):
            lmr.resolve_volume_limits(scale_mode="slayer", matthews_values=[9.0])
        with self.assertRaises(ValueError):
            lmr.resolve_occlusion_limits(occlusion_scale_mode="slayer")

    def test_slayer_compact_volume_profile(self) -> None:
        from metrics import lattice_metrics_radar as lmr  # noqa: PLC0415

        limits, meta = lmr.resolve_volume_limits(
            scale_mode="slayer-compact", matthews_values=[9.0]
        )
        # Compact SlpA bbox2 calibration.
        self.assertEqual(limits["matthews_a3_per_Da"], (7.0, 10.0))
        self.assertEqual(limits["estimated_solvent_percent"], (80.0, 90.0))
        self.assertEqual(limits["packing_density_percent"], (10.0, 15.0))
        self.assertEqual(meta["volume_profile"]["version"], 1)

    def test_volume_empirical_close_and_wide(self) -> None:
        from metrics import lattice_metrics_radar as lmr  # noqa: PLC0415

        limits, meta = lmr.resolve_volume_limits(
            scale_mode="empirical",
            matthews_values=[8.0, 9.0],
            packing_values=[12.0, 14.0],
            solvent_values=[84.0, 86.0],
        )
        self.assertEqual(limits["matthews_a3_per_Da"], (7.0, 10.0))
        self.assertEqual(limits["packing_density_percent"], (11.0, 15.0))
        self.assertEqual(limits["estimated_solvent_percent"], (83.0, 87.0))
        self.assertEqual(
            meta["empirical_bands"]["matthews_a3_per_Da"]["rule"], "close_integer_pad"
        )

        limits2, meta2 = lmr.resolve_volume_limits(
            scale_mode="empirical",
            matthews_values=[8.6, 44.2],
            packing_values=[2.4, 12.4],
            solvent_values=[85.7, 97.2],
        )
        self.assertEqual(limits2["matthews_a3_per_Da"], (8.5, 44.5))
        self.assertEqual(limits2["packing_density_percent"], (2.0, 15.0))
        self.assertEqual(limits2["estimated_solvent_percent"], (85.0, 100.0))
        self.assertEqual(meta2["empirical_bands"]["matthews_a3_per_Da"]["rule"], "wide_grid")

    def test_user_requires_matthews(self) -> None:
        from metrics import lattice_metrics_radar as lmr  # noqa: PLC0415

        with self.assertRaises(ValueError):
            lmr.resolve_volume_limits(scale_mode="user", matthews_values=[8.0])

    def test_user_derives_solvent_and_packing(self) -> None:
        from metrics import lattice_metrics_radar as lmr  # noqa: PLC0415

        limits, meta = lmr.resolve_volume_limits(
            scale_mode="user",
            matthews_values=[8.0],
            matthews_min=5.0,
            matthews_max=15.0,
        )
        self.assertEqual(meta["scale_mode"], "user")
        self.assertEqual(limits["matthews_a3_per_Da"], (5.0, 15.0))
        sol_lo, sol_hi = limits["estimated_solvent_percent"]
        self.assertGreater(sol_hi, sol_lo)
        pack_lo, pack_hi = limits["packing_density_percent"]
        self.assertGreater(pack_hi, pack_lo)
        for value in (sol_lo, sol_hi, pack_lo, pack_hi):
            self.assertEqual(value % 5, 0)
        raw_sol = meta["limits_before_round"]["estimated_solvent_percent"]
        self.assertLessEqual(sol_lo, raw_sol[0])
        self.assertGreaterEqual(sol_hi, raw_sol[1])

    def test_explicit_limits_are_not_rounded(self) -> None:
        from metrics import lattice_metrics_radar as lmr  # noqa: PLC0415

        limits, meta = lmr.resolve_volume_limits(
            scale_mode="cryst1",
            matthews_values=[8.0],
            matthews_min=5.1,
            matthews_max=13.6,
            solvent_min=78.2,
            solvent_max=94.1,
        )
        self.assertEqual(limits["matthews_a3_per_Da"], (5.1, 13.6))
        self.assertEqual(limits["estimated_solvent_percent"], (78.2, 94.1))
        self.assertEqual(
            meta["explicit_limits"],
            ["estimated_solvent_percent", "matthews_a3_per_Da"],
        )
        self.assertEqual(limits["packing_density_percent"], (20.0, 75.0))

    def test_volume_limits_used_in_scores(self) -> None:
        from metrics import lattice_metrics_radar as lmr  # noqa: PLC0415

        rows = [
            {
                "structure_stem": "a",
                "packing_density_percent": 12.0,
                "matthews_a3_per_Da": 8.0,
                "estimated_solvent_percent": 85.0,
                "reference_chain_BSA_per_kDa_reference_chain_A2": 50.0,
                "lattice_burial_fraction_percent": 10.0,
                "lattice_contact_residue_fraction_percent": 12.0,
            },
            {
                "structure_stem": "b",
                "packing_density_percent": 14.0,
                "matthews_a3_per_Da": 9.0,
                "estimated_solvent_percent": 87.0,
                "reference_chain_BSA_per_kDa_reference_chain_A2": 60.0,
                "lattice_burial_fraction_percent": 15.0,
                "lattice_contact_residue_fraction_percent": 18.0,
            },
        ]
        limits, _ = lmr.resolve_volume_limits(scale_mode="bbox", matthews_values=[8.0, 9.0])
        scores, _, kept, used = lmr.compute_scores(rows, ["a", "b"], volume_limits=limits)
        self.assertEqual(kept, ["a", "b"])
        self.assertEqual(used["matthews_a3_per_Da"], limits["matthews_a3_per_Da"])
        bsa_idx = [m.column for m in lmr.METRICS].index(
            "reference_chain_BSA_per_kDa_reference_chain_A2"
        )
        self.assertAlmostEqual(scores["a"][bsa_idx], 0.0)
        self.assertAlmostEqual(scores["b"][bsa_idx], 10.0)

    def test_occlusion_bbox_uses_profile_and_cryst1_falls_back(self) -> None:
        from metrics import lattice_metrics_radar as lmr  # noqa: PLC0415

        values = {
            "lattice_burial_fraction_percent": [9.0, 12.0],
            "lattice_contact_residue_fraction_percent": [12.0, 15.0],
            "reference_chain_BSA_per_kDa_reference_chain_A2": [43.0, 57.0],
        }
        bbox_lim, bbox_meta = lmr.resolve_occlusion_limits(
            occlusion_scale_mode="bbox", occlusion_values=values
        )
        # Fixed profile bands, not the within-run span.
        self.assertEqual(bbox_lim["lattice_burial_fraction_percent"], (2.0, 20.0))
        self.assertEqual(bbox_meta["occlusion_profile"]["profile"], "bbox")

        cryst_lim, cryst_meta = lmr.resolve_occlusion_limits(
            occlusion_scale_mode="cryst1", occlusion_values=values
        )
        self.assertEqual(cryst_lim["lattice_burial_fraction_percent"], (8.0, 13.0))
        self.assertIn("empirical", cryst_meta.get("occlusion_note", "").lower())


def _slpa_like_rows() -> list[dict[str, object]]:
    """Compact SlpA 15-molecule cohort: occlusion values close enough to collapse under min-max."""
    data = [
        ("SlpA2", 12.0, 9.0, 86.0, 56.5, 12.0, 15.0),
        ("SlpA7", 14.0, 8.0, 84.0, 56.9, 12.0, 15.0),
        ("SlpA7b", 12.0, 9.0, 86.0, 43.1, 9.0, 12.0),
        ("SlpA11", 12.0, 9.0, 86.0, 45.2, 10.0, 13.0),
        ("SlpA4d", 13.0, 8.0, 85.0, 48.0, 10.0, 12.0),
    ]
    return [dict(zip(("structure_stem", *_METRIC_COLUMNS), row)) for row in data]


class OcclusionProfileTests(unittest.TestCase):
    def test_shipped_bbox_profile_bands(self) -> None:
        from metrics import lattice_metrics_radar as lmr  # noqa: PLC0415

        limits, meta = lmr.resolve_occlusion_limits(occlusion_scale_mode="bbox")
        self.assertEqual(limits["lattice_burial_fraction_percent"], (2.0, 20.0))
        self.assertEqual(limits["lattice_contact_residue_fraction_percent"], (2.0, 25.0))
        self.assertEqual(
            limits["reference_chain_BSA_per_kDa_reference_chain_A2"], (5.0, 80.0)
        )
        profile = meta["occlusion_profile"]
        self.assertEqual(profile["profile"], "bbox")
        self.assertEqual(profile["version"], 1)
        # Provenance: observed span must sit inside the rounded band.
        obs = profile["observed"]["lattice_burial_fraction_percent"]
        self.assertGreaterEqual(obs[0], 0.0)
        self.assertLessEqual(obs[1], 20.0)
        self.assertTrue(profile["sources"])

    def test_shipped_slayer_compact_profile_bands(self) -> None:
        from metrics import lattice_metrics_radar as lmr  # noqa: PLC0415

        limits, meta = lmr.resolve_occlusion_limits(occlusion_scale_mode="slayer-compact")
        self.assertEqual(limits["lattice_burial_fraction_percent"], (5.0, 15.0))
        self.assertEqual(limits["lattice_contact_residue_fraction_percent"], (10.0, 20.0))
        self.assertEqual(
            limits["reference_chain_BSA_per_kDa_reference_chain_A2"], (40.0, 60.0)
        )
        profile = meta["occlusion_profile"]
        self.assertEqual(profile["profile"], "slayer-compact")
        self.assertEqual(profile["version"], 1)

    def test_profile_bands_sit_on_the_declared_grids(self) -> None:
        from metrics import lattice_metrics_radar as lmr  # noqa: PLC0415

        for mode in ("bbox", "slayer-compact"):
            limits, _ = lmr.resolve_occlusion_limits(occlusion_scale_mode=mode)
            for col, (lo, hi) in limits.items():
                self.assertEqual(
                    lo % lmr.occlusion_lower_round_step(col),
                    0,
                    f"{mode}/{col} lower limit off-grid",
                )
                self.assertEqual(
                    hi % lmr.occlusion_round_step(col),
                    0,
                    f"{mode}/{col} upper limit off-grid",
                )

    def test_profile_filename_slug(self) -> None:
        from metrics import lattice_metrics_radar as lmr  # noqa: PLC0415

        self.assertEqual(lmr.profile_filename("bbox", 1), "bbox_v1.json")
        self.assertEqual(lmr.profile_filename("slayer-compact", 1), "slayer_compact_v1.json")

    def test_missing_profile_errors(self) -> None:
        from metrics import lattice_metrics_radar as lmr  # noqa: PLC0415

        with tempfile.TemporaryDirectory() as td:
            with self.assertRaises(ValueError) as ctx:
                lmr.load_occlusion_profile("bbox", profile_dir=td)
            self.assertIn("calibrate_occlusion_profile", str(ctx.exception))

    def test_unknown_mode_errors(self) -> None:
        from metrics import lattice_metrics_radar as lmr  # noqa: PLC0415

        with self.assertRaises(ValueError):
            lmr.resolve_occlusion_limits(occlusion_scale_mode="nonsense")

    def test_cohort_mode_has_no_fixed_limits(self) -> None:
        from metrics import lattice_metrics_radar as lmr  # noqa: PLC0415

        limits, meta = lmr.resolve_occlusion_limits(occlusion_scale_mode="cohort")
        self.assertEqual(limits, {})
        self.assertNotIn("occlusion_profile", meta)

    def test_user_mode_requires_all_three_pairs(self) -> None:
        from metrics import lattice_metrics_radar as lmr  # noqa: PLC0415

        with self.assertRaises(ValueError) as ctx:
            lmr.resolve_occlusion_limits(
                occlusion_scale_mode="user", loi_min=0.0, loi_max=15.0
            )
        self.assertIn("--bsa-kda-min", str(ctx.exception))

    def test_user_limits_stay_exact(self) -> None:
        from metrics import lattice_metrics_radar as lmr  # noqa: PLC0415

        limits, meta = lmr.resolve_occlusion_limits(
            occlusion_scale_mode="user",
            loi_min=2.3,
            loi_max=16.4,
            contact_min=2.7,
            contact_max=20.2,
            bsa_kda_min=9.9,
            bsa_kda_max=76.8,
        )
        self.assertEqual(limits["lattice_burial_fraction_percent"], (2.3, 16.4))
        self.assertEqual(limits["lattice_contact_residue_fraction_percent"], (2.7, 20.2))
        self.assertEqual(
            limits["reference_chain_BSA_per_kDa_reference_chain_A2"], (9.9, 76.8)
        )
        self.assertEqual(len(meta["explicit_occlusion_limits"]), 3)

    def test_override_replaces_single_profile_axis(self) -> None:
        from metrics import lattice_metrics_radar as lmr  # noqa: PLC0415

        limits, meta = lmr.resolve_occlusion_limits(
            occlusion_scale_mode="bbox", loi_min=1.0, loi_max=18.0
        )
        self.assertEqual(limits["lattice_burial_fraction_percent"], (1.0, 18.0))
        # Untouched axes keep the profile band.
        self.assertEqual(limits["lattice_contact_residue_fraction_percent"], (2.0, 25.0))
        self.assertEqual(meta["explicit_occlusion_limits"], ["lattice_burial_fraction_percent"])

    def test_half_pair_and_inverted_pair_error(self) -> None:
        from metrics import lattice_metrics_radar as lmr  # noqa: PLC0415

        with self.assertRaises(ValueError) as ctx:
            lmr.resolve_occlusion_limits(occlusion_scale_mode="bbox", contact_min=5.0)
        self.assertIn("--contact-min", str(ctx.exception))
        with self.assertRaises(ValueError):
            lmr.resolve_occlusion_limits(
                occlusion_scale_mode="bbox", bsa_kda_min=80.0, bsa_kda_max=40.0
            )


class EmpiricalBandTests(unittest.TestCase):
    def test_close_span_pads_by_one_integer(self) -> None:
        from metrics import lattice_metrics_radar as lmr  # noqa: PLC0415

        lo, hi, rule = lmr.empirical_occlusion_band([9.0, 12.0, 10.0])
        self.assertEqual((lo, hi), (8.0, 13.0))
        self.assertEqual(rule, "close_integer_pad")

        lo, hi, rule = lmr.empirical_occlusion_band([9.3, 12.1])
        self.assertEqual((lo, hi), (8.0, 14.0))
        self.assertEqual(rule, "close_integer_pad")

    def test_wide_span_rounds_to_0_or_5(self) -> None:
        from metrics import lattice_metrics_radar as lmr  # noqa: PLC0415

        lo, hi, rule = lmr.empirical_occlusion_band([9.86, 76.84])
        self.assertEqual((lo, hi), (5.0, 80.0))
        self.assertEqual(rule, "wide_0_or_5")

        lo, hi, rule = lmr.empirical_occlusion_band(
            [2.3, 16.0], column="lattice_burial_fraction_percent"
        )
        self.assertEqual((lo, hi), (2.0, 20.0))
        self.assertEqual(rule, "wide_0_or_5")

    def test_boundary_span_of_5_is_close(self) -> None:
        from metrics import lattice_metrics_radar as lmr  # noqa: PLC0415

        lo, hi, rule = lmr.empirical_occlusion_band([10.0, 15.0])
        self.assertEqual(rule, "close_integer_pad")
        self.assertEqual((lo, hi), (9.0, 16.0))

    def test_negative_lower_limit_clamped(self) -> None:
        from metrics import lattice_metrics_radar as lmr  # noqa: PLC0415

        lo, hi, rule = lmr.empirical_occlusion_band([0.2, 1.5])
        self.assertEqual(rule, "close_integer_pad")
        self.assertEqual(lo, 0.0)
        self.assertEqual(hi, 3.0)

    def test_resolve_empirical_uses_per_metric_rules(self) -> None:
        from metrics import lattice_metrics_radar as lmr  # noqa: PLC0415

        values = {
            "lattice_burial_fraction_percent": [9.0, 12.0, 10.0],
            "lattice_contact_residue_fraction_percent": [12.0, 15.0],
            "reference_chain_BSA_per_kDa_reference_chain_A2": [9.9, 76.8],
        }
        limits, meta = lmr.resolve_occlusion_limits(
            occlusion_scale_mode="empirical", occlusion_values=values
        )
        self.assertEqual(limits["lattice_burial_fraction_percent"], (8.0, 13.0))
        self.assertEqual(limits["lattice_contact_residue_fraction_percent"], (11.0, 16.0))
        self.assertEqual(
            limits["reference_chain_BSA_per_kDa_reference_chain_A2"], (5.0, 80.0)
        )
        bands = meta["empirical_bands"]
        self.assertEqual(bands["lattice_burial_fraction_percent"]["rule"], "close_integer_pad")
        self.assertEqual(
            bands["reference_chain_BSA_per_kDa_reference_chain_A2"]["rule"], "wide_0_or_5"
        )

    def test_empirical_end_to_end_avoids_collapse(self) -> None:
        from metrics import lattice_metrics_radar as lmr  # noqa: PLC0415

        with tempfile.TemporaryDirectory() as td:
            csv_path = os.path.join(td, "combined_lattice_vs_ec.csv")
            _write_combined_csv(csv_path, _slpa_like_rows())
            out = os.path.join(td, "radar")
            with contextlib.redirect_stderr(io.StringIO()):
                lmr.generate_radar(
                    input_csv=csv_path,
                    output_dir=out,
                    occlusion_scale_mode="empirical",
                )
            import json  # noqa: PLC0415

            with open(
                os.path.join(out + "_vol-bbox_occ-empirical", "radar_scale_vol-bbox_occ-empirical.json"),
                encoding="utf-8",
            ) as f:
                meta = json.load(f)
            self.assertEqual(meta["occlusion_scale_mode"], "empirical")
            # SlpA LOI 9–12 → 8–13; contact 12–15 → 11–16; BSA 43.1–56.9 → 40–60
            self.assertEqual(meta["used_limits"]["lattice_burial_fraction_percent"], [8.0, 13.0])
            self.assertEqual(
                meta["used_limits"]["lattice_contact_residue_fraction_percent"], [11.0, 16.0]
            )
            self.assertEqual(
                meta["used_limits"]["reference_chain_BSA_per_kDa_reference_chain_A2"],
                [40.0, 60.0],
            )


class OcclusionScoringTests(unittest.TestCase):
    def test_bbox_avoids_0_10_collapse_that_cohort_produces(self) -> None:
        from metrics import lattice_metrics_radar as lmr  # noqa: PLC0415

        rows = _slpa_like_rows()
        stems = [str(r["structure_stem"]) for r in rows]
        occl_idx = [
            [m.column for m in lmr.METRICS].index(col) for col in lmr.OCCLUSION_COLUMNS
        ]

        cohort_scores, _, _, _ = lmr.compute_scores(rows, stems)
        cohort_vals = [cohort_scores[s][i] for s in stems for i in occl_idx]
        self.assertIn(0.0, cohort_vals)
        self.assertIn(10.0, cohort_vals)

        limits, _ = lmr.resolve_occlusion_limits(occlusion_scale_mode="bbox")
        bbox_scores, _, _, used = lmr.compute_scores(rows, stems, occlusion_limits=limits)
        for stem in stems:
            for i in occl_idx:
                score = bbox_scores[stem][i]
                self.assertGreater(score, 0.0)
                self.assertLess(score, 10.0)
        for col, band in limits.items():
            self.assertEqual(used[col], band)

    def test_score_is_clamped_to_the_ring(self) -> None:
        from metrics import lattice_metrics_radar as lmr  # noqa: PLC0415

        self.assertAlmostEqual(lmr.score_metric(30.0, 0.0, 20.0, invert=False), 10.0)
        self.assertAlmostEqual(lmr.score_metric(-5.0, 0.0, 20.0, invert=False), 0.0)
        self.assertAlmostEqual(lmr.score_metric(30.0, 0.0, 20.0, invert=True), 0.0)

    def test_clipped_values_reports_out_of_band_only(self) -> None:
        from metrics import lattice_metrics_radar as lmr  # noqa: PLC0415

        rows = _slpa_like_rows()
        stems = [str(r["structure_stem"]) for r in rows]
        limits, _ = lmr.resolve_occlusion_limits(occlusion_scale_mode="slayer-compact")
        _scores, raw_by_stem, kept, used = lmr.compute_scores(
            rows, stems, occlusion_limits=limits
        )
        # The compact cohort is what the compact band was calibrated on, so nothing clips.
        self.assertEqual(lmr.clipped_values(kept, raw_by_stem, used), [])

        # A sparse / low-occlusion lattice falls below the compact band and must be reported.
        open_row = dict(zip(("structure_stem", *_METRIC_COLUMNS),
                            ("loose_sparse", 12.0, 9.0, 86.0, 9.9, 2.3, 2.7)))
        rows2 = [*rows, open_row]
        stems2 = [*stems, "loose_sparse"]
        _s2, raw2, kept2, used2 = lmr.compute_scores(rows2, stems2, occlusion_limits=limits)
        clipped2 = lmr.clipped_values(kept2, raw2, used2)
        offenders = {(c["structure_stem"], c["metric"], c["side"]) for c in clipped2}
        self.assertIn(
            ("loose_sparse", "reference_chain_BSA_per_kDa_reference_chain_A2", "below"),
            offenders,
        )
        self.assertIn(("loose_sparse", "lattice_burial_fraction_percent", "below"), offenders)
        self.assertTrue(all(c["structure_stem"] == "loose_sparse" for c in clipped2))

    def test_default_mode_is_bbox_and_recorded_in_scale_json(self) -> None:
        import json  # noqa: PLC0415

        from metrics import lattice_metrics_radar as lmr  # noqa: PLC0415

        with tempfile.TemporaryDirectory() as td:
            csv_path = os.path.join(td, "combined_lattice_vs_ec.csv")
            _write_combined_csv(csv_path, _slpa_like_rows())
            out = os.path.join(td, "radar")
            with contextlib.redirect_stderr(io.StringIO()):
                lmr.generate_radar(input_csv=csv_path, output_dir=out)
            with open(
                os.path.join(out + "_bbox", "radar_scale_bbox.json"), encoding="utf-8"
            ) as f:
                meta = json.load(f)
            self.assertEqual(meta["scale_mode"], "bbox")
            self.assertEqual(meta["occlusion_scale_mode"], "bbox")
            self.assertEqual(meta["occlusion_profile"]["profile"], "bbox")
            self.assertEqual(meta["volume_profile"]["profile"], "bbox")
            # Sample VM=8.0 sits on the bbox volume floor (v2); occlusion stays in-band.
            occl_clip = [
                c for c in meta["clipped"] if c["metric"] in lmr.OCCLUSION_COLUMNS
            ]
            self.assertEqual(occl_clip, [])
            self.assertEqual(
                meta["used_limits"]["lattice_burial_fraction_percent"], [2.0, 20.0]
            )
            self.assertEqual(meta["used_limits"]["matthews_a3_per_Da"], [8.0, 45.0])
            self.assertEqual(meta["used_limits"]["estimated_solvent_percent"], [80.0, 100.0])

    def test_cli_accepts_occlusion_mode_and_overrides(self) -> None:
        import json  # noqa: PLC0415

        from metrics import lattice_metrics_radar as lmr  # noqa: PLC0415

        with tempfile.TemporaryDirectory() as td:
            csv_path = os.path.join(td, "combined_lattice_vs_ec.csv")
            _write_combined_csv(csv_path, _slpa_like_rows())
            out = os.path.join(td, "radar")
            with contextlib.redirect_stdout(io.StringIO()), contextlib.redirect_stderr(
                io.StringIO()
            ):
                rc = lmr.main(
                    [
                        "--input", csv_path,
                        "--output-dir", out,
                        "--occlusion-scale-mode", "cohort",
                        "--loi-min", "1.5",
                        "--loi-max", "18.5",
                    ]
                )
            self.assertEqual(rc, 0)
            with open(
                os.path.join(out + "_vol-bbox_occ-cohort", "radar_scale_vol-bbox_occ-cohort.json"),
                encoding="utf-8",
            ) as f:
                meta = json.load(f)
            self.assertEqual(meta["occlusion_scale_mode"], "cohort")
            self.assertEqual(
                meta["used_limits"]["lattice_burial_fraction_percent"], [1.5, 18.5]
            )
            # Contact % had no override and no profile, so it stayed cohort-scaled (9-12 -> 12-15).
            self.assertEqual(
                meta["used_limits"]["lattice_contact_residue_fraction_percent"], [12.0, 15.0]
            )

    def test_cli_user_mode_without_all_pairs_exits_nonzero(self) -> None:
        from metrics import lattice_metrics_radar as lmr  # noqa: PLC0415

        with tempfile.TemporaryDirectory() as td:
            csv_path = os.path.join(td, "combined_lattice_vs_ec.csv")
            _write_combined_csv(csv_path, _slpa_like_rows())
            with contextlib.redirect_stderr(io.StringIO()):
                with self.assertRaises(SystemExit):
                    lmr.main(
                        [
                            "--input", csv_path,
                            "--output-dir", os.path.join(td, "radar"),
                            "--occlusion-scale-mode", "user",
                            "--loi-min", "0",
                            "--loi-max", "20",
                        ]
                    )


class CalibrationHelperTests(unittest.TestCase):
    def test_build_profile_rounds_span_outward(self) -> None:
        from metrics import calibrate_occlusion_profile as cal  # noqa: PLC0415

        with tempfile.TemporaryDirectory() as td:
            csv_path = os.path.join(td, "combined_lattice_vs_ec.csv")
            _write_combined_csv(csv_path, _slpa_like_rows())
            doc = cal.build_profile(
                [csv_path], profile="test", version=2, population="unit test cohort"
            )
            self.assertEqual(doc["profile"], "test")
            self.assertEqual(doc["version"], 2)
            self.assertEqual(doc["n_structures"], 5)
            loi = doc["metric_limits"]["lattice_burial_fraction_percent"]
            self.assertEqual(loi["observed_min"], 9.0)
            self.assertEqual(loi["observed_max"], 12.0)
            self.assertEqual(loi["limits"], [5.0, 15.0])
            self.assertEqual(loi["lower_round_step"], 5.0)
            self.assertEqual(loi["upper_round_step"], 5.0)
            bsa = doc["metric_limits"]["reference_chain_BSA_per_kDa_reference_chain_A2"]
            self.assertEqual(bsa["limits"], [40.0, 60.0])

    def test_stems_subset_and_unknown_stem_errors(self) -> None:
        from metrics import calibrate_occlusion_profile as cal  # noqa: PLC0415

        with tempfile.TemporaryDirectory() as td:
            csv_path = os.path.join(td, "combined_lattice_vs_ec.csv")
            _write_combined_csv(csv_path, _slpa_like_rows())
            doc = cal.build_profile(
                [csv_path], profile="subset", version=1, stems=["SlpA2", "SlpA7"]
            )
            self.assertEqual(doc["n_structures"], 2)
            loi = doc["metric_limits"]["lattice_burial_fraction_percent"]
            self.assertEqual((loi["observed_min"], loi["observed_max"]), (12.0, 12.0))
            with self.assertRaises(ValueError):
                cal.build_profile([csv_path], profile="x", version=1, stems=["nope"])

    def test_p5p95_trims_the_band(self) -> None:
        from metrics import calibrate_occlusion_profile as cal  # noqa: PLC0415

        with tempfile.TemporaryDirectory() as td:
            csv_path = os.path.join(td, "combined_lattice_vs_ec.csv")
            _write_combined_csv(csv_path, _slpa_like_rows())
            doc = cal.build_profile(
                [csv_path], profile="trim", version=1, band_method="p5p95"
            )
            loi = doc["metric_limits"]["lattice_burial_fraction_percent"]
            self.assertEqual(doc["band_method"], "p5p95")
            self.assertGreaterEqual(loi["band_min_unrounded"], loi["observed_min"])
            self.assertLessEqual(loi["band_max_unrounded"], loi["observed_max"])

    def test_written_profile_is_loadable_by_the_radar_module(self) -> None:
        from metrics import calibrate_occlusion_profile as cal  # noqa: PLC0415
        from metrics import lattice_metrics_radar as lmr  # noqa: PLC0415

        with tempfile.TemporaryDirectory() as td:
            csv_path = os.path.join(td, "combined_lattice_vs_ec.csv")
            _write_combined_csv(csv_path, _slpa_like_rows())
            doc = cal.build_profile([csv_path], profile="roundtrip", version=3)
            path = cal.write_profile(doc, output_dir=td)
            self.assertEqual(os.path.basename(path), "roundtrip_v3.json")
            loaded = lmr.load_occlusion_profile("roundtrip", version=3, profile_dir=td)
            self.assertEqual(loaded["profile"], "roundtrip")
            self.assertEqual(
                loaded["metric_limits"]["lattice_burial_fraction_percent"]["limits"],
                [5.0, 15.0],
            )

    def test_shipped_profiles_match_a_fresh_recalibration(self) -> None:
        """Guard against hand-edited profile JSONs drifting from the calibration rule."""
        from metrics import calibrate_occlusion_profile as cal  # noqa: PLC0415
        from metrics import lattice_metrics_radar as lmr  # noqa: PLC0415

        for mode in ("bbox", "slayer-compact"):
            doc = lmr.load_occlusion_profile(mode)
            structures = doc["structures"]
            for col in lmr.OCCLUSION_COLUMNS:
                vals = [
                    float(s[col]) for s in structures if s.get(col) is not None
                ]
                lower_step = doc["metric_limits"][col].get(
                    "lower_round_step",
                    doc["metric_limits"][col].get("round_step"),
                )
                upper_step = doc["metric_limits"][col].get(
                    "upper_round_step",
                    doc["metric_limits"][col].get("round_step"),
                )
                expected = lmr._round_outward_asymmetric(
                    min(vals),
                    max(vals),
                    lower_step=lower_step,
                    upper_step=upper_step,
                )
                self.assertEqual(
                    tuple(doc["metric_limits"][col]["limits"]),
                    expected,
                    f"{mode}/{col} band does not match its recorded structures",
                )
            self.assertEqual(doc["band_method"], "minmax")
            self.assertTrue(doc["notes"], f"{mode} profile has no provenance notes")
            self.assertIn(doc["band_method"], cal._BAND_METHODS)


_TRANSPOSED_LABELS = (
    "Packing density (%)",
    "Matthews-like coefficient (Å3/Da)",
    "Estimated solvent (%)",
    "BSAmolA/kDa (Å²)",
    "Lattice occlusion index (LOImolA) (%)",
    "Mol. A lattice contact residues (%)",
)


def _write_transposed_csv(path: str, *, drop_labels=(), extra_rows=()) -> None:
    """Metrics as rows, structures as columns (as in metrics_summary_table.csv)."""
    rows = [["", "alpha", "beta"]]
    rows.extend(extra_rows)
    values = {
        "Packing density (%)": [30.0, 40.0],
        "Matthews-like coefficient (Å3/Da)": [2.0, 3.0],
        "Estimated solvent (%)": [40.0, 55.0],
        "BSAmolA/kDa (Å²)": [150.0, 250.0],
        "Lattice occlusion index (LOImolA) (%)": [20.0, 35.0],
        "Mol. A lattice contact residues (%)": [25.0, 45.0],
    }
    for label in _TRANSPOSED_LABELS:
        if label in set(drop_labels):
            continue
        rows.append([label, *values[label]])
    with open(path, "w", newline="", encoding="utf-8") as f:
        csv.writer(f).writerows(rows)


class OrientationTests(unittest.TestCase):
    def test_detects_and_loads_transposed_table(self) -> None:
        from metrics import lattice_metrics_radar as lmr  # noqa: PLC0415

        with tempfile.TemporaryDirectory() as td:
            p = os.path.join(td, "metrics_summary_table.csv")
            _write_transposed_csv(
                p,
                extra_rows=[
                    ["Atom density (n./1000 Å3)", 9, 10],
                    ["Total volume (1000 Å3)", 9095, 7895],
                ],
            )
            rows = lmr.load_rows(p)
            self.assertEqual([r["structure_stem"] for r in rows], ["alpha", "beta"])
            self.assertEqual(float(rows[0]["packing_density_percent"]), 30.0)
            self.assertEqual(float(rows[1]["matthews_a3_per_Da"]), 3.0)
            self.assertEqual(
                float(rows[1]["reference_chain_BSA_per_kDa_reference_chain_A2"]), 250.0
            )

    def test_both_orientations_give_identical_scores(self) -> None:
        from metrics import lattice_metrics_radar as lmr  # noqa: PLC0415

        with tempfile.TemporaryDirectory() as td:
            wide = os.path.join(td, "combined_lattice_vs_ec.csv")
            trans = os.path.join(td, "metrics_summary_table.csv")
            _write_combined_csv(wide, _sample_rows())
            _write_transposed_csv(trans)

            w_scores, _, w_kept, _ = lmr.compute_scores(lmr.load_rows(wide), ["alpha", "beta"])
            t_scores, _, t_kept, _ = lmr.compute_scores(lmr.load_rows(trans), ["alpha", "beta"])
            self.assertEqual(w_kept, t_kept)
            for stem in w_kept:
                for a, b in zip(w_scores[stem], t_scores[stem]):
                    self.assertAlmostEqual(a, b)

    def test_transposed_missing_metric_row_errors(self) -> None:
        from metrics import lattice_metrics_radar as lmr  # noqa: PLC0415

        with tempfile.TemporaryDirectory() as td:
            p = os.path.join(td, "metrics_summary_table.csv")
            _write_transposed_csv(p, drop_labels=("Lattice occlusion index (LOImolA) (%)",))
            with self.assertRaises(ValueError) as ctx:
                lmr.load_rows(p)
            self.assertIn("LOImolA", str(ctx.exception))

    def test_unrecognised_layout_errors(self) -> None:
        from metrics import lattice_metrics_radar as lmr  # noqa: PLC0415

        with tempfile.TemporaryDirectory() as td:
            p = os.path.join(td, "junk.csv")
            with open(p, "w", newline="", encoding="utf-8") as f:
                csv.writer(f).writerows([["a", "b"], ["1", "2"]])
            with self.assertRaises(ValueError):
                lmr.load_rows(p)

    def test_transposed_end_to_end_outputs(self) -> None:
        from metrics import lattice_metrics_radar as lmr  # noqa: PLC0415

        with tempfile.TemporaryDirectory() as td:
            p = os.path.join(td, "metrics_summary_table.csv")
            _write_transposed_csv(p)
            out = os.path.join(td, "radar")
            lmr.generate_radar(input_csv=p, output_dir=out)
            tagged = out + "_bbox"
            self.assertTrue(os.path.isfile(os.path.join(tagged, "radar_grid_bbox.png")))
            self.assertTrue(
                os.path.isfile(os.path.join(tagged, "radar_scores_heatmap_bbox.png"))
            )


class OutputTests(unittest.TestCase):
    def test_custom_score_range_is_written_and_rendered(self) -> None:
        import json  # noqa: PLC0415

        from metrics import lattice_metrics_radar as lmr  # noqa: PLC0415

        with tempfile.TemporaryDirectory() as td:
            csv_path = os.path.join(td, "combined_lattice_vs_ec.csv")
            _write_combined_csv(csv_path, _sample_rows())
            out = os.path.join(td, "radar")
            with contextlib.redirect_stderr(io.StringIO()):
                written = lmr.generate_radar(
                    input_csv=csv_path,
                    output_dir=out,
                    score_range=(1.0, 5.0),
                )

            tagged = out + "_bbox"
            scale_path = os.path.join(tagged, "radar_scale_bbox.json")
            scores_path = os.path.join(tagged, "radar_scores_bbox.csv")
            self.assertIn(scale_path, written)
            self.assertTrue(os.path.isfile(os.path.join(tagged, "radar_grid_bbox.png")))
            self.assertTrue(
                os.path.isfile(os.path.join(tagged, "radar_scores_heatmap_bbox.png"))
            )
            with open(scale_path, encoding="utf-8") as f:
                meta = json.load(f)
            self.assertEqual(meta["score_range"], [1.0, 5.0])
            with open(scores_path, newline="", encoding="utf-8") as f:
                score_rows = list(csv.DictReader(f))
            score_columns = [c for c in (score_rows[0].keys() if score_rows else []) if c.endswith("__score")]
            self.assertTrue(score_columns)
            for row in score_rows:
                for column in score_columns:
                    value = float(row[column])
                    self.assertGreaterEqual(value, 1.0)
                    self.assertLessEqual(value, 5.0)

    def test_cli_accepts_custom_score_range(self) -> None:
        import json  # noqa: PLC0415

        from metrics import lattice_metrics_radar as lmr  # noqa: PLC0415

        with tempfile.TemporaryDirectory() as td:
            csv_path = os.path.join(td, "combined_lattice_vs_ec.csv")
            _write_combined_csv(csv_path, _sample_rows())
            out = os.path.join(td, "radar")
            with contextlib.redirect_stdout(io.StringIO()), contextlib.redirect_stderr(
                io.StringIO()
            ):
                rc = lmr.main(
                    [
                        "--input",
                        csv_path,
                        "--output-dir",
                        out,
                        "--score-range",
                        "0",
                        "5",
                    ]
                )
            self.assertEqual(rc, 0)
            with open(
                os.path.join(out + "_bbox", "radar_scale_bbox.json"),
                encoding="utf-8",
            ) as f:
                meta = json.load(f)
            self.assertEqual(meta["score_range"], [0.0, 5.0])

    def test_grid_heatmap_and_scores_csv_written(self) -> None:
        from metrics import lattice_metrics_radar as lmr  # noqa: PLC0415

        with tempfile.TemporaryDirectory() as td:
            csv_path = os.path.join(td, "combined_lattice_vs_ec.csv")
            _write_combined_csv(csv_path, _sample_rows())
            out = os.path.join(td, "radar")
            written = lmr.generate_radar(input_csv=csv_path, output_dir=out)

            tagged = out + "_bbox"
            grid = os.path.join(tagged, "radar_grid_bbox.png")
            heat = os.path.join(tagged, "radar_scores_heatmap_bbox.png")
            scores = os.path.join(tagged, "radar_scores_bbox.csv")
            self.assertTrue(os.path.isfile(grid))
            self.assertTrue(os.path.isfile(heat))
            self.assertTrue(os.path.isfile(scores))
            self.assertIn(grid, written)
            self.assertIn(heat, written)
            # No overlay artefact should ever be produced.
            self.assertFalse(os.path.isfile(os.path.join(tagged, "radar_overlay.png")))

            with open(scores, newline="", encoding="utf-8") as f:
                rdr = csv.DictReader(f)
                header = rdr.fieldnames or []
                srows = list(rdr)
            self.assertEqual(header[0], "structure_stem")
            self.assertIn("matthews_a3_per_Da__score", header)
            self.assertIn("matthews_a3_per_Da__raw", header)
            self.assertEqual({r["structure_stem"] for r in srows}, {"alpha", "beta"})

    def test_no_heatmap_skips_score_heatmap(self) -> None:
        from metrics import lattice_metrics_radar as lmr  # noqa: PLC0415

        with tempfile.TemporaryDirectory() as td:
            csv_path = os.path.join(td, "combined_lattice_vs_ec.csv")
            _write_combined_csv(csv_path, _sample_rows())
            out = os.path.join(td, "radar")
            with contextlib.redirect_stdout(io.StringIO()), contextlib.redirect_stderr(
                io.StringIO()
            ):
                written = lmr.generate_radar(
                    input_csv=csv_path,
                    output_dir=out,
                    write_heatmap=False,
                )
                rc = lmr.main(
                    [
                        "--input",
                        csv_path,
                        "--output-dir",
                        os.path.join(td, "radar_cli"),
                        "--no-heatmap",
                        "--no-scale-tag",
                    ]
                )
            self.assertEqual(rc, 0)
            tagged = out + "_bbox"
            heat = os.path.join(tagged, "radar_scores_heatmap_bbox.png")
            self.assertFalse(os.path.isfile(heat))
            self.assertTrue(any(p.endswith("radar_grid_bbox.png") for p in written))
            self.assertFalse(any("heatmap" in os.path.basename(p) for p in written))
            self.assertFalse(
                os.path.isfile(os.path.join(td, "radar_cli", "radar_scores_heatmap.png"))
            )
            self.assertTrue(os.path.isfile(os.path.join(td, "radar_cli", "radar_grid.png")))

    def test_format_and_max_per_page(self) -> None:
        from metrics import lattice_metrics_radar as lmr  # noqa: PLC0415

        with tempfile.TemporaryDirectory() as td:
            csv_path = os.path.join(td, "combined_lattice_vs_ec.csv")
            _write_combined_csv(csv_path, _sample_rows())
            out = os.path.join(td, "radar")
            lmr.generate_radar(
                input_csv=csv_path,
                output_dir=out,
                out_format="svg",
                max_per_page=1,
            )
            tagged = out + "_bbox"
            self.assertTrue(os.path.isfile(os.path.join(tagged, "radar_grid_bbox_01.svg")))
            self.assertTrue(os.path.isfile(os.path.join(tagged, "radar_grid_bbox_02.svg")))
            self.assertTrue(
                os.path.isfile(os.path.join(tagged, "radar_scores_heatmap_bbox.svg"))
            )
            self.assertFalse(os.path.isfile(os.path.join(tagged, "radar_grid_bbox.svg")))

    def test_reference_stem_missing_errors(self) -> None:
        from metrics import lattice_metrics_radar as lmr  # noqa: PLC0415

        with tempfile.TemporaryDirectory() as td:
            csv_path = os.path.join(td, "combined_lattice_vs_ec.csv")
            _write_combined_csv(csv_path, _sample_rows())
            with self.assertRaises(ValueError):
                lmr.generate_radar(
                    input_csv=csv_path,
                    output_dir=os.path.join(td, "radar"),
                    reference_stem="does_not_exist",
                )


class GhostAndDeviationTests(unittest.TestCase):
    """Leave-one-out ghost and cohort-deviation display mode."""

    def _three_rows(self) -> list[dict[str, object]]:
        # Odd cohort so the shared median lands on a real structure, which is
        # exactly the coincidence leave-one-out is meant to avoid. Spacing is
        # deliberately uneven so the LOO mean of the extremes is not equal to
        # the middle value either.
        return [
            {
                "structure_stem": "low",
                "packing_density_percent": 10.0,
                "matthews_a3_per_Da": 10.5,
                "estimated_solvent_percent": 89.0,
                "reference_chain_BSA_per_kDa_reference_chain_A2": 40.0,
                "lattice_burial_fraction_percent": 8.0,
                "lattice_contact_residue_fraction_percent": 10.0,
            },
            {
                "structure_stem": "mid",
                "packing_density_percent": 11.0,
                "matthews_a3_per_Da": 9.0,
                "estimated_solvent_percent": 86.0,
                "reference_chain_BSA_per_kDa_reference_chain_A2": 48.0,
                "lattice_burial_fraction_percent": 9.5,
                "lattice_contact_residue_fraction_percent": 12.0,
            },
            {
                "structure_stem": "high",
                "packing_density_percent": 14.0,
                "matthews_a3_per_Da": 8.0,
                "estimated_solvent_percent": 84.0,
                "reference_chain_BSA_per_kDa_reference_chain_A2": 60.0,
                "lattice_burial_fraction_percent": 12.0,
                "lattice_contact_residue_fraction_percent": 16.0,
            },
        ]

    def test_loo_ghost_never_coincides_with_panel(self) -> None:
        from metrics import lattice_metrics_radar as lmr  # noqa: PLC0415

        rows = self._three_rows()
        scores, _raw, kept, _lim = lmr.compute_scores(rows, [r["structure_stem"] for r in rows])
        ghost_by_stem, label = lmr._ghost_values(scores, kept, None, ghost_mode="loo-median")
        self.assertIn("leave-one-out", label)
        for stem in kept:
            # Odd n: the shared cohort median equals mid's scores, but the LOO
            # ghost for mid is the mean of low and high, which differs.
            self.assertNotEqual(ghost_by_stem[stem], scores[stem])
        pack_idx = 0
        others = [scores[s][pack_idx] for s in kept if s != "mid"]
        self.assertAlmostEqual(ghost_by_stem["mid"][pack_idx], sum(others) / 2.0)

    def test_cohort_median_ghost_coincides_for_odd_n(self) -> None:
        from metrics import lattice_metrics_radar as lmr  # noqa: PLC0415

        rows = self._three_rows()
        scores, _raw, kept, _lim = lmr.compute_scores(rows, [r["structure_stem"] for r in rows])
        ghost_by_stem, label = lmr._ghost_values(scores, kept, None, ghost_mode="cohort-median")
        self.assertEqual(label, "cohort median")
        # Shared ghost equals mid's scores on every axis (odd cohort, mid is median).
        self.assertEqual(ghost_by_stem["mid"], scores["mid"])
        self.assertEqual(ghost_by_stem["low"], scores["mid"])

    def test_deviation_sd_centres_and_inverts(self) -> None:
        from metrics import lattice_metrics_radar as lmr  # noqa: PLC0415

        rows = self._three_rows()
        _scores, raw, kept, _lim = lmr.compute_scores(rows, [r["structure_stem"] for r in rows])
        dev, meta = lmr.compute_deviations(kept, raw, deviation_scale="sd")
        self.assertEqual(meta["deviation_scale"], "sd")
        # Packing (non-inverted): mean of z-scores ≈ 0, high > mid > low.
        pack = [dev[s][0] for s in ("low", "mid", "high")]
        self.assertAlmostEqual(sum(pack) / 3.0, 0.0, places=6)
        self.assertLess(pack[0], pack[1])
        self.assertLess(pack[1], pack[2])
        # Matthews (inverted): lower raw VM → higher deviation (tighter).
        matt = [dev[s][1] for s in ("low", "mid", "high")]
        self.assertLess(matt[0], matt[1])
        self.assertLess(matt[1], matt[2])
        # Sample SD of the z-scores themselves is 1.
        mu = sum(pack) / 3.0
        sd = math.sqrt(sum((z - mu) ** 2 for z in pack) / 2.0)
        self.assertAlmostEqual(sd, 1.0, places=6)

    def test_deviation_range_maps_onto_unit_interval(self) -> None:
        from metrics import lattice_metrics_radar as lmr  # noqa: PLC0415

        rows = self._three_rows()
        _scores, raw, kept, _lim = lmr.compute_scores(rows, [r["structure_stem"] for r in rows])
        dev, _meta = lmr.compute_deviations(kept, raw, deviation_scale="range")
        pack = [dev[s][0] for s in kept]
        self.assertAlmostEqual(min(pack), -1.0)
        self.assertAlmostEqual(max(pack), 1.0)

    def test_deviation_end_to_end_outputs_and_tag(self) -> None:
        from metrics import lattice_metrics_radar as lmr  # noqa: PLC0415

        with tempfile.TemporaryDirectory() as td:
            csv_path = os.path.join(td, "combined_lattice_vs_ec.csv")
            _write_combined_csv(csv_path, self._three_rows())
            out = os.path.join(td, "radar")
            written = lmr.generate_radar(
                input_csv=csv_path,
                output_dir=out,
                display_mode="deviation",
                deviation_scale="sd",
                ghost_mode="loo-median",
                out_format="png",
            )
            tagged = out + "_bbox_dev-sd"
            self.assertTrue(os.path.isdir(tagged))
            scores = os.path.join(tagged, "radar_scores_bbox_dev-sd.csv")
            scale = os.path.join(tagged, "radar_scale_bbox_dev-sd.json")
            grid = os.path.join(tagged, "radar_grid_bbox_dev-sd.png")
            heat = os.path.join(tagged, "radar_scores_heatmap_bbox_dev-sd.png")
            for path in (scores, scale, grid, heat):
                self.assertTrue(os.path.isfile(path), path)
                self.assertIn(path, written)

            with open(scores, newline="", encoding="utf-8") as f:
                header = next(csv.reader(f))
            self.assertIn("packing_density_percent__deviation", header)
            self.assertNotIn("packing_density_percent__score", header)

            import json

            with open(scale, encoding="utf-8") as f:
                meta = json.load(f)
            self.assertEqual(meta["display_mode"], "deviation")
            self.assertEqual(meta["deviation_scale"], "sd")
            self.assertEqual(meta["ghost_mode"], "loo-median")
            self.assertIn("leave-one-out", meta["ghost_label"])
            self.assertEqual(meta["plot_range"][0], -meta["plot_range"][1])
            self.assertGreater(meta["deviation_limit"], 0.0)

    def test_cli_accepts_display_and_ghost(self) -> None:
        from metrics import lattice_metrics_radar as lmr  # noqa: PLC0415

        with tempfile.TemporaryDirectory() as td:
            csv_path = os.path.join(td, "combined_lattice_vs_ec.csv")
            _write_combined_csv(csv_path, self._three_rows())
            out = os.path.join(td, "radar")
            buf = io.StringIO()
            with contextlib.redirect_stdout(buf), contextlib.redirect_stderr(io.StringIO()):
                rc = lmr.main(
                    [
                        "--input",
                        csv_path,
                        "-d",
                        out,
                        "--display",
                        "deviation",
                        "--deviation-scale",
                        "mad",
                        "--ghost",
                        "cohort-median",
                        "--no-scale-tag",
                        "--format",
                        "png",
                    ]
                )
            self.assertEqual(rc, 0)
            self.assertTrue(os.path.isfile(os.path.join(out, "radar_grid.png")))
            import json

            with open(os.path.join(out, "radar_scale.json"), encoding="utf-8") as f:
                meta = json.load(f)
            self.assertEqual(meta["display_mode"], "deviation")
            self.assertEqual(meta["deviation_scale"], "mad")
            self.assertEqual(meta["ghost_mode"], "cohort-median")

    def test_output_run_tag_composition(self) -> None:
        from metrics import lattice_metrics_radar as lmr  # noqa: PLC0415

        self.assertEqual(lmr.output_run_tag("bbox", "bbox"), "bbox")
        self.assertEqual(
            lmr.output_run_tag("bbox", "bbox", display_mode="deviation", deviation_scale="sd"),
            "bbox_dev-sd",
        )
        self.assertEqual(
            lmr.output_run_tag(
                "slayer-compact",
                "empirical",
                display_mode="deviation",
                deviation_scale="mad",
                ghost_mode="cohort-mean",
            ),
            "vol-slayer-compact_occ-empirical_dev-mad_ghost-cohort-mean",
        )
        interleaved = lmr.resolve_metric_order("interleaved")
        self.assertEqual(
            lmr.output_run_tag("bbox", "bbox", metrics=interleaved),
            "bbox_order-interleaved",
        )
        self.assertEqual(
            lmr.output_run_tag("bbox", "bbox", edge_delta=True),
            "bbox_edges-delta",
        )


class EdgeDeltaTests(unittest.TestCase):
    def test_edge_delta_styles_shared_span_mapping(self) -> None:
        from metrics import lattice_metrics_radar as lmr  # noqa: PLC0415

        # Three spokes on a 0–10 range: gaps 0, 10, 10 → fractions 0, 1, 1.
        values = [0.0, 0.0, 10.0]
        styles = lmr.edge_delta_styles(values, (0.0, 10.0))
        self.assertEqual(len(styles), 3)
        self.assertAlmostEqual(styles[0][0], lmr.EDGE_DELTA_LW_MIN)  # |0-0|
        self.assertAlmostEqual(styles[1][0], lmr.EDGE_DELTA_LW_MAX)  # |0-10|
        self.assertAlmostEqual(styles[2][0], lmr.EDGE_DELTA_LW_MAX)  # |10-0|
        self.assertAlmostEqual(styles[0][1], lmr.EDGE_DELTA_ALPHA_MIN)
        self.assertAlmostEqual(styles[1][1], lmr.EDGE_DELTA_ALPHA_MAX)
        # Mid-span gap is boosted by power < 1 (clearer than linear).
        mid = lmr.edge_delta_styles([0.0, 5.0], (0.0, 10.0))
        t = 0.5**lmr.EDGE_DELTA_POWER
        mid_lw = lmr.EDGE_DELTA_LW_MIN + (lmr.EDGE_DELTA_LW_MAX - lmr.EDGE_DELTA_LW_MIN) * t
        self.assertAlmostEqual(mid[0][0], mid_lw)
        self.assertGreater(mid_lw, lmr.EDGE_DELTA_LW_MIN + 0.5 * (lmr.EDGE_DELTA_LW_MAX - lmr.EDGE_DELTA_LW_MIN))

    def test_edge_delta_cli_writes_tagged_grid(self) -> None:
        import json  # noqa: PLC0415

        from metrics import lattice_metrics_radar as lmr  # noqa: PLC0415

        with tempfile.TemporaryDirectory() as td:
            csv_path = os.path.join(td, "combined_lattice_vs_ec.csv")
            _write_combined_csv(csv_path, _sample_rows())
            out = os.path.join(td, "radar")
            with contextlib.redirect_stdout(io.StringIO()), contextlib.redirect_stderr(
                io.StringIO()
            ):
                rc = lmr.main(
                    [
                        "--input",
                        csv_path,
                        "--output-dir",
                        out,
                        "--edge-delta",
                        "--no-heatmap",
                    ]
                )
            self.assertEqual(rc, 0)
            tagged = out + "_bbox_edges-delta"
            self.assertTrue(os.path.isfile(os.path.join(tagged, "radar_grid_bbox_edges-delta.png")))
            with open(
                os.path.join(tagged, "radar_scale_bbox_edges-delta.json"),
                encoding="utf-8",
            ) as f:
                meta = json.load(f)
            self.assertTrue(meta["edge_delta"])
            self.assertEqual(meta["edge_delta_scale"]["mode"], "plot_span")


class MetricOrderTests(unittest.TestCase):
    def test_presets_and_custom_tokens(self) -> None:
        from metrics import lattice_metrics_radar as lmr  # noqa: PLC0415

        default = lmr.resolve_metric_order("default")
        self.assertEqual([m.column for m in default], list(lmr.METRIC_ORDER_PRESETS["default"]))
        self.assertEqual(
            [m.column for m in default],
            [
                "lattice_burial_fraction_percent",
                "matthews_a3_per_Da",
                "lattice_contact_residue_fraction_percent",
                "packing_density_percent",
                "reference_chain_BSA_per_kDa_reference_chain_A2",
                "estimated_solvent_percent",
            ],
        )
        self.assertEqual(lmr.metric_order_tag(default), "")

        interleaved = lmr.resolve_metric_order("interleaved")
        cols = [m.column for m in interleaved]
        self.assertEqual(cols, list(lmr.METRIC_ORDER_PRESETS["interleaved"]))
        self.assertEqual(cols[0], "packing_density_percent")
        self.assertEqual(cols[1], "lattice_contact_residue_fraction_percent")
        self.assertEqual(cols[2], "matthews_a3_per_Da")
        self.assertEqual(cols[3], "lattice_burial_fraction_percent")
        self.assertEqual(cols[4], "estimated_solvent_percent")
        self.assertEqual(cols[5], "reference_chain_BSA_per_kDa_reference_chain_A2")
        self.assertEqual(lmr.metric_order_tag(interleaved), "order-interleaved")

        custom = lmr.resolve_metric_order("bsa,loi,contact,packing,matthews,solvent")
        self.assertEqual(
            [m.column for m in custom],
            [
                "reference_chain_BSA_per_kDa_reference_chain_A2",
                "lattice_burial_fraction_percent",
                "lattice_contact_residue_fraction_percent",
                "packing_density_percent",
                "matthews_a3_per_Da",
                "estimated_solvent_percent",
            ],
        )
        self.assertEqual(lmr.metric_order_tag(custom), "order-custom")

    def test_loi_spoke_at_top_for_default_and_interleaved(self) -> None:
        from metrics import lattice_metrics_radar as lmr  # noqa: PLC0415

        for order in ("default", "interleaved"):
            metrics = lmr.resolve_metric_order(order)
            angles = lmr.spoke_angles(metrics)
            loi_i = next(
                i
                for i, m in enumerate(metrics)
                if m.column == "lattice_burial_fraction_percent"
            )
            self.assertAlmostEqual(angles[loi_i], 0.0)
            step = 2.0 * math.pi / len(metrics)
            next_i = (loi_i + 1) % len(metrics)
            self.assertAlmostEqual(angles[next_i], step)

    def test_rejects_incomplete_or_unknown(self) -> None:
        from metrics import lattice_metrics_radar as lmr  # noqa: PLC0415

        with self.assertRaises(ValueError):
            lmr.resolve_metric_order("packing,matthews,solvent")
        with self.assertRaises(ValueError):
            lmr.resolve_metric_order("packing,packing,matthews,solvent,bsa,loi")
        with self.assertRaises(ValueError):
            lmr.resolve_metric_order("not-a-metric,matthews,solvent,bsa,loi,contact")

    def test_compute_scores_respects_order(self) -> None:
        from metrics import lattice_metrics_radar as lmr  # noqa: PLC0415

        rows = [
            {
                "structure_stem": "a",
                "packing_density_percent": "30",
                "matthews_a3_per_Da": "2",
                "estimated_solvent_percent": "40",
                "reference_chain_BSA_per_kDa_reference_chain_A2": "150",
                "lattice_burial_fraction_percent": "20",
                "lattice_contact_residue_fraction_percent": "25",
            }
        ]
        metrics = lmr.resolve_metric_order("interleaved")
        scores, raw, kept, _ = lmr.compute_scores(
            rows,
            ["a"],
            volume_limits={
                "packing_density_percent": (20.0, 40.0),
                "matthews_a3_per_Da": (1.5, 4.0),
                "estimated_solvent_percent": (25.0, 80.0),
            },
            occlusion_limits={
                "reference_chain_BSA_per_kDa_reference_chain_A2": (100.0, 200.0),
                "lattice_burial_fraction_percent": (10.0, 30.0),
                "lattice_contact_residue_fraction_percent": (15.0, 35.0),
            },
            metrics=metrics,
        )
        self.assertEqual(kept, ["a"])
        self.assertEqual(len(raw["a"]), 6)
        self.assertAlmostEqual(raw["a"][0], 30.0)  # packing
        self.assertAlmostEqual(raw["a"][1], 25.0)  # contact
        self.assertAlmostEqual(raw["a"][2], 2.0)  # Matthews
        self.assertAlmostEqual(raw["a"][5], 150.0)  # BSA


class DualLayoutRangeTests(unittest.TestCase):
    """Both CSV layouts must load and score inside default bbox bands."""

    def _rounded_summary_rows(self) -> list[list[object]]:
        # Mimic metrics_summary_table.csv display rounding (integers / 1 d.p.).
        return [
            ["", "SlpA2", "SlpA7", "SlpA7b", "SlpA11", "SlpA4dD2"],
            ["Packing density (%)", 12, 14, 12, 12, 13],
            ["Matthews-like coefficient (Å3/Da)", 9, 8, 9, 9, 8],
            ["Estimated solvent (%)", 86, 84, 86, 86, 85],
            ["BSAmolA/kDa (Å²)", 56.5, 56.9, 43.1, 45.2, 48.0],
            ["Lattice occlusion index (LOImolA) (%)", 12, 12, 9, 10, 10],
            ["Mol. A lattice contact residues (%)", 15, 15, 12, 13, 12],
        ]

    def _precise_combined_rows(self) -> list[dict[str, object]]:
        return [
            {
                "structure_stem": "s2_opt2427_15mA",
                "packing_density_percent": 11.11,
                "matthews_a3_per_Da": 9.57,
                "estimated_solvent_percent": 87.10,
                "reference_chain_BSA_per_kDa_reference_chain_A2": 56.47,
                "lattice_burial_fraction_percent": 11.8,
                "lattice_contact_residue_fraction_percent": 15.2,
            },
            {
                "structure_stem": "s7_cd630_15mA",
                "packing_density_percent": 12.36,
                "matthews_a3_per_Da": 8.62,
                "estimated_solvent_percent": 85.67,
                "reference_chain_BSA_per_kDa_reference_chain_A2": 56.94,
                "lattice_burial_fraction_percent": 12.1,
                "lattice_contact_residue_fraction_percent": 14.8,
            },
            {
                "structure_stem": "s7b_r7404_15mA",
                "packing_density_percent": 10.78,
                "matthews_a3_per_Da": 9.87,
                "estimated_solvent_percent": 87.49,
                "reference_chain_BSA_per_kDa_reference_chain_A2": 43.05,
                "lattice_burial_fraction_percent": 8.9,
                "lattice_contact_residue_fraction_percent": 11.6,
            },
            {
                "structure_stem": "s11_ox247o_15mA",
                "packing_density_percent": 10.75,
                "matthews_a3_per_Da": 9.90,
                "estimated_solvent_percent": 87.52,
                "reference_chain_BSA_per_kDa_reference_chain_A2": 45.21,
                "lattice_burial_fraction_percent": 10.1,
                "lattice_contact_residue_fraction_percent": 13.0,
            },
            {
                "structure_stem": "s4del_r2d2_15mA",
                "packing_density_percent": 11.70,
                "matthews_a3_per_Da": 9.13,
                "estimated_solvent_percent": 86.48,
                "reference_chain_BSA_per_kDa_reference_chain_A2": 47.96,
                "lattice_burial_fraction_percent": 10.2,
                "lattice_contact_residue_fraction_percent": 12.4,
            },
        ]

    def test_rounded_summary_and_precise_combined_share_bbox_bands(self) -> None:
        import json  # noqa: PLC0415

        from metrics import lattice_metrics_radar as lmr  # noqa: PLC0415

        with tempfile.TemporaryDirectory() as td:
            summary = os.path.join(td, "metrics_summary_table.csv")
            with open(summary, "w", newline="", encoding="utf-8") as f:
                csv.writer(f).writerows(self._rounded_summary_rows())
            combined = os.path.join(td, "combined_lattice_vs_ec.csv")
            _write_combined_csv(combined, self._precise_combined_rows())

            for label, path in (("summary", summary), ("combined", combined)):
                out = os.path.join(td, f"radar_{label}")
                with contextlib.redirect_stderr(io.StringIO()):
                    written = lmr.generate_radar(
                        input_csv=path,
                        output_dir=out,
                        scale_mode="bbox",
                        occlusion_scale_mode="bbox",
                    )
                tagged = out + "_bbox"
                scale_path = os.path.join(tagged, "radar_scale_bbox.json")
                self.assertTrue(os.path.isfile(scale_path), label)
                self.assertTrue(
                    any(p.endswith("radar_grid_bbox.png") for p in written), label
                )
                with open(scale_path, encoding="utf-8") as f:
                    meta = json.load(f)
                self.assertEqual(
                    meta["used_limits"]["matthews_a3_per_Da"], [8.0, 45.0], label
                )
                self.assertEqual(
                    meta["used_limits"]["estimated_solvent_percent"], [80.0, 100.0], label
                )
                self.assertEqual(
                    meta["used_limits"]["packing_density_percent"], [2.0, 15.0], label
                )
                self.assertEqual(
                    meta["used_limits"]["lattice_burial_fraction_percent"],
                    [2.0, 20.0],
                    label,
                )
                self.assertEqual(meta["n_clipped"], 0, f"{label}: {meta.get('clipped')}")

    def test_compact_mode_accepts_rounded_summary_solvent(self) -> None:
        import json  # noqa: PLC0415

        from metrics import lattice_metrics_radar as lmr  # noqa: PLC0415

        with tempfile.TemporaryDirectory() as td:
            summary = os.path.join(td, "metrics_summary_table.csv")
            with open(summary, "w", newline="", encoding="utf-8") as f:
                csv.writer(f).writerows(self._rounded_summary_rows())
            out = os.path.join(td, "radar")
            with contextlib.redirect_stderr(io.StringIO()):
                lmr.generate_radar(
                    input_csv=summary,
                    output_dir=out,
                    scale_mode="slayer-compact",
                    occlusion_scale_mode="slayer-compact",
                )
            with open(
                os.path.join(out + "_slayer-compact", "radar_scale_slayer-compact.json"),
                encoding="utf-8",
            ) as f:
                meta = json.load(f)
            self.assertEqual(meta["used_limits"]["estimated_solvent_percent"], [80.0, 90.0])
            self.assertEqual(meta["n_clipped"], 0, meta.get("clipped"))

    def test_empirical_wide_loi_keeps_nonzero_floor(self) -> None:
        from metrics import lattice_metrics_radar as lmr  # noqa: PLC0415

        limits, _ = lmr.resolve_occlusion_limits(
            occlusion_scale_mode="empirical",
            occlusion_values={
                "lattice_burial_fraction_percent": [2.3, 16.0],
                "lattice_contact_residue_fraction_percent": [2.7, 20.2],
                "reference_chain_BSA_per_kDa_reference_chain_A2": [9.9, 76.8],
            },
        )
        self.assertEqual(limits["lattice_burial_fraction_percent"], (2.0, 20.0))
        self.assertEqual(limits["lattice_contact_residue_fraction_percent"], (2.0, 25.0))
        self.assertEqual(
            limits["reference_chain_BSA_per_kDa_reference_chain_A2"], (5.0, 80.0)
        )


class ErrorTests(unittest.TestCase):
    def test_missing_metric_column_errors(self) -> None:
        from metrics import lattice_metrics_radar as lmr  # noqa: PLC0415

        with tempfile.TemporaryDirectory() as td:
            csv_path = os.path.join(td, "combined_lattice_vs_ec.csv")
            _write_combined_csv(csv_path, _sample_rows(), drop_columns=("lattice_burial_fraction_percent",))
            with self.assertRaises(ValueError):
                lmr.load_rows(csv_path)

    def test_cli_missing_column_exits_nonzero(self) -> None:
        from metrics import lattice_metrics_radar as lmr  # noqa: PLC0415

        with tempfile.TemporaryDirectory() as td:
            csv_path = os.path.join(td, "combined_lattice_vs_ec.csv")
            _write_combined_csv(csv_path, _sample_rows(), drop_columns=("matthews_a3_per_Da",))
            with contextlib.redirect_stderr(io.StringIO()):
                with self.assertRaises(SystemExit):
                    lmr.main(["--input", csv_path, "--output-dir", os.path.join(td, "radar")])


if __name__ == "__main__":
    unittest.main()
