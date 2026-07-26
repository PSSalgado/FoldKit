"""Tests for foldkit heatmap per-metric CLI helpers used by interface-matrix scripts."""

from __future__ import annotations

import pytest

from metrics.interface_matrix_heatmap_cli import METRIC_CHOICES
from utils.foldkit_heatmap import metric_str_map_from_cli


def test_metric_str_map_cmap_by_metric_parses_repeatable() -> None:
    m = metric_str_map_from_cli(
        METRIC_CHOICES,
        ["bsa=plasma", "ec=RdBu_r", "contacts=Greys"],
        allowed=None,
    )
    assert m["bsa"] == "plasma"
    assert m["ec"] == "RdBu_r"
    assert m["contacts"] == "Greys"


def test_metric_str_map_unknown_metric_raises() -> None:
    with pytest.raises(ValueError, match="Unknown heatmap metric"):
        metric_str_map_from_cli(METRIC_CHOICES, ["notametric=viridis"], allowed=None)


def test_metric_str_map_invalid_token_raises() -> None:
    with pytest.raises(ValueError, match="expected METRIC=VALUE"):
        metric_str_map_from_cli(METRIC_CHOICES, ["bsa"], allowed=None)


def test_metric_str_map_empty_value_raises() -> None:
    with pytest.raises(ValueError, match="Empty value"):
        metric_str_map_from_cli(METRIC_CHOICES, ["bsa="], allowed=None)


def test_metric_str_map_diverging_center_normalized() -> None:
    m = metric_str_map_from_cli(
        METRIC_CHOICES,
        ["ec=MEDIAN", "bsa=NONE"],
        allowed=frozenset({"none", "median"}),
        normalize_value=str.lower,
    )
    assert m["ec"] == "median"
    assert m["bsa"] == "none"
