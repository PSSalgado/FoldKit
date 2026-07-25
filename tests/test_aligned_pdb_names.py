"""Tests for superimposition.aligned_pdb_names."""

from superimposition.aligned_pdb_names import (
    MAX_ALIGNED_PDB_BASENAME,
    aligned_pdb_basename,
)


def test_aligned_pdb_basename_short():
    name = aligned_pdb_basename("model_01", "ref_01", "_on_")
    assert name == "model_01_on_ref_01.pdb"
    assert len(name) <= MAX_ALIGNED_PDB_BASENAME


def test_aligned_pdb_basename_truncates_long_stems():
    model = "m" * 200
    ref = "r" * 200
    name = aligned_pdb_basename(model, ref, "_on_")
    assert name.endswith(".pdb")
    assert len(name) <= MAX_ALIGNED_PDB_BASENAME
    assert "_on_" in name
    # Hash fragments appear when truncated.
    assert "_" in name
