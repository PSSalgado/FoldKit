"""Tests for utils.ssm_log_parse helpers."""

from utils.ssm_log_parse import ResidueKey, _normalize_icode, parse_superposition_labels


def test_residue_key_as_tuple_default_icode():
    key = ResidueKey(chain="A", resnum=42)
    assert key.as_tuple() == ("A", 42, " ")


def test_residue_key_as_tuple_explicit_icode():
    key = ResidueKey(chain="B", resnum=10, icode="A")
    assert key.as_tuple() == ("B", 10, "A")


def test_normalize_icode_empty():
    assert _normalize_icode("") == " "
    assert _normalize_icode(None) == " "
    assert _normalize_icode("   ") == " "


def test_normalize_icode_first_char():
    assert _normalize_icode("A") == "A"
    assert _normalize_icode("AB") == "A"


def test_parse_superposition_labels():
    assert parse_superposition_labels("Superposing model_01 onto ref_01") == (
        "model_01",
        "ref_01",
    )
    assert parse_superposition_labels("  Superposing model onto ref") == (
        "model",
        "ref",
    )
    assert parse_superposition_labels("not a superposition line") is None
