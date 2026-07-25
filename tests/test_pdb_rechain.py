"""Tests for file_management.pdb_rechain helpers."""

import pytest

from file_management.pdb_rechain import (
    _is_atom_line,
    _parse_merge_map,
    _set_chain_resseq,
)


def _atom_line(chain: str = "A", resseq: int = 1) -> str:
    # Minimal PDB ATOM line with chain at column 22 and resseq at 23–26 (1-based).
    # Columns: 1–6 record, 7–11 serial, 13–16 name, 18–20 resname, 22 chain, 23–26 resseq.
    return (
        f"ATOM      1  CA  ALA {chain}{resseq:4d}      "
        f"11.104  12.211  13.311  1.00 20.00           C  \n"
    )


def test_is_atom_line():
    assert _is_atom_line(_atom_line()) is True
    assert _is_atom_line("HETATM    1  O   HOH A   1      "
                         "11.104  12.211  13.311  1.00 20.00           O  \n") is True
    assert _is_atom_line("HEADER    TEST\n") is False
    assert _is_atom_line("TER\n") is False


def test_parse_merge_map():
    assert _parse_merge_map("B:A,D:C,F:E") == [("B", "A"), ("D", "C"), ("F", "E")]
    assert _parse_merge_map(" B:A , D:C ") == [("B", "A"), ("D", "C")]


def test_parse_merge_map_errors():
    with pytest.raises(ValueError, match="empty"):
        _parse_merge_map("")
    with pytest.raises(ValueError, match="FROM:TO"):
        _parse_merge_map("BA")
    with pytest.raises(ValueError, match="1 char"):
        _parse_merge_map("AB:C")


def test_set_chain_resseq():
    line = _atom_line("A", 1)
    out = _set_chain_resseq(line, "B", 99)
    assert out[21] == "B"
    assert out[22:26] == "  99"
    assert out.startswith("ATOM")
    assert out[27:] == line[27:]
