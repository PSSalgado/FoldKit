"""Tests for file_management.strip_hetatm_from_pdb."""

import sys
from pathlib import Path

import pytest

from file_management import strip_hetatm_from_pdb


def test_strip_hetatm_drops_hetatm_keeps_atom(tmp_path: Path, monkeypatch: pytest.MonkeyPatch):
    inp = tmp_path / "in.pdb"
    out = tmp_path / "out.pdb"
    inp.write_text(
        "HEADER    TEST\n"
        "ATOM      1  CA  ALA A   1      11.104  12.211  13.311  1.00 20.00           C  \n"
        "HETATM    2  O   HOH A   2      21.104  22.211  23.311  1.00 20.00           O  \n"
        "TER\n"
        "END\n",
        encoding="utf-8",
    )
    monkeypatch.setattr(
        sys,
        "argv",
        ["strip_hetatm_from_pdb.py", str(inp), str(out)],
    )
    strip_hetatm_from_pdb.main()
    text = out.read_text(encoding="utf-8")
    assert "ATOM" in text
    assert "HEADER" in text
    assert "TER" in text
    assert "END" in text
    assert "HETATM" not in text
