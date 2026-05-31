#!/usr/bin/env python3
"""
Write a new PDB containing all lines from the input **except** ``HETATM`` records.

Non-``HETATM`` lines (``ATOM``, ``TER``, ``CRYST1``, headers, ``END``, etc.) are
copied unchanged. Use this to build a protein-only (no heteroatom) coordinate
set for SASA / lattice comparisons.

Example::

  python file_management/strip_hetatm_from_pdb.py /path/to/input.pdb /path/to/output_nohet.pdb
"""

from __future__ import annotations

import argparse
import sys


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("input_pdb", help="Input PDB path")
    ap.add_argument("output_pdb", help="Output PDB path")
    args = ap.parse_args()
    n_in = n_out = n_skip = 0
    with open(args.input_pdb, encoding="utf-8", errors="replace") as fin, open(
        args.output_pdb, "w", encoding="utf-8", newline="\n"
    ) as fout:
        for line in fin:
            n_in += 1
            if line.startswith("HETATM"):
                n_skip += 1
                continue
            fout.write(line)
            n_out += 1
    print(
        f"Wrote {args.output_pdb}  (lines in: {n_in}, written: {n_out}, skipped HETATM: {n_skip})",
        file=sys.stderr,
    )


if __name__ == "__main__":
    main()
