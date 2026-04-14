#!/usr/bin/env python3
"""
Interface analyser (lattice / multi-copy)
=======================================

Pairwise chain–chain interface analysis plus detailed lattice-wide metrics for
one reference chain in a multi-copy (symmetry-expanded / supercell) model.

This is the "lattice" entrypoint: it requires --reference-chain and will print
the extended lattice-wide charge complementarity metrics.

The underlying interface calculations (pairwise contacts, BSA, etc.) are the
same as in simple mode; the difference is the additional lattice summary block.
"""

import argparse
import sys

from interface_analyser import InterfaceAnalyser, collect_structure_paths, _run_analysis


def main():
    parser = argparse.ArgumentParser(
        description="Analyse interfaces in multi-copy models and compute lattice-wide metrics for a reference chain.",
    )
    parser.add_argument(
        'input',
        nargs='+',
        help='PDB/CIF file(s), directory, or glob pattern (e.g. *.pdb).',
    )
    parser.add_argument(
        '--output', '-o',
        metavar='FILE',
        help='Output file (single run). If omitted, write to stdout.',
    )
    parser.add_argument(
        '--chains',
        metavar='IDS',
        help="Focus analysis on specific chain IDs (comma-separated). Only chain pairs where at least one chain is in this list are analysed.",
    )
    parser.add_argument(
        '--reference-chain',
        metavar='ID',
        dest='reference_chain_id',
        required=True,
        help='Reference chain ID (focal copy) for lattice-wide metrics.',
    )
    args = parser.parse_args()

    focus_chains = None
    if getattr(args, 'chains', None):
        focus_chains = [c.strip() for c in str(args.chains).split(',') if c.strip()]

    reference_chain_id = str(args.reference_chain_id).strip() or None
    if not reference_chain_id:
        print("Error: --reference-chain must be a non-empty chain ID.", file=sys.stderr)
        sys.exit(1)

    paths = collect_structure_paths(args.input)
    if not paths:
        print("No structure files found.", file=sys.stderr)
        sys.exit(1)
    if len(paths) != 1 and not args.output:
        print("Error: lattice mode expects one input structure unless you provide -o.", file=sys.stderr)
        print("Tip: use a shell loop or run per-structure in batch drivers.", file=sys.stderr)
        sys.exit(1)

    analyser = InterfaceAnalyser()
    out = sys.stdout
    if args.output:
        out = open(args.output, 'w')
    try:
        _run_analysis(
            analyser,
            paths,
            out,
            focus_chains=focus_chains,
            reference_chain_id=reference_chain_id,
        )
    finally:
        if out is not sys.stdout:
            out.close()


if __name__ == "__main__":
    main()

