#!/usr/bin/env python3
"""
Interface analyser (lattice, charge metrics)
============================================

Entry point for charge/contact-based interface analysis in symmetry-expanded /
multi-copy lattice models, reporting reference-chain lattice charge summaries.

This script does not compute electrostatic complementarity (EC).
"""

import argparse
import sys

from interface_analyser_base import collect_structure_paths, _run_analysis
from interface_analyser_base import InterfaceAnalyser
from cli_log import add_log_args, setup_log_from_args


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Analyse lattice interfaces for a reference chain and report charge-based metrics.",
    )
    parser.add_argument(
        "input",
        nargs="+",
        help="PDB/CIF file(s), directory, or glob pattern (e.g. '*.pdb').",
    )
    parser.add_argument(
        "--reference-chain",
        required=True,
        metavar="ID",
        dest="reference_chain_id",
        help="Reference chain ID (focal copy) for lattice-wide charge metrics.",
    )
    parser.add_argument(
        "--output",
        "-o",
        metavar="FILE",
        help="Output file. If omitted, write to stdout.",
    )
    parser.add_argument(
        "--chains",
        metavar="IDS",
        help="Focus analysis on specific chain IDs (comma-separated).",
    )
    add_log_args(parser)
    args = parser.parse_args()
    setup_log_from_args(
        args,
        script_path=__file__,
        inputs=list(getattr(args, "input", []) or []),
        pattern=None,
    )

    focus_chains = None
    if getattr(args, "chains", None):
        focus_chains = [c.strip() for c in str(args.chains).split(",") if c.strip()]

    reference_chain_id = str(args.reference_chain_id).strip() or None
    if not reference_chain_id:
        print("Error: --reference-chain must be a non-empty chain ID.", file=sys.stderr)
        raise SystemExit(1)

    paths = collect_structure_paths(args.input)
    if not paths:
        print("No structure files found.", file=sys.stderr)
        raise SystemExit(1)
    if len(paths) != 1 and not args.output:
        print("Error: lattice mode expects one input structure unless you provide -o.", file=sys.stderr)
        raise SystemExit(1)

    analyser = InterfaceAnalyser()
    out = sys.stdout
    if args.output:
        out = open(args.output, "w")
    try:
        out_desc = args.output if args.output else "stdout"
        print(
            f"[FoldKit] interface (lattice, charge): inputs={len(paths)} focus_chains={focus_chains or 'ALL'} "
            f"reference_chain={reference_chain_id} output={out_desc}",
            file=sys.stderr,
            flush=True,
        )
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

