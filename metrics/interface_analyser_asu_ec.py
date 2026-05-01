#!/usr/bin/env python3
"""
Interface analyser (ASU, electrostatic complementarity)
======================================================

Entry point for interface analysis in the asymmetric unit (ASU) reporting
electrostatic complementarity (EC) using the McCoy, Epa & Colman (1997)
definition.
"""

import argparse
import sys

from metrics.interface_analyser_base import InterfaceAnalyserEC, _run_analysis, collect_structure_paths

from utils.cli_log import add_log_args, setup_log_from_args


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Analyse pairwise interfaces in the ASU and report electrostatic complementarity (EC).",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples (from repository root):
  python metrics/interface_analyser_asu_ec.py model_01.pdb -o ec_results.txt
  python metrics/interface_analyser_asu_ec.py model_01.pdb --chains A,B -o ec_AB.txt
""",
    )
    parser.add_argument(
        "input",
        nargs="+",
        help="PDB/CIF file(s), directory, or glob pattern (e.g. '*.pdb').",
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
    log_setup = setup_log_from_args(
        args,
        script_path=__file__,
        inputs=list(getattr(args, "input", []) or []),
        pattern=None,
    )

    focus_chains = None
    if getattr(args, "chains", None):
        focus_chains = [c.strip() for c in str(args.chains).split(",") if c.strip()]

    paths = collect_structure_paths(args.input)
    if not paths:
        print("No structure files found.", file=sys.stderr)
        raise SystemExit(1)
    if len(paths) != 1 and not args.output:
        print("Error: ASU mode expects one input structure; specify --output (-o) when passing multiple.", file=sys.stderr)
        raise SystemExit(1)

    analyser = InterfaceAnalyserEC(ec_mode="full")
    out = sys.stdout
    if args.output:
        out = open(args.output, "w")
    try:
        out_desc = args.output if args.output else "stdout"
        print(
            f"[FoldKit] interface (ASU, EC): inputs={len(paths)} focus_chains={focus_chains or 'ALL'} "
            f"output={out_desc}",
            file=sys.stderr,
            flush=True,
        )
        if log_setup is not None:
            log_setup.task(f"Interface analysis (ASU, EC): {len(paths)} structure(s)")
            log_setup.task(f"Focus chains: {','.join(focus_chains) if focus_chains else 'ALL'}")
        _run_analysis(
            analyser,
            paths,
            out,
            focus_chains=focus_chains,
            reference_chain_id=None,
            summary_log=log_setup,
        )
    finally:
        if args.output and out is not sys.stdout:
            out.close()


if __name__ == "__main__":
    main()

