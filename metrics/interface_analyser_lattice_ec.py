#!/usr/bin/env python3
"""
Interface analyser (lattice, electrostatic complementarity)
==========================================================

Entry point for interface analysis in symmetry-expanded / multi-copy lattice
models, reporting electrostatic complementarity (EC) using the McCoy, Epa &
Colman (1997) definition.

Use ``--phase full`` (default) for one-shot SASA/BSA plus EC. Use ``--phase sasa``
then ``--phase ec`` with a JSON sidecar to split wall time (same physics; see README).
"""

import argparse
import sys

from metrics.interface_analyser_base import (
    InterfaceAnalyserEC,
    _foldkit_interface_workers,
    _run_analysis,
    collect_structure_paths,
)
from metrics.interface_lattice_ec_sidecar import (
    build_cli_signature,
    read_sidecar,
    results_from_sidecar,
    validate_cli_signature,
    validate_input_fingerprint,
    validate_interfaces_for_ec,
    write_sidecar,
)

from utils.cli_log import add_log_args, setup_log_from_args


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Analyse lattice interfaces for a reference chain and report electrostatic complementarity (EC).",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples (from repository root):
  # One-shot (default)
  python metrics/interface_analyser_lattice_ec.py expanded_assembly.pdb --reference-chain A -o lattice_ec.txt

  # Two-phase: SASA/BSA first (writes sidecar), EC later (reads sidecar + same PDB)
  python metrics/interface_analyser_lattice_ec.py asm.pdb --reference-chain A --phase sasa \\
    --sidecar-out asm_lattice_ec_sidecar.json -o lattice_sasa.txt
  python metrics/interface_analyser_lattice_ec.py asm.pdb --reference-chain A --phase ec \\
    --sidecar asm_lattice_ec_sidecar.json -o lattice_ec.txt
""",
    )
    parser.add_argument(
        "input",
        nargs="+",
        help="PDB/CIF file(s), directory, or glob pattern (e.g. '*.pdb').",
    )
    parser.add_argument(
        "--phase",
        choices=("full", "sasa", "ec"),
        default="full",
        help="full: SASA/BSA then EC (default). sasa: SASA/BSA only; requires --sidecar-out. ec: EC only from --sidecar.",
    )
    parser.add_argument(
        "--sidecar-out",
        metavar="PATH",
        help="Phase sasa: write JSON sidecar for a later --phase ec run.",
    )
    parser.add_argument(
        "--sidecar",
        metavar="PATH",
        help="Phase ec: JSON sidecar produced by --phase sasa with matching flags and PDB bytes.",
    )
    parser.add_argument(
        "--reference-chain",
        required=True,
        metavar="ID",
        dest="reference_chain_id",
        help="Reference chain ID (focal copy) for lattice-wide metrics.",
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
    parser.add_argument(
        "--skip-accessibility-sasa",
        action="store_true",
        help=(
            "Skip full-model per-residue SASA used only for contact-residue accessibility "
            "in the report (much faster on large assemblies; accessibility averages stay at "
            "defaults when SASA is missing)."
        ),
    )
    parser.add_argument(
        "--ec-max-contact-points",
        type=int,
        default=None,
        metavar="N",
        help=(
            "Cap contact heavy atoms used as EC surface samples (both chains combined). "
            "Deterministic subsampling for huge interfaces."
        ),
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

    reference_chain_id = str(args.reference_chain_id).strip() or None
    if not reference_chain_id:
        print("Error: --reference-chain must be a non-empty chain ID.", file=sys.stderr)
        raise SystemExit(1)

    phase = str(args.phase or "full").strip().lower()
    if phase == "full":
        if getattr(args, "sidecar_out", None):
            print("Error: --sidecar-out is only valid with --phase sasa.", file=sys.stderr)
            raise SystemExit(1)
        if getattr(args, "sidecar", None):
            print("Error: --sidecar is only valid with --phase ec.", file=sys.stderr)
            raise SystemExit(1)
    elif phase == "sasa":
        if not getattr(args, "sidecar_out", None):
            print("Error: --phase sasa requires --sidecar-out PATH.", file=sys.stderr)
            raise SystemExit(1)
        if getattr(args, "sidecar", None):
            print("Error: do not pass --sidecar with --phase sasa.", file=sys.stderr)
            raise SystemExit(1)
    else:
        if not getattr(args, "sidecar", None):
            print("Error: --phase ec requires --sidecar PATH.", file=sys.stderr)
            raise SystemExit(1)
        if getattr(args, "sidecar_out", None):
            print("Error: do not pass --sidecar-out with --phase ec.", file=sys.stderr)
            raise SystemExit(1)

    paths = collect_structure_paths(args.input)
    if not paths:
        print("No structure files found.", file=sys.stderr)
        raise SystemExit(1)
    if len(paths) != 1 and not args.output:
        print(
            "Error: lattice mode expects one input structure; specify --output (-o) when passing multiple.",
            file=sys.stderr,
        )
        raise SystemExit(1)
    if phase in ("sasa", "ec") and len(paths) != 1:
        print(
            "Error: phased lattice EC (--phase sasa or ec) expects exactly one structure path.",
            file=sys.stderr,
        )
        raise SystemExit(1)

    ec_cap = getattr(args, "ec_max_contact_points", None)
    analyser = InterfaceAnalyserEC(
        ec_mode="full",
        skip_accessibility_sasa=bool(getattr(args, "skip_accessibility_sasa", False)),
        ec_max_contact_points=int(ec_cap) if ec_cap is not None else None,
    )

    cli_sig = build_cli_signature(
        reference_chain_id=reference_chain_id,
        focus_chains=focus_chains or [],
        contact_distance=float(analyser.contact_distance),
        skip_accessibility_sasa=bool(analyser.skip_accessibility_sasa),
        ec_max_contact_points=getattr(analyser, "ec_max_contact_points", None),
    )

    analyse_hook = None
    after_results = None
    sidecar_obj = None

    if phase == "sasa":

        def analyse_hook(analyser_inst, pdb_file):
            return analyser_inst.analyze_interfaces_sasa_only(
                pdb_file,
                focus_chains=focus_chains,
                reference_chain_id=reference_chain_id,
            )

        sidecar_out_path = str(args.sidecar_out)

        def after_results(pdb_file, results):
            write_sidecar(
                sidecar_out_path,
                pdb_path=pdb_file,
                results=results,
                cli_signature=cli_sig,
            )

    elif phase == "ec":
        sidecar_obj = read_sidecar(str(args.sidecar))
        validate_cli_signature(sidecar_obj["cli_signature"], cli_sig)
        validate_interfaces_for_ec(sidecar_obj["interfaces"])

        def analyse_hook(analyser_inst, pdb_file):
            validate_input_fingerprint(sidecar_obj, pdb_file)
            res = results_from_sidecar(sidecar_obj)
            return analyser_inst.apply_ec_phase(pdb_file, res, reference_chain_id)

    out = sys.stdout
    if args.output:
        out = open(args.output, "w")
    try:
        out_desc = args.output if args.output else "stdout"
        iw = _foldkit_interface_workers()
        iw_note = f" parallel_workers={iw}" if iw > 1 else ""
        phase_note = f" phase={phase}"
        print(
            f"[FoldKit] interface (lattice, EC): inputs={len(paths)} focus_chains={focus_chains or 'ALL'} "
            f"reference_chain={reference_chain_id} output={out_desc}{phase_note}{iw_note}",
            file=sys.stderr,
            flush=True,
        )
        if phase == "sasa":
            print(f"[FoldKit] sidecar-out={args.sidecar_out}", file=sys.stderr, flush=True)
        elif phase == "ec":
            print(f"[FoldKit] sidecar={args.sidecar}", file=sys.stderr, flush=True)
        if log_setup is not None:
            log_setup.task(f"Interface analysis (lattice, EC): {len(paths)} structure(s)")
            log_setup.task(f"Reference chain: {reference_chain_id}")
            log_setup.task(f"Focus chains: {','.join(focus_chains) if focus_chains else 'ALL'}")
        _run_analysis(
            analyser,
            paths,
            out,
            focus_chains=focus_chains,
            reference_chain_id=reference_chain_id,
            summary_log=log_setup,
            analyse_hook=analyse_hook,
            after_results=after_results,
        )
    finally:
        if args.output and out is not sys.stdout:
            out.close()


if __name__ == "__main__":
    main()
