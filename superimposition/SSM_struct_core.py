#!/usr/bin/env python3
"""
Find the structural core from SSM all-vs-all (or multi-pair) superpositions.

Reads Coot SSM logs and/or pre-extracted ``ssm_equivalences/*.tsv`` files,
merges pairwise structural equivalences, and writes:

  ssm_structural_msa.fasta / .aln  — MSA over all merged alignment columns
  ssm_struct_core.fasta / .aln     — core only (reference residues in all pairs)
  ssm_core_equivalences.tsv        — one row per core column (residue IDs per structure)
  ssm_structural_msa_equivalences.tsv — sparse columns for full structural MSA
  ssm_msa_equivalences.tsv         — sparse columns from MAFFT alignment (when MAFFT runs)
  ssm_msa.fasta / .aln             — MAFFT alignment of full structure sequences

MAFFT is resolved from ``MAFFT_ROOT``, common install prefixes, or ``PATH`` (see ``utils/msa_external.py``).

Run on one SSM output directory or scan many subdirectories with ``--dir``.

Examples (from repository root):

  python superimposition/SSM_struct_core.py path/to/SSMaligned_all_vs_all_LABEL/
  python superimposition/SSM_struct_core.py --dir path/to/parent --subdir-glob 'SSMaligned2_*'
  python superimposition/SSM_struct_core.py --equiv-dir path/to/ssm_equivalences
"""

from __future__ import annotations

import argparse
import fnmatch
import os
import sys

_REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if _REPO_ROOT not in sys.path:
    sys.path.insert(0, _REPO_ROOT)

from utils.cli_log import add_log_args, setup_log_from_args
from utils.msa_external import (
    build_label_residue_key_lists,
    mafft_alignment_to_columns,
    write_mafft_msa,
)
from utils.ssm_log_parse import (
    ResidueKey,
    SSMSuperpositionBlock,
    _natural_sort_key,
    blocks_from_equivalence_dir,
    build_equivalence_components,
    build_reference_anchored_msa_columns,
    collect_one_letter_codes,
    discover_coot_log,
    format_clustal_like,
    format_fasta,
    parse_ssm_log,
    read_fasta,
    validate_equivalence_sequence_mapping,
    write_sparse_equivalences_tsv,
)


def _resolved(path: str) -> str:
    return os.path.abspath(os.path.expanduser(path))


def load_ssm_blocks(
    directory: str,
    *,
    log_file: str | None = None,
    equiv_dir: str | None = None,
) -> tuple[list[SSMSuperpositionBlock], str | None]:
    directory = _resolved(directory)
    blocks: list[SSMSuperpositionBlock] = []

    if equiv_dir:
        blocks.extend(blocks_from_equivalence_dir(equiv_dir))
    else:
        default_equiv = os.path.join(directory, "ssm_equivalences")
        if os.path.isdir(default_equiv):
            blocks.extend(blocks_from_equivalence_dir(default_equiv))

    if log_file:
        log_path = _resolved(log_file)
    else:
        log_path = discover_coot_log(directory)

    if log_path and os.path.isfile(log_path):
        log_blocks = parse_ssm_log(log_path)
        if blocks:
            by_pair = {(b.moving_label, b.reference_label): b for b in blocks}
            merged: list[SSMSuperpositionBlock] = []
            for lb in log_blocks:
                key = (lb.moving_label, lb.reference_label)
                if key in by_pair:
                    tsv_block = by_pair[key]
                    lb.equivalences = tsv_block.equivalences or lb.equivalences
                merged.append(lb)
            blocks = merged
        else:
            blocks = log_blocks
    else:
        log_path = None

    return blocks, log_path


def _sequences_for_columns(
    labels: list[str],
    columns: list[dict[str, ResidueKey]],
    aa_codes: dict[tuple[str, tuple], str],
) -> list[str]:
    seqs = []
    for label in labels:
        chars: list[str] = []
        for col in columns:
            if label not in col:
                chars.append("-")
                continue
            node = (label, col[label].as_tuple())
            chars.append(aa_codes.get(node, "X"))
        seqs.append("".join(chars))
    return seqs


def _write_core_tsv(
    path: str,
    labels: list[str],
    core_columns: list[dict[str, ResidueKey]],
    aa_codes: dict[tuple[str, tuple], str],
    *,
    primary_reference: str | None = None,
) -> None:
    os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
    cols = ["column"]
    for label in labels:
        cols.extend(
            [
                "{}_chain".format(label),
                "{}_resnum".format(label),
                "{}_icode".format(label),
                "{}_aa".format(label),
            ]
        )
    header = (
        "# SSM structural core: reference residues in all model→reference alignments\n"
        if primary_reference
        else "# SSM structural core: residues present in all pairwise alignments\n"
    )
    with open(path, "w", encoding="utf-8") as fh:
        fh.write(header)
        fh.write("\t".join(cols) + "\n")
        for i, mapping in enumerate(core_columns, start=1):
            row = [str(i)]
            for label in labels:
                key = mapping[label]
                aa = aa_codes.get((label, key.as_tuple()), "X")
                row.extend([key.chain, str(key.resnum), key.icode, aa])
            fh.write("\t".join(row) + "\n")


def process_directory(
    directory: str,
    *,
    log_file: str | None = None,
    equiv_dir: str | None = None,
    skip_existing: bool = False,
    force: bool = False,
    run_mafft: bool = True,
) -> bool:
    directory = _resolved(directory)
    out_msa_fasta = os.path.join(directory, "ssm_structural_msa.fasta")
    out_msa_aln = os.path.join(directory, "ssm_structural_msa.aln")
    out_core_fasta = os.path.join(directory, "ssm_struct_core.fasta")
    out_core_aln = os.path.join(directory, "ssm_struct_core.aln")
    out_core_tsv = os.path.join(directory, "ssm_core_equivalences.tsv")
    out_struct_msa_tsv = os.path.join(directory, "ssm_structural_msa_equivalences.tsv")
    out_mafft_tsv = os.path.join(directory, "ssm_msa_equivalences.tsv")
    out_mafft_fasta = os.path.join(directory, "ssm_msa.fasta")
    out_mafft_aln = os.path.join(directory, "ssm_msa.aln")

    skip_markers = [out_core_fasta, out_core_tsv, out_msa_fasta, out_struct_msa_tsv]
    if run_mafft:
        skip_markers.append(out_mafft_tsv)

    if skip_existing and not force and all(
        os.path.isfile(p) for p in skip_markers
    ):
        print("SKIP (outputs exist): {}".format(directory))
        return True

    blocks, log_path = load_ssm_blocks(directory, log_file=log_file, equiv_dir=equiv_dir)
    if not blocks:
        print("WARNING: no SSM equivalences found in {}".format(directory), file=sys.stderr)
        return False

    for block in blocks:
        for warning in validate_equivalence_sequence_mapping(block):
            print("WARNING: {}".format(warning), file=sys.stderr)

    labels, all_cols, core_cols, _pairwise, primary_reference = build_equivalence_components(
        blocks, log_file=log_path
    )
    if len(labels) < 2:
        print(
            "WARNING: need at least two structures; found {} in {}".format(
                len(labels), directory
            ),
            file=sys.stderr,
        )
        return False

    core_labels = labels
    if primary_reference and core_cols:
        core_labels = [primary_reference] + sorted(
            (lab for lab in labels if lab != primary_reference),
            key=_natural_sort_key,
        )

    aa_codes = collect_one_letter_codes(blocks)

    all_seqs = _sequences_for_columns(labels, all_cols, aa_codes)
    core_seqs = _sequences_for_columns(core_labels, core_cols, aa_codes)

    with open(out_msa_fasta, "w", encoding="utf-8") as fh:
        fh.write(format_fasta(list(zip(labels, all_seqs))))
    with open(out_msa_aln, "w", encoding="utf-8") as fh:
        fh.write(
            format_clustal_like(
                labels,
                all_seqs,
                title="SSM structural MSA ({} structures, {} columns)".format(
                    len(labels), len(all_cols)
                ),
            )
        )

    with open(out_core_fasta, "w", encoding="utf-8") as fh:
        fh.write(format_fasta(list(zip(core_labels, core_seqs))))
    with open(out_core_aln, "w", encoding="utf-8") as fh:
        fh.write(
            format_clustal_like(
                core_labels,
                core_seqs,
                title="SSM structural core ({} structures, {} core columns)".format(
                    len(core_labels), len(core_cols)
                ),
            )
        )
    _write_core_tsv(
        out_core_tsv,
        core_labels,
        core_cols,
        aa_codes,
        primary_reference=primary_reference,
    )
    write_sparse_equivalences_tsv(
        out_struct_msa_tsv,
        labels,
        all_cols,
        aa_codes,
        header_comment=(
            "SSM structural MSA: union of pairwise alignment columns "
            "({} structures, {} columns)".format(len(labels), len(all_cols))
        ),
    )

    mafft_fallback: list[tuple[str, str]] | None = None
    if primary_reference:
        ref_labels, ref_cols = build_reference_anchored_msa_columns(
            primary_reference, blocks
        )
        ref_seqs = _sequences_for_columns(ref_labels, ref_cols, aa_codes)
        mafft_fallback = list(zip(ref_labels, ref_seqs))
    elif all_cols:
        mafft_fallback = list(zip(labels, all_seqs))

    mafft_ok = False
    if run_mafft:
        mafft_ok, mafft_msg = write_mafft_msa(
            directory,
            labels,
            blocks,
            log_path,
            out_fasta=out_mafft_fasta,
            out_aln=out_mafft_aln,
            fallback_records=mafft_fallback,
        )
        if not mafft_ok:
            print("WARNING: MAFFT failed: {}".format(mafft_msg), file=sys.stderr)
        elif os.path.isfile(out_mafft_fasta):
            aligned = read_fasta(out_mafft_fasta)
            key_lists = build_label_residue_key_lists(
                directory, labels, blocks, log_file=log_path
            )
            mafft_cols = mafft_alignment_to_columns(aligned, key_lists)
            if mafft_cols:
                write_sparse_equivalences_tsv(
                    out_mafft_tsv,
                    labels,
                    mafft_cols,
                    aa_codes,
                    header_comment=(
                        "MAFFT MSA columns mapped to native residues "
                        "({} structures, {} columns)".format(
                            len(labels), len(mafft_cols)
                        )
                    ),
                )
            else:
                print(
                    "WARNING: no MAFFT equivalence columns written for {}".format(
                        directory
                    ),
                    file=sys.stderr,
                )

    print("Processed: {}".format(directory))
    if primary_reference:
        n_models = len([lab for lab in labels if lab != primary_reference])
        print("  Reference: {} ({} models)".format(primary_reference, n_models))
    print("  Structures: {} ({})".format(len(labels), ", ".join(labels)))
    print("  MSA columns (union of pairwise alignments): {}".format(len(all_cols)))
    print("  Core columns: {}".format(len(core_cols)))
    print(
        "  Wrote: {}, {}, {}, {}".format(
            out_msa_fasta, out_core_fasta, out_core_tsv, out_struct_msa_tsv
        )
    )
    if run_mafft and os.path.isfile(out_mafft_fasta):
        print("  MAFFT: {}".format(out_mafft_fasta))
        if os.path.isfile(out_mafft_tsv):
            print("  MAFFT equivalences: {}".format(out_mafft_tsv))
    if core_cols:
        print(
            "  Next: python superimposition/SSM_aligned_core.py {}".format(
                directory
            )
        )
    return True


def discover_subdirs(base_dir: str, subdir_glob: str) -> list[str]:
    base_dir = _resolved(base_dir)
    matches: list[str] = []
    for root, dirnames, _files in os.walk(base_dir):
        for name in dirnames:
            if fnmatch.fnmatch(name, subdir_glob):
                matches.append(os.path.join(root, name))
    matches.sort(key=_natural_sort_key)
    return matches


def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Build structural MSA and core sequences from SSM pairwise equivalences "
            "(Coot log and/or ssm_equivalences/*.tsv). Optional MAFFT on full sequences."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    parser.add_argument(
        "directory",
        nargs="?",
        default=None,
        help="SSM output directory (contains coot_log*.txt and/or ssm_equivalences/).",
    )
    parser.add_argument(
        "--dir",
        dest="scan_dir",
        metavar="DIR",
        default=None,
        help="Scan under DIR for subdirectories to process (see --subdir-glob).",
    )
    parser.add_argument(
        "--subdir-glob",
        default="SSMaligned*",
        help="With --dir: process subdirs matching this glob (default: SSMaligned*).",
    )
    parser.add_argument(
        "--coot-log",
        dest="log_file",
        default=None,
        help="Explicit Coot log path (single-directory mode).",
    )
    parser.add_argument(
        "--equiv-dir",
        default=None,
        help="Directory of pairwise TSV files (default: <dir>/ssm_equivalences).",
    )
    parser.add_argument(
        "--skip-existing",
        action="store_true",
        help="Skip a directory when core/MSA outputs already exist.",
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="Re-run even when output files already exist.",
    )
    parser.add_argument(
        "--no-mafft",
        action="store_true",
        help="Skip external MAFFT; write structural outputs only.",
    )
    add_log_args(parser)
    args = parser.parse_args()

    if args.scan_dir and args.directory:
        parser.error("Use either a single directory or --dir, not both.")
    if not args.scan_dir and not args.directory:
        parser.error("Provide a directory or --dir=DIR.")

    setup_log_from_args(
        args,
        script_path=__file__,
        inputs=[args.directory or args.scan_dir or ""],
        pattern=args.subdir_glob,
    )

    targets: list[str] = []
    if args.scan_dir:
        targets = discover_subdirs(args.scan_dir, args.subdir_glob)
        if not targets:
            print(
                "No subdirectories matching {!r} under {}".format(
                    args.subdir_glob, _resolved(args.scan_dir)
                ),
                file=sys.stderr,
            )
            sys.exit(1)
    else:
        targets = [_resolved(args.directory)]

    ok = 0
    failed = 0
    for target in targets:
        if process_directory(
            target,
            log_file=args.log_file,
            equiv_dir=args.equiv_dir,
            skip_existing=args.skip_existing,
            force=args.force,
            run_mafft=not args.no_mafft,
        ):
            ok += 1
        else:
            failed += 1

    print("=" * 72)
    print("Done. Processed: {}, failed: {}".format(ok, failed))
    if failed:
        sys.exit(1)


if __name__ == "__main__":
    main()
