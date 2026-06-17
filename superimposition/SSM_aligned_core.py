#!/usr/bin/env python3
"""
Trim native PDBs to the SSM structural core or MSA-defined regions, re-superimpose
in Coot SSM, and superimpose full-length models onto each aligned trimmed model.

Reads equivalence TSV from an ``SSMaligned*`` directory (after ``SSM_struct_core.py``).
Use ``--trim-source`` to select core (default), structural MSA, continuous MSA,
or MAFFT MSA columns.

Writes a sibling directory (name depends on trim source):

  ../SSM_aligned_core_<run_label>/              (--trim-source core)
  ../SSM_aligned_structural_msa_<run_label>/    (--trim-source structural_msa)
  ../SSM_aligned_continuous_msa_<run_label>/   (--trim-source continuous_msa)
  ../SSM_aligned_mafft_msa_<run_label>/         (--trim-source mafft)
    core_trim/{label}_core.pdb
    continuous_trim/{label}_cont.pdb   (continuous_msa only)
    core_ssm/               # Coot SSM: trimmed cores → anchor core
    full_ssm/               # Coot SSM: full native → own aligned core
    ssm_core_alignment_metrics.tsv
    anchor.txt
    processing_log.txt

Phase A metrics (core RMSD, % sequence identity) come from the phase A Coot log.

Requires BioPython (trim) and Coot on PATH (SSM phases).

Examples (from repository root):

  python superimposition/SSM_aligned_core.py path/to/SSMaligned_all_vs_all_<run_label>/
  python superimposition/SSM_aligned_core.py path/to/SSM/ --trim-source structural_msa
  python superimposition/SSM_aligned_core.py path/to/SSM/ --trim-source continuous_msa
  python superimposition/SSM_aligned_core.py --dir path/to/parent --subdir-glob 'SSMaligned_all_vs_all_*'
"""

from __future__ import annotations

import argparse
import fnmatch
import os
import shutil
import sys

_REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if _REPO_ROOT not in sys.path:
    sys.path.insert(0, _REPO_ROOT)

from ranking.extract_rmsd import extract_ssm_rmsd_and_seq, extract_ssm_rmsd_values
from superimposition.aligned_pdb_names import aligned_pdb_basename
from superimposition.coot_ssm_templates import (
    create_full_to_core_ssm_script,
    create_one_to_many_ssm_script,
)
from utils.cli_log import add_log_args, setup_log_from_args
from utils.coot_run import coot_available, run_coot_script
from utils.pdb_core_trim import (
    BIOPYTHON_AVAILABLE,
    count_core_ca_present,
    extract_core_structure,
    find_native_structure,
    longest_contiguous_residue_keys,
    write_structure_pdb,
)
from utils.ssm_log_parse import (
    ResidueKey,
    TRIM_SOURCE_CONTINUOUS_MSA,
    TRIM_SOURCE_CORE,
    TRIM_SOURCE_MAFFT,
    TRIM_SOURCE_STRUCTURAL_MSA,
    _natural_sort_key,
    aligned_trim_output_dir,
    audit_coot_ssm_log,
    coot_metrics_identity_summary,
    coot_metrics_rmsd_summary,
    discover_coot_log,
    metrics_rows_have_rmsd,
    pairwise_metrics_from_ssm_log,
    parse_log_models_directory,
    read_fasta,
    read_name_map,
    read_sparse_equivalences_tsv,
    resolve_models_directory,
    ssm_output_run_label,
    trim_equivalences_tsv_path,
)

SPARSE_TRIM_CA_THRESHOLD = 20

_ALIGNED_TAG = "_SSMaligned2_"
_FULL_TAG = "_SSMfull2core_"
_NO_TRIM_MSG = "no trim defined"
_NO_CORE_MSG = "no core defined"

TRIM_SOURCE_CHOICES = (
    TRIM_SOURCE_CORE,
    TRIM_SOURCE_STRUCTURAL_MSA,
    TRIM_SOURCE_CONTINUOUS_MSA,
    TRIM_SOURCE_MAFFT,
)


def _trim_stem_suffix(trim_source: str) -> str:
    if trim_source == TRIM_SOURCE_CORE:
        return "_core"
    if trim_source == TRIM_SOURCE_CONTINUOUS_MSA:
        return "_cont"
    return "_trim"


def _trim_subdir_names(trim_source: str) -> tuple[str, str]:
    if trim_source == TRIM_SOURCE_CORE:
        return "core_trim", "core_ssm"
    if trim_source == TRIM_SOURCE_CONTINUOUS_MSA:
        return "continuous_trim", "continuous_ssm"
    return "msa_trim", "msa_ssm"


def _empty_trim_message(trim_source: str) -> str:
    if trim_source == TRIM_SOURCE_CORE:
        return _NO_CORE_MSG
    if trim_source == TRIM_SOURCE_CONTINUOUS_MSA:
        return "no continuous trim defined"
    return _NO_TRIM_MSG


def _resolved(path: str) -> str:
    return os.path.abspath(os.path.expanduser(path))


def _count_pdb_ca_atoms(path: str) -> int:
    """Count Cα atoms in a PDB file (fast line scan)."""
    n = 0
    with open(path, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if not line.startswith(("ATOM  ", "HETATM")):
                continue
            if line[12:16].strip() == "CA":
                n += 1
    return n


def _trim_keys_for_label(
    label: str,
    columns: list[dict[str, ResidueKey]],
    *,
    structure_path: str | None = None,
    trim_source: str = TRIM_SOURCE_CORE,
) -> list[ResidueKey]:
    """Ordered residue keys for trim (MSA column order; deduplicated)."""
    seen: set[tuple[str, int, str]] = set()
    keys: list[ResidueKey] = []
    for col in columns:
        if label not in col:
            continue
        kt = col[label].as_tuple()
        if kt in seen:
            continue
        seen.add(kt)
        keys.append(col[label])
    if trim_source == TRIM_SOURCE_CONTINUOUS_MSA and structure_path:
        keys = longest_contiguous_residue_keys(structure_path, keys)
    return keys


def _core_keys_for_label(
    label: str,
    core_columns: list[dict[str, ResidueKey]],
) -> list[ResidueKey]:
    return _trim_keys_for_label(label, core_columns)


def _native_pdb_for_label(
    models_dir: str | None,
    label: str,
    *,
    name_map: dict[str, str] | None = None,
) -> str | None:
    if not models_dir:
        return None
    return find_native_structure(models_dir, label, name_map=name_map)


def _score_anchor_candidate(
    anchor: str,
    labels: list[str],
    core_columns: list[dict[str, ResidueKey]],
    models_dir: str | None,
    *,
    name_map: dict[str, str] | None = None,
) -> tuple[int, dict[str, str]]:
    """Return (complete_core_columns, native_source_map) for anchor candidate."""
    sources: dict[str, str | None] = {}
    for label in labels:
        sources[label] = _native_pdb_for_label(
            models_dir, label, name_map=name_map
        )

    complete = 0
    for col in core_columns:
        ok = True
        for label in labels:
            src = sources.get(label)
            if not src:
                ok = False
                break
            key = col[label]
            if count_core_ca_present(src, [key]) < 1:
                ok = False
                break
        if ok:
            complete += 1

    resolved = {k: v for k, v in sources.items() if v}
    return complete, resolved


def _score_anchor_msa_candidate(
    anchor: str,
    labels: list[str],
    columns: list[dict[str, ResidueKey]],
    models_dir: str | None,
    *,
    name_map: dict[str, str] | None = None,
    trim_source: str = TRIM_SOURCE_STRUCTURAL_MSA,
) -> tuple[int, dict[str, str]]:
    """Return (trimmed_residue_count, native_source_map) for MSA-based trim."""
    sources: dict[str, str | None] = {}
    for label in labels:
        sources[label] = _native_pdb_for_label(
            models_dir, label, name_map=name_map
        )

    src = sources.get(anchor)
    keys = _trim_keys_for_label(
        anchor,
        columns,
        structure_path=src,
        trim_source=trim_source,
    )
    if not src or not keys:
        return 0, {k: v for k, v in sources.items() if v}

    n_found = count_core_ca_present(src, keys)
    if n_found != len(keys):
        return 0, {k: v for k, v in sources.items() if v}

    return len(keys), {k: v for k, v in sources.items() if v}


def _score_anchor_for_trim_source(
    trim_source: str,
    anchor: str,
    labels: list[str],
    columns: list[dict[str, ResidueKey]],
    models_dir: str | None,
    *,
    name_map: dict[str, str] | None = None,
) -> tuple[int, dict[str, str]]:
    if trim_source == TRIM_SOURCE_CORE:
        return _score_anchor_candidate(
            anchor, labels, columns, models_dir, name_map=name_map
        )
    return _score_anchor_msa_candidate(
        anchor,
        labels,
        columns,
        models_dir,
        name_map=name_map,
        trim_source=trim_source,
    )


def choose_anchor(
    labels: list[str],
    columns: list[dict[str, ResidueKey]],
    models_dir: str | None,
    anchor_override: str | None,
    *,
    name_map: dict[str, str] | None = None,
    trim_source: str = TRIM_SOURCE_CORE,
) -> tuple[str, dict[str, str], str]:
    """Pick anchor; return (anchor, sources, reason)."""
    if trim_source == TRIM_SOURCE_CORE:
        score_label = "complete core columns"
    elif trim_source == TRIM_SOURCE_CONTINUOUS_MSA:
        score_label = "continuous trim residues"
    else:
        score_label = "trim residues"

    if anchor_override:
        if anchor_override not in labels:
            raise ValueError(
                "Anchor {!r} not among structures: {}".format(
                    anchor_override, ", ".join(labels)
                )
            )
        score, sources = _score_anchor_for_trim_source(
            trim_source,
            anchor_override,
            labels,
            columns,
            models_dir,
            name_map=name_map,
        )
        reason = "user --anchor ({} {})".format(score, score_label)
        return anchor_override, sources, reason

    best_anchor = labels[0]
    best_score = -1
    best_sources: dict[str, str] = {}
    for candidate in sorted(labels, key=_natural_sort_key):
        score, sources = _score_anchor_for_trim_source(
            trim_source,
            candidate,
            labels,
            columns,
            models_dir,
            name_map=name_map,
        )
        if score > best_score or (
            score == best_score
            and _natural_sort_key(candidate) < _natural_sort_key(best_anchor)
        ):
            best_anchor = candidate
            best_score = score
            best_sources = sources

    reason = "auto: {} {} (max over candidates)".format(best_score, score_label)
    return best_anchor, best_sources, reason


def _aligned_trim_path(
    ssm_dir: str, label: str, anchor: str, stem_suffix: str
) -> str:
    if label == anchor:
        return os.path.join(ssm_dir, "{}{}.pdb".format(label, stem_suffix))
    model_stem = "{}{}".format(label, stem_suffix)
    ref_stem = "{}{}".format(anchor, stem_suffix)
    name = aligned_pdb_basename(model_stem, ref_stem, _ALIGNED_TAG)
    return os.path.join(ssm_dir, name)


def _aligned_core_path(core_ssm_dir: str, label: str, anchor: str) -> str:
    return _aligned_trim_path(core_ssm_dir, label, anchor, "_core")


def _write_coot_metrics_tables(
    out_dir: str,
    metric_rows: list[dict[str, object]],
    anchor: str,
    *,
    no_core: bool = False,
    structure_labels: list[str] | None = None,
) -> None:
    metrics_path = os.path.join(out_dir, "ssm_core_alignment_metrics.tsv")
    pair_path = os.path.join(out_dir, "ssm_core_sequence_identity.tsv")
    summary_path = os.path.join(out_dir, "ssm_core_sequence_identity_summary.tsv")
    rmsd_summary_path = os.path.join(out_dir, "ssm_core_rmsd_summary.tsv")

    if no_core:
        _write_no_core_table_rows(
            metrics_path,
            pair_path,
            summary_path,
            rmsd_summary_path,
            structure_labels or [],
        )
        return

    with open(metrics_path, "w", encoding="utf-8") as fh:
        fh.write(
            "structure_a\tstructure_b\tn_aligned_residues\t"
            "pct_identity\tcore_rmsd_angstrom\n"
        )
        for row in metric_rows:
            fh.write(
                "{structure_a}\t{structure_b}\t{n_aligned_residues}\t"
                "{pct_identity}\t{core_rmsd_angstrom}\n".format(**row)
            )

    with open(pair_path, "w", encoding="utf-8") as fh:
        fh.write(
            "structure_a\tstructure_b\tn_aligned_residues\tpct_identity\n"
        )
        for row in metric_rows:
            fh.write(
                "{structure_a}\t{structure_b}\t{n_aligned_residues}\t"
                "{pct_identity}\n".format(
                    structure_a=row["structure_a"],
                    structure_b=row["structure_b"],
                    n_aligned_residues=row["n_aligned_residues"],
                    pct_identity=row["pct_identity"],
                )
            )

    id_sum = coot_metrics_identity_summary(metric_rows, anchor=anchor)
    with open(summary_path, "w", encoding="utf-8") as fh:
        fh.write("structure\tn_comparisons\tmean_pct_identity\n")
        for row in id_sum:
            fh.write(
                "{structure}\t{n_comparisons}\t{mean_pct_identity}\n".format(
                    **row
                )
            )

    rmsd_sum = coot_metrics_rmsd_summary(metric_rows, anchor=anchor)
    with open(rmsd_summary_path, "w", encoding="utf-8") as fh:
        fh.write("structure\tn_comparisons\tmean_core_rmsd_angstrom\n")
        for row in rmsd_sum:
            fh.write(
                "{structure}\t{n_comparisons}\t{mean_core_rmsd_angstrom}\n".format(
                    **row
                )
            )


def _no_core_row(
    structure_a: str = "",
    structure_b: str = "",
) -> dict[str, str]:
    msg = _NO_CORE_MSG
    return {
        "structure_a": structure_a,
        "structure_b": structure_b,
        "n_aligned_residues": msg,
        "pct_identity": msg,
        "core_rmsd_angstrom": msg,
    }


def _no_core_summary_row(structure: str = "") -> dict[str, str]:
    return {
        "structure": structure or _NO_CORE_MSG,
        "n_comparisons": "0",
        "mean_pct_identity": _NO_CORE_MSG,
    }


def _no_core_rmsd_summary_row(structure: str = "") -> dict[str, str]:
    return {
        "structure": structure or _NO_CORE_MSG,
        "n_comparisons": "0",
        "mean_core_rmsd_angstrom": _NO_CORE_MSG,
    }


def _write_no_core_table_rows(
    metrics_path: str,
    pair_path: str,
    summary_path: str,
    rmsd_summary_path: str,
    structure_labels: list[str],
) -> None:
    pair_row = _no_core_row()
    id_row = {
        k: v
        for k, v in pair_row.items()
        if k != "core_rmsd_angstrom"
    }

    with open(metrics_path, "w", encoding="utf-8") as fh:
        fh.write("# {}\n".format(_NO_CORE_MSG))
        fh.write(
            "structure_a\tstructure_b\tn_aligned_residues\t"
            "pct_identity\tcore_rmsd_angstrom\n"
        )
        fh.write(
            "{structure_a}\t{structure_b}\t{n_aligned_residues}\t"
            "{pct_identity}\t{core_rmsd_angstrom}\n".format(**pair_row)
        )

    with open(pair_path, "w", encoding="utf-8") as fh:
        fh.write("# {}\n".format(_NO_CORE_MSG))
        fh.write(
            "structure_a\tstructure_b\tn_aligned_residues\tpct_identity\n"
        )
        fh.write(
            "{structure_a}\t{structure_b}\t{n_aligned_residues}\t"
            "{pct_identity}\n".format(**id_row)
        )

    summary_rows = (
        [_no_core_summary_row(lab) for lab in structure_labels]
        if structure_labels
        else [_no_core_summary_row()]
    )
    with open(summary_path, "w", encoding="utf-8") as fh:
        fh.write("# {}\n".format(_NO_CORE_MSG))
        fh.write("structure\tn_comparisons\tmean_pct_identity\n")
        for row in summary_rows:
            fh.write(
                "{structure}\t{n_comparisons}\t{mean_pct_identity}\n".format(
                    **row
                )
            )

    rmsd_rows = (
        [_no_core_rmsd_summary_row(lab) for lab in structure_labels]
        if structure_labels
        else [_no_core_rmsd_summary_row()]
    )
    with open(rmsd_summary_path, "w", encoding="utf-8") as fh:
        fh.write("# {}\n".format(_NO_CORE_MSG))
        fh.write("structure\tn_comparisons\tmean_core_rmsd_angstrom\n")
        for row in rmsd_rows:
            fh.write(
                "{structure}\t{n_comparisons}\t{mean_core_rmsd_angstrom}\n".format(
                    **row
                )
            )


def _process_no_trim_directory(
    ssm_dir: str,
    out_dir: str,
    run_label: str,
    equiv_tsv: str,
    labels: list[str],
    sequences: dict[str, str],
    *,
    trim_source: str,
) -> bool:
    empty_msg = _empty_trim_message(trim_source)
    structure_labels = labels or sorted(sequences.keys(), key=_natural_sort_key)
    os.makedirs(out_dir, exist_ok=True)

    log_lines = [
        "# SSM aligned trim processing",
        "trim_source: {}".format(trim_source),
        "ssm_input_dir: {}".format(ssm_dir),
        "run_label: {}".format(run_label),
        "output_dir: {}".format(out_dir),
        "status: {}".format(empty_msg),
        "equivalences_tsv: {}".format(equiv_tsv),
        "n_msa_columns: 0",
        "",
    ]

    with open(os.path.join(out_dir, "anchor.txt"), "w", encoding="utf-8") as fh:
        fh.write("anchor: (none)\n")
        fh.write("{}\n".format(empty_msg))

    with open(
        os.path.join(out_dir, "processing_log.txt"), "w", encoding="utf-8"
    ) as fh:
        fh.write("\n".join(log_lines) + "\n")

    _write_coot_metrics_tables(
        out_dir, [], "", no_core=True, structure_labels=structure_labels
    )

    print("Processed: {} ({})".format(ssm_dir, empty_msg))
    print("  Output: {}".format(out_dir))
    return True


def _process_no_core_directory(
    ssm_dir: str,
    out_dir: str,
    run_label: str,
    core_tsv: str,
    labels: list[str],
    sequences: dict[str, str],
) -> bool:
    return _process_no_trim_directory(
        ssm_dir,
        out_dir,
        run_label,
        core_tsv,
        labels,
        sequences,
        trim_source=TRIM_SOURCE_CORE,
    )


def _aligned_output_exists(ssm_log_marker: str) -> bool:
    return os.path.isfile(ssm_log_marker)


def _clear_aligned_output(
    out_dir: str,
    trim_dir_name: str,
    ssm_dir_name: str,
) -> None:
    """Remove phase outputs so a forced re-run starts clean."""
    for sub in (trim_dir_name, ssm_dir_name, "full_ssm"):
        path = os.path.join(out_dir, sub)
        if os.path.isdir(path):
            shutil.rmtree(path)
    for name in (
        "ssm_core_alignment_metrics.tsv",
        "ssm_core_sequence_identity.tsv",
        "ssm_core_sequence_identity_summary.tsv",
        "ssm_core_rmsd_summary.tsv",
        "anchor.txt",
        "processing_log.txt",
    ):
        path = os.path.join(out_dir, name)
        if os.path.isfile(path):
            os.unlink(path)


def _phase_a_failed(
    metric_rows: list[dict[str, object]],
    expected_pairs: int,
    coot_failures: list[str],
) -> bool:
    if coot_failures:
        return True
    if expected_pairs <= 0:
        return False
    if len(metric_rows) < expected_pairs:
        return True
    return not metrics_rows_have_rmsd(metric_rows)


def _run_phase_a_coot(
    trim_dir: str,
    ssm_out_dir: str,
    anchor: str,
    labels: list[str],
    *,
    stem_suffix: str,
    not_interactive: bool,
    log_lines: list[str],
) -> tuple[int, str]:
    """Coot SSM one-to-many on trimmed structures; anchor defines the frame."""
    os.makedirs(ssm_out_dir, exist_ok=True)
    anchor_trim = os.path.join(trim_dir, "{}{}.pdb".format(anchor, stem_suffix))
    anchor_dest = os.path.join(ssm_out_dir, "{}{}.pdb".format(anchor, stem_suffix))
    shutil.copy2(anchor_trim, anchor_dest)

    moving = [
        os.path.join(trim_dir, "{}{}.pdb".format(lab, stem_suffix))
        for lab in sorted(labels, key=_natural_sort_key)
        if lab != anchor
    ]

    log_path = os.path.join(ssm_out_dir, "coot_log_core.txt")
    script = create_one_to_many_ssm_script(
        anchor_trim,
        moving,
        ssm_out_dir,
        keep_coot_open=not not_interactive,
        aligned_tag=_ALIGNED_TAG,
    )
    header = [
        "# SSM aligned trim — phase A (trimmed → anchor)",
        "# Reference: {} (native trimmed)".format(anchor_trim),
        "# Anchor: {}".format(anchor),
        "# Models: {}".format(len(moving)),
    ]
    log_lines.append("phase_A: coot_log={}".format(log_path))
    rc = run_coot_script(
        script, log_path, header_lines=header, not_interactive=not_interactive
    )
    return rc, log_path


def _run_phase_b_coot(
    ssm_out_dir: str,
    full_ssm_dir: str,
    models_dir: str,
    anchor: str,
    labels: list[str],
    *,
    stem_suffix: str,
    not_interactive: bool,
    log_lines: list[str],
    name_map: dict[str, str] | None = None,
) -> tuple[int, str]:
    """Coot SSM: each full native structure onto its own aligned trim."""
    os.makedirs(full_ssm_dir, exist_ok=True)
    pairs: list[tuple[str, str, str]] = []
    for label in sorted(labels, key=_natural_sort_key):
        ref_trim = _aligned_trim_path(ssm_out_dir, label, anchor, stem_suffix)
        full_native = find_native_structure(
            models_dir, label, name_map=name_map
        )
        if not ref_trim or not os.path.isfile(ref_trim):
            print(
                "ERROR: missing aligned trim for {}: {}".format(label, ref_trim),
                file=sys.stderr,
            )
            return 1, ""
        if not full_native or not os.path.isfile(full_native):
            print(
                "ERROR: missing native PDB for {} under {}".format(
                    label, models_dir
                ),
                file=sys.stderr,
            )
            return 1, ""
        pairs.append((ref_trim, full_native, label))

    log_path = os.path.join(full_ssm_dir, "coot_log_full.txt")
    script = create_full_to_core_ssm_script(
        pairs,
        full_ssm_dir,
        keep_coot_open=not not_interactive,
        aligned_tag=_FULL_TAG,
    )
    header = [
        "# SSM aligned trim — phase B (full native → own aligned trim)",
        "# Pairs: {}".format(len(pairs)),
    ]
    log_lines.append("phase_B: coot_log={}".format(log_path))
    rc = run_coot_script(
        script, log_path, header_lines=header, not_interactive=not_interactive
    )
    if rc == 0:
        extract_ssm_rmsd_values(
            log_path,
            rmsd_output_path=os.path.join(
                full_ssm_dir, "rmsd_SSM_values_full.txt"
            ),
        )
    return rc, log_path


def process_ssm_directory(
    ssm_dir: str,
    *,
    anchor_override: str | None = None,
    skip_existing: bool = False,
    force: bool = False,
    not_interactive: bool = True,
    skip_full_ssm: bool = False,
    trim_source: str = TRIM_SOURCE_CORE,
) -> bool:
    ssm_dir = _resolved(ssm_dir)
    out_dir = aligned_trim_output_dir(ssm_dir, trim_source)
    run_label = ssm_output_run_label(ssm_dir)
    stem_suffix = _trim_stem_suffix(trim_source)
    trim_dir_name, ssm_dir_name = _trim_subdir_names(trim_source)

    equiv_tsv = trim_equivalences_tsv_path(ssm_dir, trim_source)
    fallback_fasta = os.path.join(ssm_dir, "ssm_struct_core.fasta")
    if trim_source in (TRIM_SOURCE_STRUCTURAL_MSA, TRIM_SOURCE_CONTINUOUS_MSA):
        fallback_fasta = os.path.join(ssm_dir, "ssm_structural_msa.fasta")
    elif trim_source == TRIM_SOURCE_MAFFT:
        fallback_fasta = os.path.join(ssm_dir, "ssm_msa.fasta")

    ssm_out_dir = os.path.join(out_dir, ssm_dir_name)
    ssm_log_marker = os.path.join(ssm_out_dir, "coot_log_core.txt")

    if skip_existing and _aligned_output_exists(ssm_log_marker):
        print("SKIP (outputs exist): {}".format(out_dir))
        return True

    if _aligned_output_exists(ssm_log_marker) and not force:
        print(
            "ERROR: output exists at {} — use --force to overwrite".format(
                out_dir
            ),
            file=sys.stderr,
        )
        return False

    if force and os.path.isdir(out_dir):
        _clear_aligned_output(out_dir, trim_dir_name, ssm_dir_name)

    if not os.path.isfile(equiv_tsv):
        hint = "run SSM_struct_core.py first"
        if trim_source == TRIM_SOURCE_MAFFT:
            hint = "run SSM_struct_core.py without --no-mafft first"
        print(
            "WARNING: missing {} — {}".format(equiv_tsv, hint),
            file=sys.stderr,
        )
        return False

    if not BIOPYTHON_AVAILABLE:
        print("ERROR: BioPython required (pip install biopython)", file=sys.stderr)
        return False

    labels, trim_columns = read_sparse_equivalences_tsv(equiv_tsv)
    sequences: dict[str, str] = {}
    if os.path.isfile(fallback_fasta):
        sequences = dict(read_fasta(fallback_fasta))

    if not trim_columns:
        if not labels and sequences:
            labels = sorted(sequences.keys(), key=_natural_sort_key)
        empty_msg = _empty_trim_message(trim_source)
        print("NOTE: {} — {}".format(empty_msg, equiv_tsv), file=sys.stderr)
        return _process_no_trim_directory(
            ssm_dir,
            out_dir,
            run_label,
            equiv_tsv,
            labels,
            sequences,
            trim_source=trim_source,
        )

    if not coot_available():
        print("ERROR: coot not found on PATH", file=sys.stderr)
        return False

    if not labels:
        print("WARNING: no structure labels in {}".format(equiv_tsv), file=sys.stderr)
        return False

    log_path = discover_coot_log(ssm_dir)
    models_dir, models_dir_source = resolve_models_directory(ssm_dir, log_path)
    if not models_dir:
        print(
            "ERROR: models directory not found for {} "
            "(checked coot log, re-anchored path, _short_<run_label>, name map)".format(
                ssm_dir
            ),
            file=sys.stderr,
        )
        return False

    name_map: dict[str, str] = {}
    if models_dir_source == "name_map_parent":
        name_map_path = os.path.join(
            os.path.dirname(ssm_dir), "{}_name_map.tsv".format(run_label)
        )
        name_map = read_name_map(name_map_path)

    if models_dir_source != "coot_log":
        log_parsed = (
            parse_log_models_directory(log_path) if log_path else None
        )
        print(
            "WARNING: coot log # Directories: path not usable ({!r}); "
            "using {} ({})".format(
                log_parsed, models_dir, models_dir_source
            ),
            file=sys.stderr,
        )

    try:
        anchor, sources, anchor_reason = choose_anchor(
            labels,
            trim_columns,
            models_dir,
            anchor_override,
            name_map=name_map,
            trim_source=trim_source,
        )
    except ValueError as exc:
        print("ERROR: {}".format(exc), file=sys.stderr)
        return False

    trim_dir = os.path.join(out_dir, trim_dir_name)
    os.makedirs(out_dir, exist_ok=True)
    os.makedirs(trim_dir, exist_ok=True)

    col_label = (
        "n_core_columns"
        if trim_source == TRIM_SOURCE_CORE
        else "n_msa_columns"
    )
    log_lines: list[str] = [
        "# SSM aligned trim processing (Coot SSM re-alignment)",
        "trim_source: {}".format(trim_source),
        "ssm_input_dir: {}".format(ssm_dir),
        "run_label: {}".format(run_label),
        "output_dir: {}".format(out_dir),
        "anchor: {}".format(anchor),
        "anchor_selection: {}".format(anchor_reason),
        "models_directory: {}".format(models_dir),
        "models_directory_source: {}".format(models_dir_source),
        "{}: {}".format(col_label, len(trim_columns)),
        "equivalences_tsv: {}".format(equiv_tsv),
    ]
    if name_map:
        log_lines.append("name_map_entries: {}".format(len(name_map)))
    log_lines.extend([
        "",
        "# Phase 0: trim native PDBs",
    ])

    failed = False
    for label in sorted(labels, key=_natural_sort_key):
        src = sources.get(label) or _native_pdb_for_label(
            models_dir, label, name_map=name_map
        )
        if src and os.path.isfile(src):
            sources[label] = src
        if not src or not os.path.isfile(src):
            print("ERROR: no native PDB for {}".format(label), file=sys.stderr)
            failed = True
            continue

        n_msa_keys = len(_trim_keys_for_label(label, trim_columns))
        trim_keys = _trim_keys_for_label(
            label,
            trim_columns,
            structure_path=src,
            trim_source=trim_source,
        )
        if trim_source == TRIM_SOURCE_CONTINUOUS_MSA and n_msa_keys:
            log_lines.append(
                "{}: msa_matched={} continuous={}".format(
                    label, n_msa_keys, len(trim_keys)
                )
            )

        if not trim_keys:
            print("ERROR: no trim residues defined for {}".format(label), file=sys.stderr)
            failed = True
            continue

        n_found = count_core_ca_present(src, trim_keys)
        log_lines.append(
            "{}: native={} trim_residues={} ca_found={}".format(
                label, src, len(trim_keys), n_found
            )
        )
        if n_found != len(trim_keys):
            print(
                "ERROR: {} missing {}/{} trim Cα in {}".format(
                    label, len(trim_keys) - n_found, len(trim_keys), src
                ),
                file=sys.stderr,
            )
            failed = True
            continue

        try:
            struct = extract_core_structure(
                src, trim_keys, out_chain="A", renumber_from=1
            )
        except KeyError as exc:
            print("ERROR: {} — {}".format(label, exc), file=sys.stderr)
            failed = True
            continue

        trim_path = os.path.join(trim_dir, "{}{}.pdb".format(label, stem_suffix))
        write_structure_pdb(struct, trim_path)
        log_lines.append("  trim: {}".format(trim_path))

    if failed:
        print("ERROR: phase 0 trim failed for {}".format(ssm_dir), file=sys.stderr)
        return False

    log_lines.append("")
    log_lines.append("# Trim CA counts (sparse MSA trims may break SSM graph)")
    for label in sorted(labels, key=_natural_sort_key):
        trim_path = os.path.join(trim_dir, "{}{}.pdb".format(label, stem_suffix))
        n_ca = _count_pdb_ca_atoms(trim_path)
        log_lines.append("{}: trim_ca={}".format(label, n_ca))
        if trim_source in (
            TRIM_SOURCE_STRUCTURAL_MSA,
            TRIM_SOURCE_CONTINUOUS_MSA,
        ) and n_ca < SPARSE_TRIM_CA_THRESHOLD:
            mode = (
                "continuous_msa"
                if trim_source == TRIM_SOURCE_CONTINUOUS_MSA
                else "structural_msa"
            )
            print(
                "WARNING: {} trim has only {} CA atoms (< {}); "
                "{} Coot SSM may fail (can't make graph)".format(
                    label, n_ca, SPARSE_TRIM_CA_THRESHOLD, mode
                ),
                file=sys.stderr,
            )
            log_lines.append(
                "  warning: {} sparse trim ({} CA)".format(label, n_ca)
            )

    phase_a_failures: list[str] = []
    expected_pairs = max(0, len(labels) - 1)

    log_lines.append("")
    log_lines.append("# Phase A: Coot SSM trimmed → anchor")
    rc_a, coot_log_a = _run_phase_a_coot(
        trim_dir,
        ssm_out_dir,
        anchor,
        labels,
        stem_suffix=stem_suffix,
        not_interactive=not_interactive,
        log_lines=log_lines,
    )
    coot_failures = audit_coot_ssm_log(coot_log_a)
    if coot_failures:
        phase_a_failures.extend(coot_failures)
        for msg in coot_failures:
            print("WARNING: Coot SSM — {}".format(msg), file=sys.stderr)
    if rc_a != 0:
        print(
            "ERROR: Coot phase A failed (exit {}) — see {}".format(
                rc_a, coot_log_a
            ),
            file=sys.stderr,
        )
        with open(
            os.path.join(out_dir, "processing_log.txt"), "w", encoding="utf-8"
        ) as fh:
            fh.write("\n".join(log_lines) + "\n")
        return False

    extract_ssm_rmsd_and_seq(
        coot_log_a,
        output_dir=ssm_out_dir,
        rmsd_output_path=os.path.join(
            ssm_out_dir, "rmsd_SSM_values_core.txt"
        ),
        seq_align_output_path=os.path.join(
            ssm_out_dir, "ssm_seq_align_core.txt"
        ),
        equivalences_dir=os.path.join(ssm_out_dir, "ssm_equivalences"),
    )
    metric_rows = pairwise_metrics_from_ssm_log(coot_log_a)
    if _phase_a_failed(metric_rows, expected_pairs, coot_failures):
        print(
            "ERROR: Coot phase A incomplete or invalid — see {}".format(
                coot_log_a
            ),
            file=sys.stderr,
        )
        log_lines.append(
            "phase_A_status: failed (SSM errors or empty/incomplete metrics)"
        )
        with open(
            os.path.join(out_dir, "processing_log.txt"), "w", encoding="utf-8"
        ) as fh:
            fh.write("\n".join(log_lines) + "\n")
        return False

    _write_coot_metrics_tables(out_dir, metric_rows, anchor)

    if not skip_full_ssm:
        log_lines.append("")
        log_lines.append("# Phase B: Coot SSM full native → own aligned trim")
        full_ssm_dir = os.path.join(out_dir, "full_ssm")
        rc_b, coot_log_b = _run_phase_b_coot(
            ssm_out_dir,
            full_ssm_dir,
            models_dir,
            anchor,
            labels,
            stem_suffix=stem_suffix,
            not_interactive=not_interactive,
            log_lines=log_lines,
            name_map=name_map or None,
        )
        phase_b_failures = audit_coot_ssm_log(coot_log_b)
        for msg in phase_b_failures:
            print("WARNING: Coot SSM phase B — {}".format(msg), file=sys.stderr)
        if rc_b != 0 or phase_b_failures:
            print(
                "ERROR: Coot phase B failed (exit {}) — see {}".format(
                    rc_b, coot_log_b
                ),
                file=sys.stderr,
            )
            with open(
                os.path.join(out_dir, "processing_log.txt"), "w", encoding="utf-8"
            ) as fh:
                fh.write("\n".join(log_lines) + "\n")
            return False
    else:
        log_lines.append("")
        log_lines.append("# Phase B skipped (--skip-full-ssm)")

    with open(os.path.join(out_dir, "anchor.txt"), "w", encoding="utf-8") as fh:
        fh.write("anchor: {}\n".format(anchor))
        fh.write("{}\n".format(anchor_reason))
        fh.write("\ntrim_source: {}\n".format(trim_source))
        fh.write("\nNative models directory:\n  {}\n".format(models_dir))
        fh.write("\nPhase A aligned trims:\n  {}\n".format(ssm_out_dir))
        if not skip_full_ssm:
            fh.write("\nPhase B full aligned:\n  {}\n".format(
                os.path.join(out_dir, "full_ssm")
            ))

    with open(
        os.path.join(out_dir, "processing_log.txt"), "w", encoding="utf-8"
    ) as fh:
        fh.write("\n".join(log_lines) + "\n")

    print("Processed: {}".format(ssm_dir))
    print("  Run label: {}".format(run_label))
    print("  Trim source: {}".format(trim_source))
    print("  Output: {}".format(out_dir))
    print("  Anchor: {} ({})".format(anchor, anchor_reason))
    print("  Trim columns: {}".format(len(trim_columns)))
    print("  Phase A pairs: {}".format(len(metric_rows)))
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
            "Trim native PDBs to the SSM structural core, re-superimpose in "
            "Coot SSM (anchor frame), and superimpose full models onto aligned cores."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    parser.add_argument(
        "directory",
        nargs="?",
        default=None,
        help="SSM output directory (SSMaligned* with ssm_core_equivalences.tsv).",
    )
    parser.add_argument(
        "--dir",
        dest="scan_dir",
        metavar="DIR",
        default=None,
        help="Scan under DIR for subdirectories (see --subdir-glob).",
    )
    parser.add_argument(
        "--subdir-glob",
        default="SSMaligned*",
        help="With --dir: process subdirs matching this glob (default: SSMaligned*).",
    )
    parser.add_argument(
        "--trim-source",
        choices=TRIM_SOURCE_CHOICES,
        default=TRIM_SOURCE_CORE,
        help=(
            "Residue set for trim: core (strict intersection), "
            "structural_msa (union of SSM columns), continuous_msa (longest "
            "contiguous SSM-matched segment per structure), or mafft (MAFFT columns)."
        ),
    )
    parser.add_argument(
        "--anchor",
        default=None,
        help="Coordinate anchor structure (default: auto, max complete core columns).",
    )
    parser.add_argument(
        "--skip-existing",
        action="store_true",
        help="Skip when core_ssm/coot_log_core.txt already exists.",
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="Overwrite existing aligned output directory (requires re-run).",
    )
    parser.add_argument(
        "--skip-full-ssm",
        action="store_true",
        help="Run phase 0 and A only; skip full native → aligned core (phase B).",
    )
    parser.add_argument(
        "--interactive",
        action="store_true",
        help="Keep Coot open after alignments (default: batch exit).",
    )
    add_log_args(parser)
    args = parser.parse_args()

    if args.scan_dir and args.directory:
        parser.error("Use either a single directory or --dir, not both.")
    if not args.scan_dir and not args.directory:
        parser.error("Provide a directory or --dir=DIR.")

    not_interactive = not args.interactive

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
        if process_ssm_directory(
            target,
            anchor_override=args.anchor,
            skip_existing=args.skip_existing,
            force=args.force,
            not_interactive=not_interactive,
            skip_full_ssm=args.skip_full_ssm,
            trim_source=args.trim_source,
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
