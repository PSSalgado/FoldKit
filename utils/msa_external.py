"""
External MSA (MAFFT) and structure sequence extraction for SSM output directories.
"""

from __future__ import annotations

import os
import shutil
import subprocess
import sys
import tempfile

from utils.ssm_log_parse import (
    SSMSuperpositionBlock,
    _natural_sort_key,
    format_clustal_like,
    format_fasta,
    parse_log_chain_ids,
    parse_log_reference_path,
)

try:
    from Bio.PDB.MMCIFParser import MMCIFParser
    from Bio.PDB.PDBParser import PDBParser
    from Bio.PDB.Polypeptide import protein_letters_3to1

    BIOPYTHON_AVAILABLE = True
except ImportError:
    BIOPYTHON_AVAILABLE = False

_ALIGNED_TAGS = ("_SSMaligned2_", "_SSMaligned_")

# Prefix layout typical of RPM-style MAFFT installs (usr/bin, usr/libexec/mafft).
_MAFFT_LOCAL_SUBPATHS = (
    ("usr/bin/mafft", "usr/libexec/mafft"),
)


def _mafft_candidate_roots() -> list[str]:
    """Search order for a prefix-style MAFFT install (``MAFFT_ROOT``, then common defaults)."""
    roots: list[str] = []
    seen: set[str] = set()

    def _add(path: str | None) -> None:
        if not path:
            return
        expanded = os.path.abspath(os.path.expanduser(path))
        if expanded not in seen:
            seen.add(expanded)
            roots.append(expanded)

    _add(os.environ.get("MAFFT_ROOT"))
    _add(os.path.join(os.path.expanduser("~"), "opt", "xtal", "mafft"))
    _add(os.environ.get("FOLDKIT_XTAL"))
    if os.environ.get("FOLDKIT_XTAL"):
        _add(os.path.join(os.environ["FOLDKIT_XTAL"], "mafft"))
    return roots


def _local_mafft_at_root(root: str) -> tuple[str, str] | None:
    """Return (mafft_bin, mafft_binaries_dir) for a prefix install under root."""
    for bin_rel, libexec_rel in _MAFFT_LOCAL_SUBPATHS:
        bin_path = os.path.join(root, bin_rel)
        libexec = os.path.join(root, libexec_rel)
        if os.path.isfile(bin_path) and os.path.isdir(libexec):
            return bin_path, libexec
    return None


def resolve_mafft_installation() -> tuple[str | None, dict[str, str]]:
    """
    Locate MAFFT and return (executable path, env overrides for subprocess).

    Supports:
    - ``MAFFT_ROOT`` prefix installs (sets ``MAFFT_BINARIES``)
    - ``mafft`` already on ``PATH`` with ``MAFFT_BINARIES`` set (e.g. after ``source env.sh``)
    - system ``mafft`` on ``PATH`` (no overrides)
    """
    env_patch: dict[str, str] = {}

    for root in _mafft_candidate_roots():
        local = _local_mafft_at_root(root)
        if local:
            bin_path, libexec = local
            env_patch = {
                "MAFFT_ROOT": root,
                "MAFFT_BINARIES": libexec,
            }
            return bin_path, env_patch

    on_path = shutil.which("mafft")
    if on_path:
        if os.environ.get("MAFFT_BINARIES"):
            return on_path, {}
        # System install: helpers live under /usr/libexec/mafft
        system_libexec = "/usr/libexec/mafft"
        if os.path.isdir(system_libexec):
            return on_path, {}
        return on_path, {}

    return None, {}


def mafft_subprocess_env() -> dict[str, str]:
    """Full environment for a MAFFT subprocess (includes ``MAFFT_BINARIES`` when needed)."""
    _bin, patch = resolve_mafft_installation()
    env = os.environ.copy()
    if patch:
        env.update(patch)
        bin_dir = os.path.dirname(_bin) if _bin else ""
        if bin_dir:
            env["PATH"] = bin_dir + os.pathsep + env.get("PATH", "")
    return env


def mafft_executable() -> str | None:
    """Path to the MAFFT driver script/binary, or None if not found."""
    path, _patch = resolve_mafft_installation()
    return path


def _label_from_pdb_stem(stem: str, known_labels: set[str]) -> str | None:
    if stem in known_labels:
        return stem
    for tag in _ALIGNED_TAGS:
        if tag in stem:
            model_part = stem.split(tag)[0]
            if model_part in known_labels:
                return model_part
            ref_part = stem.split(tag, 1)[-1]
            if ref_part in known_labels:
                return ref_part
    return None


def discover_label_pdb_map(
    directory: str,
    labels: set[str],
) -> dict[str, str]:
    """Map Coot structure labels to PDB/CIF paths in an SSM output directory."""
    directory = os.path.abspath(os.path.expanduser(directory))
    mapping: dict[str, str] = {}
    if not os.path.isdir(directory):
        return mapping

    for name in sorted(os.listdir(directory)):
        lower = name.lower()
        if not (lower.endswith(".pdb") or lower.endswith(".cif")):
            continue
        stem = os.path.splitext(name)[0]
        label = _label_from_pdb_stem(stem, labels)
        if label and label not in mapping:
            mapping[label] = os.path.join(directory, name)
    return mapping


def _structure_parser():
    if not BIOPYTHON_AVAILABLE:
        return None, None
    return PDBParser(QUIET=True), MMCIFParser(QUIET=True)


def extract_chain_sequence(
    structure_path: str,
    chain_id: str | None = None,
) -> tuple[str, dict[tuple[str, int, str], str]]:
    """
    Return (one-letter sequence, (chain, resnum, icode) -> aa) for one chain.

    Residues are ordered by (chain, resnum, icode). Non-standard residues are skipped.
    """
    if not BIOPYTHON_AVAILABLE:
        return "", {}

    structure_path = os.path.abspath(os.path.expanduser(structure_path))
    pdb_parser, cif_parser = _structure_parser()
    label = os.path.splitext(os.path.basename(structure_path))[0]
    try:
        if structure_path.lower().endswith(".cif"):
            structure = cif_parser.get_structure(label, structure_path)
        else:
            structure = pdb_parser.get_structure(label, structure_path)
    except Exception:
        return "", {}

    model = structure[0]
    chains = list(model.get_chains())
    if not chains:
        return "", {}

    if chain_id:
        chain = model[chain_id]
    else:
        chain = chains[0]

    ordered: list[tuple[tuple[str, int, str], str]] = []
    for residue in chain:
        if residue.id[0] != " ":
            continue
        resnum = residue.id[1]
        icode = residue.id[2] if len(residue.id) > 2 else " "
        icode = icode if icode else " "
        try:
            aa = protein_letters_3to1[residue.resname]
        except KeyError:
            continue
        key = (chain.id, resnum, icode)
        ordered.append((key, aa.upper()))

    residue_map = {k: v for k, v in ordered}
    sequence = "".join(aa for _k, aa in ordered)
    return sequence, residue_map


def extract_chain_residue_keys(
    structure_path: str,
    chain_id: str | None = None,
) -> list[tuple]:
    """
    Return ordered (ResidueKey, one-letter aa) for one chain.

    ResidueKey is imported lazily to avoid circular imports at module load.
    """
    from utils.ssm_log_parse import ResidueKey

    if not BIOPYTHON_AVAILABLE:
        return []

    structure_path = os.path.abspath(os.path.expanduser(structure_path))
    pdb_parser, cif_parser = _structure_parser()
    label = os.path.splitext(os.path.basename(structure_path))[0]
    try:
        if structure_path.lower().endswith(".cif"):
            structure = cif_parser.get_structure(label, structure_path)
        else:
            structure = pdb_parser.get_structure(label, structure_path)
    except Exception:
        return []

    model = structure[0]
    chains = list(model.get_chains())
    if not chains:
        return []

    chain = model[chain_id] if chain_id else chains[0]

    ordered: list[tuple] = []
    for residue in chain:
        if residue.id[0] != " ":
            continue
        resnum = residue.id[1]
        icode = residue.id[2] if len(residue.id) > 2 else " "
        icode = icode if icode else " "
        try:
            aa = protein_letters_3to1[residue.resname]
        except KeyError:
            continue
        key = ResidueKey(chain.id, resnum, icode)
        ordered.append((key, aa.upper()))
    return ordered


def build_label_residue_key_lists(
    directory: str,
    labels: list[str],
    blocks: list,
    log_file: str | None = None,
) -> dict[str, list[tuple]]:
    """Map each structure label to ordered native (ResidueKey, aa) pairs for MAFFT mapping."""
    from utils.ssm_log_parse import parse_log_reference_label, parse_log_reference_path

    label_set = set(labels)
    pdb_map = discover_label_pdb_map(directory, label_set)
    ref_path = parse_log_reference_path(log_file) if log_file else None
    ref_label = parse_log_reference_label(log_file) if log_file else None
    ref_chain, model_chain = parse_log_chain_ids(log_file) if log_file else (None, None)

    out: dict[str, list[tuple]] = {}
    for label in labels:
        path = pdb_map.get(label)
        if not path and ref_path and ref_label == label and os.path.isfile(ref_path):
            path = ref_path
        if path and os.path.isfile(path):
            chain = ref_chain if label == ref_label else model_chain
            keys = extract_chain_residue_keys(path, chain_id=chain)
            if keys:
                out[label] = keys
    return out


def mafft_alignment_to_columns(
    aligned_records: list[tuple[str, str]],
    label_residue_keys: dict[str, list[tuple]],
) -> list[dict]:
    """
    Map MAFFT-aligned sequences to sparse equivalence columns (ResidueKey per label).

    Non-gap alignment columns consume the next native residue for each sequence.
    """
    from utils.ssm_log_parse import ResidueKey

    if not aligned_records:
        return []
    width = len(aligned_records[0][1])
    indices: dict[str, int] = {label: 0 for label in label_residue_keys}
    columns: list[dict] = []
    for col_i in range(width):
        mapping: dict = {}
        for label, seq in aligned_records:
            keys = label_residue_keys.get(label)
            if not keys:
                continue
            aa = seq[col_i]
            if aa == "-":
                continue
            idx = indices[label]
            if idx >= len(keys):
                continue
            key, _native_aa = keys[idx]
            indices[label] = idx + 1
            mapping[label] = key
        if mapping:
            columns.append(mapping)
    return columns


def sequence_from_blocks_fallback(
    label: str,
    blocks: list[SSMSuperpositionBlock],
) -> str:
    """Reconstruct longest ungapped sequence from Target:/Moving: alignment lines."""
    parts: list[str] = []
    for block in blocks:
        if block.reference_label == label and block.target_seq_aligned:
            seq = "".join(c for c in block.target_seq_aligned if c != "-")
            parts.append(seq)
        if block.moving_label == label and block.moving_seq_aligned:
            seq = "".join(c for c in block.moving_seq_aligned if c != "-")
            parts.append(seq)
    if not parts:
        return ""
    return max(parts, key=len)


def build_structure_fasta_records(
    directory: str,
    labels: list[str],
    blocks: list[SSMSuperpositionBlock],
    log_file: str | None = None,
) -> list[tuple[str, str]]:
    """Build one FASTA record per structure label for MAFFT input."""
    from utils.ssm_log_parse import parse_log_reference_label

    label_set = set(labels)
    pdb_map = discover_label_pdb_map(directory, label_set)
    ref_path = parse_log_reference_path(log_file) if log_file else None
    ref_label = parse_log_reference_label(log_file) if log_file else None
    ref_chain, model_chain = parse_log_chain_ids(log_file) if log_file else (None, None)

    records: list[tuple[str, str]] = []
    for label in labels:
        seq = ""
        path = pdb_map.get(label)
        if not path and ref_path and ref_label == label and os.path.isfile(ref_path):
            path = ref_path

        if path and os.path.isfile(path):
            chain = ref_chain if label == ref_label else model_chain
            seq, _res_map = extract_chain_sequence(path, chain_id=chain)
        if not seq:
            seq = sequence_from_blocks_fallback(label, blocks)
        if seq:
            records.append((label, seq))
    return records


def mafft_available() -> bool:
    return mafft_executable() is not None


def run_mafft(
    in_fasta: str,
    out_fasta: str,
) -> tuple[bool, str]:
    """
    Run MAFFT on in_fasta and write aligned FASTA to out_fasta.

    Returns (success, message).
    """
    mafft_bin = mafft_executable()
    if not mafft_bin:
        return False, (
            "MAFFT not found (set MAFFT_ROOT to a prefix install, ensure mafft "
            "is on PATH, or use --no-mafft)"
        )

    in_fasta = os.path.abspath(in_fasta)
    out_fasta = os.path.abspath(out_fasta)
    os.makedirs(os.path.dirname(out_fasta) or ".", exist_ok=True)
    env = mafft_subprocess_env()

    try:
        with open(out_fasta, "w", encoding="utf-8") as out_fh:
            result = subprocess.run(
                [mafft_bin, "--auto", in_fasta],
                stdout=out_fh,
                stderr=subprocess.PIPE,
                text=True,
                check=False,
                env=env,
            )
    except OSError as exc:
        return False, "MAFFT failed: {}".format(exc)

    if result.returncode != 0:
        err = (result.stderr or "").strip()
        return False, "MAFFT exit {}: {}".format(result.returncode, err or "unknown error")

    if not os.path.isfile(out_fasta) or os.path.getsize(out_fasta) == 0:
        return False, "MAFFT produced no output"

    return True, "MAFFT alignment written to {}".format(out_fasta)


def _read_fasta(path: str) -> list[tuple[str, str]]:
    records: list[tuple[str, str]] = []
    name: str | None = None
    parts: list[str] = []
    with open(path, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            line = line.rstrip("\n\r")
            if line.startswith(">"):
                if name is not None:
                    records.append((name, "".join(parts)))
                name = line[1:].split()[0]
                parts = []
            elif name is not None:
                parts.append(line.strip())
        if name is not None:
            records.append((name, "".join(parts)))
    return records


def write_mafft_msa(
    directory: str,
    labels: list[str],
    blocks: list[SSMSuperpositionBlock],
    log_file: str | None,
    *,
    out_fasta: str,
    out_aln: str,
    fallback_records: list[tuple[str, str]] | None = None,
) -> tuple[bool, str]:
    """
    Run MAFFT on structure sequences; on failure write fallback_records to out_fasta.
    """
    records = build_structure_fasta_records(directory, labels, blocks, log_file)
    if len(records) < 2:
        msg = "Need at least two sequences for MAFFT; found {}".format(len(records))
        if fallback_records:
            with open(out_fasta, "w", encoding="utf-8") as fh:
                fh.write(format_fasta(fallback_records))
            labels_fb = [n for n, _ in fallback_records]
            seqs_fb = [s for _, s in fallback_records]
            with open(out_aln, "w", encoding="utf-8") as fh:
                fh.write(
                    format_clustal_like(
                        labels_fb,
                        seqs_fb,
                        title="SSM structural MSA (MAFFT unavailable)",
                    )
                )
            return False, msg + "; wrote structural fallback"
        return False, msg

    with tempfile.NamedTemporaryFile(
        mode="w",
        suffix=".fasta",
        delete=False,
        encoding="utf-8",
    ) as tmp:
        tmp.write(format_fasta(records))
        tmp_path = tmp.name

    try:
        ok, msg = run_mafft(tmp_path, out_fasta)
    finally:
        try:
            os.unlink(tmp_path)
        except OSError:
            pass

    if ok:
        aligned = _read_fasta(out_fasta)
        if aligned:
            aln_labels = [n for n, _ in aligned]
            aln_seqs = [s for _, s in aligned]
            with open(out_aln, "w", encoding="utf-8") as fh:
                fh.write(
                    format_clustal_like(
                        aln_labels,
                        aln_seqs,
                        title="SSM MAFFT alignment ({} structures)".format(len(aligned)),
                    )
                )
        return True, msg

    print("WARNING: {}".format(msg), file=sys.stderr)
    if fallback_records:
        with open(out_fasta, "w", encoding="utf-8") as fh:
            fh.write(format_fasta(fallback_records))
        labels_fb = [n for n, _ in fallback_records]
        seqs_fb = [s for _, s in fallback_records]
        with open(out_aln, "w", encoding="utf-8") as fh:
            fh.write(
                format_clustal_like(
                    labels_fb,
                    seqs_fb,
                    title="SSM structural MSA (MAFFT fallback)",
                )
            )
        return False, msg + "; wrote structural fallback"
    return False, msg
