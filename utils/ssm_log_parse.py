"""
Parse Coot SSM superposition output from coot_log*.txt files.

Coot prints per-pair structural correspondence in a ``Moving Reference Distance(/A)``
table (chain/residue pairs with Cα distance) and optional horizontal one-letter
``Moving:`` / ``Target:`` sequence alignment lines before the ``INFO: core rmsd`` block.
"""

from __future__ import annotations

import os
import re
from dataclasses import dataclass, field
from typing import Iterable


def _natural_sort_key(text: str) -> list:
    return [int(part) if part.isdigit() else part.lower() for part in re.split(r"(\d+)", text)]

# One-to-many: "Superposing model onto ref"
# All-vs-all: "  Superposing model onto ref"
_RE_SUPERPOSING = re.compile(
    r"^\s*Superposing\s+(\S+)\s+onto\s+(\S+)\s*$"
)
_RE_DIST_HEADER = re.compile(r"Moving\s+Reference\s+Distance\s*\(/A\)", re.I)
_RE_DIST_PAIR = re.compile(
    r"^\s+(\S+)\s+(\d+)([^\s<]*)\s+<--->\s+(\S+)\s+(\d+)([^\s<:]*)\s+:\s+([\d.]+)\s*$"
)
_RE_MOVING_SEQ = re.compile(r"^\s*Moving:\s+(.*)\s*$")
_RE_TARGET_SEQ = re.compile(r"^\s*Target:\s+(.*)\s*$")
_RE_CORE_RMSD = re.compile(r"INFO:\s*core rmsd")
_RE_STAT = re.compile(
    r"^\s*(number of residues in reference structure|"
    r"number of residues in moving structure|"
    r"number of residues in aligned sections \(reference\)|"
    r"number of residues in aligned sections \(moving\)|"
    r"number of aligned residues|number of gaps|"
    r"number of misdirections|number of SSE combinations|"
    r"sequence identity)\s*:\s*(.+)\s*$"
)
_RE_LOG_REFERENCE = re.compile(
    r"^\s*#\s*Reference:\s+(\S+?)(?:\s+\(chain\s+\S+\))?\s*$"
)
_RE_LOG_REF_CHAIN = re.compile(r"^\s*#\s*Ref chain:\s+(\S+)", re.I)
_RE_LOG_MODEL_CHAIN = re.compile(r"^\s*#\s*Model chain:\s+(\S+)", re.I)


@dataclass(frozen=True)
class ResidueKey:
    chain: str
    resnum: int
    icode: str = " "

    def as_tuple(self) -> tuple[str, int, str]:
        ic = self.icode if self.icode else " "
        return (self.chain, self.resnum, ic)


@dataclass(frozen=True)
class AlignmentMappingEntry:
    """One structurally aligned column mapped from sequence lines and equivalences."""

    column_index: int
    moving_aa: str
    target_aa: str
    moving_key: ResidueKey
    reference_key: ResidueKey
    distance: float | None = None


@dataclass
class SSMSuperpositionBlock:
    moving_label: str
    reference_label: str
    alignment_line: str = ""
    equivalences: list[tuple[ResidueKey, ResidueKey, float | None]] = field(
        default_factory=list
    )
    moving_seq_aligned: str = ""
    target_seq_aligned: str = ""
    stats: dict[str, str] = field(default_factory=dict)


def _normalize_icode(raw: str) -> str:
    s = (raw or "").strip()
    if not s:
        return " "
    return s[0]


def _parse_residue_key(chain: str, resnum: str, icode_raw: str) -> ResidueKey:
    return ResidueKey(
        chain=chain.strip(),
        resnum=int(resnum),
        icode=_normalize_icode(icode_raw),
    )


def parse_superposition_labels(line: str) -> tuple[str, str] | None:
    m = _RE_SUPERPOSING.match(line.strip())
    if not m:
        return None
    return m.group(1), m.group(2)


def parse_ssm_log(log_file: str) -> list[SSMSuperpositionBlock]:
    """Return one block per SSM superposition found in a Coot log."""
    log_file = os.path.abspath(os.path.expanduser(log_file))
    blocks: list[SSMSuperpositionBlock] = []
    current: SSMSuperpositionBlock | None = None
    in_dist_table = False
    moving_seq_parts: list[str] = []
    target_seq_parts: list[str] = []

    def _flush() -> None:
        nonlocal current, in_dist_table, moving_seq_parts, target_seq_parts
        if current is None:
            return
        current.moving_seq_aligned = "".join(moving_seq_parts)
        current.target_seq_aligned = "".join(target_seq_parts)
        blocks.append(current)
        current = None
        in_dist_table = False
        moving_seq_parts = []
        target_seq_parts = []

    with open(log_file, encoding="utf-8", errors="replace") as fh:
        for raw in fh:
            line = raw.rstrip("\n\r")
            labels = parse_superposition_labels(line)
            if labels is not None:
                _flush()
                mov, ref = labels
                current = SSMSuperpositionBlock(
                    moving_label=mov,
                    reference_label=ref,
                    alignment_line=line.strip(),
                )
                continue

            if current is None:
                continue

            if _RE_DIST_HEADER.search(line):
                in_dist_table = True
                continue

            if in_dist_table:
                m = _RE_DIST_PAIR.match(line)
                if m:
                    mov_key = _parse_residue_key(m.group(1), m.group(2), m.group(3))
                    ref_key = _parse_residue_key(m.group(4), m.group(5), m.group(6))
                    dist = float(m.group(7))
                    current.equivalences.append((mov_key, ref_key, dist))
                    continue
                if _RE_MOVING_SEQ.match(line) or _RE_TARGET_SEQ.match(line):
                    in_dist_table = False
                elif line.strip() == "":
                    continue
                elif not line.startswith(" "):
                    in_dist_table = False

            m_mov = _RE_MOVING_SEQ.match(line)
            if m_mov:
                moving_seq_parts.append(m_mov.group(1))
                continue

            m_tgt = _RE_TARGET_SEQ.match(line)
            if m_tgt:
                target_seq_parts.append(m_tgt.group(1))
                continue

            if _RE_CORE_RMSD.search(line):
                current.stats["core rmsd line"] = line.strip()
                continue

            m_stat = _RE_STAT.match(line)
            if m_stat:
                current.stats[m_stat.group(1)] = m_stat.group(2).strip()
                continue

    _flush()
    return blocks


def canonical_pair_label(a: str, b: str) -> tuple[str, str]:
    """Undirected pair key in natural sort order."""
    ka, kb = _natural_sort_key(a), _natural_sort_key(b)
    return (a, b) if ka <= kb else (b, a)


def residue_node(structure: str, key: ResidueKey) -> tuple[str, tuple[str, int, str]]:
    return (structure, key.as_tuple())


class UnionFind:
    def __init__(self) -> None:
        self.parent: dict[tuple, tuple] = {}

    def find(self, x: tuple) -> tuple:
        if x not in self.parent:
            self.parent[x] = x
        while self.parent[x] != x:
            self.parent[x] = self.parent[self.parent[x]]
            x = self.parent[x]
        return x

    def union(self, a: tuple, b: tuple) -> None:
        ra, rb = self.find(a), self.find(b)
        if ra != rb:
            self.parent[rb] = ra


def canonical_edge(
    n1: tuple[str, tuple[str, int, str]],
    n2: tuple[str, tuple[str, int, str]],
) -> tuple[tuple[str, tuple[str, int, str]], tuple[str, tuple[str, int, str]]]:
    a = (n1[0], n1[1])
    b = (n2[0], n2[1])
    ka = (_natural_sort_key(a[0]), a[1])
    kb = (_natural_sort_key(b[0]), b[1])
    return (a, b) if ka <= kb else (b, a)


def build_pairwise_edge_set(
    blocks: Iterable[SSMSuperpositionBlock],
) -> set[tuple[tuple[str, tuple], tuple[str, tuple]]]:
    edges: set[tuple[tuple[str, tuple], tuple[str, tuple]]] = set()
    for block in blocks:
        for mov_key, ref_key, _dist in block.equivalences:
            n1 = (block.moving_label, mov_key.as_tuple())
            n2 = (block.reference_label, ref_key.as_tuple())
            edges.add(canonical_edge(n1, n2))
    return edges


def _residue_key_sort(key: ResidueKey) -> tuple:
    return (_natural_sort_key(key.chain), key.resnum, key.icode)


def map_equivalences_to_alignment(
    block: SSMSuperpositionBlock,
) -> list[AlignmentMappingEntry]:
    """
    Map distance-table equivalences to paired non-gap columns in Moving:/Target: lines.

    Returns one entry per structural equivalence row, in table order.
    """
    entries: list[AlignmentMappingEntry] = []
    mov_seq = block.moving_seq_aligned
    ref_seq = block.target_seq_aligned
    if not mov_seq or not ref_seq or len(mov_seq) != len(ref_seq):
        return entries

    eq_rows = list(block.equivalences)
    eq_i = 0
    for col_i, (mov_aa, ref_aa) in enumerate(zip(mov_seq, ref_seq)):
        if mov_aa == "-" or ref_aa == "-":
            continue
        if eq_i >= len(eq_rows):
            break
        mov_key, ref_key, dist = eq_rows[eq_i]
        entries.append(
            AlignmentMappingEntry(
                column_index=col_i,
                moving_aa=mov_aa,
                target_aa=ref_aa,
                moving_key=mov_key,
                reference_key=ref_key,
                distance=dist,
            )
        )
        eq_i += 1
    return entries


def validate_equivalence_sequence_mapping(
    block: SSMSuperpositionBlock,
) -> list[str]:
    """Return warning messages when equivalences and sequence lines disagree."""
    warnings: list[str] = []
    pair = "{} vs {}".format(block.moving_label, block.reference_label)
    mov_seq = block.moving_seq_aligned
    ref_seq = block.target_seq_aligned
    n_equiv = len(block.equivalences)

    if not mov_seq or not ref_seq:
        if n_equiv:
            warnings.append(
                "{}: {} equivalence rows but no Moving:/Target: sequence lines".format(
                    pair, n_equiv
                )
            )
        return warnings

    if len(mov_seq) != len(ref_seq):
        warnings.append(
            "{}: Moving and Target alignment lengths differ ({} vs {})".format(
                pair, len(mov_seq), len(ref_seq)
            )
        )
        return warnings

    paired_cols = sum(
        1
        for mov_aa, ref_aa in zip(mov_seq, ref_seq)
        if mov_aa != "-" and ref_aa != "-"
    )
    if paired_cols < n_equiv:
        warnings.append(
            "{}: {} equivalence rows but only {} paired alignment columns".format(
                pair, n_equiv, paired_cols
            )
        )
    elif paired_cols > n_equiv:
        # Coot often prints extra terminal pairs in Moving:/Target: lines that are
        # not listed in the distance table; TSV output uses the table only.
        pass

    mapped = map_equivalences_to_alignment(block)
    if len(mapped) < n_equiv:
        warnings.append(
            "{}: only {} of {} equivalences mapped to sequence columns".format(
                pair, len(mapped), n_equiv
            )
        )

    return warnings


def parse_log_reference_path(log_file: str) -> str | None:
    """Return the reference structure path from a FoldKit SSM log header, if present."""
    log_file = os.path.abspath(os.path.expanduser(log_file))
    try:
        with open(log_file, encoding="utf-8", errors="replace") as fh:
            for _ in range(40):
                line = fh.readline()
                if not line:
                    break
                m = _RE_LOG_REFERENCE.match(line.rstrip("\n\r"))
                if m:
                    return os.path.expanduser(m.group(1))
    except OSError:
        return None
    return None


def parse_log_reference_label(log_file: str) -> str | None:
    """Return the reference structure stem from a FoldKit SSM log header, if present."""
    ref_path = parse_log_reference_path(log_file)
    if not ref_path:
        return None
    return os.path.splitext(os.path.basename(ref_path))[0]


def parse_log_chain_ids(log_file: str) -> tuple[str | None, str | None]:
    """Return (ref_chain, model_chain) from log header lines, if recorded."""
    log_file = os.path.abspath(os.path.expanduser(log_file))
    ref_chain: str | None = None
    model_chain: str | None = None
    try:
        with open(log_file, encoding="utf-8", errors="replace") as fh:
            for _ in range(40):
                line = fh.readline()
                if not line:
                    break
                m_ref = _RE_LOG_REF_CHAIN.match(line)
                if m_ref:
                    ref_chain = m_ref.group(1)
                m_mod = _RE_LOG_MODEL_CHAIN.match(line)
                if m_mod:
                    model_chain = m_mod.group(1)
    except OSError:
        pass
    return ref_chain, model_chain


def infer_primary_reference(
    blocks: Iterable[SSMSuperpositionBlock],
    log_file: str | None = None,
) -> str | None:
    """
    Return the primary reference label for one-to-many SSM logs.

    Prefers the ``# Reference:`` log header, then a single shared reference_label
    across all blocks. Returns None for all-vs-all (multiple references).
    """
    if log_file:
        header_label = parse_log_reference_label(log_file)
        if header_label:
            return header_label

    blocks_list = list(blocks)
    if not blocks_list:
        return None

    ref_labels = {b.reference_label for b in blocks_list}
    if len(ref_labels) == 1:
        return next(iter(ref_labels))

    mov_labels = {b.moving_label for b in blocks_list}
    # One-to-many: one reference, many moving structures (reference never appears as moving).
    if len(ref_labels) == 1 and len(mov_labels & ref_labels) == 0:
        return next(iter(ref_labels))

    return None


def build_reference_anchored_core(
    reference_label: str,
    blocks: Iterable[SSMSuperpositionBlock],
) -> tuple[list[str], list[dict[str, ResidueKey]]]:
    """
    Core columns for one-to-many SSM: reference residues matched in every model block.

    Returns (labels, core_columns) where labels = [reference] + sorted moving models.
    """
    ref_blocks = [b for b in blocks if b.reference_label == reference_label]
    moving_labels = sorted(
        {b.moving_label for b in ref_blocks},
        key=_natural_sort_key,
    )
    labels = [reference_label] + moving_labels
    if not ref_blocks or not moving_labels:
        return labels, []

    ref_sets: list[set[tuple[str, int, str]]] = []
    block_maps: list[tuple[str, dict[tuple[str, int, str], ResidueKey]]] = []
    for block in ref_blocks:
        ref_to_mov: dict[tuple[str, int, str], ResidueKey] = {}
        for mov_key, ref_key, _dist in block.equivalences:
            ref_to_mov[ref_key.as_tuple()] = mov_key
        ref_sets.append(set(ref_to_mov))
        block_maps.append((block.moving_label, ref_to_mov))

    core_ref_tuples = set.intersection(*ref_sets) if ref_sets else set()
    core_columns: list[dict[str, ResidueKey]] = []
    for ref_tuple in sorted(core_ref_tuples, key=lambda t: (_natural_sort_key(t[0]), t[1], t[2])):
        ref_key = ResidueKey(ref_tuple[0], ref_tuple[1], ref_tuple[2])
        column = {reference_label: ref_key}
        for mov_label, ref_to_mov in block_maps:
            column[mov_label] = ref_to_mov[ref_tuple]
        core_columns.append(column)

    return labels, core_columns


def build_reference_anchored_msa_columns(
    reference_label: str,
    blocks: Iterable[SSMSuperpositionBlock],
) -> tuple[list[str], list[dict[str, ResidueKey]]]:
    """
    Reference-indexed structural columns: union of reference residues, gaps where unaligned.
    """
    ref_blocks = [b for b in blocks if b.reference_label == reference_label]
    moving_labels = sorted(
        {b.moving_label for b in ref_blocks},
        key=_natural_sort_key,
    )
    labels = [reference_label] + moving_labels
    if not ref_blocks:
        return labels, []

    all_ref_keys: set[tuple[str, int, str]] = set()
    per_block: list[tuple[str, dict[tuple[str, int, str], ResidueKey]]] = []
    for block in ref_blocks:
        ref_to_mov: dict[tuple[str, int, str], ResidueKey] = {}
        for mov_key, ref_key, _dist in block.equivalences:
            rt = ref_key.as_tuple()
            all_ref_keys.add(rt)
            ref_to_mov[rt] = mov_key
        per_block.append((block.moving_label, ref_to_mov))

    columns: list[dict[str, ResidueKey]] = []
    for ref_tuple in sorted(all_ref_keys, key=lambda t: (_natural_sort_key(t[0]), t[1], t[2])):
        ref_key = ResidueKey(ref_tuple[0], ref_tuple[1], ref_tuple[2])
        column: dict[str, ResidueKey] = {reference_label: ref_key}
        for mov_label, ref_to_mov in per_block:
            if ref_tuple in ref_to_mov:
                column[mov_label] = ref_to_mov[ref_tuple]
        columns.append(column)

    return labels, columns


def filter_complete_pairwise_core(
    labels: list[str],
    columns: list[dict[str, ResidueKey]],
    edge_set: set[tuple[tuple[str, tuple], tuple[str, tuple]]],
) -> list[dict[str, ResidueKey]]:
    """
    Keep columns where every structure is present and each pair is directly
    aligned in the corresponding SSM superposition (not merely transitive).
    """
    core: list[dict[str, ResidueKey]] = []
    label_set = set(labels)
    for mapping in columns:
        if set(mapping.keys()) != label_set or len(mapping) != len(label_set):
            continue
        complete = True
        for i, si in enumerate(labels):
            for sj in labels[i + 1 :]:
                n1 = (si, mapping[si].as_tuple())
                n2 = (sj, mapping[sj].as_tuple())
                if canonical_edge(n1, n2) not in edge_set:
                    complete = False
                    break
            if not complete:
                break
        if complete:
            core.append(mapping)
    return core


def build_equivalence_components(
    blocks: Iterable[SSMSuperpositionBlock],
    *,
    log_file: str | None = None,
) -> tuple[
    list[str],
    list[dict[str, ResidueKey]],
    list[dict[str, ResidueKey]],
    dict[tuple[str, str], list],
    str | None,
]:
    """
    Merge pairwise SSM equivalences into connected components (MSA columns).

    Returns (structure_labels, all_columns, core_columns, pairwise_equivs, primary_reference).
    all_columns: components with >=2 structures (union of pairwise alignments).
    core_columns: reference-anchored intersection (one-to-many) or complete pairwise core.
    """
    blocks_list = list(blocks)
    primary_reference = infer_primary_reference(blocks_list, log_file)
    structures: set[str] = set()
    pairwise: dict[tuple[str, str], list[tuple[ResidueKey, ResidueKey]]] = {}
    uf = UnionFind()

    for block in blocks_list:
        structures.add(block.moving_label)
        structures.add(block.reference_label)
        pair = canonical_pair_label(block.moving_label, block.reference_label)
        pair_equivs: list[tuple[ResidueKey, ResidueKey]] = []
        for mov_key, ref_key, _dist in block.equivalences:
            pair_equivs.append((mov_key, ref_key))
            uf.union(
                residue_node(block.moving_label, mov_key),
                residue_node(block.reference_label, ref_key),
            )
        if pair_equivs:
            pairwise.setdefault(pair, []).extend(pair_equivs)

    labels = sorted(structures, key=_natural_sort_key)

    components: dict[tuple, dict[str, ResidueKey]] = {}
    for node in uf.parent:
        root = uf.find(node)
        struct, _res_tuple = node
        _chain, _resnum, _icode = node[1]
        key = ResidueKey(_chain, _resnum, _icode)
        bucket = components.setdefault(root, {})
        if struct in bucket:
            # Inconsistent merge; keep first mapping.
            continue
        bucket[struct] = key

    core_columns: list[dict[str, ResidueKey]] = []
    edge_set = build_pairwise_edge_set(blocks_list)
    all_multi: list[dict[str, ResidueKey]] = []
    for mapping in components.values():
        if len(mapping) >= 2:
            all_multi.append(mapping)
    all_multi.sort(
        key=lambda m: (
            _natural_sort_key(labels[0]) if labels else "",
            m.get(labels[0], ResidueKey("", 0)).resnum if labels else 0,
            m.get(labels[0], ResidueKey("", 0)).icode if labels else " ",
        )
    )
    if primary_reference:
        _ref_labels, core_columns = build_reference_anchored_core(
            primary_reference, blocks_list
        )
        core_columns.sort(
            key=lambda m: _residue_key_sort(m[primary_reference])
        )
    else:
        core_columns = filter_complete_pairwise_core(labels, all_multi, edge_set)
        core_columns.sort(
            key=lambda m: (
                _natural_sort_key(labels[0]),
                m[labels[0]].resnum,
                m[labels[0]].icode,
            )
        )
    return labels, all_multi, core_columns, pairwise, primary_reference


def one_letter_from_aligned_block(
    block: SSMSuperpositionBlock,
) -> dict[tuple[str, tuple], str]:
    """
    Map (structure_label, residue_tuple) -> one-letter code from Moving:/Target: lines.

    Paired non-gap columns are matched in order to rows of the distance table.
    """
    out: dict[tuple[str, tuple], str] = {}
    mov_seq = block.moving_seq_aligned
    ref_seq = block.target_seq_aligned
    if not mov_seq or not ref_seq or len(mov_seq) != len(ref_seq):
        return out

    eq_keys = [(m, r) for m, r, _ in block.equivalences]
    eq_i = 0
    for mov_aa, ref_aa in zip(mov_seq, ref_seq):
        if mov_aa == "-" or ref_aa == "-":
            continue
        if eq_i >= len(eq_keys):
            break
        mov_key, ref_key = eq_keys[eq_i]
        if mov_aa not in "-^.":
            out[(block.moving_label, mov_key.as_tuple())] = mov_aa.upper()
        if ref_aa not in "-^.":
            out[(block.reference_label, ref_key.as_tuple())] = ref_aa.upper()
        eq_i += 1
    return out


def collect_one_letter_codes(
    blocks: Iterable[SSMSuperpositionBlock],
) -> dict[tuple[str, tuple], str]:
    """Merge one-letter codes from all blocks (later blocks do not overwrite)."""
    codes: dict[tuple[str, tuple], str] = {}
    for block in blocks:
        for key, aa in one_letter_from_aligned_block(block).items():
            codes.setdefault(key, aa)
    return codes


def write_equivalence_tsv(
    path: str,
    moving_label: str,
    reference_label: str,
    equivalences: list[tuple[ResidueKey, ResidueKey, float | None]],
    *,
    header_comment: str | None = None,
) -> None:
    out_dir = os.path.dirname(path) or "."
    os.makedirs(out_dir, exist_ok=True)
    with open(path, "w", encoding="utf-8") as fh:
        fh.write(
            "# SSM structural equivalences: moving={} reference={}\n".format(
                moving_label, reference_label
            )
        )
        if header_comment:
            fh.write("# {}\n".format(header_comment))
        fh.write(
            "moving_chain\tmoving_resnum\tmoving_icode\t"
            "reference_chain\treference_resnum\treference_icode\tdistance_A\n"
        )
        for mov_key, ref_key, dist in equivalences:
            dist_s = "" if dist is None else "{:.4f}".format(dist)
            fh.write(
                "{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(
                    mov_key.chain,
                    mov_key.resnum,
                    mov_key.icode,
                    ref_key.chain,
                    ref_key.resnum,
                    ref_key.icode,
                    dist_s,
                )
            )


def read_equivalence_tsv(path: str) -> list[tuple[ResidueKey, ResidueKey]]:
    equivs: list[tuple[ResidueKey, ResidueKey]] = []
    with open(path, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            if line.lower().startswith("moving_chain"):
                continue
            parts = line.split("\t")
            if len(parts) < 6:
                continue
            mov = ResidueKey(parts[0], int(parts[1]), parts[2])
            ref = ResidueKey(parts[3], int(parts[4]), parts[5])
            equivs.append((mov, ref))
    return equivs


def blocks_from_equivalence_dir(equiv_dir: str) -> list[SSMSuperpositionBlock]:
    """Load pairwise TSV files written by extract_rmsd.py --seq-align."""
    equiv_dir = os.path.abspath(os.path.expanduser(equiv_dir))
    blocks: list[SSMSuperpositionBlock] = []
    if not os.path.isdir(equiv_dir):
        return blocks
    for name in sorted(os.listdir(equiv_dir)):
        if not name.endswith(".tsv"):
            continue
        stem = name[:-4]
        if "_vs_" not in stem:
            continue
        mov, ref = stem.split("_vs_", 1)
        path = os.path.join(equiv_dir, name)
        equivs_raw = read_equivalence_tsv(path)
        equivs = [(m, r, None) for m, r in equivs_raw]
        blocks.append(
            SSMSuperpositionBlock(
                moving_label=mov,
                reference_label=ref,
                alignment_line="{} vs {}".format(mov, ref),
                equivalences=equivs,
            )
        )
    return blocks


def discover_coot_log(directory: str) -> str | None:
    """Return coot_log.txt or the sole coot_log_*.txt in directory, if any."""
    directory = os.path.abspath(os.path.expanduser(directory))
    default = os.path.join(directory, "coot_log.txt")
    if os.path.isfile(default):
        return default
    labeled = sorted(
        p
        for p in (
            os.path.join(directory, n)
            for n in os.listdir(directory)
            if n.startswith("coot_log_") and n.endswith(".txt")
        )
        if os.path.isfile(p)
    )
    if len(labeled) == 1:
        return labeled[0]
    return None


def _count_ssm_blocks_in_log(log_path: str) -> int:
    if not log_path or not os.path.isfile(log_path):
        return 0
    return len(parse_ssm_log(log_path))


def best_coot_log(directory: str) -> str | None:
    """
    Return the Coot log in *directory* with the most parseable SSM blocks.

    Prefer the file that actually contains superposition data over stale crash
    stubs or completion footers when multiple ``coot_log*.txt`` files exist.
    """
    directory = os.path.abspath(os.path.expanduser(directory))
    if not os.path.isdir(directory):
        return None
    candidates: list[str] = []
    for name in os.listdir(directory):
        if name == "coot_log.txt" or (
            name.startswith("coot_log_") and name.endswith(".txt")
        ):
            path = os.path.join(directory, name)
            if os.path.isfile(path):
                candidates.append(path)
    if not candidates:
        return None
    if len(candidates) == 1:
        return candidates[0]

    def sort_key(path: str) -> tuple[int, int, str]:
        return (
            _count_ssm_blocks_in_log(path),
            os.path.getsize(path),
            path,
        )

    return max(candidates, key=sort_key)


def resolve_step1_coot_log(
    ssm_dir: str,
    log_file: str | None = None,
) -> str | None:
    """Pick the best Step 1 SSM log under *ssm_dir* (most blocks wins)."""
    ssm_dir = os.path.abspath(os.path.expanduser(ssm_dir))
    best = best_coot_log(ssm_dir)
    if log_file is None:
        return best
    log_file = os.path.abspath(os.path.expanduser(log_file))
    if not os.path.isfile(log_file):
        return best
    if best and _count_ssm_blocks_in_log(best) > _count_ssm_blocks_in_log(log_file):
        return best
    return log_file


def consolidate_ssm_logs(ssm_dir: str, label: str) -> str | None:
    """
    Copy the best log in *ssm_dir* to ``coot_log_<label>.txt`` (canonical path).

    Returns the canonical log path, or None if no log exists.
    """
    ssm_dir = os.path.abspath(os.path.expanduser(ssm_dir))
    best = best_coot_log(ssm_dir)
    if not best:
        return None
    canonical = os.path.join(ssm_dir, "coot_log_{}.txt".format(label))
    if os.path.abspath(best) != os.path.abspath(canonical):
        with open(best, encoding="utf-8", errors="replace") as src:
            content = src.read()
        with open(canonical, "w", encoding="utf-8") as dst:
            dst.write(content)
    return canonical


def _graph_errors_without_blocks(
    log_file: str,
    blocks: list[SSMSuperpositionBlock],
) -> list[str]:
    block_pairs = {
        (block.moving_label, block.reference_label) for block in blocks
    }
    labels: set[str] = set()
    for mov, ref in block_pairs:
        labels.add(mov)
        labels.add(ref)

    hard: list[str] = []
    for msg in audit_coot_ssm_log(log_file):
        stem = None
        for pattern in (_RE_COOT_CANT_MAKE_GRAPH, _RE_COOT_SSM_SUPERPOSE_FAIL):
            m = pattern.search(msg)
            if m:
                stem = os.path.splitext(os.path.basename(m.group(1)))[0]
                break
        if stem is None:
            hard.append(msg)
            continue
        has_block = any(
            (stem, ref) in block_pairs for ref in labels if ref != stem
        )
        if not has_block:
            hard.append(msg)
    return hard


def block_has_valid_ssm(block: SSMSuperpositionBlock) -> bool:
    """True when a parsed block has structural equivalences or core RMSD."""
    if block.equivalences:
        return True
    if _parse_core_rmsd_angstrom(block) is not None:
        return True
    return bool(block.stats.get("core rmsd line"))


def substantive_ssm_pair_count(log_file: str) -> int:
    """
    Count superposition pairs with real SSM metrics (not just ``Superposing`` lines).

    Falls back to counting ``INFO: core rmsd achieved`` lines in the raw log.
    """
    log_file = os.path.abspath(os.path.expanduser(log_file))
    blocks = parse_ssm_log(log_file)
    count = sum(1 for block in blocks if block_has_valid_ssm(block))
    if count:
        return count
    with open(log_file, encoding="utf-8", errors="replace") as fh:
        text = fh.read()
    return len(re.findall(r"INFO:\s*core rmsd achieved", text, re.I))


def validate_ssm_extract_outputs(
    ssm_dir: str,
    label: str | None = None,
) -> tuple[bool, list[str]]:
    """Check Step 2 SSM extract artefacts; header-only equivalence TSVs are invalid."""
    issues: list[str] = []
    ssm_dir = os.path.abspath(os.path.expanduser(ssm_dir))
    if label is None:
        label = ssm_output_run_label(ssm_dir)

    rmsd_candidates = [
        os.path.join(ssm_dir, "rmsd_SSM_values_{}.txt".format(label)),
        os.path.join(ssm_dir, "rmsd_SSM_values.txt"),
    ]
    rmsd_path = next((p for p in rmsd_candidates if os.path.isfile(p)), None)
    if not rmsd_path:
        issues.append("missing rmsd_SSM_values file")
    else:
        text = open(rmsd_path, encoding="utf-8", errors="replace").read()
        if not re.search(r"core rmsd|INFO:\s*core rmsd", text, re.I):
            issues.append("RMSD file has no core rmsd values")

    equiv_dir = os.path.join(ssm_dir, "ssm_equivalences")
    if not os.path.isdir(equiv_dir):
        issues.append("missing ssm_equivalences/")
    else:
        tsvs = [f for f in os.listdir(equiv_dir) if f.endswith(".tsv")]
        if not tsvs:
            issues.append("ssm_equivalences/ has no TSV files")
        else:
            with_data = 0
            for name in tsvs:
                path = os.path.join(equiv_dir, name)
                with open(path, encoding="utf-8", errors="replace") as fh:
                    for line in fh:
                        line = line.strip()
                        if not line or line.startswith("#"):
                            continue
                        if line.startswith("moving_chain"):
                            continue
                        with_data += 1
                        break
            if with_data == 0:
                issues.append(
                    "all {} equivalence TSVs are header-only".format(len(tsvs))
                )

    return (len(issues) == 0, issues)


def format_fasta(records: list[tuple[str, str]]) -> str:
    lines: list[str] = []
    for name, seq in records:
        lines.append(">{}\n".format(name))
        for i in range(0, len(seq), 70):
            lines.append(seq[i : i + 70] + "\n")
    return "".join(lines)


_SSM_DIR_PREFIXES = (
    "SSMaligned_all_vs_all_",
    "SSMaligned2_",
    "SSMaligned_",
)
_RE_LOG_DIRECTORIES = re.compile(r"^\s*#\s*Directories:\s+(.+)\s*$")


def ssm_output_run_label(ssm_dir: str | os.PathLike) -> str:
    """
    Derive run label for ``SSM_aligned_core_<run_label>/`` from an SSM output directory.

    Example: ``SSMaligned_all_vs_all_my_run`` → ``my_run``.
    """
    ssm_dir = os.path.abspath(os.path.expanduser(str(ssm_dir)))
    name = os.path.basename(ssm_dir.rstrip("/\\"))
    for prefix in _SSM_DIR_PREFIXES:
        if name.startswith(prefix):
            suffix = name[len(prefix) :]
            if suffix:
                return suffix
    log_path = discover_coot_log(ssm_dir)
    if log_path:
        base = os.path.basename(log_path)
        if base.startswith("coot_log_") and base.endswith(".txt"):
            label = base[9:-4]
            if label:
                return label
        if base == "coot_log.txt":
            pass
    return name


def aligned_core_output_dir(ssm_dir: str | os.PathLike) -> str:
    """Sibling output directory ``SSM_aligned_core_<run_label>/``."""
    return aligned_trim_output_dir(ssm_dir, "core")


TRIM_SOURCE_CORE = "core"
TRIM_SOURCE_STRUCTURAL_MSA = "structural_msa"
TRIM_SOURCE_CONTINUOUS_MSA = "continuous_msa"
TRIM_SOURCE_MAFFT = "mafft"

TRIM_SOURCE_EQUIV_TSV = {
    TRIM_SOURCE_CORE: "ssm_core_equivalences.tsv",
    TRIM_SOURCE_STRUCTURAL_MSA: "ssm_structural_msa_equivalences.tsv",
    TRIM_SOURCE_CONTINUOUS_MSA: "ssm_structural_msa_equivalences.tsv",
    TRIM_SOURCE_MAFFT: "ssm_msa_equivalences.tsv",
}

TRIM_SOURCE_OUTPUT_PREFIX = {
    TRIM_SOURCE_CORE: "SSM_aligned_core_",
    TRIM_SOURCE_STRUCTURAL_MSA: "SSM_aligned_structural_msa_",
    TRIM_SOURCE_CONTINUOUS_MSA: "SSM_aligned_continuous_msa_",
    TRIM_SOURCE_MAFFT: "SSM_aligned_mafft_msa_",
}


def aligned_trim_output_dir(ssm_dir: str | os.PathLike, trim_source: str) -> str:
    """Sibling output directory for aligned trim (core or MSA modes)."""
    ssm_dir = os.path.abspath(os.path.expanduser(str(ssm_dir)))
    parent = os.path.dirname(ssm_dir)
    run_label = ssm_output_run_label(ssm_dir)
    prefix = TRIM_SOURCE_OUTPUT_PREFIX.get(trim_source, TRIM_SOURCE_OUTPUT_PREFIX[TRIM_SOURCE_CORE])
    return os.path.join(parent, "{}{}".format(prefix, run_label))


def trim_equivalences_tsv_path(ssm_dir: str | os.PathLike, trim_source: str) -> str:
    """Path to equivalence TSV inside an SSM output directory."""
    ssm_dir = os.path.abspath(os.path.expanduser(str(ssm_dir)))
    name = TRIM_SOURCE_EQUIV_TSV.get(trim_source, TRIM_SOURCE_EQUIV_TSV[TRIM_SOURCE_CORE])
    return os.path.join(ssm_dir, name)


def parse_log_models_directory(log_file: str) -> str | None:
    """Return model source directory from ``# Directories:`` in a Coot log header."""
    log_file = os.path.abspath(os.path.expanduser(log_file))
    try:
        with open(log_file, encoding="utf-8", errors="replace") as fh:
            for _ in range(60):
                line = fh.readline()
                if not line:
                    break
                m = _RE_LOG_DIRECTORIES.match(line.rstrip("\n\r"))
                if m:
                    raw = m.group(1).strip()
                    # AxB logs may list several directories comma-separated
                    if "," in raw:
                        raw = raw.split(",", 1)[0].strip()
                    return os.path.expanduser(raw)
    except OSError:
        return None
    return None


def read_name_map(path: str) -> dict[str, str]:
    """Parse ``{run_label}_name_map.tsv`` → short_stem → original_basename."""
    path = os.path.abspath(os.path.expanduser(path))
    mapping: dict[str, str] = {}
    try:
        with open(path, encoding="utf-8", errors="replace") as fh:
            for raw in fh:
                line = raw.strip()
                if not line or line.startswith("#"):
                    continue
                parts = line.split("\t")
                if len(parts) < 2:
                    continue
                short_stem, original = parts[0].strip(), parts[1].strip()
                if short_stem.lower() in ("short_stem", "short", "label"):
                    continue
                if short_stem and original:
                    mapping[short_stem] = original
    except OSError:
        return {}
    return mapping


def resolve_models_directory(
    ssm_dir: str,
    log_file: str | None = None,
) -> tuple[str | None, str]:
    """
    Resolve native models directory for an SSM output folder.

    Tries the Coot log ``# Directories:`` path first, then fallbacks when that
    path is stale (e.g. parent folder renamed). Returns (path, source_note).
    """
    ssm_dir = os.path.abspath(os.path.expanduser(str(ssm_dir)))
    parent = os.path.dirname(ssm_dir)
    run_label = ssm_output_run_label(ssm_dir)

    if log_file is None:
        log_file = discover_coot_log(ssm_dir)

    log_path: str | None = None
    parsed: str | None = None
    if log_file and os.path.isfile(log_file):
        log_path = os.path.abspath(os.path.expanduser(log_file))
        parsed = parse_log_models_directory(log_path)

    candidates: list[tuple[str, str]] = []
    if parsed:
        candidates.append((parsed, "coot_log"))
        reanchored = os.path.join(parent, os.path.basename(parsed.rstrip("/\\")))
        if reanchored != parsed:
            candidates.append((reanchored, "reanchored"))

    short_run = os.path.join(parent, "_short_{}".format(run_label))
    candidates.append((short_run, "_short_run_label"))

    try:
        short_dirs = sorted(
            os.path.join(parent, name)
            for name in os.listdir(parent)
            if name.startswith("_short_")
            and os.path.isdir(os.path.join(parent, name))
        )
    except OSError:
        short_dirs = []
    if len(short_dirs) == 1:
        candidates.append((short_dirs[0], "_short_glob"))
    elif len(short_dirs) > 1:
        for path in short_dirs:
            if run_label and path.endswith("_short_{}".format(run_label)):
                candidates.append((path, "_short_glob"))

    name_map_path = os.path.join(parent, "{}_name_map.tsv".format(run_label))
    if os.path.isfile(name_map_path) and read_name_map(name_map_path):
        candidates.append((parent, "name_map_parent"))

    seen: set[str] = set()
    for path, source in candidates:
        path = os.path.abspath(os.path.expanduser(path))
        if path in seen:
            continue
        seen.add(path)
        if os.path.isdir(path):
            return path, source

    return None, ""


def read_fasta(path: str) -> list[tuple[str, str]]:
    """Return (name, sequence) records from a FASTA file."""
    path = os.path.abspath(os.path.expanduser(path))
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


def _is_sparse_tsv_absent(chain: str, resnum: str) -> bool:
    chain = (chain or "").strip()
    resnum = (resnum or "").strip()
    return not chain or chain == "." or not resnum or resnum == "."


def write_sparse_equivalences_tsv(
    path: str,
    labels: list[str],
    columns: list[dict[str, ResidueKey]],
    aa_codes: dict[tuple[str, tuple], str],
    *,
    header_comment: str = "",
) -> None:
    """Write equivalence TSV; absent structures use ``.`` in all four fields."""
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
    with open(path, "w", encoding="utf-8") as fh:
        if header_comment:
            fh.write("# {}\n".format(header_comment))
        fh.write("\t".join(cols) + "\n")
        for i, mapping in enumerate(columns, start=1):
            row = [str(i)]
            for label in labels:
                if label in mapping:
                    key = mapping[label]
                    aa = aa_codes.get((label, key.as_tuple()), "X")
                    row.extend([key.chain, str(key.resnum), key.icode, aa])
                else:
                    row.extend([".", ".", ".", "."])
            fh.write("\t".join(row) + "\n")


def read_sparse_equivalences_tsv(
    path: str,
) -> tuple[list[str], list[dict[str, ResidueKey]]]:
    """
    Parse sparse equivalence TSV (structural MSA, MAFFT, or dense core).

    Returns (structure_labels, columns) where each column maps only present labels.
    """
    path = os.path.abspath(os.path.expanduser(path))
    labels: list[str] = []
    columns: list[dict[str, ResidueKey]] = []
    with open(path, encoding="utf-8", errors="replace") as fh:
        header: list[str] | None = None
        for raw in fh:
            line = raw.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if header is None:
                header = parts
                if parts[0] != "column":
                    raise ValueError(
                        "Expected 'column' as first header field in {}".format(path)
                    )
                rest = parts[1:]
                for i in range(0, len(rest), 4):
                    if i + 3 < len(rest):
                        field = rest[i]
                        if field.endswith("_chain"):
                            labels.append(field[: -len("_chain")])
                continue
            if len(parts) < 1 + 4 * len(labels):
                continue
            mapping: dict[str, ResidueKey] = {}
            base = 1
            for j, label in enumerate(labels):
                off = base + j * 4
                chain = parts[off]
                resnum_s = parts[off + 1]
                icode = parts[off + 2]
                if _is_sparse_tsv_absent(chain, resnum_s):
                    continue
                mapping[label] = ResidueKey(chain, int(resnum_s), icode)
            if mapping:
                columns.append(mapping)
    return labels, columns


def read_core_equivalences_tsv(
    path: str,
) -> tuple[list[str], list[dict[str, ResidueKey]]]:
    """
    Parse ``ssm_core_equivalences.tsv`` written by SSM_struct_core.py.

    Returns (structure_labels, core_columns).
    """
    return read_sparse_equivalences_tsv(path)


def _alignment_identity(a: str, b: str) -> tuple[int, int, float]:
    """Return (n_columns, n_identical, pct_identity) for gapped sequences."""
    if len(a) != len(b):
        n = min(len(a), len(b))
        a, b = a[:n], b[:n]
    n_cols = 0
    n_id = 0
    for ca, cb in zip(a, b):
        if ca == "-" or cb == "-":
            continue
        n_cols += 1
        if ca.upper() == cb.upper():
            n_id += 1
    pct = (100.0 * n_id / n_cols) if n_cols else 0.0
    return n_cols, n_id, pct


def pairwise_sequence_identity(
    sequences: dict[str, str],
) -> list[dict[str, object]]:
    """Pairwise % identity over columns where both sequences are non-gap."""
    labels = sorted(sequences.keys(), key=_natural_sort_key)
    rows: list[dict[str, object]] = []
    for i, la in enumerate(labels):
        for lb in labels[i:]:
            n_cols, n_id, pct = _alignment_identity(sequences[la], sequences[lb])
            rows.append(
                {
                    "structure_a": la,
                    "structure_b": lb,
                    "n_core_columns": n_cols,
                    "n_identical": n_id,
                    "pct_identity": round(pct, 2),
                }
            )
    return rows


def pairwise_core_alignment_metrics(
    sequences: dict[str, str],
    core_rmsd: dict[tuple[str, str], float],
) -> list[dict[str, object]]:
    """Merge pairwise % identity with core Cα RMSD (Å) for each structure pair."""
    rows = pairwise_sequence_identity(sequences)
    for row in rows:
        key = (str(row["structure_a"]), str(row["structure_b"]))
        rmsd = core_rmsd.get(key)
        row["core_rmsd_angstrom"] = rmsd if rmsd is not None else ""
    return rows


def core_rmsd_summary(
    core_rmsd: dict[tuple[str, str], float],
) -> list[dict[str, object]]:
    """Mean core Cα RMSD to other structures (excludes self)."""
    labels = sorted(
        {a for pair in core_rmsd for a in pair},
        key=_natural_sort_key,
    )
    rows: list[dict[str, object]] = []
    for la in labels:
        others = [lb for lb in labels if lb != la]
        if not others:
            rows.append(
                {
                    "structure": la,
                    "n_comparisons": 0,
                    "mean_core_rmsd_angstrom": "",
                }
            )
            continue
        rmsds: list[float] = []
        for lb in others:
            key = (la, lb) if _natural_sort_key(la) <= _natural_sort_key(lb) else (lb, la)
            val = core_rmsd.get(key)
            if val is not None and val == val:
                rmsds.append(float(val))
        mean_rmsd = sum(rmsds) / len(rmsds) if rmsds else float("nan")
        rows.append(
            {
                "structure": la,
                "n_comparisons": len(rmsds),
                "mean_core_rmsd_angstrom": round(mean_rmsd, 4)
                if mean_rmsd == mean_rmsd
                else "",
            }
        )
    return rows


def sequence_identity_summary(
    sequences: dict[str, str],
) -> list[dict[str, object]]:
    """Mean pairwise % identity to other structures (excludes self)."""
    labels = sorted(sequences.keys(), key=_natural_sort_key)
    rows: list[dict[str, object]] = []
    for la in labels:
        others = [lb for lb in labels if lb != la]
        if not others:
            rows.append(
                {
                    "structure": la,
                    "n_comparisons": 0,
                    "mean_pct_identity": "",
                }
            )
            continue
        pcts: list[float] = []
        for lb in others:
            _n, _id, pct = _alignment_identity(sequences[la], sequences[lb])
            pcts.append(pct)
        mean_pct = sum(pcts) / len(pcts)
        rows.append(
            {
                "structure": la,
                "n_comparisons": len(others),
                "mean_pct_identity": round(mean_pct, 2),
            }
        )
    return rows


def _parse_core_rmsd_angstrom(block: SSMSuperpositionBlock) -> float | None:
    """Parse core RMSD (Å) from an SSM superposition block."""
    line = block.stats.get("core rmsd line", "")
    if not line:
        for key, val in block.stats.items():
            if "core rmsd" in key.lower():
                line = val
                break
    m = re.search(
        r"(?:INFO:\s*)?core\s+rmsd(?:\s+achieved)?\s*:\s*([\d.]+)", line, re.I
    )
    if m:
        return float(m.group(1))
    m = re.search(r"(?:rmsd|rms\s+devi)\s*[=:]?\s*([\d.]+)", line, re.I)
    return float(m.group(1)) if m else None


def _parse_pct_identity_stat(value: str) -> float | None:
    """Parse numeric percentage from a Coot ``sequence identity`` stat line."""
    if not value:
        return None
    m = re.search(r"([\d.]+)\s*%?", value.strip())
    return float(m.group(1)) if m else None


def _normalize_coot_metric_label(label: str) -> str:
    """Strip phase suffixes so anchor ``X`` matches reference ``X_trim`` / ``X_core``."""
    for suffix in ("_trim", "_core", "_cont"):
        if label.endswith(suffix):
            return label[: -len(suffix)]
    return label


def _metric_reference_is_anchor(reference: str, anchor: str | None) -> bool:
    if not anchor:
        return True
    return _normalize_coot_metric_label(reference) == _normalize_coot_metric_label(
        anchor
    )


def pairwise_metrics_from_ssm_log(log_file: str) -> list[dict[str, object]]:
    """
    Pairwise core RMSD and % sequence identity from a Coot SSM log.

    One row per superposition block (typically moving → reference).
    """
    blocks = parse_ssm_log(log_file)
    rows: list[dict[str, object]] = []
    for block in blocks:
        rmsd = _parse_core_rmsd_angstrom(block)
        pct = _parse_pct_identity_stat(block.stats.get("sequence identity", ""))
        n_aligned = block.stats.get("number of aligned residues", "")
        if not n_aligned:
            n_aligned = block.stats.get(
                "number of residues in aligned sections (moving)", ""
            )
        rows.append(
            {
                "structure_a": block.moving_label,
                "structure_b": block.reference_label,
                "n_aligned_residues": n_aligned,
                "pct_identity": round(pct, 2) if pct is not None else "",
                "core_rmsd_angstrom": round(rmsd, 4) if rmsd is not None else "",
            }
        )
    return rows


def coot_metrics_rmsd_summary(
    metric_rows: list[dict[str, object]],
    *,
    anchor: str | None = None,
) -> list[dict[str, object]]:
    """Mean core RMSD to anchor (or all partners) per moving structure."""
    by_structure: dict[str, list[float]] = {}
    for row in metric_rows:
        mov = str(row.get("structure_a", ""))
        ref = str(row.get("structure_b", ""))
        if not _metric_reference_is_anchor(ref, anchor):
            continue
        val = row.get("core_rmsd_angstrom")
        if val == "" or val is None:
            continue
        try:
            by_structure.setdefault(mov, []).append(float(val))
        except (TypeError, ValueError):
            continue
    labels = sorted(by_structure.keys(), key=_natural_sort_key)
    out: list[dict[str, object]] = []
    for lab in labels:
        rmsds = by_structure[lab]
        mean_r = sum(rmsds) / len(rmsds) if rmsds else float("nan")
        out.append(
            {
                "structure": lab,
                "n_comparisons": len(rmsds),
                "mean_core_rmsd_angstrom": round(mean_r, 4)
                if mean_r == mean_r
                else "",
            }
        )
    return out


def coot_metrics_identity_summary(
    metric_rows: list[dict[str, object]],
    *,
    anchor: str | None = None,
) -> list[dict[str, object]]:
    """Mean % sequence identity to anchor per moving structure."""
    by_structure: dict[str, list[float]] = {}
    for row in metric_rows:
        mov = str(row.get("structure_a", ""))
        ref = str(row.get("structure_b", ""))
        if not _metric_reference_is_anchor(ref, anchor):
            continue
        val = row.get("pct_identity")
        if val == "" or val is None:
            continue
        try:
            by_structure.setdefault(mov, []).append(float(val))
        except (TypeError, ValueError):
            continue
    labels = sorted(by_structure.keys(), key=_natural_sort_key)
    out: list[dict[str, object]] = []
    for lab in labels:
        pcts = by_structure[lab]
        mean_p = sum(pcts) / len(pcts) if pcts else float("nan")
        out.append(
            {
                "structure": lab,
                "n_comparisons": len(pcts),
                "mean_pct_identity": round(mean_p, 2) if mean_p == mean_p else "",
            }
        )
    return out


def find_ssm_superposition_block(
    log_file: str,
    moving_label: str,
    reference_label: str,
):
    """Return the SSM block for ``moving_label`` onto ``reference_label``, if present."""
    for block in parse_ssm_log(log_file):
        if (
            block.moving_label == moving_label
            and block.reference_label == reference_label
        ):
            return block
    return None


_RE_COOT_CANT_MAKE_GRAPH = re.compile(r"can't make graph for\s+(\S+)", re.I)
_RE_COOT_SSM_SUPERPOSE_FAIL = re.compile(
    r"SSM superposition failed for\s+(\S+)", re.I
)


def audit_coot_ssm_log(log_file: str) -> list[str]:
    """
    Return human-readable failure lines from a Coot SSM log.

    Coot often prints ``can't make graph`` without raising a Python exception,
    so callers should audit logs even when the process exit code is 0.
    """
    log_file = os.path.abspath(os.path.expanduser(log_file))
    if not os.path.isfile(log_file):
        return []
    failures: list[str] = []
    seen: set[str] = set()
    with open(log_file, encoding="utf-8", errors="replace") as fh:
        for raw in fh:
            line = raw.rstrip("\n\r")
            for pattern in (_RE_COOT_CANT_MAKE_GRAPH, _RE_COOT_SSM_SUPERPOSE_FAIL):
                m = pattern.search(line)
                if not m:
                    continue
                msg = line.strip()
                if msg in seen:
                    continue
                seen.add(msg)
                failures.append(msg)
    return failures


def metrics_rows_have_rmsd(metric_rows: list[dict[str, object]]) -> bool:
    """True when at least one row has a numeric core RMSD."""
    for row in metric_rows:
        val = row.get("core_rmsd_angstrom")
        if val == "" or val is None:
            continue
        try:
            float(val)
            return True
        except (TypeError, ValueError):
            continue
    return False


def validate_step1_ssm_outputs(
    ssm_dir: str,
    log_file: str | None = None,
    *,
    require_equivalences: bool = True,
    min_structures: int = 2,
) -> tuple[bool, list[str]]:
    """
    Validate Step 1 SSM all-vs-all outputs (and Step 2 equivalences if present).

    Returns ``(ok, issues)`` where *issues* lists human-readable problems.
    """
    issues: list[str] = []
    ssm_dir = os.path.abspath(os.path.expanduser(ssm_dir))
    log_file = resolve_step1_coot_log(ssm_dir, log_file)
    if not log_file or not os.path.isfile(log_file):
        issues.append("missing or unreadable Coot SSM log")
        return False, issues

    blocks = parse_ssm_log(log_file)
    substantive = substantive_ssm_pair_count(log_file)
    if substantive <= 0:
        issues.append("no SSM superposition blocks with RMSD or equivalences in log")
        return False, issues

    labels: set[str] = set()
    for block in blocks:
        labels.add(block.moving_label)
        labels.add(block.reference_label)
    n = len(labels)
    expected = n * (n - 1) if n >= min_structures else 0
    if n < min_structures:
        issues.append(
            "fewer than {} structures in log (found {})".format(
                min_structures, n
            )
        )
    elif expected and substantive < expected:
        issues.append(
            "incomplete substantive SSM: {} pairs with RMSD/equiv, expected {} for {} structures".format(
                substantive, expected, n
            )
        )

    valid_blocks = [b for b in blocks if block_has_valid_ssm(b)]
    graph_failures = _graph_errors_without_blocks(log_file, valid_blocks)
    if graph_failures:
        issues.append(
            "Coot SSM graph errors with no successful block ({} message(s))".format(
                len(graph_failures)
            )
        )
        for msg in graph_failures[:5]:
            issues.append("  " + msg)
        if len(graph_failures) > 5:
            issues.append("  ... and {} more".format(len(graph_failures) - 5))

    if require_equivalences:
        equiv_dir = os.path.join(ssm_dir, "ssm_equivalences")
        if not os.path.isdir(equiv_dir):
            issues.append(
                "missing ssm_equivalences/ (run extract_rmsd with --seq-align)"
            )
        else:
            tsv_files = [
                f
                for f in os.listdir(equiv_dir)
                if f.endswith(".tsv")
            ]
            if not tsv_files:
                issues.append(
                    "ssm_equivalences/ has no TSV files "
                    "(run extract_rmsd with --seq-align)"
                )

    return (len(issues) == 0, issues)


def format_clustal_like(
    labels: list[str],
    sequences: list[str],
    *,
    title: str = "FoldKit SSM",
) -> str:
    width = max((len(s) for s in sequences), default=0)
    if width == 0:
        return ""
    name_width = max(len(l) for l in labels)
    lines = [title + "\n"]
    chunk = 60
    for start in range(0, width, chunk):
        end = min(start + chunk, width)
        for label, seq in zip(labels, sequences):
            lines.append("{:<{nw}}  {}\n".format(label, seq[start:end], nw=name_width))
        lines.append("\n")
    return "".join(lines)
