#!/usr/bin/env python3
"""
Rewrite chain IDs in PDB files (ATOM/HETATM/ANISOU/TER). Optionally merge one chain into
another and renumber residues so the source chain continues after the last residue
of the target chain (e.g. B → A with numbering after the last A residue).

Examples:
  # Rename chain B to A only (fails if chain A already has atoms):
  python pdb_rechain.py models/ --pattern '*_models.pdb' --from B --to A

  # Merge B into A: all former-B atoms become chain A, residue numbers continue
  # after the highest residue number on the current chain A:
  python pdb_rechain.py models/ --pattern '*_models.pdb' --from B --to A \\
      --merge-renumber -o rechained/

  # Single file:
  python pdb_rechain.py model_01.pdb --from B --to A --merge-renumber -o model_01_merged.pdb

  # Multiple merges in one pass (applies left-to-right):
  python pdb_rechain.py multicopy.pdb --merge-map 'B:A,D:C,F:E' -o merged.pdb

  # Multiple merges + then relabel remaining chains sequentially A..Z + renumber residues per chain:
  python pdb_rechain.py multicopy.pdb --merge-map 'B:A,D:C,F:E,H:G' --rename-sequential --renumber-per-chain -o cleaned.pdb

  # Multi-copy assemblies: control the A..Z assignment so merged pairs map predictably.
  # Example merge list:
  #   B→A, D→C, F→E, H→G, J→I, L→K, N→M, O→P
  # After these merges, the surviving chain IDs are:
  #   A, C, E, G, I, K, M, P
  # Use --chain-order so merged AB becomes new A, merged CD becomes new B, etc.
  python pdb_rechain.py multicopy.pdb \
      --merge-map 'B:A,D:C,F:E,H:G,J:I,L:K,N:M,O:P' \
      --reorder-chains --rename-sequential --chain-order 'A,C,E,G,I,K,M,P' \
      --renumber-per-chain -o multicopy_merged_ordered.pdb
"""

from __future__ import annotations

import argparse
import glob
import os
import sys
from typing import Iterable


def _is_atom_line(line: str) -> bool:
    return line.startswith(("ATOM  ", "HETATM"))


def _is_anisou_line(line: str) -> bool:
    return line.startswith("ANISOU")


def _is_coord_or_anisou_line(line: str) -> bool:
    """
    Records that carry chain+resseq in columns 22–26 and should be kept consistent.
    """
    return _is_atom_line(line) or _is_anisou_line(line)


def _is_ter_line(line: str) -> bool:
    return line.startswith("TER   ") or line.startswith("TER\n")


def _parse_atom_chain_resseq_icode(line: str) -> tuple[str, int, str] | None:
    if len(line) < 27:
        return None
    chain = line[21]
    try:
        resseq = int(line[22:26])
    except ValueError:
        return None
    icode = line[26] if len(line) > 26 else " "
    return (chain, resseq, icode)


def _set_chain_resseq(line: str, chain: str, resseq: int) -> str:
    """Rewrite columns 22–26 (chain + 4-digit resseq); keep rest of line."""
    if len(line) < 27:
        return line
    ch = chain[:1]
    rs = max(-999, min(9999, int(resseq)))
    body = f"{ch}{rs:4d}" + line[26:]
    return line[:21] + body


def _max_resseq_on_chain(lines: Iterable[str], chain: str) -> int | None:
    mx = None
    for line in lines:
        if not _is_atom_line(line):
            continue
        p = _parse_atom_chain_resseq_icode(line)
        if not p or p[0] != chain:
            continue
        if p[1] > (mx if mx is not None else p[1] - 1):
            mx = p[1]
    return mx


def _min_resseq_on_chain(lines: Iterable[str], chain: str) -> int | None:
    mn = None
    for line in lines:
        if not _is_atom_line(line):
            continue
        p = _parse_atom_chain_resseq_icode(line)
        if not p or p[0] != chain:
            continue
        r = p[1]
        if mn is None or r < mn:
            mn = r
    return mn


def _chain_has_atoms(lines: Iterable[str], chain: str) -> bool:
    for line in lines:
        if not _is_atom_line(line):
            continue
        p = _parse_atom_chain_resseq_icode(line)
        if p and p[0] == chain:
            return True
    return False


def _parse_merge_map(text: str) -> list[tuple[str, str]]:
    """
    Parse merge map like 'B:A,D:C,F:E' into [('B','A'), ...].
    """
    pairs: list[tuple[str, str]] = []
    for raw in (text or "").split(","):
        raw = raw.strip()
        if not raw:
            continue
        if ":" not in raw:
            raise ValueError(f"invalid merge-map item {raw!r} (expected FROM:TO)")
        a, b = raw.split(":", 1)
        a, b = a.strip(), b.strip()
        if len(a) != 1 or len(b) != 1:
            raise ValueError(f"invalid merge-map item {raw!r} (chain IDs must be 1 char)")
        pairs.append((a, b))
    if not pairs:
        raise ValueError("merge-map is empty")
    return pairs


def _apply_merge_map(
    lines: list[str],
    merge_map: list[tuple[str, str]],
) -> tuple[list[str] | None, str | None]:
    """
    Apply multiple --merge-renumber operations left-to-right.
    """
    cur = lines
    for cf, ct in merge_map:
        cur2, err = _merge_renumber(cur, cf, ct)
        if err:
            return None, f"merge {cf!r}->{ct!r} failed: {err}"
        assert cur2 is not None
        cur = cur2
    return cur, None


def _chains_in_order(lines: Iterable[str]) -> list[str]:
    """
    Return chain IDs in first-appearance order (ATOM/HETATM only).
    """
    seen: set[str] = set()
    order: list[str] = []
    for line in lines:
        if not _is_atom_line(line):
            continue
        p = _parse_atom_chain_resseq_icode(line)
        if not p:
            continue
        ch = p[0]
        if ch not in seen:
            seen.add(ch)
            order.append(ch)
    return order


def _parse_chain_order(text: str) -> list[str]:
    """
    Parse a comma-separated chain order like 'A,C,E,G' into ['A','C','E','G'].
    """
    out: list[str] = []
    for raw in (text or "").split(","):
        raw = raw.strip()
        if not raw:
            continue
        if len(raw) != 1:
            raise ValueError(f"invalid chain ID {raw!r} in --chain-order (must be 1 char)")
        if raw not in out:
            out.append(raw)
    if not out:
        raise ValueError("--chain-order is empty")
    return out


def _build_sequential_chain_mapping(
    present_order: list[str],
    preferred_order: list[str] | None,
) -> dict[str, str]:
    """
    Build mapping old_chain -> new_chain (A..Z) using:
    - preferred_order first (only those that are present),
    - then any remaining chains in present_order.
    """
    order: list[str] = []
    present_set = set(present_order)
    if preferred_order:
        for c in preferred_order:
            if c in present_set and c not in order:
                order.append(c)
    for c in present_order:
        if c not in order:
            order.append(c)
    if len(order) > 26:
        raise ValueError(f"rename supports at most 26 chains; found {len(order)}")
    return {old: chr(ord("A") + i) for i, old in enumerate(order)}


def _rename_chains_by_map(lines: list[str], mapping: dict[str, str]) -> list[str]:
    """
    Rename chain IDs for ATOM/HETATM/TER records by mapping.
    """
    out: list[str] = []
    for line in lines:
        if _is_coord_or_anisou_line(line) or _is_ter_line(line):
            p = _parse_atom_chain_resseq_icode(line)
            if p:
                ch, rseq, _icode = p
                new_ch = mapping.get(ch)
                if new_ch:
                    line = _set_chain_resseq(line, new_ch, rseq)
        out.append(line)
    return out


def _renumber_residues_per_chain(lines: list[str]) -> tuple[list[str] | None, str | None]:
    """
    Renumber residues for each chain separately starting at 1, in residue first-appearance order.
    Keeps insertion codes (column 27) unchanged.
    """
    new_index: dict[tuple[str, int, str], int] = {}
    next_for_chain: dict[str, int] = {}

    out: list[str] = []
    for line in lines:
        if not _is_coord_or_anisou_line(line) and not _is_ter_line(line):
            out.append(line)
            continue
        p = _parse_atom_chain_resseq_icode(line)
        if not p:
            out.append(line)
            continue
        ch, rseq, icode = p
        key = (ch, rseq, icode)
        if key not in new_index:
            nxt = next_for_chain.get(ch, 1)
            if nxt > 9999:
                return None, f"renumber-per-chain would exceed 9999 on chain {ch!r}"
            new_index[key] = nxt
            next_for_chain[ch] = nxt + 1
        line = _set_chain_resseq(line, ch, new_index[key])
        out.append(line)
    return out, None


def _reorder_coordinate_blocks_by_chain(lines: list[str]) -> list[str]:
    """
    Reorder coordinate-like records so each chain's coordinates are contiguous.

    Keeps non-coordinate records in their original relative order, but moves all
    ATOM/HETATM/ANISOU/TER blocks to the end (before END/ENDMDL if present) and
    groups them by chain in first-appearance order.
    """
    end_records = {"END", "ENDMDL"}

    header: list[str] = []
    tail_end: list[str] = []
    coord: list[str] = []

    for line in lines:
        tag = line[:6].strip()
        if tag in end_records:
            tail_end.append(line)
        elif _is_coord_or_anisou_line(line) or _is_ter_line(line):
            coord.append(line)
        else:
            header.append(line)

    chain_order = _chains_in_order(coord)
    # Include chains that might only appear on TER/ANISOU after edits
    extra: list[str] = []
    seen = set(chain_order)
    for line in coord:
        p = _parse_atom_chain_resseq_icode(line)
        if p and p[0] not in seen:
            seen.add(p[0])
            extra.append(p[0])
    chain_order.extend(extra)

    by_chain: dict[str, list[str]] = {c: [] for c in chain_order}
    other: list[str] = []
    for line in coord:
        p = _parse_atom_chain_resseq_icode(line)
        if not p:
            other.append(line)
            continue
        ch = p[0]
        if ch in by_chain:
            by_chain[ch].append(line)
        else:
            other.append(line)

    out: list[str] = []
    out.extend(header)
    out.extend(other)
    for ch in chain_order:
        out.extend(by_chain.get(ch, []))
    out.extend(tail_end)
    return out


def _rename_chain_only(
    lines: list[str],
    chain_from: str,
    chain_to: str,
) -> tuple[list[str] | None, str | None]:
    """Replace chain ID; target chain must have no atoms."""
    cf, ct = chain_from[:1], chain_to[:1]
    if _chain_has_atoms(lines, ct):
        return None, (
            f"chain {ct!r} already has atoms; use --merge-renumber to merge "
            f"{cf!r} into {ct!r} with new residue numbers, or choose an empty target chain."
        )
    out = []
    for line in lines:
        if _is_coord_or_anisou_line(line) or _is_ter_line(line):
            p = _parse_atom_chain_resseq_icode(line)
            if p and p[0] == cf:
                line = _set_chain_resseq(line, ct, p[1])
        out.append(line)
    return out, None


def _merge_renumber(
    lines: list[str],
    chain_from: str,
    chain_to: str,
) -> tuple[list[str] | None, str | None]:
    """
    Move all atoms from chain_from to chain_to; residue numbers on former chain_from
    continue after the maximum residue number on chain_to (integer resseq only).
    """
    cf, ct = chain_from[:1], chain_to[:1]
    if cf == ct:
        return None, "--from and --to must differ for --merge-renumber"
    if not _chain_has_atoms(lines, cf):
        return None, f"no ATOM/HETATM lines found for chain {cf!r}"
    if not _chain_has_atoms(lines, ct):
        return None, f"chain {ct!r} has no atoms; use rename mode without --merge-renumber"

    max_a = _max_resseq_on_chain(lines, ct)
    min_b = _min_resseq_on_chain(lines, cf)
    if max_a is None or min_b is None:
        return None, "could not determine residue ranges for merge"

    delta = max_a + 1 - min_b
    max_b = _max_resseq_on_chain(lines, cf)
    assert max_b is not None
    max_new = max_b + delta
    if max_new > 9999:
        return None, (
            f"renumbered residues would exceed 9999 (max would be {max_new}); "
            "cannot fit PDB fixed-width resseq field"
        )

    out: list[str] = []
    for line in lines:
        if _is_ter_line(line):
            # Drop TER that only terminated the moving chain; keep others loosely
            p = _parse_atom_chain_resseq_icode(line)
            if p and p[0] == cf:
                continue
            out.append(line)
            continue
        if not _is_coord_or_anisou_line(line):
            out.append(line)
            continue
        p = _parse_atom_chain_resseq_icode(line)
        if not p:
            out.append(line)
            continue
        ch, rseq, icode = p
        if ch != cf:
            out.append(line)
            continue
        new_r = rseq + delta
        line = _set_chain_resseq(line, ct, new_r)
        out.append(line)
    return out, None


def _collect_inputs(paths: list[str], pattern: str | None) -> list[str]:
    files: list[str] = []
    for p in paths:
        ap = os.path.abspath(p)
        if os.path.isfile(ap):
            files.append(ap)
            continue
        if os.path.isdir(ap):
            pat = pattern or "*.pdb"
            g = sorted(glob.glob(os.path.join(ap, pat)))
            g = [x for x in g if os.path.isfile(x)]
            if not g:
                print(f"Warning: no files matching {os.path.join(ap, pat)!r}", file=sys.stderr)
            files.extend(g)
            continue
        raise FileNotFoundError(f"not a file or directory: {p!r}")
    # de-dup preserving order
    seen = set()
    uniq = []
    for f in files:
        if f not in seen:
            seen.add(f)
            uniq.append(f)
    return uniq


def _default_out_path(in_path: str, out_dir: str | None, suffix: str) -> str:
    base = os.path.basename(in_path)
    if base.lower().endswith(".pdb"):
        stem = base[:-4]
        ext = ".pdb"
    else:
        stem, ext = base, ""
    name = f"{stem}{suffix}{ext}"
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)
        return os.path.join(out_dir, name)
    d = os.path.dirname(in_path)
    return os.path.join(d, name)


def process_file(
    in_path: str,
    chain_from: str,
    chain_to: str,
    merge_renumber: bool,
    merge_map: list[tuple[str, str]] | None,
    rename_sequential: bool,
    chain_order: list[str] | None,
    renumber_per_chain: bool,
    reorder_chains: bool,
) -> tuple[bool, str]:
    with open(in_path, encoding="utf-8", errors="replace") as f:
        lines = f.readlines()

    new_lines: list[str] | None = lines
    err: str | None = None

    if merge_map:
        new_lines, err = _apply_merge_map(new_lines, merge_map)
    else:
        if merge_renumber:
            new_lines, err = _merge_renumber(new_lines, chain_from, chain_to)
        else:
            new_lines, err = _rename_chain_only(new_lines, chain_from, chain_to)

    if err:
        return False, err
    assert new_lines is not None

    if rename_sequential:
        order = _chains_in_order(new_lines)
        try:
            mapping = _build_sequential_chain_mapping(order, chain_order)
        except ValueError as e:
            return False, str(e)
        new_lines = _rename_chains_by_map(new_lines, mapping)

    if renumber_per_chain:
        new_lines, err = _renumber_residues_per_chain(new_lines)
        if err:
            return False, err
        assert new_lines is not None

    # Optional: make each chain contiguous (prevents merged source chains from being
    # interleaved with target chains due to original record ordering).
    if reorder_chains:
        new_lines = _reorder_coordinate_blocks_by_chain(new_lines)

    return True, "".join(new_lines)


def main() -> None:
    ap = argparse.ArgumentParser(
        description="Replace PDB chain IDs; optional merge with continued residue numbering.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    ap.add_argument(
        "paths",
        nargs="+",
        help="PDB file(s) and/or directory(ies). For directories, use --pattern.",
    )
    ap.add_argument(
        "--from",
        "-f",
        dest="chain_from",
        metavar="ID",
        help="Current chain ID (one character).",
    )
    ap.add_argument(
        "--to",
        "-t",
        dest="chain_to",
        metavar="ID",
        help="New chain ID (one character).",
    )
    ap.add_argument(
        "--merge-map",
        metavar="MAP",
        help="Comma-separated list of FROM:TO merges applied left-to-right. Example: 'B:A,D:C,F:E'",
    )
    ap.add_argument(
        "--pattern",
        "-p",
        metavar="GLOB",
        help="Filename glob for each directory argument (e.g. '*_models.pdb').",
    )
    ap.add_argument(
        "--merge-renumber",
        action="store_true",
        help=(
            "Merge --from into --to: assign --to to all former --from atoms and set "
            "residue numbers to continue after the last integer resseq on --to. "
            "Requires both chains to exist in the file."
        ),
    )
    ap.add_argument(
        "--rename-sequential",
        action="store_true",
        help="After merging/renaming, relabel all chains in first-appearance order to A..Z.",
    )
    ap.add_argument(
        "--chain-order",
        metavar="IDS",
        help=(
            "When used with --rename-sequential: comma-separated preferred order of existing chain IDs "
            "(post-merge) to assign to A,B,C,... first. Example: --chain-order 'A,C,E,G,I,K,M,P' makes "
            "A->A, C->B, E->C, ... Remaining chains are appended in first-appearance order."
        ),
    )
    ap.add_argument(
        "--renumber-per-chain",
        action="store_true",
        help="After merging/renaming, renumber residues per chain starting at 1.",
    )
    ap.add_argument(
        "--reorder-chains",
        action="store_true",
        help="After merging/renaming, reorder coordinates so each chain is contiguous (helps avoid interleaving).",
    )
    ap.add_argument(
        "-o",
        "--output",
        metavar="PATH",
        help=(
            "Output file (single input only); output directory (multiple inputs); or a path template "
            "containing '{}' (multiple inputs) where '{}' is replaced by the input stem. "
            "Default: same directory as input with --suffix before .pdb."
        ),
    )
    ap.add_argument(
        "--suffix",
        default="_rechain",
        help="Suffix for output basename when -o is not set (default: %(default)s).",
    )
    ap.add_argument(
        "--dry-run",
        action="store_true",
        help="Print planned outputs only.",
    )
    args = ap.parse_args()

    merge_map = None
    if args.merge_map:
        try:
            merge_map = _parse_merge_map(args.merge_map)
        except ValueError as e:
            print(f"Error: {e}", file=sys.stderr)
            sys.exit(2)

    if merge_map is None:
        if not args.chain_from or not args.chain_to:
            print("Error: either provide --merge-map or provide both --from and --to.", file=sys.stderr)
            sys.exit(2)
        if len(args.chain_from) != 1 or len(args.chain_to) != 1:
            print("Error: --from and --to must be single-character chain IDs.", file=sys.stderr)
            sys.exit(2)

    chain_order = None
    if args.chain_order:
        try:
            chain_order = _parse_chain_order(args.chain_order)
        except ValueError as e:
            print(f"Error: {e}", file=sys.stderr)
            sys.exit(2)

    try:
        inputs = _collect_inputs(args.paths, args.pattern)
    except FileNotFoundError as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(2)

    if not inputs:
        print("Error: no input PDB files.", file=sys.stderr)
        sys.exit(2)

    for d in args.paths:
        if os.path.isdir(d) and not args.pattern:
            print(
                "Error: --pattern is required when a source is a directory "
                "(e.g. --pattern '*_batch*.pdb').",
                file=sys.stderr,
            )
            sys.exit(2)

    out_dir_multi: str | None = None
    out_template_multi: str | None = None
    if len(inputs) > 1:
        if not args.output:
            print(
                "Error: for multiple input files, -o/--output is required (directory, or a path template containing '{}').",
                file=sys.stderr,
            )
            sys.exit(2)
        out_arg = os.path.abspath(args.output)
        if "{}" in out_arg:
            out_template_multi = out_arg
            parent = os.path.dirname(out_template_multi) or "."
            os.makedirs(parent, exist_ok=True)
        else:
            out_dir_multi = out_arg
            if os.path.isfile(out_dir_multi):
                print(
                    "Error: for multiple inputs, -o must be a directory or a path template containing '{}'.",
                    file=sys.stderr,
                )
                sys.exit(2)
            os.makedirs(out_dir_multi, exist_ok=True)

    ok_n = 0
    for in_path in inputs:
        if len(inputs) == 1:
            if args.output:
                oa = os.path.abspath(args.output)
                if os.path.isdir(oa):
                    dest = os.path.join(
                        oa,
                        os.path.basename(_default_out_path(in_path, None, args.suffix)),
                    )
                else:
                    dest = oa
                    if not dest.lower().endswith(".pdb"):
                        dest = dest + ".pdb"
            else:
                dest = _default_out_path(in_path, None, args.suffix)
        else:
            if out_template_multi:
                stem = os.path.splitext(os.path.basename(in_path))[0]
                dest = out_template_multi.replace("{}", stem)
                if not dest.lower().endswith(".pdb"):
                    dest = dest + ".pdb"
            else:
                dest = _default_out_path(in_path, out_dir_multi, args.suffix)

        if args.dry_run:
            print(f"{in_path} -> {dest}")
            ok_n += 1
            continue

        ok, data = process_file(
            in_path,
            args.chain_from or "",
            args.chain_to or "",
            bool(args.merge_renumber),
            merge_map,
            bool(args.rename_sequential),
            chain_order,
            bool(args.renumber_per_chain),
            bool(args.reorder_chains),
        )
        if not ok:
            print(f"Error [{in_path}]: {data}", file=sys.stderr)
            continue
        with open(dest, "w", encoding="utf-8", newline="\n") as fo:
            fo.write(data)
        print(dest)
        ok_n += 1

    if ok_n == 0:
        sys.exit(1)


if __name__ == "__main__":
    main()
