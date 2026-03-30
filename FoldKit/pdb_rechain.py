#!/usr/bin/env python3
"""
Replace chain IDs in PDB files (ATOM/HETATM/TER). Optionally merge one chain into
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
"""

from __future__ import annotations

import argparse
import glob
import os
import sys
from typing import Iterable


def _is_atom_line(line: str) -> bool:
    return line.startswith(("ATOM  ", "HETATM"))


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
        if _is_atom_line(line) or _is_ter_line(line):
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
        if not _is_atom_line(line):
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
) -> tuple[bool, str]:
    with open(in_path, encoding="utf-8", errors="replace") as f:
        lines = f.readlines()

    if merge_renumber:
        new_lines, err = _merge_renumber(lines, chain_from, chain_to)
    else:
        new_lines, err = _rename_chain_only(lines, chain_from, chain_to)

    if err:
        return False, err
    assert new_lines is not None
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
        required=True,
        help="Current chain ID (one character).",
    )
    ap.add_argument(
        "--to",
        "-t",
        dest="chain_to",
        metavar="ID",
        required=True,
        help="New chain ID (one character).",
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
        "-o",
        "--output",
        metavar="PATH",
        help=(
            "Output file (single input only) or output directory (multiple inputs). "
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

    if len(args.chain_from) != 1 or len(args.chain_to) != 1:
        print("Error: --from and --to must be single-character chain IDs.", file=sys.stderr)
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
    if len(inputs) > 1:
        if not args.output:
            print(
                "Error: -o/--output directory is required for multiple input files.",
                file=sys.stderr,
            )
            sys.exit(2)
        out_dir_multi = os.path.abspath(args.output)
        if os.path.isfile(out_dir_multi):
            print("Error: for multiple inputs, -o must be a directory.", file=sys.stderr)
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
            dest = _default_out_path(in_path, out_dir_multi, args.suffix)

        if args.dry_run:
            print(f"{in_path} -> {dest}")
            ok_n += 1
            continue

        ok, data = process_file(
            in_path,
            args.chain_from,
            args.chain_to,
            args.merge_renumber,
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
