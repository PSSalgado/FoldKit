#!/usr/bin/env python3
"""
Open multiple structure files in Coot and set representation to C-alpha.

Two ways to select files:

  1) Directory mode: pass one or more directories. Collects *.pdb, *.cif, *.mmcif
     in each directory (non-recursive). Optional --filter SUBSTRING on basename.

  2) Pattern mode: pass glob patterns (and optional directories to search).
     Recursive search under each directory; supports PDB and mmCIF.

Examples:
  python utils/open_models_in_coot.py models/
  python utils/open_models_in_coot.py --filter set_a dir1/ dir2/
  python utils/open_models_in_coot.py '*tag*LSQaligned2*ref*pdb' models/
"""

from __future__ import annotations

import argparse
import glob
import os
import sys
from pathlib import Path
from typing import List, Optional, Set

import fnmatch

_REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if _REPO_ROOT not in sys.path:
    sys.path.insert(0, _REPO_ROOT)

from utils.cli_log import add_log_args, setup_log_from_args


def is_structure_file(file_path: Path) -> bool:
    valid_extensions = {".pdb", ".cif", ".mmcif"}
    return file_path.suffix.lower() in valid_extensions


def collect_from_directories(
    directories: List[str],
    filter_substring: Optional[str],
) -> List[str]:
    """Top-level PDB/mmCIF in each directory (non-recursive), optional basename filter."""
    found: Set[str] = set()
    for directory in directories:
        d = os.path.abspath(directory)
        if not os.path.isdir(d):
            print(f"Warning: Directory '{directory}' does not exist. Skipping.", file=sys.stderr)
            continue
        for ext in ("*.pdb", "*.cif", "*.mmcif"):
            for f in glob.glob(os.path.join(d, ext)):
                if filter_substring is None or filter_substring in os.path.basename(f):
                    found.add(os.path.abspath(f))
    return sorted(found)


def _pattern_has_glob_metachars(pattern: str) -> bool:
    return any(ch in pattern for ch in "*?[")


def find_matching_files(patterns: List[str], search_dirs: Optional[List[str]] = None) -> List[str]:
    """Find files matching patterns under search_dirs (recursive rglob).

    Glob patterns use shell-style matching (fnmatch) on each file basename.
    """
    matching_files: Set[str] = set()

    if not search_dirs:
        search_dirs = ["."]

    for search_dir in search_dirs:
        search_path = Path(search_dir)
        if not search_path.exists():
            print(f"Warning: Directory '{search_dir}' does not exist. Skipping.", file=sys.stderr)
            continue

        for pattern in patterns:
            if _pattern_has_glob_metachars(pattern):
                for file_path in search_path.rglob("*"):
                    if not file_path.is_file() or not is_structure_file(file_path):
                        continue
                    file_str = str(file_path.resolve())
                    basename = file_path.name
                    # Match basename only: fnmatch on full paths fails for patterns like *.pdb
                    # because * does not span '/' in shell-style globs.
                    if fnmatch.fnmatch(basename, pattern):
                        print(f"Found matching file: {basename}")
                        matching_files.add(file_str)
            else:
                file_path = search_path / pattern
                if file_path.is_file() and is_structure_file(file_path):
                    matching_files.add(str(file_path.resolve()))

    return sorted(matching_files)


def create_coot_script(model_files: List[str]) -> str:
    script_content = """
import os
import sys

set_nomenclature_errors_on_read("ignore")
set_show_symmetry_master(0)

for model_path in {0}:
    handle_read_draw_molecule_with_recentre(model_path, 0)
    mol_index = graphics_n_molecules() - 1
    graphics_to_ca_representation(mol_index)
    set_molecule_bonds_colour_map_rotation(mol_index, 20 * mol_index)
""".format(
        model_files,
    )
    return script_content


def main() -> None:
    ap = argparse.ArgumentParser(
        description="Open PDB/mmCIF models in Coot with CA representation.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples (from repository root):
  python utils/open_models_in_coot.py models/
  python utils/open_models_in_coot.py --filter set_a dir1/ dir2/
  python utils/open_models_in_coot.py '*model_*.cif' /path/to/search
""",
    )
    ap.add_argument(
        "--filter",
        metavar="SUBSTRING",
        help="Directory mode only: keep only files whose basename contains this substring.",
    )
    ap.add_argument(
        "paths",
        nargs="+",
        metavar="PATH",
        help="Directories (all top-level structures) and/or glob patterns.",
    )
    add_log_args(ap)
    args = ap.parse_args()
    setup_log_from_args(args, script_path=__file__, inputs=list(getattr(args, "paths", []) or []), pattern=None)

    patterns: List[str] = []
    directories: List[str] = []
    for arg in args.paths:
        if os.path.isdir(arg):
            directories.append(arg)
        else:
            patterns.append(arg)

    if not patterns and directories:
        model_files = collect_from_directories(directories, args.filter)
        if not model_files:
            flt = f" with filter '{args.filter}'" if args.filter else ""
            print(
                f"Error: No PDB/mmCIF files found in the given directories{flt}.",
                file=sys.stderr,
            )
            sys.exit(1)
        if args.filter:
            print(f"Applied filter '{args.filter}': {len(model_files)} file(s).")
    elif patterns:
        if args.filter:
            print(
                "Note: --filter applies only to directory-only mode; ignored for pattern mode.",
                file=sys.stderr,
            )
        model_files = find_matching_files(patterns, directories if directories else None)
        if not model_files:
            print(
                "Error: No valid structure files (PDB/mmCIF) found for the given patterns.",
                file=sys.stderr,
            )
            sys.exit(1)
    else:
        print("Error: No directories or patterns to process.", file=sys.stderr)
        sys.exit(1)

    print(f"Found {len(model_files)} structure file(s) to open:")
    for f in model_files:
        print(f"  - {os.path.basename(f)}")

    script_file = "temp_coot_script.py"
    with open(script_file, "w") as f:
        f.write(create_coot_script(model_files))

    print("Starting Coot...")
    os.system(f"coot --script {script_file}")

    if os.path.exists(script_file):
        os.remove(script_file)


if __name__ == "__main__":
    main()
