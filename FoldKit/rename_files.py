"""
Rename files (and optionally directories) under a path.

Legacy mode: strip embedded date stamps from PDB names matching
  {prefix}_YYYY_MM_DD_HH_MM_*.pdb  in the given directory only (non-recursive).

Regex mode: --remove and/or --replace/--with, optional recursion.
"""

from __future__ import annotations

import argparse
import glob
import os
import re
import sys


def strip_date_from_prefixed_pdbs(directory: str, name_pattern: str) -> None:
    """Rename PDBs in one directory: drop YYYY_MM_DD_HH_MM segment after prefix."""
    if not os.path.exists(directory):
        print(f"Error: Directory '{directory}' does not exist.")
        return

    pattern = os.path.join(directory, f"{name_pattern}_*_*.pdb")
    files = glob.glob(pattern)

    if not files:
        print(f"No files matching pattern '{name_pattern}_*_*.pdb' found in '{directory}'")
        return

    print(f"Found {len(files)} files to rename:")

    for old_path in files:
        old_name = os.path.basename(old_path)
        match = re.match(
            rf"({name_pattern})_\d{{4}}_\d{{2}}_\d{{2}}_\d{{2}}_\d{{2}}_(.*\.pdb)",
            old_name,
        )
        if match:
            prefix = match.group(1)
            rest = match.group(2)
            new_name = f"{prefix}_{rest}"
            new_path = os.path.join(directory, new_name)
            try:
                os.rename(old_path, new_path)
                print(f"Renamed: {old_name} -> {new_name}")
            except OSError as e:
                print(f"Error renaming {old_name}: {e}")
        else:
            print(f"Skipping {old_name}: Does not match expected pattern")


def rename_item(
    old_path: str,
    pattern_to_remove=None,
    pattern_to_replace=None,
    replacement=None,
):
    """Rename a file or directory by removing or replacing patterns in its basename."""
    parent_dir = os.path.dirname(old_path)
    old_name = os.path.basename(old_path)
    new_name = old_name

    if pattern_to_remove:
        new_name = re.sub(pattern_to_remove, "", new_name)

    if pattern_to_replace and replacement is not None:
        new_name = re.sub(pattern_to_replace, replacement, new_name)

    if new_name == old_name:
        return None

    new_path = os.path.join(parent_dir, new_name)
    try:
        os.rename(old_path, new_path)
        print(f"Renamed: {old_name} -> {new_name}")
        return new_path
    except OSError as e:
        print(f"Error renaming {old_name}: {e}")
        return None


def process_directory(
    directory: str,
    pattern_to_remove=None,
    pattern_to_replace=None,
    replacement=None,
    recursive: bool = True,
    rename_dirs: bool = True,
    root_anchor: str | None = None,
):
    """Process a directory: rename files, optionally recurse and rename dirs."""
    if root_anchor is None:
        root_anchor = os.path.abspath(directory)

    if not os.path.exists(directory):
        print(f"Error: Directory '{directory}' does not exist.")
        return

    print(f"Processing directory: {directory}")

    items = os.listdir(directory)

    for item in items:
        item_path = os.path.join(directory, item)
        if os.path.isdir(item_path):
            continue
        rename_item(item_path, pattern_to_remove, pattern_to_replace, replacement)

    if recursive:
        items = os.listdir(directory)
        for item in items:
            item_path = os.path.join(directory, item)
            if os.path.isdir(item_path):
                process_directory(
                    item_path,
                    pattern_to_remove,
                    pattern_to_replace,
                    replacement,
                    recursive,
                    rename_dirs,
                    root_anchor,
                )

    if rename_dirs and recursive:
        for item in os.listdir(directory):
            item_path = os.path.join(directory, item)
            if os.path.isdir(item_path):
                rename_item(item_path, pattern_to_remove, pattern_to_replace, replacement)

    if rename_dirs and os.path.abspath(directory) != root_anchor:
        rename_item(directory, pattern_to_remove, pattern_to_replace, replacement)


def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Rename files under a directory. "
            "Legacy: pass DIRECTORY and PREFIX to strip date segments from PREFIX_*_*.pdb. "
            "Or use --remove / --replace for regex-based renames (recursive by default)."
        )
    )
    parser.add_argument(
        "directory",
        help="Directory to operate in.",
    )
    parser.add_argument(
        "strip_prefix",
        nargs="?",
        default=None,
        metavar="PREFIX",
        help=(
            "Legacy mode only: filename prefix; renames {PREFIX}_YYYY_MM_DD_HH_MM_*.pdb "
            "to {PREFIX}_*.pdb in DIRECTORY (non-recursive)."
        ),
    )
    parser.add_argument(
        "--remove",
        metavar="PATTERN",
        help="Regex to remove from each basename.",
    )
    parser.add_argument(
        "--replace",
        metavar="PATTERN",
        help="Regex to replace in each basename.",
    )
    parser.add_argument(
        "--with",
        dest="replacement",
        metavar="STRING",
        help="Replacement string (required with --replace).",
    )
    parser.add_argument(
        "--no-recursive",
        action="store_true",
        help="Do not enter subdirectories.",
    )
    parser.add_argument(
        "--files-only",
        action="store_true",
        help="Only rename files, not directories.",
    )

    args = parser.parse_args()
    directory = os.path.abspath(args.directory)

    has_regex = args.remove or args.replace
    if args.strip_prefix is not None and has_regex:
        parser.error("Use either legacy PREFIX (positional) or --remove/--replace, not both")

    if args.strip_prefix is not None:
        strip_date_from_prefixed_pdbs(directory, args.strip_prefix)
        return

    if args.replace and args.replacement is None:
        parser.error("When using --replace, you must also specify --with")

    if not has_regex:
        parser.error(
            "Specify legacy PREFIX as second argument (date-strip PDBs), "
            "or use --remove and/or --replace for regex renaming"
        )

    process_directory(
        directory,
        pattern_to_remove=args.remove,
        pattern_to_replace=args.replace,
        replacement=args.replacement,
        recursive=not args.no_recursive,
        rename_dirs=not args.files_only,
        root_anchor=directory,
    )


if __name__ == "__main__":
    main()
