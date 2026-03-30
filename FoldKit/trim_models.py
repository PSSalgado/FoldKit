#!/usr/bin/env python3
"""
Trim PDB/mmCIF models to a common residue span (shortest structure in each batch).

Standalone script: no imports from other FoldKit modules. Optional external tool: ``gemmi``
for mmCIF → PDB conversion (same as other FoldKit structure scripts).

``trim_superimposeLSQ.py`` imports ``run_trim_job`` and helpers from this module for ``--trim`` / ``--trim-only``.
"""
from __future__ import annotations

import fnmatch
import glob
import os
import subprocess
import sys


def pattern_has_glob_chars(pattern: str) -> bool:
    return any(ch in pattern for ch in ["*", "?", "["])


def matches_filter_ci(basename: str, pattern: str) -> bool:
    """
    Case-insensitive filter match on basename (and stem).
    Glob chars (* ? [) → fnmatch; else substring match.
    """
    b = basename.lower()
    p = pattern.lower()
    stem, _ = os.path.splitext(b)
    if pattern_has_glob_chars(p):
        return fnmatch.fnmatchcase(b, p) or fnmatch.fnmatchcase(stem, p)
    return (p in b) or (p in stem)


def get_residue_range(file_path: str):
    """Return (min_res, max_res) from ATOM/HETATM records, or None."""
    min_res = float("inf")
    max_res = float("-inf")
    is_cif = file_path.lower().endswith(".cif")

    with open(file_path, "r") as f:
        for line in f:
            if is_cif:
                if line.strip().startswith("ATOM") or line.strip().startswith("HETATM"):
                    parts = line.split()
                    if len(parts) > 8:
                        try:
                            res_num = int(parts[8])
                            min_res = min(min_res, res_num)
                            max_res = max(max_res, res_num)
                        except (ValueError, IndexError):
                            continue
            else:
                if line.startswith("ATOM  ") or line.startswith("HETATM"):
                    try:
                        res_num = int(line[22:26].strip())
                        min_res = min(min_res, res_num)
                        max_res = max(max_res, res_num)
                    except ValueError:
                        continue

    if min_res == float("inf") or max_res == float("-inf"):
        return None
    return (min_res, max_res)


def trim_pdb_file(input_file: str, output_file: str, start_res: int, end_res: int) -> None:
    """Write a PDB with only residues in [start_res, end_res] (ATOM/HETATM)."""
    with open(input_file, "r") as infile, open(output_file, "w") as outfile:
        for line in infile:
            if line.startswith("ATOM  ") or line.startswith("HETATM"):
                try:
                    res_num = int(line[22:26].strip())
                    if start_res <= res_num <= end_res:
                        outfile.write(line)
                except ValueError:
                    continue
            elif line.startswith("TER") or line.startswith("END"):
                outfile.write(line)


def convert_cif_to_pdb(cif_file: str, pdb_file: str) -> bool:
    try:
        subprocess.run(
            ["gemmi", "convert", "--to=pdb", cif_file, pdb_file],
            check=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )
        return True
    except (subprocess.SubprocessError, FileNotFoundError):
        print(
            f"Warning: Could not convert {cif_file} to PDB format. Gemmi may not be installed."
        )
        return False


def generate_trimmed_filename(original_path: str, pattern: str, source_dir: str) -> str:
    original_filename = os.path.basename(original_path)
    dir_name = os.path.basename(source_dir)
    dir_name = dir_name.replace("models", "").replace("Models", "")
    if dir_name.endswith("_"):
        dir_name = dir_name[:-1]
    if dir_name.startswith("_"):
        dir_name = dir_name[1:]

    pattern_pos = original_filename.lower().find(pattern.lower())
    if pattern_pos == -1:
        suffix = original_filename
    else:
        suffix = original_filename[pattern_pos:]

    if suffix.lower().endswith(".pdb") or suffix.lower().endswith(".cif"):
        suffix = os.path.splitext(suffix)[0]

    new_filename = f"SDt_{dir_name}_{suffix}.pdb"
    return new_filename.replace("__", "_")


def collect_models_from_directories(directories: list[str]) -> list[tuple[str, str]]:
    """Return list of (absolute_path, source_directory) for top-level *.pdb / *.cif per dir."""
    all_model_files: list[tuple[str, str]] = []
    for directory in directories:
        if not os.path.exists(directory):
            print(f"Error: Directory '{directory}' does not exist. Skipping.")
            continue
        dir_models: list[str] = []
        dir_models.extend(glob.glob(os.path.join(directory, "*.pdb")))
        dir_models.extend(glob.glob(os.path.join(directory, "*.cif")))
        if not dir_models:
            print(f"Warning: No PDB or CIF files found in '{directory}'. Skipping.")
            continue
        for model in dir_models:
            all_model_files.append((model, directory))
    return all_model_files


def run_trim_job(
    all_model_files: list[tuple[str, str]],
    filter_patterns: list[str],
) -> list[str]:
    """
    Trim models to the shortest residue span among each batch.
    If filter_patterns is non-empty, one batch per pattern (output under trimmed_<pattern>/).
    If empty, single batch under trimmed_all/.
    Returns list of written output PDB paths.
    """
    trimmed_files: list[str] = []

    if filter_patterns:
        for pattern in filter_patterns:
            filtered_models = [
                (model, src_dir)
                for model, src_dir in all_model_files
                if matches_filter_ci(os.path.basename(model), pattern)
            ]
            if not filtered_models:
                print(f"Warning: No models match the filter pattern '{pattern}'. Skipping.")
                continue

            print(f"\nProcessing filter pattern: '{pattern}'")
            print(f"Found {len(filtered_models)} matching models.")

            pattern_dir = f"trimmed_{pattern}"
            if not os.path.exists(pattern_dir):
                os.makedirs(pattern_dir)

            model_ranges = {}
            for model, _ in filtered_models:
                res_range = get_residue_range(model)
                if res_range:
                    model_ranges[model] = res_range
                    print(f"Model {os.path.basename(model)}: Residues {res_range[0]}-{res_range[1]}")
                else:
                    print(
                        f"Warning: Could not determine residue range for {os.path.basename(model)}. Skipping."
                    )

            if not model_ranges:
                print(f"Error: Could not determine residue range for any models with pattern '{pattern}'.")
                continue

            shortest_model = min(model_ranges.items(), key=lambda x: x[1][1] - x[1][0] + 1)
            shortest_range = shortest_model[1]
            shortest_model_path = shortest_model[0]
            print(
                f"\nShortest model for pattern '{pattern}': {os.path.basename(shortest_model_path)} "
                f"with residues {shortest_range[0]}-{shortest_range[1]}"
            )
            print(
                f"\nTrimming models with pattern '{pattern}' to residue range "
                f"{shortest_range[0]}-{shortest_range[1]}..."
            )

            for model, src_dir in filtered_models:
                new_filename = generate_trimmed_filename(model, pattern, src_dir)
                output_file = os.path.join(pattern_dir, new_filename)

                if model.lower().endswith(".cif"):
                    temp_pdb = os.path.join(
                        pattern_dir, f"temp_{os.path.splitext(os.path.basename(model))[0]}.pdb"
                    )
                    if convert_cif_to_pdb(model, temp_pdb):
                        trim_pdb_file(temp_pdb, output_file, shortest_range[0], shortest_range[1])
                        os.remove(temp_pdb)
                        trimmed_files.append(output_file)
                        print(f"Converted and trimmed {os.path.basename(model)} -> {os.path.basename(output_file)}")
                    else:
                        print(f"Error: Could not process CIF file {os.path.basename(model)}")
                else:
                    trim_pdb_file(model, output_file, shortest_range[0], shortest_range[1])
                    trimmed_files.append(output_file)
                    print(f"Trimmed {os.path.basename(model)} -> {os.path.basename(output_file)}")

            shortest_src_dir = next(
                sd for m, sd in filtered_models if m == shortest_model_path
            )
            shortest_new_filename = generate_trimmed_filename(
                shortest_model_path, pattern, shortest_src_dir
            )
            shortest_output_file = os.path.join(pattern_dir, shortest_new_filename)
            if not os.path.exists(shortest_output_file):
                trim_pdb_file(shortest_model_path, shortest_output_file, shortest_range[0], shortest_range[1])
                trimmed_files.append(shortest_output_file)
                print(
                    f"Trimmed shortest model: {os.path.basename(shortest_model_path)} -> "
                    f"{os.path.basename(shortest_output_file)}"
                )
    else:
        trim_dir = "trimmed_all"
        if not os.path.exists(trim_dir):
            os.makedirs(trim_dir)

        model_ranges = {}
        for model, _ in all_model_files:
            res_range = get_residue_range(model)
            if res_range:
                model_ranges[model] = res_range
                print(f"Model {os.path.basename(model)}: Residues {res_range[0]}-{res_range[1]}")
            else:
                print(f"Warning: Could not determine residue range for {os.path.basename(model)}. Skipping.")

        if not model_ranges:
            print("Error: Could not determine residue range for any models.")
            sys.exit(1)

        shortest_model = min(model_ranges.items(), key=lambda x: x[1][1] - x[1][0] + 1)
        shortest_range = shortest_model[1]
        shortest_model_path = shortest_model[0]
        print(
            f"\nShortest model: {os.path.basename(shortest_model_path)} with residues "
            f"{shortest_range[0]}-{shortest_range[1]}"
        )

        shortest_src_dir = next(sd for m, sd in all_model_files if m == shortest_model_path)
        print(f"\nTrimming all models to residue range {shortest_range[0]}-{shortest_range[1]}...")

        for model, src_dir in all_model_files:
            new_filename = f"SDt_{os.path.basename(src_dir)}_{os.path.basename(model)}"
            new_filename = new_filename.replace("models", "").replace("Models", "")
            new_filename = new_filename.replace("__", "_")
            if new_filename.endswith(".cif"):
                new_filename = new_filename[:-4] + ".pdb"
            output_file = os.path.join(trim_dir, new_filename)

            if model.lower().endswith(".cif"):
                temp_pdb = os.path.join(trim_dir, f"temp_{os.path.splitext(os.path.basename(model))[0]}.pdb")
                if convert_cif_to_pdb(model, temp_pdb):
                    trim_pdb_file(temp_pdb, output_file, shortest_range[0], shortest_range[1])
                    os.remove(temp_pdb)
                    trimmed_files.append(output_file)
                    print(f"Converted and trimmed {os.path.basename(model)} -> {os.path.basename(output_file)}")
                else:
                    print(f"Error: Could not process CIF file {os.path.basename(model)}")
            else:
                trim_pdb_file(model, output_file, shortest_range[0], shortest_range[1])
                trimmed_files.append(output_file)
                print(f"Trimmed {os.path.basename(model)} -> {os.path.basename(output_file)}")

        shortest_new_filename = f"SDt_{os.path.basename(shortest_src_dir)}_{os.path.basename(shortest_model_path)}"
        shortest_new_filename = shortest_new_filename.replace("models", "").replace("Models", "")
        shortest_new_filename = shortest_new_filename.replace("__", "_")
        if shortest_new_filename.endswith(".cif"):
            shortest_new_filename = shortest_new_filename[:-4] + ".pdb"
        shortest_output_file = os.path.join(trim_dir, shortest_new_filename)
        if not os.path.exists(shortest_output_file):
            trim_pdb_file(shortest_model_path, shortest_output_file, shortest_range[0], shortest_range[1])
            trimmed_files.append(shortest_output_file)
            print(
                f"Trimmed shortest model: {os.path.basename(shortest_model_path)} -> "
                f"{os.path.basename(shortest_output_file)}"
            )

    return trimmed_files


def main() -> None:
    if len(sys.argv) < 2:
        print(
            "Usage: python trim_models.py [--filter=set_a,set_b,...] directory1 [directory2 ...]"
        )
        print("  --filter=   Optional comma-separated patterns (basename substring or glob).")
        print("              If omitted, all PDB/CIF files in the given directories are used.")
        sys.exit(1)

    filter_patterns: list[str] = []
    directories: list[str] = []
    for arg in sys.argv[1:]:
        if arg.startswith("--filter="):
            filter_patterns = [p for p in arg.split("=", 1)[1].split(",") if p]
        elif not arg.startswith("--"):
            directories.append(os.path.abspath(arg))

    if not directories:
        print("Error: No input directories specified.")
        sys.exit(1)

    all_model_files = collect_models_from_directories(directories)
    if not all_model_files:
        print("Error: No PDB or CIF files found in any of the specified directories.")
        sys.exit(1)

    out = run_trim_job(all_model_files, filter_patterns)
    if not out:
        print("Error: No trimmed outputs were produced.")
        sys.exit(1)
    print(f"\nDone. Wrote {len(out)} structure file(s).")


if __name__ == "__main__":
    main()
