"""
Shared reference/model filename pairing for superimpose_coot_LSQ.py and superimpose_coot_SSM.py
--pattern mode. Keeps SSM independent of the LSQ script.
"""

import glob
import os
import re


def _normalize_ref_model_pattern(pattern):
    """Normalise a pattern by removing special characters, converting to lowercase."""
    if pattern is None:
        return None
    return re.sub(r"[_\-\s]", "", pattern.lower())


def find_ref_model_matches(
    ref_dir,
    model_dir,
    ref_pattern,
    model_pattern,
    target_pattern=None,
    divider=None,
    strict_position=False,
    ref_file_pattern="*.pdb",
    model_file_pattern="*.cif",
):
    """Find matching reference and model files based on filename / subdirectory patterns."""
    print(f"\nLooking for reference files in: {ref_dir}")
    print(f"Using reference file pattern: {ref_file_pattern}")
    print(f"Using reference search pattern: '{ref_pattern}'")
    if target_pattern:
        print(f"Using target search pattern: '{target_pattern}'")
    if divider:
        print(f"Using divider: '{divider}'")

    norm_ref_pattern = _normalize_ref_model_pattern(ref_pattern)
    norm_model_pattern = _normalize_ref_model_pattern(model_pattern)

    ref_file_glob = os.path.join(ref_dir, ref_file_pattern)
    all_ref_files = glob.glob(ref_file_glob)
    print(f"Found {len(all_ref_files)} total files matching '{ref_file_pattern}' in reference directory")

    ref_files = []
    for ref_file in all_ref_files:
        filename = os.path.basename(ref_file)
        norm_filename = _normalize_ref_model_pattern(filename)

        if divider is None:
            if norm_ref_pattern in norm_filename:
                if target_pattern is None or _normalize_ref_model_pattern(target_pattern) in norm_filename:
                    ref_files.append(ref_file)
                    print(f"Matched reference: {filename}")
        elif divider in filename:
            norm_divider = _normalize_ref_model_pattern(divider)
            if strict_position and target_pattern is not None:
                parts = filename.split(divider)
                if len(parts) != 2:
                    continue
                first_part = _normalize_ref_model_pattern(parts[0])
                second_part = _normalize_ref_model_pattern(parts[1])
                if norm_ref_pattern in first_part and _normalize_ref_model_pattern(target_pattern) in second_part:
                    ref_files.append(ref_file)
                    print(f"Matched reference: {filename}")
            else:
                if norm_divider in norm_filename and norm_ref_pattern in norm_filename:
                    if target_pattern is None or _normalize_ref_model_pattern(target_pattern) in norm_filename:
                        ref_files.append(ref_file)
                        print(f"Matched reference: {filename}")

    if not ref_files:
        if divider is not None and strict_position and target_pattern is not None:
            print(
                f"No reference files found matching pattern '{ref_pattern}' before '{divider}' and '{target_pattern}' after"
            )
        elif divider is not None and target_pattern is not None:
            print(
                f"No reference files found containing both patterns '{ref_pattern}' and '{target_pattern}' with divider '{divider}'"
            )
        else:
            print(f"No reference files found containing pattern '{ref_pattern}'")
            if target_pattern is not None:
                print(f"and '{target_pattern}'")
        return []

    print(f"Found {len(ref_files)} matching reference files")

    print(f"\nLooking for model files in: {model_dir}")
    print(f"Using model file pattern: {model_file_pattern}")
    print(f"Using model search pattern: '{model_pattern}'")

    if not os.path.exists(model_dir):
        print(f"Error: Model directory '{model_dir}' does not exist")
        return []

    print("Contents of model directory:")
    try:
        for item in os.listdir(model_dir):
            item_path = os.path.join(model_dir, item)
            type_str = "DIR" if os.path.isdir(item_path) else "FILE"
            print(f"  {type_str} - {item}")
    except Exception as e:
        print(f"Error listing model directory: {e}")

    subdirs = [
        item
        for item in os.listdir(model_dir)
        if os.path.isdir(os.path.join(model_dir, item))
    ]
    model_dir_has_subdirs = len(subdirs) > 0
    print(f"Model directory has {len(subdirs)} subdirectories")

    if model_dir_has_subdirs:
        matching_subdirs = []
        for subdir in subdirs:
            subdir_path = os.path.join(model_dir, subdir)
            if norm_model_pattern in _normalize_ref_model_pattern(subdir):
                matching_subdirs.append(subdir_path)
                print(f"Matched subdirectory by name: {subdir}")
                continue
            has_matching_files = False
            for file_path in glob.glob(os.path.join(subdir_path, model_file_pattern)):
                if norm_model_pattern in _normalize_ref_model_pattern(os.path.basename(file_path)):
                    has_matching_files = True
                    break
            if has_matching_files:
                matching_subdirs.append(subdir_path)
                print(f"Matched subdirectory by file content: {subdir}")

        if not matching_subdirs:
            print(f"No subdirectories found matching pattern: {model_pattern}")
            return []

        print(f"Found {len(matching_subdirs)} matching subdirectories")
        model_files_by_subdir = {}
        for subdir in matching_subdirs:
            mf = glob.glob(os.path.join(subdir, model_file_pattern))
            if mf:
                model_files_by_subdir[subdir] = mf
                print(f"Found {len(mf)} model files in {os.path.basename(subdir)}")
        if not model_files_by_subdir:
            print("No model files found in matching subdirectories")
            return []

        matches = []
        for ref_file in ref_files:
            for subdir, model_files in model_files_by_subdir.items():
                matches.append((ref_file, subdir, model_files))
        return matches

    model_file_glob = os.path.join(model_dir, model_file_pattern)
    flat_models = glob.glob(model_file_glob)
    print(f"Found {len(flat_models)} total files matching '{model_file_pattern}' in model directory")

    matching_models = []
    for model_file in flat_models:
        filename = os.path.basename(model_file)
        if norm_model_pattern in _normalize_ref_model_pattern(filename):
            matching_models.append(model_file)
            print(f"Matched model file: {filename}")

    if not matching_models:
        print(f"No model files found matching pattern: {model_pattern}")
        return []

    print(f"Found {len(matching_models)} matching model files")
    matches = []
    for ref_file in ref_files:
        matches.append((ref_file, model_dir, matching_models))
    return matches
