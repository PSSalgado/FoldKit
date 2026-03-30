#!/usr/bin/env python3
"""
Batch driver for superimpose_coot_LSQ.py --pattern across multiple conditions and tags.

Layout (same spirit as superimpose_coot_LSQ.py pattern-mode examples):

  refs_base/
    <condition>/          e.g. condition_1, condition_2
      <reference_subdir>/ trimmed references (LSQ output naming)
  models_base/
    <condition>/          same condition labels as under refs_base
      *_fl_<set>*/         model subdirectories per set (set_a, set_b, …)

Override paths with --ref-base and --models-base, or edit the defaults below.
"""

import argparse
import os
import subprocess
import sys
from pathlib import Path

# Defaults: set via --ref-base / --models-base or edit here (see module docstring).
REF_STRUCTURES_BASE = "/path/to/refs"
MODELS_BASE = "/path/to/models"

# Subdirectory under each condition that holds reference PDBs (LSQ-style tree).
REFERENCE_SUBDIR = "LSQaligned2_reference_m0"

# Stem used in reference filenames (matches *LSQ*reference*.pdb style in --pattern mode).
REFERENCE_STEM = "reference"

# Condition directory name -> filename prefix token (e.g. run1 in run1_sd_set_a).
CONDITION_PREFIX = {
    "condition_1": "run1",
    "condition_2": "run2",
    "condition_3": "run3",
}

CONDITIONS = list(CONDITION_PREFIX.keys())

# Tags to process if --tags is omitted (same idea as --filter=set_a in SSM/LSQ scripts).
TAGS: list[str] = []


def get_script_path() -> str:
    return os.path.join(os.path.dirname(os.path.abspath(__file__)), "superimpose_coot_LSQ.py")


def run_superposition(
    condition: str,
    tag: str,
    script_path: str,
    ref_base: str | None = None,
    models_base: str | None = None,
) -> bool:
    """
    Run superimpose_coot_LSQ.py --pattern for one condition and tag.

    Args:
        condition: Subdirectory name under ref_base and models_base (e.g. condition_1).
        tag: Model tag, like set_a in --filter=set_a (matches *_fl_<tag>* dirs and ref glob).
        script_path: Path to superimpose_coot_LSQ.py.
        ref_base: Base directory for reference trees.
        models_base: Base directory for model trees.
    """
    ref_base = ref_base or REF_STRUCTURES_BASE
    models_base = models_base or MODELS_BASE
    prefix = CONDITION_PREFIX[condition]

    ref_dir = os.path.join(ref_base, condition, REFERENCE_SUBDIR)
    model_dir = os.path.join(models_base, condition)

    ref_pattern = f"{prefix}_sd_{tag}"
    target_pattern = f"{prefix}_sd_{REFERENCE_STEM}"
    model_pattern = f"_fl_{tag}"

    print(f"  Running superposition for {condition} — tag {tag}")
    print(f"  Reference dir: {ref_dir}")
    print(f"  Model dir: {model_dir}")
    print(f"  Reference pattern: {ref_pattern}")
    print(f"  Target pattern: {target_pattern}")
    print(f"  Model pattern: {model_pattern}")

    if not os.path.isdir(ref_dir):
        print(f"  ERROR: Reference directory does not exist: {ref_dir}")
        return False

    ref_glob = f"*{tag}*LSQ*{REFERENCE_STEM}*.pdb"
    ref_files = list(Path(ref_dir).glob(ref_glob))
    if not ref_files:
        print(f"  WARNING: No reference files matching {ref_glob!r}")
        print("  Skipping this tag/condition combination.")
        return False

    if not os.path.isdir(model_dir):
        print(f"  ERROR: Model directory does not exist: {model_dir}")
        return False

    model_dirs = list(Path(model_dir).glob(f"*_fl_{tag}*"))
    if not model_dirs:
        print(f"  WARNING: No model directories matching *_fl_{tag}*")
        print("  Skipping this tag/condition combination.")
        return False

    original_dir = os.getcwd()
    os.chdir(model_dir)

    cmd = [
        sys.executable,
        script_path,
        "--pattern",
        "--divider=LSQ_",
        f"--ref-file-pattern=*{tag}*LSQ*{REFERENCE_STEM}*.pdb",
        "--model-file-pattern=*.cif",
        ref_dir,
        model_dir,
        ref_pattern,
        model_pattern,
        target_pattern,
    ]

    print(f"  Running: {' '.join(cmd)}")

    try:
        subprocess.run(cmd, check=True, capture_output=True, text=True)
        print(f"  OK: completed {condition} — {tag}")
        success = True
    except subprocess.CalledProcessError as e:
        print(f"  ERROR: {condition} — {tag}")
        print(f"  {e}")
        print(f"  stdout: {e.stdout}")
        print(f"  stderr: {e.stderr}")
        success = False

    os.chdir(original_dir)
    return success


def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Run superimpose_coot_LSQ.py --pattern across conditions and tags "
            "(reference + models layout; same naming style as LSQ/SSM examples)."
        )
    )
    parser.add_argument(
        "--tags",
        nargs="+",
        metavar="TAG",
        help="Tags to process (e.g. set_a set_b); same role as --filter in SSM/LSQ scripts.",
    )
    parser.add_argument(
        "--conditions",
        nargs="+",
        metavar="NAME",
        help=f"Condition subdirs to process (default: {CONDITIONS})",
    )
    parser.add_argument(
        "--ref-base",
        dest="ref_base",
        help=f"Base directory for reference trees (default: {REF_STRUCTURES_BASE})",
    )
    parser.add_argument(
        "--models-base",
        dest="models_base",
        help=f"Base directory for model trees (default: {MODELS_BASE})",
    )
    args = parser.parse_args()

    ref_base = args.ref_base or REF_STRUCTURES_BASE
    models_base = args.models_base or MODELS_BASE
    tags = args.tags or TAGS
    conditions = args.conditions or CONDITIONS

    if not tags:
        parser.error(
            "No tags given. Use --tags (e.g. set_a set_b) or set TAGS in this script."
        )

    unknown_c = [c for c in conditions if c not in CONDITION_PREFIX]
    if unknown_c:
        parser.error(f"Unknown condition(s): {unknown_c!r}. Known: {list(CONDITION_PREFIX)}")

    script_path = get_script_path()

    for tag in tags:
        print("=======================================================================")
        print(f"Tag: {tag}")
        print("=======================================================================")

        for condition in conditions:
            print("-----------------------------------------------------------------------")
            print(f"  Condition: {condition}")
            print("-----------------------------------------------------------------------")
            run_superposition(condition, tag, script_path, ref_base=ref_base, models_base=models_base)
            print()

    print("Done.")
    print(
        "Outputs are written under each models_base/<condition>/ tree "
        "(superimpose_coot_LSQ.py --pattern convention)."
    )


if __name__ == "__main__":
    main()
