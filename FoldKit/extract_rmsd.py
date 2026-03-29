import os
import sys
import argparse
import re


def _peek_log_head(path: str, max_bytes: int = 65536) -> str:
    try:
        with open(path, "rb") as f:
            return f.read(max_bytes).decode("utf-8", errors="replace")
    except OSError:
        return ""


def detect_rmsd_log_format(log_file: str) -> str:
    """
    Return 'ssm' if the log looks like Coot SSM (Superposing … onto …), else 'lsq'.
    """
    sample = _peek_log_head(log_file)
    if "Superposing" in sample and " onto " in sample:
        return "ssm"
    return "lsq"


def extract_ssm_rmsd_values(log_file: str) -> None:
    """
    Extract molecule/superposition info and RMSD values from a Coot log for SSM.
    SSM logs use different wording than LSQ: scripts print "Superposing X onto …",
    and Coot may print "Aligning … to …". We capture whichever alignment line appears
    before each "INFO: core rmsd" block and write it together with the RMSD block.
    """
    if not os.path.exists(log_file):
        print(f"Error: Log file '{log_file}' not found.")
        return

    output_dir = os.path.dirname(log_file)
    base = os.path.basename(log_file)
    if base == "coot_log.txt":
        rmsd_basename = "rmsd_SSM_values.txt"
    elif base.startswith("coot_log_") and base.endswith(".txt"):
        suffix = base[9:-4]
        rmsd_basename = "rmsd_SSM_values_{}.txt".format(suffix)
    else:
        rmsd_basename = "rmsd_SSM_values.txt"
    rmsd_file = os.path.join(output_dir, rmsd_basename)

    print(f"Extracting RMSD values (SSM) from {log_file}")
    print(f"Writing to {rmsd_file}")

    current_block = []
    last_alignment_line = None
    with open(log_file, "r") as log, open(rmsd_file, "w") as rmsd:
        rmsd.write(
            "# SSM superposition: alignment line (model vs reference) then RMSD block per pair\n"
        )
        rmsd.write("# Extracted from: {}\n\n".format(log_file))
        for line in log:
            if "Superposing" in line and " onto " in line:
                last_alignment_line = line.strip()
            elif "Aligning" in line and " to " in line:
                last_alignment_line = line.strip()
            elif "INFO: core rmsd" in line:
                if current_block:
                    rmsd.write("".join(current_block))
                    rmsd.write("\n")
                current_block = []
                if last_alignment_line:
                    current_block.append("{}\n".format(last_alignment_line))
                current_block.append(line)
            elif any(
                x in line
                for x in [
                    "number of residues",
                    "number of aligned",
                    "number of gaps",
                    "number of misdirections",
                    "number of SSE",
                    "sequence identity",
                ]
            ):
                current_block.append(line)
        if current_block:
            rmsd.write("".join(current_block))
            rmsd.write("\n")


def extract_rmsd_values(
    log_file,
    aligned_pattern=None,
    reference_pattern=None,
    case_sensitive=False,
    debug=False,
):
    """
    Extract RMSD and related values from a Coot log file (LSQ-style: Aligning … to …).

    Parameters:
    - log_file: Path to the Coot log file
    - aligned_pattern: Pattern to match in the aligned model name (optional)
    - reference_pattern: Pattern to match in the reference model name (optional)
    - case_sensitive: Whether to match case sensitively
    - debug: Whether to create a debug file with all decisions
    """
    if not os.path.exists(log_file):
        print(f"Error: Log file '{log_file}' not found.")
        return

    output_dir = os.path.dirname(log_file)
    base = os.path.basename(log_file)

    suffix = ""
    if aligned_pattern and reference_pattern:
        suffix = f"_aligned_{aligned_pattern}_ref_{reference_pattern}"
    elif aligned_pattern:
        suffix = f"_aligned_{aligned_pattern}"
    elif reference_pattern:
        suffix = f"_ref_{reference_pattern}"
    else:
        if base.startswith("coot_log_") and base.endswith(".txt"):
            suffix = "_" + base[9:-4]

    rmsd_file = os.path.join(output_dir, f"rmsd_values{suffix}.txt")
    debug_file = os.path.join(output_dir, f"rmsd_debug{suffix}.txt") if debug else None

    print(f"Extracting RMSD values (LSQ) from {log_file}")
    if aligned_pattern:
        print(f"Filtering for aligned models matching: {aligned_pattern}")
    if reference_pattern:
        print(f"Filtering for reference models matching: {reference_pattern}")
    print(f"Writing to {rmsd_file}")

    all_blocks = []
    current_block = []
    current_alignment = None

    with open(log_file, "r") as log:
        for line in log:
            if "Aligning" in line and "to" in line:
                if current_block and current_alignment:
                    all_blocks.append((current_alignment, current_block))

                current_alignment = line.strip()
                current_block = [line]
            elif "INFO: core rmsd" in line or any(
                x in line
                for x in [
                    "number of residues",
                    "number of aligned",
                    "number of gaps",
                    "number of misdirections",
                    "number of SSE",
                    "sequence identity",
                ]
            ):
                if current_block:
                    current_block.append(line)

    if current_block and current_alignment:
        all_blocks.append((current_alignment, current_block))

    filtered_blocks = []
    debug_info = []

    for alignment, block in all_blocks:
        match = re.match(r"Aligning\s+(.+?)\s+to\s+(.+?)$", alignment)
        if not match:
            continue

        aligned_model = match.group(1)
        reference_model = match.group(2)

        if not case_sensitive:
            aligned_check = aligned_model.lower()
            reference_check = reference_model.lower()
            aligned_pattern_check = aligned_pattern.lower() if aligned_pattern else None
            reference_pattern_check = reference_pattern.lower() if reference_pattern else None
        else:
            aligned_check = aligned_model
            reference_check = reference_model
            aligned_pattern_check = aligned_pattern
            reference_pattern_check = reference_pattern

        aligned_match = True if aligned_pattern is None else aligned_pattern_check in aligned_check
        reference_match = (
            True if reference_pattern is None else reference_pattern_check in reference_check
        )

        if (
            reference_pattern
            and reference_pattern_check in aligned_check
            and reference_pattern_check not in reference_check
        ):
            reference_match = False

        if (
            aligned_pattern
            and aligned_pattern_check in reference_check
            and aligned_pattern_check not in aligned_check
        ):
            aligned_match = False

        include = aligned_match and reference_match

        if debug:
            decision = {
                "alignment": alignment,
                "aligned_model": aligned_model,
                "reference_model": reference_model,
                "aligned_match": aligned_match,
                "reference_match": reference_match,
                "include": include,
            }
            debug_info.append(decision)

        if include:
            filtered_blocks.append(block)

    # Normalized patterns for debug output (do not reuse loop locals — some iterations `continue` early)
    dbg_aligned_pat = (
        aligned_pattern.lower()
        if aligned_pattern and not case_sensitive
        else aligned_pattern
    )
    dbg_reference_pat = (
        reference_pattern.lower()
        if reference_pattern and not case_sensitive
        else reference_pattern
    )

    with open(rmsd_file, "w") as out:
        for block in filtered_blocks:
            out.write("".join(block))
            out.write("\n")

    if debug and debug_file:
        with open(debug_file, "w") as debug_out:
            debug_out.write("# RMSD Extraction Debug Information\n")
            for info in debug_info:
                debug_out.write(f"## {info['alignment']}\n")
                debug_out.write(f"Aligned model: {info['aligned_model']}\n")
                debug_out.write(f"Reference model: {info['reference_model']}\n")

                if aligned_pattern:
                    aligned_in_aligned = dbg_aligned_pat in (
                        info["aligned_model"].lower()
                        if not case_sensitive
                        else info["aligned_model"]
                    )
                    aligned_in_reference = dbg_aligned_pat in (
                        info["reference_model"].lower()
                        if not case_sensitive
                        else info["reference_model"]
                    )
                    debug_out.write(
                        f"Aligned pattern '{aligned_pattern}' in aligned model: {aligned_in_aligned}\n"
                    )
                    debug_out.write(
                        f"Aligned pattern '{aligned_pattern}' in reference model: {aligned_in_reference}\n"
                    )

                if reference_pattern:
                    reference_in_aligned = dbg_reference_pat in (
                        info["aligned_model"].lower()
                        if not case_sensitive
                        else info["aligned_model"]
                    )
                    reference_in_reference = dbg_reference_pat in (
                        info["reference_model"].lower()
                        if not case_sensitive
                        else info["reference_model"]
                    )
                    debug_out.write(
                        f"Reference pattern '{reference_pattern}' in aligned model: {reference_in_aligned}\n"
                    )
                    debug_out.write(
                        f"Reference pattern '{reference_pattern}' in reference model: {reference_in_reference}\n"
                    )

                debug_out.write(f"Aligned match: {info['aligned_match']}\n")
                debug_out.write(f"Reference match: {info['reference_match']}\n")
                debug_out.write(f"DECISION: {'INCLUDE' if info['include'] else 'EXCLUDE'}\n\n")

    print(
        f"Extracted {len(filtered_blocks)} of {len(all_blocks)} RMSD blocks that match the filter criteria."
    )
    if debug:
        print(f"Debug information written to {debug_file}")


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Extract RMSD blocks from Coot logs. LSQ logs use 'Aligning … to …'; "
            "SSM logs use 'Superposing … onto …'. Use --format or rely on --format auto."
        )
    )
    parser.add_argument("log_file", help="Path to the Coot log file")
    parser.add_argument(
        "--format",
        "-f",
        choices=("auto", "lsq", "ssm"),
        default="auto",
        help="Log type: lsq (filtered LSQ blocks), ssm (SSM alignment + RMSD), or auto (detect)",
    )
    parser.add_argument("--aligned", "-a", help="LSQ only: pattern for aligned model names")
    parser.add_argument("--reference", "-r", help="LSQ only: pattern for reference model names")
    parser.add_argument(
        "--case-sensitive", "-c", action="store_true", help="LSQ only: case-sensitive pattern match"
    )
    parser.add_argument(
        "--debug", "-d", action="store_true", help="LSQ only: write rmsd_debug*.txt"
    )

    args = parser.parse_args()
    log_file = os.path.abspath(args.log_file)

    fmt = args.format
    if fmt == "auto":
        fmt = detect_rmsd_log_format(log_file)
        print(f"Detected log format: {fmt}")

    if fmt == "ssm":
        if args.aligned or args.reference or args.debug:
            print(
                "Warning: --aligned, --reference, and --debug apply to LSQ mode only; ignoring.",
                file=sys.stderr,
            )
        extract_ssm_rmsd_values(log_file)
    else:
        extract_rmsd_values(
            log_file,
            args.aligned,
            args.reference,
            args.case_sensitive,
            args.debug,
        )


if __name__ == "__main__":
    main()
