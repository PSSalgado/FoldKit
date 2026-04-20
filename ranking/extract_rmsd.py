import os
import sys
import argparse
import re

_REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if _REPO_ROOT not in sys.path:
    sys.path.insert(0, _REPO_ROOT)

from cli_log import add_log_args, setup_log_from_args

_EXTRACT_EPILOG = """
Examples (from repository root):
  Single Coot log:
    python ranking/extract_rmsd.py --format ssm path/to/coot_log.txt
    python ranking/extract_rmsd.py --format lsq path/to/coot_log.txt --aligned=set_a --reference=ref_01
  Optional `file` keyword (same as a single positional log):
    python ranking/extract_rmsd.py file path/to/coot_log.txt -o /tmp/rmsd.txt
  Scan every coot_log.txt / coot_log_*.txt under a directory (-o = mirror root for outputs):
    python ranking/extract_rmsd.py --dir=/path/to/base --format ssm
    python ranking/extract_rmsd.py --dir /path/to/base -o /path/to/output_root_dir
"""


def _resolved_path(path: str) -> str:
    """Expand ~ and resolve to absolute path (Python does not expand ~ in abspath)."""
    return os.path.abspath(os.path.expanduser(path))


def _peek_log_head(path: str, max_bytes: int = 65536) -> str:
    try:
        with open(path, "rb") as f:
            return f.read(max_bytes).decode("utf-8", errors="replace")
    except OSError:
        return ""


def is_standard_coot_log_basename(name: str) -> bool:
    """
    True for Coot log names from FoldKit scripts: coot_log.txt or coot_log_<suffix>.txt.
    Excludes unrelated names such as coot_logfiles.txt or coot_logged.txt.
    """
    if name == "coot_log.txt":
        return True
    return name.startswith("coot_log_") and name.endswith(".txt")


def detect_rmsd_log_format(log_file: str) -> str:
    """
    Return 'ssm' if the log looks like Coot SSM (FoldKit format), else 'lsq'.

    FoldKit writes '# SSM alignment log' (etc.) before Coot runs; prefer that so
    --format auto works even if Coot output order differs. Also detect the
    standard script lines 'Superposing … onto …'.
    """
    sample = _peek_log_head(log_file)
    if "# SSM" in sample:
        return "ssm"
    if "# LSQ" in sample:
        return "lsq"
    if "Superposing" in sample and " onto " in sample:
        return "ssm"
    return "lsq"


def _output_kwargs_from_o_flag(output_arg):
    """
    Map --output to extractor kwargs.

    - Existing directory -> output_dir
    - Trailing path separator (e.g. out/ or out\\) -> output_dir (created if missing)
    - Otherwise -> rmsd_output_path (file path, extension optional; parent dirs created in extractors)

    To target a new directory whose name has no extension, pass a trailing separator (e.g. out/results/).
    """
    if not output_arg:
        return {}
    expanded = os.path.expanduser(output_arg)
    abs_norm = os.path.abspath(expanded)
    if os.path.isdir(abs_norm):
        return {"output_dir": abs_norm}
    stripped = expanded.rstrip("/\\")
    had_trailing_sep = len(stripped) < len(expanded)
    if had_trailing_sep:
        d = os.path.abspath(stripped) if stripped else abs_norm
        os.makedirs(d, exist_ok=True)
        return {"output_dir": d}
    return {"rmsd_output_path": abs_norm}


def extract_ssm_rmsd_values(log_file: str, output_dir=None, rmsd_output_path=None) -> None:
    """
    Extract molecule/superposition info and RMSD values from a Coot log for SSM.
    SSM logs use different wording than LSQ: scripts print "Superposing X onto …",
    and Coot may print "Aligning … to …". We capture whichever alignment line appears
    before each "INFO: core rmsd" block and write it together with the RMSD block.

    output_dir: if set, RMSD text uses the standard basename under this directory.
    rmsd_output_path: if set, write RMSD to this exact file (parent dirs created).
    When both are set, rmsd_output_path wins.
    """
    log_file = _resolved_path(log_file)
    if not os.path.exists(log_file):
        print(f"Error: Log file '{log_file}' not found.")
        return

    if os.path.isdir(log_file):
        print(
            f"Error: Expected a Coot log file, but '{log_file}' is a directory. "
            "Use --dir DIR to scan for coot_log.txt / coot_log_*.txt."
        )
        return
    if rmsd_output_path is not None:
        rmsd_file = _resolved_path(rmsd_output_path)
        os.makedirs(os.path.dirname(rmsd_file), exist_ok=True)
    else:
        if output_dir is not None:
            output_dir = _resolved_path(output_dir)
        else:
            output_dir = os.path.dirname(log_file)
        base = os.path.basename(log_file)
        if base == "coot_log.txt":
            rmsd_basename = "rmsd_SSM_values.txt"
        elif base.startswith("coot_log_") and base.endswith(".txt"):
            suffix = base[9:-4]
            rmsd_basename = "rmsd_SSM_values_{}.txt".format(suffix)
        else:
            rmsd_basename = "rmsd_SSM_values.txt"
        os.makedirs(output_dir, exist_ok=True)
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
    output_dir=None,
    rmsd_output_path=None,
):
    """
    Extract RMSD and related values from a Coot log file (LSQ-style: Aligning … to …).

    Coot LSQ logs report deviations after each pair, e.g. INFO:: LSQ matched … atoms,
    mean/rms/max/min devi — not SSM's "INFO: core rmsd". Those lines are collected here.

    Parameters:
    - log_file: Path to the Coot log file
    - aligned_pattern: Pattern to match in the aligned model name (optional)
    - reference_pattern: Pattern to match in the reference model name (optional)
    - case_sensitive: Whether to match case sensitively
    - debug: Whether to create a debug file with all decisions
    - output_dir: if set, RMSD (and debug) use standard basenames under this directory
    - rmsd_output_path: if set, write RMSD to this exact file (parent dirs created)
    When rmsd_output_path is set, debug output (if any) is placed next to that file.
    """
    log_file = _resolved_path(log_file)
    if not os.path.exists(log_file):
        print(f"Error: Log file '{log_file}' not found.")
        return

    if os.path.isdir(log_file):
        print(
            f"Error: Expected a Coot log file, but '{log_file}' is a directory. "
            "Use --dir DIR to scan for coot_log.txt / coot_log_*.txt."
        )
        return
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

    if rmsd_output_path is not None:
        rmsd_file = _resolved_path(rmsd_output_path)
        os.makedirs(os.path.dirname(rmsd_file), exist_ok=True)
        output_dir = os.path.dirname(rmsd_file)
    else:
        if output_dir is not None:
            output_dir = _resolved_path(output_dir)
        else:
            output_dir = os.path.dirname(log_file)
        os.makedirs(output_dir, exist_ok=True)
        rmsd_file = os.path.join(output_dir, f"rmsd_values{suffix}.txt")

    debug_file = os.path.join(output_dir, f"rmsd_debug{suffix}.txt") if debug else None

    print(f"Extracting RMSD values (LSQ) from {log_file}")
    if aligned_pattern:
        print(f"Filtering for aligned models matching: {aligned_pattern}")
    if reference_pattern:
        print(f"Filtering for reference models matching: {reference_pattern}")
    print(f"Writing to {rmsd_file}")

    # Lines that belong to the RMSD / deviation block after an "Aligning … to …" line.
    # SSM-style logs (sometimes mixed in): "INFO: core rmsd", residue counts.
    # Coot LSQ: "INFO:: LSQ matched", "rms devi:", etc. (often double-colon INFO::).
    _lsq_rmsd_markers = (
        "INFO: core rmsd",
        "number of residues",
        "number of aligned",
        "number of gaps",
        "number of misdirections",
        "number of SSE",
        "sequence identity",
        "LSQ matched",
        "matched atoms had",
        "Matching/moving molecule",
        "mean devi",
        "rms devi",
        "max devi",
        "min devi",
    )

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
            elif current_block and any(x in line for x in _lsq_rmsd_markers):
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

    # Normalised patterns for debug output (do not reuse loop locals — some iterations `continue` early)
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


def _run_dir_scan(base_dir: str, args) -> None:
    """Walk base_dir for standard Coot logs and extract RMSD for each."""
    print(f"Searching for Coot logs under {base_dir}")
    log_files = []
    for root, _dirs, files in os.walk(base_dir):
        for name in files:
            if is_standard_coot_log_basename(name):
                log_files.append(os.path.join(root, name))

    print(f"Found {len(log_files)} log files")
    output_root = _resolved_path(args.output) if args.output else None
    for log_file in sorted(log_files):
        fmt = args.format
        if fmt == "auto":
            fmt = detect_rmsd_log_format(log_file)
            print(f"\n[scan] {log_file}: detected format {fmt}")
        else:
            print(f"\n[scan] {log_file}: using format {fmt}")

        if output_root:
            rel = os.path.relpath(os.path.dirname(log_file), base_dir)
            out_dir = os.path.join(output_root, rel)
            os.makedirs(out_dir, exist_ok=True)
            if fmt == "ssm":
                if args.aligned or args.reference:
                    print(
                        "Warning: --aligned/--reference apply to LSQ mode only; ignoring for SSM.",
                        file=sys.stderr,
                    )
                extract_ssm_rmsd_values(log_file, output_dir=out_dir)
            else:
                extract_rmsd_values(
                    log_file,
                    args.aligned,
                    args.reference,
                    args.case_sensitive,
                    args.debug,
                    output_dir=out_dir,
                )
        else:
            if fmt == "ssm":
                if args.aligned or args.reference:
                    print(
                        "Warning: --aligned/--reference apply to LSQ mode only; ignoring for SSM.",
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


def _run_cli(argv):
    """Single parser: one log file, or --dir for recursive scan. Optional leading `file` is stripped in main."""
    parser = argparse.ArgumentParser(
        description=(
            "Extract RMSD blocks from Coot logs. LSQ logs use 'Aligning … to …'; "
            "SSM logs use 'Superposing … onto …'. Use --format or rely on --format auto. "
            "With --dir, walk the tree for coot_log.txt / coot_log_*.txt."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=_EXTRACT_EPILOG,
    )
    parser.add_argument(
        "log_file",
        nargs="?",
        default=None,
        help="Path to one Coot log (omit when using --dir).",
    )
    parser.add_argument(
        "--dir",
        dest="scan_dir",
        metavar="DIR",
        default=None,
        help=(
            "Base directory: recursively find coot_log.txt and coot_log_*.txt (not other coot_log*.txt names)."
        ),
    )
    parser.add_argument(
        "--format",
        "-f",
        choices=("auto", "lsq", "ssm"),
        default="auto",
        help="Log type: lsq, ssm, or auto (detect). With --dir, applied per file (or auto per file).",
    )
    parser.add_argument("--aligned", "-a", help="LSQ only: pattern for aligned model names.")
    parser.add_argument("--reference", "-r", help="LSQ only: pattern for reference model names.")
    parser.add_argument(
        "--case-sensitive", "-c", action="store_true", help="LSQ only: case-sensitive pattern match."
    )
    parser.add_argument(
        "--debug", "-d", action="store_true", help="LSQ only: write rmsd_debug*.txt."
    )
    parser.add_argument(
        "-o",
        "--output",
        dest="output",
        default=None,
        help=(
            "Single log: output RMSD file or directory (same rules as before). "
            "With --dir: optional root directory; RMSD files are written under it mirroring the input tree."
        ),
    )
    add_log_args(parser)

    args = parser.parse_args(argv)

    if args.scan_dir is not None and not str(args.scan_dir).strip():
        parser.error(
            "--dir was given without a path. "
            "Do not put a space after '=' (e.g. use --dir=/path or --dir /path, not --dir= /path)."
        )

    if args.scan_dir and args.log_file:
        parser.error("Use either --dir or a log file path, not both.")

    if args.scan_dir:
        scan_dir = _resolved_path(args.scan_dir)
        setup_log_from_args(
            args,
            script_path=__file__,
            inputs=[scan_dir],
            pattern=None,
        )
        _run_dir_scan(scan_dir, args)
        return

    if not args.log_file:
        parser.error("Provide a Coot log file path, or use --dir=DIR.")

    log_abs = _resolved_path(args.log_file)
    if os.path.isdir(log_abs):
        parser.error(
            "The positional path {!r} is a directory, not a Coot log file. "
            "Use directory scan: --dir {} [--format ssm] [-o ...]. "
            "If you used --dir= /path, remove the space after '='.".format(
                args.log_file,
                args.log_file,
            )
        )

    setup_log_from_args(
        args,
        script_path=__file__,
        inputs=[log_abs],
        pattern=None,
    )
    log_file = log_abs
    out_kw = _output_kwargs_from_o_flag(args.output)

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
        extract_ssm_rmsd_values(log_file, **out_kw)
    else:
        extract_rmsd_values(
            log_file,
            args.aligned,
            args.reference,
            args.case_sensitive,
            args.debug,
            **out_kw,
        )


def _use_file_subcommand(argv) -> bool:
    """`file LOG ...` — optional keyword so a log literally named 'file' is not ambiguous."""
    if not argv:
        return False
    if argv[0] != "file":
        return False
    if len(argv) >= 2 and argv[1] in ("-h", "--help"):
        return True
    if len(argv) == 1:
        return not os.path.isfile("file")
    if argv[1].startswith("-"):
        return False
    return True


def main():
    argv = sys.argv[1:]
    if _use_file_subcommand(argv):
        argv = argv[1:]
        if not argv:
            print("error: 'file' requires a Coot log path (see --help).", file=sys.stderr)
            sys.exit(2)
    _run_cli(argv)


if __name__ == "__main__":
    main()
