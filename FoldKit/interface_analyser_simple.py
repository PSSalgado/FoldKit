#!/usr/bin/env python3
"""
Interface analyser (simple / ASU-style)
=====================================

Pairwise chain–chain interface analysis for a single structure file.

This is the "simple" entrypoint: it reports pairwise interfaces and summary
statistics, but does not compute or print multi-copy lattice reference metrics.

For multi-copy / lattice-wide metrics, use `interface_analyser_lattice.py`.
"""

import os
import sys
import argparse

from interface_analyser import InterfaceAnalyser, collect_structure_paths, filter_paths_by_patterns, _run_analysis


def main():
    parser = argparse.ArgumentParser(
        description="Analyse pairwise interfaces between chains in structure files (simple mode).",
    )
    parser.add_argument(
        'input',
        nargs='+',
        help='PDB/CIF file(s), directory, or glob pattern (e.g. *.pdb).',
    )
    parser.add_argument(
        '--set', '-s',
        action='append',
        dest='sets',
        metavar='PATTERNS',
        help='Comma-separated patterns; file is included only if basename contains ALL patterns. Repeat for multiple sets.',
    )
    parser.add_argument(
        '--sets',
        dest='sets_multi',
        nargs='+',
        metavar='SET',
        help='Multiple set names in one go (one pattern per set).',
    )
    parser.add_argument(
        '--output', '-o',
        metavar='FILE',
        help='Output file. Single file: -o results.txt (all sets concatenated). Per-set: use "{}" for set label.',
    )
    parser.add_argument(
        '--per-structure', '-p',
        action='store_true',
        help='Write one output file per structure; -o must contain "{}" (replaced by file stem).',
    )
    parser.add_argument(
        '--dry-run',
        action='store_true',
        help='Print which files would be processed for each set and exit without running.',
    )
    parser.add_argument(
        '--chains',
        metavar='IDS',
        help="Focus analysis on specific chain IDs (comma-separated). Only chain pairs where at least one chain is in this list are analysed.",
    )
    args = parser.parse_args()

    focus_chains = None
    if getattr(args, 'chains', None):
        focus_chains = [c.strip() for c in str(args.chains).split(',') if c.strip()]

    paths = collect_structure_paths(args.input)

    # Per-structure mode (simple): identical to interface_analyser.py, but never passes reference_chain_id
    if args.per_structure:
        if not args.output or '{}' not in args.output:
            print("Error: --per-structure requires -o with '{}' in the path (e.g. -o '{}_interface.txt').", file=sys.stderr)
            sys.exit(1)
        paths_to_process = list(paths)
        if getattr(args, 'sets_multi', None) or args.sets:
            filtered_set = set()
            if getattr(args, 'sets_multi', None):
                for name in args.sets_multi:
                    patterns = [name.strip()] if name.strip() else []
                    if patterns:
                        filtered_set.update(filter_paths_by_patterns(paths, patterns))
            if args.sets:
                for s in args.sets:
                    patterns = [p.strip() for p in s.split(',') if p.strip()]
                    if patterns:
                        filtered_set.update(filter_paths_by_patterns(paths, patterns))
            paths_to_process = sorted(filtered_set)
        if not paths_to_process:
            print("No structure files match the filter.", file=sys.stderr)
            sys.exit(1)
        analyser = InterfaceAnalyser()
        for p in paths_to_process:
            stem = os.path.splitext(os.path.basename(p))[0]
            out_path = args.output.replace('{}', stem)
            print(f"Writing {os.path.basename(p)} -> {out_path}", file=sys.stderr)
            with open(out_path, 'w') as out:
                _run_analysis(analyser, [p], out, focus_chains=focus_chains, reference_chain_id=None)
        return

    # Build list of (label, filtered_paths). If --set/--sets not used, one set with all paths.
    set_list = []
    if getattr(args, 'sets_multi', None):
        for name in args.sets_multi:
            patterns = [name.strip()] if name.strip() else []
            if patterns:
                label = patterns[0]
                filtered = filter_paths_by_patterns(paths, patterns)
                set_list.append((label, patterns, filtered))
    if args.sets:
        for s in args.sets:
            patterns = [p.strip() for p in s.split(',') if p.strip()]
            if patterns:
                label = '_'.join(patterns)
                filtered = filter_paths_by_patterns(paths, patterns)
                set_list.append((label, patterns, filtered))
    if not set_list:
        set_list = [("all", [], list(paths))]

    if args.dry_run:
        lines = ["Dry run: no analysis will be performed.", ""]
        single_file = args.output and '{}' not in args.output
        for label, patterns, filtered in set_list:
            dest = args.output if single_file else (args.output.replace('{}', label) if args.output and '{}' in args.output else (args.output or "stdout"))
            lines.append(f"Set {label!r} (patterns: {patterns}): {len(filtered)} file(s) -> {dest}")
            for p in filtered:
                lines.append(f"  {os.path.basename(p)}")
            lines.append("")
        print("\n".join(lines), file=sys.stderr, flush=True)
        return

    if not paths:
        print("No structure files found.", file=sys.stderr)
        sys.exit(1)

    single_output_file = args.output and '{}' not in args.output
    out = None
    if single_output_file:
        out = open(args.output, 'w')
        print(f"Writing all sets to {args.output}", file=sys.stderr)

    analyser = InterfaceAnalyser()
    for label, patterns, filtered in set_list:
        if not filtered:
            print(f"No files match set {label!r} (patterns: {patterns}), skipping.", file=sys.stderr)
            continue
        if args.output and not single_output_file:
            out_path = args.output.replace('{}', label)
            out = open(out_path, 'w')
            print(f"Writing analysis for set {label!r} ({len(filtered)} file(s)) to {out_path}", file=sys.stderr)
        elif not single_output_file:
            out = sys.stdout
        try:
            _run_analysis(analyser, filtered, out, focus_chains=focus_chains, reference_chain_id=None)
        finally:
            if args.output and not single_output_file and out is not sys.stdout:
                out.close()

    if single_output_file and out is not None and out is not sys.stdout:
        out.close()


if __name__ == "__main__":
    main()

