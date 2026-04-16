#!/usr/bin/env python3
"""
CLI log helper.

Provides:
- standard log-related argparse options
- tee of stdout/stderr to a per-run log file

Default log filename is derived from:
  <script-stem>__<input-slug>[__<tag>].log
where input-slug is derived from input paths / directory / pattern text.
"""

from __future__ import annotations

import os
import re
import sys
import atexit
from dataclasses import dataclass
from datetime import datetime
from typing import TextIO


def add_log_args(parser) -> None:
    parser.add_argument(
        "--log",
        metavar="PATH",
        default=None,
        help=(
            "Write a log file (stdout/stderr tee). If PATH is a directory, a default log filename is used. "
            "If PATH contains '{}', it is replaced by the input-derived slug."
        ),
    )
    parser.add_argument(
        "--log-tag",
        metavar="TEXT",
        default="",
        help="Optional tag appended to the default log filename.",
    )
    parser.add_argument(
        "--no-log",
        action="store_true",
        help="Disable log file creation.",
    )


def strip_log_args_from_argv(argv: list[str]) -> tuple[list[str], dict[str, str | bool]]:
    """
    Remove standard logging flags from argv.

    Supports:
    - --log PATH / --log=PATH
    - --log-tag TEXT / --log-tag=TEXT
    - --no-log

    Returns (argv_without_log_flags, parsed_log_kwargs).
    """
    out: list[str] = []
    parsed: dict[str, str | bool] = {"log": None, "log_tag": "", "no_log": False}

    it = iter(range(len(argv)))
    i = 0
    while i < len(argv):
        a = argv[i]
        if a == "--no-log":
            parsed["no_log"] = True
            i += 1
            continue
        if a.startswith("--log="):
            parsed["log"] = a.split("=", 1)[1]
            i += 1
            continue
        if a.startswith("--log-tag="):
            parsed["log_tag"] = a.split("=", 1)[1]
            i += 1
            continue
        if a == "--log":
            if i + 1 < len(argv):
                parsed["log"] = argv[i + 1]
                i += 2
                continue
            # fall through: treat as normal arg
        if a == "--log-tag":
            if i + 1 < len(argv):
                parsed["log_tag"] = argv[i + 1]
                i += 2
                continue
        out.append(a)
        i += 1

    return out, parsed


class _ArgsFromDict:
    def __init__(self, d: dict[str, str | bool]):
        for k, v in d.items():
            setattr(self, k, v)


def _slugify(text: str) -> str:
    text = text.strip().replace(os.sep, "_")
    text = re.sub(r"[^A-Za-z0-9._-]+", "_", text)
    text = re.sub(r"_+", "_", text).strip("_")
    return text or "run"


def _input_slug(inputs: list[str] | None, pattern: str | None = None) -> str:
    """
    Derive a stable slug from typical CLI inputs.
    """
    bits: list[str] = []
    if inputs:
        # Prefer a single input; otherwise summarise.
        if len(inputs) == 1:
            p = inputs[0]
            p = os.path.abspath(os.path.expanduser(p))
            base = os.path.basename(p.rstrip(os.sep)) or "input"
            stem, _ext = os.path.splitext(base)
            bits.append(stem or base)
        else:
            bits.append(f"{len(inputs)}inputs")
    if pattern:
        bits.append(pattern)
    return _slugify("__".join(bits) if bits else "run")


def default_log_name(script_path: str, inputs: list[str] | None, pattern: str | None, tag: str) -> str:
    script_stem = os.path.splitext(os.path.basename(script_path))[0]
    slug = _input_slug(inputs, pattern=pattern)
    tag_part = _slugify(tag) if tag else ""
    base = f"{script_stem}__{slug}"
    if tag_part:
        base = f"{base}__{tag_part}"
    return f"{base}.log"


class _Tee(TextIO):
    def __init__(self, a: TextIO, b: TextIO):
        self._a = a
        self._b = b

    def write(self, s: str) -> int:  # type: ignore[override]
        na = self._a.write(s)
        self._b.write(s)
        return na

    def flush(self) -> None:  # type: ignore[override]
        self._a.flush()
        self._b.flush()


@dataclass
class LogSetup:
    path: str
    file: TextIO
    stdout_prev: TextIO
    stderr_prev: TextIO

    def close(self) -> None:
        sys.stdout = self.stdout_prev
        sys.stderr = self.stderr_prev
        try:
            self.file.flush()
        finally:
            self.file.close()


def setup_log_from_args(
    args,
    *,
    script_path: str,
    inputs: list[str] | None = None,
    pattern: str | None = None,
) -> LogSetup | None:
    """
    If logging is enabled, tee stdout/stderr to a log file and return a LogSetup.
    """
    if getattr(args, "no_log", False):
        return None

    log_arg = getattr(args, "log", None)
    tag = getattr(args, "log_tag", "") or ""
    name = default_log_name(script_path, inputs, pattern, tag)

    if not log_arg:
        dest = os.path.join(os.getcwd(), name)
    else:
        log_arg = os.path.abspath(os.path.expanduser(str(log_arg)))
        if "{}" in log_arg:
            slug = _input_slug(inputs, pattern=pattern)
            dest = log_arg.replace("{}", slug)
        elif os.path.isdir(log_arg) or log_arg.endswith(os.sep):
            os.makedirs(log_arg, exist_ok=True)
            dest = os.path.join(log_arg, name)
        else:
            parent = os.path.dirname(log_arg) or "."
            os.makedirs(parent, exist_ok=True)
            dest = log_arg

    # Ensure .log suffix for defaults/templates unless user provided an explicit filename with a suffix.
    if not os.path.splitext(dest)[1]:
        dest = dest + ".log"

    f = open(dest, "w", encoding="utf-8", newline="\n")
    # Header line with timestamp for quick identification.
    ts = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    f.write(f"# Log started: {ts}\n")
    f.write(f"# Command: {' '.join(map(str, sys.argv))}\n\n")
    f.flush()

    stdout_prev = sys.stdout
    stderr_prev = sys.stderr
    sys.stdout = _Tee(stdout_prev, f)  # type: ignore[assignment]
    sys.stderr = _Tee(stderr_prev, f)  # type: ignore[assignment]
    ls = LogSetup(path=dest, file=f, stdout_prev=stdout_prev, stderr_prev=stderr_prev)
    atexit.register(ls.close)
    return ls


def setup_log_from_argv(
    *,
    script_path: str,
    argv: list[str],
    inputs: list[str] | None = None,
    pattern: str | None = None,
) -> tuple[list[str], LogSetup | None]:
    """
    For scripts that do not use argparse: strip logging args from argv and set up tee logging.

    Returns (argv_without_log_flags, LogSetup|None). The caller may want to assign
    the returned argv back to sys.argv.
    """
    argv2, parsed = strip_log_args_from_argv(argv)
    args = _ArgsFromDict(parsed)
    ls = setup_log_from_args(args, script_path=script_path, inputs=inputs, pattern=pattern)
    return argv2, ls

