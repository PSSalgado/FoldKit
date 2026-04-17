#!/usr/bin/env python3
"""
CLI log helper (summary-only).

FoldKit scripts can optionally write a short run summary to a log file via `--log`.
There is **no automatic logging** and **no stdout/stderr tee**.

The summary log is intended to capture:
- basic run context (timestamp, command, cwd, inputs)
- a list of high-level tasks/steps (when scripts report them)
- any errors / uncaught exceptions (traceback captured via sys.excepthook)
"""

from __future__ import annotations

import os
import re
import sys
import atexit
import traceback
from dataclasses import dataclass
from datetime import datetime
from typing import TextIO


def add_log_args(parser) -> None:
    parser.add_argument(
        "--log",
        nargs="?",
        const="__AUTO__",
        metavar="FILE",
        default=None,
        help=(
            "Optional summary log file. If provided without a value, defaults to '<script_stem>.log'. "
            "This log is a short summary (tasks + errors); it does not mirror stdout/stderr."
        ),
    )


def strip_log_args_from_argv(argv: list[str]) -> tuple[list[str], dict[str, str | None]]:
    """
    Remove standard logging flags from argv.

    Supports:
    - --log FILE / --log=FILE
    - --log (no value) -> "__AUTO__"

    Returns (argv_without_log_flags, parsed_log_kwargs).
    """
    out: list[str] = []
    parsed: dict[str, str | None] = {"log": None}

    it = iter(range(len(argv)))
    i = 0
    while i < len(argv):
        a = argv[i]
        if a.startswith("--log="):
            parsed["log"] = a.split("=", 1)[1]
            i += 1
            continue
        if a == "--log":
            if i + 1 < len(argv):
                nxt = argv[i + 1]
                if nxt.startswith("-"):
                    parsed["log"] = "__AUTO__"
                    i += 1
                    continue
                parsed["log"] = nxt
                i += 2
                continue
            parsed["log"] = "__AUTO__"
            i += 1
            continue
        out.append(a)
        i += 1

    return out, parsed


class _ArgsFromDict:
    def __init__(self, d: dict[str, str | None]):
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


def default_log_name(script_path: str) -> str:
    script_stem = os.path.splitext(os.path.basename(script_path))[0]
    return f"{script_stem}.log"


@dataclass
class SummaryLog:
    path: str
    file: TextIO
    tasks: list[str]
    errors: list[str]
    _had_uncaught: bool = False
    _prev_excepthook: object | None = None

    def task(self, msg: str) -> None:
        msg = str(msg).strip()
        if not msg:
            return
        self.tasks.append(msg)

    def kv(self, key: str, value) -> None:
        """
        Record a key/value line in the task list.

        This keeps logs consistent across scripts while staying plain text.
        """
        k = str(key).strip()
        if not k:
            return
        v = value
        try:
            if isinstance(v, (list, tuple)) and len(v) > 10:
                v = list(v[:10]) + ["..."]
        except Exception:
            pass
        self.task(f"{k}: {v}")

    def progress(self, i: int, n: int, label: str) -> None:
        """Record a simple progress marker like '[3/12] filename'."""
        try:
            ii = int(i)
            nn = int(n)
        except Exception:
            self.task(str(label))
            return
        self.task(f"[{ii}/{nn}] {label}")

    def error(self, msg: str) -> None:
        msg = str(msg).strip()
        if not msg:
            return
        self.errors.append(msg)

    def close(self) -> None:
        try:
            self._write_footer()
            self.file.flush()
        finally:
            self.file.close()
        if self._prev_excepthook is not None:
            try:
                sys.excepthook = self._prev_excepthook  # type: ignore[assignment]
            except Exception:
                pass

    def _write_footer(self) -> None:
        self.file.write("\n")
        self.file.write("# Summary\n")
        if self.tasks:
            self.file.write("# Tasks\n")
            for t in self.tasks:
                self.file.write(f"- {t}\n")
        if self.errors:
            self.file.write("# Errors\n")
            for e in self.errors:
                self.file.write(f"- {e}\n")
        if not self.errors and not self._had_uncaught:
            self.file.write("# Status: OK\n")
        else:
            self.file.write("# Status: ERROR\n")


def setup_log_from_args(
    args,
    *,
    script_path: str,
    inputs: list[str] | None = None,
    pattern: str | None = None,
) -> SummaryLog | None:
    """
    If --log is provided, write a summary log and return a SummaryLog.
    """
    log_arg = getattr(args, "log", None)
    if not log_arg:
        return None

    if str(log_arg) == "__AUTO__":
        dest = default_log_name(script_path)
    else:
        dest = str(log_arg)

    dest = os.path.abspath(os.path.expanduser(dest))
    parent = os.path.dirname(dest) or "."
    os.makedirs(parent, exist_ok=True)
    if not os.path.splitext(dest)[1]:
        dest = dest + ".log"

    f = open(dest, "w", encoding="utf-8", newline="\n")
    ts = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    f.write(f"# Started: {ts}\n")
    f.write(f"# Command: {' '.join(map(str, sys.argv))}\n")
    f.write(f"# CWD: {os.getcwd()}\n")
    f.write(f"# Script: {os.path.abspath(script_path)}\n")
    f.write(f"# Python: {sys.executable}\n")
    if inputs:
        f.write(f"# Inputs: {inputs}\n")
    if pattern:
        f.write(f"# Pattern: {pattern}\n")
    f.write("\n")

    log = SummaryLog(path=dest, file=f, tasks=[], errors=[])

    # Capture uncaught exceptions into the summary log.
    prev_hook = getattr(sys, "excepthook", None)
    log._prev_excepthook = prev_hook  # type: ignore[assignment]

    def _hook(exc_type, exc, tb):
        log._had_uncaught = True
        try:
            log.error(f"Uncaught exception: {exc_type.__name__}: {exc}")
            f.write("\n# Uncaught exception\n")
            f.write("".join(traceback.format_exception(exc_type, exc, tb)))
            f.write("\n")
            f.flush()
        except Exception:
            pass
        if callable(prev_hook):
            prev_hook(exc_type, exc, tb)

    try:
        sys.excepthook = _hook  # type: ignore[assignment]
    except Exception:
        pass

    atexit.register(log.close)
    return log


def setup_log_from_argv(
    *,
    script_path: str,
    argv: list[str],
    inputs: list[str] | None = None,
    pattern: str | None = None,
) -> tuple[list[str], SummaryLog | None]:
    """
    For scripts that do not use argparse: strip logging args from argv and set up summary logging.

    Returns (argv_without_log_flags, SummaryLog|None). The caller may want to assign
    the returned argv back to sys.argv.
    """
    argv2, parsed = strip_log_args_from_argv(argv)
    args = _ArgsFromDict(parsed)
    ls = setup_log_from_args(args, script_path=script_path, inputs=inputs, pattern=pattern)
    return argv2, ls

