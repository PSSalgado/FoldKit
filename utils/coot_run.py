"""
Run Coot with a generated Python script and capture output to a log file.
"""

from __future__ import annotations

import os
import shutil
import sys
import tempfile
from subprocess import PIPE, Popen, STDOUT

COOT_STDOUT_ALIGNMENTS_DONE_MARKER = "FOLDKIT_ALIGNMENTS_DONE"


def coot_available() -> bool:
    """Return True if ``coot`` is on PATH."""
    return shutil.which("coot") is not None


def pipe_coot_stdout_line(line: str, log_f, *, echo_all_stdout: bool = False) -> None:
    """Append Coot stdout to the log; echo marker lines to stderr."""
    log_f.write(line)
    log_f.flush()
    stripped = line.rstrip("\n\r")
    if COOT_STDOUT_ALIGNMENTS_DONE_MARKER in line:
        print(stripped, file=sys.stderr, flush=True)
    if echo_all_stdout:
        print(stripped)


def run_coot_script(
    script_content: str,
    log_path: str,
    *,
    header_lines: list[str] | None = None,
    not_interactive: bool = True,
    echo_all_stdout: bool = False,
) -> int:
    """
    Write ``script_content`` to a temp script, run ``coot --script``, log stdout.

    Returns Coot process exit code (0 on success).
    """
    if not coot_available():
        print("ERROR: coot not found on PATH", file=sys.stderr)
        return 127

    log_path = os.path.abspath(os.path.expanduser(log_path))
    os.makedirs(os.path.dirname(log_path) or ".", exist_ok=True)

    fd, script_path = tempfile.mkstemp(suffix="_coot_script.py", prefix="foldkit_")
    try:
        with os.fdopen(fd, "w", encoding="utf-8") as fh:
            fh.write(script_content)

        with open(log_path, "w", encoding="utf-8") as log:
            if header_lines:
                for line in header_lines:
                    log.write(line.rstrip("\n") + "\n")
                log.write("\n")

            process = Popen(
                ["coot", "--script", script_path],
                stdout=PIPE,
                stderr=STDOUT,
                universal_newlines=True,
                bufsize=1,
            )
            if process.stdout is not None:
                while True:
                    output = process.stdout.readline()
                    if output == "" and process.poll() is not None:
                        break
                    if output:
                        pipe_coot_stdout_line(
                            output, log, echo_all_stdout=echo_all_stdout
                        )
            else:
                print(
                    "Warning: Coot stdout unavailable; check {}".format(log_path),
                    file=sys.stderr,
                )
            return process.wait()
    finally:
        try:
            os.unlink(script_path)
        except OSError:
            pass
