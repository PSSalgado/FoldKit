"""Tee stdout/stderr to a run log for Dali-related FoldKit CLIs."""

from __future__ import annotations

import os
import sys
from datetime import datetime
from typing import Any, IO


class _TeeTextIO:
    __slots__ = ("primary", "log_fp")

    def __init__(self, primary: Any, log_fp: IO[str]) -> None:
        self.primary = primary
        self.log_fp = log_fp

    def write(self, s: str) -> int:
        self.primary.write(s)
        self.log_fp.write(s)
        return len(s)

    def flush(self) -> None:
        self.primary.flush()
        self.log_fp.flush()

    def isatty(self) -> bool:
        fn = getattr(self.primary, "isatty", None)
        return bool(fn()) if callable(fn) else False

    def fileno(self) -> int:
        return self.primary.fileno()

    @property
    def encoding(self) -> str:
        return getattr(self.primary, "encoding", "utf-8") or "utf-8"


def install_dali_run_log(log_path: str, banner: str) -> tuple[Any, Any, IO[str]]:
    path = os.path.abspath(os.path.expanduser(log_path.strip()))
    d = os.path.dirname(path)
    if d:
        os.makedirs(d, exist_ok=True)
    log_fp = open(path, "w", encoding="utf-8")
    log_fp.write(f"# {banner} — {datetime.now().isoformat(timespec='seconds')}\n")
    log_fp.flush()
    old_out, old_err = sys.stdout, sys.stderr
    sys.stdout = _TeeTextIO(old_out, log_fp)
    sys.stderr = _TeeTextIO(old_err, log_fp)
    return old_out, old_err, log_fp


def uninstall_dali_run_log(
    state: tuple[Any, Any, IO[str]] | None, banner: str, ok: bool = True
) -> None:
    if state is None:
        return
    old_out, old_err, log_fp = state
    sys.stdout, sys.stderr = old_out, old_err
    status = "finished" if ok else "exited with error"
    try:
        log_fp.write(f"# {banner} — {status} {datetime.now().isoformat(timespec='seconds')}\n")
    finally:
        log_fp.close()


def resolve_log_path(
    explicit: str | None, default_path: str | None, no_log: bool
) -> str | None:
    if no_log:
        return None
    if explicit and explicit.strip():
        return explicit.strip()
    return default_path
