"""Tests for utils.cli_log argument stripping."""

from utils.cli_log import strip_log_args_from_argv


def test_strip_log_with_file():
    out, parsed = strip_log_args_from_argv(["script.py", "--log", "run.log", "input.pdb"])
    assert out == ["script.py", "input.pdb"]
    assert parsed == {"log": "run.log"}


def test_strip_log_equals_file():
    out, parsed = strip_log_args_from_argv(["script.py", "--log=run.log", "input.pdb"])
    assert out == ["script.py", "input.pdb"]
    assert parsed == {"log": "run.log"}


def test_strip_bare_log():
    out, parsed = strip_log_args_from_argv(["script.py", "--log", "--verbose"])
    assert out == ["script.py", "--verbose"]
    assert parsed == {"log": "__AUTO__"}


def test_strip_bare_log_at_end():
    out, parsed = strip_log_args_from_argv(["script.py", "input.pdb", "--log"])
    assert out == ["script.py", "input.pdb"]
    assert parsed == {"log": "__AUTO__"}


def test_no_log_flag():
    out, parsed = strip_log_args_from_argv(["script.py", "input.pdb"])
    assert out == ["script.py", "input.pdb"]
    assert parsed == {"log": None}
