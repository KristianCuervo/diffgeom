"""
Run every example end to end.

The examples are the library's real test suite: each one is a worked problem
whose expected answers are written into its own comments. If an example still
runs, the API it exercises still works. This caught the API drift that left a
third of them crashing.
"""
import pathlib
import subprocess
import sys

import pytest

REPO_ROOT = pathlib.Path(__file__).resolve().parent.parent
EXAMPLES = sorted((REPO_ROOT / "examples").rglob("*.py"))

# The whole suite runs in well under a minute; this is only a hang guard.
TIMEOUT_SECONDS = 600


def _example_id(path: pathlib.Path) -> str:
    return str(path.relative_to(REPO_ROOT / "examples"))


assert EXAMPLES, "no examples found -- has the layout changed?"


@pytest.mark.parametrize("example", EXAMPLES, ids=_example_id)
def test_example_runs(example: pathlib.Path):
    result = subprocess.run(
        [sys.executable, str(example)],
        capture_output=True,
        text=True,
        timeout=TIMEOUT_SECONDS,
        cwd=REPO_ROOT,
    )
    assert result.returncode == 0, (
        f"{_example_id(example)} exited {result.returncode}\n"
        f"--- stdout ---\n{result.stdout}\n"
        f"--- stderr ---\n{result.stderr}"
    )
