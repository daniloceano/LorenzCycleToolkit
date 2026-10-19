"""
Shared pytest fixtures and session checks for the LorenzCycleToolkit test suite.

The integration tests run the toolkit in a private temporary working directory,
so the personal configuration files a user keeps in ``inputs/`` of the
repository (``namelist``, ``track``, ``box_limits``) are never read or
overwritten by the suite.  A session check also fails the run if anything
under ``inputs/`` or ``samples/`` of the repository changed while the tests
were running.
"""

import hashlib
import os
import shutil
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[1]
PROTECTED_DIRS = ("inputs", "samples")


def _digest_protected(root):
    """Return {relative path: sha256} for every file under the protected dirs."""
    digests = {}
    for name in PROTECTED_DIRS:
        base = root / name
        if not base.is_dir():
            continue
        for path in sorted(base.rglob("*")):
            if path.is_file():
                rel = path.relative_to(root).as_posix()
                digests[rel] = hashlib.sha256(path.read_bytes()).hexdigest()
    return digests


def pytest_sessionstart(session):
    session.config._lectk_protected_before = _digest_protected(REPO_ROOT)


def pytest_sessionfinish(session, exitstatus):
    before = getattr(session.config, "_lectk_protected_before", None)
    if before is None:
        return
    after = _digest_protected(REPO_ROOT)
    changed = sorted(k for k in set(before) | set(after) if before.get(k) != after.get(k))
    if changed:
        session.config._lectk_protected_changed = changed
        session.exitstatus = pytest.ExitCode.TESTS_FAILED


def pytest_terminal_summary(terminalreporter, exitstatus, config):
    changed = getattr(config, "_lectk_protected_changed", None)
    if changed:
        terminalreporter.section("protected repository files changed", sep="=", red=True)
        for rel in changed:
            terminalreporter.write_line(f"modified, created or removed during the tests: {rel}")


@pytest.fixture
def lectk_workdir(tmp_path, monkeypatch):
    """
    Prepare a private working directory for one toolkit run and chdir into it.

    The returned function copies the requested templates from ``inputs/`` of
    the repository to ``inputs/namelist``, ``inputs/box_limits`` and
    ``inputs/track`` of the working directory, links ``samples/`` read-only
    by path, and returns the working directory.  Results are written to
    ``LEC_Results/`` inside it, so each test owns its inputs and outputs.
    """

    def _prepare(namelist, box_limits=None, track=None):
        workdir = tmp_path / "work"
        (workdir / "inputs").mkdir(parents=True)
        os.symlink(REPO_ROOT / "samples", workdir / "samples")
        copies = {"namelist": namelist, "box_limits": box_limits, "track": track}
        for target, template in copies.items():
            if template is not None:
                shutil.copyfile(REPO_ROOT / "inputs" / template, workdir / "inputs" / target)
        monkeypatch.chdir(workdir)
        return workdir

    return _prepare
