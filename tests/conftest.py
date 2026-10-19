"""
Shared pytest fixtures and session checks for the LorenzCycleToolkit test suite.

The integration tests run the toolkit in a private temporary working directory,
so the personal configuration files a user keeps in ``inputs/`` of the
repository (``namelist``, ``track``, ``box_limits``) are never read or
overwritten by the suite.  A session check also fails the run if anything
under ``inputs/`` or ``samples/`` of the repository changed while the tests
were running, and, after all tests have run, if the log file of any run
holds messages of another run or of a later test.

The suite runs offline: name resolution and connections are blocked for
the whole session and every attempt is recorded. An attempt from a test
not marked ``expect_network_attempts`` fails the run. The tests that draw
maps use a verified Natural Earth copy (``naturalearth_cache``).
"""

import hashlib
import logging
import os
import shutil
import socket
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[1]
PROTECTED_DIRS = ("inputs", "samples")
NATURALEARTH_MANIFEST = REPO_ROOT / "tests" / "data" / "naturalearth_manifest.tsv"

# ------------------------------------------------------------------ network
NETWORK_ATTEMPTS = []  # (test node id, call, target)
_CURRENT_TEST = {"nodeid": None}
_SOCKET_FUNCTIONS = (
    (socket.socket, "connect"),
    (socket.socket, "connect_ex"),
    (socket, "getaddrinfo"),
    (socket, "create_connection"),
)


class NetworkAccessBlocked(OSError):
    """Raised by the session network guard of the test suite."""


def _blocked(call):
    def guard(*args, **kwargs):
        target = args[1] if call.startswith("connect") and len(args) > 1 else (args[0] if args else None)
        NETWORK_ATTEMPTS.append((_CURRENT_TEST["nodeid"], call, repr(target)))
        raise NetworkAccessBlocked(f"network access blocked by the test suite: {call}({target!r})")
    return guard


def pytest_configure(config):
    config.addinivalue_line(
        "markers",
        "expect_network_attempts: the test checks that an attempt to reach the network is blocked",
    )
    config._lectk_socket_saved = [(owner, name, getattr(owner, name)) for owner, name in _SOCKET_FUNCTIONS]
    for owner, name in _SOCKET_FUNCTIONS:
        setattr(owner, name, _blocked(name))


def pytest_unconfigure(config):
    for owner, name, original in getattr(config, "_lectk_socket_saved", []):
        setattr(owner, name, original)


def pytest_collection_modifyitems(session, config, items):
    config._lectk_network_expected = {
        item.nodeid for item in items if item.get_closest_marker("expect_network_attempts")
    }


@pytest.fixture(autouse=True)
def _network_guard_scope(request):
    _CURRENT_TEST["nodeid"] = request.node.nodeid
    yield
    _CURRENT_TEST["nodeid"] = None


@pytest.fixture
def network_attempts(request):
    """The network attempts blocked during the current test."""
    class _View:
        def __iter__(self):
            return iter([a for a in NETWORK_ATTEMPTS if a[0] == request.node.nodeid])

        def __len__(self):
            return len(list(iter(self)))

    return _View()


# ------------------------------------------------------------ Natural Earth
@pytest.fixture(scope="session")
def naturalearth_cache():
    """The verified Natural Earth copy named by CARTOPY_DATA_DIR (explicit error if absent)."""
    from tests.naturalearth import cache_from_environment

    return cache_from_environment(NATURALEARTH_MANIFEST)


def _fresh_cartopy(monkeypatch, pre_existing, data_dir):
    """Point Cartopy to the given directories and forget geometries read before."""
    import cartopy
    import cartopy.feature

    monkeypatch.setitem(cartopy.config, "pre_existing_data_dir", Path(pre_existing))
    monkeypatch.setitem(cartopy.config, "data_dir", Path(data_dir))
    saved = dict(cartopy.feature._NATURAL_EARTH_GEOM_CACHE)
    cartopy.feature._NATURAL_EARTH_GEOM_CACHE.clear()
    return saved


@pytest.fixture
def naturalearth_reads(naturalearth_cache, monkeypatch, tmp_path):
    """
    Cartopy reads Natural Earth only from the verified copy (no download
    directory content, no geometries kept from earlier tests); returns the
    (resolution, category, name) sets requested during the test.
    """
    import cartopy.feature
    import cartopy.io.shapereader as shapereader

    saved = _fresh_cartopy(monkeypatch, naturalearth_cache, tmp_path / "cartopy_download_dir")
    requested = []
    original = shapereader.natural_earth

    def spy(resolution="110m", category="physical", name="coastline"):
        requested.append((resolution, category, name))
        return original(resolution=resolution, category=category, name=name)

    monkeypatch.setattr(shapereader, "natural_earth", spy)
    yield requested
    cartopy.feature._NATURAL_EARTH_GEOM_CACHE.clear()
    cartopy.feature._NATURAL_EARTH_GEOM_CACHE.update(saved)


@pytest.fixture
def naturalearth_absent(monkeypatch, tmp_path):
    """
    No Natural Earth data available: empty pre-existing and download
    directories and no geometries kept from earlier tests. Any download
    attempt meets the session network guard.
    """
    import cartopy.feature

    empty = tmp_path / "naturalearth_absent"
    empty.mkdir()
    saved = _fresh_cartopy(monkeypatch, empty, tmp_path / "cartopy_download_dir")
    yield empty
    cartopy.feature._NATURAL_EARTH_GEOM_CACHE.clear()
    cartopy.feature._NATURAL_EARTH_GEOM_CACHE.update(saved)


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


def _log_isolation_problems(workdirs, sizes_at_test_end):
    """
    Return the problems found in the run logs written under ``workdirs``.

    A run log must hold one run ("Starting LEC analysis" once), name only its
    own input file, end with the run (nothing after the last "Analysis
    complete" line) and must not grow after the test that wrote it has
    finished.  Lines added later come from later tests through the same
    logger, which the per-test checks cannot see because they run first.
    """
    problems = []
    for log, size in sorted(sizes_at_test_end.items()):
        if log.stat().st_size != size:
            grown = log.stat().st_size - size
            problems.append(f"{'/'.join(log.parts[-4:])}: {grown} bytes written after its test ended")
    for workdir in workdirs:
        for log in sorted(Path(workdir).glob("LEC_Results/*/log.*")):
            lines = log.read_text(encoding="utf-8", errors="replace").splitlines()
            where = f"{Path(workdir).parent.name}/{log.relative_to(workdir)}"
            starts = sum("Starting LEC analysis" in line for line in lines)
            if starts != 1:
                problems.append(f"{where}: {starts} runs logged, expected 1")
            stem = log.name[len("log."):]
            infiles = [line.split("infile='", 1)[1].split("'", 1)[0] for line in lines if "infile='" in line]
            foreign = sorted({f for f in infiles if Path(f).name != f"{stem}.nc"})
            if foreign:
                problems.append(f"{where}: arguments of other runs logged: {foreign}")
            complete = [i for i, line in enumerate(lines) if "Analysis complete" in line]
            if complete and complete[-1] != len(lines) - 1:
                extra = len(lines) - 1 - complete[-1]
                problems.append(f"{where}: {extra} lines written after the run completed")
    return problems


def pytest_sessionfinish(session, exitstatus):
    before = getattr(session.config, "_lectk_protected_before", None)
    if before is None:
        return
    after = _digest_protected(REPO_ROOT)
    changed = sorted(k for k in set(before) | set(after) if before.get(k) != after.get(k))
    if changed:
        session.config._lectk_protected_changed = changed
        session.exitstatus = pytest.ExitCode.TESTS_FAILED
    problems = _log_isolation_problems(
        getattr(session.config, "_lectk_workdirs", []),
        getattr(session.config, "_lectk_log_sizes", {}),
    )
    if problems:
        session.config._lectk_log_problems = problems
        session.exitstatus = pytest.ExitCode.TESTS_FAILED
    expected = getattr(session.config, "_lectk_network_expected", set())
    unexpected = [a for a in NETWORK_ATTEMPTS if a[0] not in expected]
    if unexpected:
        session.config._lectk_network_unexpected = unexpected
        session.exitstatus = pytest.ExitCode.TESTS_FAILED


def pytest_terminal_summary(terminalreporter, exitstatus, config):
    changed = getattr(config, "_lectk_protected_changed", None)
    if changed:
        terminalreporter.section("protected repository files changed", sep="=", red=True)
        for rel in changed:
            terminalreporter.write_line(f"modified, created or removed during the tests: {rel}")
    problems = getattr(config, "_lectk_log_problems", None)
    if problems:
        terminalreporter.section("run logs not isolated", sep="=", red=True)
        for problem in problems:
            terminalreporter.write_line(problem)
    unexpected = getattr(config, "_lectk_network_unexpected", None)
    if unexpected:
        terminalreporter.section("network access attempted", sep="=", red=True)
        for nodeid, call, target in unexpected:
            terminalreporter.write_line(f"{nodeid or '(outside a test)'}: {call}({target}) blocked")
    terminalreporter.write_line(
        f"network attempts blocked by the session guard: {len(NETWORK_ATTEMPTS)} "
        f"({len(unexpected or [])} unexpected)"
    )


@pytest.fixture(autouse=True)
def _close_lectk_log_handlers():
    """
    Close the log handlers of the toolkit logger after every test.

    A run configures the "lorenzcycletoolkit" logger to write to its own log
    file; without this, messages emitted by later tests through the same
    logger would be appended to the log file of the last integration case.
    """
    yield
    logger = logging.getLogger("lorenzcycletoolkit")
    for handler in list(logger.handlers):
        logger.removeHandler(handler)
        handler.close()


@pytest.fixture
def lectk_workdir(tmp_path, monkeypatch, request):
    """
    Prepare a private working directory for one toolkit run and chdir into it.

    The returned function copies the requested templates from ``inputs/`` of
    the repository to ``inputs/namelist``, ``inputs/box_limits`` and
    ``inputs/track`` of the working directory, links ``samples/`` read-only
    by path, and returns the working directory.  Results are written to
    ``LEC_Results/`` inside it, so each test owns its inputs and outputs.
    """

    workdir = tmp_path / "work"

    def _prepare(namelist, box_limits=None, track=None):
        (workdir / "inputs").mkdir(parents=True)
        if not hasattr(request.config, "_lectk_workdirs"):
            request.config._lectk_workdirs = []
        request.config._lectk_workdirs.append(workdir)
        os.symlink(REPO_ROOT / "samples", workdir / "samples")
        copies = {"namelist": namelist, "box_limits": box_limits, "track": track}
        for target, template in copies.items():
            if template is not None:
                shutil.copyfile(REPO_ROOT / "inputs" / template, workdir / "inputs" / target)
        monkeypatch.chdir(workdir)
        return workdir

    yield _prepare

    # Record the size of every run log at the end of the test; the session
    # check fails if any of them grows afterwards.
    if not hasattr(request.config, "_lectk_log_sizes"):
        request.config._lectk_log_sizes = {}
    for log in workdir.glob("LEC_Results/*/log.*"):
        request.config._lectk_log_sizes[log] = log.stat().st_size
