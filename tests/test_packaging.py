"""
Declared runtime dependencies (Wave 5: P-A1).

Every third-party module imported by the toolkit (``src/`` and
``lorenzcycletoolkit.py``) must come from a distribution listed in
``[project] dependencies`` of ``pyproject.toml``; otherwise ``pip install``
gives an installation whose command line fails (Wave 5, A1: ``No module
named 'cmocean'``). Import names are mapped to their distributions by
``importlib.metadata.packages_distributions``.
"""

import ast
import re
import sys
import tomllib
from importlib.metadata import packages_distributions, requires
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
OWN_MODULES = {"src", "lorenzcycletoolkit", "tests"}


def _third_party_imports():
    sources = sorted((REPO_ROOT / "src").rglob("*.py")) + [REPO_ROOT / "lorenzcycletoolkit.py"]
    names = set()
    for path in sources:
        for node in ast.walk(ast.parse(path.read_text(encoding="utf-8"))):
            if isinstance(node, ast.Import):
                names |= {alias.name.split(".")[0] for alias in node.names}
            elif isinstance(node, ast.ImportFrom) and node.level == 0 and node.module:
                names.add(node.module.split(".")[0])
    return sorted(n for n in names if n not in sys.stdlib_module_names and n not in OWN_MODULES)


def _normalise(name):
    return re.sub(r"[-_.]+", "-", name).lower()


def _names(dependencies):
    return {_normalise(re.split(r"[<>=!~;\[ ]", d, maxsplit=1)[0]) for d in dependencies}


def _dependencies():
    pyproject = tomllib.loads((REPO_ROOT / "pyproject.toml").read_text(encoding="utf-8"))
    return pyproject["project"]["dependencies"]


def test_every_imported_distribution_is_declared():
    declared = _names(_dependencies())
    providers = packages_distributions()
    missing = {}
    for module in _third_party_imports():
        distributions = {_normalise(d) for d in providers.get(module, [])}
        assert distributions, f"{module} is imported but no installed distribution provides it"
        if not distributions & declared:
            missing[module] = sorted(distributions)
    assert missing == {}, f"imported but not declared in pyproject.toml: {missing}"


def test_cyclophaser_declared_at_the_tested_version():
    """The compatibility layer accepts one cyclophaser version only; the metadata must ask for that one."""
    from src.utils.cyclophaser_compat import TESTED_VERSION

    pins = [d for d in _dependencies() if _normalise(d).startswith("cyclophaser")]
    assert pins == [f"cyclophaser=={TESTED_VERSION}"]


def test_installed_metadata_matches_pyproject():
    """The installed distribution declares the dependencies of pyproject.toml (a stale wheel would not)."""
    assert _names(requires("lorenzcycletoolkit") or []) == _names(_dependencies())
