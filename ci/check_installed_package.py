"""
Check, before the tests, that the installed LorenzCycleToolkit is the code of
the checkout: every Python file of the installed distribution must equal, byte
for byte, the same file of the checkout. Prints the version, the location of
the distribution and the number of files checked; any difference stops the
script with an error.

Usage: python ci/check_installed_package.py
"""

import hashlib
import sys
from importlib import metadata
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main():
    dist = metadata.distribution("lorenzcycletoolkit")
    files = [f for f in dist.files or [] if f.suffix == ".py" and ".." not in f.parts]
    problems = []
    for f in files:
        installed, checkout = Path(dist.locate_file(f)), ROOT / f
        if not checkout.is_file():
            problems.append(f"not in the checkout: {f}")
        elif sha256(installed) != sha256(checkout):
            problems.append(f"differs from the checkout: {f}")
    print(f"lorenzcycletoolkit {dist.version} installed at {dist.locate_file('')}; "
          f"{len(files)} Python files compared with the checkout {ROOT}")
    if not files or problems:
        sys.exit("The installed package is not the checked-out code: " + ("; ".join(problems) or "no files"))
    print("installed package = checkout")


if __name__ == "__main__":
    main()
