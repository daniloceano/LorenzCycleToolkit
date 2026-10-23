"""
Run the LorenzCycleToolkit command line with the network blocked (test helper).

Usage: python run_cli_offline.py <path to lorenzcycletoolkit.py> [arguments...]

Every attempt to resolve a name or open a connection is refused and reported
on stderr as "LECTK-TEST: network access blocked: ...". The script runs as
``python lorenzcycletoolkit.py`` would (its directory first on sys.path), and
the exit status is that of the command line.
"""

import os
import runpy
import socket
import sys


def _blocked(call):
    def guard(*args, **kwargs):
        target = args[1] if call.startswith("connect") and len(args) > 1 else (args[0] if args else None)
        print(f"LECTK-TEST: network access blocked: {call}({target!r})", file=sys.stderr, flush=True)
        raise OSError(f"network access blocked by the test suite: {call}({target!r})")
    return guard


for owner, name in ((socket.socket, "connect"), (socket.socket, "connect_ex"),
                    (socket, "getaddrinfo"), (socket, "create_connection")):
    setattr(owner, name, _blocked(name))

script = os.path.abspath(sys.argv[1])
sys.argv = sys.argv[1:]
sys.path.insert(0, os.path.dirname(script))
runpy.run_path(script, run_name="__main__")
