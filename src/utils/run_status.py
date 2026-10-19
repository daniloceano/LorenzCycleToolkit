"""
Status of one LorenzCycleToolkit run, written to ``run_status.json``.

The status file is written next to the results (``LEC_Results/<case>/``) by
every command-line run, with or without ``--plots``. It is an additional file:
the CSV outputs are not changed.

States of a product (``products[].status``):
    success        requested and every expected file was written by this run
    failed         requested and not (or not completely) generated
    unavailable    not applicable to these data for a known, approved reason
    not_requested  optional product not requested (no ``--plots``)
While the run is still going, requested products are ``pending``.

Global state (``status``):
    running   the run has not finished (left on disk if it is killed)
    success   scientific computation complete and every requested product generated
    partial   scientific computation complete, at least one requested figure failed
    failed    the scientific computation (or the run) failed

Exit status of the command line (``exit_code``): 0 for success, 3 for partial
(0 with ``--allow-partial-plots``; the status stays ``partial``), 1 when the
scientific computation fails (the original exception propagates).

The scientific outputs are required; the figures and the cyclophaser periods
are optional diagnostics. A figure call is the only thing guarded: an error in
the scientific computation is never turned into a partial result.
"""

import contextlib
import datetime
import json
import os
import sys
from pathlib import Path

SCHEMA_VERSION = 1

EXIT_SUCCESS = 0
EXIT_FAILED = 1
EXIT_PARTIAL = 3

STATUS_FILE = "run_status.json"

# Optional products requested with --plots, in the order they are generated.
FIXED_FIGURES = ("timeseries", "box_limits_map", "boxplots", "hovmoller", "lec_cycle")
MOVING_FIGURES = (
    "domain_maps", "track_map", "zeta_hgt_timeseries", "timeseries", "hovmoller",
    "boxplots", "periods", "lec_cycle", "lec_periods", "lps",
)
DIAGNOSTICS = ("periods",)

# Wording of the log messages, kept as in the previous fixed-framework messages.
DESCRIPTIONS = {
    "results_csv": "Main results CSV",
    "vertical_level_csvs": "Vertical-level CSVs",
    "trackfile": "System track file",
    "domain_maps": "Domain map",
    "box_limits_map": "Box limits plot",
    "track_map": "Track map",
    "zeta_hgt_timeseries": "Vorticity and geopotential height time series",
    "timeseries": "Time series plot",
    "hovmoller": "Hovmöller plot",
    "boxplots": "Boxplot terms",
    "periods": "Life-cycle periods (cyclophaser)",
    "lec_cycle": "Lorenz cycle plot",
    "lec_periods": "Lorenz cycle plots for the life-cycle periods",
    "lps": "Lorenz Phase Space plots",
}


def _now():
    return datetime.datetime.now(datetime.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def _error(exc):
    return {"type": type(exc).__name__, "message": str(exc)}


def _version():
    try:
        from importlib.metadata import version

        return version("lorenzcycletoolkit")
    except Exception:  # not installed (source tree)
        return None


class RunStatus:
    """
    In-memory record of a run, written atomically to ``run_status.json``.

    Args:
        args (argparse.Namespace): Parsed command-line arguments.
        results_directory (str): Directory of the results of this run.
        app_logger (logging.Logger): Logger of the run (optional).
        argv (list): Command line to record (defaults to ``sys.argv``).
    """

    def __init__(self, args, results_directory, app_logger=None, argv=None):
        self.results_directory = Path(results_directory)
        self.path = self.results_directory / STATUS_FILE
        self.logger = app_logger
        self.plots_requested = bool(getattr(args, "plots", False))
        self.allow_partial_plots = bool(getattr(args, "allow_partial_plots", False))
        if getattr(args, "fixed", False):
            self.method = "fixed"
        elif getattr(args, "track", False):
            self.method = "track"
        else:
            self.method = "choose"
        self.command = {
            "argv": list(sys.argv if argv is None else argv),
            "method": self.method,
            "infile": getattr(args, "infile", None),
        }
        self.started = _now()
        self.finished = None
        self.scientific = {"status": "running", "error": None}
        self.run_error = None
        self.interrupted = False
        self.log_file = None
        if app_logger is not None:
            for handler in app_logger.handlers:
                filename = getattr(handler, "baseFilename", None)
                if filename:
                    self.log_file = self._relative(filename)

        self.products = {}
        required = ["results_csv", "vertical_level_csvs"]
        if self.method != "fixed":
            required.append("trackfile")
        for name in required:
            self.products[name] = self._new(name, "scientific", True, "pending")
        figures = FIXED_FIGURES if self.method == "fixed" else MOVING_FIGURES
        for name in figures:
            kind = "diagnostic" if name in DIAGNOSTICS else "figure"
            state = "pending" if self.plots_requested else "not_requested"
            self.products[name] = self._new(name, kind, False, state)

    # ------------------------------------------------------------------ helpers
    @staticmethod
    def _new(name, kind, required, status):
        return {
            "name": name, "kind": kind, "required": required, "status": status,
            "expected": [], "written": [], "missing": [], "reason": None,
            "error": None, "depends_on": [], "consequences": [], "failures": [],
        }

    def _relative(self, path):
        path = Path(path)
        try:
            return path.resolve().relative_to(self.results_directory.resolve()).as_posix()
        except ValueError:
            return str(path)

    def _log(self, level, message):
        if self.logger is not None:
            getattr(self.logger, level)(message)

    def _scan(self, watch):
        """Return {relative path: (mtime_ns, size, inode)} of the watched files and directories."""
        state = {}
        for rel in watch:
            path = self.results_directory / rel
            if path.is_dir():
                entries = [p for p in path.iterdir() if p.is_file()]
            elif path.is_file():
                entries = [path]
            else:
                continue
            for entry in entries:
                st = entry.stat()
                state[self._relative(entry)] = (st.st_mtime_ns, st.st_size, st.st_ino)
        return state

    # ---------------------------------------------------------------- recording
    def record_scientific_outputs(self, outputs):
        """
        Record the required scientific outputs once they have been written.

        Args:
            outputs (dict): {product name: list of file paths}.

        Raises:
            RuntimeError: If a declared output file does not exist; this is a
                failure of the scientific computation, not of a figure.
        """
        for name, paths in outputs.items():
            product = self.products[name]
            rel = [self._relative(p) for p in paths]
            product["expected"] = rel
            missing = [r for r in rel if not (self.results_directory / r).is_file()]
            if missing:
                product.update(status="failed", missing=missing,
                               reason="required output not written: " + ", ".join(missing))
                raise RuntimeError(f"{DESCRIPTIONS[name]} not written: {', '.join(missing)}")
            product.update(status="success", written=rel)
        pending = [n for n, p in self.products.items() if p["required"] and p["status"] != "success"]
        if not pending:
            self.scientific = {"status": "success", "error": None}

    def run_product(self, name, func, *args, expected=(), watch=None, step=None, **kwargs):
        """
        Generate an optional product, recording whether it was really written.

        Only the call ``func(*args, **kwargs)`` is guarded. The product is a
        success when the call returns and every expected file was written by
        it (created or modified); otherwise it is failed, with the exception or
        the missing files as the reason. Calls with the same name and a
        ``step`` (one per time step) are accumulated into one product, which
        stays pending until ``close_product`` is called after the last step.

        Returns:
            bool: True if this call generated all its expected files.
        """
        product = self.products[name]
        expected = [self._relative(self.results_directory / e) for e in expected]
        if watch is None:
            watch = sorted({str(Path(e).parent) for e in expected} | set(expected))
        before = self._scan(watch)
        error = None
        try:
            func(*args, **kwargs)
        except Exception as exc:  # a figure failure must not stop the run
            error = exc
        after = self._scan(watch)
        written = sorted(p for p, s in after.items() if before.get(p) != s)
        missing = [e for e in expected if e not in written]

        product["expected"].extend(e for e in expected if e not in product["expected"])
        product["written"].extend(w for w in written if w not in product["written"])
        product["missing"].extend(m for m in missing if m not in product["missing"])
        label = f" ({step})" if step else ""
        if error is None and not missing:
            # a product made of time steps is complete only after the last one
            if product["status"] == "pending" and step is None:
                product["status"] = "success"
            self._log("info", f"  ✅ {DESCRIPTIONS[name]}{label} generated")
            return True

        if error is not None:
            cause = f"{type(error).__name__}: {error}"
        else:
            cause = "the call returned without writing " + ", ".join(missing)
        failure = {"step": step, "missing": missing, "cause": cause}
        if error is not None:
            failure["error"] = _error(error)
        product["failures"].append(failure)
        if product["error"] is None and error is not None:
            product["error"] = _error(error)
        product["status"] = "failed"
        n = len(product["failures"])
        first = product["failures"][0]["cause"]
        product["reason"] = first if n == 1 else f"{n} calls failed; first: {first}"
        self._log("error", f"  ❌ {DESCRIPTIONS[name]}{label} not generated: {cause}; "
                           f"expected {', '.join(missing) or ', '.join(expected)}")
        return False

    def close_product(self, name):
        """Mark a product made of time steps as complete after its last step."""
        product = self.products[name]
        if product["status"] == "pending":
            product["status"] = "success"

    def fail_products(self, names, exc, reason):
        """Mark requested products as failed without calling them (e.g. import error)."""
        for name in names:
            product = self.products[name]
            product.update(status="failed", error=_error(exc),
                           reason=f"{reason}: {type(exc).__name__}: {exc}")
            self._log("error", f"  ❌ {DESCRIPTIONS[name]} not generated: {product['reason']}")

    def dependency_failed(self, name, dependency, expected=()):
        """Mark a product as failed because the product it needs is not available."""
        dep = self.products[dependency]
        product = self.products[name]
        product.update(
            status="failed", depends_on=[dependency],
            expected=[self._relative(self.results_directory / e) for e in expected],
            reason=f"not generated: requires '{dependency}', which {dep['status']} ({dep['reason']})",
        )
        if name not in dep["consequences"]:
            dep["consequences"].append(name)
        self._log("error", f"  ❌ {DESCRIPTIONS[name]} not generated: requires "
                           f"{DESCRIPTIONS[dependency].lower()}, which {dep['status']}")

    def succeeded(self, name):
        return self.products[name]["status"] == "success"

    def figures_ok(self):
        """True when no requested optional product failed."""
        return not any(p["status"] == "failed" for p in self.products.values() if not p["required"])

    def failed_figures(self):
        return [n for n, p in self.products.items() if not p["required"] and p["status"] == "failed"]

    # ------------------------------------------------------------ global state
    @property
    def status(self):
        if self.run_error is not None or self.scientific["status"] == "failed":
            return "failed"
        if self.finished is None:
            return "running"
        return "success" if self.figures_ok() else "partial"

    @property
    def exit_code(self):
        status = self.status
        if status == "success":
            return EXIT_SUCCESS
        if status == "partial":
            return EXIT_SUCCESS if self.allow_partial_plots else EXIT_PARTIAL
        if status == "failed":
            return None if self.interrupted else EXIT_FAILED
        return None

    def record_exception(self, exc):
        """Record an exception that ends the run (the caller re-raises it)."""
        self.finished = _now()
        self.run_error = _error(exc)
        if not isinstance(exc, Exception):
            self.interrupted = True
        if self.scientific["status"] != "success":
            self.scientific = {"status": "failed", "error": _error(exc)}
        for product in self.products.values():
            if product["status"] == "pending":
                product["status"] = "failed"
                product["reason"] = (f"not generated: the run stopped before this product "
                                     f"({type(exc).__name__}: {exc})")

    def finish(self):
        """
        Close the record after a complete run and write it.

        Raises:
            RuntimeError: If the required scientific outputs were never
                recorded; the run is then failed, never complete.
        """
        if self.scientific["status"] != "success":
            pending = [n for n, p in self.products.items() if p["required"] and p["status"] != "success"]
            raise RuntimeError("required scientific outputs were not recorded: " + ", ".join(pending))
        self.finished = _now()
        for product in self.products.values():
            if product["status"] == "pending":
                product["status"] = "failed"
                product["reason"] = "not generated: the product was never attempted"
        self.write()

    @contextlib.contextmanager
    def guard(self):
        """Record any exception leaving the block, write the status and re-raise."""
        try:
            yield self
        except BaseException as exc:
            self.record_exception(exc)
            try:
                self.write()
            except Exception as write_error:
                self._log("error", f"❌ Could not write {self.path}: {write_error}")
            raise

    # ------------------------------------------------------------------ output
    def to_dict(self):
        products = list(self.products.values())
        summary = {s: sum(p["status"] == s for p in products)
                   for s in ("success", "failed", "unavailable", "not_requested", "pending")}
        return {
            "schema_version": SCHEMA_VERSION,
            "lectk_version": _version(),
            "status": self.status,
            "exit_code": self.exit_code,
            "allow_partial_plots": self.allow_partial_plots,
            "plots_requested": self.plots_requested,
            "interrupted": self.interrupted,
            "error": self.run_error,
            "command": self.command,
            "started_utc": self.started,
            "finished_utc": self.finished,
            "scientific": self.scientific,
            "summary": summary,
            "log_file": self.log_file,
            "products": products,
        }

    def write(self):
        """Write run_status.json atomically (temporary file, then rename)."""
        tmp = self.path.with_name(f".{STATUS_FILE}.{os.getpid()}.tmp")
        try:
            with open(tmp, "w", encoding="utf-8") as handle:
                json.dump(self.to_dict(), handle, indent=2, ensure_ascii=False)
                handle.write("\n")
                handle.flush()
                os.fsync(handle.fileno())
            os.replace(tmp, self.path)
        except BaseException:
            with contextlib.suppress(OSError):
                tmp.unlink()
            raise

    def final_message(self, framework, seconds):
        """Closing log line: the success message only when it is true."""
        if self.status == "success":
            return f"🎉 Analysis complete! {framework} framework ran in {seconds:.2f} seconds"
        failed = self.failed_figures()
        requested = [n for n, p in self.products.items()
                     if not p["required"] and p["status"] != "not_requested"]
        exit_note = ("exit status 0 because --allow-partial-plots was given"
                     if self.allow_partial_plots else f"exit status {EXIT_PARTIAL}")
        return (f"⚠️ Analysis complete with PARTIAL output: {len(failed)} of {len(requested)} requested "
                f"figure products failed ({', '.join(failed)}); see {self.path} ({exit_note}). "
                f"{framework} framework ran in {seconds:.2f} seconds")
