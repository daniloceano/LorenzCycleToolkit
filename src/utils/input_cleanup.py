"""
Optional deletion of the input NetCDF file after a run (``--delete-downloaded-input``).

By default no input file is ever deleted, whether it was given by the user,
reused or downloaded by the run itself (--cdsapi or --ncar). With
``--delete-downloaded-input`` the input file is deleted only when every
condition below holds; otherwise it is kept and ``run_status.json``
(``input_cleanup``) says why:

- it was downloaded by this run: absent before the run, and reported as
  ``downloaded`` by the run timing (for --ncar only from the toolkit's own
  retrieval; a file written by another run is ``reused``);
- the run finished with status ``success`` (not partial, failed or
  interrupted: those stop before this step or keep the file);
- the file is unchanged since the download (same inode, size and
  modification time);
- the dataset opened from it has been closed.

The provenance record of an --ncar download is kept. Temporary files of the
download itself (daily CDS files, NCAR responses) are cleaned by the
download functions, as before.
"""

import os


def fingerprint(path):
    """(inode, size, modification time in ns) of a file, or None if it does not exist."""
    try:
        st = os.stat(path)
    except OSError:
        return None
    return (st.st_ino, st.st_size, st.st_mtime_ns)


def downloaded_input(args, existed_before, timing):
    """
    Fingerprint of the input file if this run downloaded it, else None.

    Args:
        args (argparse.Namespace): Command-line arguments.
        existed_before (bool): Whether the input file existed before the run.
        timing (RunTiming): Timing of the run (origin of the input file).
    """
    if existed_before or not (getattr(args, "cdsapi", False) or getattr(args, "ncar", False)):
        return None
    if timing.input_data != "downloaded":
        return None
    return fingerprint(args.infile)


def cleanup_input(args, run_status, downloaded, data=None, logger=None):
    """
    Delete the input file only when it is safe (see the module docstring),
    record the decision in ``run_status.input_cleanup`` and write the status.

    Returns:
        dict: ``{"requested", "path", "deleted", "reason"}``.
    """
    requested = bool(getattr(args, "delete_downloaded_input", False))
    path = args.infile
    deleted, reason = False, None
    if not requested:
        reason = "kept: input files are never deleted without --delete-downloaded-input"
    elif downloaded is None:
        reason = "kept: not downloaded by this run (local, pre-existing, reused or of uncertain origin)"
    elif run_status.status != "success":
        reason = f"kept: the run is {run_status.status}, not success"
    elif fingerprint(path) != downloaded:
        reason = "kept: the file is missing or changed since it was downloaded"
    else:
        if data is not None:
            try:
                data.close()
            except Exception:  # closing is best effort; unlinking an open file is safe on POSIX
                pass
        try:
            os.remove(path)
            deleted, reason = True, "deleted: downloaded by this run, which succeeded (--delete-downloaded-input)"
        except OSError as exc:
            reason = f"kept: deletion failed ({type(exc).__name__}: {exc})"

    run_status.input_cleanup = {"requested": requested, "path": path, "deleted": deleted, "reason": reason}
    run_status.write()
    if logger is not None:
        if deleted:
            logger.info(f"🗑️ Input file {path} deleted: downloaded by this run (--delete-downloaded-input)")
        elif requested:
            logger.warning(f"⚠️ --delete-downloaded-input: {path} {reason}")
    return run_status.input_cleanup
