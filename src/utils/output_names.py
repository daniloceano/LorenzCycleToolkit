"""
Names of the results directory and of the files of one run.

Without ``-o/--outname`` every name is derived from the input file exactly as
before (``LEC_Results/<infile>_<method>/``, ``<infile>_<method>_results.csv``,
``<infile>_<method>_trackfile``, ``log.<infile>``). With ``--outname NAME``
the name replaces the input-file stem everywhere:

    LEC_Results/NAME_<method>/
        NAME_<method>_results.csv
        NAME_<method>_trackfile      (track and choose)
        log.NAME

Created by:
    Danilo Couto de Souza
    Universidade de São Paulo (USP)
    Instituto de Astornomia, Ciências Atmosféricas e Geociências
    São Paulo - Brazil

Contact:
    danilo.oceano@gmail.com
"""

import os

RESULTS_ROOT = "./LEC_Results/"


def validate_outname(name):
    """
    Check that an ``--outname`` value can be used as a directory and file stem.

    Args:
        name (str): The value given to ``-o/--outname``.

    Returns:
        str: The name, unchanged.

    Raises:
        ValueError: If the name is empty, is ``.`` or ``..``, or contains a
            path separator or a NUL character, so that it could not name a
            single directory inside ``LEC_Results/``.
    """
    if name is None or not name.strip():
        raise ValueError("the output name is empty")
    if name in (".", ".."):
        raise ValueError(f"the output name {name!r} is not allowed")
    separators = {"/", "\\", os.sep} | ({os.altsep} if os.altsep else set())
    if any(sep in name for sep in separators) or "\0" in name:
        raise ValueError(
            f"the output name {name!r} must not contain a path separator; "
            "it names one directory inside LEC_Results/"
        )
    return name


def _outname(args):
    return getattr(args, "outname", None)


def results_directory(args, method):
    """Directory of the results of this run (``LEC_Results/<stem>_<method>``)."""
    if _outname(args):
        stem = args.outname
    else:
        stem = "".join(args.infile.split("/")[-1].split(".nc"))
    return os.path.join(RESULTS_ROOT, stem + "_" + method)


def results_stem(args, method):
    """Stem of the main results CSV and of the trackfile (``<stem>_<method>``)."""
    if _outname(args):
        return f"{args.outname}_{method}"
    return f"{os.path.basename(args.infile).split('.nc')[0]}_{method}"


def log_file_name(args):
    """Name of the log file of the run (``log.<stem>``)."""
    if _outname(args):
        return f"log.{args.outname}"
    return f'log.{os.path.basename(args.infile).split(".")[0]}'
