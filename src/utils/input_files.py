# **************************************************************************** #
#                                                                              #
#                                                         :::      ::::::::    #
#    input_files.py                                     :+:      :+:    :+:    #
#                                                     +:+ +:+         +:+      #
#    By: daniloceano <danilo.oceano@gmail.com>      +#+  +:+       +#+         #
#                                                 +#+#+#+#+#+   +#+            #
#    Created: 2026/08/18 10:00:00 by daniloceano       #+#    #+#              #
#    Updated: 2026/08/18 10:00:00 by daniloceano      ###   ########.fr        #
#                                                                              #
# **************************************************************************** #

"""
Resolution and parsing of the text input files kept in ``inputs/``.

The working copies of ``inputs/namelist`` and ``inputs/box_limits`` are not
tracked by git, so that editing them for a run does not show up as a change to
the repository. Only the ``.default`` companions are distributed, and they are
used whenever no working copy is present.

This module has no intra-package imports, so both ``src/utils/select_area.py``
and ``src/utils/tools.py`` -- which import each other -- can use it.

Created by:
    Danilo Couto de Souza
    Universidade de São Paulo (USP)
    Instituto de Astornomia, Ciências Atmosféricas e Geociências
    São Paulo - Brazil

Contact:
    danilo.oceano@gmail.com
"""

import logging
import os

import pandas as pd


def resolve_input_file(path: str) -> str:
    """
    Return ``path``, or its shipped ``.default`` companion when ``path`` is absent.

    Args:
        path (str): Path to a working input file, e.g. ``inputs/namelist``.

    Returns:
        str: The path to use. When neither the working copy nor the
        ``.default`` companion exists, ``path`` is returned unchanged so that
        the caller reports the missing file the user actually asked for.
    """
    if os.path.exists(path):
        return path
    fallback = f"{path}.default"
    if os.path.exists(fallback):
        return fallback
    return path


def read_box_limits(path: str, app_logger=None):
    """
    Read a box-limits file, resolving the ``.default`` fallback first.

    Args:
        path (str): Path to the box-limits file, e.g. ``inputs/box_limits``.
        app_logger (logging.Logger, optional): Logger used to report a missing,
            empty or malformed file. Falls back to the root logger.

    Returns:
        tuple: ``(resolved_path, min_lon, max_lon, min_lat, max_lat)``, with the
        four limits as raw floats in the file's own longitude convention. They
        are returned unvalidated: :func:`src.utils.longitude.canonicalize_box`
        is the single place where the longitude pair is checked and expressed
        in the dataset's convention.

    Raises:
        FileNotFoundError: If neither the working copy nor the ``.default``
            companion exists.
        pandas.errors.EmptyDataError: If the file is empty.
        ValueError: If any of the four required rows is missing.
    """
    resolved = resolve_input_file(path)
    log_error = app_logger.error if app_logger is not None else logging.error

    try:
        dfbox = pd.read_csv(resolved, header=None, delimiter=";", index_col=0)
    except FileNotFoundError:
        log_error("❌ Box limits file not found!")
        log_error("\n" + "=" * 70)
        log_error("📁 BOX LIMITS FILE NOT FOUND")
        log_error("=" * 70)
        log_error(f"Looking for: {os.path.abspath(resolved)}")
        log_error(f"Current directory: {os.getcwd()}")
        log_error("\n💡 User Solutions:")
        log_error("   1. Create a box_limits file with the domain boundaries")
        log_error("   2. Use the shipped default: inputs/box_limits.default")
        log_error("   3. Specify a custom file with: --box_limits <path>")
        log_error("\n📝 Expected format:")
        log_error("   min_lon;-60")
        log_error("   max_lon;-30")
        log_error("   min_lat;-50")
        log_error("   max_lat;-20")
        log_error("=" * 70 + "\n")
        raise FileNotFoundError(
            f"Box limits file not found: {os.path.abspath(resolved)}. "
            f"Create one or use --box_limits to specify a path."
        )
    except pd.errors.EmptyDataError:
        log_error("❌ Box limits file is empty!")
        log_error(f"File: {os.path.abspath(resolved)}")
        raise pd.errors.EmptyDataError(f"Box limits file is empty: {resolved}")
    except Exception as e:
        log_error("❌ Error reading box_limits file!")
        log_error(f"File: {os.path.abspath(resolved)}")
        log_error(f"Error: {type(e).__name__}: {e}")
        log_error("\n💡 Check the file format (CSV with ';' delimiter)")
        raise

    required = ("min_lon", "max_lon", "min_lat", "max_lat")
    missing = [key for key in required if key not in dfbox.index]
    if missing:
        log_error("❌ Box limits file is missing required fields!")
        log_error("\n" + "=" * 70)
        log_error("📋 MISSING BOX LIMITS FIELDS")
        log_error("=" * 70)
        log_error(f"File: {resolved}")
        log_error(f"Missing fields: {missing}")
        log_error(f"Found fields: {list(dfbox.index)}")
        log_error("\n📝 Required format:")
        for key in required:
            log_error(f"   {key};<value>")
        log_error("=" * 70 + "\n")
        raise ValueError(
            f"Box limits file missing required fields: {missing}. "
            f"Found: {list(dfbox.index)}"
        )

    limits = tuple(float(dfbox.loc[key].iloc[0]) for key in required)
    return (resolved,) + limits
