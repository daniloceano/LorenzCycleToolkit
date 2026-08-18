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


def read_box_limits(path: str):
    """
    Read a box-limits file, resolving the ``.default`` fallback first.

    Args:
        path (str): Path to the box-limits file, e.g. ``inputs/box_limits``.

    Returns:
        tuple: ``(resolved_path, min_lon, max_lon, min_lat, max_lat)``, with the
        four limits as raw floats in the file's own longitude convention. They
        are returned unvalidated: :func:`src.utils.longitude.canonicalize_box`
        is the single place where the longitude pair is checked and expressed
        in the dataset's convention.

    Raises:
        KeyError: If any of the four rows is missing from the file.
    """
    resolved = resolve_input_file(path)
    dfbox = pd.read_csv(resolved, header=None, delimiter=";", index_col=0)
    limits = tuple(
        float(dfbox.loc[key].iloc[0])
        for key in ("min_lon", "max_lon", "min_lat", "max_lat")
    )
    return (resolved,) + limits
