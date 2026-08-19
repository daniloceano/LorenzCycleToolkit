# **************************************************************************** #
#                                                                              #
#                                                         :::      ::::::::    #
#    calc_budget_and_residual.py                        :+:      :+:    :+:    #
#                                                     +:+ +:+         +:+      #
#    By: daniloceano <danilo.oceano@gmail.com>      +#+  +:+       +#+         #
#                                                 +#+#+#+#+#+   +#+            #
#    Created: 2022/01/31 20:15:59 by daniloceano       #+#    #+#              #
#    Updated: 2024/01/03 22:33:44 by daniloceano      ###   ########.fr        #
#                                                                              #
# **************************************************************************** #

"""
Budget Difference and Residual Calculations for the Lorenz Energy Cycle.

Created by:
    Danilo Couto de Souza
    Universidade de São Paulo (USP)
    Instituto de Astornomia, Ciências Atmosféricas e Geociências
    São Paulo - Brazil

Contact:
    danilo.oceano@gmail.com
"""

import logging

import numpy as np
import pandas as pd


def elapsed_seconds(dates: np.ndarray, app_logger: logging.Logger = None) -> np.ndarray:
    """
    Convert an array of datetimes into elapsed seconds relative to the first
    sample, validating that the series is strictly increasing.

    Args:
        dates (Array-like): datetime64 values (or anything pandas can parse).
        app_logger (logging.Logger): optional logger.

    Returns:
        np.ndarray: elapsed time in seconds, float.

    Raises:
        ValueError: if fewer than two samples are supplied, or if the
            timestamps are not strictly increasing.
    """
    times = pd.to_datetime(np.asarray(dates))
    values = np.asarray(times, dtype="datetime64[ns]")

    if values.size < 2:
        raise ValueError(
            "At least two time steps are required to compute a tendency; "
            f"got {values.size}."
        )

    seconds = (values - values[0]) / np.timedelta64(1, "s")
    seconds = np.asarray(seconds, dtype=float)

    steps = np.diff(seconds)
    if np.any(steps <= 0):
        bad = int(np.argmin(steps))
        raise ValueError(
            "Time coordinate must be strictly increasing to compute tendencies; "
            f"found a non-positive step of {steps[bad]:.1f} s between "
            f"index {bad} ({values[bad]}) and index {bad + 1} ({values[bad + 1]})."
        )

    if app_logger is not None and steps.size:
        unique = np.unique(np.round(steps, 6))
        if unique.size > 1:
            app_logger.info(
                f"⏱️ Irregular time sampling detected: {unique.size} distinct "
                f"step lengths (min {steps.min():.1f} s, max "
                f"{steps.max():.1f} s). Tendencies use the actual time "
                "coordinate."
            )
        else:
            app_logger.debug(
                f"Regular time sampling: constant step of {steps[0]:.1f} s."
            )

    return seconds


def calc_budget_diff(df: pd.DataFrame, dates: np.ndarray, app_logger: logging.Logger):
    """
    Estimate budget values for energy terms using finite differences.

    The derivative is taken with respect to the ACTUAL time coordinate rather
    than a single leading interval, so that irregular or incomplete sampling is
    handled correctly.  For a perfectly regular series the result is identical
    to the previous behaviour to within floating-point precision.

    Args:
        df (DataFrame): DataFrame containing energy terms.
        dates (Array-like): Array of datetime objects representing time points.

    Returns:
        DataFrame: Updated DataFrame with budget values.
    """
    app_logger.debug("Estimating budget values using finite differences...")

    seconds = elapsed_seconds(dates, app_logger)
    energy_terms = ["Az", "Ae", "Kz", "Ke"]

    try:
        for term in energy_terms:
            values = np.asarray(df[term], dtype=float)
            if values.size != seconds.size:
                raise ValueError(
                    f"Length mismatch for '{term}': {values.size} values for "
                    f"{seconds.size} time steps."
                )
            df[f"∂{term}/∂t (finite diff.)"] = np.gradient(values, seconds)
    except Exception as e:
        app_logger.error(f"Error in calc_budget_diff: {e}")
        raise

    app_logger.debug("Done.")
    return df


def calc_budget_diff_4th(
    df: pd.DataFrame, time: np.ndarray, app_logger: logging.Logger
):
    """
    Estimate budget values for energy terms using 4th order finite differences.

    Args:
        df (DataFrame): DataFrame containing energy terms.
        time (Array-like): Array of datetime objects representing time points.

    Returns:
        DataFrame: Updated DataFrame with budget values.
    """
    app_logger.debug("Estimating budget values using 4th order finite differences...")
    seconds = elapsed_seconds(time, app_logger)
    steps = np.diff(seconds)
    if np.ptp(steps) > 1e-6:
        raise ValueError(
            "calc_budget_diff_4th requires a uniformly sampled time series; "
            f"step lengths range from {steps.min():.1f} s to {steps.max():.1f} s. "
            "Use calc_budget_diff instead."
        )
    dt = float(steps[0])
    energy_terms = ["Az", "Ae", "Kz", "Ke"]

    try:
        for term in energy_terms:
            df = _apply_4th_order_diff(df, term, dt)
    except Exception as e:
        app_logger.error(f"Error in calc_budget_diff_4th: {e}")
        raise

    app_logger.debug("Done.")
    return df


def _apply_4th_order_diff(df, term, dt):
    """
    Apply 4th order finite difference to a specific energy term.

    Args:
        df (DataFrame): DataFrame containing energy terms.
        term (str): The energy term to apply the finite difference.
        dt (float): Time difference in seconds.

    Returns:
        DataFrame: Updated DataFrame with the specific term calculated.
    """
    forward = (df[term].iloc[1] - df[term].iloc[0]) / dt
    central_second = (df[term].iloc[2] - df[term].iloc[0]) / (2 * dt)
    central_penultimate = (df[term].iloc[-1] - df[term].iloc[-3]) / (2 * dt)
    fourth_order = _compute_4th_order_terms(df[term], dt)
    backward = (df[term].iloc[-1] - df[term].iloc[-2]) / dt

    df[f"∂{term}/∂t (finite diff.)"] = (
        [forward, central_second] + fourth_order + [central_penultimate, backward]
    )
    return df


def _compute_4th_order_terms(series, dt):
    """
    Compute 4th order finite difference terms.

    Args:
        series (Series): Pandas Series containing the specific energy term.
        dt (float): Time difference in seconds.

    Returns:
        list: List of calculated 4th order terms.
    """
    fourth_order1 = (
        (4 / 3) * (series.iloc[3:-1].values - series.iloc[1:-3].values) / (2 * dt)
    )
    fourth_order2 = (
        (1 / 3) * (series.iloc[4:].values - series.iloc[:-4].values) / (4 * dt)
    )
    return list(fourth_order1 - fourth_order2)


def calc_residuals(df: pd.DataFrame, app_logger: logging.Logger):
    """
    Compute the residuals RGz, RKz, RGe, and RKe using estimated budget terms.

    Args:
        df (DataFrame): DataFrame containing budget and conversion terms.

    Returns:
        DataFrame: Updated DataFrame with residuals.
    """
    app_logger.debug("Estimating residuals...")

    try:
        df["RGz"] = df["∂Az/∂t (finite diff.)"] + df["Cz"] + df["Ca"] - df["BAz"]
        df["RKz"] = df["∂Kz/∂t (finite diff.)"] - df["Cz"] - df["Ck"] - df["BKz"]
        df["RGe"] = df["∂Ae/∂t (finite diff.)"] - df["Ca"] + df["Ce"] - df["BAe"]
        df["RKe"] = df["∂Ke/∂t (finite diff.)"] - df["Ce"] + df["Ck"] - df["BKe"]

    except Exception as e:
        app_logger.error(f"Error in calc_residuals: {e}")
        raise

    app_logger.debug("Done.")
    return df
