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


def calc_budget_diff(df: pd.DataFrame, dates: np.ndarray, app_logger: logging.Logger):
    """
    Estimate budget values for energy terms using finite differences method.

    Args:
        df (DataFrame): DataFrame containing energy terms.
        dates (Array-like): Array of datetime objects representing time points.

    Returns:
        DataFrame: Updated DataFrame with budget values.
    """
    app_logger.debug("Estimating budget values using finite differences...")

    dt = float((dates[1] - dates[0]) / np.timedelta64(1, "s"))
    energy_terms = ["Az", "Ae", "Kz", "Ke"]

    try:
        for term in energy_terms:
            df[f"∂{term}/∂t (finite diff.)"] = np.gradient(df[term], dt)
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
    dt = float((time[1] - time[0]) / np.timedelta64(1, "s"))
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
