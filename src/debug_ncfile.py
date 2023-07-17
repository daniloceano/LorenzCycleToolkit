# **************************************************************************** #
#                                                                              #
#                                                         :::      ::::::::    #
#    debug_ncfile.py                                    :+:      :+:    :+:    #
#                                                     +:+ +:+         +:+      #
#    By: Danilo <danilo.oceano@gmail.com>           +#+  +:+       +#+         #
#                                                 +#+#+#+#+#+   +#+            #
#    Created: 2023/07/17 14:39:44 by Danilo            #+#    #+#              #
#    Updated: 2023/07/17 17:36:54 by Danilo           ###   ########.fr        #
#                                                                              #
# **************************************************************************** #

import os
import argparse
import pandas as pd
import xarray as xr
import numpy as np
import dask
from metpy.units import units
from thermodynamics import StaticStability
import matplotlib.pyplot as plt
from Math import CalcAreaAverage, CalcZonalAverage
from metpy.constants import Re


def convert_longitude(data, lon_indexer):
    data.coords[lon_indexer] = (data.coords[lon_indexer] + 180) % 360 - 180
    data = data.sortby(data[lon_indexer])
    return data


def get_data(infile: str, varlist: str) -> xr.Dataset:
    """
    Read input data from a NetCDF file and perform necessary preprocessing.

    Args:
        infile (str): Path to the input NetCDF file.
        varlist (str): Path to the variable list file.

    Returns:
        xr.Dataset: Preprocessed input data.

    Raises:
        FileNotFoundError: If the variable list file or input NetCDF file cannot be found.
        pd.errors.EmptyDataError: If the variable list file is empty.
        Exception: If an error occurs while opening the input NetCDF file.
    """
    print(f"Variables specified by the user in: {varlist}")
    print(f"Attempting to read {varlist} file...")
    try:
        df_vars = pd.read_csv(varlist, sep=';', index_col=0, header=0)
    except FileNotFoundError:
        raise SystemExit("Error: The 'fvar' text file could not be found.")
    except pd.errors.EmptyDataError:
        raise SystemExit("Error: The 'fvar' text file is empty.")

    print("List of variables found:")
    print(df_vars)

    lon_indexer = df_vars.loc["Longitude"]["Variable"]
    lat_indexer = df_vars.loc["Latitude"]["Variable"]
    level_indexer = df_vars.loc["Vertical Level"]["Variable"]

    print("Opening input data...")
    try:
        with dask.config.set(array={"slicing": {"split_large_chunks": True}}):
            data = convert_longitude(
                xr.open_dataset(infile),
                df_vars.loc["Longitude"]["Variable"]
            )
    except FileNotFoundError:
        raise SystemExit("ERROR: Could not open file. Check if path, fvars file, and file format (.nc) are correct.")
    except Exception as e:
        raise SystemExit(f"ERROR: An error occurred while opening the file: {e}")

    print("Assigning geospatial coordinates in radians...")
    data = data.assign_coords({"rlats": np.deg2rad(data[lat_indexer])})
    data = data.assign_coords({"coslats": np.cos(np.deg2rad(data[lat_indexer]))})
    data = data.assign_coords({"rlons": np.deg2rad(data[lon_indexer])})

    levels_pa = (data[level_indexer] * units(str(data[level_indexer].units))).metpy.convert_units("Pa")
    data = data.assign_coords({level_indexer: levels_pa})

    data = data.sortby(lon_indexer).sortby(level_indexer, ascending=True).sortby(lat_indexer, ascending=True)

    print("Data opened successfully.")
    return data


def calc_ca(tair_AE, tair_ZE, v_ZE, omega_ZE, sigma_AA, PressureData, VerticalCoordIndexer, xlength, ylength):
    """
    Calculate the conversion between available potential energy terms (Ca).

    Args:
        tair_AE (xr.DataArray): Eddy component of air temperature.
        tair_ZE (xr.DataArray): Zonal average of air temperature.
        v_ZE (xr.DataArray): Zonal average of northward wind component.
        omega_ZE (xr.DataArray): Zonal average of omega velocity.
        sigma_AA (xr.DataArray): Static stability parameter.
        PressureData (xr.DataArray): Pressure data.
        VerticalCoordIndexer (str): Name of the vertical coordinate indexer.
        xlength (float): Length of the domain in the x-direction.
        ylength (xr.DataArray): Length of the domain in the y-direction.

    Returns:
        xr.DataArray: Conversion between available potential energy terms (Ca).
        xr.DataArray: First term of the integral.
        xr.DataArray: Second term of the integral.
        xr.DataArray: Derivative of air temperature with respect to pressure.
        xr.DataArray: Omega velocity multiplied by air temperature and derivative of air temperature with respect to pressure.
    """
    # First term of the integral
    DelPhi_tairAE = (tair_AE * tair_AE["coslats"]).differentiate("rlats")
    tair_v_ZE_sigma_AA = (v_ZE * tair_ZE) / (2 * Re * sigma_AA)
    term1 = tair_v_ZE_sigma_AA * DelPhi_tairAE
    term1 = CalcAreaAverage(term1, ylength, xlength=xlength)

    # Second term of the integral
    DelPres_tairAE = (tair_AE).differentiate(VerticalCoordIndexer) / PressureData.metpy.units
    omega_tair_ZE = omega_ZE * tair_ZE
    term2 = omega_tair_ZE * DelPres_tairAE
    term2 = CalcAreaAverage(term2, ylength, xlength=xlength) / sigma_AA

    function = term1 + term2

    if np.isnan(function).any():
        function = function.interpolate_na(dim=VerticalCoordIndexer) * function.metpy.units
        if np.isnan(function).any():
            function = function.dropna(dim=VerticalCoordIndexer)

    Ca = -function.integrate(VerticalCoordIndexer) * PressureData.metpy.units

    return Ca, term1, term2, DelPres_tairAE, omega_tair_ZE


def get_ca_arguments(data, df_vars, western_limit, eastern_limit, southern_limit, northern_limit, Re):
    """
    Extract necessary arguments for calculating Ca.

    Args:
        data (xr.Dataset): Input data.
        df_vars (pd.DataFrame): Variable information dataframe.
        western_limit (float): Western limit of the domain.
        eastern_limit (float): Eastern limit of the domain.
        southern_limit (float): Southern limit of the domain.
        northern_limit (float): Northern limit of the domain.
        Re (float): Earth's radius.

    Returns:
        tuple: Arguments for calculating Ca.
    """
    LonIndexer = df_vars.loc["Longitude"]["Variable"]
    LatIndexer = df_vars.loc["Latitude"]["Variable"]
    VerticalCoordIndexer = df_vars.loc["Vertical Level"]["Variable"]
    PressureData = data[VerticalCoordIndexer] * units("Pa")

    xlength = data[LonIndexer].sel({LonIndexer: eastern_limit}, method="nearest") - data[LonIndexer].sel(
        {LonIndexer: western_limit}, method="nearest"
    )
    ylength = np.sin(
        data[LatIndexer].sel({LatIndexer: northern_limit}, method="nearest")["rlats"]
    ) - np.sin(data[LatIndexer].sel({LatIndexer: southern_limit}, method="nearest")["rlats"])

    tair = (
        data[df_vars.loc["Air Temperature"]["Variable"]]
        * units(df_vars.loc["Air Temperature"]["Units"]).to("K")
    ).sel(
        **{
            LatIndexer: slice(southern_limit, northern_limit),
            LonIndexer: slice(western_limit, eastern_limit),
        }
    )
    tair_ZA = CalcZonalAverage(tair, xlength)
    tair_AA = CalcAreaAverage(tair_ZA, ylength)
    tair_ZE = tair - tair_ZA
    tair_AE = tair_ZA - tair_AA

    v = (
        data[df_vars.loc["Northward Wind Component"]["Variable"]]
        * units(df_vars.loc["Northward Wind Component"]["Units"]).to("m/s")
    ).sel(
        **{
            LatIndexer: slice(southern_limit, northern_limit),
            LonIndexer: slice(western_limit, eastern_limit),
        }
    )
    v_ZA = CalcZonalAverage(v, xlength)
    v_AA = CalcAreaAverage(v_ZA, ylength)
    v_ZE = v - v_ZA

    omega = (
        data[df_vars.loc["Omega Velocity"]["Variable"]]
        * units(df_vars.loc["Omega Velocity"]["Units"]).to("Pa/s")
    ).sel(
        **{
            LatIndexer: slice(southern_limit, northern_limit),
            LonIndexer: slice(western_limit, eastern_limit),
        }
    )
    omega_ZA = CalcZonalAverage(omega, xlength)
    omega_ZE = omega - omega_ZA

    sigma_AA = StaticStability(
        tair, PressureData, VerticalCoordIndexer, xlength, ylength
    )

    return (
        tair_AE,
        tair_ZE,
        v_ZE,
        omega_ZE,
        sigma_AA,
        PressureData,
        VerticalCoordIndexer,
        xlength,
        ylength,
    )


def plot_timeseries(times, values, title, y_label):
    flattened_values = np.array(values).flatten()
    plt.plot(times, flattened_values)
    plt.xlabel("Time")
    plt.ylabel(y_label)
    plt.title(title)
    plt.savefig(f"debug/{y_label}.png")


def main(args):
    infile = args.infile
    varlist = "../inputs/fvars_ERA5-cdsapi"
    trackfile = "/p1-nemo/danilocs/SWSA-cyclones_energetic-analysis//LEC_results-q0.999/RG3-q0.999-20070557_ERA5_track-15x15/RG3-q0.999-20070557_ERA5_track-15x15_track"

    # Indexers
    df_vars = pd.read_csv(varlist, sep=";", index_col=0, header=0)
    lon_indexer = df_vars.loc["Longitude"]["Variable"]
    lat_indexer = df_vars.loc["Latitude"]["Variable"]
    time_name = df_vars.loc["Time"]["Variable"]
    vertical_coord_indexer = df_vars.loc["Vertical Level"]["Variable"]

    # Open data
    data = get_data(infile, varlist)
    track = pd.read_csv(trackfile, parse_dates=[0], delimiter=";", index_col="time")
    times = pd.to_datetime(track.index)

    dTdt = data[df_vars.loc["Air Temperature"]["Variable"]].differentiate(
        df_vars.loc["Time"]["Variable"], datetime_unit="s"
    ) * units("K/s")

    sigma_series = []
    ca_series = []
    DelPhi_tairAE_series = []
    tair_v_ZE_sigma_AA_series = []
    DelPres_tairAE_series = []
    omega_tair_ZE_series = []
    tair_series = []

    for time in times:
        print(f"Computing for {time}")
        itrack = track.loc[time]
        width, length = itrack["width"], itrack["length"]
        western_limit, eastern_limit = itrack["Lon"] - width/2, itrack["Lon"] + width/2
        southern_limit, northern_limit = itrack["Lat"] - length/2, itrack["Lat"] + length/2

        idata = data.sel(
            {time_name: time}
        ).sel(
            {lat_indexer: slice(southern_limit, northern_limit)}
        ).sel(
            {lon_indexer: slice(western_limit, eastern_limit)}
        )

        ca_arguments = get_ca_arguments(
            idata,
            df_vars,
            western_limit,
            eastern_limit,
            southern_limit,
            northern_limit,
            Re
        )

        ylength = northern_limit - southern_limit
        xlength = eastern_limit - western_limit

        Ca, DelPhi_tairAE, tair_v_ZE_sigma_AA, DelPres_tairAE, omega_tair_ZE = calc_ca(*ca_arguments)

        ca_series.append(float(Ca))
        sigma_series.append(float(ca_arguments[4].integrate(coord=vertical_coord_indexer)))

        DelPhi_tairAE_series.append(float(DelPhi_tairAE.integrate(coord=vertical_coord_indexer)))

        tair_v_ZE_sigma_AA_series.append(float(tair_v_ZE_sigma_AA.integrate(coord=vertical_coord_indexer)))

        DelPres_tairAE_AA = CalcAreaAverage(DelPres_tairAE, ylength=ylength)
        DelPres_tairAE_series.append(float(DelPres_tairAE_AA.integrate(coord=vertical_coord_indexer)))

        omega_tair_ZE_AA = CalcAreaAverage(omega_tair_ZE, ylength=ylength, xlength=xlength)
        omega_tair_ZE_series.append(float(omega_tair_ZE_AA.integrate(coord=vertical_coord_indexer)))

        tair_AA = CalcAreaAverage(
            idata[df_vars.loc["Air Temperature"]["Variable"]], ylength=ylength, xlength=xlength)
        tair_series.append(float(tair_AA.integrate(coord=vertical_coord_indexer)))

    os.makedirs("debug", exist_ok=True)
    plot_timeseries(times, ca_series, "Timeseries of Ca", "Ca")
    plot_timeseries(times, DelPhi_tairAE_series, "Timeseries of DelPhi_tairAE", "DelPhi_tairAE")
    plot_timeseries(times, tair_v_ZE_sigma_AA_series, "Timeseries of tair_v_ZE_sigma_AA", "tair_v_ZE_sigma_AA")
    plot_timeseries(times, DelPres_tairAE_series, "Timeseries of DelPres_tairAE", "DelPres_tairAE")
    plot_timeseries(times, omega_tair_ZE_series, "Timeseries of omega_tair_ZE", "omega_tair_ZE")
    plot_timeseries(times, sigma_series, "Timeseries of Static Stability", "Static Stability")
    plot_timeseries(times, tair_series, "Timeseries of Temperature", "Temperature")





if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Lorenz Energy Cycle (LEC) program.")
    parser.add_argument("infile", help="Input .nc file with temperature, geopotential, and wind components.")
    parser.add_argument("-r", "--residuals", action="store_true", help="Compute the Dissipation and Generation terms as residuals.")
    parser.add_argument("-g", "--geopotential", action="store_true", help="Use geopotential data instead of geopotential height.")
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("-f", "--fixed", action="store_true", help="Compute the energetics for a fixed domain specified by the 'box_lims' file.")
    group.add_argument("-t", "--track", action="store_true", help="Define the box using a track file specified by the 'track' file.")
    group.add_argument("-c", "--choose", action="store_true", help="Choose the domain for each time step by clicking on the screen.")
    parser.add_argument("-z", "--zeta", action="store_true", help="Use this flag if the track file was created using vorticity.")
    parser.add_argument("-m", "--mpas", action="store_true", help="for MPAS-A data processed with MPAS-BR routines")
    parser.add_argument("-o", "--outname", type=str, help="Choose a name for saving results.")
    parser.add_argument("-v", "--verbosity", action="store_true", help="Increase output verbosity.")

    args = parser.parse_args([
        "/p1-nemo/danilocs/SWSA-cyclones_energetic-analysis/met_data/ERA5/scripts/APIs/RG3-q0.99-20070557_ERA5.nc",
        "-r", "-t", "-g", "-z"
    ])

    main(args)
