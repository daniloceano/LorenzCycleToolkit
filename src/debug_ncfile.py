# **************************************************************************** #
#                                                                              #
#                                                         :::      ::::::::    #
#    debug_ncfile.py                                    :+:      :+:    :+:    #
#                                                     +:+ +:+         +:+      #
#    By: Danilo  <danilo.oceano@gmail.com>          +#+  +:+       +#+         #
#                                                 +#+#+#+#+#+   +#+            #
#    Created: 2023/07/17 14:39:44 by Danilo            #+#    #+#              #
#    Updated: 2023/07/19 10:46:52 by Danilo           ###   ########.fr        #
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
import cartopy.crs as ccrs
import matplotlib.pyplot as plt
from cartopy.mpl.gridliner import LATITUDE_FORMATTER, LONGITUDE_FORMATTER

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
        print('interpolating')
        function = function.interpolate_na(dim=VerticalCoordIndexer) * function.metpy.units
        if np.isnan(function).any():
            print('removing nans')
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


def plot_timeseries(times, values, y_label):
    plt.close("all")
    plt.figure()
    flattened_values = np.array(values).flatten()
    plt.plot(times, flattened_values)
    plt.ylabel(y_label)
    plt.xticks(rotation=45, ha="right")  # Adjust x ticks rotation and alignment
    plt.tight_layout()  # Ensure tight layout to prevent overlapping
    plt.savefig(f"debug/{y_label}.png")

def analyse_timeseries(data, varlist, times, track, slice_flag=False):

    # Indexers
    df_vars = pd.read_csv(varlist, sep=";", index_col=0, header=0)
    lon_indexer = df_vars.loc["Longitude"]["Variable"]
    lat_indexer = df_vars.loc["Latitude"]["Variable"]
    time_name = df_vars.loc["Time"]["Variable"]
    vertical_coord_indexer = df_vars.loc["Vertical Level"]["Variable"]

    sigma_series = []
    ca_series = []
    DelPhi_tairAE_series = []
    tair_v_ZE_sigma_AA_series = []
    DelPres_tairAE_series = []
    omega_tair_ZE_series = []
    tair_series = []
    tair_ZE_series = []

    # Specify the desired dates
    start_date = pd.Timestamp("2007-09-07")
    end_date = pd.Timestamp("2007-09-11")

    # Slice the times for the specified dates
    times = times[(times >= start_date) & (times <= end_date)]

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

        tair = idata[df_vars.loc["Air Temperature"]["Variable"]]
        tair_AA = CalcAreaAverage(tair, ylength=ylength, xlength=xlength)
        tair_series.append(float(tair_AA.integrate(coord=vertical_coord_indexer)))

        tair_ZA = CalcZonalAverage(tair, xlength)
        tair_AA = CalcAreaAverage(tair_ZA, ylength)
        tair_ZE = tair - tair_ZA

        tair_ZE_AA = CalcAreaAverage(tair_ZE, ylength=ylength, xlength=xlength)
        tair_ZE_series.append(float(tair_ZE_AA.integrate(coord=vertical_coord_indexer)))

    os.makedirs("debug", exist_ok=True)

    if slice_flag == False:
        plot_timeseries(times, ca_series, "Ca")
        plot_timeseries(times, DelPhi_tairAE_series, "DelPhi_tairAE")
        plot_timeseries(times, tair_v_ZE_sigma_AA_series, "tair_v_ZE_sigma_AA")
        plot_timeseries(times, DelPres_tairAE_series, "DelPres_tairAE")
        plot_timeseries(times, omega_tair_ZE_series, "omega_tair_ZE")
        plot_timeseries(times, sigma_series, "Static Stability")
        plot_timeseries(times, tair_series, "Temperature")
        plot_timeseries(times, tair_ZE_series, "tair_ZE")
    else:
        plot_timeseries(times, ca_series, "Ca")
        plot_timeseries(times, DelPhi_tairAE_series, "DelPhi_tairAE_sliced")
        plot_timeseries(times, tair_v_ZE_sigma_AA_series, "tair_v_ZE_sigma_AA_sliced")
        plot_timeseries(times, tair_ZE_series, "tair_ZE_sliced")        

def analyse_tair(data, time, track, varlist):
    # Indexers
    df_vars = pd.read_csv(varlist, sep=";", index_col=0, header=0)
    lon_indexer = df_vars.loc["Longitude"]["Variable"]
    lat_indexer = df_vars.loc["Latitude"]["Variable"]
    time_name = df_vars.loc["Time"]["Variable"]
    vertical_coord_indexer = df_vars.loc["Vertical Level"]["Variable"]

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
    tair = idata[df_vars.loc["Air Temperature"]["Variable"]]

    # Create a meshgrid for lon, lat coordinates
    lon = idata[lon_indexer]
    lat = idata[lat_indexer]

    maps_dir = 'debug/maps'
    os.makedirs(maps_dir, exist_ok=True)

    # Loop through each vertical level
    for level in tair[vertical_coord_indexer]:

        # Plotting
        plt.close("all")
        plt.figure(figsize=(10, 6))
        ax = plt.axes(projection=ccrs.PlateCarree())
        
        level = float(level)
        level_data = tair.sel({vertical_coord_indexer: level})
        level_hPa = level/100

        # Plot contourf of tair and get the mappable object
        contour = ax.contourf(lon, lat, level_data, cmap='coolwarm', transform=ccrs.PlateCarree())
        plt.colorbar(contour, label='Temperature')

        ax.coastlines(resolution='50m', color='black',
                       linewidth=1, alpha=0.5, transform=ccrs.PlateCarree(), zorder=101)

        # Add gridlines with labels for left and bottom sides
        gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True, linewidth=1, color='gray', alpha=0.5)
        gl.top_labels = False
        gl.right_labels = False
        gl.xformatter = LONGITUDE_FORMATTER
        gl.yformatter = LATITUDE_FORMATTER

        # Save the figure for each vertical level
        plt.savefig(f'{maps_dir}/tair_level_{level_hPa}.png')

        # Clear the plot for the next level
        plt.clf()

def analyse_tair_AE(data, time, track, varlist, slice_flag=False):
    # Indexers
    df_vars = pd.read_csv(varlist, sep=";", index_col=0, header=0)
    lon_indexer = df_vars.loc["Longitude"]["Variable"]
    lat_indexer = df_vars.loc["Latitude"]["Variable"]
    time_name = df_vars.loc["Time"]["Variable"]
    vertical_coord_indexer = df_vars.loc["Vertical Level"]["Variable"]

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

    ylength = northern_limit - southern_limit
    xlength = eastern_limit - western_limit

    tair = idata[df_vars.loc["Air Temperature"]["Variable"]]
    tair_ZA = CalcZonalAverage(tair, xlength)
    tair_AA = CalcAreaAverage(tair_ZA, ylength)
    tair_AE = tair_ZA - tair_AA
    DelPres_tairAE = (tair_AE).differentiate(vertical_coord_indexer) / units('Pa')
    DelPhi_tairAE = (tair_AE * tair_AE["coslats"]).differentiate("rlats")

    tair_AE_AA = CalcAreaAverage(tair_AE, ylength)
    DelPres_tairAE_AA = CalcAreaAverage(DelPres_tairAE, ylength)
    DelPhi_tairAE_AA = CalcAreaAverage(DelPhi_tairAE, ylength)

    print(float(tair_AE_AA.integrate(coord=vertical_coord_indexer)))
    print(float(DelPres_tairAE_AA.integrate(coord=vertical_coord_indexer)))
    print(float(DelPhi_tairAE_AA.integrate(coord=vertical_coord_indexer)))

    if slice_flag == False:
        # Plot tair_AE
        plot_panel(tair_AE, lat_indexer, "debug/tair_AE")
        plot_timeseries(tair_AE.latitude, tair_AE.isel(level=0), f"tair_AE_{float(tair_AE.level[0])}hPa")
        plot_timeseries(DelPres_tairAE_AA.level, DelPres_tairAE_AA, "DelPres_tairAE_AA")
        print(f"plotting tair_AE for: {float(tair_AE.level[0])}")

        # Plot DelPres_tairAE
        plot_panel(DelPres_tairAE, lat_indexer, "debug/DelPres_tairAE")
        plot_timeseries(DelPhi_tairAE_AA.level, DelPhi_tairAE_AA, "DelPhi_tairAE_AA")

        # Plot DelPres_tairAE
        plot_panel(DelPhi_tairAE, lat_indexer, "debug/DelPhi_tairAE")
        plot_timeseries(DelPhi_tairAE_AA.level, DelPhi_tairAE_AA, "DelPhi_tairAE_AA")

    else:
        plot_timeseries(tair_AE.latitude, tair_AE.isel(level=0), f"tair_AE_{float(tair_AE.level[0])}Pa")
        plot_panel(tair_AE, lat_indexer, "debug/tair_AE_slice")
        plot_timeseries(DelPres_tairAE_AA.level, DelPres_tairAE_AA, "DelPres_tairAE_AA_slice")

        # Plot DelPres_tairAE
        plot_panel(DelPres_tairAE, lat_indexer, "debug/DelPres_tairAE_slice")
        plot_timeseries(DelPhi_tairAE_AA.level, DelPhi_tairAE_AA, "DelPhi_tairAE_AA_slice")

        # Plot DelPres_tairAE
        plot_panel(DelPhi_tairAE, lat_indexer, "debug/DelPhi_tairAE_slice")
        plot_timeseries(DelPhi_tairAE_AA.level, DelPhi_tairAE_AA, "DelPhi_tairAE_AA_slice")

def plot_panel(data, lat_indexer, title):
    # Get the number of vertical levels
    num_levels = len(data['level'])

    # Determine the number of rows and columns for subplots
    num_rows = int(num_levels / 4) + (num_levels % 4 > 0)
    num_cols = min(num_levels, 4)

    # Create the panel of subplots
    plt.close("all")
    fig, axes = plt.subplots(num_rows, num_cols, figsize=(12, 12))

    # Flatten the axes array in case there is only one row or one column
    axes = axes.flatten()

    # Loop through each vertical level and plot the data in the corresponding subplot
    for i, level in enumerate(data['level']):
        level_hPa = level/100
        ax = axes[i]
        ax.plot(data[lat_indexer], data.sel(level=level))
        ax.set_title(f"Level: {level_hPa:.2f} hPa")
        ax.set_aspect('auto')
        ax.grid(True)

    # Adjust the layout of subplots
    plt.tight_layout()

    # Save the figure
    plt.savefig(f"{title}_panel.png")

def main(args):
    infile = args.infile
    varlist = "../inputs/fvars_ERA5-cdsapi"
    trackfile = "/p1-nemo/danilocs/SWSA-cyclones_energetic-analysis//LEC_results-q0.999/RG3-q0.999-20070557_ERA5_track-15x15/RG3-q0.999-20070557_ERA5_track-15x15_track"

    # Indexers
    df_vars = pd.read_csv(varlist, sep=";", index_col=0, header=0)
    vertical_coord_indexer = df_vars.loc["Vertical Level"]["Variable"]
    
    # Open data
    data = get_data(infile, varlist)
    track = pd.read_csv(trackfile, parse_dates=[0], delimiter=";", index_col="time")
    times = pd.to_datetime(track.index)

    analyse_timeseries(data, varlist, times, track)

    # Slice the data for levels from 100000 to 100
    sliced_data = data.sel({vertical_coord_indexer: slice(1000, 100000)})
    analyse_timeseries(sliced_data, varlist, times, track, slice_flag=True)

    time = pd.Timestamp("2007-09-09 00:00")
    analyse_tair(data, time, track, varlist)

    analyse_tair_AE(data, time, track, varlist)
    analyse_tair_AE(sliced_data, time, track, varlist, slice_flag=True)


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
