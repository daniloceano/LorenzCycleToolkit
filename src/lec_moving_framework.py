# **************************************************************************** #
#                                                                              #
#                                                         :::      ::::::::    #
#    moving_framework.py                                :+:      :+:    :+:    #
#                                                     +:+ +:+         +:+      #
#    By: daniloceano <danilo.oceano@gmail.com>      +#+  +:+       +#+         #
#                                                 +#+#+#+#+#+   +#+            #
#    Created: 2023/12/19 17:32:55 by daniloceano       #+#    #+#              #
#    Updated: 2023/12/19 17:33:05 by daniloceano      ###   ########.fr        #
#                                                                              #
# **************************************************************************** #

import os
import pandas as pd

from EnergyContents import EnergyContents
from ConversionTerms import ConversionTerms
from BoundaryTerms import BoundaryTerms
from GenerationDissipationTerms import GenerationDissipationTerms
from BoxData import BoxData
from BudgetResidual import calc_budget_diff,calc_residuals

import logging

def LEC_moving(data, dfVars, dTdt, ResultsSubDirectory, FigsDirectory):
    """
    Computes the Lorenz Energy Cycle using a moving framework.
    
    Args:
    data: A Xarray Dataset containing the data to compute the energy cycle.
    
    Returns:
    None
    """
    print('Computing energetics using moving framework')

    # Indexers
    LonIndexer, LatIndexer, TimeName, VerticalCoordIndexer = (
      dfVars.loc['Longitude']['Variable'],
      dfVars.loc['Latitude']['Variable'],
      dfVars.loc['Time']['Variable'],
      dfVars.loc['Vertical Level']['Variable']
      )
    
    PressureData = data[VerticalCoordIndexer]
    
    # Create csv files for storing vertical results
    for term in ['Az', 'Ae', 'Kz', 'Ke', 'Cz', 'Ca', 'Ck', 'Ce', 'Ge', 'Gz']:
        columns = [TimeName] + [float(i) for i in PressureData.values]
        df = pd.DataFrame(columns=columns)
        file_name = term + '_' + VerticalCoordIndexer + '.csv'
        file_path = os.path.join(ResultsSubDirectory, file_name)
        df.to_csv(file_path, index=None) 
        print(file_path+' created (but still empty)')

    # Track file
    if args.track:
        trackfile = '../inputs/track'
        track = pd.read_csv(trackfile,parse_dates=[0],
                            delimiter=';',index_col='time')

    # Dictionary for saving system position and attributes
    results_keys = ['datestr', 'central_lat', 'central_lon', 'length', 'width',
                'min_max_zeta_850', 'min_hgt_850', 'max_wind_850']
    out_track = pd.DataFrame(columns=results_keys)

    # Create dict for store results
    TermsDict = {}
    energy = ['Az','Ae','Kz','Ke']
    conversion = ['Cz','Ca','Ck','Ce']
    boundary = ['BAz','BAe','BKz','BKe','BΦZ','BΦE']
    gendiss = ['Gz','Ge']
    if not args.residuals:  
        gendiss = ['Gz','Ge','Dz','De']
    for term in [*energy,*conversion,*boundary,*gendiss]:
        TermsDict[term] = []
    
    # Slice the time array so the first and the last timestep will be the same
    # as in the track file
    times = pd.to_datetime(data[TimeName].values)
    if args.track:
        times = times[(times>=track.index[0]) & (times<=track.index[-1])]
        if len(times) == 0:
            raise ValueError(
                "Mismatch between trackfile and data! Check that and try again!"
                f" file starts and ends in times {times[0]} and {times[-1]},"
                f"while track starts and ends in times {track.index[0]} and {track.index[-1]}")
        
    # times = times[:10]
    for t in times:

        idata = data.sel({TimeName: t})
        idTdt = dTdt.sel({TimeName: t})

        # Sometimes there are errors on indexing the time coordinate on the NetCDF file
        #  and I can't do anything about it. So I will have to discard one of the timesteps.
        # Any complaints can be sent directly to the original owner of the .nc file.
        if idata[TimeName].shape != ():
            idata = data.isel({TimeName: 1})
            idTdt = dTdt.isel({TimeName: 1})

        iu = (idata[dfVars.loc['Eastward Wind Component']['Variable']].compute() *
            units(dfVars.loc['Eastward Wind Component']['Units']).to('m/s'))
        iv = (idata[dfVars.loc['Northward Wind Component']['Variable']].compute() *
            units(dfVars.loc['Northward Wind Component']['Units']).to('m/s'))
        if args.geopotential:
            ihgt = ((idata[dfVars.loc['Geopotential']['Variable']].compute() *
            units(dfVars.loc['Geopotential']['Units']))/g).metpy.convert_units('gpm')
        else:
            ihgt = (idata[dfVars.loc['Geopotential Height']['Variable']].compute() *
            units(dfVars.loc['Geopotential Height']['Units']).to('gpm'))

        iu_850, iv_850, ight_850 = (
            iu.sel({VerticalCoordIndexer: 85000}),
            iv.sel({VerticalCoordIndexer: 85000}),
            ihgt.sel({VerticalCoordIndexer: 85000})
        )

        iwspd_850, izeta_850 = wind_speed(iu_850, iv_850), vorticity(iu_850, iv_850).metpy.dequantify()

        lat, lon = idata[LatIndexer], idata[LonIndexer]

        # Get current time and box limits
        itime = str(t)
        datestr = pd.to_datetime(itime).strftime('%Y-%m-%d-%H%M')

        if args.track:
            # Get current time and box limits
            if 'width'in track.columns:
                width, length = track.loc[itime]['width'],track.loc[itime]['length']
            else:
                width, length = 15, 15
                
            # Track timestep closest to the model timestep, just in case
            # the track file has a poorer temporal resolution
            closest_index = np.searchsorted(track.index, t)
            if closest_index > 0 and (closest_index == len(track.index
                                ) or abs(t - track.index[closest_index-1]
                                ) < abs(t - track.index[closest_index])):
                closest_index -= 1
            track_itime = track.index[closest_index]
                # Store system position and attributes
            central_lon, central_lat = track.loc[track_itime]['Lon'], track.loc[track_itime]['Lat']
            min_lon = central_lon-(width/2)
            max_lon = central_lon+(width/2)
            min_lat = central_lat-(length/2)
            max_lat = central_lat+(length/2)
            limits = {'min_lon':min_lon,'max_lon':max_lon,
                        'min_lat':min_lat,'max_lat':max_lat}
            
            # Slice data for defined box
            izeta_850_slice = izeta_850.sel({LatIndexer:slice(min_lat, max_lat), LonIndexer:slice(min_lon, max_lon)})
            ight_850_slice = ight_850.sel({LatIndexer:slice(min_lat, max_lat), LonIndexer:slice(min_lon, max_lon)})
            iwspd_850_slice = iwspd_850.sel({LatIndexer:slice(min_lat, max_lat), LonIndexer:slice(min_lon, max_lon)})

            # Check if 'min_max_zeta_850', 'min_hgt_850' and 'max_wind_850' columns exists in the track file.
            # If they exist, then retrieve and convert the value from the track file.
            # If they do not exist, calculate them.
            try:
                min_max_zeta = float(track.loc[track_itime]['min_max_zeta_850'])
            except KeyError:
                if args.zeta:
                    min_max_zeta_unformatted = izeta_850.sel(latitude=central_lat, longitude=central_lon, method='nearest')
                else:
                    min_max_zeta_unformatted = izeta_850_slice.min()
                if min_lat < 0:
                    min_max_zeta = float(np.nanmin(min_max_zeta_unformatted))
                else:
                    min_max_zeta = float(np.nanmax(min_max_zeta_unformatted))

            try:
                min_hgt = float(track.loc[track_itime]['min_hgt_850'])
            except KeyError:
                min_hgt = float(ight_850_slice.min())

            try:
                max_wind = float(track.loc[track_itime]['max_wind_850'])
            except KeyError:
                max_wind = float(iwspd_850_slice.max())

        elif args.choose:
            # Draw maps and ask user to specify corners for specifying the box
            limits = draw_box_map(iu_850, iv_850, izeta_850, ight_850,
                                    lat, lon, itime)
            
            # Store system position and attributes
            min_lon, max_lon = limits['min_lon'],  limits['max_lon']
            min_lat, max_lat = limits['min_lat'],  limits['max_lat']
            width, length = limits['max_lon'] - limits['min_lon'], limits['max_lat'] - limits['min_lat']
            central_lat = (limits['max_lat'] + limits['min_lat'])/2
            central_lon = (limits['max_lon'] + limits['min_lon'])/2

            # Slice data for defined box and find extremes
            izeta_850_slice = izeta_850.sel({LatIndexer:slice(min_lat, max_lat), LonIndexer:slice(min_lon, max_lon)})
            ight_850_slice = ight_850.sel({LatIndexer:slice(min_lat, max_lat), LonIndexer:slice(min_lon, max_lon)})
            iwspd_850_slice = iwspd_850.sel({LatIndexer:slice(min_lat, max_lat), LonIndexer:slice(min_lon, max_lon)})
            min_max_zeta = float(izeta_850_slice.min())
            min_hgt = float(ight_850_slice.min())
            max_wind = float(iwspd_850_slice.max())

        # Find position of the extremes
        lat_slice, lon_slice = izeta_850_slice[LatIndexer], izeta_850_slice[LonIndexer]
        min_max_zeta_lat, min_max_zeta_lon = find_extremum_coordinates(izeta_850_slice, lat_slice, lon_slice, 'min_max_zeta')
        min_hgt_lat, min_hgt_lon = find_extremum_coordinates(ight_850_slice, lat_slice, lon_slice, 'min_hgt')
        max_wind_lat, max_wind_lon = find_extremum_coordinates(iwspd_850_slice, lat_slice, lon_slice, 'max_wind')

        # Store the results in a dictionary for plotting purposes
        data850 = {
            'min_max_zeta': {
                'latitude': min_max_zeta_lat,
                'longitude': min_max_zeta_lon,
                'data': izeta_850
            },
            'min_hgt': {
                'latitude': min_hgt_lat,
                'longitude': min_hgt_lon,
                'data': ight_850
            },
            'max_wind': {
                'latitude': max_wind_lat,
                'longitude': max_wind_lon,
                'data': iwspd_850
            },
            'lat': lat,
            'lon': lon,
        }

        position = {
        'datestr': datestr,
        'central_lat': central_lat,
        'central_lon': central_lon,
        'length': length,
        'width': width,
        'min_max_zeta_850': min_max_zeta,
        'min_hgt_850': min_hgt,
        'max_wind_850': max_wind
        }

        out_track = out_track.append(position, ignore_index=True)

        plot_domain_attributes(data850, position, FigsDirectory)

        print(f'\nTime: {datestr}')
        print(f'central lat/lon: {central_lon}, {central_lat}')
        print(f'Box min_lon, max_lon: {min_lon}/{max_lon}')
        print(f'Box min_lat, max_lat: {min_lat}/{max_lat}')
        print(f'Box size (longitude): {width}')
        print(f'Box size (latitude): {length}')
        print(f'Minimum/Maximum vorticity at 850 hPa: {min_max_zeta}')
        print(f'Minimum geopotential height at 850 hPa: {min_hgt}')
        print(f'Maximum wind speed at 850 hPa: {max_wind}')
    
        # Create box object
        try:
            box_obj = BoxData(
                data=idata.compute(),
                dfVars=dfVars,
                args=args,
                western_limit=min_lon,
                eastern_limit=max_lon,
                southern_limit=min_lat,
                northern_limit=max_lat,
                output_dir=ResultsSubDirectory,
                dTdt=idTdt
            )
        except Exception as e:
            logging.exception("An exception occurred: {}".format(e))
            print('An exception occurred: {}'.format(e))
            raise SystemExit('Error creating the box for computations')
        
        # Compute energy terms
        try:
            ec_obj = EnergyContents(box_obj,method='moving')
            TermsDict['Az'].append(ec_obj.calc_az())
            TermsDict['Ae'].append(ec_obj.calc_ae())
            TermsDict['Kz'].append(ec_obj.calc_kz())
            TermsDict['Ke'].append(ec_obj.calc_ke())
        except Exception as e:
            logging.exception("An exception occurred: {}".format(e))
            print('An exception occurred: {}'.format(e))
            raise SystemExit('Error on computing Energy Contents')
        
        # Compute conversion terms
        try:
            ct_obj = ConversionTerms(box_obj,method='moving')
            TermsDict['Ca'].append(ct_obj.calc_ca())
            TermsDict['Ce'].append(ct_obj.calc_ce())
            TermsDict['Ck'].append(ct_obj.calc_ck())
            TermsDict['Cz'].append(ct_obj.calc_cz())
        except Exception as e:
            logging.exception("An exception occurred: {}".format(e))
            print('An exception occurred: {}'.format(e))
            raise SystemExit('Error on computing Conversion Terms')
        
        # Compute boundary terms
        try:
            bt_obj = BoundaryTerms(box_obj,method='moving')
            TermsDict['BAe'].append(bt_obj.calc_bae())
            TermsDict['BAz'].append(bt_obj.calc_baz())
            TermsDict['BKe'].append(bt_obj.calc_bke())
            TermsDict['BKz'].append(bt_obj.calc_bkz())
            TermsDict['BΦE'].append(bt_obj.calc_boe())
            TermsDict['BΦZ'].append(bt_obj.calc_boz())
        except Exception as e:
            logging.exception("An exception occurred: {}".format(e))
            print('An exception occurred: {}'.format(e))
            raise SystemExit('Error on computing Boundary Terms')
        
        # Compute generation/dissipation terms
        try:
            gdt_obj = GenerationDissipationTerms(box_obj,method='moving')
            TermsDict['Ge'].append(gdt_obj.calc_ge())
            TermsDict['Gz'].append(gdt_obj.calc_gz())
            if not args.residuals:
                TermsDict['De'].append(gdt_obj.calc_de())
                TermsDict['Dz'].append(gdt_obj.calc_dz())
        except Exception as e:
            logging.exception("An exception occurred: {}".format(e))
            print('An exception occurred: {}'.format(e))
            raise SystemExit('Error on computing Generation Terms')
        
    print('\nOrganising results in a Pandas DataFrame')
    df = pd.DataFrame.from_dict(TermsDict, orient ='columns',dtype = float)
    days = times.values.astype('datetime64[D]')
    hours = pd.to_datetime(times.values).hour
    df['Date'],df['Hour'] = days, hours
    
    # Print results for user
    if args.verbosity:
        for term in TermsDict.keys():
            print('\n'+term)
            print(df[term].tolist())
    
    print('\n------------------------------------------------------------------------')
    print('Estimating budget terms (∂X/∂t) using finite differences ')
    df = calc_budget_diff(df,times, args) 
    print('Ok!')
    
    print('\n------------------------------------------------------------------------')
    print('Computing residuals RGz, RKz, RGe and RKe')
    df = calc_residuals(df, args)
    print('Ok!')
    
    print('\nCreating a csv to store results...')
    outfile = ResultsSubDirectory+'/'+outfile_name+'.csv'
    df.to_csv(outfile)
    print(outfile+' created') 
    print('All done!')
    
    # Save system position as a csv file for replicability
    # df = pd.DataFrame.from_dict(position, orient='index').T
    out_track = out_track.rename(columns={'datestr':'time','central_lat':'Lat','central_lon':'Lon'})
    output_trackfile =  ResultsSubDirectory+outfile_name+'_track'
    out_track.to_csv(output_trackfile, index=False, sep=";")
    
    # 13) Make figures directory
    FigsDirectory = ResultsSubDirectory+'/Figures'
    check_create_folder(FigsDirectory)
    
    # Determine periods
    try:
        determine_periods_options = {
            "vorticity_column": 'min_max_zeta_850',
            "plot": os.path.join(ResultsSubDirectory, 'periods'),
            "plot_steps": os.path.join(ResultsSubDirectory, 'periods_didatic'),
            "export_dict": os.path.join(ResultsSubDirectory, 'periods'),
            "process_vorticity_args": {
                "use_filter": "auto",
                "use_smoothing_twice": "auto"}
        }
        determine_periods(output_trackfile, **determine_periods_options)
    except Exception as e:
        logging.exception("An exception occurred: {}".format(e))
        print('An exception occurred: {}'.format(e))
        raise SystemExit('Error on determining periods')
    
    if args.residuals:
        flag = ' -r'
    else:
        flag = ' '

    if args.plots:
        os.system("python ../plots/plot_timeseries.py "+outfile+flag)
        os.system("python ../plots/plot_vertical.py "+ResultsSubDirectory)
        os.system("python ../plots/plot_boxplot.py "+ResultsSubDirectory+flag)
        os.system("python ../plots/plot_LEC.py "+outfile)
        os.system("python ../plots/plot_LPS.py "+outfile)
        os.system("python ../plots/plot_track.py "+outfile)