# **************************************************************************** #
#                                                                              #
#                                                         :::      ::::::::    #
#    lorenz_cycle.py                                    :+:      :+:    :+:    #
#                                                     +:+ +:+         +:+      #
#    By: daniloceano <danilo.oceano@gmail.com>      +#+  +:+       +#+         #
#                                                 +#+#+#+#+#+   +#+            #
#    Created: 2022/04/20 10:05:52 by daniloceano       #+#    #+#              #
#    Updated: 2023/12/19 18:04:21 by daniloceano      ###   ########.fr        #
#                                                                              #
# **************************************************************************** #

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Main script for the Lorenz Energy Cycle (LEC) Program.

This script handles user inputs to open the data, specifies the bounding box for the calculations,
and validates user inputs for processing the Lorenz Energy Cycle using either a fixed or moving framework.

Created by Danilo Couto de Souza, Universidade de São Paulo (USP), Instituto de Astronomia, Ciências Atmosféricas e Geociências, São Paulo - Brazil
Contact: danilo.oceano@gmail.com
"""

import os
import pandas as pd
import argparse
import time
import logging
from metpy.units import units
from select_area import slice_domain
from tools import get_data
from lec_fixed_framework import lec_fixed
from lec_moving_framework import lec_moving

def main():    
    parser = argparse.ArgumentParser(description="Lorenz Energy Cycle (LEC) program.")
    parser.add_argument("infile", help="Input .nc file with temperature, geopotential, and wind components.")
    parser.add_argument("-r", "--residuals", action='store_true', help="Compute the Dissipation and Generation terms as residuals.")
    parser.add_argument("-g", "--geopotential", action='store_true', help="Use geopotential data instead of geopotential height.")
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("-f", "--fixed", action='store_true', help="Compute the energetics for a fixed domain specified by the 'box_lims' file.")
    group.add_argument("-t", "--track", action='store_true', help="Define the box using a track file specified by the 'track' file.")
    group.add_argument("-c", "--choose", action='store_true', help="Choose the domain for each time step by clicking on the screen.")
    parser.add_argument("-z", "--zeta", action='store_true', help="Use this flag if the track file was created using vorticity.")
    parser.add_argument("-m", "--mpas", action='store_true', help="for MPAS-A data processed with MPAS-BR routines")
    parser.add_argument("-p", "--plots", action='store_true', help="wether or not to make plots.")
    parser.add_argument("-o", "--outname", type=str, help="Choose a name for saving results.")
    parser.add_argument("-v", "--verbosity", action='store_true', help="Increase output verbosity.")

    args = parser.parse_args()

    # Debug:
    # args = parser.parse_args(['../samples/Reg1-Representative_NCEP-R2.nc', '-r', '-t', '-p'])
    # args = parser.parse_args(['/p1-nemo/danilocs/mpas/MPAS-BR/post_proc/py/interpolations/Catarina-2403-2903_MPAS.nc',
    # '-t', '-g', '-r'])
    # args = parser.parse_args(['/p1-nemo/danilocs/SWSA-cyclones_energetic-analysis/met_data/ERA5/DATA/10MostIntense-19830422_ERA5.nc',
    # '-c', '-g', '-r'])
    # args = parser.parse_args(['../../data_etc/netCDF_data/2022.nc', '-r', '-f', '-p'])
    # args = parser.parse_args(['../samples/Reg1-Representative_NCEP-R2.nc', '-r', '-f', '-p'])

    initialize_logging()
    data = prepare_data(args)
    
    run_lec_analysis(data, args)

def initialize_logging():
    log_file = 'error_log.txt'
    logging.basicConfig(filename=log_file, filemode='w', level=logging.ERROR, format='%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    if os.path.getsize(log_file) == 0:
        os.remove(log_file)

def prepare_data(args):
    data = get_data(args.infile, '../inputs/fvars')
    if args.mpas:
        data = data.drop_dims('standard_height')
    
    return slice_domain(data, args, '../inputs/fvars')

def run_lec_analysis(data, args):
    start_time = time.time()
    results_directory = setup_results_directory(args)

    if args.fixed:
        lec_fixed(data, results_directory, args)
        logging.info("--- %s seconds running fixed framework ---" % (time.time() - start_time))

    if args.track or args.choose:
        lec_moving(data, results_directory, args)
        logging.info("--- %s seconds running moving framework ---" % (time.time() - start_time))

def setup_results_directory(args):
    results_main_directory = '../LEC_Results'
    os.makedirs(results_main_directory, exist_ok=True)

    results_sub_directory = os.path.join(results_main_directory, construct_outfile_name(args))
    os.makedirs(results_sub_directory, exist_ok=True)

    return results_sub_directory

def construct_outfile_name(args):
    return args.outname if args.outname else ''.join(args.infile.split('/')[-1].split('.nc')) + '_' + args.method

if __name__ == "__main__":
    main()