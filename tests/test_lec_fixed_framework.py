# **************************************************************************** #
#                                                                              #
#                                                         :::      ::::::::    #
#    test_lec_fixed_framework.py                        :+:      :+:    :+:    #
#                                                     +:+ +:+         +:+      #
#    By: daniloceano <danilo.oceano@gmail.com>      +#+  +:+       +#+         #
#                                                 +#+#+#+#+#+   +#+            #
#    Created: 2023/12/26 17:44:00 by daniloceano       #+#    #+#              #
#    Updated: 2023/12/28 12:54:03 by daniloceano      ###   ########.fr        #
#                                                                              #
# **************************************************************************** #

"""
Unit tests for the LEC fixed framework.

On the top level of the project, run the tests with:

    python -m unittest discover -s tests

"""

import unittest
import shutil  
import pandas as pd
import xarray as xr
import argparse
import os
from unittest.mock import patch
from src.frameworks.lec_fixed_framework import lec_fixed


class TestLECFixedFramework(unittest.TestCase):

    def setUp(self):
        # Mock arguments
        self.args = argparse.Namespace(
            infile='samples/Reg1-Representative_NCEP-R2.nc',
            residuals=True,
            fixed=True,
            geopotential=False,
            track=False,
            choose=False,
            zeta=False,
            mpas=False,
            plots=False,
            outname=None,
            verbosity=False
        )

        # Mock DataFrame for variable_list_df
        varlist = "inputs/fvars_NCEP-R2"
        self.variable_list_df = pd.read_csv(varlist, sep=';', index_col=0, header=0)

        # Create a mock xarray Dataset for data
        self.data = xr.open_dataset(self.args.infile)

        # Create the mock results directory
        self.results_subdirectory = "../mock_results_subdirectory"
        if not os.path.exists(self.results_subdirectory):
            os.makedirs(self.results_subdirectory)

    def tearDown(self):
        # Clean up: Recursively remove the mock results directory and its contents after tests
        if os.path.exists(self.results_subdirectory):
            shutil.rmtree(self.results_subdirectory)

    def mock_exists_side_effect(self, path):
        print(f"exists called with: {path}")
        return path != self.results_subdirectory

    @patch('src.frameworks.lec_fixed_framework.os.path.exists', side_effect=mock_exists_side_effect)
    @patch('src.frameworks.lec_fixed_framework.os.path.isdir', return_value=True)
    @patch('src.frameworks.lec_fixed_framework.BoxData')
    @patch('src.frameworks.lec_fixed_framework.EnergyContents')
    @patch('src.frameworks.lec_fixed_framework.ConversionTerms')
    @patch('src.frameworks.lec_fixed_framework.BoundaryTerms')
    @patch('src.frameworks.lec_fixed_framework.GenerationDissipationTerms')
    @patch('src.frameworks.lec_fixed_framework.xr.open_dataset')
    @patch('pandas.read_csv')
    @patch("builtins.open", new_callable=unittest.mock.mock_open)
    def test_lec_fixed_functionality(self, mock_open, mock_read_csv, mock_open_dataset, mock_gdt, mock_bt, mock_ct, mock_ec, mock_box, mock_path_isdir, mock_path_exists):
        # Mocking xr.open_dataset to return a specific xarray.Dataset
        mock_open_dataset.return_value = self.data

        # Create a mock DataFrame with the correct structure
        mock_df = pd.DataFrame({'Value': [-60, -30, -42.5, -17.5]}, index=['min_lon', 'max_lon', 'min_lat', 'max_lat'])
        mock_read_csv.return_value = mock_df

        # Call lec_fixed
        lec_fixed(self.data, self.variable_list_df, self.results_subdirectory, self.args)

        # Assert that file writing was attempted
        mock_open.assert_called()

        # Additional assertions for other mocked functions
        mock_box.assert_called_once()
        mock_ec.assert_called_once()
        mock_ct.assert_called_once()
        mock_bt.assert_called_once()
        mock_gdt.assert_called_once()

if __name__ == '__main__':
    unittest.main()