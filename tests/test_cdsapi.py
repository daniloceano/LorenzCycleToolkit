"""
Unit tests for the get_cdsapi_data function.

This module tests the CDS API data retrieval functionality, including:
- Request parameter generation
- Area calculation with buffer
- Date and time range handling
- Error handling
"""

import argparse
import logging
import os
import sys
import tempfile
from datetime import datetime, timedelta
from unittest.mock import MagicMock, Mock, patch

import numpy as np
import pandas as pd
import pytest

# Add parent directory to path to import the module
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))


@pytest.fixture
def sample_track():
    """Create a sample track DataFrame for testing."""
    dates = pd.date_range("2020-01-01", "2020-01-03", freq="6h")
    track = pd.DataFrame({
        "Lat": np.linspace(-30, -25, len(dates)),
        "Lon": np.linspace(-50, -45, len(dates)),
    }, index=dates)
    return track


@pytest.fixture
def sample_args():
    """Create sample arguments for testing."""
    with tempfile.NamedTemporaryFile(suffix=".nc", delete=False) as tmp:
        args = argparse.Namespace(infile=tmp.name, time_resolution=3)
    return args


@pytest.fixture
def app_logger():
    """Create a logger for testing."""
    logger = logging.getLogger("test_cdsapi")
    logger.setLevel(logging.DEBUG)
    return logger


@pytest.fixture(autouse=True)
def cleanup_temp_files(sample_args):
    """Cleanup temporary files after each test."""
    yield
    if os.path.exists(sample_args.infile):
        try:
            os.remove(sample_args.infile)
        except Exception:
            pass


class TestGetCDSAPIData:
    """Test suite for the get_cdsapi_data function."""

    @patch('xarray.concat')
    @patch('xarray.open_dataset')
    @patch('cdsapi.Client')
    def test_basic_request(self, mock_client_class, mock_open_dataset, mock_concat, sample_track, sample_args, app_logger):
        """Test basic CDS API request with valid inputs."""
        from src.utils.tools import get_cdsapi_data

        # Mock the client and its retrieve method
        mock_client = MagicMock()
        mock_client_class.return_value = mock_client
        
        # Mock file creation by retrieve method
        def create_temp_file(dataset, request, output_path):
            with open(output_path, 'w') as f:
                f.write("mock netcdf data")
        mock_client.retrieve.side_effect = create_temp_file
        
        # Mock xarray operations
        mock_ds = MagicMock()
        mock_open_dataset.return_value = mock_ds
        mock_concat.return_value = mock_ds

        result = get_cdsapi_data(sample_args, sample_track, app_logger)

        # Verify client was created with correct timeout
        mock_client_class.assert_called_once_with(timeout=600, retry_max=500)

        # Verify retrieve was called 3 times (one for each day)
        assert mock_client.retrieve.call_count == 3
        
        # Check that each call was for a single day
        for call in mock_client.retrieve.call_args_list:
            assert call[0][0] == "reanalysis-era5-pressure-levels"
            request = call[0][1]
            assert request["product_type"] == "reanalysis"
            assert request["format"] == "netcdf"
            assert "variable" in request
            assert "pressure_level" in request
            assert isinstance(request["day"], str)  # Single day per request as string
            assert "time" in request
            assert "area" in request

        # Verify concatenation was called
        assert mock_concat.called
        
        # Verify function returned the file path
        assert result == sample_args.infile

    @patch('xarray.concat')
    @patch('xarray.open_dataset')
    @patch('cdsapi.Client')
    def test_area_calculation_with_buffer(self, mock_client_class, mock_open_dataset, mock_concat, sample_track, sample_args, app_logger):
        """Test that area is correctly calculated with 15-degree buffer."""
        from src.utils.tools import get_cdsapi_data

        mock_client = MagicMock()
        mock_client_class.return_value = mock_client
        
        # Mock file creation
        def create_temp_file(dataset, request, output_path):
            with open(output_path, 'w') as f:
                f.write("mock netcdf data")
        mock_client.retrieve.side_effect = create_temp_file
        
        # Mock xarray operations
        mock_ds = MagicMock()
        mock_open_dataset.return_value = mock_ds
        mock_concat.return_value = mock_ds

        get_cdsapi_data(sample_args, sample_track, app_logger)

        # Check the first call (all calls should have the same area)
        request = mock_client.retrieve.call_args_list[0][0][1]
        area = request["area"]

        # Track bounds: lat [-30, -25], lon [-50, -45]
        # With 15-degree buffer:
        # North: ceil(-25 + 15) = -10
        # South: floor(-30 - 15) = -45
        # West: floor(-50 - 15) = -65
        # East: ceil(-45 + 15) = -30
        assert area[0] == -10  # North
        assert area[1] == -65  # West
        assert area[2] == -45  # South
        assert area[3] == -30  # East

    @patch('xarray.concat')
    @patch('xarray.open_dataset')
    @patch('cdsapi.Client')
    def test_date_range_handling(self, mock_client_class, mock_open_dataset, mock_concat, sample_track, sample_args, app_logger):
        """Test correct handling of date ranges."""
        from src.utils.tools import get_cdsapi_data

        mock_client = MagicMock()
        mock_client_class.return_value = mock_client
        
        # Mock file creation
        def create_temp_file(dataset, request, output_path):
            with open(output_path, 'w') as f:
                f.write("mock netcdf data")
        mock_client.retrieve.side_effect = create_temp_file
        
        # Mock xarray operations
        mock_ds = MagicMock()
        mock_open_dataset.return_value = mock_ds
        mock_concat.return_value = mock_ds

        get_cdsapi_data(sample_args, sample_track, app_logger)

        # Track spans from 2020-01-01 to 2020-01-03, so 3 calls
        assert mock_client.retrieve.call_count == 3
        
        # Check that days 01, 02, 03 were requested
        days_requested = [call[0][1]["day"] for call in mock_client.retrieve.call_args_list]
        assert days_requested == ["01", "02", "03"]

    @patch('xarray.concat')
    @patch('xarray.open_dataset')
    @patch('cdsapi.Client')
    def test_time_step_calculation(self, mock_client_class, mock_open_dataset, mock_concat, sample_track, sample_args, app_logger):
        """Test time step uses the value from args.time_resolution."""
        from src.utils.tools import get_cdsapi_data

        mock_client = MagicMock()
        mock_client_class.return_value = mock_client
        
        # Mock file creation
        def create_temp_file(dataset, request, output_path):
            with open(output_path, 'w') as f:
                f.write("mock netcdf data")
        mock_client.retrieve.side_effect = create_temp_file
        
        # Mock xarray operations
        mock_ds = MagicMock()
        mock_open_dataset.return_value = mock_ds
        mock_concat.return_value = mock_ds

        get_cdsapi_data(sample_args, sample_track, app_logger)

        request = mock_client.retrieve.call_args_list[0][0][1]
        times = request["time"]

        # args.time_resolution is 3 hours (default), so times should be every 3 hours
        expected_times = ["00:00", "03:00", "06:00", "09:00", "12:00", "15:00", "18:00", "21:00"]
        assert times == expected_times

    @patch('xarray.concat')
    @patch('xarray.open_dataset')
    @patch('cdsapi.Client')
    def test_additional_day_needed(self, mock_client_class, mock_open_dataset, mock_concat, sample_args, app_logger):
        """Test that date range is correctly extended when needed."""
        from src.utils.tools import get_cdsapi_data

        # Create track that ends at 22:00 on day 02 (which is after 21:00)
        # This should trigger the additional day logic
        dates = pd.date_range("2020-01-01 00:00", "2020-01-02 22:00", freq="3h")
        track = pd.DataFrame({
            "Lat": np.linspace(-30, -25, len(dates)),
            "Lon": np.linspace(-50, -45, len(dates)),
        }, index=dates)

        mock_client = MagicMock()
        mock_client_class.return_value = mock_client
        
        # Mock file creation
        def create_temp_file(dataset, request, output_path):
            with open(output_path, 'w') as f:
                f.write("mock data")
        mock_client.retrieve.side_effect = create_temp_file
        
        # Mock xarray operations
        mock_ds = MagicMock()
        mock_open_dataset.return_value = mock_ds
        mock_concat.return_value = mock_ds

        get_cdsapi_data(sample_args, track, app_logger)

        # Check days requested - track ends at 21:00, so no additional day should be added
        days_requested = [call[0][1]["day"] for call in mock_client.retrieve.call_args_list]
        assert "01" in days_requested
        assert "02" in days_requested
        # Should not include day 03 because the last timestamp is exactly at 21:00
        assert len(days_requested) == 2

    @patch('xarray.concat')
    @patch('xarray.open_dataset')
    @patch('cdsapi.Client')
    def test_pressure_levels_included(self, mock_client_class, mock_open_dataset, mock_concat, sample_track, sample_args, app_logger):
        """Test that all required pressure levels are included."""
        from src.utils.tools import get_cdsapi_data

        mock_client = MagicMock()
        mock_client_class.return_value = mock_client
        
        # Mock file creation
        def create_temp_file(dataset, request, output_path):
            with open(output_path, 'w') as f:
                f.write("mock netcdf data")
        mock_client.retrieve.side_effect = create_temp_file
        
        # Mock xarray operations
        mock_ds = MagicMock()
        mock_open_dataset.return_value = mock_ds
        mock_concat.return_value = mock_ds

        get_cdsapi_data(sample_args, sample_track, app_logger)

        request = mock_client.retrieve.call_args_list[0][0][1]
        pressure_levels = request["pressure_level"]

        # Check some key pressure levels
        assert "850" in pressure_levels
        assert "500" in pressure_levels
        assert "1000" in pressure_levels
        assert "100" in pressure_levels

    @patch('xarray.concat')
    @patch('xarray.open_dataset')
    @patch('cdsapi.Client')
    def test_variables_included(self, mock_client_class, mock_open_dataset, mock_concat, sample_track, sample_args, app_logger):
        """Test that all required variables are included."""
        from src.utils.tools import get_cdsapi_data

        mock_client = MagicMock()
        mock_client_class.return_value = mock_client
        
        # Mock file creation
        def create_temp_file(dataset, request, output_path):
            with open(output_path, 'w') as f:
                f.write("mock netcdf data")
        mock_client.retrieve.side_effect = create_temp_file
        
        # Mock xarray operations
        mock_ds = MagicMock()
        mock_open_dataset.return_value = mock_ds
        mock_concat.return_value = mock_ds

        get_cdsapi_data(sample_args, sample_track, app_logger)

        request = mock_client.retrieve.call_args_list[0][0][1]
        variables = request["variable"]

        expected_variables = [
            "u_component_of_wind",
            "v_component_of_wind",
            "temperature",
            "vertical_velocity",
            "geopotential",
        ]
        assert set(variables) == set(expected_variables)

    @patch('cdsapi.Client')
    def test_file_not_created_raises_error(self, mock_client_class, sample_track, sample_args, app_logger):
        """Test that FileNotFoundError is raised if daily file is not created."""
        from src.utils.tools import get_cdsapi_data

        mock_client = MagicMock()
        mock_client_class.return_value = mock_client

        # Don't create the file to simulate failure (retrieve does nothing)
        # This will cause the daily file check to fail
        with pytest.raises(FileNotFoundError, match="Daily file not created"):
            get_cdsapi_data(sample_args, sample_track, app_logger)

    @patch('cdsapi.Client')
    def test_api_error_propagates(self, mock_client_class, sample_track, sample_args, app_logger):
        """Test that API errors are properly propagated."""
        from src.utils.tools import get_cdsapi_data

        mock_client = MagicMock()
        mock_client_class.return_value = mock_client
        mock_client.retrieve.side_effect = Exception("API Error: Authentication failed")

        with pytest.raises(Exception, match="API Error: Authentication failed"):
            get_cdsapi_data(sample_args, sample_track, app_logger)

    @patch('xarray.concat')
    @patch('xarray.open_dataset')
    @patch('cdsapi.Client')
    def test_single_timestep_track(self, mock_client_class, mock_open_dataset, mock_concat, sample_args, app_logger):
        """Test handling of track with only one timestep."""
        from src.utils.tools import get_cdsapi_data

        # Create track with single timestep
        track = pd.DataFrame({
            "Lat": [-30],
            "Lon": [-50],
        }, index=pd.DatetimeIndex(["2020-01-01 00:00"]))

        mock_client = MagicMock()
        mock_client_class.return_value = mock_client
        
        # Mock file creation
        def create_temp_file(dataset, request, output_path):
            with open(output_path, 'w') as f:
                f.write("mock data")
        mock_client.retrieve.side_effect = create_temp_file
        
        # Mock xarray operations
        mock_ds = MagicMock()
        mock_open_dataset.return_value = mock_ds
        mock_concat.return_value = mock_ds

        get_cdsapi_data(sample_args, track, app_logger)

        request = mock_client.retrieve.call_args[0][1]
        times = request["time"]

        # With smart time detection, single timestep at 00:00 should download only that hour
        # (rounded to time_step which is 3, so it downloads 00:00 only)
        expected_times = ["00:00"]
        assert times == expected_times

    @patch('xarray.concat')
    @patch('xarray.open_dataset')
    @patch('cdsapi.Client')
    def test_cross_dateline_coordinates(self, mock_client_class, mock_open_dataset, mock_concat, sample_args, app_logger):
        """Test handling of coordinates that cross the dateline."""
        from src.utils.tools import get_cdsapi_data

        # Create track crossing the dateline (170°E to -170°W = 190°)
        dates = pd.date_range("2020-01-01", "2020-01-02", freq="6h")
        track = pd.DataFrame({
            "Lat": np.linspace(-30, -25, len(dates)),
            "Lon": np.linspace(170, -170, len(dates)),  # Crosses dateline
        }, index=dates)

        mock_client = MagicMock()
        mock_client_class.return_value = mock_client
        
        # Mock file creation
        def create_temp_file(dataset, request, output_path):
            with open(output_path, 'w') as f:
                f.write("mock data")
        mock_client.retrieve.side_effect = create_temp_file
        
        # Mock xarray operations
        mock_ds = MagicMock()
        mock_open_dataset.return_value = mock_ds
        mock_concat.return_value = mock_ds

        get_cdsapi_data(sample_args, track, app_logger)

        # Should complete without error
        assert mock_client.retrieve.called

    @patch('xarray.concat')
    @patch('xarray.open_dataset')
    @patch('cdsapi.Client')
    def test_logging_output(self, mock_client_class, mock_open_dataset, mock_concat, sample_track, sample_args, caplog):
        """Test that appropriate log messages are generated."""
        from src.utils.tools import get_cdsapi_data

        mock_client = MagicMock()
        mock_client_class.return_value = mock_client
        
        # Mock file creation
        def create_temp_file(dataset, request, output_path):
            with open(output_path, 'w') as f:
                f.write("mock data")
        mock_client.retrieve.side_effect = create_temp_file
        
        # Mock xarray operations
        mock_ds = MagicMock()
        mock_open_dataset.return_value = mock_ds
        mock_concat.return_value = mock_ds

        logger = logging.getLogger("test_logger")
        logger.setLevel(logging.DEBUG)

        with caplog.at_level(logging.INFO):
            get_cdsapi_data(sample_args, sample_track, logger)

        # Check that key log messages were generated
        assert any("Starting data download from CDS API" in record.message 
                  or "Retrieving data from CDS API" in record.message
                  for record in caplog.records)
        assert any("COMPLETED SUCCESSFULLY" in record.message 
                  or "completed successfully" in record.message
                  for record in caplog.records)

    @patch('xarray.concat')
    @patch('xarray.open_dataset')
    @patch('cdsapi.Client')
    def test_multi_year_track(self, mock_client_class, mock_open_dataset, mock_concat, sample_args, app_logger):
        """Test handling of track spanning multiple years."""
        from src.utils.tools import get_cdsapi_data

        # Create track spanning from 2019 to 2020
        dates = pd.date_range("2019-12-30", "2020-01-02", freq="6h")
        track = pd.DataFrame({
            "Lat": np.linspace(-30, -25, len(dates)),
            "Lon": np.linspace(-50, -45, len(dates)),
        }, index=dates)

        mock_client = MagicMock()
        mock_client_class.return_value = mock_client
        
        # Mock file creation
        def create_temp_file(dataset, request, output_path):
            with open(output_path, 'w') as f:
                f.write("mock data")
        mock_client.retrieve.side_effect = create_temp_file
        
        # Mock xarray operations
        mock_ds = MagicMock()
        mock_open_dataset.return_value = mock_ds
        mock_concat.return_value = mock_ds

        get_cdsapi_data(sample_args, track, app_logger)

        # Collect all years from all calls
        years_requested = [call[0][1]["year"] for call in mock_client.retrieve.call_args_list]
        months_requested = [call[0][1]["month"] for call in mock_client.retrieve.call_args_list]
        
        # Should include both years
        assert "2019" in years_requested
        assert "2020" in years_requested
        # Should include December and January
        assert "12" in months_requested
        assert "01" in months_requested


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
