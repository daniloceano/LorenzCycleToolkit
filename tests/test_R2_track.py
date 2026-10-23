from pathlib import Path

from lorenzcycletoolkit import (create_arg_parser, initialize_logging,
                                prepare_data, run_lec_analysis,
                                setup_results_directory)

from tests import output_checks as oc

# TO DO: slice nc file for smaller time steps and improve testing


def test_different_args(monkeypatch, lectk_workdir, naturalearth_cache):
    # Prepare input data in a private working directory: the
    # repository inputs/ (user configuration) is never touched.
    lectk_workdir(namelist='namelist_NCEP-R2', track='track_testdata_NCEP-R2')

    test_args = ['lorenzcycletoolkit.py', 'samples/testdata_NCEP-R2.nc', '-r', '-t', '-p', '-v']
    monkeypatch.setattr('sys.argv', test_args)

    parser = create_arg_parser()
    args = parser.parse_args()

    assert args.residuals
    assert args.track
    assert args.plots
    assert args.verbosity
    assert args.infile == 'samples/testdata_NCEP-R2.nc'

    method = 'track'

    # Setup results directory
    results_subdirectory, figures_directory, results_subdirectory_vertical_levels = setup_results_directory(args, method)

    # Initialize logging
    app_logger = initialize_logging(results_subdirectory, args)
    app_logger.info("Starting LEC analysis")
    app_logger.info(f"Command line arguments: {args}")

    # Prepare data
    data = prepare_data(args, 'inputs/namelist', app_logger)

    # Run LEC analysis
    run_lec_analysis(data, args, results_subdirectory, figures_directory, results_subdirectory_vertical_levels, app_logger)

    # Verify the outputs (E0-08): required files, time steps, terms and finite
    # values. No approved numerical reference exists for the NCEP-R2 cases.
    results_dir = Path(results_subdirectory)
    times = oc.track_times('inputs/track')
    assert len(times) == 5
    oc.check_main_results(results_dir, 'testdata_NCEP-R2_track', times)
    oc.check_trackfile(results_dir, 'testdata_NCEP-R2_track', times)
    oc.check_vertical_levels(results_dir, 'lv_ISBL3', 'initial_time0_hours', times, 'track')
    oc.check_log(results_dir, 'testdata_NCEP-R2')

    # Run status: every requested figure, including the periods and the
    # Lorenz Phase Space, generated; success declared.
    oc.check_run_status(results_dir, 'success', 0, oc.products('track'))
    oc.check_final_messages(results_dir, 'testdata_NCEP-R2', complete=True)
    assert (results_dir / 'periods.csv').is_file()
