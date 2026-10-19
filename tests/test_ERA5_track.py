from pathlib import Path

from lorenzcycletoolkit import (create_arg_parser, initialize_logging,
                                prepare_data, run_lec_analysis,
                                setup_results_directory)

from tests import output_checks as oc

# TO DO: slice nc file for smaller time steps and improve testing


def test_different_args(monkeypatch, lectk_workdir, naturalearth_cache):
    # Prepare input data in a private working directory: the
    # repository inputs/ (user configuration) is never touched.
    lectk_workdir(namelist='namelist_ERA5', track='track_testdata_ERA5')

    test_args = ['lorenzcycletoolkit.py', 'samples/testdata_ERA5.nc', '-r', '-t', '-p', '-v']
    monkeypatch.setattr('sys.argv', test_args)

    parser = create_arg_parser()
    args = parser.parse_args()

    assert args.residuals
    assert args.track
    assert args.plots
    assert args.verbosity
    assert args.infile == 'samples/testdata_ERA5.nc'

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

    # Verify the outputs (E0-08): required files, time steps, terms, finite
    # values and agreement with the audit reference outputs.
    results_dir = Path(results_subdirectory)
    times = oc.track_times('inputs/track')
    assert len(times) == 5
    assert times == oc.input_times('samples/testdata_ERA5.nc', 'time')
    oc.check_main_results(results_dir, 'testdata_ERA5_track', times)
    oc.check_trackfile(results_dir, 'testdata_ERA5_track', times)
    oc.check_vertical_levels(results_dir, 'level', 'time', times, 'track')
    oc.check_log(results_dir, 'testdata_ERA5')
    oc.assert_matches_reference(results_dir, oc.REFERENCE_DIR / 'ERA5_track', 'testdata_ERA5_track', 'level')

    # Run status: the 5 hourly steps of this sample are fewer than the 8
    # needed by the high-resolution cyclophaser options (E0-06, P1), so the
    # periods are unavailable, cyclophaser is not called, and the Lorenz cycle
    # diagrams of the periods and the Lorenz Phase Space are unavailable too.
    # Unavailable is not a failure: the run is complete (exit status 0).
    data = oc.check_run_status(
        results_dir, 'success', 0,
        oc.products('track', periods='unavailable', lec_periods='unavailable', lps='unavailable'),
    )
    by_name = {p['name']: p for p in data['products']}
    assert by_name['periods']['reason'].startswith('5 usable time steps < minimum 8 for the high_res options')
    assert by_name['periods']['error'] is None and by_name['periods']['written'] == []
    assert by_name['lec_periods']['depends_on'] == ['periods']
    assert by_name['lps']['depends_on'] == ['periods']
    assert by_name['periods']['consequences'] == ['lec_periods', 'lps']
    assert not (results_dir / 'periods.csv').exists()
    assert not (results_dir / 'Figures' / 'Periods').exists() or not any((results_dir / 'Figures' / 'Periods').iterdir())
    oc.check_final_messages(results_dir, 'testdata_ERA5', complete=True)
    log = (results_dir / 'log.testdata_ERA5').read_text(encoding='utf-8')
    assert 'Calling Cyclophaser' not in log
    assert 'All plots generated successfully!' not in log
    assert '3 requested products unavailable for these data: periods, lec_periods, lps' in log
