import shutil

from lorenzcycletoolkit import (create_arg_parser, initialize_logging,
                                prepare_data, run_lec_analysis,
                                setup_results_directory)

# TO DO: slice nc file for smaller time steps and improve testing


def test_different_args(monkeypatch):
    # Prepare input data
    shutil.copy('inputs/track_testdata_NCEP-R2', 'inputs/track')
    shutil.copy('inputs/namelist_NCEP-R2', 'inputs/namelist')

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
