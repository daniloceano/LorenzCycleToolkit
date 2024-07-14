import shutil
from lorenz_cycle import create_arg_parser, initialize_logging, prepare_data, setup_results_directory, run_lec_analysis

def test_example_args(monkeypatch):

    # Prepare input data
    shutil.copy('inputs/box_limits_Reg1', 'inputs/box_limits')
    shutil.copy('inputs/namelist_ERA5', 'inputs/namelist')

    test_args = ['lorenz_cycle.py', 'samples/testdata_ERA5.nc', '-r', '-f', '-p', '-v']
    monkeypatch.setattr('sys.argv', test_args)
    
    parser = create_arg_parser()
    args = parser.parse_args()
    
    assert args.residuals == True
    assert args.fixed == True
    assert args.plots == True
    assert args.verbosity == True
    assert args.infile == 'samples/testdata_ERA5.nc'

    method = 'fixed'

    # Setup results directory
    results_subdirectory, figures_directory = setup_results_directory(args, method)

    # Initialize logging
    app_logger = initialize_logging(results_subdirectory, args)
    app_logger.info("Starting LEC analysis")
    app_logger.info(f"Command line arguments: {args}")

    # Prepare data
    data = prepare_data(args, 'inputs/namelist', app_logger)
    
    # Run LEC analysis
    run_lec_analysis(data, args, results_subdirectory, figures_directory, app_logger)