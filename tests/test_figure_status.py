"""
Figure failures are reported, never hidden, and never hide a scientific error.

Each test runs a sample case with --plots in a private working directory and
injects a failure: a plotting function that raises, one that returns without
writing its figure, or an error inside the scientific computation. The run
status, the log messages and the CSV outputs are then checked.
"""

from pathlib import Path

import pytest

from lorenzcycletoolkit import (create_arg_parser, initialize_logging,
                                prepare_data, run_lec_analysis,
                                setup_results_directory)

from tests import output_checks as oc


def run_case(monkeypatch, argv, method):
    monkeypatch.setattr('sys.argv', ['lorenzcycletoolkit.py', *argv])
    args = create_arg_parser().parse_args(argv)
    results, figures, levels = setup_results_directory(args, method)
    app_logger = initialize_logging(results, args)
    app_logger.info("Starting LEC analysis")
    app_logger.info(f"Command line arguments: {args}")
    data = prepare_data(args, 'inputs/namelist', app_logger)
    status = run_lec_analysis(data, args, results, figures, levels, app_logger)
    return Path(results), status


FIXED_ARGV = ['samples/testdata_ERA5.nc', '-r', '-f', '-p']


def test_fixed_figure_failures_give_partial_output(monkeypatch, lectk_workdir):
    lectk_workdir(namelist='namelist_ERA5', box_limits='box_limits_Reg1')

    def broken_hovmoller(*args, **kwargs):
        raise RuntimeError("injected Hovmöller failure")

    def silent_box_limits(*args, **kwargs):
        return None  # logs nothing and writes nothing, like the old silent returns

    monkeypatch.setattr('src.plots.plot_hovmoller.plot_hovmoller', broken_hovmoller)
    monkeypatch.setattr('src.plots.map_box_limits.plot_box_limits', silent_box_limits)

    results_dir, status = run_case(monkeypatch, FIXED_ARGV, 'fixed')

    # The scientific outputs are complete and unchanged.
    times = oc.input_times('samples/testdata_ERA5.nc', 'time')
    oc.check_main_results(results_dir, 'testdata_ERA5_fixed', times)
    oc.assert_matches_reference(results_dir, oc.REFERENCE_DIR / 'ERA5_fixed', 'testdata_ERA5_fixed', 'level')

    data = oc.check_run_status(
        results_dir, 'partial', 3,
        oc.products('fixed', hovmoller='failed', box_limits_map='failed'),
    )
    by_name = {p['name']: p for p in data['products']}
    assert by_name['hovmoller']['error'] == {'type': 'RuntimeError', 'message': 'injected Hovmöller failure'}
    assert by_name['hovmoller']['missing'] == ['Figures/hovmollers/hovmoller_energy.png']
    assert by_name['box_limits_map']['error'] is None
    assert 'returned without writing Figures/box_limits.png' in by_name['box_limits_map']['reason']
    assert data['scientific']['status'] == 'success'
    assert status.exit_code == 3
    oc.check_final_messages(results_dir, 'testdata_ERA5', complete=False)


def test_fixed_scientific_error_is_not_partial(monkeypatch, lectk_workdir):
    lectk_workdir(namelist='namelist_ERA5', box_limits='box_limits_Reg1')

    def broken_kz(self):
        raise FloatingPointError("injected error in Kz")

    monkeypatch.setattr('src.analysis.energy_contents.EnergyContents.calc_kz', broken_kz)

    with pytest.raises(FloatingPointError, match="injected error in Kz"):
        run_case(monkeypatch, FIXED_ARGV, 'fixed')

    results_dir = Path('LEC_Results/testdata_ERA5_fixed')
    assert not (results_dir / 'testdata_ERA5_fixed_results.csv').exists()
    data = oc.load_run_status(results_dir)
    assert (data['status'], data['exit_code']) == ('failed', 1)
    assert data['scientific']['error'] == {'type': 'FloatingPointError', 'message': 'injected error in Kz'}
    assert not [p['name'] for p in data['products'] if p['status'] == 'success']
    log = (results_dir / 'log.testdata_ERA5').read_text(encoding='utf-8')
    assert 'Analysis complete' not in log
