"""
Exit status and run_status.json of the real command line (DEC-02).

Each test runs ``lorenzcycletoolkit.py`` in a subprocess, in a private
working directory, with the network blocked (tests/run_cli_offline.py):

    0  computation complete, every requested product generated
    3  computation complete, a requested figure failed (partial output)
    0  the same with --allow-partial-plots (status still partial)
    1  error of the computation (never reported as partial)
    2  invalid arguments (unchanged; no results directory, no status file)
"""

import os
import subprocess
import sys
from pathlib import Path

from tests import output_checks as oc

REPO_ROOT = Path(__file__).resolve().parents[1]
LAUNCHER = Path(__file__).resolve().parent / 'run_cli_offline.py'
BLOCKED = 'LECTK-TEST: network access blocked'


def cli(*args, env=None):
    flags = ['-I'] if sys.flags.isolated else []
    command = [sys.executable, *flags, '-B', str(LAUNCHER), str(REPO_ROOT / 'lorenzcycletoolkit.py'), *args]
    environment = dict(os.environ, MPLBACKEND='Agg', **(env or {}))
    return subprocess.run(command, cwd=os.getcwd(), env=environment, capture_output=True, text=True, timeout=900)


FIXED_RESULTS = Path('LEC_Results/testdata_ERA5_fixed')
TRACK_RESULTS = Path('LEC_Results/testdata_ERA5_track')
# 5-step ERA5 sample: below the 8 steps of the high-resolution cyclophaser options.
NO_PERIODS = dict(periods='unavailable', lec_periods='unavailable', lps='unavailable')
NO_MAPS = dict(domain_maps='failed', track_map='failed', **NO_PERIODS)


def no_naturalearth(tmp_path):
    empty = tmp_path / 'no_naturalearth'
    empty.mkdir()
    return {'CARTOPY_DATA_DIR': str(empty), 'XDG_DATA_HOME': str(tmp_path / 'xdg_data')}


def test_cli_without_plots_exits_0(lectk_workdir):
    lectk_workdir(namelist='namelist_ERA5', box_limits='box_limits_Reg1')
    run = cli('samples/testdata_ERA5.nc', '-f', '-r')
    assert run.returncode == 0, run.stderr[-2000:]
    data = oc.check_run_status(FIXED_RESULTS, 'success', 0, oc.products('fixed', figures='not_requested'))
    assert data['plots_requested'] is False
    assert data['command']['argv'][1:] == ['samples/testdata_ERA5.nc', '-f', '-r']
    assert data['log_file'] == 'log.testdata_ERA5'
    oc.check_final_messages(FIXED_RESULTS, 'testdata_ERA5', complete=True)
    assert BLOCKED not in run.stderr


def test_cli_unavailable_periods_exit_0(lectk_workdir, naturalearth_cache):
    """Periods unavailable for a 5-step series are not a failure (E0-06)."""
    lectk_workdir(namelist='namelist_ERA5', track='track_testdata_ERA5')
    run = cli('samples/testdata_ERA5.nc', '-t', '-r', '-p')
    assert run.returncode == 0, run.stderr[-2000:]
    oc.check_run_status(TRACK_RESULTS, 'success', 0, oc.products('track', **NO_PERIODS))
    oc.check_final_messages(TRACK_RESULTS, 'testdata_ERA5', complete=True)
    assert BLOCKED not in run.stderr


def test_cli_partial_output_exits_3(lectk_workdir, tmp_path):
    lectk_workdir(namelist='namelist_ERA5', track='track_testdata_ERA5')
    run = cli('samples/testdata_ERA5.nc', '-t', '-r', '-p', env=no_naturalearth(tmp_path))
    assert run.returncode == 3, run.stderr[-2000:]
    data = oc.check_run_status(TRACK_RESULTS, 'partial', 3, oc.products('track', **NO_MAPS))
    assert data['allow_partial_plots'] is False
    oc.check_final_messages(TRACK_RESULTS, 'testdata_ERA5', complete=False)
    assert BLOCKED in run.stderr  # default mode: Cartopy tries the network, the guard blocks it


def test_cli_offline_maps_exit_3_without_network_attempt(lectk_workdir, tmp_path):
    lectk_workdir(namelist='namelist_ERA5', track='track_testdata_ERA5')
    run = cli('samples/testdata_ERA5.nc', '-t', '-r', '-p', '--offline-maps', env=no_naturalearth(tmp_path))
    assert run.returncode == 3, run.stderr[-2000:]
    data = oc.check_run_status(TRACK_RESULTS, 'partial', 3, oc.products('track', **NO_MAPS))
    assert data['offline_maps'] is True
    assert BLOCKED not in run.stderr
    assert 'Downloading' not in run.stderr
    oc.assert_matches_reference(TRACK_RESULTS, oc.REFERENCE_DIR / 'ERA5_track', 'testdata_ERA5_track', 'level')


def test_cli_allow_partial_plots_exits_0_and_stays_partial(lectk_workdir, tmp_path):
    lectk_workdir(namelist='namelist_ERA5', track='track_testdata_ERA5')
    run = cli('samples/testdata_ERA5.nc', '-t', '-r', '-p', '--offline-maps', '--allow-partial-plots',
              env=no_naturalearth(tmp_path))
    assert run.returncode == 0, run.stderr[-2000:]
    data = oc.check_run_status(TRACK_RESULTS, 'partial', 0, oc.products('track', **NO_MAPS))
    assert data['allow_partial_plots'] is True and data['offline_maps'] is True
    assert BLOCKED not in run.stderr
    oc.check_final_messages(TRACK_RESULTS, 'testdata_ERA5', complete=False)
    log = (TRACK_RESULTS / 'log.testdata_ERA5').read_text(encoding='utf-8')
    assert 'exit status 0 because --allow-partial-plots was given' in log


def test_cli_map_without_natural_earth_exits_3_with_results(lectk_workdir, tmp_path):
    lectk_workdir(namelist='namelist_ERA5', box_limits='box_limits_Reg1')
    empty = tmp_path / 'no_naturalearth'
    empty.mkdir()
    run = cli('samples/testdata_ERA5.nc', '-f', '-r', '-p',
              env={'CARTOPY_DATA_DIR': str(empty), 'XDG_DATA_HOME': str(tmp_path / 'xdg_data')})
    assert run.returncode == 3, run.stderr[-2000:]
    assert BLOCKED in run.stderr
    data = oc.check_run_status(FIXED_RESULTS, 'partial', 3, oc.products('fixed', box_limits_map='failed'))
    product = next(p for p in data['products'] if p['name'] == 'box_limits_map')
    assert 'network access blocked by the test suite' in product['reason']
    oc.assert_matches_reference(FIXED_RESULTS, oc.REFERENCE_DIR / 'ERA5_fixed', 'testdata_ERA5_fixed', 'level')


def test_cli_scientific_error_exits_1_and_is_failed(lectk_workdir):
    lectk_workdir(namelist='namelist_ERA5')
    Path('inputs/box_limits_invalid').write_text('min_lon;-60\nmax_lon;-30\nmin_lat;-17.5\nmax_lat;-42.5\n')
    run = cli('samples/testdata_ERA5.nc', '-f', '-r', '-p', '--box_limits', 'inputs/box_limits_invalid')
    assert run.returncode == 1
    assert 'Traceback' in run.stderr and 'ValueError' in run.stderr
    data = oc.load_run_status(FIXED_RESULTS)
    assert (data['status'], data['exit_code']) == ('failed', 1)
    assert data['scientific']['status'] == 'failed'
    assert data['scientific']['error']['type'] == 'ValueError'
    assert not [p['name'] for p in data['products'] if p['status'] == 'success']
    assert not (FIXED_RESULTS / 'testdata_ERA5_fixed_results.csv').exists()


def test_cli_argument_error_exits_2_without_results(lectk_workdir):
    lectk_workdir(namelist='namelist_ERA5')
    run = cli('samples/testdata_ERA5.nc', '-r')
    assert run.returncode == 2
    assert 'usage:' in run.stderr
    assert not Path('LEC_Results').exists()
