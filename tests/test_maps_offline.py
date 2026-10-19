"""
Maps offline (DEC-03, alternative A): verified Natural Earth copy, no download.

Scenario A - complete copy: every map is generated, Cartopy reads only the
verified copy and no network access is attempted.
Scenario B - no Natural Earth data: the maps fail with an explicit cause, the
download attempt is blocked by the session guard, and the computation is
preserved: the CSV outputs are byte-identical to the reference.
The fixture that verifies the copy fails explicitly (never skips) when the
copy is absent or differs from the manifest.
"""

import hashlib
import shutil
import warnings
from pathlib import Path

import pytest
from cartopy.io import DownloadWarning

from lorenzcycletoolkit import (create_arg_parser, initialize_logging,
                                prepare_data, run_lec_analysis,
                                setup_results_directory)

from tests import output_checks as oc
from tests.naturalearth import (NaturalEarthCacheError,
                                cache_from_environment, read_manifest,
                                verify_cache)

NATURALEARTH_MANIFEST = Path(__file__).resolve().parent / 'data' / 'naturalearth_manifest.tsv'

FIXED = (['samples/testdata_ERA5.nc', '-r', '-f', '-p'], 'fixed', 'testdata_ERA5_fixed',
         dict(namelist='namelist_ERA5', box_limits='box_limits_Reg1'))
TRACK = (['samples/testdata_ERA5.nc', '-r', '-t', '-p'], 'track', 'testdata_ERA5_track',
         dict(namelist='namelist_ERA5', track='track_testdata_ERA5'))
MAPS = {'fixed': ('box_limits_map',), 'track': ('domain_maps', 'track_map')}
# The 5-step ERA5 sample is below the 8 steps of the high-resolution
# cyclophaser options (E0-06): periods and dependent products unavailable.
NO_PERIODS = dict(periods='unavailable', lec_periods='unavailable', lps='unavailable')


def run_case(monkeypatch, lectk_workdir, case, *extra):
    argv, method, stem, inputs = case
    argv = [*argv, *extra]
    lectk_workdir(**inputs)
    monkeypatch.setattr('sys.argv', ['lorenzcycletoolkit.py', *argv])
    args = create_arg_parser().parse_args(argv)
    results, figures, levels = setup_results_directory(args, method)
    app_logger = initialize_logging(results, args)
    app_logger.info("Starting LEC analysis")
    app_logger.info(f"Command line arguments: {args}")
    data = prepare_data(args, 'inputs/namelist', app_logger)
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter('always')
        run_lec_analysis(data, args, results, figures, levels, app_logger)
    return Path(results), [w for w in caught if issubclass(w.category, DownloadWarning)]


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def assert_csv_outputs_identical_to_reference(results_dir, stem, reference):
    """The CSV outputs are byte-identical to the reference copies of the audit baseline."""
    names = [f'{stem}_results.csv'] + (['testdata_ERA5_track_trackfile'] if 'track' in stem else [])
    names += [f'results_vertical_levels/{p.name}' for p in sorted((reference / 'results_vertical_levels').iterdir())]
    assert len(names) in (24, 25)
    for name in names:
        assert sha256(results_dir / name) == sha256(reference / name), f'{name} differs from the reference'


@pytest.mark.parametrize('case', [FIXED, TRACK], ids=['fixed', 'track'])
def test_scenario_a_complete_cache_maps_generated_offline(monkeypatch, lectk_workdir, naturalearth_reads,
                                                          network_attempts, case):
    results_dir, downloads = run_case(monkeypatch, lectk_workdir, case)
    method = case[1]
    overrides = NO_PERIODS if method == 'track' else {}
    oc.check_run_status(results_dir, 'success', 0, oc.products(method, **overrides))
    assert len(network_attempts) == 0
    assert downloads == []
    # Cartopy read the maps from the verified copy only
    manifest_sets = {(r['scale'], r['category'], r['name']) for r in read_manifest(NATURALEARTH_MANIFEST)}
    assert naturalearth_reads, 'no Natural Earth data were read'
    assert set(naturalearth_reads) <= manifest_sets, set(naturalearth_reads) - manifest_sets


@pytest.mark.expect_network_attempts
@pytest.mark.parametrize('case', [FIXED, TRACK], ids=['fixed', 'track'])
def test_scenario_b_no_cache_maps_fail_computation_preserved(monkeypatch, lectk_workdir, naturalearth_absent,
                                                             network_attempts, tmp_path, case):
    results_dir, downloads = run_case(monkeypatch, lectk_workdir, case)
    argv, method, stem, _ = case

    # The computation is complete and byte-identical to the reference.
    reference = oc.REFERENCE_DIR / ('ERA5_fixed' if method == 'fixed' else 'ERA5_track')
    assert_csv_outputs_identical_to_reference(results_dir, stem, reference)

    # The maps are reported as failed, with the blocked download as the cause.
    overrides = {name: 'failed' for name in MAPS[method]}
    if method == 'track':
        overrides.update(NO_PERIODS)
    data = oc.check_run_status(results_dir, 'partial', 3, oc.products(method, **overrides))
    assert data['scientific']['status'] == 'success'
    for product in data['products']:
        if product['name'] in MAPS[method]:
            assert 'network access blocked by the test suite' in product['reason'], product['reason']
            assert product['written'] == []
    if method == 'track':
        maps = next(p for p in data['products'] if p['name'] == 'domain_maps')
        assert len(maps['failures']) == 5
    oc.check_final_messages(results_dir, 'testdata_ERA5', complete=False)

    # Every attempt to reach the network was blocked, and nothing was downloaded.
    assert len(network_attempts) > 0
    assert {call for _, call, _ in network_attempts} <= {'getaddrinfo', 'connect', 'connect_ex', 'create_connection'}
    assert downloads, 'Cartopy did not announce the download it attempted'
    assert [p for p in (tmp_path / 'cartopy_download_dir').rglob('*') if p.is_file()] == []


def test_manifest_covers_every_scale_used_by_the_code():
    """Adaptive features may need 110m, 50m or 10m; the fixed features use 50m."""
    rows = read_manifest(NATURALEARTH_MANIFEST)
    assert len(rows) == 75
    sets = {(r['category'], r['name'], r['scale']) for r in rows}
    assert len(sets) == 15
    for category, name in [('physical', 'coastline'), ('physical', 'land'), ('physical', 'ocean'),
                           ('cultural', 'admin_0_boundary_lines_land')]:
        assert {s for c, n, s in sets if (c, n) == (category, name)} == {'110m', '50m', '10m'}
    for name in ('admin_0_countries', 'admin_1_states_provinces_lines', 'populated_places'):
        assert ('cultural', name, '50m') in sets


def _small_copy(tmp_path, naturalearth_cache):
    """A copy of two manifest files and a manifest naming only them."""
    rows = [r for r in read_manifest(NATURALEARTH_MANIFEST) if r['relpath'].endswith(('.prj', '.cpg'))][:2]
    root = tmp_path / 'ne'
    for row in rows:
        (root / row['relpath']).parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(naturalearth_cache / row['relpath'], root / row['relpath'])
    manifest = tmp_path / 'manifest.tsv'
    columns = ['relpath', 'category', 'scale', 'name', 'bytes', 'sha256']
    manifest.write_text('# test\n' + '\t'.join(columns) + '\n'
                        + ''.join('\t'.join(r[c] for c in columns) + '\n' for r in rows), encoding='utf-8')
    return root, manifest, rows


def test_cache_check_accepts_the_verified_copy(naturalearth_cache):
    assert len(verify_cache(naturalearth_cache, NATURALEARTH_MANIFEST)) == 75


def test_cache_check_rejects_missing_and_altered_files(tmp_path, naturalearth_cache):
    root, manifest, rows = _small_copy(tmp_path, naturalearth_cache)
    assert len(verify_cache(root, manifest)) == 2

    (root / rows[0]['relpath']).unlink()
    with pytest.raises(NaturalEarthCacheError, match=f"missing {rows[0]['relpath']}"):
        verify_cache(root, manifest)

    shutil.copyfile(naturalearth_cache / rows[0]['relpath'], root / rows[0]['relpath'])
    altered = root / rows[1]['relpath']
    altered.write_bytes(altered.read_bytes() + b' ')
    with pytest.raises(NaturalEarthCacheError, match='differs from the manifest'):
        verify_cache(root, manifest)


def test_cache_check_requires_cartopy_data_dir(monkeypatch):
    monkeypatch.delenv('CARTOPY_DATA_DIR', raising=False)
    with pytest.raises(NaturalEarthCacheError, match='CARTOPY_DATA_DIR is not set'):
        cache_from_environment(NATURALEARTH_MANIFEST)
