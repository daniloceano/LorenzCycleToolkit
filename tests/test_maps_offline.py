"""
Maps offline (DEC-03, alternative A): verified Natural Earth copy, no download.

Scenario A - complete copy: every map is generated, Cartopy reads only the
verified copy and no network access is attempted.
Scenario B - no Natural Earth data: the maps fail with an explicit cause, the
download attempt is blocked by the session guard, and the computation is
preserved: the CSV outputs are byte-identical to those of the same case run
without figures on the same machine, and equal to the audit reference within
its registered tolerances (the reference was produced on x86_64).
The fixture that verifies the copy fails explicitly (never skips) when the
copy is absent or differs from the manifest.
"""

import hashlib
import math
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
    lectk_workdir(**inputs)
    return _run(monkeypatch, [*argv, *extra], method)


def _run(monkeypatch, argv, method):
    """Run the toolkit in the current working directory; return the results directory and the download warnings."""
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


def csv_output_names(results_dir, stem):
    """Main results, trackfile (moving framework) and every per-level CSV of a run."""
    results_dir = Path(results_dir)
    names = [f'{stem}_results.csv']
    if (results_dir / f'{stem}_trackfile').is_file():
        names.append(f'{stem}_trackfile')
    names += [f'results_vertical_levels/{p.name}' for p in sorted((results_dir / 'results_vertical_levels').glob('*.csv'))]
    return names


def assert_csv_outputs_byte_identical(results_dir, stem, other_dir, other_stem):
    """The CSV outputs of two runs are byte-identical; their file names differ only by the stem."""
    names = csv_output_names(results_dir, stem)
    assert len(names) in (24, 25)
    assert [n.replace(stem, other_stem) for n in names] == csv_output_names(other_dir, other_stem)
    for name in names:
        other = Path(other_dir) / name.replace(stem, other_stem)
        assert sha256(Path(results_dir) / name) == sha256(other), f'{name} differs from {other}'


def assert_computation_preserved(monkeypatch, results_dir, case):
    """
    The failure of the maps does not touch the computation.

    Byte for byte, the CSV outputs equal those of the same case run without
    figures on this machine. Value by value, they agree with the audit
    reference within its registered tolerances: the reference was produced
    on x86_64, where it is also byte-identical; on another platform
    (e.g. arm64) the last bits differ (tests/data/reference/PROVENANCE.md).
    """
    argv, method, stem, _ = case
    plain_dir, _ = _run(monkeypatch, [a for a in argv if a != '-p'] + ['-o', 'no_figures'], method)
    assert_csv_outputs_byte_identical(results_dir, stem, plain_dir, f'no_figures_{method}')
    reference = oc.REFERENCE_DIR / ('ERA5_fixed' if method == 'fixed' else 'ERA5_track')
    oc.assert_matches_reference(results_dir, reference, stem, 'level')


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

    # The computation is complete and unaffected by the failed maps.
    assert_computation_preserved(monkeypatch, results_dir, case)

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



def _alter_first_value(path, column, change):
    """Rewrite one cell of the first data row of a CSV in place; the rest of the file is unchanged."""
    lines = Path(path).read_text(encoding='utf-8').split('\n')
    header, cells = lines[0].split(','), lines[1].split(',')
    j = header.index(column)
    cells[j] = repr(change(float(cells[j])))
    lines[1] = ','.join(cells)
    Path(path).write_text('\n'.join(lines), encoding='utf-8')


def _copy_of_reference(tmp_path):
    reference = oc.REFERENCE_DIR / 'ERA5_fixed'
    run = tmp_path / 'run'
    shutil.copytree(reference, run)
    # positive control: an unaltered copy passes both checks
    assert_csv_outputs_byte_identical(run, 'testdata_ERA5_fixed', reference, 'testdata_ERA5_fixed')
    oc.assert_matches_reference(run, reference, 'testdata_ERA5_fixed', 'level')
    return reference, run


def test_preservation_checks_reject_a_value_beyond_the_tolerance(tmp_path):
    """Negative control: Kz changed by one part in a million fails both checks."""
    reference, run = _copy_of_reference(tmp_path)
    _alter_first_value(run / 'testdata_ERA5_fixed_results.csv', 'Kz', lambda x: x * (1 + 1e-6))
    with pytest.raises(AssertionError, match='differs'):
        assert_csv_outputs_byte_identical(run, 'testdata_ERA5_fixed', reference, 'testdata_ERA5_fixed')
    with pytest.raises(AssertionError, match='values differ from the reference'):
        oc.assert_matches_reference(run, reference, 'testdata_ERA5_fixed', 'level')


def test_byte_check_rejects_a_change_in_the_last_bit(tmp_path):
    """Negative control: one ulp passes the tolerance but not the same-machine byte check."""
    reference, run = _copy_of_reference(tmp_path)
    _alter_first_value(run / 'testdata_ERA5_fixed_results.csv', 'Kz', lambda x: math.nextafter(x, math.inf))
    oc.assert_matches_reference(run, reference, 'testdata_ERA5_fixed', 'level')
    with pytest.raises(AssertionError, match='differs'):
        assert_csv_outputs_byte_identical(run, 'testdata_ERA5_fixed', reference, 'testdata_ERA5_fixed')

@pytest.mark.parametrize('case', [FIXED, TRACK], ids=['fixed', 'track'])
def test_offline_maps_with_complete_cache(monkeypatch, lectk_workdir, naturalearth_reads, network_attempts, case):
    """--offline-maps reads the local copy: every map generated, no network attempt."""
    results_dir, downloads = run_case(monkeypatch, lectk_workdir, case, '--offline-maps')
    method = case[1]
    data = oc.check_run_status(results_dir, 'success', 0,
                               oc.products(method, **(NO_PERIODS if method == 'track' else {})))
    assert data['offline_maps'] is True
    assert len(network_attempts) == 0 and downloads == []
    assert naturalearth_reads


@pytest.mark.parametrize('case', [FIXED, TRACK], ids=['fixed', 'track'])
def test_offline_maps_without_cache_never_reach_the_network(monkeypatch, lectk_workdir, naturalearth_absent,
                                                           network_attempts, tmp_path, case):
    """
    Without data, --offline-maps fails the maps before any network access:
    the session guard records no attempt at all (it would record, and block,
    one without the option; see scenario B), nothing is announced as a
    download, and the computation is unchanged.
    """
    results_dir, downloads = run_case(monkeypatch, lectk_workdir, case, '--offline-maps')
    argv, method, stem, _ = case

    assert_computation_preserved(monkeypatch, results_dir, case)

    overrides = {name: 'failed' for name in MAPS[method]}
    if method == 'track':
        overrides.update(NO_PERIODS)
    data = oc.check_run_status(results_dir, 'partial', 3, oc.products(method, **overrides))
    for product in data['products']:
        if product['name'] in MAPS[method]:
            assert product['error']['type'] == 'NaturalEarthOfflineError', product['error']
            assert '--offline-maps forbids downloading it' in product['reason']
            assert product['written'] == []
    assert len(network_attempts) == 0
    assert downloads == []
    assert not (tmp_path / 'cartopy_download_dir').exists()


def test_offline_maps_restore_cartopy_downloaders():
    """Outside --offline-maps, and after it (also after an error), Cartopy is unchanged."""
    import cartopy
    import cartopy.io.shapereader  # noqa: F401

    from src.utils.offline_maps import NE_KEY, offline_maps

    before = dict(cartopy.config['downloaders'])
    with offline_maps(False):
        assert cartopy.config['downloaders'] == before
    with pytest.raises(KeyError):
        with offline_maps(True):
            assert cartopy.config['downloaders'][NE_KEY] is not before[NE_KEY]
            raise KeyError('error inside the run')
    assert cartopy.config['downloaders'] == before
    assert all(cartopy.config['downloaders'][k] is v for k, v in before.items())


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
