# D-Bench runner (CDS × NCAR)

Controlled benchmark of the ERA5 acquisition through the CDS API (`--cdsapi`,
guarded mode) and NCAR GDEX (`--ncar`), followed by the LEC computation in the
moving framework (`-t -r`), for the ten registered tracks of 1979–1980.

Protocol, tolerances, order, budget and stop rules are fixed **before any
result** in `docs/evolution/evidence/records/2026-10-29-dbench-preregistration.md`.
This directory is not part of the installed package.

## What a run does

`dbench.py run` makes one **cold** acquisition of one track by one provider:

- a new private directory `ROOT/runs/<case>/<provider>/r<rep>/` that must not
  exist yet (no earlier file can be taken for a new download);
- its own `inputs/namelist` and `inputs/namelist_ERA5-cdsapi`, both copies of
  the versioned CDS preset, so that the preparation, the `--ncar` output names
  and the LEC computation read the same mapping for both providers, whatever
  the working `inputs/namelist` of the checkout is;
- the track copied after checking its registered SHA-256;
- the toolkit under `/usr/bin/time -v` (peak RSS, CPU), with an external
  monotonic wall clock and `TMPDIR` inside the run directory;
- for CDS: guarded mode (deadline, few retries, one submission at a time,
  ledger before each submission, best-effort delete of the own job after a
  deadline), the credential given **only** through the environment of that
  process (`CDSAPI_URL`/`CDSAPI_KEY`), never `~/.cdsapirc`; before each CDS
  case the campaign lists the account's pending requests and submits nothing
  if any is pending or the listing fails;
- the input NetCDF is kept (no `--delete-downloaded-input`) for the
  integrity comparison and for reprocessing.

The record of each run (`ROOT/records/<case>_<provider>_r<rep>.json`, mode
600) holds the timing of `run_status.json`, the outcome (success, censored or
failed), file size and SHA-256, resources, load average and the commit. No
credential is ever written.

## Commands

```bash
python benchmarks/d_bench/dbench.py plan
python benchmarks/d_bench/dbench.py cds-pending --key-file KEYS          # read-only listing
python benchmarks/d_bench/dbench.py run --case 19800296 --provider ncar --rep 1 --root ROOT --tracks TRACKS
python benchmarks/d_bench/dbench.py run --case 19800296 --provider cds --rep 1 --root ROOT --tracks TRACKS --key-file KEYS
python benchmarks/d_bench/dbench.py compare --case 19800296 --rep 1 --root ROOT --tracks TRACKS
python benchmarks/d_bench/dbench.py campaign --root ROOT --tracks TRACKS --key-file KEYS
```

`ROOT`, `TRACKS` and `KEYS` are private paths outside Git. After a CDS stop,
`python -m src.utils.cds_guard status ROOT/private/cds_ledger.jsonl` says why;
only an operator who has checked the account runs
`resolve --confirm-no-pending --note "..."`.

Campaign limits: 12 h per case and provider; a 96 h ceiling counted from the
first start (persisted in `ROOT/records/campaign_start.json`), a run starting
only if its deadline still fits; r2 only after the ten r1 pairs succeeded; a
campaign that cannot complete is reported as partial. Limits are never raised
automatically. `summary` applies the pre-registered decision rule.

## Tests

`tests/test_dbench_runner.py`, `tests/test_cds_guard.py` and
`tests/test_input_cleanup.py` run offline: the toolkit process, the CDS
client and the job listing are fakes; no request is ever sent.
