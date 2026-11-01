Remote Data Acquisition Methods
===============================

The toolkit needs ERA5 pressure-level fields over the area of a track. They can be read from a local NetCDF
file or obtained remotely. This page describes four remote methods, what each needs, how each behaves when
something fails, and what a controlled benchmark (D-Bench, October 2026) measured for each of them.

.. list-table::
   :header-rows: 1
   :widths: 24 30 18 28

   * - Method
     - Provider
     - Command line
     - Status
   * - CDS API
     - Copernicus Climate Data Store
     - ``--cdsapi``
     - supported
   * - NCAR NCSS
     - NCAR GDEX, THREDDS NetCDF Subset Service
     - ``--ncar``
     - experimental (see :doc:`ncar_gdex_example`)
   * - NCAR kerchunk byte-range
     - NCAR GDEX, kerchunk references + HTTPS range requests
     - not available
     - benchmarked with a separate client; integration under review, not released
   * - NCAR Get a Subset
     - NCAR GDEX request API (asynchronous)
     - not available
     - benchmarked with a separate client; integration under review, not released

.. note::

   **No method is best everywhere.** The benchmark below describes ten cyclone cases of 1979–1980 requested
   from one server in Brazil on one day (2026-10-31), one run per case and method. Transfer times, queues,
   failures and the volume each method moves depend on the network, the load of each service, its queues, the
   period requested and the size of the domain. Read the results as observations under those conditions.

How to read the numbers on this page
------------------------------------

- **Measured:** observed in a stated run, on a stated date. A single measurement describes that run, not the
  service in general.
- **Estimated:** computed from sizes or documented limits, not observed.
- **Hypothesis:** a possible explanation that has not been tested.

What every method must deliver
------------------------------

Whatever the route, the toolkit runs on one NetCDF file with:

- five fields: temperature (K), geopotential (m\ :sup:`2` s\ :sup:`-2`), vertical velocity (Pa s\ :sup:`-1`),
  zonal and meridional wind (m s\ :sup:`-1`);
- the 37 ERA5 pressure levels, 1 to 1000 hPa;
- the 0.25° ERA5 grid over the extent of the track with a 15° buffer, rounded outwards to whole degrees;
- the hours of the track for the chosen ``--time-resolution``;
- variable and coordinate names that match the namelist (see :doc:`configuration`).

The amount of data a case needs is

.. math::

   V = 5 \times 37 \times n_\mathrm{lat} \times n_\mathrm{lon} \times n_\mathrm{times} \times 4\ \mathrm{bytes}

in float32. The ten benchmark cases need 0.5 to 6.7 GiB. The methods differ in **how much more than V they
transfer**, **where the area and the hours are cut**, and **when the data start to arrive**.

How each method works
---------------------

CDS API (``--cdsapi``)
~~~~~~~~~~~~~~~~~~~~~~

- **Mechanism:** one asynchronous job per day, with the five fields together. The CDS queues the job, extracts
  the area and the hours on its side, converts GRIB to NetCDF and makes the file available for download. The
  toolkit waits for each job, downloads it and merges the days.
- **Requirements and authentication:** a CDS account; the ERA5 licence (CC-BY-4.0) accepted in the CDS web
  interface (without it the request is refused with HTTP 403); a personal API key in ``~/.cdsapirc``.
- **Dates:** 1940 to a few days ago.
- **Advantages:** the reference copy of ERA5; the most recent dates; the area and hours are cut on the server
  (about V transferred).
- **Disadvantages:** queue waits that vary from minutes to hours; per-account queue and cost limits
  (``--time-resolution 1`` over large areas may exceed them); protection measures the CDS applies at times.
- **Resources:** one file per day on disk, then the merged file; memory for the merge of the daily files.
- **Failures:** stopping the client does not cancel the job: it may stay queued or running and use the
  account's quota. Check the account for pending requests before submitting again.

NCAR NCSS (``--ncar``)
~~~~~~~~~~~~~~~~~~~~~~

- **Mechanism:** synchronous HTTP requests to the THREDDS NetCDF Subset Service of GDEX (dataset d633000), one
  per field and day, two when the area crosses Greenwich. The server cuts the area and picks the hours and
  answers on the same connection. The toolkit checks every response and writes one file.
- **Requirements and authentication:** none (no account).
- **Dates:** GDEX publishes ERA5 some months after the CDS (on 2026-10-26: 1940-01-01 to 2026-07-31).
- **Advantages:** no account and no queue; about V transferred (measured: 1.00 × V); verified responses are
  kept, so an interrupted retrieval resumes.
- **Disadvantages:** one long synchronous response per request; in the benchmarks some runs failed after
  repeated timeouts or HTTP errors (see below).
- **Resources:** responses in ``<file>.parts/`` until the file is written; one response in memory at a time
  (measured peak RSS 0.5–3.2 GiB, including the energetics run in the same process).
- **Limits:** set by the toolkit (8 GiB per run, 4 attempts per request, 300 s without data). Details in
  :doc:`ncar_gdex_example`.
- **Failures (measured, 2026-10-30 and 2026-10-31):** four runs failed, two on each day, on domains of 1.9 to
  3.9 GiB: a request got no response within 300 s, or HTTP 503 or 504 from the GDEX front end, four times in a
  row. The largest domain (6.7 GiB) succeeded on 2026-10-31. The cause (service load, a property of particular
  files, or both) has **not** been established.

NCAR kerchunk byte-range
~~~~~~~~~~~~~~~~~~~~~~~~

- **Mechanism:** GDEX publishes kerchunk reference files (JSON) that map every chunk of its NetCDF-4 files to
  a byte range. A chunk is one field, one hour, all 37 levels and the **whole globe** (721 × 1440 points),
  compressed with zlib and byte shuffle. The client reads only the chunks of the hours it needs with HTTPS range
  requests, decompresses them and cuts the area **locally**.
- **Requirements and authentication:** none (no account). The reference files must be obtained first (about
  16 MB per field and decade). They name a host whose TLS certificate failed on 2026-10-30; the benchmark
  client used another GDEX host serving the same files.
- **Advantages:** no account, no queue and no server-side processing; every chunk is an independent request
  that can be retried alone; no failed range in 2,190 requests in the benchmark.
- **Disadvantages:** the whole globe is transferred for every field and hour, so the volume does not shrink
  with the domain (measured: 3.0 to 11.1 × V, more for small domains).
- **Resources:** writing chunk by chunk keeps memory small (measured peak RSS 0.70–0.74 GiB); disk only for
  the output. An earlier prototype that held all fields in memory needed 13.7 GB on the largest case.
- **Failures:** a failed range is retried; the transfer can resume from the chunks already written.

NCAR Get a Subset
~~~~~~~~~~~~~~~~~

- **Mechanism:** an asynchronous request through the GDEX API: the client submits a control document
  (dataset, date range, fields, area), polls the request until it is ``Completed``, lists its output files and
  downloads them. The server cuts the area; the output is one NetCDF-4 file per field and day.
- **Hours:** the request has no time step, so every hour of the date range is delivered. The track hours are
  picked locally and the files are merged into the single file the toolkit reads.
- **Areas crossing Greenwich:** one request is enough (measured): the files hold the longitudes of both sides
  in 0–360° order and are reassembled locally.
- **Requirements and authentication:** a UCAR/GDEX account and its personal API token. The token is sent with
  every API call: keep it out of logs, shell history and shared files.
- **Advantages:** server-side area cut with compressed files (measured: 1.22 to 1.30 × V for these tracks);
  short remote processing in the benchmark (2 to 6 minutes per request).
- **Disadvantages and limits (documented by GDEX and observed):**

  - at most 8 open requests per user; a request stays for 7 days, then it is purged;
  - **identical requests are deduplicated**: while a request is listed on the account, in any state
    (including ``Set for Purge``), submitting the same request again returns the old one instead of
    processing a new one;
  - a purge is asynchronous: after ``DELETE purge/<request>`` the request stays listed as
    ``Set for Purge`` for hours;
  - a submission can take more than two minutes to be answered, and the request may be created even when the
    client stopped waiting.

- **Resources:** the files of the request (every hour of the date range) on disk, then the merged file; with
  a field-by-field merge the memory stays small (measured peak RSS 0.43–0.58 GiB).
- **Failures:** a failed download can be repeated from the file list **without a new request**. Resubmitting
  the same request does not give a fresh run (deduplication) and uses one of the 8 slots; it must never be
  "forced" by changing the scientific parameters of the request.

Cutting on the server or locally
--------------------------------

.. list-table::
   :header-rows: 1
   :widths: 20 14 16 50

   * - Method
     - Area cut
     - Hours picked
     - Data received compared with V (measured in D-Bench Stage B)
   * - CDS API
     - server
     - server
     - not measured in Stage B (about V, estimated)
   * - NCAR NCSS
     - server
     - server
     - 1.00 × V (4 cases)
   * - NCAR Get a Subset
     - server
     - local
     - 1.22–1.30 × V (9 cases): every hour of the range, compressed to about 0.43 of float32
   * - NCAR kerchunk
     - local
     - local
     - 3.0–11.1 × V (10 cases): the whole globe for each field and hour, compressed

Benchmark results (D-Bench Stage B)
-----------------------------------

Methodology
~~~~~~~~~~~

- **Date and place:** 2026-10-31, 02:30–16:11 UTC, one Linux server at IAG-USP (São Paulo, Brazil), shared with
  other work (load average about 111 on 112 cores throughout).
- **Population:** the ten cyclone tracks of 1979–1980 used by D-Bench (19 to 71 positions every 3 hours;
  domains of 0.5 to 6.7 GiB in float32).
- **Design:** one cold run per case and method (30 planned runs), one acquisition at a time. Cases ran in a
  random order (seed recorded); within each case the three methods followed an order balanced across cases so
  that each method held each position 3 or 4 times.
- **Clients:** NCSS through the toolkit's ``--ncar`` route (frozen commit, 3-hour deadline per case);
  kerchunk and Get a Subset through instrumented benchmark clients outside the toolkit. Limits were fixed before
  the first request: per-run traffic caps, 300 s per network operation, 4 attempts, a stop rule of two
  consecutive failures per method.
- **Primary measure:** time from the first request to a complete NetCDF file that passes the common checks,
  split into remote processing, transfer, local assembly and checks.
- **Equivalence:** against the file obtained from the CDS for the same case in an earlier campaign
  (2026-10-30): the same times, levels and grid, and field differences within pre-registered tolerances
  (t 0.005 K; z 2.0 m\ :sup:`2` s\ :sup:`-2`; w 0.001 Pa s\ :sup:`-1`; u, v 0.01 m s\ :sup:`-1`). The energetics
  (moving framework with residuals, ``-t -r``) were compared column by column with tolerance
  max(10\ :sup:`-4` × scale of the column, 10\ :sup:`-9`). One case (19800318) has no CDS file and is compared
  only between the NCAR methods.
- **Records:** every run, failure and stopped method is kept as observed; nothing is imputed. The protocol was
  sealed (SHA-256) before the first request and audited afterwards from the original records.

Outcome of the 30 planned runs
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. list-table::
   :header-rows: 1
   :widths: 18 13 13 13 13 30

   * - Method
     - Cold successes
     - Failed
     - Not submitted
     - Not run
     - Notes
   * - Get a Subset
     - 9
     - 0
     - 1
     - 0
     - 19800296: an identical request from an earlier test was still listed; a submission would have been
       deduplicated (not a cold run), so nothing was submitted.
   * - kerchunk
     - 10
     - 0
     - 0
     - 0
     - 0 retries in 2,190 range requests.
   * - NCSS
     - 4
     - 2
     - 0
     - 4
     - 19800256 (no response within 300 s, 4 times) and 19800318 (HTTP 504, 4 times) failed in a row, which
       stopped the method by the pre-registered rule; the last four cases were not run. Not-run cases are not
       failures, but they leave NCSS untested on them.

Success counts describe this campaign; with 6 to 10 attempts per method they do not estimate a general
failure rate.

Time to a validated file
~~~~~~~~~~~~~~~~~~~~~~~~

.. figure:: _static/remote_data_acquisition/stageB_time.png
   :alt: Minutes to a validated NetCDF file for each case (ordered by data size) and method.
   :width: 100%

   Time to a validated NetCDF file, one cold run per case and method, 2026-10-31. Cases are ordered by their
   float32 content V. × marks a failed run, – a run not made after the NCSS stop rule, "inc" the Get a Subset
   run that was not submitted.

.. list-table:: Time to a validated file, cold successes only (seconds)
   :header-rows: 1
   :widths: 20 8 18 22 32

   * - Method
     - n
     - Median
     - Quartiles
     - Range
   * - Get a Subset
     - 9
     - 934 (15.6 min)
     - 684–1410
     - 335–2269
   * - kerchunk
     - 10
     - 2440 (40.7 min)
     - 1156–3219
     - 866–4418
   * - NCSS
     - 4
     - 1524 (25.4 min)
     - 1258–2593
     - 652–5609

The medians are not comparable across methods: they cover different cases (NCSS only three of the four smallest
domains and the largest one). The comparison in the same cases is below.

.. list-table:: Paired comparison in the same cases (ratio of times; less than 1 means the first is faster)
   :header-rows: 1
   :widths: 30 10 20 40

   * - Pair
     - Cases
     - Median ratio
     - First method faster in
   * - Get a Subset / kerchunk
     - 9
     - 0.42
     - 9 of 9 (ratios 0.29–0.61)
   * - Get a Subset / NCSS
     - 3
     - 0.33
     - 3 of 3 (ratios 0.23–0.40)
   * - kerchunk / NCSS
     - 4
     - 0.74
     - 3 of 4 (ratios 0.55–1.33)

Where the time went (cold successes): Get a Subset waited 2–6 minutes (median 4.2) for the remote processing
and downloaded at 8–11 MiB/s; kerchunk spent 74–82 % of its acquisition time in the range requests and the
rest decoding and writing; NCSS spent 84–96 % of its time acquiring the responses (its client does not separate
the server's preparation from the transfer). For all methods the common checks against the reference took
from under 1 minute (small domains) to 15–23 minutes (the largest domain); they are included in the times
above.

Data received, memory and requests
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. figure:: _static/remote_data_acquisition/stageB_transfer.png
   :alt: Bytes received divided by the float32 content V for each case and method.
   :width: 100%

   Bytes received (failed attempts included) divided by the float32 content V of each case.

.. list-table:: Totals over the cold successes
   :header-rows: 1
   :widths: 20 20 20 20 20

   * - Method
     - Received
     - HTTP requests
     - Retries
     - Peak RSS
   * - Get a Subset
     - 31.6 GiB (9 runs)
     - 25–62 per run (API + files)
     - 0
     - 0.43–0.58 GiB
   * - kerchunk
     - 125.9 GiB (10 runs)
     - 95–355 per run
     - 0
     - 0.70–0.74 GiB
   * - NCSS
     - 8.8 GiB (4 runs)
     - 15–90 per run
     - 19 (11 in one case)
     - 0.48–3.21 GiB

The bytes received by the two failed NCSS runs were not recorded by the client and are not included.

Scientific equivalence
~~~~~~~~~~~~~~~~~~~~~~

Three questions are kept apart: did the acquisition succeed, are the meteorological fields equivalent, and
are the energetics equivalent.

- **Between NCAR methods:** in every case where two or three methods succeeded, the five fields were identical
  value by value, and so were the energetics of Get a Subset and kerchunk; NCSS energetics differed from them by
  at most 2 × 10\ :sup:`-10` of the tolerance (the toolkit run inside the NCSS client).
- **Fields against the CDS:** every successful file passed the pre-registered tolerances. The largest
  differences were 2.4 × 10\ :sup:`-4` K in temperature and 0.06 m\ :sup:`2` s\ :sup:`-2` in geopotential;
  vertical velocity and both wind components were identical.
- **Energetics against the CDS:** the same verdict for every NCAR method in a case, because the fields are the
  same. **Only four of the nine cases with a CDS file are equivalent**; five are not (table below). The time
  derivative of the zonal available potential energy (∂Az/∂t) is outside the tolerance in all five, the
  residual RGz in four. The largest departure is about 0.07 % (6.7 × 10\ :sup:`-4`) of the column scale.

.. list-table:: Energetics against the CDS, column by column (tolerance: 10\ :sup:`-4` of the column scale)
   :header-rows: 1
   :widths: 16 18 20 46

   * - Case
     - Verdict
     - Largest difference / tolerance
     - Columns outside the tolerance
   * - 19790205
     - equivalent
     - 1.6 × 10\ :sup:`-5`
     - —
   * - 19790535
     - equivalent
     - 7.1 × 10\ :sup:`-5`
     - —
   * - 19790585
     - equivalent
     - 7.7 × 10\ :sup:`-5`
     - —
   * - 19800296
     - equivalent
     - 0.66
     - —
   * - 19791008
     - **not equivalent**
     - 1.24
     - ∂Az/∂t
   * - 19790644
     - **not equivalent**
     - 3.40
     - ∂Az/∂t, RGz
   * - 19800256
     - **not equivalent**
     - 3.48
     - BAz, ∂Az/∂t, ∂Ae/∂t, RGz, RGe
   * - 19790166
     - **not equivalent**
     - 5.18
     - Az, BAz, ∂Az/∂t, ∂Ae/∂t, RGz
   * - 19791150
     - **not equivalent**
     - 6.66
     - Az, Ae, Ca, BAz, Gz, ∂Az/∂t, ∂Ae/∂t, RGz, RGe
   * - 19800318
     - no CDS file
     - —
     - identical between the NCAR methods

**Case 19790644 (diagnostic, not a campaign result).** The energetics fail only in ∂Az/∂t (largest difference
0.0054 against a tolerance of 0.0016) and RGz (0.0059 against 0.0017). The GDEX temperature differs from the
CDS temperature at 2.3 % of the grid points, by at most 2.4 × 10\ :sup:`-4` K, only between 1000 and 800 hPa;
the geopotential at 2 % of the points, by at most 0.004 m\ :sup:`2` s\ :sup:`-2`. In a diagnostic run the GDEX
temperature was replaced by the CDS temperature: the energetics became equivalent to the CDS run (largest
difference 1.6 × 10\ :sup:`-5` of the tolerance). Replacing the geopotential instead left both failures
unchanged. So the departure comes from the small temperature differences, which the zonal available potential
energy and its time derivative amplify. These hybrid files were built only for this diagnosis and are not
results of the campaign.

**The other four cases outside the tolerance (hypothesis, not tested).** That the same mechanism explains
19791008, 19800256, 19790166 and 19791150 is a hypothesis: the substitution was made for 19790644 only. The
largest temperature difference alone does not predict the verdict: 19790205 and 19800296 are equivalent
although their temperatures differ from the CDS by up to 2.4 × 10\ :sup:`-4` and 2.1 × 10\ :sup:`-4` K, the
same order as in the failing cases; only 19790535 and 19790585 have identical temperatures.

**Interpretation and caveats.** The two sources hold the same ERA5 analysis, but the GDEX files are a separate
conversion of it, and their temperatures are not always bit-identical to those the CDS delivers. Differences far
below the field tolerances can still exceed the column tolerance of the energetics in tendency and residual
terms, which are small differences of large quantities. The tolerances were fixed before the campaign and were
not changed. Whether such departures matter depends on the use: they are far below the size of the energy terms
discussed in cyclone studies, but results from NCAR and CDS files should not be expected to agree to
10\ :sup:`-4` of each column's scale.

Limitations
~~~~~~~~~~~

- One run per case and method, on one day, from one place: no repetition in other windows.
- Ten cases from 1979–1980 over the South Atlantic; the Get a Subset volume depends on the ratio between the
  hours of the date range and the hours used (3-hour tracks here).
- NCSS ran with a 3-hour deadline per case instead of the 12 hours of the primary campaign, and was not run on
  four cases after its stop rule.
- The CDS was not part of Stage B; its times from the primary campaign (other day, other queue) are not compared
  here.
- The common checks of kerchunk and Get a Subset ran in the same process as their clients; NCSS also ran the
  energetics inside its client. Small differences in what each total includes are documented in the records.
- The kerchunk and Get a Subset times were measured with the benchmark clients. The integration under review in
  the toolkit reads in the same way (one request at a time) but has only been tested offline, against emulated
  services; its times have not been measured.

Choosing a method
-----------------

These are criteria, not a ranking.

.. list-table::
   :header-rows: 1
   :widths: 34 66

   * - If ...
     - consider ...
   * - the case is recent (later than the GDEX copy)
     - CDS API: GDEX publishes ERA5 months later.
   * - you need results that match a CDS-based analysis to 10\ :sup:`-4` of each column
     - CDS API: NCAR files hold the same fields within tight tolerances, but tendency and residual terms of
       the energetics can differ beyond that level (see above).
   * - you have no CDS account, or the CDS queue is long
     - NCAR NCSS (no account); Get a Subset (UCAR account) when available in the toolkit.
   * - the domain is large or the track is long
     - methods that cut on the server and do not depend on one long synchronous response per request: in Stage
       B Get a Subset completed all nine cases it ran, the largest included; NCSS failed on two of six (1.9 and
       3.7 GiB) and completed the largest (6.7 GiB). One run per case does not establish an effect of size.
   * - network volume is limited
     - methods that cut the area on the server (CDS, NCSS, Get a Subset); kerchunk transfers the whole globe.
   * - you cannot keep a personal token, or cannot wait for a remote request
     - NCSS or kerchunk (no account, no queue).

Examples
--------

**CDS API.**

.. code-block:: bash

   python lorenzcycletoolkit.py case.nc -t -r --cdsapi --trackfile inputs/track_case --time-resolution 3

**NCAR NCSS.** Show the requests first, then run:

.. code-block:: bash

   python -m src.utils.ncar_gdex case.nc --trackfile inputs/track_case --time-resolution 3 --dry-run
   python lorenzcycletoolkit.py case.nc -t -r --ncar --trackfile inputs/track_case --time-resolution 3

**NCAR kerchunk and Get a Subset** are not yet available in a release of the toolkit; examples will be added
when they are.

Reproducing the benchmark
-------------------------

The Stage B protocol and analysis rules were sealed (SHA-256) before the first request. The protocol, the raw
records (logs, ledgers, provenance records) and the audit are kept outside the repository, because they hold
account details. The table and figures of this page are produced by ``benchmarks/d_bench/stageb_report.py``
from the audit; the public table is ``_static/remote_data_acquisition/stageB_runs.csv`` (one row per planned
run, no private paths).
