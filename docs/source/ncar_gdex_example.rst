NCAR GDEX - Alternative to the CDS API (experimental)
=====================================================

``--ncar`` downloads the ERA5 pressure-level fields from the copy of ERA5 kept by the NCAR Geoscience Data
Exchange (GDEX, dataset `d633000 <https://doi.org/10.5065/BH6N-5N20>`_) instead of the Copernicus Climate Data
Store (CDS). GDEX cuts the area on its server and answers on the same connection. It needs no account and does
not go through the CDS job queue; it still takes some seconds to prepare each response, and it can be slow or
unavailable at times.

.. warning::

   This route is **experimental**. It has been validated on one two-day case: its values were bit-identical to
   an earlier download of the same GDEX data, and the toolkit reproduced that case's results byte for byte. It
   has **not** been compared with files downloaded from the CDS. Its options and messages may still change. The
   ``--cdsapi`` route and local NetCDF input are unchanged.

``--ncar`` or ``--cdsapi``?
---------------------------

Both request the same data for a track: the five fields (temperature, geopotential, vertical velocity, wind
components), the 37 pressure levels, the extent of the track with a 15-degree buffer rounded outwards to whole
degrees, and the same hours of each day for ``--time-resolution``.

.. list-table::
   :header-rows: 1
   :widths: 22 39 39

   * -
     - ``--cdsapi``
     - ``--ncar`` (experimental)
   * - Source
     - Copernicus CDS
     - NCAR GDEX (d633000)
   * - Account
     - CDS account and key
     - none
   * - Cost
     - free
     - free
   * - Dates
     - 1940 to a few days ago
     - 1940 to some months ago
   * - Waiting (1)
     - CDS job queue
     - no CDS job queue
   * - Requests (2)
     - one per day
     - one per field and day
   * - File (3)
     - CDS NetCDF file
     - NetCDF-4, float32
   * - Existing file
     - used as it is
     - used if verified
   * - Limits
     - those of the CDS
     - set by the toolkit

1. The length of the CDS queue varies. GDEX does not use it, but it still takes time to prepare each response:
   6 to 24 s per response in the validation, and the service can be slower or unavailable at times.
2. Five requests per day; two per field and day when the area crosses Greenwich.
3. The values of the responses unchanged (float32), losslessly compressed, with a provenance record next to the
   file. Without ``--ncar``, any existing NetCDF file is used as it is.

The values of GDEX and the CDS come from the same ECMWF archive, but GDEX converts the GRIB files itself; that
its values are bit-identical to the NetCDF files of the CDS has not been verified.

Dates available
---------------

GDEX publishes ERA5 some months after the CDS. On 2026-10-26 it covered 1940-01-01 to 2026-07-31; check the
current end date on the `dataset page <https://doi.org/10.5065/BH6N-5N20>`_. The last day of the track is
requested first, so a date that is not published yet fails on the first request, with a message that suggests
``--cdsapi``; nothing else is downloaded. For recent cases use ``--cdsapi``.

Costs
-----

- **Money:** none. No account or key is needed. Cite the dataset (see "Citation").
- **Volume:** about 4 bytes per value, so per day about
  5 fields × hours of the day × 37 levels × latitudes × longitudes × 4 bytes. ``--time-resolution 3`` gives 8
  hours per day, ``6`` gives 4. A two-day track at 6-hour resolution over 129 × 129 points (32° × 32°) needed 10
  requests of 9.9 MB each, 94 MiB in all.
- **Time:** in the validation (2026-10-27, from Brazil), each response waited 6 to 24 s before its first byte and
  then arrived at about 3 MB/s. One request is made at a time, so count 10 to 25 s per field and day of that size.
- **Disk:** the compressed file was 44 % of the volume downloaded (43.7 MB for 98.6 MB). During the download the
  responses are also kept in ``<file>.parts/``, removed at the end. The free disk space is checked before
  starting.

Before downloading: ``--dry-run``
---------------------------------

The standalone command shows the requests and their size without downloading anything:

.. code-block:: bash

   python -m src.utils.ncar_gdex case.nc --trackfile inputs/track_case \
       --time-resolution 6 --dry-run

For the two-day example above it printed:

.. code-block:: text

   area N-6 W-61 S-38 E-29: 129 x 129 points, 8 times, 37 levels
   10 requests, about 94.0 MiB to transfer, 94.0 MiB of values; names from ...
   https://tds.gdex.ucar.edu/thredds/ncss/grid/files/d633000/e5.oper.an.pl/201808/...
   (one line per request, cut here)

Usage
-----

.. code-block:: bash

   python lorenzcycletoolkit.py case.nc -t -r --ncar --trackfile inputs/track_case \
       --time-resolution 6

``case.nc`` does not need to exist: it is written by the download, then the energetics run on it. The track file
defines the area and the hours also in fixed mode (``-f --box_limits ...``), as with ``--cdsapi``. ``--ncar`` and
``--cdsapi`` cannot be used together.

The standalone command (the same as above without ``--dry-run``) downloads the file without running the
energetics. It also accepts:

- ``--workers``: up to 4 requests at a time;
- ``--max-gib``: largest volume received;
- ``--float64``: store the values as float64 (exact conversion);
- ``--no-compress``: write NetCDF-3 instead of compressed NetCDF-4;
- ``--keep-parts``: keep the responses after writing the file.

What is written
---------------

- ``case.nc``: the five fields, values copied unchanged (float32, no interpolation or unit conversion), with the
  variable and coordinate names of the namelist used by the run (``inputs/namelist`` or
  ``inputs/namelist.default``). The namelist must describe geopotential (m**2 s**-2), not geopotential height.
- ``case.nc.provenance.json``:

  - the URL, size and SHA-256 of every response, and when and by which version of the toolkit it was received;
  - the SHA-256 of ``case.nc`` and its numeric type;
  - how the request was derived from the track;
  - the software versions and the timings.
- ``lec_metadata.json`` of each run records the SHA-256 of its input file: a result is traced back to the
  provenance record by that SHA-256.

``case.nc`` appears only when every response has been received and checked, and the file itself has been checked.

Checks and limits
-----------------

Every response is checked before use:

- HTTP status, number of bytes and NetCDF signature;
- variable and units;
- times, levels, latitudes and longitudes exactly those requested;
- no missing or non-finite values.

Then:

- **Volume.** The default limit is 8 GiB (``--max-gib`` in the standalone command).

  - What is counted: every byte the toolkit reads from the service once decrypted (status lines, headers, chunk
    framing and bodies), failed attempts and retries included.
  - Before each request, the most its response may hold is reserved. The request is not made if that would pass
    the limit, and no connection reads more than its reservation. When the limit cannot hold every remaining
    request, this is said before the first request.
  - A response sent without Content-Length (chunked) also reserves room for its framing before its body is read.
  - When a response is abandoned before its end, bytes may have reached the computer without being read. They are
    not observed: the size of the connection's receive buffer is added to the count as an upper bound. If that
    takes the total past the limit, the program says so and makes no further request.
  - Not observed and not counted: the TLS handshake and encryption overhead, the answer of a proxy when one is
    used, and TCP/IP headers. The network traffic is therefore a little larger than the count.
- **Time.** An attempt fails after 300 s without data, or when it lasts longer than 300 s plus the time its
  largest response takes at 25 kB/s, also while bytes keep arriving slowly. Connecting (TCP and TLS handshake) is
  bounded by the 300 s without data of each step.
- **Retries.** Network failures, cut or too slow responses and HTTP 429/5xx are attempted up to 4 times in all,
  waiting 10 s, 30 s and 90 s between attempts, or the wait the service asks for (``Retry-After``) when it is
  longer. If the service asks to wait more than 300 s, the program stops and asks to run again later.
- **Stops at once.** The program stops at once on:

  - a refused request (HTTP 401, 403, 404 or another 4xx);
  - a redirect (never followed);
  - a page that is not NetCDF (e.g. maintenance);
  - a response larger than its reservation;
  - a response that is not the requested data.
- **One download per file.** While a file is downloaded, a lock file ``<file>.lock`` sits next to it. A second
  download of the same file (for instance two runs started together) stops at once and says so. Runs that only
  read an existing, verified file are not affected. If a download is killed, its lock does not block the next
  one (on Windows, where this kind of lock is not available, the message says how to remove a lock file left
  behind).

Running again, and avoiding unwanted downloads
----------------------------------------------

- **Existing file.** An existing ``case.nc`` is never overwritten.

  - It is reused only if all of these hold:

    - its provenance record matches its SHA-256;
    - it was written with the same variable names and numeric type;
    - it covers every time, level, latitude and longitude of the request (a file covering more is reused).

  - Otherwise the program stops and explains why.
  - Run without ``--ncar`` to use any NetCDF file as a plain local input.
- **Interrupted download.** The verified responses are kept in ``case.nc.parts/``: running the same command again
  downloads only what is missing. A kept response is reused only if its SHA-256 and content still match.
- **Misspelt name.** A misspelt input file name is a new file and starts a new download: check the name, or run
  ``--dry-run`` first.

Precision: float32 and float64
------------------------------

The toolkit computes in the numeric type of its input file. GDEX delivers float32, and ``--ncar`` keeps it: this
is the default.

For the same ERA5 values stored as float32 and as float64 (two-day case, 129 x 129 points, 37 levels), the results
of the fixed and moving-box runs differed by at most 2.7e-7 of the largest value of each output column, with no
change of sign. The values of the mass term ``M`` were identical. The regression references of the test suite were
computed from float64 input, so float32 input does not reproduce them to their tolerance (1e-10 relative).

``python -m src.utils.ncar_gdex ... --float64`` stores the values as float64 (exact conversion). In the case
above this cost:

- 3 % more disk for the file;
- twice the memory for the input;
- in the moving-box run, about twice the time and 40 % more memory;
- in the fixed-box run, about 5 % more of each.

A float64 file is not reused by a run that asks for float32, nor the reverse.

Citation
--------

Cite the dataset: ECMWF (2019), ERA5 Reanalysis (0.25 Degree Latitude-Longitude Grid), Research Data Archive at
the National Center for Atmospheric Research, Computational and Information Systems Laboratory,
https://doi.org/10.5065/BH6N-5N20, together with the Copernicus terms of use of ERA5 (licence CC-BY-4.0).
