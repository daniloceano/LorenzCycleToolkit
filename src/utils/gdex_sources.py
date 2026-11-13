# **************************************************************************** #
#                                                                              #
#                                                         :::      ::::::::    #
#    gdex_sources.py                                    :+:      :+:    :+:    #
#                                                     +:+ +:+         +:+      #
#    By: daniloceano <danilo.oceano@gmail.com>      +#+  +:+       +#+         #
#                                                 +#+#+#+#+#+   +#+            #
#    Created: 2026/10/31 20:00:00 by daniloceano       #+#    #+#              #
#    Updated: 2026/10/31 20:00:00 by daniloceano      ###   ########.fr        #
#                                                                              #
# **************************************************************************** #

"""
Other ways of reading the NCAR GDEX ERA5 files (experimental).

:mod:`src.utils.ncar_gdex` retrieves the ERA5 fields of a track through the
THREDDS NetCDF Subset Service (NCSS). This module adds two other access
methods to the same retrieval, for the same dataset (d633000) and the same
plan (area, 37 levels, hours of the CDS route):

``kerchunk``
    One HTTPS byte-range GET per field and hour of the original GDEX NetCDF-4
    files, located with the GDEX kerchunk reference files (JSON) and decoded
    here (zlib + shuffle, float32). Each GET returns the global chunk
    (37 levels x 721 x 1440 points, about 55-60 MiB compressed); the area is
    cut on this side, so the transfer does not shrink with the domain.
    A box crossing Greenwich has two longitude segments (two pieces per
    field and day) that share every chunk: each chunk is read once and
    both pieces are cut from it.
    The reference files must be given as a local directory
    (``--kerchunk-refs``): this version does not download them.

``subset``
    The GDEX "Get a Subset" request API (needs a personal GDEX API token):
    one request for the whole area and period is submitted, processed by
    GDEX, listed, and its files (one per field and day, every hour of the
    period) are downloaded and cut to the hours of the plan here. A box
    crossing Greenwich is never submitted as such (GDEX did not cut it
    correctly in Frente 2b): it is one request per side, with the limits of
    the plan's two longitude segments.

Both write, for every request of the plan, the same verified piece an NCSS
response gives (one field, one day, one longitude segment, the values
unchanged) into ``<file>.parts/``; :func:`src.utils.ncar_gdex.retrieve`
then checks each piece (:func:`~src.utils.ncar_gdex.check_piece`), resumes
from the pieces already verified, assembles the file piece by piece, checks
it and writes the provenance record, exactly as for NCSS. The three methods
read the same GDEX values: in the D-Bench Stage B campaign (2026-10-31) their
files were identical value by value.

The bytes received are counted against the retrieval's
:class:`~src.utils.ncar_gdex.ByteBudget` (every attempt, failed ones
included); with the toolkit's opener the counting connections apply.

Get a Subset rules, from the D-Bench campaigns:

* GDEX deduplicates identical requests while the earlier one is listed on the
  account, in any state ("Set for Purge" included). Before submitting, the
  account's requests are listed: a completed identical request is reused (and
  recorded as reused, never purged); one waiting for purge, or failed, stops
  the retrieval with an explanation. Nothing is ever resubmitted with altered
  parameters to get around the deduplication.
* An account holds at most 8 open requests: nothing is submitted when the
  requests listed plus those to submit would pass ``SUBSET_MAX_LISTED``.
* The API redirects some addresses (since 2026-11: ``status/<id>`` and the
  like, to the same path with a final slash). A redirect is followed, up to
  ``SUBSET_MAX_REDIRECTS`` times, only to the API's own scheme and host, so
  that the token is never sent elsewhere; a submission is never redirected.
* A submission may take minutes and create the request even when the client
  gave up waiting: the client waits ``SUBSET_SUBMIT_TIMEOUT`` seconds and,
  without an answer, looks for a new listed request that matches; it adopts it
  instead of submitting again (and never purges it: its origin is not proven).
* The request is recorded in ``<file>.parts/subset_request.json`` (the
  second side of a box crossing Greenwich in
  ``subset_request.segment1.json``) as soon as its number is known, so an
  interrupted retrieval resumes with the same request and the files already
  verified (SHA-256), never with a new one.
* Purge (``--gdex-purge``, off by default) only removes the requests this
  retrieval created (their number came back in the answer to the submission,
  recorded in the request state files), after the output file and its provenance record are
  written; for each, the service's answer and the request's state in the listing that
  follows are recorded. "Set for Purge" means the purge is pending, not done.

The token is read from the environment variable ``GDEX_API_TOKEN`` or from
the file ``~/.gdex_api_token`` (readable by its owner only). It is sent as
the query parameter the API expects and never written to logs, errors or
records.

Created by:
    Danilo Couto de Souza
    Universidade de São Paulo (USP)
    Instituto de Astronomia, Geofísica e Ciências Atmosféricas
    São Paulo - Brazil

Contact:
    danilo.oceano@gmail.com
"""

import datetime
import hashlib
import json
import os
import re
import stat
import threading
import time
import urllib.error
import urllib.parse
import urllib.request
import zlib

import netCDF4
import numpy as np

from . import ncar_gdex as g

# Coordinate attributes of the pieces: those of the NCSS responses, so that the
# file written from any method carries the same metadata.
PIECE_COORDINATE_ATTRS = {
    "time": {"long_name": "time", "units": g.TIME_UNITS, "calendar": "gregorian"},
    "level": {"long_name": "pressure level", "short_name": "plev", "units": "hPa", "alternate_units": "millibar"},
    "latitude": {"long_name": "latitude", "short_name": "lat", "units": "degrees_north", "qualifier": "Gaussian",
                 "standard_name": "latitude"},
    "longitude": {"long_name": "longitude", "short_name": "lon", "units": "degrees_east",
                  "standard_name": "longitude"},
}
# Field attributes never copied into a piece (storage details of the source).
_SKIP_ATTRS = {"_FillValue", "missing_value", "coordinates", "_ARRAY_DIMENSIONS",
               "QuantizeGranularBitGroomNumberOfSignificantDigits", "_QuantizeGranularBitGroomNumberOfSignificantDigits"}

# kerchunk
KERCHUNK_NAMES = {"T": "temperature", "Z": "geopotential", "W": "vertical_velocity",
                  "U": "u_component_of_wind", "V": "v_component_of_wind"}
# The references name data.rda.ucar.edu, whose TLS certificate failed on 2026-10-30; the same
# files are served by data.gdex.ucar.edu.
KERCHUNK_HOST_REWRITE = ("https://data.rda.ucar.edu/", "https://data.gdex.ucar.edu/")
KERCHUNK_LEVELS = 37
KERCHUNK_GRID_STEP = 0.25

# Get a Subset
SUBSET_API = "https://gdex.ucar.edu/api/"
SUBSET_GROUP_INDEX = "26"            # ERA5 atmospheric pressure level analysis [netCDF4]
SUBSET_MAX_LISTED = 7                # nothing is submitted with this many requests listed (GDEX limit: 8)
SUBSET_SUBMIT_TIMEOUT = 600.0
SUBSET_API_TIMEOUT = 120.0
SUBSET_API_LIMIT = 4 * 1024 ** 2     # bytes reserved for one API answer
SUBSET_REMOTE_DEADLINE = 6 * 3600.0
SUBSET_POLL = (30.0, 120.0, 900.0)   # seconds between status checks: first, later, after this long
SUBSET_COMPRESSION = 0.5             # files / float32 content of every hour, upper bound of the D-Bench ratios (0.43)
SUBSET_STATE = "subset_request.json"
SUBSET_MAX_REDIRECTS = 3             # redirects followed for one API call (same scheme and host only)
TOKEN_ENV = "GDEX_API_TOKEN"
TOKEN_FILE = os.path.join("~", ".gdex_api_token")


class _Ticket:
    """What a :class:`~src.utils.ncar_gdex.ByteBudget` needs to know of one attempt."""

    def __init__(self, key, url, max_bytes):
        self.key, self.url, self.max_bytes = key, url, int(max_bytes)


# ------------------------------------------------------------------ pieces
def write_piece(path, request, values, attributes):
    """
    Write one piece shaped like the NCSS response of ``request``.

    Args:
        path (str): File to write (NetCDF-3 64-bit offset).
        request (PieceRequest): The request the piece answers.
        values (np.ndarray): float32 values, time x level x latitude x longitude.
        attributes (dict): Attributes of the source variable (units, names...).
    """
    values = np.asarray(values)
    if values.dtype != np.float32 or values.shape != request.shape:
        raise g.NcarContentError(f"{request.key}: values of type {values.dtype} and shape {values.shape} "
                                 f"instead of float32 {request.shape}")
    fill = attributes.get("_FillValue", attributes.get("missing_value"))
    with g._NETCDF_LOCK, netCDF4.Dataset(path, "w", format="NETCDF3_64BIT_OFFSET") as out:
        for name, size in zip(g.DIMENSIONS, request.shape):
            out.createDimension(name, size)
        for name, dtype, data in (
            ("time", "i4", g._hours_since_origin(g._moments(request.day, request.hours))),
            ("level", "f8", g.PRESSURE_LEVELS_HPA),
            ("latitude", "f8", request.latitudes),
            ("longitude", "f8", request.segment.response),
        ):
            variable = out.createVariable(name, dtype, (name,))
            variable.setncatts(PIECE_COORDINATE_ATTRS[name])
            variable[:] = np.asarray(data)
        variable = out.createVariable(request.field.name, "f4", g.DIMENSIONS,
                                      fill_value=None if fill is None else np.float32(fill))
        variable.setncatts({k: v for k, v in attributes.items() if k not in _SKIP_ATTRS})
        if fill is not None:
            variable.setncattr("missing_value", np.float32(fill))
        variable.set_auto_maskandscale(False)
        variable[:] = values


def _keep_piece(request, parts_dir, values, attributes, http, attempts, elapsed, started_utc):
    """Write, check and store one piece; returns its PieceResult."""
    path = os.path.join(parts_dir, f"{request.key}.nc")
    part = f"{path}.part"
    try:
        write_piece(part, request, values, attributes)
        g.check_piece(part, request)
    except g._Retryable as error:
        raise g.NcarContentError(f"{request.key}: the piece written here cannot be read ({error})") from None
    except BaseException:
        if os.path.exists(part):
            os.remove(part)
        raise
    os.replace(part, path)
    sha256, size = g._sha256(path), os.path.getsize(path)
    received = g._iso(started_utc)
    # "url" identifies the content of the piece (the NCSS request it is equivalent to); "http" says
    # where its values came from
    g._write_json(f"{path}.json", {"url": request.url, "sha256": sha256, "bytes": size, "http": http,
                                   "received_utc": received, "lectk_commit": g._lectk_commit()})
    return g.PieceResult(request, path, sha256, size, attempts, elapsed, http, False, received, g._lectk_commit())


# ------------------------------------------------------------- transfers
def _redact(text, secret):
    text = str(text)
    if secret:
        for form in {secret, urllib.parse.quote(secret, safe=""), urllib.parse.quote_plus(secret)}:
            text = text.replace(form, "<token>")
    return text


REDIRECT_STATUS = (301, 302, 303, 307, 308)


class _Redirect(g.NcarRequestError):
    """An answer redirecting the request (raised by :func:`_transfer` when asked to report redirects)."""

    def __init__(self, message, code, location):
        super().__init__(message)
        self.code, self.location = code, location


def _transfer(budget, ticket, url, opener, timeout, cancel, sink=None, headers=None, method="GET", data=None,
              expect=(200,), secret=None, redirects=False):
    """
    One HTTP attempt, reserved in ``budget`` and charged what it received.

    The body goes to the open file ``sink``, or is returned. Bytes counted by
    the toolkit's connections are charged; with another opener a complete
    answer is charged its headers and body, any other attempt its whole
    reservation.

    The opener never follows redirects. With ``redirects`` true, a redirect
    is raised as :class:`_Redirect` with its ``Location`` so that the caller
    can decide whether to follow it; otherwise it is a refused request.

    Returns:
        tuple: (body or None, size, sha256, http dict).

    Raises:
        g._Retryable (network failure, 429/5xx, incomplete answer),
        g.NcarRequestError (other HTTP status; _Redirect for a redirect when
        ``redirects`` is true), g._Cancelled.
    """
    if cancel.is_set():
        raise g._Cancelled()
    budget.reserve(ticket)
    message = urllib.request.Request(url, data=data, method=method,
                                     headers={"User-Agent": g.USER_AGENT, **(headers or {})})
    message.lectk_limit = ticket.max_bytes
    started = time.monotonic()
    deadline = started + timeout + ticket.max_bytes / g.MIN_TRANSFER_RATE
    message.lectk_deadline = deadline
    received, header_bytes, size, complete, status = None, 0, 0, False, None
    digest, chunks = hashlib.sha256(), []
    result, failure = None, None
    try:
        try:
            response = opener(message, timeout=timeout)
        except urllib.error.HTTPError as error:
            received = getattr(getattr(error, "fp", None), "lectk_received", None)
            status = error.code
            body = _redact(g._error_text(error), secret)
            location = error.headers.get("Location") if error.headers is not None else None
            error.close()
            if redirects and error.code in REDIRECT_STATUS:
                raise _Redirect(f"{ticket.key}: HTTP {error.code} {error.reason} (redirect)", error.code, location)
            if error.code in g.RETRY_STATUS:
                raise g._Retryable(f"HTTP {error.code} {error.reason}: {body or 'no message'}",
                                   g._retry_after(error.headers))
            raise g.NcarRequestError(f"{ticket.key}: HTTP {error.code} {error.reason}: {body or 'no message'}")
        received = getattr(response, "lectk_received", None)
        with response:
            status = getattr(response, "status", None)
            header_bytes = g._header_bytes(response)
            if status not in expect:
                raise g.NcarRequestError(f"{ticket.key}: HTTP {status}, expected {' or '.join(map(str, expect))}")
            length = response.headers.get("Content-Length")
            expected = int(length) if length is not None and length.isdigit() else None
            limit = ticket.max_bytes - header_bytes - 1
            if expected is not None and expected > limit:
                raise g.NcarContentError(f"{ticket.key}: the service announced {expected} bytes, more than the "
                                         f"{limit} this request can hold")
            while True:
                if cancel.is_set():
                    raise g._Cancelled()
                if time.monotonic() > deadline:
                    raise g._Retryable(f"response too slow: {size} bytes in {time.monotonic() - started:.0f} s")
                block = response.read(min(1 << 20, limit + 1 - size))
                if not block:
                    break
                size += len(block)
                if size > limit:
                    raise g.NcarContentError(f"{ticket.key}: answer larger than the {limit} bytes it can hold")
                digest.update(block)
                if sink is None:
                    chunks.append(block)
                else:
                    sink.write(block)
            complete = expected is None or size == expected
            if not complete:
                raise g._Retryable(f"incomplete answer: {size} of {expected} bytes")
        result = (b"".join(chunks) if sink is None else None, size, digest.hexdigest(),
                  {"status": status, "header_bytes": header_bytes, "seconds": round(time.monotonic() - started, 3)})
    except (g._Retryable, g.NcarRetrievalError, g._Cancelled) as error:
        failure = error
    except g.NETWORK_ERRORS as error:
        failure = g._Retryable(_redact(f"{type(error).__name__}: {error}", secret))
    except Exception as error:   # the toolkit's connections signal the byte or time limit with their own errors
        failure = g._Retryable(_redact(f"{type(error).__name__}: {error}", secret))
    if received is not None:
        charged = received.bytes + (0 if complete else (received.window or g.WINDOW_ALLOWANCE))
    elif failure is None:
        charged = header_bytes + size
    else:
        charged = ticket.max_bytes
    budget.settle(ticket, charged, {"status": status, "seconds": round(time.monotonic() - started, 3),
                                    "result": "ok" if failure is None else _redact(g._describe(failure), secret)})
    if failure is not None:
        raise failure
    return result


def _retrying(what, attempts, backoff, cancel, logger, call):
    """Call ``call(attempt)`` until it succeeds, retrying g._Retryable failures with growing waits."""
    for attempt in range(1, attempts + 1):
        try:
            return call(attempt), attempt
        except g._Retryable as error:
            if attempt == attempts:
                raise g.NcarTransferError(f"{what}: failed {attempts} times; last error: {error}") from None
            asked = error.retry_after or 0.0
            if asked > g.MAX_WAIT:
                raise g.NcarTransferError(f"{what}: the service asked to wait {asked:.0f} s; run again later "
                                          f"(the verified pieces are kept)") from None
            delay = min(g.MAX_WAIT, max(backoff * 3 ** (attempt - 1), asked))
            logger.warning(f"⚠️ {what}: {error}; retrying in {delay:.0f} s (attempt {attempt + 1}/{attempts})")
            g._wait(delay, cancel)
    raise AssertionError("unreachable")


class _Source:
    """Common part of the access methods (see :func:`src.utils.ncar_gdex.retrieve`, ``source``)."""

    name = None
    service = None

    def __init__(self, opener=None, timeout=None, attempts=None, backoff=None):
        self.opener = opener or g.default_opener()
        self.timeout = g.DEFAULT_TIMEOUT if timeout is None else timeout
        self.attempts = g.DEFAULT_ATTEMPTS if attempts is None else attempts
        self.backoff = backoff
        self.details = {}

    @property
    def _backoff(self):
        return g.DEFAULT_BACKOFF if self.backoff is None else self.backoff

    def estimated_bytes(self, plan):
        raise NotImplementedError

    def disk_bytes(self, plan, estimate):
        """Disk space the method needs besides the output file (default: the pieces)."""
        return plan.estimated_bytes

    def prepare(self, plan, parts_dir, budget, logger, cancel):
        """Work done once before the pieces, also when every piece was kept (default: nothing)."""

    def fetch(self, request, parts_dir, budget, logger, cancel):
        raise NotImplementedError

    def finish(self, record, logger):
        """Called after the file and its provenance record are in place; returns what to add to the record."""
        return {}

    def describe(self):
        return dict(self.details, method=self.name, service=self.service)


# -------------------------------------------------------------- kerchunk
class KerchunkSource(_Source):
    """
    Byte ranges of the GDEX NetCDF-4 files, located by kerchunk references.

    Args:
        refs_dir (str): Directory with the GDEX reference files
            ``<variable>.<decade start>-<decade end>-remote.json`` (for example
            ``temperature.1970-1979-remote.json``).
        opener, timeout, attempts, backoff: As for NCSS
            (:func:`src.utils.ncar_gdex.retrieve`).
    """

    name = "kerchunk"
    service = "HTTPS byte ranges of the GDEX NetCDF-4 files located with GDEX kerchunk references"

    def __init__(self, refs_dir, **options):
        super().__init__(**options)
        if not refs_dir or not os.path.isdir(refs_dir):
            raise ValueError(f"kerchunk reference directory not found: {refs_dir!r}")
        self.refs_dir = refs_dir
        self._refs = {}
        self._chunks = None
        self._siblings = {}     # (field name, day) -> the requests of that field and day, one per segment
        self._cut = {}          # request key -> values cut for it from chunks read for another segment

    def _references(self, field_name, year):
        decade = year - year % 10
        key = (field_name, decade)
        if key not in self._refs:
            name = f"{KERCHUNK_NAMES[field_name]}.{decade}-{decade + 9}-remote.json"
            path = os.path.join(self.refs_dir, name)
            if not os.path.isfile(path):
                raise g.NcarRetrievalError(f"kerchunk reference file {name} not found in {self.refs_dir}")
            with open(path, encoding="utf-8") as handle:
                refs = json.load(handle)["refs"]
            array = refs[f"{field_name}/.zarray"]
            array = json.loads(array) if isinstance(array, str) else array
            attrs = refs.get(f"{field_name}/.zattrs") or {}
            attrs = json.loads(attrs) if isinstance(attrs, str) else attrs
            chunks = array.get("chunks") or []
            if (array.get("dtype") != "<f4" or (array.get("compressor") or {}).get("id") != "zlib"
                    or [f.get("id") for f in array.get("filters") or []] != ["shuffle"]
                    or len(chunks) != 4 or chunks[0] != 1 or chunks[1] != KERCHUNK_LEVELS
                    or chunks[3] != int(round(360 / KERCHUNK_GRID_STEP))):
                raise g.NcarContentError(f"{name}: unexpected storage of {field_name} ({array}); this version reads "
                                         f"float32 chunks of one hour, {KERCHUNK_LEVELS} levels and the global "
                                         f"0.25-degree grid, zlib + shuffle")
            self._refs[key] = {"refs": refs, "chunks": chunks, "attrs": attrs, "decade": decade,
                               "file": name, "sha256": g._sha256(path)}
            self.details.setdefault("reference_files", {})[name] = self._refs[key]["sha256"]
        return self._refs[key]

    def _chunk(self, field, moment):
        ref = self._references(field.name, moment.year)
        index = (moment - datetime.datetime(ref["decade"], 1, 1)) // datetime.timedelta(hours=1)
        try:
            url, offset, size = ref["refs"][f"{field.name}/{index}.0.0.0"]
        except (KeyError, ValueError):
            raise g.NcarRetrievalError(f"{ref['file']} has no chunk for {field.name} at {moment:%Y-%m-%d %H} UTC "
                                       f"(GDEX may not have published this date yet; --cdsapi may)") from None
        if f".{field.code}.{field.grid}.{moment:%Y%m%d}00_" not in url:
            raise g.NcarContentError(f"{ref['file']}: chunk {index} points to {url}, not to the {field.name} "
                                     f"file of {moment:%Y-%m-%d}")
        if url.startswith(KERCHUNK_HOST_REWRITE[0]):
            url = KERCHUNK_HOST_REWRITE[1] + url[len(KERCHUNK_HOST_REWRITE[0]):]
        return url, int(offset), int(size), ref

    def _plan_chunks(self, plan):
        hours = {}
        for request in plan.requests:
            for moment in g._moments(request.day, request.hours):
                hours[(request.field.name, moment)] = request.field
        return {key: self._chunk(field, key[1]) for key, field in hours.items()}

    def estimated_bytes(self, plan):
        self._chunks = self._plan_chunks(plan)
        return sum(size + g.HEADER_ALLOWANCE for _, _, size, _ in self._chunks.values())

    def prepare(self, plan, parts_dir, budget, logger, cancel):
        if self._chunks is None:
            self._chunks = self._plan_chunks(plan)
        self.details["chunks"] = len(self._chunks)
        self.details["chunk_bytes"] = sum(size for _, _, size, _ in self._chunks.values())
        self.details["host_rewrite"] = list(KERCHUNK_HOST_REWRITE)
        self._siblings, self._cut = {}, {}
        for request in plan.requests:
            self._siblings.setdefault((request.field.name, request.day), []).append(request)

    @staticmethod
    def _columns(request):
        return np.rint((np.asarray(request.segment.response) % 360.0) / KERCHUNK_GRID_STEP).astype(int) % 1440

    def fetch(self, request, parts_dir, budget, logger, cancel):
        started_utc = datetime.datetime.now(datetime.timezone.utc)
        started = time.monotonic()
        kept = self._cut.pop(request.key, None)
        if kept is not None:
            # a box crossing Greenwich has two longitude segments that share every global chunk: their values were
            # cut when the other segment's piece read the chunks, which are never read twice
            values, origins, attrs, reader = kept
            http = {"method": self.name, "status": 206, "chunks": origins, "chunks_read_with": reader}
            return _keep_piece(request, parts_dir, values, attrs, http, 1, time.monotonic() - started, started_utc)
        # the other segments of the same field and day still to obtain are cut from the same chunks
        others = [r for r in self._siblings.get((request.field.name, request.day), ())
                  if r.key != request.key and r.key not in self._cut and g._stored_piece(r, parts_dir) is None]
        targets = [request] + others
        rows = np.rint((90.0 - np.asarray(request.latitudes)) / KERCHUNK_GRID_STEP).astype(int)
        columns = [self._columns(r) for r in targets]
        cuts = [np.empty(r.shape, dtype=np.float32) for r in targets]
        values = cuts[0]
        origins, attempts, attrs = [], 1, None   # attempts: 1 + the GETs repeated for any hour
        for k, moment in enumerate(g._moments(request.day, request.hours)):
            url, offset, size, ref = self._chunks[(request.field.name, moment)]
            if rows.min() < 0 or rows.max() >= ref["chunks"][2]:
                raise g.NcarContentError(f"{request.key}: latitudes outside the chunks of {ref['file']}")
            ticket = _Ticket(f"{request.key}_{moment:%H}", f"{url}#bytes={offset}-{offset + size - 1}",
                             size + g.HEADER_ALLOWANCE)

            def get(attempt, url=url, offset=offset, size=size, ticket=ticket):
                body, got, sha, http = _transfer(budget, ticket, url, self.opener, self.timeout, cancel,
                                                 headers={"Range": f"bytes={offset}-{offset + size - 1}"},
                                                 expect=(206,))
                if got != size:
                    raise g._Retryable(f"range of {got} bytes instead of {size}")
                return body, sha, http

            (body, sha, http), used = _retrying(ticket.key, self.attempts, self._backoff, cancel, logger, get)
            attempts += used - 1
            try:
                raw = np.frombuffer(zlib.decompress(body), dtype=np.uint8)
                shape = (KERCHUNK_LEVELS, ref["chunks"][2], ref["chunks"][3])
                chunk = raw.reshape(4, -1).T.copy().view("<f4").reshape(shape)
            except (zlib.error, ValueError) as error:
                raise g.NcarContentError(f"{ticket.key}: chunk cannot be decoded ({error})") from None
            band = chunk[:, rows]
            for cut, cols in zip(cuts, columns):
                cut[k] = band[:, :, cols]
            origins.append({"url": url, "offset": offset, "bytes": size, "sha256": sha, "status": http["status"]})
            attrs = ref["attrs"]
        for other, cut in zip(others, cuts[1:]):
            self._cut[other.key] = (cut, origins, attrs, request.key)
        http = {"method": self.name, "status": 206, "chunks": origins}
        if others:
            http["chunks_also_used_by"] = [r.key for r in others]
        return _keep_piece(request, parts_dir, values, attrs, http, attempts, time.monotonic() - started,
                           started_utc)


# ------------------------------------------------------------ Get a Subset
def read_token(environ=None, path=None):
    """
    The personal GDEX API token, from GDEX_API_TOKEN or ~/.gdex_api_token.

    Raises:
        g.NcarRequestError: Without a token, or if the file can be read by
            other users (its permissions are reported, never its content).
    """
    environ = os.environ if environ is None else environ
    if environ.get(TOKEN_ENV, "").strip():
        return environ[TOKEN_ENV].strip()
    path = os.path.expanduser(path or TOKEN_FILE)
    if not os.path.isfile(path):
        raise g.NcarRequestError(f"Get a Subset needs a personal GDEX API token: set {TOKEN_ENV} or write it to "
                                 f"{path} (readable by you only)")
    mode = stat.S_IMODE(os.stat(path).st_mode)
    if mode & 0o077:
        raise g.NcarRequestError(f"{path} can be read by other users (mode {mode:o}); run chmod 600 {path}")
    with open(path, encoding="utf-8") as handle:
        token = handle.read().strip()
    if not token:
        raise g.NcarRequestError(f"{path} is empty")
    return token


def subset_control(plan):
    """The Get a Subset request of a plan: all five fields, every hour from its first to its last time."""
    return {
        "dataset": g.SOURCE["dataset"],
        "groupindex": SUBSET_GROUP_INDEX,
        "date": f"{plan.times[0]:%Y%m%d%H%M}/to/{plan.times[-1]:%Y%m%d%H%M}",
        "param": "/".join(f.name for f in g.FIELDS),
        "product": "Analysis",
        "oformat": "netCDF",
        "nlat": g._number(plan.north),
        "slat": g._number(plan.south),
        "wlon": g._number(plan.west),
        "elon": g._number(plan.east),
    }


def _plan_segments(plan):
    """The longitude segments of a plan, in their order in the output file."""
    segments = {}
    for request in plan.requests:
        segments.setdefault(request.segment_index, request.segment)
    return [segments[index] for index in sorted(segments)]


def subset_controls(plan):
    """
    The Get a Subset requests of a plan: one per longitude segment.

    A box that does not cross Greenwich is one request,
    :func:`subset_control` unchanged. A box that crosses it is one request
    per side, with the limits of the plan's longitude segments (for example
    -30 to -0.25 and 0 to 20): GDEX does not cut a box crossing 0 degrees
    correctly (Frente 2b, 2026-11), and the pieces are per segment anyway.
    """
    control = subset_control(plan)
    segments = _plan_segments(plan)
    if len(segments) == 1:
        return [control]
    return [dict(control, wlon=g._number(s.west), elon=g._number(s.east)) for s in segments]


def _state_name(segment_index):
    """Request state file of a segment (the first keeps the name of a single request)."""
    return SUBSET_STATE if segment_index == 0 else f"subset_request.segment{segment_index}.json"


def rinfo_matches(rinfo, control):
    """
    True if the request summary GDEX lists (``rinfo``) describes exactly ``control``.

    Dataset, period, product group and area must be listed and equal; the
    fields (``parameters=8!d633000:T,...``), levels (empty: all) and output
    format are compared when the summary lists them.
    """
    pairs = dict(item.split("=", 1) for item in (rinfo or "").split(";") if "=" in item)
    start, end = control["date"].split("/to/")

    def moment(text):
        return f"{text[:4]}-{text[4:6]}-{text[6:8]} {text[8:10]}:{text[10:12]}"

    try:
        if not (pairs.get("dsnum") == control["dataset"] and pairs.get("startdate") == moment(start)
                and pairs.get("enddate") == moment(end) and pairs.get("tindex") == control["groupindex"]
                and all(float(pairs[k]) == float(control[k]) for k in ("nlat", "slat", "wlon", "elon"))):
            return False
    except (KeyError, ValueError):
        return False
    if "parameters" in pairs and ({p.rsplit(":", 1)[-1] for p in pairs["parameters"].split(",")}
                                  != set(control["param"].split("/"))):
        return False
    if pairs.get("level", ""):
        return False
    return "ofmt" not in pairs or pairs["ofmt"].lower() == control["oformat"].lower()


_FILE_NAME = re.compile(r"\.(?P<code>\d+_\d+_\w)\.(?P<grid>ll025\w+)\.(?P<start>\d{10})_(?P<end>\d{10})\.nc$")


class _SubsetRequest:
    """One Get a Subset request of a retrieval (one per longitude segment)."""

    def __init__(self, index, control, single):
        self.index, self.control, self.single = index, control, single
        self.state = None
        self.details = {}

    @property
    def state_file(self):
        return _state_name(self.index)

    @property
    def label(self):
        """How the logs name the request ("" for the only one)."""
        return "" if self.single else f" (longitudes {self.control['wlon']} to {self.control['elon']})"


class SubsetSource(_Source):
    """
    One GDEX Get a Subset request per longitude segment of the plan.

    A box that does not cross Greenwich is one request; one that crosses it
    is two, one per side (:func:`subset_controls`).

    Args:
        token (str): Personal GDEX API token (default: :func:`read_token`).
        purge (bool): Purge the requests after the file is written, only
            those this retrieval created.
        api (str): API address.
        clock (callable): time.monotonic replacement (tests).
        opener, timeout, attempts, backoff: As for NCSS.
    """

    name = "subset"
    service = "GDEX Get a Subset request API"

    def __init__(self, token=None, purge=False, api=SUBSET_API, clock=None, **options):
        super().__init__(**options)
        self._token = token
        self.purge = bool(purge)
        self.api = api
        self.clock = clock or time.monotonic
        self.requests = []
        self.files = {}

    @property
    def token(self):
        if self._token is None:
            self._token = read_token()
        return self._token

    @property
    def state(self):
        """State of the first (for a box not crossing Greenwich, the only) request."""
        return self.requests[0].state if self.requests else None

    def estimated_bytes(self, plan):
        hours = int((plan.times[-1] - plan.times[0]) // datetime.timedelta(hours=1)) + 1
        return int(SUBSET_COMPRESSION * 4 * len(g.FIELDS) * hours * len(g.PRESSURE_LEVELS_HPA)
                   * len(plan.latitudes) * len(plan.longitudes))

    def disk_bytes(self, plan, estimate):
        # the downloaded files are kept in the parts directory next to the pieces
        return estimate + plan.estimated_bytes

    # --- API
    def _with_token(self, url):
        """``url`` with this account's token as its only ``token`` parameter."""
        parsed = urllib.parse.urlsplit(url)
        query = [(k, v) for k, v in urllib.parse.parse_qsl(parsed.query, keep_blank_values=True) if k != "token"]
        query.append(("token", self.token))
        return urllib.parse.urlunsplit(parsed._replace(query=urllib.parse.urlencode(query)))

    def _redirect_target(self, moved, current, method, shown):
        """
        The address a redirect of the API leads to, if it may be followed.

        Only redirects to the API's own scheme and host (and port) are
        followed, so that the token is never sent elsewhere; a submission is
        never redirected (it is never repeated, see the module notes).

        Returns:
            tuple: (URL with the token, URL shown in logs and records).
        """
        if method == "POST":
            raise g.NcarRequestError(f"GDEX API {method} {shown}: the API redirected the submission (HTTP "
                                     f"{moved.code}); a submission is never redirected nor repeated, nothing more "
                                     f"was submitted")
        if not moved.location:
            raise g.NcarRequestError(f"GDEX API {method} {shown}: HTTP {moved.code} without an address; not followed")
        target = urllib.parse.urlsplit(urllib.parse.urljoin(current, moved.location.strip()))
        api = urllib.parse.urlsplit(self.api)
        try:
            same = (target.scheme == api.scheme and target.hostname is not None and api.hostname is not None
                    and target.hostname.lower() == api.hostname.lower() and target.port == api.port
                    and target.username is None and target.password is None)
        except ValueError:          # an invalid port
            same = False
        if not same:
            where = f"{target.scheme}://{target.hostname}" if target.hostname else "an address without a host"
            raise g.NcarRequestError(f"GDEX API {method} {shown}: HTTP {moved.code} to {_redact(where, self.token)}, "
                                     f"another scheme or host than {api.scheme}://{api.hostname}; not followed (the "
                                     f"token is only sent to the API)")
        url = self._with_token(urllib.parse.urlunsplit(target))
        return url, urllib.parse.urlunsplit(target._replace(query="", fragment=""))

    def _call(self, budget, cancel, logger, method, endpoint, payload=None, timeout=None):
        url = f"{self.api}{endpoint}?{urllib.parse.urlencode({'token': self.token})}"
        shown = f"{self.api}{endpoint}"
        data = None if payload is None else json.dumps(payload).encode()
        headers = {} if payload is None else {"Content-Type": "application/json"}

        def call(attempt):
            target, where = url, shown
            for hop in range(SUBSET_MAX_REDIRECTS + 1):
                ticket = _Ticket(f"api {method} {endpoint}", where, SUBSET_API_LIMIT)
                try:
                    body, _, _, _ = _transfer(budget, ticket, target, self.opener, timeout or SUBSET_API_TIMEOUT,
                                              cancel, method=method, data=data, headers=headers, secret=self.token,
                                              redirects=True)
                    break
                except _Redirect as moved:
                    if hop == SUBSET_MAX_REDIRECTS:
                        raise g.NcarRequestError(f"GDEX API {method} {shown}: more than {SUBSET_MAX_REDIRECTS} "
                                                 f"redirects; not followed further") from None
                    target, where = self._redirect_target(moved, target, method, shown)
                    logger.debug(f"GDEX API {method} {shown}: HTTP {moved.code}, followed to {where}")
            try:
                return json.loads(body.decode("utf-8"))
            except ValueError:
                raise g.NcarContentError(f"{shown}: the API did not answer JSON") from None

        # a submission is never repeated by this loop (one attempt): a lost answer is handled by the caller
        attempts = 1 if endpoint == "submit/" else self.attempts
        answer, _ = _retrying(f"GDEX API {method} {shown}", attempts, self._backoff, cancel, logger, call)
        return answer

    def _listing(self, budget, cancel, logger):
        answer = self._call(budget, cancel, logger, "GET", "status/")
        data = answer.get("data")
        if isinstance(data, list):
            return data
        if "no requests" in " ".join(answer.get("error_messages") or []).lower():
            return []
        raise g.NcarRequestError("GDEX API: the list of requests of the account could not be read")

    def _save_state(self, parts_dir, request):
        g._write_json(os.path.join(parts_dir, request.state_file), request.state)

    def _load_state(self, parts_dir, request):
        path = os.path.join(parts_dir, request.state_file)
        if not os.path.isfile(path):
            return None
        with open(path, encoding="utf-8") as handle:
            state = json.load(handle)
        if state.get("control") != request.control:
            raise g.NcarRetrievalError(f"{path} belongs to another request; move {parts_dir} away to start again")
        return state

    def prepare(self, plan, parts_dir, budget, logger, cancel):
        controls = subset_controls(plan)
        self.requests = [_SubsetRequest(i, c, len(controls) == 1) for i, c in enumerate(controls)]
        for request in self.requests:
            request.state = self._load_state(parts_dir, request)
        # only the requests with pieces still to obtain are used; the others are only recorded
        needed = [r for r in self.requests
                  if any(g._stored_piece(p, parts_dir) is None for p in plan.requests if p.segment_index == r.index)]
        for request in self.requests:
            if request not in needed and request.state is not None:
                self._describe_state(request)
        if not needed:
            return
        self._obtain_requests([r for r in needed if r.state is None], parts_dir, budget, logger, cancel)
        for request in needed:
            self._wait_completed(request, budget, logger, cancel, parts_dir)
            self._download(request, parts_dir, budget, logger, cancel)

    def _describe_state(self, request, **extra):
        state = request.state
        request.details.update(request_id=state["request_id"], created_here=state["created_here"],
                               reused_existing_request=state["reused_existing_request"],
                               submit_answer_lost=state.get("submit_answer_lost", False), **extra)

    def describe(self):
        if len(self.requests) == 1:
            self.details.update(self.requests[0].details)
        elif self.requests:
            self.details["requests"] = [dict(segment=r.index, wlon=r.control["wlon"], elon=r.control["elon"],
                                             **r.details) for r in self.requests]
        return super().describe()

    def _obtain_requests(self, pending, parts_dir, budget, logger, cancel):
        if not pending:
            return
        listed = self._listing(budget, cancel, logger)
        to_submit = []
        for request in pending:
            same = [r for r in listed if rinfo_matches(r.get("rinfo"), request.control)]
            if not same:
                to_submit.append(request)
                continue
            found = same[0]
            status = str(found.get("status"))
            if status == "Completed" or not any(w in status.lower() for w in ("purge", "error", "fail", "cancel")):
                logger.info(f"♻️ GDEX request {found.get('request_index')} ({status}) is identical to this "
                            f"one{request.label} and is reused; it was not created by this retrieval and is never "
                            f"purged by it")
                request.state = {"control": request.control, "request_id": str(found.get("request_index")),
                                 "created_here": False, "reused_existing_request": True,
                                 "status_when_found": status, "files": {}}
                self._save_state(parts_dir, request)
                continue
            raise g.NcarRequestError(
                f"an identical GDEX request ({found.get('request_index')}) is listed as {status!r}: GDEX would "
                f"answer a new submission with it (deduplication). Wait until it leaves the list of requests "
                f"(https://gdex.ucar.edu), then run again; nothing was submitted")
        if not to_submit:
            return
        if len(listed) + len(to_submit) > SUBSET_MAX_LISTED:
            needs = (f"; this area crosses Greenwich and needs {len(to_submit)} requests, one per side"
                     if len(to_submit) > 1 else "")
            raise g.NcarRequestError(
                f"{len(listed)} GDEX requests are listed on the account (GDEX accepts 8 open requests){needs}; "
                f"nothing was submitted. Wait for requests to expire or purge your own, then run again")
        known = {str(r.get("request_index")) for r in listed}
        for request in to_submit:
            self._submit(request, known, parts_dir, budget, logger, cancel)
            known.add(request.state["request_id"])

    def _submit(self, request, known, parts_dir, budget, logger, cancel):
        control = request.control
        submitted_utc = datetime.datetime.now(datetime.timezone.utc)
        request_id, lost = None, False
        try:
            answer = self._call(budget, cancel, logger, "POST", "submit/", payload=control,
                                timeout=SUBSET_SUBMIT_TIMEOUT)
            data = answer.get("data") if isinstance(answer, dict) else None
            request_id = str(data["request_id"]) if isinstance(data, dict) and data.get("request_id") else None
            if request_id is None:
                messages = " ".join(answer.get("error_messages") or []) if isinstance(answer, dict) else ""
                raise g.NcarRequestError(f"GDEX refused the request: {_redact(messages or answer, self.token)}")
        except g.NcarTransferError:
            # no answer: the request may exist anyway; adopt a new listed request that matches, never resubmit
            lost = True
            new = [r for r in self._listing(budget, cancel, logger)
                   if str(r.get("request_index")) not in known and rinfo_matches(r.get("rinfo"), control)]
            if len(new) != 1:
                raise g.NcarTransferError(
                    "the GDEX submission got no answer and no new matching request is listed; nothing more was "
                    "submitted. Check https://gdex.ucar.edu before running again") from None
            request_id = str(new[0]["request_index"])
        reused = request_id in known
        # only a request whose number came back in the answer to this submission counts as created here (and
        # can be purged); one adopted after a lost answer is very likely, but not provably, this submission's
        created = not reused and not lost
        if reused:
            logger.info(f"♻️ GDEX answered with the existing request {request_id} (deduplicated)")
        elif lost:
            logger.info(f"📮 GDEX request {request_id} adopted after the submission got no answer; this retrieval "
                        f"never purges it (its origin is not proven)")
        request.state = {"control": control, "request_id": request_id, "created_here": created,
                         "reused_existing_request": reused, "submitted_utc": g._iso(submitted_utc),
                         "submit_answer_lost": lost, "files": {}}
        self._save_state(parts_dir, request)
        logger.info(f"📮 GDEX request {request_id}{request.label} {'reused' if reused else 'submitted'}")

    def _wait_completed(self, request, budget, logger, cancel, parts_dir):
        rid = request.state["request_id"]
        started, last, observations = self.clock(), None, []
        while True:
            answer = self._call(budget, cancel, logger, "GET", f"status/{rid}")
            data = answer.get("data") if isinstance(answer, dict) else None
            status = data.get("status") if isinstance(data, dict) else None
            elapsed = self.clock() - started
            if status != last:
                observations.append({"seconds": round(elapsed, 1), "status": status})
                logger.info(f"   GDEX request {rid}: {status}")
                last = status
            if status == "Completed":
                break
            if status is None or any(w in str(status).lower() for w in ("purge", "error", "fail", "cancel")):
                raise g.NcarRequestError(f"GDEX request {rid} is {status!r}; it cannot be downloaded. Remove "
                                         f"{os.path.join(parts_dir, request.state_file)} to submit a new request")
            if elapsed > SUBSET_REMOTE_DEADLINE:
                raise g.NcarTransferError(f"GDEX request {rid} not completed after {elapsed:.0f} s ({status}); it is "
                                          f"kept: running again resumes with it")
            g._wait(SUBSET_POLL[0] if elapsed < SUBSET_POLL[2] else SUBSET_POLL[1], cancel)
        request.details["remote_status"] = observations

    def _download(self, request, parts_dir, budget, logger, cancel):
        state = request.state
        rid = state["request_id"]
        answer = self._call(budget, cancel, logger, "GET", f"get_req_files/{rid}")
        listed = ((answer.get("data") or {}).get("web_files") or []) if isinstance(answer, dict) else []
        files = [f for f in listed if str(f.get("web_path", "")).endswith(".nc")]
        if not files:
            raise g.NcarContentError(f"GDEX request {rid} lists no NetCDF files")
        listed_bytes = sum(int(f.get("size") or 0) for f in files)
        pending_bytes = sum(int(f.get("size") or 0) for f in files
                            if os.path.basename(f["web_path"]) not in state["files"])
        room = budget.limit - budget.used - budget.reserved
        if pending_bytes > room:
            raise g.NcarTransferError(f"GDEX request {rid} lists {pending_bytes / 1024 ** 3:.2f} GiB still to "
                                      f"download, more than the {room / 1024 ** 3:.2f} GiB left by the byte limit; "
                                      f"nothing was downloaded (raise --max-gib)")
        raw_dir = os.path.join(parts_dir, f"subset_{rid}")
        os.makedirs(raw_dir, exist_ok=True)
        for item in files:
            name = os.path.basename(item["web_path"])
            size = int(item.get("size") or 0)
            path = os.path.join(raw_dir, name)
            kept = state["files"].get(name)
            if kept and os.path.isfile(path) and g._sha256(path) == kept["sha256"]:
                continue
            ticket = _Ticket(f"file {name}", item["web_path"], size + g.HEADER_ALLOWANCE)

            def get(attempt, item=item, path=path, size=size, ticket=ticket):
                with open(f"{path}.part", "wb") as sink:
                    _, got, sha, http = _transfer(budget, ticket, item["web_path"], self.opener, self.timeout,
                                                  cancel, sink=sink, secret=self.token)
                if size and got != size:
                    raise g._Retryable(f"{got} bytes instead of the {size} listed")
                os.replace(f"{path}.part", path)
                return sha, http, got

            try:
                (sha, http, got), used = _retrying(ticket.key, self.attempts, self._backoff, cancel, logger, get)
            finally:
                if os.path.exists(f"{path}.part"):
                    os.remove(f"{path}.part")
            state["files"][name] = {"url": item["web_path"], "bytes": got, "sha256": sha, "attempts": used}
            self._save_state(parts_dir, request)
            logger.info(f"  ⬇️ {name}: {got / 1024 ** 2:.1f} MiB")
        for name in state["files"]:
            match = _FILE_NAME.search(name)
            if match:
                self.files[(request.index, match["code"], match["start"][:8])] = os.path.join(raw_dir, name)
        self._describe_state(request, files=len(files), listed_bytes=listed_bytes)

    def fetch(self, request, parts_dir, budget, logger, cancel):
        started_utc = datetime.datetime.now(datetime.timezone.utc)
        started = time.monotonic()
        subset = self.requests[request.segment_index]
        path = self.files.get((request.segment_index, request.field.code, f"{request.day:%Y%m%d}"))
        if path is None:
            raise g.NcarContentError(f"{request.key}: no {request.field.name} file of {request.day} in GDEX request "
                                     f"{subset.state['request_id']}")
        with g._NETCDF_LOCK, netCDF4.Dataset(path) as source:
            variable = source.variables[request.field.name]
            variable.set_auto_maskandscale(False)
            attrs = {k: variable.getncattr(k) for k in variable.ncattrs()}
            times = np.asarray(source.variables["time"][:], dtype=np.int64)
            lats = np.asarray(source.variables["latitude"][:], dtype=np.float64)
            lons = np.asarray(source.variables["longitude"][:], dtype=np.float64) % 360.0
            levels = np.asarray(source.variables["level"][:], dtype=np.float64)
            want_t = g._hours_since_origin(g._moments(request.day, request.hours))
            index = {}
            for name, have, want in (("time", times, want_t), ("level", levels, np.asarray(g.PRESSURE_LEVELS_HPA)),
                                     ("latitude", lats, np.asarray(request.latitudes)),
                                     ("longitude", lons, np.asarray(request.segment.response) % 360.0)):
                position = {round(float(v), 6): i for i, v in enumerate(have)}
                try:
                    index[name] = [position[round(float(v), 6)] for v in want]
                except KeyError:
                    raise g.NcarContentError(f"{request.key}: the GDEX file {os.path.basename(path)} lacks "
                                             f"{name} values of the request") from None
            values = variable[index["time"]][:, index["level"]][:, :, index["latitude"]][:, :, :, index["longitude"]]
        http = {"method": self.name, "request_id": subset.state["request_id"], "file": os.path.basename(path),
                "file_sha256": subset.state["files"][os.path.basename(path)]["sha256"]}
        return _keep_piece(request, parts_dir, np.asarray(values, dtype=np.float32), attrs, http, 1,
                           time.monotonic() - started, started_utc)

    def finish(self, record, logger):
        if not self.purge:
            return {"purge": {"requested": False}}
        if len(self.requests) == 1:
            return {"purge": self._purge(self.requests[0], logger)}
        return {"purge": {"requested": True,
                          "requests": [dict(self._purge(r, logger), segment=r.index) for r in self.requests]}}

    def _purge(self, request, logger):
        """Purge one request if this retrieval created it; returns what to record."""
        if request.state is None:
            logger.info(f"GDEX request{request.label} unknown (every piece was kept and no request state was "
                        f"found): not purged")
            return {"requested": True, "done": False, "reason": "no request state in the parts directory"}
        rid = request.state["request_id"]
        if not request.state.get("created_here"):
            reason = ("adopted after a lost submission answer (origin not proven)"
                      if request.state.get("submit_answer_lost") and not request.state.get("reused_existing_request")
                      else "not created by this retrieval")
            logger.info(f"GDEX request {rid}: {reason}; not purged")
            return {"requested": True, "request_id": rid, "done": False, "reason": reason}
        budget = g.ByteBudget(16 * SUBSET_API_LIMIT)
        cancel = threading.Event()
        entry = {"requested": True, "request_id": rid}
        try:
            answer = self._call(budget, cancel, logger, "DELETE", f"purge/{rid}")
            entry["answer"] = answer.get("data") if isinstance(answer, dict) else answer
            listed = [r for r in self._listing(budget, cancel, logger) if str(r.get("request_index")) == rid]
            entry["listed_after"] = listed[0].get("status") if listed else None
            entry["removed_from_listing"] = not listed
        except Exception as error:   # the file and its record are already in place: the failure is recorded
            entry["error"] = _redact(f"{type(error).__name__}: {error}", self._token)
        state = entry.get("listed_after")
        logger.info(f"🧹 GDEX request {rid}: purge asked; "
                    + ("outcome unknown (error recorded)" if "error" in entry
                       else "no longer listed" if entry.get("removed_from_listing") else f"still listed as {state!r}"))
        return entry
