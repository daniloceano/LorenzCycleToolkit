"""
Logging isolation between toolkit runs in the same Python process.

initialize_logging attached a new file handler and a new console handler to
the "lorenzcycletoolkit" logger on every call without removing the previous
ones, so a second run in the same process also wrote its messages to the
log file of the first run and repeated every console line.
"""

import argparse
import logging

from lorenzcycletoolkit import initialize_logging


def _args(infile):
    return argparse.Namespace(verbosity=False, infile=infile)


def _flush(logger):
    for handler in logger.handlers:
        handler.flush()


def test_second_run_does_not_write_to_first_log(tmp_path):
    first, second = tmp_path / "first", tmp_path / "second"
    first.mkdir()
    second.mkdir()

    logger = initialize_logging(str(first), _args("samples/case_first.nc"))
    logger.info("message from the first run")
    logger = initialize_logging(str(second), _args("samples/case_second.nc"))
    logger.info("message from the second run")
    _flush(logger)

    first_log = (first / "log.case_first").read_text()
    second_log = (second / "log.case_second").read_text()
    assert "message from the second run" not in first_log
    assert "message from the first run" not in second_log
    assert first_log.count("message from the first run") == 1
    assert second_log.count("message from the second run") == 1


def test_handlers_do_not_accumulate(tmp_path):
    for name in ("a", "b", "c"):
        (tmp_path / name).mkdir()
        logger = initialize_logging(str(tmp_path / name), _args(f"samples/{name}.nc"))
    file_handlers = [h for h in logger.handlers if isinstance(h, logging.FileHandler)]
    console_handlers = [h for h in logger.handlers if type(h) is logging.StreamHandler]
    assert len(file_handlers) == 1
    assert file_handlers[0].baseFilename == str(tmp_path / "c" / "log.c")
    assert len(console_handlers) == 1


def test_console_messages_are_not_repeated(tmp_path, capsys):
    for name in ("a", "b"):
        (tmp_path / name).mkdir()
        logger = initialize_logging(str(tmp_path / name), _args(f"samples/{name}.nc"))
    capsys.readouterr()
    logger.info("single console message")
    _flush(logger)
    assert capsys.readouterr().err.count("single console message") == 1


def test_previous_file_handler_is_closed(tmp_path):
    (tmp_path / "a").mkdir()
    (tmp_path / "b").mkdir()
    logger = initialize_logging(str(tmp_path / "a"), _args("samples/a.nc"))
    previous = [h for h in logger.handlers if isinstance(h, logging.FileHandler)][0]
    initialize_logging(str(tmp_path / "b"), _args("samples/b.nc"))
    assert previous not in logger.handlers
    assert previous.stream is None


# The next two tests run in this order (pytest keeps the definition order
# within a module). The first configures the toolkit logger and returns
# without closing it, as an integration test does; the second must start
# with no handler attached, which holds only if the handlers of the previous
# test are closed after it (tests/conftest.py).
_CONFIGURED_LOG = {}


def test_handlers_closed_after_test_part1_configure(tmp_path, request):
    logger = initialize_logging(str(tmp_path), _args("samples/part1.nc"))
    logger.info("message of part 1")
    _flush(logger)
    log = tmp_path / "log.part1"
    _CONFIGURED_LOG["path"] = log
    # Also hand the log to the end-of-session check (tests/conftest.py),
    # which fails the run if it grows after this test.
    if not hasattr(request.config, "_lectk_log_sizes"):
        request.config._lectk_log_sizes = {}
    request.config._lectk_log_sizes[log] = log.stat().st_size


def test_handlers_closed_after_test_part2_no_inherited_handler():
    logger = logging.getLogger("lorenzcycletoolkit")
    assert logger.handlers == [], f"handlers inherited from an earlier test: {logger.handlers}"
    logger.warning("message emitted by a later test")
    if "path" in _CONFIGURED_LOG:
        assert "message emitted by a later test" not in _CONFIGURED_LOG["path"].read_text()
