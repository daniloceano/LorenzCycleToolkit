"""
Unit tests of the run status record (run_status.json).

They check the contract of the command line: the global state and the exit
status for success, partial output (with and without --allow-partial-plots),
failed scientific computation and interruption; that a figure is declared
generated only when its files were written by the call; and that the status
file is written atomically, so an interrupted or failed write never leaves a
complete-looking status behind.
"""

import json
import os

import pytest

from lorenzcycletoolkit import create_arg_parser
from src.utils import run_status as rs
from src.utils.run_status import RunStatus


def _args(*extra):
    return create_arg_parser().parse_args(["samples/case.nc", *extra])


def _record(tmp_path, *extra):
    return RunStatus(_args(*extra), tmp_path, argv=["lorenzcycletoolkit.py", "samples/case.nc", *extra])


def _complete_science(status, tmp_path):
    outputs = {}
    for name in [n for n, p in status.products.items() if p["required"]]:
        path = tmp_path / f"{name}.csv"
        path.write_text("x\n")
        outputs[name] = [path]
    status.record_scientific_outputs(outputs)


def _writer(path, text="png"):
    def write():
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(text)
    return write


def _read(tmp_path):
    return json.loads((tmp_path / "run_status.json").read_text(encoding="utf-8"))


def test_products_without_plots_are_not_requested(tmp_path):
    status = _record(tmp_path, "-f")
    assert [n for n, p in status.products.items() if p["required"]] == ["results_csv", "vertical_level_csvs"]
    assert {p["status"] for p in status.products.values() if not p["required"]} == {"not_requested"}
    track = _record(tmp_path, "-t")
    assert "trackfile" in track.products and track.products["trackfile"]["required"]
    assert set(rs.MOVING_FIGURES) <= set(track.products)


def test_success_without_plots(tmp_path):
    status = _record(tmp_path, "-f")
    status.write()
    assert _read(tmp_path)["status"] == "running"
    assert _read(tmp_path)["exit_code"] is None
    _complete_science(status, tmp_path)
    status.finish()
    data = _read(tmp_path)
    assert data["schema_version"] == 1
    assert (data["status"], data["exit_code"]) == ("success", 0)
    assert data["scientific"] == {"status": "success", "error": None}
    assert data["plots_requested"] is False
    assert data["summary"]["not_requested"] == len(rs.FIXED_FIGURES)
    assert data["summary"]["success"] == 2
    assert sum(data["summary"].values()) == len(data["products"])


def test_all_requested_figures_generated(tmp_path):
    status = _record(tmp_path, "-f", "-p")
    _complete_science(status, tmp_path)
    for name in rs.FIXED_FIGURES:
        target = tmp_path / "Figures" / f"{name}.png"
        assert status.run_product(name, _writer(target), expected=[f"Figures/{name}.png"])
    status.finish()
    data = _read(tmp_path)
    assert (data["status"], data["exit_code"]) == ("success", 0)
    for product in data["products"]:
        if not product["required"]:
            assert product["status"] == "success"
            assert product["written"] == product["expected"]
    assert status.final_message("Fixed", 1.0).startswith("🎉 Analysis complete!")


def test_failed_figure_gives_partial_and_exit_3(tmp_path):
    status = _record(tmp_path, "-f", "-p")
    _complete_science(status, tmp_path)

    def broken():
        raise RuntimeError("cartography unavailable")

    assert not status.run_product("box_limits_map", broken, expected=["Figures/box_limits.png"])
    for name in rs.FIXED_FIGURES:
        if name != "box_limits_map":
            status.run_product(name, _writer(tmp_path / "Figures" / f"{name}.png"), expected=[f"Figures/{name}.png"])
    status.finish()
    data = _read(tmp_path)
    assert (data["status"], data["exit_code"]) == ("partial", rs.EXIT_PARTIAL)
    assert data["scientific"]["status"] == "success"
    product = next(p for p in data["products"] if p["name"] == "box_limits_map")
    assert product["status"] == "failed"
    assert product["error"] == {"type": "RuntimeError", "message": "cartography unavailable"}
    assert product["missing"] == ["Figures/box_limits.png"]
    assert product["written"] == []
    message = status.final_message("Fixed", 1.0)
    assert "PARTIAL" in message and "box_limits_map" in message and "🎉" not in message


def test_allow_partial_plots_exits_0_but_stays_partial(tmp_path):
    status = _record(tmp_path, "-f", "-p", "--allow-partial-plots")
    _complete_science(status, tmp_path)
    status.run_product("timeseries", lambda: None, expected=["Figures/timeseries/a.png"])
    for name in rs.FIXED_FIGURES[1:]:
        status.run_product(name, _writer(tmp_path / "Figures" / f"{name}.png"), expected=[f"Figures/{name}.png"])
    status.finish()
    data = _read(tmp_path)
    assert (data["status"], data["exit_code"], data["allow_partial_plots"]) == ("partial", 0, True)
    assert "--allow-partial-plots" in status.final_message("Fixed", 1.0)
    assert "PARTIAL" in status.final_message("Fixed", 1.0)


def test_figure_returning_without_writing_is_failed(tmp_path):
    """A plotting function that logs an error and returns must not count as generated."""
    status = _record(tmp_path, "-t", "-p")
    assert not status.run_product("track_map", lambda: None, expected=["Figures/track.png"])
    product = status.products["track_map"]
    assert product["status"] == "failed"
    assert product["error"] is None
    assert "returned without writing Figures/track.png" in product["reason"]


def test_stale_file_from_an_earlier_run_is_not_counted(tmp_path):
    old = tmp_path / "Figures" / "track.png"
    old.parent.mkdir(parents=True)
    old.write_text("from a previous run")
    status = _record(tmp_path, "-t", "-p")
    assert not status.run_product("track_map", lambda: None, expected=["Figures/track.png"])
    assert status.products["track_map"]["missing"] == ["Figures/track.png"]
    assert status.run_product("track_map", _writer(old, "new"), expected=["Figures/track.png"]) is True
    # a later success does not hide the earlier failure of the same product
    assert status.products["track_map"]["status"] == "failed"


def test_partial_files_with_exception_are_listed_and_failed(tmp_path):
    status = _record(tmp_path, "-f", "-p")

    def half():
        _writer(tmp_path / "Figures" / "boxplots" / "a.png")()
        raise ValueError("second figure failed")

    assert not status.run_product("boxplots", half, expected=["Figures/boxplots/a.png", "Figures/boxplots/b.png"])
    product = status.products["boxplots"]
    assert product["written"] == ["Figures/boxplots/a.png"]
    assert product["missing"] == ["Figures/boxplots/b.png"]
    assert product["status"] == "failed"


def test_time_step_products_accumulate(tmp_path):
    status = _record(tmp_path, "-t", "-p")
    status.run_product("domain_maps", _writer(tmp_path / "Figures/boxes/box_1.png"),
                       expected=["Figures/boxes/box_1.png"], step="1")

    def broken():
        raise OSError("no shapefile")

    status.run_product("domain_maps", broken, expected=["Figures/boxes/box_2.png"], step="2")
    status.run_product("domain_maps", broken, expected=["Figures/boxes/box_3.png"], step="3")
    product = status.products["domain_maps"]
    assert product["status"] == "failed"
    assert product["expected"] == [f"Figures/boxes/box_{i}.png" for i in (1, 2, 3)]
    assert product["written"] == ["Figures/boxes/box_1.png"]
    assert product["missing"] == ["Figures/boxes/box_2.png", "Figures/boxes/box_3.png"]
    assert [f["step"] for f in product["failures"]] == ["2", "3"]
    assert product["reason"] == "2 calls failed; first: OSError: no shapefile"


def test_dependency_failure_is_recorded_on_both_products(tmp_path):
    status = _record(tmp_path, "-t", "-p")

    def cyclophaser():
        raise ValueError("polyorder must be less than window_length.")

    status.run_product("periods", cyclophaser, expected=["periods.csv"])
    status.dependency_failed("lps", "periods", expected=["Figures/LPS/LPS_1h.png"])
    assert status.products["lps"]["status"] == "failed"
    assert status.products["lps"]["depends_on"] == ["periods"]
    assert "requires 'periods', which failed" in status.products["lps"]["reason"]
    assert status.products["periods"]["consequences"] == ["lps"]


def test_scientific_error_is_failed_never_partial(tmp_path):
    status = _record(tmp_path, "-f", "-p")
    status.write()
    with pytest.raises(ZeroDivisionError):
        with status.guard():
            raise ZeroDivisionError("division by zero in a term")
    data = _read(tmp_path)
    assert (data["status"], data["exit_code"]) == ("failed", 1)
    assert data["scientific"]["status"] == "failed"
    assert data["scientific"]["error"]["type"] == "ZeroDivisionError"
    assert {p["status"] for p in data["products"]} == {"failed"}
    assert all("run stopped" in p["reason"] for p in data["products"])


def test_error_after_complete_science_is_still_failed(tmp_path):
    status = _record(tmp_path, "-f", "-p")
    _complete_science(status, tmp_path)
    with pytest.raises(KeyError):
        with status.guard():
            raise KeyError("unexpected")
    data = _read(tmp_path)
    assert data["status"] == "failed"
    assert data["scientific"]["status"] == "success"
    assert data["error"]["type"] == "KeyError"


def test_interruption_is_failed_without_exit_code(tmp_path):
    status = _record(tmp_path, "-t")
    with pytest.raises(KeyboardInterrupt):
        with status.guard():
            raise KeyboardInterrupt
    data = _read(tmp_path)
    assert (data["status"], data["interrupted"], data["exit_code"]) == ("failed", True, None)


def test_finish_without_scientific_outputs_is_an_error(tmp_path):
    status = _record(tmp_path, "-f")
    with pytest.raises(RuntimeError, match="required scientific outputs were not recorded"):
        with status.guard():
            status.finish()
    assert _read(tmp_path)["status"] == "failed"


def test_missing_required_output_is_a_scientific_failure(tmp_path):
    status = _record(tmp_path, "-f")
    with pytest.raises(RuntimeError, match="not written"):
        status.record_scientific_outputs({"results_csv": [tmp_path / "absent.csv"]})
    assert status.products["results_csv"]["status"] == "failed"
    assert status.scientific["status"] == "running"


def test_failed_final_write_keeps_the_running_status(tmp_path, monkeypatch):
    status = _record(tmp_path, "-f")
    status.write()
    _complete_science(status, tmp_path)

    def no_replace(src, dst):
        raise OSError("disk full")

    monkeypatch.setattr(rs.os, "replace", no_replace)
    with pytest.raises(OSError, match="disk full"):
        status.finish()
    monkeypatch.undo()
    assert _read(tmp_path)["status"] == "running"
    assert [p.name for p in tmp_path.iterdir() if p.name.endswith(".tmp")] == []


def test_write_is_atomic(tmp_path, monkeypatch):
    """The status file is replaced in one step: readers never see a half-written file."""
    calls = []
    real_replace = os.replace

    def spy(src, dst):
        calls.append((os.path.basename(src), os.path.basename(dst)))
        assert json.loads(open(src, encoding="utf-8").read())["schema_version"] == 1
        real_replace(src, dst)

    monkeypatch.setattr(rs.os, "replace", spy)
    _record(tmp_path, "-f").write()
    assert calls == [(f".run_status.json.{os.getpid()}.tmp", "run_status.json")]
