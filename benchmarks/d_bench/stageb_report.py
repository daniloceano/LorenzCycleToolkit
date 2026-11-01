"""
D-Bench Stage B: public table and figures from the audited campaign records.

Reads the audit written by the campaign's read-only audit (``audit/audit.json``) and the output
NetCDF files (for their dimensions only), and writes, without any network access:

- ``stageB_runs.csv``: one row per planned run (10 cases x 3 methods), no private paths;
- ``stageB_time.png``: time to a validated NetCDF file per case and method;
- ``stageB_transfer.png``: bytes received divided by the float32 content V of the case.

Usage:
    python benchmarks/d_bench/stageb_report.py AUDIT_JSON CAMPAIGN_ROOT OUT_DIR
"""

import csv
import json
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import netCDF4  # noqa: E402

METHODS = ("gas", "kerchunk", "ncss")
LABEL = {"gas": "Get a Subset", "kerchunk": "kerchunk", "ncss": "NCSS"}
COLOR = {"gas": "#2a78d6", "kerchunk": "#eb6834", "ncss": "#1baf7a"}    # categorical slots 1-3
MARKER = {"gas": "o", "kerchunk": "s", "ncss": "D"}
INK, MUTED, GRID = "#0b0b0b", "#52514e", "#e4e3df"
G = 2 ** 30


def content_bytes(root, case):
    """V = 5 fields x levels x lat x lon x times x 4 bytes, from the kerchunk output of the case."""
    with netCDF4.Dataset(Path(root) / f"kerchunk/runs/{case}/era5_{case}_kerchunk.nc") as nc:
        sizes = [len(d) for d in nc.dimensions.values()]
    total = 1
    for s in sizes:
        total *= s
    return 5 * 4 * total


def main(audit_json, root, out):
    out = Path(out)
    out.mkdir(parents=True, exist_ok=True)
    audit = json.loads(Path(audit_json).read_text())
    runs = audit["runs"]
    cases = list(dict.fromkeys(r["case"] for r in runs))
    volume = {c: content_bytes(root, c) for c in cases}
    rows = []
    for r in runs:
        lec = r.get("lec_vs_cds")
        rows.append({
            "case": r["case"], "cds_reference": "yes" if r["has_cds"] else "no", "method": LABEL[r["method"]],
            "position_in_block": r["position"], "class": r["cls"],
            "total_to_validated_s": "" if r.get("total_s") is None else round(r["total_s"]),
            "remote_wait_s": "" if r.get("remote_s") is None else round(r["remote_s"]),
            "transfer_s": next((round(r[k]) for k in ("download_s", "network_s") if r.get(k) is not None), ""),
            "acquisition_s": "" if r.get("acquisition_s") is None else round(r["acquisition_s"]),
            "validation_s": "" if r.get("validation_s") is None else round(r["validation_s"]),
            "received_gib": "" if not r.get("bytes") else round(r["bytes"] / G, 2),
            "received_over_V": "" if not r.get("bytes") else round(r["bytes"] / volume[r["case"]], 2),
            "V_gib": round(volume[r["case"]] / G, 2),
            "http_requests": r.get("http_requests") or "",
            "retries": "" if r.get("retries") is None else r["retries"],
            "peak_rss_gib": "" if not r.get("peak_rss") else round(r["peak_rss"] / G, 2),
            "fields_within_tolerance": "" if r["outcome"] != "success" else ("yes" if r.get("field_valid") else "no"),
            "lec_vs_cds": ("" if r["outcome"] != "success" else "no CDS reference" if not r["has_cds"]
                           else "equivalent" if lec and lec["equivalent"]
                           else f"not equivalent (max {lec['max_rel_tol']:.2f} x tolerance)" if lec else ""),
        })
    with open(out / "stageB_runs.csv", "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)

    order = sorted(cases, key=lambda c: volume[c])            # by size of the case
    plt.rcParams.update({"font.size": 9, "axes.edgecolor": MUTED, "axes.labelcolor": INK,
                         "xtick.color": MUTED, "ytick.color": MUTED})

    def base(ax, ylabel):
        ax.set_xticks(range(len(order)))
        ax.set_xticklabels([f"{c}\n{volume[c] / G:.1f} GiB" + ("" if c != "19800318" else "\n(no CDS)")
                            for c in order], fontsize=7)
        ax.set_ylabel(ylabel)
        ax.grid(axis="y", color=GRID, linewidth=0.8)
        ax.set_axisbelow(True)
        for side in ("top", "right"):
            ax.spines[side].set_visible(False)

    fig, ax = plt.subplots(figsize=(8.5, 3.8), dpi=150)
    offset = {"gas": -0.2, "kerchunk": 0.0, "ncss": 0.2}
    for m in METHODS:
        for i, c in enumerate(order):
            r = next(x for x in runs if x["case"] == c and x["method"] == m)
            x = i + offset[m]
            if r["outcome"] == "success":
                ax.plot(x, r["total_s"] / 60, MARKER[m], color=COLOR[m], markersize=7,
                        markeredgecolor="white", markeredgewidth=1.0)
            else:
                text = {"failed": "×", "not_run": "–"}.get(r["outcome"], "inc")
                ax.annotate(text, (x, 2), ha="center", va="bottom", color=MUTED, fontsize=8)
        ax.plot([], [], MARKER[m], color=COLOR[m], label=LABEL[m], markersize=7)
    base(ax, "minutes to validated NetCDF")
    ax.set_ylim(0, None)
    ax.legend(frameon=False, ncol=3, loc="upper left")
    ax.set_title("D-Bench Stage B, 2026-10-31: one cold run per case and method (× failed, – not run, "
                 "inc = not submitted)", fontsize=8.5, color=INK, loc="left")
    fig.tight_layout()
    fig.savefig(out / "stageB_time.png")
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(8.5, 3.4), dpi=150)
    for m in METHODS:
        xs, ys = [], []
        for i, c in enumerate(order):
            r = next(x for x in runs if x["case"] == c and x["method"] == m)
            if r["outcome"] == "success" and r.get("bytes"):
                xs.append(i + offset[m])
                ys.append(r["bytes"] / volume[c])
        ax.plot(xs, ys, MARKER[m], color=COLOR[m], markersize=7, markeredgecolor="white", markeredgewidth=1.0,
                label=LABEL[m])
    ax.axhline(1.0, color=MUTED, linewidth=0.8, linestyle=(0, (3, 3)))
    base(ax, "bytes received / V")
    ax.set_ylim(0, ax.get_ylim()[1] * 1.15)
    ax.legend(frameon=False, ncol=3, loc="upper right")
    ax.set_title("Data received relative to the float32 content V of each case (dashed line: V)",
                 fontsize=8.5, color=INK, loc="left")
    fig.tight_layout()
    fig.savefig(out / "stageB_transfer.png")
    plt.close(fig)


if __name__ == "__main__":
    main(*sys.argv[1:4])
