#!/usr/bin/env python3
"""
06_plot_resources.py
====================
Generate publication-quality time and disk-space comparison figures for
the three genome-similarity pipelines evaluated in this study:

    KMC (exact)         — k-mer counting + exact pairwise comparison
    FracMinHash         — sourmash manysketch + sourmash pairwise
    AlphaMaxGeomHash    — 03_alphamaxgeom_sketch.sh + 04_alphamaxgeom_pairwise.py

Two figures are produced:

  Figure 1 — Computation time (stacked bars)
      "Indexing/sketching" time + "Pairwise similarity" time for each method.
      The pairwise bar is coloured differently from the indexing bar.
      All values are reported in minutes and labelled explicitly.

  Figure 2 — Disk space (stacked bars)
      "Index/sketch storage" + "Pairwise output" size for each method.

Both figures include a secondary y-axis (right side) showing the ratio
relative to KMC pairwise computation time / KMC pairwise output size,
so the reader can immediately gauge the resource reduction.

Timing sources
--------------
  AMG sketch   : data/GTDB/alphamaxgeom_sketches/sketch_run_stats.json
  AMG pairwise : data/GTDB/alphamaxgeom_pairwise/pairwise_run_stats.json
  FMH sketch   : scripts/make_sourmash_sketches.log   (GNU time block)
  FMH pairwise : scripts/compute_sourmash_pairwise.log (GNU time block)
  KMC pairwise : scripts/02_kmc_pairwise.log           (Run Summary block)
  KMC counting : scripts/01_kmc_count.log              (GNU time block, if present)

Disk-space sources
------------------
  All sizes are measured live from the filesystem at runtime so they are
  always current.

Extensibility
-------------
To add a new method (e.g. MinHash), add one entry to METHOD_SPECS below.
Each entry declares where to find timing information and which filesystem
paths to measure for disk usage.  No other code needs to change.
"""

import argparse
import json
import logging
import os
import re
import sys
from pathlib import Path

import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%H:%M:%S",
)
log = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Matplotlib publication style
# ---------------------------------------------------------------------------
plt.rcParams.update({
    "font.family":    "DejaVu Sans",
    "font.size":      9,
    "axes.titlesize": 10,
    "axes.labelsize": 9,
    "xtick.labelsize": 8,
    "ytick.labelsize": 8,
    "savefig.dpi":    300,
    "savefig.bbox":   "tight",
    "pdf.fonttype":   42,
    "ps.fonttype":    42,
})

# ---------------------------------------------------------------------------
# Method colour palette  (extend this list to add more methods)
# ---------------------------------------------------------------------------
# Chosen to be distinguishable for deuteranopia / protanopia colour blindness.
METHOD_COLORS = {
    "KMC (exact)":     {"index": "#2166AC", "pairwise": "#92C5DE"},  # blues
    "FracMinHash":     {"index": "#D6604D", "pairwise": "#F4A582"},  # reds
    "AlphaMaxGeomHash":{"index": "#4DAC26", "pairwise": "#B8E186"},  # greens
    # Future MinHash:
    # "MinHash":       {"index": "#762A83", "pairwise": "#C2A5CF"},  # purples
}

# ---------------------------------------------------------------------------
# Parsers
# ---------------------------------------------------------------------------

def _parse_gnu_time_elapsed(log_path: Path) -> float | None:
    """
    Extract wall-clock seconds from a GNU `time` block in a log file.
    Looks for a line like:
        Elapsed (wall clock) time (h:mm:ss or m:ss): 7:18.19
    or:
        Elapsed (wall clock) time (h:mm:ss or m:ss): 2:05:10
    Returns seconds as float, or None if not found.
    """
    if not log_path.exists():
        return None
    # The GNU time format is either m:ss.fraction or h:mm:ss[.fraction].
    # Group 1 = hours-or-minutes, group 2 = minutes-or-seconds (may have
    # decimal), group 3 = seconds (present only in h:mm:ss format).
    pattern = re.compile(
        r"Elapsed \(wall clock\) time.*?:\s*(\d+):(\d+(?:\.\d+)?)(?::(\d+(?:\.\d+)?))?$"
    )
    with open(log_path) as f:
        for line in f:
            m = pattern.search(line.strip())
            if m:
                a, b, c = m.group(1), m.group(2), m.group(3)
                if c is None:
                    # m:ss format
                    return int(a) * 60 + float(b)
                else:
                    # h:mm:ss format
                    return int(a) * 3600 + int(b) * 60 + float(c)
    return None


def _parse_run_summary_hms(log_path: Path, label: str = "Wall-clock time") -> float | None:
    """
    Extract wall-clock seconds from a "Run Summary" block, e.g.:
        Wall-clock time     : 02:05:10 (hh:mm:ss)
    Returns seconds as float, or None if not found.
    """
    if not log_path.exists():
        return None
    pattern = re.compile(
        rf"{re.escape(label)}\s*:\s*(\d+):(\d+):(\d+)"
    )
    with open(log_path) as f:
        for line in f:
            m = pattern.search(line)
            if m:
                h, mn, s = int(m.group(1)), int(m.group(2)), int(m.group(3))
                return h * 3600 + mn * 60 + s
    return None


def _parse_gnu_time_ram_kb(log_path: Path) -> float | None:
    """Extract peak RAM in kB from GNU time block (Maximum resident set size)."""
    if not log_path.exists():
        return None
    pattern = re.compile(r"Maximum resident set size \(kbytes\):\s*(\d+)")
    with open(log_path) as f:
        for line in f:
            m = pattern.search(line)
            if m:
                return float(m.group(1))
    return None


def _json_stat(json_path: Path, key: str) -> float | None:
    """Read a single numeric value from a JSON stats file."""
    if not json_path.exists():
        return None
    with open(json_path) as f:
        data = json.load(f)
    return data.get(key)


def disk_bytes(*paths: Path) -> float:
    """
    Return the total size in bytes of one or more files or directories.
    Uses os.walk for directories (faster than subprocess du on large trees
    when the directory list is already cached by the OS).
    Returns 0.0 for paths that do not exist.
    """
    total = 0
    for p in paths:
        p = Path(p)
        if not p.exists():
            continue
        if p.is_file():
            total += p.stat().st_size
        elif p.is_dir():
            for dirpath, _, filenames in os.walk(p):
                for fn in filenames:
                    try:
                        total += os.path.getsize(os.path.join(dirpath, fn))
                    except OSError:
                        pass
    return float(total)


# ---------------------------------------------------------------------------
# Data collection
# ---------------------------------------------------------------------------

def collect_stats(base: Path) -> dict:
    """
    Gather timing and disk-space statistics for all three methods from the
    known file locations under *base* (the project root directory).

    Returns a nested dict:
        {
          "MethodName": {
            "index_seconds":    float or None,
            "pairwise_seconds": float or None,
            "index_bytes":      float,
            "pairwise_bytes":   float,
            "index_ram_kb":     float or None,   # peak RAM
            "pairwise_ram_kb":  float or None,
          },
          ...
        }
    """
    scripts = base / "scripts"
    data    = base / "data" / "GTDB"

    stats = {}

    # ------------------------------------------------------------------
    # KMC (exact)
    # ------------------------------------------------------------------
    kmc_index_s    = _parse_gnu_time_elapsed(scripts / "01_kmc_count.log")
    kmc_pairwise_s = _parse_run_summary_hms(scripts / "02_kmc_pairwise.log",
                                             "Wall-clock time")
    if kmc_pairwise_s is None:
        # Fallback: parse the GNU time block if the Run Summary block is absent
        kmc_pairwise_s = _parse_gnu_time_elapsed(scripts / "02_kmc_pairwise.log")

    kmc_index_ram_kb    = _parse_gnu_time_ram_kb(scripts / "01_kmc_count.log")
    kmc_pairwise_ram_kb = _parse_gnu_time_ram_kb(scripts / "02_kmc_pairwise.log")

    stats["KMC (exact)"] = {
        "index_seconds":    kmc_index_s,
        "pairwise_seconds": kmc_pairwise_s,
        "index_bytes":      disk_bytes(data / "kmc_dbs"),
        "pairwise_bytes":   disk_bytes(data / "kmc_pairwise"),
        "index_ram_kb":     kmc_index_ram_kb,
        "pairwise_ram_kb":  kmc_pairwise_ram_kb,
    }

    # ------------------------------------------------------------------
    # FracMinHash (sourmash)
    # ------------------------------------------------------------------
    fmh_index_s    = _parse_gnu_time_elapsed(scripts / "make_sourmash_sketches.log")
    fmh_pairwise_s = _parse_gnu_time_elapsed(scripts / "compute_sourmash_pairwise.log")

    fmh_index_ram_kb    = _parse_gnu_time_ram_kb(scripts / "make_sourmash_sketches.log")
    fmh_pairwise_ram_kb = _parse_gnu_time_ram_kb(scripts / "compute_sourmash_pairwise.log")

    stats["FracMinHash"] = {
        "index_seconds":    fmh_index_s,
        "pairwise_seconds": fmh_pairwise_s,
        "index_bytes":      disk_bytes(data / "gtdb_genomes_reps_r226_sketches.siz.zip"),
        "pairwise_bytes":   disk_bytes(data / "gtdb_pairwise_containment.csv"),
        "index_ram_kb":     fmh_index_ram_kb,
        "pairwise_ram_kb":  fmh_pairwise_ram_kb,
    }

    # ------------------------------------------------------------------
    # AlphaMaxGeomHash
    # ------------------------------------------------------------------
    amg_sketch_json   = data / "alphamaxgeom_sketches"  / "sketch_run_stats.json"
    amg_pairwise_json = data / "alphamaxgeom_pairwise"  / "pairwise_run_stats.json"

    amg_index_s    = _json_stat(amg_sketch_json,   "wall_clock_seconds")
    amg_pairwise_s = _json_stat(amg_pairwise_json, "wall_clock_seconds")

    stats["AlphaMaxGeomHash"] = {
        "index_seconds":    amg_index_s,
        "pairwise_seconds": amg_pairwise_s,
        "index_bytes":      disk_bytes(data / "alphamaxgeom_sketches"),
        "pairwise_bytes":   disk_bytes(data / "alphamaxgeom_pairwise"),
        "index_ram_kb":     None,   # tracked by GNU time in sketch_run_stats
        "pairwise_ram_kb":  None,
    }
    # AMG sketch RAM is stored in JSON as a human string; parse it separately
    if amg_sketch_json.exists():
        with open(amg_sketch_json) as f:
            jd = json.load(f)
        ram_str = jd.get("peak_ram", "")   # e.g. "33.27 MB"
        m = re.match(r"([\d.]+)\s*(\w+)", ram_str)
        if m:
            val, unit = float(m.group(1)), m.group(2).upper()
            multiplier = {"KB": 1, "MB": 1024, "GB": 1024**2}.get(unit, 1)
            stats["AlphaMaxGeomHash"]["index_ram_kb"] = val * multiplier

    return stats


# ---------------------------------------------------------------------------
# Plotting helpers
# ---------------------------------------------------------------------------

def _fmt_time(s: float | None) -> str:
    if s is None:
        return "N/A"
    h, rem = divmod(int(s), 3600)
    m, sc  = divmod(rem, 60)
    if h:
        return f"{h}h {m:02d}m"
    if m:
        return f"{m}m {sc:02d}s"
    return f"{sc}s"


def _fmt_bytes(b: float) -> str:
    for unit, thresh in [("TB", 1e12), ("GB", 1e9), ("MB", 1e6), ("KB", 1e3)]:
        if b >= thresh:
            return f"{b/thresh:.1f} {unit}"
    return f"{b:.0f} B"


def _bar_label(ax, bar, text: str, fontsize: int = 7):
    """Place a text label just above a bar."""
    ax.text(
        bar.get_x() + bar.get_width() / 2,
        bar.get_height() + bar.get_y() + bar.get_height() * 0.03,
        text,
        ha="center", va="bottom", fontsize=fontsize, clip_on=True,
    )


# ---------------------------------------------------------------------------
# Figure 1 — Computation time  (log y-scale)
# ---------------------------------------------------------------------------

# Small positive floor used as the visual baseline for log-scale bars so
# that bars never extend to log(0) = −∞.  0.1 min ≈ 6 seconds.
_TIME_FLOOR_MIN  = 0.1
_DISK_FLOOR_GB   = 0.001   # 1 MB

def _log_center(lo: float, hi: float) -> float:
    """Geometric mean: the visually centred position on a log axis."""
    return np.sqrt(max(lo, 1e-12) * max(hi, 1e-12))


def plot_time(stats: dict, out_dir: Path,
              method_order: list[str] | None = None):
    """
    Stacked bar chart (log y-scale) of wall-clock time:
        darker segment = indexing / k-mer counting / sketching
        lighter segment = pairwise similarity computation

    Log scale is essential: KMC (~125 min) is ~7× slower than FMH (~21 min)
    and AMG (~18 min), which are nearly indistinguishable on a linear axis.
    The ×N ratio label above each bar gives the speed-up relative to KMC
    pairwise computation time.
    """
    if method_order is None:
        method_order = list(stats.keys())

    methods = [m for m in method_order if m in stats]
    n       = len(methods)
    x       = np.arange(n)
    width   = 0.55
    floor   = _TIME_FLOOR_MIN

    fig, ax = plt.subplots(figsize=(6.5, 4.5))

    kmc_pairwise_s = stats.get("KMC (exact)", {}).get("pairwise_seconds") or 1.0

    for i, method in enumerate(methods):
        s      = stats[method]
        col    = METHOD_COLORS.get(method, {"index": "#888", "pairwise": "#BBB"})
        idx_s  = s["index_seconds"]    or 0.0
        pw_s   = s["pairwise_seconds"] or 0.0
        idx_m  = idx_s / 60
        pw_m   = pw_s  / 60

        # Stacked bars starting from the log floor.
        # For KMC (idx=0): the index segment has zero height and is invisible.
        ax.bar(i, idx_m,       width, bottom=floor,
               color=col["index"],    edgecolor="white", linewidth=0.5)
        ax.bar(i, pw_m,        width, bottom=floor + idx_m,
               color=col["pairwise"], edgecolor="white", linewidth=0.5)

        bar_top = floor + idx_m + pw_m

        # Label each segment at its geometric (log) centre if there is
        # enough vertical room (segment spans > factor 2 on log axis).
        if idx_m > 0 and (floor + idx_m) / floor > 2.0:
            cy = _log_center(floor, floor + idx_m)
            ax.text(i, cy, _fmt_time(idx_s),
                    ha="center", va="center", fontsize=7,
                    color="white", fontweight="bold")
        if pw_m > 0 and bar_top / (floor + idx_m) > 2.0:
            cy = _log_center(floor + idx_m, bar_top)
            ax.text(i, cy, _fmt_time(pw_s),
                    ha="center", va="center", fontsize=7,
                    color="white", fontweight="bold")

        # For thin segments (< factor 2), annotate outside instead.
        if idx_m > 0 and (floor + idx_m) / floor <= 2.0:
            ax.text(i, floor + idx_m * 1.05,
                    f"sketch: {_fmt_time(idx_s)}",
                    ha="center", va="bottom", fontsize=6.5, color="#333")
        if pw_m > 0 and bar_top / (floor + idx_m) <= 2.0:
            ax.text(i, bar_top * 1.05,
                    f"pairwise: {_fmt_time(pw_s)}",
                    ha="center", va="bottom", fontsize=6.5, color="#333")

        # N/A note for missing indexing time
        if s["index_seconds"] is None:
            ax.text(i, floor * 1.3, "indexing\nN/A",
                    ha="center", va="bottom", fontsize=6, color="#666",
                    style="italic")

        # ×N ratio above bar (relative to KMC pairwise)
        total_s = (s["index_seconds"] or 0) + (s["pairwise_seconds"] or 0)
        if total_s > 0:
            ratio = total_s / kmc_pairwise_s
            ax.text(i, bar_top * 1.6,
                    f"×{ratio:.2f}",
                    ha="center", va="bottom", fontsize=9, fontweight="bold")

    ax.set_yscale("log")
    ax.set_xticks(x)
    ax.set_xticklabels(methods, fontsize=9)
    ax.set_ylabel("Wall-clock time (minutes, log scale)", fontsize=9)
    # Y limits: floor/2 at the bottom, top bar × enough headroom for ×N labels
    all_tops = [
        _TIME_FLOOR_MIN + (s.get("index_seconds") or 0) / 60
                        + (s.get("pairwise_seconds") or 0) / 60
        for s in stats.values()
    ]
    ax.set_ylim(_TIME_FLOOR_MIN / 2, max(all_tops) * 5)
    ax.yaxis.grid(True, which="both", linestyle="--", linewidth=0.4, alpha=0.6)
    ax.set_axisbelow(True)

    legend_handles = [
        mpatches.Patch(color="#555555", label="Indexing / k-mer counting / sketching"),
        mpatches.Patch(color="#AAAAAA", label="Pairwise similarity computation"),
    ]
    ax.legend(handles=legend_handles, fontsize=8, loc="upper left")
    ax.annotate(
        f"×N = total time relative to KMC pairwise ({_fmt_time(kmc_pairwise_s)})",
        xy=(0.01, 0.01), xycoords="axes fraction",
        fontsize=7, va="bottom", color="#444",
    )

    for suffix in (".pdf", ".png"):
        fig.savefig(out_dir / f"resources_time{suffix}")
        log.info("Saved: %s", out_dir / f"resources_time{suffix}")
    plt.close(fig)


# ---------------------------------------------------------------------------
# Figure 2 — Disk space  (log y-scale)
# ---------------------------------------------------------------------------

def plot_disk(stats: dict, out_dir: Path,
              method_order: list[str] | None = None):
    """
    Stacked bar chart (log y-scale) of disk space:
        darker segment = index / sketch storage
        lighter segment = pairwise results storage

    Log scale is essential: KMC k-mer databases (3.9 TB) are ~880× larger
    than FMH sketches (4 GB) and ~260× larger than AMG sketches (15 GB).
    On a linear axis FMH and AMG are invisible slivers.
    """
    if method_order is None:
        method_order = list(stats.keys())

    methods = [m for m in method_order if m in stats]
    n       = len(methods)
    x       = np.arange(n)
    width   = 0.55
    floor   = _DISK_FLOOR_GB

    fig, ax = plt.subplots(figsize=(6.5, 4.5))

    kmc_pw_bytes = stats.get("KMC (exact)", {}).get("pairwise_bytes") or 1.0

    for i, method in enumerate(methods):
        s      = stats[method]
        col    = METHOD_COLORS.get(method, {"index": "#888", "pairwise": "#BBB"})
        idx_b  = s["index_bytes"]
        pw_b   = s["pairwise_bytes"]
        idx_gb = idx_b / 1e9
        pw_gb  = pw_b  / 1e9

        ax.bar(i, idx_gb, width, bottom=floor,
               color=col["index"],    edgecolor="white", linewidth=0.5)
        ax.bar(i, pw_gb,  width, bottom=floor + idx_gb,
               color=col["pairwise"], edgecolor="white", linewidth=0.5)

        bar_top = floor + idx_gb + pw_gb

        # Segment labels at log centres (only if segment spans > factor 2)
        if idx_gb > 0 and (floor + idx_gb) / floor > 2.0:
            cy = _log_center(floor, floor + idx_gb)
            ax.text(i, cy, _fmt_bytes(idx_b),
                    ha="center", va="center", fontsize=7,
                    color="white", fontweight="bold")
        if pw_gb > 0 and bar_top / (floor + idx_gb) > 2.0:
            cy = _log_center(floor + idx_gb, bar_top)
            ax.text(i, cy, _fmt_bytes(pw_b),
                    ha="center", va="center", fontsize=7,
                    color="white", fontweight="bold")

        # Thin segments: annotate outside
        if idx_gb > 0 and (floor + idx_gb) / floor <= 2.0:
            ax.text(i, floor + idx_gb * 1.05,
                    f"idx: {_fmt_bytes(idx_b)}",
                    ha="center", va="bottom", fontsize=6.5, color="#333")
        if pw_gb > 0 and bar_top / (floor + idx_gb) <= 2.0:
            ax.text(i, bar_top * 1.05,
                    f"pw: {_fmt_bytes(pw_b)}",
                    ha="center", va="bottom", fontsize=6.5, color="#333")

        # ×N ratio relative to KMC pairwise output
        total_b = idx_b + pw_b
        ratio   = total_b / kmc_pw_bytes
        ax.text(i, bar_top * 1.6,
                f"×{ratio:.0f}",
                ha="center", va="bottom", fontsize=9, fontweight="bold")

    ax.set_yscale("log")
    ax.set_xticks(x)
    ax.set_xticklabels(methods, fontsize=9)
    ax.set_ylabel("Disk space (GB, log scale)", fontsize=9)
    all_tops = [
        _DISK_FLOOR_GB + s["index_bytes"] / 1e9 + s["pairwise_bytes"] / 1e9
        for s in stats.values()
    ]
    ax.set_ylim(_DISK_FLOOR_GB / 2, max(all_tops) * 5)
    ax.yaxis.grid(True, which="both", linestyle="--", linewidth=0.4, alpha=0.6)
    ax.set_axisbelow(True)

    legend_handles = [
        mpatches.Patch(color="#555555", label="Index / sketch storage"),
        mpatches.Patch(color="#AAAAAA", label="Pairwise results storage"),
    ]
    ax.legend(handles=legend_handles, fontsize=8, loc="upper left")
    ax.annotate(
        f"×N = total storage relative to KMC pairwise output ({_fmt_bytes(kmc_pw_bytes)})",
        xy=(0.01, 0.01), xycoords="axes fraction",
        fontsize=7, va="bottom", color="#444",
    )

    for suffix in (".pdf", ".png"):
        fig.savefig(out_dir / f"resources_disk{suffix}")
        log.info("Saved: %s", out_dir / f"resources_disk{suffix}")
    plt.close(fig)


# ---------------------------------------------------------------------------
# Figure 3 — Accuracy vs. Resource trade-off scatter
# ---------------------------------------------------------------------------

def plot_tradeoff(stats: dict, sanity_json: Path, out_dir: Path,
                  method_order: list[str] | None = None):
    """
    Scatter plot of accuracy (Pearson r for Jaccard vs. KMC) against total
    computation time and total disk space — one point per method.
    This gives a quick visual summary of the accuracy/resource trade-off.

    Accuracy values are read from the sanity-check summary JSON written by
    05_sanity_check.py.
    """
    if not sanity_json.exists():
        log.warning("Sanity-check JSON not found (%s) — skipping trade-off plot.", sanity_json)
        return

    with open(sanity_json) as f:
        summary = json.load(f)

    # Build metric → pearson_r lookup
    r_by_label = {}
    for s in summary.get("error_statistics", []):
        r_by_label[s["metric"]] = s.get("pearson_r")

    amg_r = r_by_label.get("AMG vs KMC — Jaccard")
    fmh_r = r_by_label.get("FMH vs KMC — Jaccard")

    accuracy = {
        "KMC (exact)":      1.0,          # ground truth, r = 1 by definition
        "FracMinHash":      fmh_r,
        "AlphaMaxGeomHash": amg_r,
    }

    if method_order is None:
        method_order = list(stats.keys())
    methods = [m for m in method_order if m in stats]

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(9.0, 3.8))

    kmc_pw_s = stats.get("KMC (exact)", {}).get("pairwise_seconds") or 1.0
    kmc_pw_b = stats.get("KMC (exact)", {}).get("pairwise_bytes")   or 1.0

    for method in methods:
        s     = stats[method]
        col   = METHOD_COLORS.get(method, {"index": "#888"})["index"]
        r_val = accuracy.get(method)
        if r_val is None:
            continue

        total_s = ((s["index_seconds"] or 0) + (s["pairwise_seconds"] or 0)) / 60
        total_b = (s["index_bytes"] + s["pairwise_bytes"]) / 1e9

        for ax, xval, xlabel in [
            (ax1, total_s, "Total wall-clock time (minutes, log scale)"),
            (ax2, total_b, "Total disk space (GB, log scale)"),
        ]:
            ax.scatter(xval, r_val, s=120, color=col, zorder=5,
                       edgecolors="black", linewidths=0.6)
            ax.annotate(
                method, (xval, r_val),
                textcoords="offset points", xytext=(6, 4),
                fontsize=8, color=col,
            )
            ax.set_xscale("log")
            ax.set_xlabel(xlabel, fontsize=9)
            ax.set_ylabel("Pearson r  (Jaccard vs. KMC exact)", fontsize=9)
            ax.set_ylim(0.992, 1.002)
            ax.yaxis.grid(True, linestyle="--", linewidth=0.4, alpha=0.6)
            ax.xaxis.grid(True, which="both", linestyle="--", linewidth=0.4, alpha=0.4)
            ax.set_axisbelow(True)

    ax1.set_title("Accuracy vs. Computation Time", fontsize=10)
    ax2.set_title("Accuracy vs. Disk Space",        fontsize=10)
    fig.tight_layout()

    for suffix in (".pdf", ".png"):
        fig.savefig(out_dir / f"resources_tradeoff{suffix}")
        log.info("Saved: %s", out_dir / f"resources_tradeoff{suffix}")
    plt.close(fig)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def parse_args():
    p = argparse.ArgumentParser(
        description="Plot computation time and disk-space comparison figures.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("--base-dir", required=True,
                   help="Project root directory (parent of scripts/ and data/)")
    p.add_argument("--sanity-json", default="",
                   help="Path to sanity_check_summary.json from 05_sanity_check.py "
                        "(used for the accuracy vs. resource trade-off scatter plot)")
    p.add_argument("--output", required=True,
                   help="Output directory for figures")
    return p.parse_args()


def main():
    args    = parse_args()
    base    = Path(args.base_dir)
    out_dir = Path(args.output)
    out_dir.mkdir(parents=True, exist_ok=True)

    log.info("Collecting resource statistics ...")
    stats = collect_stats(base)

    # Log what was found
    for method, s in stats.items():
        idx_t = _fmt_time(s["index_seconds"])
        pw_t  = _fmt_time(s["pairwise_seconds"])
        log.info(
            "  %-22s  index=%s  pairwise=%s  "
            "index_disk=%s  pairwise_disk=%s",
            method, idx_t, pw_t,
            _fmt_bytes(s["index_bytes"]), _fmt_bytes(s["pairwise_bytes"]),
        )

    # The desired display order (most to least resource-intensive)
    method_order = ["KMC (exact)", "FracMinHash", "AlphaMaxGeomHash"]

    log.info("Generating time figure ...")
    plot_time(stats, out_dir, method_order)

    log.info("Generating disk-space figure ...")
    plot_disk(stats, out_dir, method_order)

    sanity_json = Path(args.sanity_json) if args.sanity_json else (
        base / "data" / "GTDB" / "alphamaxgeom_pairwise" / "sanity_check" /
        "sanity_check_summary.json"
    )
    log.info("Generating accuracy vs. resource trade-off figure ...")
    plot_tradeoff(stats, sanity_json, out_dir, method_order)

    # Write a machine-readable summary of all gathered numbers
    summary_path = out_dir / "resource_summary.json"
    with open(summary_path, "w") as f:
        json.dump(stats, f, indent=2)
    log.info("Resource summary written to %s", summary_path)
    log.info("Done.")


if __name__ == "__main__":
    main()
