#!/usr/bin/env python3
"""
10_plot_resources.py
====================
Generate publication-quality time and disk-space comparison figures for the
genome-similarity pipelines evaluated in this study:

    KMC (exact)               -- k-mer counting + exact pairwise comparison
    BottomK                   -- 03_bottomk_sketch.sh + 06_bottomk_pairwise.py
    FracMinHash (kmer-sketch) -- 04_fracminhash_sketch.sh + 07_fracminhash_pairwise.py
    AlphaMaxGeomHash          -- 05_alphamaxgeom_sketch.sh + 08_alphamaxgeom_pairwise.py
    Sourmash FracMinHash      -- make_sourmash_sketches.sh + compute_sourmash_pairwise.sh
                                 (included only if --include-sourmash is set)

Two figures are produced:

  Figure 1 -- Computation time (stacked bars, log scale)
      "Indexing/sketching" time + "Pairwise similarity" time for each method.

  Figure 2 -- Disk space (stacked bars, log scale)
      "Index/sketch storage" + "Pairwise output" size for each method.

Figure 3 -- Accuracy vs. resource trade-off scatter plot.

Log scale is essential: KMC (~2 h, ~3.9 TB) is orders of magnitude larger
than the sketching methods; a linear axis makes them invisible slivers.

Timing sources
--------------
  AMG sketch    : data/GTDB/alphamaxgeom_sketches/sketch_run_stats.json
  AMG pairwise  : data/GTDB/alphamaxgeom_pairwise_thr0001/pairwise_run_stats.json
  BK sketch     : data/GTDB/bottomk_sketches/sketch_run_stats.json
  BK pairwise   : data/GTDB/bottomk_pairwise_thr0001/pairwise_run_stats.json
  FMH_ks sketch : data/GTDB/fracminhash_sketches/sketch_run_stats.json
  FMH_ks pairwise: data/GTDB/fracminhash_pairwise_thr0001/pairwise_run_stats.json
  Sourmash sketch  : scripts/make_sourmash_sketches.log   (GNU time block)
  Sourmash pairwise: scripts/compute_sourmash_pairwise.log (GNU time block)
  KMC pairwise  : scripts/02_kmc_pairwise.log             (Run Summary block)
  KMC counting  : scripts/01_kmc_count.log                (GNU time block)

Disk-space sources
------------------
  All sizes are measured live from the filesystem at runtime.

Extensibility
-------------
To add a new method, add one entry to METHOD_COLORS and one block in
collect_stats().  No other code needs to change.
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
# Method colour palette  (extend this dict to add more methods)
# Chosen to be distinguishable for deuteranopia / protanopia colour blindness.
# ---------------------------------------------------------------------------
METHOD_COLORS = {
    "KMC (exact)":              {"index": "#2166AC", "pairwise": "#92C5DE"},  # blues
    "MinHash":                  {"index": "#762A83", "pairwise": "#C2A5CF"},  # purples
    "FracMinHash":              {"index": "#D6604D", "pairwise": "#F4A582"},  # reds
    "AlphaMaxGeomHash":         {"index": "#4DAC26", "pairwise": "#B8E186"},  # greens
    # Sourmash FracMinHash (off by default; add via --include-sourmash):
    "Sourmash FracMinHash":     {"index": "#8C510A", "pairwise": "#DFC27D"},  # browns
}

# ---------------------------------------------------------------------------
# Parsers
# ---------------------------------------------------------------------------

def _parse_gnu_time_elapsed(log_path):
    """
    Extract wall-clock seconds from a GNU `time` block in a log file.
    Looks for a line like:
        Elapsed (wall clock) time (h:mm:ss or m:ss): 7:18.19
    or:
        Elapsed (wall clock) time (h:mm:ss or m:ss): 2:05:10
    Returns seconds as float, or None if not found.
    """
    if not Path(log_path).exists():
        return None
    pattern = re.compile(
        r"Elapsed \(wall clock\) time.*?:\s*(\d+):(\d+(?:\.\d+)?)(?::(\d+(?:\.\d+)?))?$"
    )
    with open(log_path) as f:
        for line in f:
            m = pattern.search(line.strip())
            if m:
                a, b, c = m.group(1), m.group(2), m.group(3)
                if c is None:
                    return int(a) * 60 + float(b)
                else:
                    return int(a) * 3600 + int(b) * 60 + float(c)
    return None


def _parse_pipeline_step_elapsed(log_path, step_label):
    """
    Extract elapsed wall-clock seconds for a named step from a master pipeline
    log that uses the convention:
        === START <step_label> <date> ===
        ...
        === END   <step_label> <date> ===

    Parses the date strings with dateutil and returns the difference in seconds,
    or None if the markers are not found.
    """
    from dateutil import parser as _duparser
    if not Path(log_path).exists():
        return None
    start_pat = re.compile(rf"=== START\s+{re.escape(step_label)}\s+(.*?)\s*===")
    end_pat   = re.compile(rf"=== END\s+{re.escape(step_label)}\s+(.*?)\s*===")
    t_start = t_end = None
    with open(log_path) as f:
        for line in f:
            if t_start is None:
                m = start_pat.search(line)
                if m:
                    try:
                        t_start = _duparser.parse(m.group(1))
                    except Exception:
                        pass
            else:
                m = end_pat.search(line)
                if m:
                    try:
                        t_end = _duparser.parse(m.group(1))
                    except Exception:
                        pass
                    break
    if t_start is not None and t_end is not None:
        return (t_end - t_start).total_seconds()
    return None


def _parse_run_summary_hms(log_path, label="Wall-clock time"):
    """
    Extract wall-clock seconds from a "Run Summary" block, e.g.:
        Wall-clock time     : 02:05:10 (hh:mm:ss)
    Returns seconds as float, or None if not found.
    """
    if not Path(log_path).exists():
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


def _parse_gnu_time_ram_kb(log_path):
    """Extract peak RAM in kB from GNU time block (Maximum resident set size)."""
    if not Path(log_path).exists():
        return None
    pattern = re.compile(r"Maximum resident set size \(kbytes\):\s*(\d+)")
    with open(log_path) as f:
        for line in f:
            m = pattern.search(line)
            if m:
                return float(m.group(1))
    return None


def _json_stat(json_path, key):
    """Read a single numeric value from a JSON stats file."""
    json_path = Path(json_path)
    if not json_path.exists():
        return None
    with open(json_path) as f:
        data = json.load(f)
    return data.get(key)


def _json_cpu_seconds(json_path, time_key="wall_clock_seconds", cores_key="cores"):
    """
    Compute CPU seconds from a run_stats JSON file.
    Returns wall_clock_seconds * cores, or None if either field is missing.
    """
    json_path = Path(json_path)
    if not json_path.exists():
        return None
    with open(json_path) as f:
        data = json.load(f)
    wall_s = data.get(time_key)
    cores  = data.get(cores_key)
    if wall_s is None or cores is None:
        return None
    return float(wall_s) * float(cores)


def disk_bytes(*paths):
    """
    Return the total size in bytes of one or more files or directories.
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

def collect_stats(base, include_sourmash=False):
    """
    Gather timing and disk-space statistics for all methods.

    Returns a nested dict:
        {
          "MethodName": {
            "index_seconds":    float or None,
            "pairwise_seconds": float or None,
            "index_bytes":      float,
            "pairwise_bytes":   float,
            "index_ram_kb":     float or None,
            "pairwise_ram_kb":  float or None,
          },
          ...
        }

    Methods whose data files are absent are still included with None timings
    and 0 disk usage so that bar-chart code doesn't need to special-case them.
    Only KMC and at least one sketch method are truly expected to be present.
    """
    base    = Path(base)
    scripts = base / "scripts"
    data    = base / "data" / "GTDB"

    stats = {}

    # ------------------------------------------------------------------
    # KMC (exact)
    # No run_stats JSON exists for KMC; timing is read from log files and
    # converted to CPU seconds using the known core counts:
    #   counting  : 96 parallel jobs × 4 threads/job = 384 cores
    #   pairwise  : 192 worker processes (--cores 192 in 02_kmc_pairwise.sh)
    # The KMC counting log crashed before completion, so index time is None.
    #
    # Pairwise timing: prefer the master pipeline log (full_pipeline.log)
    # because 02_kmc_pairwise.log may contain a stale run from a previous
    # threshold setting.  Fall back to 02_kmc_pairwise.log if the pipeline
    # log does not have a matching START/END block.
    # ------------------------------------------------------------------
    _KMC_COUNT_CORES    = 384   # 96 parallel_jobs × 4 threads each
    _KMC_PAIRWISE_CORES = 192   # --cores 192 in 02_kmc_pairwise.sh

    kmc_index_wall_s    = _parse_gnu_time_elapsed(scripts / "01_kmc_count.log")
    kmc_pairwise_wall_s = _parse_pipeline_step_elapsed(
        scripts / "full_pipeline.log", "02_kmc_pairwise"
    )
    if kmc_pairwise_wall_s is None:
        kmc_pairwise_wall_s = _parse_run_summary_hms(scripts / "02_kmc_pairwise.log",
                                                      "Wall-clock time")
    if kmc_pairwise_wall_s is None:
        kmc_pairwise_wall_s = _parse_gnu_time_elapsed(scripts / "02_kmc_pairwise.log")

    kmc_index_s    = (kmc_index_wall_s    * _KMC_COUNT_CORES
                      if kmc_index_wall_s    is not None else None)
    kmc_pairwise_s = (kmc_pairwise_wall_s * _KMC_PAIRWISE_CORES
                      if kmc_pairwise_wall_s is not None else None)

    stats["KMC (exact)"] = {
        "index_seconds":    kmc_index_s,
        "pairwise_seconds": kmc_pairwise_s,
        "index_bytes":      disk_bytes(data / "kmc_dbs"),
        "pairwise_bytes":   disk_bytes(data / "kmc_pairwise_thr0001"),
        "index_ram_kb":     None,
        "pairwise_ram_kb":  None,
    }

    # ------------------------------------------------------------------
    # BottomK
    # CPU seconds = wall_clock_seconds × parallel_jobs (sketch)
    #             = wall_clock_seconds × cores          (pairwise)
    # ------------------------------------------------------------------
    bk_sketch_json   = data / "bottomk_sketches"       / "sketch_run_stats.json"
    bk_pairwise_json = data / "bottomk_pairwise_thr0001" / "pairwise_run_stats.json"

    stats["MinHash"] = {
        "index_seconds":    _json_cpu_seconds(bk_sketch_json,
                                              "wall_clock_seconds", "parallel_jobs"),
        "pairwise_seconds": _json_cpu_seconds(bk_pairwise_json,
                                              "wall_clock_seconds", "cores"),
        "index_bytes":      disk_bytes(data / "bottomk_sketches"),
        "pairwise_bytes":   disk_bytes(data / "bottomk_pairwise_thr0001"),
        "index_ram_kb":     None,
        "pairwise_ram_kb":  None,
    }

    # ------------------------------------------------------------------
    # AlphaMaxGeomHash
    # ------------------------------------------------------------------
    amg_sketch_json   = data / "alphamaxgeom_sketches"          / "sketch_run_stats.json"
    amg_pairwise_json = data / "alphamaxgeom_pairwise_thr0001"  / "pairwise_run_stats.json"

    stats["AlphaMaxGeomHash"] = {
        "index_seconds":    _json_cpu_seconds(amg_sketch_json,
                                              "wall_clock_seconds", "parallel_jobs"),
        "pairwise_seconds": _json_cpu_seconds(amg_pairwise_json,
                                              "wall_clock_seconds", "cores"),
        "index_bytes":      disk_bytes(data / "alphamaxgeom_sketches"),
        "pairwise_bytes":   disk_bytes(data / "alphamaxgeom_pairwise_thr0001"),
        "index_ram_kb":     None,
        "pairwise_ram_kb":  None,
    }

    # ------------------------------------------------------------------
    # FracMinHash (kmer-sketch binary)
    # ------------------------------------------------------------------
    fmh_ks_sketch_json   = data / "fracminhash_sketches"        / "sketch_run_stats.json"
    fmh_ks_pairwise_json = data / "fracminhash_pairwise_thr0001" / "pairwise_run_stats.json"

    stats["FracMinHash"] = {
        "index_seconds":    _json_cpu_seconds(fmh_ks_sketch_json,
                                              "wall_clock_seconds", "parallel_jobs"),
        "pairwise_seconds": _json_cpu_seconds(fmh_ks_pairwise_json,
                                              "wall_clock_seconds", "cores"),
        "index_bytes":      disk_bytes(data / "fracminhash_sketches"),
        "pairwise_bytes":   disk_bytes(data / "fracminhash_pairwise_thr0001"),
        "index_ram_kb":     None,
        "pairwise_ram_kb":  None,
    }

    # ------------------------------------------------------------------
    # Sourmash FracMinHash (optional; excluded by default)
    # Sourmash has no run_stats JSON; fall back to log-file wall-clock time.
    # The core count for sourmash is not recorded, so CPU time is unavailable
    # and we store None.
    # ------------------------------------------------------------------
    if include_sourmash:
        stats["Sourmash FracMinHash"] = {
            "index_seconds":    None,
            "pairwise_seconds": None,
            "index_bytes":  disk_bytes(data / "gtdb_genomes_reps_r226_sketches.siz.zip"),
            "pairwise_bytes": disk_bytes(data / "gtdb_pairwise_containment.csv"),
            "index_ram_kb":   None,
            "pairwise_ram_kb":  None,
        }

    return stats


# ---------------------------------------------------------------------------
# Plotting helpers
# ---------------------------------------------------------------------------

def _fmt_time(s):
    if s is None:
        return "N/A"
    h, rem = divmod(int(s), 3600)
    m, sc  = divmod(rem, 60)
    if h:
        return f"{h}h {m:02d}m"
    if m:
        return f"{m}m {sc:02d}s"
    return f"{sc}s"


def _fmt_bytes(b):
    for unit, thresh in [("TB", 1e12), ("GB", 1e9), ("MB", 1e6), ("KB", 1e3)]:
        if b >= thresh:
            return f"{b/thresh:.1f} {unit}"
    return f"{b:.0f} B"


def _bar_label(ax, bar, text, fontsize=7):
    ax.text(
        bar.get_x() + bar.get_width() / 2,
        bar.get_height() + bar.get_y() + bar.get_height() * 0.03,
        text,
        ha="center", va="bottom", fontsize=fontsize, clip_on=True,
    )


# ---------------------------------------------------------------------------
# Figure 1 -- Computation time  (log y-scale)
# ---------------------------------------------------------------------------

_TIME_FLOOR_H    = 0.001   # ~3.6 seconds — keeps log axis stable near zero
_DISK_FLOOR_GB   = 0.001   # 1 MB

def _log_center(lo, hi):
    """Geometric mean: the visually centred position on a log axis."""
    return np.sqrt(max(lo, 1e-12) * max(hi, 1e-12))


def plot_time(stats, out_dir, method_order=None):
    """
    Bar chart (log y-scale) of pairwise-similarity CPU time (hours) only.
    CPU seconds = wall_clock_seconds × number of cores used.
    Indexing / sketching time is excluded.
    """
    if method_order is None:
        method_order = list(stats.keys())

    methods = [m for m in method_order if m in stats]
    n       = len(methods)
    x       = np.arange(n)
    width   = 0.55
    floor   = _TIME_FLOOR_H

    fig, ax = plt.subplots(figsize=(max(6.5, n * 1.5), 4.5))

    kmc_pairwise_s = stats.get("KMC (exact)", {}).get("pairwise_seconds") or 1.0

    for i, method in enumerate(methods):
        s     = stats[method]
        col   = METHOD_COLORS.get(method, {"index": "#888", "pairwise": "#BBB"})
        pw_s  = s["pairwise_seconds"] or 0.0
        pw_h  = pw_s / 3600

        ax.bar(i, pw_h, width, bottom=floor,
               color=col["index"], edgecolor="white", linewidth=0.5)

        bar_top = floor + pw_h

        # Label inside bar if it's tall enough, otherwise above
        if pw_h > 0 and bar_top / floor > 2.0:
            cy = _log_center(floor, bar_top)
            ax.text(i, cy, _fmt_time(pw_s),
                    ha="center", va="center", fontsize=7,
                    color="white", fontweight="bold")
        elif pw_h > 0:
            ax.text(i, bar_top * 1.05, _fmt_time(pw_s),
                    ha="center", va="bottom", fontsize=6.5, color="#333")

        if pw_s > 0:
            ratio = pw_s / kmc_pairwise_s
            ax.text(i, bar_top * 1.6,
                    f"x{ratio:.3f}",
                    ha="center", va="bottom", fontsize=9, fontweight="bold")

    ax.set_yscale("log")
    ax.set_xticks(x)
    ax.set_xticklabels(methods, fontsize=8, rotation=15, ha="right")
    ax.set_ylabel("Pairwise similarity CPU time (hours, log scale)", fontsize=9)
    all_tops = [
        _TIME_FLOOR_H + (s.get("pairwise_seconds") or 0) / 3600
        for s in stats.values()
    ]
    ax.set_ylim(_TIME_FLOOR_H / 2, max(all_tops) * 5)
    ax.yaxis.grid(True, which="both", linestyle="--", linewidth=0.4, alpha=0.6)
    ax.set_axisbelow(True)

    ax.annotate(
        f"xN = pairwise CPU time relative to KMC pairwise ({_fmt_time(kmc_pairwise_s)})",
        xy=(0.01, 0.01), xycoords="axes fraction",
        fontsize=7, va="bottom", color="#444",
    )

    for suffix in (".pdf", ".png"):
        fig.savefig(out_dir / f"resources_time{suffix}")
        log.info("Saved: %s", out_dir / f"resources_time{suffix}")
    plt.close(fig)


# ---------------------------------------------------------------------------
# Figure 2 -- Disk space  (log y-scale)
# ---------------------------------------------------------------------------

def plot_disk(stats, out_dir, method_order=None):
    """Stacked bar chart (log y-scale) of disk space."""
    if method_order is None:
        method_order = list(stats.keys())

    methods = [m for m in method_order if m in stats]
    n       = len(methods)
    x       = np.arange(n)
    width   = 0.55
    floor   = _DISK_FLOOR_GB

    fig, ax = plt.subplots(figsize=(max(6.5, n * 1.5), 4.5))

    kmc_idx_bytes = stats.get("KMC (exact)", {}).get("index_bytes") or 1.0

    for i, method in enumerate(methods):
        s      = stats[method]
        col    = METHOD_COLORS.get(method, {"index": "#888", "pairwise": "#BBB"})
        idx_b  = s["index_bytes"]
        idx_gb = idx_b / 1e9

        ax.bar(i, idx_gb, width, bottom=floor,
               color=col["index"], edgecolor="white", linewidth=0.5)

        bar_top = floor + idx_gb

        if idx_gb > 0 and bar_top / floor > 2.0:
            cy = _log_center(floor, bar_top)
            ax.text(i, cy, _fmt_bytes(idx_b),
                    ha="center", va="center", fontsize=7,
                    color="white", fontweight="bold")
        else:
            ax.text(i, bar_top * 1.05, _fmt_bytes(idx_b),
                    ha="center", va="bottom", fontsize=6.5, color="#333")

        if idx_b > 0:
            ratio = idx_b / kmc_idx_bytes
            ax.text(i, bar_top * 1.6, f"x{ratio:.3f}",
                    ha="center", va="bottom", fontsize=9, fontweight="bold")

    ax.set_yscale("log")
    ax.set_xticks(x)
    ax.set_xticklabels(methods, fontsize=8, rotation=15, ha="right")
    ax.set_ylabel("Sketch / index size on disk (GB, log scale)", fontsize=9)
    all_tops = [
        _DISK_FLOOR_GB + s["index_bytes"] / 1e9
        for s in stats.values()
    ]
    ax.set_ylim(_DISK_FLOOR_GB / 2, max(all_tops) * 5)
    ax.yaxis.grid(True, which="both", linestyle="--", linewidth=0.4, alpha=0.6)
    ax.set_axisbelow(True)

    ax.annotate(
        f"xN = sketch size relative to KMC index ({_fmt_bytes(kmc_idx_bytes)})",
        xy=(0.01, 0.01), xycoords="axes fraction",
        fontsize=7, va="bottom", color="#444",
    )

    for suffix in (".pdf", ".png"):
        fig.savefig(out_dir / f"resources_disk{suffix}")
        log.info("Saved: %s", out_dir / f"resources_disk{suffix}")
    plt.close(fig)


# ---------------------------------------------------------------------------
# Figure 3 -- Accuracy vs. Resource trade-off scatter
# ---------------------------------------------------------------------------

def plot_tradeoff(stats, sanity_json, out_dir, method_order=None):
    """
    Scatter plot of accuracy (Pearson r for Jaccard vs. KMC) against total
    computation time and total disk space -- one point per method.

    Accuracy values are read from the sanity-check summary JSON written by
    09_sanity_check.py.  The JSON key format is "<label> vs KMC -- Jaccard".
    """
    sanity_json = Path(sanity_json)
    if not sanity_json.exists():
        log.warning("Sanity-check JSON not found (%s) -- skipping trade-off plot.",
                    sanity_json)
        return

    with open(sanity_json) as f:
        summary = json.load(f)

    # Build metric --> pearson_r lookup
    r_by_label = {}
    for s in summary.get("error_statistics", []):
        r_by_label[s["metric"]] = s.get("pearson_r")

    # Map display method name --> accuracy
    accuracy = {"KMC (exact)": 1.0}
    # 09_sanity_check.py writes labels like "AlphaMaxGeomHash vs KMC -- Jaccard"
    label_map = {
        "AlphaMaxGeomHash":          "AlphaMaxGeomHash vs KMC -- Jaccard",
        "MinHash":                   "BottomK vs KMC -- Jaccard",
        "FracMinHash":               "FracMinHash (kmer-sketch) vs KMC -- Jaccard",
        "Sourmash FracMinHash":      "Sourmash FMH vs KMC -- Jaccard",
    }
    for display_name, stat_key in label_map.items():
        if stat_key in r_by_label:
            accuracy[display_name] = r_by_label[stat_key]

    if method_order is None:
        method_order = list(stats.keys())
    methods = [m for m in method_order if m in stats]

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(9.0, 3.8))

    for method in methods:
        s     = stats[method]
        col   = METHOD_COLORS.get(method, {"index": "#888"})["index"]
        r_val = accuracy.get(method)
        if r_val is None:
            continue

        total_s = ((s["index_seconds"] or 0) + (s["pairwise_seconds"] or 0)) / 3600
        total_b = (s["index_bytes"] + s["pairwise_bytes"]) / 1e9

        for ax, xval, xlabel in [
            (ax1, total_s, "Total CPU time (hours, log scale)"),
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
            ax.set_ylim(0.990, 1.002)
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
# Figure 4 -- Pairwise accuracy  (L1 error bar chart, one per metric)
# ---------------------------------------------------------------------------

def plot_accuracy(l1_json, out_dir, metric, method_order=None):
    """
    Bar chart of mean |error| (mean of |estimate − KMC exact| over all pairs)
    for each sketching method.

    Values are read from l1_errors.json written by 10_plot_heatmaps.py.
    KMC (exact) is the reference and is excluded from this chart.

    Parameters
    ----------
    l1_json     : path to l1_errors.json (written by 10_plot_heatmaps.py)
    out_dir     : output directory
    metric      : "jaccard" or "max_containment"
    method_order: display order (KMC will be skipped automatically)
    """
    l1_json = Path(l1_json)
    if not l1_json.exists():
        log.warning("L1 errors file not found (%s) — skipping accuracy bar chart "
                    "for %s.  Run 10_plot_heatmaps.py first.", l1_json, metric)
        return

    with open(l1_json) as f:
        all_l1 = json.load(f)

    l1_for_metric = all_l1.get(metric)
    if not l1_for_metric:
        log.warning("No L1 data for metric '%s' in %s — skipping.", metric, l1_json)
        return

    # Build ordered list, skipping KMC and any method not in the L1 data
    sketch_methods = [m for m in (method_order or list(l1_for_metric.keys()))
                      if m != "KMC (exact)" and m in l1_for_metric]
    if not sketch_methods:
        log.warning("No sketch-method L1 data found for metric '%s' — skipping.", metric)
        return

    n     = len(sketch_methods)
    x     = np.arange(n)
    width = 0.55

    fig, ax = plt.subplots(figsize=(max(5.0, n * 1.5), 4.5))

    for i, method in enumerate(sketch_methods):
        col  = METHOD_COLORS.get(method, {"index": "#888"})["index"]
        l1   = l1_for_metric[method]
        bar  = ax.bar(i, l1, width, color=col, edgecolor="white", linewidth=0.5)
        ax.text(i, l1 * 1.02, f"{l1:.6f}",
                ha="center", va="bottom", fontsize=8, color="#222")

    metric_label = {
        "jaccard":         "Jaccard",
        "max_containment": "Max containment",
    }.get(metric, metric)

    ax.set_xticks(x)
    ax.set_xticklabels(sketch_methods, fontsize=8, rotation=15, ha="right")
    ax.set_ylabel(f"Mean |error|  (|estimate − KMC exact|,  {metric_label})", fontsize=9)
    ax.set_ylim(0, max(l1_for_metric[m] for m in sketch_methods) * 1.2)
    ax.yaxis.grid(True, linestyle="--", linewidth=0.4, alpha=0.6)
    ax.set_axisbelow(True)

    stem = f"accuracy_l1_{metric}"
    for suffix in (".pdf", ".png"):
        fpath = out_dir / (stem + suffix)
        fig.savefig(fpath)
        log.info("Saved: %s", fpath)
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
                   help="Path to sanity_check_summary.json from 09_sanity_check.py "
                        "(used for the accuracy vs. resource trade-off scatter plot). "
                        "If omitted, defaults to the AMG sanity_check directory.")
    p.add_argument("--include-sourmash", action="store_true",
                   help="Include Sourmash FracMinHash in the figures. Off by default.")
    p.add_argument("--output", required=True,
                   help="Output directory for figures")
    return p.parse_args()


def main():
    args    = parse_args()
    base    = Path(args.base_dir)
    out_dir = Path(args.output)
    out_dir.mkdir(parents=True, exist_ok=True)

    log.info("Collecting resource statistics ...")
    stats = collect_stats(base, include_sourmash=args.include_sourmash)

    for method, s in stats.items():
        idx_t = _fmt_time(s["index_seconds"])
        pw_t  = _fmt_time(s["pairwise_seconds"])
        log.info(
            "  %-28s  index=%s  pairwise=%s  "
            "index_disk=%s  pairwise_disk=%s",
            method, idx_t, pw_t,
            _fmt_bytes(s["index_bytes"]), _fmt_bytes(s["pairwise_bytes"]),
        )

    # Display order: KMC first (exact baseline), then sketch methods
    # left-to-right by increasing resource usage (theory prediction):
    #   MinHash < AlphaMaxGeomHash < FracMinHash
    method_order = [
        "KMC (exact)",
        "MinHash",
        "AlphaMaxGeomHash",
        "FracMinHash",
    ]
    if args.include_sourmash:
        method_order.append("Sourmash FracMinHash")

    # Exclude methods with zero data (not yet run) from the figures
    method_order = [
        m for m in method_order
        if m in stats and (
            (stats[m]["index_seconds"] or 0) > 0 or
            (stats[m]["pairwise_seconds"] or 0) > 0 or
            stats[m]["index_bytes"] > 0 or
            stats[m]["pairwise_bytes"] > 0
        )
    ]

    log.info("Generating time figure ...")
    plot_time(stats, out_dir, method_order)

    log.info("Generating disk-space figure ...")
    plot_disk(stats, out_dir, method_order)

    # Sanity JSON: default to the AMG sanity_check dir (backward-compat)
    sanity_json = Path(args.sanity_json) if args.sanity_json else (
        base / "data" / "GTDB" / "alphamaxgeom_pairwise" / "sanity_check" /
        "sanity_check_summary.json"
    )
    # Also look in a combined sanity_check dir if the AMG-specific one is absent
    if not sanity_json.exists():
        alt = base / "data" / "GTDB" / "sanity_check" / "sanity_check_summary.json"
        if alt.exists():
            sanity_json = alt

    log.info("Generating accuracy vs. resource trade-off figure ...")
    plot_tradeoff(stats, sanity_json, out_dir, method_order)

    # L1 accuracy bar charts (one per metric).
    # l1_errors.json is written by 10_plot_heatmaps.py into the same output dir.
    l1_json = out_dir / "l1_errors.json"
    sketch_method_order = [m for m in method_order if m != "KMC (exact)"]
    for metric in ("jaccard", "max_containment"):
        log.info("Generating L1 accuracy bar chart (%s) ...", metric)
        plot_accuracy(l1_json, out_dir, metric, sketch_method_order)

    summary_path = out_dir / "resource_summary.json"
    with open(summary_path, "w") as f:
        json.dump(stats, f, indent=2)
    log.info("Resource summary written to %s", summary_path)
    log.info("Done.")


if __name__ == "__main__":
    main()
