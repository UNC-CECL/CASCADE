"""
Does the northern buffer's shape matter? The first 10 buffer cells follow the real coast north of GIS 90.

    python scripts/hatteras_ms/experiments/HAT_buffer_orientation.py build
    python scripts/hatteras_ms/experiments/HAT_buffer_orientation.py run

build writes, for 1996 and 2010, the padded offset with GIS 91-100 taken from the same
dune line as GIS 1-90 (same per-domain mean, same baseline) and the Hermite closure
wrapping the last 5 right-buffer cells and the 15 left ones. GIS 1-90 are unchanged.
run drives the unchanged runner with that file swapped in, in-process: zeroBE only, one
window at a time, so the buffer is the only thing that differs from the matrix zeroBE runs.
Figures: HAT_candidate_windows_plot.py --buffer. Details: the study's README.md.

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-10-04
"""
from __future__ import annotations

import os
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd

_HERE = Path(__file__).resolve()
PROJECT_ROOT = next(_p for _p in _HERE.parents if (_p / "pyproject.toml").exists())
sys.path.insert(0, str(PROJECT_ROOT / "scripts"))
sys.path.insert(0, str(PROJECT_ROOT))

# --- CONFIG ------------------------------------------------------------------
HINDCAST = PROJECT_ROOT / "scripts" / "hatteras_ms" / "HAT_hindcast.py"
TAG = "island-offset/2026-10-04-northern-buffer-from-real-coast"
EXP_DIR = PROJECT_ROOT / "output" / "raw_runs" / "experiments" / TAG
INPUTS = EXP_DIR / "inputs"
RUN_DIR_TAG = f"{TAG}/runs"
STARTS = (1996, 2010)
N_MEASURED = 10              # right-buffer cells that follow the real coast (GIS 91..100)
BUFFERS = 15
SCENARIO = "full_management"
# -----------------------------------------------------------------------------

from site_layer.hatteras_site_config import _island_offset_file  # noqa: E402
from site_layer.hat_topo_version import BRIE_ROOT, RAW_OFFSET_DIR, dune_line_for_year  # noqa: E402


# Per-domain mean dune distance, as step2_raw_offsets_to_island_offset.calculate_relative_offset takes it
def domain_means(raw, grids):
    out = {}
    for g in grids:
        sub = raw[raw.domain_id == g]
        if sub.empty:
            continue
        first = sub.sort_values("LineID").groupby("LineID").ORIG_LEN.first()
        out[g] = float(first.mean())
    return pd.Series(out)


def offset_path(start):
    return INPUTS / f"Island_Dune_Offsets_{start}_PADDED_120_buffer{N_MEASURED}.csv"


# The base file, and the same with GIS 91-100 measured and the ring re-closed after them
def build_one(start):
    base_file = BRIE_ROOT / Path(_island_offset_file(start)).relative_to("2-brie-offset")
    base_df = pd.read_csv(base_file)              # one column, headed by the year
    base = base_df.iloc[:, 0].to_numpy(float)
    vint = dune_line_for_year(start)
    raw = pd.read_csv(base_file.parent / f"{vint}_duneline_offset_raw.csv")
    raw_ext = pd.read_csv(RAW_OFFSET_DIR / "ext" / f"{vint}_duneline_offset_raw_ext.csv")
    m_base = domain_means(raw, range(1, 91))
    m_ext = domain_means(raw_ext, range(91, 91 + N_MEASURED))
    baseline = m_base.min()
    real = (m_base - baseline).to_numpy()
    ext = (m_ext - baseline).to_numpy()
    if not np.allclose(real, base[BUFFERS:BUFFERS + 90], atol=1e-6):
        raise SystemExit(f"{start}: recomputed GIS 1-90 do not match {base_file.name}")
    measured = np.r_[real, ext]
    # Hermite from GIS 100 round the wrap to GIS 1: 5 right cells, then 15 left
    n_right = BUFFERS - N_MEASURED
    n = n_right + BUFFERS + 1
    p0, p1 = measured[-1], measured[0]
    m0, m1 = (measured[-1] - measured[-2]) * n, (measured[1] - measured[0]) * n
    t = np.arange(1, n) / n
    path = ((2 * t**3 - 3 * t**2 + 1) * p0 + (t**3 - 2 * t**2 + t) * m0
            + (-2 * t**3 + 3 * t**2) * p1 + (t**3 - t**2) * m1)
    padded = np.r_[path[n_right:], measured, path[:n_right]]
    assert len(padded) == len(base) == 120
    INPUTS.mkdir(parents=True, exist_ok=True)
    pd.DataFrame({base_df.columns[0]: padded}).to_csv(offset_path(start), index=False)
    return base, padded


def cmd_build():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from site_layer.hat_figure_style import C_1984, C_1997, INK, apply_style, record_caption, save
    apply_style()
    fig, axes = plt.subplots(2, 1, figsize=(7.2, 5.0), sharex=True, constrained_layout=True)
    for ax, start in zip(axes, STARTS):
        base, new = build_one(start)
        x = np.arange(120) - BUFFERS + 1          # GIS numbering, buffers outside 1-90
        ax.axvspan(0.5, 90.5, color="0.94", lw=0)
        ax.axvspan(90.5, 90.5 + N_MEASURED, color="#f6e3b8", lw=0)
        ax.plot(x, base, color=C_1984, lw=1.2, label="Current: Hermite closure from GIS 90")
        ax.plot(x, new, color=C_1997, lw=1.2, ls=(0, (4, 2)),
                label=f"Test: GIS 91–{90 + N_MEASURED} from the dune line, closure after")
        ax.set_ylabel(f"{start} offset (m)")
        ax.grid(axis="y")
        print(f"{start}: max change in GIS 1-90 {np.abs(new - base)[BUFFERS:BUFFERS + 90].max():.2e} m; "
              f"right buffer {base[-BUFFERS:].round(0)} -> {new[-BUFFERS:].round(0)}")
    axes[0].legend(loc="upper left", frameon=False, fontsize=8)
    axes[-1].set_xlabel("Domain (GIS numbering; buffers below 1 and above 90)")
    out = save(fig, INPUTS / "buffer_orientation_offsets", close=True)[0]
    record_caption(out, (
        "Island offset (metres, dune line distance relative to its minimum over GIS 1–90) for the two windows. "
        f"Red is the model input now: GIS 1–90 from the dune line and 15 buffer cells a side closed by a cubic "
        f"Hermite. Blue dashed is the buffer test: GIS 91–{90 + N_MEASURED} (shaded) from the same dune line and "
        "method, the closure then wrapping the remaining 5 right and 15 left cells. Grey is the modelled reach, "
        "identical in both."))
    print(out)


# Subprocess entry: point the period at the test offset, then run the unchanged runner
def _launch(start):
    import runpy
    from site_layer import hatteras_site_config as sc
    sc.HATTERAS_PERIODS[start]["island_offset_file"] = str(offset_path(start))
    # Metres mode keeps only GIS 1-90 and re-closes the ring itself, which would throw the test buffer
    # away (the first attempt did exactly that); for this file, hand Cascade the padded values as written
    import numpy as np_
    import cascade_pipeline.hindcast as hc
    original = hc.build_island_offset

    def as_written(path, geometry, mode="metres"):
        if Path(path).resolve() == offset_path(start).resolve():
            vals = np_.loadtxt(path, skiprows=1, delimiter=",")
            assert vals.size == geometry.total_domains
            print(f"buffer test: padded offset used as written ({vals.size} values)")
            return vals
        return original(path, geometry, mode)

    hc.build_island_offset = as_written
    print(f"buffer test: island offset {offset_path(start)}")
    sys.argv = [str(HINDCAST)]
    runpy.run_path(str(HINDCAST), run_name="__main__")


def run(start, preset, label, override=None):
    tag = f"{RUN_DIR_TAG}/{start}_{preset}_{label}"
    env = dict(os.environ, PYTHONIOENCODING="utf-8", PYTHONUNBUFFERED="1", HAT_IGNORE_SETTINGS="1",
               HAT_START_YEAR=str(start), HAT_SOURCE_SINK_PRESET=preset, HAT_SCENARIO=SCENARIO,
               HAT_RELOCATIONS="False", HAT_GROIN_ENABLED="False", HAT_RUN_KIND="experiment",
               HAT_RUN_TAG=tag, HAT_SAVE_MODEL_STATE="False", HAT_OVERWRITE="False")
    if override:
        env["HAT_BE_OVERRIDE"] = override
    logs = EXP_DIR / "logs"
    logs.mkdir(parents=True, exist_ok=True)
    with open(logs / f"{start}_{preset}_{label}.log", "w", encoding="utf-8") as log:
        code = subprocess.run([sys.executable, str(_HERE), "_launch", str(start)], stdout=log,
                              stderr=subprocess.STDOUT, env=env, cwd=PROJECT_ROOT).returncode
    if code:
        raise SystemExit(f"run {tag} failed ({code}); see its log")
    return result_dir(start, preset, label)


def result_dir(start, preset, label):
    base = PROJECT_ROOT / "output" / "raw_runs" / "experiments" / RUN_DIR_TAG / f"{start}_{preset}_{label}"
    hits = list(base.rglob("tables/shoreline_change_rate.csv"))
    return hits[0].parent.parent if hits else None


def cmd_run():
    # One run at a time: three at once ran the machine out of memory (2026-10-04)
    for start in STARTS:
        print(start, "zeroBE", run(start, "zeroBE", "zero"), flush=True)


def main():
    if len(sys.argv) > 2 and sys.argv[1] == "_launch":
        _launch(int(sys.argv[2]))
        return
    {"build": cmd_build, "run": cmd_run}[sys.argv[1]]()


if __name__ == "__main__":
    main()
