#!/usr/bin/env python3
"""
The groin in the calibration window, as simply as possible: the start, the observed end, and the model year by year.

    python groin_simple_gif.py [--period calibration|test|both]   ->  figures/groin_simple_<period>.gif

The domains around the Buxton groins, in the style of the runner's position_buxton_groin GIF. Fixed: the start
(the zero line) and the observed change over the window (CoastSat net-change target, the one the groin was scored
on). Moving, one frame per model year: the shoreline change since the start of the model without the groin and
with it. Calibration 1996-2009 at strength 0.16, GIS 1-15; test 2009-2025 at 0.10, off from the 2017 fill,
GIS 2-15 (GIS 1, Cape Point, gains about 160 m and would flatten the rest). Matrix runs, full management, edgeBE,
which the 10-14 refit reproduces to 0.0 m. Model values are per domain, unsmoothed.

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-10-16
"""
from __future__ import annotations

# Non-code files beside this script live in the workspace on D: (hatteras_paths.workspace_twin)
import sys as _hat_sys  # noqa: E402
from pathlib import Path as _HatPath  # noqa: E402
_HAT_DIR = next(_p for _p in _HatPath(__file__).resolve().parents if (_p / "hatteras_paths.py").exists())
if str(_HAT_DIR) not in _hat_sys.path:
    _hat_sys.path.insert(0, str(_HAT_DIR))
from hatteras_paths import workspace_twin  # noqa: E402
import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
# Hatteras/ holds the scripts; data/ and output/ live in the workspace on D: (hatteras_paths)
import sys as _hat_sys  # noqa: E402
_HAT_DIR = next(p for p in HERE.parents  if (p / "hatteras_paths.py").exists())
if str(_HAT_DIR) not in _hat_sys.path:
    _hat_sys.path.insert(0, str(_HAT_DIR))
from hatteras_paths import WORKSPACE as REPO  # noqa: E402
sys.path.insert(0, str(_HAT_DIR))
import matplotlib  # noqa: E402
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.animation import PillowWriter  # noqa: E402
from cascade_pipeline.run_registry import find_run_dir_or_expected  # noqa: E402
from site_layer import hat_figure_style as fs  # noqa: E402
from site_layer.hat_observed_rates import net_change_domain_csv  # noqa: E402

# --- CONFIG ------------------------------------------------------------------
RAW_RUNS = REPO / "output" / "raw_runs"
PERIODS = {  # name: (start, end, run name with {groin} for the groin token, strength, groin off from, GIS shown)
    "calibration": (1996, 2009, "HAT_1996-2009_managed_{groin}_edgeBE", 0.16, None, np.arange(1, 16)),
    "test": (2009, 2025, "HAT_2009-2025_managed_{groin}_edgeBE", 0.10, 2017, np.arange(2, 16)),
}
PAD = 15                                   # buffer domains before GIS 1 in the model's arrays
GROIN_FACE = 5.5
OBSERVED = dict(color=fs.C["REF"], lw=2.0, ls=(0, (5, 2)))
NOGROIN = dict(color=fs.C["ADDED"], lw=2.0)     # amber
GROIN = dict(color=fs.C["LATE"], lw=2.2)       # blue; amber and blue stay apart from the red groins and green observed
FPS = 1.5
HOLD = 6
# -----------------------------------------------------------------------------


def model_change(start, end, stem, token, gis):
    run = find_run_dir_or_expected(RAW_RUNS, stem.format(groin=token), (start, end), "edgeBE")
    x = np.load(run / f"{run.name}_shoreline_matrix.npy")
    idx = PAD + gis - 1
    return -(x[:, idx] - x[0, idx])        # the matrix counts landward +; seaward + here


def make(period):
    START, END, stem, strength, off_from, GIS = PERIODS[period]
    off, on = model_change(START, END, stem, "nogroin", GIS), model_change(START, END, stem, "groin", GIS)
    obs = pd.read_csv(net_change_domain_csv(START, END)).set_index("gis_domain")["net_change_lowess7_m"].reindex(GIS)
    years = list(range(START, START + len(on)))
    allv = np.concatenate([off.ravel(), on.ravel(), obs.to_numpy()])
    pad = 0.1 * (allv.max() - allv.min())
    ylim = (allv.min() - pad, allv.max() + pad)

    fig, ax = plt.subplots(figsize=fs.figsize("double", aspect=0.56))
    fig.subplots_adjust(left=0.1, right=0.98, top=0.88, bottom=0.25)

    def draw(k, y):
        ax.clear()
        ax.set_xlim(GIS[0] - 0.5, GIS[-1] + 0.5)
        ax.set_ylim(*ylim)
        ax.axvline(GROIN_FACE, color=fs.C["GROIN"], lw=1.4, ls=":", zorder=2)
        ax.text(GROIN_FACE + 0.12, ylim[1] - 0.04 * (ylim[1] - ylim[0]), "Buxton groins",
                color=fs.C["GROIN"], fontsize=7.5, va="top")
        ax.axhline(0, color=fs.INK_MUTED, lw=1.0, ls=(0, (4, 3)), zorder=2, label=f"{START} shoreline (start)")
        ax.plot(GIS, obs.to_numpy(), zorder=3, label=f"observed {END} (CoastSat)", **OBSERVED)
        ax.plot(GIS, off[k], marker="o", ms=3.5, zorder=4, label="model, no groin", **NOGROIN)
        ax.plot(GIS, on[k], marker="o", ms=3.5, zorder=5, label=f"model, with groin ({strength:.2f})", **GROIN)
        ax.text(0.985, 0.04, str(y), transform=ax.transAxes, ha="right", va="bottom", fontsize=18,
                fontweight="bold", color=fs.INK)
        if off_from is not None and y >= off_from:
            ax.text(0.985, 0.2, f"groin switched off from the {off_from} fill", transform=ax.transAxes,
                    ha="right", va="bottom", fontsize=8, color=fs.C["GROIN"])
        ax.set_xticks(GIS)
        ax.set_xlabel(fs.DOMAIN_AXIS_LABEL)
        ax.set_ylabel(f"Shoreline change since {START} (m)\nseaward +")
        ax.grid(axis="y")
        fs.open_frame(ax)
        ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.17), ncol=4, frameon=False, fontsize=7.5)
        ax.set_title(f"Buxton groins, {period} run {START}-{END}: with and without the groin", fontsize=9,
                     loc="left")

    OUT = workspace_twin(HERE) / "figures" / f"groin_simple_{period}.gif"
    OUT.parent.mkdir(parents=True, exist_ok=True)
    writer = PillowWriter(fps=FPS)
    with writer.saving(fig, str(OUT), dpi=150):
        for k, y in enumerate(years):
            draw(k, y)
            writer.grab_frame()
        for _ in range(HOLD):
            writer.grab_frame()
    plt.close(fig)
    i5, i6 = int(np.flatnonzero(GIS == 5)[0]), int(np.flatnonzero(GIS == 6)[0])
    print(f"wrote {OUT}  ({years[0]}-{years[-1]}; end GIS 5/6: no groin {off[-1][i5]:+.1f}/{off[-1][i6]:+.1f} m, "
          f"groin {on[-1][i5]:+.1f}/{on[-1][i6]:+.1f} m, observed {obs.loc[5]:+.1f}/{obs.loc[6]:+.1f} m)")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--period", choices=("calibration", "test", "both"), default="both")
    which = ap.parse_args().period
    fs.apply_style()
    for period in (PERIODS if which == "both" else (which,)):
        make(period)


if __name__ == "__main__":
    main()
