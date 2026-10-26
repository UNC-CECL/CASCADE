#!/usr/bin/env python3
"""
The spread groin against the runner's groin in the calibration window, as simply as possible: the start, the
observed end, and the model year by year.

    python spread_simple_gif.py [--variant spread_0_23_13|spread_12_13_16] [--period calibration|test] [--current]
        ->  figures/<variant>_simple_<period>.gif (+ last frame .png)

The style of ../2026-10-14-window-refit/groin_simple_gif.py. Fixed: the start (the zero line) and the observed change
over 1996-2009 (CoastSat net-change target). Moving, one frame per model year: the shoreline change since 1996
without the groin, with the runner's groin (GIS 6 takes all of its loss across the groin face) and with the spread
groin (shares of it on GIS 6, 7, 8: 0, 2/3, 1/3, the spread groin, or 1/2, 1/3, 1/6), each at its best gap
strength. Model values are per domain,
unsmoothed.

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-10-20
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
sys.path.insert(0, str(HERE))
import spread_relaxation as sp  # noqa: E402
import tapered_fillet as tf  # noqa: E402  (on the path through spread_relaxation)
import matplotlib  # noqa: E402
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.animation import PillowWriter  # noqa: E402
from site_layer import hat_figure_style as fs  # noqa: E402
from site_layer.hat_observed_rates import net_change_domain_csv  # noqa: E402

# --- CONFIG ------------------------------------------------------------------
# period: (start, end, GIS shown, runner's strength, groin off from). The test window shows GIS 2-15 (GIS 1, Cape
# Point, gains about 160 m and would flatten the rest), its controls are this study's own (current end rates), and
# its spread runs at the calibration strength, unchanged
WINDOWS = {"calibration": (1996, 2009, np.arange(1, 16), 0.16, None),
           "test": (2009, 2025, np.arange(2, 16), 0.10, 2017)}
SPREADS = {  # name (by its shares on GIS 6, 7, 8): (run folder, strength, legend text, title text)
    "spread_12_13_16": ("spread", 0.13, "loss spread over GIS 6-8 (1/2, 1/3, 1/6)", "spread over GIS 6-8 (1/2, 1/3, 1/6)"),
    "spread_0_23_13": ("wE", 0.10, "spread groin (0, 2/3, 1/3)", "the spread groin (0, 2/3, 1/3)"),
}
PAD = 15                                   # buffer domains before GIS 1 in the model's arrays
GROIN_FACE = 5.5
# By role: the observed reference strongest, the no-groin baseline recessive, the two groins a colour-blind-safe
# blue/orange pair (Okabe-Ito) with the new one heavier; red stays the groins' own
OBSERVED = dict(color="#1a1a1a", lw=2.4, ls=(0, (5, 2)))
NOGROIN = dict(color="#a6a6a6", lw=1.6)        # light grey
GROIN = dict(color="#0072B2", lw=2.0)          # blue
SPREAD = dict(color="#E69F00", lw=2.8)         # orange
FPS = 1.5
HOLD = 6
# -----------------------------------------------------------------------------


def model_change(rd, gis):
    x = np.load(next(rd.glob("*_shoreline_matrix.npy")))
    idx = PAD + gis - 1
    return -(x[:, idx] - x[0, idx])        # the matrix counts landward +; seaward + here


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--variant", choices=tuple(SPREADS), default="spread_0_23_13")
    ap.add_argument("--period", choices=tuple(WINDOWS), default="calibration")
    ap.add_argument("--current", action="store_true",
                    help="calibration from the 2026-10-21 rerun on today's runner (its own controls), as the test window")
    a = ap.parse_args()
    sp.CURRENT = a.current
    variant, period = a.variant, a.period
    run_folder, S_SPREAD, legend, title = SPREADS[variant]
    START, END, GIS, S_RUNNER, off_from = WINDOWS[period]
    fs.apply_style()
    if period == "calibration" and not a.current:
        off = model_change(tf.run_dir(START, None, None), GIS)
        on = model_change(tf.run_dir(START, "face", S_RUNNER), GIS)
    else:
        off = model_change(sp.run_dir(0.0, "nogroin", START), GIS)
        on = model_change(sp.run_dir(S_RUNNER, "face", START), GIS)
    spread = model_change(sp.run_dir(S_SPREAD, run_folder, START), GIS)
    obs = pd.read_csv(net_change_domain_csv(START, END)).set_index("gis_domain")["net_change_lowess7_m"].reindex(GIS)
    years = list(range(START, START + len(on)))
    allv = np.concatenate([off.ravel(), on.ravel(), spread.ravel(), obs.to_numpy()])
    pad = 0.1 * (allv.max() - allv.min())
    ylim = (allv.min() - pad, allv.max() + pad)

    fig, ax = plt.subplots(figsize=fs.figsize("double", aspect=0.56))
    fig.subplots_adjust(left=0.1, right=0.98, top=0.88, bottom=0.27)

    def draw(k, y):
        ax.clear()
        ax.set_xlim(GIS[0] - 0.5, GIS[-1] + 0.5)
        ax.set_ylim(*ylim)
        ax.axvline(GROIN_FACE, color=fs.C["GROIN"], lw=1.4, ls=":", zorder=2)
        ax.text(GROIN_FACE + 0.12, ylim[1] - 0.04 * (ylim[1] - ylim[0]), "Buxton groins",
                color=fs.C["GROIN"], fontsize=7.5, va="top")
        ax.axhline(0, color=fs.INK_MUTED, lw=0.9, ls=(0, (1, 2)), zorder=2, label=f"{START} shoreline (start)")
        ax.plot(GIS, obs.to_numpy(), zorder=7, label=f"observed {END} (CoastSat)", **OBSERVED)
        ax.plot(GIS, off[k], marker="o", ms=3.5, zorder=4, label="model, no groin", **NOGROIN)
        ax.plot(GIS, on[k], marker="o", ms=3.5, zorder=5, label=f"model, groin without the spread, strength {S_RUNNER:.2f}",
                **GROIN)
        ax.plot(GIS, spread[k], marker="o", ms=3.5, zorder=6,
                label=f"model, {legend}, strength {S_SPREAD:.2f}", **SPREAD)
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
        ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.17), ncol=3, frameon=False, fontsize=7.5)
        ax.set_title(f"Buxton groins, {period} run {START}-{END}: the groin without the spread, or {title}",
                     fontsize=9, loc="left")

    out = workspace_twin(HERE) / "figures" / f"{variant}_simple_{period}{'_current' if a.current else ''}.gif"
    out.parent.mkdir(parents=True, exist_ok=True)
    writer = PillowWriter(fps=FPS)
    with writer.saving(fig, str(out), dpi=150):
        for k, y in enumerate(years):
            draw(k, y)
            writer.grab_frame()
        for _ in range(HOLD):
            writer.grab_frame()
    fig.savefig(out.with_suffix(".png"), dpi=200)       # the last frame, for a document
    plt.close(fig)
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
