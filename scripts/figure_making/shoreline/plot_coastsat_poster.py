#!/usr/bin/env python3
"""
The LOWESS-smoothed companion to the CoastSat calibration-periods figure (7 domains, 3.5 km).

    python scripts/figure_making/shoreline/plot_coastsat_poster.py

Same periods and style as the primary, only the smoothing differs. Writes to
output/figures/2-observations/shoreline/. Details: scripts/figure_making/shoreline/README.md.

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-09-30
"""
import matplotlib
matplotlib.use("Agg")
import numpy as np
import pandas as pd
import warnings

# House style (site_layer/hat_figure_style.py), applied at import
import sys as _sys
from pathlib import Path as _P
_sys.path.insert(0, str(next(_q for _q in _P(__file__).resolve().parents
                             if (_q / "pyproject.toml").exists()) / "scripts"))
from site_layer.hat_figure_style import (apply_style, C, C_1984, C_1997, C_1984_FILL,
                              C_1997_FILL, INK_MUTED, GRID_C, figsize, save,
                              record_caption, town_bands, structures,
                              open_frame, DOMAIN_AXIS_LABEL)
apply_style()
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from statsmodels.nonparametric.smoothers_lowess import lowess

from site_layer.hatteras_site_config import HATTERAS_PERIODS, HATTERAS_ANNOTATIONS

warnings.filterwarnings("ignore")

_REPO = next(_p for _p in _P(__file__).resolve().parents
             if (_p / "pyproject.toml").exists())
from site_layer.hat_observed_rates import COASTSAT_LRR_ROOT as LRR_DIR  # noqa: E402
from site_layer import hat_figure_style as _hs  # noqa: E402

# --- CONFIG ------------------------------------------------------------------
PERIOD_STARTS = (1996, 2009)
PERIODS = [(st, HATTERAS_PERIODS[st]["end_year"]) for st in PERIOD_STARTS]
DOMAIN_COL, LRR_COL = "domain_number", "mean_lrr"
DOMAIN_MIN, DOMAIN_MAX = 1, 90
WINDOW_DOMAINS = 7                       # 3.5 km at the 500 m domain spacing; 10 until 2026-09-28
# named for its window: two_periods_10_domains until 2026-09-28
OUT = _hs.figure_dir("observations", "shoreline", f"two_periods_{WINDOW_DOMAINS}_domains")
PERIOD_COLOURS = ((C_1984, C_1984_FILL), (C_1997, C_1997_FILL))
# -----------------------------------------------------------------------------


# Domain-mean LRR for one window, from its own product
def load_period(start, end):
    path = LRR_DIR / f"{start}_{end}" / "domain_lrr_summary.csv"
    if not path.is_file():
        raise SystemExit(f"\nno CoastSat LRR product for {start}-{end} at {path}\n")
    frame = pd.read_csv(path)[[DOMAIN_COL, LRR_COL]].dropna()
    frame = frame[(frame[DOMAIN_COL] >= DOMAIN_MIN) & (frame[DOMAIN_COL] <= DOMAIN_MAX)]
    return frame.sort_values(DOMAIN_COL).reset_index(drop=True)


# Run: load and smooth each period, draw, record the caption
def main():
    frac = WINDOW_DOMAINS / (DOMAIN_MAX - DOMAIN_MIN + 1)
    km = WINDOW_DOMAINS * 0.5

    fig, ax = plt.subplots(figsize=figsize("double", height=3.6))
    fig.subplots_adjust(left=0.085, right=0.985, bottom=0.235, top=0.96)

    # THE LIMITS GO FIRST
    ax.set_xlim(DOMAIN_MIN - 0.5, DOMAIN_MAX + 0.5)

    for name, (lo, hi) in HATTERAS_ANNOTATIONS.shoal_zones.items():
        ax.axvspan(lo - 0.5, hi + 0.5, facecolor=C["ADDED"], alpha=0.13,
                   lw=0, zorder=0.5)
        # A second row, clear of the village names at 0.985
        ax.text((lo + hi) / 2, 0.925, name, transform=ax.get_xaxis_transform(),
                ha="center", va="top", fontsize=7, color="#8a620e", zorder=7)
    # NAMED, like the shoals and the structures
    town_bands(ax, shade="0.93")
    ax.axhline(0, color=INK_MUTED, lw=0.7, ls=(0, (4, 3)), zorder=2)

    for (start, end), (line_c, fill_c) in zip(PERIODS, PERIOD_COLOURS):
        frame = load_period(start, end)
        d = frame[DOMAIN_COL].values.astype(float)
        v = frame[LRR_COL].values.astype(float)
        ok = np.isfinite(v)
        sm = lowess(v[ok], d[ok], frac=frac, return_sorted=True)
        x, y = sm[:, 0], sm[:, 1]
        ax.fill_between(x, 0, y, color=fill_c, alpha=0.55, lw=0, zorder=3)
        ax.plot(x, y, color=line_c, lw=1.6, zorder=4)
        print(f"  {start}-{end}: {int(ok.sum())} domains, "
              f"smoothed island mean {y.mean():+.2f} m/yr")

    ax.set_xlabel(DOMAIN_AXIS_LABEL)
    ax.set_ylabel("shoreline change rate\n(m/yr, + seaward)")
    ax.xaxis.set_major_locator(ticker.MultipleLocator(10))
    ax.xaxis.set_minor_locator(ticker.MultipleLocator(5))
    ax.grid(axis="y", color=GRID_C, lw=0.5, zorder=1)
    open_frame(ax)
    structures(ax)

    handles = [Line2D([], [], color=c, lw=1.8,
                      label=f"Period {i}  ({st}–{en})")
               for i, ((st, en), (c, _)) in enumerate(zip(PERIODS, PERIOD_COLOURS), 1)]
    handles.append(Patch(facecolor=C["ADDED"], alpha=0.13, label="shoal influence"))
    handles.append(Patch(facecolor="0.93", label="village / community zone"))
    fig.legend(handles=handles, loc="lower center", bbox_to_anchor=(0.5, 0.005),
               ncol=4, frameon=False, fontsize=8, handlelength=1.8,
               columnspacing=2.0)

    paths = save(fig, OUT, vector=True, close=True)
    record_caption(paths[0],
        "The LOWESS-smoothed companion to coastsat_calibration_periods.png: the "
        "same CoastSat domain-mean LRR rates over the same two run periods, "
        f"{PERIODS[0][0]}–{PERIODS[0][1]} in red and {PERIODS[1][0]}–"
        f"{PERIODS[1][1]} in blue, positive seaward, smoothed with a LOWESS "
        f"window of {WINDOW_DOMAINS} domains ({km:.1f} km, frac={frac:.3f}). "
        "The smoothing is for reading the alongshore pattern only; the model "
        "is scored against the unsmoothed rates. Grey bands are the community "
        "zones, the warm bands the shoal-influence zones, the solid hairline "
        "the Buxton groin field and the dotted ones the Avon and Rodanthe "
        "piers. South to north from left to right; GIS 1 is Cape Point, GIS "
        "90 Oregon Inlet.")
    print(f"Saved: {paths[0]}")


if __name__ == "__main__":
    main()
