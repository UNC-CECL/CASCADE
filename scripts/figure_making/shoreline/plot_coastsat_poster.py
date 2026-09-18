#!/usr/bin/env python3
"""
plot_coastsat_poster.py
==============================================================================
The LOESS-SMOOTHED companion to coastsat_calibration_periods.png: the same
CoastSat rates over the same two run periods, smoothed over a 10-domain (5 km)
window so the alongshore pattern reads without the domain-to-domain scatter.

Drawn in exactly the style of the primary, so the two can be laid side by side
and only the smoothing differs.

WHAT CHANGED 2026-09-17
  * THE PERIODS ARE THE CANONICAL CHAIN, 1996 -> 2010 -> 2024, read from
    HATTERAS_PERIODS, and each curve comes from that window's own LRR product.
    It drew 1984-2004 / 2004-2024 before.
  * IT WAS 13 INCHES WIDE, a poster size, saved with a tight bbox so the
    labels rather than figsize() decided the width. At that size its 12 pt
    bold axis labels reduce to about 5 pt on a page. It is 190 mm now, with
    explicit margins and no tight bbox. Render it at a poster width
    deliberately if a poster copy is wanted.
  * THE HOUSE STYLE WAS NEVER APPLIED. The file imported record_caption at the
    BOTTOM and never called apply_style(), so it was not in the project
    typeface at all despite sitting beside figures that are.
  * OFF THE CANVAS: a two-line bold title, the S/N end labels, the place names
    (Buxton / Avon / Tri-Village / Salvo / Waves / Rodanthe), the
    Accretion/Erosion markers and a footnote paragraph. All caption material
    under 9-figures/STYLE.md, and the caption carries it now.
  * THE COLOURS ARE THE VINTAGE PAIR. It used a dark blue and a brown of its
    own; the earlier period is C_1984 red and the later C_1997 blue here, as
    in every other figure that draws two periods.
  * THE LEGEND IS OUT OF THE PANEL. It was boxed and sat on the data.
==============================================================================
"""
import matplotlib
matplotlib.use("Agg")
import numpy as np
import pandas as pd
import warnings

# HOUSE STYLE: one typeface and one palette across every figure in this
# project. See scripts/site_layer/hat_figure_style.py and 9-figures/STYLE.md. The root is
# found by searching upward (ORGANIZATION.md rule 5), so this block is
# independent of whatever this script calls its own repository variable.
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
OUT = _REPO / "output" / "figures" / "shoreline" / "two_periods_10_domains"

PERIOD_STARTS = (1996, 2010)
PERIODS = [(st, HATTERAS_PERIODS[st]["end_year"]) for st in PERIOD_STARTS]
DOMAIN_COL, LRR_COL = "domain_number", "mean_lrr"
DOMAIN_MIN, DOMAIN_MAX = 1, 90
WINDOW_DOMAINS = 10                      # 5.0 km at the 500 m domain spacing
PERIOD_COLOURS = ((C_1984, C_1984_FILL), (C_1997, C_1997_FILL))


def load_period(start, end):
    path = LRR_DIR / f"{start}_{end}" / "domain_lrr_summary.csv"
    if not path.is_file():
        raise SystemExit(f"\nno CoastSat LRR product for {start}-{end} at {path}\n")
    frame = pd.read_csv(path)[[DOMAIN_COL, LRR_COL]].dropna()
    frame = frame[(frame[DOMAIN_COL] >= DOMAIN_MIN) & (frame[DOMAIN_COL] <= DOMAIN_MAX)]
    return frame.sort_values(DOMAIN_COL).reset_index(drop=True)


def main():
    frac = WINDOW_DOMAINS / (DOMAIN_MAX - DOMAIN_MIN + 1)
    km = WINDOW_DOMAINS * 0.5

    fig, ax = plt.subplots(figsize=figsize("double", height=3.6))
    fig.subplots_adjust(left=0.085, right=0.985, bottom=0.235, top=0.96)

    # THE LIMITS GO FIRST. town_bands() skips any span outside the current
    # view and clamps a label to the visible part of its span, so calling it
    # before set_xlim silently dropped Buxton (GIS 7-8): the axes had
    # autoscaled to the shoal spans and 6.5-8.5 fell outside them. Band and
    # label both vanished, with no error (2026-09-17).
    ax.set_xlim(DOMAIN_MIN - 0.5, DOMAIN_MAX + 0.5)

    for name, (lo, hi) in HATTERAS_ANNOTATIONS.shoal_zones.items():
        ax.axvspan(lo - 0.5, hi + 0.5, facecolor=C["ADDED"], alpha=0.13,
                   lw=0, zorder=0.5)
        # a second row, clear of the village names at 0.985: Avon Shoals
        # spans Avon and Wimble Shoals spans Tri-Village, so the two sets
        # of labels overlap in x and must differ in y
        ax.text((lo + hi) / 2, 0.925, name, transform=ax.get_xaxis_transform(),
                ha="center", va="top", fontsize=7, color="#8a620e", zorder=7)
    # NAMED, like the shoals and the structures: if a band is worth
    # drawing it is worth naming (Hannah, 2026-09-17). town_bands puts
    # these at the top of the panel, and structures() already knows to
    # tuck its own labels under them.
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
        "The LOESS-smoothed companion to coastsat_calibration_periods.png: the "
        "same CoastSat domain-mean LRR rates over the same two run periods, "
        f"{PERIODS[0][0]}–{PERIODS[0][1]} in red and {PERIODS[1][0]}–"
        f"{PERIODS[1][1]} in blue, positive seaward, smoothed with a LOESS "
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
