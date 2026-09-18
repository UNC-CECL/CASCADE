#!/usr/bin/env python3
"""
plot_coastsat_calibration_periods_darkerfill.py
==============================================================================
Observed CoastSat shoreline-change rate per domain, one curve per run period.

This is the PRIMARY shoreline figure. CoastSat is the target the model is
graded against, so the DSAS version of the same plot is the independent check
and lives under supporting/ (Hannah, 2026-09-17).

WHAT CHANGED 2026-09-17, and why
  * THE PERIODS WERE THE OLD PAIR. It read the 1984-2004 and 2004-2024 LRR
    products; the canonical chain is 1996 -> 2010 -> 2024 now. That is NOT a
    relabel: the windows have their own LRR products and the rates genuinely
    differ (GIS 1 is -4.16 m/yr over 1984-2004 and +3.23 over 1996-2010), so
    the curves move with the labels. PERIOD_STARTS drives both, and the ends
    come from HATTERAS_PERIODS.
  * THE ANNOTATIONS WERE ITS OWN. Village spans, groin and pier lines were
    re-declared here as literals -- a fourth copy of what the site config
    already holds -- and drawn in this file's own style. They come from
    town_bands() and structures() now, so they match the reach figures
    exactly and cannot disagree with the config.
  * BOTH SHOAL ZONES. Only Wimble was drawn; HATTERAS_ANNOTATIONS.shoal_zones
    has Avon Shoals (GIS 24-39) as well, and both are real (Hannah,
    2026-09-17: "do both").
  * THE PLACE NAMES CAME OFF THE CANVAS. Buxton / Avon / Tri-Village / Salvo /
    Waves / Rodanthe were printed in the data area; the house style puts that
    naming in the caption, and the bands still show where they are.
  * THE DIRECTION MARKERS WENT INTO THE AXIS LABEL. "Accretion" and "Erosion"
    with arrows sat in the panel; the axis says "+ seaward" now, which is the
    same statement in the place the house style keeps it.
  * LEGEND OUT OF THE PANEL, frameless, below -- as the management figures.
==============================================================================
"""
import matplotlib
matplotlib.use("Agg")
import pandas as pd

# HOUSE STYLE: one typeface and one palette across every figure in this
# project. See scripts/site_layer/hat_figure_style.py and figure_making/STYLE.md. The root is
# found by searching upward (ORGANIZATION.md rule 5), so this block is
# independent of whatever this script calls its own repository variable.
import sys as _sys
from pathlib import Path as _P
_sys.path.insert(0, str(next(_q for _q in _P(__file__).resolve().parents
                             if (_q / "pyproject.toml").exists()) / "scripts"))
from site_layer.hat_figure_style import (apply_style, C, C_1984, C_1997, C_1984_FILL,
                              C_1997_FILL, INK, INK_MUTED, GRID_C, figsize,
                              save, record_caption, town_bands, structures,
                              open_frame, DOMAIN_AXIS_LABEL)
apply_style()
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

from site_layer.hatteras_site_config import HATTERAS_PERIODS, HATTERAS_ANNOTATIONS

# Anchored 2026-09-14: absolute into a home directory, or into a tree
# renamed since. Rule 5 of ORGANIZATION.md.
_PATH_REPO = next(_p for _p in _P(__file__).resolve().parents
                  if (_p / "pyproject.toml").exists())
from site_layer.hat_observed_rates import COASTSAT_LRR_ROOT as LRR_DIR  # noqa: E402
from site_layer import hat_figure_style as _hs  # noqa: E402
OUT = _hs.figure_dir("shoreline", "coastsat_calibration_periods")

# The canonical chain. Each end comes from the site config, and the LRR
# product for a window lives under <start>_<end>/, so naming the starts names
# the data as well as the labels.
PERIOD_STARTS = (1996, 2010)
PERIODS = [(st, HATTERAS_PERIODS[st]["end_year"]) for st in PERIOD_STARTS]

DOMAIN_COL, LRR_COL = "domain_number", "mean_lrr"
DOMAIN_MIN, DOMAIN_MAX = 1, 90

# the vintage pair: the earlier period red, the later blue, as everywhere else
PERIOD_COLOURS = ((C_1984, C_1984_FILL), (C_1997, C_1997_FILL))


def load_period(start, end):
    """Domain-mean LRR for one window, from that window's own product."""
    path = LRR_DIR / f"{start}_{end}" / "domain_lrr_summary.csv"
    if not path.is_file():
        raise SystemExit(
            f"\nno CoastSat LRR product for {start}-{end} at {path}.\n"
            f"Windows present: "
            f"{', '.join(sorted(p.name for p in LRR_DIR.iterdir() if p.is_dir()))}\n")
    frame = pd.read_csv(path)
    frame = frame[[DOMAIN_COL, LRR_COL]].dropna()
    frame = frame[(frame[DOMAIN_COL] >= DOMAIN_MIN)
                  & (frame[DOMAIN_COL] <= DOMAIN_MAX)]
    frame = frame.sort_values(DOMAIN_COL).reset_index(drop=True)
    print(f"  {start}-{end}: {len(frame)} domains, "
          f"island mean {frame[LRR_COL].mean():+.2f} m/yr")
    return frame


def main():
    fig, ax = plt.subplots(figsize=figsize("double", height=3.6))
    fig.subplots_adjust(left=0.085, right=0.985, bottom=0.235, top=0.96)

    # SHOAL ZONES, both of them, under everything. They are the one warm
    # accent on this panel; the periods own the red/blue.
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

    # village spans, from the site config, drawn as every alongshore figure
    # draws them; unlabelled, because the caption names them
    # NAMED, like the shoals and the structures: if a band is worth
    # drawing it is worth naming (Hannah, 2026-09-17). town_bands puts
    # these at the top of the panel, and structures() already knows to
    # tuck its own labels under them.
    town_bands(ax, shade="0.93")
    ax.axhline(0, color=INK_MUTED, lw=0.7, ls=(0, (4, 3)), zorder=2)

    frames = []
    for (start, end), (line_c, fill_c) in zip(PERIODS, PERIOD_COLOURS):
        frame = load_period(start, end)
        frames.append(((start, end), frame))
        x, y = frame[DOMAIN_COL].values, frame[LRR_COL].values
        ax.fill_between(x, 0, y, color=fill_c, alpha=0.55, lw=0, zorder=3)
        ax.plot(x, y, color=line_c, lw=1.4, zorder=4,
                label=f"Period {len(frames)}  ({start}–{end})")

    ax.set_xlabel(DOMAIN_AXIS_LABEL)
    ax.set_ylabel("shoreline change rate\n(m/yr, + seaward)")
    ax.xaxis.set_major_locator(ticker.MultipleLocator(10))
    ax.xaxis.set_minor_locator(ticker.MultipleLocator(5))
    ax.grid(axis="y", color=GRID_C, lw=0.5, zorder=1)
    open_frame(ax)

    # structures() measures text against the settled layout, so it goes after
    # the data and after anything that resizes the axes
    structures(ax)

    handles = [Line2D([], [], color=c, lw=1.6,
                      label=f"Period {i}  ({st}–{en})")
               for i, ((st, en), (c, _)) in enumerate(zip(PERIODS, PERIOD_COLOURS), 1)]
    handles.append(Patch(facecolor=C["ADDED"], alpha=0.13,
                         label="shoal influence"))
    handles.append(Patch(facecolor="0.93", label="village / community zone"))
    fig.legend(handles=handles, loc="lower center", bbox_to_anchor=(0.5, 0.005),
               ncol=4, frameon=False, fontsize=8, handlelength=1.8,
               columnspacing=2.0)

    paths = save(fig, OUT, vector=True, close=True)
    shoals = ", ".join(f"{n} GIS {lo}–{hi}"
                       for n, (lo, hi) in HATTERAS_ANNOTATIONS.shoal_zones.items())
    towns = ", ".join(f"{n} GIS {lo}–{hi}"
                      for n, (lo, hi) in HATTERAS_ANNOTATIONS.town_spans.items())
    record_caption(paths[0],
        "Observed shoreline change rate per model domain, CoastSat, one curve "
        f"per run period: {PERIODS[0][0]}–{PERIODS[0][1]} in red and "
        f"{PERIODS[1][0]}–{PERIODS[1][1]} in blue, positive seaward. Each "
        "curve is the mean of the CoastSat transect LRRs falling in that "
        "500 m domain, read from that window's own product under "
        "5-scr/3-rates/coastsat_lrr/; the two windows are separate fits, not one "
        "record split in two. THIS IS THE TARGET the model is graded against "
        "(see dsas_calibration_periods.png under supporting/ for the "
        f"independent DSAS check). Grey bands are the community zones "
        f"({towns}); the warm bands are the shoal-influence zones ({shoals}). "
        "The solid hairline is the Buxton groin field and the dotted ones the "
        "Avon and Rodanthe piers. South to north from left to right; GIS 1 is "
        "Cape Point, GIS 90 Oregon Inlet.")
    print(f"Saved: {paths[0]}")


if __name__ == "__main__":
    main()
