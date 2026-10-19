#!/usr/bin/env python3
"""
The two figures the model is compared against, as one: each period's end shoreline measured from its
own start line, five averaging windows each (Hannah, 2026-10-13).

    python period_target_positions.py  ->  figures/period_targets/period_target_positions_1996_2009_2025.png

(a) Calibration: the 2009 windows seaward of the calendar-1996 mean shoreline (storms 1997-2009).
(b) Test: the 2025 windows seaward of the calendar-2009 mean shoreline (storms 2010-2025).
The panels are the 2009 figure of the from_1996_calendar frame and the 2025 figure of the
from_2009_calendar frame of target_shoreline_positions.py, with its windows, colours, smoothing and
axis; the helpers are imported from there. No runs.

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-10-13
"""
from __future__ import annotations

import dataclasses

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import pandas as pd  # noqa: E402

from target_shoreline_positions import (  # noqa: E402
    DATA_END, GIS, HERE, LEGEND_ORDER, LW, PALETTE, YEARS, HATTERAS_ANNOTATIONS,
    add_geographic_annotations, calendar_positions, datum_stations, fs, key, mark_offaxis,
    mean_shoreline_csv, smooth_like_model)

# --- CONFIG ------------------------------------------------------------------
# period -> (start year, end year, storms)
PERIODS = {"Calibration": (1996, 2009, "1997-2009"), "Test": (2009, 2025, "2010-2025")}
YLIM = (-100.0, 100.0)                 # the from_<year>_calendar frames' axis
# The legend names each method without its dates, which differ by panel; the dates are in the caption
NAMES = {"calendar": "Calendar year (proposed)", "jan_6mo": "1 Jan ± 6 months", "jan_1yr": "1 Jan ± 1 yr",
         "y_6mo": "Calendar year ± 6 months", "three": "Three calendar years"}
STEM = "period_target_positions_1996_2009_2025"
# -----------------------------------------------------------------------------


def figure(st, n_obs, png):
    fs.apply_style()
    fig, axes = plt.subplots(2, 1, figsize=(10.0, 9.0), sharex=True, layout="constrained")
    lines, off = {}, []
    for i, (ax, (period, (y0, y1, storms))) in enumerate(zip(axes, PERIODS.items())):
        pos = calendar_positions(st, y0)
        add_geographic_annotations(ax, dataclasses.replace(
            HATTERAS_ANNOTATIONS, groin_label_y=0.03, groin_label_side="left"))
        fs.open_frame(ax)
        ax.grid(axis="y")
        ax.axhline(0, color=fs.INK_MUTED, lw=0.6, zorder=2)
        ax.set_ylim(*YLIM)
        for w, lab, m in YEARS[y1]["windows"]:
            sm = smooth_like_model(pos[key(w)])
            hits = mark_offaxis(ax, sm.values, *YLIM, PALETTE[m])
            if hits:
                off.append((f"({fs._letter(i)}) {lab.split(':')[0]}", hits))
            lines[m] = ax.plot(GIS, sm.values, color=PALETTE[m], lw=LW,
                               zorder=10 if m == "calendar" else 4)[0]
        ax.set_ylabel(f"Cross-shore position (m)\nrelative to the {y0} mean shoreline (+ seaward)", fontsize=9)
        fs._title(ax, i, "")
        ax.set_title(f"{period} target: {y1} shoreline relative to the calendar-year {y0} mean shoreline "
                     f"(storms {storms})", loc="center", fontsize=10)
    axes[0].tick_params(labelbottom=True)
    axes[-1].set_xlim(1, 90)
    axes[-1].set_xlabel(fs.DOMAIN_AXIS_LABEL)
    fig.legend([lines[m] for m in LEGEND_ORDER], [NAMES[m] for m in LEGEND_ORDER],
               loc="outside lower center", ncol=5, fontsize=7.5, frameon=False, columnspacing=2.0,
               title="Averaging window of the CoastSat mean shoreline at the end of each period",
               title_fontsize=8)

    def window_dates(y):
        return "; ".join(f"{NAMES[m]}, {lab.split(': ')[1].replace(' (proposed)', '')} ({n_obs[key(w)]})"
                         for w, lab, m in YEARS[y]["windows"])
    fs.caption(fig, (
        "The observed targets the model is compared against, each measured from its own period's start "
        "line: (a) calibration, the 2009 CoastSat mean shoreline seaward of the calendar-1996 mean; (b) test, "
        "the 2025 mean seaward of the calendar-2009 mean. Each line is the end shoreline averaged over one "
        "window; the black calendar-year line is the advisor's 2026-10-12 proposal and the target, the "
        "others show how much the target depends on the window (the start line is the calendar year "
        "throughout). Windows, with median satellite positions per transect: (a) " + window_dates(2009)
        + "; (b) " + window_dates(2025) + f". * The CoastSat record ends {DATA_END}; windows marked * hold "
        "data up to then. Every line is intersected with the 100 m model transects by the island-offset "
        "build's own intersection and averaged per 500 m GIS domain; minimum 5 positions per transect. "
        "Lines: 7-domain LOWESS, GIS 1-10 raw (the smoothing of the calibration and test targets). "
        "Shaded: communities; hatched: shoal zones; dash-dot: piers; dotted: Buxton groin. Same y-axis "
        f"in both panels ({YLIM[0]:g} to {YLIM[1]:g} m). The panels are the 2009 figure of "
        "from_1996_calendar and the 2025 figure of from_2009_calendar; values in their CSVs." + "".join(
            f" Off the axis, marked with a triangle at the edge: {lab} "
            + ", ".join(f"{v:+.0f} m at GIS {g}" for g, v in hits) + "." for lab, hits in off)))
    return fs.save(fig, png, close=True)


def main():
    ws = [w for y in YEARS.values() for w, _, _ in y["windows"]]
    st = datum_stations([(key(w), w, None) for w in ws])
    n_obs = {key(w): int(pd.read_csv(mean_shoreline_csv(*w)).query("included").n_obs.median()) for w in ws}
    out = figure(st, n_obs, HERE / "figures" / "period_targets" / f"{STEM}.png")
    print(f"-> {out[0]}")


if __name__ == "__main__":
    main()
