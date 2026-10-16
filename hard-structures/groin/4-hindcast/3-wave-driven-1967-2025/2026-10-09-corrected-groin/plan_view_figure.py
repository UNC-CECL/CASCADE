#!/usr/bin/env python3
"""
The corrected-groin hindcast as maps: the observed and modelled shorelines around the Buxton groins.

    python plan_view_figure.py   ->  figures/corrected_groin_5_plan_view.png

Three true-scale maps (UTM 18N, north up, 1:1) of the 1.5 km window on the groin field, at the end
of the calibration (2008), just before the 2017 fill (2016) and at the end of the test (2024). The
island is filled up to that year's observed CoastSat shoreline (the median position per transect).
Each modelled shoreline is placed per transect: the observed position over the run's start window
plus the model's shoreline change since its start in that transect's domain, along the transect's
own shore-normal. The model has one value per 500 m domain, so its lines step at domain edges. Map
layers come from ../../../1-observations/analysis/5-shoreline_maps/shoreline_plan_view_gif.py, used
unchanged.

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-10-09
"""
from __future__ import annotations

import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402
from matplotlib.patches import Patch  # noqa: E402
from shapely.geometry import box  # noqa: E402

import results_figures as rf  # noqa: E402

sys.path.insert(0, str(rf.cgm.REPO / "hard-structures" / "groin" / "1-observations" / "analysis" / "5-shoreline_maps"))
import shoreline_plan_view_gif as pv  # noqa: E402
from site_layer.hat_figure_style import C, INK_MUTED, MAP_TEXT_DARK, apply_style, figsize, record_caption, save  # noqa: E402

# --- CONFIG ------------------------------------------------------------------
SR = rf.sr
PAD = 15                                  # buffer domains before GIS 1 in the model's arrays
MAPS = ((2008, rf.CAL, "End of calibration"), (2016, rf.TEST, "Just before the 2017 fill"),
        (2024, rf.TEST, "End of the test"))
MODEL_LW = 1.6
OBS_LW = 2.4
# -----------------------------------------------------------------------------


# Each transect's observed position over a period's start window (the model's starting shoreline)
def start_positions(data, period, transects):
    ch = pd.read_csv(pv.CHAINAGE, usecols=["transect_id", "date", "chainage_m", "source"])
    ch = ch[(ch["source"] == "coastsat") & ch["transect_id"].isin(transects)]
    lo, hi = (pd.Timestamp(d) for d in SR.START_WINDOW[period])
    d = pd.to_datetime(ch["date"])
    return ch[(d >= lo) & (d <= hi)].groupby("transect_id")["chainage_m"].median()


# A modelled shoreline in year y: start position plus the model's seaward change in each transect's domain
def model_line(data, path, period, y, start, dirs):
    x = np.load(path)
    k = y - period
    seaward = -(0.5 * (x[k] + x[k + 1]) - x[0])            # mid-year, landward + in the matrix
    rows = []
    for tid, s0 in start.items():
        gis = int(data.dom.loc[tid])
        dv = seaward[PAD + gis - 1]
        o, dx, dy = dirs.loc[tid]
        rows.append((tid, o[0] + (s0 + dv) * dx, o[1] + (s0 + dv) * dy))
    df = pd.DataFrame(rows, columns=["transect_id", "px", "py"])
    return pv.shoreline(df, data.order)


# Each run's mean offset from the observed shoreline per GIS domain in year y, seaward +
def domain_offsets(data, path, period, y, start):
    x = np.load(path)
    k = y - period
    seaward = -(0.5 * (x[k] + x[k + 1]) - x[0])
    ann = data.ann[data.ann["year"] == y].set_index("transect_id")["median"]
    rows = [(int(data.dom.loc[t]), s0 + seaward[PAD + int(data.dom.loc[t]) - 1] - ann[t])
            for t, s0 in start.items() if t in ann.index]
    return pd.DataFrame(rows, columns=["gis", "off"]).groupby("gis")["off"].mean()


def main():
    apply_style()
    data = pv.Data()
    cx, cy = np.vstack(data.groins).mean(axis=0)
    x0, x1, y0, y1 = cx - pv.HALF_M, cx + pv.HALF_M, cy - pv.HALF_M, cy + pv.HALF_M
    ref = data.ref[data.ref["px"].between(x0 - 400, x1 + 400) & data.ref["py"].between(y0 - 400, y1 + 400)]
    transects = ref["transect_id"].tolist()
    origins = data.ann.drop_duplicates("transect_id").set_index("transect_id")[["origin_x", "origin_y", "dir_x", "dir_y"]]
    dirs = pd.Series({t: ((origins.loc[t, "origin_x"], origins.loc[t, "origin_y"]),
                          origins.loc[t, "dir_x"], origins.loc[t, "dir_y"]) for t in transects if t in origins.index})
    dirs = pd.DataFrame(dirs.tolist(), index=dirs.index, columns=["o", "dx", "dy"])
    starts = {p: start_positions(data, p, dirs.index) for p in (rf.CAL, rf.TEST)}
    files = {p: rf.run_files(p) for p in (rf.CAL, rf.TEST)}
    runs = [(label, col, ls) for label, (col, ls) in rf.RUNS.items() if not label.startswith(f"corrected, {rf.JOINT}")]

    fig, axes = plt.subplots(1, 3, figsize=figsize("double", height=3.9), constrained_layout=True)
    gaps, offsets = {}, {}
    for ax, (y, period, title) in zip(axes, MAPS):
        ax.set_facecolor(pv.OCEAN)
        pv.fill(ax, data.lands[y].intersection(box(x0 - 50, y0 - 50, x1 + 50, y1 + 50)), pv.SAND)
        start = starts[period].dropna()
        sx = pd.DataFrame([(t, dirs.loc[t, "o"][0] + s * dirs.loc[t, "dx"], dirs.loc[t, "o"][1] + s * dirs.loc[t, "dy"])
                           for t, s in start.items()], columns=["transect_id", "px", "py"])
        ax.plot(*pv.shoreline(sx, data.order), color=INK_MUTED, lw=0.9, ls=(0, (3, 2)), zorder=3.4)
        for label, col, ls in runs:
            ax.plot(*model_line(data, files[period][label], period, y, start, dirs), color=col, ls=ls, lw=MODEL_LW,
                    zorder=3.6)
        ax.plot(*data.lines[y], color=rf.OBS_C, lw=OBS_LW, zorder=3.5)
        pv.map_furniture(ax, data, x0, x1, y0, y1, 500, "plain", MAP_TEXT_DARK, segments=1)
        # the offsets on the map: model minus observed, GIS 5 and 6, seaward +
        ax.text(0.97, 0.39, "model $-$ observed (m)\n           GIS 5  GIS 6", transform=ax.transAxes, ha="right",
                va="bottom", fontsize=5.8, family="monospace", color=INK_MUTED, zorder=9)
        for j, (label, col, _ls) in enumerate(runs):
            off = domain_offsets(data, files[period][label], period, y, start)
            offsets[(y, label)] = off
            name = {"no groin": "no groin", "pinned groin (in the runner now)": "pinned"}.get(label, "corrected")
            ax.text(0.97, 0.35 - 0.045 * j, f"{name:<10}{off.loc[5]:+5.0f}  {off.loc[6]:+5.0f}", transform=ax.transAxes,
                    ha="right", va="bottom", fontsize=5.8, family="monospace", color=col, zorder=9,
                    fontweight="bold")
        ax.set_title(f"{title}, {y}", fontsize=8.5, loc="left")
        m = SR.model_change(files[period][f"corrected, {rf.CHOSEN} to 2016"])
        o = SR.coastsat_change(period)
        gaps[y] = (m[y - period], float(o.loc[y]) if y in o.index else np.nan)
    handles = [Line2D([], [], color=rf.OBS_C, lw=OBS_LW, label="observed shoreline (CoastSat, that year)"),
               Line2D([], [], color=INK_MUTED, lw=0.9, ls=(0, (3, 2)), label="observed at the run's start")]
    handles += [Line2D([], [], color=col, ls=ls, lw=MODEL_LW, label=label.replace(" to 2016", ", off from 2017"))
                for label, col, ls in runs]
    handles += [Line2D([], [], color=C["GROIN"], lw=2.0, label="Buxton groins"),
                Patch(facecolor=pv.SAND, label="island (to the observed shoreline)"), Patch(facecolor=pv.OCEAN, label="ocean")]
    fig.legend(handles=handles, loc="outside lower center", ncol=4, fontsize=6.5, frameon=False)
    fig.suptitle("Where are the shorelines around the Buxton groins, observed and modelled? (1.5 km windows, true scale)",
                 fontsize=9.5)
    png = rf.FIGS / "corrected_groin_5_plan_view.png"
    save(fig, png, dpi=250, close=True)
    cg_label, pin_label = f"corrected, {rf.CHOSEN} to 2016", "pinned groin (in the runner now)"
    record_caption(png, (
        "**Tests:** what the groin versions look like on a map against the observed shoreline, at the end of the "
        "calibration, just before the 2017 fill and at the end of the test. **How to read:** each map is the same "
        "1.5 km window centred on the Buxton groin field (GIS 5-7), drawn at true scale, north up (UTM 18N); "
        "white lines and labels are the GIS domain bands. Sand fill and the thick teal line: the observed CoastSat "
        "shoreline that year (the median position per transect). Dashed grey: the observed shoreline over the "
        "run's start window, where each modelled line begins (1995-1997 for the calibration run, 2008-2010 for the "
        "test run). Thin lines: the modelled shorelines, each placed per transect as that start position plus "
        "the model's shoreline change in the transect's domain, along the transect's shore-normal. The model has one "
        "value per 500 m domain, so its lines step at the domain edges; that is its resolution, not an error in "
        "the drawing. The small table on each map is each run's mean offset from the observed shoreline in GIS 5 "
        "(south of the groins) and GIS 6 (the groin field), seaward positive. **Shows:** at true scale the runs "
        "differ by tens of metres in a 1.5 km window, so the table carries the comparison. In 2008, without a groin "
        f"GIS 6 sits {-offsets[(2008, 'no groin')].loc[6]:.0f} m landward of the observed shoreline; with the "
        f"corrected groin {-offsets[(2008, cg_label)].loc[6]:.0f} m. In 2024 the pinned groin, still on, holds "
        f"GIS 6 {offsets[(2024, pin_label)].loc[6]:.0f} m seaward of the observed shoreline and GIS 5 "
        f"{-offsets[(2024, pin_label)].loc[5]:.0f} m landward of it, so its gap is too wide; the corrected groin, "
        f"off from the 2017 fill, is {offsets[(2024, cg_label)].loc[6]:+.0f} and {offsets[(2024, cg_label)].loc[5]:+.0f} m. "
        f"The corrected groin's modelled gap change since the 2009 start is {gaps[2016][0]:+.0f} m in 2016 and "
        f"{gaps[2024][0]:+.0f} m in 2024, against {gaps[2016][1]:+.0f} m and {gaps[2024][1]:+.0f} m observed. "
        "In 2024 every run sits 40-50 m seaward of the observed shoreline at GIS 7-8, north of the window's groins: "
        "the 2017 fill's sand stays where it was placed in the model, a separate misfit the groin does not touch."))
    print("wrote", png.relative_to(rf.cgm.REPO))


if __name__ == "__main__":
    main()
