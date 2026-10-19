#!/usr/bin/env python3
"""
Five averaging windows, the same methods for each of the 1996, 2009 and 2025 target shorelines, drawn like the
10-07 end-line figures but as shoreline POSITION, not net change (the advisor's request, 2026-10-12).

    python target_shoreline_positions.py  ->  figures/target_shoreline_positions_<year>_<frame>.png
                                              for year 1996, 2009, 2025 and the three frames below

Three position frames, one figure per year and frame, alongshore by GIS domain:
  coastsat_baseline   CoastSat's own cross-shore position: distance from each CoastSat transect's
                      origin, averaged per domain over the transects every window of that year
                      includes. The baseline steps between CoastSat sites, and the steps are in
                      every line alike.
  from_1996_calendar  distance seaward of the calendar-1996 mean shoreline, the same fixed line
                      for all three years (Hannah, 2026-10-12), all on the 100 m model transects
                      by the island-offset build's own intersection.
  from_2009_calendar  the same, seaward of the calendar-2009 mean shoreline, the start line of the
                      test (2009 -> 2025), so the 2025 figure reads as the change the test is
                      scored against (Hannah, 2026-10-13).
Lines are 7-domain LOWESS (GIS 1-10 raw); the raw domain means are written to CSV only. No runs.

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-10-13
"""
from __future__ import annotations

import dataclasses
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

HERE = Path(__file__).resolve().parent
REPO = next(p for p in HERE.parents if (p / "pyproject.toml").exists())
MS = REPO / "scripts" / "input_prep" / "5-scr" / "1-observations" / "mean_shoreline"
NC = REPO / "scripts" / "input_prep" / "5-scr" / "3-rates" / "coastsat" / "net_change"
sys.path[:0] = [str(REPO / "scripts"), str(MS), str(NC)]
from coastsat_mean_shoreline_compared import datum_stations  # noqa: E402
from coastsat_net_change import smooth_like_model  # noqa: E402
from cascade_pipeline.annotations import add_geographic_annotations  # noqa: E402
from site_layer import hat_figure_style as fs  # noqa: E402
from site_layer.hat_observed_rates import mean_shoreline_csv  # noqa: E402
from site_layer.hatteras_site_config import HATTERAS_ANNOTATIONS  # noqa: E402

# --- CONFIG ------------------------------------------------------------------
GIS = pd.RangeIndex(1, 91)
DATA_END = "13 Jan 2026"                       # last CoastSat position
# The same five methods applied to each target year Y (Hannah, 2026-10-12), each with one colour
# in every figure; the calendar year (the advisor's proposal) is drawn last, on top. All lines one
# width. All windows built at min 5 positions per transect.
# method -> (first, last) as ISO dates from Y, and the legend name
METHODS = {
    "y_6mo": (lambda y: (f"{y-1}-07-01", f"{y+1}-06-30"),
              "Calendar year {y} ± 6 months: 1 Jul {p} – 30 Jun {n}"),
    "three": (lambda y: (f"{y-1}-01-01", f"{y+1}-12-31"), "Three calendar years: 1 Jan {p} – 31 Dec {n}"),
    "jan_1yr": (lambda y: (f"{y-1}-01-01", f"{y+1}-01-01"), "1 Jan {y} ± 1 yr: 1 Jan {p} – 1 Jan {n}"),
    "jan_6mo": (lambda y: (f"{y-1}-07-01", f"{y}-07-01"), "1 Jan {y} ± 6 months: 1 Jul {p} – 1 Jul {y}"),
    "calendar": (lambda y: (f"{y}-01-01", f"{y}-12-31"), "Calendar year {y}: 1 Jan – 31 Dec {y} (proposed)"),
}
# Legend order: the proposal first, then the alternatives from the shortest window to the longest
LEGEND_ORDER = ("calendar", "jan_6mo", "jan_1yr", "y_6mo", "three")
PALETTE = {"y_6mo": "#2a9d8f", "three": "#8e6c8a", "jan_1yr": "#e9b949", "jan_6mo": "#e07a5f",
           "calendar": "#1d1d1d"}
LW = 1.3
ROLES = {1996: "start of calibration", 2009: "end of calibration, start of test", 2025: "end of test"}
YEARS = {
    y: dict(role=role, windows=[(f(y), lab.format(y=y, p=y - 1, n=y + 1)
                                 + ("*" if y == 2025 and f(y)[1] > "2026-01-13" else ""), m)
                                for m, (f, lab) in METHODS.items()])
    for y, role in ROLES.items()
}
FRAMES = {
    "coastsat_baseline": dict(
        ylabel="Distance seaward of the\nCoastSat transect origin (m)",
        subtitle="Cross-shore shoreline position from the CoastSat transect origin",
        zero=False,
        ylim=(0.0, 120.0),
        how=("CoastSat's own cross-shore position: the distance along each CoastSat transect from its "
             "origin, averaged per 500 m GIS domain over the transects that every window of that year "
             "includes. The CoastSat baseline is built from site segments and steps between some domains; "
             "the steps are the baseline's, not the beach's, and sit in every line alike.")),
    "from_1996_calendar": dict(
        ylabel="Cross-shore position (m)\nrelative to the 1996 mean shoreline (+ seaward)",
        subtitle="Cross-shore shoreline position relative to the calendar-year 1996 mean shoreline",
        zero=True,
        ylim=(-100.0, 100.0),            # the 10-07 net-change figures' axis
        how=("Distance seaward of the calendar-1996 mean shoreline (1 Jan-31 Dec 1996), the same fixed "
             "line in all three panels, so panel (a)'s proposed line is zero by construction: every line "
             "is intersected with the 100 m model transects by the island-offset build's own intersection "
             "(step1_line_to_raw_offsets.intersect) and averaged per 500 m GIS domain.")),
    "from_2009_calendar": dict(
        ylabel="Cross-shore position (m)\nrelative to the 2009 mean shoreline (+ seaward)",
        subtitle="Cross-shore shoreline position relative to the calendar-year 2009 mean shoreline (start of test)",
        zero=True,
        ylim=(-100.0, 100.0),            # the same as the 1996 frame
        how=("Distance seaward of the calendar-2009 mean shoreline (1 Jan-31 Dec 2009), the start line of "
             "the test period, the same fixed line in all three panels, so the 2009 figure's proposed line "
             "is zero by construction and the 2025 figure shows the 2009→2025 change the test is scored "
             "against; subtracting one fixed line shifts every window alike, so the spread between windows "
             "is the same as in the other frames. Every line is intersected with the 100 m model transects "
             "by the island-offset build's own intersection (step1_line_to_raw_offsets.intersect) and "
             "averaged per 500 m GIS domain.")),
}
# -----------------------------------------------------------------------------


def dates(w):
    if isinstance(w[0], int):
        return f"1 Jan {w[0]}–31 Dec {w[1]}"
    a, b = map(pd.Timestamp, w)
    return f"{a.day} {a:%b %Y}–{b.day} {b:%b %Y}"


def key(w):
    return f"{w[0]}_{w[1]}"


# Per domain, CoastSat's mean position from its own transect origins, over the transects every
# window of the year includes
def coastsat_positions():
    out = {}
    for y in YEARS.values():
        tabs = {key(w): pd.read_csv(mean_shoreline_csv(*w)).set_index("transect_id")
                for w, _, _ in y["windows"]}
        common = set.intersection(*(set(t.index[t.included]) for t in tabs.values()))
        for k, t in tabs.items():
            t = t.loc[sorted(common)]
            out[k] = t.groupby("domain_number")["mean_chainage_m"].mean()
    return pd.DataFrame(out).reindex(GIS)


# Per domain, distance seaward of one calendar-year mean shoreline, all on the model transects
def calendar_positions(st, ref_year):
    ref = METHODS["calendar"][0](ref_year)
    # stations grow landward, so the reference's station minus a line's is + seaward of it
    return st.rsub(st[key(ref)], axis=0).reindex(GIS)


# One year's windows in one frame, as the 10-07 end-line figures: one panel, legend below
# fs.mark_offaxis for an axis that need not be symmetric: a triangle at the edge for each value
# beyond (lo, hi), and the values back for the caption
def mark_offaxis(ax, y, lo, hi, color):
    y = np.asarray(y, dtype=float)
    out = np.isfinite(y) & ((y > hi) | (y < lo))
    if not out.any():
        return []
    up = y[out] > hi
    ax.scatter(np.asarray(GIS)[out], np.where(up, hi, lo), s=13, marker="^", zorder=14, clip_on=False,
               c=[color])
    return sorted(((int(g), float(v)) for g, v in zip(np.asarray(GIS)[out], y[out])), key=lambda t: -abs(t[1]))


def figure(pos, n_obs, frame, year, png, palette=None):
    palette = palette or PALETTE
    f, y = FRAMES[frame], YEARS[year]
    fs.apply_style()
    fig, ax = plt.subplots(figsize=(10.0, 5.6), layout="constrained")
    add_geographic_annotations(ax, dataclasses.replace(
        HATTERAS_ANNOTATIONS, groin_label_y=0.03, groin_label_side="left"))
    fs.open_frame(ax)
    ax.grid(axis="y")
    if f["zero"]:
        ax.axhline(0, color=fs.INK_MUTED, lw=0.6, zorder=2)
    lines = {}
    # The 10-07 figures' smoothing (7-domain LOWESS, GIS 1-10 raw); smoothed lines only (Hannah, 2026-10-12)
    ax.set_ylim(*f["ylim"])                     # the same in all three years (Hannah, 2026-10-12)
    off = []
    for w, lab, m in y["windows"]:
        sm = smooth_like_model(pos[key(w)])
        hits = mark_offaxis(ax, sm.values, *f["ylim"], palette[m])
        if hits:
            off.append((lab.split(":")[0], hits))
        lines[m] = ax.plot(GIS, sm.values, color=palette[m], lw=LW,
                           zorder=10 if m == "calendar" else 4, label=lab)[0]
    ax.set_xlim(1, 90)
    ax.set_ylabel(f["ylabel"], fontsize=9)
    ax.set_xlabel(fs.DOMAIN_AXIS_LABEL)
    ax.set_title(f"{year} target shoreline ({y['role']}): sensitivity to the averaging window\n"
                 f"{f['subtitle']}", loc="left", fontsize=10)
    late = any(lab.endswith("*") for _, lab, _ in y["windows"])
    fig.legend([lines[m] for m in LEGEND_ORDER], [lines[m].get_label() for m in LEGEND_ORDER],
               loc="outside lower center", ncol=2, fontsize=7.5, frameon=False, columnspacing=3.0,
               title="Averaging window of the CoastSat mean shoreline"
               + (f"   (* record ends {DATA_END})" if late else ""),
               title_fontsize=8, alignment="left")
    counts = "; ".join(f"{lab.split(':')[0]} {n_obs[key(w)]}" for w, lab, _ in y["windows"])
    star = (f" * The CoastSat record ends {DATA_END}; windows marked * hold data up to then."
            if late else "")
    fs.caption(fig, (
        f"Sensitivity of the {year} target shoreline ({y['role']}) to the averaging window. Each line is "
        "the CoastSat mean shoreline over one window, per 500 m GIS domain, south to north; the same five "
        "windows are applied to each target year (1996, 2009, 2025), each with the same colour in all "
        f"three figures. The calendar-year {year} mean is the advisor's 2026-10-12 proposal. " + f["how"]
        + f" Median satellite positions per transect: {counts}; minimum 5 per transect." + star
        + " Solid lines: 7-domain LOWESS, GIS 1-10 raw (the smoothing of the 2026-10-07 figures and of the "
        "calibration and test targets); the unsmoothed domain means are in the raw CSV. Shaded: communities; hatched: shoal zones; dash-dot: piers; dotted: Buxton "
        f"groin. Per-domain values in target_shoreline_positions_{frame}.csv. The y-axis is the same in "
        f"the 1996, 2009 and 2025 figures ({f['ylim'][0]:g} to {f['ylim'][1]:g} m)." + "".join(
            f" Off the axis, marked with a triangle at the edge: {lab} "
            + ", ".join(f"{v:+.0f} m at GIS {g}" for g, v in hits) + "." for lab, hits in off)))
    return fs.save(fig, png, close=True)


def main():
    n_obs = {key(w): int(pd.read_csv(mean_shoreline_csv(*w)).query("included").n_obs.median())
             for y in YEARS.values() for w, _, _ in y["windows"]}
    st = datum_stations([(key(w), w, None) for y in YEARS.values() for w, _, _ in y["windows"]])
    for frame, pos in (("coastsat_baseline", coastsat_positions()),
                       ("from_1996_calendar", calendar_positions(st, 1996)),
                       ("from_2009_calendar", calendar_positions(st, 2009))):
        for year in YEARS:
            out = figure(pos, n_obs, frame, year,
                         HERE / "figures" / frame / f"target_shoreline_positions_{year}_{frame}.png")
            print(f"-> {out[0].name}")
        pos.apply(smooth_like_model).round(2).to_csv(
            fs.support_dir(HERE / "figures" / frame) / f"target_shoreline_positions_{frame}_lowess7.csv",
            index_label="gis_domain")
        pos.round(2).to_csv(fs.support_dir(HERE / "figures" / frame) / f"target_shoreline_positions_{frame}.csv",
                            index_label="gis_domain")


if __name__ == "__main__":
    main()
