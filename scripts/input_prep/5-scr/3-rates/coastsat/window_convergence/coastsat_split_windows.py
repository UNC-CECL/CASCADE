"""
Shoreline position per transect, with the 1996-2024 rate and the two halves cut at a year.

    python scripts/input_prep/5-scr/3-rates/coastsat/window_convergence/coastsat_split_windows.py
    python scripts/input_prep/5-scr/3-rates/coastsat/window_convergence/coastsat_split_windows.py --site ../hatteras-shoreline-windows

Picks eight transects by behaviour (two per group), draws them cut at 2010 and
at four cutoffs, and writes every transect's record for the interactive page.
Reads the fits in 1-rate_profiles/ (run coastsat_window_profiles.py first) and
the raw CoastSat series. Details: the 4-split_windows/ README.

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-10-02
"""

import argparse
import json
import shutil
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import geopandas as gpd
from shapely.geometry import box

# Rule 5: find the root by searching upward.
_REPO = next(_p for _p in Path(__file__).resolve().parents
             if (_p / "pyproject.toml").exists())
sys.path.insert(0, str(_REPO / "scripts"))
from site_layer import hat_observed_rates as obs      # noqa: E402
from site_layer import hat_figure_style as fs         # noqa: E402
from site_layer import hat_map_layers as ml           # noqa: E402
from site_layer.hatteras_site_config import HATTERAS_ANNOTATIONS as ANN  # noqa: E402

sys.path.insert(0, str(_REPO / "scripts" / "input_prep" / "5-scr" / "lib"))
import scr_paths  # noqa: E402,F401
import coastsat_lrr as cl                             # noqa: E402

import matplotlib                                     # noqa: E402
matplotlib.use("Agg")
import matplotlib.pyplot as plt                       # noqa: E402
import matplotlib.dates                               # noqa: E402
from matplotlib.lines import Line2D                   # noqa: E402


# --- CONFIG ------------------------------------------------------------------
REF_START, REF_END = 1996, 2024
CUTOFF = 2010                          # the model's two legs
CUTOFFS = [2005, 2010, 2015, 2020]     # the sensitivity figure's columns
MIN_OBS = 10                           # as the window fits

# The four behaviour groups, in the order they are filled
GROUPS = [
    ("agree", "halves agree"),
    ("disagree", "halves disagree"),
    ("sign_flip", "sign flip"),
    ("step_2021", "2021 step"),
]
PER_GROUP = 2
AGREE_TOL_M_YR = 0.5        # both halves within this of 1996-2024
FLIP_MIN_M_YR = 0.5         # a sign flip needs both halves at least this large
EDGE_DOMAINS = (1, 90)      # the end buffers, never picked
MIN_DOMAIN_GAP = 2          # picks at least this many domains apart

# The 2021 step per transect: 2021 median minus 2019 median, as the domain table
STEP_FROM, STEP_TO = 2019, 2021

STEM_MAIN = "split_windows_{0}".format(CUTOFF)
STEM_MAIN_NO_MAP = "split_windows_{0}_no_map".format(CUTOFF)
STEM_CUTOFFS = "split_windows_cutoffs"
INTERACTIVE_DIR = "interactive"

# The locator map beside the cutoff grid
MAP_CRS = "EPSG:26918"
MAP_WIDTH = 0.9                         # the map column, in panel widths
# ONE PICK PER GROUP IN THE CUTOFF GRID, spread along the island (Hannah, 2026-10-02):
# sign flip GIS 11, 2021 step GIS 33, agree GIS 56, disagree GIS 81
CUTOFF_GRID_DOMAINS = (11, 33, 56, 81)
CUTOFF_GRID_HEIGHT = 6.6
MAP_WIDTH_MAIN = 0.30                   # the same, beside the one wide 2010 column
MAP_WATER, MAP_LAND, MAP_LAND_EDGE = "#e9eff4", "#ede9df", "0.55"
MAP_PLACES = {"Buxton": ANN.town_spans["Buxton"], "Avon": ANN.town_spans["Avon"],
              "Salvo": (ANN.village_lines["Salvo"],) * 2,
              "Waves": (ANN.village_lines["Waves"],) * 2,
              "Rodanthe": (ANN.village_lines["Rodanthe"],) * 2}
# -----------------------------------------------------------------------------


# Every fit of both nested families, one table
def load_fits():
    frames = [pd.read_csv(obs.window_profiles_dir(d, a) / "window_profiles_transects.csv")
              for d, a in (("forward", REF_START), ("backward", REF_END))]
    fits = pd.concat(frames, ignore_index=True)
    return fits.drop_duplicates(["transect_id", "window"])


# The fit for one window, by transect
def window_fit(fits, start, end):
    sub = fits[fits["window"] == "{0}_{1}".format(start, end)]
    return sub.set_index("transect_id")


# `usa_NC_0034_0054` -> its raw series, clipped to the record
def load_series(transect_id):
    site = transect_id.rsplit("_", 1)[0]
    path = obs.COASTSAT_TIMESERIES / "{0}_timeseries".format(site) / "{0}.csv".format(transect_id)
    df = cl.load_timeseries(str(path))
    return cl.filter_dates(df, "{0}-01-01".format(REF_START), "{0}-12-31".format(REF_END))


# The year-by-year median position
def annual_median(series):
    by_year = series.groupby(series["date"].dt.year)["chainage_m"].median()
    return by_year.index.to_numpy(), by_year.to_numpy()


# One row per transect: the three rates, the 2021 step, the nourished flag
def score_transects(fits, series):
    lookup = pd.read_csv(obs.TRANSECT_DOMAINS / "transect_domain_lookup.csv")
    lookup = lookup.dropna(subset=["domain_number"])
    first = window_fit(fits, REF_START, CUTOFF)["lrr_m_yr"]
    second = window_fit(fits, CUTOFF, REF_END)["lrr_m_yr"]
    whole = window_fit(fits, REF_START, REF_END)["lrr_m_yr"]
    prefill = pd.read_csv(obs.COASTSAT_TIMESERIES.parent / "detrended_position"
                          / "step_2021_by_domain_prefill.csv")
    nourished = prefill.set_index("domain_number")["nourished"]
    rows = []
    for r in lookup.sort_values(["domain_number", "transect_id"]).itertuples(index=False):
        yrs, med = annual_median(series[r.transect_id])
        med = dict(zip(yrs, med))
        step = (med[STEP_TO] - med[STEP_FROM]
                if STEP_TO in med and STEP_FROM in med else np.nan)
        rows.append({
            "transect_id": r.transect_id,
            "domain_number": int(r.domain_number),
            "lrr_first_m_yr": first.get(r.transect_id, np.nan),
            "lrr_second_m_yr": second.get(r.transect_id, np.nan),
            "lrr_whole_m_yr": whole.get(r.transect_id, np.nan),
            "step_2021_m": round(step, 2) if step == step else np.nan,
            "nourished": bool(nourished.get(int(r.domain_number), False)),
        })
    t = pd.DataFrame(rows)
    t["halves_gap_m_yr"] = (t["lrr_second_m_yr"] - t["lrr_first_m_yr"]).abs()
    t["worst_half_miss_m_yr"] = np.maximum(
        (t["lrr_first_m_yr"] - t["lrr_whole_m_yr"]).abs(),
        (t["lrr_second_m_yr"] - t["lrr_whole_m_yr"]).abs())
    t["is_sign_flip"] = ((np.sign(t["lrr_first_m_yr"]) != np.sign(t["lrr_second_m_yr"]))
                         & (t["lrr_first_m_yr"].abs() >= FLIP_MIN_M_YR)
                         & (t["lrr_second_m_yr"].abs() >= FLIP_MIN_M_YR))
    return t


# Each group's candidates, most extreme first
def ranked(t, group):
    ok = t.dropna(subset=["lrr_first_m_yr", "lrr_second_m_yr", "lrr_whole_m_yr"])
    ok = ok[~ok["domain_number"].isin(EDGE_DOMAINS)]
    if group == "agree":
        c = ok[ok["worst_half_miss_m_yr"] <= AGREE_TOL_M_YR]
        return c.sort_values("worst_half_miss_m_yr")
    if group == "disagree":
        return ok[~ok["is_sign_flip"]].sort_values("halves_gap_m_yr", ascending=False)
    if group == "sign_flip":
        return ok[ok["is_sign_flip"]].sort_values("halves_gap_m_yr", ascending=False)
    if group == "step_2021":
        c = ok[~ok["nourished"]].dropna(subset=["step_2021_m"])
        return c.sort_values("step_2021_m", ascending=False)
    raise ValueError(group)


# Two per group, greedily, never two within MIN_DOMAIN_GAP domains
def pick(t):
    taken, picks = [], []
    for group, label in GROUPS:
        n = 0
        for rank, r in enumerate(ranked(t, group).itertuples(index=False), start=1):
            if any(abs(r.domain_number - d) < MIN_DOMAIN_GAP for d in taken):
                continue
            taken.append(r.domain_number)
            row = r._asdict()
            row.update(group=group, group_label=label, rank_in_group=rank)
            picks.append(row)
            n += 1
            if n == PER_GROUP:
                break
    return pd.DataFrame(picks)


# A stored fit as the straight line it is, over the dates it saw
def fit_line(fit):
    t0, t1 = pd.Timestamp(fit["first_obs"], tz="UTC"), pd.Timestamp(fit["last_obs"], tz="UTC")
    span = (t1 - t0).total_seconds() / (86400.0 * 365.25)
    return [t0, t1], [fit["intercept_m"], fit["intercept_m"] + fit["lrr_m_yr"] * span]


# The three line colours: whole, first half, second half
LINE_STYLE = {
    "whole": (fs.C["ACCENT"], 2.0, 6),
    "first": (fs.C["EARLY"], 1.6, 7),
    "second": (fs.C["LATE"], 1.6, 7),
}


# One panel: points, the annual median, the three lines, the rates
def draw_panel(ax, series, fits, transect_id, cutoff, fontsize, inline=False):
    ax.scatter(series["date"], series["chainage_m"], s=2.0, color=fs.C["BASE"],
               alpha=0.18, lw=0, zorder=1)
    yrs, med = annual_median(series)
    ax.plot([pd.Timestamp("{0}-07-01".format(y), tz="UTC") for y in yrs], med,
            color=fs.C["INK_MUTED"], lw=0.9, marker="o", ms=2.0, mfc="white",
            mew=0.5, zorder=4)
    ax.axvline(pd.Timestamp("{0}-07-01".format(cutoff), tz="UTC"), color=fs.C["INK_MUTED"],
               lw=0.6, ls=(0, (2, 2)), zorder=2)
    lines = []
    for key, (start, end) in (("whole", (REF_START, REF_END)),
                              ("first", (REF_START, cutoff)),
                              ("second", (cutoff, REF_END))):
        fit = fits.loc[(fits["transect_id"] == transect_id)
                       & (fits["window"] == "{0}_{1}".format(start, end))]
        colour, lw, z = LINE_STYLE[key]
        if fit.empty or fit.iloc[0]["n_obs"] < MIN_OBS:
            lines.append(("{0}–{1}  no fit".format(start, end), colour))
            continue
        fit = fit.iloc[0]
        x, y = fit_line(fit)
        ax.plot(x, y, color=colour, lw=lw, zorder=z, solid_capstyle="round")
        lines.append(("{0}–{1}  {2:+.2f}".format(start, end, fit["lrr_m_yr"]), colour))
    # Headroom so the rates never sit on the data
    lo, hi = ax.get_ylim()
    ax.set_ylim(lo, hi + 0.42 * (hi - lo))
    for k, (text, colour) in enumerate(lines):
        xy = (0.01 + k * 0.16, 0.95) if inline else (0.02, 0.97 - k * 0.105)
        ax.text(xy[0], xy[1], text, transform=ax.transAxes, ha="left",
                va="top", fontsize=fontsize, color=colour, zorder=12,
                path_effects=fs._halo(2.0))
    ax.set_xlim(pd.Timestamp("{0}-01-01".format(REF_START), tz="UTC"),
                pd.Timestamp("{0}-12-31".format(REF_END), tz="UTC"))
    ax.grid(True, axis="y", alpha=0.6)


# The legend both figures share: the data on the top row, the three fits below
def legend_handles(cutoff_text):
    data = [
        Line2D([], [], color="none", marker="o", ms=3.5, mfc=fs.C["BASE"], mec="none",
               label="CoastSat position"),
        Line2D([], [], color=fs.C["INK_MUTED"], lw=0.9, marker="o", ms=2.6, mfc="white",
               mew=0.6, label="annual median"),
        Line2D([], [], color=fs.C["INK_MUTED"], lw=0.8, ls=(0, (2, 2)),
               label="cutoff year (in both windows)"),
    ]
    fits = [
        Line2D([], [], color=fs.C["ACCENT"], lw=2.0,
               label="{0}–{1}, long-term".format(REF_START, REF_END)),
        Line2D([], [], color=fs.C["EARLY"], lw=2.0,
               label="{0}–{1}, first window".format(REF_START, cutoff_text)),
        Line2D([], [], color=fs.C["LATE"], lw=2.0,
               label="{0}–{1}, second window".format(cutoff_text, REF_END)),
    ]
    # Matplotlib fills legend columns first, so interleave to get the two rows
    return [h for pair in zip(data, fits) for h in pair]


# The shared legend, two rows of three under the panels
def add_legend(fig, cutoff_text):
    fig.legend(handles=legend_handles(cutoff_text), loc="outside lower center", ncol=3,
               handlelength=2.4, columnspacing=2.2, handletextpad=0.7)


# The map column and one row per pick, north at the top; returns the panel axes
def map_and_rows(picks, ncol, width_ratios, height=fs.FIG_H_MAX):
    # ROWS RUN NORTH (top) TO SOUTH (bottom), as on the map (Hannah, 2026-10-02)
    picks = picks.sort_values(["domain_number", "transect_id"], ascending=False)
    nrow = len(picks)
    fig = plt.figure(figsize=fs.figsize("double", height=height), layout="constrained")
    gs = fig.add_gridspec(nrow, ncol + 1, width_ratios=width_ratios)
    axes = np.empty((nrow, ncol), dtype=object)
    for r in range(nrow):
        for c in range(ncol):
            axes[r, c] = fig.add_subplot(gs[r, c + 1],
                                         sharex=axes[0, 0] if (r or c) else None,
                                         sharey=axes[0, 0] if (r or c) else None)
            if r < nrow - 1:
                axes[r, c].tick_params(labelbottom=False)
            if c > 0:
                axes[r, c].tick_params(labelleft=False)
    draw_locator(fig.add_subplot(gs[:, 0]), picks, transect_points())
    return fig, axes, picks


# Row label: the map number (1 at Buxton), the group, the transect
def row_label(ax, r, nrow, p):
    ax.set_ylabel("{0} · {1}\nGIS {2} · {3}".format(
        p.map_number, p.group_label, p.domain_number, p.transect_id.replace("usa_NC_", "")),
        fontsize=6.5)


# ONE Y-AXIS FOR EVERY PANEL, so slopes compare between rows (Hannah, 2026-10-02)
def share_position_axis(axes, picks, series, headroom, year_step):
    pos = np.concatenate([series[t]["chainage_m"].to_numpy() for t in picks["transect_id"]])
    lo, hi = 25 * np.floor(pos.min() / 25), 25 * np.ceil(pos.max() / 25)
    axes[0, 0].set_xticks([pd.Timestamp("{0}-01-01".format(y), tz="UTC")
                           for y in range(REF_START, REF_END + 1, year_step)])
    axes[0, 0].xaxis.set_major_formatter(matplotlib.dates.DateFormatter("%Y"))
    for ax in axes.ravel():
        ax.set_ylim(lo, hi + headroom * (hi - lo))
        ax.set_yticks(np.arange(0, hi + 1, 50))
    # One y label for the whole grid, down the right side, centred on the rows
    nrow = axes.shape[0]
    # on a twin, so it never replaces the row label of a one-column figure
    mid = axes[nrow // 2 - 1, -1].twinx()
    mid.set_yticks([])
    mid.set_ylabel("shoreline position (m)", rotation=270, va="bottom")
    mid.yaxis.set_label_coords(1.02 if axes.shape[1] == 1 else 1.04,
                               0.0 if nrow % 2 == 0 else 0.5)


MAP_CAPTION = ("Rows run north (top) to south (bottom); the map at left numbers each "
               "row's transect on Hatteras Island, north up, 1 at Buxton to 8 at Rodanthe.")


# The eight picks, cut at CUTOFF
def draw_main(picks, fits, series, out_dir):
    fs.apply_style()
    fig, axes, picks = map_and_rows(picks, 1, [MAP_WIDTH_MAIN, 1.0])
    nrow = len(picks)
    for r, p in enumerate(picks.itertuples(index=False)):
        ax = axes[r, 0]
        draw_panel(ax, series[p.transect_id], fits, p.transect_id, CUTOFF, 6.0, inline=True)
        row_label(ax, r, nrow, p)
        ax.tick_params(labelsize=6)
    axes[-1, 0].set_xlabel("year")
    share_position_axis(axes, picks, series, 0.22, 4)
    add_legend(fig, str(CUTOFF))
    paths = fs.save(fig, Path(out_dir) / STEM_MAIN, close=True)
    fs.record_caption(paths[0],
        "Shoreline position through time at eight transects, two from each of "
        "four behaviour groups: both halves agree with the long-term rate, the "
        "halves disagree most (same sign), the halves have opposite signs, and "
        "the largest 2021 step outside nourished domains. Grey dots are every "
        "CoastSat position, {0} to {1}, and the open circles their annual "
        "median. Each straight line is an OLS fit to the raw positions in its "
        "window, drawn over the dates it saw: purple {0}–{1}, red the first "
        "window {0}–{2}, blue the second {2}–{1}; {2} belongs to both (calendar "
        "years, inclusive), and the dotted line marks it. Rates in m/yr, "
        "seaward positive. The annual median is a guide for the eye; no line "
        "is fitted to it. Position is CoastSat chainage on an origin that is "
        "arbitrary per transect; every panel shares one y-axis, so slopes "
        "compare between rows. {3} Groups and ranks are in "
        "split_windows_picks.csv.".format(REF_START, REF_END, CUTOFF, MAP_CAPTION))
    return paths[0]


# The eight picks, cut at CUTOFF, by group in a 4 x 2 grid without the map
def draw_main_no_map(picks, fits, series, out_dir):
    fs.apply_style()
    fig, axes = plt.subplots(4, 2, figsize=fs.figsize("double", height=8.8),
                             sharex=True, layout="constrained")
    for i, (ax, p) in enumerate(zip(axes.ravel(), picks.itertuples(index=False))):
        draw_panel(ax, series[p.transect_id], fits, p.transect_id, CUTOFF, 6.6)
        fs._title(ax, i, "{0} · GIS {1} · {2}".format(
            p.group_label, p.domain_number, p.transect_id.replace("usa_NC_", "")))
        if i % 2 == 0:
            ax.set_ylabel("shoreline position (m)")
        if i >= 6:
            ax.set_xlabel("year")
    add_legend(fig, str(CUTOFF))
    paths = fs.save(fig, Path(out_dir) / STEM_MAIN_NO_MAP, close=True)
    fs.record_caption(paths[0],
        "The same as split_windows_{2}.png without the locator map, grouped by "
        "behaviour rather than ordered alongshore. Shoreline position through time at eight transects, two from each of "
        "four behaviour groups: both halves agree with the long-term rate, the "
        "halves disagree most (same sign), the halves have opposite signs, and "
        "the largest 2021 step outside nourished domains. Grey dots are every "
        "CoastSat position, {0} to {1}, and the open circles their annual "
        "median. Each straight line is an OLS fit to the raw positions in its "
        "window, drawn over the dates it saw: purple {0}–{1}, red the first "
        "window {0}–{2}, blue the second {2}–{1}; {2} belongs to both (calendar "
        "years, inclusive), and the dotted line marks it. Rates in m/yr, "
        "seaward positive. The annual median is a guide for the eye; no line "
        "is fitted to it. Position is CoastSat chainage on an origin that is "
        "arbitrary per transect, so only the slopes compare between panels, "
        "and each panel is autoscaled. Groups and ranks are in "
        "split_windows_picks.csv.".format(REF_START, REF_END, CUTOFF))
    return paths[0]


# Every transect's midpoint in UTM, with its domain
def transect_points():
    lookup = pd.read_csv(obs.TRANSECT_DOMAINS / "transect_domain_lookup.csv").dropna(
        subset=["domain_number"])
    mids = transect_midpoints(lookup["transect_id"])
    pts = gpd.GeoDataFrame(
        lookup, crs="EPSG:4326",
        geometry=gpd.points_from_xy([mids[t][0] for t in lookup["transect_id"]],
                                    [mids[t][1] for t in lookup["transect_id"]]))
    pts = pts.to_crs(MAP_CRS)
    pts["e"], pts["n"] = pts.geometry.x, pts.geometry.y
    return pts.set_index("transect_id")


# The island, north up, each pick a numbered marker and the villages named
def draw_locator(ax, picks, pts):
    # Cropped to the ocean shoreline the transects sit on, so the island fills the column
    n0, n1 = pts["n"].min() - 900, pts["n"].max() + 900
    e0, e1 = pts["e"].min() - 2200, pts["e"].max() + 1300
    ax.set_facecolor(MAP_WATER)
    land = gpd.read_file(ml.ISLAND_OUTLINE).to_crs(MAP_CRS)
    land = land.geometry.intersection(box(e0 - 20000, n0 - 20000, e1 + 20000, n1 + 20000))
    gpd.GeoSeries(land[~land.is_empty], crs=MAP_CRS).plot(
        ax=ax, facecolor=MAP_LAND, edgecolor=MAP_LAND_EDGE, lw=0.4, zorder=1)
    ax.set_xlim(e0, e1)
    ax.set_ylim(n0, n1)
    ax.set_aspect("equal", adjustable="datalim")
    ax.set_xticks([]); ax.set_yticks([])
    fs.spines_for_image(ax)

    # The picks on the shoreline, their numbers stacked offshore with leaders
    sel = pts.loc[picks["transect_id"]].copy()
    # Numbered south to north, 1 at the Buxton end, as the domains are
    sel["num"] = picks["map_number"].to_numpy()
    sel = sel.sort_values("n")
    gap = 0.035 * (n1 - n0)
    ys = list(sel["n"])
    for k in range(1, len(ys)):
        ys[k] = max(ys[k], ys[k - 1] + gap)
    shift = max(0.0, ys[-1] - (n1 - gap))
    ys = [y - shift for y in ys]
    e_lab = pts["e"].max() + 750
    for (tid, r), y in zip(sel.iterrows(), ys):
        ax.plot([r["e"], e_lab], [r["n"], y], color=fs.C["INK_MUTED"], lw=0.5, zorder=3)
        ax.plot(r["e"], r["n"], "o", ms=3.0, color=fs.C["INK"], zorder=4)
        ax.text(e_lab, y, str(r["num"]), ha="center", va="center", fontsize=6.5,
                fontweight="bold", color=fs.C["INK"], zorder=6,
                bbox=dict(boxstyle="circle,pad=0.25", facecolor="white",
                          edgecolor=fs.C["INK"], lw=0.6))

    # Villages on the sound side, the south end named
    for name, (d0, d1) in MAP_PLACES.items():
        sub = pts[(pts["domain_number"] >= d0) & (pts["domain_number"] <= d1)]
        ax.text(sub["e"].min() - 450, sub["n"].mean(), name, ha="right", va="center",
                fontsize=5.8, fontstyle="italic", color=fs.C["INK"], zorder=5,
                path_effects=fs._halo(1.5))
    south = pts.loc[pts["n"].idxmin()]
    ax.text(south["e"] + 300, south["n"] - 500, "Cape Point", ha="left", va="top",
            fontsize=5.8, color=fs.C["INK"], zorder=5, path_effects=fs._halo(1.5))
    fs.north_arrow(ax, x=0.16, y=0.88, length=0.03)


# One pick per group, one row each, cut at each of CUTOFFS
def draw_cutoffs(picks, fits, series, out_dir):
    fs.apply_style()
    picks = picks[picks["domain_number"].isin(CUTOFF_GRID_DOMAINS)].copy()
    # Numbered 1-4 within this figure, still from the Buxton end (Hannah, 2026-10-02)
    picks["map_number"] = picks["domain_number"].rank(method="first").astype(int)
    fig, axes, picks = map_and_rows(picks, len(CUTOFFS), [MAP_WIDTH] + [1.0] * len(CUTOFFS),
                                    height=CUTOFF_GRID_HEIGHT)
    nrow = len(picks)
    for r, p in enumerate(picks.itertuples(index=False)):
        for c, cutoff in enumerate(CUTOFFS):
            ax = axes[r, c]
            draw_panel(ax, series[p.transect_id], fits, p.transect_id, cutoff, 6.0)
            if r == 0:
                ax.set_title("cut at {0}".format(cutoff))
            if c == 0:
                row_label(ax, r, nrow, p)
            if r == nrow - 1:
                ax.set_xlabel("year")
            ax.tick_params(labelsize=6)
    share_position_axis(axes, picks, series, 0.42, 8)
    add_legend(fig, "cut")
    paths = fs.save(fig, Path(out_dir) / STEM_CUTOFFS, close=True)
    fs.record_caption(paths[0],
        "Four of the eight transects of split_windows_{2}.png, one per behaviour "
        "group and spread along the island, numbered 1 (Buxton) to 4 (Rodanthe) "
        "within this figure; one row each, with the record cut at {3}. In each panel purple is the {0}–{1} rate, red the "
        "first window from {0} to the cut and blue the second from the cut to "
        "{1}; the cut year belongs to both, and the dotted line marks it. Rates "
        "in m/yr, seaward positive, fitted to every raw CoastSat position in the "
        "window. Every panel shares one y-axis, so slopes compare between rows. "
        "{4} Where the red and blue slopes change with the column, the rate "
        "depends on where the record is cut."
        .format(REF_START, REF_END, CUTOFF, ", ".join(str(c) for c in CUTOFFS),
                MAP_CAPTION))
    return paths[0]


# Midpoint of each NC transect line, lon/lat
def transect_midpoints(ids):
    gj = json.loads((obs.TRANSECT_DOMAINS / "CoastSat_transect_layer.geojson")
                    .read_text(encoding="utf-8"))
    want, out = set(ids), {}
    for f in gj["features"]:
        tid = f["properties"]["id"]
        if tid in want:
            xy = np.asarray(f["geometry"]["coordinates"], dtype=float)
            out[tid] = [round(float(xy[:, 0].mean()), 5), round(float(xy[:, 1].mean()), 5)]
    return out


# The island outline and village anchors for the page's map, in lon/lat
def map_layers_lonlat():
    pts = transect_points()
    win = box(pts["e"].min() - 3000, pts["n"].min() - 3000,
              pts["e"].max() + 3000, pts["n"].max() + 3000)
    land = gpd.read_file(ml.ISLAND_OUTLINE).to_crs(MAP_CRS).geometry.intersection(win)
    land = gpd.GeoSeries(land[~land.is_empty], crs=MAP_CRS).simplify(15).to_crs("EPSG:4326")
    rings = []
    for geom in land:
        for poly in getattr(geom, "geoms", [geom]):
            if poly.geom_type == "Polygon":
                rings.append([[round(x, 5), round(y, 5)] for x, y in poly.exterior.coords])
    places = []
    for name, (d0, d1) in MAP_PLACES.items():
        sub = pts[(pts["domain_number"] >= d0) & (pts["domain_number"] <= d1)]
        ll = gpd.GeoSeries(gpd.points_from_xy([sub["e"].min()], [sub["n"].mean()]),
                           crs=MAP_CRS).to_crs("EPSG:4326").iloc[0]
        places.append({"name": name, "ll": [round(ll.x, 5), round(ll.y, 5)]})
    return rings, places


# Every transect's raw record, for the interactive page
def write_interactive_data(t, picks, series, out_dir):
    mids = transect_midpoints(t["transect_id"])
    group = dict(zip(picks["transect_id"], picks["group"]))
    # Pick numbers as on the figures: 1 at the Buxton end, rising north
    order = picks.sort_values(["domain_number", "transect_id"])["transect_id"]
    number = dict((tid, k) for k, tid in enumerate(order, start=1))
    rows = []
    for r in t.itertuples(index=False):
        s = series[r.transect_id]
        # Years as compute_lrr counts them: days since 1 Jan REF_START over 365.25
        days = (s["date"] - pd.Timestamp("{0}-01-01".format(REF_START), tz="UTC")).dt.total_seconds() / 86400.0
        year = REF_START + days / 365.25
        rows.append({
            "id": r.transect_id.replace("usa_NC_", ""),
            "gis": r.domain_number,
            "ll": mids.get(r.transect_id),
            "g": group.get(r.transect_id),
            "num": number.get(r.transect_id),
            "t": [round(float(v), 4) for v in year],
            "x": [round(float(v), 1) for v in s["chainage_m"]],
        })
    rings, places = map_layers_lonlat()
    d = Path(out_dir) / INTERACTIVE_DIR
    d.mkdir(parents=True, exist_ok=True)
    path = d / "split_windows_data.json"
    path.write_text(json.dumps({
        "ref": [REF_START, REF_END], "min_obs": MIN_OBS, "cutoff": CUTOFF,
        "groups": dict(GROUPS), "outline": rings, "places": places,
        "transects": rows}, separators=(",", ":")),
        encoding="utf-8")
    return path


# Copy the page and its data into the website repo (index.html at its root)
def publish_site(out_dir, site):
    site = Path(site)
    if not (site / ".git").is_dir():
        raise SystemExit("--site {0} is not a git repo; clone or init it first".format(site))
    src = Path(out_dir) / INTERACTIVE_DIR
    # The source page is a fragment (the Claude viewer adds the document around it);
    # a plain web server sends it as-is, so give it its own document and charset
    page = (src / "split_windows_explorer.html").read_text(encoding="utf-8")
    cut = page.index("</style>") + len("</style>")
    (site / "index.html").write_text(
        "<!doctype html>\n<html lang=\"en\">\n<head>\n<meta charset=\"utf-8\">\n"
        "<meta name=\"viewport\" content=\"width=device-width, initial-scale=1\">\n"
        + page[:cut].strip() + "\n<style>body { margin: 0; }</style>\n</head>\n<body>\n"
        + page[cut:].strip() + "\n</body>\n</html>\n", encoding="utf-8")
    shutil.copyfile(src / "split_windows_data.json", site / "split_windows_data.json")
    print("site updated: {0} (commit and push there to publish)".format(site))


# Run: score, pick, draw, export
def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    ap.add_argument("--site", help="website repo to copy the page and data into "
                    "(e.g. ../hatteras-shoreline-windows)")
    ap.add_argument("--site-only", action="store_true",
                    help="only copy the existing page and data to --site, no refit")
    args = ap.parse_args()
    out = obs.split_windows_dir()
    out.mkdir(parents=True, exist_ok=True)
    if args.site_only:
        if not args.site:
            raise SystemExit("--site-only needs --site")
        publish_site(out, args.site)
        return

    # Fits and raw series
    fits = load_fits()
    lookup = pd.read_csv(obs.TRANSECT_DOMAINS / "transect_domain_lookup.csv")
    ids = lookup.dropna(subset=["domain_number"])["transect_id"]
    series = {tid: load_series(tid) for tid in ids}

    # Score every transect and pick two per group
    t = score_transects(fits, series)
    t.to_csv(fs.support_dir(out) / "split_windows_transects.csv", index=False)
    picks = pick(t)
    # Numbered south to north, 1 at the Buxton end, as the domains are
    picks["map_number"] = picks["domain_number"].rank(method="first").astype(int)
    picks.to_csv(out / "split_windows_picks.csv", index=False)
    print(picks[["group", "rank_in_group", "domain_number", "transect_id",
                 "lrr_first_m_yr", "lrr_second_m_yr", "lrr_whole_m_yr",
                 "step_2021_m"]].to_string(index=False))

    # Figures, then the page's data
    print(draw_main(picks, fits, series, out))
    print(draw_main_no_map(picks, fits, series, out))
    print(draw_cutoffs(picks, fits, series, out))
    print(write_interactive_data(t, picks, series, out))
    if args.site:
        publish_site(out, args.site)


if __name__ == "__main__":
    main()
