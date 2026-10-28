#!/usr/bin/env python3
"""
A run's modelled dune crests against every lidar survey flown in its window, domain by domain.

    python scripts/real-world-cases/Hatteras/cascade_pipeline/plotting/dune_crest_comparison.py <run_dir>

    from cascade_pipeline.plotting.dune_crest_comparison import plot_dune_crest_comparison

Reads the run's saved model (<run>.npz) and its start's survey crests
(<start>-start/survey-crests/survey_dune_crests.csv, from HAT_survey_dune_crests.py). Each survey
is set against the model state of the 1 January its model_state names. Writes, into the run folder:
figures/dune_crests/dune_crests_gis<a>_<b>.png, one per ten domains (the domain's median crest and
its 25th-75th percentile, model line and band, survey dot and bar); figures/dune_crests/
dune_crest_offset_overview.png, model minus survey median for every domain and survey; and
tables/dune_crest_vs_lidar.csv. A survey cell where the dune window is now beach counts as no dune
(the berm), as in HAT_dune_crest_vs_surveys.py; spikes and unsurveyed cells are left out. Called
by run_figures.py; skipped for a run that saved no model.

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-10-26
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

# The packages live in Hatteras/ (the folder holding hatteras_paths.py), which is not installed
_HAT_DIR = next(_p for _p in Path(__file__).resolve().parents if (_p / "hatteras_paths.py").exists())
if str(_HAT_DIR) not in sys.path:
    sys.path.insert(0, str(_HAT_DIR))

import matplotlib.patheffects as pe  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from matplotlib.colors import BoundaryNorm, ListedColormap  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402
from matplotlib.patches import Patch  # noqa: E402
from matplotlib.ticker import FixedLocator, FormatStrFormatter, MultipleLocator  # noqa: E402
from matplotlib.transforms import ScaledTranslation  # noqa: E402

from cascade_pipeline.run_layout import FIGURES, TABLES, resolve  # noqa: E402
from cascade_pipeline.saved_runs import load_saved_cascade  # noqa: E402
from hatteras_paths import WORKSPACE  # noqa: E402
from site_layer import hat_topo_version as tv  # noqa: E402
from site_layer.hat_figure_style import apply_style, mark_offaxis, record_caption  # noqa: E402
from site_layer.hat_scoring import SMOOTH_DOMAINS, smooth_like_target  # noqa: E402
from site_layer.hatteras_site_config import HATTERAS_ANNOTATIONS, HATTERAS_DOMAINS  # noqa: E402

# --- CONFIG ------------------------------------------------------------------
CHUNK = 10                       # domains per figure
OUT_DIR = "dune_crests"          # under the run's figures/
TABLE = "dune_crest_vs_lidar.csv"
MODEL_COLOR = "#8C8C8C"
SURVEY_COLOR = "#2166AC"
STORMS_MARKED = 3                # the period's largest storms by Rhigh, dotted on every panel
# Overview colour bands (m): the middle one, +/-0.5 m, reads as a match; open-ended beyond the ends
OFFSET_BINS_M = (-4, -3, -2, -1, -0.5, 0.5, 1, 2, 3, 4)
NO_SURVEY_COLOR = "#C8C8C8"
NO_DUNE = "beach_landward_of_window"
DROPPED = ("not_surveyed", "spike")
# Overview v2: shoreline position change over the window drawn on the map, on its own axis
COASTSAT_COLOR = "black"
POSITION_MODEL_COLOR = "#1B7F6B"
# v2 softens the map so the lines read over it
MAP_ALPHA = 0.7
# The offshore shoal brackets above the map, a darker gold than the config's so the labels read
SHOAL_COLOR = "#7A5800"
# Overview storm names, keyed by big_storms()'s label; the November 2006 nor'easter has no name of its own
STORM_NAMES = {"Aug 1999": "Hurricane Dennis", "Sep 2003": "Hurricane Isabel",
               "Nov 2006": "Nor'easter"}
# Surveys footnoted on the overview: unclassified lidar, whose crests may include buildings and vegetation
UNCLASSIFIED_SURVEYS = (2001,)
# v2's shoreline axis half-range (m), fixed so runs compare
POSITION_HALF_M = 70
# -----------------------------------------------------------------------------


# The start's survey crests, with each survey dated at the middle of the months flown (decimal year)
def survey_table(start_year, berm_navd88):
    t = pd.read_csv(tv.product_dir(f"{start_year}-start") / "survey-crests" / "survey_dune_crests.csv")
    t = t[~t.status.isin(DROPPED)].copy()
    t["crest"] = t.crest_m_navd88.where(t.status != NO_DUNE, berm_navd88)
    t = t.dropna(subset=["crest"])
    span = t.flown.str.split(" to ")
    first, last = pd.to_datetime(span.str[0]), pd.to_datetime(span.str[-1])
    t["when"] = (first.dt.year + first.dt.month / 12 + last.dt.year + last.dt.month / 12) / 2 - 1 / 24
    t["state_year"] = pd.to_datetime(t.model_state).dt.year
    return t


# Shoreline position change over the window by GIS domain, seaward +, the model's (its endpoint rate x
# the window's years) and CoastSat's, both smoothed as the scoring target is; None without the CoastSat table
def position_change(run_dir, window):
    f = (tv.INIT_ROOT / "5-scr" / "3-rates" / "coastsat" / "total_change" / f"{window[0]}_{window[1]}"
         / "smoothed" / "tables" / "domain_smoothed.csv")
    if not f.is_file():
        return None
    cs = pd.read_csv(f)
    observed = cs[cs.window_domains == SMOOTH_DOMAINS].set_index("domain_number").observed_m
    rate = pd.read_csv(run_dir / TABLES / "shoreline_change_rate.csv").set_index("gis_domain")
    model = smooth_like_target(rate.change_rate_m_yr * (window[1] - window[0]))
    return model, observed


# The run's crests, m NAVD88: (states, 90 GIS, cells), from its own berm and MHW
def model_crests(cascade):
    dom = HATTERAS_DOMAINS.for_total(len(cascade.barrier3d))   # the run's own buffers
    out = []
    for g in range(dom.first_gis_id, dom.last_gis_id + 1):
        b = cascade.barrier3d[dom.gis_to_pad(g)]
        out.append((np.asarray(b.DuneDomain, float).max(axis=2) + b._BermEl + b._MHW) * 10)
    return np.stack(out, axis=1)


# The period's largest storms: (decimal year, "Mon YYYY") from the run's storm summary
def big_storms(storm_file):
    summary = Path(str(storm_file).replace(".npy", "_summary.csv"))
    if not summary.is_file():
        return []
    s = pd.read_csv(summary, parse_dates=["StartTime"]).nlargest(STORMS_MARKED, "Rhigh")
    return [(t.year + (t.dayofyear - 0.5) / 365.25, t.strftime("%b %Y")) for t in s.StartTime]


# Per GIS and survey: survey and model median on the surveyed cells, and the model minus the survey
def compare(crests, surveys, first_model_year):
    rows = []
    for (g, sv), grp in surveys.groupby(["gis", "survey"]):
        state = int(grp.state_year.iloc[0]) - first_model_year
        if not 0 <= state < crests.shape[0]:
            continue
        m = crests[state, g - 1, grp.cell.values]
        rows.append(dict(gis=g, survey=sv, flown=grp.flown.iloc[0], model_state=grp.model_state.iloc[0],
                         cells=len(grp), survey_median_m=grp.crest.median(), model_median_m=np.median(m),
                         model_minus_survey_m=np.median(m) - grp.crest.median()))
    return pd.DataFrame(rows)


def _year_axis(ax, years):
    ax.set_xlim(years[0] - 1.5, years[-1] + 0.5)
    lo = int(np.ceil((years[0] - 1) / 2) * 2)
    ax.xaxis.set_major_locator(FixedLocator(range(lo, int(years[-1]) + 1, 2)))
    ax.xaxis.set_minor_locator(MultipleLocator(1))
    ax.xaxis.set_major_formatter(FormatStrFormatter("%d"))


# One figure of CHUNK domains
def _chunk_figure(crests, surveys, years, storms, gis_list, ylim, run_name, window, out):
    fig, axs = plt.subplots(2, CHUNK // 2, figsize=(13, 5.8), sharex=True, sharey=True)
    for ax, g in zip(axs.flat, gis_list):
        c = crests[:, g - 1, :]
        lo, mid, hi = np.percentile(c, [25, 50, 75], axis=1)
        ax.fill_between(years, lo, hi, color=MODEL_COLOR, alpha=0.2, lw=0)
        ax.plot(years, mid, color=MODEL_COLOR, lw=1.8)
        for _, grp in surveys[surveys.gis == g].groupby("survey"):
            q = grp.crest.quantile([0.25, 0.5, 0.75]).values
            ax.errorbar(grp.when.iloc[0], q[1], yerr=[[q[1] - q[0]], [q[2] - q[1]]], fmt="o", ms=4,
                        color=SURVEY_COLOR, elinewidth=1, capsize=0, zorder=5)
        for t, label in storms:
            ax.axvline(t, color="0.6", lw=0.8, ls=":")
            if g == gis_list[0]:
                ax.text(t, 0.03, label, rotation=90, ha="right", va="bottom", fontsize=6.5, color="0.4",
                        transform=ax.get_xaxis_transform())
        ax.set_title(f"GIS {g}", loc="left", fontsize=10, fontweight="bold")
    axs[0, 0].set_ylim(*ylim)
    _year_axis(axs[0, 0], years)
    for ax in axs[:, 0]:
        ax.set_ylabel("Dune crest (m NAVD88)")
    for ax in axs[-1]:
        ax.set_xlabel("Year")
    fig.legend(handles=[Line2D([], [], color=MODEL_COLOR, lw=2, label="Model (median; band: 25th–75th pct)"),
                        Line2D([], [], color=SURVEY_COLOR, marker="o", ls="none",
                               label="Lidar survey (median; bar: 25th–75th pct)")],
               loc="lower center", ncol=2, frameon=False)
    fig.suptitle(f"Modelled and surveyed dune crests, GIS {gis_list[0]}–{gis_list[-1]}, "
                 f"{window[0]}–{window[1]}", fontsize=12, fontweight="bold", y=0.99)
    fig.text(0.5, 0.925, run_name, ha="center", va="bottom", fontsize=8, color="0.45")
    fig.tight_layout(rect=(0, 0.05, 1, 0.92))
    fig.savefig(out, dpi=200, facecolor="white")
    plt.close(fig)


# "Aug-Oct" from "1996-08 to 1996-10"
def _months(flown):
    parts = [pd.Timestamp(p) for p in str(flown).split(" to ")]
    names = [p.strftime("%b") for p in parts]
    return names[0] if len(set(names)) == 1 else f"{names[0]}–{names[-1]}"


# Model minus survey for every domain and survey: where and when the dunes leave the lidar. Time runs
# up the page; storms are drawn between the surveys they fell between; villages and the groin on top.
# With position (model, CoastSat), the v2: the shoreline position change drawn over the softened map on a
# right axis centred on zero
def _overview_figure(table, surveys, storms, run_name, window, out, position=None):
    grid = table.pivot(index="survey", columns="gis", values="model_minus_survey_m")
    grid = grid.sort_index().reindex(columns=range(1, 91))
    # The start's own survey is the model's starting topography, a match by construction
    grid = grid[grid.index > window[0]]
    n = len(grid)
    flown = table.drop_duplicates("survey").set_index("survey").flown.reindex(grid.index)
    when = surveys.drop_duplicates("survey").set_index("survey").when.reindex(grid.index).values

    # One colour per band, white in the middle (a match), the two ends for the open-ended overflow
    half = (len(OFFSET_BINS_M) - 2) // 2
    blues = plt.get_cmap("Blues_r")(np.linspace(0.0, 0.75, half + 1))
    reds = plt.get_cmap("Reds")(np.linspace(0.25, 1.0, half + 1))
    cmap = ListedColormap(list(blues[1:]) + ["white"] + list(reds[:-1]))
    cmap.set_under(blues[0])
    cmap.set_over(reds[-1])
    cmap.set_bad(NO_SURVEY_COLOR)
    norm = BoundaryNorm(OFFSET_BINS_M, cmap.N)
    fig, ax = plt.subplots(figsize=(13, 0.38 * n + 3.0))
    im = ax.imshow(np.ma.masked_invalid(grid.values), aspect="auto", cmap=cmap, norm=norm, origin="lower",
                   extent=(0.5, 90.5, -0.5, n - 0.5), interpolation="nearest",
                   alpha=None if position is None else MAP_ALPHA)
    # Thin white rules between surveys and every tenth domain
    for y in np.arange(0.5, n - 1, 1):
        ax.axhline(y, color="white", lw=0.6)
    ax.set_yticks(range(n))
    footnoted = [s for s in grid.index if s in UNCLASSIFIED_SURVEYS]
    ax.set_yticklabels([f"{s}{'†' if s in footnoted else ''} ({_months(f)})" for s, f in flown.items()],
                       fontsize=8.5)
    ax.tick_params(axis="y", length=0)

    # Storms between the surveys they fell between; labelled once the map's top axis is known
    storm_rows = []
    for t, label in storms:
        i = int(np.searchsorted(when, t))
        if 0 < i < n:
            ax.axhline(i - 0.5, color="0.1", lw=1.2, ls=(0, (4, 2)))
            if label in STORM_NAMES:
                label = f"{STORM_NAMES[label]} ({label})"
            storm_rows.append((i - 0.5, label))

    # Villages bracketed just above the map and the offshore shoals a tier above them, offset in inches so
    # they sit the same whatever the figure's height; the groin dotted through the map
    def tier(dy):
        return ax.get_xaxis_transform() + ScaledTranslation(0, dy, fig.dpi_scale_trans)
    for spans, dy, bar_color, text_color, style in (
            (HATTERAS_ANNOTATIONS.town_spans, 0.05, "0.35", "0.25", "normal"),
            (HATTERAS_ANNOTATIONS.shoal_zones, 0.36, SHOAL_COLOR, SHOAL_COLOR, "italic")):
        for name, (lo, hi) in spans.items():
            ax.plot([lo - 0.5, hi + 0.5], [1, 1], color=bar_color, lw=2.5, transform=tier(dy), clip_on=False,
                    solid_capstyle="butt")
            ax.text((lo + hi) / 2, 1, name, transform=tier(dy + 0.04), ha="center", va="bottom", fontsize=8.5,
                    fontweight="bold", color=text_color, style=style)
    for x in HATTERAS_ANNOTATIONS.groins.values():
        ax.axvline(x, color="0.15", lw=1.0, ls=":")

    ax.xaxis.set_major_locator(MultipleLocator(10))
    ax.xaxis.set_minor_locator(MultipleLocator(1))
    ax.set_xlabel("GIS domain (south → north)")
    ax.set_ylabel("Lidar survey (months flown)")
    position_handles = []
    over = ax
    if position is not None:
        model, observed = position
        over = pos = ax.twinx()
        half = POSITION_HALF_M
        pos.set_ylim(-half, half)
        pos.axhline(0, color="0.3", lw=0.6, zorder=1)
        halo = [pe.Stroke(linewidth=3.6, foreground="white"), pe.Normal()]
        for series, color, label in ((observed, COASTSAT_COLOR, "CoastSat"),
                                     (model, POSITION_MODEL_COLOR, "Model")):
            pos.plot(series.index, series.values, color=color, lw=1.8, path_effects=halo, zorder=5)
            mark_offaxis(pos, series.index, series.values, half, color=color, size=20)
            position_handles.append(Line2D([], [], color=color, lw=1.8, label=label))
        pos.set_xlim(ax.get_xlim())
        pos.set_ylabel(f"Shoreline position change, {window[0]}–{window[1]} (m)\n"
                       f"← landward     seaward →")
        pos.tick_params(axis="y", labelsize=8)
    # The storm labels on the topmost axis, so the position lines do not cover them
    for y, label in storm_rows:
        over.text(0.7, y, f" {label} ", transform=ax.transData, ha="left", va="center", fontsize=7.5,
                  color="0.1", zorder=10, bbox=dict(boxstyle="round,pad=0.15", fc="white", ec="none", alpha=0.9))
    # Title and margins in inches, so they sit the same whatever the number of surveys; the run name bottom right
    height = fig.get_figheight()
    width = fig.get_figwidth()
    title = "Modelled minus observed median dune crest elevation by domain and lidar survey"
    if position is not None:
        title += ", with shoreline position change"
    fig.suptitle(f"{title}, {window[0]}–{window[1]}",
                 x=0.5, y=1 - 0.12 / height, ha="center", va="top", fontsize=12, fontweight="bold")
    fig.text(0.99, 0.12 / height, run_name, ha="right", va="bottom", fontsize=8, color="0.45")
    if footnoted:
        fig.text(0.01, 0.12 / height, f"† {', '.join(map(str, footnoted))}: unclassified lidar, so the crests may "
                 f"include buildings and vegetation", ha="left", va="bottom", fontsize=8, color="0.35")
    # v1 keeps v2's right margin, so the two maps sit on the same box and overlay for a slide transition
    fig.subplots_adjust(left=1.45 / width, right=1 - 0.9 / width,
                        bottom=1.55 / height, top=1 - 1.15 / height)
    box = ax.get_position()

    # Under the map, the colour bar on the left and the legend on the right, in titled groups at fixed
    # positions (inches), so v1's groups stay put when v2 adds the shoreline group
    bar_y = 0.75 / height
    cax = fig.add_axes((box.x0, bar_y, 0.45 * box.width, 0.14 / height))
    bar = fig.colorbar(im, cax=cax, orientation="horizontal", ticks=OFFSET_BINS_M, spacing="proportional",
                       extend="both", format="%g")
    bar.ax.tick_params(labelsize=8)
    bar.set_label("← model lower     Model − survey, median dune crest (m)     model higher →", fontsize=9)
    groups = [("Dune crest map", [Patch(fc="white", ec="0.6", lw=0.6, label="within ±0.5 m"),
                                  Patch(fc=NO_SURVEY_COLOR, label="no survey")], 6.7),
              ("Markers", [Line2D([], [], color="0.1", lw=1.2, ls=(0, (4, 2)), label="largest storms of the period"),
                           Line2D([], [], color="0.15", lw=1.0, ls=":", label="Buxton groin")], 8.2)]
    if position_handles:
        groups.append(("Shoreline position change", position_handles, 10.5))
    for name, handles, x in groups:
        fig.legend(handles=handles, title=name, loc="upper left", bbox_to_anchor=(x / width, 1.0 / height),
                   frameon=False, fontsize=8.5, alignment="left", borderaxespad=0,
                   title_fontproperties=dict(size=9, weight="bold"))
    fig.savefig(out, dpi=200, facecolor="white")
    plt.close(fig)


def plot_dune_crest_comparison(run_dir, cascade=None):
    """Draw the run's dune crest figures and table; returns the figure folder, or None if skipped.

    Args:
        run_dir: The run folder, holding <run>_run_metadata.json and <run>.npz.
        cascade: The run's Cascade, if the caller already loaded it.
    """
    run_dir = Path(run_dir)
    run_name = run_dir.name
    if cascade is None:
        npz = resolve(run_dir, "archive", run_name)
        if not Path(npz).is_file():
            print("  dune crests vs lidar  skipped: the run saved no model")
            return None
        cascade = load_saved_cascade(npz)
    meta = json.loads((run_dir / f"{run_name}_run_metadata.json").read_text(encoding="utf-8"))
    period = meta["period"]
    window = (int(period["start_year"]), int(period["end_year"]))
    first_model_year = int(period["first_model_year"])

    apply_style()
    crests = model_crests(cascade)
    b0 = cascade.barrier3d[0]
    berm_navd88 = (b0._BermEl + b0._MHW) * 10
    surveys = survey_table(window[0], berm_navd88)
    years = first_model_year + np.arange(crests.shape[0])
    storms = big_storms(WORKSPACE / "data" / "hatteras_init" / meta["identity"]["storm_file"])
    table = compare(crests, surveys, first_model_year)

    fig_dir = run_dir / FIGURES / OUT_DIR
    fig_dir.mkdir(parents=True, exist_ok=True)
    (run_dir / TABLES).mkdir(exist_ok=True)
    table.round(3).to_csv(run_dir / TABLES / TABLE, index=False)

    # One y axis for every chunk, so the figures compare
    lo = min(np.percentile(crests, 25, axis=2).min(), surveys.crest.quantile(0.01))
    hi = max(np.percentile(crests, 75, axis=2).max(), surveys.crest.quantile(0.99))
    ylim = (np.floor(lo) - 0.2, np.ceil(hi) + 0.2)
    surveyed = ", ".join(str(s) for s in sorted(surveys.survey.unique()))
    for a in range(1, 91, CHUNK):
        gis_list = list(range(a, a + CHUNK))
        out = fig_dir / f"dune_crests_gis{a:02d}_{a + CHUNK - 1:02d}.png"
        _chunk_figure(crests, surveys, years, storms, gis_list, ylim, run_name, window, out)
        record_caption(out, (
            f"Dune crest elevation (m NAVD88), GIS {a}-{a + CHUNK - 1}, {window[0]}-{window[1]}. Model "
            f"(grey): each annual state (1 January) of run {run_name}, the median over the domain's dune "
            f"cells (line) and the 25th-75th percentile (band). Lidar (blue): the crests of each survey "
            f"({surveyed}) on the same cells, median (dot) and 25th-75th percentile (bar), at the middle of "
            f"the months flown; a cell whose dune window is now beach counts as the berm. Dotted: the "
            f"{len(storms)} largest storms of the period by Rhigh. The 1997-2005 surveys are coarse grids "
            f"(3-9 m), so read domain values, not single cells; 2001 is unclassified."))
    out = fig_dir / "dune_crest_offset_overview.png"
    _overview_figure(table, surveys, storms, run_name, window, out)
    record_caption(out, (
        f"Model minus survey median dune crest (m) for every GIS domain (x) and lidar survey (y, earliest "
        f"at the bottom), {window[0]}-{window[1]}, run {run_name}: red, the model's dunes are higher; blue, "
        f"lower; white, within ±0.5 m; grey, no survey. Each survey is set against the model state of the "
        f"1 January its model_state names; the {window[0]} survey, the model's starting topography, is left "
        f"out. † unclassified lidar ({', '.join(map(str, UNCLASSIFIED_SURVEYS))}), whose crests may include "
        f"buildings and vegetation. Dashed: the {len(storms)} largest storms of the period by Rhigh, named, "
        f"between the surveys they fell between. Bars above: the villages and, a tier higher, the offshore shoals; dotted: the Buxton groin. The "
        f"colour bands are {', '.join(f'{b:g}' for b in OFFSET_BINS_M)} m, open-ended beyond."))
    position = position_change(run_dir, window)
    if position is None:
        print(f"  dune crest overview v2 skipped: no CoastSat total change for {window[0]}_{window[1]}")
    else:
        out = fig_dir / "dune_crest_offset_overview_v2.png"
        _overview_figure(table, surveys, storms, run_name, window, out, position)
        record_caption(out, (
            f"As dune_crest_offset_overview.png, with the map at {MAP_ALPHA:g} opacity and the shoreline position "
            f"change over {window[0]}-{window[1]} drawn over it on the right axis, zero at its middle, seaward positive: "
            f"CoastSat (black), the mean position over the last calendar year minus the first, and the model "
            f"(green), its endpoint rate x {window[1] - window[0]} years, both smoothed as the scoring target "
            f"is ({SMOOTH_DOMAINS}-domain LOWESS, GIS 1-10 raw). The axis is fixed at ±{POSITION_HALF_M} m so "
            f"runs compare; a value beyond it is marked with a triangle at the edge."))
    print(f"  dune crests vs lidar  {fig_dir}")
    return fig_dir


def main():
    import argparse
    import matplotlib
    matplotlib.use("Agg")
    ap = argparse.ArgumentParser(description=__doc__.strip().splitlines()[0])
    ap.add_argument("run_dir", type=Path, nargs="+", help="run folder(s)")
    for rd in ap.parse_args().run_dir:
        plot_dune_crest_comparison(rd)


if __name__ == "__main__":
    main()
