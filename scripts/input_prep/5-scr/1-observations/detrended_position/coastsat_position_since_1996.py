"""
The 2021 jump in plain shoreline positions: where is the shoreline each year, against 1996?

    python scripts/input_prep/5-scr/1-observations/detrended_position/coastsat_position_since_1996.py

No detrending: each transect's annual median position minus its 1996 median,
plus the raw satellite points at four transects; writes one figure and its
table beside the detrended position. Details: scripts/input_prep/5-scr/1-observations/README.md.

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-10-01
"""

import sys
from pathlib import Path

import numpy as np
import pandas as pd

_REPO = next(_p for _p in Path(__file__).resolve().parents
             if (_p / "pyproject.toml").exists())
sys.path.insert(0, str(_REPO / "scripts"))
from site_layer import hat_observed_rates as obs      # noqa: E402
from site_layer import hat_figure_style as fs         # noqa: E402

sys.path.insert(0, str(_REPO / "scripts" / "input_prep" / "5-scr" / "lib"))
import scr_paths  # noqa: E402,F401
import coastsat_lrr as cl                             # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parent))
from coastsat_detrended_position import (NOURISHED, NOURISHED_DOMAINS,  # noqa: E402
                                         FIRST_YEAR, LAST_YEAR, MIN_YEARS)

import matplotlib                                     # noqa: E402
matplotlib.use("Agg")
import matplotlib.pyplot as plt                       # noqa: E402


# --- CONFIG ------------------------------------------------------------------
# One example transect per reach, never nourished (GIS domains)
REACHES = {"Cape Point–Buxton": (1, 15), "Avon": (16, 31),
           "open reach": (32, 66), "Tri-Village–Rodanthe": (67, 90)}
JUMP = (2020, 2021)
# -----------------------------------------------------------------------------


# Every transect's file and domain
def transects():
    lookup = pd.read_csv(obs.TRANSECT_DOMAINS / "transect_domain_lookup.csv")
    lookup = lookup.dropna(subset=["domain_number"]).sort_values(
        ["domain_number", "transect_id"])
    for r in lookup.itertuples(index=False):
        site = r.transect_id.rsplit("_", 1)[0]
        path = (obs.COASTSAT_TIMESERIES / "{0}_timeseries".format(site)
                / "{0}.csv".format(r.transect_id))
        if path.is_file():
            yield r.transect_id, int(r.domain_number), path


# One transect's observations in the record window
def observations(path):
    return cl.filter_dates(cl.load_timeseries(str(path)),
                           "{0}-01-01".format(FIRST_YEAR),
                           "{0}-12-31".format(LAST_YEAR))


# Annual median position minus the 1996 median, year x transect, metres (positive = seaward)
def build_matrix():
    years = np.arange(FIRST_YEAR, LAST_YEAR + 1)
    cols, doms, paths = {}, {}, {}
    for tid, dom, path in transects():
        df = observations(path)
        med = df.groupby(df["date"].dt.year)["chainage_m"].median().reindex(years)
        if med.notna().sum() < MIN_YEARS or np.isnan(med.loc[FIRST_YEAR]):
            continue
        cols[tid] = (med - med.loc[FIRST_YEAR]).to_numpy()
        doms[tid], paths[tid] = dom, path
    P = pd.DataFrame(cols, index=years)
    P.index.name = "year"
    return P, pd.Series(doms), paths


# Per reach, the never-nourished transect whose jump is closest to the reach median
def pick_examples(P, doms):
    jump = P.loc[JUMP[1]] - P.loc[JUMP[0]]
    picks = {}
    for name, (a, b) in REACHES.items():
        sel = doms.between(a, b) & ~doms.isin(NOURISHED_DOMAINS) & jump.notna()
        j = jump[sel]
        picks[name] = (j - j.median()).abs().idxmin()
    return picks


# Island summary per year, never-nourished transects only
def year_table(P, doms):
    U = P.loc[:, ~doms.reindex(P.columns).isin(NOURISHED_DOMAINS).to_numpy()]
    med = U.median(axis=1)
    out = pd.DataFrame({
        "median_since_1996_m": med.round(2),
        "q25_m": U.quantile(0.25, axis=1).round(2),
        "q75_m": U.quantile(0.75, axis=1).round(2),
        "change_from_previous_year_m": med.diff().round(2),
        "n_transects": U.notna().sum(axis=1),
    })
    return out.reset_index()


def draw(table, P, doms, paths, picks, out_dir):
    fs.apply_style()
    fig = plt.figure(figsize=fs.figsize("double", height=8.2), layout="constrained")
    sub = fig.subfigures(3, 1, height_ratios=[1.0, 1.15, 0.8])
    yrs = table["year"].to_numpy()
    jump_band = dict(color=fs.C["INK_MUTED"], alpha=0.15, lw=0, zorder=1)

    # (a) the island: where the shoreline is each year, against 1996
    ax = sub[0].subplots()
    ax.axvspan(*JUMP, **jump_band)
    ax.fill_between(yrs, table["q25_m"], table["q75_m"], color=fs.C["LATE"],
                    alpha=0.2, lw=0, zorder=2, label="middle half of transects")
    ax.plot(yrs, table["median_since_1996_m"], color=fs.C["LATE"], lw=1.6,
            marker="o", ms=3, zorder=3, label="median transect")
    ax.axhline(0.0, color=fs.C["INK"], lw=0.8, zorder=2)
    for _name, (fy, _ds) in NOURISHED.items():
        ax.axvline(fy, color=fs.C["ADDED"], lw=0.7, ls=(0, (3, 2)), zorder=2)
    ax.text(2014.15, 0.97, "fill", transform=ax.get_xaxis_transform(),
            color=fs.C["ADDED"], fontsize=7, va="top")
    ax.text(2022.15, 0.97, "fills", transform=ax.get_xaxis_transform(),
            color=fs.C["ADDED"], fontsize=7, va="top")
    ax.set_xlim(yrs.min() - 0.5, yrs.max() + 0.5)
    ax.set_ylabel("shoreline position\nvs. 1996 (m)")
    fs._title(ax, 0, "Hatteras shoreline (never-nourished transects), relative to 1996")
    ax.legend(loc="lower left")
    ax.grid(True, axis="y", alpha=0.6)
    ax.annotate("seaward ↑\nlandward ↓", xy=(0.995, 0.05), xycoords="axes fraction",
                ha="right", va="bottom", fontsize=7, color=fs.C["INK_MUTED"])

    # (b) four real transects: every satellite image, and the yearly median
    axes = sub[1].subplots(1, 4, sharex=True)
    for k, (ax, (name, tid)) in enumerate(zip(axes, picks.items())):
        df = observations(paths[tid])
        t = df["date"].dt.year + (df["date"].dt.dayofyear - 0.5) / 365.25
        y = df["chainage_m"] - df[df["date"].dt.year == FIRST_YEAR]["chainage_m"].median()
        ax.axvspan(JUMP[0] + 0.5, JUMP[1] + 0.5, **jump_band)
        ax.plot(t, y, ls="none", marker="o", ms=1.2, color=fs.C["BASE"], alpha=0.35,
                zorder=2)
        ax.plot(yrs + 0.5, P[tid].to_numpy(), color=fs.C["LATE"], lw=1.3, zorder=3)
        ax.axhline(0.0, color=fs.C["INK"], lw=0.6, zorder=2)
        lo, hi = np.nanpercentile(y, [1, 99])
        ax.set_ylim(lo - 10, hi + 10)
        ax.set_xlim(FIRST_YEAR, LAST_YEAR + 1)
        ax.set_xticks([2000, 2010, 2020])
        ax.set_title("{0}\nGIS {1}".format(name, doms[tid]), fontsize=8)
        ax.grid(True, axis="y", alpha=0.6)
        if k == 0:
            ax.set_ylabel("vs. 1996 (m)")
    sub[1].suptitle("(b)  Four real transects: grey = each satellite image, "
                    "blue = yearly median", x=0.01, ha="left", fontsize=9)

    # (c) how 2021 compares with every other year-to-year move
    ax = sub[2].subplots()
    ch = table["change_from_previous_year_m"].to_numpy()
    colours = [fs.C["ACCENT"] if y == JUMP[1] else fs.C["BASE"] for y in yrs]
    ax.bar(yrs, ch, width=0.78, color=colours, lw=0, zorder=3)
    ax.axhline(0.0, color=fs.C["INK"], lw=0.8, zorder=4)
    i = int(np.where(yrs == JUMP[1])[0][0])
    ax.annotate("{0:+.1f} m".format(ch[i]), xy=(yrs[i], ch[i]), xytext=(0, 3),
                textcoords="offset points", ha="center", fontsize=8)
    ax.set_xlim(yrs.min() - 0.5, yrs.max() + 0.5)
    ax.set_ylim(top=ch[i] * 1.18)
    ax.set_ylabel("change from\nprevious year (m)")
    ax.set_xlabel("year")
    fs._title(ax, 2, "Year-to-year move of the median transect")
    ax.grid(True, axis="y", alpha=0.6)

    paths_out = fs.save(fig, Path(out_dir) / "shoreline_position_since_1996", close=True)
    fs.record_caption(paths_out[0],
        "The 2021 jump in plain shoreline positions, no detrending. (a) every "
        "never-nourished CoastSat transect's yearly median position minus its "
        "own 1996 median (positive = seaward): the median transect (line) and "
        "the middle half of transects (band). The median stays within about "
        "10 m of 1996 for 24 years, then moves seaward in one year (shaded, "
        "2020 to 2021) and eases back through 2024. Dashed lines are fill years; no fill transect is in "
        "this panel. (b) four never-nourished transects, one per reach, each "
        "the one whose 2020 to 2021 jump is closest to its reach's median: "
        "every satellite shoreline (grey) and the yearly median (blue). (c) "
        "the median transect's change from the year before; 2021 is the "
        "largest move in the record. Values in shoreline_position_since_1996.csv.")
    return paths_out[0]


def main():
    out_dir = obs.DETRENDED_POSITION
    P, doms, paths = build_matrix()
    table = year_table(P, doms)
    table.to_csv(out_dir / "shoreline_position_since_1996.csv", index=False)
    picks = pick_examples(P, doms)
    print(table.to_string(index=False))
    print(picks)
    print(draw(table, P, doms, paths, picks, out_dir))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
