"""
Does the dune growth rate vary alongshore? Each domain's rate fitted straight from the surveyed crests, no model runs.

    python scripts/hatteras_ms/experiments/HAT_dune_growth_alongshore_survey_fit.py
    python scripts/hatteras_ms/experiments/HAT_dune_growth_alongshore_survey_fit.py --widened
    python scripts/hatteras_ms/experiments/HAT_dune_growth_alongshore_survey_fit.py --widened --intervals

Why (2026-10-21, Hannah, option 1 of the alongshore growth-rate question): the
island-wide rate was fitted on seven stable domains by running the model. Away
from them a domain's crest change also holds storm loss, overwash and
bulldozing, so a model fit per domain would soak those up (as the cape rate
did). This takes the advisor's stable-dune idea down to the cell: a cell counts
as undisturbed when it never drops more than DROP_M between consecutive
surveys (about three times the 0.17 m survey-to-survey noise), is surveyed at
the end and is not on the beach or the landward window edge there.

Per undisturbed cell, the rate r that carries its starting crest to its final
surveyed crest under Barrier3D's own yearly step,
    h(k+1) = h(k) + r h(k) (1 - h(k)/Dmax),   h above the berm (1.7 m NAVD88), Dmax 10.5 m,
over the model years between them (start 1 January 1997 -> 2009 lidar on
1 January 2010, 13 steps; start 2010 -> 2019 lidar on 1 January 2020, 10
steps). The start is the model's own starting dunes (dune-topo arrays); the
drop test runs on the surveys only, so the coarse 1996-2005 grids are compared
with each other and not with a 1 m DEM. Per domain: median and middle half of
the cell rates, where at least MIN_CELLS cells are undisturbed.

--widened (Hannah, same day): 1,148 of the 2009-2019 cells ended with their
crest on the landward edge of the dune window, emptying whole domains. Where a
survey's crest sits on that edge, its crest from the window grown 2 cells
(20 m) each side is used instead (the survey table's widened column), in the
drop test and at the end, and the edge exclusion is dropped. Written beside the
first run with a "_widened_at_edge" suffix.

--intervals (Hannah, same day: the calibration window has nine surveys, use
them): logistic growth is a straight line in log-odds, d/dt ln(h/(Dmax-h)) = r,
so per cell r is the summed log-odds gain over its storm-free survey intervals
divided by the years they span. An interval is skipped when the crest drops
more than DROP_M across it, either end is on the beach, unsurveyed or below
MIN_H0_M or within CEILING_GAP_M of Dmax (buildings in the unclassified 2001
returns; log-odds runs away there), or a fill was placed in that domain during it; growth restarts at the
next survey. A cell needs MIN_QUIET_YEARS of kept intervals. Times are the
surveys' flown midpoints, the start dunes dated to the start survey. Continuous
r is within ~1.5% of Barrier3D's yearly-step r at 0.03. Where every interval of
a cell is kept the sum telescopes to start -> end, so the coarse grids'
smoothing cancels except across a skipped interval.

Check: the stable domains (GIS 20-22, 59-62) should come out near the model-fit
0.02 (1996-2009). Domains a fill reached inside the window are flagged; the
roadway managers' dune rebuilding is not in any record and is not removed.

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-10-21
"""
from __future__ import annotations

import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

_HERE = Path(__file__).resolve()
# Hatteras/ holds the scripts; data/ and output/ live in the workspace on D: (hatteras_paths)
import sys as _hat_sys  # noqa: E402
_HAT_DIR = next(_p for _p in _HERE.parents  if (_p / "hatteras_paths.py").exists())
if str(_HAT_DIR) not in _hat_sys.path:
    _hat_sys.path.insert(0, str(_HAT_DIR))
from hatteras_paths import WORKSPACE as PROJECT_ROOT  # noqa: E402
sys.path[:0] = [str(_HAT_DIR), str(_HERE.parent)]
import HAT_dune_crest_vs_surveys as V  # noqa: E402
from site_layer import hat_figure_style as fs  # noqa: E402
from site_layer import hat_topo_version as tv  # noqa: E402
from site_layer.hatteras_site_config import HATTERAS_NOURISHMENT_PROJECTS  # noqa: E402

# --- CONFIG ------------------------------------------------------------------
TAG = "dunes/2026-10-21-alongshore-dune-growth-survey-fit"
EXP_DIR = PROJECT_ROOT / "output" / "raw_runs" / "experiments" / TAG
BERM_NAVD88 = 1.7
DMAX_M = 10.5                     # Dmaxel 12.2 m NAVD88 minus the berm
DROP_M = 0.5                      # largest survey-to-survey drop an undisturbed cell may have
MIN_H0_M = 0.5                    # a cell starting this close to the berm has no logistic rate to speak of
MIN_CELLS = 10                    # undisturbed cells a domain needs for a rate
# start -> (topography product, first model year, last survey used, adopted range, last year of the window)
WINDOWS = {1996: ("1996-start", 1997, 2009, (0.01, 0.03), 2009),
           2009: ("2009-start", 2010, 2019, (0.02, 0.04), 2025)}
STABLE = (20, 21, 22, 59, 60, 61, 62)
CELLS = 50
WIDENED = "--widened" in sys.argv
INTERVALS = "--intervals" in sys.argv
SUFFIX = ("_widened_at_edge" if WIDENED else "") + ("_intervals" if INTERVALS else "")
CEILING_GAP_M = 0.5               # --intervals: a crest this close to Dmax is a building (2001 is unclassified), not a dune
MIN_QUIET_YEARS = 3.0             # --intervals: storm-free survey time a cell needs
START_FLOWN = {1996: "1996-08 to 1996-10", 2009: "2009-08"}     # when the start dunes were flown
FILL_SPAN_YR = 1.25               # a fill's placement counted from 1 January of its year (Buxton 2017 ran to Feb 2018)
# Why a cell has no rate, checked in this order (the figure's lower strips)
REASONS = {"used": ("Rate fitted", "#f1a340"),
           "beach": ("Became beach in a survey", "#fee391"),
           "drop": (f"Lost > {DROP_M} m between surveys (storm)", "#7f3b08"),
           "unsurveyed": ("Not surveyed at the end", "#e0e0e0"),
           "edge": ("End crest at the landward window edge", "#bdbdbd"),
           "low": (f"Started < {MIN_H0_M} m above the berm", "#fdd49e"),
           "out_of_range": ("Rate outside -0.3 to 0.6", "#525252"),
           "short": (f"Under {MIN_QUIET_YEARS:.0f} storm-free survey years", "#7f3b08")}
# Sand palette for the rate panels (Hannah, 2026-10-21: ambers, oranges and yellows for dunes)
SAND = {"median": "#b35806", "running": "#543005", "iqr": "#d8a25e", "adopted": "#fff0b3", "fill": "#d9d9d9"}
# -----------------------------------------------------------------------------


# Barrier3D's yearly logistic step, vectorised over cells
def grow(h0, r, steps):
    h = h0.copy()
    for _ in range(steps):
        h = h + r * h * (1 - h / DMAX_M)
    return h


# The r carrying h0 to h1 in `steps` years, by bisection (growth rises with r); NaN outside [-0.3, 0.6]
def implied_rate(h0, h1, steps):
    lo, hi = np.full_like(h0, -0.3), np.full_like(h0, 0.6)
    ok = (grow(h0, lo, steps) <= h1) & (grow(h0, hi, steps) >= h1)
    for _ in range(60):
        mid = (lo + hi) / 2
        up = grow(h0, mid, steps) < h1
        lo, hi = np.where(up, mid, lo), np.where(up, hi, mid)
    return np.where(ok, (lo + hi) / 2, np.nan)


def start_heights(product):
    return np.concatenate([np.load(tv.array_path("dune", g, product)).astype(float) * 10 for g in range(1, 91)])


def filled(start, end):
    return {g for p in HATTERAS_NOURISHMENT_PROJECTS if start < p.year <= end for g in p.gis_domains}


def cell_rates(start):
    product, first, last, _, end_year = WINDOWS[start]
    s = V.surveys(start).sort_values(["survey", "gis", "cell"])
    if WIDENED:
        at_edge = s.at_landward_edge.fillna(False).astype(bool) & (s.status == "ok")
        s["crest"] = s.crest.where(~at_edge, s.crest_widened_window_m_navd88)
    years = sorted(s.survey.unique())
    grid = {y: g.set_index(["gis", "cell"]) for y, g in s.groupby("survey")}
    idx = pd.MultiIndex.from_product([range(1, 91), range(CELLS)], names=["gis", "cell"])
    crest = np.stack([grid[y].crest.reindex(idx).to_numpy() for y in years])          # surveys x cells, m NAVD88
    status = np.stack([grid[y].status.reindex(idx).to_numpy() for y in years])
    edge = grid[last].at_landward_edge.reindex(idx).fillna(False).to_numpy(bool) & (not WIDENED)

    # Largest drop between consecutive surveyed values, per cell
    worst = np.zeros(crest.shape[1])
    prev = crest[0].copy()
    for row in crest[1:]:
        have = np.isfinite(row) & np.isfinite(prev)
        worst = np.where(have, np.maximum(worst, prev - row), worst)
        prev = np.where(np.isfinite(row), row, prev)
    on_beach = (status == "beach_landward_of_window").any(axis=0)

    h0 = start_heights(product)
    h1 = crest[years.index(last)] - BERM_NAVD88
    steps = int(grid[last].model_state.iloc[0][:4]) - first
    keep = (worst <= DROP_M) & ~on_beach & ~edge & np.isfinite(h1) & (h0 >= MIN_H0_M)
    r = np.where(keep, implied_rate(h0, np.where(np.isfinite(h1), h1, 0.0), steps), np.nan)
    reason = np.select([on_beach, worst > DROP_M, ~np.isfinite(h1), edge, h0 < MIN_H0_M, np.isnan(r)],
                       ["beach", "drop", "unsurveyed", "edge", "low", "out_of_range"], "used")
    return pd.DataFrame(dict(gis=idx.get_level_values(0), cell=idx.get_level_values(1), h0_m=h0, h1_m=h1,
                             worst_drop_m=worst, on_beach=on_beach, landward_edge=edge, undisturbed=keep,
                             r=r, reason=reason)), steps, years


# Decimal year of a flown string's midpoint: "1996-08 to 1996-10", "1998-09"
def flown_mid(text):
    ends = [int(t[:4]) + (int(t[5:7]) - 0.5) / 12 for t in text.split(" to ")]
    return sum(ends) / len(ends)


def logit(h):
    return np.log(h / (DMAX_M - h))


def interval_cell_rates(start):
    product, first, last, _, end_year = WINDOWS[start]
    s = V.surveys(start)
    s = s[s.survey != start].sort_values(["survey", "gis", "cell"])          # the start dunes stand in for it
    if WIDENED:
        at_edge = s.at_landward_edge.fillna(False).astype(bool) & (s.status == "ok")
        s["crest"] = s.crest.where(~at_edge, s.crest_widened_window_m_navd88)
    years = sorted(s.survey.unique())
    grid = {y: g.set_index(["gis", "cell"]) for y, g in s.groupby("survey")}
    idx = pd.MultiIndex.from_product([range(1, 91), range(CELLS)], names=["gis", "cell"])
    h0 = start_heights(product)
    h = np.vstack([h0] + [grid[y].crest.reindex(idx).to_numpy() - BERM_NAVD88 for y in years])
    beach = np.vstack([np.zeros(h0.size, bool)] +
                      [(grid[y].status.reindex(idx) == "beach_landward_of_window").to_numpy() for y in years])
    t = np.array([flown_mid(START_FLOWN[start])] + [flown_mid(grid[y].flown.iloc[0]) for y in years])
    gis = idx.get_level_values(0).to_numpy()
    fills = [(np.isin(gis, p.gis_domains), p.year) for p in HATTERAS_NOURISHMENT_PROJECTS
             if t[0] < p.year + FILL_SPAN_YR and p.year < t[-1]]

    gain, quiet, worst = np.zeros(h0.size), np.zeros(h0.size), np.zeros(h0.size)
    rows = []
    for a in range(len(t) - 1):
        ha, hb = h[a], h[a + 1]
        drop = ha - hb
        worst = np.where(np.isfinite(drop), np.maximum(worst, drop), worst)
        ok = (np.isfinite(ha) & np.isfinite(hb) & ~beach[a] & ~beach[a + 1]
              & (ha >= MIN_H0_M) & (hb >= MIN_H0_M) & (drop <= DROP_M)
              & (ha < DMAX_M - CEILING_GAP_M) & (hb < DMAX_M - CEILING_GAP_M))
        for where, y in fills:
            if t[a] < y + FILL_SPAN_YR and y < t[a + 1]:
                ok &= ~where
        dl = np.where(ok, logit(np.clip(hb, 1e-6, DMAX_M - 1e-6)) - logit(np.clip(ha, 1e-6, DMAX_M - 1e-6)), 0.0)
        gain += dl
        quiet += np.where(ok, t[a + 1] - t[a], 0.0)
        label = f"{'start' if a == 0 else years[a - 1]}-{years[a]}"
        rows.append(dict(interval=label, years=round(t[a + 1] - t[a], 2), cells_kept=int(ok.sum()),
                         cells_dropped=int((np.isfinite(drop) & (drop > DROP_M)).sum()),
                         island_median_r=round(float(np.median(dl[ok] / (t[a + 1] - t[a]))), 4) if ok.any() else np.nan))
    enough = quiet >= MIN_QUIET_YEARS
    r = np.where(enough, gain / np.where(quiet > 0, quiet, 1), np.nan)
    pd.DataFrame(rows).to_csv(EXP_DIR / "tables" / f"interval_island_rates_{start}_start{SUFFIX}.csv", index=False)
    cells = pd.DataFrame(dict(gis=gis, cell=idx.get_level_values(1), h0_m=h0, h1_m=h[-1], worst_drop_m=worst,
                              quiet_years=np.round(quiet, 2), undisturbed=enough, r=r,
                              reason=np.where(enough, "used", "short")))
    return cells, round(t[-1] - t[0], 2), years


def domain_table(start, cells):
    fills = filled(start, WINDOWS[start][4])
    rows = []
    for g, d in cells.groupby("gis"):
        rr = d.r.dropna()
        enough = len(rr) >= MIN_CELLS
        rows.append(dict(gis=g, cells_undisturbed=int(d.undisturbed.sum()), cells_with_rate=len(rr),
                         r_median=round(float(rr.median()), 4) if enough else np.nan,
                         r_p25=round(float(rr.quantile(.25)), 4) if enough else np.nan,
                         r_p75=round(float(rr.quantile(.75)), 4) if enough else np.nan,
                         h0_median_m=round(float(d.h0_m.median()), 2),
                         stable_domain=g in STABLE, fill_in_window=g in fills))
    return pd.DataFrame(rows)


def figure(tables, cells):
    from matplotlib.patches import Patch
    fs.apply_style()
    fig, grid = plt.subplots(2 * len(tables), 1, figsize=fs.figsize("double", height=7.8), sharex=True,
                             layout="constrained", gridspec_kw={"height_ratios": [3, 1.3] * len(tables)})
    axes, strips = grid[0::2], grid[1::2]
    rate_items, cell_items = {}, {}
    for ax, sx, (start, t) in zip(axes, strips, tables.items()):
        _, _, last, adopted, _ = WINDOWS[start]
        ax.axhspan(*adopted, color=SAND["adopted"], zorder=0)
        rate_items["Adopted island-wide range"] = Patch(color=SAND["adopted"])
        fills = t.gis[t.fill_in_window]
        for g in fills:
            for a in (ax, sx):
                a.axvspan(g - 0.5, g + 0.5, facecolor="none", edgecolor=SAND["fill"], hatch="////", lw=0, zorder=0)
        if len(fills):
            rate_items["Fill placed inside the window"] = Patch(facecolor="white", edgecolor=SAND["fill"], hatch="////")
        ok = t.r_median.notna()
        rate_items["Middle half of the cells"] = ax.vlines(t.gis[ok], t.r_p25[ok], t.r_p75[ok], color=SAND["iqr"], lw=1.0)
        rate_items["Domain median"], = ax.plot(t.gis[ok], t.r_median[ok], "o", ms=4, color=SAND["median"], ls="none")
        roll = t.r_median.rolling(7, center=True, min_periods=4).median()
        rate_items["7-domain running median"], = ax.plot(t.gis, roll, color=SAND["running"], lw=1.6)
        ax.axhline(0, color=fs.INK_MUTED, lw=0.6, ls="--")
        ax.set_ylabel(f"Survey-fit r, {start}–{last}")
        how = "storm-free survey intervals" if INTERVALS else f"{start} start dunes → {last} lidar"
        ax.text(0.01, 0.95, f"{how}; adopted {adopted[0]}–{adopted[1]}", transform=ax.transAxes, fontsize=7, va="top")
        ax.grid(axis="y", color=fs.GRID_C, lw=0.5)
        fs.open_frame(ax)
        # No point: fewer than MIN_CELLS cells left; the strip below says what happened to them
        none = t.gis[~ok]
        if len(none):                    # an empty unclipped line breaks the constrained layout
            rate_items[f"No point: fewer than {MIN_CELLS} of {CELLS} cells usable"], = ax.plot(
                none, np.full(len(none), 0.0), "x", color=fs.INK_MUTED, ms=4, mew=1.0, ls="none",
                transform=ax.get_xaxis_transform(), clip_on=False)

        counts = cells[start].groupby(["gis", "reason"]).size().unstack(fill_value=0)
        counts = counts.reindex(index=range(1, 91), columns=list(REASONS), fill_value=0)
        bottom = np.zeros(90)
        for key, (lab, col) in REASONS.items():
            sx.bar(counts.index, counts[key], bottom=bottom, width=0.85, lw=0, color=col)
            if counts[key].sum():
                cell_items[lab] = Patch(color=col)
            bottom += counts[key].to_numpy()
        cell_items[f"{MIN_CELLS} cells: fewest for a point"] = sx.axhline(MIN_CELLS, color=fs.INK, lw=0.8, ls="--")
        sx.set_ylim(0, CELLS)
        sx.set_yticks([0, MIN_CELLS, CELLS])
        sx.set_ylabel("Cells per\ndomain", fontsize=8)
        fs.open_frame(sx)
    strips[-1].set_xlabel(fs.DOMAIN_AXIS_LABEL)
    strips[-1].set_xlim(0.5, 90.5)
    for items, loc, title in ((rate_items, "outside lower left", "Growth rate panels"),
                              (cell_items, "outside lower right", f"Cell strips: what happened to each domain's {CELLS} cells")):
        leg = fig.legend(list(items.values()), list(items), loc=loc, ncol=2, frameon=False, fontsize=7,
                         title=title, title_fontsize=7.5, alignment="left")
        leg.get_title().set_fontweight("bold")
    out = EXP_DIR / "figures" / f"dune_growth_rate_alongshore_survey_fit{SUFFIX}.png"
    fs.save(fig, out, close=True)
    return out


def main():
    (EXP_DIR / "tables").mkdir(parents=True, exist_ok=True)
    tables, summary, allcells = {}, [], {}
    for start in WINDOWS:
        cells, steps, years = (interval_cell_rates if INTERVALS else cell_rates)(start)
        cells.to_csv(EXP_DIR / "tables" / f"cell_rates_{start}_start{SUFFIX}.csv", index=False)
        t = domain_table(start, cells)
        t.to_csv(EXP_DIR / "tables" / f"domain_rates_{start}_start{SUFFIX}.csv", index=False)
        tables[start] = t
        allcells[start] = cells
        r = cells.r.dropna()
        st = cells[cells.gis.isin(STABLE)].r.dropna()
        summary.append(dict(start=start, model_steps=steps, surveys=" ".join(map(str, years)),
                            cells=len(cells), undisturbed=int(cells.undisturbed.sum()), with_rate=len(r),
                            domains_with_rate=int(t.r_median.notna().sum()),
                            island_median_r=round(float(r.median()), 4),
                            island_iqr=f"{r.quantile(.25):.4f} to {r.quantile(.75):.4f}",
                            stable_median_r=round(float(st.median()), 4), stable_cells=len(st),
                            domain_median_spread=f"{t.r_median.quantile(.1):.4f} to {t.r_median.quantile(.9):.4f} (P10-P90)"))
    s = pd.DataFrame(summary)
    s.to_csv(EXP_DIR / "tables" / f"summary{SUFFIX}.csv", index=False)
    with pd.option_context("display.width", 250, "display.max_columns", 30, "display.max_rows", 200):
        print(s.T.to_string())
        for start, t in tables.items():
            print(f"\n{start}")
            print(t.to_string(index=False))
    print(figure(tables, allcells))


if __name__ == "__main__":
    main()
