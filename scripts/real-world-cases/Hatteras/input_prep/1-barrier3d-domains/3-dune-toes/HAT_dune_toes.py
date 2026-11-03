"""
Dune toes and beach widths per 10 m dune cell, from a start's lidar surveys: the dune-collision inputs.

    python scripts/real-world-cases/Hatteras/input_prep/1-barrier3d-domains/3-dune-toes/HAT_dune_toes.py --start 1996-start
    python scripts/real-world-cases/Hatteras/input_prep/1-barrier3d-domains/3-dune-toes/HAT_dune_toes.py --start 2009-start
    python scripts/real-world-cases/Hatteras/input_prep/1-barrier3d-domains/3-dune-toes/HAT_dune_toes.py --figures

The dune-collision module (cascade/dune_collision.py) lowers a dune the water
reaches but does not overtop, down to its toe, and less on a wide beach. It needs
each cell's toe and beach width. Step 3 of the stage, after 2-survey-crests:
both are taken from the same profiles HAT_survey_dune_crests.py reads crests from (the start's own arrays, then every
survey of its window, on the start's frame and dune windows): the crest is the
highest point of the dune window, the beach start the first cell above the beach
threshold, the toe the point of the beach-to-crest profile furthest below the
straight line between them (the beach start on a straight or convex face), and
the beach width the dry beach, the beach start to the toe (0 where the crest stands within one cell of the
beach start). Until 2026-11-03 the beach width ran to the crest, the dune face included (kept as
beach_to_crest_m): it made a broad dune face read as a wide beach, and the toe is where the swash meets the dune.

Each cell then gets one value, fixed for the run (2026-10-29: kept simple):
the median toe over the surveys whose toe is above TOE_MIN_M_NAVD88 and at least
TOE_BELOW_CREST_M below that survey's crest, and the median beach width over the
surveys; only survey cells with crest status "ok" in survey_dune_crests.csv count
(spikes out). A cell with no value takes its domain's median; a domain with none
stays empty (NaN), and the module leaves it alone. Surveys in TOE_SURVEYS_EXCLUDED
are left out of both: the 1996 start's 1998 grid puts 38% of its toes on the wet
beach (median 1.7 m, some below 0 m NAVD88).

Ported from dune-collision/4-hindcast/1-calibration-1996-2009/HAT_dune_collision.py (toes, cell_toes, the
per-domain fill in _launch), where the 1996-start values were calibrated; this
script reproduces that experiment's dune_toes.csv for the 1996 start.

Reads <start>/survey-crests/survey_dune_crests.csv (step 2) for the spike
status. Writes data/hatteras_init/1-barrier3d-domains/<start>/dune-toes/, beside
dune-topo/ and dune-growth-rate/ (the other dune inputs; survey-crests/ is what
the model's dunes are scored against, so the inputs are kept out of it):
survey_dune_toes.csv (one row per survey, domain and cell) and
dune_toe_cells.csv (one row per domain and cell: the values the runner uses).
--figures draws, from the saved tables of every start, into
output/figures/3-model-inputs/1-domains/dune_toes/: dune_toes_alongshore.png
(toe, crest above toe and beach width per domain), dune_toes_survey_spread.png
(each survey's toe against the fixed value) and dune_toes_method.png (how a toe is
read off one profile and what the module does with it, on one typical cell, EXAMPLE).

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-10-29
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

_HERE = Path(__file__).resolve()
_HAT_DIR = next(_p for _p in _HERE.parents if (_p / "hatteras_paths.py").exists())
for _path in (_HAT_DIR, _HERE.parents[1] / "2-survey-crests"):
    if str(_path) not in sys.path:
        sys.path.insert(0, str(_path))
import HAT_survey_dune_crests as S  # noqa: E402

# --- CONFIG ------------------------------------------------------------------
MHW_NAVD88_M = 0.36              # the survey profiles are m MHW; the tables are m NAVD88
MIN_BEACH_CELLS = 2              # a toe needs a beach start at least this many cells seaward of the crest
TOE_MIN_M_NAVD88 = 1.0           # above the beach start (0.86 m NAVD88)
TOE_BELOW_CREST_M = 0.5          # and at least this far below the crest
FOREDUNE_DROP_M = 1.0            # a peak with a dip this deep behind it is the first dune: the toe is picked under
                                 # a chord to it, not to the crest (2026-11-03, after the 1997 imagery check)
TOE_SURVEYS_EXCLUDED = {"1996-start": (1998,), "2009-start": ()}
CELLS = 50                       # 10 m dune cells per domain
GIS = range(1, 91)
# -----------------------------------------------------------------------------


# The seaward-most dune on a beach-start-to-crest profile (index into seg): the first peak, at least
# MIN_BEACH_CELLS from the beach start, with a dip of FOREDUNE_DROP_M or more behind it before the profile climbs
# above it again; the crest (the last point) where no peak qualifies, a single foredune
def first_dune(seg):
    for p in range(MIN_BEACH_CELLS, len(seg) - 1):
        if seg[p] < seg[p - 1] or seg[p] < seg[p + 1]:
            continue
        behind = seg[p + 1:]
        higher = np.flatnonzero(behind > seg[p])
        upto = behind[:higher[0]] if higher.size else behind
        if upto.size and seg[p] - upto.min() >= FOREDUNE_DROP_M:
            return p
    return len(seg) - 1


# Per profile: crest, beach start, toe and their cell indices on one domain's straightened grid
def crest_and_toe(z, i0, i1, X):
    n = z.shape[0]
    out = {k: np.full(n, np.nan) for k in ("crest", "toe", "beach_start_idx", "crest_idx", "toe_idx",
                                           "foredune_idx")}
    above = z > X.BEACH_START_THR_M
    sb = np.where(above.any(axis=1), above.argmax(axis=1), -1)
    for i in range(n):
        if sb[i] < 0:
            continue
        a = max(i0, int(sb[i])) if X.CLIP_WINDOW_TO_BEACH else i0
        b = min(i1, z.shape[1])
        if b <= a:
            continue
        win = z[i, a:b]
        if (win <= X.NODATA_SENTINEL_M + 1e-9).any():
            continue
        valid = win > X.SENTINEL_WATER_M + 1e-9
        if not valid.any():
            continue
        c = a + int(np.argmax(np.where(valid, win, -np.inf)))
        out["crest"][i], out["crest_idx"][i], out["beach_start_idx"][i] = z[i, c], c, sb[i]
        s = int(sb[i])
        if c - s < MIN_BEACH_CELLS:
            continue
        seg = z[i, s:c + 1]
        if (seg <= X.NODATA_SENTINEL_M + 1e-9).any():
            continue
        # The chord runs to the first dune, not the highest: on a field of several ridges the highest is a
        # landward one, and the deepest dip under a chord to it is a swale inside the field (2026-11-03)
        f = first_dune(seg)
        out["foredune_idx"][i] = s + f
        seg = seg[:f + 1]
        line = seg[0] + (seg[-1] - seg[0]) * np.arange(len(seg)) / (len(seg) - 1)
        k = int(np.argmax(line - seg))
        if k in (0, len(seg) - 1):
            k = 0 if seg[0] >= seg[-1] else k      # a straight or convex face: the toe is the beach start
        out["toe"][i], out["toe_idx"][i] = seg[k], s + k
    return out


# (survey, gis, profiles m MHW, dune window i0, i1) for the start's own arrays ("start"), then every survey
def domain_profiles(start, domains=GIS):
    import rasterio
    tv = S.tv
    prod, surveys = S.STARTS[start]
    windows = json.loads(tv.picks_dir(start).joinpath(
        f"HAT_dune_search_windows_{tv.topo_dirs(start)[2]}.json").read_text(encoding="utf-8"))
    frames = {g: S.frame(start, g, windows) for g in domains}
    clips = S.product(prod).gapfill_1m
    # Every domain's grid: a raster's unit is read off the whole reach
    grids = {g: S.clip_grid(clips / f"clip_domain_{g}_filled.tif") for g in GIS}
    offsets = S.SURVEY_VERTICAL_OFFSET_M.get(start, {})
    # The start's own arrays: the dunes the model starts from
    for g in domains:
        shear, c0, i0, i1 = frames[g]
        yield "start", g, S.profiles(np.load(tv.npy_dirs(start)[0] / f"domain_{g}.npy").astype(float), shear, c0), i0, i1
    print("start arrays done", flush=True)
    for year, (kind, path, flown, state, note) in surveys.items():
        files = sorted(path.parent.glob(path.name)) if "*" in path.name else [path]
        srcs = [rasterio.open(f) for f in files]
        try:
            unit = "m" if kind == "mosaic" else S.unit_of(srcs, start, grids)[0]
            scale = S.US_FT if unit == "ft" else 1.0
            for g in domains:
                shear, c0, i0, i1 = frames[g]
                on_clip = (S.mosaic_on_clip(srcs[0], grids[g]) if kind == "mosaic"
                           else S.raster_on_clip(srcs, grids[g]))
                yield year, g, S.profiles(S.downsample(on_clip * scale + offsets.get(year, 0.0)), shear, c0), i0, i1
        finally:
            for s in srcs:
                s.close()
        print(f"{year} done", flush=True)


# Every survey's (and the start arrays') toe, crest and beach width per cell, with the crest script's status
def survey_toes(start):
    X, tv = S.X, S.tv
    rows = []
    for survey, g, z, i0, i1 in domain_profiles(start):
        r = crest_and_toe(z, i0, i1, X)
        for cell in range(CELLS):
            rows.append(dict(survey=survey, gis=g, cell=cell,
                             crest_m_navd88=r["crest"][cell] + MHW_NAVD88_M,
                             toe_m_navd88=r["toe"][cell] + MHW_NAVD88_M,
                             beach_to_crest_m=(r["crest_idx"][cell] - r["beach_start_idx"][cell]) * 10,
                             toe_to_crest_m=(r["crest_idx"][cell] - r["toe_idx"][cell]) * 10))
    t = pd.DataFrame(rows)
    # The beach width the module uses (2026-11-03): the dry beach, beach start to toe, no part of the dune. Where
    # the crest stands within MIN_BEACH_CELLS of the beach start no toe is picked, and there is no dry beach: 0
    t["beach_width_m"] = np.where(t.toe_to_crest_m.notna(), t.beach_to_crest_m - t.toe_to_crest_m,
                                  np.where(t.beach_to_crest_m < MIN_BEACH_CELLS * 10, 0.0, np.nan))
    # The crest script's spike screen: its status beside each survey cell, and its crest as a check
    sc = pd.read_csv(tv.product_dir(start) / "survey-crests" / "survey_dune_crests.csv")[["survey", "gis", "cell", "status", "crest_m_navd88"]]
    sc["survey"] = sc.survey.astype(str)
    t["survey"] = t.survey.astype(str)
    t = t.merge(sc.rename(columns={"crest_m_navd88": "crest_script_m_navd88"}),
                on=["survey", "gis", "cell"], how="left")
    t.loc[t.survey == "start", "status"] = "ok"
    same = (t.crest_m_navd88 - t.crest_script_m_navd88).abs()
    print(f"crest check against survey_dune_crests.csv: max |difference| {same.max():.4f} m")
    return t


# One toe and beach width per cell (m NAVD88, m): survey medians, the domain median where a cell has none
def cell_table(t, start):
    st = t[(t.survey != "start") & (t.status == "ok")].copy()
    st["survey"] = st.survey.astype(int)
    st = st[~st.survey.isin(TOE_SURVEYS_EXCLUDED[start])]
    good = st[(st.toe_m_navd88 > TOE_MIN_M_NAVD88)
              & (st.toe_m_navd88 < st.crest_m_navd88 - TOE_BELOW_CREST_M)]
    columns = {"toe_m_navd88": good.groupby(["gis", "cell"]).toe_m_navd88.median(),
               "beach_width_m": st.groupby(["gis", "cell"]).beach_width_m.median(),
               "beach_to_crest_m": st.groupby(["gis", "cell"]).beach_to_crest_m.median()}
    full = pd.MultiIndex.from_product([GIS, range(CELLS)], names=["gis", "cell"])
    out = pd.DataFrame(index=full)
    for name, series in columns.items():
        own = series.reindex(full)
        out[name] = own.fillna(series.groupby(level="gis").median().reindex(full.get_level_values("gis")).set_axis(full))
        out[name.replace("_m_navd88", "").replace("_m", "") + "_source"] = np.where(
            own.notna(), "cell", np.where(out[name].notna(), "domain median", "none"))
    out["surveys_with_toe"] = good.groupby(["gis", "cell"]).size().reindex(full).fillna(0).astype(int)
    return out.reset_index()


# Figures

# The two starts in the house window colours (1996 red, 2009 amber, as HAT_stable_dune_figures.py)
def _start_colours():
    from site_layer import hat_figure_style as fs
    return {"1996-start": fs.C_1984, "2009-start": fs.C["ADDED"]}


# Both starts' tables, where built: {start: (survey rows, cell table)}
def _tables():
    out = {}
    for start in sorted(S.STARTS, key=lambda s: int(s[:4])):     # calibration start first
        d = S.tv.product_dir(start) / "dune-toes"
        if (d / "dune_toe_cells.csv").exists():
            out[start] = (pd.read_csv(d / "survey_dune_toes.csv"), pd.read_csv(d / "dune_toe_cells.csv"))
    if not out:
        raise SystemExit("no dune_toe_cells.csv yet: run with --start first")
    return out


# Per GIS: the median and quartiles of a cell column across the domain's 50 cells
def _by_domain(values):
    g = values.groupby(level="gis")
    return g.median(), g.quantile(0.25), g.quantile(0.75)


def figures(beach_width_limit_m=25.0):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    from matplotlib.patches import Patch
    from site_layer import hat_figure_style as fs
    fs.apply_style()
    colours = _start_colours()
    tables = _tables()
    out_dir = fs.figure_dir("inputs", "1-domains", "dune_toes")
    span = " and ".join(f"{s} ({', '.join(str(y) for y in S.STARTS[s][1] if y not in TOE_SURVEYS_EXCLUDED[s])})"
                        for s in tables)

    # 1. Alongshore: toe, dune height above the toe, beach width
    fig, axes = plt.subplots(3, 1, sharex=True, figsize=fs.figsize("double", height=7.0),
                             layout="constrained")
    for start, (t, cells) in tables.items():
        c = colours[start]
        cells = cells.set_index(["gis", "cell"])
        st = t[(t.survey != "start") & (t.status == "ok")].copy()
        st = st[~st.survey.astype(int).isin(TOE_SURVEYS_EXCLUDED[start])]
        crest = st.groupby(["gis", "cell"]).crest_m_navd88.median().reindex(cells.index)
        series = (cells.toe_m_navd88, crest - cells.toe_m_navd88, cells.beach_width_m)
        for ax, values in zip(axes, series):
            mid, lo, hi = _by_domain(values)
            ax.fill_between(mid.index, lo, hi, color=c, alpha=0.18, lw=0, step="mid", zorder=2)
            ax.plot(mid.index, mid, color=c, lw=1.6, drawstyle="steps-mid", zorder=3)
        # Domains whose toe is mostly the domain median, and the ones with none
        share = (cells.toe_source != "cell").groupby(level="gis").mean()
        filled = share[(share >= 0.5) & (cells.toe_source == "none").groupby(level="gis").mean().lt(1)].index
        empty = (cells.toe_source == "none").groupby(level="gis").all()
        empty = empty[empty].index
        mid = _by_domain(cells.toe_m_navd88)[0]
        axes[0].plot(filled, mid.reindex(filled), ls="none", marker="o", ms=5, mfc="white", mec=c, mew=1.2, zorder=4)
        for g in empty:
            axes[0].plot(g, 1.15, ls="none", marker="x", ms=6, mew=1.4, color=c, zorder=4)
            axes[0].annotate(f"GIS {g}: no toe", (g, 1.15), xytext=(6, 0), textcoords="offset points",
                             ha="left", va="center", fontsize=7, color=fs.INK_MUTED)
    axes[2].axhline(beach_width_limit_m, color=fs.INK_MUTED, lw=0.8, ls="--", zorder=2)
    labels = ("Dune toe (m NAVD88)", "Crest above toe (m)", "Beach width (m)")
    for k, (ax, label) in enumerate(zip(axes, labels)):
        ax.set_xlim(0.5, 90.5)
        ax.set_ylabel(label)
        ax.grid(axis="y")
        fs.open_frame(ax)
        fs.town_bands(ax, label=(k == 0))
        fs.structures(ax, label=(k == len(axes) - 1))
    axes[0].set_ylim(1.0, None)
    axes[2].set_ylim(0, None)
    axes[-1].set_xlabel(fs.DOMAIN_AXIS_LABEL)
    handles = [Line2D([], [], color=colours[s], lw=1.6, label=f"{s}: domain median") for s in tables]
    handles += [Patch(color=fs.INK_MUTED, alpha=0.25, label="middle half of the domain's 50 cells"),
                Line2D([], [], ls="none", marker="o", ms=5, mfc="white", mec=fs.INK, label="half or more of the cells take the domain median"),
                Line2D([], [], color=fs.INK_MUTED, lw=0.8, ls="--",
                       label=f"beach width {beach_width_limit_m:g} m: no collision loss at or beyond")]
    fig.legend(handles=handles, loc="outside lower center", ncol=2, fontsize=7, frameon=False)
    png = out_dir / "dune_toes_alongshore.png"
    fs.save(fig, png, close=True)
    fs.record_caption(png, (
        "The dune-collision inputs per GIS domain, as the runner reads them "
        "(input_prep/1-barrier3d-domains/3-dune-toes): top, each 10 m dune cell's toe, the median over the "
        f"surveys whose toe is above {TOE_MIN_M_NAVD88:g} m NAVD88 and at least {TOE_BELOW_CREST_M:g} m below "
        "that survey's crest; middle, the cell's median crest over the same surveys less its toe, the height "
        "a storm can take off; bottom, the cell's median beach width (the dry beach, beach start to toe), with the width "
        "beyond which the module takes nothing. Lines: the domain median; bands: the middle half of the "
        "domain's cells. A cell without a value takes its domain's median (open circles where half or more "
        "do); a domain with none is left out of the collision (cross). Surveys: " + span + "; spike-screened "
        "cells left out."))
    print(f"wrote {png}")

    # 2. How much a fixed toe hides: each survey's domain-median toe against the value used
    fig, axes = plt.subplots(len(tables), 1, sharex=True, sharey=True, layout="constrained",
                             figsize=fs.figsize("double", height=2.6 * len(tables)), squeeze=False)
    for k, (ax, (start, (t, cells))) in enumerate(zip(axes[:, 0], tables.items())):
        c = colours[start]
        st = t[(t.survey != "start") & (t.status == "ok")].copy()
        st["survey"] = st.survey.astype(int)
        st = st[~st.survey.isin(TOE_SURVEYS_EXCLUDED[start])]
        good = st[(st.toe_m_navd88 > TOE_MIN_M_NAVD88) & (st.toe_m_navd88 < st.crest_m_navd88 - TOE_BELOW_CREST_M)]
        per_survey = good.groupby(["survey", "gis"]).toe_m_navd88.median().reset_index()
        ax.scatter(per_survey.gis, per_survey.toe_m_navd88, s=9, color=c, alpha=0.45, lw=0, zorder=3)
        used = cells.groupby("gis").toe_m_navd88.median()
        ax.plot(used.index, used, ls="none", marker="_", ms=7, mew=1.6, color=fs.INK, zorder=4)
        spread = per_survey.groupby("gis").toe_m_navd88.agg(lambda x: x.max() - x.min())
        ax.set_title(f"{start}: {per_survey.survey.nunique()} surveys", loc="left")
        ax.set_title(f"median survey-to-survey range {spread.median():.2f} m, widest {spread.max():.2f} m "
                     f"(GIS {spread.idxmax()})", loc="right", fontsize=7, color=fs.INK_MUTED)
        ax.set_xlim(0.5, 90.5)
        ax.set_ylim(1.0, 4.9)          # headroom so no gridline runs through the town names
        ax.set_ylabel("Dune toe (m NAVD88)")
        ax.grid(axis="y")
        fs.open_frame(ax)
        fs.town_bands(ax, label=(k == 0))
        fs.structures(ax, label=(k == len(tables) - 1))
    axes[-1, 0].set_xlabel(fs.DOMAIN_AXIS_LABEL)
    fig.legend(handles=[Line2D([], [], ls="none", marker="o", ms=4, color=fs.INK_MUTED, alpha=0.6,
                               label="one survey's domain-median toe"),
                        Line2D([], [], ls="none", marker="_", ms=7, mew=1.6, color=fs.INK,
                               label="the fixed toe the runner uses (domain median of the cell values)")],
               loc="outside lower center", ncol=2, fontsize=7, frameon=False)
    png = out_dir / "dune_toes_survey_spread.png"
    fs.save(fig, png, close=True)
    fs.record_caption(png, (
        "How much the fixed dune toe hides. Per GIS domain, each survey's median toe over the domain's "
        "screened cells (dots), against the one value the runner holds for the whole run (tick: the domain "
        "median of the per-cell survey medians). The spread between surveys is what the fixed toe leaves "
        f"out, a choice kept for simplicity (2026-10-29). Surveys: {span}."))
    print(f"wrote {png}")


# The method figure's example: a typical cell (toe 2.93 m NAVD88, beach 35 m, crest 3.25 m above the toe,
# a toe in all 8 surveys), and the storm series of its window
EXAMPLE = ("1996-start", 22, 32)
HEADROOM_M = 2.3                 # above the crest on the profile panels, room for the notes
EXAMPLE_STORMS = ("data", "hatteras_init", "3-env-forcings", "3-storms", "hindcast_storms", "1997_2009",
                  "1997_2009_storms_v3_240_summary.csv")


def method_figure(beach_width_limit_m=25.0):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    from hatteras_paths import WORKSPACE
    from site_layer import hat_figure_style as fs
    fs.apply_style()
    colours = _start_colours()
    tables = _tables()
    start, g, cell = EXAMPLE
    X = S.X
    c = colours[start]
    excluded = TOE_SURVEYS_EXCLUDED[start]
    profiles = {survey: (z[cell], crest_and_toe(z, i0, i1, X))
                for survey, _, z, i0, i1 in domain_profiles(start, [g])}
    fixed = tables[start][1].set_index(["gis", "cell"]).loc[(g, cell)]
    storms = pd.read_csv(WORKSPACE.joinpath(*EXAMPLE_STORMS))
    peak = storms.Rhigh * 10 + MHW_NAVD88_M                          # dam MHW to m NAVD88

    fig = plt.figure(figsize=fs.figsize("double", height=6.4), layout="constrained")
    grid = fig.add_gridspec(2, 2)
    ax_a, ax_b, ax_d = fig.add_subplot(grid[0, 0]), fig.add_subplot(grid[0, 1]), fig.add_subplot(grid[1, 1])
    ax_c_top, ax_c = grid[1, 0].subgridspec(2, 1, height_ratios=(1, 2), hspace=0.05).subplots(sharex=True)

    # One survey's profile, the anatomy drawn on it: the survey whose toe is nearest the cell's fixed one
    survey = min((s for s in profiles if s != "start" and int(s) not in excluded
                  and np.isfinite(profiles[s][1]["toe"][cell])),
                 key=lambda s: abs(profiles[s][1]["toe"][cell] + MHW_NAVD88_M - fixed.toe_m_navd88))
    z, r = profiles[survey]
    sb, ci, ti = (int(r[k][cell]) for k in ("beach_start_idx", "crest_idx", "toe_idx"))
    lo, hi = max(sb - 4, 0), min(ci + 6, len(z))
    x = (np.arange(lo, hi) - sb) * 10.0
    zz = z[lo:hi] + MHW_NAVD88_M
    toe, crest, base = zz[ti - lo], zz[ci - lo], zz[sb - lo]
    for ax in (ax_a, ax_b):
        ax.plot(x, zz, color=fs.INK, lw=1.6, zorder=4)
        ax.set_xlim(x[0], x[-1])
        ax.set_ylim(0, crest + HEADROOM_M)
        ax.set_xlabel("Distance landward of the beach start (m)")
        ax.set_ylabel("Elevation (m NAVD88)")
        fs.open_frame(ax)
    ax_a.fill_between(x, zz, 0, color=fs.C["BASE_FILL"], lw=0, zorder=1)
    ax_a.plot([0, (ci - sb) * 10], [base, crest], color=fs.INK_MUTED, lw=0.9, ls="--", zorder=3)
    for (xi, zi, name, dx, dy, ha) in ((0, base, "beach start", 8, -12, "left"),
                                        ((ti - sb) * 10, toe, "toe", 6, 6, "left"),
                                        ((ci - sb) * 10, crest, "crest", -6, 6, "right")):
        ax_a.plot(xi, zi, "o", ms=6, color=c if name == "toe" else fs.INK, zorder=5)
        ax_a.annotate(f"{name} {zi:.2f} m", (xi, zi), xytext=(dx, dy), textcoords="offset points", ha=ha,
                      fontsize=7.5, color=fs.INK)
    # The dry beach: beach start to toe, no part of the dune face
    ax_a.annotate("", (0, 0.35), ((ti - sb) * 10, 0.35), arrowprops=dict(arrowstyle="<->", color=fs.INK_MUTED, lw=0.8))
    ax_a.text((ti - sb) * 5, 0.5, f"beach width {(ti - sb) * 10} m", ha="center", fontsize=7, color=fs.INK_MUTED)
    # The height collision can remove, measured off on the open water seaward of the beach
    xc = x[0] + 12
    for level, xi in ((toe, (ti - sb) * 10), (crest, (ci - sb) * 10)):
        ax_a.plot([xc - 4, xi], [level, level], color=c, lw=0.6, ls=":", zorder=2)
    ax_a.annotate("", (xc, toe), (xc, crest), arrowprops=dict(arrowstyle="<->", color=c, lw=0.9))
    ax_a.text(xc + 3, (toe + crest) / 2, f"{crest - toe:.2f} m that\ncollision can\nremove", fontsize=7,
              color=fs.INK, va="center")
    ax_a.text(0.98, 0.97, f"GIS {g}, cell {cell}, {survey} survey\ntoe: furthest below the dashed chord\n"
              f"the runner keeps: toe {fixed.toe_m_navd88:.2f} m, beach {fixed.beach_width_m:g} m (survey medians)",
              transform=ax_a.transAxes, fontsize=7, color=fs.INK_MUTED, va="top", ha="right")
    fs._title(ax_a, 0, "Reading one profile")

    # The regimes the toe and crest divide the storm water into (Sallenger 2000), and the window's storms
    bands = ((0, toe, fs.C["WATER"], "swash: below the toe, nothing"),
             (toe, crest, fs.C["ACCENT_FILL"], "collision: lowered, never\nbelow the toe (this module)"),
             (crest, crest + HEADROOM_M, fs.C["ADDED_FILL"], "overwash: over the crest\n(Barrier3D)"))
    for b0, b1, colour, label in bands:
        ax_b.axhspan(b0, b1, color=colour, alpha=0.55, lw=0, zorder=0)
        n = int(((peak >= b0) & (peak < b1)).sum())
        ax_b.text(x[0] + 4, (b0 + b1) / 2, f"{label}\n{n} storm{'s' if n != 1 else ''}", fontsize=7,
                  va="center", color=fs.INK, zorder=6,
                  bbox=dict(facecolor="white", alpha=0.8, edgecolor="none", boxstyle="square,pad=0.2"))
    jitter = np.random.default_rng(0).uniform(-4, 4, len(peak))
    ax_b.scatter(x[-1] - 10 + jitter, peak, s=8, color=fs.INK, alpha=0.6, lw=0, zorder=5, clip_on=False)
    ax_b.text(x[-1] - 2, crest + HEADROOM_M - 0.15, "storm peaks", ha="right", fontsize=7, color=fs.INK_MUTED, va="top")
    fs._title(ax_b, 1, "What the toe and crest decide")

    # The beach factor, and how many cells it lets collision reach
    w = np.linspace(0, 160, 200)
    ax_c_top.plot(w, np.clip(1 - w / beach_width_limit_m, 0, None), color=fs.INK, lw=1.6)
    ax_c_top.set_ylabel("factor")
    ax_c_top.set_ylim(0, 1.1)
    ax_c_top.set_yticks([0, 1])
    ax_c_top.text(beach_width_limit_m + 3, 0.55, f"1 − W / {beach_width_limit_m:g} m", fontsize=7.5,
                  color=fs.INK)
    shares = []
    for s, (_, cells) in tables.items():
        widths = cells.beach_width_m.dropna()
        ax_c.hist(widths, bins=np.arange(0, 165, 10), histtype="step", lw=1.6, color=colours[s], label=s)
        shares.append(f"{s[:4]} {100 * (widths < beach_width_limit_m).mean():.0f}%")
    for ax in (ax_c_top, ax_c):
        ax.axvline(beach_width_limit_m, color=fs.INK_MUTED, lw=0.8, ls="--")
        ax.set_xlim(0, 160)
        fs.open_frame(ax)
    ax_c.set_xlabel("Beach width W (m)")
    ax_c.set_ylabel("Dune cells")
    ax_c.text(0.98, 0.95, "cells under " + f"{beach_width_limit_m:g} m:\n" + ", ".join(shares),
              transform=ax_c.transAxes, ha="right", va="top", fontsize=7, color=fs.INK)
    ax_c.legend(loc="center right", fontsize=7, frameon=False)
    fs._title(ax_c_top, 2, "A wide beach shields the dune")

    # The same cell in every survey: each one's toe, and the one value the runner keeps
    for s, (zs, rs) in profiles.items():
        if s == "start" or int(s) in excluded or not np.isfinite(rs["toe"][cell]):
            continue
        ax_d.plot((np.arange(lo, hi) - sb) * 10.0, zs[lo:hi] + MHW_NAVD88_M, color=c, lw=0.9, alpha=0.45)
        ax_d.plot((rs["toe_idx"][cell] - sb) * 10, rs["toe"][cell] + MHW_NAVD88_M, "o", ms=4.5, color=c, zorder=4)
    ax_d.axhline(fixed.toe_m_navd88, color=fs.INK, lw=1.2, ls="--", zorder=3)
    ax_d.text(x[0] + 3, fixed.toe_m_navd88 + 0.12, f"fixed toe {fixed.toe_m_navd88:.2f} m (median)",
              fontsize=7.5, color=fs.INK)
    ax_d.set_xlim(x[0], x[-1])
    ax_d.set_ylim(0, crest + HEADROOM_M)
    ax_d.set_xlabel("Distance landward of the beach start (m)")
    ax_d.set_ylabel("Elevation (m NAVD88)")
    fs.open_frame(ax_d)
    used = [s for s in profiles if s != "start" and int(s) not in excluded]
    fs._title(ax_d, 3, f"One value from {len(used)} surveys")
    ax_d.legend(handles=[Line2D([], [], color=c, lw=0.9, alpha=0.6, label="one survey's profile"),
                         Line2D([], [], ls="none", marker="o", ms=4.5, color=c, label="its toe")],
                loc="upper left", fontsize=7, frameon=False)

    png = fs.figure_dir("inputs", "1-domains", "dune_toes") / "dune_toes_method.png"
    fs.save(fig, png, close=True)
    fs.record_caption(png, (
        f"How the dune-collision inputs are read, on one typical cell (GIS {g}, cell {cell}, {start}). "
        "(a) A surveyed cross-shore profile: the beach start is the first point above the beach threshold, "
        "the crest the highest point in the dune window, and the toe the point furthest below the straight "
        "chord between them (dashed); the beach width is the dry beach, from the beach start to the toe, and the crest "
        "above the toe is what collision can remove. (b) The toe and crest split a storm's peak water level "
        "(R_high, Stockdon run-up on the tide gauge) into the regimes of Sallenger (2000): below the toe the "
        "dune is untouched, between toe and crest the module lowers it storm by storm (Larson et al. 2004) but "
        "never below the toe, and over the crest Barrier3D's overwash takes over; dots are the peak water "
        f"levels of the {len(peak)} storms in the 1997-2009 series, counted per band at this cell. (c) Collision "
        f"is scaled by max(0, 1 − W / {beach_width_limit_m:g} m), so a beach {beach_width_limit_m:g} m or wider "
        "shields its dune completely; below, the beach widths of every dune cell for both starts. (d) The same "
        "cell in every survey that found a toe there, with each survey's toe; the runner holds one value per "
        "cell for the whole run, the median (dashed)."))
    print(f"wrote {png}")


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[1],
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--start", choices=list(S.STARTS))
    ap.add_argument("--figures", action="store_true",
                    help="draw the figures from the saved tables of every start built so far")
    args = ap.parse_args()
    if args.figures:
        figures()
        method_figure()
        return
    if not args.start:
        ap.error("--start is required unless --figures")
    t = survey_toes(args.start)
    cells = cell_table(t, args.start)
    out = S.tv.product_dir(args.start) / "dune-toes"
    out.mkdir(parents=True, exist_ok=True)
    t.round(3).to_csv(out / "survey_dune_toes.csv", index=False)
    cells.round(3).to_csv(out / "dune_toe_cells.csv", index=False)
    print(cells.groupby("toe_source").size().to_string())
    print(cells[["toe_m_navd88", "beach_width_m"]].describe().round(2).to_string())
    print(f"wrote {out / 'survey_dune_toes.csv'}\nwrote {out / 'dune_toe_cells.csv'}")


if __name__ == "__main__":
    main()
