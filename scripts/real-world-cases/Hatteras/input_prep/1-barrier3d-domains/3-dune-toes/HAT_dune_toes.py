"""
Dune toes and beach widths per 10 m dune cell, from a start's lidar surveys: the dune-collision inputs.

    python scripts/real-world-cases/Hatteras/input_prep/1-barrier3d-domains/3-dune-toes/HAT_dune_toes.py --start 1996-start
    python scripts/real-world-cases/Hatteras/input_prep/1-barrier3d-domains/3-dune-toes/HAT_dune_toes.py --start 2009-start
    python scripts/real-world-cases/Hatteras/input_prep/1-barrier3d-domains/3-dune-toes/HAT_dune_toes.py --start 2009-start --cells
    python scripts/real-world-cases/Hatteras/input_prep/1-barrier3d-domains/3-dune-toes/HAT_dune_toes.py --figures

The dune-collision module (cascade/dune_collision.py) lowers a dune the water
reaches but does not overtop, down to its toe, and less on a wide beach. It needs
each cell's toe and beach width. Step 3 of the stage, after 2-survey-crests:
both are taken from the same profiles HAT_survey_dune_crests.py reads crests from (the start's own arrays, then every
survey of its window, on the start's frame and dune windows): the crest is the
highest point of the dune window, the beach start the first cell above the beach
threshold, the toe the point of the profile furthest below the straight line from the beach start to the
first dune (the beach start on a straight or convex face), and
the beach width the dry beach, the beach start to the toe (0 where the crest stands within one cell of the
beach start). Until 2026-11-03 the beach width ran to the crest, the dune face included (kept as
beach_to_crest_m): it made a broad dune face read as a wide beach, and the toe is where the swash meets the dune.
The first dune is the seaward-most peak with a dip of FOREDUNE_DROP_M behind it, the crest where none has one.
Until 2026-11-03 the line ran to the crest: on a field of several ridges (GIS ~32-35, 63-67 in 1997) the toe then
fell in a swale 35-90 m inside the field, against the dune lines traced from imagery
(dune-collision/1-observations/2-toes-vs-dunelines/); the earlier tables are kept in
9-dune-toes/0-previous-methods/<start>/deepest-dip-until-2026-11-03/.

Each domain gets one toe, fixed for the run (since 2026-11-04): the mean of its screened cells over the lidar
surveys nearest the start (TOE_DOMAIN_SURVEYS: 1996 and 1997; 2004, 2005 and 2009), the toe screen as below. Until
then each cell had its own (2026-10-29), kept as toe_cell_median_m_navd88:
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
read off one profile, made one per domain and used by the module, on one typical cell, EXAMPLE).

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-11-04
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
# The toe the runner uses since 2026-11-04: one per domain, the mean of its screened cells over the lidar surveys
# nearest the start, (start whose survey table holds them, surveys). No lidar before 1996, so 1996 and 1997; for 2009
# the surveys up to it, read on the 1996 start's windows (nothing after the start). Per-cell survey medians before
# (kept as toe_cell_median_m_navd88): in the 2009-2025 test they took the toes from 2014-2019, after the start, and a
# single cell's low toe over-erodes (the erosion goes with the square of the water above it)
BEACH_WIDTH_LIMIT_M = 25.0       # the module's beach_width_limit (dry beach), for the domain beach factor
TOE_DOMAIN_SURVEYS = {"1996-start": ("1996-start", (1996, 1997)), "2009-start": ("1996-start", (2004, 2005, 2009))}
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
    # One window for every profile, or one per profile (the survey windows that follow the dune, since 2026-11-03)
    i0s, i1s = np.broadcast_to(np.asarray(i0, int), (n,)), np.broadcast_to(np.asarray(i1, int), (n,))
    above = z > X.BEACH_START_THR_M
    sb = np.where(above.any(axis=1), above.argmax(axis=1), -1)
    for i in range(n):
        if sb[i] < 0:
            continue
        i0, i1 = max(int(i0s[i]), 0), int(i1s[i])
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
                # The survey's windows follow the dune, as the crest script moved them (S.survey_windows, 2026-11-03)
                w0, w1 = S.survey_windows(start, year, g, i0, i1)
                on_clip = (S.mosaic_on_clip(srcs[0], grids[g]) if kind == "mosaic"
                           else S.raster_on_clip(srcs, grids[g]))
                yield year, g, S.profiles(S.downsample(on_clip * scale + offsets.get(year, 0.0)), shear, c0), w0, w1
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
    # The toe the runner reads: the domain mean over the surveys nearest the start (TOE_DOMAIN_SURVEYS)
    src, years = TOE_DOMAIN_SURVEYS[start]
    n = t if src == start else pd.read_csv(S.tv.product_dir(src) / "dune-toes" / "survey_dune_toes.csv")
    n = n[(n.survey.astype(str).isin([str(y) for y in years])) & (n.status == "ok")]
    n = n[(n.toe_m_navd88 > TOE_MIN_M_NAVD88) & (n.toe_m_navd88 < n.crest_m_navd88 - TOE_BELOW_CREST_M)]
    dm = n.groupby("gis").toe_m_navd88.mean()
    out = out.rename(columns={"toe_m_navd88": "toe_cell_median_m_navd88"})
    out["toe_m_navd88"] = dm.reindex(full.get_level_values("gis")).to_numpy()
    # The beach the same way since 2026-11-04, one value per domain: the domain mean of the module's beach factor,
    # max(0, 1 - W / BEACH_WIDTH_LIMIT_M), over every ok cell of those surveys (beach_factor). The mean of the
    # factor, not of the width: the factor is 0 beyond the limit, so a mean width over it switches the module off on
    # a domain whose dune is half behind a narrow beach. beach_width_m is the width that gives that factor,
    # BEACH_WIDTH_LIMIT_M (1 - factor), what the module reads; the per-cell medians kept as beach_width_cell_median_m
    nw = t if src == start else pd.read_csv(S.tv.product_dir(src) / "dune-toes" / "survey_dune_toes.csv")
    nw = nw[(nw.survey.astype(str).isin([str(y) for y in years])) & (nw.status == "ok")]
    bf = np.clip(1 - nw.beach_width_m / BEACH_WIDTH_LIMIT_M, 0, None).groupby(nw.gis).mean()
    out = out.rename(columns={"beach_width_m": "beach_width_cell_median_m"})
    out["beach_factor"] = bf.reindex(full.get_level_values("gis")).to_numpy()
    out["beach_width_m"] = BEACH_WIDTH_LIMIT_M * (1 - out.beach_factor)
    out["beach_width_source"] = np.where(out.beach_width_m.notna(), out.toe_source, "none")
    print(f"beach factor: domain mean over {src} {years}: island {bf.mean():.2f}, range {bf.min():.2f}-{bf.max():.2f}")
    out["toe_source"] = np.where(out.toe_m_navd88.notna(), f"domain mean {src} {'+'.join(map(str, years))}", "none")
    print(f"toe: domain mean over {src} {years}: island {dm.mean():.2f} m NAVD88, range {dm.min():.2f}-{dm.max():.2f}")
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


# The method figure's example: a typical cell (GIS 22, cell 32) and its domain, and the storm series of the
# calibration window
EXAMPLE = ("1996-start", 22, 32)
HEADROOM_M = 2.3                 # above the crest on the profile panels, room for the notes
EXAMPLE_STORMS = ("data", "hatteras_init", "3-env-forcings", "3-storms", "hindcast_storms", "1997_2009",
                  "1997_2009_storms_v3_240_summary.csv")
PERIOD_LABEL = {"1996-start": "1996–2009", "2009-start": "2009–2025"}


def method_figure():
    """How the toe is measured and used (redrawn 2026-11-04 for one toe per domain, no beach factor): (a) one
    profile, (b) the regimes the domain toe and the crest decide, (c) one domain, every cell's toe in its start
    surveys and the one domain toe, (d) the domain toes along the island for both periods."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from hatteras_paths import WORKSPACE
    from site_layer import hat_figure_style as fs
    fs.apply_style()
    colours = _start_colours()
    tables = _tables()
    start, g, cell = EXAMPLE
    X = S.X
    c = colours[start]
    years = [str(y) for y in TOE_DOMAIN_SURVEYS[start][1]]
    profiles = {survey: (z[cell], crest_and_toe(z, i0, i1, X))
                for survey, _, z, i0, i1 in domain_profiles(start, [g])}
    dom_toe = float(tables[start][1].groupby("gis").toe_m_navd88.first()[g])
    storms = pd.read_csv(WORKSPACE.joinpath(*EXAMPLE_STORMS))
    peak = storms.Rhigh * 10 + MHW_NAVD88_M                          # dam MHW to m NAVD88

    fig, ((ax_a, ax_b), (ax_c, ax_d)) = plt.subplots(2, 2, figsize=fs.figsize("double", height=6.4),
                                                     layout="constrained")

    # (a) One start survey's profile at the example cell: the survey whose toe is nearest the domain toe
    cands = [s for s in profiles if str(s) in years and np.isfinite(profiles[s][1]["toe"][cell])]
    survey = min(cands, key=lambda s: abs(profiles[s][1]["toe"][cell] + MHW_NAVD88_M - dom_toe))
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
    fi = int(r["foredune_idx"][cell])
    top = z[fi] + MHW_NAVD88_M
    ax_a.plot([0, (fi - sb) * 10], [base, top], color=fs.INK_MUTED, lw=0.9, ls="--", zorder=3)
    for (xi, zi, name, dx, dy, ha) in ((0, base, "beach start", 8, -12, "left"),
                                        ((ti - sb) * 10, toe, "toe", 6, -12, "left"),
                                        ((ci - sb) * 10, crest, "crest", -6, 6, "right")):
        ax_a.plot(xi, zi, "o", ms=6, color=c if name == "toe" else fs.INK, zorder=5)
        ax_a.annotate(f"{name} {zi:.2f} m", (xi, zi), xytext=(dx, dy), textcoords="offset points", ha=ha,
                      fontsize=7.5, color=fs.INK)
    # The largest gap under the line, the toe's definition
    k = (ti - sb) * 10
    chord_at = base + (top - base) * k / max((fi - sb) * 10, 1)
    if chord_at - toe > 0.15:
        ax_a.annotate("", (k, toe), (k, chord_at), arrowprops=dict(arrowstyle="<->", color=c, lw=0.8))
    ax_a.text(0.98, 0.97, f"GIS {g}, cell {cell}, {survey} survey\n"
              "toe: the point furthest below the dashed line\nfrom the beach start to the first dune",
              transform=ax_a.transAxes, fontsize=7, color=fs.INK_MUTED, va="top", ha="right")
    fs._title(ax_a, 0, "Reading the toe off one profile")

    # (b) The regimes the domain toe and this crest split the storms into (Sallenger 2000)
    bands = ((0, dom_toe, fs.C["WATER"], "swash: below the toe, nothing"),
             (dom_toe, crest, fs.C["ACCENT_FILL"], "collision: Larson lowers the\ndune, never below the toe"),
             (crest, crest + HEADROOM_M, fs.C["ADDED_FILL"], "overwash: over the crest\n(Barrier3D)"))
    for b0, b1, colour, label in bands:
        ax_b.axhspan(b0, b1, color=colour, alpha=0.55, lw=0, zorder=0)
        n = int(((peak >= b0) & (peak < b1)).sum())
        ax_b.text(x[0] + 4, (b0 + b1) / 2, f"{label}\n{n} storm{'s' if n != 1 else ''}", fontsize=7,
                  va="center", color=fs.INK, zorder=6,
                  bbox=dict(facecolor="white", alpha=0.8, edgecolor="none", boxstyle="square,pad=0.2"))
    ax_b.axhline(dom_toe, color=c, lw=1.2, ls="--", zorder=3)
    ax_b.text(x[-1] - 22, dom_toe - 0.12, f"domain toe {dom_toe:.2f} m", ha="right", va="top", fontsize=7,
              color=c, zorder=6)
    jitter = np.random.default_rng(0).uniform(-4, 4, len(peak))
    ax_b.scatter(x[-1] - 10 + jitter, peak, s=8, color=fs.INK, alpha=0.6, lw=0, zorder=5, clip_on=False)
    ax_b.text(x[-1] - 2, crest + HEADROOM_M - 0.15, "storm peaks, 1997–2009", ha="right", fontsize=7,
              color=fs.INK_MUTED, va="top")
    fs._title(ax_b, 1, "What the toe decides")

    # (c) One domain: every cell's toe in the start surveys, screened, and the one domain toe (their mean)
    t = tables[start][0]
    d = t[(t.gis == g) & t.survey.astype(str).isin(years) & (t.status == "ok") & t.toe_m_navd88.notna()].copy()
    ok = (d.toe_m_navd88 > TOE_MIN_M_NAVD88) & (d.toe_m_navd88 < d.crest_m_navd88 - TOE_BELOW_CREST_M)
    crest_med = d.groupby("cell").crest_m_navd88.median()
    ax_c.plot(crest_med.index * 10 + 5, crest_med, color=fs.INK_MUTED, lw=1.0, drawstyle="steps-mid",
              label="crest (survey median)")
    fills = [c, "white", fs.INK]
    for y, mfc in zip(years, fills):
        e = d[ok & (d.survey.astype(str) == y)]
        ax_c.plot(e.cell * 10 + 5, e.toe_m_navd88, ls="none", marker="o", ms=4.5, mew=1.0, mfc=mfc, mec=c,
                  label=f"toe, {y} survey")
    bad = d[~ok]
    if len(bad):
        ax_c.plot(bad.cell * 10 + 5, bad.toe_m_navd88, ls="none", marker="x", ms=4.5, color=fs.INK_MUTED,
                  label="screened out")
    ax_c.axhline(dom_toe, color=c, lw=1.6, ls="--", label=f"domain toe {dom_toe:.2f} m (mean)")
    ax_c.axvline(cell * 10 + 5, color=fs.INK_MUTED, lw=0.6, ls=":")
    ax_c.text(cell * 10 + 9, 0.6, f"cell {cell}, panel (a)", fontsize=6.5, color=fs.INK_MUTED)
    ax_c.set_xlim(0, CELLS * 10)
    ax_c.set_ylim(0.5, 9.0)
    ax_c.set_xlabel(f"Alongshore position in GIS {g}, south to north (m)")
    ax_c.set_ylabel("Elevation (m NAVD88)")
    ax_c.legend(loc="upper left", ncol=2, fontsize=6.5, frameon=False)
    fs.open_frame(ax_c)
    fs._title(ax_c, 2, f"One toe for GIS {g}")

    # (d) The domain toes along the island, both periods
    ax_d.axvspan(0.5, 5.5, color=fs.INK_MUTED, alpha=0.12, lw=0)
    ax_d.text(3, 1.08, "module\noff", ha="center", va="bottom", fontsize=6.5, color=fs.INK_MUTED)
    for s, (_, cells) in tables.items():
        v = cells.groupby("gis").toe_m_navd88.first()
        src_years = ", ".join(str(y) for y in TOE_DOMAIN_SURVEYS[s][1])
        ax_d.plot(v.index, v, color=colours[s], lw=1.6, drawstyle="steps-mid",
                  label=f"{PERIOD_LABEL.get(s, s)} run ({src_years} surveys)")
    ax_d.plot(g, dom_toe, "o", ms=6, mfc="white", mec=c, mew=1.4, zorder=5)
    ax_d.annotate(f"GIS {g}", (g, dom_toe), xytext=(-6, -16), textcoords="offset points", ha="right", fontsize=7, color=fs.INK)
    ax_d.set_xlim(0.5, 90.5)
    ax_d.set_ylim(1.0, 4.5)
    ax_d.set_xlabel(fs.DOMAIN_AXIS_LABEL)
    ax_d.set_ylabel("Domain toe (m NAVD88)")
    ax_d.grid(axis="y")
    ax_d.legend(loc="upper right", fontsize=7, frameon=False)
    fs.open_frame(ax_d)
    fs._title(ax_d, 3, "One toe per domain, along the island")

    png = fs.figure_dir("inputs", "1-domains", "dune_toes") / "dune_toes_method.png"
    fs.save(fig, png, close=True)
    d.assign(kept=ok).round(3).to_csv(png.parent / "supporting" / "dune_toes_method_domain.csv", index=False)
    fs.record_caption(png, (
        "How the dune toe is measured and used by the dune-collision module. (a) One lidar profile (GIS "
        f"{g}, cell {cell}, {survey}): the beach start is the first point above the beach threshold, the first "
        "dune the most seaward peak with a 1 m dip behind it (the crest on a single foredune), and the toe the "
        "point furthest below the straight line between them (dashed). (b) The domain toe and the crest split "
        "each storm's peak water level into the regimes of Sallenger (2000): below the toe nothing happens; "
        "between toe and crest the module lowers the dune by Larson et al. (2004), never below the toe; over "
        "the crest Barrier3D's overwash takes over. Dots: the peak water levels of the 1997-2009 storm series. "
        f"(c) Every 10 m cell's toe in GIS {g} in the {' and '.join(years)} surveys; toes below "
        f"{TOE_MIN_M_NAVD88:g} m NAVD88 or within {TOE_BELOW_CREST_M:g} m of the crest are screened out (x), and "
        "the mean of the rest is the domain's one toe (dashed), held for the run. (d) The domain toes along the "
        "island: 1996-2009 from the 1996 and 1997 surveys, 2009-2025 from the 2004, 2005 and 2009 surveys "
        "(nothing after the start). The module is off at Cape Point (GIS 1-5)."))
    print(f"wrote {png}")
    print(f"GIS {g}: {int(ok.sum())} toes kept, {int((~ok).sum())} screened out, domain toe {dom_toe:.2f} m, "
          f"kept range {d[ok].toe_m_navd88.min():.2f}-{d[ok].toe_m_navd88.max():.2f} m")


def transect_figure():
    """One lidar transect across the example domain and the toe picked on it (2026-11-04): the simple version of
    the method figure's panel (a), with the domain toe the model uses."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from site_layer import hat_figure_style as fs
    fs.apply_style()
    start, g, cell = EXAMPLE
    c = _start_colours()[start]
    years = [str(y) for y in TOE_DOMAIN_SURVEYS[start][1]]
    dom_toe = float(_tables()[start][1].groupby("gis").toe_m_navd88.first()[g])
    survey = years[-1]
    z, r = next((z[cell], crest_and_toe(z, i0, i1, S.X)) for s, _, z, i0, i1 in domain_profiles(start, [g])
                if str(s) == survey)
    sb, ci, ti, fi = (int(r[k][cell]) for k in ("beach_start_idx", "crest_idx", "toe_idx", "foredune_idx"))
    lo, hi = max(sb - 4, 0), min(ci + 6, len(z))
    x = (np.arange(lo, hi) - sb) * 10.0
    zz = z[lo:hi] + MHW_NAVD88_M
    toe, crest, base, top = zz[ti - lo], zz[ci - lo], zz[sb - lo], z[fi] + MHW_NAVD88_M

    fig, ax = plt.subplots(figsize=fs.figsize("single", height=3.0), layout="constrained")
    ax.fill_between(x, zz, 0, color=fs.C["BASE_FILL"], lw=0, zorder=1)
    ax.plot(x, zz, color=fs.INK, lw=1.6, zorder=4)
    ax.plot([0, (fi - sb) * 10], [base, top], color=fs.INK_MUTED, lw=0.9, ls="--", zorder=3)
    ax.plot((ti - sb) * 10, toe, "o", ms=8, color=c, zorder=5)
    ax.annotate(f"dune toe {toe:.2f} m", ((ti - sb) * 10, toe), xytext=(8, -14), textcoords="offset points",
                fontsize=8, color=c, fontweight="bold")
    ax.plot(0, base, "o", ms=5, color=fs.INK, zorder=5)
    ax.annotate("beach start", (0, base), xytext=(6, -12), textcoords="offset points", fontsize=7.5, color=fs.INK)
    ax.plot((ci - sb) * 10, crest, "o", ms=5, color=fs.INK, zorder=5)
    ax.annotate(f"crest {crest:.2f} m", ((ci - sb) * 10, crest), xytext=(-6, 6), textcoords="offset points",
                ha="right", fontsize=7.5, color=fs.INK)
    ax.text(x[0] + 2, 0.25, "ocean", fontsize=7.5, style="italic", color="#4F7CA6")
    ax.set_xlim(x[0], x[-1])
    ax.set_ylim(0, crest + 1.0)
    ax.set_xlabel("Distance landward of the beach start (m)")
    ax.set_ylabel("Elevation (m NAVD88)")
    fs.open_frame(ax)
    png = fs.figure_dir("inputs", "1-domains", "dune_toes") / "dune_toe_transect.png"
    fs.save(fig, png, close=True)
    fs.record_caption(png, (
        f"The dune toe on one lidar transect across GIS {g} (cell {cell}, {survey} survey). The toe (red) is the "
        "point on the profile furthest below the straight line (dashed) from the beach start to the dune: the "
        "break in slope where the beach turns up into the dune face. The model uses one toe per domain, the mean "
        f"of every cell's toe in the {' and '.join(years)} surveys ({dom_toe:.2f} m NAVD88 for GIS {g})."))
    print(f"wrote {png}: toe {toe:.2f} m, crest {crest:.2f} m, domain toe {dom_toe:.2f} m")


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[1],
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--start", choices=list(S.STARTS))
    ap.add_argument("--cells", action="store_true",
                    help="rewrite dune_toe_cells.csv from the saved survey_dune_toes.csv (no profiles read)")
    ap.add_argument("--figures", action="store_true",
                    help="draw the figures from the saved tables of every start built so far")
    args = ap.parse_args()
    if args.figures:
        figures()
        method_figure()
        transect_figure()
        return
    if not args.start:
        ap.error("--start is required unless --figures")
    out = S.tv.product_dir(args.start) / "dune-toes"
    if args.cells:
        t = pd.read_csv(out / "survey_dune_toes.csv")
        t["survey"] = t.survey.astype(str)
    else:
        t = survey_toes(args.start)
    cells = cell_table(t, args.start)
    out.mkdir(parents=True, exist_ok=True)
    if not args.cells:
        t.round(3).to_csv(out / "survey_dune_toes.csv", index=False)
    cells.round(3).to_csv(out / "dune_toe_cells.csv", index=False)
    print(cells.groupby("toe_source").size().to_string())
    print(cells[["toe_m_navd88", "beach_width_m"]].describe().round(2).to_string())
    print(f"wrote {out / 'survey_dune_toes.csv'}\nwrote {out / 'dune_toe_cells.csv'}")


if __name__ == "__main__":
    main()
