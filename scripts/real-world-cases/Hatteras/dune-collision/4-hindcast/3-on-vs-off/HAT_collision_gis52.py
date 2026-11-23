"""
Why GIS 52 does not drop (2026-11-19): the dune crest through 1996-2009 at GIS 50-54, lidar against the model
without the collision module and with it (Lower front row, equal rows, Cs 2e-6).

    python scripts/real-world-cases/Hatteras/dune-collision/4-hindcast/3-on-vs-off/HAT_collision_gis52.py

Writes 3-collision-design/0-all-designs/gis50-54_crest_through_time_1996_2009.png. Each panel prints how far the
module lowered the crest below no collision by the end, so the module's effect is the gap between the two lines.

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-11-19
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

_HERE = Path(__file__).resolve()
sys.path.insert(0, str(next(p for p in _HERE.parents if (p / "hatteras_paths.py").exists())))
sys.path.insert(0, str(_HERE.parent))
import collision_designs as CD  # noqa: E402
import HAT_collision_design_visuals as V  # noqa: E402

H = V.H
GIS = (50, 51, 52, 53, 54)


def main():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from site_layer import hat_figure_style as S
    S.apply_style()
    runs = (("No collision", "cal_step2_growth_noscrape_vdefault", CD.EROSION["nocollision"][2]),
            ("Lower front row, Cs 2e-6", H.DESIGN_SWEEP[("lowerfront_equal", 2)], CD.EROSION["lowerfront"][2]))
    sv = pd.read_csv(H.SURVEY_CRESTS)
    sv = sv[(sv.status == "ok") & sv.gis.isin(GIS)].copy()
    t = pd.to_datetime(sv.flown.astype(str).str[:7], errors="coerce")
    sv["yr"] = t.dt.year + (t.dt.month - 0.5) / 12
    lid = sv.groupby(["gis", "survey"]).agg(z=("crest_m_navd88", "median"), yr=("yr", "first")).reset_index()
    fig, axs = plt.subplots(1, len(GIS), figsize=(S.FIG_W_DOUBLE, 3.0), sharey=True, layout="constrained")
    ends = {}
    for lab, v, col in runs:
        c, d, first = V.load(v)
        for ax, g in zip(axs, GIS):
            b = c.barrier3d[d.gis_to_pad(g)]
            z = (np.asarray(b.DuneDomain, float).max(axis=2) + b._BermEl).mean(axis=1) * 10 + H.MHW_NAVD88_M
            ax.plot(first + np.arange(z.size), z, color=col, lw=2.2, label=lab)
            ends[(lab, g)] = z[-1]
    for ax, g in zip(axs, GIS):
        L = lid[lid.gis == g].sort_values("yr")
        ax.plot(L.yr, L.z, "o", color="0.1", ms=5, label="Lidar crest (domain median)")
        gap = ends[("Lower front row, Cs 2e-6", g)] - ends[("No collision", g)]
        obs = L.z.iloc[-1] - L.z.iloc[0]
        ax.text(0.04, 0.04, f"module: {gap:+.1f} m\nlidar 1996→2009: {obs:+.1f} m", transform=ax.transAxes,
                fontsize=7.5, va="bottom", color="0.15")
        ax.axvline(2003.7, color="0.6", lw=0.7, ls=":")
        ax.set_title(f"GIS {g}", fontsize=9.5, fontweight="bold")
        ax.set_xlim(1996, 2010.5)
        ax.set_xticks([1998, 2002, 2006, 2010])
        ax.tick_params(labelsize=7.5)
        ax.grid(axis="y", color="0.92", lw=0.5)
        ax.spines[["top", "right"]].set_visible(False)
    axs[0].set_ylabel("Dune crest (m NAVD88)", fontsize=8.5)
    axs[0].text(2003.8, 1.02, "Isabel", transform=axs[0].get_xaxis_transform(), fontsize=7, color="0.4")
    h, l = axs[0].get_legend_handles_labels()
    fig.legend(h, l, loc="outside lower center", ncol=3, fontsize=8, frameon=False)
    o = CD.ALL_DIR / "gis50-54_crest_through_time_1996_2009.png"
    S.save(fig, o, bbox_inches="tight", close=True)
    S.record_caption(o, (
        "Dune crest (the higher model row, domain mean, m NAVD88) through 1996-2009 at GIS 50-54: no collision "
        f"(run {runs[0][1]}) and Lower front row, equal rows, Cs 2e-6 (run {runs[1][1]}), against the lidar domain-"
        "median crest at each survey. Text: the module's effect at the end (with minus without) and the lidar's "
        "change from the 1996 to the 2009 survey. Dotted: Isabel."))
    print(f"wrote {o}")


def observed(runs=()):
    """gis50-54_crest_observed_1996_2009.png: the lidar alone, each survey's domain-median crest joined in time,
    the middle half of the domain's cells as a bar, the 1996-2009 change printed. With runs ((label, run, colour),
    ...): each run's domain-median crest (the higher row) added as a line, written as
    gis50-54_crest_observed_and_model_1996_2009.png."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from site_layer import hat_figure_style as S
    S.apply_style()
    sv = pd.read_csv(H.SURVEY_CRESTS)
    sv = sv[(sv.status == "ok") & sv.gis.isin(GIS)].copy()
    t = pd.to_datetime(sv.flown.astype(str).str[:7], errors="coerce")
    sv["yr"] = t.dt.year + (t.dt.month - 0.5) / 12
    lid = sv.groupby(["gis", "survey"]).agg(z=("crest_m_navd88", "median"), lo=("crest_m_navd88", lambda x: x.quantile(.25)),
                                            hi=("crest_m_navd88", lambda x: x.quantile(.75)), yr=("yr", "first")).reset_index()
    fig, axs = plt.subplots(1, len(GIS), figsize=(S.FIG_W_DOUBLE, 3.0), sharey=True, layout="constrained")
    for lab, v, col in runs:
        c, d, first = V.load(v)
        for ax, g in zip(axs, GIS):
            b = c.barrier3d[d.gis_to_pad(g)]
            z = np.median(np.asarray(b.DuneDomain, float).max(axis=2) + b._BermEl, axis=1) * 10 + H.MHW_NAVD88_M
            ax.plot(first + np.arange(z.size), z, color=col, lw=2.0, label=lab, zorder=3)
    for ax, g in zip(axs, GIS):
        L = lid[lid.gis == g].sort_values("yr")
        ax.vlines(L.yr, L.lo, L.hi, color="0.75", lw=3, zorder=1)
        ax.plot(L.yr, L.z, "-o", color="0.1", ms=5, lw=1.2, zorder=4, label="Lidar crest (domain median)")
        ax.text(0.97, 0.97, f"1996→2009: {L.z.iloc[-1] - L.z.iloc[0]:+.1f} m", transform=ax.transAxes, fontsize=7.5,
                color="0.15", ha="right", va="top",
                bbox=dict(boxstyle="round,pad=0.2", fc="white", ec="none", alpha=0.9))
        ax.axvline(2003.7, color="0.6", lw=0.7, ls=":")
        ax.set_title(f"GIS {g}", fontsize=9.5, fontweight="bold")
        ax.set_xlim(1996, 2010.5)
        ax.set_xticks([1998, 2002, 2006, 2010])
        ax.tick_params(labelsize=7.5)
        ax.grid(axis="y", color="0.92", lw=0.5)
        ax.spines[["top", "right"]].set_visible(False)
    axs[0].set_ylabel("Lidar dune crest (m NAVD88)", fontsize=8.5)
    axs[0].text(2003.8, 1.02, "Isabel", transform=axs[0].get_xaxis_transform(), fontsize=7, color="0.4")
    if runs:
        axs[0].set_ylim(top=axs[0].get_ylim()[1] + 0.6)     # headroom for the change label above the lines
        h, l = axs[0].get_legend_handles_labels()
        fig.legend(h, l, loc="outside lower center", ncol=len(h), fontsize=8, frameon=False)
    o = CD.ALL_DIR / ("gis50-54_crest_observed_and_model_1996_2009.png" if runs
                      else "gis50-54_crest_observed_1996_2009.png")
    S.save(fig, o, bbox_inches="tight", close=True)
    S.record_caption(o, (
        "Lidar dune crest at GIS 50-54, 1996-2009: each survey's domain median (dots, joined in time) and the middle "
        "half of the domain's 10 m cells (grey bar, 25th-75th percentile), the highest point of each profile's "
        "dune window (survey_dune_crests.csv). The 2001 survey is unclassified lidar. Dotted: Isabel."
        + ("" if not runs else " Lines: the model's domain-median crest (the higher dune row) at 1 January of each "
           "year, " + "; ".join(f"{lab} (run {v})" for lab, v, _ in runs) + ".")))
    print(f"wrote {o}")


SHORE_WINDOWS = ((1996, 1997), (1998, 1999), (2000, 2001), (2002, 2003), (2004, 2005), (2006, 2007), (2008, 2009))
SHORE_COLOURS = ("#D7301F", "#F28E2B", "#F2C14E", "#1A9850", "#2C7FB8", "#3F2A8C", "#8E44AD")   # ROYGBIV in time
PHOTO_YEARS = (1998, 2004)
ZOOM_GIS = (51, 53)                    # the imagery figure zooms on the drop and its two neighbours
MAP_CRS = "EPSG:26918"
MIN_POSITIONS = 3                      # satellite positions a transect needs in a window for its mean


def imagery():
    """gis50-54_shorelines_on_imagery.png: two-year CoastSat mean shorelines at GIS 50-54 (1996-97 to 2008-09,
    ROYGBIV in time) on the 1998 and 2004 aerial photographs, the nearest flights before and after the dune loss
    at GIS 52 (none 1999-2003); domain boxes outlined. Prints each domain's mean position per window."""
    import geopandas as gpd
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.patheffects as pe
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    _hat = next(p for p in _HERE.parents if (p / "hatteras_paths.py").exists())
    sys.path.insert(0, str(_hat / "input_prep" / "5-scr" / "1-observations" / "mean_shoreline"))
    sys.path.insert(0, str(_hat / "input_prep" / "5-scr" / "lib"))
    import coastsat_mean_shoreline as cms
    from coastsat_lrr import load_timeseries
    from site_layer import hat_figure_style as S
    from site_layer.hat_aerial_imagery import Imagery
    from site_layer.hat_observed_rates import DOMAIN_BOXES, TRANSECT_DOMAINS
    S.apply_style()
    boxes = gpd.read_file(DOMAIN_BOXES).to_crs(MAP_CRS)
    boxes["gis"] = np.arange(1, len(boxes) + 1)          # file order is south -> north, GIS 1-90
    lk = pd.read_csv(TRANSECT_DOMAINS / "transect_domain_lookup.csv").dropna(subset=["domain_number"])
    lk = lk[lk.domain_number.between(min(GIS) - 1, max(GIS) + 1)]
    geom = cms.transect_geometry(set(lk.transect_id))
    rows = []
    for t in lk.itertuples():
        if t.transect_id not in geom:
            continue
        f = cms.timeseries_file(t.transect_id)
        if not f.exists():
            continue
        ts = load_timeseries(str(f)).dropna(subset=["chainage_m"])
        x0, y0, ux, uy = geom[t.transect_id][:4]
        for lo, hi in SHORE_WINDOWS:
            w = ts[ts.date.dt.year.between(lo, hi)]
            if len(w) >= MIN_POSITIONS:
                ch = float(w.chainage_m.mean())
                rows.append(dict(transect_id=t.transect_id, gis=int(t.domain_number), window=f"{lo}–{hi}",
                                 n=len(w), chainage_m=ch, x=x0 + ch * ux, y=y0 + ch * uy))
    pos = pd.DataFrame(rows)
    first = pos[pos.window == f"{SHORE_WINDOWS[0][0]}–{SHORE_WINDOWS[0][1]}"].set_index("transect_id").chainage_m
    pos["change_m"] = pos.chainage_m - pos.transect_id.map(first)
    seg = boxes[boxes.gis.between(*ZOOM_GIS)]
    ylo, yhi = seg.total_bounds[1], seg.total_bounds[3]
    p = pos[pos.y.between(ylo, yhi)]
    b = (p.x.min() - 170.0, ylo, p.x.max() + 90.0, yhi)
    fig, axs = plt.subplots(1, len(PHOTO_YEARS) + 1, figsize=(S.FIG_W_DOUBLE, 7.4), layout="constrained",
                            gridspec_kw={"width_ratios": [1] * len(PHOTO_YEARS) + [1.25]})
    # the change panel: each window's position minus 1996-97, transect by transect, on the maps' alongshore axis
    ac = axs[-1]
    for (lo, hi), col in zip(SHORE_WINDOWS[1:], SHORE_COLOURS[1:]):
        s = pos[(pos.window == f"{lo}–{hi}") & pos.y.between(ylo, yhi)].sort_values("y")
        ac.plot(s.change_m, s.y, color=col, lw=1.8)
    ac.axvline(0, color=SHORE_COLOURS[0], lw=1.8)
    for g, poly in zip(seg.gis, seg.geometry):
        ac.axhline(poly.bounds[1], color="0.75", lw=0.6)
    ac.set_ylim(ylo, yhi)
    ac.set_yticks([(poly.bounds[1] + poly.bounds[3]) / 2 for poly in seg.geometry])
    ac.set_yticklabels([f"GIS {g}" for g in seg.gis], fontsize=8)
    ac.set_xlabel("Shoreline change since 1996–97 (m)\n← landward    seaward →", fontsize=8.5)
    ac.grid(axis="x", color="0.92", lw=0.5)
    ac.spines[["top", "right"]].set_visible(False)
    ac.set_title("Change by window", fontsize=9.5, loc="left")
    for ax, yr in zip(axs[:-1], PHOTO_YEARS):
        img = Imagery(yr, MAP_CRS, strict=False).read(b[0], b[2], b[1], b[3], 1.0, MAP_CRS).copy()
        img[img.max(axis=2) == 0] = 255
        ax.imshow(img, extent=(b[0], b[2], b[1], b[3]), origin="upper", interpolation="bilinear", zorder=0)
        for (lo, hi), col in zip(SHORE_WINDOWS, SHORE_COLOURS):
            s = pos[pos.window == f"{lo}–{hi}"].sort_values("y")
            ax.plot(s.x, s.y, color=col, lw=1.6, zorder=3,
                    path_effects=[pe.withStroke(linewidth=3.0, foreground="white")])
        for g, poly in zip(seg.gis, seg.geometry):
            xs, ys = poly.exterior.xy
            ax.plot(xs, ys, color="white", lw=0.8, zorder=2,
                    path_effects=[pe.withStroke(linewidth=1.8, foreground="0.35")])
            ax.text(b[0] + 0.05 * (b[2] - b[0]), (poly.bounds[1] + poly.bounds[3]) / 2, f"GIS {g}", fontsize=8,
                    va="center", ha="left", zorder=5,
                    bbox=dict(facecolor="white", alpha=0.85, edgecolor="none", boxstyle="square,pad=0.2"))
        ax.set_xlim(b[0], b[2])
        ax.set_ylim(b[1], b[3])
        ax.set_aspect("equal")
        ax.set_xticks([])
        ax.set_yticks([])
        S.spines_for_image(ax)
        ax.set_title(f"Aerial photograph, {yr}", fontsize=9.5, loc="left")
    S._scalebar(axs[0], 100.0, show_cells=False)
    S._north_arrow(axs[0], x=0.82, y=0.06)
    h = [Line2D([], [], color=c, lw=2.0, path_effects=[pe.withStroke(linewidth=3.4, foreground="white")])
         for c in SHORE_COLOURS]
    fig.legend(h, [f"{lo}–{hi}" for lo, hi in SHORE_WINDOWS], loc="outside lower center", ncol=4, fontsize=8,
               frameon=False, title="CoastSat mean shoreline", title_fontsize=8.5)
    o = CD.ALL_DIR / f"gis{ZOOM_GIS[0]}-{ZOOM_GIS[1]}_shorelines_on_imagery.png"
    S.save(fig, o, vector=False, close=True)
    S.record_caption(o, (
        "CoastSat mean shoreline over each two-year window, 1996-97 to 2008-09 (ROYGBIV, red earliest), per "
        f"transect from at least {MIN_POSITIONS} satellite positions, geolocated along each transect, at GIS "
        f"{ZOOM_GIS[0]}-{ZOOM_GIS[1]}, on the {PHOTO_YEARS[0]} and {PHOTO_YEARS[1]} aerial photographs (the nearest flights before and "
        "after the dune loss at GIS 52; none 1999-2003). White outlines: the Barrier3D domains. North up. "
        "Right: each window's position minus 1996-97 per transect, on the maps' alongshore axis. 1996-1999 "
        "positions are from fewer satellite images (Landsat 5 era)."))
    pos.to_csv(o.parent / "supporting" / o.name.replace(".png", ".csv"), index=False)
    print(f"wrote {o}")
    tab = pos[pos.gis.between(min(GIS), max(GIS))].pivot_table(index="window", columns="gis", values="chainage_m")
    print("Domain-mean shoreline position, m seaward along the transects (change from 1996-97):")
    print((tab - tab.iloc[0]).round(1).to_string())
    print(pos[pos.gis.between(min(GIS), max(GIS))].groupby("window").n.median().to_string())


PROFILE_GIS = (51, 52, 53)
PROFILE_SURVEYS = (1996, 1998, 1999, 2000, 2001, 2004, 2009)          # ROYGBIV in time
PROFILE_XLIM = (60, 330)                                              # cross-shore range shown, m from the frame's ocean edge
PROFILE_WATER_M = 0.0                                                 # below this (m NAVD88) is water or swash: not drawn


SIMPLE_GIS = 52
SIMPLE_SURVEYS = (1996, 2000, 2001, 2004, 2009)
SIMPLE_COLOURS = ("#D7301F", "#1A9850", "#2C7FB8", "#3F2A8C", "#8E44AD")


def profiles(simple=False):
    """gis51-53_lidar_profiles_1996_2009.png: the lidar cross-shore profile at each survey, the median over each
    domain's 50 alongshore 10 m cells, ocean on the left, in the 1996 start's frame (the same straightened grid
    the crests are picked on), so a retreating beach moves right. simple: one domain (SIMPLE_GIS), five surveys,
    the beach, foredune and landward ridge labelled: gis52_lidar_profile_simple_1996_2009.png."""
    global PROFILE_GIS, PROFILE_SURVEYS, SHORE_COLOURS
    if simple:
        PROFILE_GIS, PROFILE_SURVEYS, SHORE_COLOURS = (SIMPLE_GIS,), SIMPLE_SURVEYS, SIMPLE_COLOURS
    import rasterio
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    _hat = next(p for p in _HERE.parents if (p / "hatteras_paths.py").exists())
    sys.path.insert(0, str(_hat / "input_prep" / "1-barrier3d-domains" / "2-survey-crests"))
    import HAT_survey_dune_crests as SC
    import json
    from site_layer import hat_figure_style as S
    S.apply_style()
    X, tv, start = SC.X, SC.tv, "1996-start"
    prod, surveys = SC.STARTS[start]
    windows = json.loads(tv.picks_dir(start).joinpath(
        f"HAT_dune_search_windows_{tv.topo_dirs(start)[2]}.json").read_text(encoding="utf-8"))
    frames = {g: SC.frame(start, g, windows) for g in PROFILE_GIS}
    grids = {g: SC.clip_grid(SC.product(prod).gapfill_1m / f"clip_domain_{g}_filled.tif") for g in range(1, 91)}
    prof = {}
    for year in PROFILE_SURVEYS:
        kind, path = surveys[year][:2]
        files = sorted(path.parent.glob(path.name)) if "*" in path.name else [path]
        srcs = [rasterio.open(f) for f in files]
        try:
            unit = "m" if kind == "mosaic" else SC.unit_of(srcs, start, grids)[0]
            scale = SC.US_FT if unit == "ft" else 1.0
            shift = SC.SURVEY_VERTICAL_OFFSET_M.get(start, {}).get(year, 0.0)
            for g in PROFILE_GIS:
                shear, c0, i0, i1 = frames[g]
                on = SC.mosaic_on_clip(srcs[0], grids[g]) if kind == "mosaic" else SC.raster_on_clip(srcs, grids[g])
                z = SC.profiles(SC.downsample(on * scale + shift), shear, c0)[:50]
                z = np.where(z <= X.NODATA_SENTINEL_M + 1e-9, np.nan, z)
                z = np.where(z <= X.SENTINEL_WATER_M + 1e-9, np.nan, z) + X.MHW_M          # m NAVD88
                prof[(year, g)] = np.nanmedian(z, axis=0)
        finally:
            for s in srcs:
                s.close()
        print(f"{year} read", flush=True)
    fig, axs = plt.subplots(1, len(PROFILE_GIS), figsize=(S.FIG_W_SINGLE * 1.5, 3.6) if simple else
                            (S.FIG_W_DOUBLE, 3.2), sharey=True, layout="constrained")
    axs = np.atleast_1d(axs)
    for ax, g in zip(axs, PROFILE_GIS):
        for yr, col in zip(PROFILE_SURVEYS, SHORE_COLOURS):
            z = prof[(yr, g)][: PROFILE_XLIM[1] // 10 + 1]
            z = np.where(z < PROFILE_WATER_M, np.nan, z)
            ax.plot(np.arange(z.size) * 10, z, color=col, lw=2.0 if g == 52 else 1.6, label=str(yr))
        ax.axhline(X.MHW_M, color="0.45", lw=0.7, ls="--")
        ax.axhline(X.BERM_ELEV_NAVD_M, color="0.45", lw=0.7, ls=":")
        ax.set_title(f"GIS {g}", fontsize=9.5, fontweight="bold")
        ax.set_xlabel("Cross-shore distance (m), ocean → island", fontsize=8.5)
        ax.set_xlim(*((90, 300) if simple else PROFILE_XLIM))
        ax.grid(color="0.93", lw=0.5)
        if simple:                  # name the parts of the profile
            note = dict(fontsize=8, color="0.15", ha="center")
            arrow = dict(arrowstyle="-", color="0.4", lw=0.7)
            ax.annotate("foredune\n6 m in 1996,\ngone by 2001", xy=(150, 5.9), xytext=(195, 5.6), va="center",
                        arrowprops=arrow, **note)
            ax.annotate("landward ridge\n(unchanged to 2004)", xy=(240, 3.45), xytext=(272, 4.6), va="center",
                        arrowprops=arrow, **note)
            ax.annotate("beach", xy=(112, 0.9), xytext=(112, 2.6), va="bottom", arrowprops=arrow, **note)
            ax.annotate("", xy=(132, 0.55), xytext=(108, 0.55),
                        arrowprops=dict(arrowstyle="->", color="0.15", lw=1.1))
            ax.text(137, 0.55, "beach moved ~30 m landward", fontsize=7.5, color="0.15", ha="left", va="center")
        ax.spines[["top", "right"]].set_visible(False)
    xr = 300 if simple else PROFILE_XLIM[1]
    axs[0].text(xr - 2, X.MHW_M + 0.1, "MHW", fontsize=7, color="0.4", ha="right", va="bottom")
    axs[0].text(xr - 2, X.BERM_ELEV_NAVD_M + 0.1, "berm / model toe", fontsize=7, color="0.4", ha="right",
                va="bottom")
    axs[0].set_ylabel("Elevation (m NAVD88)", fontsize=8.5)
    h, l = axs[0].get_legend_handles_labels()
    fig.legend(h, l, loc="outside lower center", ncol=len(h), fontsize=8, frameon=False, title="Lidar survey",
               title_fontsize=8.5)
    o = CD.ALL_DIR / (f"gis{SIMPLE_GIS}_lidar_profile_simple_1996_2009.png" if simple
                      else "gis51-53_lidar_profiles_1996_2009.png")
    S.save(fig, o, bbox_inches="tight", close=True)
    S.record_caption(o, (
        f"Lidar cross-shore profiles at GIS {PROFILE_GIS[0]}-{PROFILE_GIS[-1]}, the median over each domain's 50 "
        "alongshore 10 m cells at each survey (ROYGBIV in time), resampled to 10 m and placed in the 1996 start's "
        "straightened frame (the grid the dune crests are picked on), ocean on the left; water, swash below 0 m "
        "NAVD88 and unsurveyed cells left out. Dashed: MHW (0.36 m); dotted: the berm, the model's dune toe (1.70 m). The 1996-2005 surveys are "
        "coarse grids that smooth sharp crests; 2001 is unclassified lidar."))
    print(f"wrote {o}")


BEACH_SURVEYS = (1996, 2000, 2001, 2004, 2009)
BEACH_COLOURS = ("#D7301F", "#1A9850", "#2C7FB8", "#3F2A8C", "#8E44AD")   # the profile figure's colours
SMOOTH_CELLS = 5                                                          # alongshore running median of the lines


def imagery_beach():
    """gis51-53_beach_on_imagery.png: on the 1998 and 2004 photographs, each survey's foredune crest (solid; the
    highest point of the survey's dune window, as the crests are picked) and shoreline (dashed; the first cell above
    0.5 m MHW, the beach start the windows are clipped to), per 10 m alongshore cell, placed on the map by sending
    the clip's own coordinates through the same orient-straighten-trim steps as the elevations."""
    import rasterio
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.patheffects as pe
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    import json
    _hat = next(p for p in _HERE.parents if (p / "hatteras_paths.py").exists())
    sys.path.insert(0, str(_hat / "input_prep" / "1-barrier3d-domains" / "2-survey-crests"))
    import HAT_survey_dune_crests as SC
    from site_layer import hat_figure_style as S
    from site_layer.hat_aerial_imagery import Imagery
    S.apply_style()
    X, tv, start = SC.X, SC.tv, "1996-start"
    prod, surveys = SC.STARTS[start]
    windows = json.loads(tv.picks_dir(start).joinpath(
        f"HAT_dune_search_windows_{tv.topo_dirs(start)[2]}.json").read_text(encoding="utf-8"))
    frames = {g: SC.frame(start, g, windows) for g in PROFILE_GIS}
    grids = {g: SC.clip_grid(SC.product(prod).gapfill_1m / f"clip_domain_{g}_filled.tif") for g in range(1, 91)}
    crs = grids[PROFILE_GIS[0]][2].to_string()
    # map coordinates of every profile cell, through the elevations' own steps
    xy = {}
    for g in PROFILE_GIS:
        t, (h, w), _ = grids[g]
        cols, rows = np.meshgrid(np.arange(w) + 0.5, np.arange(h) + 0.5)
        shear, c0, *_ = frames[g]
        xy[g] = tuple(SC.profiles(SC.downsample(t.c + cols * t.a + rows * t.b), shear, c0)[:50] if k == 0 else
                      SC.profiles(SC.downsample(t.f + cols * t.d + rows * t.e), shear, c0)[:50] for k in (0, 1))
    lines = {}
    for year in BEACH_SURVEYS:
        kind, path = surveys[year][:2]
        files = sorted(path.parent.glob(path.name)) if "*" in path.name else [path]
        srcs = [rasterio.open(f) for f in files]
        try:
            unit = "m" if kind == "mosaic" else SC.unit_of(srcs, start, grids)[0]
            scale = SC.US_FT if unit == "ft" else 1.0
            shift = SC.SURVEY_VERTICAL_OFFSET_M.get(start, {}).get(year, 0.0)
            for g in PROFILE_GIS:
                shear, c0, i0, i1 = frames[g]
                on = SC.mosaic_on_clip(srcs[0], grids[g]) if kind == "mosaic" else SC.raster_on_clip(srcs, grids[g])
                z = SC.profiles(SC.downsample(on * scale + shift), shear, c0)[:50]
                i0s, i1s = SC.survey_windows(start, year, g, i0, i1)
                _, _, status, idx, _ = SC.crests(z, i0s, i1s, with_index=True)
                above = z > X.BEACH_START_THR_M
                sb = np.where(above.any(axis=1), above.argmax(axis=1), -1)
                gx, gy = xy[g]
                for name, k, ok in (("crest", idx, (status == "ok") & (idx >= 0)), ("shore", sb, sb >= 0)):
                    kk = np.clip(k, 0, gx.shape[1] - 1)
                    px = np.where(ok, gx[np.arange(50), kk], np.nan)
                    py = np.where(ok, gy[np.arange(50), kk], np.nan)
                    lines.setdefault((year, name), []).append(pd.DataFrame({"x": px, "y": py}))
        finally:
            for s in srcs:
                s.close()
        print(f"{year} read", flush=True)
    for k, v in list(lines.items()):          # one line per survey and feature, south to north, smoothed
        d = pd.concat(v, ignore_index=True).dropna().sort_values("y")
        d = d.rolling(SMOOTH_CELLS, center=True, min_periods=2).median()          # drop single-cell jumps
        lines[k] = d.rolling(SMOOTH_CELLS + 4, center=True, min_periods=2).mean().dropna()   # soften 10 m steps
    allp = pd.concat(lines.values())
    b = (allp.x.min() - 120.0, allp.y.min(), allp.x.max() + 90.0, allp.y.max())
    fig, axs = plt.subplots(1, len(PHOTO_YEARS), figsize=(S.FIG_W_DOUBLE * 0.72, 7.6), layout="constrained")
    # before / after: each photo with the survey nearest it, the other drawn faintly on the after panel
    pairs = ((PHOTO_YEARS[0], (BEACH_SURVEYS[0],)), (PHOTO_YEARS[1], (BEACH_SURVEYS[0], PHOTO_YEARS[1])))
    colour = dict(zip(BEACH_SURVEYS, BEACH_COLOURS))
    for ax, (yr, show) in zip(axs, pairs):
        img = Imagery(yr, crs, strict=False).read(b[0], b[2], b[1], b[3], 0.8, crs).copy()
        img[img.max(axis=2) == 0] = 255
        ax.imshow(img, extent=(b[0], b[2], b[1], b[3]), origin="upper", interpolation="bilinear", zorder=0)
        for year in show:
            faint = len(show) > 1 and year == show[0]
            for name, ls in (("shore", (0, (4, 2))), ("crest", "-")):
                d = lines[(year, name)]
                ax.plot(d.x, d.y, color=colour[year], lw=1.3 if faint else 2.2, ls=ls, zorder=2 if faint else 3,
                        alpha=0.55 if faint else 1.0,
                        path_effects=[pe.withStroke(linewidth=2.6 if faint else 3.8, foreground="white")])
        ax.set_xlim(b[0], b[2])
        ax.set_ylim(b[1], b[3])
        ax.set_aspect("equal")
        ax.set_xticks([])
        ax.set_yticks([])
        S.spines_for_image(ax)
        ax.set_title(f"{'Before' if yr == PHOTO_YEARS[0] else 'After'}: {yr} photo, "
                     + " + ".join(f"{y} lidar" for y in show[::-1]), fontsize=9, loc="left")
    S._scalebar(axs[0], 100.0, show_cells=False)
    S._north_arrow(axs[0], x=0.82, y=0.06)
    halo = [pe.withStroke(linewidth=3.4, foreground="white")]
    used = sorted({y for _, s in pairs for y in s})
    h = [Line2D([], [], color=colour[y], lw=2.2, path_effects=halo) for y in used]
    h += [Line2D([], [], color="0.2", lw=2.0, ls="-"), Line2D([], [], color="0.2", lw=2.0, ls=(0, (4, 2)))]
    fig.legend(h, [f"{y} lidar" for y in used] + ["foredune crest", "shoreline (where the beach starts)"],
               loc="outside lower center", ncol=2, fontsize=8, frameon=False)
    fig.suptitle(f"GIS {PROFILE_GIS[0]}–{PROFILE_GIS[-1]}", fontsize=10, fontweight="bold")
    o = CD.ALL_DIR / f"gis{PROFILE_GIS[0]}-{PROFILE_GIS[-1]}_beach_on_imagery.png"
    S.save(fig, o, vector=False, close=True)
    S.record_caption(o, (
        f"Lidar foredune crest (solid; the highest point in each survey's dune window, as survey_dune_crests.csv "
        f"picks it) and shoreline (dashed; the first 10 m cell above 0.5 m MHW, where the beach starts) at GIS "
        f"{PROFILE_GIS[0]}-{PROFILE_GIS[-1]} for the {', '.join(map(str, BEACH_SURVEYS))} surveys, a "
        f"{SMOOTH_CELLS}-cell alongshore running median, on the {PHOTO_YEARS[0]} and {PHOTO_YEARS[1]} aerial "
        "photographs. North up; positions on the 10 m profile grid, so accurate to about 10 m."))
    print(f"wrote {o}")
    for name in ("shore", "crest"):
        base = lines[(BEACH_SURVEYS[0], name)]
        print(name, {y: round(float(np.median(lines[(y, name)].x) - np.median(base.x)), 1) for y in BEACH_SURVEYS},
              "(median x shift vs 1996, m east = seaward)")


MODEL_RUNS = (("No collision, equal rows", "cal_step2_growth_noscrape_vdefault", CD.EROSION["nocollision"][2]),
              ("Lower front row, equal rows, Cs 2e-5", H.DESIGN_SWEEP[("lowerfront_equal", 20)],
               CD.EROSION["lowerfront"][2]))


if __name__ == "__main__":
    a = sys.argv[1:]
    if a == ["observed"]:
        observed()
    elif a == ["observed_model"]:
        observed(MODEL_RUNS)
    elif a == ["imagery"]:
        imagery()
    elif a == ["profiles"]:
        profiles()
    elif a == ["profile_simple"]:
        profiles(simple=True)
    elif a == ["beach"]:
        imagery_beach()
    else:
        main()
