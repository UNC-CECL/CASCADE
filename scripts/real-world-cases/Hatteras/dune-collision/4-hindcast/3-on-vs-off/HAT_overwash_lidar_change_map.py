"""
Lidar elevation change on the island map beside the model's overwash (2026-11-19): the DEM difference between two
surveys, at 10 m, next to a run's overwash per GIS domain over the same model years, both on the island map of
HAT_overwash_maps_side_by_side.py.

    python scripts/real-world-cases/Hatteras/dune-collision/4-hindcast/3-on-vs-off/HAT_overwash_lidar_change_map.py
    python .../HAT_overwash_lidar_change_map.py grid_faceretreat_equal_cs3_vdefault
    python .../HAT_overwash_lidar_change_map.py model [run]   # the same map for a run: model_change_map_<run>_1996_2009.png
    python .../HAT_overwash_lidar_change_map.py observed      # the lidar alone, 1996 -> 2009: lidar_change_map_observed_1996_2009.png

Lidar: each survey read onto every GIS domain's 1 m clip grid (as HAT_survey_dune_crests.py reads it, with its
vertical offsets), brought to 10 m the same way, and differenced; water (below WATER_NAVD_M in either survey) left
out. The island-wide median change (the datum offset between the surveys) is taken off, and changes smaller than
SHOW_M are left clear, so the land shows through. Blue: ground gained (sand added, or vegetation grown); red:
ground lost. Model: the run's overwash (QowTS) summed over the model years between the two surveys' model states.
Pairs: 2000 -> 2004 (Isabel) and 2005 -> 2009. The lidar grids are cached in the scratch folder below.

Writes 3-collision-design/0-all-designs/overwash_lidar_change_map_<run>_<y0>_<y1>.png.

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-11-19
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

_HERE = Path(__file__).resolve()
_HAT = next(p for p in _HERE.parents if (p / "hatteras_paths.py").exists())
for p in (_HAT, _HERE.parent, _HAT / "input_prep" / "1-barrier3d-domains" / "2-survey-crests",
          _HAT / "input_prep" / "8-overwash-analysis" / "2-record",
          _HAT / "input_prep" / "8-overwash-analysis" / "1-observations"):
    sys.path.insert(0, str(p))
import collision_designs as CD  # noqa: E402
import HAT_overwash_lidar_deposition as LD  # noqa: E402

H = LD.H

# --- CONFIG ------------------------------------------------------------------
START = "1996-start"
PAIRS = ((2000, 2004), (2005, 2009))
DEFAULT_RUN = "grid_lowerfront_equal_cs200_vdefault"
WATER_NAVD_M = 0.3                   # ground below this in either survey: water or wet sand, left out
SHOW_M = 0.15                        # changes smaller than this left clear
LIMIT_M = 1.0                        # colour scale, +/- m
MODEL_MAX = 1000.0                   # m3/m, top of the model scale (as the side-by-side maps)
CACHE = CD.ALL_DIR / "supporting" / "lidar_change_10m"
CAL_SURVEYS = (1996, 1997, 1998, 1999, 2000, 2001, 2004, 2005, 2009)   # the calibration period's surveys
ZOOMS = ((1, 12), (30, 40), (76, 88))  # GIS reaches drawn close up
ZOOM_NAMES = {(1, 12): "Cape Point–Buxton", (30, 40): "North of Avon", (76, 88): "Rodanthe–S-curves"}
# -----------------------------------------------------------------------------


def survey_grids(year):
    """{gis: (z at 10 m, m NAVD88, extent)} for one survey on the domains' clip grids; cached."""
    f = CACHE / f"survey_{year}_10m.npz"
    if f.exists():
        d = np.load(f)
        return {g: (d[f"z{g}"], tuple(d[f"e{g}"])) for g in range(1, 91)}
    import rasterio
    import HAT_survey_dune_crests as SC
    prod, surveys = SC.STARTS[START]
    kind, path = surveys[year][:2]
    files = sorted(path.parent.glob(path.name)) if "*" in path.name else [path]
    grids = {g: SC.clip_grid(SC.product(prod).gapfill_1m / f"clip_domain_{g}_filled.tif") for g in range(1, 91)}
    srcs = [rasterio.open(x) for x in files]
    out = {}
    try:
        unit = "m" if kind == "mosaic" else SC.unit_of(srcs, START, grids)[0]
        scale = SC.US_FT if unit == "ft" else 1.0
        shift = SC.SURVEY_VERTICAL_OFFSET_M.get(START, {}).get(year, 0.0)
        for g in range(1, 91):
            t, (h, w), _ = grids[g]
            on = SC.mosaic_on_clip(srcs[0], grids[g]) if kind == "mosaic" else SC.raster_on_clip(srcs, grids[g])
            h10, w10 = h // SC.BLOCK * SC.BLOCK, w // SC.BLOCK * SC.BLOCK
            z = SC.downsample((on * scale + shift)[:h10, :w10])
            ext = (t.c, t.c + w10 * t.a, t.f + h10 * t.e, t.f)       # left, right, bottom, top
            out[g] = (z, ext)
    finally:
        for s in srcs:
            s.close()
    CACHE.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(f, **{f"z{g}": v[0] for g, v in out.items()}, **{f"e{g}": np.array(v[1]) for g, v in out.items()})
    print(f"  {year} read", flush=True)
    return out


def main(run=DEFAULT_RUN):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.cm import ScalarMappable
    from matplotlib.colors import ListedColormap, Normalize, PowerNorm
    from overwash_map_periods import draw_island, load_geometry, scalebar_and_north
    from HAT_collision_design_cs_sweep import cs_label
    from site_layer import hat_figure_style as S
    S.apply_style()

    label = next((d.name for d in CD.DESIGNS if d.run == run), run)
    for (key, c), v in H.DESIGN_SWEEP.items():
        if v == run:
            e, s = key.split("_")
            label = f"{CD.EROSION[e][0]}, {CD.START[s][0].lower()}, Cs {cs_label(c)}"
    dom, land, road, bounds = load_geometry()
    cmap_m = ListedColormap(plt.get_cmap("YlOrRd")(np.linspace(0.10, 0.80, 256)))
    norm_m = PowerNorm(0.5, vmin=0.0, vmax=MODEL_MAX, clip=True)
    # red loss, blue gain
    lo = plt.get_cmap("Reds")(np.linspace(0.85, 0.20, 128))     # red loss
    hi = plt.get_cmap("Blues")(np.linspace(0.20, 0.85, 128))    # blue gain
    cmap_l = ListedColormap(np.vstack([lo, hi]))
    norm_l = Normalize(-LIMIT_M, LIMIT_M)
    grids = {y: survey_grids(y) for y in sorted({y for p in PAIRS for y in p})}

    for y0, y1 in PAIRS:
        diff = {}
        for g in range(1, 91):
            (a, ext), (b, _) = grids[y0][g], grids[y1][g]
            d = b - a
            d[(a < WATER_NAVD_M) | (b < WATER_NAVD_M)] = np.nan
            diff[g] = (d, ext)
        offset = float(np.nanmedian(np.concatenate([d.ravel() for d, _ in diff.values()])))
        offset = 0.0 if abs(offset) < 0.005 else offset
        ow, span = LD.model_overwash(run, y0, y1)
        fill_m = {g: (cmap_m(norm_m(v)) if v > 0 else "#ececec") for g, v in zip(range(1, 91), ow)}

        fig, (ax_m, ax_l) = plt.subplots(1, 2, figsize=(S.FIG_W_DOUBLE * 0.8, 8.4), layout="constrained")
        draw_island(ax_m, dom, land, road, bounds, fill_m, "a", f"Model, {span[0]}–{span[1]}",
                    reach_labels=False, pad_w=1500.0)
        draw_island(ax_l, dom, land, road, bounds, {}, "b", f"Lidar, {y0}→{y1}",
                    reach_labels=False, pad_w=1500.0)
        for d, (l, r, bt, tp) in diff.values():
            d = d - offset
            d = np.where(np.abs(d) < SHOW_M, np.nan, d)
            ax_l.imshow(d, extent=(l, r, bt, tp), cmap=cmap_l, norm=norm_l, interpolation="nearest",
                        zorder=3.5)
        road.plot(ax=ax_l, color=S.C["ROAD"], lw=0.5, zorder=4)
        for a, t in ((ax_m, f"(a) Model, {span[0]}–{span[1]}"), (ax_l, f"(b) Lidar, {y0} → {y1}")):
            a.set_title("", loc="center")
            a.set_title(t, loc="left", fontsize=10)
        x0, x1 = ax_m.get_xlim()
        ax_m.set_xlim(x0 + 0.44 * (x1 - x0), x1 + 0.44 * (x1 - x0))
        scalebar_and_north(ax_m)
        ax_m.set_xlim(x0, x1)
        cb = fig.colorbar(ScalarMappable(norm=norm_m, cmap=cmap_m), ax=ax_m, orientation="horizontal",
                          location="bottom", shrink=0.85, pad=0.01, aspect=22)
        cb.set_ticks([0, 50, 100, 250, 500, 1000])
        cb.set_label("m³ per m of shore (square-root scale); grey: none", fontsize=7.5)
        cb.ax.tick_params(labelsize=7)
        cb2 = fig.colorbar(ScalarMappable(norm=norm_l, cmap=cmap_l), ax=ax_l, orientation="horizontal",
                           location="bottom", shrink=0.85, pad=0.01, aspect=22, extend="both")
        cb2.set_label(f"ground change, m (blue gain, red loss; clear < {SHOW_M:g} m)", fontsize=7.5)
        cb2.ax.tick_params(labelsize=7)
        fig.suptitle(f"{label}\nlidar datum offset {offset:+.2f} m removed", fontsize=9.5)
        o = CD.ALL_DIR / f"overwash_lidar_change_map_{run}_{y0}_{y1}.png"
        S.save(fig, o, vector=False, close=True)
        S.record_caption(o, (
            f"(a) The run's overwash per GIS domain over model years {span[0]}-{span[1]} (run {run}), yellow to red "
            f"on a square-root scale, grey where none. (b) The lidar change in ground height from the {y0} to the "
            f"{y1} survey at 10 m on the domains' clip grids, after taking off the island-wide median change "
            f"({offset:+.2f} m, the datum offset between the surveys); blue gain, red loss, clear where the "
            f"change is under {SHOW_M:g} m or either survey is below {WATER_NAVD_M:g} m NAVD88. The change also "
            "holds vegetation growth and loss, shoreline and dune change, road clearing and building."))
        print(f"{y0} -> {y1}: offset {offset:+.3f} m; wrote {o}")

        # close up: each reach as a model | lidar pair, colour bars in their own row
        fig = plt.figure(figsize=(S.FIG_W_DOUBLE, 3.9), layout="constrained")
        gs = fig.add_gridspec(2, 2 * len(ZOOMS), height_ratios=[1, 0.06])
        for k, (g0, g1) in enumerate(ZOOMS):
            am, al = fig.add_subplot(gs[0, 2 * k]), fig.add_subplot(gs[0, 2 * k + 1])
            draw_island(am, dom, land, road, bounds, fill_m, "a", "", reach_labels=False, pad_w=0.0, label_every=2)
            draw_island(al, dom, land, road, bounds, {}, "a", "", reach_labels=False, pad_w=0.0, label_every=2)
            for g in range(max(1, g0 - 1), min(90, g1 + 1) + 1):
                d, (l, r, bt, tp) = diff[g]
                d = d - offset
                d = np.where(np.abs(d) < SHOW_M, np.nan, d)
                al.imshow(d, extent=(l, r, bt, tp), cmap=cmap_l, norm=norm_l, interpolation="nearest", zorder=3.5)
            road.plot(ax=al, color=S.C["ROAD"], lw=0.6, zorder=4)
            sub = dom[dom.domain_id.between(g0, g1)].total_bounds
            for i, (a, t) in enumerate(((am, f"Model\nGIS {g0}–{g1}"), (al, f"Lidar\nGIS {g0}–{g1}"))):
                a.set_xlim(sub[0] - 300, sub[2] + 1100)
                a.set_ylim(sub[1] - 200, sub[3] + 200)
                a.set_title("", loc="left")
                a.set_title(f"({chr(97 + 2 * k + i)}) {t}", fontsize=8, loc="left")
                for txt in a.texts:
                    txt.set_clip_on(True)
        cm_ax = fig.add_subplot(gs[1, :len(ZOOMS)])
        cl_ax = fig.add_subplot(gs[1, len(ZOOMS):])
        cb = fig.colorbar(ScalarMappable(norm=norm_m, cmap=cmap_m), cax=cm_ax, orientation="horizontal",
                          ticks=[0, 50, 100, 250, 500, 1000])
        cb.set_label("model overwash (a, c, e), m³/m, square-root scale", fontsize=7.5)
        cb2 = fig.colorbar(ScalarMappable(norm=norm_l, cmap=cmap_l), cax=cl_ax, orientation="horizontal",
                           extend="both")
        cb2.set_label(f"lidar change (b, d, f), m; blue gain, red loss, clear < {SHOW_M:g} m", fontsize=7.5)
        for c_ in (cb, cb2):
            c_.ax.tick_params(labelsize=7)
        fig.suptitle(f"{label}: model years {span[0]}–{span[1]}; lidar {y0} → {y1}, "
                     f"datum offset {offset:+.2f} m removed", fontsize=9)
        o = CD.ALL_DIR / f"overwash_lidar_change_map_{run}_{y0}_{y1}_zoom.png"
        S.save(fig, o, vector=False, close=True)
        S.record_caption(o, "Close-ups of the map of the same name without _zoom: Cape Point to Buxton (GIS 1-12), "
                            "Avon (30-40) and Rodanthe (76-88), each as the model's overwash beside the lidar change.")
        print(f"wrote {o}")


def lidar_change(y0, y1):
    """{gis: (change m, extent)} from survey y0 to y1 with the island-wide median (the datum offset) taken off and
    water left out (not yet thresholded), and the offset."""
    g0s, g1s = survey_grids(y0), survey_grids(y1)
    diff = {}
    for g in range(1, 91):
        (a, ext), (b, _) = g0s[g], g1s[g]
        d = b - a
        d[(a < WATER_NAVD_M) | (b < WATER_NAVD_M)] = np.nan
        diff[g] = (d, ext)
    off = float(np.nanmedian(np.concatenate([d.ravel() for d, _ in diff.values()])))
    off = 0.0 if abs(off) < 0.005 else off
    return {g: (d - off, e) for g, (d, e) in diff.items()}, off


def model_states(run, y0, y1):
    """{gis: (model ground at survey y0's state, at survey y1's state)} in m NAVD88 on the domains' 10 m clip grids.
    The model's grid goes back on the map through the start's profile frame (HAT_survey_dune_crests.frame): its
    column i is the domain's profile i, its interior row 0 the cell behind the 1996 dune crest on that profile
    (HAT_dune_topo_extractor.build_interior), its two dune rows the crest cell and the one seaward of it, and the
    whole column steps landward with the run's shoreline (x_s_TS, one cell per dam, as the design visuals do).
    The cells the dunes leave behind as the shoreline retreats are beach, at the berm."""
    import pandas as pd
    import HAT_survey_dune_crests as SC
    import HAT_collision_design_visuals as V
    X, tv = SC.X, SC.tv
    c, d, first = V.load(run)
    st = pd.read_csv(H.SURVEY_CRESTS, usecols=["survey", "model_state"]).drop_duplicates("survey").set_index("survey")
    ts = [pd.Timestamp(st.model_state[y]).year - first for y in (y0, y1)]     # DomainTS[k]: after model year first+k-1
    windows = json.loads(tv.picks_dir(START).joinpath(
        f"HAT_dune_search_windows_{tv.topo_dirs(START)[2]}.json").read_text(encoding="utf-8"))
    out = {}
    for g in range(1, 91):
        shear, c0, i0, i1 = SC.frame(START, g, windows)
        raw = np.load(tv.npy_dirs(START)[0] / f"domain_{g}.npy").astype(float)
        P = SC.profiles(raw, shear, c0)
        above = P > X.BEACH_START_THR_M
        sb = np.where(above.any(axis=1), above.argmax(axis=1), -1)
        _, dune_loc = X.find_dunes(P, sb, i0, i1)
        b = c.barrier3d[d.gis_to_pad(g)]
        xs = np.asarray(b.x_s_TS, float)
        berm = b._BermEl * 10 + H.MHW_NAVD88_M
        n, m = raw.shape
        zs = []
        for t in ts:
            col = np.vstack([(np.asarray(b.DuneDomain, float)[t] + b._BermEl).T,
                             np.asarray(b.DomainTS[t], float)]) * 10 + H.MHW_NAVD88_M     # (rows, 50) m NAVD88
            k = int(round(xs[t] - xs[0]))
            z = np.full((n, m), np.nan)
            for i in range(min(n, col.shape[1])):
                if dune_loc[i] < 0:
                    continue
                j0 = dune_loc[i] - 1
                js = np.arange(j0 + k, j0 + k + col.shape[0])
                vals = col[:, i]
                if k > 0:                                     # the beach left behind
                    js = np.concatenate([np.arange(j0, j0 + k), js])
                    vals = np.concatenate([np.full(k, berm), vals])
                cols_ = (m - 1) - (js + c0 + int(shear[i]))   # frame -> the 10 m clip grid (ocean right, flipped)
                ok = (cols_ >= 0) & (cols_ < m)
                z[n - 1 - i, cols_[ok]] = vals[ok]
            z[z < WATER_NAVD_M] = np.nan
            zs.append(z)
        out[g] = tuple(zs)
    return out


def summary_map(diff, title, cbar_label, o, caption, crop=None, hatch_villages=True):
    """One island map of a change grid with the three close-ups beside it; crop: {gis: grid} whose data extent sets
    each close-up's width (default the change itself), so two maps can share their frames."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import matplotlib.patches as mpatches
    from matplotlib.cm import ScalarMappable
    from matplotlib.colors import ListedColormap, Normalize
    from overwash_map_periods import draw_island, load_geometry, scalebar_and_north
    from overwash_data import SECTIONS
    from site_layer import hat_figure_style as S
    S.apply_style()
    dom, land, road, bounds = load_geometry()
    lo = plt.get_cmap("Reds")(np.linspace(0.85, 0.20, 128))     # red loss
    hi = plt.get_cmap("Blues")(np.linspace(0.20, 0.85, 128))    # blue gain
    cmap_l = ListedColormap(np.vstack([lo, hi]))
    norm_l = Normalize(-LIMIT_M, LIMIT_M)
    crop = crop or {g: d for g, (d, _) in diff.items()}
    shown = {g: (np.where(np.abs(d) < SHOW_M, np.nan, d), e) for g, (d, e) in diff.items()}
    villages = [(lo_, hi_) for _, lo_, hi_, kind in SECTIONS if kind == "village"]
    plt.rcParams["hatch.linewidth"] = 0.6
    developed = []                              # the village domains' boxes, clipped to land
    if hatch_villages:
        for _, r in dom.iterrows():
            if any(lo_ <= int(r.domain_id) <= hi_ for lo_, hi_ in villages):
                part = r.geometry.intersection(land)
                developed += [q for q in getattr(part, "geoms", [part])
                              if q.geom_type == "Polygon" and not q.is_empty]

    def paint(ax, gis=range(1, 91), lw=0.5):
        for g in gis:
            d, (l, r, bt, tp) = shown[g]
            ax.imshow(d, extent=(l, r, bt, tp), cmap=cmap_l, norm=norm_l, interpolation="nearest", zorder=3.5)
        road.plot(ax=ax, color=S.C["ROAD"], lw=lw, zorder=4)
        for q in developed:                     # developed: change there is mostly buildings and roads
            ax.add_patch(mpatches.Polygon(np.asarray(q.exterior.coords), closed=True, fill=False, hatch="////",
                                          edgecolor="0.35", lw=0.0, zorder=5))

    fig = plt.figure(figsize=(S.FIG_W_DOUBLE, 7.6), layout="constrained")
    gs = fig.add_gridspec(2, 1 + len(ZOOMS), height_ratios=[1, 0.035], width_ratios=[1.15] + [1] * len(ZOOMS))
    ax = fig.add_subplot(gs[0, 0])
    draw_island(ax, dom, land, road, bounds, {}, "a", "", reach_labels=False, pad_w=1500.0)
    paint(ax)
    ax.set_title("", loc="center")
    ax.set_title("(a) Hatteras Island", loc="left", fontsize=9)
    if hatch_villages:
        ax.legend(handles=[mpatches.Patch(facecolor="none", edgecolor="0.35", hatch="////", lw=0.0,
                                          label="village\n(developed)")],
                  loc="center left", bbox_to_anchor=(0.0, 0.62), fontsize=7, frameon=True, framealpha=0.9,
                  edgecolor="none", handlelength=1.6, handleheight=1.6)
    for g0, g1 in ZOOMS:                        # the close-up reaches outlined on the island
        sub = dom[dom.domain_id.between(g0, g1)].total_bounds
        ax.add_patch(plt.Rectangle((sub[0] - 300, sub[1] - 200), sub[2] - sub[0] + 600, sub[3] - sub[1] + 400,
                                   fill=False, ec="0.25", lw=0.8, ls="--", zorder=6))
    x0, x1 = ax.get_xlim()
    ax.set_xlim(x0 + 0.44 * (x1 - x0), x1 + 0.44 * (x1 - x0))
    scalebar_and_north(ax)
    ax.set_xlim(x0, x1)
    for k, (g0, g1) in enumerate(ZOOMS):         # south to north, as the GIS numbers run
        a = fig.add_subplot(gs[0, 1 + k])
        draw_island(a, dom, land, road, bounds, {}, "a", "", reach_labels=False, pad_w=0.0, label_every=2)
        paint(a, range(max(1, g0 - 1), min(90, g1 + 1) + 1), lw=0.7)
        sub = dom[dom.domain_id.between(g0, g1)].total_bounds
        xs = []                                  # cropped to the data, room for the domain numbers east of it
        for g in range(g0, g1 + 1):
            (l, r, bt, tp), dd = diff[g][1], crop[g]
            cols = np.where(np.isfinite(dd).any(axis=0))[0]
            if len(cols):
                w = (r - l) / dd.shape[1]
                xs += [l + cols[0] * w, l + (cols[-1] + 1) * w]
        a.set_xlim(min(xs) - 150, max(xs) + 900)
        a.set_ylim(sub[1] - 200, sub[3] + 200)
        a.set_title("", loc="center")
        a.set_title(f"({chr(98 + k)}) {ZOOM_NAMES.get((g0, g1), '')}\nGIS {g0}–{g1}", loc="left", fontsize=8.5)
        for t in a.texts:
            t.set_clip_on(True)
    cb = fig.colorbar(ScalarMappable(norm=norm_l, cmap=cmap_l), cax=fig.add_subplot(gs[1, :]),
                      orientation="horizontal", extend="both")
    cb.set_label(cbar_label, fontsize=8)
    cb.ax.tick_params(labelsize=7.5)
    fig.suptitle(title, fontsize=10)
    S.save(fig, o, vector=False, close=True)
    S.record_caption(o, caption + (" Hatched: the village domains (" + ", ".join(
        f"GIS {a_}-{b_}" for a_, b_ in villages) + "; overwash_data.SECTIONS), where buildings, pavement and "
        "cleared sand make the change hard to read." if hatch_villages else ""))
    print(f"wrote {o}")


def observed():
    """The lidar change alone from the calibration period's first survey to its last (1996 -> 2009)."""
    y0, y1 = CAL_SURVEYS[0], CAL_SURVEYS[-1]
    diff, off = lidar_change(y0, y1)
    print(f"  {y0} -> {y1}: datum offset {off:+.3f} m", flush=True)
    summary_map(
        diff, f"Lidar ground change over the calibration period, {y0} → {y1}",
        f"ground change {y0} → {y1}, m (blue gain, red loss; clear < {SHOW_M:g} m or no data)",
        CD.ALL_DIR / f"lidar_change_map_observed_{y0}_{y1}.png",
        f"The lidar change in ground height from the {y0} survey (the starting dunes' survey, on the 2009 datum as "
        f"the survey crests are) to the {y1} survey, at 10 m on the GIS domains' clip grids: (a) the island, dashed "
        "boxes the close-ups; (b-d) Cape Point to Buxton, north of Avon and Rodanthe to the S-curves, south to "
        f"north. The island-wide median change ({off:+.2f} m) is taken off. Blue gain, red loss; clear where the "
        f"change is under {SHOW_M:g} m, either survey is below {WATER_NAVD_M:g} m NAVD88, or either survey is "
        "missing. The change holds washover, dune and shoreline change, vegetation, road work (NC-12 at Buxton "
        "moved in 1999) and building.")


def modelled(run=DEFAULT_RUN):
    """The same map for a run: its ground change between the model states of the 1996 and 2009 surveys, shown only
    where the lidar has both surveys, so the two maps cover the same ground."""
    from HAT_collision_design_cs_sweep import cs_label
    y0, y1 = CAL_SURVEYS[0], CAL_SURVEYS[-1]
    lid, _ = lidar_change(y0, y1)
    states = model_states(run, y0, y1)
    diff = {g: (np.where(np.isfinite(lid[g][0]), states[g][1] - states[g][0], np.nan), lid[g][1])
            for g in range(1, 91)}
    label = next((x.name for x in CD.DESIGNS if x.run == run), run)
    for (key, c), v in H.DESIGN_SWEEP.items():
        if v == run:
            e, s = key.split("_")
            label = f"{CD.EROSION[e][0]}, {CD.START[s][0].lower()}, Cs {cs_label(c)}"
    summary_map(
        diff, f"Model ground change over the calibration period, {y0} → {y1}\n{label}",
        f"model ground change {y0} → {y1}, m (blue gain, red loss; clear < {SHOW_M:g} m or no lidar)",
        CD.ALL_DIR / f"model_change_map_{run}_{y0}_{y1}.png",
        f"The model's change in ground height (run {run}) from its state for the {y0} survey (its start) to its "
        f"state for the {y1} survey, on the GIS domains' 10 m clip grids, drawn as lidar_change_map_observed_"
        f"{y0}_{y1}.png and only where that map has both surveys. The model's grid goes back on the map through "
        "the start's profile frame: column i on profile i, the two dune rows at and seaward of the 1996 crest, "
        "the interior behind, the whole column moved landward with the run's shoreline; the cells left behind "
        f"are beach at the berm. Blue gain, red loss; clear under {SHOW_M:g} m.",
        crop={g: lid[g][0] for g in lid}, hatch_villages=False)


if __name__ == "__main__":
    a = sys.argv[1:]
    if a and a[0] == "observed":
        observed()
    elif a and a[0] == "model":
        modelled(a[1] if len(a) > 1 else DEFAULT_RUN)
    else:
        main(a[0] if a else DEFAULT_RUN)
