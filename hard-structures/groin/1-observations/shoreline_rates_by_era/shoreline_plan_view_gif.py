#!/usr/bin/env python3
"""
The shoreline around the Buxton groins as a true map, one frame per year, from CoastSat.

    python shoreline_plan_view_gif.py [--only whole|closeup|tiles|static] [--basemap plain|satellite]
        plain:     figures/shoreline_map_whole_reach.gif, figures/shoreline_map_groin_closeup.gif,
                   figures/shoreline_map_tiles/shoreline_map_tile_<nn>_GIS<a>-<b>.gif,
                   figures/shoreline_rates_by_era_5_groin_window_four_years.png (--only static)
        satellite: the same GIFs over satellite imagery, under figures/v2_satellite/

Every CoastSat observation is placed at its real position (UTM 18N): the transect's origin plus
its chainage along that transect's own shore-normal direction, from the CoastSat transect layer.
The maps are drawn at 1:1, north up, with no stretching. Each frame is one calendar year: that
year's shoreline (the median position per transect, at least MIN_OBS images), the earlier years
coloured by year, the 1984-1986 shoreline dashed, the four groins, NC-12, the GIS domain bands and
the place names. On the plain basemap the island is filled in sand up to that year's shoreline;
on the satellite basemap the imagery (Esri World Imagery, one recent date) shows the land instead.
The 1.5 km views add a locator and the window's median change since 1984-1986, with the beach
fills that touched the window marked.

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-10-08
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import geopandas as gpd  # noqa: E402
import matplotlib.patheffects as pe  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from matplotlib.animation import PillowWriter  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402
from matplotlib.patches import Patch, Rectangle  # noqa: E402
from shapely.geometry import Polygon, box  # noqa: E402
from shapely.ops import unary_union  # noqa: E402

HERE = Path(__file__).resolve().parent
REPO = next(p for p in HERE.parents if (p / "pyproject.toml").exists())
sys.path.insert(0, str(REPO / "scripts"))
from site_layer import hat_observed_rates as obs  # noqa: E402
from site_layer.hat_figure_style import (  # noqa: E402
    C, C_1997, INK, INK_MUTED, MAP_TEXT, MAP_TEXT_DARK, apply_style, figsize, north_dart, open_frame,
    record_caption, save, scale_bar_km, spines_for_image)
from site_layer.hat_map_layers import ISLAND_OUTLINE  # noqa: E402
from site_layer.hat_topo_version import ROAD_LINE_CURRENT  # noqa: E402
from site_layer.hatteras_site_config import HATTERAS_ANNOTATIONS  # noqa: E402

# --- CONFIG ------------------------------------------------------------------
CHAINAGE = HERE / "output" / "groin_analysis_chainage_all.csv"
GROINS = HERE.parent / "gis_data" / "groins_hatteras.geojson"
DOMAINS = REPO / "data" / "hatteras_init" / "5-scr" / "2-transect-frame" / "transect_domains" / "HAT_domains.json"
FILLS = (REPO / "data" / "hatteras_init" / "4-mgmt-forcing" / "nourishment" / "2-model-input"
         / "nourishment_projects.csv")
FIGS = HERE / "figures"
UTM = "EPSG:32618"
YEARS = (1984, 2025)
REF_YEARS = (1984, 1986)              # the dashed reference shoreline
MIN_OBS = 3                           # images per transect per year
GAP_M = 300.0                         # no line across a stretch without transects
REACH = (1, 65)                       # GIS domains covered by the maps
SAND, OCEAN = "#e6d5ae", C["WATER"]
ROAD_PLAIN, ROAD_SAT = "0.35", "#ffd23f"   # NC-12: grey on the plain map, road yellow on imagery
FILL_C = C["ACCENT"]                  # beach fills on the change series
# the current NC-12 line in the repo starts ~1 km north of the groins (northing 3,902,564), so the groin
# window shows no road; the captions say so
OUTLINE_GROW_M = 250.0                # the outline widened (square corners) to reach every shoreline in the
                                      # maps: they sit at most ~210 m seaward of it (GIS 1); GIS 90 is off-map
HALF_M = 750.0                        # the tight views: 1.5 km squares
LOCATOR_HALF_M = 6000.0               # the locator: 12 km around the window
TILE_STEP_M = 1400.0                  # tile centres along the 1984-1986 shoreline; 100 m overlap
STATIC_YEARS = (1985, 1998, 2010, 2025)
FPS = 3
HOLD = 6                              # repeats of the last frame
RAMPS = {"plain": matplotlib.colors.ListedColormap(plt.get_cmap("Blues")(np.linspace(0.3, 0.95, 256))),
         "satellite": matplotlib.colors.ListedColormap(plt.get_cmap("YlOrRd")(np.linspace(0.2, 0.95, 256)))}
# -----------------------------------------------------------------------------


# Unit shore-normal direction per transect, from the CoastSat transect layer
def transect_directions(ids, origins):
    t = gpd.read_file(obs.TRANSECT_LAYER)
    t = t[t["id"].isin(ids)].to_crs(UTM)
    rows = []
    for tid, g in zip(t["id"], t.geometry):
        c = np.asarray(g.coords)
        a, b = c[0], c[-1]
        o = origins.loc[tid].to_numpy()
        if np.hypot(*(b - o)) < np.hypot(*(a - o)):      # the line runs toward the origin: flip it
            a, b = b, a
        d = (b - a) / np.hypot(*(b - a))
        rows.append((tid, d[0], d[1], np.hypot(*(a - o))))
    d = pd.DataFrame(rows, columns=["transect_id", "dir_x", "dir_y", "origin_offset_m"]).set_index("transect_id")
    assert d["origin_offset_m"].max() < 5.0, "transect lines do not start at the chainage origins"
    return d


# One year's shoreline as map coordinates in alongshore order, broken where transects are missing
def shoreline(df, order):
    df = df.assign(o=df["transect_id"].map(order)).sort_values("o")
    x, y = df["px"].to_numpy(), df["py"].to_numpy()
    br = np.where(np.hypot(np.diff(x), np.diff(y)) > GAP_M)[0] + 1
    return np.insert(x, br, np.nan), np.insert(y, br, np.nan)


# The island in one year: the widened outline minus everything seaward of that year's shoreline
def land(df, island):
    pts = df[["px", "py"]].to_numpy()
    sea = pts + 1500.0 * df[["dir_x", "dir_y"]].to_numpy()
    return island.difference(Polygon(np.vstack([pts, sea[::-1]])).buffer(0))


def fill(ax, geom, colour, z=1):
    for poly in getattr(geom, "geoms", [geom]):
        if not poly.is_empty and poly.geom_type == "Polygon":
            ax.fill(*poly.exterior.xy, color=colour, lw=0, zorder=z)
            for hole in poly.interiors:
                ax.fill(*hole.xy, color=OCEAN, lw=0, zorder=z)


class Data:
    """Everything the maps draw, computed once."""

    def __init__(self):
        ch = pd.read_csv(CHAINAGE, usecols=["transect_id", "date", "chainage_m", "source", "domain",
                                            "alongshore_m", "origin_x", "origin_y"])
        ch = ch[ch["source"] == "coastsat"]
        ch["year"] = pd.to_datetime(ch["date"]).dt.year
        origins = ch.groupby("transect_id")[["origin_x", "origin_y"]].first()
        self.dom = ch.groupby("transect_id")["domain"].first()
        self.order = ch.groupby("transect_id")["alongshore_m"].first()
        dirs = transect_directions(origins.index, origins)
        print(f"{len(dirs)} transects with directions")
        ann = (ch[ch["year"].between(*YEARS)].groupby(["transect_id", "year"])["chainage_m"]
               .agg(["median", "size"]).reset_index())
        ann = ann[ann["size"] >= MIN_OBS].join(origins, on="transect_id").join(dirs, on="transect_id")
        ann = ann.dropna(subset=["dir_x"])
        ann["px"] = ann["origin_x"] + ann["median"] * ann["dir_x"]
        ann["py"] = ann["origin_y"] + ann["median"] * ann["dir_y"]
        ref = (ch[ch["year"].between(*REF_YEARS)].groupby("transect_id")["chainage_m"].median()
               .rename("median").to_frame().join(origins).join(dirs).dropna().reset_index())
        ref["px"] = ref["origin_x"] + ref["median"] * ref["dir_x"]
        ref["py"] = ref["origin_y"] + ref["median"] * ref["dir_y"]
        # change since 1984-1986 along each transect, seaward positive
        ann = ann.join(ref.set_index("transect_id")["median"].rename("ref_median"), on="transect_id")
        ann["change_m"] = ann["median"] - ann["ref_median"]
        self.ann, self.ref = ann, ref
        self.years = sorted(ann["year"].unique())
        in_reach = self.dom[self.dom.between(*REACH)].index
        self.lines = {y: shoreline(ann[(ann["year"] == y) & ann["transect_id"].isin(in_reach)], self.order)
                      for y in self.years}
        self.ref_line = shoreline(ref[ref["transect_id"].isin(in_reach)], self.order)
        self.outline = unary_union(gpd.read_file(ISLAND_OUTLINE).to_crs(UTM).geometry)
        island = self.outline.buffer(OUTLINE_GROW_M, join_style=2)
        base = ref.set_index("transect_id")[["px", "py", "dir_x", "dir_y"]]
        self.lands = {}
        for y in self.years:
            cur = ann[ann["year"] == y].set_index("transect_id")[["px", "py", "dir_x", "dir_y"]]
            full = cur.combine_first(base).reset_index()
            full = full.assign(o=full["transect_id"].map(self.order)).sort_values("o")
            self.lands[y] = land(full, island)
        g = json.load(open(GROINS))
        self.groins = [np.asarray(f["geometry"]["coordinates"]) for f in g["features"]]
        d = json.load(open(DOMAINS))
        self.domain_bands = {f["properties"]["domain_id"]: (min(c[1] for c in f["geometry"]["coordinates"][0]),
                                                            max(c[1] for c in f["geometry"]["coordinates"][0]))
                             for f in d["features"]}
        self.road = np.asarray(gpd.read_file(ROAD_LINE_CURRENT).to_crs(UTM).geometry.iloc[0].coords)
        f = pd.read_csv(FILLS)
        f = f[f["source"].str.startswith("model input") & (f["first_gis"] <= REACH[1])]
        self.fills = [dict(name=r["name"], year=int(r["year"]), lo=int(r["first_gis"]), hi=int(r["last_gis"]),
                           m3=float(r["volume_m3"])) for _, r in f.iterrows()]
        # places: each village at the middle of its domains, the piers on the 1984-1986 shoreline
        rd = ref.assign(d=ref["transect_id"].map(self.dom))
        dom_xy = rd.groupby("d")[["px", "py"]].mean()
        self.villages = {t: dom_xy.loc[a:b].mean().to_numpy()
                         for t, (a, b) in HATTERAS_ANNOTATIONS.town_spans.items() if a <= REACH[1]}
        self.piers = {p: dom_xy.loc[d].to_numpy() for p, (d, _) in HATTERAS_ANNOTATIONS.piers.items()
                      if d in dom_xy.index}


def map_furniture(ax, data, x0, x1, y0, y1, bar, style, text, segments=2):
    """Domain bands, NC-12, piers, villages, groins, scale bar and north arrow on a north-up map."""
    line_c = "white" if style == "plain" else "0.85"
    for gid, (b0, b1) in data.domain_bands.items():
        if b1 < y0 or b0 > y1:
            continue
        for b in (b0, b1):
            if y0 < b < y1:
                ax.axhline(b, color=line_c, lw=0.8, alpha=0.9 if style == "plain" else 0.6, zorder=2.5)
        if (min(b1, y1) - max(b0, y0)) > 0.12 * (y1 - y0):
            ax.text(x0 + 0.03 * (x1 - x0), (max(b0, y0) + min(b1, y1)) / 2, f"GIS {gid}", ha="left", va="center",
                    **{**text, "fontsize": 7.5, "zorder": 8})
    halo = [pe.withStroke(linewidth=2.6, foreground="white" if style == "plain" else "0.1")]
    ax.plot(data.road[:, 0], data.road[:, 1], color=ROAD_PLAIN if style == "plain" else ROAD_SAT, lw=1.2,
            path_effects=halo, zorder=3.2)
    for name, (px, py) in data.piers.items():
        if x0 < px < x1 and y0 < py < y1:
            ax.plot(px, py, marker="s", ms=4, color=INK if style == "plain" else "white", zorder=6)
            ax.text(px + 0.03 * (x1 - x0), py, name, ha="left", va="center", **{**text, "fontsize": 7, "zorder": 8})
    for name, (vx, vy) in data.villages.items():
        if x0 < vx < x1 and y0 < vy < y1:
            ax.text(vx - 0.06 * (x1 - x0), vy, name, ha="right", va="center", fontstyle="italic",
                    **{**text, "fontsize": 9, "zorder": 8})
    for gl in data.groins:
        ax.plot(gl[:, 0], gl[:, 1], color=C["GROIN"], lw=2.0, solid_capstyle="butt", zorder=5)
    ax.set_xlim(x0, x1)
    ax.set_ylim(y0, y1)
    ax.set_aspect("equal")
    ax.set_xticks([])
    ax.set_yticks([])
    spines_for_image(ax)
    scale_bar_km(ax, length_m=bar, segments=segments, unit="km" if bar >= 1000 else "m", x=0.05, y=0.05,
                 text_kw=text)
    north_dart(ax, (x0 + 0.9 * (x1 - x0), y0 + 0.1 * (y1 - y0)),   # bottom right, clear of the scale bar
               arrow_m=min(0.07 * (y1 - y0), 0.12 * (x1 - x0)), text_kw=text)


def draw_shorelines(ax, data, k, y, x0, x1, y0, y1, style, image=None):
    ramp = RAMPS[style]
    if style == "plain":
        ax.set_facecolor(OCEAN)
        fill(ax, data.lands[y].intersection(box(x0 - 50, y0 - 50, x1 + 50, y1 + 50)), SAND)
    else:
        img, ext = image
        ax.imshow(img, extent=ext, zorder=0, interpolation="bilinear")
    for yy in data.years[:k]:
        ax.plot(*data.lines[yy], color=ramp((yy - data.years[0]) / max(1, data.years[-1] - data.years[0])),
                lw=0.7, alpha=0.9, zorder=3)
    dark = style == "plain"
    ax.plot(*data.ref_line, color=INK_MUTED if dark else "white", lw=1.0, ls="--", zorder=3.5)
    ax.plot(*data.lines[y], color=INK if dark else "white", lw=1.6, zorder=4,
            path_effects=None if dark else [pe.withStroke(linewidth=2.8, foreground="0.05")])


def legend_handles(style, n_tile=True):
    dark = style == "plain"
    h = [Line2D([], [], color=INK if dark else "white", lw=1.6, label="this year's shoreline",
                path_effects=None if dark else [pe.withStroke(linewidth=2.8, foreground="0.05")]),
         Line2D([], [], color=RAMPS[style](0.6), lw=0.8, label="earlier years (by year)"),
         Line2D([], [], color=INK_MUTED if dark else "0.4", lw=1.0, ls="--",
                label=f"{REF_YEARS[0]}-{REF_YEARS[1]} shoreline"),
         Line2D([], [], color=C["GROIN"], lw=2.0, label="Buxton groins"),
         Line2D([], [], color=ROAD_PLAIN if dark else ROAD_SAT, lw=1.2, label="NC-12 (current)",
                path_effects=[pe.withStroke(linewidth=2.6, foreground="white" if dark else "0.1")])]
    if n_tile:
        h.append(Patch(facecolor=FILL_C, alpha=0.3, label="beach fill in this window"))
    if dark:
        h += [Patch(facecolor=SAND, label="island"), Patch(facecolor=OCEAN, label="ocean")]
    return h


def satellite(x0, x1, y0, y1):
    import contextily as ctx
    from pyproj import Transformer
    tr = Transformer.from_crs(UTM, "EPSG:3857", always_xy=True)
    (mx0, mx1), (my0, my1) = tr.transform([x0, x1], [y0, y1])
    img, ext = ctx.bounds2img(mx0, my0, mx1, my1, source=ctx.providers.Esri.WorldImagery)
    img, ext = ctx.warp_tiles(img, ext, t_crs=UTM)
    return img, (ext[0], ext[1], ext[2], ext[3])


def window_stats(data, x0, x1, y0, y1):
    inside = data.ref[(data.ref["px"].between(x0, x1)) & (data.ref["py"].between(y0, y1))]["transect_id"]
    sub = data.ann[data.ann["transect_id"].isin(inside)]
    series = sub.groupby("year")["change_m"].median().reindex(data.years)
    n_max = int(sub.groupby("year")["transect_id"].nunique().max()) if len(sub) else 0
    doms = data.dom.loc[inside]
    lo, hi = (int(doms.min()), int(doms.max())) if len(doms) else (0, 0)
    fills = [f for f in data.fills if f["lo"] <= hi and f["hi"] >= lo]
    return series, n_max, (lo, hi), fills


def change_panel(axt, data, series, n_max, fills, y=None, style="plain"):
    lim = np.nanmax(np.abs(series.to_numpy())) if series.notna().any() else 10
    lim = np.ceil(max(lim, 10) / 10) * 10 * 1.15
    axt.axhline(0, color=INK_MUTED, lw=0.6)
    for f in fills:
        axt.axvspan(f["year"] - 0.45, f["year"] + 0.45, color=FILL_C, alpha=0.3, lw=0, zorder=0)
        axt.text(f["year"], lim * 0.97, f"{f['name'].split()[0]} fill", rotation=90, ha="right", va="top",
                 fontsize=6, color=INK_MUTED)
    if y is None:
        axt.plot(series.index, series.to_numpy(), color=C_1997, lw=1.6, marker="o", ms=2.2)
    else:
        axt.plot(series.index, series.to_numpy(), color="0.85", lw=1.0)
        upto = series.loc[:y]
        axt.plot(upto.index, upto.to_numpy(), color=C_1997, lw=1.6, marker="o", ms=2.2)
        if np.isfinite(series.get(y, np.nan)):
            axt.plot(y, series[y], marker="o", ms=6, color=INK, zorder=5)
        axt.axvline(y, color=INK_MUTED, lw=0.6, ls=":")
    axt.set_xlim(data.years[0] - 1, data.years[-1] + 1)
    axt.set_ylim(-lim, lim)
    axt.set_ylabel(f"Change since {REF_YEARS[0]}-{REF_YEARS[1]}\n(m, seaward +)", fontsize=7)
    axt.tick_params(labelsize=6.5)
    axt.grid(axis="y")
    open_frame(axt)
    axt.set_title(f"How far has this window's shoreline moved?\nmedian of {n_max} transects", fontsize=8)


def locator(axl, data, x0, x1, y0, y1, style):
    """The window in its surroundings: about 12 km of island around it, villages, NC-12, the groins."""
    cx, cy = (x0 + x1) / 2, (y0 + y1) / 2
    lx0, lx1, ly0, ly1 = cx - LOCATOR_HALF_M, cx + LOCATOR_HALF_M, cy - LOCATOR_HALF_M, cy + LOCATOR_HALF_M
    ink = INK if style == "plain" else "white"
    if style == "plain":
        gpd.GeoSeries([data.outline]).plot(ax=axl, facecolor=SAND, edgecolor="0.55", lw=0.3)
        axl.set_facecolor(OCEAN)
    else:
        gpd.GeoSeries([data.outline]).plot(ax=axl, facecolor="0.55", edgecolor="0.3", lw=0.3)
        axl.set_facecolor("0.12")
    axl.plot(data.road[:, 0], data.road[:, 1], color=ROAD_PLAIN if style == "plain" else ROAD_SAT, lw=0.8, zorder=3)
    for gl in data.groins:
        axl.plot(gl[:, 0], gl[:, 1], color=C["GROIN"], lw=1.5, zorder=4)
    axl.add_patch(Rectangle((x0, y0), x1 - x0, y1 - y0, fill=False, edgecolor=C["LOCATOR"], lw=1.6, zorder=6))
    for name, (vx, vy) in data.villages.items():
        if lx0 < vx < lx1 and ly0 < vy < ly1:
            axl.text(vx - 900, vy, name, ha="right", va="center", fontsize=7, fontstyle="italic", color=ink, zorder=7)
    bx = data.outline.bounds
    cape = (bx[0] + 12000, bx[1] + 1500)
    if lx0 < cape[0] < lx1 and ly0 < cape[1] < ly1:
        axl.text(*cape, "Cape Point", ha="center", fontsize=7, fontstyle="italic", color=ink, zorder=7)
    axl.set_xlim(lx0, lx1)
    axl.set_ylim(ly0, ly1)
    axl.set_aspect("equal")
    axl.set_xticks([])
    axl.set_yticks([])
    spines_for_image(axl)
    scale_bar_km(axl, length_m=2000, segments=1, unit="km", x=0.06, y=0.06,
                 text_kw={**(MAP_TEXT_DARK if style == "plain" else MAP_TEXT), "fontsize": 6.5})
    axl.set_title(f"Where: {2 * LOCATOR_HALF_M / 1000:.0f} km around the window", fontsize=8)


# A 1.5 km map with a locator and the window's mean change, one frame per year
def tight_gif(data, cx, cy, out, label, style):
    x0, x1, y0, y1 = cx - HALF_M, cx + HALF_M, cy - HALF_M, cy + HALF_M
    series, n_max, (lo, hi), fills = window_stats(data, x0, x1, y0, y1)
    gis = f"GIS {lo}-{hi}" if n_max else "no transects"
    text = MAP_TEXT_DARK if style == "plain" else MAP_TEXT
    image = satellite(x0 - 60, x1 + 60, y0 - 60, y1 + 60) if style == "satellite" else None

    fig = plt.figure(figsize=figsize("double", height=5.4))
    gs = fig.add_gridspec(2, 2, width_ratios=(1.75, 1), height_ratios=(1.2, 0.95), left=0.01, right=0.985,
                          top=0.92, bottom=0.16, wspace=0.16, hspace=0.42)
    ax, axl, axt = fig.add_subplot(gs[:, 0]), fig.add_subplot(gs[0, 1]), fig.add_subplot(gs[1, 1])
    fig.legend(handles=legend_handles(style), loc="lower center", ncol=4, fontsize=6.6, frameon=False)
    locator(axl, data, x0, x1, y0, y1, style)
    writer = PillowWriter(fps=FPS)
    with writer.saving(fig, str(out), dpi=120):
        for k, y in enumerate(data.years):
            ax.clear()
            draw_shorelines(ax, data, k, y, x0, x1, y0, y1, style, image)
            map_furniture(ax, data, x0, x1, y0, y1, 500, style, text)
            for f in fills:
                if f["year"] == y:
                    ax.text(0.03, 0.97, f"{y}: {f['name']}, GIS {f['lo']}-{f['hi']}\n"
                                        f"{f['m3'] / 1e6:.1f} million m³ placed", transform=ax.transAxes,
                            ha="left", va="top", fontsize=7, zorder=9,
                            bbox=dict(facecolor="white", alpha=0.85, edgecolor=FILL_C, boxstyle="round,pad=0.3"))
            axt.clear()
            change_panel(axt, data, series, n_max, fills, y, style)
            fig.suptitle(f"The shoreline in {y}: {label}, {gis} (1.5 km across)", fontsize=9.5)
            writer.grab_frame()
        for _ in range(HOLD):
            writer.grab_frame()
    plt.close(fig)
    return gis


def whole_gif(data, out, style):
    r = data.ref[data.ref["transect_id"].isin(data.dom[data.dom.between(*REACH)].index)]
    x0, x1 = r["px"].min() - 150, r["px"].max() + 150
    y0, y1 = r["py"].min() - 150, r["py"].max() + 150
    text = MAP_TEXT_DARK if style == "plain" else MAP_TEXT
    image = satellite(x0, x1, y0, y1) if style == "satellite" else None
    fig, ax = plt.subplots(figsize=figsize("single", height=6.8), constrained_layout=True)
    fig.legend(handles=legend_handles(style, n_tile=False), loc="outside lower center", ncol=2, fontsize=6.6,
               frameon=False)
    writer = PillowWriter(fps=FPS)
    with writer.saving(fig, str(out), dpi=120):
        for k, y in enumerate(data.years):
            ax.clear()
            draw_shorelines(ax, data, k, y, x0, x1, y0, y1, style, image)
            for gl in data.groins:
                ax.plot(gl[:, 0], gl[:, 1], color=C["GROIN"], lw=2.0, zorder=5)
            halo = [pe.withStroke(linewidth=2.6, foreground="white" if style == "plain" else "0.1")]
            ax.plot(data.road[:, 0], data.road[:, 1], color=ROAD_PLAIN if style == "plain" else ROAD_SAT, lw=1.0,
                    path_effects=halo, zorder=3.2)
            for name, (vx, vy) in data.villages.items():
                if y0 < vy < y1:
                    ax.text(vx - 900, vy, name, ha="right", va="center", fontstyle="italic",
                            **{**text, "fontsize": 8, "zorder": 8})
            ax.set_xlim(x0, x1)
            ax.set_ylim(y0, y1)
            ax.set_aspect("equal")
            ax.set_xticks([])
            ax.set_yticks([])
            spines_for_image(ax)
            scale_bar_km(ax, length_m=2000, segments=2, unit="km", x=0.05, y=0.04, text_kw=text)
            north_dart(ax, (x0 + 0.85 * (x1 - x0), y0 + 0.06 * (y1 - y0)), arrow_m=0.12 * (x1 - x0), text_kw=text)
            fig.suptitle(f"The shoreline in {y}, GIS {REACH[0]}-{REACH[1]}", fontsize=9)
            writer.grab_frame()
        for _ in range(HOLD):
            writer.grab_frame()
    plt.close(fig)
    print(f"wrote {out.relative_to(REPO)}")


# Tile centres every TILE_STEP_M along the 1984-1986 shoreline, GIS 1-65
def tile_centres(data):
    r = data.ref[data.ref["transect_id"].isin(data.dom[data.dom.between(*REACH)].index)]
    r = r.assign(o=r["transect_id"].map(data.order)).sort_values("o")
    xy = r[["px", "py"]].to_numpy()
    s = np.r_[0, np.cumsum(np.hypot(*np.diff(xy, axis=0).T))]
    at = np.arange(HALF_M, s[-1] - HALF_M + TILE_STEP_M / 2, TILE_STEP_M)
    return [(np.interp(a, s, xy[:, 0]), np.interp(a, s, xy[:, 1])) for a in at]


# The groin window at four years and its change, as one static figure for papers and slides
def static_figure(data):
    cx, cy = np.vstack(data.groins).mean(axis=0)
    x0, x1, y0, y1 = cx - HALF_M, cx + HALF_M, cy - HALF_M, cy + HALF_M
    series, n_max, (lo, hi), fills = window_stats(data, x0, x1, y0, y1)
    fig = plt.figure(figsize=figsize("double", height=4.9))
    gs = fig.add_gridspec(2, 4, height_ratios=(1.0, 0.75), left=0.08, right=0.985, top=0.94, bottom=0.13,
                          wspace=0.06, hspace=0.32)
    letters = "abcde"
    for i, y in enumerate(STATIC_YEARS):
        ax = fig.add_subplot(gs[0, i])
        k = data.years.index(y)
        ax.set_facecolor(OCEAN)
        fill(ax, data.lands[y].intersection(box(x0 - 50, y0 - 50, x1 + 50, y1 + 50)), SAND)
        ax.plot(*data.ref_line, color=INK_MUTED, lw=0.9, ls="--", zorder=3.5)
        ax.plot(*data.lines[y], color=INK, lw=1.4, zorder=4)
        map_furniture(ax, data, x0, x1, y0, y1, 500, "plain", MAP_TEXT_DARK, segments=1)
        ax.set_title(f"({letters[i]}) {y}", fontsize=9, loc="left")
        del k
    axt = fig.add_subplot(gs[1, :])
    change_panel(axt, data, series, n_max, fills, None)
    axt.set_title("")
    axt.set_title(f"({letters[4]}) How far has the shoreline around the groins moved? (median of {n_max} transects)",
                  fontsize=9, loc="left")
    for y in STATIC_YEARS:
        axt.axvline(y, color=INK_MUTED, lw=0.6, ls=":")
    h = legend_handles("plain")
    fig.legend(handles=[h[0], h[2], h[3], h[5], h[6], h[7]], loc="lower center", ncol=6, fontsize=6.6,
               frameon=False)
    png = FIGS / "shoreline_rates_by_era_5_groin_window_four_years.png"
    save(fig, png, close=True)
    fl = "; ".join(f"{f['name']} {f['year']} (GIS {f['lo']}-{f['hi']}, {f['m3'] / 1e6:.1f} million m³)" for f in fills)
    s = series.dropna()
    record_caption(png, (
        "**Tests:** how far, and when, the shoreline around the Buxton groins retreated or advanced, seen as maps "
        "rather than rates. **How to read:** (a-d) the same 1.5 km window centred on the groin field "
        f"({'GIS %d-%d' % (lo, hi)}) at four years, drawn at true scale, north up: the island in sand up to that "
        "year's CoastSat shoreline (black; the median position per transect that year), the 1984-1986 shoreline "
        "dashed, the four groins in red, GIS domain bands in white (the current NC-12 line on file starts about "
        "1 km north of the groins, so no road falls in this window). (e) The median change since "
        f"1984-1986 of the {n_max} transects in the window, seaward positive; shading marks the beach fills whose "
        f"footprints reach the window ({fl}); dotted lines mark the four mapped years. **Shows:** the window's "
        f"shoreline stood {s.loc[1985]:+.0f} m from its 1984-1986 position in 1985, {s.loc[1998]:+.0f} m in 1998, "
        f"{s.loc[2010]:+.0f} m in 2010 and {s.loc[2025]:+.0f} m in 2025; its lowest was {s.min():+.0f} m in "
        f"{int(s.idxmin())}. Each position is the CoastSat transect origin plus its chainage along that transect's "
        "own shore-normal direction."))
    print(f"wrote {png.relative_to(REPO)}")


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--only", choices=("whole", "closeup", "tiles", "static"))
    ap.add_argument("--basemap", choices=("plain", "satellite"), default="plain")
    a = ap.parse_args()
    apply_style()
    data = Data()
    style = a.basemap
    root = FIGS if style == "plain" else FIGS / "v2_satellite"
    tiles = root / "shoreline_map_tiles"
    root.mkdir(parents=True, exist_ok=True)
    if a.only in (None, "static") and style == "plain":
        static_figure(data)
    if a.only in (None, "whole"):
        whole_gif(data, root / "shoreline_map_whole_reach.gif", style)
    if a.only in (None, "closeup"):
        cx, cy = np.vstack(data.groins).mean(axis=0)
        gis = tight_gif(data, cx, cy, root / "shoreline_map_groin_closeup.gif", "around the Buxton groins", style)
        print(f"wrote groin close-up ({gis})")
    if a.only in (None, "tiles"):
        tiles.mkdir(parents=True, exist_ok=True)
        for old in tiles.glob("*.gif"):
            old.unlink()
        for i, (cx, cy) in enumerate(tile_centres(data), 1):
            tmp = tiles / f"_tile_{i:02d}.gif"
            gis = tight_gif(data, cx, cy, tmp, f"tile {i}", style)
            tmp.rename(tiles / f"shoreline_map_tile_{i:02d}_{gis.replace(' ', '')}.gif")
            print(f"wrote tile {i:02d} ({gis})")


if __name__ == "__main__":
    main()
