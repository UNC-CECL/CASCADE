"""
How built up each village is, from building footprints: the evidence for a per-village overwash filter.

    python scripts/real-world-cases/Hatteras/input_prep/4-mgmt-forcings/village_building_density/HAT_village_building_density.py

The beach/dune manager returns a share of each year's overwash in a village to the shoreface (overwash_filter;
one 40% for every village now, the residential end of Rogers et al.'s 2015 40-90%, residential to commercial).
Whether the villages should differ is a question of how much each blocks overwash, so this measures it rather
than fitting it to the shoreline (2026-11-03, Hannah). Per GIS domain, on its 500 m box: the land (2009 lidar
above MHW, 0.36 m NAVD88), the shoreline per 1 m row (the seaward-most land pixel), and the building
footprints rasterized at 1 m. Building cover (footprint share of the land) in the oceanfront band (the first
OCEANFRONT_M landward of the shoreline, where village overwash lands) and on all the land; building count; and
the count whose footprint reaches into the first FRONT_ROW_M (the oceanfront row). Per village (VILLAGES, the
community zones split by village) the domains pooled. A measurement only: the runner still reads one filter.

Footprints: Microsoft US Building Footprints v2 (ODbL; machine-learned from Bing imagery, mostly 2014-2019),
gis_sources/buildings/ms_usbuildings_v2/, clipped here to the island once (hatteras_buildings.geojson). Not the
Dare County tax parcels: they carry owners' names.

Writes data/hatteras_init/4-mgmt-forcing/village_building_density/building_density_by_domain.csv and
building_density_by_village.csv, and output/figures/3-model-inputs/4-management/village_building_density/
village_building_density.png (PDF and caption under supporting/).

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-11-03
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

_HERE = Path(__file__).resolve()
_HAT_DIR = next(_p for _p in _HERE.parents if (_p / "hatteras_paths.py").exists())
if str(_HAT_DIR) not in sys.path:
    sys.path.insert(0, str(_HAT_DIR))
from hatteras_paths import WORKSPACE as PROJECT_ROOT  # noqa: E402
from site_layer.hat_gis_root import GIS_ROOT  # noqa: E402
import site_layer.hat_figure_style as fs  # noqa: E402

# --- CONFIG ------------------------------------------------------------------
FOOTPRINTS = GIS_ROOT / "buildings" / "ms_usbuildings_v2"
FOOTPRINT_ZIP = FOOTPRINTS / "NorthCarolina.geojson.zip"
FOOTPRINT_CLIP = FOOTPRINTS / "hatteras_buildings.geojson"
BBOX_LONLAT = (-75.75, 35.18, -75.40, 35.82)       # Hatteras Island with a margin
DOMAINS = PROJECT_ROOT / "data" / "hatteras_init" / "5-scr" / "2-transect-frame" / "transect_domains" / "HAT_domains.json"
DEM_2009 = GIS_ROOT / "elevation" / "dem" / "2009" / "noaa_usace2009_nc_dem_Job1076020" / "2009_full.tif"
MHW_NAVD88_M = 0.36
OCEANFRONT_M = 250        # the band landward of the shoreline where village overwash lands
FRONT_ROW_M = 100         # a building reaching this close to the shoreline is in the oceanfront row
VILLAGES = {"Buxton": (7, 8), "Avon": (21, 31), "Salvo": (68, 71), "Waves": (72, 76), "Rodanthe": (77, 83)}
OUT_DIR = PROJECT_ROOT / "data" / "hatteras_init" / "4-mgmt-forcing" / "village_building_density"
FIG_DIR = PROJECT_ROOT / "output" / "figures" / "3-model-inputs" / "4-management" / "village_building_density"
# Roads (2026-11-03): where overwash is cleared, NC-12 by NCDOT and the village streets by the county and owners, and
# much of it put back on the beach - the model's village rule (the share returned to the shoreface), where the
# buildings are what blocks and traps it. NCDOT's route inventory (gis_sources/roads/NCRoutes.gdb) by RouteClass,
# 3 NC routes (NC-12), 4 state secondary roads, 5-6 non-system; it has no widths, so the pavement takes typical ones:
# NC-12 two 12 ft lanes (7.3 m), secondary roads 20 ft (6.1 m), the rest 16 ft (5.0 m). Road length needs none
ROADS = GIS_ROOT / "roads" / "NCRoutes.gdb"
ROAD_WIDTH_M = {3: 7.3, 4: 6.1, 5: 5.0, 6: 5.0}
FRONT_ROW_CUTOFFS_M = (75, 100, 150)      # the front-row count at each, as the cutoff matters (Avon's row is ~110 m)
# -----------------------------------------------------------------------------


def footprints():
    import geopandas as gpd
    if not FOOTPRINT_CLIP.exists():
        b = gpd.read_file(f"/vsizip/{FOOTPRINT_ZIP.as_posix()}/NorthCarolina.geojson", bbox=BBOX_LONLAT)
        b.to_file(FOOTPRINT_CLIP, driver="GeoJSON")
        print(f"clipped {len(b)} buildings to {FOOTPRINT_CLIP}", flush=True)
    return gpd.read_file(FOOTPRINT_CLIP)


def roads():
    import geopandas as gpd
    import pyogrio
    from pyproj import Transformer
    crs = pyogrio.read_info(ROADS)["crs"]
    tr = Transformer.from_crs(4326, crs, always_xy=True)
    (x0, y0), (x1, y1) = tr.transform(*BBOX_LONLAT[:2]), tr.transform(*BBOX_LONLAT[2:])
    r = gpd.read_file(ROADS, bbox=(x0, y0, x1, y1))
    r["geometry"] = r.geometry.force_2d()
    r = r.set_crs(2264, allow_override=True).to_crs(26918)    # the NC State Plane (ftUS) part of its compound CRS
    r["route_class"] = r.RouteClass.astype(int)
    r["width_m"] = r.route_class.map(ROAD_WIDTH_M).fillna(5.0)
    r["ncdot"] = r.RouteMaintCode.str.strip().eq("System")     # NCDOT-maintained; "Non-System" otherwise
    return r


def domain_boxes():
    d = json.loads(DOMAINS.read_text(encoding="utf-8"))
    out = {}
    for f in d["features"]:
        xs, ys = zip(*f["geometry"]["coordinates"][0])
        out[int(f["properties"]["domain_id"])] = (min(xs), min(ys), max(xs), max(ys))
    return out


def measure():
    import rasterio
    from rasterio.crs import CRS
    from rasterio.features import rasterize
    from rasterio.transform import from_origin
    from rasterio.warp import reproject, Resampling
    from shapely.geometry import box
    utm = CRS.from_epsg(26918)     # the boxes' EPSG:3725 is the same UTM 18N grid to within centimetres
    b = footprints().to_crs(utm)
    rd = roads()
    rows = []
    with rasterio.open(DEM_2009) as dem:
        for g, (x0, y0, x1, y1) in sorted(domain_boxes().items()):
            w, h = int(round(x1 - x0)), int(round(y1 - y0))
            t = from_origin(x0, y1, 1.0, 1.0)
            z = np.full((h, w), np.nan, np.float32)
            reproject(rasterio.band(dem, 1), z, src_crs=dem.crs or utm, src_nodata=dem.nodata, dst_transform=t,
                      dst_crs=utm, dst_nodata=np.nan, resampling=Resampling.bilinear)
            land = np.isfinite(z) & (z > MHW_NAVD88_M) & (z < 50)
            # each row's shoreline: the seaward-most (east) land pixel; distance landward of it per pixel
            has = land.any(axis=1)
            shore = np.where(has, w - 1 - np.argmax(land[:, ::-1], axis=1), -1)
            dist = (shore[:, None] - np.arange(w)[None, :]).astype(float)
            dist[~has] = np.nan
            band = land & (dist >= 0) & (dist <= OCEANFRONT_M)
            sel = b[b.intersects(box(x0, y0, x1, y1))]
            built = (rasterize([(geom, 1) for geom in sel.geometry], out_shape=(h, w), transform=t, fill=0)
                     .astype(bool) if len(sel) else np.zeros((h, w), bool))
            # buildings reaching the oceanfront row: any of their pixels within each cutoff of the shoreline
            front = {m: 0 for m in FRONT_ROW_CUTOFFS_M}
            if len(sel):
                ids = rasterize([(geom, i + 1) for i, geom in enumerate(sel.geometry)], out_shape=(h, w),
                                transform=t, fill=0)
                for m in FRONT_ROW_CUTOFFS_M:
                    near = (dist >= 0) & (dist <= m)
                    front[m] = len(np.unique(ids[near & (ids > 0)]))
            # roads: the pavement (each centreline buffered to its typical width), the length in the band, and how
            # far NC-12 runs from the shoreline (the median over the rows it crosses)
            rsel = rd[rd.intersects(box(x0, y0, x1, y1))]
            def pave(sub):
                return (rasterize([(geom.buffer(wd / 2, cap_style="flat"), 1) for geom, wd in zip(sub.geometry, sub.width_m)],
                                  out_shape=(h, w), transform=t, fill=0).astype(bool)
                        if len(sub) else np.zeros((h, w), bool))

            def centre(sub):
                return (rasterize([(geom, 1) for geom in sub.geometry], out_shape=(h, w), transform=t, fill=0,
                                  all_touched=True).astype(bool) if len(sub) else np.zeros((h, w), bool))

            # who clears it: NCDOT's own roads (RouteMaintCode "System": NC-12 and the state secondary roads) and the
            # rest (private, HOA or county)
            state = rsel.ncdot
            paved_ncdot, paved_other = pave(rsel[state]), pave(rsel[~state])
            paved = paved_ncdot | paved_other
            line_ncdot, line_other = centre(rsel[state]), centre(rsel[~state])
            line = line_ncdot | line_other
            nc12 = rsel[rsel.route_class == 3]
            nc12_px = (rasterize([(geom, 1) for geom in nc12.geometry], out_shape=(h, w), transform=t, fill=0,
                                 all_touched=True).astype(bool) if len(nc12) else np.zeros((h, w), bool))
            d12 = np.where(nc12_px, dist, np.nan)
            d12_rows = np.nanmin(np.where(np.isfinite(d12), d12, np.inf), axis=1)
            d12_rows = d12_rows[np.isfinite(d12_rows)]
            hard = built | paved
            rows.append(dict(gis=g, buildings=len(sel), buildings_front_row=front[FRONT_ROW_M],
                             **{f"buildings_within_{m}m": front[m] for m in FRONT_ROW_CUTOFFS_M},
                             land_m2=int(land.sum()), band_land_m2=int(band.sum()),
                             cover_land=float(built[land].mean()) if land.any() else np.nan,
                             cover_oceanfront=float(built[band].mean()) if band.any() else np.nan,
                             paved_land=float(paved[land].mean()) if land.any() else np.nan,
                             paved_oceanfront=float(paved[band].mean()) if band.any() else np.nan,
                             paved_ncdot_oceanfront=float((paved_ncdot & ~paved_other)[band].mean()) if band.any() else np.nan,
                             paved_other_oceanfront=float(paved_other[band].mean()) if band.any() else np.nan,
                             hard_oceanfront=float(hard[band].mean()) if band.any() else np.nan,
                             # 1 m pixels along the centrelines: about their length in m (a little over on diagonals)
                             road_m_band=float((line & band).sum()),
                             road_ncdot_m_band=float((line_ncdot & band).sum()),
                             road_other_m_band=float((line_other & band).sum()),
                             nc12_from_shore_m=float(np.median(d12_rows)) if d12_rows.size else np.nan))
            if g % 10 == 0:
                print(f"GIS {g} done", flush=True)
    return pd.DataFrame(rows)


def village_of(g):
    return next((n for n, (lo, hi) in VILLAGES.items() if lo <= g <= hi), "")


def summarise(t):
    t = t.assign(village=t.gis.map(village_of))
    out = []
    for name, d in [*((n, t[t.village == n]) for n in VILLAGES), ("outside the villages", t[t.village == ""])]:
        out.append(dict(village=name, gis=f"{d.gis.min()}-{d.gis.max()}" if name in VILLAGES else "",
                        domains=len(d), buildings=int(d.buildings.sum()),
                        buildings_per_km_shore=d.buildings.sum() / (0.5 * len(d)),
                        front_row_per_km=d.buildings_front_row.sum() / (0.5 * len(d)),
                        cover_oceanfront=(d.cover_oceanfront * d.band_land_m2).sum() / d.band_land_m2.sum(),
                        cover_land=(d.cover_land * d.land_m2).sum() / d.land_m2.sum(),
                        **{f"within_{m}m_per_km": d[f"buildings_within_{m}m"].sum() / (0.5 * len(d))
                           for m in FRONT_ROW_CUTOFFS_M},
                        paved_oceanfront=(d.paved_oceanfront * d.band_land_m2).sum() / d.band_land_m2.sum(),
                        hard_oceanfront=(d.hard_oceanfront * d.band_land_m2).sum() / d.band_land_m2.sum(),
                        road_km_per_km_band=d.road_m_band.sum() / (0.5 * len(d)) / 1000,
                        paved_ncdot_oceanfront=(d.paved_ncdot_oceanfront * d.band_land_m2).sum() / d.band_land_m2.sum(),
                        paved_other_oceanfront=(d.paved_other_oceanfront * d.band_land_m2).sum() / d.band_land_m2.sum(),
                        road_ncdot_km_per_km_band=d.road_ncdot_m_band.sum() / (0.5 * len(d)) / 1000,
                        road_other_km_per_km_band=d.road_other_m_band.sum() / (0.5 * len(d)) / 1000,
                        nc12_from_shore_m=d.nc12_from_shore_m.median()))
    return t, pd.DataFrame(out)


# The plan-view windows: two domains (1 km of shore) near each village's middle, all at one scale
MAP_DOMAINS = {"Buxton": (7, 8), "Avon": (24, 25), "Salvo": (69, 70), "Waves": (73, 74), "Rodanthe": (78, 79)}
MAP_LANDWARD_M, MAP_SEAWARD_M = 700, 80
LAND_C, BUILDING_C = "#e9e5dc", "0.12"
NCDOT_ROAD_C, OTHER_ROAD_C = "0.38", "0.70"


# One village window: land and water from the 2009 lidar on a 2 m grid, and each row's shoreline (easting by
# northing, a running median over 50 m)
def window(name, boxes):
    import rasterio
    from rasterio.crs import CRS
    from rasterio.transform import from_origin
    from rasterio.warp import reproject, Resampling
    lo, hi = MAP_DOMAINS[name]
    y0, y1 = boxes[lo][1], boxes[hi][3]
    xa, xb = min(boxes[lo][0], boxes[hi][0]), max(boxes[lo][2], boxes[hi][2])
    res = 2.0
    w, h = int((xb - xa) / res), int((y1 - y0) / res)
    utm = CRS.from_epsg(26918)
    z = np.full((h, w), np.nan, np.float32)
    with rasterio.open(DEM_2009) as dem:
        reproject(rasterio.band(dem, 1), z, src_crs=dem.crs or utm, src_nodata=dem.nodata,
                  dst_transform=from_origin(xa, y1, res, res), dst_crs=utm, dst_nodata=np.nan,
                  resampling=Resampling.bilinear)
    land = np.isfinite(z) & (z > MHW_NAVD88_M) & (z < 50)
    xs = xa + (np.arange(w) + 0.5) * res
    ys = y1 - (np.arange(h) + 0.5) * res
    shore = np.where(land.any(axis=1), xs[w - 1 - np.argmax(land[:, ::-1], axis=1)], np.nan)
    shore = pd.Series(shore).rolling(25, center=True, min_periods=5).median().to_numpy()
    xm = np.nanmedian(shore)
    return land, (xa, xb, y0, y1), ys, shore, (xm - MAP_LANDWARD_M, xm + MAP_SEAWARD_M, y0, y1)


def _panel_title(ax, i, text):
    ax.set_title(f"({chr(ord('a') + i)})  {text}", loc="left", fontsize=9, fontweight="bold", pad=4)


def figure(t, v):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import ListedColormap
    from matplotlib.lines import Line2D
    from matplotlib.patches import Patch
    from rasterio.crs import CRS
    from shapely.geometry import box
    fs.apply_style()
    b = footprints().to_crs(CRS.from_epsg(26918))
    rd = roads()
    boxes = domain_boxes()
    vv = v[v.village != "outside the villages"].set_index("village")

    fig = plt.figure(figsize=fs.figsize("double", height=7.4), layout="constrained")
    gs = fig.add_gridspec(3, 1, height_ratios=[2.05, 2.0, 1.85])
    top = gs[0].subgridspec(1, len(MAP_DOMAINS))
    land_cmap = ListedColormap([fs.C["WATER"], LAND_C])
    for i, name in enumerate(MAP_DOMAINS):
        ax = fig.add_subplot(top[0, i])
        land, (xa, xb, y0, y1), ys, shore, ext = window(name, boxes)
        ax.imshow(land.astype(int), extent=(xa, xb, y0, y1), cmap=land_cmap, vmin=0, vmax=1,
                  interpolation="nearest", zorder=0)
        win = box(ext[0], ext[2], ext[1], ext[3])
        rw = rd[rd.intersects(win)]
        rw[~rw.ncdot].plot(ax=ax, color=OTHER_ROAD_C, lw=0.6, zorder=1.4)
        rw[rw.ncdot & (rw.route_class != 3)].plot(ax=ax, color=NCDOT_ROAD_C, lw=1.0, zorder=1.5)
        rw[rw.route_class == 3].plot(ax=ax, color=fs.C["ROAD"], lw=1.8, zorder=1.6)
        b[b.intersects(win)].plot(ax=ax, color=BUILDING_C, lw=0, zorder=2)
        for off, ls in ((FRONT_ROW_M, (0, (3, 2))), (OCEANFRONT_M, (0, (1, 1.5)))):
            ax.plot(shore - off, ys, color=fs.C["ACCENT"], lw=0.9, ls=ls, zorder=3)
        ax.set_xlim(ext[0], ext[1])
        ax.set_ylim(ext[2], ext[3])
        ax.set_aspect("equal")
        ax.set_xticks([])
        ax.set_yticks([])
        for sp in ax.spines.values():
            sp.set_linewidth(0.6)
        lo, hi = MAP_DOMAINS[name]
        _panel_title(ax, i, name)
        ax.text(0.04, 0.985, f"GIS {lo}–{hi}", transform=ax.transAxes, ha="left", va="top", fontsize=7,
                color=fs.INK_MUTED, bbox=dict(facecolor="white", alpha=0.85, edgecolor="none", pad=1.2))
        r = vv.loc[name]
        ax.text(0.04, 0.02, f"{100 * r.hard_oceanfront:.1f}% built or paved\nNC-12 {r.nc12_from_shore_m:.0f} m back",
                transform=ax.transAxes, ha="left", va="bottom", fontsize=7, color=fs.INK,
                bbox=dict(facecolor="white", alpha=0.85, edgecolor="none", pad=1.5))
        if i == len(MAP_DOMAINS) - 1:
            bx, by = ext[1] - 40 - 200, ext[3] - 60
            ax.plot([bx, bx + 200], [by, by], color=fs.INK, lw=2.0, solid_capstyle="butt", zorder=5)
            ax.text(bx + 100, by + 15, "200 m", ha="center", va="bottom", fontsize=7, zorder=5,
                    bbox=dict(facecolor="white", alpha=0.85, edgecolor="none", pad=0.8))
            ax.annotate("N", xy=(0.88, 0.80), xytext=(0.88, 0.70), xycoords="axes fraction",
                        ha="center", va="center", fontsize=7, fontweight="bold",
                        arrowprops=dict(arrowstyle="-|>", lw=0.8, color=fs.INK))
    handles = [Patch(facecolor=BUILDING_C, label="Building footprint"),
               Line2D([], [], color=fs.C["ROAD"], lw=1.8, label="NC-12"),
               Line2D([], [], color=NCDOT_ROAD_C, lw=1.0, label="NCDOT secondary roads"),
               Line2D([], [], color=OTHER_ROAD_C, lw=0.6, label="Private, HOA or county roads"),
               Patch(facecolor=LAND_C, edgecolor="0.6", lw=0.4, label="Land (2009 lidar, above MHW)"),
               Patch(facecolor=fs.C["WATER"], label="Water"),
               Line2D([], [], color=fs.C["ACCENT"], lw=0.9, ls=(0, (3, 2)), label=f"{FRONT_ROW_M} m from the shoreline"),
               Line2D([], [], color=fs.C["ACCENT"], lw=0.9, ls=(0, (1, 1.5)),
                      label=f"{OCEANFRONT_M} m from the shoreline")]
    fig.legend(handles=handles, loc="outside upper center", ncols=4, fontsize=7.5, frameon=False)

    # Along the island
    ax = fig.add_subplot(gs[1])
    ax.bar(t.gis, 100 * t.cover_oceanfront, width=0.8, color=fs.C["ACCENT"], lw=0, zorder=2,
           label=f"Within {OCEANFRONT_M} m of the shoreline")
    ax.plot(t.gis, 100 * t.cover_land, color=fs.INK, lw=1.0, zorder=3, label="All land in the domain")
    fs.town_bands(ax, spans=VILLAGES, fontsize=7, strip=0.08, shade="0.93")
    ax.set_xlim(0.5, 90.5)
    ax.set_ylim(0, 15)
    ax.set_yticks(range(0, 13, 2))
    ax.yaxis.grid(True, zorder=0)
    ax.set_xlabel(fs.DOMAIN_AXIS_LABEL)
    ax.set_ylabel("Building cover (% of land)")
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, 0.8), ncols=2, fontsize=7.5, frameon=False)
    _panel_title(ax, len(MAP_DOMAINS), "Building cover along the island")

    # Per village: what blocks the overwash, the front row, and how close NC-12 runs to the shore
    bot = gs[2].subgridspec(1, 3)
    names = list(vv.index)
    x = np.arange(len(names))
    ax = fig.add_subplot(bot[0, 0])
    bld = 100 * vv.cover_oceanfront.to_numpy()
    pav_dot, pav_oth = 100 * vv.paved_ncdot_oceanfront.to_numpy(), 100 * vv.paved_other_oceanfront.to_numpy()
    pav = pav_dot + pav_oth
    ax.bar(x, bld, 0.62, color=fs.C["ACCENT"], lw=0, zorder=2, label="Buildings")
    ax.bar(x, pav_dot, 0.62, bottom=bld, color=NCDOT_ROAD_C, lw=0, zorder=2, label="NCDOT roads")
    ax.bar(x, pav_oth, 0.62, bottom=bld + pav_dot, color=OTHER_ROAD_C, lw=0, zorder=2, label="Other roads")
    for xi, tot in zip(x, bld + pav):
        ax.text(xi, tot, f"{tot:.1f}", ha="center", va="bottom", fontsize=7, color=fs.INK)
    ax.set_xticks(x, names, rotation=30, ha="right", rotation_mode="anchor")
    ax.set_ylim(0, (bld + pav).max() * 1.25)
    ax.yaxis.grid(True, zorder=0)
    ax.set_ylabel(f"Cover within {OCEANFRONT_M} m (% of land)", fontsize=8)
    ax.legend(loc="upper right", fontsize=7, frameon=False)
    _panel_title(ax, len(MAP_DOMAINS) + 1, "Buildings and roads")
    ax = fig.add_subplot(bot[0, 1])
    shades = (fs.C["BASE_FILL"], fs.C["BASE"], fs.INK)
    wbar = 0.26
    for k, (m, c) in enumerate(zip(FRONT_ROW_CUTOFFS_M, shades)):
        ax.bar(x + (k - 1) * wbar, vv[f"within_{m}m_per_km"], wbar, color=c, lw=0, zorder=2, label=f"{m} m")
    ax.set_xticks(x, names, rotation=30, ha="right", rotation_mode="anchor")
    ax.yaxis.grid(True, zorder=0)
    ax.set_ylabel("Buildings per km of shore", fontsize=8)
    ax.legend(title="Within", title_fontsize=7, loc="upper right", fontsize=7, frameon=False)
    _panel_title(ax, len(MAP_DOMAINS) + 2, "Oceanfront-row buildings")
    ax = fig.add_subplot(bot[0, 2])
    d12 = vv.nc12_from_shore_m.to_numpy()
    bars = ax.bar(x, d12, 0.62, color=fs.C["ROAD"], lw=0, zorder=2)
    for rect, val in zip(bars, d12):
        ax.text(rect.get_x() + rect.get_width() / 2, val, f"{val:.0f}", ha="center", va="bottom", fontsize=7)
    ax.set_xticks(x, names, rotation=30, ha="right", rotation_mode="anchor")
    ax.set_ylim(0, np.nanmax(d12) * 1.18)
    ax.yaxis.grid(True, zorder=0)
    ax.set_ylabel("NC-12 from the shoreline (m)", fontsize=8)
    _panel_title(ax, len(MAP_DOMAINS) + 3, "NC-12 setback")
    png = FIG_DIR / "village_building_density.png"
    fs.save(fig, png)
    fs.record_caption(png, (
        f"Buildings and roads in the Hatteras villages: building footprints from the Microsoft US Building Footprints "
        f"(v2, ODbL; mostly 2014-2019 imagery), roads from NCDOT's route inventory (NC-12 and the state and "
        f"non-system roads). (a-e) Plan view of 1 km of shore (two GIS domains) near the middle of each village, all "
        f"at one scale and north up, on the land of the 2009 lidar (above MHW, 0.36 m NAVD88), with lines "
        f"{FRONT_ROW_M} m and {OCEANFRONT_M} m landward of the shoreline (the seaward edge of the land, smoothed "
        f"over 50 m); labels: the village's built-or-paved share of the land within {OCEANFRONT_M} m (g) and how far "
        f"NC-12 runs from the shoreline (i). (f) Per GIS domain, the building footprint share of the land within "
        f"{OCEANFRONT_M} m of the shoreline, where village overwash lands (bars), and of all the domain's land "
        f"(line); villages shaded. (g) Per village, the share of the land within {OCEANFRONT_M} m under buildings, under "
        f"NCDOT-maintained pavement (NC-12 and the state secondary roads) and under private, HOA or county pavement (centrelines at typical widths: NC-12 7.3 m, state secondary roads 6.1 m, others 5.0 m; the "
        f"inventory has none). (h) Buildings reaching within 75, 100 and 150 m of the shoreline, per km of shore. "
        f"(i) NC-12's distance from the shoreline, the median over the domain rows it crosses, then over the "
        f"village's domains. Every village is low-density residential (2-15% built or paved); they differ in their "
        f"oceanfront row (dense at Buxton and Rodanthe, none at Salvo and Waves within 100 m) and in NC-12, which "
        f"runs 140 m behind the shoreline at Buxton and about 600 m behind it in the Tri-Village. The beach/dune "
        f"manager returns 40% of village overwash to the shoreface everywhere now, the residential end of Rogers "
        f"et al.'s (2015) 40-90%."))
    print(f"wrote {png}")


def main():
    if "--figure" in sys.argv:
        # redraw from the saved tables
        figure(pd.read_csv(OUT_DIR / "building_density_by_domain.csv"),
               pd.read_csv(OUT_DIR / "building_density_by_village.csv"))
        return
    t = measure()
    t, v = summarise(t)
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    t.round(4).to_csv(OUT_DIR / "building_density_by_domain.csv", index=False)
    v.round(4).to_csv(OUT_DIR / "building_density_by_village.csv", index=False)
    pd.set_option("display.width", 200)
    print(v.round(3).to_string(index=False))
    figure(t, v)


if __name__ == "__main__":
    main()
