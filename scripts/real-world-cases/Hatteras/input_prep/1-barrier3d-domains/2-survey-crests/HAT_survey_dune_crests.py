"""
Dune crests from the lidar surveys flown inside a hindcast window, on that start's dune windows.

    python scripts/input_prep/1-barrier3d-domains/2-survey-crests/HAT_survey_dune_crests.py                    # 2009-start
    python scripts/input_prep/1-barrier3d-domains/2-survey-crests/HAT_survey_dune_crests.py --start 1996-start

For each survey of the start's window, per GIS domain: put the survey on the
start product's 1 m domain clip, resample to 10 m as HAT_dem_resample_clip.py
does, orient and straighten in the start's own frame (its shear and water trim),
and take each alongshore profile's highest point inside the start's picked dune
window, clipped to the survey's own beach start, as find_dunes does: the same
crest definition the model's starting dunes have. A crest is kept only where
every cell of its window was surveyed.

Two kinds of source. "mosaic": a 1 m DEM on the grid the clips share (2014-2019),
read by pixel offset as HAT_dem_gap_fill.py does. "raster": the 1996-2005 lidar,
downloaded from NOAA's Data Access Viewer as grids of 3-9 m or ft in their own
projections, reprojected bilinearly onto the clip; their vertical unit is read
from the CRS where stated and otherwise from the ratio of their land elevations
to the start's arrays, and converted to metres. Coarse grids smooth sharp
crests, so these read crests somewhat low; the 1996 grid against the 1 m 1996
survey the starting dunes came from measures by how much.

Vertical offsets: a survey can be shifted by SURVEY_VERTICAL_OFFSET_M before
anything else (since 2026-10-20: the 1996 ALACE grid by -0.24 m, as the
1996-start DEM is, 1997 and 1999 by -0.34 and -0.32 m, and the 2014-2019 surveys
by their own -0.10 to -0.21 m: each survey's road offset from the 2009 lidar). Every survey is checked on NC-12 pavement against the 2009 lidar:
survey_sources.csv gets the offset applied and the median remaining
survey-minus-2009 difference along the 2008 road line (GIS 9-13 left out, the
road moved there in 1999), which should be near zero for a corrected survey.

Before any survey, the start's own arrays go through this code and must
reproduce its dune files. A second column widens the window by WIDEN_CELLS each
side, as a sensitivity only.

Lidar spikes (since 2026-10-26): once every survey is in, a crest more than
SPIKE_M above the same cell's crest both before (the previous survey, or the
start's own crest for the first) and after (the next survey) is blanked with
status "spike"; the value it had stays in crest_raw_m_navd88. 179 cells on the
1996 start, mostly 1996 and 1998, up to 60 m; none on the 2009 start.

Writes data/hatteras_init/1-barrier3d-domains/<start>/survey-crests/:
survey_dune_crests.csv (one row per survey, domain and 10 m cell),
survey_coverage.csv (per survey and domain, with cells_spike) and
survey_sources.csv (per survey).

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-10-26
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import geopandas as gpd  # noqa: E402
import rasterio  # noqa: E402
from rasterio.features import rasterize  # noqa: E402
import rasterio.windows  # noqa: E402
from rasterio.crs import CRS  # noqa: E402
from rasterio.enums import Resampling  # noqa: E402
from rasterio.warp import reproject, transform_bounds  # noqa: E402
from rasterio.windows import Window  # noqa: E402

# Hatteras/ holds the scripts; data/ and output/ live in the workspace on D: (hatteras_paths)
import sys as _hat_sys  # noqa: E402
_HAT_DIR = next(p for p in Path(__file__).resolve().parents  if (p / "hatteras_paths.py").exists())
if str(_HAT_DIR) not in _hat_sys.path:
    _hat_sys.path.insert(0, str(_HAT_DIR))
from hatteras_paths import WORKSPACE as REPO  # noqa: E402
sys.path.insert(0, str(_HAT_DIR))
sys.path.insert(0, str(_HAT_DIR / "input_prep" / "1-barrier3d-domains" / "1-extraction"))
import HAT_dune_topo_extractor as X  # noqa: E402
from site_layer import hat_topo_version as tv  # noqa: E402
from site_layer.hat_elevation_products import product  # noqa: E402
from site_layer.hat_gis_root import GIS_ROOT  # noqa: E402

# --- CONFIG ------------------------------------------------------------------
POLYGONS = GIS_ROOT / "elevation" / "dem"
POINTS = GIS_ROOT / "elevation" / "lidar"
# Per start: the product whose 1 m clips set the grid, and its surveys:
# year -> (kind, file or tile glob, flight dates from the metadata, the model state compared, note)
STARTS = {
    "2009-start": ("2009-2014", {
        2014: ("mosaic", POLYGONS / "2014/noaa_2014_NOAA_Post_Sandy_DEM_Job1076021/2014_full.tif",
               "2013-11 to 2014-06", "2014-01-01", ""),
        2016: ("mosaic", POLYGONS / "2016/noaa_2016_USACE_post_matthew_Job1076022/2016_full.tif",
               "2016-10 to 2016-12", "2017-01-01", ""),
        2017: ("mosaic", POLYGONS / "2017/noaa_usace2017_east_coast_dem_Job1076023/2017_full.tif",
               "2017-06 to 2017-08", "2017-01-01", ""),
        2018: ("mosaic", POLYGONS / "2018/noaa_usace2018_post_florence_dem_Job1076024/2018_full.tif",
               "2018-10", "2019-01-01", ""),
        2019: ("mosaic", POLYGONS / "2019/noaa_nc2019_dunex_Job1076025/2019_full.tif",
               "2019-09 to 2019-10", "2020-01-01", ""),
    }),
    "1996-start": ("2009-2014-1996", {
        1996: ("raster", POINTS / "1996/noaa_beta_1996_FallEC_J114106/1996_full.tif",
               "1996-08 to 1996-10", "1997-01-01", "coarse grid of the survey the starting dunes came from: bias check"),
        1997: ("raster", POINTS / "1997/noaa_1997_FallEC_Job1111523/1997_FallEC_Job1111523.tif", "1997-09 to 1997-10", "1998-01-01", ""),
        1998: ("raster", POINTS / "1998/noaa_beta_1998_FallEC_J114104/*.tif", "1998-09", "1999-01-01", ""),
        1999: ("raster", POINTS / "1999/noaa_beta_1999_SC-NC-VA_Post-Dennis_J114103/*.tif", "1999-09",
               "2000-01-01", "post-Dennis, flown before Floyd; the model state includes Floyd"),
        2000: ("raster", POINTS / "2000/noaa_beta_2000_Spring_NC_ATM_J136368/2000_full.tif", "2000-02 to 2000-03", "2000-01-01", ""),
        2001: ("raster", POINTS / "2001/noaa_2001_NC_PhaseI_Job1135397/2001_NC_PhaseI_Job1135397.tif",
               "2001-01 to 2001-03", "2001-01-01", "unclassified returns: buildings or vegetation can raise a crest"),
        2004: ("raster", POINTS / "2004/noaa_beta_2004_PreIvan_J114100/*.tif", "2004-04 to 2004-09", "2004-01-01",
               "after Isabel (2003)"),
        2005: ("raster", POINTS / "2005/noaa_beta_2005_Topobathy_J114099/*.tif", "2005-08 to 2005-11", "2006-01-01", ""),
        2009: ("mosaic", POLYGONS / "2009/noaa_usace2009_nc_dem_Job1076020/2009_full.tif", "2009-08", "2010-01-01",
               "the calibration end, here on the 1996 windows"),
    }),
}
# Added to a survey's elevations (m) before anything else: per start, per survey year. The 1996 ALACE
# grid reads NC-12 0.24 m above the 2009 lidar, as the ALACE tiles in the 1996-start DEM do
# (0-elevation/1-source-selection/HAT_alace_vertical_offset.py); the same campaign, the same shift
# 1997 and 1999 read NC-12 +0.34 and +0.32 m above 2009 on this script's road check (2026-10-20), so they
# are shifted by those amounts (Hannah); 2000 (+0.10) and 2001 (-0.10, unclassified) are left as measured
# The 2009-start surveys all read NC-12 above 2009 (2014 +0.13, 2016 +0.21, 2017 +0.10, 2018 +0.12, 2019 +0.12;
# road check 2026-10-20): consistent across three agencies, so each is shifted to the 2009 reference (Hannah,
# 2026-10-21), keeping every survey on the datum of the model's starting DEM
SURVEY_VERTICAL_OFFSET_M = {"1996-start": {1996: -0.24, 1997: -0.34, 1999: -0.32},
                            "2009-start": {2014: -0.13, 2016: -0.21, 2017: -0.10, 2018: -0.12, 2019: -0.12}}
# The road check: every survey against the 2009 lidar on the 2008 NC-12 line
REF_2009 = POLYGONS / "2009/noaa_usace2009_nc_dem_Job1076020/2009_full.tif"
ROAD_LINE = tv.MGMT_ROOT / "road_offset" / "raw_offset" / "2008" / "nc12_2008.geojson"
RELOCATED_1999 = range(9, 14)            # GIS 9-13: NC-12 moved in 1999, not the same pavement
# A raster without a CRS (the 2000 mosaic) sits on the clips' UTM 18N grid in metres
ASSUMED_CRS = CRS.from_epsg(26918)
US_FT = 1200 / 3937
WIDEN_CELLS = 2                  # sensitivity column only: the window grown this many 10 m cells each side
BLOCK = 10                       # 1 m -> 10 m
# Lidar spikes (2026-10-26): the 1996-1998 ATM grids carry returns up to 60 m (a run of 25-60 m at GIS 38 in
# 1998, 5-6 m there in 1997 and 1999), which the window's highest point picks up. A crest more than SPIKE_M above
# the cell's crest both before and after it in time is a spike (drop_spikes). Not a ceiling, which would cut the
# dunes that stand 11-13 m survey after survey (GIS 18, 26-27); not the cell's median over the surveys, which also
# caught real dunes that were later lost (GIS 81-82 in 1999, 2009-start GIS 59 and 75 after 2014)
SPIKE_M = 5.0
# -----------------------------------------------------------------------------


# The resample rule of HAT_dem_resample_clip.py: arcgis_bilinear, all four centre pixels required
def downsample(arr):
    h, w = arr.shape
    b = arr.reshape(h // BLOCK, BLOCK, w // BLOCK, BLOCK)
    lo, hi = BLOCK // 2 - 1, BLOCK // 2
    core = np.stack([b[:, lo, :, lo], b[:, lo, :, hi], b[:, hi, :, lo], b[:, hi, :, hi]])
    return np.where(np.isfinite(core).all(axis=0), core.mean(axis=0), np.nan)


def clip_grid(clip_path):
    with rasterio.open(clip_path) as c:
        return c.transform, c.shape, c.crs


def clean(arr, nodata):
    arr = arr.astype(float)
    if nodata is not None:
        arr[arr == nodata] = np.nan
    arr[(arr < -100) | (arr > 200)] = np.nan          # float max / -999999 sentinels, in m or ft
    return arr


# A 1 m mosaic on the clips' shared grid, by pixel offset; m NAVD88
def mosaic_on_clip(src, grid):
    t, (h, w), _ = grid
    col0 = int(round((t.c - src.transform.c) / src.transform.a))
    row0 = int(round((src.transform.f - t.f) / -src.transform.e))
    return clean(src.read(1, window=Window(col0, row0, w, h), boundless=True, fill_value=np.nan), src.nodata)


# A coarse grid (or its tiles) reprojected bilinearly onto the clip; in the source's own vertical unit
def raster_on_clip(sources, grid):
    t, (h, w), dst_crs = grid
    hz = CRS.from_wkt(dst_crs.to_wkt())
    out = np.full((h, w), np.nan)
    bounds = (t.c, t.f - h, t.c + w, t.f)
    for src in sources:
        s_crs = src.crs or ASSUMED_CRS
        sb = transform_bounds(dst_crs, s_crs, *bounds)
        if sb[2] < src.bounds.left or sb[0] > src.bounds.right or sb[3] < src.bounds.bottom or sb[1] > src.bounds.top:
            continue
        # Only the patch under the clip, padded by a few source cells for the bilinear kernel
        win = rasterio.windows.from_bounds(*sb, transform=src.transform).round_offsets().round_lengths()
        win = Window(win.col_off - 3, win.row_off - 3, win.width + 6, win.height + 6)
        a = np.full((h, w), np.nan)
        patch = src.read(1, window=win, boundless=True, fill_value=np.nan)
        reproject(clean(patch, src.nodata), a, src_transform=src.window_transform(win), src_crs=s_crs,
                  src_nodata=np.nan, dst_transform=t, dst_crs=hz, dst_nodata=np.nan,
                  resampling=Resampling.bilinear)
        out = np.where(np.isfinite(out), out, a)
    return out


# Vertical unit from a compound CRS's vertical part, when it names one
def stated_unit(src):
    if src.crs is None:
        return None
    wkt = src.crs.to_wkt()
    if "VERT_CS" not in wkt and "VERTCRS" not in wkt:
        return None
    vert = wkt[wkt.find("VERT"):]
    return "ft" if "foot" in vert.lower() else "m"


# The start's frame for one domain: shear, water trim and windows, from its own arrays
def frame(start, gis, windows):
    raw = np.load(tv.npy_dirs(start)[0] / f"domain_{gis}.npy").astype(float)
    a = X.orient_ocean_right(raw, X.OCEAN_LOC)[:, ::-1]
    z = a - X.MHW_M
    z[z < X.WATER_CLAMP_M] = X.SENTINEL_WATER_M
    z[a <= X.RAW_NODATA_MAX_NAVD] = X.NODATA_SENTINEL_M
    above = z > X.BEACH_START_THR_M
    sb = np.where(above.any(axis=1), above.argmax(axis=1), -1)
    zs, _, shear, _ = X.straighten_profiles(z, sb)
    c0, _ = X.water_col_bounds(zs, X.SENTINEL_WATER_M)
    w = windows[f"domain_{gis}"]
    if w.get("trim_offset_c0") is not None and int(w["trim_offset_c0"]) != c0:
        raise SystemExit(f"GIS {gis}: water trim {c0} differs from the picks' {w['trim_offset_c0']}")
    return shear, c0, int(w["i0"]), int(w["i1"])


# Profiles in the start's frame, ocean first, m MHW; nodata at the NODATA sentinel
def profiles(raw_navd, shear, c0):
    a = X.orient_ocean_right(raw_navd, X.OCEAN_LOC)[:, ::-1]
    z = a - X.MHW_M
    z[z < X.WATER_CLAMP_M] = X.SENTINEL_WATER_M
    z[~np.isfinite(a) | (a <= X.RAW_NODATA_MAX_NAVD)] = X.NODATA_SENTINEL_M
    n_along, n_cross = z.shape
    zs = np.full_like(z, X.NODATA_SENTINEL_M)
    for i in range(n_along):
        k = int(shear[i])
        zs[i, :n_cross - k] = z[i, k:]
    return zs[:, c0:]


# find_dunes, with the crest kept only where the whole window was surveyed; m MHW, and why a crest is missing
def crests(z, i0, i1):
    above = z > X.BEACH_START_THR_M
    sb = np.where(above.any(axis=1), above.argmax(axis=1), -1)
    out = np.full(z.shape[0], np.nan)
    edge = np.zeros(z.shape[0], bool)
    status = np.full(z.shape[0], "ok", dtype=object)
    for i in range(z.shape[0]):
        a = max(i0, int(sb[i])) if (X.CLIP_WINDOW_TO_BEACH and sb[i] >= 0) else i0
        b = min(i1, z.shape[1])
        if b <= a:
            status[i] = "beach_landward_of_window"    # the starting dune position is now beach
            continue
        win = z[i, a:b]
        if (win <= X.NODATA_SENTINEL_M + 1e-9).any():
            status[i] = "not_surveyed"
            continue
        valid = win > X.SENTINEL_WATER_M + 1e-9
        if not valid.any():
            status[i] = "water"
            continue
        k = int(np.argmax(np.where(valid, win, -np.inf)))
        out[i], edge[i] = win[k], (a + k == b - 1)
    return out, edge, status


# The check: the start's arrays through this code, exact windows, must give its dune files. Returns the start's
# crests (m NAVD88, domain x cell), the reference the first survey's spikes are judged against
def check_start(start, frames):
    worst = 0.0
    out = np.full((91, 50), np.nan)
    for g in range(1, 91):
        shear, c0, i0, i1 = frames[g]
        z = profiles(np.load(tv.npy_dirs(start)[0] / f"domain_{g}.npy").astype(float), shear, c0)
        above = z > X.BEACH_START_THR_M
        sb = np.where(above.any(axis=1), above.argmax(axis=1), -1)
        e, _ = X.find_dunes(z, sb, i0, i1)
        out[g] = np.asarray(e[:50], float) + X.MHW_M
        h = e - (X.BERM_ELEV_NAVD_M - X.MHW_M)
        h[np.isfinite(h) & (h < 0)] = X.MIN_DUNE_H_M
        h[~np.isfinite(h)] = X.MIN_DUNE_H_M
        ref = np.load(tv.topo_dirs(start)[1] / tv.array_name("dune", g)).ravel()[:50] * 10
        worst = max(worst, float(np.abs(h[:50] - ref).max()))
    if worst > 1e-3:
        raise SystemExit(f"{start} check failed: max |difference| from the dune files {worst:.4f} m")
    print(f"{start} check passed: this code reproduces its dune files (max |difference| {worst:.1e} m)")
    return out


# Vertical unit of a raster survey: stated in its CRS, else from its land elevations against the start's arrays
def unit_of(sources, start, grids):
    stated = {stated_unit(s) for s in sources} - {None}
    if len(stated) == 1:
        return stated.pop(), "stated in the CRS"
    ratios = []
    for g in (10, 30, 50, 70):
        a = downsample(raster_on_clip(sources, grids[g]))
        ref = np.load(tv.npy_dirs(start)[0] / f"domain_{g}.npy").astype(float)
        m = np.isfinite(a) & (ref > 1.0) & (a > 0)
        if m.sum() > 200:
            ratios.append(float(np.median(a[m]) / np.median(ref[m])))
    r = float(np.median(ratios))
    if not (0.7 < r < 1.4 or 2.5 < r < 4.2):
        raise SystemExit(f"cannot tell the vertical unit: land elevation ratio {r:.2f} to the start's arrays")
    return ("ft" if r > 2 else "m"), f"land elevation ratio {r:.2f} to the start's arrays"


# Spikes, once every survey is in: a crest more than SPIKE_M above the cell's crest both before and after it in
# time, an up-and-straight-back-down no dune does. Before: the previous survey with a crest, or for the first survey
# the start's own crest (start_crests, from check_start); after: the next. A survey with nothing on one side is
# judged on the other. A dune that persists, or is lost and stays lost, fails one side. Blanked, status "spike",
# the value kept in crest_raw_m_navd88; the widened column is judged the same way on its own values
def spike_mask(t, col, start_crests):
    w = t.pivot_table(index=["gis", "cell"], columns="survey", values=col, dropna=False).sort_index(axis=1)
    # The start's crest as a column ahead of the first survey, so it stands in until a survey has a crest
    with_start = pd.concat([pd.Series(start_crests[w.index.get_level_values("gis"), w.index.get_level_values("cell")],
                                      index=w.index, name="start"), w], axis=1)
    prev = with_start.ffill(axis=1).shift(1, axis=1)[w.columns]
    nxt = w.bfill(axis=1).shift(-1, axis=1)
    up_prev, up_next = w - prev > SPIKE_M, w - nxt > SPIKE_M
    spike = (up_prev | prev.isna()) & (up_next | nxt.isna()) & (prev.notna() | nxt.notna()) & w.notna()
    s = spike.stack().rename("spike")
    return t.join(s, on=["gis", "cell", "survey"]).spike.fillna(False).astype(bool)


def drop_spikes(t, cover, start_crests):
    t["crest_raw_m_navd88"] = t.crest_m_navd88
    spikes = {col: spike_mask(t, col, start_crests) for col in ("crest_m_navd88", "crest_widened_window_m_navd88")}
    for col, spike in spikes.items():
        t.loc[spike, col] = np.nan
    main_spike = spikes["crest_m_navd88"]
    t.loc[main_spike, "status"] = "spike"
    t.loc[main_spike, "at_landward_edge"] = False
    # Coverage recounted from the cleaned table
    by = t.assign(crest=t.crest_m_navd88.notna(), spike=main_spike, edge=t.at_landward_edge & t.crest_m_navd88.notna())         .groupby(["survey", "gis"])[["crest", "spike", "edge"]].sum()
    key = pd.MultiIndex.from_frame(cover[["survey", "gis"]])
    cover["cells_with_crest"] = by.crest.reindex(key).to_numpy()
    cover["landward_edge_hits"] = by.edge.reindex(key).to_numpy()
    cover["cells_spike"] = by.spike.reindex(key).to_numpy()
    print(f"spikes (> {SPIKE_M:g} m above the crest before and after): "
          f"{t[main_spike].groupby('survey').size().to_dict()}; widened window: "
          f"{t[spikes['crest_widened_window_m_navd88']].groupby('survey').size().to_dict()}", flush=True)
    return t, cover


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    ap.add_argument("--start", choices=list(STARTS), default="2009-start")
    start = ap.parse_args().start
    prod, surveys = STARTS[start]
    windows = json.loads(tv.picks_dir(start).joinpath(
        f"HAT_dune_search_windows_{tv.topo_dirs(start)[2]}.json").read_text(encoding="utf-8"))
    frames = {g: frame(start, g, windows) for g in range(1, 91)}
    start_crests = check_start(start, frames)
    clips = product(prod).gapfill_1m
    grids = {g: clip_grid(clips / f"clip_domain_{g}_filled.tif") for g in range(1, 91)}
    # 1 m road pixels per domain and the 2009 lidar on them, for the road check
    line = gpd.read_file(ROAD_LINE)
    roads, ref09 = {}, {}
    with rasterio.open(REF_2009) as ref:
        for g in range(1, 91):
            t, shape, crs = grids[g]
            roads[g] = rasterize([(x, 1) for x in line.to_crs(crs).geometry], out_shape=shape, transform=t,
                                 fill=0).astype(bool)
            ref09[g] = mosaic_on_clip(ref, grids[g])
    offsets = SURVEY_VERTICAL_OFFSET_M.get(start, {})
    rows, cover, sources_tab = [], [], []
    for year, (kind, path, flown, state, note) in surveys.items():
        files = sorted(path.parent.glob(path.name)) if "*" in path.name else [path]
        srcs = [rasterio.open(f) for f in files]
        try:
            if kind == "mosaic":
                unit, how = "m", "1 m DEM in m NAVD88"
            else:
                unit, how = unit_of(srcs, start, grids)
            scale = US_FT if unit == "ft" else 1.0
            shift = offsets.get(year, 0.0)
            road_diffs = []
            res = srcs[0].res[0]
            sources_tab.append(dict(survey=year, kind=kind, files=len(files), source=str(path),
                                    resolution=f"{res:g} {'ft' if srcs[0].crs and srcs[0].crs.linear_units.lower().startswith(('foot', 'us survey')) else 'm'}",
                                    vertical_unit=unit, unit_from=how, vertical_offset_applied_m=shift,
                                    flown=flown, model_state=state, note=note))
            print(f"{year}: {kind}, {len(files)} file(s), vertical {unit} ({how})", flush=True)
            for g in range(1, 91):
                shear, c0, i0, i1 = frames[g]
                on_clip = mosaic_on_clip(srcs[0], grids[g]) if kind == "mosaic" else raster_on_clip(srcs, grids[g])
                on_m = on_clip * scale + shift
                if g not in RELOCATED_1999:
                    d = (on_m - ref09[g])[roads[g]]
                    road_diffs.append(d[np.isfinite(d)])
                z = profiles(downsample(on_m), shear, c0)
                c, edge, status = crests(z, i0, i1)
                c_wide, _, _ = crests(z, max(i0 - WIDEN_CELLS, 0), i1 + WIDEN_CELLS)
                for cell in range(50):
                    rows.append(dict(survey=year, flown=flown, model_state=state, gis=g, cell=cell,
                                     crest_m_navd88=c[cell] + X.MHW_M,
                                     crest_widened_window_m_navd88=c_wide[cell] + X.MHW_M,
                                     at_landward_edge=bool(edge[cell]), status=status[cell]))
                cover.append(dict(survey=year, gis=g, cells_with_crest=int(np.isfinite(c[:50]).sum()),
                                  cells_not_surveyed=int((status[:50] == "not_surveyed").sum()),
                                  cells_beach_landward_of_window=int((status[:50] == "beach_landward_of_window").sum()),
                                  landward_edge_hits=int(edge[:50][np.isfinite(c[:50])].sum())))
            d = np.concatenate(road_diffs) if road_diffs else np.array([])
            sources_tab[-1]["road_offset_vs_2009_m"] = round(float(np.median(d)), 3) if d.size else np.nan
            sources_tab[-1]["road_px"] = int(d.size)
            print(f"{year}: shifted {shift:+.2f} m; on NC-12 it now sits "
                  f"{sources_tab[-1]['road_offset_vs_2009_m']:+.3f} m from the 2009 lidar ({d.size} px)", flush=True)
        finally:
            for s in srcs:
                s.close()
        n = sum(r["cells_with_crest"] for r in cover if r["survey"] == year)
        print(f"{year}: crest in {n} of 4500 cells ({n / 45:.0f}%)", flush=True)
    t, cover = drop_spikes(pd.DataFrame(rows), pd.DataFrame(cover), start_crests)
    out = tv.product_dir(start) / "survey-crests"
    out.mkdir(parents=True, exist_ok=True)
    t.round(3).to_csv(out / "survey_dune_crests.csv", index=False)
    cover.to_csv(out / "survey_coverage.csv", index=False)
    pd.DataFrame(sources_tab).to_csv(out / "survey_sources.csv", index=False)
    print(out)


if __name__ == "__main__":
    main()
