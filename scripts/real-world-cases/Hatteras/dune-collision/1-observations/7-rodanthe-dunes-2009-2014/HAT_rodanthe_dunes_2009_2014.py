"""
Rodanthe village (GIS 77-81): the 2009, 2011 and 2014 lidar against the nearest aerial imagery.

    python scripts/real-world-cases/Hatteras/dune-collision/1-observations/7-rodanthe-dunes-2009-2014/HAT_rodanthe_dunes_2009_2014.py [GIS,...]

Asked 2026-11-03 while finding why the model over-erodes the Rodanthe village shoreline in 2009-2025: did the dunes
the model starts from match the ground, and were they rebuilt between the 2009 and 2011 surveys? One row per domain:
the March 2008 NOAA mosaic, Hannah's May 2009 Google Earth captures, the 2009, 2011 and 2014 lidar (each on the 2009
datum by HAT_survey_dune_crests.SURVEY_VERTICAL_OFFSET_M: 2011 -0.79 m, 2014 -0.13 m) with a 3 m contour, the January-
April 2014 NOAA mosaic, and the three profiles (median and interquartile range of the domain's 500 east-west 1 m rows,
against the distance from the 2009 shoreline). Hannah's 2009 dune line is drawn on every map. There is no 2010-2012
imagery on disk.

Finding (2026-11-03): the crests barely change from 2009 to 2014 (no rebuild; the model starts from the right dunes);
the beach is cut back 20-40 m in the post-Irene 2011 survey at GIS 79-81 and partly back by 2014.

Writes, in the workspace's dune-collision/1-observations/7-rodanthe-dunes-2009-2014/: rodanthe_dunes_2009_2014.png,
and supporting/rodanthe_dunes_2009_2014.pdf with its caption in supporting/CAPTIONS.md.

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-11-03
"""
import glob, json, sys
from pathlib import Path
import numpy as np
import rasterio
from affine import Affine
from rasterio.crs import CRS
from rasterio.transform import from_origin
from rasterio.warp import reproject, Resampling
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import geopandas as gpd

_HAT_DIR = next(p for p in Path(__file__).resolve().parents if (p / "hatteras_paths.py").exists())
sys.path.insert(0, str(_HAT_DIR))
sys.path.insert(0, str(_HAT_DIR / "input_prep" / "1-barrier3d-domains" / "2-survey-crests"))
from hatteras_paths import WORKSPACE  # noqa: E402
from site_layer.hat_gis_root import GIS_ROOT as _GIS_ROOT  # noqa: E402
import HAT_survey_dune_crests as S  # noqa: E402

UTM = CRS.from_epsg(26918)
GIS_ROOT = str(_GIS_ROOT).replace("\\", "/")
OUT_DIR = WORKSPACE / "dune-collision" / "1-observations" / "7-rodanthe-dunes-2009-2014"
NAME = "rodanthe_dunes_2009_2014"
DEM = {2009: f"{GIS_ROOT}/elevation/dem/2009/noaa_usace2009_nc_dem_Job1076020/2009_full.tif",
       2011: f"{GIS_ROOT}/elevation/dem/2011/noaa_2011_post_irene_NC_J1437890_grid1m/2011_full.tif",
       2014: f"{GIS_ROOT}/elevation/dem/2014/noaa_2014_NOAA_Post_Sandy_DEM_Job1076021/2014_full.tif"}
IMG = {"Mar 2008 (NOAA)": [f"{GIS_ROOT}/aerial_imagery/2008/2008_full_aerial.tif"],
       "May 2009 (Google Earth)": "GE",
       "Jan-Apr 2014 (NOAA)": sorted(glob.glob(f"{GIS_ROOT}/aerial_imagery/2014/noaa_2014_4BandImagery_J1472573/*.tif"))}
DOMS = json.load(open(WORKSPACE / "data" / "hatteras_init" / "5-scr" / "2-transect-frame" / "transect_domains" / "HAT_domains.json"))
LINES = WORKSPACE / "data" / "hatteras_init" / "2-brie-offset" / "dunelines"
GIS = [int(a) for a in sys.argv[1].split(",")] if len(sys.argv) > 1 else [77, 78, 79, 80, 81]
RES = 1.0
# Each survey on the 2009 datum, as the crest script reads it (road check against the 2009 lidar)
OFFSET = {y: S.SURVEY_VERTICAL_OFFSET_M["2009-start"][y] for y in (2011, 2014)}


def box(g):
    f = next(f for f in DOMS["features"] if f["properties"]["domain_id"] == g)
    xs, ys = zip(*f["geometry"]["coordinates"][0])
    return min(xs), min(ys), max(xs), max(ys)


def read_dem(path, bounds):
    x0, y0, x1, y1 = bounds
    w, h = int((x1 - x0) / RES), int((y1 - y0) / RES)
    out = np.full((h, w), np.nan, np.float32)
    with rasterio.open(path) as s:
        reproject(rasterio.band(s, 1), out, src_crs=s.crs or UTM, src_nodata=s.nodata, dst_transform=from_origin(x0, y1, RES, RES),
                  dst_crs=UTM, dst_nodata=np.nan, resampling=Resampling.bilinear)
    out[out < -50] = np.nan
    return out


def ge_files():
    res = []
    for jpg in glob.glob(f"{GIS_ROOT}/aerial_imagery/2009/google_earth_captures/raw_ge/HAT_2009_*.jpg"):
        try:
            a, d, b, e, c, f = [float(v) for v in open(jpg[:-4] + ".jgwx").read().split()]
        except FileNotFoundError:
            continue
        res.append((jpg, Affine(a, b, c, d, e, f)))
    return res


def image(key, bounds, res=0.5):
    x0, y0, x1, y1 = bounds
    w, h = int((x1 - x0) / res), int((y1 - y0) / res)
    t = from_origin(x0, y1, res, res)
    out = np.zeros((3, h, w), np.uint8)
    srcs = ge_files() if IMG[key] == "GE" else [(f, None) for f in IMG[key]]
    for f, aff in srcs:
        with rasterio.open(f) as s:
            a = np.zeros_like(out)
            tr = aff or s.transform
            crs = UTM if aff else s.crs
            L, T_ = tr.c, tr.f
            R, B = tr * (s.width, s.height)
            if max(L, R) < x0 or min(L, R) > x1 or max(T_, B) < y0 or min(T_, B) > y1:
                continue
            for i in range(3):
                reproject(s.read(i + 1), a[i], src_transform=tr, src_crs=crs, src_nodata=0, dst_transform=t, dst_crs=UTM,
                          dst_nodata=0, resampling=Resampling.bilinear)
            empty = (out == 0).all(axis=0)
            out[:, empty] = a[:, empty]
    return out if (out > 0).any() else None


line09 = gpd.read_file(f"{LINES}/duneline_2009.geojson").to_crs(UTM)
ncols = 7
fig, axes = plt.subplots(len(GIS), ncols, figsize=(2.15 * ncols, 3.0 * len(GIS)), constrained_layout=True,
                         gridspec_kw=dict(width_ratios=[1] * 6 + [2.2]))
for r, g in enumerate(GIS):
    bx = box(g)
    full = {y: read_dem(p, bx) + OFFSET.get(y, 0.0) for y, p in DEM.items()}
    med09 = np.nanmedian(full[2009], axis=0)
    xs = bx[0] + (np.arange(full[2009].shape[1]) + 0.5) * RES
    shore = xs[np.flatnonzero(med09 > 0.5).max()]
    win = (shore - 230, bx[1], shore + 50, bx[3])
    ext = (win[0], win[2], win[1], win[3])
    i0, i1 = int((win[0] - bx[0]) / RES), int((win[2] - bx[0]) / RES)
    panels = [("img", "Mar 2008 (NOAA)"), ("img", "May 2009 (Google Earth)"), ("dem", 2009), ("dem", 2011), ("dem", 2014),
              ("img", "Jan-Apr 2014 (NOAA)")]
    for c, (kind, key) in enumerate(panels):
        ax = axes[r, c]
        if kind == "img":
            im = image(key, win)
            if im is None:
                ax.text(0.5, 0.5, "no image", ha="center", transform=ax.transAxes)
            else:
                ax.imshow(np.moveaxis(im, 0, -1), extent=ext)
        else:
            z = full[key][:, i0:i1]
            hm = ax.imshow(z, extent=ext, cmap="terrain", vmin=-1, vmax=7)
            ax.contour(np.linspace(win[0], win[2], z.shape[1]), np.linspace(win[3], win[1], z.shape[0]), z, levels=[3.0],
                       colors="k", linewidths=0.6)
        line09.plot(ax=ax, color="#ff2fd0", lw=1.0)
        ax.set_xlim(win[0], win[2]); ax.set_ylim(win[1], win[3])
        ax.set_xticks([]); ax.set_yticks([])
        if r == 0:
            ax.set_title((f"{key} lidar" + (f" ({OFFSET[key]:+.2f} m)" if key in OFFSET else "")) if kind == "dem" else key, fontsize=8)
        if c == 0:
            ax.set_ylabel(f"GIS {g}", fontsize=9)
    ax = axes[r, -1]
    dist = xs[i0:i1] - shore
    for y, col in ((2009, "#b2182b"), (2011, "#2166ac"), (2014, "#5a5a5a")):
        z = full[y][:, i0:i1]
        q = np.nanpercentile(z, [25, 50, 75], axis=0)
        ax.fill_between(dist, q[0], q[2], color=col, alpha=0.15, lw=0)
        ax.plot(dist, q[1], color=col, lw=1.4, label=f"{y} lidar" + (f" {OFFSET[y]:+.2f} m" if y in OFFSET else "") + " (median, IQR)")
    ax.axhline(1.7, color="k", ls=":", lw=0.8)
    ax.text(dist[0] + 3, 1.75, "model berm 1.7 m", fontsize=6, va="bottom")
    ax.set_xlim(dist[0], dist[-1]); ax.set_ylim(-1, 8)
    ax.set_xlabel("m from 2009 shoreline (+ seaward)", fontsize=7); ax.set_ylabel("Elevation (m NAVD88)", fontsize=7)
    ax.tick_params(labelsize=6)
    if r == 0:
        ax.legend(fontsize=6, frameon=False, loc="upper right")
cb = fig.colorbar(hm, ax=axes[:, 2:5], orientation="horizontal", fraction=0.02, pad=0.01)
cb.set_label("Lidar elevation (m NAVD88); black contour 3 m; pink = Hannah's 2009 dune line", fontsize=7)
(OUT_DIR / "supporting").mkdir(parents=True, exist_ok=True)
fig.savefig(OUT_DIR / f"{NAME}.png", dpi=170)
fig.savefig(OUT_DIR / "supporting" / f"{NAME}.pdf")
print(OUT_DIR / f"{NAME}.png")
