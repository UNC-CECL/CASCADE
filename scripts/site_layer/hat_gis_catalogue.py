"""
Measure every dataset in the GIS drive's catalogue, and how much of each domain's land it covers.

    python scripts/site_layer/hat_gis_catalogue.py               # every row
    python scripts/site_layer/hat_gis_catalogue.py --ids A1996   # just these rows

Reads catalogue.csv at the drive root (one row per dataset: the source facts, typed in by hand) and
fills its measured columns (files, size, CRS, resolution, extent, domains covered). Writes
coverage_by_domain.csv beside it: per dataset and domain, the share of the domain's land (2009 DEM
above MHW, on a 100 m grid) where the dataset has data. Point clouds count their order footprint and
shoreline layers count a domain they cross. Details: D:\\hatteras_ch1\\gis_sources\\README.md.

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-10-21
"""
from __future__ import annotations

import argparse
import re
import sys
import warnings
from pathlib import Path

import geopandas as gpd
import numpy as np
import pandas as pd
import pyogrio
import rasterio
from affine import Affine
from pyproj import CRS, Transformer
from rasterio.windows import Window
from shapely.geometry import Point
from shapely.ops import unary_union

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from site_layer.hat_gis_root import GIS_ROOT  # noqa: E402

warnings.filterwarnings("ignore", message="Measured")
# JPG frames georeferenced in ArcGIS: georef() reads their world files instead
warnings.filterwarnings("ignore", category=rasterio.errors.NotGeoreferencedWarning)

# --- CONFIG ------------------------------------------------------------------
CATALOGUE = GIS_ROOT / "catalogue.csv"
COVERAGE = GIS_ROOT / "coverage_by_domain.csv"
DOMAINS = GIS_ROOT / "domains.geojson"
LAND_DEM = GIS_ROOT / "elevation" / "dem" / "2009" / "noaa_usace2009_nc_dem_Job1076020" / "2009_full.tif"
MHW_NAVD88 = 0.36            # m; land is above this in the 2009 DEM
GRID_M = 100.0               # sample spacing inside each domain
MIN_LAND_POINTS = 5          # fewer land points than this: coverage left blank
COVERED_SHARE = 0.5          # a domain counts as covered at this share of its land
RASTER_EXT = {".tif", ".tiff", ".jpg", ".jpeg", ".png"}
WORLD_EXT = (".jgwx", ".jgw", ".pgwx", ".pgw", ".tfwx", ".tfw", ".wld")
SKIP_PARTS = ("thumbnail", "_reorg", "frames_zipped")
MEASURED = ["files", "size_gb", "crs", "resolution_m", "bands", "lon_min", "lat_min", "lon_max", "lat_max",
            "domains_covered", "coverage_method"]
# -----------------------------------------------------------------------------


# Land sample points of every domain, in the land DEM's CRS: {domain_id: (x array, y array)}
def land_points():
    dom = gpd.read_file(DOMAINS)
    with rasterio.open(LAND_DEM) as dem:
        dom = dom.to_crs(dem.crs)
        pts, crs = {}, dem.crs
        for _, r in dom.iterrows():
            x0, y0, x1, y1 = r.geometry.bounds
            xs, ys = np.meshgrid(np.arange(x0 + GRID_M / 2, x1, GRID_M), np.arange(y0 + GRID_M / 2, y1, GRID_M))
            inside = [(x, y) for x, y in zip(xs.ravel(), ys.ravel()) if r.geometry.contains(Point(x, y))]
            z = np.array([v[0] for v in dem.sample(inside)], dtype=float) if inside else np.array([])
            nod = dem.nodata
            land = [(x, y) for (x, y), v in zip(inside, z) if np.isfinite(v) and v != nod and v > MHW_NAVD88]
            pts[int(r["domain_id"])] = (np.array([p[0] for p in land]), np.array([p[1] for p in land]))
    return pts, crs, dom


# A raster's transform and CRS, reading an ArcGIS world file and .aux.xml where GDAL does not
def georef(path: Path, src):
    if src.crs is not None and not src.transform.is_identity:
        return src.transform, src.crs
    for ext in WORLD_EXT:
        w = path.with_suffix(ext)
        if w.exists():
            a, d, b, e, c, f = [float(v) for v in w.read_text().split()[:6]]
            aux = path.with_name(path.name + ".aux.xml")
            m = re.search(r"<(?:WKT|SRS)>(.*?)</(?:WKT|SRS)>", aux.read_text(errors="ignore"), re.S) if aux.exists() else None
            crs = CRS.from_wkt(m.group(1).replace("&quot;", '"')) if m else src.crs
            return Affine(a, b, c, d, e, f), crs
    return None, None


# Short CRS label: EPSG code where there is one (horizontal part of a compound CRS), else its name
def crs_label(c: CRS) -> str:
    h = c.sub_crs_list[0] if c.sub_crs_list else c
    e = h.to_epsg(min_confidence=70)
    return f"EPSG:{e}" if e else h.name


# Data files of a dataset path (a file, or a folder searched recursively)
def data_files(p: Path, kind: str):
    if p.is_file():
        return [p]
    if kind in ("raster", "elevation"):
        return sorted(f for f in p.rglob("*") if f.suffix.lower() in RASTER_EXT
                      and not any(s in str(f).lower() for s in SKIP_PARTS) and not f.name.startswith("._"))
    return []


# Share of each domain's land points where any of the rasters has data
def raster_coverage(files, kind, pts, pts_crs):
    opened, meta = [], {"crs": "", "res": [], "bands": 0, "bounds": []}
    for f in files:
        src = rasterio.open(f)
        tr, crs = georef(f, src)
        if tr is None or crs is None:
            src.close()
            continue
        c = CRS.from_user_input(crs)
        unit = 1.0
        if c.is_geographic:
            unit = 111_320.0 * np.cos(np.radians(35.4))
        elif any(w in c.axis_info[0].unit_name.lower() for w in ("foot", "feet")):
            unit = 0.3048006
        meta["res"].append(np.hypot(tr.a, tr.d) * unit)     # pixel width, also for rotated frames
        meta["crs"] = meta["crs"] or crs_label(c)
        meta["bands"] = max(meta["bands"], src.count)
        b = rasterio.transform.array_bounds(src.height, src.width, tr)      # west, south, east, north
        meta["bounds"].append(Transformer.from_crs(crs, 4326, always_xy=True).transform_bounds(*b))
        opened.append((src, tr, crs, b))
    cov = {}
    trans = {}
    for did, (xs, ys) in pts.items():
        if len(xs) < MIN_LAND_POINTS:
            cov[did] = np.nan
            continue
        hit = np.zeros(len(xs), bool)
        for src, tr, crs, (b_w, b_s, b_e, b_n) in opened:
            key = CRS.from_user_input(crs).to_string()
            if key not in trans:
                trans[key] = Transformer.from_crs(pts_crs, crs, always_xy=True)
            fx, fy = trans[key].transform(xs, ys)
            todo = np.where(~hit & (fx > b_w) & (fx < b_e) & (fy > b_s) & (fy < b_n))[0]
            if not len(todo):
                continue
            inv = ~tr
            nb = min(3, src.count)
            for i in todo:
                c, r = inv * (fx[i], fy[i])
                c, r = int(c), int(r)
                if not (0 <= r < src.height and 0 <= c < src.width):
                    continue
                v = src.read(list(range(1, nb + 1)), window=Window(c, r, 1, 1))[:, 0, 0].astype(float)
                if kind == "elevation":
                    ok = np.isfinite(v[0]) and v[0] != src.nodata and v[0] > -1e5
                else:
                    ok = not (src.nodata is not None and np.all(v == src.nodata)) and np.any(v > 0)
                hit[i] = ok
            if hit.all():
                break
        cov[did] = hit.mean()
    for src, *_ in opened:
        src.close()
    return cov, meta


# Share of each domain's land points inside a point-cloud order's footprint (tile index, clipped)
def footprint_coverage(p: Path, pts, pts_crs):
    tiles = list(p.glob("*_tileindex.shp"))
    clip = list(p.glob("*_clippoly.shp"))
    if not tiles:
        return {d: np.nan for d in pts}, {}
    poly = unary_union(gpd.read_file(tiles[0]).to_crs(pts_crs).geometry)
    if clip:
        poly = poly.intersection(unary_union(gpd.read_file(clip[0]).to_crs(pts_crs).geometry))
    cov = {d: (np.mean([poly.contains(Point(x, y)) for x, y in zip(xs, ys)]) if len(xs) >= MIN_LAND_POINTS else np.nan)
           for d, (xs, ys) in pts.items()}
    lon0, lat0, lon1, lat1 = gpd.GeoSeries([poly], crs=pts_crs).to_crs(4326).total_bounds
    return cov, {"crs": "point cloud", "bounds": [(lon0, lat0, lon1, lat1)]}


# 1 where a shoreline layer crosses the domain, else 0 (every layer of a file or folder)
def line_coverage(p: Path, dom):
    files = [p] if p.is_file() else sorted(list(p.glob("*.shp")) + list(p.glob("*.gpkg")))
    geoms, crs_s = [], ""
    for f in files:
        for layer, _ in pyogrio.list_layers(f):
            g = gpd.read_file(f, layer=layer)
            crs_s = crs_s or crs_label(CRS.from_user_input(g.crs))
            geoms.append(unary_union(g.to_crs(dom.crs).geometry))
    allg = unary_union(geoms)
    cov = {int(r["domain_id"]): float(allg.intersects(r.geometry)) for _, r in dom.iterrows()}
    lon0, lat0, lon1, lat1 = gpd.GeoSeries([allg], crs=dom.crs).to_crs(4326).total_bounds
    return cov, {"crs": crs_s, "bounds": [(lon0, lat0, lon1, lat1)]}


# Domain list as ranges: 1-15, 20, 24-90
def ranges(ids):
    ids, out = sorted(ids), []
    for i in ids:
        if out and i == out[-1][1] + 1:
            out[-1][1] = i
        else:
            out.append([i, i])
    return ", ".join(f"{a}-{b}" if a != b else f"{a}" for a, b in out)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ids", nargs="*", help="catalogue ids to measure (default: all)")
    args = ap.parse_args()

    cat = pd.read_csv(CATALOGUE, dtype=str, keep_default_na=False)
    for c in MEASURED:
        if c not in cat.columns:
            cat[c] = ""
    old = pd.read_csv(COVERAGE, dtype={"id": str}).set_index("id") if COVERAGE.exists() else pd.DataFrame()
    pts, pts_crs, dom = land_points()
    print(f"land points: {sum(len(v[0]) for v in pts.values())} in {len(pts)} domains")

    # Measure each row
    rows = {}
    for i, r in cat.iterrows():
        if args.ids and r["id"] not in args.ids:
            continue
        p = GIS_ROOT / r["path"]
        if not p.exists():
            print(f"  {r['id']}: MISSING {p}")
            cat.loc[i, "coverage_method"] = "path missing"
            continue
        allf = [p] if p.is_file() else [f for f in p.rglob("*") if f.is_file()]
        cat.loc[i, "files"] = str(len(allf))
        cat.loc[i, "size_gb"] = f"{sum(f.stat().st_size for f in allf) / 1e9:.2f}"
        kind = r["kind"]
        cov, meta, method = None, {}, ""
        if kind in ("raster", "elevation"):
            files = data_files(p, kind)
            cov, meta = raster_coverage(files, kind, pts, pts_crs)
            method = "pixels with data, share of land points" if meta.get("res") else "no georeferenced rasters"
            if not meta.get("res"):
                cov = None
        elif kind == "pointcloud":
            cov, meta = footprint_coverage(p, pts, pts_crs)
            method = "order footprint (tile index x clip polygon), share of land points"
        elif kind == "line":
            cov, meta = line_coverage(p, dom)
            method = "1 if a line crosses the domain"
        if meta.get("crs"):
            cat.loc[i, "crs"] = meta["crs"]
        if meta.get("res"):
            cat.loc[i, "resolution_m"] = f"{np.median(meta['res']):.2f}"
            cat.loc[i, "bands"] = str(meta["bands"])
        if meta.get("bounds"):
            b = np.array(meta["bounds"])
            for k, v in zip(("lon_min", "lat_min", "lon_max", "lat_max"),
                            (b[:, 0].min(), b[:, 1].min(), b[:, 2].max(), b[:, 3].max())):
                cat.loc[i, k] = f"{v:.4f}"
        cat.loc[i, "coverage_method"] = method
        if cov is not None:
            rows[r["id"]] = cov
            cat.loc[i, "domains_covered"] = ranges([d for d, s in cov.items() if np.isfinite(s) and s >= COVERED_SHARE])
        print(f"  {r['id']:<14} {kind:<10} {cat.loc[i, 'domains_covered'] or '-'}")

    # Write both tables (coverage rows of unmeasured ids are kept)
    cat.to_csv(CATALOGUE, index=False)
    new = pd.DataFrame.from_dict(rows, orient="index").round(2)
    new.columns = [f"d{c:02d}" for c in new.columns]
    out = pd.concat([old.drop(index=new.index, errors="ignore"), new])
    out = out.reindex([i for i in cat["id"] if i in out.index])
    out.index.name = "id"
    out.to_csv(COVERAGE)
    print(f"wrote {CATALOGUE} and {COVERAGE}")


if __name__ == "__main__":
    main()
