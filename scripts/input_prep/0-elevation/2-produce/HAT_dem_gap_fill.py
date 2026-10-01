"""
Step 1 of 3 (2009-start DEM): clip the 2009 DEM per domain and fill its gaps from the 2014 NOAA Post-Sandy DEM.

    python scripts/input_prep/0-elevation/2-produce/HAT_dem_gap_fill.py

Fills only nodata, under four selection rules (coverage, floor, connectivity,
no value changes). Writes a 1 m filled clip and a survey-year raster per domain
under data/hatteras_init/0-elevation/<product>/1-gapfill-1m/; next is
HAT_dem_resample_clip.py. Details: scripts/input_prep/0-elevation/README.md.

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-09-18
"""

from __future__ import annotations

import csv
import sys
import os
from pathlib import Path

import numpy as np
import geopandas as gpd
import rasterio
from pyproj import CRS
from rasterio.enums import Resampling
from rasterio.vrt import WarpedVRT
from rasterio.windows import Window, from_bounds
from scipy.interpolate import griddata
from scipy.ndimage import binary_dilation, distance_transform_edt, label


# Walk up until a directory holds data/hatteras_init
def _find_project_root(start: Path) -> Path:
    for p in [start, *start.parents]:
        if (p / "data" / "hatteras_init").is_dir():
            return p
    raise SystemExit(f"cannot find data/hatteras_init above {start}")


PROJECT_ROOT = _find_project_root(Path(__file__).resolve())
import sys as _elsys
from pathlib import Path as _ELP
_elsys.path.insert(0, str(next(_q for _q in _ELP(__file__).resolve().parents
                               if (_q / "pyproject.toml").exists()) / "scripts"))
from site_layer import hat_elevation_products as _el  # noqa: E402
# --- CONFIG ------------------------------------------------------------------
ELEVATION_DIR = _el.ELEVATION_ROOT
GIS_ROOT = Path(r"D:\Hatteras_GIS")

BASE_DEM_PATH = (GIS_ROOT / "Elevation" / "Polygons" / "2009"
                 / "usace2009_nc_dem_Job1076020" / "2009_full.tif")

# Fill source: 2014 NOAA Post-Sandy DEM, chosen by measured dry-gap coverage (see README)
FILL_DEM_PATH = (GIS_ROOT / "Elevation" / "Polygons" / "2014"
                 / "2014_NOAA_Post_Sandy_DEM_Job1076021" / "2014_full.tif")
FILL_SOURCE_TAG = "2014_NOAA_PostSandy"
FILL_SOURCE_YEAR = 2014

DOMAIN_FILE = GIS_ROOT / "domains.geojson"
DOMAIN_ID_FIELD = "domain_id"

# The product this script builds
PRODUCT_TAG = "2009-2014"
# Paths from site_layer/hat_elevation_products.py, never built by hand
sys.path.insert(0, str(PROJECT_ROOT / "scripts"))
from site_layer.hat_elevation_products import product as _product  # noqa: E402

OUTPUT_DIR = _product(PRODUCT_TAG, check=False).gapfill_1m
AUDIT_CSV = "gapfill_audit.csv"

GRID_SIZE_M = 10.0    # the eventual Barrier3D cell; the clip must divide by it
EXPECTED_CLIP = (500, 2000)   # rows, cols at 1 m

# vertical datum (matches HAT_dune_topo_extractor / hindcast config)
MHW_ELEVATION = 0.36          # m NAVD88, Duck NC gauge 8651370

# Rule 3: an elevation floor at -2.64 m NAVD88, a guard matching the extractor's water clamp
APPLY_ELEV_FLOOR = True
FILL_MIN_ELEV_NAVD = MHW_ELEVATION - 3.0   # -2.64 m NAVD88, the extractor's clamp

# RULE 2: connectivity
REQUIRE_ISLAND_CONNECTION = True
CONTEXT_BUFFER_M = 200.0      # window padding for connectivity + boundary context

# Bridge gaps up to 20 m (a tidal creek) before the connectivity test; measured, see README
GAP_BRIDGE_M = 20.0

# Rule 4: vertical reconciliation; both value-changing steps off, both still reported
APPLY_BIAS_CORRECTION = False
APPLY_FEATHER = False

OVERLAP_RING_PX = 10          # collar used for the (reported) ring estimate
FEATHER_WIDTH_PX = 5          # blend distance, only if APPLY_FEATHER
MIN_RING_CELLS = 25           # below this the ring estimate is not trustworthy

# How the 2008 ground returns become a 1 m raster
GRID_STAT = "median"

NODATA_OUT = -9999.0

# The _survey raster stores the year each cell's elevation came from, so it needs no legend
SURVEY_2009, SURVEY_NONE = 2009, 0
SURVEY_FILL = FILL_SOURCE_YEAR   # the year written for a filled cell
SURVEY_NODATA = 65535   # unused sentinel; 0 is a real value here, not nodata
# -----------------------------------------------------------------------------


# Geometry - the domain window, snapped to the DEM's own grid

# Polygon bounds as an integer window on the source grid, trimmed to whole blocks; returns the shift
def snap_window(bounds, transform, res_x, res_y, block):
    minx, miny, maxx, maxy = bounds
    left, top = transform.c, transform.f

    col0 = int(round((minx - left) / res_x))
    row0 = int(round((top - maxy) / res_y))
    col1 = int(round((maxx - left) / res_x))
    row1 = int(round((top - miny) / res_y))

    width, height = col1 - col0, row1 - row0
    trim_w, trim_h = width % block, height % block
    width -= trim_w
    height -= trim_h

    adj = {"snap_x_m": abs((left + col0 * res_x) - minx),
           "snap_y_m": abs((top - row0 * res_y) - maxy),
           "trim_x_m": trim_w * res_x, "trim_y_m": trim_h * res_y}
    return Window(col0, row0, width, height), adj


# The window grown by pad_px on every side
def pad_window(win, pad_px):
    return Window(win.col_off - pad_px, win.row_off - pad_px,
                  win.width + 2 * pad_px, win.height + 2 * pad_px)


# Read a domain window boundless, so a domain off the DEM edge gives nodata rather than an error
def read_window(src, win, nodata_in):
    fill = nodata_in if nodata_in is not None else NODATA_OUT
    arr = src.read(1, window=win, boundless=True, fill_value=fill).astype(np.float64)
    if nodata_in is not None and not np.isnan(nodata_in):
        arr = np.where(arr == nodata_in, np.nan, arr)
    return arr


# The fill source DEM

# The candidate DEM, read onto the base DEM's grid on demand
class FillSource:

    def __init__(self, path, dst_crs):
        self.src = rasterio.open(path)
        self.vrt = WarpedVRT(self.src, crs=dst_crs,
                             resampling=Resampling.nearest)
        crs = CRS.from_wkt(self.src.crs.to_wkt()) if self.src.crs else None
        hz = crs.sub_crs_list[0] if (crs is not None and crs.is_compound) else crs
        self.epsg = hz.to_epsg() if hz is not None else None
        self.res = self.src.transform.a
        self.nodata = self.src.nodata

    # Read onto the base grid; WarpedVRT rejects boundless reads, so overlap is pasted at its offset
    def read_on_grid(self, bounds, shape):
        bx0, by0, bx1, by1 = bounds
        out = np.full(shape, np.nan, np.float32)
        vb = self.vrt.bounds
        ix0, iy0 = max(bx0, vb.left), max(by0, vb.bottom)
        ix1, iy1 = min(bx1, vb.right), min(by1, vb.top)
        if ix1 <= ix0 or iy1 <= iy0:
            return out
        h, w = int(round(iy1 - iy0)), int(round(ix1 - ix0))
        if h <= 0 or w <= 0:
            return out
        d = self.vrt.read(1, window=from_bounds(ix0, iy0, ix1, iy1,
                                                transform=self.vrt.transform),
                          out_shape=(h, w), masked=True,
                          resampling=Resampling.nearest)
        d = np.ma.filled(d.astype(np.float32), np.nan)
        r0, c0 = max(int(round(by1 - iy1)), 0), max(int(round(ix0 - bx0)), 0)
        hh, ww = min(h, shape[0] - r0), min(w, shape[1] - c0)
        if hh > 0 and ww > 0:
            out[r0:r0 + hh, c0:c0 + ww] = d[:hh, :ww]
        return out

    def close(self):
        self.vrt.close()
        self.src.close()


# The fill

# Candidate cells reachable from the island through valid ground or other candidates
def island_connected(valid, candidate, bridge_px=0):
    if not valid.any():
        return np.zeros_like(candidate)
    lab_v, n_v = label(valid, structure=np.ones((3, 3)))
    if n_v == 0:
        return np.zeros_like(candidate)
    sizes = np.bincount(lab_v.ravel()); sizes[0] = 0
    island = lab_v == int(np.argmax(sizes))

    land = valid | candidate
    if bridge_px > 0:
        land = binary_dilation(land, iterations=int(bridge_px))
    lab_a, _ = label(land, structure=np.ones((3, 3)))
    keep = set(np.unique(lab_a[island])) - {0}
    return candidate & np.isin(lab_a, list(keep))


# z[a] - z[b] over every 4-connected pair with a in mask_a, b in mask_b
def _neighbour_diffs(z, mask_a, mask_b):
    out = []
    for sl_a, sl_b in ((np.s_[:, :-1], np.s_[:, 1:]),
                       (np.s_[:, 1:], np.s_[:, :-1]),
                       (np.s_[:-1, :], np.s_[1:, :]),
                       (np.s_[1:, :], np.s_[:-1, :])):
        m = mask_a[sl_a] & mask_b[sl_b]
        if m.any():
            out.append(z[sl_a][m] - z[sl_b][m])
    return np.concatenate(out) if out else np.empty(0)


# How big a step does the fill create where it meets measured 2009 ground? Reported against a CONTROL
def seam_check(z, measured, filled):
    seam = _neighbour_diffs(z, filled, measured)
    ctrl = _neighbour_diffs(z, measured, measured)
    # Every key present every time, so the audit CSV keeps stable columns
    res = {"seam_n": int(seam.size), "ctrl_n": int(ctrl.size),
           "seam_median_signed_m": None, "seam_median_abs_m": None,
           "seam_p90_abs_m": None, "ctrl_median_abs_m": None,
           "ctrl_p90_abs_m": None, "seam_over_ctrl": None}
    if seam.size:
        res["seam_median_signed_m"] = float(np.median(seam))
        res["seam_median_abs_m"] = float(np.median(np.abs(seam)))
        res["seam_p90_abs_m"] = float(np.percentile(np.abs(seam), 90))
    if ctrl.size:
        res["ctrl_median_abs_m"] = float(np.median(np.abs(ctrl)))
        res["ctrl_p90_abs_m"] = float(np.percentile(np.abs(ctrl), 90))
    if res["seam_median_abs_m"] is not None and res["ctrl_median_abs_m"]:
        res["seam_over_ctrl"] = round(
            res["seam_median_abs_m"] / res["ctrl_median_abs_m"], 3)
    return res


# Median base - fill over a ring around the fill, and how many cells it used
def estimate_bias(base, fill, fill_mask, ring_px):
    ring = binary_dilation(fill_mask, iterations=ring_px) & ~fill_mask
    both = ring & ~np.isnan(base) & ~np.isnan(fill)
    n = int(both.sum())
    if n < MIN_RING_CELLS:
        return 0.0, n
    return float(np.nanmedian(base[both] - fill[both])), n


# Continuation of the 2009 surface into the fill area: the blend target at the seam, so no hard step
def boundary_extrapolation(base, fill_mask, buffer_px):
    dil = binary_dilation(fill_mask, iterations=buffer_px)
    border = dil & ~fill_mask & ~np.isnan(base)
    br, bc = np.where(border)
    fr, fc = np.where(fill_mask)
    out = np.full(base.shape, np.nan)
    if len(br) < 10 or len(fr) == 0:
        return out
    vals = griddata((br, bc), base[br, bc], (fr, fc), method="linear")
    bad = np.isnan(vals)
    if bad.any():
        vals[bad] = griddata((br, bc), base[br, bc], (fr[bad], fc[bad]),
                             method="nearest")
    out[fr, fc] = vals
    return out


# Blend the fill into the base over feather_px from the fill edge
def feathered_merge(base, fill_vals, fill_mask, extrap, feather_px):
    dist = distance_transform_edt(fill_mask)
    w = np.clip(dist / max(feather_px, 1e-9), 0, 1)
    blended = np.where(np.isnan(extrap), fill_vals,
                       (1 - w) * extrap + w * fill_vals)
    out = base.copy()
    out[fill_mask] = blended[fill_mask]
    return out


# Io

# Write one single-band GeoTIFF
def write_raster(arr, transform, crs, path, dtype, nodata):
    profile = {"driver": "GTiff", "height": arr.shape[0], "width": arr.shape[1],
               "count": 1, "dtype": dtype, "crs": crs, "transform": transform,
               "nodata": nodata, "compress": "deflate"}
    with rasterio.open(path, "w", **profile) as dst:
        dst.write(arr.astype(dtype), 1)


# The DEM's compound CRS compared properly, since a plain `!=` reports a no-op reprojection
def resolve_crs(src, gdf):
    if src.crs is None:
        print("  WARNING: DEM has no CRS; assuming domains already match.")
        return gdf
    dem_crs = CRS.from_wkt(src.crs.to_wkt())
    hz = dem_crs.sub_crs_list[0] if dem_crs.is_compound else dem_crs
    if dem_crs.is_compound:
        print(f"  DEM CRS compound: horizontal EPSG:{hz.to_epsg()}, "
              f"vertical {dem_crs.sub_crs_list[1].name}")
    if gdf.crs is None:
        print("  WARNING: domains have no CRS; assuming they match the DEM.")
    elif gdf.crs.to_epsg() != hz.to_epsg():
        print(f"  Reprojecting domains {gdf.crs} -> EPSG:{hz.to_epsg()}")
        gdf = gdf.to_crs(hz)
    else:
        print(f"  Domains already EPSG:{hz.to_epsg()} - no reprojection")
    return gdf


# Run: per domain, clip, fill, check the seams, write the rasters, then the audit
def main():
    for label_, path in (("base DEM", BASE_DEM_PATH),
                         ("domain file", DOMAIN_FILE),
                         ("fill source", FILL_DEM_PATH)):
        if not path.exists():
            raise FileNotFoundError(
                f"{label_} not found: {path}\n"
                f"  Fix the paths at the top of this script.")

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    src = rasterio.open(BASE_DEM_PATH)
    res_x, res_y = src.transform.a, -src.transform.e
    print(f"\n2009 DEM: {src.height} x {src.width}, {res_x} m, nodata={src.nodata}")
    block = int(round(GRID_SIZE_M / res_x))
    pad_px = int(round(CONTEXT_BUFFER_M / res_x))

    gdf = gpd.read_file(DOMAIN_FILE)
    print(f"\n{len(gdf)} domains from {DOMAIN_FILE.name}")
    gdf = resolve_crs(src, gdf)

    fill = FillSource(FILL_DEM_PATH, src.crs)
    print(f"\nfill source: {FILL_DEM_PATH.name}  ({FILL_SOURCE_TAG})")
    print(f"  EPSG:{fill.epsg}  {fill.res} m  nodata={fill.nodata}"
          f"{'  -> reprojected on read' if fill.epsg else ''}")

    print(f"\nRules: coverage={FILL_SOURCE_YEAR} DEM | "
          f"connectivity={REQUIRE_ISLAND_CONNECTION} (bridge {GAP_BRIDGE_M:g} m)"
          f" | floor={'%.2f m NAVD88' % FILL_MIN_ELEV_NAVD if APPLY_ELEV_FLOOR else 'off'}"
          f" | bias={APPLY_BIAS_CORRECTION} feather={APPLY_FEATHER}")

    audit = []
    for _, row in gdf.iterrows():
        dom = row[DOMAIN_ID_FIELD]
        try:
            dom = int(dom)
        except (TypeError, ValueError):
            pass

        win, adj = snap_window(row.geometry.bounds, src.transform, res_x, res_y, block)
        big = pad_window(win, pad_px)

        base_b = read_window(src, big, src.nodata)
        t_big = src.window_transform(big)

        bx0, by1 = t_big.c, t_big.f
        bx1 = bx0 + big.width * res_x
        by0 = by1 - big.height * res_y
        g2008 = fill.read_on_grid((bx0, by0, bx1, by1), base_b.shape)

        gap = np.isnan(base_b)
        valid = ~gap
        has08 = ~np.isnan(g2008)

        cand = gap & has08
        n_cand = int(cand.sum())

        if REQUIRE_ISLAND_CONNECTION:
            # bridging a gap of GAP_BRIDGE_M needs half that dilation per side
            bridge_px = int(round(GAP_BRIDGE_M / res_x / 2.0))
            conn = island_connected(valid, cand, bridge_px)
        else:
            conn = cand
        n_conn_drop = n_cand - int(conn.sum())

        # Both bias estimates reported every run, applied or not
        bias_ring, n_ring = estimate_bias(base_b, g2008, conn, OVERLAP_RING_PX)
        ov = valid & has08
        n_ov = int(ov.sum())
        bias_overlap = (float(np.nanmedian(base_b[ov] - g2008[ov]))
                        if n_ov >= MIN_RING_CELLS else 0.0)

        bias = bias_ring if APPLY_BIAS_CORRECTION else 0.0
        fill_vals = g2008 + bias

        n_floor_drop = 0
        # Counted on the cropped domain window, to compare with `filled`
        _c = np.s_[pad_px:pad_px + win.height, pad_px:pad_px + win.width]
        n_sub_mhw = int((conn[_c] & (fill_vals[_c] < MHW_ELEVATION)).sum())
        if APPLY_ELEV_FLOOR:
            too_low = conn & (fill_vals < FILL_MIN_ELEV_NAVD)
            n_floor_drop = int(too_low.sum())
            conn = conn & ~too_low

        n_fill = int(conn.sum())
        if n_fill and APPLY_FEATHER:
            extrap = boundary_extrapolation(base_b, conn, pad_px)
            merged = feathered_merge(base_b, fill_vals, conn, extrap, FEATHER_WIDTH_PX)
        elif n_fill:
            merged = base_b.copy()
            merged[conn] = fill_vals[conn]
        else:
            merged = base_b

        seam = seam_check(merged, valid, conn)

        # back to the exact domain window
        r0 = pad_px; c0 = pad_px
        out = merged[r0:r0 + win.height, c0:c0 + win.width]
        filled_mask = conn[r0:r0 + win.height, c0:c0 + win.width]
        gap_out = np.isnan(out)

        survey = np.full(out.shape, SURVEY_2009, np.uint16)
        survey[filled_mask] = SURVEY_FILL
        survey[gap_out] = SURVEY_NONE

        if out.shape != EXPECTED_CLIP:
            print(f"  domain {dom}: ERROR clip {out.shape}, expected {EXPECTED_CLIP}")

        t_out = src.window_transform(win)
        write_raster(np.where(np.isnan(out), NODATA_OUT, out), t_out, src.crs,
                     OUTPUT_DIR / f"clip_domain_{dom}_filled.tif",
                     "float32", NODATA_OUT)
        write_raster(survey, t_out, src.crs,
                     OUTPUT_DIR / f"clip_domain_{dom}_survey.tif",
                     "uint16", SURVEY_NODATA)

        gap_before = int(np.isnan(base_b[r0:r0 + win.height, c0:c0 + win.width]).sum())
        n_filled_out = int(filled_mask.sum())
        print(f"  domain {dom:>3}: nodata {gap_before:6d} -> {int(gap_out.sum()):6d}  "
              f"filled {n_filled_out:5d}  "
              f"seam {seam.get('seam_median_abs_m', float('nan')):.3f} vs ctrl "
              f"{seam.get('ctrl_median_abs_m', float('nan')):.3f} m "
              f"(x{seam.get('seam_over_ctrl', float('nan')):.2f})  "
              f"bias ring {bias_ring:+.3f} / overlap {bias_overlap:+.3f}  "
              f"drop: conn {n_conn_drop}, floor {n_floor_drop}", flush=True)

        audit.append({
            "domain": dom, "rows": out.shape[0], "cols": out.shape[1],
            "nodata_before": gap_before, "nodata_after": int(gap_out.sum()),
            "filled": n_filled_out,
            "cand_2008_coverage": n_cand,
            "dropped_connectivity": n_conn_drop,
            "dropped_elev_floor": n_floor_drop,
            "cand_below_mhw": n_sub_mhw,
            "bias_applied_m": round(bias, 4),
            "bias_ring_m": round(bias_ring, 4), "bias_ring_cells": n_ring,
            "bias_overlap_m": round(bias_overlap, 4), "bias_overlap_cells": n_ov,
            **{k: (round(v, 4) if isinstance(v, float) else v)
               for k, v in seam.items()},
            "fill_cells_in_window": int(has08.sum()),
            "snap_x_m": round(adj["snap_x_m"], 4), "snap_y_m": round(adj["snap_y_m"], 4),
            "trim_x_m": adj["trim_x_m"], "trim_y_m": adj["trim_y_m"],
        })

    src.close()
    fill.close()

    path = OUTPUT_DIR / AUDIT_CSV
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(audit[0].keys()))
        w.writeheader(); w.writerows(audit)

    tot_before = sum(a["nodata_before"] for a in audit)
    tot_after = sum(a["nodata_after"] for a in audit)
    tot_fill = sum(a["filled"] for a in audit)
    tot_conn = sum(a["dropped_connectivity"] for a in audit)
    tot_floor = sum(a["dropped_elev_floor"] for a in audit)
    print(f"\n{'=' * 70}\n{len(audit)} domains")
    print(f"  nodata cells   {tot_before:,} -> {tot_after:,}  "
          f"({100 * (tot_before - tot_after) / max(tot_before, 1):.1f}% recovered)")
    print(f"  filled         {tot_fill:,}")
    print(f"  rejected by connectivity {tot_conn:,}   by elevation floor {tot_floor:,}")
    print(f"  audit: {path}")
    print("\nNext: HAT_dem_resample_clip.py resamples these 1 m clips to 10 m.")


if __name__ == "__main__":
    main()
