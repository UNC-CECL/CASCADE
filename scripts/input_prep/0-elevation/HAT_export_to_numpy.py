"""
HAT_export_to_numpy.py

Step 3 of 3: converts each 10 m domain raster into the .npy array that
HAT_dune_topo_extractor.py reads, matching the documented ArcGIS
RasterToNumPyArray convention:
    - nodata cells filled with -10 (not NaN)
    - no unit conversion - stays in METRES NAVD88
    - no axis transpose - rasterio's row-major read uses the same
      north-at-top / west-at-left convention as arcpy.RasterToNumPyArray

THE CONTRACT THIS HAS TO SATISFY
---------------------------------
Read out of HAT_dune_topo_extractor.py rather than assumed:

  LOAD_PATH (line 142)   INIT_ROOT/1-barrier3d-domains/{YEAR}-raw/
                         {YEAR}-npy-arrays/{DEM_NAME}
  filenames (line 2557)  startswith("domain_") and endswith(".npy")
  load      (line 994)   np.load(...).astype(float), must be 2D
  nodata    (line 1015)  raw <= RAW_NODATA_MAX_NAVD (-9.0); raw nodata is
                         exactly -10.0 m NAVD88
  units     (line 1017)  z = raw - MHW_M, so raw must be m NAVD88
  shape                  ALONG_COLS=50 alongshore, TOPO_ROWS=200 cross-shore,
                         OCEAN_LOC="right" -> ocean at the HIGH column index

Our 10 m rasters are 50 rows (alongshore, 500 m) x 200 cols (cross-shore,
2000 m) with east at the high column index, and east is the ocean side. So the
arrays go through as-is: no transpose, no flip, no unit change.

THE FILENAME IS NOT domain_N_topography_2009.npy
-------------------------------------------------
That is what the extractor WRITES. What it READS is domain_<N>.npy. Getting
this backwards produces a folder the extractor silently finds zero domains in.

WHY A SEPARATE DEM_NAME FOLDER
-------------------------------
Output goes to {YEAR}-npy-arrays/2009_pea_hatteras_filled/, not over
2009_pea_hatteras. Dune picks are keyed per run/version, so leaving the
un-filled arrays in place keeps the existing 2009_v4 picks valid and lets you
A/B the road-drowning result. Point DEM_NAME at the new folder for a fresh run.

THE SURVEY ARRAYS GO IN A SIBLING FOLDER, DELIBERATELY
-------------------------------------------------------
The extractor globs domain_*.npy. A survey array named domain_5_survey.npy would
match that glob and be loaded as if it were a domain. So the survey arrays go to
a separate directory that the extractor never looks at, under the same
domain_<N>.npy name.

INPUTS  (data/hatteras_init/0-elevation/2-resampled-10m/)
    resampled_domain_<N>_filled.tif
    resampled_domain_<N>_survey.tif

OUTPUTS
    data/hatteras_init/1-barrier3d-domains/2009-raw/2009-npy-arrays/
        2009_pea_hatteras_filled/domain_<N>.npy          m NAVD88, -10 nodata
        2009_pea_hatteras_filled_survey/domain_<N>.npy   2009 / fill year / 0

Requires: rasterio, numpy
"""

import csv
import re
from pathlib import Path

import numpy as np
import rasterio

# =============================================================================
# CONFIG
# =============================================================================

PROJECT_ROOT = Path(__file__).resolve().parents[3]
INIT_ROOT = PROJECT_ROOT / "data" / "hatteras_init"
ELEVATION_DIR = INIT_ROOT / "0-elevation"

SOURCE_TAG = "2014_NOAA_PostSandy"   # match HAT_dem_gap_fill.py
INPUT_DIR = ELEVATION_DIR / "2-resampled-10m" / SOURCE_TAG
INPUT_GLOB = "resampled_domain_*_filled.tif"
ID_PATTERN = re.compile(r"resampled_domain_(\w+)_filled\.tif$")

DEM_YEAR = "2009"
DEM_NAME = "2009_pea_hatteras_filled"   # set DEM_NAME in the extractor to match

NPY_ROOT = INIT_ROOT / "1-barrier3d-domains" / f"{DEM_YEAR}-raw" / f"{DEM_YEAR}-npy-arrays"
OUTPUT_DIR = NPY_ROOT / DEM_NAME
SURVEY_DIR = NPY_ROOT / f"{DEM_NAME}_survey"

NODATA_FILL = -10.0   # the extractor detects nodata as raw <= -9.0
EXPECTED_SHAPE = (50, 200)

AUDIT_CSV = "export_audit.csv"

SURVEY_2009, SURVEY_NONE = 2009, 0
# Keep in sync with FILL_SOURCE_YEAR in HAT_dem_gap_fill.py.
SURVEY_FILL = 2014


def main():
    if not INPUT_DIR.exists():
        raise FileNotFoundError(
            f"{INPUT_DIR} not found - run HAT_dem_resample_clip.py first.")
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    SURVEY_DIR.mkdir(parents=True, exist_ok=True)

    paths = sorted(INPUT_DIR.glob(INPUT_GLOB))
    if not paths:
        raise FileNotFoundError(f"no {INPUT_GLOB} in {INPUT_DIR}")
    print(f"{len(paths)} domain raster(s) in {INPUT_DIR}")
    print(f"-> {OUTPUT_DIR}")

    audit, n_bad = [], 0
    for p in paths:
        m = ID_PATTERN.search(p.name)
        if not m:
            print(f"  WARNING: cannot parse domain id from {p.name} - skipping")
            continue
        dom = m.group(1)

        with rasterio.open(p) as src:
            arr = src.read(1).astype(np.float32)
            nd = src.nodata
        if nd is not None and not np.isnan(nd):
            arr = np.where(arr == nd, np.nan, arr)
        out = np.where(np.isnan(arr), NODATA_FILL, arr).astype(np.float32)

        if EXPECTED_SHAPE is not None and out.shape != tuple(EXPECTED_SHAPE):
            n_bad += 1
            print(f"  domain {dom}: ERROR shape {out.shape}, expected "
                  f"{tuple(EXPECTED_SHAPE)}")

        np.save(OUTPUT_DIR / f"domain_{dom}.npy", out)

        n_filled = -1
        survey_tif = INPUT_DIR / f"resampled_domain_{dom}_survey.tif"
        if survey_tif.exists():
            with rasterio.open(survey_tif) as src:
                survey = src.read(1).astype(np.uint16)
            np.save(SURVEY_DIR / f"domain_{dom}.npy", survey)
            n_filled = int((survey == SURVEY_FILL).sum())

        n_nodata = int((out == NODATA_FILL).sum())
        valid = out[out != NODATA_FILL]
        if len(valid):
            print(f"  domain {dom:>3}: {out.shape}  min={valid.min():6.2f} "
                  f"mean={valid.mean():6.2f} max={valid.max():6.2f}  "
                  f"nodata {n_nodata:5d}  filled {n_filled}")
        else:
            print(f"  domain {dom:>3}: {out.shape}  ALL NODATA")

        audit.append({
            "domain": dom, "rows": out.shape[0], "cols": out.shape[1],
            "nodata_cells": n_nodata, "filled_cells": n_filled,
            "min_m_navd88": float(valid.min()) if len(valid) else None,
            "mean_m_navd88": float(valid.mean()) if len(valid) else None,
            "max_m_navd88": float(valid.max()) if len(valid) else None,
        })

    path = OUTPUT_DIR / AUDIT_CSV
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(audit[0].keys()))
        w.writeheader(); w.writerows(audit)

    print(f"\n{len(audit)} array(s) written")
    if n_bad:
        print(f"*** {n_bad} domain(s) have the wrong shape ***")
    print(f"""
TO USE THIS SET, in HAT_dune_topo_extractor.py:
    DEM_NAME = "{DEM_NAME}"
    VERSION  = <a new version, e.g. "v5">     # picks are per-version
Your 2009_v4 picks stay valid against the un-filled 2009_pea_hatteras set, so
you can run both and compare the road-drowning outcome.""")


if __name__ == "__main__":
    main()
