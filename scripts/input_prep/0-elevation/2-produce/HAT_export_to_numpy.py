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

  LOAD_PATH              INIT_ROOT/1-barrier3d-domains/{TOPO_PRODUCT}/
                         npy-arrays
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

ONE FOLDER PER START PERIOD
----------------------------
Output goes to 1-barrier3d-domains/<TOPO_TARGET>/npy-arrays/, where TOPO_TARGET
is the period the arrays are for - "1984-start" or "2004-start". The two periods
start from different DEMs:

    --product 2009-2014-1996  ->  1984-start
    --product 2009-2014       ->  2004-start

Before 2026-08-25 the tree was keyed on the DEM year and both periods read one
set of arrays. Dune picks are keyed per version WITHIN a product, so a new
version starts from defaults rather than inheriting another version's windows.

THE SURVEY ARRAYS GO IN A SIBLING FOLDER, DELIBERATELY
-------------------------------------------------------
The extractor globs domain_*.npy. A survey array named domain_5_survey.npy would
match that glob and be loaded as if it were a domain. So the survey arrays go to
a separate directory that the extractor never looks at, under the same
domain_<N>.npy name.

INPUTS  (data/hatteras_init/0-elevation/<PRODUCT>/2-resampled-10m/)
    resampled_domain_<N>_filled.tif
    resampled_domain_<N>_survey.tif

OUTPUTS
    data/hatteras_init/1-barrier3d-domains/<TOPO_TARGET>/
        npy-arrays/domain_<N>.npy          m NAVD88, -10 nodata
        npy-arrays_survey/domain_<N>.npy   provenance codes

    python HAT_export_to_numpy.py --product 2009-2014-1996   # -> 1984-start
    python HAT_export_to_numpy.py                            # -> 2004-start

Requires: rasterio, numpy
"""

import csv
import re
import sys
from pathlib import Path

import numpy as np
import rasterio

# =============================================================================
# CONFIG
# =============================================================================

def _find_project_root(start: Path) -> Path:
    """
    Walk up until a directory holds data/hatteras_init.

    NOT parents[N]. This file moved into 2-produce/ on 2026-08-25, and the
    old parents[3] then resolved to input_prep/ rather than the project root.
    That raises nothing - it just makes every path below it wrong, silently,
    until some glob comes back empty. Same helper and same reason as
    4-mgmt-forcings/road_offset/2-audit/HAT_road_setback_audit.py.
    """
    for p in [start, *start.parents]:
        if (p / "data" / "hatteras_init").is_dir():
            return p
    raise SystemExit(f"cannot find data/hatteras_init above {start}")


PROJECT_ROOT = _find_project_root(Path(__file__).resolve())
INIT_ROOT = PROJECT_ROOT / "data" / "hatteras_init"
ELEVATION_DIR = INIT_ROOT / "0-elevation"

# Which elevation product to export. "2009-2014" is the baseline;
# "2009-2014-1996" is the 1984-start DEM. Resolved through
# scripts/hat_elevation_products.py so a layout change cannot leave this
# pointing at a directory that is no longer there.
SOURCE_TAG = "2009-2014"
if "--product" in sys.argv:
    SOURCE_TAG = sys.argv[sys.argv.index("--product") + 1]

sys.path.insert(0, str(PROJECT_ROOT / "scripts"))
from hat_elevation_products import product as _product  # noqa: E402

INPUT_DIR = _product(SOURCE_TAG).resampled_10m
INPUT_GLOB = "resampled_domain_*_filled.tif"
ID_PATTERN = re.compile(r"resampled_domain_(\w+)_filled\.tif$")

DEM_YEAR = "2009"
# WHICH PERIOD PRODUCT these arrays are for. Must match TOPO_PRODUCT in
# HAT_dune_topo_extractor.py - that is the script that reads them.
#     --product 2009-2014-1996  ->  TOPO_TARGET "1984-start"
#     --product 2009-2014       ->  TOPO_TARGET "2004-start"
TOPO_TARGET = "1984-start" if SOURCE_TAG == "2009-2014-1996" else "2004-start"
if "--target" in sys.argv:
    TOPO_TARGET = sys.argv[sys.argv.index("--target") + 1]

DEM_NAME = SOURCE_TAG   # recorded in the audit; no longer a path segment

# Period-first (2026-08-25). Was {DEM_YEAR}-raw/{DEM_YEAR}-npy-arrays/{DEM_NAME}.
NPY_ROOT = INIT_ROOT / "1-barrier3d-domains" / TOPO_TARGET
OUTPUT_DIR = NPY_ROOT / "npy-arrays"
SURVEY_DIR = NPY_ROOT / "npy-arrays_survey"

NODATA_FILL = -10.0   # the extractor detects nodata as raw <= -9.0
EXPECTED_SHAPE = (50, 200)

AUDIT_CSV = "export_audit.csv"

SURVEY_2009, SURVEY_NONE = 2009, 0
# Every non-base code this product's survey rasters may carry. Taken from the
# resolver rather than hardcoded to 2014: the 1984 product also carries 1996,
# and a hardcoded 2014 reported its fill count as if the 1996 graft were not
# there.
from hat_elevation_products import fill_codes as _fill_codes  # noqa: E402
SURVEY_FILL_CODES = list(_fill_codes(SOURCE_TAG)) or [2014]


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
            n_filled = int(np.isin(survey, SURVEY_FILL_CODES).sum())

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
    TOPO_PRODUCT = "{TOPO_TARGET}"
    VERSION      = <a new version, e.g. "v1">   # picks are PER VERSION
Picks live at 1-barrier3d-domains/{TOPO_TARGET}/picks/ and save_windows() writes
back after every domain, so a VERSION whose pick file does not exist starts from
defaults - it does not inherit another version's windows.""")


if __name__ == "__main__":
    main()
