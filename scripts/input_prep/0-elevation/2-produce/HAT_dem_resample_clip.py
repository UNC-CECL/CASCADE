"""
Step 2 of 3: resample each 1 m domain clip to the 10 m Barrier3D grid (50 x 200 per domain).

    python scripts/input_prep/0-elevation/2-produce/HAT_dem_resample_clip.py
    python scripts/input_prep/0-elevation/2-produce/HAT_dem_resample_clip.py --product 2009-2014-1996

Reads a product's 1-gapfill-1m/ clips and writes 2-resampled-10m/ rasters and
resample_audit.csv; next is HAT_export_to_numpy.py. Details: scripts/input_prep/0-elevation/README.md.

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-09-18
"""

import csv
import os
import re
import sys
from pathlib import Path

import numpy as np
import rasterio
from rasterio.transform import Affine


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
ELEVATION_DIR = _el.ELEVATION_ROOT

# --- CONFIG ------------------------------------------------------------------
# Which product to resample: matches the producing script's tag; chosen with --product
DEFAULT_SOURCE = "2009-2014"

SOURCE_TAG = DEFAULT_SOURCE
for _flag in ("--product", "--source"):      # --source kept as an alias
    if _flag in sys.argv:
        SOURCE_TAG = sys.argv[sys.argv.index(_flag) + 1]
        break

sys.path.insert(0, str(PROJECT_ROOT / "scripts"))
from site_layer.hat_elevation_products import product as _product, fill_codes  # noqa: E402

_P = _product(SOURCE_TAG)
INPUT_DIR = _P.gapfill_1m
INPUT_GLOB = "clip_domain_*_filled.tif"
ID_PATTERN = re.compile(r"clip_domain_(\w+)_filled\.tif$")

OUTPUT_DIR = _P.resampled_10m
AUDIT_CSV = "resample_audit.csv"

GRID_SIZE_M = 10.0
EXPECTED_SHAPE = (50, 200)   # rows, cols at 10 m; None to skip the check

# Aggregation: arcgis_bilinear reproduces the existing domains exactly (options in README)
AGGREGATION = "arcgis_bilinear"

# ArcGIS emitted a partial-weight value at some nodata edges and nodata at others
BILINEAR_REQUIRE_ALL_FOUR = True

NODATA_OUT = -9999.0
SURVEY_2009, SURVEY_NONE = 2009, 0
SURVEY_NODATA = 65535

# Every non-base code a survey raster may carry, MOST SPECIFIC FIRST
SURVEY_FILL_CODES = list(fill_codes(SOURCE_TAG))

# Kept because the audit and the console line report "filled cells" against one code
SURVEY_FILL = 2014
# -----------------------------------------------------------------------------


# Resampling - exact block reduction

# Reduces a (block*R, block*C) array to (R, C)
def downsample(arr, block, method):
    h, w = arr.shape
    if h % block or w % block:
        raise ValueError(f"clip {arr.shape} is not a whole multiple of {block}")
    b = arr.reshape(h // block, block, w // block, block)

    if method == "mean":
        return np.nanmean(b, axis=(1, 3)), 0
    if method == "nearest":
        i = block // 2
        return b[:, i, :, i], 0
    if method == "arcgis_bilinear":
        if block % 2:
            raise ValueError("arcgis_bilinear assumes an even block factor")
        lo, hi = block // 2 - 1, block // 2
        core = np.stack([b[:, lo, :, lo], b[:, lo, :, hi],
                         b[:, hi, :, lo], b[:, hi, :, hi]])
        n_valid = np.sum(~np.isnan(core), axis=0)
        partial = int(((n_valid > 0) & (n_valid < 4)).sum())
        if BILINEAR_REQUIRE_ALL_FOUR:
            out = np.where(n_valid == 4, np.nanmean(core, axis=0), np.nan)
        else:
            out = np.nanmean(core, axis=0)
        return out, partial
    raise ValueError(f"unknown AGGREGATION: {method!r}")


# Survey year for the same four cells bilinear reads, so the flag describes the value written
def downsample_survey(survey, block):
    h, w = survey.shape
    b = survey.reshape(h // block, block, w // block, block)
    lo, hi = block // 2 - 1, block // 2
    core = np.stack([b[:, lo, :, lo], b[:, lo, :, hi],
                     b[:, hi, :, lo], b[:, hi, :, hi]])
    out = np.full(core.shape[1:], SURVEY_2009, np.uint16)
    # Reverse precedence so the most specific code is assigned LAST and wins.
    present = [(core == c).any(axis=0) for c in SURVEY_FILL_CODES]
    for c, m in zip(reversed(SURVEY_FILL_CODES), reversed(present)):
        out[m] = c
    n_mixed = int((sum(m.astype(np.int8) for m in present) > 1).sum())         if len(present) > 1 else 0
    out[(core == SURVEY_NONE).any(axis=0)] = SURVEY_NONE
    return out, n_mixed


# One band, with its transform, CRS and nodata
def read_raster(path):
    with rasterio.open(path) as s:
        arr = s.read(1)
        nd = s.nodata
        return arr, s.transform, s.crs, nd


# Write one single-band GeoTIFF
def write_raster(arr, transform, crs, path, dtype, nodata):
    profile = {"driver": "GTiff", "height": arr.shape[0], "width": arr.shape[1],
               "count": 1, "dtype": dtype, "crs": crs, "transform": transform,
               "nodata": nodata, "compress": "deflate"}
    with rasterio.open(path, "w", **profile) as dst:
        dst.write(arr.astype(dtype), 1)


# Run: per domain, resample elevation and provenance, then the audit
def main():
    if not INPUT_DIR.exists():
        raise FileNotFoundError(
            f"{INPUT_DIR} not found - run HAT_dem_gap_fill.py first.")
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    paths = sorted(INPUT_DIR.glob(INPUT_GLOB))
    if not paths:
        raise FileNotFoundError(f"no {INPUT_GLOB} in {INPUT_DIR}")
    print(f"{len(paths)} domain clip(s) in {INPUT_DIR}")
    print(f"Aggregation: {AGGREGATION}")

    audit, n_bad = [], 0
    for p in paths:
        m = ID_PATTERN.search(p.name)
        if not m:
            print(f"  WARNING: cannot parse domain id from {p.name} - skipping")
            continue
        dom = m.group(1)

        arr, t, crs, nd = read_raster(p)
        arr = arr.astype(np.float64)
        if nd is not None and not np.isnan(nd):
            arr = np.where(arr == nd, np.nan, arr)

        block = int(round(GRID_SIZE_M / t.a))
        out, n_partial = downsample(arr, block, AGGREGATION)
        t_out = Affine(GRID_SIZE_M, 0, t.c, 0, -GRID_SIZE_M, t.f)

        ok = EXPECTED_SHAPE is None or out.shape == tuple(EXPECTED_SHAPE)
        if not ok:
            n_bad += 1
            print(f"  domain {dom}: ERROR {out.shape}, expected "
                  f"{tuple(EXPECTED_SHAPE)} - Barrier3D will reject this")

        survey_path = INPUT_DIR / f"clip_domain_{dom}_survey.tif"
        n_filled_10m = -1
        per_code, n_mixed = {}, 0
        if survey_path.exists():
            survey, _, _, _ = read_raster(survey_path)
            survey_out, n_mixed = downsample_survey(survey, block)
            write_raster(survey_out, t_out, crs,
                         OUTPUT_DIR / f"resampled_domain_{dom}_survey.tif",
                         "uint16", SURVEY_NODATA)
            per_code = {c: int((survey_out == c).sum())
                        for c in SURVEY_FILL_CODES}
            n_filled_10m = sum(per_code.values())

        write_raster(np.where(np.isnan(out), NODATA_OUT, out), t_out, crs,
                     OUTPUT_DIR / f"resampled_domain_{dom}_filled.tif",
                     "float32", NODATA_OUT)

        valid = out[~np.isnan(out)]
        if len(valid):
            print(f"  domain {dom:>3}: {out.shape}  min={valid.min():6.2f} "
                  f"mean={valid.mean():6.2f} max={valid.max():6.2f}  "
                  f"nodata {100 * (1 - len(valid) / out.size):5.1f}%  "
                  f"filled {', '.join(f'{c}:{per_code.get(c, 0)}' for c in SURVEY_FILL_CODES)}"
                  f"{f'  mixed {n_mixed}' if n_mixed else ''}")
        else:
            print(f"  domain {dom:>3}: {out.shape}  ALL NODATA")

        audit.append({
            "domain": dom, "rows": out.shape[0], "cols": out.shape[1],
            "shape_ok": ok, "origin_x": t.c, "origin_y": t.f,
            "nodata_frac": round(1 - len(valid) / out.size, 4),
            "filled_cells_10m": n_filled_10m,
            **{f"cells_10m_{c}": per_code.get(c, 0) for c in SURVEY_FILL_CODES},
            "mixed_source_cells": n_mixed,
            "partial_edge_cells": n_partial,
            "min_m": float(valid.min()) if len(valid) else None,
            "mean_m": float(valid.mean()) if len(valid) else None,
            "max_m": float(valid.max()) if len(valid) else None,
        })

    path = OUTPUT_DIR / AUDIT_CSV
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(audit[0].keys()))
        w.writeheader(); w.writerows(audit)

    print(f"\n{len(audit)} domain(s) -> {OUTPUT_DIR}")
    print(f"audit: {path}")
    if n_bad:
        print(f"*** {n_bad} domain(s) are NOT {tuple(EXPECTED_SHAPE)} - see "
              f"shape_ok before running step 3. ***")
    print("\nNext: HAT_export_to_numpy.py writes the .npy arrays the "
          "dune/topo extractor reads.")


if __name__ == "__main__":
    main()
