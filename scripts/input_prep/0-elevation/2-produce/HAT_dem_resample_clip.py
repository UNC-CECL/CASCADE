"""
HAT_dem_resample_clip.py

Step 2 of 3: resamples each gap-filled 1 m domain clip (output of
HAT_dem_gap_fill.py) to the 10 m Barrier3D grid, 50 x 200 per domain.
HAT_export_to_numpy.py converts these to .npy in the final step.

CLIPPING NOW HAPPENS IN STEP 1 - AND THE ORDER MATTERS
-------------------------------------------------------
Barrier3D needs each domain to be an exact 50 x 200 array of true 10 x 10 m
cells. That only works if the 10 m grid is built from the domain's own corner,
so each 10 m cell is an exact 10 x 10 block of 1 m source cells.

Resampling the whole island first and clipping second cannot deliver that. The
domain boxes do not sit on a global 10 m grid - their origins land at arbitrary
sub-10 m offsets (450439.120, 454507.796, ...). Cutting them from a grid whose
cell edges fall on multiples of 10 snaps outward to the enclosing cells, giving
51 x 201 for most domains and shifting every domain up to half a cell off its
polygon.

Confirmed against the existing ArcGIS outputs: clip_domain_N.tif and
resampled_domain_N.tif share a byte-identical origin for every domain checked,
so the 10 m grid was built inside the clip. There is also no global grid worth
preserving - domains 50 and 51 share an x origin but their y origins differ by
505.126 m against a 500 m extent, so the boxes were never a tiling of one grid.

Clipping therefore lives in step 1, which needs the 1 m window anyway to do the
fill. This step only reduces 1 m -> 10 m, which keeps the window definition in
exactly one place.

RESAMPLE METHOD - SETTLED EMPIRICALLY
--------------------------------------
Reconstructed resampled_domain_N.tif from clip_domain_N.tif for domains 1, 2, 3,
25, 50, 75, 100 and 120. ArcGIS used BILINEAR, not the Nearest Neighbor default
the old docstring assumed:

    method                          agreement with existing files
    nearest (any of the 100 cells)  0.1% of cells (chance)
    block mean (all 100 cells)      max diff 1.16 m
    bilinear                        exact, max diff 2.4e-07 (float32 epsilon)

At an exact 10x reduction the 10 m cell center lands on the boundary between
source cells 4 and 5 on both axes, so bilinear collapses to a tie: the plain
mean of the central 2 x 2 source cells. It uses 4 of the 100 cells under each
output cell and ignores the other 96. That reproduces your existing domains, so
it is the default - see AGGREGATION for the alternative.

INPUTS  (data/hatteras_init/0-elevation/1-gapfill-1m/)
    clip_domain_<N>_filled.tif
    clip_domain_<N>_survey.tif

OUTPUTS (data/hatteras_init/0-elevation/2-resampled-10m/)
    resampled_domain_<N>_filled.tif     the 50 x 200 Barrier3D grid
    resampled_domain_<N>_survey.tif     2009 / fill year / 0 per cell
    resample_audit.csv

Requires: rasterio, numpy
"""

import csv
import os
import re
import sys
from pathlib import Path

import numpy as np
import rasterio
from rasterio.transform import Affine

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
ELEVATION_DIR = PROJECT_ROOT / "data" / "hatteras_init" / "0-elevation"

# Must match FILL_SOURCE_TAG in HAT_dem_gap_fill.py, or PRODUCT_TAG in
# HAT_dem_1984_mosaic.py - each source keeps its own subfolder so a re-run
# cannot clobber another source's rasters or its audit CSV.
#
# Selected on the command line rather than by editing, because there are now
# two live products and hand-editing a constant per run is how the figures
# ended up labelled with the wrong source once already:
#
#     python HAT_dem_resample_clip.py                              # baseline
#     python HAT_dem_resample_clip.py --product 2009-2014-1996     # 1984 start
DEFAULT_SOURCE = "2009-2014"

SOURCE_TAG = DEFAULT_SOURCE
for _flag in ("--product", "--source"):      # --source kept as an alias
    if _flag in sys.argv:
        SOURCE_TAG = sys.argv[sys.argv.index(_flag) + 1]
        break

sys.path.insert(0, str(PROJECT_ROOT / "scripts"))
from hat_elevation_products import product as _product, fill_codes  # noqa: E402

_P = _product(SOURCE_TAG)
INPUT_DIR = _P.gapfill_1m
INPUT_GLOB = "clip_domain_*_filled.tif"
ID_PATTERN = re.compile(r"clip_domain_(\w+)_filled\.tif$")

OUTPUT_DIR = _P.resampled_10m
AUDIT_CSV = "resample_audit.csv"

GRID_SIZE_M = 10.0
EXPECTED_SHAPE = (50, 200)   # rows, cols at 10 m; None to skip the check

# AGGREGATION
#   "arcgis_bilinear"  mean of the central 2 x 2 source cells - reproduces your
#                      existing domains exactly. Uses 4 of 100 cells.
#   "mean"             mean of all 100 cells. More defensible for elevation, but
#                      will NOT reproduce the existing files (up to ~1.2 m at
#                      dune crests, where sampling 4 cells is least
#                      representative).
#   "nearest"          single source cell. NOT what produced the existing files.
AGGREGATION = "arcgis_bilinear"

# ArcGIS emitted a partial-weight value at some nodata edges and nodata at
# others. Strict (all 4 required) matched it exactly in the interior and
# differed only at <= 13 cells per domain, all on nodata margins. Strict does
# not invent elevation at the water edge, so it is the default.
BILINEAR_REQUIRE_ALL_FOUR = True

NODATA_OUT = -9999.0
SURVEY_2009, SURVEY_NONE = 2009, 0
SURVEY_NODATA = 65535

# Every non-base code a survey raster may carry, MOST SPECIFIC FIRST. This is
# the precedence downsample_survey resolves a mixed 2 x 2 block with, so the
# order is a decision, not a list:
#
#   1996 first  it is the override that defines the 1984-start product. A 10 m
#               cell that drew any of its four read cells from ALACE should say
#               so - that is the flag a reader uses to find the graft.
#   2014 next   the gap fill, and the only code the 2009-start product has.
#
# SURVEY_NONE outranks both: an unsurveyed cell in the core is also when the
# elevation output is nodata under BILINEAR_REQUIRE_ALL_FOUR, so the two agree
# by construction. Blocks that mix two fill codes are counted per domain in the
# audit as `mixed_source_cells`, so the precedence never hides how often it had
# to choose.
SURVEY_FILL_CODES = list(fill_codes(SOURCE_TAG))

# Kept because the audit and the console line report "filled cells" against one
# code. For the 1984 product that is the 1996 count; the 2014 count is reported
# beside it.
SURVEY_FILL = 2014


# =============================================================================
# RESAMPLING - exact block reduction
# =============================================================================

def downsample(arr, block, method):
    """
    Reduces a (block*R, block*C) array to (R, C). Every output cell is an exact
    block x block window of source cells - that is what makes the 10 m cells
    true 10 x 10 m cells rather than resampled approximations.
    """
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


def downsample_survey(survey, block):
    """
    Survey year for the SAME four cells bilinear actually reads, so the flag
    describes the value that was written rather than the whole block. A 10 m
    cell reads a fill code if any of the central 2 x 2 came from that fill, and
    0 if any of them was unsurveyed - which is also when the elevation output is
    nodata under BILINEAR_REQUIRE_ALL_FOUR, so the two agree by construction.

    Returns (survey_10m, n_mixed) where n_mixed counts blocks whose four read
    cells carried more than one fill code and SURVEY_FILL_CODES had to break
    the tie.
    """
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


def read_raster(path):
    with rasterio.open(path) as s:
        arr = s.read(1)
        nd = s.nodata
        return arr, s.transform, s.crs, nd


def write_raster(arr, transform, crs, path, dtype, nodata):
    profile = {"driver": "GTiff", "height": arr.shape[0], "width": arr.shape[1],
               "count": 1, "dtype": dtype, "crs": crs, "transform": transform,
               "nodata": nodata, "compress": "deflate"}
    with rasterio.open(path, "w", **profile) as dst:
        dst.write(arr.astype(dtype), 1)


# =============================================================================
# MAIN
# =============================================================================

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
