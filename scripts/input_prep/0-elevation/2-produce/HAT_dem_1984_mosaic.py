"""
Build the 1984-start DEM: 2009 + 2014, with the 1996 ALACE beach and foredune laid over them.

    python scripts/input_prep/0-elevation/2-produce/HAT_dem_1984_mosaic.py

Wherever 1996 has data it wins (1996 > 2009 > 2014), with a split floor, a 12 m
ceiling and island connectivity; there is no road boundary. Writes per-domain
1 m clips, a survey-provenance raster and mosaic_1984_audit.csv under
data/hatteras_init/0-elevation/2009-2014-1996/1-gapfill-1m/. Details: scripts/input_prep/0-elevation/README.md.

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-09-18
"""

from __future__ import annotations

import csv
import sys
from pathlib import Path

import numpy as np
import geopandas as gpd
import rasterio
from rasterio.features import rasterize

# Shared helpers come from the 2009 script rather than being copied
sys.path.insert(0, str(Path(__file__).resolve().parent))
import HAT_dem_gap_fill as gf


# Named for its composition, not its start period; resolved through hat_elevation_products.py
PRODUCT_TAG = "2009-2014-1996"
sys.path.insert(0, str(gf.PROJECT_ROOT / "scripts"))
from site_layer.hat_elevation_products import product as _product  # noqa: E402

# --- CONFIG ------------------------------------------------------------------
OUTPUT_DIR = _product(PRODUCT_TAG, check=False).gapfill_1m
AUDIT_CSV = "mosaic_1984_audit.csv"

# the override source
OVERRIDE_DIR = (gf.GIS_ROOT / "Elevation" / "Polygons" / "1996"
                / "1996_FallEC_J1441002")
OVERRIDE_GLOB = "*_FallEC_*t*.tif"     # the two data tiles, not the clippoly
OVERRIDE_TAG = "1996_NASA_ALACE"
OVERRIDE_YEAR = 1996

# The road lines are diagnostic only since 2026-08-26: read and reported, gating nothing
from site_layer.hat_topo_version import road_line_file, road_line_for_year  # noqa: E402
ROAD_LINES = {y: road_line_file(road_line_for_year(y)) for y in (1984, 2004)}

# The split floor: two floors for two questions (see README)
OVERRIDE_FLOOR_INTO_GAP = gf.FILL_MIN_ELEV_NAVD   # -2.64 m NAVD88
OVERRIDE_FLOOR_TO_REPLACE = gf.MHW_ELEVATION      #  0.36 m NAVD88

# The ceiling: reject ALACE returns above 12 m NAVD88, the uncorrected-return tail; nothing is clamped
APPLY_OVERRIDE_CEILING = True
OVERRIDE_MAX_ELEV_NAVD = 12.0     # m NAVD88

# everything below is inherited so the two products cannot drift apart
REQUIRE_ISLAND_CONNECTION = gf.REQUIRE_ISLAND_CONNECTION
GAP_BRIDGE_M = gf.GAP_BRIDGE_M
CONTEXT_BUFFER_M = gf.CONTEXT_BUFFER_M
MHW_ELEVATION = gf.MHW_ELEVATION
NODATA_OUT = gf.NODATA_OUT
EXPECTED_CLIP = gf.EXPECTED_CLIP
GRID_SIZE_M = gf.GRID_SIZE_M

SURVEY_NONE, SURVEY_2009 = 0, 2009
SURVEY_1996, SURVEY_2014 = OVERRIDE_YEAR, gf.FILL_SOURCE_YEAR

# Reproduces HAT_dune_topo_extractor.py, so the audit's shift column matches the extractor
BEACH_START_THR_M = 0.50    # m MHW, strict '>' - the extractor's BEACH_START_THR_M
WATER_CLAMP_M = -3.0        # m MHW - the extractor's WATER_CLAMP_M
# -----------------------------------------------------------------------------


# The override source

# Several single-band tiles read as one raster on the base grid
class TiledSource:

    def __init__(self, directory, pattern, dst_crs):
        paths = sorted(Path(directory).glob(pattern))
        if not paths:
            raise FileNotFoundError(f"no {pattern} in {directory}")
        self.tiles = [gf.FillSource(p, dst_crs) for p in paths]
        self.paths = paths
        first = self.tiles[0]
        self.epsg, self.res, self.nodata = first.epsg, first.res, first.nodata

    def read_on_grid(self, bounds, shape):
        out = np.full(shape, np.nan, np.float32)
        for t in self.tiles:
            d = t.read_on_grid(bounds, shape)
            take = np.isnan(out) & np.isfinite(d)
            out[take] = d[take]
        return out

    def close(self):
        for t in self.tiles:
            t.close()


# The ocean-side boundary

# True ocean-side (east) of the road, per raster row
def ocean_side_mask(road_geom, shape, transform, min_rows=10):
    rl = rasterize([(road_geom, 1)], out_shape=shape, transform=transform,
                   fill=0, all_touched=True).astype(bool)

    col = np.full(shape[0], np.nan)
    for r in range(shape[0]):
        cs = np.where(rl[r])[0]
        if cs.size:
            col[r] = cs.max()

    have = np.isfinite(col)
    n_have = int(have.sum())
    if n_have < min_rows:
        return np.zeros(shape, bool), n_have

    idx = np.arange(shape[0])
    col = np.interp(idx, idx[have], col[have])   # flat-held past both ends

    mask = np.zeros(shape, bool)
    for r in range(shape[0]):
        c = int(round(col[r])) + 1
        if c < shape[1]:
            mask[r, c:] = True
    return mask, n_have


# The diagnostic the window origin depends on

# Median of the extractor's start_beach over the profiles, in m from the window's ocean edge
def start_beach_median_m(arr):
    z = arr[:, ::-1] - MHW_ELEVATION
    z = np.where(np.isnan(z), WATER_CLAMP_M, z)
    z[z < WATER_CLAMP_M] = WATER_CLAMP_M
    above = z > BEACH_START_THR_M
    sb = np.where(above.any(axis=1), above.argmax(axis=1), -1)
    v = sb[sb >= 0]
    return float(np.median(v)) if v.size else float("nan")


# One fill stage

# Connectivity, shared by both stages
def select(base_valid, cand, res_x):
    if not REQUIRE_ISLAND_CONNECTION:
        return cand, 0
    bridge_px = int(round(GAP_BRIDGE_M / res_x / 2.0))
    keep = gf.island_connected(base_valid, cand, bridge_px)
    return keep, int(cand.sum()) - int(keep.sum())


# Run: per domain, build the mosaic, write the clip and provenance, then the audit
def main():
    for label_, path in ([("base DEM", gf.BASE_DEM_PATH),
                          ("domain file", gf.DOMAIN_FILE),
                          ("2014 fill", gf.FILL_DEM_PATH),
                          ("1996 override", OVERRIDE_DIR)]
                         + [(f"{y} road line (diagnostic)", q)
                            for y, q in ROAD_LINES.items()]):
        if not path.exists():
            raise FileNotFoundError(
                f"{label_} not found: {path}\n"
                f"  Fix the paths at the top of this script.")

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    src = rasterio.open(gf.BASE_DEM_PATH)
    res_x, res_y = src.transform.a, -src.transform.e
    block = int(round(GRID_SIZE_M / res_x))
    pad_px = int(round(CONTEXT_BUFFER_M / res_x))
    print(f"\n2009 DEM: {src.height} x {src.width}, {res_x} m, "
          f"nodata={src.nodata}")

    gdf = gpd.read_file(gf.DOMAIN_FILE)
    gdf = gf.resolve_crs(src, gdf)
    print(f"{len(gdf)} domains from {gf.DOMAIN_FILE.name}")

    dem_crs = gf.CRS.from_wkt(src.crs.to_wkt())
    hz = dem_crs.sub_crs_list[0] if dem_crs.is_compound else dem_crs
    road_geoms = {}
    for yr, q in ROAD_LINES.items():
        rl = gpd.read_file(q)
        print(f"{yr} NC-12: {q.name}, {rl.crs} -> EPSG:{hz.to_epsg()}"
              f"   DIAGNOSTIC ONLY - gates nothing")
        road_geoms[yr] = rl.to_crs(hz).union_all()

    fill14 = gf.FillSource(gf.FILL_DEM_PATH, src.crs)
    print(f"gap fill : {gf.FILL_DEM_PATH.name}  ({gf.FILL_SOURCE_TAG})  "
          f"EPSG:{fill14.epsg}  {fill14.res} m")

    ov96 = TiledSource(OVERRIDE_DIR, OVERRIDE_GLOB, src.crs)
    print(f"override : {OVERRIDE_DIR.name}  ({OVERRIDE_TAG})  "
          f"EPSG:{ov96.epsg}  {ov96.res} m  nodata={ov96.nodata}")
    for p in ov96.paths:
        print(f"             {p.name}")

    print(f"\nPriority : {OVERRIDE_YEAR} > 2009 > {SURVEY_2014}, EVERYWHERE 1996 has data")
    print(f"           no road boundary - the landward limit is the ALACE "
          f"swath edge itself, not NC-12; the road lines are reported, "
          f"not applied")
    _ceil = (f"{OVERRIDE_MAX_ELEV_NAVD:+.2f} m NAVD88"
             if APPLY_OVERRIDE_CEILING else "OFF")
    print(f"Ceiling  : 1996 {_ceil}"
          f"   - the ALACE return tail reaches +256.65 m")
    print(f"Floors   : 1996 into gap {OVERRIDE_FLOOR_INTO_GAP:+.2f} m NAVD88   "
          f"1996 over measured {OVERRIDE_FLOOR_TO_REPLACE:+.2f} m NAVD88")
    print(f"           2014 into gap {gf.FILL_MIN_ELEV_NAVD:+.2f} m NAVD88")
    print(f"Vertical : bias correction OFF, feathering OFF - a written cell is "
          f"its own survey's measurement\n")

    audit = []
    for _, row in gdf.iterrows():
        dom = row[gf.DOMAIN_ID_FIELD]
        try:
            dom = int(dom)
        except (TypeError, ValueError):
            pass

        win, adj = gf.snap_window(row.geometry.bounds, src.transform,
                                  res_x, res_y, block)
        big = gf.pad_window(win, pad_px)
        t_big = src.window_transform(big)

        base = gf.read_window(src, big, src.nodata)
        bx0, by1 = t_big.c, t_big.f
        bx1 = bx0 + big.width * res_x
        by0 = by1 - big.height * res_y
        bounds = (bx0, by0, bx1, by1)

        g96 = ov96.read_on_grid(bounds, base.shape)
        g14 = fill14.read_on_grid(bounds, base.shape)

        valid09 = np.isfinite(base)
        has96, has14 = np.isfinite(g96), np.isfinite(g14)

        # Built for the AUDIT ONLY - neither is applied
        oceans = {yr: ocean_side_mask(g, base.shape, t_big)
                  for yr, g in road_geoms.items()}
        ocean, n_road_rows = oceans[1984]
        ocean04 = oceans[2004][0]
        has_road = n_road_rows > 0

        # A gap means no other survey, 2009 or 2014, saw the cell
        covered_by_other = valid09 | has14
        # NO ROAD BOUNDARY - the one substantive change of 2026-08-26
        cand96_cov = has96.copy()
        if APPLY_OVERRIDE_CEILING:
            too_high = cand96_cov & (g96 > OVERRIDE_MAX_ELEV_NAVD)
            n96_ceiling = int(too_high.sum())
            cand96_cov = cand96_cov & ~too_high
        else:
            n96_ceiling = 0
        into_gap = (cand96_cov & ~covered_by_other
                    & (g96 >= OVERRIDE_FLOOR_INTO_GAP))
        replace = (cand96_cov & covered_by_other
                   & (g96 >= OVERRIDE_FLOOR_TO_REPLACE))
        n96_cov = int(cand96_cov.sum())
        n96_floor_gap = (int((cand96_cov & ~covered_by_other).sum())
                         - int(into_gap.sum()))
        n96_floor_rep = (int((cand96_cov & covered_by_other).sum())
                         - int(replace.sum()))

        cand96 = into_gap | replace
        cand96, n96_conn = select(valid09, cand96, res_x)

        merged = base.copy()
        if cand96.any():
            merged[cand96] = g96[cand96]

        # Bias reported every run, applied or not; negative means 1996 sits above 2009
        ov = valid09 & has96
        bias96 = (float(np.nanmedian(base[ov] - g96[ov]))
                  if int(ov.sum()) >= gf.MIN_RING_CELLS else 0.0)
        # The ocean-side-of-1984 figure the docstring's datum argument was built on
        ov84 = ov & ocean
        bias96_o84 = (float(np.nanmedian(base[ov84] - g96[ov84]))
                      if int(ov84.sum()) >= gf.MIN_RING_CELLS else 0.0)
        seam96 = gf.seam_check(merged, valid09 & ~cand96, cand96)

        # Run against the POST-1996 surface
        valid_now = np.isfinite(merged)
        cand14 = ~valid_now & has14 & (g14 >= gf.FILL_MIN_ELEV_NAVD)
        n14_cov = int((~valid_now & has14).sum())
        n14_floor = n14_cov - int(cand14.sum())
        cand14, n14_conn = select(valid_now, cand14, res_x)
        if cand14.any():
            merged[cand14] = g14[cand14]

        # 'Before' rebuilt through the same code path, so a stale product cannot change the answer
        before = base.copy()
        b_cand = ~valid09 & has14 & (g14 >= gf.FILL_MIN_ELEV_NAVD)
        b_cand, _ = select(valid09, b_cand, res_x)
        if b_cand.any():
            before[b_cand] = g14[b_cand]

        r0 = c0 = pad_px
        crop = np.s_[r0:r0 + win.height, c0:c0 + win.width]
        sb_before = start_beach_median_m(before[crop])
        sb_after = start_beach_median_m(merged[crop])
        sb_shift = sb_before - sb_after       # +ve = moved seaward

        # write
        out = merged[crop]
        gap_out = np.isnan(out)
        survey = np.full(out.shape, SURVEY_2009, np.uint16)
        survey[cand14[crop]] = SURVEY_2014
        survey[cand96[crop]] = SURVEY_1996
        survey[gap_out] = SURVEY_NONE

        if out.shape != EXPECTED_CLIP:
            print(f"  domain {dom}: ERROR clip {out.shape}, "
                  f"expected {EXPECTED_CLIP}")

        t_out = src.window_transform(win)
        gf.write_raster(np.where(gap_out, NODATA_OUT, out), t_out, src.crs,
                        OUTPUT_DIR / f"clip_domain_{dom}_filled.tif",
                        "float32", NODATA_OUT)
        gf.write_raster(survey, t_out, src.crs,
                        OUTPUT_DIR / f"clip_domain_{dom}_survey.tif",
                        "uint16", gf.SURVEY_NODATA)

        n96_out = int(cand96[crop].sum())
        n96_new = int((cand96 & ~covered_by_other)[crop].sum())
        n96_rep = int((cand96 & covered_by_other)[crop].sum())
        n14_out = int(cand14[crop].sum())
        gap_before = int(np.isnan(base[crop]).sum())

        print(f"  domain {dom:>3}: {'road' if has_road else ' -- '}  "
              f"1996 {n96_out:6d} (new {n96_new:6d} / over {n96_rep:6d})  "
              f"2014 {n14_out:6d}  nodata {gap_before:6d} -> "
              f"{int(gap_out.sum()):6d}  "
              f"start_beach {sb_before:6.1f} -> {sb_after:6.1f} m "
              f"({sb_shift:+6.1f} m, {sb_shift / 10:+.1f} cell)  "
              f"bias96 {bias96:+.3f}", flush=True)

        audit.append({
            "domain": dom, "rows": out.shape[0], "cols": out.shape[1],
            # DIAGNOSTIC ONLY from 2026-08-26 - neither alignment gates the override any more
            "road_line": has_road, "road_rows": n_road_rows,
            "ocean_cells": int(ocean[crop].sum()),
            "road_rows_2004": int(oceans[2004][1]),
            "ocean_cells_2004": int(ocean04[crop].sum()),
            # The numbers that justify dropping the boundary, per domain.
            "n1996_landward_of_1984_road": int((cand96 & ~ocean)[crop].sum()),
            "n1996_landward_of_2004_road": int((cand96 & ~ocean04)[crop].sum()),

            "start_beach_before_m": round(sb_before, 1),
            "start_beach_after_m": round(sb_after, 1),
            "start_beach_shift_m": round(sb_shift, 1),
            "start_beach_shift_cells10": round(sb_shift / 10.0, 2),

            "n1996_written": n96_out,
            "n1996_new_land": n96_new,
            "n1996_overwrote_2009": n96_rep,
            "cand1996_coverage": n96_cov,
            "dropped_1996_ceiling": n96_ceiling,
            "dropped_1996_floor_gap": n96_floor_gap,
            "dropped_1996_floor_replace": n96_floor_rep,
            "dropped_1996_connectivity": n96_conn,
            "bias_2009_minus_1996_m": round(bias96, 4),
            "bias_2009_minus_1996_cells": int(ov.sum()),
            "bias_2009_minus_1996_ocean1984_m": round(bias96_o84, 4),
            "bias_2009_minus_1996_ocean1984_cells": int(ov84.sum()),

            "n2014_written": n14_out,
            "cand2014_coverage": n14_cov,
            "dropped_2014_floor": n14_floor,
            "dropped_2014_connectivity": n14_conn,

            "max_elev_m": (round(float(np.nanmax(out)), 2)
                           if np.isfinite(out).any()
                           else float("nan")),
            "nodata_before": gap_before,
            "nodata_after": int(gap_out.sum()),
            **{f"seam96_{k}": (round(v, 4) if isinstance(v, float) else v)
               for k, v in seam96.items()},
            "snap_x_m": round(adj["snap_x_m"], 4),
            "snap_y_m": round(adj["snap_y_m"], 4),
            "trim_x_m": adj["trim_x_m"], "trim_y_m": adj["trim_y_m"],
        })

    src.close()
    fill14.close()
    ov96.close()

    path = OUTPUT_DIR / AUDIT_CSV
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(audit[0].keys()))
        w.writeheader()
        w.writerows(audit)

    with_road = [a for a in audit if a["road_line"]]
    # EVERY domain now receives 1996, so the shift statistics run over all of them
    shifts = np.array([a["start_beach_shift_m"] for a in audit])
    print(f"\n{'=' * 78}\n{len(audit)} domains  "
          f"(1996 admitted in ALL of them - no road boundary; "
          f"{len(with_road)} carry a 1984 road line for the diagnostic, "
          f"{len(audit) - len(with_road)} do not)")
    print(f"  1996 written   {sum(a['n1996_written'] for a in audit):,}  "
          f"(new land {sum(a['n1996_new_land'] for a in audit):,}, "
          f"over 2009 {sum(a['n1996_overwrote_2009'] for a in audit):,})")
    print(f"  2014 written   {sum(a['n2014_written'] for a in audit):,}")
    print(f"  1996 rejected  ceiling "
          f"{sum(a['dropped_1996_ceiling'] for a in audit):,}   "
          f"floor-in-gap "
          f"{sum(a['dropped_1996_floor_gap'] for a in audit):,}   "
          f"floor-to-replace "
          f"{sum(a['dropped_1996_floor_replace'] for a in audit):,}   "
          f"connectivity "
          f"{sum(a['dropped_1996_connectivity'] for a in audit):,}")
    tb = sum(a["nodata_before"] for a in audit)
    ta = sum(a["nodata_after"] for a in audit)
    print(f"  nodata         {tb:,} -> {ta:,} "
          f"({100 * (tb - ta) / max(tb, 1):.1f}% recovered)")
    lw84 = sum(a["n1996_landward_of_1984_road"] for a in audit)
    lw04 = sum(a["n1996_landward_of_2004_road"] for a in audit)
    tot96 = max(sum(a["n1996_written"] for a in audit), 1)
    print(f"  1996 landward of NC-12   "
          f"1984 line {lw84:,} ({100 * lw84 / tot96:.1f}% of what 1996 wrote)"
          f"   2004 line {lw04:,} ({100 * lw04 / tot96:.1f}%)")
    print(f"    That is what the ocean-side boundary used to discard. It is "
          f"backdune, not sound - the ALACE swath ends on its own well short "
          f"of the bay shore in every domain measured.")
    print(f"\n  start_beach shift, m (+ve = seaward), {len(audit)} domains:")
    print(f"    median {np.median(shifts):+.1f}   "
          f"mean {shifts.mean():+.1f}   "
          f"min {shifts.min():+.1f}   max {shifts.max():+.1f}")
    print(f"    seaward >= 1 cell (10 m)  {int((shifts >= 10).sum())}"
          f"    unchanged  {int((shifts == 0).sum())}"
          f"    LANDWARD  {int((shifts < 0).sum())}")
    # The threshold that matters is one Barrier3D cell, not one metre
    hard = [(a["domain"], a["start_beach_shift_m"])
            for a in audit if a["start_beach_shift_m"] <= -GRID_SIZE_M]
    soft = [(a["domain"], a["start_beach_shift_m"])
            for a in audit if -GRID_SIZE_M < a["start_beach_shift_m"] < 0]
    if soft:
        print(f"    landward, SUB-CELL (< {GRID_SIZE_M:g} m, "
              f"cannot move a 10 m cell by more than one): {soft}")
    if hard:
        print(f"\n    *** LANDWARD BY >= ONE MODEL CELL: {hard}")
        print(f"    *** The split floor is supposed to take these to zero. A "
              f"non-empty list here means 1996 wet-edge returns are still "
              f"displacing measured ground - re-check OVERRIDE_FLOOR_TO_REPLACE "
              f"and the `covered_by_other` test before using this product.")
    peak = max(a["max_elev_m"] for a in audit)
    print(f"\n  highest cell in the product  {peak:.2f} m NAVD88"
          f"   (2009+2014 island-wide max is 10.18 m)")
    if APPLY_OVERRIDE_CEILING and peak > OVERRIDE_MAX_ELEV_NAVD + 0.01:
        print(f"    *** ABOVE THE CEILING. Either it is not being "
              f"applied, or the 2009 / 2014 surface itself carries "
              f"a spike. Find it before using this product.")
    print(f"\n  audit: {path}")
    print(f"\nNext: HAT_dem_resample_clip.py with SOURCE_TAG = "
          f"'{PRODUCT_TAG}', then HAT_plot_gapfill.py --source {PRODUCT_TAG}.")


if __name__ == "__main__":
    main()
