"""
HAT_dem_1984_mosaic.py

Step 1 of 3 for the 1984-START topography. The 2009-start product is made by
HAT_dem_gap_fill.py and is NOT touched by this script; the two write to
different tags and neither can overwrite the other.

WHAT THIS DOES THAT HAT_dem_gap_fill.py DOES NOT
------------------------------------------------
HAT_dem_gap_fill.py only ever writes into NODATA. Its four rules SELECT which
measured cells are allowed into a hole; a cell the 2009 survey saw is never
changed.

This script adds a fifth stage that OVERWRITES measured 2009 ground with the
1996 NOAA/NASA ALACE survey. That is a different operation and it needs its own
justification, which is below.

    wherever 1996 has data :   1996  >  2009  >  2014
    everywhere else        :           2009  >  2014

There is NO road boundary. Until 2026-08-26 the override was confined to the
ocean side of the 1984 NC-12 alignment; it no longer is. The landward limit is
the ALACE swath's own edge. Why, and what it cost, is the next section.

WHY THERE IS NO ROAD BOUNDARY
-----------------------------
The original design used "ocean side of NC-12" as the landward limit, on the
reasoning that a rule you can state is better than a survey artefact. Choosing
WHICH NC-12 - the 1984 or the 2004 alignment - forced the question of what the
boundary was actually doing, and it turned out to be doing very little.

All three options were measured across all 90 domains by running stage 1 three
times with only the mask changed, every other rule identical:

    boundary        1996 cells written      new land      overwrote measured
    1984 line             11,983,935       2,084,834               9,899,101
    2004 line             12,251,459       2,086,884              10,164,575
    none                  14,707,324       2,320,397              12,386,927

Read the middle column. Moving 1984 -> 2004 buys 267,524 cells of which only
2,050 are ground no other survey saw. 99.2% of that gain is 1996 replacing a
2009 or 2014 measurement in the backdune - the part of the ALACE swath with
the WORST coverage. It re-vintages a strip; it recovers almost nothing.

Dropping the boundary buys 2,723,389 cells, 235,563 of them new land - 115x
the recovered ground of the 2004 line - and it is the only option that gives
domains 1-7 any 1996 at all.

The two objections to dropping it both failed on measurement:

  * IT DOES NOT REACH THE SOUND. Cross-shore reach of what 1996 writes, metres
    from the ocean edge, median over profiles, against the island's own extent:

        domain      island      road 84 / 04      1996 reach, no boundary
            11        1999          416 / 499                          560
            60        1274          621 / 621                          979
            77        1627          717 / 717                          429
            85        1782          185 / 305                          432
            90        1999          382 / 382                          441

    ALACE stops itself at 429-979 m against an island 1274-1999 m wide. The
    road was never what held it back. Domain 77 is the proof: the reach is
    429 m under every option while the road sits at 717 m, so there the
    boundary was not even binding.

  * IT ADMITS NO EXTRA JUNK. Island-wide maximum stays 12.00 m under all three
    options - the ceiling, not the mask, is what holds the return tail.

And the window origin is untouched by the choice. start_beach moves in ZERO
domains between the 1984 and 2004 masks (max |delta| 0.00 m) and in 6 domains
between 1984 and none, max 26 m. All three agree at the seaward edge; every
difference is landward, behind the beach start.

WHAT WAS GIVEN UP. A stated rule was traded for a survey footprint. "Ocean
side of NC-12" is defensible in a methods section; "wherever ALACE flew" is a
survey artefact standing in for a landform boundary. The honest description of
this product is that its 1996 extent is the ALACE swath, and the audit reports
n1996_landward_of_1984_road / _2004_road per domain so the discarded band is a
number rather than a claim.

Both alignments are still read and still rasterized every run. They gate
nothing. They are reported.

WHAT THE 1996 SOURCE ACTUALLY IS
--------------------------------
1996 Fall East Coast NOAA/NASA Airborne LiDAR Assessment of Coastal Erosion
(ALACE), tiles J1441002 R0_C0 and R1_C0. 3 m grid, EPSG:6347 (NAD83(2011) /
UTM 18N), nodata -999999, +/-15 cm stated vertical accuracy.

Its abstract is the whole methodological problem in one sentence: the aircraft
surveyed "from the low water line to the landward base of the sand dunes."

Measured consequences, over the ocean-side band in 13 sampled domains:

    coverage of the ocean-side band   30-84%, median ~53%
    coverage LANDWARD of the road     0-16%   - it stops before NC-12
                                      ^ this is now ADMITTED, not discarded,
                                        and is 17.2% of what 1996 writes
    seaward reach vs 2009 waterline   +70 to +133 m
    median(1996 - 2009) on overlap    +0.11 to +0.36 m in the developed reach

So this is NOT "the ocean side becomes 1996". It is a 1996 beach and foredune
grafted onto a 2009 backdune and interior, and the seam lands at the ALACE
swath's landward edge - near the dune toe - not at the road. Which is the same
place HAT_dune_topo_extractor.py picks its dune windows. Read a dune-crest
number near that line knowing it sits beside a survey boundary.

VERTICAL DATUM IS AN INFERENCE, NOT A DECLARATION
-------------------------------------------------
The .prj carries the horizontal CRS only. The sole evidence for the vertical is
`geoid18` in the NOAA Digital Coast S3 path the metadata cites, i.e. NAVD88 via
GEOID18, metres. The empirical check is what actually supports it: median
(1996 - 2009) over co-measured cells is +0.11 to +0.36 m in the developed
reach. A NGVD29 or ellipsoid-height mistake would show METRES of disagreement,
not decimetres, and the sign is the one 12 years of erosion predicts. Both
per-domain bias estimates are written to the audit every run so this stays
checkable rather than asserted.

Domains 11 and 86 report +1.64 and +1.66 m. Both have 4% and 22% 2009 coverage
in the band, so the overlap is a sliver at the waterline. Do not read those two
as a datum finding.

THE SPLIT FLOOR, AND WHY IT IS NOT THE SAME NUMBER TWICE
--------------------------------------------------------
The first build used the existing -2.64 m NAVD88 floor for both cases. It moved
33 of 83 domains' beach start LANDWARD, by up to 27 m, while 1996 measured
HIGHER than 2009 in the overlap - the opposite of what erosion predicts.

The cause is the ALACE swath running to the low water line. It carries wet
swash and beach-face returns between -2.64 m and MHW. Ocean-side of the road
those cells overwrote 2009's dry beach and pushed the first cell above
BEACH_START_THR_M (0.50 m MHW) landward. Measured on the worst offenders,
median shift of start_beach in metres, positive = seaward:

    domain   floor -2.64, overwrite   floor MHW, overwrite   1996 into gaps only
      73             -27.5                    0.0                    0.0
      42             -23.0                    0.0                    0.0
      58             -22.0                    0.0                    0.5
      90             -16.0                   +1.0                   +1.0
      79             +52.5                  +52.5                  +51.5
      80             +71.0                  +71.0                  +71.0
      77             +59.0                  +59.0                  +59.0
      11             +57.0                  +57.0                  +57.0

Two things fall out of that table and both matter:

  * Raising the floor to MHW *for the overwrite case only* takes every landward
    mover to zero and leaves every seaward mover untouched.
  * The seaward signal is identical when 1996 is only allowed into nodata. The
    real gain is where the 2009 survey saw NOTHING, not where the two surveys
    disagree. The overwrite is still kept, because it is what carries the 1996
    beach and foredune ELEVATIONS - it just does not move the window.

Hence two floors, for two different questions:

    NO other survey    ->  admit 1996 at -2.64 m, the extractor's WATER_CLAMP_M,
    saw the cell           exactly as 2014 is admitted. Nothing is lost, and the
                           sub-MHW seaward gain is kept.
    2009 OR 2014 saw   ->  admit 1996 only above MHW. One survey may replace
    the cell               another's measurement, but not with a wet return.

    The first build tested only `~valid09` for the first case and shipped 16
    landward movers, D73 at -21 m and D90 at -16 m. The 2009 holes are not
    empty in the product - 2014 fills them - so that let a -2 m 1996 swash
    return displace a +1 m 2014 measurement. The guard at the end of main()
    is what caught it, and it stays there for the next time.

No value is modified by either. Only the selection differs, and every cell each
floor rejects is counted per domain in the audit.

Resampling 3 m -> 1 m was tested and is irrelevant here: nearest and bilinear
agree to within 1 m of start_beach in every domain checked. Nearest is used, as
in HAT_dem_gap_fill.FillSource, because it invents no values at survey edges.

DOMAINS 1-7 - RESOLVED 2026-08-26
---------------------------------
Both NC-12 exports end at domain 8. Domains 1-7 sit 265 m to 3.3 km beyond the
terminus, so there was no line to define an ocean side for them, and under the
old boundary they received no 1996 at all. They carried a 2009 beach while
8-90 carried a 1996 one, leaving an alongshore discontinuity in the initial
condition at the 7/8 seam. The documented alternative was to invent a boundary
by extending domain 8's alignment south across Hatteras village.

Removing the boundary dissolves the problem instead of trading it. 1-7 now
receive 1996 on exactly the same terms as every other domain - 77k to 175k
cells each, 6k to 60k of it ground no other survey saw - and nothing is
invented. The seam is gone.

They are still flagged `road_line=False` in the audit. That flag is now purely
descriptive: it says these domains have no NC-12 export, not that they were
treated differently.

WHAT MOVES DOWNSTREAM, AND WHAT DOES NOT
-----------------------------------------
MOVES.  Each domain's cross-shore window origin. HAT_dune_topo_extractor.py
        sets start_beach to the first cell above 0.50 m MHW, so adding 1996
        land above that threshold slides the window seaward. Measured with the
        split floor: 35 of 83 domains move >= 1 model cell seaward, 6 are
        unchanged, none move landward. The developed reach 76-81 moves 47-72 m,
        i.e. 5-7 Barrier3D cells. Every domain's before/after is in the audit as
        start_beach_before_m / start_beach_after_m / start_beach_shift_m, so
        this is verifiable per domain rather than asserted island-wide.
        A NEW PICK SET IS THEREFORE REQUIRED - the 2009_v5 dune windows were
        picked against a different origin.

DOES NOT.  `shoreline_offset`. That comes from
        2-brie-offset/1984/Island_Dune_Offsets_1984_CASCADE_Input.csv,
        is measured independently of any DEM, already exists for 1984, and is
        passed straight to Cascade(). Nothing here touches it. It was a
        deliberate choice to leave the two independent rather than re-derive one
        from the other.

INPUTS
    D:/Hatteras_GIS/.../2009_full.tif                base, 1 m, EPSG:3725
    D:/Hatteras_GIS/.../2014_full.tif                gap fill, 1 m, EPSG:6347
    D:/Hatteras_GIS/.../1996_FallEC_J1441002/*.tif   override, 3 m, EPSG:6347
    D:/Hatteras_GIS/domains.geojson                  90 boxes, 2000 x 500 m
    4-mgmt-forcing/.../1978/nc12_1978.geojson        EPSG:2264, US survey FEET

OUTPUTS (data/hatteras_init/0-elevation/2009-2014-1996/1-gapfill-1m/)
    clip_domain_<N>_filled.tif   the mosaic, m NAVD88
    clip_domain_<N>_survey.tif   provenance per cell:
                                     1996 = ALACE, wherever it has data
                                     2009 = measured by the 2009 DEM
                                     2014 = filled from 2014 NOAA Post-Sandy
                                        0 = no survey saw it
    mosaic_1984_audit.csv        per-domain counts for every rule above

Filenames match HAT_dem_gap_fill.py's so HAT_dem_resample_clip.py needs only
its SOURCE_TAG changed. The survey raster now carries FOUR codes rather than
three, which its downsample_survey() has been extended to handle.

Requires: rasterio, geopandas, numpy, scipy

    python HAT_dem_1984_mosaic.py
"""

from __future__ import annotations

import csv
import sys
from pathlib import Path

import numpy as np
import geopandas as gpd
import rasterio
from rasterio.features import rasterize

# Shared helpers come from the 2009 script rather than being copied. Same
# directory, and importing it runs no IO - its module level is paths and
# constants only, and main() is guarded.
sys.path.insert(0, str(Path(__file__).resolve().parent))
import HAT_dem_gap_fill as gf


# =============================================================================
# CONFIG
# =============================================================================

# Named for its COMPOSITION rather than the start period - see the note on
# PRODUCT_TAG in HAT_dem_gap_fill.py. Resolved through
# scripts/site_layer/hat_elevation_products.py so the layout lives in one place.
PRODUCT_TAG = "2009-2014-1996"
sys.path.insert(0, str(gf.PROJECT_ROOT / "scripts"))
from site_layer.hat_elevation_products import product as _product  # noqa: E402

OUTPUT_DIR = _product(PRODUCT_TAG, check=False).gapfill_1m
AUDIT_CSV = "mosaic_1984_audit.csv"

# --- the override source ---
OVERRIDE_DIR = (gf.GIS_ROOT / "Elevation" / "Polygons" / "1996"
                / "1996_FallEC_J1441002")
OVERRIDE_GLOB = "*_FallEC_*t*.tif"     # the two data tiles, not the clippoly
OVERRIDE_TAG = "1996_NASA_ALACE"
OVERRIDE_YEAR = 1996

# --- THE ROAD LINES ARE DIAGNOSTIC ONLY (2026-08-26) ------------------------
#
# They no longer gate anything. 1996 is admitted wherever it has data that
# clears the floors, the ceiling and connectivity, and the landward limit is
# the ALACE swath's own edge. See "WHY THERE IS NO ROAD BOUNDARY" above.
#
# Both vintages are still read and still reported per domain, because "how far
# landward of NC-12 did 1996 actually write?" is the first question a reader
# will put to this product, and the audit should answer it rather than leave
# it to be re-derived. NO_ROAD_POLICY is gone with the boundary it configured.
# Keyed by PERIOD; each period's LINE comes from ROAD_LINE_FOR_YEAR.
from site_layer.hat_topo_version import road_line_file, road_line_for_year  # noqa: E402
ROAD_LINES = {y: road_line_file(road_line_for_year(y)) for y in (1984, 2004)}

# --- THE SPLIT FLOOR. See the module docstring; these are not the same number
#     twice by accident, and collapsing them re-introduces the 33 landward
#     movers. ---
OVERRIDE_FLOOR_INTO_GAP = gf.FILL_MIN_ELEV_NAVD   # -2.64 m NAVD88
OVERRIDE_FLOOR_TO_REPLACE = gf.MHW_ELEVATION      #  0.36 m NAVD88

# --- THE CEILING. 2009 and 2014 are clean gridded DEMs and need none; this is
#     specific to ALACE and it is not optional. ---
#
# The 1996 grid runs from -16.56 to +256.65 m NAVD88. Its 99th percentile is
# 7.61 m, so everything above roughly 12 m is the uncorrected-return tail that
# ALACE-era ATM data is known for - cloud, bird and aircraft returns the
# vendor's "visual inspection" filter did not catch.
#
# Shipped without a ceiling, the first build put 250 m spikes into 66 of 90
# domains, against a 2009+2014 island-wide max of 10.18 m. The floors only ever
# guarded the low side; nothing guarded this one.
#
# 12 m, from the distribution of the 12,021,270 cells 1996 actually wrote,
# 1 m bins:
#
#      5 m 459288   6 m 222927   7 m 86532   8 m 34755   9 m 18645
#     10 m  15561  11 m  13362  12 m 10458  13 m  4881  14 m  2328
#     15 m   1923  16 m   1752  17 m   885  18 m  1224  19 m   675
#     20 m    594  ... a flat few-hundred-per-bin tail to 58 m, then to 256 m
#
# The terrain population decays steeply to 12 m and HALVES at 13 m; past ~15 m
# the histogram is flat, which is an artifact population, not a landform. 12 m
# also sits 18% above the 10.18 m island-wide max of the 2009+2014 product -
# enough headroom for a 1996 dune genuinely taller than anything the 2009
# survey measured, which is the erosion signal this product exists to capture,
# without admitting the tail.
#
# Rejects 37,335 cells, 0.31% of what 1996 wrote. A rejected cell is NOT lost:
# it falls through to 2009, then to 2014, so the ceiling costs coverage only
# where no other survey saw the cell at all. Counted per domain in the audit as
# dropped_1996_ceiling.
#
# NOTHING IS CLAMPED. A cell over the ceiling is rejected, not pulled down to
# it - the same posture as every other rule here.
APPLY_OVERRIDE_CEILING = True
OVERRIDE_MAX_ELEV_NAVD = 12.0     # m NAVD88

# --- everything below is inherited so the two products cannot drift apart ---
REQUIRE_ISLAND_CONNECTION = gf.REQUIRE_ISLAND_CONNECTION
GAP_BRIDGE_M = gf.GAP_BRIDGE_M
CONTEXT_BUFFER_M = gf.CONTEXT_BUFFER_M
MHW_ELEVATION = gf.MHW_ELEVATION
NODATA_OUT = gf.NODATA_OUT
EXPECTED_CLIP = gf.EXPECTED_CLIP
GRID_SIZE_M = gf.GRID_SIZE_M

SURVEY_NONE, SURVEY_2009 = 0, 2009
SURVEY_1996, SURVEY_2014 = OVERRIDE_YEAR, gf.FILL_SOURCE_YEAR

# Reproduces HAT_dune_topo_extractor.py so the audit's shift column means the
# same thing the extractor will do. Changing either without the other makes the
# diagnostic quietly wrong.
BEACH_START_THR_M = 0.50    # m MHW, strict '>' - extractor line 297
WATER_CLAMP_M = -3.0        # m MHW - extractor line 297


# =============================================================================
# THE OVERRIDE SOURCE
# =============================================================================

class TiledSource:
    """
    Several single-band tiles read as one raster on the base grid.

    ALACE ships J1441002 as two row-tiles that abut at y=3947502 rather than as
    one mosaic. They are read independently and pasted in order; a later tile
    never overwrites a cell an earlier one already filled, so the shared edge
    resolves deterministically instead of by read order.

    Each tile goes through gf.FillSource, so the CRS difference (EPSG:6347 ->
    EPSG:3725) and the 3 m -> 1 m resolution difference are handled exactly as
    they are for the 2014 fill, by the same code, with the same nearest-
    neighbour choice.
    """

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


# =============================================================================
# THE OCEAN-SIDE BOUNDARY
# =============================================================================

def ocean_side_mask(road_geom, shape, transform, min_rows=10):
    """
    True ocean-side (east) of the road, per raster row.

    OCEAN_LOC is "right" throughout this pipeline - the domains are 2000 m
    cross-shore in x by 500 m alongshore in y, and the Atlantic is at increasing
    easting. So the boundary is one column index per row.

    The MOST SEAWARD road cell in a row is used, not the first. Where the
    alignment runs diagonally through a row it occupies several columns, and
    taking the max yields the SMALLEST ocean region - the choice that cannot
    accidentally place a road cell on the ocean side of itself.

    Rows the road does not reach are interpolated from the rows it does, then
    held flat past the ends. Without that, the 200 m context pad and any row
    where the line steps outside the window would punch holes in the mask, and
    the override would be ragged along the domain edges for no physical reason.

    Returns (mask, n_rows_with_road). An all-False mask with 0 rows means this
    window has no road at all - the caller decides what that means.
    """
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


# =============================================================================
# THE DIAGNOSTIC THE WINDOW ORIGIN DEPENDS ON
# =============================================================================

def start_beach_median_m(arr):
    """
    Median over the alongshore profiles of the extractor's start_beach, in
    metres from the ocean edge of the window. Smaller = further seaward.

    Mirrors HAT_dune_topo_extractor.load_domain exactly:
        ocean-first (arr[:, ::-1], since OCEAN_LOC="right")
        z = raw - MHW ; z below WATER_CLAMP_M pinned to the sentinel
        first index where z > BEACH_START_THR_M
    then default_window()'s median over profiles. Nodata is treated as the water
    sentinel here, which is what the clamp does to it downstream.

    This is computed at 1 m. The extractor works at 10 m, so divide by 10 for
    Barrier3D cells - the audit carries both.
    """
    z = arr[:, ::-1] - MHW_ELEVATION
    z = np.where(np.isnan(z), WATER_CLAMP_M, z)
    z[z < WATER_CLAMP_M] = WATER_CLAMP_M
    above = z > BEACH_START_THR_M
    sb = np.where(above.any(axis=1), above.argmax(axis=1), -1)
    v = sb[sb >= 0]
    return float(np.median(v)) if v.size else float("nan")


# =============================================================================
# ONE FILL STAGE
# =============================================================================

def select(base_valid, cand, res_x):
    """Connectivity, shared by both stages. `cand` is already coverage- and
    floor-filtered; this only drops what the island cannot reach."""
    if not REQUIRE_ISLAND_CONNECTION:
        return cand, 0
    bridge_px = int(round(GAP_BRIDGE_M / res_x / 2.0))
    keep = gf.island_connected(base_valid, cand, bridge_px)
    return keep, int(cand.sum()) - int(keep.sum())


# =============================================================================
# MAIN
# =============================================================================

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

        # Built for the AUDIT ONLY - neither is applied. They answer "how far
        # landward of each alignment did 1996 write?", which is what a reader
        # needs in order to judge this product now that the survey's own swath
        # edge is the only landward limit.
        oceans = {yr: ocean_side_mask(g, base.shape, t_big)
                  for yr, g in road_geoms.items()}
        ocean, n_road_rows = oceans[1984]
        ocean04 = oceans[2004][0]
        has_road = n_road_rows > 0

        # --- STAGE 1: the 1996 override, ocean side only -------------------
        # "Gap" means NO OTHER SURVEY SAW THIS CELL - 2009 or 2014. Not just
        # 2009. The first build tested `~valid09` and the guard below caught
        # it: 16 domains still moved landward, D73 by 21 m and D90 by 16 m.
        # The 2009 holes are not empty in the product, 2014 fills them, so
        # `~valid09` let 1996 put a -2 m swash return over a +1 m 2014
        # measurement - precisely the replacement the MHW floor exists to
        # stop, just wearing a different survey's name. The stated principle
        # was always "one survey may replace another's measurement, but not
        # with a wet return"; ANOTHER means any, so the test is coverage by
        # any other survey.
        #
        # has14 is used raw, before 2014's own floor and connectivity. That is
        # deliberate and conservative in the right direction: where 2014 has
        # any value at all, 1996 must clear MHW to displace it.
        covered_by_other = valid09 | has14
        # NO ROAD BOUNDARY - the one substantive change of 2026-08-26. `ocean`
        # is built above and deliberately NOT applied here. The landward limit
        # is the ALACE swath edge, measured at 429-979 m from the ocean edge
        # against an island extent of 1274-1999 m: the survey stops far short
        # of the sound unaided, so the road was never what held it back. In
        # domain 77 the 1996 reach is 429 m while the road sits at 717 m - the
        # boundary was not even binding there. Every rule below is unchanged.
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

        # Reported every run whether or not anything is applied, so "we chose
        # not to correct" stays a checkable claim rather than a comment.
        # base - fill, the sign convention of gf.estimate_bias. NEGATIVE
        # means 1996 sits ABOVE 2009, which is what erosion predicts and
        # what the run reports (about -0.25 m through the developed reach).
        ov = valid09 & has96
        bias96 = (float(np.nanmedian(base[ov] - g96[ov]))
                  if int(ov.sum()) >= gf.MIN_RING_CELLS else 0.0)
        # The ocean-side-of-1984 figure the docstring's datum argument was
        # built on. Kept so that argument stays checkable on its ORIGINAL
        # footprint now that the written population is a fifth larger, rather
        # than being quietly restated against a different set of cells.
        ov84 = ov & ocean
        bias96_o84 = (float(np.nanmedian(base[ov84] - g96[ov84]))
                      if int(ov84.sum()) >= gf.MIN_RING_CELLS else 0.0)
        seam96 = gf.seam_check(merged, valid09 & ~cand96, cand96)

        # --- STAGE 2: the 2014 gap fill ------------------------------------
        # Run against the POST-1996 surface: 1996 cells are measured ground and
        # are legitimate connectivity anchors. Landward of the road this is
        # identical to HAT_dem_gap_fill.py, because nothing changed there.
        valid_now = np.isfinite(merged)
        cand14 = ~valid_now & has14 & (g14 >= gf.FILL_MIN_ELEV_NAVD)
        n14_cov = int((~valid_now & has14).sum())
        n14_floor = n14_cov - int(cand14.sum())
        cand14, n14_conn = select(valid_now, cand14, res_x)
        if cand14.any():
            merged[cand14] = g14[cand14]

        # --- the before/after the window origin depends on ------------------
        # "Before" is rebuilt here rather than read from the 2014 product, so
        # the two sides of the comparison come from one code path and a stale
        # product on disk cannot silently change the answer.
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

        # --- write ----------------------------------------------------------
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
            # DIAGNOSTIC ONLY from 2026-08-26 - neither alignment gates the
            # override any more. road_line is retained because readers of this
            # audit key off it, and because "does this domain have a road line
            # at all" remains a real property of the domain.
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
    # EVERY domain now receives 1996, so the shift statistics run over all of
    # them. Before 2026-08-26 they ran over `with_road` only, because domains
    # with no road line got no override at all and their shift was
    # structurally zero - averaging those in would have diluted the number.
    # That no longer applies. `with_road` survives to report the split.
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
    # The threshold that matters is one Barrier3D cell, not one metre. This
    # diagnostic runs at 1 m; the extractor works at 10 m, so a shift under
    # 10 m cannot move the window by more than a single cell and usually moves
    # it by none. Both bands are printed - the sub-cell one is residual, the
    # whole-cell one is the wet-edge failure the split floor exists to prevent
    # and means the rule needs re-checking before this product is used.
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
