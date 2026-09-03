"""
HAT_dem_duneline_coverage.py

Does the 1996 ALACE swath in the `2009-2014-1996` DEM actually reach the 1984
dune line?

WHY THIS EXISTS
---------------
`2009-2014-1996` is the 1984-start DEM: a 1996 beach and foredune grafted onto
a 2009 backdune, with NO road boundary - the landward limit is the ALACE
swath's own edge. ALACE surveyed "from the low water line to the landward base
of the sand dunes", so the graft seam lands at the dune toe.

That is fine as long as the 1984 dune is INSIDE the swath. Where the swath
stops seaward of the 1984 dune line, the model's t=0 dune is a 2009 dune
wearing a 1996 beach, and no amount of care in the pick pass can recover the
1984 crest from a surface that does not contain it.

This script measures that, per profile, at 1 m, over all 90 domains. It
CHANGES NOTHING. It writes no elevation, edits no product, and proposes no
correction - it reports where the question has a bad answer.

THE REFERENCE LINES
-------------------
    duneline_1984.geojson   495 vertices   the initial condition being tested
    duneline_1997.geojson   581 vertices   the CONTEMPORANEOUS control

1997 is one year after the ALACE flight, so it is where the 1996 DEM's OWN
dune sits. Reporting the reach against both separates two things that a single
number confounds:

    "1996 never flew that far landward"     <- fails against BOTH lines
    "the dune moved between the dates"      <- fails against 1984 only

The 1984-1997 separation is reported per domain so the second term is a
measured quantity rather than an assumption. Island-wide it is small - the
`raw-duneline-geojson/README.md` splits the naive 1984-vs-row-0 offset as
feature +16.2 m, date +0.8 m - but it is not small everywhere, and the
per-domain column is the point.

BOTH FILES ARE UTM 18N IN METRES. That README lists 1984 as EPSG:26918 and
1997 as EPSG:3725; those are the NAD83 and NAD83(NSRS2007) realizations of the
same projection and the transform between them at Hatteras is 0.000 m. They
are reprojected on load anyway, but no datum shift is being absorbed silently.

WHAT IS MEASURED, AND THE SIGN CONVENTION
-----------------------------------------
Every cross-shore distance in the outputs is METRES LANDWARD FROM THE OCEAN
EDGE of the 2000 m domain window. Larger = further landward. This is the same
origin `start_beach_median_m` in HAT_dem_1984_mosaic.py uses and the same one
HAT_dune_topo_extractor.py works in (OCEAN_LOC="right", so arrays are read
ocean-first).

Per profile (raster row; 500 per domain at 1 m):

    d1984_m, d1997_m      the dune lines' own positions
    reach_contig_m        walk landward from the FIRST 1996 cell and stop at
                          the first cell that is not 1996. The solid swath.
    reach_max_m           the landward-most 1996 cell anywhere in the profile,
                          holes ignored. The swath's true extent.
    gap1984_contig_m      d1984_m - reach_contig_m
    gap1984_max_m         d1984_m - reach_max_m    (and the same two for 1997)

    A POSITIVE GAP MEANS 1996 STOPS SEAWARD OF THE LINE - the coverage is
    MISSING there. Negative means 1996 reaches past it.

Both reach rules are reported because ALACE coverage in the band is patchy
(30-84%, median ~53%). Where the two diverge widely, what sits landward of the
contiguous swath is speckle rather than surface, and that divergence
(`reach_spread_m`) is itself a finding. Neither rule is applied to the other's
exclusion, and no hole-bridging tolerance is invented.

Profiles with no 1996 at all are given reach 0 - the swath reaches the ocean
edge and no further - rather than being dropped. `n_rows_no_1996` counts them
so a domain whose median rests on empty profiles is visible.

ABSENT IS NOT THE SAME AS REJECTED
----------------------------------
A cell can lack 1996 because ALACE never flew it, or because ALACE flew it and
this pipeline threw the return away. Those are different problems with
different fixes, and `clip_domain_*_survey.tif` cannot tell them apart - it
records only the winner.

So stage 1 of HAT_dem_1984_mosaic.py is RE-RUN here, in the padded window, with
its guards imported rather than restated, and every cell is classified:

    0  absent              ALACE has no data for this cell
    1  written             1996 won; this is what the DEM carries
    2  rej_ceiling         above 12.00 m NAVD88, the uncorrected-return tail
    3  rej_floor_gap       below -2.64 m NAVD88 where NO other survey saw it
    4  rej_floor_replace   below MHW where another survey did - a wet swash
                           return that would have displaced dry measured beach
    5  rej_connectivity    passed the floors, unreachable from the island

The recomputed `written` mask is checked cell for cell against the shipped
`clip_domain_<N>_survey.tif`. A nonzero mismatch means this diagnostic and the
product on disk have drifted apart; it is printed per domain and totalled at
the end. It should be zero.

THE BAND METRIC
---------------
Separately from the reach test, the composition of THE EXTRACTOR'S OWN WINDOW
is reported: 0-80 m landward of each profile's own beach start, the default
search band HAT_dune_topo_extractor.py picks dune crests in. That window is
anchored to beach start, not to the dune line, so it answers a different
question - "will the pick pass be picking in 1996 or in 2009?" - and is kept in
its own columns rather than blended with the reach numbers.

THE PER-DOMAIN FLAG
-------------------
    dune84_carried_by_1996 = median over profiles of gap1984_contig_m <= 0

Position alone, against the contiguous swath. No coverage threshold is
involved: a threshold would be a number picked out of the air, and the question
asked was where the fill stops relative to the 1984 line. The band fractions
sit beside the flag for a reader who wants to weigh it, and every input to it
is in the profile CSV.

INPUTS
    data/.../0-elevation/2009-2014-1996/1-gapfill-1m/clip_domain_*_survey.tif
    D:/Hatteras_GIS/.../2009_full.tif, 2014_full.tif, 1996_FallEC_J1441002/
    D:/Hatteras_GIS/domains.geojson
    data/.../1-barrier3d-domains/raw-duneline-geojson/duneline_1984.geojson
    data/.../1-barrier3d-domains/raw-duneline-geojson/duneline_1997.geojson

OUTPUTS (data/hatteras_init/0-elevation/2009-2014-1996-duneline/)
    duneline_coverage_domains.csv     90 rows, one per domain
    duneline_coverage_profiles.csv    45,000 rows, one per 1 m profile
    1-alace-class-10m/clip_domain_<N>_alaceclass.tif
                                      the six-code classification at 10 m,
                                      modal over each 10 x 10 block
    figures/                          drawn by
                                      3-figures/HAT_plot_duneline_coverage.py

Requires: rasterio, geopandas, numpy, scipy

    python HAT_dem_duneline_coverage.py
    python HAT_dem_duneline_coverage.py --domains 8,9,77   # a subset, to check

`--domains` is for checking the code path on a few windows. It writes the CSVs
for that subset only, so a subset run OVERWRITES the full ones - re-run without
it before reading anything.
"""

from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

import numpy as np
import geopandas as gpd
import rasterio
from rasterio.features import rasterize

sys.path.insert(0, str(Path(__file__).resolve().parent))
import HAT_dem_gap_fill as gf
import HAT_dem_1984_mosaic as m84

sys.path.insert(0, str(gf.PROJECT_ROOT / "scripts"))
from hat_elevation_products import ELEVATION_ROOT, product as _product  # noqa: E402


# =============================================================================
# CONFIG
# =============================================================================

SOURCE_TAG = "2009-2014-1996"

# A DIAGNOSTIC SIBLING, NOT A PRODUCT. It holds no elevation raster and forks
# nothing: the 178 MB of 1 m tifs stay where they are and are read in place, so
# this folder and the product it describes cannot drift apart on disk. It is
# deliberately NOT registered in hat_elevation_products.PRODUCTS - product()
# resolves things that have gapfill_1m and resampled_10m stages, and this has
# neither.
OUT_DIR = ELEVATION_ROOT / f"{SOURCE_TAG}-duneline"
CLASS_DIR = OUT_DIR / "1-alace-class-10m"
DOMAIN_CSV = "duneline_coverage_domains.csv"
PROFILE_CSV = "duneline_coverage_profiles.csv"

DUNE_DIR = (gf.PROJECT_ROOT / "data" / "hatteras_init" / "1-barrier3d-domains"
            / "raw-duneline-geojson")
DUNE_LINES = {1984: DUNE_DIR / "duneline_1984.geojson",
              1997: DUNE_DIR / "duneline_1997.geojson"}
TARGET_YEAR = 1984      # the line the flag is about
CONTROL_YEAR = 1997     # the contemporaneous control

# The extractor's default search band, 0-80 m landward of beach start. Changing
# this without changing HAT_dune_topo_extractor.py makes the band columns mean
# something the pick pass does not do.
BAND_START_M = 0.0
BAND_END_M = 80.0

# Classification codes. The order matters twice: it is the tie-break precedence
# for the 10 m downsample below, and it is the column order in the CSVs.
CLS_ABSENT = 0
CLS_WRITTEN = 1
CLS_CEILING = 2
CLS_FLOOR_GAP = 3
CLS_FLOOR_REPLACE = 4
CLS_CONNECTIVITY = 5
CLS_NAMES = {CLS_ABSENT: "absent", CLS_WRITTEN: "written",
             CLS_CEILING: "rej_ceiling", CLS_FLOOR_GAP: "rej_floor_gap",
             CLS_FLOOR_REPLACE: "rej_floor_replace",
             CLS_CONNECTIVITY: "rej_connectivity"}
N_CLS = len(CLS_NAMES)

# Most informative FIRST. A 10 x 10 block split evenly between "written" and a
# rejection is drawn as the rejection, because the figure exists to show where
# the fill fails and a tie that hid the failure would defeat it. Ties are rare
# and the count-based CSV is unaffected either way.
CLS_TIE_ORDER = [CLS_CONNECTIVITY, CLS_FLOOR_REPLACE, CLS_FLOOR_GAP,
                 CLS_CEILING, CLS_ABSENT, CLS_WRITTEN]

CLASS_NODATA = 255
GRID_SIZE_M = gf.GRID_SIZE_M      # 10 m, the Barrier3D cell


# =============================================================================
# THE REFERENCE LINES
# =============================================================================

def line_col_per_row(geom, shape, transform, min_rows=10):
    """
    Column index of a digitized line, one per raster row.

    Deliberately DIFFERENT from HAT_dem_1984_mosaic.ocean_side_mask, which takes
    the MAX column. That is right for a boundary - it yields the smallest ocean
    region and cannot place a road cell on the ocean side of itself - but this
    is not a boundary. It is a reference POSITION, and where the line runs
    diagonally through a row it occupies several columns with no reason to
    prefer either end. The MEAN is used, and the within-row span comes back
    alongside it so a reader can see how diagonal the line is where a profile's
    number looks odd.

    Rows the line does not reach are interpolated from the rows it does and held
    flat past the ends, exactly as ocean_side_mask does, so a line that steps
    briefly outside the padded window does not punch a hole in the series.

    Returns (col_per_row, span_per_row, n_rows_hit). n_rows_hit below min_rows
    means the line does not meaningfully cross this window and the column series
    comes back all-NaN - the caller decides what that means.
    """
    rl = rasterize([(geom, 1)], out_shape=shape, transform=transform,
                   fill=0, all_touched=True).astype(bool)

    col = np.full(shape[0], np.nan)
    span = np.full(shape[0], np.nan)
    for r in range(shape[0]):
        cs = np.where(rl[r])[0]
        if cs.size:
            col[r] = cs.mean()
            span[r] = cs.max() - cs.min()

    have = np.isfinite(col)
    n_have = int(have.sum())
    if n_have < min_rows:
        return np.full(shape[0], np.nan), span, n_have

    idx = np.arange(shape[0])
    col = np.interp(idx, idx[have], col[have])
    return col, span, n_have


# =============================================================================
# THE REACH RULES
# =============================================================================

def reach_indices(written_ocean_first):
    """
    Per profile, in ocean-first column indices:

        first    the first 1996 cell walking landward from the ocean edge
        contig   the last cell of the CONTIGUOUS run that starts there
        far      the landward-most 1996 cell anywhere in the profile

    A profile with no 1996 gets -1 in all three; the caller maps that to a reach
    of 0 m, meaning the swath reaches the ocean edge and no further. That is a
    measurement, not a fill - those profiles are counted separately as
    n_rows_no_1996 so a median resting on them is visible rather than implied.
    """
    w = written_ocean_first
    _, width = w.shape
    any96 = w.any(axis=1)

    first = np.where(any96, w.argmax(axis=1), -1)

    cols = np.arange(width)[None, :]
    gap = (~w) & (cols >= first[:, None])
    has_gap = gap.any(axis=1)
    first_false = np.where(has_gap, gap.argmax(axis=1), width)
    contig = np.where(any96, first_false - 1, -1)

    far = np.where(any96, width - 1 - w[:, ::-1].argmax(axis=1), -1)
    return first, contig, far


def start_beach_per_row(arr_crop):
    """
    HAT_dune_topo_extractor's start_beach, per profile rather than as a median.

    Mirrors HAT_dem_1984_mosaic.start_beach_median_m cell for cell - ocean-first,
    minus MHW, clamped at WATER_CLAMP_M, first index strictly above
    BEACH_START_THR_M - and differs from it only in not taking the median. The
    constants are imported from that module rather than restated, so the two
    cannot drift.

    Returns ocean-first indices, -1 where the profile never clears the threshold.
    """
    z = arr_crop[:, ::-1] - m84.MHW_ELEVATION
    z = np.where(np.isnan(z), m84.WATER_CLAMP_M, z)
    z[z < m84.WATER_CLAMP_M] = m84.WATER_CLAMP_M
    above = z > m84.BEACH_START_THR_M
    return np.where(above.any(axis=1), above.argmax(axis=1), -1)


# =============================================================================
# THE 10 m CLASSIFICATION
# =============================================================================

def downsample_class(cls_1m, block):
    """
    Modal class over each block x block cell, ties broken by CLS_TIE_ORDER.

    NOT the same rule as HAT_dem_resample_clip.downsample_survey, and the
    difference is deliberate. That function reports the provenance of the four
    cells bilinear actually reads, because its output has to describe the
    elevation value written beside it. Nothing here writes an elevation. This
    raster exists to be looked at, so it reports what MOST of the block is, and
    the tie-break only decides the rare even split.
    """
    h, w = cls_1m.shape
    b = cls_1m.reshape(h // block, block, w // block, block)
    counts = np.stack([(b == c).sum(axis=(1, 3)) for c in range(N_CLS)])
    # A weight strictly smaller than one cell, so it orders ties and nothing
    # else.
    bonus = np.zeros(N_CLS)
    for rank, c in enumerate(CLS_TIE_ORDER):
        bonus[c] = (len(CLS_TIE_ORDER) - rank) / (block * block * 10.0)
    return (counts + bonus[:, None, None]).argmax(axis=0).astype(np.uint8)


# =============================================================================
# MAIN
# =============================================================================

def main(only_domains=None):
    src_product = _product(SOURCE_TAG)
    shipped_1m = src_product.gapfill_1m

    checks = [("base DEM", gf.BASE_DEM_PATH),
              ("domain file", gf.DOMAIN_FILE),
              ("2014 fill", gf.FILL_DEM_PATH),
              ("1996 override", m84.OVERRIDE_DIR),
              (f"shipped {SOURCE_TAG}", shipped_1m)]
    checks += [(f"{y} dune line", q) for y, q in DUNE_LINES.items()]
    for label_, path in checks:
        if not path.exists():
            raise SystemExit(
                f"{label_} not found: {path}\n"
                f"  Fix the paths at the top of this script.")

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    CLASS_DIR.mkdir(parents=True, exist_ok=True)

    src = rasterio.open(gf.BASE_DEM_PATH)
    res_x, res_y = src.transform.a, -src.transform.e
    block = int(round(GRID_SIZE_M / res_x))
    pad_px = int(round(m84.CONTEXT_BUFFER_M / res_x))

    gdf = gpd.read_file(gf.DOMAIN_FILE)
    gdf = gf.resolve_crs(src, gdf)

    dem_crs = gf.CRS.from_wkt(src.crs.to_wkt())
    hz = dem_crs.sub_crs_list[0] if dem_crs.is_compound else dem_crs
    dune_geoms = {}
    for yr, q in DUNE_LINES.items():
        dl = gpd.read_file(q)
        n_v = sum(len(g.coords) for g in dl.geometry
                  if g is not None and g.geom_type == "LineString")
        print(f"{yr} dune line: {q.name}, {dl.crs} -> EPSG:{hz.to_epsg()}, "
              f"{n_v} vertices")
        dune_geoms[yr] = dl.to_crs(hz).union_all()

    fill14 = gf.FillSource(gf.FILL_DEM_PATH, src.crs)
    ov96 = m84.TiledSource(m84.OVERRIDE_DIR, m84.OVERRIDE_GLOB, src.crs)
    print(f"\n2009 DEM : {src.height} x {src.width}, {res_x} m")
    print(f"override : {m84.OVERRIDE_DIR.name}  EPSG:{ov96.epsg}  "
          f"{ov96.res} m")
    print(f"shipped  : {shipped_1m}")
    print(f"\nDistances are metres LANDWARD from the ocean edge of the 2000 m "
          f"window.\nA POSITIVE GAP means 1996 stops SEAWARD of the dune line "
          f"- coverage is MISSING.\n")

    domain_rows, profile_rows = [], []
    total_mismatch = 0

    for _, row in gdf.iterrows():
        dom = row[gf.DOMAIN_ID_FIELD]
        try:
            dom = int(dom)
        except (TypeError, ValueError):
            pass
        if only_domains is not None and dom not in only_domains:
            continue

        win, _adj = gf.snap_window(row.geometry.bounds, src.transform,
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

        # --- STAGE 1 OF THE MOSAIC, RE-RUN FOR ITS REJECTS ------------------
        # Every constant comes from m84. If that script's guards change, this
        # diagnostic changes with them, and the shipped-raster check below is
        # what proves the two are still in step.
        covered_by_other = valid09 | has14
        cand96_cov = has96.copy()
        if m84.APPLY_OVERRIDE_CEILING:
            too_high = cand96_cov & (g96 > m84.OVERRIDE_MAX_ELEV_NAVD)
            cand96_cov = cand96_cov & ~too_high
        else:
            too_high = np.zeros_like(has96)

        into_gap = (cand96_cov & ~covered_by_other
                    & (g96 >= m84.OVERRIDE_FLOOR_INTO_GAP))
        replace = (cand96_cov & covered_by_other
                   & (g96 >= m84.OVERRIDE_FLOOR_TO_REPLACE))
        floor_gap_rej = cand96_cov & ~covered_by_other & ~into_gap
        floor_rep_rej = cand96_cov & covered_by_other & ~replace

        passed = into_gap | replace
        cand96, _ = m84.select(valid09, passed, res_x)
        conn_rej = passed & ~cand96

        merged = base.copy()
        if cand96.any():
            merged[cand96] = g96[cand96]
        valid_now = np.isfinite(merged)
        cand14 = ~valid_now & has14 & (g14 >= gf.FILL_MIN_ELEV_NAVD)
        cand14, _ = m84.select(valid_now, cand14, res_x)
        if cand14.any():
            merged[cand14] = g14[cand14]

        r0 = c0 = pad_px
        crop = np.s_[r0:r0 + win.height, c0:c0 + win.width]

        cls = np.full(base.shape, CLS_ABSENT, np.uint8)
        cls[too_high] = CLS_CEILING
        cls[floor_gap_rej] = CLS_FLOOR_GAP
        cls[floor_rep_rej] = CLS_FLOOR_REPLACE
        cls[conn_rej] = CLS_CONNECTIVITY
        cls[cand96] = CLS_WRITTEN
        cls_c = cls[crop]

        # --- IS THIS STILL THE PRODUCT ON DISK? -----------------------------
        # The recomputed winners against the shipped provenance raster, cell for
        # cell. Nonzero means this diagnostic is describing a DEM that is not
        # the one in 1-gapfill-1m, and every number below it is suspect.
        sp = shipped_1m / f"clip_domain_{dom}_survey.tif"
        with rasterio.open(sp) as s:
            shipped = s.read(1)
        mismatch = int(((shipped == m84.SURVEY_1996)
                        != (cls_c == CLS_WRITTEN)).sum())
        total_mismatch += mismatch

        # --- REACH AND LINES, OCEAN-FIRST -----------------------------------
        width = cls_c.shape[1]
        w96 = (cls_c == CLS_WRITTEN)[:, ::-1]
        first_i, contig_i, far_i = reach_indices(w96)
        any96 = first_i >= 0
        reach_contig_m = np.where(any96, contig_i, 0) * res_x
        reach_max_m = np.where(any96, far_i, 0) * res_x

        sb_i = start_beach_per_row(merged[crop])
        sb_m = np.where(sb_i >= 0, sb_i, np.nan) * res_x

        d_m, span_m, n_hit = {}, {}, {}
        for yr, geom in dune_geoms.items():
            col, span, n = line_col_per_row(geom, base.shape, t_big)
            col_c = col[r0:r0 + win.height] - c0
            d_m[yr] = (width - 1 - col_c) * res_x     # ocean-first metres
            span_m[yr] = span[r0:r0 + win.height] * res_x
            n_hit[yr] = n

        # --- THE EXTRACTOR'S OWN WINDOW -------------------------------------
        cols = np.arange(width)[None, :]
        sb0 = np.where(sb_i >= 0, sb_i, 0)[:, None]
        band = ((cols >= sb0 + int(BAND_START_M / res_x))
                & (cols < sb0 + int(BAND_END_M / res_x))
                & (sb_i >= 0)[:, None])
        cls_of = cls_c[:, ::-1]
        band_n = band.sum(axis=1)
        band_frac = {}
        for c in range(N_CLS):
            hit = (band & (cls_of == c)).sum(axis=1)
            band_frac[c] = np.where(band_n > 0, hit / np.maximum(band_n, 1),
                                    np.nan)

        # --- PROFILE ROWS ---------------------------------------------------
        gaps = {}
        for yr in DUNE_LINES:
            gaps[(yr, "contig")] = d_m[yr] - reach_contig_m
            gaps[(yr, "max")] = d_m[yr] - reach_max_m
        sep = d_m[TARGET_YEAR] - d_m[CONTROL_YEAR]

        def _r(v, nd=2):
            return round(float(v), nd) if np.isfinite(v) else ""

        for r in range(cls_c.shape[0]):
            rec = {
                "domain": dom, "row": r,
                "start_beach_m": _r(sb_m[r]),
                "d1984_m": _r(d_m[1984][r]),
                "d1997_m": _r(d_m[1997][r]),
                "sep_1984_1997_m": _r(sep[r]),
                "line_span_1984_m": _r(span_m[1984][r]),
                "has_1996": bool(any96[r]),
                "reach_contig_m": _r(reach_contig_m[r]),
                "reach_max_m": _r(reach_max_m[r]),
                "reach_spread_m": _r(reach_max_m[r] - reach_contig_m[r]),
                "gap1984_contig_m": _r(gaps[(1984, "contig")][r]),
                "gap1984_max_m": _r(gaps[(1984, "max")][r]),
                "gap1997_contig_m": _r(gaps[(1997, "contig")][r]),
                "gap1997_max_m": _r(gaps[(1997, "max")][r]),
                "band_n": int(band_n[r]),
            }
            for c in range(N_CLS):
                rec[f"band_frac_{CLS_NAMES[c]}"] = _r(band_frac[c][r], 4)
            profile_rows.append(rec)

        # --- DOMAIN ROW -----------------------------------------------------
        def pct(a, p):
            a = np.asarray(a, float)
            a = a[np.isfinite(a)]
            return float(np.percentile(a, p)) if a.size else float("nan")

        med_gap84 = pct(gaps[(1984, "contig")], 50)
        carried = bool(med_gap84 <= 0)

        d = {
            "domain": dom,
            "n_rows": int(cls_c.shape[0]),
            "n_rows_no_1996": int((~any96).sum()),
            "n_rows_no_beach_start": int((sb_i < 0).sum()),
            "dune_rows_1984": n_hit[1984], "dune_rows_1997": n_hit[1997],
            "shipped_survey_mismatch": mismatch,

            "d1984_med_m": round(pct(d_m[1984], 50), 2),
            "d1997_med_m": round(pct(d_m[1997], 50), 2),
            "sep_1984_1997_med_m": round(pct(sep, 50), 2),
            "sep_1984_1997_p25_m": round(pct(sep, 25), 2),
            "sep_1984_1997_p75_m": round(pct(sep, 75), 2),

            "start_beach_med_m": round(pct(sb_m, 50), 2),
            "reach_contig_med_m": round(pct(reach_contig_m, 50), 2),
            "reach_max_med_m": round(pct(reach_max_m, 50), 2),
            "reach_spread_med_m": round(pct(reach_max_m - reach_contig_m, 50),
                                        2),

            "gap1984_contig_med_m": round(med_gap84, 2),
            "gap1984_contig_p25_m": round(pct(gaps[(1984, "contig")], 25), 2),
            "gap1984_contig_p75_m": round(pct(gaps[(1984, "contig")], 75), 2),
            "gap1984_max_med_m": round(pct(gaps[(1984, "max")], 50), 2),
            "gap1997_contig_med_m": round(pct(gaps[(1997, "contig")], 50), 2),
            "gap1997_max_med_m": round(pct(gaps[(1997, "max")], 50), 2),

            "gap1984_contig_med_cells": round(med_gap84 / GRID_SIZE_M, 2),
            "frac_rows_gap1984_positive":
                round(float(np.mean(gaps[(1984, "contig")] > 0)), 4),
            "dune84_carried_by_1996": carried,
        }
        for c in range(N_CLS):
            d[f"n_{CLS_NAMES[c]}"] = int((cls_c == c).sum())
            d[f"band_frac_{CLS_NAMES[c]}"] = round(
                float(np.nanmean(band_frac[c])), 4)
        domain_rows.append(d)

        # --- THE 10 m CLASSIFICATION ----------------------------------------
        t_out = src.window_transform(win)
        cls10 = downsample_class(cls_c, block)
        gf.write_raster(
            cls10,
            rasterio.Affine(GRID_SIZE_M, 0, t_out.c, 0, -GRID_SIZE_M, t_out.f),
            src.crs, CLASS_DIR / f"clip_domain_{dom}_alaceclass.tif",
            "uint8", CLASS_NODATA)

        flag = "1996" if carried else "2009  <-- 1984 dune NOT in the swath"
        warn = f"  MISMATCH {mismatch}" if mismatch else ""
        print(f"  domain {dom:>3}: d1984 {d['d1984_med_m']:>7.1f} m   "
              f"reach {d['reach_contig_med_m']:>7.1f} m   "
              f"gap {med_gap84:>+8.1f} m ({med_gap84 / GRID_SIZE_M:>+5.1f} "
              f"cells)   band 1996 {d['band_frac_written'] * 100:>5.1f}%   "
              f"{flag}{warn}")

    ov96.close()
    fill14.close()
    src.close()

    dp = OUT_DIR / DOMAIN_CSV
    with open(dp, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(domain_rows[0].keys()))
        w.writeheader()
        w.writerows(domain_rows)
    pp = OUT_DIR / PROFILE_CSV
    with open(pp, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(profile_rows[0].keys()))
        w.writeheader()
        w.writerows(profile_rows)

    # --- SUMMARY ------------------------------------------------------------
    n = len(domain_rows)
    carried = [r for r in domain_rows if r["dune84_carried_by_1996"]]
    missed = [r for r in domain_rows if not r["dune84_carried_by_1996"]]
    g = np.array([r["gap1984_contig_med_m"] for r in domain_rows])
    gm = np.array([r["gap1984_max_med_m"] for r in domain_rows])
    sep = np.array([r["sep_1984_1997_med_m"] for r in domain_rows])

    print(f"\n{'=' * 78}\n{n} domains")
    print(f"  1984 dune inside the 1996 swath : {len(carried)}")
    print(f"  1984 dune NOT inside it         : {len(missed)}")
    if missed:
        print(f"      {[r['domain'] for r in missed]}")
    print(f"\n  gap to the 1984 line, contiguous reach, m "
          f"(+ve = 1996 stops seaward):")
    print(f"      min {g.min():+.1f}   p25 {np.percentile(g, 25):+.1f}   "
          f"median {np.median(g):+.1f}   p75 {np.percentile(g, 75):+.1f}   "
          f"max {g.max():+.1f}")
    print(f"      in 10 m cells: median {np.median(g) / GRID_SIZE_M:+.1f}, "
          f"worst {g.max() / GRID_SIZE_M:+.1f}")
    print(f"  same, landward-most reach: median {np.median(gm):+.1f} m   "
          f"(speckle beyond the solid swath buys "
          f"{np.median(g) - np.median(gm):.1f} m)")
    print(f"\n  1984 - 1997 line separation, m (+ve = 1984 line landward of "
          f"1997): median {np.median(sep):+.1f}, "
          f"p25 {np.percentile(sep, 25):+.1f}, "
          f"p75 {np.percentile(sep, 75):+.1f}")

    tot = {c: sum(r[f"n_{CLS_NAMES[c]}"] for r in domain_rows)
           for c in range(N_CLS)}
    allc = max(sum(tot.values()), 1)
    print(f"\n  why 1996 is or is not there, island-wide, all cells:")
    for c in range(N_CLS):
        print(f"      {CLS_NAMES[c]:<18} {tot[c]:>12,}  "
              f"{tot[c] / allc * 100:>5.2f}%")
    print(f"  same, inside the extractor's {BAND_START_M:.0f}-{BAND_END_M:.0f} "
          f"m window, mean of per-domain fractions:")
    for c in range(N_CLS):
        v = np.nanmean([r[f"band_frac_{CLS_NAMES[c]}"] for r in domain_rows])
        print(f"      {CLS_NAMES[c]:<18} {v * 100:>7.2f}%")

    if total_mismatch:
        print(f"\n  *** {total_mismatch:,} cells disagree with the shipped "
              f"survey raster. This diagnostic is NOT describing the product "
              f"in {shipped_1m.name}. Do not use these numbers. ***")
    else:
        print(f"\n  shipped-raster check: 0 mismatched cells - the recomputed "
              f"1996 winners are exactly the product on disk.")

    print(f"\n  domains : {dp}")
    print(f"  profiles: {pp}")
    print(f"  class   : {CLASS_DIR}")


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[2])
    ap.add_argument("--domains", default=None,
                    help="comma-separated subset, for checking the code path. "
                         "Overwrites the full CSVs - re-run without it.")
    a = ap.parse_args()
    subset = ({int(x) for x in a.domains.split(",")} if a.domains else None)
    main(subset)
