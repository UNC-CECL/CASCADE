# ==============================================================================
# HAT_road_offset_from_dune_start.py
#
# NC-12 road setback and road elevation per Barrier3D domain, measured from the
# SAME reference CASCADE indexes against: interior row 0 of the extracted
# topography, which is one cell landward of the picked dune crest.
#
# WHY THIS EXISTS
#   `roadway_manager.bulldoze` places the road with
#
#       road_start = int(road_setback / dy)                 # roadway_manager.py:99
#       old_road_domain = xyz_interior_grid[road_start:road_end, :]
#
#   so `road_setback` is metres landward of interior row 0, and
#   `cascade_pipeline/roadway.py:147` already documents it as "metres landward of
#   the dune line". This script is the first measurement that actually honours
#   that convention: it re-derives the dune line with the same picked windows and
#   the same straightening the topography was built with, then measures the road
#   against it, per alongshore profile.
#
#   It also removes a large error source. The clip boxes are north-up while the
#   island trends NNE, so NC-12 crosses each 500 m domain diagonally: on GIS 11
#   the road spans raw cross-shore cells 151-161, i.e. 110 m of cross-shore
#   extent for a 20 m road. Any method that reduces that to a raw cross-shore
#   median inherits the smear. Shearing the road mask with the SAME per-profile
#   shear as the topography (`shear_like`) collapses it.
#
# TWO REFERENCE FRAMES, BOTH REPORTED
#   The topography is 2009; the road_offset are 1984 and 2004. Those disagree, and the
#   disagreement is the subject of old_method_offset/RoadSetback_audit.md. This
#   script does not pick a winner:
#     * setback_2009_m   measured directly here, road-vs-2009-dune. Internally
#                        consistent with the grid CASCADE actually runs.
#     * setback_sameyear_m  read from the EXISTING RoadSetback_<year>.csv, which
#                        is a road-vs-same-year-dune measurement. No
#                        extrapolation is performed to produce it.
#     * delta_m          the difference, which should track the dune-line
#                        retreat between <year> and 2009. The 2-brie-offset
#                        files are used ONLY to validate that, never as an input.
#
# OUTPUTS (nothing existing is overwritten; new tree under dunestart_offset\)
#   dunestart_offset\<year>\RoadSetback_<year>_dunestart.csv    2-row, model-facing
#   dunestart_offset\<year>\RoadElevation_<year>_dunestart.csv  2-row, m MHW
#   dunestart_offset\<year>\RoadOffset_<year>_domains.csv       per-domain detail
#   dunestart_offset\<year>\RoadOffset_<year>_profiles.csv      per-profile detail
#   dunestart_offset\RoadOffset_dunestart_audit.md              the write-up
#
# WHERE THIS DEPARTS FROM "MEASURE, DON'T CORRECT" -- TWO PLACES, BOTH FLAGGED
#   Both act ONLY on the model-facing CSV. `setback_2009_m` in the _domains.csv
#   is always the measurement, and every adjusted domain carries a flag naming
#   what was done to it. (1) is the negative floor, below. (2) is the seaward
#   relocation of roadways that drown at initialisation, documented at
#   RELOCATE_DROWNING.
#
# (1) NEGATIVE SETBACKS ARE FLOORED
#   A 1984 road measured against the 2009 dune line can land SEAWARD of interior
#   row 0, giving a negative setback. A negative cannot be written to the
#   model-facing file, because `int(-50/10) = -5` and
#   `xyz_interior_grid[-5:-3, :]` is valid Python that indexes from the LANDWARD
#   end -- the road would be bulldozed into the bay, silently. Nor is 0.0 a
#   "no road" sentinel: `build_roadway_management_on` decides which domains are
#   managed from the road SPAN, not the value, so a 0.0 inside the span still
#   means "managed, road at row 0".
#
#   So the model-facing CSV carries max(setback, 0.0) and the true signed value
#   is preserved in the _domains.csv and the audit, flagged NEGATIVE. That is a
#   floor, and it is called one everywhere it appears.
#
# USAGE
#     python HAT_road_offset_from_dune_start.py
# ==============================================================================

from __future__ import annotations

import importlib.util
import json
from datetime import datetime
from pathlib import Path

import numpy as np

import matplotlib

# Only the array helpers are used, never the picker.
matplotlib.use("Agg")

# ==============================================================================
# CONFIG
# ==============================================================================

PROJECT_ROOT = Path(r"C:\Users\hanna\PycharmProjects\CASCADE")
INIT_ROOT = PROJECT_ROOT / "data" / "hatteras_init"

# IMPORTANT: the topography_dunes copy, NOT the sibling copy in this folder.
# There are four copies of HAT_dune_topo_extractor.py in the repo and only this
# one has ALONGSHORE_FLIP = True. Importing the local copy would measure the road
# in the mirrored alongshore frame -- the exact mismatch this whole exercise is
# meant to eliminate.
EXTRACTOR = (PROJECT_ROOT / "scripts" / "input_prep" / "1-barrier3d-domains"
             / "topography_dunes" / "HAT_dune_topo_extractor.py")

ROADS_ROOT = INIT_ROOT / "4-mgmt-forcing" / "road_offset"
MASK_DIR_FMT = ROADS_ROOT / "raster" / "{year}" / "masks"
MASK_NAME_FMT = "domain_{domain}_road_{year}.npy"

# The existing same-year measurements, used as the second reference frame.
EXISTING_SETBACK_FMT = ROADS_ROOT / "old_method_offset" / "{year}" / "RoadSetback_{year}.csv"

# Offset files, used ONLY to validate delta_m against measured retreat.
OFFSET_FMT = (INIT_ROOT / "2-brie-offset" / "hindcast_{year}"
              / "Island_Dune_Offsets_{year}_CASCADE_Input.csv")

OUT_ROOT = ROADS_ROOT / "dunestart_offset"

YEARS = [1984, 2004]
DOMAINS = list(range(1, 91))        # D1 = Cape Point (south) -> D90 = Rodanthe

# --- ROAD SPAN ----------------------------------------------------------
# Domains written to the model-facing CSVs. Matches the legacy files and the
# rasterizer's own FIRST_ROAD_DOMAIN/LAST_ROAD_DOMAIN = 9, 90.
#
# D8 is measured and reported but EXCLUDED, on evidence rather than convention:
# at the Buxton bend NC-12 turns to run nearly east-west, PARALLEL to the raster
# rows, so it crosses D8's corner rather than crossing the domain shore-normal.
# HAT_check_geojson_vs_mask.py measures the consequence -- the road spans a
# median of 75 cells per row there (102 max) against 2 cells in a normal domain,
# and it touches only 25 of 50 profiles. A cross-shore "distance landward of the
# dune" is not a meaningful quantity for a road running alongshore, so a scalar
# setback for D8 would be a number without a physical reading.
#
# Domains outside this span are still measured, still land in the _domains.csv,
# and are flagged EXCLUDED_FROM_SPAN so the exclusion is visible rather than
# silent.
ROAD_SPAN = (9, 90)

# --- UNSTRAIGHTENED CONTROL PASS ----------------------------------------
# A second measurement with STRAIGHTEN = False, so the old-vs-new difference can
# be attributed instead of just reported. Holding the dune feature, the road
# mask, the aggregation and this code fixed and changing ONLY the frame splits
# the total change exactly:
#
#     new_straightened - legacy  =  (new_straightened - new_raw)   <- frame
#                                 + (new_raw          - legacy)    <- dune feature
#
# The control needs the pre-straightening pick set, because a window picked in
# one frame is a valid index range in the other that points at different cells --
# which is what the extractor's `straightened` guard exists to catch. That file
# predates the flag (straighten: null), so the guard reads it as False and lets
# the control through, which is correct.
#
# Residual caveat, stated rather than hidden: the two passes necessarily use two
# different pick files, so a small part of the "frame" component is really
# window-choice. The windows were drawn to bracket the same dune in each frame,
# so it is second-order, but it is not zero.
#
# ALONGSHORE_FLIP is deliberately NOT varied. The setback is one scalar per
# domain, collapsed by median over profiles; reversing the profile order leaves
# that multiset unchanged, so the flip cannot move the setback at all. It decides
# where the road sits WITHIN a domain, which is a different question.
CONTROL_UNSTRAIGHTENED = True
CONTROL_WINDOW_JSON = (INIT_ROOT / "1-barrier3d-domains" / "2009-dune-topo"
                       / "picks" / "HAT_dune_search_windows_2009_pea_hatteras.json")
CONTROL_SUFFIX = "_rawframe"

# Assumed roadway width. RoadwayConfig.road_width_m is a single global 20.0
# (cascade_pipeline/roadway.py:45), so a per-domain width is not consumable
# today -- but cascade.py:43 does accept an array, so the MEASURED width is
# reported as a diagnostic in case that changes.
ASSUMED_ROAD_WIDTH_M = 20.0

# Flag thresholds. None of these alter a value; they only label it.
SCATTER_SETBACK_SPREAD_M = 40.0     # p90-p10 of the per-profile setback
SCATTER_ELEV_STD_M = 0.50           # std of road-cell elevation
MIN_ROAD_PROFILES = 25              # of 50; below this the domain is PARTIAL
WATER_FRAC_FLAG = 0.20              # fraction of road cells at sentinel water

# --- SEAWARD RELOCATION OF DROWNING ROADWAYS ----------------------------
# The second departure from "measure, don't correct". Decided 2026-08-18.
#
# WHAT IT DOES
#   A roadway whose flanking rows are >20% water at t=0 is width-drowned by
#   `roadway_manager.bulldoze` on the first call: RoadwayManager sets
#   _drown_break, cascade.py never calls update() again, and that domain spends
#   the entire hindcast as an UNMANAGED barrier wearing a road label -- no
#   overwash removal, no dune rebuilding, no relocation. Rather than lose the
#   domain, the road is moved to the nearest row SEAWARD that survives the test.
#
# WHAT "NEAREST VIABLE" MEANS
#   The largest road_start < current for which bulldoze's own flanking test
#   passes -- `drown_test` below is that test, not an approximation of it. The
#   road's own band is NOT required to be dry, because bulldoze never looks at
#   it; requiring it would move roads further seaward than the model needs.
#   There is NO cap on the distance: the nearest viable row is taken however far
#   it is. See the caveat below for what that costs on GIS 79.
#
# THE ASSUMPTION THIS RESTS ON, STATED PLAINLY
#   On the 2009_v3 topography the domains this fires on do NOT drown on measured
#   water. They drown on LiDAR coverage gaps. Across the six flanking rows that
#   fail at GIS 78/79/80 there are 106 wet cells, of which 105 were NEVER
#   SURVEYED and 1 is genuinely measured wet. The extractor writes no-data back
#   as SENTINEL_WATER_M because Barrier3D has no representation for "unknown",
#   and `wet_fraction` counts every cell (see the audit script -- that is
#   deliberate, so the figure reports what the run does).
#
#   So this relocation is NOT "the road was in water, we moved it out". It is
#   "the 2009 DEM has no data there, CASCADE reads no-data as water, so the road
#   is moved onto surveyed ground to keep the domain managed". Anyone reading a
#   managed-vs-unmanaged result at GIS 78-80 needs that sentence.
#
# WHAT IT COSTS, PER DOMAIN
#   GIS 78  490 -> 470 m ( 20 m)  lands on the measured per-profile MINIMUM
#   GIS 80  490 -> 450 m ( 40 m)  inside the measured spread, ~p10 (already
#                                 flagged SCATTER_SETBACK(71m))
#   GIS 79  510 -> 400 m (110 m)  90 m SEAWARD OF THE SEAWARD-MOST PROFILE ever
#                                 measured in that domain. This one is not a
#                                 re-pick of the alongshore statistic; it is a
#                                 position no profile showed. Accepted knowingly
#                                 (no cap), and flagged BEYOND_MEASURED so it can
#                                 never be mistaken for a measurement.
#
# Set RELOCATE_DROWNING = False to write the measured setbacks unchanged and let
# the domains drown; everything is still measured and flagged either way.
RELOCATE_DROWNING = True

# bulldoze's own constants. Not tuning knobs -- roadway_manager.py:50-51,125-165;
# RoadwayManager re-asserts 0.2 at :531; drown_threshold = 0 at :766.
DROWN_THRESHOLD_M = 0.0
DROWN_PCT = 0.2
ROAD_WIDTH_CELLS = 2                # int(road_width 20 m / dx 10 m)


# ==============================================================================
# EXTRACTOR IMPORT
# ==============================================================================

def load_extractor():
    """Import the corrected extractor and assert it is the flipped one."""
    spec = importlib.util.spec_from_file_location("hat_extractor", EXTRACTOR)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)

    if not module.ALONGSHORE_FLIP:
        raise SystemExit(
            f"\n{EXTRACTOR} has ALONGSHORE_FLIP = False.\n"
            "The road offsets must be measured in the same alongshore frame as\n"
            "the topography. Either point EXTRACTOR at the corrected copy or\n"
            "set ALONGSHORE_FLIP = True there.\n"
        )
    if not module.STRAIGHTEN:
        print("[warn] extractor has STRAIGHTEN = False; the road mask will not "
              "be sheared, so obliquity stays in the measurement.")
    return module


# ==============================================================================
# FRAME ALIGNMENT
#
# align_mask_to_topography() and interior_row0_line() USED TO LIVE HERE. They
# were private copies that re-derived the extractor's frame from outside it --
# two definitions of one chain, kept in step by hand. They now live in
# HAT_dune_topo_extractor.py next to shear_like(), which is the only place that
# knows the frame, and this script calls them through `ext`:
#
#     ext.align_mask_to_topography(raw_mask, dom)     (no `ext` first arg)
#     ext.interior_row0_line(prof_arr, dune_loc)
#
# The extractor needed them anyway, to draw NC-12 on the picker and to report
# road distances from interior row 0 in its settings sheet. Moving rather than
# duplicating is what keeps this script's setback and the extractor's diagnostic
# columns the same measurement.
#
# TWO LATENT BUGS WERE FIXED IN THE MOVE, both no-ops on the current settings:
# interior_row0_line's all-water test now matches remove_water_rows (`<=`, so a
# leading row of pure no-data trims like the real thing instead of being kept),
# and lead_trim is only applied when TRIM_INTERIOR_ROWS is True. lead_trim is 0
# on all 90 domains either way, which is why the outputs did not move -- verified
# by diffing the CSVs before and after the refactor.
# ==============================================================================


# ==============================================================================
# GEOLOCATION  (diagnostics only -- never feeds a setback)
#
# Put the measured per-profile cells back on the map, so RoadOffset_*_profiles.csv
# can be loaded into GIS and eyeballed against the road instead of being read as
# bare cell indices. Borrowed in spirit from
# roya_files/road_setback_from_road_json_to_interior_start_REVISED.py, which works
# natively in map coordinates -- but NOT its method: that script reconstructs the
# interior point by walking (c0 + interior_start + 0.5) * cell metres from the
# domain polygon edge along a PCA axis, a formula with no shear term, which would
# put back the obliquity this pipeline exists to remove. Here the affine comes
# straight from the same resampled_domain_*.tif the mask was burned onto, and the
# index chain is inverted exactly.
# ==============================================================================

TIF_FMT = (INIT_ROOT / "1-barrier3d-domains" / "2009-raw"
           / "2009-domain-clipresample" / "domain_{domain}"
           / "resampled_domain_{domain}.tif")

WRITE_PROFILE_COORDS = True
_geo_warned: set[str] = set()


def domain_georeference(ext, domain: int):
    """(transform, crs, n_rows, n_cols) for one domain, or None if unavailable.

    Returns None rather than raising: the coordinates are a QC convenience, and a
    missing tif or a missing rasterio must never stop a setback being measured.
    """
    if not WRITE_PROFILE_COORDS:
        return None
    try:
        import rasterio
    except ImportError:
        if "rasterio" not in _geo_warned:
            _geo_warned.add("rasterio")
            print("  [coords] rasterio not installed; profile x/y left blank")
        return None

    # Only OCEAN_LOC = "right" is invertible with the simple chain below. The
    # "top"/"bottom" branches apply np.rot90, which swaps the axes -- inverting
    # that needs its own case, and getting it silently wrong would put points in
    # the wrong place on a map, which is worse than leaving them out.
    if ext.OCEAN_LOC != "right":
        if "ocean_loc" not in _geo_warned:
            _geo_warned.add("ocean_loc")
            print(f"  [coords] OCEAN_LOC={ext.OCEAN_LOC!r} is not invertible here; "
                  f"profile x/y left blank")
        return None

    path = Path(str(TIF_FMT).format(domain=domain))
    if not path.is_file():
        if "tif" not in _geo_warned:
            _geo_warned.add("tif")
            print(f"  [coords] no resampled tif under {path.parent.parent}; "
                  f"profile x/y left blank")
        return None
    with rasterio.open(path) as src:
        return src.transform, src.crs, src.shape[0], src.shape[1]


def cell_to_map(geo, ext, dom: dict, profile: int, cross_cell: int):
    """Aligned-frame (profile, cross-shore cell) -> (x, y) in the tif's CRS.

    Inverts load_profiles' chain, in reverse order:

        aligned cell k  ->  ocean-first column   k + c0 + shear[profile]
        ocean-first     ->  oriented column      (n_cols - 1) - j
        ALONGSHORE_FLIP ->  original row         (n_rows - 1) - profile
        affine          ->  x, y at the CELL CENTRE (+0.5, +0.5)

    The +0.5 is the same convention HAT_check_geojson_vs_mask.py documents: cell k
    spans [k, k+1), so its centre is k + 0.5. Dropping it shifts every point half
    a cell (5 m) northwest.
    """
    if geo is None or cross_cell < 0:
        return "", ""
    transform, _crs, n_rows, n_cols = geo

    j = int(cross_cell) + int(dom["c0"]) + int(dom["shear"][profile])
    col = (n_cols - 1) - j
    row = (n_rows - 1) - profile if ext.ALONGSHORE_FLIP else profile
    if not (0 <= col < n_cols and 0 <= row < n_rows):
        return "", ""

    x, y = transform * (col + 0.5, row + 0.5)
    return round(float(x), 2), round(float(y), 2)


# ==============================================================================
# PER-DOMAIN MEASUREMENT
# ==============================================================================

def measure_domain(ext, domain: int, year: int, windows: dict) -> tuple[dict, list[dict]]:
    """Measure one domain's road setback and elevation against its dune start."""
    stem = f"domain_{domain}"
    dem_path = ext.LOAD_PATH / f"{stem}.npy"
    mask_path = Path(str(MASK_DIR_FMT).format(year=year)) / \
        MASK_NAME_FMT.format(domain=domain, year=year)

    record = {
        "year": year, "domain": domain,
        "section": ext.section_for(stem),
        "setback_2009_m": np.nan,
        "setback_2009_floored_m": 0.0,
        # What the model-facing CSV actually carries, and why it differs.
        # Defaults live here rather than being added later so every record --
        # in-span, excluded, no-road, and the control pass -- has the same
        # fieldnames for write_csv's DictWriter.
        "setback_model_m": np.nan,
        "drowns_at_init": "",
        "relocated_seaward_m": 0.0,
        "n_road_profiles": 0,
        "setback_p10_m": np.nan, "setback_p90_m": np.nan,
        "setback_min_m": np.nan, "setback_max_m": np.nan,
        "setback_center_m": np.nan,
        "measured_road_width_m": np.nan,
        "road_elev_mhw_median": np.nan,
        "road_elev_mhw_mean": np.nan,
        "road_elev_std_m": np.nan,
        "road_elev_navd_median": np.nan,
        "n_road_cells": 0,
        "water_frac": np.nan,
        "obliquity_deg": np.nan,
        "shear_max_cells": 0,
        "lead_trim_rows": 0,
        "window": "",
        "flags": "",
    }
    flags: list[str] = []

    if not dem_path.is_file():
        record["flags"] = "NO_DEM"
        return record, []
    if not mask_path.is_file():
        record["flags"] = "NO_MASK"
        return record, []

    dom = ext.load_profiles(dem_path)
    prof_arr = dom["z"]
    record["obliquity_deg"] = dom["obliquity_deg"]
    record["shear_max_cells"] = int(np.max(dom["shear"]))

    # Same window the topography was extracted with, and the same frame guard the
    # extractor's run pass applies -- a window picked unstraightened is a valid
    # index range that points at different cells.
    w = windows.get(stem)
    if w is None:
        i0, i1 = ext.default_window(prof_arr, dom["start_beach"])
        record["window"] = f"default[{i0},{i1}]"
        flags.append("DEFAULT_WINDOW")
    else:
        if bool(w.get("straightened", False)) != bool(ext.STRAIGHTEN):
            record["flags"] = "WINDOW_FRAME_MISMATCH"
            return record, []
        i0, i1 = int(w["i0"]), int(w["i1"])
        record["window"] = f"[{i0},{i1}]"

    dune_elev, dune_loc = ext.find_dunes(prof_arr, dom["start_beach"], i0, i1)
    if not np.any(dune_loc >= 0):
        record["flags"] = "NO_DUNE"
        return record, []

    row0, lead_trim = ext.interior_row0_line(prof_arr, dune_loc)
    record["lead_trim_rows"] = lead_trim

    mask = ext.align_mask_to_topography(np.load(mask_path), dom)
    if not mask.any():
        record["flags"] = "NO_ROAD"
        return record, []

    n_along = min(ext.ALONG_COLS, prof_arr.shape[0])
    cell = ext.CELL_SIZE_M
    geo = domain_georeference(ext, domain)   # None -> x/y columns stay blank

    profiles: list[dict] = []
    setback_cells: list[float] = []
    center_cells: list[float] = []
    width_cells: list[float] = []
    elev_mhw: list[float] = []
    n_water = 0
    n_cells = 0

    for i in range(n_along):
        road_cells = np.flatnonzero(mask[i])
        if road_cells.size == 0 or row0[i] < 0:
            continue

        # road_setback locates the SEAWARD edge of the road block, because
        # bulldoze indexes [road_start : road_start + road_width]. Ocean-first
        # cross-shore means the seaward edge is the minimum index.
        seaward = int(road_cells.min())
        landward = int(road_cells.max())
        center = float(road_cells.mean())

        sb_cells = seaward - int(row0[i])
        setback_cells.append(sb_cells)
        center_cells.append(center - int(row0[i]))
        width_cells.append(landward - seaward + 1)

        profile_elev = prof_arr[i, road_cells]
        wet = ~(profile_elev > ext.SENTINEL_WATER_M + 1e-9)
        n_water += int(wet.sum())
        n_cells += int(road_cells.size)
        elev_mhw.extend(profile_elev[~wet].tolist())

        # Map coordinates of the two cells the setback is measured between, so
        # the pair can be plotted in GIS against the road. Diagnostics only --
        # setback_m above is computed from the indices, never from these.
        interior_x, interior_y = cell_to_map(geo, ext, dom, i, int(row0[i]))
        road_x, road_y = cell_to_map(geo, ext, dom, i, seaward)

        profiles.append({
            "year": year, "domain": domain, "profile": i,
            "dune_crest_cell": int(dune_loc[i]),
            "interior_row0_cell": int(row0[i]),
            "road_seaward_cell": seaward,
            "road_landward_cell": landward,
            "road_center_cell": round(center, 2),
            "setback_m": round(sb_cells * cell, 1),
            "road_width_m": round((landward - seaward + 1) * cell, 1),
            "n_road_cells": int(road_cells.size),
            "n_wet_road_cells": int(wet.sum()),
            "interior_x": interior_x, "interior_y": interior_y,
            "road_x": road_x, "road_y": road_y,
        })

    if not setback_cells:
        record["flags"] = "NO_ROAD"
        return record, []

    sb = np.asarray(setback_cells, dtype=float) * cell
    record["n_road_profiles"] = len(sb)
    record["setback_2009_m"] = round(float(np.median(sb)), 1)
    record["setback_2009_floored_m"] = round(max(float(np.median(sb)), 0.0), 1)
    record["setback_p10_m"] = round(float(np.percentile(sb, 10)), 1)
    record["setback_p90_m"] = round(float(np.percentile(sb, 90)), 1)
    record["setback_min_m"] = round(float(sb.min()), 1)
    record["setback_max_m"] = round(float(sb.max()), 1)
    record["setback_center_m"] = round(
        float(np.median(np.asarray(center_cells) * cell)), 1)
    record["measured_road_width_m"] = round(
        float(np.median(np.asarray(width_cells) * cell)), 1)

    record["n_road_cells"] = n_cells
    record["water_frac"] = round(n_water / n_cells, 3) if n_cells else np.nan
    if elev_mhw:
        arr = np.asarray(elev_mhw, dtype=float)
        record["road_elev_mhw_median"] = round(float(np.median(arr)), 3)
        record["road_elev_mhw_mean"] = round(float(arr.mean()), 3)
        record["road_elev_std_m"] = round(float(arr.std()), 3)
        record["road_elev_navd_median"] = round(
            float(np.median(arr)) + ext.MHW_M, 3)
    else:
        flags.append("ALL_ROAD_CELLS_WET")

    # --- flags: labels only, no value is altered ---
    if record["setback_2009_m"] < 0:
        flags.append(f"NEGATIVE({record['setback_2009_m']:.0f}->floored 0)")
    if record["n_road_profiles"] < MIN_ROAD_PROFILES:
        flags.append(f"PARTIAL({record['n_road_profiles']}/{n_along})")
    spread = record["setback_p90_m"] - record["setback_p10_m"]
    if spread > SCATTER_SETBACK_SPREAD_M:
        flags.append(f"SCATTER_SETBACK({spread:.0f}m)")
    if (np.isfinite(record["road_elev_std_m"])
            and record["road_elev_std_m"] > SCATTER_ELEV_STD_M):
        flags.append(f"SCATTER_ELEV({record['road_elev_std_m']:.2f})")
    if np.isfinite(record["water_frac"]) and record["water_frac"] > WATER_FRAC_FLAG:
        flags.append(f"IN_WATER({record['water_frac']:.0%})")
    if abs(record["measured_road_width_m"] - ASSUMED_ROAD_WIDTH_M) > 20.0:
        flags.append(f"WIDTH({record['measured_road_width_m']:.0f}m)")

    record["flags"] = ",".join(flags)
    return record, profiles


# ==============================================================================
# SEAWARD RELOCATION -- see RELOCATE_DROWNING for the decision and its cost
# ==============================================================================

def load_saved_interior(ext, domain: int) -> np.ndarray | None:
    """
    The interior array CASCADE actually initialises with.

    Deliberately the SAVED .npy, not a re-derived interior: the drown test has
    to run on the same bytes `HAT_road_placement_on_domains.py` and
    `HAT_road_setback_audit.py` read, or the three scripts can disagree about
    which domains drown. Paths come off the extractor rather than being
    hardcoded, so bumping VERSION does not silently point this at stale arrays.
    """
    p = Path(ext.TOPO_SAVE_PATH) / f"domain_{domain}_topography_{ext.TAG}.npy"
    return np.load(p) if p.is_file() else None


def drown_test(interior: np.ndarray, road_start: int) -> tuple | None:
    """
    bulldoze's width-drown test, transcribed. Returns (seaside, bayside, drowns).

    The rows tested are the NEIGHBOURS of the bulldozed band -- road_start - 1
    and road_end + 1 -- never the band itself. Interior values are decametres
    MHW and bulldoze compares `grid * dz` against the threshold in metres, so
    the * 10 is the model's, not a display convenience. Every cell counts,
    no-data included: the extractor stores no-data as the water sentinel and
    bulldoze reads the literal array.

    None means the placement is not testable -- bulldoze indexes road_end + 1
    with no bounds check, so that is an IndexError at t=0, not a drowning.
    """
    n = interior.shape[0]
    end = road_start + ROAD_WIDTH_CELLS
    if road_start < 0 or end + 1 >= n:
        return None
    wet = lambda row: float((row * 10.0 <= DROWN_THRESHOLD_M).mean())
    seaside = wet(interior[road_start - 1, :]) if road_start > 0 else 0.0
    bayside = wet(interior[end + 1, :])
    return seaside, bayside, bool(seaside > DROWN_PCT or bayside > DROWN_PCT)


def nearest_viable_seaward(interior: np.ndarray, road_start: int) -> int | None:
    """
    Largest road_start strictly seaward of the current one that does not drown.

    Scans landward-to-seaward and returns the first pass, so the road moves the
    shortest distance that survives initialisation. None means no row seaward of
    here is viable -- the domain keeps its measured setback and drowns.
    """
    for start in range(road_start - 1, -1, -1):
        t = drown_test(interior, start)
        if t is not None and not t[2]:
            return start
    return None


def relocate_drowning(ext, in_span: list[dict], cell: float) -> None:
    """
    Fill `setback_model_m` for every in-span domain, moving the drowned ones.

    Mutates the records in place. `setback_2009_m` is never touched -- the
    measurement stays recoverable from the _domains.csv, which is the whole
    point of doing this here rather than in the measurement.
    """
    for r in in_span:
        floored = float(r["setback_2009_floored_m"])
        r["setback_model_m"] = round(floored, 1)

        interior = load_saved_interior(ext, r["domain"])
        if interior is None:
            r["flags"] = ",".join(filter(None, [r["flags"], "NO_SAVED_TOPO"]))
            continue

        start = int(floored / cell)
        t = drown_test(interior, start)
        if t is None:
            # Past the end of the array: an IndexError mid-run, not a drowning.
            # Left alone here; the audit script is where that is adjudicated.
            r["flags"] = ",".join(filter(None, [r["flags"], "OVERRUN"]))
            continue

        sea, bay, drowns = t
        r["drowns_at_init"] = "yes" if drowns else "no"
        if not drowns or not RELOCATE_DROWNING:
            continue

        best = nearest_viable_seaward(interior, start)
        if best is None:
            r["flags"] = ",".join(filter(None, [
                r["flags"],
                f"DROWNS_NO_VIABLE_ROW(sea{sea:.0%},bay{bay:.0%})"]))
            continue

        moved_m = (start - best) * cell
        r["setback_model_m"] = round(best * cell, 1)
        r["relocated_seaward_m"] = round(moved_m, 1)
        note = [f"MOVED_SEAWARD({moved_m:.0f}m,sea{sea:.0%},bay{bay:.0%})"]

        # A move inside the domain's own per-profile spread is a re-pick of the
        # alongshore statistic. A move beyond it is a position no profile
        # showed, which is a different claim and has to say so.
        if (np.isfinite(r["setback_min_m"])
                and r["setback_model_m"] < r["setback_min_m"]):
            note.append(f"BEYOND_MEASURED(min {r['setback_min_m']:.0f}m)")
        r["flags"] = ",".join(filter(None, [r["flags"]] + note))


# ==============================================================================
# EXISTING SAME-YEAR SETBACKS + OFFSET VALIDATION
# ==============================================================================

def read_two_row_csv(path: Path) -> dict[int, float]:
    """Read a 2-row (GIS IDs, values) CASCADE forcing file."""
    if not path.is_file():
        return {}
    raw = np.loadtxt(path, delimiter=",")
    if raw.ndim != 2 or raw.shape[0] != 2:
        raise ValueError(f"{path}: expected 2 rows, got shape {raw.shape}")
    return {int(k): float(v) for k, v in zip(raw[0], raw[1])}


def write_two_row_csv(path: Path, values: dict[int, float]) -> None:
    """Write the 2-row format load_padded_series expects."""
    path.parent.mkdir(parents=True, exist_ok=True)
    ids = sorted(values)
    with open(path, "w", newline="") as f:
        f.write(",".join(f"{i}.000" for i in ids) + "\n")
        f.write(",".join(f"{values[i]:.3f}" for i in ids) + "\n")


def load_offsets(year: int) -> np.ndarray | None:
    """Per-domain dune-line offset, seaward-positive metres. 90 rows expected."""
    path = Path(str(OFFSET_FMT).format(year=year))
    if not path.is_file():
        return None
    values = np.loadtxt(path, delimiter=",", skiprows=1)
    values = np.asarray(values, dtype=float).reshape(-1)
    if values.size == 120:            # 15 + 90 + 15 padding
        values = values[15:105]
    return values if values.size == len(DOMAINS) else None


# ==============================================================================
# MAIN
# ==============================================================================

def main() -> None:
    ext = load_extractor()
    windows = json.load(open(ext.WINDOW_JSON))

    print("=" * 84)
    print("NC-12 road offset measured from the extracted dune start")
    print("=" * 84)
    print(f"  extractor      : {EXTRACTOR.name}  "
          f"(ALONGSHORE_FLIP={ext.ALONGSHORE_FLIP}, STRAIGHTEN={ext.STRAIGHTEN})")
    print(f"  DEMs           : {ext.LOAD_PATH}")
    print(f"  picked windows : {ext.WINDOW_JSON.name} "
          f"({sum(1 for k in windows if k != '_meta')} domains)")
    print(f"  reference      : interior row 0 = dune crest + 1 cell")
    print(f"  output root    : {OUT_ROOT}")

    audit: dict[int, dict] = {}

    for year in YEARS:
        print(f"\n--- {year} " + "-" * 68)
        mask_dir = Path(str(MASK_DIR_FMT).format(year=year))
        if not mask_dir.is_dir():
            print(f"  [skip] no mask folder: {mask_dir}")
            continue

        records, all_profiles = [], []
        for domain in DOMAINS:
            record, profiles = measure_domain(ext, domain, year, windows)
            records.append(record)
            all_profiles.extend(profiles)

        with_road = [r for r in records if r["n_road_profiles"] > 0]
        if not with_road:
            print("  [skip] no domain carried road")
            continue

        # Measured everywhere, written only within ROAD_SPAN. Flag the excluded
        # ones in the diagnostics so the gap between "has road" and "is forced"
        # is on the record.
        in_span = [r for r in with_road
                   if ROAD_SPAN[0] <= r["domain"] <= ROAD_SPAN[1]]
        excluded = [r for r in with_road if r not in in_span]
        for r in excluded:
            r["flags"] = ",".join(filter(None, [r["flags"], "EXCLUDED_FROM_SPAN"]))
        if not in_span:
            print(f"  [skip] no domain with road inside ROAD_SPAN {ROAD_SPAN}")
            continue

        road_ids = [r["domain"] for r in in_span]
        first_gis, last_gis = min(road_ids), max(road_ids)
        if excluded:
            print(f"  excluded from span   : "
                  f"{[r['domain'] for r in excluded]} "
                  f"(measured and kept in _domains.csv, not written to the "
                  f"model-facing files)")

        # --- second reference frame: the existing same-year measurement ---
        existing = read_two_row_csv(Path(str(EXISTING_SETBACK_FMT).format(year=year)))
        for r in records:
            same = existing.get(r["domain"], np.nan)
            r["setback_legacy_m"] = same
            r["delta_vs_legacy_m"] = (
                round(same - r["setback_2009_m"], 1)
                if np.isfinite(same) and np.isfinite(r["setback_2009_m"])
                else np.nan)

        out_dir = OUT_ROOT / str(year)

        # --- seaward relocation of roadways that drown at initialisation ---
        # Runs on the FLOORED value, because that is what CASCADE would index
        # with, and fills setback_model_m for every in-span domain.
        relocate_drowning(ext, in_span, ext.CELL_SIZE_M)

        # --- model-facing files ---
        write_two_row_csv(
            out_dir / f"RoadSetback_{year}_dunestart.csv",
            {r["domain"]: r["setback_model_m"] for r in in_span},
        )
        elev = {r["domain"]: r["road_elev_mhw_median"] for r in in_span
                if np.isfinite(r["road_elev_mhw_median"])}
        write_two_row_csv(out_dir / f"RoadElevation_{year}_dunestart.csv", elev)

        # --- diagnostics: every measured domain, in-span or not ---
        write_csv(out_dir / f"RoadOffset_{year}_domains.csv", records)
        write_csv(out_dir / f"RoadOffset_{year}_profiles.csv",
                  [p for p in all_profiles
                   if ROAD_SPAN[0] <= p["domain"] <= ROAD_SPAN[1]])

        audit[year] = summarize(year, records, in_span, first_gis, last_gis)
        report(audit[year], records)

    if CONTROL_UNSTRAIGHTENED:
        run_control(ext)

    if audit:
        write_audit(audit, ext)
        print(f"\n[audit] {OUT_ROOT / 'RoadOffset_dunestart_audit.md'}")


def run_control(ext) -> None:
    """Re-measure with STRAIGHTEN = False so the frame effect is separable."""
    if not CONTROL_WINDOW_JSON.is_file():
        print(f"\n[control] skipped: no unstraightened picks at "
              f"{CONTROL_WINDOW_JSON}")
        return

    print(f"\n--- UNSTRAIGHTENED CONTROL " + "-" * 52)
    print(f"  picks: {CONTROL_WINDOW_JSON.name}")

    windows = json.load(open(CONTROL_WINDOW_JSON))
    saved = ext.STRAIGHTEN
    ext.STRAIGHTEN = False
    try:
        for year in YEARS:
            mask_dir = Path(str(MASK_DIR_FMT).format(year=year))
            if not mask_dir.is_dir():
                continue
            records = [measure_domain(ext, d, year, windows)[0] for d in DOMAINS]
            for r in records:
                r.pop("setback_legacy_m", None)
                r.pop("delta_vs_legacy_m", None)
            out = (OUT_ROOT / str(year)
                   / f"RoadOffset_{year}_domains{CONTROL_SUFFIX}.csv")
            write_csv(out, records)
            got = [r for r in records if r["n_road_profiles"] > 0]
            sb = np.array([r["setback_2009_m"] for r in got], dtype=float)
            print(f"  {year}: {len(got)} domains, setback median "
                  f"{np.nanmedian(sb):.0f} m, range {np.nanmin(sb):.0f} to "
                  f"{np.nanmax(sb):.0f} m  ->  {out.name}")
    finally:
        ext.STRAIGHTEN = saved


def write_csv(path: Path, rows: list[dict]) -> None:
    import csv
    if not rows:
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def summarize(year, records, with_road, first_gis, last_gis) -> dict:
    sb = np.array([r["setback_2009_m"] for r in with_road], dtype=float)
    delta = np.array([r["delta_vs_legacy_m"] for r in with_road], dtype=float)
    negative = [r for r in with_road if r["setback_2009_m"] < 0]

    offsets_year = load_offsets(year)
    offsets_2004 = load_offsets(2004)
    retreat = None
    corr = np.nan
    if offsets_year is not None and offsets_2004 is not None and year != 2004:
        retreat = offsets_year - offsets_2004     # seaward-positive: >0 = retreat
        # Does delta_vs_legacy actually behave like retreat? If the legacy file
        # and this one were measuring the same thing in different years, this
        # correlation would be strongly POSITIVE. It is not -- see the audit.
        by_domain = {d: retreat[d - 1] for d in DOMAINS}
        pairs = [(r["delta_vs_legacy_m"], by_domain[r["domain"]])
                 for r in with_road
                 if np.isfinite(r["delta_vs_legacy_m"])
                 and r["domain"] in by_domain]
        if len(pairs) > 2:
            a, b = np.asarray(pairs, dtype=float).T
            corr = float(np.corrcoef(a, b)[0, 1])

    return {
        "year": year,
        "first_gis": first_gis, "last_gis": last_gis,
        "n_with_road": len(with_road),
        "n_no_road": sum(1 for r in records if r["flags"].startswith("NO_ROAD")),
        "setback_min": float(np.nanmin(sb)), "setback_max": float(np.nanmax(sb)),
        "setback_median": float(np.nanmedian(sb)),
        "n_negative": len(negative),
        "negative_domains": [(r["domain"], r["setback_2009_m"]) for r in negative],
        "n_drowned": sum(1 for r in with_road if r["drowns_at_init"] == "yes"),
        "relocated": [(r["domain"], r["setback_2009_floored_m"],
                       r["setback_model_m"], r["relocated_seaward_m"],
                       "BEYOND_MEASURED" in r["flags"])
                      for r in with_road if r["relocated_seaward_m"]],
        "unfixable": [r["domain"] for r in with_road
                      if "DROWNS_NO_VIABLE_ROW" in r["flags"]],
        "delta_median": float(np.nanmedian(delta)) if np.isfinite(delta).any() else np.nan,
        "delta_min": float(np.nanmin(delta)) if np.isfinite(delta).any() else np.nan,
        "delta_max": float(np.nanmax(delta)) if np.isfinite(delta).any() else np.nan,
        "retreat_median": float(np.median(retreat)) if retreat is not None else np.nan,
        "delta_retreat_corr": corr,
        "flagged": [(r["domain"], r["flags"]) for r in records if r["flags"]],
    }


def report(s: dict, records: list[dict]) -> None:
    print(f"  road span            : GIS {s['first_gis']}-{s['last_gis']} "
          f"({s['n_with_road']} domains with road, {s['n_no_road']} without)")
    print(f"  setback vs 2009 dune : median {s['setback_median']:.0f} m, "
          f"range {s['setback_min']:.0f} to {s['setback_max']:.0f} m")
    print(f"  delta vs LEGACY file : median {s['delta_median']:.0f} m, "
          f"range {s['delta_min']:.0f} to {s['delta_max']:.0f} m")
    if np.isfinite(s["retreat_median"]):
        print(f"  offset-derived retreat: median {s['retreat_median']:.0f} m")
        print(f"  corr(delta, retreat) : {s['delta_retreat_corr']:+.3f}  "
              f"<- would be near +1 if the two files measured the same feature")
    print(f"  NEGATIVE (floored)   : {s['n_negative']} domain(s)"
          + (f" -> {[d for d, _ in s['negative_domains']]}" if s["n_negative"] else ""))
    print(f"  DROWNS at init       : {s['n_drowned']} domain(s)"
          + ("" if RELOCATE_DROWNING else "   [RELOCATE_DROWNING=False]"))
    if s["relocated"]:
        print(f"  moved seaward to the nearest viable row:")
        for d, was, now, moved, beyond in s["relocated"]:
            print(f"      GIS {d:>2}: {was:>5.0f} -> {now:>5.0f} m "
                  f"({moved:>4.0f} m seaward)"
                  + ("   BEYOND the measured per-profile range" if beyond
                     else "   within the measured range"))
    if s["unfixable"]:
        print(f"  [warn] no viable row seaward, still drowning: {s['unfixable']}")


def write_audit(audit: dict, ext) -> None:
    OUT_ROOT.mkdir(parents=True, exist_ok=True)
    path = OUT_ROOT / "RoadOffset_dunestart_audit.md"

    # Does anything actually drown on THIS extraction? The no-data explanation
    # below is a measurement taken on 2009_v4, not a recomputation, so it must
    # not be printed under a version where nothing drowns -- it would assert a
    # cell count that was never measured on the arrays in front of it. The
    # gap-filled DEM (2009_v5) is exactly that case: 0 drowns, both years.
    total_drowned = sum(s["n_drowned"] for s in audit.values())

    lines = [
        "# NC-12 road offset measured from the extracted dune start",
        "",
        f"Generated {datetime.now().isoformat(timespec='seconds')} by "
        "`HAT_road_offset_from_dune_start.py`.",
        "",
        "| | |",
        "|---|---|",
        f"| Extractor | `{EXTRACTOR.name}` "
        f"(ALONGSHORE_FLIP={ext.ALONGSHORE_FLIP}, STRAIGHTEN={ext.STRAIGHTEN}) |",
        f"| DEMs | `{ext.LOAD_PATH}` |",
        f"| Picked windows | `{ext.WINDOW_JSON.name}` |",
        f"| Reference | interior row 0 = picked dune crest + 1 cell |",
        f"| Road reference | seaward-most road cell per profile |",
        f"| Alongshore collapse | median over profiles with road |",
        f"| Elevation | median of non-water road cells, m MHW |",
        f"| Assumed road width | {ASSUMED_ROAD_WIDTH_M:.0f} m |",
        "",
        "`road_setback` is metres landward of interior row 0 "
        "(`roadway_manager.py:99`). The road mask is put through the same "
        "orient / alongshore-flip / shear / trim chain as the topography, so "
        "the road and the dune line are compared on identical ground.",
        "",
        "## Negative setbacks are floored, not measured away",
        "",
        "A negative setback means the road falls SEAWARD of interior row 0. It "
        "cannot go into the model-facing CSV: `int(-50/10) = -5` and "
        "`xyz_interior_grid[-5:-3, :]` indexes from the landward end, so the "
        "road would be bulldozed into the bay with no error. `0.0` is not a "
        "\"no road\" sentinel either -- `build_roadway_management_on` uses the "
        "road span, not the value. The model-facing file therefore carries "
        "`max(setback, 0)`; the true signed value is in "
        "`RoadOffset_<year>_domains.csv` under `setback_2009_m`.",
        "",
        "## Roadways that drown at initialisation are moved seaward",
        "",
        "A roadway whose flanking rows are >20% water at t=0 is width-drowned "
        "by `roadway_manager.bulldoze` on the first call. `RoadwayManager` sets "
        "`_drown_break`, `cascade.py` never calls `update()` again, and that "
        "domain spends the whole hindcast as an **unmanaged barrier wearing a "
        "road label** -- no overwash removal, no dune rebuilding, no "
        "relocation. To keep the domain managed, the model-facing setback is "
        "moved to the nearest row seaward that passes bulldoze's own test.",
        "",
        "The test transcribed here is bulldoze's: the rows checked are the "
        "NEIGHBOURS of the bulldozed band (`road_start - 1`, `road_end + 1`), "
        "never the band itself, and every cell counts. There is **no cap** on "
        "the distance moved -- the nearest viable row is taken however far it "
        "is. `setback_2009_m` is never altered; only "
        "`RoadSetback_<year>_dunestart.csv` carries the moved value, and every "
        "moved domain is flagged `MOVED_SEAWARD`.",
        "",
        "### The assumption this rests on",
        "",
    ]

    if total_drowned == 0:
        lines += [
            f"**Nothing drowns at initialisation on {ext.RUN_NAME}** -- 0 "
            "domains in either year -- so no setback was moved and the rest of "
            "this section does not apply to this run.",
            "",
            "That is the point of the gap-filled DEM. On `2009_v4` three "
            "domains per year drowned, and the audit for that version showed "
            "they drowned on **LiDAR coverage gaps, not measured water**: "
            "across the six flanking rows failing at GIS 78/79/80 there were "
            "106 wet cells, of which 105 had never been surveyed and 1 was "
            "genuinely wet. The extractor writes no-data back as the water "
            "sentinel because Barrier3D has no representation for \"unknown\", "
            "so unsurveyed ground read as ocean and drowned the roadway. "
            "Filling those holes from the 2014 NOAA Post-Sandy DEM removes the "
            "cause, and the drown count goes to zero.",
            "",
            "Those figures are quoted from the `2009_v4` audit and were "
            "measured there, not recomputed here -- see "
            "`dunestart_offset_ARCHIVE_2009_v4/`.",
            "",
        ]
    else:
        lines += [
            f"On the {ext.RUN_NAME} topography these domains **do not drown on "
            "measured water -- they drown on LiDAR coverage gaps.** Across the "
            "six flanking rows that fail at GIS 78/79/80 there are 106 wet "
            "cells, of which **105 were never surveyed and 1 is genuinely "
            "measured wet.** The extractor writes no-data back as the water "
            "sentinel because Barrier3D has no representation for \"unknown\", "
            "and the drown test counts every cell -- deliberately, so that what "
            "is reported is what the run does.",
            "",
            "NOTE: that 106/105/1 split was measured on `2009_v4` and is not "
            "recomputed per run. Confirm it still holds before quoting it on a "
            "different extraction.",
            "",
            "So this is **not** \"the road was in water, we moved it out\". It "
            "is \"the 2009 DEM has no data there, CASCADE reads no-data as "
            "water, so the road is moved onto surveyed ground to keep the "
            "domain managed\". Any managed-vs-unmanaged result at GIS 78-80 "
            "needs that sentence attached to it.",
            "",
        ]

    lines += [
        "A move that lands inside the domain's own per-profile spread is a "
        "re-pick of the alongshore statistic. A move beyond it is a position no "
        "profile showed -- a different claim, flagged `BEYOND_MEASURED`.",
        "",
        "## `delta_vs_legacy_m` is NOT the dune-line retreat",
        "",
        "It is tempting to read the difference against the legacy "
        "`RoadSetback_<year>.csv` as the retreat between `<year>` and 2009. It "
        "is not, and the reported `corr(delta, retreat)` shows it: a like-for-"
        "like pair of measurements taken in two different years would correlate "
        "near +1 with the offset-derived retreat, and the measured correlation "
        "is strongly negative.",
        "",
        "Three reasons the two files are not commensurable, from "
        "`HAT_setback_from_lines.py` (retired 2026-08-17, git blob "
        "`c37ec03e`):",
        "",
        "1. **Different dune feature.** The legacy file computes "
        "`setback_i = road_cell_i - dune_cell_i` against a digitized same-year "
        "*dune-line geojson*. This script measures against the *DEM dune crest* "
        "found inside the picked search window. Those are different features, "
        "not the same feature at two dates.",
        "2. **Different frame.** The legacy measurement is taken with both lines "
        "\"raw, ocean-first\" (`HAT_setback_from_lines.py:254`) -- i.e. "
        "unstraightened, so it still carries the obliquity smear this script "
        "removes.",
        "3. **The legacy file is already floored.** It prints "
        "`FLOORED(x) raw setback was negative`, so at the domains where the two "
        "disagree most, the legacy value is itself a clamp, not a measurement.",
        "",
        "So `delta_vs_legacy_m` is a *migration diagnostic* -- how much each "
        "domain's forcing changes if you adopt this method -- and nothing more. "
        "Closing the same-year comparison properly needs a same-year dune crest "
        "measured the same way, which needs a same-year DEM.",
        "",
        "## The road reference is a buffered mask edge, not the road",
        "",
        "`road_setback` is measured to the seaward-most cell of the "
        "**rasterized mask**, which `HAT_rasterize_road_to_domains.py` burns with "
        "`ROAD_BUFFER_M = 6.0` and `ALL_TOUCHED = True` -- roughly a 24 m mask "
        "for an ~8 m road. So the reference is not NC-12's centerline, and not "
        "its pavement edge either. `HAT_road_buffer_bias.py` measures the "
        "difference against the source geojson, in the original raster frame "
        "(orient / flip / shear / trim all cancel out of it, since interior row 0 "
        "is common to both):",
        "",
        "| | 1984 | 2004 |",
        "|---|---|---|",
        "| median bias | -6.8 m | -6.8 m |",
        "| p10 / p90 | -10.8 / -2.7 m | -10.8 / -2.8 m |",
        "",
        "Negative means the mask edge sits **seaward** of the centerline, so "
        "every setback here runs about 7 m smaller than a centerline measurement "
        "would.",
        "",
        "### That is close to right, and correcting it would be worse",
        "",
        f"`bulldoze` lays a {ASSUMED_ROAD_WIDTH_M:.0f} m block LANDWARD from "
        "`road_start`, so a block centred on the real road wants "
        f"`road_start = centerline - {ASSUMED_ROAD_WIDTH_M / 2:.0f} m`. The "
        "buffer supplies `centerline - 6.9 m`, leaving a residual misplacement "
        "of **3.1 m -- about a third of a cell**, and an order of magnitude "
        "below the 20-40 m p90 placement error the alongshore collapse to one "
        "scalar already causes (`HAT_method_comparison_figures.py`).",
        "",
        "It is deliberately left alone. Per-profile setbacks are integer cell "
        "differences times the cell size, so domain medians land almost exactly "
        "on cell boundaries, and `bulldoze` truncates with `int()`. Applying the "
        "3.1 m adjustment moves `road_start` a **full cell (10 m) seaward in 83% "
        "of 1984 domains and 90% of 2004 domains** -- trading a 3 m error for a "
        "7 m error in the opposite direction. The buffer is accidentally doing "
        "very nearly the right job; the quantization is what would bite.",
        "",
        "Two independent paths agree on the number: the raster-frame column "
        "comparison above, and the `road_x`/`road_y` columns in "
        "`RoadOffset_<year>_profiles.csv`, which are produced by inverting the "
        "full index chain back to map coordinates and whose shapely distance to "
        "the same geojson is 6.62-6.68 m median (max 12.6 m, no point beyond "
        "60 m). That agreement also validates the index inversion -- an "
        "alongshore-flip error would mirror points within each 500 m domain and "
        "show up as distances in the hundreds of metres.",
        "",
        "D8 is the one extreme outlier (-353 m in 1984). That is the Buxton bend, "
        "where NC-12 runs parallel to the raster rows so a \"seaward-most cell\" "
        "has no physical reading -- the same reason D8 is already "
        "`EXCLUDED_FROM_SPAN`. Its appearance here is a consistency check "
        "passing, not a new problem.",
        "",
    ]
    for year, s in audit.items():
        lines += [
            f"## {year}",
            "",
            f"- Road span: GIS {s['first_gis']}-{s['last_gis']} "
            f"({s['n_with_road']} with road, {s['n_no_road']} without)",
            f"- Setback vs 2009 dune: median {s['setback_median']:.0f} m, "
            f"range {s['setback_min']:.0f} to {s['setback_max']:.0f} m",
            f"- Delta vs the LEGACY `RoadSetback_{year}.csv`: median "
            f"{s['delta_median']:.0f} m, range {s['delta_min']:.0f} to "
            f"{s['delta_max']:.0f} m",
        ]
        if np.isfinite(s["retreat_median"]):
            lines += [
                f"- Offset-derived dune-line retreat {year}->2004: median "
                f"{s['retreat_median']:.0f} m",
                f"- **corr(delta, retreat) = {s['delta_retreat_corr']:+.3f}** "
                f"-- see the caveat below.",
            ]
        lines.append(f"- NEGATIVE, floored to 0: {s['n_negative']} domain(s)")
        if s["negative_domains"]:
            lines += ["", "| GIS | true setback (m) |", "|---:|---:|"]
            lines += [f"| {d} | {v:.0f} |" for d, v in s["negative_domains"]]
        lines.append(f"- DROWNS at initialisation: {s['n_drowned']} domain(s)"
                     + ("" if RELOCATE_DROWNING
                        else "  (RELOCATE_DROWNING = False, none moved)"))
        if s["relocated"]:
            lines += [
                "",
                "| GIS | measured (m) | written (m) | moved seaward (m) | "
                "inside the measured per-profile range? |",
                "|---:|---:|---:|---:|---|",
            ]
            lines += [
                f"| {d} | {was:.0f} | {now:.0f} | {moved:.0f} | "
                + ("**no -- beyond it**" if beyond else "yes") + " |"
                for d, was, now, moved, beyond in s["relocated"]
            ]
        if s["unfixable"]:
            lines.append(f"- No viable row seaward, still drowning: "
                         f"{s['unfixable']}")
        if s["flagged"]:
            lines += ["", "### Flags", "", "| GIS | flags |", "|---:|---|"]
            lines += [f"| {d} | {f} |" for d, f in s["flagged"]]
        lines.append("")

    path.write_text("\n".join(lines), encoding="utf-8")


if __name__ == "__main__":
    main()
