# ==============================================================================
# HAT_dune_topo_extractor.py
#
# Hatteras CASCADE dune & interior topography extractor with per-domain,
# interactively selected dune search windows.
#
# INPUT  : domain_#.npy, shape = (alongshore_rows, cross_shore_cols),
#          elevation in m NAVD88
# OUTPUT : interior topography and dune height arrays, in decameters (dam)
#          + a JSON of the per-domain dune search windows (re-runnable)
#
# WORKFLOW
#   1. MODE = "pick"          -> step through domains, drag a dune search
#                                window on the profile stack, saves to JSON
#                                after every domain (safe to quit and resume)
#   2. MODE = "run"           -> extract using the saved JSON windows
#   3. MODE = "pick_and_run"  -> both in one pass
#
# PICKER LAYOUT / KEYS
#   The picker draws the domain with the OCEAN AT THE BOTTOM: cross-shore runs
#   vertically (cell 0 = ocean, landward upward), alongshore runs left-right.
#   Left panel = elevation map, right panel = profile stack, shared y-axis.
#   This is display only; i0/i1 are still cross-shore indices from the ocean.
#
#   NC-12 is drawn on both panels (v4): filled road cells and a dashed centre
#   line on the map, a shaded cross-shore envelope on the profile stack, one
#   colour per road vintage. It is there to stop the window being dragged onto
#   the road embankment, which a dune-crest argmax will happily lock onto. The
#   road constrains nothing in code -- see ROAD OVERLAY in CONFIG.
#
#   drag vertically on either panel : set the cross-shore search window
#   enter / close window            : accept current window
#   r                               : reset to the default window
#   s                               : skip this domain (use DEFAULT_WINDOW_PX)
#   q / esc                         : quit picking, keep everything saved so far
# ==============================================================================

from __future__ import annotations

import csv
import json
import os
import re
from datetime import datetime
from pathlib import Path

import numpy as np

import matplotlib
try:
    # needed for a real, blocking, interactive picker window. Only force it if
    # tkinter actually exists, otherwise matplotlib fails later with a confusing
    # error at figure-creation time instead of here.
    import tkinter  # noqa: F401
    matplotlib.use("TkAgg")
except Exception:
    print("[warn] TkAgg unavailable; the picker needs an interactive backend.")
import matplotlib.pyplot as plt
from matplotlib.colors import FuncNorm, ListedColormap
from matplotlib.transforms import blended_transform_factory
from matplotlib.widgets import SpanSelector

plt.rcParams["font.size"] = 11

# ==============================================================================
# CONFIG
# ==============================================================================

# --- MODE ---------------------------------------------------------------
# v4 RE-PICKS every window with the road drawn on the picker (see ROAD OVERLAY).
# The v3 windows were drawn without knowing where NC-12 sits, so a window could
# sit landward of the road and take the road embankment for a dune crest. This
# pass is an ADJUSTMENT, not a blind redraw: the picker opens on each domain's
# saved v4 window (seeded from v3), so "r" resets to the v3 pick and accepting
# without dragging keeps it.
MODE = "pick_and_run"      # "pick" | "run" | "pick_and_run"
# Domains to process. This filters BOTH passes -- picking AND running -- so it
# is the modelled set, not just a picking subset.
#
#   list(range(1, 91))      the 90 domains CASCADE runs (D1 = Cape Point ->
#                           D90 = Rodanthe). The DEM folder holds 131; 91-131
#                           are north of the study area and are not modelled,
#                           so picking or extracting them is wasted work.
#   None                    every domain_*.npy found
#   [1, 5, 11, 33, 67, 74]  a trial set. Worth doing before committing to 90
#                           hand-picks: 33 (25.6 deg) and 67 (22.8 deg) are the
#                           worst obliquity on the island and two of the four
#                           SCATTER flags; 11 (8.0 deg) is the domain the road
#                           setbacks kept failing on; 74 (5.7 deg) is a
#                           near-square control where panels 1 and 2 should look
#                           almost identical; 1 and 5 are Cape Point, where a
#                           LINEAR fit to the shoreline is most likely to fail.
PICK_DOMAINS = list(range(1, 91))
# True for the v4 road re-pick. The v4 picks file is SEEDED from v3, so every
# domain is already present and False would skip all 90 -- nothing would be
# re-picked. Set back to False once the re-pick is finished, so a later run in
# "pick_and_run" resumes instead of starting over.
REPICK_EXISTING = True     # False = skip domains already present in the JSON
SAVE_QC_FIGS = True        # per-domain dune-detection QC figure
SAVE_COMPARISON_FIGS = True  # per-domain raw GIS vs processed CASCADE input figure
SAVE_SETTINGS_SHEET = True   # per-domain settings/results sheet (csv + xlsx)

# --- PATHS --------------------------------------------------------------
# Everything for one settings variant lands in ONE run folder, so comparing
# versions means comparing two directories:
#
#   data\hatteras_init\dune_topo\
#       picks\
#           HAT_dune_search_windows_2009_pea_hatteras.json  <- your picks (see below)
#       2009_v1\
#           RUN_MANIFEST.txt                     <- every setting that made this folder
#           HAT_dune_topo_settings_2009_v1.xlsx  <- per-domain sheet (+ .csv)
#           HAT_dune_topo_summary_2009_v1.png    <- all domains on one page
#           topography\  domain_7_topography_2009.npy   <- CASCADE reads these two
#           dunes\       domain_7_dune_2009.npy
#           figures\
#               gis_vs_processed\  domain_007_gis_vs_processed.png
#               qc\                domain_007_qc.png
#       2009_v2\   ... same shape, nothing shared, nothing overwritten
#
# NOTE: topography\ and dunes\ moved. Point the hindcast runner at
#       RUN_DIR\topography and RUN_DIR\dunes for whichever version you're using.
# v3 = v2 settings + ALONGSHORE_FLIP = True (see GEOMETRY below). Bumped rather
# than reused so 2009_v2 survives as the unflipped reference to diff against.
# v4 = v3 settings + the NC-12 overlay and a full re-pick of the dune windows
# with the road visible. No processing setting changed, so a v4 run in "run"
# mode on the v3 windows reproduces the v3 arrays byte for byte -- what makes
# v4 different is the WINDOWS, which is exactly why it gets its own folder and
# its own picks file rather than overwriting v3.
VERSION = "v4"             # bump this per settings variant, like Lexi's sheet
DEM_YEAR = "2009"
DEM_NAME = "2009_pea_hatteras"
RUN_NAME = f"{DEM_YEAR}_{VERSION}"

PROJECT_ROOT = Path(r"C:\Users\hanna\PycharmProjects\CASCADE")
INIT_ROOT = PROJECT_ROOT / "data" / "hatteras_init"
# These two moved when data\hatteras_init was reorganized into the numbered
# 1-barrier3d-domains \ 2-brie-offset \ ... tree. They previously read
# INIT_ROOT/"elevations"/DEM_NAME and INIT_ROOT/"dune_topo", neither of which
# exists any more, so the script could not find its own inputs or the 2009_v2
# run folder it wrote.
LOAD_PATH = (INIT_ROOT / "1-barrier3d-domains" / f"{DEM_YEAR}-raw"
             / f"{DEM_YEAR}-npy-arrays" / DEM_NAME)

DUNE_TOPO_ROOT = INIT_ROOT / "1-barrier3d-domains" / f"{DEM_YEAR}-dune-topo"
RUN_DIR = DUNE_TOPO_ROOT / RUN_NAME
TOPO_SAVE_PATH = RUN_DIR / "topography"
DUNE_SAVE_PATH = RUN_DIR / "dunes"
FIG_DIR_CMP = RUN_DIR / "figures" / "gis_vs_processed"
FIG_DIR_QC = RUN_DIR / "figures" / "qc"
SHEET_SAVE_PATH = RUN_DIR / f"HAT_dune_topo_settings_{RUN_NAME}"   # .csv + .xlsx
SUMMARY_FIG_PATH = RUN_DIR / f"HAT_dune_topo_summary_{RUN_NAME}.png"
ISLAND_FIG_PATH = RUN_DIR / f"HAT_dune_topo_island_offsets_{RUN_NAME}.png"
PLAN_FIG_STEM = f"HAT_dune_topo_island_planview_{RUN_NAME}"  # + _{year}.png
MANIFEST_PATH = RUN_DIR / "RUN_MANIFEST.txt"

# Picks live OUTSIDE the run folder. They are the only artifact here you cannot
# regenerate, and they describe where the dune sits in the DEM -- not which
# settings variant you're testing -- so v2, v3... reuse them by default. Set
# PICK_SET = RUN_NAME instead if you want a version to carry its own picks.
# The picks are FRAME-DEPENDENT. A window picked on a straightened array is a
# valid index range on an unstraightened one -- it just points at different
# cells -- so the two frames cannot share a file. save_windows() stamps
# "straightened" on every entry and the run pass refuses on a mismatch, but it
# writes back to WINDOW_JSON on every domain: pointing this at the v1 set would
# overwrite v1's unstraightened windows as you re-pick, and 2009_v1 would stop
# being reproducible.
# Set this by hand to match STRAIGHTEN below (it is defined further down, so it
# cannot be referenced here):
#     STRAIGHTEN = True   ->  f"{DEM_NAME}_straight"
#     STRAIGHTEN = False  ->  DEM_NAME
#
# v4 CARRIES ITS OWN PICKS -- PICK_SET = RUN_NAME, not the shared straight set.
# This is the case the warning above describes. v4 re-picks all 90 windows, and
# save_windows() writes back to WINDOW_JSON after EVERY domain, so pointing this
# at f"{DEM_NAME}_straight" would destroy the v3 picks as you worked and 2009_v3
# would stop being reproducible. The v4 file was seeded by copying the v3 one, so
# every domain opens on its v3 window and an unchanged domain stays unchanged:
#
#   picks\HAT_dune_search_windows_2009_pea_hatteras_straight.json   v1-v3, FROZEN
#   picks\HAT_dune_search_windows_2009_v4.json                      v4, re-picked
PICK_SET = RUN_NAME
PICKS_DIR = DUNE_TOPO_ROOT / "picks"
WINDOW_JSON = PICKS_DIR / f"HAT_dune_search_windows_{PICK_SET}.json"

TAG = DEM_YEAR             # appended to the .npy filenames CASCADE reads

# --- ISLAND OFFSETS -----------------------------------------------------
# Measured per-domain dune offsets used to place domains in a common cross-shore
# frame. One value per GIS domain, header row = year.
#
# CONVENTION (both inferred from the data -- change if your pipeline says otherwise):
#   OFFSET_ROW_ORDER = "D1_first"    row 0 of the CSV is domain 1 (Cape Point).
#       Check: this puts the largest raw_offset (6301 m) at Cape Point, which is what
#       a seaward-protruding headland should look like. Reversed puts the max at
#       Pea Island instead.
#   OFFSET_SEAWARD_POSITIVE = True   larger value = further seaward.
#       Check: 75/90 domains go negative 1984->2004, mean -2.07 m/yr. Seaward-
#       positive reads that as island-wide retreat at ~2 m/yr (right for Hatteras);
#       the other sign reads it as island-wide accretion (wrong).
SAVE_ISLAND_FIG = True
OFFSET_DIR = INIT_ROOT / "2-brie-offset"
OFFSET_FILES = {
    1984: OFFSET_DIR / "hindcast_1984" / "Island_Dune_Offsets_1984_CASCADE_Input.csv",
    2004: OFFSET_DIR / "hindcast_2004" / "Island_Dune_Offsets_2004_CASCADE_Input.csv",
}
CELL_SIZE_M = 10.0             # DEM cell size, cross-shore and alongshore (DAM_TO_M)
NUM_REAL_DOMAINS = 90
N_BUFFER_DOMAINS = 15          # raw_offset files may be padded to 15+90+15 = 120 rows
OFFSET_COLUMN = 0              # for multi-year raw_offset files, which column to read

# Plan-view canvas, reproducing plot_initialization_poster_no_border.py exactly:
#   offset_cells = round(offset_m / 10); each domain's topo row 0 (ocean side)
#   lands on canvas row = offset_cells; alongshore flipped with np.fliplr.
SAVE_ISLAND_PLAN_FIG = True
# ISLAND_FLIP_ALONGSHORE was REMOVED (2026-08-17). It applied np.fliplr to each
# domain INSIDE the placement loop, so it did not mirror the island -- it reversed
# the 50 cells of every 500 m block against the ascending domain order. That was a
# workaround for the source arrays having the within-domain alongshore order
# backwards, which ALONGSHORE_FLIP now fixes at load. Keeping both double-flipped:
# measured seam/inner discontinuity ratio on the plotted canvas was
#
#     v2 arrays + flip  1.97 (right)      v2 arrays, no flip  21.15
#     v3 arrays + flip 21.15 (WRONG)      v3 arrays, no flip   1.97 (right)
#
# so v2's plan view was only correct because two errors cancelled. The flip is
# gone rather than defaulted False, per the decision to keep this code lean.
# CONSEQUENCE: plotting a PRE-CORRECTION run (2009_v2 or earlier) through this
# function now sawtooths, and _assert_alongshore_continuity below will say so.
# To regenerate a v2 plan view, re-extract it with ALONGSHORE_FLIP = True instead.
ISLAND_INCLUDE_DUNE = True     # write the dune crest into canvas row raw_offset-1
ISLAND_ELEV_MIN_M = -1.0       # poster value. Raising this flattens land contrast:
ISLAND_ELEV_MAX_M = 4.0        #   land maps to ramp 0.35-1.0, so max=4 puts a 2 m
                               #   cell at tan, max=5 at pale yellow-green.
ISLAND_SEA_LEVEL_POS = 0.35    # colormap position for 0 m, as in the poster
ISLAND_SENTINEL_AS_OCEAN = False  # False = poster behaviour: sentinel (-3 m) water
                                  #   cells clip to the dark navy bottom of terrain,
                                  #   so the model's cross-shore extent stays visible
                                  #   and only outside-canvas NaN is light blue.
                                  # True = sentinel also renders light blue.
ISLAND_OCEAN_COLOR = "#b0cfe8"    # poster set_bad colour
# Display only, does not touch the .npy files. ONE FIGURE PER YEAR PER MODE, so a
# single run gives you both versions to compare:
#   "trimmed" = each domain only as tall as its own island, as stored
#   "padded"  = every domain given the same cross-shore extent (ISLAND_PAD_ROWS),
#               sentinel-filled landward, so the water behind each domain shows as
#               the navy wedge in the poster figure
# If TRIM_INTERIOR_ROWS = False the arrays are already padded and the two match.
ISLAND_CROSS_SHORE_MODES = ["trimmed", "padded"]
ISLAND_PAD_ROWS = 200               # cells of constant cross-shore extent when padded.
                                    # 200 = 2000 m, equal to TOPO_ROWS, so nothing is
                                    #       cropped (the poster look)
                                    # 100 = 1000 m, tighter; the script warns if that
                                    #       crops real land off the sound side
OFFSET_ROW_ORDER = "D1_first"  # "D1_first" | "D90_first"
OFFSET_SEAWARD_POSITIVE = True

# --- ISLAND SECTIONS ----------------------------------------------------
# D1 = Cape Point (south) -> D90 = near Rodanthe (north). Labels the sheet, the
# per-domain figures and the summary figure. Set to [] to disable.
SECTIONS = [
    ((1, 6),   "Cape Point"),
    ((7, 8),   "Buxton"),
    ((9, 20),  "Buxton-Avon"),
    ((21, 31), "Avon"),
    ((32, 67), "Avon-Tri-Village / Wimble Shoals"),
    ((68, 83), "Tri-Village"),
    ((84, 90), "Pea Island / N. Rodanthe"),
]

# --- DATUMS / THRESHOLDS ------------------------------------------------
MHW_M = 0.36               # m NAVD88
BERM_ELEV_NAVD_M = 1.70    # m NAVD88 (matches the 1.7 m storm/collision threshold)
BEACH_START_THR_M = 0.50   # m MHW-relative, strict '>' comparison
WATER_CLAMP_M = -3.0       # m MHW-relative; below this -> sentinel.
                           #   -3.0 keeps back-barrier marsh cells (Lexi's v3 edit)
                           #   -1.0 was the original dune_topo_extractor_from_GIS behavior
SENTINEL_WATER_M = -3.0    # m MHW-relative
MIN_DUNE_H_M = 0.1         # m, floor on dune height above berm

# --- NO DATA IS NOT WATER -----------------------------------------------
# The source LiDAR carries -10.0 m NAVD88 where it has no return. Clamping used
# to fold that into SENTINEL_WATER_M, so a cell the survey never saw became
# indistinguishable from a cell measured below MHW. That conflation is not
# cosmetic: roadway_manager.bulldoze drowns a roadway when >20% of the cells
# BORDERING it sit at or below 0 m MHW, and a no-data cell satisfies that test.
# In GIS 78/79/80 the row landward of NC-12 is 17-25 no-data cells and ZERO
# genuinely wet ones, so all three roadways "width-drowned" at t=0 on the
# strength of missing survey coverage. The giveaway is where the data stops:
# those profiles end while the ground is still 0.5-0.7 m ABOVE MHW, whereas
# their neighbours grade down through zero the way a real sound margin does.
#
# So no-data is tracked separately and written to its own array. The topography
# CASCADE reads is UNCHANGED -- no-data still lands on SENTINEL_WATER_M there,
# because Barrier3D has no representation for "unknown" and inventing an
# elevation would be worse. What changes is that the information now survives,
# so any consumer that cares can ask.
#
#   <stem>_nodata_<TAG>.npy   bool, same shape as the topography array
#                             True = this cell was never surveyed
RAW_NODATA_MAX_NAVD = -9.0   # raw NoData is exactly -10.0 m NAVD88
NODATA_SENTINEL_M = -99.0    # internal only; never written to the topography
WRITE_NODATA_MASK = True

# --- GEOMETRY -----------------------------------------------------------
TOPO_ROWS = 200            # max interior rows written
ALONG_COLS = 50            # alongshore profiles per domain (500 m / 10 m)
OCEAN_LOC = "right"        # "right", "left", "top", "bottom" in the RAW array

# True reverses alongshore order after orienting. REQUIRED for Hatteras, because
# the GIS row order and the domain numbering run in OPPOSITE directions:
#
#   1. Every resampled_domain_*.tfw has pixel Y size = -10 and zero rotation, so
#      the rasters are north-up and array ROW 0 IS THE NORTHERNMOST row. With
#      OCEAN_LOC="right" no rotation is applied, so axis 0 stays alongshore:
#      within a domain, index 0 = north, index 49 = south.
#   2. The .tfw upper-left northing increases monotonically with domain number
#      (D1 = 3,899,274 -> D90 = 3,944,002 m, 502.6 m per step), so DOMAIN NUMBER
#      INCREASES NORTHWARD. D1 = Cape Point / south, as OFFSET_ROW_ORDER and
#      SECTIONS already assume.
#
# Unflipped, the assembled island is a 90-tooth sawtooth: each 500 m block is
# internally mirrored against its neighbours. Measured on the 2009_v2 output,
# mean jump at domain seams vs. mean jump within a domain:
#
#                        as saved   flipped
#     island width        21.1x      1.97x     (134 m seam jumps -> 12.5 m)
#     dune height          3.45x     1.54x
#     mean interior elev   4.87x     1.41x
#
# WHAT THIS DOES AND DOES NOT AFFECT. BRIE resolves ONE node per 500 m Barrier3D
# domain (brie_coupler.py: dy=500, alongshore_section_count=ny) and exchanges a
# scalar x_s, so it never sees within-domain cells -- the sawtooth does not
# corrupt alongshore transport. Barrier3D's own 50-cell axis is near
# mirror-symmetric (Q1/Q3 weighting and the i>0 / i<BarrierLength-1 boundaries
# are symmetric; the router writes only into row d+1 so sweep order is
# immaterial). The one asymmetry is a bug: DiffuseDunes loops
# `range(2, BarrierLength)`, so alongshore cell 0 never exchanges sand with cell
# 1 while cell 49 does.
#
# So this flag is about GEOGRAPHIC FIDELITY AND CROSS-INPUT CONSISTENCY, which is
# where the real exposure is: NC-12 masks, community/nourishment zones and
# setback CSVs are all built in GIS row order and must share ONE frame with the
# topography, or the road sits at the mirrored alongshore position inside every
# domain. Flipping HERE -- at load, before the shear -- is what keeps the picker,
# the road masks, shear_like(), the QC figures and the saved .npy in that one
# frame. Do NOT reverse only the arrays written to disk: that is what
# PEA_dune_topo_extractor_..._alongshore_corrected.py does
# (CORRECT_SAVED_ALONGSHORE), and it leaves every figure and mask mirrored
# against the files CASCADE reads.
#
# THE SAVED PICKS ARE UNAFFECTED. i0/i1 are scalar CROSS-SHORE indices, and the
# flip commutes with straighten + water-trim: reversing alongshore negates the
# polyfit slope but leaves ref[i] -> ref[n-1-i], hence shear[i] -> shear[n-1-i]
# and an identical per-profile cross-shore frame. Verified on all 90 domains
# (z_flipped == z_unflipped[::-1] exactly, c0 and n_cross_trimmed unchanged) by
# HAT_alongshore_frame_check.py. No re-picking is needed.
ALONGSHORE_FLIP = True

# --- ROAD OVERLAY -------------------------------------------------------
# NC-12 drawn on the picker and every per-domain figure. DISPLAY AND DIAGNOSTICS
# ONLY: the road never enters find_dunes, build_interior or straighten_profiles,
# so the arrays CASCADE reads are byte-identical with SHOW_ROAD either way. What
# it changes is what YOU see while picking -- a search window that sits landward
# of the road is picking the road embankment, not a dune crest.
#
# THE MASKS ARE NOT MADE HERE. HAT_rasterize_road_to_domains.py burns the road
# geojson onto each domain's resampled_*.tif affine, so they are cell-for-cell
# aligned with the DEM .npy by construction. This script only reads them, and
# refuses anything whose shape disagrees.
#
# TWO VINTAGES ON A 2009 DEM. The road lines are 1984 and 2004; the topography is
# 2009. That mismatch IS the subject of RoadOffset_dunestart_audit.md, so both are
# drawn, in different colours, with the year in every label -- never read the
# 1984 line as 1984 topography.
SHOW_ROAD = True
ROAD_YEARS = [1984, 2004]
REQUIRE_ROAD_MASKS = True   # missing file or shape mismatch -> hard error.
                            # All 90 domains have masks for both years, so this
                            # only ever fires if the raster tree moved or a
                            # re-export changed a grid.
ROAD_RASTER_ROOT = INIT_ROOT / "4-mgmt-forcing" / "road_offset" / "raster"
ROAD_MASK_DIR_FMT = "{year}/masks"
ROAD_MASK_NAME_FMT = "domain_{domain}_road_{year}.npy"

# D1-D7 (Cape Point) have ZERO road cells in both vintages -- NC-12 does not
# reach the point. An empty mask is normal and silent; only a missing file or a
# shape mismatch is an error.
ROAD_COLORS = {1984: "#6A1B9A", 2004: "#111111"}   # purple 1984, near-black 2004
ROAD_EDGE_COLOR = "#FFFFFF"   # thin outline so the road reads on dark water AND
                              # light land, which a single fill colour cannot
ROAD_PLAN_ALPHA = 0.55        # filled road cells on a map panel
ROAD_ENVELOPE_ALPHA = 0.10    # cross-shore envelope band on a profile panel
SHOW_ROAD_CENTER_LINE = True  # dashed alongshore line through the road centre

# --- DUNE SEARCH --------------------------------------------------------
DEFAULT_WINDOW_PX = 8      # fallback window length (px landward of beach start)
CLIP_WINDOW_TO_BEACH = True  # per profile, start the search at max(i0, beach_start)
                             # so a picked window can't wander onto the wet beach
                             # where the shoreline curves within a domain

# --- INTERIOR -----------------------------------------------------------
USE_CONST_INTERIOR = False  # interior starts one cell landward of the MOST landward
                           # dune in the domain -> alongshore alignment preserved.
                           # False = each profile starts behind its own dune (v3's
                           # original behavior; breaks alongshore alignment).
FILL_MISSING_DUNE = True   # profiles with no dune found get MIN_DUNE_H_M instead of
                           # the -3.0 m sentinel (a negative dune height is not a
                           # valid Barrier3D input)
TRIM_INTERIOR_ROWS = True  # CHANGES THE .npy CASCADE READS, not just the figures.
                           # True  = Lexi's v3: drop all-water rows, so each domain's
                           #         interior array is only as tall as its own island.
                           # False = your dune_topo_extractor_from_GIS.py: every domain
                           #         is exactly (TOPO_ROWS, ALONG_COLS), padded landward
                           #         with sentinel. This is what produced the 2009_v1
                           #         arrays behind the poster figure.

# --- STRAIGHTEN ---------------------------------------------------------
# Shear each alongshore profile so the shoreline runs HORIZONTALLY, before the
# dune window is picked.
#
# Why: the clip boxes are north-up (rot = 0.00 on all 131 rasters) while
# Hatteras trends NNW, so cross-shore is due east-west and the shoreline crosses
# each 500 m domain diagonally. Domain 11's dune runs cell 3 -> 13: 8 deg of
# obliquity, 100 m of drift. Two consequences, and they are separate:
#
#   1. THE WINDOW has to be wide enough to span the diagonal. Domain 11's is
#      [3, 15] = 120 m -- also wide enough to catch a back-dune, a wooded ridge,
#      or a house. The two worst-obliquity domains (GIS 33 at 25.6 deg, GIS 67
#      at 22.8 deg) are two of the four SCATTER flags in the road setbacks.
#
#   2. THE WEDGE. USE_CONST_INTERIOR cuts horizontally at max(dune_loc) + 1,
#      throwing away everything seaward of that on every other profile.
#
# Straightening fixes (1). USE_CONST_INTERIOR = False fixes (2). Set both --
# even straightened, residual dune variability is real and the const cut still
# costs 20-40 m of it.
#
# NEITHER fixes the distance inflation: the profiles are still due east-west, so
# cross-shore AND alongshore distances stay long by 1/cos(theta) -- 1% at 8 deg,
# 4% at 16 deg, 11% at 26 deg. A "500 m" domain spans 500/cos(theta) m of
# shoreline. Only re-clipping with rotated boxes fixes that.
#
# THE PICKS BECOME FRAME-DEPENDENT. A window picked straightened is a valid
# index range on an unstraightened array; it just points at different cells.
# save_windows records STRAIGHTEN and the run pass refuses on a mismatch. Use a
# NEW WINDOW_JSON and a NEW RUN_NAME rather than overwriting a picked set.
STRAIGHTEN = True

# What to align on. "beach" = the first cell above BEACH_START_THR_M, which is
# computed without a window -- no chicken-and-egg with the dune pick.
STRAIGHTEN_REF = "beach"

# "linear"  fit a straight line to start_beach, shear by that. Over 500 m the
#           island does not curve, so the obliquity IS linear: the fit removes
#           exactly the diagonal and leaves real alongshore variability in the
#           array. Immune to a few bad profiles.
# "raw"     shear by start_beach itself. Flattens real structure too and folds
#           every noisy pixel into the geometry. Diagnostic only.
STRAIGHTEN_FIT = "linear"

# Below this many profiles with a beach, skip straightening and say so.
STRAIGHTEN_MIN_PROFILES = 10


# ==============================================================================
# ARRAY HELPERS
# ==============================================================================

def water_col_bounds(domain_array: np.ndarray, w_elev: float) -> tuple[int, int]:
    """
    First/last (inclusive) columns that are not entirely water.

    `<=` rather than `==`: no-data now carries NODATA_SENTINEL_M, which is below
    w_elev. With `==` a column of pure no-data would read as "not water" and
    survive trimming, silently changing every array's shape. Water and no-data
    are both "nothing to model here" for trimming purposes, so both trim.
    """
    keep = [c for c in range(domain_array.shape[1])
            if not np.all(domain_array[:, c] <= w_elev + 1e-9)]
    if not keep:
        raise ValueError("domain is entirely water after clamping")
    return min(keep), max(keep)


def remove_water_rows(domain_array: np.ndarray, w_elev: float) -> np.ndarray:
    """Trim leading/trailing rows that are entirely water. `<=` for the reason
    in water_col_bounds: no-data must trim like water or shapes change."""
    keep = [r for r in range(domain_array.shape[0])
            if not np.all(domain_array[r, :] <= w_elev + 1e-9)]
    if not keep:
        return domain_array[:0, :]
    return domain_array[min(keep):max(keep) + 1, :]


def orient_ocean_right(arr: np.ndarray, ocean_loc: str) -> np.ndarray:
    """
    Return an (alongshore, cross_shore) array with the ocean in the LAST column.

    NOTE: v3's "top"/"left" branches used np.flip(arr), which flips BOTH axes and
    silently reverses the alongshore order. These branches do not.
    """
    if ocean_loc == "right":
        out = arr
    elif ocean_loc == "left":
        out = arr[:, ::-1]
    elif ocean_loc == "bottom":
        out = np.rot90(arr)        # last row -> last column
    elif ocean_loc == "top":
        out = np.rot90(arr, -1)    # first row -> last column
    else:
        raise ValueError(f"OCEAN_LOC must be right/left/top/bottom, got {ocean_loc!r}")
    if ALONGSHORE_FLIP:
        out = out[::-1, :]
    return np.ascontiguousarray(out)


def straighten_profiles(z: np.ndarray, start_beach: np.ndarray):
    """
    Shear each alongshore profile so the shoreline is horizontal.

    Returns (z_sheared, start_beach_sheared, shear, obliquity_deg).

    shear[i] = cells dropped from the SEAWARD end of profile i. Anything that
    indexes the same grid afterwards -- an NC-12 mask, a dune-line mask -- has
    to be sheared with this same array via shear_like(), or it points at
    different ground than the topography does.

    The fit's slope is cells cross-shore per cell alongshore. Both axes are
    CELL_SIZE_M, so obliquity = atan(slope): the shoreline's angle to the grid's
    east-west axis, measured directly rather than inferred from dune_loc's span.
    """
    n_along, n_cross = z.shape
    shear = np.zeros(n_along, dtype=int)

    if not STRAIGHTEN:
        return z, start_beach, shear, 0.0

    sb = np.asarray(start_beach)
    ok = sb >= 0
    if int(ok.sum()) < STRAIGHTEN_MIN_PROFILES:
        print(f"       [warn] only {int(ok.sum())} profiles have a beach; "
              f"not straightening this domain")
        return z, start_beach, shear, 0.0

    x = np.arange(n_along)
    slope, intercept = np.polyfit(x[ok], sb[ok].astype(float), 1)
    obliq = float(np.degrees(np.arctan(abs(slope))))

    if STRAIGHTEN_FIT == "linear":
        ref = slope * x + intercept
    elif STRAIGHTEN_FIT == "raw":
        ref = sb.astype(float).copy()
        if (~ok).any():
            ref[~ok] = np.interp(x[~ok], x[ok], sb[ok].astype(float))
    else:
        raise ValueError(f"STRAIGHTEN_FIT must be linear/raw, "
                         f"got {STRAIGHTEN_FIT!r}")

    shear = np.round(ref - np.nanmin(ref)).astype(int)
    shear = np.clip(shear, 0, n_cross - 1)

    # Cells shifted in from beyond the original array were never surveyed
    # either, so they are no-data rather than water. They trim identically
    # (see water_col_bounds), so this changes the mask, not the topography.
    zs = np.full_like(z, NODATA_SENTINEL_M)
    for i in range(n_along):
        k = int(shear[i])
        if k:
            zs[i, :n_cross - k] = z[i, k:]
        else:
            zs[i, :] = z[i, :]

    sb_new = np.where(ok, np.maximum(sb - shear, 0), -1)
    return zs, sb_new, shear, obliq


def shear_like(arr: np.ndarray, shear: np.ndarray) -> np.ndarray:
    """
    Apply an existing shear to another array on the same grid.

    Must be the SAME shear straighten_profiles() returned for this domain. Used
    by HAT_road_setback_extract.py to put the NC-12 and dune-line masks in the
    frame the topography was saved in.
    """
    out = np.zeros_like(arr)
    n_along, n_cross = arr.shape
    for i in range(min(n_along, len(shear))):
        k = int(shear[i])
        if k:
            out[i, :n_cross - k] = arr[i, k:]
        else:
            out[i, :] = arr[i, :]
    return out


def align_mask_to_topography(raw_mask: np.ndarray, dom: dict) -> np.ndarray:
    """Put a raw GIS mask into the frame the topography was saved in.

    MOVED HERE from HAT_road_offset_from_dune_start.py (2026-08-18). It lived
    there as a private copy that re-derived this script's frame from outside it,
    which meant two definitions of the same chain that had to be kept in step by
    hand. It is now defined once, next to shear_like, and the setback script
    imports it. Any change to load_profiles' ordering has to be reflected here or
    the assert at the bottom fires.

    The chain must match ``load_profiles`` operation for operation, or the mask
    indexes different ground than ``dune_loc`` does:

        orient_ocean_right   same OCEAN_LOC and the same ALONGSHORE_FLIP
        [:, ::-1]            ocean-first, as load_profiles does to build ``raw``
        shear_like(shear)    the SAME per-profile shear, not a re-fit one
        [:, c0:c0+n_cross]   the SAME water-trim window

    ``load_profiles`` returns c0 but not c1; the trimmed width of ``z`` supplies
    the rest, which is also a check that the two arrays end up the same shape.
    """
    mask = np.squeeze(np.asarray(raw_mask))
    if mask.ndim != 2:
        raise ValueError(f"expected a 2-D mask, got shape {mask.shape}")

    binary = np.isfinite(mask) & (mask > 0)

    oriented = orient_ocean_right(binary, OCEAN_LOC)
    ocean_first = np.ascontiguousarray(oriented[:, ::-1])

    # shear_like fills with np.zeros_like, so on a bool array the cells shifted
    # in from beyond the seaward end become False -- correct for a mask, which
    # is why the bool dtype has to survive this call.
    sheared = shear_like(ocean_first, dom["shear"])

    c0 = int(dom["c0"])
    n_cross = dom["z"].shape[1]
    trimmed = sheared[:, c0:c0 + n_cross]

    if trimmed.shape != dom["z"].shape:
        raise ValueError(
            f"aligned mask {trimmed.shape} does not match topography "
            f"{dom['z'].shape}; the mask was not exported on the same grid"
        )
    return np.ascontiguousarray(trimmed, dtype=bool)


def interior_row0_line(prof_arr: np.ndarray,
                       dune_loc: np.ndarray) -> tuple[np.ndarray, int]:
    """Source cross-shore cell that becomes SAVED interior row 0, per profile.

    MOVED HERE from HAT_road_offset_from_dune_start.py alongside
    align_mask_to_topography, and for the same reason: the setback measures from
    interior row 0, this script's figures and road columns measure from interior
    row 0, and there must be exactly one definition of where that is.

    ``build_interior`` with USE_CONST_INTERIOR = False fills each column from
    ``prof_arr[i, dune_loc[i] + 1:]``, so interior row 0 is the cell one landward
    of the crest. But TRIM_INTERIOR_ROWS = True then runs ``remove_water_rows``,
    which drops leading AND trailing all-water rows -- so if interior row 0 were
    all-water across every profile, the SAVED row 0 would be a different cell
    and every setback would be off by that shift.

    It is currently zero on all 90 domains, but it is computed rather than
    assumed, because a change to the dune window or the water clamp could make
    it nonzero without any other visible symptom.

    TWO FIXES APPLIED IN THE MOVE, both of which are no-ops on the current
    settings and both of which were latent bugs in the original:
      1. the all-water test is now `<= SENTINEL_WATER_M + 1e-9`, matching
         remove_water_rows. The original used `== SENTINEL_WATER_M`, which does
         NOT catch a leading row of pure no-data (NODATA_SENTINEL_M = -99 is
         below the sentinel, so `==` kept a row that the real trim dropped).
      2. lead_trim is only applied when TRIM_INTERIOR_ROWS is True. The original
         assumed it, so with TRIM_INTERIOR_ROWS = False it would have shifted
         row 0 by a trim that never happened.
    """
    topo, start_island = build_interior(prof_arr, dune_loc)

    if TRIM_INTERIOR_ROWS:
        keep = [r for r in range(topo.shape[0])
                if not np.all(topo[r, :] <= SENTINEL_WATER_M + 1e-9)]
        lead_trim = int(min(keep)) if keep else 0
    else:
        lead_trim = 0

    if start_island is not None:
        # USE_CONST_INTERIOR: the cut is horizontal, so row 0 is the same cell on
        # every profile whether or not that profile found a dune.
        row0 = np.full(len(dune_loc), start_island + lead_trim, dtype=int)
    else:
        row0 = np.where(dune_loc >= 0, dune_loc + 1 + lead_trim, -1)
    return row0, lead_trim


def natural_key(name: str) -> tuple:
    m = re.search(r"domain_(\d+)", name)
    return (int(m.group(1)) if m else 10**9, name)


def domain_number(stem: str) -> int | None:
    m = re.search(r"domain_(\d+)", stem)
    return int(m.group(1)) if m else None


def section_for(stem: str) -> str:
    """Island section label for a domain, per the D1=Cape Point convention."""
    n = domain_number(stem)
    if n is None:
        return ""
    for (lo, hi), label in SECTIONS:
        if lo <= n <= hi:
            return label
    return ""


def fig_stem(stem: str) -> str:
    """Zero-padded figure name so 90+ domains sort correctly in Explorer."""
    n = domain_number(stem)
    return f"domain_{n:03d}" if n is not None else stem


def fig_title(stem: str) -> str:
    sec = section_for(stem)
    return f"{stem}  ({sec})" if sec else stem


# ==============================================================================
# ROAD MASKS
#
# Read-only consumers of HAT_rasterize_road_to_domains.py. Nothing here writes a
# mask, and nothing here influences the dune search or the saved arrays.
# ==============================================================================

def road_mask_path(year: int, domain_id: int) -> Path:
    """Where HAT_rasterize_road_to_domains.py put one domain's mask."""
    return (Path(ROAD_RASTER_ROOT) / ROAD_MASK_DIR_FMT.format(year=year)
            / ROAD_MASK_NAME_FMT.format(domain=domain_id, year=year))


def load_road_masks(dem_path: Path, dom: dict) -> tuple[dict, dict]:
    """Load every ROAD_YEARS mask for one domain, in both frames.

    Returns (road_raw, road_aligned), each {year: bool array}:

        road_raw      (n_along, n_raw)   ocean-first, UNSHEARED and UNTRIMMED --
                                        the frame panel 1 of the comparison
                                        figure draws, so the road's real diagonal
                                        across the domain stays visible
        road_aligned  (n_along, n_cross) the frame the topography and dune_loc
                                        live in, via align_mask_to_topography

    Empty dicts if SHOW_ROAD is False. A domain with no road cells still gets an
    all-False entry rather than being omitted -- D1-D7 are that case every run,
    and callers should not have to distinguish "no road here" from "not loaded".
    """
    if not SHOW_ROAD:
        return {}, {}

    domain_id = domain_number(Path(dem_path).stem)
    if domain_id is None:
        raise ValueError(f"cannot determine domain ID for road mask: {dem_path.name}")

    raw_out, aligned_out = {}, {}
    for year in ROAD_YEARS:
        path = road_mask_path(year, domain_id)
        if not path.is_file():
            msg = f"no {year} road mask for domain {domain_id}: {path}"
            if REQUIRE_ROAD_MASKS:
                raise FileNotFoundError(msg)
            print(f"       [road warn] {msg}")
            continue

        mask = np.squeeze(np.load(path))
        if mask.ndim != 2:
            raise ValueError(f"{path.name}: expected a 2-D mask, got {mask.shape}")
        # The mask is checked against the RAW DEM shape, before any orienting, so
        # a grid mismatch is reported against the thing the rasterizer actually
        # snapped to. No transpose or resize: that would hide a misregistration.
        if tuple(mask.shape) != tuple(dom["raw_shape_unoriented"]):
            raise ValueError(
                f"{path.name}: road shape {mask.shape} does not match DEM shape "
                f"{tuple(dom['raw_shape_unoriented'])}. Re-run "
                f"HAT_rasterize_road_to_domains.py for {year}."
            )

        binary = np.isfinite(mask) & (mask > 0)
        oriented = orient_ocean_right(binary, OCEAN_LOC)
        raw_out[year] = np.ascontiguousarray(oriented[:, ::-1], dtype=bool)
        aligned_out[year] = align_mask_to_topography(mask, dom)

    return raw_out, aligned_out


def road_profile_positions(mask: np.ndarray | None) -> tuple[np.ndarray, np.ndarray,
                                                             np.ndarray]:
    """Seaward edge, landward edge and centre cell of the road, per profile.

    NaN on profiles the road does not cross, which is a real state -- the road
    leaves the domain, or the domain has no road at all.
    """
    if mask is None:
        return (np.array([], dtype=float),) * 3
    n_along = mask.shape[0]
    seaward = np.full(n_along, np.nan)
    landward = np.full(n_along, np.nan)
    center = np.full(n_along, np.nan)
    for i in range(n_along):
        cells = np.flatnonzero(mask[i])
        if cells.size:
            seaward[i] = float(cells.min())
            landward[i] = float(cells.max())
            center[i] = float(np.mean(cells))
    return seaward, landward, center


def processed_road_grid(mask: np.ndarray | None, row0_line: np.ndarray,
                        n_rows: int, n_cols: int) -> np.ndarray:
    """Re-index road cells into the SAVED interior grid.

    The interior is cut per profile at interior row 0, so a road cell's row in
    the saved array is (source cell - row0[i]). This is what shows whether NC-12
    is inside the array CASCADE actually reads, and at which interior row --
    which is the same quantity roadway_manager.bulldoze indexes with
    int(road_setback / dy).
    """
    out = np.zeros((n_rows, n_cols), dtype=bool)
    if mask is None:
        return out
    n_fill = min(n_cols, mask.shape[0], len(row0_line))
    for i in range(n_fill):
        start = int(row0_line[i])
        if start < 0:
            continue
        dest = np.flatnonzero(mask[i]) - start
        keep = (dest >= 0) & (dest < n_rows)
        out[dest[keep], i] = True
    return out


def road_offset_stats(mask: np.ndarray | None,
                      row0_line: np.ndarray) -> dict:
    """Per-domain road geometry and setback from SAVED interior row 0.

    `setback_median_m` IS `setback_2009_m` from RoadOffset_<year>_domains.csv --
    same reference row, same frame, same shear, same edge, same statistic. It is
    computed here independently so the two can be diffed; if they disagree, one
    of the two frames has drifted.

    MEASURED FROM THE SEAWARD EDGE OF THE ROAD BLOCK, NOT ITS CENTRE, and
    reported as a MEDIAN over profiles. Both of those match
    HAT_road_offset_from_dune_start.py, and neither is arbitrary:

      * the seaward edge is what `roadway_manager.bulldoze` indexes --
        `[road_start : road_start + road_width]` starts at the seaward edge, so
        that is the cell `int(road_setback / dy)` has to land on. Measuring the
        centre instead reads high by half the mask width, which on a ~24 m mask
        (ROAD_BUFFER_M = 6 plus all_touched on 10 m cells) is a systematic ~10 m.
      * the median, because NC-12 leaves some domains diagonally: the handful of
        profiles where the road clips a corner drag a mean by up to ~150 m while
        the median stays on the road proper.

    The centre and the width are still reported, as geometry -- they are what
    tells you the mask is a fat buffer around an 8 m road rather than the road.

    Sign: POSITIVE = road LANDWARD of interior row 0 (the normal case, and the
    only one Barrier3D can represent). Negative means the road sits seaward of
    the dune line, which is what the audit's NEGATIVE floor is about.
    """
    empty = {
        "road_profiles": 0, "road_cells": 0,
        "road_span_cells": np.nan, "road_center_cell": np.nan,
        "road_width_cells": np.nan,
        "setback_median_m": np.nan, "setback_mean_m": np.nan,
        "setback_min_m": np.nan, "setback_max_m": np.nan,
        "center_median_m": np.nan, "n_seaward": 0,
    }
    if mask is None or not np.any(mask):
        return empty

    seaward, landward, center = road_profile_positions(mask)

    # Cap at ALONG_COLS, as the setback script does: profiles beyond the 50
    # CASCADE keeps are not part of the measurement.
    n = min(len(center), len(row0_line), ALONG_COLS)
    row0 = np.asarray(row0_line[:n], dtype=float)
    valid = np.isfinite(center[:n]) & (row0 >= 0)

    out = dict(empty)
    out["road_profiles"] = int(valid.sum())
    out["road_cells"] = int(np.count_nonzero(mask[:n]))
    if np.isfinite(center[:n]).any():
        out["road_span_cells"] = float(np.nanmax(landward[:n])
                                      - np.nanmin(seaward[:n]) + 1)
        out["road_center_cell"] = float(np.nanmedian(center[:n]))

    if valid.any():
        sb = (seaward[:n][valid] - row0[valid]) * CELL_SIZE_M
        ctr = (center[:n][valid] - row0[valid]) * CELL_SIZE_M
        width = (landward[:n][valid] - seaward[:n][valid] + 1)
        out.update({
            "setback_median_m": float(np.median(sb)),
            "setback_mean_m": float(np.mean(sb)),
            "setback_min_m": float(np.min(sb)),
            "setback_max_m": float(np.max(sb)),
            "center_median_m": float(np.median(ctr)),
            "road_width_cells": float(np.median(width)),
            # counted on the SEAWARD edge, because that is the value that gets
            # floored before it reaches the model
            "n_seaward": int(np.count_nonzero(sb < 0)),
        })
    return out


def add_road_plan_overlay(ax, mask: np.ndarray | None, year: int,
                          *, zorder: float = 5.0, label: bool = True) -> None:
    """Draw exact road cells on an (alongshore x, cross-shore y) map panel.

    The white contour is not decoration: the fill sits on a terrain colormap that
    runs from near-black water to pale land, and no single fill colour is legible
    against both. The outline is.
    """
    if mask is None or not np.any(mask):
        return
    color = ROAD_COLORS.get(year, "#111111")
    n_along, n_cross = mask.shape
    display = np.ma.masked_where(~mask.T, np.ones((n_cross, n_along)))
    ax.imshow(display, aspect="auto", origin="lower",
              extent=[-0.5, n_along - 0.5, -0.5, n_cross - 0.5],
              interpolation="nearest", cmap=ListedColormap([color]),
              vmin=0.0, vmax=1.0, alpha=ROAD_PLAN_ALPHA, zorder=zorder)
    if np.any(~mask):
        ax.contour(np.arange(n_along), np.arange(n_cross), mask.T.astype(float),
                   levels=[0.5], colors=ROAD_EDGE_COLOR, linewidths=0.7,
                   zorder=zorder + 0.2)
    if SHOW_ROAD_CENTER_LINE:
        _, _, center = road_profile_positions(mask)
        if np.isfinite(center).any():
            ax.plot(np.arange(n_along), center, color=color, lw=1.2, ls="--",
                    zorder=zorder + 0.4)
    if label:
        ax.plot([], [], color=color, lw=5, alpha=ROAD_PLAN_ALPHA,
                label=f"NC-12 {year}")


def add_road_envelope(ax, mask: np.ndarray | None, year: int,
                      *, label: bool = True) -> None:
    """Road's cross-shore envelope on a profile-stack panel (elevation x, cell y).

    A profile stack has no alongshore axis, so the road cannot be drawn cell by
    cell -- what it can show is the band of cross-shore cells the road occupies
    anywhere in the domain, which is what you compare the search window against.
    """
    if mask is None or not np.any(mask):
        return
    color = ROAD_COLORS.get(year, "#111111")
    seaward, landward, center = road_profile_positions(mask)
    if not np.isfinite(center).any():
        return
    ax.axhspan(float(np.nanmin(seaward)) - 0.5, float(np.nanmax(landward)) + 0.5,
               color=color, alpha=ROAD_ENVELOPE_ALPHA, zorder=0,
               label=f"NC-12 {year} envelope" if label else None)
    ax.axhline(float(np.nanmedian(center)), color=color, lw=1.1, ls="--",
               zorder=1)


# ==============================================================================
# LOAD / PREP
# ==============================================================================

def load_profiles(in_path: Path) -> dict:
    """
    Load one domain. Returns a dict:
        raw         : (n_along, n_cross_raw) RAW GIS elevation in m NAVD88,
                      oriented OCEAN-FIRST, untrimmed. For the comparison figure.
        z           : (n_along, n_cross) MHW-relative, clamped, water-trimmed,
                      OCEAN-FIRST (index 0 = ocean). This is what gets processed.
        start_beach : (n_along,) first index in z where z > BEACH_START_THR_M,
                      -1 if none
        c0          : cross-shore raw_offset of z within raw. WITHOUT
                      straightening, z[:, k] == raw column k + c0. WITH
                      straightening the mapping is per profile:
                          z[i, k] == raw[i, k + c0 + shear[i]]
        shear       : (n_along,) cells dropped from the seaward end of each
                      profile to make the shoreline horizontal. Zeros if
                      STRAIGHTEN is False. Any mask that has to index the same
                      grid must be put through shear_like() with THIS array.
        obliquity_deg : the shoreline's angle to the grid's east-west axis,
                      from the slope of the start_beach fit. 0.0 if not
                      straightened.
        road_raw    : {year: bool (n_along, n_cross_raw)} NC-12, ocean-first,
                      unsheared and untrimmed. Empty if SHOW_ROAD is False.
        road_masks  : {year: bool (n_along, n_cross)} NC-12 in the SAME frame as
                      z, via align_mask_to_topography. Empty if SHOW_ROAD is
                      False. Display and diagnostics only -- nothing downstream
                      of here lets the road affect the dune search or the arrays.
    """
    arr = np.load(in_path).astype(float, copy=False)
    if arr.ndim != 2:
        raise ValueError(f"expected 2D array, got {arr.ndim}D")

    raw_shape_unoriented = arr.shape   # the shape the road rasterizer snapped to
    arr = orient_ocean_right(arr, OCEAN_LOC)

    # sanity check that OCEAN_LOC is actually right (the original script's AUTO_ORIENT,
    # demoted to a warning so orientation stays an explicit, documented choice)
    edge = max(1, min(5, arr.shape[1] // 20))
    left_q = np.nanpercentile(arr[:, :edge], 25)
    right_q = np.nanpercentile(arr[:, -edge:], 25)
    if right_q > left_q:
        print(f"[warn] {in_path.name}: after orienting, the LEFT edge is lower "
              f"(p25 left={left_q:.2f}, right={right_q:.2f}). Check OCEAN_LOC.")

    raw = np.ascontiguousarray(arr[:, ::-1])  # NAVD88, ocean first, untrimmed

    # No-data is identified on the RAW array, before the clamp folds it into the
    # water sentinel. It then rides the same shear/trim/slice path as z, marked
    # by a value far below the clamp, and is separated out again at save time.
    nodata = raw <= RAW_NODATA_MAX_NAVD

    z = raw - MHW_M
    z[z < WATER_CLAMP_M] = SENTINEL_WATER_M
    z[nodata] = NODATA_SENTINEL_M

    n_along = z.shape[0]
    if n_along < ALONG_COLS:
        print(f"[warn] {in_path.name}: alongshore={n_along} < {ALONG_COLS}; "
              f"trailing output cols remain sentinel.")
    elif n_along > ALONG_COLS:
        print(f"[warn] {in_path.name}: alongshore={n_along} > {ALONG_COLS}; "
              f"only first {ALONG_COLS} profiles used.")

    # ORDER MATTERS: start_beach -> straighten -> water trim.
    # start_beach is found on the untrimmed array because the shear is what
    # defines the frame; c0 must then be measured on the array the window is
    # actually picked in, or the two disagree.
    above = z > BEACH_START_THR_M
    start_beach = np.where(above.any(axis=1), above.argmax(axis=1), -1)

    z, start_beach, shear, obliq = straighten_profiles(z, start_beach)

    c0, c1 = water_col_bounds(z, SENTINEL_WATER_M)
    z = z[:, c0:c1 + 1]
    start_beach = np.where(start_beach >= 0,
                           np.maximum(start_beach - c0, 0), -1)

    dom = {"raw": raw, "z": z, "start_beach": start_beach, "c0": int(c0),
           "shear": shear, "obliquity_deg": round(float(obliq), 2),
           "name": in_path.name,
           "raw_shape_unoriented": raw_shape_unoriented}

    # LAST, because align_mask_to_topography needs the finished z, c0 and shear.
    # Loading the road cannot change any of them -- if it ever appears to, the
    # shape assert inside align_mask_to_topography is what will say so.
    dom["road_raw"], dom["road_masks"] = load_road_masks(in_path, dom)
    return dom


def masked_profiles(prof_arr: np.ndarray) -> np.ndarray:
    """NaN out sentinel water cells for plotting/statistics."""
    return np.where(prof_arr <= SENTINEL_WATER_M + 1e-9, np.nan, prof_arr)


def default_window(prof_arr: np.ndarray, start_beach: np.ndarray) -> tuple[int, int]:
    valid = start_beach[start_beach >= 0]
    base = int(np.median(valid)) if valid.size else 0
    return base, min(base + DEFAULT_WINDOW_PX, prof_arr.shape[1])


# ==============================================================================
# INTERACTIVE PICKER
# ==============================================================================

def _span_selector(ax, on_select):
    """SpanSelector with a props/rectprops fallback for matplotlib < 3.5."""
    try:
        return SpanSelector(ax, on_select, "vertical", useblit=True,
                            props=dict(alpha=0.2, facecolor="#FF8C00"))
    except TypeError:
        return SpanSelector(ax, on_select, "vertical", useblit=True,
                            rectprops=dict(alpha=0.2, facecolor="#FF8C00"))


def pick_window(stem: str, prof_arr: np.ndarray, start_beach: np.ndarray,
                init: tuple[int, int],
                road_masks: dict | None = None) -> tuple[str, int, int]:
    """
    Show the domain ocean-at-bottom and let the user drag a dune search window.

    Display only: the cross-shore axis runs vertically with cell 0 (ocean) at the
    bottom, landward upward. i0/i1 are still cross-shore indices from the ocean.

    NC-12 is drawn on both panels (v4). It is there to tell you when a window has
    wandered onto the road embankment: the road is a hard, flat-topped ridge a
    dune-crest argmax will happily lock onto, and the v3 windows were picked
    without being able to see it. The road never constrains the window in code --
    only your eye.
    """
    n_along, n_cross = prof_arr.shape
    zm = masked_profiles(prof_arr)
    state = {"i0": init[0], "i1": init[1], "action": None}

    fig, (ax_map, ax_prof) = plt.subplots(
        1, 2, figsize=(13, 9), sharey=True,
        gridspec_kw={"width_ratios": [1.3, 1.0], "wspace": 0.06},
    )
    try:
        fig.canvas.manager.set_window_title(f"dune search window - {stem}")
    except Exception:
        pass

    finite = zm[np.isfinite(zm)]
    vmax = float(np.percentile(finite, 99)) if finite.size else 3.0

    # --- map panel: alongshore on x, cross-shore on y, ocean at the bottom ---
    im = ax_map.imshow(
        np.ma.masked_invalid(zm.T), aspect="auto", origin="lower",
        extent=[-0.5, n_along - 0.5, -0.5, n_cross - 0.5],
        cmap="terrain", vmin=-1.0, vmax=max(vmax, 2.0),
    )
    sb = np.where(start_beach >= 0, start_beach, np.nan)
    ax_map.plot(np.arange(n_along), sb, color="k", lw=1.2, label="beach start")
    for _yr, _m in (road_masks or {}).items():
        add_road_plan_overlay(ax_map, _m, _yr)
    ax_map.set_xlabel("alongshore cell")
    ax_map.set_ylabel("cross-shore cell  (0 = ocean, landward up)")
    ax_map.legend(loc="upper right", fontsize=9, framealpha=0.9)

    # --- profile panel: elevation on x, cross-shore on y ---
    y = np.arange(n_cross)
    ax_prof.plot(zm.T, y, color="0.75", lw=0.6)
    med = np.nanmedian(zm, axis=0)
    ax_prof.plot(med, y, color="k", lw=2.0, label="median profile")
    ax_prof.axvline(BEACH_START_THR_M, color="#1565C0", ls="--", lw=1.0,
                    label=f"beach thr {BEACH_START_THR_M} m")
    ax_prof.axvline(BERM_ELEV_NAVD_M - MHW_M, color="#B71C1C", ls=":", lw=1.2,
                    label=f"berm {BERM_ELEV_NAVD_M} m NAVD88")
    for _yr, _m in (road_masks or {}).items():
        add_road_envelope(ax_prof, _m, _yr)
    ax_prof.set_xlabel("elev (m MHW)")
    ax_prof.set_xlim(-1.2, max(vmax, 2.0) + 0.5)
    ax_prof.set_ylim(-0.5, n_cross - 0.5)
    ax_prof.legend(loc="upper right", fontsize=9, framealpha=0.9)
    plt.setp(ax_prof.get_yticklabels(), visible=False)

    try:
        fig.colorbar(im, ax=[ax_map, ax_prof], location="right",
                     fraction=0.035, pad=0.02, label="elev (m MHW)")
    except (TypeError, ValueError):
        fig.colorbar(im, ax=ax_prof, label="elev (m MHW)")

    spans = [None, None]

    def redraw():
        for k, ax in enumerate((ax_map, ax_prof)):
            if spans[k] is not None:
                spans[k].remove()
            spans[k] = ax.axhspan(state["i0"], state["i1"], color="#FF8C00",
                                  alpha=0.25, zorder=0)
        fig.suptitle(
            f"{stem}   search window = [{state['i0']}, {state['i1']}]  "
            f"({state['i1'] - state['i0']} cells)\n"
            "drag vertically on either panel | enter = accept | r = reset | "
            "s = skip | q = quit",
            fontsize=11,
        )
        fig.canvas.draw_idle()

    redraw()

    def on_select(ymin, ymax):
        i0 = int(np.clip(round(ymin), 0, n_cross - 2))
        i1 = int(np.clip(round(ymax), i0 + 1, n_cross))
        state["i0"], state["i1"] = i0, i1
        redraw()

    def on_key(event):
        if event.key in ("enter", "return"):
            state["action"] = "accept"
            plt.close(fig)
        elif event.key == "r":
            state["i0"], state["i1"] = init
            redraw()
        elif event.key == "s":
            state["action"] = "skip"
            plt.close(fig)
        elif event.key in ("q", "escape"):
            state["action"] = "quit"
            plt.close(fig)

    selectors = [_span_selector(ax_map, on_select), _span_selector(ax_prof, on_select)]

    fig.canvas.mpl_connect("key_press_event", on_key)
    plt.show()  # blocks
    del selectors

    if state["action"] is None:  # window closed with the X
        state["action"] = "accept"
    return state["action"], state["i0"], state["i1"]


# ==============================================================================
# EXTRACTION
# ==============================================================================

def find_dunes(prof_arr, start_beach, i0, i1):
    """Dune elevation and location per profile within the picked window."""
    n_along, n_cross = prof_arr.shape
    dune_elev = np.full(n_along, np.nan)
    dune_loc = np.full(n_along, -1, dtype=int)

    for i in range(n_along):
        a = i0
        if CLIP_WINDOW_TO_BEACH and start_beach[i] >= 0:
            a = max(i0, int(start_beach[i]))
        b = min(i1, n_cross)
        if b <= a:
            continue
        w = prof_arr[i, a:b]
        valid = w > SENTINEL_WATER_M + 1e-9
        if not valid.any():
            continue
        # argmax INSIDE the window. v3 used np.where(prof == dune_elev)[0][0] over
        # the whole profile, which can snap the dune onto an earlier cell of equal
        # elevation (common on a quantized DEM).
        k = int(np.argmax(np.where(valid, w, -np.inf)))
        dune_elev[i] = float(w[k])
        dune_loc[i] = a + k

    return dune_elev, dune_loc


def build_interior(prof_arr, dune_loc):
    """Return interior topography, (cross_shore_rows, alongshore_cols), ocean-first."""
    n_along, n_cross = prof_arr.shape
    topo = np.full((TOPO_ROWS, ALONG_COLS), SENTINEL_WATER_M, dtype=float)
    n_fill = min(ALONG_COLS, n_along)

    if USE_CONST_INTERIOR:
        valid = dune_loc[dune_loc >= 0]
        if valid.size == 0:
            return topo, None
        start_island = int(valid.max()) + 1
        block = prof_arr[:n_fill, start_island:].T  # (cross, along)
        rows = min(TOPO_ROWS, block.shape[0])
        topo[:rows, :n_fill] = block[:rows, :]
        if block.shape[0] > TOPO_ROWS:
            print(f"       [warn] interior truncated: {block.shape[0]} -> {TOPO_ROWS} rows")
        return topo, start_island

    for i in range(n_fill):
        if dune_loc[i] < 0:
            continue
        use = prof_arr[i, dune_loc[i] + 1:-1]
        if use.size:
            rows = min(TOPO_ROWS, use.size)
            topo[:rows, i] = use[:rows]
    return topo, None


def extract_domain(stem, prof_arr, start_beach, i0, i1, topo_dir, dune_dir,
                   shear=None, obliquity_deg=0.0, road_masks=None):
    n_along, n_cross = prof_arr.shape
    dune_elev, dune_loc = find_dunes(prof_arr, start_beach, i0, i1)

    n_found = int(np.isfinite(dune_elev).sum())
    if n_found == 0:
        print(f"[skip] {stem}: no dune found in window [{i0}, {i1}]")
        return None

    dune_h = dune_elev - (BERM_ELEV_NAVD_M - MHW_M)  # both MHW-relative
    dune_h[np.isfinite(dune_h) & (dune_h < 0.0)] = MIN_DUNE_H_M

    dune_m = np.full(ALONG_COLS, SENTINEL_WATER_M, dtype=float)
    n_fill = min(ALONG_COLS, n_along)
    dune_m[:n_fill] = dune_h[:n_fill]
    missing = ~np.isfinite(dune_m[:n_fill])
    if missing.any():
        idx = np.where(missing)[0].tolist()
        if FILL_MISSING_DUNE:
            dune_m[:n_fill][missing] = MIN_DUNE_H_M
            print(f"       [warn] {stem}: no dune at profiles {idx} -> "
                  f"filled with {MIN_DUNE_H_M} m")
        else:
            dune_m[:n_fill][missing] = SENTINEL_WATER_M
            print(f"       [warn] {stem}: no dune at profiles {idx} -> sentinel")

    topo_m, start_island = build_interior(prof_arr, dune_loc)
    if TRIM_INTERIOR_ROWS:
        topo_m = remove_water_rows(topo_m, SENTINEL_WATER_M)

    # Where SAVED interior row 0 sits on each profile, and the road measured
    # against it. Diagnostics only -- computed from dune_loc and the same
    # build_interior call the arrays came from, and written to nothing but the
    # settings sheet and the figures.
    row0_line, lead_trim = interior_row0_line(prof_arr, dune_loc)
    road_stats = {yr: road_offset_stats(m, row0_line)
                  for yr, m in (road_masks or {}).items()}

    # Split no-data back out. The topography written is byte-identical to what
    # this script produced before the mask existed: every no-data cell goes back
    # to SENTINEL_WATER_M, because Barrier3D has no representation for
    # "unknown". The mask is what carries the distinction forward.
    topo_nodata = topo_m <= NODATA_SENTINEL_M + 1e-9
    dune_nodata = dune_m <= NODATA_SENTINEL_M + 1e-9
    topo_m = np.where(topo_nodata, SENTINEL_WATER_M, topo_m)
    dune_m = np.where(dune_nodata, SENTINEL_WATER_M, dune_m)

    topo_dm = topo_m * 0.1
    dune_dm = dune_m * 0.1

    topo_out = Path(topo_dir) / f"{stem}_topography_{TAG}.npy"
    dune_out = Path(dune_dir) / f"{stem}_dune_{TAG}.npy"
    np.save(topo_out, topo_dm)
    np.save(dune_out, dune_dm)
    if WRITE_NODATA_MASK:
        np.save(Path(topo_dir) / f"{stem}_nodata_{TAG}.npy", topo_nodata)
        if topo_nodata.any():
            print(f"       [nodata] {stem}: {int(topo_nodata.sum())} of "
                  f"{topo_nodata.size} interior cells never surveyed "
                  f"({topo_nodata.mean():.1%})")
    print(f"[ok] {stem}: window [{i0}, {i1}], {n_found}/{n_along} dunes, "
          f"interior {topo_dm.shape}, mean dune h = "
          f"{np.nanmean(dune_h):.2f} m")

    return {
        "stem": stem, "i0": i0, "i1": i1, "n_found": n_found, "n_along": n_along,
        "mean_dune_h_m": float(np.nanmean(dune_h)),
        "min_dune_h_m": float(np.nanmin(dune_h)),
        "max_dune_h_m": float(np.nanmax(dune_h)),
        "n_filled": int(missing.sum()),
        "interior_rows": int(topo_dm.shape[0]),
        "interior_cols": int(topo_dm.shape[1]),
        "mean_interior_elev_m": float(np.mean(topo_m[topo_m > SENTINEL_WATER_M + 1e-9]))
        if np.any(topo_m > SENTINEL_WATER_M + 1e-9) else np.nan,
        "start_island": start_island,
        # the saved arrays are in the straightened frame; nothing about a .npy
        # says so, and a window or a mask from the other frame is silently wrong
        "straightened": bool(STRAIGHTEN),
        "obliquity_deg": obliquity_deg,
        "shear_max_cells": int(np.max(shear)) if shear is not None else 0,
        "shear": shear,
        "row0_line": row0_line, "lead_trim_rows": lead_trim,
        "road_stats": road_stats,
        "dune_elev": dune_elev, "dune_loc": dune_loc, "dune_h": dune_h,
        "topo_dm": topo_dm, "dune_dm": dune_dm,
        "topo_file": topo_out.name, "dune_file": dune_out.name,
    }


# ==============================================================================
# QC FIGURE
# ==============================================================================

def qc_figure(stem, prof_arr, start_beach, res, fig_dir: Path,
              road_masks: dict | None = None):
    """QC figure in the same ocean-at-bottom orientation as the picker."""
    n_along, n_cross = prof_arr.shape
    zm = masked_profiles(prof_arr)
    i0, i1 = res["i0"], res["i1"]
    dl = np.where(res["dune_loc"] >= 0, res["dune_loc"], np.nan)

    fig = plt.figure(figsize=(12.5, 9.5))
    gs = fig.add_gridspec(2, 2, width_ratios=[1.3, 1.0], height_ratios=[1.0, 0.35],
                          wspace=0.06, hspace=0.12)
    ax_map = fig.add_subplot(gs[0, 0])
    ax_prof = fig.add_subplot(gs[0, 1], sharey=ax_map)
    ax_h = fig.add_subplot(gs[1, 0], sharex=ax_map)

    finite = zm[np.isfinite(zm)]
    vmax = float(np.percentile(finite, 99)) if finite.size else 3.0

    # --- map: alongshore x, cross-shore y (ocean at bottom) ---
    im = ax_map.imshow(np.ma.masked_invalid(zm.T), aspect="auto", origin="lower",
                       extent=[-0.5, n_along - 0.5, -0.5, n_cross - 0.5],
                       cmap="terrain", vmin=-1.0, vmax=max(vmax, 2.0))
    ax_map.axhspan(i0, i1, color="#FF8C00", alpha=0.20, zorder=0)
    ax_map.plot(np.arange(n_along), np.where(start_beach >= 0, start_beach, np.nan),
                color="k", lw=1.0, label="beach start")
    ax_map.plot(np.arange(n_along), dl, color="#B71C1C", lw=0, marker=".", ms=4,
                label="dune crest")
    if res["start_island"] is not None:
        ax_map.axhline(res["start_island"], color="w", lw=1.2, ls="--",
                       label="interior start")
    for _yr, _m in (road_masks or {}).items():
        add_road_plan_overlay(ax_map, _m, _yr)
    ax_map.set_ylabel("cross-shore cell  (0 = ocean, landward up)")
    ax_map.legend(loc="upper right", fontsize=8, framealpha=0.9)
    plt.setp(ax_map.get_xticklabels(), visible=False)

    # --- profiles: elevation x, cross-shore y ---
    y = np.arange(n_cross)
    ax_prof.plot(zm.T, y, color="0.75", lw=0.5)
    ax_prof.plot(np.nanmedian(zm, axis=0), y, color="k", lw=2.0)
    ax_prof.axhspan(i0, i1, color="#FF8C00", alpha=0.20, zorder=0)
    ax_prof.axvline(BERM_ELEV_NAVD_M - MHW_M, color="#B71C1C", ls=":", lw=1.2)
    ax_prof.plot(res["dune_elev"], dl, lw=0, marker=".", ms=4, color="#B71C1C")
    for _yr, _m in (road_masks or {}).items():
        add_road_envelope(ax_prof, _m, _yr, label=False)
    ax_prof.set_xlabel("elev (m MHW)")
    ax_prof.set_xlim(-1.2, max(vmax, 2.0) + 0.5)
    ax_prof.set_ylim(-0.5, n_cross - 0.5)
    plt.setp(ax_prof.get_yticklabels(), visible=False)

    # --- dune height alongshore, aligned under the map ---
    ax_h.plot(np.arange(n_along), res["dune_h"], color="#FF8C00", lw=1.5,
              marker="o", ms=3)
    ax_h.axhline(MIN_DUNE_H_M, color="0.5", ls="--", lw=1.0)
    ax_h.set_xlabel("alongshore cell")
    ax_h.set_ylabel("dune height\nabove berm (m)")
    ax_h.set_xlim(-0.5, n_along - 0.5)

    try:
        fig.colorbar(im, ax=[ax_map, ax_prof], location="right",
                     fraction=0.035, pad=0.02, label="elev (m MHW)")
    except (TypeError, ValueError):
        fig.colorbar(im, ax=ax_prof, label="elev (m MHW)")

    road_note = ""
    for _yr in ROAD_YEARS:
        st = (res.get("road_stats") or {}).get(_yr)
        if st and np.isfinite(st["setback_median_m"]):
            road_note += (f"   |  NC-12 {_yr}: {st['setback_median_m']:+.0f} m "
                          f"from interior row 0")
        elif st is not None:
            road_note += f"   |  NC-12 {_yr}: not in domain"
    fig.suptitle(f"{fig_title(stem)}  |  search window [{i0}, {i1}]  |  "
                 f"{res['n_found']}/{n_along} dunes found{road_note}",
                 fontsize=12)

    fig_dir.mkdir(parents=True, exist_ok=True)
    fig.savefig(fig_dir / f"{fig_stem(stem)}_qc.png", dpi=130,
                bbox_inches="tight")
    plt.close(fig)


# ==============================================================================
# GIS vs PROCESSED COMPARISON FIGURE
# ==============================================================================

def comparison_figure(stem, dom, res, fig_dir: Path):
    """
    The whole chain, left to right: what came in, what you picked on, what
    CASCADE gets.

    RAW           the DEM as GIS exported it (m NAVD88, untrimmed). The clip
                  boxes are north-up while the island trends NNW, so the
                  shoreline crosses each domain diagonally -- and so does every
                  overlay, because they are mapped back through the per-profile
                  shear: raw[i, k + c0 + shear[i]]. The diagonal you see here is
                  what the shear removes.

    STRAIGHTENED  the array the picker showed you and the window was drawn on
                  (m MHW, clamped, sheared, water-trimmed). Same overlays, now
                  horizontal. If the beach-start line still slopes here, the
                  linear fit did not capture that domain's shoreline.

    PROCESSED     the .npy files CASCADE reads, converted back to m for display.

    With STRAIGHTEN = False the first two panels are the same picture, minus the
    trim and the datum shift. That is the point of showing both.
    """
    raw = dom["raw"]
    zs = dom["z"]
    c0 = dom["c0"]
    n_along, n_raw = raw.shape
    n_cross = zs.shape[1]
    i0, i1 = res["i0"], res["i1"]
    xs = np.arange(n_along)

    shear = np.asarray(dom.get("shear", np.zeros(n_along, dtype=int)))
    off = c0 + shear.astype(float)          # processed cell k -> raw column k+off
    obliq = dom.get("obliquity_deg", 0.0)

    dl = np.where(res["dune_loc"] >= 0, res["dune_loc"].astype(float), np.nan)
    sb = np.where(dom["start_beach"] >= 0, dom["start_beach"].astype(float),
                  np.nan)

    topo_m = res["topo_dm"] * 10.0          # back to m MHW
    dune_m = res["dune_dm"] * 10.0          # m above berm
    topo_disp = np.where(topo_m <= SENTINEL_WATER_M + 1e-9, np.nan, topo_m)
    n_int_rows, n_int_cols = topo_disp.shape

    fig = plt.figure(figsize=(19, 9))
    gs = fig.add_gridspec(2, 3, width_ratios=[1.0, 1.0, 1.0],
                          height_ratios=[1.0, 0.14], wspace=0.22, hspace=0.06)
    ax_raw = fig.add_subplot(gs[:, 0])
    ax_str = fig.add_subplot(gs[:, 1])
    ax_int = fig.add_subplot(gs[0, 2])
    ax_dun = fig.add_subplot(gs[1, 2])

    # Shared elevation span so the three panels are comparable by eye -- just
    # labelled NAVD88 vs MHW. Scaling each to its own percentile lets the -10 m
    # shoreface swamp the ramp and the island reads as one flat colour.
    land = raw[raw > MHW_M]
    hi_navd = float(np.percentile(land, 99)) if land.size else MHW_M + 4.0
    vmin_r = MHW_M - 1.0
    vmax_r = max(hi_navd, MHW_M + 2.0)

    # ---------------- 1: raw GIS ----------------
    im_r = ax_raw.imshow(raw.T, aspect="auto", origin="lower",
                         extent=[-0.5, n_along - 0.5, -0.5, n_raw - 0.5],
                         cmap="terrain", vmin=vmin_r, vmax=vmax_r)
    try:
        fig.colorbar(im_r, ax=ax_raw, location="left", pad=0.02, fraction=0.04,
                     label="elev (m NAVD88)")
    except (TypeError, ValueError):
        fig.colorbar(im_r, ax=ax_raw, pad=0.01, fraction=0.04,
                     label="elev (m NAVD88)")
    ax_raw.contour(xs, np.arange(n_raw), raw.T, levels=[MHW_M],
                   colors="#1565C0", linewidths=1.0)
    ax_raw.fill_between(xs, i0 + off, i1 + off, color="#FF8C00", alpha=0.22,
                        zorder=0, label="dune search window")
    ax_raw.plot(xs, sb + off, color="k", lw=1.0, label="beach start")
    ax_raw.plot(xs, dl + off, color="#B71C1C", lw=0, marker=".", ms=4,
                label="picked dune crest")
    if res["start_island"] is not None:
        ax_raw.plot(xs, res["start_island"] + off, color="w", lw=1.4, ls="--",
                    label="interior start")
    ax_raw.plot(xs, off - 0.5, color="0.3", lw=1.0, ls="-.",
                label="water trim edge")

    # The road in the RAW frame: unsheared and untrimmed, so this is NC-12's real
    # diagonal across the north-up clip box. Compare it with the same road on
    # panel 2 -- that difference is what the shear removes, and it is the error
    # any raw cross-shore median of the road inherits.
    for _yr, _m in (dom.get("road_raw") or {}).items():
        add_road_plan_overlay(ax_raw, _m, _yr)

    # crop to the island so the diagonal is legible next to the straightened
    # panel; the full raw is mostly sound and shoreface
    lo_r = max(float(np.nanmin(off)) - 4.0, -0.5)
    hi_r = min(float(np.nanmax(off)) + n_cross + 4.0, n_raw - 0.5)
    ax_raw.set_ylim(lo_r, hi_r)
    ax_raw.set_xlabel("alongshore cell")
    ax_raw.set_ylabel("cross-shore cell  (0 = ocean, landward up)")
    ax_raw.set_title(f"1. RAW GIS input   {raw.shape}   m NAVD88\n"
                     f"{dom['name']}", fontsize=11)
    ax_raw.legend(loc="upper right", fontsize=8, framealpha=0.9)

    # ---------------- 2: straightened (the picking frame) ----------------
    zm = masked_profiles(zs)
    im_s = ax_str.imshow(np.ma.masked_invalid(zm.T), aspect="auto",
                         origin="lower",
                         extent=[-0.5, n_along - 0.5, -0.5, n_cross - 0.5],
                         cmap="terrain", vmin=vmin_r - MHW_M,
                         vmax=vmax_r - MHW_M)
    fig.colorbar(im_s, ax=ax_str, pad=0.01, fraction=0.04, label="elev (m MHW)")
    ax_str.axhspan(i0, i1, color="#FF8C00", alpha=0.22, zorder=0,
                   label="dune search window")
    ax_str.plot(xs, sb, color="k", lw=1.0, label="beach start")
    ax_str.plot(xs, dl, color="#B71C1C", lw=0, marker=".", ms=4,
                label="picked dune crest")
    if res["start_island"] is not None:
        ax_str.axhline(res["start_island"], color="w", lw=1.4, ls="--",
                       label="interior start")
    for _yr, _m in (dom.get("road_masks") or {}).items():
        add_road_plan_overlay(ax_str, _m, _yr)
    ax_str.set_xlabel("alongshore cell")
    ax_str.set_ylabel("cross-shore cell  (straightened frame)")
    if STRAIGHTEN:
        t2 = (f"2. STRAIGHTENED — the picking frame   {zs.shape}\n"
              f"{STRAIGHTEN_FIT} fit on {STRAIGHTEN_REF} start  |  "
              f"obliquity {obliq:.1f}°  |  shear 0–{int(np.max(shear))} cells")
    else:
        t2 = (f"2. PROCESSED PROFILES — the picking frame   {zs.shape}\n"
              f"STRAIGHTEN = False, so this is the raw frame minus the trim")
    ax_str.set_title(t2, fontsize=11)
    ax_str.legend(loc="upper right", fontsize=8, framealpha=0.9)

    # ---------------- 3: what CASCADE reads ----------------
    im_i = ax_int.imshow(np.ma.masked_invalid(topo_disp), aspect="auto",
                         origin="lower",
                         extent=[-0.5, n_int_cols - 0.5, -0.5, n_int_rows - 0.5],
                         cmap="terrain", vmin=vmin_r - MHW_M, vmax=vmax_r - MHW_M)
    fig.colorbar(im_i, ax=ax_int, pad=0.01, fraction=0.04, label="elev (m MHW)")

    # The road re-indexed into the SAVED interior grid. This is the panel that
    # answers the question the other two cannot: is NC-12 inside the array
    # CASCADE reads, and at which interior row -- the same row
    # roadway_manager.bulldoze lands on with int(road_setback / dy). A road that
    # falls off the seaward edge here has a negative setback.
    _row0 = np.asarray(res.get("row0_line", []), dtype=int)
    for _yr, _m in (dom.get("road_masks") or {}).items():
        if _row0.size:
            add_road_plan_overlay(
                ax_int,
                processed_road_grid(_m, _row0, n_int_rows, n_int_cols).T,
                _yr, label=False)
    ax_int.set_ylabel("interior row\n(0 = ocean side)")
    ax_int.set_title(f"3. PROCESSED CASCADE input\ninterior {res['topo_file']}  "
                     f"{res['topo_dm'].shape}  (dam on disk)", fontsize=11)
    plt.setp(ax_int.get_xticklabels(), visible=False)

    dune_disp = np.where(dune_m <= SENTINEL_WATER_M + 1e-9, np.nan, dune_m)
    im_d = ax_dun.imshow(np.ma.masked_invalid(dune_disp[None, :]), aspect="auto",
                         origin="lower", extent=[-0.5, ALONG_COLS - 0.5, -0.5, 0.5],
                         cmap="YlOrBr", vmin=0.0,
                         vmax=max(float(np.nanmax(dune_disp)), 0.5))
    fig.colorbar(im_d, ax=ax_dun, pad=0.01, fraction=0.04,
                 label="dune h above berm (m)")
    ax_dun.set_yticks([0])
    ax_dun.set_yticklabels(["dune"])
    ax_dun.set_xlabel("alongshore cell")

    interior = ("sheared per profile" if res["start_island"] is None
                else f"const from row {res['start_island']}")
    road_note = ""
    for _yr in ROAD_YEARS:
        st = (res.get("road_stats") or {}).get(_yr)
        if st and np.isfinite(st["setback_median_m"]):
            road_note += f"   |   NC-12 {_yr} {st['setback_median_m']:+.0f} m"
    fig.suptitle(
        f"{fig_title(stem)}  —  {RUN_NAME}   |   window [{i0}, {i1}] "
        f"({i1 - i0} cells)   |   {res['n_found']}/{res['n_along']} dunes   |   "
        f"mean dune h {res['mean_dune_h_m']:.2f} m   |   "
        f"obliquity {obliq:.1f}°   |   interior {interior}{road_note}",
        fontsize=12,
    )

    fig_dir.mkdir(parents=True, exist_ok=True)
    fig.savefig(fig_dir / f"{fig_stem(stem)}_gis_vs_processed.png", dpi=130,
                bbox_inches="tight")
    plt.close(fig)


# ==============================================================================
# SETTINGS SHEET & JSON I/O
# ==============================================================================

def global_config() -> dict:
    """Every global knob, for the config sheet / manifest / provenance."""
    return {
        "run name": RUN_NAME,
        "version": VERSION,
        "run_time": datetime.now().isoformat(timespec="seconds"),
        "script": Path(__file__).name,
        "root GIS domains": LOAD_PATH.name,
        "run folder": str(RUN_DIR),
        "load_path_full": str(LOAD_PATH),
        "topo_path_full": str(TOPO_SAVE_PATH),
        "dune_path_full": str(DUNE_SAVE_PATH),
        "window_json": str(WINDOW_JSON),
        "MHW_M": MHW_M,
        "BERM_ELEV_NAVD_M": BERM_ELEV_NAVD_M,
        "beach start": BEACH_START_THR_M,
        "WATER_CLAMP_M": WATER_CLAMP_M,
        "SENTINEL_WATER_M": SENTINEL_WATER_M,
        "MIN_DUNE_H_M": MIN_DUNE_H_M,
        "TOPO_ROWS": TOPO_ROWS,
        "ALONG_COLS": ALONG_COLS,
        "OCEAN_LOC": OCEAN_LOC,
        "ALONGSHORE_FLIP": ALONGSHORE_FLIP,
        "DEFAULT_WINDOW_PX": DEFAULT_WINDOW_PX,
        "CLIP_WINDOW_TO_BEACH": CLIP_WINDOW_TO_BEACH,
        "shift interior": False,  # not implemented in this version (see header)
        "constant interior row": USE_CONST_INTERIOR,
        "FILL_MISSING_DUNE": FILL_MISSING_DUNE,
        "STRAIGHTEN": STRAIGHTEN,
        "STRAIGHTEN_REF": STRAIGHTEN_REF if STRAIGHTEN else "",
        "STRAIGHTEN_FIT": STRAIGHTEN_FIT if STRAIGHTEN else "",
        "TRIM_INTERIOR_ROWS": TRIM_INTERIOR_ROWS,
        # Road overlay. Recorded because the road columns in the sheet are
        # meaningless without knowing which vintage and which mask tree produced
        # them -- but note that NONE of these affect the saved arrays.
        "SHOW_ROAD": SHOW_ROAD,
        "ROAD_YEARS": str(ROAD_YEARS) if SHOW_ROAD else "",
        "ROAD_RASTER_ROOT": str(ROAD_RASTER_ROOT) if SHOW_ROAD else "",
        "REQUIRE_ROAD_MASKS": REQUIRE_ROAD_MASKS if SHOW_ROAD else "",
    }


def road_columns(res: dict) -> dict:
    """Per-year NC-12 columns for the settings sheet.

    `road setback <year> (m)` is the cross-check column: seaward edge, median
    over profiles, metres from SAVED interior row 0, positive = landward. It
    should equal `setback_2009_m` in RoadOffset_<year>_domains.csv to the metre.
    See road_offset_stats for why the edge and the statistic are what they are.

    Blank rather than 0 where the road is absent: D1-D7 have no NC-12 at all, and
    a 0 there would read as "road exactly at interior row 0".
    """
    def num(v, nd=1):
        return round(v, nd) if np.isfinite(v) else ""

    out = {}
    for year in ROAD_YEARS:
        st = (res.get("road_stats") or {}).get(year)
        if st is None:
            continue
        out[f"road profiles {year}"] = st["road_profiles"]
        out[f"road cells {year}"] = st["road_cells"]
        out[f"road span cells {year}"] = num(st["road_span_cells"])
        out[f"road width cells {year}"] = num(st["road_width_cells"])
        out[f"road center cell {year}"] = num(st["road_center_cell"])
        # the comparable one, named for what it is
        out[f"road setback {year} (m)"] = num(st["setback_median_m"])
        out[f"road setback mean {year} (m)"] = num(st["setback_mean_m"])
        out[f"road setback min {year} (m)"] = num(st["setback_min_m"])
        out[f"road setback max {year} (m)"] = num(st["setback_max_m"])
        # centre-referenced, for continuity with the roya-style dune-to-road
        # number; NOT the value bulldoze indexes
        out[f"road center offset {year} (m)"] = num(st["center_median_m"])
        # Profiles where the road is SEAWARD of interior row 0. Barrier3D cannot
        # represent that (int(negative/dy) indexes from the landward end), so it
        # is the flag the setback audit's NEGATIVE floor exists to handle.
        out[f"road seaward profiles {year}"] = st["n_seaward"]
    return out


def settings_row(res: dict, w: dict | None) -> dict:
    """One sheet row per domain: settings used + what came out."""
    return {
        "domain": domain_number(res["stem"]),
        "stem": res["stem"],
        "section": section_for(res["stem"]),
        # --- settings, mirroring the tracking sheet ---
        "root GIS domains": LOAD_PATH.name,
        "root dunes/topo": RUN_NAME,
        "figure name": f"{fig_stem(res['stem'])}_gis_vs_processed.png",
        "beach start": BEACH_START_THR_M,
        "dune window start": res["i0"],
        "dune window end": res["i1"],
        "dune window": res["i1"] - res["i0"],
        "window source": "picked" if w else "default",
        "clip window to beach": CLIP_WINDOW_TO_BEACH,
        "shift interior": False,
        "constant interior row": USE_CONST_INTERIOR,
        "straightened": res.get("straightened", ""),
        "straighten fit": STRAIGHTEN_FIT if STRAIGHTEN else "",
        "obliquity (deg)": res.get("obliquity_deg", ""),
        "shear max (cells)": res.get("shear_max_cells", ""),
        "shear max (m)": (round(res.get("shear_max_cells", 0) * CELL_SIZE_M, 1)
                          if res.get("shear_max_cells") else ""),
        "MHW (m NAVD88)": MHW_M,
        "berm (m NAVD88)": BERM_ELEV_NAVD_M,
        "water clamp (m MHW)": WATER_CLAMP_M,
        # --- results ---
        "interior start row": res["start_island"],
        "lead trim rows": res.get("lead_trim_rows", ""),
        "dunes found": res["n_found"],
        "profiles": res["n_along"],
        "dunes filled": res["n_filled"],
        "mean dune h (m)": round(res["mean_dune_h_m"], 3),
        "min dune h (m)": round(res["min_dune_h_m"], 3),
        "max dune h (m)": round(res["max_dune_h_m"], 3),
        "interior rows": res["interior_rows"],
        "interior cols": res["interior_cols"],
        "mean interior elev (m MHW)": round(res["mean_interior_elev_m"], 3),
        "topo file": res["topo_file"],
        "dune file": res["dune_file"],
        "picked": (w or {}).get("picked", ""),
        **road_columns(res),
    }


def write_settings_sheet(rows: list, base_path: Path) -> None:
    """Write the per-domain settings sheet as CSV, plus XLSX if pandas is around."""
    if not rows:
        return
    base_path.parent.mkdir(parents=True, exist_ok=True)

    # Union of keys in first-seen order, not rows[0].keys(). With
    # REQUIRE_ROAD_MASKS = False a domain whose mask is missing contributes no
    # road columns, and DictWriter raises on any row holding a key the header
    # does not -- so keying off the first row would turn one absent mask into a
    # crash after all 90 domains had been processed.
    fieldnames = list(dict.fromkeys(k for r in rows for k in r))
    csv_path = base_path.with_suffix(".csv")
    with open(csv_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, restval="")
        writer.writeheader()
        writer.writerows(rows)
    print(f"[sheet] {csv_path}")

    try:
        import pandas as pd
        cfg = global_config()
        cfg_df = pd.DataFrame({"setting": list(cfg.keys()),
                               "value": [str(v) for v in cfg.values()]})
        xlsx_path = base_path.with_suffix(".xlsx")
        with pd.ExcelWriter(xlsx_path, engine="openpyxl") as xw:
            pd.DataFrame(rows).to_excel(xw, sheet_name="domains", index=False)
            cfg_df.to_excel(xw, sheet_name="global_config", index=False)
        print(f"[sheet] {xlsx_path}")
    except ImportError:
        print("[info] pandas/openpyxl not available; wrote CSV only")


def write_manifest(rows: list, path: Path) -> None:
    """Plain-text record so a run folder found later explains itself."""
    path.parent.mkdir(parents=True, exist_ok=True)
    cfg = global_config()
    width = max(len(k) for k in cfg)
    lines = [
        "=" * 78,
        f"HAT dune & topography extraction  --  run {RUN_NAME}",
        "=" * 78,
        "",
        "CONTENTS",
        "  topography\\            interior elevation arrays (dam), CASCADE input",
        "  dunes\\                 dune height above berm (dam), CASCADE input",
        "  figures\\gis_vs_processed\\   raw DEM vs what CASCADE receives, per domain",
        "  figures\\qc\\            dune detection QC, per domain",
        f"  {SHEET_SAVE_PATH.name}.xlsx / .csv   per-domain settings + results",
        f"  {SUMMARY_FIG_PATH.name}   all domains on one page",
        f"  {PLAN_FIG_STEM}_<year>_<mode>.png   plan view, per year per mode",
        "",
        "PICKS (not in this folder -- not regenerable)",
        f"  {WINDOW_JSON}",
    ]
    if PICK_SET == RUN_NAME:
        lines += [
            "  ^ THIS VERSION CARRIES ITS OWN PICKS. It was seeded by copying the",
            "    shared straight set, then re-picked with NC-12 drawn on the",
            "    picker, so the shared set is untouched and earlier versions stay",
            "    reproducible.",
        ]
    else:
        lines += ["  ^ shared across versions -- re-picking here rewrites it for "
                  "every version that points at it"]
    if SHOW_ROAD:
        lines += [
            "",
            "ROAD OVERLAY (display + diagnostics only)",
            f"  masks read from : {ROAD_RASTER_ROOT}",
            f"  vintages        : {', '.join(str(y) for y in ROAD_YEARS)}",
            "  The road does NOT enter find_dunes, build_interior or",
            "  straighten_profiles. topography\\, dunes\\ and the nodata masks are",
            "  byte-identical with SHOW_ROAD either way; only the figures and the",
            "  road columns of the settings sheet change.",
            "  Road distances are metres from SAVED interior row 0, positive =",
            "  landward, matching setback_2009_m in RoadOffset_<year>_domains.csv.",
        ]
    lines += [
        "",
        "SETTINGS",
    ]
    lines += [f"  {k:<{width}} : {v}" for k, v in cfg.items()]
    lines += ["", f"DOMAINS WRITTEN: {len(rows)}"]
    if rows:
        nums = [r["domain"] for r in rows if r["domain"] is not None]
        if nums:
            lines.append(f"  range: {min(nums)} - {max(nums)}")
        defaults = [r["stem"] for r in rows if r["window source"] == "default"]
        if defaults:
            lines.append(f"  !! fallback default window (never picked): "
                         f"{', '.join(defaults)}")
        filled = [r["stem"] for r in rows if r["dunes filled"] > 0]
        if filled:
            lines.append(f"  !! profiles with no dune found, filled with "
                         f"{MIN_DUNE_H_M} m: {', '.join(filled)}")
    lines.append("")
    path.write_text("\n".join(lines))
    print(f"[manifest] {path}")


def summary_figure(rows: list, path: Path) -> None:
    """Every domain on one page: windows, dune heights, interior extent."""
    if not rows:
        return
    rows = sorted([r for r in rows if r["domain"] is not None],
                  key=lambda r: r["domain"])
    if not rows:
        return
    d = np.array([r["domain"] for r in rows])

    fig, (ax0, ax1, ax2) = plt.subplots(3, 1, figsize=(13, 9), sharex=True,
                                        gridspec_kw={"hspace": 0.12})

    # section bands + labels, matching the HAT_hindcast annotation style
    for k, ((lo, hi), label) in enumerate(SECTIONS):
        if not ((d >= lo) & (d <= hi)).any():
            continue
        for ax in (ax0, ax1, ax2):
            ax.axvspan(lo - 0.5, hi + 0.5, color="#90AFC5",
                       alpha=0.18 if k % 2 == 0 else 0.08, lw=0, zorder=0)
        trans = blended_transform_factory(ax0.transData, ax0.transAxes)
        ax0.text((lo + hi) / 2, 0.96, label, transform=trans, ha="center",
                 va="top", fontsize=8, rotation=0,
                 bbox=dict(boxstyle="round,pad=0.2", fc="white", ec="none",
                           alpha=0.85))

    # 1) search window and interior start
    ax0.fill_between(d, [r["dune window start"] for r in rows],
                     [r["dune window end"] for r in rows],
                     color="#FF8C00", alpha=0.45, step="mid",
                     label="dune search window")
    ax0.step(d, [r["interior start row"] for r in rows], where="mid",
             color="#B71C1C", lw=1.4, label="interior start row")

    # NC-12's median cross-shore cell per domain, on the same axis as the search
    # window. Where the road line dips INTO the orange band, that domain's window
    # and the road overlap -- the crest argmax could be locking onto the road.
    for _yr in ROAD_YEARS:
        key = f"road center cell {_yr}"
        if not any(isinstance(r.get(key), (int, float)) for r in rows):
            continue
        vals = np.array([r.get(key) if isinstance(r.get(key), (int, float))
                         else np.nan for r in rows], dtype=float)
        if np.isfinite(vals).any():
            ax0.plot(d, vals, color=ROAD_COLORS.get(_yr, "#111111"), lw=1.2,
                     ls="--", label=f"NC-12 {_yr} centre")
    ax0.set_ylabel("cross-shore cell\n(0 = ocean)")
    ax0.legend(loc="lower right", fontsize=8, framealpha=0.9)
    ax0.set_title(f"HAT dune & topography extraction  —  {RUN_NAME}  —  "
                  f"{len(rows)} domains", fontsize=12)

    # 2) dune height
    ax1.fill_between(d, [r["min dune h (m)"] for r in rows],
                     [r["max dune h (m)"] for r in rows],
                     color="#FF8C00", alpha=0.25, label="min-max")
    ax1.plot(d, [r["mean dune h (m)"] for r in rows], color="#FF8C00", lw=1.6,
             marker="o", ms=3, label="mean")
    ax1.axhline(MIN_DUNE_H_M, color="0.4", ls="--", lw=1.0,
                label=f"floor {MIN_DUNE_H_M} m")
    ax1.set_ylabel("dune height\nabove berm (m)")
    ax1.legend(loc="upper right", fontsize=8, framealpha=0.9)

    # 3) interior extent + data-quality flags
    ax2.plot(d, [r["interior rows"] for r in rows], color="#1565C0", lw=1.6,
             marker="o", ms=3, label="interior rows")
    ax2.set_ylabel("interior rows")
    ax2.set_xlabel("domain  (1 = Cape Point / south  →  90 = Rodanthe / north)")

    flag = np.array([r["dunes filled"] > 0 or r["window source"] == "default"
                     for r in rows])
    if flag.any():
        trans2 = blended_transform_factory(ax2.transData, ax2.transAxes)
        ax2.plot(d[flag], np.full(flag.sum(), 0.04), transform=trans2, lw=0,
                 marker="v", ms=6, color="#B71C1C", clip_on=False,
                 label="default window / filled dunes")
    ax2.legend(loc="upper right", fontsize=8, framealpha=0.9)
    ax2.set_xlim(d.min() - 0.5, d.max() + 0.5)

    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=140, bbox_inches="tight")
    plt.close(fig)
    print(f"[summary] {path}")


def resolve_offset_file(year: int, configured) -> Path | None:
    """
    Return the raw_offset CSV for a year. If the configured path is wrong, search
    OFFSET_DIR rather than silently skipping -- the hindcast_* subfolder naming
    is a convention this script only partly knows.
    """
    configured = Path(configured)
    if configured.exists():
        return configured
    root = Path(OFFSET_DIR)
    if not root.is_dir():
        return None
    hits = sorted(p for p in root.rglob("Island_Dune_Offsets*.csv")
                  if str(year) in p.name)
    if len(hits) == 1:
        print(f"[path] offsets {year}: configured path not found, found instead")
        print(f"       {hits[0]}")
        print(f"       -> update OFFSET_FILES[{year}] to match")
        return hits[0]
    if len(hits) > 1:
        print(f"[path] offsets {year}: configured path not found; "
              f"{len(hits)} candidates, none picked (ambiguous):")
        for h in hits:
            print(f"       {h}")
    return None


def check_paths() -> bool:
    """Preflight every path before doing any work. False = don't run."""
    print("=" * 78)
    print(f"PATH CHECK  —  run {RUN_NAME}")
    print("=" * 78)
    ok = True

    def row(tag, label, path, note=""):
        print(f"  [{tag:^7}] {label:<20} {path}{note}")

    for label, pth in (("project root", PROJECT_ROOT), ("hatteras_init", INIT_ROOT)):
        if Path(pth).is_dir():
            row("ok", label, pth)
        else:
            row("MISSING", label, pth)
            ok = False

    if Path(LOAD_PATH).is_dir():
        n = len([f for f in os.listdir(LOAD_PATH)
                 if f.startswith("domain_") and f.endswith(".npy")])
        if n:
            row("ok", "DEM domains", LOAD_PATH, f"   ({n} domain_*.npy)")
        else:
            row("EMPTY", "DEM domains", LOAD_PATH, "   no domain_*.npy here")
            ok = False
    else:
        row("MISSING", "DEM domains", LOAD_PATH)
        ok = False

    for year in sorted(OFFSET_FILES):
        found = resolve_offset_file(year, OFFSET_FILES[year])
        if found is None:
            row("MISSING", f"offsets {year}", OFFSET_FILES[year],
                "   (figures will skip this year)")
        elif Path(found) != Path(OFFSET_FILES[year]):
            row("FOUND", f"offsets {year}", found, "   (not where configured)")
        else:
            row("ok", f"offsets {year}", found)

    if SHOW_ROAD:
        for year in ROAD_YEARS:
            mdir = Path(ROAD_RASTER_ROOT) / ROAD_MASK_DIR_FMT.format(year=year)
            if not mdir.is_dir():
                row("MISSING", f"road masks {year}", mdir,
                    "   run HAT_rasterize_road_to_domains.py")
                if REQUIRE_ROAD_MASKS:
                    ok = False
                continue
            wanted = PICK_DOMAINS if PICK_DOMAINS is not None \
                else list(range(1, NUM_REAL_DOMAINS + 1))
            absent = [n for n in wanted if not road_mask_path(year, n).is_file()]
            if absent:
                row("MISSING" if REQUIRE_ROAD_MASKS else "WARN",
                    f"road masks {year}", mdir,
                    f"   {len(wanted) - len(absent)}/{len(wanted)} present, "
                    f"missing {absent[:8]}{' ...' if len(absent) > 8 else ''}")
                if REQUIRE_ROAD_MASKS:
                    ok = False
            else:
                row("ok", f"road masks {year}", mdir,
                    f"   ({len(wanted)} masks)")

    wj = Path(WINDOW_JSON)
    if wj.exists():
        try:
            n = len([k for k in json.loads(wj.read_text()) if not k.startswith("_")])
            row("ok", "picks", wj, f"   ({n} domains picked)")
        except Exception as e:
            row("BAD", "picks", wj, f"   unreadable: {e}")
            ok = False
    else:
        row("new", "picks", wj, "   (created on first pick)")

    # save_windows() rewrites WINDOW_JSON after EVERY domain, so a picking run
    # against a SHARED pick set destroys the picks of every version pointing at
    # it, one domain at a time, with no prompt. That is how 2009_v3 would have
    # been lost to the v4 re-pick. Loud, but not fatal: sharing a pick set is
    # legitimate when you are only topping up domains that were never picked.
    if MODE in ("pick", "pick_and_run") and PICK_SET != RUN_NAME:
        row("WARN", "pick set", wj,
            f"\n            ^ MODE = {MODE!r} will REWRITE this shared pick set "
            f"({PICK_SET!r}).\n"
            f"              Every version that reads it changes with it. To give "
            f"{RUN_NAME} its own,\n"
            f"              set PICK_SET = RUN_NAME and copy the old file to "
            f"HAT_dune_search_windows_{RUN_NAME}.json first.")

    if MODE in ("pick", "pick_and_run") and not REPICK_EXISTING and wj.exists():
        try:
            n_saved = len([k for k in json.loads(wj.read_text())
                           if not k.startswith("_")])
            n_want = (len(PICK_DOMAINS) if PICK_DOMAINS is not None else n_saved)
            if n_saved >= n_want:
                row("WARN", "pick pass", wj,
                    f"\n            ^ REPICK_EXISTING = False and all {n_want} "
                    f"requested domains are already saved,\n"
                    f"              so the pick pass will skip every one of them. "
                    f"Set REPICK_EXISTING = True to re-draw.")
        except Exception:
            pass

    for legacy in (Path(RUN_DIR) / f"{PLAN_FIG_STEM}.png",):
        if legacy.exists():
            row("STALE", "original plan view", legacy,
                "\n            ^ pre-dates the per-year/per-mode naming, nothing "
                "overwrites it — delete it")

    rd = Path(RUN_DIR)
    row("exists" if rd.exists() else "new", "run folder", rd,
        "   (will be overwritten)" if rd.exists() else "   (created on run)")

    print("=" * 78)
    if not ok:
        print("[stop] required inputs missing — fix the paths above and re-run.\n")
    return ok


def load_offsets() -> dict:
    """{year: (domain_numbers, offset_m)} from the CASCADE island-raw_offset CSVs."""
    out = {}
    for year, configured in OFFSET_FILES.items():
        path = resolve_offset_file(year, configured)
        if path is None:
            print(f"[warn] raw_offset file not found, skipping {year}: {configured}")
            continue
        v = np.loadtxt(path, skiprows=1, delimiter=",", ndmin=2).astype(float)
        if v.shape[1] > 1:
            print(f"[raw_offset] {year}: {v.shape[1]} columns, using column {OFFSET_COLUMN}")
            v = v[:, OFFSET_COLUMN]
        else:
            v = v[:, 0]
        padded = NUM_REAL_DOMAINS + 2 * N_BUFFER_DOMAINS
        if v.size == padded:
            print(f"[raw_offset] {year}: padded file ({padded} rows), stripping "
                  f"{N_BUFFER_DOMAINS} buffer domains from each end")
            v = v[N_BUFFER_DOMAINS:N_BUFFER_DOMAINS + NUM_REAL_DOMAINS]
        n = v.size
        if OFFSET_ROW_ORDER == "D1_first":
            dom = np.arange(n) + 1
        else:
            dom = n - np.arange(n)
        order = np.argsort(dom)
        out[year] = (dom[order], v[order])
        print(f"[raw_offset] {year}: {n} domains, {v.min():.0f}-{v.max():.0f} m "
              f"({path.name})")
    return out


def _island_norm():
    """Poster colormap: terrain with 0 m pinned to colormap position 0.35."""
    lo, hi, pos = ISLAND_ELEV_MIN_M, ISLAND_ELEV_MAX_M, ISLAND_SEA_LEVEL_POS

    def fwd(x):
        out = np.where(x < 0.0, pos * (x - lo) / (0.0 - lo),
                       pos + (1.0 - pos) * x / hi)
        return np.where(np.isnan(x), np.nan, out)

    def inv(x):
        return np.where(x < pos, lo + (x / pos) * (0.0 - lo),
                        (x - pos) / (1.0 - pos) * hi)

    cmap = plt.cm.terrain.copy()
    cmap.set_bad(color=ISLAND_OCEAN_COLOR)
    return cmap, FuncNorm((fwd, inv), vmin=lo, vmax=hi)


def _assert_alongshore_continuity(grids, label: str, warn_ratio: float = 5.0):
    """Warn if the assembled alongshore axis is discontinuous at domain seams.

    A per-domain alongshore reversal is almost invisible in a 45 km plan view --
    it reads as roughness -- but it puts every 500 m block backwards. The signal
    is unmistakable in numbers: compare the mean jump ACROSS domain seams with the
    mean jump WITHIN a domain. A continuous island sits near 1; a per-domain
    reversal drove this to 21 on the 2009_v3 arrays when the legacy plotting flip
    was still applied, which is the bug this guard exists to catch.

    Returns the ratio, or nan when it cannot be computed.
    """
    series = []
    for g in grids:
        land = (g > SENTINEL_WATER_M + 1e-9).sum(axis=0).astype(float)
        if land.size:
            series.append(land)
    if len(series) < 3:
        return float("nan")

    width = min(s.size for s in series)
    stack = np.array([s[:width] for s in series])
    seam = np.abs(stack[:-1, -1] - stack[1:, 0])
    inner = np.abs(np.diff(stack, axis=1))
    if not inner.size or inner.mean() <= 0:
        return float("nan")

    ratio = float(seam.mean() / inner.mean())
    if ratio > warn_ratio:
        print(f"       [ALONGSHORE WARNING] {label}: seam/inner discontinuity "
              f"ratio {ratio:.1f} (expected ~1-2).")
        print(f"       The alongshore axis is reversed WITHIN each domain "
              f"relative to the domain order, so every 500 m block is backwards.")
        print(f"       Either the arrays were extracted with ALONGSHORE_FLIP = "
              f"False, or a per-domain np.fliplr has been reintroduced here.")
    return ratio


def _build_island_canvas(recs, offset_m_by_domain, mode):
    """Stitch processed domains onto one plan-view canvas at their dune offsets."""
    use = [(n, r) for n, r in recs if n in offset_m_by_domain]
    if not use:
        return None, None, None, None
    off_cells = [int(round(offset_m_by_domain[n] / CELL_SIZE_M)) for n, _ in use]

    grids, dunes, cropped = [], [], []
    roads = {yr: [] for yr in ROAD_YEARS}
    for n_dom, r in use:
        g = r["topo_dm"] * CELL_SIZE_M                 # dam -> m MHW

        # The road in the SAVED interior frame, so it is padded, cropped and
        # placed by exactly the same rules as the topography it sits on. Built
        # before the pad/crop below so it goes through both with the grid.
        _row0 = np.asarray(r.get("row0_line", []), dtype=int)
        for yr in ROAD_YEARS:
            m = (r.get("road_masks") or {}).get(yr)
            roads[yr].append(
                processed_road_grid(m, _row0, g.shape[0], g.shape[1])
                if (m is not None and _row0.size)
                else np.zeros(g.shape, dtype=bool))

        if mode == "padded":
            # pad landward, matching where dune_topo_extractor_from_GIS.py left its
            # sentinel: interior row 0 is the ocean side, rows increase landward
            n = g.shape[0]
            if n < ISLAND_PAD_ROWS:
                g = np.vstack([g, np.full((ISLAND_PAD_ROWS - n, g.shape[1]),
                                          SENTINEL_WATER_M)])
                for yr in ROAD_YEARS:
                    roads[yr][-1] = np.vstack([
                        roads[yr][-1],
                        np.zeros((ISLAND_PAD_ROWS - n, g.shape[1]), dtype=bool)])
            elif n > ISLAND_PAD_ROWS:
                n_land = int(np.sum(g[ISLAND_PAD_ROWS:] > SENTINEL_WATER_M + 1e-9))
                if n_land:
                    cropped.append((n_dom, n, n_land))
                g = g[:ISLAND_PAD_ROWS]
                # The road can be cropped off entirely here: NC-12 sits well
                # landward on the wide domains, so ISLAND_PAD_ROWS = 100 loses it
                # where it also loses real land. That is the same loss the
                # `cropped` warning already reports, not a separate bug.
                for yr in ROAD_YEARS:
                    roads[yr][-1] = roads[yr][-1][:ISLAND_PAD_ROWS]
        if ISLAND_SENTINEL_AS_OCEAN:
            g = np.where(g <= SENTINEL_WATER_M + 1e-9, np.nan, g)
        d = r["dune_dm"] * CELL_SIZE_M + (BERM_ELEV_NAVD_M - MHW_M)   # -> m MHW
        # No per-domain flip here: the arrays already run south -> north within a
        # domain, matching the ascending domain order. See the removal note at
        # ISLAND_INCLUDE_DUNE.
        grids.append(g)
        dunes.append(d)

    if cropped:
        print(f"[planview] ISLAND_PAD_ROWS={ISLAND_PAD_ROWS} "
              f"({ISLAND_PAD_ROWS * CELL_SIZE_M:.0f} m) crops real land from "
              f"{len(cropped)} domain(s):")
        for n_dom, n_rows, n_land in cropped[:8]:
            print(f"           D{n_dom}: {n_rows} rows stored, {n_land} land cells lost")
        if len(cropped) > 8:
            print(f"           ... and {len(cropped) - 8} more")

    _assert_alongshore_continuity(grids, f"planview canvas ({mode})")

    max_rows = max(g.shape[0] for g in grids)
    canvas_rows = max(off_cells) + max_rows + 5
    total_cols = sum(g.shape[1] for g in grids)
    canvas = np.full((canvas_rows, total_cols), np.nan)

    road_canvas = {yr: np.zeros((canvas_rows, total_cols), dtype=bool)
                   for yr in ROAD_YEARS}

    col = 0
    starts = []
    for k, (g, d) in enumerate(zip(grids, dunes)):
        n_rows, n_cols = g.shape
        starts.append(col)
        origin = off_cells[k]
        end = min(origin + n_rows, canvas_rows)
        canvas[origin:end, col:col + n_cols] = g[:end - origin, :]
        if ISLAND_INCLUDE_DUNE and origin >= 1:
            canvas[origin - 1, col:col + min(n_cols, d.size)] = d[:n_cols]
        # Same origin, same columns, same clip as the grid above -- the road is
        # placed by the topography's rule, not its own.
        for yr in ROAD_YEARS:
            rg = roads[yr][k]
            if rg.shape[0] >= end - origin:
                road_canvas[yr][origin:end, col:col + n_cols] = rg[:end - origin, :]
        col += n_cols

    return canvas, np.array(starts), [n for n, _ in use], road_canvas


def island_plan_figure(summary: list, offsets: dict, run_dir: Path) -> None:
    """
    Plan view of the processed dune + interior for domains 1-90 at the measured
    offsets. Writes ONE FIGURE PER OFFSET YEAR PER CROSS-SHORE MODE, styled to
    match plot_initialization_poster_no_border.py.
    """
    recs = sorted([(domain_number(r["stem"]), r) for r in summary
                   if domain_number(r["stem"]) is not None])
    if not recs or not offsets:
        return
    cmap, norm = _island_norm()

    for year in sorted(offsets):
        dom, v = offsets[year]
        omap = {int(a): float(b) for a, b in zip(dom, v)}
        for mode in ISLAND_CROSS_SHORE_MODES:
            canvas, starts, used, road_canvas = _build_island_canvas(
                recs, omap, mode)
            if canvas is None:
                print(f"[planview] {year}/{mode}: no domains overlap the offsets, skipped")
                continue

            n_cs, n_al = canvas.shape
            fig_w = 20.0
            fig_h = min(max(4.5, fig_w * (n_cs / n_al) * 1.8), 7.5)   # poster aspect
            fig = plt.figure(figsize=(fig_w, fig_h), facecolor="white")
            ax = fig.add_axes([0.06, 0.18, 0.88, 0.68])
            ax.set_facecolor(ISLAND_OCEAN_COLOR)

            im = ax.pcolormesh(np.ma.masked_invalid(canvas), cmap=cmap, norm=norm,
                               shading="auto", rasterized=True)

            # NC-12 across the whole island, in the offset frame. At 45 km wide a
            # single 20 m road is around one pixel, so it is drawn with
            # pcolormesh on the same canvas rather than as a line: that way it
            # cannot drift relative to the topography it was placed against.
            for _yr in ROAD_YEARS:
                rc = (road_canvas or {}).get(_yr)
                if rc is None or not rc.any():
                    continue
                ax.pcolormesh(np.ma.masked_where(~rc, rc.astype(float)),
                              cmap=ListedColormap([ROAD_COLORS.get(_yr, "#111111")]),
                              vmin=0.0, vmax=1.0, shading="auto",
                              rasterized=True, zorder=4)
                ax.plot([], [], color=ROAD_COLORS.get(_yr, "#111111"), lw=3,
                        label=f"NC-12 {_yr}")
            if any((road_canvas or {}).get(y) is not None
                   and (road_canvas or {})[y].any() for y in ROAD_YEARS):
                ax.legend(loc="upper right", fontsize=9, framealpha=0.9)

            ax.set_xlim(0, n_al)
            ax.set_ylim(0, n_cs)

            cax = fig.add_axes([0.955, 0.18, 0.013, 0.68])
            cbar = plt.colorbar(im, cax=cax)
            cbar.set_label("Elevation (m MHW)", fontsize=12, color="#1a1a2e",
                           labelpad=10, rotation=270)
            cbar.ax.yaxis.set_tick_params(color="#1a1a2e", labelcolor="#1a1a2e")
            cbar.outline.set_edgecolor("#cccccc")
            cbar.set_ticks([-1, 0, 1, 2, 3, 4])

            ticks, labels = [], []
            for k, n in enumerate(used):
                if n % 5 == 0 or n == 1:
                    ticks.append(starts[k] + ALONG_COLS // 2)
                    labels.append(str(n))
            ax.set_xticks(ticks)
            ax.set_xticklabels(labels, fontsize=9)
            ax.set_xlabel("Domain (S → N,  Cape Hatteras to Rodanthe)", fontsize=12,
                          labelpad=8)
            ax.set_ylabel("Cross-shore cell (raw_offset frame)", fontsize=12)
            for k, n in enumerate(used):
                if n % 10 == 0:
                    ax.axvline(starts[k] - 0.5, color="#aaaaaa", lw=0.4, alpha=0.5,
                               zorder=2)
            for sp in ("top", "right"):
                ax.spines[sp].set_visible(False)
            for sp in ("bottom", "left"):
                ax.spines[sp].set_color("#999999")

            what = "dune + interior" if ISLAND_INCLUDE_DUNE else "interior"
            if mode == "padded":
                extent = f"{ISLAND_PAD_ROWS} cells / {ISLAND_PAD_ROWS * CELL_SIZE_M:.0f} m"
                note = f"{what}, every domain padded to {extent} cross-shore"
            else:
                note = f"{what}, each domain trimmed to its own island"
            ax.set_title(f"Hatteras Island — CASCADE Initialization  |  {year} offsets  "
                         f"|  {DEM_YEAR} extracted {note}  ({len(used)} domains)",
                         fontsize=14, fontweight="bold", color="#1a1a2e", pad=12)

            path = Path(run_dir) / f"{PLAN_FIG_STEM}_{year}_{mode}.png"
            path.parent.mkdir(parents=True, exist_ok=True)
            fig.savefig(path, dpi=200, bbox_inches="tight", facecolor="white")
            plt.close(fig)
            print(f"[planview] {path}")


def island_figure(summary: list, offsets: dict, path: Path) -> None:
    """
    All domains together: measured dune offsets vs the crest this run extracted.

    The measured offsets and the extracted crest live in DIFFERENT frames -- the
    offsets are in the model's common cross-shore frame, the extracted crest is in
    the per-domain raw DEM array frame (m landward of that array's cell 0). They
    are NOT differenced here. They're plotted on separate panels, and the script
    reports the correlation between them so the frame relationship is testable
    rather than assumed.
    """
    recs = sorted([(domain_number(r["stem"]), r) for r in summary
                   if domain_number(r["stem"]) is not None])
    if not recs:
        return
    d = np.array([n for n, _ in recs])

    pos_mean, pos_lo, pos_hi = [], [], []
    for _, r in recs:
        loc = r["dune_loc"].astype(float)
        loc[loc < 0] = np.nan
        # back to the raw cross-shore axis: k + c0 + shear[i] (shear is zeros
        # when STRAIGHTEN is False, so this is the original `+ c0`)
        _sh = np.asarray(r.get("shear", 0))
        pm = (loc + r["c0"] + _sh) * CELL_SIZE_M
        pos_mean.append(np.nanmean(pm))
        pos_lo.append(np.nanmin(pm))
        pos_hi.append(np.nanmax(pm))
    pos_mean = np.array(pos_mean)

    fig, (ax0, ax1, ax2) = plt.subplots(3, 1, figsize=(13, 10), sharex=True,
                                        gridspec_kw={"hspace": 0.12})

    for k, ((lo, hi), label) in enumerate(SECTIONS):
        for ax in (ax0, ax1, ax2):
            ax.axvspan(lo - 0.5, hi + 0.5, color="#90AFC5",
                       alpha=0.18 if k % 2 == 0 else 0.08, lw=0, zorder=0)
        trans = blended_transform_factory(ax0.transData, ax0.transAxes)
        ax0.text((lo + hi) / 2, 0.97, label, transform=trans, ha="center",
                 va="top", fontsize=8,
                 bbox=dict(boxstyle="round,pad=0.2", fc="white", ec="none",
                           alpha=0.85))

    # 1) measured offsets, all domains
    colors = {1984: "#1565C0", 2004: "#B71C1C"}
    for year in sorted(offsets):
        dom, v = offsets[year]
        ax0.plot(dom, v, color=colors.get(year, "0.4"), lw=1.6,
                 label=f"measured dune raw_offset {year}")
    ax0.set_ylabel("measured raw_offset\n(m, common frame)")
    ax0.legend(loc="upper right", fontsize=8, framealpha=0.9)
    sign = "seaward +" if OFFSET_SEAWARD_POSITIVE else "landward +"
    ax0.set_title(f"Island dune offsets vs extracted crest  —  {RUN_NAME}  "
                  f"({len(recs)} domains, row0={OFFSET_ROW_ORDER}, {sign})",
                  fontsize=12)

    # 2) crest extracted from the DEM this run, raw-array frame
    ax1.fill_between(d, pos_lo, pos_hi, color="#FF8C00", alpha=0.25,
                     label="within-domain min-max")
    ax1.plot(d, pos_mean, color="#FF8C00", lw=1.6, marker="o", ms=3,
             label=f"extracted crest, {DEM_YEAR} DEM")
    ax1.set_ylabel("extracted crest\n(m from raw cell 0)")
    ax1.legend(loc="upper right", fontsize=8, framealpha=0.9)

    # correlation against each measured year -- diagnostic for the frame relation
    notes = []
    for year in sorted(offsets):
        dom, v = offsets[year]
        common = np.intersect1d(d, dom)
        if common.size < 3:
            continue
        a = pos_mean[np.isin(d, common)]
        b = v[np.isin(dom, common)]
        ok = np.isfinite(a) & np.isfinite(b)
        if ok.sum() < 3:
            continue
        r = float(np.corrcoef(a[ok], b[ok])[0, 1])
        notes.append(f"r(extracted, {year}) = {r:+.2f}  (n={ok.sum()})")
    if notes:
        trans = blended_transform_factory(ax1.transAxes, ax1.transAxes)
        ax1.text(0.01, 0.06, "   |   ".join(notes), transform=trans, fontsize=8,
                 va="bottom", ha="left",
                 bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="0.7",
                           alpha=0.9))

    # 3) measured change
    if 1984 in offsets and 2004 in offsets:
        dom_a, va = offsets[1984]
        dom_b, vb = offsets[2004]
        common = np.intersect1d(dom_a, dom_b)
        ch = (vb[np.isin(dom_b, common)] - va[np.isin(dom_a, common)]) / 20.0
        ax2.axhline(0, color="0.5", lw=1.0)
        ax2.plot(common, ch, color="#B71C1C", lw=1.4, marker="o", ms=3)
        ax2.set_ylabel("measured change\n1984→2004 (m/yr)")
    else:
        ax2.text(0.5, 0.5, "need both 1984 and 2004 offsets", ha="center",
                 va="center", transform=ax2.transAxes, fontsize=10, color="0.4")
    ax2.set_xlabel("domain  (1 = Cape Point / south  →  90 = Rodanthe / north)")
    ax2.set_xlim(min(d.min(), 1) - 0.5, max(d.max(), 90) + 0.5)

    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=140, bbox_inches="tight")
    plt.close(fig)
    print(f"[island] {path}")
    for n in notes:
        print(f"         {n}")


def load_windows(path: Path) -> dict:
    if path.exists():
        with open(path) as f:
            return json.load(f)
    return {"_meta": {}}


def save_windows(path: Path, windows: dict) -> None:
    windows["_meta"] = {
        "updated": datetime.now().isoformat(timespec="seconds"),
        "load_path": str(LOAD_PATH),
        "mhw_m": MHW_M,
        "berm_elev_navd_m": BERM_ELEV_NAVD_M,
        "beach_start_thr_m": BEACH_START_THR_M,
        "water_clamp_m": WATER_CLAMP_M,
        "ocean_loc": OCEAN_LOC,
        "alongshore_flip": ALONGSHORE_FLIP,
        "straighten": STRAIGHTEN,
        "straighten_fit": STRAIGHTEN_FIT if STRAIGHTEN else "",
        "note": "i0/i1 are cross-shore indices in the ocean-first, "
                "water-trimmed profile array. If straighten is true they are "
                "in the STRAIGHTENED frame and do not apply to an "
                "unstraightened array -- the index range is still valid, it "
                "just points at different cells. alongshore_flip is recorded "
                "for provenance only: the flip reverses the alongshore axis "
                "but leaves every profile's cross-shore frame bit-identical "
                "(shear[i] -> shear[n-1-i], c0 unchanged), so these windows "
                "are valid under either setting. Verified for all 90 domains "
                "by HAT_alongshore_frame_check.py.",
    }
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w") as f:
        json.dump(windows, f, indent=2, sort_keys=True)


# ==============================================================================
# MAIN
# ==============================================================================

def main():
    print("=" * 86)
    print(f"HAT dune / topography extraction  |  {RUN_NAME}  |  MODE = {MODE}")
    print("=" * 86)
    if STRAIGHTEN:
        print(f"  STRAIGHTEN = True ({STRAIGHTEN_FIT} fit on "
              f"{STRAIGHTEN_REF} start)")
        print(f"    Each profile is sheared so the shoreline is horizontal")
        print(f"    BEFORE the window is picked. Windows should come out narrow")
        print(f"    (~4 cells) instead of spanning the diagonal (~12).")
        print(f"    Picks are frame-dependent: this WINDOW_JSON must have been")
        print(f"    picked with STRAIGHTEN=True or the run pass will refuse.")
        print(f"      {WINDOW_JSON}")
    else:
        print(f"  STRAIGHTEN = False -- the shoreline crosses each domain")
        print(f"    diagonally, so windows must span it and dune picks get")
        print(f"    loose in the high-obliquity domains.")
    if USE_CONST_INTERIOR:
        print(f"  USE_CONST_INTERIOR = True -- the interior is cut horizontally")
        print(f"    at max(dune_loc)+1, which throws away everything seaward of")
        print(f"    that on every other profile. Straightening shrinks the loss")
        print(f"    but does not remove it; residual dune variability is real.")
        print(f"    Consider False.")
    else:
        print(f"  USE_CONST_INTERIOR = False -- interior sheared per profile,")
        print(f"    no wedge lost. InteriorDomain[0] is each profile's dune.")
    if SHOW_ROAD:
        print(f"  SHOW_ROAD = True -- NC-12 "
              f"({', '.join(str(y) for y in ROAD_YEARS)}) drawn on the picker and")
        print(f"    every per-domain figure, in the SAME sheared/trimmed frame as")
        print(f"    the topography. Display and diagnostics only: the arrays")
        print(f"    CASCADE reads are byte-identical with SHOW_ROAD off.")
        print(f"    A search window sitting landward of the road is picking the")
        print(f"    road embankment, not a dune crest.")
    if MODE in ("pick", "pick_and_run"):
        print(f"  PICKS -> {WINDOW_JSON.name}")
        print(f"    rewritten after every domain. REPICK_EXISTING="
              f"{REPICK_EXISTING}.")
    print()
    if not check_paths():
        return
    load_dir = Path(LOAD_PATH)
    topo_dir = Path(TOPO_SAVE_PATH)
    dune_dir = Path(DUNE_SAVE_PATH)
    topo_dir.mkdir(parents=True, exist_ok=True)
    dune_dir.mkdir(parents=True, exist_ok=True)

    names = sorted(
        [n for n in os.listdir(load_dir)
         if n.endswith(".npy") and n.startswith("domain_")],
        key=natural_key,
    )
    if PICK_DOMAINS is not None:
        wanted = set(PICK_DOMAINS)
        names = [n for n in names if domain_number(Path(n).stem) in wanted]
    print(f"[info] Found {len(names)} domain file(s) in {load_dir}")

    windows = load_windows(Path(WINDOW_JSON))

    # ---------------- PICK PASS ----------------
    if MODE in ("pick", "pick_and_run"):
        for name in names:
            stem = Path(name).stem
            if not REPICK_EXISTING and stem in windows:
                w = windows[stem]
                if bool(w.get("straightened", False)) != bool(STRAIGHTEN):
                    print(f"[repick] {stem}: saved window is from the "
                          f"STRAIGHTEN={bool(w.get('straightened', False))} "
                          f"frame, this run is STRAIGHTEN={STRAIGHTEN} -- "
                          f"re-picking")
                else:
                    print(f"[skip pick] {stem}: window [{w['i0']}, {w['i1']}] "
                          f"already saved")
                    continue
            try:
                dom = load_profiles(load_dir / name)
            except Exception as e:
                print(f"[skip] {name}: {e}")
                continue

            # Open on the SAVED window when there is one, so a re-pick is an
            # adjustment: the v4 file was seeded from v3, so each domain shows
            # its v3 window, "r" resets to it, and accepting without dragging
            # keeps it. Falling back to default_window here -- as this did before
            # v4 -- would have made every one of the 90 re-picks a blind redraw.
            saved = windows.get(stem)
            if saved and bool(saved.get("straightened", False)) == bool(STRAIGHTEN):
                init = (int(saved["i0"]), int(saved["i1"]))
            else:
                init = default_window(dom["z"], dom["start_beach"])
            action, i0, i1 = pick_window(stem, dom["z"], dom["start_beach"], init,
                                         road_masks=dom.get("road_masks"))

            if action == "quit":
                print("[info] quit picking; saved windows kept")
                break
            if action == "skip":
                print(f"[pick] {stem}: skipped -> default window {init}")
                continue

            windows[stem] = {
                "i0": int(i0), "i1": int(i1),
                "n_cross_trimmed": int(dom["z"].shape[1]),
                "n_along": int(dom["z"].shape[0]),
                "trim_offset_c0": int(dom["c0"]),
                # the frame this window was picked in. A window picked
                # straightened is a valid index range on an unstraightened
                # array; it just points at different cells. The run pass
                # refuses on a mismatch rather than quietly using it.
                "straightened": bool(STRAIGHTEN),
                "obliquity_deg": dom["obliquity_deg"],
                "shear_max_cells": int(np.max(dom["shear"])),
                "picked": datetime.now().isoformat(timespec="seconds"),
                # What this window replaced, so the v3 -> v4 change is auditable
                # from the picks file alone. Equal values mean the road showed
                # nothing wrong with the old window and it was accepted as-is,
                # which is a result worth being able to see.
                "prev_i0": init[0], "prev_i1": init[1],
                "changed": bool((int(i0), int(i1)) != (init[0], init[1])),
            }
            save_windows(Path(WINDOW_JSON), windows)  # save every domain
            moved = "" if (int(i0), int(i1)) == (init[0], init[1]) else \
                f"  (was [{init[0]}, {init[1]}])"
            print(f"[pick] {stem}: window [{i0}, {i1}] saved{moved}")

    # ---------------- RUN PASS ----------------
    if MODE in ("run", "pick_and_run"):
        summary, rows = [], []
        for name in names:
            stem = Path(name).stem
            try:
                dom = load_profiles(load_dir / name)
            except Exception as e:
                print(f"[skip] {name}: {e}")
                continue
            prof_arr, start_beach = dom["z"], dom["start_beach"]

            w = windows.get(stem)
            if w is None:
                i0, i1 = default_window(prof_arr, start_beach)
                print(f"[warn] {stem}: no picked window, using default [{i0}, {i1}]")
            else:
                # absent == picked before straightening existed == False.
                # NOT "unknown, proceed": that would silently apply an
                # unstraightened window to a straightened array, which is a
                # valid index range pointing at the wrong cells.
                w_str = bool(w.get("straightened", False))
                if w_str != bool(STRAIGHTEN):
                    print(f"[skip] {stem}: window was picked with "
                          f"STRAIGHTEN={w_str} but this run has "
                          f"STRAIGHTEN={STRAIGHTEN}. The frames differ by the "
                          f"shear, so the window points at different cells. "
                          f"Re-pick, or point WINDOW_JSON at the matching set.")
                    continue
                if w.get("n_cross_trimmed") != prof_arr.shape[1]:
                    print(f"[warn] {stem}: trimmed width changed "
                          f"({w.get('n_cross_trimmed')} -> {prof_arr.shape[1]}); "
                          f"the saved window may no longer line up. Re-pick.")
                i0, i1 = int(w["i0"]), int(w["i1"])

            res = extract_domain(stem, prof_arr, start_beach, i0, i1,
                                 topo_dir, dune_dir,
                                 shear=dom["shear"],
                                 obliquity_deg=dom["obliquity_deg"],
                                 road_masks=dom.get("road_masks"))
            if res is None:
                continue
            res["c0"] = dom["c0"]
            # carried for island_plan_figure, which needs the road in the saved
            # interior frame and so needs the masks and the row0 line together
            res["road_masks"] = dom.get("road_masks")
            if SAVE_QC_FIGS:
                qc_figure(stem, prof_arr, start_beach, res, Path(FIG_DIR_QC),
                          road_masks=dom.get("road_masks"))
            if SAVE_COMPARISON_FIGS:
                comparison_figure(stem, dom, res, Path(FIG_DIR_CMP))
            summary.append(res)
            rows.append(settings_row(res, w))

        if SAVE_SETTINGS_SHEET:
            write_settings_sheet(rows, Path(SHEET_SAVE_PATH))
        summary_figure(rows, Path(SUMMARY_FIG_PATH))
        if SAVE_ISLAND_FIG:
            _off = load_offsets()
            island_figure(summary, _off, Path(ISLAND_FIG_PATH))
            if SAVE_ISLAND_PLAN_FIG:
                island_plan_figure(summary, _off, Path(RUN_DIR))
        write_manifest(rows, Path(MANIFEST_PATH))

        print("\n" + "=" * 86)
        print(f"{'domain':<12}{'section':<26}{'window':>12}{'dunes':>8}"
              f"{'mean h (m)':>12}{'interior rows':>16}")
        print("-" * 86)
        for r in rows:
            win = "[{},{}]".format(r["dune window start"], r["dune window end"])
            print(f"{r['stem']:<12}{r['section'][:25]:<26}{win:>12}"
                  f"{r['dunes found']:>8}{r['mean dune h (m)']:>12.2f}"
                  f"{r['interior rows']:>16}")
        print("=" * 86)

        n_default = sum(1 for r in rows if r["window source"] == "default")
        n_filled = sum(1 for r in rows if r["dunes filled"] > 0)
        if n_default:
            print(f"[!] {n_default} domain(s) used the fallback default window "
                  f"— re-run MODE='pick' for those")
        if n_filled:
            print(f"[!] {n_filled} domain(s) had profiles with no dune found")

        if SHOW_ROAD:
            for year in ROAD_YEARS:
                have = [r for r in rows if r.get(f"road profiles {year}", 0)]
                seaward = [r["stem"] for r in rows
                           if isinstance(r.get(f"road seaward profiles {year}"), int)
                           and r[f"road seaward profiles {year}"] > 0]
                print(f"[road {year}] NC-12 present in {len(have)}/{len(rows)} "
                      f"domain(s)")
                if seaward:
                    print(f"[road {year}] road SEAWARD of interior row 0 on some "
                          f"profiles in {len(seaward)} domain(s): "
                          f"{', '.join(seaward[:10])}"
                          f"{' ...' if len(seaward) > 10 else ''}")
                    print(f"            Barrier3D cannot place a road at a "
                          f"negative setback; see RoadOffset_dunestart_audit.md.")
        print(f"\n[done] {len(summary)} domain(s). Everything for this run is in:"
              f"\n       {RUN_DIR}\n"
              f"       start with RUN_MANIFEST.txt and "
              f"{SUMMARY_FIG_PATH.name}\n"
              f"       CASCADE inputs: {topo_dir}\n"
              f"                       {dune_dir}")


if __name__ == "__main__":
    main()
