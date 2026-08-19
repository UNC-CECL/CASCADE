#!/usr/bin/env python3
"""
HATTERAS ISLAND — CASCADE hindcast runner (1984-2004 / 2004-2024)

Modeled shoreline change rate:
  change_rate = (shoreline_m[-1, :] - shoreline_m[0, :]) / (END_YEAR - START_YEAR)  [m/yr]

CoastSat overlays: mean_lrr (m/yr) from domain_lrr_summary.csv, GIS 1-90 mapped to
padded indices 15-104 (see Section 1). Active period solid line, secondary dashed.

Overwash filter: developed communities only, set to 0.4; everything else (incl.
buffers) is 0 -- Buxton GIS 7-8 (pad 21-22), Avon GIS 21-31 (pad 35-45),
Salvo/Waves/Rodanthe GIS 68-83 (pad 82-97).

CASCADE IS RUN AS PUBLISHED
  This script does not patch cascade.roadway_manager or cascade.barrier3d. An
  earlier "Section 0" monkey-patched bulldoze() to guard an unchecked index and
  to label drowning messages by domain; both are gone. The guard's job is now
  done upstream of the model -- HAT_road_setback_audit.py refuses any setback
  that would overrun -- and domain attribution comes from CASCADE's own
  post-run state (see summarise_road_management). The one known consequence is
  recorded in the audit doc under "Unguarded index in bulldoze()".

HOW THIS FILE IS ORGANIZED
  SECTION 1   Domain configuration (fixed geometry, not period-dependent)
  SECTION 2   Display option
  SECTION 3   File paths (static directories + filenames)
  SECTION 4   PERIOD CONFIGURATION <-- everything that changes between 1984-2004
              and 2004-2024: SLR rate, storm file, dune offsets, road setback
              file, nourishment flags, background-erosion rates. Flip START_YEAR
              below and every period-dependent value resolves automatically.
  SECTION 5   CoastSat datasets + LOESS smoothing config
  SECTION 6   Geographic annotation styling (publication figure)
  SECTION 7   Simulation parameters that do NOT vary by period
  SECTION 8   Historical management events (road relocations, bridges,
              nourishment) - calendar-year-keyed, period-filtered automatically
  SECTION 9   Load input data (files resolved in Section 4/8)
  SECTION 10  Elevation + dune file lists
  HELPER FUNCTIONS / CASCADE RUNNER / MAIN
"""

import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from matplotlib.transforms import blended_transform_factory
from statsmodels.nonparametric.smoothers_lowess import lowess
# While TESTING the groin, use the SANDBOX copy of cascade.py (cascade_groin.py)
# with the 3-line pre-AST hook added. Real cascade.py stays untouched until the
# groin is proven; then fold the hook in and set USE_SANDBOX_CASCADE=False.
USE_SANDBOX_CASCADE = True

if USE_SANDBOX_CASCADE:
    from cascade.cascade_groin import Cascade   # hooked sandbox copy
else:
    from cascade.cascade import Cascade         # real package (hook folded in)
import cascade.roadway_manager as _rm

# =============================================================================
# SECTION 1: DOMAIN CONFIGURATION
# =============================================================================
# Fixed geometry of the model grid. Nothing here depends on START_YEAR.

NUM_REAL_DOMAINS   = 90
NUM_BUFFER_DOMAINS = 15

FIRST_FILE_NUMBER = 1    # GIS domain IDs: 1 .. 90
LAST_FILE_NUMBER  = FIRST_FILE_NUMBER + NUM_REAL_DOMAINS - 1  # = 90

TOTAL_DOMAINS    = NUM_BUFFER_DOMAINS + NUM_REAL_DOMAINS + NUM_BUFFER_DOMAINS  # 120
START_REAL_INDEX = NUM_BUFFER_DOMAINS          # = 15  (first real domain in padded array)
END_REAL_INDEX   = START_REAL_INDEX + NUM_REAL_DOMAINS  # = 105

# Roads occupy GIS domains 9-90  ->  padded indices 23-104
FIRST_ROAD_DOMAIN = 9
LAST_ROAD_DOMAIN  = 90
START_ROAD_INDEX  = (FIRST_ROAD_DOMAIN - 1) + NUM_BUFFER_DOMAINS   # = 23
END_ROAD_INDEX    = (LAST_ROAD_DOMAIN  - 1) + NUM_BUFFER_DOMAINS + 1  # = 105


def _gis_to_pad(gis_id):
    """Convert a 1-based GIS domain ID to a CASCADE padded array index."""
    return START_REAL_INDEX + (gis_id - FIRST_FILE_NUMBER)


# =============================================================================
# SECTION 1B: GROIN CONFIGURATION (Buxton Groin Field)
# =============================================================================
# Fixed real-world facts about the structure -- NOT period-dependent, so this
# lives here rather than in Section 4. install_year is the true historical
# construction date; deterioration is specified as a delay relative to it, so
# the same config resolves correctly whichever period (START_YEAR) is active:
#   Period 1 (1984-2004): install year predates the run, so the groin is
#     already active at year 0; the 1996->2003 deterioration plays out
#     mid-run (years 12-19) -- this is the only run that exercises the ramp.
#   Period 2 (2004-2024): deterioration completed before this period starts,
#     so the groin is already at its floor fraction from year 0 onward.
#
# GROIN_ENABLED is the master switch -- flip it and rerun, same pattern as
# flipping START_YEAR, rather than a run-matrix loop.
GROIN_ENABLED = False   # TEMP: set by end-to-end test run; restore to True

GROIN_UPDRIFT_GIS   = 6     # source: updrift accretes (dipole +)
GROIN_DOWNDRIFT_GIS = 5     # sink:   downdrift erodes (dipole -)
GROIN_INSTALL_YEAR  = 1969  # confirmed historical construction date

GROIN_TRAPPING_RATE_M_YR = 60.0   # M -- placeholder from the 1967-1997 test; retune against Period 1 CoastSat LRR

# Deterioration: last repair 1996 -> major storm damage 2003. Modeled as a
# linear ramp bridging the two documented dates (maintenance stops in 1996,
# unrepaired damage accumulates, the 2003 storm locks in the new steady
# state) rather than treating either date alone as the trigger.
GROIN_DETERIORATION_DELAY_YEARS = 1996 - GROIN_INSTALL_YEAR   # = 27
GROIN_DETERIORATION_MODE        = "linear_ramp"
GROIN_DETERIORATION_RAMP_YEARS  = 2003 - 1996                  # = 7
GROIN_DETERIORATION_FRACTION    = 0.2   # placeholder -- solve from era LRR ratio at D5/D6

# Module import path -- update if HAT_groin_module.py has moved since the
# 1967-1997 test.
GROIN_MODULE_IMPORT_PATH = "scripts.groin.HAT-hindcast-groin-test.version_control.HAT_groin_module"


# =============================================================================
# SECTION 2: DISPLAY OPTION
# =============================================================================

# Set to True  -> x-axis shows only the 90 real domains (no buffers)
# Set to False -> x-axis shows all 120 padded domains (buffers shaded)
PLOT_REAL_DOMAINS_ONLY = True

# Pops up the annotated publication figure (*_annotated.png) at run end.
# It's always saved to disk either way - this only controls the on-screen
# plt.show() popup, so turning it off just quiets the run without losing
# the file. Kept off by default so a run only pops up the main rate plot.
SHOW_ANNOTATED_FIGURE = False

# -----------------------------------------------------------------------------
# SECTION 2a: SHORELINE ANIMATION (GIF)
# -----------------------------------------------------------------------------
# Plan-view animation of the shoreline, one frame/model year, saved to the
# run's comparison folder. Blue = seaward of year-0 shoreline (accretion), red =
# landward (erosion). Each GIF_JOBS entry produces one GIF.
#
# "range":
#   "real"       90 real GIS domains (1-90), full geographic annotation
#   "all"        all 120 padded domains (buffers shaded, GIS ids on top axis)
#   "groin"      auto-window around a groin (+/- job's "pad"), read from
#                ANN_GROINS (Section 6). Fans out over every ANN_GROINS entry,
#                so adding a structure there needs no GIF_JOBS change. Restrict
#                with "which": a name, list of names, or "all" (default).
#   "groin_span" one window enclosing every ANN_GROINS groin (+/- pad); use
#                once structures are close enough to want a single frame.
#   (lo, hi)     explicit inclusive GIS domain window, e.g. (1, 30)
#
# "mode" (y-axis):
#   "position"     cross-shore position relative to year-0 alongshore mean;
#                  preserves planform shape, best at narrow ranges.
#   "displacement" cumulative change from each domain's own year-0 position
#                  (every domain starts at 0); preferred for "real"/"all".
#   "difference"   this run minus a baseline run, in displacement terms:
#                  (run - run_yr0) - (baseline - baseline_yr0). Blue = this
#                  run seaward of baseline. Needs GIF_BASELINE_NPY, else skipped.
#
# "pad" only read by "groin"/"groin_span"; "which" only by "groin".
MAKE_SHORELINE_GIF = True
GIF_JOBS = [
    dict(range="real",  mode="displacement"),                       # whole island overview
    dict(range="real",  mode="position", auto_open=True),           # whole island, planform position - shown at run end
    dict(range="groin", mode="position",   pad=9),      # zoom, one per groin
    dict(range="groin", mode="difference", pad=9),      # groin effect vs baseline
    # dict(range="groin_span", mode="difference", pad=9),   # all groins in one frame
    # dict(range="groin", mode="difference", pad=9, which="Buxton Groin"),  # just one
]

# Previous run's *_shoreline_matrix.npy to difference against (set to the
# no-groin run when testing the groin module). None skips "difference" jobs.
# Must match domain count + run length of the current run, e.g.:
#   os.path.join(OUTPUT_BASE_DIR, "HAT_1984_2004_noGroin_Hs1.0",
#                "HAT_1984_2004_noGroin_Hs1.0_shoreline_matrix.npy")
GIF_BASELINE_NPY   = None
GIF_BASELINE_LABEL = "no-groin baseline"   # shown in the difference-GIF title

GIF_FPS            = 3               # frames per second
GIF_YEAR_STRIDE    = 1               # 1 = every model year; 2 = every other, etc.
GIF_ANNOTATE       = True            # geographic annotation layer (GIS-axis modes)
GIF_AUTO_OPEN      = False           # pop the GIF open in the default viewer when done
GIF_KEEP_FRAMES    = False           # also write the individual frame PNGs to disk
GIF_SAVE_MATRIX    = True            # save shoreline_m.npy (needed as a future baseline)
GIF_OCEAN_AT_BOTTOM = True           # y-axis orientation: True -> ocean/seaward at the
                                      # bottom of the plot, sound/landward at the top
                                      # (matches Hatteras' real cross-shore layout).
                                      # False -> original convention, seaward at the top.

# =============================================================================
# SECTION 3: FILE PATHS
# =============================================================================

# Derived from this file's location (scripts/hatteras_ms/) so the script runs on
# any checkout without editing a hardcoded path.
PROJECT_BASE_DIR   = str(Path(__file__).resolve().parents[2])
HATTERAS_DATA_BASE = os.path.join(PROJECT_BASE_DIR, "data", "hatteras_init")
print(f"[DEBUG] PROJECT_BASE_DIR   = {PROJECT_BASE_DIR!r}")
print(f"[DEBUG] HATTERAS_DATA_BASE = {HATTERAS_DATA_BASE!r}")
OUTPUT_BASE_DIR    = os.path.join(PROJECT_BASE_DIR, "comparison", "raw_runs")
COASTSAT_BASE_DIR  = os.path.join(
    PROJECT_BASE_DIR, "scripts", "input_prep", "CoastSat"
)

# Filename only - CASCADE resolves it relative to datadir (hatteras_init)
PARAMETER_FILE = "Hatteras-CASCADE-parameters.yaml"

# --- Topo/dune init year ---
# TOPO_DUNE_INIT_YEAR : year label used in the *filename* (e.g. domain_7_topography_2009.npy)
# TOPO_DUNE_VERSION   : version folder holding the dunes/ and topography/ subdirs.
#                       Change it to switch between extractor versions (2009_v2,
#                       etc.) without touching the filename year.
# NOTE: this is the same initial condition for BOTH periods (it is the surface
# the model starts from at whichever START_YEAR is selected) - it is not part
# of the period/calibration switch in Section 4.
TOPO_DUNE_INIT_YEAR = "2009"      # used in filenames - matches extractor comparison
TOPO_DUNE_VERSION   = "2009_v2"   # version folder under 1-barrier3d-domains/2009-dune-topo/

# Layout: 1-barrier3d-domains/2009-dune-topo/<version>/{dunes,topography}/domain_<n>_*_<year>.npy
BARRIER3D_DIR = os.path.join(HATTERAS_DATA_BASE, "1-barrier3d-domains")
DUNE_TOPO_DIR = os.path.join(BARRIER3D_DIR, "2009-dune-topo", TOPO_DUNE_VERSION)
BUFFER_DIR    = os.path.join(BARRIER3D_DIR, "2009-buffer")

os.chdir(PROJECT_BASE_DIR)
os.makedirs(OUTPUT_BASE_DIR, exist_ok=True)

# =============================================================================
# SECTION 4: PERIOD CONFIGURATION  <-- everything that changes with START_YEAR
# =============================================================================
# Flip START_YEAR below (1984 or 2004) and every period-dependent value
# resolves automatically: end year, SLR rate, storm file, island offset file,
# road setback file, nourishment on/off + volume default, and the
# background-erosion rates loaded into DOMAIN_BE_RATES.
# First place to check "did I get 1984 vs 2004 right?"
#
# To enable nourishment for 2004: set enable_nourishment=True and update
# nourishment_volume (m^3/m) in the 2004 entry below (already on).

RUN_NAME_SUFFIX = "full_calibrated"   # <- edit this to label each experiment

# <- Flip between 1984 and 2004 to switch periods; END_YEAR is set automatically
START_YEAR = 1984
# ── SOURCE_SINK_PRESET ───────────────────────────────────────────────────────
# Which background-erosion value set is used, independent of period:
#   "base"       only edge boundary corrections (GIS 1 and 90); interior = 0.0
#   "calibrated" full corrections
# Both presets are period-keyed below (Section 4a), so the correct GIS-1/GIS-90
# edge values (and, for "calibrated", the full interior fit) load automatically.
SOURCE_SINK_PRESET = "calibrated"

# -----------------------------------------------------------------------------
# SECTION 4a: BACKGROUND EROSION (SOURCE/SINK) RATES — keyed by period + preset
# -----------------------------------------------------------------------------
# Units: m/yr

DOMAIN_BE_RATES_BASE_BY_PERIOD = {
    1984: {
        1:  -40, #-40
        90:  15, # 15
    },
    2004: {
        1:   35, #35
        90:  35, #35
    },
}

# "Calibrated" preset
DOMAIN_BE_RATES_CALIBRATED_BY_PERIOD = {
    1984: {
          1: -40.0,  # LOCKED — use your solved value, not 0.0
          2: -0.6,  # Cape Point / Shoal Dynamics
          3: -0.7,  # Cape Point / Shoal Dynamics
          4: 0.0,
          5: +1.2,  # Cape Point / Shoal Dynamics
          6: +2.0,  # Cape Point / Shoal Dynamics
          7: +0.9,  # Cape Point / Shoal Dynamics
          8: +1.1,  # Cape Point / Shoal Dynamics
          9: 0.0,
         10: -0.9,  # Cape Point / Shoal Dynamics
         11: -1.3,  # Buxton–Avon Transition
         12: -0.9,  # Buxton–Avon Transition
         13: -0.6,  # Buxton–Avon Transition
         14: 0.0,
         15: 0.0,
         16: 0.0,
         17: +0.6,  # Buxton–Avon Transition
         18: +0.6,  # Buxton–Avon Transition
         19: +0.6,  # Buxton–Avon Transition
         20: 0.0,
         21: 0.0,
         22: 0.0,
         23: 0.0,
         24: +0.3,  # Avon
         25: +0.3,  # Avon
         26: +0.8,  # Avon
         27: +1.4,  # Avon
         28: +1.9,  # Avon
         29: +2.3,  # Avon
         30: +2.2,  # Avon
         31: +2.2,  # Avon
         32: +1.9,  # Mid-island
         33: +1.4,  # Mid-island
         34: +1.0,  # Mid-island
         35: 0.0,
         36: 0.0,
         37: 0.0,
         38: 0.0,
         39: 0.0,
         40: 0.0,
         41: 0.0,
         42: 0.0,
         43: 0.0,
         44: +0.4,  # Mid-island
         45: 0.0,
         46: 0.0,
         47: 0.0,
         48: -0.3,  # Mid-island
         49: -1.0,  # Mid-island
         50: -1.0,  # Mid-island
         51: -1.2,  # Mid-island
         52: -1.3,  # Mid-island
         53: -1.3,  # Mid-island
         54: -1.2,  # Mid-island
         55: -1.0,  # Mid-island
         56: -1.0,  # Mid-island
         57: -0.3,  # Mid-island
         58: 0.0,
         59: 0.0,
         60: 0.0,
         61: 0.0,
         62: +0.4,  # Wimble Shoals Influence
         63: 0.0,
         64: 0.0,
         65: 0.0,
         66: 0.0,
         67: 0.0,
         68: +0.7,  # Wimble Shoals Influence
         69: +1.2,  # Wimble Shoals Influence
         70: +1.6,  # Wimble Shoals Influence
         71: +1.9,  # Wimble Shoals Influence
         72: +2.2,  # Wimble Shoals Influence
         73: +2.4,  # Wimble Shoals Influence
         74: +2.6,  # Wimble Shoals Influence
         75: +2.5,  # Tri-Village / Rodanthe
         76: +1.7,  # Tri-Village / Rodanthe
         77: +1.3,  # Tri-Village / Rodanthe
         78: +0.7,  # Tri-Village / Rodanthe
         79: 0.0,
         80: -1.3,  # Tri-Village / Rodanthe
         81: -1.9,  # Tri-Village / Rodanthe
         82: -2.1,  # Tri-Village / Rodanthe
         83: -2.3,  # Tri-Village / Rodanthe
         84: -2.5,  # Pea Island NWR
         85: -2.5,  # Pea Island NWR
         86: -2.4,  # Pea Island NWR
         87: -2.0,  # Pea Island NWR
         88: -1.6,  # Pea Island NWR
         89: -1.0,  # Pea Island NWR
         90: 15.0,  # LOCKED — use your solved value, not 0.0
    },

    2004: {
          1: 35.0,  # LOCKED — use your solved value, not 0.0
          2: +1.0,  # Cape Point / Shoal Dynamics
          3: +1.6,  # Cape Point / Shoal Dynamics
          4: +1.6,  # Cape Point / Shoal Dynamics
          5: +1.2,  # Cape Point / Shoal Dynamics
          6: -1.5,  # Cape Point / Shoal Dynamics
          7: 0.0,
          8: +1.1,  # Cape Point / Shoal Dynamics
          9: +2.2,  # Cape Point / Shoal Dynamics
         10: +3.1,  # Cape Point / Shoal Dynamics
         11: +3.4,  # Buxton–Avon Transition
         12: +3.4,  # Buxton–Avon Transition
         13: +3.3,  # Buxton–Avon Transition
         14: +3.2,  # Buxton–Avon Transition
         15: +3.0,  # Buxton–Avon Transition
         16: +3.1,  # Buxton–Avon Transition
         17: +3.0,  # Buxton–Avon Transition
         18: +2.9,  # Buxton–Avon Transition
         19: +2.6,  # Buxton–Avon Transition
         20: +2.2,  # Buxton–Avon Transition
         21: +1.9,  # Avon
         22: +1.4,  # Avon
         23: +0.8,  # Avon
         24: +0.3,  # Avon
         25: +0.3,  # Avon
         26: +0.8,  # Avon
         27: +1.4,  # Avon
         28: +1.9,  # Avon
         29: +2.3,  # Avon
         30: +3.1,  # Avon
         31: +3.6,  # Avon
         32: +3.8,  # Mid-island
         33: +3.7,  # Mid-island
         34: +3.6,  # Mid-island
         35: +3.5,  # Mid-island
         36: +3.3,  # Mid-island
         37: +3.1,  # Mid-island
         38: +2.7,  # Mid-island
         39: +2.4,  # Mid-island
         40: +2.1,  # Mid-island
         41: +1.8,  # Mid-island
         42: +1.5,  # Mid-island
         43: +1.1,  # Mid-island
         44: +0.4,  # Mid-island
         45: 0.0,
         46: 0.0,
         47: 0.0,
         48: -0.3,  # Mid-island
         49: 0.0,
         50: -1.0,  # Mid-island
         51: -1.2,  # Mid-island
         52: -1.3,  # Mid-island
         53: -1.3,  # Mid-island
         54: -1.2,  # Mid-island
         55: -1.0,  # Mid-island
         56: 0.0,
         57: -0.3,  # Mid-island
         58: 0.0,
         59: 0.0,
         60: 0.0,
         61: 0.0,
         62: +0.4,  # Wimble Shoals Influence
         63: +1.1,  # Wimble Shoals Influence
         64: +1.6,  # Wimble Shoals Influence
         65: +2.0,  # Wimble Shoals Influence
         66: +2.4,  # Wimble Shoals Influence
         67: +2.8,  # Wimble Shoals Influence
         68: +3.3,  # Wimble Shoals Influence
         69: +3.7,  # Wimble Shoals Influence
         70: +4.0,  # Wimble Shoals Influence
         71: +4.2,  # Wimble Shoals Influence
         72: +4.1,  # Wimble Shoals Influence
         73: +4.0,  # Wimble Shoals Influence
         74: +3.7,  # Wimble Shoals Influence
         75: +3.4,  # Tri-Village / Rodanthe
         76: +3.0,  # Tri-Village / Rodanthe
         77: +2.7,  # Tri-Village / Rodanthe
         78: +2.4,  # Tri-Village / Rodanthe
         79: +2.3,  # Tri-Village / Rodanthe
         80: +2.0,  # Tri-Village / Rodanthe
         81: +1.7,  # Tri-Village / Rodanthe
         82: +1.4,  # Tri-Village / Rodanthe
         83: +1.0,  # Tri-Village / Rodanthe
         84: +0.7,  # Pea Island NWR
         85: 0.0,
         86: +0.5,  # Pea Island NWR
         87: +0.6,  # Pea Island NWR
         88: +0.7,  # Pea Island NWR
         89: +1.2,  # Pea Island NWR
         90: 35.0,  # LOCKED — use your solved value, not 0.0
    },  #
}


if DOMAIN_BE_RATES_CALIBRATED_BY_PERIOD.get(2004) is None:
    DOMAIN_BE_RATES_CALIBRATED_BY_PERIOD[2004] = dict(DOMAIN_BE_RATES_CALIBRATED_BY_PERIOD[1984])
    DOMAIN_BE_RATES_CALIBRATED_2004_IS_PLACEHOLDER = True
else:
    DOMAIN_BE_RATES_CALIBRATED_2004_IS_PLACEHOLDER = False

# -----------------------------------------------------------------------------
# SECTION 4b: PERIOD CONFIG — file paths, SLR rate, nourishment, BE preset
# -----------------------------------------------------------------------------

PERIOD_CONFIG = {
    1984: dict(
        end_year             = 2004,
        sea_level_rise_rate  = 0.004,   # m/yr - from duck_rslr_analysis.py
        storm_file           = os.path.join(
            HATTERAS_DATA_BASE, "3-env-forcings", "storms", "hindcast_storms", "1984_2004",
            "1984_2004_storms_v3_72.npy",
        ),
        island_offset_file     = os.path.join(
            HATTERAS_DATA_BASE, "2-brie-offset", "hindcast_1984",
            f"Island_Dune_Offsets_1984_PADDED_{TOTAL_DOMAINS}.csv",
        ),
        road_setback_file    = os.path.join(
            HATTERAS_DATA_BASE, "4-mgmt-forcing", "road_offset", "old_method_offset", "1984", "RoadSetback_1984.csv",
        ),
        enable_nourishment   = False,
        nourishment_volume   = 0,
        be_rates_base        = DOMAIN_BE_RATES_BASE_BY_PERIOD[1984],
        be_rates_calibrated  = DOMAIN_BE_RATES_CALIBRATED_BY_PERIOD[1984],
    ),
    2004: dict(
        end_year             = 2024,
        sea_level_rise_rate  = 0.006,   # m/yr - from duck_rslr_analysis.py
        storm_file           = os.path.join(
            HATTERAS_DATA_BASE, "3-env-forcings", "storms", "hindcast_storms", "2004_2024",
            "2004_2024_storms_v3_72.npy",
        ),
        island_offset_file     = os.path.join(
            HATTERAS_DATA_BASE, "2-brie-offset", "hindcast_2004",
            f"Island_Dune_Offsets_2004_PADDED_{TOTAL_DOMAINS}.csv",
        ),
        road_setback_file    = os.path.join(
            HATTERAS_DATA_BASE, "4-mgmt-forcing", "road_offset", "old_method_offset", "2004", "RoadSetback_2004.csv",
        ),  # UPDATE path once file is generated, if not already final
        enable_nourishment   = True,    # historical BN injected per-year in time loop
        nourishment_volume   = 100,     # m^3/m default passed to Cascade init (threshold unused)
        be_rates_base        = DOMAIN_BE_RATES_BASE_BY_PERIOD[2004],
        be_rates_calibrated  = DOMAIN_BE_RATES_CALIBRATED_BY_PERIOD[2004],
    ),
}

if START_YEAR not in PERIOD_CONFIG:
    print(f"ERROR: Invalid START_YEAR {START_YEAR}. Must be one of {list(PERIOD_CONFIG)}.")
    sys.exit(1)

_cfg = PERIOD_CONFIG[START_YEAR]

END_YEAR             = _cfg["end_year"]
RUN_YEARS            = END_YEAR - START_YEAR    # passed to CASCADE as time_step_count
SEA_LEVEL_RISE_RATE  = _cfg["sea_level_rise_rate"]
STORM_FILE           = _cfg["storm_file"]
ISLAND_OFFSET_FILE     = _cfg["island_offset_file"]
ROAD_SETBACK_FILE    = _cfg["road_setback_file"]
ENABLE_NOURISHMENT   = _cfg["enable_nourishment"]
NOURISHMENT_VOLUME   = _cfg["nourishment_volume"]
RUN_NAME_BASE        = f"HAT_{START_YEAR}_{END_YEAR}_{RUN_NAME_SUFFIX}"

# ── Resolve active BE preset for the active period ──────────────────────────
_BE_PRESETS_THIS_PERIOD = {
    "base":       _cfg["be_rates_base"],
    "calibrated": _cfg["be_rates_calibrated"],
}
if SOURCE_SINK_PRESET not in _BE_PRESETS_THIS_PERIOD:
    print(f"ERROR: SOURCE_SINK_PRESET '{SOURCE_SINK_PRESET}' is not valid. "
          f"Must be one of {list(_BE_PRESETS_THIS_PERIOD)}.")
    sys.exit(1)
DOMAIN_BE_RATES = _BE_PRESETS_THIS_PERIOD[SOURCE_SINK_PRESET]
if DOMAIN_BE_RATES is None:
    print(f"ERROR: No '{SOURCE_SINK_PRESET}' background-erosion rates defined "
          f"for START_YEAR={START_YEAR}. See Section 4a.")
    sys.exit(1)

# Loud, impossible-to-miss warning: catches the exact mistake that caused
# "editing the 2004 calibrated values doesn't appear to do anything" - the
# 2004 calibrated dict is independent now (see Section 4a fix), but it may
# still BE an unverified copy of the 1984 fit rather than a real Period-2
# calibration, depending on whether DOMAIN_BE_RATES_CALIBRATED_2004_IS_PLACEHOLDER
# has been set to False after pasting in real solved values.
if (START_YEAR == 2004 and SOURCE_SINK_PRESET == "calibrated"
        and DOMAIN_BE_RATES_CALIBRATED_2004_IS_PLACEHOLDER):
    print("=" * 80)
    print("WARNING: START_YEAR=2004, SOURCE_SINK_PRESET='calibrated' is still")
    print("         using the 1984 PLACEHOLDER interior shape (see Section 4a).")
    print("         Edits to this dict are independent of 1984 now, but no real")
    print("         Period-2 differential_evolution fit has been confirmed yet.")
    print("         Set DOMAIN_BE_RATES_CALIBRATED_2004_IS_PLACEHOLDER = False")
    print("         once you've pasted in real solved Period-2 values.")
    print("=" * 80)

SEA_LEVEL_CONSTANT = True
TO_METERS = True

print("=" * 80)
print("HATTERAS ISLAND CASCADE - PERIOD CONFIGURATION")
print("=" * 80)
print(f"Period:               {START_YEAR}-{END_YEAR}  ({RUN_YEARS} years)")
print(f"SLR rate:             {SEA_LEVEL_RISE_RATE * 1000:.1f} mm/yr")
print(f"SOURCE_SINK_PRESET:   '{SOURCE_SINK_PRESET}'  "
      f"({len(DOMAIN_BE_RATES)} non-default domain(s))")
print(f"Storm file:           {os.path.basename(STORM_FILE)}")
print(f"Island offset file:  {os.path.basename(ISLAND_OFFSET_FILE)}")
print(f"Road setback file:    {os.path.basename(ROAD_SETBACK_FILE)}")
print(f"Nourishment enabled:  {ENABLE_NOURISHMENT}  (default volume: {NOURISHMENT_VOLUME} m^3/m)")
print("=" * 80 + "\n")

# =============================================================================
# SECTION 5: COASTSAT DATASETS (transect-level) + LOESS SMOOTHING CONFIG
# =============================================================================
# Each entry points to transect_lrr_full.csv - one row per CoastSat transect
# (~50 m spacing, hundreds of transects total). LOESS is applied at transect
# resolution then aggregated to domain resolution for CASCADE comparison.
#
# period_start controls active vs reference styling:
#   active period (matches START_YEAR) -> full opacity scatter + solid LOESS lines
#   reference period                   -> faded scatter + faded LOESS lines
#
# NOTE: both CoastSat periods are always loaded (for the active/reference
# overlay distinction above) - this is observational data, not simulation
# config, so it intentionally is NOT folded into PERIOD_CONFIG in Section 4.

COASTSAT_DATASETS = [
    dict(
        label           = "CoastSat LRR (1984-2004)",
        period_start    = 1984,
        csv             = os.path.join(COASTSAT_BASE_DIR, "1984_2004", "transect_lrr_full.csv"),
        domain_col      = "domain_number",
        rate_col        = "lrr_m_yr",
        transect_id_col = "transect_id",
    ),
    dict(
        label           = "CoastSat LRR (2004-2024)",
        period_start    = 2004,
        csv             = os.path.join(
            COASTSAT_BASE_DIR, "2004_2024", "transect_lrr_full.csv"
        ),
        domain_col      = "domain_number",
        rate_col        = "lrr_m_yr",
        transect_id_col = "transect_id",
    ),
]

# --- LOESS smoothing applied to CoastSat overlay ---
# List one or two window sizes (domain units; 1 domain ~= 500 m).
# LOESS is applied at transect resolution (physical along-coast distance as x),
# then aggregated to domain resolution for the CASCADE comparison.
#   10 domains -> ~5.0 km  <- primary reference
#    7 domains -> ~3.5 km  <- narrower comparison window
LOESS_WINDOW_DOMAINS = [7, 10]   # list of 1 or 2 window sizes

# Styling per LOESS window (matched by list position).
# Tuple: (linewidth, linestyle, alpha_factor_for_active_period)
LOESS_WINDOW_STYLES = [
    (1.8, "-", 1.00),   # narrower window (7-dom): solid - distinguished by color
    (2.0, "-", 1.00),   # wider window   (10-dom): solid, full opacity, primary reference
]

# Which LOESS window to use as the reference curve in residual calculations.
# Must be one of the values in LOESS_WINDOW_DOMAINS.
RESIDUALS_LOESS_WINDOW = 10

# Southernmost GIS domains (1 through this value) shown as raw per-domain mean
# LRR instead of LOESS-smoothed - Oregon Inlet boundary effects dominate this
# zone and smoothing can obscure the sharp gradient there.
#   - Widest LOESS window: raw means (d1-N) stitched onto the LOESS line
#     (d(N+1)-90) so the curve stays visually continuous.
#   - Narrower windows: the line simply starts at domain N+1.
# Set to 0 to use LOESS for all domains.
LOESS_SKIP_SOUTHERN_DOMAINS = 10

# =============================================================================
# COLOUR PALETTE REFERENCE
# =============================================================================
# CoastSat - cool blue family. Three layers, light -> dark = less -> more smoothed:
#   transect scatter  #9ECAE1   very light blue    individual ~50 m transect LRR (dots)
#    7-domain LOESS   #6BAED6   medium sky blue    transect-based LOESS, narrow window
#   10-domain LOESS   #08519C   deep ocean blue    transect-based LOESS, primary reference
#
# Model line  ->  ANN_MODEL_COLOR  (warm orange)
#
# Design rationale:
#   Cool blue = observations (CoastSat)   Warm orange = model comparison (CASCADE)
#   Blue + orange is the most colorblind-safe pairing (deuteranopia / protanopia)
#   Within the blue family, lighter -> darker encodes less -> more smoothed
#   Period distinction (1984-2004 vs 2004-2024) is carried by opacity, not color
# =============================================================================

# CoastSat LOESS colors keyed by window size (domain count).
CS_WINDOW_COLORS = {
     7: "#6BAED6",   # medium sky blue  - 7-domain LOESS
    10: "#08519C",   # deep ocean blue  - 10-domain LOESS
}
CS_WINDOW_COLOR_DEFAULT = "#4A7C8E"   # fallback for any unlisted window size

# Individual transect scatter - plotted at lowest zorder as context.
CS_RAW_COLOR            = "#5BA3C9"    # medium blue - darkened from original #9ECAE1 for legibility
PLOT_RAW_LRR            = True         # set False to hide transect scatter entirely
RAW_LRR_SOUTHERN_ONLY   = True         # True  -> scatter only for D1-LOESS_SKIP_SOUTHERN_DOMAINS
                                        #          (the zone where LOESS is suppressed)
                                        # False -> scatter for all domains D1-90
PLOT_REFERENCE_PERIOD   = False        # set True to also show the non-active CoastSat period (faded)
RAW_LRR_SCATTER_SIZE    = 6            # marker area in points^2 - bumped from 4 for visibility
RAW_LRR_SCATTER_ALPHA   = 0.60         # opacity for active period (up from 0.35); x0.35 for reference

# =============================================================================
# SECTION 6: GEOGRAPHIC ANNOTATION STYLING (annotated publication figure)
# =============================================================================

ANN_TOWN_SPANS = {
    "Buxton":      (7,  8),
    "Avon":        (21, 31),
    "Tri-Village": (68, 83),   # Salvo / Waves / Rodanthe
}
ANN_VILLAGE_LINES = {"Salvo": 69, "Waves": 74, "Rodanthe": 80}

ANN_PIER_LABEL_Y  = 0.76   # default rotated label y for any pier (0=bottom, 1=top axes fraction)
ANN_GROIN_LABEL_Y = 0.68   # rotated label y for groin lines
ANN_PIERS = {
    "Avon Pier":     (26, ANN_PIER_LABEL_Y),   # (domain, label_y) - adjust per pier
    "Rodanthe Pier": (79, ANN_PIER_LABEL_Y),
}
ANN_GROINS        = {"Buxton Groin": 5.5}      # boundary between domains 5 and 6
ANN_WIMBLE_SHOALS = (60, 74)
ANN_AVON_SHOALS   = (24, 39)   # Avon Shoals influence zone

# Accretion / Erosion side labels - set to None for auto-computed midpoint,
# or a 0-1 axes fraction to pin to a fixed position.
LABEL_ACCRETION_Y = None
LABEL_EROSION_Y   = None

ANN_C_TOWN_SPAN    = "#90AFC5"
ANN_C_WIMBLE       = "#E0A800"   # amber - both shoal zones share this color
ANN_C_AVON_SHOALS  = "#E0A800"   # same amber as Wimble Shoals (same feature type)
ANN_C_VILLAGE_LINE = "0.40"
ANN_C_PIER         = "#1565C0"
ANN_C_GROIN        = "#B71C1C"

ANN_MODEL_COLOR = "#FF8C00"   # warm orange - modeled shoreline change rate

# =============================================================================
# SECTION 7: SIMULATION PARAMETERS (shared across both periods)
# =============================================================================
# Anything here is intentionally the SAME regardless of START_YEAR. If a
# value needs to differ by period, it belongs in PERIOD_CONFIG (Section 4)
# instead, not here.

# --- Site-specific vertical datum parameters ---
# Override CASCADE/Barrier3D defaults (berm_elevation=1.9, MHW=0.46), which are
# set for a different site. Must match dune_topo_extractor and HAT_create_storms.
BERM_ELEVATION = 1.7    # m NAVD88 - Hatteras Island, NCDOT-derived via NC State
MHW_ELEVATION  = 0.36   # m NAVD88 - Duck NC gauge (NOAA 8651370)

NUM_CORES = 1   # >1 causes crashes; leave at 1 for Hatteras runs

# --- Wave climate ---
# WAVE_HEIGHTS_TO_TEST: use a single value for production runs.
# Multiple values trigger a calibration sweep (each gets its own run folder).
# Wave parameters are always recorded in the run metadata file.
WAVE_HEIGHTS_TO_TEST           = [2.5]   # m - Hs calibration value(s)
FIXED_WAVE_PERIOD              = 8       # s
FIXED_WAVE_ASYMMETRY           = 0.7
FIXED_WAVE_ANGLE_HIGH_FRACTION = 0.1

# --- Storm / morphodynamics ---
DUNE_REBUILD_HEIGHT    = 3.0    # m
REBUILD_ELEV_THRESHOLD = 0.01   # dam
OVERWASH_TO_DUNE       = 9
SANDBAG_ELEV           = 0

# --- Management flags (shared across both periods) ---
ENABLE_ROADWAY_MANAGEMENT = True
ENABLE_SANDBAG_PLACEMENT  = False
# Note: ENABLE_NOURISHMENT and NOURISHMENT_VOLUME ARE period-specific -> Section 4

# --- Road relocation domain control ---
# Villages have fixed infrastructure that prevents NC-12 from relocating
# landward; the manager still runs overwash removal + dune rebuild there, just
# not relocation. Relocation BLOCKED in: Buxton (GIS 7-8), Avon (GIS 21-31),
# Salvo/Waves/Rodanthe (GIS 68-83). Allowed everywhere else (GIS 9-20, 32-67, 84-90).
VILLAGE_GIS_RANGES_NO_RELOCATION = [
    (7,  8),   # Buxton
    (21, 31),  # Avon
    (68, 83),  # Salvo / Waves / Rodanthe (Tri-Village)
]

DOMAIN_TICK_STEP    = 5
DOMAIN_SPACING_M    = 500   # metres per CASCADE domain (used to convert window_domains -> km)
DOMAIN_LENGTH_M     = 500   # same value - explicit alias used in nourishment volume conversions

FLIP_SIGN_MODEL = True  # flips only the modeled sign (no alongshore reversal)

# -----------------------------------------------------------------------------
# OVERWASH FILTER - per community
# 0 = no filtering (undeveloped / buffer domains)
# 0.4 is based on Laura's article
# -----------------------------------------------------------------------------
OVERWASH_FILTER_DEFAULT         = 0.0   # undeveloped domains + all buffer domains
OVERWASH_FILTER_BUXTON          = 0.4   # GIS  7- 8  (south end of Buxton)
OVERWASH_FILTER_AVON            = 0.4   # GIS 21-31
OVERWASH_FILTER_SALVO_WAVES_ROD = 0.4   # GIS 68-83  (Salvo / Waves / Rodanthe)

# ── Build the full 120-domain background-erosion array from DOMAIN_BE_RATES ──
# (DOMAIN_BE_RATES itself was resolved per-period in Section 4.)
BACKGROUND_EROSION_RATES = [0.0] * TOTAL_DOMAINS
for _gis_id, _rate in DOMAIN_BE_RATES.items():
    _pad_idx = START_REAL_INDEX + (_gis_id - 1)   # GIS 1 -> pad 15, GIS 90 -> pad 104
    if 0 <= _pad_idx < TOTAL_DOMAINS:
        BACKGROUND_EROSION_RATES[_pad_idx] = _rate
    else:
        print(f"WARNING: DOMAIN_BE_RATES GIS ID {_gis_id} -> pad index {_pad_idx} "
              f"out of range (0-{TOTAL_DOMAINS - 1}). Skipped.")

USE_BACKGROUND_EROSION = any(r != 0.0 for r in BACKGROUND_EROSION_RATES)

assert len(BACKGROUND_EROSION_RATES) == TOTAL_DOMAINS, (
    f"BACKGROUND_EROSION_RATES has {len(BACKGROUND_EROSION_RATES)} entries, "
    f"expected {TOTAL_DOMAINS}."
)

print(f"Background erosion array: {len(BACKGROUND_EROSION_RATES)} domains "
      f"| non-zero: {sum(1 for v in BACKGROUND_EROSION_RATES if v != 0)}  "
      f"| USE_BACKGROUND_EROSION={USE_BACKGROUND_EROSION}")
print(f"Simulation configuration: WAVE_HEIGHTS={WAVE_HEIGHTS_TO_TEST}  "
      f"FLIP_SIGN_MODEL={FLIP_SIGN_MODEL}  PLOT_REAL_DOMAINS_ONLY={PLOT_REAL_DOMAINS_ONLY}")
print("=" * 80 + "\n")

# =============================================================================
# SECTION 8: HISTORICAL MANAGEMENT EVENTS (road + nourishment)
# =============================================================================
# Calendar-year-keyed rather than period-keyed: each event fires only if its
# year falls inside [START_YEAR, END_YEAR], so the same event list is safe
# across periods - Period 1 (1984) skips anything dated 2004+, and vice versa.
# That's why these live here rather than inside PERIOD_CONFIG.

# --- Historical road management events ---
# Two event types, applied inside the time step loop:
#   "relocate" - road physically moved landward at the known year; ADDS
#                relocation_displacement_m to the setback each listed GIS domain
#                is currently carrying, then management continues normally from
#                there. A displacement, not an absolute setback -- see the long
#                note above HISTORICAL_ROAD_EVENTS for why that distinction is
#                what kept NC-12 out of Pamlico Sound.
#   "bridge"   - bridge/alternate route replaced the road surface; disables
#                RoadwayManager for the listed domains from that year onward
#                (fully natural dynamics). Applied in ALL periods/scenarios
#                as a fixed boundary condition (the bridge is already built).
#
# Master toggle: set False to skip ALL "relocate"-type events below (e.g. if
# a relocated setback crashes the run and you want to isolate whether that's
# the cause, or run a no-relocation baseline). "bridge" events are unaffected
# by this flag — Jug Handle Bridge is a fixed boundary condition, not a manual
# relocation, and carries no setback value that could trigger the same issue.
ENABLE_HISTORICAL_ROAD_RELOCATIONS = False

# Format: list of dicts with keys:
#   year      : calendar year the event occurs
#   gis_start : first GIS domain affected (inclusive)
#   gis_end   : last GIS domain affected (inclusive)
#   type      : "relocate" | "bridge"
#   note      : human-readable description (for metadata log)
#   enabled   : optional, defaults to True. Set False to skip just this one
#               event (e.g. while re-checking its relocation_displacement_m
#               values) without disabling the others via the master toggle.

# RELOCATIONS ARE A DISPLACEMENT, NOT A SETBACK
# ---------------------------------------------
# These events used to carry post_relocation_setback_m -- an ABSOLUTE setback
# built as (1984 setback) + (1978->1997 road displacement). That double-counted
# the shoreline retreat, and it is what drove NC-12 into Pamlico Sound:
#
#   The 1984 setback is measured from the 1984 DUNE LINE. It is spent against a
#   grid whose row 0 is the 2009 DEM dune crest, which sits landward of the 1984
#   dune by the 1984->2009 retreat R. So the origin has ALREADY moved landward
#   by R. Adding the physical displacement on top of a 1984-referenced setback
#   counts R a second time, and the road is placed R metres behind where NCDOT
#   actually put it.
#
#   Checked against RoadSetback_2004.csv -- a same-year (2004 road vs 2004 dune)
#   measurement taken AFTER both events, so it already IS the post-relocation
#   position. The old absolute values overshot it by a median of +35 m in the
#   1999 block and +96 m at Pea Island, where retreat was worst:
#
#       GIS 11: 129 prescribed vs  81 measured   (+48)
#       GIS 84: 163 prescribed vs  50 measured  (+113)
#       GIS 86: 205 prescribed vs  88 measured  (+117)
#
#   At GIS 84 the measured 2004 setback (50 m) is SMALLER than the 1984 setback
#   (93 m) even though the road was moved 70 m landward in 1989 -- the shoreline
#   overtook the road faster than NCDOT could move it. That is the real Pea
#   Island story, and the absolute value erased it.
#
# So the events now carry relocation_displacement_m: how far the road physically
# moved, and nothing else. CASCADE already decrements road_setback by
# dune_migrated every year (roadway_manager.road_relocation_checks), so by the
# event year the model's setback has absorbed the modelled retreat on its own.
# Adding the displacement to the CURRENT setback counts the retreat exactly
# once, from the model's own dune migration.
#
#   new_setback = rm._road_setback + relocation_displacement_m[gis]
#
# The displacements below are the 1978->1997 cross-shore offsets measured in
# ArcGIS Pro -- unchanged, they were never the problem.
HISTORICAL_ROAD_EVENTS = [
    # 1989: Pea Island / N. Rodanthe relocation
    dict(
        year=1989, gis_start=84, gis_end=87, type="relocate",
        note="NC-12 relocated landward 1989, Pea Island (GIS 84-87)",
        enabled=True,
        relocation_displacement_m={
            84:  70.0,
            85: 120.0,
            86: 120.0,
            87:  25.0,
        },
    ),
    # 1999: Inter-village south relocation
    # NOTE: GIS 11 is the domain whose old 129.0 m absolute setback pushed the
    # road off the back of a narrow interior and triggered the bulldoze()
    # IndexError. Under the displacement form it moves +80 m from wherever the
    # model has the road in 1999, which is the physical statement.
    dict(
        year=1999, gis_start=9, gis_end=15, type="relocate",
        note="NC-12 relocated landward 1999, inter-village south (GIS 9-15)",
        enabled=True,
        relocation_displacement_m={
             9:  30.0,
            10:  55.0,
            11:  80.0,
            12:  70.0,
            13:  55.0,
            14:  30.0,
            15:   5.0,
        },
    ),
    # 2022: Jug Handle Bridge - road surface removed, natural dynamics restored
    dict(
        year=2022, gis_start=82, gis_end=88, type="bridge",
        note="Jug Handle Bridge installed 2022; road removed GIS 82-88 -> unmanaged",
    ),
]

# Independent check on a corrected relocation: RoadSetback_2004.csv for the same
# domain. Both events precede 2004, so a corrected setback should land near the
# 2004 same-year measurement. Printed alongside each relocation when available.
RELOCATION_CHECK_2004 = {
     9: 89.0, 10: 83.0, 11: 81.0, 12: 89.0, 13: 87.0, 14: 93.0, 15: 71.0,
    84: 50.0, 85: 85.0, 86: 88.0, 87: 40.0,
}

# --- Historical beach nourishment data (events fall within 2004-2024 only) ---
# Source: Hatteras_Management_Timelines.xlsx -> Nourishment_Timeline sheet
#
# Three nourishment projects fall within the 2004-2024 hindcast window:
#   2014: Rodanthe emergency fill  - GIS 85-88  (4 domains,  1,620,000 cy)
#   2022: Avon shore protection    - GIS 23-26  (4 domains,  2,200,000 cy)
#   2022: Buxton shore protection  - GIS  6-15  (10 domains, 1,200,000 cy)
#
# Volume conversion: 1 cy = 0.764555 m^3; per-domain m^3 = (project cy /
# n_domains) x 0.764555 (computed inline below); per-domain m^3/m = that /
# DOMAIN_LENGTH_M (500 m). GIS -> pad uses the Section 1 formula.
#
# HAT_BN_YEARS must stay sorted chronologically. Every domain's volume list
# needs exactly len(HAT_BN_YEARS) entries (0 = not nourished that year). The
# build function (HELPER FUNCTIONS, below) skips years outside the active
# period, so Period 1 (1984) returns all-zero arrays with no code change.

_CY_TO_M3 = 0.764555   # cubic yards -> cubic metres

# Calendar years of nourishment events - SORTED, one entry per distinct year
HAT_BN_YEARS = [2014, 2022]
#                 ^       ^
#              Rodanthe  Avon + Buxton

# HAT_BN_VOLUME_BY_DOMAIN
# Key   : GIS domain ID (1-90)
# Value : [volume_for_2014_m3, volume_for_2022_m3]
#         0 = not nourished that year; non-zero = total m^3 for that domain
HAT_BN_VOLUME_BY_DOMAIN = {
    # --- Buxton: GIS 6-15 (10 domains, 1,200,000 cy in 2022) ---
    #   2014 volume = 0 (no Buxton event)
    #   2022 volume = 1,200,000 / 10 x 0.764555 = 91,746.6 m^3/domain
     6: [0, round(1_200_000 / 10 * _CY_TO_M3, 1)],
     7: [0, round(1_200_000 / 10 * _CY_TO_M3, 1)],
     8: [0, round(1_200_000 / 10 * _CY_TO_M3, 1)],
     9: [0, round(1_200_000 / 10 * _CY_TO_M3, 1)],
    10: [0, round(1_200_000 / 10 * _CY_TO_M3, 1)],
    11: [0, round(1_200_000 / 10 * _CY_TO_M3, 1)],
    12: [0, round(1_200_000 / 10 * _CY_TO_M3, 1)],
    13: [0, round(1_200_000 / 10 * _CY_TO_M3, 1)],
    14: [0, round(1_200_000 / 10 * _CY_TO_M3, 1)],
    15: [0, round(1_200_000 / 10 * _CY_TO_M3, 1)],
    # --- Avon: GIS 23-26 (4 domains, 2,200,000 cy in 2022) ---
    #   2014 volume = 0 (no Avon event)
    #   2022 volume = 2,200,000 / 4 x 0.764555 = 420,505.3 m^3/domain
    23: [0, round(2_200_000 / 4 * _CY_TO_M3, 1)],
    24: [0, round(2_200_000 / 4 * _CY_TO_M3, 1)],
    25: [0, round(2_200_000 / 4 * _CY_TO_M3, 1)],
    26: [0, round(2_200_000 / 4 * _CY_TO_M3, 1)],
    # --- Rodanthe: GIS 85-88 (4 domains, 1,620,000 cy in 2014) ---
    #   2014 volume = 1,620,000 / 4 x 0.764555 = 309,644.8 m^3/domain
    #   2022 volume = 0 (no Rodanthe event)
    85: [round(1_620_000 / 4 * _CY_TO_M3, 1), 0],
    86: [round(1_620_000 / 4 * _CY_TO_M3, 1), 0],
    87: [round(1_620_000 / 4 * _CY_TO_M3, 1), 0],
    88: [round(1_620_000 / 4 * _CY_TO_M3, 1), 0],
}

# =============================================================================
# SECTION 9: LOAD INPUT DATA
# =============================================================================

print("Loading input data...")

for _label, _path in (("ISLAND_OFFSET_FILE", ISLAND_OFFSET_FILE),
                      ("ROAD_SETBACK_FILE", ROAD_SETBACK_FILE)):
    if not os.path.isfile(_path):
        print(f"CRITICAL ERROR: Missing data file ({_label})")
        print(f"  as given:  {_path}")
        print(f"  abspath:   {os.path.abspath(_path)}")
        sys.exit(1)

try:
    island_offset_all = np.loadtxt(ISLAND_OFFSET_FILE, skiprows=1, delimiter=",")
    island_offset_dam = island_offset_all / 10.0  # m -> dam (single-column file, 1 value per domain)
    print(f"  Loaded island offsets: {island_offset_dam.size} domains (dam)")

    road_setbacks_raw = np.loadtxt(ROAD_SETBACK_FILE, skiprows=1, delimiter=",")
    print(f"  Loaded road setbacks: {road_setbacks_raw.size} values")

except Exception as e:
    print(f"CRITICAL ERROR loading data: {e}")
    sys.exit(1)

# Road setbacks: padded array (zeros outside road span)
road_setbacks_full = np.zeros(TOTAL_DOMAINS)
num_road_values = min(len(road_setbacks_raw), END_ROAD_INDEX - START_ROAD_INDEX)
road_setbacks_full[START_ROAD_INDEX : START_ROAD_INDEX + num_road_values] = \
    road_setbacks_raw[:num_road_values]
print(f"  Road setback array prepared ({TOTAL_DOMAINS} domains)")

# =============================================================================
# PER-DOMAIN ROAD ELEVATION
# =============================================================================
# NOT period-dependent, and deliberately so.
#
# road_ele is a PRESCRIPTION, not an observation: bulldoze() writes it into the
# interior grid every year. So the question it answers is "what grade is NC-12
# built to here", and the answer has to be measured on the surface the model
# actually runs on -- which is the 2009 LiDAR, for BOTH hindcast periods. There
# is one topography, so there is one road elevation.
#
# The 2004 alignment is used for both periods because it is the only digitised
# line that exists on the 2009 surface. Sampled along the 1984 line, GIS 10
# returns 3.98 m NAVD88 -- a dune that migrated over the alignment NCDOT
# abandoned in 1999. The two lines agree to a mean |difference| of 0.00 m over
# the 71 domains where the road never moved, so using 2004 everywhere costs
# nothing where it does not matter. Where it does matter, the 1984 line reads
# HIGHER, not lower: mean 2.37 m against 1.54 m, sd up to 1.70 m, because the
# dune migrated over the abandoned corridor. Both candidates sit above the
# neighbouring un-relocated grade (1.26 m), so 2004 is the lower of the two and
# the only one measured on a graded surface.
#
# Produced by HAT_road_elevation.py: MEAN of the 1 m native LiDAR clip under a
# 3.5 m buffer on nc12_2004.geojson (~7 m corridor, one carriageway).
#
# The ArcGIS elevation_2009 column this was once validated against (mean
# |difference| 0.03 m) NO LONGER EXISTS -- nc12_2004.csv now carries avg_elev_m,
# z_mean, z_max, z_min, road_z and relief_m, every one of them all-zero or
# all-empty. That check cannot be re-run from the repo. What replaces it is
# internal: a corridor-width sweep, a 1 m vs 10 m comparison, an alongshore
# continuity check and the relocation bracket. See RoadElevation_audit.md.
#
# ALREADY MHW-RELATIVE, which is the frame bulldoze() compares against (the
# extractor subtracts MHW_M = 0.36 before anything else). Do NOT subtract again.
# Toggle, overridable from the shell so a paired sensitivity run needs no edit:
#     HAT_PER_DOMAIN_ROAD_ELE=0 python HAT_hindcast_1984_2024.py
# The uniform run is tagged in RUN_NAME_BASE below, so the two do not collide in
# comparison/raw_runs/.
USE_PER_DOMAIN_ROAD_ELEVATION = (
    os.environ.get("HAT_PER_DOMAIN_ROAD_ELE", "1").lower()
    not in ("0", "false", "no")
)

# Buffer domains only. VERIFIED INERT: all 38 padded indices outside the road
# span (0-22, 105-119) have ROADWAY_MANAGEMENT_ON = False, so no RoadwayManager
# is ever constructed for them and this value cannot reach the model. The 27
# village indices inside the road span (GIS 21-31, 68-83) do receive a measured
# elevation from the file but are likewise unmanaged, so it is unused there too.
# 55 of the 82 road-span domains are actually managed.
ROAD_ELEVATION_FALLBACK = 1.45
# One file for both periods -- one 2009 DEM, so one elevation set. No year in
# the name by design; see road_elevation/RoadElevation_audit.md.
ROAD_ELEVATION_FILE = os.path.join(
    HATTERAS_DATA_BASE, "4-mgmt-forcing", "road_elevation",
    "RoadElevation.csv")

road_elevation_full = np.full(TOTAL_DOMAINS, ROAD_ELEVATION_FALLBACK)
if USE_PER_DOMAIN_ROAD_ELEVATION:
    if not os.path.exists(ROAD_ELEVATION_FILE):
        print(f"\n  [warn] road elevation file not found:\n"
              f"         {ROAD_ELEVATION_FILE}\n"
              f"         falling back to the scalar "
              f"{ROAD_ELEVATION_FALLBACK} m for every domain.\n"
              f"         Run HAT_road_elevation.py to create it.")
    else:
        _re_raw = np.loadtxt(ROAD_ELEVATION_FILE, skiprows=1, delimiter=",")
        _n_re = min(len(_re_raw), END_ROAD_INDEX - START_ROAD_INDEX)
        road_elevation_full[START_ROAD_INDEX : START_ROAD_INDEX + _n_re] = \
            _re_raw[:_n_re]
        _road_only = road_elevation_full[START_ROAD_INDEX:START_ROAD_INDEX + _n_re]
        print(f"  Road elevation: per-domain from "
              f"{os.path.basename(ROAD_ELEVATION_FILE)} ({_n_re} domains)")
        print(f"    m MHW-relative: min {_road_only.min():.2f} | "
              f"median {np.median(_road_only):.2f} | max {_road_only.max():.2f}")
        print(f"    m NAVD88      : min {_road_only.min() + MHW_ELEVATION:.2f} | "
              f"median {np.median(_road_only) + MHW_ELEVATION:.2f} | "
              f"max {_road_only.max() + MHW_ELEVATION:.2f}")
        print(f"    vs the old scalar {ROAD_ELEVATION_FALLBACK}: median change "
              f"{np.median(_road_only) - ROAD_ELEVATION_FALLBACK:+.2f} m, "
              f"lower in {int((_road_only < ROAD_ELEVATION_FALLBACK).sum())} of "
              f"{_n_re} domains")
        print(f"    bulldoze() writes this into the interior at t=0, so a value "
              f"near the real surface")
        print(f"    perturbs the initial grid far less than a uniform one did.")
else:
    print(f"  Road elevation: uniform {ROAD_ELEVATION_FALLBACK} m "
          f"(USE_PER_DOMAIN_ROAD_ELEVATION = False)")
    # Tag the uniform run so a paired sensitivity pair lands in two folders.
    RUN_NAME_BASE = f"{RUN_NAME_BASE}_uniformRoadEle"
    print(f"    run name tagged -> {RUN_NAME_BASE}")

# =============================================================================
# UNITS / DATUM ASSERTIONS
# =============================================================================
# CASCADE takes TWO conventions side by side, and the difference is invisible
# in a value that merely looks reasonable:
#
#   CONVERTED FOR YOU  supply in m NAVD88; load_input.py does /10 - MHW
#       MHW (:227), BermEl (:241), Dmaxel (:304), ShrubEl_* (:346-7)
#
#   ALREADY CONVERTED  supply in the model's own frame; nothing touches it
#       elevation_file  dam, MHW-relative   (load_elevation() is just np.load;
#                       configuration.py documents "[dam x dam x dam MHW]")
#       dune_file       dam, ABOVE BERM     (load_input.py:249 assigns straight
#                       into DuneDomain)
#       road_ele        m, MHW-relative     (bulldoze() only does road_ele/dz)
#
# barrier3d.py:1197 pops MHW into self._MHW and _MHW appears NOWHERE else in
# the file -- the model never applies it to a grid. Anything grid-shaped has to
# arrive pre-converted, and a slip produces a plausible-looking run rather than
# an error. These assertions make it fail at launch instead.
#
# Ranges are deliberately loose: they discriminate a 10x unit error or a 0.36 m
# datum shift, not whether the values are physically ideal. Run
# scripts/input_prep/4-mgmt-forcings/road_offset/HAT_units_datum_check.py for
# the full trace and the plausibility checks.
def _assert_units_and_datum():
    _topo_sentinel_dam = -0.30      # SENTINEL_WATER_M (-3.0 m MHW) in dam

    _t = np.load(ELEVATION_FILE_PATHS[START_REAL_INDEX], mmap_mode="r")
    assert -0.35 <= float(np.nanmin(_t)) <= 0.0, (
        f"topography min {float(np.nanmin(_t)):.3f} dam: expected dam "
        f"MHW-relative with a {_topo_sentinel_dam} dam water sentinel. "
        f"A metres array would be ~10x.")
    assert 0.05 < float(np.nanmax(_t)) < 2.0, (
        f"topography max {float(np.nanmax(_t)):.3f} dam "
        f"= {float(np.nanmax(_t)) * 10:.1f} m MHW: expected dam. A metres "
        f"array would read ~{float(np.nanmax(_t)) * 10:.0f} here.")

    _d = np.load(DUNE_FILE_PATHS[START_REAL_INDEX], mmap_mode="r")
    _dv = np.asarray(_d, dtype=float).ravel()
    _dv = _dv[_dv > _topo_sentinel_dam + 1e-9]
    if _dv.size:
        assert float(np.nanmax(_dv)) < 1.5, (
            f"dune max {float(np.nanmax(_dv)):.3f} dam "
            f"= {float(np.nanmax(_dv)) * 10:.1f} m above berm: expected dam "
            f"ABOVE BERM. A metres array would be ~10x.")

    _r = road_elevation_full[START_ROAD_INDEX:END_ROAD_INDEX]
    assert 0.0 < float(_r.min()) and float(_r.max()) < 3.0, (
        f"road_ele range {_r.min():.2f} to {_r.max():.2f}: expected metres "
        f"MHW-relative. NAVD88 would sit ~{MHW_ELEVATION:.2f} m higher; "
        f"decametres would be ~10x smaller.")

    # The extractor's berm and Barrier3D's must be the same number, or the dune
    # heights are measured from a different datum than the model rebuilds to.
    _berm_b3d = (BERM_ELEVATION / 10.0 - MHW_ELEVATION / 10.0) * 10.0
    assert abs(_berm_b3d - (1.70 - MHW_ELEVATION)) < 1e-9, (
        f"berm mismatch: load_input gives {_berm_b3d:.3f} m MHW, the extractor "
        f"used BERM_ELEV_NAVD_M - MHW_M = {1.70 - MHW_ELEVATION:.3f}. "
        f"Dune heights are measured above the extractor's berm.")
    print("  Units/datum assertions passed "
          "(topography dam MHW, dunes dam above berm, road_ele m MHW)")


# NOTE: called AFTER ELEVATION_FILE_PATHS / DUNE_FILE_PATHS are populated
# (see the call site further down); defining it here keeps it next to the
# convention it documents.

# ROADWAY_MANAGEMENT_ON: per-domain flag passed directly to Cascade.
# Non-village road domains (GIS 9-90, excl. villages): management ON
#   -> RoadwayManager runs each step: bulldoze, rebuild dunes, relocate.
# Village domains (Buxton GIS 7-8, Avon GIS 21-31, Tri-Village GIS 68-83):
#   management OFF -> no RoadwayManager, no relocation possible.
# Buffer domains: always OFF.
# ENABLE_ROADWAY_MANAGEMENT acts as a global on/off switch.
ROADWAY_MANAGEMENT_ON = [False] * TOTAL_DOMAINS
if ENABLE_ROADWAY_MANAGEMENT:
    for _pad in range(START_ROAD_INDEX, END_ROAD_INDEX):
        ROADWAY_MANAGEMENT_ON[_pad] = True   # road span default: on
    for _gis_start, _gis_end in VILLAGE_GIS_RANGES_NO_RELOCATION:
        for _gis in range(_gis_start, _gis_end + 1):
            _pad = START_REAL_INDEX + (_gis - FIRST_FILE_NUMBER)
            if 0 <= _pad < TOTAL_DOMAINS:
                ROADWAY_MANAGEMENT_ON[_pad] = False  # village: off

_n_road_managed = sum(ROADWAY_MANAGEMENT_ON)
_n_village_off  = sum(
    1 for _gs, _ge in VILLAGE_GIS_RANGES_NO_RELOCATION
    for _g in range(_gs, _ge + 1)
    if 0 <= START_REAL_INDEX + (_g - FIRST_FILE_NUMBER) < TOTAL_DOMAINS
)
print(
    f"  ROADWAY_MANAGEMENT_ON: {_n_road_managed} managed domains - "
    f"village domains off ({_n_village_off}: no bulldoze/rebuild/relocation)"
)

SANDBAG_MANAGEMENT_ON     = [ENABLE_SANDBAG_PLACEMENT]  * TOTAL_DOMAINS
NOURISHMENT_MANAGEMENT_ON = [ENABLE_NOURISHMENT]        * TOTAL_DOMAINS

# =============================================================================
# SECTION 10: ELEVATION + DUNE FILE LISTS
# =============================================================================

print("Generating elevation + dune profile file paths...")

ELEVATION_FILE_PATHS = []
DUNE_FILE_PATHS      = []

for _ in range(START_REAL_INDEX):  # left buffers
    DUNE_FILE_PATHS.append(os.path.join(BUFFER_DIR, "sample_1_dune.npy"))
    ELEVATION_FILE_PATHS.append(os.path.join(BUFFER_DIR, "sample_1_topography.npy"))

for i_list in range(START_REAL_INDEX, END_REAL_INDEX):  # real domains
    file_num = FIRST_FILE_NUMBER + (i_list - START_REAL_INDEX)
    DUNE_FILE_PATHS.append(
        os.path.join(DUNE_TOPO_DIR, "dunes",
                     f"domain_{file_num}_dune_{TOPO_DUNE_INIT_YEAR}.npy"))
    ELEVATION_FILE_PATHS.append(
        os.path.join(DUNE_TOPO_DIR, "topography",
                     f"domain_{file_num}_topography_{TOPO_DUNE_INIT_YEAR}.npy"))

for _ in range(END_REAL_INDEX, TOTAL_DOMAINS):  # right buffers
    DUNE_FILE_PATHS.append(os.path.join(BUFFER_DIR, "sample_1_dune.npy"))
    ELEVATION_FILE_PATHS.append(os.path.join(BUFFER_DIR, "sample_1_topography.npy"))

# Fail here rather than deep inside Barrier3D's init with an opaque traceback.
_missing = [p for p in ELEVATION_FILE_PATHS + DUNE_FILE_PATHS if not os.path.exists(p)]
if _missing:
    raise FileNotFoundError(
        f"{len(_missing)} of {2 * TOTAL_DOMAINS} init files missing - check "
        f"TOPO_DUNE_VERSION ({TOPO_DUNE_VERSION!r}) and DUNE_TOPO_DIR.\n"
        f"  first missing: {_missing[0]}"
    )

print(f"  Generated {len(ELEVATION_FILE_PATHS)} elevation file paths")
print(f"  Generated {len(DUNE_FILE_PATHS)} dune file paths")

# Units/datum gate. Defined above next to the convention it documents; called
# here because it reads the init arrays. A datum or unit slip in any of these
# produces a plausible-looking run rather than an error, so fail at launch.
try:
    _assert_units_and_datum()
except AssertionError as _e:
    print(f"\n  [FATAL] units/datum check failed:\n         {_e}\n")
    raise
except (IndexError, FileNotFoundError, OSError) as _e:
    print(f"  [warn] units/datum assertions skipped ({type(_e).__name__}: {_e})")

print("=" * 80 + "\n")

# =============================================================================
# HELPER FUNCTIONS
# =============================================================================

def build_nourishment_arrays_from_manual_inputs():
    """
    Build per-year nourishment-on and volume arrays for the CASCADE time loop,
    from HAT_BN_YEARS + HAT_BN_VOLUME_BY_DOMAIN (Section 8). Years outside
    [START_YEAR, END_YEAR] are silently skipped, so Period 1 returns all-zero
    arrays with no code change needed.

    Returns
    -------
    nourishment_on_by_year     : dict {year: np.ndarray[TOTAL_DOMAINS]}, 1/0
    nourishment_volume_by_year : dict {year: list[TOTAL_DOMAINS]}, m^3/m
    """
    nourishment_on_by_year     = {}
    nourishment_volume_by_year = {}

    for year in range(START_YEAR, END_YEAR + 1):
        nourishment_on_by_year[year]     = np.zeros(TOTAL_DOMAINS)
        nourishment_volume_by_year[year] = [0.0] * TOTAL_DOMAINS

    for gis_id, volumes_m3 in HAT_BN_VOLUME_BY_DOMAIN.items():
        if len(HAT_BN_YEARS) != len(volumes_m3):
            raise ValueError(
                f"GIS domain {gis_id}: HAT_BN_YEARS and volume list must have "
                f"the same length ({len(HAT_BN_YEARS)} vs {len(volumes_m3)})."
            )

        pad_idx = _gis_to_pad(gis_id)
        if not (0 <= pad_idx < TOTAL_DOMAINS):
            print(f"  WARNING: GIS {gis_id} -> pad {pad_idx} out of range - skipped.")
            continue

        for year, total_m3 in zip(HAT_BN_YEARS, volumes_m3):
            if year < START_YEAR or year > END_YEAR:
                continue   # event outside this period - skip silently

            volume_m3_per_m = float(total_m3) / DOMAIN_LENGTH_M
            nourishment_on_by_year[year][pad_idx]     = 1
            nourishment_volume_by_year[year][pad_idx] = volume_m3_per_m

    # Print schedule summary
    has_events = False
    for year in range(START_YEAR, END_YEAR + 1):
        active_pad = np.where(nourishment_on_by_year[year] == 1)[0]
        if len(active_pad) > 0:
            has_events = True
            active_gis = [
                FIRST_FILE_NUMBER + (idx - START_REAL_INDEX)
                for idx in active_pad
                if START_REAL_INDEX <= idx < END_REAL_INDEX
            ]
            total_vol = (
                np.sum(np.asarray(nourishment_volume_by_year[year], dtype=float))
                * DOMAIN_LENGTH_M
            )
            print(f"  {year}: GIS domains {active_gis}  |  total = {total_vol:,.0f} m^3")

    if not has_events:
        print("  (no nourishment events in this period's date range)")

    return nourishment_on_by_year, nourishment_volume_by_year


def get_x_s_TS(b3d):
    """Extract shoreline time series from a Barrier3D object."""
    if hasattr(b3d, "x_s_TS"):
        return np.asarray(b3d.x_s_TS, dtype=float)
    if hasattr(b3d, "_x_s_TS"):
        return np.asarray(b3d._x_s_TS, dtype=float)
    raise AttributeError(
        "No shoreline time series found (x_s_TS / _x_s_TS) on Barrier3D object."
    )


def build_shoreline_matrix(cascade, to_meters=True):
    """Build [time x domain] shoreline matrix from a CASCADE object."""
    b3d_list = cascade.barrier3d
    ndom = len(b3d_list)
    nt   = len(get_x_s_TS(b3d_list[0]))
    shoreline = np.zeros((nt, ndom), dtype=float)
    for j in range(ndom):
        shoreline[:, j] = get_x_s_TS(b3d_list[j])
    if to_meters:
        shoreline *= 10.0  # dam -> m
    return shoreline


def estimate_transect_spacing(along_coast_m):
    """Median spacing between consecutive transects in metres (positive diffs only)."""
    arr   = np.sort(along_coast_m)
    diffs = np.diff(arr)
    pos   = diffs[diffs > 0]
    return float(np.median(pos)) if len(pos) else 50.0


def load_transect_data(ds):
    """
    Load individual transect LRR values from transect_lrr_full.csv.
    Derives along-coast distance by spreading each domain's transects evenly
    across its 500 m band (mirrors coastsat_smoothed_method_comparison.py).

    Returns
    -------
    domain_ids    : int array   - GIS domain ID (1-90) for each transect
    lrr_values    : float array - LRR (m/yr) for each transect
    along_coast_m : float array - cumulative along-coast distance (m)
    All None on load failure.
    """
    csv_path   = ds["csv"]
    domain_col = ds["domain_col"]
    rate_col   = ds["rate_col"]
    id_col     = ds.get("transect_id_col", "transect_id")

    if not os.path.exists(csv_path):
        print(f"  WARNING: Transect CSV not found: {csv_path}")
        return None, None, None

    df = pd.read_csv(csv_path)
    df.columns = [c.split(".csv")[-1] if ".csv" in c else c for c in df.columns]

    for col in [domain_col, rate_col]:
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors="coerce")

    df = df.dropna(subset=[domain_col, rate_col])
    df[domain_col] = df[domain_col].astype(int)
    df = df[(df[domain_col] >= FIRST_FILE_NUMBER) & (df[domain_col] <= LAST_FILE_NUMBER)]

    sort_cols = [domain_col, id_col] if id_col in df.columns else [domain_col]
    df = df.sort_values(sort_cols).reset_index(drop=True)

    # Derive along_coast_m: spread each domain's transects evenly across its 500 m band.
    # Vectorised - avoids groupby.apply and the associated pandas deprecation warning.
    df["_rank"]         = df.groupby(domain_col).cumcount()
    df["_n"]            = df.groupby(domain_col)[domain_col].transform("count")
    df["along_coast_m"] = (
        (df[domain_col] - 1) * DOMAIN_SPACING_M
        + (df["_rank"] + 0.5) * (DOMAIN_SPACING_M / df["_n"])
    )
    df = df.drop(columns=["_rank", "_n"])

    domain_ids    = df[domain_col].values.astype(int)
    lrr_values    = df[rate_col].values.astype(float)
    along_coast_m = df["along_coast_m"].values.astype(float)

    spacing = estimate_transect_spacing(along_coast_m)
    print(f"  {ds['label']}: {len(df)} transects  "
          f"est. spacing {spacing:.0f} m  "
          f"LRR range {np.nanmin(lrr_values):+.2f}-{np.nanmax(lrr_values):+.2f} m/yr")
    return domain_ids, lrr_values, along_coast_m


def loess_smooth_transect_to_domains(along_coast_m, lrr, domain_ids, window_domains):
    """
    Apply LOESS at transect resolution using physical along-coast distance (m),
    then aggregate smoothed values to GIS domain resolution (mean per domain).

    Returns
    -------
    gis_x    : int array   - GIS domain IDs (1-90) that have at least one transect
    smoothed : float array - domain-averaged smoothed LRR (m/yr)
    frac     : float       - LOESS frac used (for logging)
    """
    window_km = window_domains * DOMAIN_SPACING_M / 1000.0
    spacing_m = estimate_transect_spacing(along_coast_m)
    n         = len(along_coast_m)
    frac      = float(np.clip((window_km * 1000.0 / spacing_m) / n, 0.02, 1.0))

    valid = np.isfinite(lrr)
    if valid.sum() < 5:
        print(f"  WARNING: Too few valid transects ({valid.sum()}) - skipping LOESS")
        return None, None, frac

    result            = lowess(lrr[valid], along_coast_m[valid], frac=frac, return_sorted=True)
    smoothed_t        = np.full(n, np.nan)
    smoothed_t[valid] = np.interp(along_coast_m[valid], result[:, 0], result[:, 1])

    dom_agg = (pd.DataFrame({"domain": domain_ids, "smoothed": smoothed_t})
                 .groupby("domain")["smoothed"].mean()
                 .dropna())

    return dom_agg.index.values.astype(int), dom_agg.values, frac


def compute_domain_means(domain_ids, lrr_values, gis_min, gis_max):
    """
    Mean LRR per GIS domain for domains in [gis_min, gis_max]. Used to
    substitute raw per-domain averages for LOESS smoothing in the
    southernmost domains, where boundary effects dominate.

    Returns
    -------
    gis_x : int array   - GIS domain IDs with at least one transect
    means : float array - mean LRR per domain
    """
    df  = pd.DataFrame({"domain": domain_ids, "lrr": lrr_values})
    sub = df[(df["domain"] >= gis_min) & (df["domain"] <= gis_max) & df["lrr"].notna()]
    if len(sub) == 0:
        return np.array([], dtype=int), np.array([], dtype=float)
    agg = sub.groupby("domain")["lrr"].mean()
    return agg.index.values.astype(int), agg.values


def splice_loess_with_raw_south(win_gis_x, win_smoothed,
                                transect_domain_ids, transect_lrr_values,
                                skip_n=LOESS_SKIP_SOUTHERN_DOMAINS,
                                is_widest_window=False):
    """
    Return (plot_x, plot_y) for a single LOESS window with optional southern splice.

    Widest window (is_widest_window=True) with skip_n > 0:
      domains 1-skip_n -> per-domain raw means; domains skip_n+1+ -> LOESS,
      concatenated so the plotted line is continuous.
    Narrower windows or skip_n == 0:
      domains 1-skip_n omitted entirely (line starts at skip_n+1).

    Returns
    -------
    plot_x : int/float array - GIS domain IDs to plot
    plot_y : float array     - LRR values to plot
    """
    if skip_n == 0:
        return win_gis_x, win_smoothed

    # LOESS portion: domains strictly north of skip_n
    mask   = win_gis_x > skip_n
    lx, ly = win_gis_x[mask], win_smoothed[mask]

    # All windows: LOESS line starts at domain skip_n+1.
    # Domains 1-skip_n show raw scatter only - no smoothed line.
    plot_x, plot_y = lx, ly

    return plot_x, plot_y


# =============================================================================
# ANNOTATED FIGURE HELPERS
# =============================================================================

def add_geographic_annotations(ax):
    """
    Add all geographic reference annotations to an axis (bottom -> top):
      1. Wimble Shoals influence zone (hatched amber fill, bottom label)
      2. Community shaded spans       (steel-blue fill, top labels)
      3. Village center lines         (dashed gray,  y=0.84)
      4. Pier lines                   (dash-dot blue, y=0.76, rotated)
      5. Groin lines                  (dotted red,    y=0.76, rotated)

    Label y-positions use blended axes-fraction coords (fixed regardless of
    data range). X-axis must be in GIS domain IDs (1-90).
    """
    trans = blended_transform_factory(ax.transData, ax.transAxes)

    # 1a. Avon Shoals influence zone (drawn first, beneath Wimble Shoals and everything else)
    alo, ahi = ANN_AVON_SHOALS
    ax.axvspan(alo - 0.5, ahi + 0.5,
               color=ANN_C_AVON_SHOALS, alpha=0.10, zorder=0,
               hatch="///", edgecolor=ANN_C_AVON_SHOALS, linewidth=0)
    ax.text((alo + ahi) / 2.0, 0.04,
            "Avon Shoals\nPosition", transform=trans,
            ha="center", va="bottom", fontsize=7, color="#7A5800", style="italic",
            bbox=dict(boxstyle="round,pad=0.2", fc="white", ec="none", alpha=0.80))

    # 1b. Wimble Shoals influence zone
    wlo, whi = ANN_WIMBLE_SHOALS
    ax.axvspan(wlo - 0.5, whi + 0.5,
               color=ANN_C_WIMBLE, alpha=0.10, zorder=0,
               hatch="///", edgecolor=ANN_C_WIMBLE, linewidth=0)
    ax.text((wlo + whi) / 2.0, 0.04,
            "Wimble Shoals\nPosition", transform=trans,
            ha="center", va="bottom", fontsize=7, color="#7A5800", style="italic",
            bbox=dict(boxstyle="round,pad=0.2", fc="white", ec="none", alpha=0.80))

    # 2. Community / town spans
    for span_label, (d_lo, d_hi) in ANN_TOWN_SPANS.items():
        ax.axvspan(d_lo - 0.5, d_hi + 0.5,
                   color=ANN_C_TOWN_SPAN, alpha=0.14, zorder=0)
        ax.text((d_lo + d_hi) / 2.0, 0.90,
                span_label, transform=trans,
                ha="center", va="top", fontsize=8, color="0.25", fontweight="bold",
                bbox=dict(boxstyle="round,pad=0.2", fc="white", ec="none", alpha=0.85))

    # 3. Village center lines (within Tri-Village span)
    for vname, dom in ANN_VILLAGE_LINES.items():
        ax.axvline(dom, color=ANN_C_VILLAGE_LINE, lw=0.9, ls="--", alpha=0.65, zorder=1)
        ax.text(dom, 0.84, vname, transform=trans,
                ha="center", va="top", fontsize=7.5, color="0.30",
                bbox=dict(boxstyle="round,pad=0.15", fc="white", ec="none", alpha=0.80))

    # 4. Pier lines - per-pier label y from (domain, label_y) tuple
    for pname, (dom, lbl_y) in ANN_PIERS.items():
        ax.axvline(dom, color=ANN_C_PIER, lw=1.0, ls="-.", alpha=0.80, zorder=2)
        ax.text(dom, lbl_y, pname, transform=trans,
                ha="center", va="top", fontsize=7, color=ANN_C_PIER, rotation=90,
                bbox=dict(boxstyle="round,pad=0.15", fc="white", ec="none", alpha=0.80))

    # 5. Groin lines
    for gname, dom in ANN_GROINS.items():
        ax.axvline(dom, color=ANN_C_GROIN, lw=1.1, ls=":", alpha=0.85, zorder=2)
        ax.text(dom, ANN_GROIN_LABEL_Y, gname, transform=trans,
                ha="center", va="top", fontsize=7, color=ANN_C_GROIN, rotation=90,
                bbox=dict(boxstyle="round,pad=0.15", fc="white", ec="none", alpha=0.80))


def annotation_legend_handles():
    """Return proxy legend artists for the geographic annotation layers."""
    return [
        Patch(fc=ANN_C_TOWN_SPAN, alpha=0.30, label="Community"),
        Patch(fc=ANN_C_WIMBLE, alpha=0.25, hatch="///",
              edgecolor=ANN_C_WIMBLE, linewidth=0, label="Shoals position (Avon / Wimble)"),
        Line2D([0], [0], color=ANN_C_VILLAGE_LINE, lw=0.9, ls="--", label="Village center"),
        Line2D([0], [0], color=ANN_C_PIER, lw=1.0, ls="-.", label="Pier"),
        Line2D([0], [0], color=ANN_C_GROIN, lw=1.1, ls=":", label="Groin"),
    ]


# =============================================================================
# SHORELINE ANIMATION (GIF)
# =============================================================================

def _open_file(path):
    """Open a file in the OS default application (cross-platform, non-fatal)."""
    import subprocess
    try:
        if sys.platform.startswith("win"):
            os.startfile(path)  # type: ignore[attr-defined]  # noqa
        elif sys.platform == "darwin":
            subprocess.run(["open", path], check=False)
        else:
            subprocess.run(["xdg-open", path], check=False)
    except Exception as e:
        print(f"  [GIF] could not auto-open ({e}); open it manually: {path}")


def _slug(text):
    """Filename-safe short label (e.g. 'Buxton Groin' -> 'BuxtonGroin')."""
    return "".join(ch for ch in str(text) if ch.isalnum())


def _groin_zoom_window(gdom, pad=9):
    """
    Inclusive GIS window centred on one groin, clipped to the real domains.
    ANN_GROINS values sit on a domain boundary (e.g. 5.5 = the 5/6 interface),
    so the window is built around that interface, not a cell centre.
    """
    lo = int(max(FIRST_FILE_NUMBER, np.floor(gdom - pad)))
    hi = int(min(LAST_FILE_NUMBER,  np.ceil(gdom + pad)))
    return lo, hi


def _select_groins(which="all"):
    """
    Resolve a job's "which" key into [(name, domain), ...] from ANN_GROINS.

    which : "all" (default, every groin — new structures are picked up with no
            config change), a single name, or a list of names.
    """
    if not ANN_GROINS:
        raise ValueError("range='groin' requires at least one entry in ANN_GROINS.")

    if isinstance(which, str) and which.strip().lower() == "all":
        return list(ANN_GROINS.items())

    names = [which] if isinstance(which, str) else list(which)
    missing = [n for n in names if n not in ANN_GROINS]
    if missing:
        raise ValueError(
            f"groin(s) {missing} not in ANN_GROINS (have: {list(ANN_GROINS)})."
        )
    return [(n, ANN_GROINS[n]) for n in names]


def _resolve_gif_domain_range(spec, pad=9, groin=None):
    """
    Turn a GIF_JOBS "range" value into concrete plotting coordinates.

    groin : (name, domain) tuple, required when spec == "groin". Supplied by
            make_all_shoreline_gifs, which fans one "groin" job out into one
            call per structure.

    Returns
    -------
    pad_lo, pad_hi : padded-array slice bounds (pad_hi exclusive)
    axis_kind      : "gis" (x-axis in GIS domain IDs) or "pad" (padded index)
    x_lo, x_hi     : inclusive x-axis limits in whichever space axis_kind names
    tag            : short string for the comparison filename
    """
    if isinstance(spec, str):
        key = spec.strip().lower()
        if key == "real":
            return (START_REAL_INDEX, END_REAL_INDEX, "gis",
                    FIRST_FILE_NUMBER, LAST_FILE_NUMBER,
                    f"D{FIRST_FILE_NUMBER}-{LAST_FILE_NUMBER}")
        if key == "all":
            return (0, TOTAL_DOMAINS, "pad",
                    0, TOTAL_DOMAINS - 1,
                    f"ALL{TOTAL_DOMAINS}pad")
        if key == "groin":
            if groin is None:
                raise ValueError(
                    "range='groin' must be resolved per-structure; call it via "
                    "make_all_shoreline_gifs, or pass groin=(name, domain)."
                )
            gname, gdom = groin
            gis_lo, gis_hi = _groin_zoom_window(gdom, pad)
            return (_gis_to_pad(gis_lo), _gis_to_pad(gis_hi) + 1, "gis",
                    gis_lo, gis_hi, f"groinZoom_{_slug(gname)}_D{gis_lo}-{gis_hi}")
        if key == "groin_span":
            # One window enclosing EVERY groin — useful once structures interact.
            doms = list(ANN_GROINS.values())
            if not doms:
                raise ValueError("range='groin_span' requires entries in ANN_GROINS.")
            gis_lo, _ = _groin_zoom_window(min(doms), pad)
            _, gis_hi = _groin_zoom_window(max(doms), pad)
            return (_gis_to_pad(gis_lo), _gis_to_pad(gis_hi) + 1, "gis",
                    gis_lo, gis_hi, f"groinSpan_D{gis_lo}-{gis_hi}")
        raise ValueError(
            f"GIF range={spec!r} not understood; use 'real', 'all', 'groin', "
            f"'groin_span', or (lo, hi)."
        )

    gis_lo, gis_hi = int(spec[0]), int(spec[1])
    if gis_lo > gis_hi:
        gis_lo, gis_hi = gis_hi, gis_lo
    if gis_lo < FIRST_FILE_NUMBER or gis_hi > LAST_FILE_NUMBER:
        raise ValueError(
            f"GIF range=({gis_lo}, {gis_hi}) falls outside the real domains "
            f"{FIRST_FILE_NUMBER}-{LAST_FILE_NUMBER}. Use 'all' to include buffers."
        )
    return (_gis_to_pad(gis_lo), _gis_to_pad(gis_hi) + 1, "gis",
            gis_lo, gis_hi, f"D{gis_lo}-{gis_hi}")


def make_shoreline_gif(shoreline_m, run_dir, run_name, Hs=None,
                       domain_range="real", mode="displacement", pad=9, groin=None,
                       baseline_m=None, baseline_label=None,
                       fps=None, stride=None,
                       annotate=None, auto_open=None, keep_frames=None):
    """
    Plan-view animation of the shoreline, one frame per model year. Current
    shoreline drawn over a year-0 reference (dashed grey); shaded blue where
    seaward of the reference, red where landward. Axis orientation is set by
    GIF_OCEAN_AT_BOTTOM (Section 2a): ocean at bottom / sound at top by default,
    matching Hatteras' real cross-shore layout.

    SIGN CONVENTION: build_shoreline_matrix() returns raw x_s_TS (x10 for m),
    NO sign flip applied -- in BRIE/Barrier3D, x_s_TS INCREASES as the shoreline
    retreats landward. FLIP_SIGN_MODEL is applied internally here so larger =
    more seaward, matching change_rate's sign handling in the runner. Do not
    pre-flip shoreline_m before passing it in.

    MODE:
      "displacement" (default) - each domain's change from its OWN year-0
        position; all start at 0, so a few m/yr stays legible across all 90
        domains. Use for the full island.
      "position" - cross-shore position relative to year-0 alongshore mean.
        Keeps planform shape, but x_s_TS is referenced to each domain's own
        local grid, so the alongshore "shape" is partly an artefact of how the
        initial DEMs were extracted, and ~80 m of 40-yr change can be nearly
        invisible against a planform range of several hundred m. Fine for a
        narrow groin-zone window like (1, 30); use "displacement" for the
        full island.
      "difference" - this run minus baseline_m, in displacement terms:
        (run - run_yr0) - (base - base_yr0). Differencing displacements (not
        raw positions) cancels any per-domain initial-condition offset between
        runs, so year 0 is exactly zero. With a no-groin baseline this isolates
        the groin's effect: blue = seaward of baseline (updrift fillet),
        red = landward (downdrift notch).

    Parameters
    ----------
    shoreline_m : 2-D array [n_years, TOTAL_DOMAINS] from
        build_shoreline_matrix(cascade, to_meters=True); raw x_s_TS convention.
    domain_range : "real" | "all" | "groin" | (gis_lo, gis_hi).
    pad : half-width in domains; only read when domain_range="groin".
    baseline_m : same-shape matrix from a previous run (same raw x_s_TS
        convention); required for mode="difference".
    Everything else defaults to the Section 2a config constants.
    """
    import io

    try:
        from PIL import Image
    except Exception as e:      # Pillow ships with matplotlib, so this is unlikely
        print(f"  [GIF] Pillow unavailable ({e}); skipping animation.")
        return None

    fps         = GIF_FPS         if fps         is None else fps
    stride      = GIF_YEAR_STRIDE if stride      is None else stride
    annotate    = GIF_ANNOTATE    if annotate    is None else annotate
    auto_open   = GIF_AUTO_OPEN   if auto_open   is None else auto_open
    keep_frames = GIF_KEEP_FRAMES if keep_frames is None else keep_frames
    baseline_label = GIF_BASELINE_LABEL if baseline_label is None else baseline_label

    mode = str(mode).strip().lower()
    if mode not in ("position", "displacement", "difference"):
        raise ValueError(
            f"GIF mode={mode!r} must be 'position', 'displacement', or 'difference'."
        )
    if mode == "difference" and baseline_m is None:
        print("  [GIF] mode='difference' needs a baseline run "
              "(set GIF_BASELINE_NPY); skipping this job.")
        return None

    pad_lo, pad_hi, axis_kind, x_lo, x_hi, range_tag = _resolve_gif_domain_range(
        domain_range, pad=pad, groin=groin
    )

    raw = np.asarray(shoreline_m, dtype=float)
    if raw.ndim != 2 or raw.shape[1] < pad_hi:
        print(f"  [GIF] shoreline matrix shape {raw.shape} too small for "
              f"padded slice [{pad_lo}:{pad_hi}]; skipping animation.")
        return None

    # x_s_TS increases landward -> flip so that larger = more seaward.
    flip = -1.0 if FLIP_SIGN_MODEL else 1.0
    pos  = (flip * raw)[:, pad_lo:pad_hi]
    n_years = pos.shape[0]
    if n_years < 2:
        print("  [GIF] fewer than 2 model years; skipping animation.")
        return None

    x    = np.arange(x_lo, x_hi + 1)
    init = pos[0].copy()

    if mode == "difference":
        base = np.asarray(baseline_m, dtype=float)
        if base.shape != raw.shape:
            print(f"  [GIF] baseline shape {base.shape} != run shape {raw.shape}; "
                  f"the two runs must share domain count and length. Skipping.")
            return None
        base_pos = (flip * base)[:, pad_lo:pad_hi]
        # Difference the DISPLACEMENTS so any initial-condition raw_offset cancels.
        series = (pos - init[None, :]) - (base_pos - base_pos[0][None, :])
        _up_word = "landward" if GIF_OCEAN_AT_BOTTOM else "seaward"
        ylabel = f"\u0394 shoreline vs {baseline_label} (m)\n{_up_word} \u25b2"
        lbl_up, lbl_dn = (f"Seaward of {baseline_label}", f"Landward of {baseline_label}")
        lbl_ref = "No difference"
    elif mode == "displacement":
        series = pos - init[None, :]
        _up_word = "landward" if GIF_OCEAN_AT_BOTTOM else "seaward"
        ylabel = f"Shoreline displacement since year 0 (m)\n{_up_word} \u25b2"
        lbl_up, lbl_dn = "Accretion (seaward of year 0)", "Erosion (landward of year 0)"
        lbl_ref = f"Year 0 ({START_YEAR}) shoreline"
    else:
        series = pos - np.nanmean(init)
        _up_word = "landward" if GIF_OCEAN_AT_BOTTOM else "seaward"
        ylabel = f"Cross-shore position (m, rel. year-0 mean)\n{_up_word} \u25b2"
        lbl_up, lbl_dn = "Accretion (seaward of year 0)", "Erosion (landward of year 0)"
        lbl_ref = f"Year 0 ({START_YEAR}) shoreline"

    ref = series[0]   # year-0 reference; identically zero except in "position"

    # Fixed axes across all frames so the shoreline doesn't jitter frame-to-frame.
    ymin, ymax = float(np.nanmin(series)), float(np.nanmax(series))
    ypad = (ymax - ymin) * 0.10 or 1.0
    # GIF_OCEAN_AT_BOTTOM reverses the (bottom, top) order passed to ax.set_ylim(),
    # which is matplotlib's documented way to invert an axis - larger/seaward
    # values then land at the bottom instead of the top. The underlying series,
    # shading comparison (cur >= ref), and accretion/erosion colors are untouched;
    # only the visual orientation changes.
    if GIF_OCEAN_AT_BOTTOM:
        ylim = (ymax + ypad, ymin - ypad)
    else:
        ylim = (ymin - ypad, ymax + ypad)

    n_dom   = len(x)
    figsize = (float(np.clip(6.0 + 0.115 * n_dom, 9.0, 18.0)), 5.0)

    frames_dir = os.path.join(run_dir, "gif_frames")
    if keep_frames:
        os.makedirs(frames_dir, exist_ok=True)

    year_idx = list(range(0, n_years, max(int(stride), 1)))
    if year_idx[-1] != n_years - 1:
        year_idx.append(n_years - 1)   # always land on the final year

    be_lbl = "on" if USE_BACKGROUND_EROSION else "off"
    hs_lbl = f"Hs={Hs} m  |  " if Hs is not None else ""

    frames = []
    for t in year_idx:
        cur = series[t]

        fig, ax = plt.subplots(figsize=figsize, dpi=110)
        # Fixed margins (NOT bbox_inches="tight") so every frame is byte-identical
        # in size — mismatched frame dimensions break GIF assembly.
        fig.subplots_adjust(left=0.085, right=0.985, top=0.90, bottom=0.235)
        fig.patch.set_facecolor("white")
        ax.set_facecolor("white")

        # -- Geographic context ------------------------------------------------
        if axis_kind == "gis" and annotate:
            add_geographic_annotations(ax)
        elif axis_kind == "pad":
            ax.axvspan(-0.5, START_REAL_INDEX - 0.5, color="red", alpha=0.10, zorder=0)
            ax.axvspan(END_REAL_INDEX - 0.5, TOTAL_DOMAINS - 0.5, color="red", alpha=0.10, zorder=0)
            ax.axvline(START_REAL_INDEX - 0.5, color="k", ls="--", lw=1.0, alpha=0.5, zorder=2)
            ax.axvline(END_REAL_INDEX - 0.5,   color="k", ls="--", lw=1.0, alpha=0.5, zorder=2)
            if annotate:
                trans_pad = blended_transform_factory(ax.transData, ax.transAxes)
                for span_label, (d_lo, d_hi) in ANN_TOWN_SPANS.items():
                    ax.axvspan(_gis_to_pad(d_lo) - 0.5, _gis_to_pad(d_hi) + 0.5,
                               color=ANN_C_TOWN_SPAN, alpha=0.14, zorder=0)
                    ax.text((_gis_to_pad(d_lo) + _gis_to_pad(d_hi)) / 2.0, 0.90,
                            span_label, transform=trans_pad, ha="center", va="top",
                            fontsize=8, color="0.25", fontweight="bold",
                            bbox=dict(boxstyle="round,pad=0.2", fc="white",
                                      ec="none", alpha=0.85))
                for gname, dom in ANN_GROINS.items():
                    ax.axvline(_gis_to_pad(dom), color=ANN_C_GROIN, lw=1.1, ls=":",
                               alpha=0.85, zorder=2)
                    ax.text(_gis_to_pad(dom), ANN_GROIN_LABEL_Y, gname,
                            transform=trans_pad, ha="center", va="top", fontsize=7,
                            color=ANN_C_GROIN, rotation=90,
                            bbox=dict(boxstyle="round,pad=0.15", fc="white",
                                      ec="none", alpha=0.80))

        # -- Shoreline + erosion/accretion shading -----------------------------
        ax.fill_between(x, cur, ref, where=(cur >= ref), interpolate=True,
                        color="#1565C0", alpha=0.28, zorder=1, label=lbl_up)
        ax.fill_between(x, cur, ref, where=(cur < ref), interpolate=True,
                        color="#B71C1C", alpha=0.22, zorder=1, label=lbl_dn)
        ax.plot(x, ref, color="0.55", ls="--", lw=1.2, zorder=3, label=lbl_ref)
        ax.plot(x, cur, color="#1a2a3a", lw=2.0, zorder=4, label="Shoreline")

        # -- Axes --------------------------------------------------------------
        ax.set_xlim(x_lo - 0.5, x_hi + 0.5)
        ax.set_ylim(*ylim)
        ax.xaxis.set_major_locator(ticker.MultipleLocator(10 if n_dom > 40 else 5))
        ax.xaxis.set_minor_locator(ticker.MultipleLocator(5 if n_dom > 40 else 1))
        ax.tick_params(axis="both", which="major", labelsize=9, direction="in", length=5)
        ax.tick_params(axis="both", which="minor", direction="in", length=3)
        ax.grid(True, which="major", ls=":", lw=0.6, alpha=0.4, color="gray")
        ax.spines[["top", "right"]].set_visible(False)

        if axis_kind == "pad":
            ax.set_xlabel(
                f"CASCADE domain index (buffers included, 0\u2013{TOTAL_DOMAINS - 1})"
                "   |   \u2190 S (Cape Point)      N (Pea Island) \u2192",
                fontsize=11, fontweight="bold")
            top_ax = ax.secondary_xaxis("top")
            tp, tl = [], []
            for gid in range(FIRST_FILE_NUMBER, LAST_FILE_NUMBER + 1, 10):
                tp.append(_gis_to_pad(gid))
                tl.append(str(gid))
            top_ax.set_xticks(tp)
            top_ax.set_xticklabels(tl, fontsize=8)
            top_ax.set_xlabel(f"GIS Domain ID ({FIRST_FILE_NUMBER}\u2013{LAST_FILE_NUMBER})",
                              fontsize=9)
        else:
            ax.set_xlabel(
                "CASCADE Model Domain (500 m alongshore)"
                "   |   \u2190 S (Cape Point)      N (Pea Island) \u2192",
                fontsize=11, fontweight="bold")

        ax.set_ylabel(ylabel, fontsize=10.5, fontweight="bold")

        # -- Year stamp + title + legend ---------------------------------------
        ax.text(0.985, 0.045, f"{START_YEAR + t}  (year {t})",
                transform=ax.transAxes, ha="right", va="bottom",
                fontsize=13, fontweight="bold",
                bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="0.7"))
        ax.set_title(
            (f"Shoreline difference \u2014 {run_name} minus {baseline_label}  |  "
             f"{START_YEAR}\u2013{END_YEAR}  |  {hs_lbl}BE={be_lbl}"
             if mode == "difference" else
             f"Shoreline evolution \u2014 Hatteras Island  |  {START_YEAR}\u2013{END_YEAR}  |  "
             f"{hs_lbl}BE={be_lbl}  |  {run_name}"),
            fontsize=11, fontweight="bold", color="#1a2a3a", pad=(14 if axis_kind == "gis" else 26))
        # Legend outside (below) the axes so it never covers the town / shoal
        # labels, which sit at fixed axes fractions inside the panel.
        ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.115),
                  bbox_transform=ax.transAxes, fontsize=8.5, framealpha=0.95,
                  edgecolor="#cccccc", frameon=True, ncol=4)

        # -- Render to an in-memory PNG (constant size, no disk churn) ----------
        buf = io.BytesIO()
        fig.savefig(buf, format="png", facecolor="white")
        buf.seek(0)
        frames.append(Image.open(buf).convert("RGB"))
        if keep_frames:
            fig.savefig(os.path.join(frames_dir, f"frame_{t:03d}.png"),
                        facecolor="white")
        plt.close(fig)

    sizes = {f.size for f in frames}
    if len(sizes) > 1:
        print(f"  [GIF] frame sizes differ ({sizes}); aborting animation.")
        return None

    gif_path = os.path.join(run_dir, f"{run_name}_shoreline_{mode}_{range_tag}.gif")
    frames[0].save(gif_path, save_all=True, append_images=frames[1:],
                   duration=int(1000.0 / max(float(fps), 0.1)), loop=0, optimize=True)

    print(f"  [GIF] saved shoreline animation ({len(frames)} frames, {mode}, "
          f"{range_tag}): {gif_path}")
    if keep_frames:
        print(f"  [GIF] frames kept in: {frames_dir}")
    if auto_open:
        _open_file(gif_path)
    return gif_path


def make_all_shoreline_gifs(shoreline_m, run_dir, run_name, Hs=None):
    """
    Run every entry in GIF_JOBS for one completed run, plus save the shoreline
    matrix so this run can serve as the baseline for a later difference GIF.

    Each job is independent: one failing (bad range, missing baseline) logs and
    is skipped rather than taking the others — or the run — down with it.
    """
    # Save the raw matrix first: it is the thing a future "difference" job needs,
    # and it is cheap (120 domains x ~40 years). Raw x_s_TS convention, metres,
    # NO sign flip — make_shoreline_gif expects it that way.
    if GIF_SAVE_MATRIX:
        npy_path = os.path.join(run_dir, f"{run_name}_shoreline_matrix.npy")
        np.save(npy_path, np.asarray(shoreline_m, dtype=float))
        print(f"  Saved shoreline matrix: {npy_path}")

    # Load the baseline once, not per job.
    baseline_m = None
    wants_diff = any(str(j.get("mode", "")).lower() == "difference" for j in GIF_JOBS)
    if wants_diff:
        if not GIF_BASELINE_NPY:
            print("  [GIF] GIF_BASELINE_NPY is None; 'difference' jobs will be skipped.")
        elif not os.path.exists(GIF_BASELINE_NPY):
            print(f"  [GIF] baseline not found: {GIF_BASELINE_NPY}; "
                  f"'difference' jobs will be skipped.")
        else:
            baseline_m = np.load(GIF_BASELINE_NPY)
            print(f"  [GIF] baseline loaded ({baseline_m.shape}): {GIF_BASELINE_NPY}")

    made = []
    for job in GIF_JOBS:
        # A "groin" job fans out into one GIF per structure in ANN_GROINS, so
        # adding a second groin needs no change here or in GIF_JOBS. Use
        # "which" to restrict, or range="groin_span" for one enclosing window.
        if str(job.get("range", "")).strip().lower() == "groin":
            try:
                targets = [{"groin": g} for g in _select_groins(job.get("which", "all"))]
            except Exception as e:
                print(f"  [GIF] job {job} skipped ({e}).")
                continue
        else:
            targets = [{"groin": None}]

        for extra in targets:
            try:
                path = make_shoreline_gif(
                    shoreline_m, run_dir, run_name, Hs=Hs,
                    domain_range=job.get("range", "real"),
                    mode=job.get("mode", "displacement"),
                    pad=job.get("pad", 9),
                    groin=extra["groin"],
                    baseline_m=baseline_m,
                    fps=job.get("fps"), stride=job.get("stride"),
                    annotate=job.get("annotate"), auto_open=job.get("auto_open"),
                    keep_frames=job.get("keep_frames"),
                )
                if path:
                    made.append(path)
            except Exception as e:
                label = f"{job} [{extra['groin'][0]}]" if extra["groin"] else str(job)
                print(f"  [GIF] job {label} failed ({e}); continuing.")
    return made


# =============================================================================
# RUN METADATA
# =============================================================================

def write_run_metadata(run_dir, run_name, Hs):
    """
    Write every run parameter to {run_name}_run_metadata.txt - the single
    source of truth for reproducing/comparing a run without reading the
    script itself. Call right after run_dir is created (params are fixed at
    script launch, so timing relative to the cascade run doesn't matter).
    """
    import datetime

    # Collect non-zero background erosion rates (GIS ID -> rate)
    be_nonzero = {
        gis_id: rate
        for gis_id, rate in DOMAIN_BE_RATES.items()
        if rate != 0.0
    }

    lines = [
        "# CASCADE Run Metadata",
        "# Generated automatically - do not edit by hand.",
        "# Re-run the script to regenerate after any parameter change.",
        "",
        "[Run Identity]",
        f"run_name              = {run_name}",
        f"timestamp             = {datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')}",
        f"script                = HAT_hindcast_1984_2024.py",
        f"run_name_suffix       = {RUN_NAME_SUFFIX}",
        "",
        "[Period]",
        f"start_year            = {START_YEAR}",
        f"end_year              = {END_YEAR}",
        f"run_years             = {RUN_YEARS}",
        "",
        "[Wave Climate]",
        f"wave_height_m         = {Hs}",
        f"wave_period_s         = {FIXED_WAVE_PERIOD}",
        f"wave_asymmetry        = {FIXED_WAVE_ASYMMETRY}",
        f"wave_angle_high_frac  = {FIXED_WAVE_ANGLE_HIGH_FRACTION}",
        "",
        "[Sea Level Rise]",
        f"slr_rate_mm_per_yr    = {SEA_LEVEL_RISE_RATE * 1000:.2f}",
        f"slr_constant          = {SEA_LEVEL_CONSTANT}",
        "",
        "[Vertical Datum]",
        f"berm_elevation_m      = {BERM_ELEVATION}   # m NAVD88",
        f"mhw_elevation_m       = {MHW_ELEVATION}    # m NAVD88",
        "",
        "[Storm]",
        f"storm_file            = {os.path.basename(STORM_FILE)}",
        f"storm_file_path       = {STORM_FILE}",
        "",
        "[Topo / Dune Initialisation]",
        f"topo_dune_init_year   = {TOPO_DUNE_INIT_YEAR}",
        f"topo_dune_version     = {TOPO_DUNE_VERSION}",
        f"topo_dune_dir         = {DUNE_TOPO_DIR}",
        f"island_offset_file      = {os.path.basename(ISLAND_OFFSET_FILE)}",
        "",
        "[Morphodynamics]",
        f"dune_rebuild_height_m = {DUNE_REBUILD_HEIGHT}",
        f"rebuild_elev_thresh   = {REBUILD_ELEV_THRESHOLD}  # dam",
        f"overwash_to_dune      = {OVERWASH_TO_DUNE}",
        f"flip_sign_model       = {FLIP_SIGN_MODEL}",
        "",
        "[Background Erosion]",
        f"source_sink_preset    = {SOURCE_SINK_PRESET}",
        f"non_zero_domains      = {len(be_nonzero)}",
    ]

    if be_nonzero:
        for gis_id, rate in sorted(be_nonzero.items()):
            lines.append(f"  GIS_{gis_id:02d}             = {rate:.4f}  # dam/yr")
    else:
        lines.append("  (all domains = 0.0)")

    lines += [
        "",
        "[Overwash Filter]",
        f"default               = {OVERWASH_FILTER_DEFAULT}  (undeveloped + buffers)",
        f"buxton                = {OVERWASH_FILTER_BUXTON}   (GIS 7-8)",
        f"avon                  = {OVERWASH_FILTER_AVON}   (GIS 21-31)",
        f"tri_village           = {OVERWASH_FILTER_SALVO_WAVES_ROD}   (GIS 68-83)",
        "",
        "[Management]",
        f"roadway_management    = {ENABLE_ROADWAY_MANAGEMENT}",
        f"  relocation_blocked_villages = {VILLAGE_GIS_RANGES_NO_RELOCATION}",
        f"  ENABLE_HISTORICAL_ROAD_RELOCATIONS = {ENABLE_HISTORICAL_ROAD_RELOCATIONS}",
        f"  historical_road_events      = "
        f"{[(e['note'], 'enabled' if e.get('enabled', True) else 'DISABLED') for e in HISTORICAL_ROAD_EVENTS]}",
        f"road_setback_file     = {ROAD_SETBACK_FILE}",
        f"road_elevation        = "
        f"{'per-domain' if USE_PER_DOMAIN_ROAD_ELEVATION else 'uniform'}"
        f" ({ROAD_ELEVATION_FALLBACK} m fallback), m MHW-relative",
        f"  road_elevation_file         = {ROAD_ELEVATION_FILE}",
        f"  road_elevation_span_mhw     = "
        f"{road_elevation_full[START_ROAD_INDEX:END_ROAD_INDEX].min():.2f} to "
        f"{road_elevation_full[START_ROAD_INDEX:END_ROAD_INDEX].max():.2f} "
        f"(median {np.median(road_elevation_full[START_ROAD_INDEX:END_ROAD_INDEX]):.2f})",
        f"nourishment           = {ENABLE_NOURISHMENT}",
        f"nourishment_volume    = {NOURISHMENT_VOLUME}  # m3/m (Cascade init default; historical injected per-year)",
        f"nourishment_bn_years  = {HAT_BN_YEARS}",
        f"nourishment_gis_doms  = {sorted(HAT_BN_VOLUME_BY_DOMAIN.keys())}",
        f"sandbag               = {ENABLE_SANDBAG_PLACEMENT}",
        "",
        "[Domain Configuration]",
        f"num_real_domains      = {NUM_REAL_DOMAINS}",
        f"num_buffer_domains    = {NUM_BUFFER_DOMAINS}",
        f"total_domains         = {TOTAL_DOMAINS}",
        f"first_gis_domain      = {FIRST_FILE_NUMBER}",
        f"last_gis_domain       = {LAST_FILE_NUMBER}",
        f"start_real_index      = {START_REAL_INDEX}   # padded",
        f"end_real_index        = {END_REAL_INDEX}  # padded, exclusive",
        "",
        "[Dune Growth]",
        f"rmin                  = 0.55",
        f"rmax                  = 0.95",
        "",
        "[Compute]",
        f"num_cores             = {NUM_CORES}",
    ]

    out_path = os.path.join(run_dir, f"{run_name}_run_metadata.txt")
    with open(out_path, "w", encoding="utf-8") as f:
        f.write("\n".join(lines) + "\n")

    print(f"  Saved run metadata: {out_path}")
    return out_path


# =============================================================================
# CASCADE RUNNER
# =============================================================================
# Note: this function references three module-level globals intentionally -
#   HATTERAS_DATA_BASE : root data directory passed to Cascade as datadir
#   PARAMETER_FILE     : YAML filename resolved by CASCADE relative to datadir
#   OUTPUT_BASE_DIR    : root comparison directory for cascade.save()

def summarise_road_management(cascade, run_dir=None):
    """
    Which domains kept their road, and which lost it — read from CASCADE's own
    post-run state.

    This replaces the domain labelling the old Section 0 monkey-patch provided.
    That version wrapped RoadwayManager.__init__/.update to stamp a label onto
    each drowning print as it happened. It worked, but it meant running a model
    that differed from the package. Everything it told you is recoverable
    afterwards from state CASCADE already exposes, with no patching:

        cascade.road_break[i]              1 once management stopped
        roadways[i].drown_break            road drowned (width or SLR)
        roadways[i].relocation_break       no room left to relocate
        roadways[i]._road_ele_TS           written every managed year, so the
                                           last non-zero entry is the last year
                                           the road was managed

    The one thing it cannot recover is WHY, at the moment it happened, when a
    domain both drowned and was labelled in the same print. drown_break vs
    relocation_break still distinguishes the two causes, which is the part that
    matters for reporting.
    """
    roadways = getattr(cascade, "roadways", None)
    if roadways is None:
        print("\n  [road summary] no roadway managers on this run")
        return []

    # CASCADE constructs a RoadwayManager for EVERY domain -- cascade.py:428,
    # "always initialize just in case we want to add a road ... during the
    # simulation" -- so the existence of one says nothing about whether the road
    # is managed. Only domains whose roadway_management_module flag is True ever
    # have update() called. Village domains (GIS 21-31, 68-83) have a manager
    # object that is never touched, and would otherwise be miscounted here.
    mgmt = getattr(cascade, "roadway_management_module", None)

    rows = []
    for pad, rm in enumerate(roadways):
        if rm is None:
            continue
        gis = FIRST_FILE_NUMBER + (pad - START_REAL_INDEX)
        if not (FIRST_ROAD_DOMAIN <= gis <= LAST_ROAD_DOMAIN):
            continue
        if mgmt is not None and not mgmt[pad]:
            continue    # manager exists but is never updated (e.g. a village)

        ele_ts = np.asarray(getattr(rm, "_road_ele_TS", []), dtype=float)
        nz = np.nonzero(ele_ts)[0]
        last_managed_year = int(nz[-1]) + 1 if nz.size else 0

        drowned = bool(getattr(rm, "drown_break", 0))
        no_room = bool(getattr(rm, "relocation_break", 0))
        rows.append(dict(
            gis=gis, pad=pad,
            drowned=drowned, relocation_blocked=no_room,
            last_managed_year=last_managed_year,
            reason=("drowned" if drowned else
                    "no room to relocate" if no_room else "managed throughout"),
            overwash_removed_m3=float(np.sum(
                getattr(rm, "_road_overwash_volume", [0.0]))),
            dunes_rebuilt=int(np.sum(getattr(rm, "_dunes_rebuilt_TS", [0]))),
            relocations=int(np.sum(getattr(rm, "_road_relocated_TS", [0]))),
        ))

    lost = [r for r in rows if r["drowned"] or r["relocation_blocked"]]
    print("\n" + "=" * 78)
    print("ROAD MANAGEMENT SUMMARY")
    print("=" * 78)
    print(f"  {len(rows)} managed domain(s); {len(lost)} lost their road")
    if lost:
        print(f"\n  {'GIS':>4} {'pad':>4} {'last yr':>8}  reason")
        for r in sorted(lost, key=lambda x: x["last_managed_year"]):
            print(f"  {r['gis']:>4} {r['pad']:>4} {r['last_managed_year']:>8}  "
                  f"{r['reason']}")
        print("\n  Once a road is lost, RoadwayManager.update() returns "
              "immediately every")
        print("  later year: no overwash removal, no dune rebuilding, no "
              "relocation. Those")
        print("  domains are unmanaged barriers for the rest of the run and "
              "must be excluded")
        print("  by name from any managed-vs-unmanaged comparison.")
    else:
        print("  Every managed domain kept its road for the whole run.")

    if run_dir and rows:
        os.makedirs(run_dir, exist_ok=True)
        p = os.path.join(run_dir, "road_management_summary.csv")
        pd.DataFrame(rows).to_csv(p, index=False)
        print(f"\n  [out] {p}")
    return rows


def run_cascade_simulation(
    nt, name, storm_file, alongshore_section_count, num_cores,
    rmin, rmax, elevation_file, dune_file,
    dune_design_elevation, dune_minimum_elevation,
    road_ele, road_width, road_setback,
    historical_road_events,
    overwash_filter, overwash_to_dune,
    nourishment_volume, background_erosion,
    roadway_management_on, beach_dune_manager_on,
    sea_level_rise_rate, sea_level_constant,
    sandbag_management_on, sandbag_elevation,
    enable_shoreline_offset, shoreline_offset,
    wave_height, wave_period, wave_asymmetry, wave_angle_high_fraction,
    berm_elevation, MHW,
    historical_nourishment_on_by_year=None,
    historical_nourishment_volume_by_year=None,
):
    cascade = Cascade(
        HATTERAS_DATA_BASE,
        name,
        storm_file=storm_file,
        elevation_file=elevation_file,
        dune_file=dune_file,
        parameter_file=PARAMETER_FILE,

        berm_elevation=berm_elevation,
        MHW=MHW,

        wave_height=wave_height,
        wave_period=wave_period,
        wave_asymmetry=wave_asymmetry,
        wave_angle_high_fraction=wave_angle_high_fraction,

        sea_level_rise_rate=sea_level_rise_rate,
        sea_level_rise_constant=sea_level_constant,

        background_erosion=background_erosion,
        alongshore_section_count=alongshore_section_count,
        time_step_count=nt,

        min_dune_growth_rate=rmin,
        max_dune_growth_rate=rmax,
        num_cores=num_cores,

        roadway_management_module=roadway_management_on,
        beach_nourishment_module=beach_dune_manager_on,
        sandbag_management_on=sandbag_management_on,
        alongshore_transport_module=True,
        community_economics_module=False,

        road_ele=road_ele,
        road_width=road_width,
        road_setback=road_setback,

        dune_design_elevation=dune_design_elevation,
        dune_minimum_elevation=dune_minimum_elevation,
        sandbag_elevation=sandbag_elevation,

        overwash_filter=overwash_filter,
        overwash_to_dune=overwash_to_dune,

        enable_shoreline_offset=enable_shoreline_offset,
        shoreline_offset=shoreline_offset,

        nourishment_volume=nourishment_volume,
        nourishment_interval=None,
    )

    # ── GROIN (Buxton Groin Field) ──────────────────────────────────────────
    # Inert unless GROIN_ENABLED; see Section 1B. Deferred import so a run
    # with GROIN_ENABLED=False never needs HAT_groin_module.py importable.
    groin_cb = None
    if GROIN_ENABLED:
        try:
            import importlib
            _groin_mod = importlib.import_module(GROIN_MODULE_IMPORT_PATH)
            GroinCallback = _groin_mod.GroinCallback
        except ImportError as e:
            sys.exit(f"ERROR: GROIN_ENABLED=True needs HAT_groin_module.py importable: {e}")

        groin_cb = GroinCallback(
            updrift_pad=_gis_to_pad(GROIN_UPDRIFT_GIS),
            downdrift_pad=_gis_to_pad(GROIN_DOWNDRIFT_GIS),
            trapping_rate_m_yr=GROIN_TRAPPING_RATE_M_YR,
            start_year=START_YEAR,             # active period's start (1984 or 2004)
            install_year=GROIN_INSTALL_YEAR,   # always 1969, regardless of period
            n_domains=alongshore_section_count,
            deterioration_delay_years=GROIN_DETERIORATION_DELAY_YEARS,
            deterioration_mode=GROIN_DETERIORATION_MODE,
            deterioration_ramp_years=GROIN_DETERIORATION_RAMP_YEARS,
            deterioration_fraction=GROIN_DETERIORATION_FRACTION,
        )
        cascade._groin_callback = groin_cb
        print(f"  Groin attached: updrift D{GROIN_UPDRIFT_GIS} "
              f"(pad {_gis_to_pad(GROIN_UPDRIFT_GIS)}), "
              f"downdrift D{GROIN_DOWNDRIFT_GIS} (pad {_gis_to_pad(GROIN_DOWNDRIFT_GIS)}), "
              f"M={GROIN_TRAPPING_RATE_M_YR} m/yr, install {GROIN_INSTALL_YEAR}, "
              f"deterioration @ +{GROIN_DETERIORATION_DELAY_YEARS}yr "
              f"({GROIN_DETERIORATION_MODE}, floor={GROIN_DETERIORATION_FRACTION})")

    historical_nourishment_log = []

    # Tracks padded indices where roadway management has been permanently
    # terminated mid-run (e.g. Jug Handle Bridge 2022 removes GIS 82-88).
    # We track these locally rather than nulling cascade.roadways[pad] because
    # CASCADE's own update() loop calls self._roadways[iB3D].drown_break on every
    # domain unconditionally - setting the object to None causes an AttributeError.
    # By keeping the CASCADE roadway objects intact and skipping affected pads in
    # our management code, the model continues without crashing.
    terminated_road_pads: set = set()

    # Village domains have ROADWAY_MANAGEMENT_ON=False so no RoadwayManager
    # is instantiated for them - relocation is impossible by construction.
    # Non-village road domains run full management (bulldoze + rebuild + relocate).

    for time_step in range(nt - 1):
        current_year = START_YEAR + time_step

        # Reset per-step nourishment flags so nothing carries over from prior year.
        cascade.nourish_now = np.zeros(alongshore_section_count)

        # ── HISTORICAL BEACH NOURISHMENT ────────────────────────────────────
        # If this year appears in the historical schedule, set nourish_now and
        # overwrite the per-domain volume on the corresponding Cascade object.
        if historical_nourishment_on_by_year is not None:
            if current_year in historical_nourishment_on_by_year:
                nourish_now = np.asarray(
                    historical_nourishment_on_by_year[current_year], dtype=float
                )
                nourish_vol = np.asarray(
                    historical_nourishment_volume_by_year[current_year], dtype=float
                )

                if np.any(nourish_now == 1):
                    print(f"\n  -> Applying historical BN in {current_year}:")
                    cascade.nourish_now = nourish_now

                    for iB3D in range(alongshore_section_count):
                        if nourish_now[iB3D] != 1:
                            continue

                        # Update nourishment volume for this domain.
                        # Try the private attribute first (matches Benton/Eve pattern),
                        # fall back to the public attribute if not found.
                        nourishment_obj = cascade.nourishments[iB3D]
                        if hasattr(nourishment_obj, "_nourishment_volume"):
                            nourishment_obj._nourishment_volume = float(nourish_vol[iB3D])
                        elif hasattr(nourishment_obj, "nourishment_volume"):
                            nourishment_obj.nourishment_volume = float(nourish_vol[iB3D])
                        else:
                            raise AttributeError(
                                f"Cannot set nourishment volume on cascade.nourishments[{iB3D}] "
                                f"- no '_nourishment_volume' or 'nourishment_volume' attribute."
                            )

                        gis_id = FIRST_FILE_NUMBER + (iB3D - START_REAL_INDEX)
                        vol_m3_per_m = float(nourish_vol[iB3D])
                        print(
                            f"    GIS {gis_id:3d} (pad {iB3D:3d}): "
                            f"{vol_m3_per_m:.1f} m^3/m  |  "
                            f"{vol_m3_per_m * DOMAIN_LENGTH_M:,.0f} m^3 total"
                        )
                        historical_nourishment_log.append(dict(
                            run_name                  = name,
                            model_year                = current_year,
                            time_step                 = time_step,
                            padded_index              = iB3D,
                            gis_domain                = gis_id,
                            nourishment_volume_m3_per_m = vol_m3_per_m,
                            nourishment_volume_m3_total = vol_m3_per_m * DOMAIN_LENGTH_M,
                        ))

        print(f"\rYear {time_step + 1}/{nt}", end="", flush=True)

        # ── HISTORICAL ROAD MANAGEMENT EVENTS (see Section 8 for event types) ──
        if historical_road_events and hasattr(cascade, "roadways") and cascade.roadways is not None:
            for event in historical_road_events:
                if current_year != event["year"]:
                    continue

                gis_range = range(event["gis_start"], event["gis_end"] + 1)
                pad_indices = [_gis_to_pad(g) for g in gis_range
                               if 0 <= _gis_to_pad(g) < alongshore_section_count]

                if event["type"] == "relocate":
                    if not ENABLE_HISTORICAL_ROAD_RELOCATIONS or not event.get("enabled", True):
                        print(f"\n  -> Road relocation event in {current_year} SKIPPED "
                              f"(disabled by toggle): {event['note']}")
                        continue
                    print(f"\n  -> Road relocation event in {current_year}: "
                          f"{event['note']}")
                    per_domain_disp = event.get("relocation_displacement_m", {})
                    b3d_all = cascade.barrier3d
                    for pad in pad_indices:
                        rm = cascade.roadways[pad] if cascade.roadways is not None else None
                        if rm is None or not hasattr(rm, "_road_setback"):
                            continue
                        gis_id = FIRST_FILE_NUMBER + (pad - START_REAL_INDEX)
                        old_sb = rm._road_setback

                        # DISPLACEMENT, not an absolute setback. rm._road_setback
                        # has been decremented by dune_migrated every year since
                        # t=0, so it already carries the modelled retreat; adding
                        # the measured 1978->1997 offset counts that retreat
                        # exactly once. See the note on HISTORICAL_ROAD_EVENTS.
                        if gis_id in per_domain_disp:
                            new_sb = old_sb + float(per_domain_disp[gis_id])
                        else:
                            new_sb = rm._road_relocation_setback

                        # ---- CASCADE's own relocation guards -----------------
                        # The model-driven relocation path runs both of these.
                        # This prescribed path does NOT, because these events
                        # actually happened and refusing them would be
                        # historically wrong. So evaluate and report, then
                        # apply anyway -- the disagreement goes on the record
                        # rather than being invisible.
                        b3d = b3d_all[pad] if b3d_all is not None else None
                        warnings_out = []
                        if b3d is not None:
                            try:
                                probe_ele, probe_drown = _rm.get_road_relocation_elevation(
                                    time_index=rm._time_index,
                                    xyz_interior_grid=b3d.InteriorDomain,
                                    road_setback=new_sb,
                                    road_width=rm._road_relocation_width,
                                    dx=10, dy=10, dz=10,
                                )
                                if probe_drown:
                                    warnings_out.append(
                                        f"get_road_relocation_elevation REFUSES "
                                        f"(road would be {probe_ele:.2f} m MSL)")
                            except IndexError:
                                warnings_out.append(
                                    "get_road_relocation_elevation: setback is "
                                    "past the end of the interior grid")
                            avg_w = float(b3d.InteriorWidth_AvgTS[-1]) * 10.0
                            need = new_sb + 2 * rm._road_relocation_width
                            if need > avg_w:
                                warnings_out.append(
                                    f"road_relocation_checks REFUSES "
                                    f"(needs {need:.0f} m, island is {avg_w:.0f} m)")

                        # road_ele is DELIBERATELY not touched here.
                        #
                        # A real relocation rebuilds the road at grade on new
                        # ground, so it looks like road_ele should be reset. It
                        # is already correct, and resetting it would be an exact
                        # no-op:
                        #
                        #   road_ele is initialised from RoadElevation.csv,
                        #   measured on the 2004 alignment -- which IS the
                        #   post-relocation road, since both events precede
                        #   2004. roadway_manager then decrements it in the
                        #   Lagrangian frame, so at year t it holds
                        #       measured_2004 - sum(RSLR[0:t]) * 10
                        #   and rebuilding at grade on that same alignment in
                        #   the model's frame gives exactly the same number.
                        #
                        # It would only bite if the post-relocation alignment
                        # had a different measured elevation than the initial
                        # one. It does not, by construction. Do not "fix" this.
                        chk = RELOCATION_CHECK_2004.get(gis_id)
                        chk_txt = (f" | 2004 measured {chk:.0f} m"
                                   if chk is not None else "")
                        rm._road_setback = new_sb
                        rm._road_relocation_setback = new_sb  # keep in sync
                        rm._road_setback_TS[rm._time_index - 1] = new_sb
                        print(f"    GIS {gis_id:3d} (pad {pad:3d}): "
                              f"setback {old_sb:.1f} m + "
                              f"{per_domain_disp.get(gis_id, 0.0):.0f} m "
                              f"displacement -> {new_sb:.1f} m{chk_txt}")
                        for w in warnings_out:
                            print(f"                       [warn] {w}")
                        if warnings_out:
                            print(f"                       -> applying anyway "
                                  f"(prescribed historical event)")

                elif event["type"] == "bridge":
                    print(f"\n  -> Bridge/road removal event in {current_year}: "
                          f"{event['note']}")
                    for pad in pad_indices:
                        rm = cascade.roadways[pad] if cascade.roadways is not None else None
                        if rm is None:
                            continue
                        # DO NOT set cascade.roadways[pad] = None here.
                        # CASCADE's update() loop calls self._roadways[iB3D].drown_break
                        # on every domain with no None check - nulling the object
                        # causes an AttributeError on the very next cascade.update() call.
                        # Instead, record the pad in terminated_road_pads so any future
                        # per-step management code can skip these domains.
                        terminated_road_pads.add(pad)
                        gis_id = FIRST_FILE_NUMBER + (pad - START_REAL_INDEX)
                        print(f"    GIS {gis_id:3d} (pad {pad:3d}): "
                              f"roadway management terminated (bridge/road removed)")
        cascade.update()

        if getattr(cascade, "b3d_break", False):
            print(f"\nModel stopped at year {time_step + 1} (b3d_break)")
            break

    if groin_cb is not None and len(groin_cb.year_TS) == 0:
        print("\n" + "!" * 78)
        print("WARNING: groin callback was never called. The pre-AST hook in")
        print("cascade_groin.py is missing, so this run is identical to no_groin.")
        print("!" * 78)

    save_path = os.path.join(OUTPUT_BASE_DIR, name)
    os.makedirs(save_path, exist_ok=True)
    cascade.save(save_path)
    print(f"\n  Saved: {save_path}")

    # ── Save historical BN log ───────────────────────────────────────────────
    if len(historical_nourishment_log) > 0:
        bn_df = pd.DataFrame(historical_nourishment_log)
        bn_csv = os.path.join(save_path, f"{name}_historical_BN_log.csv")
        bn_df.to_csv(bn_csv, index=False)
        print(f"  Saved BN log ({len(bn_df)} events): {bn_csv}")
    else:
        print(f"  (no historical BN events applied in this run)")

    # ── Save groin diagnostics ────────────────────────────────────────────────
    if groin_cb is not None and len(groin_cb.year_TS) > 0:
        groin_df = pd.DataFrame(groin_cb.diagnostics_frame())
        groin_csv = os.path.join(save_path, f"{name}_groin_diagnostics.csv")
        groin_df.to_csv(groin_csv, index=False)
        print(f"  Saved groin diagnostics ({len(groin_cb.year_TS)} yrs): {groin_csv}")

    return cascade


# =============================================================================
# MAIN
# =============================================================================

def main():
    RMIN = [0.55] * TOTAL_DOMAINS
    RMAX = [0.95] * TOTAL_DOMAINS
    DUNE_DESIGN_ELEVATION  = [DUNE_REBUILD_HEIGHT]    * TOTAL_DOMAINS
    DUNE_MINIMUM_ELEVATION = [REBUILD_ELEV_THRESHOLD] * TOTAL_DOMAINS

    # Per-domain, m MHW-relative, built at module scope from
    # RoadElevation.csv (see PER-DOMAIN ROAD ELEVATION above). CASCADE
    # accepts an array here -- cascade.py:47 keeps it as-is when np.size > 1 and
    # hands self._road_ele[iB3D] to each RoadwayManager. Falls back to a uniform
    # ROAD_ELEVATION_FALLBACK if the file is missing or the toggle is off.
    ROAD_ELEVATION = road_elevation_full
    ROAD_WIDTH     = 20.0

    # ── Per-domain overwash filter ───────────────────────────────────────────
    OVERWASH_FILTER = np.full(TOTAL_DOMAINS, OVERWASH_FILTER_DEFAULT)

    # Buxton: GIS 7-8 -> padded 21-22
    OVERWASH_FILTER[_gis_to_pad(7)  : _gis_to_pad(8)  + 1] = OVERWASH_FILTER_BUXTON
    # Avon: GIS 21-31 -> padded 35-45
    OVERWASH_FILTER[_gis_to_pad(21) : _gis_to_pad(31) + 1] = OVERWASH_FILTER_AVON
    # Salvo/Waves/Rodanthe: GIS 68-83 -> padded 82-97
    OVERWASH_FILTER[_gis_to_pad(68) : _gis_to_pad(83) + 1] = OVERWASH_FILTER_SALVO_WAVES_ROD

    OVERWASH_FILTER = list(OVERWASH_FILTER)

    print("Overwash filter configuration:")
    print(f"  Domains with active filter (> 0): "
          f"{sum(1 for v in OVERWASH_FILTER if v > 0)} of {TOTAL_DOMAINS}")
    for label, gis_start, gis_end, val in [
        ("Buxton",              7,  8,  OVERWASH_FILTER_BUXTON),
        ("Avon",               21, 31,  OVERWASH_FILTER_AVON),
        ("Salvo/Waves/Rod.",   68, 83,  OVERWASH_FILTER_SALVO_WAVES_ROD),
    ]:
        p0, p1 = _gis_to_pad(gis_start), _gis_to_pad(gis_end)
        print(f"  {label:<22s} GIS {gis_start:2d}-{gis_end:2d} "
              f"-> padded {p0:3d}-{p1:3d}  filter={val}")
    print()

    time_span_years = END_YEAR - START_YEAR
    if time_span_years == 0:
        time_span_years = None

    # ── Build historical nourishment arrays ──────────────────────────────────
    # Returns per-year dicts for every year in [START_YEAR, END_YEAR].
    # Period 1 (1984): no events in range -> all-zero arrays returned.
    # Period 2 (2004): injects 2014 Rodanthe and 2022 Avon + Buxton events.
    print("\nBuilding historical nourishment schedule...")
    HIST_NOURISH_ON, HIST_NOURISH_VOL = build_nourishment_arrays_from_manual_inputs()

    # ── Per-domain beach/dune manager flag ───────────────────────────────────
    # Enable only for domains that actually receive nourishment; leave all
    # others False so CASCADE doesn't run the nourishment module unnecessarily.
    # For Period 1 (no events) this collapses to all-False.
    if ENABLE_NOURISHMENT and any(
        np.any(arr == 1) for arr in HIST_NOURISH_ON.values()
    ):
        NOURISHMENT_MANAGEMENT_ON = [False] * TOTAL_DOMAINS
        for gis_id in HAT_BN_VOLUME_BY_DOMAIN:
            pad_idx = _gis_to_pad(gis_id)
            if 0 <= pad_idx < TOTAL_DOMAINS:
                NOURISHMENT_MANAGEMENT_ON[pad_idx] = True
        n_managed = sum(NOURISHMENT_MANAGEMENT_ON)
        print(f"  NOURISHMENT_MANAGEMENT_ON: True for {n_managed} domains "
              f"(GIS {sorted(HAT_BN_VOLUME_BY_DOMAIN.keys())})")
    else:
        # Period 1 or nourishment disabled: all False
        NOURISHMENT_MANAGEMENT_ON = [False] * TOTAL_DOMAINS
        print("  NOURISHMENT_MANAGEMENT_ON: all False (Period 1 or nourishment disabled)")

    # ── Load CoastSat transects + apply LOESS at transect resolution ──────────
    print("\nLoading CoastSat transect data...")
    cs_series = []
    for ds in COASTSAT_DATASETS:
        domain_ids, lrr_values, along_coast_m = load_transect_data(ds)
        if domain_ids is None:
            continue
        windows = []
        for w in LOESS_WINDOW_DOMAINS:
            gis_x, smoothed, frac = loess_smooth_transect_to_domains(
                along_coast_m, lrr_values, domain_ids, w
            )
            if gis_x is None:
                continue
            print(f"  LOESS applied: window={w} domains "
                  f"({w * DOMAIN_SPACING_M / 1000.0:.1f} km)  "
                  f"frac={frac:.3f}  ({ds['label']})")
            windows.append(dict(window=w, gis_x=gis_x, smoothed=smoothed, frac=frac))
        cs_series.append(dict(
            label                = ds["label"],
            period_start         = ds["period_start"],
            active               = (ds["period_start"] == START_YEAR),
            transect_domains     = domain_ids,
            transect_rates       = lrr_values,
            transect_along_coast = along_coast_m,
            windows              = windows,
        ))

    # ── Run CASCADE + compute modeled change rates ───────────────────────────
    rate_profiles = {}

    for Hs in WAVE_HEIGHTS_TO_TEST:
        # Single-value run: name = RUN_NAME_BASE (Hs recorded in metadata only)
        # Multi-value sweep: append _Hs{X} to distinguish runs within the sweep
        if len(WAVE_HEIGHTS_TO_TEST) == 1:
            run_name_hs = RUN_NAME_BASE
        else:
            run_name_hs = f"{RUN_NAME_BASE}_Hs{str(Hs).replace('.', 'p')}"

        cascade = run_cascade_simulation(
            nt=RUN_YEARS,
            name=run_name_hs,
            storm_file=STORM_FILE,
            alongshore_section_count=TOTAL_DOMAINS,
            num_cores=NUM_CORES,

            rmin=RMIN, rmax=RMAX,

            elevation_file=ELEVATION_FILE_PATHS,
            dune_file=DUNE_FILE_PATHS,

            dune_design_elevation=DUNE_DESIGN_ELEVATION,
            dune_minimum_elevation=DUNE_MINIMUM_ELEVATION,

            road_ele=ROAD_ELEVATION,
            road_width=ROAD_WIDTH,
            road_setback=road_setbacks_full,
            historical_road_events=HISTORICAL_ROAD_EVENTS,

            overwash_filter=OVERWASH_FILTER,
            overwash_to_dune=OVERWASH_TO_DUNE,

            nourishment_volume=NOURISHMENT_VOLUME,
            background_erosion=BACKGROUND_EROSION_RATES,

            roadway_management_on=ROADWAY_MANAGEMENT_ON,
            beach_dune_manager_on=NOURISHMENT_MANAGEMENT_ON,

            sea_level_rise_rate=SEA_LEVEL_RISE_RATE,
            sea_level_constant=SEA_LEVEL_CONSTANT,

            sandbag_management_on=SANDBAG_MANAGEMENT_ON,
            sandbag_elevation=SANDBAG_ELEV,

            enable_shoreline_offset=True,
            shoreline_offset=island_offset_dam,

            wave_height=Hs,
            wave_period=FIXED_WAVE_PERIOD,
            wave_asymmetry=FIXED_WAVE_ASYMMETRY,
            wave_angle_high_fraction=FIXED_WAVE_ANGLE_HIGH_FRACTION,
            berm_elevation=BERM_ELEVATION,
            MHW=MHW_ELEVATION,

            historical_nourishment_on_by_year=(
                HIST_NOURISH_ON if ENABLE_NOURISHMENT else None
            ),
            historical_nourishment_volume_by_year=(
                HIST_NOURISH_VOL if ENABLE_NOURISHMENT else None
            ),
        )

        shoreline_m = build_shoreline_matrix(cascade, to_meters=TO_METERS)
        nt_actual, ndom = shoreline_m.shape

        denom = time_span_years if time_span_years is not None else max(nt_actual - 1, 1)
        total_change = shoreline_m[-1, :] - shoreline_m[0, :]
        change_rate  = total_change / float(denom)

        if FLIP_SIGN_MODEL:
            change_rate *= -1.0

        rate_profiles[Hs] = change_rate  # length = TOTAL_DOMAINS (padded)

        # ── Per-run comparison folder ────────────────────────────────────────────
        run_dir = os.path.join(OUTPUT_BASE_DIR, run_name_hs)
        os.makedirs(run_dir, exist_ok=True)

        # Which domains kept their road. Replaces the old Section 0 domain
        # labelling, using CASCADE's own post-run state instead of patching it.
        summarise_road_management(cascade, run_dir=run_dir)

        # ── Save CSV of modeled rates ────────────────────────────────────────
        csv_path = os.path.join(run_dir, f"{run_name_hs}_shoreline_change_rate.csv")
        pd.DataFrame({
            "cascade_padded_index": np.arange(TOTAL_DOMAINS),
            "gis_domain_id": (
                [None] * START_REAL_INDEX
                + list(range(FIRST_FILE_NUMBER, LAST_FILE_NUMBER + 1))
                + [None] * (TOTAL_DOMAINS - END_REAL_INDEX)
            ),
            "model_rate_m_per_yr": change_rate,
        }).to_csv(csv_path, index=False)
        print(f"  Saved rate CSV: {csv_path}")

        # Write full parameter record alongside the model comparison
        write_run_metadata(run_dir, run_name_hs, Hs)

        # ====================================================================
        # PLOT
        # ====================================================================
        # CoastSat overlay styling:
        #   active period (matches START_YEAR) -> solid line, full opacity
        #   other period                       -> dashed line, lighter

        if PLOT_REAL_DOMAINS_ONLY:
            # ── REAL DOMAINS ONLY ────────────────────────────────────────────
            gis_ids = np.arange(FIRST_FILE_NUMBER, LAST_FILE_NUMBER + 1)  # 1..90

            fig, ax = plt.subplots(figsize=(14, 5), constrained_layout=True)

            real_rate = change_rate[START_REAL_INDEX:END_REAL_INDEX]
            ax.plot(gis_ids, real_rate, color=ANN_MODEL_COLOR, linewidth=2,
                    label=f"Model Hs={Hs} m", zorder=6)

            widest_window = max(LOESS_WINDOW_DOMAINS)
            for cs in cs_series:
                is_active = cs["active"]
                if not is_active and not PLOT_REFERENCE_PERIOD:
                    continue
                scatter_x = cs["transect_along_coast"] / DOMAIN_SPACING_M + FIRST_FILE_NUMBER
                if PLOT_RAW_LRR:
                    raw_alpha = RAW_LRR_SCATTER_ALPHA if is_active else RAW_LRR_SCATTER_ALPHA * 0.35
                    ax.scatter(scatter_x, cs["transect_rates"],
                               color=CS_RAW_COLOR, s=RAW_LRR_SCATTER_SIZE,
                               alpha=raw_alpha, zorder=1, linewidths=0)
                widest_win = max(LOESS_WINDOW_DOMAINS)
                for idx, win in enumerate(cs["windows"]):
                    cs_color = CS_WINDOW_COLORS.get(win["window"], CS_WINDOW_COLOR_DEFAULT)
                    lw_base, ls, alpha_factor = (
                        LOESS_WINDOW_STYLES[idx] if idx < len(LOESS_WINDOW_STYLES)
                        else (1.5, "-", 0.80)
                    )
                    is_widest = (win["window"] == widest_win)
                    plot_x, plot_y = splice_loess_with_raw_south(
                        win["gis_x"], win["smoothed"],
                        cs["transect_domains"], cs["transect_rates"],
                        is_widest_window=is_widest,
                    )
                    lbl = f"{cs['label']} — LOESS {win['window']}-dom"
                    if is_active:
                        ax.plot(plot_x, plot_y,
                                color=cs_color, lw=lw_base, ls=ls,
                                alpha=alpha_factor, zorder=4,
                                label=lbl)
                    else:
                        ax.plot(plot_x, plot_y,
                                color=cs_color, lw=lw_base * 0.85, ls=ls,
                                alpha=0.40 * alpha_factor, zorder=3,
                                label=lbl + " (ref)")

            community_spans = [
                (7,  8,  "Buxton"),
                (21, 31, "Avon"),
                (68, 83, "Salvo/Waves/Rod."),
            ]
            for gis_lo, gis_hi, comm_label in community_spans:
                ax.axvspan(gis_lo - 0.5, gis_hi + 0.5,
                           alpha=0.10, color="steelblue", zorder=0)
                ax.text((gis_lo + gis_hi) / 2, ax.get_ylim()[0],
                        comm_label, ha="center", va="bottom",
                        fontsize=7, style="italic", color="steelblue")

            ax.axhline(0.0, linestyle="--", linewidth=1, color="gray", alpha=0.7)

            be_label = "on" if USE_BACKGROUND_EROSION else "off"
            xticks = np.arange(FIRST_FILE_NUMBER, LAST_FILE_NUMBER + 1, DOMAIN_TICK_STEP)
            ax.set_xticks(xticks)
            ax.set_xticklabels([str(i) for i in xticks], rotation=45, ha="right", fontsize=9)

            ax.set_xlabel("GIS Domain ID (1–90)")
            ax.set_ylabel("Shoreline change rate (m/yr)")
            ax.set_title(
                f"Modeled vs CoastSat Shoreline Change Rate – Hatteras Island | "
                f"Real domains only | SLR={SEA_LEVEL_RISE_RATE * 1000:.1f} mm/yr | "
                f"{START_YEAR}–{END_YEAR} | Hs={Hs} m | BE={be_label}"
            )
            ax.grid(alpha=0.3)
            ax.legend()

            fig_suffix = "REAL_DOMAINS_ONLY"

        else:
            # ── ALL DOMAINS (BUFFERS INCLUDED) ───────────────────────────────
            domain_numbers = np.arange(TOTAL_DOMAINS)

            fig, ax = plt.subplots(figsize=(22, 6), constrained_layout=True)

            ax.axvspan(0, START_REAL_INDEX - 0.5, alpha=0.12, color="red", label="Buffer")
            ax.axvspan(END_REAL_INDEX - 0.5, TOTAL_DOMAINS - 1, alpha=0.12, color="red")

            community_spans_pad = [
                (_gis_to_pad(7),  _gis_to_pad(8),  "Buxton"),
                (_gis_to_pad(21), _gis_to_pad(31), "Avon"),
                (_gis_to_pad(68), _gis_to_pad(83), "Salvo/Waves/Rod."),
            ]
            for pad_lo, pad_hi, comm_label in community_spans_pad:
                ax.axvspan(pad_lo - 0.5, pad_hi + 0.5,
                           alpha=0.10, color="steelblue", zorder=0)

            ax.plot(domain_numbers, change_rate, color=ANN_MODEL_COLOR,
                    linewidth=2, label=f"Model Hs={Hs} m", zorder=6)

            for cs in cs_series:
                is_active = cs["active"]
                if not is_active and not PLOT_REFERENCE_PERIOD:
                    continue
                # Convert along-coast distance to fractional padded index
                pad_scatter_x = cs["transect_along_coast"] / DOMAIN_SPACING_M + START_REAL_INDEX
                if PLOT_RAW_LRR:
                    raw_alpha = RAW_LRR_SCATTER_ALPHA if is_active else RAW_LRR_SCATTER_ALPHA * 0.35
                    ax.scatter(pad_scatter_x, cs["transect_rates"],
                               color=CS_RAW_COLOR, s=RAW_LRR_SCATTER_SIZE,
                               alpha=raw_alpha, zorder=1, linewidths=0)
                widest_win = max(LOESS_WINDOW_DOMAINS)
                for idx, win in enumerate(cs["windows"]):
                    cs_color = CS_WINDOW_COLORS.get(win["window"], CS_WINDOW_COLOR_DEFAULT)
                    lw_base, ls, alpha_factor = (
                        LOESS_WINDOW_STYLES[idx] if idx < len(LOESS_WINDOW_STYLES)
                        else (1.5, "-", 0.80)
                    )
                    is_widest = (win["window"] == widest_win)
                    # Splice in GIS space, then convert to padded indices
                    plot_gis_x, plot_y = splice_loess_with_raw_south(
                        win["gis_x"], win["smoothed"],
                        cs["transect_domains"], cs["transect_rates"],
                        is_widest_window=is_widest,
                    )
                    pad_x = plot_gis_x - FIRST_FILE_NUMBER + START_REAL_INDEX
                    lbl   = f"{cs['label']} — LOESS {win['window']}-dom"
                    if is_active:
                        ax.plot(pad_x, plot_y,
                                color=cs_color, lw=lw_base, ls=ls,
                                alpha=alpha_factor, zorder=4, label=lbl)
                    else:
                        ax.plot(pad_x, plot_y,
                                color=cs_color, lw=lw_base * 0.85, ls=ls,
                                alpha=0.40 * alpha_factor, zorder=3,
                                label=lbl + " (ref)")

            ax.axvline(START_REAL_INDEX, linestyle="--", linewidth=1, color="k", alpha=0.5)
            ax.axvline(END_REAL_INDEX - 1, linestyle="--", linewidth=1, color="k", alpha=0.5)
            ax.axhline(0.0, linestyle="--", linewidth=1, color="gray", alpha=0.7)

            y_top = ax.get_ylim()[1] if ax.get_ylim()[1] != 1.0 else 5.0
            ax.text((0 + START_REAL_INDEX) / 2, 0, "Left\nbuffer",
                    ha="center", va="center", fontsize=8, style="italic", alpha=0.6)
            ax.text((START_REAL_INDEX + END_REAL_INDEX) / 2, 0, "Real island (GIS 1–90)",
                    ha="center", va="center", fontsize=9, fontweight="bold", alpha=0.5)
            ax.text((END_REAL_INDEX + TOTAL_DOMAINS) / 2, 0, "Right\nbuffer",
                    ha="center", va="center", fontsize=8, style="italic", alpha=0.6)

            xticks = np.arange(0, TOTAL_DOMAINS, DOMAIN_TICK_STEP)
            ax.set_xticks(xticks)
            ax.set_xticklabels([str(i) for i in xticks], rotation=45, ha="right", fontsize=9)
            ax.set_xlabel("CASCADE domain index (including buffers, 0–119)")

            top_ax = ax.secondary_xaxis("top")
            top_positions = []
            top_labels    = []
            for gis_id in range(FIRST_FILE_NUMBER, LAST_FILE_NUMBER + 1, DOMAIN_TICK_STEP):
                j = START_REAL_INDEX + (gis_id - FIRST_FILE_NUMBER)
                top_positions.append(j)
                top_labels.append(str(gis_id))
            top_ax.set_xticks(top_positions)
            top_ax.set_xticklabels(top_labels, fontsize=9)
            top_ax.set_xlabel(f"GIS Domain ID ({FIRST_FILE_NUMBER}–{LAST_FILE_NUMBER})")

            be_label = "on" if USE_BACKGROUND_EROSION else "off"
            ax.set_ylabel("Shoreline change rate (m/yr)")
            ax.set_title(
                f"Modeled vs CoastSat Shoreline Change Rate – Hatteras Island | "
                f"All domains (buffers included) | SLR={SEA_LEVEL_RISE_RATE * 1000:.1f} mm/yr | "
                f"{START_YEAR}–{END_YEAR} | Hs={Hs} m | BE={be_label}"
            )
            ax.grid(alpha=0.3)
            ax.legend()

            fig_suffix = "ALL_DOMAINS_WITH_BUFFERS"

        fig_out = os.path.join(run_dir, f"{run_name_hs}_shoreline_change_rate_{fig_suffix}.png")
        fig.savefig(fig_out, dpi=300, bbox_inches="tight")
        print(f"  Saved plot: {fig_out}")
        plt.show()

        # ====================================================================
        # ANNOTATED PUBLICATION FIGURE
        # ====================================================================
        # Second figure: modeled rate + both CoastSat periods on GIS domain
        # IDs 1–90, with full geographic annotation layer (communities,
        # village centers, piers, groin, Wimble Shoals influence zone).
        # Always uses the real-domains-only (GIS 1–90) x-axis regardless of
        # PLOT_REAL_DOMAINS_ONLY, since this is a publication/poster figure.

        gis_ids_ann   = np.arange(FIRST_FILE_NUMBER, LAST_FILE_NUMBER + 1)  # 1..90
        real_rate_ann = change_rate[START_REAL_INDEX:END_REAL_INDEX]

        fig2, ax2 = plt.subplots(figsize=(13, 7.0), constrained_layout=True)
        fig2.patch.set_facecolor("white")
        ax2.set_facecolor("white")

        # -- Geographic annotations (drawn first so data renders on top) ------
        add_geographic_annotations(ax2)

        # -- CoastSat scatter + multi-window LOESS lines ----------------------
        data_handles2  = []
        widest_window2 = max(LOESS_WINDOW_DOMAINS)
        for cs in cs_series:
            is_active = cs["active"]
            if not is_active and not PLOT_REFERENCE_PERIOD:
                continue

            # Scatter: respect RAW_LRR_SOUTHERN_ONLY toggle (set in Section 5).
            #   True  -> dots only for D1-LOESS_SKIP_SOUTHERN_DOMAINS (the no-LOESS zone)
            #   False -> dots for all domains D1-90
            scatter_x = cs["transect_along_coast"] / DOMAIN_SPACING_M + FIRST_FILE_NUMBER
            if PLOT_RAW_LRR:
                if RAW_LRR_SOUTHERN_ONLY:
                    south_mask = cs["transect_domains"] <= LOESS_SKIP_SOUTHERN_DOMAINS
                    scatter_x_plot = scatter_x[south_mask]
                    scatter_y_plot = cs["transect_rates"][south_mask]
                    raw_lbl = (f"{cs['label']} — transect LRR (D1-{LOESS_SKIP_SOUTHERN_DOMAINS})"
                               if is_active else None)
                else:
                    scatter_x_plot = scatter_x
                    scatter_y_plot = cs["transect_rates"]
                    raw_lbl = f"{cs['label']} — transect LRR" if is_active else None
                raw_alpha = RAW_LRR_SCATTER_ALPHA if is_active else RAW_LRR_SCATTER_ALPHA * 0.35
                ax2.scatter(scatter_x_plot, scatter_y_plot,
                            color=CS_RAW_COLOR, s=RAW_LRR_SCATTER_SIZE,
                            alpha=raw_alpha, zorder=1, linewidths=0, label=raw_lbl)
                if is_active:
                    data_handles2.append(
                        Line2D([0], [0], color=CS_RAW_COLOR, marker=".", ms=5,
                               ls="none", alpha=RAW_LRR_SCATTER_ALPHA, label=raw_lbl)
                    )
            for idx, win in enumerate(cs["windows"]):
                cs_color = CS_WINDOW_COLORS.get(win["window"], CS_WINDOW_COLOR_DEFAULT)
                lw_base, ls, alpha_factor = (
                    LOESS_WINDOW_STYLES[idx] if idx < len(LOESS_WINDOW_STYLES)
                    else (1.5, "-", 0.80)
                )
                is_widest = (win["window"] == widest_window2)
                gis_x, rate = splice_loess_with_raw_south(
                    win["gis_x"], win["smoothed"],
                    cs["transect_domains"], cs["transect_rates"],
                    is_widest_window=is_widest,
                )
                lbl = f"{cs['label']} — LOESS {win['window']}-dom"
                if is_active:
                    if is_widest:
                        ax2.fill_between(gis_x, rate, 0,
                                         where=(rate < 0),  alpha=0.14, color=cs_color,
                                         interpolate=True)
                        ax2.fill_between(gis_x, rate, 0,
                                         where=(rate >= 0), alpha=0.10, color=cs_color,
                                         interpolate=True)
                    ax2.plot(gis_x, rate, color=cs_color, lw=lw_base,
                             ls=ls, alpha=alpha_factor, zorder=5, label=lbl)
                    data_handles2.append(
                        Line2D([0], [0], color=cs_color, lw=lw_base, ls=ls,
                               alpha=alpha_factor, label=lbl)
                    )
                else:
                    ax2.plot(gis_x, rate, color=cs_color, lw=lw_base * 0.85,
                             ls=ls, alpha=0.40 * alpha_factor, zorder=4)
                    data_handles2.append(
                        Line2D([0], [0], color=cs_color, lw=lw_base * 0.85, ls=ls,
                               alpha=0.40 * alpha_factor, label=lbl + " (ref)")
                    )

        # -- Model line -------------------------------------------------------
        ax2.plot(gis_ids_ann, real_rate_ann,
                 color=ANN_MODEL_COLOR, linewidth=2.4, zorder=6,
                 label=f"Model Hs={Hs} m")
        data_handles2.append(
            Line2D([0], [0], color=ANN_MODEL_COLOR, lw=2.4,
                   label=f"Model Hs={Hs} m")
        )

        # -- Zero reference line ----------------------------------------------
        ax2.axhline(0, color="#2c2c2c", linewidth=1.2, linestyle="--",
                    alpha=0.65, zorder=3)

        # -- Scatter / LOESS transition marker --------------------------------
        # Only drawn when RAW_LRR_SOUTHERN_ONLY=True; marks where raw dots end
        # and LOESS lines begin.
        if PLOT_RAW_LRR and RAW_LRR_SOUTHERN_ONLY and LOESS_SKIP_SOUTHERN_DOMAINS > 0:
            ax2.axvline(LOESS_SKIP_SOUTHERN_DOMAINS + 0.5,
                        color="0.60", lw=0.8, ls=(0, (4, 4)), alpha=0.55, zorder=2)

        # -- Axes formatting --------------------------------------------------
        ax2.set_xlim(FIRST_FILE_NUMBER - 0.5, LAST_FILE_NUMBER + 0.5)
        ax2.xaxis.set_major_locator(ticker.MultipleLocator(10))
        ax2.xaxis.set_minor_locator(ticker.MultipleLocator(5))
        ax2.yaxis.set_major_locator(ticker.MultipleLocator(1))
        ax2.tick_params(axis="both", which="major", labelsize=11,
                        direction="in", length=5)
        ax2.tick_params(axis="both", which="minor", direction="in", length=3)
        ax2.grid(True, which="major", linestyle=":", linewidth=0.6,
                 alpha=0.4, color="gray")
        ax2.spines[["top", "right"]].set_visible(False)
        ax2.spines[["left", "bottom"]].set_linewidth(1.2)

        # Lock ylim to data range, then place accretion / erosion side labels
        all_vals_ann = np.concatenate(
            [real_rate_ann] +
            [w["smoothed"][np.isfinite(w["smoothed"])]
             for cs in cs_series for w in cs["windows"]]
        )
        ymin_ann = all_vals_ann.min()
        ymax_ann = all_vals_ann.max()
        ypad_ann = (ymax_ann - ymin_ann) * 0.06
        ax2.set_ylim(ymin_ann - ypad_ann, ymax_ann + ypad_ann)

        ybot2, ytop2 = ax2.get_ylim()
        zero_frac2   = (0 - ybot2) / (ytop2 - ybot2)
        acc_y = (LABEL_ACCRETION_Y if LABEL_ACCRETION_Y is not None
                 else zero_frac2 + (1 - zero_frac2) / 2)
        ero_y = (LABEL_EROSION_Y   if LABEL_EROSION_Y   is not None
                 else zero_frac2 / 2)
        ax2.text(1.0, acc_y, "Accretion ▲", transform=ax2.transAxes,
                 fontsize=9.5, color="#555555", ha="right", va="center",
                 style="italic")
        ax2.text(1.0, ero_y, "Erosion ▼", transform=ax2.transAxes,
                 fontsize=9.5, color="#555555", ha="right", va="center",
                 style="italic")

        # -- Axis labels & orientation text -----------------------------------
        ax2.set_xlabel("CASCADE Model Domain (500 m alongshore)",
                       fontsize=12, fontweight="bold", labelpad=4)
        ax2.set_ylabel("Shoreline Change Rate (m/yr)",
                       fontsize=12, fontweight="bold", labelpad=8)
        ax2.text(0.0, 1.01, "← S  |  Cape Hatteras",
                 transform=ax2.transAxes, fontsize=9, color="#444444",
                 ha="left", va="bottom", style="italic", clip_on=False)
        ax2.text(1.0, 1.01, "Pea Island  |  N →",
                 transform=ax2.transAxes, fontsize=9, color="#444444",
                 ha="right", va="bottom", style="italic", clip_on=False)

        # -- Title ------------------------------------------------------------
        be_label2 = "on" if USE_BACKGROUND_EROSION else "off"
        ax2.set_title(
            f"Modeled vs Observed Shoreline Change — Hatteras Island, NC  |  "
            f"{START_YEAR}–{END_YEAR}  |  Hs={Hs} m  |  BE={be_label2}",
            fontsize=12, fontweight="bold", pad=12, color="#1a2a3a"
        )

        # -- Legend: data series + annotation proxies -------------------------
        # Placed outside the axes (below) so it never overlaps data or lines.
        ax2.legend(
            handles=data_handles2 + annotation_legend_handles(),
            loc="upper center",
            bbox_to_anchor=(0.5, -0.10),
            bbox_transform=ax2.transAxes,
            fontsize=9.0, framealpha=0.95, edgecolor="#cccccc",
            frameon=True, ncol=5,
        )

        # -- Caption ----------------------------------------------------------
        ax2.annotate(
            f"Model: CASCADE  |  Obs: CoastSat LRR per 500-m domain  |  "
            f"SLR={SEA_LEVEL_RISE_RATE * 1000:.1f} mm/yr  |  Run: {run_name_hs}",
            xy=(0, 0), xycoords="axes fraction",
            xytext=(0, -0.22), textcoords="axes fraction",
            fontsize=7.5, color="#666666", ha="left", va="top", style="italic",
            annotation_clip=False,
        )

        fig2_out = os.path.join(run_dir, f"{run_name_hs}_annotated.png")
        fig2.savefig(fig2_out, dpi=300, bbox_inches="tight", facecolor="white")
        print(f"  Saved annotated plot: {fig2_out}")
        if SHOW_ANNOTATED_FIGURE:
            plt.show()

        # ====================================================================
        # SHORELINE ANIMATION (GIF)
        # ====================================================================
        # Non-fatal: a failed animation must never lose a completed run.
        if MAKE_SHORELINE_GIF:
            try:
                make_all_shoreline_gifs(shoreline_m, run_dir, run_name_hs, Hs=Hs)
            except Exception as e:
                print(f"  [GIF] generation failed ({e}); continuing with the run.")

    # ── end of for Hs loop ───────────────────────────────────────────────────
    if not rate_profiles:
        print("No successful runs; nothing to plot.")
        sys.exit(1)


if __name__ == "__main__":
    main()

