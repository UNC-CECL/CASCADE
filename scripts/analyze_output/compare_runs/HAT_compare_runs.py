#!/usr/bin/env python3
"""
HATTERAS ISLAND — CASCADE Run Comparison Plot
=============================================
Loads up to 4 pre-computed CASCADE runs from their saved shoreline change
rate CSVs and overlays them on one figure for visual comparison.  Adds
LOESS-smoothed CoastSat LRR for reference.

Each run must have already been executed by the hindcast runner, which saves
a rate CSV automatically inside the run directory:

  <run_dir>/{run_name}_shoreline_change_rate.csv
  Columns: gis_domain | change_rate_m_yr | lrr_m_yr | lrr_r2
  lrr_m_yr is the one read -- see RUN_DOMAIN_COL / RUN_RATE_COL.

<run_dir> is resolved from (run_name, period, preset, arm) by
cascade_pipeline.run_registry.find_run_dir -- never joined by hand.

Outputs (saved to COMPARISON_ROOT_DIR/{COMPARISON_NAME}, i.e. under
output/comparisons/):
  {COMPARISON_NAME}_diagnostic.png      — quick multi-run diagnostic
  {COMPARISON_NAME}_annotated.png       — publication figure with geographic annotations
  {COMPARISON_NAME}_residuals.png       — optional panel: each model minus active CoastSat
"""

import os
import pathlib
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import matplotlib.colors as mcolors
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from matplotlib.transforms import blended_transform_factory
from statsmodels.nonparametric.smoothers_lowess import lowess

# =============================================================================
# SECTION 1: DOMAIN CONFIGURATION
# (must match the values used in HAT_hindcast_1984_2024_old version.py)
# =============================================================================

NUM_REAL_DOMAINS   = 90
NUM_BUFFER_DOMAINS = 15

FIRST_FILE_NUMBER = 1     # GIS domain IDs: 1–90
LAST_FILE_NUMBER  = FIRST_FILE_NUMBER + NUM_REAL_DOMAINS - 1   # = 90

TOTAL_DOMAINS    = NUM_BUFFER_DOMAINS + NUM_REAL_DOMAINS + NUM_BUFFER_DOMAINS  # 120
START_REAL_INDEX = NUM_BUFFER_DOMAINS             # = 15
END_REAL_INDEX   = START_REAL_INDEX + NUM_REAL_DOMAINS  # = 105

DOMAIN_TICK_STEP    = 5
DOMAIN_SPACING_M    = 500   # metres per CASCADE domain (used to convert window_domains → km)

# =============================================================================
# SECTION 2: FILE PATHS
# =============================================================================

# ANCHORED, NOT TYPED. These were absolute literals on one machine, and two
# of the three pointed at folders that do not exist: "comparison/raw_runs"
# (the tree is output/raw_runs) and "input_prep/CoastSat" (it is under
# 5-scr/). Anchoring on the pyproject.toml at the repo root makes them follow
# the checkout and survive this file changing depth.
PROJECT_BASE_DIR = next(
    p for p in pathlib.Path(__file__).resolve().parents
    if (p / "pyproject.toml").exists()
)
RAW_RUNS = PROJECT_BASE_DIR / "output" / "raw_runs"

COASTSAT_BASE_DIR = os.path.join(
    PROJECT_BASE_DIR, "scripts", "input_prep", "5-scr", "CoastSat"
)

# Where comparison figures are saved. Products belong under output/, never
# beside the script -- see output/README.md, which names comparisons/ as the
# home for cross-run figures. A subfolder named COMPARISON_NAME is created
# automatically.
COMPARISON_ROOT_DIR = os.path.join(PROJECT_BASE_DIR, "output", "comparisons")

# The rate CSV's schema. WHICH RATE COLUMN IS READ IS A METHOD CHOICE, NOT A
# SPELLING. The file carries two: lrr_m_yr, an OLS slope through the annual
# shoreline positions, and change_rate_m_yr, the endpoint rate. This script
# reads lrr_m_yr because the CoastSat target it is plotted against is one too
# (see rate_col in COASTSAT_DATASETS below) -- putting an endpoint rate up
# against an OLS one moves the difference between two estimators into the
# residual panel, where it reads as model error.
RUN_DOMAIN_COL = "gis_domain"
RUN_RATE_COL   = "lrr_m_yr"

sys.path.insert(0, str(PROJECT_BASE_DIR / "scripts"))

from cascade_pipeline.run_registry import find_run_dir   # noqa: E402

# =============================================================================
# SECTION 3: RUNS TO COMPARE
# =============================================================================
#
# run_name   : folder name AND filename prefix for this run's comparison -
#              must match RUN_NAME_BASE used when the run was executed in
#              the hindcast script. The folder itself is RESOLVED from
#              (run_name, period, preset, arm) by
#              cascade_pipeline.run_registry.find_run_dir; it is never joined
#              by hand here. If run_dir (below) is set, run_name is still
#              used as the filename prefix but the resolver is bypassed.
# period     : "1984_2004" or "2004_2024" - the run's period directory.
# preset     : source/sink preset the run was made under ("calibBE",
#              "edgeBE", "zeroBE"), i.e. the directory below the period.
# arm        : (optional) forcing arm. Omit for the calibration arm, which is
#              where every run made before arms existed sits. Naming no arm
#              never silently searches: a run forced off the calibration wave
#              climate must be asked for by name.
# label      : short display label for the legend.
# start_year : 1984 or 2004. Determines which CoastSat period is drawn solid
#              for this run, AND which panel it appears in for the two-panel
#              period-comparison figure (see plot_two_period below).
# sort_key   : (optional) numeric value used to order this run along the
#              color gradient - typically the swept parameter (Hs, a BE
#              multiplier, etc). If omitted, the run's position in this list
#              is used instead. Lower sort_key -> lighter color, higher -> darker.
# color      : (optional) explicit matplotlib color string. If given, this
#              OVERRIDES the auto-generated gradient color for this run only -
#              useful for pinning one run to a fixed color while letting the
#              rest auto-generate. Leave unset (or None) for normal use.
# run_dir    : (optional) ESCAPE HATCH - full path to the folder containing
#              this run's CSV, for a run that is not in the raw_runs tree at
#              all (another drive, a folder shared with a collaborator, an
#              archive). If given it OVERRIDES the resolver for just this
#              run; other entries still resolve through find_run_dir. Prefer
#              period/preset: a hand-typed path is how a figure ends up
#              drawing a run other than the one it names. The CSV inside must
#              still be named {run_name}_shoreline_change_rate.csv.
#
# Auto-color gradient: each run's color is drawn from RUN_COLORMAP (Section 5)
# at a position determined by its rank along sort_key (or list order). This
# means reordering runs, adding new ones, or changing how many you're
# comparing never requires manually re-picking hex codes.

# THE FOUR RUNS THIS LIST HELD ARE GONE. HAT_1984_2004_SStest,
# HAT_1984_2004_FinalSS, HAT_2004_2024_SStest and HAT_2004_2024_FinalSS were
# addressed by absolute paths into output/raw_runs/source&sink_tests/ and at
# the flat raw_runs/<name> level; neither exists in the tree any more, and
# none of the four names appears in run_index.csv or in
# archived_output_20260828/. The figures they produced are kept at
# output/comparisons/source_sink_zones/ but cannot be regenerated as-is.
# Name live runs below before running this script. (Checked 2026-09-02.)
RUNS_TO_COMPARE = [
    # dict(
    #     run_name   = "HAT_1984_2004_calibBE_road_bdm_groin",
    #     period     = "1984_2004",
    #     preset     = "calibBE",
    #     label      = "Calibrated (1984)",
    #     start_year = 1984,
    #     sort_key   = 2.5,
    # ),
    # dict(
    #     run_name   = "HAT_2004_2024_calibBE_road_bdm_nourish_groin",
    #     period     = "2004_2024",
    #     preset     = "calibBE",
    #     label      = "Calibrated (2004)",
    #     start_year = 2004,
    #     sort_key   = 3.0,
    # ),
    # A run outside the raw_runs tree entirely -- see run_dir in Section 3.
    # dict(
    #     run_name   = "HAT_2004_2024_L7_Hs2p5",
    #     label      = "Hs2.5 (2004)",
    #     start_year = 2004,
    #     sort_key   = 2.5,
    #     run_dir    = r"D:\shared\HAT_2004_2024_L7_Hs2p5",
    # ),
]

# =============================================================================
# >>> NAME THIS ANALYSIS <<<  (controls the comparison folder + filenames)
# =============================================================================
COMPARISON_NAME = "source_sink_zones"   # <-- EDIT THIS to name your comparison folder

# =============================================================================
# SECTION 4: COASTSAT DATASETS
# =============================================================================
# Each entry points to a transect_lrr_full.csv — one row per CoastSat transect
# (~50 m spacing, hundreds of transects total).  LOESS is applied at transect
# resolution then aggregated to domain resolution for the CASCADE comparison.
#
# period_start controls "active vs reference" styling:
#   active period  → full opacity scatter + solid LOESS lines
#   reference period → faded scatter + faded LOESS lines
#
# The active period is inferred from the majority of RUNS_TO_COMPARE start years.
# Override ACTIVE_PERIOD_START manually if your runs span different periods.

COASTSAT_DATASETS = [
    dict(
        label           = "CoastSat (1984–2004)",
        period_start    = 1984,
        csv             = os.path.join(COASTSAT_BASE_DIR, "1984_2004", "transect_lrr_full.csv"),
        domain_col      = "domain_number",
        rate_col        = "lrr_m_yr",
        transect_id_col = "transect_id",
    ),
    dict(
        label           = "CoastSat (2004–2024)",
        period_start    = 2004,
        # *** BUG FIX: this previously pointed to "2004_2024_specific_dates",
        # *** which doesn't match the folder HAT_hindcast_1984_2024_old version.py
        # *** actually writes/reads ("2004_2024"). Since that path didn't
        # *** exist, the load silently failed and cs_series only ever had
        # *** the 1984 entry - every 2004 run was then being compared
        # *** against 1984 CoastSat data in every figure (diagnostic,
        # *** annotated, two-period, residuals), which is why residuals for
        # *** the 2004 runs were reaching +9 m/yr instead of a realistic
        # *** +-1-2 m/yr, and why the right-side 2004 panel had no transect
        # *** scatter at all in the two-period figure.
        csv             = os.path.join(
            COASTSAT_BASE_DIR, "2004_2024", "transect_lrr_full.csv"
        ),
        domain_col      = "domain_number",
        rate_col        = "lrr_m_yr",
        transect_id_col = "transect_id",
    ),
]

# *** NOTE: TRANSECT_DATASETS below is NOT referenced anywhere else in this
# *** file - COASTSAT_DATASETS above already carries transect-level data
# *** (read via the same domain_col/rate_col pattern, see load_transect_data).
# *** This block appears to be left over from an earlier version of the
# *** script's structure. Left here with the same path fix applied so it's
# *** at least not actively wrong if something starts using it again, but
# *** it's currently dead code - safe to delete if you don't need it.
TRANSECT_DATASETS = [
    dict(
        period_start = 1984,
        csv          = os.path.join(COASTSAT_BASE_DIR, "1984_2004", "transect_lrr_full.csv"),
        domain_col   = "domain_number",   # CASCADE domain ID column (1–90)
        lrr_col      = "lrr_m_yr",        # individual transect LRR column
    ),
    dict(
        period_start = 2004,
        csv          = os.path.join(
            COASTSAT_BASE_DIR, "2004_2024", "transect_lrr_full.csv"
        ),
        domain_col   = "domain_number",
        lrr_col      = "lrr_m_yr",
    ),
]

# Which CoastSat period is drawn solid.  Inferred from runs if None.
ACTIVE_PERIOD_START = None   # 1984 or 2004, or None for auto

# --- LOESS smoothing applied to CoastSat overlay ---
# List one or two window sizes (domain units; 1 domain ≈ 500 m).
# Two windows are drawn with distinct line styles so you can compare
# smoothing bandwidth side-by-side.  Set to a single-element list to
# revert to the original single-curve behaviour.
#   10 domains → frac ≈ 0.111 → ~5 km  ← recommended primary
#    7 domains → frac ≈ 0.078 → ~3.5 km ← narrower reference
LOESS_WINDOW_DOMAINS = [7, 10]   # list of 1 or 2 window sizes (domains)

# Styling for each entry in LOESS_WINDOW_DOMAINS (matched by list position).
# Tuple: (linewidth, linestyle, alpha_factor_for_active_period)
# The fill in plot_annotated is drawn only for windows with linestyle "-".
LOESS_WINDOW_STYLES = [
    (1.8, "-",  1.00),   # narrower window (7-dom): solid — distinguished from 10-dom by color
    (2.0, "-",  1.00),   # wider window   (10-dom): solid, full opacity, primary reference
]

# Which LOESS window (domain count) to use as the reference curve in the
# residuals plot.  Must be one of the values in LOESS_WINDOW_DOMAINS.
RESIDUALS_LOESS_WINDOW = 10

# =============================================================================
# SECTION 5: PLOT OPTIONS
# =============================================================================

# Show the residuals figure (model − CoastSat for each run)?
PLOT_RESIDUALS = True

# Show the two-panel period-comparison figure? Left panel = all runs whose
# start_year is 1984, right panel = all runs whose start_year is 2004, each
# with its own active CoastSat overlay, sharing one y-axis and one combined
# legend. Useful when RUNS_TO_COMPARE mixes runs from both periods and you
# want everything side-by-side in a single figure rather than picking one
# ACTIVE_PERIOD_START. Has no effect if all runs share the same start_year
# (the panel for the missing period is simply skipped).
PLOT_TWO_PERIOD = True

# --- Annotation label y-positions (axes fraction: 0.0 = bottom, 1.0 = top) ---
# Pier and groin labels sit on the vertical lines; adjust to avoid overlap.
ANN_PIER_LABEL_Y  = 0.80   # default rotated label y for any pier not given its own override
ANN_PIERS = {
    # *** BUG FIX: label_y values below were 85 and 70 (presumably intended
    # *** as percentages), but the rotated pier labels are drawn with
    # *** blended_transform_factory(ax.transData, ax.transAxes) - meaning y
    # *** is in AXES-FRACTION coordinates (valid range 0.0-1.0), same as
    # *** ANN_PIER_LABEL_Y above. A label_y of 85 placed the label 85x the
    # *** axes height above the plot, which made matplotlib's tight-bbox
    # *** calculation balloon to ~300+ inches on save and crash with an
    # *** out-of-memory error. Fixed to 0.0-1.0 values matching the
    # *** documented default and the equivalent hindcast-script values.
    "Avon Pier":     (26, 0.85),   # (domain, label_y) - label_y is axes-fraction [0,1]
    "Rodanthe Pier": (79, 0.70),
}
ANN_GROIN_LABEL_Y = 0.65

# Accretion / Erosion side labels on the annotated figure.
# Set to None to use the auto-computed midpoint between zero and the plot edge.
# Override with a 0–1 axes fraction to pin the label to a fixed position.
LABEL_ACCRETION_Y = None   # e.g. 0.80 to pin near the top
LABEL_EROSION_Y   = None   # e.g. 0.15 to pin near the bottom

# =============================================================================
# COLOUR PALETTE REFERENCE
# =============================================================================
# CoastSat — cool blue family.  Four layers, light → dark = less → more processed:
#   transect scatter  #9ECAE1   very light blue    individual ~50 m transect LRR (dots)
#   domain avg. line  #9ECAE1   very light blue    domain-averaged LRR (dotted line)
#    7-domain LOESS   #6BAED6   medium sky blue    LOESS narrow window, solid
#   10-domain LOESS   #08519C   deep ocean blue    LOESS primary reference, solid + fill
#
# CASCADE runs — warm orange-red gradient (light → dark = low → high parameter),
# generated automatically from RUN_COLORMAP below rather than hardcoded hex
# values, so any number of runs can be compared without manually assigning
# colors. A run's `color` field (if set) overrides the auto color for just
# that run.
#
# Design rationale:
#   Cool blue = observations (CoastSat)   Warm orange-red = model comparison (CASCADE)
#   Blue + orange-red is the most colorblind-safe pairing (deuteranopia / protanopia)
#   Within each family, lighter → darker encodes low → high smoothing / parameter value
#   Period distinction (1984–2004 vs 2004–2024) is carried by linestyle, not color
# =============================================================================

# Colormap used to auto-generate run colors, sampled light->dark across the
# rank-ordered sort_key (or list order if sort_key is omitted). "YlOrRd" and
# "Oranges" both stay in the warm orange-red family used by prior versions of
# this script; "plasma" or "inferno" work well for >5 runs since they keep
# more contrast between adjacent steps.
RUN_COLORMAP = "YlOrRd"

# The auto gradient is sampled from this range of the colormap (0=lightest
# end, 1=darkest end). Avoiding the very ends keeps the lightest run visible
# against a white background and the darkest run distinguishable from black text.
RUN_COLORMAP_RANGE = (0.35, 0.95)

# CoastSat LOESS colors keyed by window size (domain count).
# Add entries here if you add new window sizes to LOESS_WINDOW_DOMAINS.
CS_WINDOW_COLORS = {
     7: "#6BAED6",   # medium sky blue  — 7-domain LOESS
    10: "#08519C",   # deep ocean blue  — 10-domain LOESS
}
CS_WINDOW_COLOR_DEFAULT = "#4A7C8E"   # fallback for any unlisted window size

# Individual transect LRR scatter — plotted at lowest zorder as context.
# Styling matches HAT_hindcast_1984_2024_old version.py for visual consistency between
# the two scripts.
CS_RAW_COLOR            = "#5BA3C9"    # medium blue
PLOT_RAW_LRR            = True         # set False to hide transect scatter from all figures
RAW_LRR_SOUTHERN_ONLY   = True         # True  -> scatter only D1-LOESS_SKIP_SOUTHERN_DOMAINS
                                        #          (the zone where LOESS is suppressed)
                                        # False -> scatter for all domains D1-90
RAW_LRR_SCATTER_SIZE    = 6            # marker area in points²
RAW_LRR_SCATTER_ALPHA   = 0.60         # opacity for active period; ×0.35 for reference period

# Geographic annotation colors (shared with hindcast script)
ANN_TOWN_SPANS = {
    "Buxton":      (7,   8),
    "Avon":        (21, 31),
    "Tri-Village": (68, 83),
}
ANN_VILLAGE_LINES = {"Salvo": 69, "Waves": 74, "Rodanthe": 80}
ANN_GROINS        = {"Buxton Groin": 5.5}   # boundary between domains 5 and 6
ANN_WIMBLE_SHOALS = (60, 74)
ANN_AVON_SHOALS   = (24, 39)   # Avon Shoals influence zone (same feature type as Wimble Shoals)

ANN_C_TOWN_SPAN    = "#90AFC5"
ANN_C_WIMBLE       = "#E0A800"   # amber - both shoal zones share this color
ANN_C_AVON_SHOALS  = "#E0A800"   # same amber as Wimble Shoals (same feature type)
ANN_C_VILLAGE_LINE = "0.40"
ANN_C_PIER         = "#1565C0"
ANN_C_GROIN        = "#B71C1C"

# Southernmost domains (1 through this value) for which raw per-domain scatter
# is shown instead of LOESS smoothing - Oregon Inlet boundary effects dominate
# this zone and LOESS smoothing there can obscure the sharp gradient. The raw
# scatter toggle below restricts dots to just this zone so domain 11+ is
# represented only by the LOESS lines (matches HAT_hindcast_1984_2024_old version.py).
LOESS_SKIP_SOUTHERN_DOMAINS = 10

# =============================================================================
# HELPER FUNCTIONS — domain utilities
# =============================================================================

def _gis_to_pad(gis_id):
    """Convert a 1-based GIS domain ID to a CASCADE padded array index."""
    return START_REAL_INDEX + (gis_id - FIRST_FILE_NUMBER)


def assign_run_colors(runs):
    """
    Auto-generate a light->dark color gradient for a list of run configs.

    Runs are ranked by `sort_key` if every run provides one; otherwise by
    their position in the list (so simply listing runs low-to-high parameter
    value works without setting sort_key explicitly). Colors are sampled
    evenly across RUN_COLORMAP_RANGE of RUN_COLORMAP, lightest first.

    Any run with an explicit, non-None `color` field keeps that color
    untouched and is excluded from the rank-based assignment - useful for
    pinning one run (e.g. a baseline) to a fixed color while the rest of the
    sweep auto-generates.

    Parameters
    ----------
    runs : list of dict - entries from RUNS_TO_COMPARE (or equivalent)

    Returns
    -------
    list of dict - same entries, each with a resolved 'color' key (hex string)
    """
    cmap = plt.get_cmap(RUN_COLORMAP)
    lo, hi = RUN_COLORMAP_RANGE

    # Split into runs that need an auto color vs. runs with an explicit override
    auto_runs   = [r for r in runs if r.get("color") is None]
    fixed_runs  = [r for r in runs if r.get("color") is not None]

    if auto_runs:
        have_all_keys = all(r.get("sort_key") is not None for r in auto_runs)
        if have_all_keys:
            ranked = sorted(auto_runs, key=lambda r: r["sort_key"])
        else:
            ranked = auto_runs   # fall back to list order

        n = len(ranked)
        for i, r in enumerate(ranked):
            # n==1 -> sample at the dark end so a single run isn't washed out
            frac = (lo + hi) / 2.0 if n == 1 else lo + (hi - lo) * (i / (n - 1))
            r["color"] = mcolors.to_hex(cmap(frac))

    return fixed_runs + auto_runs if fixed_runs else auto_runs


# =============================================================================
# HELPER FUNCTIONS — data loading
# =============================================================================

def load_run_rates(run_name, period=None, preset=None, arm=None, run_dir=None):
    """
    Load the shoreline change rate CSV produced by HAT_hindcast_1984_2024_old version.py.

    Parameters
    ----------
    run_name : str       — used to build the default path AND the expected
                            CSV filename ({run_name}_shoreline_change_rate.csv),
                            which is unchanged regardless of how the folder
                            was resolved.
    period   : str, optional — "1984_2004" or "2004_2024".
    preset   : str, optional — source/sink preset the run was made under.
    arm      : str, optional — forcing arm; None means the calibration arm.
    run_dir  : str, optional — ESCAPE HATCH. Full path to the folder holding
                            that CSV, for a run outside the raw_runs tree.
                            If given it OVERRIDES the resolver entirely.
                            Otherwise period and preset are both required and
                            find_run_dir locates the run, raising with what
                            IS on disk when it is absent.

    Returns
    -------
    gis_ids   : int array — GIS domain IDs, READ FROM the CSV's gis_domain
                            column rather than assumed. Normally 1-90; a
                            short array means the run wrote fewer rows, which
                            is warned about rather than padded over.
    rates_myr : float array, same shape — the run's LRR in m/yr, from
                            RUN_RATE_COL, aligned to gis_ids by construction.
    run_dir   : str                                   — full path to run folder (resolved)
    """
    # RESOLVED, NOT JOINED BY HAND. A run lives at
    # raw_runs/[<arm>/]<period>/<preset>/<run_name>; this function used to
    # join a two-level path that had no slot for the preset or the arm, so
    # every run read as missing. find_run_dir raises naming the arms a run IS
    # under, which reads as "it is over there" rather than "it never existed".
    if run_dir is None:
        if period is None or preset is None:
            raise ValueError(
                f"run '{run_name}' names neither (period, preset) nor an "
                f"explicit run_dir; one of the two is required."
            )
        kwargs = {"arm": arm} if arm else {}
        run_dir = str(find_run_dir(RAW_RUNS, run_name, period, preset, **kwargs))
    csv_path = os.path.join(run_dir, f"{run_name}_shoreline_change_rate.csv")

    if not os.path.exists(csv_path):
        raise FileNotFoundError(
            f"Rate CSV not found for run '{run_name}'.\n"
            f"Expected: {csv_path}\n"
            f"Run HAT_hindcast_1984_2024_old version.py first to generate this file."
        )

    df = pd.read_csv(csv_path)
    required = {RUN_DOMAIN_COL, RUN_RATE_COL}
    if not required.issubset(df.columns):
        raise ValueError(
            f"Rate CSV for '{run_name}' is missing required columns.\n"
            f"Need: {sorted(required)}  |  Found: {sorted(df.columns)}\n"
            f"  {csv_path}"
        )

    # KEYED ON THE FILE'S OWN DOMAIN COLUMN, NOT ON ROW ORDER. This used to
    # slice padded indices 15-104 out of a 120-row CSV and pair them
    # POSITIONALLY with a hand-built arange(1, 91), so the domain ids were
    # asserted rather than read: a run written in the other alongshore order
    # would have shifted every rate against its label with nothing raised.
    # The CSV now carries the 90 real domains keyed on gis_domain, and
    # reading the ids from the file is what makes that failure impossible.
    real = df[df[RUN_DOMAIN_COL].between(FIRST_FILE_NUMBER, LAST_FILE_NUMBER)]
    real = real.sort_values(RUN_DOMAIN_COL)

    if len(real) != NUM_REAL_DOMAINS:
        print(f"  ⚠️  '{run_name}': expected {NUM_REAL_DOMAINS} real-domain rows, "
              f"got {len(real)} — results may be incomplete.")

    gis_ids   = real[RUN_DOMAIN_COL].values.astype(int)
    rates_myr = real[RUN_RATE_COL].values.astype(float)

    return gis_ids, rates_myr, run_dir


def estimate_transect_spacing(along_coast_m):
    """Median spacing between consecutive transects in metres (positive diffs only)."""
    arr   = np.sort(along_coast_m)
    diffs = np.diff(arr)
    pos   = diffs[diffs > 0]
    return float(np.median(pos)) if len(pos) else 50.0


def load_transect_data(ds):
    """
    Load individual transect LRR values from transect_lrr_full.csv and derive
    along-coast distance by spreading each domain's transects evenly across its
    500 m band (mirrors 6-scr-smooth/HAT_loess_method_comparison.py: load_transect_csv).

    Returns
    -------
    domain_ids    : int array   — CASCADE domain ID for each transect
    lrr_values    : float array — LRR (m/yr) for each transect
    along_coast_m : float array — cumulative along-coast distance (m)
    All three arrays share the same length (one entry per transect).
    Returns (None, None, None) on load failure.
    """
    csv_path   = ds["csv"]
    domain_col = ds["domain_col"]
    rate_col   = ds["rate_col"]
    id_col     = ds.get("transect_id_col", "transect_id")

    if not os.path.exists(csv_path):
        print(f"  ⚠️  Transect CSV not found: {csv_path}")
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

    # Spread each domain's transects evenly across its 500 m band so that
    # physical spacing can be estimated correctly for the LOESS frac.
    def _spread(grp):
        n         = len(grp)
        base      = (grp[domain_col].iloc[0] - 1) * DOMAIN_SPACING_M
        offsets   = (np.arange(n) + 0.5) * (DOMAIN_SPACING_M / n)
        grp       = grp.copy()
        grp["along_coast_m"] = base + offsets
        return grp

    df = df.groupby(domain_col, group_keys=False).apply(_spread)

    domain_ids    = df[domain_col].values.astype(int)
    lrr_values    = df[rate_col].values.astype(float)
    along_coast_m = df["along_coast_m"].values.astype(float)

    spacing = estimate_transect_spacing(along_coast_m)
    print(f"  ✓ {ds['label']}: {len(df)} transects  "
          f"est. spacing {spacing:.0f} m  "
          f"LRR range {np.nanmin(lrr_values):+.2f}–{np.nanmax(lrr_values):+.2f} m/yr")
    return domain_ids, lrr_values, along_coast_m


def loess_smooth_transect_to_domains(along_coast_m, lrr, domain_ids, window_domains):
    """
    Apply LOESS at transect resolution using physical along-coast distance (m) as x,
    then aggregate smoothed values to CASCADE domain resolution by averaging within
    each domain.  Mirrors smooth_transect_df() + aggregate_to_domains() from
    6-scr-smooth/HAT_loess_method_comparison.py.

    window_domains is converted to km (× DOMAIN_SPACING_M) so the physical window
    is consistent regardless of transect density.

    Returns
    -------
    gis_x    : int array   — domain IDs that have at least one transect
    smoothed : float array — domain-averaged smoothed LRR (m/yr), same length as gis_x
    frac     : float       — LOESS frac used (for logging)
    """
    window_km = window_domains * DOMAIN_SPACING_M / 1000.0
    spacing_m = estimate_transect_spacing(along_coast_m)
    n         = len(along_coast_m)
    frac      = float(np.clip((window_km * 1000.0 / spacing_m) / n, 0.02, 1.0))

    valid = np.isfinite(lrr)
    if valid.sum() < 5:
        print(f"  ⚠️  Too few valid transects ({valid.sum()}) for LOESS — skipping")
        return None, None, frac

    result            = lowess(lrr[valid], along_coast_m[valid], frac=frac, return_sorted=True)
    smoothed_t        = np.full(n, np.nan)
    smoothed_t[valid] = np.interp(along_coast_m[valid], result[:, 0], result[:, 1])

    # Average smoothed transect values within each CASCADE domain
    dom_agg = (pd.DataFrame({"domain": domain_ids, "smoothed": smoothed_t})
                 .groupby("domain")["smoothed"].mean()
                 .dropna())

    return dom_agg.index.values.astype(int), dom_agg.values, frac


def splice_loess_with_raw_south(win_gis_x, win_smoothed, skip_n=None):
    """
    Trim a LOESS curve so it starts north of the southernmost `skip_n`
    domains, leaving that southern zone to show raw transect scatter only.

    Ported from HAT_hindcast_1984_2024_old version.py's function of the same name for
    visual consistency between the two scripts - domains 1-skip_n are
    boundary-affected (Oregon Inlet dynamics) and LOESS smoothing there can
    obscure the sharp gradient rather than clarify it, so the LOESS line is
    simply not drawn there; the raw scatter (already restricted to this same
    zone via RAW_LRR_SOUTHERN_ONLY) carries the signal instead.

    Parameters
    ----------
    win_gis_x    : int array   - GIS domain IDs from the LOESS result
    win_smoothed : float array - LOESS-smoothed LRR (m/yr)
    skip_n       : int         - domains 1..skip_n are excluded from the
                                  returned line. Defaults to
                                  LOESS_SKIP_SOUTHERN_DOMAINS.

    Returns
    -------
    plot_x, plot_y : arrays - the LOESS curve restricted to domains > skip_n
    """
    if skip_n is None:
        skip_n = LOESS_SKIP_SOUTHERN_DOMAINS
    if skip_n <= 0:
        return win_gis_x, win_smoothed
    mask = win_gis_x > skip_n
    return win_gis_x[mask], win_smoothed[mask]


# =============================================================================
# HELPER FUNCTIONS — geographic annotations
# =============================================================================

def add_geographic_annotations(ax):
    """
    Draw the standard Hatteras geographic annotation layer onto ax.
    X-axis must be in GIS domain IDs (1–90).
    """
    trans = blended_transform_factory(ax.transData, ax.transAxes)

    # 1a. Avon Shoals influence zone (same feature type as Wimble Shoals)
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

    # 2. Community spans
    for span_label, (d_lo, d_hi) in ANN_TOWN_SPANS.items():
        ax.axvspan(d_lo - 0.5, d_hi + 0.5,
                   color=ANN_C_TOWN_SPAN, alpha=0.14, zorder=0)
        ax.text((d_lo + d_hi) / 2.0, 0.90,
                span_label, transform=trans,
                ha="center", va="top", fontsize=8, color="0.25", fontweight="bold",
                bbox=dict(boxstyle="round,pad=0.2", fc="white", ec="none", alpha=0.85))

    # 3. Village center lines
    for vname, dom in ANN_VILLAGE_LINES.items():
        ax.axvline(dom, color=ANN_C_VILLAGE_LINE, lw=0.9, ls="--", alpha=0.65, zorder=1)
        ax.text(dom, 0.84, vname, transform=trans,
                ha="center", va="top", fontsize=7.5, color="0.30",
                bbox=dict(boxstyle="round,pad=0.15", fc="white", ec="none", alpha=0.80))

    # 4. Pier lines
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
        Line2D([0], [0], color=ANN_C_PIER,         lw=1.0, ls="-.", label="Pier"),
        Line2D([0], [0], color=ANN_C_GROIN,        lw=1.1, ls=":",  label="Groin"),
    ]


def _style_ax(ax, ylabel="Shoreline change rate (m/yr)"):
    """Apply shared axis styling."""
    ax.set_xlim(FIRST_FILE_NUMBER - 0.5, LAST_FILE_NUMBER + 0.5)
    ax.axhline(0.0, color="#2c2c2c", linewidth=1.0, linestyle="--", alpha=0.55, zorder=3)
    ax.xaxis.set_major_locator(ticker.MultipleLocator(10))
    ax.xaxis.set_minor_locator(ticker.MultipleLocator(5))
    ax.tick_params(axis="both", which="major", labelsize=10, direction="in", length=5)
    ax.tick_params(axis="both", which="minor", direction="in", length=3)
    ax.grid(True, which="major", linestyle=":", linewidth=0.6, alpha=0.35, color="gray")
    ax.spines[["top", "right"]].set_visible(False)
    ax.spines[["left", "bottom"]].set_linewidth(1.1)
    ax.set_ylabel(ylabel, fontsize=11, fontweight="bold", labelpad=8)


# =============================================================================
# PLOTTING
# =============================================================================

def plot_diagnostic(run_data, cs_series, active_period, out_path, comparison_name):
    """
    Quick diagnostic plot — all runs + CoastSat + geographic annotations on
    one panel. Uses the same drawing helper as plot_annotated/plot_two_period
    so all comparison figures stay visually consistent.

    Previously this plot had NO geographic annotations and used
    loc="best" for the legend, which with 4+ runs landed on top of the
    data (see uploaded screenshot) - both fixed below.
    """
    fig, ax = plt.subplots(figsize=(15, 7))
    fig.patch.set_facecolor("white")
    ax.set_facecolor("white")

    model_handles, cs_handles = _draw_comparison_panel(ax, run_data, cs_series, active_period)

    ax.set_xlabel("CASCADE Model Domain (500 m alongshore)",
                  fontsize=11, fontweight="bold", labelpad=4)
    ax.set_title(
        f"CASCADE Run Comparison — Hatteras Island, NC  |  {comparison_name}",
        fontsize=12, fontweight="bold", pad=12, color="#1a2a3a"
    )

    # Reserve bottom margin BEFORE placing the legend so it stays in
    # on-canvas figure-fraction coordinates (see plot_annotated for the
    # full explanation of why this matters - a previous off-canvas
    # bbox_to_anchor caused an out-of-memory crash on save in this script).
    fig.subplots_adjust(bottom=0.26)

    all_handles = model_handles + cs_handles + annotation_legend_handles()
    fig.legend(handles=all_handles,
               loc="lower center",
               bbox_to_anchor=(0.5, 0.04),
               fontsize=9, framealpha=0.95, edgecolor="#cccccc",
               frameon=True, ncol=4)

    fig.savefig(out_path, dpi=300, facecolor="white")
    plt.close(fig)
    print(f"✓ Saved diagnostic:  {os.path.basename(out_path)}")


def _draw_comparison_panel(ax, run_data, cs_series, active_period,
                            panel_title=None, show_reference_period=True):
    """
    Draw geographic annotations + CoastSat curves + model run lines onto a
    single axis. Shared by plot_annotated (one panel) and plot_two_period
    (two panels, called once per period) so both stay visually identical.

    Parameters
    ----------
    show_reference_period : bool
        True  -> the inactive CoastSat period is still drawn, faded, with a
                 "(ref)" label (useful in plot_annotated's single-panel view,
                 where showing the other period for context is informative).
        False -> only the active period (cs["period_start"] == active_period)
                 is drawn at all; the inactive period is skipped entirely.
                 This is what plot_two_period uses, since each panel is
                 already dedicated to one period - showing the other period
                 there just duplicates what the OTHER panel is for.

    Returns
    -------
    model_handles, cs_handles : lists of Line2D proxies for the legend
    """
    add_geographic_annotations(ax)

    cs_handles = []
    widest_window = max(LOESS_WINDOW_DOMAINS)   # fill drawn only for the widest window
    for cs in cs_series:
        is_active = cs["period_start"] == active_period
        if not is_active and not show_reference_period:
            continue   # skip the inactive period entirely for this panel
        # --- individual transect scatter (lowest zorder — contextual background) ---
        scatter_x = cs["transect_along_coast"] / DOMAIN_SPACING_M + FIRST_FILE_NUMBER
        if PLOT_RAW_LRR:
            if RAW_LRR_SOUTHERN_ONLY:
                south_mask     = cs["transect_domains"] <= LOESS_SKIP_SOUTHERN_DOMAINS
                scatter_x_plot = scatter_x[south_mask]
                scatter_y_plot = cs["transect_rates"][south_mask]
                raw_lbl = (f"{cs['label']} — transect LRR (D1-{LOESS_SKIP_SOUTHERN_DOMAINS})"
                           if is_active else None)
            else:
                scatter_x_plot = scatter_x
                scatter_y_plot = cs["transect_rates"]
                raw_lbl = f"{cs['label']} — transect LRR" if is_active else None
            raw_alpha = RAW_LRR_SCATTER_ALPHA if is_active else RAW_LRR_SCATTER_ALPHA * 0.35
            ax.scatter(scatter_x_plot, scatter_y_plot,
                       color=CS_RAW_COLOR, s=RAW_LRR_SCATTER_SIZE,
                       alpha=raw_alpha, zorder=1, linewidths=0, label=raw_lbl)
            if is_active:
                cs_handles.append(
                    Line2D([0], [0], color=CS_RAW_COLOR, marker=".", ms=5,
                           ls="none", alpha=RAW_LRR_SCATTER_ALPHA, label=raw_lbl)
                )
        # --- LOESS smoothed curves (transect-based, aggregated to domain resolution) ---
        for idx, win in enumerate(cs["windows"]):
            cs_color  = CS_WINDOW_COLORS.get(win["window"], CS_WINDOW_COLOR_DEFAULT)
            lw_base, ls, alpha_factor = (
                LOESS_WINDOW_STYLES[idx] if idx < len(LOESS_WINDOW_STYLES)
                else (1.5, "--", 0.80)
            )
            # Trim the LOESS line to start north of LOESS_SKIP_SOUTHERN_DOMAINS
            # when RAW_LRR_SOUTHERN_ONLY is on, leaving D1-10 to show only the
            # raw transect scatter (matches HAT_hindcast_1984_2024_old version.py styling).
            if RAW_LRR_SOUTHERN_ONLY:
                w_gis_x, rate = splice_loess_with_raw_south(win["gis_x"], win["smoothed"])
            else:
                w_gis_x = win["gis_x"]
                rate     = win["smoothed"]
            w_lbl    = f"LOESS {win['window']}-dom"
            lbl      = f"{cs['label']} — {w_lbl}"
            if is_active:
                if win["window"] == widest_window:
                    ax.fill_between(w_gis_x, rate, 0,
                                    where=(rate < 0),  alpha=0.14, color=cs_color,
                                    interpolate=True)
                    ax.fill_between(w_gis_x, rate, 0,
                                    where=(rate >= 0), alpha=0.10, color=cs_color,
                                    interpolate=True)
                ax.plot(w_gis_x, rate, color=cs_color, linewidth=lw_base,
                        linestyle=ls, alpha=alpha_factor, zorder=4, label=lbl)
                cs_handles.append(
                    Line2D([0], [0], color=cs_color, lw=lw_base, ls=ls,
                           alpha=alpha_factor, label=lbl)
                )
            else:
                ax.plot(w_gis_x, rate, color=cs_color, linewidth=lw_base * 0.85,
                        linestyle=ls, alpha=0.40 * alpha_factor, zorder=3)
                cs_handles.append(
                    Line2D([0], [0], color=cs_color, lw=lw_base * 0.85, ls=ls,
                           alpha=0.40 * alpha_factor, label=lbl + " (ref)")
                )

    model_handles = []
    for run in run_data:
        ax.plot(run["gis_ids"], run["rates"],
                color=run["color"], linewidth=2.4, zorder=5, label=run["label"])
        model_handles.append(
            Line2D([0], [0], color=run["color"], lw=2.4, label=run["label"])
        )

    _style_ax(ax)

    # Compass / orientation labels
    ax.text(0.0, 1.01, "← S  |  Cape Hatteras",
            transform=ax.transAxes, fontsize=9, color="#444444",
            ha="left", va="bottom", style="italic", clip_on=False)
    ax.text(1.0, 1.01, "Pea Island  |  N →",
            transform=ax.transAxes, fontsize=9, color="#444444",
            ha="right", va="bottom", style="italic", clip_on=False)

    if panel_title:
        ax.set_title(panel_title, fontsize=11, fontweight="bold", pad=10, color="#1a2a3a")

    return model_handles, cs_handles


def plot_annotated(run_data, cs_series, active_period, out_path, comparison_name):
    """
    Publication-quality figure with full geographic annotation layer.
    """
    fig, ax = plt.subplots(figsize=(14, 7.5))
    fig.patch.set_facecolor("white")
    ax.set_facecolor("white")

    model_handles, cs_handles = _draw_comparison_panel(ax, run_data, cs_series, active_period)

    ax.set_xlabel("CASCADE Model Domain (500 m alongshore)",
                  fontsize=12, fontweight="bold", labelpad=4)

    # Accretion / erosion side labels
    all_vals = np.concatenate(
        [r["rates"] for r in run_data] +
        [w["smoothed"][np.isfinite(w["smoothed"])]
         for cs in cs_series for w in cs["windows"]]
    )
    ymin = all_vals.min(); ymax = all_vals.max()
    ypad = (ymax - ymin) * 0.07
    ax.set_ylim(ymin - ypad, ymax + ypad)
    ybot, ytop = ax.get_ylim()
    zero_frac  = (0 - ybot) / (ytop - ybot)
    acc_y = LABEL_ACCRETION_Y if LABEL_ACCRETION_Y is not None else zero_frac + (1 - zero_frac) / 2
    ero_y = LABEL_EROSION_Y   if LABEL_EROSION_Y   is not None else zero_frac / 2
    ax.text(1.0, acc_y, "Accretion ▲",
            transform=ax.transAxes, fontsize=9, color="#555555",
            ha="right", va="center", style="italic")
    ax.text(1.0, ero_y, "Erosion ▼",
            transform=ax.transAxes, fontsize=9, color="#555555",
            ha="right", va="center", style="italic")

    ax.set_title(
        f"CASCADE Run Comparison — Hatteras Island, NC  |  {comparison_name}",
        fontsize=12, fontweight="bold", pad=12, color="#1a2a3a"
    )

    # Reserve bottom margin BEFORE placing the legend/caption so both stay
    # in on-canvas figure-fraction coordinates ([0, 1]). A prior version
    # placed these via ax.legend(bbox_to_anchor=(...), bbox_transform=
    # ax.transAxes) with a negative y - in this matplotlib version that
    # combination produced a degenerate tight-bbox computation (~490 inches
    # tall) on save, causing an out-of-memory crash. fig.legend/fig.text in
    # plain figure-fraction coordinates with a reserved margin avoids that
    # failure mode entirely.
    fig.subplots_adjust(bottom=0.26)

    # Legend: model runs | CoastSat | geographic annotations
    all_handles = model_handles + cs_handles + annotation_legend_handles()
    fig.legend(handles=all_handles,
               loc="lower center",
               bbox_to_anchor=(0.5, 0.07),
               fontsize=9, framealpha=0.95, edgecolor="#cccccc",
               frameon=True, ncol=4)

    # Caption — figure-fraction, sits below the legend, fully on-canvas.
    fig.text(
        0.5, 0.01,
        f"Model: CASCADE  |  Observed: CoastSat LRR "
        f"(LOESS {'/'.join(str(w) for w in LOESS_WINDOW_DOMAINS)}-domain windows)  |  "
        f"Comparison: {comparison_name}",
        fontsize=7.5, color="#666666", ha="center", va="bottom", style="italic",
    )

    fig.savefig(out_path, dpi=300, facecolor="white")
    plt.close(fig)
    print(f"✓ Saved annotated:   {os.path.basename(out_path)}")


def plot_two_period(run_data, cs_series, out_path, comparison_name):
    """
    Two-panel figure: left = all runs with start_year=1984 (vs. 1984-2004
    CoastSat), right = all runs with start_year=2004 (vs. 2004-2024 CoastSat).
    Shares one y-axis range and one combined legend below both panels.

    If every run shares the same start_year, the empty panel is skipped and
    a single-panel figure is produced instead (so this is always safe to call
    regardless of what's in RUNS_TO_COMPARE).
    """
    runs_1984 = [r for r in run_data if r["start_year"] == 1984]
    runs_2004 = [r for r in run_data if r["start_year"] == 2004]

    panels = []
    if runs_1984:
        panels.append((1984, runs_1984, "1984-2004"))
    if runs_2004:
        panels.append((2004, runs_2004, "2004-2024"))

    if not panels:
        print("  ⚠️  No runs with a recognized start_year — skipping two-period figure.")
        return

    n_panels = len(panels)
    fig, axes = plt.subplots(
        1, n_panels, figsize=(11 * n_panels, 8), sharey=True,
    )
    if n_panels == 1:
        axes = [axes]
    fig.patch.set_facecolor("white")

    all_model_handles = []
    all_cs_handles_by_period = {}
    all_vals_chunks = []

    for ax, (period_start, panel_runs, period_label) in zip(axes, panels):
        ax.set_facecolor("white")
        model_handles, cs_handles = _draw_comparison_panel(
            ax, panel_runs, cs_series, period_start,
            panel_title=f"{period_label}",
            show_reference_period=False,   # each panel shows ONLY its own period's CoastSat
        )
        all_model_handles.extend(model_handles)
        # Keep cs_handles per-period so duplicate "active" labels across
        # panels (e.g. both panels showing their own active LOESS curve)
        # don't collide - keyed by period so panel 2's handles don't get
        # silently dropped as "duplicates" of panel 1's.
        all_cs_handles_by_period[period_start] = cs_handles
        all_vals_chunks.append(np.concatenate(
            [r["rates"] for r in panel_runs] +
            [w["smoothed"][np.isfinite(w["smoothed"])]
             for cs in cs_series if cs["period_start"] == period_start
             for w in cs["windows"]]
        ) if panel_runs else np.array([]))
        ax.set_xlabel("CASCADE Model Domain (500 m alongshore)",
                      fontsize=11, fontweight="bold", labelpad=4)

    # Shared y-limits across both panels, computed from ALL data in either panel
    all_vals = np.concatenate([c for c in all_vals_chunks if c.size])
    ymin = all_vals.min(); ymax = all_vals.max()
    ypad = (ymax - ymin) * 0.07
    for ax in axes:
        ax.set_ylim(ymin - ypad, ymax + ypad)
        ybot, ytop = ax.get_ylim()
        zero_frac  = (0 - ybot) / (ytop - ybot)
        acc_y = LABEL_ACCRETION_Y if LABEL_ACCRETION_Y is not None else zero_frac + (1 - zero_frac) / 2
        ero_y = LABEL_EROSION_Y   if LABEL_EROSION_Y   is not None else zero_frac / 2
        ax.text(1.0, acc_y, "Accretion ▲",
                transform=ax.transAxes, fontsize=8.5, color="#555555",
                ha="right", va="center", style="italic")
        ax.text(1.0, ero_y, "Erosion ▼",
                transform=ax.transAxes, fontsize=8.5, color="#555555",
                ha="right", va="center", style="italic")

    # Only the leftmost panel gets the y-axis label (sharey hides the rest)
    axes[0].set_ylabel("Shoreline change rate (m/yr)", fontsize=11,
                        fontweight="bold", labelpad=8)

    # Reserve real margin at the top (suptitle) and bottom (legend + caption)
    # BEFORE placing those artists, so every coordinate below stays inside
    # the canvas (figure-fraction range [0, 1]). Placing fig.legend/fig.text
    # at negative y (i.e. below the canvas) forces bbox_inches="tight" to
    # expand the saved bounding box to include that off-canvas content,
    # which on a wide multi-panel figure produced a runaway canvas size and
    # an out-of-memory crash on save - this keeps everything on-canvas instead.
    # bottom=0.20 leaves room for: xlabel (~0.04) + legend (now deduplicated
    # by visual style across the two panels - typically ~12 handles at
    # ncol=4 -> 3 rows, ~0.10) + caption (~0.02). Previously this was 0.30,
    # sized for the legend BEFORE the cross-panel CoastSat dedup fix (which
    # used to list every CoastSat entry twice, once per period); the smaller
    # deduplicated legend left a large unused gap with the original margin.
    fig.subplots_adjust(top=0.88, bottom=0.20, wspace=0.06)

    fig.suptitle(
        f"CASCADE Run Comparison by Period — Hatteras Island, NC  |  {comparison_name}",
        fontsize=13, fontweight="bold", color="#1a2a3a", y=0.97,
    )

    # Combined legend: model runs from both panels (deduplicated by label,
    # since the same run could theoretically appear via panel_runs once each)
    # + one CoastSat entry set per period shown + annotation proxies.
    seen_labels = set()
    dedup_model_handles = []
    for h in all_model_handles:
        if h.get_label() not in seen_labels:
            dedup_model_handles.append(h)
            seen_labels.add(h.get_label())

    combined_cs_handles = []
    for period_start, _, period_label in panels:
        combined_cs_handles.extend(all_cs_handles_by_period.get(period_start, []))

    # De-duplicate CoastSat legend entries by VISUAL STYLE rather than label
    # text. Both panels draw their own period's CoastSat in the same colors
    # (CS_WINDOW_COLORS is keyed by window size only, not by period; CS_RAW_COLOR
    # is a single global color) - so "CoastSat (1984-2004) - LOESS 7-dom" and
    # "CoastSat (2004-2024) - LOESS 7-dom" render as the literal same line
    # style. Previously the legend listed both anyway since they were
    # deduplicated by exact label text, which never matched (different period
    # in the label). Since each panel's title already says which period it
    # is, one shared legend entry per visual style is enough here - relabeled
    # to be period-agnostic ("LOESS 7-dom" instead of "CoastSat (1984-2004)
    # - LOESS 7-dom").
    def _style_key(h):
        return (h.get_color(), h.get_linestyle(), round(h.get_linewidth(), 2),
                h.get_marker(), round(h.get_alpha() or 1.0, 2))

    seen_styles = set()
    dedup_cs_handles = []
    for h in combined_cs_handles:
        key = _style_key(h)
        if key in seen_styles:
            continue
        seen_styles.add(key)
        # Strip the period prefix ("CoastSat (1984-2004) - ") from the label,
        # leaving just the part that describes the line itself.
        label = h.get_label()
        if " — " in label:
            label = label.split(" — ", 1)[1]
        h.set_label(label)
        dedup_cs_handles.append(h)

    all_handles = dedup_model_handles + dedup_cs_handles + annotation_legend_handles()

    fig.legend(handles=all_handles,
               loc="lower center",
               bbox_to_anchor=(0.5, 0.055),
               fontsize=8.5, framealpha=0.95, edgecolor="#cccccc",
               frameon=True, ncol=4)

    fig.text(
        0.5, 0.005,
        f"Model: CASCADE  |  Observed: CoastSat LRR "
        f"(LOESS {'/'.join(str(w) for w in LOESS_WINDOW_DOMAINS)}-domain windows)  |  "
        f"Comparison: {comparison_name}  |  Left/right panels each show their own active period",
        fontsize=7.5, color="#666666", ha="center", va="bottom", style="italic",
    )

    # bbox_inches=None (figure's own bounds) rather than "tight": with the
    # margins already reserved above, the figure's own bounding box is
    # correct, and "tight" recomputation is what caused the runaway size.
    fig.savefig(out_path, dpi=300, facecolor="white")
    plt.close(fig)
    print(f"✓ Saved two-period:  {os.path.basename(out_path)}")


def plot_residuals(run_data, cs_series, active_period, out_path, comparison_name):
    """
    Residual panel — each model run minus the active CoastSat LOESS curve.
    Helps identify where each run over- or under-predicts relative to observations.
    Only produced if PLOT_RESIDUALS = True.
    """
    # Find the active CoastSat series and the designated residuals window
    active_cs = next(
        (cs for cs in cs_series if cs["period_start"] == active_period), None
    )
    if active_cs is None:
        print("  ⚠️  No active CoastSat series found — skipping residuals plot.")
        return

    # Select the configured residuals window; fall back to the last (widest) if missing
    active_win = next(
        (w for w in active_cs["windows"] if w["window"] == RESIDUALS_LOESS_WINDOW),
        active_cs["windows"][-1],
    )
    cs_gis_x    = active_win["gis_x"]
    cs_smoothed = active_win["smoothed"]
    all_gis_ids = np.arange(FIRST_FILE_NUMBER, LAST_FILE_NUMBER + 1, dtype=int)

    # Interpolate CoastSat to full 1–90 grid (fills any gaps)
    cs_interp = np.interp(
        all_gis_ids.astype(float),
        cs_gis_x.astype(float),
        cs_smoothed,
    )

    fig, ax = plt.subplots(figsize=(15, 6))
    fig.patch.set_facecolor("white")
    ax.set_facecolor("white")
    add_geographic_annotations(ax)

    # Per-run summary stats, shown in the legend so the single most useful
    # diagnostic number (typical fit quality) is visible at a glance instead
    # of requiring a visual estimate from the curve alone.
    run_stats = []
    for run in run_data:
        residual = run["rates"] - cs_interp
        mean_abs = np.nanmean(np.abs(residual))
        rmse     = np.sqrt(np.nanmean(residual ** 2))
        run_stats.append((run, residual, mean_abs, rmse))
        ax.plot(run["gis_ids"], residual,
                color=run["color"], linewidth=2.0, zorder=4,
                label=f"{run['label']}  (MAE={mean_abs:.2f}, RMSE={rmse:.2f})")
        ax.fill_between(run["gis_ids"], residual, 0,
                        where=(residual > 0), alpha=0.08, color=run["color"])
        ax.fill_between(run["gis_ids"], residual, 0,
                        where=(residual < 0), alpha=0.08, color=run["color"])

    _style_ax(ax, ylabel="Model − CoastSat (m/yr)")
    ax.set_xlabel("CASCADE Model Domain (500 m alongshore)",
                  fontsize=11, fontweight="bold", labelpad=4)
    ax.set_title(
        f"Residuals: Model − CoastSat (active period)  |  {comparison_name}",
        fontsize=12, fontweight="bold", pad=12, color="#1a2a3a"
    )

    run_handles = [Line2D([0], [0], color=r["color"], lw=2.0,
                          label=f"{r['label']}  (MAE={mae:.2f}, RMSE={rmse:.2f})")
                   for r, _, mae, rmse in run_stats]
    cs_label_handle = Line2D([0], [0], color="gray", lw=1.0, ls="--",
                              label=f"Reference: {active_cs['label']} LOESS {active_win['window']}-dom")

    # Reserve bottom margin BEFORE placing the legend so it stays in
    # on-canvas figure-fraction coordinates (see plot_annotated for the
    # full explanation - loc="best" previously placed the legend directly
    # on top of data, see uploaded screenshot).
    fig.subplots_adjust(bottom=0.24)
    fig.legend(handles=run_handles + [cs_label_handle],
               loc="lower center",
               bbox_to_anchor=(0.5, 0.03),
               fontsize=9, framealpha=0.95, edgecolor="#cccccc",
               frameon=True, ncol=2)

    fig.savefig(out_path, dpi=300, facecolor="white")
    plt.close(fig)
    print(f"✓ Saved residuals:   {os.path.basename(out_path)}")
    for run, _, mae, rmse in run_stats:
        print(f"    {run['label']:<20s} MAE={mae:.3f} m/yr   RMSE={rmse:.3f} m/yr")


# =============================================================================
# MAIN
# =============================================================================

def main():
    # ── Resolve comparison name ───────────────────────────────────────────────
    global COMPARISON_NAME
    if COMPARISON_NAME is None:
        COMPARISON_NAME = "_vs_".join(r["label"].replace(" ", "") for r in RUNS_TO_COMPARE)
    comparison_name = COMPARISON_NAME

    # ── Resolve run colors (auto-gradient unless a run overrides 'color') ────
    RUNS_TO_COMPARE[:] = assign_run_colors(RUNS_TO_COMPARE)

    # ── Resolve active CoastSat period ────────────────────────────────────────
    global ACTIVE_PERIOD_START
    if ACTIVE_PERIOD_START is None:
        years = [r["start_year"] for r in RUNS_TO_COMPARE]
        ACTIVE_PERIOD_START = max(set(years), key=years.count)
        print(f"ACTIVE_PERIOD_START inferred from runs: {ACTIVE_PERIOD_START}")
    active_period = ACTIVE_PERIOD_START

    # ── Create comparison folder ──────────────────────────────────────────────────
    out_dir = os.path.join(COMPARISON_ROOT_DIR, comparison_name)
    os.makedirs(out_dir, exist_ok=True)
    print(f"\nComparison:   {comparison_name}")
    print(f"Output dir:   {out_dir}")
    print("=" * 70)

    # ── Load model runs ───────────────────────────────────────────────────────
    print("\nLoading CASCADE run rate CSVs...")
    run_data = []
    for run_cfg in RUNS_TO_COMPARE:
        try:
            gis_ids, rates, run_dir = load_run_rates(
                run_cfg["run_name"],
                period  = run_cfg.get("period"),
                preset  = run_cfg.get("preset"),
                arm     = run_cfg.get("arm"),
                run_dir = run_cfg.get("run_dir"),
            )
            run_data.append(dict(
                run_name   = run_cfg["run_name"],
                label      = run_cfg["label"],
                color      = run_cfg["color"],
                start_year = run_cfg["start_year"],
                gis_ids    = gis_ids,
                rates      = rates,
                run_dir    = run_dir,
            ))
            print(f"  ✓ {run_cfg['label']:<20s}  "
                  f"rate range {rates.min():.2f}–{rates.max():.2f} m/yr  "
                  f"color={run_cfg['color']}  "
                  f"({run_dir})")
        except FileNotFoundError as e:
            print(f"  ❌ SKIPPED '{run_cfg['run_name']}': {e}")

    if not run_data:
        print("\n❌ No valid runs loaded — RUNS_TO_COMPARE is empty "
              "or every entry failed to resolve. See Section 3.")
        sys.exit(1)

    print(f"\n  {len(run_data)} run(s) loaded successfully.")

    # ── Load CoastSat transects + apply LOESS at transect resolution ─────────
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
            print(f"  ✓ LOESS applied: window={w} domains "
                  f"({w * DOMAIN_SPACING_M / 1000.0:.1f} km)  "
                  f"frac={frac:.3f}  ({ds['label']})")
            windows.append(dict(window=w, gis_x=gis_x, smoothed=smoothed, frac=frac))
        cs_series.append(dict(
            label                 = ds["label"],
            period_start          = ds["period_start"],
            transect_domains      = domain_ids,       # one entry per transect (domain ID)
            transect_rates        = lrr_values,        # one entry per transect (LRR y)
            transect_along_coast  = along_coast_m,    # one entry per transect (physical x)
            windows               = windows,           # LOESS curves at domain resolution
        ))

    if not cs_series:
        print("  ⚠️  No CoastSat data loaded — plots will show model lines only.")
    else:
        # Catch the case that caused real confusion before: one period's
        # CoastSat CSV failed to load (e.g. wrong path) while another
        # period's loaded fine. cs_series then silently has only one
        # period's data, so any run expecting the missing period gets
        # compared against the wrong period's CoastSat with no error - only
        # a one-line warning printed earlier, easy to miss.
        loaded_periods = {cs["period_start"] for cs in cs_series}
        needed_periods = {r["start_year"] for r in RUNS_TO_COMPARE}
        missing_periods = needed_periods - loaded_periods
        if missing_periods:
            print("=" * 70)
            print(f"WARNING: CoastSat data for period(s) {sorted(missing_periods)} "
                  f"did NOT load (see 'Transect CSV not found' warning above for "
                  f"the exact path checked).")
            print(f"         Runs with start_year in {sorted(missing_periods)} will "
                  f"be compared against the WRONG period's CoastSat data, or show "
                  f"no CoastSat overlay at all, in figures that need it.")
            print(f"         Fix the path in COASTSAT_DATASETS (Section 4) before "
                  f"trusting any figure that includes these runs.")
            print("=" * 70)

    # ── Produce figures ───────────────────────────────────────────────────────
    print("\nGenerating figures...")

    plot_diagnostic(
        run_data, cs_series, active_period,
        os.path.join(out_dir, f"{comparison_name}_diagnostic.png"),
        comparison_name,
    )

    plot_annotated(
        run_data, cs_series, active_period,
        os.path.join(out_dir, f"{comparison_name}_annotated.png"),
        comparison_name,
    )

    if PLOT_TWO_PERIOD and cs_series:
        plot_two_period(
            run_data, cs_series,
            os.path.join(out_dir, f"{comparison_name}_two_period.png"),
            comparison_name,
        )

    if PLOT_RESIDUALS and cs_series:
        plot_residuals(
            run_data, cs_series, active_period,
            os.path.join(out_dir, f"{comparison_name}_residuals.png"),
            comparison_name,
        )

    print(f"\n✓ All figures saved to:\n  {out_dir}")
    print("=" * 70)


if __name__ == "__main__":
    main()
