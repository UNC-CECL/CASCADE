"""
HAT_be_zone_analysis.py
========================
Derives defensible background erosion (BE) source/sink corrections for CASCADE
from the residual between a LOESS-smoothed CoastSat observed shoreline change
rate and the CASCADE base-run LRR.

Philosophy
----------
Corrections are applied only where THREE conditions are simultaneously met:
  1. The residual (smoothed observed rate minus CASCADE LRR) exceeds a
     significance threshold (not just noise)
  2. The signal is spatially coherent across multiple adjacent domains
  3. A physical mechanism can be named for that zone

This minimises the number of free parameters and maximises scientific
defensibility — every correction has a name and a reason.

Workflow
--------
  1. Load CoastSat domain-averaged LRR (P1 and P2) — the observed shoreline
     change rate
  2. Load CASCADE base-run LRR from NPZ (P1 and P2, management included, BE=0)
  3. LOESS-smooth the observed shoreline change rate (10-domain window),
     excluding domains 1-GROIN_EXCLUDE_THROUGH_DOMAIN (Buxton groin influence
     zone) from the fit entirely — those domains pass through with their raw,
     unsmoothed rate
  4. Compute residual = smoothed CoastSat rate - CASCADE LRR, per domain per
     period (the fully raw, unsmoothed residual is also retained, for
     diagnostic comparison only — it no longer drives any decision)
  5. Identify significant zones (|residual| > SIGNIFICANCE_THRESHOLD,
     spatially coherent over >= MIN_ZONE_WIDTH adjacent domains)
  6. Classify each domain: stable correction vs shifting between periods
  7. Output:
       - Diagnostic figures showing raw/smoothed rate, residual, and zone
         identification
       - DOMAIN_BE_RATES dicts for P1, P2, and three forecast scenarios
       - Summary CSV with all metrics and physical zone assignments

Forecast scenarios (for domains where P1 ≠ P2 correction)
----------------------------------------------------------
  "continue" : use P2 correction (current trajectory continues into future)
  "revert"   : use P1 correction (system returns to pre-2004 state)
  "neutral"  : use mean(P1, P2) (no prior on future state)

Usage
-----
  python HAT_be_zone_analysis.py

Dependencies
------------
  pip install pandas numpy matplotlib scipy statsmodels tqdm
"""

import json
import os
import re
import sys
from pathlib import Path


def _never_die_on_a_print():
    """Stop a console encoding from killing a finished computation.

    This script prints arrows, ellipses and plus-minus signs. A Windows
    console is cp1252, which cannot encode any of them, so `print` raises
    UnicodeEncodeError -- and on 2026-08-28 that happened at the very last
    status line, AFTER the whole calibration had been computed and BEFORE
    DOMAIN_BE_RATES.txt was written. The run looked like a crash, the numbers
    were gone, and the stale file left behind still carried the previous
    pass's date, so nothing downstream would have noticed it was old.

    Reconfiguring is preferred to ASCII-ifying every print: the next arrow
    someone types would reintroduce the bug, and the failure mode is silent
    data loss rather than a wrong character. If UTF-8 cannot be set, fall
    back to errors="replace" so an unencodable character degrades to "?"
    instead of raising.
    """
    for stream in (sys.stdout, sys.stderr):
        try:
            stream.reconfigure(encoding="utf-8", errors="replace")
        except Exception:
            try:
                stream.reconfigure(errors="replace")
            except Exception:
                pass


_never_die_on_a_print()

import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.nonparametric.smoothers_lowess import lowess

# The pipeline's own code builds both inputs below. Imported, never copied:
# this file used to reimplement the observed smoothing and the modelled LRR
# extraction, and both had drifted from what the model is scored against.
_HERE = Path(__file__).resolve()
PROJECT_BASE_DIR = _HERE.parents[4]
if not (PROJECT_BASE_DIR / "pyproject.toml").exists():
    raise RuntimeError(
        f"CASCADE repo root not found: {PROJECT_BASE_DIR} has no "
        f"pyproject.toml. This file expects to live in "
        f"scripts/input_prep/7-source-sink/loess_smooth/.")
SCRIPTS_DIR = PROJECT_BASE_DIR / "scripts"
if str(SCRIPTS_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPTS_DIR))

from cascade_pipeline.coastsat_loess import (        # noqa: E402
    CoastSatDataset,
    LoessConfig,
    build_coastsat_series,
    compute_domain_means,
)
from cascade_pipeline.hindcast import build_target_table   # noqa: E402
from cascade_pipeline.run_registry import (                # noqa: E402
    CALIBRATION_ARM, preset_dir_for)
from hatteras_site_config import HATTERAS_DOMAINS          # noqa: E402
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.gridspec as gridspec
from tqdm import tqdm

# ============================================================
# CONFIG — edit paths and thresholds here
# ============================================================

# CoastSat TRANSECT-level LRR. Not the domain-averaged summary: the target
# is smoothed at transect resolution and only then averaged, and doing it in
# the other order gives a measurably different curve.
COASTSAT_BASE = SCRIPTS_DIR / "input_prep" / "5-scr" / "CoastSat"
P1_COASTSAT_CSV = str(COASTSAT_BASE / "1984_2004" / "transect_lrr_full.csv")
P2_COASTSAT_CSV = str(COASTSAT_BASE / "2004_2024" / "transect_lrr_full.csv")

# The base run, per period. edgeBE is the preset that is ZERO at every domain
# being solved (2-89) while carrying the independently solved values at the
# locked ends -- so an interior residual is the whole correction rather than
# an increment on an existing one, and the open-boundary artifact at GIS 1/90
# is absorbed there instead of diffusing alongshore into interior terms.
#
# full_management because the observed CoastSat rate is from a real island
# that WAS managed. nogroin because domains 1-10 are already excluded from
# the LOESS fit for groin influence, and running the base with the groin on
# would push its signal into the residual and double-count it against the
# separate M/f sweep.
# ITERATION. The calibration is a ONE-SHOT solve: it measures the residual of a
# base run and imposes it as the BE field. That is exact only if imposing X m/yr
# moves the domain's LRR by X m/yr, and it does not -- BRIE diffuses an imposed
# rate alongshore, so a domain keeps only a fraction g of what it is given.
# Measured on 2026-08-24, g is not a constant: a contiguous same-signed block of
# corrections passes at g ~ 0.8-1.2, while a pattern that alternates sign at the
# grid scale is damped to g ~ 0.1. One pass therefore closes 42% (P1) and 57%
# (P2) of the misfit rather than all of it, and the shortfall is uneven.
#
# The fix is to iterate rather than to guess g: point this at the CURRENT
# calibBE runs, measure what residual is left, and ADD it to the field already
# in place (apply_be_fit.py --add). Each pass closes fraction g of whatever
# remains, so it converges whatever g turns out to be, and it needs no estimate
# of g at all. Amplifying by 1/g instead was rejected: g varies by an order of
# magnitude with wavelength, so narrow features would be amplified ~10x into
# rates that are indefensible read as sediment fluxes.
#
#   pass 0   HAT_BE_BASE_PRESET unset -> edgeBE   ->  apply_be_fit.py
#   pass 1+  HAT_BE_BASE_PRESET=calibBE           ->  apply_be_fit.py --add
#
# Stop when no zone clears SIGNIFICANCE_THRESHOLD, which is then the tolerance
# the field is converged to.
BASE_PRESET   = os.environ.get("HAT_BE_BASE_PRESET", "edgeBE").strip() or "edgeBE"
BASE_SCENARIO = "road_bdm_nogroin"
RAW_RUNS_DIR  = PROJECT_BASE_DIR / "output" / "raw_runs"
# Where a CONCLUDED experiment's forcing arms are kept. A live calibration
# probe is still written into raw_runs/ under HAT_ARM_TAG -- that is what stops
# it overwriting the base run it probes -- and is moved here when the question
# it was answering is settled, so raw_runs/ stays the production matrix.
ARM_RUNS_DIR  = PROJECT_BASE_DIR / "output" / "hs_experiment" / "runs"

# The section 8 settings, matching the runner. TARGET_WINDOW is the widest
# window; `rate_comparison` resolves the reference the same way.
LOESS_CONFIG  = LoessConfig(window_domains=(7, 10), skip_southern_domains=10)
TARGET_WINDOW = 10

# HAT_BE_OUTPUT_DIR redirects every output -- be_zone_metrics.csv,
# DOMAIN_BE_RATES*.txt, convergence_history.json, the figures. A what-if pass
# (a different Hs, a trial base run) MUST set it: the production directory
# holds a converged calibration whose stopping point is a recorded scientific
# claim, and an exploratory pass silently overwriting it would destroy the
# provenance without anyone noticing.
OUTPUT_DIR = os.environ.get("HAT_BE_OUTPUT_DIR", "").strip()     or str(_HERE.parent / "output")

# ── Column names in CoastSat CSVs ─────────────────────────────────────────────
LRR_COL    = "median_lrr"   # use median — more robust to outlier transects
DOMAIN_COL = "domain_number"

# ── CASCADE structure ─────────────────────────────────────────────────────────
NUM_REAL_DOMAINS = 90
START_REAL_INDEX = 15        # buffer domains before domain 1
CASCADE_SIGN     = -1        # x_s_TS increases landward = erosion → flip to standard
P1_START, P1_END = 1984, 2004
P2_START, P2_END = 2004, 2024

# ── Correction thresholds ─────────────────────────────────────────────────────
# Minimum smoothed residual magnitude to warrant any correction at all.
# Below this = within model noise/uncertainty, leave at zero.
SIGNIFICANCE_THRESHOLD = 0.5   # m/yr

# Minimum number of adjacent domains with |residual| > threshold to form a zone.
# Prevents correcting isolated noisy domains.
# WHICH MODEL COLUMN THE RESIDUAL IS BUILT FROM. Must be the same
# estimator as the observed target, which is a per-transect OLS slope;
# see load_model_lrr. "change_rate_m_yr" reproduces a pre-2026-08-22
# calibration and nothing else.
RATE_COLUMN = "lrr_m_yr"

# WHICH BASE RUN THE RESIDUAL IS DERIVED FROM.
#   True   prefer the groin-ON base run carrying the joint fit's (M, f),
#          so the source/sink carries only what the MODULES could not
#          explain. Falls back to the no-groin run, loudly, when the
#          joint fit has not run or no matching run exists.
#   False  always the no-groin run. Reproduces every calibration before
#          2026-08-22, and hands the source/sink the groin's D5/D6
#          signal to absorb.
#
# ORDERING. Fit the groin against zeroBE/edgeBE, which impose nothing at
# D5/D6, then recalibrate the source/sink against a run carrying the
# fitted groin. Running this with the switch on BEFORE the joint fit is
# harmless -- it warns and falls back -- but the result is a pre-groin
# calibration and should not be quoted as a post-groin one.
GROIN_AWARE_BASE_RUN = True

MIN_ZONE_WIDTH = 3   # domains

# If |P1_correction - P2_correction| exceeds this, the zone is "shifting"
# and needs period-specific values + forecast scenarios.
SHIFT_THRESHOLD = 0.75   # m/yr

# ── LOESS smoothing ───────────────────────────────────────────────────────────
# Fraction of data used for each local regression (larger = smoother).
# 10-domain window matches the CoastSat LOESS calibration window used
# throughout the dissertation. As of this version, smoothing is applied to
# the OBSERVED SHORELINE CHANGE RATE itself (before differencing against
# CASCADE) — not to the residual. LOESS_FRAC is calibrated for the full
# 90-domain array; see smooth_shoreline_rate() for how it's re-derived once
# the groin zone is excluded from the fit.
LOESS_FRAC = 0.111  # exactly 10 domains at 90 total
LOESS_WINDOW_DOMAINS = 10  # the actual window width LOESS_FRAC is calibrated to hit

# Domains 1 through this value are excluded ENTIRELY from the shoreline-rate
# LOESS fit (Buxton groin influence zone) — the groin's localized signal
# would otherwise bleed into the smoothed estimate at neighbouring domains.
# These domains always keep their raw, unsmoothed CoastSat rate.
GROIN_EXCLUDE_THROUGH_DOMAIN = 10

# ── Manual overrides ─────────────────────────────────────────────────────────
# Domain-level corrections that override the LOESS-derived value.
# Use sparingly — only where the smoothing window demonstrably under/over-corrects
# and you have a clear physical justification for the different value.
# Format: domain → (p1_override, p2_override, reason)
# Set either value to None to keep the LOESS-derived value for that period.
# Forecast scenarios for overridden domains use the same logic as normal:
#   continue = p2, revert = p1, neutral = mean(p1, p2)
MANUAL_OVERRIDES = {
    # No overrides — pure LOESS-derived values against the current baseline.
    # This script is for DISCOVERY: see what the data says before any
    # manual intervention. Once you have chosen final values (informed by
    # this comparison plus your own judgement), enter them in
    # HAT_be_apply_final_rates.py to generate the final figures and
    # forecast scenarios.
}

# ── Locked domains ────────────────────────────────────────────────────────────
# Domains where the BE rate has already been solved independently (e.g. D1 and
# D90, found by reproducing the buffer-cell shoreline change rate directly).
# These are forced to 0.0 in ALL DOMAIN_BE_RATES outputs and excluded from the
# significance/strategy calculation entirely — the script will not suggest a
# correction for these domains, since you will always supply your own value.
# Add a short note here for your own records of what each locked value is.
# ── Frozen zone set ───────────────────────────────────────────────────────────
# THE ZONES ARE THE SCIENCE; THE MAGNITUDE IS THE ARITHMETIC. Zone membership
# says "this stretch of coast has a real sediment-budget deficit, and here is
# the process". Magnitude says "deliver the amount you diagnosed, given that
# BRIE diffuses roughly half of it away". Iterating BOTH lets the second
# quietly rewrite the first: each pass re-derives zones from a NEW residual, so
# as the coherent features are satisfied, progressively less coherent ones
# cross the 0.5 m/yr threshold and get corrected. Worse, adding BE at a domain
# pushes sediment into its neighbours and changes THEIR residuals -- so later
# passes partly correct the alongshore spillover of earlier passes. That is
# bookkeeping, not geomorphology, and it never terminates.
#
# Measured 2026-08-24: two unmasked passes corrected 19 domains in 1984-2004
# and 12 in 2004-2024 that the pass-0 zone identification never selected --
# including D5-D7 in period 2, the groin's own footprint.
#
# So zones are identified ONCE, from the edgeBE residual under the ordinary
# significance and width rules, and then held fixed. Everything outside stays
# at 0.0 no matter what its residual does; that residual remains in the results
# as honest unexplained variance, which is the defensible number anyway.
#
# Regenerate ONLY by re-running pass 0 against edgeBE and re-deriving the set;
# do not edit a domain in here to chase a residual.
FROZEN_ZONE_DOMAINS = {
    1984: (
        5, 6, 7, 8, 9, 10, 11, 12, 13, 22, 27, 28, 29, 30, 31, 32, 33,
        34, 44, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 62, 68, 69, 70,
        71, 72, 73, 74, 75, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88,
        89
    ),
    2004: (
        8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 27,
        28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43,
        44, 48, 50, 51, 52, 53, 54, 55, 57, 62, 63, 64, 65, 66, 67, 68,
        69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 83, 84, 85, 86, 87,
        88, 89
    ),
}

# ── Domains reserved for the groin module ─────────────────────────────────────
# D5-D7 is the Buxton groin's own footprint, and the residual left there is the
# GROIN's residual, not a source/sink term. Measured 2026-08-24, after two
# iteration passes:
#
#     period 1   D6 = +1.72 m/yr   observed is more seaward than modelled --
#                                  M = 60 does not build enough fillet
#     period 2   D6 = -2.21 m/yr   modelled is more seaward than observed --
#                                  the fillet does not release, which a module
#                                  with trapping bounded at >= 0 CANNOT do
#
# Letting the BE field absorb those is exactly the double-count that
# GROIN_AWARE_BASE_RUN exists to prevent: the source/sink term would quietly do
# the groin's job and the groin would look better calibrated than it is.
#
# THIS WAS ALREADY HAPPENING, BY ACCIDENT. D5-D6 (period 1) and D6-D7 (period 2)
# each form a significant run only 2 domains wide, and MIN_ZONE_WIDTH = 3 means
# no zone forms, so no correction is ever applied. That is the right outcome
# reached by a rule that knows nothing about the groin -- and it would silently
# reverse if MIN_ZONE_WIDTH were ever retuned. Naming the domains here makes the
# behaviour intentional and survives that.
#
# FREEZE, NOT ZERO. These emit 0.0, which under `apply_be_fit.py --add` means
# "add nothing" -- the value already in the config is kept. That matters: at D5
# only ~16% of the standing correction is the groin, the other ~84% being Cape
# Point background measured with the groin switched OFF. Zeroing would discard
# it. Under the pass-0 replace path 0.0 would zero them, so reserve domains only
# once the field they should keep is already in place.
GROIN_RESERVED_DOMAINS = (5, 6, 7)

LOCKED_DOMAINS = {
    1:  "Solved directly via buffer-cell reproduction (see GIS 1 value)",
    90: "Solved directly via buffer-cell reproduction (see GIS 90 value)",
}

# ── Physical zone definitions ─────────────────────────────────────────────────
# These are your prior hypotheses about where physical mechanisms operate.
# The script will test whether the residual data supports them.
# Format: zone_name → (domain_start, domain_end, mechanism_description)
PHYSICAL_ZONES = {
    "Cape Point / Shoal Dynamics":  (1,  10,  "Cape/shoal attachment-detachment cycle + post-Isabel recovery"),
    "Buxton–Avon Transition":       (9,  20,  "Post-Isabel geomorphic recovery, background SLR erosion"),
    "Avon":                         (21, 31,  "Pier-induced sediment shadow, nourishment interactions"),
    "Mid-island":                   (32, 59,  "Background SLR-driven erosion, minimal local forcing"),
    "Wimble Shoals Influence":      (60, 74,  "Offshore shoal migration — sediment delivery to nearshore"),
    "Tri-Village / Rodanthe":       (75, 83,  "Persistent erosion hotspot, event-driven, infrastructure effects"),
    "Pea Island NWR":               (84, 90,  "Oregon Inlet dynamics, northern Wimble Shoals influence"),
}

# Display-only shortenings for the in-place zone strip labels (fig_be_rates).
# The canonical PHYSICAL_ZONES names above are kept everywhere else — CSV
# comparison, DOMAIN_BE_RATES comments, mechanism lookups — since the longer
# names carry more information there. Add more entries here if you want
# other zones shortened on the chart too.
ZONE_DISPLAY_NAMES = {
    "Cape Point / Shoal Dynamics":  "Cape Hatteras",
    "Buxton–Avon Transition":       "Buxton-Avon",
    "Wimble Shoals Influence":      "Wimble Shoals",
    "Tri-Village / Rodanthe":       "Tri-Village",
}

# ── Annotation colours (publication style) ───────────────────────────────────
ANN_TOWN_SPANS = {
    "Buxton":      (7,  8),
    "Avon":        (21, 31),
    "Tri-Village": (68, 83),
}
ANN_WIMBLE_SHOALS = (60, 74)
ANN_PIERS   = {"Avon Pier": 26, "Rodanthe Pier": 79}
ANN_GROINS  = {"Buxton Groin": 5.5}
ANN_C_TOWN  = "#90AFC5"
ANN_C_WIMBLE = "#E0A800"
ANN_C_PIER  = "#1565C0"
ANN_C_GROIN = "#B71C1C"

# ── Figure font sizes ─────────────────────────────────────────────────────────
# Centralised here so the whole figure can be resized for readability (e.g.
# for sharing with Laura / printing) by editing these six values instead of
# hunting through every plot function.
FONT_SUPTITLE = 18   # main figure title (fig.suptitle)
FONT_TITLE    = 13   # per-panel titles
FONT_LABEL    = 12   # axis labels (x/y)
FONT_TICK     = 11   # tick labels
FONT_LEGEND   = 11   # legend text
FONT_ANNOT    = 10   # small annotation text (town/pier/groin labels in annotate_ax)

# ============================================================
# CASCADE LOADER
# ============================================================

JOINT_FIT_JSON = (PROJECT_BASE_DIR / "output" / "groin_sweep"
                  / "joint_fit.json")


def _fitted_groin(preset=BASE_PRESET):
    """The (M, fraction) the joint fit settled on, or None.

    Args:
        preset: Which preset's fit to read. The base run is BASE_PRESET,
            so its own fit is the one that describes it.

    Returns:
        An (M, fraction) tuple, or None if the joint fit has not run or
        holds no entry for this preset.
    """
    if not JOINT_FIT_JSON.exists():
        return None
    try:
        fits = json.loads(JOINT_FIT_JSON.read_text())
    except (OSError, ValueError):
        return None
    fit = fits.get(preset)
    if not fit or "M" not in fit or "fraction" not in fit:
        return None
    return float(fit["M"]), float(fit["fraction"])


def _groin_run_at(period_dir, stem, fitted, tolerance=1e-6):
    """The groin base run carrying exactly the fitted (M, f), if present.

    Matched on the run's own metadata rather than on its directory name:
    the name records that a groin ran, not which one.

    Args:
        period_dir: <period>/<preset> directory to search.
        stem: Invariant leading part of the run name.
        fitted: (M, fraction) to match.
        tolerance: Absolute tolerance on both values.

    Returns:
        The matching Path, or None.
    """
    if not period_dir.exists():
        return None
    want_m, want_f = fitted
    for run in sorted(period_dir.glob(f"{stem}*_groin")):
        if not run.is_dir() or run.name.endswith("_nogroin"):
            continue
        if "reloc" in run.name or "nonourish" in run.name:
            continue
        meta = run / f"{run.name}_run_metadata.json"
        if not meta.exists():
            continue
        try:
            groin = json.loads(meta.read_text()).get("groin", {})
        except (OSError, ValueError):
            continue
        got_m = groin.get("trapping_rate_m_yr")
        # The floor is not stored as its own field. The runner writes the
        # deterioration as a HUMAN-READABLE STRING -- "linear_ramp, floor 0.6"
        # -- so `deterioration_fraction` is always None and, before this,
        # every candidate was skipped and the analysis fell back to the
        # no-groin base while reporting "no groin base run at the fitted
        # M, f". The runs were correct; only the lookup was wrong, and the
        # fallback is a warning rather than an error, so it produced a
        # complete pre-groin calibration that LOOKED like a post-groin one.
        got_f = groin.get("deterioration_fraction")
        if got_f is None:
            match = re.search(r"floor\s*([\d.]+)", str(groin.get("deterioration", "")))
            got_f = match.group(1) if match else None
        if got_m is None or got_f is None:
            continue
        if (abs(float(got_m) - want_m) <= tolerance
                and abs(float(got_f) - want_f) <= tolerance):
            return run
    return None


def _wave_arm():
    """The forcing arm HAT_BE_HS selects, as run_registry spells arms.

    CALIBRATION_ARM at the calibration wave climate, so every run made before
    arms existed resolves exactly where it always did. The token comes from
    wave_climate_token -- the same function the runner derives its own arm tag
    from -- rather than being spelled here, so the two cannot disagree about
    which directory a run was filed in.

    Returns:
        An arm name for preset_dir_for.
    """
    # hatteras_ms is not on the path for this script the way SCRIPTS_DIR is.
    # Added inside the function, mirroring HAT_groin_sweep_config._wave_scope,
    # so the module does not gain an import-time dependency on the runner.
    _ms = str(SCRIPTS_DIR / "hatteras_ms")
    if _ms not in sys.path:
        sys.path.insert(0, _ms)
    from cascade_pipeline.hindcast import wave_climate_token
    from HAT_hindcast_config import field_default

    fields = ("hs", "wave_period_s", "wave_asymmetry",
              "wave_angle_high_fraction")
    defaults = {name: field_default(name) for name in fields}
    values = dict(defaults)
    raw = os.environ.get("HAT_BE_HS", "").strip()
    if raw:
        values["hs"] = float(raw)
    return wave_climate_token(values, defaults) or CALIBRATION_ARM


def base_run_dir(period_start, period_end):
    """The base run directory for one period, resolved from what is on disk.

    The scenario tokens are NOT the same in both periods: period 2 has
    nourishment scheduled, so its full_management run carries a `nourish`
    token that period 1 has no reason to. Globbing the invariant part and
    excluding the arms that must not be picked is what keeps this from
    silently resolving to the wrong run when the token set changes again.

    Raises:
        FileNotFoundError: If no base run is present. Loud rather than
            falling back to another preset -- a calibration silently derived
            against the wrong base would look entirely normal downstream.
        RuntimeError: If more than one candidate matches, rather than
            guessing which run the calibration should rest on.
    """
    # Runs are filed [<forcing arm>/]<period>/<preset>/, the arm being absent
    # at the calibration climate. HAT_BE_HS points this at the tree for a
    # different wave height, which is what lets the SAME calibration method be
    # run against a differently-forced model and the two compared.
    #
    # BOTH THE ARM'S SPELLING AND THE LAYOUT COME FROM SHARED CODE. This block
    # previously spelled "waveHs3" itself and joined the path itself, so it
    # held a second copy of two rules the runner also implements -- and a
    # second copy is how the two drift into reading different directories.
    #
    # WHICH ROOT, THOUGH. The calibration arm lives in raw_runs/; an
    # off-calibration arm does NOT. The Hs 3.0 arms were moved to
    # hs_experiment/runs/ on 2026-09-02, beside the DECISION.md they are the
    # evidence for, so that raw_runs/ holds only the production matrix and its
    # sensitivity cells and no run name appears there twice. The layout INSIDE
    # each root is identical, which is why one preset_dir_for call serves both.
    arm = _wave_arm()
    root = RAW_RUNS_DIR if arm == CALIBRATION_ARM else ARM_RUNS_DIR
    period_dir = preset_dir_for(root, (period_start, period_end),
                                BASE_PRESET, arm=arm)
    stem = f"HAT_{period_start}_{period_end}_{BASE_PRESET}_road_bdm"

    # GROIN-ON BASE RUN, WHEN ONE EXISTS AT THE FITTED (M, f).
    #
    # The source/sink field is meant to carry what the MODULES could not
    # explain. Deriving it from a groin-OFF run hands it the groin's whole
    # signal at D5/D6, and a later calibBE-plus-groin run then applies both
    # -- the double count section 7 of HAT_hindcast_methods.md warns about.
    # It is not hypothetical: the 2026-08-22 calibration put -1.4 m/yr at
    # D6 in 2004-2024, absorbing exactly the fillet relaxation the groin is
    # now known to be able to produce.
    #
    # THE (M, f) MATCH IS THE POINT, not merely finding a groin run. A seed
    # run exists at PROVISIONAL values (M = 50, f = 0.9) purely to give the
    # sweep a drift-guard reference; calibrating against that would fit the
    # source/sink to a groin nobody has fitted yet. So the run must carry
    # the values in joint_fit.json, and until the joint fit has run there
    # is nothing to match and this falls back -- loudly.
    if GROIN_AWARE_BASE_RUN:
        fitted = _fitted_groin()
        if fitted is None:
            print("  WARNING: GROIN_AWARE_BASE_RUN is on but "
                  f"{JOINT_FIT_JSON.name} does not exist yet -- falling "
                  "back to the no-groin base run. The source/sink field "
                  "will absorb the groin's D5/D6 signal.")
        else:
            matched = _groin_run_at(period_dir, stem, fitted)
            if matched is not None:
                return matched
            print(f"  WARNING: no groin base run at the fitted "
                  f"M={fitted[0]:g}, f={fitted[1]:g} under {period_dir} "
                  "-- falling back to the no-groin base run. Run stage 6 "
                  "first if the source/sink should carry only what the "
                  "groin could not.")

    hits = sorted(
        d for d in period_dir.glob(f"{stem}*_nogroin")
        if d.is_dir()
        and "reloc" not in d.name          # the prescribed-relocation arm
        and "nonourish" not in d.name)     # the fills-off contrast
    if len(hits) == 1:
        return hits[0]
    if not hits:
        raise FileNotFoundError(
            f"base run not found under {period_dir}\n"
            f"  Expected the {BASE_PRESET} full_management run (groin off, "
            f"relocations off) for {period_start}-{period_end}.\n"
            f"  Run it with:\n"
            f"    python scripts/hatteras_ms/HAT_run_all.py --stages 2 "
            f"--presets {BASE_PRESET} --scenarios full_management")
    raise RuntimeError(
        f"{len(hits)} candidate base runs in {period_dir}: "
        f"{[d.name for d in hits]}. Refusing to guess which one the "
        f"calibration should rest on.")


def load_model_lrr(period_start, period_end):
    """Per-GIS-domain modelled LRR, m/yr, (+) seaward.

    Read from the run's own `*_shoreline_change_rate.csv` rather than
    re-derived from the .npz. The pipeline writes that file from the same
    array section 12 scores, so this cannot disagree with the model about
    sign, units, or padded-index alignment.

    Takes `lrr_m_yr`, the OLS slope through the run's annual states --
    NOT `change_rate_m_yr`, which is (x[-1] - x[0]) / span. The residual
    this calibration turns into a background-erosion rate is model minus
    observed, and the observed side is a per-transect OLS slope, so the
    model side has to be the same estimator or the residual carries the
    difference between two estimators as if it were a sediment budget.
    Every BE preset before 2026-08-22 was fit on the endpoint column;
    this function is the reason those values are not reproducible from
    the current pipeline without setting RATE_COLUMN back.
    """
    run_dir = base_run_dir(period_start, period_end)
    hits = sorted(run_dir.glob("*_shoreline_change_rate.csv"))
    if not hits:
        raise FileNotFoundError(
            f"no *_shoreline_change_rate.csv in {run_dir}")
    print(f"  model  {hits[0].parent.name}  [{RATE_COLUMN}]")
    frame = pd.read_csv(hits[0])
    if RATE_COLUMN not in frame.columns:
        raise KeyError(
            f"{hits[0].name} has no {RATE_COLUMN!r} column. It predates the "
            f"LRR estimator; re-run the base run, or backfill it with "
            f"scripts/input_prep/7-source-sink/backfill_lrr.py.")
    return frame.set_index("gis_domain")[RATE_COLUMN]


def load_observed(period_start, csv_path):
    """(raw_per_domain_mean, target) for one period, m/yr, (+) seaward.

    `target` is the curve the runner's section 8 builds and section 12 grades
    against: LOESS at transect resolution over along-coast distance, averaged
    to domains, with GIS 1..skip_southern_domains spliced in as raw means.
    `raw` is the unsmoothed per-domain mean, kept for the diagnostic residual
    only -- it drives nothing.
    """
    series = build_coastsat_series(
        [CoastSatDataset(label=f"CoastSat {period_start}",
                         period_start=period_start, csv_path=csv_path)],
        period_start, LOESS_CONFIG)
    if not series:
        raise FileNotFoundError(f"CoastSat transects failed to load: "
                                f"{csv_path}")
    cs = series[0]
    table = build_target_table(cs, LOESS_CONFIG, HATTERAS_DOMAINS,
                               TARGET_WINDOW)
    target = pd.Series(np.asarray(table["target_lrr_m_yr"], dtype=float),
                       index=np.asarray(table["gis_domain"], dtype=int))

    gis_x, means = compute_domain_means(
        cs["transect_domains"], cs["transect_rates"],
        HATTERAS_DOMAINS.first_gis_id, HATTERAS_DOMAINS.last_gis_id)
    raw = pd.Series(np.asarray(means, dtype=float),
                    index=np.asarray(gis_x, dtype=int))
    return raw, target


def load_cascade_lrr(npz_path, start_year, end_year):
    """Extract per-domain LRR (m/yr) from CASCADE base-run NPZ."""
    print(f"  Loading: {os.path.basename(npz_path)}")
    data    = np.load(npz_path, allow_pickle=True)
    cascade = data["cascade"][0]
    b3d     = cascade.barrier3d
    n_years = end_year - start_year
    years   = np.arange(start_year, end_year + 1)

    lrr = {}
    for dom in range(1, NUM_REAL_DOMAINS + 1):
        idx   = START_REAL_INDEX + (dom - 1)
        b3d_i = b3d[idx]
        if hasattr(b3d_i, "x_s_TS"):
            xs = np.array(b3d_i.x_s_TS, dtype=float)
        elif hasattr(b3d_i, "_x_s_TS"):
            xs = np.array(b3d_i._x_s_TS, dtype=float)
        else:
            lrr[dom] = np.nan; continue

        xs_m = xs * 10.0 * CASCADE_SIGN   # dam → m, flip sign
        nt   = min(len(xs_m), n_years + 1)
        if nt < 4:
            lrr[dom] = np.nan; continue

        slope, *_ = stats.linregress(years[:nt], xs_m[:nt])
        lrr[dom]  = slope

    return pd.Series(lrr, name="cascade_lrr")


# ============================================================
# SMOOTHING
# ============================================================

def loess_smooth(values, frac=LOESS_FRAC):
    """Apply LOESS smoothing to a domain-indexed array, handling NaNs."""
    domains = np.arange(1, len(values) + 1, dtype=float)
    mask    = ~np.isnan(values)
    if mask.sum() < 5:
        return values.copy()
    smoothed_valid = lowess(values[mask], domains[mask],
                            frac=frac, return_sorted=False)
    out = np.full_like(values, np.nan)
    out[mask] = smoothed_valid
    return out


def smooth_shoreline_rate(raw_rate, exclude_through=GROIN_EXCLUDE_THROUGH_DOMAIN,
                          window_domains=LOESS_WINDOW_DOMAINS):
    """
    LOESS-smooth a domain-indexed OBSERVED shoreline change rate, excluding
    domains 1..exclude_through entirely from the fit (Buxton groin influence
    zone) and passing those domains through unchanged with their raw rate.

    Domains > exclude_through are smoothed using ONLY data from domains
    > exclude_through — the groin-zone values never enter the regression at
    all, so they cannot bleed into the smoothed estimate near the zone
    boundary (this is stronger than merely overwriting the comparison for
    domains 1..exclude_through after smoothing over the full array).

    frac is re-derived here rather than reusing LOESS_FRAC directly: LOESS_FRAC
    (0.111) is calibrated to give a 10-domain window when fit over all 90
    domains. Once the groin zone is excluded, only 80 domains remain in the
    fit, so reusing 0.111 unchanged would narrow the window to ~8.9 domains.
    Recomputing frac = window_domains / n_valid preserves the true ~10-domain
    (~5 km) window this dissertation uses everywhere else.
    """
    n_valid = len(raw_rate) - exclude_through
    frac    = window_domains / n_valid

    fit_input = raw_rate.copy()
    fit_input[:exclude_through] = np.nan                       # groin zone never enters the fit
    smoothed  = loess_smooth(fit_input, frac=frac)              # NaN-aware LOESS
    smoothed[:exclude_through] = raw_rate[:exclude_through]     # restore raw, unsmoothed values
    return smoothed


# ============================================================
# ZONE IDENTIFICATION
# ============================================================

def identify_correction_zones(smoothed_residual, min_width=MIN_ZONE_WIDTH,
                               threshold=SIGNIFICANCE_THRESHOLD):
    """
    Find contiguous runs of domains where |smoothed_residual| > threshold
    and the run is at least min_width domains wide.

    Returns array of booleans: True = correction warranted.
    """
    domains  = np.arange(1, NUM_REAL_DOMAINS + 1)
    sig      = np.abs(smoothed_residual) > threshold
    warranted = np.zeros(NUM_REAL_DOMAINS, dtype=bool)

    # Find contiguous runs
    i = 0
    while i < NUM_REAL_DOMAINS:
        if sig[i]:
            j = i
            while j < NUM_REAL_DOMAINS and sig[j]:
                j += 1
            run_length = j - i
            if run_length >= min_width:
                warranted[i:j] = True
            i = j
        else:
            i += 1
    return warranted


def assign_physical_zone(domain):
    """Return the physical zone name for a given domain number."""
    for zone_name, (d0, d1, _) in PHYSICAL_ZONES.items():
        if d0 <= domain <= d1:
            return zone_name
    return "Unassigned"


# ============================================================
# BE RATE COMPUTATION
# ============================================================

def compute_be_rates(raw_p1, raw_p2, smooth_p1, smooth_p2):
    """
    For each domain, determine the appropriate BE correction strategy.

    Rules:
      - If smoothed residual not significant in either period → BE = 0
      - If significant in one or both periods:
          - If |correction_P1 - correction_P2| < SHIFT_THRESHOLD → stable → use mean
          - Otherwise → shifting → flag for scenario treatment, provide P1/P2/mean
    """
    domains = np.arange(1, NUM_REAL_DOMAINS + 1)
    rows    = []

    sig_p1 = identify_correction_zones(smooth_p1)
    sig_p2 = identify_correction_zones(smooth_p2)

    for i, dom in enumerate(domains):
        r1_raw  = raw_p1[i]
        r2_raw  = raw_p2[i]
        r1_sm   = smooth_p1[i]
        r2_sm   = smooth_p2[i]
        s1      = sig_p1[i]
        s2      = sig_p2[i]

        # Groin-reserved domains: emit 0.0 and skip the strategy logic, so an
        # iteration pass adds nothing here and the standing value is kept.
        # See GROIN_RESERVED_DOMAINS for why this residual is not ours.
        if dom in GROIN_RESERVED_DOMAINS:
            zone = assign_physical_zone(dom)
            rows.append({
                "domain":               dom,
                "raw_residual_p1":      r1_raw,
                "raw_residual_p2":      r2_raw,
                "smooth_residual_p1":   r1_sm,
                "smooth_residual_p2":   r2_sm,
                "correction_warranted_p1": False,
                "correction_warranted_p2": False,
                "strategy":             "groin-reserved",
                "be_hindcast_p1":       0.0,
                "be_hindcast_p2":       0.0,
                "be_forecast_continue": 0.0,
                "be_forecast_revert":   0.0,
                "be_forecast_neutral":  0.0,
                "physical_zone":        zone,
                "mechanism":            "RESERVED — the Buxton groin's own "
                                        "residual; absorbing it here would "
                                        "double-count against the M/f fit",
            })
            continue

        # Locked domains: force to 0.0, skip all significance/strategy logic.
        # These domains already have an independently-solved BE rate that you
        # will always supply yourself — the script should not suggest a value.
        if dom in LOCKED_DOMAINS:
            strategy = "locked"
            be_hindcast_p1   = 0.0
            be_hindcast_p2   = 0.0
            be_forecast_continue = 0.0
            be_forecast_revert   = 0.0
            be_forecast_neutral  = 0.0
            zone = assign_physical_zone(dom)
            mechanism = f"LOCKED — {LOCKED_DOMAINS[dom]}"
            rows.append({
                "domain":               dom,
                "raw_residual_p1":      r1_raw,
                "raw_residual_p2":      r2_raw,
                "smooth_residual_p1":   r1_sm,
                "smooth_residual_p2":   r2_sm,
                "correction_warranted_p1": False,
                "correction_warranted_p2": False,
                "strategy":             strategy,
                "be_hindcast_p1":       be_hindcast_p1,
                "be_hindcast_p2":       be_hindcast_p2,
                "be_forecast_continue": be_forecast_continue,
                "be_forecast_revert":   be_forecast_revert,
                "be_forecast_neutral":  be_forecast_neutral,
                "physical_zone":        zone,
                "mechanism":            mechanism,
            })
            continue

        # Correction value to apply: use smoothed residual (spatially coherent)
        # Only apply where warranted
        corr_p1 = r1_sm if s1 else 0.0
        corr_p2 = r2_sm if s2 else 0.0

        # Is either period significant?
        any_sig = s1 or s2

        if not any_sig:
            strategy = "zero"
            be_hindcast_p1   = 0.0
            be_hindcast_p2   = 0.0
            be_forecast_continue = 0.0
            be_forecast_revert   = 0.0
            be_forecast_neutral  = 0.0
        else:
            delta = abs(corr_p1 - corr_p2)
            if delta < SHIFT_THRESHOLD:
                strategy = "stable"
                be_mean = float(np.nanmean([corr_p1, corr_p2]))
                be_hindcast_p1   = be_mean
                be_hindcast_p2   = be_mean
                be_forecast_continue = be_mean
                be_forecast_revert   = be_mean
                be_forecast_neutral  = be_mean
            else:
                strategy = "shifting"
                be_hindcast_p1   = corr_p1
                be_hindcast_p2   = corr_p2
                # Forecast scenarios
                be_forecast_continue = corr_p2          # current state continues
                be_forecast_revert   = corr_p1          # reverts to pre-2004 state
                be_forecast_neutral  = float(np.nanmean([corr_p1, corr_p2]))

        zone = assign_physical_zone(dom)
        _, _, mechanism = PHYSICAL_ZONES.get(zone, (None, None, "Unknown"))

        # Apply manual overrides if defined for this domain
        if dom in MANUAL_OVERRIDES:
            ov_p1, ov_p2, ov_reason = MANUAL_OVERRIDES[dom]
            if ov_p1 is not None:
                be_hindcast_p1       = ov_p1
                be_forecast_revert   = ov_p1
            if ov_p2 is not None:
                be_hindcast_p2       = ov_p2
                be_forecast_continue = ov_p2
            # Recalculate neutral as mean of (possibly overridden) P1 and P2
            be_forecast_neutral = float(np.nanmean([be_hindcast_p1, be_hindcast_p2]))
            strategy = strategy + "*"  # flag as manually overridden in comparison
            mechanism = ov_reason

        rows.append({
            "domain":               dom,
            "raw_residual_p1":      r1_raw,
            "raw_residual_p2":      r2_raw,
            "smooth_residual_p1":   r1_sm,
            "smooth_residual_p2":   r2_sm,
            "correction_warranted_p1": s1,
            "correction_warranted_p2": s2,
            "strategy":             strategy,
            "be_hindcast_p1":       be_hindcast_p1,
            "be_hindcast_p2":       be_hindcast_p2,
            "be_forecast_continue": be_forecast_continue,
            "be_forecast_revert":   be_forecast_revert,
            "be_forecast_neutral":  be_forecast_neutral,
            "physical_zone":        zone,
            "mechanism":            mechanism,
        })

    frame = pd.DataFrame(rows).set_index("domain")

    # Hold the zone set fixed -- see FROZEN_ZONE_DOMAINS. Applied here rather
    # than inside the loop so the metrics CSV still records the residual and
    # the significance verdict for every domain: the diagnosis stays visible,
    # only the correction is withheld.
    for period, column in ((1984, "be_hindcast_p1"), (2004, "be_hindcast_p2")):
        outside = [d for d in frame.index if d not in FROZEN_ZONE_DOMAINS[period]]
        frame.loc[outside, column] = 0.0
    inside_either = set(FROZEN_ZONE_DOMAINS[1984]) | set(FROZEN_ZONE_DOMAINS[2004])
    outside_both = [d for d in frame.index if d not in inside_either]
    for column in ("be_forecast_continue", "be_forecast_revert",
                   "be_forecast_neutral"):
        frame.loc[outside_both, column] = 0.0
    print(f"  Frozen zone set: {len(FROZEN_ZONE_DOMAINS[1984])} domains P1, "
          f"{len(FROZEN_ZONE_DOMAINS[2004])} P2; corrections outside withheld")
    return frame


# ============================================================
# ANNOTATION HELPER
# ============================================================

def annotate_ax(ax, ylim):
    ymin, ymax = ylim; yspan = ymax - ymin
    ax.axvspan(ANN_WIMBLE_SHOALS[0]-0.5, ANN_WIMBLE_SHOALS[1]+0.5,
               color=ANN_C_WIMBLE, alpha=0.12, zorder=0)
    for name, (d0, d1) in ANN_TOWN_SPANS.items():
        ax.axvspan(d0-0.5, d1+0.5, color=ANN_C_TOWN, alpha=0.18, zorder=0)
        ax.text((d0+d1)/2, ymax - 0.02*yspan, name, ha="center", va="top",
                fontsize=FONT_ANNOT, color=ANN_C_TOWN, fontweight="bold")
    for pname, pdom in ANN_PIERS.items():
        ax.axvline(pdom, color=ANN_C_PIER, lw=0.9, ls="--", zorder=2)
        ax.text(pdom+0.3, ymin+0.72*yspan, pname, rotation=90,
                va="bottom", fontsize=FONT_ANNOT, color=ANN_C_PIER)
    for gname, gdom in ANN_GROINS.items():
        ax.axvline(gdom, color=ANN_C_GROIN, lw=0.9, ls="--", zorder=2)
        ax.text(gdom+0.3, ymin+0.62*yspan, gname, rotation=90,
                va="bottom", fontsize=FONT_ANNOT, color=ANN_C_GROIN)
    ax.axhline(0, color="black", lw=0.8, ls="--", alpha=0.5)
    ax.axhline( SIGNIFICANCE_THRESHOLD, color="#888", lw=0.6, ls=":", alpha=0.7)
    ax.axhline(-SIGNIFICANCE_THRESHOLD, color="#888", lw=0.6, ls=":", alpha=0.7)
    ax.set_xlim(0.5, NUM_REAL_DOMAINS + 0.5)


def find_zone_runs(zone_series, domains):
    """Collapse a per-domain physical_zone Series into contiguous
    (start_domain, end_domain, zone_name) runs, in domain order."""
    runs = []
    current_zone = None
    run_start = None
    for dom in domains:
        z = zone_series.loc[dom]
        if z != current_zone:
            if current_zone is not None:
                runs.append((run_start, prev_dom, current_zone))
            current_zone, run_start = z, dom
        prev_dom = dom
    runs.append((run_start, prev_dom, current_zone))
    return runs


def label_zone_runs(ax, fig, runs, y=0.5, fontsize_start=FONT_LABEL,
                    fontsize_min=6.5, pad_frac=0.90):
    """
    Place one centered text label per zone run directly on the strip,
    shrinking that label's own font (and only that one) until it fits
    within its own zone's width — so a narrow zone with a long name
    (e.g. "Cape Point / Shoal Dynamics") never bleeds into its neighbour.
    """
    renderer = fig.canvas.get_renderer()
    for d0, d1, name in runs:
        if name is None or name == "Unassigned":
            continue
        display_name = ZONE_DISPLAY_NAMES.get(name, name)
        xc = (d0 + d1) / 2.0
        x0_disp = ax.transData.transform((d0 - 0.5, y))[0]
        x1_disp = ax.transData.transform((d1 + 0.5, y))[0]
        avail_px = abs(x1_disp - x0_disp) * pad_frac

        fs = fontsize_start
        txt = ax.text(xc, y, display_name, ha="center", va="center",
                      fontsize=fs, color="black", zorder=5, clip_on=False,
                      bbox=dict(facecolor="white", alpha=0.78,
                                edgecolor="none", boxstyle="round,pad=0.15"))
        fig.canvas.draw()
        bbox = txt.get_window_extent(renderer=renderer)
        while bbox.width > avail_px and fs > fontsize_min:
            fs -= 0.5
            txt.set_fontsize(fs)
            fig.canvas.draw()
            bbox = txt.get_window_extent(renderer=renderer)


# ============================================================
# FIGURE 1 — Diagnostic: raw residual, smoothed, zone identification
# ============================================================

def plot_diagnostic(cs_p1, cs_p2, casc_p1, casc_p2,
                    cs_p1_smooth, cs_p2_smooth,
                    raw_p1, raw_p2, smooth_p1, smooth_p2,
                    results, out_path):

    domains = np.arange(1, NUM_REAL_DOMAINS + 1)
    fig, axes = plt.subplots(5, 1, figsize=(18, 18), sharex=True,
                             gridspec_kw={"height_ratios": [3, 3, 3, 3, 1.5]})

    # ── Panel 1: CoastSat vs CASCADE LRR — Period 1 ───────────────────────────
    ax = axes[0]
    ax.plot(domains, cs_p1,   "o-", ms=3, lw=1.0, color="#b2182b",
            label="CoastSat LRR (raw)", zorder=3, alpha=0.55)
    ax.plot(domains, cs_p1_smooth, "-", lw=1.8, color="#762a83",
            label=f"CoastSat LOESS-smoothed (D{GROIN_EXCLUDE_THROUGH_DOMAIN+1}+)",
            zorder=4)
    ax.plot(domains, casc_p1, "s-", ms=3, lw=1.0, color="#2166ac",
            label="CASCADE base LRR", zorder=3, alpha=0.8)
    ax.axvline(GROIN_EXCLUDE_THROUGH_DOMAIN + 0.5, color="#762a83",
               lw=0.8, ls=":", alpha=0.7, zorder=2)
    ax.set_ylabel("LRR  (m/yr)", fontsize=FONT_LABEL)
    ax.set_title(f"Period 1  (1984–2004)  |  CoastSat vs CASCADE base run  |  "
                 f"rate smoothing excludes D1-{GROIN_EXCLUDE_THROUGH_DOMAIN} (groin)",
                 fontsize=FONT_TITLE, loc="left", pad=3)
    ylim = (min(np.nanmin([cs_p1, cs_p1_smooth, casc_p1])*1.2, -3),
            max(np.nanmax([cs_p1, cs_p1_smooth, casc_p1])*1.2,  3))
    ax.set_ylim(ylim); ax.legend(fontsize=FONT_LEGEND, loc="upper center")
    ax.tick_params(labelsize=FONT_TICK)
    annotate_ax(ax, ylim)

    # ── Panel 2: CoastSat vs CASCADE LRR — Period 2 ───────────────────────────
    ax = axes[1]
    ax.plot(domains, cs_p2,   "o-", ms=3, lw=1.0, color="#b2182b",
            label="CoastSat LRR (raw)", zorder=3, alpha=0.55)
    ax.plot(domains, cs_p2_smooth, "-", lw=1.8, color="#762a83",
            label=f"CoastSat LOESS-smoothed (D{GROIN_EXCLUDE_THROUGH_DOMAIN+1}+)",
            zorder=4)
    ax.plot(domains, casc_p2, "s-", ms=3, lw=1.0, color="#2166ac",
            label="CASCADE base LRR", zorder=3, alpha=0.8)
    ax.axvline(GROIN_EXCLUDE_THROUGH_DOMAIN + 0.5, color="#762a83",
               lw=0.8, ls=":", alpha=0.7, zorder=2)
    ax.set_ylabel("LRR  (m/yr)", fontsize=FONT_LABEL)
    ax.set_title(f"Period 2  (2004–2024)  |  CoastSat vs CASCADE base run  |  "
                 f"rate smoothing excludes D1-{GROIN_EXCLUDE_THROUGH_DOMAIN} (groin)",
                 fontsize=FONT_TITLE, loc="left", pad=3)
    ylim = (min(np.nanmin([cs_p2, cs_p2_smooth, casc_p2])*1.2, -3),
            max(np.nanmax([cs_p2, cs_p2_smooth, casc_p2])*1.2,  3))
    ax.set_ylim(ylim); ax.legend(fontsize=FONT_LEGEND, loc="upper center")
    ax.tick_params(labelsize=FONT_TICK)
    annotate_ax(ax, ylim)

    # ── Panel 3: Raw vs smoothed-rate residual — Period 1 ────────────────────
    ax = axes[2]
    ax.bar(domains, raw_p1, width=0.7, color="#d9d9d9", alpha=0.6,
           label="Raw residual (unsmoothed both sides)", zorder=2)
    ax.plot(domains, smooth_p1, "-", lw=2.0, color="#b2182b",
            label="Residual (from smoothed rate)", zorder=3)
    # Shade warranted zones, distinguishing APPLIED from WITHHELD.
    # correction_warranted_* is deliberately left unmasked so the metrics CSV
    # keeps the diagnosis for every domain -- but a domain outside
    # FROZEN_ZONE_DOMAINS is diagnosed and NOT corrected, and shading the two
    # alike would tell the reader a correction was applied where none was.
    sig = results["correction_warranted_p1"].values
    for i, dom in enumerate(domains):
        if not sig[i]:
            continue
        if dom in FROZEN_ZONE_DOMAINS[1984] and dom not in GROIN_RESERVED_DOMAINS:
            ax.axvspan(dom-0.5, dom+0.5, color="#b2182b", alpha=0.12, zorder=0)
        else:
            ax.axvspan(dom-0.5, dom+0.5, facecolor="none", edgecolor="#777777",
                       hatch="///", linewidth=0.0, alpha=0.55, zorder=0)
    ax.set_ylabel("Residual  (m/yr)\nCoastSat − CASCADE", fontsize=FONT_LABEL)
    ax.set_title(f"Period 1 residual  |  solid = corrected;  hatched = diagnosed but WITHHELD (outside the frozen zone set, or reserved for the groin)  "
                 f"(|residual| > {SIGNIFICANCE_THRESHOLD} m/yr, ≥{MIN_ZONE_WIDTH} domains wide)",
                 fontsize=FONT_TITLE, loc="left", pad=3)
    all_vals = np.concatenate([raw_p1[~np.isnan(raw_p1)],
                               smooth_p1[~np.isnan(smooth_p1)]])
    ylim = (min(np.nanmin(all_vals)*1.2, -3), max(np.nanmax(all_vals)*1.2, 3))
    ax.set_ylim(ylim); ax.legend(fontsize=FONT_LEGEND, loc="upper center")
    ax.tick_params(labelsize=FONT_TICK)
    annotate_ax(ax, ylim)

    # ── Panel 4: Raw vs smoothed-rate residual — Period 2 ────────────────────
    ax = axes[3]
    ax.bar(domains, raw_p2, width=0.7, color="#d9d9d9", alpha=0.6,
           label="Raw residual (unsmoothed both sides)", zorder=2)
    ax.plot(domains, smooth_p2, "-", lw=2.0, color="#2166ac",
            label="Residual (from smoothed rate)", zorder=3)
    # Shade warranted zones, distinguishing APPLIED from WITHHELD.
    # correction_warranted_* is deliberately left unmasked so the metrics CSV
    # keeps the diagnosis for every domain -- but a domain outside
    # FROZEN_ZONE_DOMAINS is diagnosed and NOT corrected, and shading the two
    # alike would tell the reader a correction was applied where none was.
    sig = results["correction_warranted_p2"].values
    for i, dom in enumerate(domains):
        if not sig[i]:
            continue
        if dom in FROZEN_ZONE_DOMAINS[2004] and dom not in GROIN_RESERVED_DOMAINS:
            ax.axvspan(dom-0.5, dom+0.5, color="#2166ac", alpha=0.12, zorder=0)
        else:
            ax.axvspan(dom-0.5, dom+0.5, facecolor="none", edgecolor="#777777",
                       hatch="///", linewidth=0.0, alpha=0.55, zorder=0)
    ax.set_ylabel("Residual  (m/yr)\nCoastSat − CASCADE", fontsize=FONT_LABEL)
    ax.set_title(f"Period 2 residual  |  solid = corrected;  hatched = diagnosed but WITHHELD (outside the frozen zone set, or reserved for the groin)",
                 fontsize=FONT_TITLE, loc="left", pad=3)
    all_vals = np.concatenate([raw_p2[~np.isnan(raw_p2)],
                               smooth_p2[~np.isnan(smooth_p2)]])
    ylim = (min(np.nanmin(all_vals)*1.2, -3), max(np.nanmax(all_vals)*1.2, 3))
    ax.set_ylim(ylim); ax.legend(fontsize=FONT_LEGEND, loc="upper center")
    ax.tick_params(labelsize=FONT_TICK)
    annotate_ax(ax, ylim)

    # ── Panel 5: Strategy strip ───────────────────────────────────────────────
    ax = axes[4]
    # Use .get() with strip of * so any starred variant falls back gracefully
    strat_colors = {"zero": "#ffffff", "stable": "#4dac26",
                    "shifting": "#d73027", "locked": "#999999"}
    for i, dom in enumerate(domains):
        strat = results.loc[dom, "strategy"]
        ax.bar(dom, 1, width=0.9, color=strat_colors.get(strat.rstrip("*"), "#cccccc"),
               edgecolor="#cccccc", linewidth=0.3, zorder=2)
    ax.set_yticks([])
    ax.set_ylabel("BE strategy", fontsize=FONT_LABEL)
    ax.set_xlabel("CASCADE domain  (1 = Buxton  →  90 = Rodanthe)", fontsize=FONT_LABEL)
    ax.set_ylim(0, 1); ax.tick_params(labelsize=FONT_TICK)
    ax.axvspan(ANN_WIMBLE_SHOALS[0]-0.5, ANN_WIMBLE_SHOALS[1]+0.5,
               color=ANN_C_WIMBLE, alpha=0.12, zorder=0)
    for name, (d0, d1) in ANN_TOWN_SPANS.items():
        ax.axvspan(d0-0.5, d1+0.5, color=ANN_C_TOWN, alpha=0.18, zorder=0)

    patches = [
        mpatches.Patch(color="#ffffff", ec="#888", label="No correction (within noise)"),
        mpatches.Patch(color="#4dac26", label="Stable correction (single value P1=P2)"),
        mpatches.Patch(color="#d73027", label="Shifting (period-specific + forecast scenarios)"),
    ]
    ax.legend(handles=patches, fontsize=FONT_LEGEND, loc="upper right", framealpha=0.9)

    fig.suptitle(
        "Hatteras Island  |  BE source/sink zone identification\n"
        f"Significance threshold: ±{SIGNIFICANCE_THRESHOLD} m/yr  ·  "
        f"Min zone width: {MIN_ZONE_WIDTH} domains  ·  "
        f"Shift threshold: {SHIFT_THRESHOLD} m/yr",
        fontsize=FONT_SUPTITLE, fontweight="bold")
    fig.tight_layout(rect=[0, 0, 1, 0.97])
    fig.savefig(out_path, dpi=180, bbox_inches="tight")
    plt.close(fig)
    print(f"  Diagnostic figure saved → {out_path}")


# ============================================================
# FIGURE 2 — Final BE rates: hindcast + forecast scenarios
# ============================================================

def plot_be_rates(results, out_path):
    domains = np.arange(1, NUM_REAL_DOMAINS + 1)
    fig, axes = plt.subplots(3, 1, figsize=(18, 12), sharex=True,
                             gridspec_kw={"height_ratios": [4, 4, 1.5]})

    # ── Panel 1: Hindcast BE rates ────────────────────────────────────────────
    ax = axes[0]
    be_p1 = results["be_hindcast_p1"].values
    be_p2 = results["be_hindcast_p2"].values
    w = 0.38
    ax.bar(domains - w/2, be_p1, width=w, color="#4575b4", alpha=0.85,
           label="P1 hindcast BE  (1984–2004)", zorder=2)
    ax.bar(domains + w/2, be_p2, width=w, color="#d73027", alpha=0.85,
           label="P2 hindcast BE  (2004–2024)", zorder=2)
    ylim = (min(np.nanmin([be_p1, be_p2])*1.3, -2),
            max(np.nanmax([be_p1, be_p2])*1.3,  2))
    ax.set_ylim(ylim)
    ax.legend(fontsize=FONT_LEGEND, loc="lower center", ncol=2,
              frameon=True, framealpha=0.9)
    ax.set_ylabel("BE rate  (m/yr)\n+ = accretion source\n− = erosion sink", fontsize=FONT_LABEL)
    ax.set_title("Hindcast BE rates  |  only applied where physically warranted",
                 fontsize=FONT_TITLE, loc="left", pad=3)
    ax.tick_params(labelsize=FONT_TICK)
    annotate_ax(ax, ylim)

    # Mark shifting zones
    for i, dom in enumerate(domains):
        if results.loc[dom, "strategy"] == "shifting":
            ax.axvspan(dom-0.5, dom+0.5, color="#ff7f00", alpha=0.10, zorder=0)

    # ── Panel 2: Forecast scenario BE rates ───────────────────────────────────
    ax = axes[1]
    be_cont = results["be_forecast_continue"].values
    be_rev  = results["be_forecast_revert"].values
    be_neut = results["be_forecast_neutral"].values

    ax.fill_between(domains,
                    np.minimum(be_rev, be_cont),
                    np.maximum(be_rev, be_cont),
                    color="#aaaaaa", alpha=0.30, label="Scenario range (revert↔continue)")
    ax.plot(domains, be_cont, "o-", ms=3, lw=1.2, color="#d73027",
            label='Scenario "continue" (P2 state)', zorder=3)
    ax.plot(domains, be_rev,  "s-", ms=3, lw=1.2, color="#4575b4",
            label='Scenario "revert" (P1 state)', zorder=3, alpha=0.8)
    ax.plot(domains, be_neut, "^-", ms=3, lw=1.0, color="#333333",
            ls="--", label='Scenario "neutral" (mean)', zorder=3, alpha=0.7)

    ylim = (min(np.nanmin([be_cont, be_rev])*1.3, -2),
            max(np.nanmax([be_cont, be_rev])*1.3,  2))
    ax.set_ylim(ylim)
    ax.legend(fontsize=FONT_LEGEND, loc="lower center", ncol=4,
              frameon=True, framealpha=0.9)
    ax.set_ylabel("BE rate  (m/yr)", fontsize=FONT_LABEL)
    ax.set_title("Forecast scenario BE rates  |  grey band = physical uncertainty range",
                 fontsize=FONT_TITLE, loc="left", pad=3)
    ax.tick_params(labelsize=FONT_TICK)
    annotate_ax(ax, ylim)

    # ── Panel 3: Zone identification strip ────────────────────────────────────
    ax = axes[2]
    zone_palette = {
        "Cape Point / Shoal Dynamics":  "#b2182b",
        "Buxton–Avon Transition":       "#f4a582",
        "Avon":                         "#4393c3",
        "Mid-island":                   "#d9d9d9",
        "Wimble Shoals Influence":      "#E0A800",
        "Tri-Village / Rodanthe":       "#762a83",
        "Pea Island NWR":               "#4dac26",
    }
    for i, dom in enumerate(domains):
        zone  = results.loc[dom, "physical_zone"]
        strat = results.loc[dom, "strategy"]
        color = zone_palette.get(zone, "#eeeeee")
        alpha = 0.85 if strat != "zero" else 0.30
        ax.bar(dom, 1, width=0.9, color=color, alpha=alpha,
               edgecolor="white", linewidth=0.3, zorder=2)

    ax.set_yticks([])
    ax.set_ylabel("Physical zone", fontsize=FONT_LABEL)
    ax.set_xlabel("CASCADE domain  (1 = Buxton  →  90 = Rodanthe)", fontsize=FONT_LABEL)
    ax.set_ylim(0, 1); ax.tick_params(labelsize=FONT_TICK)

    ax.legend(handles=[mpatches.Patch(color="#aaaaaa", alpha=0.3,
                                      label="Faded = no correction applied")],
              fontsize=FONT_ANNOT, loc="lower right", framealpha=0.85)

    fig.suptitle(
        "Hatteras Island  |  Final BE source/sink rates\n"
        "Corrections applied only where residual is significant and physically motivated",
        fontsize=FONT_SUPTITLE, fontweight="bold")
    fig.tight_layout(rect=[0, 0, 1, 0.97])

    # Direct in-place zone labels, fitted to each zone's own width — done
    # AFTER tight_layout so the pixel-width measurements used to size each
    # label reflect the figure's final axes geometry, not a stale pre-layout
    # one.
    zone_runs = find_zone_runs(results["physical_zone"], domains)
    label_zone_runs(ax, fig, zone_runs)

    fig.savefig(out_path, dpi=180, bbox_inches="tight")
    plt.close(fig)
    print(f"  BE rates figure saved → {out_path}")


# ============================================================
# PRINT DOMAIN_BE_RATES DICTS
# ============================================================

def print_be_dicts(results, txt_path=None):
    """Print ready-to-paste DOMAIN_BE_RATES dicts for all scenarios
    and optionally write to a txt file."""
    scenarios = {
        "P1 hindcast":         "be_hindcast_p1",
        "P2 hindcast":         "be_hindcast_p2",
        "Forecast: continue":  "be_forecast_continue",
        "Forecast: revert":    "be_forecast_revert",
        "Forecast: neutral":   "be_forecast_neutral",
    }

    lines = []
    lines.append("=" * 70)
    lines.append("DOMAIN_BE_RATES — ready to paste into CASCADE run script")
    lines.append("=" * 70)

    for label, col in scenarios.items():
        nonzero = (results[col] != 0).sum()
        lines.append(f"")
        lines.append(f"# {label}  ({nonzero} domains with non-zero correction)")
        lines.append("DOMAIN_BE_RATES = {")
        for dom in results.index:
            val = results.loc[dom, col]
            strat = results.loc[dom, "strategy"]
            zone  = results.loc[dom, "physical_zone"]
            if strat == "locked":
                lines.append(f"    {dom:3d}: 0.0,  # LOCKED — use your solved value, not 0.0")
            elif val != 0.0:
                lines.append(f"    {dom:3d}: {val:+.1f},  # {zone}")
            else:
                lines.append(f"    {dom:3d}: 0.0,")
        lines.append("}")

    # Print to console
    print("\n" + "\n".join(lines))

    # Write to txt file
    if txt_path is not None:
        with open(txt_path, "w") as f:
            f.write("\n".join(lines) + "\n")
        print(f"\n  BE rates txt saved → {txt_path}")


# ============================================================
# MAIN
# ============================================================

def main():
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    def as_array(series):
        return np.array([series.get(d, np.nan)
                         for d in range(1, NUM_REAL_DOMAINS + 1)])

    # ── Observed: the section 8 target, not the domain-averaged CSV ──────────
    print("Loading CoastSat transects and building the section 8 target …")
    raw_p1_ser, tgt_p1_ser = load_observed(P1_START, P1_COASTSAT_CSV)
    raw_p2_ser, tgt_p2_ser = load_observed(P2_START, P2_COASTSAT_CSV)
    cs_p1, cs_p2 = as_array(raw_p1_ser), as_array(raw_p2_ser)
    cs_p1_smooth, cs_p2_smooth = as_array(tgt_p1_ser), as_array(tgt_p2_ser)
    print(f"  P1: {np.sum(~np.isnan(cs_p1_smooth))} target domains")
    print(f"  P2: {np.sum(~np.isnan(cs_p2_smooth))} target domains")

    # ── Modelled: the base run's own rate CSV ───────────────────────────────
    print(f"\nLoading {BASE_PRESET} base-run rates …")
    casc_p1 = as_array(load_model_lrr(P1_START, P1_END))
    casc_p2 = as_array(load_model_lrr(P2_START, P2_END))

    domains = np.arange(1, NUM_REAL_DOMAINS + 1)
    pd.DataFrame({"domain": domains, "casc_p1": casc_p1, "casc_p2": casc_p2}
                 ).to_csv(os.path.join(OUTPUT_DIR, "cascade_base_lrr.csv"), index=False)

    # The observed curve is already smoothed -- `build_target_table` did it at
    # transect resolution. `smooth_shoreline_rate` is deliberately NOT called
    # here any more: running it would LOESS an already-LOESSed curve, and the
    # second pass would flatten exactly the coherent zones this script exists
    # to detect. The function is kept for reference by the diagnostic figure.

    # ── Raw residual — fully unsmoothed on both sides, kept for diagnostic
    # comparison only. It no longer drives any decision below. ────────────────
    raw_p1 = cs_p1 - casc_p1
    raw_p2 = cs_p2 - casc_p2
    print(f"  Mean |raw residual| P1: {np.nanmean(np.abs(raw_p1)):.2f} m/yr")
    print(f"  Mean |raw residual| P2: {np.nanmean(np.abs(raw_p2)):.2f} m/yr")

    # ── Residual from the smoothed observed rate — THIS is what drives zone
    # identification, significance testing, and BE corrections from here on ──
    print("Computing residual from smoothed shoreline rate …")
    smooth_p1 = cs_p1_smooth - casc_p1
    smooth_p2 = cs_p2_smooth - casc_p2

    # ── Compute BE corrections ────────────────────────────────────────────────
    print("Computing BE corrections …")
    results = compute_be_rates(raw_p1, raw_p2, smooth_p1, smooth_p2)

    # Summary
    n_zero     = (results["strategy"] == "zero").sum()
    n_stable   = (results["strategy"] == "stable").sum()
    n_shifting = (results["strategy"] == "shifting").sum()
    print(f"\n  No correction needed:          {n_zero:3d} domains")
    print(f"  Stable correction (single BE): {n_stable:3d} domains")
    print(f"  Shifting (period-specific):    {n_shifting:3d} domains")

    # ── Save CSV ──────────────────────────────────────────────────────────────
    csv_out = os.path.join(OUTPUT_DIR, "be_zone_metrics.csv")
    results.to_csv(csv_out)
    print(f"\n  Metrics CSV saved → {csv_out}")

    # ── Figures ───────────────────────────────────────────────────────────────
    txt_out = os.path.join(OUTPUT_DIR, "DOMAIN_BE_RATES.txt")
    print_be_dicts(results, txt_path=txt_out)

    print("\nGenerating figures …")
    plot_diagnostic(
        cs_p1, cs_p2, casc_p1, casc_p2,
        cs_p1_smooth, cs_p2_smooth,
        raw_p1, raw_p2, smooth_p1, smooth_p2,
        results,
        os.path.join(OUTPUT_DIR, "fig_be_diagnostic.png"))

    plot_be_rates(
        results,
        os.path.join(OUTPUT_DIR, "fig_be_rates.png"))

    # ── Print dicts ───────────────────────────────────────────────────────────
    print("\nDone.")


if __name__ == "__main__":
    main()
