"""
Derive the source/sink (background erosion) field from the residual between smoothed CoastSat and a base run.

    python scripts/input_prep/7-source-sink/2-calibrate/be_zone_residual_fit.py
    python scripts/input_prep/7-source-sink/2-calibrate/be_zone_residual_fit.py   # HAT_BE_BASE_PRESET=calibBE for an iteration pass

Zones fixed at pass 0, magnitudes iterated; writes the metrics, a paste-ready
rates block and diagnostic figures for be_apply_fit_to_config.py. Details: scripts/input_prep/7-source-sink/README.md.

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-09-30
"""

import json
import os
import re
import sys
from pathlib import Path


# Stop a console encoding from killing a finished computation
def _never_die_on_a_print():
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

# The pipeline's own code builds both inputs below
_HERE = Path(__file__).resolve()
PROJECT_BASE_DIR = next(_p for _p in _HERE.parents
                        if (_p / "pyproject.toml").exists())
if not (PROJECT_BASE_DIR / "pyproject.toml").exists():
    raise RuntimeError(
        f"CASCADE repo root not found: {PROJECT_BASE_DIR} has no "
        f"pyproject.toml. This file expects to live in "
        f"scripts/input_prep/7-source-sink/2-calibrate/.")
SCRIPTS_DIR = PROJECT_BASE_DIR / "scripts"
if str(SCRIPTS_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPTS_DIR))

from cascade_pipeline.coastsat_lowess import (        # noqa: E402
    CoastSatDataset,
    LowessConfig,
    build_coastsat_series,
    compute_domain_means,
)
from cascade_pipeline.hindcast import build_target_table   # noqa: E402
from cascade_pipeline.run_layout import resolve as resolve_run_file  # noqa: E402
from cascade_pipeline.run_registry import (                # noqa: E402
    CALIBRATION_ARM, preset_dir_for)
from site_layer.hatteras_site_config import HATTERAS_DOMAINS          # noqa: E402
from site_layer.hatteras_site_config import HATTERAS_BE_RATES_CALIBRATED  # noqa: E402
from site_layer.hatteras_site_config import HATTERAS_PERIODS              # noqa: E402
from site_layer.hat_figure_style import (                             # noqa: E402
    apply_style, figsize, save, caption, town_bands, open_frame,
    DOMAIN_AXIS_LABEL, C, C_1984, C_1997, C_1984_FILL, C_1997_FILL,
    INK, INK_MUTED, halo, _title)
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.gridspec as gridspec
from tqdm import tqdm

apply_style()

# Config — edit paths and thresholds here

# CoastSat transect-level LRR: smoothed at transect resolution, then averaged
from site_layer.hat_observed_rates import COASTSAT_LRR_ROOT as COASTSAT_BASE  # noqa: E402
# --- CONFIG ------------------------------------------------------------------
# The pair being fitted: HAT_BE_PERIODS names two starts; P1/P2 mean earlier/later
_PERIOD_ENV = os.environ.get("HAT_BE_PERIODS", "").strip()
PERIOD_STARTS = (tuple(int(x) for x in _PERIOD_ENV.split(","))
                 if _PERIOD_ENV else (1984, 2004))   # see DEFAULT_PERIOD_STARTS
if len(PERIOD_STARTS) != 2:
    raise SystemExit(
        f"\nHAT_BE_PERIODS must name exactly two period starts, got "
        f"{PERIOD_STARTS!r}. The fit compares an earlier window with a "
        f"later one.\n")
for _st in PERIOD_STARTS:
    if _st not in HATTERAS_PERIODS:
        raise SystemExit(
            f"\nperiod start {_st} is not in HATTERAS_PERIODS "
            f"({sorted(HATTERAS_PERIODS)}).\n")

DEFAULT_PERIOD_STARTS = (1984, 2004)
P1_START, P2_START = PERIOD_STARTS
P1_END = HATTERAS_PERIODS[P1_START]["end_year"]
P2_END = HATTERAS_PERIODS[P2_START]["end_year"]
P1_LABEL = f"{P1_START}\u2013{P1_END}"
P2_LABEL = f"{P2_START}\u2013{P2_END}"
P1_COASTSAT_CSV = str(COASTSAT_BASE / f"{P1_START}_{P1_END}"
                      / "transect_lrr_full.csv")
P2_COASTSAT_CSV = str(COASTSAT_BASE / f"{P2_START}_{P2_END}"
                      / "transect_lrr_full.csv")

# The base run per period: edgeBE, full management, no groin; iterate with calibBE (README)
BASE_PRESET   = os.environ.get("HAT_BE_BASE_PRESET", "edgeBE").strip() or "edgeBE"
# BASE_SCENARIO was here until 2026-09-18
RAW_RUNS_DIR  = PROJECT_BASE_DIR / "output" / "raw_runs"
# Where a CONCLUDED experiment's forcing arms are kept
ARM_RUNS_DIR  = PROJECT_BASE_DIR / "output" / "calibration" / "hs" / "runs"

# The section 8 settings, matching the runner
LOWESS_CONFIG  = LowessConfig(window_domains=(7,), skip_southern_domains=10)
TARGET_WINDOW = 7   # 10 until 2026-09-28, with the runner

# HAT_BE_OUTPUT_DIR redirects every output
from site_layer import hat_source_sink as _be  # noqa: E402
_PAIR_TAG = _be.pair_tag(P1_START, P1_END, P2_START, P2_END)
OUTPUT_DIR = (os.environ.get("HAT_BE_OUTPUT_DIR", "").strip()
              or str(_be.calibrate_dir(_PAIR_TAG)))

# Figures are read out of the data tree, not out of scripts/
FIG_DIR = (os.environ.get("HAT_BE_OUTPUT_DIR", "").strip()
           or str(_be.figures_dir(_PAIR_TAG)))

# Column names in CoastSat CSVs
LRR_COL    = "median_lrr"   # use median — more robust to outlier transects
DOMAIN_COL = "domain_number"

# CASCADE structure
NUM_REAL_DOMAINS = 90
START_REAL_INDEX = 15        # buffer domains before domain 1
CASCADE_SIGN     = -1        # x_s_TS increases landward = erosion → flip to standard
# P1_START/P1_END and P2_START/P2_END are derived from PERIOD_STARTS above.

# Correction thresholds

# Minimum smoothed residual magnitude to warrant any correction at all
SIGNIFICANCE_THRESHOLD = 0.5   # m/yr

# Minimum number of adjacent domains with |residual| > threshold to form a zone
RATE_COLUMN = "lrr_m_yr"

# WHICH BASE RUN THE RESIDUAL IS DERIVED FROM
GROIN_AWARE_BASE_RUN = True

MIN_ZONE_WIDTH = 3   # domains

# If |P1_correction - P2_correction| exceeds this
SHIFT_THRESHOLD = 0.75   # m/yr

# LOWESS smoothing

# Fraction of data used for each local regression (larger = smoother)
LOWESS_FRAC = 7 / 90  # exactly 7 domains at 90 total
LOWESS_WINDOW_DOMAINS = 7  # the actual window width LOWESS_FRAC is calibrated to hit

# Domains 1..N excluded from the LOWESS fit (Buxton groin influence)
GROIN_EXCLUDE_THROUGH_DOMAIN = 10
# -----------------------------------------------------------------------------

# Manual overrides

# Domain-level corrections that override the LOWESS-derived value
MANUAL_OVERRIDES = {
    # No overrides — pure LOWESS-derived values against the current baseline
}

# Locked domains

# Solved independently (D1, D90): forced to 0.0 and kept out of the significance tests

# Frozen zone set

# The zones are fixed at pass 0; only the magnitudes iterate (history in README)
def frozen_zones(period_start):
    try:
        return FROZEN_ZONE_DOMAINS[period_start]
    except KeyError:
        raise SystemExit(
            f"\nno frozen zone set for period {period_start}. Solved "
            f"periods: {sorted(FROZEN_ZONE_DOMAINS)}.\n\n"
            f"This is not a missing file -- it is a judgement that has not "
            f"been made. Each zone in this table carries a named physical "
            f"mechanism, and inventing one for {period_start} by reusing "
            f"another period's zones would assert that the same mechanisms "
            f"act over a different window.\n\n"
            f"Run with HAT_BE_FREEZE=off to see the candidate zones for "
            f"{period_start} without applying them, then add the ones you "
            f"accept to FROZEN_ZONE_DOMAINS in\n  {__file__}\n") from None


FROZEN_ZONE_DOMAINS = {
    1984: (
        5, 6, 7, 8, 10, 11, 12, 13, 27, 28, 29, 30, 31, 32, 33,
        34, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 68, 69, 70,
        71, 72, 73, 74, 75, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88,
        89
    ),
    2004: (
        9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 27,
        28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43,
        44, 50, 51, 52, 53, 54, 55, 62, 63, 64, 65, 66, 67, 68,
        69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 83, 84, 85, 86, 87,
        88, 89
    ),
}

# Domains reserved for the groin module

# D5-D7 is the Buxton groin's own footprint, and the residual left there is the GROIN's residual
GROIN_RESERVED_DOMAINS = (5, 6, 7)

LOCKED_DOMAINS = {
    1:  "Solved directly via buffer-cell reproduction (see GIS 1 value)",
    90: "Solved directly via buffer-cell reproduction (see GIS 90 value)",
}

# Physical zone definitions

# These are your prior hypotheses about where physical mechanisms operate
PHYSICAL_ZONES = {
    "Cape Point / Shoal Dynamics":  (1,  10,  "Cape/shoal attachment-detachment cycle + post-Isabel recovery"),
    "Buxton–Avon Transition":       (9,  20,  "Post-Isabel geomorphic recovery, background SLR erosion"),
    "Avon":                         (21, 31,  "Pier-induced sediment shadow, nourishment interactions"),
    "Mid-island":                   (32, 59,  "Background SLR-driven erosion, minimal local forcing"),
    "Wimble Shoals Influence":      (60, 74,  "Offshore shoal migration — sediment delivery to nearshore"),
    "Tri-Village / Rodanthe":       (75, 83,  "Persistent erosion hotspot, event-driven, infrastructure effects"),
    "Pea Island NWR":               (84, 90,  "Oregon Inlet dynamics, northern Wimble Shoals influence"),
}

# Display-only shortenings for the in-place zone strip labels (fig_be_rates)
ZONE_DISPLAY_NAMES = {
    "Cape Point / Shoal Dynamics":  "Cape Point",
    "Buxton–Avon Transition":       "Buxton-Avon",
    "Wimble Shoals Influence":      "Wimble Shoals",
    "Tri-Village / Rodanthe":       "Tri-Village",
}

# Alongshore annotation

# The village spans come from the site config through `town_bands()`, so this file cannot disagree with it
ANN_WIMBLE_SHOALS = (60, 74)
ANN_PIERS   = {"Avon pier": 26, "Rodanthe pier": 79}
ANN_GROINS  = {"Buxton groin": 5.5}

# Type comes from hat_figure_style.apply_style()
FONT_ANNOT  = 7.0    # structure names written inside a panel
FONT_STRIP  = 7.0    # names written on a one-line strip

# CASCADE loader

JOINT_FIT_JSON = (PROJECT_BASE_DIR / "output" / "calibration" / "groin"
                  / "joint_fit.json")


# The (M, fraction) the joint fit settled on, or None
def _fitted_groin(preset=BASE_PRESET):
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


# The groin base run carrying exactly the fitted (M, f), if present
def _groin_run_at(period_dir, stem, fitted, tolerance=1e-6):
    if not period_dir.exists():
        return None
    want_m, want_f = fitted
    for run in sorted(period_dir.glob(f"{stem}*_groin")):
        if not run.is_dir() or run.name.endswith("_nogroin"):
            continue
        if "reloc" in run.name or "nonourish" in run.name:
            continue
        meta = resolve_run_file(run, "metadata_json", run.name)
        if not meta.exists():
            continue
        try:
            groin = json.loads(meta.read_text()).get("groin", {})
        except (OSError, ValueError):
            continue
        got_m = groin.get("trapping_rate_m_yr")
        # The floor is not stored as its own field
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


# The forcing arm HAT_BE_HS selects, as run_registry spells arms
def _wave_arm():
    # Hatteras_ms is not on the path for this script the way SCRIPTS_DIR is
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


# The base run directory for one period, resolved from what is on disk
def base_run_dir(period_start, period_end):
    # Runs are filed [<forcing arm>/]<period>/<preset>/, the arm being absent at the calibration climate
    arm = _wave_arm()
    if arm == CALIBRATION_ARM:
        # The matrix, wherever the registry files it (raw_runs/matrix/ since 2026-09-16
        period_dir = preset_dir_for(RAW_RUNS_DIR, (period_start, period_end),
                                    BASE_PRESET)
    else:
        # Hs_experiment/runs/ keeps its 2026-09-02 shape, <arm>/<period>/ <preset>/, and is a closed experiment
        period_dir = (Path(ARM_RUNS_DIR) / arm
                      / f"{period_start}_{period_end}" / BASE_PRESET)
    stem = f"HAT_{period_start}_{period_end}_{BASE_PRESET}_road_bdm"

    # GROIN-ON BASE RUN, WHEN ONE EXISTS AT THE FITTED (M, f)
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


# Per-GIS-domain modelled LRR, m/yr, (+) seaward
def load_model_lrr(period_start, period_end):
    run_dir = base_run_dir(period_start, period_end)
    # Resolved, not joined: the rate CSV's path depends on the run layout
    csv_path = resolve_run_file(run_dir, "rate_csv", run_dir.name)
    if not csv_path.is_file():
        raise FileNotFoundError(
            f"no shoreline change rate CSV in {run_dir}")
    print(f"  model  {run_dir.name}  [{RATE_COLUMN}]")
    frame = pd.read_csv(csv_path)
    if RATE_COLUMN not in frame.columns:
        raise KeyError(
            f"{csv_path.name} has no {RATE_COLUMN!r} column. It predates the "
            f"LRR estimator; re-run the base run, or backfill it with "
            f"scripts/input_prep/7-source-sink/1-prepare/"
            f"backfill_run_lrr.py.")
    return frame.set_index("gis_domain")[RATE_COLUMN]


# (raw_per_domain_mean, target) for one period, m/yr, (+) seaward
def load_observed(period_start, csv_path):
    series = build_coastsat_series(
        [CoastSatDataset(label=f"CoastSat {period_start}",
                         period_start=period_start, csv_path=csv_path)],
        period_start, LOWESS_CONFIG)
    if not series:
        raise FileNotFoundError(f"CoastSat transects failed to load: "
                                f"{csv_path}")
    cs = series[0]
    table = build_target_table(cs, LOWESS_CONFIG, HATTERAS_DOMAINS,
                               TARGET_WINDOW)
    target = pd.Series(np.asarray(table["target_lrr_m_yr"], dtype=float),
                       index=np.asarray(table["gis_domain"], dtype=int))

    gis_x, means = compute_domain_means(
        cs["transect_domains"], cs["transect_rates"],
        HATTERAS_DOMAINS.first_gis_id, HATTERAS_DOMAINS.last_gis_id)
    raw = pd.Series(np.asarray(means, dtype=float),
                    index=np.asarray(gis_x, dtype=int))
    return raw, target


# Extract per-domain LRR (m/yr) from CASCADE base-run NPZ
def load_cascade_lrr(npz_path, start_year, end_year):
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


# Smoothing

# Apply LOWESS smoothing to a domain-indexed array, handling NaNs
def lowess_smooth(values, frac=LOWESS_FRAC):
    domains = np.arange(1, len(values) + 1, dtype=float)
    mask    = ~np.isnan(values)
    if mask.sum() < 5:
        return values.copy()
    smoothed_valid = lowess(values[mask], domains[mask],
                            frac=frac, return_sorted=False)
    out = np.full_like(values, np.nan)
    out[mask] = smoothed_valid
    return out


# LOWESS-smooth an observed rate, the groin zone (1..exclude_through) left out and passed through raw
def smooth_shoreline_rate(raw_rate, exclude_through=GROIN_EXCLUDE_THROUGH_DOMAIN,
                          window_domains=LOWESS_WINDOW_DOMAINS):
    n_valid = len(raw_rate) - exclude_through
    frac    = window_domains / n_valid

    fit_input = raw_rate.copy()
    fit_input[:exclude_through] = np.nan                       # groin zone never enters the fit
    smoothed  = lowess_smooth(fit_input, frac=frac)              # NaN-aware LOWESS
    smoothed[:exclude_through] = raw_rate[:exclude_through]     # restore raw, unsmoothed values
    return smoothed


# Zone identification

# Contiguous runs of domains with |smoothed_residual| > threshold, at least min_width wide
def identify_correction_zones(smoothed_residual, min_width=MIN_ZONE_WIDTH,
                               threshold=SIGNIFICANCE_THRESHOLD):
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


# Return the physical zone name for a given domain number
def assign_physical_zone(domain):
    for zone_name, (d0, d1, _) in PHYSICAL_ZONES.items():
        if d0 <= domain <= d1:
            return zone_name
    return "Unassigned"


# BE rate computation

# For each domain, determine the appropriate BE correction strategy
def compute_be_rates(raw_p1, raw_p2, smooth_p1, smooth_p2):
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

        # Groin-reserved domains emit 0.0, so an iteration adds nothing here
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

        # Locked domains are forced to 0.0: their rates are solved independently
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

        # The smoothed residual, applied only where warranted
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

    # Exploratory pass (HAT_BE_FREEZE=off): a diagnosis, never a field to apply
    if os.environ.get("HAT_BE_FREEZE", "").strip().lower() == "off":
        kept = [int(d) for d in frame.index
                if frame.loc[d, "be_hindcast_p1"] or frame.loc[d, "be_hindcast_p2"]]
        print("  FREEZE OFF: zone set NOT applied. Candidate domains with a "
              f"correction: {kept}")
        print("  These are CANDIDATES, not a field. Name a mechanism for each "
              "zone you accept, add it to FROZEN_ZONE_DOMAINS, and re-run "
              "without HAT_BE_FREEZE to produce a field that can be applied.")
        return frame

    # Hold the zone set fixed -- see FROZEN_ZONE_DOMAINS
    for period, column in ((P1_START, "be_hindcast_p1"),
                           (P2_START, "be_hindcast_p2")):
        zones = frozen_zones(period)
        outside = [d for d in frame.index if d not in zones]
        frame.loc[outside, column] = 0.0
    inside_either = set(frozen_zones(P1_START)) | set(frozen_zones(P2_START))
    outside_both = [d for d in frame.index if d not in inside_either]
    for column in ("be_forecast_continue", "be_forecast_revert",
                   "be_forecast_neutral"):
        frame.loc[outside_both, column] = 0.0
    print(f"  Frozen zone set: {len(frozen_zones(P1_START))} domains "
          f"{P1_LABEL}, {len(frozen_zones(P2_START))} {P2_LABEL}; "
          f"corrections outside withheld")
    return frame


# Annotation helper

# The alongshore furniture every panel shares
def annotate_ax(ax, ylim, villages=True, wimble=True, thresholds=True,
                label_at="bottom"):
    ymin, ymax = ylim
    yspan = ymax - ymin
    ax.set_xlim(0.5, NUM_REAL_DOMAINS + 0.5)
    if villages:
        town_bands(ax, where="top", strip=0.05, fontsize=FONT_ANNOT)
    if wimble:
        town_bands(ax, where="bottom", strip=0.05, fontsize=FONT_ANNOT,
                   shade="0.90",
                   spans={"Wimble Shoals": ANN_WIMBLE_SHOALS})
    ty, tva = ((0.93, "top") if label_at == "top" else (0.02, "bottom"))
    for name, dom in list(ANN_PIERS.items()) + list(ANN_GROINS.items()):
        ax.axvline(dom, color=INK_MUTED, lw=0.7, ls=(0, (3, 3)), zorder=2)
        ax.text(dom + 0.8, ymin + ty * yspan, name, rotation=90,
                ha="left", va=tva, fontsize=FONT_ANNOT, color=INK_MUTED,
                zorder=6, path_effects=halo(2.2))
    ax.axhline(0, color=INK, lw=0.7, zorder=2)
    if thresholds:
        ax.axhline(SIGNIFICANCE_THRESHOLD, color=INK_MUTED, lw=0.6, ls=":",
                   zorder=2)
        ax.axhline(-SIGNIFICANCE_THRESHOLD, color=INK_MUTED, lw=0.6, ls=":",
                   zorder=2)


# Collapse a per-domain physical_zone Series into (start_domain, end_domain, zone_name) runs
def find_zone_runs(zone_series, domains):
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


# Collapse a sorted list of domains into (first, last) runs
def _contiguous(domains):
    runs = []
    for d in sorted(domains):
        if runs and d == runs[-1][1] + 1:
            runs[-1] = (runs[-1][0], d)
        else:
            runs.append((d, d))
    return runs


# One centred label per zone run on the strip, its font shrunk until it fits the zone's width
def label_zone_runs(ax, fig, runs, y=0.5, fontsize_start=FONT_STRIP,
                    fontsize_min=6.5, pad_frac=0.90):
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


# Figure 1 — Diagnostic: raw residual, smoothed, zone identification

# Five stacked panels
def plot_diagnostic(cs_p1, cs_p2, casc_p1, casc_p2,
                    cs_p1_smooth, cs_p2_smooth,
                    raw_p1, raw_p2, smooth_p1, smooth_p2,
                    results, out_path):
    domains = np.arange(1, NUM_REAL_DOMAINS + 1)
    fig, axes = plt.subplots(5, 1, figsize=figsize("double", height=9.4),
                             sharex=True, constrained_layout=True,
                             gridspec_kw={"height_ratios": [3, 3, 3, 3, 1.0]})

    rate_panels = (
        (0, axes[0], P1_LABEL, cs_p1, cs_p1_smooth, casc_p1, C_1984),
        (1, axes[1], P2_LABEL, cs_p2, cs_p2_smooth, casc_p2, C_1997),
    )
    for i, ax, label, raw, smooth, model, colour in rate_panels:
        ax.plot(domains, raw, "o-", ms=2.2, lw=0.7, color=colour, alpha=0.55,
                zorder=3,
                label="CoastSat linear regression rate, unsmoothed")
        ax.plot(domains, smooth, "-", lw=1.6, color=C["REF"], zorder=5,
                label="CoastSat rate, LOWESS-smoothed (the target)")
        ax.plot(domains, model, "-", lw=1.0, color=C["ACCENT"], alpha=0.9,
                zorder=4, label="model base run, linear regression rate")
        ax.axvline(GROIN_EXCLUDE_THROUGH_DOMAIN + 0.5, color=INK_MUTED,
                   lw=0.7, ls=":", zorder=2)
        ax.set_ylabel("shoreline rate (m/yr)")
        _title(ax, i, f"shoreline rate, {label}")
        lo = min(np.nanmin([raw, smooth, model]) * 1.2, -3)
        hi = max(np.nanmax([raw, smooth, model]) * 1.2, 3)
        ylim = (lo, hi + 0.30 * (hi - lo))    # headroom for the legend
        ax.set_ylim(ylim)
        ax.legend(loc="upper center", bbox_to_anchor=(0.5, 0.94), ncol=2,
                  fontsize=6.5, framealpha=1.0, handlelength=1.5,
                  columnspacing=1.0).set_zorder(10)
        open_frame(ax)
        annotate_ax(ax, ylim, thresholds=False)

    resid_panels = (
        (2, axes[2], "1984\u20132004", raw_p1, smooth_p1, C_1984, C_1984_FILL,
         "correction_warranted_p1", FROZEN_ZONE_DOMAINS[1984]),
        (3, axes[3], "2004\u20132024", raw_p2, smooth_p2, C_1997, C_1997_FILL,
         "correction_warranted_p2", FROZEN_ZONE_DOMAINS[2004]),
    )
    for i, ax, label, raw, smooth, colour, fill, column, frozen in resid_panels:
        ax.bar(domains, raw, width=0.7, color=C["BASE_FILL"], zorder=2,
               label="residual from the unsmoothed rates")
        ax.plot(domains, smooth, "-", lw=1.6, color=colour, zorder=3,
                label="residual that drives the correction")
        # Shade warranted zones, distinguishing APPLIED from WITHHELD
        sig = results[column].values
        for k, dom in enumerate(domains):
            if not sig[k]:
                continue
            if dom in frozen and dom not in GROIN_RESERVED_DOMAINS:
                ax.axvspan(dom - 0.5, dom + 0.5, color=fill, alpha=0.5, lw=0,
                           zorder=0)
            else:
                ax.axvspan(dom - 0.5, dom + 0.5, facecolor="none",
                           edgecolor=C["BASE"], hatch="///", linewidth=0.0,
                           alpha=0.6, zorder=0)
        ax.set_ylabel("residual (m/yr)\nobserved \u2212 modelled")
        _title(ax, i, f"residual, {label}")
        all_vals = np.concatenate([raw[~np.isnan(raw)],
                                   smooth[~np.isnan(smooth)]])
        lo = min(np.nanmin(all_vals) * 1.2, -3)
        hi = max(np.nanmax(all_vals) * 1.2, 3)
        ylim = (lo, hi + 0.22 * (hi - lo))    # headroom for the legend
        ax.set_ylim(ylim)
        handles, _labels = ax.get_legend_handles_labels()
        handles += [mpatches.Patch(facecolor=fill, alpha=0.5, edgecolor="none",
                                   label="correction applied"),
                    mpatches.Patch(facecolor="none", edgecolor=C["BASE"],
                                   hatch="///", label="diagnosed, withheld")]
        ax.legend(handles=handles, loc="upper center",
                  bbox_to_anchor=(0.5, 0.94), ncol=2, fontsize=6.5,
                  framealpha=1.0, handlelength=1.5,
                  columnspacing=1.0).set_zorder(10)
        open_frame(ax)
        annotate_ax(ax, ylim)

    # Panel 5: what each domain was given
    ax = axes[4]
    strat_colors = {"zero": "white", "stable": C["REF"],
                    "shifting": C["ACCENT"], "locked": C["BASE"],
                    "groin-reserved": C["BASE_FILL"]}
    for dom in domains:
        strat = results.loc[dom, "strategy"]
        ax.bar(dom, 1, width=0.94,
               color=strat_colors.get(strat.rstrip("*"), C["BASE_FILL"]),
               edgecolor=C["GRID"], linewidth=0.3, zorder=2)
    ax.set_yticks([])
    ax.set_ylabel("strategy")
    ax.set_xlabel(DOMAIN_AXIS_LABEL)
    ax.set_ylim(0, 1)
    ax.set_xlim(0.5, NUM_REAL_DOMAINS + 0.5)
    town_bands(ax, where="top", strip=0.30, fontsize=FONT_STRIP)
    _title(ax, 4, "correction strategy")
    patches = [
        mpatches.Patch(facecolor="white", ec=C["BASE"], label="none (within noise)"),
        mpatches.Patch(color=C["REF"], label="one value for both periods"),
        mpatches.Patch(color=C["ACCENT"], label="period-specific"),
        mpatches.Patch(color=C["BASE"], label="boundary domain, solved separately"),
        mpatches.Patch(color=C["BASE_FILL"], label="reserved for the groin"),
    ]
    ax.legend(handles=patches, loc="upper center",
              bbox_to_anchor=(0.5, -0.70), ncol=5, fontsize=7, frameon=False)

    caption(fig, (
        "Background-erosion zone identification along Hatteras Island; domain 1 "
        "is at Cape Point and domain 90 at Pea Island, 500 m per domain. "
        "(a, b) the CoastSat linear regression rate at each domain, unsmoothed "
        "and after the LOWESS smoothing applied north of domain "
        f"{GROIN_EXCLUDE_THROUGH_DOMAIN} (domains 1\u2013"
        f"{GROIN_EXCLUDE_THROUGH_DOMAIN} pass through unsmoothed because the "
        "Buxton groin dominates them), against the CASCADE base run. "
        "(c, d) the residual, observed minus modelled: the grey bars use the "
        "unsmoothed rates on both sides and are shown for comparison only, "
        "while the coloured line is the residual the calibration acts on. "
        "Shaded domains met the significance test "
        f"(|residual| > {SIGNIFICANCE_THRESHOLD} m/yr over at least "
        f"{MIN_ZONE_WIDTH} adjacent domains); hatched domains met it but lie "
        "outside the frozen zone set or inside the groin's reserved footprint "
        "(domains 5\u20137) and were left uncorrected. "
        "(e) what each domain was given: one value for both periods where the "
        f"two differ by less than {SHIFT_THRESHOLD} m/yr, period-specific "
        "values otherwise. Village spans are shaded along the top edge and the "
        "Wimble Shoals reach along the bottom; the piers and the groin are "
        "marked with dashed rulers. A diagnostic figure, not a manuscript one."))

    save(fig, out_path)
    plt.close(fig)
    print(f"  Diagnostic figure saved \u2192 {out_path}")


# Figure 2 — Final BE rates: hindcast + forecast scenarios

# The CALIBRATED FIELD as the model actually carries it, from the config
def _field_from_config(domains):
    out = {}
    for period, tag in ((1984, "p1"), (2004, "p2")):
        table = HATTERAS_BE_RATES_CALIBRATED[period]
        out[tag] = np.array([0.0 if d in LOCKED_DOMAINS else table.get(d, 0.0)
                             for d in domains], dtype=float)
    return out["p1"], out["p2"]


# The field per domain for both periods, with the residuals
def plot_be_rates(results, out_path):
    domains = np.arange(1, NUM_REAL_DOMAINS + 1)
    field_p1, field_p2 = _field_from_config(domains)
    fig, axes = plt.subplots(3, 1, figsize=figsize("double", height=6.6),
                             sharex=True, constrained_layout=True,
                             gridspec_kw={"height_ratios": [4, 4, 1.2]})

    # Panel 1: the hindcast field, both periods
    ax = axes[0]
    be_p1 = field_p1
    be_p2 = field_p2
    w = 0.40
    ax.bar(domains - w / 2, be_p1, width=w, color=C_1984, zorder=3,
           label="1984\u20132004")
    ax.bar(domains + w / 2, be_p2, width=w, color=C_1997, zorder=3,
           label="2004\u20132024")
    ylim = (min(np.nanmin([be_p1, be_p2]) * 1.3, -2),
            max(np.nanmax([be_p1, be_p2]) * 1.3, 2))
    ax.set_ylim(ylim)
    ax.set_ylabel("background erosion\nrate (m/yr)")
    _title(ax, 0, "hindcast field")
    open_frame(ax)
    # Wimble Shoals is named in panel (c)
    annotate_ax(ax, ylim, wimble=False, thresholds=False, label_at="top")

    # Domains whose correction differs between the periods, as a strip at the foot
    shifting = [d for d, a, b in zip(domains, field_p1, field_p2)
                if abs(a - b) >= SHIFT_THRESHOLD]
    town_bands(ax, where="bottom", strip=0.035, label=False,
               shade=C["ACCENT_FILL"],
               spans={f"_{d0}": (d0, d1) for d0, d1 in _contiguous(shifting)})
    handles, _lbl = ax.get_legend_handles_labels()
    handles.append(mpatches.Patch(color=C["ACCENT_FILL"],
                                  label="value differs between the periods"))
    ax.legend(handles=handles, loc="lower center", ncol=3, fontsize=7,
              framealpha=1.0).set_zorder(10)

    # Panel 2: the three forecast scenarios
    ax = axes[1]
    # The three scenarios are DEFINITIONS over the two hindcast fields, not separate fits
    be_cont = field_p2
    be_rev = field_p1
    be_neut = (field_p1 + field_p2) / 2.0

    ax.fill_between(domains, np.minimum(be_rev, be_cont),
                    np.maximum(be_rev, be_cont), color=C["BASE_FILL"],
                    lw=0, zorder=1, label="range spanned by the two")
    ax.plot(domains, be_cont, "-", lw=1.3, color=C_1997, zorder=3,
            label="2004\u20132024 field carried forward")
    ax.plot(domains, be_rev, "-", lw=1.3, color=C_1984, zorder=3,
            label="1984\u20132004 field restored")
    ax.plot(domains, be_neut, ls=(0, (4, 3)), lw=1.0, color=INK_MUTED, zorder=4,
            label="mean of the two")
    ylim = (min(np.nanmin([be_cont, be_rev]) * 1.3, -2),
            max(np.nanmax([be_cont, be_rev]) * 1.3, 2))
    ax.set_ylim(ylim)
    ax.set_ylabel("background erosion\nrate (m/yr)")
    _title(ax, 1, "forecast scenarios")
    open_frame(ax)
    annotate_ax(ax, ylim, wimble=False, thresholds=False, label_at="top")
    ax.legend(loc="lower center", ncol=4, fontsize=7,
              framealpha=1.0).set_zorder(10)

    # Panel 3: the physical zones, and where a correction was applied

    # The zone strip is not a categorical colour scale
    ax = axes[2]
    # 'Correction applied' means the field is nonzero in either period
    for dom, a, b in zip(domains, field_p1, field_p2):
        corrected = bool(a) or bool(b)
        ax.bar(dom, 1, width=1.0,
               color=C["ACCENT_FILL"] if corrected else C["BASE_FILL"],
               edgecolor="white", linewidth=0.25, zorder=2)
    for _d0, d1, _name in find_zone_runs(results["physical_zone"], domains):
        ax.axvline(d1 + 0.5, color="white", lw=1.4, zorder=3)
    ax.set_yticks([])
    ax.set_ylabel("physical zone")
    ax.set_xlabel(DOMAIN_AXIS_LABEL)
    ax.set_ylim(0, 1)
    ax.set_xlim(0.5, NUM_REAL_DOMAINS + 0.5)
    _title(ax, 2, "physical zones")
    ax.legend(handles=[
        mpatches.Patch(color=C["ACCENT_FILL"], label="correction applied"),
        mpatches.Patch(color=C["BASE_FILL"], label="left at zero")],
        loc="upper center", bbox_to_anchor=(0.5, -0.70), ncol=2, fontsize=7,
        frameon=False)

    caption(fig, (
        "The calibrated background-erosion field and the forecast scenarios "
        "built from it; domain 1 is at Cape Point and domain 90 at Pea Island, "
        "500 m per domain. A positive rate is a sediment source, a negative "
        "one a sink. (a) the hindcast field for each period, applied only where "
        "the residual was significant, spatially coherent and attributable to a "
        "named process; the purple wash marks the domains whose value differs "
        f"between the periods by more than {SHIFT_THRESHOLD} m/yr. (b) three "
        "ways of carrying the field into the future -- the 2004\u20132024 field "
        "continued, the 1984\u20132004 field restored, and the mean of the two "
        "-- with the grey band showing the spread between the first two, which "
        "is the physical uncertainty the choice carries. (c) the physical zones "
        "named once each in place, shaded where a correction was applied and "
        "left pale where the domain keeps a rate of zero. Village spans are "
        "shaded along the top edge of the upper panels and the Wimble Shoals "
        "reach along the bottom; the Avon and Rodanthe piers and the Buxton "
        "groin are marked with dashed rulers."))

    # Direct in-place zone labels, fitted to each zone's own width
    fig.canvas.draw()
    zone_runs = find_zone_runs(results["physical_zone"], domains)
    label_zone_runs(ax, fig, zone_runs)

    save(fig, out_path)
    plt.close(fig)
    print(f"  BE rates figure saved \u2192 {out_path}")


# Print DOMAIN_BE_RATES dicts

# Print ready-to-paste DOMAIN_BE_RATES dicts for all scenarios and optionally write to a txt file
def print_be_dicts(results, txt_path=None):
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


# Run: residuals, zones, corrections, then the tables and figures
def main():
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    os.makedirs(FIG_DIR, exist_ok=True)

    def as_array(series):
        return np.array([series.get(d, np.nan)
                         for d in range(1, NUM_REAL_DOMAINS + 1)])

    # Observed: the section 8 target, not the domain-averaged CSV
    print("Loading CoastSat transects and building the section 8 target …")
    raw_p1_ser, tgt_p1_ser = load_observed(P1_START, P1_COASTSAT_CSV)
    raw_p2_ser, tgt_p2_ser = load_observed(P2_START, P2_COASTSAT_CSV)
    cs_p1, cs_p2 = as_array(raw_p1_ser), as_array(raw_p2_ser)
    cs_p1_smooth, cs_p2_smooth = as_array(tgt_p1_ser), as_array(tgt_p2_ser)
    print(f"  P1: {np.sum(~np.isnan(cs_p1_smooth))} target domains")
    print(f"  P2: {np.sum(~np.isnan(cs_p2_smooth))} target domains")

    # Modelled: the base run's own rate CSV
    print(f"\nLoading {BASE_PRESET} base-run rates …")
    casc_p1 = as_array(load_model_lrr(P1_START, P1_END))
    casc_p2 = as_array(load_model_lrr(P2_START, P2_END))

    domains = np.arange(1, NUM_REAL_DOMAINS + 1)
    pd.DataFrame({"domain": domains, "casc_p1": casc_p1, "casc_p2": casc_p2}
                 ).to_csv(os.path.join(OUTPUT_DIR, "cascade_base_lrr.csv"), index=False)

    # The observed curve is already smoothed -- `build_target_table` did it at transect resolution

    # Raw residual, unsmoothed on both sides: diagnostic only, it drives nothing below
    raw_p1 = cs_p1 - casc_p1
    raw_p2 = cs_p2 - casc_p2
    print(f"  Mean |raw residual| P1: {np.nanmean(np.abs(raw_p1)):.2f} m/yr")
    print(f"  Mean |raw residual| P2: {np.nanmean(np.abs(raw_p2)):.2f} m/yr")

    # Residual from the smoothed observed rate: this drives everything from here on
    print("Computing residual from smoothed shoreline rate …")
    smooth_p1 = cs_p1_smooth - casc_p1
    smooth_p2 = cs_p2_smooth - casc_p2

    # Compute BE corrections
    print("Computing BE corrections …")
    results = compute_be_rates(raw_p1, raw_p2, smooth_p1, smooth_p2)

    # Summary
    n_zero     = (results["strategy"] == "zero").sum()
    n_stable   = (results["strategy"] == "stable").sum()
    n_shifting = (results["strategy"] == "shifting").sum()
    print(f"\n  No correction needed:          {n_zero:3d} domains")
    print(f"  Stable correction (single BE): {n_stable:3d} domains")
    print(f"  Shifting (period-specific):    {n_shifting:3d} domains")

    # Save CSV
    csv_out = os.path.join(OUTPUT_DIR, "be_zone_metrics.csv")
    results.to_csv(csv_out)
    print(f"\n  Metrics CSV saved → {csv_out}")

    # Figures
    txt_out = os.path.join(OUTPUT_DIR, "DOMAIN_BE_RATES.txt")
    print_be_dicts(results, txt_path=txt_out)

    print("\nGenerating figures …")
    plot_diagnostic(
        cs_p1, cs_p2, casc_p1, casc_p2,
        cs_p1_smooth, cs_p2_smooth,
        raw_p1, raw_p2, smooth_p1, smooth_p2,
        results,
        os.path.join(FIG_DIR, "1-field", "fig_be_diagnostic.png"))

    plot_be_rates(
        results,
        os.path.join(FIG_DIR, "1-field", "fig_be_rates.png"))

    # Print dicts
    print("\nDone.")


if __name__ == "__main__":
    main()
