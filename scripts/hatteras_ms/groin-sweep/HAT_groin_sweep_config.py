#!/usr/bin/env python3
"""
Shared constants for the groin / background-erosion sweep: target, fit window and cell naming.

    python scripts/hatteras_ms/groin-sweep/HAT_groin_sweep_config.py   # imported by the orchestrator, the worker and the figures

The orchestrator and the worker must agree on these, so they live here once. Details: scripts/hatteras_ms/README.md.

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-09-30
"""

from __future__ import annotations

import re
import sys
import os
from pathlib import Path

import numpy as np
import pandas as pd

_HERE = Path(__file__).resolve()
# parents[3]: this file is in hatteras_ms/groin-sweep/; the guard makes a move fail loudly here
PROJECT_BASE_DIR = next(_p for _p in _HERE.parents
                        if (_p / "pyproject.toml").exists())
if not (PROJECT_BASE_DIR / "pyproject.toml").exists():
    raise RuntimeError(
        f"CASCADE repo root not found: {PROJECT_BASE_DIR} has no "
        f"pyproject.toml. This file expects to live in "
        f"scripts/hatteras_ms/groin-sweep/.")
SCRIPTS_DIR = PROJECT_BASE_DIR / "scripts"
# --- CONFIG ------------------------------------------------------------------
# Every groin-sweep product lives under here
GROIN_SWEEP_ROOT = PROJECT_BASE_DIR / "output" / "calibration" / "groin"
if str(SCRIPTS_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPTS_DIR))

from cascade_pipeline.coastsat_lowess import compute_domain_means  # noqa: E402
from site_layer.hatteras_site_config import (  # noqa: E402
    HATTERAS_BE_RATES_EDGE,
    HATTERAS_PERIODS,
    last_model_year,
)

# Periods and conventions

# THE CANONICAL CHAIN IS 1996 -> 2009 -> 2025, DEM to DEM (2026-10-05; was 1996 -> 2010 -> 2024 from 2026-09-17)
PERIODS = (1996, 2009)
PRESETS = ("edgeBE", "zeroBE")

DAM_TO_M = 10.0              # Barrier3D works in decameters
FLIP_SIGN_MODEL = True       # x_s_TS increases landward; flip so + = seaward

for _period in PERIODS:
    if _period not in HATTERAS_PERIODS:
        raise ValueError(
            f"sweep period {_period} is not in HATTERAS_PERIODS "
            f"({sorted(HATTERAS_PERIODS)}); the sweep and the site config "
            f"disagree about which periods exist.")

END_YEAR = {p: v["end_year"] for p, v in HATTERAS_PERIODS.items()}
# The last calendar year each period runs (MODEL_YEARS.md); END_YEAR is only the label
LAST_MODEL_YEAR = {p: last_model_year(p) for p in HATTERAS_PERIODS}


# The groin -- everything except the two swept knobs

# FOUR STRUCTURES, ONE DIPOLE -- BY DESIGN, NOT BY OVERSIGHT
GROIN_UPDRIFT_GIS = 6        # source: accretes -- also holds the field itself
GROIN_DOWNDRIFT_GIS = 5      # sink:   erodes
GROIN_INSTALL_YEAR = 1969    # confirmed construction date
GROIN_LAST_REPAIR_YEAR = 1996
GROIN_STORM_YEAR = 2003      # storm damage; end of the deterioration ramp

GROIN_DETERIORATION_DELAY_YEARS = GROIN_LAST_REPAIR_YEAR - GROIN_INSTALL_YEAR
GROIN_DETERIORATION_MODE = "linear_ramp"
GROIN_DETERIORATION_RAMP_YEARS = GROIN_STORM_YEAR - GROIN_LAST_REPAIR_YEAR

# The SAME schedule is used in both periods, deliberately

# Fraction of the peak effect that defines the fillet's edge
GROIN_EXTENT_THRESHOLD_FRAC = 0.10


# The grid

# M = 0 needs a no-groin run at the same be1; the ceiling is 110 (README)
M_VALUES = [0.0, 40.0, 50.0, 60.0, 70.0, 80.0, 95.0, 110.0, 125.0, 140.0, 160.0]
F_VALUES = [0.0, 0.2, 0.4, 0.6, 0.8, 1.0]

# be1 swept in period 1 only; -42.6 is the calibrated value itself (README)
BE1_VALUES_1984 = [-10.0, -16.0, -22.0, -28.0, -34.0, -40.0, -42.6, -46.0]
# -----------------------------------------------------------------------------


# Background-erosion values to sweep at GIS 1 for one period/preset
def be1_values(period, preset):
    if preset == "zeroBE":
        return [None]
    if period == 1984:
        return list(BE1_VALUES_1984)
    return [None]


# The fixed north-end rate for one period, from the site config
def be_gis90(period):
    return float(HATTERAS_BE_RATES_EDGE[period][90])


# The site-config be1 for a period whose be1 is not swept
def be_gis1_default(period):
    return float(HATTERAS_BE_RATES_EDGE[period][1])


# The observational target

# Raw per-domain transect means over D1-D12, not the LOWESS-smoothed table
FIT_GIS_MIN, FIT_GIS_MAX = 1, 12
FIT_DOMAINS_GIS = tuple(range(FIT_GIS_MIN, FIT_GIS_MAX + 1))

# Resolved through hat_observed_rates.py (2026-09-18), not typed.
from site_layer.hat_observed_rates import COASTSAT_LRR_ROOT as COASTSAT_DIR  # noqa: E402

# The values the M-only sweep carried as a literal table
_OBSERVED_LRR_1984_PUBLISHED = {
    1: -4.16, 2: -4.80, 3: -4.77, 4: -3.98, 5: -2.37, 6: -1.39,
    7: -2.19, 8: -1.70, 9: -2.13, 10: -2.74, 11: -2.96, 12: -2.35,
}


# Computes the per-domain observed LRR over the fit window
def _load_observed_lrr(period):
    path = COASTSAT_DIR / f"{period}_{END_YEAR[period]}" / "transect_lrr_full.csv"
    if not path.exists():
        raise FileNotFoundError(
            f"CoastSat transect file for {period}-{END_YEAR[period]} not "
            f"found at {path}. The sweep cannot score without a target.")
    frame = pd.read_csv(path)
    gis_x, means = compute_domain_means(
        frame["domain_number"].values, frame["lrr_m_yr"].values,
        FIT_GIS_MIN, FIT_GIS_MAX)
    observed = {int(g): float(v) for g, v in zip(gis_x, means)}
    missing = [g for g in FIT_DOMAINS_GIS if g not in observed]
    if missing:
        raise ValueError(
            f"{period}: no CoastSat transects fall in domains {missing}, so "
            f"the fit window is incomplete.")
    return observed


OBSERVED_LRR = {p: _load_observed_lrr(p) for p in PERIODS}

# The computed 1984 target must reproduce the published table to 2 dp
_OBSERVED_LRR_1984 = OBSERVED_LRR.get(1984) or _load_observed_lrr(1984)
for _gis, _published in _OBSERVED_LRR_1984_PUBLISHED.items():
    _computed = _OBSERVED_LRR_1984[_gis]
    if abs(_computed - _published) > 5e-3:
        raise ValueError(
            f"observed 1984-2004 LRR at D{_gis} computed as {_computed:.4f} "
            f"but the published sweep table says {_published:.2f}. The target "
            f"has moved; every period-1 sweep number predates the change.")

# The ranking target: observed updrift minus downdrift
OBSERVED_DIFFERENTIAL = {
    p: OBSERVED_LRR[p][GROIN_UPDRIFT_GIS] - OBSERVED_LRR[p][GROIN_DOWNDRIFT_GIS]
    for p in PERIODS
}

# WHAT THE DIFFERENTIAL ACTUALLY MEASURES
OBSERVED_DIFFERENTIAL_IS_BUILDING = {
    p: OBSERVED_DIFFERENTIAL[p] >= 0.0 for p in PERIODS
}

# Old name kept so nothing importing it breaks
PERIOD_DIFFERENTIAL_IS_REACHABLE = OBSERVED_DIFFERENTIAL_IS_BUILDING


# The validation reference

# The worker duplicates build_cascade and run_cascade_simulation from the hindcast runner
VALIDATION_TOLERANCE_M_YR = 5e-3
VALIDATION_REQUIRED_PRESET = "edgeBE"


# Scenario tokens that NEGATE a module
_NEGATED_SCENARIO_TOKENS = frozenset(
    {"nogroin", "nonourish", "noroad", "nobdm"})


# The raw_runs path component for this sweep's wave climate
def _wave_scope():
    from cascade_pipeline.hindcast import wave_climate_token
    # Hatteras_ms is not on the path for this package the way SCRIPTS_DIR is
    _ms = str(PROJECT_BASE_DIR / "scripts" / "hatteras_ms")
    if _ms not in sys.path:
        sys.path.insert(0, _ms)
    from HAT_hindcast_config import field_default

    fields = ("hs", "wave_period_s", "wave_asymmetry",
              "wave_angle_high_fraction")
    defaults = {name: field_default(name) for name in fields}
    values = dict(defaults)

    raw = os.environ.get("HAT_SWEEP_HS", "").strip()
    if raw:
        values["hs"] = float(raw)
    return wave_climate_token(values, defaults) or ""


# The island-offset mode a sweep runs under
def sweep_offset_mode():
    from cascade_pipeline.hindcast import ISLAND_OFFSET_MODES
    _ms = str(PROJECT_BASE_DIR / "scripts" / "hatteras_ms")
    if _ms not in sys.path:
        sys.path.insert(0, _ms)
    from HAT_hindcast_config import field_default

    mode = (os.environ.get("HAT_SWEEP_OFFSET_MODE", "").strip()
            or field_default("offset_mode"))
    if mode not in ISLAND_OFFSET_MODES:
        raise ValueError(f"HAT_SWEEP_OFFSET_MODE={mode!r} is not one of "
                         f"{ISLAND_OFFSET_MODES}")
    return mode


# "offset<mode>" for any mode but asrun, else ""
def _offset_token():
    mode = sweep_offset_mode()
    return "" if mode == "asrun" else f"offset{mode}"


# Directory of the matrix run this period's sweep validates against
def validation_run_dir(period):
    # Runs are filed <wave scope>/<period>/<preset>/
    from cascade_pipeline.run_registry import preset_dir_for
    base = preset_dir_for(PROJECT_BASE_DIR / "output" / "raw_runs",
                          (period, END_YEAR[period]), "edgeBE",
                          arm=_wave_scope() or "calibration")
    if not base.exists():
        return None, f"no run directory for {period}-{END_YEAR[period]}: {base}"

    # The offset token sits after the preset in a run name (HAT_1996_2010_edgeBE_offsetmetres_road_bdm_...)
    token = _offset_token()
    stem = (f"HAT_{period}_{END_YEAR[period]}_edgeBE"
            + (f"_{token}" if token else "") + "_road_bdm")
    matches = sorted(
        path for path in base.iterdir()
        if path.is_dir() and path.name.startswith(stem)
        and path.name.endswith("_groin")
        and not _NEGATED_SCENARIO_TOKENS & set(path.name.split("_")))

    if not matches:
        return None, (
            f"no edgeBE / roadway+beach-dune / groin run under {base}. The "
            f"sweep validates its duplicated model code against one; run the "
            f"seed stage first, or pass --skip-validation.")
    if len(matches) > 1:
        return None, (
            f"{len(matches)} candidate reference runs under {base}: "
            f"{[p.name for p in matches]}. Exactly one is expected; the "
            f"sweep will not guess which one it should match.")
    return matches[0], None


# Naming

# Directory / row name for one combination
def combo_dir_name(M, be1, fraction):
    be_token = "NA" if be1 is None else f"{be1:g}"
    if M == 0:
        return f"M0_be{be_token}"
    return f"M{M:g}_be{be_token}_f{fraction:.2f}"


# Where one period/preset sweep writes its results
def sweep_output_dir(period, preset):
    stem = f"{period}_{END_YEAR[period]}_{preset}"
    # A sweep run at a non-default Hs gets its OWN directory
    raw = os.environ.get("HAT_SWEEP_HS", "").strip()
    if raw and float(raw) != 2.5:
        stem += "_Hs" + f"{float(raw):g}".replace(".", "p")
    # The offset mode likewise (2026-09-24)
    if _offset_token():
        stem += "_" + _offset_token()
    return GROIN_SWEEP_ROOT / stem


# Where the joint fit writes, scoped by wave climate
def joint_fit_paths():
    out = GROIN_SWEEP_ROOT
    raw = os.environ.get("HAT_SWEEP_HS", "").strip()
    suffix = ""
    if raw and float(raw) != 2.5:
        suffix = "_Hs" + f"{float(raw):g}".replace(".", "p")
    if _offset_token():
        suffix += "_" + _offset_token()
    return out / f"joint_fit{suffix}.json", out / f"joint_fit{suffix}.csv"


# Every combination for one period/preset sweep
def build_grid(period, preset):
    combos = []
    for be1 in be1_values(period, preset):
        combos.append((0.0, be1, F_VALUES[0]))     # the paired baseline
        for M in M_VALUES:
            if M == 0:
                continue
            for fraction in F_VALUES:
                combos.append((M, be1, fraction))
    return combos


# Extent (measured, never fit)

# Fillet size -- the ranking metric

# Why size and not slope (README)


# Fillet size at the end of a run, relative to a paired baseline
def measure_fillet(shoreline_m, baseline_m, geometry,
                   updrift_gis=None, downdrift_gis=None):
    up = geometry.gis_to_pad(GROIN_UPDRIFT_GIS if updrift_gis is None
                             else updrift_gis)
    down = geometry.gis_to_pad(GROIN_DOWNDRIFT_GIS if downdrift_gis is None
                               else downdrift_gis)
    run = float(shoreline_m[-1, down] - shoreline_m[-1, up])
    ref = float(baseline_m[-1, down] - baseline_m[-1, up])
    return run - ref


# THE GROIN IS A SUB-GRID FEATURE
GROIN_TRANSECT_BOUNDARY = 76.5     # groin sits between transects ...0076 / ...0077
FILLET_HALFWIDTH_TRANSECTS = 6     # ~370 m each side at ~62 m transect spacing
FILLET_TREND_DOMAINS = (3, 9)      # window the regional trend is fitted over
FILLET_TREND_ORDER = 2


WETDRY_CHANGE_TABLE = (
    PROJECT_BASE_DIR / "hard-structures" / "groin" / "HAT-groin-buxton-output"
    / "shoreline_position_output" / "Change_from_wetdry_1967_D2_D12.csv")


# Fillet change over one period, from the fixed 1967 wet/dry datum
def observed_fillet_m(period, trend_order=1):
    if not WETDRY_CHANGE_TABLE.exists():
        raise FileNotFoundError(
            f"wet/dry change table not found at {WETDRY_CHANGE_TABLE}; the "
            f"fillet target cannot be built without it.")
    frame = pd.read_csv(WETDRY_CHANGE_TABLE).set_index("Domain_ID")

    series = {}
    for column in frame.columns:
        match = re.match(r"change_from_wetdry_1967_wetdry_(\d{4})", column)
        if not match:
            continue
        up = frame.loc[GROIN_UPDRIFT_GIS, column]
        down = frame.loc[GROIN_DOWNDRIFT_GIS, column]
        if pd.isna(up) or pd.isna(down):
            continue
        series.setdefault(int(match.group(1)), []).append(float(down - up))

    end = END_YEAR[period]
    years = np.array(sorted(y for y in series if period <= y <= end), dtype=float)
    if len(years) < trend_order + 2:
        raise ValueError(
            f"{period}-{end}: only {len(years)} dated wet/dry observations in "
            f"the window; too few to fit an order-{trend_order} trend.")
    values = np.array([np.mean(series[int(y)]) for y in years])
    slope = np.polyfit(years, values, trend_order)[0]
    return float(slope * (end - period))


OBSERVED_FILLET_M = {p: observed_fillet_m(p) for p in PERIODS}


# Alongshore extent of the groin's effect, from a paired baseline run
def measure_groin_extent(shoreline_m, baseline_m, geometry, updrift_gis,
                         downdrift_gis, threshold_frac):
    effect = -(np.asarray(shoreline_m)[-1] - np.asarray(baseline_m)[-1])
    peak = float(np.nanmax(np.abs(effect))) if effect.size else 0.0
    threshold = threshold_frac * peak

    def span(start_gis, step):
        # A zero peak means the two runs are identical, so the threshold is also zero
        if not np.isfinite(peak) or peak <= 0.0:
            return 0
        count, gis = 0, start_gis
        while geometry.first_gis_id <= gis <= geometry.last_gis_id:
            value = effect[geometry.gis_to_pad(gis)]
            if not np.isfinite(value) or abs(value) < threshold:
                break
            count += 1
            gis += step
        return count

    up = span(updrift_gis, +1)
    down = span(downdrift_gis, -1)
    return dict(
        peak_m=peak, threshold_m=threshold,
        updrift_domains=up, updrift_m=up * geometry.domain_spacing_m,
        downdrift_domains=down, downdrift_m=down * geometry.domain_spacing_m,
    )


if __name__ == "__main__":
    print("sweep configuration")
    for p in PERIODS:
        print(f"\n  period {p}-{END_YEAR[p]}")
        print(f"    be90                 {be_gis90(p):+g} m/yr")
        print(f"    observed D6-D5       {OBSERVED_DIFFERENTIAL[p]:+.2f} m/yr"
              + ("" if PERIOD_DIFFERENTIAL_IS_REACHABLE[p]
                 else "   NEGATIVE -- unreachable at M >= 0, reports a bound"))
        for preset in PRESETS:
            print(f"    grid {preset:<8}        {len(build_grid(p, preset)):>4} combinations")
    total = sum(len(build_grid(p, s)) for p in PERIODS for s in PRESETS)
    print(f"\n  total sweep combinations  {total}")
