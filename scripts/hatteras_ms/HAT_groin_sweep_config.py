#!/usr/bin/env python3
"""Shared constants for the groin / background-erosion sweep, both periods.

The orchestrator and the worker must agree on the target, the fit window and
the combination naming, or the orchestrator will rank rows the worker scored
against something else. Both import from here so there is one copy of each.

WHAT CHANGED FROM THE 1984-2004-ONLY VERSION
    Three things, all forced by fitting M and f jointly across both periods
    rather than fitting M alone in one period:

    1. Everything period-specific is now keyed by start year. There is no
       module-level START_YEAR; a caller states which period it means.
    2. The deterioration fraction f is a SWEPT AXIS, not a two-point bracket.
       The old version's argument for bracketing it -- that over 1984-2004 the
       cumulative trapping is M*(16 + 4f), so f barely separates from M -- is
       still true and is why f cannot be fit from period 1 alone. It is fit
       from the two periods together (see JOINT IDENTIFIABILITY below).
    3. The observed targets are COMPUTED from the CoastSat transect files
       rather than typed in. The 1984-2004 values are asserted against the
       numbers the old version carried, so the change is provably a
       refactor and not a new target.

JOINT IDENTIFIABILITY -- why two periods are needed for two knobs
    Cumulative trapping per unit M, measured from `cascade.groin.GroinCallback`
    with the documented 1969 install / 1996 onset / 7-year ramp schedule:

        1984-2004    16 + 4f      (f moves this only 16.0 -> 20.0)
        2004-2024    20f          (the run sits entirely past the 2003 ramp,
                                   so only the product M*f is identifiable)

    Period 1 alone leaves f nearly free: sliding f across its whole range needs
    only a 25% change in M to hold cumulative trapping fixed. Period 2 alone
    cannot separate M from f at all. Together they intersect:

        f = 4*B / (5*A - B),   A = period-1 constraint, B = period-2

    Period 2's observed differential is NEGATIVE (see below), which no positive
    M can produce, so B is pinned near zero and the intersection drives f
    toward 0 -- i.e. the structure stopped trapping after the 2003 storm. That
    is a RESULT, not a fitting artifact, but it does mean the period-2 leg
    reports a bound rather than an optimum, and the orchestrator says so.

Author: Hannah A. Henry, UNC CECL
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

_HERE = Path(__file__).resolve()
PROJECT_BASE_DIR = _HERE.parents[2]
SCRIPTS_DIR = PROJECT_BASE_DIR / "scripts"
if str(SCRIPTS_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPTS_DIR))

from cascade_pipeline.coastsat_loess import compute_domain_means  # noqa: E402
from hatteras_site_config import (  # noqa: E402
    HATTERAS_BE_RATES_EDGE,
    HATTERAS_PERIODS,
)

# =============================================================================
# PERIODS AND CONVENTIONS
# =============================================================================
PERIODS = (1984, 2004)
PRESETS = ("edgeBE", "zeroBE")

DAM_TO_M = 10.0              # Barrier3D works in decameters
FLIP_SIGN_MODEL = True       # x_s_TS increases landward; flip so + = seaward

for _period in PERIODS:
    if _period not in HATTERAS_PERIODS:
        raise ValueError(
            f"sweep period {_period} is not in HATTERAS_PERIODS "
            f"({sorted(HATTERAS_PERIODS)}); the sweep and the site config "
            f"disagree about which periods exist.")

END_YEAR = {p: HATTERAS_PERIODS[p]["end_year"] for p in PERIODS}


# =============================================================================
# THE GROIN -- everything except the two swept knobs
# =============================================================================
GROIN_UPDRIFT_GIS = 6        # source: accretes
GROIN_DOWNDRIFT_GIS = 5      # sink:   erodes
GROIN_INSTALL_YEAR = 1969    # confirmed construction date
GROIN_LAST_REPAIR_YEAR = 1996
GROIN_STORM_YEAR = 2003      # storm damage; end of the deterioration ramp

GROIN_DETERIORATION_DELAY_YEARS = GROIN_LAST_REPAIR_YEAR - GROIN_INSTALL_YEAR
GROIN_DETERIORATION_MODE = "linear_ramp"
GROIN_DETERIORATION_RAMP_YEARS = GROIN_STORM_YEAR - GROIN_LAST_REPAIR_YEAR

# The SAME schedule is used in both periods, deliberately. M and f are one
# structure-level pair, so period 2 must not be given a special "no
# deterioration, M already reduced" configuration -- that would hard-code the
# collapse to M*f instead of letting it emerge, and would make the period-2
# fit a different parameter from the period-1 one.

# Fraction of the peak effect that defines the fillet's edge. Matches the
# hindcast runner so a sweep extent and a matrix extent mean the same thing.
GROIN_EXTENT_THRESHOLD_FRAC = 0.10


# =============================================================================
# THE GRID
# =============================================================================
# M = 0 is the paired-baseline column, not a weak groin: measure_groin_extent
# needs a no-groin run at the SAME be1, and a baseline at a different be1
# would report the background-erosion difference as fillet. f is meaningless
# when M = 0, so the grid builder collapses those cells to one per be1.
#
# The ceiling is 110, not the 80 the M-only sweep used. Freeing f forces it:
# period-1 cumulative trapping is M*(16 + 4f), so reaching at f = 0 what
# M = 80 delivers at f = 0.9 needs M ~= 98. An 80 ceiling would clip the
# ridge exactly where the joint solution is heading.
M_VALUES = [0.0, 40.0, 50.0, 60.0, 70.0, 80.0, 95.0, 110.0]
F_VALUES = [0.0, 0.2, 0.4, 0.6, 0.8, 1.0]

# be1 is swept in period 1 only. The 2004 edgeBE values are taken as given --
# there is no prior fit behind a 2004 bracket, and inventing one would add an
# axis that resolves nothing. See the orchestrator docstring.
BE1_VALUES_1984 = [-22.0, -28.0, -34.0, -40.0, -46.0]


def be1_values(period, preset):
    """Background-erosion values to sweep at GIS 1 for one period/preset.

    Args:
        period: 1984 or 2004.
        preset: "edgeBE" or "zeroBE".

    Returns:
        A list of be1 values in m/yr, or [None] when the preset has no free
        background-erosion knob. None means "use the preset as it stands",
        which the worker turns into an all-zero field for zeroBE and into the
        site-config table values for a period whose be1 is not swept.
    """
    if preset == "zeroBE":
        return [None]
    if period == 1984:
        return list(BE1_VALUES_1984)
    return [None]


def be_gis90(period):
    """The fixed north-end rate for one period, from the site config.

    Read from HATTERAS_BE_RATES_EDGE rather than pinned here. The old version
    pinned 15.0, the table has since been re-solved to 10.0 for 1984, and the
    orchestrator's drift guard aborts on exactly that mismatch -- a pinned
    copy is a stale copy waiting to happen.

    Args:
        period: 1984 or 2004.

    Returns:
        The rate at GIS 90 in m/yr.
    """
    return float(HATTERAS_BE_RATES_EDGE[period][90])


def be_gis1_default(period):
    """The site-config be1 for a period whose be1 is not swept."""
    return float(HATTERAS_BE_RATES_EDGE[period][1])


# =============================================================================
# THE OBSERVATIONAL TARGET
# =============================================================================
# Raw per-domain transect means over D1-D12, not the LOESS-smoothed table:
# LoessConfig's skip_southern_domains is 10, so D1-D10 are raw in
# COASTSAT_TARGET anyway, and taking D11-D12 raw as well makes the whole
# window one construction instead of two. Sign: + is seaward/accreting.
FIT_GIS_MIN, FIT_GIS_MAX = 1, 12
FIT_DOMAINS_GIS = tuple(range(FIT_GIS_MIN, FIT_GIS_MAX + 1))

COASTSAT_DIR = SCRIPTS_DIR / "input_prep" / "5-scr" / "CoastSat"

# The values the M-only sweep carried as a literal table. Kept ONLY as an
# assertion target: if computing them from the transect file no longer
# reproduces these, the target moved and every published period-1 number
# needs re-checking, so that must fail loudly rather than quietly re-fit.
_OBSERVED_LRR_1984_PUBLISHED = {
    1: -4.16, 2: -4.80, 3: -4.77, 4: -3.98, 5: -2.37, 6: -1.39,
    7: -2.19, 8: -1.70, 9: -2.13, 10: -2.74, 11: -2.96, 12: -2.35,
}


def _load_observed_lrr(period):
    """Computes the per-domain observed LRR over the fit window.

    Args:
        period: 1984 or 2004.

    Returns:
        A dict of GIS domain ID to mean LRR in m/yr.

    Raises:
        FileNotFoundError: If the period's transect file is absent.
        ValueError: If any domain in the fit window has no transects.
    """
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

# Provable refactor: the computed 1984 target must reproduce the table the
# M-only sweep was scored against, to the 2 dp it was written at.
for _gis, _published in _OBSERVED_LRR_1984_PUBLISHED.items():
    _computed = OBSERVED_LRR[1984][_gis]
    if abs(_computed - _published) > 5e-3:
        raise ValueError(
            f"observed 1984-2004 LRR at D{_gis} computed as {_computed:.4f} "
            f"but the published sweep table says {_published:.2f}. The target "
            f"has moved; every period-1 sweep number predates the change.")

# The ranking metric's target: observed updrift minus downdrift. Derived from
# OBSERVED_LRR rather than retyped, so the two cannot disagree.
OBSERVED_DIFFERENTIAL = {
    p: OBSERVED_LRR[p][GROIN_UPDRIFT_GIS] - OBSERVED_LRR[p][GROIN_DOWNDRIFT_GIS]
    for p in PERIODS
}

# A positive differential is the signature a functioning groin leaves. Period 2
# is NEGATIVE, which the source/sink pair cannot produce at any M >= 0: the
# callback adds -M updrift and +M downdrift, so the modelled differential is
# non-negative by construction. The period-2 leg therefore reports a bound
# (the fit rails to the smallest cumulative trapping on the grid), and the
# orchestrator labels it as such instead of quoting it as an optimum.
PERIOD_DIFFERENTIAL_IS_REACHABLE = {
    p: OBSERVED_DIFFERENTIAL[p] >= 0.0 for p in PERIODS
}


# =============================================================================
# THE VALIDATION REFERENCE
# =============================================================================
# The worker duplicates build_cascade and run_cascade_simulation from the
# hindcast runner, which is itself a mirror of the notebook -- three copies
# that can drift. The guard runs the reference matrix run's own configuration
# through the worker and differences the two rate curves, so drift is caught
# numerically rather than by hashing source text.
#
# The reference's PARAMETERS ARE NOT PINNED HERE. They are read from that
# run's metadata JSON at sweep time, because the matrix run gets re-run with
# new values as the fit is refined -- and a pinned M or be1 here would report
# "drift" the first time it did, when nothing had drifted.
VALIDATION_TOLERANCE_M_YR = 1e-6
VALIDATION_REQUIRED_PRESET = "edgeBE"


def validation_run_dir(period):
    """Directory of the matrix run this period's sweep validates against.

    RESOLVED BY LOOKUP, NOT CONSTRUCTED. The runner's RUN_NAME is derived in
    its section 7.5 from what sections 5-6 actually built, and it carries a
    `nourish` token only when the period has fill scheduled -- so the
    reference is `..._edgeBE_road_bdm_groin` in 1984-2004 but
    `..._edgeBE_road_bdm_nourish_groin` in 2004-2024. Building the name here
    means reimplementing that derivation, and getting it wrong is silent: the
    sweep reports "reference missing" for a run that is sitting on disk under
    a name one token different.

    Args:
        period: 1984 or 2004.

    Returns:
        (path, None) for exactly one match, else (None, message).
    """
    base = PROJECT_BASE_DIR / "output" / "raw_runs" / f"{period}_{END_YEAR[period]}"
    if not base.exists():
        return None, f"no run directory for {period}-{END_YEAR[period]}: {base}"

    stem = f"HAT_{period}_{END_YEAR[period]}_edgeBE_road_bdm"
    matches = sorted(
        path for path in base.iterdir()
        if path.is_dir() and path.name.startswith(stem)
        and path.name.endswith("_groin")
        and not path.name.endswith("_nogroin"))

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


# =============================================================================
# NAMING
# =============================================================================

def combo_dir_name(M, be1, fraction):
    """Directory / row name for one combination.

    Sweep combinations must not use the matrix's RUN_NAME scheme: that name
    carries a preset token ("edgeBE") but no be1 or f value, so every
    combination would resolve to the same directory and silently overwrite
    the last.

    M = 0 collapses the f axis -- with no groin attached there is nothing to
    deteriorate, so every f would produce a byte-identical run. Those cells
    are named without an f token so the grid builder cannot enqueue six
    duplicates of the same baseline.

    Args:
        M: Groin trapping rate, m/yr.
        be1: Background erosion at GIS 1 in m/yr, or None when the period /
            preset has no swept be1.
        fraction: Groin deterioration floor.

    Returns:
        A filesystem-safe name, e.g. "M60_be-34_f0.40" or "M0_beNA".
    """
    be_token = "NA" if be1 is None else f"{be1:g}"
    if M == 0:
        return f"M0_be{be_token}"
    return f"M{M:g}_be{be_token}_f{fraction:.2f}"


def sweep_output_dir(period, preset):
    """Where one period/preset sweep writes its results."""
    return (PROJECT_BASE_DIR / "output" / "groin_sweep"
            / f"{period}_{END_YEAR[period]}_{preset}")


def build_grid(period, preset):
    """Every combination for one period/preset sweep.

    The M = 0 baseline is emitted once per be1 rather than once per (be1, f):
    with no groin attached f cannot apply, so the extra cells would be exact
    duplicates and would inflate the run count by five per be1 for nothing.

    Args:
        period: 1984 or 2004.
        preset: "edgeBE" or "zeroBE".

    Returns:
        A list of (M, be1, fraction) tuples.
    """
    combos = []
    for be1 in be1_values(period, preset):
        combos.append((0.0, be1, F_VALUES[0]))     # the paired baseline
        for M in M_VALUES:
            if M == 0:
                continue
            for fraction in F_VALUES:
                combos.append((M, be1, fraction))
    return combos


# =============================================================================
# EXTENT (measured, never fit)
# =============================================================================

def measure_groin_extent(shoreline_m, baseline_m, geometry, updrift_gis,
                         downdrift_gis, threshold_frac):
    """Alongshore extent of the groin's effect, from a paired baseline run.

    Copied from the hindcast runner's section 12.3. It lives here rather than
    in the worker because the orchestrator needs it too -- pairing a groin run
    with its baseline is a property of the grid, not of one combination -- and
    importing the worker would pull CASCADE into the orchestrator process for
    what is a few lines of numpy.

    The baseline MUST be the M = 0 run at the same be1. Pairing across
    different be1 would report the background-erosion difference as fillet.

    Args:
        shoreline_m: [time x domain] matrix from the groin run.
        baseline_m: [time x domain] matrix from the paired M = 0 run.
        geometry: DomainGeometry describing the padded array.
        updrift_gis, downdrift_gis: The groin's flanking domains.
        threshold_frac: Fraction of the peak effect defining the edge.

    Returns:
        A dict with peak_m, threshold_m, and the updrift/downdrift extents in
        domains and meters.
    """
    effect = -(np.asarray(shoreline_m)[-1] - np.asarray(baseline_m)[-1])
    peak = float(np.nanmax(np.abs(effect))) if effect.size else 0.0
    threshold = threshold_frac * peak

    def span(start_gis, step):
        # A zero peak means the two runs are identical, so the threshold is
        # also zero, every "abs(value) < 0" test is False, and the walk would
        # report the whole island as fillet.
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
