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

    THIS BLOCK USED TO SAY that period 2's negative differential pinned B near
    zero and drove f to 0 -- "the structure stopped trapping after the 2003
    storm" -- and called that a RESULT rather than a fitting artifact. That was
    wrong, and it was wrong because of the AGGREGATION, not the physics.
    Corrected 2026-08-23 against the CoastSat transects:

        the fillet is ~190 m wide; one model domain is 500 m, so the whole
        dipole fits inside a single cell. Domain averaging blends the fillet's
        peak with its own far field, and in 2004-2024 that does not merely
        weaken the signal -- it INVERTS it. Transect-scale the dipole is
        +0.48 m/yr (updrift accreting, groin-like); the domain-mean D6-D5
        difference is -2.47 m/yr (anti-groin).

    Scored against the domain-mean target the sweep drove f to 0, because the
    honest way to match a negative fillet is to trap nothing. Scored against
    the transects, the SAME runs show the groin holding ~44% of its period-1
    strength, which is what the imagery shows: deteriorated, still working.

    So `observed_fillet_m` is now built from transects (see below), and the
    period-2 leg constrains the PRODUCT M*f rather than reporting a bound at
    zero. `differential`, being a domain-mean quantity, carries the same
    aggregation defect and is retained for reporting only -- never ranked on.

Author: Hannah A. Henry, UNC CECL
"""

from __future__ import annotations

import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

_HERE = Path(__file__).resolve()
# parents[3], not [2]: this file lives in scripts/hatteras_ms/groin-sweep/.
# The guard below is what makes a future move fail here, loudly, instead of
# resolving to scripts/scripts and surfacing as a missing data file several
# imports deeper.
PROJECT_BASE_DIR = _HERE.parents[3]
if not (PROJECT_BASE_DIR / "pyproject.toml").exists():
    raise RuntimeError(
        f"CASCADE repo root not found: {PROJECT_BASE_DIR} has no "
        f"pyproject.toml. This file expects to live in "
        f"scripts/hatteras_ms/groin-sweep/.")
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
# FOUR STRUCTURES, ONE DIPOLE -- BY DESIGN, NOT BY OVERSIGHT.
# The Buxton groin field is four groins spanning northing 3901373-3901789
# (HAT_groin_shoreline_analysis_v2.py metadata, n_groin_features = 4). D6 spans
# 3901298-3901798, so the ENTIRE field falls inside a single model domain.
#
# The module therefore represents the field's CUMULATIVE effect as one
# source/sink pair rather than four. There is nothing to gain from splitting
# them: four dipoles inside one 500 m cell would sum to exactly the one dipole
# the cell can express, and `GroinCallback` requires adjacent domains anyway,
# so a per-structure representation is not expressible on this grid.
#
# WHAT THIS MAKES M. M is the trapping rate of the FIELD AS A WHOLE, not of a
# groin. It cannot be divided by four to get a per-structure rate, and it
# cannot be compared against a published single-groin trapping rate. Combined
# with the resolution note in `observed_fillet_m` -- the field's fillet is
# ~190 m wide inside a 500 m cell -- M is an effective, grid-specific,
# field-aggregate quantity. The defensible structure-level statement is the
# RATIO of trapping between periods, which is what f carries.
GROIN_UPDRIFT_GIS = 6        # source: accretes -- also holds the field itself
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
# EXTENDED ABOVE 110 ON 2026-08-24, to document the ridge rather than argue it.
# The 1984 edgeBE fit landed on M = 110, f = 0.4 -- the grid maximum -- with a
# fillet error of 0.0026 m. That is NOT "the grid ran out before the optimum":
# fillet is monotonically increasing in M at every f, so M > 110 at f = 0.4
# OVERSHOOTS the 22.11 m target. The at_grid_bound flag is mechanical.
#
# What the grid edge hides is a RIDGE. Reading the be1 = -34 surface, the target
# is matched at roughly (M 46, f 1.0), (M 57, f 0.8), (M 75, f 0.6) and
# (M 110, f 0.4) -- every one an equally good endpoint fit. Extending to 160
# adds the f = 0.2 arm of that ridge (which needs M ~ 200 to reach the target),
# so the figure shows a ridge running off the grid instead of a peak sitting on
# its edge. It resolves nothing on its own; the trajectory metric is what
# separates these cells, because they reach the same fillet along very
# different paths.
M_VALUES = [0.0, 40.0, 50.0, 60.0, 70.0, 80.0, 95.0, 110.0, 125.0, 140.0, 160.0]
F_VALUES = [0.0, 0.2, 0.4, 0.6, 0.8, 1.0]

# be1 is swept in period 1 only. The 2004 edgeBE values are taken as given --
# there is no prior fit behind a 2004 bracket, and inventing one would add an
# axis that resolves nothing. See the orchestrator docstring.
#
# THE SHALLOW END EXISTS TO BRACKET, NOT TO SEARCH. The site config's own
# edgeBE value at GIS 1 for 1984 is -24.0, which sat between the old
# bracket's first two rungs and within one step of its shallow edge. A fit
# landing on -22.0 could then mean either 'the optimum is -22' or 'the
# optimum is shallower than -22 and the grid ran out', and joint_fit's
# at_grid_bound check cannot tell those apart -- it can only say the value
# railed. -16 and -10 were added on 2026-08-22 so the configured value is
# bracketed on BOTH sides and a shallow optimum is a result rather than a
# bound. They cost 2 x 43 = 86 extra cells in the 1984 edgeBE sweep.
BE1_VALUES_1984 = [-10.0, -16.0, -22.0, -28.0, -34.0, -40.0, -46.0]


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

# WHAT THE DIFFERENTIAL ACTUALLY MEASURES. LRR[D6] - LRR[D5] is algebraically
# the OLS slope of the fillet, x_s[D5] - x_s[D6] -- its TREND across the
# window, not its size. Three consequences, all counter-intuitive, all
# measured on the 2026-08-22 1984 zeroBE sweep:
#
#   - the fillet SATURATES (~18 m at M = 40, reached by 1990, where BRIE's
#     diffusion balances the trapping), so a healthy CONSTANT groin scores
#     near ZERO;
#   - a DEGRADING groin scores NEGATIVE, because its fillet is relaxing;
#   - sign therefore reports building-vs-failing, not strong-vs-weak, and two
#     very different groins can score identically.
#
# THIS BLOCK USED TO CLAIM a negative differential was unreachable at M >= 0,
# on the grounds that the callback adds -M updrift and +M downdrift. That
# holds for the fillet's SIZE, not its slope. Falsified twice:
#
#   1. cell M40_beNA_f0.00 scores -0.713, MORE negative than the no-groin
#      baseline's -0.229 -- a groin at M = 40 driving the differential below
#      no groin at all;
#   2. the period-2 seed run inherits a 25.3 m fillet at t = 0 (the real
#      Buxton fillet is in the 2004 initial shoreline) and relaxes to 21.9 m,
#      scoring -0.168 m/yr at M = 50, f = 0.9. Negative, groin attached.
#
# So the observed 2004-2024 differential of -2.47 m/yr IS reachable, and is
# what a fillet relaxing after the 2003 storm damage looks like. Period 2
# failing to match it is informative rather than guaranteed, which is why it
# is now used out-of-sample rather than summed into the objective.
#
# The flag is kept, renamed to say what it really tests: whether the OBSERVED
# differential carries the sign of a groin still BUILDING its fillet across
# the window. Useful to report; not a statement about what the model can do.
OBSERVED_DIFFERENTIAL_IS_BUILDING = {
    p: OBSERVED_DIFFERENTIAL[p] >= 0.0 for p in PERIODS
}

# Old name kept so nothing importing it breaks. It no longer means 'the model
# cannot reach this'.
PERIOD_DIFFERENTIAL_IS_REACHABLE = OBSERVED_DIFFERENTIAL_IS_BUILDING


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
# SET FROM THE MODEL'S MEASURED NOISE FLOOR, NOT FROM ZERO.
#
# This was 1e-6, which assumes the model is deterministic. It is not. Running
# the SAME worker combination (2004 edgeBE, M = 50, f = 0.9, be1 = 50.3) twice
# on 2026-08-23 gave two rate curves differing by 7.7e-4 m/yr -- the same
# magnitude, and at the same domain, as the "drift" the guard was reporting
# against the published matrix run. Reference-vs-run2 (1.6e-4) came out
# SMALLER than run1-vs-run2, which is run-to-run scatter, not a code offset.
#
# WHERE THE SCATTER COMES FROM. The difference is float noise (median 1.8e-8
# across the reach) everywhere except a smooth bump centred on D16, which is
# the first domain OUTSIDE the beach-dune manager's footprint -- the 2022
# Buxton fill covers D6-D15, so D15/D16 is the sharpest discontinuity in the
# model. A threshold decision there turns last-bit differences into something
# a diffusion solver then spreads over ~10 domains.
#
# So a 1e-6 guard can never pass however well-synced the code is, and every
# sweep aborts on physics rather than on drift. 5e-3 sits ~6x above the
# measured floor and still catches real drift by orders of magnitude: a
# genuine divergence in build_cascade or the update loop moves rates by
# tenths of a m/yr, not thousandths.
#
# 7.4e-4 m/yr inside the fit window is ~1.5 cm of shoreline over 20 years,
# far below the resolution of a fit whose observed target cannot even
# reproduce the SHAPE of the groin pair (see observed_fillet_m).
VALIDATION_TOLERANCE_M_YR = 5e-3
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
    # Runs are filed <period>/<preset>/. The stem below already pins the
    # preset to edgeBE, so the directory pins it to the same thing.
    base = (PROJECT_BASE_DIR / "output" / "raw_runs"
            / f"{period}_{END_YEAR[period]}" / "edgeBE")
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
    """Where one period/preset sweep writes its results.

    Args:
        period: 1984 or 2004.
        preset: "edgeBE" or "zeroBE".

    Returns:
        The output directory Path.
    """
    stem = f"{period}_{END_YEAR[period]}_{preset}"
    return PROJECT_BASE_DIR / "output" / "groin_sweep" / stem


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

# =============================================================================
# FILLET SIZE -- the ranking metric
# =============================================================================
# WHY SIZE AND NOT SLOPE. `differential` (LRR[D6] - LRR[D5]) is algebraically
# the fillet's OLS SLOPE. The fillet saturates -- ~18 m at M = 40, by 1990 --
# so its slope over a 20-year window is near zero and carries almost no
# information about M, while its LEVEL carries all of it. Scored on slope the
# 1984 zeroBE pilot railed at M = 110, f = 1.0, both grid edges, at 115% of
# the reach sediment budget and still 35% short. Rescored on size, the SAME 43
# runs put the optimum at M = 50, f = 0.6 -- interior on both axes, 0.1 m
# error, and an f consistent with the 1996 repair and 2003 storm damage.


def measure_fillet(shoreline_m, baseline_m, geometry,
                   updrift_gis=None, downdrift_gis=None):
    """Fillet size at the end of a run, relative to a paired baseline.

    The fillet is `x_s[downdrift] - x_s[updrift]`: positive when the updrift
    domain sits seaward of the downdrift one, which is what a groin builds.
    Differenced against the M = 0 run at the SAME be1, so what is reported is
    the groin's contribution and not the coast's natural shape -- D5/D6 sit on
    a curving shoreline at Cape Point and carry an offset with no groin at all.

    Args:
        shoreline_m: [time x padded domain] array for the groin run.
        baseline_m: The same for its paired M = 0 run.
        geometry: DomainGeometry describing the padded array.
        updrift_gis, downdrift_gis: The groin pair; module defaults if None.

    Returns:
        Fillet size in metres at the final year, groin minus baseline.
    """
    up = geometry.gis_to_pad(GROIN_UPDRIFT_GIS if updrift_gis is None
                             else updrift_gis)
    down = geometry.gis_to_pad(GROIN_DOWNDRIFT_GIS if downdrift_gis is None
                               else downdrift_gis)
    run = float(shoreline_m[-1, down] - shoreline_m[-1, up])
    ref = float(baseline_m[-1, down] - baseline_m[-1, up])
    return run - ref


# THE GROIN IS A SUB-GRID FEATURE. Measured from the CoastSat transects on
# 2026-08-23, the fillet's observable width is ~190 m while one model domain
# is 500 m: the whole dipole fits inside a single cell.
#
# That is why this target is built from TRANSECTS and not from OBSERVED_LRR.
# The domain means average the fillet's peak against its own far field and
# what survives is the regional gradient, not the groin. In 2004-2024 that
# does not merely weaken the signal, it INVERTS it: the transect-scale dipole
# is +0.48 m/yr (groin-like, updrift accreting) while the domain-mean
# D6-D5 difference is -2.47 m/yr (anti-groin). Scored against the domain-mean
# version the sweep drove f to 0 -- "the groin stopped trapping" -- because
# the honest way to match a negative fillet is to trap nothing. Scored against
# the transects the same runs show the groin retaining ~44% of its period-1
# strength, which is what the imagery shows.
GROIN_TRANSECT_BOUNDARY = 76.5     # groin sits between transects ...0076 / ...0077
FILLET_HALFWIDTH_TRANSECTS = 6     # ~370 m each side at ~62 m transect spacing
FILLET_TREND_DOMAINS = (3, 9)      # window the regional trend is fitted over
FILLET_TREND_ORDER = 2


WETDRY_CHANGE_TABLE = (
    PROJECT_BASE_DIR / "hard-structures" / "groin" / "HAT-groin-test-output"
    / "shoreline_position_output" / "Change_from_wetdry_1967_D2_D12.csv")


def observed_fillet_m(period, trend_order=1):
    """Fillet change over one period, from the fixed 1967 wet/dry datum.

    THIRD TARGET IN TWO DAYS, and the reasons the first two failed are the
    reason this one is built from a different dataset entirely:

      1. Domain-mean LRR anomaly (-43.2 m for 2004-2024). The fillet is ~190 m
         wide inside a 500 m cell, so domain averaging blends its peak with its
         own far field.
      2. Transect LRR peak-to-trough (+22.1 / +9.7 m). `max(one side) -
         min(other side)` is biased POSITIVE: simulated at the observed scatter
         the bias alone is +27 m and +42 m, exceeding the signal it claimed to
         measure. It inverted period 2's sign.

    `Change_from_wetdry_1967_D2_D12.csv` avoids both. It is a purpose-built
    table of shoreline change against a FIXED 1967 datum, per domain, at 24
    dates from 1967 to 2023 -- so the fillet at D5/D6 is differenced directly,
    with no smoothing window to blend it and no order statistic to bias it.

    A TREND, NOT TWO ENDPOINTS. The year-to-year scatter is 11-14 m, so
    differencing the two end years throws away eight intermediate observations
    and inherits the noise of both. An OLS fit across the period gives
    +3.46 +/- 0.70 m/yr for 1984-2004 and -3.85 +/- 0.76 m/yr for 2004-2024 --
    4.9 and 5.1 sigma. Both periods carry a real, well-determined signal, which
    neither earlier target could show.

    SIGN matches `measure_fillet`: the table is landward-positive (+ = erosion),
    so D5 minus D6 is positive when the downdrift domain has retreated further
    than the updrift one -- which is what a groin builds.

    Args:
        period: 1984 or 2004.
        trend_order: Polynomial order for the fit across the period.

    Returns:
        Fillet CHANGE in metres over the period.

    Raises:
        FileNotFoundError: If the wet/dry change table is absent.
        ValueError: If a period has too few dated observations to fit.
    """
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
