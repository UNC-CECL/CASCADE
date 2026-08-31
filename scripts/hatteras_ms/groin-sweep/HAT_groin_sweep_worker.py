#!/usr/bin/env python3
"""One combination of the groin / background-erosion sweep.

Runs a SINGLE (period, preset, M, deterioration_fraction, be1) combination in
its own fresh process, scores it against that period's CoastSat target, prints
one machine-readable line, and exits. Launched as a subprocess by
`HAT_groin_sweep.py` -- never in a loop inside one interpreter.

WHY A SUBPROCESS PER COMBINATION
    An earlier in-process sweep (all combinations back-to-back in one Python
    process) died partway through with a Windows access violation
    (0xC0000005): the signature of state accumulating across many repeated
    Cascade/Barrier3D constructions, not a bug in any single run. Process
    isolation makes the OS reclaim everything between combinations.

WHAT IS SWEPT, AND WHAT IS HELD FIXED
    Swept:  M (trapping rate) and f (deterioration floor), jointly across both
            periods, plus be1 in period 1 under edgeBE only.
    Fixed:  scenario    full_management (roadway + beach_dune + fills)
            Hs          2.5 m
            groin       GIS 6 updrift / GIS 5 downdrift, installed 1969,
                        linear_ramp deterioration, onset 1996, ramp 7 yr --
                        the SAME schedule in both periods, so period 2's
                        collapse to M*f emerges rather than being hard-coded
            be90        from HATTERAS_BE_RATES_EDGE, not pinned here
    See the orchestrator's docstring for why M and f need both periods.

WHAT IS SHARED, AND WHAT GUARDS THE REST
    `build_cascade` is imported from `cascade_pipeline.hindcast` -- the same
    definition the notebook and `HAT_hindcast_1984_2024.py` use, so a model
    built here is built exactly as a published run is. Only
    `run_cascade_simulation` below is this file's own, and it differs on
    purpose (see below).

    That difference is still guarded numerically, not by a source hash: the
    orchestrator reads the period's published matrix run's own M, f and be1
    out of its metadata, runs that combination through this worker, and
    refuses to report anything if the two rate curves differ. The reference's
    parameters are read rather than pinned, because the reference gets re-run
    as the fit is refined.

    Everything else that CAN come from the shared packages does: the presets,
    road events, nourishment projects and community zones from
    `hatteras_site_config`; the loaders, audits and shoreline maths from
    `cascade_pipeline`.

TWO DELIBERATE DIFFERENCES FROM THE HINDCAST RUNNER
    - `cascade.save()` is NOT called. It writes a ~140 MB npz per run; the
      344 cells of the full sweep would be ~48 GB of state nothing downstream
      reads. The shoreline matrix and the rate curve are written instead
      (~50 KB).
    - No figures, no GIFs, no run_index.csv row. A sweep combination is not a
      run of the matrix and must not be filed as one.

Usage:  HAT_groin_sweep_worker.py <period> <preset> <M> <fraction> <be1|none>
                                  <out_dir>
Prints: RESULT_JSON={...} on success. Non-zero exit means the combination
        failed; the orchestrator records it and moves on.

Author: Hannah A. Henry, UNC CECL
"""

from __future__ import annotations

import contextlib
import json
import os
import shutil
import sys
import time
from pathlib import Path

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
if not (PROJECT_BASE_DIR / "pyproject.toml").exists():
    raise RuntimeError(
        f"CASCADE repo root not found: {PROJECT_BASE_DIR} has no "
        f"pyproject.toml. This file expects to live in scripts/hatteras_ms/.")
if str(SCRIPTS_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPTS_DIR))

import numpy as np
import pandas as pd

from cascade.groin import GroinCallback, predict_fillet

from cascade_pipeline import nourishment, roadway
from cascade_pipeline import roadway as roadway_module
# The pre-AST groin hook lives only in the sandbox copy. Attaching a callback
# to stock Cascade silently does nothing, so which Cascade is imported is not
# optional -- cascade_pipeline.hindcast makes that choice once, for the
# notebook, the runner and this worker alike, via USE_SANDBOX_CASCADE.
from cascade_pipeline.hindcast import build_cascade
from cascade_pipeline.run_registry import git_provenance, values_digest
from cascade_pipeline.shoreline import (build_shoreline_matrix,
                                        compute_change_rate, compute_lrr)

from hatteras_site_config import (
    HATTERAS_BEACH_DUNE,
    HATTERAS_COMMUNITY_ZONES,
    HATTERAS_DOMAINS,
    HATTERAS_FIRST_ROAD_DOMAIN,
    HATTERAS_LAST_ROAD_DOMAIN,
    HATTERAS_NOURISHMENT_PROJECTS,
    HATTERAS_PERIODS,
    HATTERAS_RELOCATION_CHECK_2004,
    HATTERAS_ROAD_ELEVATION_FILE,
    HATTERAS_ROAD_EVENTS,
)

from HAT_groin_sweep_config import (
    DAM_TO_M,
    FIT_DOMAINS_GIS,
    FLIP_SIGN_MODEL,
    GROIN_DOWNDRIFT_GIS,
    GROIN_EXTENT_THRESHOLD_FRAC,
    GROIN_INSTALL_YEAR,
    GROIN_LAST_REPAIR_YEAR,
    GROIN_STORM_YEAR,
    GROIN_UPDRIFT_GIS,
    OBSERVED_LRR,
    OBSERVED_DIFFERENTIAL,
    PERIODS,
    PRESETS,
    be_gis1_default,
    be_gis90,
    combo_dir_name,
)


# =============================================================================
# 0. WHICH RUN THIS IS -- parsed before anything derived from it
# =============================================================================
# The period cannot come from a module constant any more: one worker serves
# both. It is parsed here, at the top, because every path and forcing file in
# section 1 is derived from it -- reading it later would mean two sources of
# truth for the same value for the length of the module.
#
# Parsed at import rather than inside main() for the same reason, and made
# fatal rather than defaulted: a worker that silently fell back to 1984 would
# fill a 2004 sweep directory with period-1 results that look entirely valid.

_USAGE = (f"Usage: {Path(sys.argv[0]).name} <period> <preset> <M> <fraction> "
          f"<be1|none> <out_dir>\n"
          f"  period  one of {list(PERIODS)}\n"
          f"  preset  one of {list(PRESETS)}\n"
          f"  be1     background erosion at GIS 1 in m/yr, or 'none' to take "
          f"the site-config value for the period")


def _parse_cli(argv):
    """Reads the six positional arguments, or exits with usage.

    Returns:
        (period, preset, M, fraction, be1, out_dir). be1 is None when the
        caller passed 'none', meaning "use the preset as it stands".
    """
    if len(argv) != 7:
        sys.exit(_USAGE)
    try:
        period = int(argv[1])
        preset = argv[2]
        M = float(argv[3])
        fraction = float(argv[4])
        be1 = None if argv[5].lower() == "none" else float(argv[5])
    except ValueError as exc:
        sys.exit(f"{_USAGE}\n\ncould not read arguments: {exc}")
    if period not in PERIODS:
        sys.exit(f"period must be one of {list(PERIODS)}, got {period}")
    if preset not in PRESETS:
        sys.exit(f"preset must be one of {list(PRESETS)}, got {preset!r}")
    if preset == "zeroBE" and be1 is not None:
        sys.exit("zeroBE has no background-erosion knob; pass be1 as 'none'")
    return period, preset, M, fraction, be1, argv[6]


(START_YEAR, SOURCE_SINK_PRESET, _CLI_M, _CLI_FRACTION,
 _CLI_BE1, _CLI_OUT_DIR) = _parse_cli(sys.argv)


# =============================================================================
# 1. FIXED CONFIGURATION -- period, paths, scenario
# =============================================================================

# Resolved from the extractor, not pinned -- the same source the hindcast
# runner and the road setbacks use, so the sweep cannot drift out of step with
# them. See scripts/hat_topo_version.py.
from hat_topo_version import topo_dirs, product_for_year  # scripts/, on path

# THE PRODUCT MUST BE NAMED. topo_dirs() with no argument resolves
# DEFAULT_PRODUCT, which is "2004-start" -- so a 1984 sweep built its model on
# the 2004 island while its reference matrix run used 1984-start. All 90
# domains differ between the two products and 65 differ in interior SHAPE, so
# the sweep was fitting M and f against a different barrier. The orchestrator's
# drift guard caught it on 2026-08-29 at 0.126 m/yr, 25x its tolerance, and
# blamed "section 3" -- which was in sync all along.
#
# This line used to carry a comment saying the bare call "still resolves
# 2004-start, which is what this sweep read before the 2026-08-25
# restructure". That was true and it was the bug: preserving pre-restructure
# behaviour is exactly wrong once the periods stopped sharing a topography.
TOPO_PRODUCT = product_for_year(START_YEAR)

_TOPO_DIR, _DUNE_DIR, TOPO_DUNE_VERSION = topo_dirs(TOPO_PRODUCT)

# No array-name constant here any more - see build_domain_file_paths below.

HATTERAS_DATA_BASE = PROJECT_BASE_DIR / "data" / "hatteras_init"
PARAMETER_FILE = "Hatteras-CASCADE-parameters.yaml"   # resolved by CASCADE
BARRIER3D_DIR = HATTERAS_DATA_BASE / "1-barrier3d-domains"
# Taken from what topo_dirs() RETURNED rather than re-joined from parts -
# re-joining is how a resolver gets bypassed without anyone noticing.
from hat_topo_version import BUFFER_DIR   # noqa: E402
DUNE_TOPO_DIR = _TOPO_DIR.parent

# Cascade resolves its parameter file relative to cwd, exactly as the runner
# does. Every path above is already absolute, so this only affects CASCADE.
os.chdir(PROJECT_BASE_DIR)

PERIOD = HATTERAS_PERIODS[START_YEAR]

# CONTINUOUS-WINDOW OVERRIDES. Unset, these change nothing and the worker runs
# exactly the period the site config defines. Set, they let ONE run span both
# hindcast periods -- 1984-2024 rather than 1984-2004 -- which is what makes M
# and f separable: the fillet builds in the first half and declines in the
# second, so the two knobs are constrained by different stretches of one curve
# instead of trading off inside a single 20-year window.
#
# They are environment variables and not a new HATTERAS_PERIODS entry on
# purpose. The site config is keyed by START YEAR, so a 1984-2024 entry would
# collide with 1984-2004, and every consumer that loops over HATTERAS_PERIODS
# -- the run matrix above all -- would silently gain a period it was never
# meant to run.
_end_override = os.environ.get("HAT_SWEEP_END_YEAR", "").strip()
END_YEAR = int(_end_override) if _end_override else PERIOD["end_year"]
CONTINUOUS_WINDOW = bool(_end_override) and END_YEAR != PERIOD["end_year"]

RUN_YEARS = END_YEAR - START_YEAR
SEA_LEVEL_RISE_RATE = PERIOD["sea_level_rise_rate"]
ISLAND_OFFSET_FILE = HATTERAS_DATA_BASE / PERIOD["island_offset_file"]

_storm_override = os.environ.get("HAT_SWEEP_STORM_FILE", "").strip()
STORM_FILE = (HATTERAS_DATA_BASE / _storm_override if _storm_override
              else HATTERAS_DATA_BASE / PERIOD["storm_file"])
ROAD_SETBACK_FILE = HATTERAS_DATA_BASE / PERIOD["road_setback_file"]

# SCENARIO = "full_management", expanded. Written out rather than looked up so
# a sweep combination cannot inherit a scenario-table edit made for the matrix.
ENABLE_ROADWAY_MANAGEMENT = True
ENABLE_BEACH_DUNE_MANAGEMENT = True
ENABLE_NOURISHMENT_FILLS = True
ENABLE_HISTORICAL_ROAD_RELOCATIONS = False

# ASYMMETRIC SOURCE/SINK. 1.0 is the volume-neutral pair the module has always
# used; 0.0 makes the groin a pure updrift source, the physical statement that
# the trapped sand comes from outside the pair rather than from the downdrift
# beach. Read from the environment so a variant is a separate sweep in its own
# directory, never mixed into the symmetric results.
_sink_raw = os.environ.get("HAT_GROIN_SINK_FRACTION", "").strip()
SINK_FRACTION = float(_sink_raw) if _sink_raw else 1.0

NUM_CORES = 1        # >1 has crashed on this configuration; leave at 1
RMIN_VALUE, RMAX_VALUE = 0.55, 0.95
DUNE_DESIGN_ELEVATION_M = 3.0
DUNE_MINIMUM_ELEVATION_M = 0.0
ROAD_WIDTH = 20.0
# Hs IS A CALIBRATION VALUE, NOT AN OBSERVATION. The runner records it as
# "calibration value 2.5" and notes that a sweep over it belongs in a separate
# script -- so overriding it here is using the knob as intended, not bending a
# measured constant.
#
# It matters for the groin because BRIE's alongshore diffusivity scales roughly
# as Hs^2.5, and diffusivity sets two things that pull in opposite directions:
#
#     fillet DECAY rate  ~ diffusivity        (too slow at Hs = 2.5: the
#                                              modelled fillet relaxes -15.3 m
#                                              over 1984-2024 against -24.4 m
#                                              observed, even with no groin)
#     fillet EXTENT      ~ sqrt(diffusivity)  (currently 2,000 m modelled
#                                              against 2,250 m observed, i.e.
#                                              already slightly narrow)
#
# So raising Hs speeds the decay toward the observations while widening the
# fillet away from them. The diagnostic exists to measure that trade rather
# than argue it. Unset, this changes nothing.
_hs_raw = os.environ.get("HAT_SWEEP_HS", "").strip()
Hs = float(_hs_raw) if _hs_raw else 2.5
FIXED_WAVE_PERIOD = 8
FIXED_WAVE_ASYMMETRY = 0.7
FIXED_WAVE_ANGLE_HIGH_FRACTION = 0.1
BERM_ELEVATION = 1.7     # m NAVD88
MHW_ELEVATION = 0.36     # m NAVD88
ENABLE_SANDBAG_PLACEMENT = False
SANDBAG_ELEVATION = 0
SEA_LEVEL_CONSTANT = True


# =============================================================================
# 2. INPUT ASSEMBLY  (mirrors runner sections 2-6, prints suppressed)
# =============================================================================

# build_domain_file_paths() used to be COPIED here from the hindcast runner,
# names and all. Both copies spelled f"domain_{gis_id}_topography_{init_year}"
# and both had to be edited in step. It now comes from cascade_pipeline, which
# itself delegates to hat_topo_version.domain_arrays() - one definition of the
# directory and one of the filename, for every reader in the repo.
from cascade_pipeline.hindcast import build_domain_file_paths  # noqa: E402


def load_island_offset_dam(offset_path, geometry):
    """Loads a padded BRIE island-offset file and converts it to decameters.

    Copied from the hindcast runner. See its section 3.1.

    Args:
        offset_path: Path to a single-column padded offset CSV, in meters.
        geometry: DomainGeometry the file must match in length.

    Returns:
        A 1-D array of offsets in decameters, one per padded domain.

    Raises:
        ValueError: If the file is not one value per padded domain.
    """
    offset_m = np.loadtxt(offset_path, skiprows=1, delimiter=",")
    if offset_m.ndim != 1 or offset_m.size != geometry.total_domains:
        raise ValueError(f"{offset_path.name}: expected "
                         f"{geometry.total_domains} values, got shape "
                         f"{offset_m.shape}")
    return offset_m / DAM_TO_M


def build_background_erosion(be_rates, geometry):
    """Expands sparse per-GIS background-erosion rates onto the padded array.

    Copied from the hindcast runner. See its section 4.3.

    Args:
        be_rates: Mapping of GIS domain ID to rate in m/yr. Domains absent
            from the mapping get 0.0.
        geometry: DomainGeometry describing the padded array.

    Returns:
        A list of geometry.total_domains rates, ready to pass to Cascade().

    Raises:
        ValueError: If a GIS ID falls outside the padded array.
    """
    rates = [0.0] * geometry.total_domains
    for gis_id, rate in be_rates.items():
        pad_index = geometry.gis_to_pad(gis_id)
        if not 0 <= pad_index < geometry.total_domains:
            raise ValueError(f"GIS {gis_id} -> pad index {pad_index}, outside "
                             f"0-{geometry.total_domains - 1}")
        rates[pad_index] = float(rate)
    return rates


def assemble_forcing(be1):
    """Builds every per-domain forcing array the model needs.

    Everything here is independent of the groin, so it is identical across
    combinations that share `be1`. It is rebuilt per process anyway -- the
    whole point of process isolation is that nothing is shared.

    Args:
        be1: Background-erosion rate at GIS 1, m/yr, or None to take the
            site-config value for this period. Ignored under zeroBE.

    Returns:
        A dict of the arrays build_cascade needs, plus the schedule and
        management masks the run loop and the audits use.
    """
    geometry = HATTERAS_DOMAINS
    # No product argument: this worker sweeps the DEFAULT_PRODUCT, the same one
    # topo_dirs() resolved above for DUNE_TOPO_DIR.
    # The PRODUCT is passed, exactly as the runner does at its
    # section 3. Omitting it silently selects DEFAULT_PRODUCT.
    elevation_paths, dune_paths = build_domain_file_paths(
        geometry, TOPO_PRODUCT)

    missing = [p for p in elevation_paths + dune_paths if not Path(p).exists()]
    if missing:
        raise FileNotFoundError(
            f"{len(missing)} init files missing; first: {missing[0]}")

    # --- source/sink --------------------------------------------------------
    # Built directly rather than through resolve_be_preset(): where be1 is
    # swept the sweep varies a value the presets hold fixed, so calling it
    # "edgeBE" would file a swept run under a preset name whose numbers it
    # does not use. The SHAPE is the preset's and nothing else --
    #   zeroBE  nothing imposed anywhere, so there is no knob to sweep and
    #           the groin has to account for the southern deficit alone.
    #   edgeBE  two end domains, interior untouched. be1 is swept in period 1
    #           and taken from the site config in period 2, where no prior fit
    #           stands behind a bracket.
    if SOURCE_SINK_PRESET == "zeroBE":
        be_rates = {}
    else:
        be_rates = {
            1: float(be1) if be1 is not None else be_gis1_default(START_YEAR),
            90: be_gis90(START_YEAR),
        }
    background_erosion = build_background_erosion(be_rates, geometry)

    # --- roadway ------------------------------------------------------------
    road_span = (HATTERAS_FIRST_ROAD_DOMAIN, HATTERAS_LAST_ROAD_DOMAIN)
    road_config = roadway.RoadwayConfig()
    road_setbacks_full, _ = roadway.load_road_setbacks(
        ROAD_SETBACK_FILE, geometry, *road_span)
    road_elevation_full, _ = roadway.load_road_elevations(
        HATTERAS_DATA_BASE / HATTERAS_ROAD_ELEVATION_FILE, geometry,
        *road_span, config=road_config)
    roadway_management_on = roadway.build_roadway_management_on(
        geometry, *road_span, community_zones=HATTERAS_COMMUNITY_ZONES,
        enabled=ENABLE_ROADWAY_MANAGEMENT)

    # --- beach/dune ---------------------------------------------------------
    bn_schedule = nourishment.build_schedule(
        HATTERAS_NOURISHMENT_PROJECTS, geometry, START_YEAR, END_YEAR)
    bn_schedule_applied = bn_schedule if ENABLE_NOURISHMENT_FILLS else (
        nourishment.build_schedule([], geometry, START_YEAR, END_YEAR))
    overwash_filter = (
        nourishment.build_overwash_filter(
            geometry, HATTERAS_COMMUNITY_ZONES, config=HATTERAS_BEACH_DUNE)
        if ENABLE_BEACH_DUNE_MANAGEMENT
        else [0.0] * geometry.total_domains)
    beach_dune_management_on = nourishment.build_beach_dune_management_on(
        geometry, HATTERAS_COMMUNITY_ZONES, bn_schedule.nourished_gis,
        enabled=ENABLE_BEACH_DUNE_MANAGEMENT)

    return dict(
        elevation_paths=elevation_paths,
        dune_paths=dune_paths,
        island_offset_dam=load_island_offset_dam(ISLAND_OFFSET_FILE, geometry),
        background_erosion=background_erosion,
        be_rates=be_rates,
        road_setbacks_full=road_setbacks_full,
        road_elevation_full=road_elevation_full,
        roadway_management_on=roadway_management_on,
        overwash_filter=overwash_filter,
        overwash_to_dune=HATTERAS_BEACH_DUNE.overwash_to_dune_pct,
        beach_dune_management_on=beach_dune_management_on,
        bn_schedule_applied=bn_schedule_applied,
    )


# =============================================================================
# 3. run_cascade_simulation
# =============================================================================
# `build_cascade` is NOT here: it is imported from cascade_pipeline.hindcast,
# the same definition the runner and the notebook build from, so the three can
# no longer disagree about how a model is constructed.
#
# `run_cascade_simulation` below is a deliberately different function, not a
# copy -- it collects the per-year events instead of printing them, writes no
# artifacts, and RAISES on a missing groin hook where the runner warns. The
# module docstring states why.


def run_cascade_simulation(
    cascade, run_years, name, start_year, geometry,
    alongshore_section_count,
    historical_road_events=(), relocations_enabled=True, setback_check=None,
    nourishment_schedule=None, groin_callback=None,
):
    """Steps a built Cascade through its period.

    Differs from the runner's copy only in what it writes: no artifacts are
    saved here (the caller writes the two small files it needs) and the
    per-year event log is collected rather than printed.

    Args:
        cascade: A Cascade from build_cascade, not yet stepped.
        run_years: Annual transitions to simulate.
        name: Run name, recorded on nourishment log rows.
        start_year: Calendar year of the run's first state.
        geometry: DomainGeometry, for GIS <-> pad translation.
        alongshore_section_count: Padded domain count.
        historical_road_events: RelocationEvent / BridgeEvent sequence.
        relocations_enabled: Global toggle for relocation events.
        setback_check: {gis: measured_setback_m} reported beside relocations.
        nourishment_schedule: A NourishmentSchedule, or None for no fills.
        groin_callback: The attached GroinCallback, for the drive check.

    Returns:
        A (cascade, events) tuple. `events` is a list of dicts recording the
        nourishment and roadway actions the loop took.

    Raises:
        RuntimeError: If a groin was attached but never called -- the
            pre-AST hook is missing and the run is silently a no-groin run.
    """
    events = []

    for time_step in range(run_years):
        current_year = start_year + time_step

        if nourishment_schedule is not None:
            applied = nourishment_schedule.apply_to_cascade(
                cascade, current_year)
            for row in applied:
                # `row` already carries `year` (apply_to_cascade puts it
                # there), so passing year=current_year as well is a duplicate
                # keyword and raises. The runner spreads the row without
                # re-supplying year; this now matches it.
                events.append(dict(kind="nourishment", run_name=name,
                                   time_step=time_step, **row))
        else:
            cascade.nourish_now = np.zeros(alongshore_section_count)

        for event in historical_road_events or ():
            if current_year != event.year:
                continue
            for row in roadway_module.apply_historical_event(
                    cascade, event, geometry,
                    relocations_enabled=relocations_enabled,
                    setback_check=setback_check):
                # `row` carries its own `kind` ("bridge", "relocation", ...) --
                # the runner switches on it. A dict literal lets the row win
                # rather than colliding, which `dict(kind=..., **row)` did.
                events.append({"year": current_year, "run_name": name, **row})

        cascade.update()

        if getattr(cascade, "b3d_break", False):
            events.append(dict(kind="b3d_break", year=current_year,
                               time_step=time_step))
            break

    # The failure mode this guards is silent: with the hook missing, a groin
    # run produces a valid-looking no-groin result and the sweep would fit M
    # against a model M never touched.
    if groin_callback is not None and not groin_callback.year_TS:
        raise RuntimeError(
            "the groin callback was never called -- the pre-AST hook in "
            "cascade/cascade_groin.py is missing, so this run is identical "
            "to a no-groin run despite a groin being attached.")

    return cascade, events


# =============================================================================
# 4. SCORING
# =============================================================================

def score_combo(model_lrr, geometry):
    """Scores one finished run against the CoastSat 1984-2004 target.

    Four metrics, all computed, one ranked on. The differential is the ranking
    metric because it is the only one that identifies M: over D1-D12 the
    profile RMSE moves 7% while M moves 4x, whereas the differential moves 5x.
    The profile numbers are diagnostics, and the extent is neither -- it is
    the pre-registered check from `cascade.groin`, reported so it stays a test
    rather than a target.

    Args:
        model_lrr: Padded model LRR array, m/yr, seaward-positive -- the
            OLS slope through every annual state, the same estimator
            OBSERVED_LRR is. NOT the endpoint difference.
        geometry: DomainGeometry describing the padded array.

    Returns:
        A dict of the four metrics plus the modelled per-domain rates over
        FIT_DOMAINS_GIS.
    """
    model = {gis: float(model_lrr[geometry.gis_to_pad(gis)])
             for gis in FIT_DOMAINS_GIS}
    differential = model[GROIN_UPDRIFT_GIS] - model[GROIN_DOWNDRIFT_GIS]

    # IN A CONTINUOUS WINDOW THESE TARGETS DO NOT APPLY. OBSERVED_LRR and
    # OBSERVED_DIFFERENTIAL are keyed by START_YEAR and describe that period's
    # own 20-year window. Scoring a 40-year run against them would silently
    # compare a 1984-2024 rate curve to a 1984-2004 target and write the result
    # into the CSV as if it meant something. The modelled rates are still
    # returned -- they are a property of the run, not of any target -- and the
    # continuous sweep scores the per-domain CHANGE PROFILE itself.
    #
    # NOTE FOR ANY FUTURE CALLER: HAT_groin_sweep.py treats a NaN
    # differential_err as "this cell failed", so the production orchestrator
    # must NOT be pointed at a continuous window. Use the continuous driver.
    if CONTINUOUS_WINDOW:
        nan = float("nan")
        return dict(
            differential_m_yr=differential, differential_err=nan,
            rmse_pair=nan, rmse_window=nan, bias_window=nan,
            model_rates={f"rate_D{g}": model[g] for g in FIT_DOMAINS_GIS},
        )

    observed = OBSERVED_LRR[START_YEAR]

    def rmse(domains):
        errors = [model[g] - observed[g] for g in domains]
        return float(np.sqrt(np.mean(np.square(errors))))

    def bias(domains):
        return float(np.mean([model[g] - observed[g] for g in domains]))

    pair = (GROIN_DOWNDRIFT_GIS, GROIN_UPDRIFT_GIS)
    return dict(
        differential_m_yr=differential,
        differential_err=abs(differential - OBSERVED_DIFFERENTIAL[START_YEAR]),
        rmse_pair=rmse(pair),
        rmse_window=rmse(FIT_DOMAINS_GIS),
        bias_window=bias(FIT_DOMAINS_GIS),
        model_rates={f"rate_D{g}": model[g] for g in FIT_DOMAINS_GIS},
    )


# Extent is measured by the orchestrator, not here: pairing a groin run with
# its M = 0 baseline is a property of the grid, and the baseline may not have
# finished when this combination does. See measure_groin_extent in
# HAT_groin_sweep_config.


# =============================================================================
# 4.5 CONSTRUCTION LOCK
# =============================================================================
# CASCADE's Cascade.__init__ reaches brie_coupler.initialize_equal, which calls
# set_yaml("TMAX", ...) -- and set_yaml opens the SHARED parameter file for
# reading, then reopens it for writing:
#
#     with open(file_name) as f: doc = full_load(f)
#     doc[var_name] = new_vals
#     with open(file_name, "w") as f: dump(doc, f)
#
# Every worker points at the same data/hatteras_init/Hatteras-CASCADE-
# parameters.yaml. Run two at once and one truncates the file while the other
# is reading it: full_load returns None and the next line raises
# "TypeError: 'NoneType' object does not support item assignment". It is
# invisible in a serial run and fires on essentially every cell in a parallel
# one -- 215 of 215 in the first attempt at this sweep.
#
# Serialising CONSTRUCTION alone fixes it. Construction is seconds; the
# 20-year update loop it precedes is minutes and stays fully parallel, so the
# pool keeps almost all of its speedup.
#
# The alternative -- giving each worker a private copy of the datadir -- was
# rejected as the bigger change: the yaml sits beside the init surfaces the
# runs read, and duplicating that per worker trades a lock for a much larger
# surface of things that can go stale.

CONSTRUCT_LOCK = PROJECT_BASE_DIR / "output" / "groin_sweep" / ".construct.lock"
CONSTRUCT_LOCK_TIMEOUT_S = 900
CONSTRUCT_LOCK_STALE_S = 300

# Windows refuses to unlink a file another process still has open, and both
# the release and the stale-steal below race exactly that. Unguarded, a single
# WinError 32 killed the worker AND left the lock behind, after which every
# later cell blocked on it -- 21 of 43 cells lost to one transient collision,
# and 17 hours of wall clock, on 2026-08-22. Retry briefly, then give up
# quietly: an orphan is self-healing via CONSTRUCT_LOCK_STALE_S, whereas an
# exception here is not.
CONSTRUCT_LOCK_UNLINK_TRIES = 20
CONSTRUCT_LOCK_UNLINK_PAUSE_S = 0.1


def _release_construct_lock():
    """Removes the construction lock, tolerating Windows sharing errors.

    Returns:
        True if the lock is gone, False if it could not be removed. False is
        not an error: the stale-steal path reclaims it after
        CONSTRUCT_LOCK_STALE_S, which is strictly better than raising and
        taking the worker down with the lock still held.
    """
    for _ in range(CONSTRUCT_LOCK_UNLINK_TRIES):
        try:
            CONSTRUCT_LOCK.unlink(missing_ok=True)
            return True
        except PermissionError:
            time.sleep(CONSTRUCT_LOCK_UNLINK_PAUSE_S)
        except OSError:
            return False
    return False


@contextlib.contextmanager
def cascade_construction_lock():
    """Serialises Cascade construction across worker processes.

    A stale lock is stolen rather than waited on: a worker killed mid-
    construction (the 0xC0000005 this design already expects) would otherwise
    block every later cell until the sweep's own timeout.

    Raises:
        RuntimeError: If the lock cannot be taken within the timeout.
    """
    CONSTRUCT_LOCK.parent.mkdir(parents=True, exist_ok=True)
    deadline = time.time() + CONSTRUCT_LOCK_TIMEOUT_S
    handle = None
    while handle is None:
        try:
            handle = os.open(str(CONSTRUCT_LOCK),
                             os.O_CREAT | os.O_EXCL | os.O_WRONLY)
        except FileExistsError:
            try:
                age = time.time() - CONSTRUCT_LOCK.stat().st_mtime
            except FileNotFoundError:
                continue          # released between the failure and the stat
            if age > CONSTRUCT_LOCK_STALE_S:
                # A failed steal is not fatal -- another worker may be
                # stealing the same lock this instant. Fall through and let
                # the loop retry rather than dying on the collision.
                _release_construct_lock()
                continue
            if time.time() > deadline:
                raise RuntimeError(
                    f"could not take the Cascade construction lock within "
                    f"{CONSTRUCT_LOCK_TIMEOUT_S}s ({CONSTRUCT_LOCK} held for "
                    f"{age:.0f}s). Delete it if no sweep is running.")
            time.sleep(0.25)

    try:
        os.write(handle, f"{os.getpid()} {time.time():.0f}".encode())
        yield
    finally:
        os.close(handle)
        _release_construct_lock()


PRISTINE_PARAMETERS = (PROJECT_BASE_DIR / "output" / "groin_sweep"
                       / ".parameters_pristine.yaml")


def ensure_parameter_file_intact():
    """Repairs the shared parameter file if a previous run left it broken.

    The lock stops two workers interleaving their writes, but it cannot help
    if a worker is KILLED between set_yaml's truncating open and its dump --
    and 0xC0000005 does exactly that on this configuration. The file is left
    half-written, and because it is shared, every remaining cell of the sweep
    then dies on a YAML scanner error. That is how a single transient crash
    turns into an overnight run with zero results.

    Every value CASCADE writes here (TMAX, the shoreface terms, the three
    file paths) is recomputed at the next construction, so restoring a
    snapshot loses nothing: only the hand-set parameters matter and those are
    identical in every version of this file.

    Must be called while holding the construction lock.
    """
    PRISTINE_PARAMETERS.parent.mkdir(parents=True, exist_ok=True)
    live = HATTERAS_DATA_BASE / PARAMETER_FILE

    try:
        import yaml
        parsed = yaml.full_load(live.read_text())
        readable = isinstance(parsed, dict) and bool(parsed)
    except Exception:
        readable = False

    if readable:
        # First healthy worker leaves a snapshot for the rest of the sweep.
        if not PRISTINE_PARAMETERS.exists():
            shutil.copy2(live, PRISTINE_PARAMETERS)
        return

    if not PRISTINE_PARAMETERS.exists():
        raise RuntimeError(
            f"{live} is unreadable and no pristine snapshot exists at "
            f"{PRISTINE_PARAMETERS}. Restore it from git "
            f"(git checkout -- {live.relative_to(PROJECT_BASE_DIR)}) before "
            f"sweeping.")
    shutil.copy2(PRISTINE_PARAMETERS, live)
    print(f"  repaired {live.name} from the pristine snapshot", flush=True)


def build_cascade_locked(**kwargs):
    """`build_cascade` under the construction lock.

    A wrapper rather than a lock inside `build_cascade`, so the shared
    definition in cascade_pipeline.hindcast stays free of sweep-only
    concurrency machinery.
    """
    with cascade_construction_lock():
        ensure_parameter_file_intact()
        return build_cascade(data_base=HATTERAS_DATA_BASE,
                             parameter_file=PARAMETER_FILE, **kwargs)


# =============================================================================
# 5. MAIN
# =============================================================================

def run_combo(M, be1, fraction, out_dir):
    """Builds, runs and scores one combination.

    Args:
        M: Groin trapping rate, m/yr. 0.0 attaches no groin at all -- the
            paired baseline column.
        be1: Background erosion at GIS 1, m/yr.
        fraction: Groin deterioration floor, in [0, 1].
        out_dir: Directory for this combination's two output files.

    Returns:
        A result dict, written to disk as result.json and printed as
        RESULT_JSON= for the orchestrator to parse.
    """
    geometry = HATTERAS_DOMAINS
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    name = combo_dir_name(M, be1, fraction)

    forcing = assemble_forcing(be1)

    # M = 0 means no groin, not a groin of zero strength: GroinCallback would
    # still be attached and still append diagnostics, and a "0 m/yr groin run"
    # in the record is a run someone will later mistake for a failed one.
    groin_callback = None
    if M > 0:
        groin_callback = GroinCallback(
            updrift_pad=geometry.gis_to_pad(GROIN_UPDRIFT_GIS),
            downdrift_pad=geometry.gis_to_pad(GROIN_DOWNDRIFT_GIS),
            trapping_rate_m_yr=M,
            start_year=START_YEAR,
            install_year=GROIN_INSTALL_YEAR,
            n_domains=geometry.total_domains,
            deterioration_delay_years=GROIN_LAST_REPAIR_YEAR - GROIN_INSTALL_YEAR,
            deterioration_mode="linear_ramp",
            deterioration_ramp_years=GROIN_STORM_YEAR - GROIN_LAST_REPAIR_YEAR,
            deterioration_fraction=fraction,
            sink_fraction=SINK_FRACTION,
        )

    total = geometry.total_domains
    t0 = time.perf_counter()
    cascade = build_cascade_locked(
        run_years=RUN_YEARS,
        name=name,
        storm_file=str(STORM_FILE),
        alongshore_section_count=total,
        num_cores=NUM_CORES,
        rmin=[RMIN_VALUE] * total, rmax=[RMAX_VALUE] * total,
        elevation_file=forcing["elevation_paths"],
        dune_file=forcing["dune_paths"],
        dune_design_elevation=[DUNE_DESIGN_ELEVATION_M] * total,
        dune_minimum_elevation=[DUNE_MINIMUM_ELEVATION_M] * total,
        road_ele=forcing["road_elevation_full"],
        road_width=ROAD_WIDTH,
        road_setback=forcing["road_setbacks_full"],
        overwash_filter=forcing["overwash_filter"],
        overwash_to_dune=forcing["overwash_to_dune"],
        nourishment_volume=[0.0] * total,
        background_erosion=forcing["background_erosion"],
        roadway_management_on=forcing["roadway_management_on"],
        beach_dune_manager_on=forcing["beach_dune_management_on"],
        sea_level_rise_rate=SEA_LEVEL_RISE_RATE,
        sea_level_constant=SEA_LEVEL_CONSTANT,
        sandbag_management_on=[ENABLE_SANDBAG_PLACEMENT] * total,
        sandbag_elevation=SANDBAG_ELEVATION,
        enable_shoreline_offset=True,
        shoreline_offset=forcing["island_offset_dam"],
        wave_height=Hs,
        wave_period=FIXED_WAVE_PERIOD,
        wave_asymmetry=FIXED_WAVE_ASYMMETRY,
        wave_angle_high_fraction=FIXED_WAVE_ANGLE_HIGH_FRACTION,
        berm_elevation=BERM_ELEVATION,
        MHW=MHW_ELEVATION,
        groin_callback=groin_callback,
    )

    # BRIE's diffusion number at the initial state, for the analytic fillet
    # prediction. Recorded per combination because it is what makes the
    # emergent extent checkable against something written down beforehand.
    brie = cascade._brie_coupler._brie
    r_ipl = float(brie._coast_diff[int(np.clip(round(90), 1, brie._wave_climl))]
                  * brie._dt / 2.0 / brie._dy ** 2)
    predicted_amplitude_m = predicted_extent_m = float("nan")
    if groin_callback is not None:
        predicted_amplitude_m, _, predicted_extent_m = predict_fillet(
            trapping_rate_m_yr=M, r_ipl=r_ipl, run_years=RUN_YEARS,
            dy_m=geometry.domain_spacing_m)

    cascade, events = run_cascade_simulation(
        cascade=cascade,
        run_years=RUN_YEARS,
        name=name,
        start_year=START_YEAR,
        geometry=geometry,
        alongshore_section_count=total,
        historical_road_events=HATTERAS_ROAD_EVENTS,
        relocations_enabled=ENABLE_HISTORICAL_ROAD_RELOCATIONS,
        setback_check=HATTERAS_RELOCATION_CHECK_2004,
        nourishment_schedule=forcing["bn_schedule_applied"],
        groin_callback=groin_callback,
    )
    run_seconds = time.perf_counter() - t0

    shoreline_m = build_shoreline_matrix(cascade)
    states = shoreline_m.shape[0]
    if states - 1 != RUN_YEARS:
        raise RuntimeError(
            f"span mismatch: {states} states implies {states - 1} years "
            f"elapsed, but RUN_YEARS is {RUN_YEARS} -- the run ended early "
            f"(b3d_break or drowning) and its rate denominator would be wrong.")

    # Both estimators, and the sweep is SCORED on the LRR. OBSERVED_LRR is a
    # per-transect OLS slope, so fitting M against an endpoint difference
    # would fit the groin to the gap between two estimators as well as to
    # the groin. The endpoint column is still written: the drift guard in
    # HAT_groin_sweep.py compares it against the published matrix run.
    change_rate = compute_change_rate(
        shoreline_m, span_years=RUN_YEARS, flip_sign=FLIP_SIGN_MODEL)
    model_lrr, model_lrr_r2 = compute_lrr(
        shoreline_m, span_years=RUN_YEARS, flip_sign=FLIP_SIGN_MODEL)

    np.save(out_dir / "shoreline_matrix.npy", shoreline_m)
    _real = slice(geometry.start_real_index, geometry.end_real_index)
    pd.DataFrame({
        "gis_domain": np.arange(geometry.first_gis_id, geometry.last_gis_id + 1),
        "change_rate_m_yr": change_rate[_real],
        "lrr_m_yr": model_lrr[_real],
        "lrr_r2": model_lrr_r2[_real],
    }).to_csv(out_dir / "shoreline_change_rate.csv", index=False)

    # Road drowning happens inside CASCADE's roadway_manager, not through the
    # historical-event list, so it has to be read back off the finished model
    # -- the same source the hindcast runner's section 12.2 uses. Counting
    # `events` here would report 0 for a run that drowned eight road_offset.
    road_summary = roadway.summarise_road_management(
        cascade, geometry, HATTERAS_FIRST_ROAD_DOMAIN,
        HATTERAS_LAST_ROAD_DOMAIN)

    result = dict(
        M=float(M),
        be1=None if be1 is None else float(be1),
        fraction=float(fraction),
        be90=(None if SOURCE_SINK_PRESET == "zeroBE"
              else float(be_gis90(START_YEAR))),
        period=int(START_YEAR),
        preset=SOURCE_SINK_PRESET,
        # WHAT THIS CELL WAS BUILT ON. Added 2026-08-31, and the reason is
        # specific: result.json recorded M, be1, f, period, preset and scores
        # and NOTHING about its inputs, so a sweep cell could not be told apart
        # from a cell of the same combination built on different topography.
        #
        # That is not hypothetical. The worker resolved topo_dirs() without
        # naming a product until 2026-08-30, so DEFAULT_PRODUCT ("2004-start")
        # answered and a 1984 sweep built on the 2004 barrier. When the
        # question "are these cells stale?" was asked on 2026-08-31, the
        # production runs answered in seconds from their own topo_product and
        # be_values_digest columns; the sweep could only be answered by
        # RE-RUNNING CELLS and diffing, at ~6 minutes each.
        #
        # These four fields are named to match run_index.csv exactly, so a
        # sweep cell and a matrix run can be compared without a translation
        # step.
        topo_product=TOPO_PRODUCT,
        topo_dune_version=TOPO_DUNE_VERSION,
        be_values_digest=values_digest(forcing["be_rates"]),
        **{f"git_{key}": value
           for key, value in git_provenance(PROJECT_BASE_DIR).items()
           if key in ("commit", "dirty")},
        # The effective rate the groin actually applied, summed over the run
        # and divided by its length. This is the quantity the two periods have
        # in common: period 1 delivers M*(16 + 4f)/20 and period 2 delivers
        # M*f. Recorded so the joint fit does not have to re-derive it from
        # the schedule, and so a period-2 row states plainly that M and f
        # reached it only as a product.
        mean_applied_M_m_yr=(
            float(np.mean(groin_callback.trapping_rate_applied_TS))
            if groin_callback is not None else 0.0),
        combo=name,
        run_seconds=round(run_seconds, 1),
        annual_states=int(states),
        r_ipl_t0=r_ipl,
        predicted_amplitude_m=predicted_amplitude_m,
        predicted_extent_m=predicted_extent_m,
        roads_drowned=sum(1 for r in road_summary if r["drowned"]),
        n_events=len(events),
        **score_combo(model_lrr, geometry),
    )
    # Nested dict flattened so the orchestrator's DataFrame gets one column
    # per domain rather than a column of dicts.
    result.update(result.pop("model_rates"))

    (out_dir / "result.json").write_text(json.dumps(result, indent=2))
    return result


def main():
    # argv was parsed and validated at import (section 0), because the period
    # decides every path in section 1. Nothing is re-read here.
    result = run_combo(_CLI_M, _CLI_BE1, _CLI_FRACTION, _CLI_OUT_DIR)
    print("RESULT_JSON=" + json.dumps(result))


if __name__ == "__main__":
    main()
