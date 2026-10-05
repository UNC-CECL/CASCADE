#!/usr/bin/env python3
"""
One cell of the groin / background-erosion sweep, in its own process.

    python scripts/hatteras_ms/groin-sweep/HAT_groin_sweep_worker.py   # launched by HAT_groin_sweep.py with the cell's arguments

Builds and runs the model for one (period, preset, M, f, be1), scores it and
writes result.json with what it was built on. Details: scripts/hatteras_ms/README.md.

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-09-24
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
# parents[3]: this file is in hatteras_ms/groin-sweep/; the guard makes a move fail loudly here
PROJECT_BASE_DIR = next(_p for _p in _HERE.parents
                        if (_p / "pyproject.toml").exists())
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
# The pre-AST groin hook lives only in the sandbox copy
from cascade_pipeline.hindcast import build_cascade
from cascade_pipeline.run_registry import git_provenance, values_digest
from cascade_pipeline.shoreline import (build_shoreline_matrix,
                                        compute_change_rate, compute_lrr)

from site_layer.hatteras_site_config import (
    run_years,
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
    sweep_offset_mode,
    GROIN_SWEEP_ROOT,
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


# 0. which run this is -- parsed before anything derived from it

# The period cannot come from a module constant any more

_USAGE = (f"Usage: {Path(sys.argv[0]).name} <period> <preset> <M> <fraction> "
          f"<be1|none> <out_dir>\n"
          f"  period  one of {list(PERIODS)}\n"
          f"  preset  one of {list(PRESETS)}\n"
          f"  be1     background erosion at GIS 1 in m/yr, or 'none' to take "
          f"the site-config value for the period")


# Reads the six positional arguments, or exits with usage
def _parse_cli(argv):
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


# --- CONFIG ------------------------------------------------------------------
(START_YEAR, SOURCE_SINK_PRESET, _CLI_M, _CLI_FRACTION,
 _CLI_BE1, _CLI_OUT_DIR) = _parse_cli(sys.argv)


# 1. fixed configuration -- period, paths, scenario

# Resolved from the extractor, not pinned -- the same source the hindcast runner and the road setbacks use
from site_layer.hat_topo_version import topo_dirs, product_for_year  # scripts/, on path

# The product must be named, or a 1984 sweep builds on the 2004 island
TOPO_PRODUCT = product_for_year(START_YEAR)

_TOPO_DIR, _DUNE_DIR, TOPO_DUNE_VERSION = topo_dirs(TOPO_PRODUCT)

# No array-name constant here any more - see build_domain_file_paths below.

HATTERAS_DATA_BASE = PROJECT_BASE_DIR / "data" / "hatteras_init"
PARAMETER_FILE = "Hatteras-CASCADE-parameters.yaml"   # resolved by CASCADE
from site_layer.hat_topo_version import DOMAIN_ROOT as BARRIER3D_DIR  # noqa: E402
# Taken from what topo_dirs() RETURNED rather than re-joined from parts
from site_layer.hat_topo_version import BUFFER_DIR   # noqa: E402
DUNE_TOPO_DIR = _TOPO_DIR.parent

# Cascade resolves its parameter file relative to cwd, exactly as the runner does
os.chdir(PROJECT_BASE_DIR)

PERIOD = HATTERAS_PERIODS[START_YEAR]

# Continuous-window overrides (environment): unset, the period runs as configured
_end_override = os.environ.get("HAT_SWEEP_END_YEAR", "").strip()
END_YEAR = int(_end_override) if _end_override else PERIOD["end_year"]
CONTINUOUS_WINDOW = bool(_end_override) and END_YEAR != PERIOD["end_year"]

# Configured period: start..last_model_year (MODEL_YEARS.md); an override end stays an exclusive boundary
RUN_YEARS = END_YEAR - START_YEAR if CONTINUOUS_WINDOW else run_years(START_YEAR)
LAST_MODEL_YEAR = START_YEAR + RUN_YEARS - 1
SEA_LEVEL_RISE_RATE = PERIOD["sea_level_rise_rate"]
ISLAND_OFFSET_FILE = HATTERAS_DATA_BASE / PERIOD["island_offset_file"]

_storm_override = os.environ.get("HAT_SWEEP_STORM_FILE", "").strip()
STORM_FILE = (HATTERAS_DATA_BASE / _storm_override if _storm_override
              else HATTERAS_DATA_BASE / PERIOD["storm_file"])
ROAD_SETBACK_FILE = HATTERAS_DATA_BASE / PERIOD["road_setback_file"]

# SCENARIO = "full_management", written out so a matrix-table edit cannot reach the sweep
ENABLE_ROADWAY_MANAGEMENT = True
ENABLE_BEACH_DUNE_MANAGEMENT = True
ENABLE_NOURISHMENT_FILLS = True
ENABLE_HISTORICAL_ROAD_RELOCATIONS = False

# Asymmetric source/sink: 1.0 volume-neutral, 0.0 a pure updrift source
_sink_raw = os.environ.get("HAT_GROIN_SINK_FRACTION", "").strip()
SINK_FRACTION = float(_sink_raw) if _sink_raw else 1.0

NUM_CORES = 1        # >1 has crashed on this configuration; leave at 1
RMIN_VALUE, RMAX_VALUE = 0.55, 0.95
DUNE_DESIGN_ELEVATION_M = 3.0
DUNE_MINIMUM_ELEVATION_M = 0.0
ROAD_WIDTH = 20.0
# Hs IS A CALIBRATION VALUE, NOT AN OBSERVATION
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
# -----------------------------------------------------------------------------


# 2. input assembly (mirrors runner sections 2-6, prints suppressed)

# Domain file paths from cascade_pipeline, one definition for every reader
from cascade_pipeline.hindcast import build_domain_file_paths  # noqa: E402


# The island offset comes from cascade_pipeline.hindcast.build_island_offset, the runner's own builder
from cascade_pipeline.hindcast import build_island_offset  # noqa: E402


# Expands sparse per-GIS background-erosion rates onto the padded array
def build_background_erosion(be_rates, geometry):
    rates = [0.0] * geometry.total_domains
    for gis_id, rate in be_rates.items():
        pad_index = geometry.gis_to_pad(gis_id)
        if not 0 <= pad_index < geometry.total_domains:
            raise ValueError(f"GIS {gis_id} -> pad index {pad_index}, outside "
                             f"0-{geometry.total_domains - 1}")
        rates[pad_index] = float(rate)
    return rates


# Builds every per-domain forcing array the model needs
def assemble_forcing(be1):
    geometry = HATTERAS_DOMAINS
    # The PRODUCT is passed, exactly as the runner does at its section 3
    elevation_paths, dune_paths = build_domain_file_paths(
        geometry, TOPO_PRODUCT)

    missing = [p for p in elevation_paths + dune_paths if not Path(p).exists()]
    if missing:
        raise FileNotFoundError(
            f"{len(missing)} init files missing; first: {missing[0]}")

    # Source/sink

    # Built directly rather than through resolve_be_preset()
    if SOURCE_SINK_PRESET == "zeroBE":
        be_rates = {}
    else:
        be_rates = {
            1: float(be1) if be1 is not None else be_gis1_default(START_YEAR),
            90: be_gis90(START_YEAR),
        }
    background_erosion = build_background_erosion(be_rates, geometry)

    # Roadway
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

    # Beach/dune
    bn_schedule = nourishment.build_schedule(
        HATTERAS_NOURISHMENT_PROJECTS, geometry, START_YEAR, LAST_MODEL_YEAR)
    bn_schedule_applied = bn_schedule if ENABLE_NOURISHMENT_FILLS else (
        nourishment.build_schedule([], geometry, START_YEAR, LAST_MODEL_YEAR))
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
        island_offset=build_island_offset(ISLAND_OFFSET_FILE, geometry,
                                          mode=sweep_offset_mode()),
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


# 3. run_cascade_simulation

# build_cascade comes from cascade_pipeline; run_cascade_simulation below is different on purpose


# Steps a built Cascade through its period
def run_cascade_simulation(
    cascade, run_years, name, start_year, geometry,
    alongshore_section_count,
    historical_road_events=(), relocations_enabled=True, setback_check=None,
    nourishment_schedule=None, groin_callback=None,
):
    events = []

    for time_step in range(run_years):
        current_year = start_year + time_step

        if nourishment_schedule is not None:
            applied = nourishment_schedule.apply_to_cascade(
                cascade, current_year)
            for row in applied:
                # `row` already carries year
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
                # `row` carries its own kind, so the row wins
                events.append({"year": current_year, "run_name": name, **row})

        cascade.update()

        if getattr(cascade, "b3d_break", False):
            events.append(dict(kind="b3d_break", year=current_year,
                               time_step=time_step))
            break

    # The failure mode this guards is silent
    if groin_callback is not None and not groin_callback.year_TS:
        raise RuntimeError(
            "the groin callback was never called -- the pre-AST hook in "
            "cascade/cascade_groin.py is missing, so this run is identical "
            "to a no-groin run despite a groin being attached.")

    return cascade, events


# 4. scoring

# Scores one finished run against the CoastSat 1984-2004 target
def score_combo(model_lrr, geometry):
    model = {gis: float(model_lrr[geometry.gis_to_pad(gis)])
             for gis in FIT_DOMAINS_GIS}
    differential = model[GROIN_UPDRIFT_GIS] - model[GROIN_DOWNDRIFT_GIS]

    # IN A CONTINUOUS WINDOW THESE TARGETS DO NOT APPLY
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


# Extent is measured by the orchestrator, not here


# 4.5 construction lock

# Reading, then reopens it for writing

CONSTRUCT_LOCK = GROIN_SWEEP_ROOT / ".construct.lock"
CONSTRUCT_LOCK_TIMEOUT_S = 900
CONSTRUCT_LOCK_STALE_S = 300

# Windows refuses to unlink a file another process still has open
CONSTRUCT_LOCK_UNLINK_TRIES = 20
CONSTRUCT_LOCK_UNLINK_PAUSE_S = 0.1


# Removes the construction lock, tolerating Windows sharing errors
def _release_construct_lock():
    for _ in range(CONSTRUCT_LOCK_UNLINK_TRIES):
        try:
            CONSTRUCT_LOCK.unlink(missing_ok=True)
            return True
        except PermissionError:
            time.sleep(CONSTRUCT_LOCK_UNLINK_PAUSE_S)
        except OSError:
            return False
    return False


# Serialises Cascade construction across worker processes
@contextlib.contextmanager
def cascade_construction_lock():
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
                # A failed steal is not fatal -- another worker may be stealing the same lock this instant
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


PRISTINE_PARAMETERS = (GROIN_SWEEP_ROOT
                       / ".parameters_pristine.yaml")


# Repairs the shared parameter file if a previous run left it broken
def ensure_parameter_file_intact():
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


# `build_cascade` under the construction lock
def build_cascade_locked(**kwargs):
    with cascade_construction_lock():
        ensure_parameter_file_intact()
        return build_cascade(data_base=HATTERAS_DATA_BASE,
                             parameter_file=PARAMETER_FILE, **kwargs)


# 5. main

# Builds, runs and scores one combination
def run_combo(M, be1, fraction, out_dir):
    geometry = HATTERAS_DOMAINS
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    name = combo_dir_name(M, be1, fraction)

    forcing = assemble_forcing(be1)

    # Still be attached and still append diagnostics
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
        shoreline_offset=forcing["island_offset"],   # metres, as Cascade wants
        wave_height=Hs,
        wave_period=FIXED_WAVE_PERIOD,
        wave_asymmetry=FIXED_WAVE_ASYMMETRY,
        wave_angle_high_fraction=FIXED_WAVE_ANGLE_HIGH_FRACTION,
        berm_elevation=BERM_ELEVATION,
        MHW=MHW_ELEVATION,
        groin_callback=groin_callback,
    )

    # BRIE's diffusion number at the initial state, for the analytic fillet prediction
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

    # Both estimators, and the sweep is SCORED on the LRR
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

    # Road drowning happens inside CASCADE's roadway_manager, not through the historical-event list
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
        # What this cell was built on: named to match run_index.csv (README)
        topo_product=TOPO_PRODUCT,
        topo_dune_version=TOPO_DUNE_VERSION,
        be_values_digest=values_digest(forcing["be_rates"]),
        # The offset mode, since 2026-09-24
        offset_mode=sweep_offset_mode(),
        **{f"git_{key}": value
           for key, value in git_provenance(PROJECT_BASE_DIR).items()
           if key in ("commit", "dirty")},
        # The effective rate the groin actually applied, summed over the run and divided by its length
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
    # Nested dict flattened so the orchestrator's DataFrame gets one column per domain rather than a column of dicts
    result.update(result.pop("model_rates"))

    (out_dir / "result.json").write_text(json.dumps(result, indent=2))
    return result


# Run: build, run, score and write the cell
def main():
    # Argv was parsed and validated at import (section 0), because the period decides every path in section 1
    result = run_combo(_CLI_M, _CLI_BE1, _CLI_FRACTION, _CLI_OUT_DIR)
    print("RESULT_JSON=" + json.dumps(result))


if __name__ == "__main__":
    main()
