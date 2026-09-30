#!/usr/bin/env python3
"""
Hatteras Island CASCADE hindcast runner: the headless twin of HAT_hindcast_1984_2024.ipynb.

    python scripts/hatteras_ms/HAT_hindcast_1984_2024.py   # settings from hat_run.yaml, or HAT_* variables

Same sections, same order and same code as the notebook (the source of truth),
minus its display-only QC plots; the run is chosen in hat_run.yaml. Details: scripts/hatteras_ms/README.md.

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-09-30
"""

# 1. Imports

import datetime
import os
import shutil
import sys
import time
from pathlib import Path

# The packages live in scripts/, which is not installed: found from this file's location
_HERE = Path(__file__).resolve()
PROJECT_BASE_DIR = next(_p for _p in _HERE.parents
                        if (_p / "pyproject.toml").exists())
SCRIPTS_DIR = PROJECT_BASE_DIR / "scripts"
# HAT_hindcast_config sits beside this file; added explicitly so the notebook reaches it the same way
HATTERAS_MS_DIR = _HERE.parent
if not (PROJECT_BASE_DIR / "pyproject.toml").exists():
    raise RuntimeError(
        f"CASCADE repo root not found: {PROJECT_BASE_DIR} has no "
        f"pyproject.toml. This file expects to live in scripts/hatteras_ms/.")
for _path in (SCRIPTS_DIR, HATTERAS_MS_DIR):
    if str(_path) not in sys.path:
        sys.path.insert(0, str(_path))

# The settings file, read before anything it selects

# Read before the imports it selects: use_sandbox_cascade and show_figures are settled at import
from HAT_hindcast_config import load_run_config          # noqa: E402

_BOOT_CONFIG = load_run_config()
os.environ["CASCADE_USE_SANDBOX"] = (
    "1" if _BOOT_CONFIG.use_sandbox_cascade else "0")

import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt


from cascade.groin import BlockingGroinCallback, GroinCallback, predict_fillet

from cascade_pipeline import nourishment, roadway
from cascade_pipeline import reports
from cascade_pipeline.hindcast import (
    DAM_TO_M,
    USE_SANDBOX_CASCADE,
    brie_r_ipl,
    build_background_erosion,
    build_cascade,
    build_domain_file_paths,
    build_shoreline_target,
    build_target_table,
    load_barrier3d_contract,
    build_island_offset,
    island_offset_tilts,
    load_island_offset_dam,
    load_storm_series,
    measure_groin_extent,
    run_cascade_simulation,
    scenario_run_name,
    wave_climate_token,
    relocation_setback_token,
)
from cascade_pipeline.coastsat_lowess import (
    CoastSatDataset,
    LowessConfig,
    build_coastsat_series,
)
from cascade_pipeline.plotting.rate_comparison import (
    DEFAULT_RATE_COMPARISON,
    plot_annotated_rate_comparison,
    plot_rate_comparison,
)
from cascade_pipeline.plotting.shoreline_gif import GifConfig, make_all_shoreline_gifs
from cascade_pipeline.run_info import RunInfo
from cascade_pipeline.run_layout import resolve, write_path
from cascade_pipeline.run_registry import (
    INDEX_SECTION,
    barrier3d_provenance,
    rebuild_run_index,
    sweep_family,
    git_provenance,
    guard_run_dir,
    preset_dir_for,
    run_dir_contents,
    skill_vs_target,
    values_digest,
    timestamp,
    CALIBRATION_ARM,
    RUN_INDEX_FILENAME,
    write_run_metadata,
)
from cascade_pipeline.shoreline import (build_shoreline_matrix,
                                        compute_change_rate, compute_lrr)

from site_layer.hatteras_site_config import (
    HATTERAS_ANNOTATIONS,
    HATTERAS_BEACH_DUNE,
    HATTERAS_BE_PRESETS,
    be_rates,
    HATTERAS_BE_EDGE_DOMAINS,
    HATTERAS_BE_RATES_2004_IS_PLACEHOLDER,
    HATTERAS_COMMUNITY_ZONES,
    HATTERAS_DOMAINS,
    HATTERAS_FIRST_ROAD_DOMAIN,
    HATTERAS_LAST_ROAD_DOMAIN,
    HATTERAS_NOURISHMENT_PROJECTS,
    HATTERAS_PERIODS,
    HATTERAS_RELOCATION_CHECK_2004,
    HATTERAS_ROAD_ELEVATION_FILE,
    HATTERAS_ROAD_EVENTS,
    island_offset_version,
    resolve_be_preset,
    HATTERAS_GEOMETRY,
    HATTERAS_GEOMETRY_EXTENDED,
    HATTERAS_OFFSET_SOURCE,
    SCORE_INTERIOR_GIS,
)
from site_layer.hat_topo_version import DEFAULT_OFFSET_SOURCE  # noqa: E402
from cascade_pipeline.domains import DEFAULT_DOMAINS  # the surveyed reach, GIS 1-90

# Saved, not shown: a headless run cannot display (show_figures: null means False here)
SHOW_FIGURES = (False if _BOOT_CONFIG.show_figures is None
                else _BOOT_CONFIG.show_figures)
if not SHOW_FIGURES:
    matplotlib.use("Agg")

try:
    from tqdm.auto import tqdm
except ImportError:                                     # optional dependency
    tqdm = None

print(f"Imports OK from {SCRIPTS_DIR}")
print(f"USE_SANDBOX_CASCADE = {USE_SANDBOX_CASCADE}")
print(f"HATTERAS_DOMAINS.total_domains = {HATTERAS_DOMAINS.total_domains}  "
      f"(geometry {HATTERAS_GEOMETRY!r}: GIS {HATTERAS_DOMAINS.first_gis_id} "
      f"to {HATTERAS_DOMAINS.last_gis_id})")


# 2. Dune/topo -- per period

# The topography product is per period (HATTERAS_PERIODS); its version is resolved, never pinned

# 2.1 Project paths

# No per-period array names here: build_domain_file_paths() delegates to the resolver

# Which extraction: resolved by topo_dirs(); an old one only by override, e.g. topo_dirs("2004-start", override="v3")
from site_layer.hat_topo_version import topo_dirs, current_topo_versions  # scripts/, on sys.path above
from site_layer.hat_topo_version import BUFFER_DIR as _BUFFER_DIR

# _BOOT_CONFIG, not RUN_CONFIG: the period is needed before section 3; checked again there
TOPO_PRODUCT = HATTERAS_PERIODS[_BOOT_CONFIG.start_year]["topo_product"]

_TOPO_DIR, _DUNE_DIR, TOPO_DUNE_VERSION = topo_dirs(TOPO_PRODUCT)

print(f"topography            {TOPO_PRODUCT} / {TOPO_DUNE_VERSION}  "
      f"(product from the period, version resolved)")

HATTERAS_DATA_BASE = PROJECT_BASE_DIR / "data" / "hatteras_init"
OUTPUT_ROOT = PROJECT_BASE_DIR / "output" / "raw_runs"
# Where the rate fits live is hat_observed_rates' to say
from site_layer.hat_observed_rates import COASTSAT_LRR_ROOT as COASTSAT_BASE_DIR  # noqa: E402
PARAMETER_FILE = "Hatteras-CASCADE-parameters.yaml"  # resolved by CASCADE

from site_layer.hat_topo_version import DOMAIN_ROOT as BARRIER3D_DIR  # noqa: E402
# Taken from what topo_dirs() returned, never re-joined from parts
DUNE_TOPO_DIR = _TOPO_DIR.parent
BUFFER_DIR = _BUFFER_DIR

os.chdir(PROJECT_BASE_DIR)
# Runs are filed per period (OUTPUT_BASE_DIR, section 3); run_index.csv stays at the root
OUTPUT_ROOT.mkdir(parents=True, exist_ok=True)

reports.path_inventory([
    ("PROJECT_BASE_DIR", PROJECT_BASE_DIR),
    ("HATTERAS_DATA_BASE", HATTERAS_DATA_BASE),
    ("DUNE_TOPO_DIR", DUNE_TOPO_DIR),
    ("BUFFER_DIR", BUFFER_DIR),
    ("COASTSAT_BASE_DIR", COASTSAT_BASE_DIR),
    ("OUTPUT_ROOT", OUTPUT_ROOT),
])


# 2.2 Build the padded file lists


# The product is passed; hat_topo_version resolves directory and filename together
ELEVATION_FILE_PATHS, DUNE_FILE_PATHS = build_domain_file_paths(
    HATTERAS_DOMAINS, TOPO_PRODUCT)

print(f"\n{len(ELEVATION_FILE_PATHS)} elevation + {len(DUNE_FILE_PATHS)} dune "
      f"paths (expect {HATTERAS_DOMAINS.total_domains} each)")


# 2.3 Verify every file exists

# A stale version or a moved data folder fails here, with a count and the first offender

_expected_files = 2 * HATTERAS_DOMAINS.total_domains
_missing = [path for path in ELEVATION_FILE_PATHS + DUNE_FILE_PATHS
            if not Path(path).exists()]

if _missing:
    raise FileNotFoundError(
        f"{len(_missing)} of {_expected_files} init files missing. Check "
        f"TOPO_DUNE_VERSION ({TOPO_DUNE_VERSION!r}) and DUNE_TOPO_DIR.\n"
        f"  First missing: {_missing[0]}")

print(f"All {_expected_files} init files present.")


# 2.4 Units check against Barrier3D's input contract

# The .npy arrays must already be decameters MHW (a metres file runs 10x too tall); checked raw


BARRIER3D_CONTRACT = load_barrier3d_contract(HATTERAS_DATA_BASE / PARAMETER_FILE)

print(f"\nContract from {PARAMETER_FILE}:")
print(f"  BarrierLength -> {BARRIER3D_CONTRACT['barrier_length_cells']} "
      f"alongshore cells")
print(f"  MHW           -> {BARRIER3D_CONTRACT['mhw_dam']:.3f} dam")
print(f"  BermEl        -> {BARRIER3D_CONTRACT['berm_el_dam']:.3f} dam "
      f"above MHW\n")

reports.run_units_check(ELEVATION_FILE_PATHS, DUNE_FILE_PATHS,
                        BARRIER3D_CONTRACT, HATTERAS_DOMAINS)


# 3. Island orientation -- set START_YEAR

# START_YEAR selects the period; values come from hat_run.yaml, re-read here, and can be typed over

from HAT_hindcast_config import (                     # noqa: E402
    load_run_config, describe as _describe_run_config, preflight as _preflight,
    field_default as _field_default)

RUN_CONFIG = load_run_config()

# The two values section 1 already spent: a yaml edited since then raises
_BOOT_DRIFT = {
    name: (getattr(_BOOT_CONFIG, name), getattr(RUN_CONFIG, name))
    for name in ("use_sandbox_cascade", "show_figures")
    if getattr(_BOOT_CONFIG, name) != getattr(RUN_CONFIG, name)
}
if _BOOT_DRIFT:
    raise RuntimeError(
        "hat_run.yaml changed after section 1 imported on the old values:\n"
        + "\n".join(f"  {name}: section 1 used {was!r}, the file now says "
                    f"{now!r}" for name, (was, now) in _BOOT_DRIFT.items())
        + "\nThese are settled at import. Re-run section 1 (in a notebook, "
          "restart the kernel and Run All).")

# --- CONFIG ------------------------------------------------------------------
START_YEAR = RUN_CONFIG.start_year   # 1984, 1996, 2004 or 2010

# Source/sink preset: "zeroBE" nothing, "edgeBE" the two end domains, "calibBE" the full fit
SOURCE_SINK_PRESET = RUN_CONFIG.source_sink_preset

# Scenario -- the management combination this run simulates

# Named scenarios: each sets the roadway, beach_dune, fills and relocations switches
SCENARIOS = {
    # nothing human acts on the island: the counterfactual
    "natural": dict(roadway=False, beach_dune=False,
                    fills=False, relocations=False),
    # NC-12 defended, villages left to behave naturally
    "roadway_only": dict(roadway=True, beach_dune=False,
                         fills=False, relocations=False),
    # villages managed and nourished, the road passive
    "beachdune_only": dict(roadway=False, beach_dune=True,
                           fills=True, relocations=False),
    # everything: the status-quo hindcast
    "full_management": dict(roadway=True, beach_dune=True,
                            fills=True, relocations=False),
    # everything except the sand -- isolates the fills against full_management
    "full_no_fill": dict(roadway=True, beach_dune=True,
                         fills=False, relocations=False),
}

SCENARIO = RUN_CONFIG.scenario

# Read here: the run-name preview below depends on it
OFFSET_MODE = RUN_CONFIG.offset_mode

# The reach: must match the HAT_GEOMETRY hatteras_site_config was built from ("base" = GIS 1-90)
GEOMETRY = RUN_CONFIG.geometry
if GEOMETRY != HATTERAS_GEOMETRY:
    raise SystemExit(
        f"\n[stop] geometry mismatch: the config says {GEOMETRY!r} but "
        f"hatteras_site_config built {HATTERAS_GEOMETRY!r} from the "
        f"environment. Set HAT_GEOMETRY={GEOMETRY} in the environment "
        f"(the yaml alone cannot select a reach).\n")


# The groin is not part of the scenario: every scenario runs with and without it (12.3)
GROIN_ENABLED = RUN_CONFIG.groin_enabled

# Which groin: "dipole" (+/-M a year) or "blocking" (a fraction b of the transport); own name tokens
GROIN_KIND = RUN_CONFIG.groin_kind
if GROIN_KIND not in ("dipole", "blocking"):
    raise ValueError(f"groin kind {GROIN_KIND!r} must be 'dipole' or 'blocking'")
GROIN_TOKEN = (("groin" if GROIN_KIND == "dipole" else "groinblock")
               if GROIN_ENABLED else "nogroin")

# Where a relocated road is rebuilt, metres behind the dune line; "measured" = the t=0 setback
RELOCATION_SETBACK_M = RUN_CONFIG.relocation_setback_m
# -----------------------------------------------------------------------------

# Sensitivity tokens

# Forcing tokens keep a sensitivity cell out of the matrix run's name; None at the calibration values
_WAVE_VALUES = {
    "hs": RUN_CONFIG.hs,
    "wave_period_s": RUN_CONFIG.wave_period_s,
    "wave_asymmetry": RUN_CONFIG.wave_asymmetry,
    "wave_angle_high_fraction": RUN_CONFIG.wave_angle_high_fraction,
}
_WAVE_DEFAULTS = {name: _field_default(name) for name in _WAVE_VALUES}
# A sensitivity cell is its baseline's name plus this token (filed under sensitivity/<axis>/)
WAVE_TOKEN = wave_climate_token(_WAVE_VALUES, _WAVE_DEFAULTS)
RELOCATION_SETBACK_TOKEN = relocation_setback_token(
    RELOCATION_SETBACK_M, _field_default("relocation_setback_m"))

print("\n" + _describe_run_config())

# Expand
if SCENARIO not in SCENARIOS:
    raise ValueError(f"SCENARIO must be one of {sorted(SCENARIOS)}, "
                     f"got {SCENARIO!r}")
_SCENARIO_PRESET = SCENARIOS[SCENARIO]
ENABLE_ROADWAY_MANAGEMENT = _SCENARIO_PRESET["roadway"]
ENABLE_BEACH_DUNE_MANAGEMENT = _SCENARIO_PRESET["beach_dune"]
ENABLE_NOURISHMENT_FILLS = _SCENARIO_PRESET["fills"]
# HAT_RELOCATIONS overrides the scenario when set; None leaves the preset in charge
ENABLE_HISTORICAL_ROAD_RELOCATIONS = (
    _SCENARIO_PRESET["relocations"] if RUN_CONFIG.relocations is None
    else RUN_CONFIG.relocations)

# One-off overrides

# Uncomment to depart from the named scenario for one run; the departure is detected and named below
# ENABLE_ROADWAY_MANAGEMENT = False
# ENABLE_BEACH_DUNE_MANAGEMENT = False
# ENABLE_NOURISHMENT_FILLS = False
# ENABLE_HISTORICAL_ROAD_RELOCATIONS = True

_SCENARIO_DEPARTURES = {
    key: (want, got) for key, want, got in (
        ("roadway", _SCENARIO_PRESET["roadway"], ENABLE_ROADWAY_MANAGEMENT),
        ("beach_dune", _SCENARIO_PRESET["beach_dune"],
         ENABLE_BEACH_DUNE_MANAGEMENT),
        ("fills", _SCENARIO_PRESET["fills"], ENABLE_NOURISHMENT_FILLS),
        ("relocations", _SCENARIO_PRESET["relocations"],
         ENABLE_HISTORICAL_ROAD_RELOCATIONS),
    ) if want != got
}
if START_YEAR not in HATTERAS_PERIODS:
    raise ValueError(f"START_YEAR must be one of {sorted(HATTERAS_PERIODS)}, "
                     f"got {START_YEAR}")

# Normalised to the canonical preset name, so an alias never reaches RUN_NAME
SOURCE_SINK_PRESET, _PRESET_BY_PERIOD = resolve_be_preset(SOURCE_SINK_PRESET)

PERIOD = HATTERAS_PERIODS[START_YEAR]

# Section 2 chose the topography from the boot config: checked against the period that runs
if PERIOD["topo_product"] != TOPO_PRODUCT:
    raise SystemExit(
        f"\n[stop] topography product mismatch.\n"
        f"  section 2 loaded : {TOPO_PRODUCT}  (from boot start_year "
        f"{_BOOT_CONFIG.start_year})\n"
        f"  this period wants: {PERIOD['topo_product']}  (START_YEAR "
        f"{START_YEAR})\n"
        f"The domain arrays already in memory are the wrong ones. Re-run "
        f"with a consistent start_year.\n")
END_YEAR = PERIOD["end_year"]
RUN_YEARS = END_YEAR - START_YEAR

SEA_LEVEL_RISE_RATE = PERIOD["sea_level_rise_rate"]
ENABLE_NOURISHMENT = PERIOD["enable_nourishment"]
NOURISHMENT_VOLUME = PERIOD["nourishment_volume"]
ISLAND_OFFSET_FILE = HATTERAS_DATA_BASE / PERIOD["island_offset_file"]
# Which version folder that resolved to (v1, v2, or flat); see the json write.
ISLAND_OFFSET_VERSION = island_offset_version(START_YEAR)
STORM_FILE = HATTERAS_DATA_BASE / PERIOD["storm_file"]
ROAD_SETBACK_FILE = HATTERAS_DATA_BASE / PERIOD["road_setback_file"]

# Resolve the combinations that cannot both be true

# Contradictory switches (fill without a manager, relocation without a road) are resolved and announced
_FILLS_FORCED_OFF = (ENABLE_NOURISHMENT_FILLS
                     and not ENABLE_BEACH_DUNE_MANAGEMENT)
if _FILLS_FORCED_OFF:
    ENABLE_NOURISHMENT_FILLS = False
_RELOCATIONS_FORCED_OFF = (ENABLE_HISTORICAL_ROAD_RELOCATIONS
                           and not ENABLE_ROADWAY_MANAGEMENT)
if _RELOCATIONS_FORCED_OFF:
    ENABLE_HISTORICAL_ROAD_RELOCATIONS = False

# Run name: the period stem; the scenario suffix is derived from the switches in 7.5
RUN_NAME_STEM = f"HAT_{START_YEAR}_{END_YEAR}"

# Run directories are filed by period, so every later path (RUN_DIR, the 12.3 baseline) is too
PERIOD_TAG = f"{START_YEAR}_{END_YEAR}"
# Filed by purpose: matrix/, sensitivity/<axis>/, experiments/<tag>/, versions/<tag>/ (run_registry)
RUN_KIND = RUN_CONFIG.run_kind
RUN_TAG = RUN_CONFIG.run_tag
# HAT_ARM_TAG, the pre-09-16 spelling, is read as an experiment tag
_LEGACY_ARM = os.environ.get("HAT_ARM_TAG", "").strip()
if _LEGACY_ARM and not RUN_TAG:
    print(f"HAT_ARM_TAG={_LEGACY_ARM!r} is the pre-2026-09-16 spelling; read as "
          f"HAT_RUN_KIND=experiment HAT_RUN_TAG={_LEGACY_ARM!r}")
    RUN_KIND, RUN_TAG = "experiment", _LEGACY_ARM
# A run off the calibration waves cannot be filed in matrix/: refused
if RUN_KIND == "matrix" and WAVE_TOKEN:
    raise ValueError(
        f"wave climate is off calibration ({WAVE_TOKEN}) but HAT_RUN_KIND is "
        f"matrix. A forced run is a sensitivity cell (HAT_RUN_KIND=sensitivity) "
        f"or an experiment (HAT_RUN_KIND=experiment HAT_RUN_TAG=<name>).")
# Nor can a run on a non-default offset source (its name would collide): refused
if RUN_KIND == "matrix" and HATTERAS_OFFSET_SOURCE != DEFAULT_OFFSET_SOURCE:
    raise ValueError(
        f"island offset source is {HATTERAS_OFFSET_SOURCE!r} (not "
        f"{DEFAULT_OFFSET_SOURCE!r}) but HAT_RUN_KIND is matrix. The run name "
        f"carries no offset token, so this would overwrite the matrix run it is "
        f"meant to be compared against. Use HAT_RUN_KIND=experiment with "
        f"HAT_RUN_TAG=<name>.")
if RUN_KIND == "sensitivity" and not (WAVE_TOKEN or RELOCATION_SETBACK_TOKEN):
    raise ValueError(
        "HAT_RUN_KIND=sensitivity but every swept forcing is at its default; "
        "this cell would be the matrix run under another name.")
if RUN_KIND == "sensitivity" and not RUN_TAG:
    RUN_TAG = sweep_family(f"x_{WAVE_TOKEN or RELOCATION_SETBACK_TOKEN}")
# Model state (~99% of a run's size) kept for matrix runs unless output.save_model_state says otherwise
SAVE_MODEL_STATE = (RUN_CONFIG.save_model_state
                    if (RUN_KIND == "matrix"
                        or RUN_CONFIG.origins["save_model_state"] != "default")
                    else False)
# Built by run_registry, so the writer and the readers agree on the layout
OUTPUT_BASE_DIR = preset_dir_for(OUTPUT_ROOT, PERIOD_TAG, SOURCE_SINK_PRESET,
                                 kind=RUN_KIND, tag=RUN_TAG)
OUTPUT_BASE_DIR.mkdir(parents=True, exist_ok=True)

print(f"\nSTART_YEAR = {START_YEAR}  ->  {START_YEAR}-{END_YEAR}, "
      f"{RUN_YEARS} model years")
print(f"RUN_NAME_STEM = {RUN_NAME_STEM!r}"
      "   (scenario suffix derived in 7.5)")
print(f"SOURCE_SINK_PRESET = {SOURCE_SINK_PRESET!r}")
print(f"OUTPUT_BASE_DIR = {OUTPUT_BASE_DIR}")
print(f"RUN_KIND = {RUN_KIND!r}   RUN_TAG = {RUN_TAG!r}   SAVE_MODEL_STATE = {SAVE_MODEL_STATE}")

# Which Barrier3D: its checked-out branch is the model; recorded, warned about if not the fix branch
_B3D = barrier3d_provenance()
import cascade.beach_dune_manager as _bdm_module  # noqa: E402  (dune-cap provenance)
print(f"BARRIER3D = {_B3D['branch']}@{str(_B3D['commit'])[:7]}   "
      f"route_overwash fix: {_B3D['route_overwash_fix']}"
      + ("   (tree dirty)" if _B3D["dirty"] else ""))
print(f"            overwash gap/momentum fixes: {_B3D.get('gap_momentum_fix')}   "
      f"per-cell dune ceilings: {_B3D.get('per_cell_ceiling')}")
# The adopted model needs per-cell dune ceilings: a Barrier3D without them is refused
if _B3D.get("per_cell_ceiling") is not True:
    raise RuntimeError(
        f"Barrier3D {_B3D['branch']}@{str(_B3D['commit'])[:7]} has no per-cell dune "
        "ceilings, which the parameter template asks for. Check out hatteras/adopted "
        "in the Barrier3D repository.")
if _B3D.get("gap_momentum_fix") is not True:
    print("  WARNING: this Barrier3D does not carry the 2026-09-28 overwash gap and "
          "momentum fixes (hatteras/adopted has them).")
if _B3D["route_overwash_fix"] is not True:
    print("  WARNING: this Barrier3D does NOT carry the route_overwash index fix. "
          "Runs read the wrong cells in overwash and can crash silently. Check out "
          "fix/route-overwash-axis-swap in the Barrier3D repository unless this run "
          "reproduces an old one on purpose.")

# The name this scenario will produce, predicted from the switches

# Advisory preview of the run name; 7.5 derives the real one and raises if they differ
_PERIOD_HAS_FILL = bool(nourishment.build_schedule(
    HATTERAS_NOURISHMENT_PROJECTS, HATTERAS_DOMAINS,
    START_YEAR, END_YEAR).projects)
_PREVIEW_TOKENS = [
    SOURCE_SINK_PRESET,
    None if OFFSET_MODE == "asrun" else f"offset{OFFSET_MODE}",
    "road" if ENABLE_ROADWAY_MANAGEMENT else "noroad",
    "reloc" if ENABLE_HISTORICAL_ROAD_RELOCATIONS else None,
    "bdm" if ENABLE_BEACH_DUNE_MANAGEMENT else "nobdm",
    ("nourish" if (ENABLE_NOURISHMENT_FILLS and _PERIOD_HAS_FILL)
     else ("nonourish" if _PERIOD_HAS_FILL and ENABLE_BEACH_DUNE_MANAGEMENT
           else None)),
    GROIN_TOKEN,
    RELOCATION_SETBACK_TOKEN,
    WAVE_TOKEN,
]
RUN_NAME_PREVIEW = (f"{RUN_NAME_STEM}_"
                    + "_".join(t for t in _PREVIEW_TOKENS if t))

# Name, destination, collision and rough run time, reported before sections 4-10 run
print("\n" + _preflight(RUN_NAME_PREVIEW,
                        OUTPUT_BASE_DIR / RUN_NAME_PREVIEW,
                        config=RUN_CONFIG))

reports.scenario_report(
    scenario=SCENARIO, scenarios=SCENARIOS, departures=_SCENARIO_DEPARTURES,
    roadway_on=ENABLE_ROADWAY_MANAGEMENT,
    relocations_on=ENABLE_HISTORICAL_ROAD_RELOCATIONS,
    relocations_forced_off=_RELOCATIONS_FORCED_OFF,
    beach_dune_on=ENABLE_BEACH_DUNE_MANAGEMENT,
    fills_on=ENABLE_NOURISHMENT_FILLS, fills_forced_off=_FILLS_FORCED_OFF,
    period_expects_nourishment=ENABLE_NOURISHMENT,
    groin_enabled=GROIN_ENABLED, run_name_preview=RUN_NAME_PREVIEW,
    input_files=[("island offset", ISLAND_OFFSET_FILE),
                 ("storms", STORM_FILE),
                 ("road setback", ROAD_SETBACK_FILE)])


# 3.1 Island offsets

# Each padded domain's cross-shore starting position: shoreline_offset on the Cascade() call


# OFFSET_MODE picks the variant: "metres" (default since 2026-09-24) or "asrun" (the old /10)
island_offset = build_island_offset(
    ISLAND_OFFSET_FILE, HATTERAS_DOMAINS, mode=OFFSET_MODE)
OFFSET_TILTS = island_offset_tilts(island_offset, HATTERAS_DOMAINS)

_real = slice(HATTERAS_DOMAINS.start_real_index, HATTERAS_DOMAINS.end_real_index)
# Metres, except asrun (the file / 10)
_offset_m = island_offset * (DAM_TO_M if OFFSET_MODE == "asrun" else 1.0)
print(f"\n{START_YEAR} offsets ({OFFSET_MODE}): {island_offset.size} padded domains | "
      f"file span {_offset_m[_real].min():.0f}-{_offset_m[_real].max():.0f} m | "
      f"handed to Cascade {island_offset[_real].min():.0f}-"
      f"{island_offset[_real].max():.0f} m")


# 4. Period forcings -- RSLR, storms, SOURCE/SINK

# Everything here is resolved by the START_YEAR set in section 3.

# 4.1 Relative sea level rise

print(f"\nSEA_LEVEL_RISE_RATE = {SEA_LEVEL_RISE_RATE} m/yr")
print(f"  over {RUN_YEARS} years -> "
      f"{SEA_LEVEL_RISE_RATE * RUN_YEARS:.3f} m total rise")


# 4.2 Storm series

# One row per storm: time step, Rhigh, Rlow (dam), period, duration (h)


STORM_SERIES = load_storm_series(STORM_FILE)

reports.storm_report(storms=STORM_SERIES, storm_file=STORM_FILE,
                     run_years=RUN_YEARS)


# 4.3 Source/sink (background erosion)

# Per-domain rate, m/yr, as Barrier3D's Rat: (-) erosion, (+) accretion; absent domains 0.0


# Names the missing fit; an unsolved period's first edge probe starts empty when HAT_BE_OVERRIDE is set
try:
    DOMAIN_BE_RATES = be_rates(SOURCE_SINK_PRESET, START_YEAR)
except ValueError:
    if (SOURCE_SINK_PRESET == "edgeBE"
            and os.environ.get("HAT_BE_OVERRIDE", "").strip()):
        DOMAIN_BE_RATES = {}
    else:
        raise

# HAT_BE_OVERRIDE: a solve step without editing the config; the BE digest still tells it apart
_BE_OVERRIDE_RAW = os.environ.get("HAT_BE_OVERRIDE", "").strip()
if _BE_OVERRIDE_RAW:
    # Copied: the preset dict is the config's own
    DOMAIN_BE_RATES = dict(DOMAIN_BE_RATES)
    _BE_OVERRIDES = {}
    for _pair in _BE_OVERRIDE_RAW.split(","):
        if not _pair.strip():
            continue
        _gis, _, _rate = _pair.partition("=")
        if not _:
            raise ValueError(
                f"HAT_BE_OVERRIDE entry {_pair!r} is not 'gis=rate'. "
                f"Expected e.g. '1=-45.2,90=11.8'.")
        _BE_OVERRIDES[int(_gis)] = float(_rate)

    _unknown = sorted(g for g in _BE_OVERRIDES
                      if not HATTERAS_DOMAINS.first_gis_id <= g
                      <= HATTERAS_DOMAINS.last_gis_id)
    if _unknown:
        raise ValueError(
            f"HAT_BE_OVERRIDE names domain(s) {_unknown} outside the modelled "
            f"reach GIS {HATTERAS_DOMAINS.first_gis_id}-"
            f"{HATTERAS_DOMAINS.last_gis_id}.")

    print(f"\nHAT_BE_OVERRIDE       {len(_BE_OVERRIDES)} domain(s) forced off "
          f"preset {SOURCE_SINK_PRESET!r}")
    for _gis, _rate in sorted(_BE_OVERRIDES.items()):
        _was = DOMAIN_BE_RATES.get(_gis, 0.0)
        print(f"  GIS {_gis:<3}           {_was:+.4f} -> {_rate:+.4f} m/yr")
        DOMAIN_BE_RATES[_gis] = _rate

# An edgeBE end with no solved value and no override is refused
if SOURCE_SINK_PRESET == "edgeBE":
    _no_edge = [g for g in HATTERAS_BE_EDGE_DOMAINS
                if not DOMAIN_BE_RATES.get(g)]
    if _no_edge:
        raise SystemExit(
            f"\n[stop] edgeBE has no nonzero rate at end domain(s) "
            f"{_no_edge} of geometry {HATTERAS_GEOMETRY!r}. Supply one "
            f"with HAT_BE_OVERRIDE (e.g. '115=12.0'), or run zeroBE.\n")

BACKGROUND_EROSION_RATES = build_background_erosion(
    DOMAIN_BE_RATES, HATTERAS_DOMAINS)
USE_BACKGROUND_EROSION = any(rate != 0.0 for rate in BACKGROUND_EROSION_RATES)

# Derived from the preset, and checked against it
_EXPECT_BE_ON = SOURCE_SINK_PRESET != "zeroBE"
if USE_BACKGROUND_EROSION != _EXPECT_BE_ON:
    raise ValueError(
        f"preset {SOURCE_SINK_PRESET!r} implies "
        f"USE_BACKGROUND_EROSION={_EXPECT_BE_ON}, but the expanded rates give "
        f"{USE_BACKGROUND_EROSION}. The preset in hatteras_site_config.py does "
        f"not match its name.")

reports.background_erosion_report(
    preset=SOURCE_SINK_PRESET, start_year=START_YEAR,
    domain_rates=DOMAIN_BE_RATES, rates=BACKGROUND_EROSION_RATES,
    use_background_erosion=USE_BACKGROUND_EROSION, geometry=HATTERAS_DOMAINS,
    rates_2004_are_placeholder=HATTERAS_BE_RATES_2004_IS_PLACEHOLDER)


# 5. roadway_manager -- setbacks, per-domain elevation, historical events

# NC-12 forcing: setback (per period), elevation (m MHW), managed domains; events carry displacements

# Loaded and audited even with management off; it just never reaches a RoadwayManager

ROADWAY = roadway.RoadwayConfig()
_road_span = (HATTERAS_FIRST_ROAD_DOMAIN, HATTERAS_LAST_ROAD_DOMAIN)

# Setback: by period
road_setbacks_full, _missing_setbacks = roadway.load_road_setbacks(
    ROAD_SETBACK_FILE, HATTERAS_DOMAINS, *_road_span)

# Elevation: one set for every period
ROAD_ELEVATION_FILE = HATTERAS_DATA_BASE / HATTERAS_ROAD_ELEVATION_FILE
road_elevation_full, _missing_elevations = roadway.load_road_elevations(
    ROAD_ELEVATION_FILE, HATTERAS_DOMAINS, *_road_span, config=ROADWAY)

# Which domains CASCADE actually manages
ROADWAY_MANAGEMENT_ON = roadway.build_roadway_management_on(
    HATTERAS_DOMAINS, *_road_span,
    community_zones=HATTERAS_COMMUNITY_ZONES,
    enabled=ENABLE_ROADWAY_MANAGEMENT)

_road_slice = slice(HATTERAS_DOMAINS.gis_to_pad(HATTERAS_FIRST_ROAD_DOMAIN),
                    HATTERAS_DOMAINS.gis_to_pad(HATTERAS_LAST_ROAD_DOMAIN) + 1)
reports.roadway_report(
    setback_file=ROAD_SETBACK_FILE, setbacks=road_setbacks_full,
    missing_setbacks=_missing_setbacks,
    elevation_file=ROAD_ELEVATION_FILE, elevations=road_elevation_full,
    missing_elevations=_missing_elevations, road_slice=_road_slice,
    config=ROADWAY, roadway_on=ENABLE_ROADWAY_MANAGEMENT,
    management_on=ROADWAY_MANAGEMENT_ON,
    first_road_gis=HATTERAS_FIRST_ROAD_DOMAIN,
    last_road_gis=HATTERAS_LAST_ROAD_DOMAIN,
    community_zones=HATTERAS_COMMUNITY_ZONES,
    road_events=HATTERAS_ROAD_EVENTS,
    relocations_enabled=ENABLE_HISTORICAL_ROAD_RELOCATIONS)


# 5.1 Pre-flight audit: which road_offset will not survive year one

# A road whose flanking rows are >20% water drowns, and the domain is unmanaged from then on

road_audit = roadway.audit_setbacks(
    ELEVATION_FILE_PATHS, road_setbacks_full, HATTERAS_DOMAINS, *_road_span,
    management_on=ROADWAY_MANAGEMENT_ON, config=ROADWAY)
audit_summary = roadway.summarise_audit(road_audit)

reports.road_audit_report(audit=road_audit, summary=audit_summary)


# 6. beach_dune_manager -- nourishment schedule + overwash filter

# beach_dune_manager: always-on filter + fixed dune line, event-driven fills; overwash_filter is a PERCENT

# Schedule: one project list, period-filtered

# Every Hatteras project is in 2004-2024, so a 1984 run builds an empty schedule from this list
BN_SCHEDULE = nourishment.build_schedule(
    HATTERAS_NOURISHMENT_PROJECTS, HATTERAS_DOMAINS, START_YEAR, END_YEAR)

# The schedule the model is actually driven by

# Fills off = an empty schedule, not a skipped call; the footprint is unchanged
BN_SCHEDULE_APPLIED = BN_SCHEDULE if ENABLE_NOURISHMENT_FILLS else (
    nourishment.build_schedule([], HATTERAS_DOMAINS, START_YEAR, END_YEAR))

# What CASCADE is handed

# Zeros with the module off: the array states what the run does
OVERWASH_FILTER = (
    nourishment.build_overwash_filter(
        HATTERAS_DOMAINS, HATTERAS_COMMUNITY_ZONES, config=HATTERAS_BEACH_DUNE)
    if ENABLE_BEACH_DUNE_MANAGEMENT
    else [0.0] * HATTERAS_DOMAINS.total_domains)
OVERWASH_TO_DUNE = HATTERAS_BEACH_DUNE.overwash_to_dune_pct
BEACH_DUNE_MANAGEMENT_ON = nourishment.build_beach_dune_management_on(
    HATTERAS_DOMAINS, HATTERAS_COMMUNITY_ZONES, BN_SCHEDULE.nourished_gis,
    enabled=ENABLE_BEACH_DUNE_MANAGEMENT)

# Placeholder: rewritten every year by the schedule, so a missed schedule nourishes nothing
NOURISHMENT_VOLUME_INIT = [0.0] * HATTERAS_DOMAINS.total_domains

DOUBLE_MANAGED_GIS = nourishment.find_double_managed(
    BEACH_DUNE_MANAGEMENT_ON, ROADWAY_MANAGEMENT_ON, HATTERAS_DOMAINS)
BN_AUDIT = nourishment.audit_schedule(
    BN_SCHEDULE_APPLIED, BEACH_DUNE_MANAGEMENT_ON, config=HATTERAS_BEACH_DUNE)

reports.beach_dune_report(
    start_year=START_YEAR, end_year=END_YEAR,
    beach_dune_enabled=ENABLE_BEACH_DUNE_MANAGEMENT,
    fills_enabled=ENABLE_NOURISHMENT_FILLS,
    roadway_enabled=ENABLE_ROADWAY_MANAGEMENT,
    config=HATTERAS_BEACH_DUNE, management_on=BEACH_DUNE_MANAGEMENT_ON,
    overwash_to_dune=OVERWASH_TO_DUNE,
    community_zones=HATTERAS_COMMUNITY_ZONES,
    schedule=BN_SCHEDULE, schedule_applied=BN_SCHEDULE_APPLIED,
    audit=BN_AUDIT, double_managed=DOUBLE_MANAGED_GIS,
    geometry=HATTERAS_DOMAINS)


# 7. Hard structures / groin -- Buxton groin field

# The groin callback adds -M updrift, +M downdrift each year before the transport solve; the fillet is emergent

# 7.1 Switches, structure, sediment-budget reference

# GROIN_ENABLED is set in section 3; the guard stays here

# Only the sandbox Cascade has the pre-AST hook; elsewhere the groin silently does nothing
if GROIN_ENABLED and not USE_SANDBOX_CASCADE:
    raise RuntimeError(
        "GROIN_ENABLED=True requires USE_SANDBOX_CASCADE=True (section 1). "
        "cascade.cascade.Cascade.update() has no _groin_callback hook, so the "
        "groin would silently do nothing.\n"
        "USE_SANDBOX_CASCADE is pinned True and has no line in hat_run.yaml, "
        "so reaching this means HAT_USE_SANDBOX_CASCADE=0 is set in the "
        "environment. Unset it, or turn the groin off.")

# Net transport is southward, so updrift = north; the two domains share the boundary at GIS 5.5
GROIN_UPDRIFT_GIS = 6       # source: accretes
GROIN_DOWNDRIFT_GIS = 5     # sink:   erodes
GROIN_INSTALL_YEAR = 1969   # confirmed construction date

# Amplitude and deterioration floor: the two tunable knobs

# M and f come from the config (the sweep sets them per run); fit jointly across both periods
GROIN_TRAPPING_RATE_M_YR = RUN_CONFIG.groin_trapping_rate_m_yr
# b, the blocking fraction: read only when GROIN_KIND is "blocking"; f scales it as it scales M
GROIN_BLOCKING_FRACTION = RUN_CONFIG.groin_blocking_fraction
GROIN_M_PROVENANCE = ("joint two-period fit against the CoastSat D6-D5 "
                      "differential; see output/calibration/groin/ for the M-f "
                      "ridge and which grid bounds the solution touches")

# Deterioration: intact until the 2003 storm, failed from 2004

# Instant failure at the 2003 storm (from the 2004 step), as the observed D5-D6 gap shows
GROIN_DETERIORATION_DELAY_YEARS = 2004 - GROIN_INSTALL_YEAR   # = 35
GROIN_DETERIORATION_MODE = "instant"
GROIN_DETERIORATION_RAMP_YEARS = 0.0
GROIN_DETERIORATION_FRACTION = RUN_CONFIG.groin_deterioration_fraction

# Sediment-budget reference
REACH_TRANSPORT_LOSS_M3_YR = 5.9e5
REACH_TRANSPORT_CITATION = ("Inman & Dolan (1989), via Moore et al. (2010), "
                            "doi:10.1029/2009JF001299")
REACH_TRANSPORT_CAVEAT = (
    "reach-integrated transport-gradient LOSS, Oregon Inlet to Cape Hatteras "
    "(~60 km) -- a divergence, not a gross flux at Buxton, so this bounds "
    "order of magnitude only")

# Active profile height: both candidates printed; section 11 resolves it from the model
GROIN_PROFILE_HEIGHT_CANDIDATES_M = (12.0, 24.0)


# 7.2 Build the callback

# GROIN_CB is built either way for the 7.4 report; GROIN_CALLBACK is the one attached

_GROIN_COMMON = dict(
    updrift_pad=HATTERAS_DOMAINS.gis_to_pad(GROIN_UPDRIFT_GIS),
    downdrift_pad=HATTERAS_DOMAINS.gis_to_pad(GROIN_DOWNDRIFT_GIS),
    start_year=START_YEAR,
    install_year=GROIN_INSTALL_YEAR,
    n_domains=HATTERAS_DOMAINS.total_domains,
    deterioration_delay_years=GROIN_DETERIORATION_DELAY_YEARS,
    deterioration_mode=GROIN_DETERIORATION_MODE,
    deterioration_ramp_years=GROIN_DETERIORATION_RAMP_YEARS,
    deterioration_fraction=GROIN_DETERIORATION_FRACTION,
)
if GROIN_KIND == "blocking":
    GROIN_CB = BlockingGroinCallback(
        blocking_fraction=GROIN_BLOCKING_FRACTION, **_GROIN_COMMON)
else:
    GROIN_CB = GroinCallback(
        trapping_rate_m_yr=GROIN_TRAPPING_RATE_M_YR, **_GROIN_COMMON)
GROIN_CALLBACK = GROIN_CB if GROIN_ENABLED else None


# 7.4 Report

reports.groin_report(
    enabled=GROIN_ENABLED, callback=GROIN_CB,
    updrift_gis=GROIN_UPDRIFT_GIS, downdrift_gis=GROIN_DOWNDRIFT_GIS,
    install_year=GROIN_INSTALL_YEAR, start_year=START_YEAR, end_year=END_YEAR,
    geometry=HATTERAS_DOMAINS,
    trapping_rate_m_yr=GROIN_TRAPPING_RATE_M_YR,
    m_provenance=GROIN_M_PROVENANCE,
    groin_kind=GROIN_KIND, blocking_fraction=GROIN_BLOCKING_FRACTION,
    deterioration_mode=GROIN_DETERIORATION_MODE,
    deterioration_delay_years=GROIN_DETERIORATION_DELAY_YEARS,
    deterioration_ramp_years=GROIN_DETERIORATION_RAMP_YEARS,
    deterioration_fraction=GROIN_DETERIORATION_FRACTION,
    profile_height_candidates_m=GROIN_PROFILE_HEIGHT_CANDIDATES_M,
    reach_transport_loss_m3_yr=REACH_TRANSPORT_LOSS_M3_YR,
    reach_transport_citation=REACH_TRANSPORT_CITATION,
    reach_transport_caveat=REACH_TRANSPORT_CAVEAT,
    source_sink_preset=SOURCE_SINK_PRESET, domain_be_rates=DOMAIN_BE_RATES)


# 7.5 Scenario summary, and the run name derived from it

# The run-name suffix is derived from the switches, never typed

SCENARIO_SWITCHES = [
    ("period", f"{START_YEAR}-{END_YEAR} ({RUN_YEARS} yr)", None),
    ("source/sink preset", SOURCE_SINK_PRESET, SOURCE_SINK_PRESET),
    # No token for "asrun", so the older runs keep their names
    ("shoreline offset", OFFSET_MODE,
     None if OFFSET_MODE == "asrun" else f"offset{OFFSET_MODE}"),
    # No token of its own: implied by the preset, and checked against it in 4.3
    ("background erosion", USE_BACKGROUND_EROSION, None),
    ("roadway management", f"{sum(ROADWAY_MANAGEMENT_ON)} domains"
     if ENABLE_ROADWAY_MANAGEMENT else "off",
     "road" if ENABLE_ROADWAY_MANAGEMENT else "noroad"),
    ("historical relocations", ENABLE_HISTORICAL_ROAD_RELOCATIONS,
     "reloc" if ENABLE_HISTORICAL_ROAD_RELOCATIONS else None),
    ("beach/dune manager", f"{sum(BEACH_DUNE_MANAGEMENT_ON)} domains"
     if ENABLE_BEACH_DUNE_MANAGEMENT else "off",
     "bdm" if ENABLE_BEACH_DUNE_MANAGEMENT else "nobdm"),
    # "nonourish" only when there was fill to withhold and a module to withhold it in
    ("nourishment fills", f"{len(BN_SCHEDULE_APPLIED.projects)} applied"
     if ENABLE_NOURISHMENT_FILLS
     else ("suppressed" if BN_SCHEDULE.projects else "none in period"),
     "nourish" if BN_SCHEDULE_APPLIED.projects
     else ("nonourish" if BN_SCHEDULE.projects and ENABLE_BEACH_DUNE_MANAGEMENT
           else None)),
    ("groin", ("on" if GROIN_KIND == "dipole" else "on (blocking)")
     if GROIN_ENABLED else "off", GROIN_TOKEN),
    ("relocation target",
     "each domain's measured offset" if RELOCATION_SETBACK_M is None
     else f"{RELOCATION_SETBACK_M:g} m behind the dune line",
     RELOCATION_SETBACK_TOKEN),
    # Tokened only off the calibration value; last, so sweep_family reads the axis off it
    ("wave climate",
     f"Hs {RUN_CONFIG.hs} m, Tp {RUN_CONFIG.wave_period_s} s, "
     f"asym {RUN_CONFIG.wave_asymmetry}, "
     f"high-angle {RUN_CONFIG.wave_angle_high_fraction}"
     + ("" if WAVE_TOKEN is None else "  (off calibration)"),
     WAVE_TOKEN),
]

RUN_NAME_SUFFIX = "_".join(
    token for _, _, token in SCENARIO_SWITCHES if token)
RUN_NAME_BASE = f"{RUN_NAME_STEM}_{RUN_NAME_SUFFIX}"

# A switch that did not reach its module is raised, not warned
if RUN_NAME_BASE != RUN_NAME_PREVIEW:
    raise AssertionError(
        f"run name disagrees with the section 3 preview\n"
        f"  section 3    {RUN_NAME_PREVIEW}\n"
        f"  section 7.5  {RUN_NAME_BASE}\n"
        f"The preview follows the switches; this follows the built modules.")

reports.scenario_summary_report(
    scenario=SCENARIO, departures=_SCENARIO_DEPARTURES,
    switches=SCENARIO_SWITCHES, run_name_base=RUN_NAME_BASE,
    double_managed=DOUBLE_MANAGED_GIS, groin_enabled=GROIN_ENABLED,
    updrift_gis=GROIN_UPDRIFT_GIS,
    beach_dune_on_updrift=BEACH_DUNE_MANAGEMENT_ON[
        HATTERAS_DOMAINS.gis_to_pad(GROIN_UPDRIFT_GIS)])


# 8. CoastSat target rates -- LOWESS windows

# The CoastSat target: LOWESS at transect resolution, averaged to domains; GIS 1-10 unsmoothed (display only)

# 8.1 Load both periods

# Both periods always load: section 12 draws the other as reference

COASTSAT_DATASETS = [
    CoastSatDataset(
        label="CoastSat LRR (1984-2004)",
        period_start=1984,
        csv_path=str(COASTSAT_BASE_DIR / "1984_2004" / "transect_lrr_full.csv"),
    ),
    CoastSatDataset(
        label="CoastSat LRR (2004-2024)",
        period_start=2004,
        csv_path=str(COASTSAT_BASE_DIR / "2004_2024" / "transect_lrr_full.csv"),
    ),
    # Every window loads on every run
    CoastSatDataset(
        label="CoastSat LRR (1996-2010)",
        period_start=1996,
        csv_path=str(COASTSAT_BASE_DIR / "1996_2010" / "transect_lrr_full.csv"),
    ),
    CoastSatDataset(
        label="CoastSat LRR (2010-2024)",
        period_start=2010,
        csv_path=str(COASTSAT_BASE_DIR / "2010_2024" / "transect_lrr_full.csv"),
    ),
]

# 7-domain smoothing, the group's range (2026-09-28)
LOWESS_CONFIG = LowessConfig(window_domains=(7,), skip_southern_domains=10)

# The reference window, named: rate_comparison uses max(window_domains)
TARGET_WINDOW = 7
if TARGET_WINDOW != max(LOWESS_CONFIG.window_domains):
    raise ValueError(
        f"TARGET_WINDOW={TARGET_WINDOW} but rate_comparison will use "
        f"max(window_domains)={max(LOWESS_CONFIG.window_domains)} as the "
        f"reference curve -- section 12 would compare against a different "
        f"curve than the one reported here.")

# An extended geometry appends its own CoastSat rows beyond GIS 90
if HATTERAS_GEOMETRY_EXTENDED:
    COASTSAT_DATASETS = [
        CoastSatDataset(label=ds.label + " + extension",
                        period_start=ds.period_start,
                        csv_path=str(Path(ds.csv_path).parent / "ext"
                                     / "transect_lrr_with_base.csv"))
        if ds.period_start == START_YEAR else ds
        for ds in COASTSAT_DATASETS]

cs_series = build_coastsat_series(
    COASTSAT_DATASETS, active_period_start=START_YEAR,
    lowess_config=LOWESS_CONFIG, domains=HATTERAS_DOMAINS)

CS_ACTIVE = next((cs for cs in cs_series if cs["active"]), None)
if CS_ACTIVE is None:
    raise RuntimeError(
        f"no CoastSat dataset has period_start == START_YEAR ({START_YEAR}); "
        f"loaded: {[cs['period_start'] for cs in cs_series]}")


# 8.2 The target, as a table


COASTSAT_TARGET = build_target_table(
    CS_ACTIVE, LOWESS_CONFIG, HATTERAS_DOMAINS, TARGET_WINDOW)

# The surveyed reach's own target: what the interior score is graded against in every geometry
COASTSAT_TARGET_BASE = COASTSAT_TARGET
if HATTERAS_GEOMETRY_EXTENDED:
    _cs_base = build_coastsat_series(
        [CoastSatDataset(
            label=f"CoastSat LRR ({START_YEAR}-{END_YEAR}) surveyed",
            period_start=START_YEAR,
            csv_path=str(COASTSAT_BASE_DIR / f"{START_YEAR}_{END_YEAR}"
                         / "transect_lrr_full.csv"))],
        active_period_start=START_YEAR, lowess_config=LOWESS_CONFIG,
        domains=DEFAULT_DOMAINS)
    COASTSAT_TARGET_BASE = build_target_table(
        _cs_base[0], LOWESS_CONFIG, DEFAULT_DOMAINS, TARGET_WINDOW)


# 8.4 Report

reports.coastsat_report(
    target=COASTSAT_TARGET, active=CS_ACTIVE, target_window=TARGET_WINDOW,
    lowess_config=LOWESS_CONFIG, geometry=HATTERAS_DOMAINS, cs_series=cs_series,
    updrift_gis=GROIN_UPDRIFT_GIS, downdrift_gis=GROIN_DOWNDRIFT_GIS)


# 9. Figure configuration

# Keyword bundles for every figure call: omit the annotations and a figure loses its geography silently

# 9.1 Site config every figure needs

RATE_FIG_KWARGS = dict(
    domains=HATTERAS_DOMAINS,
    annotations=HATTERAS_ANNOTATIONS,
    lowess_config=LOWESS_CONFIG,        # section 8's, not DEFAULT_LOWESS
    config=DEFAULT_RATE_COMPARISON,
)

gif_config = GifConfig(
    fps=3,
    year_stride=1,
    annotate=True,
    auto_open=False,
    keep_frames=False,
    save_matrix=True,          # lets a later run difference against this one
    ocean_at_bottom=True,      # Hatteras' real cross-shore layout
    baseline_label="no-groin baseline",
    target_label=f"observed {END_YEAR} dune line",   # section 9.4's target
)

GIF_KWARGS = dict(
    domains=HATTERAS_DOMAINS,
    annotations=HATTERAS_ANNOTATIONS,
    gif_config=gif_config,
)

# Figure-level conventions the section 12 calls read.
FLIP_SIGN_MODEL = True          # x_s_TS increases landward; flip so up = seaward
PLOT_REAL_DOMAINS_ONLY = True   # GIS 1-90 axis; False adds the buffer domains

# The estimator the figures draw: "lrr" matches the CoastSat target; "endpoint" only redraws old figures
RATE_ESTIMATOR = "lrr"
if RATE_ESTIMATOR not in ("lrr", "endpoint"):
    raise ValueError(f"RATE_ESTIMATOR must be 'lrr' or 'endpoint', "
                     f"got {RATE_ESTIMATOR!r}")

# Reporting threshold only: nothing is dropped or flagged below it
LRR_R2_FLOOR = 0.50


# 9.2 Animation jobs

# range: "real" | "all" | "groin" | "groin_span" | (lo, hi); mode: position | displacement | difference

# make_gifs false empties the job list; the call still writes the shoreline matrix
GIF_JOBS = [
    dict(range="real", mode="displacement"),
    dict(range="real", mode="position"),
    dict(range="groin", mode="position", pad=9),
    dict(range="groin", mode="difference", pad=9),
] if RUN_CONFIG.make_gifs else []


# 9.3 Baseline for difference jobs


GIF_BASELINE_NAME = None
GIF_BASELINE_NPY = None
if GROIN_ENABLED:
    GIF_BASELINE_NAME = scenario_run_name(
        SCENARIO_SWITCHES, RUN_NAME_STEM, groin="nogroin")
    _baseline = resolve(OUTPUT_BASE_DIR / GIF_BASELINE_NAME, "matrix",
                        GIF_BASELINE_NAME)
    GIF_BASELINE_NPY = str(_baseline) if _baseline.exists() else None


# 9.4 Validation target -- the surveyed island position in the end year

# From the raw transect CSVs (fixed datum), not the padded offsets
from site_layer.hat_topo_version import RAW_OFFSET_DIR  # noqa: E402  (2026-09-18)


# 9.5 Report, and the annotation guard

# An empty AnnotationConfig would strip every figure's geography: fail here

_ann_populated = any([HATTERAS_ANNOTATIONS.town_spans,
                      HATTERAS_ANNOTATIONS.village_lines,
                      HATTERAS_ANNOTATIONS.piers,
                      HATTERAS_ANNOTATIONS.groins,
                      HATTERAS_ANNOTATIONS.shoal_zones])
if not _ann_populated:
    raise RuntimeError(
        "HATTERAS_ANNOTATIONS is empty -- every figure would render with no "
        "geographic layer and no error. Check the import in section 1.")

reports.figure_config_report(
    annotations=HATTERAS_ANNOTATIONS, lowess_config=LOWESS_CONFIG,
    gif_config=gif_config, gif_jobs=GIF_JOBS,
    flip_sign_model=FLIP_SIGN_MODEL,
    real_domains_only=PLOT_REAL_DOMAINS_ONLY,
    groin_enabled=GROIN_ENABLED, baseline_npy=GIF_BASELINE_NPY,
    baseline_name=GIF_BASELINE_NAME, output_base_dir=OUTPUT_BASE_DIR)


# 10. build_cascade + run_cascade_simulation

# Both live in cascade_pipeline/hindcast.py; run_years counts transitions, not states


print("\nbuild_cascade + run_cascade_simulation defined")
print(f"  run_years -> {RUN_YEARS} transitions, "
      f"time_step_count={RUN_YEARS + 1}, {RUN_YEARS + 1} annual states "
      f"({START_YEAR}-{END_YEAR})")
print("  nourishment via BN_SCHEDULE.apply_to_cascade -> "
      "cascade.nourishment_volume")
print("  road events via cascade_pipeline.roadway.apply_historical_event")


# 11. Initialize CASCADE -- single config, no sweep

# 11.1 Parameters sections 2-9 do not produce

NUM_CORES = 1        # >1 has crashed on this configuration; leave at 1

# Dune growth (Barrier3D logistic growth bounds)
RMIN = [0.55] * HATTERAS_DOMAINS.total_domains
RMAX = [0.95] * HATTERAS_DOMAINS.total_domains

# Dune rebuild thresholds, m MHW

# In metres, the documented unit
DUNE_DESIGN_ELEVATION_M = 3.0    # rebuild target
DUNE_MINIMUM_ELEVATION_M = 0.0   # rebuild trigger: CASCADE's berm floor governs, explicitly
DUNE_DESIGN_ELEVATION = [DUNE_DESIGN_ELEVATION_M] * HATTERAS_DOMAINS.total_domains
DUNE_MINIMUM_ELEVATION = [DUNE_MINIMUM_ELEVATION_M] * HATTERAS_DOMAINS.total_domains

# Roadway
ROAD_ELEVATION = road_elevation_full   # per-domain, m MHW, from section 5
ROAD_WIDTH = 20.0

# Wave climate: one configuration, no sweep

# From hat_run.yaml's physics block, so a sweep never edits this file
Hs = RUN_CONFIG.hs                                       # m, calibration 2.5
FIXED_WAVE_PERIOD = RUN_CONFIG.wave_period_s             # s, calibration 8
FIXED_WAVE_ASYMMETRY = RUN_CONFIG.wave_asymmetry         # calibration 0.7
FIXED_WAVE_ANGLE_HIGH_FRACTION = RUN_CONFIG.wave_angle_high_fraction  # 0.1

# The run was named for these in 7.5: they must still match
assert _WAVE_VALUES == {
    "hs": Hs, "wave_period_s": FIXED_WAVE_PERIOD,
    "wave_asymmetry": FIXED_WAVE_ASYMMETRY,
    "wave_angle_high_fraction": FIXED_WAVE_ANGLE_HIGH_FRACTION,
}, "section 11's wave climate is not the one section 3 named the run for"

# Datums
BERM_ELEVATION = 1.7    # m NAVD88, Hatteras Island, NCDOT-derived via NC State
MHW_ELEVATION = 0.36    # m NAVD88, Duck NC gauge (NOAA 8651370)

# Sandbags: off for the hindcast

# From hat_run.yaml (top-level sandbags); reaches run_index.csv as sandbags_on
ENABLE_SANDBAG_PLACEMENT = RUN_CONFIG.sandbags
SANDBAG_MANAGEMENT_ON = [ENABLE_SANDBAG_PLACEMENT] * HATTERAS_DOMAINS.total_domains
SANDBAG_ELEVATION = 0

SEA_LEVEL_CONSTANT = True

# Run identity, from section 7.5's derived name

# OVERWRITE: False stops on an existing result; True empties the directory and reuses it
OVERWRITE = RUN_CONFIG.overwrite

RUN_NAME = RUN_NAME_BASE
RUN_DIR = str(OUTPUT_BASE_DIR / RUN_NAME)
# Read before the guard, which empties the directory when OVERWRITE is set
_replacing = run_dir_contents(RUN_DIR)
guard_run_dir(RUN_DIR, overwrite=OVERWRITE)
print(f"\nRUN_DIR               {RUN_DIR}"
      + ("   (OVERWRITE=True)" if OVERWRITE else ""))
if _replacing:
    print(f"  replaced            {len(_replacing)} file(s) from the previous run")


# 11.2 Build

# Each run gets its own parameter file (CASCADE rewrites it while it builds); the template stays read-only
RUN_PARAMETER_FILE = Path(RUN_DIR) / f"{RUN_NAME}-parameters.yaml"
shutil.copyfile(HATTERAS_DATA_BASE / PARAMETER_FILE, RUN_PARAMETER_FILE)

cascade = build_cascade(
    run_years=RUN_YEARS,
    name=RUN_NAME,
    storm_file=str(STORM_FILE),
    alongshore_section_count=HATTERAS_DOMAINS.total_domains,
    num_cores=NUM_CORES,
    rmin=RMIN, rmax=RMAX,
    elevation_file=ELEVATION_FILE_PATHS,
    dune_file=DUNE_FILE_PATHS,
    dune_design_elevation=DUNE_DESIGN_ELEVATION,
    dune_minimum_elevation=DUNE_MINIMUM_ELEVATION,
    road_ele=ROAD_ELEVATION,
    road_width=ROAD_WIDTH,
    road_setback=road_setbacks_full,
    overwash_filter=OVERWASH_FILTER,
    overwash_to_dune=OVERWASH_TO_DUNE,
    nourishment_volume=NOURISHMENT_VOLUME_INIT,
    background_erosion=BACKGROUND_EROSION_RATES,
    roadway_management_on=ROADWAY_MANAGEMENT_ON,
    beach_dune_manager_on=BEACH_DUNE_MANAGEMENT_ON,
    sea_level_rise_rate=SEA_LEVEL_RISE_RATE,
    sea_level_constant=SEA_LEVEL_CONSTANT,
    sandbag_management_on=SANDBAG_MANAGEMENT_ON,
    sandbag_elevation=SANDBAG_ELEVATION,
    enable_shoreline_offset=True,
    shoreline_offset=island_offset,
    wave_height=Hs,
    wave_period=FIXED_WAVE_PERIOD,
    wave_asymmetry=FIXED_WAVE_ASYMMETRY,
    wave_angle_high_fraction=FIXED_WAVE_ANGLE_HIGH_FRACTION,
    berm_elevation=BERM_ELEVATION,
    MHW=MHW_ELEVATION,
    data_base=HATTERAS_DATA_BASE,
    parameter_file=str(RUN_PARAMETER_FILE),
    groin_callback=GROIN_CALLBACK,
    relocation_setback_m=RELOCATION_SETBACK_M,
)


# 11.3 Pre-run diagnostics, and the groin prediction

# Read off the built model; the fillet prediction is written down before the run


# r_ipl at the groin cell's starting angle (shore-normal is negative under option A)
if GROIN_ENABLED:
    _up = HATTERAS_DOMAINS.gis_to_pad(GROIN_UPDRIFT_GIS)
    R_IPL_THETA_DEG = float(np.degrees(np.arctan2(
        island_offset[_up + 1] - island_offset[_up],
        HATTERAS_DOMAINS.domain_spacing_m)))
else:
    R_IPL_THETA_DEG = 0.0
R_IPL = brie_r_ipl(cascade, theta_deg=R_IPL_THETA_DEG)
_brie = cascade._brie_coupler._brie
_d_sf_m = float(_brie.d_sf)
_h_b_m = float(cascade.barrier3d[0].h_b_TS[0]) * DAM_TO_M
_profile_height_m = _d_sf_m + _h_b_m
_berm_floor_m = float(cascade.barrier3d[0].BermEl) * DAM_TO_M

# The M / (4 r_ipl) prediction is the dipole's; a blocking groin has no a-priori amplitude
if GROIN_CALLBACK is not None and GROIN_KIND == "dipole":
    (GROIN_PREDICTED_AMPLITUDE_M,
     GROIN_PREDICTED_EXTENT_DOMAINS,
     GROIN_PREDICTED_EXTENT_M) = predict_fillet(
        trapping_rate_m_yr=GROIN_CALLBACK.M,
        r_ipl=R_IPL,
        run_years=RUN_YEARS,
        dy_m=HATTERAS_DOMAINS.domain_spacing_m,
    )
else:
    GROIN_PREDICTED_AMPLITUDE_M = None
    GROIN_PREDICTED_EXTENT_DOMAINS = None
    GROIN_PREDICTED_EXTENT_M = None

reports.pre_run_report(
    run_name=RUN_NAME, run_dir=RUN_DIR, run_years=RUN_YEARS,
    start_year=START_YEAR, end_year=END_YEAR, geometry=HATTERAS_DOMAINS,
    roadway_on=ROADWAY_MANAGEMENT_ON, beach_dune_on=BEACH_DUNE_MANAGEMENT_ON,
    wave_height=Hs, d_sf_m=_d_sf_m, h_b_m=_h_b_m,
    profile_height_m=_profile_height_m, berm_floor_m=_berm_floor_m,
    dune_design_elevation_m=DUNE_DESIGN_ELEVATION_M,
    dune_minimum_elevation_m=DUNE_MINIMUM_ELEVATION_M,
    groin_callback=GROIN_CALLBACK, groin_enabled=GROIN_ENABLED, r_ipl=R_IPL,
    predicted_amplitude_m=GROIN_PREDICTED_AMPLITUDE_M,
    predicted_extent_domains=GROIN_PREDICTED_EXTENT_DOMAINS,
    predicted_extent_m=GROIN_PREDICTED_EXTENT_M,
    reach_transport_loss_m3_yr=REACH_TRANSPORT_LOSS_M3_YR)


# 12. Run the loop, verify, then figures

# Verification before figures, and the nourishment check asserts

GROIN_EXTENT_THRESHOLD_FRAC = 0.10   # fraction of peak effect defining "extent"

# 12.1 Run

# A model already stepped would run past TMAX

if len(cascade.barrier3d[0].x_s_TS) > 1:
    raise RuntimeError(
        f"this cascade has already been stepped "
        f"({len(cascade.barrier3d[0].x_s_TS)} states). Rebuild the model "
        f"(section 11) before running section 12 again.")

_t0 = time.perf_counter()
_run_kwargs = dict(
    cascade=cascade,
    run_years=RUN_YEARS,
    name=RUN_NAME,
    run_dir=RUN_DIR,
    start_year=START_YEAR,
    geometry=HATTERAS_DOMAINS,
    alongshore_section_count=HATTERAS_DOMAINS.total_domains,
    historical_road_events=HATTERAS_ROAD_EVENTS,
    relocations_enabled=ENABLE_HISTORICAL_ROAD_RELOCATIONS,
    setback_check=HATTERAS_RELOCATION_CHECK_2004,
    nourishment_schedule=BN_SCHEDULE_APPLIED,
    groin_callback=GROIN_CALLBACK,
    save_model_state=SAVE_MODEL_STATE,
)
if tqdm is not None:
    with tqdm(total=RUN_YEARS, desc=f"{RUN_NAME}", unit="yr") as _bar:
        cascade = run_cascade_simulation(progress=_bar, **_run_kwargs)
else:
    cascade = run_cascade_simulation(progress=None, **_run_kwargs)
RUN_SECONDS = time.perf_counter() - _t0
print(f"\nruntime               {RUN_SECONDS / 60:.1f} min "
      f"({RUN_SECONDS / RUN_YEARS:.1f} s per model year)")


# 12.2 Verify

shoreline_m = build_shoreline_matrix(cascade)
_states, _ = shoreline_m.shape

reports.run_length_report(states=_states, run_years=RUN_YEARS)

# The denominator behind the old 5% bias: checked, not trusted
assert _states - 1 == RUN_YEARS, (
    f"span mismatch: {_states} states implies {_states - 1} years elapsed, "
    f"but RUN_YEARS is {RUN_YEARS}")
# Two estimators: change_rate (endpoint, conserves net movement) and model_lrr (OLS, matches the target)
change_rate = compute_change_rate(
    shoreline_m, span_years=RUN_YEARS, flip_sign=FLIP_SIGN_MODEL)
model_lrr, model_lrr_r2 = compute_lrr(
    shoreline_m, span_years=RUN_YEARS, flip_sign=FLIP_SIGN_MODEL)

# The estimator section 12.5 draws, resolved once
PLOTTED_RATE = model_lrr if RATE_ESTIMATOR == "lrr" else change_rate

# Nourishment: the model's own record, not the schedule's intent
BN_REPORT = nourishment.verify_nourishment(
    cascade, BN_SCHEDULE_APPLIED, BEACH_DUNE_MANAGEMENT_ON)
reports.nourishment_report(report=BN_REPORT, run_years=RUN_YEARS)
assert BN_REPORT["ok"], "nourishment did not reach the model as scheduled"

# The double-management consequence section 6 predicted
reports.frozen_setbacks_report(
    double_managed=DOUBLE_MANAGED_GIS,
    rows=(nourishment.verify_setbacks_frozen(
              cascade, DOUBLE_MANAGED_GIS, HATTERAS_DOMAINS)
          if DOUBLE_MANAGED_GIS else ()))

# Which road_offset survived
ROAD_SUMMARY = roadway.summarise_road_management(
    cascade, HATTERAS_DOMAINS, HATTERAS_FIRST_ROAD_DOMAIN,
    HATTERAS_LAST_ROAD_DOMAIN)
_drowned, _blocked = reports.roadway_outcome_report(summary=ROAD_SUMMARY)


# 12.3 Groin: the pre-registered extent check


GROIN_EXTENT = None
if GROIN_CALLBACK is None:
    print(f"\nGROIN EXTENT          skipped: groin not enabled in this run")
elif not GIF_BASELINE_NPY:
    print(f"\nGROIN EXTENT          skipped: no paired no-groin baseline")
    print(f"  expected            {GIF_BASELINE_NAME}")
    print(f"  Run once with GROIN_ENABLED = False to create it, then re-run.")
else:
    _baseline_m = np.load(GIF_BASELINE_NPY)
    GROIN_EXTENT = measure_groin_extent(
        shoreline_m, _baseline_m, HATTERAS_DOMAINS,
        GROIN_UPDRIFT_GIS, GROIN_DOWNDRIFT_GIS, GROIN_EXTENT_THRESHOLD_FRAC)
    reports.groin_extent_report(
        extent=GROIN_EXTENT, threshold_frac=GROIN_EXTENT_THRESHOLD_FRAC,
        baseline_name=GIF_BASELINE_NAME,
        predicted_extent_domains=GROIN_PREDICTED_EXTENT_DOMAINS,
        predicted_extent_m=GROIN_PREDICTED_EXTENT_M)


# 12.4 Write

_shoreface_depth_m = float(cascade._brie_coupler._brie.d_sf)

run = RunInfo(
    run_name=RUN_NAME, run_dir=RUN_DIR,
    start_year=START_YEAR, end_year=END_YEAR, Hs=Hs,
    flip_sign_model=FLIP_SIGN_MODEL,
    background_erosion_on=USE_BACKGROUND_EROSION,
)

_rate_csv = str(write_path(RUN_DIR, "rate_csv", RUN_NAME))
# Both estimators ship, with lrr_r2 beside them
pd.DataFrame({
    "gis_domain": np.arange(HATTERAS_DOMAINS.first_gis_id,
                            HATTERAS_DOMAINS.last_gis_id + 1),
    "change_rate_m_yr": change_rate[_real],
    "lrr_m_yr": model_lrr[_real],
    "lrr_r2": model_lrr_r2[_real],
}).to_csv(_rate_csv, index=False)
print(f"\nwrote                 {os.path.basename(_rate_csv)}")

if ROAD_SUMMARY:
    _road_csv = str(write_path(RUN_DIR, "road_csv", RUN_NAME))
    pd.DataFrame(ROAD_SUMMARY).to_csv(_road_csv, index=False)
    print(f"                      {os.path.basename(_road_csv)}")

# Skill against the section 8 target

# Island-wide and interior spans both reported; SKILL is the LRR one, SKILL_ENDPOINT kept beside it
SKILL = skill_vs_target(model_lrr, COASTSAT_TARGET, HATTERAS_DOMAINS)
SKILL_ENDPOINT = skill_vs_target(change_rate, COASTSAT_TARGET,
                                 HATTERAS_DOMAINS)
# Interior = GIS 2-89 against the surveyed target in every geometry
for _skill, _rate in ((SKILL, model_lrr), (SKILL_ENDPOINT, change_rate)):
    _skill["mean_bias_reach_interior_m_yr"] = _skill["mean_bias_interior_m_yr"]
    _skill["rmse_reach_interior_m_yr"] = _skill["rmse_interior_m_yr"]
    _base = skill_vs_target(_rate, COASTSAT_TARGET_BASE, HATTERAS_DOMAINS,
                            interior_gis=SCORE_INTERIOR_GIS)
    for _key in ("mean_bias_interior_m_yr", "rmse_interior_m_yr",
                 "n_domains_interior"):
        _skill[_key] = _base[_key]
print(f"\nSKILL vs CoastSat     model LRR - target LRR, m/yr")
print(f"  island-wide         bias {SKILL['mean_bias_m_yr']:+.3f}   "
      f"RMSE {SKILL['rmse_m_yr']:.3f}   (n={SKILL['n_domains']})")
print(f"  interior (GIS 2-89) bias {SKILL['mean_bias_interior_m_yr']:+.3f}   "
      f"RMSE {SKILL['rmse_interior_m_yr']:.3f}   "
      f"(n={SKILL['n_domains_interior']})")
print(f"  endpoint estimator  bias "
      f"{SKILL_ENDPOINT['mean_bias_interior_m_yr']:+.3f}   "
      f"RMSE {SKILL_ENDPOINT['rmse_interior_m_yr']:.3f}   "
      f"(interior; the pre-LRR metric)")

# Where a straight line summarises the trajectory poorly (nourished or barely moving domains)
_lrr_r2_real = model_lrr_r2[_real]
_poor_fit = [int(_gis) for _gis, _r2 in zip(
    range(HATTERAS_DOMAINS.first_gis_id, HATTERAS_DOMAINS.last_gis_id + 1),
    _lrr_r2_real) if np.isfinite(_r2) and _r2 < LRR_R2_FLOOR]
print(f"\nLRR FIT QUALITY       median r2 "
      f"{np.nanmedian(_lrr_r2_real):.3f}   "
      f"{len(_poor_fit)} of {_lrr_r2_real.size} domains below "
      f"{LRR_R2_FLOOR:.2f}")
if _poor_fit:
    print(f"  a step, or no trend GIS "
          f"{', '.join(str(_g) for _g in _poor_fit[:12])}"
          f"{' ...' if len(_poor_fit) > 12 else ''}")

# Run metadata: the scenario, plus what distinguishes this run

# One structure renders both the .txt and the .json
_GIT = git_provenance(PROJECT_BASE_DIR)
_TIMESTAMP = timestamp()
GENERATED_BY = "HAT_hindcast_1984_2024.py"

# What the source/sink preset actually was

# The preset's name and its numbers are recorded separately
_BE_NONZERO = sum(1 for _rate in DOMAIN_BE_RATES.values() if _rate)
_BE_DIGEST = values_digest(DOMAIN_BE_RATES)
_BE_EDGE_RATES = {f"rate_gis{_gis}_m_yr": DOMAIN_BE_RATES.get(_gis, 0.0)
                  for _gis in HATTERAS_BE_EDGE_DOMAINS}

_META = {
    "identity": {
        "run_name": RUN_NAME,
        "timestamp": _TIMESTAMP,
        "generated_by": GENERATED_BY,
        "runtime_min": f"{RUN_SECONDS / 60:.1f}",
        # The three that drift across a batch: Cascade class, extractor version, commit
        "use_sandbox_cascade": (USE_SANDBOX_CASCADE,
                                "cascade.cascade_groin, not cascade.cascade"),
        # Product and version both: version numbers restart per product
        "topo_product": TOPO_PRODUCT,
        "topo_dune_version": TOPO_DUNE_VERSION,
        # The island offset version, which the run name does not carry
        "island_offset_version": ISLAND_OFFSET_VERSION,
        # The reach: "base" is GIS 1-90; an extended one is not a matrix run
        "geometry": (HATTERAS_GEOMETRY,
                     f"GIS {HATTERAS_DOMAINS.first_gis_id} to "
                     f"{HATTERAS_DOMAINS.last_gis_id}"),
        "run_kind": RUN_KIND,
        "run_tag": RUN_TAG,
        "parameter_file": (RUN_PARAMETER_FILE.name,
                           f"this run's copy of {PARAMETER_FILE}; CASCADE "
                           f"rewrites the copy, not the template"),
        "save_model_state": SAVE_MODEL_STATE,
        "git_commit": _GIT["commit"],
        "git_branch": _GIT["branch"],
        "git_dirty": (_GIT["dirty"],
                      "True: the commit alone does not reproduce this run"),
        # Barrier3D's branch is part of the model
        "barrier3d_branch": _B3D["branch"],
        "barrier3d_commit": _B3D["commit"],
        "barrier3d_dirty": _B3D["dirty"],
        "barrier3d_route_overwash_fix": (_B3D["route_overwash_fix"],
                                         "True: the loaded Barrier3D has the "
                                         "2026-09-24 route_overwash index fix"),
        "barrier3d_gap_momentum_fix": (_B3D.get("gap_momentum_fix"),
                                       "True: DuneGaps, gap discharge slice and "
                                       "inundation momentum fixed (2026-09-28)"),
        "barrier3d_per_cell_ceiling": _B3D.get("per_cell_ceiling"),
        # Read off the constructed model, not the template
        "dune_ceiling": ("per-cell, from the starting dunes"
                         if getattr(cascade.barrier3d[0], "_DuneCeilingFromStart", False)
                         else f"uniform Dmaxel {cascade.barrier3d[0].Dmaxel * 10 + MHW_ELEVATION:.2f} m NAVD88"),
        "storm_file": str(STORM_FILE.relative_to(HATTERAS_DATA_BASE)),
        # beach_dune_manager's 4 m bulldozer cap: what it clips (2026-09-28).
        "bdm_dune_cap_applies_to": getattr(_bdm_module, "DUNE_CAP_APPLIES_TO", "whole dune cell"),
    },
    "scenario": {label: value for label, value, _token in SCENARIO_SWITCHES},
    "period": {
        "start_year": START_YEAR,
        "end_year": END_YEAR,
        "run_years": (RUN_YEARS, "transitions"),
        "annual_states": _states,
    },
    "wave climate": {
        "wave_height_m": Hs,
        "wave_period_s": FIXED_WAVE_PERIOD,
        "wave_asymmetry": FIXED_WAVE_ASYMMETRY,
        "wave_angle_high_frac": FIXED_WAVE_ANGLE_HIGH_FRACTION,
        "shoreface_depth_m": (f"{_shoreface_depth_m:.2f}", "8.9 * Hs"),
    },
    "sea level": {
        "rslr_m_yr": SEA_LEVEL_RISE_RATE,
        "rslr_constant": SEA_LEVEL_CONSTANT,
    },
    "dunes and roadway": {
        "rmin / rmax": f"{RMIN[0]} / {RMAX[0]}",
        "dune_design_ele_m_mhw": (DUNE_DESIGN_ELEVATION_M,
                                  "floored by roadway_manager"),
        "dune_min_ele_m_mhw": (DUNE_MINIMUM_ELEVATION_M,
                               "floored by roadway_manager"),
        "road_width_m": ROAD_WIDTH,
        "relocations_enabled": ENABLE_HISTORICAL_ROAD_RELOCATIONS,
    },
    "source/sink": {
        "preset": SOURCE_SINK_PRESET,
        "background_erosion_on": (USE_BACKGROUND_EROSION,
                                  "implied by the preset, checked in 4.3"),
        "domains_specified": len(DOMAIN_BE_RATES),
        "nonzero_domains": _BE_NONZERO,
        "values_digest": (_BE_DIGEST,
                          "changes if any rate changed, interior included"),
        **_BE_EDGE_RATES,
    },
    "groin": {"enabled": GROIN_ENABLED},
    "skill": {
        "target": (f"CoastSat LOWESS {TARGET_WINDOW}-domain",
                   "raw means over GIS 1-"
                   f"{LOWESS_CONFIG.skip_southern_domains}"),
        "estimator": (RATE_ESTIMATOR,
                      "OLS slope through every annual state, matching the "
                      "target's own definition"),
        "mean_bias_m_yr": f"{SKILL['mean_bias_m_yr']:+.4f}",
        "rmse_m_yr": f"{SKILL['rmse_m_yr']:.4f}",
        "mean_bias_interior_m_yr": (
            f"{SKILL['mean_bias_interior_m_yr']:+.4f}",
            "GIS 2-89: the locked end domains excluded"),
        "rmse_interior_m_yr": f"{SKILL['rmse_interior_m_yr']:.4f}",
        "lrr_r2_median": (f"{np.nanmedian(_lrr_r2_real):.4f}",
                          "how well a line describes the modelled trajectory"),
        "lrr_r2_below_floor": (f"{len(_poor_fit)}",
                               f"domains under r2 {LRR_R2_FLOOR:.2f}"),
        "endpoint_mean_bias_interior_m_yr": (
            f"{SKILL_ENDPOINT['mean_bias_interior_m_yr']:+.4f}",
            "the pre-LRR estimator; calibBE and groin M were fit on this"),
        "endpoint_rmse_interior_m_yr":
            f"{SKILL_ENDPOINT['rmse_interior_m_yr']:.4f}",
    },
    "verification": {
        "nourishment_ok": BN_REPORT["ok"],
        "roads_drowned": len(_drowned),
        "roads_reloc_blocked": len(_blocked),
    },
}

if GROIN_CALLBACK is not None:
    _META["groin"].update({
        "kind": GROIN_KIND,
        **({"trapping_rate_m_yr": GROIN_CALLBACK.M} if GROIN_KIND == "dipole"
           else {"blocking_fraction": GROIN_CALLBACK.blocking_fraction,
                 "mean_trapping_rate_m_yr":
                     f"{GROIN_CALLBACK.mean_trapping_rate_m_yr:.2f}"}),
        "updrift / downdrift": f"GIS {GROIN_UPDRIFT_GIS} / {GROIN_DOWNDRIFT_GIS}",
        "install_year": GROIN_CALLBACK.install_year,
        "deterioration": f"{GROIN_CALLBACK.deterioration_mode}, "
                         f"floor {GROIN_CALLBACK.deterioration_fraction}",
        "r_ipl_t0": f"{R_IPL:.4f}",
        "r_ipl_theta_deg": f"{R_IPL_THETA_DEG:.1f}",
        "predicted_extent_m": ("n/a" if GROIN_PREDICTED_EXTENT_M is None
                               else f"{GROIN_PREDICTED_EXTENT_M:.0f}"),
    })
    if GROIN_EXTENT is not None:
        _META["groin"]["measured_extent_m"] = (
            f"{GROIN_EXTENT['updrift_m']:.0f} updrift / "
            f"{GROIN_EXTENT['downdrift_m']:.0f} downdrift")

# Cross-run index: one row per run, for comparing the matrix

# One row per run, rebuilt from disk; at OUTPUT_ROOT so one table covers both periods
RUN_INDEX_PATH = OUTPUT_ROOT / RUN_INDEX_FILENAME
_index_row = {
    "run_name": RUN_NAME,
    "timestamp": _TIMESTAMP,
    "start_year": START_YEAR,
    "end_year": END_YEAR,
    "source_sink_preset": SOURCE_SINK_PRESET,
    "be_nonzero_domains": _BE_NONZERO,
    "be_values_digest": _BE_DIGEST,
    **{f"be_{_key}": _value for _key, _value in _BE_EDGE_RATES.items()},
    "scenario": SCENARIO,
    "scenario_overridden": bool(_SCENARIO_DEPARTURES),
    "groin_enabled": GROIN_ENABLED,
    "groin_kind": GROIN_KIND if GROIN_CALLBACK is not None else "",
    # A blocking groin's M is emergent: the mean rate it actually applied.
    "groin_trapping_m_yr": (np.nan if GROIN_CALLBACK is None
                            else GROIN_CALLBACK.M if GROIN_KIND == "dipole"
                            else GROIN_CALLBACK.mean_trapping_rate_m_yr),
    "groin_blocking_b": (GROIN_CALLBACK.blocking_fraction
                         if GROIN_CALLBACK is not None and GROIN_KIND == "blocking"
                         else np.nan),
    # f, recorded since 2026-08-31 (the seed runs used 0.9); be1 is be_rate_gis1_m_yr
    "groin_deterioration_f": (GROIN_CALLBACK.deterioration_fraction
                              if GROIN_CALLBACK is not None else np.nan),
    "roadway_management": ENABLE_ROADWAY_MANAGEMENT,
    "relocations_enabled": ENABLE_HISTORICAL_ROAD_RELOCATIONS,
    "beach_dune_management": ENABLE_BEACH_DUNE_MANAGEMENT,
    "nourishment_fills": ENABLE_NOURISHMENT_FILLS,
    "bdm_domains": int(sum(BEACH_DUNE_MANAGEMENT_ON)),
    "nourishment_projects": len(BN_SCHEDULE_APPLIED.projects),
    "Hs_m": Hs,
    # Where the run is filed, spelled out
    "kind": RUN_KIND,
    "tag": RUN_TAG,
    "sandbags_on": ENABLE_SANDBAG_PLACEMENT,
    "rslr_m_yr": SEA_LEVEL_RISE_RATE,
    "annual_states": _states,
    "run_complete": _states == RUN_YEARS + 1,
    "rate_estimator": RATE_ESTIMATOR,
    "mean_bias_m_yr": SKILL["mean_bias_m_yr"],
    "rmse_m_yr": SKILL["rmse_m_yr"],
    "mean_bias_interior_m_yr": SKILL["mean_bias_interior_m_yr"],
    "rmse_interior_m_yr": SKILL["rmse_interior_m_yr"],
    "lrr_r2_median": float(np.nanmedian(_lrr_r2_real)),
    "lrr_r2_below_floor": len(_poor_fit),
    "endpoint_mean_bias_interior_m_yr":
        SKILL_ENDPOINT["mean_bias_interior_m_yr"],
    "endpoint_rmse_interior_m_yr": SKILL_ENDPOINT["rmse_interior_m_yr"],
    "nourishment_ok": BN_REPORT["ok"],
    "roads_drowned": len(_drowned),
    "roads_reloc_blocked": len(_blocked),
    "groin_extent_updrift_m": (GROIN_EXTENT["updrift_m"]
                               if GROIN_EXTENT is not None else np.nan),
    "groin_extent_downdrift_m": (GROIN_EXTENT["downdrift_m"]
                                 if GROIN_EXTENT is not None else np.nan),
    "runtime_min": round(RUN_SECONDS / 60, 1),
    "use_sandbox_cascade": USE_SANDBOX_CASCADE,
    "topo_product": TOPO_PRODUCT,          # see the note at the json write
    "topo_dune_version": TOPO_DUNE_VERSION,
    "island_offset_version": ISLAND_OFFSET_VERSION,
    "geometry": HATTERAS_GEOMETRY,
    "rmse_reach_interior_m_yr": SKILL["rmse_reach_interior_m_yr"],
    "git_commit": _GIT["commit"][:12],
    "git_dirty": _GIT["dirty"],
    "barrier3d_commit": str(_B3D["commit"])[:12],
    "barrier3d_gap_momentum_fix": _B3D.get("gap_momentum_fix"),
    "dune_ceiling": ("per-cell" if getattr(cascade.barrier3d[0], "_DuneCeilingFromStart", False)
                     else "uniform"),
    "storm_file": STORM_FILE.name,
    "bdm_dune_cap": getattr(_bdm_module, "DUNE_CAP_APPLIES_TO", "whole dune cell"),
    "barrier3d_route_overwash_fix": _B3D["route_overwash_fix"],
}
# The row goes into the metadata; run_index.csv is rebuilt from every run's metadata
_META[INDEX_SECTION] = _index_row
_meta_txt, _meta_json = write_run_metadata(
    RUN_DIR, RUN_NAME, _META,
    header=[f"CASCADE run metadata -- generated by {GENERATED_BY}, "
            f"section 12.",
            "Companion .json holds the same values, machine-readable."])
print(f"                      {_meta_txt.name}")
print(f"                      {_meta_json.name}")
RUN_INDEX = rebuild_run_index(OUTPUT_ROOT, RUN_INDEX_PATH,
                              current_versions=current_topo_versions())
print(f"                      {RUN_INDEX_FILENAME}  "
      f"({len(RUN_INDEX)} runs indexed; rebuilt from every run's metadata)")


# 12.5 Figures

# Section 9's bundles, or the figures lose their geography

plot_rate_comparison(
    PLOTTED_RATE, cs_series, run,
    real_domains_only=PLOT_REAL_DOMAINS_ONLY, estimator=RATE_ESTIMATOR,
    sea_level_rise_rate_m_yr=SEA_LEVEL_RISE_RATE,
    save_path=str(write_path(
        RUN_DIR,
        "figure_rate" if PLOT_REAL_DOMAINS_ONLY else "figure_rate_buffers",
        RUN_NAME)),
    show=SHOW_FIGURES, **RATE_FIG_KWARGS)

plot_annotated_rate_comparison(
    PLOTTED_RATE, cs_series, run,
    estimator=RATE_ESTIMATOR,
    sea_level_rise_rate_m_yr=SEA_LEVEL_RISE_RATE,
    save_path=str(write_path(RUN_DIR, "figure_rate_buffers", RUN_NAME)),
    show=SHOW_FIGURES, **RATE_FIG_KWARGS)

# Section 9.4's target, now that the run has a year 0; buffers are NaN
SHORELINE_TARGET_M, OBSERVED_CHANGE_M = build_shoreline_target(
    shoreline_m[0], START_YEAR, END_YEAR, HATTERAS_DOMAINS, RAW_OFFSET_DIR)

reports.target_misfit_report(
    target_m=SHORELINE_TARGET_M, observed_change_m=OBSERVED_CHANGE_M,
    shoreline_m=shoreline_m, end_year=END_YEAR, geometry=HATTERAS_DOMAINS,
    raw_offset_dir=RAW_OFFSET_DIR)

# Road relocations, so the GIFs mark the year the road moved; None draws no markers
_RELOCATION_EVENTS = None
if getattr(cascade, "_roadways", None):
    _RELOCATION_EVENTS = np.zeros(
        (shoreline_m.shape[0], HATTERAS_DOMAINS.total_domains), dtype=bool)
    for _pad, _roadway in enumerate(cascade._roadways):
        # Not `getattr(...) or []`: the attribute is a numpy array
        _raw = getattr(_roadway, "_road_relocated_TS", None)
        _series = (np.asarray([], dtype=float) if _raw is None
                   else np.asarray(_raw, dtype=float))
        if _series.size:
            _n = min(_series.size, _RELOCATION_EVENTS.shape[0])
            _RELOCATION_EVENTS[:_n, _pad] = _series[:_n] > 0
    print(f"  roadway relocations   {int(_RELOCATION_EVENTS.sum())} events "
          f"across {int(_RELOCATION_EVENTS.any(axis=0).sum())} domains")

GIF_PATHS = make_all_shoreline_gifs(
    shoreline_m, run, GIF_JOBS,
    baseline_npy=GIF_BASELINE_NPY,
    target_m=SHORELINE_TARGET_M,
    relocations=_RELOCATION_EVENTS, **GIF_KWARGS)

plt.close("all")
print(f"\ndone                  {RUN_DIR}")
