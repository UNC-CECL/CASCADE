#!/usr/bin/env python3
"""HATTERAS ISLAND -- CASCADE hindcast runner (1984-2004 / 2004-2024).

Headless twin of `HAT_hindcast_1984_2024.ipynb`. The notebook is the source of
truth: this file mirrors it section for section, in the same order, calling the
same `cascade_pipeline` / `hatteras_site_config` code with the same values, so
the two produce the same run.

WHAT THAT MEANS IN PRACTICE

  - Nothing is duplicated that the packages already own. The background-erosion
    presets, road events, nourishment projects, community zones and annotations
    all come from `hatteras_site_config`; the loaders, audits, verifiers and
    plotting come from `cascade_pipeline`. An earlier version of this file kept
    its own copies of all of it, and they drifted.
  - The functions the notebook defines in its own cells (`build_cascade`,
    `run_cascade_simulation`, `build_target_table`, `brie_r_ipl`,
    `measure_groin_extent`, ...) are copied here verbatim. Edit one, edit both.
  - Executed top to bottom at module level, which is what a notebook does.
    Nothing imports this file.

THE TWO DELIBERATE DIFFERENCES FROM THE NOTEBOOK

  - The notebook's QC plots are not here: island orientation, initialization
    plan view, RSLR, storms, source/sink, roadway plan view, beach/dune
    footprints, groin setting, CoastSat target. They display, they do not feed
    the run. Every *check* that sits beside them is kept -- the Barrier3D units
    contract, the setback audit, the nourishment audit, the annotation guard.
  - `SHOW_FIGURES` (below) defaults to False. The notebook renders its final
    figures inline with `show=True`; a headless run cannot. The files written
    to disk are identical either way.

WHERE THE SWITCHES LIVE  (same sections as the notebook)

  section 1   USE_SANDBOX_CASCADE, SHOW_FIGURES
  section 3   START_YEAR, SOURCE_SINK_PRESET, SCENARIO, GROIN_ENABLED
  section 11  OVERWRITE
              SCENARIO expands to ENABLE_ROADWAY_MANAGEMENT,
              ENABLE_BEACH_DUNE_MANAGEMENT, ENABLE_NOURISHMENT_FILLS and
              ENABLE_HISTORICAL_ROAD_RELOCATIONS; commented override lines
              sit directly below the table
  section 7   GROIN_TRAPPING_RATE_M_YR, deterioration settings
  section 9   FLIP_SIGN_MODEL, PLOT_REAL_DOMAINS_ONLY, GIF_JOBS
  section 11  Hs, dune thresholds, datums, ENABLE_SANDBAG_PLACEMENT

RUN_NAME is derived from those switches in section 7.5 -- it is not typed by
hand, so the output directory cannot disagree with what was simulated.
"""

# =============================================================================
# 1. IMPORTS
# =============================================================================

import datetime
import os
import sys
import time
from pathlib import Path

# cascade_pipeline and hatteras_site_config live in scripts/, which isn't
# installed. The notebook walks up from cwd to find pyproject.toml; a script
# knows where it is.
_HERE = Path(__file__).resolve()
PROJECT_BASE_DIR = _HERE.parents[2]
SCRIPTS_DIR = PROJECT_BASE_DIR / "scripts"
# HAT_hindcast_config sits beside this file. Running the .py puts that
# directory on sys.path automatically; importing it from the notebook does
# not, so it is added explicitly and both files reach the module the same way.
HATTERAS_MS_DIR = _HERE.parent
if not (PROJECT_BASE_DIR / "pyproject.toml").exists():
    raise RuntimeError(
        f"CASCADE repo root not found: {PROJECT_BASE_DIR} has no "
        f"pyproject.toml. This file expects to live in scripts/hatteras_ms/.")
for _path in (SCRIPTS_DIR, HATTERAS_MS_DIR):
    if str(_path) not in sys.path:
        sys.path.insert(0, str(_path))

import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import yaml

# True: sandbox copy with the pre-AST groin hook. False: real package, hook
# folded in.
USE_SANDBOX_CASCADE = True
if USE_SANDBOX_CASCADE:
    from cascade.cascade_groin import Cascade
else:
    from cascade.cascade import Cascade

from cascade.groin import GroinCallback, predict_fillet

from cascade_pipeline import nourishment, roadway
from cascade_pipeline import roadway as roadway_module
from cascade_pipeline.coastsat_loess import (
    CoastSatDataset,
    LoessConfig,
    build_coastsat_series,
    compute_domain_means,
)
from cascade_pipeline.plotting.rate_comparison import (
    DEFAULT_RATE_COMPARISON,
    plot_annotated_rate_comparison,
    plot_rate_comparison,
)
from cascade_pipeline.plotting.shoreline_gif import GifConfig, make_all_shoreline_gifs
from cascade_pipeline.run_info import RunInfo
from cascade_pipeline.run_registry import (
    append_run_index,
    git_provenance,
    guard_run_dir,
    run_dir_contents,
    skill_vs_target,
    values_digest,
    timestamp,
    RUN_INDEX_FILENAME,
    write_run_metadata,
)
from cascade_pipeline.shoreline import build_shoreline_matrix, compute_change_rate

from hatteras_site_config import (
    HATTERAS_ANNOTATIONS,
    HATTERAS_BEACH_DUNE,
    HATTERAS_BE_PRESETS,
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
    resolve_be_preset,
)

# The notebook draws its final figures inline. A headless run cannot, so the
# figures are saved and not shown. This is the only behavioural difference in
# the output path; the files on disk are identical.
SHOW_FIGURES = False
if not SHOW_FIGURES:
    matplotlib.use("Agg")

try:
    from tqdm.auto import tqdm
except ImportError:                                     # optional dependency
    tqdm = None

print(f"Imports OK from {SCRIPTS_DIR}")
print(f"USE_SANDBOX_CASCADE = {USE_SANDBOX_CASCADE}")
print(f"HATTERAS_DOMAINS.total_domains = {HATTERAS_DOMAINS.total_domains}")


# =============================================================================
# 2. FIXED DUNE/TOPO (period-independent)
# =============================================================================
# Same 2009 init surface for both periods, doesn't depend on START_YEAR.

# --- 2.1 project paths -------------------------------------------------------

TOPO_DUNE_INIT_YEAR = "2009"     # year label inside the filenames
TOPO_DUNE_VERSION = "2009_v3"    # extractor version folder
#   v1/v2 were moved to 2009-dune-topo/incorrect/ on 2026-08-17, so pinning
#   either fails the init-file check below. v3 is also the extraction the
#   dune-start road setbacks are measured against -- keep the two in step.

HATTERAS_DATA_BASE = PROJECT_BASE_DIR / "data" / "hatteras_init"
OUTPUT_ROOT = PROJECT_BASE_DIR / "output" / "raw_runs"
COASTSAT_BASE_DIR = (PROJECT_BASE_DIR / "scripts" / "input_prep"
                     / "5-scr" / "CoastSat")
PARAMETER_FILE = "Hatteras-CASCADE-parameters.yaml"  # resolved by CASCADE

BARRIER3D_DIR = HATTERAS_DATA_BASE / "1-barrier3d-domains"
DUNE_TOPO_DIR = BARRIER3D_DIR / "2009-dune-topo" / TOPO_DUNE_VERSION
BUFFER_DIR = BARRIER3D_DIR / "2009-buffer"

os.chdir(PROJECT_BASE_DIR)
# Runs are filed per period: section 3 builds OUTPUT_BASE_DIR =
# OUTPUT_ROOT / "<start>_<end>" once START_YEAR is expanded, so a
# 1984-2004 run and a 2004-2024 run of one scenario cannot land beside
# each other. run_index.csv stays at the root, covering both periods.
OUTPUT_ROOT.mkdir(parents=True, exist_ok=True)

for _name, _path in [
    ("PROJECT_BASE_DIR", PROJECT_BASE_DIR),
    ("HATTERAS_DATA_BASE", HATTERAS_DATA_BASE),
    ("DUNE_TOPO_DIR", DUNE_TOPO_DIR),
    ("BUFFER_DIR", BUFFER_DIR),
    ("COASTSAT_BASE_DIR", COASTSAT_BASE_DIR),
    ("OUTPUT_ROOT", OUTPUT_ROOT),
]:
    print(f"{_name:<20} {'ok' if _path.exists() else 'MISSING':<8} {_path}")


# --- 2.2 build the padded file lists -----------------------------------------

def build_domain_file_paths(geometry, dune_topo_dir, buffer_dir, init_year):
    """Builds one elevation and one dune file path per padded domain.

    Args:
        geometry: DomainGeometry describing the padded domain array.
        dune_topo_dir: Directory holding the topography/ and dunes/ subdirs.
        buffer_dir: Directory holding the buffer sample profiles.
        init_year: Year label embedded in the real-domain filenames.

    Returns:
        An (elevation_paths, dune_paths) tuple of string lists. Each list is
        geometry.total_domains long and index-aligned with the padded array.
    """
    buffer_elevation = str(buffer_dir / "sample_1_topography.npy")
    buffer_dune = str(buffer_dir / "sample_1_dune.npy")

    elevation_paths = [buffer_elevation] * geometry.num_buffer_domains
    dune_paths = [buffer_dune] * geometry.num_buffer_domains

    for gis_id in range(geometry.first_gis_id, geometry.last_gis_id + 1):
        elevation_paths.append(str(
            dune_topo_dir / "topography"
            / f"domain_{gis_id}_topography_{init_year}.npy"))
        dune_paths.append(str(
            dune_topo_dir / "dunes"
            / f"domain_{gis_id}_dune_{init_year}.npy"))

    elevation_paths += [buffer_elevation] * geometry.num_buffer_domains
    dune_paths += [buffer_dune] * geometry.num_buffer_domains

    return elevation_paths, dune_paths


ELEVATION_FILE_PATHS, DUNE_FILE_PATHS = build_domain_file_paths(
    HATTERAS_DOMAINS, DUNE_TOPO_DIR, BUFFER_DIR, TOPO_DUNE_INIT_YEAR)

print(f"\n{len(ELEVATION_FILE_PATHS)} elevation + {len(DUNE_FILE_PATHS)} dune "
      f"paths (expect {HATTERAS_DOMAINS.total_domains} each)")


# --- 2.3 verify every file exists --------------------------------------------
# A stale TOPO_DUNE_VERSION or a moved data folder fails here with a count and
# the first offender, rather than as an opaque traceback inside Barrier3D.

_expected_files = 2 * HATTERAS_DOMAINS.total_domains
_missing = [path for path in ELEVATION_FILE_PATHS + DUNE_FILE_PATHS
            if not Path(path).exists()]

if _missing:
    raise FileNotFoundError(
        f"{len(_missing)} of {_expected_files} init files missing. Check "
        f"TOPO_DUNE_VERSION ({TOPO_DUNE_VERSION!r}) and DUNE_TOPO_DIR.\n"
        f"  First missing: {_missing[0]}")

print(f"All {_expected_files} init files present.")


# --- 2.4 units check against Barrier3D's input contract ----------------------
# load_input.py converts the SCALAR yaml parameters but loads the elevation and
# dune .npy files verbatim, so the arrays on disk must already be in decameters
# relative to MHW. A file written in metres would run without error and model an
# island 10x too tall. Checked on the RAW arrays, before any DAM_TO_M scaling.

DAM_TO_M = 10.0          # Barrier3D works in decameters
WATER_CLAMP_DAM = -0.3   # SENTINEL_WATER_M / WATER_CLAMP_M from RUN_MANIFEST.txt
MAX_PLAUSIBLE_DAM = 2.0  # 20 m; a barrier island in metres would blow past this


def load_barrier3d_contract(parameter_path):
    """Reads the unit-relevant Barrier3D parameters from the CASCADE YAML.

    Mirrors the conversions in barrier3d/load_input.py so the values can be
    compared against the raw arrays on disk.

    Args:
        parameter_path: Path to the CASCADE parameter YAML.

    Returns:
        A dict with barrier_length_cells, mhw_dam and berm_el_dam.
    """
    with open(parameter_path, encoding="utf-8") as handle:
        params = yaml.safe_load(handle)

    mhw_dam = params["MHW"] / 10.0
    return {
        "barrier_length_cells": int(params["BarrierLength"] / 10.0),
        "mhw_dam": mhw_dam,
        "berm_el_dam": params["BermEl"] / 10.0 - mhw_dam,
    }


def check_domain_units(elevation_dam, dune_dam, contract):
    """Checks one domain's raw arrays against the Barrier3D input contract.

    Args:
        elevation_dam: Raw elevation array as stored on disk.
        dune_dam: Raw dune array as stored on disk.
        contract: Mapping from load_barrier3d_contract.

    Returns:
        A dict of check name -> bool, True when the array satisfies the check.
    """
    alongshore_cells = (elevation_dam.shape[1] if elevation_dam.ndim == 2
                        else None)
    return {
        "elevation is 2-D (cross_shore, alongshore)": elevation_dam.ndim == 2,
        "dune is 1-D (one crest per alongshore cell)": dune_dam.ndim == 1,
        "alongshore cells match dune length":
            alongshore_cells == dune_dam.size,
        f"alongshore >= BarrierLength ({contract['barrier_length_cells']})":
            alongshore_cells is not None
            and alongshore_cells >= contract["barrier_length_cells"],
        "dtype is floating point":
            np.issubdtype(elevation_dam.dtype, np.floating)
            and np.issubdtype(dune_dam.dtype, np.floating),
        "magnitudes are decameters, not metres":
            np.abs(elevation_dam).max() <= MAX_PLAUSIBLE_DAM,
        f"seaward floor at water clamp ({WATER_CLAMP_DAM} dam)":
            np.isclose(elevation_dam.min(), WATER_CLAMP_DAM),
        "dune heights positive (above berm)": dune_dam.min() > 0,
    }


BARRIER3D_CONTRACT = load_barrier3d_contract(HATTERAS_DATA_BASE / PARAMETER_FILE)

print(f"\nContract from {PARAMETER_FILE}:")
print(f"  BarrierLength -> {BARRIER3D_CONTRACT['barrier_length_cells']} "
      f"alongshore cells")
print(f"  MHW           -> {BARRIER3D_CONTRACT['mhw_dam']:.3f} dam")
print(f"  BermEl        -> {BARRIER3D_CONTRACT['berm_el_dam']:.3f} dam "
      f"above MHW\n")

_results = {}
for _gis_id in range(HATTERAS_DOMAINS.first_gis_id,
                     HATTERAS_DOMAINS.last_gis_id + 1):
    _pad = HATTERAS_DOMAINS.gis_to_pad(_gis_id)
    _checks = check_domain_units(np.load(ELEVATION_FILE_PATHS[_pad]),
                                 np.load(DUNE_FILE_PATHS[_pad]),
                                 BARRIER3D_CONTRACT)
    for _name, _passed in _checks.items():
        _results.setdefault(_name, []).append((_gis_id, _passed))

print(f"{'Check':<46} {'domains passing':>15}")
_failed_any = False
for _name, _outcomes in _results.items():
    _failures = [gis_id for gis_id, passed in _outcomes if not passed]
    _failed_any = _failed_any or bool(_failures)
    print(f"  {_name:<44} {len(_outcomes) - len(_failures):>6}"
          f"/{len(_outcomes)}"
          + (f"   FAILED: {_failures[:8]}" if _failures else ""))

print("\n" + ("SOME CHECKS FAILED - do not run the model on these arrays"
               if _failed_any else
               "All units and shapes match what Barrier3D expects."))


# =============================================================================
# 3. ISLAND ORIENTATION -- SET START_YEAR
# =============================================================================
# START_YEAR is the one flip. It selects a period from HATTERAS_PERIODS and
# everything in section 4 follows from it: run length, RSLR rate, storm series,
# background-erosion preset, road-setback file, nourishment settings.
#
# The run-selecting values in this section come from HAT_hindcast_config, so a
# driver can set them through the environment without editing tracked source.
# Interactively, type over any assignment below -- the imported value is only
# the default, and the assignment is still the last word.

from HAT_hindcast_config import RUN_CONFIG            # noqa: E402
from HAT_hindcast_config import describe as _describe_run_config  # noqa: E402

START_YEAR = RUN_CONFIG.start_year   # 1984 or 2004

# The source/sink axis of the run matrix. Each name states a hypothesis about
# where the alongshore sediment budget is unresolved:
#   "zeroBE"   nothing imposed anywhere
#   "edgeBE"   only the two end domains, absorbing the open-boundary artifact
#   "calibBE"  the full per-domain fit against the CoastSat LRR
SOURCE_SINK_PRESET = RUN_CONFIG.source_sink_preset

# The management axis of the run matrix. Every switch that decides whether a
# CASCADE management module runs at all lives here, in one block, because a
# switch kept beside the section that uses it is how a run ends up half-managed
# and named for a scenario it did not simulate. Sections 5 and 6 read these;
# they do not define their own.
#
#   ENABLE_ROADWAY_MANAGEMENT     roadway_manager: bulldozing, dune rebuild to
#                                 the design elevation, setback tracking, road
#                                 drowning. Off leaves the road as forcing that
#                                 nothing acts on.
#   ENABLE_BEACH_DUNE_MANAGEMENT  beach_dune_manager, module and all: overwash
#                                 filtering, the fixed dune line
#                                 (dune_migration_on = False), the 50 m
#                                 community-width drowning check, and fills.
#                                 Off is the only way to get natural shoreline
#                                 behaviour in the village domains.
#   ENABLE_NOURISHMENT_FILLS      Within an enabled beach_dune_manager, whether
#                                 historical fill is actually spent. False
#                                 leaves the module and its footprint exactly
#                                 as they were, so a fills-on / fills-off pair
#                                 differs in the fill and nothing else.
#
# Both modules off is the natural-dynamics run. The groin is deliberately not
# in this block: it is a structure, not a management module, and keeps its own
# switch in 7.1 so groin-only scenarios stay reachable.
# =============================================================================
# SCENARIO -- the management combination this run simulates
# =============================================================================
# The management axis of the run matrix, named rather than typed as four
# booleans, so a scenario is chosen in one word and cannot be assembled wrong
# by accident. What each switch controls:
#
#   roadway      roadway_manager: bulldozing, dune rebuild to the design
#                elevation, setback tracking, road drowning. Off leaves NC-12
#                as forcing that nothing acts on -- still loaded, still
#                audited in 5.1, never handed to a RoadwayManager.
#   beach_dune   beach_dune_manager, module and all: overwash filtering, the
#                fixed dune line (dune_migration_on = False), the 50 m
#                community-width drowning check, and fills. Off is the ONLY
#                way to get natural shoreline behaviour in the village
#                domains -- see the section 6 markdown on what is always-on.
#   fills        Within an enabled beach_dune_manager, whether the historical
#                fill is actually spent. False leaves the module and its
#                footprint exactly as they are, so a fills-on / fills-off pair
#                differs in the fill and nothing else.
#   relocations  Historical NC-12 relocation events. False in every named
#                scenario; read 5.1 before overriding it on.
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

# The groin is NOT part of the scenario, deliberately. 12.3 measures its effect
# against a paired no-groin baseline identical in every other token, so every
# scenario is run twice -- False first to create the baseline, then True.
# Folding it into the table would double the table to say the same thing.
GROIN_ENABLED = RUN_CONFIG.groin_enabled

print("\n" + _describe_run_config())

# --- expand ------------------------------------------------------------------
if SCENARIO not in SCENARIOS:
    raise ValueError(f"SCENARIO must be one of {sorted(SCENARIOS)}, "
                     f"got {SCENARIO!r}")
_SCENARIO_PRESET = SCENARIOS[SCENARIO]
ENABLE_ROADWAY_MANAGEMENT = _SCENARIO_PRESET["roadway"]
ENABLE_BEACH_DUNE_MANAGEMENT = _SCENARIO_PRESET["beach_dune"]
ENABLE_NOURISHMENT_FILLS = _SCENARIO_PRESET["fills"]
ENABLE_HISTORICAL_ROAD_RELOCATIONS = _SCENARIO_PRESET["relocations"]

# --- one-off overrides -------------------------------------------------------
# Uncomment to depart from the named scenario for a single run. The departure
# is detected and printed below, and the run name is still derived from the
# switches rather than from SCENARIO, so an overridden run cannot be filed
# under the scenario label it departed from.
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

# Normalised to the canonical key. The deprecated aliases ("base",
# "calibrated") still run, but it is the canonical name that reaches RUN_NAME
# in 7.5 -- an alias can never put a stale token in a directory name.
SOURCE_SINK_PRESET, _PRESET_BY_PERIOD = resolve_be_preset(SOURCE_SINK_PRESET)

PERIOD = HATTERAS_PERIODS[START_YEAR]
END_YEAR = PERIOD["end_year"]
RUN_YEARS = END_YEAR - START_YEAR

SEA_LEVEL_RISE_RATE = PERIOD["sea_level_rise_rate"]
ENABLE_NOURISHMENT = PERIOD["enable_nourishment"]
NOURISHMENT_VOLUME = PERIOD["nourishment_volume"]
ISLAND_OFFSET_FILE = HATTERAS_DATA_BASE / PERIOD["island_offset_file"]
STORM_FILE = HATTERAS_DATA_BASE / PERIOD["storm_file"]
ROAD_SETBACK_FILE = HATTERAS_DATA_BASE / PERIOD["road_setback_file"]

# --- resolve the combinations that cannot both be true -----------------------
# A fill cannot land in a domain with no BeachDuneManager: 6.3 would report it
# as dropped and the run would nourish nothing. A relocation event cannot move
# a setback nothing reads. Both are resolved here and announced below, rather
# than left as a contradiction for a later section to trip over.
_FILLS_FORCED_OFF = (ENABLE_NOURISHMENT_FILLS
                     and not ENABLE_BEACH_DUNE_MANAGEMENT)
if _FILLS_FORCED_OFF:
    ENABLE_NOURISHMENT_FILLS = False
_RELOCATIONS_FORCED_OFF = (ENABLE_HISTORICAL_ROAD_RELOCATIONS
                           and not ENABLE_ROADWAY_MANAGEMENT)
if _RELOCATIONS_FORCED_OFF:
    ENABLE_HISTORICAL_ROAD_RELOCATIONS = False

# Run name: the period stem only. The scenario suffix is derived from the active
# management switches in 7.5 -- a hand-typed label is how a groin-off run ends
# up in a directory named for a groin-on one.
RUN_NAME_STEM = f"HAT_{START_YEAR}_{END_YEAR}"

# Run directories are filed by period. Resolved here, not in section 1,
# because the period is not known until START_YEAR is expanded above.
# Every later section derives its paths from OUTPUT_BASE_DIR -- RUN_DIR
# in 9, the paired groin baseline in 12.3 -- so scoping it here scopes
# all of them, and the baseline lookup can no longer resolve to a run
# from the other period.
PERIOD_TAG = f"{START_YEAR}_{END_YEAR}"
OUTPUT_BASE_DIR = OUTPUT_ROOT / PERIOD_TAG
OUTPUT_BASE_DIR.mkdir(parents=True, exist_ok=True)

print(f"\nSTART_YEAR = {START_YEAR}  ->  {START_YEAR}-{END_YEAR}, "
      f"{RUN_YEARS} model years")
print(f"RUN_NAME_STEM = {RUN_NAME_STEM!r}"
      "   (scenario suffix derived in 7.5)")
print(f"SOURCE_SINK_PRESET = {SOURCE_SINK_PRESET!r}")
print(f"OUTPUT_BASE_DIR = {OUTPUT_BASE_DIR}")

# --- the name this scenario will produce, predicted from the switches --------
# Advisory only. 7.5 derives the authoritative RUN_NAME_BASE from what sections
# 5 and 6 actually built and raises if the two disagree, so this preview can
# never quietly become the thing that names the directory. Token order matches
# SCENARIO_SWITCHES in 7.5 exactly -- that is what makes the comparison valid.
# The period's fill is resolved with build_schedule, the same call 6 makes,
# rather than by re-implementing the date filter here.
_PERIOD_HAS_FILL = bool(nourishment.build_schedule(
    HATTERAS_NOURISHMENT_PROJECTS, HATTERAS_DOMAINS,
    START_YEAR, END_YEAR).projects)
_PREVIEW_TOKENS = [
    SOURCE_SINK_PRESET,
    "road" if ENABLE_ROADWAY_MANAGEMENT else "noroad",
    "reloc" if ENABLE_HISTORICAL_ROAD_RELOCATIONS else None,
    "bdm" if ENABLE_BEACH_DUNE_MANAGEMENT else "nobdm",
    ("nourish" if (ENABLE_NOURISHMENT_FILLS and _PERIOD_HAS_FILL)
     else ("nonourish" if _PERIOD_HAS_FILL and ENABLE_BEACH_DUNE_MANAGEMENT
           else None)),
    "groin" if GROIN_ENABLED else "nogroin",
]
RUN_NAME_PREVIEW = (f"{RUN_NAME_STEM}_"
                    + "_".join(t for t in _PREVIEW_TOKENS if t))

print(f"\nSCENARIO              {SCENARIO!r}"
      + (f"   OVERRIDDEN: {', '.join(sorted(_SCENARIO_DEPARTURES))}"
         if _SCENARIO_DEPARTURES else ""))
for _key, _preset in SCENARIOS.items():
    print(f"  {'->' if _key == SCENARIO else '  '} {_key:<18}"
          f"road={str(_preset['roadway']):<5} "
          f"bdm={str(_preset['beach_dune']):<5} "
          f"fills={str(_preset['fills'])}")
for _key, (_want, _got) in sorted(_SCENARIO_DEPARTURES.items()):
    print(f"  override            {_key}: scenario says {_want}, "
          f"this run uses {_got}")

print("\nMANAGEMENT            what the scenario above expands to")
print(f"  roadway_manager     "
      f"{'on' if ENABLE_ROADWAY_MANAGEMENT else 'OFF'}")
print(f"    relocations       "
      f"{'on' if ENABLE_HISTORICAL_ROAD_RELOCATIONS else 'off'}"
      + ("   FORCED off: roadway_manager is off"
         if _RELOCATIONS_FORCED_OFF else ""))
print(f"  beach_dune_manager  "
      f"{'on' if ENABLE_BEACH_DUNE_MANAGEMENT else 'OFF'}")
print(f"    nourishment fills "
      f"{'on' if ENABLE_NOURISHMENT_FILLS else 'off'}"
      + ("   FORCED off: beach_dune_manager is off"
         if _FILLS_FORCED_OFF else ""))
# PERIOD["enable_nourishment"] describes the period; nothing reads it. It is
# printed beside the switch that does reach the model so the two cannot quietly
# disagree -- the volume there is a legacy default that section 6 overwrites.
# Flagged in one direction only: a period that expects fill and is not getting
# it is a scenario choice worth restating, whereas fills left on in a period
# with no projects withholds nothing. 6.3 reports what was actually applied.
print(f"    period expects    {ENABLE_NOURISHMENT}"
      + ("   <- the period has fill and this run suppresses it"
         if ENABLE_NOURISHMENT and not ENABLE_NOURISHMENT_FILLS else ""))
print(f"  groin               {'on' if GROIN_ENABLED else 'off'}"
      f"   (paired baseline: run each scenario with False first)")

print(f"\nRUN_NAME (predicted)  {RUN_NAME_PREVIEW!r}")
print(f"                      confirmed against the built modules in 7.5")
print()
for _name, _path in [("island offset", ISLAND_OFFSET_FILE),
                     ("storms", STORM_FILE),
                     ("road setback", ROAD_SETBACK_FILE)]:
    print(f"  {_name:<14} {'ok' if _path.exists() else 'MISSING':<8} {_path.name}")


# --- 3.1 island offsets ------------------------------------------------------
# One value per padded domain giving that domain's cross-shore starting
# position. Becomes `shoreline_offset` on the Cascade() call.

def load_island_offset_dam(offset_path, geometry):
    """Loads a padded BRIE island-offset file and converts it to decameters.

    Args:
        offset_path: Path to a single-column padded offset CSV, in meters.
        geometry: DomainGeometry the file must match in length.

    Returns:
        A 1-D array of offsets in decameters, one per padded domain, ready to
        pass to Cascade() as shoreline_offset.

    Raises:
        ValueError: If the file is not one value per padded domain.
    """
    offset_m = np.loadtxt(offset_path, skiprows=1, delimiter=",")
    if offset_m.ndim != 1 or offset_m.size != geometry.total_domains:
        raise ValueError(f"{offset_path.name}: expected "
                         f"{geometry.total_domains} values, got shape "
                         f"{offset_m.shape}")
    return offset_m / DAM_TO_M


island_offset_dam = load_island_offset_dam(ISLAND_OFFSET_FILE, HATTERAS_DOMAINS)

_real = slice(HATTERAS_DOMAINS.start_real_index, HATTERAS_DOMAINS.end_real_index)
print(f"\n{START_YEAR} offsets: {island_offset_dam.size} padded domains | "
      f"real span {island_offset_dam[_real].min() * DAM_TO_M:.0f}-"
      f"{island_offset_dam[_real].max() * DAM_TO_M:.0f} m")


# =============================================================================
# 4. PERIOD FORCINGS -- RSLR, STORMS, SOURCE/SINK
# =============================================================================
# Everything here is resolved by the START_YEAR set in section 3.

# --- 4.1 relative sea level rise ---------------------------------------------

print(f"\nSEA_LEVEL_RISE_RATE = {SEA_LEVEL_RISE_RATE} m/yr")
print(f"  over {RUN_YEARS} years -> "
      f"{SEA_LEVEL_RISE_RATE * RUN_YEARS:.3f} m total rise")


# --- 4.2 storm series --------------------------------------------------------
# One row per storm: time (1-based model step), Rhigh, Rlow (decameters),
# period, duration (hours).

STORM_COLUMNS = ("time", "Rhigh", "Rlow", "period", "duration")


def load_storm_series(storm_path):
    """Loads a Barrier3D storm series into a DataFrame.

    Args:
        storm_path: Path to the storm .npy file.

    Returns:
        A DataFrame with the raw columns plus Rhigh_m and Rlow_m in meters.

    Raises:
        ValueError: If the file does not have the expected column count.
    """
    storms = np.load(storm_path)
    if storms.ndim != 2 or storms.shape[1] != len(STORM_COLUMNS):
        raise ValueError(f"{storm_path.name}: expected "
                         f"(n, {len(STORM_COLUMNS)}), got {storms.shape}")

    series = pd.DataFrame(storms, columns=STORM_COLUMNS)
    series["Rhigh_m"] = series["Rhigh"] * DAM_TO_M
    series["Rlow_m"] = series["Rlow"] * DAM_TO_M
    return series


STORM_SERIES = load_storm_series(STORM_FILE)

print(f"\n{STORM_FILE.name}: {len(STORM_SERIES)} storms")
print(f"  time steps    {STORM_SERIES['time'].min():.0f} to "
      f"{STORM_SERIES['time'].max():.0f}  (run is {RUN_YEARS} years)")
print(f"  Rhigh         {STORM_SERIES['Rhigh_m'].min():.2f} to "
      f"{STORM_SERIES['Rhigh_m'].max():.2f} m")
print(f"  duration      {STORM_SERIES['duration'].min():.0f} to "
      f"{STORM_SERIES['duration'].max():.0f} hours")
print(f"  storms/year   {len(STORM_SERIES) / RUN_YEARS:.1f} mean")

_beyond_run = (STORM_SERIES["time"] > RUN_YEARS + 1).sum()
if _beyond_run:
    print(f"  NOTE: {_beyond_run} storms land past model year {RUN_YEARS + 1}")


# --- 4.3 source/sink (background erosion) ------------------------------------
# Per-domain rate in m/yr, passed to Barrier3D as `Rat`. Sign convention from
# cascade/brie_coupler.py: (-) = erosion, (+) = accretion. Presets are sparse --
# a domain absent from a preset gets 0.0 m/yr.

def build_background_erosion(be_rates, geometry):
    """Expands sparse per-GIS background-erosion rates onto the padded array.

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


DOMAIN_BE_RATES = HATTERAS_BE_PRESETS[SOURCE_SINK_PRESET][START_YEAR]
BACKGROUND_EROSION_RATES = build_background_erosion(
    DOMAIN_BE_RATES, HATTERAS_DOMAINS)
USE_BACKGROUND_EROSION = any(rate != 0.0 for rate in BACKGROUND_EROSION_RATES)

# A derived fact, not a switch: the preset decides it. Checked because the two
# used to be separate knobs that could disagree, and the run name carried both.
_EXPECT_BE_ON = SOURCE_SINK_PRESET != "zeroBE"
if USE_BACKGROUND_EROSION != _EXPECT_BE_ON:
    raise ValueError(
        f"preset {SOURCE_SINK_PRESET!r} implies "
        f"USE_BACKGROUND_EROSION={_EXPECT_BE_ON}, but the expanded rates give "
        f"{USE_BACKGROUND_EROSION}. The preset in hatteras_site_config.py does "
        f"not match its name.")

print(f"\npreset {SOURCE_SINK_PRESET!r}, period {START_YEAR}: "
      f"{len(DOMAIN_BE_RATES)} domains specified, "
      f"{sum(1 for r in BACKGROUND_EROSION_RATES if r != 0)} non-zero of "
      f"{HATTERAS_DOMAINS.total_domains}")
if DOMAIN_BE_RATES:
    print(f"  rate range     {min(DOMAIN_BE_RATES.values()):+.1f} to "
          f"{max(DOMAIN_BE_RATES.values()):+.1f} m/yr")
else:
    print("  rate range     no domains specified -- 0.0 m/yr everywhere")
print(f"  USE_BACKGROUND_EROSION = {USE_BACKGROUND_EROSION}")
if (START_YEAR == 2004 and SOURCE_SINK_PRESET == "calibBE"
        and HATTERAS_BE_RATES_2004_IS_PLACEHOLDER):
    print("  WARNING: 2004 calibrated rates are a copy of the 1984 fit")


# =============================================================================
# 5. roadway_manager -- SETBACKS, PER-DOMAIN ELEVATION, HISTORICAL EVENTS
# =============================================================================
# NC-12's forcing is three things: where the road sits (setback, by period), how
# high it is (elevation, period-independent), and which domains are managed.
#
# Two conventions, both Barrier3D's:
#   - road_ele is metres MHW-relative, not NAVD88 -- bulldoze writes it straight
#     into the interior grid, which the extractor stores MHW-relative.
#   - relocation events carry a DISPLACEMENT, not an absolute setback. CASCADE
#     already decrements the setback by dune migration each year, so adding the
#     measured displacement counts the retreat once; an absolute setback
#     referenced to an older dune line counts it twice.

# ENABLE_ROADWAY_MANAGEMENT and ENABLE_HISTORICAL_ROAD_RELOCATIONS are set in
# section 3, with the other management switches. The forcing below is loaded
# either way: with management off the setbacks and elevations are still read
# and audited in 5.1, they simply never reach a RoadwayManager.

ROADWAY = roadway.RoadwayConfig()
_road_span = (HATTERAS_FIRST_ROAD_DOMAIN, HATTERAS_LAST_ROAD_DOMAIN)

# --- setback: by period ------------------------------------------------------
road_setbacks_full, _missing_setbacks = roadway.load_road_setbacks(
    ROAD_SETBACK_FILE, HATTERAS_DOMAINS, *_road_span)

# --- elevation: one set for every period -------------------------------------
ROAD_ELEVATION_FILE = HATTERAS_DATA_BASE / HATTERAS_ROAD_ELEVATION_FILE
road_elevation_full, _missing_elevations = roadway.load_road_elevations(
    ROAD_ELEVATION_FILE, HATTERAS_DOMAINS, *_road_span, config=ROADWAY)

# --- which domains CASCADE actually manages ----------------------------------
ROADWAY_MANAGEMENT_ON = roadway.build_roadway_management_on(
    HATTERAS_DOMAINS, *_road_span,
    community_zones=HATTERAS_COMMUNITY_ZONES,
    enabled=ENABLE_ROADWAY_MANAGEMENT)

_road_slice = slice(HATTERAS_DOMAINS.gis_to_pad(HATTERAS_FIRST_ROAD_DOMAIN),
                    HATTERAS_DOMAINS.gis_to_pad(HATTERAS_LAST_ROAD_DOMAIN) + 1)
print(f"\nsetback file    {ROAD_SETBACK_FILE.name}")
print(f"  {road_setbacks_full[_road_slice].min():.0f}"
      f"-{road_setbacks_full[_road_slice].max():.0f} m"
      + (f"   MISSING {_missing_setbacks}" if _missing_setbacks else ""))
print(f"elevation file  {ROAD_ELEVATION_FILE.name}  (period-independent)")
print(f"  {road_elevation_full[_road_slice].min():.2f}"
      f"-{road_elevation_full[_road_slice].max():.2f} m MHW"
      f"   buffers {ROADWAY.elevation_fallback_m} m"
      + (f"   MISSING {_missing_elevations}" if _missing_elevations else ""))
if ENABLE_ROADWAY_MANAGEMENT:
    print(f"managed         {sum(ROADWAY_MANAGEMENT_ON)} domains "
          f"({HATTERAS_LAST_ROAD_DOMAIN - HATTERAS_FIRST_ROAD_DOMAIN + 1} "
          f"carry road, minus community zones "
          f"{list(HATTERAS_COMMUNITY_ZONES)})")
else:
    print(f"managed         0 domains -- ENABLE_ROADWAY_MANAGEMENT is False "
          f"(section 3)")
    print(f"                no bulldozing, no dune rebuild, no setback "
          f"tracking, no road drowning")
for _event in HATTERAS_ROAD_EVENTS:
    print(f"event {_event.year}      {type(_event).__name__:<15} {_event.note}")
print(f"ENABLE_HISTORICAL_ROAD_RELOCATIONS = "
      f"{ENABLE_HISTORICAL_ROAD_RELOCATIONS}")


# --- 5.1 pre-flight audit: which road_offset will not survive year one -------------
# bulldoze tests the two rows FLANKING the road -- never the road's own cells --
# and drowns it when either flank is more than 20% water. A drowned road is not
# a warning: roadway_manager sets _drown_break and returns immediately on every
# later year, so the domain becomes an unmanaged barrier wearing a road label.

road_audit = roadway.audit_setbacks(
    ELEVATION_FILE_PATHS, road_setbacks_full, HATTERAS_DOMAINS, *_road_span,
    management_on=ROADWAY_MANAGEMENT_ON, config=ROADWAY)
audit_summary = roadway.summarise_audit(road_audit)

print(f"\n{audit_summary['n_domains']} road domains, "
      f"{audit_summary['n_managed']} managed\n")
if audit_summary["blocking"]:
    print("BLOCKING - these would corrupt the run rather than drown it:")
    for _gis, _wall in audit_summary["blocking"]:
        print(f"  GIS {_gis}: {_wall}")
    print()
print(f"drowns at t=0, managed   : {audit_summary['drowning']}")
print(f"drowns at t=0, unmanaged : {audit_summary['drowning_unmanaged']}"
      f"   (community zones; no manager runs there)\n")

_bad = [r for r in road_audit if r["drowns"] and r["managed"]]
if _bad:
    print(f"{'GIS':>4} {'setback':>9} {'road rows':>11} {'%sea':>6} "
          f"{'%bay':>6} {'%road wet':>10} {'land ends min/med':>18}")
    for _r in _bad:
        _rows = f"{_r['road_start']}-{_r['road_end'] - 1}"
        _ends = f"{_r['width_min']} / {_r['width_median']}"
        print(f"{_r['gis']:>4} {_r['setback_m']:>8.0f}m {_rows:>11} "
              f"{_r['sea_water'] * 100:>5.0f}% {_r['bay_water'] * 100:>5.0f}% "
              f"{_r['road_cells_water'] * 100:>9.0f}% {_ends:>18}")


# =============================================================================
# 6. beach_dune_manager -- NOURISHMENT SCHEDULE + OVERWASH FILTER
# =============================================================================
# Two different things arrive bundled in one CASCADE module.
#
# ALWAYS-ON wherever the module is enabled, every year: a percentage of overwash
# deposition is removed from the interior and returned to the shoreface, the
# dune line is held fixed (dune_migration_on = False), and the community is
# abandoned if the average interior width falls below 50 m.
#
# EVENT-DRIVEN, only where and when nourish_now says so: sand is added to the
# shoreface. There is no way to get the second without the first.
#
# Three conventions, all places this pipeline has been wrong before:
#   - overwash_filter is a PERCENT, not a fraction. filter_overwash divides by
#     100. A value of 0.4 filters 0.4% of overwash, which is indistinguishable
#     from no filtering. BeachDuneConfig now refuses the fraction scale.
#   - The per-year volume goes to cascade.nourishment_volume, which is where
#     CASCADE reads it from. A volume written onto the BeachDuneManager instance
#     lands on the attribute CASCADE overwrites one line before the manager
#     reads it, and the fill quietly spends the Cascade() init default.
#   - The manager's time series are offset by one; NourishmentSchedule.time_index
#     is the single place that conversion lives.
#
# The module footprint is the UNION of the permanent community zones and every
# domain that receives fill. Where that overlaps the roadway footprint, both
# managers run on the same domain -- reported here, verified in section 12.

# --- schedule: one project list, period-filtered ------------------------------
# Every Hatteras project falls in 2004-2024, so a 1984 run builds an empty
# schedule from this same list rather than needing a period-keyed one.
BN_SCHEDULE = nourishment.build_schedule(
    HATTERAS_NOURISHMENT_PROJECTS, HATTERAS_DOMAINS, START_YEAR, END_YEAR)

# --- the schedule the model is actually driven by -----------------------------
# ENABLE_NOURISHMENT_FILLS (section 3). Suppressing fills builds an EMPTY
# schedule for the same period rather than skipping apply_to_cascade: the loop
# still rewrites nourish_now and nourishment_volume to zero every year, and the
# audit below and the verification in 12.2 both check against what was driven
# rather than what was intended. BN_SCHEDULE itself is left intact and still
# defines the module footprint below, so turning fills off changes the fill and
# not the footprint.
BN_SCHEDULE_APPLIED = BN_SCHEDULE if ENABLE_NOURISHMENT_FILLS else (
    nourishment.build_schedule([], HATTERAS_DOMAINS, START_YEAR, END_YEAR))

# --- what CASCADE is handed ---------------------------------------------------
# The filter is inert with the module off -- only BeachDuneManager applies it --
# so it is built as zeros there rather than left at its community values. The
# array handed to CASCADE should state what the run does, not what it would do.
OVERWASH_FILTER = (
    nourishment.build_overwash_filter(
        HATTERAS_DOMAINS, HATTERAS_COMMUNITY_ZONES, config=HATTERAS_BEACH_DUNE)
    if ENABLE_BEACH_DUNE_MANAGEMENT
    else [0.0] * HATTERAS_DOMAINS.total_domains)
OVERWASH_TO_DUNE = HATTERAS_BEACH_DUNE.overwash_to_dune_pct
BEACH_DUNE_MANAGEMENT_ON = nourishment.build_beach_dune_management_on(
    HATTERAS_DOMAINS, HATTERAS_COMMUNITY_ZONES, BN_SCHEDULE.nourished_gis,
    enabled=ENABLE_BEACH_DUNE_MANAGEMENT)

# Placeholder for the Cascade() call. Every value is rewritten each model year
# by BN_SCHEDULE.apply_to_cascade(), so this is 0 rather than
# PERIOD["nourishment_volume"]: if the schedule ever fails to reach the model, a
# year should nourish nothing rather than quietly nourish the default.
NOURISHMENT_VOLUME_INIT = [0.0] * HATTERAS_DOMAINS.total_domains

DOUBLE_MANAGED_GIS = nourishment.find_double_managed(
    BEACH_DUNE_MANAGEMENT_ON, ROADWAY_MANAGEMENT_ON, HATTERAS_DOMAINS)
BN_AUDIT = nourishment.audit_schedule(
    BN_SCHEDULE_APPLIED, BEACH_DUNE_MANAGEMENT_ON, config=HATTERAS_BEACH_DUNE)

print(f"\nperiod                {START_YEAR}-{END_YEAR}")
if ENABLE_BEACH_DUNE_MANAGEMENT:
    print(f"overwash_filter       "
          f"{HATTERAS_BEACH_DUNE.community_overwash_filter_pct:.0f}% in "
          f"community zones, "
          f"{HATTERAS_BEACH_DUNE.default_overwash_filter_pct:.0f}% "
          f"elsewhere   (PERCENT, not a fraction)")
    print(f"overwash_to_dune      {OVERWASH_TO_DUNE:.0f}%")
    print(f"module on             {sum(BEACH_DUNE_MANAGEMENT_ON)} domains "
          f"= community zones {list(HATTERAS_COMMUNITY_ZONES)} "
          f"UNION nourished {list(BN_SCHEDULE.nourished_gis)}")
    if not ENABLE_NOURISHMENT_FILLS and BN_SCHEDULE.projects:
        print(f"                      the footprint still carries the project "
              f"extents: the fills are off, the module is not")
else:
    print(f"module on             0 domains -- "
          f"ENABLE_BEACH_DUNE_MANAGEMENT is False (section 3)")
    print(f"                      no overwash filtering, no fixed dune line, "
          f"no width-drowning check, no fill")
    print(f"                      overwash_filter handed to CASCADE is all "
          f"zero, so the array says what the run does")

print(f"\nnourishment schedule  "
      f"{len(BN_SCHEDULE_APPLIED.projects)} project(s) applied"
      + (f"   {len(BN_SCHEDULE.projects)} in period SUPPRESSED by "
         f"ENABLE_NOURISHMENT_FILLS=False"
         if not ENABLE_NOURISHMENT_FILLS and BN_SCHEDULE.projects else ""))
for _project in BN_SCHEDULE_APPLIED.projects:
    _per_m = _project.volume_m3_per_m(HATTERAS_DOMAINS.domain_spacing_m)
    print(f"  {_project.year}  {_project.name:<26} "
          f"GIS {_project.gis_domains[0]}-{_project.gis_domains[-1]:<3} "
          f"{_project.volume_cubic_yards:>9,.0f} cy -> "
          f"{_per_m:>7.1f} m^3/m per domain   "
          f"(TS index {BN_SCHEDULE_APPLIED.time_index(_project.year)})")
for _project, _reason in BN_SCHEDULE.skipped:
    print(f"  --    {_project.name:<26} skipped: {_reason}")

print(f"\npre-flight audit      {'ok' if BN_AUDIT['ok'] else 'PROBLEMS'}")
if BN_AUDIT["module_off"]:
    print(f"  scheduled but module OFF (fill would be dropped): "
          f"{BN_AUDIT['module_off']}")
for _row in BN_AUDIT["implausible_volume"]:
    print(f"  implausible volume: GIS {_row['gis']} {_row['year']} "
          f"{_row['volume_m3_per_m']:.1f} m^3/m")
for _note in BN_AUDIT["notes"]:
    print(f"  note: {_note}")

if DOUBLE_MANAGED_GIS:
    print(f"\nBOTH MODULES ON       GIS {list(DOUBLE_MANAGED_GIS)} "
          f"({len(DOUBLE_MANAGED_GIS)} domains)")
    print("  Overwash is removed twice over: RoadwayManager bulldozes and")
    print("  rewrites DomainTS, then BeachDuneManager filters what survived.")
    print("  BeachDuneManager also holds dune_migration_on False, so")
    print("  ShorelineChangeTS stays 0 -- the value RoadwayManager reads to")
    print("  decrement the setback. NC-12 stops retreating here.")
    print("  Run as published; verified against the finished run in section 12.")
elif not (ENABLE_ROADWAY_MANAGEMENT and ENABLE_BEACH_DUNE_MANAGEMENT):
    print("\nBOTH MODULES ON       none - at least one module is off, so the "
          "overlap cannot arise")
else:
    print("\nBOTH MODULES ON       none - the two footprints are disjoint "
          "this period")


# =============================================================================
# 7. hard_structures / GROIN -- BUXTON GROIN FIELD
# =============================================================================
# cascade/groin.py attaches through cascade._groin_callback and is called once
# per model year from inside Cascade.update(), immediately before the alongshore
# transport solve. Each active year it adds -M to the updrift domain and +M to
# the downdrift domain of x_s_dt. BRIE's implicit diffusion solve spreads that
# dipole in the same step, so the fillet's taper and extent are EMERGENT -- only
# its amplitude is imposed.
#
# It is a forcing, not a barrier: nothing here blocks alongshore transport. The
# pair is volume-neutral by construction, trapping never saturates on state, and
# install_year is inert (1969 precedes both periods).
#
# Two independent estimates of M disagree by an order of magnitude -- the
# shoreline-position fit says 50 m/yr, the sediment budget says that is
# unaffordable. Configured at the sweep's value; the breach is REPORTED, not
# corrected. Section 12 reports the resulting misfit.

# --- 7.1 switches, structure, sediment-budget reference ----------------------

# GROIN_ENABLED is set in section 3, beside the scenario, so one cell decides
# what a run simulates. The guard below stays here: this is where a wrongly
# configured groin would silently do nothing.

# The pre-AST hook exists ONLY in cascade/cascade_groin.py. Attaching a callback
# to a Cascade built from the real package is a silent no-op: the run succeeds,
# the groin does nothing, and the output looks like a valid groin-on run.
if GROIN_ENABLED and not USE_SANDBOX_CASCADE:
    raise RuntimeError(
        "GROIN_ENABLED=True requires USE_SANDBOX_CASCADE=True (section 1). "
        "cascade.cascade.Cascade.update() has no _groin_callback hook, so the "
        "groin would silently do nothing.")

# Net alongshore transport on Hatteras is southward, so updrift = north. The two
# domains must be adjacent -- they share the blocked boundary at GIS 5.5.
GROIN_UPDRIFT_GIS = 6       # source: accretes
GROIN_DOWNDRIFT_GIS = 5     # sink:   erodes
GROIN_INSTALL_YEAR = 1969   # confirmed construction date

# --- amplitude and deterioration floor: the two tunable knobs ----------------
# Both come from HAT_hindcast_config so the sweep driver can set them per run.
# They are fit JOINTLY across the two periods, not independently: over
# 1984-2004 the cumulative trapping is M*(16 + 4f), so f barely separates from
# M there, while over 2004-2024 the run sits entirely past the 2003 ramp and
# only the product M*f is identifiable. Neither window pins both on its own.
GROIN_TRAPPING_RATE_M_YR = RUN_CONFIG.groin_trapping_rate_m_yr
GROIN_M_PROVENANCE = ("joint two-period fit against the CoastSat D6-D5 "
                      "differential; see output/groin_sweep/ for the M-f "
                      "ridge and which grid bounds the solution touches")

# --- deterioration: 1996 last repair -> 2003 storm damage --------------------
GROIN_DETERIORATION_DELAY_YEARS = 1996 - GROIN_INSTALL_YEAR   # = 27
GROIN_DETERIORATION_MODE = "linear_ramp"
GROIN_DETERIORATION_RAMP_YEARS = 2003 - 1996                  # = 7
GROIN_DETERIORATION_FRACTION = RUN_CONFIG.groin_deterioration_fraction

# --- sediment-budget reference -----------------------------------------------
REACH_TRANSPORT_LOSS_M3_YR = 5.9e5
REACH_TRANSPORT_CITATION = ("Inman & Dolan (1989), via Moore et al. (2010), "
                            "doi:10.1029/2009JF001299")
REACH_TRANSPORT_CAVEAT = (
    "reach-integrated transport-gradient LOSS, Oregon Inlet to Cape Hatteras "
    "(~60 km) -- a divergence, not a gross flux at Buxton, so this bounds "
    "order of magnitude only")

# Active profile height, converting shoreline displacement to volume. Both
# candidates print here rather than one being asserted; section 11 resolves it
# from the constructed model.
GROIN_PROFILE_HEIGHT_CANDIDATES_M = (12.0, 24.0)


# --- 7.2 build the callback ---------------------------------------------------
# GROIN_CB is built unconditionally so 7.4's report renders either way.
# GROIN_CALLBACK is the one attached to the model in section 11.

GROIN_CB = GroinCallback(
    updrift_pad=HATTERAS_DOMAINS.gis_to_pad(GROIN_UPDRIFT_GIS),
    downdrift_pad=HATTERAS_DOMAINS.gis_to_pad(GROIN_DOWNDRIFT_GIS),
    trapping_rate_m_yr=GROIN_TRAPPING_RATE_M_YR,
    start_year=START_YEAR,
    install_year=GROIN_INSTALL_YEAR,
    n_domains=HATTERAS_DOMAINS.total_domains,
    deterioration_delay_years=GROIN_DETERIORATION_DELAY_YEARS,
    deterioration_mode=GROIN_DETERIORATION_MODE,
    deterioration_ramp_years=GROIN_DETERIORATION_RAMP_YEARS,
    deterioration_fraction=GROIN_DETERIORATION_FRACTION,
)
GROIN_CALLBACK = GROIN_CB if GROIN_ENABLED else None


# --- 7.3 diagnostics ----------------------------------------------------------

def groin_trapping_schedule(callback, start_year, end_year):
    """Effective trapping rate for every model year of a period.

    Reads the callback's deterioration curve without running it.
    `_effective_trapping_rate` is a pure function of the year -- it does not
    touch `_call_count` -- so querying it here leaves the callback's year
    counter at zero for the actual run.

    Args:
        callback: A GroinCallback.
        start_year: First model year.
        end_year: Last model year, exclusive.

    Returns:
        A (years, M_eff) tuple of 1-D arrays.
    """
    years = np.arange(start_year, end_year)
    m_eff = np.array([callback._effective_trapping_rate(int(y)) for y in years])
    return years, m_eff


def implied_interception_m3_yr(m_per_yr, profile_height_m, geometry):
    """Sediment volume the dipole transfers across the groin each year.

    A shoreline displacement of `m_per_yr` over one domain implies a volume of
    `m_per_yr * dy * (d_sf + h_b)`. The dipole moves that from the downdrift
    side to the updrift side; the transfer is volume-neutral overall.

    Args:
        m_per_yr: Trapping rate M, in meters per year.
        profile_height_m: Active profile height (d_sf + h_b), in meters.
        geometry: DomainGeometry supplying the alongshore domain width.

    Returns:
        The implied transfer in cubic meters per year.
    """
    return m_per_yr * geometry.domain_spacing_m * profile_height_m


# --- 7.4 report ---------------------------------------------------------------

_up_pad = HATTERAS_DOMAINS.gis_to_pad(GROIN_UPDRIFT_GIS)
_down_pad = HATTERAS_DOMAINS.gis_to_pad(GROIN_DOWNDRIFT_GIS)
_years, _m_eff = groin_trapping_schedule(GROIN_CB, START_YEAR, END_YEAR)

print(f"\ngroin                 {'ENABLED' if GROIN_ENABLED else 'disabled'}"
      f"   (module cascade.groin, hook cascade_groin.py pre-AST)")
print(f"structure             updrift GIS {GROIN_UPDRIFT_GIS} (pad {_up_pad}) "
      f"/ downdrift GIS {GROIN_DOWNDRIFT_GIS} (pad {_down_pad}), "
      f"installed {GROIN_INSTALL_YEAR}")
print(f"trapping rate M       {GROIN_TRAPPING_RATE_M_YR:.0f} m/yr")
print(f"  provenance          {GROIN_M_PROVENANCE}")
print(f"deterioration         {GROIN_DETERIORATION_MODE}, "
      f"onset {GROIN_CB.deterioration_year:.0f} "
      f"(+{GROIN_DETERIORATION_DELAY_YEARS:.0f} yr), "
      f"ramp {GROIN_DETERIORATION_RAMP_YEARS:.0f} yr, "
      f"floor {GROIN_DETERIORATION_FRACTION:.2f}")
print(f"  M over {START_YEAR}-{END_YEAR}    "
      f"{_m_eff[0]:.1f} -> {_m_eff[-1]:.1f} m/yr   "
      f"(mean {_m_eff.mean():.1f})")
if GROIN_INSTALL_YEAR <= START_YEAR:
    print(f"  note                install {GROIN_INSTALL_YEAR} predates the "
          f"run: active every step, and {START_YEAR - GROIN_INSTALL_YEAR} yr "
          f"of fillet is already in the initial shoreline")
if GROIN_DETERIORATION_FRACTION >= 0.5:
    print(f"  TENSION             floor {GROIN_DETERIORATION_FRACTION:.2f} "
          f"leaves the structure at "
          f"{GROIN_DETERIORATION_FRACTION * 100:.0f}% strength after failure; "
          f"the 1996/2003 record does not support that")

print(f"\nSEDIMENT BUDGET       reference {REACH_TRANSPORT_LOSS_M3_YR:,.0f} "
      f"m3/yr")
print(f"  source              {REACH_TRANSPORT_CITATION}")
print(f"  CAVEAT              {REACH_TRANSPORT_CAVEAT}")
print(f"  implied transfer at M = {GROIN_TRAPPING_RATE_M_YR:.0f} m/yr over "
      f"dy = {HATTERAS_DOMAINS.domain_spacing_m:.0f} m:")
for _h in GROIN_PROFILE_HEIGHT_CANDIDATES_M:
    _vol = implied_interception_m3_yr(
        GROIN_TRAPPING_RATE_M_YR, _h, HATTERAS_DOMAINS)
    _pct = 100.0 * _vol / REACH_TRANSPORT_LOSS_M3_YR
    print(f"    profile height {_h:5.1f} m ->  {_vol:>9,.0f} m3/yr "
          f"= {_pct:5.1f}% of the reach budget"
          + ("   BREACH" if _pct > 50 else ""))
print("  Not corrected. Section 11 resolves the profile height from the "
      "constructed model")
print("  rather than the yaml, and section 12 reports the resulting misfit.")

print(f"\nSOURCE/SINK OVERLAP   preset {SOURCE_SINK_PRESET!r}, "
      f"period {START_YEAR}   ((-) erosion / (+) accretion)")
for _label, _gis in (("updrift  ", GROIN_UPDRIFT_GIS),
                     ("downdrift", GROIN_DOWNDRIFT_GIS)):
    _be = DOMAIN_BE_RATES.get(_gis, 0.0)
    _dipole = -GROIN_TRAPPING_RATE_M_YR if _gis == GROIN_UPDRIFT_GIS \
        else +GROIN_TRAPPING_RATE_M_YR
    print(f"  {_label} GIS {_gis:<3} BE {_be:+6.2f} m/yr    "
          f"groin x_s_dt {_dipole:+7.1f} m/yr")
print("  The calibrated BE rates were fit to observed CoastSat LRR spanning "
      "the")
print("  functional-groin era, so the groin's signature is plausibly counted "
      "twice.")
print("  Reported pending the source/sink re-analysis; not corrected here.")


# --- 7.5 scenario summary, and the run name derived from it ------------------
# RUN_NAME_SUFFIX used to be typed by hand: run with the groin off, forget to
# retype the label, and the output lands in a directory named for a different
# experiment. It is now derived from the switches themselves.

SCENARIO_SWITCHES = [
    ("period", f"{START_YEAR}-{END_YEAR} ({RUN_YEARS} yr)", None),
    ("source/sink preset", SOURCE_SINK_PRESET, SOURCE_SINK_PRESET),
    # No token of its own: it is implied by the preset, and checked against it
    # in 4.3. Emitting both produced names like "..._base_noBE_..." that said
    # the same thing twice without saying which zero it was.
    ("background erosion", USE_BACKGROUND_EROSION, None),
    ("roadway management", f"{sum(ROADWAY_MANAGEMENT_ON)} domains"
     if ENABLE_ROADWAY_MANAGEMENT else "off",
     "road" if ENABLE_ROADWAY_MANAGEMENT else "noroad"),
    ("historical relocations", ENABLE_HISTORICAL_ROAD_RELOCATIONS,
     "reloc" if ENABLE_HISTORICAL_ROAD_RELOCATIONS else None),
    ("beach/dune manager", f"{sum(BEACH_DUNE_MANAGEMENT_ON)} domains"
     if ENABLE_BEACH_DUNE_MANAGEMENT else "off",
     "bdm" if ENABLE_BEACH_DUNE_MANAGEMENT else "nobdm"),
    # "nonourish" only when there was fill to withhold and a module to
    # withhold it in: with beach_dune_manager off, "nobdm" already says no
    # fill, and 1984 has no project to suppress. The project count is
    # deliberately not in the name -- it is a property of the period, not of
    # the scenario, and run_index.csv carries it as nourishment_projects.
    ("nourishment fills", f"{len(BN_SCHEDULE_APPLIED.projects)} applied"
     if ENABLE_NOURISHMENT_FILLS
     else ("suppressed" if BN_SCHEDULE.projects else "none in period"),
     "nourish" if BN_SCHEDULE_APPLIED.projects
     else ("nonourish" if BN_SCHEDULE.projects and ENABLE_BEACH_DUNE_MANAGEMENT
           else None)),
    ("groin", "on" if GROIN_ENABLED else "off",
     "groin" if GROIN_ENABLED else "nogroin"),
]

RUN_NAME_SUFFIX = "_".join(
    token for _, _, token in SCENARIO_SWITCHES if token)
RUN_NAME_BASE = f"{RUN_NAME_STEM}_{RUN_NAME_SUFFIX}"

# The section 3 preview was predicted from the switches; this name is derived
# from what sections 5 and 6 actually built. A difference means a switch did
# not reach the module it names -- the failure that produces a run filed under
# a scenario it did not simulate. Raised, not warned.
if RUN_NAME_BASE != RUN_NAME_PREVIEW:
    raise AssertionError(
        f"run name disagrees with the section 3 preview\n"
        f"  section 3    {RUN_NAME_PREVIEW}\n"
        f"  section 7.5  {RUN_NAME_BASE}\n"
        f"The preview follows the switches; this follows the built modules.")

print(f"\nSCENARIO              {SCENARIO!r}"
      + (f"   OVERRIDDEN: {', '.join(sorted(_SCENARIO_DEPARTURES))}"
         if _SCENARIO_DEPARTURES else ""))
for _label, _value, _token in SCENARIO_SWITCHES:
    print(f"  {_label:<24} {str(_value):<22} "
          f"{'-> ' + _token if _token else ''}")
print(f"\nRUN_NAME_BASE         {RUN_NAME_BASE!r}")
print(f"                      matches the section 3 preview (asserted)")
if DOUBLE_MANAGED_GIS:
    print(f"  both managers on    GIS {list(DOUBLE_MANAGED_GIS)} "
          f"-- see section 6")
# Keyed on the module flag, not the schedule: dune_migration_on is held False
# wherever beach_dune_manager runs, fill year or not.
if GROIN_ENABLED and BEACH_DUNE_MANAGEMENT_ON[
        HATTERAS_DOMAINS.gis_to_pad(GROIN_UPDRIFT_GIS)]:
    print(f"  groin + beach/dune  updrift GIS {GROIN_UPDRIFT_GIS} is inside "
          f"the beach_dune_manager footprint;")
    print(f"                      beach_dune_manager holds dune_migration_on "
          f"False there")


# =============================================================================
# 8. COASTSAT TARGET RATES -- LOESS WINDOWS
# =============================================================================
# The observational target. The "calibBE" source/sink preset was fit against
# the curve produced here, and section 12 draws the model against it.
#
# LOESS runs at TRANSECT resolution using along-coast distance as x, then
# averages to domain resolution. The widest window is the target, implicitly:
# rate_comparison picks max(window_domains). TARGET_WINDOW names that choice.
#
# skip_southern_domains = 10 suppresses LOESS across GIS 1-10. That is
# DISPLAY-ONLY -- LOESS still fits over all transects and only the result is
# truncated, so the southern transects still pull the values just north of the
# cut. Left as built, because the calibrated preset was fit against this curve.
#
# The Buxton groin sits at GIS 5.5, so both its flanking domains fall in the
# unsmoothed zone: the groin's observational target is a raw mean over a handful
# of transects, not a point on the reference curve.

# --- 8.1 load both periods ----------------------------------------------------
# Both periods always load. This is observational data, not simulation forcing,
# and section 12 needs the non-active period for reference styling -- so it is
# deliberately NOT folded into HATTERAS_PERIODS.

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
]

LOESS_CONFIG = LoessConfig(window_domains=(7, 10), skip_southern_domains=10)

# Named here rather than left implicit. rate_comparison resolves the reference
# window as max(window_domains); this makes that choice visible, and the
# assertion below catches a window list whose maximum is not what was intended.
TARGET_WINDOW = 10
if TARGET_WINDOW != max(LOESS_CONFIG.window_domains):
    raise ValueError(
        f"TARGET_WINDOW={TARGET_WINDOW} but rate_comparison will use "
        f"max(window_domains)={max(LOESS_CONFIG.window_domains)} as the "
        f"reference curve -- section 12 would compare against a different "
        f"curve than the one reported here.")

cs_series = build_coastsat_series(
    COASTSAT_DATASETS, active_period_start=START_YEAR,
    loess_config=LOESS_CONFIG, domains=HATTERAS_DOMAINS)

CS_ACTIVE = next((cs for cs in cs_series if cs["active"]), None)
if CS_ACTIVE is None:
    raise RuntimeError(
        f"no CoastSat dataset has period_start == START_YEAR ({START_YEAR}); "
        f"loaded: {[cs['period_start'] for cs in cs_series]}")


# --- 8.2 the target, as a table ----------------------------------------------

def build_target_table(cs, loess_config, geometry, window):
    """Per-domain target rate, labelled with where each value came from.

    The target is not one curve: GIS 1..skip_southern_domains are raw
    per-domain means (LOESS is suppressed there), and the rest is the
    LOESS-smoothed reference window. Every row records which.

    Args:
        cs: One entry from build_coastsat_series.
        loess_config: LoessConfig used to build it.
        geometry: DomainGeometry describing the real-domain span.
        window: Reference window width, in domains.

    Returns:
        A DataFrame with gis_domain, target_lrr_m_yr, source, n_transects.

    Raises:
        ValueError: If `window` was not computed for this dataset.
    """
    match = [w for w in cs["windows"] if w["window"] == window]
    if not match:
        raise ValueError(
            f"window {window} not in {[w['window'] for w in cs['windows']]}")
    smoothed = dict(zip(match[0]["gis_x"], match[0]["smoothed"]))

    skip = loess_config.skip_southern_domains
    raw_x, raw_y = compute_domain_means(
        cs["transect_domains"], cs["transect_rates"],
        geometry.first_gis_id, skip)
    raw = dict(zip(raw_x, raw_y))
    counts = pd.Series(cs["transect_domains"]).value_counts()

    rows = []
    for gis in range(geometry.first_gis_id, geometry.last_gis_id + 1):
        if gis <= skip:
            value, source = raw.get(gis, np.nan), f"raw mean (D1-{skip})"
        else:
            value, source = smoothed.get(gis, np.nan), f"LOESS {window}-dom"
        rows.append(dict(gis_domain=gis, target_lrr_m_yr=value, source=source,
                         n_transects=int(counts.get(gis, 0))))
    return pd.DataFrame(rows)


def groin_differential(cs, updrift_gis, downdrift_gis):
    """Observed updrift-minus-downdrift mean LRR for one CoastSat period.

    The observational check on the dipole section 7 imposes: a positive
    differential means the updrift domain is retreating more slowly (or
    advancing faster) than the downdrift one, which is the signature a
    functioning groin should leave.

    Args:
        cs: One entry from build_coastsat_series.
        updrift_gis: GIS domain on the updrift side.
        downdrift_gis: GIS domain on the downdrift side.

    Returns:
        A dict of label, updrift, downdrift, differential, and per-domain
        transect counts. Rates are NaN where a domain has no transects.
    """
    lo, hi = min(updrift_gis, downdrift_gis), max(updrift_gis, downdrift_gis)
    gis_x, means = compute_domain_means(
        cs["transect_domains"], cs["transect_rates"], lo, hi)
    by_domain = dict(zip(gis_x, means))
    counts = pd.Series(cs["transect_domains"]).value_counts()
    up = by_domain.get(updrift_gis, np.nan)
    down = by_domain.get(downdrift_gis, np.nan)
    return dict(
        label=cs["label"], period_start=cs["period_start"],
        active=cs["active"], updrift=up, downdrift=down,
        differential=up - down,
        n_updrift=int(counts.get(updrift_gis, 0)),
        n_downdrift=int(counts.get(downdrift_gis, 0)),
    )


COASTSAT_TARGET = build_target_table(
    CS_ACTIVE, LOESS_CONFIG, HATTERAS_DOMAINS, TARGET_WINDOW)


# --- 8.4 report ---------------------------------------------------------------

_n_loess = int((COASTSAT_TARGET["source"].str.startswith("LOESS")).sum())
_n_raw = len(COASTSAT_TARGET) - _n_loess
_missing = COASTSAT_TARGET["target_lrr_m_yr"].isna().sum()

print(f"\nactive dataset        {CS_ACTIVE['label']}")
print(f"target window         LOESS {TARGET_WINDOW}-domain "
      f"({TARGET_WINDOW * HATTERAS_DOMAINS.domain_spacing_m / 1000:.1f} km)"
      f"   -- rate_comparison uses max(window_domains), asserted above")
print(f"COASTSAT_TARGET       {len(COASTSAT_TARGET)} domains: "
      f"{_n_raw} raw mean (D1-{LOESS_CONFIG.skip_southern_domains}), "
      f"{_n_loess} LOESS"
      + (f", {_missing} with no transects" if _missing else ""))
print(f"  range               "
      f"{COASTSAT_TARGET['target_lrr_m_yr'].min():+.2f} to "
      f"{COASTSAT_TARGET['target_lrr_m_yr'].max():+.2f} m/yr")
print("  This is the curve the 'calibBE' source/sink preset was fit "
      "against (section 4.3).")

print(f"\nUNSMOOTHED ZONE       D1-{LOESS_CONFIG.skip_southern_domains}, "
      f"raw per-domain means -- no LOESS line here")
for _row in COASTSAT_TARGET[
        COASTSAT_TARGET["gis_domain"]
        <= LOESS_CONFIG.skip_southern_domains].itertuples():
    _mark = ""
    if _row.gis_domain == GROIN_UPDRIFT_GIS:
        _mark = "  <- groin updrift"
    elif _row.gis_domain == GROIN_DOWNDRIFT_GIS:
        _mark = "  <- groin downdrift"
    print(f"    D{_row.gis_domain:<3} {_row.target_lrr_m_yr:+7.2f} m/yr   "
          f"n={_row.n_transects:<3}{_mark}")

print(f"\nGROIN DIFFERENTIAL    observed updrift D{GROIN_UPDRIFT_GIS} minus "
      f"downdrift D{GROIN_DOWNDRIFT_GIS}")
print("  positive = updrift doing better than downdrift = groin-like")
for _cs in cs_series:
    _d = groin_differential(_cs, GROIN_UPDRIFT_GIS, GROIN_DOWNDRIFT_GIS)
    _verdict = ("supports the dipole direction" if _d["differential"] > 0
                else "OPPOSITE to the imposed dipole")
    print(f"  {_d['label']:<28} {'(active)' if _d['active'] else '        '} "
          f"up {_d['updrift']:+6.2f}  down {_d['downdrift']:+6.2f}  "
          f"diff {_d['differential']:+6.2f} m/yr")
    print(f"    {'':<28}          n = {_d['n_updrift']} / "
          f"{_d['n_downdrift']} transects   {_verdict}")


# =============================================================================
# 9. FIGURE CONFIGURATION
# =============================================================================
# Not the plotting functions -- the configuration they need, gathered before
# section 12 calls them.
#
# Every plotting entry point defaults to the PACKAGE defaults, and
# DEFAULT_ANNOTATIONS is an AnnotationConfig() with every field empty. Omit
# `annotations=HATTERAS_ANNOTATIONS` at any one call site and that figure comes
# out with no villages, no piers, no groin line and no shoal zones -- nothing
# raises. Hence the keyword bundles below, splatted into every call.

# --- 9.1 site config every figure needs --------------------------------------

RATE_FIG_KWARGS = dict(
    domains=HATTERAS_DOMAINS,
    annotations=HATTERAS_ANNOTATIONS,
    loess_config=LOESS_CONFIG,        # section 8's, not DEFAULT_LOESS
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


# --- 9.2 animation jobs -------------------------------------------------------
# range: "real" | "all" | "groin" | "groin_span" | (gis_lo, gis_hi)
# mode:  "position" | "displacement" | "difference"
# pad:   half-width in domains, read only by "groin" / "groin_span"
#
# "groin" fans out into one GIF per structure in annotations.groins, so new
# structures are picked up without editing this list.

GIF_JOBS = [
    dict(range="real", mode="displacement"),
    dict(range="real", mode="position"),
    dict(range="groin", mode="position", pad=9),
    dict(range="groin", mode="difference", pad=9),
]


# --- 9.3 baseline for difference jobs ----------------------------------------

def scenario_run_name(switches, stem, **overrides):
    """Run name for a variant of the current scenario.

    Rebuilds the name from the same switch tokens section 7.5 used, with named
    switches overridden. Built from the token list rather than by editing the
    name string, so flipping "groin" cannot accidentally match inside "nogroin".

    Args:
        switches: SCENARIO_SWITCHES, as (label, value, token) triples.
        stem: RUN_NAME_STEM, the period prefix.
        **overrides: label -> replacement token. A token of None or "" drops
            that switch from the name, matching how 7.5 omits default states.

    Returns:
        The run name for that variant.

    Raises:
        KeyError: If an override names a switch that does not exist -- a typo
            would otherwise silently produce the unmodified name.
    """
    labels = {label for label, _, _ in switches}
    unknown = set(overrides) - labels
    if unknown:
        raise KeyError(f"no such switch: {sorted(unknown)}; have {sorted(labels)}")
    tokens = [overrides.get(label, token) for label, _, token in switches]
    return f"{stem}_{'_'.join(t for t in tokens if t)}"


GIF_BASELINE_NAME = None
GIF_BASELINE_NPY = None
if GROIN_ENABLED:
    GIF_BASELINE_NAME = scenario_run_name(
        SCENARIO_SWITCHES, RUN_NAME_STEM, groin="nogroin")
    _baseline = (OUTPUT_BASE_DIR / GIF_BASELINE_NAME
                 / f"{GIF_BASELINE_NAME}_shoreline_matrix.npy")
    GIF_BASELINE_NPY = str(_baseline) if _baseline.exists() else None


# --- 9.4 validation target -- the surveyed island position in the end year ----
# Read from the RAW transect CSVs, not the padded offset files: the padded files
# are each zeroed on their own most-seaward domain, so differencing two years
# subtracts a constant and flips the sign of the mean, which would put the
# target on the wrong side of the model. The raw CSVs share a fixed offshore
# datum, so differencing them is a real shoreline change.
RAW_OFFSET_DIR = HATTERAS_DATA_BASE / "2-brie-offset" / "raw_offsets"

# Column names in the raw dune-line intersection CSVs, matching COL_MAP in
# scripts/input_prep/2-brie-offset/island_offset_hybrid.py:50.
RAW_OFFSET_COLUMNS = dict(domain="domain_id", distance="ORIG_LEN", transect="LineID")


def load_absolute_dune_distance(year, geometry, raw_dir=RAW_OFFSET_DIR,
                                columns=RAW_OFFSET_COLUMNS):
    """Mean distance from the offshore datum to the dune line, per GIS domain.

    The absolute quantity the padded offset files are built from, BEFORE
    island_offset_hybrid.py subtracts that year's own minimum. Absolute
    distances share a fixed datum across years, so differencing two of them is
    a real shoreline change; differencing the padded files is not.

    One row per transect (`LineID`) is kept before averaging, mirroring
    island_offset_hybrid.py:108-113, so a transect sampled at many points does
    not outweigh one sampled at few.

    Args:
        year: Survey year; reads <raw_dir>/<year>_duneline_offset_raw.csv.
        geometry: DomainGeometry supplying the real-domain GIS range.
        raw_dir: Directory holding the raw per-year CSVs.
        columns: Mapping with "domain", "distance" and "transect" keys.

    Returns:
        1-D array of length geometry.num_real_domains, in meters, increasing
        LANDWARD, indexed by GIS domain. NaN where a domain has no transects.

    Raises:
        FileNotFoundError: If that year has no raw file.
        KeyError: If an expected column is missing.
    """
    path = raw_dir / f"{year}_duneline_offset_raw.csv"
    raw = pd.read_csv(path)
    missing = [c for c in columns.values() if c not in raw.columns]
    if missing:
        raise KeyError(f"{path.name}: missing column(s) {missing}; "
                       f"have {list(raw.columns)[:12]}...")

    per_transect = raw.drop_duplicates(
        subset=[columns["domain"], columns["transect"]])
    means = per_transect.groupby(columns["domain"])[columns["distance"]].mean()

    gis_ids = np.arange(geometry.first_gis_id, geometry.last_gis_id + 1)
    return means.reindex(gis_ids).to_numpy(dtype=float)


def build_shoreline_target(model_year0_m, start_year, end_year, geometry):
    """Surveyed end-year shoreline position, in the model's own x_s frame.

    Takes the model's year-0 position and adds the OBSERVED start->end change,
    so the two share a base and the gap between them is the model's misfit.
    Returns positions rather than a change because that is what
    make_shoreline_gif's `target_m` expects.

    Args:
        model_year0_m: shoreline_m[0], the run's year-0 position (raw x_s_TS
            convention, meters, increasing landward), one value per PADDED
            domain.
        start_year: Run start year; must have a raw offset file.
        end_year: Run end year; returns None if it has no raw offset file.
        geometry: DomainGeometry.

    Returns:
        (target_m, observed_change_m):
            target_m: padded array of target positions, NaN in the buffers.
            observed_change_m: real-domain change, + = landward.
        Both None if the end year was never surveyed.
    """
    if not (RAW_OFFSET_DIR / f"{end_year}_duneline_offset_raw.csv").exists():
        return None, None

    start_dist = load_absolute_dune_distance(start_year, geometry)
    end_dist = load_absolute_dune_distance(end_year, geometry)
    observed_change_m = end_dist - start_dist   # + = landward, as x_s_TS is

    target_m = np.full(geometry.total_domains, np.nan)
    real = slice(geometry.start_real_index, geometry.end_real_index)
    target_m[real] = np.asarray(model_year0_m)[real] + observed_change_m
    return target_m, observed_change_m


# --- 9.5 report, and the annotation guard ------------------------------------
# A swapped-in empty AnnotationConfig would strip every figure's geography
# silently; fail here instead, where the cause is one line away.

_ann_populated = any([HATTERAS_ANNOTATIONS.town_spans,
                      HATTERAS_ANNOTATIONS.village_lines,
                      HATTERAS_ANNOTATIONS.piers,
                      HATTERAS_ANNOTATIONS.groins,
                      HATTERAS_ANNOTATIONS.shoal_zones])
if not _ann_populated:
    raise RuntimeError(
        "HATTERAS_ANNOTATIONS is empty -- every figure would render with no "
        "geographic layer and no error. Check the import in section 1.")

print(f"\nannotations           "
      f"{len(HATTERAS_ANNOTATIONS.town_spans)} towns, "
      f"{len(HATTERAS_ANNOTATIONS.piers)} piers, "
      f"{len(HATTERAS_ANNOTATIONS.groins)} groin(s), "
      f"{len(HATTERAS_ANNOTATIONS.shoal_zones)} shoal zones")
print(f"loess_config          windows {LOESS_CONFIG.window_domains}, "
      f"skip D1-{LOESS_CONFIG.skip_southern_domains}   (section 8's)")
print(f"gif_config            {gif_config.fps} fps, stride "
      f"{gif_config.year_stride}, ocean at "
      f"{'bottom' if gif_config.ocean_at_bottom else 'top'}, "
      f"save_matrix={gif_config.save_matrix}")
print(f"sign convention       flip_sign_model={FLIP_SIGN_MODEL} "
      f"(up = seaward), real domains only={PLOT_REAL_DOMAINS_ONLY}")

print(f"\nGIF_JOBS              {len(GIF_JOBS)} job(s)")
_diff_jobs = [j for j in GIF_JOBS if j.get("mode") == "difference"]
for _job in GIF_JOBS:
    print(f"  {_job['range']:<11} {_job['mode']:<13}"
          + (f" pad={_job['pad']}" if "pad" in _job else ""))

print("\nbaseline              ", end="")
if not GROIN_ENABLED:
    print("n/a - this run has the groin off, so it IS the baseline")
    if _diff_jobs:
        print(f"  {len(_diff_jobs)} difference job(s) will be skipped by "
              f"make_all_shoreline_gifs")
elif GIF_BASELINE_NPY:
    print(f"found\n  {GIF_BASELINE_NAME}\n  {GIF_BASELINE_NPY}")
else:
    print(f"NOT FOUND\n  looked for {GIF_BASELINE_NAME}")
    print(f"  under {OUTPUT_BASE_DIR}")
    if _diff_jobs:
        print(f"  {len(_diff_jobs)} difference job(s) will be SKIPPED. "
              f"Run once with GROIN_ENABLED=False to create it.")


# =============================================================================
# 10. build_cascade + run_cascade_simulation
# =============================================================================
# The split is what lets section 11 hold a built-but-unstepped Cascade. BRIE's
# diffusivity and the groin's fillet prediction are only meaningful as initial
# conditions, and a prediction printed after the run is not one.
#
# `run_years` is TRANSITIONS, not states. Barrier3D seeds _x_s_TS = [x_s] at
# init and appends one entry per update, so N updates produce N+1 annual states
# spanning N years. The original signature took `nt` and looped `range(nt - 1)`
# with `nt = END_YEAR - START_YEAR`, which ran 19 updates for a 20-year period
# while dividing by 20 -- every rate came out low by 19/20, and storm years 20
# and 21 were never applied. Here time_step_count = run_years + 1 and the loop
# runs exactly run_years updates.
#
# Nourishment goes to cascade.nourishment_volume, which is where CASCADE reads
# it. Writing it onto the manager instead hits the attribute CASCADE overwrites
# one line before the manager reads it, so the fill spends the init default.


def build_cascade(
    run_years, name, storm_file, alongshore_section_count, num_cores,
    rmin, rmax, elevation_file, dune_file,
    dune_design_elevation, dune_minimum_elevation,
    road_ele, road_width, road_setback,
    overwash_filter, overwash_to_dune,
    nourishment_volume, background_erosion,
    roadway_management_on, beach_dune_manager_on,
    sea_level_rise_rate, sea_level_constant,
    sandbag_management_on, sandbag_elevation,
    enable_shoreline_offset, shoreline_offset,
    wave_height, wave_period, wave_asymmetry, wave_angle_high_fraction,
    berm_elevation, MHW, groin_callback=None,
):
    """Constructs a Cascade and attaches the groin, without stepping it.

    Separate from run_cascade_simulation so section 11 can inspect the model
    before any time step runs: BRIE's diffusivity and the shoreface depth are
    only meaningful as initial conditions, and the fillet prediction they feed
    is only a prediction if it is printed before the run.

    Args:
        run_years: Annual transitions the run will simulate. The model is built
            with time_step_count = run_years + 1, giving run_years + 1 annual
            states. See the section 10 comment on why this is not run_years.
        name: Run name, passed to Cascade.
        storm_file, elevation_file, dune_file: Barrier3D input paths.
        alongshore_section_count: Padded domain count.
        num_cores: Cores for the parallel Barrier3D step.
        rmin, rmax: Dune growth rate bounds, per domain.
        dune_design_elevation: Rebuild target, m MHW. roadway_manager raises
            this to berm + 1.0 m on the first step if it is lower.
        dune_minimum_elevation: Rebuild trigger, m MHW. roadway_manager raises
            this to berm + 0.3 m on the first step if it is lower.
        road_ele, road_width, road_setback: Padded roadway forcing.
        overwash_filter, overwash_to_dune: beach_dune_manager forcing.
        nourishment_volume: Per-domain init volume. Every scheduled year
            overwrites it via the schedule; see section 6.
        background_erosion: Padded source/sink rates.
        roadway_management_on, beach_dune_manager_on: Per-domain module flags.
        sea_level_rise_rate, sea_level_constant: RSLR forcing.
        sandbag_management_on, sandbag_elevation: Sandbag forcing.
        enable_shoreline_offset, shoreline_offset: Island orientation, in dam.
        wave_height, wave_period, wave_asymmetry, wave_angle_high_fraction:
            Wave climate. wave_height also sets BRIE's shoreface depth,
            d_sf = 8.9 * Hs (brie.py:270).
        berm_elevation, MHW: Barrier3D datums.
        groin_callback: A GroinCallback to attach, or None.

    Returns:
        The constructed Cascade, before any update().
    """
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

        # run_years transitions need run_years + 1 states. TMAX is set from
        # this (brie_coupler.py:117), and the loop runs exactly run_years
        # updates, so the last write lands on the final valid index.
        time_step_count=run_years + 1,

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

    if groin_callback is not None:
        cascade._groin_callback = groin_callback

    return cascade


def run_cascade_simulation(
    cascade, run_years, name, run_dir, start_year, geometry,
    alongshore_section_count,
    historical_road_events=(), relocations_enabled=True, setback_check=None,
    nourishment_schedule=None, groin_callback=None, progress=None,
):
    """Steps a built Cascade through its period and writes the run artifacts.

    Args:
        cascade: A Cascade from build_cascade, not yet stepped.
        run_years: Annual transitions to simulate; the loop runs exactly this
            many updates.
        name: Run name, used for output filenames.
        run_dir: Directory for the saved model and logs.
        start_year: Calendar year of the run's first state.
        geometry: DomainGeometry, for GIS <-> pad translation.
        alongshore_section_count: Padded domain count.
        historical_road_events: RelocationEvent / BridgeEvent sequence.
        relocations_enabled: Global toggle for relocation events.
        setback_check: {gis: measured_setback_m} reported beside relocations.
        nourishment_schedule: A NourishmentSchedule, or None for no fills.
        groin_callback: The attached GroinCallback, for diagnostics output.
        progress: Optional tqdm-like object with update() and write().
            When given, the year counter advances the bar and every event
            message is written through it, so the messages scroll above a
            live bar instead of shredding it. None falls back to a plain
            carriage-return counter, so callers without tqdm still work.

    Returns:
        The same Cascade, after the run.
    """
    nourishment_log = []
    # One emitter for everything the loop says, so the display mechanism is
    # the caller's choice and this function does not care which it is.
    emit = progress.write if progress is not None else print

    for time_step in range(run_years):
        current_year = start_year + time_step

        # --- historical beach nourishment -----------------------------------
        # apply_to_cascade rewrites BOTH nourish_now and nourishment_volume in
        # full every year, so nothing carries over. It writes the volume to
        # cascade.nourishment_volume, which stock Cascade.update() copies into
        # each BeachDuneManager before calling it -- see section 6.
        if nourishment_schedule is not None:
            applied = nourishment_schedule.apply_to_cascade(
                cascade, current_year)
            if applied:
                print(f"\n  -> nourishment {current_year}:")
                for _row in applied:
                    print(f"       GIS {_row['gis']:>3} (pad {_row['pad']:>3})  "
                          f"{_row['volume_m3_per_m']:.1f} m^3/m")
                    nourishment_log.append(dict(
                        run_name=name, time_step=time_step, **_row))
        else:
            cascade.nourish_now = np.zeros(alongshore_section_count)

        # --- historical roadway events --------------------------------------
        for _event in historical_road_events or ():
            if current_year != _event.year:
                continue
            rows = roadway_module.apply_historical_event(
                cascade, _event, geometry,
                relocations_enabled=relocations_enabled,
                setback_check=setback_check)
            if not rows:
                print(f"\n  -> {current_year} event skipped: {_event.note}")
                continue
            print(f"\n  -> {current_year}: {_event.note}")
            for _row in rows:
                if _row["kind"] == "bridge":
                    emit(f"       GIS {_row['gis']:>3} (pad {_row['pad']:>3})  "
                         f"roadway management OFF")
                    continue
                _chk = (f" | 2004 measured {_row['check_m']:.0f} m"
                        if _row["check_m"] is not None else "")
                emit(f"       GIS {_row['gis']:>3} (pad {_row['pad']:>3})  "
                     f"setback {_row['old_setback_m']:.1f} + "
                     f"{_row['displacement_m']:.0f} -> "
                     f"{_row['new_setback_m']:.1f} m{_chk}")
                for _warn in _row["warnings"]:
                    emit(f"                   [warn] {_warn}")
                if _row["warnings"]:
                    emit(f"                   -> applied anyway "
                         f"(prescribed historical event)")

        cascade.update()

        if progress is not None:
            progress.update(1)
        else:
            print(f"\rYear {time_step + 1}/{run_years}", end="", flush=True)

        if getattr(cascade, "b3d_break", False):
            print(f"\nModel stopped at year {time_step + 1} (b3d_break)")
            break

    # --- did the run do what it was configured to do? ------------------------
    _states = len(cascade.barrier3d[0].x_s_TS)
    if _states != run_years + 1:
        print(f"\n  NOTE: {_states} annual states for {run_years} run_years "
              f"(expected {run_years + 1}) -- the run ended early")
    if groin_callback is not None and len(groin_callback.year_TS) == 0:
        print("\n" + "!" * 74)
        print("WARNING: the groin callback was never called. The pre-AST hook")
        print("in cascade_groin.py is missing, so this run is identical to a")
        print("no-groin run despite GROIN_ENABLED being True.")
        print("!" * 74)

    # --- artifacts -----------------------------------------------------------
    os.makedirs(run_dir, exist_ok=True)
    # The .npz model pickle is ~160 MB and is what lets a figure be re-derived
    # without re-running, so it is written by default. Everything downstream
    # (rate CSV, shoreline matrix, metadata) is written regardless, so a run
    # skipped here is still a complete row in run_index.csv -- just one that
    # cannot be re-plotted from state.
    if RUN_CONFIG.save_model_state:
        cascade.save(run_dir)
        print(f"\n  saved: {run_dir}")
    else:
        print(f"\n  model state NOT saved (HAT_SAVE_MODEL_STATE=false)")

    if nourishment_log:
        _bn_csv = os.path.join(run_dir, f"{name}_nourishment_log.csv")
        pd.DataFrame(nourishment_log).to_csv(_bn_csv, index=False)
        print(f"  nourishment log ({len(nourishment_log)} events): {_bn_csv}")

    if groin_callback is not None and groin_callback.year_TS:
        _groin_csv = os.path.join(run_dir, f"{name}_groin_diagnostics.csv")
        pd.DataFrame(groin_callback.diagnostics_frame()).to_csv(
            _groin_csv, index=False)
        print(f"  groin diagnostics ({len(groin_callback.year_TS)} yrs): "
              f"{_groin_csv}")

    return cascade


print("\nbuild_cascade + run_cascade_simulation defined")
print(f"  run_years -> {RUN_YEARS} transitions, "
      f"time_step_count={RUN_YEARS + 1}, {RUN_YEARS + 1} annual states "
      f"({START_YEAR}-{END_YEAR})")
print("  nourishment via BN_SCHEDULE.apply_to_cascade -> "
      "cascade.nourishment_volume")
print("  road events via cascade_pipeline.roadway.apply_historical_event")


# =============================================================================
# 11. INITIALIZE CASCADE -- SINGLE CONFIG, NO SWEEP
# =============================================================================

# --- 11.1 parameters sections 2-9 do not produce -----------------------------

NUM_CORES = 1        # >1 has crashed on this configuration; leave at 1

# --- dune growth (Barrier3D logistic growth bounds) --------------------------
RMIN = [0.55] * HATTERAS_DOMAINS.total_domains
RMAX = [0.95] * HATTERAS_DOMAINS.total_domains

# --- dune rebuild thresholds, m MHW ------------------------------------------
# Both are floored by roadway_manager on the first step:
#   dune_design_elevation  = max(passed, BermEl * 10 + 1.0)
#   dune_minimum_elevation = max(passed, BermEl * 10 + 0.3)
# Stated in the documented unit rather than the original's 0.01 "# dam", which
# was the wrong unit for the argument and far below the floor anyway.
DUNE_DESIGN_ELEVATION_M = 3.0    # rebuild target
DUNE_MINIMUM_ELEVATION_M = 0.0   # rebuild trigger: let CASCADE's berm floor
                                 # govern, explicitly rather than by accident
DUNE_DESIGN_ELEVATION = [DUNE_DESIGN_ELEVATION_M] * HATTERAS_DOMAINS.total_domains
DUNE_MINIMUM_ELEVATION = [DUNE_MINIMUM_ELEVATION_M] * HATTERAS_DOMAINS.total_domains

# --- roadway ------------------------------------------------------------------
ROAD_ELEVATION = road_elevation_full   # per-domain, m MHW, from section 5
ROAD_WIDTH = 20.0

# --- wave climate: one Hs, no sweep ------------------------------------------
# A sweep over Hs is a separate script; this runs one configuration.
Hs = 2.5                              # m, calibration value
FIXED_WAVE_PERIOD = 8                 # s
FIXED_WAVE_ASYMMETRY = 0.7
FIXED_WAVE_ANGLE_HIGH_FRACTION = 0.1

# --- datums -------------------------------------------------------------------
BERM_ELEVATION = 1.7    # m NAVD88, Hatteras Island, NCDOT-derived via NC State
MHW_ELEVATION = 0.36    # m NAVD88, Duck NC gauge (NOAA 8651370)

# --- sandbags: off for the hindcast ------------------------------------------
ENABLE_SANDBAG_PLACEMENT = False
SANDBAG_MANAGEMENT_ON = [ENABLE_SANDBAG_PLACEMENT] * HATTERAS_DOMAINS.total_domains
SANDBAG_ELEVATION = 0

SEA_LEVEL_CONSTANT = True

# --- run identity, from section 7.5's derived name ----------------------------
# OVERWRITE is checked here, before the model is built, so a name collision
# costs nothing either way.
#
#   False  the matrix default. A directory that already holds a result
#          stops the run. 12.3 resolves the groin baseline by directory
#          name, so replacing one run's output silently redefines the
#          paired run's answer -- which is why this refuses rather than
#          asks.
#   True   iterating on one scenario: tweak a value, re-run, read the
#          figures, tweak again. The directory is EMPTIED and reused, so
#          it never mixes two trials' files, and run_index.csv's row for
#          this name is replaced. The previous trial is gone. Set it back
#          to False before running the matrix.
OVERWRITE = RUN_CONFIG.overwrite

RUN_NAME = RUN_NAME_BASE
RUN_DIR = str(OUTPUT_BASE_DIR / RUN_NAME)
# Read before the guard runs: with OVERWRITE=True the guard empties the
# directory, and a silent wipe is how a trial gets mistaken for the run it
# replaced.
_replacing = run_dir_contents(RUN_DIR)
guard_run_dir(RUN_DIR, overwrite=OVERWRITE)
print(f"\nRUN_DIR               {RUN_DIR}"
      + ("   (OVERWRITE=True)" if OVERWRITE else ""))
if _replacing:
    print(f"  replaced            {len(_replacing)} file(s) from the previous run")


# --- 11.2 build ---------------------------------------------------------------

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
    shoreline_offset=island_offset_dam,
    wave_height=Hs,
    wave_period=FIXED_WAVE_PERIOD,
    wave_asymmetry=FIXED_WAVE_ASYMMETRY,
    wave_angle_high_fraction=FIXED_WAVE_ANGLE_HIGH_FRACTION,
    berm_elevation=BERM_ELEVATION,
    MHW=MHW_ELEVATION,
    groin_callback=GROIN_CALLBACK,
)


# --- 11.3 pre-run diagnostics, and the groin prediction ----------------------
# Reads the constructed model rather than the yaml, which resolves the active
# profile height section 7 had to leave open. The predicted fillet amplitude and
# extent are written down BEFORE the run; section 12 checks the emergent extent
# against them. Amplitude was tuned, extent was not.

def brie_r_ipl(cascade, theta_deg=0.0):
    """BRIE's diffusion number at the freshly built model's initial state.

    Reproduces `brie.py:1294`, which computes r_ipl as a local variable inside
    update() and never stores it:

        r_ipl = coast_diff[clip(round(90 - theta))] * dt / 2 / dy**2

    `_coast_diff` is the wave-climate-averaged shoreline diffusivity, built once
    in BRIE's __init__ and deleted by its finalize(). Nothing in CASCADE calls
    finalize, but the angle-dependent index changes as the shoreline evolves, so
    this is only the initial-condition value.

    Args:
        cascade: A built Cascade.
        theta_deg: Shoreline angle to evaluate at. 0 is shore-normal, the
            reference used for the fillet scaling.

    Returns:
        The dimensionless diffusion number.
    """
    brie = cascade._brie_coupler._brie
    index = int(np.clip(round(90 - theta_deg), 1, brie._wave_climl))
    return float(brie._coast_diff[index] * brie._dt / 2.0 / brie._dy ** 2)


R_IPL = brie_r_ipl(cascade)
_brie = cascade._brie_coupler._brie
_d_sf_m = float(_brie.d_sf)
_h_b_m = float(cascade.barrier3d[0].h_b_TS[0]) * DAM_TO_M
_profile_height_m = _d_sf_m + _h_b_m
_berm_floor_m = float(cascade.barrier3d[0].BermEl) * DAM_TO_M

print(f"\nbuilt                 {RUN_NAME}")
print(f"  run_dir             {RUN_DIR}")
print(f"  {RUN_YEARS} transitions -> {RUN_YEARS + 1} annual states "
      f"({START_YEAR}-{END_YEAR})")
print(f"  domains             {HATTERAS_DOMAINS.total_domains} padded, "
      f"{sum(ROADWAY_MANAGEMENT_ON)} roadway, "
      f"{sum(BEACH_DUNE_MANAGEMENT_ON)} beach/dune")

print(f"\nSHOREFACE             Hs {Hs} m -> d_sf = 8.9*Hs = {_d_sf_m:.2f} m "
      f"(brie.py:270)")
print(f"  h_b at t=0          {_h_b_m:.2f} m")
print(f"  active profile      {_profile_height_m:.2f} m  "
      f"(resolves the section 7 volume diagnostic)")

print(f"\nDUNE THRESHOLDS       passed -> floor roadway_manager will apply")
print(f"  design              {DUNE_DESIGN_ELEVATION_M:.2f} -> "
      f"{max(DUNE_DESIGN_ELEVATION_M, _berm_floor_m + 1.0):.2f} m MHW")
print(f"  minimum             {DUNE_MINIMUM_ELEVATION_M:.2f} -> "
      f"{max(DUNE_MINIMUM_ELEVATION_M, _berm_floor_m + 0.3):.2f} m MHW")
print(f"  BermEl * 10         {_berm_floor_m:.2f} m   "
      + ("(consistent with the metre reading)" if _berm_floor_m < 5
         else "SUSPECT -- see the notebook markdown on BermEl units"))

if GROIN_CALLBACK is not None:
    (GROIN_PREDICTED_AMPLITUDE_M,
     GROIN_PREDICTED_EXTENT_DOMAINS,
     GROIN_PREDICTED_EXTENT_M) = predict_fillet(
        trapping_rate_m_yr=GROIN_CALLBACK.M,
        r_ipl=R_IPL,
        run_years=RUN_YEARS,
        dy_m=HATTERAS_DOMAINS.domain_spacing_m,
    )
    _vol = implied_interception_m3_yr(
        GROIN_CALLBACK.M, _profile_height_m, HATTERAS_DOMAINS)
    _pct = 100.0 * _vol / REACH_TRANSPORT_LOSS_M3_YR

    print(f"\nGROIN                 attached, M = {GROIN_CALLBACK.M:.0f} m/yr, "
          f"updrift pad {GROIN_CALLBACK.updrift_pad} / "
          f"downdrift pad {GROIN_CALLBACK.downdrift_pad}")
    print(f"  r_ipl (t=0)         {R_IPL:.4f}  (shore-normal, brie.py:1294)")
    print(f"  PREDICTED amplitude {GROIN_PREDICTED_AMPLITUDE_M:.1f} m       = M / (4 * r_ipl)")
    print(f"  PREDICTED extent    {GROIN_PREDICTED_EXTENT_DOMAINS:.1f} domains "
          f"({GROIN_PREDICTED_EXTENT_M:.0f} m)   = dy * sqrt(2 * r_ipl * t)")
    print(f"  Extent does not depend on M. Section 12 checks the emergent")
    print(f"  extent against this number, written down before the run.")
    print(f"  implied transfer    {_vol:,.0f} m3/yr = {_pct:.0f}% of the "
          f"{REACH_TRANSPORT_LOSS_M3_YR:,.0f} m3/yr reach budget")
    if _pct > 50:
        print(f"                      BREACH, reported not corrected "
              f"(section 7)")
else:
    print(f"\nGROIN                 not attached (GROIN_ENABLED = "
          f"{GROIN_ENABLED})")
    print(f"  r_ipl (t=0)         {R_IPL:.4f}")


# =============================================================================
# 12. RUN THE LOOP, VERIFY, THEN FIGURES
# =============================================================================
# Verification comes before figures, deliberately, and the nourishment check
# ASSERTS. A fill that never reached the model invalidates the run, and that
# exact failure went unnoticed for a long time precisely because it produced
# plausible-looking output. Nothing is lost when it fires --
# run_cascade_simulation saves the model and its logs before returning.

GROIN_EXTENT_THRESHOLD_FRAC = 0.10   # fraction of peak effect defining "extent"

# --- 12.1 run ----------------------------------------------------------------
# Barrier3D seeds x_s_TS with one entry and appends one per update, so a model
# that has already been stepped has more than one. Stepping it again would run
# past TMAX and raise from deep inside Barrier3D.

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
)
if tqdm is not None:
    with tqdm(total=RUN_YEARS, desc=f"{RUN_NAME}", unit="yr") as _bar:
        cascade = run_cascade_simulation(progress=_bar, **_run_kwargs)
else:
    cascade = run_cascade_simulation(progress=None, **_run_kwargs)
RUN_SECONDS = time.perf_counter() - _t0
print(f"\nruntime               {RUN_SECONDS / 60:.1f} min "
      f"({RUN_SECONDS / RUN_YEARS:.1f} s per model year)")


# --- 12.2 verify --------------------------------------------------------------

shoreline_m = build_shoreline_matrix(cascade)
_states, _ = shoreline_m.shape

print(f"\nRUN LENGTH            {_states} annual states for {RUN_YEARS} "
      f"transitions")
if _states != RUN_YEARS + 1:
    print(f"  INCOMPLETE          expected {RUN_YEARS + 1}; the run stopped "
          f"early (b3d_break or drowning)")

# The denominator that produced a 5% bias in every earlier run. Checked, not
# trusted -- see the section 10 comment.
assert _states - 1 == RUN_YEARS, (
    f"span mismatch: {_states} states implies {_states - 1} years elapsed, "
    f"but RUN_YEARS is {RUN_YEARS}")
change_rate = compute_change_rate(
    shoreline_m, span_years=RUN_YEARS, flip_sign=FLIP_SIGN_MODEL)
print(f"  denominator         {RUN_YEARS} yr == states - 1   (asserted)")

# --- nourishment: the model's own record, not the schedule's intent ---------
BN_REPORT = nourishment.verify_nourishment(
    cascade, BN_SCHEDULE_APPLIED, BEACH_DUNE_MANAGEMENT_ON)
print(f"\nNOURISHMENT           {'OK' if BN_REPORT['ok'] else 'FAILED'}   "
      f"(read from each manager's _nourishment_volume_TS)")
for _row in BN_REPORT["confirmed"]:
    print(f"  ok       GIS {_row['gis']:>3} {_row['year']}  "
          f"{_row['actual_m3_per_m']:>7.1f} m^3/m")
for _row in BN_REPORT["wrong_volume"]:
    print(f"  WRONG    GIS {_row['gis']:>3} {_row['year']}  expected "
          f"{_row['expected_m3_per_m']:.1f}, got "
          f"{_row['actual_m3_per_m']:.1f} m^3/m")
for _row in BN_REPORT["missing"]:
    print(f"  MISSING  GIS {_row['gis']:>3} {_row['year']}  expected "
          f"{_row['expected_m3_per_m']:.1f} m^3/m, never applied")
for _row in BN_REPORT["unexpected"]:
    print(f"  EXTRA    GIS {_row['gis']:>3} {_row['year']}  "
          f"{_row['volume_m3_per_m']:.1f} m^3/m, not scheduled")
for _note in BN_REPORT["notes"]:
    print(f"  note: {_note}")
assert BN_REPORT["ok"], "nourishment did not reach the model as scheduled"

# --- the double-management consequence section 6 predicted ------------------
if DOUBLE_MANAGED_GIS:
    print(f"\nFROZEN SETBACKS       predicted in section 6 for GIS "
          f"{list(DOUBLE_MANAGED_GIS)}")
    for _row in nourishment.verify_setbacks_frozen(
            cascade, DOUBLE_MANAGED_GIS, HATTERAS_DOMAINS):
        print(f"  GIS {_row['gis']:>3}  {_row['setback_start_m']:.0f} -> "
              f"{_row['setback_end_m']:.0f} m"
              f"{'   FROZEN as predicted' if _row['frozen'] else ''}")

# --- which road_offset survived ----------------------------------------------------
ROAD_SUMMARY = roadway.summarise_road_management(
    cascade, HATTERAS_DOMAINS, HATTERAS_FIRST_ROAD_DOMAIN,
    HATTERAS_LAST_ROAD_DOMAIN)
_drowned = [r for r in ROAD_SUMMARY if r["drowned"]]
_blocked = [r for r in ROAD_SUMMARY if r["relocation_blocked"]]
print(f"\nROADWAY               {len(ROAD_SUMMARY)} managed domains, "
      f"{len(_drowned)} drowned, {len(_blocked)} relocation-blocked")
for _row in _drowned[:10]:
    print(f"  drowned  GIS {_row['gis']:>3}  last managed "
          f"{_row['last_managed_year']}  {_row['reason']}")
if len(_drowned) > 10:
    print(f"  ... and {len(_drowned) - 10} more (full list in the CSV)")


# --- 12.3 groin: the pre-registered extent check -----------------------------

def measure_groin_extent(shoreline_m, baseline_m, geometry, updrift_gis,
                         downdrift_gis, threshold_frac):
    """Alongshore extent of the groin's effect, from a paired baseline run.

    The effect is this run's final shoreline minus the no-groin baseline's,
    sign-flipped so positive means seaward of the baseline. Extent is the
    contiguous run of domains, outward from the groin, where the effect holds
    at or above `threshold_frac` of its own peak magnitude.

    Args:
        shoreline_m: [time x domain] matrix from this run.
        baseline_m: [time x domain] matrix from the no-groin run.
        geometry: DomainGeometry describing the padded array.
        updrift_gis, downdrift_gis: The groin's flanking domains.
        threshold_frac: Fraction of the peak effect defining the edge.

    Returns:
        A dict with effect (padded array, m), peak_m, threshold_m, and the
        updrift/downdrift extents in domains and meters.
    """
    effect = -(np.asarray(shoreline_m)[-1] - np.asarray(baseline_m)[-1])
    peak = float(np.nanmax(np.abs(effect))) if effect.size else 0.0
    threshold = threshold_frac * peak

    def span(start_gis, step):
        # A zero peak means the two runs are identical, so the threshold is
        # also zero and every "abs(value) < 0" test is False -- the walk would
        # run to the end of the grid and report the whole island as fillet.
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
        effect=effect, peak_m=peak, threshold_m=threshold,
        updrift_domains=up, updrift_m=up * geometry.domain_spacing_m,
        downdrift_domains=down, downdrift_m=down * geometry.domain_spacing_m,
    )


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
    print(f"\nGROIN EXTENT          measured against {GIF_BASELINE_NAME}")
    print(f"  peak effect         {GROIN_EXTENT['peak_m']:.1f} m   "
          f"(threshold {GROIN_EXTENT_THRESHOLD_FRAC:.0%} = "
          f"{GROIN_EXTENT['threshold_m']:.1f} m)")
    print(f"  updrift             {GROIN_EXTENT['updrift_domains']} domains "
          f"({GROIN_EXTENT['updrift_m']:.0f} m)")
    print(f"  downdrift           {GROIN_EXTENT['downdrift_domains']} domains "
          f"({GROIN_EXTENT['downdrift_m']:.0f} m)")
    print(f"  PREDICTED (sec 11)  "
          f"{GROIN_PREDICTED_EXTENT_DOMAINS:.1f} domains "
          f"({GROIN_PREDICTED_EXTENT_M:.0f} m)")
    print(f"  Amplitude was tuned; extent was not. This is the part that "
          f"could have failed.")


# --- 12.4 write ---------------------------------------------------------------

_shoreface_depth_m = float(cascade._brie_coupler._brie.d_sf)

run = RunInfo(
    run_name=RUN_NAME, run_dir=RUN_DIR,
    start_year=START_YEAR, end_year=END_YEAR, Hs=Hs,
    flip_sign_model=FLIP_SIGN_MODEL,
    background_erosion_on=USE_BACKGROUND_EROSION,
)

_rate_csv = os.path.join(RUN_DIR, f"{RUN_NAME}_shoreline_change_rate.csv")
pd.DataFrame({
    "gis_domain": np.arange(HATTERAS_DOMAINS.first_gis_id,
                            HATTERAS_DOMAINS.last_gis_id + 1),
    "change_rate_m_yr": change_rate[HATTERAS_DOMAINS.start_real_index:
                                    HATTERAS_DOMAINS.end_real_index],
}).to_csv(_rate_csv, index=False)
print(f"\nwrote                 {os.path.basename(_rate_csv)}")

if ROAD_SUMMARY:
    _road_csv = os.path.join(RUN_DIR, "road_management_summary.csv")
    pd.DataFrame(ROAD_SUMMARY).to_csv(_road_csv, index=False)
    print(f"                      {os.path.basename(_road_csv)}")

# --- skill against the section 8 target --------------------------------------
# Both spans are reported. The end domains carry the locked source/sink values
# (tens of m/yr), so an island-wide RMSE for a calibBE or edgeBE run is
# dominated by two domains that were pinned rather than predicted -- and a
# zeroBE run has no such term. Ranking the presets on the island-wide number
# alone would mostly rank their boundary treatment.
SKILL = skill_vs_target(change_rate, COASTSAT_TARGET, HATTERAS_DOMAINS)
print(f"\nSKILL vs CoastSat     model - target, m/yr")
print(f"  island-wide         bias {SKILL['mean_bias_m_yr']:+.3f}   "
      f"RMSE {SKILL['rmse_m_yr']:.3f}   (n={SKILL['n_domains']})")
print(f"  interior (GIS 2-89) bias {SKILL['mean_bias_interior_m_yr']:+.3f}   "
      f"RMSE {SKILL['rmse_interior_m_yr']:.3f}   "
      f"(n={SKILL['n_domains_interior']})")

# --- run metadata: the scenario, plus what distinguishes this run -----------
# One structure renders both files: the .txt to read, the .json to parse.
# Built together so they cannot disagree -- the earlier version emitted only
# prose, so anything downstream had to re-parse it.
_GIT = git_provenance(PROJECT_BASE_DIR)
_TIMESTAMP = timestamp()
GENERATED_BY = "HAT_hindcast_1984_2024.py"

# --- what the source/sink preset actually was -------------------------
# The preset's name and its numbers are separate facts.
# HATTERAS_BE_RATES_EDGE is a slice of HATTERAS_BE_RATES_CALIBRATED at
# the end domains, so editing an end domain to test edgeBE changes
# calibBE too -- and both keep their names. Recording the values is what
# makes two trials of one preset tell themselves apart afterwards.
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
        # The three that drift silently across a multi-day batch: which
        # Cascade class ran, which extractor version built the init surface,
        # and which commit produced both.
        "use_sandbox_cascade": (USE_SANDBOX_CASCADE,
                                "cascade.cascade_groin, not cascade.cascade"),
        "topo_dune_version": TOPO_DUNE_VERSION,
        "parameter_file": PARAMETER_FILE,
        "git_commit": _GIT["commit"],
        "git_branch": _GIT["branch"],
        "git_dirty": (_GIT["dirty"],
                      "True: the commit alone does not reproduce this run"),
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
        "target": (f"CoastSat LOESS {TARGET_WINDOW}-domain",
                   "raw means over GIS 1-"
                   f"{LOESS_CONFIG.skip_southern_domains}"),
        "mean_bias_m_yr": f"{SKILL['mean_bias_m_yr']:+.4f}",
        "rmse_m_yr": f"{SKILL['rmse_m_yr']:.4f}",
        "mean_bias_interior_m_yr": (
            f"{SKILL['mean_bias_interior_m_yr']:+.4f}",
            "GIS 2-89: the locked end domains excluded"),
        "rmse_interior_m_yr": f"{SKILL['rmse_interior_m_yr']:.4f}",
    },
    "verification": {
        "nourishment_ok": BN_REPORT["ok"],
        "roads_drowned": len(_drowned),
        "roads_reloc_blocked": len(_blocked),
    },
}

if GROIN_CALLBACK is not None:
    _META["groin"].update({
        "trapping_rate_m_yr": GROIN_CALLBACK.M,
        "updrift / downdrift": f"GIS {GROIN_UPDRIFT_GIS} / {GROIN_DOWNDRIFT_GIS}",
        "install_year": GROIN_CALLBACK.install_year,
        "deterioration": f"{GROIN_CALLBACK.deterioration_mode}, "
                         f"floor {GROIN_CALLBACK.deterioration_fraction}",
        "r_ipl_t0": f"{R_IPL:.4f}",
        "predicted_extent_m": f"{GROIN_PREDICTED_EXTENT_M:.0f}",
    })
    if GROIN_EXTENT is not None:
        _META["groin"]["measured_extent_m"] = (
            f"{GROIN_EXTENT['updrift_m']:.0f} updrift / "
            f"{GROIN_EXTENT['downdrift_m']:.0f} downdrift")

_meta_txt, _meta_json = write_run_metadata(
    RUN_DIR, RUN_NAME, _META,
    header=[f"CASCADE run metadata -- generated by {GENERATED_BY}, "
            f"section 12.",
            "Companion .json holds the same values, machine-readable."])
print(f"                      {_meta_txt.name}")
print(f"                      {_meta_json.name}")

# --- cross-run index: one row per run, for comparing the matrix --------------
# Keyed on run_name, so a re-run replaces its row rather than adding a second.
# This is what makes 12 runs comparable without opening 12 metadata files.
# At OUTPUT_ROOT rather than the period directory: run_name already
# carries the period stem, so one file covers the whole matrix and the
# two periods stay comparable in a single table.
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
    "groin_trapping_m_yr": (GROIN_CALLBACK.M if GROIN_CALLBACK is not None
                            else np.nan),
    "roadway_management": ENABLE_ROADWAY_MANAGEMENT,
    "relocations_enabled": ENABLE_HISTORICAL_ROAD_RELOCATIONS,
    "beach_dune_management": ENABLE_BEACH_DUNE_MANAGEMENT,
    "nourishment_fills": ENABLE_NOURISHMENT_FILLS,
    "bdm_domains": int(sum(BEACH_DUNE_MANAGEMENT_ON)),
    "nourishment_projects": len(BN_SCHEDULE_APPLIED.projects),
    "Hs_m": Hs,
    "rslr_m_yr": SEA_LEVEL_RISE_RATE,
    "annual_states": _states,
    "run_complete": _states == RUN_YEARS + 1,
    "mean_bias_m_yr": SKILL["mean_bias_m_yr"],
    "rmse_m_yr": SKILL["rmse_m_yr"],
    "mean_bias_interior_m_yr": SKILL["mean_bias_interior_m_yr"],
    "rmse_interior_m_yr": SKILL["rmse_interior_m_yr"],
    "nourishment_ok": BN_REPORT["ok"],
    "roads_drowned": len(_drowned),
    "roads_reloc_blocked": len(_blocked),
    "groin_extent_updrift_m": (GROIN_EXTENT["updrift_m"]
                               if GROIN_EXTENT is not None else np.nan),
    "groin_extent_downdrift_m": (GROIN_EXTENT["downdrift_m"]
                                 if GROIN_EXTENT is not None else np.nan),
    "runtime_min": round(RUN_SECONDS / 60, 1),
    "use_sandbox_cascade": USE_SANDBOX_CASCADE,
    "topo_dune_version": TOPO_DUNE_VERSION,
    "git_commit": _GIT["commit"][:12],
    "git_dirty": _GIT["dirty"],
}
RUN_INDEX = append_run_index(RUN_INDEX_PATH, _index_row)
print(f"                      {RUN_INDEX_FILENAME}  "
      f"({len(RUN_INDEX)} runs indexed)")


# --- 12.5 figures -------------------------------------------------------------
# Site config arrives via section 9's bundles; omitting them would silently
# strip the geographic layer (DEFAULT_ANNOTATIONS is empty).

plot_rate_comparison(
    change_rate, cs_series, run,
    real_domains_only=PLOT_REAL_DOMAINS_ONLY,
    sea_level_rise_rate_m_yr=SEA_LEVEL_RISE_RATE,
    save_path=os.path.join(
        RUN_DIR, f"{RUN_NAME}_shoreline_change_rate"
        f"{'_REAL_DOMAINS_ONLY' if PLOT_REAL_DOMAINS_ONLY else ''}.png"),
    show=SHOW_FIGURES, **RATE_FIG_KWARGS)

plot_annotated_rate_comparison(
    change_rate, cs_series, run,
    sea_level_rise_rate_m_yr=SEA_LEVEL_RISE_RATE,
    save_path=os.path.join(RUN_DIR, f"{RUN_NAME}_annotated.png"),
    show=SHOW_FIGURES, **RATE_FIG_KWARGS)

# Section 9.4's validation target, resolved now that the run has a year 0.
# Buffers are NaN, so the comparison below is over the real domains only.
SHORELINE_TARGET_M, OBSERVED_CHANGE_M = build_shoreline_target(
    shoreline_m[0], START_YEAR, END_YEAR, HATTERAS_DOMAINS)

print()
if SHORELINE_TARGET_M is None:
    print(f"target                none -- no {END_YEAR} dune-line survey in "
          f"{RAW_OFFSET_DIR.name}/; GIFs draw the model alone")
else:
    _real = slice(HATTERAS_DOMAINS.start_real_index, HATTERAS_DOMAINS.end_real_index)
    _modeled = (shoreline_m[-1] - shoreline_m[0])[_real]
    _misfit = _modeled - OBSERVED_CHANGE_M
    print(f"target                observed {END_YEAR} dune line "
          f"({np.count_nonzero(~np.isnan(OBSERVED_CHANGE_M))} real domains)")
    print(f"  observed change     mean {np.nanmean(OBSERVED_CHANGE_M):+7.1f} m   "
          f"({np.nanmin(OBSERVED_CHANGE_M):+.1f} to "
          f"{np.nanmax(OBSERVED_CHANGE_M):+.1f})   + = landward")
    print(f"  modeled change      mean {_modeled.mean():+7.1f} m   "
          f"({_modeled.min():+.1f} to {_modeled.max():+.1f})")
    print(f"  misfit              mean {np.nanmean(_misfit):+7.1f} m   "
          f"RMSE {np.sqrt(np.nanmean(_misfit ** 2)):.1f} m")

GIF_PATHS = make_all_shoreline_gifs(
    shoreline_m, run, GIF_JOBS,
    baseline_npy=GIF_BASELINE_NPY,
    target_m=SHORELINE_TARGET_M, **GIF_KWARGS)

plt.close("all")
print(f"\ndone                  {RUN_DIR}")
