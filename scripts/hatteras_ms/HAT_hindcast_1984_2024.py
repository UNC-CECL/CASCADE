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

WHERE THE SWITCHES LIVE

  `hat_run.yaml`, beside this file. Edit it, save it, run. That file is the
  interface; the sections below only read it, through HAT_hindcast_config:

    section 1   output.show_figures -- read at import, because it selects a
                matplotlib backend
    section 3   start_year, source_sink, scenario, relocations, offset_mode
                (scenario expands to ENABLE_ROADWAY_MANAGEMENT,
                ENABLE_BEACH_DUNE_MANAGEMENT, ENABLE_NOURISHMENT_FILLS and
                ENABLE_HISTORICAL_ROAD_RELOCATIONS; commented override lines
                sit directly below the table)
    section 7   groin.enabled, groin.trapping_M, groin.deterioration_f
    section 9   output.make_gifs
    section 11  physics.wave_height_Hs, sandbags, output.overwrite,
                output.save_model_state

  An environment variable beats the file, so `HAT_run_all.py` can drive a
  matrix without editing anything; it also sets HAT_IGNORE_SETTINGS=1, so a
  half-finished experiment left in the yaml cannot reach a batch run.

  STILL TYPED IN THIS FILE, because they are properties of the study rather
  than of a run: dune thresholds and datums (section 11), FLIP_SIGN_MODEL and
  PLOT_REAL_DOMAINS_ONLY (section 9), the groin's geometry and deterioration
  schedule (section 7), NUM_CORES (section 11).

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

# --- the settings file, read before anything it selects ----------------------
# `hat_run.yaml` is where a run is chosen; see HAT_hindcast_config for the
# precedence rules. It is read HERE rather than in section 3 because two of
# its values are settled at import time and cannot be changed afterwards:
#
#   use_sandbox_cascade  decides which Cascade class cascade_pipeline.hindcast
#                        binds, which happens the moment it is imported below.
#                        It is pinned True and has no line in the yaml -- see
#                        the "not settable here" block at the foot of that
#                        file for why it must not follow groin.enabled.
#   show_figures         decides the matplotlib backend, which is fixed by the
#                        first figure created.
#
# Everything else is re-read in section 3, so an interactive edit to the yaml
# applies on a re-run of that cell without restarting the kernel. Section 3
# also checks these two against the file as it stands then, and raises if they
# have changed -- a stale sandbox flag is the failure where the groin silently
# does nothing.
from HAT_hindcast_config import load_run_config          # noqa: E402

_BOOT_CONFIG = load_run_config()
os.environ["CASCADE_USE_SANDBOX"] = (
    "1" if _BOOT_CONFIG.use_sandbox_cascade else "0")

import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt


from cascade.groin import GroinCallback, predict_fillet

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
)
from cascade_pipeline.coastsat_loess import (
    CoastSatDataset,
    LoessConfig,
    build_coastsat_series,
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
from cascade_pipeline.shoreline import (build_shoreline_matrix,
                                        compute_change_rate, compute_lrr)

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
#
# `output.show_figures: null` in hat_run.yaml means "whichever file is being
# run decides", which here is False. Setting it true in the yaml forces
# figures open in a headless run, which is why it is not the default.
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
print(f"HATTERAS_DOMAINS.total_domains = {HATTERAS_DOMAINS.total_domains}")


# =============================================================================
# 2. DUNE/TOPO -- PER PERIOD
# =============================================================================
# NO LONGER period-independent (changed 2026-08-25). Until then both periods
# read one topography and this section said so. They now start from different
# DEMs, so the product is selected from HATTERAS_PERIODS:
#
#     1984  ->  1984-start   from DEM 2009-2014-1996
#     2004  ->  2004-start   from DEM 2009-2014
#
# The VERSION within a product is still resolved, never pinned.

# --- 2.1 project paths -------------------------------------------------------

# There is no TOPO_DUNE_INIT_YEAR any more. The arrays carry no year - the
# period is the PRODUCT DIRECTORY - and this runner no longer builds their
# names at all; build_domain_file_paths() delegates to the resolver. See the
# note at the top of hat_topo_version.py for why a per-period tag was tried
# and reverted.

# WHICH EXTRACTION -- resolved, not pinned. topo_dirs() reads VERSION out of
# HAT_dune_topo_extractor.py, so the runner, the dune-start road setbacks and
# the audits always describe the same extraction, and a version missing from
# disk is a loud error rather than a silently stale read.
#
# Pinning it by hand is what let this runner sit on 2009_v3 while the setbacks
# were rebuilt on v4 (2026-08-19) and then v5 (2026-08-20). The setback is
# measured from interior row 0, so that combination places the road against a
# row that does not exist on the grid being run.
#
# To reproduce an older run deliberately:
#     topo_dirs("2004-start", override="v3").
from hat_topo_version import topo_dirs  # scripts/, on sys.path above
from hat_topo_version import BUFFER_DIR as _BUFFER_DIR

# _BOOT_CONFIG, not RUN_CONFIG: the period must be known HERE, and RUN_CONFIG is
# not loaded until section 3. The two are compared a few lines below section 3's
# reload, and section 3.0 re-asserts that this product matches the period that
# actually ran, so a boot/run divergence cannot silently pick the wrong barrier.
TOPO_PRODUCT = HATTERAS_PERIODS[_BOOT_CONFIG.start_year]["topo_product"]

_TOPO_DIR, _DUNE_DIR, TOPO_DUNE_VERSION = topo_dirs(TOPO_PRODUCT)

print(f"topography            {TOPO_PRODUCT} / {TOPO_DUNE_VERSION}  "
      f"(product from the period, version resolved)")

HATTERAS_DATA_BASE = PROJECT_BASE_DIR / "data" / "hatteras_init"
OUTPUT_ROOT = PROJECT_BASE_DIR / "output" / "raw_runs"
COASTSAT_BASE_DIR = (PROJECT_BASE_DIR / "scripts" / "input_prep"
                     / "5-scr" / "CoastSat")
PARAMETER_FILE = "Hatteras-CASCADE-parameters.yaml"  # resolved by CASCADE

BARRIER3D_DIR = HATTERAS_DATA_BASE / "1-barrier3d-domains"
# Taken from what topo_dirs() RETURNED rather than re-joined from parts. The
# old line rebuilt the path independently, which is how a resolver gets bypassed
# without anyone noticing - the same failure mode as HAT_road_elevation.py.
DUNE_TOPO_DIR = _TOPO_DIR.parent
BUFFER_DIR = _BUFFER_DIR

os.chdir(PROJECT_BASE_DIR)
# Runs are filed per period: section 3 builds OUTPUT_BASE_DIR =
# OUTPUT_ROOT / "<start>_<end>" once START_YEAR is expanded, so a
# 1984-2004 run and a 2004-2024 run of one scenario cannot land beside
# each other. run_index.csv stays at the root, covering both periods.
OUTPUT_ROOT.mkdir(parents=True, exist_ok=True)

reports.path_inventory([
    ("PROJECT_BASE_DIR", PROJECT_BASE_DIR),
    ("HATTERAS_DATA_BASE", HATTERAS_DATA_BASE),
    ("DUNE_TOPO_DIR", DUNE_TOPO_DIR),
    ("BUFFER_DIR", BUFFER_DIR),
    ("COASTSAT_BASE_DIR", COASTSAT_BASE_DIR),
    ("OUTPUT_ROOT", OUTPUT_ROOT),
])


# --- 2.2 build the padded file lists -----------------------------------------


# The PRODUCT is passed, not the directories: build_domain_file_paths now
# delegates to hat_topo_version.domain_arrays(), which resolves the directory
# and the filename together. Nothing here spells an array name.
ELEVATION_FILE_PATHS, DUNE_FILE_PATHS = build_domain_file_paths(
    HATTERAS_DOMAINS, TOPO_PRODUCT)

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


BARRIER3D_CONTRACT = load_barrier3d_contract(HATTERAS_DATA_BASE / PARAMETER_FILE)

print(f"\nContract from {PARAMETER_FILE}:")
print(f"  BarrierLength -> {BARRIER3D_CONTRACT['barrier_length_cells']} "
      f"alongshore cells")
print(f"  MHW           -> {BARRIER3D_CONTRACT['mhw_dam']:.3f} dam")
print(f"  BermEl        -> {BARRIER3D_CONTRACT['berm_el_dam']:.3f} dam "
      f"above MHW\n")

reports.run_units_check(ELEVATION_FILE_PATHS, DUNE_FILE_PATHS,
                        BARRIER3D_CONTRACT, HATTERAS_DOMAINS)


# =============================================================================
# 3. ISLAND ORIENTATION -- SET START_YEAR
# =============================================================================
# START_YEAR is the one flip. It selects a period from HATTERAS_PERIODS and
# everything in section 4 follows from it: run length, RSLR rate, storm series,
# background-erosion preset, road-setback file, nourishment settings.
#
# The run-selecting values in this section come from `hat_run.yaml`, read
# through HAT_hindcast_config: edit that file, save it, run. A driver can also
# set them through the environment without editing anything, and
# interactively you can still type over any assignment below -- the loaded
# value is only the default, and the assignment is still the last word.
#
# load_run_config() RE-READS the yaml on every call rather than reusing the
# object section 1 built. In a notebook the module is cached by the import
# system, so a fresh read is the only way an edit to the yaml applies without
# restarting the kernel.

from HAT_hindcast_config import (                     # noqa: E402
    load_run_config, describe as _describe_run_config, preflight as _preflight)

RUN_CONFIG = load_run_config()

# The two values section 1 already spent. They select an import and a
# matplotlib backend, so a yaml edited after section 1 ran would be reported
# in the block below while the process is still running the old choice --
# and for the sandbox flag that means a groin-on run whose groin silently
# does nothing. Raised, not warned.
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

# Read here rather than beside the offset load in section 4: the run-name
# preview below needs it, and a value the name depends on must be settled
# before the name is predicted.
OFFSET_MODE = RUN_CONFIG.offset_mode

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
# HAT_RELOCATIONS overrides the scenario when set; None leaves the preset
# in charge. Either way the departure is detected just below and the
# `reloc` token still comes from the switch, not from SCENARIO.
ENABLE_HISTORICAL_ROAD_RELOCATIONS = (
    _SCENARIO_PRESET["relocations"] if RUN_CONFIG.relocations is None
    else RUN_CONFIG.relocations)

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

# Section 2 picked the topography from _BOOT_CONFIG.start_year, before
# RUN_CONFIG existed. If those disagree the run would model the wrong barrier
# entirely, so it is checked rather than assumed.
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
# Filed by period, then by source/sink preset. The preset directory is
# redundant with the preset token in RUN_NAME on purpose: the token is what
# run_index.csv, the logs and the figure captions all key on, and the
# directory is only there so the three presets of one period can be read
# side by side instead of interleaved in one listing of thirty-odd runs.
OUTPUT_BASE_DIR = OUTPUT_ROOT / PERIOD_TAG / SOURCE_SINK_PRESET
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
    None if OFFSET_MODE == "asrun" else f"offset{OFFSET_MODE}",
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

# What this run will be called, where it will land, whether something already
# lives there, and roughly how long it will take -- reported HERE rather than
# at the section 11 guard so a name collision is visible before sections 4-10
# do their work. Advisory: `guard_run_dir` in 11 is still the authority on
# what a collision does, and this does not duplicate that decision.
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


# --- 3.1 island offsets ------------------------------------------------------
# One value per padded domain giving that domain's cross-shore starting
# position. Becomes `shoreline_offset` on the Cascade() call.


# OFFSET_MODE selects which shoreline_offset variant is built; see
# cascade_pipeline.hindcast.build_island_offset. The default, "asrun",
# reproduces the historical unit error (the file is METRES and Cascade
# wants metres, but load_island_offset_dam divided by 10), so previously
# published runs stay reproducible until the correction is adopted.
island_offset = build_island_offset(
    ISLAND_OFFSET_FILE, HATTERAS_DOMAINS, mode=OFFSET_MODE)
OFFSET_TILTS = island_offset_tilts(island_offset, HATTERAS_DOMAINS)

_real = slice(HATTERAS_DOMAINS.start_real_index, HATTERAS_DOMAINS.end_real_index)
print(f"\n{START_YEAR} offsets: {island_offset.size} padded domains | "
      f"real span {island_offset[_real].min() * DAM_TO_M:.0f}-"
      f"{island_offset[_real].max() * DAM_TO_M:.0f} m")


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


STORM_SERIES = load_storm_series(STORM_FILE)

reports.storm_report(storms=STORM_SERIES, storm_file=STORM_FILE,
                     run_years=RUN_YEARS)


# --- 4.3 source/sink (background erosion) ------------------------------------
# Per-domain rate in m/yr, passed to Barrier3D as `Rat`. Sign convention from
# cascade/brie_coupler.py: (-) = erosion, (+) = accretion. Presets are sparse --
# a domain absent from a preset gets 0.0 m/yr.


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

reports.background_erosion_report(
    preset=SOURCE_SINK_PRESET, start_year=START_YEAR,
    domain_rates=DOMAIN_BE_RATES, rates=BACKGROUND_EROSION_RATES,
    use_background_erosion=USE_BACKGROUND_EROSION, geometry=HATTERAS_DOMAINS,
    rates_2004_are_placeholder=HATTERAS_BE_RATES_2004_IS_PLACEHOLDER)


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


# --- 5.1 pre-flight audit: which road_offset will not survive year one -------------
# bulldoze tests the two rows FLANKING the road -- never the road's own cells --
# and drowns it when either flank is more than 20% water. A drowned road is not
# a warning: roadway_manager sets _drown_break and returns immediately on every
# later year, so the domain becomes an unmanaged barrier wearing a road label.

road_audit = roadway.audit_setbacks(
    ELEVATION_FILE_PATHS, road_setbacks_full, HATTERAS_DOMAINS, *_road_span,
    management_on=ROADWAY_MANAGEMENT_ON, config=ROADWAY)
audit_summary = roadway.summarise_audit(road_audit)

reports.road_audit_report(audit=road_audit, summary=audit_summary)


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


# =============================================================================
# 7. hard-structures / GROIN -- BUXTON GROIN FIELD
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
        "groin would silently do nothing.\n"
        "USE_SANDBOX_CASCADE is pinned True and has no line in hat_run.yaml, "
        "so reaching this means HAT_USE_SANDBOX_CASCADE=0 is set in the "
        "environment. Unset it, or turn the groin off.")

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


# --- 7.4 report ---------------------------------------------------------------

reports.groin_report(
    enabled=GROIN_ENABLED, callback=GROIN_CB,
    updrift_gis=GROIN_UPDRIFT_GIS, downdrift_gis=GROIN_DOWNDRIFT_GIS,
    install_year=GROIN_INSTALL_YEAR, start_year=START_YEAR, end_year=END_YEAR,
    geometry=HATTERAS_DOMAINS,
    trapping_rate_m_yr=GROIN_TRAPPING_RATE_M_YR,
    m_provenance=GROIN_M_PROVENANCE,
    deterioration_mode=GROIN_DETERIORATION_MODE,
    deterioration_delay_years=GROIN_DETERIORATION_DELAY_YEARS,
    deterioration_ramp_years=GROIN_DETERIORATION_RAMP_YEARS,
    deterioration_fraction=GROIN_DETERIORATION_FRACTION,
    profile_height_candidates_m=GROIN_PROFILE_HEIGHT_CANDIDATES_M,
    reach_transport_loss_m3_yr=REACH_TRANSPORT_LOSS_M3_YR,
    reach_transport_citation=REACH_TRANSPORT_CITATION,
    reach_transport_caveat=REACH_TRANSPORT_CAVEAT,
    source_sink_preset=SOURCE_SINK_PRESET, domain_be_rates=DOMAIN_BE_RATES)


# --- 7.5 scenario summary, and the run name derived from it ------------------
# RUN_NAME_SUFFIX used to be typed by hand: run with the groin off, forget to
# retype the label, and the output lands in a directory named for a different
# experiment. It is now derived from the switches themselves.

SCENARIO_SWITCHES = [
    ("period", f"{START_YEAR}-{END_YEAR} ({RUN_YEARS} yr)", None),
    ("source/sink preset", SOURCE_SINK_PRESET, SOURCE_SINK_PRESET),
    # No token when "asrun": every run predating the shoreline_offset
    # unit finding used it, and adding a token would rename them all.
    ("shoreline offset", OFFSET_MODE,
     None if OFFSET_MODE == "asrun" else f"offset{OFFSET_MODE}"),
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

reports.scenario_summary_report(
    scenario=SCENARIO, departures=_SCENARIO_DEPARTURES,
    switches=SCENARIO_SWITCHES, run_name_base=RUN_NAME_BASE,
    double_managed=DOUBLE_MANAGED_GIS, groin_enabled=GROIN_ENABLED,
    updrift_gis=GROIN_UPDRIFT_GIS,
    beach_dune_on_updrift=BEACH_DUNE_MANAGEMENT_ON[
        HATTERAS_DOMAINS.gis_to_pad(GROIN_UPDRIFT_GIS)])


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


COASTSAT_TARGET = build_target_table(
    CS_ACTIVE, LOESS_CONFIG, HATTERAS_DOMAINS, TARGET_WINDOW)


# --- 8.4 report ---------------------------------------------------------------

reports.coastsat_report(
    target=COASTSAT_TARGET, active=CS_ACTIVE, target_window=TARGET_WINDOW,
    loess_config=LOESS_CONFIG, geometry=HATTERAS_DOMAINS, cs_series=cs_series,
    updrift_gis=GROIN_UPDRIFT_GIS, downdrift_gis=GROIN_DOWNDRIFT_GIS)


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

# WHICH ESTIMATOR THE FIGURES DRAW. "lrr" is the OLS slope through every
# annual state; "endpoint" is (x[-1] - x[0]) / RUN_YEARS. The observational
# target is an LRR -- CoastSat's transect_lrr_full.csv is a per-transect OLS
# slope with r_squared and unc_m_yr beside it -- so "lrr" is the only setting
# that puts the same quantity on both sides of the comparison. "endpoint"
# exists to redraw a pre-2026-08-22 figure, not as an alternative.
RATE_ESTIMATOR = "lrr"
if RATE_ESTIMATOR not in ("lrr", "endpoint"):
    raise ValueError(f"RATE_ESTIMATOR must be 'lrr' or 'endpoint', "
                     f"got {RATE_ESTIMATOR!r}")

# Below this, a domain's LRR is reported as a poor summary of its trajectory
# rather than passed over. 0.50 is a reporting threshold and nothing else --
# no value is dropped, filtered, or flagged in the CSV on account of it.
LRR_R2_FLOOR = 0.50


# --- 9.2 animation jobs -------------------------------------------------------
# range: "real" | "all" | "groin" | "groin_span" | (gis_lo, gis_hi)
# mode:  "position" | "displacement" | "difference"
# pad:   half-width in domains, read only by "groin" / "groin_span"
#
# "groin" fans out into one GIF per structure in annotations.groins, so new
# structures are picked up without editing this list.

# `output.make_gifs: false` in hat_run.yaml empties this list rather than
# skipping the section 12 call. The shoreline matrix .npy is written by that
# same call, OUTSIDE the job loop, and section 12.3's paired groin baseline
# and HAT_scenario_grid.py both read it -- so short-circuiting the call would
# cost the run its matrix, while an empty job list costs only the animations.
GIF_JOBS = [
    dict(range="real", mode="displacement"),
    dict(range="real", mode="position"),
    dict(range="groin", mode="position", pad=9),
    dict(range="groin", mode="difference", pad=9),
] if RUN_CONFIG.make_gifs else []


# --- 9.3 baseline for difference jobs ----------------------------------------


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

reports.figure_config_report(
    annotations=HATTERAS_ANNOTATIONS, loess_config=LOESS_CONFIG,
    gif_config=gif_config, gif_jobs=GIF_JOBS,
    flip_sign_model=FLIP_SIGN_MODEL,
    real_domains_only=PLOT_REAL_DOMAINS_ONLY,
    groin_enabled=GROIN_ENABLED, baseline_npy=GIF_BASELINE_NPY,
    baseline_name=GIF_BASELINE_NAME, output_base_dir=OUTPUT_BASE_DIR)


# =============================================================================
# 10. build_cascade + run_cascade_simulation
# =============================================================================
# Both live in cascade_pipeline/hindcast.py, imported in section 1, and are
# shared with HAT_groin_sweep_worker.py -- which used to carry its own copy.
#
# The split is what lets section 11 hold a built-but-unstepped Cascade: BRIE's
# diffusivity and the groin's fillet prediction are only meaningful as initial
# conditions, and a prediction printed after the run is not one.
#
# `run_years` is TRANSITIONS, not states. See the module comment there for the
# off-by-one this replaced, which ran 19 updates for a 20-year period while
# dividing by 20.


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
# From hat_run.yaml (physics.wave_height_Hs); reaches run_index.csv as Hs_m,
# so a changed value is recoverable from the index and not only from the file.
Hs = RUN_CONFIG.hs                    # m, calibration value 2.5
FIXED_WAVE_PERIOD = 8                 # s
FIXED_WAVE_ASYMMETRY = 0.7
FIXED_WAVE_ANGLE_HIGH_FRACTION = 0.1

# --- datums -------------------------------------------------------------------
BERM_ELEVATION = 1.7    # m NAVD88, Hatteras Island, NCDOT-derived via NC State
MHW_ELEVATION = 0.36    # m NAVD88, Duck NC gauge (NOAA 8651370)

# --- sandbags: off for the hindcast ------------------------------------------
# From hat_run.yaml (top-level `sandbags`: a management decision, not a
# property of the coast); reaches run_index.csv as sandbags_on.
ENABLE_SANDBAG_PLACEMENT = RUN_CONFIG.sandbags
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
    shoreline_offset=island_offset,
    wave_height=Hs,
    wave_period=FIXED_WAVE_PERIOD,
    wave_asymmetry=FIXED_WAVE_ASYMMETRY,
    wave_angle_high_fraction=FIXED_WAVE_ANGLE_HIGH_FRACTION,
    berm_elevation=BERM_ELEVATION,
    MHW=MHW_ELEVATION,
    data_base=HATTERAS_DATA_BASE,
    parameter_file=PARAMETER_FILE,
    groin_callback=GROIN_CALLBACK,
)


# --- 11.3 pre-run diagnostics, and the groin prediction ----------------------
# Reads the constructed model rather than the yaml, which resolves the active
# profile height section 7 had to leave open. The predicted fillet amplitude and
# extent are written down BEFORE the run; section 12 checks the emergent extent
# against them. Amplitude was tuned, extent was not.


R_IPL = brie_r_ipl(cascade)
_brie = cascade._brie_coupler._brie
_d_sf_m = float(_brie.d_sf)
_h_b_m = float(cascade.barrier3d[0].h_b_TS[0]) * DAM_TO_M
_profile_height_m = _d_sf_m + _h_b_m
_berm_floor_m = float(cascade.barrier3d[0].BermEl) * DAM_TO_M

if GROIN_CALLBACK is not None:
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
    save_model_state=RUN_CONFIG.save_model_state,
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

reports.run_length_report(states=_states, run_years=RUN_YEARS)

# The denominator that produced a 5% bias in every earlier run. Checked, not
# trusted -- see the run_years comment in cascade_pipeline/hindcast.py.
assert _states - 1 == RUN_YEARS, (
    f"span mismatch: {_states} states implies {_states - 1} years elapsed, "
    f"but RUN_YEARS is {RUN_YEARS}")
# TWO ESTIMATORS OF THE SAME THING, and they are not interchangeable.
#
#   change_rate  (x[-1] - x[0]) / RUN_YEARS. A net displacement over a span.
#                Reads two of the RUN_YEARS+1 states, so a single-year
#                excursion still present in the final state arrives at full
#                amplitude. Kept because it is the only column that exactly
#                conserves the period's net shoreline movement, and because
#                every run before 2026-08-22 was scored on it.
#   model_lrr    OLS slope through every state. This is what section 8's
#                target IS -- CoastSat's transect_lrr_full.csv holds a
#                per-transect OLS slope -- so this is the column the figures
#                and the skill metrics use.
#
# The distinction is not cosmetic here. A nourishment fill enters as an
# instantaneous step in x_s, and BRIE's Crank-Nicolson alongshore solve
# answers a step with a grid-scale mode that alternates sign along the coast
# and decays only ~39%/yr (see compute_lrr). The 2022 Avon and Buxton fills
# are two years from the end of the 2004-2024 run, so an endpoint rate
# reports that ringing as a +-1.7 m/yr sawtooth through GIS 6-15 and 22-28.
# The LRR reports about a quarter of it, and the residue is the real fill
# edge rather than the solver.
change_rate = compute_change_rate(
    shoreline_m, span_years=RUN_YEARS, flip_sign=FLIP_SIGN_MODEL)
model_lrr, model_lrr_r2 = compute_lrr(
    shoreline_m, span_years=RUN_YEARS, flip_sign=FLIP_SIGN_MODEL)

# One name for whichever estimator section 12.5 draws, resolved once here so
# no figure call picks its own.
PLOTTED_RATE = model_lrr if RATE_ESTIMATOR == "lrr" else change_rate

# --- nourishment: the model's own record, not the schedule's intent ---------
BN_REPORT = nourishment.verify_nourishment(
    cascade, BN_SCHEDULE_APPLIED, BEACH_DUNE_MANAGEMENT_ON)
reports.nourishment_report(report=BN_REPORT, run_years=RUN_YEARS)
assert BN_REPORT["ok"], "nourishment did not reach the model as scheduled"

# --- the double-management consequence section 6 predicted ------------------
reports.frozen_setbacks_report(
    double_managed=DOUBLE_MANAGED_GIS,
    rows=(nourishment.verify_setbacks_frozen(
              cascade, DOUBLE_MANAGED_GIS, HATTERAS_DOMAINS)
          if DOUBLE_MANAGED_GIS else ()))

# --- which road_offset survived ---------------------------------------------
ROAD_SUMMARY = roadway.summarise_road_management(
    cascade, HATTERAS_DOMAINS, HATTERAS_FIRST_ROAD_DOMAIN,
    HATTERAS_LAST_ROAD_DOMAIN)
_drowned, _blocked = reports.roadway_outcome_report(summary=ROAD_SUMMARY)


# --- 12.3 groin: the pre-registered extent check -----------------------------


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


# --- 12.4 write ---------------------------------------------------------------

_shoreface_depth_m = float(cascade._brie_coupler._brie.d_sf)

run = RunInfo(
    run_name=RUN_NAME, run_dir=RUN_DIR,
    start_year=START_YEAR, end_year=END_YEAR, Hs=Hs,
    flip_sign_model=FLIP_SIGN_MODEL,
    background_erosion_on=USE_BACKGROUND_EROSION,
)

_rate_csv = os.path.join(RUN_DIR, f"{RUN_NAME}_shoreline_change_rate.csv")
# Both estimators ship, so a consumer states which one it wants rather than
# inheriting whichever the pipeline happened to write. lrr_r2 rides along
# because a slope through a domain that stepped rather than trended is a
# summary worth flagging at the point of use -- the nourished domains come
# out near 0.00.
pd.DataFrame({
    "gis_domain": np.arange(HATTERAS_DOMAINS.first_gis_id,
                            HATTERAS_DOMAINS.last_gis_id + 1),
    "change_rate_m_yr": change_rate[_real],
    "lrr_m_yr": model_lrr[_real],
    "lrr_r2": model_lrr_r2[_real],
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
# SKILL is the LRR one: the target is an LRR, so this is the only pairing
# that compares like with like. SKILL_ENDPOINT is kept beside it because
# calibBE and the groin M were fit before this distinction was drawn, and a
# preset's provenance is unreadable once the metric it was fit on stops being
# recorded. Expect SKILL to be slightly WORSE than SKILL_ENDPOINT on most
# runs -- modelled erosion decelerates across a period, so an all-years slope
# is more erosive than an endpoint difference. That is the estimator changing,
# not the model.
SKILL = skill_vs_target(model_lrr, COASTSAT_TARGET, HATTERAS_DOMAINS)
SKILL_ENDPOINT = skill_vs_target(change_rate, COASTSAT_TARGET,
                                 HATTERAS_DOMAINS)
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

# Where a straight line is a poor summary of the modelled trajectory, so the
# LRR is read with that in mind rather than silently. A nourished domain sits
# flat for most of the period and then steps, which no slope describes well.
# r2 is a variance ratio, so a domain that barely moved lands here too -- the
# two cases are told apart by whether the domain took a fill.
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
        # BOTH are recorded, and the product is not optional. Version numbers
        # restart at v1 per product (2026-08-26), so "v1" alone no longer
        # identifies a topography - 1984-start/v1 and 2004-start/v1 are
        # different surfaces. A run that logs only the version cannot be traced
        # back to the arrays it read.
        "topo_product": TOPO_PRODUCT,
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
    # f was absent from this index until 2026-08-31, so no groin run before
    # that date records which deterioration fraction it used -- and the pair
    # is quoted as (M, f), not as M alone. The seed runs turned out to be
    # f = 0.9, recoverable only from a comment in
    # HAT_be_zone_LOESS_analysis.py. be1 needs no column: it is already here
    # as be_rate_gis1_m_yr.
    "groin_deterioration_f": (GROIN_CALLBACK.deterioration_fraction
                              if GROIN_CALLBACK is not None else np.nan),
    "roadway_management": ENABLE_ROADWAY_MANAGEMENT,
    "relocations_enabled": ENABLE_HISTORICAL_ROAD_RELOCATIONS,
    "beach_dune_management": ENABLE_BEACH_DUNE_MANAGEMENT,
    "nourishment_fills": ENABLE_NOURISHMENT_FILLS,
    "bdm_domains": int(sum(BEACH_DUNE_MANAGEMENT_ON)),
    "nourishment_projects": len(BN_SCHEDULE_APPLIED.projects),
    "Hs_m": Hs,
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
    PLOTTED_RATE, cs_series, run,
    real_domains_only=PLOT_REAL_DOMAINS_ONLY, estimator=RATE_ESTIMATOR,
    sea_level_rise_rate_m_yr=SEA_LEVEL_RISE_RATE,
    save_path=os.path.join(
        RUN_DIR, f"{RUN_NAME}_shoreline_change_rate"
        f"{'_REAL_DOMAINS_ONLY' if PLOT_REAL_DOMAINS_ONLY else ''}.png"),
    show=SHOW_FIGURES, **RATE_FIG_KWARGS)

plot_annotated_rate_comparison(
    PLOTTED_RATE, cs_series, run,
    estimator=RATE_ESTIMATOR,
    sea_level_rise_rate_m_yr=SEA_LEVEL_RISE_RATE,
    save_path=os.path.join(RUN_DIR, f"{RUN_NAME}_annotated.png"),
    show=SHOW_FIGURES, **RATE_FIG_KWARGS)

# Section 9.4's validation target, resolved now that the run has a year 0.
# Buffers are NaN, so the comparison below is over the real domains only.
SHORELINE_TARGET_M, OBSERVED_CHANGE_M = build_shoreline_target(
    shoreline_m[0], START_YEAR, END_YEAR, HATTERAS_DOMAINS, RAW_OFFSET_DIR)

reports.target_misfit_report(
    target_m=SHORELINE_TARGET_M, observed_change_m=OBSERVED_CHANGE_M,
    shoreline_m=shoreline_m, end_year=END_YEAR, geometry=HATTERAS_DOMAINS,
    raw_offset_dir=RAW_OFFSET_DIR)

# Roadway relocations, so the shoreline GIFs mark the year the road moved.
# Assembled here rather than inside the plotting module: the module takes an
# array and should not know that a relocation lives on a RoadwayManager.
# None when the run has no roadways, which draws no markers and leaves the
# GIFs byte-for-byte as they were.
_RELOCATION_EVENTS = None
if getattr(cascade, "_roadways", None):
    _RELOCATION_EVENTS = np.zeros(
        (shoreline_m.shape[0], HATTERAS_DOMAINS.total_domains), dtype=bool)
    for _pad, _roadway in enumerate(cascade._roadways):
        # NOT `getattr(...) or []`: the attribute is a numpy array and
        # `array or []` raises "truth value of an array is ambiguous".
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
