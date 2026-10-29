#!/usr/bin/env python3
"""
Everything a Hatteras hindcast run needs before the model is built, as one step per section.

    from HAT_hindcast_setup import prepare_run
    r = prepare_run(RUN_CONFIG)          # every step in order; r.RUN_NAME_BASE, r.ELEVATION_FILE_PATHS, ...

Each setup_* function is one section of the runner (2 to 11.1), reads what the earlier ones put on `r`,
adds its own, and prints its report, so HAT_hindcast.ipynb can call them one at a time. Moved out of
HAT_hindcast.py on 2026-10-23, unchanged. Import it after HAT_hindcast.py's section 1 (the CASCADE
patches and the sandbox switch are settled there). Details: scripts/real-world-cases/Hatteras/hatteras_ms/README.md.

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-10-23
"""

import os
from pathlib import Path
from types import SimpleNamespace


from hatteras_paths import WORKSPACE as PROJECT_BASE_DIR

from cascade_pipeline import nourishment, roadway
from cascade_pipeline import reports
from cascade_pipeline.hindcast import (
    build_background_erosion,
    build_domain_file_paths,
    build_target_table,
    load_barrier3d_contract,
    build_island_offset,
    island_offset_tilts,
    load_storm_series,
    scenario_run_name,
    wave_climate_token,
    relocation_setback_token,
)
from cascade_pipeline.coastsat_lowess import (
    CoastSatDataset,
    LowessConfig,
    build_coastsat_series,
)
from cascade_pipeline.run_layout import resolve
from cascade_pipeline.run_names import (
    SITE as RUN_NAME_SITE, build_run_name, other_spellings, scenario_field)
from cascade_pipeline.run_registry import (
    barrier3d_provenance,
    sweep_family,
    preset_dir_for,
)

from site_layer.hatteras_site_config import (
    first_model_year,
    last_model_year,
    run_years,
    HATTERAS_BEACH_DUNE,
    be_rates,
    HATTERAS_BE_EDGE_DOMAINS,
    HATTERAS_BE_RATES_2004_IS_PLACEHOLDER,
    HATTERAS_COMMUNITY_ZONES,
    HATTERAS_DOMAINS,
    HATTERAS_FIRST_ROAD_DOMAIN,
    HATTERAS_LAST_ROAD_DOMAIN,
    HATTERAS_NOURISHMENT_PROJECTS,
    HATTERAS_PERIODS,
    HATTERAS_ROAD_ELEVATION_FILE,
    HATTERAS_ROAD_EVENTS,
    island_offset_version,
    resolve_be_preset,
    HATTERAS_GEOMETRY,
    HATTERAS_GEOMETRY_EXTENDED,
    HATTERAS_OFFSET_SOURCE,
)
from site_layer.hat_topo_version import RUN_OFFSET_SOURCE  # noqa: E402
from cascade_pipeline.domains import DEFAULT_DOMAINS  # the surveyed reach, GIS 1-90
# Which extraction: resolved by topo_dirs(); an old one only by override, e.g. topo_dirs("2009-start", override="v3")
from site_layer.hat_topo_version import topo_dirs
from site_layer.hat_topo_version import BUFFER_DIR as _BUFFER_DIR
# Where the rate fits live is hat_observed_rates' to say
from site_layer.hat_observed_rates import COASTSAT_LRR_ROOT as COASTSAT_BASE_DIR
from HAT_hindcast_config import (
    describe as _describe_run_config, preflight as _preflight,
    field_default as _field_default, groin_strength_for)
import inspect
import cascade_pipeline.hindcast as _hindcast_module


# Named scenarios: each sets the roadway, beach_dune, fills and relocations switches.
# Literal dict(...) calls at module level: HAT_run_matrix.py reads this table from the source
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


# What the steps build, one attribute per value, named as the runner named it
class RunSetup(SimpleNamespace):
    pass


# Section 2: the topography product, its padded file lists, and the units check
def setup_topography(r):
    # 2. Dune/topo -- per period

    # The topography product is per period (HATTERAS_PERIODS); its version is resolved, never pinned

    # 2.1 Project paths

    # No per-period array names here: build_domain_file_paths() delegates to the resolver

    r.TOPO_PRODUCT = HATTERAS_PERIODS[r.RUN_CONFIG.start_year]["topo_product"]

    _TOPO_DIR, _DUNE_DIR, r.TOPO_DUNE_VERSION = topo_dirs(r.TOPO_PRODUCT)

    print(f"topography            {r.TOPO_PRODUCT} / {r.TOPO_DUNE_VERSION}  "
          f"(product from the period, version resolved)")

    r.HATTERAS_DATA_BASE = PROJECT_BASE_DIR / "data" / "hatteras_init"
    r.OUTPUT_ROOT = PROJECT_BASE_DIR / "output" / "raw_runs"
    r.PARAMETER_FILE = "Hatteras-CASCADE-parameters.yaml"  # resolved by CASCADE

    # Taken from what topo_dirs() returned, never re-joined from parts
    r.DUNE_TOPO_DIR = _TOPO_DIR.parent
    r.BUFFER_DIR = _BUFFER_DIR

    os.chdir(PROJECT_BASE_DIR)
    # Runs are filed per period (OUTPUT_BASE_DIR, section 3)
    r.OUTPUT_ROOT.mkdir(parents=True, exist_ok=True)

    reports.path_inventory([
        ("PROJECT_BASE_DIR", PROJECT_BASE_DIR),
        ("HATTERAS_DATA_BASE", r.HATTERAS_DATA_BASE),
        ("DUNE_TOPO_DIR", r.DUNE_TOPO_DIR),
        ("BUFFER_DIR", r.BUFFER_DIR),
        ("COASTSAT_BASE_DIR", COASTSAT_BASE_DIR),
        ("OUTPUT_ROOT", r.OUTPUT_ROOT),
    ])

    # 2.2 Build the padded file lists

    # The product is passed; hat_topo_version resolves directory and filename together
    r.ELEVATION_FILE_PATHS, r.DUNE_FILE_PATHS = build_domain_file_paths(
        HATTERAS_DOMAINS, r.TOPO_PRODUCT)

    print(f"\n{len(r.ELEVATION_FILE_PATHS)} elevation + {len(r.DUNE_FILE_PATHS)} dune "
          f"paths (expect {HATTERAS_DOMAINS.total_domains} each)")

    # 2.3 Verify every file exists

    # A stale version or a moved data folder fails here, with a count and the first offender

    _expected_files = 2 * HATTERAS_DOMAINS.total_domains
    _missing = [path for path in r.ELEVATION_FILE_PATHS + r.DUNE_FILE_PATHS
                if not Path(path).exists()]

    if _missing:
        raise FileNotFoundError(
            f"{len(_missing)} of {_expected_files} init files missing. Check "
            f"TOPO_DUNE_VERSION ({r.TOPO_DUNE_VERSION!r}) and DUNE_TOPO_DIR.\n"
            f"  First missing: {_missing[0]}")

    print(f"All {_expected_files} init files present.")

    # 2.4 Units check against Barrier3D's input contract

    # The .npy arrays must already be decameters MHW (a metres file runs 10x too tall); checked raw

    r.BARRIER3D_CONTRACT = load_barrier3d_contract(r.HATTERAS_DATA_BASE / r.PARAMETER_FILE)

    print(f"\nContract from {r.PARAMETER_FILE}:")
    print(f"  BarrierLength -> {r.BARRIER3D_CONTRACT['barrier_length_cells']} "
          f"alongshore cells")
    print(f"  MHW           -> {r.BARRIER3D_CONTRACT['mhw_dam']:.3f} dam")
    print(f"  BermEl        -> {r.BARRIER3D_CONTRACT['berm_el_dam']:.3f} dam "
          f"above MHW\n")

    reports.run_units_check(r.ELEVATION_FILE_PATHS, r.DUNE_FILE_PATHS,
                            r.BARRIER3D_CONTRACT, HATTERAS_DOMAINS)


# Section 3: the period, the scenario switches, where the run is filed, and the run-name preview
def setup_scenario(r):
    # 3. Island orientation -- set START_YEAR

    # START_YEAR selects the period; values come from hat_run.yaml, re-read here, and can be typed over

    # --- CONFIG ------------------------------------------------------------------
    r.START_YEAR = r.RUN_CONFIG.start_year   # 1984, 1996, 2004 or 2009

    # Source/sink preset: "zeroBE" nothing, "edgeBE" the two end domains, "calibBE" the full fit
    r.SOURCE_SINK_PRESET = r.RUN_CONFIG.source_sink_preset

    # Scenario -- the management combination this run simulates

    r.SCENARIOS = SCENARIOS

    r.SCENARIO = r.RUN_CONFIG.scenario

    # The reach: must match the HAT_GEOMETRY hatteras_site_config was built from ("base" = GIS 1-90)
    r.GEOMETRY = r.RUN_CONFIG.geometry
    if r.GEOMETRY != HATTERAS_GEOMETRY:
        raise SystemExit(
            f"\n[stop] geometry mismatch: the config says {r.GEOMETRY!r} but "
            f"hatteras_site_config built {HATTERAS_GEOMETRY!r} from the "
            f"environment. Set HAT_GEOMETRY={r.GEOMETRY} in the environment "
            f"(the yaml alone cannot select a reach).\n")

    # The groin is not part of the scenario: every scenario runs with and without it (12.3)
    r.GROIN_ENABLED = r.RUN_CONFIG.groin_enabled

    # One groin since 2026-10-14, the corrected groin; its runs carry the plain token
    r.GROIN_TOKEN = "groin" if r.GROIN_ENABLED else "nogroin"

    # Where a relocated road is rebuilt, metres behind the dune line; "measured" = the t=0 setback
    r.RELOCATION_SETBACK_M = r.RUN_CONFIG.relocation_setback_m

    # The road manager rebuilds dunes only where the road is within this many metres of the dune (obx-cascade)
    r.ROAD_DUNE_REBUILD_TRIGGER_M = r.RUN_CONFIG.road_dune_rebuild_trigger_m

    # The beach/dune manager clips village dunes to this every year, m NAVD88 (lidar-measured, 2026-10-29)
    r.VILLAGE_DUNE_CAP_NAVD88_M = r.RUN_CONFIG.village_dune_cap_navd88_m
    # -----------------------------------------------------------------------------

    # Sensitivity tokens

    # Forcing tokens keep a sensitivity cell out of the matrix run's name; None at the calibration values
    r.WAVE_VALUES = {
        "hs": r.RUN_CONFIG.hs,
        "wave_period_s": r.RUN_CONFIG.wave_period_s,
        "wave_asymmetry": r.RUN_CONFIG.wave_asymmetry,
        "wave_angle_high_fraction": r.RUN_CONFIG.wave_angle_high_fraction,
    }
    _WAVE_DEFAULTS = {name: _field_default(name) for name in r.WAVE_VALUES}
    # A sensitivity cell is its baseline's name plus this token (filed under sensitivity/<axis>/)
    r.WAVE_TOKEN = wave_climate_token(r.WAVE_VALUES, _WAVE_DEFAULTS)
    r.RELOCATION_SETBACK_TOKEN = relocation_setback_token(
        r.RELOCATION_SETBACK_M, _field_default("relocation_setback_m"))

    print("\n" + _describe_run_config())

    # Expand
    if r.SCENARIO not in r.SCENARIOS:
        raise ValueError(f"SCENARIO must be one of {sorted(r.SCENARIOS)}, "
                         f"got {r.SCENARIO!r}")
    _SCENARIO_PRESET = r.SCENARIOS[r.SCENARIO]
    r.ENABLE_ROADWAY_MANAGEMENT = _SCENARIO_PRESET["roadway"]
    r.ENABLE_BEACH_DUNE_MANAGEMENT = _SCENARIO_PRESET["beach_dune"]
    r.ENABLE_NOURISHMENT_FILLS = _SCENARIO_PRESET["fills"]
    # HAT_RELOCATIONS overrides the scenario when set; None leaves the preset in charge
    r.ENABLE_HISTORICAL_ROAD_RELOCATIONS = (
        _SCENARIO_PRESET["relocations"] if r.RUN_CONFIG.relocations is None
        else r.RUN_CONFIG.relocations)

    # One-off overrides

    # Uncomment to depart from the named scenario for one run; the departure is detected and named below
    # ENABLE_ROADWAY_MANAGEMENT = False
    # ENABLE_BEACH_DUNE_MANAGEMENT = False
    # ENABLE_NOURISHMENT_FILLS = False
    # ENABLE_HISTORICAL_ROAD_RELOCATIONS = True

    r.SCENARIO_DEPARTURES = {
        key: (want, got) for key, want, got in (
            ("roadway", _SCENARIO_PRESET["roadway"], r.ENABLE_ROADWAY_MANAGEMENT),
            ("beach_dune", _SCENARIO_PRESET["beach_dune"],
             r.ENABLE_BEACH_DUNE_MANAGEMENT),
            ("fills", _SCENARIO_PRESET["fills"], r.ENABLE_NOURISHMENT_FILLS),
            ("relocations", _SCENARIO_PRESET["relocations"],
             r.ENABLE_HISTORICAL_ROAD_RELOCATIONS),
        ) if want != got
    }
    if r.START_YEAR not in HATTERAS_PERIODS:
        raise ValueError(f"START_YEAR must be one of {sorted(HATTERAS_PERIODS)}, "
                         f"got {r.START_YEAR}")

    # Normalised to the canonical preset name, so an alias never reaches RUN_NAME
    r.SOURCE_SINK_PRESET, _PRESET_BY_PERIOD = resolve_be_preset(r.SOURCE_SINK_PRESET)

    r.PERIOD = HATTERAS_PERIODS[r.START_YEAR]

    r.END_YEAR = r.PERIOD["end_year"]
    # end_year is the window label; the run steps start..last_model_year (scripts/real-world-cases/Hatteras/hatteras_ms/docs/MODEL_YEARS.md)
    r.LAST_MODEL_YEAR = last_model_year(r.START_YEAR)
    # The period key is the DEM year; the model clock starts at the first storms (since 2026-10-12)
    r.FIRST_MODEL_YEAR = first_model_year(r.START_YEAR)
    r.RUN_YEARS = run_years(r.START_YEAR)

    r.SEA_LEVEL_RISE_RATE = r.PERIOD["sea_level_rise_rate"]
    r.ENABLE_NOURISHMENT = r.PERIOD["enable_nourishment"]
    r.NOURISHMENT_VOLUME = r.PERIOD["nourishment_volume"]
    r.ISLAND_OFFSET_FILE = r.HATTERAS_DATA_BASE / r.PERIOD["island_offset_file"]
    # Which version folder that resolved to (v1, v2, or flat); see the json write.
    r.ISLAND_OFFSET_VERSION = island_offset_version(r.START_YEAR)
    r.STORM_FILE = r.HATTERAS_DATA_BASE / r.PERIOD["storm_file"]
    r.ROAD_SETBACK_FILE = r.HATTERAS_DATA_BASE / r.PERIOD["road_setback_file"]

    # Resolve the combinations that cannot both be true

    # Contradictory switches (fill without a manager, relocation without a road) are resolved and announced
    _FILLS_FORCED_OFF = (r.ENABLE_NOURISHMENT_FILLS
                         and not r.ENABLE_BEACH_DUNE_MANAGEMENT)
    if _FILLS_FORCED_OFF:
        r.ENABLE_NOURISHMENT_FILLS = False
    _RELOCATIONS_FORCED_OFF = (r.ENABLE_HISTORICAL_ROAD_RELOCATIONS
                               and not r.ENABLE_ROADWAY_MANAGEMENT)
    if _RELOCATIONS_FORCED_OFF:
        r.ENABLE_HISTORICAL_ROAD_RELOCATIONS = False

    # Run name: the period stem; the scenario suffix is derived from the switches in 7.5
    r.RUN_NAME_STEM = f"{RUN_NAME_SITE}_{r.START_YEAR}-{r.END_YEAR}"

    # Run directories are filed by period, so every later path (RUN_DIR, the 12.3 baseline) is too
    r.PERIOD_TAG = f"{r.START_YEAR}_{r.END_YEAR}"
    # Filed by purpose: matrix/, sensitivity/<axis>/, experiments/<tag>/, versions/<tag>/ (run_registry)
    r.RUN_KIND = r.RUN_CONFIG.run_kind
    r.RUN_TAG = r.RUN_CONFIG.run_tag
    # HAT_ARM_TAG, the pre-09-16 spelling, is read as an experiment tag
    _LEGACY_ARM = os.environ.get("HAT_ARM_TAG", "").strip()
    if _LEGACY_ARM and not r.RUN_TAG:
        print(f"HAT_ARM_TAG={_LEGACY_ARM!r} is the pre-2026-09-16 spelling; read as "
              f"HAT_RUN_KIND=experiment HAT_RUN_TAG={_LEGACY_ARM!r}")
        r.RUN_KIND, r.RUN_TAG = "experiment", _LEGACY_ARM
    # A run off the calibration waves cannot be filed in matrix/: refused
    if r.RUN_KIND == "matrix" and r.WAVE_TOKEN:
        raise ValueError(
            f"wave climate is off calibration ({r.WAVE_TOKEN}) but HAT_RUN_KIND is "
            f"matrix. A forced run is a sensitivity cell (HAT_RUN_KIND=sensitivity) "
            f"or an experiment (HAT_RUN_KIND=experiment HAT_RUN_TAG=<name>).")
    # Nor can a run on a non-default offset source (its name would collide): refused
    if r.RUN_KIND == "matrix" and HATTERAS_OFFSET_SOURCE != RUN_OFFSET_SOURCE:
        raise ValueError(
            f"island offset source is {HATTERAS_OFFSET_SOURCE!r} (not "
            f"{RUN_OFFSET_SOURCE!r}) but HAT_RUN_KIND is matrix. The run name "
            f"carries no offset token, so this would overwrite the matrix run it is "
            f"meant to be compared against. Use HAT_RUN_KIND=experiment with "
            f"HAT_RUN_TAG=<name>.")
    if r.RUN_KIND == "sensitivity" and not (r.WAVE_TOKEN or r.RELOCATION_SETBACK_TOKEN):
        raise ValueError(
            "HAT_RUN_KIND=sensitivity but every swept forcing is at its default; "
            "this cell would be the matrix run under another name.")
    if r.RUN_KIND == "sensitivity" and not r.RUN_TAG:
        r.RUN_TAG = sweep_family(f"x_{r.WAVE_TOKEN or r.RELOCATION_SETBACK_TOKEN}")
    # Model state (~99% of a run's size) kept for matrix runs unless output.save_model_state says otherwise
    r.SAVE_MODEL_STATE = (r.RUN_CONFIG.save_model_state
                        if (r.RUN_KIND == "matrix"
                            or r.RUN_CONFIG.origins["save_model_state"] != "default")
                        else False)
    # Built by run_registry, so the writer and the readers agree on the layout
    r.OUTPUT_BASE_DIR = preset_dir_for(r.OUTPUT_ROOT, r.PERIOD_TAG, r.SOURCE_SINK_PRESET,
                                     kind=r.RUN_KIND, tag=r.RUN_TAG)
    r.OUTPUT_BASE_DIR.mkdir(parents=True, exist_ok=True)

    print(f"\nSTART_YEAR = {r.START_YEAR}  ->  {r.START_YEAR}-{r.END_YEAR}, "
          f"{r.RUN_YEARS} model years, storms {r.FIRST_MODEL_YEAR}-{r.LAST_MODEL_YEAR}")
    print(f"RUN_NAME_STEM = {r.RUN_NAME_STEM!r}"
          "   (scenario suffix derived in 7.5)")
    print(f"SOURCE_SINK_PRESET = {r.SOURCE_SINK_PRESET!r}")
    print(f"OUTPUT_BASE_DIR = {r.OUTPUT_BASE_DIR}")
    print(f"RUN_KIND = {r.RUN_KIND!r}   RUN_TAG = {r.RUN_TAG!r}   SAVE_MODEL_STATE = {r.SAVE_MODEL_STATE}")

    # Which Barrier3D: its checked-out branch is the model; recorded, warned about if not the fix branch
    r.B3D = barrier3d_provenance()
    print(f"BARRIER3D = {r.B3D['branch']}@{str(r.B3D['commit'])[:7]}   "
          f"route_overwash fix: {r.B3D['route_overwash_fix']}"
          + ("   (tree dirty)" if r.B3D["dirty"] else ""))
    print(f"            overwash gap/momentum fixes: {r.B3D.get('gap_momentum_fix')}   "
          f"per-cell dune ceilings: {r.B3D.get('per_cell_ceiling')}")
    # One island-wide Dmaxel since 2026-10-13, so a Barrier3D without per-cell ceilings runs (it was refused)
    if r.B3D.get("gap_momentum_fix") is not True:
        print("  WARNING: this Barrier3D does not carry the 2026-09-28 overwash gap and "
              "momentum fixes (hatteras/adopted has them).")
    if r.B3D["route_overwash_fix"] is not True:
        print("  WARNING: this Barrier3D does NOT carry the route_overwash index fix. "
              "Runs read the wrong cells in overwash and can crash silently. Check out "
              "fix/route-overwash-axis-swap in the Barrier3D repository unless this run "
              "reproduces an old one on purpose.")

    # The name this scenario will produce, predicted from the switches

    # Advisory preview of the run name; 7.5 derives the real one and raises if they differ
    _PERIOD_HAS_FILL = bool(nourishment.build_schedule(
        HATTERAS_NOURISHMENT_PROJECTS, HATTERAS_DOMAINS,
        r.FIRST_MODEL_YEAR, r.LAST_MODEL_YEAR).projects)
    # Grammar in cascade_pipeline/run_names.py: HAT_<start>-<end>_<scenario>_<groin>_<preset>[_<perturbation>]
    r.RUN_NAME_PREVIEW = build_run_name(
        r.START_YEAR, r.END_YEAR,
        roadway=r.ENABLE_ROADWAY_MANAGEMENT,
        beach_dune=r.ENABLE_BEACH_DUNE_MANAGEMENT,
        groin=r.GROIN_ENABLED, preset=r.SOURCE_SINK_PRESET,
        fills_withheld=(_PERIOD_HAS_FILL and r.ENABLE_BEACH_DUNE_MANAGEMENT
                        and not r.ENABLE_NOURISHMENT_FILLS),
        relocations=r.ENABLE_HISTORICAL_ROAD_RELOCATIONS,
        perturbations=(r.RELOCATION_SETBACK_TOKEN, r.WAVE_TOKEN))

    # Name, destination, collision and rough run time, reported before sections 4-10 run
    print("\n" + _preflight(r.RUN_NAME_PREVIEW,
                            r.OUTPUT_BASE_DIR / r.RUN_NAME_PREVIEW,
                            config=r.RUN_CONFIG))

    reports.scenario_report(
        scenario=r.SCENARIO, scenarios=r.SCENARIOS, departures=r.SCENARIO_DEPARTURES,
        roadway_on=r.ENABLE_ROADWAY_MANAGEMENT,
        relocations_on=r.ENABLE_HISTORICAL_ROAD_RELOCATIONS,
        relocations_forced_off=_RELOCATIONS_FORCED_OFF,
        beach_dune_on=r.ENABLE_BEACH_DUNE_MANAGEMENT,
        fills_on=r.ENABLE_NOURISHMENT_FILLS, fills_forced_off=_FILLS_FORCED_OFF,
        period_expects_nourishment=r.ENABLE_NOURISHMENT,
        groin_enabled=r.GROIN_ENABLED, run_name_preview=r.RUN_NAME_PREVIEW,
        input_files=[("island offset", r.ISLAND_OFFSET_FILE),
                     ("storms", r.STORM_FILE),
                     ("road setback", r.ROAD_SETBACK_FILE)])


# Section 3.1: each padded domain's cross-shore starting position
def setup_island_offset(r):
    # 3.1 Island offsets

    # Each padded domain's cross-shore starting position: shoreline_offset on the Cascade() call

    # Metres, as the file holds them (always since 2026-10-23)
    r.island_offset = build_island_offset(r.ISLAND_OFFSET_FILE, HATTERAS_DOMAINS)
    r.OFFSET_TILTS = island_offset_tilts(r.island_offset, HATTERAS_DOMAINS)

    r.real = slice(HATTERAS_DOMAINS.start_real_index, HATTERAS_DOMAINS.end_real_index)
    print(f"\n{r.START_YEAR} offsets (metres): {r.island_offset.size} padded domains | "
          f"span {r.island_offset[r.real].min():.0f}-{r.island_offset[r.real].max():.0f} m")


# Section 4: sea-level rise, the storm series and the source/sink rates
def setup_forcings(r):
    # 4. Period forcings -- RSLR, storms, SOURCE/SINK

    # Everything here is resolved by the START_YEAR set in section 3.

    # 4.1 Relative sea level rise

    print(f"\nSEA_LEVEL_RISE_RATE = {r.SEA_LEVEL_RISE_RATE} m/yr")
    print(f"  over {r.RUN_YEARS} years -> "
          f"{r.SEA_LEVEL_RISE_RATE * r.RUN_YEARS:.3f} m total rise")

    # 4.2 Storm series

    # One row per storm: time step, Rhigh, Rlow (dam), period, duration (h)

    r.STORM_SERIES = load_storm_series(r.STORM_FILE)

    reports.storm_report(storms=r.STORM_SERIES, storm_file=r.STORM_FILE,
                         run_years=r.RUN_YEARS)

    # 4.3 Source/sink (background erosion)

    # Per-domain rate, m/yr, as Barrier3D's Rat: (-) erosion, (+) accretion; absent domains 0.0

    # Names the missing fit; an unsolved period's first edge probe starts empty when HAT_BE_OVERRIDE is set
    try:
        r.DOMAIN_BE_RATES = be_rates(r.SOURCE_SINK_PRESET, r.START_YEAR)
    except ValueError:
        if (r.SOURCE_SINK_PRESET == "edgeBE"
                and os.environ.get("HAT_BE_OVERRIDE", "").strip()):
            r.DOMAIN_BE_RATES = {}
        else:
            raise

    # HAT_BE_OVERRIDE: a solve step without editing the config; the BE digest still tells it apart
    _BE_OVERRIDE_RAW = os.environ.get("HAT_BE_OVERRIDE", "").strip()
    if _BE_OVERRIDE_RAW:
        # Copied: the preset dict is the config's own
        r.DOMAIN_BE_RATES = dict(r.DOMAIN_BE_RATES)
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
              f"preset {r.SOURCE_SINK_PRESET!r}")
        for _gis, _rate in sorted(_BE_OVERRIDES.items()):
            _was = r.DOMAIN_BE_RATES.get(_gis, 0.0)
            print(f"  GIS {_gis:<3}           {_was:+.4f} -> {_rate:+.4f} m/yr")
            r.DOMAIN_BE_RATES[_gis] = _rate

    # An edgeBE end with no solved value and no override is refused
    if r.SOURCE_SINK_PRESET == "edgeBE":
        _no_edge = [g for g in HATTERAS_BE_EDGE_DOMAINS
                    if not r.DOMAIN_BE_RATES.get(g)]
        if _no_edge:
            raise SystemExit(
                f"\n[stop] edgeBE has no nonzero rate at end domain(s) "
                f"{_no_edge} of geometry {HATTERAS_GEOMETRY!r}. Supply one "
                f"with HAT_BE_OVERRIDE (e.g. '115=12.0'), or run zeroBE.\n")

    r.BACKGROUND_EROSION_RATES = build_background_erosion(
        r.DOMAIN_BE_RATES, HATTERAS_DOMAINS)
    r.USE_BACKGROUND_EROSION = any(rate != 0.0 for rate in r.BACKGROUND_EROSION_RATES)

    # Derived from the preset, and checked against it
    _EXPECT_BE_ON = r.SOURCE_SINK_PRESET != "zeroBE"
    if r.USE_BACKGROUND_EROSION != _EXPECT_BE_ON:
        raise ValueError(
            f"preset {r.SOURCE_SINK_PRESET!r} implies "
            f"USE_BACKGROUND_EROSION={_EXPECT_BE_ON}, but the expanded rates give "
            f"{r.USE_BACKGROUND_EROSION}. The preset in hatteras_site_config.py does "
            f"not match its name.")

    reports.background_erosion_report(
        preset=r.SOURCE_SINK_PRESET, start_year=r.START_YEAR,
        domain_rates=r.DOMAIN_BE_RATES, rates=r.BACKGROUND_EROSION_RATES,
        use_background_erosion=r.USE_BACKGROUND_EROSION, geometry=HATTERAS_DOMAINS,
        rates_2004_are_placeholder=HATTERAS_BE_RATES_2004_IS_PLACEHOLDER)


# Section 5: NC-12 setbacks, elevations, managed domains, and the year-one audit
def setup_roadway(r):
    # 5. roadway_manager -- setbacks, per-domain elevation, historical events

    # NC-12 forcing: setback (per period), elevation (m MHW), managed domains; events carry displacements

    # Loaded and audited even with management off; it just never reaches a RoadwayManager

    r.ROADWAY = roadway.RoadwayConfig()
    _road_span = (HATTERAS_FIRST_ROAD_DOMAIN, HATTERAS_LAST_ROAD_DOMAIN)

    # Setback: by period
    r.road_setbacks_full, _missing_setbacks = roadway.load_road_setbacks(
        r.ROAD_SETBACK_FILE, HATTERAS_DOMAINS, *_road_span)

    # Elevation: one set for every period
    r.ROAD_ELEVATION_FILE = r.HATTERAS_DATA_BASE / HATTERAS_ROAD_ELEVATION_FILE
    r.road_elevation_full, _missing_elevations = roadway.load_road_elevations(
        r.ROAD_ELEVATION_FILE, HATTERAS_DOMAINS, *_road_span, config=r.ROADWAY)

    # Which domains CASCADE actually manages
    r.ROADWAY_MANAGEMENT_ON = roadway.build_roadway_management_on(
        HATTERAS_DOMAINS, *_road_span,
        community_zones=HATTERAS_COMMUNITY_ZONES,
        enabled=r.ENABLE_ROADWAY_MANAGEMENT)

    _road_slice = slice(HATTERAS_DOMAINS.gis_to_pad(HATTERAS_FIRST_ROAD_DOMAIN),
                        HATTERAS_DOMAINS.gis_to_pad(HATTERAS_LAST_ROAD_DOMAIN) + 1)
    reports.roadway_report(
        setback_file=r.ROAD_SETBACK_FILE, setbacks=r.road_setbacks_full,
        missing_setbacks=_missing_setbacks,
        elevation_file=r.ROAD_ELEVATION_FILE, elevations=r.road_elevation_full,
        missing_elevations=_missing_elevations, road_slice=_road_slice,
        config=r.ROADWAY, roadway_on=r.ENABLE_ROADWAY_MANAGEMENT,
        management_on=r.ROADWAY_MANAGEMENT_ON,
        first_road_gis=HATTERAS_FIRST_ROAD_DOMAIN,
        last_road_gis=HATTERAS_LAST_ROAD_DOMAIN,
        community_zones=HATTERAS_COMMUNITY_ZONES,
        road_events=HATTERAS_ROAD_EVENTS,
        relocations_enabled=r.ENABLE_HISTORICAL_ROAD_RELOCATIONS)

    # 5.1 Pre-flight audit: which road_offset will not survive year one

    # A road whose flanking rows are >20% water drowns, and the domain is unmanaged from then on

    r.road_audit = roadway.audit_setbacks(
        r.ELEVATION_FILE_PATHS, r.road_setbacks_full, HATTERAS_DOMAINS, *_road_span,
        management_on=r.ROADWAY_MANAGEMENT_ON, config=r.ROADWAY)
    r.audit_summary = roadway.summarise_audit(r.road_audit)

    reports.road_audit_report(audit=r.road_audit, summary=r.audit_summary)


# Section 6: the nourishment schedule, overwash filter and beach/dune managed domains
def setup_beach_dune(r):
    # 6. beach_dune_manager -- nourishment schedule + overwash filter

    # beach_dune_manager: always-on filter + fixed dune line, event-driven fills; overwash_filter is a PERCENT

    # Schedule: one project list, period-filtered

    # Every Hatteras project is in 2004-2024, so a 1984 run builds an empty schedule from this list
    r.BN_SCHEDULE = nourishment.build_schedule(
        HATTERAS_NOURISHMENT_PROJECTS, HATTERAS_DOMAINS, r.FIRST_MODEL_YEAR, r.LAST_MODEL_YEAR)

    # The schedule the model is actually driven by

    # Fills off = an empty schedule, not a skipped call; the footprint is unchanged
    r.BN_SCHEDULE_APPLIED = r.BN_SCHEDULE if r.ENABLE_NOURISHMENT_FILLS else (
        nourishment.build_schedule([], HATTERAS_DOMAINS, r.FIRST_MODEL_YEAR, r.LAST_MODEL_YEAR))

    # What CASCADE is handed

    # Zeros with the module off: the array states what the run does
    r.OVERWASH_FILTER = (
        nourishment.build_overwash_filter(
            HATTERAS_DOMAINS, HATTERAS_COMMUNITY_ZONES, config=HATTERAS_BEACH_DUNE)
        if r.ENABLE_BEACH_DUNE_MANAGEMENT
        else [0.0] * HATTERAS_DOMAINS.total_domains)
    r.OVERWASH_TO_DUNE = HATTERAS_BEACH_DUNE.overwash_to_dune_pct
    # One manager per domain (obx-cascade refuses both): road domains are left to RoadwayManager, which nourishes them
    r.BEACH_DUNE_MANAGEMENT_ON = nourishment.build_beach_dune_management_on(
        HATTERAS_DOMAINS, HATTERAS_COMMUNITY_ZONES, r.BN_SCHEDULE.nourished_gis,
        enabled=r.ENABLE_BEACH_DUNE_MANAGEMENT, roadway_on=r.ROADWAY_MANAGEMENT_ON)

    # Placeholder: rewritten every year by the schedule, so a missed schedule nourishes nothing
    r.NOURISHMENT_VOLUME_INIT = [0.0] * HATTERAS_DOMAINS.total_domains

    r.DOUBLE_MANAGED_GIS = nourishment.find_double_managed(
        r.BEACH_DUNE_MANAGEMENT_ON, r.ROADWAY_MANAGEMENT_ON, HATTERAS_DOMAINS)
    r.BN_AUDIT = nourishment.audit_schedule(
        r.BN_SCHEDULE_APPLIED, r.BEACH_DUNE_MANAGEMENT_ON, config=HATTERAS_BEACH_DUNE,
        roadway_on=r.ROADWAY_MANAGEMENT_ON)

    reports.beach_dune_report(
        start_year=r.START_YEAR, end_year=r.END_YEAR,
        beach_dune_enabled=r.ENABLE_BEACH_DUNE_MANAGEMENT,
        fills_enabled=r.ENABLE_NOURISHMENT_FILLS,
        roadway_enabled=r.ENABLE_ROADWAY_MANAGEMENT,
        config=HATTERAS_BEACH_DUNE, management_on=r.BEACH_DUNE_MANAGEMENT_ON,
        overwash_to_dune=r.OVERWASH_TO_DUNE,
        community_zones=HATTERAS_COMMUNITY_ZONES,
        schedule=r.BN_SCHEDULE, schedule_applied=r.BN_SCHEDULE_APPLIED,
        audit=r.BN_AUDIT, double_managed=r.DOUBLE_MANAGED_GIS,
        geometry=HATTERAS_DOMAINS)


# Section 7: the Buxton groin's settings for the GroinManager
def setup_groin(r):
    # 7. Hard structures / groin -- Buxton groin field

    # The corrected groin moves a share of the drift arriving at GIS 5|6 updrift each year, before the transport
    # solve, and cancels BRIE's leak of it; the fillet is emergent

    # 7.1 Switches and structure

    # GROIN_ENABLED is set in section 3; the guard stays here

    # obx-cascade's Cascade has the GroinManager module (Cascade(groin_module=True, ...)); without it there is no groin
    if r.GROIN_ENABLED and "groin_module" not in inspect.signature(_hindcast_module.Cascade).parameters:
        raise RuntimeError(
            f"GROIN_ENABLED=True but {_hindcast_module.Cascade.__module__}.Cascade has no groin_module "
            "argument, so it has no GroinManager. "
            "Turn the groin off, or run against a CASCADE that has the GroinManager.")

    # Net transport is southward, so updrift = north; the two domains share the boundary at GIS 5.5
    r.GROIN_UPDRIFT_GIS = 6       # source: accretes
    r.GROIN_DOWNDRIFT_GIS = 5     # sink:   erodes
    r.GROIN_INSTALL_YEAR = 1969   # confirmed construction date

    # The groin's settings (hat_run.yaml): the share of the arriving drift it holds back, per window; the
    # first year it no longer works (the 2017 Buxton fill); and the spread, the shares of GIS 6's loss across
    # the face taken from GIS 6, 7, 8 (0, 2/3, 1/3: the spread groin, fitted with its strength on the
    # calibration window in hard-structures/groin/4-hindcast/3-wave-driven-1967-2025/2026-10-20-spread-
    # relaxation/). Another start year needs groin.strength set
    r.GROIN_STRENGTH = groin_strength_for(r.START_YEAR, r.RUN_CONFIG)
    if r.GROIN_STRENGTH is None:
        if r.GROIN_ENABLED:
            raise ValueError(
                f"no groin strength for a {r.START_YEAR} start: hat_run.yaml has groin.strength_1996 and "
                f"groin.strength_2009 only. Set groin.strength (or HAT_GROIN_STRENGTH) to run it with the groin on.")
        r.GROIN_STRENGTH = 0.0   # the groin is off; the strength history is still drawn up for the 7.4 report
    r.GROIN_END_YEAR = r.RUN_CONFIG.groin_end_year
    r.GROIN_SPREAD_WEIGHTS = r.RUN_CONFIG.groin_spread_weights

    # 7.2 The GroinManager's settings

    # Cascade builds its GroinManager from these when GROIN_ENABLED (section 11.2: groin_module=GROIN_ENABLED);
    # the strength history is drawn up either way for the 7.4 report

    r.GROIN_UPDRIFT_PAD = HATTERAS_DOMAINS.gis_to_pad(r.GROIN_UPDRIFT_GIS)
    r.GROIN_DOWNDRIFT_PAD = HATTERAS_DOMAINS.gis_to_pad(r.GROIN_DOWNDRIFT_GIS)
    r.GROIN_STRENGTH_HISTORY = [(r.GROIN_INSTALL_YEAR, r.GROIN_STRENGTH), (r.GROIN_END_YEAR, 0.0)]

    # 7.4 Report

    reports.groin_report(
        enabled=r.GROIN_ENABLED, strength_history=r.GROIN_STRENGTH_HISTORY,
        updrift_gis=r.GROIN_UPDRIFT_GIS, downdrift_gis=r.GROIN_DOWNDRIFT_GIS,
        install_year=r.GROIN_INSTALL_YEAR, start_year=r.FIRST_MODEL_YEAR, end_year=r.LAST_MODEL_YEAR + 1,
        geometry=HATTERAS_DOMAINS, strength=r.GROIN_STRENGTH, off_year=r.GROIN_END_YEAR,
        source_sink_preset=r.SOURCE_SINK_PRESET, domain_be_rates=r.DOMAIN_BE_RATES)


# Section 7.5: the run name, derived from the built switches and checked against the preview
def setup_run_name(r):
    # 7.5 Scenario summary, and the run name derived from it

    # The run-name suffix is derived from the switches, never typed

    # Tokens join in the run_names grammar's order; the management switches feed the one scenario slot
    # "-nofill" only when there was fill to withhold and a module to withhold it in
    _FILLS_WITHHELD = bool(not r.BN_SCHEDULE_APPLIED.projects and r.BN_SCHEDULE.projects
                           and r.ENABLE_BEACH_DUNE_MANAGEMENT)
    r.SCENARIO_SWITCHES = [
        ("period", f"{r.START_YEAR}-{r.END_YEAR} ({r.RUN_YEARS} yr)", None),
        ("scenario", r.SCENARIO, scenario_field(
            r.ENABLE_ROADWAY_MANAGEMENT, r.ENABLE_BEACH_DUNE_MANAGEMENT,
            _FILLS_WITHHELD, r.ENABLE_HISTORICAL_ROAD_RELOCATIONS)),
        ("roadway management", f"{sum(r.ROADWAY_MANAGEMENT_ON)} domains"
         if r.ENABLE_ROADWAY_MANAGEMENT else "off", None),
        ("historical relocations", r.ENABLE_HISTORICAL_ROAD_RELOCATIONS, None),
        ("beach/dune manager", f"{sum(r.BEACH_DUNE_MANAGEMENT_ON)} domains"
         if r.ENABLE_BEACH_DUNE_MANAGEMENT else "off", None),
        # Whether a period has fills follows from the period, so applied fills carry no token
        ("nourishment fills", f"{len(r.BN_SCHEDULE_APPLIED.projects)} applied"
         if r.ENABLE_NOURISHMENT_FILLS
         else ("suppressed" if r.BN_SCHEDULE.projects else "none in period"), None),
        ("groin", "on" if r.GROIN_ENABLED else "off", r.GROIN_TOKEN),
        ("source/sink preset", r.SOURCE_SINK_PRESET, r.SOURCE_SINK_PRESET),
        # Always metres since 2026-10-23; no token since 2026-10-29
        ("shoreline offset", "metres", None),
        # No token of its own: implied by the preset, and checked against it in 4.3
        ("background erosion", r.USE_BACKGROUND_EROSION, None),
        ("relocation target",
         "each domain's measured offset" if r.RELOCATION_SETBACK_M is None
         else f"{r.RELOCATION_SETBACK_M:g} m behind the dune line",
         r.RELOCATION_SETBACK_TOKEN),
        # Tokened only off the calibration value; last, so sweep_family reads the axis off it
        ("wave climate",
         f"Hs {r.RUN_CONFIG.hs} m, Tp {r.RUN_CONFIG.wave_period_s} s, "
         f"asym {r.RUN_CONFIG.wave_asymmetry}, "
         f"high-angle {r.RUN_CONFIG.wave_angle_high_fraction}"
         + ("" if r.WAVE_TOKEN is None else "  (off calibration)"),
         r.WAVE_TOKEN),
    ]

    r.RUN_NAME_SUFFIX = "_".join(
        token for _, _, token in r.SCENARIO_SWITCHES if token)
    r.RUN_NAME_BASE = f"{r.RUN_NAME_STEM}_{r.RUN_NAME_SUFFIX}"

    # A switch that did not reach its module is raised, not warned
    if r.RUN_NAME_BASE != r.RUN_NAME_PREVIEW:
        raise AssertionError(
            f"run name disagrees with the section 3 preview\n"
            f"  section 3    {r.RUN_NAME_PREVIEW}\n"
            f"  section 7.5  {r.RUN_NAME_BASE}\n"
            f"The preview follows the switches; this follows the built modules.")

    reports.scenario_summary_report(
        scenario=r.SCENARIO, departures=r.SCENARIO_DEPARTURES,
        switches=r.SCENARIO_SWITCHES, run_name_base=r.RUN_NAME_BASE,
        double_managed=r.DOUBLE_MANAGED_GIS, groin_enabled=r.GROIN_ENABLED,
        updrift_gis=r.GROIN_UPDRIFT_GIS,
        beach_dune_on_updrift=r.BEACH_DUNE_MANAGEMENT_ON[
            HATTERAS_DOMAINS.gis_to_pad(r.GROIN_UPDRIFT_GIS)])


# Section 8: the CoastSat target the run is scored against
def setup_coastsat_target(r):
    # 8. CoastSat target rates -- LOWESS windows

    # The CoastSat target: LOWESS at transect resolution, averaged to domains; GIS 1-10 unsmoothed (display only)

    # 8.1 Load both periods

    # Both periods always load: section 12 draws the other as reference

    r.COASTSAT_DATASETS = [
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
            label=f"CoastSat LRR (1996-{HATTERAS_PERIODS[1996]['end_year']})",
            period_start=1996,
            csv_path=str(COASTSAT_BASE_DIR / f"1996_{HATTERAS_PERIODS[1996]['end_year']}"
                         / "transect_lrr_full.csv"),
        ),
        CoastSatDataset(
            label=f"CoastSat LRR (2009-{HATTERAS_PERIODS[2009]['end_year']})",
            period_start=2009,
            csv_path=str(COASTSAT_BASE_DIR / f"2009_{HATTERAS_PERIODS[2009]['end_year']}"
                         / "transect_lrr_full.csv"),
        ),
    ]

    # 7-domain smoothing, the group's range (2026-09-28)
    r.LOWESS_CONFIG = LowessConfig(window_domains=(7,), skip_southern_domains=10)

    # The reference window, named: rate_comparison uses max(window_domains)
    r.TARGET_WINDOW = 7
    if r.TARGET_WINDOW != max(r.LOWESS_CONFIG.window_domains):
        raise ValueError(
            f"TARGET_WINDOW={r.TARGET_WINDOW} but rate_comparison will use "
            f"max(window_domains)={max(r.LOWESS_CONFIG.window_domains)} as the "
            f"reference curve -- section 12 would compare against a different "
            f"curve than the one reported here.")

    # An extended geometry appends its own CoastSat rows beyond GIS 90
    if HATTERAS_GEOMETRY_EXTENDED:
        r.COASTSAT_DATASETS = [
            CoastSatDataset(label=ds.label + " + extension",
                            period_start=ds.period_start,
                            csv_path=str(Path(ds.csv_path).parent / "ext"
                                         / "transect_lrr_with_base.csv"))
            if ds.period_start == r.START_YEAR else ds
            for ds in r.COASTSAT_DATASETS]

    r.cs_series = build_coastsat_series(
        r.COASTSAT_DATASETS, active_period_start=r.START_YEAR,
        lowess_config=r.LOWESS_CONFIG, domains=HATTERAS_DOMAINS)

    r.CS_ACTIVE = next((cs for cs in r.cs_series if cs["active"]), None)
    if r.CS_ACTIVE is None:
        raise RuntimeError(
            f"no CoastSat dataset has period_start == START_YEAR ({r.START_YEAR}); "
            f"loaded: {[cs['period_start'] for cs in r.cs_series]}")

    # 8.2 The target, as a table

    r.COASTSAT_TARGET = build_target_table(
        r.CS_ACTIVE, r.LOWESS_CONFIG, HATTERAS_DOMAINS, r.TARGET_WINDOW)

    # The surveyed reach's own target: what the interior score is graded against in every geometry
    r.COASTSAT_TARGET_BASE = r.COASTSAT_TARGET
    if HATTERAS_GEOMETRY_EXTENDED:
        _cs_base = build_coastsat_series(
            [CoastSatDataset(
                label=f"CoastSat LRR ({r.START_YEAR}-{r.END_YEAR}) surveyed",
                period_start=r.START_YEAR,
                csv_path=str(COASTSAT_BASE_DIR / f"{r.START_YEAR}_{r.END_YEAR}"
                             / "transect_lrr_full.csv"))],
            active_period_start=r.START_YEAR, lowess_config=r.LOWESS_CONFIG,
            domains=DEFAULT_DOMAINS)
        r.COASTSAT_TARGET_BASE = build_target_table(
            _cs_base[0], r.LOWESS_CONFIG, DEFAULT_DOMAINS, r.TARGET_WINDOW)

    # 8.4 Report

    reports.coastsat_report(
        target=r.COASTSAT_TARGET, active=r.CS_ACTIVE, target_window=r.TARGET_WINDOW,
        lowess_config=r.LOWESS_CONFIG, geometry=HATTERAS_DOMAINS, cs_series=r.cs_series,
        updrift_gis=r.GROIN_UPDRIFT_GIS, downdrift_gis=r.GROIN_DOWNDRIFT_GIS)


# Sections 9-10: rate conventions and the groin's paired no-groin baseline
def setup_conventions(r):
    # 9. Rate conventions, and the groin's paired baseline

    # The figures are drawn after the run from its folder (figure_making/model_output/run_figures.py)

    r.FLIP_SIGN_MODEL = True          # x_s_TS increases landward; flip so up = seaward

    # The estimator scored and recorded: "lrr" matches the CoastSat target; "endpoint" only reproduces old runs
    r.RATE_ESTIMATOR = "lrr"
    if r.RATE_ESTIMATOR not in ("lrr", "endpoint"):
        raise ValueError(f"RATE_ESTIMATOR must be 'lrr' or 'endpoint', "
                         f"got {r.RATE_ESTIMATOR!r}")

    # Reporting threshold only: nothing is dropped or flagged below it
    r.LRR_R2_FLOOR = 0.50

    # A groin run is measured against its no-groin twin (12.3), and its difference animation drawn from it
    r.GROIN_BASELINE_NAME = None
    r.GROIN_BASELINE_NPY = None
    if r.GROIN_ENABLED:
        r.GROIN_BASELINE_NAME = scenario_run_name(
            r.SCENARIO_SWITCHES, r.RUN_NAME_STEM, groin="nogroin")
        # A twin made before the 2026-10-29 rename sits under its old spelling
        for _twin in [r.GROIN_BASELINE_NAME] + other_spellings(r.GROIN_BASELINE_NAME):
            _baseline = resolve(r.OUTPUT_BASE_DIR / _twin, "matrix", _twin)
            if _baseline.exists():
                r.GROIN_BASELINE_NPY = str(_baseline)
                break

    # 10. build_cascade + run_cascade_simulation

    # Both live in cascade_pipeline/hindcast.py; run_years counts transitions, not states

    print("\nbuild_cascade + run_cascade_simulation defined")
    print(f"  run_years -> {r.RUN_YEARS} transitions, "
          f"time_step_count={r.RUN_YEARS + 1}, {r.RUN_YEARS + 1} annual states "
          f"({r.FIRST_MODEL_YEAR}-{r.LAST_MODEL_YEAR}, final state 1 Jan {r.LAST_MODEL_YEAR + 1})")
    print("  nourishment via BN_SCHEDULE.apply_to_cascade -> "
          "cascade.nourishment_volume")
    print("  road events via cascade_pipeline.roadway.apply_historical_event")


# Section 11.1: the model parameters sections 2-9 do not produce
def setup_parameters(r):
    # 11. Initialize CASCADE -- single config, no sweep

    # 11.1 Parameters sections 2-9 do not produce

    r.NUM_CORES = 1        # >1 has crashed on this configuration; leave at 1

    # Dune growth (Barrier3D logistic growth bounds). 0.01-0.03 since 2026-10-20, paired with
    # Dmax 10.5 m above the berm (Dmaxel 12.2 m, parameter file): on the 1996 start with the ALACE
    # offset corrected (dune-topo v3) its mean, 0.02, reproduces the 2009 crest of the stable dunes
    # (GIS 59-62, 20-22), and the range doubles the cell-to-cell spread of a fixed rate
    # (data/.../1996-start/dune-growth-rate/README.md); was 0.55-0.95
    # Per period since 2026-10-21 (Hannah): 2009-2025 takes 0.02-0.04, the rate that tracks the 2014-2019
    # lidar's median crest within 0.16 m (0.01-0.03 falls 0.14-0.39 m behind); fitted to the test window's own
    # surveys, so that window tests overwash and shoreline, not dune growth
    # (data/.../2009-start/dune-growth-rate/README.md). 1984/2004 keep the calibration rate.
    r.DUNE_GROWTH_RATES = {1996: (0.01, 0.03), 2009: (0.02, 0.04)}
    _RLO, _RHI = r.DUNE_GROWTH_RATES.get(r.START_YEAR, r.DUNE_GROWTH_RATES[1996])
    r.RMIN = [_RLO] * HATTERAS_DOMAINS.total_domains
    r.RMAX = [_RHI] * HATTERAS_DOMAINS.total_domains
    # Cape Point no longer gets its own rate (Hannah, 2026-10-21): 0.10-0.20 at GIS 1-5 matched the 2014-2017
    # cape surveys but built the 1996-2009 cape dunes 3-4 m above the surveys and barely moved GIS 1-2's
    # shoreline (dunes/2026-10-21-cape-growth-rate); one rate per period island-wide, the cape's missing
    # sand handled by the source/sink instead.
    # HAT_GROWTH_OVERRIDE="gis=rmin:rmax,...": per-domain bounds for an experiment without editing this file
    # (2026-10-21, the alongshore growth sensitivity test); unset, every domain keeps the period's rate
    r.GROWTH_OVERRIDE_RAW = os.environ.get("HAT_GROWTH_OVERRIDE", "").strip()
    for _pair in filter(None, r.GROWTH_OVERRIDE_RAW.split(",")):
        _g, _bounds = _pair.split("=")
        _lo, _hi = (float(_v) for _v in _bounds.split(":"))
        if not 0 <= _lo <= _hi:
            raise SystemExit(f"\nHAT_GROWTH_OVERRIDE entry {_pair!r}: need 0 <= rmin <= rmax.\n")
        r.RMIN[HATTERAS_DOMAINS.gis_to_pad(int(_g))] = _lo
        r.RMAX[HATTERAS_DOMAINS.gis_to_pad(int(_g))] = _hi
    if r.GROWTH_OVERRIDE_RAW:
        print(f"\nHAT_GROWTH_OVERRIDE   {len(r.GROWTH_OVERRIDE_RAW.split(','))} domain(s) given their own growth bounds")

    # Dune rebuild thresholds, m MHW

    # In metres, the documented unit
    # Rebuild target: the NC-12 safe dune crest, 4.3 m NAVD88 (Velasquez-Montoya et al. 2021,
    # Natural Hazards Review 22(2)), less MHW 0.36 m. Was 3.0, unsourced, until 2026-10-13
    r.DUNE_DESIGN_ELEVATION_M = 3.94   # rebuild target
    r.DUNE_MINIMUM_ELEVATION_M = 0.0   # rebuild trigger: CASCADE's berm floor governs, explicitly
    r.DUNE_DESIGN_ELEVATION = [r.DUNE_DESIGN_ELEVATION_M] * HATTERAS_DOMAINS.total_domains
    r.DUNE_MINIMUM_ELEVATION = [r.DUNE_MINIMUM_ELEVATION_M] * HATTERAS_DOMAINS.total_domains

    # Roadway
    r.ROAD_ELEVATION = r.road_elevation_full   # per-domain, m MHW, from section 5
    r.ROAD_WIDTH = 20.0

    # Wave climate: one configuration, no sweep

    # From hat_run.yaml's physics block, so a sweep never edits this file
    r.Hs = r.RUN_CONFIG.hs                                       # m, calibration 2.5
    r.FIXED_WAVE_PERIOD = r.RUN_CONFIG.wave_period_s             # s, calibration 8
    r.FIXED_WAVE_ASYMMETRY = r.RUN_CONFIG.wave_asymmetry         # calibration 0.7
    r.FIXED_WAVE_ANGLE_HIGH_FRACTION = r.RUN_CONFIG.wave_angle_high_fraction  # 0.1

    # The run was named for these in 7.5: they must still match
    assert r.WAVE_VALUES == {
        "hs": r.Hs, "wave_period_s": r.FIXED_WAVE_PERIOD,
        "wave_asymmetry": r.FIXED_WAVE_ASYMMETRY,
        "wave_angle_high_fraction": r.FIXED_WAVE_ANGLE_HIGH_FRACTION,
    }, "section 11's wave climate is not the one section 3 named the run for"

    # Datums
    r.BERM_ELEVATION = 1.7    # m NAVD88, Hatteras Island, NCDOT-derived via NC State
    r.MHW_ELEVATION = 0.36    # m NAVD88, Duck NC gauge (NOAA 8651370)
    r.BEACH_SLOPE = 0.06      # Barrier3D beta; the storm series' run-up slope (since 2026-10-20)

    # Sandbags: off for the hindcast

    # From hat_run.yaml (top-level sandbags); reaches the run metadata as sandbags_on
    r.ENABLE_SANDBAG_PLACEMENT = r.RUN_CONFIG.sandbags
    r.SANDBAG_MANAGEMENT_ON = [r.ENABLE_SANDBAG_PLACEMENT] * HATTERAS_DOMAINS.total_domains
    r.SANDBAG_ELEVATION = 0

    r.SEA_LEVEL_CONSTANT = True


STEPS = (setup_topography, setup_scenario, setup_island_offset, setup_forcings, setup_roadway, setup_beach_dune, setup_groin, setup_run_name, setup_coastsat_target, setup_conventions, setup_parameters)


# Every step in order, for one run's settings
def prepare_run(config):
    r = RunSetup(RUN_CONFIG=config)
    for step in STEPS:
        step(r)
    return r
