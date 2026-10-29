#!/usr/bin/env python3
"""
The shortest Hatteras CASCADE run: load the inputs, build the model, run it, save the shoreline.

    python scripts/real-world-cases/Hatteras/hatteras_ms/HAT_run_simple.py

Every input file is named below, and every value the model gets is written out. Edit the settings
and run. No checks, reports or figures: for those use HAT_hindcast.py, which gives the model the
same inputs. The historical NC-12 relocations are not applied (they are off in every hindcast scenario).

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-10-23
"""

import os
import shutil
import sys
from pathlib import Path

import numpy as np

# --- SETTINGS ----------------------------------------------------------------
START_YEAR = 1996               # 1996 (storms 1997-2009) or 2009 (storms 2010-2025)
ROAD_MANAGEMENT = True          # NC-12 defended
BEACH_DUNE_MANAGEMENT = True    # the villages managed
NOURISHMENT = True              # the historical beach fills (needs BEACH_DUNE_MANAGEMENT)
GROIN = False                   # the Buxton groin
SAVE_MODEL_STATE = True         # the ~160 MB .npz of the whole model

# Source/sink rates, m/yr by GIS domain ((-) erosion, (+) accretion); {} is zeroBE.
# edgeBE is {1: 3.9368, 90: 1.2077} for 1996 and {1: 64.2148, 90: 7.547} for 2009;
# domainBE has a rate for every domain (HATTERAS_BE_RATES in site_layer/hatteras_site_config.py)
BACKGROUND_EROSION = {}
# -----------------------------------------------------------------------------

# --- INPUT FILES ---------------------------------------------------------------
WORKSPACE = Path(r"D:\hatteras_ch1\model_runs")
DATA = WORKSPACE / "data" / "hatteras_init"
DOMAINS = DATA / "1-barrier3d-domains"
OFFSETS = DATA / "2-brie-offset"
STORMS = DATA / "3-env-forcings" / "3-storms" / "hindcast_storms"
ROADS = DATA / "4-mgmt-forcing"

PERIODS = {
    1996: dict(
        first_year=1997, last_year=2009,                       # the storm years
        topography=DOMAINS / "1996-start" / "dune-topo" / "v3",
        island_offset=OFFSETS / "1996" / "shoreline" / "v4" / "Island_Shoreline_Offsets_1996_PADDED_120.csv",
        storms=STORMS / "1997_2009" / "1997_2009_storms_v3_240.npy",
        road_setback=ROADS / "road_offset" / "dunestart_offset" / "derived" / "1996" / "RoadSetback_1996_dunestart.csv",
        sea_level_rise=0.002,                                   # m/yr
        dune_growth=(0.01, 0.03),                               # Barrier3D rmin, rmax
    ),
    2009: dict(
        first_year=2010, last_year=2025,
        topography=DOMAINS / "2009-start" / "dune-topo" / "v1",
        island_offset=OFFSETS / "2009" / "shoreline" / "v4" / "Island_Shoreline_Offsets_2009_PADDED_120.csv",
        storms=STORMS / "2010_2025" / "2010_2025_storms_v3_240.npy",
        road_setback=ROADS / "road_offset" / "dunestart_offset" / "derived" / "2009" / "RoadSetback_2009_dunestart.csv",
        sea_level_rise=0.006,
        dune_growth=(0.02, 0.04),
    ),
}
BUFFER_TOPOGRAPHY = DOMAINS / "shared" / "buffer" / "sample_1_topography.npy"
BUFFER_DUNES = DOMAINS / "shared" / "buffer" / "sample_1_dune.npy"
ROAD_ELEVATION = ROADS / "road_elevation" / "RoadElevation.csv"
PARAMETERS = DATA / "Hatteras-CASCADE-parameters.yaml"
# -----------------------------------------------------------------------------

# --- THE ISLAND ------------------------------------------------------------------
# 90 real domains (GIS 1 at Cape Point to GIS 90 at Pea Island), 500 m each, with 30 buffer domains
# on each end: 150 in all. Domain GIS g sits at position g - 1 + 30 in every list below.
N_REAL, N_BUFFER, DOMAIN_WIDTH_M = 90, 30, 500.0
N = N_REAL + 2 * N_BUFFER


def pad(gis):
    return gis - 1 + N_BUFFER


VILLAGES = [(7, 8), (21, 31), (68, 83)]       # GIS ranges: Buxton, Avon, Rodanthe-Waves-Salvo
ROAD_GIS = range(9, 91)                        # where NC-12 runs
FILLS = [                                      # (year, GIS domains, cubic yards)
    (2014, range(82, 89), 1_620_000),          # Rodanthe emergency fill
    (2017, range(6, 17), 2_600_000),           # Buxton beach nourishment
    (2022, range(6, 17), 1_200_000),           # Buxton shore protection
    (2022, range(21, 29), 1_000_000),          # Avon shore protection
]
BRIDGE = (2022, range(82, 89))                 # Jug Handle Bridge: NC-12 off the island, no longer managed
# -----------------------------------------------------------------------------

# Fixes to CASCADE the lab has not merged yet, applied before the model is built
HATTERAS = next(p for p in Path(__file__).resolve().parents if (p / "hatteras_paths.py").exists())
sys.path.insert(0, str(HATTERAS))
import cascade_patch  # noqa: E402
cascade_patch.apply()
from cascade.cascade import Cascade  # noqa: E402
from cascade_pipeline.hindcast import build_island_offset  # noqa: E402
from cascade_pipeline.domains import DomainGeometry  # noqa: E402
from cascade_pipeline.run_names import scenario_field  # noqa: E402

period = PERIODS[START_YEAR]
first_year, last_year = period["first_year"], period["last_year"]
years = last_year - first_year + 1

# The hindcast's grammar (cascade_pipeline/run_names.py) with its own site word and no preset
name = (f"simple_{START_YEAR}-{last_year}_{scenario_field(ROAD_MANAGEMENT, BEACH_DUNE_MANAGEMENT)}"
        f"_{'groin' if GROIN else 'nogroin'}")
run_dir = WORKSPACE / "output" / "raw_runs" / "simple" / name
run_dir.mkdir(parents=True, exist_ok=True)
os.chdir(WORKSPACE)

# --- 1. Load the inputs --------------------------------------------------------
topography = period["topography"]
elevation_files = ([str(BUFFER_TOPOGRAPHY)] * N_BUFFER
                   + [str(topography / "topography" / f"domain_{g}_topography.npy") for g in range(1, N_REAL + 1)]
                   + [str(BUFFER_TOPOGRAPHY)] * N_BUFFER)
dune_files = ([str(BUFFER_DUNES)] * N_BUFFER
              + [str(topography / "dunes" / f"domain_{g}_dune.npy") for g in range(1, N_REAL + 1)]
              + [str(BUFFER_DUNES)] * N_BUFFER)

# Where each domain's shoreline starts, in metres. The file's buffers are replaced with a smooth curve
# that joins the two ends of the island, because BRIE treats the coast as a loop.
island_offset = build_island_offset(
    period["island_offset"],
    DomainGeometry(num_real_domains=N_REAL, num_buffer_domains=N_BUFFER, first_gis_id=1,
                   domain_spacing_m=DOMAIN_WIDTH_M))

background_erosion = [0.0] * N
for gis, rate in BACKGROUND_EROSION.items():
    background_erosion[pad(gis)] = rate


def by_domain(path, fill):
    """A 2-row file (GIS domains, values) as one value per domain; `fill` where it has none."""
    gis_ids, values = np.loadtxt(path, delimiter=",")
    out = np.full(N, fill)
    for gis, value in zip(gis_ids.astype(int), values):
        if gis in ROAD_GIS:
            out[pad(gis)] = value
    return out


road_setback = by_domain(period["road_setback"], fill=0.0)          # m behind the dune line
road_elevation = by_domain(ROAD_ELEVATION, fill=1.45)               # m above MHW

in_village = [False] * N
for first, last in VILLAGES:
    for gis in range(first, last + 1):
        in_village[pad(gis)] = True
road_on = [ROAD_MANAGEMENT and (i - N_BUFFER + 1) in ROAD_GIS and not in_village[i] for i in range(N)]

# The fills in this period, as m^3 per m of beach in each domain
fills = {}
for year, gis_domains, cubic_yards in FILLS:
    if NOURISHMENT and BEACH_DUNE_MANAGEMENT and first_year <= year <= last_year:
        per_m = cubic_yards * 0.764555 / len(gis_domains) / DOMAIN_WIDTH_M
        for gis in gis_domains:
            fills.setdefault(year, {})
            fills[year][pad(gis)] = fills[year].get(pad(gis), 0.0) + per_m

# The beach/dune manager runs in the villages and wherever a fill lands, except where the road manager runs
filled = {pad(gis) for year, gis_domains, _ in FILLS if first_year <= year <= last_year for gis in gis_domains}
beach_dune_on = [BEACH_DUNE_MANAGEMENT and (in_village[i] or i in filled) and not road_on[i] for i in range(N)]
overwash_filter = [40.0 if (BEACH_DUNE_MANAGEMENT and in_village[i]) else 0.0 for i in range(N)]   # % removed

# Each run gets its own copy of the parameter file, because CASCADE edits it while it builds
parameter_file = run_dir / f"{name}-parameters.yaml"
shutil.copyfile(PARAMETERS, parameter_file)

# --- 2. Build the model --------------------------------------------------------
rmin, rmax = period["dune_growth"]

cascade = Cascade(
    DATA, name,
    storm_file=str(period["storms"]),
    elevation_file=elevation_files,
    dune_file=dune_files,
    parameter_file=str(parameter_file),
    time_step_count=years + 1,
    alongshore_section_count=N,
    num_cores=1,

    # the island
    enable_shoreline_offset=True,
    shoreline_offset=island_offset,
    min_dune_growth_rate=[rmin] * N,
    max_dune_growth_rate=[rmax] * N,
    berm_elevation=1.7,                 # m NAVD88
    MHW=0.36,                           # m NAVD88 (Duck gauge)
    beta=0.06,                          # beach slope

    # the forcing
    wave_height=2.0,                    # m
    wave_period=7.5,                    # s
    wave_asymmetry=0.6,
    wave_angle_high_fraction=0.5,
    sea_level_rise_rate=period["sea_level_rise"],
    sea_level_rise_constant=True,
    background_erosion=background_erosion,
    alongshore_transport_module=True,
    outwash_module=False,
    community_economics_module=False,

    # NC-12
    roadway_management_module=road_on,
    road_setback=road_setback,
    road_ele=road_elevation,
    road_width=20.0,                    # m
    dune_design_elevation=[3.94] * N,   # m above MHW: the NC-12 safe dune crest
    dune_minimum_elevation=[0.0] * N,
    road_relocation_setback=20.0,       # m behind the dune, where a relocated road is rebuilt
    dune_mgmt_road_setback_trigger=20.0,

    # the villages
    beach_nourishment_module=beach_dune_on,
    overwash_filter=overwash_filter,
    overwash_to_dune=9.0,               # % of the removed overwash put back on the dune
    nourishment_volume=[0.0] * N,       # set each year from the fills below
    nourishment_interval=None,
    sandbag_management_on=[False] * N,
    sandbag_elevation=0,

    # the Buxton groin, between GIS 5 and 6: holds back 10% of the drift from 1969 until the 2017 fill
    groin_module=GROIN,
    groin_updrift_pad=pad(6),
    groin_downdrift_pad=pad(5),
    groin_strength_history=[(1969, 0.10), (2017, 0.0)],
    groin_start_year=first_year,
    groin_spread_weights=(0.0, 2 / 3, 1 / 3),
)

# --- 3. Run it, one year at a time ---------------------------------------------
for step in range(years):
    year = first_year + step

    # This year's beach fills (none in most years)
    this_year = fills.get(year, {})
    cascade.nourish_now = np.array([1.0 if i in this_year else 0.0 for i in range(N)])
    cascade.nourishment_volume = [this_year.get(i, 0.0) for i in range(N)]

    # The Jug Handle Bridge: the road there is no longer managed
    if year == BRIDGE[0]:
        flags = list(cascade.roadway_management_module)
        for gis in BRIDGE[1]:
            if cascade.roadways[pad(gis)] is not None:
                flags[pad(gis)] = False
        cascade.roadway_management_module = flags

    cascade.update()

    # A domain whose road drowned is handed to the beach/dune manager
    road_flags = list(cascade.roadway_management_module)
    beach_flags = list(cascade.beach_nourishment_module)
    drowned = [i for i in range(N)
               if road_flags[i] and cascade.road_break[i] and cascade.roadways[i].drown_break]
    if drowned:
        for i in drowned:
            road_flags[i], beach_flags[i] = False, True
        cascade.roadway_management_module = road_flags
        cascade.beach_nourishment_module = beach_flags

    print(f"year {step + 1}/{years} ({year})")
    if getattr(cascade, "b3d_break", False):
        print("the model stopped early (a barrier drowned)")
        break

# --- 4. Save -------------------------------------------------------------------
shoreline_m = np.array([b3d.x_s_TS for b3d in cascade.barrier3d]).T * 10.0     # year x domain, metres
np.save(run_dir / f"{name}_shoreline_matrix.npy", shoreline_m)
if SAVE_MODEL_STATE:
    cascade.save(str(run_dir))
print(f"done  {run_dir}")
