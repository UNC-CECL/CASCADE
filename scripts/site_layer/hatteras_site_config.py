"""
Hatteras Island site config: the real domains, periods, forcing files, presets and management events.

    from site_layer.hatteras_site_config import HATTERAS_DOMAINS, HATTERAS_PERIODS, HATTERAS_ANNOTATIONS

Fills cascade_pipeline's generic dataclasses with Hatteras content; a new site writes a sibling
module in this shape. The calibrated source/sink table is edited as text by 7-source-sink/2-calibrate. Details: scripts/site_layer/README.md.

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-09-30
"""

import csv
import os
import re

from cascade_pipeline.annotations import AnnotationConfig
from cascade_pipeline.domains import DomainGeometry
from cascade_pipeline.nourishment import BeachDuneConfig, NourishmentProject
from cascade_pipeline.roadway import (
    BridgeEvent, RelocationEvent, load_road_setbacks)

# hat_topo_version owns the data root and the year -> input pairings (topography, road line, setback file)
from site_layer.hat_topo_version import (INIT_ROOT, ROAD_LINE_FOR_YEAR,  # noqa: F401
                              YEAR_PRODUCT, road_setback_relpath)
# The other 4-mgmt-forcing paths (road elevation, relocation) come from hat_topo_version too
from site_layer import hat_topo_version as _tv_mgmt  # noqa: E402
# The storm series, from 3-env-forcings (2026-09-18).
from site_layer import hat_env_forcings as _env  # noqa: E402
from site_layer.hat_extension_domains import (BASE_GEOMETRY, gis_bounds,  # noqa: F401
                                   geometry_label, is_extended)

# GIS 1-90 plus 15 buffers a side; HAT_GEOMETRY picks an extended reach (unset = base)
HATTERAS_GEOMETRY = (os.environ.get("HAT_GEOMETRY", "").strip() or BASE_GEOMETRY)

# Island offset source for this run: "shoreline" (default since 2026-10-05) or "duneline"; off-default is not a matrix run
HATTERAS_OFFSET_SOURCE = (os.environ.get("HAT_ISLAND_OFFSET_SOURCE", "").strip()
                          or _tv_mgmt.RUN_OFFSET_SOURCE)
if HATTERAS_OFFSET_SOURCE not in _tv_mgmt.OFFSET_SOURCES:
    raise SystemExit(
        f"HAT_ISLAND_OFFSET_SOURCE={HATTERAS_OFFSET_SOURCE!r} is not a known "
        f"island-offset source; known: {', '.join(_tv_mgmt.OFFSET_SOURCES)}")
HATTERAS_GEOMETRY_EXTENDED = is_extended(HATTERAS_GEOMETRY)
_FIRST_GIS, _LAST_GIS = gis_bounds(HATTERAS_GEOMETRY)
HATTERAS_DOMAINS = DomainGeometry(
    num_real_domains=_LAST_GIS - _FIRST_GIS + 1,
    num_buffer_domains=15,
    first_gis_id=_FIRST_GIS,
    domain_spacing_m=500.0,
)

# The interior score is GIS 2-89 against the surveyed CoastSat table in every geometry
SCORE_INTERIOR_GIS = (2, 89)

HATTERAS_ANNOTATIONS = AnnotationConfig(
    town_spans={
        "Buxton": (7, 8),
        "Avon": (21, 31),
        "Tri-Village": (68, 83),  # Salvo / Waves / Rodanthe
    },
    village_lines={
        "Salvo": 69, "Waves": 74, "Rodanthe": 80,
    },
    piers={
        "Avon Pier": (26, 0.76),
        "Rodanthe Pier": (79, 0.76),
    },
    groins={
        "Buxton Groin": 5.5,  # boundary between domains 5 and 6
    },
    shoal_zones={
        "Avon Shoals": (24, 39),
        "Wimble Shoals": (60, 74),
    },
    region_name="Hatteras Island",
    low_end_label="S | Cape Point",
    high_end_label="Pea Island | N",
    obs_source_name="CoastSat",
)


# Period forcing config

# Periods: the start year resolves every period-dependent forcing; paths are relative to data/hatteras_init

# "topo_product" is the period's Barrier3D domain folder (from YEAR_PRODUCT); offset builds resolved below
def _island_offset_file(start_year):
    """Path of the padded offset file, relative to INIT_ROOT.

    120 domains in the base geometry. An extended geometry reads its own
    build under <year>/ext/<geometry>/ (island_offset_hybrid.py --geometry),
    which is not a version and does not move CURRENT; the path is returned
    unchecked because every period resolves here at import and only the
    period being run needs the file to exist -- the runner checks that.
    """
    # Every part of this path comes from hat_topo_version
    source = HATTERAS_OFFSET_SOURCE
    base = _tv_mgmt.offset_start_dir(start_year, source)
    if HATTERAS_GEOMETRY_EXTENDED:
        name = _tv_mgmt.offset_basename(start_year, source)
        return _tv_mgmt.init_relpath(
            base / "ext" / HATTERAS_GEOMETRY
            / f"{name}_PADDED_{HATTERAS_DOMAINS.total_domains}.csv")
    # The version choice lives in hat_topo_version.offset_version
    return _tv_mgmt.init_relpath(_tv_mgmt.offset_file(start_year, "padded", 120,
                                                      source=source))


def island_offset_version(start_year):
    """The version segment of the offset file this period resolves to.

    "duneline/v2" since 2026-09-22, when every build moved under the SOURCE it
    was measured from; "v2" before that, and "flat" for the unversioned layout
    older still. Recorded in run metadata and run_index.csv (2026-09-15)
    because a v1 run and a v2 run are otherwise identical on disk: the run name
    carries no offset token and the file name is the same in every version
    folder. Resolves through _island_offset_file so it can never disagree with
    the file that was read.

    READING AN OLD RUN: a token with no "/" is from before the sources were
    split, and every build then was dune-derived -- "v1" means "duneline/v1".
    """
    if HATTERAS_GEOMETRY_EXTENDED:
        return f"ext/{HATTERAS_GEOMETRY}"
    # <year>/<source>/<version>/<file>, relative to 2-brie-offset.
    parts = _island_offset_file(start_year).split("/")
    source, version = parts[2], parts[3]
    return f"{source}/{version}" if re.fullmatch(r"v\d+", version) else source


# MODEL YEARS (2026-10-04). "end_year" is the window LABEL: it names the storm file, the CoastSat
# target folder and the run folders (1996_2009, 2009_2025). "last_model_year" is the last calendar
# year the model steps through. A run covers start_year .. last_model_year inclusive, one transition
# per calendar year, so run_years(start) = last_model_year - start_year + 1, and the final saved
# state is 1 January of last_model_year + 1. Full explanation: scripts/hatteras_ms/MODEL_YEARS.md
HATTERAS_PERIODS = {
    1984: {
        "end_year": 2004,
        "last_model_year": 2003,
        # 0.00391 m/yr fitted over 1984-2004 (Duck gauge), stored to 0.001
        "sea_level_rise_rate": 0.004,
        "storm_file": _env.init_relpath(_env.storm_series_file(1984, 2004)),
        "island_offset_file": _island_offset_file(1984),
        # Measured on the 1978 NC-12 line against row 0 of 1984-start; paired with the topography version
        "road_setback_file": road_setback_relpath(1984),
        "topo_product": YEAR_PRODUCT[1984],
        "enable_nourishment": False,
        "nourishment_volume": 0,  # m^3/m
    },
    2004: {
        "end_year": 2024,
        "last_model_year": 2023,
        # 0.00639 m/yr fitted over 2004-2024; see rslr/fits/duck_rslr_rates.csv.
        "sea_level_rise_rate": 0.006,
        "storm_file": _env.init_relpath(_env.storm_series_file(2004, 2024)),
        "island_offset_file": _island_offset_file(2004),
        # Measured on the 2008 NC-12 line against row 0 of 2004-start
        "road_setback_file": road_setback_relpath(2004),
        "topo_product": YEAR_PRODUCT[2004],
        "enable_nourishment": True,  # historical BN injected per-year in the time loop
        "nourishment_volume": 100,  # m^3/m passed to Cascade init
    },

    # DEM to DEM (2026-10-05): calibrate 1996 -> 2009, test 2009 -> 2025; see MODEL YEARS above
    1996: {
        # 1996-2009 since 2026-10-05: ALACE 1996 to the USACE 2009 DEM (was 1996-2015, before it 1996-2010)
        "end_year": 2009,
        # Exclusive label: 13 years, 1996-2008, ending 1 Jan 2009; the DEM was flown 2009-08-17
        "last_model_year": 2008,
        # 0.00252 m/yr fitted over 1996-2009 (Duck gauge; CI +/-0.0025). 1996-2015 was 0.00376
        "sea_level_rise_rate": 0.003,
        "storm_file": _env.init_relpath(_env.storm_series_file(1996, 2009)),
        # Derived, not surveyed: built from the 1997 dune line
        "island_offset_file": _island_offset_file(1996),
        # Derived: the 1984 setbacks with the 1989 Pea Island relocation applied
        "road_setback_file": road_setback_relpath(1996),
        "topo_product": YEAR_PRODUCT[1996],
        # No fill falls in 1996-2008; Rodanthe 2014 is the first
        "enable_nourishment": False,
        "nourishment_volume": 100,  # m^3/m passed to Cascade init
    },
    2009: {
        # Key was 2010 until 2026-10-05: the period starts in its DEM's year (USACE 2009, flown 2009-08-10 to 08-24)
        "end_year": 2025,
        # Exclusive label: 16 years, 2009-2024, ending 1 Jan 2025; the target is 2025-08-17 +/- 6 months
        "last_model_year": 2024,
        # 0.00509 m/yr fitted over 2009-2025 (Duck gauge, record through 2025)
        "sea_level_rise_rate": 0.005,
        "storm_file": _env.init_relpath(_env.storm_series_file(2009, 2025)),
        # Derived, not surveyed: the 2009 dune line, or the CoastSat mean over 2008-08-17 to 2010-08-17
        "island_offset_file": _island_offset_file(2009),
        # A copy of the 2004 file: same topography, same road line, no relocation between
        "road_setback_file": road_setback_relpath(2009),
        "topo_product": YEAR_PRODUCT[2009],
        # Rodanthe 2014, Buxton 2017 and both 2022 projects fall inside 2009-2024
        "enable_nourishment": True,
        "nourishment_volume": 100,  # m^3/m passed to Cascade init
    },
}

# The last calendar year a period's run steps through; raises if end_year was changed without it
def last_model_year(start_year):
    p = HATTERAS_PERIODS[start_year]
    last = p["last_model_year"]
    if last not in (p["end_year"] - 1, p["end_year"]):
        raise ValueError(
            f"period {start_year}: last_model_year {last} does not sit at end_year {p['end_year']} "
            f"or the year before; set both when changing a window (scripts/hatteras_ms/MODEL_YEARS.md)")
    return last


# Annual transitions a period's run makes: start_year .. last_model_year inclusive
def run_years(start_year):
    return last_model_year(start_year) - start_year + 1


# Background erosion (source/sink) rates, m/yr by GIS domain

# Source/sink presets, m/yr, sparse (absent = 0.0), (-) erosion: zeroBE, edgeBE, calibBE

# "zeroBE": empty, since an absent domain is already 0.0
HATTERAS_BE_RATES_ZERO = {
    1984: {},
    2004: {},
    1996: {},
    2009: {},
}

# The edge-only domains: the first and last real domains (buffers stay 0.0)
HATTERAS_BE_EDGE_DOMAINS = (HATTERAS_DOMAINS.first_gis_id,
                            HATTERAS_DOMAINS.last_gis_id)

# "calibBE": the per-domain fit against the CoastSat LRR; edgeBE slices its GIS 1 from it
# Refit 2026-08-28 — both periods, five passes, on the current topography

# Every value below re-solved 2026-08-28 on the current topography (method and RMSE in the README)
HATTERAS_BE_RATES_CALIBRATED = {
    1984: {
          1: -42.6,  # LOCKED — end domain, LRR-solved; see the end-domain note above
          2: +0.0,  # Cape Point / Shoal Dynamics
          3: +0.0,  # Cape Point / Shoal Dynamics
          4: +0.0,  # Cape Point / Shoal Dynamics
          5: +0.0,  # Cape Point / Shoal Dynamics
          6: +0.0,  # Cape Point / Shoal Dynamics
          7: +0.0,  # Cape Point / Shoal Dynamics
          8: +3.7,  # Cape Point / Shoal Dynamics
          9: +0.0,  # Cape Point / Shoal Dynamics
         10: -3.0,  # Cape Point / Shoal Dynamics
         11: -1.4,  # Buxton–Avon Transition
         12: -1.4,  # Buxton–Avon Transition
         13: -0.7,  # Buxton–Avon Transition
         14: +0.0,  # Buxton–Avon Transition
         15: +0.0,  # Buxton–Avon Transition
         16: +0.0,  # Buxton–Avon Transition
         17: +0.0,  # Buxton–Avon Transition
         18: +0.0,  # Buxton–Avon Transition
         19: +0.0,  # Buxton–Avon Transition
         20: +0.0,  # Buxton–Avon Transition
         21: +0.0,  # Avon
         22: +0.0,  # Avon
         23: +0.0,  # Avon
         24: +0.0,  # Avon
         25: +0.0,  # Avon
         26: +0.0,  # Avon
         27: +0.9,  # Avon
         28: +1.2,  # Avon
         29: +2.4,  # Avon
         30: +3.3,  # Avon
         31: +4.1,  # Avon
         32: +3.6,  # Mid-island
         33: +2.0,  # Mid-island
         34: +0.9,  # Mid-island
         35: +0.0,  # Mid-island
         36: +0.0,  # Mid-island
         37: +0.0,  # Mid-island
         38: +0.0,  # Mid-island
         39: +0.0,  # Mid-island
         40: +0.0,  # Mid-island
         41: +0.0,  # Mid-island
         42: +0.0,  # Mid-island
         43: +0.0,  # Mid-island
         44: +0.0,  # Mid-island
         45: +0.0,  # Mid-island
         46: +0.0,  # Mid-island
         47: +0.0,  # Mid-island
         48: -1.1,  # Mid-island
         49: -1.4,  # Mid-island
         50: -2.1,  # Mid-island
         51: -2.4,  # Mid-island
         52: -2.2,  # Mid-island
         53: -2.0,  # Mid-island
         54: -1.8,  # Mid-island
         55: -1.0,  # Mid-island
         56: -1.0,  # Mid-island
         57: -0.3,  # Mid-island
         58: +0.0,  # Mid-island
         59: +0.0,  # Mid-island
         60: +0.0,  # Wimble Shoals Influence
         61: +0.0,  # Wimble Shoals Influence
         62: +0.0,  # Wimble Shoals Influence
         63: +0.0,  # Wimble Shoals Influence
         64: +0.0,  # Wimble Shoals Influence
         65: +0.0,  # Wimble Shoals Influence
         66: +0.0,  # Wimble Shoals Influence
         67: +0.0,  # Wimble Shoals Influence
         68: +0.5,  # Wimble Shoals Influence
         69: +1.6,  # Wimble Shoals Influence
         70: +2.2,  # Wimble Shoals Influence
         71: +3.4,  # Wimble Shoals Influence
         72: +4.0,  # Wimble Shoals Influence
         73: +4.2,  # Wimble Shoals Influence
         74: +3.8,  # Wimble Shoals Influence
         75: +1.9,  # Tri-Village / Rodanthe
         76: +0.0,  # Tri-Village / Rodanthe
         77: +0.0,  # Tri-Village / Rodanthe
         78: -1.2,  # Tri-Village / Rodanthe
         79: -2.1,  # Tri-Village / Rodanthe
         80: -3.5,  # Tri-Village / Rodanthe
         81: -3.8,  # Tri-Village / Rodanthe
         82: -4.2,  # Tri-Village / Rodanthe
         83: -4.8,  # Tri-Village / Rodanthe
         84: -5.3,  # Pea Island NWR
         85: -5.2,  # Pea Island NWR
         86: -4.8,  # Pea Island NWR
         87: -3.9,  # Pea Island NWR
         88: -3.0,  # Pea Island NWR
         89: -1.8,  # Pea Island NWR
         90: +32.8,  # LOCKED — end domain, LRR-solved; see the end-domain note above
    },
    2004: {
          1: +50.3,  # LOCKED — end domain, LRR-solved; see the end-domain note above
          2: +0.0,  # Cape Point / Shoal Dynamics
          3: +0.0,  # Cape Point / Shoal Dynamics
          4: +0.0,  # Cape Point / Shoal Dynamics
          5: +0.0,  # Cape Point / Shoal Dynamics
          6: +0.0,  # Cape Point / Shoal Dynamics
          7: +0.0,  # Cape Point / Shoal Dynamics
          8: +0.0,  # Cape Point / Shoal Dynamics
          9: +0.3,  # Cape Point / Shoal Dynamics
         10: +1.5,  # Cape Point / Shoal Dynamics
         11: +1.2,  # Buxton–Avon Transition
         12: +1.5,  # Buxton–Avon Transition
         13: +1.9,  # Buxton–Avon Transition
         14: +2.0,  # Buxton–Avon Transition
         15: +2.1,  # Buxton–Avon Transition
         16: +2.3,  # Buxton–Avon Transition
         17: +2.4,  # Buxton–Avon Transition
         18: +2.4,  # Buxton–Avon Transition
         19: +2.1,  # Buxton–Avon Transition
         20: +1.7,  # Buxton–Avon Transition
         21: +0.4,  # Avon
         22: -1.9,  # Avon
         23: +0.0,  # Avon
         24: +0.0,  # Avon
         25: +0.0,  # Avon
         26: +0.0,  # Avon
         27: -1.1,  # Avon
         28: +1.2,  # Avon
         29: +2.4,  # Avon
         30: +3.2,  # Avon
         31: +4.0,  # Avon
         32: +4.2,  # Mid-island
         33: +3.8,  # Mid-island
         34: +3.5,  # Mid-island
         35: +3.0,  # Mid-island
         36: +2.6,  # Mid-island
         37: +2.3,  # Mid-island
         38: +2.0,  # Mid-island
         39: +2.3,  # Mid-island
         40: +2.0,  # Mid-island
         41: +1.7,  # Mid-island
         42: +1.4,  # Mid-island
         43: +1.0,  # Mid-island
         44: +0.6,  # Mid-island
         45: +0.0,  # Mid-island
         46: +0.0,  # Mid-island
         47: +0.0,  # Mid-island
         48: +0.0,  # Mid-island
         49: +0.0,  # Mid-island
         50: -1.3,  # Mid-island
         51: -1.6,  # Mid-island
         52: -2.2,  # Mid-island
         53: -2.0,  # Mid-island
         54: -1.8,  # Mid-island
         55: -1.0,  # Mid-island
         56: +0.0,  # Mid-island
         57: +0.0,  # Mid-island
         58: +0.0,  # Mid-island
         59: +0.0,  # Mid-island
         60: +0.0,  # Wimble Shoals Influence
         61: +0.0,  # Wimble Shoals Influence
         62: +0.0,  # Wimble Shoals Influence
         63: +0.7,  # Wimble Shoals Influence
         64: +1.5,  # Wimble Shoals Influence
         65: +1.9,  # Wimble Shoals Influence
         66: +2.3,  # Wimble Shoals Influence
         67: +2.3,  # Wimble Shoals Influence
         68: +3.0,  # Wimble Shoals Influence
         69: +3.8,  # Wimble Shoals Influence
         70: +4.1,  # Wimble Shoals Influence
         71: +4.3,  # Wimble Shoals Influence
         72: +4.1,  # Wimble Shoals Influence
         73: +2.8,  # Wimble Shoals Influence
         74: +2.4,  # Wimble Shoals Influence
         75: +2.4,  # Tri-Village / Rodanthe
         76: +1.5,  # Tri-Village / Rodanthe
         77: +1.2,  # Tri-Village / Rodanthe
         78: +0.9,  # Tri-Village / Rodanthe
         79: +0.4,  # Tri-Village / Rodanthe
         80: +0.0,  # Tri-Village / Rodanthe
         81: +0.0,  # Tri-Village / Rodanthe
         82: +0.0,  # Tri-Village / Rodanthe
         83: -2.6,  # Tri-Village / Rodanthe
         84: -3.5,  # Pea Island NWR
         85: -3.8,  # Pea Island NWR
         86: -3.7,  # Pea Island NWR
         87: -3.9,  # Pea Island NWR
         88: -3.0,  # Pea Island NWR
         89: -1.8,  # Pea Island NWR
         90: +57.9,  # LOCKED — end domain, LRR-solved; see the end-domain note above
    },
}

# Falls back to the 1984 fit only if 2004 was never solved separately
if HATTERAS_BE_RATES_CALIBRATED.get(2004) is None:
    HATTERAS_BE_RATES_CALIBRATED[2004] = dict(HATTERAS_BE_RATES_CALIBRATED[1984])
    HATTERAS_BE_RATES_2004_IS_PLACEHOLDER = True
else:
    HATTERAS_BE_RATES_2004_IS_PLACEHOLDER = False

# GIS 90 is the one end value the two presets do not share; GIS 1 is shared
HATTERAS_BE_EDGE_SHARED_DOMAINS = (1,)
HATTERAS_BE_EDGE_SPLIT_DOMAINS = (90,)

# GIS 90 under edgeBE, solved on the edgeBE road_bdm base run
HATTERAS_BE_EDGE_D90 = {
    1984: +13.0,
    2004: +46.7,
}

HATTERAS_BE_RATES_EDGE = {
    period: {**{gis: rates[gis]
                for gis in HATTERAS_BE_EDGE_SHARED_DOMAINS if gis in rates},
             90: HATTERAS_BE_EDGE_D90[period]}
    for period, rates in HATTERAS_BE_RATES_CALIBRATED.items()
}

# Periods with no calibrated fit carry their two end values explicitly
HATTERAS_BE_EDGE_ONLY = {
    # (GIS 1, GIS 90), m/yr: the adopted model with split12 storms (solve history in the README)
    # Lowess-7 option a values on the pre-adoption model, superseded 2026-09-28

    # Pre-adoption option A values on the LOWESS-7 target
    # Lowess-10 option a values, superseded 2026-09-28

    # Option A values on the LOWESS-10 target
    # /10-offset solve, 1996, superseded 2026-09-27

    # 1996: 1996-2015, five Newton steps (candidate-windows experiment, 2026-10-02)
    # 1996-2009: solved 2026-10-05 on the DEM-to-DEM net change (raw GIS 1, LOWESS-7 GIS 90), five secant steps from zeroBE
    # (experiments/end-domain-boundaries/2026-10-05-ends-solved-on-net-change-1996_2009). 1996-2015: (+3.39, +37.60)
    1996: (+1.4981, +10.5659),
    # 1996-2010 before it: (+4.3888, +19.0935), split12 storms, 2026-09-29; trim24 (+4.3509, +19.0935); adopted model, 2026-09-28; pre-adoption LOWESS-7 (+4.8394, +18.2545); LOWESS-10 +17.545; /10 (+32.2, +10.0)

    # /10-offset solve, 2010, superseded 2026-09-27

    # 2010: 2010-2026, five Newton steps (candidate-windows experiment, 2026-10-02); GIS 1 gain fell to ~0.01 near the top
    # GIS 90 re-solved 2026-10-04 after the run-length fix and the Rodanthe 82-88 / Buxton 6-16 footprints: 24.08 -> 38.3
    # (probes 28.1 -0.655, 39.0 +0.088, 37.7 -0.065, 38.3 -0.051 m/yr; the response is noise-limited near here). 1996 kept at
    # +37.60 (-0.144): probes 32.0-40.7 all scored worse, no trend (experiments/end-domain-boundaries/2026-10-04-gis90-runlength-footprints)
    # 2009-2025: the TEST period carries the calibration ends unchanged (2026-10-05, the DEM-to-DEM plan); not solved here.
    # Its residual becomes the second source/sink set. 2010-2026: (+172.89, +38.3)
    2009: (+1.4981, +10.5659),
    # before 2026-10-04: (+172.89, +24.08); 2010-2024 before it: (+8.0405, +21.2582), split12 storms, 2026-09-29; trim24 (+8.0, +21.2582) after the dune-cap fix, 2026-09-28; adopted before it (+8.0, +22.4937); pre-adoption LOWESS-7 (+18.8657, +24.2358); LOWESS-10 (+18.8, +24.535); /10 (+72.6, +31.3)
}

# Option B (2010-2024 at Hs 2.5 with its own ends): recorded, not wired
HATTERAS_WAVE_OPTION_B = {
    2010: {
        "hs": 2.5,
        "wave_period_s": 7.5,
        "wave_asymmetry": 0.6,
        "wave_angle_high_fraction": 0.5,
        "be_edge": (+8.0, +40.399),   # (GIS 1, GIS 90), m/yr
    },
}

for _period, (_d1, _d90) in HATTERAS_BE_EDGE_ONLY.items():
    if _period in HATTERAS_BE_RATES_EDGE:
        raise ValueError(
            f"period {_period} has both a calibrated fit and an entry in "
            f"HATTERAS_BE_EDGE_ONLY, so GIS 1 has two homes. Delete the "
            f"edge-only entry -- the calibrated preset is the one the "
            f"comprehension above keeps in step.")
    HATTERAS_BE_RATES_EDGE[_period] = {1: _d1, 90: _d90}

# An extended geometry keeps standing values only at domains that are still ends
if HATTERAS_GEOMETRY_EXTENDED:
    HATTERAS_BE_RATES_EDGE = {
        _period: {gis: rate for gis, rate in _rates.items()
                  if gis in HATTERAS_BE_EDGE_DOMAINS}
        for _period, _rates in HATTERAS_BE_RATES_EDGE.items()}

for _period, _rates in (() if HATTERAS_GEOMETRY_EXTENDED
                        else HATTERAS_BE_RATES_EDGE.items()):
    _absent = [gis for gis in HATTERAS_BE_EDGE_DOMAINS if gis not in _rates]
    if _absent:
        raise ValueError(
            f"edgeBE {_period}: end domains {_absent} are not in the "
            f"calibrated preset, so the edge preset would silently model them "
            f"at 0.0 m/yr. Add them to HATTERAS_BE_RATES_CALIBRATED.")

    # An edge domain present but zero is refused: edgeBE would equal zeroBE there
    _zeroed = [gis for gis in HATTERAS_BE_EDGE_DOMAINS
               if gis in _rates and _rates[gis] == 0.0]
    if _zeroed:
        raise ValueError(
            f"edgeBE {_period}: end domains {_zeroed} are present but 0.0, so "
            f"edgeBE is identical to zeroBE there and the edge secant has no "
            f"second bracket. Give them a nonzero starting value, or run "
            f"zeroBE if that is what you actually want.")

# The GIS 90 split must stay earned: each side solved against its own preset
for _period in HATTERAS_BE_RATES_CALIBRATED:
    if _period not in HATTERAS_BE_EDGE_D90:
        raise ValueError(
            f"edgeBE {_period}: GIS 90 is a SPLIT domain but has no entry in "
            f"HATTERAS_BE_EDGE_D90, so the edge preset has no value solved "
            f"against its own base run.")

# Canonical presets: the keys are the run-name tokens
HATTERAS_BE_PRESETS = {
    "zeroBE": HATTERAS_BE_RATES_ZERO,
    "edgeBE": HATTERAS_BE_RATES_EDGE,
    "calibBE": HATTERAS_BE_RATES_CALIBRATED,
}

# Deprecated spellings; resolve_be_preset() maps them to the canonical key
HATTERAS_BE_PRESET_ALIASES = {
    "base": "zeroBE",          # "base" was all-zeros, with the edge values commented out beside them
    "calibrated": "calibBE",
}


def resolve_be_preset(name):
    """Resolves a source/sink preset name to its canonical key and rates.

    Args:
        name: A canonical preset key or a deprecated alias.

    Returns:
        A (canonical_name, rates_by_period) tuple. rates_by_period maps start
        year to a sparse {gis_id: rate_m_yr} dict.

    Raises:
        ValueError: If the name is neither a canonical key nor an alias.
    """
    canonical = HATTERAS_BE_PRESET_ALIASES.get(name, name)
    if canonical not in HATTERAS_BE_PRESETS:
        raise ValueError(
            f"unknown source/sink preset {name!r}; expected one of "
            f"{sorted(HATTERAS_BE_PRESETS)} "
            f"(deprecated aliases: {sorted(HATTERAS_BE_PRESET_ALIASES)})")
    return canonical, HATTERAS_BE_PRESETS[canonical]


def be_rates(name, start_year):
    """The per-domain rates one preset supplies for one period.

    WHY THIS IS NOT A PLAIN DICT LOOKUP. A period can be WIRED before it is
    CALIBRATED, and since 2026-09-11 two of them are: all four periods resolve
    their forcing files, but only 1984 and 2004 have solved edge and calibrated
    presets. The other two carry zeroBE alone until their end domains are
    re-solved against their OWN CoastSat target.

    Indexing the preset dict directly gives a bare KeyError on the year, which
    reads like a typo in the settings file rather than what it is -- a fit that
    has not been done yet. This says which, and what would fix it.

    Args:
        name: A canonical preset key or a deprecated alias.
        start_year: The period's start year.

    Returns:
        A sparse {gis_id: rate_m_yr} dict. An absent domain is 0.0 m/yr.

    Raises:
        ValueError: If the preset is unknown, or is not solved for that period.
    """
    canonical, by_period = resolve_be_preset(name)
    if start_year not in by_period:
        raise ValueError(
            f"the {canonical} preset has no rates for the {start_year} "
            f"period; it is solved for {sorted(by_period)}. These are fitted "
            f"per period against that period's own observed rates, so another "
            f"period's numbers cannot stand in for them. Run {start_year} "
            f"under zeroBE, or solve its end domains first -- see "
            f"HATTERAS_BE_EDGE_D90 and the end-domain note above it.")
    return by_period[start_year]

# Nc-12 roadway

# GIS domains carrying NC-12 (none on Cape Point, GIS 1-8)
HATTERAS_FIRST_ROAD_DOMAIN = 9
HATTERAS_LAST_ROAD_DOMAIN = 90

# Permanent community zones: roadway management off inside the villages
HATTERAS_COMMUNITY_ZONES = (
    (7, 8),     # Buxton
    (21, 31),   # Avon
    (68, 83),   # Salvo / Waves / Rodanthe (Tri-Village)
)

# Per-domain road elevation, m MHW-relative (2009 lidar under the 2004 line); one set for every period
HATTERAS_ROAD_ELEVATION_FILE = _tv_mgmt.init_relpath(_tv_mgmt.ROAD_ELEVATION_FILE)

# Historical NC-12 events: relocations carry a displacement, not an absolute setback

# Per-domain measurement from HAT_road_relocation_distance.py, named for its line vintages
_RELOCATION_MEASUREMENT_FILE = _tv_mgmt.init_relpath(
    _tv_mgmt.road_relocation_file(1978, 2008))

# The 1999 event stops at GIS 14: GIS 15's lines cross, so its displacement is undefined


def _measured_displacements(gis_domains):
    """Reads the measured landward displacement for one relocation event.

    Args:
        gis_domains: The GIS domains the event moves. Hand-specified per event,
            NOT taken from the file: the measurement classifies ~40 domains
            'relocated' across the whole island, and which of them belong to
            which historical event is a fact about NC-12, not about the lines.

    Returns:
        A {gis: displacement_m} dict, holding mean_signed_landward_m from the
        measurement rounded to the nearest CELL_M (see ROUNDED TO WHOLE CELLS
        above); `measured_displacement_m` returns the unrounded value.

    Raises:
        ValueError: If a domain is absent from the measurement, or the
            measurement classifies it 'no_edit' / 'redigitized'. Either means
            the CSV cannot support a relocation there, and forcing one anyway
            would put a number the data does not contain into the model.
    """
    path = INIT_ROOT / _RELOCATION_MEASUREMENT_FILE
    with open(path, newline="", encoding="utf-8") as handle:
        measured = {int(row["domain"]): row for row in csv.DictReader(handle)}

    absent = [gis for gis in gis_domains if gis not in measured]
    if absent:
        raise ValueError(
            f"relocated domains {absent} are absent from {path.name}")

    unmoved = {gis: measured[gis]["classification"] for gis in gis_domains
               if measured[gis]["classification"] != "relocated"}
    if unmoved:
        raise ValueError(
            f"{path.name} classifies {unmoved} -- the two digitised lines are "
            f"copied or re-traced there, so they cannot measure a relocation")

    return {gis: round_to_cell(float(measured[gis]["mean_signed_landward_m"]))
            for gis in gis_domains}


# Rounding the displacements to whole cells, and what it cost at GIS 11

CELL_M = 10.0   # Barrier3D cell, m: a prescribed move is a whole number of these


def round_to_cell(displacement_m):
    """A displacement rounded to the nearest whole cell (half-up, sign-aware)."""
    import math
    sign = -1.0 if displacement_m < 0 else 1.0
    return sign * math.floor(abs(displacement_m) / CELL_M + 0.5) * CELL_M


def measured_displacement_m(gis):
    """The UNROUNDED measured displacement for one domain, for labels and
    reports that want to show the measurement beside the forcing."""
    path = INIT_ROOT / _RELOCATION_MEASUREMENT_FILE
    with open(path, newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            if int(row["domain"]) == gis:
                return float(row["mean_signed_landward_m"])
    raise KeyError(gis)


HATTERAS_ROAD_EVENTS = (
    RelocationEvent(
        year=1989,
        displacement_m=_measured_displacements((84, 85, 86, 87)),
        note="NC-12 relocated landward 1989, Pea Island (GIS 84-87)",
    ),
    RelocationEvent(
        year=1999,
        displacement_m=_measured_displacements((9, 10, 11, 12, 13, 14)),
        note="NC-12 relocated landward 1999, inter-village south (GIS 9-14)",
    ),
    BridgeEvent(
        year=2022,
        gis_domains=tuple(range(82, 89)),
        note="Jug Handle Bridge 2022; road removed GIS 82-88 -> unmanaged",
    ),
)

# Report-only check: a relocated setback against the measured 2004 setback


def _relocation_check(period_start=2004):
    """Reads the measured setback for every domain a relocation event moves.

    The domain list is taken from HATTERAS_ROAD_EVENTS rather than retyped, so
    adding a relocation adds its cross-check automatically.

    Args:
        period_start: Start year whose road_setback_file supplies the measured
            values. Both relocation events precede 2004, so period 2 is the
            first same-year measurement that postdates them.

    Returns:
        A {gis: measured_setback_m} dict.

    Raises:
        ValueError: If a relocated domain is absent from the setback file.
            load_road_setbacks fills an absent domain with 0.0, which would
            read here as "the road sits on the dune line" -- a plausible-
            looking number for a measurement that was never made.
    """
    path = INIT_ROOT / HATTERAS_PERIODS[period_start]["road_setback_file"]
    setbacks, missing = load_road_setbacks(
        path, HATTERAS_DOMAINS,
        HATTERAS_FIRST_ROAD_DOMAIN, HATTERAS_LAST_ROAD_DOMAIN)

    relocated = sorted({gis for event in HATTERAS_ROAD_EVENTS
                        if isinstance(event, RelocationEvent)
                        for gis in event.displacement_m})
    absent = sorted(set(relocated) & set(missing))
    if absent:
        raise ValueError(
            f"relocated domains {absent} are not in {path.name}, so their "
            f"cross-check would report 0.0 m as a measured setback")

    return {gis: float(setbacks[HATTERAS_DOMAINS.gis_to_pad(gis)])
            for gis in relocated}


HATTERAS_RELOCATION_CHECK_2004 = _relocation_check()

# Beach nourishment

# Nourishment projects (Hatteras_Management_Timelines.xlsx); volumes in yd3, spread over the domains
HATTERAS_NOURISHMENT_PROJECTS = (
    NourishmentProject(
        name="Rodanthe emergency fill",
        year=2014,
        # GIS 82-88 since 2026-10-04 (was 84-89), from the 2013 USACE notice: 2.13 mi "from 1.5 miles north of the
        # Pea Island NWR border into the Mirlo Beach community to just north of the Rodanthe pier"; the north limit
        # falls ~380 m into GIS 88, the south at the GIS 82 south edge. The CoastSat change agrees (nourishment/reported_extent)
        gis_domains=tuple(range(82, 89)),
        volume_cubic_yards=1_620_000,
        note="Mirlo Beach S-curves, 2.13 mi from 1.5 mi north of the Pea Island refuge border south into Mirlo Beach",
    ),
    NourishmentProject(
        name="Buxton beach nourishment",
        year=2017,
        # Placed 2017-06-21 to 2018-02-27 (~46% by Nov 2017); fired in the start year. Same 2.9 mi as 2022
        # GIS 6-16 since 2026-10-04: southernmost groin ~200 m into GIS 6, Haulover Day Use Area ~150 m into GIS 16
        # Observed sand from the GIS 6 end moved south past the groin to GIS 4-5 within ~1 yr; the model keeps it at
        # GIS 6. Kept as reported, the mismatch reported (2026-10-05; 4-mgmt-forcing/README.md)
        gis_domains=tuple(range(6, 17)),
        volume_cubic_yards=2_600_000,
        note="Haulover Day Use Area to the lighthouse groin, 2.9 mi; Outer Banks Voice 2018-03-01",
    ),
    NourishmentProject(
        name="Buxton shore protection",
        year=2022,
        # Haulover Day Use Area to the lighthouse groin field, the 2017 footprint: GIS 6-16 since 2026-10-04 (was 6-15)
        gis_domains=tuple(range(6, 17)),
        volume_cubic_yards=1_200_000,
        note="Extends north out of Buxton village (7-8) into the road corridor",
    ),
    NourishmentProject(
        name="Avon shore protection",
        year=2022,
        # Avon: 3,000 ft north of the pier south to the village boundary, GIS 21-28
        gis_domains=tuple(range(21, 29)),
        volume_cubic_yards=1_000_000,
        note="Due East Rd to Askins Creek N Dr; inside the Avon zone (21-31)",
    ),
)

# Overwash filter for developed ground, a PERCENT: 40, the residential end of Rogers et al. (2015)
HATTERAS_BEACH_DUNE = BeachDuneConfig(
    community_overwash_filter_pct=40.0,
    default_overwash_filter_pct=0.0,
    overwash_to_dune_pct=9.0,
)
