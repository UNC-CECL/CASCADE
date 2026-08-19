"""Hatteras Island site config: the real place names, spans, and labels.

This is application config, not library code -- it imports cascade_pipeline's
generic dataclasses and fills them in with Hatteras-specific content, the
same way any other CASCADE study site would. Nothing in cascade_pipeline itself
knows these values exist; a different site (Ocracoke, etc.) would write its
own sibling module in this same shape and never touch the package.

Import these presets from your run script / notebook:
    from hatteras_site_config import HATTERAS_DOMAINS, HATTERAS_ANNOTATIONS
"""

from cascade_pipeline.annotations import AnnotationConfig
from cascade_pipeline.domains import DomainGeometry
from cascade_pipeline.nourishment import BeachDuneConfig, NourishmentProject
from cascade_pipeline.roadway import BridgeEvent, RelocationEvent

# Real domains GIS 1-90 (500 m each, south to north, Cape Point to Pea
# Island), padded by 15 buffer domains on each side. These happen to match
# DomainGeometry's own field defaults, but naming the instance explicitly
# here (rather than relying on cascade_pipeline.domains.DEFAULT_DOMAINS) keeps
# "this is Hatteras' geometry" visible at the call site.
HATTERAS_DOMAINS = DomainGeometry(
    num_real_domains=90,
    num_buffer_domains=15,
    first_gis_id=1,
    domain_spacing_m=500.0,
)

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
    low_end_label="S | Cape Hatteras",
    high_end_label="Pea Island | N",
    obs_source_name="CoastSat",
)


# =============================================================================
# Period forcing config
# =============================================================================
# The hindcast runs as one of two periods. Picking the start year resolves
# every period-dependent forcing: run length, RSLR rate, storm series, the
# BRIE island-offset file that sets the starting shoreline, the road setback
# file, nourishment defaults, and which background-erosion preset applies.
#
# Paths are relative to data/hatteras_init/ so this module stays independent
# of where the repo is checked out; join them onto your own data base dir.
#
# ROAD SETBACK: switched to the DUNE-START method on 2026-08-18.
#
# The setback is metres landward of interior row 0 -- the row
# roadway_manager.py:99 indexes against -- and the dune-start method is the
# first one that measures it there. The legacy files measured against the
# same-year digitised dune line instead, a different feature from a different
# year, which put the road a median +40 m (1984) / +23 m (2004) landward of the
# road actually rasterized onto the model grid, reaching 130 m.
#
# TWO THINGS THAT MUST MOVE TOGETHER WITH THIS:
#   1. TOPO_DUNE_VERSION must be 2009_v3. These setbacks are measured against
#      2009_v3's interior row 0; spend them on a different extraction and the
#      reference they were measured from no longer exists.
#   2. The legacy files are still on disk under old_method_offset/ for the
#      method comparison. They are NOT interchangeable -- see
#      scripts/input_prep/4-mgmt-forcings/road_offset/README.md.
#
# Both are 2 rows x 82 cols, GIS 9-90, so the swap itself is a drop-in.

HATTERAS_PERIODS = {
    1984: {
        "end_year": 2004,
        "sea_level_rise_rate": 0.004,  # m/yr, from duck_rslr_analysis.py
        "storm_file": (
            "3-env-forcings/storms/hindcast_storms/1984_2004/"
            "1984_2004_storms_v3_72.npy"),
        "island_offset_file": (
            "2-brie-offset/hindcast_1984/Island_Dune_Offsets_1984_PADDED_120.csv"),
        "road_setback_file": (
            "4-mgmt-forcing/road_offset/dunestart_offset/1984/"
            "RoadSetback_1984_dunestart.csv"),
        "enable_nourishment": False,
        "nourishment_volume": 0,  # m^3/m
    },
    2004: {
        "end_year": 2024,
        "sea_level_rise_rate": 0.006,  # m/yr, from duck_rslr_analysis.py
        "storm_file": (
            "3-env-forcings/storms/hindcast_storms/2004_2024/"
            "2004_2024_storms_v3_72.npy"),
        "island_offset_file": (
            "2-brie-offset/hindcast_2004/Island_Dune_Offsets_2004_PADDED_120.csv"),
        "road_setback_file": (
            "4-mgmt-forcing/road_offset/dunestart_offset/2004/"
            "RoadSetback_2004_dunestart.csv"),
        "enable_nourishment": True,  # historical BN injected per-year in the time loop
        "nourishment_volume": 100,  # m^3/m passed to Cascade init
    },
}

# =============================================================================
# Background erosion (source/sink) rates, m/yr by GIS domain
# =============================================================================
# Sparse by design: a domain absent from a preset gets 0.0 m/yr. Sign follows
# cascade/brie_coupler.py: (-) = erosion, (+) = accretion, in m/yr.
#
# Three presets, one per hypothesis about where the alongshore sediment budget
# is unresolved -- they are the source/sink axis of the run matrix:
#
#   zeroBE   no source/sink anywhere. Whatever the shoreline does is what
#            Barrier3D + BRIE produce unaided.
#   edgeBE   only the two end domains, which absorb the open-boundary artifact
#            at the ends of the modelled reach. Nothing is imposed on the
#            interior, so an interior misfit is the model's own.
#   calibBE  the full per-domain fit against the CoastSat LRR target rates.
#
# Moved here from HAT_hindcast_1984_2024.py so the notebook and the
# run script read the same numbers. Note data/hatteras_init/7-source-sink/
# holds older partial copies of these -- they disagree with these values and
# the 2004 one is truncated mid-dict; these are the ones that have been run.

# "zeroBE" preset -- empty rather than {1: 0, 90: 0}, because the sparse
# contract already gives an absent domain 0.0 and an explicit zero reads like a
# value that was solved for.
HATTERAS_BE_RATES_ZERO = {
    1984: {},
    2004: {},
}

# GIS domains carrying the edge-only preset: the first and last REAL domains.
# Not the padded buffers, which stay 0.0 in every preset.
HATTERAS_BE_EDGE_DOMAINS = (1, 90)

# "calibBE" preset -- the per-domain fit. HATTERAS_BE_RATES_EDGE is derived
# from this below, so the two presets cannot disagree about the end domains.
HATTERAS_BE_RATES_CALIBRATED = {
    1984: {
          1: -24.0,  # LOCKED — use your solved value, not 0.0 # -24
          2: -0.6,  # Cape Point / Shoal Dynamics
          3: -0.7,  # Cape Point / Shoal Dynamics
          4: 0.0,
          5: +1.2,  # Cape Point / Shoal Dynamics
          6: +2.0,  # Cape Point / Shoal Dynamics
          7: +0.9,  # Cape Point / Shoal Dynamics
          8: +1.1,  # Cape Point / Shoal Dynamics
          9: 0.0,
         10: -0.9,  # Cape Point / Shoal Dynamics
         11: -1.3,  # Buxton–Avon Transition
         12: -0.9,  # Buxton–Avon Transition
         13: -0.6,  # Buxton–Avon Transition
         14: 0.0,
         15: 0.0,
         16: 0.0,
         17: +0.6,  # Buxton–Avon Transition
         18: +0.6,  # Buxton–Avon Transition
         19: +0.6,  # Buxton–Avon Transition
         20: 0.0,
         21: 0.0,
         22: 0.0,
         23: 0.0,
         24: +0.3,  # Avon
         25: +0.3,  # Avon
         26: +0.8,  # Avon
         27: +1.4,  # Avon
         28: +1.9,  # Avon
         29: +2.3,  # Avon
         30: +2.2,  # Avon
         31: +2.2,  # Avon
         32: +1.9,  # Mid-island
         33: +1.4,  # Mid-island
         34: +1.0,  # Mid-island
         35: 0.0,
         36: 0.0,
         37: 0.0,
         38: 0.0,
         39: 0.0,
         40: 0.0,
         41: 0.0,
         42: 0.0,
         43: 0.0,
         44: +0.4,  # Mid-island
         45: 0.0,
         46: 0.0,
         47: 0.0,
         48: -0.3,  # Mid-island
         49: -1.0,  # Mid-island
         50: -1.0,  # Mid-island
         51: -1.2,  # Mid-island
         52: -1.3,  # Mid-island
         53: -1.3,  # Mid-island
         54: -1.2,  # Mid-island
         55: -1.0,  # Mid-island
         56: -1.0,  # Mid-island
         57: -0.3,  # Mid-island
         58: 0.0,
         59: 0.0,
         60: 0.0,
         61: 0.0,
         62: +0.4,  # Wimble Shoals Influence
         63: 0.0,
         64: 0.0,
         65: 0.0,
         66: 0.0,
         67: 0.0,
         68: +0.7,  # Wimble Shoals Influence
         69: +1.2,  # Wimble Shoals Influence
         70: +1.6,  # Wimble Shoals Influence
         71: +1.9,  # Wimble Shoals Influence
         72: +2.2,  # Wimble Shoals Influence
         73: +2.4,  # Wimble Shoals Influence
         74: +2.6,  # Wimble Shoals Influence
         75: +2.5,  # Tri-Village / Rodanthe
         76: +1.7,  # Tri-Village / Rodanthe
         77: +1.3,  # Tri-Village / Rodanthe
         78: +0.7,  # Tri-Village / Rodanthe
         79: 0.0,
         80: -1.3,  # Tri-Village / Rodanthe
         81: -1.9,  # Tri-Village / Rodanthe
         82: -2.1,  # Tri-Village / Rodanthe
         83: -2.3,  # Tri-Village / Rodanthe
         84: -2.5,  # Pea Island NWR
         85: -2.5,  # Pea Island NWR
         86: -2.4,  # Pea Island NWR
         87: -2.0,  # Pea Island NWR
         88: -1.6,  # Pea Island NWR
         89: -1.0,  # Pea Island NWR
         90: 10.0,  # LOCKED — use your solved value, not 0.0 # 10
    },

    2004: {
          1: 35.0,  # LOCKED — use your solved value, not 0.0
          2: +1.0,  # Cape Point / Shoal Dynamics
          3: +1.6,  # Cape Point / Shoal Dynamics
          4: +1.6,  # Cape Point / Shoal Dynamics
          5: +1.2,  # Cape Point / Shoal Dynamics
          6: -1.5,  # Cape Point / Shoal Dynamics
          7: 0.0,
          8: +1.1,  # Cape Point / Shoal Dynamics
          9: +2.2,  # Cape Point / Shoal Dynamics
         10: +3.1,  # Cape Point / Shoal Dynamics
         11: +3.4,  # Buxton–Avon Transition
         12: +3.4,  # Buxton–Avon Transition
         13: +3.3,  # Buxton–Avon Transition
         14: +3.2,  # Buxton–Avon Transition
         15: +3.0,  # Buxton–Avon Transition
         16: +3.1,  # Buxton–Avon Transition
         17: +3.0,  # Buxton–Avon Transition
         18: +2.9,  # Buxton–Avon Transition
         19: +2.6,  # Buxton–Avon Transition
         20: +2.2,  # Buxton–Avon Transition
         21: +1.9,  # Avon
         22: +1.4,  # Avon
         23: +0.8,  # Avon
         24: +0.3,  # Avon
         25: +0.3,  # Avon
         26: +0.8,  # Avon
         27: +1.4,  # Avon
         28: +1.9,  # Avon
         29: +2.3,  # Avon
         30: +3.1,  # Avon
         31: +3.6,  # Avon
         32: +3.8,  # Mid-island
         33: +3.7,  # Mid-island
         34: +3.6,  # Mid-island
         35: +3.5,  # Mid-island
         36: +3.3,  # Mid-island
         37: +3.1,  # Mid-island
         38: +2.7,  # Mid-island
         39: +2.4,  # Mid-island
         40: +2.1,  # Mid-island
         41: +1.8,  # Mid-island
         42: +1.5,  # Mid-island
         43: +1.1,  # Mid-island
         44: +0.4,  # Mid-island
         45: 0.0,
         46: 0.0,
         47: 0.0,
         48: -0.3,  # Mid-island
         49: 0.0,
         50: -1.0,  # Mid-island
         51: -1.2,  # Mid-island
         52: -1.3,  # Mid-island
         53: -1.3,  # Mid-island
         54: -1.2,  # Mid-island
         55: -1.0,  # Mid-island
         56: 0.0,
         57: -0.3,  # Mid-island
         58: 0.0,
         59: 0.0,
         60: 0.0,
         61: 0.0,
         62: +0.4,  # Wimble Shoals Influence
         63: +1.1,  # Wimble Shoals Influence
         64: +1.6,  # Wimble Shoals Influence
         65: +2.0,  # Wimble Shoals Influence
         66: +2.4,  # Wimble Shoals Influence
         67: +2.8,  # Wimble Shoals Influence
         68: +3.3,  # Wimble Shoals Influence
         69: +3.7,  # Wimble Shoals Influence
         70: +4.0,  # Wimble Shoals Influence
         71: +4.2,  # Wimble Shoals Influence
         72: +4.1,  # Wimble Shoals Influence
         73: +4.0,  # Wimble Shoals Influence
         74: +3.7,  # Wimble Shoals Influence
         75: +3.4,  # Tri-Village / Rodanthe
         76: +3.0,  # Tri-Village / Rodanthe
         77: +2.7,  # Tri-Village / Rodanthe
         78: +2.4,  # Tri-Village / Rodanthe
         79: +2.3,  # Tri-Village / Rodanthe
         80: +2.0,  # Tri-Village / Rodanthe
         81: +1.7,  # Tri-Village / Rodanthe
         82: +1.4,  # Tri-Village / Rodanthe
         83: +1.0,  # Tri-Village / Rodanthe
         84: +0.7,  # Pea Island NWR
         85: 0.0,
         86: +0.5,  # Pea Island NWR
         87: +0.6,  # Pea Island NWR
         88: +0.7,  # Pea Island NWR
         89: +1.2,  # Pea Island NWR
         90: 35.0,  # LOCKED — use your solved value, not 0.0
    },  #
}

# The 2004 calibrated preset falls back to the 1984 fit if it was never
# solved separately. It has been, so this stays False -- but the flag is
# what tells you whether a "calibrated 2004" run is really calibrated.
if HATTERAS_BE_RATES_CALIBRATED.get(2004) is None:
    HATTERAS_BE_RATES_CALIBRATED[2004] = dict(HATTERAS_BE_RATES_CALIBRATED[1984])
    HATTERAS_BE_RATES_2004_IS_PLACEHOLDER = True
else:
    HATTERAS_BE_RATES_2004_IS_PLACEHOLDER = False

# The edge-only preset is SLICED from the calibrated fit rather than retyped,
# so re-solving an end domain updates both presets at once. Retyping is how the
# old "base" preset drifted into holding commented-out copies of these numbers.
HATTERAS_BE_RATES_EDGE = {
    period: {gis: rates[gis]
             for gis in HATTERAS_BE_EDGE_DOMAINS if gis in rates}
    for period, rates in HATTERAS_BE_RATES_CALIBRATED.items()
}

for _period, _rates in HATTERAS_BE_RATES_EDGE.items():
    _absent = [gis for gis in HATTERAS_BE_EDGE_DOMAINS if gis not in _rates]
    if _absent:
        raise ValueError(
            f"edgeBE {_period}: end domains {_absent} are not in the "
            f"calibrated preset, so the edge preset would silently model them "
            f"at 0.0 m/yr. Add them to HATTERAS_BE_RATES_CALIBRATED.")

# Canonical presets. Keys are the tokens that land in run-directory names, so
# each one states which hypothesis was run -- "base" did not.
HATTERAS_BE_PRESETS = {
    "zeroBE": HATTERAS_BE_RATES_ZERO,
    "edgeBE": HATTERAS_BE_RATES_EDGE,
    "calibBE": HATTERAS_BE_RATES_CALIBRATED,
}

# Deprecated spellings, kept so older scripts keep running. resolve_be_preset()
# maps these to the canonical key, and it is the canonical key that reaches the
# run name -- an alias can never put a stale token in a directory name.
HATTERAS_BE_PRESET_ALIASES = {
    "base": "zeroBE",          # "base" was all-zeros, with the edge values
                               # commented out beside them
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

# =============================================================================
# NC-12 roadway
# =============================================================================
# GIS domains carrying NC-12. Domains 1-8 (Cape Point) have no road in the
# modelled span.
HATTERAS_FIRST_ROAD_DOMAIN = 9
HATTERAS_LAST_ROAD_DOMAIN = 90

# Permanent community zones. Roadway management is OFF here: inside a village
# the road is a street network that is maintained rather than relocated, so
# CASCADE's relocate-or-abandon logic does not describe it.
#
# These are the PERMANENT settlement footprints only, deliberately narrower
# than the beach-nourishment footprints, which extend past the villages.
HATTERAS_COMMUNITY_ZONES = (
    (7, 8),     # Buxton
    (21, 31),   # Avon
    (68, 83),   # Salvo / Waves / Rodanthe (Tri-Village)
)

# Per-domain road elevation, meters MHW-RELATIVE: the MEAN of the 2009 LiDAR
# 1 m cells in a 7 m corridor under the digitised 2004 NC-12 alignment.
# Written by scripts/input_prep/4-mgmt-forcings/road_elevation/
# HAT_road_elevation.py; see RoadElevation_audit.md beside the file.
#
# NOT period-dependent, and that is deliberate: road_ele is a property of the
# surveyed surface, and there is one topography (2009) for every period, so one
# elevation set serves both. Hence no year in the filename -- a per-year name
# would imply a measured change in roadbed height that no data supports.
#
# The 2004 alignment is used rather than 1984 because it is the only digitised
# line that lands on roadbed everywhere on the 2009 surface. Sampled along the
# 1984 line, the relocated domains (GIS 9-15, 84-87) return 2.37 m mean with a
# within-domain sigma up to 1.70 m -- the abandoned corridor now lies UNDER the
# foredune, so that sample is dune, not road. The 2004 value is both the better
# measurement and the lower of the two.
HATTERAS_ROAD_ELEVATION_FILE = "4-mgmt-forcing/road_elevation/RoadElevation.csv"

# Historical NC-12 management events.
#
# Relocations carry a DISPLACEMENT, not an absolute setback. The displacement
# is the 1978->1997 cross-shore offset measured between the two digitised road
# lines in ArcGIS Pro. CASCADE adds it to whatever setback the model is
# carrying at the event year, which already reflects the modelled shoreline
# retreat -- so the retreat is counted once. An absolute setback referenced to
# the 1984 dune line would count it twice, because the topography's own dune
# line has already moved landward by that amount.
HATTERAS_ROAD_EVENTS = (
    RelocationEvent(
        year=1989,
        displacement_m={84: 70.0, 85: 120.0, 86: 120.0, 87: 25.0},
        note="NC-12 relocated landward 1989, Pea Island (GIS 84-87)",
    ),
    RelocationEvent(
        year=1999,
        displacement_m={9: 30.0, 10: 55.0, 11: 80.0, 12: 70.0, 13: 55.0,
                        14: 30.0, 15: 5.0},
        note="NC-12 relocated landward 1999, inter-village south (GIS 9-15)",
    ),
    BridgeEvent(
        year=2022,
        gis_domains=tuple(range(82, 89)),
        note="Jug Handle Bridge 2022; road removed GIS 82-88 -> unmanaged",
    ),
)

# Independent check on a relocated setback: the measured RoadSetback_2004.csv
# value for the same domain. Both relocation events precede 2004, so a
# correctly displaced setback should land near the 2004 same-year measurement.
# Reporting only -- nothing reads this to decide anything.
HATTERAS_RELOCATION_CHECK_2004 = {
     9: 89.0, 10: 83.0, 11: 81.0, 12: 89.0, 13: 87.0, 14: 93.0, 15: 71.0,
    84: 50.0, 85: 85.0, 86: 88.0, 87: 40.0,
}

# =============================================================================
# Beach nourishment
# =============================================================================
# Source: Hatteras_Management_Timelines.xlsx -> Nourishment_Timeline sheet.
# Volumes are the reported project totals in cubic yards; cascade_pipeline
# converts to m^3/m and spreads them evenly across each project's domains.
#
# All three projects fall in 2004-2024, so a 1984-2004 run builds an empty
# schedule from this same list -- no period-keying needed here.
#
# The project extents are wider than HATTERAS_COMMUNITY_ZONES on purpose:
# these are engineering footprints, and both Buxton and Rodanthe extend past
# the settlement into the NC-12 corridor. That is what puts GIS 9-15 and
# 85-88 under both the roadway and beach-dune managers.
HATTERAS_NOURISHMENT_PROJECTS = (
    NourishmentProject(
        name="Rodanthe emergency fill",
        year=2014,
        gis_domains=tuple(range(85, 89)),
        volume_cubic_yards=1_620_000,
        note="Emergency fill north of Tri-Village; GIS 85-88 carry NC-12",
    ),
    NourishmentProject(
        name="Buxton shore protection",
        year=2022,
        gis_domains=tuple(range(6, 16)),
        volume_cubic_yards=1_200_000,
        note="Extends north out of Buxton village (7-8) into the road corridor",
    ),
    NourishmentProject(
        name="Avon shore protection",
        year=2022,
        gis_domains=tuple(range(23, 27)),
        volume_cubic_yards=2_200_000,
        note="Entirely inside the Avon community zone (21-31)",
    ),
)

# Overwash filtering percent for developed ground. CASCADE's overwash_filter
# is a PERCENT -- filter_overwash divides by 100 -- and Rogers et al. (2015)
# give 40-90%, residential to commercial. Hatteras' villages are residential,
# so 40 is the low end of that range.
#
# Earlier versions of this pipeline passed 0.4 here, on the assumption it was
# a fraction. That filtered 0.4% of overwash, which is indistinguishable from
# no filtering at all; BeachDuneConfig now rejects the fraction scale.
HATTERAS_BEACH_DUNE = BeachDuneConfig(
    community_overwash_filter_pct=40.0,
    default_overwash_filter_pct=0.0,
    overwash_to_dune_pct=9.0,
)
