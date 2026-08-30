"""Hatteras Island site config: the real place names, spans, and labels.

This is application config, not library code -- it imports cascade_pipeline's
generic dataclasses and fills them in with Hatteras-specific content, the
same way any other CASCADE study site would. Nothing in cascade_pipeline itself
knows these values exist; a different site (Ocracoke, etc.) would write its
own sibling module in this same shape and never touch the package.

Import these presets from your run script / notebook:
    from hatteras_site_config import HATTERAS_DOMAINS, HATTERAS_ANNOTATIONS
"""

import csv

from cascade_pipeline.annotations import AnnotationConfig
from cascade_pipeline.domains import DomainGeometry
from cascade_pipeline.nourishment import BeachDuneConfig, NourishmentProject
from cascade_pipeline.roadway import (
    BridgeEvent, RelocationEvent, load_road_setbacks)

# Sibling module in scripts/, which owns "where is data/hatteras_init". The
# relocation cross-check below reads the period-2 setback file rather than
# carrying a copy of its numbers, so this module needs the data root.
from hat_topo_version import INIT_ROOT, YEAR_PRODUCT

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
#   1. The topography must be the extraction these setbacks were measured
#      against -- they are metres landward of ITS interior row 0, so spending
#      them on another extraction measures from a row that does not exist.
#      This no longer has to be remembered: the runner, the sweep worker and
#      the road scripts all resolve the version through
#      scripts/hat_topo_version.py, which reads VERSION out of the extractor.
#      Bump it there and everything moves together. (This comment used to
#      pin 2009_v3 by hand and went stale the day the setbacks moved to v4,
#      then again when it still said "Current: 2009_v5" after the tree went
#      period-first. It no longer names a version at all -- ask topo_dirs().)
#      There are now TWO extractions, one per period; see "topo_product".
#   2. The legacy files are still on disk under old_method_offset/ for the
#      method comparison. They are NOT interchangeable -- see
#      scripts/input_prep/4-mgmt-forcings/road_offset/README.md.
#
# Both are 2 rows x 82 cols, GIS 9-90, so the swap itself is a drop-in.

# "topo_product" names the folder under
# data/hatteras_init/1-barrier3d-domains/ that this period's Barrier3D domains
# come from. It sits beside storm_file / island_offset_file / road_setback_file
# because it is the same kind of thing: a per-period input the run ingests.
#
# ADDED 2026-08-25. Before that the runner hardcoded ONE topography
# ("2009-dune-topo" / <version>) and BOTH periods read it, so a 1984 run and a
# 2004 run started from the same barrier. They no longer do:
#
#     1984  <- DEM 2009-2014-1996  (1996 ALACE overwriting measured ground
#                                   wherever ALACE has data; the landward limit
#                                   is the swath edge, near the dune toe. The
#                                   "ocean-side of the 1984 NC-12 line"
#                                   boundary this comment used to name was
#                                   dropped 2026-08-26 -- see that product's
#                                   README for the three-way measurement.)
#     2004  <- DEM 2009-2014       (the baseline gap-filled DEM)
#
# The version WITHIN a product is still resolved, never pinned - see
# scripts/hat_topo_version.py.
#
# NOT A LITERAL ANY MORE (2026-08-26). The value comes from YEAR_PRODUCT in
# hat_topo_version.py, which is the same mapping every road script in
# 4-mgmt-forcings now resolves through. It was spelled out in four places and
# omitted in three, and the three that omitted it - the placement figure, the
# method diagnostic and the per-domain views - gave BOTH vintages 2004-start
# interiors. Since 65 of 90 domains have a different interior shape between the
# products, that is a different island, not a rounding difference. One mapping,
# imported, so the runner and the forcing that feeds it cannot disagree.
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
        "topo_product": YEAR_PRODUCT[1984],
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
        "topo_product": YEAR_PRODUCT[2004],
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
#
# FIT AGAINST THE LRR, NOT THE ENDPOINT DIFFERENCE. These rates are the
# LOESS-smoothed residual of an edgeBE base run against the CoastSat
# target, and both sides of that residual are now the same estimator: the
# target is a per-transect OLS slope (transect_lrr_full.csv, lrr_m_yr) and
# the model side is read from the run's lrr_m_yr, its own OLS slope
# through 21 annual states. Before 2026-08-22 the model side was
# change_rate_m_yr -- (x[-1] - x[0]) / span -- so the residual carried the
# gap between two estimators as if it were a sediment budget.
#
# The refit moved 36 of 88 interior domains in 1984-2004 (mean 0.21, max
# 0.50 m/yr) and 65 of 88 in 2004-2024 (mean 0.39, max 1.60 m/yr). Period
# 2 moves further because that is the period with nourishment: a fill is
# an instantaneous step in x_s, BRIE answers a step with a slowly-decaying
# alongshore grid mode, and the 2022 Buxton and Avon fills are two years
# from the end of the run -- so the endpoint difference was reading solver
# ringing into the residual and calling it background erosion.
#
# Regenerate with scripts/input_prep/7-source-sink/loess_smooth/
# HAT_be_zone_LOESS_analysis.py; GIS 1 and 90 are NOT taken from its
# output, which writes them as 0.0 -- they stay the separately solved
# buffer-cell values below.
#
# THE END DOMAINS WERE RE-SOLVED AGAINST THE LRR ON 2026-08-23, for the
# same reason the interior was: the values they replace were solved when
# the model side was the endpoint difference, so they were fit to a
# different estimator than the one now plotted. Previous values were
# GIS 1 / 90 = -24.0 / +10.0 (1984) and +35.0 / +35.0 (2004).
#
# HOW THEY WERE SOLVED. Two Newton steps against the base run below: a
# first from the zeroBE/edgeBE secant, then a second on the local secant
# through the two bracketing runs. The gain is SMALL -- d(LRR)/d(BE) is
# 0.092 to 0.123 across the four cases, i.e. only about a TENTH of an
# imposed edge rate survives in that domain's own shoreline, the rest
# being diffused alongshore by BRIE within a few domains. So each value
# here is roughly
# ten times the misfit it is there to close, and is a boundary-artifact
# absorber, not a sediment budget: read as a real flux over the 22.25 m
# shoreface, GIS 1 in 1984-2004 would be ~5e5 m^3/yr across one 500 m
# domain. It is not capped for plausibility, because capping it would
# just move the open-boundary artifact back into the figure.
#
# THE GAIN IS NOT CONSTANT, which is why one step was not enough. It
# STIFFENS with the imposed rate at GIS 1 in 1984-2004 (0.098 over
# BE 0 to -24, but 0.123 over -24 to -46.4) and SOFTENS in 2004-2024
# (0.097 to 0.092). At -42 m/yr for 20 years GIS 1 is being pushed
# ~840 m landward, far outside the range where a background rate acts
# as a small perturbation, so a single global slope should not be
# assumed if these are ever re-solved.
#
# WHAT THEY WERE FIT TO. The value the rate-comparison figure DRAWS at
# each end, so fit and figure cannot disagree:
#   GIS 1   raw per-domain transect mean -- LoessConfig.skip_southern_
#           domains is 10, so D1-D10 are drawn raw, not smoothed.
#   GIS 90  the LOESS-10 value, the primary window, which is what is
#           drawn everywhere north of D10.
# The two ends therefore use different estimators. That is deliberate:
# it mirrors the splice the figure already makes.
#
# BASE RUN. edgeBE / road_bdm / groin off, per period -- the same run
# HAT_be_zone_LOESS_analysis.py derives the interior residual from.
#
# WHY GROIN OFF. This USED to say "because output/groin_sweep/joint_fit.json
# does not exist yet, so _fitted_groin() returns None and the analysis falls
# back". THAT IS NO LONGER TRUE -- the joint fit landed 2026-08-24 and pins
# M = 60, f = 0.6 for every preset (calibBE copies edgeBE's pair; the two that
# were fitted independently both landed on it). So groin-off is now a CHOICE,
# not a fallback, and it is the right one: the groin is a structure whose
# trapping is fitted against the same period-1 shoreline these BE rates are
# fitted against, and letting both absorb the same misfit makes neither
# identifiable. The BE fit goes first, groin off; the sweep then pins be1 at
# the production edgeBE value and fits the groin on top.
#
# The groin sits at GIS 5/6 and does not reach GIS 90 at all (on the pre-refit
# runs D90 was identical to the 3 dp reported, groin-on vs groin-off, in
# both periods); it moves GIS 1 by about 0.4 m/yr, so GIS 1 is worth
# re-checking against the fitted groin.
#
# STALE AS OF 2026-08-26 FOR PERIOD 1. Every rate below was derived from a
# base run on the pre-restructure shared topography (now 2004-start/v1). The
# 1984-start product is a different surface, so the 1984 rates must be refit
# before any calibBE or edgeBE run on it. Period 2 is unaffected -- its base
# run read 2004-start, which has not changed.
#
# D1 <-> D90 CROSS-TALK IS NEGLIGIBLE, so solving the two ends
# independently is safe: with 15 buffer domains a side they are 30
# domains -- 15 km -- apart around the ring, against a BRIE diffusion
# length of sqrt(D*t) ~ 3.2 km over a 20 year period.
# ============================================================================
# INTERIOR ZEROED 2026-08-28 — THE CALIBRATED FIT HAS NOT BEEN REDONE
# ============================================================================
# GIS 2-89 were set to +0.0 for BOTH periods on 2026-08-28. They are not a fit
# result; they are the absence of one. The previous values are in git history
# (`git show HEAD:scripts/hatteras_site_config.py`) and in
# `output/archived_output_20260828/`.
#
# WHY. Period 1's rates were derived on a base run against the pre-restructure
# shared topography, and `1984-start` is now a different surface with
# re-measured road setbacks (15 of 83 road domains moved on 2026-08-28). Rather
# than refit period 1 and leave period 2 on an older sitting, both are being
# re-solved together — the stage-5 groin joint fit intersects BOTH periods'
# surfaces, so mixed calibration vintages would make M and f mean two things.
#
# GIS 1 AND 90 ARE DELIBERATELY NOT ZEROED. They are `LOCKED_GIS` in
# `apply_be_fit.py`: solved by Newton steps on the edge secant, not by the
# LOESS interior residual, and they are the SECOND BRACKET that edge solve
# needs. `HATTERAS_BE_RATES_EDGE` slices GIS 1 straight out of this table, so
# zeroing it would make the edgeBE bracket run identical to zeroBE at GIS 1 and
# collapse the secant to 0/0 — silently, because the guard below only tests for
# ABSENT domains, not zeroed ones. The values standing here are the previous
# solution, kept as a starting bracket and expected to be overwritten.
#
# UNTIL THE REFIT LANDS, calibBE IS EFFECTIVELY zeroBE PLUS TWO EDGE CELLS.
# Do not read a calibBE run made in this window as a calibrated result.
HATTERAS_BE_RATES_CALIBRATED = {
    1984: {
          1: -42.6,  # LOCKED — end domain, LRR-solved; see the end-domain note above
          2: +0.0,  # Cape Point / Shoal Dynamics
          3: +0.0,  # Cape Point / Shoal Dynamics
          4: +0.0,  # Cape Point / Shoal Dynamics
          5: +0.0,  # Cape Point / Shoal Dynamics
          6: +0.0,  # Cape Point / Shoal Dynamics
          7: +0.0,  # Cape Point / Shoal Dynamics
          8: +1.2,  # Cape Point / Shoal Dynamics
          9: +0.0,  # Cape Point / Shoal Dynamics
         10: -1.1,  # Cape Point / Shoal Dynamics
         11: -0.7,  # Buxton–Avon Transition
         12: -0.7,  # Buxton–Avon Transition
         13: -0.6,  # Buxton–Avon Transition
         14: +0.0,  # Buxton–Avon Transition
         15: +0.0,  # Buxton–Avon Transition
         16: +0.0,  # Buxton–Avon Transition
         17: +0.0,  # Buxton–Avon Transition
         18: +0.0,  # Buxton–Avon Transition
         19: +0.0,  # Buxton–Avon Transition
         20: +0.0,  # Buxton–Avon Transition
         21: +0.0,  # Avon
         22: +0.3,  # Avon
         23: +0.0,  # Avon
         24: +0.0,  # Avon
         25: +0.0,  # Avon
         26: +0.0,  # Avon
         27: +0.9,  # Avon
         28: +1.2,  # Avon
         29: +2.3,  # Avon
         30: +2.9,  # Avon
         31: +3.1,  # Avon
         32: +2.7,  # Mid-island
         33: +1.4,  # Mid-island
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
         44: +0.3,  # Mid-island
         45: +0.0,  # Mid-island
         46: +0.0,  # Mid-island
         47: +0.0,  # Mid-island
         48: -0.4,  # Mid-island
         49: -1.1,  # Mid-island
         50: -2.1,  # Mid-island
         51: -2.4,  # Mid-island
         52: -2.2,  # Mid-island
         53: -2.0,  # Mid-island
         54: -1.8,  # Mid-island
         55: -1.3,  # Mid-island
         56: -1.0,  # Mid-island
         57: -0.3,  # Mid-island
         58: +0.0,  # Mid-island
         59: +0.0,  # Mid-island
         60: +0.0,  # Wimble Shoals Influence
         61: +0.0,  # Wimble Shoals Influence
         62: +0.3,  # Wimble Shoals Influence
         63: +0.0,  # Wimble Shoals Influence
         64: +0.0,  # Wimble Shoals Influence
         65: +0.0,  # Wimble Shoals Influence
         66: +0.0,  # Wimble Shoals Influence
         67: +0.0,  # Wimble Shoals Influence
         68: +0.5,  # Wimble Shoals Influence
         69: +1.6,  # Wimble Shoals Influence
         70: +2.2,  # Wimble Shoals Influence
         71: +2.6,  # Wimble Shoals Influence
         72: +2.7,  # Wimble Shoals Influence
         73: +3.0,  # Wimble Shoals Influence
         74: +2.4,  # Wimble Shoals Influence
         75: +1.6,  # Tri-Village / Rodanthe
         76: +0.0,  # Tri-Village / Rodanthe
         77: +0.0,  # Tri-Village / Rodanthe
         78: -1.2,  # Tri-Village / Rodanthe
         79: -2.1,  # Tri-Village / Rodanthe
         80: -3.4,  # Tri-Village / Rodanthe
         81: -3.7,  # Tri-Village / Rodanthe
         82: -4.0,  # Tri-Village / Rodanthe
         83: -4.3,  # Tri-Village / Rodanthe
         84: -4.5,  # Pea Island NWR
         85: -4.3,  # Pea Island NWR
         86: -3.8,  # Pea Island NWR
         87: -2.9,  # Pea Island NWR
         88: -2.2,  # Pea Island NWR
         89: -1.3,  # Pea Island NWR
         90: +26.8,  # LOCKED — end domain, LRR-solved; see the end-domain note above
    },
    2004: {
          1: +50.3,  # LOCKED — end domain, LRR-solved; see the end-domain note above
          2: +0.0,  # Cape Point / Shoal Dynamics
          3: +0.0,  # Cape Point / Shoal Dynamics
          4: +0.0,  # Cape Point / Shoal Dynamics
          5: +0.0,  # Cape Point / Shoal Dynamics
          6: +0.0,  # Cape Point / Shoal Dynamics
          7: +0.0,  # Cape Point / Shoal Dynamics
          8: -0.8,  # Cape Point / Shoal Dynamics
          9: +0.8,  # Cape Point / Shoal Dynamics
         10: +1.6,  # Cape Point / Shoal Dynamics
         11: +1.7,  # Buxton–Avon Transition
         12: +1.9,  # Buxton–Avon Transition
         13: +2.1,  # Buxton–Avon Transition
         14: +2.1,  # Buxton–Avon Transition
         15: +2.2,  # Buxton–Avon Transition
         16: +2.4,  # Buxton–Avon Transition
         17: +2.5,  # Buxton–Avon Transition
         18: +2.4,  # Buxton–Avon Transition
         19: +2.1,  # Buxton–Avon Transition
         20: +1.7,  # Buxton–Avon Transition
         21: +0.3,  # Avon
         22: -0.9,  # Avon
         23: +0.0,  # Avon
         24: +0.0,  # Avon
         25: +0.0,  # Avon
         26: +0.0,  # Avon
         27: -0.2,  # Avon
         28: +1.2,  # Avon
         29: +2.3,  # Avon
         30: +2.9,  # Avon
         31: +3.3,  # Avon
         32: +3.6,  # Mid-island
         33: +3.5,  # Mid-island
         34: +3.5,  # Mid-island
         35: +3.3,  # Mid-island
         36: +2.9,  # Mid-island
         37: +2.6,  # Mid-island
         38: +2.3,  # Mid-island
         39: +2.3,  # Mid-island
         40: +2.0,  # Mid-island
         41: +1.7,  # Mid-island
         42: +1.4,  # Mid-island
         43: +1.0,  # Mid-island
         44: +0.3,  # Mid-island
         45: +0.0,  # Mid-island
         46: +0.0,  # Mid-island
         47: +0.0,  # Mid-island
         48: -0.4,  # Mid-island
         49: +0.0,  # Mid-island
         50: -1.3,  # Mid-island
         51: -1.6,  # Mid-island
         52: -2.2,  # Mid-island
         53: -2.0,  # Mid-island
         54: -1.8,  # Mid-island
         55: -1.3,  # Mid-island
         56: +0.0,  # Mid-island
         57: -0.3,  # Mid-island
         58: +0.0,  # Mid-island
         59: +0.0,  # Mid-island
         60: +0.0,  # Wimble Shoals Influence
         61: +0.0,  # Wimble Shoals Influence
         62: +0.3,  # Wimble Shoals Influence
         63: +1.0,  # Wimble Shoals Influence
         64: +1.5,  # Wimble Shoals Influence
         65: +1.9,  # Wimble Shoals Influence
         66: +2.3,  # Wimble Shoals Influence
         67: +2.6,  # Wimble Shoals Influence
         68: +3.0,  # Wimble Shoals Influence
         69: +3.8,  # Wimble Shoals Influence
         70: +4.1,  # Wimble Shoals Influence
         71: +4.3,  # Wimble Shoals Influence
         72: +4.1,  # Wimble Shoals Influence
         73: +3.3,  # Wimble Shoals Influence
         74: +2.4,  # Wimble Shoals Influence
         75: +2.0,  # Tri-Village / Rodanthe
         76: +1.5,  # Tri-Village / Rodanthe
         77: +1.2,  # Tri-Village / Rodanthe
         78: +0.9,  # Tri-Village / Rodanthe
         79: +0.4,  # Tri-Village / Rodanthe
         80: +0.0,  # Tri-Village / Rodanthe
         81: +0.0,  # Tri-Village / Rodanthe
         82: +0.0,  # Tri-Village / Rodanthe
         83: -2.1,  # Tri-Village / Rodanthe
         84: -2.7,  # Pea Island NWR
         85: -2.9,  # Pea Island NWR
         86: -2.7,  # Pea Island NWR
         87: -2.9,  # Pea Island NWR
         88: -2.2,  # Pea Island NWR
         89: -1.3,  # Pea Island NWR
         90: +52.1,  # LOCKED — end domain, LRR-solved; see the end-domain note above
    },
}

# The 2004 calibrated preset falls back to the 1984 fit if it was never
# solved separately. It has been, so this stays False -- but the flag is
# what tells you whether a "calibrated 2004" run is really calibrated.
if HATTERAS_BE_RATES_CALIBRATED.get(2004) is None:
    HATTERAS_BE_RATES_CALIBRATED[2004] = dict(HATTERAS_BE_RATES_CALIBRATED[1984])
    HATTERAS_BE_RATES_2004_IS_PLACEHOLDER = True
else:
    HATTERAS_BE_RATES_2004_IS_PLACEHOLDER = False

# GIS 90 IS THE ONE VALUE THE TWO PRESETS DO NOT SHARE.
#
# GIS 1 is still sliced out of the calibrated fit below, so re-solving it
# updates both presets at once -- that slicing is what stopped the old "base"
# preset drifting into commented-out copies of these numbers, and it still
# holds everywhere it can. GIS 1 can be shared because it lands within
# 0.22 m/yr of target in BOTH presets: nothing is forced near it (the
# calibrated fit is 0.0 at GIS 2-4 in 1984-2004 and at GIS 2-7 in 2004-2024),
# so its neighbourhood is the same either way.
#
# GIS 90 CANNOT. Measured on the 2026-08-23 matrix, one shared value put
# calibBE 1.136 m/yr (1984-2004) and 0.498 m/yr (2004-2024) off target there,
# against 0.003 and 0.013 for edgeBE at the same number. The asymmetry is
# structural, not a fitting error:
#
#   an ISOLATED forced cell is diffused away by its unforced neighbours, so
#   its gain d(LRR)/d(BE) is only about 0.1 and it needs roughly ten times
#   the misfit it closes. A forced ZONE moves together and is not diffused
#   away, so interior corrections run at a gain near 1.
#
# GIS 90 sits on that join. In edgeBE it is isolated and the 0.1 gain holds.
# In calibBE the fit puts a coherent erosive zone at GIS 84-89 (down to
# -3.0 m/yr in 1984-2004) hard against it, and that zone drags GIS 90 down
# by more than a metre a year -- an effect the per-domain LOESS residual
# never sees, because it fits each domain against a base run in which those
# neighbours were unforced. So the two presets genuinely need different
# numbers at GIS 90, and forcing one on them would mean one of them is
# always wrong there.
HATTERAS_BE_EDGE_SHARED_DOMAINS = (1,)
HATTERAS_BE_EDGE_SPLIT_DOMAINS = (90,)

# GIS 90 under edgeBE, solved on the edgeBE road_bdm base run. The calibBE
# value for the same domain lives in HATTERAS_BE_RATES_CALIBRATED.
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

for _period, _rates in HATTERAS_BE_RATES_EDGE.items():
    _absent = [gis for gis in HATTERAS_BE_EDGE_DOMAINS if gis not in _rates]
    if _absent:
        raise ValueError(
            f"edgeBE {_period}: end domains {_absent} are not in the "
            f"calibrated preset, so the edge preset would silently model them "
            f"at 0.0 m/yr. Add them to HATTERAS_BE_RATES_CALIBRATED.")

    # An edge domain present but ZERO is the same failure the check above
    # exists to prevent, and the absence test does not catch it: edgeBE would
    # be byte-identical to zeroBE at that domain, which makes the two presets
    # indistinguishable there and collapses the edge secant to 0/0 with no
    # error. Found on 2026-08-28 while zeroing the interior of the calibrated
    # table -- GIS 1 is SLICED from that table, so emptying it would have
    # silently disarmed the edge preset.
    _zeroed = [gis for gis in HATTERAS_BE_EDGE_DOMAINS
               if gis in _rates and _rates[gis] == 0.0]
    if _zeroed:
        raise ValueError(
            f"edgeBE {_period}: end domains {_zeroed} are present but 0.0, so "
            f"edgeBE is identical to zeroBE there and the edge secant has no "
            f"second bracket. Give them a nonzero starting value, or run "
            f"zeroBE if that is what you actually want.")

# The split is only defensible while each side is actually solved against its
# own preset. If a period ever gains an edgeBE D90 with no calibBE counterpart
# -- or the two silently converge back to one number, which would mean the
# split has stopped earning its complexity -- say so rather than let it pass.
for _period in HATTERAS_BE_RATES_CALIBRATED:
    if _period not in HATTERAS_BE_EDGE_D90:
        raise ValueError(
            f"edgeBE {_period}: GIS 90 is a SPLIT domain but has no entry in "
            f"HATTERAS_BE_EDGE_D90, so the edge preset has no value solved "
            f"against its own base run.")

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
# NOT period-dependent, and that is deliberate -- but the REASON changed on
# 2026-08-26 and the old one is worth recording because it is now false.
#
# It used to read: "there is one topography (2009) for every period, so one
# elevation set serves both". There are two products now, and they do NOT agree
# under the road. Sampling the 2004 alignment in a 3.5 m corridor on both:
#
#     2009-2014-1996 minus 2009-2014, corridor mean:  median +0.222 m
#     54 of 82 domains move more than 0.05 m; cell counts identical
#
# Identical cell counts means ALACE REPLACED measured 2009 pavement rather
# than filling holes in it. And +0.222 m is not a road: it is the island-wide
# 1996-vs-2009 survey offset, which mosaic_1984_audit.csv reports per domain at
# median +0.255 m (p10 +0.14, p90 +0.33) and HAT_dem_1984_mosaic.py leaves
# UNCORRECTED on purpose ("bias correction OFF, feathering OFF").
#
# So a per-period road elevation built from each period's own DEM would push
# 1984 road_ele up ~0.22 m island-wide, and that increment would be the survey
# offset, not a measured roadbed. A higher road is buried by overwash less
# often, so it would reach the model. One file, built on the 2009-2014
# baseline, keeps that offset out of a forcing.
#
# Hence still no year in the filename -- a per-year name would imply a measured
# change in roadbed height that no data supports.
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
# Relocations carry a DISPLACEMENT, not an absolute setback. CASCADE adds it to
# whatever setback the model is carrying at the event year, which already
# reflects the modelled shoreline retreat -- so the retreat is counted once. An
# absolute setback referenced to the 1984 dune line would count it twice,
# because the topography's own dune line has already moved landward by that
# amount.
#
# WHERE THE DISPLACEMENTS COME FROM -- measured, not typed.
#
# Until 2026-08-20 these were eleven hand-entered literals attributed to a
# 1978->1997 cross-shore offset digitised in ArcGIS Pro. The 1997 line is not in
# the repo -- only nc12_1984.geojson and nc12_2004.geojson are -- so those
# numbers could not be re-derived, checked, or corrected. They are now read from
# HAT_road_relocation_distance.py's per-domain measurement of the two lines that
# ARE on disk.
#
# WHICH COLUMN, AND WHY THE SIGNED ONE
#   mean_relocation_m        unsigned nearest-distance old line -> new line.
#   mean_signed_landward_m   the same displacement vectors projected onto the
#                            landward normal. THIS ONE.
# The unsigned column is dragged toward zero wherever the two lines share
# vertices: the 2004 line was digitised by editing a copy of the 1984 one, so an
# unedited stretch contributes an exact 0.000 m that is an editing artefact, not
# a measurement. On GIS 9 and 87 that is 23% and 33% of samples, and the signed
# mean is correspondingly LARGER (18.0 vs 13.8 m, 17.5 vs 11.8 m). Sign also
# matters on its own terms: a relocation is landward by definition, so a column
# that cannot express direction cannot contradict that claim.
#
# THE VINTAGE GAP. The lines were digitised off 1978 and 2008 imagery, so the
# measured interval brackets BOTH events rather than either one. That is safe
# only because the two events are disjoint in space -- 1989 moves GIS 84-87,
# 1999 moves GIS 9-15 -- so no domain's displacement is claimed twice. Do not
# add a third relocation overlapping either span without revisiting this.

# Per-domain measurement written by
# scripts/input_prep/4-mgmt-forcings/road_relocation/HAT_road_relocation_distance.py
_RELOCATION_MEASUREMENT_FILE = ("4-mgmt-forcing/road_relocation/1984_2004/"
                                "road_relocation_1984_2004.csv")

# WHY THE 1999 EVENT STOPS AT GIS 14. The measurement classifies GIS 15
# 'relocated', but cannot say by how much or in which direction: the two
# digitised lines CROSS inside that domain, so 56% of samples read landward and
# 44% seaward (sign_agreement 0.56), and the mean (+1.8 m) and median (-1.0 m)
# disagree in SIGN -- that mean is the residue of two opposing populations, not
# a displacement. Its mean magnitude, 4.0 m, is below the 5 m re-digitising
# threshold; only its 12.2 m maximum is above, which is the sole reason it
# classified 'relocated' rather than 'redigitized' at all. GIS 15 is therefore
# left out of the 1999 event entirely rather than carried with a forced 0.0 m.


def _measured_displacements(gis_domains):
    """Reads the measured landward displacement for one relocation event.

    Args:
        gis_domains: The GIS domains the event moves. Hand-specified per event,
            NOT taken from the file: the measurement classifies ~40 domains
            'relocated' across the whole island, and which of them belong to
            which historical event is a fact about NC-12, not about the lines.

    Returns:
        A {gis: displacement_m} dict, holding mean_signed_landward_m from the
        measurement.

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

    return {gis: float(measured[gis]["mean_signed_landward_m"])
            for gis in gis_domains}


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

# Independent check on a relocated setback: the MEASURED 2004 setback for the
# same domain. Both relocation events precede 2004, so a correctly displaced
# setback should land near the 2004 same-year measurement.
# Reporting only -- nothing reads this to decide anything.
#
# READ FROM THE FILE PERIOD 2 ACTUALLY RUNS, not typed out here. This used to
# hold eleven literals copied from the LEGACY old-method RoadSetback_2004.csv
# (GIS 9 = 89 m, 10 = 83 m, 11 = 81 m ...). That copy survived the move to the
# dune-start method on 2026-08-18 and both topography bumps after it, so the
# printed check was comparing a dune-start setback against a number measured to
# a different feature -- median +23 m landward, and on a topography the run no
# longer uses. The same domains read 50 / 40 / 20 m in the file period 2
# actually loads. A check that disagrees by construction is worse than none:
# it invites a correct relocation to be read as a failed one.
#
# Deriving it has the same motive as slicing HATTERAS_BE_RATES_EDGE out of the
# calibrated preset -- the two can no longer drift apart, and bumping the
# topography version moves this with the setbacks it is checking.


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
# VOLUMES AND EXTENTS ARE FROM THE PROJECT RECORD, and both halves matter:
# CASCADE applies volume/length, so a right volume on a short footprint is as
# wrong as a wrong volume. Checked against the record 2026-08-22; Buxton was
# already correct, Avon and Rodanthe were not.
#
#   project    was              now              record
#   Rodanthe   619 m^3/m        408 m^3/m        ~380 m^3/m
#   Buxton     183 m^3/m        183 m^3/m        ~197 m^3/m   (unchanged)
#   Avon       841 m^3/m        191 m^3/m        190-216 m^3/m
#
# Avon was wrong on BOTH axes -- 2.2x the volume placed, over 55% of the
# footprint -- which compounded to 3.9-4.4x the real fill density and a 68 m
# instantaneous shoreline step. That step is what BRIE's Crank-Nicolson solve
# rings on; see compute_lrr and HAT_hindcast_methods.md section 12.
HATTERAS_NOURISHMENT_PROJECTS = (
    NourishmentProject(
        name="Rodanthe emergency fill",
        year=2014,
        # USACE/NCDOT emergency fill at the Mirlo Beach "S-curves", ~2 miles
        # (3.2 km) immediately NORTH of Rodanthe village, whose community zone
        # ends at GIS 83. Six domains is 3.0 km (1.86 mi) rather than the
        # 6.4 domains 2 miles would take: GIS 90 is the locked north-end
        # domain carrying the source/sink boundary value (tens of m/yr), and
        # putting fill into a domain whose rate was pinned rather than
        # modelled would make the fill and the boundary condition
        # indistinguishable. The 7% shortfall is the cost of that exclusion
        # and is stated rather than absorbed into the volume.
        gis_domains=tuple(range(84, 90)),
        volume_cubic_yards=1_600_000,
        note="Mirlo Beach S-curves, ~2 mi N of Rodanthe; GIS 84-89 carry NC-12",
    ),
    NourishmentProject(
        name="Buxton shore protection",
        year=2022,
        # 2.9 mi (4.7 km) north from the oceanfront groin at the lighthouse,
        # which sits at GIS 5.5 -- so 6-15, which is what was already here.
        # Left untouched: it is the one project whose configured density
        # (183 m^3/m) already matched the record (197 m^3/m).
        gis_domains=tuple(range(6, 16)),
        volume_cubic_yards=1_200_000,
        note="Extends north out of Buxton village (7-8) into the road corridor",
    ),
    NourishmentProject(
        name="Avon shore protection",
        year=2022,
        # Placed from 3,000 ft north of Avon Pier (Due East Road) south to
        # Askins Creek North Drive at the village's southern boundary. The
        # pier is GIS 26 (HATTERAS_ANNOTATIONS.piers), 3,000 ft is 1.8
        # domains, and the southern village boundary is GIS 21 -- so 21-28,
        # 4.0 km, matching the ~2.5 mi the record gives. The previous 23-26
        # was the MIDDLE of this footprint, missing ~2 km to the south and
        # ~1 km to the north. Still entirely inside the Avon community zone
        # (21-31), so the beach_dune_manager footprint is unchanged by this.
        gis_domains=tuple(range(21, 29)),
        volume_cubic_yards=1_000_000,
        note="Due East Rd to Askins Creek N Dr; inside the Avon zone (21-31)",
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
