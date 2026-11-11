"""
Dune toes from the lidar, and a collision-regime dune erosion rule calibrated on the 1996-2009 surveys.

    python scripts/real-world-cases/Hatteras/dune-collision/4-hindcast/1-calibration-1996-2009/HAT_dune_collision.py toes
    python scripts/real-world-cases/Hatteras/dune-collision/4-hindcast/1-calibration-1996-2009/HAT_dune_collision.py twl
    python scripts/real-world-cases/Hatteras/dune-collision/4-hindcast/1-calibration-1996-2009/HAT_dune_collision.py calibrate
    python scripts/real-world-cases/Hatteras/dune-collision/4-hindcast/1-calibration-1996-2009/HAT_dune_collision.py empirical
    python scripts/real-world-cases/Hatteras/dune-collision/4-hindcast/1-calibration-1996-2009/HAT_dune_collision.py run
    python scripts/real-world-cases/Hatteras/dune-collision/4-hindcast/1-calibration-1996-2009/HAT_dune_collision.py score
    python scripts/real-world-cases/Hatteras/dune-collision/4-hindcast/1-calibration-1996-2009/HAT_dune_collision.py plot

Barrier3D lowers a dune only when Rhigh overtops it. Sallenger (2000) puts a storm in the collision regime
when the water reaches the dune toe but not the crest (D_low < R_high < D_high), where the real dunes lose
sand to scarping. The magnitude here is the wave-impact rule of Larson et al. (2004),
dV = 4 Cs (R - z_toe)^2 t / T, summed hour by hour over the hours the total water level stands above the
toe, and taken as a crest lowering over the model's two-cell (20 m) dune: dH = dV / 20 m.

toes       per survey, domain and 10 m cell, on the 1996-start dune windows: the crest as
           HAT_survey_dune_crests.py finds it, the beach start (first cell above 0.5 m MHW) and the toe,
           the point between them farthest below the straight line from the beach start to the crest
           (the maximum-distance-below-line definition). Also the start's own arrays, the toe the model
           would use. Writes dune_toes.csv
twl        the hourly total water level the storm series is built from (Duck gauge + Stockdon R2% on WIS
           ST63228, slope 0.06, as historical_storm_creation_v3_HAT.py), checked against the series' peaks.
           Writes twl_hourly_1996_2009.csv
calibrate  for each pair of consecutive surveys, the cells the water never overtopped between them: the
           observed crest change against the exposure sum((TWL - toe)^2 dt / Tp) from the cell's toe (the
           median over the surveys that pass the toe screen), and the Cs that fits. Writes
           collision_calibration_cells.csv, collision_calibration.csv and collision_calibration.png
empirical  the same cells in the terms of Long et al. (2014): crest lowering over pre-storm crest (dz/D_high)
           against R_high/D_high and against beach width, and the fit dz/D_high = a max(0, 1 - W / W0) on the
           storm intervals. Writes collision_empirical.csv, collision_empirical_fit.csv and collision_empirical.png
run        the managed edgeBE 1996-2009 hindcast with the groin, as the headline matrix run, with collision
           erosion added after each model year, in-process only: every dune cell loses
           Cs * 4 sum((TWL - toe)^2 dt / Tp) / 20 m over that calendar year's hours with the water between its
           toe and its crest, down to no lower than the toe; both dune rows scaled together, as the forced
           crests are. Optionally the eroded sand goes onto the beach, moving the shoreline seaward as
           CASCADE's shoreface nourishment does (dx_s = 2 V / (2 h_b + D_sf)), passed on to BRIE:
               larson_lost     Cs 1.7e-4 (Dennis fit, = Larson's US East Coast storm value), sand lost
               larson_beach    Cs 1.7e-4, sand to the beach
               strong_lost     Cs 4e-4 (what the Dennis hotspots, GIS 11 and 86, need), sand lost
               strong_beach    Cs 4e-4, sand to the beach
           or by the beach-width rule from empirical: once in each calendar year whose highest water reaches
           the cell's toe but not its crest, the crest loses a max(0, 1 - W / W0) of its height (m NAVD88),
           W the cell's beach width in the last survey flown before that year's highest water:
               width_lost          sand lost
               width_beach         sand to the beach
               width_storms_beach  sand to the beach, only in years whose highest water reaches 3.1 m NAVD88,
                                   the lowest R_high of the two storm intervals the rule is fit on
           and the dune growth rate refit under width_storms_beach, which takes the 2009 median crest 1.2 m
           below the lidar's at the calibrated 0.01-0.03: the same rule with Barrier3D's growth bounds on
           every domain at mean r +- 0.01:
               width_storms_r04 ... width_storms_r10   mean r 0.04, 0.06, 0.08, 0.10
           or, process-based, by cascade/dune_collision.py inside Barrier3D's storm loop: each storm lowers the
           cells it leaves in collision (toe < its peak Rhigh < crest) by Larson's rule per storm,
           Cs 4 (Rhigh - toe)^2 duration / T / 20 m, before the next storm; nothing is changed after the year.
           Cs from calibrate's per-storm fit (island 2.0e-5, Dennis 4.2e-5); the eroded sand lost:
               storm_cs2e5, storm_cs4e5               the calibrated growth rate, 0.01-0.03
               storm_cs2e5_r06, storm_cs4e5_r06       growth 0.05-0.07
           A threshold rule (no loss below an exposure, steeper above) fit Dennis no better than the linear
           one (r 0.45 against 0.50 on domain medians; best threshold 0), so it is not run. The control is
           the forced-dune-crests experiment's control run (same code, nothing added)
score      island skill, per-GIS LRR against CoastSat, and the model's crests against the lidar's in each
           survey's model state: the RMSE of the island's median crest over the surveys, and the stable
           dunes the growth rate was calibrated on (GIS 20-22, 59-62) in 2009. Writes collision_runs_summary.csv, collision_runs_by_gis.csv and
           collision_runs_crests.csv
plot       every run's end-minus-start shoreline position change against the CoastSat net change, in the
           layout of the run's shoreline_position_change figure, for the island, Buxton and Rodanthe; the
           runs named after it (and the control; the files then carry their names), or every run.
           Writes collision_runs_position_<view>.png

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-10-26
"""
from __future__ import annotations

import argparse
import json
import os
import subprocess
import sys
import time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np
import pandas as pd

_HERE = Path(__file__).resolve()
# Hatteras/ holds the scripts; data/ and output/ live in the workspace on D: (hatteras_paths)
_HAT_DIR = next(_p for _p in _HERE.parents if (_p / "hatteras_paths.py").exists())
if str(_HAT_DIR) not in sys.path:
    sys.path.insert(0, str(_HAT_DIR))
from hatteras_paths import WORKSPACE as PROJECT_ROOT  # noqa: E402

# --- CONFIG ------------------------------------------------------------------
TAG = "2026-10-26-dune-collision"
EXP_DIR = PROJECT_ROOT / "output" / "raw_runs" / "dune-storm-response" / TAG
START = "1996-start"
SURVEY_CRESTS = (PROJECT_ROOT / "data" / "hatteras_init" / "1-barrier3d-domains" / START
                 / "survey-crests" / "survey_dune_crests.csv")
STORM_SUMMARY = (PROJECT_ROOT / "data" / "hatteras_init" / "3-env-forcings" / "3-storms" / "hindcast_storms"
                 / "1997_2009" / "1997_2009_storms_v3_240_summary.csv")
MHW_NAVD88_M = 0.36
TWL_START, TWL_END = "1996-01-01 00:00:00", "2009-12-31 23:00:00"
BEACH_SLOPE = 0.06               # as the storm series
DUNE_WIDTH_M = 20.0              # Barrier3D's two 10 m dune cells
MIN_BEACH_CELLS = 2              # toe needs a beach start at least this many cells seaward of the crest
# The calibration's toe: per cell, the median over the surveys of the toes that pass the screen. Not each
# survey's own: the 1998 grid puts 38% of its toes on the wet beach (median 1.7 m, some below 0 m NAVD88),
# and it is the survey before Dennis
TOE_SURVEYS_EXCLUDED = (1998,)
TOE_MIN_M_NAVD88 = 1.0           # above the beach start (0.86 m NAVD88)
TOE_BELOW_CREST_M = 0.5          # and at least this far below the crest
# Fit on the storm intervals without the post-Isabel repairs (2001-2004) inside them; all intervals as well
FIT_INTERVALS = ("1996-1997", "1997-1998", "1998-1999", "1999-2000", "2000-2001", "2004-2005", "2005-2009")
# Empirical scaling: fit on the intervals with a storm signal in the crests (Dennis; 2004-05), each less its
# median change on beaches wider than this
EMPIRICAL_INTERVALS = ("1998-1999", "2004-2005")
BACKGROUND_BEACH_WIDTH_M = 60.0
# Model runs
HINDCAST = _HAT_DIR / "hatteras_ms" / "HAT_hindcast.py"
CONTROL_DIR = (PROJECT_ROOT / "output" / "raw_runs" / "dune-storm-response" / "2026-10-26-forced-dune-crests"
               / "runs" / "control")
RUN_START, RUN_END, PRESET, SCENARIO = 1996, 2009, "edgeBE", "full_management"
# The rule, its Cs (Larson) or the lowest year-maximum TWL it acts in (beach width, m NAVD88; None: any year
# the water reaches the toe), whether the eroded sand goes onto the beach, and the dune growth bounds
# (rmin, rmax) on every domain (None: the runner's calibrated 0.01-0.03)
VARIANTS = {"larson_lost": ("larson", 1.7e-4, False, None),
            "larson_beach": ("larson", 1.7e-4, True, None),
            "strong_lost": ("larson", 4e-4, False, None),
            "strong_beach": ("larson", 4e-4, True, None),
            "width_lost": ("width", None, False, None),
            "width_beach": ("width", None, True, None),
            "width_storms_beach": ("width", 3.1, True, None),
            # In-model, storm by storm (cascade/dune_collision.py): Cs from the per-storm calibration, island fit
            # 2.0e-5 or Dennis 4.2e-5; the eroded sand lost
            "storm_cs2e5": ("storm", 2.0e-5, False, None),
            "storm_cs4e5": ("storm", 4.2e-5, False, None),
            "storm_cs2e5_r06": ("storm", 2.0e-5, False, (0.05, 0.07)),
            "storm_cs4e5_r06": ("storm", 4.2e-5, False, (0.05, 0.07)),
            **{f"width_storms_r{round(r * 100):02d}": ("width", 3.1, True, (r - 0.01, r + 0.01))
               for r in (0.04, 0.06, 0.08, 0.10)},
            # Committee test: a source/sink where the model's dunes do not erode, with and without collision
            "sink": (None, None, False, None),
            "storm_cs2e5_r06_sink": ("storm", 2.0e-5, False, (0.05, 0.07)),
            # Per storm with the beach: Larson scaled by max(0, 1 - W / 45 m), W each cell's median beach width
            # over the surveys (1998 left out); Cs refit with the factor (domain medians, r 0.48, Dennis 0.86)
            "storm_w45": ("storm", 1.24e-4, False, None),
            "storm_w45_r04": ("storm", 1.24e-4, False, (0.03, 0.05)),
            "storm_w45_r06": ("storm", 1.24e-4, False, (0.05, 0.07)),
            # From 2026-10-28 the parameter file sets StormStart: 1 (the 1997 storms run); every run above had 2
            "control_ss1": (None, None, False, None),
            "storm_w45_r04_ss1": ("storm", 1.24e-4, False, (0.03, 0.05)),
            # Advisors' Hs 3.0 m (CASCADE's wave_height: BRIE's alongshore transport and shoreface depth; the
            # storm series and so the collision forcing are unchanged); the edgeBE ends were solved at 2.0
            "control_hs3": (None, None, False, None),
            "storm_w45_r04_hs3": ("storm", 1.24e-4, False, (0.03, 0.05)),
            # The same on the ends solved at Hs 3.0 (ENDS_HS3)
            "control_hs3_ends": (None, None, False, None),
            "storm_w45_r04_hs3_ends": ("storm", 1.24e-4, False, (0.03, 0.05)),
            # Hs 2.5 on the Hs 2.0 ends: is 2.5 better than 3.0 (the ends moved the Hs 3.0 scores <= 0.07 m/yr)
            "control_hs2p5": (None, None, False, None),
            "storm_w45_r04_hs2p5": ("storm", 1.24e-4, False, (0.03, 0.05)),
            # Groin strength sweep at Hs 2.5 and 3.0 (the 0.10 runs are storm_w45_r04_hs2p5 / _hs3_ends)
            **{f"storm_w45_r04_{h}_g{round(g * 100):02d}": ("storm", 1.24e-4, False, (0.03, 0.05))
               for h in ("hs2p5", "hs3") for g in (0.04, 0.06, 0.08, 0.13, 0.16, 0.19)},
            # Cs 1.45e-4: the per-storm beach fit with the module's toe floor (HAT_collision_rule_comparison.py,
            # 2026-10-29); 1.24e-4 above was fit without it
            "storm_cs145_r04_hs3_ends": ("storm", 1.45e-4, False, (0.03, 0.05)),
            "storm_cs145_r04_hs2p5": ("storm", 1.45e-4, False, (0.03, 0.05)),
            # Calibration on the 2026-10-29 defaults (village dune cap 8.0 m NAVD88, Hs 3.0): a new control, and
            # the module (Cs 1.24e-4, W0 45 m) over the growth rates, Hs 3.0 on ENDS_HS3; compare within this set
            "cal_control": (None, None, False, None),
            **{f"cal_r{round(r * 100):02d}": ("storm", 1.24e-4, False, (round(r - 0.01, 3), round(r + 0.01, 3)))
               for r in (0.02, 0.03, 0.04, 0.05, 0.06)},
            # Zoned growth (HAT_collision_calibration.py growth): 0.02-0.04 everywhere but the Buxton-Avon road
            # (GIS 8-20), which still sat 0.35 m above the lidar at 0.01-0.03
            **{f"cal_zone_bx{lo * 100:.0f}{hi * 100:.0f}": ("storm", 1.24e-4, False,
                                              {g: ((lo, hi) if 8 <= g <= 20 else (0.02, 0.04)) for g in range(1, 91)})
               for lo, hi in ((0.0, 0.02), (0.0, 0.01))},
            # Hs sweep: cal_r03 (growth 0.02-0.04) at Hs 2.0-4.0, each on its own ends (ENDS_BY_HS); cal_r03 is Hs 3.0
            **{f"cal_hs{h:g}".replace(".", "p"): ("storm", 1.24e-4, False, (0.02, 0.04)) for h in (2.0, 2.5, 3.5, 4.0)},
            # Groin strength sweep on cal_r03 (Hs 3.0, growth 0.02-0.04; 0.10 is cal_r03 itself)
            **{f"cal_g{round(g * 100):02d}": ("storm", 1.24e-4, False, (0.02, 0.04)) for g in (0.13, 0.16, 0.19, 0.22)},
            # The control at the chosen settings (Hs 3.0, groin 0.16), against cal_g16
            "cal_control_g16": (None, None, False, None),
            # A quick look at Hs 2.75 on interpolated ends (ENDS_BY_HS), groin 0.16
            "cal_hs2p75": ("storm", 1.24e-4, False, (0.02, 0.04)),
            # cal_g16 with the module off at Cape Point, south of the Buxton groin (GIS 1-5), where it lowers
            # the lidar's crests too far and the toes cannot be picked reliably (berm and runnel)
            "cal_g16_nocape": ("storm", 1.24e-4, False, (0.02, 0.04)),
            # cal_g16_nocape with a dune cap per village, then also no bulldozed overwash on Rodanthe's dunes
            "cal_g16_nocape_vcap": ("storm", 1.24e-4, False, (0.02, 0.04)),
            "cal_g16_nocape_vcap_ow0": ("storm", 1.24e-4, False, (0.02, 0.04)),
            # cal_g16_nocape on toes and beach widths from the 1996 start profile alone (TOES_FROM_START), not the
            # median over the 1996-2009 surveys: a Cs sweep, the start's toes taking nothing from after it
            **{f"cal_g16_nocape_ts{round(cs * 1e6):d}": ("storm", cs, False, (0.02, 0.04))
               for cs in (1.24e-4, 2.5e-4, 5e-4, 1e-3)},
            # cal_g16_nocape with the beach following the modelled shoreline (DuneCollisionManager
            # beach_follows_shoreline): from the median widths, then from the 1996 start profile's widths
            "cal_g16_nocape_dynw": ("storm", 1.24e-4, False, (0.02, 0.04)),
            "cal_g16_nocape_dynw_sw": ("storm", 1.24e-4, False, (0.02, 0.04)),
            # cal_g16_nocape with the road's dune rebuild set from the lidar rebuilds (ROAD_DUNE_REBUILD): the
            # setback gate alone, then with the rebuild trigger height too
            "cal_g16_nocape_rd150": ("storm", 1.24e-4, False, (0.02, 0.04)),
            "cal_g16_nocape_rd150h": ("storm", 1.24e-4, False, (0.02, 0.04)),
            # ... and with a 198 m shoreline-to-road gate, the FEIS's 650 ft monitoring warning sign, not an NCDOT
            # dune maintenance rule (SHORELINE_GATE)
            "cal_g16_nocape_rd198": ("storm", 1.24e-4, False, (0.02, 0.04)),
            # ... with the beach acting on the water instead of the volume (RUNUP_LOSS): R_high less 0.02 m per m
            # of beach at the dune, Cs 1.64e-4 its offline fit with the toe floor; and the ramp at its own such fit
            "cal_g16_nocape_rd198_rl02": ("storm", 1.64e-4, False, (0.02, 0.04)),
            "cal_g16_nocape_rd198_cs145": ("storm", 1.45e-4, False, (0.02, 0.04)),
            # ... with the beach factor on the dry beach (beach start to toe, 2026-11-03) at W0 25 m: Cs 1.02e-4 its
            # offline fit with the toe floor, and 0.87e-4 (scaled as the crest-based 1.45e-4 -> 1.24e-4 in-model)
            "cal_g16_nocape_rd198_dry": ("storm", 1.02e-4, False, (0.02, 0.04)),
            # the module-off control with the 198 m road rule too: against cal_g16_nocape_rd198_dry, the module alone
            "cal_control_g16_rd198": (None, None, False, None),
            # The model-change steps (2026-11-04), each adding one change to the one before: the base; + the dune
            # growth paired with the module; + the module (dry beach, W0 25 m, Cs 1.02e-4, off at Cape Point);
            # + the 198 m road rule; + Rodanthe 90% (the default since 2026-11-04); + toes and beach widths
            # from the 1996 start profile. Steps 1-4 hold Rodanthe at 40% (RODANTHE_40)
            # the base at the growth bounds the calibration had before the module (0.01-0.03); the 1996 default
            # has been 0.02-0.04 since 2026-10-29, paired with the module
            "cal_step1_base": (None, None, False, (0.01, 0.03)),
            "cal_step2_growth": (None, None, False, (0.02, 0.04)),
            "cal_step3_module": ("storm", 1.02e-4, False, (0.02, 0.04)),
            "cal_step4_road": ("storm", 1.02e-4, False, (0.02, 0.04)),
            "cal_step5_rodanthe90": ("storm", 1.02e-4, False, (0.02, 0.04)),
            "cal_step6_starttoes": ("storm", 1.02e-4, False, (0.02, 0.04)),
            # cal_step5_rodanthe90 with one toe per domain, the mean of its screened cells over the two surveys
            # nearest the start, 1996 and 1997 (TOE_DOMAIN_AVG; 2026-11-04)
            "cal_step6_toeavg": ("storm", 1.02e-4, False, (0.02, 0.04)),
            # steps 3 and 4 again on the domain toes, the default since 2026-11-04 (step 5 is cal_step6_toeavg)
            "cal_step3_module_dt": ("storm", 1.02e-4, False, (0.02, 0.04)),
            "cal_step4_road_dt": ("storm", 1.02e-4, False, (0.02, 0.04)),
            "cal_step6_toeavg_cw": ("storm", 1.02e-4, False, (0.02, 0.04)),
            # cal_step6_toeavg with Larson alone: the beach factor 1 everywhere (width 0), Cs 1.9e-5 its offline fit
            # without the factor (rule comparison, larson_storm, re-read crests; 2.0e-5 on the old), 2026-11-04
            "cal_step6_toeavg_la": ("storm", 1.9e-5, False, (0.02, 0.04)),
            "cal_g16_nocape_rd198_dry_cs087": ("storm", 0.87e-4, False, (0.02, 0.04)),
            # The alternative to cal_g16_nocape: the module everywhere, with faster dune recovery at Cape Point
            # (GIS 1-5; 0.10-0.20 is the cape's own rate of the 2026-10-20 calibration, dropped 2026-10-21)
            **{f"cal_g16_capeg{lo * 100:02.0f}": ("storm", 1.24e-4, False,
                                                 {g: ((lo, hi) if g <= 5 else (0.02, 0.04)) for g in range(1, 91)})
               for lo, hi in ((0.05, 0.10), (0.10, 0.20))}}
# The village overwash filter (2026-11-03): the share of each year's overwash the villages (Buxton, Avon, the
# Tri-Village; HATTERAS_COMMUNITY_ZONES) return to the shoreface, one setting for all of them, at points in Rogers et
# al.'s (2015) 40-90% (residential to commercial); 40 is the hindcast's. The 9% bulldozed onto the dune is kept.
# At Rodanthe the village rule drove most of the 2009-2025 shoreline miss (2-test-2009-2025, control_rd198_ow*).
# Chosen here and checked on 2009-2025 (HAT_collision_test_2009.py VILLAGE_FILTER). _vf40 reruns its base on the
# survey windows that follow the dune (toes regenerated the same day), so the three compare
VILLAGE_FILTER = {f"{b}_vf{p}": p for b in ("cal_control_g16", "cal_g16_nocape_rd198_dry") for p in (40, 65, 90)}
VILLAGE_FILTER.pop("cal_control_g16_vf40")       # the module off reads no toes: cal_control_g16 is its base as run
for _v in VILLAGE_FILTER:
    VARIANTS[_v] = VARIANTS[_v.rsplit("_vf", 1)[0]]
# The GIS the module runs on per variant; GIS 1-90 otherwise
COLLISION_GIS = {"cal_g16_nocape": range(6, 91), "cal_g16_nocape_vcap": range(6, 91),
                 "cal_g16_nocape_vcap_ow0": range(6, 91),
                 **{v: range(6, 91) for v in VARIANTS if v.startswith(("cal_g16_nocape_ts", "cal_g16_nocape_dynw",
                                                                         "cal_g16_nocape_rd"))}}
TOES_FROM_START = tuple(v for v in VARIANTS if v.startswith("cal_g16_nocape_ts")) + ("cal_step6_starttoes",)
CAL_STEPS = ("cal_step1_base", "cal_step2_growth", "cal_step3_module", "cal_step4_road", "cal_step5_rodanthe90",
             "cal_step6_starttoes")
COLLISION_GIS.update({v: range(6, 91) for v in CAL_STEPS[2:]})
# Rodanthe village (GIS 77-83) held at the 40% the other villages keep, as before 2026-11-04
RODANTHE_40 = ("cal_step1_base", "cal_step2_growth", "cal_step3_module", "cal_step4_road")
# One toe per domain: the mean of its screened cells (toe above TOE_MIN_M_NAVD88 and TOE_BELOW_CREST_M under the
# crest) over these surveys, every cell of the domain given it
TOE_DOMAIN_AVG = {v: (1996, 1997) for v in ("cal_step6_toeavg", "cal_step3_module_dt", "cal_step4_road_dt",
                                             "cal_step6_toeavg_cw", "cal_step6_toeavg_la")}
LARSON_ALONE = ("cal_step6_toeavg_la",)
# Larson alone's Cs, in the model (2026-11-04): cal_step6_toeavg_la at other Cs, chosen on 1996-2009 alone
LA_CS_SWEEP = {f"cal_la_cs{round(cs * 1e6):d}": cs for cs in (1e-5, 3e-5, 5e-5, 8e-5)}
for _v, _cs in LA_CS_SWEEP.items():
    VARIANTS[_v] = ("storm", _cs, False, (0.02, 0.04))
    TOE_DOMAIN_AVG[_v] = (1996, 1997)
LARSON_ALONE = LARSON_ALONE + tuple(LA_CS_SWEEP)
# The model-change steps with Larson alone, the default since 2026-11-04 (Cs 1.9e-5, no beach factor): steps 3
# and 4 again (Rodanthe 40%); step 5 is cal_step6_toeavg_la
for _v in ("cal_step3_module_la", "cal_step4_road_la"):
    VARIANTS[_v] = ("storm", 1.9e-5, False, (0.02, 0.04))
    TOE_DOMAIN_AVG[_v] = (1996, 1997)
    RODANTHE_40 = RODANTHE_40 + (_v,)
LARSON_ALONE = LARSON_ALONE + ("cal_step3_module_la", "cal_step4_road_la")
# cal_step6_toeavg_la with the 198 m gate as one island-wide road trigger, 198 m less the island's median beach
# (~50 m) = 150 m, the roadway manager's own setting (road_dune_rebuild_trigger_m), in place of the per-domain gates
# (2026-11-04); set through ROAD_DUNE_REBUILD below
VARIANTS["cal_final_rd150u"] = ("storm", 1.9e-5, False, (0.02, 0.04))
TOE_DOMAIN_AVG["cal_final_rd150u"] = (1996, 1997)
LARSON_ALONE = LARSON_ALONE + ("cal_final_rd150u",)
# cal_step3_module_la with the domain toe averaged over every calibration survey (1998 left out, as the toe screen),
# not only the two nearest the start: the toe after Bonnie and Dennis lowered the beach, at the cost of using
# surveys from after the start (2026-11-04, option 1 for Buxton)
VARIANTS["cal_step3_module_la_alltoes"] = ("storm", 1.9e-5, False, (0.02, 0.04))
TOE_DOMAIN_AVG["cal_step3_module_la_alltoes"] = (1996, 1997, 1999, 2000, 2001, 2004, 2005, 2009)
RODANTHE_40 = RODANTHE_40 + ("cal_step3_module_la_alltoes",)
LARSON_ALONE = LARSON_ALONE + ("cal_step3_module_la_alltoes",)
# ... and at stronger Cs: 5e-5, and 1e-4, about the older beach-factor runs' Cs at Buxton's narrow beaches
for _v, _cs in (("cal_step3_module_la_alltoes_cs50", 5e-5), ("cal_step3_module_la_alltoes_cs100", 1e-4),
                 ("cal_step3_module_la_alltoes_cs10", 1e-5), ("cal_step3_module_la_alltoes_cs30", 3e-5),
                 # the hourly offline fit's Cs (larson_hourly), in the per-storm module
                 ("cal_step3_module_la_alltoes_cs110", 1.1e-4)):
    VARIANTS[_v] = ("storm", _cs, False, (0.02, 0.04))
    TOE_DOMAIN_AVG[_v] = TOE_DOMAIN_AVG["cal_step3_module_la_alltoes"]
    RODANTHE_40 = RODANTHE_40 + (_v,)
    LARSON_ALONE = LARSON_ALONE + (_v,)
# cal_step3_module_la_alltoes_cs110 with the later steps: + the 150 m road rule (Rodanthe 40%), then + Rodanthe 90%
for _v, _rod40 in (("cal_alltoes_cs110_road", True), ("cal_alltoes_cs110_final", False)):
    VARIANTS[_v] = ("storm", 1.1e-4, False, (0.02, 0.04))
    TOE_DOMAIN_AVG[_v] = TOE_DOMAIN_AVG["cal_step3_module_la_alltoes"]
    LARSON_ALONE = LARSON_ALONE + (_v,)
    if _rod40:
        RODANTHE_40 = RODANTHE_40 + (_v,)
# The same two Cs on the default toes (1996 and 1997), to set against the all-survey toes (no road rule, Rodanthe 40%,
# as cal_step3_module_la; the cal_la_cs* sweep carries the road rule and Rodanthe 90%)
for _v, _cs in (("cal_step3_module_la_cs50", 5e-5), ("cal_step3_module_la_cs100", 1e-4)):
    VARIANTS[_v] = ("storm", _cs, False, (0.02, 0.04))
    TOE_DOMAIN_AVG[_v] = (1996, 1997)
    RODANTHE_40 = RODANTHE_40 + (_v,)
    LARSON_ALONE = LARSON_ALONE + (_v,)
# The 1996 transect toes (input_prep/9-dune-toes/HAT_toe_method_line.py: one per domain, the median of five
# 100 m transects of the 1996 lidar, front dune inside the model window, no lower than the berm; 2026-11-05), at
# the offline refit's Cs (8e-6), 1e-5 and the current 1.9e-5; and the current toes at 8e-6 to set against them. As
# cal_step3_module_la: Larson alone, dune growth 0.02-0.04, no road rule, Rodanthe 40%
TRANSECT_TOES_FILE = (PROJECT_ROOT / "data" / "hatteras_init" / "9-dune-toes" / "1-transects-1996" / "1-toes"
                      / "domain_toes_1996.csv")
TRANSECT_TOES = ("cal_transect_cs8", "cal_transect_cs10", "cal_transect_cs19")
for _v, _cs in (("cal_transect_cs8", 8e-6), ("cal_transect_cs10", 1e-5), ("cal_transect_cs19", 1.9e-5),
                ("cal_step3_module_la_cs8", 8e-6)):
    VARIANTS[_v] = ("storm", _cs, False, (0.02, 0.04))
    TOE_DOMAIN_AVG[_v] = (1996, 1997)            # overridden by the transect toes for TRANSECT_TOES
    RODANTHE_40 = RODANTHE_40 + (_v,)
    LARSON_ALONE = LARSON_ALONE + (_v,)
# The hand-checked 1996 Automorph toes (input_prep/9-dune-toes/HAT_toe_automorph_final.py, 2026-11-09), one per domain,
# at about their offline refit's Cs (9.7e-6) and at the current 1.9e-5: otherwise cal_step6_toeavg_la (road rule, Rodanthe 90%)
AUTOMORPH_TOES_FILE = TRANSECT_TOES_FILE.with_name("automorph_domain_toes_1996.csv")
# cal_automorph_module_cs10: the module alone, cal_step2_growth with the module on (no road rule, Rodanthe 40%)
AUTOMORPH_TOES = ("cal_automorph_cs10", "cal_automorph_cs19", "cal_automorph_module_cs10")
for _v, _cs in (("cal_automorph_cs10", 1e-5), ("cal_automorph_cs19", 1.9e-5), ("cal_automorph_module_cs10", 1e-5)):
    VARIANTS[_v] = ("storm", _cs, False, (0.02, 0.04))
    TOE_DOMAIN_AVG[_v] = (1996, 1997)            # overridden by the Automorph toes
    LARSON_ALONE = LARSON_ALONE + (_v,)
RODANTHE_40 = RODANTHE_40 + ("cal_automorph_module_cs10",)
# cal_step2_growth and cal_automorph_module_cs10 with the road's first bulldozing kept off the dunes (2026-11-09): in
# year 1 it clears the starting topography above road grade, up to 17 m3/m, and put that on the dune
_NOSCRAPE_CS = {f"cal_automorph_module_cs{_c}_noscrape": _c * 1e-6 for _c in (10, 30, 50, 100, 200)}
NO_FIRST_SCRAPE = ("cal_step2_growth_noscrape", *_NOSCRAPE_CS)
VARIANTS["cal_step2_growth_noscrape"] = (None, None, False, (0.02, 0.04))
for _v, _cs in _NOSCRAPE_CS.items():
    VARIANTS[_v] = ("storm", _cs, False, (0.02, 0.04))
    TOE_DOMAIN_AVG[_v] = (1996, 1997)            # overridden by the Automorph toes
    LARSON_ALONE = LARSON_ALONE + (_v,)
RODANTHE_40 = RODANTHE_40 + NO_FIRST_SCRAPE
AUTOMORPH_TOES = AUTOMORPH_TOES + tuple(_NOSCRAPE_CS)
# The toe at the model's own berm (Barrier3D's dune stands on it; Laura, 2026-11-09), every cell, with the road fix, at
# lower Cs: the water above the toe is about 1-2 m more than with the lidar toes, and it enters squared
TOE_AT_BERM_M = 1.7                              # Barrier3D's berm, every domain [m NAVD88]
TOE_AT_BERM = {f"cal_toeberm_cs{_c}": _c * 1e-6 for _c in (1, 2, 5, 10, 30, 100)}
TOE_AT_BERM.update({"cal_toeberm_cs0p1": 1e-7, "cal_toeberm_cs0p3": 3e-7})   # the low end of the sensitivity sweep
# The same with lateral erosion (COLLISION_ROWS below; 2026-11-11): the front row keeps its height until the face has
# retreated through it, so the Cs needed is higher
TOE_AT_BERM.update({f"cal_toeberm_lat_cs{_c}": _c * 1e-6 for _c in (1, 3, 10, 20, 30, 50, 70, 100, 150, 200, 250)})
for _v, _cs in TOE_AT_BERM.items():
    VARIANTS[_v] = ("storm", _cs, False, (0.02, 0.04))
    TOE_DOMAIN_AVG[_v] = (1996, 1997)            # overridden by the berm
    LARSON_ALONE = LARSON_ALONE + (_v,)
RODANTHE_40 = RODANTHE_40 + tuple(TOE_AT_BERM)
NO_FIRST_SCRAPE = NO_FIRST_SCRAPE + tuple(TOE_AT_BERM)
# Which dune rows lose the collision volume, per variant (cascade/dune_collision.lower_dunes): "proportional" (both,
# in proportion to their height; every variant not listed), "front" (the seaward row only), "equal" (both alike) or
# "lateral" (the face retreats through the front row, which keeps its height until its sand is gone)
COLLISION_ROWS = {v: "lateral" for v in TOE_AT_BERM if "_lat_" in v}
# cal_automorph_module_cs10 at stronger Cs, island-wide, read on GIS 9-15 where the dunes stay too high (2026-11-09, test)
AUTOMORPH_CS_SWEEP = {f"cal_automorph_module_cs{round(cs * 1e6):d}": cs for cs in (3e-5, 5e-5, 1e-4, 2e-4)}
for _v, _cs in AUTOMORPH_CS_SWEEP.items():
    VARIANTS[_v] = ("storm", _cs, False, (0.02, 0.04))
    TOE_DOMAIN_AVG[_v] = (1996, 1997)            # overridden by the Automorph toes
    LARSON_ALONE = LARSON_ALONE + (_v,)
    RODANTHE_40 = RODANTHE_40 + (_v,)
AUTOMORPH_TOES = AUTOMORPH_TOES + tuple(AUTOMORPH_CS_SWEEP)
# cal_automorph_module_cs10 with automatic Automorph toes on GIS 9-15(HAT_toe_automorph_cells.py, 2026-11-09, test): each
# cell's own toe, or the median of the same cells for the whole domain; the other domains keep the hand-checked toes
AUTOMORPH_CELL_FILE = TRANSECT_TOES_FILE.with_name("automorph_cell_toes_1996.csv")
AUTOMORPH_CELL_TOES = {"cal_automorph_cells_9_15": "toe_m_navd88", "cal_automorph_cellmed_9_15": "toe_domain_median_m_navd88"}
for _v in AUTOMORPH_CELL_TOES:
    VARIANTS[_v] = ("storm", 1e-5, False, (0.02, 0.04))
    TOE_DOMAIN_AVG[_v] = (1996, 1997)            # overridden by the Automorph toes
    LARSON_ALONE = LARSON_ALONE + (_v,)
    RODANTHE_40 = RODANTHE_40 + (_v,)
AUTOMORPH_TOES = AUTOMORPH_TOES + tuple(AUTOMORPH_CELL_TOES)
# The automatic Automorph toe of every cell of every domain, at the Cs of the domain-toe sweep (2026-11-09, test): as
# cal_automorph_module_cs10 and AUTOMORPH_CS_SWEEP, only per cell
AUTOMORPH_CELLS_ALL_FILE = TRANSECT_TOES_FILE.with_name("automorph_cell_toes_all_1996.csv")
AUTOMORPH_CELLS_ALL = {f"cal_automorph_cellsall_cs{_c}{_s}": _c * 1e-6 for _c in (10, 30, 50, 100, 200)
                       for _s in ("", "_noscrape")}       # _noscrape: with the road fix (NO_FIRST_SCRAPE)
NO_FIRST_SCRAPE = NO_FIRST_SCRAPE + tuple(v for v in AUTOMORPH_CELLS_ALL if v.endswith("_noscrape"))
for _v, _cs in AUTOMORPH_CELLS_ALL.items():
    VARIANTS[_v] = ("storm", _cs, False, (0.02, 0.04))
    TOE_DOMAIN_AVG[_v] = (1996, 1997)            # overridden by the cell toes
    LARSON_ALONE = LARSON_ALONE + (_v,)
    RODANTHE_40 = RODANTHE_40 + (_v,)
# Domain toes with the per-cell beach widths (survey medians), for the comparison with the domain beach factor
CELL_WIDTHS = ("cal_step6_toeavg_cw",)
RODANTHE_40 = RODANTHE_40 + ("cal_step3_module_dt", "cal_step4_road_dt")
CAL_STEPS = CAL_STEPS + tuple(TOE_DOMAIN_AVG)
COLLISION_GIS.update({v: range(6, 91) for v in TOE_DOMAIN_AVG})
BEACH_FOLLOWS_SHORELINE = ("cal_g16_nocape_dynw", "cal_g16_nocape_dynw_sw")
WIDTH_FROM_START = TOES_FROM_START + ("cal_g16_nocape_dynw_sw",)
# The road's dune rebuild per variant, (setback gate m, rebuild trigger m MHW or None for the hindcast's berm floor),
# from the lidar rebuilds along NC-12 (dune-storm-response/2026-10-29-collision-test-2009-2025, tables/
# observed_dune_rebuilds.csv): rebuilt up to 140 m behind the dune; crest before a rebuild median 2.66 m NAVD88
ROAD_DUNE_REBUILD = {"cal_g16_nocape_rd150": (150.0, None), "cal_g16_nocape_rd150h": (150.0, 2.66 - MHW_NAVD88_M),
                     "cal_final_rd150u": (150.0, None),
                     "cal_alltoes_cs110_road": (150.0, None), "cal_alltoes_cs110_final": (150.0, None)}
# A 198 m shoreline-to-road gate: the Bonner Bridge Replacement FEIS (NCDOT/FHWA 2008, p. 2-128) lists a shoreline
# within 650 ft (198 m) of the road as a monitoring warning sign for bridge phasing (the landward reach of Isabel's
# washover fans), not as a dune maintenance rule; NCDOT's criterion is 230 ft (70 m) from MHW to the pavement
# (p. 2-64; input_prep/4-mgmt-forcings/road_dune_rebuild/NCDOT_road_dune_criteria.txt). The road manager measures
# from the dune, so each domain's gate is 198 m less its median beach width (the toe table's, as the module's)
SHORELINE_GATE, SHORELINE_GATE_M = ("cal_g16_nocape_rd198", "cal_g16_nocape_rd198_rl02", "cal_g16_nocape_rd198_cs145",
                                    "cal_g16_nocape_rd198_dry", "cal_g16_nocape_rd198_dry_cs087"), 198.0
# The collision beach width on the dry beach, beach start to toe (2026-11-03; HAT_dune_toes.py since then); the
# road gate above keeps the width to the crest, so that only the beach factor changes
DRY_BEACH = ("cal_g16_nocape_rd198_dry", "cal_g16_nocape_rd198_dry_cs087")
_VF_ON = tuple(v for v in VILLAGE_FILTER if v.startswith("cal_g16_nocape_rd198_dry"))
SHORELINE_GATE, DRY_BEACH = SHORELINE_GATE + _VF_ON, DRY_BEACH + _VF_ON
SHORELINE_GATE = SHORELINE_GATE + ("cal_control_g16_rd198", "cal_step4_road", "cal_step5_rodanthe90",
                                   "cal_step6_starttoes", "cal_step6_toeavg", "cal_step4_road_dt",
                                   "cal_step6_toeavg_cw", "cal_step6_toeavg_la", *LA_CS_SWEEP, "cal_step4_road_la",
                                   "cal_automorph_cs10", "cal_automorph_cs19")
DRY_BEACH = DRY_BEACH + CAL_STEPS[2:]
# The beach lowering the water at the dune, m per m of beach (DuneCollisionManager beach_runup_loss), per variant
RUNUP_LOSS = {"cal_g16_nocape_rd198_rl02": 0.02}
# Per-village dune caps, m NAVD88: the lidar 95th percentile of each village's crest cells
# (input_prep/4-mgmt-forcings/village_dune_cap/HAT_village_dune_cap.py), in place of the one 8.0 m
VILLAGE_DOMAINS = {"Buxton": (7, 8), "Avon": (21, 31), "Salvo": (68, 71), "Waves": (72, 76), "Rodanthe": (77, 83)}
VILLAGE_CAPS_M = {"Buxton": 7.6, "Avon": 7.9, "Salvo": 6.7, "Waves": 8.2, "Rodanthe": 6.8}
VILLAGE_CAP_VARIANTS = ("cal_g16_nocape_vcap", "cal_g16_nocape_vcap_ow0")
# No bulldozed overwash onto the Rodanthe village dunes (the manager's overwash_to_dune share 0)
RODANTHE_NO_OW_TO_DUNE = ("cal_g16_nocape_vcap_ow0",)
# Groin strength per sweep variant; every other variant runs at GROIN_STRENGTH_BASE, the runner's default when
# these runs were made (hat_run.yaml moved to 0.16 on 2026-10-29)
GROIN_STRENGTH_BASE = 0.10
GROIN_STRENGTH = {f"storm_w45_r04_{h}_g{round(g * 100):02d}": g for h in ("hs2p5", "hs3") for g in (0.04, 0.06, 0.08, 0.13, 0.16, 0.19)}
GROIN_STRENGTH.update({f"cal_g{round(g * 100):02d}": g for g in (0.13, 0.16, 0.19, 0.22)})
GROIN_STRENGTH["cal_control_g16"] = GROIN_STRENGTH["cal_hs2p75"] = GROIN_STRENGTH["cal_g16_nocape"] = 0.16
GROIN_STRENGTH.update({"cal_g16_capeg05": 0.16, "cal_g16_capeg10": 0.16, "cal_g16_nocape_vcap": 0.16,
                       "cal_g16_nocape_vcap_ow0": 0.16, **{v: 0.16 for v in TOES_FROM_START + BEACH_FOLLOWS_SHORELINE + tuple(ROAD_DUNE_REBUILD) + SHORELINE_GATE}})
GROIN_STRENGTH.update({v: 0.16 for v in VILLAGE_FILTER})
GROIN_STRENGTH.update({v: 0.16 for v in CAL_STEPS})
GROIN_STRENGTH["cal_step2_growth_noscrape"] = 0.16
# CASCADE's wave height (m) per variant off the runner's 2.0
HS_M = {"control_hs3": 3.0, "storm_w45_r04_hs3": 3.0, "control_hs3_ends": 3.0, "storm_w45_r04_hs3_ends": 3.0,
        "control_hs2p5": 2.5, "storm_w45_r04_hs2p5": 2.5, "storm_cs145_r04_hs3_ends": 3.0, "storm_cs145_r04_hs2p5": 2.5,
        **{v: (2.5 if "_hs2p5_" in v else 3.0) for v in GROIN_STRENGTH},
        **{v: 3.0 for v in VARIANTS if v.startswith("cal_")},
        **{f"cal_hs{h:g}".replace(".", "p"): h for h in (2.0, 2.5, 2.75, 3.5, 4.0)}}
# The end rates (m/yr) solved at Hs 3.0 with StormStart 1 (input_prep/7-source-sink/2-calibrate/
# be_end_solve_net_change.py --hs 3, 2026-10-28, step 8), in place of edgeBE's Hs 2.0 values
ENDS_HS3 = {1: 3.0688, 90: 2.0647}
# The same solve at each wave height on the 2026-10-29 defaults (village dune cap 8.0 m; --hs 2 ... 4, last step)
ENDS_BY_HS = {2.0: {1: 3.9595, 90: 1.2255}, 2.5: {1: 3.4009, 90: 1.5958}, 3.0: ENDS_HS3,
              3.5: {1: 2.8201, 90: 2.6656}, 4.0: {1: 2.6813, 90: 3.5346}}
# Hs 2.75 not solved: halfway between the 2.5 and 3.0 ends, for a quick look (cal_hs2p75)
ENDS_BY_HS[2.75] = {g: round((ENDS_BY_HS[2.5][g] + ENDS_BY_HS[3.0][g]) / 2, 4) for g in (1, 90)}
ENDS = {"control_hs3_ends": ENDS_HS3, "storm_w45_r04_hs3_ends": ENDS_HS3, "storm_cs145_r04_hs3_ends": ENDS_HS3,
        **{v: ENDS_HS3 for v in GROIN_STRENGTH if "_hs3_" in v},
        **{v: ENDS_BY_HS[HS_M[v]] for v in VARIANTS if v.startswith("cal_")}}
# The beach width W0 (m) beyond which a dune loses nothing to collision, per variant using the beach
BEACH_W0_M = {"storm_w45": 45.0, "storm_w45_r04": 45.0, "storm_w45_r06": 45.0, "storm_w45_r04_ss1": 45.0,
              "storm_w45_r04_hs3": 45.0, "storm_w45_r04_hs3_ends": 45.0, "storm_w45_r04_hs2p5": 45.0,
              "storm_cs145_r04_hs3_ends": 45.0, "storm_cs145_r04_hs2p5": 45.0,
              **{v: 45.0 for v in GROIN_STRENGTH if VARIANTS[v][0]},
              **{v: 45.0 for v in VARIANTS if v.startswith(("cal_r", "cal_zone", "cal_hs", "cal_g"))}}
BEACH_W0_M.update({v: 25.0 for v in DRY_BEACH})
# The sink, m/yr as Barrier3D's Rat ((-) erosion): CoastSat LRR less the control's at GIS 9-14, where the
# control erodes too little and edgeBE has no rate
BUXTON_SINK = {9: -0.61, 10: -0.56, 11: -0.35, 12: -1.77, 13: -1.92, 14: -0.89}
SINK_VARIANTS = ("sink", "storm_cs2e5_r06_sink")
# The stable dunes the growth rate was calibrated on (HAT_stable_dune_growth_rate.py)
STABLE_GIS = (20, 21, 22, 59, 60, 61, 62)
WORKERS = 4
REPORT_GIS = range(5, 21)
# -----------------------------------------------------------------------------


# Toes

# Beach start, crest and toe per profile (m MHW, ocean first): the input step's picker
# (input_prep/1-barrier3d-domains/3-dune-toes/HAT_dune_toes.py), so the calibration and the runner read one rule.
# Until 2026-11-03 this was a copy of it, with the chord to the highest crest (deepest dip); since, the chord runs
# to the first dune (FOREDUNE_DROP_M there)
def crest_and_toe(z, i0, i1, X):
    p = str(_HAT_DIR / "input_prep" / "1-barrier3d-domains" / "3-dune-toes")
    if p not in sys.path:
        sys.path.insert(0, p)
    import HAT_dune_toes
    return HAT_dune_toes.crest_and_toe(z, i0, i1, X)


def toes():
    sys.path.insert(0, str(_HAT_DIR / "input_prep" / "1-barrier3d-domains" / "2-survey-crests"))
    import HAT_survey_dune_crests as S
    import json
    import rasterio
    X, tv = S.X, S.tv
    prod, surveys = S.STARTS[START]
    windows = json.loads(tv.picks_dir(START).joinpath(
        f"HAT_dune_search_windows_{tv.topo_dirs(START)[2]}.json").read_text(encoding="utf-8"))
    frames = {g: S.frame(START, g, windows) for g in range(1, 91)}
    clips = S.product(prod).gapfill_1m
    grids = {g: S.clip_grid(clips / f"clip_domain_{g}_filled.tif") for g in range(1, 91)}
    offsets = S.SURVEY_VERTICAL_OFFSET_M.get(START, {})
    rows = []

    def add(survey, g, z, i0, i1):
        r = crest_and_toe(z, i0, i1, X)
        for cell in range(50):
            rows.append(dict(survey=survey, gis=g, cell=cell,
                             crest_m_navd88=r["crest"][cell] + MHW_NAVD88_M,
                             toe_m_navd88=r["toe"][cell] + MHW_NAVD88_M,
                             beach_width_m=(r["crest_idx"][cell] - r["beach_start_idx"][cell]) * 10,
                             toe_to_crest_m=(r["crest_idx"][cell] - r["toe_idx"][cell]) * 10))

    # The start's own arrays: the dunes the model starts from
    for g in range(1, 91):
        shear, c0, i0, i1 = frames[g]
        z = S.profiles(np.load(tv.npy_dirs(START)[0] / f"domain_{g}.npy").astype(float), shear, c0)
        add("start", g, z, i0, i1)
    print("start arrays done", flush=True)
    for year, (kind, path, flown, state, note) in surveys.items():
        files = sorted(path.parent.glob(path.name)) if "*" in path.name else [path]
        srcs = [rasterio.open(f) for f in files]
        try:
            unit = "m" if kind == "mosaic" else S.unit_of(srcs, START, grids)[0]
            scale = S.US_FT if unit == "ft" else 1.0
            for g in range(1, 91):
                shear, c0, i0, i1 = frames[g]
                on_clip = (S.mosaic_on_clip(srcs[0], grids[g]) if kind == "mosaic"
                           else S.raster_on_clip(srcs, grids[g]))
                z = S.profiles(S.downsample(on_clip * scale + offsets.get(year, 0.0)), shear, c0)
                # The survey's windows follow the dune, as the crest script moved them (2026-11-03)
                add(year, g, z, *S.survey_windows(START, year, g, i0, i1))
        finally:
            for s in srcs:
                s.close()
        print(f"{year} done", flush=True)
    t = pd.DataFrame(rows)
    # The crest script's spike screen: keep its status beside each survey cell
    sc = pd.read_csv(SURVEY_CRESTS)[["survey", "gis", "cell", "status", "crest_m_navd88"]]
    sc["survey"] = sc.survey.astype(str)
    t["survey"] = t.survey.astype(str)
    t = t.merge(sc.rename(columns={"crest_m_navd88": "crest_script_m_navd88"}), on=["survey", "gis", "cell"],
                how="left")
    t.loc[t.survey == "start", "status"] = "ok"
    EXP_DIR.mkdir(parents=True, exist_ok=True)
    t.round(3).to_csv(EXP_DIR / "dune_toes.csv", index=False)
    same = (t.crest_m_navd88 - t.crest_script_m_navd88).abs()
    print(f"crest check against survey_dune_crests.csv: max |difference| {same.max():.4f} m")
    ok = t[(t.status == "ok") & t.toe_m_navd88.notna()]
    print(ok.groupby("survey")[["toe_m_navd88", "crest_m_navd88", "beach_width_m"]].median().round(2))
    print(f"wrote {EXP_DIR / 'dune_toes.csv'}")


# Hourly total water level

def twl():
    from site_layer import hat_env_forcings as env
    w = pd.read_csv(env.DUCK_GAUGE_FILE, index_col="t", parse_dates=True)["v"]
    wis = pd.read_csv(env.WIS_FILE, index_col="time", parse_dates=True)[["waveHs", "waveTp"]]
    idx = pd.date_range(TWL_START, TWL_END, freq="h")
    df = pd.DataFrame({"water_level": w.reindex(idx), "Hs": wis.waveHs.reindex(idx),
                       "Tp": wis.waveTp.reindex(idx)}).dropna()
    # Stockdon et al. (2006) R2%, as historical_storm_creation_v3_HAT.calculate_r2_percent
    L0 = 9.81 * df.Tp ** 2 / (2 * np.pi)
    setup = 0.35 * BEACH_SLOPE * np.sqrt(df.Hs * L0)
    swash = np.sqrt((0.75 * BEACH_SLOPE * np.sqrt(df.Hs * L0)) ** 2 + (0.06 * np.sqrt(df.Hs * L0)) ** 2)
    df["R2"] = 1.1 * (setup + swash / 2)
    df["TWL"] = df.water_level + df.R2
    # Check: each storm's peak in the series is the peak of these hours
    ev = pd.read_csv(STORM_SUMMARY, parse_dates=["StartTime", "EndTime"])
    peak = np.array([df.TWL.loc[a:b].max() for a, b in zip(ev.StartTime, ev.EndTime)])
    diff = peak - (ev.Rhigh * 10 + MHW_NAVD88_M)
    print(f"{len(df)} hours; storm peaks against the series: max |difference| {np.nanmax(np.abs(diff)):.4f} m "
          f"over {len(ev)} events")
    EXP_DIR.mkdir(parents=True, exist_ok=True)
    df.round(4).to_csv(EXP_DIR / "twl_hourly_1996_2009.csv", index_label="time")
    print(f"wrote {EXP_DIR / 'twl_hourly_1996_2009.csv'}")


# Calibration

# The survey rows that passed the crest script's screen
def survey_table():
    t = pd.read_csv(EXP_DIR / "dune_toes.csv")
    t = t[(t.survey != "start") & (t.status == "ok")].copy()
    t["survey"] = t.survey.astype(int)
    return t


# Per cell, the median toe over the surveys whose toe passes the screen (m NAVD88)
def cell_toes(t):
    good = t[~t.survey.isin(TOE_SURVEYS_EXCLUDED) & (t.toe_m_navd88 > TOE_MIN_M_NAVD88)
             & (t.toe_m_navd88 < t.crest_m_navd88 - TOE_BELOW_CREST_M)]
    return good.groupby(["gis", "cell"]).toe_m_navd88.median().rename("toe_cell_m_navd88")


# The model's storm series: peak Rhigh (m NAVD88), period (s) and duration (h) per storm
def storm_table():
    st = pd.read_csv(STORM_SUMMARY, parse_dates=["StartTime", "EndTime"])
    return st.assign(rhigh_m_navd88=st.Rhigh * 10 + MHW_NAVD88_M)


# Per cell, the Larson exposure sum((Rhigh - toe)^2 duration / T) over the given storms (m^2)
def storm_exposure(storms, toe):
    excess = np.clip(storms.rhigh_m_navd88.to_numpy()[None, :] - toe[:, None], 0, None)
    return (excess ** 2 * (storms.duration * 3600.0 / storms.period).to_numpy()[None, :]).sum(axis=1)


# First and last day flown, from "YYYY-MM" or "YYYY-MM to YYYY-MM"
def flown_span(flown):
    parts = flown.split(" to ")
    first = pd.Timestamp(parts[0] + "-01")
    last = pd.Timestamp(parts[-1] + "-01") + pd.offsets.MonthEnd(0) + pd.Timedelta(hours=23)
    return first, last


def calibrate():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from site_layer.hat_figure_style import apply_style, record_caption
    apply_style()

    t = survey_table()
    toe_cell = cell_toes(t)
    print(f"per-cell toe: {len(toe_cell)} cells, median {toe_cell.median():.2f}, "
          f"10-90% {toe_cell.quantile(0.1):.2f}-{toe_cell.quantile(0.9):.2f} m NAVD88")
    flown = pd.read_csv(SURVEY_CRESTS).groupby("survey").flown.first()
    h = pd.read_csv(EXP_DIR / "twl_hourly_1996_2009.csv", index_col="time", parse_dates=True)
    storms = storm_table()
    years = sorted(t.survey.unique())
    cells, summary = [], []
    for a, b in zip(years[:-1], years[1:]):
        t0, t1 = flown_span(flown[a])[1], flown_span(flown[b])[0]
        hh = h.loc[t0:t1]
        pa = t[t.survey == a].set_index(["gis", "cell"])
        pb = t[t.survey == b].set_index(["gis", "cell"])
        j = pa.join(pb[["crest_m_navd88"]], rsuffix="_b", how="inner").join(toe_cell, how="inner").dropna(
            subset=["toe_cell_m_navd88", "crest_m_navd88", "crest_m_navd88_b"])
        twl_max = hh.TWL.max()
        # Not overtopped: the water never reached the earlier crest
        j = j[j.crest_m_navd88 > twl_max].copy()
        # Exposure, m^2: sum over hours of (TWL - toe)^2 * 3600 s / Tp
        toe = j.toe_cell_m_navd88.to_numpy()[:, None]
        excess = np.clip(hh.TWL.to_numpy()[None, :] - toe, 0, None)
        j["exposure_m2"] = (excess ** 2 * 3600.0 / hh.Tp.to_numpy()[None, :]).sum(axis=1)
        j["hours_above_toe"] = (excess > 0).sum(axis=1)
        # The same per storm, as the model sees it: each storm's peak Rhigh held over its duration
        ss = storms[(storms.StartTime >= t0) & (storms.StartTime <= t1)]
        j["exposure_storm_m2"] = storm_exposure(ss, j.toe_cell_m_navd88.to_numpy())
        j["crest_loss_m"] = j.crest_m_navd88 - j.crest_m_navd88_b
        j["interval"] = f"{a}-{b}"
        j["years"] = (t1 - t0).days / 365.25
        cells.append(j.reset_index())
        summary.append(dict(interval=f"{a}-{b}", years=round(j.years.iloc[0], 2), twl_max_m_navd88=round(twl_max, 2),
                            cells=len(j), toe_median=j.toe_cell_m_navd88.median(),
                            hours_above_toe_median=j.hours_above_toe.median(),
                            exposure_median_m2=j.exposure_m2.median(),
                            crest_loss_median_m=j.crest_loss_m.median()))
    c = pd.concat(cells, ignore_index=True)
    s = pd.DataFrame(summary).set_index("interval")

    # Fit dH = (4 Cs / W) E - g * years: Cs the collision coefficient, g the net growth between surveys
    A = np.column_stack([4 * c.exposure_m2 / DUNE_WIDTH_M, -c.years])
    cs_all, g_all = np.linalg.lstsq(A, c.crest_loss_m, rcond=None)[0]
    # The same on domain medians per interval, which damps the lidar's cell-to-cell noise
    d = c.groupby(["interval", "gis"]).agg(exposure_m2=("exposure_m2", "median"), years=("years", "first"),
                                           crest_loss_m=("crest_loss_m", "median")).reset_index()
    A = np.column_stack([4 * d.exposure_m2 / DUNE_WIDTH_M, -d.years])
    cs_every, _ = np.linalg.lstsq(A, d.crest_loss_m, rcond=None)[0]
    fit = d.interval.isin(FIT_INTERVALS).to_numpy()
    cs_dom, g_dom = np.linalg.lstsq(A[fit], d.crest_loss_m[fit], rcond=None)[0]
    # Per interval, the growth term held at the island fit
    s["cs_interval"] = [
        float(np.sum(4 * x.exposure_m2 / DUNE_WIDTH_M * (x.crest_loss_m + g_dom * x.years))
              / np.sum((4 * x.exposure_m2 / DUNE_WIDTH_M) ** 2)) if (x.exposure_m2 > 0).any() else np.nan
        for _, x in d.groupby("interval", sort=False)]
    c["predicted_loss_m"] = 4 * cs_dom * c.exposure_m2 / DUNE_WIDTH_M - g_dom * c.years
    d["predicted_loss_m"] = 4 * cs_dom * d.exposure_m2 / DUNE_WIDTH_M - g_dom * d.years
    r_cells = np.corrcoef(c.predicted_loss_m, c.crest_loss_m)[0, 1]
    r_dom = np.corrcoef(d.predicted_loss_m, d.crest_loss_m)[0, 1]
    # Per-storm form, fit the same way: the Cs the model's module takes
    ds = c.groupby(["interval", "gis"]).agg(exposure_m2=("exposure_storm_m2", "median"), years=("years", "first"),
                                            crest_loss_m=("crest_loss_m", "median")).reset_index()
    A = np.column_stack([4 * ds.exposure_m2 / DUNE_WIDTH_M, -ds.years])
    cs_storm, g_storm = np.linalg.lstsq(A[fit], ds.crest_loss_m[fit], rcond=None)[0]
    r_storm = np.corrcoef(A[fit] @ [cs_storm, g_storm], ds.crest_loss_m[fit])[0, 1]
    s["cs_storm_interval"] = [
        float(np.sum(4 * x.exposure_m2 / DUNE_WIDTH_M * (x.crest_loss_m + g_storm * x.years))
              / np.sum((4 * x.exposure_m2 / DUNE_WIDTH_M) ** 2)) if (x.exposure_m2 > 0).any() else np.nan
        for _, x in ds.groupby("interval", sort=False)]
    s["exposure_storm_median_m2"] = c.groupby("interval", sort=False).exposure_storm_m2.median()

    pd.set_option("display.width", 250)
    pd.set_option("display.max_columns", 30)
    print(s.round(3))
    print(f"\ncells:   Cs {cs_all:.2e}  growth {g_all:+.3f} m/yr   (n={len(c)})")
    print(f"domains: Cs {cs_dom:.2e}  growth {g_dom:+.3f} m/yr   r {r_dom:.2f} (cells r {r_cells:.2f}, n={len(d)}); "
          f"every interval, 2001-2004 too: Cs {cs_every:.2e}")
    print(f"per storm: Cs {cs_storm:.2e}  growth {g_storm:+.3f} m/yr   r {r_storm:.2f} (domain medians)")
    print(f"Larson et al. (2004) laboratory Cs: 1.7e-04")
    c.round(4).to_csv(EXP_DIR / "collision_calibration_cells.csv", index=False)
    s.assign(cs_domains=cs_dom, growth_domains_m_yr=g_dom, r_domains=r_dom, cs_storm_domains=cs_storm,
             growth_storm_m_yr=g_storm, r_storm_domains=r_storm).round(6).to_csv(
        EXP_DIR / "collision_calibration.csv")

    fig, axs = plt.subplots(1, 2, figsize=(10, 4.2))
    ax = axs[0]
    for name, x in d.groupby("interval", sort=False):
        ax.plot(4 * x.exposure_m2 / DUNE_WIDTH_M, x.crest_loss_m, "o", ms=3.5, alpha=0.8, label=name)
    xs = np.linspace(0, (4 * d.exposure_m2 / DUNE_WIDTH_M).max(), 50)
    ax.plot(xs, cs_dom * xs, color="k", lw=1.2, label=f"fit, Cs = {cs_dom:.1e}")
    ax.axhline(0, color="0.6", lw=0.8)
    ax.set_xlabel("4 Σ (TWL − toe)² Δt / T / 20 m   (m)")
    ax.set_ylabel("Crest lowering between surveys (m)")
    ax.set_title("Domain medians, cells never overtopped", loc="left")
    ax.legend(fontsize=7, ncol=2)
    ax = axs[1]
    toe0 = t[t.survey == years[0]].groupby("gis").toe_m_navd88.median()
    crest0 = t[t.survey == years[0]].groupby("gis").crest_m_navd88.median()
    ax.plot(toe0.index, crest0, color="0.5", lw=1.2, label=f"crest, {years[0]} survey")
    ax.plot(toe0.index, toe0, color="#2166AC", lw=1.6, label=f"toe, {years[0]} survey")
    st = pd.read_csv(EXP_DIR / "dune_toes.csv")
    st = st[st.survey == "start"].groupby("gis").toe_m_navd88.median()
    ax.plot(st.index, st, color="#2166AC", lw=1.0, ls="--", label="toe, model start arrays")
    ax.axhline(1.7, color="0.3", lw=0.8, ls=":", label="Barrier3D berm, 1.7 m")
    ax.set_xlabel("GIS domain")
    ax.set_ylabel("m NAVD88")
    ax.set_title("Dune toe and crest by domain (median)", loc="left")
    ax.legend(fontsize=7)
    fig.tight_layout()
    out = EXP_DIR / "collision_calibration.png"
    fig.savefig(out, dpi=200, facecolor="white")
    record_caption(out, (
        "Left: per domain and survey interval, the median crest lowering of the 10 m cells the total water "
        "level never overtopped between the two surveys, against the Larson et al. (2004) wave-impact exposure "
        "4 Σ (TWL − z_toe)² Δt / T over the hours above the earlier survey's toe, per 20 m of dune width; line: "
        f"the least-squares Cs with a net growth term ({g_dom:+.3f} m/yr) held for every interval. "
        "Right: domain-median toe and crest, the toe as the point farthest below the line from the beach start "
        f"(0.5 m MHW) to the crest. Experiment {TAG}."))
    print(f"wrote {out}")


# Empirical scaling, in the terms of Long et al. (2014): dz / D_high against R_high / D_high, and the beach width

def empirical():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from scipy.optimize import curve_fit
    from site_layer.hat_figure_style import apply_style, record_caption
    apply_style()

    c = pd.read_csv(EXP_DIR / "collision_calibration_cells.csv")
    r_high = pd.read_csv(EXP_DIR / "collision_calibration.csv").set_index("interval").twl_max_m_navd88
    c = c[c.interval != "2001-2004"].copy()
    c["r_high"] = c.interval.map(r_high)
    c["x"] = c.r_high / c.crest_m_navd88
    c["y"] = c.crest_loss_m / c.crest_m_navd88
    # Less each interval's background (growth, repair, lidar offset): the median on the wide beaches
    c["y_net"] = c.y - c.interval.map(c[c.beach_width_m > BACKGROUND_BEACH_WIDTH_M].groupby("interval").y.median())
    st = c[c.interval.isin(EMPIRICAL_INTERVALS) & (c.y.abs() < 1.5)]

    def rule(w, a, w0):
        return a * np.clip(1 - w / w0, 0, None)

    (a, w0), _ = curve_fit(rule, st.beach_width_m, st.y_net, p0=[0.4, 40.0])
    r_cells = np.corrcoef(rule(st.beach_width_m, a, w0), st.y_net)[0, 1]
    dm = st.assign(p=rule(st.beach_width_m, a, w0)).groupby(["interval", "gis"])[["p", "y_net", "x"]].median()
    r_dom, r_dom_x = np.corrcoef(dm.p, dm.y_net)[0, 1], np.corrcoef(dm.x, dm.y_net)[0, 1]
    print(f"dz/D_high = {a:.3f} max(0, 1 - W / {w0:.1f} m)   r {r_cells:.2f} cells, {r_dom:.2f} domain medians "
          f"(R_high/D_high alone: {r_dom_x:.2f}); fit on {', '.join(EMPIRICAL_INTERVALS)}, n={len(st)}")

    widths = [0, 20, 30, 40, 60, 400]
    tab = c.groupby(["interval", pd.cut(c.beach_width_m, widths), pd.cut(c.x, [0, .6, .7, .8, .9, 1])],
                    observed=True).y.agg(["size", "median"]).reset_index()
    tab.columns = ["interval", "beach_width_m", "r_high_over_d_high", "cells", "dz_over_d_high_median"]
    tab.round(4).to_csv(EXP_DIR / "collision_empirical.csv", index=False)
    pd.DataFrame([dict(a=a, w0_m=w0, r_cells=r_cells, r_domains=r_dom, r_domains_rhigh_only=r_dom_x,
                       intervals=" ".join(EMPIRICAL_INTERVALS), cells=len(st))]).round(4).to_csv(
        EXP_DIR / "collision_empirical_fit.csv", index=False)

    fig, axs = plt.subplots(1, 2, figsize=(10, 4.2))
    shades = ["#7F2704", "#D94801", "#FD8D3C", "#FDBE85", "#BDBDBD"]
    xb = np.arange(0.45, 1.0, 0.1)
    ax = axs[0]
    for (lo, hi), col in zip(zip(widths[:-1], widths[1:]), shades):
        s = st[(st.beach_width_m > lo) & (st.beach_width_m <= hi)]
        m = s.groupby(pd.cut(s.x, np.arange(0.4, 1.01, 0.1)), observed=False).y
        ax.plot(xb, m.median(), "o-", color=col, lw=1.6, ms=5,
                label=f"beach {lo}–{hi} m" if hi < 400 else f"beach > {lo} m")
    ax.axhspan(-0.1, 0.1, color="0.9", zorder=0, label="Long et al. (2014): within 10%")
    ax.axhline(0, color="0.6", lw=0.8)
    ax.set_xlabel("R_high / D_high")
    ax.set_ylabel("Crest lowering / D_high")
    ax.set_title("Dennis and 2004–05, cells in collision", loc="left")
    ax.legend(fontsize=7)
    ax = axs[1]
    wb = np.arange(5, 100, 10)
    for iv, x in c.groupby("interval"):
        storm = iv in EMPIRICAL_INTERVALS
        m = x.groupby(pd.cut(x.beach_width_m, np.arange(0, 101, 10)), observed=False).y_net.median()
        ax.plot(wb, m, "o-" if storm else "-", color="#D94801" if storm else "0.7", lw=1.6 if storm else 1.0,
                ms=4, label=f"{iv} (R_high {r_high[iv]:.2f} m)", zorder=3 if storm else 2)
    ws = np.linspace(0, 100, 200)
    ax.plot(ws, rule(ws, a, w0), color="k", lw=1.4, label=f"{a:.2f} max(0, 1 − W / {w0:.0f} m)")
    ax.axhline(0, color="0.6", lw=0.8)
    ax.set_xlabel("Beach width, 0.5 m MHW contour to the crest (m)")
    ax.set_ylabel("Crest lowering / D_high, less background")
    ax.set_title("Every survey interval", loc="left")
    ax.legend(fontsize=6.5)
    fig.tight_layout()
    out = EXP_DIR / "collision_empirical.png"
    fig.savefig(out, dpi=200, facecolor="white")
    record_caption(out, (
        "Dune crest change across the 1996–2009 lidar surveys in the normalisation of Long et al. (2014), for the "
        "10 m cells the total water level (Duck gauge + Stockdon R2%, slope 0.06) never overtopped between two "
        "surveys. Left: median crest lowering over pre-storm crest elevation (NAVD88) against R_high / D_high, "
        "per beach-width class, for the two storm intervals (1998–99, Dennis; 2004–05); shaded, the collision-"
        "regime bound Long et al. report. Right: the same, less each interval's median on beaches wider than "
        f"{BACKGROUND_BEACH_WIDTH_M:g} m, against beach width (beach start at 0.5 m MHW to the crest); storm "
        f"intervals orange, the rest grey; line, the least-squares fit on the storm intervals (r {r_cells:.2f} "
        f"cells, {r_dom:.2f} domain medians). Experiment {TAG}."))
    print(f"wrote {out}")


# Model runs

# Crest lowering (m) by Larson et al. (2004) over the hours in the collision regime: above the toe, below the crest
def larson_loss(hours, cs):
    twl, tp = hours.TWL.to_numpy(), hours.Tp.to_numpy()

    def loss(g, crest, z_toe):
        reach = (twl[None, :] > z_toe[:, None]) & (twl[None, :] < crest[:, None])
        excess = np.where(reach, twl[None, :] - z_toe[:, None], 0.0)
        return cs * 4 * (excess ** 2 * 3600.0 / tp[None, :]).sum(axis=1) / DUNE_WIDTH_M
    return loss


# Crest lowering (m) by the beach-width rule, once, where the year's highest water is between the toe and the crest
def width_loss(hours, width, a, w0, min_twl):
    r_high = hours.TWL.max()

    def loss(g, crest, z_toe):
        if min_twl is not None and r_high < min_twl:
            return np.zeros_like(crest)
        w = width[g][: len(crest)]
        hit = (r_high > z_toe) & (r_high < crest)
        return np.where(hit, a * np.clip(1 - w / w0, 0, None) * crest, 0.0)
    return loss


# Collision erosion over one calendar year for every real domain, applied to the state just written, the
# lowering per cell from loss(gis, crest, toe); returns per GIS the mean crest lowering (m), the sand removed
# (m^3/m) and the shoreline move (m, + seaward)
def collide(cascade, dom, toe, loss, to_beach):
    out = {}
    for g in range(1, 91):
        b = cascade.barrier3d[dom.gis_to_pad(g)]
        t = b._time_index - 1
        rows = b._DuneDomain[t]                                       # dam above the berm, (cells, 2)
        top = rows.max(axis=1)
        crest = (top + b._BermEl) * 10 + MHW_NAVD88_M                 # m NAVD88
        z_toe = toe[g][: len(crest)]
        lowering = np.nan_to_num(loss(g, crest, z_toe))               # m
        # No lower than the toe; a crest already below its toe stays where it is
        new_crest = np.maximum(crest - lowering, np.minimum(z_toe, crest))
        new_top = np.maximum((new_crest - MHW_NAVD88_M) / 10 - b._BermEl, b._DuneRestart)
        safe_top = np.where(top > 0, top, 1.0)
        new_rows = np.where(top[:, None] > 0,
                            np.maximum(rows * (new_top / safe_top)[:, None], b._DuneRestart), rows)
        new_rows = np.minimum(new_rows, rows)                         # erosion only
        removed = rows - new_rows                                     # dam, per cell and row
        b._DuneDomain[t] = new_rows
        # m^3 per m of shore: each row is 10 m cross-shore, averaged over the domain's cells
        vol = float(removed.sum(axis=1).mean()) * 10 * 10
        dx = 0.0
        if to_beach and vol > 0:
            # CASCADE's shoreface_nourishment, in metres: dx_s = 2 V / (2 h_b + D_sf)
            dx = 2 * vol / (2 * b.h_b_TS[-1] * 10 + b._DShoreface * 10)
            b.x_s = b.x_s - dx / 10
            b.x_s_TS[-1] = b.x_s
        out[g] = (float(np.mean((top - new_rows.max(axis=1)) * 10)), vol, dx)
    return out


# Child process: collision erosion after each model year, then the unchanged hindcast
def _launch(variant):
    import runpy
    import cascade_pipeline.hindcast as H
    from site_layer.hatteras_site_config import HATTERAS_DOMAINS
    rule, cs, to_beach, _ = VARIANTS[variant]
    st = survey_table()

    # Per GIS, 50 cells from a (gis, cell) series, the domain median where a cell has none
    def per_gis(sr):
        dom_median = sr.groupby(level="gis").median()
        out = {}
        for g in range(1, 91):
            v = np.full(50, dom_median.get(g, np.nan))
            if g in sr.index.get_level_values("gis"):
                s = sr.loc[g]
                v[s.index.to_numpy()] = s.to_numpy()
            out[g] = v
        return out

    toe = per_gis(cell_toes(st))
    if variant in TOE_DOMAIN_AVG:
        s_ = st[st.survey.isin(TOE_DOMAIN_AVG[variant])]
        s_ = s_[(s_.toe_m_navd88 > TOE_MIN_M_NAVD88) & (s_.toe_m_navd88 < s_.crest_m_navd88 - TOE_BELOW_CREST_M)]
        dm = s_.groupby("gis").toe_m_navd88.mean()
        toe = {g: np.full(50, dm.get(g, dm.median())) for g in range(1, 91)}
        print()
        print(f"EXPERIMENT one toe per domain, the mean over the {TOE_DOMAIN_AVG[variant]} surveys: island "
              f"{dm.mean():.2f} m NAVD88, range {dm.min():.2f}-{dm.max():.2f}", flush=True)
    if variant in TRANSECT_TOES:
        tt = pd.read_csv(TRANSECT_TOES_FILE).set_index("gis").toe_m_navd88
        toe = {g: np.full(50, tt[g]) for g in range(1, 91)}
        print(f"EXPERIMENT one toe per domain from the 1996 transects ({TRANSECT_TOES_FILE.name}): island "
              f"{tt.mean():.2f} m NAVD88, range {tt.min():.2f}-{tt.max():.2f}", flush=True)
    if variant in AUTOMORPH_TOES:
        tt = pd.read_csv(AUTOMORPH_TOES_FILE).set_index("gis").toe_m_navd88
        toe = {g: np.full(50, tt[g]) for g in range(1, 91)}
        print(f"EXPERIMENT one toe per domain, Automorph hand-checked ({AUTOMORPH_TOES_FILE.name}): island "
              f"{tt.mean():.2f} m NAVD88, range {tt.min():.2f}-{tt.max():.2f}", flush=True)
    if variant in AUTOMORPH_CELL_TOES:
        ct = pd.read_csv(AUTOMORPH_CELL_FILE).sort_values(["gis", "cell"])
        for g, c_ in ct.groupby("gis"):
            toe[g] = c_[AUTOMORPH_CELL_TOES[variant]].to_numpy(float)
        print(f"EXPERIMENT automatic Automorph toes on GIS {ct.gis.min()}-{ct.gis.max()} "
              f"({AUTOMORPH_CELL_TOES[variant]}), range {ct[AUTOMORPH_CELL_TOES[variant]].min():.2f}-"
              f"{ct[AUTOMORPH_CELL_TOES[variant]].max():.2f} m NAVD88", flush=True)
    if variant in TOE_AT_BERM:
        toe = {g: np.full(50, TOE_AT_BERM_M) for g in range(1, 91)}
        print(f"EXPERIMENT the toe at the berm, {TOE_AT_BERM_M:g} m NAVD88, every cell", flush=True)
    if variant in AUTOMORPH_CELLS_ALL:
        ct = pd.read_csv(AUTOMORPH_CELLS_ALL_FILE).sort_values(["gis", "cell"])
        toe = {g: c_.toe_m_navd88.to_numpy(float) for g, c_ in ct.groupby("gis")}
        print(f"EXPERIMENT automatic Automorph toe per cell on every domain ({len(ct)} cells): island "
              f"{ct.toe_m_navd88.mean():.2f} m NAVD88, range {ct.toe_m_navd88.min():.2f}-{ct.toe_m_navd88.max():.2f}",
              flush=True)
    if variant in WIDTH_FROM_START:
        sp = pd.read_csv(EXP_DIR / "dune_toes.csv")
        sp = sp[sp.survey == "start"]
        start_width = per_gis(sp.set_index(["gis", "cell"]).beach_width_m)
    if variant in TOES_FROM_START:
        good = sp[(sp.toe_m_navd88 > TOE_MIN_M_NAVD88) & (sp.toe_m_navd88 < sp.crest_m_navd88 - TOE_BELOW_CREST_M)]
        toe = per_gis(good.set_index(["gis", "cell"]).toe_m_navd88)
    hours = pd.read_csv(EXP_DIR / "twl_hourly_1996_2009.csv", index_col="time", parse_dates=True)
    if rule == "width":
        fit = pd.read_csv(EXP_DIR / "collision_empirical_fit.csv").iloc[0]
        flown = pd.read_csv(SURVEY_CRESTS).groupby("survey").flown.first()
        widths = {y: per_gis(st[st.survey == y].set_index(["gis", "cell"]).beach_width_m)
                  for y in sorted(st.survey.unique())}

        # The beach widths of the last survey flown before the year's highest water (the first survey if none)
        def width_for(yh):
            peak = yh.TWL.idxmax()
            before = [y for y in widths if flown_span(flown[y])[1] < peak]
            return widths[before[-1] if before else min(widths)]

        def loss_for(yh):
            return width_loss(yh, width_for(yh), fit.a, fit.w0_m, cs)
    else:
        def loss_for(yh):
            return larson_loss(yh, cs)
    log = []
    orig_build = H.build_cascade

    # The village overwash filter (VILLAGE_FILTER), on every managed village domain's beach/dune manager
    def village_filter(c):
        p = VILLAGE_FILTER.get(variant)
        if p is None:
            return
        n = 0
        for i, mgr in enumerate(c.nourishments):
            if mgr is not None and c._beach_nourishment_module[i]:
                mgr._overwash_filter = p
                n += 1
        print(f"\nEXPERIMENT village overwash filter {p:g}% to the shoreface on {n} village domains", flush=True)

    # the 198 m (650 ft) shoreline-to-road gate (SHORELINE_GATE) as each road manager's gate behind the dune; was applied
    # only with the module on until 2026-11-04 (cal_control_g16_rd198 ran without it)
    def road_gate(c, dom, width):
        if variant not in SHORELINE_GATE:
            return
        w = pd.Series({g: np.nanmedian(width[g]) for g in range(1, 91)})
        w = w.fillna(w.median())
        for g in range(1, 91):
            r = c.roadways[dom.gis_to_pad(g)]
            if r is not None:
                r._dune_mgmt_road_setback_trigger = max(SHORELINE_GATE_M - w[g], 0.0)
        print(f"\nEXPERIMENT road dune rebuild where NC 12 is within {SHORELINE_GATE_M:g} m of the shoreline: "
              f"gates median {SHORELINE_GATE_M - w.median():.0f} m, range {SHORELINE_GATE_M - w.max():.0f}-"
              f"{SHORELINE_GATE_M - w.min():.0f} m", flush=True)

    # Rodanthe (GIS 77-83) held at 40% (RODANTHE_40), and the road's first bulldozing kept off the dunes (NO_FIRST_SCRAPE)
    def rodanthe_40(c, dom):
        if variant in NO_FIRST_SCRAPE:
            n = sum(1 for r in c.roadways if r is not None)
            for r in c.roadways:
                if r is not None:
                    r._discard_first_bulldoze = True
            print(f"\nEXPERIMENT the first road bulldozing kept off the dunes on {n} road domains", flush=True)
        if variant not in RODANTHE_40:
            return
        n = 0
        for g in range(77, 84):
            i = dom.gis_to_pad(g)
            if c.nourishments[i] is not None and c._beach_nourishment_module[i]:
                c.nourishments[i]._overwash_filter = 40.0
                n += 1
        print(f"\nEXPERIMENT Rodanthe overwash filter held at 40% on {n} domains (GIS 77-83)", flush=True)

    managers = {}

    def build_cascade(*a, **k):
        if rule == "storm":
            # The process itself: Cascade's dune_collision_module (cascade/dune_collision.py) inside each
            # domain's storm loop, on the variant's GIS (COLLISION_GIS; the buffers left out), nothing after the year
            dom = HATTERAS_DOMAINS.for_total(k["alongshore_section_count"])
            on = COLLISION_GIS.get(variant, range(1, 91))
            w0 = BEACH_W0_M.get(variant)
            toe_m = [None] * dom.total_domains
            width_m = [None] * dom.total_domains
            if w0:
                # Each cell's median beach width over the surveys, 1998 left out as for the toe
                width = per_gis(st[~st.survey.isin(TOE_SURVEYS_EXCLUDED)].groupby(["gis", "cell"]).beach_width_m.median())
                if variant in WIDTH_FROM_START:
                    width = start_width
                width_collision = width
                if variant in DRY_BEACH:
                    # the dry beach, as HAT_dune_toes.py: 0 where the crest is within one cell of the beach start
                    s0 = sp if variant in WIDTH_FROM_START else st[~st.survey.isin(TOE_SURVEYS_EXCLUDED)]
                    dry = np.where(s0.toe_to_crest_m.notna(), s0.beach_width_m - s0.toe_to_crest_m,
                                   np.where(s0.beach_width_m < 2 * 10, 0.0, np.nan))
                    width_collision = per_gis(s0.assign(dry=dry).groupby(["gis", "cell"]).dry.median())
                    if variant in TOE_DOMAIN_AVG and variant not in CELL_WIDTHS:
                        # one beach factor per domain, the mean of max(0, 1 - W / W0) over its cells in the toe's
                        # surveys, given to the module as the width W0 (1 - factor) (2026-11-04; HAT_dune_toes.py)
                        sw = st[st.survey.isin(TOE_DOMAIN_AVG[variant])]
                        dw = pd.Series(np.where(sw.toe_to_crest_m.notna(), sw.beach_width_m - sw.toe_to_crest_m,
                                                np.where(sw.beach_width_m < 2 * 10, 0.0, np.nan)), index=sw.index)
                        bf = np.clip(1 - dw / w0, 0, None).groupby(sw.gis).mean()
                        width_collision = {g: np.full(50, w0 * (1 - bf.get(g, bf.median()))) for g in range(1, 91)}
                        print()
                        print(f"EXPERIMENT one beach factor per domain, the mean over the {TOE_DOMAIN_AVG[variant]} "
                              f"surveys: island {bf.mean():.2f}", flush=True)
                    if variant in LARSON_ALONE:
                        width_collision = {g: np.zeros(50) for g in range(1, 91)}
                        print(f"EXPERIMENT Larson alone: beach factor 1 everywhere, Cs {cs:.1e}", flush=True)
            for g in on:
                toe_m[dom.gis_to_pad(g)] = toe[g] - MHW_NAVD88_M
                if w0:
                    width_m[dom.gis_to_pad(g)] = width_collision[g]
            # Larson-alone variants zero the widths above (beach factor 1), so they run as before without them
            if (w0 and variant not in LARSON_ALONE) or variant in BEACH_FOLLOWS_SHORELINE or variant in RUNUP_LOSS:
                raise SystemExit(f"{variant}: the module's beach options were removed from cascade/dune_collision.py on "
                                 "2026-11-04; check out git tag collision-with-beach-factor to rerun it")
            k.update(dune_collision_module=True, dune_collision_toe_elevation=toe_m, dune_collision_cs=cs,
                     dune_collision_rows=COLLISION_ROWS.get(variant, "proportional"))
            if variant in ROAD_DUNE_REBUILD:
                gate, trigger = ROAD_DUNE_REBUILD[variant]
                k["road_dune_rebuild_trigger_m"] = gate
                if trigger is not None:
                    k["dune_minimum_elevation"] = [trigger] * len(k["dune_minimum_elevation"])
                print(f"\nEXPERIMENT road dune rebuild: setback gate {gate:g} m, trigger "
                      + ("berm floor" if trigger is None else f"{trigger:.2f} m MHW"), flush=True)
            c = orig_build(*a, **k)
            road_gate(c, dom, width)
            rodanthe_40(c, dom)
            # Per-village dune caps and Rodanthe's overwash-to-dune share, on the village managers (m above the berm)
            for name, (lo, hi) in VILLAGE_DOMAINS.items():
                for g in range(lo, hi + 1):
                    mgr = c.nourishments[dom.gis_to_pad(g)]
                    if mgr is None or not c._beach_nourishment_module[dom.gis_to_pad(g)]:
                        continue
                    if variant in VILLAGE_CAP_VARIANTS:
                        mgr._artificial_maximum_dune_height = VILLAGE_CAPS_M[name] - k["berm_elevation"]
                    if variant in RODANTHE_NO_OW_TO_DUNE and name == "Rodanthe":
                        mgr._overwash_to_dune = 0
            village_filter(c)
            managers.update({g: dom.gis_to_pad(g) for g in on}, cascade=c)
            print(f"\nEXPERIMENT collision: dune_collision_module on GIS {on[0]}-{on[-1]}, Cs {cs:.1e}"
                  + (f", beach W0 {w0:g} m" if w0 else ""), flush=True)
            return c
        c = orig_build(*a, **k)
        village_filter(c)
        dom = HATTERAS_DOMAINS.for_total(k["alongshore_section_count"])
        road_gate(c, dom, per_gis(st[~st.survey.isin(TOE_SURVEYS_EXCLUDED)].groupby(["gis", "cell"]).beach_width_m.median()))
        rodanthe_40(c, dom)
        if rule is None:
            return c
        cls = type(c)
        orig_update = cls.update
        calls = [0]

        # Patched on the class, so the saved model still pickles
        def update(self):
            orig_update(self)
            calls[0] += 1
            year = RUN_START + calls[0]                     # update k runs calendar year 1996 + k
            dom = HATTERAS_DOMAINS.for_total(len(self.barrier3d))
            res = collide(self, dom, toe, loss_for(hours.loc[str(year)]), to_beach)
            if to_beach and self._alongshore_transport_module:
                b3 = self.barrier3d
                n = len(b3)
                self._brie_coupler.update_brie_for_human_modifications(
                    [b3[i].x_t_TS[-1] for i in range(n)], [b3[i].x_s_TS[-1] for i in range(n)],
                    [b3[i].x_b_TS[-1] for i in range(n)], [b3[i].h_b_TS[-1] for i in range(n)],
                    [b3[i].s_sf_TS[-1] for i in range(n)])
            yr = pd.DataFrame([dict(year=year, gis=g, crest_lowering_m=dh, volume_m3_m=vol, shoreline_seaward_m=dx)
                               for g, (dh, vol, dx) in res.items()])
            log.append(yr)
            top = yr.loc[yr.crest_lowering_m.idxmax()]
            print(f"\nEXPERIMENT collision {year}: mean crest lowering {yr.crest_lowering_m.mean():.3f} m, "
                  f"max {top.crest_lowering_m:.2f} m (GIS {int(top.gis)}), sand {yr.volume_m3_m.mean():.2f} m3/m"
                  + (f", shoreline {yr.shoreline_seaward_m.mean():.2f} m seaward" if to_beach else ", lost"),
                  flush=True)

        cls.update = update
        return c

    H.build_cascade = build_cascade
    sys.path.insert(0, str(HINDCAST.parent))
    sys.argv = [str(HINDCAST)]
    try:
        runpy.run_path(str(HINDCAST), run_name="__main__")
    finally:
        if managers:
            # Each domain's manager, read off the model at the end (Cascade.dune_collision)
            collision = managers.pop("cascade").dune_collision
            log = [pd.DataFrame([dict(year=RUN_START + 1 + i, gis=g, crest_lowering_m=dh * 10, volume_m3_m=vol * 1000,
                                      shoreline_seaward_m=0.0)
                                 for g, pad in managers.items()
                                 for i, (dh, vol) in enumerate(zip(collision[pad].crest_lowering_TS,
                                                                   collision[pad].volume_TS))])]
        if log:
            pd.concat(log).round(4).to_csv(EXP_DIR / "logs" / f"collision_{variant}.csv", index=False)
        if managers:
            collision_cell_table(collision, managers, variant)


# The managers' per-cell series (DuneCollisionManager: collision_volume_TS, collision_storms_TS, dune_row_loss_TS,
# collision_deficit_TS) as one row per domain and model year, beside the run's log and in its tables/ folder
def collision_cell_table(collision, pads, variant):
    rows = []
    for g, pad in pads.items():
        m = collision[pad]
        vol = getattr(m, "collision_volume_TS", None)
        if vol is None:
            continue
        storms, loss, deficit = m.collision_storms_TS, m.dune_row_loss_TS, m.collision_deficit_TS
        for t in range(1, len(m.crest_lowering_TS) + 1):
            r = dict(year=RUN_START + t, gis=g, rows=m.rows,
                     volume_m3_m=float(vol[t].mean() * 100),
                     collision_storms_per_cell=float(storms[t].mean()),
                     cells_in_collision=int((storms[t] > 0).sum()),
                     front_rows_collapsed=int(loss[t, :, 0].sum()),
                     back_rows_collapsed=int(loss[t, :, 1:].sum()))
            if deficit is not None:
                r["front_deficit_m3_m"] = float(deficit[t, :, 0].mean() * 100)
                r["back_deficit_m3_m"] = float(deficit[t, :, 1:].mean() * 100)
            rows.append(r)
    if not rows:
        return
    t = pd.DataFrame(rows).round(4)
    t.to_csv(EXP_DIR / "logs" / f"collision_cells_{variant}.csv", index=False)
    rd = run_dir(variant)
    if rd is not None:
        t.to_csv(rd / "tables" / "collision_cells.csv", index=False)
    print(f"EXPERIMENT collision per cell: {t.volume_m3_m.sum() / t.gis.nunique():.1f} m3/m eroded per domain, "
          f"{t.front_rows_collapsed.sum()} front and {t.back_rows_collapsed.sum()} back rows collapsed", flush=True)


# Where and when the dune cells collapse, and the collision volume, one figure per run (action "collapse")
def plot_collapse(variants):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from site_layer import hat_figure_style as S
    S.apply_style()
    for v in variants:
        f = EXP_DIR / "logs" / f"collision_cells_{v}.csv"
        if not f.exists():
            print(f"  no {f.name}: run {v} first")
            continue
        t = pd.read_csv(f)
        years = sorted(t.year.unique())
        fig, axs = plt.subplots(1, 2, figsize=(S.FIG_W_DOUBLE, 4.6), sharey=True, layout="constrained")
        for ax, col, cmap, lab in ((axs[0], "volume_m3_m", "Purples", "Collision volume eroded (m³/m)"),
                                   (axs[1], "front_rows_collapsed", "Oranges", "Front-row cells collapsed (of 50)")):
            grid = t.pivot_table(index="year", columns="gis", values=col).reindex(index=years, columns=range(1, 91))
            im = ax.imshow(grid.values, aspect="auto", origin="lower", cmap=cmap, interpolation="nearest",
                           extent=(0.5, 90.5, years[0] - 0.5, years[-1] + 0.5), vmin=0)
            fig.colorbar(im, ax=ax, orientation="horizontal", location="bottom", shrink=0.8, label=lab)
            ax.set_xlabel(S.DOMAIN_AXIS_LABEL)
        axs[0].set_ylabel("Model year")
        for k, ax in enumerate(axs):
            ax.text(-0.02, 1.02, f"({'ab'[k]})", transform=ax.transAxes, fontsize=10, fontweight="bold", ha="right")
        out = EXP_DIR / f"collision_collapse_{v}.png"
        S.save(fig, out, bbox_inches="tight", close=True)
        S.record_caption(out, (f"Collision-regime dune erosion per GIS domain and model year, run {v} (rows mode "
                               f"{t.rows.iloc[0]}): (a) the volume the collision rule eroded, domain mean (m3/m); "
                               "(b) front-row dune cells eroded through to the toe (collapsed), of the domain's 50."))
        print(f"wrote {out}")


def run_dir(variant):
    root = CONTROL_DIR if variant == "control" else EXP_DIR / "runs" / variant
    base = root / f"{RUN_START}_{RUN_END}" / PRESET
    # the newest run: run names changed on 2026-10-29, so a variant's folder can hold both
    meta = sorted(base.glob("HAT_*/*_run_metadata.json"), key=lambda f: f.stat().st_mtime) if base.exists() else []
    return meta[-1].parent if meta else None


# One run in its own process, logged, with its launch record
def launch(variant):
    logs = EXP_DIR / "logs"
    logs.mkdir(parents=True, exist_ok=True)
    env = dict(os.environ)
    env.update({"HAT_IGNORE_SETTINGS": "1", "HAT_START_YEAR": str(RUN_START),
                "HAT_SOURCE_SINK_PRESET": PRESET, "HAT_SCENARIO": SCENARIO, "HAT_GROIN_ENABLED": "1",
                "HAT_RUN_KIND": "dune-storm-response", "HAT_RUN_TAG": f"{TAG}/runs/{variant}",
                "HAT_SAVE_MODEL_STATE": "1", "HAT_MAKE_GIFS": "0", "HAT_SHOW_FIGURES": "0",
                "HAT_OVERWRITE": "1", "PYTHONIOENCODING": "utf-8"})
    rule, cs, to_beach, growth = VARIANTS[variant]
    # The hindcast turns the module on by default since 2026-10-29; here only the "storm" variants have it, through
    # their own toes and Cs (_launch), and every other variant runs without it, as it was made
    env["HAT_DUNE_COLLISION"] = "1" if rule == "storm" else "0"
    if growth:
        per = growth if isinstance(growth, dict) else {g: growth for g in range(1, 91)}
        env["HAT_GROWTH_OVERRIDE"] = ",".join(f"{g}={lo:.3f}:{hi:.3f}" for g, (lo, hi) in sorted(per.items()))
    sink = BUXTON_SINK if variant in SINK_VARIANTS else None
    overrides = {**(sink or {}), **ENDS.get(variant, {})}
    if overrides:
        env["HAT_BE_OVERRIDE"] = ",".join(f"{g}={r}" for g, r in overrides.items())
    if variant in HS_M:
        env["HAT_HS"] = str(HS_M[variant])
    env["HAT_GROIN_STRENGTH"] = str(GROIN_STRENGTH.get(variant, GROIN_STRENGTH_BASE))
    t0 = time.time()
    with open(logs / f"{variant}.log", "w", encoding="utf-8") as fh:
        p = subprocess.run([sys.executable, str(_HERE), "_launch", variant],
                           env=env, cwd=HINDCAST.parent, stdout=fh, stderr=subprocess.STDOUT)
    rec = dict(variant=variant, rule=rule, **({"min_twl_m": cs} if rule == "width" else {"cs": cs}),
               sand_to_beach=to_beach, growth=growth, sink=sink, beach_w0_m=BEACH_W0_M.get(variant),
               hs_m=HS_M.get(variant), ends=ENDS.get(variant), groin_strength=GROIN_STRENGTH.get(variant, GROIN_STRENGTH_BASE),
               returncode=p.returncode, minutes=round((time.time() - t0) / 60, 1))
    with open(logs / "launches.jsonl", "a", encoding="utf-8") as fh:
        fh.write(json.dumps(rec) + "\n")
    print(f"{variant:14s} exit {p.returncode}  {rec['minutes']} min", flush=True)


# Island skill, per-GIS LRR against CoastSat, and the modelled crests against each survey's
def score():
    from cascade_pipeline.coastsat_lowess import build_coastsat_series
    from cascade_pipeline.hindcast import build_target_table
    from cascade_pipeline.saved_runs import load_saved_cascade
    from site_layer.hatteras_site_config import HATTERAS_DOMAINS as DOM
    sys.path.insert(0, str(_HAT_DIR / "figure_making" / "model_output"))
    from rerender_run_figures import COASTSAT_DATASETS, LOWESS_CONFIG
    from run_figures import TARGET_WINDOW

    cs = build_coastsat_series(COASTSAT_DATASETS, active_period_start=RUN_START,
                               lowess_config=LOWESS_CONFIG, domains=DOM)
    target = build_target_table(next(c for c in cs if c["active"]), LOWESS_CONFIG, DOM, TARGET_WINDOW)
    target = target.set_index("gis_domain")["target_lrr_m_yr"]
    sv = pd.read_csv(SURVEY_CRESTS)
    sv = sv[sv.status == "ok"].copy()
    sv["row"] = pd.to_datetime(sv.model_state).dt.year - (RUN_START + 1)   # row 0: 1 Jan of the first model year

    frames, summary, crests = {"coastsat": pd.DataFrame({"lrr": target})}, [], []
    for v in ["control", *VARIANTS]:
        rd = run_dir(v)
        if rd is None:
            print(f"  missing {v}")
            continue
        c = load_saved_cascade(rd / f"{rd.name}.npz")
        dom = DOM.for_total(len(c.barrier3d))
        rate = pd.read_csv(rd / "tables" / "shoreline_change_rate.csv").set_index("gis_domain")
        ow = np.load(next(rd.glob("*_overwash.npy")))
        meta = json.loads((rd / f"{rd.name}_run_metadata.json").read_text(encoding="utf-8"))
        model = {}
        for g in range(1, 91):
            b = c.barrier3d[dom.gis_to_pad(g)]
            model[g] = (np.asarray(b.DuneDomain, float).max(axis=2) + b._BermEl) * 10 + MHW_NAVD88_M
        frames[v] = pd.DataFrame([dict(gis=g, lrr=rate.lrr_m_yr[g], overwash=ow[1:, dom.gis_to_pad(g)].sum())
                                  for g in range(1, 91)]).set_index("gis")
        # Crests in each survey's model state, cell for cell
        for (survey, row), grp in sv.groupby(["survey", "row"]):
            if row < 0 or row >= len(model[1]):
                continue
            m = np.array([model[g][row][cell] for g, cell in zip(grp.gis, grp.cell)])
            o = grp.crest_m_navd88.to_numpy()
            mg = pd.Series(m).groupby(grp.gis.to_numpy()).median()
            og = pd.Series(o).groupby(grp.gis.to_numpy()).median()
            crests.append(dict(run=v, survey=survey, model_median=np.median(m), lidar_median=np.median(o),
                               model_p10=np.percentile(m, 10), lidar_p10=np.percentile(o, 10),
                               model_below_4m=np.mean(m < 4), lidar_below_4m=np.mean(o < 4),
                               model_stable_median=np.median(m[grp.gis.isin(STABLE_GIS).to_numpy()]),
                               lidar_stable_median=np.median(o[grp.gis.isin(STABLE_GIS).to_numpy()]),
                               r_domains=np.corrcoef(mg, og)[0, 1]))
        lrr = frames[v].lrr
        sk = meta["skill"]
        summary.append(dict(run=v, bias=float(sk["mean_bias_interior_m_yr"]), rmse=float(sk["rmse_interior_m_yr"]),
                            rmse_gis_8_15=float(np.sqrt(((lrr.loc[8:15] - target.loc[8:15]) ** 2).mean())),
                            rmse_gis_80_89=float(np.sqrt(((lrr.loc[80:89] - target.loc[80:89]) ** 2).mean())),
                            island_overwash=frames[v].overwash.loc[2:89].sum()))
    pd.set_option("display.width", 250)
    pd.set_option("display.max_columns", 40)
    crests = pd.DataFrame(crests)
    # Dunes: the island's median crest over the surveys after the start, and the stable dunes in 2009
    later = crests[crests.survey != crests.survey.min()]
    dunes = pd.DataFrame({
        "crest_median_rmse": later.groupby("run").apply(
            lambda x: np.sqrt(((x.model_median - x.lidar_median) ** 2).mean()), include_groups=False),
        "stable_2009_minus_lidar": crests[crests.survey == crests.survey.max()].set_index("run").eval(
            "model_stable_median - lidar_stable_median")})
    summary = pd.DataFrame(summary).set_index("run").join(dunes)
    out = pd.concat(frames, axis=1)
    print(summary.round(3))
    print(out.xs("lrr", axis=1, level=1).loc[list(REPORT_GIS) + list(range(80, 90))].round(2))
    print(crests.pivot(index="survey", columns="run", values="r_domains").round(2))
    print(crests.pivot(index="survey", columns="run", values="model_p10").round(2).assign(
        lidar=crests.groupby("survey").lidar_p10.first().round(2)))
    summary.to_csv(EXP_DIR / "collision_runs_summary.csv")
    out.to_csv(EXP_DIR / "collision_runs_by_gis.csv")
    crests.round(4).to_csv(EXP_DIR / "collision_runs_crests.csv", index=False)
    print(f"\nwrote {EXP_DIR / 'collision_runs_summary.csv'}")


# Figures

# Label, colour and line width per run, all solid: hue for the rule (purple Larson, orange strong, green beach width; blue is
# CoastSat), light for the sand
# lost and dark for the sand on the beach; drawn in this order, so later runs sit on top where they coincide
RUN_STYLES = {
    "control": ("No collision (control)", "#8C8C8C", 1.6),
    "larson_lost": ("Cs 1.7e-4, sand lost", "#BCBDDC", 1.6),
    "larson_beach": ("Cs 1.7e-4, sand to beach", "#6A51A3", 1.8),
    "strong_lost": ("Cs 4e-4, sand lost", "#FDAE6B", 1.6),
    "strong_beach": ("Cs 4e-4, sand to beach", "#D94801", 1.8),
    "width_lost": ("Beach width, sand lost", "#A1D99B", 1.6),
    "width_beach": ("Beach width, sand to beach", "#238B45", 1.8),
    "width_storms_beach": ("Beach width, storm years, sand to beach", "#00441B", 1.8),
    "width_storms_r04": ("Storm years, growth r 0.03-0.05", "#FCC5C0", 1.6),
    "width_storms_r06": ("Storm years, growth r 0.05-0.07", "#F768A1", 1.7),
    "width_storms_r08": ("Storm years, growth r 0.07-0.09", "#AE017E", 1.8),
    "width_storms_r10": ("Storm years, growth r 0.09-0.11", "#49006A", 1.8),
    "storm_cs2e5": ("Per storm, Cs 2.0e-5", "#80CDC1", 1.6),
    "storm_cs4e5": ("Per storm, Cs 4.2e-5", "#DFC27D", 1.6),
    "storm_cs2e5_r06": ("Per storm, Cs 2.0e-5, growth r 0.05-0.07", "#01665E", 1.8),
    "storm_cs4e5_r06": ("Per storm, Cs 4.2e-5, growth r 0.05-0.07", "#8C510A", 1.8),
    "sink": ("No collision, sink at GIS 9-14", "#F4A582", 1.6),
    "storm_cs2e5_r06_sink": ("Per storm, Cs 2.0e-5, growth r 0.05-0.07, sink at GIS 9-14", "#B2182B", 1.8),
    "storm_w45": ("Per storm with beach, Cs 1.24e-4", "#C2A5CF", 1.6),
    "storm_w45_r04": ("Per storm with beach, Cs 1.24e-4, growth r 0.03-0.05", "#9970AB", 1.7),
    "storm_w45_r06": ("Per storm with beach, Cs 1.24e-4, growth r 0.05-0.07", "#762A83", 1.8),
    "control_ss1": ("No collision, 1997 storms on", "#4D4D4D", 1.6),
    "storm_w45_r04_ss1": ("Per storm with beach, growth r 0.03-0.05, 1997 storms on", "#E7298A", 1.8),
    "control_hs3": ("No collision, Hs 3.0 m", "#FDB863", 1.6),
    "storm_w45_r04_hs3": ("Per storm with beach, growth r 0.03-0.05, Hs 3.0 m", "#B35806", 1.8),
    "control_hs3_ends": ("No collision, Hs 3.0 m, ends re-solved", "#BABABA", 1.6),
    "storm_w45_r04_hs3_ends": ("Per storm with beach, growth r 0.03-0.05, Hs 3.0 m, ends re-solved", "#542788", 1.8),
    "control_hs2p5": ("No collision, Hs 2.5 m", "#D9D9D9", 1.6),
    "storm_w45_r04_hs2p5": ("Per storm with beach, growth r 0.03-0.05, Hs 2.5 m", "#1B7837", 1.8),
    "storm_w45_r04_hs2p5_g04": ("Hs 2.5, groin 0.04", "#D9F0D3", 1.6),
    "storm_w45_r04_hs2p5_g06": ("Hs 2.5, groin 0.06", "#A6DBA0", 1.6),
    "storm_w45_r04_hs2p5_g08": ("Hs 2.5, groin 0.08", "#5AAE61", 1.7),
    "storm_w45_r04_hs3_g04": ("Hs 3.0, groin 0.04", "#E7D4E8", 1.6),
    "storm_w45_r04_hs3_g06": ("Hs 3.0, groin 0.06", "#C2A5CF", 1.6),
    "storm_w45_r04_hs3_g08": ("Hs 3.0, groin 0.08", "#9970AB", 1.7),
    "storm_w45_r04_hs2p5_g13": ("Hs 2.5, groin 0.13", "#FDAE61", 1.6),
    "storm_w45_r04_hs2p5_g16": ("Hs 2.5, groin 0.16", "#F46D43", 1.7),
    "storm_w45_r04_hs2p5_g19": ("Hs 2.5, groin 0.19", "#A50026", 1.8),
    "storm_w45_r04_hs3_g13": ("Hs 3.0, groin 0.13", "#80CDC1", 1.6),
    "storm_w45_r04_hs3_g16": ("Hs 3.0, groin 0.16", "#35978F", 1.7),
    "storm_w45_r04_hs3_g19": ("Hs 3.0, groin 0.19", "#01665E", 1.8),
    "storm_cs145_r04_hs3_ends": ("Cs 1.45e-4, Hs 3.0 m, ends re-solved", "#E08214", 1.8),
    "storm_cs145_r04_hs2p5": ("Cs 1.45e-4, Hs 2.5 m", "#C51B7D", 1.8),
    # The 2026-10-29 calibration (village dune cap 8.0 m); drawn against cal_control, Hs in purples off
    # CoastSat's blue
    "cal_control": ("No collision, Hs 3.0 m", "#8C8C8C", 1.6),
    "cal_automorph_cs10": ("Collision, Automorph toes, Cs 1e-5", "#7B3F99", 1.8),
    "cal_automorph_module_cs10": ("Collision module on (Automorph toes, Cs 1e-5)", "#7B3F99", 1.8),
    "cal_step2_growth": ("Collision module off", "#8C8C8C", 1.6),
    "cal_automorph_module_cs30": ("Collision module on, Cs 3e-5", "#74A9CF", 1.6),
    "cal_automorph_module_cs50": ("Collision module on, Cs 5e-5", "#31A354", 1.6),
    "cal_automorph_module_cs100": ("Collision module on, Cs 1e-4", "#FD8D3C", 1.6),
    "cal_automorph_module_cs200": ("Collision module on, Cs 2e-4", "#D7301F", 1.6),
    "cal_automorph_cellsall_cs10": ("Toe per cell, Cs 1e-5", "#7B3F99", 1.8),
    "cal_automorph_cellsall_cs30": ("Toe per cell, Cs 3e-5", "#74A9CF", 1.6),
    "cal_automorph_cellsall_cs50": ("Toe per cell, Cs 5e-5", "#31A354", 1.6),
    "cal_automorph_cellsall_cs100": ("Toe per cell, Cs 1e-4", "#FD8D3C", 1.6),
    "cal_automorph_cellsall_cs200": ("Toe per cell, Cs 2e-4", "#D7301F", 1.6),
    "cal_automorph_cells_9_15": ("Collision module on, a toe per cell, GIS 9-15", "#D7301F", 1.6),
    "cal_automorph_cellmed_9_15": ("Collision module on, the cells' domain median, GIS 9-15", "#2B8CBE", 1.6),
    "cal_hs2": ("Hs 2.0 m", "#CBC9E2", 1.6),
    "cal_hs2p5": ("Hs 2.5 m", "#9E9AC8", 1.6),
    "cal_r03": ("Hs 3.0 m (growth 0.02-0.04, groin 0.10)", "#756BB1", 1.8),
    "cal_hs3p5": ("Hs 3.5 m", "#54278F", 1.7),
    "cal_hs4": ("Hs 4.0 m", "#3F007D", 1.8),
    "cal_g13": ("Hs 3.0, groin 0.13", "#99D8C9", 1.6),
    "cal_g16": ("Hs 3.0, groin 0.16", "#41AE76", 1.7),
    "cal_g19": ("Hs 3.0, groin 0.19", "#006D2C", 1.8),
    "cal_g22": ("Hs 3.0, groin 0.22", "#00441B", 1.8),
    "cal_control_g16": ("No collision", "#7f7f7f", 1.6),
    "cal_hs2p75": ("Hs 2.75 m", "#E7298A", 1.8),
    "cal_g16_nocape": ("Current calibration (collision, off at Cape Point)", "#7b3294", 2.0),
    "cal_g16_capeg05": ("Collision, cape growth 0.05-0.10", "#E7298A", 1.7),
    "cal_g16_capeg10": ("Collision, cape growth 0.10-0.20", "#7570B3", 1.9),
    "cal_g16_nocape_vcap": ("Per-village dune caps", "#7B3294", 1.8),
    "cal_g16_nocape_vcap_ow0": ("Per-village caps, no overwash onto Rodanthe dunes", "#E66101", 1.8),
}
# (file tag, GIS range or None for the whole island, y half-range in m)
VIEWS = (("island", None, 100.0), ("buxton", (1, 25), 60.0), ("avon_salvo_road", (45, 60), 40.0), ("rodanthe", (70, 90), 60.0))
# Map features with nothing to point at in a zoomed view
ZOOM_HIDDEN_FEATURES = {"Shoal influence zone", "Village center", "Pier"}
# Per view where it differs: the Avon-Salvo road has no village, pier or groin in view
ZOOM_HIDDEN_BY_VIEW = {"avon_salvo_road": {"Shoal influence zone", "Village center", "Pier", "Groin", "Community"}}


# Every run's end-minus-start position change against the CoastSat net change, in the layout of the run's
# shoreline_position_change figure (net_change_comparison.py), as the forced-dune-crests plot does
def plot(variants=None, base="control"):
    import dataclasses
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import matplotlib.ticker as ticker
    from matplotlib.lines import Line2D
    from cascade_pipeline.coastsat_lowess import CoastSatDataset, build_coastsat_series
    from cascade_pipeline.plotting.comparison_style import restyle
    from cascade_pipeline.plotting.net_change_comparison import LOWESS_CONFIG, smooth_model
    from cascade_pipeline.plotting.rate_comparison import (
        DEFAULT_RATE_COMPARISON, plot_annotated_rate_comparison)
    from cascade_pipeline.run_info import RunInfo
    from cascade_pipeline.run_layout import resolve
    from cascade_pipeline.shoreline import compute_change_rate
    from site_layer.hat_figure_style import mark_offaxis, offaxis_clause, record_caption
    from site_layer.hat_observed_rates import NET_CHANGE_WINDOWS, WINDOW_ROLE, net_change_dir
    from site_layer.hatteras_site_config import HATTERAS_ANNOTATIONS, HATTERAS_DOMAINS

    # End minus start (m, seaward positive), smoothed as the run figure smooths it
    runs = {}
    for v in [base, *(variants or [v for v in RUN_STYLES if v != base and not v.startswith("cal_")])]:
        rd = run_dir(v)
        if rd is None:
            print(f"  missing {v}")
            continue
        shoreline_m = np.load(resolve(rd, "matrix", rd.name))
        runs[v] = (rd, smooth_model(compute_change_rate(shoreline_m, span_years=1, flip_sign=True)))
    if base not in runs:
        print("  the control run is needed to draw the figures")
        return

    tag = f"{RUN_START}_{RUN_END}"
    cs_series = build_coastsat_series(
        [CoastSatDataset(label=f"CoastSat net change ({RUN_START}–{RUN_END})", period_start=RUN_START,
                         csv_path=str(net_change_dir(RUN_START, RUN_END) / f"transect_net_change_{tag}.csv"),
                         rate_col="net_change_m")],
        active_period_start=RUN_START, lowess_config=LOWESS_CONFIG, domains=HATTERAS_DOMAINS)

    control_dir, control = runs[base]
    label_c, color_c, _ = RUN_STYLES[base]
    gis = np.arange(HATTERAS_DOMAINS.first_gis_id, HATTERAS_DOMAINS.last_gis_id + 1)
    geometry = HATTERAS_DOMAINS.for_total(len(control))     # the runs' own buffers, not the config's
    real = slice(geometry.start_real_index, geometry.end_real_index)
    role = WINDOW_ROLE.get((RUN_START, RUN_END), "").replace(" (held out)", "")

    for view, gis_range, y_half in VIEWS:
        where = f"GIS {gis_range[0]}–{gis_range[1]}" if gis_range else "whole island"
        subtitle = " · ".join(p for p in (role, where, "full management + groin, edgeBE",
                                          "collision-regime dune erosion") if p)
        config = dataclasses.replace(
            DEFAULT_RATE_COMPARISON, quantity="position", publication_text=True,
            show_features=True, plot_domain_means=False, raw_lrr_southern_only=True,
            smoothing_marker=False, observed_label="CoastSat (7-domain LOWESS)", model_label=label_c,
            title=f"Shoreline position change with dune collision erosion, {RUN_START}–{RUN_END}",
            ylim=(-y_half, y_half), subtitle=subtitle)
        run = RunInfo(run_name=control_dir.name, run_dir=str(control_dir), start_year=RUN_START,
                      end_year=RUN_END, Hs=2.0, flip_sign_model=True)
        # The control is drawn by the run figure's own call, in the control's colour
        fig, ax = plot_annotated_rate_comparison(
            control, cs_series, run, estimator=None, domains=HATTERAS_DOMAINS,
            annotations=dataclasses.replace(HATTERAS_ANNOTATIONS, model_color=color_c),
            lowess_config=LOWESS_CONFIG, config=config)

        in_view = np.ones_like(gis, bool) if gis_range is None else (
            (gis >= gis_range[0]) & (gis <= gis_range[1]))
        named = [(label_c, mark_offaxis(ax, gis[in_view], control[real][in_view], y_half, color=color_c))]
        handles = []
        for v, (rd, change) in runs.items():
            label, color, lw = RUN_STYLES[v]
            handles.append(Line2D([0], [0], color=color, lw=lw + 0.4, label=label))
            if v == base:
                # Drawn by the run figure's call; given its width here, under the collision runs
                for line in ax.get_lines():
                    if line.get_label() == label_c:
                        line.set_linewidth(lw)
                continue
            ax.plot(gis, change[real], color=color, lw=lw, zorder=6)
            named.append((label, mark_offaxis(ax, gis[in_view], change[real][in_view], y_half, color=color)))

        if gis_range is not None:
            # Zoom: hide the labels of features outside the view; centre the transects on their domain
            ax.set_xlim(gis_range[0] - 0.5, gis_range[1] + 0.5)
            ax.xaxis.set_major_locator(ticker.MultipleLocator(5))
            ax.xaxis.set_minor_locator(ticker.MultipleLocator(1))
            for txt in ax.texts:
                x_is_data = txt.get_transform() not in (ax.transAxes, fig.transFigure)
                if x_is_data and not gis_range[0] - 0.5 <= txt.get_position()[0] <= gis_range[1] + 0.5:
                    txt.set_visible(False)
            for coll in ax.collections:
                if "transect" in str(coll.get_label()):
                    coll.set_offsets(coll.get_offsets() - [0.5, 0.0])

        # The figure's legend with every run in place of the one model line, then the house finish
        old = fig.legends.pop()
        hidden = {label_c} | (set() if gis_range is None else ZOOM_HIDDEN_BY_VIEW.get(view, ZOOM_HIDDEN_FEATURES))
        entries = [(h, t.get_text()) for h, t in zip(old.legend_handles, old.get_texts())
                   if t.get_text() not in hidden]
        fig.legend(handles + [h for h, _ in entries],
                   [h.get_label() for h in handles] + [t for _, t in entries],
                   loc="outside lower center")
        restyle(fig, ax, subtitle)

        out = EXP_DIR / f"collision_runs_position_{view}{'__' + '_'.join(variants) if variants else ''}.png"
        fig.savefig(out, dpi=300, bbox_inches="tight", facecolor="white")
        plt.close("all")
        (s0, s1), (e0, e1) = NET_CHANGE_WINDOWS[(RUN_START, RUN_END)]
        record_caption(out, (
            f"Shoreline position change, {RUN_START}–{RUN_END}, {where}; positive is seaward. Modelled "
            f"(CASCADE): the shoreline at the end of each run minus the start, smoothed with a 7-domain "
            f"LOWESS north of domain 10 as in the run figure, for the full-management edgeBE run with the "
            f"groin without collision erosion (grey) and with every dune cell lowered after each model year "
            f"by Cs·4Σ(TWL − toe)²Δt/T / 20 m over the hours the total water level stood between its lidar toe "
            f"and its crest (Larson et al. 2004): Cs 1.7e-4, the Dennis fit (purple), or 4e-4, what the Dennis "
            f"hotspots need (orange); or, green, in each year whose highest water stood between the toe and the "
            f"crest, by a·max(0, 1 − W/W0) of the crest elevation, W the last survey's beach width (the fit on "
            f"Dennis and 2004–05; darkest, only in years whose highest water reached 3.1 m NAVD88); light, the "
            f"eroded sand lost; dark, placed on the beach as CASCADE's shoreface nourishment places fill. Runs "
            f"drawn: {', '.join(RUN_STYLES[v][0] for v in runs)}. Where runs coincide the last drawn sits on top. "
            f"Observed "
            f"(CoastSat, blue dots and line): per transect, the mean position over {e0} to {e1} minus the "
            f"mean over {s0} to {s1}, smoothed with a 7-domain LOWESS north of domain 10; dots south of it. "
            f"The y axis is fixed at ±{y_half:g} m."
            + offaxis_clause([(k, n) for k, n in named if n], y_half)
            + f" Experiment {TAG}."))
        print(f"wrote {out}")


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("action", choices=["toes", "twl", "calibrate", "empirical", "run", "score", "plot", "collapse",
                                       "_launch"])
    ap.add_argument("variant", nargs="*")
    args = ap.parse_args()
    if args.action == "_launch":
        _launch(args.variant[0])
    elif args.action == "run":
        with ThreadPoolExecutor(WORKERS) as ex:
            list(ex.map(launch, args.variant or list(VARIANTS)))
    elif args.action == "collapse":
        plot_collapse(args.variant)
    elif args.action == "plot":
        # The calibration runs (cal_*) are on the 2026-10-29 defaults; their baseline is cal_control
        cal = any(v.startswith("cal_") for v in args.variant)
        plot([v for v in args.variant if v != "cal_control"] or None, base="cal_control" if cal else "control")
    else:
        {"toes": toes, "twl": twl, "calibrate": calibrate, "empirical": empirical, "score": score}[args.action]()


if __name__ == "__main__":
    main()
