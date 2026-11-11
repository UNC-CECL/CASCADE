"""
Does the 1996-2009 calibration of the dune collision module hold in 2009-2025, with nothing refit?

    python scripts/real-world-cases/Hatteras/dune-collision/4-hindcast/2-test-2009-2025/HAT_collision_test_2009.py run [variant ...]
    python scripts/real-world-cases/Hatteras/dune-collision/4-hindcast/2-test-2009-2025/HAT_collision_test_2009.py score

The test window run on the hindcast's own defaults (managed edgeBE, groin 0.16 until the 2017 Buxton fill, Hs 3.0
on the 2009 ends solved at Hs 3.0, dune growth 0.02-0.04, village dune cap 8.0 m NAVD88, the 2014, 2017 and 2022
nourishments) with the module as calibrated on 1996-2009 (Cs 1.24e-4, beach limit 45 m; 2026-10-29-collision-
calibration), on the 2009 start's own toes and beach widths. Variants:

    control      the module off
    collision    the module on everywhere (the hindcast default)
    nocape       the module off at Cape Point, south of the Buxton groin (GIS 1-5), as in the calibration
    nocape_rNN   nocape with faster dune growth, NN its mean (r +/- 0.01): do the 2009-2025 dunes grow faster?
    *_hsX        control and nocape at another wave height (alongshore diffusivity), on approximate ends (ENDS_BY_HS)
    nocape_dynw  nocape with the beach following the modelled shoreline; nocape_dynw_sw from the start widths
    nocape_dry, nocape_rd198_dry  nocape and nocape_rd198 on the dry beach (beach start to toe), Cs 1.02e-4, limit
                 25 m (2026-11-03). The variants before them ran on the beach to the crest, Cs 1.24e-4, limit 45 m;
                 rerun now they would pick up the dry inputs and the new defaults
    nocape_toestart  nocape on toes and beach widths from the 2009 start profile alone (the default median is of
                 the 2014-2019 surveys, all after Irene and Sandy)

run     each variant in its own process (HAT_hindcast.py, its run figures included), under
        output/raw_runs/dune-storm-response/<TAG>/runs/<variant>/
score   shoreline LRR RMSE against CoastSat by zone and island-wide, and the dune crests against the 2009-start lidar
        surveys (share of domain-survey cells within 0.5 m, and over 1 m too high). Writes tables/scores.csv

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-10-29
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
_HAT_DIR = next(_p for _p in _HERE.parents if (_p / "hatteras_paths.py").exists())
if str(_HAT_DIR) not in sys.path:
    sys.path.insert(0, str(_HAT_DIR))
from hatteras_paths import WORKSPACE as PROJECT_ROOT  # noqa: E402

# --- CONFIG ------------------------------------------------------------------
TAG = "2026-10-29-collision-test-2009-2025"
EXP_DIR = PROJECT_ROOT / "output" / "raw_runs" / "dune-storm-response" / TAG
HINDCAST = _HAT_DIR / "hatteras_ms" / "HAT_hindcast.py"
RUN_START = 2009
PRESET, SCENARIO = "edgeBE", "full_management"
# The module on (True) or off, and the GIS it leaves out
VARIANTS = {"control": (False, ()), "collision": (True, ()), "nocape": (True, tuple(range(1, 6))),
            **{f"nocape_r{round(r * 100):02d}": (True, tuple(range(1, 6))) for r in (0.04, 0.05, 0.06, 0.08)}}
# Wave height per variant (the hindcast's 3.0 otherwise), each on ends for that Hs. Only Hs 3.0 is solved on the
# current settings (+51.0016, +10.3255). Hs 2.0 is the 2026-10-21 solve on the earlier settings, 2.5 halfway
# between, and 1.5, 3.5 and 4.0 linear extensions of the 2.0-3.0 line: a quick look, not solved ends
ENDS_BY_HS = {1.5: (70.8214, 6.1578), 2.0: (64.2148, 7.5470), 2.5: (57.6082, 8.9363), 3.5: (37.7884, 13.1040),
              4.0: (24.5752, 15.8825)}
for _hs in ENDS_BY_HS:
    for _base in ("control", "nocape"):
        VARIANTS[f"{_base}_hs{_hs:g}".replace(".", "p")] = VARIANTS[_base]
HS_OF = {f"{b}_hs{h:g}".replace(".", "p"): h for h in ENDS_BY_HS for b in ("control", "nocape")}
# Dune growth bounds per variant, island-wide; the hindcast's 0.02-0.04 otherwise
GROWTH = {f"nocape_r{round(r * 100):02d}": (round(r - 0.01, 3), round(r + 0.01, 3)) for r in (0.04, 0.05, 0.06, 0.08)}
# Toes and beach widths from the 2009 start profile alone, not the median of the 2014-2019 surveys (all after
# Irene and Sandy): the same screens, the domain median where a cell has none
VARIANTS["nocape_toestart"] = VARIANTS["nocape_toestart_cs250"] = VARIANTS["nocape"]
TOES_FROM_START = ("nocape_toestart", "nocape_toestart_cs250")
# The beach following the modelled shoreline (DuneCollisionManager beach_follows_shoreline), on the default toes:
# from the default median widths, then from the 2009 start profile's widths
VARIANTS["nocape_dynw"] = VARIANTS["nocape_dynw_sw"] = VARIANTS["nocape"]
BEACH_FOLLOWS_SHORELINE = ("nocape_dynw", "nocape_dynw_sw")
# The control without the 2014, 2017 and 2022 fills (the hindcast's full_no_fill scenario): against control, where
# the fill sand goes
VARIANTS["control_nofill"] = VARIANTS["control"]
SCENARIO_OF = {"control_nofill": "full_no_fill"}
# The road's dune rebuild from the lidar rebuilds (tables/observed_dune_rebuilds.csv), on nocape: the setback gate
# (rebuilt up to 140 m behind the dune), then with the trigger height too (crest before a rebuild, median 2.66 m
# NAVD88); (gate m, trigger m MHW or None for the hindcast's berm floor)
VARIANTS["nocape_rd150"] = VARIANTS["nocape_rd150h"] = VARIANTS["nocape"]
ROAD_DUNE_REBUILD = {"nocape_rd150": (150.0, None), "nocape_rd150h": (150.0, 2.66 - 0.36)}   # 0.36: MHW, m NAVD88
# A 198 m shoreline-to-road gate in place of the fitted gate: the Bonner Bridge Replacement FEIS (NCDOT/FHWA 2008,
# p. 2-128) lists a shoreline within 650 ft (198 m) of the road as a monitoring warning sign for bridge phasing (from
# Isabel's washover fans), not as a dune maintenance rule; NCDOT's criterion is 230 ft (70 m) from MHW to the
# pavement (p. 2-64; input_prep/4-mgmt-forcings/road_dune_rebuild/NCDOT_road_dune_criteria.txt). The road manager
# measures the road from the dune, so each domain's gate is 198 m less its median beach width
VARIANTS["nocape_rd198"] = VARIANTS["nocape"]
# ... with the beach acting on the water instead of the volume (RUNUP_LOSS, Cs 1.64e-4), and the ramp at Cs 1.45e-4:
# each at its offline fit with the toe floor (HAT_collision_rule_comparison.py)
VARIANTS["nocape_rd198_rl02"] = VARIANTS["nocape_rd198_cs145"] = VARIANTS["nocape"]
# The module on the dry beach (beach start to toe, HAT_dune_toes.py since 2026-11-03), Cs 1.02e-4 and limit 25 m as
# refit on 1996-2009 (rule comparison, then cal_g16_nocape_rd198_dry); the variants above ran before the switch,
# on the beach to the crest, Cs 1.24e-4 and 45 m. The road gate stays on the beach to the crest
VARIANTS["nocape_dry"] = VARIANTS["nocape_rd198_dry"] = VARIANTS["nocape"]
DRY_BEACH, DRY_CS, DRY_LIMIT_M = ("nocape_dry", "nocape_rd198_dry"), 1.02e-4, 25.0
OLD_CS, OLD_LIMIT_M = 1.24e-4, 45.0
# nocape_rd198_dry with the eroded dune sand put on the beach and shoreface, not deleted (2026-11-03): after each
# storm's lowering, the domain's eroded volume is placed as CASCADE places a nourishment (shoreface_nourishment,
# Ashton & Lorenzo-Trueba 2018), the shoreline moving seaward by 2V / (2 h_b + D_shoreface). An experiment here, not
# yet a DuneCollisionManager option
VARIANTS["nocape_rd198_dry_sand"] = VARIANTS["nocape"]
SAND_TO_BEACH = {"nocape_rd198_dry_sand": 1.0}   # share of the eroded volume put back
DRY_BEACH = DRY_BEACH + ("nocape_rd198_dry_sand",)
# The control with the 198 m road rule and the module off: against nocape_rd198_dry, the module alone (2026-11-03)
VARIANTS["control_rd198"] = VARIANTS["control"]
# Rodanthe village dune upkeep (2026-11-03), on control_rd198, at GIS 77-83 only. The model over-erodes the village
# shoreline 2009-2018 (about 4-5 m/yr against CoastSat's ~1) while its dunes there (GIS 78, 82) drop 0.7-1 m in
# 2010-11 and overwash every year; the lidar crests held within a few tenths of a metre 2009-2014
# (1-observations/7-rodanthe-dunes-2009-2014). Two ways the village could keep its dunes:
#   option 1, OVERWASH_TO_DUNE: more of each year's overwash pushed back onto the dune (the village default bulldozes
#     9% and returns 40% to the shoreface, HATTERAS_BEACH_DUNE); the sand the storm moved, nothing added
#   option 2, DUNE_UPKEEP: a dune cell more than the trigger below its 2009 starting height is rebuilt to that height
#     after the year's storms, the volume taken off the beach (the shoreline moves landward by 2V / (2 h_b + D_sf),
#     shoreface_nourishment run backwards)
UPKEEP_GIS = tuple(range(77, 84))
OVERWASH_TO_DUNE = {"control_rd198_ow30": 30.0, "control_rd198_ow59": 59.0}   # percent; filter 40 + this < 100
DUNE_UPKEEP = {"control_rd198_upkeep": 0.5}                                     # trigger, m below the 2009 crest
for _v in (*OVERWASH_TO_DUNE, *DUNE_UPKEEP):
    VARIANTS[_v] = VARIANTS["control"]
# The village overwash filter (2026-11-03): the share of each year's overwash every village (Buxton, Avon, the
# Tri-Village) returns to the shoreface, one setting for all, at points in Rogers et al.'s (2015) 40-90%; 40 is the
# hindcast's, and the 9% onto the dune is kept. Chosen on 1996-2009 (HAT_dune_collision.py VILLAGE_FILTER), checked
# here. _vf40 reruns each base on the toes from the survey windows that follow the dune (regenerated the same day)
VILLAGE_FILTER = {f"{b}_vf{p}": p for b in ("control_rd198", "nocape_rd198_dry") for p in (40, 65, 90)}
for _v in VILLAGE_FILTER:
    VARIANTS[_v] = VARIANTS[_v.rsplit("_vf", 1)[0]]
DRY_BEACH = DRY_BEACH + tuple(v for v in VILLAGE_FILTER if v.startswith("nocape_rd198_dry"))
# nocape_rd198_dry on toes and beach widths from the 2009 start profile alone (TOES_FROM_START), no survey after
# the start (2026-11-04): the clean hindcast against the default median of the 2014-2019 surveys
VARIANTS["nocape_rd198_dry_ts"] = VARIANTS["nocape"]
TOES_FROM_START = TOES_FROM_START + ("nocape_rd198_dry_ts",)
DRY_BEACH = DRY_BEACH + ("nocape_rd198_dry_ts",)
# The model-change steps (2026-11-04), each adding one change to the one before: the base; + the module (dry beach,
# W0 25 m, Cs 1.02e-4, off at Cape Point); + the 198 m road rule; then nocape_rd198_dry (+ Rodanthe 90%, the
# default since 2026-11-04) and nocape_rd198_dry_ts (+ 2009 start toes). Steps 1-3 hold Rodanthe at 40% (RODANTHE_40)
VARIANTS["step1_base"] = VARIANTS["control"]
VARIANTS["step2_module"] = VARIANTS["step3_road"] = VARIANTS["nocape"]
DRY_BEACH = DRY_BEACH + ("step2_module", "step3_road")
RODANTHE_40 = ("step1_base", "step2_module", "step3_road")
# nocape_rd198_dry with one toe per domain, the mean of its screened cells over the lidar up to the start, 2004, 2005
# and 2009 (read on the 1996 start's windows moved with the shoreline; nothing after 2009), 2026-11-04
VARIANTS["nocape_rd198_dry_toeavg"] = VARIANTS["nocape"]
TOE_DOMAIN_AVG = {"nocape_rd198_dry_toeavg": (2004, 2005, 2009)}
TOE_DOMAIN_AVG_TABLE = PROJECT_ROOT / "output" / "raw_runs" / "dune-storm-response" / "2026-10-26-dune-collision" / "dune_toes.csv"
DRY_BEACH = DRY_BEACH + tuple(TOE_DOMAIN_AVG)
# Steps 2 and 3 again on the domain toes, which dune_toe_cells.csv holds since 2026-11-04 (HAT_dune_toes.py
# TOE_DOMAIN_SURVEYS; step 4 is nocape_rd198_dry_toeavg, the same toes)
VARIANTS["step2_module_dt"] = VARIANTS["step3_road_dt"] = VARIANTS["nocape"]
TOE_DOMAIN_AVG.update({"step2_module_dt": (2004, 2005, 2009), "step3_road_dt": (2004, 2005, 2009)})
# Domain toes with the per-cell beach widths (the 2011-2020 survey medians, beach_width_cell_median_m), for the
# comparison with the domain beach factor
VARIANTS["nocape_rd198_dry_toeavg_cw"] = VARIANTS["nocape"]
TOE_DOMAIN_AVG["nocape_rd198_dry_toeavg_cw"] = (2004, 2005, 2009)
CELL_WIDTHS = ("nocape_rd198_dry_toeavg_cw",)
DRY_BEACH = DRY_BEACH + CELL_WIDTHS
DRY_BEACH = DRY_BEACH + ("step2_module_dt", "step3_road_dt")
# nocape_rd198_dry_toeavg with Larson alone: the beach factor 1 everywhere (width 0), Cs 1.9e-5 its offline fit
# without the factor (rule comparison, larson_storm, re-read crests; 2.0e-5 on the old), 2026-11-04
LARSON_ALONE = {"nocape_rd198_dry_toeavg_la": 1.9e-5}
VARIANTS["nocape_rd198_dry_toeavg_la"] = VARIANTS["nocape"]
TOE_DOMAIN_AVG["nocape_rd198_dry_toeavg_la"] = (2004, 2005, 2009)
DRY_BEACH = DRY_BEACH + tuple(LARSON_ALONE)
RODANTHE_40 = RODANTHE_40 + ("step2_module_dt", "step3_road_dt")
SHORELINE_GATE = ("nocape_rd198", "nocape_rd198_rl02", "nocape_rd198_cs145", "nocape_rd198_dry", "control_rd198", "nocape_rd198_dry_sand", "nocape_rd198_dry_ts", "step3_road", "nocape_rd198_dry_toeavg", "step3_road_dt", "nocape_rd198_dry_toeavg_cw", "nocape_rd198_dry_toeavg_la",
                  *OVERWASH_TO_DUNE, *DUNE_UPKEEP, *VILLAGE_FILTER)
SHORELINE_GATE_M = 198.0
# Larson alone's Cs, in the model (2026-11-04): the calibration's sweep (HAT_dune_collision.py LA_CS_SWEEP) run
# here to show the sensitivity; Cs is chosen on 1996-2009
for _cs in (1e-5, 3e-5, 5e-5, 8e-5):
    _v = f"la_cs{round(_cs * 1e6):d}"
    LARSON_ALONE[_v] = _cs
    VARIANTS[_v] = VARIANTS["nocape"]
    TOE_DOMAIN_AVG[_v] = (2004, 2005, 2009)
    DRY_BEACH = DRY_BEACH + (_v,)
    SHORELINE_GATE = SHORELINE_GATE + (_v,)
# The model-change steps with Larson alone, the default since 2026-11-04 (Cs 1.9e-5, no beach factor): steps 2 and
# 3 again (Rodanthe 40%); step 4 is nocape_rd198_dry_toeavg_la
for _v in ("step2_module_la", "step3_road_la"):
    LARSON_ALONE[_v] = 1.9e-5
    VARIANTS[_v] = VARIANTS["nocape"]
    TOE_DOMAIN_AVG[_v] = (2004, 2005, 2009)
    DRY_BEACH = DRY_BEACH + (_v,)
    RODANTHE_40 = RODANTHE_40 + (_v,)
SHORELINE_GATE = SHORELINE_GATE + ("step3_road_la",)
# nocape_rd198_dry_toeavg_la with the 198 m gate as one island-wide road trigger, 198 m less the island's median beach
# (~50 m) = 150 m, the roadway manager's own setting (road_dune_rebuild_trigger_m), in place of the per-domain
# gates (2026-11-04)
ROAD_TRIGGER_M = {"final_rd150u": 150.0}
LARSON_ALONE["final_rd150u"] = 1.9e-5
VARIANTS["final_rd150u"] = VARIANTS["nocape"]
TOE_DOMAIN_AVG["final_rd150u"] = (2004, 2005, 2009)
DRY_BEACH = DRY_BEACH + ("final_rd150u",)
# The beach lowering the water at the dune, m per m of beach (DuneCollisionManager beach_runup_loss)
RUNUP_LOSS = {"nocape_rd198_rl02": 0.02}
WIDTH_FROM_START = ("nocape_dynw_sw",)
# Cs per variant, the hindcast's 1.24e-4 otherwise: 2.5e-4 the best island crest fit of the 1996-2009 sweep on
# start-profile toes (HAT_dune_collision.py cal_g16_nocape_ts250)
CS = {"nocape_toestart_cs250": 2.5e-4, "nocape_rd198_rl02": 1.64e-4, "nocape_rd198_cs145": 1.45e-4}
TOE_SURVEYS = (PROJECT_ROOT / "data" / "hatteras_init" / "1-barrier3d-domains" / "2009-start" / "dune-toes"
               / "survey_dune_toes.csv")
TOE_MIN_M_NAVD88, TOE_BELOW_CREST_M, MHW_M_NAVD88 = 1.0, 0.5, 0.36   # as HAT_dune_toes.py and the hindcast setup
# (name, first GIS, last GIS), as the calibration's zones
ZONES = (("cape", 1, 7), ("buxton_avon_road", 8, 20), ("avon", 21, 31), ("avon_salvo_road", 32, 67),
         ("tri_village", 68, 83), ("rodanthe_north", 84, 90))
WORKERS = 3
# The surveys the dune crests are scored on: the runs before 2026-11-03 compared these five; later runs also
# carry 2011, 2012 and 2020 (added to the crest surveys that day), left out so every run is scored alike
SCORE_SURVEYS = (2014, 2016, 2017, 2018, 2019)
# -----------------------------------------------------------------------------


# The newest run folder of a variant
def run_dir(variant):
    base = EXP_DIR / "runs" / variant
    meta = sorted(base.glob("**/HAT_*/*_run_metadata.json"), key=lambda f: f.stat().st_mtime) if base.exists() else []
    return meta[-1].parent if meta else None


# Each road manager's rebuild gate (m behind the dune) set to SHORELINE_GATE_M less its domain's median beach width
def shoreline_gate(c, toe_cells):
    from site_layer.hatteras_site_config import HATTERAS_DOMAINS
    dom = HATTERAS_DOMAINS.for_total(len(c.barrier3d))
    t = pd.read_csv(toe_cells)
    col = "beach_to_crest_m" if "beach_to_crest_m" in t else "beach_width_m"   # the road is measured from the crest
    w = t.groupby("gis")[col].median().reindex(range(1, 91))
    w = w.fillna(w.median())
    for g in range(1, 91):
        r = c.roadways[dom.gis_to_pad(g)]
        if r is not None:
            r._dune_mgmt_road_setback_trigger = max(SHORELINE_GATE_M - w[g], 0.0)
    print(f"EXPERIMENT road gates {SHORELINE_GATE_M:g} m less the beach: median "
          f"{max(SHORELINE_GATE_M - w.median(), 0):.0f} m, range {SHORELINE_GATE_M - w.max():.0f}-"
          f"{SHORELINE_GATE_M - w.min():.0f} m", flush=True)


# Child process: the hindcast, with the module left out of the variant's GIS
def _launch(variant):
    import runpy
    import cascade_pipeline.hindcast as hc
    from site_layer.hatteras_site_config import HATTERAS_DOMAINS
    _, off = VARIANTS[variant]
    orig_build = hc.build_cascade

    def build_cascade(*a, **k):
        if k.get("dune_collision_module") and k.get("dune_collision_beach_width") is None:
            # hat_run.yaml's default since 2026-11-04 is Larson alone (beach: none), which sends no widths; these
            # variants were made with the beach factor, so give them the toe table's widths as before (Larson-alone
            # variants zero them below)
            dom = HATTERAS_DOMAINS.for_total(k["alongshore_section_count"])
            cw = pd.read_csv(TOE_SURVEYS.with_name("dune_toe_cells.csv")).sort_values(["gis", "cell"])
            widths = [None] * dom.total_domains
            for gis, d in cw.groupby("gis"):
                widths[dom.gis_to_pad(int(gis))] = d.beach_width_m.to_numpy()
            k["dune_collision_beach_width"] = widths
        if variant in TOE_DOMAIN_AVG and k.get("dune_collision_module"):
            t = pd.read_csv(TOE_DOMAIN_AVG_TABLE)
            t = t[(t.survey.astype(str).isin([str(y) for y in TOE_DOMAIN_AVG[variant]])) & (t.status == "ok")]
            t = t[(t.toe_m_navd88 > TOE_MIN_M_NAVD88) & (t.toe_m_navd88 < t.crest_m_navd88 - TOE_BELOW_CREST_M)]
            dm = t.groupby("gis").toe_m_navd88.mean()
            dom = HATTERAS_DOMAINS.for_total(k["alongshore_section_count"])
            toes = list(k["dune_collision_toe_elevation"])
            for gis in range(1, 91):
                i = dom.gis_to_pad(gis)
                if toes[i] is not None:
                    toes[i] = np.full(len(toes[i]), dm.get(gis, dm.median()) - MHW_M_NAVD88)
            k["dune_collision_toe_elevation"] = toes
            # one beach factor per domain, the mean of max(0, 1 - W / DRY_LIMIT_M) over every ok cell of the same
            # surveys, given to the module as the width DRY_LIMIT_M (1 - factor) (2026-11-04; HAT_dune_toes.py)
            w = pd.read_csv(TOE_DOMAIN_AVG_TABLE)
            w = w[(w.survey.astype(str).isin([str(y) for y in TOE_DOMAIN_AVG[variant]])) & (w.status == "ok")]
            dw = pd.Series(np.where(w.toe_to_crest_m.notna(), w.beach_width_m - w.toe_to_crest_m,
                                    np.where(w.beach_width_m < 2 * 10, 0.0, np.nan)), index=w.index)
            bf = np.clip(1 - dw / DRY_LIMIT_M, 0, None).groupby(w.gis).mean()
            wm = DRY_LIMIT_M * (1 - bf)
            widths = list(k["dune_collision_beach_width"])
            for gis in range(1, 91):
                i = dom.gis_to_pad(gis)
                if widths[i] is not None:
                    widths[i] = np.full(len(widths[i]), wm.get(gis, wm.median()))
            k["dune_collision_beach_width"] = widths
            print(f"EXPERIMENT one beach factor per domain, the mean over the same surveys: island {bf.mean():.2f}",
                  flush=True)
            if variant in CELL_WIDTHS:
                cw = pd.read_csv(TOE_SURVEYS.with_name("dune_toe_cells.csv")).sort_values(["gis", "cell"])
                widths = list(k["dune_collision_beach_width"])
                for gis, d in cw.groupby("gis"):
                    i = dom.gis_to_pad(int(gis))
                    if widths[i] is not None:
                        widths[i] = d.beach_width_cell_median_m.to_numpy()
                k["dune_collision_beach_width"] = widths
                print("EXPERIMENT ... replaced by the per-cell beach widths (beach_width_cell_median_m)", flush=True)
            print()
            print(f"EXPERIMENT one toe per domain, the mean over the {TOE_DOMAIN_AVG[variant]} surveys: island "
                  f"{dm.mean():.2f} m NAVD88, range {dm.min():.2f}-{dm.max():.2f}", flush=True)
        if variant in TOES_FROM_START and k.get("dune_collision_module"):
            t = pd.read_csv(TOE_SURVEYS)
            t = t[t.survey.astype(str) == "start"]
            good = t[(t.toe_m_navd88 > TOE_MIN_M_NAVD88) & (t.toe_m_navd88 < t.crest_m_navd88 - TOE_BELOW_CREST_M)]
            dom = HATTERAS_DOMAINS.for_total(k["alongshore_section_count"])
            toes, widths = list(k["dune_collision_toe_elevation"]), list(k["dune_collision_beach_width"])
            for gis, d in t.groupby("gis"):
                cells = d.sort_values("cell")
                toe = good[good.gis == gis].set_index("cell").toe_m_navd88.reindex(cells.cell)
                toe = toe.fillna(toe.median())
                if toe.notna().any():
                    toes[dom.gis_to_pad(int(gis))] = (toe - MHW_M_NAVD88).to_numpy()
                    widths[dom.gis_to_pad(int(gis))] = cells.beach_width_m.to_numpy()
            k["dune_collision_toe_elevation"], k["dune_collision_beach_width"] = toes, widths
            print("\nEXPERIMENT toes and beach widths from the 2009 start profile", flush=True)
        if variant in WIDTH_FROM_START and k.get("dune_collision_module"):
            t = pd.read_csv(TOE_SURVEYS)
            t = t[t.survey.astype(str) == "start"]
            dom = HATTERAS_DOMAINS.for_total(k["alongshore_section_count"])
            widths = list(k["dune_collision_beach_width"])
            for gis, d in t.groupby("gis"):
                if widths[dom.gis_to_pad(int(gis))] is not None:
                    widths[dom.gis_to_pad(int(gis))] = d.sort_values("cell").beach_width_m.to_numpy()
            k["dune_collision_beach_width"] = widths
            print("\nEXPERIMENT beach widths from the 2009 start profile", flush=True)
        if variant in BEACH_FOLLOWS_SHORELINE or variant in RUNUP_LOSS:
            raise SystemExit(f"{variant}: beach_follows_shoreline and beach_runup_loss were removed from "
                             "cascade/dune_collision.py on 2026-11-04; check out an earlier commit to rerun it")
        if variant in ROAD_DUNE_REBUILD:
            gate, trigger = ROAD_DUNE_REBUILD[variant]
            k["road_dune_rebuild_trigger_m"] = gate
            if trigger is not None:
                k["dune_minimum_elevation"] = [trigger] * len(k["dune_minimum_elevation"])
            print(f"\nEXPERIMENT road dune rebuild: setback gate {gate:g} m, trigger "
                  + ("berm floor" if trigger is None else f"{trigger:.2f} m MHW"), flush=True)
        if variant in SHORELINE_GATE:
            print(f"\nEXPERIMENT road dune rebuild where NC 12 is within {SHORELINE_GATE_M:g} m of the shoreline", flush=True)
        if variant in DRY_BEACH and k.get("dune_collision_module"):
            k["dune_collision_cs"], k["dune_collision_beach_width_limit"] = DRY_CS, DRY_LIMIT_M
            print(f"\nEXPERIMENT dry beach (start to toe), Cs {DRY_CS:.2e}, limit {DRY_LIMIT_M:g} m", flush=True)
        if variant in CS and k.get("dune_collision_module"):
            k["dune_collision_cs"] = CS[variant]
            print(f"\nEXPERIMENT Cs {CS[variant]:.2e}", flush=True)
        if variant in OVERWASH_TO_DUNE:
            dom = HATTERAS_DOMAINS.for_total(k["alongshore_section_count"])
            share = k["overwash_to_dune"]
            share = list(share) if np.size(share) > 1 else [share] * k["alongshore_section_count"]
            for g in UPKEEP_GIS:
                share[dom.gis_to_pad(g)] = OVERWASH_TO_DUNE[variant]
            k["overwash_to_dune"] = share
            print(f"\nEXPERIMENT {OVERWASH_TO_DUNE[variant]:g}% of the overwash bulldozed onto the dune at GIS "
                  f"{UPKEEP_GIS[0]}-{UPKEEP_GIS[-1]}", flush=True)
        if variant in ROAD_TRIGGER_M:
            k["road_dune_rebuild_trigger_m"] = ROAD_TRIGGER_M[variant]
            print(f"\nEXPERIMENT road dune rebuild trigger {ROAD_TRIGGER_M[variant]:g} m island-wide", flush=True)
        if variant in LARSON_ALONE and k.get("dune_collision_module"):
            k["dune_collision_cs"] = LARSON_ALONE[variant]
            k["dune_collision_beach_width"] = [None if w is None else np.zeros(len(w))
                                               for w in k["dune_collision_beach_width"]]
            print()
            print(f"EXPERIMENT Larson alone: beach factor 1 everywhere, Cs {LARSON_ALONE[variant]:.1e}", flush=True)
        if off and k.get("dune_collision_module"):
            dom = HATTERAS_DOMAINS.for_total(k["alongshore_section_count"])
            for key in ("dune_collision_toe_elevation", "dune_collision_beach_width"):
                if k.get(key) is not None:
                    k[key] = list(k[key])
                    for g in off:
                        k[key][dom.gis_to_pad(g)] = None
            print(f"\nEXPERIMENT collision off at GIS {off[0]}-{off[-1]}", flush=True)
        # The module has no beach since 2026-11-04: Larson-alone variants (widths zeroed above) run as before
        widths = k.pop("dune_collision_beach_width", None)
        k.pop("dune_collision_beach_width_limit", None)
        if (k.get("dune_collision_module") and variant not in LARSON_ALONE
                and widths is not None and any(w is not None for w in widths)):
            raise SystemExit(f"{variant}: the module's beach options were removed from cascade/dune_collision.py on "
                             "2026-11-04; check out git tag collision-with-beach-factor to rerun it")
        c = orig_build(*a, **k)
        if variant in SHORELINE_GATE:
            shoreline_gate(c, TOE_SURVEYS.with_name("dune_toe_cells.csv"))
        if variant in RODANTHE_40:
            dom = HATTERAS_DOMAINS.for_total(k["alongshore_section_count"])
            n = 0
            for g in range(77, 84):
                i = dom.gis_to_pad(g)
                if c.nourishments[i] is not None and c._beach_nourishment_module[i]:
                    c.nourishments[i]._overwash_filter = 40.0
                    n += 1
            print(f"\nEXPERIMENT Rodanthe overwash filter held at 40% on {n} domains (GIS 77-83)", flush=True)
        if variant in VILLAGE_FILTER:
            n = 0
            for i, mgr in enumerate(c.nourishments):
                if mgr is not None and c._beach_nourishment_module[i]:
                    mgr._overwash_filter = VILLAGE_FILTER[variant]
                    n += 1
            print(f"\nEXPERIMENT village overwash filter {VILLAGE_FILTER[variant]:g}% to the shoreface on {n} "
                  f"village domains", flush=True)
        if variant in DUNE_UPKEEP:
            dom = HATTERAS_DOMAINS.for_total(len(c.barrier3d))
            for g in UPKEEP_GIS:
                p = dom.gis_to_pad(g)
                # the 2009 starting dune (dam above the berm), the height each cell is kept to
                c.nourishments[p]._upkeep_target = np.asarray(c.barrier3d[p].DuneDomain[0], float).copy()
                c.nourishments[p]._upkeep_gis = g
        return c

    if variant in DUNE_UPKEEP:
        from cascade.beach_dune_manager import BeachDuneManager, shoreface_nourishment
        trigger_dam, orig_bdm_update = DUNE_UPKEEP[variant] / 10.0, BeachDuneManager.update
        upkeep_log = []

        def bdm_update(self, barrier3d, nourish_now, rebuild_dune_now, nourishment_interval):
            out = orig_bdm_update(self, barrier3d, nourish_now, rebuild_dune_now, nourishment_interval)
            target = getattr(self, "_upkeep_target", None)
            if target is None or self.narrow_break:
                return out
            t = barrier3d.time_index - 1
            dunes = np.asarray(barrier3d.DuneDomain[t], float)
            low = dunes.max(axis=1) < target.max(axis=1) - trigger_dam
            if not low.any():
                return out
            new = dunes.copy()
            new[low] = np.maximum(dunes[low], target[low])
            v = float((new - dunes).sum()) / barrier3d.BarrierLength   # dam^3 per dam of shore
            barrier3d.DuneDomain[t] = new
            # the sand comes off the beach: the nourishment formula with the volume taken away
            x_s, s_sf, _ = shoreface_nourishment(barrier3d.x_s, barrier3d.x_t, -v, barrier3d.h_b_TS[-1],
                                                 barrier3d.DShoreface, 0.0)
            dx_m = (x_s - barrier3d.x_s) * 10
            barrier3d.x_s, barrier3d.x_s_TS[-1], barrier3d.s_sf_TS[-1] = x_s, x_s, s_sf
            self._beach_width[t] = max(self._beach_width[t] - dx_m, 0.0)
            upkeep_log.append(dict(year=RUN_START + t, gis=self._upkeep_gis, cells=int(low.sum()),
                                   volume_m3_per_m=v * 100, shoreline_landward_m=dx_m))
            return out

        BeachDuneManager.update = bdm_update
        import atexit

        def write_log():
            if upkeep_log:
                rd = run_dir(variant)
                if rd is not None:
                    (rd / "tables").mkdir(exist_ok=True)
                    pd.DataFrame(upkeep_log).to_csv(rd / "tables" / "dune_upkeep_log.csv", index=False)
        atexit.register(write_log)
        print(f"\nEXPERIMENT village dunes at GIS {UPKEEP_GIS[0]}-{UPKEEP_GIS[-1]} kept to their 2009 height "
              f"(rebuilt when {DUNE_UPKEEP[variant]:g} m below), the sand taken off the beach", flush=True)

    if variant in SAND_TO_BEACH:
        from cascade import dune_collision as dc
        from cascade.beach_dune_manager import shoreface_nourishment
        share, orig_update = SAND_TO_BEACH[variant], dc.DuneCollisionManager.update

        def update(self, barrier3d, DuneDomain, t):
            before = self._volume_TS[-1] if self._volume_TS else 0.0
            orig_update(self, barrier3d, DuneDomain, t)
            v = (self._volume_TS[-1] - before) * share   # dam^3 per dam of shore, this storm
            if v > 0:
                barrier3d._x_s = shoreface_nourishment(barrier3d._x_s, barrier3d._x_t, v, barrier3d._h_b_TS[-1],
                                                       barrier3d._DShoreface, 0.0)[0]
                self.sand_to_beach_dam3 = getattr(self, "sand_to_beach_dam3", 0.0) + v

        dc.DuneCollisionManager.update = update
        print(f"\nEXPERIMENT {share:.0%} of the collision-eroded dune sand onto the beach and shoreface", flush=True)
    hc.build_cascade = build_cascade
    sys.path.insert(0, str(HINDCAST.parent))
    sys.argv = [str(HINDCAST)]
    if variant in SAND_TO_BEACH:
        # Cascade's domains run in joblib worker processes, which import the module afresh and would not see the
        # patched update: threads keep them in this process
        import joblib
        with joblib.parallel_config(backend="threading"):
            runpy.run_path(str(HINDCAST), run_name="__main__")
    else:
        runpy.run_path(str(HINDCAST), run_name="__main__")


# One variant in its own process, logged
def launch(variant):
    logs = EXP_DIR / "logs"
    logs.mkdir(parents=True, exist_ok=True)
    on, _ = VARIANTS[variant]
    env = dict(os.environ, HAT_IGNORE_SETTINGS="1", HAT_START_YEAR=str(RUN_START), HAT_SOURCE_SINK_PRESET=PRESET,
               HAT_SCENARIO=SCENARIO_OF.get(variant, SCENARIO), HAT_GROIN_ENABLED="1", HAT_DUNE_COLLISION="1" if on else "0",
               HAT_RUN_KIND="dune-storm-response", HAT_RUN_TAG=f"{TAG}/runs/{variant}", HAT_SAVE_MODEL_STATE="1",
               HAT_MAKE_GIFS="0", HAT_SHOW_FIGURES="0", HAT_OVERWRITE="1", PYTHONIOENCODING="utf-8")
    if variant in HS_OF:
        hs = HS_OF[variant]
        env["HAT_HS"] = repr(hs)
        env["HAT_BE_OVERRIDE"] = f"1={ENDS_BY_HS[hs][0]},90={ENDS_BY_HS[hs][1]}"
    if variant in GROWTH:
        lo, hi = GROWTH[variant]
        env["HAT_GROWTH_OVERRIDE"] = ",".join(f"{g}={lo:.3f}:{hi:.3f}" for g in range(1, 91))
    t0 = time.time()
    with open(logs / f"{variant}.log", "w", encoding="utf-8") as fh:
        p = subprocess.run([sys.executable, str(_HERE), "_launch", variant], env=env, cwd=HINDCAST.parent,
                           stdout=fh, stderr=subprocess.STDOUT)
    print(f"{variant:10s} exit {p.returncode}  {(time.time() - t0) / 60:.1f} min", flush=True)


def score():
    from cascade_pipeline.coastsat_lowess import build_coastsat_series
    from cascade_pipeline.hindcast import build_target_table
    from cascade_pipeline.plotting.dune_crest_comparison import TABLE, TABLES
    from site_layer.hatteras_site_config import HATTERAS_DOMAINS as DOM
    sys.path.insert(0, str(_HAT_DIR / "figure_making" / "model_output"))
    from rerender_run_figures import COASTSAT_DATASETS, LOWESS_CONFIG
    from run_figures import TARGET_WINDOW
    cs = build_coastsat_series(COASTSAT_DATASETS, active_period_start=RUN_START, lowess_config=LOWESS_CONFIG,
                               domains=DOM)
    target = build_target_table(next(c for c in cs if c["active"]), LOWESS_CONFIG, DOM, TARGET_WINDOW)
    target = target.set_index("gis_domain")["target_lrr_m_yr"]
    rows = []
    for v in VARIANTS:
        rd = run_dir(v)
        if rd is None:
            print(f"missing {v}")
            continue
        rate = pd.read_csv(rd / "tables" / "shoreline_change_rate.csv").set_index("gis_domain").lrr_m_yr
        meta = json.loads((rd / f"{rd.name}_run_metadata.json").read_text(encoding="utf-8"))
        err = (rate - target).loc[2:89]
        row = dict(run=v, rmse=float(meta["skill"]["rmse_interior_m_yr"]),
                   bias=float(meta["skill"]["mean_bias_interior_m_yr"]))
        for z, lo, hi in ZONES:
            e = err.loc[max(lo, 2):min(hi, 89)]
            row[f"lrr_rmse_{z}"] = float(np.sqrt((e ** 2).mean()))
        crest_table = rd / TABLES / TABLE
        if crest_table.exists():
            d = pd.read_csv(crest_table)
            d = d[d.survey.isin(SCORE_SURVEYS)].model_minus_survey_m.dropna()
            row.update(dune_cells=len(d), dune_within_05=float((d.abs() <= 0.5).mean()),
                       dune_over_1m_high=float((d > 1).mean()), dune_over_1m_low=float((d < -1).mean()))
        rows.append(row)
    t = pd.DataFrame(rows).set_index("run")
    (EXP_DIR / "tables").mkdir(parents=True, exist_ok=True)
    t.round(4).to_csv(EXP_DIR / "tables" / "scores.csv")
    pd.set_option("display.width", 200)
    print(t.round(3).T.to_string())


def main():
    ap = argparse.ArgumentParser(description=__doc__.strip().splitlines()[0])
    ap.add_argument("action", choices=["run", "score", "_launch"])
    ap.add_argument("variant", nargs="*")
    a = ap.parse_args()
    if a.action == "_launch":
        _launch(a.variant[0])
    elif a.action == "run":
        with ThreadPoolExecutor(WORKERS) as ex:
            list(ex.map(launch, a.variant or list(VARIANTS)))
    else:
        score()


if __name__ == "__main__":
    main()
