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
# NCDOT's own criterion in place of the fitted gate: dune maintenance where NC 12 is within 650 ft (198 m) of the
# shoreline (Bonner Bridge Replacement FEIS, NCDOT/FHWA 2008, p. 2-128, from Isabel's washover fans). The road
# manager measures the road from the dune, so each domain's gate is 198 m less its median beach width
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
SHORELINE_GATE = ("nocape_rd198", "nocape_rd198_rl02", "nocape_rd198_cs145", "nocape_rd198_dry")
SHORELINE_GATE_M = 198.0
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
        if off and k.get("dune_collision_module"):
            dom = HATTERAS_DOMAINS.for_total(k["alongshore_section_count"])
            for key in ("dune_collision_toe_elevation", "dune_collision_beach_width"):
                if k.get(key) is not None:
                    k[key] = list(k[key])
                    for g in off:
                        k[key][dom.gis_to_pad(g)] = None
            print(f"\nEXPERIMENT collision off at GIS {off[0]}-{off[-1]}", flush=True)
        c = orig_build(*a, **k)
        if variant in SHORELINE_GATE:
            shoreline_gate(c, TOE_SURVEYS.with_name("dune_toe_cells.csv"))
        return c

    hc.build_cascade = build_cascade
    sys.path.insert(0, str(HINDCAST.parent))
    sys.argv = [str(HINDCAST)]
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
