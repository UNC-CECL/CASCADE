"""
The 1996-2010 window on today's setup: run through December 2010, both ends re-solved.

    python scripts/hatteras_ms/experiments/HAT_window_1996_2010.py

Drives the unchanged runner with the 1996 period patched in-process to end in 2010
(last model year 2010, the 1996_2010 storms, graded on CoastSat 1996_2010): zeroBE once,
then edgeBE with GIS 1 and GIS 90 solved together by secant steps. No fill falls in the
window. Figures: HAT_candidate_windows_plot.py --matrix --calib 1996_2010 (into this study's figures/).

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-10-04
"""
from __future__ import annotations

import json
import os
import subprocess
import sys
from pathlib import Path

import pandas as pd

_HERE = Path(__file__).resolve()
PROJECT_ROOT = next(_p for _p in _HERE.parents if (_p / "pyproject.toml").exists())
sys.path.insert(0, str(PROJECT_ROOT / "scripts"))
sys.path.insert(0, str(PROJECT_ROOT))

# --- CONFIG ------------------------------------------------------------------
HINDCAST = PROJECT_ROOT / "scripts" / "hatteras_ms" / "HAT_hindcast_1984_2024.py"
TAG = "candidate-windows/2026-10-04-window-1996-2010"
EXP_DIR = PROJECT_ROOT / "output" / "raw_runs" / "experiments" / TAG
RUN_DIR_TAG = f"{TAG}/runs"
START, END = 1996, 2010
START_ENDS = {1: 4.3888, 90: 19.0935}    # the 2026-09-29 1996-2010 solve, the first probe
FIRST_GAIN = {1: 0.70, 90: 0.15}         # d(LRR)/d(BE) at each end, from the 2026-10-04 solves
TOL = 0.05
MAX_STEPS = 6
SCENARIO = "full_management"
# -----------------------------------------------------------------------------


# Subprocess entry: end the 1996 period in 2010, then run the unchanged runner
def _launch():
    import runpy
    from site_layer import hatteras_site_config as sc
    from site_layer import hat_env_forcings as env
    p = sc.HATTERAS_PERIODS[START]
    p["end_year"] = END
    p["last_model_year"] = END            # CoastSat 1996_2010 and its storm file run through Dec 2010
    p["storm_file"] = env.init_relpath(env.storm_series_file(START, END))
    print(f"1996-2010 window: storms {p['storm_file']}, RSLR {p['sea_level_rise_rate']}")
    sys.argv = [str(HINDCAST)]
    runpy.run_path(str(HINDCAST), run_name="__main__")


def run(preset, label, override=None):
    tag = f"{RUN_DIR_TAG}/{preset}_{label}"
    e = dict(os.environ, PYTHONIOENCODING="utf-8", PYTHONUNBUFFERED="1", HAT_IGNORE_SETTINGS="1",
             HAT_START_YEAR=str(START), HAT_SOURCE_SINK_PRESET=preset, HAT_SCENARIO=SCENARIO,
             HAT_RELOCATIONS="False", HAT_GROIN_ENABLED="False", HAT_RUN_KIND="experiment",
             HAT_RUN_TAG=tag, HAT_SAVE_MODEL_STATE="False", HAT_OVERWRITE="False")
    if override:
        e["HAT_BE_OVERRIDE"] = override
    logs = EXP_DIR / "logs"
    logs.mkdir(parents=True, exist_ok=True)
    with open(logs / f"{preset}_{label}.log", "w", encoding="utf-8") as log:
        code = subprocess.run([sys.executable, str(_HERE), "_launch"], stdout=log, stderr=subprocess.STDOUT,
                              env=e, cwd=PROJECT_ROOT).returncode
    if code:
        raise SystemExit(f"run {tag} failed ({code}); see its log")
    base = PROJECT_ROOT / "output" / "raw_runs" / "experiments" / RUN_DIR_TAG / f"{preset}_{label}"
    return next(base.rglob("tables/shoreline_change_rate.csv")).parent.parent


# The targets at both ends, as the runner scores them (LOWESS-7, raw at GIS 1-10)
def targets():
    from cascade_pipeline.hindcast import build_target_table
    from cascade_pipeline.coastsat_lowess import CoastSatDataset, LowessConfig, build_coastsat_series
    from site_layer.hatteras_site_config import HATTERAS_DOMAINS
    from site_layer.hat_observed_rates import lrr_csv
    cfg = LowessConfig(window_domains=(7,), skip_southern_domains=10)
    s = build_coastsat_series([CoastSatDataset(label="c", period_start=START, csv_path=str(lrr_csv(START, END)))],
                              START, cfg, domains=HATTERAS_DOMAINS)
    t = build_target_table(s[0], cfg, HATTERAS_DOMAINS, 7).set_index("gis_domain")["target_lrr_m_yr"]
    return {1: float(t.loc[1]), 90: float(t.loc[90])}


def main():
    if len(sys.argv) > 1 and sys.argv[1] == "_launch":
        _launch()
        return
    run("zeroBE", "zero")
    tgt = targets()
    be, hist = dict(START_ENDS), []
    for k in range(MAX_STEPS):
        d = run("edgeBE", f"step{k}", override=f"1={be[1]:.3f},90={be[90]:.3f}")
        lrr = pd.read_csv(d / "tables" / "shoreline_change_rate.csv").set_index("gis_domain")["lrr_m_yr"]
        res = {g: float(lrr.loc[g]) - tgt[g] for g in (1, 90)}
        hist.append(dict(step=k, gis1=be[1], gis90=be[90], res1=res[1], res90=res[90]))
        print(f"step{k}  GIS 1 {be[1]:+.3f} -> {res[1]:+.3f}   GIS 90 {be[90]:+.3f} -> {res[90]:+.3f}", flush=True)
        if all(abs(res[g]) < TOL for g in (1, 90)):
            break
        for g, rk, bk in ((1, "res1", "gis1"), (90, "res90", "gis90")):
            if abs(res[g]) < TOL:
                continue
            gain = FIRST_GAIN[g]
            if len(hist) >= 2 and hist[-1][bk] != hist[-2][bk]:
                sec = (hist[-1][rk] - hist[-2][rk]) / (hist[-1][bk] - hist[-2][bk])
                gain = sec if sec > 0.02 else gain
            be[g] = be[g] - res[g] / gain
    best = min(hist, key=lambda h: abs(h["res1"]) + abs(h["res90"]))
    out = EXP_DIR / "runs"
    out.mkdir(parents=True, exist_ok=True)
    (out / "solve_summary.json").write_text(json.dumps(dict(history=hist, best=best), indent=1), encoding="utf-8")
    print("best", best)


if __name__ == "__main__":
    main()
