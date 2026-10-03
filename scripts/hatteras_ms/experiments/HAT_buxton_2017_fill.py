"""
Does adding the 2017 Buxton fill improve the 2010-2026 hindcast?

    python scripts/hatteras_ms/experiments/HAT_buxton_2017_fill.py run
    python scripts/hatteras_ms/experiments/HAT_buxton_2017_fill.py score

Runs 2010-2026 edgeBE full_management twice, with the fill list as committed
and with Buxton 2017 removed in-process; scored against the CoastSat LRR target
island-wide and over the fill footprint, GIS 6-15. Details: scripts/hatteras_ms/experiments/README.md.

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-10-03
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
PROJECT_ROOT = next(_p for _p in _HERE.parents if (_p / "pyproject.toml").exists())
sys.path.insert(0, str(PROJECT_ROOT / "scripts"))

# --- CONFIG ------------------------------------------------------------------
TAG = "management/2026-10-03-buxton-2017-fill"
EXP_DIR = PROJECT_ROOT / "output" / "raw_runs" / "experiments" / TAG
HINDCAST = PROJECT_ROOT / "scripts" / "hatteras_ms" / "HAT_hindcast_1984_2024.py"
MATRIX_RUN = (PROJECT_ROOT / "output" / "raw_runs" / "matrix" / "2010_2026" / "edgeBE"
              / "HAT_2010_2026_edgeBE_offsetmetres_road_bdm_nourish_nogroin")
WINDOW = (2010, 2026)
VARIANTS = ("with2017", "without2017")
DROP = ("Buxton beach nourishment", 2017)
FOOTPRINT_GIS = tuple(range(6, 16))
TARGET_WINDOW = 7
# -----------------------------------------------------------------------------


# Child process: drop the fill if asked, then run the unchanged hindcast
def _launch(variant):
    import runpy
    sys.path.insert(0, str(HINDCAST.parent))
    from site_layer import hatteras_site_config as sc
    if variant == "without2017":
        kept = tuple(p for p in sc.HATTERAS_NOURISHMENT_PROJECTS if (p.name, p.year) != DROP)
        assert len(kept) == len(sc.HATTERAS_NOURISHMENT_PROJECTS) - 1
        sc.HATTERAS_NOURISHMENT_PROJECTS = kept
    sys.argv = [str(HINDCAST)]
    runpy.run_path(str(HINDCAST), run_name="__main__")


def run_dir(v):
    base = EXP_DIR / "runs" / v / f"{WINDOW[0]}_{WINDOW[1]}" / "edgeBE"
    hits = sorted(d for d in base.glob("HAT_*") if any(d.glob("*_run_metadata.json"))) if base.exists() else []
    return hits[-1] if hits else None


def launch(v):
    logs = EXP_DIR / "logs"
    logs.mkdir(parents=True, exist_ok=True)
    env = dict(os.environ)
    env.update({"HAT_IGNORE_SETTINGS": "1", "HAT_START_YEAR": str(WINDOW[0]),
                "HAT_SOURCE_SINK_PRESET": "edgeBE", "HAT_SCENARIO": "full_management",
                "HAT_RUN_KIND": "experiment", "HAT_RUN_TAG": f"{TAG}/runs/{v}",
                "HAT_SAVE_MODEL_STATE": "1", "HAT_MAKE_GIFS": "0", "HAT_SHOW_FIGURES": "0",
                "HAT_OVERWRITE": "1", "PYTHONIOENCODING": "utf-8"})
    t0 = time.time()
    log = logs / f"{v}.log"
    with open(log, "w", encoding="utf-8") as fh:
        p = subprocess.run([sys.executable, str(_HERE), "_launch", v], env=env, cwd=HINDCAST.parent,
                           stdout=fh, stderr=subprocess.STDOUT)
    rec = dict(variant=v, returncode=p.returncode, minutes=round((time.time() - t0) / 60, 1))
    with open(logs / "launches.jsonl", "a", encoding="utf-8") as fh:
        fh.write(json.dumps(rec) + "\n")
    print(f"{v:12s} exit {p.returncode} {rec['minutes']} min", flush=True)


# CoastSat target as the runner scores it, plus the raw per-domain mean
def targets():
    from cascade_pipeline.hindcast import build_target_table
    from cascade_pipeline.coastsat_lowess import CoastSatDataset, LowessConfig, build_coastsat_series
    from site_layer.hatteras_site_config import HATTERAS_DOMAINS
    from site_layer.hat_observed_rates import lrr_csv
    cfg = LowessConfig(window_domains=(TARGET_WINDOW,), skip_southern_domains=10)
    series = build_coastsat_series([CoastSatDataset(label=f"CoastSat {WINDOW[0]}", period_start=WINDOW[0],
                                                    csv_path=str(lrr_csv(*WINDOW)))], WINDOW[0], cfg,
                                   domains=HATTERAS_DOMAINS)
    t = build_target_table(series[0], cfg, HATTERAS_DOMAINS, TARGET_WINDOW).set_index("gis_domain")
    raw = pd.read_csv(lrr_csv(*WINDOW)).groupby("domain_number")["lrr_m_yr"].mean()
    t["raw_mean_lrr_m_yr"] = raw.reindex(t.index).values
    return t


def model_lrr(rd):
    return pd.read_csv(rd / "tables" / "shoreline_change_rate.csv").set_index("gis_domain")["lrr_m_yr"]


def stats(resid):
    r = resid.dropna()
    return float(r.mean()), float(np.sqrt((r ** 2).mean()))


def score():
    t = targets()
    runs = {v: run_dir(v) for v in VARIANTS}
    runs["matrix"] = MATRIX_RUN
    lrr = {k: model_lrr(rd) for k, rd in runs.items() if rd is not None}
    if "with2017" in lrr:
        d = (lrr["with2017"] - lrr["matrix"]).abs().max()
        print(f"with2017 vs matrix run: max |LRR difference| {d:.2e} m/yr")

    rows = []
    for v in VARIANTS:
        if v not in lrr:
            print(f"  missing {v}")
            continue
        md = json.loads(next(runs[v].glob("*_run_metadata.json")).read_text())["skill"]
        m = lrr[v]
        fp = list(FOOTPRINT_GIS)
        b_t, r_t = stats(m.loc[fp] - t.loc[fp, "target_lrr_m_yr"])
        b_r, r_r = stats(m.loc[fp] - t.loc[fp, "raw_mean_lrr_m_yr"])
        rows.append(dict(variant=v,
                         interior_bias_m_yr=float(md["mean_bias_interior_m_yr"]),
                         interior_rmse_m_yr=float(md["rmse_interior_m_yr"]),
                         gis6_15_bias_vs_target=b_t, gis6_15_rmse_vs_target=r_t,
                         gis6_15_bias_vs_raw=b_r, gis6_15_rmse_vs_raw=r_r,
                         gis6_15_model_mean_lrr=float(m.loc[fp].mean())))
    s = pd.DataFrame(rows)
    per = pd.DataFrame({"target_lrr_m_yr": t["target_lrr_m_yr"], "raw_mean_lrr_m_yr": t["raw_mean_lrr_m_yr"],
                        **{f"model_{k}": v for k, v in lrr.items() if k in VARIANTS}})
    if set(VARIANTS) <= set(lrr):
        per["with_minus_without"] = per["model_with2017"] - per["model_without2017"]
    (EXP_DIR / "tables").mkdir(parents=True, exist_ok=True)
    s.to_csv(EXP_DIR / "tables" / "scores.csv", index=False)
    per.to_csv(EXP_DIR / "tables" / "per_domain_lrr.csv")
    with pd.option_context("display.width", 220, "display.max_columns", 20):
        print(s.round(3).to_string(index=False))
        print(per.loc[1:25].round(2).to_string())


def main():
    if len(sys.argv) > 1 and sys.argv[1] == "_launch":
        _launch(sys.argv[2])
        return
    ap = argparse.ArgumentParser()
    ap.add_argument("action", choices=["run", "score"])
    a = ap.parse_args()
    if a.action == "run":
        jobs = [v for v in VARIANTS if run_dir(v) is None]
        with ThreadPoolExecutor(len(jobs) or 1) as ex:
            list(ex.map(launch, jobs))
    else:
        score()


if __name__ == "__main__":
    main()
