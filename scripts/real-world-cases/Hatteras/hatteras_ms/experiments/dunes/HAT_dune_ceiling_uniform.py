"""
Does one island-wide dune ceiling do as well as the per-cell ceilings, and at what height?

    python scripts/real-world-cases/Hatteras/hatteras_ms/experiments/dunes/HAT_dune_ceiling_uniform.py run
    python scripts/real-world-cases/Hatteras/hatteras_ms/experiments/dunes/HAT_dune_ceiling_uniform.py score

Replaces every domain's per-cell ceilings with one value, in-process only, at
5.0, 5.5, 6.0 and 6.5 m NAVD88; managed edgeBE runs of both periods against the
per-cell control; scored on overwash skill, dune crest, low-dune recovery,
dune rebuilds and shoreline. Details: scripts/real-world-cases/Hatteras/hatteras_ms/experiments/dunes/README.md.

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-10-13
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
import sys as _hat_sys  # noqa: E402
_HAT_DIR = next(_p for _p in _HERE.parents  if (_p / "hatteras_paths.py").exists())
if str(_HAT_DIR) not in _hat_sys.path:
    _hat_sys.path.insert(0, str(_HAT_DIR))
from hatteras_paths import WORKSPACE as PROJECT_ROOT  # noqa: E402
from cascade_pipeline.saved_runs import load_saved_cascade  # noqa: E402
sys.path.insert(0, str(_HAT_DIR))
sys.path.insert(0, str(_HAT_DIR / "input_prep" / "8-overwash-analysis" / "4-vs-model"))

# --- CONFIG ------------------------------------------------------------------
TAG = "dunes/2026-10-13-dune-ceiling-uniform"
EXP_DIR = PROJECT_ROOT / "output" / "raw_runs" / "experiments" / TAG
HINDCAST = _HAT_DIR / "hatteras_ms" / "HAT_hindcast.py"
MATRIX = PROJECT_ROOT / "output" / "raw_runs" / "matrix"

PERIODS = {1996: (2009, (1997, 2009)), 2009: (2025, (2010, 2025))}   # start: (end, storm years)
PRESET, SCENARIO = "edgeBE", "full_management"
MHW_NAVD88_M = 0.36
VARIANTS = {"percell": None, "uni50": 5.0, "uni55": 5.5, "uni60": 6.0, "uni65": 6.5}   # m NAVD88
LOW_QUARTILE = 25            # percentile of starting crests that defines a "low" dune cell
FLAT_M = 0.5                 # m above the berm: a cell this low counts as flattened
WORKERS = 4
# -----------------------------------------------------------------------------


# Every domain grows toward one ceiling: Barrier3D's scalar path, as with DuneCeilingFromStart false
def set_uniform(cascade, ceiling_navd88_m):
    for b in cascade.barrier3d:
        dmaxel = ceiling_navd88_m / 10.0 - MHW_NAVD88_M / 10.0        # dam MHW, as load_input stores Dmaxel
        b._DuneCeiling = None
        b._Dmaxel = dmaxel
        b._Dmax = dmaxel - b._BermEl


# Child process: wrap build_cascade to set the ceiling, then run the unchanged hindcast
def _launch(variant):
    import runpy
    import cascade_pipeline.hindcast as H
    value = VARIANTS[variant]
    if value is not None:
        orig = H.build_cascade

        def build_cascade(*a, **k):
            c = orig(*a, **k)
            set_uniform(c, value)
            print(f"EXPERIMENT uniform dune ceiling = {value:.2f} m NAVD88 (per-cell ceilings removed)", flush=True)
            return c

        H.build_cascade = build_cascade
    sys.path.insert(0, str(HINDCAST.parent))
    sys.argv = [str(HINDCAST)]
    runpy.run_path(str(HINDCAST), run_name="__main__")


def run_dir(start, variant):
    base = EXP_DIR / "runs" / variant / f"{start}_{PERIODS[start][0]}" / PRESET
    hits = sorted(d for d in base.glob("HAT_*") if any(d.glob("*_run_metadata.json"))) if base.exists() else []
    return hits[-1] if hits else None


# One run in its own process, logged, with its launch record
def launch(job):
    start, variant = job
    logs = EXP_DIR / "logs"
    logs.mkdir(parents=True, exist_ok=True)
    env = dict(os.environ)
    env.update({"HAT_IGNORE_SETTINGS": "1", "HAT_START_YEAR": str(start),
                "HAT_SOURCE_SINK_PRESET": PRESET, "HAT_SCENARIO": SCENARIO,
                "HAT_RUN_KIND": "experiment", "HAT_RUN_TAG": f"{TAG}/runs/{variant}",
                "HAT_SAVE_MODEL_STATE": "1", "HAT_MAKE_GIFS": "0", "HAT_SHOW_FIGURES": "0",
                "HAT_OVERWRITE": "1", "PYTHONIOENCODING": "utf-8"})
    t0 = time.time()
    log = logs / f"{variant}_{start}.log"
    with open(log, "w", encoding="utf-8") as fh:
        p = subprocess.run([sys.executable, str(_HERE), "_launch", variant],
                           env=env, cwd=HINDCAST.parent, stdout=fh, stderr=subprocess.STDOUT)
    rec = dict(start=start, variant=variant, ceiling_m_navd88=VARIANTS[variant], returncode=p.returncode,
               minutes=round((time.time() - t0) / 60, 1))
    with open(logs / "launches.jsonl", "a", encoding="utf-8") as fh:
        fh.write(json.dumps(rec) + "\n")
    print(f"{start} {variant:8s} exit {p.returncode}  {rec['minutes']} min", flush=True)


# Scores

def load(rd):
    return load_saved_cascade(rd / f"{rd.name}.npz")


# The geometry a run was made with: dom with the run's own buffer count. Runs before 2026-10-21 have 15 a side,
# the site config 30 since, so dom alone reads 15 domains off
def run_dom(c, dom):
    return dom.for_total(len(c.barrier3d))


# Per GIS domain, per cell, per saved year: crest height above the berm (m), the higher dune row
def cell_crests(c, dom):
    dom = run_dom(c, dom)
    return [np.asarray(c.barrier3d[dom.gis_to_pad(g)].DuneDomain, float).max(axis=2) * 10 for g in range(1, 91)]


def to_navd88(c, dom, above_berm):
    return above_berm + c.barrier3d[run_dom(c, dom).gis_to_pad(1)]._BermEl * 10 + MHW_NAVD88_M


# Observed against modelled overwash per image and domain, overwash_vs_model's rule on this period's storms
def overwash_skill(c, obs, dom, storm_years):
    import overwash_vs_model as ovm
    import storms_vs_overwash as SV
    dom = run_dom(c, dom)
    from site_layer import hat_env_forcings as env
    q = np.array([np.asarray(b.QowTS) for b in c.barrier3d]).T
    s = pd.read_csv(env.storm_summary_file(*storm_years), parse_dates=["StartTime", "EndTime"])
    s["dated"] = s.EndTime - ovm.GRACE
    largest = s.loc[s.groupby("time").Rhigh.idxmax()].set_index("time").dated
    imgs = ovm.images_in((storm_years[0], storm_years[1] + 1), obs)
    rows = []
    for _, im in imgs.iterrows():
        o = obs[obs.Obs_ID == im.Obs_ID].set_index("domain").overwash
        for gis in range(1, 91):
            p = dom.gis_to_pad(gis)
            yrs = [t for t in range(1, q.shape[0]) if q[t, p] > ovm.HEADLINE_THR]
            m = any(im["from"] < largest[t] <= im.date for t in yrs if t in largest)
            rows.append(dict(gis=gis, observed=o.get(gis, np.nan), model=int(m)))
    d = pd.DataFrame(rows).dropna(subset=["observed"])
    o, m = d.observed.astype(int), d.model.astype(int)
    d["outcome"] = np.select([(o == 1) & (m == 1), (o == 1) & (m == 0), (o == 0) & (m == 1)],
                             ["hit", "miss", "false_alarm"], "neither")
    return SV.skill(d)


def score():
    from site_layer import hat_overwash as ow
    from site_layer.hatteras_site_config import HATTERAS_DOMAINS as DOM
    obs = pd.read_csv(ow.OBSERVATIONS / "overwash_observations.csv")
    lidar_rd = run_dir(2009, "percell")
    lidar = np.array([np.median(cc[0]) for cc in cell_crests(load(lidar_rd), DOM)]) if lidar_rd else None
    rows = []
    for start, (end, storm_years) in PERIODS.items():
        for v, value in VARIANTS.items():
            rd = run_dir(start, v)
            if rd is None:
                print(f"  missing {start} {v}")
                continue
            c = load(rd)
            meta = json.loads((rd / f"{rd.name}_run_metadata.json").read_text(encoding="utf-8"))
            cc = cell_crests(c, DOM)
            first = np.concatenate([x[0] for x in cc])
            last = np.concatenate([x[-1] for x in cc])
            every = np.concatenate([x[1:].ravel() for x in cc])
            low = first <= np.percentile(first, LOW_QUARTILE)
            med_end = np.array([np.median(x[-1]) for x in cc])
            sk = overwash_skill(c, obs, DOM, storm_years)
            rm = pd.read_csv(rd / "tables" / "road_management.csv")
            row = dict(period=f"{start}_{end}", variant=v, ceiling_m_navd88=value,
                       skill=round(sk["skill"], 3), hit_rate=round(sk["hit_rate"], 3),
                       false_alarm_rate=round(sk["false_alarm_rate"], 3),
                       hits=sk["hits"], misses=sk["misses"], false_alarms=sk["false_alarms"],
                       crest_start_median_m_navd88=round(float(to_navd88(c, DOM, np.median(first))), 2),
                       crest_end_median_m_navd88=round(float(to_navd88(c, DOM, np.median(last))), 2),
                       low_quartile_start_m_navd88=round(float(to_navd88(c, DOM, np.median(first[low]))), 2),
                       low_quartile_end_m_navd88=round(float(to_navd88(c, DOM, np.median(last[low]))), 2),
                       share_cell_years_flat=round(float((every < FLAT_M).mean()), 3),
                       alongshore_crest_sd_end_m=round(float(np.std(med_end)), 2),
                       dunes_rebuilt=int(rm.dunes_rebuilt.sum()),
                       rmse_interior_m_yr=float(meta["skill"]["rmse_interior_m_yr"]))
            if start == 1996 and lidar is not None:
                row.update(crest_end_minus_lidar_m=round(float(np.median(med_end - lidar)), 2),
                           crest_end_vs_lidar_r=round(float(np.corrcoef(med_end, lidar)[0, 1]), 2))
            rows.append(row)
            print(f"  scored {start} {v}", flush=True)
    t = pd.DataFrame(rows)
    (EXP_DIR / "tables").mkdir(parents=True, exist_ok=True)
    t.to_csv(EXP_DIR / "tables" / "scores.csv", index=False)
    with pd.option_context("display.width", 300, "display.max_columns", 40):
        print(t.to_string(index=False))


def main():
    if len(sys.argv) > 1 and sys.argv[1] == "_launch":
        _launch(sys.argv[2])
        return
    ap = argparse.ArgumentParser()
    ap.add_argument("action", choices=["run", "score"])
    ap.add_argument("--workers", type=int, default=WORKERS)
    a = ap.parse_args()
    if a.action == "run":
        jobs = [(s, v) for s in PERIODS for v in VARIANTS if run_dir(s, v) is None]
        with ThreadPoolExecutor(a.workers) as ex:
            list(ex.map(launch, jobs))
    else:
        score()


if __name__ == "__main__":
    main()
