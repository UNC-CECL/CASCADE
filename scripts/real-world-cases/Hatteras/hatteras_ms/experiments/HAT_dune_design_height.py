"""
Does raising the dune design height from 3.0 to 3.94 m MHW (4.3 m NAVD88) change the calibration run?

    python scripts/hatteras_ms/experiments/HAT_dune_design_height.py run
    python scripts/hatteras_ms/experiments/HAT_dune_design_height.py score

Runs the unchanged hindcast for the 1996-2009 managed edgeBE run twice, the
control and the design height set in-process only; scored on overwash skill,
dune crest against the 2009 lidar, dune rebuilds and shoreline. Details:
scripts/hatteras_ms/experiments/README.md.

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
sys.path.insert(0, str(_HAT_DIR))
sys.path.insert(0, str(_HAT_DIR / "input_prep" / "8-overwash-analysis" / "4-vs-model"))

# --- CONFIG ------------------------------------------------------------------
TAG = "dunes/2026-10-13-dune-design-height"
EXP_DIR = PROJECT_ROOT / "output" / "raw_runs" / "experiments" / TAG
HINDCAST = _HAT_DIR / "hatteras_ms" / "HAT_hindcast.py"
MATRIX = PROJECT_ROOT / "output" / "raw_runs" / "matrix"

START, END = 1996, 2009                 # the calibration period
STORM_YEARS = (1997, 2009)              # its storm series
TEST_START, TEST_END = 2009, 2025       # whose starting dunes are the 2009 lidar
PRESET, SCENARIO = "edgeBE", "full_management"
MHW_NAVD88_M = 0.36
VARIANTS = {"control": None,                              # the runner's own 3.0 m MHW
            "design394": 4.3 - MHW_NAVD88_M}              # Velasquez-Montoya et al. (2021) NC-12 safe crest
WORKERS = 4
# -----------------------------------------------------------------------------


# Child process: override the design height on the build_cascade call, then run the unchanged hindcast
def _launch(variant):
    import runpy
    import cascade_pipeline.hindcast as H
    value = VARIANTS[variant]
    if value is not None:
        orig = H.build_cascade

        def build_cascade(*a, **k):
            k["dune_design_elevation"] = [value] * len(k["dune_design_elevation"])
            print(f"EXPERIMENT dune_design_elevation = {value:.2f} m MHW (runner constant overridden)", flush=True)
            return orig(*a, **k)

        H.build_cascade = build_cascade
    sys.path.insert(0, str(HINDCAST.parent))
    sys.argv = [str(HINDCAST)]
    runpy.run_path(str(HINDCAST), run_name="__main__")


def run_dir(variant):
    base = EXP_DIR / "runs" / variant / f"{START}_{END}" / PRESET
    hits = sorted(d for d in base.glob("HAT_*") if any(d.glob("*_run_metadata.json"))) if base.exists() else []
    return hits[-1] if hits else None


def matrix_dir():
    return next((MATRIX / f"{START}_{END}" / PRESET).glob("HAT_*_road_bdm_*nogroin"))


# One run in its own process, logged, with its launch record
def launch(variant):
    logs = EXP_DIR / "logs"
    logs.mkdir(parents=True, exist_ok=True)
    env = dict(os.environ)
    env.update({"HAT_IGNORE_SETTINGS": "1", "HAT_START_YEAR": str(START),
                "HAT_SOURCE_SINK_PRESET": PRESET, "HAT_SCENARIO": SCENARIO,
                "HAT_RUN_KIND": "experiment", "HAT_RUN_TAG": f"{TAG}/runs/{variant}",
                "HAT_SAVE_MODEL_STATE": "1", "HAT_MAKE_GIFS": "0", "HAT_SHOW_FIGURES": "0",
                "HAT_OVERWRITE": "1", "PYTHONIOENCODING": "utf-8"})
    t0 = time.time()
    log = logs / f"{variant}.log"
    with open(log, "w", encoding="utf-8") as fh:
        p = subprocess.run([sys.executable, str(_HERE), "_launch", variant],
                           env=env, cwd=HINDCAST.parent, stdout=fh, stderr=subprocess.STDOUT)
    rec = dict(variant=variant, design_m_mhw=VARIANTS[variant], returncode=p.returncode,
               minutes=round((time.time() - t0) / 60, 1))
    with open(logs / "launches.jsonl", "a", encoding="utf-8") as fh:
        fh.write(json.dumps(rec) + "\n")
    print(f"{variant:10s} exit {p.returncode}  {rec['minutes']} min", flush=True)


# Scores

def load(rd):
    return np.load(rd / f"{rd.name}.npz", allow_pickle=True)["cascade"][0]


# Median crest per GIS domain (m MHW) for every saved year: highest dune row, median over the cells
def crests(c, dom):
    out = []
    for gis in range(1, 91):
        b = c.barrier3d[dom.gis_to_pad(gis)]
        dd = np.asarray(b.DuneDomain, float)
        out.append(np.median(dd.max(axis=2), axis=1) * 10 + b._BermEl * 10)
    return np.array(out)                                   # (90, years)


# Observed against modelled overwash per image and domain, overwash_vs_model's rule on this period's storms
def overwash_cells(c, obs, dom):
    import overwash_vs_model as ovm
    import storms_vs_overwash as SV
    from site_layer import hat_env_forcings as env
    q = np.array([np.asarray(b.QowTS) for b in c.barrier3d]).T
    s = pd.read_csv(env.storm_summary_file(*STORM_YEARS), parse_dates=["StartTime", "EndTime"])
    s["dated"] = s.EndTime - ovm.GRACE
    largest = s.loc[s.groupby("time").Rhigh.idxmax()].set_index("time").dated
    imgs = ovm.images_in((STORM_YEARS[0], STORM_YEARS[1] + 1), obs)
    rows = []
    for _, im in imgs.iterrows():
        o = obs[obs.Obs_ID == im.Obs_ID].set_index("domain").overwash
        for gis in range(1, 91):
            p = dom.gis_to_pad(gis)
            yrs = [t for t in range(1, q.shape[0]) if q[t, p] > ovm.HEADLINE_THR]
            m = any(im["from"] < largest[t] <= im.date for t in yrs if t in largest)
            rows.append(dict(obs_id=im.Obs_ID, gis=gis, observed=o.get(gis, np.nan), model=int(m)))
    d = pd.DataFrame(rows).dropna(subset=["observed"])
    o, m = d.observed.astype(int), d.model.astype(int)
    d["outcome"] = np.select([(o == 1) & (m == 1), (o == 1) & (m == 0), (o == 0) & (m == 1)],
                             ["hit", "miss", "false_alarm"], "neither")
    return SV.skill(d)


def score():
    from site_layer import hat_overwash as ow
    from site_layer.hatteras_site_config import HATTERAS_DOMAINS as DOM
    obs = pd.read_csv(ow.OBSERVATIONS / "overwash_observations.csv")
    test = next((MATRIX / f"{TEST_START}_{TEST_END}" / PRESET).glob("HAT_*_road_bdm_*nogroin"))
    lidar = crests(load(test), DOM)[:, 0]                 # the 2009 run's starting crests = 2009 lidar
    runs = {"matrix": matrix_dir(), **{v: run_dir(v) for v in VARIANTS}}
    rows, shore = [], {}
    for name, rd in runs.items():
        if rd is None:
            print(f"  missing {name}")
            continue
        c = load(rd)
        meta = json.loads((rd / f"{rd.name}_run_metadata.json").read_text(encoding="utf-8"))
        cr = crests(c, DOM)
        rm = pd.read_csv(rd / "tables" / "road_management.csv")
        design = [r._dune_design_elevation_TS[0] for r in c.roadways if r is not None]
        shore[name] = np.load(rd / f"{rd.name}_shoreline_matrix.npy")
        rows.append(dict(run=name, design_applied_m_mhw=round(float(np.median(design)), 2),
                         **{k: round(float(v), 3) for k, v in overwash_cells(c, obs, DOM).items()},
                         crest_end_median_m_mhw=round(float(np.median(cr[:, -1])), 2),
                         crest_end_minus_lidar_m=round(float(np.median(cr[:, -1] - lidar)), 2),
                         crest_end_vs_lidar_r=round(float(np.corrcoef(cr[:, -1], lidar)[0, 1]), 2),
                         dunes_rebuilt=int(rm.dunes_rebuilt.sum()),
                         rmse_interior_m_yr=meta.get("skill", {}).get("rmse_interior_m_yr"),
                         endpoint_rmse_interior_m_yr=meta.get("skill", {}).get("endpoint_rmse_interior_m_yr")))
        print(f"  scored {name}", flush=True)
    t = pd.DataFrame(rows)
    if "control" in shore:
        for name, m in shore.items():
            t.loc[t.run == name, "max_abs_shore_diff_vs_control_m"] = round(float(np.nanmax(np.abs(m - shore["control"]))), 3)
    (EXP_DIR / "tables").mkdir(parents=True, exist_ok=True)
    t.to_csv(EXP_DIR / "tables" / "scores.csv", index=False)
    with pd.option_context("display.width", 260, "display.max_columns", 30):
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
        jobs = [v for v in VARIANTS if run_dir(v) is None]
        with ThreadPoolExecutor(a.workers) as ex:
            list(ex.map(launch, jobs))
    else:
        score()


if __name__ == "__main__":
    main()
