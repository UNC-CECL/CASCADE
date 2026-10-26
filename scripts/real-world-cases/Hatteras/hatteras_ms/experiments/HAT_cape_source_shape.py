"""
Which shape should the Cape Point sediment source take: all at GIS 1, or tapered over GIS 1-5?

    python scripts/hatteras_ms/experiments/HAT_cape_source_shape.py run
    python scripts/hatteras_ms/experiments/HAT_cape_source_shape.py score

Why (2026-10-21): the cape's own dune growth rate (0.10-0.20) was dropped. It
fitted the 2014-2017 cape surveys but overbuilt the 1996-2009 cape dunes by
3-4 m and barely moved GIS 1-2 (dunes/2026-10-21-cape-growth-rate). Cape Point
accretes because sand converges there, which the periodic BRIE domain cannot
carry, so the cape gets its sand from the source/sink instead. The end solve
puts all of it at GIS 1; this asks whether spreading it over the cape fits
better.

Shapes, each a total source Q (m/yr, summed over the cape domains) split by
fixed weights:
    gis1_only           1                   (today's end-solve shape)
    linear_5_4_3_2_1    5:4:3:2:1 over GIS 1-5
    halving_16_8_4_2_1  16:8:4:2:1 over GIS 1-5
Per shape and window, Q is solved by secant so the mean net change over GIS
1-5 matches CoastSat's (domain means, raw). Same criterion for every shape, so
they are compared at equal cape-zone sand balance. GIS 90 keeps the config's
solved rate. edgeBE through HAT_BE_OVERRIDE (no main-code change), full
management, groin on, relocations off, the runner's island-wide growth rate.
Scores: per-domain change at GIS 1-5 against CoastSat, cape RMSE, GIS 6-10,
interior RMSE and bias.

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-10-21
"""
from __future__ import annotations

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
_HAT_DIR = next(p for p in _HERE.parents  if (p / "hatteras_paths.py").exists())
if str(_HAT_DIR) not in _hat_sys.path:
    _hat_sys.path.insert(0, str(_HAT_DIR))
from hatteras_paths import WORKSPACE as PROJECT_ROOT  # noqa: E402
HINDCAST = _HAT_DIR / "hatteras_ms" / "HAT_hindcast.py"
sys.path.insert(0, str(_HAT_DIR))
from site_layer.hat_observed_rates import net_change_dir  # noqa: E402
from site_layer.hatteras_site_config import HATTERAS_BE_EDGE_ONLY, HATTERAS_DOMAINS as D, HATTERAS_PERIODS  # noqa: E402

# --- CONFIG ------------------------------------------------------------------
TAG = "source-sink/2026-10-21-cape-source-shape"
EXP_DIR = PROJECT_ROOT / "output" / "raw_runs" / "experiments" / TAG
SHAPES = {"gis1_only": (1,), "linear_5_4_3_2_1": (5, 4, 3, 2, 1), "halving_16_8_4_2_1": (16, 8, 4, 2, 1)}
CAPE = (1, 2, 3, 4, 5)
STARTS = (1996, 2009)
Q_START = {1996: (5.0, 15.0), 2009: (50.0, 100.0)}     # first two probes of the total source, m/yr
TOL_M = 1.0                                             # cape-mean residual that counts as solved
MAX_STEPS = 7
WORKERS = 4
# -----------------------------------------------------------------------------


def window(start):
    return f"{start}_{HATTERAS_PERIODS[start]['end_year']}"


def rates(shape, q):
    w = np.array(SHAPES[shape], float)
    return {g: q * wi / w.sum() for g, wi in zip(CAPE, w)}


def override(shape, q, start):
    r = rates(shape, q)
    r[90] = HATTERAS_BE_EDGE_ONLY[start][1]
    return ",".join(f"{g}={v:.4f}" for g, v in r.items())


def run_dir(shape, start, step):
    base = EXP_DIR / f"{shape}_step{step}" / window(start) / "edgeBE"   # tag stays within 4 parts
    hits = sorted(d for d in base.glob("HAT_*") if any(d.glob("*_run_metadata.json"))) if base.exists() else []
    return hits[-1] if hits else None


def launch(shape, start, step, q):
    if run_dir(shape, start, step) is not None:
        return run_dir(shape, start, step)
    logs = EXP_DIR / "logs"
    logs.mkdir(parents=True, exist_ok=True)
    env = {k: v for k, v in os.environ.items() if not k.startswith("HAT_")}
    env.update({"HAT_IGNORE_SETTINGS": "1", "HAT_START_YEAR": str(start), "HAT_SOURCE_SINK_PRESET": "edgeBE",
                "HAT_BE_OVERRIDE": override(shape, q, start), "HAT_SCENARIO": "full_management",
                "HAT_RELOCATIONS": "False", "HAT_GROIN_ENABLED": "True",
                "HAT_RUN_KIND": "experiment", "HAT_RUN_TAG": f"{TAG}/{shape}_step{step}",
                "HAT_SAVE_MODEL_STATE": "0", "HAT_MAKE_GIFS": "0", "HAT_SHOW_FIGURES": "0",
                "HAT_OVERWRITE": "1", "PYTHONIOENCODING": "utf-8", "MPLBACKEND": "Agg"})
    with open(logs / f"{shape}_{start}_step{step}.log", "w", encoding="utf-8") as fh:
        code = subprocess.run([sys.executable, str(HINDCAST)], env=env, cwd=HINDCAST.parent,
                              stdout=fh, stderr=subprocess.STDOUT).returncode
    if code:
        raise RuntimeError(f"{shape} {start} step{step} failed ({code})")
    return run_dir(shape, start, step)


# Observed net change per GIS domain: the mean of its transects (raw)
def observed(start):
    w = window(start)
    t = pd.read_csv(net_change_dir(*map(int, w.split("_"))) / f"transect_net_change_{w}.csv")
    return t.groupby("domain_number").net_change_m.mean()


# Model net change per GIS domain 1-90, seaward +
def model_change(rd):
    m = np.load(next(rd.glob("*_shoreline_matrix.npy")))
    chg = -(m[-1] - m[0])
    return pd.Series([chg[D.gis_to_pad(g)] for g in range(1, 91)], index=range(1, 91))


# Secant on Q until the GIS 1-5 mean closes
def solve(job):
    shape, start = job
    target = observed(start).loc[list(CAPE)].mean()
    hist = []
    for step in range(MAX_STEPS):
        if step < 2:
            q = Q_START[start][step]
        else:
            (q0, r0), (q1, r1) = hist[-2], hist[-1]
            q = q1 - r1 * (q1 - q0) / (r1 - r0) if r1 != r0 else q1
        t0 = time.time()
        rd = launch(shape, start, step, q)
        resid = model_change(rd).loc[list(CAPE)].mean() - target
        hist.append((q, resid))
        line = f"{start} {shape} step{step} Q {q:+.3f}  cape-mean residual {resid:+.2f} m  ({(time.time() - t0) / 60:.1f} min)"
        print(line, flush=True)
        with open(EXP_DIR / "solve.txt", "a", encoding="utf-8") as fh:
            fh.write(line + "\n")
        if abs(resid) < TOL_M:
            break
    return shape, start, hist


def best(shape, start):
    rows = []
    for step in range(MAX_STEPS):
        rd = run_dir(shape, start, step)
        if rd is None:
            continue
        meta = json.loads(next(rd.glob("*_run_metadata.json")).read_text(encoding="utf-8"))
        q = sum(v for g, v in meta_rates(rd).items() if g in CAPE)
        rows.append((abs(model_change(rd).loc[list(CAPE)].mean() - observed(start).loc[list(CAPE)].mean()), q, rd, meta))
    return min(rows, key=lambda r: r[0]) if rows else None


def meta_rates(rd):
    import re
    txt = next(rd.glob("*_run_metadata.txt")).read_text(encoding="utf-8", errors="replace")
    return {int(g): float(v) for g, v in re.findall(r"be_rate_gis(\d+)_m_yr\s+=\s+([-\d.]+)", txt)}


def score():
    obs = {s: observed(s) for s in STARTS}
    rows = []
    for s in STARTS:
        for shape in SHAPES:
            b = best(shape, s)
            if b is None:
                continue
            resid, q, rd, meta = b
            c = model_change(rd)
            err = c.loc[list(CAPE)] - obs[s].loc[list(CAPE)]
            rows.append(dict(window=window(s), shape=shape, total_source_m_yr=round(q, 2),
                             cape_mean_resid_m=round(resid, 2),
                             cape_change_m=" ".join(f"{c.loc[g]:+.0f}" for g in CAPE),
                             cape_rmse_m=round(float(np.sqrt((err ** 2).mean())), 1),
                             gis6_10_change_m=" ".join(f"{c.loc[g]:+.0f}" for g in range(6, 11)),
                             shoreline_rmse=round(float(meta["skill"]["rmse_interior_m_yr"]), 3),
                             shoreline_bias=round(float(meta["skill"]["mean_bias_interior_m_yr"]), 3),
                             run=str(rd.relative_to(EXP_DIR))))
        rows.append(dict(window=window(s), shape="CoastSat",
                         cape_change_m=" ".join(f"{obs[s].loc[g]:+.0f}" for g in CAPE),
                         gis6_10_change_m=" ".join(f"{obs[s].loc[g]:+.0f}" for g in range(6, 11))))
    t = pd.DataFrame(rows)
    (EXP_DIR / "tables").mkdir(parents=True, exist_ok=True)
    t.to_csv(EXP_DIR / "tables" / "cape_source_shape_scores.csv", index=False)
    with pd.option_context("display.width", 320, "display.max_columns", 30, "display.max_colwidth", 60):
        print(t.drop(columns="run").to_string(index=False))


def main():
    if len(sys.argv) > 1 and sys.argv[1] == "run":
        EXP_DIR.mkdir(parents=True, exist_ok=True)
        jobs = [(shape, s) for s in (2009, 1996) for shape in SHAPES]
        with ThreadPoolExecutor(WORKERS) as ex:
            list(ex.map(solve, jobs))
    score()


if __name__ == "__main__":
    main()
