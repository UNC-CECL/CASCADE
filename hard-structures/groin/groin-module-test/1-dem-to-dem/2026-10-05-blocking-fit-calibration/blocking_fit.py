#!/usr/bin/env python3
"""
Fit the blocking groin (b, f) on the 1996-2009 calibration period, then check it on 2009-2025.

    python blocking_fit.py grid  [--streams 4]       # the calibration grid
    python blocking_fit.py score                     # date RMSE per cell -> grid_scores.csv
    python blocking_fit.py test --b 0.6 --f 0.4      # the chosen cell on the test period

Every run: full management, the solved edgeBE ends, relocations off, failure instant from
the 2004 step. Scored by the D5-D6 gap change at the wet/dry photo dates inside the window
(the 2026-09-29 date RMSE). Runs file under experiments/groin/2026-10-05-blocking-fit-dem-to-dem/.

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-10-05
"""
from __future__ import annotations

import argparse
import os
import queue
import subprocess
import sys
import threading
import time
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
REPO = next(p for p in HERE.parents if (p / "pyproject.toml").exists())
sys.path[:0] = [str(HERE.parents[1] / "0-solver-audit" / "2026-09-29-option-a-real-planform"),
                str(REPO / "scripts" / "hatteras_ms" / "groin-sweep"),
                str(REPO / "scripts" / "hatteras_ms"), str(REPO / "scripts")]
from score_instant_grid import observed_series  # noqa: E402
from site_layer.hatteras_site_config import HATTERAS_PERIODS, run_years  # noqa: E402

# --- CONFIG ------------------------------------------------------------------
RUNNER = REPO / "scripts" / "hatteras_ms" / "HAT_hindcast_1984_2024.py"
STUDY = "groin/2026-10-05-blocking-fit-dem-to-dem"
RAW = REPO / "output" / "raw_runs"
LOGS = HERE / "logs"
B_VALUES = (0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9)
F_VALUES = (0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.8)
CALIBRATION, TEST = 1996, 2009
UP, DOWN = 15 + 6 - 1, 15 + 5 - 1      # GIS 6 / GIS 5 in the padded array
# -----------------------------------------------------------------------------


def member(b, f):
    return f"b{b:.2f}_f{f:.1f}"


def window(period):
    return f"{period}_{HATTERAS_PERIODS[period]['end_year']}"


# The matrix files of a cell, or the no-groin edgeBE baseline when b is None
def shoreline_file(period, b=None, f=None):
    if b is None:
        fill = "_nourish" if HATTERAS_PERIODS[period]["enable_nourishment"] else ""
        name = f"HAT_{window(period)}_edgeBE_offsetmetres_road_bdm{fill}_nogroin"
        return RAW / "matrix" / window(period) / "edgeBE" / name / f"{name}_shoreline_matrix.npy"
    hits = sorted((RAW / "experiments" / STUDY / member(b, f)).glob(
        f"{window(period)}/edgeBE/*/*_shoreline_matrix.npy"))
    return hits[-1] if hits else None


def run(period, b, f):
    env = {k: v for k, v in os.environ.items() if not k.startswith("HAT_")}
    env.update({
        "HAT_IGNORE_SETTINGS": "1", "HAT_START_YEAR": str(period),
        "HAT_SOURCE_SINK_PRESET": "edgeBE", "HAT_SCENARIO": "full_management",
        "HAT_RELOCATIONS": "false", "HAT_GROIN_ENABLED": "true", "HAT_GROIN_KIND": "blocking",
        "HAT_GROIN_BLOCKING_FRACTION": str(b), "HAT_GROIN_DETERIORATION_FRACTION": str(f),
        "HAT_RUN_KIND": "experiment", "HAT_RUN_TAG": f"{STUDY}/{member(b, f)}",
        "HAT_SAVE_MODEL_STATE": "false", "HAT_SHOW_FIGURES": "false", "HAT_MAKE_GIFS": "false",
        "HAT_OVERWRITE": "false", "MPLBACKEND": "Agg", "PYTHONIOENCODING": "utf-8",
    })
    LOGS.mkdir(exist_ok=True)
    log = LOGS / f"{window(period)}_{member(b, f)}.log"
    t0 = time.time()
    with open(log, "w", encoding="utf-8") as fh:
        code = subprocess.call([sys.executable, str(RUNNER)], cwd=REPO, env=env,
                               stdout=fh, stderr=subprocess.STDOUT)
    print(f"{time.strftime('%H:%M:%S')}  {window(period)} {member(b, f)}: exit {code} "
          f"({(time.time() - t0) / 60:.1f} min)", flush=True)
    return code


def run_all(cells, streams):
    todo = [c for c in cells if shoreline_file(*c) is None]
    print(f"{len(cells) - len(todo)} of {len(cells)} done; running {len(todo)} in {streams} streams",
          flush=True)
    jobs, failed = queue.Queue(), []
    for c in todo:
        jobs.put(c)

    def worker():
        while True:
            try:
                c = jobs.get_nowait()
            except queue.Empty:
                return
            if run(*c):
                failed.append(c)

    threads = [threading.Thread(target=worker) for _ in range(streams)]
    for t in threads:
        t.start()
    for t in threads:
        t.join()
    print(f"DONE: {len(todo) - len(failed)}/{len(todo)} clean; failed {failed}", flush=True)
    return not failed


# Gap change at the photo dates inside the window, model against observed
def score(path, period, obs):
    x = np.load(path)
    g = x[:, DOWN] - x[:, UP]
    g = g - g[0]
    n = run_years(period)
    yrs = [y for y in obs.index if period < y <= period + n]
    g0 = np.interp(period, obs.index, obs.values)
    o = np.array([obs[y] for y in yrs]) - g0
    m = np.array([g[y - period] for y in yrs])
    return dict(n_dates=len(yrs), rmse_dates=float(np.sqrt(np.mean((m - o) ** 2))),
                dates=" ".join(map(str, yrs)),
                model_at_dates=" ".join(f"{v:+.0f}" for v in m),
                observed_at_dates=" ".join(f"{v:+.0f}" for v in o))


def cmd_score(_a):
    obs = observed_series()
    rows = []
    for period in (CALIBRATION, TEST):
        base = shoreline_file(period)
        if base.is_file():
            rows.append(dict(period=period, b=0.0, f=np.nan, **score(base, period, obs)))
        for b in B_VALUES:
            for f in F_VALUES:
                p = shoreline_file(period, b, f)
                if p is not None:
                    rows.append(dict(period=period, b=b, f=f, **score(p, period, obs)))
    t = pd.DataFrame(rows)
    t.to_csv(HERE / "grid_scores.csv", index=False)
    pd.set_option("display.width", 200)
    for period, sub in t.groupby("period"):
        print(f"\n=== {window(period)}: dates {sub.dates.iloc[0]}, observed {sub.observed_at_dates.iloc[0]} ===")
        b0 = sub[sub.b == 0.0]
        if len(b0):
            print(f"no groin: {b0.model_at_dates.iloc[0]}  RMSE {b0.rmse_dates.iloc[0]:.1f}")
        g = sub[sub.b > 0]
        if len(g) > 1:
            print(g.pivot(index="b", columns="f", values="rmse_dates").to_string(float_format="%.1f"))
        for r in g.nsmallest(5, "rmse_dates").itertuples():
            print(f"  b {r.b:.2f} f {r.f:.1f}: RMSE {r.rmse_dates:.1f}  model {r.model_at_dates}")


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    sub = ap.add_subparsers(dest="cmd", required=True)
    g = sub.add_parser("grid")
    g.add_argument("--streams", type=int, default=4)
    sub.add_parser("score")
    t = sub.add_parser("test")
    t.add_argument("--b", type=float, required=True)
    t.add_argument("--f", type=float, required=True)
    a = ap.parse_args()
    if a.cmd == "grid":
        ok = run_all([(CALIBRATION, b, f) for b in B_VALUES for f in F_VALUES], a.streams)
        sys.exit(0 if ok else 1)
    if a.cmd == "test":
        sys.exit(0 if run_all([(TEST, a.b, a.f)], 1) else 1)
    cmd_score(a)


if __name__ == "__main__":
    main()
