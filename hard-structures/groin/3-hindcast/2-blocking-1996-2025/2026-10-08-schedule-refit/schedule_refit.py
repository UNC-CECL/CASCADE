#!/usr/bin/env python3
"""
Refit the blocking groin (b, f) on 1996-2009 under three failure schedules, scored on the annual CoastSat gap.

    python schedule_refit.py grid  [--streams 4]                 # 3 schedules x the b-f grid
    python schedule_refit.py score                               # -> grid_scores.csv
    python schedule_refit.py test --schedule ramp1996 --b 0.6 --f 0.6   # one cell on 2009-2025

Schedules: the runner's instant failure at the 2004 step; instant from 1996, the first
year after the 1995 last repair; a linear ramp from full strength in 1995 to the floor in
2003, so 1996 is the first weakened year. The unchanged runner is driven in a child
process with BlockingGroinCallback patched to the schedule (rule: experiments don't
touch main code); the runner's own report still prints its default schedule, so each
run gets a schedule.json. Every run: full management, the solved edgeBE ends,
relocations off. Primary score: RMSE of the modelled GIS 5|6 gap change against the
annual CoastSat gap change (1-observations/gap_across_groins), model mid-year vs the
calendar-year mean, both relative to the start (model t = 0, CoastSat over the
DEM-centred start window). Secondary: the 2026-10-05 photo-date RMSE.

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-10-08
"""
from __future__ import annotations

import argparse
import json
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
GROIN = REPO / "hard-structures" / "groin"
sys.path[:0] = [str(GROIN / "1-observations" / "gap_across_groins"),
                str(GROIN / "2-module-tests" / "3-real-planform"),
                str(REPO / "scripts" / "hatteras_ms" / "groin-sweep"),
                str(REPO / "scripts" / "hatteras_ms"), str(REPO / "scripts")]

# --- CONFIG ------------------------------------------------------------------
RUNNER = REPO / "scripts" / "hatteras_ms" / "HAT_hindcast_1984_2024.py"
STUDY = "groin/2026-10-08-schedule-refit"
RAW = REPO / "output" / "raw_runs"
LOGS = HERE / "logs"
INSTALL = 1969
# name: (deterioration_delay_years from install, mode, ramp years)
SCHEDULES = {
    "instant2004": (2004 - INSTALL, "instant", 0.0),       # the runner's schedule now
    "instant1996": (1996 - INSTALL, "instant", 0.0),       # first year after the 1995 repair
    "ramp1996": (1995 - INSTALL, "linear_ramp", 8.0),      # full in 1995, floor in 2003
}
B_VALUES = (0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9)
F_VALUES = (0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.8)
CALIBRATION, TEST = 1996, 2009
START_WINDOW = {1996: ("1995-10-12", "1997-10-12"), 2009: ("2008-08-17", "2010-08-17")}
UP, DOWN = 15 + 6 - 1, 15 + 5 - 1      # GIS 6 / GIS 5 in the padded array
# -----------------------------------------------------------------------------


def member(b, f):
    return f"b{b:.2f}_f{f:.1f}"


def window(period):
    from site_layer.hatteras_site_config import HATTERAS_PERIODS
    return f"{period}_{HATTERAS_PERIODS[period]['end_year']}"


def run_dir(schedule, period, b, f):
    hits = sorted((RAW / "experiments" / STUDY / schedule / member(b, f)).glob(
        f"{window(period)}/edgeBE/*/"))
    return hits[-1] if hits else None


def shoreline_file(schedule, period, b, f):
    d = run_dir(schedule, period, b, f)
    hits = sorted(d.glob("*_shoreline_matrix.npy")) if d else []
    return hits[-1] if hits else None


def baseline_file(period):
    from site_layer.hatteras_site_config import HATTERAS_PERIODS
    fill = "_nourish" if HATTERAS_PERIODS[period]["enable_nourishment"] else ""
    name = f"HAT_{window(period)}_edgeBE_offsetmetres_road_bdm{fill}_nogroin"
    return RAW / "matrix" / window(period) / "edgeBE" / name / f"{name}_shoreline_matrix.npy"


# Child process: patch the callback to the schedule, then run the unchanged runner
def _launch(schedule):
    import runpy
    import cascade.groin as cg
    delay, mode, ramp = SCHEDULES[schedule]

    class Scheduled(cg.BlockingGroinCallback):
        def __init__(self, *a, **kw):
            kw.update(deterioration_delay_years=delay, deterioration_mode=mode,
                      deterioration_ramp_years=ramp)
            super().__init__(*a, **kw)

    cg.BlockingGroinCallback = Scheduled
    sys.path.insert(0, str(RUNNER.parent))
    sys.argv = [str(RUNNER)]
    runpy.run_path(str(RUNNER), run_name="__main__")


def run(schedule, period, b, f):
    env = {k: v for k, v in os.environ.items() if not k.startswith("HAT_")}
    env.update({
        "HAT_IGNORE_SETTINGS": "1", "HAT_START_YEAR": str(period),
        "HAT_SOURCE_SINK_PRESET": "edgeBE", "HAT_SCENARIO": "full_management",
        "HAT_RELOCATIONS": "false", "HAT_GROIN_ENABLED": "true", "HAT_GROIN_KIND": "blocking",
        "HAT_GROIN_BLOCKING_FRACTION": str(b), "HAT_GROIN_DETERIORATION_FRACTION": str(f),
        "HAT_RUN_KIND": "experiment", "HAT_RUN_TAG": f"{STUDY}/{schedule}/{member(b, f)}",
        "HAT_SAVE_MODEL_STATE": "false", "HAT_SHOW_FIGURES": "false", "HAT_MAKE_GIFS": "false",
        "HAT_OVERWRITE": "false", "MPLBACKEND": "Agg", "PYTHONIOENCODING": "utf-8",
    })
    LOGS.mkdir(exist_ok=True)
    log = LOGS / f"{schedule}_{window(period)}_{member(b, f)}.log"
    t0 = time.time()
    with open(log, "w", encoding="utf-8") as fh:
        code = subprocess.call([sys.executable, str(Path(__file__)), "_launch", schedule], cwd=REPO,
                               env=env, stdout=fh, stderr=subprocess.STDOUT)
    d = run_dir(schedule, period, b, f)
    if code == 0 and d is not None:
        delay, mode, ramp = SCHEDULES[schedule]
        (d / "schedule.json").write_text(json.dumps(dict(
            schedule=schedule, install_year=INSTALL, deterioration_delay_years=delay,
            deterioration_year=INSTALL + delay, deterioration_mode=mode,
            deterioration_ramp_years=ramp, blocking_b=b, deterioration_f=f,
            note="patched in-process by schedule_refit.py; the runner's report prints its default"),
            indent=2))
    print(f"{time.strftime('%H:%M:%S')}  {schedule} {window(period)} {member(b, f)}: exit {code} "
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


# The observed CoastSat GIS 6 minus GIS 5 gap change: annual means minus the start-window mean
def coastsat_change(period):
    import coastsat_gap_across_groins as cg
    up, down = cg.sides()["domain"]
    lo, hi = (pd.Timestamp(d, tz="UTC") for d in START_WINDOW[period])

    def side(tids):
        ann, win = {}, {}
        for t in tids:
            d = pd.read_csv(cg.SERIES_DIR / f"{cg.SITE}_timeseries" / f"{t}.csv", header=0)
            d.columns = ["date", "chainage_m"] + list(d.columns[2:])
            d["date"] = pd.to_datetime(d["date"], utc=True)
            d["chainage_m"] = pd.to_numeric(d["chainage_m"], errors="coerce")
            ref = cg.annual(t).loc[cg.FIT[0]:cg.FIT[1]].mean()
            ann[t] = cg.annual(t) - ref
            w = d[(d["date"] >= lo) & (d["date"] <= hi)]["chainage_m"]
            win[t] = w.mean() - ref
        a = pd.DataFrame(ann)
        ok = a.notna().sum(axis=1) >= np.ceil(cg.MIN_SHARE * len(tids))
        return a[ok].mean(axis=1), pd.Series(win).mean()

    su, wu = side(up)
    sd, wd = side(down)
    return (su - sd).dropna() - (wu - wd)


# Model gap change (seaward-positive sense: x5 - x6 in the landward-positive frame), mid-year
def model_change(path):
    x = np.load(path)
    g = x[:, DOWN] - x[:, UP]
    g = g - g[0]
    return 0.5 * (g[:-1] + g[1:])            # year k: mean of 1 Jan k and 1 Jan k+1


def score_annual(path, period, obs):
    from site_layer.hatteras_site_config import run_years
    m = model_change(path)
    yrs = [y for y in range(period, period + run_years(period)) if y in obs.index]
    mv = np.array([m[y - period] for y in yrs])
    ov = obs.loc[yrs].to_numpy()
    return dict(n_years=len(yrs), rmse_annual=float(np.sqrt(np.mean((mv - ov) ** 2))),
                bias_annual=float(np.mean(mv - ov)), end_model=float(mv[-1]), end_obs=float(ov[-1]))


def score_photos(path, period, photos):
    from site_layer.hatteras_site_config import run_years
    x = np.load(path)
    g = x[:, DOWN] - x[:, UP]
    g = g - g[0]
    n = run_years(period)
    yrs = [y for y in photos.index if period < y <= period + n]
    g0 = np.interp(period, photos.index, photos.values)
    o = np.array([photos[y] for y in yrs]) - g0
    m = np.array([g[y - period] for y in yrs])
    return dict(rmse_photos=float(np.sqrt(np.mean((m - o) ** 2))),
                photo_dates=" ".join(map(str, yrs)))


def cmd_score(_a):
    from score_instant_grid import observed_series
    photos = observed_series()
    rows, obs_rows = [], []
    for period in (CALIBRATION, TEST):
        obs = coastsat_change(period)
        obs_rows += [dict(period=period, year=int(y), coastsat_change_m=v) for y, v in obs.items()]
        base = baseline_file(period)
        if base.is_file():
            rows.append(dict(schedule="none", period=period, b=0.0, f=np.nan,
                             **score_annual(base, period, obs), **score_photos(base, period, photos)))
        for s in SCHEDULES:
            for b in B_VALUES:
                for f in F_VALUES:
                    p = shoreline_file(s, period, b, f)
                    if p is not None:
                        rows.append(dict(schedule=s, period=period, b=b, f=f,
                                         **score_annual(p, period, obs),
                                         **score_photos(p, period, photos)))
    t = pd.DataFrame(rows)
    t.to_csv(HERE / "grid_scores.csv", index=False, float_format="%.3f")
    pd.DataFrame(obs_rows).to_csv(HERE / "coastsat_gap_change.csv", index=False, float_format="%.2f")
    pd.set_option("display.width", 220)
    cal = t[t.period == CALIBRATION]
    print(cal[cal.schedule == "none"].to_string())
    for s, sub in cal[cal.schedule != "none"].groupby("schedule"):
        print(f"\n=== {s}: annual CoastSat RMSE (m) ===")
        print(sub.pivot(index="b", columns="f", values="rmse_annual").to_string(float_format="%.1f"))
        for r in sub.nsmallest(3, "rmse_annual").itertuples():
            print(f"  b {r.b:.2f} f {r.f:.1f}: annual {r.rmse_annual:.1f} (bias {r.bias_annual:+.1f}), "
                  f"photos {r.rmse_photos:.1f}")


def main():
    if len(sys.argv) > 2 and sys.argv[1] == "_launch":
        return _launch(sys.argv[2])
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    sub = ap.add_subparsers(dest="cmd", required=True)
    g = sub.add_parser("grid")
    g.add_argument("--streams", type=int, default=4)
    g.add_argument("--schedules", nargs="+", default=list(SCHEDULES))
    sub.add_parser("score")
    t = sub.add_parser("test")
    t.add_argument("--schedule", required=True, choices=list(SCHEDULES))
    t.add_argument("--b", type=float, required=True)
    t.add_argument("--f", type=float, required=True)
    a = ap.parse_args()
    if a.cmd == "grid":
        cells = [(s, CALIBRATION, b, f) for s in a.schedules for b in B_VALUES for f in F_VALUES]
        sys.exit(0 if run_all(cells, a.streams) else 1)
    if a.cmd == "test":
        sys.exit(0 if run_all([(a.schedule, TEST, a.b, a.f)], 1) else 1)
    cmd_score(a)


if __name__ == "__main__":
    main()
