#!/usr/bin/env python3
"""
Refit the blocking groin (b, f) on 1996-2009 with a conserving face exchange, and compare with the 10-05 pin.

    python conserving_groin.py grid  [--streams 4]        # the calibration grid
    python conserving_groin.py score                      # -> grid_scores.csv, budget.csv
    python conserving_groin.py test --b 0.6 --f 0.6       # one cell on 2009-2025

The pinned BlockingGroinCallback cancels b of each cell's own coupling to the face, with
each cell's own diffusion number, so GIS 5 loses about four times what GIS 6 gains. Here
both sides use one coefficient, r_face: BRIE's diffusion number for the shoreline angle
across the GIS 5|6 face (the forward difference brie.py assigns to the lower cell), so
dx_up = -dx_down every year. The unchanged runner is driven in a child process with
_r_ipl patched; everything else is the 10-05 setup: full management, the solved edgeBE
ends, relocations off, instant failure at the 2004 step. Scores: the 10-05 photo-date
RMSE and the 10-08 annual CoastSat RMSE of the GIS 5|6 gap change.

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-10-08
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
GROIN = REPO / "hard-structures" / "groin"
sys.path[:0] = [str(HERE.parent / "2026-10-08-schedule-refit"),
                str(GROIN / "1-observations" / "gap_across_groins"),
                str(GROIN / "2-module-tests" / "3-real-planform"),
                str(REPO / "scripts" / "hatteras_ms" / "groin-sweep"),
                str(REPO / "scripts" / "hatteras_ms"), str(REPO / "scripts")]

# --- CONFIG ------------------------------------------------------------------
RUNNER = REPO / "scripts" / "hatteras_ms" / "HAT_hindcast_1984_2024.py"
STUDY = "groin/2026-10-08-conserving-groin"
PIN_STUDY = "groin/2026-10-05-blocking-fit-dem-to-dem"
RAW = REPO / "output" / "raw_runs"
LOGS = HERE / "logs"
B_VALUES = (0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.4, 0.5, 0.6)   # b 0.6 overshoots 2004 by ~100 m
F_VALUES = (0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.8)
CALIBRATION, TEST = 1996, 2009
UP, DOWN = 15 + 6 - 1, 15 + 5 - 1      # GIS 6 / GIS 5 in the padded array
NEAR = {"GIS 3": 15 + 3 - 1, "GIS 4": 15 + 4 - 1, "GIS 7": 15 + 7 - 1, "GIS 8": 15 + 8 - 1}
# -----------------------------------------------------------------------------


def member(b, f):
    return f"b{b:.2f}_f{f:.1f}"


def window(period):
    from site_layer.hatteras_site_config import HATTERAS_PERIODS
    return f"{period}_{HATTERAS_PERIODS[period]['end_year']}"


def run_dir(study, period, b, f):
    hits = sorted((RAW / "experiments" / study / member(b, f)).glob(f"{window(period)}/edgeBE/*/"))
    return hits[-1] if hits else None


def shoreline_file(study, period, b, f):
    d = run_dir(study, period, b, f)
    hits = sorted(d.glob("*_shoreline_matrix.npy")) if d else []
    return hits[-1] if hits else None


# Child process: one face coefficient for both sides, then the unchanged runner
def _launch():
    import runpy
    import cascade.groin as cg
    pinned_r = cg.BlockingGroinCallback._r_ipl

    class Conserving(cg.BlockingGroinCallback):
        kind = "blocking"

        def _r_ipl(self, brie, x_s, i):
            return pinned_r(brie, x_s, self._lo)

    cg.BlockingGroinCallback = Conserving
    sys.path.insert(0, str(RUNNER.parent))
    sys.argv = [str(RUNNER)]
    runpy.run_path(str(RUNNER), run_name="__main__")


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
        code = subprocess.call([sys.executable, str(Path(__file__)), "_launch"], cwd=REPO,
                               env=env, stdout=fh, stderr=subprocess.STDOUT)
    d = run_dir(STUDY, period, b, f)
    if code == 0 and d is not None:
        (d / "conserving.txt").write_text(
            "patched in-process by conserving_groin.py: both sides of the face use r_face "
            "(BRIE's diffusion number of the lower cell); the runner's report does not show it\n")
    print(f"{time.strftime('%H:%M:%S')}  {window(period)} {member(b, f)}: exit {code} "
          f"({(time.time() - t0) / 60:.1f} min)", flush=True)
    return code


def run_all(cells, streams):
    todo = [c for c in cells if shoreline_file(STUDY, *c) is None]
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


# The groin's applied shoreline changes, from the run's own diagnostics
def budget(path):
    d = pd.read_csv(path.parent / "tables" / "groin_diagnostics.csv")
    a = d[d.groin_active]
    up, down = a.applied_dx_updrift_m.sum(), a.applied_dx_downdrift_m.sum()
    return dict(applied_up_m=up, applied_down_m=down, applied_net_m=up + down)


# Landward-positive x_s, so seaward change is start minus end
def near_change(path):
    x = np.load(path)
    return {f"{k}_seaward_m": float(x[0, i] - x[-1, i]) for k, i in NEAR.items()}


def cmd_score(_a):
    from schedule_refit import baseline_file, coastsat_change, score_annual, score_photos
    from score_instant_grid import observed_series
    photos = observed_series()
    rows = []
    for period in (CALIBRATION, TEST):
        obs = coastsat_change(period)
        base = baseline_file(period)
        if base.is_file():
            rows.append(dict(arm="none", period=period, b=0.0, f=np.nan,
                             **score_annual(base, period, obs), **score_photos(base, period, photos),
                             **near_change(base)))
        for arm, study in (("conserving", STUDY), ("pinned", PIN_STUDY)):
            for b in B_VALUES:
                for f in F_VALUES:
                    p = shoreline_file(study, period, b, f)
                    if p is not None:
                        rows.append(dict(arm=arm, period=period, b=b, f=f,
                                         **score_annual(p, period, obs),
                                         **score_photos(p, period, photos),
                                         **budget(p), **near_change(p)))
    t = pd.DataFrame(rows)
    t.to_csv(HERE / "grid_scores.csv", index=False, float_format="%.3f")
    pd.set_option("display.width", 220)
    cal = t[(t.period == CALIBRATION) & (t.arm == "conserving")]
    for score in ("rmse_photos", "rmse_annual"):
        print(f"\n=== conserving, calibration: {score} (m) ===")
        print(cal.pivot(index="b", columns="f", values=score).to_string(float_format="%.1f"))
    print(f"\nlargest |applied net| in the conserving arm: {cal.applied_net_m.abs().max():.2e} m")


def main():
    if len(sys.argv) > 1 and sys.argv[1] == "_launch":
        return _launch()
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
        cells = [(CALIBRATION, b, f) for b in B_VALUES for f in F_VALUES]
        sys.exit(0 if run_all(cells, a.streams) else 1)
    if a.cmd == "test":
        sys.exit(0 if run_all([(TEST, a.b, a.f)], 1) else 1)
    cmd_score(a)


if __name__ == "__main__":
    main()
