#!/usr/bin/env python3
"""
Refit the corrected groin's strength for each window, run from the observed start shoreline, on today's runner.

    python window_refit.py run   [--streams 4]   # the strength grid in both windows, through the unchanged runner
    python window_refit.py score                 # -> window_scores.csv

The hindcast matrix starts each window from the observed 1996 or 2009 shoreline, which already holds the fillet
the groins built. The 1967-2025 fit (0.25, 2026-10-13-pregroin-planform/) describes the groin over its whole life;
started from an observed shoreline, 0.25 overshoots (10-12 window scores, 2026-10-09-corrected-groin/). This
redoes the 10-12 window fit on the runner as it is now (corrected groin built in; dune design 3.94 m and the
island-wide Dmaxel since 10-13), with the setup of that fit: edgeBE, full_management, no relocation, off 2017.
Scored on the gap across the groins (GIS 5|6) against annual CoastSat and, for 1996-2009, the photo dates, with
the scoring functions of corrected_groin.py.

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-10-14
"""

from __future__ import annotations

import argparse
import os
import re
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
sys.path.insert(0, str(HERE.parent / "2026-10-09-corrected-groin"))
import corrected_groin as cg  # noqa: E402  (scoring only; its run code needs commit 4d937b79)

# --- CONFIG ------------------------------------------------------------------
RUNNER = REPO / "scripts" / "hatteras_ms" / "HAT_hindcast_1984_2024.py"
STUDY = "groin/2026-10-14-window-refit"
RAW = REPO / "output" / "raw_runs" / "experiments" / STUDY
LOGS = HERE / "logs"
END_YEAR = 2017
# None = groin off (the no-groin run on the same setup)
STRENGTHS = {1996: (None, 0.10, 0.13, 0.16, 0.19, 0.22, 0.25),
             2009: (None, 0.04, 0.06, 0.08, 0.10, 0.12, 0.16, 0.25)}
WINDOW = {1996: "1996_2009", 2009: "2009_2025"}
WETDRY_CHANGE_TABLE = (REPO / "hard-structures" / "groin" / "1-observations" / "analysis"
                       / "1-aerial_photo_gap" / "Change_from_wetdry_1967_D2_D12.csv")
# -----------------------------------------------------------------------------


def member(s):
    return "nogroin" if s is None else f"s{s:.2f}"


def shoreline_file(period, s):
    hits = sorted((RAW / member(s) / WINDOW[period]).glob("*/*/*_shoreline_matrix.npy"))
    return hits[0] if hits else None


def run(period, s):
    env = {k: v for k, v in os.environ.items() if not k.startswith(("HAT_", "CORRECTED_GROIN_"))}
    env.update({
        "HAT_IGNORE_SETTINGS": "1", "HAT_START_YEAR": str(period),
        "HAT_SOURCE_SINK_PRESET": "edgeBE", "HAT_SCENARIO": "full_management", "HAT_RELOCATIONS": "false",
        "HAT_GROIN_ENABLED": "false" if s is None else "true",
        "HAT_GROIN_STRENGTH": "0" if s is None else str(s), "HAT_GROIN_END_YEAR": str(END_YEAR),
        "HAT_RUN_KIND": "experiment", "HAT_RUN_TAG": f"{STUDY}/{member(s)}",
        "HAT_SAVE_MODEL_STATE": "false", "HAT_SHOW_FIGURES": "false", "HAT_MAKE_GIFS": "false",
        "HAT_OVERWRITE": "false", "MPLBACKEND": "Agg", "PYTHONIOENCODING": "utf-8",
    })
    LOGS.mkdir(exist_ok=True)
    log = LOGS / f"{WINDOW[period]}_{member(s)}.log"
    t0 = time.time()
    with open(log, "w", encoding="utf-8") as fh:
        code = subprocess.call([sys.executable, str(RUNNER)], cwd=REPO, env=env,
                               stdout=fh, stderr=subprocess.STDOUT)
    print(f"{time.strftime('%H:%M:%S')}  {WINDOW[period]} {member(s)}: exit {code} "
          f"({(time.time() - t0) / 60:.1f} min)", flush=True)
    return code


def cmd_run(a):
    cells = [(p, s) for p, ss in STRENGTHS.items() for s in ss]
    todo = [c for c in cells if shoreline_file(*c) is None]
    print(f"{len(cells) - len(todo)} of {len(cells)} done; running {len(todo)} in {a.streams} streams", flush=True)
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

    threads = [threading.Thread(target=worker) for _ in range(a.streams)]
    for t in threads:
        t.start()
    for t in threads:
        t.join()
    print(f"DONE: {len(todo) - len(failed)}/{len(todo)} clean; failed {failed}", flush=True)


# The observed gap across the groins (GIS 5 minus GIS 6) at each photo date, from the 1967 wet/dry line;
# as score_instant_grid.observed_series, whose config module went with groin-sweep/ on 2026-10-14
def observed_series():
    tab = pd.read_csv(WETDRY_CHANGE_TABLE).set_index("Domain_ID")
    obs = {}
    for c in tab.columns:
        m = re.match(r"change_from_wetdry_1967_wetdry_(\d{4})", c)
        if m:
            obs.setdefault(int(m.group(1)), []).append(tab.loc[5, c] - tab.loc[6, c])
    return pd.Series({y: np.mean(v) for y, v in obs.items()}).sort_index()


def cmd_score(_a):
    photos = observed_series()
    rows = []
    for period, ss in STRENGTHS.items():
        obs = cg.observed_change(period)
        for s in ss:
            p = shoreline_file(period, s)
            if p is None:
                continue
            row = dict(period=period, strength=0.0 if s is None else s, groin=s is not None,
                       **cg.score_annual(p, period, obs), **cg.split_scores(p, period, obs))
            if period == cg.CALIBRATION:
                row.update(cg.score_photos(p, period, photos))
            rows.append(row)
    t = pd.DataFrame(rows)
    t.to_csv(HERE / "window_scores.csv", index=False, float_format="%.3f")
    pd.set_option("display.width", 200)
    print(t.drop(columns=["photo_dates"], errors="ignore").round(1).to_string(index=False))
    for period in STRENGTHS:
        g = t[(t.period == period) & t.groin]
        if len(g):
            best = g.loc[g.rmse_annual.idxmin()]
            print(f"{WINDOW[period]}: best annual at strength {best.strength:.2f} ({best.rmse_annual:.1f} m)")


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    sub = ap.add_subparsers(dest="cmd", required=True)
    r = sub.add_parser("run")
    r.add_argument("--streams", type=int, default=4)
    sub.add_parser("score")
    a = ap.parse_args()
    {"run": cmd_run, "score": cmd_score}[a.cmd](a)


if __name__ == "__main__":
    main()
