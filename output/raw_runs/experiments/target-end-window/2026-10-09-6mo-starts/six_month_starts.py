#!/usr/bin/env python3
"""
Does a +/-6-month start shoreline change the runs, compared with the +/-1-yr one?

    python six_month_starts.py run --arm matrix   ->  the 6 production runs on the current config (+/-1 yr)
    python six_month_starts.py run --arm 6mo      ->  the same 6 runs on shoreline offset v3 (+/-6 months)
    python six_month_starts.py score              ->  scores.csv, domain_net_change.csv

Two arms, three presets, both periods, the unchanged runner. The 6mo arm pins
HAT_OFFSET_VERSION_<year>_SHORELINE=v3 and files as an experiment. Each arm is
scored against the target that shares its start line: the 6mo arm's target
swaps the start window for the +/-6-month one and keeps the end window.

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-10-09
"""
from __future__ import annotations

import argparse
import os
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
REPO = next(p for p in HERE.parents if (p / "pyproject.toml").exists())
NC = REPO / "scripts" / "input_prep" / "5-scr" / "3-rates" / "coastsat" / "net_change"
sys.path[:0] = [str(REPO / "scripts"), str(NC)]
from coastsat_net_change import domain_net_change, smooth_like_model, window_means  # noqa: E402
from site_layer.hat_observed_rates import NET_CHANGE_WINDOWS  # noqa: E402
from site_layer.hatteras_site_config import HATTERAS_DOMAINS as D  # noqa: E402

# --- CONFIG ------------------------------------------------------------------
HINDCAST = REPO / "scripts" / "hatteras_ms" / "HAT_hindcast_1984_2024.py"
RAW = REPO / "output" / "raw_runs"
TAG = "target-end-window/2026-10-09-6mo-starts"
LOGS = HERE / "logs"
WORKERS = 4
GIS = pd.RangeIndex(1, 91)
INTERIOR = slice(2, 89)
PERIODS = {1996: (1996, 2009), 2009: (2009, 2025)}
# The +/-6-month start lines (offset v3) and the +/-1-yr ones the config uses (offset v2)
START_6MO = {1996: ("1996-04-12", "1997-04-12"), 2009: ("2009-02-17", "2010-02-17")}
# (preset, groin) per period; the fill token follows the period
RUNS = [("zeroBE", False), ("edgeBE", False), ("domainBE", True)]
# -----------------------------------------------------------------------------


def run_name(start, preset, groin):
    fill = "_nourish" if start == 2009 else ""
    return (f"HAT_{start}_{PERIODS[start][1]}_{preset}_offsetmetres_road_bdm{fill}_"
            f"{'groinblock' if groin else 'nogroin'}")


def run_dir(arm, start, preset, groin):
    root = RAW / "matrix" if arm == "matrix" else RAW / "experiments" / TAG
    return root / f"{start}_{PERIODS[start][1]}" / preset / run_name(start, preset, groin)


def launch(arm, start, preset, groin):
    env = dict(os.environ, PYTHONIOENCODING="utf-8", PYTHONUNBUFFERED="1", MPLBACKEND="Agg",
               HAT_IGNORE_SETTINGS="1", HAT_START_YEAR=str(start), HAT_SOURCE_SINK_PRESET=preset,
               HAT_SCENARIO="full_management", HAT_RELOCATIONS="False",
               HAT_GROIN_ENABLED=str(groin), HAT_OVERWRITE="False", HAT_MAKE_GIFS="False",
               HAT_SHOW_FIGURES="False")
    if arm == "matrix":
        env["HAT_RUN_KIND"] = "matrix"
    else:
        env.update(HAT_RUN_KIND="experiment", HAT_RUN_TAG=TAG, HAT_SAVE_MODEL_STATE="False",
                   **{f"HAT_OFFSET_VERSION_{start}_SHORELINE": "v3"})
    LOGS.mkdir(exist_ok=True)
    log = LOGS / f"{arm}_{start}_{preset}.log"
    with open(log, "w", encoding="utf-8") as f:
        code = subprocess.run([sys.executable, str(HINDCAST)], stdout=f, stderr=subprocess.STDOUT,
                              env=env, cwd=REPO).returncode
    print(f"{arm:6s} {start} {preset:8s} exit {code}  ({log.name})", flush=True)
    return code


def model_change(rd):
    m = np.load(next(Path(rd).glob("*_shoreline_matrix.npy")))
    return pd.Series(-(m[-1] - m[0])[D.start_real_index:D.end_real_index], index=GIS)


def target(start_w, end_w):
    s, e = window_means(start_w), window_means(end_w)
    both = s.index.intersection(e.index)
    tr = pd.DataFrame({"domain_number": s.loc[both, "domain_number"].astype(int),
                       "net_change_m": e.loc[both, "mean_chainage_m"] - s.loc[both, "mean_chainage_m"],
                       "se_net_change_m": np.hypot(s.loc[both, "se_chainage_m"],
                                                   e.loc[both, "se_chainage_m"])})
    return domain_net_change(tr)["net_change_lowess7_m"]


def score(model, obs):
    m, o = smooth_like_model(model).loc[INTERIOR], obs.loc[INTERIOR]
    d = m - o
    return dict(bias_m=d.mean(), rmse_m=float(np.sqrt((d ** 2).mean())),
                r=float(np.corrcoef(m, o)[0, 1]), obs_mean_m=o.mean(), model_mean_m=m.mean())


def cmd_run(arm):
    jobs = [(arm, s, p, g) for s in PERIODS for p, g in RUNS]
    with ThreadPoolExecutor(WORKERS) as ex:
        codes = list(ex.map(lambda j: launch(*j), jobs))
    sys.exit(max(codes))


def cmd_score():
    out, rows = pd.DataFrame(index=GIS), []
    for start, period in PERIODS.items():
        start_1yr, end_w = NET_CHANGE_WINDOWS[period]
        obs = {"matrix": target(start_1yr, end_w), "6mo": target(START_6MO[start], end_w)}
        for arm, o in obs.items():
            out[f"{start}_target_{arm}_lowess7_m"] = o
        for preset, groin in RUNS:
            mod = {arm: model_change(run_dir(arm, start, preset, groin)) for arm in obs}
            for arm in obs:
                out[f"{start}_{preset}_{arm}_lowess7_m"] = smooth_like_model(mod[arm])
                rows.append(dict(period=f"{period[0]}_{period[1]}", run=preset, arm=arm,
                                 start_window=" to ".join(START_6MO[start] if arm == "6mo" else start_1yr),
                                 **score(mod[arm], obs[arm])))
            d = (smooth_like_model(mod["6mo"]) - smooth_like_model(mod["matrix"])).loc[INTERIOR]
            print(f"{start} {preset:8s} model change, 6mo minus matrix (interior LOWESS-7): "
                  f"mean {d.mean():+.2f} m, RMS {np.sqrt((d ** 2).mean()):.2f} m")
        t = (obs["6mo"] - obs["matrix"]).loc[INTERIOR]
        print(f"{start} target, 6mo minus matrix: mean {t.mean():+.2f} m, RMS {np.sqrt((t ** 2).mean()):.2f} m")
    out.round(3).to_csv(HERE / "domain_net_change.csv", index_label="gis_domain")
    t = pd.DataFrame(rows).round(2)
    t.to_csv(HERE / "scores.csv", index=False)
    pd.set_option("display.width", 200)
    print("\n" + t.to_string(index=False))


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    sub = ap.add_subparsers(dest="cmd", required=True)
    r = sub.add_parser("run")
    r.add_argument("--arm", choices=["matrix", "6mo"], required=True)
    sub.add_parser("score")
    a = ap.parse_args()
    cmd_run(a.arm) if a.cmd == "run" else cmd_score()


if __name__ == "__main__":
    main()
