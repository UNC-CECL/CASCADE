"""
Does running the storms by the notebook's own rule (no split, no trim) change the model's skill?

    python scripts/hatteras_ms/experiments/HAT_storm_notebook_rule.py run      # run, then score
    python scripts/hatteras_ms/experiments/HAT_storm_notebook_rule.py score    # score the runs on disk
    python scripts/hatteras_ms/experiments/HAT_storm_notebook_rule.py figure   # storms by duration, both series

Two storm series, both windows (1996-2009 calibration, 2009-2025 test):
current  v3_split12_trim24, the series every run reads now
v3_240   historical_storm_creation_v3.ipynb unchanged: storms grouped when <24 h
         apart, not split, not trimmed, dropped only above 240 h (none are)
Managed edgeBE, full_management, groin off, adopted dunes (5.5 m, r 0.55-0.95).
Nothing in the main code changes: the period's storm file is swapped in the
run's own process before the unchanged runner executes.
Scores: overwash skill against the imagery (as the dune grid), interior
shoreline RMSE and bias, total overwash volume, and for 1996 the 2009 crest
against the lidar.

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-10-20
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
sys.path.insert(0, str(_HERE.parent))
import HAT_dune_ceiling_uniform as U  # noqa: E402
import HAT_dune_ceiling_dem_growth as DG  # noqa: E402
import HAT_dune_ceiling_growth_grid as G  # noqa: E402

# --- CONFIG ------------------------------------------------------------------
TAG = "storms-overwash/2026-10-20-storm-notebook-rule"
EXP_DIR = U.PROJECT_ROOT / "output" / "raw_runs" / "experiments" / TAG
VARIANTS = {"current": "v3_split12_trim24", "v3_240": "v3_240"}
JOBS = [(start, v) for start in U.PERIODS for v in VARIANTS]
WORKERS = 4


def storm_file(start, variant):
    from site_layer import hat_env_forcings as env
    return env.storm_series_file(*U.PERIODS[start][1], variant=VARIANTS[variant])


# Child process: point the period at the storm file, then run the unchanged hindcast
def _launch(start, variant):
    import runpy
    from site_layer import hatteras_site_config as sc
    from site_layer import hat_env_forcings as env
    sc.HATTERAS_PERIODS[start]["storm_file"] = env.init_relpath(storm_file(start, variant))
    print(f"EXPERIMENT storm file {sc.HATTERAS_PERIODS[start]['storm_file']}", flush=True)
    sys.path.insert(0, str(U.HINDCAST.parent))
    sys.argv = [str(U.HINDCAST)]
    runpy.run_path(str(U.HINDCAST), run_name="__main__")


def run_dir(start, variant):
    base = EXP_DIR / "runs" / variant / f"{start}_{U.PERIODS[start][0]}" / U.PRESET
    hits = sorted(d for d in base.glob("HAT_*") if any(d.glob("*_run_metadata.json"))) if base.exists() else []
    return hits[-1] if hits else None


def launch(job):
    start, variant = job
    logs = EXP_DIR / "logs"
    logs.mkdir(parents=True, exist_ok=True)
    env = {k: v for k, v in os.environ.items() if not k.startswith("HAT_")}
    env.update({"HAT_IGNORE_SETTINGS": "1", "HAT_START_YEAR": str(start),
                "HAT_SOURCE_SINK_PRESET": U.PRESET, "HAT_SCENARIO": U.SCENARIO,
                "HAT_GROIN_ENABLED": "False",
                "HAT_RUN_KIND": "experiment", "HAT_RUN_TAG": f"{TAG}/runs/{variant}",
                "HAT_SAVE_MODEL_STATE": "1", "HAT_MAKE_GIFS": "0", "HAT_SHOW_FIGURES": "0",
                "HAT_OVERWRITE": "1", "PYTHONIOENCODING": "utf-8", "MPLBACKEND": "Agg"})
    t0 = time.time()
    with open(logs / f"{variant}_{start}.log", "w", encoding="utf-8") as fh:
        p = subprocess.run([sys.executable, str(_HERE), "_launch", str(start), variant],
                           env=env, cwd=U.HINDCAST.parent, stdout=fh, stderr=subprocess.STDOUT)
    rec = dict(start=start, variant=variant, returncode=p.returncode, minutes=round((time.time() - t0) / 60, 1))
    with open(logs / "launches.jsonl", "a", encoding="utf-8") as fh:
        fh.write(json.dumps(rec) + "\n")
    print(f"{start} {variant:8s} exit {p.returncode}  {rec['minutes']} min", flush=True)


def metrics(start, variant):
    from site_layer import hat_overwash as ow
    from site_layer.hatteras_site_config import HATTERAS_DOMAINS as DOM
    rd = run_dir(start, variant)
    if rd is None:
        return None
    c = U.load(rd)
    obs = pd.read_csv(ow.OBSERVATIONS / "overwash_observations.csv")
    sk = DG.overwash_skill(c, obs, DOM, U.PERIODS[start][1])
    meta = json.loads((rd / f"{rd.name}_run_metadata.json").read_text(encoding="utf-8"))
    storms = np.load(storm_file(start, variant))
    qow = np.array([np.asarray(b.QowTS) for b in c.barrier3d])
    row = dict(period=f"{start}_{U.PERIODS[start][0]}", series=variant,
               events=len(storms), storm_hours=int(storms[:, 4].sum()),
               skill=round(sk["skill"], 3), hits=sk["hits"], misses=sk["misses"], false_alarms=sk["false_alarms"],
               rmse_interior_m_yr=float(meta["skill"]["rmse_interior_m_yr"]),
               bias_interior_m_yr=float(meta["skill"]["mean_bias_interior_m_yr"]),
               overwash_m3_per_m=round(float(qow.sum(axis=1).mean()), 1),
               dunes_rebuilt=int(pd.read_csv(rd / "tables" / "road_management.csv").dunes_rebuilt.sum()))
    if start == 1996:
        med_end = np.array([np.median(x[-1]) for x in U.cell_crests(c, DOM)])
        row.update(crest_minus_lidar_m=round(float(np.median(med_end - G.lidar())), 3),
                   crest_vs_lidar_r=round(float(np.corrcoef(med_end, G.lidar())[0, 1]), 3))
    return row


# Storms by duration, calibration against test, one panel per series
def figure():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from site_layer.hat_figure_style import (apply_style, C_1984, C_1997, figsize, open_frame, panel_title,
                                             record_caption, save)
    apply_style()
    labels = {"current": "Current series (split, trimmed to 24 h)", "v3_240": "Notebook rule (not split or trimmed)"}
    colours = {1996: C_1984, 2009: C_1997}
    dur = {(s, v): np.load(storm_file(s, v))[:, 4] for s in U.PERIODS for v in VARIANTS}
    edges = np.arange(0, 8 * np.ceil(max(d.max() for d in dur.values()) / 8) + 8, 8)
    fig, axes = plt.subplots(2, 1, figsize=figsize("double", height=4.6), sharex=True, sharey=True,
                             constrained_layout=True)
    w = 3.4
    for i, (ax, v) in enumerate(zip(axes, VARIANTS)):
        for k, s in enumerate(U.PERIODS):
            d = dur[(s, v)]
            # right-closed bins, so a 24 h storm counts in 16-24
            n = np.array([((d > a) & (d <= b)).sum() for a, b in zip(edges[:-1], edges[1:])])
            yrs = U.PERIODS[s][1]
            ax.bar(edges[:-1] + 4 + (k - 0.5) * w, n, width=w, color=colours[s], lw=0,
                   label=f"{'calibration' if s == 1996 else 'test'} {yrs[0]}-{yrs[1]}: {len(d)} storms, "
                         f"{int(d.sum()):,} h")
        ax.set_ylabel("number of storms")
        ax.text(0.99, 0.95, panel_title("ab"[i], labels[v]), transform=ax.transAxes, ha="right", va="top")
        ax.legend(frameon=False, fontsize=7.5, loc="center right")
        open_frame(ax)
    axes[-1].set_xlabel("storm duration the model receives (h, 8 h bins)")
    axes[-1].set_xticks(np.arange(0, edges[-1] + 1, 24))
    axes[-1].set_xlim(0, edges[-1])
    out = save(fig, EXP_DIR / "figures" / "storm_duration_histogram_current_vs_notebook_rule.png")
    plt.close(fig)
    record_caption(out[0],
        "Storms by the duration the model receives, calibration (red, 1997-2009) against test (blue, 2010-2025), "
        "in 8 h bins (a bin holds durations above its left edge up to and including its right edge). "
        "(a) The current series, v3_split12_trim24: storms less than 24 h apart are grouped, split where the "
        "water stays below the berm 12 h or more, and cut to the 24 h around the peak, so every storm longer "
        "than 24 h sits in the 16-24 h bin. (b) The notebook's rule, v3_240 "
        "(historical_storm_creation_v3.ipynb unchanged): grouped, not split, not trimmed; only storms over "
        "240 h would be dropped and none is. Duration counts hours with total water level above the 1.7 m "
        "NAVD88 berm. Both series from the Duck and WIS records extended through 2025.")
    print("\n".join(str(p) for p in out))


def main():
    if len(sys.argv) > 1 and sys.argv[1] == "_launch":
        _launch(int(sys.argv[2]), sys.argv[3])
        return
    if len(sys.argv) > 1 and sys.argv[1] == "figure":
        figure()
        return
    if len(sys.argv) > 1 and sys.argv[1] == "run":
        with ThreadPoolExecutor(WORKERS) as ex:
            list(ex.map(launch, [j for j in JOBS if run_dir(*j) is None]))
    df = pd.DataFrame([m for m in (metrics(*j) for j in JOBS) if m])
    (EXP_DIR / "tables").mkdir(parents=True, exist_ok=True)
    df.to_csv(EXP_DIR / "tables" / "scores.csv", index=False)
    with pd.option_context("display.width", 250, "display.max_columns", 30):
        print(df.to_string(index=False), flush=True)


if __name__ == "__main__":
    main()
