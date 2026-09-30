"""
Do a Hatteras dune ceiling and a rebuild that never lowers dunes fix the excess overwash?

    python scripts/hatteras_ms/experiments/HAT_dune_ceiling_rebuild.py run --workers 6
    python scripts/hatteras_ms/experiments/HAT_dune_ceiling_rebuild.py score

Dmaxel 3.4-9.0 m against three rebuild rules, both storm series and
windows, set in each run's own process. Details: scripts/hatteras_ms/experiments/README.md.

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-09-29
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
sys.path.insert(0, str(_HERE.parent))
import HAT_storm_length_selection as S  # noqa: E402

# --- CONFIG ------------------------------------------------------------------
TAG = "storms-and-overwash/2026-09-28-dune-ceiling-and-rebuild"
EXP_DIR = PROJECT_ROOT / "output" / "raw_runs" / "experiments" / TAG
DMAX = (3.4, 5.5, 7.5, 9.0)
RULES = ("current", "nolower", "nolower43")
STORMS = ("drop72", "trim24")
BERM_MHW_M = 1.34
DESIGN43_MHW_M = 4.3 - 0.36
# -----------------------------------------------------------------------------


# A Dmaxel's name token
def dtok(d):
    return f"dmax{str(d).replace('.', 'p')}"


# A job's group name
def group(storm, dmax, rule, scenario):
    return f"{storm}_{dtok(dmax)}_{rule}_{scenario}"


# Every job not already covered by a reused run
def jobs():
    out = []
    for w in S.WINDOWS:
        for st in STORMS:
            for d in DMAX:
                for r in RULES:
                    if not (d == 3.4 and r == "current"):
                        out.append((w, st, d, r, "full_management"))
                if d != 3.4:
                    out.append((w, st, d, "current", "natural"))
    return out


# The reused run for the current Dmaxel and rule, or None
def existing(w, st, d, r, scen):
    if d == 3.4 and r == "current":
        if st == "drop72":
            return S.MATRIX / S.wtag(w) / "edgeBE" / S.CONTROLS[(w[0], scen)]
        return S.run_path("trim24", scen, w)
    return None


# A job's run folder (reused or new), or None
def run_dir(w, st, d, r, scen):
    e = existing(w, st, d, r, scen)
    if e is not None:
        return e
    base = EXP_DIR / "runs" / group(st, d, r, scen) / S.wtag(w) / "edgeBE"
    hits = sorted(base.glob("HAT_*")) if base.exists() else []
    return hits[-1] if hits else None


# The in-process changes

# Child process: set Dmaxel and wrap rebuild_dunes, then run the hindcast
def _launch(start, storm_path, dmax, rule):
    import cascade.brie_coupler as bc
    import cascade.roadway_manager as rm
    import cascade.beach_dune_manager as bdm

    dmax = float(dmax)
    orig_set = bc.set_yaml

    def set_yaml(var_name, new_vals, file_name):
        orig_set(var_name, new_vals, file_name)
        if var_name != "Dmaxel":
            orig_set("Dmaxel", dmax, file_name)          # m NAVD88, as load_input expects

    bc.set_yaml = set_yaml

    if rule != "current":
        orig_rebuild = rm.rebuild_dunes

        def rebuild_dunes(yxz_dune_grid, max_dune_height=3.0, min_dune_height=2.4, dz=10, rng=True):
            if rule == "nolower43":
                max_dune_height = min_dune_height = DESIGN43_MHW_M - BERM_MHW_M
            new, _ = orig_rebuild(yxz_dune_grid, max_dune_height, min_dune_height, dz, rng)
            new = np.maximum(new, yxz_dune_grid)              # never lower a dune cell
            return new, float(np.sum(new - yxz_dune_grid))

        rm.rebuild_dunes = rebuild_dunes
        bdm.rebuild_dunes = rebuild_dunes

    sys.argv = [sys.argv[0]]
    S.MD._launch(start, storm_path)


# One run in its own process, logged, with its launch record
def launch(job):
    w, st, d, r, scen = job
    g = group(st, d, r, scen)
    logs = EXP_DIR / "logs"
    logs.mkdir(parents=True, exist_ok=True)
    log = logs / f"{g}_{S.wtag(w)}.log"
    env = dict(os.environ)
    env.update({"HAT_IGNORE_SETTINGS": "1", "HAT_START_YEAR": str(w[0]),
                "HAT_SOURCE_SINK_PRESET": "edgeBE", "HAT_SCENARIO": scen,
                "HAT_RUN_KIND": "experiment", "HAT_RUN_TAG": f"{TAG}/runs/{g}",
                "HAT_SAVE_MODEL_STATE": "1", "HAT_MAKE_GIFS": "0", "HAT_SHOW_FIGURES": "0",
                "HAT_OVERWRITE": "1", "PYTHONIOENCODING": "utf-8"})
    t0 = time.time()
    with open(log, "w", encoding="utf-8") as fh:
        p = subprocess.run([sys.executable, str(_HERE), "_launch", str(w[0]), str(S.storm_file(w, st)), str(d), r],
                           env=env, cwd=S.MD.HINDCAST.parent, stdout=fh, stderr=subprocess.STDOUT)
    rec = dict(window=S.wtag(w), storm=st, dmaxel=d, rule=r, scenario=scen, returncode=p.returncode,
               minutes=round((time.time() - t0) / 60, 1))
    with open(logs / "launches.jsonl", "a", encoding="utf-8") as fh:
        fh.write(json.dumps(rec) + "\n")
    print(f"{S.wtag(w)} {g:44s} exit {p.returncode} {rec['minutes']} min", flush=True)
    return rec


# Launch every job that has no run folder yet, in parallel
def run(workers):
    with ThreadPoolExecutor(workers) as ex:
        list(ex.map(launch, [j for j in jobs() if run_dir(*j) is None]))


# Scores

# Median dune crest per domain per year, m
def crest_series(c, pads):
    return np.array([[np.median((np.asarray(c.barrier3d[p].DuneDomain[t]).max(axis=1)
                                 + c.barrier3d[p].BermEl) * 10) for p in pads]
                     for t in range(len(c.barrier3d[0].x_s_TS))])


# Overwash, dune crest and shoreline scores for every job
def score():
    from site_layer.hatteras_site_config import HATTERAS_DOMAINS as DOM
    from site_layer import hat_overwash as ow
    ovm = S.overwash_module()
    obs = pd.read_csv(ow.OBSERVATIONS / "overwash_observations.csv")
    pads = [DOM.gis_to_pad(g) for g in range(1, 91)]
    lidar_c = S.load_state(S.MATRIX / "2010_2024/edgeBE" / S.CONTROLS[(2010, "full_management")])
    lidar = crest_series(lidar_c, pads)[0]                       # the 2009 lidar crests
    cells_all = []
    everything = [(w, st, d, r, sc) for w in S.WINDOWS for st in STORMS for d in DMAX
                  for r in RULES for sc in ("full_management", "natural")
                  if sc == "full_management" or r == "current"]
    rows = []
    for w, st, d, r, sc in everything:
        rd = run_dir(w, st, d, r, sc)
        if rd is None:
            print(f"  missing {S.wtag(w)} {group(st, d, r, sc)}")
            continue
        c = S.load_state(rd)
        b0 = c.barrier3d[pads[0]]
        got_dmax = round(float(b0._Dmaxel * 10 + b0._MHW * 10), 2) if hasattr(b0, "_MHW") else np.nan
        summ = pd.read_csv(S.summary_csv(w, st))
        cells = S.overwash_cells(c, w, summ, obs, ovm)
        cells_all.append(cells.assign(window=S.wtag(w), storm=st, dmaxel=d, rule=r, scenario=sc))
        cs = crest_series(c, pads)
        row = dict(window=S.wtag(w), storm=st, dmaxel=d, rule=r, scenario=sc, dmaxel_loaded_m=got_dmax,
                   crest_start_m=float(np.median(cs[0])), crest_end_m=float(np.median(cs[-1])),
                   crest_min_year_m=float(np.median(cs, axis=1).min()),
                   **S.overwash_scores(cells, 0.0), **S.shoreline_scores(rd))
        if w[0] == 1996:
            dl = cs[-1] - lidar
            row.update(crest_2010_minus_lidar_m=float(np.median(dl)), domains_lower_1m=int((dl < -1).sum()))
        rows.append(row)
        print(f"  scored {S.wtag(w)} {group(st, d, r, sc)}", flush=True)
    t = pd.DataFrame(rows)
    (EXP_DIR / "tables").mkdir(parents=True, exist_ok=True)
    t.to_csv(EXP_DIR / "tables" / "scores.csv", index=False)
    pd.concat(cells_all).to_csv(EXP_DIR / "tables" / "cells.csv", index=False)
    with pd.option_context("display.width", 260, "display.max_columns", 30):
        print(t[["window", "scenario", "storm", "dmaxel", "rule", "dmaxel_loaded_m", "crest_end_m",
                 "crest_2010_minus_lidar_m", "POD", "POFD", "PSS", "timing_r", "space_r",
                 "rmse_interior_m_yr", "bias_interior_m_yr"]].round(2).to_string(index=False))


# Run: the action asked for; `_launch` is the child process's entry point
def main():
    if len(sys.argv) > 1 and sys.argv[1] == "_launch":
        _launch(*sys.argv[2:6])
        return
    ap = argparse.ArgumentParser()
    ap.add_argument("action", choices=["run", "score"])
    ap.add_argument("--workers", type=int, default=6)
    a = ap.parse_args()
    run(a.workers) if a.action == "run" else score()


if __name__ == "__main__":
    main()
