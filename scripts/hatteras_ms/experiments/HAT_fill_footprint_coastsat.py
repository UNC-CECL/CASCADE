"""
Does placing each fill on its CoastSat-observed footprint improve the 2009-2025 test run?

    python scripts/hatteras_ms/experiments/HAT_fill_footprint_coastsat.py run
    python scripts/hatteras_ms/experiments/HAT_fill_footprint_coastsat.py score

Runs the test period (domainBE, blocking groin, full management, relocations off) twice:
with the reported footprints as committed, and with each fill's domains swapped in-process
for the observed range in 4-extent-checks/coastsat/nourishment_extent_coastsat_summary.csv. Volumes
are unchanged, so the narrower footprints carry more sand per metre. Scored on net change
in metres against the test target, island-wide (both sides 7-domain LOWESS) and per fill
footprint (raw domain values). Details: the study NOTE.md.

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-10-06
"""
from __future__ import annotations

import argparse
import dataclasses
import json
import os
import subprocess
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

_HERE = Path(__file__).resolve()
PROJECT_ROOT = next(_p for _p in _HERE.parents if (_p / "pyproject.toml").exists())
sys.path[:0] = [str(PROJECT_ROOT / "scripts"),
                str(PROJECT_ROOT / "scripts" / "input_prep" / "5-scr" / "3-rates" / "coastsat" / "net_change")]

# --- CONFIG ------------------------------------------------------------------
TAG = "management/2026-10-06-fill-footprint-coastsat"
EXP_DIR = PROJECT_ROOT / "output" / "raw_runs" / "experiments" / TAG
HINDCAST = PROJECT_ROOT / "scripts" / "hatteras_ms" / "HAT_hindcast_1984_2024.py"
MATRIX_RUN = (PROJECT_ROOT / "output" / "raw_runs" / "matrix" / "2009_2025" / "domainBE"
              / "HAT_2009_2025_domainBE_offsetmetres_road_bdm_nourish_groinblock")
OBSERVED_EXTENT = (PROJECT_ROOT / "data" / "hatteras_init" / "4-mgmt-forcing" / "nourishment"
                   / "4-extent-checks" / "coastsat" / "nourishment_extent_coastsat_summary.csv")
WINDOW = (2009, 2025)
VARIANTS = ("reported", "coastsat")
# -----------------------------------------------------------------------------


def observed_ranges():
    t = pd.read_csv(OBSERVED_EXTENT)
    return {(r.project, int(r.model_year)): tuple(range(int(r.observed_first_gis), int(r.observed_last_gis) + 1))
            for r in t.itertuples()}


# Child process: swap the footprints if asked, then run the unchanged hindcast
def _launch(variant):
    import runpy
    sys.path.insert(0, str(HINDCAST.parent))
    from site_layer import hatteras_site_config as sc
    if variant == "coastsat":
        ranges = observed_ranges()
        assert set(ranges) == {(p.name, p.year) for p in sc.HATTERAS_NOURISHMENT_PROJECTS}
        sc.HATTERAS_NOURISHMENT_PROJECTS = tuple(
            dataclasses.replace(p, gis_domains=ranges[(p.name, p.year)]) for p in sc.HATTERAS_NOURISHMENT_PROJECTS)
    sys.argv = [str(HINDCAST)]
    runpy.run_path(str(HINDCAST), run_name="__main__")


def run_dir(v):
    hits = sorted((EXP_DIR / "runs" / v).glob("*/*/*/*_shoreline_matrix.npy"))
    return hits[-1].parent if hits else None


def launch(v):
    logs = EXP_DIR / "logs"
    logs.mkdir(parents=True, exist_ok=True)
    env = {k: val for k, val in os.environ.items() if not k.startswith("HAT_")}
    env.update(PYTHONIOENCODING="utf-8", PYTHONUNBUFFERED="1", MPLBACKEND="Agg",
               HAT_IGNORE_SETTINGS="1", HAT_START_YEAR=str(WINDOW[0]),
               HAT_SOURCE_SINK_PRESET="domainBE", HAT_SCENARIO="full_management", HAT_RELOCATIONS="False",
               HAT_GROIN_ENABLED="True", HAT_OVERWRITE="False", HAT_MAKE_GIFS="False",
               HAT_SHOW_FIGURES="False", HAT_RUN_KIND="experiment", HAT_RUN_TAG=f"{TAG}/runs/{v}",
               HAT_SAVE_MODEL_STATE="False")
    t0 = time.time()
    with open(logs / f"{v}.log", "w", encoding="utf-8") as fh:
        p = subprocess.run([sys.executable, str(_HERE), "_launch", v], env=env, cwd=PROJECT_ROOT,
                           stdout=fh, stderr=subprocess.STDOUT)
    rec = dict(variant=v, returncode=p.returncode, minutes=round((time.time() - t0) / 60, 1))
    with open(logs / "launches.jsonl", "a", encoding="utf-8") as fh:
        fh.write(json.dumps(rec) + "\n")
    print(f"{v:9s} exit {p.returncode} {rec['minutes']} min", flush=True)
    if p.returncode:
        raise SystemExit(f"{v} failed; see {logs / f'{v}.log'}")


def model_net(rd):
    from site_layer.hatteras_site_config import HATTERAS_DOMAINS as D
    m = np.load(next(Path(rd).glob("*_shoreline_matrix.npy")))
    return pd.Series(-(m[-1] - m[0])[D.start_real_index:D.end_real_index], index=pd.RangeIndex(1, 91))


def score():
    from coastsat_net_change import smooth_like_model
    from site_layer.hat_observed_rates import net_change_domain_csv
    from site_layer.hatteras_site_config import HATTERAS_NOURISHMENT_PROJECTS
    obs = pd.read_csv(net_change_domain_csv(*WINDOW), index_col=0)
    runs = {v: run_dir(v) for v in VARIANTS}
    missing = [v for v, rd in runs.items() if rd is None]
    if missing:
        raise SystemExit(f"missing runs: {missing}")
    net = {v: model_net(rd) for v, rd in runs.items()}
    d = float((net["reported"] - model_net(MATRIX_RUN)).abs().max())
    print(f"reported vs matrix run: max |net change difference| {d:.2e} m")

    rows = []
    for v in VARIANTS:
        m, o = smooth_like_model(net[v]).loc[2:89], obs["net_change_lowess7_m"].loc[2:89]
        rows.append(dict(variant=v, area="interior GIS 2-89, smoothed", bias_m=(m - o).mean(),
                         rmse_m=float(np.sqrt(((m - o) ** 2).mean())), r=float(np.corrcoef(m, o)[0, 1])))
    observed = observed_ranges()
    for p in HATTERAS_NOURISHMENT_PROJECTS:
        span = sorted(set(p.gis_domains) | set(observed[(p.name, p.year)]))
        o = obs["net_change_m"].loc[span]
        for v in VARIANTS:
            dm = net[v].loc[span] - o
            rows.append(dict(variant=v, area=f"{p.name} {p.year}, GIS {span[0]}-{span[-1]} raw",
                             bias_m=dm.mean(), rmse_m=float(np.sqrt((dm ** 2).mean())), r=np.nan))
    s = pd.DataFrame(rows)
    per = pd.DataFrame({"observed_m": obs["net_change_m"], "observed_lowess7_m": obs["net_change_lowess7_m"],
                        **{f"model_{v}_m": net[v] for v in VARIANTS}})
    per["coastsat_minus_reported_m"] = per["model_coastsat_m"] - per["model_reported_m"]
    (EXP_DIR / "tables").mkdir(parents=True, exist_ok=True)
    s.round(3).to_csv(EXP_DIR / "tables" / "scores.csv", index=False)
    per.round(3).to_csv(EXP_DIR / "tables" / "per_domain_net_change.csv")
    with pd.option_context("display.width", 200, "display.max_columns", 20):
        print(s.round(2).to_string(index=False))


def main():
    if len(sys.argv) > 1 and sys.argv[1] == "_launch":
        _launch(sys.argv[2])
        return
    ap = argparse.ArgumentParser()
    ap.add_argument("action", choices=["run", "score"])
    a = ap.parse_args()
    if a.action == "run":
        for v in VARIANTS:
            if run_dir(v) is None:
                launch(v)
    else:
        score()


if __name__ == "__main__":
    main()
