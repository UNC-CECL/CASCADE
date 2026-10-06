"""
Solve a source/sink rate for every domain so the calibration run's net change matches CoastSat.

    python scripts/hatteras_ms/HAT_be_domain_solve_net_change.py solve --period 1996
    python scripts/hatteras_ms/HAT_be_domain_solve_net_change.py report --period 1996

Starts from the pinned calibration run (solved end rates, blocking groin, 0 elsewhere).
Each pass compares the run's raw end-minus-start shoreline with the net-change target
(7-domain LOWESS, GIS 1-10 raw) and moves each domain's rate by RELAX x residual / gain,
where the gain is metres of change per m/yr over the run. Stops when the interior
(GIS 2-89) RMSE is under 0.5 m and no domain misses by more than 1 m, or after
MAX_STEPS passes. Full management, groin on, relocations off. Probes are experiment
runs under source-sink/<date>-be-domain-solve-<window>/; the field per pass is written
beside them.

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-10-05
"""
from __future__ import annotations

import argparse
import os
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd

_HERE = Path(__file__).resolve()
PROJECT_ROOT = next(p for p in _HERE.parents if (p / "pyproject.toml").exists())
sys.path.insert(0, str(PROJECT_ROOT / "scripts"))

from site_layer.hat_observed_rates import net_change_domain_csv  # noqa: E402
from site_layer.hatteras_site_config import (  # noqa: E402
    HATTERAS_BE_RATES_EDGE, HATTERAS_DOMAINS, HATTERAS_PERIODS, run_years)

# --- CONFIG ------------------------------------------------------------------
HINDCAST = _HERE.parent / "HAT_hindcast_1984_2024.py"
SCENARIO = "full_management"
SOLVE_DATE = "2026-10-05"
MAX_STEPS = 10
RELAX = 0.8
RMSE_TOL_M, MAX_TOL_M = 0.5, 1.0
GAIN_OVERRIDE = {90: 4.4}          # m per m/yr measured in the end solve; elsewhere the run length
RAW = PROJECT_ROOT / "output" / "raw_runs"
EXP_ROOT = RAW / "experiments" / "source-sink"
LOG_ROOT = PROJECT_ROOT / "output" / "logs" / "driver" / "dem_to_dem"
# The pinned calibration run every solve starts from (end rates + blocking groin b 0.6 f 0.6)
START_RUNS = {1996: RAW / "experiments" / "groin" / "2026-10-05-blocking-fit-dem-to-dem" / "b0.60_f0.6"
              / "1996_2009" / "edgeBE" / "HAT_1996_2009_edgeBE_offsetmetres_road_bdm_groinblock"}
GIS = np.arange(HATTERAS_DOMAINS.first_gis_id, HATTERAS_DOMAINS.last_gis_id + 1)
# -----------------------------------------------------------------------------


def window(period):
    return f"{period}_{HATTERAS_PERIODS[period]['end_year']}"


def study(period):
    return EXP_ROOT / f"{SOLVE_DATE}-be-domain-solve-{window(period)}"


def tag(period, k):
    return f"source-sink/{study(period).name}/runs/step{k}"


def model_net(run_dir):
    m = np.load(next(Path(run_dir).glob("*_shoreline_matrix.npy")))
    D = HATTERAS_DOMAINS
    return pd.Series(-(m[-1] - m[0])[D.start_real_index:D.end_real_index], index=GIS)


def target(period):
    return pd.read_csv(net_change_domain_csv(period, HATTERAS_PERIODS[period]["end_year"]),
                       index_col=0)["net_change_lowess7_m"].reindex(GIS)


def step_run_dir(period, k):
    if k == 0:
        return START_RUNS[period]
    hits = sorted((study(period) / "runs" / f"step{k}").glob("*/*/*/*_shoreline_matrix.npy"))
    return hits[-1].parent if hits else None


def field_file(period, k):
    return study(period) / "fields" / f"be_field_step{k}.csv"


def run(period, k, field):
    override = ",".join(f"{g}={r:.4f}" for g, r in field.items())
    env = {key: v for key, v in os.environ.items() if not key.startswith("HAT_")}
    env.update(PYTHONIOENCODING="utf-8", PYTHONUNBUFFERED="1", MPLBACKEND="Agg",
               HAT_IGNORE_SETTINGS="1", HAT_START_YEAR=str(period),
               HAT_SOURCE_SINK_PRESET="edgeBE", HAT_SCENARIO=SCENARIO, HAT_RELOCATIONS="False",
               HAT_GROIN_ENABLED="True", HAT_OVERWRITE="False", HAT_MAKE_GIFS="False",
               HAT_SHOW_FIGURES="False", HAT_RUN_KIND="experiment", HAT_RUN_TAG=tag(period, k),
               HAT_SAVE_MODEL_STATE="False", HAT_BE_OVERRIDE=override)
    log = LOG_ROOT / f"be_domain_solve_{window(period)}_step{k}.log"
    log.parent.mkdir(parents=True, exist_ok=True)
    print(f"start step{k} -> {log.relative_to(PROJECT_ROOT)}", flush=True)
    with open(log, "w", encoding="utf-8") as f:
        code = subprocess.run([sys.executable, str(HINDCAST)], stdout=f, stderr=subprocess.STDOUT,
                              env=env, cwd=PROJECT_ROOT).returncode
    if code:
        raise SystemExit(f"step{k} failed ({code}); see {log}")


def score(resid):
    interior = resid.loc[2:89]
    return dict(bias_m=float(interior.mean()), rmse_m=float(np.sqrt((interior ** 2).mean())),
                max_abs_m=float(resid.abs().max()), worst_gis=int(resid.abs().idxmax()))


def cmd_solve(a):
    period = a.period
    years = run_years(period)
    want = target(period)
    (study(period) / "fields").mkdir(parents=True, exist_ok=True)
    field0 = pd.Series(0.0, index=GIS)
    for g, r in HATTERAS_BE_RATES_EDGE[period].items():
        field0[g] = r
    if not field_file(period, 0).is_file():
        field0.rename("be_m_yr").to_csv(field_file(period, 0))
    log_rows = []
    k = 0
    while True:
        rd = step_run_dir(period, k)
        field = pd.read_csv(field_file(period, k), index_col=0)["be_m_yr"].reindex(GIS)
        if rd is None:
            run(period, k, field)
            rd = step_run_dir(period, k)
        resid = want - model_net(rd)
        s = score(resid)
        log_rows.append(dict(step=k, run_dir=str(Path(rd).relative_to(PROJECT_ROOT)), **s,
                             be_min=float(field.min()), be_max=float(field.max())))
        pd.DataFrame(log_rows).round(4).to_csv(study(period) / "solve_log.csv", index=False)
        pd.DataFrame({"target_m": want, "model_m": model_net(rd), "residual_m": resid,
                      "be_m_yr": field}).round(4).to_csv(study(period) / "fields" / f"residual_step{k}.csv")
        print(f"step{k}: interior bias {s['bias_m']:+.2f} m, RMSE {s['rmse_m']:.2f} m, "
              f"max |miss| {s['max_abs_m']:.2f} m at GIS {s['worst_gis']}; "
              f"field {field.min():+.2f} to {field.max():+.2f} m/yr", flush=True)
        if s["rmse_m"] < RMSE_TOL_M and s["max_abs_m"] < MAX_TOL_M:
            print(f"converged at step{k}")
            break
        if k >= MAX_STEPS:
            print(f"stopped at MAX_STEPS ({MAX_STEPS})")
            break
        gain = pd.Series(float(years), index=GIS)
        for g, v in GAIN_OVERRIDE.items():
            gain[g] = v
        nxt = field + RELAX * resid / gain
        k += 1
        if not field_file(period, k).is_file():
            nxt.rename("be_m_yr").to_csv(field_file(period, k))
    print(f"field: {field_file(period, k).relative_to(PROJECT_ROOT)}")


def cmd_report(a):
    print(pd.read_csv(study(a.period) / "solve_log.csv").to_string(index=False))


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n", 2)[1])
    sub = ap.add_subparsers(dest="cmd", required=True)
    for name in ("solve", "report"):
        p = sub.add_parser(name)
        p.add_argument("--period", type=int, default=1996, choices=sorted(START_RUNS))
    a = ap.parse_args(argv)
    {"solve": cmd_solve, "report": cmd_report}[a.cmd](a)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
