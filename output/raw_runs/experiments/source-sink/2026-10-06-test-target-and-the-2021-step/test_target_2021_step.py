#!/usr/bin/env python3
"""
How much of the 2009-2025 test misfit is the 2020-2021 CoastSat step?

    python test_target_2021_step.py   ->  scores.csv, domain_residuals.csv (beside this script)

Scores the existing test runs three ways, interior GIS 2-89, both sides 7-domain LOWESS:
  full      the 2025 target (2025-08-17 +/-6 mo minus the 2009 start mean) vs the 1 Jan 2025 state
  pre_step  2019-08-17 +/-1 yr minus the 2009 start mean, vs the 1 Jan 2019 state (10 model years)
  no_step   the full target minus each domain's measured 2020->2021 step (pre-fill table)

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-10-06
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
REPO = next(p for p in HERE.parents if (p / "pyproject.toml").exists())
sys.path[:0] = [str(REPO / "scripts"),
                str(REPO / "scripts" / "input_prep" / "5-scr" / "3-rates" / "coastsat" / "net_change")]
from coastsat_net_change import domain_net_change, smooth_like_model, window_means  # noqa: E402
from site_layer.hat_observed_rates import (  # noqa: E402
    NET_CHANGE_WINDOWS, OBSERVATIONS, net_change_domain_csv)
from site_layer.hatteras_site_config import HATTERAS_DOMAINS as D  # noqa: E402

# --- CONFIG ------------------------------------------------------------------
RAW = REPO / "output" / "raw_runs"
START_WINDOW = NET_CHANGE_WINDOWS[(2009, 2025)][0]
PRE_STEP_WINDOW = ("2018-08-17", "2020-08-17")
PRE_STEP_INDEX = 10                     # 1 Jan 2019, the same 7.5-month offset as every target
STEP_TABLE = OBSERVATIONS / "detrended_position" / "step_2021_by_domain_prefill.csv"
GIS = pd.RangeIndex(1, 91)
RUNS = {
    "zeroBE, no groin": RAW / "matrix/2009_2025/zeroBE/HAT_2009_2025_zeroBE_offsetmetres_road_bdm_nourish_nogroin",
    "edgeBE, no groin": RAW / "matrix/2009_2025/edgeBE/HAT_2009_2025_edgeBE_offsetmetres_road_bdm_nourish_nogroin",
    "edgeBE + groin": RAW / ("experiments/groin/2026-10-05-blocking-fit-dem-to-dem/b0.60_f0.6/2009_2025/edgeBE/"
                             "HAT_2009_2025_edgeBE_offsetmetres_road_bdm_nourish_groinblock"),
    "domainBE + groin": RAW / "matrix/2009_2025/domainBE/HAT_2009_2025_domainBE_offsetmetres_road_bdm_nourish_groinblock",
}
# -----------------------------------------------------------------------------


# Observed change between two window means, per domain: raw and smoothed
def observed(start_w, end_w):
    s, e = window_means(start_w), window_means(end_w)
    both = s.index.intersection(e.index)
    tr = pd.DataFrame({"domain_number": s.loc[both, "domain_number"].astype(int),
                       "net_change_m": e.loc[both, "mean_chainage_m"] - s.loc[both, "mean_chainage_m"],
                       "se_net_change_m": np.hypot(s.loc[both, "se_chainage_m"],
                                                   e.loc[both, "se_chainage_m"])})
    return domain_net_change(tr)["net_change_lowess7_m"]


def model_change(rd, index):
    m = np.load(next(Path(rd).glob("*_shoreline_matrix.npy")))
    return pd.Series(-(m[index] - m[0])[D.start_real_index:D.end_real_index], index=GIS)


def score(model, target):
    m, o = smooth_like_model(model).loc[2:89], target.loc[2:89]
    d = m - o
    return dict(bias_m=d.mean(), rmse_m=float(np.sqrt((d ** 2).mean())),
                r=float(np.corrcoef(m, o)[0, 1]), obs_mean_m=o.mean(), model_mean_m=m.mean())


def main():
    full = pd.read_csv(net_change_domain_csv(2009, 2025), index_col=0)["net_change_lowess7_m"]
    pre = observed(START_WINDOW, PRE_STEP_WINDOW)
    step = pd.read_csv(STEP_TABLE).set_index("domain_number")["step_2020_2021_m"].reindex(GIS)
    no_step = full - smooth_like_model(step)
    targets = {"full": (full, -1), "pre_step": (pre, PRE_STEP_INDEX), "no_step": (no_step, -1)}
    rows, resid = [], {}
    for run, rd in RUNS.items():
        for name, (tgt, idx) in targets.items():
            mod = model_change(rd, idx)
            rows.append(dict(run=run, target=name, **score(mod, tgt)))
            if run == "domainBE + groin":
                resid[f"{name}_residual_m"] = tgt - smooth_like_model(mod)
    t = pd.DataFrame(rows).round(2)
    t.to_csv(HERE / "scores.csv", index=False)
    pd.DataFrame(dict(resid, step_2020_2021_m=step, pre_step_target_m=pre, full_target_m=full,
                      no_step_target_m=no_step)).round(3).to_csv(HERE / "domain_residuals.csv")
    pd.set_option("display.width", 160)
    print(t.to_string(index=False))
    print(f"\nisland step (interior mean, smoothed): {smooth_like_model(step).loc[2:89].mean():+.1f} m")


if __name__ == "__main__":
    main()
