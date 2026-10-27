"""
Which dune growth rate, down to none, best follows the surveyed dunes through 1996-2009?

    python scripts/real-world-cases/Hatteras/hatteras_ms/experiments/dunes/HAT_dune_growth_low_rates.py run      # all stages, in order
    python scripts/real-world-cases/Hatteras/hatteras_ms/experiments/dunes/HAT_dune_growth_low_rates.py score    # rank from runs on disk

Stage 1: rmin-rmax 0 (no growth), 0.01-0.03, 0.03-0.07, 0.05-0.15, 0.10-0.20 x
Dmaxel 5.0 / 5.5 m NAVD88, plus the control (5.5 m, 0.55-0.95), 1996-2009 only.
Each run is scored against the nine calibration surveys (1996-start/survey-crests:
the 1997-2005 lidar grids and the 2009 lidar, on the 1996 dune windows), each at
the model state its flight dates point to.

Ranking agreed with Hannah (2026-10-19), dune trajectory first: rank by the sum of
two ranks, the RMSE over the surveys of the low-dune 10th percentile (model minus
survey, m) and of the share of cells below 4.3 m (model minus survey); ties
broken by the mean alongshore domain-median r. Guards: overwash skill at most
GUARD_SKILL below the control, shoreline RMSE at most GUARD_RMSE above it.
Stage 2: the top three passing the guards and the control re-run with three
other growth-rate seeds; the ranking is redone on the seed means. Stage 3: the
winner and the control on 2009-2025, groin on, scored the same way against the
2014-2019 surveys, plus the winner on 1996-2009 with the groin on.
Managed edgeBE, in-process only, through HAT_dune_ceiling_growth_grid.py's launcher.
Details: scripts/real-world-cases/Hatteras/hatteras_ms/experiments/dunes/README.md.

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-10-19
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

_HERE = Path(__file__).resolve()
sys.path.insert(0, str(next(_p for _p in _HERE.parents if (_p / "experiment_paths.py").exists())))
import experiment_paths  # noqa: E402,F401  (every experiments/<theme>/ folder onto sys.path)
import HAT_dune_ceiling_uniform as U  # noqa: E402
import HAT_dune_ceiling_growth_grid as G  # noqa: E402
import HAT_dune_crest_vs_surveys as V  # noqa: E402

# --- CONFIG ------------------------------------------------------------------
TAG = "dunes/2026-10-19-dune-growth-low-rates"
CEILINGS = (5.0, 5.5)
RATES = ((0.0, 0.0), (0.01, 0.03), (0.03, 0.07), (0.05, 0.15), (0.10, 0.20))
TOP_N = 3
GUARD_SKILL = 0.016          # two seed SDs of overwash skill (grid study: 0.008)
GUARD_RMSE = 0.005           # m/yr
LOW_PCT = 10
DESIGN_NAVD88 = 4.3
# -----------------------------------------------------------------------------

# The grid study's launcher and run folders, re-pointed at this study
G.TAG, G.EXP_DIR = TAG, U.PROJECT_ROOT / "output" / "raw_runs" / "experiments" / TAG
EXP_DIR = G.EXP_DIR


# Model against every survey of the window: low-dune P10, share below 4.3 m, domain-median r
def trajectory(start, variant):
    prod, first = V.PERIODS[start]
    if G.run_dir(start, variant) is None:
        return None
    mod = V.model_crests(start, variant)
    rows = []
    for year, g in V.surveys(start).groupby("survey"):
        k = int(g.model_state.iloc[0][:4]) - first
        if k <= 0 or k >= mod.shape[0]:
            continue
        g = g.sort_values(["gis", "cell"])
        obs, m = g.crest.to_numpy(), mod[k]
        ok = np.isfinite(obs)
        c = V.compare(obs, m)
        rows.append(dict(survey=year, p10_err=np.percentile(m[ok], LOW_PCT) - np.percentile(obs[ok], LOW_PCT),
                         below_err=float((m[ok] < DESIGN_NAVD88).mean() - (obs[ok] < DESIGN_NAVD88).mean()),
                         domain_r=c["domain_median_r"]))
    t = pd.DataFrame(rows)
    return dict(p10_rmse=float(np.sqrt((t.p10_err ** 2).mean())), p10_bias=float(t.p10_err.mean()),
                below_rmse=float(np.sqrt((t.below_err ** 2).mean())), below_bias=float(t.below_err.mean()),
                domain_r_mean=float(t.domain_r.mean()), n_surveys=len(t))


def score(jobs):
    rows = []
    for start, variant in jobs:
        m = G.metrics(start, variant)
        if m is None:
            continue
        m.update(trajectory(start, variant))
        rows.append(m)
    return pd.DataFrame(rows)


def rank(df):
    d = df.copy()
    d["combo"] = d.variant.str.replace(r"_s\d+(_groin)?$", "", regex=True)
    agg = d.groupby("combo").agg(ceiling_m_navd88=("ceiling_m_navd88", "first"), rmin=("rmin", "first"),
                                 rmax=("rmax", "first"), n_seeds=("seed", "size"),
                                 p10_rmse=("p10_rmse", "mean"), p10_bias=("p10_bias", "mean"),
                                 below_rmse=("below_rmse", "mean"), below_bias=("below_bias", "mean"),
                                 domain_r_mean=("domain_r_mean", "mean"), skill=("skill", "mean"),
                                 skill_sd=("skill", "std"), hits=("hits", "mean"),
                                 false_alarms=("false_alarms", "mean"), rmse=("rmse_interior_m_yr", "mean"))
    ctrl = agg.loc["control"]
    agg["passes_guards"] = (agg.skill >= ctrl.skill - GUARD_SKILL) & (agg.rmse <= ctrl.rmse + GUARD_RMSE)
    agg["rank_sum"] = agg.p10_rmse.rank() + agg.below_rmse.rank()
    return agg.sort_values(["passes_guards", "rank_sum", "domain_r_mean"], ascending=[False, True, False]).round(4)


def write(df, name, index=False):
    (EXP_DIR / "tables").mkdir(parents=True, exist_ok=True)
    df.to_csv(EXP_DIR / "tables" / name, index=index)


def stages(run):
    with pd.option_context("display.width", 260, "display.max_columns", 30):
        stage1 = [(1996, G.name(c, r)) for c in CEILINGS for r in RATES] + [(1996, G.name(*G.CONTROL))]
        if run:
            G.run_all(stage1)
        s1 = score(stage1)
        write(s1, "stage1_1996.csv")
        r1 = rank(s1)
        write(r1, "stage1_ranking_1996.csv", index=True)
        print(r1.to_string(), flush=True)

        top = [c for c in r1[r1.passes_guards].index if c != "control"][:TOP_N]
        if not top:
            print("no candidate passes the guards: the dune trajectory and overwash skill disagree; stopping "
                  "before the seed and test stages", flush=True)
            return
        combos = [G.parse(c + "_s0")[:2] for c in top] + [G.CONTROL]
        stage2 = [(1996, G.name(c, r, s)) for c, r in combos for s in G.SEEDS]
        if run:
            G.run_all(stage2)
        s2 = score(stage2)
        write(s2, "stage2_seeds_1996.csv")
        r2 = rank(pd.concat([s1[s1.variant.str.replace(r"_s0$", "", regex=True).isin(top + ["control"])], s2]))
        write(r2, "stage2_ranking_1996.csv", index=True)
        print(r2.to_string(), flush=True)
        winner = next(c for c in r2[r2.passes_guards].index if c != "control")
        wc, wr = G.parse(winner + "_s0")[:2]
        print(f"winner {winner}", flush=True)

        stage3 = [(2009, G.name(wc, wr, groin=True)), (2009, G.name(*G.CONTROL, groin=True)),
                  (1996, G.name(wc, wr, groin=True)), (1996, G.name(*G.CONTROL, groin=True))]
        if run:
            G.run_all(stage3)
        s3 = score(stage3)
        write(s3, "stage3_winner_vs_control.csv")
        print(s3.to_string(index=False), flush=True)


def main():
    stages(run=len(sys.argv) > 1 and sys.argv[1] == "run")


if __name__ == "__main__":
    main()
