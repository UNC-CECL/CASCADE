"""
Which dune growth rate reproduces the stable dunes of 1996-2009 under Dmax 10.5 m above the berm?

    python scripts/real-world-cases/Hatteras/hatteras_ms/experiments/dunes/HAT_stable_dune_growth_rate.py run
    python scripts/real-world-cases/Hatteras/hatteras_ms/experiments/dunes/HAT_stable_dune_growth_rate.py score
    python scripts/real-world-cases/Hatteras/hatteras_ms/experiments/dunes/HAT_stable_dune_growth_rate.py check   # island-wide + 2009-2025

Hannah's advisor's method (2026-10-20): set r by matching a stable dune between 1996
and 2009 with the new ceiling (Dmaxel 12.2 m NAVD88 = berm 1.7 + 10.5). Stable
domains, picked from the nine 1996-2009 lidar surveys on the model's 10 m dune
cells and CoastSat: no overwash in any 1996-2009 image, shoreline LRR within
+/-0.5 m/yr, the dune present in every survey and changing gently without dips:
GIS 59-62 (Mid-island / Wimble Shoals) and 20-22 (Buxton-Avon / Avon). Their
crests (6.0-7.3 m) sit above Isabel's run-up (5.58 m), so the model's storms
should rarely touch them and growth alone sets their change.
Runs: the 1996 start, managed edgeBE, groin on, through the grid study's launcher,
one per rate (RATES). Score per domain group: model minus lidar domain-median
crest at each survey's model state, on the cells the survey covers; the headline
is 1 January 2010 against the 2009 lidar. Figures: HAT_stable_dune_figures.py.
"check" runs CHECK_RATES on both starts and scores them island-wide as
HAT_low_growth_test_period.py does (dune trajectory, major-event images, overwash
skill, shoreline RMSE), to compare 0.005-0.01 with the adopted 0-0.01.

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-10-20
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
from site_layer import hat_figure_style as fs  # noqa: E402

# --- CONFIG ------------------------------------------------------------------
# One study folder per 1996-start topography: v2 (before the 2026-10-20 ALACE offset correction) keeps the
# original folder, later versions get their own so runs on different starting DEMs never mix
from site_layer.hat_topo_version import topo_dirs as _topo_dirs  # noqa: E402
TOPO = _topo_dirs("1996-start")[2]
TAG = "dunes/2026-10-20-stable-dune-growth-rate" + ("" if TOPO == "v2" else f"-topo-{TOPO}")
BERM_NAVD88 = 1.7            # m; Barrier3D BermEl, the same in every domain
CEILING = round(BERM_NAVD88 + 10.5, 2)     # Dmax 10.5 m above the berm
RATES = ((0.0, 0.0), (0.0, 0.01), (0.005, 0.01), (0.01, 0.01), (0.01, 0.02), (0.02, 0.02), (0.02, 0.03), (0.03, 0.03),
         (0.04, 0.04), (0.01, 0.03), (0.005, 0.035))   # last two: wide ranges about 0.02, for the cell spread
GROUPS = {"GIS 59-62": (59, 60, 61, 62), "GIS 20-22": (20, 21, 22)}
GROUP_COLOURS = {"GIS 59-62": fs.C["ACCENT"], "GIS 20-22": "#e08214"}
CELLS = 50                   # dune cells per domain
START = 1996
CHECK_RATES = ((0.01, 0.02), (0.02, 0.02), (0.01, 0.03))
# -----------------------------------------------------------------------------

G.TAG, G.EXP_DIR = TAG, U.PROJECT_ROOT / "output" / "raw_runs" / "experiments" / TAG
EXP_DIR = G.EXP_DIR
# Tables and figures: the study set the runner's 1996 growth rate, so since 2026-10-21 they sit with the 1996
# start's inputs; the runs stay in EXP_DIR. v2 (before the ALACE correction) keeps its experiment folder
from site_layer import hat_topo_version as _tv  # noqa: E402
OUT_DIR = EXP_DIR if TOPO == "v2" else _tv.product_dir("1996-start") / "dune-growth-rate"


def jobs():
    return [(START, G.name(CEILING, r, groin=True)) for r in RATES]


# Model minus lidar domain-median crest, per survey and domain group
def score(variant):
    mod = V.model_crests(START, variant)
    first = V.PERIODS[START][1]
    rows = []
    for year, g in V.surveys(START).groupby("survey"):
        k = int(g.model_state.iloc[0][:4]) - first
        if k < 0 or k >= mod.shape[0]:
            continue
        for label, doms in GROUPS.items():
            s = g[g.gis.isin(doms)]
            ok = np.isfinite(s.crest.to_numpy())
            idx = ((s.gis - 1) * CELLS + s.cell).to_numpy()[ok]
            obs, m = s.crest.to_numpy()[ok], mod[k][idx]
            rows.append(dict(survey=year, model_state=g.model_state.iloc[0], group=label,
                             lidar_median=np.median(obs), model_median=np.median(m),
                             error=np.median(m) - np.median(obs), cells=int(ok.sum())))
    return pd.DataFrame(rows)


def check():
    import HAT_low_growth_test_period as T            # importing it re-points the launcher; point it back
    G.TAG, G.EXP_DIR = TAG, EXP_DIR
    L = T.L
    todo = [(s, G.name(CEILING, r, groin=True), r) for s in (1996, 2009) for r in CHECK_RATES]
    G.run_all([(s, v) for s, v, _ in todo])
    rows = []
    for start, variant, r in todo:
        m, tr, d = G.metrics(start, variant), L.trajectory(start, variant), T.cells(start, variant)
        row = dict(window=f"{start}_{U.PERIODS[start][0]}", rate=f"{r[0]:.3f}-{r[1]:.3f}",
                   low_dune_p10_rmse=round(tr["p10_rmse"], 3), low_dune_p10_bias=round(tr["p10_bias"], 3),
                   share_below_4p3_rmse=round(tr["below_rmse"], 3), domain_r=round(tr["domain_r_mean"], 3),
                   overall_skill=m["skill"], hits=m["hits"], false_alarms=m["false_alarms"],
                   shoreline_rmse=round(m["rmse_interior_m_yr"], 4))
        for image in T.MAJOR[start]:
            row[f"image_{image}"] = T.major(d, image)
        rows.append(row)
    t = pd.DataFrame(rows)
    (OUT_DIR / "tables").mkdir(parents=True, exist_ok=True)
    t.to_csv(OUT_DIR / "tables" / "island_and_test_check.csv", index=False)
    with pd.option_context("display.width", 300, "display.max_columns", 30):
        for w, g in t.groupby("window"):
            print(g.dropna(axis=1, how="all").to_string(index=False), "\n")


def main():
    if len(sys.argv) > 1 and sys.argv[1] == "check":
        return check()
    if len(sys.argv) > 1 and sys.argv[1] == "run":
        G.run_all(jobs())
    out = []
    for (start, variant), r in zip(jobs(), RATES):
        if G.run_dir(start, variant) is None:
            print(f"missing {variant}")
            continue
        t = score(variant)
        t.insert(0, "rate", f"{r[0]:.3f}-{r[1]:.3f}")
        out.append(t)
    t = pd.concat(out)
    (OUT_DIR / "tables").mkdir(parents=True, exist_ok=True)
    t.round(3).to_csv(OUT_DIR / "tables" / "stable_dune_scores.csv", index=False)
    end = t[t.survey == 2009].pivot(index="rate", columns="group", values="error")
    traj = t[t.survey > START].groupby(["rate", "group"]).error.apply(lambda e: np.sqrt((e ** 2).mean())).unstack()
    summary = pd.concat({"2010 error vs 2009 lidar (m)": end, "RMSE over the surveys (m)": traj}, axis=1)
    summary.round(3).to_csv(OUT_DIR / "tables" / "stable_dune_summary.csv")
    print(t[t.survey == 2009].groupby("group")[["lidar_median"]].first().round(2).to_string(), "\n")
    print(summary.round(2).to_string())


if __name__ == "__main__":
    main()
