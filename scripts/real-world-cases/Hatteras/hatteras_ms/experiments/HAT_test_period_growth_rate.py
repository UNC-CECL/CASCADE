"""
Did the dunes build faster after 2009? Which growth rate reproduces the 2009-2025 dunes under Dmax 10.5 m?

    python scripts/hatteras_ms/experiments/HAT_test_period_growth_rate.py run
    python scripts/hatteras_ms/experiments/HAT_test_period_growth_rate.py score

Asked by Hannah 2026-10-20: the island-wide median crest rose ~0.06 m/yr over
2009-2019 against ~0 over 1996-2009, so test whether the test window wants a
faster rate than the calibrated 0.01-0.03. Base Barrier3D growth, only rmin-rmax
changes. The 2009 start, Dmaxel 12.2 m, groin on, today's setup and the ends
re-solved 2026-10-20 (for 0.01-0.03), through the grid study's launcher. The
0.01-0.03 reference is the same setup's run in
dunes/2026-10-20-shoreline-rescore-new-ends. Scored against the 2014-2019 lidar
on the 2009 dune windows (island-wide median and low-dune 10th-percentile
crest per survey, as the 2009-start survey-crest table holds them), and on
overwash (Irene, after Sandy, March 2024) and the interior shoreline.

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-10-20
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

_HERE = Path(__file__).resolve()
sys.path.insert(0, str(_HERE.parent))
import HAT_dune_ceiling_uniform as U  # noqa: E402
import HAT_dune_ceiling_growth_grid as G  # noqa: E402

# --- CONFIG ------------------------------------------------------------------
TAG = "dunes/2026-10-20-test-period-growth-rate"
EXP_DIR = U.PROJECT_ROOT / "output" / "raw_runs" / "experiments" / TAG
# The table: the study set the runner's 2009 growth rate, so since 2026-10-21 it sits with the 2009 start's
# inputs; the runs stay in EXP_DIR
from site_layer import hat_topo_version as _tv  # noqa: E402
OUT_DIR = _tv.product_dir("2009-start") / "dune-growth-rate"
REFERENCE = ("dunes/2026-10-20-shoreline-rescore-new-ends",
             U.PROJECT_ROOT / "output" / "raw_runs" / "experiments" / "dunes" / "2026-10-20-shoreline-rescore-new-ends")
CEILING = 12.2
ADOPTED = (0.01, 0.03)
RATES = ((0.02, 0.04), (0.03, 0.07), (0.05, 0.15), (0.10, 0.20), (0.25, 0.45))
START = 2009
# -----------------------------------------------------------------------------


def point(tag, exp_dir):
    G.TAG, G.EXP_DIR = tag, exp_dir


def jobs():
    point(TAG, EXP_DIR)
    return [(START, G.name(CEILING, r, groin=True)) for r in RATES]


# Model against each 2014-2019 survey on the cells it covers: island-wide median and 10th percentile
def survey_errors(variant):
    import HAT_dune_crest_vs_surveys as V
    mod = V.model_crests(START, variant)
    first = V.PERIODS[START][1]
    rows = []
    for year, g in V.surveys(START).groupby("survey"):
        k = int(g.model_state.iloc[0][:4]) - first
        g = g.sort_values(["gis", "cell"])
        obs, m = g.crest.to_numpy(), mod[k]
        ok = np.isfinite(obs)
        rows.append(dict(survey=int(year), median_err=np.median(m[ok]) - np.median(obs[ok]),
                         p10_err=np.percentile(m[ok], 10) - np.percentile(obs[ok], 10)))
    return pd.DataFrame(rows)


def score(tag, exp_dir, variant, label):
    point(tag, exp_dir)
    import HAT_low_growth_test_period as T          # importing re-points the launcher; point back
    point(tag, exp_dir)
    e = survey_errors(variant)
    d = T.cells(START, variant)
    o, mo = d.observed.astype(int), d.model.astype(int)
    meta = json.loads(next(G.run_dir(START, variant).glob("*_run_metadata.json")).read_text(encoding="utf-8"))
    row = dict(rate=label,
               median_err_2014_2019=" / ".join(f"{v:+.2f}" for v in e.median_err),
               median_rmse=round(float(np.sqrt((e.median_err ** 2).mean())), 2),
               p10_err_2014_2019=" / ".join(f"{v:+.2f}" for v in e.p10_err),
               p10_rmse=round(float(np.sqrt((e.p10_err ** 2).mean())), 2),
               hits=int(((o == 1) & (mo == 1)).sum()), observed=int(o.sum()),
               model_only=int(((o == 0) & (mo == 1)).sum()),
               shoreline_rmse=round(float(meta["skill"]["rmse_interior_m_yr"]), 3),
               shoreline_bias=round(float(meta["skill"]["mean_bias_interior_m_yr"]), 3))
    for image in T.MAJOR[START]:
        row[f"image_{image}"] = T.major(d, image)
    return row


def main():
    if len(sys.argv) > 1 and sys.argv[1] == "run":
        G.run_all(jobs())
    rows = [score(*REFERENCE, G.name(CEILING, ADOPTED, groin=True), f"{ADOPTED[0]:.2f}-{ADOPTED[1]:.2f} (adopted)")]
    for (start, variant), r in zip(jobs(), RATES):
        point(TAG, EXP_DIR)
        if G.run_dir(start, variant) is None:
            print(f"missing {variant}")
            continue
        rows.append(score(TAG, EXP_DIR, variant, f"{r[0]:.2f}-{r[1]:.2f}"))
    t = pd.DataFrame(rows)
    (OUT_DIR / "tables").mkdir(parents=True, exist_ok=True)
    t.to_csv(OUT_DIR / "tables" / "test_period_growth_scores.csv", index=False)
    with pd.option_context("display.width", 300, "display.max_columns", 30):
        print(t.to_string(index=False))


if __name__ == "__main__":
    main()
