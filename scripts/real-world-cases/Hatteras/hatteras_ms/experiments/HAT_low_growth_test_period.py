"""
Do no or very slow dune growth hold up in the 2009-2025 test window, scored on dunes, major events and shoreline?

    python scripts/hatteras_ms/experiments/HAT_low_growth_test_period.py run
    python scripts/hatteras_ms/experiments/HAT_low_growth_test_period.py score

Runs no growth and rmin-rmax 0.03-0.07 (both Dmaxel 5.0 m NAVD88) and the control
(5.5 m, 0.55-0.95), groin on, on 1996-2009 and 2009-2025, in the low growth-rate
study's folder through the grid study's launcher. Scores, per window: the dune
trajectory against that window's surveys (low-dune 10th percentile, share of
cells below 4.3 m, domain-median r); the major-event images, as
overwash_vs_model scores them (1996-2009: 2004-05-25, Isabel; 2009-2025:
2011-08-27, Irene; 2013-04-06, after Sandy; 2024-03-14); overall overwash skill;
and the shoreline RMSE. Major events first, Hannah 2026-10-19: one image a year
makes a false alarm weak evidence.

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
sys.path.insert(0, str(_HERE.parent))
import HAT_dune_growth_low_rates as L  # noqa: E402
G, U = L.G, L.U
import overwash_vs_model as ovm  # noqa: E402
from site_layer import hat_overwash as ow  # noqa: E402

# --- CONFIG ------------------------------------------------------------------
SETTINGS = {"no growth, 5.0 m": (5.0, (0.0, 0.0)), "0.03-0.07, 5.0 m": (5.0, (0.03, 0.07)),
            "control, 5.5 m, 0.55-0.95": G.CONTROL}
MAJOR = {1996: ("2004-05-25",), 2009: ("2011-08-27", "2013-04-06", "2024-03-14")}
OUT = L.EXP_DIR / "tables" / "test_period_check.csv"
# -----------------------------------------------------------------------------


def jobs():
    return [(s, G.name(c, r, groin=True)) for s in (1996, 2009) for c, r in SETTINGS.values()]


# Image x domain outcomes for one run, overwash_vs_model's rule
def cells(start, variant):
    c = U.load(G.run_dir(start, variant))
    q = np.array([np.asarray(b.QowTS) for b in c.barrier3d]).T
    orig = ovm.load_qow
    ovm.load_qow = lambda w, arm: (q, variant)
    try:
        obs = pd.read_csv(ow.OBSERVATIONS / "overwash_observations.csv")
        return ovm.compare((start, U.PERIODS[start][0]), "managed", obs, ovm.HEADLINE_THR).dropna(subset=["observed"])
    finally:
        ovm.load_qow = orig


def major(d, image):
    g = d[d.image.astype(str) == image]
    o, m = g.observed.astype(int), g.model.astype(int)
    return f"{int(((o == 1) & (m == 1)).sum())}/{int((o == 1).sum())} hit, {int(((o == 0) & (m == 1)).sum())} FA"


def main():
    if len(sys.argv) > 1 and sys.argv[1] == "run":
        G.run_all(jobs())
    rows = []
    for (start, variant), label in zip(jobs(), list(SETTINGS) * 2):
        if G.run_dir(start, variant) is None:
            print(f"missing {start} {variant}")
            continue
        m, tr, d = G.metrics(start, variant), L.trajectory(start, variant), cells(start, variant)
        row = dict(window=f"{start}_{U.PERIODS[start][0]}", setting=label, variant=variant,
                   low_dune_p10_rmse=round(tr["p10_rmse"], 3), low_dune_p10_bias=round(tr["p10_bias"], 3),
                   share_below_4p3_rmse=round(tr["below_rmse"], 3), domain_r=round(tr["domain_r_mean"], 3),
                   overall_skill=m["skill"], hits=m["hits"], false_alarms=m["false_alarms"],
                   shoreline_rmse=round(m["rmse_interior_m_yr"], 4))
        for image in MAJOR[start]:
            row[f"image_{image}"] = major(d, image)
        rows.append(row)
    t = pd.DataFrame(rows)
    OUT.parent.mkdir(parents=True, exist_ok=True)
    t.to_csv(OUT, index=False)
    with pd.option_context("display.width", 300, "display.max_columns", 30):
        for w, g in t.groupby("window"):
            print(g.dropna(axis=1, how="all").drop(columns=["variant"]).to_string(index=False), "\n")


if __name__ == "__main__":
    main()
