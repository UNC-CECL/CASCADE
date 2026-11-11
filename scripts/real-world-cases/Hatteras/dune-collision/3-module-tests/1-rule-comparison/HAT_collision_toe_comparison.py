"""
Which dune toe should the collision module read? The rule comparison's fit, rerun on six toe choices.

    python scripts/real-world-cases/Hatteras/dune-collision/3-module-tests/1-rule-comparison/HAT_collision_toe_comparison.py

The module's erosion grows with the square of the water above the toe (Larson et al. 2004), so a toe too low
costs more than one too high saves, and a noisy single profile over-erodes (the 1996 start toes, after Bertha and
Fran, put the calibration crest RMSE from 0.71 to 1.33 m in the model, 2026-11-04). Here, offline on the same
calibration cells and survey intervals as HAT_collision_rule_comparison.py, each rule is refit with the toe:
    cell_median     per cell, the median over the 1996-2009 surveys (1998 left out): the module's now
    domain_median   per domain, the median of those
    island          one toe for the island, the median of those
    start_cell      per cell, the 1996 start profile's (screened as the module's; the domain median where none)
    start_domain    per domain, the median of the start profile's
    survey          per cell, each interval's own starting survey (screened; the cell median where none)
The beach widths and everything else are as the rule comparison's. Rules: no erosion, Larson per storm alone,
Larson per storm x the dry-beach factor (W0 25 m, the module; and W0 fit) and Long's 10% per storm x the
dry-beach factor.

Writes, in the workspace's dune-collision/3-module-tests/1-rule-comparison/tables/: toe_comparison.csv.

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-11-04
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

_HERE = Path(__file__).resolve()
sys.path.insert(0, str(_HERE.parent))
import HAT_collision_rule_comparison as R  # noqa: E402

H = R.H
TOES = ("cell_median", "domain_median", "island", "start_cell", "start_domain", "survey")
RULE_NAMES = ("null", "larson_storm", R.MODULE_RULE, "larson_storm_dry", "long10_storm_dry", "long10_storm_beach")


# Each calibration cell's toe (m NAVD88) under one choice
def toe_column(c, choice):
    med = c.toe_cell_m_navd88
    if choice == "cell_median":
        return med
    if choice == "domain_median":
        return med.groupby(c.gis).transform("median")
    if choice == "island":
        return pd.Series(med.median(), index=c.index)
    ok = lambda t: (t.toe_m_navd88 > H.TOE_MIN_M_NAVD88) & (t.toe_m_navd88 < t.crest_m_navd88 - H.TOE_BELOW_CREST_M)
    if choice == "survey":
        own = c.toe_m_navd88.where(ok(c))
        return own.fillna(med)
    sp = pd.read_csv(H.EXP_DIR / "dune_toes.csv")
    sp = sp[sp.survey.astype(str) == "start"]
    sp = sp[ok(sp)].set_index(["gis", "cell"]).toe_m_navd88
    if choice == "start_domain":
        v = c.gis.map(sp.groupby(level="gis").median())
    else:
        v = pd.Series(sp.reindex(pd.MultiIndex.from_arrays([c.gis, c.cell])).to_numpy(), index=c.index)
        v = v.fillna(c.gis.map(sp.groupby(level="gis").median()))
    return v.fillna(sp.median())


def main():
    orig = pd.read_csv
    rules = [r for r in R.RULES if r["name"] in RULE_NAMES]
    rows = []
    for choice in TOES:
        def read_csv(path, *a, **k):
            t = orig(path, *a, **k)
            if Path(str(path)).name == "collision_calibration_cells.csv":
                t["toe_cell_m_navd88"] = toe_column(t, choice)
            return t
        pd.read_csv = read_csv
        try:
            c, ramp = R.features()
        finally:
            pd.read_csv = orig
        spread = c.groupby("gis").toe_cell_m_navd88.std().median()
        print(f"\n== {choice}: toe median {c.toe_cell_m_navd88.median():.2f} m, within-domain spread {spread:.2f} m",
              flush=True)
        for rule in rules:
            row, _, _ = R.score(rule, c, ramp)
            rows.append(dict(toe=choice, toe_median_m=c.toe_cell_m_navd88.median(), within_domain_sd_m=spread, **row))
    t = pd.DataFrame(rows)
    out = R.OUT / "tables" / "toe_comparison.csv"
    out.parent.mkdir(parents=True, exist_ok=True)
    t.round(6).to_csv(out, index=False)
    pd.set_option("display.width", 200)
    for col in ("cv_rmse", "r", "dennis_r", "dennis_rmse"):
        print(f"\n{col}")
        print(t.pivot_table(index="rule", columns="toe", values=col, sort=False)[list(TOES)].round(3).to_string())
    print(f"\nwrote {out}")


if __name__ == "__main__":
    main()
