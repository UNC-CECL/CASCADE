"""
Cs refit with the 1996 transect toes: the collision rule comparison's offline fit (Larson et al. 2004 alone, per
storm, against the lidar crest changes of the calibration intervals), rerun with each domain's toe from
HAT_toe_method_line.py, as HAT_collision_toe_comparison.py does for its toe choices.

    python scripts/real-world-cases/Hatteras/input_prep/9-dune-toes/HAT_cs_refit_transect_toes.py

Toes (one per domain unless named):
    cell_median        per cell, median over the 1996-2009 surveys: the toe Cs = 1.9e-5 was first fit with
    current            the module's now: the mean of the 1996 and 1997 surveys' screened cells
    transect           the 1996 transects, front dune inside the model window, no lower than the berm (used)
    transect_window    the 1996 transects, the window's highest point as the crest
    transect_free      the 1996 transects, the front dune within 120 m of the shoreline
    automorph          the 1996 Automorph toes you checked by hand (HAT_toe_automorph_final.py)
Rules: no erosion, and Larson per storm alone (the module).

Writes, in the workspace's data/hatteras_init/9-dune-toes/1-transects-1996/7-cs-refit/ (D:): cs_refit_1996.csv.

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-11-05
"""
from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd

_HERE = Path(__file__).resolve()
_HAT = next(p for p in _HERE.parents if (p / "hatteras_paths.py").exists())
sys.path.insert(0, str(_HAT / "dune-collision" / "3-module-tests" / "1-rule-comparison"))
sys.path.insert(0, str(_HERE.parent))
import HAT_collision_rule_comparison as R  # noqa: E402
import HAT_toe_common as M  # noqa: E402

# --- CONFIG ------------------------------------------------------------------
TOES = {"cell_median": None, "current": "current_toe_m_navd88", "transect": "toe_m_navd88",
        "transect_window": "toe_window_m_navd88", "transect_free": "toe_free_m_navd88",
        "automorph": "automorph_toe_m_navd88"}
RULE_NAMES = ("null", "larson_storm")
OUT = M.OUT / "7-cs-refit" / f"cs_refit_{M.YEAR}.csv"
# -----------------------------------------------------------------------------


def main():
    dom = pd.read_csv(M.TOES / f"domain_toes_{M.YEAR}.csv").set_index("gis")
    dom["automorph_toe_m_navd88"] = pd.read_csv(M.TOES / f"automorph_domain_toes_{M.YEAR}.csv").set_index(
        "gis").toe_m_navd88
    orig = pd.read_csv
    rules = [r for r in R.RULES if r["name"] in RULE_NAMES]
    rows = []
    for name, col in TOES.items():
        def read_csv(path, *a, **k):
            t = orig(path, *a, **k)
            if col and Path(str(path)).name == "collision_calibration_cells.csv":
                t["toe_cell_m_navd88"] = t.gis.map(dom[col])
            return t
        pd.read_csv = read_csv
        try:
            c, ramp = R.features()
        finally:
            pd.read_csv = orig
        print(f"\n== {name}: toe median {c.toe_cell_m_navd88.median():.2f} m over {len(c)} calibration cells",
              flush=True)
        for rule in rules:
            row, _, _ = R.score(rule, c, ramp)
            rows.append(dict(toe=name, toe_median_m=c.toe_cell_m_navd88.median(), **row))
    t = pd.DataFrame(rows)
    t.round(6).to_csv(OUT, index=False)
    print()
    print(t[t.rule == "larson_storm"][["toe", "toe_median_m", "cs", "growth_m_yr", "cv_rmse", "rmse", "r",
                                       "dennis_r", "dennis_rmse"]].round(4).to_string(index=False))
    print(f"no erosion: cv {t[t.rule == 'null'].cv_rmse.iloc[0]:.3f}")
    print(f"wrote {OUT}")


if __name__ == "__main__":
    main()
