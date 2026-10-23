"""
Is the shoreline worse because of the new dune settings, or because of the rest of the 2026-10-20 setup?

    python scripts/hatteras_ms/experiments/HAT_shoreline_rescore_new_ends.py

The adopted dunes (Dmaxel 12.2 m, rmin-rmax 0.01-0.03) and the previous ones (5.5 m,
0.55-0.95), both on today's setup (ALACE-corrected 1996 start, dune-topo v3; the
v3_240 storms and beach slope 0.06 from another session) with the end rates
re-solved 2026-10-20 for the adopted dunes, groin on, on 1996-2009 and 2009-2025,
through the grid study's launcher. Scores the interior shoreline (LRR RMSE and mean
bias, GIS 2-89) and overwash skill. If the previous dunes fit the shoreline on this
setup and the adopted ones do not, the dunes are the cause; if neither does, the
rest of the setup is. The previous dunes keep the adopted dunes' ends, which
favours the adopted ones slightly near GIS 1 and 90 only.

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-10-20
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import pandas as pd

_HERE = Path(__file__).resolve()
sys.path.insert(0, str(_HERE.parent))
import HAT_dune_ceiling_uniform as U  # noqa: E402
import HAT_dune_ceiling_growth_grid as G  # noqa: E402

# --- CONFIG ------------------------------------------------------------------
TAG = "dunes/2026-10-20-shoreline-rescore-new-ends"
SETTINGS = {"adopted: 12.2 m, 0.01-0.03": (12.2, (0.01, 0.03)),
            "previous: 5.5 m, 0.55-0.95": (5.5, (0.55, 0.95))}
# -----------------------------------------------------------------------------

G.TAG, G.EXP_DIR = TAG, U.PROJECT_ROOT / "output" / "raw_runs" / "experiments" / TAG
EXP_DIR = G.EXP_DIR


def jobs():
    return [(s, G.name(c, r, groin=True)) for s in (1996, 2009) for c, r in SETTINGS.values()]


# The runner's interior mean bias, from the run's own metadata
def bias(start, variant):
    meta = next(G.run_dir(start, variant).glob("*_run_metadata.json"))
    return float(json.loads(meta.read_text(encoding="utf-8"))["skill"]["mean_bias_interior_m_yr"])


def main():
    G.run_all(jobs())
    rows = []
    for (start, variant), label in zip(jobs(), list(SETTINGS) * 2):
        m = G.metrics(start, variant)
        rows.append(dict(window=f"{start}_{U.PERIODS[start][0]}", dunes=label, variant=variant,
                         shoreline_rmse=round(m["rmse_interior_m_yr"], 3),
                         shoreline_bias=round(bias(start, variant), 3),
                         overwash_skill=m["skill"], hits=m["hits"], false_alarms=m["false_alarms"]))
    t = pd.DataFrame(rows)
    (EXP_DIR / "tables").mkdir(parents=True, exist_ok=True)
    t.to_csv(EXP_DIR / "tables" / "shoreline_rescore.csv", index=False)
    print(t.drop(columns="variant").to_string(index=False))


if __name__ == "__main__":
    main()
