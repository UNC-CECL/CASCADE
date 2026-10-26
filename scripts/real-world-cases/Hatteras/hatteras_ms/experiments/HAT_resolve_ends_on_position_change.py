"""
Solve the two end domains on position change, not LRR.

    python scripts/hatteras_ms/experiments/HAT_resolve_ends_on_position_change.py

The metres solver with its residual swapped for (model - observed) change / 14 yr,
at option A, full management; the solver and the config are unchanged. Details: scripts/hatteras_ms/experiments/README.md.

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-09-30
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import pandas as pd

_HERE = Path(__file__).resolve()
sys.path.insert(0, str(_HERE.parent))
import HAT_resolve_ends_metres as R  # noqa: E402
import HAT_wave_recommendation_figures as W  # noqa: E402

E, common = R.E, R.common
# --- CONFIG ------------------------------------------------------------------
TAG = "end-domain-boundaries/2026-09-28-ends-solved-on-position-change"
YEARS = 14
SEED = "1996=4.8394,18.2545;2010=18.8657,24.2358"
OPTION_A = ["--hs", "2.0", "--tp", "7.5", "--asym", "0.6", "--ahf", "0.5"]
OBS = {p: W.observed_change_smoothed(p) for p in (1996, 2010)}
# -----------------------------------------------------------------------------


# (model - observed) end-minus-start change, per year
def residuals(run_dir, period, targets):
    t = pd.read_csv(Path(run_dir) / "tables" / "shoreline_change_rate.csv").set_index("gis_domain")
    return {g: float(t.change_rate_m_yr[g] - OBS[period][g] / YEARS) for g in E.ENDS}


# Run: the metres solver on position change, then note the target in ends.json
def main():
    E.residuals = residuals
    sys.argv = [sys.argv[0], "--tag", TAG, "--seed", SEED, *OPTION_A, *sys.argv[1:]]
    try:
        code = R.main()
    finally:
        f = R.grid.RAW_RUNS / "experiments" / TAG / "tables" / "ends.json"
        if f.is_file():
            j = json.loads(f.read_text(encoding="utf-8"))
            j["target"] = ("CoastSat end-minus-start position change per window "
                           f"(GIS 1 raw, GIS 90 LOWESS {common.SMOOTH_DOMAINS}); residuals in m/yr = m / {YEARS}")
            j["observed_change_m"] = {str(p): {"1": round(float(OBS[p][1]), 3), "90": round(float(OBS[p][90]), 3)}
                                      for p in OBS}
            f.write_text(json.dumps(j, indent=1), encoding="utf-8")
    return code


if __name__ == "__main__":
    sys.exit(main())
