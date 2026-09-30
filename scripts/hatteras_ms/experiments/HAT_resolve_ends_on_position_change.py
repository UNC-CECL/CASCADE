"""Solve the two end domains on position change, not LRR (2026-09-28).

Hannah, 2026-09-28: "calculate the edge source sink values we should be
matching based on the position change instead" -- option 1 of three: an
experiment, the solver and the config unchanged. The LRR-solved ends match
the rates to 0.15 m/yr but overshoot the observed end-minus-start change at
GIS 1 1996-2010 (+38 m) and GIS 90 2010-2024 (+21 m), because the observed
ends do not move linearly and the model does.

    target     the CoastSat total change, mean of the end year's images minus
               mean of the start year's (5-scr/3-rates/coastsat/total_change/
               <window>/smoothed/tables/domain_smoothed.csv): GIS 1 raw, GIS 90
               LOWESS over common.SMOOTH_DOMAINS (7) -- the same treatment the
               LRR target gives each end (Hannah: "use the smoothed value at
               GIS 90"; built by HAT_wave_recommendation_figures.
               observed_change_smoothed)
    model      end-minus-start position, the run's change_rate_m_yr x 14
    residual   (model - observed) / 14 years, in m/yr, so HAT_resolve_ends_
               metres' solver, gains and 0.02 m/yr tolerance (0.28 m of
               change) apply as they are
    waves      option A (Hs 2.0, Tp 7.5, asymmetry 0.6, high-angle 0.5),
               full management, as the LRR solve
    seed       the LOWESS-7 LRR-solved ends (2026-09-28-ends-resolved-lowess7)

WHERE: output/raw_runs/experiments/end-domain-boundaries/2026-09-28-ends-solved-on-position-change/

    python scripts/hatteras_ms/experiments/HAT_resolve_ends_on_position_change.py
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
TAG = "end-domain-boundaries/2026-09-28-ends-solved-on-position-change"
YEARS = 14
SEED = "1996=4.8394,18.2545;2010=18.8657,24.2358"
OPTION_A = ["--hs", "2.0", "--tp", "7.5", "--asym", "0.6", "--ahf", "0.5"]
OBS = {p: W.observed_change_smoothed(p) for p in (1996, 2010)}


def residuals(run_dir, period, targets):
    """(model - observed) end-minus-start change, per year. `targets` (the
    LRR target R.main builds) is not used."""
    t = pd.read_csv(Path(run_dir) / "tables" / "shoreline_change_rate.csv").set_index("gis_domain")
    return {g: float(t.change_rate_m_yr[g] - OBS[period][g] / YEARS) for g in E.ENDS}


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
