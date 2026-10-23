"""
How well does each dune setting reproduce where Isabel overwashed, beside how well it follows the surveyed dunes?

    python scripts/hatteras_ms/experiments/HAT_isabel_pattern_scores.py

Re-scores the 1996-2009 runs of the ceiling x growth grid and the low growth-rate
study (no new runs). Isabel: model overwash in model year 2003, per GIS domain,
against the first image after Isabel (2004-05-25, overwash in 48 of 90 domains).
The imagery is about one image a year, so a model overwash the image does not
show is weak evidence; the hits on the major event are the strong evidence
(Hannah, 2026-10-19). Dune trajectory: the low growth-rate study's RMSE of the
low-dune 10th percentile and of the share of cells below 4.3 m over the eight
calibration surveys. Read-only.

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
from site_layer import hat_overwash as ow  # noqa: E402
from site_layer.hatteras_site_config import HATTERAS_DOMAINS as DOM  # noqa: E402

# --- CONFIG ------------------------------------------------------------------
TAG = "storms-overwash/2026-10-19-isabel-pattern-scores"
EXP_DIR = U.PROJECT_ROOT / "output" / "raw_runs" / "experiments" / TAG
STUDIES = {"grid": "dunes/2026-10-19-dune-ceiling-growth-grid",
           "low rates": "dunes/2026-10-19-dune-growth-low-rates"}
EVENT_YEAR, FIRST_YEAR = 2003, 1997      # Isabel; QowTS row t is model year FIRST_YEAR + t - 1
IMAGE = "2004-05-25"                     # the first image after Isabel
# -----------------------------------------------------------------------------


def observed():
    o = pd.read_csv(ow.OBSERVATIONS / "overwash_observations.csv")
    return o[o.Imagery_Date == IMAGE].set_index("domain").overwash.reindex(range(1, 91))


def isabel(c, obs):
    t = EVENT_YEAR - FIRST_YEAR + 1
    q = np.array([float(np.asarray(c.barrier3d[DOM.gis_to_pad(g)].QowTS)[t]) for g in range(1, 91)])
    m = pd.Series((q > 0).astype(int), index=range(1, 91))
    ok = obs.notna()
    o, m = obs[ok].astype(int), m[ok]
    hits, miss = int(((o == 1) & (m == 1)).sum()), int(((o == 1) & (m == 0)).sum())
    fa, cn = int(((o == 0) & (m == 1)).sum()), int(((o == 0) & (m == 0)).sum())
    pod, pofd = hits / max(hits + miss, 1), fa / max(fa + cn, 1)
    return dict(isabel_hits=hits, isabel_misses=miss, isabel_false_alarms=fa,
                isabel_hit_rate=round(pod, 3), isabel_skill=round(pod - pofd, 3),
                isabel_domains_modelled=int(m.sum()), isabel_agreement=round((hits + cn) / len(o), 3))


def main():
    obs = observed()
    rows = []
    for study, tag in STUDIES.items():
        G.TAG, G.EXP_DIR = tag, U.PROJECT_ROOT / "output" / "raw_runs" / "experiments" / tag
        base = G.EXP_DIR / "runs"
        for d in sorted(base.glob("*_s0")):
            variant = d.name
            if G.run_dir(1996, variant) is None:
                continue
            ceiling, rates, _, _ = G.parse(variant)
            c = U.load(G.run_dir(1996, variant))
            row = dict(study=study, variant=variant, ceiling_m_navd88=ceiling, rmin=rates[0], rmax=rates[1])
            row.update(isabel(c, obs))
            row.update({k: round(v, 3) for k, v in L.trajectory(1996, variant).items()})
            m = G.metrics(1996, variant)
            row.update(overall_skill=m["skill"], overall_false_alarms=m["false_alarms"],
                       shoreline_rmse=round(m["rmse_interior_m_yr"], 4))
            rows.append(row)
    t = pd.DataFrame(rows).drop_duplicates(subset=["variant"]).sort_values(["rmin", "ceiling_m_navd88"])
    (EXP_DIR / "tables").mkdir(parents=True, exist_ok=True)
    t.to_csv(EXP_DIR / "tables" / "isabel_and_dune_trajectory.csv", index=False)
    cols = ["variant", "isabel_hits", "isabel_misses", "isabel_false_alarms", "isabel_skill", "p10_rmse",
            "below_rmse", "domain_r_mean", "overall_skill", "shoreline_rmse"]
    with pd.option_context("display.width", 250, "display.max_columns", 30):
        print(f"Isabel image {IMAGE}: overwash observed in {int(obs.sum())} of {int(obs.notna().sum())} domains")
        print(t[cols].to_string(index=False))


if __name__ == "__main__":
    main()
