"""
How much do the overwash scores depend on how strictly modelled overwash is matched to the imagery?

    python scripts/hatteras_ms/experiments/HAT_overwash_scoring_tolerance.py

Re-scores existing runs (no new runs) under four matching rules, from the
current one to a tolerant one, and checks each rule still separates a model
from the same output shifted to the wrong domains. Details:
scripts/hatteras_ms/experiments/README.md.

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-10-13
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

_HERE = Path(__file__).resolve()
sys.path.insert(0, str(_HERE.parent))
import HAT_dune_ceiling_uniform as U  # noqa: E402

# --- CONFIG ------------------------------------------------------------------
OUT = U.PROJECT_ROOT / "output" / "raw_runs" / "experiments" / "storms-overwash-and-dunes" / \
    "2026-10-13-overwash-scoring-tolerance"
EXP = U.PROJECT_ROOT / "output" / "raw_runs" / "experiments" / "storms-overwash-and-dunes"
RUNS = {"single 5.5 m (adopted)": "2026-10-13-dune-ceiling-range/runs/c550",
        "per-cell (fork)": "2026-10-13-dune-ceiling-uniform/runs/percell",
        "per-domain median": "2026-10-13-dune-ceiling-domain-value/runs/dom50"}
RULES = {"strict (current)": dict(any_storm=False, dt=0, dx=0),
         "any storm in the model year": dict(any_storm=True, dt=0, dx=0),
         "+ one image either side": dict(any_storm=True, dt=1, dx=0),
         "+ one domain either side": dict(any_storm=True, dt=1, dx=1)}
SHIFTS = range(10, 81, 10)          # alongshore offsets, in domains, for the wrong-place check
# -----------------------------------------------------------------------------


# Observed and modelled overwash as (image, GIS 1-90) grids; model: per image, did a model year
# whose storm (largest only, or any) is dated in that image's window overwash the domain
def grids(c, obs, dom, storm_years, any_storm):
    import overwash_vs_model as ovm
    from site_layer import hat_env_forcings as env
    q = np.array([np.asarray(b.QowTS) for b in c.barrier3d]).T
    s = pd.read_csv(env.storm_summary_file(*storm_years), parse_dates=["StartTime", "EndTime"])
    s["dated"] = s.EndTime - ovm.GRACE
    largest = s.loc[s.groupby("time").Rhigh.idxmax()].set_index("time").dated
    every = s.groupby("time").dated.apply(list)
    imgs = ovm.images_in((storm_years[0], storm_years[1] + 1), obs)
    O = np.full((len(imgs), 90), np.nan)
    M = np.zeros((len(imgs), 90))
    for i, im in imgs.iterrows():
        o = obs[obs.Obs_ID == im.Obs_ID].set_index("domain").overwash
        for gis in range(1, 91):
            O[i, gis - 1] = o.get(gis, np.nan)
            p = dom.gis_to_pad(gis)
            yrs = [t for t in range(1, q.shape[0]) if q[t, p] > ovm.HEADLINE_THR]
            if any_storm:
                M[i, gis - 1] = any(im["from"] < d <= im.date for t in yrs if t in every for d in every[t])
            else:
                M[i, gis - 1] = any(im["from"] < largest[t] <= im.date for t in yrs if t in largest)
    return O, M


# Widen a 0/1 grid by dt images and dx domains
def widen(G, dt, dx):
    W = np.nan_to_num(G).copy()
    out = W.copy()
    for a in range(-dt, dt + 1):
        for b in range(-dx, dx + 1):
            out = np.maximum(out, np.roll(np.roll(W, a, axis=0), b, axis=1) *
                             _valid(W.shape, a, b))
    return out


def _valid(shape, a, b):
    v = np.ones(shape)
    if a > 0:
        v[:a] = 0
    elif a < 0:
        v[a:] = 0
    if b > 0:
        v[:, :b] = 0
    elif b < 0:
        v[:, b:] = 0
    return v


# Symmetric tolerance: a hit is observed overwash with modelled overwash nearby; a false alarm is
# modelled overwash with no observed overwash nearby; skill = hit rate - false-alarm rate
def score(O, M, dt, dx):
    ok = ~np.isnan(O)
    o, m = np.nan_to_num(O) == 1, M == 1
    m_near, o_near = widen(M, dt, dx) == 1, widen(np.where(ok, O, 0), dt, dx) == 1
    hits = (ok & o & m_near).sum()
    misses = (ok & o & ~m_near).sum()
    fa = (ok & m & ~o_near).sum()
    cn = (ok & ~m & ~o).sum()
    pod, pofd = hits / max(hits + misses, 1), fa / max(fa + cn, 1)
    return pod, pofd, pod - pofd


def main():
    from site_layer import hat_overwash as ow
    from site_layer.hatteras_site_config import HATTERAS_DOMAINS as DOM
    obs = pd.read_csv(ow.OBSERVATIONS / "overwash_observations.csv")
    rows = []
    for name, rel in RUNS.items():
        for start, (end, storm_years) in U.PERIODS.items():
            rd = sorted((EXP / rel / f"{start}_{end}" / U.PRESET).glob("HAT_*"))[-1]
            c = U.load(rd)
            cache = {a: grids(c, obs, DOM, storm_years, a) for a in (False, True)}
            for rule, r in RULES.items():
                O, M = cache[r["any_storm"]]
                pod, pofd, pss = score(O, M, r["dt"], r["dx"])
                wrong = [score(O, np.roll(M, k, axis=1), r["dt"], r["dx"])[2] for k in SHIFTS]
                rows.append(dict(run=name, period=f"{start}_{end}", rule=rule, hit_rate=round(pod, 3),
                                 false_alarm_rate=round(pofd, 3), skill=round(pss, 3),
                                 skill_wrong_place=round(float(np.mean(wrong)), 3),
                                 skill_above_wrong_place=round(pss - float(np.mean(wrong)), 3)))
            print(f"  scored {name} {start}", flush=True)
    t = pd.DataFrame(rows)
    OUT.mkdir(parents=True, exist_ok=True)
    t.to_csv(OUT / "scores.csv", index=False)
    with pd.option_context("display.width", 250, "display.max_columns", 20):
        print(t.to_string(index=False))


if __name__ == "__main__":
    main()
