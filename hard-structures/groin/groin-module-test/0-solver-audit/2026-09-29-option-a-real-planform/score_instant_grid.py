#!/usr/bin/env python3
"""Score the full-model groin grids (instant 2004 failure) on the observed gap.

For every run under raw_runs/experiments/groin/2026-09-29-instant-2004-grid/
(blocking b<b>_f<f>, dipole M<M>_f<f>) and the matrix no-groin baseline:

  date RMSE    the D5 - D6 gap change since the window start, sampled at the
               wet/dry table's own dates inside the window, against the
               observed change (2010 start interpolated between 2008 and 2014).
               The primary score: an OLS trend alone is matched by
               build-then-collapse shapes the data do not show.
  OLS change   trend x 14 yr, as observed_fillet_m builds its target.

Both windows, both kinds; joint = RMS of the two windows' date RMSEs.
Sign: landward-positive, + = downdrift sits further landward of updrift.

    python score_instant_grid.py  ->  instant_grid_scores.csv + printed tables

Author: Hannah A. Henry, UNC CECL
"""

from __future__ import annotations

import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[4]
sys.path[:0] = [str(REPO / "scripts" / "hatteras_ms" / "groin-sweep"),
                str(REPO / "scripts" / "hatteras_ms"),
                str(REPO / "scripts" / "site_layer"), str(REPO / "scripts")]
from HAT_groin_sweep_config import WETDRY_CHANGE_TABLE, observed_fillet_m  # noqa: E402

RAW = REPO / "output" / "raw_runs"
GRID = RAW / "experiments" / "groin" / "2026-09-29-instant-2004-grid"
BASELINES = {
    1996: RAW / "matrix/1996_2010/edgeBE/HAT_1996_2010_edgeBE_offsetmetres_road_bdm_nogroin",
    2010: RAW / "matrix/2010_2024/edgeBE/HAT_2010_2024_edgeBE_offsetmetres_road_bdm_nourish_nogroin",
}
UP, DOWN = 15 + 6 - 1, 15 + 5 - 1
YEARS = 14


def observed_series():
    tab = pd.read_csv(WETDRY_CHANGE_TABLE).set_index("Domain_ID")
    obs = {}
    for c in tab.columns:
        m = re.match(r"change_from_wetdry_1967_wetdry_(\d{4})", c)
        if m:
            obs.setdefault(int(m.group(1)), []).append(tab.loc[5, c] - tab.loc[6, c])
    return pd.Series({y: np.mean(v) for y, v in obs.items()}).sort_index()


def score(x, start, obs):
    g = x[:, DOWN] - x[:, UP]
    g = g - g[0]
    yrs = [y for y in obs.index if start < y <= start + YEARS]
    g0 = np.interp(start, obs.index, obs.values)
    o = np.array([obs[y] for y in yrs]) - g0
    m = np.array([g[y - start] for y in yrs])
    t = np.arange(g.size)
    return dict(rmse_dates=float(np.sqrt(np.mean((m - o) ** 2))),
                ols_change=float(np.polyfit(t, g, 1)[0] * YEARS),
                observed_ols=observed_fillet_m(start),
                model_at_dates=" ".join(f"{v:+.0f}" for v in m),
                observed_at_dates=" ".join(f"{v:+.0f}" for v in o))


def main():
    obs = observed_series()
    rows = []
    for start, bdir in BASELINES.items():
        x = np.load(bdir / f"{bdir.name}_shoreline_matrix.npy")
        rows.append(dict(kind="none", strength=0.0, f=np.nan, start=start,
                         **score(x, start, obs)))
        for f in sorted(GRID.glob(f"*_f*/{start}_*/edgeBE/*/*_shoreline_matrix.npy")):
            member = f.parents[3].name
            kind = "blocking" if member.startswith("b") else "dipole"
            strength = float(member.split("_f")[0][1:])
            rows.append(dict(kind=kind, strength=strength,
                             f=float(member.split("_f")[1]), start=start,
                             **score(np.load(f), start, obs)))
    t = pd.DataFrame(rows)
    t.to_csv(HERE / "instant_grid_scores.csv", index=False)

    pd.set_option("display.width", 220, "display.max_colwidth", 60)
    counts = t.groupby(["kind", "start"]).size()
    print("runs scored:\n" + counts.to_string())
    for start in BASELINES:
        sub = t[t.start == start]
        print(f"\n=== {start}: observed at dates {sub.observed_at_dates.iloc[0]} "
              f"(OLS {sub.observed_ols.iloc[0]:+.1f}) ===")
        base = sub[sub.kind == "none"].iloc[0]
        print(f"no groin: {base.model_at_dates}  rmse {base.rmse_dates:.1f}")
        for kind in ("blocking", "dipole"):
            k = sub[sub.kind == kind]
            if len(k):
                print(f"\n{kind} date-RMSE (rows strength, cols f)")
                print(k.pivot(index="strength", columns="f", values="rmse_dates")
                      .to_string(float_format="%.1f"))
    j = t[t.kind != "none"].pivot_table(index=["kind", "strength", "f"],
                                        columns="start", values="rmse_dates")
    base = t[t.kind == "none"].set_index("start").rmse_dates
    if {1996, 2010} <= set(j.columns):
        j = j.dropna()
        j["joint"] = np.sqrt((j[1996] ** 2 + j[2010] ** 2) / 2)
        print(f"\nno groin joint: {np.sqrt((base[1996]**2 + base[2010]**2) / 2):.1f}")
        for kind in ("blocking", "dipole"):
            if kind in j.index.get_level_values(0):
                print(f"\n{kind}: best joint date-RMSE (m)")
                print(j.loc[kind].sort_values("joint").head(6).round(1).to_string())


if __name__ == "__main__":
    main()
