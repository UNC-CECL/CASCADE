#!/usr/bin/env python3
"""
Score the option A groin grid against the observed Buxton fillet change.

    python score_groin_grid.py

Two readings per cell, against observed_fillet_m(period): groin_contribution_m
(end-year gap of the groin run minus its baseline, the M = 60 fit's metric) and
total_change_m (the run's own OLS gap change, built like the observation).
M = 0 is the adopted matrix no-groin run. Writes grid_scores.csv.
Details: README.md beside this script.

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-10-01
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
REPO = next(p for p in Path(__file__).resolve().parents if (p / "pyproject.toml").exists())
sys.path[:0] = [str(REPO / "scripts" / "hatteras_ms" / "groin-sweep"),
                str(REPO / "scripts" / "hatteras_ms"),
                str(REPO / "scripts" / "site_layer"), str(REPO / "scripts")]
from HAT_groin_sweep_config import observed_fillet_m  # noqa: E402

# --- CONFIG ------------------------------------------------------------------
RAW = REPO / "output" / "raw_runs"
GRID = RAW / "experiments" / "groin" / "2026-09-29-option-a-grid"
BASELINES = {
    1996: RAW / "matrix/1996_2010/edgeBE/HAT_1996_2010_edgeBE_offsetmetres_road_bdm_nogroin",
    2010: RAW / "matrix/2010_2024/edgeBE/HAT_2010_2024_edgeBE_offsetmetres_road_bdm_nourish_nogroin",
}
UP, DOWN = 15 + 6 - 1, 15 + 5 - 1      # GIS 6 / GIS 5 in the padded array
# -----------------------------------------------------------------------------


# A run's shoreline matrix (years x padded domains)
def load(run_dir):
    return np.load(run_dir / f"{run_dir.name}_shoreline_matrix.npy")


# D5 - D6 gap per year, landward-positive
def fillet_series(x):
    return x[:, DOWN] - x[:, UP]


# Gap change over the run, OLS slope x run length
def total_change(x):
    s = fillet_series(x)
    t = np.arange(s.size)
    return float(np.polyfit(t, s, 1)[0] * (s.size - 1))


# Run: score the baselines and every grid cell, write the CSV, print the tables
def main():
    rows = []
    for period, base_dir in BASELINES.items():
        base = load(base_dir)
        target = observed_fillet_m(period)
        rows.append(dict(period=period, M=0.0, f=np.nan,
                         groin_contribution_m=0.0,
                         total_change_m=total_change(base),
                         observed_change_m=target))
        for run_dir in sorted(GRID.glob(f"M*_f*/{period}_*/edgeBE/*")):
            member = run_dir.parents[2].name
            M = float(member.split("_")[0][1:])
            f = float(member.split("_f")[1])
            x = load(run_dir)
            rows.append(dict(
                period=period, M=M, f=f,
                groin_contribution_m=float(fillet_series(x)[-1] - fillet_series(base)[-1]),
                total_change_m=total_change(x),
                observed_change_m=target))
    t = pd.DataFrame(rows)
    t["err_total_m"] = t.total_change_m - t.observed_change_m
    t["err_contribution_m"] = t.groin_contribution_m - t.observed_change_m
    t.to_csv(HERE / "grid_scores.csv", index=False)

    pd.set_option("display.width", 200)
    for period, sub in t.groupby("period"):
        print(f"\n=== {period}: observed fillet change {sub.observed_change_m.iloc[0]:+.1f} m ===")
        print(f"no-groin total change {sub.loc[sub.M == 0, 'total_change_m'].iloc[0]:+.1f} m")
        for col in ("total_change_m", "groin_contribution_m"):
            piv = sub[sub.M > 0].pivot(index="M", columns="f", values=col)
            print(f"\n{col} (rows M, cols f)")
            print(piv.to_string(float_format="%+.1f"))
        best = sub.loc[sub.err_total_m.abs().idxmin()]
        print(f"\nbest on total change: M={best.M:g} f={best.f} "
              f"-> {best.total_change_m:+.1f} m (err {best.err_total_m:+.1f})")


if __name__ == "__main__":
    main()
