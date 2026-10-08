#!/usr/bin/env python3
"""
The observed shoreline-change profile for a continuous 1984-2024 window.

    python scripts/hatteras_ms/groin-sweep/HAT_fullperiod_target.py

A profile, not a scalar: every scalar target produced a ridge of equal fits. Details: scripts/hatteras_ms/README.md.

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-09-14
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

_HERE = Path(__file__).resolve()
PROJECT_BASE_DIR = next(_p for _p in _HERE.parents
                        if (_p / "pyproject.toml").exists())
# --- CONFIG ------------------------------------------------------------------
CHAINAGE_CSV = (PROJECT_BASE_DIR / "hard-structures" / "groin"
                / "1-observations" / "coastsat_shoreline" / "shoreline_output_coastsat"
                / "groin_analysis_chainage_all.csv")

START_YEAR, END_YEAR = 1984, 2024
FIT_DOMAINS_GIS = tuple(range(1, 13))     # D1-D12, the groin neighbourhood
MIN_OBS_PER_DOMAIN = 20
# -----------------------------------------------------------------------------


# Observed shoreline change per domain over the window, in metres
def observed_change_profile(start=START_YEAR, end=END_YEAR,
                            domains=FIT_DOMAINS_GIS):
    if not CHAINAGE_CSV.exists():
        raise FileNotFoundError(
            f"CoastSat chainage not found at {CHAINAGE_CSV}; it is produced by "
            f"HAT_groin_shoreline_analysis_v2.py.")
    frame = pd.read_csv(CHAINAGE_CSV,
                        usecols=["chainage_m", "source", "domain", "decimal_year"])
    frame = frame[(frame["source"] == "coastsat")
                  & (frame["decimal_year"] >= start)
                  & (frame["decimal_year"] <= end)]

    out, thin = {}, []
    for domain in domains:
        rows = frame[frame["domain"] == domain]
        if len(rows) < MIN_OBS_PER_DOMAIN:
            thin.append(domain)
            continue
        slope, _ = np.polyfit(rows["decimal_year"], rows["chainage_m"], 1)
        out[int(domain)] = float(slope * (end - start))
    if thin:
        raise ValueError(
            f"domains {thin} have fewer than {MIN_OBS_PER_DOMAIN} CoastSat "
            f"observations in {start}-{end}; the profile would be undefined "
            f"there.")
    return out


# Modelled shoreline change per domain, in metres, seaward-positive
def model_change_profile(shoreline_m, geometry, domains=FIT_DOMAINS_GIS):
    matrix = np.asarray(shoreline_m, dtype=float)
    # x_s increases LANDWARD; negate so + is seaward, matching chainage.
    change = -(matrix[-1] - matrix[0])
    return {int(d): float(change[geometry.gis_to_pad(d)]) for d in domains}


# RMSE between two per-domain change profiles, over shared domains
def profile_rmse(model, observed):
    shared = sorted(set(model) & set(observed))
    if not shared:
        raise ValueError("model and observed profiles share no domains")
    errors = np.array([model[d] - observed[d] for d in shared], dtype=float)
    return float(np.sqrt(np.mean(errors ** 2)))


if __name__ == "__main__":
    obs = observed_change_profile()
    print(f"observed change {START_YEAR}->{END_YEAR}, seaward-positive (m)\n")
    for d in sorted(obs):
        print(f"  D{d:<3} {obs[d]:+8.1f}")
    values = np.array(list(obs.values()))
    print(f"\n  range {values.min():+.1f} to {values.max():+.1f} m, "
          f"mean {values.mean():+.1f}")
