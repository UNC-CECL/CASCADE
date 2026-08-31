#!/usr/bin/env python3
"""The observed shoreline-change profile for a continuous 1984-2024 window.

WHY A PROFILE AND NOT A SCALAR
    Every scalar target tried on this calibration produced a RIDGE: many (M, f)
    pairs matched one number equally well, and the ranking then picked whichever
    cell sat at a grid edge. A per-domain CHANGE PROFILE constrains shape as
    well as magnitude, which is what separated the knobs on the 1967 rig -- the
    one configuration that produced an interior optimum in both.

WHY A CONTINUOUS WINDOW
    The fillet BUILDS through the first half of 1984-2024 (+52 m by 2004, from
    the fixed-datum wet/dry surveys) and DECLINES through the second (-76 m by
    2023). A 20-year window sees only one of those, so M and f trade off inside
    it. One 40-year run sees both, and they constrain different combinations.

HOW THE TARGET IS BUILT
    From `groin_analysis_chainage_all.csv` -- 793,259 CoastSat shoreline
    POSITIONS in the window, every one of the 90 domains carrying at least 20.
    For each domain an OLS line is fitted to chainage against decimal year and
    evaluated at both ends; the difference is that domain's change.

    FITTED, NOT DIFFERENCED. Two year-bins differenced give a noisy endpoint
    estimate -- at D1 that route disagreed with the published LRR by a factor
    of two. The OLS fit is also the estimator the rest of the project uses, so
    the target and the model rates are the same kind of number.

SIGN
    Chainage is SEAWARD-POSITIVE, verified by correlating its 1984->2004 rate
    against the published CoastSat LRR (+0.74). Barrier3D's x_s is
    landward-positive, so the model side is negated before comparison.

Author: Hannah A. Henry, UNC CECL
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

_HERE = Path(__file__).resolve()
PROJECT_BASE_DIR = _HERE.parents[3]

CHAINAGE_CSV = (PROJECT_BASE_DIR / "hard-structures" / "groin"
                / "HAT-groin-gis-analysis" / "shoreline_output_coastsat"
                / "groin_analysis_chainage_all.csv")

START_YEAR, END_YEAR = 1984, 2024
FIT_DOMAINS_GIS = tuple(range(1, 13))     # D1-D12, the groin neighbourhood
MIN_OBS_PER_DOMAIN = 20


def observed_change_profile(start=START_YEAR, end=END_YEAR,
                            domains=FIT_DOMAINS_GIS):
    """Observed shoreline change per domain over the window, in metres.

    Args:
        start, end: Window bounds, inclusive.
        domains: GIS domain ids to return.

    Returns:
        A dict of domain -> change in metres, SEAWARD-POSITIVE.

    Raises:
        FileNotFoundError: If the chainage table is absent.
        ValueError: If any requested domain has too few observations.
    """
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


def model_change_profile(shoreline_m, geometry, domains=FIT_DOMAINS_GIS):
    """Modelled shoreline change per domain, in metres, seaward-positive.

    Args:
        shoreline_m: [state x padded domain] array, metres, landward-positive.
        geometry: DomainGeometry for GIS -> pad translation.
        domains: GIS domain ids to return.

    Returns:
        A dict of domain -> change in metres.
    """
    matrix = np.asarray(shoreline_m, dtype=float)
    # x_s increases LANDWARD; negate so + is seaward, matching chainage.
    change = -(matrix[-1] - matrix[0])
    return {int(d): float(change[geometry.gis_to_pad(d)]) for d in domains}


def profile_rmse(model, observed):
    """RMSE between two per-domain change profiles, over shared domains."""
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
