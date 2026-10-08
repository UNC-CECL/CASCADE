#!/usr/bin/env python3
"""
The observed groin-fillet trajectory, 1967-2023, and how to score a model run against it.

    python HAT_groin_trajectory_target.py

Reads Change_from_wetdry_1967_D2_D12.csv; the fillet is D5 minus D6 change
from the fixed 1967 datum, landward-positive. Run directly it prints the
trajectory and its build / plateau / decline phases; import it for
observed_trajectory(), model_trajectory() and score().
Details: README.md beside this script.

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-10-01
"""

from __future__ import annotations

import pathlib
import re

import numpy as np
import pandas as pd


# --- CONFIG ------------------------------------------------------------------
# Repo root, found by searching upward (ORGANIZATION.md rule 5)
PROJECT_BASE_DIR = next(
    p for p in pathlib.Path(__file__).resolve().parents
    if (p / "pyproject.toml").exists())

WETDRY_CHANGE_TABLE = (
    PROJECT_BASE_DIR / "hard-structures" / "groin" / "1-observations"
    / "wetdry_photo_positions" / "Change_from_wetdry_1967_D2_D12.csv")

DATUM_YEAR = 1967
UPDRIFT_GIS = 6
DOWNDRIFT_GIS = 5
# -----------------------------------------------------------------------------


# Observed fillet per dated survey against the 1967 datum: year -> fillet_m
def observed_trajectory(table=WETDRY_CHANGE_TABLE):
    if not table.exists():
        raise FileNotFoundError(
            f"wet/dry change table not found at {table}. It is produced by the "
            f"shoreline-position prep in 3-hindcast/1-dipole-1967-2017/inputs/input_prep/.")
    frame = pd.read_csv(table).set_index("Domain_ID")

    rows = {}
    for column in frame.columns:
        match = re.match(r"change_from_wetdry_1967_wetdry_(\d{4})", column)
        if not match:
            continue
        up, down = frame.loc[UPDRIFT_GIS, column], frame.loc[DOWNDRIFT_GIS, column]
        if pd.isna(up) or pd.isna(down):
            continue
        rows.setdefault(int(match.group(1)), []).append(float(down - up))
    if not rows:
        raise ValueError(
            f"{table.name} has no dated wetdry column with both D{UPDRIFT_GIS} "
            f"and D{DOWNDRIFT_GIS} populated.")

    series = pd.DataFrame(
        {"year": sorted(rows), "fillet_m": [float(np.mean(rows[y])) for y in sorted(rows)]}
    ).set_index("year")
    # The datum year is 0 by construction; carried so both trajectories start at 0
    if DATUM_YEAR not in series.index:
        series.loc[DATUM_YEAR] = 0.0
        series = series.sort_index()
    return series


# Modelled fillet per simulated year from a [state x padded domain] matrix, referenced to year 0
def model_trajectory(shoreline_m, start_year, updrift_pad, downdrift_pad):
    matrix = np.asarray(shoreline_m, dtype=float)
    diff = matrix[:, downdrift_pad] - matrix[:, updrift_pad]
    fillet = diff - diff[0]
    years = start_year + np.arange(matrix.shape[0])
    return pd.DataFrame({"year": years, "fillet_m": fillet}).set_index("year")


# RMSE and bias on surveyed years only, never interpolated (README)
def score(model, observed):
    joined = observed.join(model, how="inner", lsuffix="_obs", rsuffix="_mod")
    if joined.empty:
        raise ValueError(
            "model and observed trajectories share no years; check start_year "
            "and the run length.")
    residual = joined["fillet_m_mod"] - joined["fillet_m_obs"]
    return dict(
        rmse_m=float(np.sqrt(np.mean(residual ** 2))),
        bias_m=float(residual.mean()),
        n_years=int(len(joined)),
        paired=joined,
    )


if __name__ == "__main__":
    obs = observed_trajectory()
    print(f"observed fillet trajectory -- {len(obs)} dated surveys "
          f"{obs.index.min()}-{obs.index.max()}\n")
    for year, value in obs["fillet_m"].items():
        print(f"  {year}  {value:8.1f} m")
    print("\nphases:")
    for a, b, label in ((1967, 1978, "build  "), (1978, 2004, "plateau"),
                        (2004, 2023, "decline")):
        va, vb = obs["fillet_m"].get(a), obs["fillet_m"].get(b)
        if va is not None and vb is not None:
            print(f"  {label} {a}-{b}: {va:6.1f} -> {vb:6.1f} m  ({vb - va:+.1f})")
