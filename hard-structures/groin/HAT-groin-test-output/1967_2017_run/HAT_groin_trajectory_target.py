#!/usr/bin/env python3
"""The observed groin-fillet TRAJECTORY, 1967-2023, and how to score against it.

WHY A TRAJECTORY AND NOT AN END STATE
    M and f are not separable from a single fillet measurement. A groin with a
    high trapping rate and a low deterioration floor, and one with a low rate
    and a high floor, reach a similar fillet by 2024 along completely different
    paths -- the first builds fast and then relaxes, the second builds slowly
    and holds. Every scalar target tried on the 1984-2004 and 2004-2024 windows
    produced a RIDGE of equally good (M, f) pairs for exactly this reason.

    The path is what separates them, and the path is measured: 24 dated wet/dry
    surveys between 1967 and 2023.

WHY 1967 AND NOT 1984
    The fillet's whole life, measured from `Change_from_wetdry_1967_D2_D12.csv`:

        build    1967-1978     0 -> 117 m
        plateau  1978-2004   117 -> 150 m
        decline  2004-2023   150 ->  74 m

    Both hindcast windows begin AFTER the build. A run starting in 1984
    inherits the fillet in its initial condition rather than predicting it, so
    its fillet CHANGE carries almost no information about M. The 1967 window is
    the only one containing the growth that M controls, while the post-2004
    decline is what constrains f. One run, both knobs, separated by different
    stretches of the same curve.

NO PAIRED BASELINE IS NEEDED
    Elsewhere the modelled fillet is differenced against an M = 0 run, because
    a run starting in 1984 begins with the real fillet already in its initial
    shoreline and the baseline is what removes it. Here the observations are
    themselves changes from a fixed 1967 datum, and the model starts before the
    groin existed, so both sides are already referenced to the same zero. The
    fillet is simply the change in (D5 - D6) since year 0.

SIGN
    The wet/dry table is landward-positive (+ = erosion), matching Barrier3D's
    x_s. The fillet is D5 minus D6: positive when the downdrift domain has
    retreated further than the updrift one, which is what a groin builds.

Author: Hannah A. Henry, UNC CECL
"""

from __future__ import annotations

import pathlib
import re

import numpy as np
import pandas as pd

PROJECT_BASE_DIR = next(
    p for p in pathlib.Path(__file__).resolve().parents
    if (p / "pyproject.toml").exists())

WETDRY_CHANGE_TABLE = (
    PROJECT_BASE_DIR / "hard-structures" / "groin" / "HAT-groin-test-output"
    / "shoreline_position_output" / "Change_from_wetdry_1967_D2_D12.csv")

DATUM_YEAR = 1967
UPDRIFT_GIS = 6
DOWNDRIFT_GIS = 5


def observed_trajectory(table=WETDRY_CHANGE_TABLE):
    """Observed fillet against the fixed 1967 datum, per dated survey.

    Args:
        table: The wet/dry change table.

    Returns:
        A DataFrame indexed by year with a `fillet_m` column, sorted by year.

    Raises:
        FileNotFoundError: If the table is absent.
        ValueError: If neither groin domain has any dated column.
    """
    if not table.exists():
        raise FileNotFoundError(
            f"wet/dry change table not found at {table}. It is produced by the "
            f"shoreline-position prep in HAT-groin-test-input/input_prep/.")
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
    # The datum year is 0 by construction; carry it explicitly so a model
    # trajectory that starts at 0 is compared against a target that does too.
    if DATUM_YEAR not in series.index:
        series.loc[DATUM_YEAR] = 0.0
        series = series.sort_index()
    return series


def model_trajectory(shoreline_m, start_year, updrift_pad, downdrift_pad):
    """Modelled fillet per simulated year, referenced to year 0.

    Args:
        shoreline_m: [state x padded domain] array, metres, landward-positive.
        start_year: Calendar year of state 0.
        updrift_pad, downdrift_pad: Padded indices of the groin pair.

    Returns:
        A DataFrame indexed by year with a `fillet_m` column.
    """
    matrix = np.asarray(shoreline_m, dtype=float)
    diff = matrix[:, downdrift_pad] - matrix[:, updrift_pad]
    fillet = diff - diff[0]
    years = start_year + np.arange(matrix.shape[0])
    return pd.DataFrame({"year": years, "fillet_m": fillet}).set_index("year")


def score(model, observed):
    """RMSE between a modelled and observed trajectory, on shared years.

    Only years the survey actually sampled are compared -- the model has a
    value every year, the observations do not, and interpolating the survey
    onto the model's grid would invent data and weight the well-sampled
    2010s the same as the sparse 1970s.

    Args:
        model: From `model_trajectory`.
        observed: From `observed_trajectory`.

    Returns:
        A dict with rmse_m, n_years, bias_m, and the paired frame.
    """
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
