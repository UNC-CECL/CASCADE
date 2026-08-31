#!/usr/bin/env python3
"""Does the chosen groin hold up THROUGH TIME, not just at the end year?

Every other figure here compares a single end state. That cannot distinguish a
groin that tracks the observations year by year from one that wanders and
happens to arrive in the right place -- and it cannot show WHEN a fit starts to
fail. This one plots the fillet's whole trajectory against the surveys.

WHAT IS PLOTTED
    Three curves per panel:

      observed   the surveyed fillet from `Change_from_wetdry_1967_D2_D12.csv`,
                 re-referenced to the period's start year so it begins at zero
                 like the model does. Markers only on the years actually
                 surveyed -- the record is irregular and joining it with a line
                 would imply samples that do not exist.
      groin ON   the chosen cell, M = 60, f = 0.6.
      groin OFF  the paired M = 0 run at the same be1. The gap between the two
                 model curves is the groin's contribution; the gap from ON to
                 observed is what the source/sink field is left to absorb.

WHY THE MODEL IS RE-REFERENCED TOO
    A run starting in 1984 or 2004 inherits the real fillet in its initial
    shoreline, so its absolute D5-D6 offset is not comparable to a survey
    measured from a 1967 datum. Differencing both sides against their own start
    year removes the inherited part and leaves the CHANGE, which is the only
    quantity the two share.

WHAT THIS FIGURE IS EXPECTED TO SHOW
    A shortfall, and a documented one. M = 60 / f = 0.6 was not chosen to match
    the fillet -- no admissible M can, on this grid -- but by a DIRECT FIT to
    the period-1 D4-D8 change profile (demeaned RMSE 11.69 m against 15.58 m
    with no groin), bounded above by affordability (719,000 m3/yr against a
    5-7e5 littoral drift). The stability bounds an earlier version of this text
    cited -- "M >= 70 unstable, M >= 100 drowns" -- were measured on the
    41-domain rig and do NOT transfer: every production cell through M = 160
    ran clean.
    The residual this figure shows IS the quantity that calibration absorbs,
    together with the Cape Point dynamics the dipole does not represent. Read
    it as the split between the two, not as a failed fit.

Usage:
    python HAT_groin_timeseries_check.py
    python HAT_groin_timeseries_check.py --M 50 --fraction 0.6

Writes output/groin_sweep/figures/timeseries_check.png

Author: Hannah A. Henry, UNC CECL
"""

from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

_HERE = Path(__file__).resolve()
PROJECT_BASE_DIR = _HERE.parents[3]
if not (PROJECT_BASE_DIR / "pyproject.toml").exists():
    raise RuntimeError(
        f"CASCADE repo root not found: {PROJECT_BASE_DIR} has no pyproject.toml.")
SCRIPTS_DIR = PROJECT_BASE_DIR / "scripts"
for _path in (SCRIPTS_DIR, _HERE.parent):
    if str(_path) not in sys.path:
        sys.path.insert(0, str(_path))

from hatteras_site_config import HATTERAS_DOMAINS as GEOMETRY  # noqa: E402

from HAT_groin_sweep_config import (  # noqa: E402
    END_YEAR,
    GROIN_DOWNDRIFT_GIS,
    GROIN_UPDRIFT_GIS,
    PERIODS,
    PRESETS,
    WETDRY_CHANGE_TABLE,
    combo_dir_name,
    sweep_output_dir,
)

FIGURE_DIR = PROJECT_BASE_DIR / "output" / "groin_sweep" / "figures"
OBSERVED_COLOR, ON_COLOR, OFF_COLOR = "#1A1A1A", "#FF8C00", "#777777"

# be1 is swept only in the 1984 edgeBE sweep. -40 is the grid value nearest
# production's -41.8, and is the be1 the period-1 D4-D8 fit was pinned at, so
# this figure and the M = 60 choice rest on the same cell. (An earlier version
# used -34, the reach-RMSE minimiser, which is no longer what M is fitted on.)
EDGE_BE1_1984 = -40.0


def observed_fillet_by_year():
    """Surveyed fillet against the fixed 1967 datum, {year: metres}."""
    frame = pd.read_csv(WETDRY_CHANGE_TABLE).set_index("Domain_ID")
    out = {}
    for column in frame.columns:
        # The SECOND year is the survey year; the first is the 1967 datum.
        match = re.match(r"change_from_wetdry_1967_wetdry_(\d{4})", column)
        if not match:
            continue
        up, down = frame.loc[GROIN_UPDRIFT_GIS, column], frame.loc[GROIN_DOWNDRIFT_GIS, column]
        if pd.isna(up) or pd.isna(down):
            continue
        out.setdefault(int(match.group(1)), []).append(float(down - up))
    return {year: float(np.mean(values)) for year, values in sorted(out.items())}


def model_fillet_series(period, preset, combo):
    """Modelled fillet per year, referenced to the run's own year 0."""
    path = sweep_output_dir(period, preset) / combo / "shoreline_matrix.npy"
    if not path.exists():
        return None, None
    matrix = np.load(path)
    up, down = GEOMETRY.gis_to_pad(GROIN_UPDRIFT_GIS), GEOMETRY.gis_to_pad(GROIN_DOWNDRIFT_GIS)
    offset = matrix[:, down] - matrix[:, up]
    return period + np.arange(matrix.shape[0]), offset - offset[0]


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--M", type=float, default=60.0)
    parser.add_argument("--fraction", type=float, default=0.6)
    args = parser.parse_args()

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    observed = observed_fillet_by_year()
    FIGURE_DIR.mkdir(parents=True, exist_ok=True)
    figure, axes = plt.subplots(len(PERIODS), len(PRESETS), figsize=(14, 8.5),
                               sharex="row")

    for row, period in enumerate(PERIODS):
        end = END_YEAR[period]
        # Re-reference the survey to this period's start: the model begins at
        # zero, so the observations must too.
        in_window = {y: v for y, v in observed.items() if period <= y <= end}
        if not in_window:
            continue
        base_year = min(in_window)
        obs_years = np.array(sorted(in_window))
        obs_vals = np.array([in_window[y] - in_window[base_year] for y in obs_years])

        for col, preset in enumerate(PRESETS):
            axis = axes[row][col]
            be1 = EDGE_BE1_1984 if (preset == "edgeBE" and period == 1984) else None
            on_combo = combo_dir_name(args.M, be1, args.fraction)
            off_combo = combo_dir_name(0.0, be1, 0.0)

            years_on, fil_on = model_fillet_series(period, preset, on_combo)
            years_off, fil_off = model_fillet_series(period, preset, off_combo)

            if years_off is not None:
                axis.plot(years_off, fil_off, color=OFF_COLOR, linestyle=":",
                          linewidth=1.8, label="model, groin OFF", zorder=3)
            if years_on is not None:
                axis.plot(years_on, fil_on, color=ON_COLOR, linewidth=2.4,
                          label=f"model, M={args.M:g} f={args.fraction:g}", zorder=4)
            axis.plot(obs_years, obs_vals, marker="o", markersize=6,
                      linestyle="none", color=OBSERVED_COLOR,
                      label=f"surveyed ({len(obs_years)} dates)", zorder=5)

            axis.axhline(0.0, color="#BBBBBB", linewidth=0.8, zorder=1)
            gap = ""
            if years_on is not None:
                modelled_end = float(fil_on[-1])
                observed_end = float(obs_vals[-1])
                gap = (f"   end: model {modelled_end:+.0f} m vs surveyed "
                       f"{observed_end:+.0f} m  (residual {observed_end - modelled_end:+.0f} m)")
            axis.set_title(f"{period}-{end}  {preset}{gap}", fontsize=10)
            axis.set_ylabel("fillet change since start (m)")
            axis.set_xlabel("year")
            axis.grid(alpha=0.25)
            axis.legend(loc="best", fontsize=8)

    figure.suptitle(
        f"Groin fillet through time -- M={args.M:g}, f={args.fraction:g} "
        f"vs surveyed shoreline", fontsize=13)
    figure.tight_layout(rect=(0, 0.06, 1, 0.965))
    figure.text(
        0.01, 0.01,
        "Both sides are differenced against their own start year: a run beginning in 1984 or 2004 inherits "
        "the real fillet in its initial shoreline, so only the CHANGE is comparable. The gap between the two "
        "model curves is the groin's contribution; the gap from the solid curve to the markers is the residual "
        "left for the source/sink calibration, together with the Cape Point dynamics this dipole cannot represent. "
        "M=60/f=0.6 was chosen by a direct fit to the period-1 D4-D8 change profile, bounded by affordability -- not to "
        "match the fillet, which no "
        "admissible M can on this grid.",
        fontsize=7.5, color="#444444", wrap=True)

    path = FIGURE_DIR / "timeseries_check.png"
    figure.savefig(path, dpi=150, facecolor="white")
    plt.close(figure)
    print(f"wrote {path}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
