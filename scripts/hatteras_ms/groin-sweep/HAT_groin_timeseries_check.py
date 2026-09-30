#!/usr/bin/env python3
"""
Does the chosen groin hold up through time, not just at the end year?

    python scripts/hatteras_ms/groin-sweep/HAT_groin_timeseries_check.py
    python scripts/hatteras_ms/groin-sweep/HAT_groin_timeseries_check.py --M 50 --fraction 0.6

The fillet year by year against the surveys, both windows. Details: scripts/hatteras_ms/README.md.

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-09-18
"""

from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

_HERE = Path(__file__).resolve()
PROJECT_BASE_DIR = next(_p for _p in _HERE.parents
                        if (_p / "pyproject.toml").exists())
if not (PROJECT_BASE_DIR / "pyproject.toml").exists():
    raise RuntimeError(
        f"CASCADE repo root not found: {PROJECT_BASE_DIR} has no pyproject.toml.")
SCRIPTS_DIR = PROJECT_BASE_DIR / "scripts"
for _path in (SCRIPTS_DIR, _HERE.parent):
    if str(_path) not in sys.path:
        sys.path.insert(0, str(_path))

from site_layer.hatteras_site_config import HATTERAS_DOMAINS as GEOMETRY  # noqa: E402

from site_layer.hat_figure_style import (apply_style, C, INK, INK_MUTED,  # noqa: E402
                              caption, figsize, open_frame, save, _title)
from HAT_groin_sweep_config import (  # noqa: E402
    GROIN_SWEEP_ROOT,
    END_YEAR,
    GROIN_DOWNDRIFT_GIS,
    GROIN_UPDRIFT_GIS,
    PERIODS,
    PRESETS,
    WETDRY_CHANGE_TABLE,
    combo_dir_name,
    sweep_output_dir,
)

# --- CONFIG ------------------------------------------------------------------
FIGURE_DIR = GROIN_SWEEP_ROOT / "figures"
# House colours (2026-09-11)
OBSERVED_COLOR, ON_COLOR, OFF_COLOR = INK, C["ACCENT"], C["BASE"]

# Be1 is swept only in the 1984 edgeBE sweep
EDGE_BE1_1984 = -40.0
# -----------------------------------------------------------------------------


# Surveyed fillet against the fixed 1967 datum, {year
def observed_fillet_by_year():
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


# Modelled fillet per year, referenced to the run's own year 0
def model_fillet_series(period, preset, combo):
    path = sweep_output_dir(period, preset) / combo / "shoreline_matrix.npy"
    if not path.exists():
        return None, None
    matrix = np.load(path)
    up, down = GEOMETRY.gis_to_pad(GROIN_UPDRIFT_GIS), GEOMETRY.gis_to_pad(GROIN_DOWNDRIFT_GIS)
    offset = matrix[:, down] - matrix[:, up]
    return period + np.arange(matrix.shape[0]), offset - offset[0]


# Run: the figure for the chosen pair
def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    parser.add_argument("--M", type=float, default=60.0)
    parser.add_argument("--fraction", type=float, default=0.6)
    args = parser.parse_args()

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    observed = observed_fillet_by_year()
    FIGURE_DIR.mkdir(parents=True, exist_ok=True)
    apply_style()
    figure, axes = plt.subplots(len(PERIODS), len(PRESETS),
                                figsize=figsize("double", aspect=0.62),
                                sharex="row", constrained_layout=True)

    residuals = []
    for row, period in enumerate(PERIODS):
        end = END_YEAR[period]
        # Re-reference the survey to this period's start
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
                          linewidth=1.4, label="modelled, groin off", zorder=3)
            if years_on is not None:
                axis.plot(years_on, fil_on, color=ON_COLOR, linewidth=1.8,
                          label=f"modelled, M {args.M:g}, f {args.fraction:g}",
                          zorder=4)
            axis.plot(obs_years, obs_vals, marker="o", markersize=3.4,
                      linestyle="none", color=OBSERVED_COLOR,
                      label=f"surveyed, {len(obs_years)} dates", zorder=5)

            axis.axhline(0.0, color=INK_MUTED, linewidth=0.8,
                         linestyle=(0, (4, 3)), zorder=1)
            if years_on is not None:
                residuals.append((period, preset, float(fil_on[-1]),
                                  float(obs_vals[-1])))
            _title(axis, row * len(PRESETS) + col,
                   f"{period} to {end}, {preset}")
            axis.set_ylabel("fillet change since start (m)")
            axis.set_xlabel("year")
            # Whole years
            axis.xaxis.set_major_locator(
                plt.matplotlib.ticker.MultipleLocator(5))
            axis.xaxis.set_major_formatter(
                plt.matplotlib.ticker.FuncFormatter(lambda v, _: f"{v:.0f}"))
            axis.grid(axis="y")
            axis.set_axisbelow(True)
            open_frame(axis)

    handles, labels = axes[0][0].get_legend_handles_labels()
    figure.legend(handles, labels, loc="outside lower center", ncol=3,
                  frameon=False)

    ends = "; ".join(
        "{} {}, modelled {:+.0f} m against {:+.0f} m surveyed, residual "
        "{:+.0f} m".format(pr, ps, mod, obs, obs - mod)
        for pr, ps, mod, obs in residuals)
    caption(figure,
            "The modelled fillet through time against the surveys, for both "
            "hindcast periods and both source/sink presets, at M = {M:g} and "
            "f = {f:g}. Both sides are differenced against their own start "
            "year: a run beginning in 1984 or 2004 inherits the real fillet in "
            "its initial shoreline, so only the CHANGE is comparable. The gap "
            "between the two modelled curves is the groin's contribution; the "
            "gap from the solid curve to the markers is the residual left for "
            "the source/sink calibration, together with the Cape Point "
            "dynamics this dipole cannot represent. At the end of each window: "
            "{ends}. The pair was chosen by a direct fit to the period-1 "
            "D4\u2013D8 change profile, bounded by affordability, and NOT to "
            "match the fillet, which no admissible M can on this grid."
            .format(M=args.M, f=args.fraction, ends=ends))

    written = save(figure, FIGURE_DIR / "timeseries_check.png", close=True)
    print(f"wrote {written[0]}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
