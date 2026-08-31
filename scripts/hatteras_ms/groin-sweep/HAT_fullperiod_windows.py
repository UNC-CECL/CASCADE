#!/usr/bin/env python3
"""Does a groin signal survive in a NARROW window? Rescores the sweep several ways.

THE QUESTION
    On the full D1-D12 window the no-groin cell fits best and RMSE rises
    monotonically with trapping. But that window is dominated by two features
    the groin module cannot produce:

        D2-D4   observed ACCRETES (+3 to +7 m over 40 years); the model erodes
                there. Cape Point, which the parameterisation does not
                represent.
        D6-D7   observed erodes 48 and 63 m, a trough peaking one domain NORTH
                of the structure. A groin pushes D6 seaward, i.e. the wrong
                way, so raising M makes the single largest misfit worse.

    Neither is a groin signal, and together they swamp one. The question this
    file answers is whether a groin signal exists underneath, in the few
    domains the structure actually reaches.

TWO SCORES PER WINDOW, AND WHY BOTH ARE NEEDED
    raw         RMSE of modelled against observed change. As the window
                narrows this is increasingly dominated by a constant OFFSET:
                if the model sits 20 m low across D4-D8, raw RMSE is ~20 m
                whatever M does, and the groin's contribution is invisible.

    detrended   A straight line is fitted across the window and removed from
                BOTH profiles first, leaving only SHAPE. A groin's dipole is a
                shape -- a step across the structure -- so this is the score
                that can actually see it. A cell that gets the local shape
                right while sitting on a biased background shows up here and
                nowhere else.

    Reported together on purpose: a groin that improves the detrended score
    while leaving the raw one unchanged has explained the local pattern
    without fixing the level, which is a precise and reportable statement
    rather than a pass or a fail.

WHAT WOULD COUNT AS A SIGNAL
    An INTERIOR optimum at M > 0 that beats the M = 0 baseline. If every
    window at both scores still prefers M = 0, that is a much stronger
    negative result than one window alone -- it says the groin's effect is not
    merely swamped at reach scale but absent at the scale it acts on.

Usage:
    python HAT_fullperiod_windows.py

Reads  output/groin_sweep/fullperiod_1984_2024/results.csv
Writes fit_windows.csv beside it, and prints the comparison.

Author: Hannah A. Henry, UNC CECL
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

_HERE = Path(__file__).resolve()
PROJECT_BASE_DIR = _HERE.parents[3]
for _path in (PROJECT_BASE_DIR / "scripts", _HERE.parent):
    if str(_path) not in sys.path:
        sys.path.insert(0, str(_path))

from HAT_fullperiod_target import observed_change_profile  # noqa: E402

OUT_ROOT = PROJECT_BASE_DIR / "output" / "groin_sweep" / "fullperiod_1984_2024"
RESULTS_CSV = OUT_ROOT / "results.csv"
WINDOWS_CSV = OUT_ROOT / "fit_windows.csv"

# The groin field occupies D6; the pair is D5/D6. Measured influence is 2.25 km
# updrift (~4.5 domains) and ZERO downdrift, so windows widen mostly northward
# in spirit even though they are drawn symmetrically for simplicity.
WINDOWS = {
    "D5-D7  (pair +1)": (5, 7),
    "D4-D8  (pair +/-2)": (4, 8),
    "D3-D9  (pair +/-3)": (3, 9),
    "D1-D12 (reach)": (1, 12),
}


def demean(values):
    """Removes a CONSTANT offset, preserving every gradient and step.

    This is the score that matches how the model is actually applied: a
    per-domain source/sink correction is fitted afterwards, so a uniform level
    error in the groin's neighbourhood is absorbed downstream and is not the
    groin's job to fix. What the groin must get right is the SHAPE.

    Deliberately weaker than `detrend`. Removing a linear trend across a short
    window also removes the gradient that a dipole PRODUCES -- fit a line
    through D4-D8 and the step across D5/D6 is partly absorbed into it, so a
    working groin can be scored as no better than none. Subtracting the mean
    cannot do that.
    """
    values = np.asarray(values, dtype=float)
    return values - values.mean()


def detrend(values, x):
    """Removes a straight line, leaving shape only."""
    values = np.asarray(values, dtype=float)
    x = np.asarray(x, dtype=float)
    if len(values) < 3:
        # Two points define the line exactly, so detrending would zero them
        # and every cell would score identically. Return as-is and let the
        # caller see the raw score instead of a meaningless zero.
        return values - values.mean()
    slope, intercept = np.polyfit(x, values, 1)
    return values - (slope * x + intercept)


def score(frame, observed, lo, hi):
    """raw and detrended RMSE for every cell over one window."""
    domains = [d for d in range(lo, hi + 1) if f"change_D{d}" in frame.columns]
    x = np.array(domains, dtype=float)
    obs = np.array([observed[d] for d in domains], dtype=float)
    obs_shape = detrend(obs, x)

    obs_demeaned = demean(obs)

    raw, centred, shaped = [], [], []
    for _, row in frame.iterrows():
        model = np.array([float(row[f"change_D{d}"]) for d in domains])
        raw.append(float(np.sqrt(np.mean((model - obs) ** 2))))
        centred.append(float(np.sqrt(np.mean((demean(model) - obs_demeaned) ** 2))))
        shaped.append(float(np.sqrt(np.mean((detrend(model, x) - obs_shape) ** 2))))
    out = frame[["combo", "M", "fraction"]].copy()
    out["raw_rmse"] = raw
    out["demeaned_rmse"] = centred        # <- the criterion that matters
    out["detrended_rmse"] = shaped
    return out


def main():
    if not RESULTS_CSV.exists():
        raise SystemExit(f"{RESULTS_CSV} not found -- run the sweep first.")
    frame = pd.read_csv(RESULTS_CSV)
    if frame.empty:
        raise SystemExit("results.csv has no scored cells.")

    collected = []
    print(f"{len(frame)} cells scored over {len(WINDOWS)} windows\n")
    header = f"{'window':<20}{'score':<12}{'best cell':<20}{'RMSE':>8}{'no-groin':>10}{'verdict':>26}"
    print(header)
    print("-" * len(header))

    for label, (lo, hi) in WINDOWS.items():
        scored = score(frame, observed_change_profile(domains=range(lo, hi + 1)), lo, hi)
        scored["window"] = label
        collected.append(scored)
        for column in ("raw_rmse", "demeaned_rmse", "detrended_rmse"):
            ranked = scored.sort_values(column)
            best = ranked.iloc[0]
            baseline = scored[scored.M == 0]
            base = float(baseline.iloc[0][column]) if not baseline.empty else np.nan
            groin_best = scored[scored.M > 0].sort_values(column).iloc[0]
            beats = np.isfinite(base) and groin_best[column] < base
            verdict = (f"groin WINS by {base - groin_best[column]:.2f} m"
                       if beats else "no groin still best")
            name = f"M={best.M:g}, f={best.fraction:.2f}" if best.M > 0 else "M=0 (no groin)"
            print(f"{label:<20}{column.replace('_rmse',''):<12}{name:<20}"
                  f"{best[column]:>8.2f}{base:>10.2f}{verdict:>26}")
        print()

    pd.concat(collected, ignore_index=True).to_csv(WINDOWS_CSV, index=False)
    print(f"written {WINDOWS_CSV}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
