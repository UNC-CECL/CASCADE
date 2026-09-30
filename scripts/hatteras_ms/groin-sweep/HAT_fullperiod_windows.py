#!/usr/bin/env python3
"""
Does a groin signal survive in a narrow window? The full-period sweep rescored several ways.

    python scripts/hatteras_ms/groin-sweep/HAT_fullperiod_windows.py

Writes the windows table. Details: scripts/hatteras_ms/README.md.

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-09-18
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

_HERE = Path(__file__).resolve()
PROJECT_BASE_DIR = next(_p for _p in _HERE.parents
                        if (_p / "pyproject.toml").exists())
for _path in (PROJECT_BASE_DIR / "scripts", _HERE.parent):
    if str(_path) not in sys.path:
        sys.path.insert(0, str(_path))

from HAT_groin_sweep_config import GROIN_SWEEP_ROOT  # noqa: E402
from HAT_fullperiod_target import observed_change_profile  # noqa: E402

# --- CONFIG ------------------------------------------------------------------
OUT_ROOT = GROIN_SWEEP_ROOT / "fullperiod_1984_2024"
RESULTS_CSV = OUT_ROOT / "results.csv"
WINDOWS_CSV = OUT_ROOT / "fit_windows.csv"

# The groin field is D6 (pair D5/D6); influence 2.25 km updrift, none downdrift
WINDOWS = {
    "D5-D7  (pair +1)": (5, 7),
    "D4-D8  (pair +/-2)": (4, 8),
    "D3-D9  (pair +/-3)": (3, 9),
    "D1-D12 (reach)": (1, 12),
}
# -----------------------------------------------------------------------------


# Removes a CONSTANT offset, preserving every gradient and step
def demean(values):
    values = np.asarray(values, dtype=float)
    return values - values.mean()


# Removes a straight line, leaving shape only
def detrend(values, x):
    values = np.asarray(values, dtype=float)
    x = np.asarray(x, dtype=float)
    if len(values) < 3:
        # Two points define the line exactly, so detrending would zero them and every cell would score identically
        return values - values.mean()
    slope, intercept = np.polyfit(x, values, 1)
    return values - (slope * x + intercept)


# raw and detrended RMSE for every cell over one window
def score(frame, observed, lo, hi):
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


# Run: rescore every window
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
