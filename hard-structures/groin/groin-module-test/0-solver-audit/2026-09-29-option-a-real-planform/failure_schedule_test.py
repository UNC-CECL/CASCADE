#!/usr/bin/env python3
"""
Does an instant groin failure at the 2003 storm let one groin fit both windows?

    python failure_schedule_test.py

Blocking (blocking_groin_emulator.py, callback) and dipole groins at full
strength until the failure year and x f after it, scored as the adjusted OLS
D5-D6 gap change against observed_fillet_m. Writes failure_schedule_scores.csv.
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
sys.path[:0] = [str(HERE), str(HERE.parent)]
import blocking_groin_emulator as bg  # noqa: E402

# --- CONFIG ------------------------------------------------------------------
FAIL_YEARS = (2003, 2004)
B_VALUES = np.round(np.arange(0.0, 1.01, 0.05), 2)
M_VALUES = (0, 1, 2, 3, 4, 5, 6, 7, 8, 10, 12, 15)
F_VALUES = (0.0, 0.05, 0.1, 0.2, 0.3)
# -----------------------------------------------------------------------------


# Score every failure year x window x f x strength, blocking and dipole
def run_all():
    rows = []
    planforms = {s: bg.load_planform(s) for s in (1996, 2010)}
    offsets = {s: bg.FULL_NOGROIN[s] - bg.gap_change(bg.solve(planforms[s], s, 0.0, 1.0, "exact"))
               for s in planforms}
    for fail in FAIL_YEARS:
        # Instant schedule: swap the module's ramp for a step at `fail`.
        bg.b_of_year = lambda b0, f, year, fail=fail: b0 if year < fail else b0 * f
        for start, x0 in planforms.items():
            for f in F_VALUES:
                for b0 in B_VALUES:
                    ch = bg.gap_change(bg.solve(x0, start, b0, f, "callback"))
                    rows.append(dict(kind="blocking", fail=fail, start=start,
                                     strength=b0, f=f, adjusted_m=ch + offsets[start]))
                for M in M_VALUES:
                    ch = bg.gap_change(dipole(x0, start, M, f, fail))
                    rows.append(dict(kind="dipole", fail=fail, start=start,
                                     strength=M, f=f, adjusted_m=ch + offsets[start]))
    t = pd.DataFrame(rows)
    t["observed_m"] = t.start.map(bg.OBSERVED)
    t["err_m"] = t.adjusted_m - t.observed_m
    return t


# GroinCallback's +/-M dipole on the same solve (x_s_dt only), instant failure at `fail`
def dipole(x0, start, M, f, fail):
    ny = x0.size
    cd, _, _ = bg.brie_diffusivity(2.0, 7.5, 0.6, 0.5, ny)
    x = x0.copy()
    states = [x.copy()]
    for t in range(1, bg.YEARS + 1):
        m = M if (start + t - 1) < fail else M * f
        theta = 180.0 * np.arctan2(x[np.r_[1:ny, 0]] - x, bg.DY_M) / np.pi
        r = np.maximum(0.0, cd[np.clip(np.round(90 - theta).astype(int), 1, 179)]
                       * bg.DT_YR / 2.0 / bg.DY_M ** 2)
        x_s_dt = np.zeros(ny)
        x_s_dt[bg.U] -= m
        x_s_dt[bg.D] += m
        rows = np.r_[np.arange(ny), np.arange(ny), np.arange(ny)]
        cols = np.r_[np.arange(ny), np.r_[1:ny, 0], np.r_[ny - 1, 0:ny - 1]]
        vals = np.r_[1.0 + 2 * r, -r, -r]
        rhs = (x + r * (x[np.r_[1:ny, 0]] - x) + r * (x[np.r_[ny - 1, 0:ny - 1]] - x)
               + x_s_dt)
        x = bg.spsolve(bg.csr_matrix((vals, (rows, cols)), shape=(ny, ny)), rhs)
        states.append(x.copy())
    return np.array(states)


# Run: score, write the CSV, print the best joint cells
def main():
    t = run_all()
    t.to_csv(HERE / "failure_schedule_scores.csv", index=False)
    pd.set_option("display.width", 200)
    for (kind, fail), sub in t.groupby(["kind", "fail"]):
        w = sub.pivot_table(index=["strength", "f"], columns="start", values="err_m")
        w["rms"] = np.sqrt((w[1996] ** 2 + w[2010] ** 2) / 2)
        print(f"\n=== {kind}, instant failure {fail}: best joint (err m; obs 1996 -4.3, 2010 -60.4) ===")
        print(w.sort_values("rms").head(6).round(1).to_string())


if __name__ == "__main__":
    main()
