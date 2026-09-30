#!/usr/bin/env python3
"""Does an INSTANT failure at the 2003 storm let one groin fit both windows?

The observed D5-D6 gap (wet/dry table, 24 dates) grows 1967-1995, holds at
134-155 m through 2004, then falls to 125 (2008), 104 (2016), 63-74 (2019-23).
That is an intact groin failing after the 2003 storm, NOT the linear
1996 -> 2003 wear-down both groin emulators were given -- which puts the decline
inside the 1996 window, where the data show none.

Tests both representations under GroinCallback's existing "instant" mode (full
strength until the failure year, x f from then on; no new code needed for the
schedule):
  blocking  approach 1, callback implementation (blocking_groin_emulator.py)
  dipole    today's GroinCallback, M m/yr

Scored as there: adjusted OLS gap change vs observed_fillet_m.

Author: Hannah A. Henry, UNC CECL
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path[:0] = [str(HERE), str(HERE.parent)]
import blocking_groin_emulator as bg  # noqa: E402

FAIL_YEARS = (2003, 2004)
B_VALUES = np.round(np.arange(0.0, 1.01, 0.05), 2)
M_VALUES = (0, 1, 2, 3, 4, 5, 6, 7, 8, 10, 12, 15)
F_VALUES = (0.0, 0.05, 0.1, 0.2, 0.3)


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


def dipole(x0, start, M, f, fail):
    """GroinCallback's +/-M dipole on the same solve (x_s_dt only)."""
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
