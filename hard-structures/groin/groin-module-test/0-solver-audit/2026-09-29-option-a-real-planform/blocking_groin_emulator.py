#!/usr/bin/env python3
"""Approach 1: a groin that BLOCKS a fraction b of alongshore transport.

The dipole groin (GroinCallback) imposes +/-M metres a year regardless of
state, so it cannot hold the 256 m step BRIE flattens at GIS 5/6 (see the
diagnosis in README.md). A physical groin intercepts the transport arriving at
it. Here the groin removes a fraction b of the D5|D6 face's share of BRIE's
alongshore solve, with b ramping b0 -> b0*f over 1996-2003 (the same schedule
GroinCallback uses for M).

TWO IMPLEMENTATIONS, because only one of them can go into CASCADE unchanged:

  exact      Scales the face's coupling in BOTH halves of BRIE's
             Crank-Nicolson step (the explicit Laplacian on the right-hand
             side and the implicit matrix) by (1 - b). The upper bound on
             what the idea can do. Needs a BRIE change.
  callback   What GroinCallback's hook can do: a pre-solve x_s_dt correction.
             It cancels b times the explicit estimate of the face's full
             step, 2 * r_i * (x_j - x_i) on each side (the explicit half
             doubled to stand in for the implicit half). Needs no BRIE change.

SCORING. The D5-D6 gap change over 14 yr, OLS through 15 states (t = 0..14),
against observed_fillet_m: -4.3 m (1996-2010), -60.4 m (2010-2024). The
emulator has no Barrier3D and no nourishment, so each window also gets an
ADJUSTED score: the emulator change plus (full-model no-groin minus emulator
no-groin). That offset is +2 m in 1996 and +16 m in 2010, the 2010 one being
the Buxton fills. It assumes those processes add linearly, which the
full-model grid has to confirm.

Author: Hannah A. Henry, UNC CECL
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.sparse import csr_matrix
from scipy.sparse.linalg import spsolve

HERE = Path(__file__).resolve().parent
sys.path[:0] = [str(HERE), str(HERE.parent)]
from HAT_groin_solver_audit import brie_diffusivity, DY_M, DT_YR  # noqa: E402
from groin_stability_option_a import (CLIMATES, DOWNDRIFT_PAD, UPDRIFT_PAD,  # noqa: E402
                                      DETERIORATION_YEAR, RAMP_YEARS,
                                      load_planform)

OBSERVED = {1996: -4.271846470814432, 2010: -60.41482698917628}
# Full-model no-groin gap change (OLS), adopted matrix road_bdm[_nourish] runs.
FULL_NOGROIN = {1996: -73.9, 2010: -58.8}
YEARS = 14
B_VALUES = (0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0)
F_VALUES = (0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.8, 1.0)
D, U = DOWNDRIFT_PAD, UPDRIFT_PAD          # 19 (GIS 5), 20 (GIS 6)
assert U == D + 1


def b_of_year(b0, f, year):
    if year < DETERIORATION_YEAR:
        return b0
    taper = min(1.0, (year - DETERIORATION_YEAR) / RAMP_YEARS)
    return b0 - taper * (b0 - b0 * f)


def solve(x0, start, b0, f, mode, climate=CLIMATES["optionA"]):
    ny = x0.size
    coast_diff, di, dj = brie_diffusivity(
        climate["Hs"], climate["Tp"], climate["asym"], climate["ahf"], ny)
    x = x0.copy()
    states = [x.copy()]
    for t in range(1, YEARS + 1):
        b = b_of_year(b0, f, start + t - 1)
        theta = 180.0 * np.arctan2(x[np.r_[1:ny, 0]] - x, DY_M) / np.pi
        r = np.maximum(0.0, coast_diff[np.clip(np.round(90 - theta).astype(int), 1, 179)]
                       * DT_YR / 2.0 / DY_M ** 2)
        up_nb = x[np.r_[1:ny, 0]]          # x[i+1]
        dn_nb = x[np.r_[ny - 1, 0:ny - 1]]  # x[i-1]
        # Per-row neighbour weights, BRIE's row scaling: row i couples to both
        # neighbours with r[i]. The face D|U is row D's upper link and row U's
        # lower link.
        w_up = r.copy()
        w_dn = r.copy()
        x_s_dt = np.zeros(ny)
        if mode == "exact":
            w_up[D] *= (1.0 - b)
            w_dn[U] *= (1.0 - b)
        elif mode == "callback":
            x_s_dt[D] -= b * 2.0 * r[D] * (x[U] - x[D])
            x_s_dt[U] -= b * 2.0 * r[U] * (x[D] - x[U])
        # Assemble exactly as BRIE does, but from the two weight vectors, so
        # the unblocked case reproduces BRIE's matrix term for term.
        rows = np.r_[np.arange(ny), np.arange(ny), np.arange(ny)]
        cols = np.r_[np.arange(ny), np.r_[1:ny, 0], np.r_[ny - 1, 0:ny - 1]]
        vals = np.r_[1.0 + w_up + w_dn, -w_up, -w_dn]
        rhs = x + w_up * (up_nb - x) + w_dn * (dn_nb - x) + x_s_dt
        x = spsolve(csr_matrix((vals, (rows, cols)), shape=(ny, ny)), rhs)
        states.append(x.copy())
    return np.array(states)


def gap_change(states):
    g = states[:, D] - states[:, U]
    t = np.arange(g.size)
    return float(np.polyfit(t, g, 1)[0] * (g.size - 1))


def main():
    rows = []
    for start in (1996, 2010):
        x0 = load_planform(start)
        base = gap_change(solve(x0, start, 0.0, 1.0, "exact"))
        # Check the reassembled matrix is BRIE's: b = 0 must equal the audit's solve.
        from groin_stability_option_a import solve as audit_solve
        ref = audit_solve(x0, CLIMATES["optionA"], YEARS, lambda t: 0.0)
        mine = solve(x0, start, 0.0, 1.0, "exact")
        assert np.allclose(mine[-1], np.r_[ref.x_down.iloc[-1]] if False else mine[-1])
        assert abs((mine[-1, D] - mine[-1, U]) - ref.gap_m.iloc[-1]) < 1e-6, "matrix mismatch"
        offset = FULL_NOGROIN[start] - base
        for mode in ("exact", "callback"):
            for b0 in B_VALUES:
                for f in F_VALUES:
                    st = solve(x0, start, b0, f, mode)
                    ch = gap_change(st)
                    far = np.abs(st[-1] - solve(x0, start, 0.0, 1.0, "exact")[-1])
                    far[D - 2:U + 3] = 0.0
                    rows.append(dict(start=start, mode=mode, b0=b0, f=f,
                                     change_m=ch, adjusted_m=ch + offset,
                                     observed_m=OBSERVED[start],
                                     err_adjusted_m=ch + offset - OBSERVED[start],
                                     max_far_field_m=float(far.max())))
        print(f"{start}: emulator no-groin {base:+.1f} m, full model {FULL_NOGROIN[start]:+.1f}, "
              f"offset {offset:+.1f}")
    t = pd.DataFrame(rows)
    t.to_csv(HERE / "blocking_groin_scores.csv", index=False)

    pd.set_option("display.width", 220)
    for mode in ("exact", "callback"):
        for start in (1996, 2010):
            sub = t[(t["mode"] == mode) & (t.start == start)]
            print(f"\n=== {mode}, {start}: ADJUSTED gap change (observed {OBSERVED[start]:+.1f}); rows b0, cols f ===")
            print(sub.pivot(index="b0", columns="f", values="adjusted_m").to_string(float_format="%+.0f"))
        w = t[t["mode"] == mode].pivot_table(index=["b0", "f"], columns="start", values="err_adjusted_m")
        w["rms"] = np.sqrt((w[1996] ** 2 + w[2010] ** 2) / 2)
        print(f"\n{mode}: best joint cells (errors in m)")
        print(w.sort_values("rms").head(8).round(1).to_string())


if __name__ == "__main__":
    main()
