#!/usr/bin/env python3
"""
Is the drift groin's sand leak BRIE's solve? The drift module rerun with a conservative solve.

    python straight_coast_conservative_solve_check.py [--years 50]   ->  comparison/tables/conservative_solve_check.csv

The same straight-coast emulator, with BRIE's row-scaled diffusion (each cell's row uses that cell's
own r for both neighbours) swapped for a face-centred one (each face's r used by both cells it joins),
so sand that leaves one cell arrives in the next. The tilt's flux is kept explicitly, since it no longer
cancels where r varies along the coast. A diagnostic only: BRIE itself is unchanged.

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-10-08
"""
from __future__ import annotations

import argparse

import numpy as np
import pandas as pd
from scipy.sparse import diags
from scipy.sparse.linalg import spsolve

import straight_coast_test as sc

OUT = sc.HERE / "comparison" / "tables" / "conservative_solve_check.csv"


# One run on a face-centred Crank-Nicolson solve: module "none" or "drift"
def run_face_centred(t, theta0, module, years, blocking_fraction=sc.BLOCKING_FRACTION):
    ny, lo, hi, dy = sc.NY, sc.LO, sc.HI, sc.DY_M
    slope = np.tan(np.deg2rad(theta0))
    eta = np.zeros(ny)
    up = hi if float(sc.drift(t, theta0)[0]) < 0 else lo
    for _ in range(years):
        step = slope * dy + (np.roll(eta, -1) - eta)
        theta = 180.0 * np.arctan2(step, dy) / np.pi
        rf = np.maximum(0.0, t["coast_diff"][sc.angle_index(theta, t["climl"])] * sc.DT_YR / 2 / dy ** 2)
        rm = np.roll(rf, 1)                                   # the face on each cell's lower side
        dx = np.zeros(ny)
        if module == "drift":
            approach = theta[hi] if up == hi else theta[lo - 1]
            q = blocking_fraction * float(sc.drift(t, approach)[0]) * sc.DT_YR / (dy * t["depth"])
            dx[lo] -= q
            dx[hi] += q
        lap = rf * (np.roll(eta, -1) - eta) - rm * (eta - np.roll(eta, 1))
        tilt = 2 * slope * dy * (rf - rm)
        a = diags([1 + rf + rm, -rf[:-1], -rm[1:], [-rf[-1]], [-rm[0]]],
                  [0, 1, -1, -(ny - 1), ny - 1], format="csr")
        eta = spsolve(a, eta + lap + tilt + dx)
    return float(-eta[up]), float(-eta.sum())


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--years", type=int, default=50)
    a = ap.parse_args()
    t = sc.brie_tables()
    summ = pd.read_csv(sc.HERE / "comparison" / "tables" / "straight_coast_summary_year50.csv")
    brie = summ[summ.module == "drift"].set_index("theta0_deg")
    rows = []
    for th in sc.ORIENTATIONS:
        _, none_total = run_face_centred(t, th, "none", a.years)
        up, total = run_face_centred(t, th, "drift", a.years)
        rows.append(dict(theta0_deg=th, no_groin_total_m=none_total, drift_updrift_m=up, drift_total_m=total,
                         brie_drift_updrift_m=brie.loc[th, "updrift_seaward_m"],
                         brie_drift_total_m=brie.loc[th, "reach_change_m"]))
    frame = pd.DataFrame(rows)
    frame.to_csv(OUT, index=False, float_format="%.2f")
    print(frame.round(1).to_string(index=False))


if __name__ == "__main__":
    main()
