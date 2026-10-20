#!/usr/bin/env python3
"""
Mirror check: does the runner's groin handle drift from either direction the same way?

    python mirror_check.py   ->  tables/mirror_one_year.csv, tables/mirror_runs.csv, tables/brie_table_symmetry.csv

Three parts, on the straight-coast rig of ../6-stress/ (61 domains, the runner's GroinCallback):

    A  one call of the callback on a shoreline and on its exact reversal (domain order flipped, so theta0 becomes
       -theta0, the groin face LO|LO+1 becomes NY-2-LO, and the waves come from the other side, with BRIE's tables
       reflected exactly). Its x_s_dt must come back exactly reversed: the groin's side, sign and approach angle.
    B  whole 100-year runs and their mirrors, with the groin and without it. BRIE's own step is not mirror-
       symmetric (each row takes r from its forward face), so no run mirrors exactly; the question is whether
       the groin adds to BRIE's asymmetry.
    C  how far BRIE's tables for asymmetry 1 - a are from an exact reflection of those for a (its 180 angle bins
       are 1.006 deg wide but are looked up as 1 deg).

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-10-15
"""
from __future__ import annotations

import sys
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "6-stress"))
import stress_test as st  # noqa: E402
from cascade.groin import GroinCallback  # noqa: E402

# --- CONFIG ------------------------------------------------------------------
NY, LO = st.NY, st.LO
MIRROR_LO = NY - 2 - LO                       # face LO|LO+1 reflected
STRENGTHS = (0.25, 1.0)
NOTCHES_M = (0.0, -150.0, 150.0)
TOLERANCE_M = 1e-9                            # what "mirrors exactly" means for one call
TABLES = HERE / "tables"
# -----------------------------------------------------------------------------

_drift_table = GroinCallback._drift_table


# The drift table the callback reads: a stand-in's own reflected table if it carries one, else BRIE's
def drift_table(brie):
    return brie._q_table if hasattr(brie, "_q_table") else _drift_table(brie)


# A stand-in carrying BRIE's tables for `climate`, reflected exactly: index k (theta = 90 - k) <-> 180 - k
def reflected(climate):
    b = st.brie_like(climate)
    q = _drift_table(b)
    k = np.arange(1, b._wave_climl)
    d2, q2 = b._coast_diff.copy(), q.copy()
    d2[k], q2[k] = b._coast_diff[180 - k], -q[180 - k]
    return SimpleNamespace(**{**vars(b), "_coast_diff": d2, "_q_table": q2, "x_s": None})


# Register a stand-in under a climate key the rig can name
def register(name, b):
    climate = dict(st.OPTION_A, asym=name)
    st._BRIE[tuple(sorted(climate.items()))] = b
    return climate


# Part A: one callback call on the shoreline and on its reversal; how far the reversed x_s_dt is from the first
def one_year(theta0, strength, notch_m, b, b_mirror):
    slope = np.tan(np.deg2rad(theta0))
    origin = (LO + 1 + NY // 2) % NY
    side = np.where(np.arange(NY) <= LO, 1.0, -1.0) * (1.0 if st.drift0(b, theta0) > 0 else -1.0)
    dist = np.minimum(np.abs(np.arange(NY) - (LO + 0.5)), NY - np.abs(np.arange(NY) - (LO + 0.5))) - 0.5
    x = slope * st.DY_M * ((np.arange(NY) - origin) % NY) - 0.5 * notch_m * side * np.exp(-dist / st.NOTCH_DECAY_CELLS)
    out = []
    for brie, shore, lo, th in ((b, x, LO, theta0), (b_mirror, x[::-1].copy(), MIRROR_LO, -theta0)):
        up, down = (lo, lo + 1) if st.drift0(brie, th) > 0 else (lo + 1, lo)
        cb = GroinCallback(up, down, [(1, strength)], start_year=1, n_domains=NY)
        brie.x_s = shore
        out.append(np.asarray(cb(SimpleNamespace(_brie_coupler=SimpleNamespace(_brie=brie)), list(np.zeros(NY))),
                              dtype=float))
    return dict(theta0=theta0, strength=strength, notch_m=notch_m, held_m=float(np.abs(out[0]).max()),
                max_diff_m=float(np.abs(out[0] - out[1][::-1]).max()))


# Part B: a whole run and its mirror; how far the mirror's final shoreline is from the reflected run
def pair(theta0, strength, notch_m, mirror_climate, groin=True):
    a = st.run(theta0, strength=strength, notch_m=notch_m, groin=groin)
    m = st.run(-theta0, strength=strength, notch_m=notch_m, groin=groin, climate=mirror_climate, lo=MIRROR_LO)
    return float(np.abs(a.attrs["eta"] - m.attrs["eta"][::-1]).max()), float(np.abs(a.attrs["eta"]).max())


# Run: parts C, A and B
def main():
    GroinCallback._drift_table = staticmethod(drift_table)   # in-process only: reads the stand-in's table
    TABLES.mkdir(parents=True, exist_ok=True)
    pd.set_option("display.width", 200)

    # Part C: BRIE's own 1 - a tables against an exact reflection of its a tables
    sym = []
    for asym in (0.55, 0.6, 0.65, 0.7):
        b, r = st.brie_like(dict(st.OPTION_A, asym=asym)), reflected(dict(st.OPTION_A, asym=asym))
        own = st.brie_like(dict(st.OPTION_A, asym=round(1 - asym, 2)))
        k = np.arange(1, b._wave_climl)
        sym.append(dict(asymmetry=asym,
                        diffusivity_max_rel_diff=float(np.abs(own._coast_diff[k] - r._coast_diff[k]).max()
                                                       / np.abs(b._coast_diff).max()),
                        drift_max_rel_diff=float(np.abs(_drift_table(own)[k] - r._q_table[k]).max()
                                                 / np.abs(r._q_table).max())))
    sym = pd.DataFrame(sym)
    sym.to_csv(TABLES / "brie_table_symmetry.csv", index=False, float_format="%.5f")
    print("C. BRIE's 1 - a tables against an exact reflection:")
    print(sym.round(4).to_string(index=False))

    exact_b = reflected(st.OPTION_A)
    exact = register("reflected-0.6", exact_b)
    b = st.brie_like(st.OPTION_A)

    # Part A
    one = pd.DataFrame([one_year(th, s, n, b, exact_b) for th in st.ORIENTATIONS for s in STRENGTHS
                        for n in NOTCHES_M])
    one["mirrors"] = one.max_diff_m <= TOLERANCE_M
    one.to_csv(TABLES / "mirror_one_year.csv", index=False, float_format="%.3e")
    print(f"\nA. one call: {int(one.mirrors.sum())} of {len(one)} mirror to {TOLERANCE_M:g} m "
          f"(largest difference {one.max_diff_m.max():.2e} m; the groin moves up to {one.held_m.max():.1f} m)")

    # Part B
    rows = []
    for th in st.ORIENTATIONS:
        for notch in NOTCHES_M:
            no_groin, no_groin_size = pair(th, 0.0, notch, exact, groin=False)
            for s in STRENGTHS:
                with_groin, size = pair(th, s, notch, exact)
                rows.append(dict(theta0=th, notch_m=notch, strength=s, mirror_diff_groin_m=with_groin,
                                 largest_change_groin_m=size, mirror_diff_no_groin_m=no_groin,
                                 largest_change_no_groin_m=no_groin_size))
    runs = pd.DataFrame(rows)
    runs["groin_adds_m"] = runs.mirror_diff_groin_m - runs.mirror_diff_no_groin_m
    runs.to_csv(TABLES / "mirror_runs.csv", index=False, float_format="%.3f")
    print("\nB. whole runs, final shoreline (m):")
    print(runs.round(2).to_string(index=False))


if __name__ == "__main__":
    main()
