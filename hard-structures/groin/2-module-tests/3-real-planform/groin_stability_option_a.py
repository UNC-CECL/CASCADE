#!/usr/bin/env python3
"""
Is the groin dipole stable under the option A wave climate, on the real Hatteras coast?

    python groin_stability_option_a.py

The solver audit's emulator (BRIE's alongshore solve plus the dipole; no
Barrier3D, source/sink or storms), started from the adopted matrix runs' BRIE
shoreline after year 1, under the old /10 climate and option A. Writes
stability_summary.csv and stability_traces.csv. Needs brie, numpy, pandas, scipy.
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
from scipy.sparse import csr_matrix
from scipy.sparse.linalg import spsolve

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "2-solver-audit"))
from HAT_groin_solver_audit import brie_diffusivity, DY_M, DT_YR  # noqa: E402

# --- CONFIG ------------------------------------------------------------------
REPO = next(p for p in Path(__file__).resolve().parents if (p / "pyproject.toml").exists())
MATRIX = REPO / "output" / "raw_runs" / "matrix"
PLANFORM_RUNS = {                  # the adopted road_nobdm edgeBE runs (README)
    1996: MATRIX / "1996_2010/edgeBE/HAT_1996_2010_edgeBE_offsetmetres_road_nobdm_nogroin",
    2010: MATRIX / "2010_2024/edgeBE/HAT_2010_2024_edgeBE_offsetmetres_road_nobdm_nogroin",
}
RUN_YEARS = 14                     # 1996-2009 and 2010-2023, as the matrix runs them

NUM_BUFFER = 15                    # GroinCallback wiring as in HAT_groin_sweep_config
UPDRIFT_PAD, DOWNDRIFT_PAD = NUM_BUFFER + 6 - 1, NUM_BUFFER + 5 - 1   # GIS 6 updrift, GIS 5 downdrift
REAL = slice(NUM_BUFFER, NUM_BUFFER + 90)
INSTALL_YEAR, DETERIORATION_YEAR, RAMP_YEARS = 1969, 1996, 7

CLIMATES = {
    "optionA": dict(Hs=2.0, Tp=7.5, asym=0.6, ahf=0.5),
    "old_div10": dict(Hs=2.5, Tp=8.0, asym=0.7, ahf=0.1),
}
M_VALUES = (1, 2, 3, 5, 10, 20, 40, 60, 80, 120)
F_VALUES = (0.6, 1.0)
LONG_YEARS = 200                   # constant-M run for the stability boundary
FILLET_TARGET_M = 22.11            # the period-1 fillet M = 60 was fitted to
# -----------------------------------------------------------------------------


# BRIE x_s (metres, landward +) after the first model year of a matrix run
def load_planform(start):
    run = PLANFORM_RUNS[start]
    npz = np.load(run / f"{run.name}.npz", allow_pickle=True)
    brie = npz["cascade"][0].brie
    return np.asarray(brie._x_s_save, float)[:, 1].copy()


# Emulated BRIE alongshore solve with one dipole; one row per year
def solve(x0, climate, years, M_of_year):
    ny = x0.size
    coast_diff, di, dj = brie_diffusivity(
        climate["Hs"], climate["Tp"], climate["asym"], climate["ahf"], ny)
    x_s = x0.copy()
    rows = []
    for t in range(1, years + 1):
        theta = 180.0 * np.arctan2(x_s[np.r_[1:ny, 0]] - x_s, DY_M) / np.pi
        r_ipl = np.maximum(
            0.0,
            coast_diff[np.clip(np.round(90 - theta).astype(int), 1, 179)]
            * DT_YR / 2.0 / DY_M ** 2)
        M = M_of_year(t)
        x_s_dt = np.zeros(ny)
        x_s_dt[UPDRIFT_PAD] -= M
        x_s_dt[DOWNDRIFT_PAD] += M
        values = np.r_[-r_ipl[-1], -r_ipl[1:], 1 + 2 * r_ipl, -r_ipl[0:-1], -r_ipl[0]]
        rhs = (x_s + r_ipl * (x_s[np.r_[1:ny, 0]] - 2 * x_s
                              + x_s[np.r_[ny - 1, 0:ny - 1]]) + x_s_dt)
        x_s = spsolve(csr_matrix((values, (di, dj))), rhs)
        rows.append(dict(
            t=t, M_applied=M,
            r_ipl_up=r_ipl[UPDRIFT_PAD], r_ipl_down=r_ipl[DOWNDRIFT_PAD],
            theta_up=theta[UPDRIFT_PAD], theta_down=theta[DOWNDRIFT_PAD],
            n_real_shut=int((r_ipl[REAL] <= 0).sum()),
            n_ring_shut=int((r_ipl <= 0).sum()),
            gap_m=x_s[DOWNDRIFT_PAD] - x_s[UPDRIFT_PAD],
            x_up=x_s[UPDRIFT_PAD], x_down=x_s[DOWNDRIFT_PAD],
            mean_real=x_s[REAL].mean()))
    return pd.DataFrame(rows)


# GroinCallback._effective_trapping_rate, linear_ramp mode
def effective_M(M, f, year):
    if year < DETERIORATION_YEAR:
        return M
    taper = min(1.0, (year - DETERIORATION_YEAR) / RAMP_YEARS)
    return M - taper * (M - M * f)


# Groin run minus the no-groin run from the same planform
def run_case(x0, climate, years, M_of_year, baseline):
    g = solve(x0, climate, years, M_of_year)
    g["fillet_m"] = g.gap_m - baseline.gap_m.values
    g["updrift_advance_m"] = baseline.x_up.values - g.x_up
    g["downdrift_retreat_m"] = g.x_down - baseline.x_down.values
    g["new_real_shut"] = g.n_real_shut - baseline.n_real_shut.values
    return g


# One row per case: end fillet, whether anything shut down, and when
def summarise(g):
    shut = g.loc[g.new_real_shut > 0, "t"]
    return dict(
        fillet_end_m=float(g.fillet_m.iloc[-1]),
        updrift_advance_end_m=float(g.updrift_advance_m.iloc[-1]),
        downdrift_retreat_end_m=float(g.downdrift_retreat_m.iloc[-1]),
        first_shutdown_yr=int(shut.iloc[0]) if len(shut) else None,
        max_new_real_shut=int(g.new_real_shut.max()),
        min_r_ipl_groin=float(np.minimum(g.r_ipl_up, g.r_ipl_down).min()),
        theta_down_end=float(g.theta_down.iloc[-1]),
        theta_up_end=float(g.theta_up.iloc[-1]),
        fillet_still_growing_m_yr=float(g.fillet_m.diff().iloc[-1]),
    )


# Run: every start x climate x window x M x f, then the tables
def main():
    rows, traces = [], []
    for start in PLANFORM_RUNS:
        x_metres = load_planform(start)
        for cname, climate in CLIMATES.items():
            # The old climate is judged on the /10 planform it was calibrated on
            x0 = x_metres / 10.0 if cname == "old_div10" else x_metres
            coast_diff, _, _ = brie_diffusivity(
                climate["Hs"], climate["Tp"], climate["asym"], climate["ahf"], x0.size)
            theta0 = 180 * np.arctan2(np.roll(x0, -1) - x0, DY_M) / np.pi
            D0 = coast_diff[np.clip(np.round(90 - theta0).astype(int), 1, 179)]
            print(f"\n{start} {cname}: groin cell angles "
                  f"up {theta0[UPDRIFT_PAD]:+.1f} / down {theta0[DOWNDRIFT_PAD]:+.1f} deg, "
                  f"D {D0[UPDRIFT_PAD]:.0f} / {D0[DOWNDRIFT_PAD]:.0f} m2/yr; "
                  f"real-reach median D {np.median(D0[REAL]):.0f}")

            for window, years, schedule in (
                    ("hindcast", RUN_YEARS, "ramp"),
                    ("long_constant", LONG_YEARS, "constant")):
                base = solve(x0, climate, years, lambda t: 0.0)
                for M in M_VALUES:
                    for f in F_VALUES:
                        if schedule == "constant":
                            if f != 1.0:
                                continue
                            M_of_year = lambda t, M=M: M
                        else:
                            M_of_year = (lambda t, M=M, f=f, s=start:
                                         effective_M(M, f, s + t - 1))
                        g = run_case(x0, climate, years, M_of_year, base)
                        rows.append(dict(start=start, climate=cname, window=window,
                                         years=years, M=M, f=f if schedule == "ramp" else None,
                                         **summarise(g)))
                        g.insert(0, "M", M)
                        g.insert(0, "f", f)
                        g.insert(0, "window", window)
                        g.insert(0, "climate", cname)
                        g.insert(0, "start", start)
                        traces.append(g)

    # Write both tables and print each window
    table = pd.DataFrame(rows)
    table.to_csv(HERE / "stability_summary.csv", index=False)
    pd.concat(traces).to_csv(HERE / "stability_traces.csv", index=False)

    pd.set_option("display.width", 200)
    cols = ["start", "climate", "M", "f", "fillet_end_m", "first_shutdown_yr",
            "max_new_real_shut", "theta_down_end", "fillet_still_growing_m_yr"]
    for window in ("hindcast", "long_constant"):
        print(f"\n=== {window} ===")
        print(table.loc[table.window == window, cols].to_string(index=False,
                                                                 float_format="%.2f"))

    # M giving the period-1 fillet, by interpolation over the stable f = 0.6 cells
    print(f"\nM giving a {FILLET_TARGET_M} m fillet over 14 yr at f = 0.6 "
          "(linear interpolation, stable cells only):")
    for (start, cname), sub in table[(table.window == "hindcast") & (table.f == 0.6)
                                     & table.first_shutdown_yr.isna()].groupby(["start", "climate"]):
        sub = sub.sort_values("fillet_end_m")
        m = np.interp(FILLET_TARGET_M, sub.fillet_end_m, sub.M, left=np.nan, right=np.nan)
        print(f"  {start} {cname:10s} M ~ {m:.1f}")


if __name__ == "__main__":
    main()
