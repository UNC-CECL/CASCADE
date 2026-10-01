#!/usr/bin/env python3
"""
Stage 0 of the groin module test: BRIE's alongshore solve under a groin dipole, without CASCADE.

    python HAT_groin_solver_audit.py [--years 200] [--outdir DIR]

An emulator of BRIE's implicit shoreline-diffusion solve plus GroinCallback's
+/-M dipole, built from BRIE's own coast_diff table and sparse indices. Six
audits: diffusivity, fillet vs M, closure vs f, domain count, groin fields, the
chosen rig. Writes solver_audit_<audit>.csv. Needs brie, numpy, pandas, scipy.
Details: README.md beside this script.

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-10-01
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.sparse import csr_matrix
from scipy.sparse.linalg import spsolve

from brie import Brie

# --- CONFIG ------------------------------------------------------------------
HERE = Path(__file__).resolve().parent
DY_M = 500.0                 # fixed by the CASCADE coupler ("do not change")
DT_YR = 1.0                  # fixed by the CASCADE coupler ("do not change")
DEFAULT_CLIMATE = dict(Hs=1.0, Tp=7.0, asym=0.8, ahf=0.2)    # CASCADE / Barrier3D default
HATTERAS_CLIMATE = dict(Hs=2.5, Tp=9.0, asym=0.8, ahf=0.2)   # not interchangeable with the default (README)
# -----------------------------------------------------------------------------

# Diffusivity tables already built, by (Hs, Tp, asym, ahf, ny)
_CD_CACHE: dict = {}


# (coast_diff, di, dj) from a real Brie instance, so the audit cannot drift from BRIE
def brie_diffusivity(Hs, Tp, asym, ahf, ny):
    key = (Hs, Tp, asym, ahf, ny)
    if key not in _CD_CACHE:
        brie = Brie(
            ast_model=True,
            barrier_model=False,
            inlet_model=False,
            b3d=True,
            wave_height=Hs,
            wave_period=Tp,
            wave_asymmetry=asym,
            wave_angle_high_fraction=ahf,
            alongshore_section_length=DY_M,
            alongshore_section_count=ny,
            time_step=DT_YR,
            time_step_count=10,
            save_spacing=1,
        )
        _CD_CACHE[key] = (brie._coast_diff.copy(), brie._di.copy(), brie._dj.copy())
    return _CD_CACHE[key]


# Shoreline angles (negative, positive cutoff, deg) bounding BRIE's positive-diffusivity band
def shutdown_angle_deg(coast_diff):
    theta = np.arange(-89, 90)
    idx = np.clip(np.round(90 - theta).astype(int), 1, 179)
    positive = theta[coast_diff[idx] > 0]
    return float(positive.min()), float(positive.max())


# BRIE's shoreline diffusion with groin dipoles only, from a straight coast: (frame, shutdown_year, profiles)
def solve_reach(ny, groins, years, climate, record_years=()):
    coast_diff, di, dj = brie_diffusivity(
        climate["Hs"], climate["Tp"], climate["asym"], climate["ahf"], ny
    )
    for updrift, downdrift, _, _ in groins:
        if abs(updrift - downdrift) != 1:
            raise ValueError(
                "groin domains must be adjacent, as GroinCallback requires; "
                f"got updrift={updrift}, downdrift={downdrift}."
            )

    x_s = np.zeros(ny)
    # Net source the field injects per year, seaward-negative (zero when every f = 1)
    net_source = -sum(M * (1.0 - f) for _, _, M, f in groins)

    shutdown_year = None
    rows, profiles = [], {}

    for year in range(1, years + 1):
        # BRIE's forward-difference angle and the diffusion number it selects, clipped at zero as BRIE does
        theta = 180.0 * np.arctan2(x_s[np.r_[1:ny, 0]] - x_s, DY_M) / np.pi
        r_ipl = np.maximum(
            0.0,
            coast_diff[np.clip(np.round(90 - theta).astype(int), 1, 179)]
            * DT_YR / 2.0 / DY_M**2,
        )
        if shutdown_year is None and np.any(r_ipl <= 0):
            shutdown_year = year

        x_s_dt = np.zeros(ny)
        for updrift, downdrift, M, f in groins:
            x_s_dt[updrift] -= M
            x_s_dt[downdrift] += M * f

        values = np.r_[-r_ipl[-1], -r_ipl[1:], 1 + 2 * r_ipl, -r_ipl[0:-1], -r_ipl[0]]
        rhs = (
            x_s
            + r_ipl * (x_s[np.r_[1:ny, 0]] - 2 * x_s + x_s[np.r_[ny - 1, 0:ny - 1]])
            + x_s_dt
        )
        x_s = spsolve(csr_matrix((values, (di, dj))), rhs)

        mean_if_conservative = year * net_source / ny
        row = dict(
            year=year,
            mean_m=float(x_s.mean()),
            mean_if_conservative_m=float(mean_if_conservative),
            closure_error_m=float(x_s.mean() - mean_if_conservative),
            min_r_ipl=float(r_ipl.min()),
            max_abs_offset_m=float(np.abs(np.diff(np.r_[x_s, x_s[0]])).max()),
        )
        for i, (updrift, downdrift, _, _) in enumerate(groins):
            row[f"fillet_{i}_m"] = float(x_s[downdrift] - x_s[updrift])
        row["fillet_m"] = row["fillet_0_m"]
        rows.append(row)
        if year in record_years:
            profiles[year] = x_s.copy()

    return pd.DataFrame(rows), shutdown_year, profiles


# A single structure at the middle of the reach, drift from high index
def one_groin(ny, M, f):
    updrift = ny // 2
    return [(updrift, updrift - 1, M, f)]


# The audits, in the order main() runs them

# Audit 1: diffusivity and shutdown angle for each wave climate
def audit_diffusivity():
    print("\n=== 1. DIFFUSIVITY AND SHUTDOWN ANGLE ===")
    print("The shutdown angle is a property of the ANGULAR wave distribution, not")
    print("of wave height: Hs scales the diffusivity, asym and ahf set its sign.")
    print("'offset/cell' is the shoreline offset across one 500 m domain at which")
    print("diffusion stops, which is the fillet a groin must not exceed.\n")
    header = (f"{'Hs':>5} {'Tp':>5} {'asym':>6} {'ahf':>5} {'D(0) m2/yr':>12} "
              f"{'r_ipl(0)':>9} {'shutdown':>10} {'offset/cell':>12}")
    rows = []

    # Print and record one climate's row
    def line(Hs, Tp, asym, ahf):
        cd, _, _ = brie_diffusivity(Hs, Tp, asym, ahf, 41)
        lo, hi = shutdown_angle_deg(cd)
        r0 = cd[90] * DT_YR / 2.0 / DY_M**2
        offset = np.tan(np.deg2rad(lo)) * DY_M
        print(f"{Hs:5.1f} {Tp:5.1f} {asym:6.2f} {ahf:5.2f} {cd[90]:12.0f} "
              f"{r0:9.3f} {lo:9.0f}d {offset:11.0f}m")
        rows.append(dict(Hs=Hs, Tp=Tp, asym=asym, ahf=ahf, D0=float(cd[90]),
                         r_ipl_0=r0, shutdown_angle_neg_deg=lo,
                         shutdown_angle_pos_deg=hi, shutdown_offset_m=offset))

    print("wave height, at the default angular distribution")
    print(header)
    for Hs in (1.0, 1.5, 2.0, 2.5, 3.0, 3.5):
        line(Hs, 9.0, 0.8, 0.2)

    print("\nangular distribution, at the Hatteras wave height")
    print(header)
    for ahf in (0.1, 0.2, 0.3, 0.4, 0.5):
        line(2.5, 9.0, 0.8, ahf)
    for asym in (0.5, 0.6, 0.7, 0.9):
        line(2.5, 9.0, asym, 0.2)

    return pd.DataFrame(rows)


# Audit 2: fillet against M across wave heights, flagging the runaway boundary
def audit_amplitude_and_runaway(years, ny=41, f=0.6):
    print(f"\n=== 2. FILLET vs M, AND THE RUNAWAY BOUNDARY ({ny} domains, "
          f"{years} yr, f={f}) ===")
    print("'RUN@yr' is the year the diffusion number first hit zero. The fillet")
    print("for those cells is meaningless: it is still growing when the run ends.\n")
    M_values = (10, 20, 40, 60, 80, 120)
    print(f"{'Hs':>5} {'r_ipl':>7} " + " ".join(f"{'M=' + str(M):>10}" for M in M_values))
    rows = []
    for Hs in (1.0, 1.5, 2.0, 2.5, 3.0, 3.5):
        climate = dict(Hs=Hs, Tp=9.0, asym=0.8, ahf=0.2)
        cd, _, _ = brie_diffusivity(Hs, 9.0, 0.8, 0.2, ny)
        r0 = cd[90] * DT_YR / 2.0 / DY_M**2
        cells = []
        for M in M_values:
            frame, shutdown, _ = solve_reach(ny, one_groin(ny, M, f), years, climate)
            fillet = float(frame.fillet_m.iloc[-1])
            cells.append(f"{'RUN@' + str(shutdown):>10}" if shutdown
                         else f"{fillet:9.0f}m")
            rows.append(dict(Hs=Hs, r_ipl_0=r0, M=M, f=f, ny=ny, years=years,
                             fillet_m=fillet, shutdown_year=shutdown,
                             fillet_over_M_over_r=fillet / (M / r0)))
        print(f"{Hs:5.1f} {r0:7.3f} " + " ".join(cells))

    frame = pd.DataFrame(rows)
    ratio = frame.loc[frame.shutdown_year.isna(), "fillet_over_M_over_r"]
    print(f"\nStable cells collapse onto fillet = C * M / r_ipl, "
          f"C = {ratio.mean():.3f} +/- {ratio.std():.3f} (n = {len(ratio)}).")
    print("So M and the wave climate are NOT independent knobs: the fillet is")
    print("bought by the ratio M / r_ipl, and r_ipl scales as Hs^2.4. A fitted M")
    print("is a statement about this wave climate and no other.")
    return frame


# Audit 3: fillet and volume closure against f, at fixed M
def audit_sink_fraction(years, ny=41, M=60):
    print(f"\n=== 3. FILLET AND VOLUME CLOSURE vs f ({ny} domains, {years} yr, "
          f"M={M}, Hatteras climate) ===")
    print("f scales the downdrift sink only, so f < 1 makes the groin a NET")
    print("SOURCE and the reach mean must advance seaward. 'closure_err' is the")
    print("part of the drift the scheme invents: it should be zero and is not.\n")
    print(f"{'f':>5} {'fillet_m':>10} {'mean_m':>10} {'conservative':>13} "
          f"{'closure_err':>12} {'err_m_per_yr':>13}")
    rows = []
    for f in (0.0, 0.2, 0.4, 0.6, 0.8, 1.0):
        frame, shutdown, _ = solve_reach(ny, one_groin(ny, M, f), years,
                                         HATTERAS_CLIMATE)
        last = frame.iloc[-1]
        print(f"{f:5.1f} {last.fillet_m:10.1f} {last.mean_m:10.1f} "
              f"{last.mean_if_conservative_m:13.1f} {last.closure_error_m:12.1f} "
              f"{last.closure_error_m / years:13.3f}")
        rows.append(dict(M=M, f=f, ny=ny, years=years, shutdown_year=shutdown,
                         **last.to_dict()))
    return pd.DataFrame(rows)


# Audit 4: how the fillet and the mean drift depend on reach length
def audit_domain_count(years, M=60, f=0.6):
    print(f"\n=== 4. DOMAIN-COUNT CONVERGENCE ({years} yr, M={M}, f={f}, "
          f"Hatteras climate) ===")
    print("The solve is PERIODIC. A short reach wraps the fillet into its own")
    print("downdrift notch, and spreads the net source over fewer cells, so the")
    print("spurious mean drift is inflated as 1/ny. This is the whole argument")
    print("for the rig's width, and it owes nothing to Hatteras.\n")
    print(f"{'ny':>5} {'reach_km':>9} {'fillet_m':>10} {'vs ny=121':>10} "
          f"{'mean_m':>10} {'mean_m_per_yr':>14}")
    rows = []
    for ny in (5, 7, 11, 21, 31, 40, 61, 81, 121):
        frame, shutdown, _ = solve_reach(ny, one_groin(ny, M, f), years,
                                         HATTERAS_CLIMATE)
        last = frame.iloc[-1]
        rows.append(dict(ny=ny, reach_km=ny * DY_M / 1000.0, years=years, M=M, f=f,
                         shutdown_year=shutdown, **last.to_dict()))
    reference = rows[-1]["fillet_m"]
    for row in rows:
        bias = 100.0 * (row["fillet_m"] - reference) / reference
        print(f"{row['ny']:5d} {row['reach_km']:9.1f} {row['fillet_m']:10.1f} "
              f"{bias:9.1f}% {row['mean_m']:10.1f} "
              f"{row['mean_m'] / years:14.3f}")
    return pd.DataFrame(rows)


# Audit 5: several structures, stacked, spaced and opposed: do they superpose, and when do they merge?
def audit_groin_field(years, ny=121, M=60, f=0.6):
    print(f"\n=== 5. GROIN FIELDS ({ny} domains, {years} yr, M={M}, f={f}, "
          f"Hatteras climate) ===")
    rows = []

    print("\nSTACKED -- N dipoles on one domain pair, against one dipole of N*M.")
    print("Exact superposition would make these two columns equal.\n")
    print(f"{'N':>4} {'N*M':>6} {'stacked fillet':>15} {'single M=N*M':>14} "
          f"{'departure':>10}")
    for n in (1, 2, 4, 8):
        updrift = ny // 2
        stacked = [(updrift, updrift - 1, M, f)] * n
        frame_s, shut_s, _ = solve_reach(ny, stacked, years, HATTERAS_CLIMATE)
        frame_1, shut_1, _ = solve_reach(ny, one_groin(ny, n * M, f), years,
                                        HATTERAS_CLIMATE)
        a = float(frame_s.fillet_m.iloc[-1])
        b = float(frame_1.fillet_m.iloc[-1])
        flag = "  (runaway)" if (shut_s or shut_1) else ""
        print(f"{n:4d} {n * M:6d} {a:15.1f} {b:14.1f} "
              f"{100.0 * (a - b) / b:9.2f}%{flag}")
        rows.append(dict(case="stacked", n_groins=n, spacing=0, M=M, f=f, ny=ny,
                         years=years, fillet_m=a, equivalent_single_fillet_m=b,
                         shutdown_year=shut_s))

    print("\nSPACED -- N dipoles every `spacing` domains, all the same sense.")
    print("'end' is the most updrift structure, 'interior' the middle one. When")
    print("interior/end is near zero the field holds ONE fillet and the inner")
    print("structures are inert; near one, each structure holds its own.\n")
    print(f"{'N':>4} {'spacing':>8} {'km apart':>9} {'end fillet':>11} "
          f"{'interior':>10} {'ratio':>7} {'reach mean':>11}")
    for n_groins in (3, 5):
        for spacing in (1, 2, 4, 8, 16):
            first = ny // 2 - (n_groins // 2) * spacing
            field = [(first + i * spacing, first + i * spacing - 1, M, f)
                     for i in range(n_groins)]
            if min(d for _, d, _, _ in field) < 0 or max(u for u, _, _, _ in field) >= ny:
                continue
            frame, shutdown, _ = solve_reach(ny, field, years, HATTERAS_CLIMATE)
            last = frame.iloc[-1]
            end = float(last["fillet_0_m"])
            interior = float(last[f"fillet_{n_groins // 2}_m"])
            print(f"{n_groins:4d} {spacing:8d} {spacing * DY_M / 1000:9.1f} "
                  f"{end:11.1f} {interior:10.1f} {interior / end:7.2f} "
                  f"{last.mean_m:11.1f}")
            rows.append(dict(case="spaced", n_groins=n_groins, spacing=spacing,
                             M=M, f=f, ny=ny, years=years, fillet_m=end,
                             interior_fillet_m=interior,
                             interior_over_end=interior / end,
                             mean_m=float(last.mean_m),
                             shutdown_year=shutdown))

    print("\nOPPOSED -- two structures facing each other (the second dipole")
    print("reversed), which is what a field looks like either side of a")
    print("divergence in the drift.\n")
    print(f"{'spacing':>8} {'km apart':>9} {'fillet A':>10} {'fillet B':>10} "
          f"{'reach mean':>11}")
    for spacing in (2, 4, 8, 16):
        a_up = ny // 2
        b_up = a_up + spacing
        if b_up >= ny:
            continue
        field = [(a_up, a_up - 1, M, f), (b_up, b_up + 1, M, f)]
        frame, shutdown, _ = solve_reach(ny, field, years, HATTERAS_CLIMATE)
        last = frame.iloc[-1]
        print(f"{spacing:8d} {spacing * DY_M / 1000:9.1f} "
              f"{last['fillet_0_m']:10.1f} {last['fillet_1_m']:10.1f} "
              f"{last.mean_m:11.1f}")
        rows.append(dict(case="opposed", n_groins=2, spacing=spacing, M=M, f=f,
                         ny=ny, years=years, fillet_m=float(last["fillet_0_m"]),
                         interior_fillet_m=float(last["fillet_1_m"]),
                         mean_m=float(last.mean_m), shutdown_year=shutdown))

    return pd.DataFrame(rows)


# Audit 6: the rig as designed (20 working domains, n_buffer per side): buffer size and runaway boundary
def audit_chosen_rig(years=200, ny=40, n_buffer=10, f=0.6):
    print(f"\n=== 6. THE CHOSEN RIG: {ny} domains "
          f"({ny - 2 * n_buffer} working + {n_buffer} buffer per side), "
          f"{years} yr, f={f} ===")
    print("'reach' is the dipole's diffusive reach in domains at the end of the")
    print("run, which depends on the wave climate and NOT on M. 'buffer?' is")
    print("whether that reach stays inside the buffer -- where it does not, the")
    print("fillet's tail has wrapped, and the extent test needs a wider grid even")
    print("though the fillet itself is still good to a few percent.\n")
    M_values = (10, 20, 40, 60, 80)
    print(f"{'Hs':>5} {'r_ipl':>7} {'reach':>7} {'buffer?':>8} "
          + " ".join(f"{'M=' + str(M):>10}" for M in M_values))
    rows = []
    for Hs in (1.0, 1.5, 2.0, 2.5, 3.0):
        climate = dict(Hs=Hs, Tp=9.0, asym=0.8, ahf=0.2)
        cd, _, _ = brie_diffusivity(Hs, 9.0, 0.8, 0.2, ny)
        r0 = cd[90] * DT_YR / 2.0 / DY_M**2
        reach = np.sqrt(2.0 * r0 * years)
        ok = "yes" if reach <= n_buffer else "NO"
        cells = []
        for M in M_values:
            frame, shutdown, _ = solve_reach(ny, one_groin(ny, M, f), years, climate)
            fillet = float(frame.fillet_m.iloc[-1])
            cells.append(f"{'RUN@' + str(shutdown):>10}" if shutdown
                         else f"{fillet:9.0f}m")
            rows.append(dict(Hs=Hs, r_ipl_0=r0, M=M, f=f, ny=ny,
                             n_buffer=n_buffer, years=years,
                             M_over_r=M / r0, fillet_m=fillet,
                             reach_domains=reach, buffer_ok=(reach <= n_buffer),
                             shutdown_year=shutdown))
        print(f"{Hs:5.1f} {r0:7.3f} {reach:7.1f} {ok:>8} " + " ".join(cells))
    print("\nThe buffer is exceeded from the wave height where reach > "
          f"{n_buffer}; at 200 years that is the time")
    print("sqrt(2*r_ipl*t) crosses the buffer width, so it is a property of the")
    print("wave climate and the run length, never of the groin.")
    return pd.DataFrame(rows)


# Run: all six audits, then one CSV per audit
def main():
    parser = argparse.ArgumentParser(
        description="Audit BRIE's alongshore solve under a groin dipole.")
    parser.add_argument("--years", type=int, default=200,
                        help="run length for the sweeps (default 200)")
    parser.add_argument("--outdir", type=Path, default=HERE,
                        help="where to write the CSV tables")
    args = parser.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)

    # Run every audit
    tables = {
        "diffusivity": audit_diffusivity(),
        "amplitude_runaway": audit_amplitude_and_runaway(args.years),
        "sink_fraction": audit_sink_fraction(args.years),
        "domain_count": audit_domain_count(args.years),
        "groin_field": audit_groin_field(args.years),
        "chosen_rig": audit_chosen_rig(args.years),
    }
    # Write one CSV per audit
    print()
    for name, frame in tables.items():
        path = args.outdir / f"solver_audit_{name}.csv"
        frame.to_csv(path, index=False)
        print(f"wrote {path}")


if __name__ == "__main__":
    main()
