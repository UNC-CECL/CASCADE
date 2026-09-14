#!/usr/bin/env python3
"""Stage 0 of the groin module test: audit BRIE's solve without running CASCADE.

WHY THIS EXISTS
    Every claim about what a groin "does" in CASCADE is a claim about what
    BRIE's implicit alongshore-diffusion solve does with the +/-M dipole that
    `cascade.groin.GroinCallback` injects into `x_s_dt`. That solve is cheap:
    one sparse tridiagonal-plus-corners system per year. Reproducing it here
    WITHOUT Barrier3D means a 1,000-year, many-axis experiment costs seconds
    instead of days, and -- more importantly -- it separates the groin's
    alongshore behaviour from the cross-shore behaviour Barrier3D adds on top.
    Anything the full rig shows that this does not is a Barrier3D feedback.

    This is an EMULATOR, not a reimplementation for production use. It mirrors
    `brie.brie.Brie.update()` (the `if self._ast_model_on:` block, lines
    ~1290-1330) exactly: the same `coast_diff` table, taken from a real Brie
    instance rather than recomputed, the same row-scaled matrix assembly from
    `_di`/`_dj`, the same periodic wrap, the same `np.maximum(0, ...)` clip on
    the diffusion number. If BRIE changes, this must be re-checked against it.

WHAT IT MEASURES
    1. DIFFUSIVITY AND SHUTDOWN ANGLE -- BRIE's diffusion number is clipped at
                                         zero, and the wave-climate-averaged
                                         diffusivity goes NEGATIVE past a
                                         critical shoreline angle. Past that
                                         angle a cell stops exchanging sand
                                         with its neighbours and the dipole
                                         accumulates without limit.
    2. FILLET AMPLITUDE vs M          -- is the response linear in M, as
                                         `groin.predict_fillet` asserts, and
                                         where is the runaway boundary?
    3. VOLUME CLOSURE vs f            -- the matrix is row-scaled by a
                                         per-domain diffusion number, so its
                                         column sums are not unity and the
                                         scheme is NOT exactly conservative
                                         once the shoreline is not straight.
                                         This reports the spurious mean drift
                                         against the drift the injected volume
                                         actually implies.
    4. DOMAIN-COUNT CONVERGENCE       -- the solve is PERIODIC in the
                                         alongshore, so a short reach wraps the
                                         fillet into its own downdrift notch.
                                         This reports how few domains the
                                         fillet tolerates and how badly the
                                         mean drift is inflated by a short one.
    5. GROIN FIELDS                   -- several dipoles at a given spacing.
                                         Two questions: do they superpose (a
                                         field of N reads as one groin of N*M),
                                         and at what spacing does each
                                         structure still hold its own fillet
                                         rather than the field holding one?

HOW TO READ "RUNAWAY"
    A cell whose diffusion number has been clipped to zero keeps receiving its
    share of the dipole and has no way to pass it on, so the reach translates at
    roughly M metres per year indefinitely. In the full model this presents as a
    barrier that migrates absurdly or drowns; in BRIE alone it eventually
    presents as `IndexError: index 180 is out of bounds` from brie.py, because
    the diffusivity lookup indexes `coast_diff` (length 180) with `90 - theta`
    clipped to 180 rather than 179. That needs a 57 km offset across one 500 m
    cell, which a runaway reaches in roughly 1,000 years at M = 60. It is the
    runaway surfacing, not a separate coding mistake.

Author: Hannah A. Henry, UNC CECL
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.sparse import csr_matrix
from scipy.sparse.linalg import spsolve

from brie import Brie

HERE = Path(__file__).resolve().parent

# BRIE fixes both of these in the CASCADE coupler, whose comment says "do not
# change", so they are constants here rather than arguments.
DY_M = 500.0
DT_YR = 1.0

# Wave climate. The first is the CASCADE / Barrier3D default, the second is the
# Hatteras hindcast setting. They are NOT interchangeable for this test: see
# the diffusivity table, where they differ by a factor of four in the restoring
# rate and therefore in the M a given fillet costs.
DEFAULT_CLIMATE = dict(Hs=1.0, Tp=7.0, asym=0.8, ahf=0.2)
HATTERAS_CLIMATE = dict(Hs=2.5, Tp=9.0, asym=0.8, ahf=0.2)

_CD_CACHE: dict = {}


def brie_diffusivity(Hs, Tp, asym, ahf, ny):
    """Return (coast_diff, di, dj) from a real Brie instance.

    The diffusivity table and the sparse index arrays are taken from BRIE
    rather than recomputed, so this audit cannot quietly drift away from the
    model it is auditing.
    """
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


def shutdown_angle_deg(coast_diff):
    """Shoreline angles bounding the band where BRIE's diffusivity is positive.

    Outside that band the diffusion number is clipped to zero and the cell
    decouples from its neighbours. Returns (negative_cutoff, positive_cutoff)
    in degrees.
    """
    theta = np.arange(-89, 90)
    idx = np.clip(np.round(90 - theta).astype(int), 1, 179)
    positive = theta[coast_diff[idx] > 0]
    return float(positive.min()), float(positive.max())


def solve_reach(ny, groins, years, climate, record_years=()):
    """Integrate BRIE's shoreline diffusion with groin dipoles, nothing else.

    Starts from a perfectly straight shoreline at x_s = 0, so every metre of
    structure in the answer came from a dipole.

    Parameters
    ----------
    ny : int
        Domain count. The reach is PERIODIC, so this is a circumference.
    groins : list of (updrift, downdrift, M, f)
        One tuple per structure. `updrift` and `downdrift` are domain indices
        and must be adjacent, matching `GroinCallback`'s own check. Sign
        follows `cascade.groin`: the updrift cell gets -M (seaward advance),
        the downdrift cell gets +M*f (landward retreat).
    years : int
        Run length.
    climate : dict
        Keys Hs, Tp, asym, ahf.
    record_years : iterable of int
        Years whose full shoreline profile to keep.

    Returns
    -------
    frame : DataFrame
        One row per year: fillet at each structure, reach mean, the mean a
        conservative scheme would give, the closure error, and the minimum
        diffusion number anywhere.
    shutdown_year : int or None
        First year the diffusion number hit zero anywhere.
    profiles : dict
        Shoreline profiles for `record_years`.
    """
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
    # The net source the field injects per year, seaward-negative. A field with
    # every f = 1 is volume neutral and must leave the reach mean alone.
    net_source = -sum(M * (1.0 - f) for _, _, M, f in groins)

    shutdown_year = None
    rows, profiles = [], {}

    for year in range(1, years + 1):
        # BRIE's forward-difference shoreline angle, and the diffusion number it
        # selects. The clip at zero is BRIE's, and it is what lets the scheme
        # stop diffusing rather than go unstable.
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


def one_groin(ny, M, f):
    """A single structure at the middle of the reach, drift from high index."""
    updrift = ny // 2
    return [(updrift, updrift - 1, M, f)]


# ---------------------------------------------------------------------------
# The five audits
# ---------------------------------------------------------------------------
def audit_diffusivity():
    """Diffusivity and shutdown angle for each wave climate."""
    print("\n=== 1. DIFFUSIVITY AND SHUTDOWN ANGLE ===")
    print("The shutdown angle is a property of the ANGULAR wave distribution, not")
    print("of wave height: Hs scales the diffusivity, asym and ahf set its sign.")
    print("'offset/cell' is the shoreline offset across one 500 m domain at which")
    print("diffusion stops, which is the fillet a groin must not exceed.\n")
    header = (f"{'Hs':>5} {'Tp':>5} {'asym':>6} {'ahf':>5} {'D(0) m2/yr':>12} "
              f"{'r_ipl(0)':>9} {'shutdown':>10} {'offset/cell':>12}")
    rows = []

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


def audit_amplitude_and_runaway(years, ny=41, f=0.6):
    """Fillet against M across wave heights; flag the runaway boundary."""
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


def audit_sink_fraction(years, ny=41, M=60):
    """Fillet and volume closure against f, at fixed M."""
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


def audit_domain_count(years, M=60, f=0.6):
    """How the fillet and the mean drift depend on reach length."""
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


def audit_groin_field(years, ny=121, M=60, f=0.6):
    """Several structures: do they superpose, and at what spacing do they merge?

    Two separate questions, deliberately kept apart.

    STACKED. Several dipoles on the SAME pair of domains is the case
    GROIN_PLAN calls "four groins, one dipole, deliberately" -- the real Buxton
    field fits inside one 500 m cell. In a LINEAR solve N stacked dipoles of
    amplitude M are exactly one dipole of amplitude N*M. The solve is not
    linear, because the diffusion number depends on the shoreline angle, so
    this measures how far from N*M the answer actually lands.

    SPACED. Dipoles every `spacing` domains. A field whose structures are far
    apart holds one fillet each; a field whose structures are close holds one
    fillet for the whole field, with the interior ones doing nothing. The
    interior-to-end fillet ratio says which regime a spacing is in.
    """
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


def audit_chosen_rig(years=200, ny=40, n_buffer=10, f=0.6):
    """The rig as designed: 20 working domains, `n_buffer` buffers per side.

    Two things this has to establish, because the rig's width rests on them.

    WHETHER THE BUFFER IS BIG ENOUGH. It is NOT an absorbing boundary -- BRIE's
    solve is periodic, so a buffer separates the groin from the wrap, it does
    not soak anything up. The diffusive reach of the dipole grows as
    dy*sqrt(2*r_ipl*t) with no dependence on M, so the buffer is exceeded after
    a time that depends only on the wave climate. Past that point the fillet's
    own tail has come round the back. Section 4 of this audit says that costs
    about 3% of the fillet, and rather more of the reach mean and of any attempt
    to measure the alongshore extent.

    WHERE THE RUNAWAY BOUNDARY FALLS ACROSS THE WAVE-HEIGHT RANGE. The fillet is
    bought by M / r_ipl, so sweeping wave height IS sweeping the cost of M. This
    prints the grid the CASCADE sweep should cover and flags which cells are
    past the shutdown.
    """
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


def main():
    parser = argparse.ArgumentParser(
        description="Audit BRIE's alongshore solve under a groin dipole.")
    parser.add_argument("--years", type=int, default=200,
                        help="run length for the sweeps (default 200)")
    parser.add_argument("--outdir", type=Path, default=HERE,
                        help="where to write the CSV tables")
    args = parser.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)

    tables = {
        "diffusivity": audit_diffusivity(),
        "amplitude_runaway": audit_amplitude_and_runaway(args.years),
        "sink_fraction": audit_sink_fraction(args.years),
        "domain_count": audit_domain_count(args.years),
        "groin_field": audit_groin_field(args.years),
        "chosen_rig": audit_chosen_rig(args.years),
    }
    print()
    for name, frame in tables.items():
        path = args.outdir / f"solver_audit_{name}.csv"
        frame.to_csv(path, index=False)
        print(f"wrote {path}")


if __name__ == "__main__":
    main()
