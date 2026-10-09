#!/usr/bin/env python3
"""
A Buxton-sized groin field on a 25 m one-line model, against the four groin modules on BRIE's 500 m grid.

    python fine_grid_reference.py [--years 50]   ->  tables/, figures/

The reference: a one-line shoreline model in flux form on 25 m cells, driven by the same wave-climate
drift BRIE uses (its diffusivity, clipped at zero as BRIE clips it, integrated into a drift curve
anchored at the real net drift of the undisturbed coast). Three groins 200 m apart, the downdrift one
on the coarse groin face and the field inside the updrift domain, as at Buxton. Each groin blocks
transport GENESIS-style: it starts by holding the modules' blocking fraction of the drift and passes
more sand round its tip as the fillet grows toward it, all of it once the fillet reaches the tip.
The reference is averaged to 500 m domains and compared with each module run with one groin on the
coarse face (../1-straight-coast/straight_coast_test.py and ../4-groin-field/groin_field_test.py,
imported unchanged). Details: README.md beside this script.

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-10-09
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path[:0] = [str(HERE.parent / "1-straight-coast"), str(HERE.parent / "4-groin-field")]
import straight_coast_test as sc  # noqa: E402
import groin_field_test as gf  # noqa: E402

# --- CONFIG ------------------------------------------------------------------
NY = 61                                   # coarse domains, as the straight-coast rig
FINE_DX = 25.0                            # m
SUB = int(sc.DY_M / FINE_DX)              # fine cells per domain
LO = sc.LO                                # the coarse groin face, between domains LO and LO + 1
GROIN_OFFSETS_M = (0.0, 200.0, 400.0)     # the field's groins updrift of the coarse face (Buxton: ~200 m apart)
GROIN_LENGTH_M = 110.0                    # Buxton's digitized groins are 90-155 m
LENGTH_SENSITIVITY_M = (55.0, 220.0)
ORIENTATIONS = (-30, -20, -10, 0, 10, 20, 30)
MODULES = ("dipole", "pinned", "drift", "corrected")
RECORD_YEARS = tuple(range(1, 51))
CFL = 0.4                                 # explicit step: dt = CFL dx^2 / max diffusivity
WINDOW = 10                               # domains either side of the face in the profile error
# -----------------------------------------------------------------------------


# The drift curve the reference uses: BRIE's diffusivity clipped at zero, integrated, anchored at theta0
def drift_curve(t, theta0):
    th = np.arange(-89.0, 89.01, 0.1)
    d = np.maximum(0.0, t["coast_diff"][sc.angle_index(th, t["climl"])])
    q = np.concatenate([[0.0], np.cumsum(0.5 * (d[1:] + d[:-1]) * np.diff(np.deg2rad(th)))]) * t["depth"]
    q += float(sc.drift(t, theta0)[0]) - np.interp(theta0, th, q)
    return th, q


# One reference run; returns {year: fine shoreline, seaward +}
def run_reference(t, theta0, years, length_m=GROIN_LENGTH_M, groins=True):
    th_tab, q_tab = drift_curve(t, theta0)
    n = NY * SUB
    slope = np.tan(np.deg2rad(theta0))
    x = np.zeros(n)                                       # landward +, as BRIE's eta
    up_is_hi = float(sc.drift(t, theta0)[0]) < 0
    face0 = (LO + 1) * SUB - 1                            # fine face on the coarse face
    sign = 1 if up_is_hi else -1
    faces = [face0 + sign * int(round(o / FINE_DX)) for o in GROIN_OFFSETS_M] if groins else []
    dt = CFL * FINE_DX ** 2 / t["coast_diff"].max()
    steps_per_year = int(np.ceil(1.0 / dt))
    dt = 1.0 / steps_per_year
    out = {}
    for year in range(1, years + 1):
        for _ in range(steps_per_year):
            step = slope * FINE_DX + (np.roll(x, -1) - x)
            theta = np.degrees(np.arctan2(step, FINE_DX))
            q = np.interp(theta, th_tab, q_tab)          # flux through face j toward higher index
            for j in faces:
                cell_up = j + 1 if up_is_hi else j      # the cell just updrift of the groin
                approach = q[j + 1] if up_is_hi else q[j - 1]
                advance = max(0.0, -x[cell_up])          # how far the fillet has built toward the tip
                blocked = sc.BLOCKING_FRACTION * max(0.0, 1.0 - advance / length_m) ** (2.0 / 3.0)
                q[j] = (1.0 - blocked) * approach
            x -= (np.roll(q, 1) - q) * dt / (FINE_DX * t["depth"])
        if year in RECORD_YEARS:
            out[year] = -x.copy()
    return out


# Fine cells averaged to domains
def to_domains(fine):
    return fine.reshape(NY, SUB).mean(axis=1)


# Run: the reference at three groin lengths, the four modules, and the comparison tables
def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--years", type=int, default=50)
    a = ap.parse_args()
    sc.NY = gf.NY = NY
    t = sc.brie_tables()
    years = [y for y in RECORD_YEARS if y <= a.years]
    up, down = LO + 1, LO                                 # drift runs toward lower index at every orientation here

    rows, prof = [], []
    for th in ORIENTATIONS:
        assert float(sc.drift(t, th)[0]) < 0
        for length in (GROIN_LENGTH_M,) + LENGTH_SENSITIVITY_M:
            ref = run_reference(t, th, a.years, length)
            for y in years:
                c = to_domains(ref[y])
                rows.append(dict(source="reference", groin_length_m=length, theta0_deg=th, year=y,
                                 updrift_m=c[up], downdrift_m=c[down], gap_m=c[up] - c[down],
                                 net_change_m=float(ref[y].sum() / SUB)))
            if length == GROIN_LENGTH_M:
                prof.append(pd.DataFrame(dict(source="reference_fine", theta0_deg=th,
                                              domain=np.arange(NY * SUB) / SUB, seaward_m=ref[a.years])))
                prof.append(pd.DataFrame(dict(source="reference", theta0_deg=th, domain=np.arange(NY) + 0.5,
                                              seaward_m=to_domains(ref[a.years]))))
        for m in MODULES:
            for y in years:
                c = gf.run_field(t, th, m, [LO], y)
                rows.append(dict(source=m, groin_length_m=np.nan, theta0_deg=th, year=y,
                                 updrift_m=c[up], downdrift_m=c[down], gap_m=c[up] - c[down],
                                 net_change_m=float(c.sum())))
                if y == a.years:
                    prof.append(pd.DataFrame(dict(source=m, theta0_deg=th, domain=np.arange(NY) + 0.5,
                                                  seaward_m=c)))
        print(f"{th:+d} deg done")
    frame = pd.DataFrame(rows)
    prof = pd.concat(prof, ignore_index=True)

    # How far each module sits from the reference: the gap over time and the final profile
    ref = frame[(frame.source == "reference") & (frame.groin_length_m == GROIN_LENGTH_M)].set_index(["theta0_deg", "year"])
    ref_prof = prof[prof.source == "reference"].set_index(["theta0_deg", "domain"]).seaward_m
    near = np.arange(LO - WINDOW, LO + WINDOW + 2) + 0.5
    score = []
    for (m, th), s_ in frame[frame.source.isin(MODULES)].groupby(["source", "theta0_deg"]):
        s_ = s_.set_index("year")
        r_ = ref.loc[th]
        mp = prof[(prof.source == m) & (prof.theta0_deg == th)].set_index("domain").seaward_m
        score.append(dict(module=m, theta0_deg=th,
                          ref_gap_end_m=r_.gap_m.iloc[-1], module_gap_end_m=s_.gap_m.iloc[-1],
                          gap_rmse_m=float(np.sqrt(((s_.gap_m - r_.gap_m) ** 2).mean())),
                          profile_rmse_m=float(np.sqrt(((mp.loc[near] - ref_prof.loc[th].loc[near]) ** 2).mean())),
                          ref_net_m=r_.net_change_m.iloc[-1], module_net_m=s_.net_change_m.iloc[-1]))
    score = pd.DataFrame(score)

    out = HERE / "tables"
    out.mkdir(parents=True, exist_ok=True)
    frame.to_csv(out / "fine_reference_by_year.csv", index=False, float_format="%.3f")
    prof.to_csv(out / "fine_reference_profiles.csv", index=False, float_format="%.3f")
    score.to_csv(out / "fine_reference_scores.csv", index=False, float_format="%.3f")

    pd.set_option("display.width", 220)
    print(score.pivot(index="theta0_deg", columns="module", values="gap_rmse_m").round(1).to_string())
    print(score.pivot(index="theta0_deg", columns="module", values="profile_rmse_m").round(1).to_string())
    print(score.pivot(index="theta0_deg", columns="module", values="module_gap_end_m").round(1)
          .join(score.groupby("theta0_deg").ref_gap_end_m.first().round(1)).to_string())
    print(frame[(frame.source == "reference") & (frame.year.isin([5, 10, 25, 50]))]
          .pivot_table(index=["theta0_deg", "groin_length_m"], columns="year", values="gap_m").round(1).to_string())
    print("reference net change, max |m|:", frame[frame.source == "reference"].net_change_m.abs().max())


if __name__ == "__main__":
    main()
