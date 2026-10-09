#!/usr/bin/env python3
"""
How the groin modules hold up over 200 years: one groin on a straight coast, against the real-sized field.

    python long_term_test.py [--years 200] [--skip-reference]   ->  tables/

The single-groin straight-coast rig (61 domains, option A waves, BRIE's alongshore solve alone) run
for --years with each module, recording every year: the gap across the groin, the net sand change
over the coast, and the steepest shoreline angle at any face (BRIE's diffusivity table ends at +-90
degrees, where a runaway stops it). The 25 m reference field with 110 m groins
(../5-fine-grid-reference/fine_grid_reference.py) is run alongside for the gap. All three scripts it
imports are used unchanged. Figures: long_term_figures.py. Details: README.md beside this script.

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
sys.path[:0] = [str(HERE.parent / d) for d in ("1-straight-coast", "4-groin-field", "5-fine-grid-reference")]
import straight_coast_test as sc  # noqa: E402
import groin_field_test as gf  # noqa: E402
import fine_grid_reference as fr  # noqa: E402

# --- CONFIG ------------------------------------------------------------------
NY = fr.NY                                # 61 domains, as the straight-coast rig
LO = fr.LO
ORIENTATIONS = fr.ORIENTATIONS            # -30 to +30 degrees
MODULES = ("dipole", "pinned", "drift", "corrected", "length")
# -----------------------------------------------------------------------------


# Steepest shoreline angle at any face, degrees, from a seaward-positive shoreline
def max_face_angle(theta0, x):
    eta = -x
    step = np.tan(np.deg2rad(theta0)) * sc.DY_M + (np.roll(eta, -1) - eta)
    return float(np.abs(np.degrees(np.arctan2(step, sc.DY_M))).max())


# Run: every module and orientation for --years, the reference alongside, one row per year
def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--years", type=int, default=200)
    ap.add_argument("--skip-reference", action="store_true")
    a = ap.parse_args()
    sc.NY = gf.NY = fr.NY = NY
    t = sc.brie_tables()
    years = range(1, a.years + 1)
    up, down = LO + 1, LO

    rows = []
    for th in ORIENTATIONS:
        for m in MODULES:
            try:
                snaps = gf.run_field(t, th, m, [LO], a.years, record_years=set(years))
            except Exception as err:                       # a runaway stops the run; record where
                rows.append(dict(source=m, theta0_deg=th, year=np.nan, failed=repr(err)))
                continue
            for y, x in snaps.items():
                rows.append(dict(source=m, theta0_deg=th, year=y, gap_m=x[up] - x[down], updrift_m=x[up],
                                 downdrift_m=x[down], net_change_m=float(x.sum()),
                                 max_angle_deg=max_face_angle(th, x), finite=bool(np.isfinite(x).all())))
        if not a.skip_reference:
            fr.RECORD_YEARS = tuple(years)
            ref = fr.run_reference(t, th, a.years, fr.GROIN_LENGTH_M)
            for y, xf in ref.items():
                x = fr.to_domains(xf)
                rows.append(dict(source="reference", theta0_deg=th, year=y, gap_m=x[up] - x[down], updrift_m=x[up],
                                 downdrift_m=x[down], net_change_m=float(xf.sum() / fr.SUB),
                                 max_angle_deg=np.nan, finite=bool(np.isfinite(xf).all())))
        print(f"{th:+d} deg done", flush=True)
    frame = pd.DataFrame(rows)
    out = HERE / "tables"
    out.mkdir(parents=True, exist_ok=True)
    frame.to_csv(out / "long_term_by_year.csv", index=False, float_format="%.3f")

    # Summary: the gap at checkpoints, its growth over the last quarter, net change, steepest angle
    f = frame.dropna(subset=["year"])
    last, q = a.years, a.years - a.years // 4
    piv = f[f.year.isin([50, 100, last])].pivot_table(index=["theta0_deg", "source"], columns="year", values="gap_m")
    g = f.set_index(["theta0_deg", "source", "year"]).gap_m
    summ = piv.copy()
    summ.columns = [f"gap_{c}_m" for c in summ.columns]
    summ["gap_growth_last_quarter_m"] = [g.loc[(th, s, last)] - g.loc[(th, s, q)] for th, s in summ.index]
    summ["net_change_end_m"] = f[f.year == last].set_index(["theta0_deg", "source"]).net_change_m
    summ["max_angle_deg"] = f.groupby(["theta0_deg", "source"]).max_angle_deg.max()
    summ["all_finite"] = f.groupby(["theta0_deg", "source"]).finite.all()
    summ = summ.reset_index()
    summ.to_csv(out / "long_term_summary.csv", index=False, float_format="%.3f")
    pd.set_option("display.width", 220)
    print(summ.round(1).to_string(index=False))
    if "failed" in frame:
        print(frame[frame.failed.notna()][["source", "theta0_deg", "failed"]].to_string(index=False))


if __name__ == "__main__":
    main()
