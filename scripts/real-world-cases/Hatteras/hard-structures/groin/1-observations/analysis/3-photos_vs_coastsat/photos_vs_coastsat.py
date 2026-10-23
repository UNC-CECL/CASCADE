#!/usr/bin/env python3
"""
The aerial-photo gap against the CoastSat gap across the Buxton groin.

    python photos_vs_coastsat.py   ->  tables/photo_vs_coastsat.csv, tables/coastsat_sides_on_1967_datum.csv,
                                       tables/boundary_gaps_on_1967_datum.csv, tables/placement_targets.csv,
                                       tables/photo_minus_coastsat_by_side.csv, tables/coastsat_one_day_spread.csv,
                                       tables/era_rates_both.csv  (figures: photos_vs_coastsat_figures.py)

Reads the tables the two single-source analyses write: ../1-aerial_photo_gap/tables/
(run photo_gap.py first) and ../2-coastsat_gap/tables/ (run
coastsat_gap_across_groins.py first). The photos measure change since 1967 against
a fixed offshore datum; CoastSat each transect about its own 1984-2016 mean. Their
datums differ, so the photo gap is shifted onto the CoastSat whole-domain gap by the
mean difference over the shared pre-fill survey years, and only the shapes compare.
The reverse puts CoastSat on the photos' 1967 datum: each CoastSat side is shifted
onto the same side's photo positions over those years. CoastSat starts in 1984, 15
years into the groin's life, so only the photos can say how big the gap was.
Era rates are OLS slopes over the same eras for both records.

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-10-14
"""
from __future__ import annotations

# Non-code files beside this script live in the workspace on D: (hatteras_paths.workspace_twin)
import sys as _hat_sys  # noqa: E402
from pathlib import Path as _HatPath  # noqa: E402
_HAT_DIR = next(_p for _p in _HatPath(__file__).resolve().parents if (_p / "hatteras_paths.py").exists())
if str(_HAT_DIR) not in _hat_sys.path:
    _hat_sys.path.insert(0, str(_HAT_DIR))
from hatteras_paths import workspace_twin  # noqa: E402
import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "2-coastsat_gap"))
import coastsat_gap_across_groins as cg  # noqa: E402
from coastsat_gap_across_groins import BOUNDARIES, ERAS, FIT, era_rate  # noqa: E402

# --- CONFIG ------------------------------------------------------------------
PHOTO_TABLES = workspace_twin(HERE.parent) / "1-aerial_photo_gap" / "tables"
COASTSAT_TABLES = workspace_twin(HERE.parent) / "2-coastsat_gap" / "tables"
TABLES = workspace_twin(HERE) / "tables"
KEYS = (1996, 2009, 2025)            # calibration start and end, test end (calendar years)
# -----------------------------------------------------------------------------


# How far one day's CoastSat side position sits from that year's mean: the scatter a single-day photo
# would show against an annual mean even if the two records agreed perfectly. Per image date, each side
# is the mean of its transects' positions (each about the same reference as the annual series); per
# year, the standard deviation of those dates about their mean.
def one_day_spread():
    up, down = cg.sides()["domain"]
    days = {}
    for side, tids in (("north", up), ("south", down)):
        cols = {}
        for t in tids:
            d = pd.read_csv(cg.SERIES_DIR / f"{cg.SITE}_timeseries" / f"{t}.csv", header=0)
            d.columns = ["date", "chainage_m"] + list(d.columns[2:])
            d["day"] = pd.to_datetime(d["date"], utc=True).dt.normalize()
            d["chainage_m"] = pd.to_numeric(d["chainage_m"], errors="coerce")
            ref = cg.annual(t).loc[FIT[0]:FIT[1]].mean()
            cols[t] = d.dropna(subset=["chainage_m"]).groupby("day")["chainage_m"].mean() - ref
        a = pd.DataFrame(cols)
        ok = a.notna().sum(axis=1) >= np.ceil(cg.MIN_SHARE * len(tids))
        days[side] = a[ok].mean(axis=1)
    d = pd.DataFrame(days).dropna()
    d["gap"] = d["north"] - d["south"]
    d["year"] = d.index.year
    g = d.groupby("year")
    out = g[["north", "south", "gap"]].std().add_suffix("_one_day_sd_m")
    out["n_images"] = g.size()
    return out[out["n_images"] >= cg.MIN_OBS_YEAR]


# The gap across each candidate boundary from both records, CoastSat moved onto the photos' 1967 datum,
# and the numbers a placement run would be scored on
def placement():
    pg = pd.read_csv(workspace_twin(PHOTO_TABLES) / "photo_boundary_gaps.csv", index_col="year")
    cg_ = pd.read_csv(workspace_twin(COASTSAT_TABLES) / "boundary_gaps.csv", index_col="year")
    on, rows = pd.DataFrame(index=cg_.index), []
    for s, n in BOUNDARIES:
        k = f"gap_{s}_{n}_m"
        shared = [y for y in pg.index if FIT[0] <= y <= FIT[1] and y in cg_.index]
        shift = (pg.loc[shared, k] - cg_.loc[shared, k]).mean()
        c = cg_[k] + shift
        on[f"coastsat_{k}"] = c
        diff = (pg.loc[shared, k] - c.loc[shared])
        row = dict(boundary=f"{s}/{n}", shift_m=shift, photo_minus_coastsat_rms_m=float(np.sqrt((diff ** 2).mean())))
        for y in KEYS:
            row[f"coastsat_{y}_m"] = c.get(y, np.nan)
            row[f"photo_{y}_m"] = pg[k].get(y, np.nan)
        row[f"coastsat_change_{KEYS[0]}_{KEYS[1]}_m"] = c.get(KEYS[1], np.nan) - c.get(KEYS[0], np.nan)
        row[f"coastsat_change_{KEYS[1]}_{KEYS[2]}_m"] = c.get(KEYS[2], np.nan) - c.get(KEYS[1], np.nan)
        for lo, hi in ERAS:
            r, se, _ = era_rate(cg_[k], lo, hi)
            row[f"coastsat_rate_{lo}_{hi}_m_yr"], row[f"coastsat_rate_{lo}_{hi}_se"] = r, se
        row["photo_peak_m"], row["photo_peak_year"] = pg[k].max(), int(pg[k].idxmax())
        rows.append(row)
    return on.rename_axis("year"), pd.DataFrame(rows)


def main():
    TABLES.mkdir(exist_ok=True)
    on, targets = placement()
    on.to_csv(workspace_twin(TABLES) / "boundary_gaps_on_1967_datum.csv", float_format="%.2f")
    targets.set_index("boundary").round(2).T.to_csv(workspace_twin(TABLES) / "placement_targets.csv")
    photos = pd.read_csv(workspace_twin(PHOTO_TABLES) / "photo_gap_annual.csv", index_col="year")["gap_m"]
    a = pd.read_csv(workspace_twin(COASTSAT_TABLES) / "gap_annual.csv")
    coastsat = a[a["gap"] == "domain"].set_index("year")["gap_m"]

    # The photo gap shifted onto the CoastSat whole-domain gap over the shared pre-fill years
    shared = [y for y in photos.index if FIT[0] <= y <= FIT[1] and y in coastsat.index]
    shift = coastsat.loc[shared].mean() - photos.loc[shared].mean()
    ph = pd.DataFrame({"year": photos.index, "gap_change_since_1967_m": photos.values,
                       "shifted_m": photos.values + shift})
    ph["coastsat_domain_m"] = ph["year"].map(coastsat)
    ph["photo_minus_coastsat_m"] = ph["shifted_m"] - ph["coastsat_domain_m"]
    ph.to_csv(workspace_twin(TABLES) / "photo_vs_coastsat.csv", index=False, float_format="%.2f")

    # Each CoastSat side shifted onto the same side's photo positions over the shared pre-fill years,
    # so both sides (and their difference, the gap) read on the photos' 1967 datum
    ps = pd.read_csv(workspace_twin(PHOTO_TABLES) / "photo_gap_annual.csv", index_col="year")
    d = a[a["gap"] == "domain"].set_index("year")
    sides = pd.DataFrame(index=d.index.rename("year"))
    for side, col in (("north", "updrift_anom_m"), ("south", "downdrift_anom_m")):
        sides[f"{side}_m"] = d[col] + (ps.loc[shared, f"{side}_m"] - d.loc[shared, col]).mean()
    sides["gap_m"] = sides["north_m"] - sides["south_m"]
    sides.to_csv(workspace_twin(TABLES) / "coastsat_sides_on_1967_datum.csv", float_format="%.2f")

    # Photo minus CoastSat per side and for the gap, on the 1967 datum, in each photo year CoastSat covers
    res = pd.DataFrame({f"{k}_photo_minus_coastsat_m": ps[f"{k}_m"] - sides[f"{k}_m"]
                        for k in ("north", "south", "gap")}).dropna().rename_axis("year")
    res.to_csv(workspace_twin(TABLES) / "photo_minus_coastsat_by_side.csv", float_format="%.2f")
    one_day_spread().to_csv(workspace_twin(TABLES) / "coastsat_one_day_spread.csv", float_format="%.2f")

    # Both records' rates over the same eras (CoastSat's from its own table, the photos' here)
    e = pd.read_csv(workspace_twin(COASTSAT_TABLES) / "era_rates.csv")
    eras = e[e["source"] == "domain"].to_dict("records")
    for lo, hi in ERAS:
        r, se, k = era_rate(photos, lo, hi)
        eras.append(dict(source="photos", era=f"{lo}-{hi}", rate_m_yr=r, se_m_yr=se, n_years=k))
    pd.DataFrame(eras).to_csv(workspace_twin(TABLES) / "era_rates_both.csv", index=False, float_format="%.3f")
    print(f"shift {shift:+.2f} m over {len(shared)} shared years")
    print(ph.round(1).to_string(index=False))
    print(pd.DataFrame(eras).round(2).to_string(index=False))


if __name__ == "__main__":
    main()
