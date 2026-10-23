#!/usr/bin/env python3
"""
The gap across the Buxton groin from the aerial photos alone, 1967-2023.

    python photo_gap.py   ->  tables/photo_gap_annual.csv, tables/photo_era_rates.csv,
                              tables/photo_boundary_gaps.csv, tables/photo_alongshore_change.csv
                              (figures: photo_gap_figures.py)

Each side's position is its GIS domain's row in Change_from_wetdry_1967_D2_D12.csv
(the 24 dated wet/dry surveys, change since 1967, made by
HAT_geometric_distance_sanity_check.py), negated to seaward positive. The gap is
north (GIS 6, updrift) minus south (GIS 5, downdrift), so a rising gap is the groin
holding the north side. Two surveys in one year are averaged. Era rates are OLS
slopes over the surveys inside each era.

For the groin-placement question (does the model groin belong at the GIS 5/6 or the
6/7 boundary?) it also writes the gap across both boundaries, and each survey's change
since 1967 in 100 m strips of northing across GIS 4-8, measured as the position table
is (mean distance to the offshore datum line) but within each strip instead of each
domain.

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
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent

# --- CONFIG ------------------------------------------------------------------
PHOTOS = workspace_twin(HERE) / "Change_from_wetdry_1967_D2_D12.csv"
TABLES = workspace_twin(HERE) / "tables"
UPDRIFT, DOWNDRIFT = 6, 5            # GIS domains north and south of the groin
# The groin's stages in the photo record: building, slow growth to the 1995 peak, after the last repair
ERAS = ((1967, 1978), (1978, 1995), (1996, 2016))
BOUNDARIES = ((5, 6), (6, 7))        # the two candidate places for the model groin, (south, north) domain
STRIP_DOMAINS = (4, 8)               # the alongshore profile runs from GIS 4 to GIS 8
STRIP_M = 100.0                      # strip width in northing
STRIP_SAMPLES = 10                   # points along each strip's piece of line
# -----------------------------------------------------------------------------


# Each side's position per survey year, change since 1967, seaward positive, and the gap
def photo_sides():
    t = pd.read_csv(PHOTOS).set_index("Domain_ID")
    obs = {}
    for c in t.columns:
        m = re.match(r"change_from_wetdry_1967_wetdry_(\d{4})", c)
        if m:
            obs.setdefault(int(m.group(1)), []).append(-t.loc[[UPDRIFT, DOWNDRIFT], c].to_numpy())
    s = pd.DataFrame({y: np.mean(v, axis=0) for y, v in obs.items()}, index=["north", "south"]).T.sort_index()
    s["gap"] = s["north"] - s["south"]
    s["n_surveys"] = pd.Series({y: len(v) for y, v in obs.items()})
    return s


# OLS slope, its standard error and the count over one era
def era_rate(s, lo, hi):
    w = s.loc[lo:hi].dropna()
    t, y = w.index.to_numpy(float), w.to_numpy(float)
    if len(t) < 3:
        return np.nan, np.nan, len(t)
    X = np.column_stack([np.ones_like(t), t])
    beta, *_ = np.linalg.lstsq(X, y, rcond=None)
    r = y - X @ beta
    return beta[1], np.sqrt((r @ r) / (len(t) - 2) / ((t - t.mean()) ** 2).sum()), len(t)


# The gap across each candidate boundary: north domain minus south domain, seaward positive
def boundary_gaps():
    t = pd.read_csv(PHOTOS).set_index("Domain_ID")
    obs = {}
    for c in t.columns:
        m = re.match(r"change_from_wetdry_1967_wetdry_(\d{4})", c)
        if m:
            obs.setdefault(int(m.group(1)), []).append(
                {f"gap_{s}_{n}_m": t.loc[s, c] - t.loc[n, c] for s, n in BOUNDARIES})
    return pd.DataFrame({y: pd.DataFrame(v).mean() for y, v in obs.items()}).T.sort_index().rename_axis("year")


# Each survey's change since 1967 in 100 m strips of northing, the position table's method per strip
def alongshore_change():
    from shapely.geometry import box
    sys.path.insert(0, str(HERE))
    import HAT_geometric_distance_sanity_check as h
    doms = h.load_domain_polygons(h.DOMAIN_POLYGONS_FILE)
    datum = h.load_single_geom(h.DATUM_LINE_FILE)
    lines = h.load_shorelines(str(Path(h.INPUT_DIR) / "wet_dry_shorelines_groin.geojson"))
    y0, y1 = doms[STRIP_DOMAINS[0]].bounds[1], doms[STRIP_DOMAINS[1]].bounds[3]
    edges = np.arange(y0, y1 + 1, STRIP_M)
    x0, x1 = datum.bounds[0] - 20000, datum.bounds[0]
    rows = []
    for key, (g, sort_key) in lines.items():
        for lo, hi in zip(edges[:-1], edges[1:]):
            clip = g.intersection(box(x0, lo, x1, hi))
            parts = [p for p in getattr(clip, "geoms", [clip]) if p.geom_type == "LineString" and p.length > 0]
            if not parts:
                continue
            d = [datum.distance(p.interpolate(i / (STRIP_SAMPLES - 1), normalized=True))
                 for p in parts for i in range(STRIP_SAMPLES)]
            rows.append(dict(survey=key, year=int(sort_key), y_mid=(lo + hi) / 2, distance_m=float(np.mean(d))))
    r = pd.DataFrame(rows)
    ref = r[r["survey"] == "wetdry_1967"].set_index("y_mid")["distance_m"]
    r["change_m"] = ref.reindex(r["y_mid"]).to_numpy() - r["distance_m"]    # seaward positive
    r["domain"] = [next((d for d, p in doms.items() if p.bounds[1] <= y < p.bounds[3]), np.nan) for y in r["y_mid"]]
    return r.groupby(["year", "y_mid", "domain"], as_index=False)[["distance_m", "change_m"]].mean()


def main():
    TABLES.mkdir(exist_ok=True)
    boundary_gaps().to_csv(workspace_twin(TABLES) / "photo_boundary_gaps.csv", float_format="%.2f")
    alongshore_change().to_csv(workspace_twin(TABLES) / "photo_alongshore_change.csv", index=False, float_format="%.2f")
    s = photo_sides()
    s.rename_axis("year").rename(columns={"north": "north_m", "south": "south_m", "gap": "gap_m"}).to_csv(
        workspace_twin(TABLES) / "photo_gap_annual.csv", float_format="%.2f")
    eras = []
    for lo, hi in ERAS:
        r, se, k = era_rate(s["gap"], lo, hi)
        eras.append(dict(era=f"{lo}-{hi}", rate_m_yr=r, se_m_yr=se, n_years=k))
    pd.DataFrame(eras).to_csv(workspace_twin(TABLES) / "photo_era_rates.csv", index=False, float_format="%.3f")
    print(s.round(1).to_string())
    print(pd.DataFrame(eras).round(2).to_string())


if __name__ == "__main__":
    main()
