#!/usr/bin/env python3
"""
The CoastSat gap across the Buxton groin, 1984-2025, and the year it stopped widening.

    python coastsat_gap_across_groins.py   ->  tables/*.csv  (figures: coastsat_gap_figures.py)

Annual means per transect, anomaly about each transect's own mean over the fit
window, then two gaps, seaward positive (updrift minus downdrift, so a rising gap
is the groin holding the north side): GIS 6 minus GIS 5 domain means (the fit
quantity), and transects within NEAR_M north and south of the groin field (the
condition check). A continuous one-break hinge is fitted on the pre-fill years;
the break year's interval comes from a residual bootstrap. A two-break hinge
and a level step at the model's 2004 failure are compared by BIC. CoastSat only:
the aerial photos are in ../1-aerial_photo_gap/, the two compared in ../3-photos_vs_coastsat/.

For the groin-placement question (GIS 5/6 or 6/7) it also writes the gap across both
boundaries from the GIS 5, 6 and 7 side means, and every GIS 4-8 transect's rate in each
era, located by northing from the southern groin.

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-10-14
"""
from __future__ import annotations

import json
from itertools import combinations
from pathlib import Path

import numpy as np
import pandas as pd
from pyproj import Transformer

HERE = Path(__file__).resolve().parent
REPO = next(p for p in HERE.parents if (p / "pyproject.toml").exists())

# --- CONFIG ------------------------------------------------------------------
FRAME = REPO / "data" / "hatteras_init" / "5-scr" / "2-transect-frame" / "transect_domains"
SERIES_DIR = REPO / "data" / "hatteras_init" / "5-scr" / "1-observations" / "coastsat_timeseries"
GROINS = REPO / "hard-structures" / "groin" / "1-observations" / "data" / "gis" / "groins_hatteras.geojson"
TABLES = HERE / "tables"
SITE = "usa_NC_0032"                 # the CoastSat site holding Buxton
UPDRIFT, DOWNDRIFT = 6, 5            # GIS domains north and south of the groin
NEAR_M = 500.0                       # near-field band beyond each end of the groin field
MIN_OBS_YEAR = 3                     # observations a transect needs in a year
MIN_SHARE = 0.5                      # share of a side's transects a year needs
FIT = (1984, 2016)                   # pre-fill: the Buxton 2017 fill lands on GIS 6
BREAK_RANGE = (1988, 2012)           # candidate break years, >= 4 yr from each end
MIN_SEG = 4                          # shortest segment in the two-break search, years
N_BOOT = 2000
N_BOOT_CHECK = 500                   # per robustness check
SEED = 20261008
STEP_YEAR = 2004                     # the model's failure step, after Isabel
ERAS = ((1984, 1995), (1996, 2003), (2004, 2016))   # to the last repair, after it, after Isabel
BOUNDARIES = ((5, 6), (6, 7))        # the two candidate places for the model groin, (south, north) domain
PROFILE_DOMAINS = (4, 5, 6, 7, 8)    # the alongshore profile
PROFILE_ERAS = ERAS + ((2017, 2025),)                # the three eras, and after the fill
# Robustness: (label, fit window, years left out)
CHECKS = (("as fitted", FIT, ()), ("from 1988", (1988, FIT[1]), ()),
          ("1995 left out", FIT, (1995,)), ("to 2013", (FIT[0], 2013), ()))
# -----------------------------------------------------------------------------


# Northing of each end of the groin field, from the four groin lines
def groin_field():
    g = json.load(open(GROINS))
    ys = [c[1] for f in g["features"] for c in f["geometry"]["coordinates"]]
    return min(ys), max(ys)


# Both ends (UTM 18N) of every transect at the site
def transect_lines():
    g = json.load(open(FRAME / "CoastSat_transect_layer.geojson"))
    to_utm = Transformer.from_crs("EPSG:4326", "EPSG:26918", always_xy=True)
    rows = {}
    for f in g["features"]:
        tid = f["properties"]["id"].replace("-", "_")
        if tid.startswith(SITE):
            (x0, y0), (x1, y1) = (to_utm.transform(*c) for c in f["geometry"]["coordinates"][:2])
            rows[tid] = dict(x0=x0, y0=y0, x1=x1, y1=y1)
    return pd.DataFrame(rows).T


# The transects on each side, for both gap definitions
def sides():
    lut = pd.read_csv(FRAME / "transect_domain_lookup.csv")
    dom = {d: lut.loc[lut["domain_number"] == d, "transect_id"].tolist() for d in (UPDRIFT, DOWNDRIFT)}
    south, north = groin_field()
    y = transect_lines()["y0"]
    near_up = y[(y > north) & (y <= north + NEAR_M)].index.tolist()
    near_down = y[(y < south) & (y >= south - NEAR_M)].index.tolist()
    return {"domain": (dom[UPDRIFT], dom[DOWNDRIFT]), "near_field": (near_up, near_down)}


# One transect's image positions: date and chainage, one row per image
def images(tid):
    d = pd.read_csv(SERIES_DIR / f"{SITE}_timeseries" / f"{tid}.csv", header=0)
    d.columns = ["date", "chainage_m"] + list(d.columns[2:])
    d["date"] = pd.to_datetime(d["date"], utc=True)
    d["year"] = d["date"].dt.year
    d["chainage_m"] = pd.to_numeric(d["chainage_m"], errors="coerce")
    return d.dropna(subset=["chainage_m"])[["date", "year", "chainage_m"]]


# One transect's calendar-year means and image counts, every year
def annual_counts(tid):
    g = images(tid).groupby("year")["chainage_m"]
    return pd.DataFrame({"mean": g.mean(), "n": g.size()})


# One transect's calendar-year means, years with at least MIN_OBS_YEAR images
def annual(tid):
    out = annual_counts(tid)
    return out[out["n"] >= MIN_OBS_YEAR]["mean"]


# Side mean of per-transect anomalies (each about its own fit-window mean), years with enough transects
def side_series(tids):
    a = pd.DataFrame({t: annual(t) for t in tids})
    a = a - a.loc[FIT[0]:FIT[1]].mean()
    ok = a.notna().sum(axis=1) >= np.ceil(MIN_SHARE * len(tids))
    return a[ok].mean(axis=1), a.notna().sum(axis=1)[ok]


# Design matrix: continuous hinge with breaks at `brk`, plus level steps at `steps`
def design(t, brk, steps=()):
    return np.column_stack([np.ones_like(t), t] + [np.clip(t - b, 0, None) for b in brk]
                           + [(t >= s).astype(float) for s in steps])


# Least squares on that design
def hinge(t, y, brk, steps=()):
    X = design(t, brk, steps)
    beta, *_ = np.linalg.lstsq(X, y, rcond=None)
    fit = X @ beta
    return beta, fit, float(((y - fit) ** 2).sum())


# Profile of SSE over one break year, and the best break
def one_break(t, y):
    cand = np.arange(BREAK_RANGE[0], BREAK_RANGE[1] + 1, dtype=float)
    sse = np.array([hinge(t, y, [b])[2] for b in cand])
    return cand, sse, cand[np.argmin(sse)]


# Best pair of breaks
def two_breaks(t, y):
    cand = np.arange(t.min() + MIN_SEG, t.max() - MIN_SEG + 1, dtype=float)
    best = min(((hinge(t, y, [a, b])[2], a, b) for a, b in combinations(cand, 2) if b - a >= MIN_SEG))
    return best[1], best[2], best[0]


# BIC of a least-squares fit with k parameters (a break year counts as one)
def bic(sse, n, k):
    return n * np.log(sse / n) + k * np.log(n)


# Residual bootstrap of the break year
def boot_break(t, y, brk, n_boot):
    _, fit, _ = hinge(t, y, [brk])
    res = y - fit
    rng = np.random.default_rng(SEED)
    return np.array([one_break(t, fit + rng.choice(res, res.size, replace=True))[2]
                     for _ in range(n_boot)])


# Does a level drop at STEP_YEAR earn its place beside the free break? Size, SE and BIC
def step_test(t, y):
    n = len(t)
    cand = np.arange(BREAK_RANGE[0], BREAK_RANGE[1] + 1, dtype=float)
    sse, b = min((hinge(t, y, [b], [STEP_YEAR])[2], b) for b in cand if b != STEP_YEAR)
    X = design(t, [b], [STEP_YEAR])
    beta, *_ = np.linalg.lstsq(X, y, rcond=None)
    cov = sse / (n - X.shape[1]) * np.linalg.inv(X.T @ X)
    return dict(step_break_year=int(b), step_m=float(beta[-1]), step_se_m=float(np.sqrt(cov[-1, -1])),
                bic_line_plus_step=bic(hinge(t, y, [], [STEP_YEAR])[2], n, 3),
                bic_break_plus_step=bic(sse, n, 5))


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


# The gap across each candidate boundary: north side mean minus south side mean, each about its own mean
def boundary_gaps():
    lut = pd.read_csv(FRAME / "transect_domain_lookup.csv")
    side = {d: side_series(lut.loc[lut["domain_number"] == d, "transect_id"].tolist())[0]
            for d in sorted({d for b in BOUNDARIES for d in b})}
    out = pd.DataFrame({f"gis{d}_anom_m": s for d, s in side.items()})
    for s, n in BOUNDARIES:
        out[f"gap_{s}_{n}_m"] = out[f"gis{n}_anom_m"] - out[f"gis{s}_anom_m"]
    return out.rename_axis("year")


# Every GIS 4-8 transect's least-squares rate in each era, and where it sits relative to the groins
def alongshore_rates():
    lut = pd.read_csv(FRAME / "transect_domain_lookup.csv")
    lut = lut[lut["domain_number"].isin(PROFILE_DOMAINS)]
    ends = transect_lines()
    south, north = groin_field()
    rows = []
    for _, r in lut.iterrows():
        a = annual(r["transect_id"])
        for lo, hi in PROFILE_ERAS:
            rate, se, k = era_rate(a, lo, hi)
            rows.append(dict(transect_id=r["transect_id"], domain=int(r["domain_number"]),
                             northing=ends.loc[r["transect_id"], "y0"],
                             from_south_groin_m=ends.loc[r["transect_id"], "y0"] - south,
                             era=f"{lo}-{hi}", rate_m_yr=rate, se_m_yr=se, n_years=k))
    return pd.DataFrame(rows).sort_values(["era", "northing"])


def main():
    TABLES.mkdir(exist_ok=True)
    boundary_gaps().to_csv(TABLES / "boundary_gaps.csv", float_format="%.2f")
    alongshore_rates().to_csv(TABLES / "alongshore_rates.csv", index=False, float_format="%.3f")
    groups = sides()
    lines = transect_lines()
    used = [dict(transect_id=t, gap=name, side=side, **lines.loc[t].to_dict())
            for name, (up, down) in groups.items()
            for side, tids in (("updrift", up), ("downdrift", down)) for t in tids]
    pd.DataFrame(used).round(2).to_csv(TABLES / "transects_used.csv", index=False)

    # Each transect's annual mean, image count and anomaly about its own fit-window mean (figures 2-3)
    per = []
    for name, (up, down) in groups.items():
        for side, tids in (("updrift", up), ("downdrift", down)):
            for t in tids:
                a = annual_counts(t)
                ref = a.loc[a["n"] >= MIN_OBS_YEAR, "mean"].loc[FIT[0]:FIT[1]].mean()
                a = a.assign(transect_id=t, gap=name, side=side, used=a["n"] >= MIN_OBS_YEAR,
                             anom_m=a["mean"] - ref, fit_mean_m=ref)
                per.append(a.rename_axis("year").reset_index())
    pd.concat(per)[["gap", "side", "transect_id", "year", "n", "used", "mean", "fit_mean_m", "anom_m"]].rename(
        columns={"n": "n_images", "mean": "annual_mean_m"}).to_csv(
        TABLES / "transect_annual.csv", index=False, float_format="%.2f")

    series, rows, fits, checks, profiles, boots_rows, hinge_rows, eras = {}, [], [], [], [], [], [], []
    for name, (up, down) in groups.items():
        su, nu = side_series(up)
        sd, nd = side_series(down)
        gap = (su - sd).dropna()
        series[name] = gap
        rows += [dict(gap=name, year=int(yr), gap_m=gap[yr], updrift_anom_m=su[yr],
                      downdrift_anom_m=sd[yr], n_updrift=int(nu[yr]), n_downdrift=int(nd[yr]))
                 for yr in gap.index]
        w = gap.loc[FIT[0]:FIT[1]]
        t, y = w.index.to_numpy(float), w.to_numpy(float)
        cand, sse, brk = one_break(t, y)
        boots = boot_break(t, y, brk, N_BOOT)
        beta, fit, sse1 = hinge(t, y, [brk])
        b1, b2, sse2 = two_breaks(t, y)
        beta2, _, _ = hinge(t, y, [b1, b2])
        n = len(t)
        fits.append(dict(
            gap=name, n_updrift=len(up), n_downdrift=len(down), n_years=n,
            break_year=int(brk), break_ci90_lo=int(np.percentile(boots, 5)),
            break_ci90_hi=int(np.percentile(boots, 95)),
            boot_share_2002_2005=float(np.mean((boots >= 2002) & (boots <= 2005))),
            boot_share_le_1997=float(np.mean(boots <= 1997)),
            rate_before_m_yr=float(beta[1]), rate_after_m_yr=float(beta[1] + beta[2]),
            two_break_years=f"{int(b1)} {int(b2)}",
            two_break_rates_m_yr=" ".join(f"{v:+.2f}" for v in np.cumsum(beta2[1:4])),
            bic_line=bic(hinge(t, y, [])[2], n, 2), bic_one_break=bic(sse1, n, 4),
            bic_two_breaks=bic(sse2, n, 6), **step_test(t, y)))
        profiles += [dict(gap=name, break_year=int(c), sse_over_best=s / sse.min())
                     for c, s in zip(cand, sse)]
        boots_rows += [dict(gap=name, break_year=int(b)) for b in boots]
        hinge_rows += [dict(gap=name, year=int(a), fit_m=f) for a, f in zip(t, fit)]
        for lo, hi in ERAS:
            r, se, k = era_rate(gap, lo, hi)
            eras.append(dict(source=name, era=f"{lo}-{hi}", rate_m_yr=r, se_m_yr=se, n_years=k))
        for label, (lo, hi), drop in CHECKS:
            w = gap.loc[lo:hi].drop(list(drop), errors="ignore")
            tt, yy = w.index.to_numpy(float), w.to_numpy(float)
            bb = one_break(tt, yy)[2]
            bs = boot_break(tt, yy, bb, N_BOOT_CHECK)
            checks.append(dict(gap=name, check=label, n_years=len(tt), break_year=int(bb),
                               bic_one_break=bic(hinge(tt, yy, [bb])[2], len(tt), 4),
                               break_ci90_lo=int(np.percentile(bs, 5)),
                               break_ci90_hi=int(np.percentile(bs, 95)), **step_test(tt, yy)))

    pd.DataFrame(rows).to_csv(TABLES / "gap_annual.csv", index=False, float_format="%.2f")
    pd.DataFrame(fits).to_csv(TABLES / "gap_breakpoint.csv", index=False, float_format="%.3f")
    pd.DataFrame(checks).to_csv(TABLES / "gap_robustness.csv", index=False, float_format="%.3f")
    pd.DataFrame(profiles).to_csv(TABLES / "break_profile.csv", index=False, float_format="%.4f")
    pd.DataFrame(boots_rows).to_csv(TABLES / "break_bootstrap.csv", index=False)
    pd.DataFrame(hinge_rows).to_csv(TABLES / "hinge_fit.csv", index=False, float_format="%.2f")
    pd.DataFrame(eras).to_csv(TABLES / "era_rates.csv", index=False, float_format="%.3f")
    print(pd.DataFrame(fits).round(2).T.to_string())
    print(pd.DataFrame(checks).round(2).to_string())
    print(pd.DataFrame(eras).round(2).to_string())


if __name__ == "__main__":
    main()
