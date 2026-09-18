r"""
duneline_lrr.py
==============================================================================
A linear regression rate through the digitised DUNE LINES, per transect, in
the layout of the CoastSat LRR product, so the same LOESS target builder and
the same scoring code read it.

WHY  (Hannah, 2026-09-16: "score it using the other dune line geojsons I
      have measured to make LRR rates")
    CASCADE's shoreline is a dune line with a fixed berm in front of it, so
    the dune line is the feature the model represents, and CoastSat's
    waterline is not (the two are close to uncorrelated in 1996-2010 and
    2004-2024, see 5-scr/4-comparisons/duneline_vs_coastsat/README.md). A two-survey
    endpoint is the weakest estimator of a rate; with every island-wide line
    inside a window in the fit, one bad survey has less leverage. This is
    that fit.

WHAT IS FITTED
    Island-wide lines only: 1984, 1997, 2004, 2009, 2023, all built by
    duneline_to_raw_offsets.py on the same 450 transects (LineID 12-463,
    five per 500 m domain). A window takes every one of them from its start
    vintage to its end vintage inclusive (hat_topo_version
    .DUNE_LINE_FOR_YEAR: 1996 starts on the 1997 line, 2010 on the 2009 line,
    2024 ends on the 2023 line), so

        1984-2004   1984, 1997, 2004          three surveys
        1996-2010   1997, 2004, 2009          three
        2004-2024   2004, 2009, 2023          three
        2010-2024   2009, 2023                two: the OLS IS the endpoint

    The Buxton-only clips (1967, 2017; GIS 2-12, ArcGIS exports about a metre
    landward of the shapely build) and the 1978 line (island-wide, but before
    every window) are NOT used, by Hannah's choice, so every domain in a
    window has the same design and every line the same method.

    Per transect: ordinary least squares of seaward position (-ORIG_LEN, the
    station from the offshore datum, sign flipped so seaward is positive)
    against decimal survey date. Survey dates from duneline_vs_coastsat
    .KNOWN_SURVEY_DATES by vintage; the 2023 NOAA set has no known flight
    date and is centred on 2023-07-01, flagged in PROVENANCE.md. With three
    points the slope's standard error is reported; with two it is NaN.

OUTPUT   data/hatteras_init/5-scr/3-rates/duneline_lrr/<start>_<end>/
             transect_lrr_full.csv     transect_id, domain_number, match_method,
                                       lrr_m_yr, r_squared, p_value, unc_m_yr,
                                       n_obs, start_date, end_date
                                       -- the CoastSat columns, so
                                       CoastSatDataset reads it unchanged
             domain_lrr_summary.csv    per-domain mean/median/std, as CoastSat's
             PROVENANCE.md             the surveys, their dates, what is assumed

USAGE
    python scripts/input_prep/5-scr/duneline_lrr/duneline_lrr.py            # all four windows
    python scripts/input_prep/5-scr/duneline_lrr/duneline_lrr.py --start-year 2004 --end-year 2024

Author: Hannah A. Henry, UNC CECL
==============================================================================
"""
from __future__ import annotations

import argparse
import datetime as _dt
import importlib.util
import sys
from pathlib import Path

import numpy as np
import pandas as pd

PROJECT_ROOT = next(_p for _p in Path(__file__).resolve().parents
                    if (_p / "pyproject.toml").exists())
sys.path.insert(0, str(PROJECT_ROOT / "scripts"))

from site_layer.hat_topo_version import dune_line_for_year, dune_raw_file  # noqa: E402
from site_layer.hat_observed_rates import DUNELINE_LRR_ROOT, TRANSECT_FILE, DOMAIN_FILE  # noqa: E402

ISLAND_WIDE = (1984, 1997, 2004, 2009, 2023)
WINDOWS = [(1984, 2004), (1996, 2010), (2004, 2024), (2010, 2024)]
GIS_FIRST, GIS_LAST = 1, 90


def _dune_module():
    path = (PROJECT_ROOT / "scripts" / "input_prep" / "5-scr"
            / "duneline_vs_coastsat" / "duneline_vs_coastsat.py")
    spec = importlib.util.spec_from_file_location("duneline_vs_coastsat", path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def survey_date(dates, vintage):
    known = dates.get(vintage)
    if known:
        return _dt.datetime.strptime(known, "%Y-%m-%d"), False
    return _dt.datetime(vintage, 7, 1), True


def decimal_year(d):
    y0 = _dt.datetime(d.year, 1, 1)
    y1 = _dt.datetime(d.year + 1, 1, 1)
    return d.year + (d - y0).total_seconds() / (y1 - y0).total_seconds()


def stations(vintage):
    """One station per (domain, transect) for one island-wide line."""
    raw = pd.read_csv(dune_raw_file(vintage), encoding="utf-8-sig")
    t = raw.drop_duplicates(subset=["domain_id", "LineID"])
    t = t[t["domain_id"].between(GIS_FIRST, GIS_LAST)]
    return t[["domain_id", "LineID", "ORIG_LEN"]].rename(
        columns={"ORIG_LEN": f"p{vintage}"})


def fit_window(start, end, dates):
    v0, v1 = dune_line_for_year(start), dune_line_for_year(end)
    vintages = [v for v in ISLAND_WIDE if v0 <= v <= v1]
    assert vintages[0] == v0 and vintages[-1] == v1, (vintages, v0, v1)
    meta = []
    table = None
    for v in vintages:
        d, assumed = survey_date(dates, v)
        meta.append(dict(vintage=v, date=d.date().isoformat(), assumed=assumed,
                         t=decimal_year(d)))
        s = stations(v)
        table = s if table is None else table.merge(s, on=["domain_id", "LineID"], how="inner")
    t = np.array([m["t"] for m in meta])
    # seaward positive: ORIG_LEN grows landward
    y = -table[[f"p{m['vintage']}" for m in meta]].to_numpy(dtype=float)
    n = len(t)
    tm = t - t.mean()
    sxx = float((tm ** 2).sum())
    slope = (y * tm).sum(axis=1) / sxx
    intercept = y.mean(axis=1) - slope * t.mean()
    yhat = intercept[:, None] + slope[:, None] * t[None, :]
    ss_res = ((y - yhat) ** 2).sum(axis=1)
    ss_tot = ((y - y.mean(axis=1, keepdims=True)) ** 2).sum(axis=1)
    with np.errstate(divide="ignore", invalid="ignore"):
        r2 = np.where(ss_tot > 0, 1.0 - ss_res / ss_tot, np.nan)
        if n > 2:
            se = np.sqrt(ss_res / (n - 2) / sxx)
            from scipy import stats
            tstat = slope / se
            p = 2 * stats.t.sf(np.abs(tstat), df=n - 2)
        else:
            se = np.full(len(slope), np.nan)
            p = np.full(len(slope), np.nan)
            r2 = np.full(len(slope), np.nan)   # two points: no residual
    out = pd.DataFrame({
        "transect_id": ["dune_{0}".format(int(i)) for i in table["LineID"]],
        "domain_number": table["domain_id"].astype(float).to_numpy(),
        "match_method": "duneline_transect",
        "lrr_m_yr": np.round(slope, 4),
        "r_squared": np.round(r2, 4),
        "p_value": np.round(p, 4),
        "unc_m_yr": np.round(se, 4),
        "n_obs": n,
        "start_date": meta[0]["date"],
        "end_date": meta[-1]["date"],
    })
    return out, meta


def domain_summary(full):
    g = full.groupby("domain_number")["lrr_m_yr"]
    s = pd.DataFrame({
        "n_valid": g.count(),
        "mean_lrr": g.mean().round(3),
        "median_lrr": g.median().round(3),
        "std_lrr": g.std().round(3),
        "min_lrr": g.min().round(3),
        "max_lrr": g.max().round(3),
        "pct_eroding": (g.apply(lambda v: 100.0 * (v < 0).mean())).round(1),
        "n_transects": g.size(),
    }).reset_index()
    return s


def provenance(start, end, meta, full):
    lines = [f"# duneline_lrr/{start}_{end} - provenance", "",
             f"Written {_dt.date.today().isoformat()} by "
             "scripts/input_prep/5-scr/duneline_lrr/duneline_lrr.py.", "",
             "## Surveys in the fit", "",
             "| vintage | file | survey date | |", "|---|---|---|---|"]
    for m in meta:
        lines.append(f"| {m['vintage']} | {dune_raw_file(m['vintage']).name} | {m['date']} | "
                     f"{'ASSUMED mid-year, flight date unknown' if m['assumed'] else ''} |")
    n = len(meta)
    lines += ["",
              f"{n} island-wide surveys on {len(full)} transects; ordinary least squares of "
              "seaward position (-ORIG_LEN) against decimal survey date per transect. "
              + ("With two surveys the slope is the endpoint rate and r_squared, p_value "
                 "and unc_m_yr are NaN." if n == 2 else
                 f"unc_m_yr is the slope's standard error with {n - 2} degree(s) of freedom."),
              "",
              "Buxton-only clips (1967, 2017) and the 1978 line are not in the fit; see the "
              "script header.", "",
              "## Island summary", "",
              f"mean of transect rates {full['lrr_m_yr'].mean():+.3f} m/yr, "
              f"median {full['lrr_m_yr'].median():+.3f}, "
              f"{100 * (full['lrr_m_yr'] < 0).mean():.0f}% of transects landward.", ""]
    return "\n".join(lines)


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--start-year", type=int)
    ap.add_argument("--end-year", type=int)
    args = ap.parse_args(argv)
    windows = ([(args.start_year, args.end_year)]
               if args.start_year and args.end_year else WINDOWS)
    dates = _dune_module().KNOWN_SURVEY_DATES
    for start, end in windows:
        full, meta = fit_window(start, end, dates)
        out = DUNELINE_LRR_ROOT / f"{start}_{end}"
        out.mkdir(parents=True, exist_ok=True)
        full.to_csv(out / TRANSECT_FILE, index=False)
        domain_summary(full).to_csv(out / DOMAIN_FILE, index=False)
        (out / "PROVENANCE.md").write_text(provenance(start, end, meta, full),
                                          encoding="utf-8")
        surveys = ", ".join(f"{m['vintage']} ({m['date']}{'*' if m['assumed'] else ''})"
                            for m in meta)
        print(f"{start}-{end}  {len(meta)} surveys: {surveys}  "
              f"{len(full)} transects  mean {full['lrr_m_yr'].mean():+.3f} m/yr  "
              f"-> {out.relative_to(PROJECT_ROOT)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
