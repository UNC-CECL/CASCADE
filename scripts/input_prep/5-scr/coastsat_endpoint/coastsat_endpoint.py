"""
coastsat_endpoint.py
==============================================================================
The stored CoastSat net shoreline change: for each CoastSat transect, the mean
shoreline position in a +/-6-month window about the window's END dune-line
survey date minus the mean about its START survey date, in METRES and m/yr,
per transect and per GIS domain. Built 2026-09-18 (Hannah, by interview) as
the counterpart of 3-rates/duneline/endpoint/, so the shoreline and the dune
line difference like for like: same windows, same survey dates, same sign.

WHY THE DUNE DATES
    A dune line is a survey at a moment. Centring the CoastSat windows on the
    same moments (1997-10-12, 2009-05-30, 2023-07-01 assumed; the 1984 and
    2004 lines for the older windows) means a gap between the two changes is
    beach-width change, not a date mismatch. Each end averages one full
    seasonal cycle of satellite positions. The window means are the ones
    duneline_vs_coastsat.py uses (window_mean / endpoint_by_transect,
    imported, not copied), so its CoastSat endpoint and this agree.

    change_m = end mean - start mean. CoastSat chainage grows SEAWARD, so
    SEAWARD IS POSITIVE, as in the dune product. rate_m_yr divides by the
    survey interval and inherits any assumed date; change_m does not depend on
    the interval (the windows are still centred on the assumed date).

OUTPUT   data/hatteras_init/5-scr/3-rates/coastsat/endpoint/<start>_<end>/
    transect_endpoint.csv         per CoastSat transect: both window means,
                                  how many positions fell in each, the first
                                  and last date inside each, change_m,
                                  rate_m_yr, the survey dates
    domain_endpoint_summary.csv   per domain: n, mean/std/min/max change_m,
                                  mean rate, pct_landward, the median
                                  positions per end window, and how many of
                                  its transects had an empty end window
    PROVENANCE.md
    Read through hat_observed_rates.coastsat_endpoint_csv(start, end, level).

USAGE
    python scripts/input_prep/5-scr/coastsat_endpoint/coastsat_endpoint.py
    python scripts/input_prep/5-scr/coastsat_endpoint/coastsat_endpoint.py --windows 1996_2024
==============================================================================
"""

from __future__ import annotations

import argparse
import datetime as dt
import sys
from pathlib import Path

import numpy as np
import pandas as pd

_REPO = next(_p for _p in Path(__file__).resolve().parents
             if (_p / "pyproject.toml").exists())
sys.path.insert(0, str(_REPO / "scripts"))
for _sub in ("duneline_vs_coastsat", "duneline_endpoint"):
    sys.path.insert(0, str(_REPO / "scripts" / "input_prep" / "5-scr" / _sub))

from duneline_vs_coastsat import (DAYS_PER_YEAR, SIX_MONTHS_DAYS,  # noqa: E402
                                  endpoint_by_transect)
from duneline_endpoint import WINDOWS, survey_date  # noqa: E402
from site_layer.hat_observed_rates import (  # noqa: E402
    COASTSAT_ENDPOINT_ROOT, ENDPOINT_DOMAIN_FILE, ENDPOINT_TRANSECT_FILE,
    transect_lookup,
)
from site_layer.hat_topo_version import dune_line_for_year  # noqa: E402

N_DOMAINS = 90


def build(start: int, end: int, cache: dict) -> dict:
    v0, v1 = dune_line_for_year(start), dune_line_for_year(end)
    d0, a0 = survey_date(v0)
    d1, a1 = survey_date(v1)
    t0 = dt.datetime(d0.year, d0.month, d0.day, tzinfo=dt.timezone.utc)
    t1 = dt.datetime(d1.year, d1.month, d1.day, tzinfo=dt.timezone.utc)
    years = (d1 - d0).days / DAYS_PER_YEAR

    lookup = pd.read_csv(transect_lookup())
    lookup = lookup[lookup["domain_number"].between(1, N_DOMAINS)]
    ep = endpoint_by_transect(lookup, t0, t1, SIX_MONTHS_DAYS, cache)
    for c in ("first_obs_start", "last_obs_start", "first_obs_end", "last_obs_end"):
        ep[c] = pd.to_datetime(ep[c], utc=True).dt.strftime("%Y-%m-%d")
    ep = ep.rename(columns={"mean_start_m": "position_start_m",
                            "mean_end_m": "position_end_m",
                            "endpoint_rate_m_yr": "rate_m_yr"})
    ep = ep.assign(start_vintage=v0, end_vintage=v1,
                   start_date=d0.isoformat(), end_date=d1.isoformat(),
                   start_date_assumed=a0, end_date_assumed=a1,
                   interval_yr=round(years, 4),
                   half_window_days=SIX_MONTHS_DAYS)

    ok = ep[ep["change_m"].notna()]
    g = ok.groupby("domain_number")
    empty = ep.assign(e=(ep["n_start"] == 0) | (ep["n_end"] == 0)).groupby("domain_number")["e"].sum()
    dom = pd.DataFrame({
        "n_transects": g.size(),
        "mean_change_m": g["change_m"].mean(), "std_change_m": g["change_m"].std(),
        "min_change_m": g["change_m"].min(), "max_change_m": g["change_m"].max(),
        "mean_rate_m_yr": g["rate_m_yr"].mean(), "std_rate_m_yr": g["rate_m_yr"].std(),
        "pct_landward": g["change_m"].apply(lambda s: 100.0 * (s < 0).mean()),
        "median_n_start": g["n_start"].median(), "median_n_end": g["n_end"].median(),
        "n_transects_empty_window": empty,
    }).round(3).reindex(range(1, N_DOMAINS + 1)).rename_axis("domain_number").reset_index()

    out = COASTSAT_ENDPOINT_ROOT / f"{start}_{end}"
    out.mkdir(parents=True, exist_ok=True)
    ep.to_csv(out / ENDPOINT_TRANSECT_FILE, index=False, float_format="%.4f")
    dom.to_csv(out / ENDPOINT_DOMAIN_FILE, index=False)

    n_empty = int(((ep["n_start"] == 0) | (ep["n_end"] == 0)).sum())
    assumed = [v for v, a in ((v0, a0), (v1, a1)) if a]
    (out / "PROVENANCE.md").write_text("\n".join([
        f"# 3-rates/coastsat/endpoint/{start}_{end} - provenance",
        "",
        f"Written {dt.datetime.now():%Y-%m-%d %H:%M} by "
        "scripts/input_prep/5-scr/coastsat_endpoint/coastsat_endpoint.py.",
        "",
        "Per CoastSat transect: the mean shoreline position within "
        f"±{SIX_MONTHS_DAYS:g} days of each survey date, end minus start, seaward "
        "positive; per domain the mean over its transects. The survey dates are "
        f"the dune lines' ({start} and {end} read the {v0} and {v1} lines "
        "through DUNE_LINE_FOR_YEAR), so this differences like for like with "
        "3-rates/duneline/endpoint/.",
        "",
        "| end | survey date | median positions per transect | median span inside the window |",
        "|---|---|---|---|",
        f"| start ({v0}) | {d0.isoformat()}{' **ASSUMED**' if a0 else ''} | "
        f"{ep['n_start'].median():.0f} | {pd.to_datetime(ep['first_obs_start']).median().date()} "
        f"to {pd.to_datetime(ep['last_obs_start']).median().date()} |",
        f"| end ({v1}) | {d1.isoformat()}{' **ASSUMED**' if a1 else ''} | "
        f"{ep['n_end'].median():.0f} | {pd.to_datetime(ep['first_obs_end']).median().date()} "
        f"to {pd.to_datetime(ep['last_obs_end']).median().date()} |",
        "",
        f"Interval {years:.2f} yr."
        + (f" The {', '.join(map(str, assumed))} flight date is not known and is "
           "assumed to be 1 July; the end window is centred on it." if assumed else ""),
        f" {len(ep)} transects, {n_empty} with an empty end window (no change "
        "computed).",
        "",
        "## Island summary",
        "",
        f"Domain mean change {dom['mean_change_m'].mean():+.1f} m "
        f"({dom['mean_rate_m_yr'].mean():+.2f} m/yr); "
        f"{int((dom['mean_change_m'] < 0).sum())} of {N_DOMAINS} domains landward; "
        f"range {dom['mean_change_m'].min():+.1f} to {dom['mean_change_m'].max():+.1f} m.",
        "",
    ]), encoding="utf-8")
    return dict(window=f"{start}_{end}", lines=f"{v0}->{v1}", years=years,
                n=len(ep), empty=n_empty, mean_m=dom["mean_change_m"].mean(),
                mean_rate=dom["mean_rate_m_yr"].mean(),
                landward=int((dom["mean_change_m"] < 0).sum()),
                n0=ep["n_start"].median(), n1=ep["n_end"].median())


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description="Net CoastSat shoreline change per window.")
    ap.add_argument("--windows", nargs="+", metavar="START_END")
    a = ap.parse_args(argv)
    wins = ([tuple(int(x) for x in w.split("_")) for w in a.windows]
            if a.windows else WINDOWS)
    cache: dict = {}
    for s, e in wins:
        r = build(s, e, cache)
        print(f"{r['window']}  dune dates {r['lines']}  {r['years']:5.2f} yr  "
              f"{r['n']} transects ({r['empty']} empty)  positions/window "
              f"{r['n0']:.0f}/{r['n1']:.0f}  mean {r['mean_m']:+6.1f} m "
              f"({r['mean_rate']:+.2f} m/yr)  landward {r['landward']}/90")
    return 0


if __name__ == "__main__":
    sys.exit(main())
