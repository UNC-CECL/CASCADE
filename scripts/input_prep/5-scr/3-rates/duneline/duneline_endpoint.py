"""
Net dune-line change between the two dune lines that bound a window, per transect and per GIS domain.

    python scripts/input_prep/5-scr/3-rates/duneline/duneline_endpoint.py
    python scripts/input_prep/5-scr/3-rates/duneline/duneline_endpoint.py --windows 1996_2010

In metres and m/yr; replaced the dune-line LRR product. Details: scripts/input_prep/5-scr/3-rates/README.md.

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-09-22
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
sys.path.insert(0, str(_REPO / "scripts" / "input_prep" / "5-scr" / "lib"))
import scr_paths  # noqa: E402,F401  (5-scr sibling modules onto sys.path)

from coastsat_vs_duneline import DAYS_PER_YEAR, KNOWN_SURVEY_DATES  # noqa: E402
from site_layer.hat_observed_rates import (  # noqa: E402
    DUNE_ENDPOINT_DOMAIN_FILE, DUNE_ENDPOINT_TRANSECT_FILE, DUNELINE_ENDPOINT_ROOT,
)
from site_layer.hat_topo_version import dune_line_for_year, dune_raw_file  # noqa: E402

# --- CONFIG ------------------------------------------------------------------
# The four model windows and the long context window.
WINDOWS = [(1984, 2004), (1996, 2010), (2004, 2024), (2010, 2024), (1996, 2024)]
ASSUMED_MONTH_DAY = "07-01"
N_DOMAINS = 90
# -----------------------------------------------------------------------------


# (date, assumed) for one line vintage
def survey_date(vintage: int):
    known = KNOWN_SURVEY_DATES.get(vintage)
    if known:
        return dt.date.fromisoformat(known), False
    return dt.date.fromisoformat(f"{vintage}-{ASSUMED_MONTH_DAY}"), True


# First row per transect, as the hindcast loader reads it
def stations(vintage: int) -> pd.DataFrame:
    raw = pd.read_csv(dune_raw_file(vintage), encoding="utf-8-sig")
    return (raw.drop_duplicates(subset=["domain_id", "LineID"])
               [["domain_id", "LineID", "ORIG_LEN"]])


# One window's per-transect and per-domain tables
def build(start: int, end: int) -> dict:
    v0, v1 = dune_line_for_year(start), dune_line_for_year(end)
    d0, a0 = survey_date(v0)
    d1, a1 = survey_date(v1)
    years = (d1 - d0).days / DAYS_PER_YEAR

    t = (stations(v0).rename(columns={"ORIG_LEN": "p0"})
         .merge(stations(v1).rename(columns={"ORIG_LEN": "p1"}),
                on=["domain_id", "LineID"], how="inner"))
    t = t[t["domain_id"].between(1, N_DOMAINS)].sort_values(["domain_id", "LineID"])
    tr = pd.DataFrame({
        "transect_id": "dune_" + t["LineID"].astype(str),
        "domain_number": t["domain_id"].astype(int),
        "line_id": t["LineID"].astype(int),
        f"position_{v0}_m": t["p0"].round(3),
        f"position_{v1}_m": t["p1"].round(3),
        "change_m": (t["p0"] - t["p1"]).round(3),        # seaward positive
        "rate_m_yr": ((t["p0"] - t["p1"]) / years).round(4),
        "start_vintage": v0, "end_vintage": v1,
        "start_date": d0.isoformat(), "end_date": d1.isoformat(),
        "start_date_assumed": a0, "end_date_assumed": a1,
        "interval_yr": round(years, 4),
    })

    g = tr.groupby("domain_number")
    dom = pd.DataFrame({
        "n_transects": g.size(),
        "mean_change_m": g["change_m"].mean(), "std_change_m": g["change_m"].std(),
        "min_change_m": g["change_m"].min(), "max_change_m": g["change_m"].max(),
        "mean_rate_m_yr": g["rate_m_yr"].mean(), "std_rate_m_yr": g["rate_m_yr"].std(),
        "pct_landward": g["change_m"].apply(lambda s: 100.0 * (s < 0).mean()),
    }).round(3).reindex(range(1, N_DOMAINS + 1)).rename_axis("domain_number").reset_index()

    out = DUNELINE_ENDPOINT_ROOT / f"{start}_{end}"
    out.mkdir(parents=True, exist_ok=True)
    tr.to_csv(out / DUNE_ENDPOINT_TRANSECT_FILE, index=False)
    dom.to_csv(out / DUNE_ENDPOINT_DOMAIN_FILE, index=False)

    def row(v, d, a):
        f = dune_raw_file(v)
        built = dt.datetime.fromtimestamp(f.stat().st_mtime).strftime("%Y-%m-%d %H:%M")
        return (f"| {v} | `{f.name}` (written {built}) | {d.isoformat()}"
                f"{' **ASSUMED**' if a else ''} |")

    assumed = [v for v, a in ((v0, a0), (v1, a1)) if a]
    (out / "PROVENANCE.md").write_text("\n".join([
        f"# 3-rates/duneline/endpoint/{start}_{end} - provenance",
        "",
        f"Written {dt.datetime.now():%Y-%m-%d %H:%M} by "
        "scripts/input_prep/5-scr/3-rates/duneline/duneline_endpoint.py.",
        "",
        "Net change between the two dune lines that bound the window, per 100 m "
        "transect, seaward positive; per domain the mean over its transects. "
        "`change_m` needs no date; `rate_m_yr` is `change_m` over the survey "
        "interval.",
        "",
        "| line | raw stations | survey date |",
        "|---|---|---|",
        row(v0, d0, a0),
        row(v1, d1, a1),
        "",
        f"Interval {years:.2f} yr ({start} and {end} read the {v0} and {v1} lines "
        "through DUNE_LINE_FOR_YEAR)."
        + (f" The {', '.join(map(str, assumed))} flight date is not known and is "
           f"centred on 1 July; `rate_m_yr` inherits that assumption, `change_m` "
           f"does not." if assumed else ""),
        "",
        "## Island summary",
        "",
        f"{len(tr)} transects in {int(dom['n_transects'].notna().sum())} domains. "
        f"Domain mean change {dom['mean_change_m'].mean():+.1f} m "
        f"({dom['mean_rate_m_yr'].mean():+.2f} m/yr); "
        f"{int((dom['mean_change_m'] < 0).sum())} of {N_DOMAINS} domains landward; "
        f"range {dom['mean_change_m'].min():+.1f} to {dom['mean_change_m'].max():+.1f} m.",
        "",
    ]), encoding="utf-8")
    return dict(window=f"{start}_{end}", lines=f"{v0}->{v1}", years=years,
                assumed=assumed, n=len(tr),
                mean_m=dom["mean_change_m"].mean(), mean_rate=dom["mean_rate_m_yr"].mean(),
                landward=int((dom["mean_change_m"] < 0).sum()))


# Run: every window, or those given
def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description="Net dune-line change per window.")
    ap.add_argument("--windows", nargs="+", metavar="START_END")
    a = ap.parse_args(argv)
    wins = ([tuple(int(x) for x in w.split("_")) for w in a.windows]
            if a.windows else WINDOWS)
    for s, e in wins:
        r = build(s, e)
        print(f"{r['window']}  lines {r['lines']}  {r['years']:5.2f} yr"
              f"{'  (date assumed: ' + ','.join(map(str, r['assumed'])) + ')' if r['assumed'] else ''}"
              f"  {r['n']} transects  mean {r['mean_m']:+6.1f} m ({r['mean_rate']:+.2f} m/yr)"
              f"  landward {r['landward']}/90")
    return 0


if __name__ == "__main__":
    sys.exit(main())
