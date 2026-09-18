"""
CoastSat rates for the coast BEYOND the 90 surveyed domains
============================================================
The Pea Island extension experiment (2026-09-16) models GIS 1-115 (and
0-115) instead of 1-90, and solves the end-domain source/sink at the new
ends against CoastSat, as the matrix does at GIS 1 and 90. The committed
rate tables stop at GIS 90 only because the transect-to-domain lookup was a
polygon join onto 90 polygons; the CoastSat record itself runs to Oregon
Inlet (263 transects on disk between GIS 90 and 115, ~10 per domain, a
median 260 observations each over 1996-2010).

This script numbers those transects the way the surveyed ones were
numbered -- the transect's origin point within a 500 m domain polygon, now
Hannah's whole-island polygons (hat_extension_domains.join_origins) -- and
fits them with the same LRR as coastsat_domain_lrr_fixed.py, for one window.
A transect no polygon covers is left out, as the surveyed mapping leaves
them out. It writes, beside the surveyed products and never into them:

    5-scr/2-transect-frame/transect_domains/transect_domain_lookup_ext.csv
    5-scr/3-rates/coastsat_lrr/<start>_<end>/ext/transect_lrr_full.csv
    5-scr/3-rates/coastsat_lrr/<start>_<end>/ext/domain_lrr_summary.csv
    5-scr/3-rates/coastsat_lrr/<start>_<end>/ext/transect_lrr_with_base.csv

The last is the surveyed table with the extension rows appended: what an
extended-geometry run loads as its active dataset. The window's own
transect_lrr_full.csv must already exist (coastsat_domain_lrr_fixed.py).

Usage
-----
    python coastsat_extension_lrr.py --start-year 1996 --end-year 2010
"""
from __future__ import annotations

import argparse
import glob
import os
import sys
from pathlib import Path

import geopandas as gpd
import numpy as np
import pandas as pd

_REPO = next(_p for _p in Path(__file__).resolve().parents
             if (_p / "pyproject.toml").exists())
sys.path.insert(0, str(_REPO / "scripts"))
sys.path.insert(0, str(Path(__file__).resolve().parent))

from site_layer.hat_extension_domains import DOMAIN_CRS, SURVEYED_GIS, join_origins  # noqa: E402
from site_layer.hat_observed_rates import (COASTSAT_TIMESERIES, EXT_DIR,  # noqa: E402
                                TRANSECT_DOMAINS, WITH_BASE_FILE, lrr_csv,
                                transect_lookup, window_dir)
from coastsat_lrr_analysis import (_empty_lrr, compute_lrr,  # noqa: E402
                                   filter_dates, load_timeseries)

TRANSECT_LAYER = TRANSECT_DOMAINS / "CoastSat_transect_layer.geojson"
MIN_OBS = 3  # as coastsat_domain_lrr_fixed.py


def extension_lookup():
    """Every on-disk CoastSat transect beyond GIS 1-90, numbered by its polygon."""
    on_disk = {os.path.splitext(os.path.basename(f))[0]
               for f in glob.glob(str(COASTSAT_TIMESERIES / "*" / "*.csv"))}
    layer = gpd.read_file(TRANSECT_LAYER)
    layer.columns = [c.split(".")[-1] for c in layer.columns]
    layer = layer[layer["id"].isin(on_disk)].to_crs(DOMAIN_CRS)
    surveyed = set(pd.read_csv(transect_lookup())["transect_id"])
    layer = layer[~layer["id"].isin(surveyed)].reset_index(drop=True)
    # The origin point within a polygon, the rule of
    # coastsat_domain_mapping.py, onto Hannah's whole-island polygons
    # (2026-09-16). A transect no polygon covers is left out.
    rows = []
    for (_, tr), gis in zip(layer.iterrows(), join_origins(layer)):
        if gis is None or np.isnan(gis):
            continue
        rows.append({"transect_id": tr["id"], "domain_number": int(gis),
                     "distance_m": 0.0, "match_method": "polygon_join",
                     "northing_m": round(float(tr.geometry.coords[0][1]), 3)})
    lookup = pd.DataFrame(rows).sort_values(["domain_number", "transect_id"])
    return lookup.reset_index(drop=True)


def fit_window(lookup, start_year, end_year):
    """coastsat_domain_lrr_fixed.compute_all_lrr, for the extension rows.
    That script parses its arguments at import, so the ten lines are
    repeated here rather than imported."""
    csv_map = {os.path.splitext(os.path.basename(f))[0]: f
               for f in glob.glob(str(COASTSAT_TIMESERIES / "*" / "*.csv"))}
    start, end = f"{start_year}-01-01", f"{end_year}-12-31"
    records = []
    for _, row in lookup.iterrows():
        tid = str(row["transect_id"])
        if tid not in csv_map:
            result = _empty_lrr(0)
        else:
            series = filter_dates(load_timeseries(csv_map[tid]), start, end)
            result = (_empty_lrr(len(series)) if len(series) < MIN_OBS
                      else compute_lrr(series))
        records.append({"transect_id": tid, "domain_number": int(row["domain_number"]),
                        "match_method": row["match_method"], **result})
    return pd.DataFrame(records)


def summarise(transects):
    valid = transects[transects["lrr_m_yr"].notna()]
    out = (valid.groupby("domain_number")["lrr_m_yr"]
           .agg(n_valid="count", mean_lrr="mean", median_lrr="median",
                std_lrr="std", min_lrr="min", max_lrr="max").reset_index())
    pct = valid.groupby("domain_number")["lrr_m_yr"].apply(lambda g: (g < 0).mean() * 100)
    out["pct_eroding"] = out["domain_number"].map(pct)
    total = transects.groupby("domain_number")["transect_id"].count()
    out["n_transects"] = out["domain_number"].map(total)
    for col in ("mean_lrr", "median_lrr", "std_lrr", "min_lrr", "max_lrr", "pct_eroding"):
        out[col] = out[col].round(3)
    return out


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--start-year", type=int, required=True)
    ap.add_argument("--end-year", type=int, required=True)
    a = ap.parse_args(argv)

    base_csv = lrr_csv(a.start_year, a.end_year)   # raises if the window is unbuilt
    lookup = extension_lookup()
    lookup_path = TRANSECT_DOMAINS / "transect_domain_lookup_ext.csv"
    lookup.to_csv(lookup_path, index=False)
    per_dom = lookup.groupby("domain_number").size()
    print(f"extension transects: {len(lookup)}  domains "
          f"{lookup.domain_number.min()}..{lookup.domain_number.max()}  "
          f"({int((lookup.domain_number < SURVEYED_GIS[0]).sum())} south, "
          f"{int((lookup.domain_number > SURVEYED_GIS[1]).sum())} north)")
    print(f"  per domain: min {per_dom.min()} median {per_dom.median():.0f} max {per_dom.max()}")
    print(f"  wrote {lookup_path}")

    fitted = fit_window(lookup, a.start_year, a.end_year)
    out_dir = window_dir(a.start_year, a.end_year) / EXT_DIR
    out_dir.mkdir(parents=True, exist_ok=True)
    fitted.to_csv(out_dir / "transect_lrr_full.csv", index=False)
    summary = summarise(fitted)
    summary.to_csv(out_dir / "domain_lrr_summary.csv", index=False)
    base = pd.read_csv(base_csv)
    with_base = pd.concat([base, fitted[[c for c in base.columns if c in fitted.columns]]],
                          ignore_index=True)
    with_base.to_csv(out_dir / WITH_BASE_FILE, index=False)
    n_valid = int(fitted["lrr_m_yr"].notna().sum())
    print(f"\n{a.start_year}-{a.end_year}: {n_valid} of {len(fitted)} extension "
          f"transects fitted; median n_obs {fitted['n_obs'].median():.0f}")
    print(summary[["domain_number", "n_valid", "mean_lrr", "median_lrr"]].to_string(index=False))
    print(f"\nwrote {out_dir / 'transect_lrr_full.csv'}\n      {out_dir / 'domain_lrr_summary.csv'}"
          f"\n      {out_dir / WITH_BASE_FILE}  ({len(base)} surveyed + {len(fitted)} extension rows)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
