"""
Step 3: the mean shoreline over a window, as one line (geojson) you can measure like a dune line.

    python mean_shoreline_template.py --start-year 1995 --end-year 1997
    python mean_shoreline_template.py --from 1995-10-12 --to 1997-10-12

Per transect: mean chainage inside the window, placed on the map as
origin + mean chainage * unit vector along the transect, then the points
joined in alongshore (id) order. Window as shoreline_rates_template.py.
Input: the CoastSat transect layer, and one CSV per transect (col 1 date,
col 2 chainage in m from the transect origin), filename = transect id.
Does not need the step-1 lookup. Needs geopandas, pandas, numpy, matplotlib.

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-10-28
"""

from __future__ import annotations

import argparse
from pathlib import Path

import geopandas as gpd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from shapely.geometry import LineString

# --- CONFIG ------------------------------------------------------------------
TRANSECTS_FILE  = Path("data/transects.geojson")
TRANSECT_ID_COL = "id"                 # CoastSat id = time-series filename
TIMESERIES_DIR  = Path("data/timeseries")
OUTPUT_DIR      = Path("output")
START_YEAR, END_YEAR = 1995, 1997
TARGET_CRS      = "EPSG:26918"         # metres; use the CRS of the line you will compare against

MIN_OBS = 10
# -----------------------------------------------------------------------------


# One spelling for an id: CoastSat layers write "usa_NC_0032-0001" or "usa_NC_0032_0001"
def norm_id(tid) -> str:
    return str(tid).strip().replace("-", "_")


# Read one transect CSV: date + chainage, by column order
def load_timeseries(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path)
    df = df.rename(columns={df.columns[0]: "date", df.columns[1]: "chainage_m"})
    df["date"] = pd.to_datetime(df["date"], utc=True, errors="coerce")
    df["chainage_m"] = pd.to_numeric(df["chainage_m"], errors="coerce")
    return df.dropna(subset=["date", "chainage_m"]).sort_values("date")


# Origin and seaward unit vector of every transect, in TARGET_CRS
def load_transects(path: Path) -> tuple[gpd.GeoDataFrame, dict]:
    gdf = gpd.read_file(path)
    if gdf.crs is None:
        raise SystemExit(f"{path} has no CRS; set it at the source, don't guess")
    if TRANSECT_ID_COL not in gdf.columns:
        raise SystemExit(f"{path} has no column '{TRANSECT_ID_COL}': {list(gdf.columns)}")
    gdf = gdf.to_crs(TARGET_CRS)
    print(f"  {len(gdf):,} transects from {path.name}, reprojected to {TARGET_CRS}")
    geom = {}
    for tid, line in zip(gdf[TRANSECT_ID_COL], gdf.geometry):
        if line.geom_type == "MultiLineString":
            line = max(line.geoms, key=lambda g: g.length)
        (x0, y0), (x1, y1) = line.coords[0], line.coords[-1]   # first vertex = origin
        length = np.hypot(x1 - x0, y1 - y0)
        if length > 0:
            geom[norm_id(tid)] = (x0, y0, (x1 - x0) / length, (y1 - y0) / length)
    return gdf, geom


# Window from the command line: whole calendar years, or first and last day.
# Returns (lo, hi, tag), lo inclusive, hi exclusive (the day after the last day)
def choose_window(a: argparse.Namespace):
    if bool(a.date_from) != bool(a.date_to):
        raise SystemExit("give both --from and --to, or neither")
    if a.date_from:
        if a.start_year is not None or a.end_year is not None:
            raise SystemExit("give --from/--to or --start-year/--end-year, not both")
        lo = pd.Timestamp(a.date_from, tz="UTC")
        last = pd.Timestamp(a.date_to, tz="UTC")
        if last < lo:
            raise SystemExit("--to is before --from")
        return lo, last + pd.Timedelta(days=1), f"{lo.date()}_{last.date()}"
    y0 = START_YEAR if a.start_year is None else a.start_year
    y1 = END_YEAR if a.end_year is None else a.end_year
    if y1 < y0:
        raise SystemExit("--end-year is before --start-year")
    return (pd.Timestamp(f"{y0}-01-01", tz="UTC"),
            pd.Timestamp(f"{y1 + 1}-01-01", tz="UTC"), f"{y0}_{y1}")


# Mean chainage in the window, and where that puts the shoreline
def mean_position(df: pd.DataFrame, geom) -> dict:
    n = len(df)
    if n == 0:
        return {"n_obs": 0}
    ch = df["chainage_m"].to_numpy(dtype=float)
    x0, y0, ux, uy = geom
    m = ch.mean()
    sd = ch.std(ddof=1) if n > 1 else np.nan
    return {"n_obs": n, "mean_chainage_m": m, "sd_chainage_m": sd,
            "se_chainage_m": sd / np.sqrt(n) if n > 1 else np.nan,
            "first_date": df["date"].min().date().isoformat(),
            "last_date": df["date"].max().date().isoformat(),
            "x": x0 + m * ux, "y": y0 + m * uy}


# Map of the transects and the line, to check the alongshore order
def draw_map(transects, points, tag, path):
    fig, ax = plt.subplots(figsize=(9, 9))
    transects.plot(ax=ax, color="0.7", lw=0.6)
    ax.plot(points["x"], points["y"], "-", color="#2166ac", lw=1.2, zorder=3,
            label=f"mean shoreline ({len(points)} transects)")
    ax.plot(points["x"].iloc[[0, -1]], points["y"].iloc[[0, -1]], "o",
            color="#b2182b", ms=4, zorder=4, label="first and last vertex")
    ax.set_title(f"Mean shoreline, {tag.replace('_', ' to ')}")
    ax.set_aspect("equal")
    ax.legend(fontsize=8)
    ax.set_axis_off()
    fig.tight_layout()
    fig.savefig(path, dpi=200)
    plt.close(fig)


# Run: window, average, place on the map, join, write
def main() -> None:
    ap = argparse.ArgumentParser(description="Mean shoreline over a window, as a line.")
    ap.add_argument("--start-year", type=int)
    ap.add_argument("--end-year", type=int)
    ap.add_argument("--from", dest="date_from", metavar="YYYY-MM-DD")
    ap.add_argument("--to", dest="date_to", metavar="YYYY-MM-DD")
    ap.add_argument("--transects", type=Path, default=TRANSECTS_FILE)
    ap.add_argument("--timeseries", type=Path, default=TIMESERIES_DIR)
    ap.add_argument("--out", type=Path, default=OUTPUT_DIR)
    a = ap.parse_args()
    lo, hi, tag = choose_window(a)

    # Find the transect files and the geometry
    files = sorted(a.timeseries.rglob("*.csv"))
    if not files:
        raise SystemExit(f"no CSVs in {a.timeseries.resolve()}")
    transects, geom = load_transects(a.transects)

    # Mean position of every transect inside the window
    rows = []
    for f in files:
        tid = norm_id(f.stem)
        if tid not in geom:
            rows.append({"transect_id": tid, "n_obs": 0, "included": False,
                         "excluded_because": "not in the transect layer"})
            continue
        df = load_timeseries(f)
        df = df[(df["date"] >= lo) & (df["date"] < hi)]
        r = {"transect_id": tid, **mean_position(df, geom[tid])}
        r["included"] = r["n_obs"] >= MIN_OBS
        r["excluded_because"] = "" if r["included"] else f"n_obs {r['n_obs']} < {MIN_OBS}"
        rows.append(r)
    # Alongshore order = id order (CoastSat ids are zero-padded site_number); check the map
    per_transect = pd.DataFrame(rows).sort_values("transect_id").reset_index(drop=True)
    points = per_transect[per_transect["included"]].reset_index(drop=True)
    print(f"{len(files)} transects, {lo.date()} to {(hi - pd.Timedelta(days=1)).date()}: "
          f"{len(points)} used, {len(per_transect) - len(points)} excluded")
    if len(points) < 2:
        raise SystemExit("fewer than two usable transects; nothing to join into a line")

    # Write the per-transect table and the line
    a.out.mkdir(parents=True, exist_ok=True)
    per_transect.to_csv(a.out / f"mean_shoreline_points_{tag}.csv", index=False)
    line = gpd.GeoDataFrame(
        {"feature_type": ["Shoreline (CoastSat window mean)"],
         "window": [f"{lo.date()}/{(hi - pd.Timedelta(days=1)).date()}"],
         "n_transects": [len(points)],
         "median_n_obs": [int(points["n_obs"].median())]},
        geometry=[LineString(zip(points["x"], points["y"]))], crs=TARGET_CRS)
    line.to_file(a.out / f"mean_shoreline_{tag}.geojson", driver="GeoJSON")

    # Draw the check map
    draw_map(transects, points, tag, a.out / f"mean_shoreline_{tag}.png")

    print(f"  positions per transect: median {int(points['n_obs'].median())}, "
          f"min {int(points['n_obs'].min())}; standard error of a mean: "
          f"median {points['se_chainage_m'].median():.1f} m")
    print(f"  wrote line, table and map to {a.out.resolve()} -- check the map")


if __name__ == "__main__":
    main()
