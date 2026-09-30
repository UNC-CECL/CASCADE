"""
transect_zone_join_template.py
==============================================================================
Which transect belongs to which zone: CoastSat transect lines + your zone
polygons in, a two-column lookup table out.

This is a TEMPLATE, not a product. It is the join the 5-scr folder runs
(2-transect-frame/coastsat_domain_mapping.py), stripped of everything specific
to Hatteras Island. Its output is exactly the LOOKUP_CSV that
shoreline_rates_template.py and shoreline_endpoint_template.py read, so it is
step 0 of both.

It imports nothing from this repository. It needs geopandas and matplotlib.

THE PROCESS
-----------
    1  LOAD      the transect lines and the zone polygons. Both must carry a
                 CRS; a file without one is an error here, never a guess.
    2  POINT     reduce each transect line to ONE point -- by default its
                 origin (the landward end CoastSat measures chainage from).
    3  JOIN      point-in-polygon. A point inside exactly one zone is matched.
    4  REFUSE    a point inside NO zone is left out. A point inside TWO zones
                 (overlapping polygons) is left out too. Neither is snapped to
                 the nearest zone: see the note below.
    5  WRITE     the lookup, a problems table, and a map to check by eye.

WHY THERE IS NO NEAREST-NEIGHBOUR SNAPPING
------------------------------------------
Snapping an unmatched transect to the closest zone looks like tidying up and
is the one step in this file that makes a WRONG answer invisible. On a curved
coast the nearest polygon is often across the curve, and once the transect is
in the lookup nothing downstream can tell: it just pulls that zone's mean
toward its neighbour's. A transect left out is visible in problems.csv and on
the map; a transect put in the wrong zone is visible nowhere. Fix the polygons
or the join point instead.

WHAT TO CHANGE
--------------
CONFIG, and nothing else unless you mean to change the method. If most of your
transects come back unmatched, look at the map before touching anything: the
usual causes are polygons that stop short of the transect origins (try
JOIN_POINT = "midpoint") or two files in different places entirely.

USAGE
    python transect_zone_join_template.py
    python transect_zone_join_template.py --transects t.geojson --zones z.geojson

INPUT this expects
    TRANSECTS_FILE   lines, one per transect, with an id column. The CoastSat
                     transect GeoJSON (coastsat.space) is already this shape;
                     its "id" is the same string as the time-series filename
                     (usa_NC_0032_0021 <-> usa_NC_0032_0021.csv), which is
                     what the other two templates match on.
    ZONES_FILE       polygons, one per zone, with an id column. Any format
                     geopandas reads (GeoJSON, shapefile, GeoPackage).

OUTPUT (in OUTPUT_DIR)
    transect_zones.csv            transect_id, zone_id   -- the lookup
    transect_zones_problems.csv   transect_id, problem, zone_ids
    transect_zones_map.png        points coloured matched / unmatched / ambiguous
==============================================================================
"""

from __future__ import annotations

import argparse
from pathlib import Path

import geopandas as gpd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
from shapely.geometry import Point

# =============================================================================
# CONFIG  -- the only part you should need to edit
# =============================================================================

TRANSECTS_FILE = Path("data/transects.geojson")
TRANSECT_ID_COL = "id"          # CoastSat's transect layer calls it "id"

ZONES_FILE = Path("data/zones.geojson")
ZONE_ID_COL = "zone_id"

OUTPUT_DIR = Path("data")       # the lookup lands where the rates template
                                # looks for it: data/transect_zones.csv

# Which single point of each transect line decides its zone.
#   "origin"    first vertex, the landward end CoastSat measures from
#   "midpoint"  halfway along the line; use it if your polygons are drawn
#               over the beach and stop short of the transect origins
JOIN_POINT = "origin"

# =============================================================================
# 1  LOAD
# =============================================================================

def load_layer(path: Path, id_col: str, what: str) -> gpd.GeoDataFrame:
    gdf = gpd.read_file(path)
    if gdf.crs is None:
        raise SystemExit(
            f"{what} file {path} has no CRS. Set it where the file was made "
            f"(or with gdf.set_crs) -- guessing one here would put every "
            f"point in the wrong place without an error.")
    if id_col not in gdf.columns:
        raise SystemExit(f"{what} file {path} has no column '{id_col}'. "
                         f"Columns: {list(gdf.columns)}")
    gdf = gdf[[id_col, gdf.geometry.name]].copy()
    gdf[id_col] = gdf[id_col].astype(str).str.strip()
    print(f"  {len(gdf):,} {what.lower()} from {path.name}  ({gdf.crs})")
    return gdf


# =============================================================================
# 2  POINT
# =============================================================================

def join_point(line, how: str) -> Point:
    if line.geom_type == "MultiLineString":
        line = max(line.geoms, key=lambda g: g.length)
    if how == "origin":
        return Point(line.coords[0])
    if how == "midpoint":
        return line.interpolate(0.5, normalized=True)
    raise SystemExit(f"JOIN_POINT must be 'origin' or 'midpoint', not {how!r}")


# =============================================================================
# the run
# =============================================================================

def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[2])
    ap.add_argument("--transects", type=Path, default=TRANSECTS_FILE)
    ap.add_argument("--zones", type=Path, default=ZONES_FILE)
    ap.add_argument("--out", type=Path, default=OUTPUT_DIR)
    ap.add_argument("--join-point", default=JOIN_POINT,
                    choices=("origin", "midpoint"))
    args = ap.parse_args()

    transects = load_layer(args.transects, TRANSECT_ID_COL, "Transects")
    zones = load_layer(args.zones, ZONE_ID_COL, "Zones")

    # Reproject to a metric CRS for the join. The result is the same in any
    # CRS; a metric one keeps the map honest.
    crs = zones.estimate_utm_crs()
    zones = zones.to_crs(crs).rename(columns={ZONE_ID_COL: "zone_id"})
    points = gpd.GeoDataFrame(
        {"transect_id": transects[TRANSECT_ID_COL].values},
        geometry=[join_point(g, args.join_point) for g in transects.geometry],
        crs=transects.crs).to_crs(crs)

    # A global CoastSat layer has hundreds of thousands of transects; keep the
    # ones near the zones before the join. 5 km of slack is only for speed --
    # the join itself decides membership.
    points = points[points.intersects(zones.union_all().envelope.buffer(5000))]
    print(f"  {len(points):,} transects near the zones")

    # 3  JOIN
    hits = gpd.sjoin(points, zones, how="left", predicate="within")
    per_transect = (hits.groupby("transect_id")["zone_id"]
                        .agg(lambda s: sorted(s.dropna().astype(str)))
                        .reset_index(name="zone_ids"))

    # 4  REFUSE -- exactly one zone, or not in the lookup at all
    n_zones = per_transect["zone_ids"].str.len()
    matched = per_transect[n_zones == 1].copy()
    matched["zone_id"] = matched["zone_ids"].str[0]
    problems = per_transect[n_zones != 1].copy()
    problems["problem"] = n_zones[n_zones != 1].map(
        lambda n: "no zone" if n == 0 else "more than one zone")
    problems["zone_ids"] = problems["zone_ids"].str.join(";")

    # 5  WRITE
    args.out.mkdir(parents=True, exist_ok=True)
    matched[["transect_id", "zone_id"]].to_csv(args.out / "transect_zones.csv",
                                               index=False)
    problems[["transect_id", "problem", "zone_ids"]].to_csv(
        args.out / "transect_zones_problems.csv", index=False)

    counts = matched["zone_id"].value_counts()
    empty = sorted(set(zones["zone_id"]) - set(counts.index))
    print(f"  matched {len(matched):,}   no zone "
          f"{(problems['problem'] == 'no zone').sum():,}   more than one zone "
          f"{(problems['problem'] == 'more than one zone').sum():,}")
    print(f"  transects per zone: min {counts.min() if len(counts) else 0}, "
          f"median {counts.median() if len(counts) else 0:.0f}, "
          f"max {counts.max() if len(counts) else 0}")
    if empty:
        print(f"  {len(empty)} zone(s) with NO transect: "
              f"{', '.join(empty[:10])}{' ...' if len(empty) > 10 else ''}")

    status = pd.Series("matched", index=per_transect["transect_id"])
    status[problems["transect_id"].values] = problems["problem"].values
    points["status"] = points["transect_id"].map(status)

    fig, ax = plt.subplots(figsize=(9, 9))
    zones.plot(ax=ax, facecolor="#f2efe6", edgecolor="0.55", lw=0.6)
    for zid, geom in zip(zones["zone_id"], zones.geometry):
        c = geom.representative_point()
        ax.annotate(zid, (c.x, c.y), fontsize=6, ha="center", color="0.35")
    style = {"matched": dict(color="#2166ac", markersize=4, marker="o"),
             "no zone": dict(color="#b2182b", markersize=14, marker="x"),
             "more than one zone": dict(color="#e08214", markersize=14, marker="D")}
    for label, kw in style.items():
        sel = points[points["status"] == label]
        if len(sel):
            sel.plot(ax=ax, label=f"{label} ({len(sel)})", zorder=3, **kw)
    ax.set_title(f"Transect {args.join_point} -> zone, point in polygon only")
    ax.legend(loc="best", fontsize=8)
    ax.set_axis_off()
    fig.tight_layout()
    fig.savefig(args.out / "transect_zones_map.png", dpi=200)
    plt.close(fig)

    print(f"  wrote transect_zones.csv, transect_zones_problems.csv and "
          f"transect_zones_map.png to {args.out.resolve()}")
    print("  Look at the map before using the lookup.")


if __name__ == "__main__":
    main()
