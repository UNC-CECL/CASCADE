"""
Step 1: assign each transect to one zone (point in polygon).

Writes the lookup the step-2 scripts read. A transect inside no zone, or
inside two, is left out and listed in the problems file -- never snapped.

    python transect_zone_join_template.py
    python transect_zone_join_template.py --transects t.geojson --zones z.geojson

Needs geopandas, matplotlib.

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-09-30
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

# --- CONFIG ------------------------------------------------------------------
TRANSECTS_FILE  = Path("example/data/masonboro_transects.geojson")
TRANSECT_ID_COL = "id"                 # CoastSat id = time-series filename
ZONES_FILE      = Path("example/data/masonboro_domains.geojson")
ZONE_ID_COL     = "domainID"
OUTPUT_DIR      = Path("example/data")
JOIN_POINT      = "origin"             # or "midpoint"
# -----------------------------------------------------------------------------


# Read a layer; refuse a missing CRS or id column
def load_layer(path: Path, id_col: str, what: str) -> gpd.GeoDataFrame:
    gdf = gpd.read_file(path)
    if gdf.crs is None:
        raise SystemExit(f"{path} has no CRS; set it at the source, don't guess")
    if id_col not in gdf.columns:
        raise SystemExit(f"{path} has no column '{id_col}': {list(gdf.columns)}")
    gdf = gdf[[id_col, gdf.geometry.name]].copy()
    gdf[id_col] = gdf[id_col].astype(str).str.strip()
    print(f"  {len(gdf):,} {what} from {path.name} ({gdf.crs})")
    return gdf


# One point per transect line: origin or midpoint
def join_point(line, how: str) -> Point:
    if line.geom_type == "MultiLineString":
        line = max(line.geoms, key=lambda g: g.length)
    return Point(line.coords[0]) if how == "origin" else line.interpolate(0.5, normalized=True)


# Map of zones and transect points, coloured by match status
def draw_map(zones, points, how, path):
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
    ax.set_title(f"Transect {how} -> zone")
    ax.legend(fontsize=8)
    ax.set_axis_off()
    fig.tight_layout()
    fig.savefig(path, dpi=200)
    plt.close(fig)


# Run: load, join, sort matches from problems, write
def main() -> None:
    ap = argparse.ArgumentParser(description="Assign transects to zones.")
    ap.add_argument("--transects", type=Path, default=TRANSECTS_FILE)
    ap.add_argument("--zones", type=Path, default=ZONES_FILE)
    ap.add_argument("--out", type=Path, default=OUTPUT_DIR)
    ap.add_argument("--join-point", default=JOIN_POINT, choices=("origin", "midpoint"))
    args = ap.parse_args()

    # Load both layers into one metric CRS
    transects = load_layer(args.transects, TRANSECT_ID_COL, "transects")
    zones = load_layer(args.zones, ZONE_ID_COL, "zones")
    crs = zones.estimate_utm_crs()
    zones = zones.to_crs(crs).rename(columns={ZONE_ID_COL: "zone_id"})
    points = gpd.GeoDataFrame(
        {"transect_id": transects[TRANSECT_ID_COL].values},
        geometry=[join_point(g, args.join_point) for g in transects.geometry],
        crs=transects.crs).to_crs(crs)
    points = points[points.intersects(zones.union_all().envelope.buffer(5000))]  # speed only

    # Point in polygon; keep transects in exactly one zone
    hits = gpd.sjoin(points, zones, how="left", predicate="within")
    per_transect = (hits.groupby("transect_id")["zone_id"]
                        .agg(lambda s: sorted(s.dropna().astype(str)))
                        .reset_index(name="zone_ids"))
    n_zones = per_transect["zone_ids"].str.len()
    matched = per_transect[n_zones == 1].assign(zone_id=lambda d: d["zone_ids"].str[0])
    problems = per_transect[n_zones != 1].assign(
        problem=n_zones[n_zones != 1].map(lambda n: "no zone" if n == 0 else "more than one zone"),
        zone_ids=lambda d: d["zone_ids"].str.join(";"))

    # Write the lookup and the problems list
    args.out.mkdir(parents=True, exist_ok=True)
    matched[["transect_id", "zone_id"]].to_csv(args.out / "transect_zones.csv", index=False)
    problems[["transect_id", "problem", "zone_ids"]].to_csv(
        args.out / "transect_zones_problems.csv", index=False)

    # Draw the check map
    status = pd.Series("matched", index=per_transect["transect_id"])
    status[problems["transect_id"].values] = problems["problem"].values
    points["status"] = points["transect_id"].map(status)
    draw_map(zones, points, args.join_point, args.out / "transect_zones_map.png")

    # Report counts and zones left empty
    empty = sorted(set(zones["zone_id"]) - set(matched["zone_id"]))
    print(f"  matched {len(matched):,}, "
          + ", ".join(f"{k} {v}" for k, v in problems["problem"].value_counts().items()))
    if empty:
        print(f"  zones with no transect: {', '.join(empty)}")
    print(f"  wrote lookup, problems and map to {args.out.resolve()} -- check the map")


if __name__ == "__main__":
    main()
