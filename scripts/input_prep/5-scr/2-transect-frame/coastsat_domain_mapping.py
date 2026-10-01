"""
Step 1 of 2: join CoastSat transects to the CASCADE domains, the lookup every rate script uses.

    python scripts/input_prep/5-scr/2-transect-frame/coastsat_domain_mapping.py

Writes transect_domain_lookup.csv and a map of the join. Details: scripts/input_prep/5-scr/2-transect-frame/README.md.

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-09-22
"""

# --- CONFIG ------------------------------------------------------------------
# CoastSat transect geometry

# Downloaded from coastsat.space; clipped to the study area automatically
import sys as _sys
from pathlib import Path as _RP
_sys.path.insert(0, str(next(_q for _q in _RP(__file__).resolve().parents
                             if (_q / "pyproject.toml").exists()) / "scripts"))
from site_layer import hat_observed_rates as _obs  # noqa: E402
TRANSECT_GEOM_PATH = str(_obs.TRANSECT_LAYER)

# Column in the transect file that holds the transect ID From the global CoastSat GeoJSON this is typically "id"
TRANSECT_ID_COL = "id"

# CASCADE domain geometry

# Exported from ArcGIS
DOMAIN_GEOM_PATH = str(_obs.DOMAIN_BOXES)

# Column in the domain file that holds the domain number From your attribute table this is "domain_id"
DOMAIN_ID_COL = "domain_id"

# Output
OUTPUT_DIR = str(_obs.TRANSECT_DOMAINS)
LOOKUP_CSV = "transect_domain_lookup.csv"
MAP_PNG    = "transect_domain_map.png"

# Coordinate Reference System for distance calculations
PROJECTED_CRS = "EPSG:32618"

# Buffer (in degrees) added around the domain bounding box when pre-filtering the global transect file
BBOX_BUFFER_DEG = 0.1

# Nearest-feature snapping is disabled — Hatteras Island's curvature causes cross-curve mismatches
MAX_SNAP_DISTANCE_M = None  # disabled
# -----------------------------------------------------------------------------

import os
import warnings
import numpy as np
import pandas as pd
import geopandas as gpd
from shapely.geometry import Point, box
import matplotlib.pyplot as plt
warnings.filterwarnings("ignore")


# Validate the CRS of a GeoDataFrame
def fix_crs(gdf: gpd.GeoDataFrame, fallback_crs: str = "EPSG:32618") -> gpd.GeoDataFrame:
    try:
        if gdf.crs is None:
            raise ValueError("No CRS defined")
        # Check CRS is resolvable by trying to get its EPSG code
        gdf.crs.to_authority()
        return gdf
    except Exception:
        print(f"  WARNING: CRS '{gdf.crs}' is unrecognised or invalid.")
        print(f"           Assuming {fallback_crs} (UTM Zone 18N).")
        print(f"           Verify this is correct before trusting the join.\n")
        return gdf.set_crs(fallback_crs, allow_override=True)


# Load CoastSat transect geometry
def load_transects(path: str, id_col: str) -> gpd.GeoDataFrame:
    gdf = gpd.read_file(path)
    print(f"Loaded {len(gdf):,} transects from: {os.path.basename(path)}")
    print(f"  Columns   : {list(gdf.columns)}")
    print(f"  CRS       : {gdf.crs}")
    print(f"  Geom types: {gdf.geom_type.value_counts().to_dict()}")

    gdf = fix_crs(gdf, fallback_crs="EPSG:4326")

    if id_col not in gdf.columns:
        raise ValueError(
            f"Column '{id_col}' not found in transect file.\n"
            f"Available columns: {list(gdf.columns)}"
        )

    # The transect origin as the only geometry column
    if gdf.geom_type.isin(["LineString", "MultiLineString"]).any():
        print("  → LineString detected; extracting origin point as active geometry.")
        origin_pts = gdf.geometry.apply(
            lambda g: Point(g.coords[0]) if g.geom_type == "LineString"
                      else Point(g.geoms[0].coords[0])
        )
        crs = gdf.crs
        # Keep only non-geometry attribute columns + rebuild as clean GeoDataFrame
        attr_cols = [c for c in gdf.columns if c != gdf.geometry.name]
        gdf = gdf[attr_cols].copy()
        gdf = gpd.GeoDataFrame(gdf, geometry=origin_pts, crs=crs)

    print()
    return gdf


# Load CASCADE domain geometry from GeoJSON or shapefile
def load_domains(path: str, id_col: str) -> gpd.GeoDataFrame:
    gdf = gpd.read_file(path)
    print(f"Loaded {len(gdf)} domains from: {os.path.basename(path)}")
    print(f"  Columns   : {list(gdf.columns)}")
    print(f"  CRS       : {gdf.crs}")
    print(f"  Geom types: {gdf.geom_type.value_counts().to_dict()}")

    # Replace the invalid EPSG:3725 with UTM 18N (EPSG:32618)
    gdf = fix_crs(gdf, fallback_crs=PROJECTED_CRS)

    if id_col not in gdf.columns:
        raise ValueError(
            f"Column '{id_col}' not found in domain file.\n"
            f"Available columns: {list(gdf.columns)}"
        )
    print()
    return gdf


# Transects within a buffered box around the domain extent, so the global file joins in reasonable time
def clip_transects_to_study_area(transects: gpd.GeoDataFrame,
                                  domains: gpd.GeoDataFrame,
                                  buffer_deg: float) -> gpd.GeoDataFrame:
    domains_wgs84 = domains.to_crs("EPSG:4326")
    minx, miny, maxx, maxy = domains_wgs84.total_bounds
    study_bbox = box(minx - buffer_deg, miny - buffer_deg,
                     maxx + buffer_deg, maxy + buffer_deg)

    transects_wgs84 = transects.to_crs("EPSG:4326")
    in_bbox = transects_wgs84.within(study_bbox)
    clipped = transects[in_bbox].copy().reset_index(drop=True)

    print(f"  Study area bbox: ({minx:.3f}, {miny:.3f}) → ({maxx:.3f}, {maxy:.3f})")
    print(f"  Clipped {len(transects):,} → {len(clipped):,} transects\n")

    if len(clipped) == 0:
        raise RuntimeError(
            "No transects found within the study area bounding box.\n"
            "Check that TRANSECT_GEOM_PATH and DOMAIN_GEOM_PATH are correct,\n"
            "and that both cover the same geographic area."
        )
    return clipped


# Join transect origin points to domain polygons — point-in-polygon only
def spatial_join_polygon(transects: gpd.GeoDataFrame,
                          domains: gpd.GeoDataFrame,
                          t_id_col: str,
                          d_id_col: str,
                          proj_crs: str,
                          max_snap_m: float = None) -> pd.DataFrame:
    # Reproject to projected CRS for accurate geometry operations
    t_proj = gpd.GeoDataFrame(
        {t_id_col: transects[t_id_col]},
        geometry=transects.geometry,
        crs=transects.crs
    ).to_crs(proj_crs)

    d_proj = gpd.GeoDataFrame(
        {d_id_col: domains[d_id_col]},
        geometry=domains.geometry,
        crs=domains.crs
    ).to_crs(proj_crs)

    # Point-in-polygon join — only exact matches kept
    joined = gpd.sjoin(t_proj, d_proj, how="left", predicate="within")
    joined = joined.rename(columns={d_id_col: "domain_number"})

    matched   = joined[joined["domain_number"].notna()].copy()
    unmatched = joined[joined["domain_number"].isna()].copy()

    print(f"  Point-in-polygon: {len(matched):,} matched")
    if len(unmatched) > 0:
        print(f"  Unmatched (excluded): {len(unmatched):,} transects")
        print(f"  These transects do not fall within any domain polygon and")
        print(f"  will not appear in the lookup table or downstream analysis.")

    # Build lookup — matched transects only
    lookup = matched[[t_id_col, "domain_number"]].copy()
    lookup["distance_m"]   = 0.0
    lookup["match_method"] = "point_in_polygon"
    lookup = lookup.reset_index(drop=True)

    return lookup.rename(columns={t_id_col: "transect_id"})


# Quick-look map
def make_verification_map(transects: gpd.GeoDataFrame,
                           domains: gpd.GeoDataFrame,
                           lookup: pd.DataFrame,
                           t_id_col: str,
                           d_id_col: str,
                           out_path: str):
    fig, ax = plt.subplots(1, 1, figsize=(12, 10))

    # Domain polygons
    d_plot = domains.to_crs("EPSG:4326")
    d_plot.plot(ax=ax, color="lightyellow", edgecolor="grey",
                linewidth=0.8, alpha=0.7, zorder=1)
    for _, row in d_plot.iterrows():
        c = row.geometry.centroid
        ax.annotate(str(row[d_id_col]), xy=(c.x, c.y),
                    fontsize=7, ha="center", color="dimgrey", zorder=2)

    # Transect origin points coloured by match method
    t_plot = gpd.GeoDataFrame(
        {t_id_col: transects[t_id_col]},
        geometry=transects.geometry,
        crs=transects.crs
    ).to_crs("EPSG:4326")
    t_plot = t_plot.merge(
        lookup[["transect_id", "domain_number", "match_method"]],
        left_on=t_id_col, right_on="transect_id", how="inner"
    )

    style = {
        "point_in_polygon": dict(color="steelblue", markersize=6,  marker="o", label="Point-in-polygon"),
        "nearest_snap"    : dict(color="orange",    markersize=8,  marker="o", label="Nearest snap"),
        "unmatched"       : dict(color="red",       markersize=10, marker="x", label="Unmatched"),
    }
    for method, kwargs in style.items():
        mask = t_plot["match_method"] == method
        if mask.any():
            t_plot[mask].plot(ax=ax, zorder=3, **kwargs)

    ax.set_title("CoastSat Transects → CASCADE Domain Assignment", fontsize=12)
    ax.set_xlabel("Longitude")
    ax.set_ylabel("Latitude")
    ax.legend(fontsize=9, loc="upper left")
    ax.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(out_path, dpi=150)
    plt.show()
    print(f"Verification map saved: {out_path}")


# Run: load both layers, join, write the lookup and the map
def main():
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    # Load
    print("=" * 55)
    print("Loading CoastSat transects...")
    transects = load_transects(TRANSECT_GEOM_PATH, TRANSECT_ID_COL)

    print("Loading CASCADE domains...")
    domains = load_domains(DOMAIN_GEOM_PATH, DOMAIN_ID_COL)

    # Clip global file to study area
    print("Clipping transects to study area bounding box...")
    transects = clip_transects_to_study_area(transects, domains, BBOX_BUFFER_DEG)

    # Spatial join
    print("Running spatial join...")
    lookup = spatial_join_polygon(transects, domains,
                                  TRANSECT_ID_COL, DOMAIN_ID_COL,
                                  PROJECTED_CRS, MAX_SNAP_DISTANCE_M)

    # Summary
    print(f"\n{'='*55}")
    print(f"  Total transects    : {len(lookup):,}")
    print(f"  Matched            : {lookup['domain_number'].notna().sum():,}")
    print(f"  Unmatched          : {lookup['domain_number'].isna().sum():,}")
    print(f"  Unique domains hit : {lookup['domain_number'].nunique()}")
    print(f"\n  Transects per domain:")
    counts = (lookup[lookup["domain_number"].notna()]
              .groupby("domain_number")["transect_id"]
              .count().rename("n_transects").sort_index())
    print(counts.to_string())
    print(f"{'='*55}\n")

    # Save
    out_csv = os.path.join(OUTPUT_DIR, LOOKUP_CSV)
    lookup.to_csv(out_csv, index=False)
    print(f"Lookup table saved: {out_csv}")

    # Verification map
    out_map = os.path.join(OUTPUT_DIR, MAP_PNG)
    make_verification_map(transects, domains, lookup,
                          TRANSECT_ID_COL, DOMAIN_ID_COL, out_map)

    return lookup


if __name__ == "__main__":
    lookup = main()
