"""
CoastSat Transect → CASCADE Domain Spatial Mapping
===================================================
Step 1 of 2 in the domain-level LRR workflow.

This script spatially joins CoastSat transect geometry to your CASCADE
domain geometry, producing a lookup table:

    transect_id  |  domain_number  |  distance_to_domain_m  |  match_method

That lookup table is then used by coastsat_domain_lrr.py to compute
domain-level LRR summaries.

Inputs
------
  - CoastSat transects GeoJSON (downloaded from coastsat.space)
  - CASCADE domain polygons GeoJSON or shapefile (exported from ArcGIS)

Outputs
-------
  transect_domain_lookup.csv  –  saved to OUTPUT_DIR
  transect_domain_map.png     –  quick-look map to verify the join

Dependencies
------------
  pip install geopandas pandas numpy matplotlib shapely pyproj
"""

# ============================================================
# CONFIG  –  edit these paths and column names before running
# ============================================================

# --- CoastSat transect geometry ---
# GeoJSON downloaded from coastsat.space → "transects" button
# This may be the full global file — the script will automatically
# clip it to your study area using the domain bounding box.
TRANSECT_GEOM_PATH = r"C:\Users\hanna\PycharmProjects\CASCADE\scripts\input_preperation\CoastSat\CoastSat_transect_layer.geojson"

# Column in the transect file that holds the transect ID
# From the global CoastSat GeoJSON this is typically "id"
TRANSECT_ID_COL = "id"

# --- CASCADE domain geometry ---
# GeoJSON exported from ArcGIS (right-click layer → Data → Export Features → GeoJSON)
DOMAIN_GEOM_PATH = r"C:\Users\hanna\PycharmProjects\CASCADE\scripts\input_preperation\CoastSat\HAT_domains.json"

# Column in the domain file that holds the domain number
# From your attribute table this is "domain_id"
DOMAIN_ID_COL = "domain_id"

# --- Output ---
OUTPUT_DIR = r"."
LOOKUP_CSV = "transect_domain_lookup.csv"
MAP_PNG    = "transect_domain_map.png"

# Coordinate Reference System for distance calculations.
# UTM Zone 18N is correct for the NC Outer Banks.
PROJECTED_CRS = "EPSG:32618"

# Buffer (in degrees) added around the domain bounding box when
# pre-filtering the global transect file. 0.1 deg ~ 10 km.
BBOX_BUFFER_DEG = 0.1

# Maximum distance (m) to snap a transect to a domain if no
# point-in-polygon match is found.
MAX_SNAP_DISTANCE_M = 2000

# ============================================================
# IMPORTS
# ============================================================
import os
import warnings
import numpy as np
import pandas as pd
import geopandas as gpd
from shapely.geometry import Point, box
import matplotlib.pyplot as plt
warnings.filterwarnings("ignore")


# ============================================================
# FUNCTIONS
# ============================================================

def fix_crs(gdf: gpd.GeoDataFrame, fallback_crs: str = "EPSG:32618") -> gpd.GeoDataFrame:
    """
    Validate the CRS of a GeoDataFrame. If the CRS is missing or
    unrecognised (e.g. EPSG:3725 which is not a real EPSG code),
    replace it with fallback_crs and warn the user.

    This commonly happens when ArcGIS exports GeoJSON with a non-standard
    CRS code stored internally that geopandas cannot resolve.
    """
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


def load_transects(path: str, id_col: str) -> gpd.GeoDataFrame:
    """
    Load CoastSat transect geometry.
    If geometries are LineStrings, extracts the first point (transect origin)
    and sets it as the active geometry for point-in-polygon joining.
    The old LineString geometry column is dropped to avoid confusion.
    """
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

    # Extract transect origin point from LineString and make it the
    # ONLY geometry column (avoids the 'active geometry not present' error)
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


def load_domains(path: str, id_col: str) -> gpd.GeoDataFrame:
    """Load CASCADE domain geometry from GeoJSON or shapefile."""
    gdf = gpd.read_file(path)
    print(f"Loaded {len(gdf)} domains from: {os.path.basename(path)}")
    print(f"  Columns   : {list(gdf.columns)}")
    print(f"  CRS       : {gdf.crs}")
    print(f"  Geom types: {gdf.geom_type.value_counts().to_dict()}")

    # Fix invalid CRS — EPSG:3725 is not a real code; domains exported
    # from ArcGIS in UTM 18N should be EPSG:32618
    gdf = fix_crs(gdf, fallback_crs=PROJECTED_CRS)

    if id_col not in gdf.columns:
        raise ValueError(
            f"Column '{id_col}' not found in domain file.\n"
            f"Available columns: {list(gdf.columns)}"
        )
    print()
    return gdf


def clip_transects_to_study_area(transects: gpd.GeoDataFrame,
                                  domains: gpd.GeoDataFrame,
                                  buffer_deg: float) -> gpd.GeoDataFrame:
    """
    Pre-filter the global transect dataset to only those within a
    buffered bounding box around the domain extent.

    Essential when working with the full global CoastSat GeoJSON
    (230,000+ transects) — without this the spatial join is very slow.
    """
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


def spatial_join_polygon(transects: gpd.GeoDataFrame,
                          domains: gpd.GeoDataFrame,
                          t_id_col: str,
                          d_id_col: str,
                          proj_crs: str,
                          max_snap_m: float) -> pd.DataFrame:
    """
    Join transect origin points to domain polygons.

    Strategy:
      1. Reproject both layers to the projected CRS.
      2. Point-in-polygon exact match.
      3. For any unmatched transects, nearest-feature snap up to max_snap_m.
    """
    # Reproject — build clean GeoDataFrames with just the columns we need
    # so there's no ambiguity about which geometry column is active
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

    # --- Step 1: point-in-polygon ---
    joined = gpd.sjoin(t_proj, d_proj, how="left", predicate="within")
    joined = joined.rename(columns={d_id_col: "domain_number"})

    matched   = joined[joined["domain_number"].notna()].copy()
    unmatched = joined[joined["domain_number"].isna()][[t_id_col, "geometry"]].copy()
    print(f"  Point-in-polygon: {len(matched):,} matched, {len(unmatched):,} unmatched")

    # --- Step 2: nearest snap for unmatched ---
    snapped_rows = []
    if len(unmatched) > 0:
        for _, row in unmatched.iterrows():
            distances   = d_proj.geometry.distance(row.geometry)
            nearest_idx = distances.idxmin()
            dist        = distances[nearest_idx]
            matched_domain = d_proj.loc[nearest_idx, d_id_col] if dist <= max_snap_m else None
            snapped_rows.append({
                t_id_col       : row[t_id_col],
                "domain_number": matched_domain,
                "distance_m"   : round(dist, 1),
                "match_method" : "nearest_snap" if dist <= max_snap_m else "unmatched",
            })

        n_snap = sum(1 for r in snapped_rows if r["match_method"] == "nearest_snap")
        n_un   = sum(1 for r in snapped_rows if r["match_method"] == "unmatched")
        print(f"  Nearest snap:    {n_snap:,} additional matched")
        if n_un:
            print(f"  WARNING: {n_un:,} transects unmatched (>{max_snap_m} m from any domain)")

    # Build final lookup table
    matched = matched[[t_id_col, "domain_number"]].copy()
    matched["distance_m"]   = 0.0
    matched["match_method"] = "point_in_polygon"
    matched = matched.reset_index(drop=True)

    if snapped_rows:
        lookup = pd.concat([matched, pd.DataFrame(snapped_rows)], ignore_index=True)
    else:
        lookup = matched

    return lookup.rename(columns={t_id_col: "transect_id"})


def make_verification_map(transects: gpd.GeoDataFrame,
                           domains: gpd.GeoDataFrame,
                           lookup: pd.DataFrame,
                           t_id_col: str,
                           d_id_col: str,
                           out_path: str):
    """
    Quick-look map: transect origins coloured by match method,
    domain polygons shown in light yellow with domain ID labels.
    Inspect this carefully before using the lookup table.
    """
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


# ============================================================
# MAIN
# ============================================================

def main():
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    # ---- Load ----
    print("=" * 55)
    print("Loading CoastSat transects...")
    transects = load_transects(TRANSECT_GEOM_PATH, TRANSECT_ID_COL)

    print("Loading CASCADE domains...")
    domains = load_domains(DOMAIN_GEOM_PATH, DOMAIN_ID_COL)

    # ---- Clip global file to study area ----
    print("Clipping transects to study area bounding box...")
    transects = clip_transects_to_study_area(transects, domains, BBOX_BUFFER_DEG)

    # ---- Spatial join ----
    print("Running spatial join...")
    lookup = spatial_join_polygon(transects, domains,
                                  TRANSECT_ID_COL, DOMAIN_ID_COL,
                                  PROJECTED_CRS, MAX_SNAP_DISTANCE_M)

    # ---- Summary ----
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

    # ---- Save ----
    out_csv = os.path.join(OUTPUT_DIR, LOOKUP_CSV)
    lookup.to_csv(out_csv, index=False)
    print(f"Lookup table saved: {out_csv}")

    # ---- Verification map ----
    out_map = os.path.join(OUTPUT_DIR, MAP_PNG)
    make_verification_map(transects, domains, lookup,
                          TRANSECT_ID_COL, DOMAIN_ID_COL, out_map)

    return lookup


if __name__ == "__main__":
    lookup = main()
