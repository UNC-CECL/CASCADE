#!/usr/bin/env python3

"""
PEA ISLAND ROAD SETBACK CALCULATOR

Calculates NC-12 road setback for Pea Island CASCADE domains D80-D119 using
ONLY:
  1) NC-12 road GeoJSON
  2) domain polygons GeoJSON
  3) interior_start_by_along + trim_offset_c0 from the dune/topography picker JSON

No reference/previous road-setback values are used.
No road-mask values are used.
No Barrier3D/CASCADE source files are modified.

For each of the 50 alongshore profiles in a domain:
  - reconstruct the source-grid ocean-side cell centers in map coordinates,
  - locate the CASCADE interior-row-0 point from trim_offset_c0 + interior_start_by_along,
  - find the nearest point on the NC-12 road geometry,
  - measure the planar distance in meters.

The domain road setback is the mean of the valid profile distances.

A simulation-ready CSV is also written with exactly two columns:
    domain_id, road_setback_m
"""

from pathlib import Path
import json

import numpy as np
import pandas as pd
import geopandas as gpd
import matplotlib.pyplot as plt

from shapely.geometry import LineString, Point
from shapely.ops import nearest_points


# =============================================================================
# CONFIGURATION
# =============================================================================

PROJECT_ROOT = Path(
    "/Users/rsahrae/PycharmProjects/PeaIsland_Hindcast/CASCADE"
)

ROAD_GEOJSON = Path(
    "/Users/rsahrae/PycharmProjects/PeaIsland_Hindcast/"
    "CASCADE/data/PeaIsland_init/Roads/"
    "nc12road_1992__FeaturesToJSO.json"
)

DOMAIN_GEOJSON = Path(
    "/Users/rsahrae/PycharmProjects/PeaIsland_Hindcast/"
    "CASCADE/scripts/input_preperation/CoastSat/PeaIsland_domains.json"
)

INTERIOR_START_JSON = Path(
    "/Users/rsahrae/PycharmProjects/PeaIsland_Hindcast/"
    "CASCADE/data/PeaIsland_init/1996/dune_topo_new/picks/"
    "PEA_dune_search_windows_Hybrid_DEM.json"
)

OUTPUT_FOLDER = PROJECT_ROOT / "data" / "PeaIsland_init" / "road_setback"
OUTPUT_FOLDER.mkdir(parents=True, exist_ok=True)

FIRST_DOMAIN = 80
LAST_DOMAIN = 119
DOMAIN_IDS = list(range(FIRST_DOMAIN, LAST_DOMAIN + 1))
DOMAIN_ID_FIELD = "ID"

CELL_SIZE_M = 10.0
ALONGSHORE_CELLS = 50
DOMAIN_ALONGSHORE_LENGTH_M = ALONGSHORE_CELLS * CELL_SIZE_M

# Your extractor uses OCEAN_LOC="right", then reverses the cross-shore axis so
# index 0 becomes ocean-first. The interior JSON indices are therefore measured
# landward from the ocean side after trimming.
OCEAN_FIRST = True

# ArcGIS RasterToNumPyArray convention: array row 0 is the top/north side.
# The dune/topography workflow keeps ALONGSHORE_FLIP=False, so profile 0 is
# treated as the north-most 10-m profile within each domain.
PROFILE_ZERO_IS_NORTH = True


# =============================================================================
# HELPERS
# =============================================================================

def load_interior_json():
    if not INTERIOR_START_JSON.is_file():
        raise FileNotFoundError(
            f"Interior-start JSON not found:\n{INTERIOR_START_JSON}"
        )

    with open(INTERIOR_START_JSON, "r") as f:
        return json.load(f)


def polygon_principal_axes(polygon):
    """Return centroid, alongshore unit vector, and cross-shore unit vector.

    Uses PCA on polygon boundary coordinates instead of Shapely's oriented
    envelope. For these ~500 m x ~2000 m domain rectangles, the long principal
    axis is cross-shore and the short principal axis is alongshore.
    """
    if polygon is None or polygon.is_empty:
        raise ValueError("Empty domain polygon")

    coords = np.asarray(polygon.exterior.coords[:-1], dtype=float)
    if coords.shape[0] < 4:
        raise ValueError("Domain polygon has too few vertices")

    center = np.array([polygon.centroid.x, polygon.centroid.y], dtype=float)
    centered = coords - center

    cov = np.cov(centered.T)
    eigvals, eigvecs = np.linalg.eigh(cov)

    order = np.argsort(eigvals)
    alongshore_unit = eigvecs[:, order[0]]
    crossshore_unit = eigvecs[:, order[-1]]

    alongshore_unit = alongshore_unit / np.linalg.norm(alongshore_unit)
    crossshore_unit = crossshore_unit / np.linalg.norm(crossshore_unit)

    return center, alongshore_unit, crossshore_unit


def orient_axes(alongshore_unit, crossshore_unit):
    """Orient axes to match the Pea Island DEM/CASCADE extraction.

    IMPORTANT
    ---------
    The road is NOT used to decide which way is landward.

    For the source Pea Island GIS arrays used by the extractor:
        OCEAN_LOC = "right"

    Therefore, in map coordinates for these domains:
        east/right  = ocean
        west/left   = landward/sound

    CASCADE's ocean-first cross-shore index increases LANDWARD, so the positive
    cross-shore unit vector must point approximately west (negative x/Easting).

    Alongshore profile 0 is kept on the north/top side because
    ALONGSHORE_FLIP=False in the extractor.
    """

    # Positive cross-shore direction = landward = approximately west.
    if crossshore_unit[0] > 0:
        crossshore_unit = -crossshore_unit

    # Profile 0 -> 49 should run north -> south.
    if PROFILE_ZERO_IS_NORTH and alongshore_unit[1] > 0:
        alongshore_unit = -alongshore_unit

    return alongshore_unit, crossshore_unit


def build_profile_line(center, alongshore_unit, crossshore_unit, profile_index,
                       extension_m=4000.0):
    """Build the cross-shore line through one 10-m alongshore cell center."""
    offset_m = (
        (profile_index + 0.5) * CELL_SIZE_M
        - DOMAIN_ALONGSHORE_LENGTH_M / 2.0
    )

    profile_center = center + alongshore_unit * offset_m

    p0 = profile_center - crossshore_unit * extension_m
    p1 = profile_center + crossshore_unit * extension_m

    return LineString([tuple(p0), tuple(p1)])


def _geometry_endpoints(geom):
    """Collect useful endpoint Points from line/polygon intersection geometry."""
    points = []

    if geom.is_empty:
        return points

    gtype = geom.geom_type

    if gtype == "Point":
        points.append(geom)

    elif gtype == "MultiPoint":
        points.extend(list(geom.geoms))

    elif gtype == "LineString":
        coords = list(geom.coords)
        if coords:
            points.extend([Point(coords[0]), Point(coords[-1])])

    elif gtype == "MultiLineString":
        for part in geom.geoms:
            coords = list(part.coords)
            if coords:
                points.extend([Point(coords[0]), Point(coords[-1])])

    elif gtype == "GeometryCollection":
        for part in geom.geoms:
            points.extend(_geometry_endpoints(part))

    return points


def ocean_side_domain_entry(profile_line, polygon, crossshore_unit):
    """Return the ocean/seaward domain-boundary point on this profile.

    crossshore_unit points landward, so the boundary point with the smallest
    projection along this vector is the ocean-side edge.
    """
    intersection = profile_line.intersection(polygon)
    points = _geometry_endpoints(intersection)

    if not points:
        return None

    projections = []
    for p in points:
        xy = np.array([p.x, p.y], dtype=float)
        projections.append(float(np.dot(xy, crossshore_unit)))

    return points[int(np.argmin(projections))]


def reconstruct_interior_start_point(
    seaward_boundary_point,
    crossshore_unit,
    trim_offset_c0,
    interior_start_cell,
):
    """Convert JSON grid indices to the map-coordinate center of interior row 0.

    The extractor first makes the array ocean-first, then trims c0 ocean-side
    columns, and finally stores interior_start_by_along in the trimmed array.

    Therefore source ocean-first cell index = c0 + interior_start_cell.
    Add 0.5 cell because the polygon boundary is a cell EDGE while the raster
    index describes a cell CENTER.
    """
    source_cell_index = float(trim_offset_c0) + float(interior_start_cell)
    distance_from_ocean_edge_m = (source_cell_index + 0.5) * CELL_SIZE_M

    xy0 = np.array(
        [seaward_boundary_point.x, seaward_boundary_point.y],
        dtype=float,
    )

    xy = xy0 + crossshore_unit * distance_from_ocean_edge_m
    return Point(float(xy[0]), float(xy[1]))


# =============================================================================
# MAIN
# =============================================================================

def main():
    print("\n" + "=" * 90)
    print("PEA ISLAND ROAD SETBACK: ROAD JSON -> CASCADE INTERIOR START")
    print("=" * 90)

    for path, label in [
        (ROAD_GEOJSON, "Road GeoJSON"),
        (DOMAIN_GEOJSON, "Domain GeoJSON"),
        (INTERIOR_START_JSON, "Interior-start JSON"),
    ]:
        if not path.is_file():
            raise FileNotFoundError(f"{label} not found:\n{path}")

    road_gdf = gpd.read_file(ROAD_GEOJSON)
    domains_gdf = gpd.read_file(DOMAIN_GEOJSON)

    print("Road CRS:  ", road_gdf.crs)
    print("Domain CRS:", domains_gdf.crs)

    if road_gdf.crs is None:
        raise ValueError("Road GeoJSON has no CRS")
    if domains_gdf.crs is None:
        raise ValueError("Domain GeoJSON has no CRS")

    if road_gdf.crs != domains_gdf.crs:
        print("Reprojecting road to domain CRS...")
        road_gdf = road_gdf.to_crs(domains_gdf.crs)

    if domains_gdf.crs.is_geographic:
        raise ValueError(
            "Domain CRS is geographic (degrees). A projected meter CRS is required."
        )

    if DOMAIN_ID_FIELD not in domains_gdf.columns:
        raise KeyError(
            f"'{DOMAIN_ID_FIELD}' not found in domain polygons. "
            f"Available columns: {list(domains_gdf.columns)}"
        )

    road_geometry = road_gdf.geometry.union_all()
    interior_json = load_interior_json()

    profile_results = []
    domain_results = []

    for domain_id in DOMAIN_IDS:
        print("\n" + "-" * 90)
        print(f"Processing D{domain_id}")

        selection = domains_gdf[domains_gdf[DOMAIN_ID_FIELD] == domain_id]
        if selection.empty:
            raise ValueError(f"No domain polygon found for D{domain_id}")

        polygon = selection.geometry.iloc[0]
        if not polygon.is_valid:
            polygon = polygon.buffer(0)

        stem = f"domain_{domain_id}"
        if stem not in interior_json:
            raise KeyError(f"{stem} missing from interior-start JSON")

        entry = interior_json[stem]

        if "interior_start_by_along" not in entry:
            raise KeyError(
                f"{stem} does not contain 'interior_start_by_along'"
            )

        if "trim_offset_c0" not in entry:
            raise KeyError(
                f"{stem} does not contain 'trim_offset_c0'. "
                "This is required to geolocate the trimmed-array indices."
            )

        interior_start = np.asarray(
            entry["interior_start_by_along"], dtype=float
        )
        trim_offset_c0 = float(entry["trim_offset_c0"])

        if len(interior_start) != ALONGSHORE_CELLS:
            raise ValueError(
                f"D{domain_id}: expected {ALONGSHORE_CELLS} interior-start "
                f"values, got {len(interior_start)}"
            )

        center, alongshore_unit, crossshore_unit = polygon_principal_axes(polygon)
        alongshore_unit, crossshore_unit = orient_axes(
            alongshore_unit,
            crossshore_unit,
        )

        distances = []

        for profile in range(ALONGSHORE_CELLS):
            interior_cell = interior_start[profile]

            if not np.isfinite(interior_cell) or interior_cell < 0:
                profile_results.append(
                    {
                        "domain": domain_id,
                        "profile": profile,
                        "trim_offset_c0": trim_offset_c0,
                        "interior_start_cell": interior_cell,
                        "setback_m": np.nan,
                        "valid": False,
                        "reason": "invalid interior-start index",
                    }
                )
                continue

            profile_line = build_profile_line(
                center,
                alongshore_unit,
                crossshore_unit,
                profile,
            )

            seaward_point = ocean_side_domain_entry(
                profile_line,
                polygon,
                crossshore_unit,
            )

            if seaward_point is None:
                profile_results.append(
                    {
                        "domain": domain_id,
                        "profile": profile,
                        "trim_offset_c0": trim_offset_c0,
                        "interior_start_cell": interior_cell,
                        "setback_m": np.nan,
                        "valid": False,
                        "reason": "no profile/domain intersection",
                    }
                )
                continue

            interior_point = reconstruct_interior_start_point(
                seaward_point,
                crossshore_unit,
                trim_offset_c0,
                interior_cell,
            )

            # Pure geometric distance from the reconstructed interior-row-0 point
            # to the NC-12 road JSON. No previous/reference setback values used.
            nearest_road_point = nearest_points(interior_point, road_geometry)[1]
            setback_m = float(interior_point.distance(nearest_road_point))

            # Also save the signed cross-shore component for QC only.
            road_vector = np.array(
                [
                    nearest_road_point.x - interior_point.x,
                    nearest_road_point.y - interior_point.y,
                ],
                dtype=float,
            )
            crossshore_component_m = float(
                np.dot(road_vector, crossshore_unit)
            )

            distances.append(setback_m)

            profile_results.append(
                {
                    "domain": domain_id,
                    "profile": profile,
                    "trim_offset_c0": trim_offset_c0,
                    "interior_start_cell": interior_cell,
                    "source_oceanfirst_cell": trim_offset_c0 + interior_cell,
                    "interior_x": interior_point.x,
                    "interior_y": interior_point.y,
                    "road_x": nearest_road_point.x,
                    "road_y": nearest_road_point.y,
                    "setback_m": setback_m,
                    "crossshore_component_m": crossshore_component_m,
                    "valid": True,
                    "reason": "",
                }
            )

        distances = np.asarray(distances, dtype=float)
        valid = np.isfinite(distances)
        n_valid = int(valid.sum())

        if n_valid:
            mean_setback = float(np.mean(distances[valid]))
            median_setback = float(np.median(distances[valid]))
            min_setback = float(np.min(distances[valid]))
            max_setback = float(np.max(distances[valid]))
            std_setback = float(np.std(distances[valid]))
        else:
            mean_setback = np.nan
            median_setback = np.nan
            min_setback = np.nan
            max_setback = np.nan
            std_setback = np.nan

        domain_results.append(
            {
                "domain": domain_id,
                "valid_profiles": n_valid,
                "total_profiles": ALONGSHORE_CELLS,
                "trim_offset_c0": trim_offset_c0,
                "mean_setback_m": mean_setback,
                "median_setback_m": median_setback,
                "min_setback_m": min_setback,
                "max_setback_m": max_setback,
                "std_setback_m": std_setback,
            }
        )

        print(f"  valid profiles : {n_valid}/{ALONGSHORE_CELLS}")
        if np.isfinite(mean_setback):
            print(f"  mean setback   : {mean_setback:.1f} m")
            print(f"  median         : {median_setback:.1f} m")
            print(f"  range          : {min_setback:.1f} - {max_setback:.1f} m")
            print(f"  trim_offset_c0 : {trim_offset_c0:.0f} cells")

    # =========================================================================
    # SAVE
    # =========================================================================

    profile_df = pd.DataFrame(profile_results)
    domain_df = pd.DataFrame(domain_results)

    profile_csv = OUTPUT_FOLDER / "road_setback_profiles_from_road_json_D80_D119.csv"
    domain_csv = OUTPUT_FOLDER / "road_setback_mean_from_road_json_D80_D119.csv"

    # Clean two-column file intended for direct use in simulation setup.
    simulation_csv = OUTPUT_FOLDER / "road_setback_SIMULATION_D80_D119.csv"

    npy_file = OUTPUT_FOLDER / "road_setback_mean_from_road_json_m_D80_D119.npy"
    plot_file = OUTPUT_FOLDER / "road_setback_from_road_json_D80_D119.png"

    profile_df.to_csv(profile_csv, index=False)
    domain_df.to_csv(domain_csv, index=False)

    simulation_df = (
        domain_df[["domain", "mean_setback_m"]]
        .rename(
            columns={
                "domain": "domain_id",
                "mean_setback_m": "road_setback_m",
            }
        )
        .sort_values("domain_id")
        .reset_index(drop=True)
    )

    # Keep useful numerical precision while remaining easy to inspect/edit.
    simulation_df.to_csv(
        simulation_csv,
        index=False,
        float_format="%.3f",
    )

    setback_array = simulation_df["road_setback_m"].to_numpy(dtype=float)
    np.save(npy_file, setback_array)

    fig, ax = plt.subplots(figsize=(15, 6))
    ax.plot(
        domain_df["domain"],
        domain_df["mean_setback_m"],
        marker="o",
        linewidth=1.5,
    )
    ax.set_xlabel("Pea Island domain")
    ax.set_ylabel("Mean distance: interior row 0 to NC-12 (m)")
    ax.set_title("NC-12 Road Setback from Reconstructed CASCADE Interior Start")
    ax.set_xticks(DOMAIN_IDS)
    ax.tick_params(axis="x", rotation=90)
    ax.grid(alpha=0.25)
    plt.tight_layout()
    fig.savefig(plot_file, dpi=200, bbox_inches="tight")
    plt.close(fig)

    print("\n" + "=" * 90)
    print("CALCULATED ROAD SETBACKS — NO REFERENCE VALUES USED")
    print("=" * 90)

    setback_dict = {
        int(row["domain"]): (
            round(float(row["mean_setback_m"]), 1)
            if np.isfinite(row["mean_setback_m"])
            else None
        )
        for _, row in domain_df.iterrows()
    }

    print("\nDictionary:")
    print(setback_dict)

    print("\nList D80 -> D119:")
    print(setback_array.tolist())

    print("\nSaved:")
    print(profile_csv)
    print(domain_csv)
    print(simulation_csv)
    print(npy_file)
    print(plot_file)

    print("\nSimulation-ready CSV preview:")
    print(simulation_df.to_string(index=False))


if __name__ == "__main__":
    main()
