# =============================================================================
# Calculate 1992–1996 road relocation distance by CASCADE domain
#
# Inputs:
#   1. Complete 1992 road GeoJSON
#   2. Complete 1996 road GeoJSON
#   3. CASCADE domain polygon GeoJSON
#
# The script intersects the 1992 road with each domain in memory, samples
# points along that segment, and measures distance to the complete 1996 road.
# =============================================================================

from pathlib import Path

import geopandas as gpd
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from shapely.geometry import (
    GeometryCollection,
    LineString,
    MultiLineString,
)
from shapely.ops import unary_union


# =============================================================================
# USER SETTINGS
# =============================================================================

ROAD_1992_FILE = Path(
    r"/Users/rsahrae/PycharmProjects/PeaIsland_Hindcast/CASCADE/data/PeaIsland_init/Roads/nc12road_1992.geojson"
)

ROAD_1996_FILE = Path(
    r"/Users/rsahrae/PycharmProjects/PeaIsland_Hindcast/CASCADE/data/PeaIsland_init/Roads/nc12-hat_present.geojson"
)

DOMAIN_FILE = Path(
    r"/Users/rsahrae/PycharmProjects/PeaIsland_Hindcast/CASCADE/scripts/input_preperation/CoastSat/PeaIsland_domains.json"
)

OUTPUT_CSV = Path(
    r"/Users/rsahrae/PycharmProjects/PeaIsland_Hindcast/CASCADE/data/PeaIsland_init/road_relocation_dis/road_relocation_1992_1996.csv"
)

OUTPUT_SAMPLE_POINTS = Path(
    r"/Users/rsahrae/path/to/road_relocation_sample_points.geojson"
)

OUTPUT_FIGURE = Path(
    r"/Users/rsahrae/PycharmProjects/PeaIsland_Hindcast/CASCADE/data/PeaIsland_init/road_relocation_dis/road_relocation_1992_1996.png"
)

# This field must be in the DOMAIN polygon GeoJSON, not the road file.
DOMAIN_ID_FIELD = "ID"

# Domains where you want to calculate relocation.
# Use None to analyze all domains.
DOMAIN_NUMBERS = [
    80,
    81,
    82,
    83,
    84,
    85,
    86,
    87,
    88,
    89,
    90,
    100,
    101,
    102,
    103,
    104,
    105,
    106,
    107,
    108,
    109,
    110,
    111,
    112,
    113,
    114,
    115,
    116,
    117,
    118,
    119
]

# Distance between points sampled along the 1992 road
SAMPLE_SPACING_M = 5.0

# NAD83(2011) / UTM Zone 18N
TARGET_CRS = "EPSG:6347"


# =============================================================================
# FUNCTIONS
# =============================================================================

def extract_lines(geometry):
    """
    Extract LineString objects from different geometry types.
    """

    lines = []

    if geometry is None or geometry.is_empty:
        return lines

    if isinstance(geometry, LineString):
        lines.append(geometry)

    elif isinstance(geometry, MultiLineString):
        lines.extend(list(geometry.geoms))

    elif isinstance(geometry, GeometryCollection):
        for part in geometry.geoms:
            lines.extend(extract_lines(part))

    return lines


def sample_line_geometry(geometry, spacing):
    """
    Generate regularly spaced sample points along line geometry.
    """

    sample_points = []

    for line in extract_lines(geometry):

        if line.length == 0:
            continue

        sample_distances = np.arange(
            0,
            line.length,
            spacing,
        )

        # Include the final endpoint
        if (
            len(sample_distances) == 0
            or sample_distances[-1] < line.length
        ):
            sample_distances = np.append(
                sample_distances,
                line.length,
            )

        for distance in sample_distances:
            sample_points.append(
                line.interpolate(distance)
            )

    return sample_points


# =============================================================================
# CREATE OUTPUT DIRECTORIES
# =============================================================================

OUTPUT_CSV.parent.mkdir(parents=True, exist_ok=True)
OUTPUT_SAMPLE_POINTS.parent.mkdir(parents=True, exist_ok=True)
OUTPUT_FIGURE.parent.mkdir(parents=True, exist_ok=True)


# =============================================================================
# READ INPUT FILES
# =============================================================================

road_1992 = gpd.read_file(ROAD_1992_FILE)
road_1996 = gpd.read_file(ROAD_1996_FILE)
domains = gpd.read_file(DOMAIN_FILE)

print("=" * 70)
print("INPUT FILE INFORMATION")
print("=" * 70)

print("\n1992 road columns:")
print(road_1992.columns.tolist())

print("\n1996 road columns:")
print(road_1996.columns.tolist())

print("\nDomain polygon columns:")
print(domains.columns.tolist())

print("\nOriginal coordinate systems:")
print(f"1992 road: {road_1992.crs}")
print(f"1996 road: {road_1996.crs}")
print(f"Domains:   {domains.crs}")

print("\nNumber of features:")
print(f"1992 road features: {len(road_1992)}")
print(f"1996 road features: {len(road_1996)}")
print(f"Domain features:    {len(domains)}")


# =============================================================================
# CHECK INPUT DATA
# =============================================================================

if road_1992.empty:
    raise ValueError("The 1992 road file is empty.")

if road_1996.empty:
    raise ValueError("The 1996 road file is empty.")

if domains.empty:
    raise ValueError("The domain polygon file is empty.")

if road_1992.crs is None:
    raise ValueError("The 1992 road file does not have a defined CRS.")

if road_1996.crs is None:
    raise ValueError("The 1996 road file does not have a defined CRS.")

if domains.crs is None:
    raise ValueError("The domain file does not have a defined CRS.")

if DOMAIN_ID_FIELD not in domains.columns:
    raise ValueError(
        f"\nThe field '{DOMAIN_ID_FIELD}' was not found in "
        f"the DOMAIN file.\n\n"
        f"Available domain fields are:\n"
        f"{domains.columns.tolist()}\n\n"
        f"Change DOMAIN_ID_FIELD to the field containing "
        f"numbers 80–119."
    )


# =============================================================================
# REPROJECT EVERYTHING TO THE SAME CRS
# =============================================================================

road_1992 = road_1992.to_crs(TARGET_CRS)
road_1996 = road_1996.to_crs(TARGET_CRS)
domains = domains.to_crs(TARGET_CRS)

print("\nProjected coordinate systems:")
print(f"1992 road: {road_1992.crs}")
print(f"1996 road: {road_1996.crs}")
print(f"Domains:   {domains.crs}")

print("\nDistances will be calculated in meters.")


# =============================================================================
# CLEAN DOMAIN NUMBERS
# =============================================================================

domains[DOMAIN_ID_FIELD] = pd.to_numeric(
    domains[DOMAIN_ID_FIELD],
    errors="coerce",
)

invalid_domains = domains[DOMAIN_ID_FIELD].isna().sum()

if invalid_domains > 0:
    print(
        f"\nRemoving {invalid_domains} domain polygon(s) "
        "with invalid domain IDs."
    )

domains = domains.dropna(
    subset=[DOMAIN_ID_FIELD]
).copy()

domains[DOMAIN_ID_FIELD] = (
    domains[DOMAIN_ID_FIELD].astype(int)
)

available_domains = sorted(
    domains[DOMAIN_ID_FIELD].unique().tolist()
)

print("\nAvailable domain numbers:")
print(available_domains)


# =============================================================================
# SELECT DOMAINS
# =============================================================================

if DOMAIN_NUMBERS is not None:

    missing_domains = sorted(
        set(DOMAIN_NUMBERS) - set(available_domains)
    )

    if missing_domains:
        print(
            "\nWarning: these requested domain numbers were "
            f"not found:\n{missing_domains}"
        )

    selected_domains = domains[
        domains[DOMAIN_ID_FIELD].isin(DOMAIN_NUMBERS)
    ].copy()

else:
    selected_domains = domains.copy()


if selected_domains.empty:
    raise ValueError(
        "No domains were selected. Check DOMAIN_ID_FIELD "
        "and DOMAIN_NUMBERS."
    )


# =============================================================================
# COMBINE ROAD FEATURES
# =============================================================================

# Your road files contain multiple features, so combine them into
# one road geometry for each year.
road_1992_geometry = unary_union(
    road_1992.geometry.dropna().tolist()
)

road_1996_geometry = unary_union(
    road_1996.geometry.dropna().tolist()
)

if road_1992_geometry.is_empty:
    raise ValueError("The combined 1992 road geometry is empty.")

if road_1996_geometry.is_empty:
    raise ValueError("The combined 1996 road geometry is empty.")


# =============================================================================
# CALCULATE RELOCATION DISTANCE BY DOMAIN
# =============================================================================

results = []
sample_records = []

print("\n" + "=" * 70)
print("ROAD RELOCATION RESULTS")
print("=" * 70)

for _, domain_row in selected_domains.iterrows():

    domain_number = domain_row[DOMAIN_ID_FIELD]
    domain_geometry = domain_row.geometry

    # Intersect only the 1992 road with the current domain.
    #
    # This determines which portion of the old road belongs to this domain.
    road_1992_in_domain = road_1992_geometry.intersection(
        domain_geometry
    )

    if road_1992_in_domain.is_empty:
        print(
            f"Domain {domain_number} | "
            "1992 road does not intersect this domain."
        )
        continue

    # Generate points along the 1992 road inside this domain
    sample_points = sample_line_geometry(
        geometry=road_1992_in_domain,
        spacing=SAMPLE_SPACING_M,
    )

    if not sample_points:
        print(
            f"Domain {domain_number} | "
            "no sample points could be generated."
        )
        continue

    # Measure each old-road point to the complete 1996 road
    relocation_distances = np.array(
        [
            point.distance(road_1996_geometry)
            for point in sample_points
        ],
        dtype=float,
    )

    mean_distance = np.mean(relocation_distances)
    median_distance = np.median(relocation_distances)
    minimum_distance = np.min(relocation_distances)
    maximum_distance = np.max(relocation_distances)
    standard_deviation = np.std(relocation_distances)

    results.append(
        {
            "domain": domain_number,
            "number_of_samples": len(relocation_distances),
            "road_1992_length_m": road_1992_in_domain.length,
            "mean_relocation_m": mean_distance,
            "median_relocation_m": median_distance,
            "minimum_relocation_m": minimum_distance,
            "maximum_relocation_m": maximum_distance,
            "std_relocation_m": standard_deviation,
        }
    )

    for sample_number, (point, distance) in enumerate(
        zip(sample_points, relocation_distances),
        start=1,
    ):
        sample_records.append(
            {
                "domain": domain_number,
                "sample_number": sample_number,
                "distance_m": distance,
                "geometry": point,
            }
        )

    print(
        f"Domain {domain_number} | "
        f"samples={len(relocation_distances)} | "
        f"road length={road_1992_in_domain.length:.2f} m | "
        f"mean={mean_distance:.2f} m | "
        f"median={median_distance:.2f} m | "
        f"min={minimum_distance:.2f} m | "
        f"max={maximum_distance:.2f} m"
    )


# =============================================================================
# SAVE SUMMARY TABLE
# =============================================================================

results_df = pd.DataFrame(results)

if results_df.empty:
    raise ValueError(
        "No distances were calculated. Check the road and "
        "domain locations and their coordinate systems."
    )

results_df = (
    results_df
    .sort_values("domain")
    .reset_index(drop=True)
)

results_df.to_csv(
    OUTPUT_CSV,
    index=False,
)

print("\n" + "=" * 70)
print("FINAL RESULTS")
print("=" * 70)

print(results_df.round(2).to_string(index=False))

print(f"\nSaved CSV:\n{OUTPUT_CSV}")


# =============================================================================
# SAVE SAMPLE POINTS
# =============================================================================

sample_gdf = gpd.GeoDataFrame(
    sample_records,
    geometry="geometry",
    crs=TARGET_CRS,
)

sample_gdf.to_file(
    OUTPUT_SAMPLE_POINTS,
    driver="GeoJSON",
)

print(f"\nSaved sample points:\n{OUTPUT_SAMPLE_POINTS}")


# =============================================================================
# PLOT RESULTS
# =============================================================================

fig, axes = plt.subplots(
    nrows=1,
    ncols=2,
    figsize=(18, 9),
    constrained_layout=True,
)


# -----------------------------------------------------------------------------
# Panel 1: Map
# -----------------------------------------------------------------------------

selected_domains.boundary.plot(
    ax=axes[0],
    color="black",
    linewidth=1,
    zorder=1,
)

road_1992.plot(
    ax=axes[0],
    color="blue",
    linewidth=2.5,
    label="1992 road",
    zorder=2,
)

road_1996.plot(
    ax=axes[0],
    color="red",
    linewidth=2.5,
    label="1996 road",
    zorder=3,
)

sample_gdf.plot(
    ax=axes[0],
    column="distance_m",
    cmap="viridis",
    markersize=18,
    legend=True,
    legend_kwds={
        "label": "Distance from 1992 to 1996 road (m)",
        "shrink": 0.7,
    },
    zorder=4,
)

# Add domain labels
for _, domain_row in selected_domains.iterrows():

    label_point = domain_row.geometry.representative_point()

    axes[0].annotate(
        text=str(domain_row[DOMAIN_ID_FIELD]),
        xy=(label_point.x, label_point.y),
        ha="center",
        va="center",
        fontsize=9,
        fontweight="bold",
        bbox={
            "boxstyle": "round,pad=0.2",
            "facecolor": "white",
            "edgecolor": "black",
            "alpha": 0.8,
        },
        zorder=5,
    )

axes[0].set_title(
    "1992–1996 Road Relocation by Domain",
    fontsize=14,
)

axes[0].set_xlabel("Easting (m)")
axes[0].set_ylabel("Northing (m)")
axes[0].set_aspect("equal")
axes[0].legend()


# -----------------------------------------------------------------------------
# Panel 2: Mean relocation distance
# -----------------------------------------------------------------------------

x_positions = np.arange(len(results_df))

axes[1].bar(
    x_positions,
    results_df["mean_relocation_m"],
    color="darkorange",
    edgecolor="black",
)

axes[1].set_xticks(x_positions)

axes[1].set_xticklabels(
    results_df["domain"].astype(str),
    rotation=90,
)

axes[1].set_title(
    "Mean Road Relocation Distance by Domain",
    fontsize=14,
)

axes[1].set_xlabel("CASCADE domain")
axes[1].set_ylabel("Mean relocation distance (m)")
axes[1].grid(axis="y", alpha=0.3)

for position, value in zip(
    x_positions,
    results_df["mean_relocation_m"],
):
    axes[1].text(
        position,
        value,
        f"{value:.1f}",
        ha="center",
        va="bottom",
        fontsize=8,
        rotation=90,
    )


# =============================================================================
# SAVE FIGURE
# =============================================================================

plt.savefig(
    OUTPUT_FIGURE,
    dpi=300,
    bbox_inches="tight",
)

plt.show()

print(f"\nSaved figure:\n{OUTPUT_FIGURE}")