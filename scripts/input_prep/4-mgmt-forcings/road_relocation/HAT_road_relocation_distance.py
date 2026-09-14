r"""
HAT_road_relocation_distance.py
===============================================================================
How far NC-12 moved between two digitised vintages, per Barrier3D domain.

Hatteras port of roya_files/road_relocation_dis.py. Same measurement: sample
the OLD road inside each domain, measure each sample point to the WHOLE new
road, and summarise per domain. Everything Hatteras-specific -- the road
vintages, the 90-polygon domain file, the CRS chain -- is in CONFIG below.

    python HAT_road_relocation_distance.py
    HAT_RELOC_FROM=1984 HAT_RELOC_TO=2004 python HAT_road_relocation_distance.py

WHAT THE NUMBER IS
------------------
`mean_relocation_m` is an UNSIGNED nearest-distance from the old road to the
new one. It says how far the road moved, not which way, and because it takes
the *nearest* point on the new road it is a lower bound on cross-shore
movement: where the two roads cross obliquely the nearest point is diagonally
alongshore, not straight across.

So this script also reports a SIGNED column, `mean_signed_landward_m`, built
from the same displacement vectors:

  landward = the new road is LEFT of the old road heading south -> north
           = west, i.e. away from the ocean         (see OCEAN_ON_RIGHT)

Positive is landward, negative is seaward. Read the signed column when you want
direction and the unsigned one when you want magnitude; `sign_agreement` tells
you what fraction of a domain's samples agreed on the direction, and anything
below ~0.9 means the two roads cross inside that domain and the domain mean is
averaging two directions.

That convention has one place it breaks. Around the Cape Point bend the island
turns east-west, the ocean is no longer on the right of a northward road, and
the sign becomes meaningless -- so domains where the road runs more than 60
degrees off north are flagged OBLIQUE_SIGN. On the 1984-2004 pair that is
domain 8 and nothing else. Magnitude is unaffected: it never used the tangent.

READ THIS BEFORE READING THE ZEROS
----------------------------------
The two Hatteras road files SHARE MOST OF THEIR GEOMETRY. 558 of the 1984
line's 791 vertices are identical to a 2004 vertex to the millimetre: the 2004
line was digitised by editing a copy of the 1984 one, and only the stretches
that visibly moved were re-drawn. Roughly 72% of the old line is therefore
exactly 0.000 m from the new line by construction.

A 0.00 m domain in this table means NOBODY EDITED THAT STRETCH. It is not
evidence the road held still. Real movement can only be claimed where the
lines actually diverge, so the script measures the shared fraction itself,
reports it per domain as `coincident_fraction`, and flags those domains
NO_EDIT rather than letting them read as a measurement.

A second artefact sits behind the first: some stretches WERE re-drawn, but
only re-traced -- the new centreline wanders a metre or two and never leaves
the old road's own footprint. That is one road digitised twice, not a road
that moved, so domains whose largest displacement anywhere stays under
REDIGITIZE_MAX_M are flagged REDIGITIZED.

Every domain therefore carries a `classification`:

    no_edit      the line was copied through unedited
    redigitized  re-traced within the road's own width
    relocated    the road actually moved  <- the only measured subset

Nothing is dropped from the CSV, but only `relocated` domains are summarised
and drawn. Any island-wide average over the full table is meaningless -- it
averages real movement against copy and re-tracing alike.

WHAT THIS IS, AND IS NOT
------------------------
Not a setback. `road_setback` comes from ../road_offset/, measured against
interior row 0 of the model grid; this is a GIS-frame observation of how far
apart two digitised lines lie.

It IS a CASCADE forcing as of 2026-08-20. `HATTERAS_ROAD_EVENTS` in
scripts/hatteras_site_config.py reads `mean_signed_landward_m` out of the
1984_2004 CSV for the domains its two historical relocation events move -- GIS
84-87 (1989) and GIS 9-15 (1999). It replaced eleven hand-entered literals
attributed to a 1978->1997 ArcGIS measurement whose 1997 line is not in the
repo. So re-running this script with those vintages CHANGES WHAT THE MODEL IS
FORCED WITH; it is no longer a diagnostic that can be regenerated freely. The
config refuses any domain this script classifies `no_edit` or `redigitized`,
so a re-run that reclassifies one of those eleven raises at import rather than
quietly forcing a number the lines do not support.

A NOTE ON THE VINTAGES
----------------------
nc12_1984.geojson and nc12_2004.geojson were digitised off 1978 and 2008
imagery -- the nearest usable coverage to the two hindcast period starts. The
labels are the periods they stand in for, not the photo dates, so the interval
measured here is really ~30 years, not 20.

REQUIREMENTS
------------
  geopandas, shapely, numpy, pandas, matplotlib
===============================================================================
"""

import os
import sys
from pathlib import Path

import geopandas as gpd
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from shapely.geometry import (
    GeometryCollection,
    LineString,
    MultiLineString,
    Point,
)
from shapely.affinity import rotate as shapely_rotate
from shapely.ops import nearest_points, unary_union


# =============================================================================
# CONFIG
# =============================================================================

PROJECT_ROOT = next(_p for _p in Path(__file__).resolve().parents
                    if (_p / "pyproject.toml").exists())
DATA_DIR = PROJECT_ROOT / "data" / "hatteras_init"

# Which pair of road vintages to compare. Override from the shell to compare a
# different pair without editing the file.
YEAR_FROM = int(os.environ.get("HAT_RELOC_FROM", 1984))
YEAR_TO = int(os.environ.get("HAT_RELOC_TO", 2004))


def road_file(year):
    """The digitised NC-12 centreline for one vintage."""

    return (
        DATA_DIR
        / "4-mgmt-forcing"
        / "road_offset"
        / "raw_offset"
        / str(year)
        / f"nc12_{year}.geojson"
    )


ROAD_FROM_FILE = road_file(YEAR_FROM)
ROAD_TO_FILE = road_file(YEAR_TO)

# The same 90-polygon domain file the shoreline and groin work uses.
DOMAIN_FILE = (
    PROJECT_ROOT
    / "scripts"
    / "input_prep"
    / "5-scr"
    / "CoastSat"
    / "HAT_domains.json"
)

OUTPUT_DIR = (
    DATA_DIR
    / "4-mgmt-forcing"
    / "road_relocation"
    / f"{YEAR_FROM}_{YEAR_TO}"
)

OUTPUT_CSV = OUTPUT_DIR / f"road_relocation_{YEAR_FROM}_{YEAR_TO}.csv"
OUTPUT_SAMPLE_POINTS = (
    OUTPUT_DIR / f"road_relocation_{YEAR_FROM}_{YEAR_TO}_sample_points.geojson"
)
# Three figures, each answering one question, rather than one page trying to
# answer all three at incompatible scales.
_FIGURE_STEM = f"road_relocation_{YEAR_FROM}_{YEAR_TO}"

OUTPUT_FIGURE_ALONGSHORE = OUTPUT_DIR / f"{_FIGURE_STEM}_alongshore.png"
OUTPUT_FIGURE_SITES = OUTPUT_DIR / f"{_FIGURE_STEM}_sites.png"
OUTPUT_FIGURE_DOMAIN_MAP = OUTPUT_DIR / f"{_FIGURE_STEM}_domain_map.png"

# The field carrying the domain number in the DOMAIN file (not the road file).
DOMAIN_ID_FIELD = "domain_id"

# Domains to analyse. None means every domain the old road touches.
# Barrier3D runs GIS 9-90; the road file reaches 8-90.
DOMAIN_NUMBERS = None

# Spacing of the points sampled along the old road, metres.
SAMPLE_SPACING_M = 5.0

# NAD83(2011) / UTM Zone 18N. The road files are EPSG:2264 (NC State Plane,
# US survey feet) and the domains EPSG:3725, so everything is reprojected here
# and every distance below is metres.
TARGET_CRS = "EPSG:6347"

# NC-12 runs south -> north up Hatteras with the Atlantic to the east, so the
# ocean is on the RIGHT of the direction of travel and landward is LEFT. This
# is what turns an unsigned distance into a signed one; flip it for a site
# where the ocean sits on the other hand.
OCEAN_ON_RIGHT = True

# Town and village spans, used only to name the relocation sites on the
# figure. Taken from hatteras_site_config rather than restated here, so a
# domain number means the same place in this figure as in every other one.
sys.path.insert(0, str(PROJECT_ROOT / "scripts"))

try:
    from hatteras_site_config import HATTERAS_ANNOTATIONS, HATTERAS_ROAD_EVENTS

    TOWN_SPANS = HATTERAS_ANNOTATIONS.town_spans
    VILLAGE_LINES = HATTERAS_ANNOTATIONS.village_lines
    # The domains the model's two historical events actually move. The site
    # figure labels each of these with the displacement it is forced with, so
    # a reader sees the number that reaches CASCADE, not only the colour; the
    # domains a site spans but the events leave alone (GIS 8, 15) say so.
    PRESCRIBED_DOMAINS = {
        gis
        for event in HATTERAS_ROAD_EVENTS
        for gis in getattr(event, "displacement_m", {}) or {}
    }
    # The forcing is the measurement rounded to the nearest cell (see
    # hatteras_site_config, ROUNDED TO WHOLE CELLS); label both.
    from hatteras_site_config import round_to_cell as _round_to_cell

except ImportError:
    # The measurement does not need them; only the labels do.
    print("Note: hatteras_site_config unavailable -- sites go unnamed.")
    TOWN_SPANS = {}
    VILLAGE_LINES = {}
    PRESCRIBED_DOMAINS = set()
    _round_to_cell = None

# Below this fraction of samples agreeing on direction, the two roads cross
# inside the domain and the signed mean is mixing landward with seaward.
SIGN_AGREEMENT_FLAG = 0.90

# A sample closer than this to the new road is sitting on SHARED geometry --
# a vertex the digitiser never edited -- not on a measured non-movement. Well
# below any plausible digitising precision, so nothing real is discarded.
# These samples carry no direction and are excluded from the signed statistics.
COINCIDENT_TOLERANCE_M = 0.5

# Above this coincident fraction, the domain is unedited copy rather than a
# measurement, and is flagged NO_EDIT.
NO_EDIT_FLAG = 0.90

# A domain whose LARGEST displacement anywhere still falls below this was
# re-traced, not relocated: the new centreline never leaves the old road's own
# footprint. NC-12 is ~8 m wide, so 5 m is about half a road width -- two
# digitisings of one road off different photos disagree by that much from
# georeferencing alone.
#
# The Hatteras data splits cleanly here and the exact value does not matter:
# the re-traced domains top out at 3.3 m and the real relocations start at
# 12.2 m, so anything from ~3.5 to ~12 m gives the same answer. Corroborated
# by direction -- re-traced domains are internally coherent but flip sign
# arbitrarily between neighbours (16-18 seaward, 19-22 landward, 24-27
# seaward), which is what a smoothed line does and a moved road does not.
REDIGITIZE_MAX_M = 5.0

# Northward component of the road's local tangent below which the landward /
# seaward convention stops being safe: 0.5 is a bearing more than 60 degrees
# off north. Around the Cape Point bend the island turns east-west and the
# ocean is no longer on the right of a northward road, so those domains are
# flagged OBLIQUE_SIGN and their sign should not be trusted. Magnitude is
# unaffected -- it never depended on the tangent.
OBLIQUE_TANGENT = 0.5


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

    Returns (point, tangent) pairs. The tangent is the local direction of the
    sampled line, oriented northward, and is what gives the distance a sign.
    """

    samples = []

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
            samples.append(
                (
                    line.interpolate(distance),
                    local_tangent(line, distance),
                )
            )

    return samples


def local_tangent(line, distance, step=1.0):
    """
    Unit direction of `line` at `distance` along it, oriented so it points
    north. Returns None where the tangent cannot be resolved.
    """

    back = line.interpolate(max(0.0, distance - step))
    forward = line.interpolate(min(line.length, distance + step))

    dx = forward.x - back.x
    dy = forward.y - back.y

    magnitude = np.hypot(dx, dy)

    if magnitude == 0:
        return None

    dx /= magnitude
    dy /= magnitude

    # Orient northward so "left" means the same thing everywhere, whichever
    # way the digitised line happens to run.
    if dy < 0:
        dx, dy = -dx, -dy

    return dx, dy


def signed_relocation(point, tangent, target_geometry):
    """
    Distance from `point` to `target_geometry`, signed positive landward.

    The sign is the side of the old road the new road falls on: the z of the
    cross product of the northward tangent with the displacement vector is
    positive to the left, and left is landward when the ocean is on the right.
    """

    _, nearest = nearest_points(point, target_geometry)

    dx = nearest.x - point.x
    dy = nearest.y - point.y

    distance = float(np.hypot(dx, dy))

    if tangent is None or distance == 0:
        return distance, np.nan, nearest

    cross = tangent[0] * dy - tangent[1] * dx

    landward = cross > 0 if OCEAN_ON_RIGHT else cross < 0

    return distance, distance if landward else -distance, nearest


# =============================================================================
# CREATE OUTPUT DIRECTORIES
# =============================================================================

OUTPUT_DIR.mkdir(parents=True, exist_ok=True)


# =============================================================================
# READ INPUT FILES
# =============================================================================

for label, path in (
    (f"{YEAR_FROM} road", ROAD_FROM_FILE),
    (f"{YEAR_TO} road", ROAD_TO_FILE),
    ("domains", DOMAIN_FILE),
):
    if not path.exists():
        raise FileNotFoundError(f"{label} not found:\n{path}")

road_from = gpd.read_file(ROAD_FROM_FILE)
road_to = gpd.read_file(ROAD_TO_FILE)
domains = gpd.read_file(DOMAIN_FILE)

print("=" * 70)
print("INPUT FILE INFORMATION")
print("=" * 70)

print(f"\n{YEAR_FROM} road: {ROAD_FROM_FILE}")
print(f"{YEAR_TO} road: {ROAD_TO_FILE}")
print(f"Domains:   {DOMAIN_FILE}")

print(f"\n{YEAR_FROM} road columns:")
print(road_from.columns.tolist())

print(f"\n{YEAR_TO} road columns:")
print(road_to.columns.tolist())

print("\nDomain polygon columns:")
print(domains.columns.tolist())

print("\nOriginal coordinate systems:")
print(f"{YEAR_FROM} road: {road_from.crs}")
print(f"{YEAR_TO} road: {road_to.crs}")
print(f"Domains:   {domains.crs}")

print("\nNumber of features:")
print(f"{YEAR_FROM} road features: {len(road_from)}")
print(f"{YEAR_TO} road features: {len(road_to)}")
print(f"Domain features:    {len(domains)}")


# =============================================================================
# CHECK INPUT DATA
# =============================================================================

if road_from.empty:
    raise ValueError(f"The {YEAR_FROM} road file is empty.")

if road_to.empty:
    raise ValueError(f"The {YEAR_TO} road file is empty.")

if domains.empty:
    raise ValueError("The domain polygon file is empty.")

if road_from.crs is None:
    raise ValueError(f"The {YEAR_FROM} road file does not have a defined CRS.")

if road_to.crs is None:
    raise ValueError(f"The {YEAR_TO} road file does not have a defined CRS.")

if domains.crs is None:
    raise ValueError("The domain file does not have a defined CRS.")

if DOMAIN_ID_FIELD not in domains.columns:
    raise ValueError(
        f"\nThe field '{DOMAIN_ID_FIELD}' was not found in "
        f"the DOMAIN file.\n\n"
        f"Available domain fields are:\n"
        f"{domains.columns.tolist()}\n\n"
        f"Change DOMAIN_ID_FIELD to the field carrying GIS domain numbers."
    )


# =============================================================================
# REPROJECT EVERYTHING TO THE SAME CRS
# =============================================================================

road_from = road_from.to_crs(TARGET_CRS)
road_to = road_to.to_crs(TARGET_CRS)
domains = domains.to_crs(TARGET_CRS)

print("\nProjected coordinate systems:")
print(f"{YEAR_FROM} road: {road_from.crs}")
print(f"{YEAR_TO} road: {road_to.crs}")
print(f"Domains:   {domains.crs}")

print("\nDistances will be calculated in metres.")


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
print(f"{available_domains[0]}-{available_domains[-1]} "
      f"({len(available_domains)} polygons)")


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

selected_domains = selected_domains.sort_values(DOMAIN_ID_FIELD)


# =============================================================================
# COMBINE ROAD FEATURES
# =============================================================================

# Each vintage may be one feature or many; combine into a single geometry so
# the nearest-point search sees the whole road.
road_from_geometry = unary_union(
    road_from.geometry.dropna().tolist()
)

road_to_geometry = unary_union(
    road_to.geometry.dropna().tolist()
)

if road_from_geometry.is_empty:
    raise ValueError(f"The combined {YEAR_FROM} road geometry is empty.")

if road_to_geometry.is_empty:
    raise ValueError(f"The combined {YEAR_TO} road geometry is empty.")


# =============================================================================
# HOW MUCH GEOMETRY DO THE TWO VINTAGES SHARE?
# =============================================================================
# Run before anything is measured, because it decides how the whole table
# should be read. Two lines that share vertices were not digitised
# independently, and every shared vertex contributes a 0.000 m "relocation"
# that is an artefact of the editing workflow.

shared_vertex_fraction = np.nan

from_vertices = [
    np.asarray(line.coords) for line in extract_lines(road_from_geometry)
]
to_vertices = [
    np.asarray(line.coords) for line in extract_lines(road_to_geometry)
]

if from_vertices and to_vertices:

    from_stack = np.round(np.vstack(from_vertices)[:, :2], 3)
    to_stack = np.round(np.vstack(to_vertices)[:, :2], 3)

    to_set = set(map(tuple, to_stack))
    shared = sum(1 for xy in map(tuple, from_stack) if xy in to_set)

    shared_vertex_fraction = shared / len(from_stack)

    print("\n" + "=" * 70)
    print("SHARED GEOMETRY CHECK")
    print("=" * 70)
    print(
        f"{YEAR_FROM} vertices: {len(from_stack)}\n"
        f"{YEAR_TO} vertices: {len(to_stack)}\n"
        f"Identical to the millimetre: {shared} "
        f"({shared_vertex_fraction:.1%} of the {YEAR_FROM} line)"
    )

    if shared_vertex_fraction > 0.05:
        print(
            f"\n  WARNING: the two vintages are not independent digitisations.\n"
            f"  {shared_vertex_fraction:.0%} of the {YEAR_FROM} line was carried\n"
            f"  into {YEAR_TO} unedited, and every one of those vertices\n"
            f"  measures 0.000 m by construction. Domains built mostly from\n"
            f"  them are flagged NO_EDIT below. Do NOT read those zeros as\n"
            f"  'the road did not move', and do not average across the table."
        )


# =============================================================================
# CALCULATE RELOCATION DISTANCE BY DOMAIN
# =============================================================================

results = []
sample_records = []
skipped_domains = []

print("\n" + "=" * 70)
print(f"ROAD RELOCATION RESULTS  {YEAR_FROM} -> {YEAR_TO}")
print("=" * 70)

for _, domain_row in selected_domains.iterrows():

    domain_number = domain_row[DOMAIN_ID_FIELD]
    domain_geometry = domain_row.geometry

    # Intersect only the OLD road with the current domain. This determines
    # which portion of the old road belongs to this domain.
    road_from_in_domain = road_from_geometry.intersection(
        domain_geometry
    )

    if road_from_in_domain.is_empty:
        skipped_domains.append(domain_number)
        continue

    # Generate points along the old road inside this domain
    samples = sample_line_geometry(
        geometry=road_from_in_domain,
        spacing=SAMPLE_SPACING_M,
    )

    if not samples:
        print(
            f"Domain {domain_number} | "
            "no sample points could be generated."
        )
        skipped_domains.append(domain_number)
        continue

    # Measure each old-road point to the complete new road
    measured = [
        signed_relocation(point, tangent, road_to_geometry)
        for point, tangent in samples
    ]

    relocation_distances = np.array(
        [item[0] for item in measured],
        dtype=float,
    )

    signed_distances = np.array(
        [item[1] for item in measured],
        dtype=float,
    )

    mean_distance = np.mean(relocation_distances)
    median_distance = np.median(relocation_distances)
    minimum_distance = np.min(relocation_distances)
    maximum_distance = np.max(relocation_distances)
    standard_deviation = np.std(relocation_distances)

    # Samples sitting on shared, unedited geometry carry no direction: their
    # sign is numerical noise on a sub-millimetre displacement. Keep them out
    # of every signed statistic, and count them instead.
    moved = relocation_distances >= COINCIDENT_TOLERANCE_M

    coincident_fraction = float(np.mean(~moved))

    # Where the road runs east-west rather than north-south -- around the
    # Cape Point bend -- "ocean on the right of northward travel" stops
    # holding, and the sign of the displacement cannot be trusted. The
    # northward component of the local tangent is what detects it.
    northward = np.array(
        [abs(tangent[1]) if tangent else np.nan for _, tangent in samples],
        dtype=float,
    )

    with np.errstate(invalid="ignore"):
        oblique_fraction = float(np.nanmean(northward < OBLIQUE_TANGENT))

    finite_signed = signed_distances[moved & np.isfinite(signed_distances)]

    if finite_signed.size:
        mean_signed = float(np.mean(finite_signed))
        median_signed = float(np.median(finite_signed))
        landward_fraction = float(np.mean(finite_signed > 0))
        sign_agreement = max(landward_fraction, 1.0 - landward_fraction)
    else:
        mean_signed = np.nan
        median_signed = np.nan
        landward_fraction = np.nan
        sign_agreement = np.nan

    # Three-way classification. Only `relocated` domains carry a measurement
    # of the road moving; the other two are artefacts of how the lines were
    # drawn, and are kept in the table but out of every summary.
    if coincident_fraction >= NO_EDIT_FLAG:
        classification = "no_edit"

    elif maximum_distance < REDIGITIZE_MAX_M:
        classification = "redigitized"

    else:
        classification = "relocated"

    flags = []

    if classification == "no_edit":
        flags.append("NO_EDIT")

    elif classification == "redigitized":
        flags.append("REDIGITIZED")

    elif np.isfinite(sign_agreement) and sign_agreement < SIGN_AGREEMENT_FLAG:
        # Only meaningful once there is real movement to have a direction.
        flags.append("ROADS_CROSS")

    if oblique_fraction > 0.5:
        flags.append("OBLIQUE_SIGN")

    if road_from_in_domain.length < SAMPLE_SPACING_M:
        flags.append("SHORT_SEGMENT")

    results.append(
        {
            "domain": domain_number,
            "number_of_samples": len(relocation_distances),
            f"road_{YEAR_FROM}_length_m": road_from_in_domain.length,
            "mean_relocation_m": mean_distance,
            "median_relocation_m": median_distance,
            "minimum_relocation_m": minimum_distance,
            "maximum_relocation_m": maximum_distance,
            "std_relocation_m": standard_deviation,
            "mean_signed_landward_m": mean_signed,
            "median_signed_landward_m": median_signed,
            "landward_fraction": landward_fraction,
            "sign_agreement": sign_agreement,
            "coincident_fraction": coincident_fraction,
            "oblique_fraction": oblique_fraction,
            "classification": classification,
            "flags": ";".join(flags),
        }
    )

    for sample_number, ((point, _tangent), (distance, signed, _nearest)) in (
        enumerate(zip(samples, measured), start=1)
    ):
        sample_records.append(
            {
                "domain": domain_number,
                "sample_number": sample_number,
                "distance_m": distance,
                "signed_landward_m": (
                    signed
                    if distance >= COINCIDENT_TOLERANCE_M
                    else np.nan
                ),
                "coincident": bool(distance < COINCIDENT_TOLERANCE_M),
                "classification": classification,
                "geometry": point,
            }
        )

    print(
        f"Domain {domain_number:>3} | "
        f"samples={len(relocation_distances):>4} | "
        f"road length={road_from_in_domain.length:7.2f} m | "
        f"mean={mean_distance:7.2f} m | "
        f"median={median_distance:7.2f} m | "
        f"max={maximum_distance:7.2f} m | "
        f"edited={1.0 - coincident_fraction:5.1%} | "
        + (
            f"signed={mean_signed:8.2f} m"
            if np.isfinite(mean_signed)
            else "signed=      -- "
        )
        + (f" | {';'.join(flags)}" if flags else "")
    )


if skipped_domains:
    print(
        f"\n{len(skipped_domains)} domain(s) carry no {YEAR_FROM} road "
        f"and were skipped:\n{skipped_domains}"
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

no_edit = results_df[results_df["classification"] == "no_edit"]
redigitized = results_df[results_df["classification"] == "redigitized"]
relocated = results_df[results_df["classification"] == "relocated"]
crossing = relocated[relocated["flags"].str.contains("ROADS_CROSS")]

print("\n" + "-" * 70)
print(
    f"Domains carrying the {YEAR_FROM} road: {len(results_df)}\n"
    f"  never edited (0 m by construction): {len(no_edit)}\n"
    f"  re-traced, not relocated:           {len(redigitized)}  "
    f"{redigitized['domain'].tolist()}\n"
    f"  RELOCATED:                          {len(relocated)}  "
    f"{relocated['domain'].tolist()}"
)

if not redigitized.empty:
    print(
        f"\nRe-traced domains: largest displacement anywhere is "
        f"{redigitized['maximum_relocation_m'].max():.2f} m, below the "
        f"{REDIGITIZE_MAX_M:.1f} m\n"
        f"cut, so the new centreline never leaves the old road's footprint.\n"
        f"Excluded from the numbers below."
    )

if not relocated.empty:
    print(
        f"\nAcross the {len(relocated)} RELOCATED domains only:\n"
        f"  mean relocation:       "
        f"{relocated['mean_relocation_m'].mean():.2f} m\n"
        f"  median relocation:     "
        f"{relocated['median_relocation_m'].median():.2f} m\n"
        f"  largest domain mean:   "
        f"{relocated['mean_relocation_m'].max():.2f} m "
        f"(domain "
        f"{relocated.loc[relocated['mean_relocation_m'].idxmax(), 'domain']})\n"
        f"  mostly landward:       "
        f"{int((relocated['mean_signed_landward_m'] > 0).sum())} domains\n"
        f"  mostly seaward:        "
        f"{int((relocated['mean_signed_landward_m'] < 0).sum())} domains\n"
        f"  roads cross inside:    "
        f"{len(crossing)} domains {crossing['domain'].tolist()}"
    )

print(
    "\nNo island-wide mean is reported: it would average real relocation\n"
    "against unedited copy and re-traced line. Every domain is still in the\n"
    "CSV, labelled in the `classification` column -- nothing was dropped."
)
print("-" * 70)

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
# Three figures, because they work at three incompatible scales and one page
# cannot serve all of them:
#
#   1. alongshore   WHERE along the island the road moved      (domain axis)
#   2. sites        WHAT the move looked like                  (~2 km, true scale)
#   3. domain map   WHICH domains carry a relocation           (45 km, true scale)
#
# Hatteras is 8 km wide and 45 km long -- aspect 0.18 -- so any true-scale map
# of the whole island is a hairline in a column of white space. Figure 3 gets
# round that by rotating the island to run left-right; figure 2 sidesteps it by
# only ever drawing ~2 km at a time.
#
# HOUSE STYLE (Hannah, 2026-09-10: "more professional and academic styled").
# The dune-line figures' style block is reused, not restated: Arial 8-10 pt,
# thin dark-grey axes, the ColorBrewer RdBu poles (earlier vintage red, later
# blue), panel letters, north arrow and scale bar on tickless maps, frameless
# legends outside the axes, and NO title sentences, statistics lines or
# footnotes on the canvas. That text is in CAPTIONS.md beside the figures,
# written at the end of this section with the numbers filled from the table,
# so it cannot go stale against the picture without the CSV going stale too.

from hat_figure_style import (  # noqa: E402   scripts/ is on sys.path since the site-config import
    DOMAIN_AXIS_LABEL, INK, INK_MUTED, C_1984 as C_EARLY, C_1997 as C_LATE,
    C_1984_FILL as C_EARLY_FILL, C_1997_FILL as C_LATE_FILL,
    FIG_H_MAX, FIG_W_DOUBLE, apply_style, figsize, save,
    _north_arrow, _title as _panel_title, _halo, open_frame, town_bands,
)
apply_style()
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

CELL_M = 10.0

# Legend wording: no working vocabulary (Hannah, 2026-09-08). A reader outside
# the project has to be able to take each entry literally.
LABEL_NO_ROAD = f"no {YEAR_FROM} road in domain"
LABEL_NO_EDIT = "centreline unchanged between surveys"
LABEL_REDIG = f"centreline re-digitised, no relocation (< {REDIGITIZE_MAX_M:.0f} m)"
LABEL_NOT_PRESCRIBED = "no prescribed move in the model"
SHADE_NO_EDIT, SHADE_REDIG = "0.93", "0.84"


def contiguous_runs(numbers, max_gap=1):
    """
    Group sorted domain numbers into runs, allowing gaps up to `max_gap` so a
    single unmeasured domain does not split one relocation into two sites.
    """

    runs = []

    for number in sorted(numbers):

        if runs and number - runs[-1][-1] <= max_gap:
            runs[-1].append(number)
        else:
            runs.append([number])

    return runs


def site_label(first, last, neighbour_range=6):
    """
    Name a run of domains after the place it overlaps.

    Falling back to the nearest village matters here: the northern relocation
    sits at domains 84-87, past the end of every span in the site config, and
    would otherwise go unnamed on the figure. It is the stretch north of
    Rodanthe -- so say that, rather than claim it IS Rodanthe.
    """

    names = [
        name
        for name, (start, end) in TOWN_SPANS.items()
        if first <= end and last >= start
    ]

    names += [
        name
        for name, position in VILLAGE_LINES.items()
        if first <= position <= last
    ]

    if names:
        return " / ".join(dict.fromkeys(names))

    if not VILLAGE_LINES:
        return ""

    nearest, position = min(
        VILLAGE_LINES.items(),
        key=lambda item: min(abs(item[1] - first), abs(item[1] - last)),
    )

    if min(abs(position - first), abs(position - last)) > neighbour_range:
        return ""

    return f"N of {nearest}" if position < first else f"S of {nearest}"


def add_scale_bar(ax, length_m=500, fraction_x=0.06, fraction_y=0.06):
    """A capped bar in data units, so it scales with the panel and cannot
    disagree with it. White halo so it reads on any backdrop."""

    x0, x1 = ax.get_xlim()
    y0, y1 = ax.get_ylim()

    x_start = x0 + fraction_x * (x1 - x0)
    y_start = y0 + fraction_y * (y1 - y0)
    tick = 0.008 * (y1 - y0)

    ax.plot([x_start, x_start + length_m], [y_start, y_start], color=INK, linewidth=1.6,
            solid_capstyle="butt", zorder=12, path_effects=_halo(3.6))
    for x_tick in (x_start, x_start + length_m):
        ax.plot([x_tick, x_tick], [y_start - tick, y_start + tick], color=INK, linewidth=1.0,
                solid_capstyle="butt", zorder=12, path_effects=_halo(2.6))
    ax.text(x_start + length_m / 2, y_start + 1.8 * tick,
            f"{length_m:,} m" if length_m < 1000 else f"{length_m / 1000:g} km",
            ha="center", va="bottom", fontsize=8, color=INK, zorder=12,
            bbox=dict(facecolor="white", alpha=0.8, edgecolor="none", boxstyle="square,pad=0.15"))


def _east_arrow(ax, x=0.955, y=0.10, length=0.03):
    """North arrow for the rotated strip map, where north is to the RIGHT."""
    ax.annotate("", xy=(x + length, y), xytext=(x, y), xycoords="axes fraction",
                textcoords="axes fraction", zorder=20,
                arrowprops=dict(arrowstyle="-|>", color=INK, lw=1.0, shrinkA=0, shrinkB=0,
                                mutation_scale=11))
    ax.text(x + length + 0.004, y, "N", transform=ax.transAxes, ha="left", va="center",
            fontsize=8.5, fontweight="bold", color=INK, zorder=20)


relocation_sites = contiguous_runs(relocated["domain"].tolist())

print("\nRelocation sites found: " + ", ".join(
    f"{run[0]}-{run[-1]}"
    + (f" ({site_label(run[0], run[-1])})"
       if site_label(run[0], run[-1]) else "")
    for run in relocation_sites
))

domain_positions = results_df["domain"].to_numpy()

is_relocated = (results_df["classification"] == "relocated").to_numpy()


def measured_only(column):
    """The column, blanked wherever the domain is not a relocation."""

    values = results_df[column].to_numpy(dtype=float).copy()
    values[~is_relocated] = np.nan

    return values


# One colour scale shared by figures 2 and 3, so a colour means the same
# distance in both. Viridis: sequential, print- and CVD-safe.
colour_norm = plt.Normalize(
    vmin=0,
    vmax=float(relocated["maximum_relocation_m"].max()),
)
COLOUR_MAP = "viridis"


# =============================================================================
# FIGURE 1 -- ALONGSHORE: where along the island the road moved
# =============================================================================

fig_alongshore, overview_ax = plt.subplots(
    figsize=figsize("double", aspect=0.40),
    constrained_layout=True,
)

# Shade the domains that carry no measurement of the road moving. Both stay
# lighter than the data; re-traced is the darker of the two because it is the
# one a reader is likelier to mistake for a small relocation.
for frame, shade in ((no_edit, SHADE_NO_EDIT), (redigitized, SHADE_REDIG)):
    for _, row in frame.iterrows():
        overview_ax.axvspan(row["domain"] - 0.5, row["domain"] + 0.5,
                            color=shade, linewidth=0, zorder=0)

signed_mean = measured_only("mean_signed_landward_m")
signed_direction = np.sign(np.nan_to_num(signed_mean))
landward = signed_mean >= 0

# Pale bar = the largest displacement anywhere in the domain, solid bar = the
# domain mean. The gap between them is how much of the domain actually moved.
# Landward carries the later vintage's blue, seaward the earlier one's red,
# so the sign reads off the colour before the axis is consulted.
overview_ax.bar(
    domain_positions,
    measured_only("maximum_relocation_m") * signed_direction,
    color=np.where(landward, C_LATE_FILL, C_EARLY_FILL),
    width=0.9,
    zorder=2,
)

overview_ax.bar(
    domain_positions,
    signed_mean,
    color=np.where(landward, C_LATE, C_EARLY),
    width=0.9,
    zorder=3,
)

overview_ax.axhline(0, color=INK, linewidth=0.8, zorder=4)

# Bracket and name each relocation site.
label_height = float(np.nanmax(measured_only("maximum_relocation_m"))) * 1.10

for run in relocation_sites:

    name = site_label(run[0], run[-1])

    overview_ax.annotate(
        "",
        xy=(run[0] - 0.5, label_height),
        xytext=(run[-1] + 0.5, label_height),
        arrowprops={"arrowstyle": "|-|", "linewidth": 1.0, "color": INK,
                    "mutation_scale": 3},
        annotation_clip=False,
    )

    overview_ax.text(
        (run[0] + run[-1]) / 2,
        label_height * 1.04,
        (f"{name}\ndomains {run[0]}–{run[-1]}" if name
         else f"domains {run[0]}–{run[-1]}"),
        ha="center",
        va="bottom",
        fontsize=8.5,
        color=INK,
    )

# The villages as the house bands behind the axis, named once, so a domain
# number means a place. `town_bands` reads the same spans from the site config.
town_bands(overview_ax, where="bottom", strip=0.07, shade="0.72")
# A STRIP, not a wash: this panel already shades two data classes
# ("centreline unchanged", "re-digitised") full height in grey, and a
# third full-height grey for the villages cannot be told from them.

# room at the left for the Buxton band's name, which sits at the very end
overview_ax.set_xlim(domain_positions.min() - 4, domain_positions.max() + 1)
overview_ax.set_ylim(-0.10 * label_height, label_height * 1.30)

overview_ax.set_yticks(
    np.arange(0, label_height, 25 if label_height > 100 else 10)
)

overview_ax.set_xlabel(DOMAIN_AXIS_LABEL)
overview_ax.set_ylabel(f"Displacement of NC-12 centreline,\n{YEAR_FROM}–{YEAR_TO} (m, landward positive)")
overview_ax.grid(axis="y", zorder=1)
open_frame(overview_ax)

alongshore_handles = [
    Patch(facecolor=C_LATE, label="mean displacement in domain, landward"),
    Patch(facecolor=C_LATE_FILL, label="largest displacement within domain"),
    Patch(facecolor=SHADE_NO_EDIT, label=LABEL_NO_EDIT),
    Patch(facecolor=SHADE_REDIG, label=LABEL_REDIG),
]
if bool(np.any(signed_mean < 0)):
    alongshore_handles.insert(1, Patch(facecolor=C_EARLY, label="mean displacement in domain, seaward"))

# two rows, so the key cannot push the figure past its printed width
fig_alongshore.legend(handles=alongshore_handles, loc="outside upper center",
                      ncol=2, frameon=False)

save(fig_alongshore, OUTPUT_FIGURE_ALONGSHORE, bbox_inches="tight")

print(f"\nSaved alongshore figure:\n{OUTPUT_FIGURE_ALONGSHORE} (+ .pdf)")


# =============================================================================
# FIGURE 2 -- SITES: true-scale zoom on each stretch that actually moved
# =============================================================================

# ONE SCALE ACROSS ALL PANELS. Buxton spans 8 domains and Rodanthe 4, so
# framing each on its own extent renders them at different metres-per-inch and
# the shorter site looks like the bigger relocation. Every panel gets the same
# vertical span (the tallest site, padded); its width follows its own site, so
# a centimetre means the same distance in both and no panel is mostly white.
site_windows = [
    domains[domains[DOMAIN_ID_FIELD].isin(run)].total_bounds
    for run in relocation_sites
]

SITE_PAD = 1.10
common_height = SITE_PAD * max(b[3] - b[1] for b in site_windows)
panel_widths = [max(SITE_PAD * (b[2] - b[0]), 0.55 * common_height) for b in site_windows]

# Drawn at the printed width: the panels share what the colour bar and the
# margins leave of a double column, and their common height follows.
SITE_CHROME_IN = 1.6                     # colour bar, labels and margins
panel_height_in = min(
    FIG_H_MAX - 0.9,
    (FIG_W_DOUBLE - SITE_CHROME_IN) / (sum(panel_widths) / common_height),
)

fig_sites, site_axes = plt.subplots(
    nrows=1,
    ncols=max(len(relocation_sites), 1),
    figsize=figsize("double", height=panel_height_in + 0.9),
    gridspec_kw={"width_ratios": panel_widths},
    constrained_layout=True,
)

site_axes = np.atleast_1d(site_axes)

site_signed = relocated.set_index("domain")["mean_signed_landward_m"]

for column, run in enumerate(relocation_sites):

    ax = site_axes[column]

    # One domain of context on each side, so the road is seen entering and
    # leaving the relocation rather than starting at the panel edge.
    context = domains[
        domains[DOMAIN_ID_FIELD].between(run[0] - 1, run[-1] + 1)
    ]

    context.plot(ax=ax, facecolor="0.96", edgecolor="0.75", linewidth=0.5, zorder=1)

    site_domains = domains[domains[DOMAIN_ID_FIELD].isin(run)]

    site_domains.plot(ax=ax, facecolor="none", edgecolor=INK, linewidth=0.8, zorder=2)

    road_to.plot(ax=ax, color=C_LATE, linewidth=2.0, zorder=4)

    site_samples = sample_gdf[sample_gdf["domain"].isin(run)]

    unmoved_here = site_samples[site_samples["coincident"]]

    if not unmoved_here.empty:
        unmoved_here.plot(ax=ax, color="0.6", markersize=3, zorder=5)

    moved_here = site_samples[~site_samples["coincident"]]

    if not moved_here.empty:
        moved_here.plot(ax=ax, column="distance_m", cmap=COLOUR_MAP, norm=colour_norm,
                        markersize=22, zorder=6)

    # The old road goes ON TOP of its own sample points, as a thin dashed
    # spine. Drawn underneath it is invisible -- the samples sit exactly on it
    # -- and the reader loses the one thing this panel exists to show: the two
    # alignments pulling apart.
    road_from.plot(ax=ax, color=C_EARLY, linewidth=0.9, linestyle=(0, (4, 3)), zorder=7)

    # Same vertical span on every site, centred on this one; width its own.
    minx, miny, maxx, maxy = site_windows[column]

    centre_x = (minx + maxx) / 2
    centre_y = (miny + maxy) / 2

    ax.set_xlim(centre_x - panel_widths[column] / 2, centre_x + panel_widths[column] / 2)
    ax.set_ylim(centre_y - common_height / 2, centre_y + common_height / 2)
    ax.set_aspect("equal")

    # Domain number and, under it, the measured displacement beside the value
    # CASCADE is forced with (that mean rounded to the nearest cell; see
    # hatteras_site_config, ROUNDED TO WHOLE CELLS). Anchored low in the left
    # of each box: the road runs up the right side of every site box, and
    # crosses mid-height only in the Buxton loop (GIS 8), so the lower-left
    # corner is clear in all of them. A domain the site spans but no event
    # moves says so, rather than being read as a prescribed move from its
    # colour.
    for _, domain_row in site_domains.iterrows():

        gis = int(domain_row[DOMAIN_ID_FIELD])
        bminx, bminy, bmaxx, bmaxy = domain_row.geometry.bounds
        anchor = (bminx + 0.03 * (bmaxx - bminx), bminy + 0.10 * (bmaxy - bminy))
        label_va = "bottom"
        # A box whose foot sits in the panel's bottom band shares it with the
        # scale bar; that one is labelled from its top-left corner instead.
        y_lo, y_hi = ax.get_ylim()
        if anchor[1] < y_lo + 0.16 * (y_hi - y_lo):
            anchor = (anchor[0], bmaxy - 0.10 * (bmaxy - bminy))
            label_va = "top"

        if gis in PRESCRIBED_DOMAINS and gis in site_signed.index and _round_to_cell is not None:
            measured_m = float(site_signed[gis])
            forced = f"{measured_m:+.1f} m measured, {_round_to_cell(measured_m):+.0f} m prescribed"
            forced_colour = INK
        else:
            forced = LABEL_NOT_PRESCRIBED
            forced_colour = INK_MUTED

        ax.annotate(
            text=f"{gis}\n{forced}",
            xy=anchor,
            ha="left",
            va=label_va,
            fontsize=7,
            linespacing=1.35,
            color=forced_colour,
            zorder=8,
            bbox=dict(boxstyle="square,pad=0.2", facecolor="white", edgecolor="none", alpha=0.8),
        )

    name = site_label(run[0], run[-1])

    _panel_title(ax, column, (f"{name}, " if name else "") + f"domains {run[0]}–{run[-1]}")

    ax.set_xticks([])
    ax.set_yticks([])

    # The whole signed metric rests on which side the ocean is, so the map has
    # to say. North-up, so the ocean side of a northward road is the right.
    ax.text(
        0.985 if OCEAN_ON_RIGHT else 0.015,
        0.5,
        "Atlantic Ocean",
        transform=ax.transAxes,
        rotation=90 if OCEAN_ON_RIGHT else -90,
        ha="right" if OCEAN_ON_RIGHT else "left",
        va="center",
        fontsize=8.5,
        color=INK_MUTED,
        style="italic",
    )

    add_scale_bar(ax, length_m=500)
    _north_arrow(ax, x=0.93, y=0.06)

sites_handles = [
    Line2D([], [], color=C_LATE, linewidth=2.0, label=f"{YEAR_TO} road centreline"),
    Line2D([], [], color=C_EARLY, linewidth=0.9, linestyle=(0, (4, 3)), label=f"{YEAR_FROM} road centreline"),
    Line2D([], [], marker="o", linestyle="none", markersize=5,
           markerfacecolor=plt.get_cmap(COLOUR_MAP)(0.6), markeredgecolor="none",
           label=f"{YEAR_FROM} centreline samples, coloured by measured displacement"),
    Line2D([], [], marker="o", linestyle="none", markersize=3, markerfacecolor="0.6",
           markeredgecolor="none", label="samples where the two centrelines coincide"),
]
# two rows: a single row of these four is wider than the printed page
fig_sites.legend(handles=sites_handles, loc="outside lower center", ncol=2,
                 frameon=False)

fig_sites.colorbar(
    plt.cm.ScalarMappable(norm=colour_norm, cmap=COLOUR_MAP),
    ax=list(site_axes),
    label=f"Measured landward displacement of the NC-12 centreline, {YEAR_FROM}–{YEAR_TO} (m)",
    shrink=0.8,
    pad=0.02,
)

save(fig_sites, OUTPUT_FIGURE_SITES, bbox_inches="tight")

print(f"Saved site figure:\n{OUTPUT_FIGURE_SITES} (+ .pdf)")


# =============================================================================
# FIGURE 3 -- DOMAIN MAP: which domains carry a relocation
# =============================================================================
# Every domain in the file, in its real place, coloured by what it carries.
#
# Drawn north-up this is 8 km across and 45 km tall and nothing is legible, so
# the island is rotated 90 degrees clockwise to run south (left) -> north
# (right). That puts the ocean at the bottom and matches the domain axis of
# figure 1, so the two figures read the same way round. Rotation preserves
# distance, so the scale bar is still honest.

ROTATION_DEGREES = -90

rotation_origin = tuple(
    unary_union(domains.geometry.tolist()).centroid.coords[0]
)


def to_strip(frame):
    """Rotate a GeoDataFrame so the island runs left-right."""

    return frame.set_geometry(
        frame.geometry.apply(
            lambda geometry: shapely_rotate(
                geometry,
                ROTATION_DEGREES,
                origin=rotation_origin,
            )
        ),
        crs=frame.crs,
    )


# Carry the classification onto the domain polygons. Domains the road never
# reaches get their own category rather than being lumped in with "no edit":
# there is nothing there to have edited.
domain_map = domains.merge(
    results_df[["domain", "classification", "mean_relocation_m"]],
    left_on=DOMAIN_ID_FIELD,
    right_on="domain",
    how="left",
)

domain_map["classification"] = domain_map["classification"].fillna("no_road")

strip_domains = to_strip(domain_map)
strip_road_from = to_strip(road_from)
strip_road_to = to_strip(road_to)

fig_map, map_ax = plt.subplots(
    figsize=figsize("double", aspect=0.36),
    constrained_layout=True,
)

BACKDROP = {
    "no_road": ("white", LABEL_NO_ROAD),
    "no_edit": (SHADE_NO_EDIT, LABEL_NO_EDIT),
    "redigitized": (SHADE_REDIG, LABEL_REDIG),
}

for key, (colour, _) in BACKDROP.items():

    subset = strip_domains[strip_domains["classification"] == key]

    if subset.empty:
        continue

    subset.plot(ax=map_ax, facecolor=colour, edgecolor="0.6", linewidth=0.35, zorder=1)

# The highlight: relocated domains filled by how far the road moved, so the
# figure says which domains AND how much in one read.
strip_relocated = strip_domains[
    strip_domains["classification"] == "relocated"
]

strip_relocated.plot(ax=map_ax, column="mean_relocation_m", cmap=COLOUR_MAP, norm=colour_norm,
                     edgecolor=INK, linewidth=0.9, zorder=3)

strip_road_from.plot(ax=map_ax, color=C_EARLY, linewidth=0.7, linestyle=(0, (4, 3)), zorder=4)
strip_road_to.plot(ax=map_ax, color=C_LATE, linewidth=0.8, zorder=5)

minx, miny, maxx, maxy = strip_domains.total_bounds

x_pad = 0.01 * (maxx - minx)
y_pad = 0.55 * (maxy - miny)

map_ax.set_xlim(minx - x_pad, maxx + x_pad)
map_ax.set_ylim(miny - y_pad, maxy + y_pad)
map_ax.set_aspect("equal")
map_ax.set_xticks([])
map_ax.set_yticks([])

for spine in map_ax.spines.values():
    spine.set_visible(False)

# Label every tenth domain along the strip, plus the ENDS of each relocation
# run: at the printed width, numbering every relocated domain runs the labels
# into each other, and the bracket above already spans the run.
label_domains = set(range(10, 91, 10)) | {
    run[edge] for run in relocation_sites for edge in (0, -1)
}

for _, domain_row in strip_domains.iterrows():

    if domain_row[DOMAIN_ID_FIELD] not in label_domains:
        continue

    label_point = domain_row.geometry.representative_point()

    is_site = domain_row["classification"] == "relocated"

    map_ax.annotate(
        text=str(domain_row[DOMAIN_ID_FIELD]),
        xy=(label_point.x, maxy + 0.04 * (maxy - miny)),
        ha="center",
        va="bottom",
        fontsize=7,
        color=INK if is_site else INK_MUTED,
        fontweight="bold" if is_site else "normal",
        rotation=90,
        zorder=6,
    )

# Name each relocation site above the strip it belongs to. The statistics
# are in the caption.
for run in relocation_sites:

    run_geometry = strip_domains[
        strip_domains[DOMAIN_ID_FIELD].isin(run)
    ].total_bounds

    name = site_label(run[0], run[-1])

    map_ax.annotate(
        "",
        xy=(run_geometry[0], maxy + 0.30 * (maxy - miny)),
        xytext=(run_geometry[2], maxy + 0.30 * (maxy - miny)),
        arrowprops={"arrowstyle": "|-|", "linewidth": 1.0, "color": INK, "mutation_scale": 3},
        annotation_clip=False,
    )

    map_ax.text(
        (run_geometry[0] + run_geometry[2]) / 2,
        maxy + 0.34 * (maxy - miny),
        (f"{name}\n" if name else "") + f"domains {run[0]}–{run[-1]}",
        ha="center",
        va="bottom",
        fontsize=8.5,
        color=INK,
    )

# Which way is which, now that the map is rotated off north: the ends named,
# the ocean named, and a north arrow pointing along the strip.
map_ax.text(minx, miny - 0.16 * (maxy - miny), "Cape Hatteras (south)",
            ha="left", va="top", fontsize=8.5, color=INK_MUTED)
map_ax.text(maxx, miny - 0.16 * (maxy - miny), "Pea Island (north)",
            ha="right", va="top", fontsize=8.5, color=INK_MUTED)
map_ax.text((minx + maxx) / 2, miny - 0.10 * (maxy - miny), "Atlantic Ocean",
            ha="center", va="top", fontsize=8.5, color=INK_MUTED, style="italic")

map_handles = [
    Patch(facecolor=colour, edgecolor="0.6", linewidth=0.5, label=text)
    for colour, text in BACKDROP.values()
]
map_handles.append(
    Patch(facecolor=plt.get_cmap(COLOUR_MAP)(0.6), edgecolor=INK, linewidth=1.0,
          label=f"relocated ({len(relocated)} domains), filled by mean displacement")
)
map_handles += [
    Line2D([], [], color=C_EARLY, linewidth=0.9, linestyle=(0, (4, 3)), label=f"{YEAR_FROM} road centreline"),
    Line2D([], [], color=C_LATE, linewidth=1.2, label=f"{YEAR_TO} road centreline"),
]

fig_map.legend(handles=map_handles, loc="outside lower center", ncol=2, frameon=False)

add_scale_bar(map_ax, length_m=5000, fraction_x=0.70, fraction_y=0.08)
_east_arrow(map_ax, x=0.94, y=0.11)

fig_map.colorbar(
    plt.cm.ScalarMappable(norm=colour_norm, cmap=COLOUR_MAP),
    ax=map_ax,
    label="Mean landward displacement\nof the NC-12 centreline (m)",
    shrink=0.7,
    pad=0.01,
)

save(fig_map, OUTPUT_FIGURE_DOMAIN_MAP, bbox_inches="tight")

print(f"Saved domain map:\n{OUTPUT_FIGURE_DOMAIN_MAP} (+ .pdf)")


# =============================================================================
# CAPTIONS: the words that used to be on the canvas
# =============================================================================
# One caption per figure, numbers filled from the table this run wrote, so a
# document takes its caption from here and the picture stays a picture.

def write_captions(path):
    site_lines = []
    for column, run in enumerate(relocation_sites):
        stats = relocated[relocated["domain"].isin(run)]
        name = site_label(run[0], run[-1])
        site_lines.append(
            f"({chr(97 + column)}) {name + ', ' if name else ''}domains {run[0]}–{run[-1]}: "
            f"{len(stats)} domains, mean displacement {stats['mean_relocation_m'].mean():.0f} m, "
            f"largest {stats['maximum_relocation_m'].max():.0f} m"
        )
    prescribed = sorted(g for g in PRESCRIBED_DOMAINS if g in site_signed.index)
    in_sites = sorted(g for run in relocation_sites for g in run)
    unprescribed = [g for g in in_sites if g not in PRESCRIBED_DOMAINS]
    unpres_txt = ""
    if unprescribed:
        agree = results_df.set_index("domain")["sign_agreement"]
        parts = [f"{g} (sign agreement {agree.get(g, float('nan')):.2f})" for g in unprescribed]
        unpres_txt = (f" Domains {', '.join(parts)} lie within the digitised change but carry no "
                      f"prescribed move: the two centrelines cross there, so the samples do not agree "
                      f"on a direction and the mean is not a displacement.")
    forced_txt = ", ".join(
        f"GIS {g} {float(site_signed[g]):.1f} → {_round_to_cell(float(site_signed[g])):.0f} m"
        for g in prescribed) if (_round_to_cell is not None and prescribed) else ""
    max_row = relocated.loc[relocated["maximum_relocation_m"].idxmax()]
    lines = [
        f"# Captions — road_relocation_{YEAR_FROM}_{YEAR_TO} figures",
        "",
        f"Written by `HAT_road_relocation_distance.py` on {pd.Timestamp.now():%Y-%m-%d %H:%M}; every number "
        f"is from `{OUTPUT_CSV.name}` beside the figures. The two NC-12 centrelines were digitised from "
        f"1978 and 2008 imagery and stand in for the {YEAR_FROM} and {YEAR_TO} hindcast starts; the "
        f"interval measured is therefore about 30 years. Displacement is the shortest distance from a "
        f"point on the {YEAR_FROM} centreline to the {YEAR_TO} centreline, projected onto the landward "
        f"normal (landward positive), sampled every {SAMPLE_SPACING_M:g} m.",
        "",
        f"**Figure 1 — `{OUTPUT_FIGURE_ALONGSHORE.name}`.** Displacement of the NC-12 centreline, "
        f"{YEAR_FROM}–{YEAR_TO}, per GIS domain from south (left) to north (right). Dark bars are the "
        f"mean signed displacement over the samples in the domain; pale bars the largest displacement "
        f"within it. Grey bands mark domains where the centreline is unchanged between the two surveys "
        f"(light) or was re-digitised without moving, under {REDIGITIZE_MAX_M:.0f} m (darker). "
        f"{len(relocated)} of {len(results_df)} road-carrying domains relocated, all landward; the largest "
        f"displacement is {max_row['maximum_relocation_m']:.0f} m in GIS {int(max_row['domain'])}. "
        f"Brackets name the two relocation sites; the named bands along the foot of the "
        f"panel are the villages.",
        "",
        f"**Figure 2 — `{OUTPUT_FIGURE_SITES.name}`.** The two relocations at true scale, north up, both "
        f"panels at one scale. " + "; ".join(site_lines) + f". The {YEAR_FROM} centreline is drawn dashed "
        f"under its sample points, coloured by the measured landward displacement to the {YEAR_TO} "
        f"centreline (solid); grey points are samples where the two centrelines coincide. Each domain is "
        f"labelled with its mean measured displacement and the value the model is forced with, that mean "
        f"rounded to the nearest {CELL_M:.0f} m Barrier3D cell"
        + (f" ({forced_txt})" if forced_txt else "") + f".{unpres_txt}",
        "",
        f"**Figure 3 — `{OUTPUT_FIGURE_DOMAIN_MAP.name}`.** All {len(domains)} domains in place, the island "
        f"rotated 90° clockwise so that south is to the left, north to the right and the Atlantic at the "
        f"bottom; rotation preserves distance, so the scale bar holds. Relocated domains are outlined and "
        f"filled by their mean displacement on the colour scale shared with Figure 2; the grey classes are "
        f"those of Figure 1, and white domains carry no {YEAR_FROM} road. Every tenth domain and every "
        f"relocated domain is numbered along the top.",
        "",
    ]
    path.write_text("\n".join(lines), encoding="utf-8")
    return path


captions_path = write_captions(OUTPUT_DIR / "CAPTIONS.md")
print(f"Captions:\n{captions_path}")

plt.show()
