"""
Shoreline Data Inventory — Hatteras Island study area
======================================================
Cross-source inventory of shoreline observations available across the
full Hatteras Island study area (all CASCADE domains, Cape Point to
north of Rodanthe). Reports what data you have from each source, how
they overlap in time, and where the gaps are.

The study area is defined by a spatial filter file (see CONFIG) so this
tool can be reused for narrower analyses later (e.g., point it at the
Buxton transects shapefile to focus on the groin study area).

Sources handled
---------------
1. User-digitized wet-dry lines           (GeoJSON/shapefile with date attr)
2. NC Coastal Management historical lines (GeoJSON/shapefile with DATE_ attr)
3. CoastSat satellite-derived shorelines  (folder of per-transect CSVs)

Outputs
-------
  shoreline_inventory_by_year.csv    — one row per year, count per source
  shoreline_inventory_observations.csv — one row per (source, date) observation
  shoreline_inventory_timeline.png   — three-panel visual timeline

Usage
-----
Edit the CONFIG section below, then run:
    python HAT_shoreline_inventory.py
"""

# ============================================================
# CONFIG
# ============================================================

# --- Study area spatial filter ---
# Path to a shapefile or GeoJSON defining the study area. Can be either:
#   • Polygon features (e.g., bounding box, CASCADE domains): merged and
#     optionally buffered.
#   • Line features (e.g., transects): merged and buffered by
#     STUDY_AREA_BUFFER_M to create a filter polygon.
# The script auto-detects which case you're using.
STUDY_AREA_FILTER_PATH = r"/scripts/input_prep/shoreline_inventory/cascade_area.geojson"

# Buffer distance (m) applied around the filter geometry.
# For a bounding box: 500 m to catch drifted historic shorelines.
# For CASCADE domain polygons: 500 m.
# For transect lines: 1000+ m for a corridor filter.
STUDY_AREA_BUFFER_M = 500

# --- Source 1: user-digitized wet-dry lines ---
WET_DRY_PATH     = r"/scripts/input_prep/shoreline_inventory/wet_dry_groin.geojson"
WET_DRY_DATE_COL = "date"

# --- Source 2: NC Coastal Management historical shorelines ---
NC_STATE_PATH     = r"/scripts/input_prep/shoreline_inventory/nc_shorelines.geojson"
NC_STATE_DATE_COL = "DATE_"

# --- Source 3: CoastSat time-series CSVs ---
# GeoJSON with CoastSat transect geometry (used to spatially filter)
COASTSAT_TRANSECT_GEOM = r"C:\Users\hanna\PycharmProjects\CASCADE\scripts\input_preperation\CoastSat\CoastSat_transect_layer.geojson"
# ID column in the CoastSat geometry file
COASTSAT_TRANSECT_ID_COL = "id"
# Root folder containing per-transect CSVs in site subfolders
COASTSAT_ROOT_DIR = r"/scripts/input_prep/CoastSat/coastsat_timeseries"
# A CoastSat date is "well-covered" if this many study-area transects
# reported an observation on that date
COASTSAT_WELL_COVERED_MIN_TRANSECTS = 10

# --- Visualization ---
# Whether to produce spatial and alongshore visualizations
MAKE_VISUALIZATIONS = True

# Colormap for year-based shoreline coloring.
# 'viridis' is perceptually uniform, colorblind-safe, and reads intuitively
# as "dark = original, bright = new". Other good options: 'plasma', 'cividis'.
SHORELINE_COLORMAP = "viridis"

# Number of representative CoastSat snapshot dates to reconstruct.
# The script picks these to span the full temporal record (evenly-spaced bins)
# and picks the most-complete date within each bin (highest transect coverage).
N_COASTSAT_SNAPSHOTS = 10

# Basemap for the spatial map view. Requires the 'contextily' library:
#   pip install contextily
# Options:
#   "none"           — plain background (default, no dependencies)
#   "esri_satellite" — ESRI World Imagery (aerial photography)
#   "esri_topo"      — ESRI Topographic (labeled road_offset, terrain)
#   "carto_light"    — Carto positron (clean, high-contrast for shorelines)
# If contextily is unavailable, the script falls back to "none" with a warning.
BASEMAP = "carto_light"

# --- Output ---
OUTPUT_DIR = r"/scripts/input_prep/shoreline_inventory/comparison"

# CRS for spatial operations (UTM 18N covers NC Outer Banks)
PROJECTED_CRS = "EPSG:32618"

# Geographic labels drawn on the right edge of each map panel.
# Values are (name, northing_meters_UTM18N). Approximate positions of major
# Hatteras Island communities. TO FIX A LABEL POSITION:
#   1. Open your data in ArcGIS Pro or QGIS
#   2. Click on the actual location of the town on the map
#   3. Read the y-coordinate (northing) from the status bar (should be a
#      value between ~3,890,000 and ~3,960,000 for Hatteras)
#   4. Update the tuple below with that value
PLACE_LABELS = [
    ("Cape Point", 3_897_500),   # southern hook, Cape Hatteras
    ("Buxton",     3_902_000),   # village north of the cape
    ("Avon",       3_911_000),   # mid-island community
    ("Salvo",      3_933_000),   # tri-village south
    ("Waves",      3_936_000),   # tri-village middle
    ("Rodanthe",   3_939_000),   # tri-village north
    ("Pea Island", 3_950_000),   # NWR north of Rodanthe
]


# ============================================================
# IMPORTS
# ============================================================
import os
import glob
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import geopandas as gpd
from shapely.geometry import Point
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import matplotlib.patches as mpatches

warnings.filterwarnings("ignore")

# Consistent source labels and colors used everywhere
SOURCE_LABELS = {
    "wet_dry":  "Aerial wet-dry lines",
    "nc_state": "NC Coastal Mgmt",
    "coastsat": "CoastSat (satellite)",
}
SOURCE_COLORS = {
    "wet_dry":  "#9C4D51",   # muted brick red
    "nc_state": "#3E5F7E",   # steel blue
    "coastsat": "#C27B49",   # warm terracotta
}
SOURCE_ORDER = ["wet_dry", "nc_state", "coastsat"]


# ============================================================
# HELPERS
# ============================================================

def fix_crs(gdf: gpd.GeoDataFrame, fallback: str = PROJECTED_CRS) -> gpd.GeoDataFrame:
    """Repair GeoDataFrames whose CRS is missing or unrecognized."""
    try:
        if gdf.crs is None:
            raise ValueError("No CRS defined")
        gdf.crs.to_authority()
        return gdf
    except Exception:
        print(f"    WARNING: CRS '{gdf.crs}' unrecognized, assuming {fallback}")
        return gdf.set_crs(fallback, allow_override=True)


def load_study_area_filter(path: str, buffer_m: float):
    """Load a shapefile/GeoJSON defining the study area and return a
    single buffered polygon used as the spatial filter for every source.

    Auto-detects whether the input contains polygon features (CASCADE
    domains) or line features (transects), and handles each accordingly:
      - Polygons: unioned, then buffered by `buffer_m`.
      - Lines:    buffered by `buffer_m`, then unioned.
      - Points:   buffered by `buffer_m`, then unioned.
    """
    print(f"Loading study area filter: {os.path.basename(path)}")
    gdf = gpd.read_file(path)
    gdf = fix_crs(gdf)
    print(f"  {len(gdf)} features, CRS: {gdf.crs}")

    geom_types = gdf.geom_type.value_counts().to_dict()
    print(f"  Geometry types: {geom_types}")

    # Reproject to projected CRS so buffer is in meters
    gdf_proj = gdf.to_crs(PROJECTED_CRS)
    merged = gdf_proj.geometry.unary_union

    # Buffer if we have any non-polygon features, or if buffer_m > 0 for polygons
    has_polygons = any(t in geom_types for t in ("Polygon", "MultiPolygon"))
    has_lines    = any(t in geom_types for t in ("LineString", "MultiLineString"))
    has_points   = any(t in geom_types for t in ("Point", "MultiPoint"))

    if has_lines or has_points:
        # Need to buffer to create an area
        filter_geom = merged.buffer(buffer_m)
        print(f"  Detected line/point features → buffered by {buffer_m} m")
    elif has_polygons and buffer_m > 0:
        # Polygons with optional buffer
        filter_geom = merged.buffer(buffer_m)
        print(f"  Detected polygon features → buffered by {buffer_m} m")
    else:
        # Polygons, no buffer
        filter_geom = merged
        print(f"  Detected polygon features → used as-is (no buffer)")

    filter_gdf = gpd.GeoDataFrame(geometry=[filter_geom], crs=PROJECTED_CRS)
    print(f"  Filter polygon area: {filter_geom.area / 1e6:.2f} km²")
    # bbox in km for context
    minx, miny, maxx, maxy = filter_geom.bounds
    print(f"  Extent (projected): {(maxx - minx) / 1000:.1f} km × "
          f"{(maxy - miny) / 1000:.1f} km\n")
    return filter_gdf


# ============================================================
# LOADER: shapefile / GeoJSON with date attribute
# ============================================================

def parse_date_column(series: pd.Series, col_name: str = "") -> pd.Series:
    """Parse a date column that could be in one of several formats:

    - String dates: '1972-07-30', '7/22/1998', 'July 30, 1972', etc.
    - Numeric milliseconds since Unix epoch: 81302400000 → 1972-07-30
      (This format is produced by some ArcGIS / QGIS exports.)
    - Numeric year-only integers: 1972 → 1972-01-01

    Returns tz-naive datetime series with NaT for unparseable values.
    """
    # Case 1: string dtype — use pandas' flexible parser
    if series.dtype == object or str(series.dtype).startswith("string"):
        return pd.to_datetime(series, errors="coerce")

    # Case 2: numeric dtype — decide between ms-since-epoch and year-integers
    if series.dtype.kind in "iuf":
        s = series.dropna()
        if len(s) == 0:
            return pd.to_datetime(series, errors="coerce")

        abs_max = s.abs().max()

        # Year integers are small (typically 1800-2200 range)
        if abs_max < 3000:
            print(f"  Date column '{col_name}' looks like year integers "
                  f"(range: {int(s.min())}-{int(s.max())})")
            return pd.to_datetime(series.astype("Int64").astype(str),
                                   format="%Y", errors="coerce")

        # Milliseconds since epoch: |values| roughly 10^10 to 10^13
        if abs_max > 1e9:
            print(f"  Date column '{col_name}' looks like Unix milliseconds "
                  f"since epoch → parsing with unit='ms'")
            return pd.to_datetime(series, unit="ms", errors="coerce")

        # Fallback: treat as generic numeric
        return pd.to_datetime(series, errors="coerce")

    # Case 3: already datetime
    if pd.api.types.is_datetime64_any_dtype(series):
        return series

    # Last resort
    return pd.to_datetime(series, errors="coerce")


def load_shapefile_dates(path: str,
                          date_col: str,
                          filter_polygon: gpd.GeoDataFrame,
                          source_key: str) -> pd.DataFrame:
    """Load a shoreline shapefile/GeoJSON, spatially filter to the study area,
    and extract observation dates.

    Returns a DataFrame with columns:
      source, date, year, extent_m
    where `extent_m` is the length of the shoreline within the study area
    (in meters), giving a proxy for alongshore coverage per feature.
    """
    print(f"[{SOURCE_LABELS[source_key]}] loading: {os.path.basename(path)}")

    if not os.path.exists(path):
        print(f"  FILE NOT FOUND — skipping this source\n")
        return pd.DataFrame(columns=["source", "date", "year", "extent_m"])

    gdf = gpd.read_file(path)
    gdf = fix_crs(gdf)
    print(f"  {len(gdf)} raw features, CRS: {gdf.crs}")

    if date_col not in gdf.columns:
        print(f"  ERROR: date column '{date_col}' not in file.")
        print(f"  Available columns: {list(gdf.columns)}\n")
        return pd.DataFrame(columns=["source", "date", "year", "extent_m"])

    # Reproject to projected CRS for spatial operations
    gdf_proj = gdf.to_crs(PROJECTED_CRS)

    # Spatial filter: keep features that intersect the study area
    filter_geom = filter_polygon.geometry.iloc[0]
    in_area = gdf_proj.geometry.intersects(filter_geom)
    gdf_proj = gdf_proj[in_area].copy()
    print(f"  After spatial filter: {len(gdf_proj)} features in study area")

    if len(gdf_proj) == 0:
        return pd.DataFrame(columns=["source", "date", "year", "extent_m"])

    # Parse dates using the format-aware parser
    parsed = parse_date_column(gdf_proj[date_col], col_name=date_col)
    n_bad = parsed.isna().sum()
    if n_bad > 0:
        print(f"  WARNING: {n_bad} features had unparseable dates — excluded")
    gdf_proj = gdf_proj.assign(_parsed_date=parsed).dropna(subset=["_parsed_date"])

    # Compute the length of each shoreline feature within the study area
    # (proxy for alongshore extent of that observation)
    def clip_length(geom):
        try:
            clipped = geom.intersection(filter_geom)
            if clipped.is_empty:
                return 0.0
            return float(clipped.length)
        except Exception:
            return 0.0

    gdf_proj["extent_m"] = gdf_proj.geometry.apply(clip_length)

    # Aggregate by unique date — many shoreline shapefiles store single
    # shorelines as multiple line segments (topologically clipped, split at
    # grid edges, etc.). Counting each segment as a separate observation
    # would inflate high-density years like 1980 (~140 segments for one
    # shoreline). We collapse to one row per date and sum the segment lengths.
    n_features_raw = len(gdf_proj)
    per_date = gdf_proj.groupby("_parsed_date").agg(
        extent_m=("extent_m", "sum"),
        n_segments=("extent_m", "count"),
    ).reset_index().rename(columns={"_parsed_date": "date"})
    per_date["source"] = source_key
    per_date["year"] = per_date["date"].dt.year
    per_date = per_date[["source", "date", "year", "extent_m", "n_segments"]]

    print(f"  → {n_features_raw} raw segments collapsed to "
          f"{len(per_date)} unique-date observations, "
          f"{per_date['year'].nunique()} unique years")
    if len(per_date) > 0:
        print(f"  → total shoreline length: "
              f"{per_date['extent_m'].sum() / 1000:.1f} km "
              f"(mean per observation: "
              f"{per_date['extent_m'].mean() / 1000:.2f} km, "
              f"mean segments per observation: "
              f"{per_date['n_segments'].mean():.1f})\n")
    return per_date


# ============================================================
# LOADER: CoastSat CSV time-series
# ============================================================

def collect_csv_map(root_dir: str) -> dict:
    """Walk one level of subfolders under root_dir and return
    {csv_stem: full_filepath} for every CSV found."""
    csv_map = {}
    if not os.path.isdir(root_dir):
        print(f"  ROOT NOT FOUND: {root_dir}")
        return csv_map
    subfolders = [
        os.path.join(root_dir, d)
        for d in sorted(os.listdir(root_dir))
        if os.path.isdir(os.path.join(root_dir, d))
    ]
    for sf in subfolders:
        for fpath in glob.glob(os.path.join(sf, "*.csv")):
            stem = os.path.splitext(os.path.basename(fpath))[0]
            csv_map[stem] = fpath
    return csv_map


def load_coastsat_dates(transect_geom_path: str,
                        transect_id_col: str,
                        root_dir: str,
                        filter_polygon: gpd.GeoDataFrame,
                        min_transects_well_covered: int,
                        source_key: str = "coastsat"):
    """Load CoastSat observations, filter to study area, and aggregate.

    Returns:
        per_date_summary : DataFrame with columns source/date/year/n_transects_covered
        raw_observations : DataFrame with columns transect_id/date/chainage_m
                           (used later for shoreline reconstruction; empty if none)
        transect_geoms   : GeoDataFrame of study-area CoastSat transects
                           with columns transect_id/geometry (LineString in PROJECTED_CRS)
                           (used later for reconstructing shoreline points)
    """
    print(f"[{SOURCE_LABELS[source_key]}] processing CoastSat data")

    empty_summary = pd.DataFrame(columns=["source", "date", "year",
                                          "n_transects_covered"])
    empty_raw     = pd.DataFrame(columns=["transect_id", "date", "chainage_m"])
    empty_geoms   = gpd.GeoDataFrame(columns=["transect_id", "geometry"],
                                     geometry="geometry", crs=PROJECTED_CRS)

    if not os.path.exists(transect_geom_path):
        print(f"  Transect geometry not found: {transect_geom_path}")
        print(f"  Skipping CoastSat source\n")
        return empty_summary, empty_raw, empty_geoms

    # Load CoastSat transect geometry (LineStrings)
    gdf = gpd.read_file(transect_geom_path)
    gdf = fix_crs(gdf, fallback="EPSG:4326")
    print(f"  Loaded {len(gdf):,} CoastSat transects globally")

    if transect_id_col not in gdf.columns:
        print(f"  ERROR: ID column '{transect_id_col}' not found. "
              f"Available: {list(gdf.columns)}\n")
        return empty_summary, empty_raw, empty_geoms

    # Reproject transects to projected CRS
    gdf_proj = gdf.to_crs(PROJECTED_CRS)

    # Extract origin point of each transect for spatial filtering
    def to_origin(g):
        if g.geom_type == "LineString":
            return Point(g.coords[0])
        if g.geom_type == "MultiLineString":
            return Point(g.geoms[0].coords[0])
        return g

    gdf_proj = gdf_proj.copy()
    gdf_proj["origin"] = gdf_proj.geometry.apply(to_origin)

    # Spatial filter: keep transects whose origin is in the study area
    filter_geom = filter_polygon.geometry.iloc[0]
    in_area = gpd.GeoSeries(gdf_proj["origin"], crs=PROJECTED_CRS).within(filter_geom)
    gdf_proj = gdf_proj[in_area].copy()
    print(f"  In study area (within {STUDY_AREA_BUFFER_M} m of filter): "
          f"{len(gdf_proj):,} transects")

    if len(gdf_proj) == 0:
        return empty_summary, empty_raw, empty_geoms

    # Build the transect_geoms comparison (kept for reconstruction)
    transect_geoms = gdf_proj[[transect_id_col, "geometry"]].rename(
        columns={transect_id_col: "transect_id"}
    ).reset_index(drop=True)

    # Look up the CSV for each study-area transect
    csv_map = collect_csv_map(root_dir)
    print(f"  Found {len(csv_map):,} CSVs in the data folder")

    study_ids = set(gdf_proj[transect_id_col].astype(str).values)
    matched = [tid for tid in study_ids if tid in csv_map]
    missing = study_ids - set(matched)
    print(f"  With CSV data available: {len(matched):,} / {len(study_ids):,}")
    if missing:
        print(f"  Missing CSVs: {len(missing):,} transects "
              f"(satellite covered but not downloaded)")

    # Read date + chainage columns from each CSV
    all_obs = []
    for i, tid in enumerate(matched):
        if (i + 1) % 25 == 0 or (i + 1) == len(matched):
            print(f"    Reading CSVs: {i + 1}/{len(matched)}", end="\r")
        try:
            df = pd.read_csv(csv_map[tid], usecols=[0, 1])
            df.columns = ["date", "chainage_m"]
            df["date"] = pd.to_datetime(df["date"], utc=True, errors="coerce")
            df["chainage_m"] = pd.to_numeric(df["chainage_m"], errors="coerce")
            df = df.dropna(subset=["date", "chainage_m"])
            df["transect_id"] = tid
            all_obs.append(df)
        except Exception as e:
            print(f"\n    ERROR reading {tid}: {e}")

    print()
    if not all_obs:
        return empty_summary, empty_raw, transect_geoms

    full = pd.concat(all_obs, ignore_index=True)
    # Floor to day using UTC to collapse same-day satellite passes, then strip
    # the timezone so these dates can be concat'd/sorted alongside tz-naive
    # dates from the shapefile sources.
    full["date_key"] = (full["date"]
                        .dt.tz_convert("UTC")
                        .dt.floor("D")
                        .dt.tz_localize(None))

    # Per-date summary (for inventory)
    per_date = (full.groupby("date_key")["transect_id"]
                    .nunique()
                    .reset_index()
                    .rename(columns={"transect_id": "n_transects_covered"}))
    per_date["source"] = source_key
    per_date["date"] = pd.to_datetime(per_date["date_key"])
    per_date["year"] = per_date["date"].dt.year
    per_date = per_date[["source", "date", "year", "n_transects_covered"]]

    # Raw observations (kept for visualization)
    raw = full[["transect_id", "date_key", "chainage_m"]].copy()
    raw = raw.rename(columns={"date_key": "date"})
    raw["date"] = pd.to_datetime(raw["date"])

    print(f"  → {len(per_date):,} unique observation dates spanning "
          f"{per_date['year'].min()}–{per_date['year'].max()}")
    n_well = (per_date["n_transects_covered"] >= min_transects_well_covered).sum()
    print(f"  → {n_well:,} well-covered dates "
          f"(≥{min_transects_well_covered} transects)")
    print(f"  → {len(raw):,} raw per-transect observations retained\n")
    return per_date, raw, transect_geoms


# ============================================================
# AGGREGATION
# ============================================================

def build_year_summary(observations: pd.DataFrame,
                        min_transects_well_covered: int) -> pd.DataFrame:
    """One row per year with per-source counts."""
    if len(observations) == 0:
        return pd.DataFrame()

    # Total observations per source per year
    pivot = (observations.groupby(["year", "source"])
                         .size()
                         .unstack(fill_value=0)
                         .reindex(columns=SOURCE_ORDER, fill_value=0))

    # For CoastSat also compute well-covered count
    cs = observations[observations["source"] == "coastsat"].copy()
    if len(cs) > 0:
        cs_well = (cs[cs["n_transects_covered"] >= min_transects_well_covered]
                      .groupby("year").size())
        pivot["coastsat_well_covered"] = cs_well.reindex(pivot.index,
                                                           fill_value=0)
    else:
        pivot["coastsat_well_covered"] = 0

    pivot["total_all_sources"] = pivot[SOURCE_ORDER].sum(axis=1)
    pivot["distinct_sources"]  = (pivot[SOURCE_ORDER] > 0).sum(axis=1)

    # Reindex to include every year in the span (fills years with 0s)
    full_range = range(int(pivot.index.min()), int(pivot.index.max()) + 1)
    pivot = pivot.reindex(full_range, fill_value=0)
    pivot.index.name = "year"

    # Rename columns to be human-friendly
    pivot = pivot.rename(columns={
        "wet_dry":  "wet_dry",
        "nc_state": "nc_state",
        "coastsat": "coastsat_all",
    })
    # Reorder
    col_order = ["wet_dry", "nc_state", "coastsat_all", "coastsat_well_covered",
                 "total_all_sources", "distinct_sources"]
    return pivot[col_order].reset_index()


# ============================================================
# PLOTTING
# ============================================================

def plot_inventory(observations: pd.DataFrame,
                    year_summary: pd.DataFrame,
                    min_transects_well_covered: int,
                    output_path: str):
    """Three-panel visual inventory of shoreline availability."""

    from matplotlib.gridspec import GridSpec

    fig = plt.figure(figsize=(15, 10))
    gs = GridSpec(3, 1, height_ratios=[1.3, 1.0, 1.0],
                  hspace=0.45, top=0.94, bottom=0.06, left=0.07, right=0.98)

    # ─── Panel 1: Stacked bar of observations per year ────────────
    ax1 = fig.add_subplot(gs[0])
    if len(year_summary) > 0:
        years = year_summary["year"].values
        bottom = np.zeros(len(years))
        for src in SOURCE_ORDER:
            col_name = "coastsat_all" if src == "coastsat" else src
            vals = year_summary[col_name].values
            ax1.bar(years, vals, bottom=bottom,
                    color=SOURCE_COLORS[src], label=SOURCE_LABELS[src],
                    width=0.9, edgecolor="none")
            bottom += vals

        ax1.set_xlabel("Year", fontsize=10)
        ax1.set_ylabel("Number of shoreline observations", fontsize=10)
        ax1.set_title("Shoreline observations per year, by source", fontsize=11)
        ax1.legend(loc="upper left", fontsize=9)
        ax1.grid(True, axis="y", alpha=0.3)
        ax1.set_xlim(years.min() - 1, years.max() + 1)

    # ─── Panel 2: Individual observation timeline ─────────────────
    ax2 = fig.add_subplot(gs[1], sharex=None)
    y_pos = {src: i for i, src in enumerate(SOURCE_ORDER)}
    for src in SOURCE_ORDER:
        obs = observations[observations["source"] == src]
        if len(obs) == 0:
            continue
        # For CoastSat, distinguish well-covered vs. sparse
        if src == "coastsat":
            well = obs[obs["n_transects_covered"] >= min_transects_well_covered]
            sparse = obs[obs["n_transects_covered"] < min_transects_well_covered]
            ax2.scatter(sparse["date"], [y_pos[src]] * len(sparse),
                        marker="|", s=80, color=SOURCE_COLORS[src], alpha=0.25,
                        linewidths=0.8)
            ax2.scatter(well["date"], [y_pos[src]] * len(well),
                        marker="|", s=100, color=SOURCE_COLORS[src], alpha=0.9,
                        linewidths=1.2)
        else:
            ax2.scatter(obs["date"], [y_pos[src]] * len(obs),
                        marker="|", s=140, color=SOURCE_COLORS[src], alpha=0.9,
                        linewidths=1.8)

    ax2.set_yticks(list(y_pos.values()))
    ax2.set_yticklabels([SOURCE_LABELS[s] for s in SOURCE_ORDER], fontsize=9)
    ax2.set_ylim(-0.6, len(SOURCE_ORDER) - 0.4)
    ax2.set_xlabel("Date", fontsize=10)
    ax2.set_title(
        f"Every observation as a tick mark  "
        f"(faded CoastSat = sparse, <{min_transects_well_covered} transects)",
        fontsize=11,
    )
    ax2.grid(True, axis="x", alpha=0.3)

    # ─── Panel 3: Cumulative observations over time ───────────────
    ax3 = fig.add_subplot(gs[2])
    if len(observations) > 0:
        obs_sorted = observations.sort_values("date").copy()
        for src in SOURCE_ORDER:
            src_obs = obs_sorted[obs_sorted["source"] == src].copy()
            if len(src_obs) == 0:
                continue
            src_obs["cumulative"] = np.arange(1, len(src_obs) + 1)
            ax3.step(src_obs["date"], src_obs["cumulative"],
                     where="post", color=SOURCE_COLORS[src], lw=1.8,
                     label=SOURCE_LABELS[src])

        ax3.set_xlabel("Date", fontsize=10)
        ax3.set_ylabel("Cumulative observations", fontsize=10)
        ax3.set_title("Cumulative shorelines available over time",
                      fontsize=11)
        ax3.legend(loc="upper left", fontsize=9)
        ax3.grid(True, alpha=0.3)

    fig.suptitle(
        "Shoreline data inventory — Buxton study area",
        fontsize=13, fontweight="bold", y=0.99,
    )
    plt.savefig(output_path, dpi=150, bbox_inches="tight")
    print(f"→ saved {output_path}")
    plt.close(fig)


# ============================================================
# VISUALIZATION: shoreline reconstruction & spatial plots
# ============================================================

def pick_representative_coastsat_dates(per_date_summary: pd.DataFrame,
                                       n_snapshots: int,
                                       min_transects: int) -> pd.DataFrame:
    """Pick N CoastSat dates that (a) span the full temporal record and
    (b) each have high transect coverage.

    Method: divide the time span into N evenly-spaced bins. In each bin,
    pick the date with the highest n_transects_covered.
    """
    if len(per_date_summary) == 0:
        return per_date_summary.iloc[:0].copy()

    candidates = per_date_summary[
        per_date_summary["n_transects_covered"] >= min_transects
    ].copy()
    if len(candidates) == 0:
        # No dates meet the strict threshold — relax and use all
        candidates = per_date_summary.copy()

    # Sort by date, then bin
    candidates = candidates.sort_values("date").reset_index(drop=True)
    t_min = candidates["date"].min()
    t_max = candidates["date"].max()

    # Handle the degenerate case
    if t_max == t_min:
        return candidates.head(n_snapshots).copy()

    bin_edges = pd.date_range(t_min, t_max, periods=n_snapshots + 1)
    picks = []
    for i in range(n_snapshots):
        lo, hi = bin_edges[i], bin_edges[i + 1]
        mask = (candidates["date"] >= lo) & (candidates["date"] <= hi)
        bin_dates = candidates[mask]
        if len(bin_dates) == 0:
            continue
        best = bin_dates.loc[bin_dates["n_transects_covered"].idxmax()]
        picks.append(best)

    result = pd.DataFrame(picks).drop_duplicates(subset=["date"])
    result = result.sort_values("date").reset_index(drop=True)
    return result


def get_transect_direction(geom):
    """Return unit direction vector (dx, dy) from a transect LineString's
    origin (first vertex) to its endpoint (last vertex)."""
    if geom.geom_type == "MultiLineString":
        geom = geom.geoms[0]
    coords = list(geom.coords)
    x0, y0 = coords[0][:2]
    x1, y1 = coords[-1][:2]
    dx, dy = x1 - x0, y1 - y0
    length = np.hypot(dx, dy)
    if length == 0:
        return 0.0, 0.0
    return dx / length, dy / length


def reconstruct_coastsat_shoreline_points(raw_obs: pd.DataFrame,
                                          transect_geoms: gpd.GeoDataFrame,
                                          target_dates: pd.Series) -> pd.DataFrame:
    """For each (date, transect_id) in target_dates, compute the shoreline
    point (x, y) in the projected CRS by walking along the transect from
    its origin for `chainage_m` meters.

    Returns DataFrame with columns:
      date, transect_id, x, y, chainage_m, origin_x, origin_y
    """
    if len(raw_obs) == 0 or len(transect_geoms) == 0:
        return pd.DataFrame(columns=["date", "transect_id", "x", "y",
                                     "chainage_m", "origin_x", "origin_y"])

    # Pre-compute origin and direction per transect
    geom_info = {}
    for _, row in transect_geoms.iterrows():
        g = row["geometry"]
        if g.geom_type == "MultiLineString":
            g = g.geoms[0]
        coords = list(g.coords)
        origin_x, origin_y = coords[0][:2]
        dx, dy = get_transect_direction(g)
        geom_info[row["transect_id"]] = (origin_x, origin_y, dx, dy)

    # Filter observations to target dates
    target_dates_set = set(pd.to_datetime(target_dates).values)
    obs = raw_obs[raw_obs["date"].isin(target_dates_set)].copy()
    if len(obs) == 0:
        return pd.DataFrame(columns=["date", "transect_id", "x", "y",
                                     "chainage_m", "origin_x", "origin_y"])

    xs, ys, ox, oy = [], [], [], []
    for tid, chn in zip(obs["transect_id"].values, obs["chainage_m"].values):
        info = geom_info.get(tid)
        if info is None:
            xs.append(np.nan); ys.append(np.nan)
            ox.append(np.nan); oy.append(np.nan)
            continue
        origin_x, origin_y, dx, dy = info
        xs.append(origin_x + chn * dx)
        ys.append(origin_y + chn * dy)
        ox.append(origin_x)
        oy.append(origin_y)

    obs["x"], obs["y"] = xs, ys
    obs["origin_x"], obs["origin_y"] = ox, oy
    obs = obs.dropna(subset=["x", "y"])
    return obs[["date", "transect_id", "x", "y",
                "chainage_m", "origin_x", "origin_y"]]


def compute_alongshore_position(transect_geoms: gpd.GeoDataFrame) -> pd.DataFrame:
    """Assign each CoastSat transect an alongshore position (m) — cumulative
    distance from the southernmost transect origin, walking north via
    nearest-neighbor.

    Returns DataFrame: transect_id, origin_x, origin_y, alongshore_m
    """
    if len(transect_geoms) == 0:
        return pd.DataFrame(columns=["transect_id", "origin_x", "origin_y",
                                     "alongshore_m"])

    rows = []
    for _, row in transect_geoms.iterrows():
        g = row["geometry"]
        if g.geom_type == "MultiLineString":
            g = g.geoms[0]
        coords = list(g.coords)
        rows.append({
            "transect_id": row["transect_id"],
            "origin_x":    coords[0][0],
            "origin_y":    coords[0][1],
        })
    df = pd.DataFrame(rows)

    # Nearest-neighbor path starting from the southernmost origin
    remaining = df.copy().reset_index(drop=True)
    start_idx = remaining["origin_y"].idxmin()
    ordered = [remaining.loc[start_idx].to_dict()]
    ordered[0]["alongshore_m"] = 0.0
    remaining = remaining.drop(start_idx).reset_index(drop=True)

    cum = 0.0
    cur_x = ordered[-1]["origin_x"]
    cur_y = ordered[-1]["origin_y"]

    while len(remaining) > 0:
        d2 = ((remaining["origin_x"] - cur_x) ** 2 +
              (remaining["origin_y"] - cur_y) ** 2)
        nn_idx = d2.idxmin()
        dist = float(np.sqrt(d2.loc[nn_idx]))
        cum += dist
        next_row = remaining.loc[nn_idx].to_dict()
        next_row["alongshore_m"] = cum
        ordered.append(next_row)
        cur_x = next_row["origin_x"]
        cur_y = next_row["origin_y"]
        remaining = remaining.drop(nn_idx).reset_index(drop=True)

    return pd.DataFrame(ordered)


# ============================================================
# VISUALIZATION: map view
# ============================================================

def _add_basemap_if_requested(ax, basemap_kind: str):
    """Attempt to overlay a basemap; fall back silently if unavailable."""
    if basemap_kind == "none":
        return
    try:
        import contextily as cx
    except ImportError:
        print(f"    contextily not installed; skipping '{basemap_kind}' basemap "
              f"(pip install contextily to enable)")
        return

    providers = {
        "esri_satellite": cx.providers.Esri.WorldImagery,
        "esri_topo":      cx.providers.Esri.WorldTopoMap,
        "carto_light":    cx.providers.CartoDB.Positron,
    }
    source = providers.get(basemap_kind)
    if source is None:
        print(f"    unknown basemap '{basemap_kind}', skipping")
        return
    try:
        cx.add_basemap(ax, crs=PROJECTED_CRS, source=source,
                       attribution="", zorder=0.5)
    except Exception as e:
        print(f"    basemap fetch failed ({e}); continuing without basemap")


def _compute_data_bounds(*gdfs_and_frames, filter_polygon):
    """Return zoom bounds from the union of all shoreline data, with a small
    buffer. Falls back to the filter polygon bounds if no data provided."""
    xmins, ymins, xmaxs, ymaxs = [], [], [], []
    for obj in gdfs_and_frames:
        if obj is None or len(obj) == 0:
            continue
        if isinstance(obj, gpd.GeoDataFrame):
            b = obj.total_bounds
            xmins.append(b[0]); ymins.append(b[1])
            xmaxs.append(b[2]); ymaxs.append(b[3])
        elif isinstance(obj, pd.DataFrame) and {"x", "y"}.issubset(obj.columns):
            xmins.append(obj["x"].min()); ymins.append(obj["y"].min())
            xmaxs.append(obj["x"].max()); ymaxs.append(obj["y"].max())
    if not xmins:
        return filter_polygon.geometry.iloc[0].bounds
    xmin, ymin = min(xmins), min(ymins)
    xmax, ymax = max(xmaxs), max(ymaxs)
    # Expand horizontally more than vertically because coast is narrow strip
    pad_x = max((xmax - xmin) * 0.20, 500)
    pad_y = (ymax - ymin) * 0.02
    return xmin - pad_x, ymin - pad_y, xmax + pad_x, ymax + pad_y


def _draw_scale_bar(ax, xmin, ymin, xmax, ymax, length_km=5):
    """Small horizontal scale bar in the lower-left of an axis."""
    length_m = length_km * 1000
    # Bottom-left inset position
    x0 = xmin + (xmax - xmin) * 0.05
    y0 = ymin + (ymax - ymin) * 0.03
    ax.plot([x0, x0 + length_m], [y0, y0], color="black", lw=2.5,
            solid_capstyle="butt", zorder=10)
    ax.text(x0 + length_m / 2, y0 + (ymax - ymin) * 0.008,
            f"{length_km} km", fontsize=8, ha="center", va="bottom",
            color="black", zorder=10)


def plot_shoreline_map(wet_dry_gdf: gpd.GeoDataFrame,
                        nc_state_gdf: gpd.GeoDataFrame,
                        coastsat_recon: pd.DataFrame,
                        coastsat_alongshore: pd.DataFrame,
                        filter_polygon: gpd.GeoDataFrame,
                        colormap: str,
                        basemap_kind: str,
                        output_path: str):
    """Four-panel horizontal map view of shoreline coverage over time.

    Design decisions:
      - Panels share y-axis (only leftmost shows northing labels)
      - Named locations labeled as secondary y-tick labels on leftmost panel
      - Axes in km (raw UTM /1000) rather than meters with 1e6 raw_offset
      - Rank-based year coloring so historic (1849) and modern (2024) years are
        both distinguishable rather than compressed into narrow color bands
      - Zoom to actual data extent + horizontal buffer rather than the full
        study area polygon (island is much narrower than the polygon)
      - Optional satellite / topo / positron basemap via contextily
    """
    from matplotlib.gridspec import GridSpec
    from matplotlib.colors import BoundaryNorm, ListedColormap
    from matplotlib.cm import ScalarMappable
    from matplotlib.lines import Line2D
    from matplotlib.ticker import FuncFormatter
    from matplotlib.patches import Rectangle

    # ─── Collect unique observation years for rank-based coloring ─────────
    all_years = []
    if len(wet_dry_gdf) > 0:
        all_years.extend(wet_dry_gdf["_year"].tolist())
    if len(nc_state_gdf) > 0:
        all_years.extend(nc_state_gdf["_year"].tolist())
    if len(coastsat_recon) > 0:
        all_years.extend(pd.to_datetime(coastsat_recon["date"]).dt.year.tolist())

    unique_years = sorted(set(int(y) for y in all_years))
    n_years = len(unique_years)
    year_to_rank = {yr: i for i, yr in enumerate(unique_years)}
    cmap_obj = plt.get_cmap(colormap)

    def year_color(yr):
        rank = year_to_rank.get(int(yr))
        if rank is None:
            return (0.5, 0.5, 0.5, 1.0)
        return cmap_obj(rank / max(n_years - 1, 1))

    # ─── Figure and axes ──────────────────────────────────────────────────
    # Independent y-axes on each panel (sharey with equal-aspect + geopandas
    # is a known incompatibility — the axes have to be independent).
    fig = plt.figure(figsize=(18, 12))
    gs = GridSpec(1, 4, wspace=0.02,
                  top=0.93, bottom=0.20, left=0.04, right=0.99)
    axes = [fig.add_subplot(gs[0, i]) for i in range(4)]
    titles = ["Aerial wet-dry lines", "NC Coastal Mgmt",
              "CoastSat (satellite)", "All sources overlaid"]

    if n_years == 0:
        for ax, title in zip(axes, titles):
            ax.text(0.5, 0.5, "No data available",
                    transform=ax.transAxes, ha="center", va="center")
            ax.set_title(title)
        plt.savefig(output_path, dpi=200, bbox_inches="tight")
        plt.close(fig)
        return

    # ─── Use the study area polygon bounds so every panel shows the same
    # extent as the input bounding box. Extra right padding gives labels
    # room outside the island's east edge, and importantly makes the plot
    # aspect wide enough that panels fill their grid cells rather than
    # shrinking to a thin strip (which would create empty gaps between them).
    filter_geom = filter_polygon.geometry.iloc[0]
    fxmin, fymin, fxmax, fymax = filter_geom.bounds
    pad_x_left  = (fxmax - fxmin) * 0.03
    pad_x_right = (fxmax - fxmin) * 0.18
    pad_y = (fymax - fymin) * 0.02
    xmin, ymin = fxmin - pad_x_left,  fymin - pad_y
    xmax, ymax = fxmax + pad_x_right, fymax + pad_y

    # Use the PLACE_LABELS defined at the top of the file (config)
    place_labels = PLACE_LABELS

    def format_km(x, _pos):
        return f"{x/1000:.0f}"

    # ─── Style each axis ──────────────────────────────────────────────────
    for i, (ax, title) in enumerate(zip(axes, titles)):
        # No polygon backdrop is drawn — the plot extent alone reflects the
        # bounding-box limits set by the user's cascade_area.geojson.

        ax.set_xlim(xmin, xmax)
        ax.set_ylim(ymin, ymax)
        ax.set_aspect("equal", adjustable="box")

        # Axis formatting: km instead of meters
        ax.xaxis.set_major_formatter(FuncFormatter(format_km))
        ax.yaxis.set_major_formatter(FuncFormatter(format_km))
        ax.tick_params(labelsize=7)

        # Title with source color (pulls from SOURCE_COLORS for consistency
        # with the timeline plot)
        title_colors = [
            SOURCE_COLORS["wet_dry"],
            SOURCE_COLORS["nc_state"],
            SOURCE_COLORS["coastsat"],
            "#333333",
        ]
        ax.set_title(title, fontsize=11, pad=8,
                     color=title_colors[i], fontweight="bold")

        ax.set_xlabel("Easting (km, UTM 18N)", fontsize=8)
        ax.set_ylabel("Northing (km, UTM 18N)", fontsize=8)

        # Grid — subtle
        ax.grid(True, alpha=0.15, linewidth=0.4)

        # Location labels on the right side of every panel, in the reserved
        # right padding area (outside the island's footprint). Subtle white
        # background with rounded corners keeps text readable if any part
        # happens to fall over a shoreline feature.
        label_x = fxmax + pad_x_right * 0.35
        for name, y in place_labels:
            if ymin <= y <= ymax:
                # Small horizontal tick from island edge toward the label
                ax.plot([fxmax + pad_x_right * 0.05,
                         label_x - (fxmax - fxmin) * 0.01],
                        [y, y],
                        color="#666", lw=0.5, alpha=0.7, zorder=6)
                ax.text(label_x, y, name,
                        fontsize=7, ha="left", va="center",
                        color="#333", zorder=7,
                        bbox=dict(boxstyle="round,pad=0.20",
                                  facecolor="white", edgecolor="none",
                                  alpha=0.80))

    # ─── Basemap ──────────────────────────────────────────────────────────
    if basemap_kind != "none":
        print(f"    Adding basemap: {basemap_kind}")
        for ax in axes:
            _add_basemap_if_requested(ax, basemap_kind)

    # Helper to batch-plot a GeoDataFrame by unique year in one geopandas
    # call per year (much faster and avoids many aspect resets).
    def plot_gdf_by_year(gdf, ax, linewidth, linestyle="-", alpha=0.85,
                          zorder=3):
        if len(gdf) == 0:
            return
        for yr in sorted(gdf["_year"].unique()):
            subset = gdf[gdf["_year"] == yr]
            subset.plot(ax=ax, color=year_color(yr),
                        linewidth=linewidth, linestyle=linestyle,
                        alpha=alpha, zorder=zorder, aspect=None)

    # ─── Panel 1: wet-dry ─────────────────────────────────────────────────
    plot_gdf_by_year(wet_dry_gdf, axes[0],
                      linewidth=2.5, alpha=0.90, zorder=3)

    # ─── Panel 2: NC state (draw oldest first so newest on top) ───────────
    plot_gdf_by_year(nc_state_gdf, axes[1],
                      linewidth=1.6, alpha=0.80, zorder=3)

    # ─── Panel 3: CoastSat reconstructed polylines ────────────────────────
    if len(coastsat_recon) > 0 and len(coastsat_alongshore) > 0:
        alongshore_map = dict(zip(coastsat_alongshore["transect_id"],
                                   coastsat_alongshore["alongshore_m"]))
        recon = coastsat_recon.copy()
        recon["alongshore_m"] = recon["transect_id"].map(alongshore_map)
        recon = recon.dropna(subset=["alongshore_m"])
        for date, grp in sorted(recon.groupby("date"), key=lambda x: x[0]):
            year = pd.Timestamp(date).year
            grp_sorted = grp.sort_values("alongshore_m")
            axes[2].plot(grp_sorted["x"], grp_sorted["y"],
                          color=year_color(year), linewidth=1.8,
                          alpha=0.85, zorder=3)

    # ─── Panel 4: All sources overlaid, by linestyle ──────────────────────
    plot_gdf_by_year(nc_state_gdf, axes[3],
                      linewidth=1.3, linestyle="--", alpha=0.60, zorder=3)
    if len(coastsat_recon) > 0 and len(coastsat_alongshore) > 0:
        alongshore_map = dict(zip(coastsat_alongshore["transect_id"],
                                   coastsat_alongshore["alongshore_m"]))
        recon = coastsat_recon.copy()
        recon["alongshore_m"] = recon["transect_id"].map(alongshore_map)
        recon = recon.dropna(subset=["alongshore_m"])
        for date, grp in sorted(recon.groupby("date"), key=lambda x: x[0]):
            year = pd.Timestamp(date).year
            grp_sorted = grp.sort_values("alongshore_m")
            axes[3].plot(grp_sorted["x"], grp_sorted["y"],
                          color=year_color(year), linewidth=1.2,
                          linestyle=":", alpha=0.70, zorder=2)
    # Wet-dry drawn last (topmost) since it's the highest-precision source
    plot_gdf_by_year(wet_dry_gdf, axes[3],
                      linewidth=2.4, linestyle="-", alpha=0.90, zorder=4)

    style_legend = [
        Line2D([0], [0], color="#333", lw=2.4, linestyle="-",
               label="Wet-dry"),
        Line2D([0], [0], color="#333", lw=1.3, linestyle="--",
               label="NC state"),
        Line2D([0], [0], color="#333", lw=1.2, linestyle=":",
               label="CoastSat"),
    ]
    axes[3].legend(handles=style_legend, loc="upper left",
                    fontsize=8, framealpha=0.92,
                    edgecolor="#888", fancybox=False)

    # ─── Scale bar on rightmost panel ─────────────────────────────────────
    _draw_scale_bar(axes[3], xmin, ymin, xmax, ymax, length_km=5)

    # ─── Continuous year colorbar with tick marks at observed years ───────
    # Uses the rank-based color mapping (each unique year gets equal color
    # space) but presented as a smooth continuous bar rather than swatches.
    # Ticks are placed at every N-th observed year to avoid overlap.
    from matplotlib.colors import Normalize
    from matplotlib.cm import ScalarMappable

    cbar_ax = fig.add_axes([0.08, 0.07, 0.88, 0.020])
    norm = Normalize(vmin=0, vmax=max(n_years - 1, 1))
    sm = ScalarMappable(norm=norm, cmap=cmap_obj)
    sm.set_array([])
    cbar = fig.colorbar(sm, cax=cbar_ax, orientation="horizontal")
    cbar.outline.set_linewidth(0.4)
    cbar.outline.set_edgecolor("#888")

    # Choose how many year labels to show — avoid overlap
    n_desired_labels = 12
    step = max(1, n_years // n_desired_labels)
    tick_ranks  = list(range(0, n_years, step))
    if tick_ranks[-1] != n_years - 1:
        tick_ranks.append(n_years - 1)
    tick_labels = [str(unique_years[i]) for i in tick_ranks]

    cbar.set_ticks(tick_ranks)
    cbar.set_ticklabels(tick_labels)
    cbar.ax.tick_params(labelsize=8, length=4, color="#666")

    # Section title above colorbar
    fig.text(0.5, 0.115,
              f"Observation year  ({n_years} unique years, "
              f"rank-ordered on {colormap} colormap)",
              ha="center", va="bottom", fontsize=10, fontweight="bold")

    # ─── Suptitle ─────────────────────────────────────────────────────────
    fig.suptitle("Shoreline data coverage — Hatteras Island",
                 fontsize=14, fontweight="bold", y=0.965)

    plt.savefig(output_path, dpi=200, bbox_inches="tight")
    print(f"→ saved {output_path}")
    plt.close(fig)


# ============================================================
# VISUALIZATION: alongshore profile (CoastSat)
# ============================================================

def plot_shoreline_map_detailed(wet_dry_gdf: gpd.GeoDataFrame,
                                nc_state_gdf: gpd.GeoDataFrame,
                                coastsat_recon: pd.DataFrame,
                                coastsat_alongshore: pd.DataFrame,
                                filter_polygon: gpd.GeoDataFrame,
                                colormap: str,
                                basemap_kind: str,
                                output_path: str):
    """Standalone large-format all-sources overlaid map for detailed
    inspection. Same data as the 4-panel view but in one big panel where
    the shorelines are individually visible.

    Sources distinguished by line style + width:
      wet-dry   → thick solid   (highest precision, user-digitized)
      NC state  → medium dashed (historical survey data)
      CoastSat  → thin dotted   (satellite-derived, densely sampled)

    Year encoded by rank-based viridis coloring.
    """
    from matplotlib.gridspec import GridSpec
    from matplotlib.colors import Normalize
    from matplotlib.cm import ScalarMappable
    from matplotlib.lines import Line2D
    from matplotlib.ticker import FuncFormatter

    # Rank-based coloring across all sources
    all_years = []
    if len(wet_dry_gdf) > 0:
        all_years.extend(wet_dry_gdf["_year"].tolist())
    if len(nc_state_gdf) > 0:
        all_years.extend(nc_state_gdf["_year"].tolist())
    if len(coastsat_recon) > 0:
        all_years.extend(pd.to_datetime(coastsat_recon["date"]).dt.year.tolist())

    unique_years = sorted(set(int(y) for y in all_years))
    n_years = len(unique_years)
    if n_years == 0:
        return
    year_to_rank = {yr: i for i, yr in enumerate(unique_years)}
    cmap_obj = plt.get_cmap(colormap)

    def year_color(yr):
        rank = year_to_rank.get(int(yr))
        if rank is None:
            return (0.5, 0.5, 0.5, 1.0)
        return cmap_obj(rank / max(n_years - 1, 1))

    # Bounds — same as the 4-panel view (from the filter polygon)
    filter_geom = filter_polygon.geometry.iloc[0]
    fxmin, fymin, fxmax, fymax = filter_geom.bounds
    pad_x_right = (fxmax - fxmin) * 0.18
    pad_x_left  = (fxmax - fxmin) * 0.03
    pad_y = (fymax - fymin) * 0.02
    xmin, ymin = fxmin - pad_x_left,  fymin - pad_y
    xmax, ymax = fxmax + pad_x_right, fymax + pad_y

    # ─── Figure: tall single panel + right-side colorbar column ───────────
    fig = plt.figure(figsize=(10, 16))
    gs = GridSpec(1, 2, width_ratios=[10, 0.6], wspace=0.02,
                  top=0.94, bottom=0.06, left=0.10, right=0.95)
    ax = fig.add_subplot(gs[0, 0])
    cbar_ax = fig.add_subplot(gs[0, 1])

    def format_km(x, _pos):
        return f"{x/1000:.0f}"

    ax.set_xlim(xmin, xmax)
    ax.set_ylim(ymin, ymax)
    ax.set_aspect("equal", adjustable="box")
    ax.xaxis.set_major_formatter(FuncFormatter(format_km))
    ax.yaxis.set_major_formatter(FuncFormatter(format_km))
    ax.set_xlabel("Easting (km, UTM 18N)", fontsize=10)
    ax.set_ylabel("Northing (km, UTM 18N)", fontsize=10)
    ax.tick_params(labelsize=8)
    ax.grid(True, alpha=0.15, linewidth=0.4)

    # Basemap (if requested)
    if basemap_kind != "none":
        print(f"    Adding basemap: {basemap_kind}")
        _add_basemap_if_requested(ax, basemap_kind)

    # ─── Plot each source with distinct style ─────────────────────────────
    # NC state — medium dashed (drawn oldest first so newer on top)
    if len(nc_state_gdf) > 0:
        for yr in sorted(nc_state_gdf["_year"].unique()):
            subset = nc_state_gdf[nc_state_gdf["_year"] == yr]
            subset.plot(ax=ax, color=year_color(yr),
                        linewidth=1.4, linestyle="--",
                        alpha=0.75, zorder=3, aspect=None)

    # CoastSat — thin dotted lines
    if len(coastsat_recon) > 0 and len(coastsat_alongshore) > 0:
        alongshore_map = dict(zip(coastsat_alongshore["transect_id"],
                                   coastsat_alongshore["alongshore_m"]))
        recon = coastsat_recon.copy()
        recon["alongshore_m"] = recon["transect_id"].map(alongshore_map)
        recon = recon.dropna(subset=["alongshore_m"])
        for date, grp in sorted(recon.groupby("date"), key=lambda x: x[0]):
            year = pd.Timestamp(date).year
            grp_sorted = grp.sort_values("alongshore_m")
            ax.plot(grp_sorted["x"], grp_sorted["y"],
                    color=year_color(year), linewidth=1.1,
                    linestyle=":", alpha=0.75, zorder=2)

    # Wet-dry — thick solid, drawn on top (highest precision source)
    if len(wet_dry_gdf) > 0:
        for yr in sorted(wet_dry_gdf["_year"].unique()):
            subset = wet_dry_gdf[wet_dry_gdf["_year"] == yr]
            subset.plot(ax=ax, color=year_color(yr),
                        linewidth=2.8, linestyle="-",
                        alpha=0.95, zorder=4, aspect=None)

    # ─── Location labels on the right side ────────────────────────────────
    label_x = fxmax + pad_x_right * 0.35
    for name, y in PLACE_LABELS:
        if ymin <= y <= ymax:
            ax.plot([fxmax + pad_x_right * 0.05,
                     label_x - (fxmax - fxmin) * 0.01],
                    [y, y],
                    color="#666", lw=0.6, alpha=0.7, zorder=6)
            ax.text(label_x, y, name,
                    fontsize=9, ha="left", va="center",
                    color="#222", zorder=7,
                    bbox=dict(boxstyle="round,pad=0.25",
                              facecolor="white", edgecolor="none",
                              alpha=0.85))

    # ─── Scale bar in lower-left ──────────────────────────────────────────
    _draw_scale_bar(ax, xmin, ymin, xmax, ymax, length_km=5)

    # ─── Source legend (top-left) ─────────────────────────────────────────
    source_legend = [
        Line2D([0], [0], color="#333", lw=2.8, linestyle="-",
               label="Aerial wet-dry lines"),
        Line2D([0], [0], color="#333", lw=1.4, linestyle="--",
               label="NC Coastal Mgmt"),
        Line2D([0], [0], color="#333", lw=1.1, linestyle=":",
               label="CoastSat (satellite)"),
    ]
    ax.legend(handles=source_legend, loc="upper left",
              fontsize=10, framealpha=0.92,
              edgecolor="#888", fancybox=False,
              title="Source", title_fontsize=10)

    # ─── Vertical continuous colorbar ─────────────────────────────────────
    norm = Normalize(vmin=0, vmax=max(n_years - 1, 1))
    sm = ScalarMappable(norm=norm, cmap=cmap_obj)
    sm.set_array([])
    cbar = fig.colorbar(sm, cax=cbar_ax, orientation="vertical")
    cbar.outline.set_linewidth(0.4)
    cbar.outline.set_edgecolor("#888")

    # Ticks at every N-th observed year to avoid clutter
    n_desired_labels = 15
    step = max(1, n_years // n_desired_labels)
    tick_ranks = list(range(0, n_years, step))
    if tick_ranks[-1] != n_years - 1:
        tick_ranks.append(n_years - 1)
    cbar.set_ticks(tick_ranks)
    cbar.set_ticklabels([str(unique_years[i]) for i in tick_ranks])
    cbar.ax.tick_params(labelsize=8, length=3, color="#666")
    cbar.set_label(f"Year  (rank-ordered, {n_years} unique years)",
                    fontsize=10, labelpad=8)

    # ─── Suptitle ─────────────────────────────────────────────────────────
    fig.suptitle("All shoreline sources — Hatteras Island (detailed view)",
                 fontsize=13, fontweight="bold", y=0.98)

    plt.savefig(output_path, dpi=200, bbox_inches="tight")
    print(f"→ saved {output_path}")
    plt.close(fig)


def plot_alongshore_profile(coastsat_recon: pd.DataFrame,
                             coastsat_alongshore: pd.DataFrame,
                             colormap: str,
                             output_path: str):
    """CoastSat shoreline chainage vs. alongshore position, colored by year.

    X-axis: alongshore position (km, south to north)
    Y-axis: chainage (m from transect origin, positive = seaward)
    """
    from matplotlib.colors import Normalize
    from matplotlib.cm import ScalarMappable

    fig, ax = plt.subplots(figsize=(15, 7))

    if len(coastsat_recon) == 0 or len(coastsat_alongshore) == 0:
        ax.text(0.5, 0.5, "No CoastSat reconstructed shorelines available",
                transform=ax.transAxes, ha="center", va="center")
        plt.savefig(output_path, dpi=150, bbox_inches="tight")
        plt.close(fig)
        return

    alongshore_map = dict(zip(coastsat_alongshore["transect_id"],
                               coastsat_alongshore["alongshore_m"]))
    recon = coastsat_recon.copy()
    recon["alongshore_m"] = recon["transect_id"].map(alongshore_map)
    recon = recon.dropna(subset=["alongshore_m"])
    recon["alongshore_km"] = recon["alongshore_m"] / 1000.0

    years = pd.to_datetime(recon["date"]).dt.year
    y_min, y_max = years.min(), years.max()
    norm = Normalize(vmin=y_min, vmax=y_max)
    cmap = plt.get_cmap(colormap)

    # Plot each date as a line, sorted by alongshore position
    for date, grp in recon.groupby("date"):
        year = pd.Timestamp(date).year
        color = cmap(norm(year))
        grp_sorted = grp.sort_values("alongshore_m")
        ax.plot(grp_sorted["alongshore_km"], grp_sorted["chainage_m"],
                color=color, lw=1.4, alpha=0.85,
                label=pd.Timestamp(date).strftime("%Y-%m-%d"))

    ax.set_xlabel("Alongshore position (km, south → north)", fontsize=11)
    ax.set_ylabel("Chainage from transect origin (m, positive = seaward)",
                  fontsize=11)
    ax.set_title(
        f"CoastSat reconstructed shorelines by alongshore position  "
        f"({len(recon.groupby('date'))} snapshot dates)",
        fontsize=12,
    )
    ax.grid(True, alpha=0.3)

    # Colorbar
    sm = ScalarMappable(norm=norm, cmap=cmap)
    sm.set_array([])
    cbar = fig.colorbar(sm, ax=ax, orientation="vertical", pad=0.02,
                        fraction=0.03)
    cbar.set_label("Year", fontsize=10)

    # Small legend of the snapshot dates
    ax.legend(loc="upper right", fontsize=7, ncol=2, framealpha=0.85,
              title="Snapshot dates")

    plt.tight_layout()
    plt.savefig(output_path, dpi=150, bbox_inches="tight")
    print(f"→ saved {output_path}")
    plt.close(fig)


# ============================================================
# VISUALIZATION: helper to reload shapefile geoms for plotting
# ============================================================

def load_shapefile_geoms_for_plot(path: str, date_col: str,
                                   filter_polygon: gpd.GeoDataFrame) -> gpd.GeoDataFrame:
    """Reload a shapefile and return the filtered geometries + parsed year.
    Used by the map-view plotting step (separate from the inventory step)."""
    if not os.path.exists(path):
        return gpd.GeoDataFrame(columns=["_year", "geometry"],
                                 geometry="geometry", crs=PROJECTED_CRS)
    gdf = gpd.read_file(path)
    gdf = fix_crs(gdf)
    if date_col not in gdf.columns:
        return gpd.GeoDataFrame(columns=["_year", "geometry"],
                                 geometry="geometry", crs=PROJECTED_CRS)
    gdf = gdf.to_crs(PROJECTED_CRS)
    filter_geom = filter_polygon.geometry.iloc[0]
    gdf = gdf[gdf.geometry.intersects(filter_geom)].copy()
    if len(gdf) == 0:
        return gpd.GeoDataFrame(columns=["_year", "geometry"],
                                 geometry="geometry", crs=PROJECTED_CRS)
    gdf["_parsed_date"] = parse_date_column(gdf[date_col], col_name=date_col)
    gdf = gdf.dropna(subset=["_parsed_date"])
    gdf["_year"] = gdf["_parsed_date"].dt.year
    return gdf


# ============================================================
# MAIN
# ============================================================

def main():
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    print("=" * 72)
    print("SHORELINE DATA INVENTORY — Hatteras Island / Buxton study area")
    print("=" * 72 + "\n")

    # ─── Load spatial filter ───────────────────────────────────────
    filter_polygon = load_study_area_filter(
        STUDY_AREA_FILTER_PATH, STUDY_AREA_BUFFER_M,
    )

    # ─── Load each source ──────────────────────────────────────────
    wet_dry_obs  = load_shapefile_dates(
        WET_DRY_PATH, WET_DRY_DATE_COL,
        filter_polygon, source_key="wet_dry",
    )
    nc_state_obs = load_shapefile_dates(
        NC_STATE_PATH, NC_STATE_DATE_COL,
        filter_polygon, source_key="nc_state",
    )
    coastsat_obs, coastsat_raw, coastsat_transects = load_coastsat_dates(
        COASTSAT_TRANSECT_GEOM, COASTSAT_TRANSECT_ID_COL,
        COASTSAT_ROOT_DIR, filter_polygon,
        COASTSAT_WELL_COVERED_MIN_TRANSECTS,
        source_key="coastsat",
    )

    all_obs = pd.concat([wet_dry_obs, nc_state_obs, coastsat_obs],
                        ignore_index=True)

    if len(all_obs) == 0:
        print("No observations found across any source. Check config paths.")
        return

    print("=" * 72)
    print("AGGREGATE SUMMARY")
    print("=" * 72)
    for src in SOURCE_ORDER:
        n = (all_obs["source"] == src).sum()
        n_years = all_obs.loc[all_obs["source"] == src, "year"].nunique()
        if n > 0:
            span = (all_obs.loc[all_obs["source"] == src, "year"].min(),
                    all_obs.loc[all_obs["source"] == src, "year"].max())
            print(f"  {SOURCE_LABELS[src]:22s}: {n:5d} observations, "
                  f"{n_years:3d} unique years  ({span[0]}–{span[1]})")
        else:
            print(f"  {SOURCE_LABELS[src]:22s}: 0 observations")

    print(f"\n  Total observations across sources: {len(all_obs)}")
    print(f"  Full span: {all_obs['year'].min()}–{all_obs['year'].max()}")

    # ─── Build year summary ───────────────────────────────────────
    year_summary = build_year_summary(all_obs,
                                       COASTSAT_WELL_COVERED_MIN_TRANSECTS)

    # ─── Save outputs ─────────────────────────────────────────────
    obs_out  = os.path.join(OUTPUT_DIR, "shoreline_inventory_observations.csv")
    year_out = os.path.join(OUTPUT_DIR, "shoreline_inventory_by_year.csv")

    all_obs_sorted = all_obs.sort_values(["date", "source"]).reset_index(drop=True)
    all_obs_sorted.to_csv(obs_out, index=False)
    year_summary.to_csv(year_out, index=False)
    print(f"\n→ saved {obs_out}")
    print(f"→ saved {year_out}")

    plot_out = os.path.join(OUTPUT_DIR, "shoreline_inventory_timeline.png")
    plot_inventory(all_obs, year_summary,
                    COASTSAT_WELL_COVERED_MIN_TRANSECTS, plot_out)

    # ─── Spatial visualizations ───────────────────────────────────
    if MAKE_VISUALIZATIONS:
        print("\n" + "=" * 72)
        print("VISUALIZATION")
        print("=" * 72)

        # Reload shapefile geometries with parsed year (used by the map view).
        # Cheap since these files are small.
        print("\nReloading shapefile geometries for map plotting...")
        wet_dry_gdf  = load_shapefile_geoms_for_plot(
            WET_DRY_PATH, WET_DRY_DATE_COL, filter_polygon,
        )
        nc_state_gdf = load_shapefile_geoms_for_plot(
            NC_STATE_PATH, NC_STATE_DATE_COL, filter_polygon,
        )
        print(f"  Wet-dry features to plot: {len(wet_dry_gdf)}")
        print(f"  NC state features to plot: {len(nc_state_gdf)}")

        # Pick representative CoastSat snapshot dates
        print(f"\nPicking {N_COASTSAT_SNAPSHOTS} representative CoastSat dates...")
        snapshot_dates = pick_representative_coastsat_dates(
            coastsat_obs, N_COASTSAT_SNAPSHOTS,
            COASTSAT_WELL_COVERED_MIN_TRANSECTS,
        )
        print(f"  Selected {len(snapshot_dates)} snapshot dates:")
        for _, row in snapshot_dates.iterrows():
            print(f"    {pd.Timestamp(row['date']).strftime('%Y-%m-%d')} "
                  f"— {int(row['n_transects_covered']):3d} transects")

        # Reconstruct spatial coordinates of CoastSat shorelines at those dates
        print(f"\nReconstructing CoastSat shoreline points...")
        coastsat_recon = reconstruct_coastsat_shoreline_points(
            coastsat_raw, coastsat_transects, snapshot_dates["date"],
        )
        print(f"  Reconstructed {len(coastsat_recon)} points across "
              f"{coastsat_recon['date'].nunique()} dates")

        # Compute alongshore position for each transect
        print("\nComputing alongshore positions for CoastSat transects...")
        coastsat_alongshore = compute_alongshore_position(coastsat_transects)
        if len(coastsat_alongshore) > 0:
            print(f"  Alongshore extent: 0 → "
                  f"{coastsat_alongshore['alongshore_m'].max() / 1000:.1f} km")

        # Map view
        print("\nGenerating map view...")
        map_out = os.path.join(OUTPUT_DIR, "shoreline_inventory_map.png")
        plot_shoreline_map(
            wet_dry_gdf, nc_state_gdf,
            coastsat_recon, coastsat_alongshore,
            filter_polygon, SHORELINE_COLORMAP, BASEMAP, map_out,
        )

        # Detailed all-sources overlaid map (standalone large-format figure)
        print("\nGenerating detailed all-sources map...")
        detailed_out = os.path.join(OUTPUT_DIR,
                                      "shoreline_inventory_map_detailed.png")
        plot_shoreline_map_detailed(
            wet_dry_gdf, nc_state_gdf,
            coastsat_recon, coastsat_alongshore,
            filter_polygon, SHORELINE_COLORMAP, BASEMAP, detailed_out,
        )

        # Alongshore profile (CoastSat only)
        print("\nGenerating alongshore profile view (CoastSat)...")
        profile_out = os.path.join(OUTPUT_DIR,
                                    "shoreline_inventory_alongshore_profile.png")
        plot_alongshore_profile(
            coastsat_recon, coastsat_alongshore,
            SHORELINE_COLORMAP, profile_out,
        )

    # ─── Print the year summary to console ────────────────────────
    print("\n" + "=" * 72)
    print("YEAR-BY-YEAR SUMMARY (non-zero years only)")
    print("=" * 72)
    non_zero = year_summary[year_summary["total_all_sources"] > 0].copy()
    with pd.option_context("display.max_rows", None):
        print(non_zero.to_string(index=False))

    # ─── Highlight gaps ───────────────────────────────────────────
    print("\n" + "=" * 72)
    print("GAP ANALYSIS")
    print("=" * 72)
    if len(non_zero) > 0:
        years = sorted(non_zero["year"].values)
        gaps = []
        for i in range(len(years) - 1):
            gap = years[i + 1] - years[i]
            if gap > 3:
                gaps.append((years[i], years[i + 1], gap))
        if gaps:
            print(f"Gaps > 3 years with no data from any source:")
            for start, end, gap in gaps:
                print(f"  {start} → {end}  ({gap} years, no observations)")
        else:
            print("No gaps > 3 years — good temporal coverage across sources")

    print("\n" + "=" * 72)
    print("Done. Files saved to:")
    print(f"  {OUTPUT_DIR}")
    print("=" * 72)


if __name__ == "__main__":
    main()
