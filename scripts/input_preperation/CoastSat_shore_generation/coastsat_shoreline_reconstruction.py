"""
CoastSat Shoreline Reconstruction
===================================
Generates interpolated shorelines at regular time steps from CoastSat
time-series data and exports them as shapefiles and GeoJSON for GIS.

Method
------
For each target date:
  1. For each transect, find the two nearest observations that bracket
     the target date (one before, one after).
  2. Linearly interpolate the cross-shore distance to the target date.
  3. Convert that distance to a geographic coordinate using the transect
     origin point and bearing from the CoastSat GeoJSON.
  4. Connect all valid transect points into a polyline.
  5. Export as shapefile + GeoJSON — one file per shoreline date.

Transects with no observations bracketing a given date are skipped
(no extrapolation). A summary of coverage per shoreline is printed.

After the run, two diagnostic plots open automatically:
  1. Coverage bar chart — % of transects successfully interpolated per date
  2. Shoreline map — all written shorelines coloured by date

You can also call plot_transect_timeseries() interactively from the
PyCharm Python console to inspect any individual transect.

Dependencies
------------
  pip install geopandas pandas numpy shapely pyproj matplotlib
"""

# ============================================================
# CONFIG  -  edit before running
# ============================================================

# CoastSat transects GeoJSON - needed for origin coordinates and bearings
TRANSECT_GEOM_PATH = r"C:\Users\hanna\PycharmProjects\CASCADE\scripts\input_preperation\CoastSat\CoastSat_transect_layer.geojson"
TRANSECT_ID_COL    = "id"

# Time-series CSVs - same setup as coastsat_domain_lrr.py
ROOT_DATA_DIR = r"C:\Users\hanna\PycharmProjects\CASCADE\scripts\input_preperation\CoastSat\coastsat_timeseries"
SITE_FILTER   = "usa_NC"

# Lookup table from coastsat_domain_mapping.py
# Used to restrict to Hatteras transects and order them along-shore.
# Set to None to use all transects found in ROOT_DATA_DIR.
LOOKUP_CSV = r"C:\Users\hanna\PycharmProjects\CASCADE\scripts\input_preperation\CoastSat\transect_domain_lookup.csv"

# Date range for shoreline generation (inclusive)
START_DATE = "2024-01-01"
END_DATE   = "2024-06-01"

# Time step between generated shorelines
# "MS" = monthly, "SMS" = semi-monthly (~every 2 weeks), "YS" = yearly
# SMS recommended for modern data (post-2015) given ~15 day median observation gap
DATE_FREQ = "SMS"

# Maximum gap allowed between bracketing observations (days).
# Recommended: 30 days for modern data (SMS), 365 days for early Landsat era
MAX_INTERP_GAP_DAYS = 30

# Minimum number of transects required to write a shoreline.
# ~900 total Hatteras transects; 400 = ~44% minimum coverage
MIN_TRANSECTS_PER_SHORELINE = 400

# Domain range - only include transects from your 90 real domains.
# Set to None to include all matched transects.
DOMAIN_MIN = 1
DOMAIN_MAX = 90

# Output directory
OUTPUT_DIR = r"C:\Users\hanna\PycharmProjects\CASCADE\scripts\input_preperation\CoastSat_shore_generation\CoastSat_shoreline_outputs"

# Coordinate Reference System for output shapefiles
# UTM Zone 18N is appropriate for the NC Outer Banks
OUTPUT_CRS = "EPSG:32618"

# ============================================================
# IMPORTS
# ============================================================
import os
import glob
import warnings
import numpy as np
import pandas as pd
import geopandas as gpd
from shapely.geometry import LineString
import matplotlib.pyplot as plt
import matplotlib.cm as cm
warnings.filterwarnings("ignore")


# ============================================================
# STEP 1: LOAD TRANSECT GEOMETRY
# ============================================================

def load_transect_geometry(geojson_path, id_col):
    """
    Load transect origin coordinates and bearings from CoastSat GeoJSON.
    Returns a DataFrame indexed by transect ID with columns:
      origin_lon, origin_lat, bearing_deg
    """
    gdf = gpd.read_file(geojson_path)
    print(f"Loaded {len(gdf):,} transects from GeoJSON")

    if gdf.crs is None or gdf.crs.to_epsg() != 4326:
        gdf = gdf.set_crs("EPSG:4326", allow_override=True)

    records = []
    for _, row in gdf.iterrows():
        geom = row.geometry
        if geom is None or geom.is_empty:
            continue
        if geom.geom_type == "LineString":
            coords = list(geom.coords)
        elif geom.geom_type == "MultiLineString":
            coords = list(geom.geoms[0].coords)
        else:
            continue
        if len(coords) < 2:
            continue

        ox, oy = coords[0]
        ex, ey = coords[-1]

        d_lon   = np.radians(ex - ox)
        lat1    = np.radians(oy)
        lat2    = np.radians(ey)
        x       = np.sin(d_lon) * np.cos(lat2)
        y       = np.cos(lat1) * np.sin(lat2) - np.sin(lat1) * np.cos(lat2) * np.cos(d_lon)
        bearing = (np.degrees(np.arctan2(x, y)) + 360) % 360

        records.append({
            "transect_id": str(row[id_col]),
            "origin_lon" : ox,
            "origin_lat" : oy,
            "bearing_deg": bearing,
        })

    geom_df = pd.DataFrame(records).set_index("transect_id")
    print(f"  Extracted geometry for {len(geom_df):,} transects")
    return geom_df


# ============================================================
# STEP 2: LOAD TIME-SERIES CSVs
# ============================================================

def collect_csv_map(root_dir, site_filter=""):
    """Auto-discover all time-series CSVs one level below root_dir."""
    csv_map = {}
    if not os.path.isdir(root_dir):
        print(f"WARNING: ROOT_DATA_DIR not found: {root_dir}")
        return csv_map

    subfolders = [
        os.path.join(root_dir, d)
        for d in sorted(os.listdir(root_dir))
        if os.path.isdir(os.path.join(root_dir, d))
        and (site_filter == "" or site_filter in d)
    ]
    for sf in subfolders:
        for fpath in glob.glob(os.path.join(sf, "*.csv")):
            stem = os.path.splitext(os.path.basename(fpath))[0]
            csv_map[stem] = fpath

    print(f"Found {len(csv_map):,} time-series CSVs across {len(subfolders)} site folders")
    return csv_map


def load_timeseries(filepath):
    """Load a single CoastSat time-series CSV into a DataFrame."""
    df = pd.read_csv(filepath, header=0)
    df.columns = [c.strip() for c in df.columns]
    df = df.rename(columns={df.columns[0]: "date", df.columns[1]: "chainage_m"})
    df["date"]       = pd.to_datetime(df["date"], utc=True)
    df["chainage_m"] = pd.to_numeric(df["chainage_m"], errors="coerce")
    df = df.dropna().sort_values("date").reset_index(drop=True)
    return df


def load_all_timeseries(csv_map, valid_ids):
    """Load time-series for all transects in valid_ids."""
    ts_data = {}
    missing = []
    for tid in sorted(valid_ids):
        if tid not in csv_map:
            missing.append(tid)
            continue
        try:
            ts_data[tid] = load_timeseries(csv_map[tid])
        except Exception as e:
            print(f"  ERROR loading {tid}: {e}")
    if missing:
        print(f"  WARNING: {len(missing)} transects in lookup not found in CSVs")
    print(f"  Loaded time-series for {len(ts_data):,} transects")
    return ts_data


# ============================================================
# STEP 3: INTERPOLATE POSITION AT TARGET DATE
# ============================================================

def interpolate_chainage(ts, target, max_gap_days):
    """
    Linearly interpolate cross-shore chainage for a target date.
    Returns None if no bracketing observations exist or gap exceeds limit.
    """
    before = ts[ts["date"] <= target]
    after  = ts[ts["date"] >  target]

    if before.empty or after.empty:
        return None

    t_before = before.iloc[-1]
    t_after  = after.iloc[0]

    gap_days = (t_after["date"] - t_before["date"]).days
    if gap_days > max_gap_days:
        return None

    frac = ((target - t_before["date"]).total_seconds() /
            (t_after["date"] - t_before["date"]).total_seconds())
    return t_before["chainage_m"] + frac * (t_after["chainage_m"] - t_before["chainage_m"])


# ============================================================
# STEP 4: CONVERT CHAINAGE TO GEOGRAPHIC COORDINATE
# ============================================================

def chainage_to_point(origin_lon, origin_lat, bearing_deg, distance_m):
    """Convert transect origin + bearing + distance to (lon, lat)."""
    R           = 6371000.0
    bearing_rad = np.radians(bearing_deg)
    lat1        = np.radians(origin_lat)
    lon1        = np.radians(origin_lon)
    d_R         = distance_m / R

    lat2 = np.arcsin(np.sin(lat1) * np.cos(d_R) +
                     np.cos(lat1) * np.sin(d_R) * np.cos(bearing_rad))
    lon2 = lon1 + np.arctan2(
        np.sin(bearing_rad) * np.sin(d_R) * np.cos(lat1),
        np.cos(d_R) - np.sin(lat1) * np.sin(lat2)
    )
    return (np.degrees(lon2), np.degrees(lat2))


# ============================================================
# STEP 5: BUILD AND EXPORT SHORELINES
# ============================================================

def build_shoreline_for_date(target_date, ordered_transects, ts_data,
                              geom_df, max_gap_days):
    """Interpolate shoreline position for all transects at a target date."""
    points  = []
    n_valid = 0
    n_skip  = 0

    for tid in ordered_transects:
        if tid not in ts_data or tid not in geom_df.index:
            n_skip += 1
            continue
        chainage = interpolate_chainage(ts_data[tid], target_date, max_gap_days)
        if chainage is None:
            n_skip += 1
            continue
        row = geom_df.loc[tid]
        lon, lat = chainage_to_point(
            row["origin_lon"], row["origin_lat"], row["bearing_deg"], chainage
        )
        points.append((lon, lat))
        n_valid += 1

    return points, n_valid, n_skip


def export_shoreline(points, target_date, shp_dir, geojson_dir, output_crs):
    """
    Export shoreline as shapefile (projected) and GeoJSON (WGS84).
    GeoJSON is directly loadable into DSAS.
    """
    if len(points) < 2:
        return

    fname_stem = "shoreline_" + target_date.strftime("%Y_%m_%d")
    attrs = {
        "date"    : [target_date.strftime("%Y-%m-%d")],
        "year"    : [target_date.year],
        "month"   : [target_date.month],
        "geometry": [LineString(points)],
    }

    gpd.GeoDataFrame(attrs, crs="EPSG:4326").to_crs(output_crs).to_file(
        os.path.join(shp_dir, fname_stem + ".shp")
    )
    gpd.GeoDataFrame(attrs, crs="EPSG:4326").to_file(
        os.path.join(geojson_dir, fname_stem + ".geojson"),
        driver="GeoJSON"
    )


# ============================================================
# VISUALISATION
# ============================================================

def plot_shoreline_coverage(coverage_csv):
    """
    Bar chart of transect coverage % per date.
    Blue >50%, orange 25-50%, red <25%.
    """
    df = pd.read_csv(coverage_csv)
    df["date"] = pd.to_datetime(df["date"])

    colors = [
        "steelblue" if p >= 50 else "#f4a261" if p >= 25 else "#e76f51"
        for p in df["pct_coverage"]
    ]

    fig, ax = plt.subplots(figsize=(16, 4))
    ax.bar(range(len(df)), df["pct_coverage"], color=colors, edgecolor="none")
    ax.axhline(50, color="black", lw=1, ls="--", alpha=0.5, label="50% threshold")

    tick_step = max(1, len(df) // 20)
    ax.set_xticks(range(0, len(df), tick_step))
    ax.set_xticklabels(
        [df["date"].iloc[i].strftime("%Y-%m") for i in range(0, len(df), tick_step)],
        rotation=45, ha="right", fontsize=8
    )
    ax.set_ylabel("% Transects Covered", fontsize=11)
    ax.set_title("CoastSat Shoreline Coverage by Date", fontsize=13, fontweight="bold")
    ax.legend(fontsize=9)
    ax.set_ylim(0, 105)
    ax.grid(True, axis="y", alpha=0.3)
    plt.tight_layout()
    plt.show()


def plot_shorelines_map(shoreline_records, max_lines=24):
    """
    Map of generated shorelines coloured by date (plasma colormap).
    shoreline_records: list of (target_date, points) tuples from main().
    """
    if not shoreline_records:
        print("No shoreline records to plot.")
        return

    step   = max(1, len(shoreline_records) // max_lines)
    subset = shoreline_records[::step]
    n      = len(subset)
    cmap   = cm.plasma

    fig, ax = plt.subplots(figsize=(6, 14))
    for i, (date, points) in enumerate(subset):
        if len(points) < 2:
            continue
        lons = [p[0] for p in points]
        lats = [p[1] for p in points]
        ax.plot(lons, lats, color=cmap(i / max(n - 1, 1)), lw=1.2, alpha=0.8)

    sm = cm.ScalarMappable(cmap=cmap, norm=plt.Normalize(0, n - 1))
    sm.set_array([])
    cbar = plt.colorbar(sm, ax=ax, fraction=0.03, pad=0.02)
    cbar.set_ticks([0, n - 1])
    cbar.set_ticklabels([
        subset[0][0].strftime("%Y-%m"),
        subset[-1][0].strftime("%Y-%m"),
    ])
    cbar.set_label("Date", fontsize=10)

    ax.set_xlabel("Longitude", fontsize=10)
    ax.set_ylabel("Latitude", fontsize=10)
    ax.set_title(
        "Reconstructed Shorelines\n"
        + subset[0][0].strftime("%Y-%m") + " to " + subset[-1][0].strftime("%Y-%m"),
        fontsize=12, fontweight="bold"
    )
    ax.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.show()


def plot_transect_timeseries(transect_id, ts_data, start=None, end=None):
    """
    Plot the raw CoastSat time series for a single transect.
    Call interactively from PyCharm Python console after main() finishes:

        plot_transect_timeseries("usa_NC_0033_0045", ts_data,
                                  start="2020-01-01", end="2025-01-01")
    """
    if transect_id not in ts_data:
        print(f"Transect {transect_id} not found in loaded data.")
        return

    df = ts_data[transect_id].copy()
    if start:
        df = df[df["date"] >= pd.Timestamp(start, tz="UTC")]
    if end:
        df = df[df["date"] <= pd.Timestamp(end, tz="UTC")]

    fig, ax = plt.subplots(figsize=(12, 4))
    ax.scatter(df["date"], df["chainage_m"], color="steelblue", s=20, zorder=3,
               label="Observations")
    ax.plot(df["date"], df["chainage_m"], color="steelblue", lw=0.8, alpha=0.5)
    ax.axhline(df["chainage_m"].mean(), color="grey", lw=1, ls="--",
               label=f"Mean ({df['chainage_m'].mean():.1f} m)")
    ax.set_xlabel("Date", fontsize=11)
    ax.set_ylabel("Cross-shore distance (m)", fontsize=11)
    ax.set_title(f"CoastSat Time Series - {transect_id}",
                 fontsize=12, fontweight="bold")
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.show()


# ============================================================
# MAIN
# ============================================================

def main():
    print("=" * 65)
    print("CoastSat Shoreline Reconstruction")
    print("=" * 65)

    os.makedirs(OUTPUT_DIR, exist_ok=True)
    shp_dir     = os.path.join(OUTPUT_DIR, "shorelines")
    geojson_dir = os.path.join(OUTPUT_DIR, "shorelines_geojson")
    os.makedirs(shp_dir,     exist_ok=True)
    os.makedirs(geojson_dir, exist_ok=True)

    # Load transect geometry
    print("\nLoading transect geometry...")
    geom_df = load_transect_geometry(TRANSECT_GEOM_PATH, TRANSECT_ID_COL)

    # Load lookup table
    if LOOKUP_CSV and os.path.exists(LOOKUP_CSV):
        print("Loading domain lookup table...")
        lookup = pd.read_csv(LOOKUP_CSV)
        lookup["transect_id"]   = lookup["transect_id"].astype(str)
        lookup["domain_number"] = pd.to_numeric(lookup["domain_number"], errors="coerce")
        if DOMAIN_MIN is not None and DOMAIN_MAX is not None:
            lookup = lookup[(lookup["domain_number"] >= DOMAIN_MIN) &
                            (lookup["domain_number"] <= DOMAIN_MAX)]
        lookup            = lookup.sort_values("domain_number")
        valid_ids         = set(lookup["transect_id"].tolist())
        ordered_transects = lookup["transect_id"].tolist()
        print(f"  Using {len(valid_ids):,} transects from domains {DOMAIN_MIN}-{DOMAIN_MAX}")
    else:
        print("No lookup table - using all transects found in CSVs")
        csv_map_temp      = collect_csv_map(ROOT_DATA_DIR, SITE_FILTER)
        valid_ids         = set(csv_map_temp.keys())
        ordered_transects = sorted(valid_ids)

    # Load time-series
    print("\nLoading time-series CSVs...")
    csv_map = collect_csv_map(ROOT_DATA_DIR, SITE_FILTER)
    ts_data = load_all_timeseries(csv_map, valid_ids)

    # Generate target dates
    target_dates = pd.date_range(
        start=START_DATE, end=END_DATE, freq=DATE_FREQ, tz="UTC"
    )
    print(f"\nGenerating {len(target_dates)} shorelines "
          f"({START_DATE} to {END_DATE}, freq={DATE_FREQ})")
    print(f"  Max interpolation gap : {MAX_INTERP_GAP_DAYS} days")
    print(f"  Min transects         : {MIN_TRANSECTS_PER_SHORELINE}")
    print("-" * 65)

    # Generate and export
    coverage_records  = []
    shoreline_records = []
    n_written = 0
    n_skipped = 0

    for target in target_dates:
        points, n_valid, n_skip = build_shoreline_for_date(
            target, ordered_transects, ts_data, geom_df, MAX_INTERP_GAP_DAYS
        )

        coverage_records.append({
            "date"        : target.strftime("%Y-%m-%d"),
            "n_valid"     : n_valid,
            "n_skipped"   : n_skip,
            "pct_coverage": round(n_valid / max(n_valid + n_skip, 1) * 100, 1),
        })

        if n_valid < MIN_TRANSECTS_PER_SHORELINE:
            print(f"  SKIP  {target.strftime('%Y-%m-%d')}  "
                  f"({n_valid} transects < minimum {MIN_TRANSECTS_PER_SHORELINE})")
            n_skipped += 1
            continue

        export_shoreline(points, target, shp_dir, geojson_dir, OUTPUT_CRS)
        shoreline_records.append((target, points))
        print(f"  OK    {target.strftime('%Y-%m-%d')}  "
              f"({n_valid}/{n_valid + n_skip} transects, "
              f"{round(n_valid / (n_valid + n_skip) * 100)}% coverage)")
        n_written += 1

    # Save coverage summary
    cov_path = os.path.join(OUTPUT_DIR, "shoreline_coverage.csv")
    pd.DataFrame(coverage_records).to_csv(cov_path, index=False)

    print(f"\n{'=' * 65}")
    print(f"  Shorelines written : {n_written}")
    print(f"  Dates skipped      : {n_skipped} (too few transects)")
    print(f"  Shapefiles         : {shp_dir}")
    print(f"  GeoJSON files      : {geojson_dir}")
    print(f"  Coverage summary   : {cov_path}")
    print(f"{'=' * 65}")

    # Diagnostic plots
    print("\nGenerating diagnostic plots...")
    plot_shoreline_coverage(cov_path)
    if shoreline_records:
        plot_shorelines_map(shoreline_records)
    else:
        print("  No shorelines written - nothing to map.")

    # ts_data returned for interactive use in PyCharm console:
    #   plot_transect_timeseries("usa_NC_0033_0045", ts_data)
    return ts_data


if __name__ == "__main__":
    main()
