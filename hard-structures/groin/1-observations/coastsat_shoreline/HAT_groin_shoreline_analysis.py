"""
Buxton groin field: shoreline change rates around the groin, from three shoreline sources on CoastSat's transects.

    python hard-structures/groin/1-observations/coastsat_shoreline/HAT_groin_shoreline_analysis.py

Inputs: aerial wet-dry lines, NC Coastal Management shorelines and CoastSat,
all measured on CoastSat's own transects, plus the groin and domain layers in
gis_data/. Writes CSVs, profile figures and two GIFs to OUTPUT_DIR. Needs
geopandas, shapely, matplotlib and Pillow. Details: README.md beside this script.

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-10-01
"""

from pathlib import Path as _AnchorPath
import sys as _sys
from pathlib import Path as _RP
_sys.path.insert(0, str(next(_q for _q in _RP(__file__).resolve().parents
                             if (_q / "pyproject.toml").exists()) / "scripts"))
from site_layer import hat_observed_rates as _obs  # noqa: E402
import os
import glob
import json
import warnings
from pathlib import Path
from datetime import datetime

import numpy as np
import pandas as pd
import geopandas as gpd
from shapely.geometry import LineString, Point, MultiPoint, MultiLineString
from shapely.strtree import STRtree

import matplotlib.pyplot as plt
import matplotlib.collections
import matplotlib.dates as mdates
from matplotlib.ticker import FuncFormatter
from matplotlib.colors import Normalize, to_rgb
from matplotlib.cm import ScalarMappable
from matplotlib.gridspec import GridSpec
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from matplotlib.transforms import blended_transform_factory

warnings.filterwarnings("ignore")


# --- CONFIG ------------------------------------------------------------------
_PATH_REPO = next(_p for _p in _AnchorPath(__file__).resolve().parents
                  if (_p / "pyproject.toml").exists())

# Paths

STUDY_AREA_FILTER_PATH = str(_PATH_REPO / "hard-structures" / "groin" / "1-observations" / "gis_data" / "cascade_area.geojson")
STUDY_AREA_BUFFER_M    = 500

WET_DRY_PATH     = str(_PATH_REPO / "hard-structures" / "groin" / "1-observations" / "gis_data" / "wet_dry_groin.geojson")
WET_DRY_DATE_COL = "date"

# NC Coastal Management shorelines, read from the data tree (README)
NC_STATE_PATH     = str(_obs.SHORELINE_INVENTORY / "nc_shorelines.geojson")
NC_STATE_DATE_COL = "DATE_"

# CoastSat transects and time series, from the data tree
COASTSAT_TRANSECT_GEOM   = str(_obs.TRANSECT_LAYER)
COASTSAT_TRANSECT_ID_COL = "id"
COASTSAT_ROOT_DIR        = str(_obs.COASTSAT_TIMESERIES)

# CASCADE domain reference

# The 90 real domain boxes, D1-D90; replaced a fixed-spacing formula (README)
DOMAINS_JSON_PATH = str(_PATH_REPO / "hard-structures" / "groin" / "1-observations" / "gis_data" / "HAT_domains.json")

OUTPUT_DIR = str(_PATH_REPO / "hard-structures" / "groin" / "1-observations" / "coastsat_shoreline" / "shoreline_output_grid100m")

# CRS for all spatial operations (UTM 18N covers NC Outer Banks)
PROJECTED_CRS = "EPSG:32618"

# Plot x-axis limit

# Profile plots run from domain 1 (downdrift) to PLOT_UPDRIFT_MAX_KM updrift
PLOT_UPDRIFT_MAX_KM = 30
PLOT_X_MAX_KM = 60   # generic safety cap; kept as the GIF/other fallback default

# Groin location

# Groin features: the northernmost is the 0 km origin, the full extent is shaded
GROIN_GEOJSON_PATH = str(_PATH_REPO / "hard-structures" / "groin" / "1-observations" / "gis_data" / "groins_hatteras.geojson")

# Used only if the groin GeoJSON is missing; None to fail instead
GROIN_NORTHING_FALLBACK = 3_901_580.0    # from groins_hatteras.geojson center

# Net transport is southward, so updrift is north: distance is + north, - south
UPDRIFT_DIRECTION = "north"    # "north" or "south"

# Groin installation

# Breakpoint reference, and the cutoff for the pre-install baseline
GROIN_INSTALLATION_YEAR = 1970

# Pre-installation baseline

# Every observation before installation, pooled across sources (README)
PRE_INSTALLATION_YEAR_CUTOFF = GROIN_INSTALLATION_YEAR

# Nourishment zones (name, year, y_min, y_max): flagged, not dropped (README)
NOURISHMENT_EXCLUSIONS = [
    ("Buxton 1966",   1966, 3_899_000, 3_907_000),
    ("Buxton 1966 → 1967 shoreline", 1967, 3_899_000, 3_907_000),  # 1966 nourishment biases 1967
]

# Analysis eras

# Three eras from the groin's maintenance history, CSE 2013 (README)
ERAS = [
    ("Pre-install",      1849, PRE_INSTALLATION_YEAR_CUTOFF - 1),
    ("Functional groin",  1970, 1995),
    ("Deteriorated",      1996, 2024),
]

# Events drawn as annotated lines on the decadal plots, not era boundaries
NOTABLE_EVENTS = [
    (1970, "Groins built"),
    (1994, "Gordon damage"),
    (1995, "Last repair (S. groin)"),
    (2003, "Isabel"),
    (2022, "Jug Handle Bridge → mgmt ends"),
]

# Decadal windows

# Fixed, non-overlapping bins; replaced a 1-yr-step sliding window (README)
DECADE_START_YEAR      = 1960
DECADE_LENGTH_YEARS    = 10
DECADE_MIN_OBSERVATIONS = 5   # minimum observations in a decade to fit LRR

# Signal-extent detection

# Contiguous run outward from the groin where |anomaly| exceeds the threshold (README)
SIGNAL_ANOMALY_THRESHOLD_M_YR = 1.0     # m/yr above baseline required
SIGNAL_MAX_SEARCH_DISTANCE_M  = 15_000  # cap on search distance from groin
SIGNAL_EXTENT_BIN_WIDTH_M     = 500     # bin width for contiguity check (= CASCADE domain length)
SIGNAL_EXTENT_MAX_GAP_BINS    = 1       # below-threshold bins allowed in a row (absorbs single-bin noise)

# Distance bands (m) for the signal-through-time series, applied on both sides
DISTANCE_BAND_EDGES_M = [0, 1000, 3000, 6000, 10000, 15000]

# Piecewise breakpoint search

# Search range for the one breakpoint year; the fit itself uses every observation
BREAKPOINT_SEARCH_START = 1970
BREAKPOINT_SEARCH_END   = 2020
BREAKPOINT_MIN_POINTS_PER_SEGMENT = 5

# Transect chainage extraction

# Transect length used for intersections, longer than CoastSat's ~500 m
TRANSECT_INTERSECTION_LENGTH_M = 800

# Chainages beyond this are treated as spurious intersections
CHAINAGE_MAX_ABS_M = 700

# Minimum observations per transect to include it in analysis
MIN_OBSERVATIONS_PER_TRANSECT = 8

# Shoreline evolution GIF (groin area)

# One frame per observation date, or per calendar year pooled (README)
GIF_FRAME_MODE            = "date"    # "date" or "year"
# A frame needs both this fraction and this count of valid transects (README)
GIF_MIN_TRANSECT_FRACTION = 0.75   # min fraction of in-window transects valid to draw a frame
GIF_MIN_TRANSECT_ABS      = 3      # AND at least this many actual transects
# Before this year the fraction bar is waived: every sparse historical frame is kept
GIF_PRE_COASTSAT_CUTOFF_YEAR = 1984
# Second GIF, zoomed to domain 1 through this domain
GIF_ZOOM_WINDOW_DOMAIN_MAX = 25
GIF_FRAME_DURATION_S      = 0.15   # seconds per frame (many frames now -- keep this brisk)
GIF_TRAIL_FRAMES          = 5      # preceding frames shown as a fading trail
GIF_DPI                   = 130
# Last pre-install shoreline, drawn as a fixed reference on every later frame
GIF_REFERENCE_YEAR        = 1967
# Water seaward of the shoreline, island landward
GIF_WATER_COLOR = "#8FC1E3"
GIF_LAND_COLOR  = "#E3D5A8"
# Requires Pillow: pip install pillow

# GIF transect source

# coastsat, hybrid or grid100m; grid100m reintroduces a known artefact (README)
GIF_TRANSECT_SOURCE = "grid100m"   # "coastsat", "hybrid", or "grid100m"
TRANSECTS_100M_PATH = str(_PATH_REPO / "hard-structures" / "groin" / "1-observations" / "gis_data" / "transects_100m.geojson")

# CASCADE domain assignment

# ~500 m domains, D1 in the south to D90 in the north; boundaries from DOMAINS_JSON_PATH
NUM_REAL_DOMAINS   = 90
DOMAIN_LENGTH_M    = 500

# Plot styling

SOURCE_COLORS = {
    "wet_dry":  "#9C4D51",   # muted brick red
    "nc_state": "#3E5F7E",   # steel blue
    "coastsat": "#C27B49",   # warm terracotta
}
SOURCE_LABELS = {
    "wet_dry":  "Aerial wet-dry lines",
    "nc_state": "NC Coastal Mgmt",
    "coastsat": "CoastSat (satellite)",
}

ERA_COLORS = {
    "Pre-install":       "#5A5A5A",
    "Functional groin":  "#2E7D32",
    "Deteriorated":      "#C62828",
}

# Display-only sub-period: Deteriorated before the 2022 nourishments, not an official era
PRE_NOURISHMENT_PERIOD = ("1996-2021 (pre-2022 nourishment)", 1996, 2021)
PRE_NOURISHMENT_COLOR = "#F57C00"   # orange -- distinguishable from Deteriorated's red

# Decade-increment LRR plot

# LRR profile per fixed window since installation, one plot per increment
DECADE_PLOT_START_YEAR      = GROIN_INSTALLATION_YEAR   # 1970
DECADE_PLOT_INCREMENTS_YEARS = [5, 10]   # generates one plot per increment in this list
DECADE_PLOT_END_YEAR        = 2024
DECADE_PLOT_COLORMAP        = "Greens"   # monochromatic ramp -- color = chronological order

# Zone-panel profile plot

# One panel per era with shaded zones, adapted from HAT_groin_zone_investigation.py
ZONE_PANEL_DOWNDRIFT_DOMAINS = (1, 4)     # downdrift analysis zone
ZONE_PANEL_UPDRIFT_DOMAINS   = (7, 20)    # updrift analysis zone
ZONE_PANEL_X_MAX_DOMAIN      = 22         # a bit past the updrift zone, for padding
ZONE_PANEL_COLOR_DOWNDRIFT   = "#E67E22"   # orange
ZONE_PANEL_COLOR_UPDRIFT     = "#1565C0"   # blue
# Extra rows: the first N years of the Functional era; [] to skip
ZONE_PANEL_FUNCTIONAL_EARLY_WINDOWS_YEARS = [5, 10]

# Colour for difference lines, distinct from every era colour
DELTA_COLOR = "#6A3D9A"   # purple

# LOWESS fraction for the smoothed overlays, a visual aid only
SMOOTHED_OVERLAY_LOWESS_FRAC = 0.04

# Lower LOWESS fraction for the era profile, to keep the sharp jumps at the groin
ERA_PROFILE_LOWESS_FRAC = 0.02

# Text sizing (all plots)

FONT_TITLE      = 15
FONT_AXIS_LABEL = 13
FONT_TICK       = 11
FONT_LEGEND     = 10.5
FONT_ANNOTATION = 10
FONT_LEGEND_GIF = 8.5   # the GIF's legend is denser (data + Ocean/Island + all annotations)
FONT_LEGEND_GIF_ZOOMED = 7   # zoomed GIF is narrower, so the same legend needs to be smaller still

# Geographic annotations, positioned by CASCADE domain and converted to km when drawn

ANN_TOWN_SPANS = {
    "Buxton":      (7,  8),
    "Avon":        (21, 31),
    # Tri-Village removed: it sat on the edge of the 30 km window
}
ANN_VILLAGE_LINES = {}   # Salvo, Waves, Rodanthe all removed -- see notes above

# Annotations past this domain are dropped; spans crossing it are clipped
ANN_MAX_DOMAIN = 70

ANN_PIER_LABEL_Y  = 0.76   # default rotated label y for any pier (0=bottom, 1=top axes fraction)
ANN_GROIN_LABEL_Y = 0.68   # rotated label y for groin lines
ANN_PIERS = {
    "Avon Pier":     (26, ANN_PIER_LABEL_Y),   # (domain, label_y) - adjust per pier
    # Rodanthe Pier removed -- same edge-of-window collision as above.
}
# The groin is not listed: every plot already marks it
ANN_WIMBLE_SHOALS = (60, 74)
ANN_AVON_SHOALS   = (24, 39)   # Avon Shoals influence zone

# Accretion/erosion label heights: None for the midpoint, or a 0-1 axes fraction
LABEL_ACCRETION_Y = None
LABEL_EROSION_Y   = None

ANN_C_TOWN_SPAN    = "#90AFC5"
ANN_C_WIMBLE       = "#E0A800"   # amber - both shoal zones share this color
ANN_C_AVON_SHOALS  = "#E0A800"   # same amber as Wimble Shoals (same feature type)
ANN_C_VILLAGE_LINE = "0.40"
ANN_C_PIER         = "#1565C0"
ANN_C_GROIN        = "#B71C1C"
# Darker outlines for town and shoal spans, which can vanish against the GIF background
ANN_C_TOWN_SPAN_EDGE = "#2F5A73"
ANN_C_SHOAL_EDGE     = "#8A6800"

ANN_MODEL_COLOR = "#FF8C00"   # warm orange - modeled shoreline change rate
# -----------------------------------------------------------------------------


# Utilities

# Stop with a clear message if a configured file is missing (GDAL's own error is opaque)
def _require_file(path: str, description: str):
    if not os.path.isfile(path):
        raise FileNotFoundError(
            f"\n\n  Could not find {description} at:\n    {path}\n"
            f"  This path is set in the CONFIG section near the top of "
            f"this script. Either save the file to that exact location, "
            f"or update the config variable to point at wherever you "
            f"actually saved it.\n")


_CRS_LOG = []   # populated by fix_crs(); printed via print_crs_summary()


# Reproject to PROJECTED_CRS, fixing a missing CRS, and log it for the summary
def fix_crs(gdf: gpd.GeoDataFrame, source_name: str) -> gpd.GeoDataFrame:
    native_crs = gdf.crs
    if native_crs is None:
        print(f"  ! {source_name}: no CRS defined, assuming EPSG:4326")
        gdf = gdf.set_crs("EPSG:4326")
        native_crs = gdf.crs
    try:
        reprojected = gdf.to_crs(PROJECTED_CRS)
        bounds = reprojected.total_bounds if len(reprojected) else None
        _CRS_LOG.append({
            "source": source_name, "native_crs": str(native_crs),
            "target_crs": PROJECTED_CRS, "ok": True, "bounds": bounds,
        })
        return reprojected
    except Exception as e:
        print(f"  ! {source_name}: CRS reprojection failed ({e}), "
              f"falling back to EPSG:32618")
        fallback = gdf.set_crs("EPSG:32618", allow_override=True)
        _CRS_LOG.append({
            "source": source_name, "native_crs": str(native_crs),
            "target_crs": "EPSG:32618 (fallback -- reprojection failed)",
            "ok": False, "bounds": fallback.total_bounds if len(fallback) else None,
        })
        return fallback


# Print every loaded file's native CRS and reprojected bounds
def print_crs_summary():
    print(f"\n{'='*72}\nSpatial reference system check")
    print(f"  Working CRS for all analysis: {PROJECTED_CRS}")
    if not _CRS_LOG:
        print("  ! no files logged yet")
        return
    for entry in _CRS_LOG:
        status = "OK" if entry["ok"] else "FAILED -- used fallback, verify manually"
        b = entry["bounds"]
        bounds_str = (f"E {b[0]:,.0f}-{b[2]:,.0f}  N {b[1]:,.0f}-{b[3]:,.0f}"
                      if b is not None else "n/a")
        print(f"  {entry['source']:<30} native={entry['native_crs']:<18} "
              f"reprojected bounds: {bounds_str}  [{status}]")
    print(f"  (All bounds should fall in a similar range if every layer "
          f"is genuinely covering the same stretch of coast -- a layer "
          f"with wildly different Easting/Northing here means its "
          f"reprojection is suspect, check that file's native CRS.)")


# Parse a date column stored as text, Unix milliseconds or a bare year
def parse_date_column(series: pd.Series, column_name: str) -> pd.Series:
    if pd.api.types.is_numeric_dtype(series):
        max_val = series.abs().max()
        if max_val > 1e9:
            # Unix milliseconds since epoch
            print(f"  '{column_name}': Unix milliseconds → parsing with unit='ms'")
            return pd.to_datetime(series, unit="ms", errors="coerce")
        elif max_val < 3000:
            # Bare year integer
            print(f"  '{column_name}': bare year integer → constructing YYYY-01-01")
            return pd.to_datetime(series.astype("Int64").astype(str) + "-01-01",
                                  errors="coerce")
    # String dates (assume ISO or common US formats)
    return pd.to_datetime(series, errors="coerce")


# Datetime to decimal year (1985-07-02 -> ~1985.50)
def decimal_year(dt: pd.Series) -> pd.Series:
    dt = pd.to_datetime(dt)
    years = dt.dt.year
    year_starts = pd.to_datetime(years.astype(str) + "-01-01")
    year_ends = pd.to_datetime((years + 1).astype(str) + "-01-01")
    frac = (dt - year_starts) / (year_ends - year_starts)
    return years + frac


# The 90 CASCADE domain polygons, with their midpoint northings
def load_domain_reference() -> pd.DataFrame:
    print(f"\nLoading CASCADE domain reference: "
          f"{os.path.basename(DOMAINS_JSON_PATH)}")
    _require_file(DOMAINS_JSON_PATH,
                  "the CASCADE domain reference (DOMAINS_JSON_PATH)")
    gdf = gpd.read_file(DOMAINS_JSON_PATH)
    print(f"  {len(gdf)} domain polygons, native CRS: {gdf.crs}")
    gdf = fix_crs(gdf, "CASCADE domains")

    id_col = "domain_id" if "domain_id" in gdf.columns else gdf.columns[0]
    records = []
    for _, row in gdf.iterrows():
        minx, miny, maxx, maxy = row.geometry.bounds
        records.append({"domain_id": int(row[id_col]), "y_min": miny, "y_max": maxy})
    domain_table = pd.DataFrame(records).sort_values("domain_id").reset_index(drop=True)
    domain_table["y_mid"] = (domain_table["y_min"] + domain_table["y_max"]) / 2.0

    print(f"  Domain range: D{domain_table['domain_id'].min()} to "
          f"D{domain_table['domain_id'].max()}")
    print(f"  Northing range: {domain_table['y_min'].min():,.0f} to "
          f"{domain_table['y_max'].max():,.0f}")
    return domain_table


# Northing to CASCADE domain, by nearest domain midpoint
def assign_domain_from_northing(y, domain_table: pd.DataFrame):
    y_arr = np.atleast_1d(np.asarray(y, dtype=float))
    mids = domain_table["y_mid"].values
    ids = domain_table["domain_id"].values
    order = np.argsort(mids)
    mids_sorted = mids[order]
    ids_sorted = ids[order]
    idx = np.searchsorted(mids_sorted, y_arr)
    idx = np.clip(idx, 1, len(mids_sorted) - 1)
    left = idx - 1
    choose_left = np.abs(y_arr - mids_sorted[left]) <= np.abs(y_arr - mids_sorted[idx])
    nearest = np.where(choose_left, left, idx)
    result = ids_sorted[nearest].astype(int)
    return int(result[0]) if np.isscalar(y) or np.ndim(y) == 0 else result


# Mean distance from the groin (km) of a domain's transects, or None
def domain_dist_km(transects: gpd.GeoDataFrame, domain_number: int) -> float:
    sub = transects[transects["domain"] == domain_number]
    if len(sub) == 0:
        print(f"  ! domain_dist_km: no transects assigned to domain "
              f"D{domain_number} -- falling back to auto-detected x-limit")
        return None
    pos_km = float(sub["dist_from_groin_m"].mean()) / 1000.0
    print(f"  domain_dist_km: D{domain_number} -> {pos_km:+.2f} km "
          f"({len(sub)} transect(s): "
          f"{sorted(sub['transect_id'].astype(str).tolist())[:5]}"
          f"{'...' if len(sub) > 5 else ''})")
    return pos_km


# Domain number <-> distance from the groin (km), by interpolation
def build_domain_to_dist_km(transects: gpd.GeoDataFrame):
    agg = (transects.groupby("domain")["dist_from_groin_m"]
                     .mean().sort_index())
    domains = agg.index.values.astype(float)
    dist_km = agg.values / 1000.0
    order = np.argsort(dist_km)
    dist_sorted = dist_km[order]
    dom_sorted = domains[order]

    def dist_to_domain(x):
        return np.interp(x, dist_sorted, dom_sorted)

    def domain_to_dist(d):
        return np.interp(d, domains, dist_km)

    return domain_to_dist, dist_to_domain


# Domain number on the bottom axis, distance from the groin (km) on the top
def set_domain_primary_axis(ax, transects: gpd.GeoDataFrame):
    agg = (transects.groupby("domain")["dist_from_groin_m"]
                     .mean().sort_index())
    if len(agg) < 2:
        return None
    dist_km = agg.values / 1000.0
    domain_to_dist, dist_to_domain = build_domain_to_dist_km(transects)

    x_lo, x_hi = ax.get_xlim()
    d_lo, d_hi = float(dist_to_domain(x_lo)), float(dist_to_domain(x_hi))
    d_lo, d_hi = min(d_lo, d_hi), max(d_lo, d_hi)
    step = 5
    domain_ticks = np.arange(np.floor(d_lo / step) * step, d_hi + step, step)
    domain_ticks = domain_ticks[domain_ticks >= 1]
    tick_positions_km = domain_to_dist(domain_ticks)
    # Skip ticks within a small margin of the plot edges, where labels collide
    edge_margin_km = 0.6
    visible = ((tick_positions_km >= x_lo + edge_margin_km) &
                (tick_positions_km <= x_hi - edge_margin_km))

    # Always show D65 when it is truly in range; domain_to_dist clips, so check first
    domain_min, domain_max = agg.index.min(), agg.index.max()
    if 65 not in domain_ticks[visible] and domain_min <= 65 <= domain_max:
        km_65 = float(domain_to_dist(65))
        if x_lo <= km_65 <= x_hi:
            domain_ticks = np.append(domain_ticks, 65.0)
            tick_positions_km = domain_to_dist(domain_ticks)
            visible = ((tick_positions_km >= x_lo + edge_margin_km) &
                        (tick_positions_km <= x_hi - edge_margin_km))
            visible = visible | (domain_ticks == 65.0)

    ax.set_xticks(tick_positions_km[visible])
    ax.set_xticklabels([f"D{int(t)}" for t in domain_ticks[visible]], fontsize=FONT_TICK)
    ax.set_xlabel("CASCADE domain", fontsize=FONT_AXIS_LABEL)
    ax.tick_params(axis="y", labelsize=FONT_TICK)

    # Thin line at every domain position, not just the labelled ticks
    all_positions_km = dist_km  # one position per domain (order doesn't matter here)
    all_visible = (all_positions_km >= x_lo) & (all_positions_km <= x_hi)
    for pos in all_positions_km[all_visible]:
        ax.axvline(pos, color="#bbbbbb", lw=0.4, alpha=0.6, zorder=0)

    secax = ax.secondary_xaxis("top", functions=(lambda x: x, lambda x: x))
    secax.set_xlabel("Distance from groin (km, alongshore)  [+ = updrift]",
                       fontsize=FONT_AXIS_LABEL)
    secax.tick_params(labelsize=FONT_TICK)
    return secax


# A lighter shade of a hex colour (0 = unchanged, 1 = white)
def lighten_color(hex_color: str, amount: float = 0.45):
    rgb = np.array(to_rgb(hex_color))
    return tuple(rgb + (np.array([1.0, 1.0, 1.0]) - rgb) * amount)


# Towns, piers and shoals drawn on a distance-from-groin plot
def add_geographic_annotations(ax, transects: gpd.GeoDataFrame,
                                gif_mode: bool = False):
    domain_to_dist, _ = build_domain_to_dist_km(transects)
    x_lo, x_hi = ax.get_xlim()
    trans = blended_transform_factory(ax.transData, ax.transAxes)

    def visible(km):
        return x_lo <= km <= x_hi

    span_alpha = 0.30 if gif_mode else 0.14
    shoal_alpha = 0.30 if gif_mode else 0.10
    base_zorder = 2 if gif_mode else 0   # above the GIF's water/island fill (zorder=1)
    bbox = dict(boxstyle="round,pad=0.2", fc="white", ec="none", alpha=0.82)

    # ─── Shoal zones (hatched fill + label) -- drawn first/underneath ──
    for name, (d_lo, d_hi), color in [
        ("Wimble Shoals", ANN_WIMBLE_SHOALS, ANN_C_WIMBLE),
        ("Avon Shoals", ANN_AVON_SHOALS, ANN_C_AVON_SHOALS),
    ]:
        if d_lo > ANN_MAX_DOMAIN:
            continue
        d_hi = min(d_hi, ANN_MAX_DOMAIN)
        km_lo, km_hi = float(domain_to_dist(d_lo)), float(domain_to_dist(d_hi))
        km_lo, km_hi = min(km_lo, km_hi), max(km_lo, km_hi)
        if km_hi < x_lo or km_lo > x_hi:
            continue
        clipped_lo, clipped_hi = max(km_lo, x_lo), min(km_hi, x_hi)
        ax.axvspan(clipped_lo, clipped_hi, color=color, alpha=shoal_alpha,
                    zorder=base_zorder, hatch="///", edgecolor=color, linewidth=0)
        mid_km = (clipped_lo + clipped_hi) / 2.0
        ax.text(mid_km, 0.04, name, transform=trans,
                 ha="center", va="bottom", fontsize=FONT_ANNOTATION,
                 color="#7A5800", style="italic", bbox=bbox, zorder=base_zorder + 10)

    # Town spans: shaded band and label, nudged off-centre where a pier would collide
    for name, (d_lo, d_hi) in ANN_TOWN_SPANS.items():
        if d_lo > ANN_MAX_DOMAIN:
            continue
        d_hi = min(d_hi, ANN_MAX_DOMAIN)
        km_lo, km_hi = float(domain_to_dist(d_lo)), float(domain_to_dist(d_hi))
        km_lo, km_hi = min(km_lo, km_hi), max(km_lo, km_hi)
        if km_hi < x_lo or km_lo > x_hi:
            continue
        clipped_lo, clipped_hi = max(km_lo, x_lo), min(km_hi, x_hi)
        ax.axvspan(clipped_lo, clipped_hi, color=ANN_C_TOWN_SPAN,
                    alpha=span_alpha, zorder=base_zorder)
        mid_domain = (d_lo + d_hi) / 2.0
        pier_conflict = any(abs(pd_ - mid_domain) < 0.15 * (d_hi - d_lo)
                             for pd_, _ in ANN_PIERS.values())
        frac = 0.25 if pier_conflict else 0.5
        label_km = clipped_lo + frac * (clipped_hi - clipped_lo)
        ax.text(label_km, 0.90, name, transform=trans,
                 ha="center", va="top", fontsize=FONT_ANNOTATION, color="0.25",
                 fontweight="bold", bbox=bbox, zorder=base_zorder + 10)

    # Village center lines
    for name, d in ANN_VILLAGE_LINES.items():
        if d > ANN_MAX_DOMAIN:
            continue
        km = float(domain_to_dist(d))
        if not visible(km):
            continue
        ax.axvline(km, color=ANN_C_VILLAGE_LINE, lw=0.9, linestyle="--",
                    alpha=0.65, zorder=base_zorder + 1)
        ax.text(km, 0.84, name, transform=trans,
                 ha="center", va="top", fontsize=FONT_ANNOTATION, color="0.30",
                 bbox=bbox, zorder=base_zorder + 10)

    # Piers
    for name, (d, label_y_frac) in ANN_PIERS.items():
        if d > ANN_MAX_DOMAIN:
            continue
        km = float(domain_to_dist(d))
        if not visible(km):
            continue
        ax.axvline(km, color=ANN_C_PIER, lw=1.0, linestyle="-.",
                    alpha=0.80, zorder=base_zorder + 2)
        ax.text(km, label_y_frac, name, transform=trans, rotation=90,
                 ha="center", va="top", fontsize=FONT_ANNOTATION,
                 color=ANN_C_PIER, bbox=bbox, zorder=base_zorder + 10)


# One legend entry per annotation layer type
def annotation_legend_handles():
    handles = [Patch(fc=ANN_C_TOWN_SPAN, alpha=0.30, label="Community")]
    if ANN_WIMBLE_SHOALS[0] <= ANN_MAX_DOMAIN or ANN_AVON_SHOALS[0] <= ANN_MAX_DOMAIN:
        handles.append(Patch(fc=ANN_C_WIMBLE, alpha=0.25, hatch="///",
                               edgecolor=ANN_C_WIMBLE, linewidth=0,
                               label="Shoals position (Avon / Wimble)"))
    if ANN_VILLAGE_LINES:
        handles.append(Line2D([0], [0], color=ANN_C_VILLAGE_LINE, lw=0.9,
                                ls="--", label="Village center"))
    if ANN_PIERS:
        handles.append(Line2D([0], [0], color=ANN_C_PIER, lw=1.0, ls="-.",
                                label="Pier"))
    return handles


# (x, y) of every real data artist on an axis, reference lines skipped
def _iter_axes_real_xy(ax):
    for line in ax.get_lines():
        xd, yd = np.asarray(line.get_xdata()), np.asarray(line.get_ydata())
        if len(xd) == 0:
            continue
        if len(xd) <= 2 and (np.ptp(xd) == 0 or np.ptp(yd) == 0):
            continue
        yield xd, yd
    for coll in ax.collections:   # scatter points + fill_between IQR bands
        if isinstance(coll, matplotlib.collections.PathCollection):
            # scatter: real data lives in offsets, not the marker path
            offsets = np.asarray(coll.get_offsets())
            if offsets.ndim == 2 and offsets.shape[0] > 0:
                yield offsets[:, 0], offsets[:, 1]
        else:
            # fill_between / PolyCollection: the real data is the polygon vertices
            verts = np.concatenate([p.vertices for p in coll.get_paths()], axis=0) \
                    if coll.get_paths() else np.empty((0, 2))
            if len(verts) > 0:
                yield verts[:, 0], verts[:, 1]


# Clip x to the window and fit y to the data visible in it
def apply_distance_xlim_and_autoscale_y(ax, x_max_km: float = None,
                                         x_min_km: float = None,
                                         y_margin_frac: float = 0.10,
                                         x_data_margin_km: float = 1.0,
                                         y_percentile: tuple = (1, 99)):
    if x_max_km is None:
        x_max_km = PLOT_X_MAX_KM

    if x_min_km is None:
        data_x_min = np.inf
        for xd, _ in _iter_axes_real_xy(ax):
            finite = np.isfinite(xd)
            if finite.any():
                data_x_min = min(data_x_min, float(xd[finite].min()))
        x_min_km = (max(-x_max_km, data_x_min - x_data_margin_km)
                    if np.isfinite(data_x_min) else -x_max_km)

    ax.set_xlim(x_min_km, x_max_km)

    all_y_vis = []
    for xd, yd in _iter_axes_real_xy(ax):
        mask = (xd >= x_min_km) & (xd <= x_max_km)
        yd_vis = yd[mask]
        finite = np.isfinite(yd_vis)
        if finite.any():
            all_y_vis.append(yd_vis[finite])

    if all_y_vis:
        all_y_vis = np.concatenate(all_y_vis)
        y_min, y_max = np.percentile(all_y_vis, y_percentile)
        if y_max > y_min:
            pad = (y_max - y_min) * y_margin_frac
            ax.set_ylim(y_min - pad, y_max + pad)
    # else leave matplotlib's own autoscale in place


# Groin geometry

# Groin features: north-end origin, field extent, per-groin reference points
def load_groin_geometry(domain_table: pd.DataFrame = None) -> dict:
    print(f"\nLoading groin geometry: {os.path.basename(GROIN_GEOJSON_PATH)}")
    if not os.path.exists(GROIN_GEOJSON_PATH):
        if GROIN_NORTHING_FALLBACK is None:
            raise FileNotFoundError(
                f"Groin file not found: {GROIN_GEOJSON_PATH}")
        print(f"  ! Not found — using GROIN_NORTHING_FALLBACK = "
              f"{GROIN_NORTHING_FALLBACK:,.0f}")
        return dict(reference_x=None, reference_y=GROIN_NORTHING_FALLBACK,
                    y_min=GROIN_NORTHING_FALLBACK,
                    y_max=GROIN_NORTHING_FALLBACK,
                    southernmost_x=None, southernmost_y=GROIN_NORTHING_FALLBACK,
                    middle_x=None, middle_y=GROIN_NORTHING_FALLBACK,
                    geometry=None,
                    n_features=0)

    gdf = gpd.read_file(GROIN_GEOJSON_PATH)
    print(f"  {len(gdf)} feature(s), native CRS: {gdf.crs}")
    gdf = fix_crs(gdf, "groin")
    geom_types = sorted(set(gdf.geometry.geom_type))
    print(f"  Geometry types: {geom_types}")
    if not any("Line" in gt for gt in geom_types):
        print(f"  ! WARNING: expected LineString feature(s) for a groin "
              f"structure, got {geom_types} instead. This usually means "
              f"GROIN_GEOJSON_PATH points at the wrong file -- double "
              f"check it's the groin shapefile, not e.g. the domain "
              f"reference or study area filter.")
    if len(gdf) > 20:
        print(f"  ! WARNING: {len(gdf)} features is a lot for a single "
              f"groin field ({GROIN_GEOJSON_PATH}) -- verify this is "
              f"really the groin file and not something else entirely.")

    merged = gdf.unary_union
    minx, miny, maxx, maxy = merged.bounds

    # Per-feature centroids: the northernmost is the 0 km origin; the full span is kept for shading
    centroids = gdf.geometry.centroid
    northernmost_idx = centroids.y.idxmax()
    reference_x = float(centroids.x.loc[northernmost_idx])
    reference_y = float(centroids.y.loc[northernmost_idx])

    print(f"  Groin field extent (northing): {miny:,.0f} → {maxy:,.0f} "
          f"({(maxy - miny):.0f} m span across {len(gdf)} feature(s))")
    print(f"  Northernmost groin feature centroid (distance-from-groin "
          f"origin): x={reference_x:,.0f}, y={reference_y:,.0f}")

    if domain_table is not None:
        ref_domain = assign_domain_from_northing(reference_y, domain_table)
        span_domains = sorted(set(
            int(d) for d in assign_domain_from_northing(
                np.array([miny, maxy]), domain_table)))
        print(f"  → Northernmost groin feature is in domain D{ref_domain}; "
              f"the whole groin field spans domain(s) "
              f"D{span_domains[0]}-D{span_domains[-1]}"
              if len(span_domains) > 1 else
              f"  → Northernmost groin feature is in domain D{ref_domain}; "
              f"the whole groin field is within domain D{span_domains[0]}")
        print(f"  ^ if this doesn't match where you expect the groin to "
              f"be, GROIN_GEOJSON_PATH most likely points at the wrong "
              f"file.")

    # Southernmost and middle groins too, for lighter reference lines
    order = np.argsort(centroids.y.values)
    sorted_x = centroids.x.values[order]
    sorted_y = centroids.y.values[order]
    southernmost_x, southernmost_y = float(sorted_x[0]), float(sorted_y[0])
    mid_idx = len(sorted_y) // 2
    middle_x, middle_y = float(sorted_x[mid_idx]), float(sorted_y[mid_idx])

    return dict(reference_x=reference_x, reference_y=reference_y,
                y_min=miny, y_max=maxy,
                southernmost_x=southernmost_x, southernmost_y=southernmost_y,
                middle_x=middle_x, middle_y=middle_y,
                geometry=merged, n_features=len(gdf))


# Signed alongshore distance from the groin (m), + updrift
def assign_distance_from_groin(transects: gpd.GeoDataFrame,
                                groin: dict) -> gpd.GeoDataFrame:
    sign = 1 if UPDRIFT_DIRECTION.lower().startswith("n") else -1
    tx = transects.copy()

    # The groin's alongshore position: the transect whose shore point is nearest the northernmost groin
    if "shore_x" in tx.columns:
        ref_x, ref_y = tx["shore_x"], tx["shore_y"]
    else:
        ref_x, ref_y = tx["origin_x"], tx["origin_y"]

    gx = groin.get("reference_x")
    gy = groin["reference_y"]
    if gx is not None:
        d2 = (ref_x - gx) ** 2 + (ref_y - gy) ** 2
    else:
        # No geometry available -- fall back to nearest-by-northing
        d2 = (ref_y - gy) ** 2
    groin_idx = d2.idxmin()
    groin_alongshore_m = float(tx.loc[groin_idx, "alongshore_m"])

    tx["dist_from_groin_m"] = sign * (tx["alongshore_m"] - groin_alongshore_m)
    print(f"  Groin alongshore position: {groin_alongshore_m/1000:.2f} km "
          f"(nearest transect id={tx.loc[groin_idx, 'transect_id']}, "
          f"{d2.loc[groin_idx]**0.5:.0f} m away)")
    return tx


# Study area filter

# Study-area polygon, buffered
def load_study_area_filter() -> gpd.GeoDataFrame:
    print(f"\nLoading study area filter: {os.path.basename(STUDY_AREA_FILTER_PATH)}")
    _require_file(STUDY_AREA_FILTER_PATH, "the study area filter (STUDY_AREA_FILTER_PATH)")
    gdf = gpd.read_file(STUDY_AREA_FILTER_PATH)
    gdf = fix_crs(gdf, "study area filter")

    geom_types = set(gdf.geometry.geom_type)
    if geom_types & {"Polygon", "MultiPolygon"}:
        # Merge all polygons, then buffer by STUDY_AREA_BUFFER_M
        merged = gdf.unary_union.buffer(STUDY_AREA_BUFFER_M)
    else:
        # Lines: buffer to create corridor
        merged = gdf.unary_union.buffer(STUDY_AREA_BUFFER_M)

    filter_gdf = gpd.GeoDataFrame(geometry=[merged], crs=PROJECTED_CRS)
    bounds = filter_gdf.total_bounds
    print(f"  Filter polygon area: {filter_gdf.geometry.area.sum() / 1e6:.2f} km²")
    print(f"  Extent: {(bounds[2]-bounds[0])/1000:.1f} × "
          f"{(bounds[3]-bounds[1])/1000:.1f} km")
    return filter_gdf


# Loader: wet-dry and NC state shorelines

# Shoreline lines from a GeoJSON: projected, clipped to the study area, dated
def load_shoreline_lines(path: str, date_col: str, source_key: str,
                         filter_gdf: gpd.GeoDataFrame) -> gpd.GeoDataFrame:
    print(f"\n[{SOURCE_LABELS[source_key]}] loading: {os.path.basename(path)}")
    if not os.path.exists(path):
        print(f"  ! File not found: {path}")
        return gpd.GeoDataFrame(columns=["geometry", "date", "year", "source"],
                                 geometry="geometry", crs=PROJECTED_CRS)

    gdf = gpd.read_file(path)
    print(f"  {len(gdf)} raw features, CRS: {gdf.crs}")
    gdf = fix_crs(gdf, source_key)

    # Spatial filter
    filter_geom = filter_gdf.geometry.iloc[0]
    keep_mask = gdf.geometry.intersects(filter_geom)
    gdf = gdf.loc[keep_mask].copy()
    print(f"  After spatial filter: {len(gdf)} features")

    if len(gdf) == 0:
        return gpd.GeoDataFrame(columns=["geometry", "date", "year", "source"],
                                 geometry="geometry", crs=PROJECTED_CRS)

    # Parse dates
    if date_col not in gdf.columns:
        available = ", ".join(gdf.columns[:15])
        raise ValueError(f"Date column '{date_col}' not found. "
                         f"Available: {available}")

    gdf["_parsed_date"] = parse_date_column(gdf[date_col], date_col)
    gdf = gdf.dropna(subset=["_parsed_date"])
    gdf["date"] = gdf["_parsed_date"]
    gdf["year"] = gdf["_parsed_date"].dt.year
    gdf["source"] = source_key

    print(f"  → {len(gdf)} valid features, "
          f"{gdf['year'].nunique()} unique years, "
          f"{gdf['date'].nunique()} unique dates")

    return gdf[["geometry", "date", "year", "source"]].reset_index(drop=True)


# Loader: CoastSat transect network (the analysis backbone)

# CoastSat's transects, the backbone for every source: geometry, order, domain
def load_coastsat_transects(filter_gdf: gpd.GeoDataFrame,
                             domain_table: pd.DataFrame) -> gpd.GeoDataFrame:
    print(f"\n[CoastSat transects] loading: "
          f"{os.path.basename(COASTSAT_TRANSECT_GEOM)}")
    _require_file(COASTSAT_TRANSECT_GEOM, "the CoastSat transect layer (COASTSAT_TRANSECT_GEOM)")
    gdf = gpd.read_file(COASTSAT_TRANSECT_GEOM)
    print(f"  {len(gdf)} global transects, native CRS: {gdf.crs}")
    gdf = fix_crs(gdf, "CoastSat transects")

    # Spatial filter using transect origin (start point)
    def get_origin(geom):
        coords = list(geom.coords)
        return Point(coords[0])

    gdf["origin"] = gdf.geometry.apply(get_origin)
    filter_geom = filter_gdf.geometry.iloc[0]
    keep_mask = gpd.GeoSeries(gdf["origin"], crs=PROJECTED_CRS).within(filter_geom)
    gdf = gdf.loc[keep_mask].copy()
    print(f"  In study area: {len(gdf)} transects")

    if len(gdf) == 0:
        raise RuntimeError("No CoastSat transects fell within study area — "
                           "check STUDY_AREA_FILTER_PATH and transect CRS.")

    # Compute origin and unit direction from geometry
    def unpack_geom(geom):
        coords = list(geom.coords)
        x0, y0 = coords[0]
        x1, y1 = coords[-1]
        dx, dy = x1 - x0, y1 - y0
        length = np.hypot(dx, dy)
        if length == 0:
            return x0, y0, 1.0, 0.0
        return x0, y0, dx / length, dy / length

    unpacked = gdf.geometry.apply(unpack_geom)
    gdf["origin_x"] = unpacked.apply(lambda t: t[0])
    gdf["origin_y"] = unpacked.apply(lambda t: t[1])
    gdf["dir_x"]    = unpacked.apply(lambda t: t[2])
    gdf["dir_y"]    = unpacked.apply(lambda t: t[3])
    gdf["transect_id"] = gdf[COASTSAT_TRANSECT_ID_COL].astype(str)
    # Shore-side reference point for groin matching; the origin is already near shore here
    gdf["shore_x"] = gdf["origin_x"]
    gdf["shore_y"] = gdf["origin_y"]

    # Alongshore position: trusted ID order if possible, else NN fallback
    gdf = _compute_alongshore_positions(gdf, id_col=COASTSAT_TRANSECT_ID_COL)

    gdf["domain"] = assign_domain_from_northing(gdf["origin_y"].values, domain_table)

    print(f"  → alongshore extent: 0 → "
          f"{gdf['alongshore_m'].max()/1000:.1f} km")
    print(f"  → domain range: D{gdf['domain'].min()} to D{gdf['domain'].max()}")

    return gdf[[
        "transect_id", "geometry",
        "origin_x", "origin_y", "dir_x", "dir_y",
        "shore_x", "shore_y", "alongshore_m", "domain",
    ]].reset_index(drop=True)


# The parallel 100 m transect grid, used only by the GIF's hybrid and grid100m modes
def load_100m_grid_transects(domain_table: pd.DataFrame) -> gpd.GeoDataFrame:
    print(f"\n[100m grid transects] loading: "
          f"{os.path.basename(TRANSECTS_100M_PATH)}")
    _require_file(TRANSECTS_100M_PATH,
                  "the 100 m transect grid (TRANSECTS_100M_PATH) -- "
                  "only needed for GIF_TRANSECT_SOURCE != 'coastsat'")
    gdf = gpd.read_file(TRANSECTS_100M_PATH)
    print(f"  {len(gdf)} transects, native CRS: {gdf.crs}")
    gdf = fix_crs(gdf, "100m grid transects")

    id_col = ("Transects_100m.LineID" if "Transects_100m.LineID" in gdf.columns
               else gdf.columns[0])
    gdf["transect_id"] = ("grid100m_" +
        pd.to_numeric(gdf[id_col], errors="coerce").astype("Int64").astype(str))

    def unpack(geom):
        coords = list(geom.coords)
        sx, sy = coords[0]    # seaward end (on the straight reference line)
        lx, ly = coords[-1]   # landward end
        # Direction points landward to seaward, as in CoastSat, so + chainage = accretion (README)
        dx, dy = sx - lx, sy - ly
        length = np.hypot(dx, dy)
        if length == 0:
            return sx, sy, lx, ly, 1.0, 0.0
        return sx, sy, lx, ly, dx / length, dy / length

    unpacked = gdf.geometry.apply(unpack)
    gdf["shore_x"]  = unpacked.apply(lambda t: t[0])
    gdf["shore_y"]  = unpacked.apply(lambda t: t[1])
    gdf["origin_x"] = unpacked.apply(lambda t: t[2])
    gdf["origin_y"] = unpacked.apply(lambda t: t[3])
    gdf["dir_x"]    = unpacked.apply(lambda t: t[4])
    gdf["dir_y"]    = unpacked.apply(lambda t: t[5])

    gdf = gdf.sort_values("shore_y").reset_index(drop=True)
    gdf["alongshore_m"] = gdf["shore_y"] - gdf["shore_y"].min()
    gdf["domain"] = assign_domain_from_northing(gdf["shore_y"].values, domain_table)

    print(f"  → alongshore extent: 0 → "
          f"{gdf['alongshore_m'].max() / 1000:.1f} km")
    return gdf[["transect_id", "geometry", "shore_x", "shore_y",
                "origin_x", "origin_y", "dir_x", "dir_y",
                "alongshore_m", "domain"]].reset_index(drop=True)


# Hybrid GIF mode: each CoastSat transect's position on the 100 m grid
def build_hybrid_position_map(transects: gpd.GeoDataFrame,
                               grid100m: gpd.GeoDataFrame,
                               groin: dict) -> pd.Series:
    gx, gy = groin.get("reference_x"), groin["reference_y"]
    if gx is not None:
        d2 = (grid100m["shore_x"] - gx) ** 2 + (grid100m["shore_y"] - gy) ** 2
    else:
        d2 = (grid100m["shore_y"] - gy) ** 2
    groin_alongshore_100m = float(grid100m.loc[d2.idxmin(), "alongshore_m"])
    sign = 1 if UPDRIFT_DIRECTION.lower().startswith("n") else -1

    grid_xy = grid100m[["shore_x", "shore_y"]].values
    grid_along = grid100m["alongshore_m"].values
    cs_xy = transects[["shore_x", "shore_y"]].values
    cs_ids = transects["transect_id"].astype(str).values

    results = {}
    for i, tid in enumerate(cs_ids):
        dists = np.hypot(grid_xy[:, 0] - cs_xy[i, 0], grid_xy[:, 1] - cs_xy[i, 1])
        nearest = int(np.argmin(dists))
        results[tid] = sign * (grid_along[nearest] - groin_alongshore_100m)
    return pd.Series(results, name="dist_from_groin_m_100mgrid")


# grid100m GIF mode: rebuild each point in x, y and re-measure it on the 100 m grid
def reconstruct_and_remeasure_on_grid100m(chainage_sub: pd.DataFrame,
                                           coastsat_transects: gpd.GeoDataFrame,
                                           grid100m: gpd.GeoDataFrame
                                           ) -> pd.DataFrame:
    tx = coastsat_transects.set_index("transect_id")[
        ["origin_x", "origin_y", "dir_x", "dir_y"]].copy()
    tx.index = tx.index.astype(str)

    sub = chainage_sub.copy()
    sub["transect_id"] = sub["transect_id"].astype(str)
    sub = sub[sub["transect_id"].isin(tx.index)].reset_index(drop=True)
    if len(sub) == 0:
        return sub

    row_info = tx.loc[sub["transect_id"]]
    ox, oy = row_info["origin_x"].values, row_info["origin_y"].values
    dx, dy = row_info["dir_x"].values, row_info["dir_y"].values
    ch = sub["chainage_m"].values
    px = ox + ch * dx
    py = oy + ch * dy

    grid_lines = grid100m.geometry.tolist()
    tree = STRtree(grid_lines)
    grid_origin = grid100m[["origin_x", "origin_y"]].values
    grid_dir = grid100m[["dir_x", "dir_y"]].values
    grid_ids = grid100m["transect_id"].values

    new_chainage = np.full(len(sub), np.nan)
    new_tid = np.full(len(sub), None, dtype=object)
    match_dist = np.full(len(sub), np.nan)
    for i in range(len(sub)):
        pt = Point(px[i], py[i])
        # True nearest-neighbour query with no distance cap, so far points are never dropped (README)
        best_idx = int(tree.nearest(pt))
        match_dist[i] = grid_lines[best_idx].distance(pt)
        gox, goy = grid_origin[best_idx]
        gdx, gdy = grid_dir[best_idx]
        new_chainage[i] = (px[i] - gox) * gdx + (py[i] - goy) * gdy
        new_tid[i] = grid_ids[best_idx]

    if len(match_dist) > 0:
        print(f"  grid100m match distance (m): median {np.median(match_dist):.0f}, "
              f"95th pct {np.percentile(match_dist, 95):.0f}, "
              f"max {np.max(match_dist):.0f} "
              f"-- large values mean the reconstructed point sat far "
              f"from any 100m-grid transect (most likely near sharp "
              f"coastline curvature), not that data is missing")

    out = sub.copy()
    out["transect_id"] = new_tid
    out["chainage_m"] = new_chainage
    return out.dropna(subset=["chainage_m", "transect_id"]).reset_index(drop=True)


# Cumulative alongshore distance per transect, in transect-ID order
def _compute_alongshore_positions(gdf: gpd.GeoDataFrame,
                                   id_col: str = None) -> gpd.GeoDataFrame:
    gdf = gdf.copy()
    order = None
    if id_col is not None and id_col in gdf.columns:
        numeric_id = pd.to_numeric(gdf[id_col], errors="coerce")
        if numeric_id.notna().all() and numeric_id.nunique() == len(gdf):
            order = np.argsort(numeric_id.values)
            print(f"  Alongshore order: sorted by '{id_col}' "
                  f"(trusted sequential ID order)")
        else:
            # Compound IDs '{site}-{index}': sort by (site, number) to recover alongshore order
            id_str = gdf[id_col].astype(str).reset_index(drop=True)
            suffix_num = pd.to_numeric(
                id_str.str.extract(r"(\d+)$")[0], errors="coerce")
            prefix = id_str.str.replace(r"[-_]?\d+$", "", regex=True)
            if suffix_num.notna().all():
                sort_df = pd.DataFrame({"prefix": prefix, "suffix_num": suffix_num})
                order = sort_df.sort_values(["prefix", "suffix_num"]).index.values
                print(f"  Alongshore order: sorted by '{id_col}' site prefix + "
                      f"numeric suffix (trusted compound ID order, "
                      f"{prefix.nunique()} site chunk(s): "
                      f"{sorted(prefix.unique().tolist())})")

    coords = np.column_stack([gdf["origin_x"].values, gdf["origin_y"].values])
    n = len(coords)

    if order is None:
        print(f"  ! '{id_col}' didn't parse as a clean, unique numeric "
              f"sequence -- falling back to greedy nearest-neighbor path")
        start_idx = int(np.argmin(coords[:, 1]))
        order = [start_idx]
        visited = np.zeros(n, dtype=bool)
        visited[start_idx] = True
        for _ in range(n - 1):
            current = coords[order[-1]]
            dists = np.hypot(coords[:, 0] - current[0], coords[:, 1] - current[1])
            dists[visited] = np.inf
            nxt = int(np.argmin(dists))
            order.append(nxt)
            visited[nxt] = True
        order = np.array(order)

    ordered_coords = coords[order]
    segment_lengths = np.hypot(np.diff(ordered_coords[:, 0]),
                                np.diff(ordered_coords[:, 1]))

    # Flag suspiciously large jumps in the run log
    if len(segment_lengths) > 0:
        typical = float(np.median(segment_lengths))
        bad = np.where(segment_lengths > max(20 * typical, 500))[0]
        if len(bad) > 0:
            print(f"  ! {len(bad)} suspiciously large alongshore jump(s) "
                  f"(>20x median spacing of {typical:.0f} m) -- inspect "
                  f"these transects before trusting far-field distances: "
                  f"{[str(gdf.iloc[order[b+1]].get('transect_id', order[b+1])) for b in bad[:10]]}"
                  f"{'...' if len(bad) > 10 else ''}")

    cumdist = np.concatenate([[0], np.cumsum(segment_lengths)])
    alongshore_map = dict(zip(order, cumdist))
    gdf["alongshore_m"] = [alongshore_map[i] for i in range(n)]
    return gdf


# Loader: CoastSat time-series CSVs

# {csv stem: path} for every CSV one folder level under root_dir
def collect_csv_map(root_dir: str) -> dict:
    csv_map = {}
    if not os.path.isdir(root_dir):
        print(f"  ! CoastSat root not found: {root_dir}")
        return csv_map
    for subfolder in sorted(os.listdir(root_dir)):
        sub_path = os.path.join(root_dir, subfolder)
        if not os.path.isdir(sub_path):
            continue
        for csv_file in sorted(glob.glob(os.path.join(sub_path, "*.csv"))):
            stem = Path(csv_file).stem
            csv_map[stem] = csv_file
    return csv_map


# CoastSat's date and chainage column names, whichever version wrote them
def _detect_coastsat_columns(df: pd.DataFrame) -> tuple:
    date_col = next((c for c in ("dates", "date", "time", "datetime")
                       if c in df.columns), None)
    chainage_col = next((c for c in ("chainage", "chainage_m",
                                       "shoreline_position", "distance",
                                       "distance_m") if c in df.columns), None)
    if date_col is None:
        date_col = df.columns[0]
    if chainage_col is None:
        chainage_col = df.columns[1] if len(df.columns) > 1 else None
    return date_col, chainage_col


# Every CoastSat per-transect CSV, as long-format chainage by transect ID
def load_coastsat_chainage(transects: gpd.GeoDataFrame,
                            root_dir: str) -> pd.DataFrame:
    print(f"\n[CoastSat chainage] reading CSVs: {root_dir}")
    csv_map = collect_csv_map(root_dir)
    print(f"  Found {len(csv_map)} CSVs in the data folder")
    tids = transects["transect_id"].astype(str).tolist()

    records = []
    n_found = 0
    n_missing = 0
    for i, tid in enumerate(tids):
        if (i + 1) % 100 == 0 or i == len(tids) - 1:
            print(f"    Reading CSVs: {i+1}/{len(tids)}", end="\r")
        csv_path = csv_map.get(tid)
        if csv_path is None:
            n_missing += 1
            continue
        try:
            df = pd.read_csv(csv_path)
        except Exception:
            n_missing += 1
            continue
        n_found += 1

        date_col, chainage_col = _detect_coastsat_columns(df)
        if chainage_col is None:
            continue
        df["_dt"] = pd.to_datetime(df[date_col], errors="coerce",
                                     utc=True).dt.tz_localize(None)
        df = df.dropna(subset=["_dt", chainage_col])
        for _, r in df.iterrows():
            records.append({
                "transect_id": tid,
                "date": r["_dt"],
                "chainage_m": float(r[chainage_col]),
                "source": "coastsat",
            })
    print()
    print(f"  Read {n_found} CSVs, {n_missing} missing")
    out = pd.DataFrame(records) if records else pd.DataFrame(
        columns=["transect_id", "date", "chainage_m", "source"])
    print(f"  → {len(out)} chainage observations")
    if len(out) > 0:
        print(f"  → date range: {out['date'].min().date()} → "
              f"{out['date'].max().date()}")
    return out


# Chainage by geometric intersection

# Chainage where each shoreline crosses each transect
def extract_chainage_by_intersection(shorelines: gpd.GeoDataFrame,
                                      transects: gpd.GeoDataFrame,
                                      source_key: str,
                                      transect_length_m: float = None,
                                      valid_range_m: tuple = None) -> pd.DataFrame:
    if transect_length_m is None:
        transect_length_m = TRANSECT_INTERSECTION_LENGTH_M
    if valid_range_m is None:
        valid_range_m = (-CHAINAGE_MAX_ABS_M, CHAINAGE_MAX_ABS_M)

    print(f"\n[{SOURCE_LABELS[source_key]} chainage] extracting via intersection...")
    print(f"  {len(shorelines)} shoreline features × {len(transects)} transects, "
          f"transect length={transect_length_m:.0f} m, "
          f"valid chainage range={valid_range_m}")

    if len(shorelines) == 0 or len(transects) == 0:
        return pd.DataFrame(columns=["transect_id", "date", "chainage_m", "source"])

    # Build transect LineStrings + STRtree for fast bbox pre-filter
    transect_lines = []
    transect_meta = []
    for _, row in transects.iterrows():
        ox, oy = row["origin_x"], row["origin_y"]
        dx, dy = row["dir_x"], row["dir_y"]
        end_x = ox + transect_length_m * dx
        end_y = oy + transect_length_m * dy
        line = LineString([(ox, oy), (end_x, end_y)])
        transect_lines.append(line)
        transect_meta.append({
            "transect_id": row["transect_id"],
            "origin_x":     ox,
            "origin_y":     oy,
            "dir_x":        dx,
            "dir_y":        dy,
        })

    tree = STRtree(transect_lines)
    # Map from id(line) → index into transect_meta
    line_to_idx = {id(ln): i for i, ln in enumerate(transect_lines)}

    records = []
    n_dropped_range = 0
    n_shorelines = len(shorelines)
    for i, (_, sh_row) in enumerate(shorelines.iterrows()):
        if (i + 1) % 20 == 0 or i == n_shorelines - 1:
            print(f"    Processing shoreline {i+1}/{n_shorelines}", end="\r")
        shore_geom = sh_row.geometry
        shore_date = sh_row["date"]

        # STRtree.query returns candidate line geometries
        candidates = tree.query(shore_geom)
        # Handle both API styles (indices in newer shapely, geoms in older)
        for cand in candidates:
            if isinstance(cand, (int, np.integer)):
                idx = int(cand)
                t_line = transect_lines[idx]
            else:
                idx = line_to_idx.get(id(cand))
                t_line = cand
                if idx is None:
                    continue

            try:
                intersect = shore_geom.intersection(t_line)
            except Exception:
                continue
            if intersect.is_empty:
                continue

            meta = transect_meta[idx]
            pt = _choose_intersection_point(intersect, meta)
            if pt is None:
                continue
            # Chainage = signed dot product along direction
            chainage = ((pt.x - meta["origin_x"]) * meta["dir_x"]
                        + (pt.y - meta["origin_y"]) * meta["dir_y"])
            if not (valid_range_m[0] <= chainage <= valid_range_m[1]):
                n_dropped_range += 1
                continue
            records.append({
                "transect_id": meta["transect_id"],
                "date":         shore_date,
                "chainage_m":   float(chainage),
                "source":       source_key,
            })
    print()
    print(f"  {n_dropped_range} intersections dropped as outside valid range")

    result = pd.DataFrame(records)
    if len(result) > 0:
        pct = np.percentile(result["chainage_m"].values, [1, 25, 50, 75, 99])
        print(f"  Chainage percentiles [1,25,50,75,99]: "
              f"{np.round(pct, 1).tolist()}")
        # Several features of one date on one transect (segmented NC lines): take the mean
        result = (result.groupby(["transect_id", "date", "source"])
                          .agg(chainage_m=("chainage_m", "mean"))
                          .reset_index())
    print(f"  → {len(result)} per-(transect, date) chainage observations")
    return result


# One crossing point from whatever geometry an intersection returns
def _choose_intersection_point(intersect, meta):
    origin = Point(meta["origin_x"], meta["origin_y"])
    if isinstance(intersect, Point):
        return intersect
    if isinstance(intersect, MultiPoint):
        pts = list(intersect.geoms)
        pts.sort(key=lambda p: origin.distance(p))
        return pts[0]
    if isinstance(intersect, LineString):
        return intersect.interpolate(0.5, normalized=True)
    if isinstance(intersect, MultiLineString):
        best = None
        best_dist = np.inf
        for ln in intersect.geoms:
            mid = ln.interpolate(0.5, normalized=True)
            d = origin.distance(mid)
            if d < best_dist:
                best_dist = d
                best = mid
        return best
    # GeometryCollection or other — try first Point-like element
    if hasattr(intersect, "geoms"):
        for g in intersect.geoms:
            if isinstance(g, Point):
                return g
        for g in intersect.geoms:
            if isinstance(g, (LineString, MultiLineString)):
                return _choose_intersection_point(g, meta)
    return None


# Unified chainage table

# All three sources in one table, with position, domain and decimal year
def build_unified_chainage_table(
    coastsat_chainage: pd.DataFrame,
    wet_dry_chainage:  pd.DataFrame,
    nc_state_chainage: pd.DataFrame,
    transects: gpd.GeoDataFrame,
) -> pd.DataFrame:
    all_chainage = pd.concat(
        [coastsat_chainage, wet_dry_chainage, nc_state_chainage],
        ignore_index=True,
    )
    if len(all_chainage) == 0:
        return all_chainage

    # Ensure transect_id is string on both sides
    all_chainage["transect_id"] = all_chainage["transect_id"].astype(str)
    tx = transects[["transect_id", "alongshore_m", "domain",
                    "origin_x", "origin_y"]].copy()
    tx["transect_id"] = tx["transect_id"].astype(str)

    merged = all_chainage.merge(tx, on="transect_id", how="inner")
    merged["decimal_year"] = decimal_year(merged["date"])
    merged["year"] = merged["date"].dt.year

    # Drop transects with too few observations
    counts = merged.groupby("transect_id").size()
    good = counts[counts >= MIN_OBSERVATIONS_PER_TRANSECT].index
    n_before = merged["transect_id"].nunique()
    merged = merged[merged["transect_id"].isin(good)].copy()
    n_after = merged["transect_id"].nunique()
    print(f"\n[merged] {len(merged)} total obs, "
          f"{n_after} transects kept "
          f"(dropped {n_before - n_after} with <{MIN_OBSERVATIONS_PER_TRANSECT} obs)")
    return merged


# Analysis: linear regression

# Unweighted OLS: slope (m/yr), intercept, r squared, n, slope SE
def regress_lrr(dec_years: np.ndarray, chainages: np.ndarray) -> dict:
    x = np.asarray(dec_years, dtype=float)
    y = np.asarray(chainages, dtype=float)
    mask = np.isfinite(x) & np.isfinite(y)
    x = x[mask]
    y = y[mask]
    n = len(x)
    if n < 2:
        return dict(slope=np.nan, intercept=np.nan, r_squared=np.nan,
                    n=n, se_slope=np.nan)
    mx = x.mean()
    my = y.mean()
    dx = x - mx
    dy = y - my
    denom = (dx * dx).sum()
    if denom == 0:
        return dict(slope=np.nan, intercept=np.nan, r_squared=np.nan,
                    n=n, se_slope=np.nan)
    slope = (dx * dy).sum() / denom
    intercept = my - slope * mx
    y_pred = slope * x + intercept
    ss_res = ((y - y_pred) ** 2).sum()
    ss_tot = ((y - my) ** 2).sum()
    r_squared = 1 - ss_res / ss_tot if ss_tot > 0 else np.nan
    if n > 2 and ss_res > 0:
        residual_var = ss_res / (n - 2)
        se_slope = np.sqrt(residual_var / denom)
    else:
        se_slope = np.nan
    return dict(slope=slope, intercept=intercept, r_squared=r_squared,
                n=n, se_slope=se_slope)


# Analysis: pre-installation LRR

# Pre-install LRR per transect; nourishment observations flagged, not dropped
def compute_preinstall_lrr(chainage: pd.DataFrame,
                            transects: gpd.GeoDataFrame) -> pd.DataFrame:
    print(f"\n{'='*72}\nPre-installation LRR")
    print(f"  Cutoff: all shorelines dated before {PRE_INSTALLATION_YEAR_CUTOFF}")

    # Everything before the cutoff, whichever pre-1970 years exist
    sub = chainage[chainage["year"] < PRE_INSTALLATION_YEAR_CUTOFF].copy()
    preinstall_years = sorted(sub["year"].unique().tolist())
    print(f"  Years present in data: {preinstall_years}")
    print(f"  → {len(sub)} observations across {sub['transect_id'].nunique()} "
          f"transects")
    print(f"  Per-source breakdown of pre-install observations:")
    for src, cnt in sub["source"].value_counts().items():
        print(f"    {SOURCE_LABELS.get(src, src):24s} {cnt}")

    # Flag, do not drop, observations in a documented nourishment zone
    def in_any_nourishment_zone(y):
        for name, year, y_min, y_max in NOURISHMENT_EXCLUSIONS:
            if year < PRE_INSTALLATION_YEAR_CUTOFF and y_min <= y <= y_max:
                return True
        return False

    for name, year, y_min, y_max in NOURISHMENT_EXCLUSIONS:
        if year >= PRE_INSTALLATION_YEAR_CUTOFF:
            continue
        n_flagged = int(((sub["year"] == year) &
                          (sub["origin_y"] >= y_min) &
                          (sub["origin_y"] <= y_max)).sum())
        print(f"  '{name}': {n_flagged} observations FLAGGED "
              f"(kept in the regression, marked nourishment_zone=True)")

    # Regress per transect (all data retained)
    tx_lookup = transects.set_index("transect_id")[["alongshore_m", "domain",
                                                      "origin_y",
                                                      "dist_from_groin_m"]].copy()
    tx_lookup.index = tx_lookup.index.astype(str)

    records = []
    for tid, grp in sub.groupby("transect_id"):
        fit = regress_lrr(grp["decimal_year"].values,
                          grp["chainage_m"].values)
        try:
            tx_row = tx_lookup.loc[str(tid)]
        except KeyError:
            continue
        records.append({
            "transect_id": tid,
            "alongshore_m": tx_row["alongshore_m"],
            "domain":       tx_row["domain"],
            "origin_y":     tx_row["origin_y"],
            "dist_from_groin_m": tx_row["dist_from_groin_m"],
            "n_obs":        fit["n"],
            "slope_m_yr":   fit["slope"],
            "intercept":    fit["intercept"],
            "r_squared":    fit["r_squared"],
            "se_slope":     fit["se_slope"],
            "nourishment_zone": in_any_nourishment_zone(tx_row["origin_y"]),
        })
    out = pd.DataFrame(records).sort_values("alongshore_m").reset_index(drop=True)
    n_good = int((out["n_obs"] >= 3).sum())
    print(f"  → {n_good}/{len(out)} transects have ≥3 pre-install observations "
          f"({int(out['nourishment_zone'].sum())} flagged nourishment-zone "
          f"transects included among them)")

    # Noise floor: far-field, unflagged pre-install spread, to compare with the signal threshold
    control = out[(out["n_obs"] >= 3) & (~out["nourishment_zone"]) &
                  (out["dist_from_groin_m"].abs() > SIGNAL_MAX_SEARCH_DISTANCE_M)]
    if len(control) >= 5:
        noise_sd = float(control["slope_m_yr"].std())
        print(f"\n  [noise floor] far-field pre-install slope std dev: "
              f"{noise_sd:.2f} m/yr  (n={len(control)} transects, "
              f"beyond {SIGNAL_MAX_SEARCH_DISTANCE_M/1000:.0f} km)")
        print(f"  → SIGNAL_ANOMALY_THRESHOLD_M_YR is currently "
              f"{SIGNAL_ANOMALY_THRESHOLD_M_YR:.2f} m/yr = "
              f"{SIGNAL_ANOMALY_THRESHOLD_M_YR / noise_sd:.2f}x that noise floor "
              f"(a common rule of thumb: use 1-2x the far-field std dev, "
              f"so ordinary variability isn't flagged as groin signal)")
    else:
        print(f"\n  [noise floor] too few far-field control transects "
              f"({len(control)}) to estimate -- widen the check or use a "
              f"different control definition")

    out.attrs["preinstall_years"] = preinstall_years
    return out


# Post-install LRR per transect (1970 on)
def compute_post_install_lrr(chainage: pd.DataFrame,
                              transects: gpd.GeoDataFrame) -> pd.DataFrame:
    print(f"\n{'='*72}\nPost-installation LRR")
    print(f"  Cutoff: all shorelines dated {GROIN_INSTALLATION_YEAR} or later")

    sub = chainage[chainage["year"] >= GROIN_INSTALLATION_YEAR].copy()
    print(f"  → {len(sub)} observations across {sub['transect_id'].nunique()} "
          f"transects")

    tx_lookup = transects.set_index("transect_id")[["alongshore_m", "domain",
                                                      "dist_from_groin_m"]].copy()
    tx_lookup.index = tx_lookup.index.astype(str)

    records = []
    for tid, grp in sub.groupby("transect_id"):
        fit = regress_lrr(grp["decimal_year"].values, grp["chainage_m"].values)
        try:
            tx_row = tx_lookup.loc[str(tid)]
        except KeyError:
            continue
        records.append({
            "transect_id": tid,
            "alongshore_m": tx_row["alongshore_m"],
            "domain":       tx_row["domain"],
            "dist_from_groin_m": tx_row["dist_from_groin_m"],
            "n_obs":        fit["n"],
            "slope_m_yr":   fit["slope"],
            "intercept":    fit["intercept"],
            "r_squared":    fit["r_squared"],
            "se_slope":     fit["se_slope"],
        })
    out = pd.DataFrame(records).sort_values("alongshore_m").reset_index(drop=True)
    n_good = int((out["n_obs"] >= 3).sum())
    print(f"  → {n_good}/{len(out)} transects have ≥3 post-install observations")
    return out


# Analysis: full-period LRR and piecewise breakpoint

# Full-period LRR per transect, and the best two-segment breakpoint
def compute_fullperiod_lrr(chainage: pd.DataFrame,
                            transects: gpd.GeoDataFrame) -> pd.DataFrame:
    print(f"\n{'='*72}\nFull-period LRR + piecewise breakpoint")

    tx_lookup = transects.set_index("transect_id")[["alongshore_m", "domain",
                                                      "origin_y",
                                                      "dist_from_groin_m"]].copy()
    tx_lookup.index = tx_lookup.index.astype(str)

    records = []
    all_tids = sorted(chainage["transect_id"].unique())
    for i, tid in enumerate(all_tids):
        if (i + 1) % 100 == 0 or i == len(all_tids) - 1:
            print(f"    Regressing transect {i+1}/{len(all_tids)}", end="\r")
        grp = chainage[chainage["transect_id"] == tid]
        if len(grp) < 4:
            continue
        # Single-slope fit
        full = regress_lrr(grp["decimal_year"].values,
                           grp["chainage_m"].values)
        # Piecewise fit
        pw = _piecewise_breakpoint(grp["decimal_year"].values,
                                    grp["chainage_m"].values)
        try:
            tx_row = tx_lookup.loc[str(tid)]
        except KeyError:
            continue
        records.append({
            "transect_id": tid,
            "alongshore_m": tx_row["alongshore_m"],
            "domain":       tx_row["domain"],
            "origin_y":     tx_row["origin_y"],
            "dist_from_groin_m": tx_row["dist_from_groin_m"],
            "n_obs":        full["n"],
            "slope_full":   full["slope"],
            "r2_full":      full["r_squared"],
            "breakpoint_year":    pw["breakpoint_year"],
            "slope_before":       pw["slope_before"],
            "slope_after":        pw["slope_after"],
            "rss_full":     pw["rss_full"],
            "rss_piecewise": pw["rss_piecewise"],
            "rss_reduction_frac": pw["rss_reduction_frac"],
        })
    print()
    return pd.DataFrame(records).sort_values("alongshore_m").reset_index(drop=True)


# Breakpoint year minimising the total RSS of a two-segment fit
def _piecewise_breakpoint(dec_years: np.ndarray, chainages: np.ndarray) -> dict:
    x = np.asarray(dec_years, dtype=float)
    y = np.asarray(chainages, dtype=float)
    mask = np.isfinite(x) & np.isfinite(y)
    x = x[mask]
    y = y[mask]
    if len(x) < 2 * BREAKPOINT_MIN_POINTS_PER_SEGMENT:
        return dict(breakpoint_year=np.nan, slope_before=np.nan,
                    slope_after=np.nan, rss_full=np.nan,
                    rss_piecewise=np.nan, rss_reduction_frac=np.nan)
    # Single-slope RSS as baseline
    full_fit = regress_lrr(x, y)
    y_full_pred = full_fit["slope"] * x + full_fit["intercept"]
    rss_full = float(((y - y_full_pred) ** 2).sum())

    best_rss = np.inf
    best = dict(breakpoint_year=np.nan, slope_before=np.nan,
                slope_after=np.nan)
    for bp in range(BREAKPOINT_SEARCH_START, BREAKPOINT_SEARCH_END + 1):
        m_before = x < bp
        m_after  = x >= bp
        if m_before.sum() < BREAKPOINT_MIN_POINTS_PER_SEGMENT:
            continue
        if m_after.sum() < BREAKPOINT_MIN_POINTS_PER_SEGMENT:
            continue
        fit_b = regress_lrr(x[m_before], y[m_before])
        fit_a = regress_lrr(x[m_after],  y[m_after])
        if not (np.isfinite(fit_b["slope"]) and np.isfinite(fit_a["slope"])):
            continue
        pred_b = fit_b["slope"] * x[m_before] + fit_b["intercept"]
        pred_a = fit_a["slope"] * x[m_after]  + fit_a["intercept"]
        rss = float(((y[m_before] - pred_b) ** 2).sum()
                    + ((y[m_after]  - pred_a) ** 2).sum())
        if rss < best_rss:
            best_rss = rss
            best = dict(breakpoint_year=bp,
                        slope_before=fit_b["slope"],
                        slope_after=fit_a["slope"])
    if not np.isfinite(best_rss):
        return dict(breakpoint_year=np.nan, slope_before=np.nan,
                    slope_after=np.nan, rss_full=rss_full,
                    rss_piecewise=np.nan, rss_reduction_frac=np.nan)
    frac_reduction = (rss_full - best_rss) / rss_full if rss_full > 0 else np.nan
    return dict(
        breakpoint_year=best["breakpoint_year"],
        slope_before=best["slope_before"],
        slope_after=best["slope_after"],
        rss_full=rss_full,
        rss_piecewise=best_rss,
        rss_reduction_frac=frac_reduction,
    )


# Analysis: decadal LRR (fixed 10-yr bins)

# LRR per transect in fixed, non-overlapping decades
def compute_decadal_lrr(chainage: pd.DataFrame) -> pd.DataFrame:
    print(f"\n{'='*72}\nDecadal LRR ({DECADE_LENGTH_YEARS}-yr non-overlapping bins)")
    last_year = int(np.ceil(chainage["year"].max())) if len(chainage) else DECADE_START_YEAR
    decade_starts = list(range(DECADE_START_YEAR, last_year + 1, DECADE_LENGTH_YEARS))
    print(f"  Decades: {decade_starts[0]}s … {decade_starts[-1]}s "
          f"({DECADE_LENGTH_YEARS} yr each, non-overlapping), "
          f"min obs {DECADE_MIN_OBSERVATIONS}")

    records = []
    all_tids = sorted(chainage["transect_id"].unique())
    for i, tid in enumerate(all_tids):
        if (i + 1) % 100 == 0 or i == len(all_tids) - 1:
            print(f"    Transect {i+1}/{len(all_tids)}", end="\r")
        grp = chainage[chainage["transect_id"] == tid]
        dy = grp["decimal_year"].values
        yr = grp["year"].values
        ch = grp["chainage_m"].values
        for d0 in decade_starts:
            d1 = d0 + DECADE_LENGTH_YEARS
            mask = (yr >= d0) & (yr < d1)
            if mask.sum() < DECADE_MIN_OBSERVATIONS:
                continue
            fit = regress_lrr(dy[mask], ch[mask])
            records.append({
                "transect_id":  tid,
                "decade_start":  d0,          # decade START year
                "decade_label": f"{d0}s",
                "n_obs":        fit["n"],
                "slope_m_yr":   fit["slope"],
                "r_squared":    fit["r_squared"],
            })
    print()
    result = pd.DataFrame(records)
    print(f"  → {len(result)} (transect, decade) fits produced across "
          f"{result['decade_start'].nunique() if len(result) else 0} decades")
    return result


# Analysis: per-era LRR

# LRR per transect per era
def compute_era_lrrs(chainage: pd.DataFrame,
                      transects: gpd.GeoDataFrame,
                      extra_periods: list = None) -> pd.DataFrame:
    print(f"\n{'='*72}\nPer-era LRR (one rate per transect per era)")
    tx = transects.set_index("transect_id")[["alongshore_m", "domain",
                                                "origin_y",
                                                "dist_from_groin_m"]]
    tx.index = tx.index.astype(str)

    periods_to_compute = list(ERAS) + list(extra_periods or [])

    records = []
    all_tids = sorted(chainage["transect_id"].unique())
    for i, tid in enumerate(all_tids):
        if (i + 1) % 200 == 0 or i == len(all_tids) - 1:
            print(f"    Transect {i+1}/{len(all_tids)}", end="\r")
        grp = chainage[chainage["transect_id"] == tid]
        try:
            tx_row = tx.loc[str(tid)]
        except KeyError:
            continue
        for era_name, y_lo, y_hi in periods_to_compute:
            sub = grp[(grp["year"] >= y_lo) & (grp["year"] <= y_hi)]
            if len(sub) < 3:
                continue
            fit = regress_lrr(sub["decimal_year"].values,
                               sub["chainage_m"].values)
            records.append({
                "transect_id":       tid,
                "alongshore_m":      tx_row["alongshore_m"],
                "dist_from_groin_m": tx_row["dist_from_groin_m"],
                "domain":            tx_row["domain"],
                "origin_y":          tx_row["origin_y"],
                "era":               era_name,
                "y_lo":              y_lo,
                "y_hi":              y_hi,
                "n_obs":             fit["n"],
                "slope_m_yr":        fit["slope"],
                "intercept":         fit["intercept"],
                "r_squared":         fit["r_squared"],
                "se_slope":          fit["se_slope"],
            })
    print()
    out = pd.DataFrame(records)
    if len(out) > 0:
        summary = out.groupby("era")["transect_id"].nunique()
        for era_name, n_tx in summary.items():
            n_total = len(out[out["era"] == era_name])
            print(f"  {era_name:>15s}: {n_tx} transects, {n_total} fits")
    return out


# Regional pre-install baseline: per-transect values joined by linear interpolation
def interp_preinstall_baseline(preinstall_lrr: pd.DataFrame) -> callable:
    good = preinstall_lrr[(preinstall_lrr["n_obs"] >= 3) &
                           preinstall_lrr["slope_m_yr"].notna() &
                           (~preinstall_lrr["nourishment_zone"])].copy()
    if len(good) < 2:
        print(f"  ! Only {len(good)} transects with usable pre-install LRR; "
              f"baseline will be a flat constant (median)")
        median = good["slope_m_yr"].median() if len(good) else 0.0
        return lambda x: np.full_like(np.asarray(x, dtype=float),
                                        float(median))

    good = good.sort_values("dist_from_groin_m")
    xs = good["dist_from_groin_m"].values
    ys = good["slope_m_yr"].values

    def baseline(query_dist_from_groin_m):
        q = np.asarray(query_dist_from_groin_m, dtype=float)
        return np.interp(q, xs, ys, left=ys[0], right=ys[-1])

    print(f"  Built connect-the-dots baseline (linear interpolation, "
          f"NOT binned/smoothed) over {len(good)} pre-install transects, "
          f"median = {np.median(ys):+.2f} m/yr")
    return baseline


# Analysis: decadal anomaly against the baseline

# Decadal LRR minus the pre-install baseline, with distance from the groin
def compute_decadal_anomaly(decadal_lrr: pd.DataFrame,
                             transects: gpd.GeoDataFrame,
                             baseline_fn) -> pd.DataFrame:
    print(f"\n{'='*72}\nDecadal anomaly relative to baseline")
    tx = transects[["transect_id", "alongshore_m", "dist_from_groin_m",
                     "origin_y", "domain"]].copy()
    tx["transect_id"] = tx["transect_id"].astype(str)
    merged = decadal_lrr.merge(tx, on="transect_id", how="inner")
    merged["baseline_lrr"] = baseline_fn(merged["dist_from_groin_m"].values)
    merged["anomaly_m_yr"] = merged["slope_m_yr"] - merged["baseline_lrr"]
    print(f"  {len(merged)} decadal anomaly rows across "
          f"{merged['transect_id'].nunique()} transects and "
          f"{merged['decade_start'].nunique()} decades")
    return merged


# Analysis: signal extent over time

# How far the signal runs outward from the groin on one side
def _contiguous_extent_one_side(dist_abs: np.ndarray, anomaly: np.ndarray,
                                 thresh: float, bin_w: float, max_d: float,
                                 max_gap_bins: int, direction: int) -> tuple:
    n_bins = int(np.ceil(max_d / bin_w))
    if len(dist_abs) == 0:
        return 0.0, np.nan
    bin_idx = np.clip((dist_abs // bin_w).astype(int), 0, n_bins - 1)
    bin_medians = np.full(n_bins, np.nan)
    for b in range(n_bins):
        vals = anomaly[bin_idx == b]
        if len(vals) > 0:
            bin_medians[b] = np.median(vals)

    extent_bins = 0
    gap = 0
    for b in range(n_bins):
        v = bin_medians[b]
        exceeds = np.isfinite(v) and (direction * v > thresh)
        if exceeds:
            extent_bins = b + 1
            gap = 0
        else:
            gap += 1
            if gap > max_gap_bins:
                break
    extent_m = extent_bins * bin_w

    if extent_bins > 0:
        within = bin_medians[:extent_bins]
        within = within[np.isfinite(within)]
        peak = float(within.max() if direction > 0 else within.min()) \
               if len(within) else np.nan
    else:
        peak = np.nan
    return extent_m, peak


# Signal extent updrift and downdrift, per decade
def compute_signal_extent_over_time(decadal_anomaly: pd.DataFrame) -> pd.DataFrame:
    print(f"\n{'='*72}\nSignal extent over time")
    print(f"  Threshold: |anomaly| ≥ {SIGNAL_ANOMALY_THRESHOLD_M_YR:.2f} m/yr, "
          f"contiguous from groin, bin={SIGNAL_EXTENT_BIN_WIDTH_M:.0f} m, "
          f"max gap={SIGNAL_EXTENT_MAX_GAP_BINS} bin(s)")
    print(f"  Max search distance: {SIGNAL_MAX_SEARCH_DISTANCE_M/1000:.1f} km "
          f"either side of groin")

    records = []
    thresh = SIGNAL_ANOMALY_THRESHOLD_M_YR
    max_d = SIGNAL_MAX_SEARCH_DISTANCE_M
    bin_w = SIGNAL_EXTENT_BIN_WIDTH_M
    max_gap = SIGNAL_EXTENT_MAX_GAP_BINS
    for c, grp in decadal_anomaly.groupby("decade_start"):
        grp = grp[grp["dist_from_groin_m"].abs() <= max_d]
        up = grp[grp["dist_from_groin_m"] > 0]
        dn = grp[grp["dist_from_groin_m"] < 0]

        updrift_extent, updrift_peak = _contiguous_extent_one_side(
            up["dist_from_groin_m"].values, up["anomaly_m_yr"].values,
            thresh, bin_w, max_d, max_gap, direction=+1)
        downdrift_extent, downdrift_peak = _contiguous_extent_one_side(
            dn["dist_from_groin_m"].abs().values, dn["anomaly_m_yr"].values,
            thresh, bin_w, max_d, max_gap, direction=-1)

        records.append({
            "decade_start":            int(c),
            "updrift_extent_m":       updrift_extent,
            "updrift_peak_anomaly":   updrift_peak,
            "downdrift_extent_m":     downdrift_extent,
            "downdrift_peak_anomaly": downdrift_peak,
            "n_updrift_transects":   int(len(up)),
            "n_downdrift_transects": int(len(dn)),
        })
    out = pd.DataFrame(records).sort_values("decade_start").reset_index(drop=True)
    return out


# Plots

# Figure: pre-install and full-period LRR against distance from the groin
def plot_alongshore_lrr_profile(preinstall_lrr: pd.DataFrame,
                                 fullperiod_lrr: pd.DataFrame,
                                 postinstall_lrr: pd.DataFrame,
                                 transects: gpd.GeoDataFrame,
                                 groin: dict,
                                 output_path: str):
    print(f"\n[plot] alongshore LRR profile → {output_path}")

    fig, ax = plt.subplots(figsize=(14, 7))

    # Full-period per-transect points, connected in order
    fs = fullperiod_lrr[fullperiod_lrr["slope_full"].notna()].copy()
    fs = fs.sort_values("dist_from_groin_m")
    if len(fs) > 0:
        ax.plot(fs["dist_from_groin_m"] / 1000, fs["slope_full"],
                 color="#333", lw=0.7, alpha=0.45, marker=".", ms=4,
                 markeredgewidth=0, zorder=2,
                 label=f"Full-period per-transect LRR  "
                        f"(n={len(fs)} transects)")

    # ─── Post-install (1970-2024) per-transect points + connecting line ───
    ps = postinstall_lrr[postinstall_lrr["slope_m_yr"].notna()].copy()
    ps = ps.sort_values("dist_from_groin_m")
    if len(ps) > 0:
        ax.plot(ps["dist_from_groin_m"] / 1000, ps["slope_m_yr"],
                 color="#1565C0", lw=0.7, alpha=0.45, marker=".", ms=4,
                 markeredgewidth=0, zorder=2,
                 label=f"Post-install ({GROIN_INSTALLATION_YEAR}-present) "
                        f"per-transect LRR  (n={len(ps)} transects)")

    # Pre-install per-transect points; nourishment-flagged ones get the same marker
    ok_pre = preinstall_lrr[(preinstall_lrr["n_obs"] >= 3) &
                             preinstall_lrr["slope_m_yr"].notna()].sort_values(
        "dist_from_groin_m")
    preinstall_years = preinstall_lrr.attrs.get("preinstall_years", [])
    yr_range_str = (f"{min(preinstall_years)}–{max(preinstall_years)}"
                    if preinstall_years else "no years found")
    if len(ok_pre) > 0:
        ax.plot(ok_pre["dist_from_groin_m"] / 1000, ok_pre["slope_m_yr"],
                 color=ERA_COLORS["Pre-install"], lw=1.0, alpha=0.6,
                 marker=".", ms=6, markeredgewidth=0, zorder=5,
                 label=f"Pre-install baseline  (per-transect LRR, "
                        f"n={len(ok_pre)}, {yr_range_str}, "
                        f"years={preinstall_years})")

    # Reference lines
    ax.axhline(0, color="#999", lw=0.5, linestyle=":")
    ax.axvline(0, color="black", lw=1.0, alpha=0.7,
                label=f"Groin  (updrift = "
                       f"{'north' if UPDRIFT_DIRECTION == 'north' else 'south'})")

    # Groin field band, south end to north end; not centred on 0, which is the north end
    if groin.get("geometry") is not None:
        sign = 1 if UPDRIFT_DIRECTION.lower().startswith("n") else -1
        g_a = sign * (groin["y_min"] - groin["reference_y"]) / 1000
        g_b = sign * (groin["y_max"] - groin["reference_y"]) / 1000
        ax.axvspan(min(g_a, g_b), max(g_a, g_b),
                    color="black", alpha=0.12, zorder=0,
                    label="Groin field  (first to last groin)")

    # x-axis from domain 1 to PLOT_UPDRIFT_MAX_KM updrift
    d1_km = domain_dist_km(transects, 1)
    apply_distance_xlim_and_autoscale_y(ax, x_max_km=PLOT_UPDRIFT_MAX_KM,
                                          x_min_km=d1_km)
    x_lo, x_hi = ax.get_xlim()
    set_domain_primary_axis(ax, transects)
    add_geographic_annotations(ax, transects)

    ax.set_ylabel("LRR  (m/yr, + = accretion)", fontsize=FONT_AXIS_LABEL)
    ax.set_title(
        "Shoreline Change Rate vs. Distance from Groin",
        fontsize=FONT_TITLE, fontweight="bold", pad=10,
    )
    ax.grid(True, alpha=0.25, linewidth=0.4)
    handles, labels = ax.get_legend_handles_labels()
    handles += annotation_legend_handles()
    ax.legend(handles=handles, loc="lower right", fontsize=FONT_LEGEND, framealpha=0.92)

    plt.tight_layout()
    plt.savefig(output_path, dpi=200, bbox_inches="tight")
    plt.close(fig)


_LOWESS_WARNED = False   # print the "statsmodels missing" warning only once


# LOWESS-smoothed (x, y) for an overlay, or (None, None) without statsmodels
def _lowess_overlay(x: np.ndarray, y: np.ndarray, frac: float):
    global _LOWESS_WARNED
    try:
        from statsmodels.nonparametric.smoothers_lowess import lowess
    except ImportError:
        if not _LOWESS_WARNED:
            print("  ! statsmodels not installed (pip install statsmodels) "
                  "-- skipping smoothed overlay lines")
            _LOWESS_WARNED = True
        return None, None
    mask = np.isfinite(x) & np.isfinite(y)
    x, y = x[mask], y[mask]
    if len(x) < 10:
        return None, None
    order = np.argsort(x)
    smoothed = lowess(y[order], x[order], frac=frac, return_sorted=True)
    return smoothed[:, 0], smoothed[:, 1]


# Figure: LRR profile, one colour per era, with smoothed overlays
def plot_era_lrr_profile(era_lrrs: pd.DataFrame,
                          preinstall_lrr: pd.DataFrame,
                          transects: gpd.GeoDataFrame,
                          groin: dict,
                          output_path: str):
    print(f"\n[plot] per-era LRR profile → {output_path}")
    if len(era_lrrs) == 0:
        print("  ! no era LRR data")
        return

    fig, ax = plt.subplots(figsize=(14, 7.5))

    # Pre-install baseline as reference
    ok_pre = preinstall_lrr[(preinstall_lrr["n_obs"] >= 3) &
                             preinstall_lrr["slope_m_yr"].notna()].sort_values(
        "dist_from_groin_m")
    preinstall_years = preinstall_lrr.attrs.get("preinstall_years", [])
    yr_range_str = (f"{min(preinstall_years)}–{max(preinstall_years)}"
                    if preinstall_years else "no years found")
    pre_color = ERA_COLORS["Pre-install"]
    if len(ok_pre) > 0:
        ax.plot(ok_pre["dist_from_groin_m"] / 1000, ok_pre["slope_m_yr"],
                 color=pre_color, lw=0.6, alpha=0.35, marker=".", ms=4,
                 markeredgewidth=0, zorder=3,
                 label=f"Pre-install baseline  (per-transect LRR, "
                        f"n={len(ok_pre)}, {yr_range_str})")
        xs, ys = _lowess_overlay(ok_pre["dist_from_groin_m"].values / 1000,
                                   ok_pre["slope_m_yr"].values,
                                   ERA_PROFILE_LOWESS_FRAC)
        if xs is not None:
            ax.plot(xs, ys, color=pre_color, lw=2.8, zorder=7)

    # Post-install eras and the pre-nourishment sub-period: faint points, bold smoothed line
    periods_to_plot = [(n, lo, hi) for n, lo, hi in ERAS if n != "Pre-install"]
    periods_to_plot.append(PRE_NOURISHMENT_PERIOD)
    for era_name, y_lo, y_hi in periods_to_plot:
        sub = era_lrrs[era_lrrs["era"] == era_name].sort_values("dist_from_groin_m")
        if len(sub) < 5:
            continue
        is_pre_nourishment = (era_name == PRE_NOURISHMENT_PERIOD[0])
        color = PRE_NOURISHMENT_COLOR if is_pre_nourishment else ERA_COLORS.get(era_name, "#333")

        ax.plot(sub["dist_from_groin_m"] / 1000, sub["slope_m_yr"],
                 color=color, lw=0.6, alpha=0.35, marker=".", ms=3.5,
                 markeredgewidth=0, zorder=3,
                 label=f"{era_name}  ({y_lo}–{y_hi}, n={len(sub)})")
        xs, ys = _lowess_overlay(sub["dist_from_groin_m"].values / 1000,
                                   sub["slope_m_yr"].values,
                                   ERA_PROFILE_LOWESS_FRAC)
        if xs is not None:
            # Dash the pre-nourishment line so it does not blur into Deteriorated
            ax.plot(xs, ys, color=color, lw=2.8, zorder=6,
                     linestyle="--" if is_pre_nourishment else "-")

    # Reference lines
    ax.axhline(0, color="#999", lw=0.5, linestyle=":")
    ax.axvline(0, color="black", lw=1.0, alpha=0.7,
                label=f"Groin  (updrift = "
                       f"{'north' if UPDRIFT_DIRECTION == 'north' else 'south'})")
    if groin.get("geometry") is not None:
        sign = 1 if UPDRIFT_DIRECTION.lower().startswith("n") else -1
        g_a = sign * (groin["y_min"] - groin["reference_y"]) / 1000
        g_b = sign * (groin["y_max"] - groin["reference_y"]) / 1000
        ax.axvspan(min(g_a, g_b), max(g_a, g_b),
                    color="black", alpha=0.12, zorder=0,
                    label="Groin field  (first to last groin)")

    d1_km = domain_dist_km(transects, 1)
    apply_distance_xlim_and_autoscale_y(ax, x_max_km=PLOT_UPDRIFT_MAX_KM,
                                          x_min_km=d1_km)
    x_lo, x_hi = ax.get_xlim()
    set_domain_primary_axis(ax, transects)
    add_geographic_annotations(ax, transects)
    ax.set_ylabel("LRR  (m/yr, + = accretion)", fontsize=FONT_AXIS_LABEL)
    ax.set_title(
        "Shoreline Change Rate by Structural Era",
        fontsize=FONT_TITLE, fontweight="bold", pad=10,
    )
    ax.grid(True, alpha=0.25, linewidth=0.4)
    handles, labels = ax.get_legend_handles_labels()
    handles += annotation_legend_handles()
    ax.legend(handles=handles, loc="best", fontsize=FONT_LEGEND, framealpha=0.92, ncol=1)

    plt.tight_layout()
    plt.savefig(output_path, dpi=200, bbox_inches="tight")
    plt.close(fig)


# Non-overlapping (name, start, end) windows from the start to the end year
def _build_decade_periods(start_year: int, increment_years: int,
                           end_year: int) -> list:
    periods = []
    y = start_year
    while y < end_year:
        y_hi = min(y + increment_years - 1, end_year)
        periods.append((f"{y}-{y_hi}", y, y_hi))
        y += increment_years
    return periods


# Figure: LRR profile per fixed window since installation
def plot_decade_lrr_profile(decade_period_lrrs: pd.DataFrame,
                             transects: gpd.GeoDataFrame,
                             groin: dict,
                             periods: list,
                             increment_years: int,
                             output_path: str):
    print(f"\n[plot] decade-increment LRR profile ({increment_years}-yr) → {output_path}")
    fig, ax = plt.subplots(figsize=(14, 7.5))

    n = len(periods)
    cmap = plt.get_cmap(DECADE_PLOT_COLORMAP)
    # Sample the colormap from 0.45 so the earliest line still reads on white
    colors = [cmap(0.45 + 0.55 * (i / max(n - 1, 1))) for i in range(n)]

    any_plotted = False
    for i, (period_name, y_lo, y_hi) in enumerate(periods):
        sub = decade_period_lrrs[decade_period_lrrs["era"] == period_name] \
            .sort_values("dist_from_groin_m")
        if len(sub) < 5:
            print(f"  ! {period_name}: only {len(sub)} transects -- skipping line")
            continue
        color = colors[i]
        xs, ys = _lowess_overlay(sub["dist_from_groin_m"].values / 1000,
                                   sub["slope_m_yr"].values,
                                   SMOOTHED_OVERLAY_LOWESS_FRAC)
        if xs is not None:
            ax.plot(xs, ys, color=color, lw=1.5, zorder=6,
                     label=f"{period_name}  (n={len(sub)})")
        any_plotted = True

    if not any_plotted:
        print("  ! no period had enough data -- skipping plot")
        plt.close(fig)
        return

    ax.axhline(0, color="#999", lw=0.5, linestyle=":")
    ax.axvline(0, color="black", lw=1.0, alpha=0.7,
                label=f"Groin  (updrift = "
                       f"{'north' if UPDRIFT_DIRECTION == 'north' else 'south'})")
    if groin.get("geometry") is not None:
        sign = 1 if UPDRIFT_DIRECTION.lower().startswith("n") else -1
        g_a = sign * (groin["y_min"] - groin["reference_y"]) / 1000
        g_b = sign * (groin["y_max"] - groin["reference_y"]) / 1000
        ax.axvspan(min(g_a, g_b), max(g_a, g_b),
                    color="black", alpha=0.12, zorder=0,
                    label="Groin field  (first to last groin)")

    d1_km = domain_dist_km(transects, 1)
    apply_distance_xlim_and_autoscale_y(ax, x_max_km=PLOT_UPDRIFT_MAX_KM,
                                          x_min_km=d1_km)
    x_lo, x_hi = ax.get_xlim()
    set_domain_primary_axis(ax, transects)
    add_geographic_annotations(ax, transects)
    ax.set_ylabel("LRR  (m/yr, + = accretion)", fontsize=FONT_AXIS_LABEL)
    ax.set_title(
        f"Shoreline Change Rate — {increment_years}-Year Windows\n"
        f"Color: light → dark = {periods[0][0].split('-')[0]} → present",
        fontsize=FONT_TITLE, fontweight="bold", pad=10,
    )
    ax.grid(True, alpha=0.25, linewidth=0.4)
    handles, labels = ax.get_legend_handles_labels()
    handles += annotation_legend_handles()
    ax.legend(handles=handles, loc="best", fontsize=FONT_LEGEND, framealpha=0.92, ncol=2)

    plt.tight_layout()
    plt.savefig(output_path, dpi=200, bbox_inches="tight")
    plt.close(fig)


# A zone's mean rate as a heavy segment with a label
def _zone_mean_and_label(ax, df: pd.DataFrame, zone_domains: tuple,
                          color: str, domain_to_dist):
    d_lo, d_hi = zone_domains
    km_lo, km_hi = float(domain_to_dist(d_lo)), float(domain_to_dist(d_hi))
    km_lo, km_hi = min(km_lo, km_hi), max(km_lo, km_hi)
    mask = ((df["dist_from_groin_m"] / 1000 >= km_lo) &
            (df["dist_from_groin_m"] / 1000 <= km_hi))
    if mask.sum() == 0:
        return
    zone_mean = df.loc[mask, "slope_m_yr"].mean()
    ax.hlines(zone_mean, km_lo, km_hi, colors=color, lw=2.5, alpha=0.85, zorder=6)
    ax.annotate(f"{zone_mean:+.1f}", xy=((km_lo + km_hi) / 2, zone_mean),
                 xytext=(0, 6 if zone_mean >= 0 else -12),
                 textcoords="raw_offset points", ha="center",
                 fontsize=8, color=color, fontweight="bold", zorder=7)


# Figure: one panel per era, analysis zones shaded with mean rates
def plot_zone_panels(preinstall_lrr: pd.DataFrame, era_lrrs: pd.DataFrame,
                      transects: gpd.GeoDataFrame, groin: dict,
                      output_path: str,
                      extra_functional_panels: list = None):
    print(f"\n[plot] zone panels → {output_path}")
    domain_to_dist, _ = build_domain_to_dist_km(transects)

    ok_pre = preinstall_lrr[(preinstall_lrr["n_obs"] >= 3) &
                             preinstall_lrr["slope_m_yr"].notna()]
    panels = [
        ("Pre-install", ok_pre),
        ("Functional groin", era_lrrs[era_lrrs["era"] == "Functional groin"]),
    ]
    panels.extend(extra_functional_panels or [])
    panels.append(("Deteriorated", era_lrrs[era_lrrs["era"] == "Deteriorated"]))
    n_panels = len(panels)

    d1_km = domain_dist_km(transects, 1)
    x_min_km = d1_km if d1_km is not None else -3.0
    x_max_km = float(domain_to_dist(ZONE_PANEL_X_MAX_DOMAIN))

    fig, axes = plt.subplots(n_panels, 1, figsize=(12, 3.4 * n_panels),
                               sharex=True)
    if n_panels == 1:
        axes = [axes]

    preinstall_df = panels[0][1].sort_values("dist_from_groin_m")

    zones = [
        ("Downdrift zone", ZONE_PANEL_DOWNDRIFT_DOMAINS, ZONE_PANEL_COLOR_DOWNDRIFT),
        ("Updrift zone", ZONE_PANEL_UPDRIFT_DOMAINS, ZONE_PANEL_COLOR_UPDRIFT),
    ]

    for i, (era_name, df) in enumerate(panels):
        ax = axes[i]
        if era_name in ERA_COLORS:
            color = ERA_COLORS[era_name]
        elif "Functional" in era_name:
            color = lighten_color(ERA_COLORS["Functional groin"], amount=0.35)
        else:
            color = "#333"
        if len(df) == 0:
            ax.text(0.5, 0.5, f"{era_name}: no data", transform=ax.transAxes,
                     ha="center", va="center")
            continue
        df_sorted = df.sort_values("dist_from_groin_m")

        for zname, zdomains, zcolor in zones:
            zd_lo, zd_hi = zdomains
            zkm_lo, zkm_hi = float(domain_to_dist(zd_lo)), float(domain_to_dist(zd_hi))
            zkm_lo, zkm_hi = min(zkm_lo, zkm_hi), max(zkm_lo, zkm_hi)
            ax.axvspan(zkm_lo, zkm_hi, color=zcolor, alpha=0.12, zorder=0,
                        label=zname if i == 0 else None)

        ax.axhline(0, color="k", lw=0.5, zorder=2)
        ax.axvline(0, color="black", lw=1.0, alpha=0.7, zorder=2,
                    label="Groin" if i == 0 else None)

        # Pre-install reference line (faint), on every panel except its own
        if era_name != "Pre-install" and len(preinstall_df) > 0:
            xs_pc, ys_pc = _lowess_overlay(
                preinstall_df["dist_from_groin_m"].values / 1000,
                preinstall_df["slope_m_yr"].values,
                ERA_PROFILE_LOWESS_FRAC)
            if xs_pc is not None:
                ax.plot(xs_pc, ys_pc, color=ERA_COLORS["Pre-install"], lw=1.2,
                         ls="--", alpha=0.7, zorder=3,
                         label="Pre-install (smoothed)" if i == 1 else None)

        ax.plot(df_sorted["dist_from_groin_m"] / 1000, df_sorted["slope_m_yr"],
                 marker="o", ms=2.5, lw=0.8, color=color, alpha=0.6, zorder=4)

        xs, ys = _lowess_overlay(df_sorted["dist_from_groin_m"].values / 1000,
                                   df_sorted["slope_m_yr"].values,
                                   ERA_PROFILE_LOWESS_FRAC)
        if xs is not None:
            ax.plot(xs, ys, color=color, lw=2.2, alpha=0.95, zorder=5,
                     label="Smoothed" if i == 0 else None)

        ax.set_xlim(x_min_km, x_max_km)
        ax.set_ylabel(f"{era_name}\n(m/yr)", fontsize=FONT_ANNOTATION)
        ax.grid(True, alpha=0.3)

        for _, zdomains, zcolor in zones:
            _zone_mean_and_label(ax, df_sorted, zdomains, zcolor, domain_to_dist)

        if i == 0:
            domain_ticks = np.arange(1, ZONE_PANEL_X_MAX_DOMAIN + 1, 2)
            tick_km = domain_to_dist(domain_ticks)
            secax = ax.secondary_xaxis("top")
            secax.set_xticks(tick_km)
            secax.set_xticklabels([f"D{int(d)}" for d in domain_ticks],
                                   fontsize=FONT_TICK - 2)

    axes[-1].set_xlabel("Distance from groin (km, alongshore)  [+ = updrift]",
                          fontsize=FONT_AXIS_LABEL)
    fig.suptitle("Alongshore Shoreline Change Rate by Era, with Analysis Zones",
                  fontsize=FONT_TITLE, fontweight="bold", y=1.0)

    handles, labels = axes[0].get_legend_handles_labels()
    if len(axes) > 1:
        h2, l2 = axes[1].get_legend_handles_labels()
        for h, l in zip(h2, l2):
            if l not in labels:
                handles.append(h)
                labels.append(l)
    fig.legend(handles, labels, loc="upper center", bbox_to_anchor=(0.5, 0.955),
                ncol=len(handles), fontsize=FONT_LEGEND - 1, frameon=True)

    plt.tight_layout(rect=[0, 0, 1, 0.90])
    plt.savefig(output_path, dpi=200, bbox_inches="tight")
    print(f"  → saved {output_path}")
    plt.close(fig)


# Figure: change in LRR between consecutive eras, per transect
def plot_era_difference_profile(era_lrrs: pd.DataFrame,
                                 preinstall_lrr: pd.DataFrame,
                                 transects: gpd.GeoDataFrame,
                                 groin: dict,
                                 output_path: str):
    print(f"\n[plot] era difference profile → {output_path}")

    pre = preinstall_lrr[(preinstall_lrr["n_obs"] >= 3) &
                          preinstall_lrr["slope_m_yr"].notna()][
        ["transect_id", "dist_from_groin_m", "slope_m_yr"]].copy()
    func = era_lrrs[era_lrrs["era"] == "Functional groin"][
        ["transect_id", "dist_from_groin_m", "slope_m_yr"]].copy()
    det = era_lrrs[era_lrrs["era"] == "Deteriorated"][
        ["transect_id", "dist_from_groin_m", "slope_m_yr"]].copy()

    if len(pre) == 0 or len(func) == 0 or len(det) == 0:
        print("  ! missing one of pre-install/functional/deteriorated -- skipping")
        return

    pre["transect_id"] = pre["transect_id"].astype(str)
    func["transect_id"] = func["transect_id"].astype(str)
    det["transect_id"] = det["transect_id"].astype(str)

    print("  [diagnostic] Pre-install vs Functional groin:")
    _diagnose_era_overlap("Pre-install", pre["transect_id"],
                           "Functional groin", func["transect_id"])
    print("  [diagnostic] Functional groin vs Deteriorated:")
    _diagnose_era_overlap("Functional groin", func["transect_id"],
                           "Deteriorated", det["transect_id"])

    # Functional minus Pre-install, transects in both; drop func's duplicate distance column first
    m1 = pre.merge(func.drop(columns=["dist_from_groin_m"]),
                    on="transect_id", suffixes=("_pre", "_func"))
    m1["diff"] = m1["slope_m_yr_func"] - m1["slope_m_yr_pre"]

    # Deteriorated − Functional: only transects present in BOTH
    m2 = func.merge(det.drop(columns=["dist_from_groin_m"]),
                     on="transect_id", suffixes=("_func", "_det"))
    m2["diff"] = m2["slope_m_yr_det"] - m2["slope_m_yr_func"]

    print(f"  Functional − Pre-install: {len(m1)} transects with a valid "
          f"rate in both eras")
    print(f"  Deteriorated − Functional: {len(m2)} transects with a valid "
          f"rate in both eras")
    _diagnose_alongshore_gaps(m1, "dist_from_groin_m",
                               "Pre-install", pre, "Functional groin", func)
    _diagnose_alongshore_gaps(m2, "dist_from_groin_m",
                               "Functional groin", func, "Deteriorated", det)

    fig, ax = plt.subplots(figsize=(14, 7))

    m1 = m1.sort_values("dist_from_groin_m")
    m2 = m2.sort_values("dist_from_groin_m")

    ax.plot(m1["dist_from_groin_m"] / 1000, m1["diff"],
             color=ERA_COLORS["Functional groin"], lw=1.1, alpha=0.55, zorder=4)
    ax.scatter(m1["dist_from_groin_m"] / 1000, m1["diff"],
                s=14, color=ERA_COLORS["Functional groin"], alpha=0.45,
                edgecolor="none",
                label=f"Functional groin − Pre-install baseline  "
                       f"(n={len(m1)}; did the groin change the rate?)")
    ax.plot(m2["dist_from_groin_m"] / 1000, m2["diff"],
             color=ERA_COLORS["Deteriorated"], lw=1.1, alpha=0.55, zorder=4)
    ax.scatter(m2["dist_from_groin_m"] / 1000, m2["diff"],
                s=14, color=ERA_COLORS["Deteriorated"], alpha=0.45,
                edgecolor="none",
                label=f"Deteriorated − Functional groin  "
                       f"(n={len(m2)}; positive = rate higher after "
                       f"deterioration)")

    ax.axhline(0, color="#999", lw=0.7, linestyle="-")
    ax.axvline(0, color="black", lw=1.0, alpha=0.7,
                label=f"Groin  (updrift = "
                       f"{'north' if UPDRIFT_DIRECTION == 'north' else 'south'})")
    if groin.get("geometry") is not None:
        sign = 1 if UPDRIFT_DIRECTION.lower().startswith("n") else -1
        g_a = sign * (groin["y_min"] - groin["reference_y"]) / 1000
        g_b = sign * (groin["y_max"] - groin["reference_y"]) / 1000
        ax.axvspan(min(g_a, g_b), max(g_a, g_b),
                    color="black", alpha=0.12, zorder=0,
                    label="Groin field  (first to last groin)")

    d1_km = domain_dist_km(transects, 1)
    apply_distance_xlim_and_autoscale_y(ax, x_max_km=PLOT_UPDRIFT_MAX_KM,
                                          x_min_km=d1_km)
    x_lo, x_hi = ax.get_xlim()
    set_domain_primary_axis(ax, transects)
    add_geographic_annotations(ax, transects)

    ax.set_ylabel("Δ LRR  (m/yr)", fontsize=FONT_AXIS_LABEL)
    ax.set_title(
        "Era-to-Era Difference in Shoreline Change Rate",
        fontsize=FONT_TITLE, fontweight="bold", pad=10,
    )
    ax.grid(True, alpha=0.25, linewidth=0.4)
    handles, labels = ax.get_legend_handles_labels()
    handles += annotation_legend_handles()
    ax.legend(handles=handles, loc="best", fontsize=FONT_LEGEND, framealpha=0.92, ncol=2)

    plt.tight_layout()
    plt.savefig(output_path, dpi=200, bbox_inches="tight")
    plt.close(fig)


# Print how many transects two era sets share before differencing
def _diagnose_era_overlap(name_a: str, ids_a, name_b: str, ids_b):
    set_a, set_b = set(ids_a), set(ids_b)
    overlap = set_a & set_b
    print(f"  {name_a}: {len(set_a)} unique transects with a valid rate")
    print(f"  {name_b}: {len(set_b)} unique transects with a valid rate")
    print(f"  In both: {len(overlap)} transects")
    smaller = min(len(set_a), len(set_b))
    if smaller > 0 and len(overlap) < smaller * 0.9:
        only_a = sorted(set_a - set_b)[:5]
        only_b = sorted(set_b - set_a)[:5]
        print(f"  ! overlap ({len(overlap)}) is well below the smaller "
              f"individual count ({smaller}) -- if {name_a} and {name_b} "
              f"are supposed to be nearly the same transects at "
              f"different times, this is worth investigating rather "
              f"than assumed correct. Examples only in {name_a}: "
              f"{only_a}; only in {name_b}: {only_b}")


# Print alongshore stretches with no matched data after a merge
def _diagnose_alongshore_gaps(merged_df: pd.DataFrame, dist_col: str,
                               name_a: str, full_a: pd.DataFrame,
                               name_b: str, full_b: pd.DataFrame,
                               gap_threshold_km: float = 1.5):
    if len(merged_df) < 2:
        return
    d = np.sort(merged_df[dist_col].values) / 1000.0
    gaps = np.diff(d)
    gap_idx = np.where(gaps > gap_threshold_km)[0]
    if len(gap_idx) == 0:
        return
    print(f"  [diagnostic] {len(gap_idx)} alongshore gap(s) > "
          f"{gap_threshold_km} km in the {name_a}/{name_b} match "
          f"(shows up as a straight-line bridge on the plot):")
    for i in gap_idx:
        lo, hi = float(d[i]), float(d[i + 1])
        # Strictly inside the gap: lo and hi exist on both sides by construction
        a_has = int(((full_a[dist_col] / 1000 > lo) &
                      (full_a[dist_col] / 1000 < hi)).sum())
        b_has = int(((full_b[dist_col] / 1000 > lo) &
                      (full_b[dist_col] / 1000 < hi)).sum())
        if a_has == 0 and b_has > 0:
            cause = f"{name_a} has no coverage there -- likely a real data gap in that era"
        elif b_has == 0 and a_has > 0:
            cause = f"{name_b} has no coverage there -- likely a real data gap in that era"
        elif a_has == 0 and b_has == 0:
            cause = "neither era has coverage there -- real gap in both"
        else:
            cause = (f"BOTH sides have data there ({name_a}={a_has}, "
                      f"{name_b}={b_has}) but it isn't matching by "
                      f"transect_id -- worth checking transect_id "
                      f"consistency between these two tables")
        print(f"    {lo:+.2f} to {hi:+.2f} km: {name_a} alone has "
              f"{a_has} transect(s) there, {name_b} alone has {b_has} "
              f"-- {cause}")


# Shared drawing for the single-comparison difference figures
def _plot_single_difference(dist_km, diff_values, n: int, color: str,
                             series_label: str, question_text: str,
                             title_prefix: str, transects: gpd.GeoDataFrame,
                             groin: dict, output_path: str,
                             curve_a_values=None, curve_a_label: str = None,
                             curve_a_color: str = None,
                             curve_b_values=None, curve_b_label: str = None,
                             curve_b_color: str = None):
    fig, ax = plt.subplots(figsize=(14, 7))

    # Both original curves and the shaded gap, drawn first so the delta stays on top
    if curve_a_values is not None and curve_b_values is not None:
        ax.fill_between(dist_km, curve_a_values, curve_b_values,
                          color=color, alpha=0.12, zorder=1)
        ax.plot(dist_km, curve_a_values, color=curve_a_color, lw=1.3,
                 alpha=0.75, zorder=2, label=curve_a_label)
        ax.plot(dist_km, curve_b_values, color=curve_b_color, lw=1.3,
                 alpha=0.75, zorder=2, label=curve_b_label)

    ax.plot(dist_km, diff_values, color=color, lw=0.7, alpha=0.35, zorder=4)
    ax.scatter(dist_km, diff_values, s=8, color=color, alpha=0.35,
                edgecolor="none", label=f"{series_label}  (n={n})")
    xs, ys = _lowess_overlay(np.asarray(dist_km), np.asarray(diff_values),
                               SMOOTHED_OVERLAY_LOWESS_FRAC)
    if xs is not None:
        ax.plot(xs, ys, color=color, lw=2.8, zorder=6,
                 label=f"{series_label}  (smoothed)")

    ax.axhline(0, color="#999", lw=0.7, linestyle="-")
    ax.axvline(0, color="black", lw=1.0, alpha=0.7,
                label=f"Groin  (updrift = "
                       f"{'north' if UPDRIFT_DIRECTION == 'north' else 'south'})")
    if groin.get("geometry") is not None:
        sign = 1 if UPDRIFT_DIRECTION.lower().startswith("n") else -1
        g_a = sign * (groin["y_min"] - groin["reference_y"]) / 1000
        g_b = sign * (groin["y_max"] - groin["reference_y"]) / 1000
        ax.axvspan(min(g_a, g_b), max(g_a, g_b), color="black", alpha=0.12,
                    zorder=0, label="Groin field  (first to last groin)")

    d1_km = domain_dist_km(transects, 1)
    apply_distance_xlim_and_autoscale_y(ax, x_max_km=PLOT_UPDRIFT_MAX_KM,
                                          x_min_km=d1_km)
    x_lo, x_hi = ax.get_xlim()
    set_domain_primary_axis(ax, transects)
    add_geographic_annotations(ax, transects)

    ax.set_ylabel("LRR (m/yr)  /  Δ LRR (m/yr)", fontsize=FONT_AXIS_LABEL)
    ax.set_title(
        title_prefix,
        fontsize=FONT_TITLE, fontweight="bold", pad=10,
    )
    ax.grid(True, alpha=0.25, linewidth=0.4)
    handles, labels = ax.get_legend_handles_labels()
    handles += annotation_legend_handles()
    ax.legend(handles=handles, loc="best", fontsize=FONT_LEGEND, framealpha=0.92, ncol=2)

    plt.tight_layout()
    plt.savefig(output_path, dpi=200, bbox_inches="tight")
    plt.close(fig)


# Figure: Functional minus Pre-install LRR
def plot_pre_to_functional_difference(era_lrrs: pd.DataFrame,
                                       preinstall_lrr: pd.DataFrame,
                                       transects: gpd.GeoDataFrame,
                                       groin: dict, output_path: str):
    print(f"\n[plot] pre-to-functional difference profile → {output_path}")
    pre = preinstall_lrr[(preinstall_lrr["n_obs"] >= 3) &
                          preinstall_lrr["slope_m_yr"].notna()][
        ["transect_id", "dist_from_groin_m", "slope_m_yr"]].copy()
    func = era_lrrs[era_lrrs["era"] == "Functional groin"][
        ["transect_id", "dist_from_groin_m", "slope_m_yr"]].copy()
    if len(pre) == 0 or len(func) == 0:
        print("  ! missing pre-install or functional era data -- skipping")
        return
    pre["transect_id"] = pre["transect_id"].astype(str)
    func["transect_id"] = func["transect_id"].astype(str)
    _diagnose_era_overlap("Pre-install", pre["transect_id"],
                           "Functional groin", func["transect_id"])
    m = pre.merge(func.drop(columns=["dist_from_groin_m"]),
                   on="transect_id", suffixes=("_pre", "_func"))
    m["diff"] = m["slope_m_yr_func"] - m["slope_m_yr_pre"]
    m = m.sort_values("dist_from_groin_m")
    print(f"  Functional − Pre-install: {len(m)} transects with a valid "
          f"rate in both eras")
    _diagnose_alongshore_gaps(m, "dist_from_groin_m",
                               "Pre-install", pre, "Functional groin", func)

    _plot_single_difference(
        m["dist_from_groin_m"] / 1000, m["diff"], len(m),
        DELTA_COLOR,
        "Functional groin − Pre-install baseline",
        "did the groin change the rate relative to background?",
        "Did the groin change the shoreline change rate?",
        transects, groin, output_path,
        curve_a_values=m["slope_m_yr_pre"],
        curve_a_label=f"Pre-install baseline  (n={len(m)})",
        curve_a_color=ERA_COLORS["Pre-install"],
        curve_b_values=m["slope_m_yr_func"],
        curve_b_label=f"Functional groin  (n={len(m)})",
        curve_b_color=ERA_COLORS["Functional groin"])


# Figure: Deteriorated minus Pre-install LRR
def plot_deteriorated_to_preinstall_difference(era_lrrs: pd.DataFrame,
                                                preinstall_lrr: pd.DataFrame,
                                                transects: gpd.GeoDataFrame,
                                                groin: dict, output_path: str):
    print(f"\n[plot] deteriorated-to-preinstall difference profile → {output_path}")
    pre = preinstall_lrr[(preinstall_lrr["n_obs"] >= 3) &
                          preinstall_lrr["slope_m_yr"].notna()][
        ["transect_id", "dist_from_groin_m", "slope_m_yr"]].copy()
    det = era_lrrs[era_lrrs["era"] == "Deteriorated"][
        ["transect_id", "dist_from_groin_m", "slope_m_yr"]].copy()
    if len(pre) == 0 or len(det) == 0:
        print("  ! missing pre-install or deteriorated era data -- skipping")
        return
    pre["transect_id"] = pre["transect_id"].astype(str)
    det["transect_id"] = det["transect_id"].astype(str)
    _diagnose_era_overlap("Pre-install", pre["transect_id"],
                           "Deteriorated", det["transect_id"])
    m = pre.merge(det.drop(columns=["dist_from_groin_m"]),
                   on="transect_id", suffixes=("_pre", "_det"))
    m["diff"] = m["slope_m_yr_det"] - m["slope_m_yr_pre"]
    m = m.sort_values("dist_from_groin_m")
    print(f"  Deteriorated − Pre-install: {len(m)} transects with a valid "
          f"rate in both eras")
    _diagnose_alongshore_gaps(m, "dist_from_groin_m",
                               "Pre-install", pre, "Deteriorated", det)

    _plot_single_difference(
        m["dist_from_groin_m"] / 1000, m["diff"], len(m),
        DELTA_COLOR,
        "Deteriorated − Pre-install baseline",
        "positive = the rate today is higher than it was before the "
        "groin was ever built",
        "How has the rate changed across the whole historical record?",
        transects, groin, output_path,
        curve_a_values=m["slope_m_yr_pre"],
        curve_a_label=f"Pre-install baseline  (n={len(m)})",
        curve_a_color=ERA_COLORS["Pre-install"],
        curve_b_values=m["slope_m_yr_det"],
        curve_b_label=f"Deteriorated  (n={len(m)})",
        curve_b_color=ERA_COLORS["Deteriorated"])


# Figure: Deteriorated minus Functional LRR
def plot_functional_to_deteriorated_difference(era_lrrs: pd.DataFrame,
                                                transects: gpd.GeoDataFrame,
                                                groin: dict, output_path: str):
    print(f"\n[plot] functional-to-deteriorated difference profile → {output_path}")
    func = era_lrrs[era_lrrs["era"] == "Functional groin"][
        ["transect_id", "dist_from_groin_m", "slope_m_yr"]].copy()
    det = era_lrrs[era_lrrs["era"] == "Deteriorated"][
        ["transect_id", "dist_from_groin_m", "slope_m_yr"]].copy()
    if len(func) == 0 or len(det) == 0:
        print("  ! missing functional or deteriorated era data -- skipping")
        return
    func["transect_id"] = func["transect_id"].astype(str)
    det["transect_id"] = det["transect_id"].astype(str)
    _diagnose_era_overlap("Functional groin", func["transect_id"],
                           "Deteriorated", det["transect_id"])
    m = func.merge(det.drop(columns=["dist_from_groin_m"]),
                    on="transect_id", suffixes=("_func", "_det"))
    m["diff"] = m["slope_m_yr_det"] - m["slope_m_yr_func"]
    m = m.sort_values("dist_from_groin_m")
    print(f"  Deteriorated − Functional: {len(m)} transects with a valid "
          f"rate in both eras")
    _diagnose_alongshore_gaps(m, "dist_from_groin_m",
                               "Functional groin", func, "Deteriorated", det)

    _plot_single_difference(
        m["dist_from_groin_m"] / 1000, m["diff"], len(m),
        DELTA_COLOR,
        "Deteriorated − Functional groin",
        "positive = rate was higher after deterioration than during "
        "the functional era",
        "Did the rate change again as the groin deteriorated?",
        transects, groin, output_path,
        curve_a_values=m["slope_m_yr_func"],
        curve_a_label=f"Functional groin  (n={len(m)})",
        curve_a_color=ERA_COLORS["Functional groin"],
        curve_b_values=m["slope_m_yr_det"],
        curve_b_label=f"Deteriorated  (n={len(m)})",
        curve_b_color=ERA_COLORS["Deteriorated"])


# Median value per distance band per decade, both sides (CSV only)
def compute_distance_band_series(decadal_df: pd.DataFrame,
                                  value_col: str) -> pd.DataFrame:
    edges = DISTANCE_BAND_EDGES_M
    df = decadal_df.copy()
    df["side"] = np.where(df["dist_from_groin_m"] >= 0, "updrift", "downdrift")
    df["abs_dist_m"] = df["dist_from_groin_m"].abs()

    records = []
    for lo, hi in zip(edges[:-1], edges[1:]):
        label = f"{lo/1000:.0f}\u2013{hi/1000:.0f} km"
        band = df[(df["abs_dist_m"] >= lo) & (df["abs_dist_m"] < hi)]
        for (side, cy), sgrp in band.groupby(["side", "decade_start"]):
            vals = sgrp[value_col].dropna()
            if len(vals) == 0:
                continue
            records.append({
                "decade_start": cy, "side": side, "band_label": label,
                "band_lo_m": lo, "band_hi_m": hi,
                "n": len(vals), "median": float(vals.median()),
            })
    return pd.DataFrame(records)


# Shoreline evolution GIF (groin area)

# Animated raw shoreline position through time around the groin
def create_groin_evolution_gif(chainage: pd.DataFrame,
                                transects: gpd.GeoDataFrame,
                                groin: dict,
                                output_dir: str,
                                domain_table: pd.DataFrame = None,
                                grid100m: gpd.GeoDataFrame = None,
                                window_domain_max: int = None,
                                output_suffix: str = ""):
    try:
        from PIL import Image
    except ImportError:
        print("\n[gif] Pillow not installed (pip install pillow) -- "
              "skipping shoreline evolution GIF")
        return None

    if GIF_TRANSECT_SOURCE not in ("coastsat", "hybrid", "grid100m"):
        raise ValueError(f"GIF_TRANSECT_SOURCE must be 'coastsat', "
                          f"'hybrid', or 'grid100m', got "
                          f"{GIF_TRANSECT_SOURCE!r}")

    # x_max_km: PLOT_UPDRIFT_MAX_KM, or the position of window_domain_max when given
    if window_domain_max is not None:
        x_max_km = domain_dist_km(transects, window_domain_max)
        if x_max_km is None:
            print(f"  ! domain {window_domain_max} has no transects -- "
                  f"falling back to PLOT_UPDRIFT_MAX_KM for the window")
            x_max_km = PLOT_UPDRIFT_MAX_KM
        window_desc = f"domain 1 to domain {window_domain_max} (~{x_max_km:+.1f} km)"
    else:
        x_max_km = PLOT_UPDRIFT_MAX_KM
        window_desc = f"domain 1 to +{PLOT_UPDRIFT_MAX_KM} km updrift"

    print(f"\n{'='*72}\nShoreline evolution GIF (groin area, {window_desc})")
    print(f"  Transect source: {GIF_TRANSECT_SOURCE}")

    frame_dir = os.path.join(output_dir, f"gif_frames_groin_area{output_suffix}")
    os.makedirs(frame_dir, exist_ok=True)

    effective_source = GIF_TRANSECT_SOURCE
    if GIF_TRANSECT_SOURCE in ("hybrid", "grid100m") and grid100m is None:
        print(f"  ! GIF_TRANSECT_SOURCE={GIF_TRANSECT_SOURCE!r} but no 100m "
              f"grid was loaded -- falling back to 'coastsat'")
        effective_source = "coastsat"

    # Restrict to domain 1 through x_max_km, in CoastSat's coordinates
    d1_km = domain_dist_km(transects, 1)
    x_min_km = d1_km if d1_km is not None else -PLOT_UPDRIFT_MAX_KM
    tx = transects[(transects["dist_from_groin_m"] / 1000 >= x_min_km) &
                   (transects["dist_from_groin_m"] / 1000 <= x_max_km)].copy()
    tx = tx.sort_values("dist_from_groin_m").reset_index(drop=True)
    print(f"  {len(tx)} CoastSat transects from domain 1 ({x_min_km:+.2f} km) "
          f"to {x_max_km:+.1f} km updrift of the groin")
    if len(tx) < 5:
        print("  ! too few transects in window -- skipping GIF")
        return None

    sub = chainage[chainage["transect_id"].astype(str).isin(
        tx["transect_id"].astype(str))].copy()
    sub["transect_id"] = sub["transect_id"].astype(str)
    sub["date"] = pd.to_datetime(sub["date"])
    if len(sub) == 0:
        print("  ! no chainage data in the groin-area window -- skipping GIF")
        return None

    if effective_source == "coastsat":
        tids_ordered = tx["transect_id"].astype(str).tolist()
        dist_km = tx["dist_from_groin_m"].values / 1000.0

    elif effective_source == "hybrid":
        # CoastSat's chainage is untouched; only each transect's x position changes
        print("  Repositioning CoastSat transects onto the 100 m grid's "
              "alongshore coordinates (measurements unchanged)...")
        tx["transect_id"] = tx["transect_id"].astype(str)
        pos_map = build_hybrid_position_map(tx, grid100m, groin)
        tx["dist_from_groin_m"] = tx["transect_id"].map(pos_map)
        tx = tx.dropna(subset=["dist_from_groin_m"])
        tx = tx.sort_values("dist_from_groin_m").reset_index(drop=True)
        tids_ordered = tx["transect_id"].tolist()
        dist_km = tx["dist_from_groin_m"].values / 1000.0
        print(f"  → repositioned {len(tx)} transects "
              f"({dist_km.min():+.2f} to {dist_km.max():+.2f} km)")

    else:   # "grid100m"
        print("  ! reconstructing shoreline points and re-measuring "
              "against the 100 m grid's own (parallel) transects -- "
              "this reintroduces a known geometric artifact, see "
              "GIF_TRANSECT_SOURCE config comment")
        sub = reconstruct_and_remeasure_on_grid100m(sub, tx, grid100m)
        if len(sub) == 0:
            print("  ! no chainage survived re-measurement -- skipping GIF")
            return None

        # Window and order the 100 m-grid transects the same way, by northing
        gx, gy = groin.get("reference_x"), groin["reference_y"]
        if gx is not None:
            d2 = (grid100m["shore_x"] - gx) ** 2 + (grid100m["shore_y"] - gy) ** 2
        else:
            d2 = (grid100m["shore_y"] - gy) ** 2
        groin_along_100m = float(grid100m.loc[d2.idxmin(), "alongshore_m"])
        sign = 1 if UPDRIFT_DIRECTION.lower().startswith("n") else -1
        grid100m = grid100m.copy()
        grid100m["dist_from_groin_m"] = sign * (grid100m["alongshore_m"] - groin_along_100m)

        x_min_km_100m = -PLOT_UPDRIFT_MAX_KM
        if domain_table is not None:
            d1_row = domain_table[domain_table["domain_id"] == 1]
            if len(d1_row) > 0:
                d1_y_min = float(d1_row["y_min"].iloc[0])
                grid_y_min = float(grid100m["shore_y"].min())
                d1_along_100m = d1_y_min - grid_y_min
                x_min_km_100m = sign * (d1_along_100m - groin_along_100m) / 1000

        grid_window = grid100m[
            (grid100m["dist_from_groin_m"] / 1000 >= x_min_km_100m) &
            (grid100m["dist_from_groin_m"] / 1000 <= x_max_km)
        ].sort_values("dist_from_groin_m").reset_index(drop=True)
        tids_ordered = grid_window["transect_id"].astype(str).tolist()
        dist_km = grid_window["dist_from_groin_m"].values / 1000.0
        tx = grid_window   # used downstream by set_domain_primary_axis / add_geographic_annotations
        print(f"  → {len(tids_ordered)} 100m-grid transects in window "
              f"({x_min_km_100m:+.2f} to {x_max_km:.0f} km)")

    # One row of raw chainage per frame (date or year)
    if GIF_FRAME_MODE == "year":
        frame_keys_all = sorted(sub["year"].unique())
        print(f"  {len(frame_keys_all)} unique years in the window "
              f"({frame_keys_all[0]} → {frame_keys_all[-1]}), "
              f"pooling all sources/dates within each year")
    elif GIF_FRAME_MODE == "date":
        frame_keys_all = sorted(sub["date"].unique())
        print(f"  {len(frame_keys_all)} unique observation dates in the "
              f"window ({frame_keys_all[0].date()} → "
              f"{frame_keys_all[-1].date()})")
    else:
        raise ValueError(f"GIF_FRAME_MODE must be 'year' or 'date', "
                          f"got {GIF_FRAME_MODE!r}")

    frames_data = {}
    n_transects_per_key = {}
    for k in frame_keys_all:
        day_data = sub[sub["year"] == k] if GIF_FRAME_MODE == "year" \
                   else sub[sub["date"] == k]
        agg = day_data.groupby("transect_id")["chainage_m"].median()
        row = np.full(len(tids_ordered), np.nan)
        for i, tid in enumerate(tids_ordered):
            if tid in agg.index:
                row[i] = agg.loc[tid]
        n_present = int(np.sum(~np.isnan(row)))
        n_transects_per_key[k] = n_present
        frame_year = k if GIF_FRAME_MODE == "year" else k.year
        required_fraction = (0.0 if frame_year < GIF_PRE_COASTSAT_CUTOFF_YEAR
                             else GIF_MIN_TRANSECT_FRACTION)
        if (n_present >= required_fraction * len(tids_ordered)
                and n_present >= GIF_MIN_TRANSECT_ABS):
            frames_data[k] = row

    valid_keys = sorted(frames_data.keys())
    unit = "years" if GIF_FRAME_MODE == "year" else "dates"
    n_pre_coastsat_kept = sum(1 for k in valid_keys
                                if (k if GIF_FRAME_MODE == "year" else k.year)
                                < GIF_PRE_COASTSAT_CUTOFF_YEAR)
    print(f"  {len(valid_keys)}/{len(frame_keys_all)} {unit} have enough "
          f"transect coverage to render a frame: before "
          f"{GIF_PRE_COASTSAT_CUTOFF_YEAR} every {unit[:-1]} with >= "
          f"{GIF_MIN_TRANSECT_ABS} transects is kept regardless of "
          f"coverage fraction ({n_pre_coastsat_kept} such {unit} kept); "
          f"{GIF_PRE_COASTSAT_CUTOFF_YEAR} onward requires "
          f">= {GIF_MIN_TRANSECT_FRACTION:.0%} of "
          f"{len(tids_ordered)} transects")
    if len(frame_keys_all) > 0:
        counts = np.array(list(n_transects_per_key.values()))
        print(f"  Coverage per {unit[:-1]}: median {int(np.median(counts))}, "
              f"a few dates would additionally qualify at lower "
              f"thresholds -- e.g. >=1: {int((counts >= 1).sum())}, "
              f">=5: {int((counts >= 5).sum())}, "
              f">=20: {int((counts >= 20).sum())} "
              f"(adjust GIF_MIN_TRANSECT_FRACTION / GIF_MIN_TRANSECT_ABS "
              f"if you want more or fewer frames)")
    if len(valid_keys) < 3:
        print(f"  ! too few valid {unit} -- skipping GIF "
              "(try lowering GIF_MIN_TRANSECT_FRACTION / GIF_MIN_TRANSECT_ABS)")
        return None

    # Reference baseline: the GIF_REFERENCE_YEAR shoreline, built like any frame (README)
    ref_data = sub[sub["date"].dt.year == GIF_REFERENCE_YEAR]
    ref_row = np.full(len(tids_ordered), np.nan)
    if len(ref_data) > 0:
        ref_agg = ref_data.groupby("transect_id")["chainage_m"].median()
        for i, tid in enumerate(tids_ordered):
            if tid in ref_agg.index:
                ref_row[i] = ref_agg.loc[tid]
    n_ref = int(np.sum(~np.isnan(ref_row)))
    print(f"  {GIF_REFERENCE_YEAR} reference baseline: {n_ref}/{len(tx)} "
          f"transects in window have a {GIF_REFERENCE_YEAR} observation"
          + ("" if n_ref > 0 else f"  ! no {GIF_REFERENCE_YEAR} data in "
                                    f"this window -- reference line will "
                                    f"be blank"))

    # Diagnostic: observations per transect, to explain gaps that recur across frames
    coverage_counts = sub.groupby("transect_id").size()
    coverage_df = pd.DataFrame({
        "transect_id": tids_ordered,
        "dist_km": dist_km,
        "n_obs": [coverage_counts.get(tid, 0) for tid in tids_ordered],
    })
    sparsest = coverage_df.nsmallest(10, "n_obs")
    print(f"  Per-transect total observation count in this window: "
          f"median {int(coverage_df['n_obs'].median())}; "
          f"sparsest 10 transects (total obs across the whole record):")
    for _, row in sparsest.iterrows():
        print(f"    {row['transect_id']}  ({row['dist_km']:+.2f} km)  "
              f"n_obs={row['n_obs']}")

    # Shared y-axis: true min and max over every frame and the reference row
    all_vals = np.concatenate(
        [v[~np.isnan(v)] for v in frames_data.values()] + [ref_row[~np.isnan(ref_row)]])
    y_lo, y_hi = float(np.nanmin(all_vals)), float(np.nanmax(all_vals))
    pad = (y_hi - y_lo) * 0.08 if y_hi > y_lo else 1.0
    ylim = (y_lo - pad, y_hi + pad)

    def era_for_year(yr):
        for name, lo, hi in ERAS:
            if lo <= yr <= hi:
                return name
        return None

    x_lo_plot, x_hi_plot = x_min_km, x_max_km

    # Groin field footprint, relative to the north-end origin
    groin_span = None
    if groin.get("geometry") is not None:
        sign = 1 if UPDRIFT_DIRECTION.lower().startswith("n") else -1
        g_a = sign * (groin["y_min"] - groin["reference_y"]) / 1000
        g_b = sign * (groin["y_max"] - groin["reference_y"]) / 1000
        groin_span = (min(g_a, g_b), max(g_a, g_b))

    print(f"  Rendering {len(valid_keys)} frames...")
    frame_paths = []
    for fi, k in enumerate(valid_keys):
        if (fi + 1) % 50 == 0 or fi == len(valid_keys) - 1:
            print(f"    {fi+1}/{len(valid_keys)}", end="\r")

        if GIF_FRAME_MODE == "year":
            frame_year, label_str = int(k), str(int(k))
        else:
            frame_year, label_str = k.year, str(k.date())

        fig, ax = plt.subplots(figsize=(12, 5.5))

        # Reference baseline on every frame after the reference year
        if n_ref > 0 and frame_year > GIF_REFERENCE_YEAR:
            ax.plot(dist_km, ref_row, color="#222", lw=0.9, linestyle="--",
                     alpha=0.6, zorder=3,
                     label=f"{GIF_REFERENCE_YEAR} reference baseline")

        # Fading trail of preceding frames
        for tr in range(1, GIF_TRAIL_FRAMES + 1):
            j = fi - tr
            if j < 0:
                continue
            alpha = 0.30 * (1 - tr / (GIF_TRAIL_FRAMES + 1))
            ax.plot(dist_km, frames_data[valid_keys[j]], color="#999",
                     lw=1.0, alpha=alpha, zorder=2)

        row = frames_data[k]
        era_name = era_for_year(frame_year)
        color = ERA_COLORS.get(era_name, "#1565C0")

        # Water/island shading from an interpolated row; the line itself keeps the raw gaps
        row_series = pd.Series(row)
        if row_series.notna().sum() >= 2:
            row_filled = row_series.interpolate(
                limit_direction="both").to_numpy()
        else:
            row_filled = row
        ax.fill_between(dist_km, row_filled, ylim[1], color=GIF_WATER_COLOR,
                          alpha=0.55, zorder=1, linewidth=0)
        ax.fill_between(dist_km, ylim[0], row_filled, color=GIF_LAND_COLOR,
                          alpha=0.65, zorder=1, linewidth=0)

        ax.plot(dist_km, row, color=color, lw=1.0, marker="o", ms=1.3,
                 zorder=5,
                 label=f"{label_str}" + (f"  ({era_name})" if era_name else ""))

        ax.axvline(0, color="black", lw=1.0, alpha=0.7,
                    label=f"Groin  (updrift = "
                           f"{'north' if UPDRIFT_DIRECTION == 'north' else 'south'})")
        if groin_span is not None:
            ax.axvspan(groin_span[0], groin_span[1], color="black",
                        alpha=0.12, zorder=0,
                        label="Groin field  (first to last groin)")

        ax.set_xlim(x_lo_plot, x_hi_plot)
        ax.set_ylim(ylim[1], ylim[0])   # inverted: island (small chainage) on top, ocean (large chainage) on bottom
        set_domain_primary_axis(ax, tx)
        add_geographic_annotations(ax, tx, gif_mode=True)

        ax.set_ylabel("Shoreline position\n(chainage, m)", fontsize=FONT_AXIS_LABEL)
        ax.set_title(f"Groin-area shoreline position  —  {label_str}",
                       fontsize=FONT_TITLE, fontweight="bold")
        handles, labels = ax.get_legend_handles_labels()
        handles += [Patch(facecolor=GIF_WATER_COLOR, alpha=0.55, label="Ocean"),
                     Patch(facecolor=GIF_LAND_COLOR, alpha=0.65, label="Island")]
        handles += annotation_legend_handles()
        if window_domain_max is not None:
            ax.legend(handles=handles, loc="upper center",
                       fontsize=FONT_LEGEND_GIF_ZOOMED, ncol=2)
        else:
            ax.legend(handles=handles, loc="upper right",
                       fontsize=FONT_LEGEND_GIF, ncol=2)
        ax.grid(True, alpha=0.25, linewidth=0.4)

        plt.tight_layout()
        frame_path = os.path.join(frame_dir, f"frame_{fi:04d}_{label_str}.png")
        fig.savefig(frame_path, dpi=GIF_DPI, bbox_inches="tight")
        plt.close(fig)
        frame_paths.append(frame_path)
    print()

    gif_path = os.path.join(
        output_dir, f"groin_analysis_shoreline_evolution{output_suffix}.gif")
    print(f"  Assembling {len(frame_paths)}-frame GIF...")
    frames = [Image.open(fp).convert("RGBA") for fp in frame_paths]
    frames[0].save(
        gif_path, save_all=True, append_images=frames[1:],
        duration=int(GIF_FRAME_DURATION_S * 500),  # *500 corrects a Pillow duration bug
        loop=0,
    )
    print(f"  → saved {gif_path}")
    return gif_path


# Metadata

# Plain-text run metadata
def write_metadata(output_path: str, stats: dict):
    with open(output_path, "w", encoding="utf-8") as f:
        f.write("HAT_groin_shoreline_analysis_v2.py — run metadata\n")
        f.write("=" * 60 + "\n\n")
        f.write(f"Run at: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n\n")
        f.write("--- CONFIG ---\n")
        f.write(f"Pre-install year cutoff:  < {PRE_INSTALLATION_YEAR_CUTOFF} "
                f"(all shorelines before this year used)\n")
        f.write(f"Groin geojson: {GROIN_GEOJSON_PATH}\n")
        f.write(f"Updrift direction: {UPDRIFT_DIRECTION}\n")
        f.write(f"Signal anomaly threshold: "
                f"{SIGNAL_ANOMALY_THRESHOLD_M_YR} m/yr\n")
        f.write(f"Signal max search distance: "
                f"{SIGNAL_MAX_SEARCH_DISTANCE_M/1000:.1f} km\n")
        f.write(f"Signal extent bin width: "
                f"{SIGNAL_EXTENT_BIN_WIDTH_M:.0f} m, "
                f"max gap {SIGNAL_EXTENT_MAX_GAP_BINS} bin(s), contiguous from groin\n")
        f.write(f"Baseline method: connect-the-dots linear interpolation "
                f"across per-transect pre-install rates (no binning/smoothing)\n")
        f.write(f"Decadal bins: {DECADE_LENGTH_YEARS} yr, "
                f"start {DECADE_START_YEAR}, "
                f"min obs {DECADE_MIN_OBSERVATIONS}\n")
        f.write(f"Breakpoint search: {BREAKPOINT_SEARCH_START}–"
                f"{BREAKPOINT_SEARCH_END}, "
                f"min pts/segment {BREAKPOINT_MIN_POINTS_PER_SEGMENT}\n")
        f.write(f"Shoreline evolution GIF window: domain 1 to "
                f"+{PLOT_UPDRIFT_MAX_KM} km updrift, "
                f"{GIF_FRAME_DURATION_S}s/frame, mode={GIF_FRAME_MODE}\n")
        f.write(f"CoastSat transect layer: {COASTSAT_TRANSECT_GEOM}\n")
        f.write(f"Transect length: {TRANSECT_INTERSECTION_LENGTH_M} m, "
                f"chainage max abs {CHAINAGE_MAX_ABS_M} m\n")
        f.write(f"CASCADE domain reference: {DOMAINS_JSON_PATH}\n")
        f.write(f"Plot x-axis clip: ±{PLOT_X_MAX_KM} km\n")
        f.write(f"Min obs per transect: {MIN_OBSERVATIONS_PER_TRANSECT}\n\n")
        f.write("--- STATS ---\n")
        for k, v in stats.items():
            f.write(f"{k}: {v}\n")


# Run: load, measure, fit, plot, animate, write
def main():
    print("=" * 72)
    print("HAT GROIN SHORELINE ANALYSIS — Hatteras Island")
    print("=" * 72)

    os.makedirs(OUTPUT_DIR, exist_ok=True)

    # 1. Study area filter
    filter_gdf = load_study_area_filter()

    # 2. CASCADE domain reference, loaded first so the groin can be checked against it
    domain_table = load_domain_reference()

    # 3. Load groin geometry (sets the origin for distance-from-groin)
    groin = load_groin_geometry(domain_table)

    # 4. Load shoreline sources
    wet_dry_gdf = load_shoreline_lines(WET_DRY_PATH, WET_DRY_DATE_COL,
                                        "wet_dry", filter_gdf)
    nc_state_gdf = load_shoreline_lines(NC_STATE_PATH, NC_STATE_DATE_COL,
                                         "nc_state", filter_gdf)

    # 5. CoastSat's transect network, the backbone for every source
    transects = load_coastsat_transects(filter_gdf, domain_table)
    transects = assign_distance_from_groin(transects, groin)
    print(f"  → distance-from-groin range: "
          f"{transects['dist_from_groin_m'].min()/1000:+.2f} to "
          f"{transects['dist_from_groin_m'].max()/1000:+.2f} km")
    n_updrift = int((transects["dist_from_groin_m"] > 0).sum())
    n_downdrift = int((transects["dist_from_groin_m"] < 0).sum())
    print(f"  → {n_updrift} updrift transects, {n_downdrift} downdrift "
          f"transects")

    # 5b. Drop transects south of domain 1, by alongshore position (README)
    domain1_row = domain_table[domain_table["domain_id"] == 1]
    if len(domain1_row) > 0:
        domain1_y_min = float(domain1_row["y_min"].iloc[0])
        slope, intercept = np.polyfit(transects["origin_y"].values,
                                        transects["dist_from_groin_m"].values, 1)
        dist_cutoff_m = slope * domain1_y_min + intercept
        n_before = len(transects)
        transects = transects[transects["dist_from_groin_m"] >= dist_cutoff_m].reset_index(drop=True)
        n_dropped = n_before - len(transects)
        print(f"  Domain 1's true northing boundary ({domain1_y_min:,.0f}) "
              f"maps to {dist_cutoff_m/1000:+.2f} km alongshore (via a "
              f"linear fit across all transects, robust to Cape Point's "
              f"local northing-vs-alongshore disagreement)")
        print(f"  Excluded {n_dropped}/{n_before} transect(s) south of "
              f"that alongshore position -- too close to Cape Point to "
              f"be reliable for this analysis")
        if len(transects) > 0:
            print(f"  → distance-from-groin range after exclusion: "
                  f"{transects['dist_from_groin_m'].min()/1000:+.2f} to "
                  f"{transects['dist_from_groin_m'].max()/1000:+.2f} km")
    else:
        print("  ! domain 1 not found in domain reference -- skipping "
              "Cape Point exclusion")

    # 6. CoastSat chainage, read directly by transect ID
    coastsat_chainage = load_coastsat_chainage(transects, COASTSAT_ROOT_DIR)

    # 6. Wet-dry and NC state chainage, by intersection with the same transects
    wet_dry_chainage = extract_chainage_by_intersection(
        wet_dry_gdf, transects, "wet_dry")
    nc_state_chainage = extract_chainage_by_intersection(
        nc_state_gdf, transects, "nc_state")

    # Sanity check: wet-dry coverage on each side of the groin
    wd_check = wet_dry_chainage.merge(
        transects[["transect_id", "dist_from_groin_m"]], on="transect_id", how="left")
    print(f"\n  [sanity check] wet-dry observations: "
          f"{int((wd_check['dist_from_groin_m'] > 0).sum())} updrift, "
          f"{int((wd_check['dist_from_groin_m'] < 0).sum())} downdrift "
          f"of the groin")

    # Spatial reference of every file loaded so far
    print_crs_summary()

    # 7. Merge into unified per-(transect, date) table
    chainage = build_unified_chainage_table(
        coastsat_chainage, wet_dry_chainage, nc_state_chainage, transects)
    chainage.to_csv(
        os.path.join(OUTPUT_DIR, "groin_analysis_chainage_all.csv"),
        index=False)
    print(f"→ saved chainage_all.csv ({len(chainage)} rows)")

    # 8. Pre-installation LRR (per transect where computable)
    preinstall_lrr = compute_preinstall_lrr(chainage, transects)
    preinstall_lrr.to_csv(
        os.path.join(OUTPUT_DIR, "groin_analysis_preinstall_lrr.csv"),
        index=False)

    # 9. Regional baseline for the anomaly: per-transect pre-install values, interpolated
    print(f"\n{'='*72}\nRegional baseline (connect-the-dots, not binned/smoothed)")
    baseline_fn = interp_preinstall_baseline(preinstall_lrr)

    # 10. Full-period LRR + piecewise breakpoint
    fullperiod_lrr = compute_fullperiod_lrr(chainage, transects)
    fullperiod_lrr.to_csv(
        os.path.join(OUTPUT_DIR, "groin_analysis_full_period_lrr.csv"),
        index=False)

    # 10b. Post-install-only LRR, 1970 on
    postinstall_lrr = compute_post_install_lrr(chainage, transects)
    postinstall_lrr.to_csv(
        os.path.join(OUTPUT_DIR, "groin_analysis_post_install_lrr.csv"),
        index=False)

    # 11. Decadal LRR (fixed 10-yr non-overlapping bins, starting 1960)
    decadal_lrr = compute_decadal_lrr(chainage)
    decadal_lrr.to_csv(
        os.path.join(OUTPUT_DIR, "groin_analysis_decadal_lrr.csv"),
        index=False)

    # 11b. Per-era LRR (one rate per transect per era)
    era_lrrs = compute_era_lrrs(chainage, transects,
                                  extra_periods=[PRE_NOURISHMENT_PERIOD])
    era_lrrs.to_csv(
        os.path.join(OUTPUT_DIR, "groin_analysis_era_lrrs.csv"),
        index=False)

    # 12. Decadal anomaly, with distance from the groin
    decadal_anomaly = compute_decadal_anomaly(decadal_lrr, transects,
                                                baseline_fn)
    decadal_anomaly.to_csv(
        os.path.join(OUTPUT_DIR, "groin_analysis_decadal_anomaly.csv"),
        index=False)

    # 13. Signal extent over time
    signal_extent = compute_signal_extent_over_time(decadal_anomaly)
    signal_extent.to_csv(
        os.path.join(OUTPUT_DIR, "groin_analysis_signal_extent.csv"),
        index=False)

    # 13b. Distance-band time series (CSV only)
    lrr_band_series = compute_distance_band_series(decadal_lrr.merge(
        transects[["transect_id", "dist_from_groin_m"]]
            .assign(transect_id=lambda d: d["transect_id"].astype(str)),
        on="transect_id", how="inner"), "slope_m_yr")
    lrr_band_series.to_csv(
        os.path.join(OUTPUT_DIR, "groin_analysis_lrr_distance_band_series.csv"),
        index=False)
    anomaly_band_series = compute_distance_band_series(
        decadal_anomaly, "anomaly_m_yr")
    anomaly_band_series.to_csv(
        os.path.join(OUTPUT_DIR, "groin_analysis_anomaly_distance_band_series.csv"),
        index=False)

    # 14. Plots (the dropped plots' CSVs are still written)
    print("\n" + "=" * 72)
    print("GENERATING PLOTS")
    print("=" * 72)
    plot_alongshore_lrr_profile(
        preinstall_lrr, fullperiod_lrr, postinstall_lrr, transects, groin,
        os.path.join(OUTPUT_DIR, "groin_analysis_alongshore_profile.png"))
    plot_era_lrr_profile(
        era_lrrs, preinstall_lrr, transects, groin,
        os.path.join(OUTPUT_DIR, "groin_analysis_era_lrr_profile.png"))
    plot_era_difference_profile(
        era_lrrs, preinstall_lrr, transects, groin,
        os.path.join(OUTPUT_DIR, "groin_analysis_era_difference_profile.png"))
    plot_pre_to_functional_difference(
        era_lrrs, preinstall_lrr, transects, groin,
        os.path.join(OUTPUT_DIR, "groin_analysis_diff_pre_to_functional.png"))
    plot_functional_to_deteriorated_difference(
        era_lrrs, transects, groin,
        os.path.join(OUTPUT_DIR, "groin_analysis_diff_functional_to_deteriorated.png"))
    plot_deteriorated_to_preinstall_difference(
        era_lrrs, preinstall_lrr, transects, groin,
        os.path.join(OUTPUT_DIR, "groin_analysis_diff_deteriorated_to_preinstall.png"))

    # Decade-increment LRR profiles, one per increment
    for increment_years in DECADE_PLOT_INCREMENTS_YEARS:
        decade_plot_periods = _build_decade_periods(
            DECADE_PLOT_START_YEAR, increment_years, DECADE_PLOT_END_YEAR)
        decade_period_lrrs = compute_era_lrrs(chainage, transects,
                                               extra_periods=decade_plot_periods)
        decade_period_lrrs.to_csv(
            os.path.join(OUTPUT_DIR,
                          f"groin_analysis_decade_increment_lrrs_{increment_years}yr.csv"),
            index=False)
        plot_decade_lrr_profile(
            decade_period_lrrs, transects, groin, decade_plot_periods,
            increment_years,
            os.path.join(OUTPUT_DIR,
                          f"groin_analysis_decade_increment_profile_{increment_years}yr.png"))

    # Zone panels: one per era, plus the early-functional windows
    early_functional_periods = [
        (f"Functional groin (first {n}yr)", GROIN_INSTALLATION_YEAR,
         GROIN_INSTALLATION_YEAR + n - 1)
        for n in ZONE_PANEL_FUNCTIONAL_EARLY_WINDOWS_YEARS
    ]
    early_functional_lrrs = compute_era_lrrs(chainage, transects,
                                              extra_periods=early_functional_periods)
    extra_functional_panels = [
        (name, early_functional_lrrs[early_functional_lrrs["era"] == name])
        for name, _, _ in early_functional_periods
    ]
    plot_zone_panels(
        preinstall_lrr, era_lrrs, transects, groin,
        os.path.join(OUTPUT_DIR, "groin_analysis_zone_panels.png"),
        extra_functional_panels=extra_functional_panels)

    # 14b. Shoreline evolution GIF
    grid100m = (load_100m_grid_transects(domain_table)
                if GIF_TRANSECT_SOURCE in ("hybrid", "grid100m") else None)
    create_groin_evolution_gif(chainage, transects, groin, OUTPUT_DIR,
                                domain_table=domain_table, grid100m=grid100m)
    create_groin_evolution_gif(chainage, transects, groin, OUTPUT_DIR,
                                domain_table=domain_table, grid100m=grid100m,
                                window_domain_max=GIF_ZOOM_WINDOW_DOMAIN_MAX,
                                output_suffix="_zoomed")

    # 15. Metadata
    stats = {
        "groin_reference_northing (northernmost groin)": groin["reference_y"],
        "groin_field_extent_min":  groin["y_min"],
        "groin_field_extent_max":  groin["y_max"],
        "n_groin_features":        groin["n_features"],
        "n_coastsat_transects":    len(transects),
        "n_transects_analyzed":    chainage["transect_id"].nunique(),
        "n_total_observations":    len(chainage),
        "n_wetdry_obs_updrift":    int((wd_check["dist_from_groin_m"] > 0).sum()),
        "n_wetdry_obs_downdrift":  int((wd_check["dist_from_groin_m"] < 0).sum()),
        "n_preinstall_transects":  int((preinstall_lrr["n_obs"] >= 3).sum()),
        "preinstall_years_used":   preinstall_lrr.attrs.get("preinstall_years", []),
        "n_piecewise_bp_found":    int(fullperiod_lrr["breakpoint_year"].notna().sum()),
        "median_pre_install_LRR_m_yr":
            float(preinstall_lrr["slope_m_yr"].median()),
        "median_full_period_LRR_m_yr":
            float(fullperiod_lrr["slope_full"].median()),
        "median_updrift_extent_km":
            float(signal_extent["updrift_extent_m"].median() / 1000),
        "median_downdrift_extent_km":
            float(signal_extent["downdrift_extent_m"].median() / 1000),
        "max_updrift_peak_anomaly_m_yr":
            float(signal_extent["updrift_peak_anomaly"].max()),
        "max_downdrift_peak_anomaly_m_yr":
            float(signal_extent["downdrift_peak_anomaly"].min()),
    }
    write_metadata(
        os.path.join(OUTPUT_DIR, "groin_analysis_metadata.txt"), stats)

    print("\n" + "=" * 72)
    print("DONE.")
    print("=" * 72)


if __name__ == "__main__":
    main()
