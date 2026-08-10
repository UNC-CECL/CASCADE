#!/usr/bin/env python3
"""
HATTERAS ISLAND: Shoreline Change Analysis from CASCADE NPZ Output
==================================================================
Loads pre-saved CASCADE NPZ comparison and produces:
  1. Yearly relative shoreline change + BN bar panel → GIF
  2. Yearly absolute shoreline position + BN bar panel → GIF
  3. Publication-quality rate profile vs CoastSat LOESS

Annotation system, color palette, LOESS CoastSat pipeline, and geographic
annotation data all match HAT_hindcast_1984_2024_old version.py exactly.

Usage:
  1. Set NPZ_PATHS_BY_LABEL  (Section 2) — one entry per run to compare.
  2. Set START_YEAR / END_YEAR            — must match the loaded run period.
  3. Set SOURCE_SINK_PRESET and Hs labels as needed for plot titles.
  4. Run from PyCharm or command line.
"""

import os
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from matplotlib.transforms import blended_transform_factory
from statsmodels.nonparametric.smoothers_lowess import lowess


# =============================================================================
# SECTION 1: DOMAIN CONFIGURATION  (matches HAT_hindcast_1984_2024_old version.py)
# =============================================================================

NUM_REAL_DOMAINS   = 90
NUM_BUFFER_DOMAINS = 15

FIRST_FILE_NUMBER = 1      # GIS domain IDs 1 .. 90
LAST_FILE_NUMBER  = FIRST_FILE_NUMBER + NUM_REAL_DOMAINS - 1   # = 90

TOTAL_DOMAINS    = NUM_BUFFER_DOMAINS + NUM_REAL_DOMAINS + NUM_BUFFER_DOMAINS  # 120
START_REAL_INDEX = NUM_BUFFER_DOMAINS           # = 15  (first real domain in padded array)
END_REAL_INDEX   = START_REAL_INDEX + NUM_REAL_DOMAINS  # = 105

DOMAIN_SPACING_M = 500
DOMAIN_LENGTH_M  = 500

print("=" * 80)
print("HATTERAS ISLAND   –   NPZ ANALYSIS SCRIPT")
print("=" * 80)
print(f"Real Domains:   {NUM_REAL_DOMAINS}  (GIS {FIRST_FILE_NUMBER}–{LAST_FILE_NUMBER})")
print(f"Buffers:        {NUM_BUFFER_DOMAINS} each side  |  TOTAL: {TOTAL_DOMAINS}")
print(f"Padded range:   [{START_REAL_INDEX}..{END_REAL_INDEX - 1}]")
print("=" * 80 + "\n")


def _gis_to_pad(gis_id):
    """1-based GIS domain ID → CASCADE padded array index."""
    return START_REAL_INDEX + (gis_id - FIRST_FILE_NUMBER)


def _pad_to_gis(pad_idx):
    """Padded index → GIS domain 1–90, or None if outside real range."""
    if START_REAL_INDEX <= pad_idx < END_REAL_INDEX:
        return FIRST_FILE_NUMBER + (pad_idx - START_REAL_INDEX)
    return None


# =============================================================================
# SECTION 2: FILE PATHS
# =============================================================================

PROJECT_BASE_DIR  = r"C:\Users\hanna\PycharmProjects\CASCADE"
OUTPUT_BASE_DIR   = os.path.join(PROJECT_BASE_DIR, "scripts", "figure_making", "position_plot", "comparison")
COASTSAT_BASE_DIR = os.path.join(
    PROJECT_BASE_DIR, "scripts", "input_prep", "CoastSat"
)

# ---------------------------------------------------------------------------
# NPZ paths — one entry per run you want to compare.
# Keys are used as legend labels in every plot.
# ---------------------------------------------------------------------------
NPZ_PATHS_BY_LABEL = {
    "Hindcast 1984–2004": (
        r"C:\Users\hanna\PycharmProjects\CASCADE\output\raw_runs\HAT_2004_2024_base_newbufferv3\HAT_2004_2024_base_newbufferv3.npz"
    ),
    # "Hindcast 2004–2024": (
    #     r"...\HAT_2004_2024_base_newbufferv3\cascade.npz"
    # ),
}

# CoastSat transect-level LRR datasets — matches COASTSAT_DATASETS in run script
COASTSAT_DATASETS = [
    dict(
        label           = "CoastSat LRR (1984–2004)",
        period_start    = 1984,
        csv             = os.path.join(COASTSAT_BASE_DIR, "1984_2004", "transect_lrr_full.csv"),
        domain_col      = "domain_number",
        rate_col        = "lrr_m_yr",
        transect_id_col = "transect_id",
    ),
    dict(
        label           = "CoastSat LRR (2004–2024)",
        period_start    = 2004,
        csv             = os.path.join(
            COASTSAT_BASE_DIR, "2004_2024", "transect_lrr_full.csv"
        ),
        domain_col      = "domain_number",
        rate_col        = "lrr_m_yr",
        transect_id_col = "transect_id",
    ),
]

os.makedirs(OUTPUT_BASE_DIR, exist_ok=True)


# =============================================================================
# SECTION 3: PLOT PARAMETERS
# =============================================================================

# ← Flip to 1984 or 2004 to match the period(s) in NPZ_PATHS_BY_LABEL
START_YEAR = 2004
END_YEAR   = 2024
RUN_YEARS  = END_YEAR - START_YEAR

SEA_LEVEL_RISE_RATE = 0.006   # m/yr — used in plot titles; 1984=0.004, 2004=0.006
FLIP_SIGN_MODEL     = True    # CASCADE x_s_TS increases with retreat → flip

SOURCE_SINK_PRESET = "base"   # used in plot title only (no calculation here)
Hs_LABEL           = 2.5      # used in plot title / legend

# LOESS windows to overlay (domain count) — matches run script
LOESS_WINDOW_DOMAINS = [7, 10]
LOESS_WINDOW_STYLES  = [
    (1.8, "-", 1.00),   # 7-domain: solid
    (2.0, "-", 1.00),   # 10-domain: solid, primary reference
]
RESIDUALS_LOESS_WINDOW      = 10
LOESS_SKIP_SOUTHERN_DOMAINS = 10   # raw means for domains 1–10; LOESS from 11+

PLOT_RAW_LRR          = True
PLOT_REFERENCE_PERIOD = False
RAW_LRR_SCATTER_SIZE  = 6
RAW_LRR_SCATTER_ALPHA = 0.60

MAKE_YEARLY_BN_GIF    = True
GIF_DURATION_SECONDS  = 4
DOMAIN_TICK_STEP      = 5


# =============================================================================
# SECTION 4: BEACH NOURISHMENT DISPLAY ARRAYS
# =============================================================================
# Matches Section 6b of HAT_hindcast_1984_2024_old version.py exactly.
# Used only for BN bar panels in yearly plots; actual nourishment was applied
# during the hindcast run, not here.

_CY_TO_M3 = 0.764555

HAT_BN_YEARS = [2014, 2022]

HAT_BN_VOLUME_BY_DOMAIN = {
    # Buxton: GIS 6–15, 1,200,000 cy in 2022
     6: [0, round(1_200_000 / 10 * _CY_TO_M3, 1)],
     7: [0, round(1_200_000 / 10 * _CY_TO_M3, 1)],
     8: [0, round(1_200_000 / 10 * _CY_TO_M3, 1)],
     9: [0, round(1_200_000 / 10 * _CY_TO_M3, 1)],
    10: [0, round(1_200_000 / 10 * _CY_TO_M3, 1)],
    11: [0, round(1_200_000 / 10 * _CY_TO_M3, 1)],
    12: [0, round(1_200_000 / 10 * _CY_TO_M3, 1)],
    13: [0, round(1_200_000 / 10 * _CY_TO_M3, 1)],
    14: [0, round(1_200_000 / 10 * _CY_TO_M3, 1)],
    15: [0, round(1_200_000 / 10 * _CY_TO_M3, 1)],
    # Avon: GIS 23–26, 2,200,000 cy in 2022
    23: [0, round(2_200_000 / 4 * _CY_TO_M3, 1)],
    24: [0, round(2_200_000 / 4 * _CY_TO_M3, 1)],
    25: [0, round(2_200_000 / 4 * _CY_TO_M3, 1)],
    26: [0, round(2_200_000 / 4 * _CY_TO_M3, 1)],
    # Rodanthe: GIS 85–88, 1,620,000 cy in 2014
    85: [round(1_620_000 / 4 * _CY_TO_M3, 1), 0],
    86: [round(1_620_000 / 4 * _CY_TO_M3, 1), 0],
    87: [round(1_620_000 / 4 * _CY_TO_M3, 1), 0],
    88: [round(1_620_000 / 4 * _CY_TO_M3, 1), 0],
}


# =============================================================================
# SECTION 5: COLOUR PALETTE  (matches run script exactly)
# =============================================================================

# Model line
ANN_MODEL_COLOR = "#FF8C00"   # warm orange

# Community / geographic annotation colours
ANN_C_TOWN_SPAN    = "#90AFC5"
ANN_C_WIMBLE       = "#E0A800"
ANN_C_VILLAGE_LINE = "0.40"
ANN_C_PIER         = "#1565C0"
ANN_C_GROIN        = "#B71C1C"

# CoastSat LOESS lines keyed by window size
CS_WINDOW_COLORS = {
     7: "#6BAED6",   # medium sky blue  — 7-domain LOESS
    10: "#08519C",   # deep ocean blue  — 10-domain LOESS
}
CS_WINDOW_COLOR_DEFAULT = "#4A7C8E"

CS_RAW_COLOR = "#5BA3C9"   # individual transect scatter


# =============================================================================
# SECTION 6: GEOGRAPHIC ANNOTATION DATA  (matches run script exactly)
# =============================================================================

ANN_TOWN_SPANS    = {
    "Buxton":      (7,  8),
    "Avon":        (21, 31),
    "Tri-Village": (68, 83),   # Salvo / Waves / Rodanthe
}
ANN_VILLAGE_LINES = {"Salvo": 69, "Waves": 74, "Rodanthe": 80}

ANN_PIER_LABEL_Y  = 0.76
ANN_GROIN_LABEL_Y = 0.68

ANN_PIERS  = {
    "Avon Pier":     (26, ANN_PIER_LABEL_Y),
    "Rodanthe Pier": (79, ANN_PIER_LABEL_Y),
}
ANN_GROINS        = {"Buxton Groin": 5.5}
ANN_WIMBLE_SHOALS = (60, 74)

LABEL_ACCRETION_Y = None   # None = auto-computed mid-point of positive side
LABEL_EROSION_Y   = None   # None = auto-computed mid-point of negative side


# =============================================================================
# HELPER FUNCTIONS
# =============================================================================

def load_cascade_from_npz(npz_path):
    """
    Load a saved CASCADE object from an NPZ file.

    Standard CASCADE save format:
        np.load(path, allow_pickle=True)["cascade.npy"].item()
    """
    if not os.path.exists(npz_path):
        raise FileNotFoundError(f"NPZ not found: {npz_path}")

    print(f"  Loading: {os.path.basename(os.path.dirname(npz_path))}/")
    data = np.load(npz_path, allow_pickle=True)

    if "cascade.npy" in data.files:
        cascade = data["cascade.npy"].item()
    elif len(data.files) == 1:
        cascade = data[data.files[0]].item()
        print(f"  Note: used NPZ key '{data.files[0]}' (expected 'cascade.npy')")
    else:
        raise KeyError(
            f"Cannot find 'cascade.npy' in {npz_path}. "
            f"Available keys: {data.files}"
        )
    data.close()
    return cascade


def get_x_s_TS(b3d):
    """Extract shoreline time series from a Barrier3D object."""
    if hasattr(b3d, "x_s_TS"):
        return np.asarray(b3d.x_s_TS, dtype=float)
    if hasattr(b3d, "_x_s_TS"):
        return np.asarray(b3d._x_s_TS, dtype=float)
    raise AttributeError(
        "No shoreline time series found (x_s_TS / _x_s_TS) on Barrier3D object."
    )


def build_shoreline_matrix(cascade, to_meters=True):
    """Build [time × domain] shoreline matrix (matches run script exactly)."""
    try:
        b3d_list = cascade.barrier3d
    except AttributeError:
        b3d_list = cascade._barrier3d

    ndom = len(b3d_list)
    nt   = len(get_x_s_TS(b3d_list[0]))
    shoreline = np.zeros((nt, ndom), dtype=float)
    for j in range(ndom):
        shoreline[:, j] = get_x_s_TS(b3d_list[j])
    if to_meters:
        shoreline *= 10.0   # dam → m
    return shoreline


def build_relative_shoreline_change_matrix(cascade, to_meters=True, flip_sign=True):
    """
    [time × domain] relative shoreline change from t=0.
    flip_sign=True: positive = accretion, negative = erosion.
    """
    sl_m = build_shoreline_matrix(cascade, to_meters=to_meters)
    change = sl_m - sl_m[0, :]
    if flip_sign:
        change *= -1.0
    return change


def build_nourishment_arrays():
    """
    Build year-keyed nourishment on/volume dicts for BN bar panels.
    Any event year outside [START_YEAR, END_YEAR] is silently skipped.
    """
    on_by_year  = {yr: np.zeros(TOTAL_DOMAINS)  for yr in range(START_YEAR, END_YEAR + 1)}
    vol_by_year = {yr: [0.0] * TOTAL_DOMAINS    for yr in range(START_YEAR, END_YEAR + 1)}

    for gis_id, volumes_m3 in HAT_BN_VOLUME_BY_DOMAIN.items():
        pad_idx = _gis_to_pad(gis_id)
        if not (0 <= pad_idx < TOTAL_DOMAINS):
            continue
        for year, total_m3 in zip(HAT_BN_YEARS, volumes_m3):
            if START_YEAR <= year <= END_YEAR and total_m3 > 0:
                on_by_year[year][pad_idx]  = 1
                vol_by_year[year][pad_idx] = float(total_m3) / DOMAIN_LENGTH_M

    active = [yr for yr in range(START_YEAR, END_YEAR + 1) if np.any(on_by_year[yr] == 1)]
    if active:
        print("Nourishment events in this period:")
        for yr in active:
            doms = [_pad_to_gis(i) for i in np.where(on_by_year[yr] == 1)[0]
                    if _pad_to_gis(i) is not None]
            print(f"  {yr}: GIS domains {doms}")
    else:
        print("No nourishment events in this period.")

    return on_by_year, vol_by_year


def _bn_group_labels(active_gis, bn_real):
    """
    Group consecutive nourished GIS domains; return one label descriptor per group.

    Instead of a cramped label above every individual bar (which overlaps badly
    when 10 Buxton or 4 Avon domains are all adjacent), this detects runs of
    consecutive nourished domains and returns a single centered annotation per run.

    Parameters
    ----------
    active_gis : list[int]   GIS domain IDs with non-zero BN this year (sorted)
    bn_real    : array[90]   BN volume (m³) per real domain (index = gis_d - 1)

    Returns
    -------
    list of dicts with keys:
        x_center  float   — GIS domain x position to center the label over
        top_vol   float   — max volume in group (drives the label's y position)
        label     str     — e.g. "D6–15\n91.7k/dom" or "D85\n309.6k"
    """
    if not active_gis:
        return []

    sorted_gis = sorted(active_gis)
    groups = []
    current = [sorted_gis[0]]

    for gis_d in sorted_gis[1:]:
        if gis_d == current[-1] + 1:         # consecutive → extend current run
            current.append(gis_d)
        else:                                  # gap → close run, start new one
            groups.append(current[:])
            current = [gis_d]
    groups.append(current)

    result = []
    for grp in groups:
        vols     = [bn_real[d - FIRST_FILE_NUMBER] for d in grp]
        top_vol  = max(vols)
        x_center = (grp[0] + grp[-1]) / 2.0
        vol_k    = vols[0] / 1000.0            # all domains in a project share the same volume

        if len(grp) == 1:
            label = f"D{grp[0]}\n{vol_k:.0f}k"
        else:
            label = f"D{grp[0]}–{grp[-1]}\n{vol_k:.0f}k/dom"

        result.append(dict(x_center=x_center, top_vol=top_vol, label=label))

    return result


# =============================================================================
# COASTSAT LOESS PIPELINE  (copied from HAT_hindcast_1984_2024_old version.py)
# =============================================================================

def estimate_transect_spacing(along_coast_m):
    """Median spacing between consecutive transects in metres (positive diffs only)."""
    arr   = np.sort(along_coast_m)
    diffs = np.diff(arr)
    pos   = diffs[diffs > 0]
    return float(np.median(pos)) if len(pos) else 50.0


def load_transect_data(ds):
    """
    Load individual transect LRR values from transect_lrr_full.csv.
    Derives along-coast distance by spreading each domain's transects evenly
    across its 500 m band.

    Returns domain_ids, lrr_values, along_coast_m — all None on load failure.
    """
    csv_path   = ds["csv"]
    domain_col = ds["domain_col"]
    rate_col   = ds["rate_col"]
    id_col     = ds.get("transect_id_col", "transect_id")

    if not os.path.exists(csv_path):
        print(f"  ⚠️  Transect CSV not found: {csv_path}")
        return None, None, None

    df = pd.read_csv(csv_path)
    df.columns = [c.split(".csv")[-1] if ".csv" in c else c for c in df.columns]

    for col in [domain_col, rate_col]:
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors="coerce")

    df = df.dropna(subset=[domain_col, rate_col])
    df[domain_col] = df[domain_col].astype(int)
    df = df[(df[domain_col] >= FIRST_FILE_NUMBER) & (df[domain_col] <= LAST_FILE_NUMBER)]

    sort_cols = [domain_col, id_col] if id_col in df.columns else [domain_col]
    df = df.sort_values(sort_cols).reset_index(drop=True)

    df["_rank"] = df.groupby(domain_col).cumcount()
    df["_n"]    = df.groupby(domain_col)[domain_col].transform("count")
    df["along_coast_m"] = (
        (df[domain_col] - 1) * DOMAIN_SPACING_M
        + (df["_rank"] + 0.5) * (DOMAIN_SPACING_M / df["_n"])
    )
    df = df.drop(columns=["_rank", "_n"])

    domain_ids    = df[domain_col].values.astype(int)
    lrr_values    = df[rate_col].values.astype(float)
    along_coast_m = df["along_coast_m"].values.astype(float)

    spacing = estimate_transect_spacing(along_coast_m)
    print(f"  ✓ {ds['label']}: {len(df)} transects  "
          f"est. spacing {spacing:.0f} m  "
          f"LRR range {np.nanmin(lrr_values):+.2f}–{np.nanmax(lrr_values):+.2f} m/yr")
    return domain_ids, lrr_values, along_coast_m


def loess_smooth_transect_to_domains(along_coast_m, lrr, domain_ids, window_domains):
    """
    Apply LOESS at transect resolution using physical along-coast distance (m),
    then aggregate smoothed values to GIS domain resolution.

    Returns gis_x, smoothed, frac — all None on failure.
    """
    window_km = window_domains * DOMAIN_SPACING_M / 1000.0
    spacing_m = estimate_transect_spacing(along_coast_m)
    n         = len(along_coast_m)
    frac      = float(np.clip((window_km * 1000.0 / spacing_m) / n, 0.02, 1.0))

    valid = np.isfinite(lrr)
    if valid.sum() < 5:
        print(f"  ⚠️  Too few valid transects ({valid.sum()}) — skipping LOESS")
        return None, None, frac

    result            = lowess(lrr[valid], along_coast_m[valid], frac=frac, return_sorted=True)
    smoothed_t        = np.full(n, np.nan)
    smoothed_t[valid] = np.interp(along_coast_m[valid], result[:, 0], result[:, 1])

    dom_agg = (pd.DataFrame({"domain": domain_ids, "smoothed": smoothed_t})
                 .groupby("domain")["smoothed"].mean()
                 .dropna())

    return dom_agg.index.values.astype(int), dom_agg.values, frac


def splice_loess_with_raw_south(win_gis_x, win_smoothed,
                                transect_domain_ids, transect_lrr_values,
                                skip_n=LOESS_SKIP_SOUTHERN_DOMAINS,
                                is_widest_window=False):
    """
    Return (plot_x, plot_y) for a LOESS window with optional southern splice.

    Widest window + skip_n > 0: domains 1–skip_n use raw per-domain means.
    All other windows: line simply starts at domain skip_n+1.
    """
    if skip_n == 0:
        return win_gis_x, win_smoothed

    # LOESS portion: domains strictly north of skip_n
    mask   = win_gis_x > skip_n
    lx, ly = win_gis_x[mask], win_smoothed[mask]
    return lx, ly


def load_all_coastsat(active_start_year):
    """
    Load all COASTSAT_DATASETS, apply LOESS at transect resolution.

    Returns a list of cs_series dicts matching the run-script structure.
    active_start_year controls which dataset gets full-opacity 'active' styling.
    """
    print("\nLoading CoastSat transect data...")
    cs_series = []
    for ds in COASTSAT_DATASETS:
        domain_ids, lrr_values, along_coast_m = load_transect_data(ds)
        if domain_ids is None:
            continue
        windows = []
        for w in LOESS_WINDOW_DOMAINS:
            gis_x, smoothed, frac = loess_smooth_transect_to_domains(
                along_coast_m, lrr_values, domain_ids, w
            )
            if gis_x is None:
                continue
            print(f"  ✓ LOESS window={w} dom "
                  f"({w * DOMAIN_SPACING_M / 1000.0:.1f} km)  "
                  f"frac={frac:.3f}  ({ds['label']})")
            windows.append(dict(window=w, gis_x=gis_x, smoothed=smoothed, frac=frac))
        cs_series.append(dict(
            label                = ds["label"],
            period_start         = ds["period_start"],
            active               = (ds["period_start"] == active_start_year),
            transect_domains     = domain_ids,
            transect_rates       = lrr_values,
            transect_along_coast = along_coast_m,
            windows              = windows,
        ))
    return cs_series


# =============================================================================
# GEOGRAPHIC ANNOTATION FUNCTIONS  (matches run script exactly)
# =============================================================================

def add_geographic_annotations(ax):
    """
    Add all geographic reference annotations to an axis.

    Layer order (bottom → top):
      1. Wimble Shoals influence zone  (hatched amber fill, bottom label)
      2. Community shaded spans        (steel-blue fill, top labels)
      3. Village center lines          (dashed gray,    y=0.84)
      4. Pier lines                    (dash-dot blue,  y=ANN_PIER_LABEL_Y)
      5. Groin lines                   (dotted red,     y=ANN_GROIN_LABEL_Y)

    X-axis must be in GIS domain IDs (1–90).
    Y-axis label positions use blended axes-fraction coordinates (data x, axes y).
    """
    trans = blended_transform_factory(ax.transData, ax.transAxes)

    # 1. Wimble Shoals influence zone
    wlo, whi = ANN_WIMBLE_SHOALS
    ax.axvspan(wlo - 0.5, whi + 0.5,
               color=ANN_C_WIMBLE, alpha=0.10, zorder=0,
               hatch="///", edgecolor=ANN_C_WIMBLE, linewidth=0)
    ax.text((wlo + whi) / 2.0, 0.04,
            "Wimble Shoals\nPosition", transform=trans,
            ha="center", va="bottom", fontsize=7, color="#7A5800", style="italic",
            bbox=dict(boxstyle="round,pad=0.2", fc="white", ec="none", alpha=0.80))

    # 2. Community / town spans
    for span_label, (d_lo, d_hi) in ANN_TOWN_SPANS.items():
        ax.axvspan(d_lo - 0.5, d_hi + 0.5,
                   color=ANN_C_TOWN_SPAN, alpha=0.14, zorder=0)
        ax.text((d_lo + d_hi) / 2.0, 0.90,
                span_label, transform=trans,
                ha="center", va="top", fontsize=8, color="0.25", fontweight="bold",
                bbox=dict(boxstyle="round,pad=0.2", fc="white", ec="none", alpha=0.85))

    # 3. Village center lines (within Tri-Village span)
    for vname, dom in ANN_VILLAGE_LINES.items():
        ax.axvline(dom, color=ANN_C_VILLAGE_LINE, lw=0.9, ls="--", alpha=0.65, zorder=1)
        ax.text(dom, 0.84, vname, transform=trans,
                ha="center", va="top", fontsize=7.5, color="0.30",
                bbox=dict(boxstyle="round,pad=0.15", fc="white", ec="none", alpha=0.80))

    # 4. Pier lines
    for pname, (dom, lbl_y) in ANN_PIERS.items():
        ax.axvline(dom, color=ANN_C_PIER, lw=1.0, ls="-.", alpha=0.80, zorder=2)
        ax.text(dom, lbl_y, pname, transform=trans,
                ha="center", va="top", fontsize=7, color=ANN_C_PIER, rotation=90,
                bbox=dict(boxstyle="round,pad=0.15", fc="white", ec="none", alpha=0.80))

    # 5. Groin lines
    for gname, dom in ANN_GROINS.items():
        ax.axvline(dom, color=ANN_C_GROIN, lw=1.1, ls=":", alpha=0.85, zorder=2)
        ax.text(dom, ANN_GROIN_LABEL_Y, gname, transform=trans,
                ha="center", va="top", fontsize=7, color=ANN_C_GROIN, rotation=90,
                bbox=dict(boxstyle="round,pad=0.15", fc="white", ec="none", alpha=0.80))


def annotation_legend_handles():
    """Return proxy legend artists for the geographic annotation layers."""
    return [
        Patch(fc=ANN_C_TOWN_SPAN, alpha=0.30, label="Community"),
        Patch(fc=ANN_C_WIMBLE, alpha=0.25, hatch="///",
              edgecolor=ANN_C_WIMBLE, linewidth=0, label="Wimble Shoals position"),
        Line2D([0], [0], color=ANN_C_VILLAGE_LINE, lw=0.9, ls="--", label="Village center"),
        Line2D([0], [0], color=ANN_C_PIER, lw=1.0, ls="-.", label="Pier"),
        Line2D([0], [0], color=ANN_C_GROIN, lw=1.1, ls=":", label="Groin"),
    ]


# =============================================================================
# YEARLY RELATIVE SHORELINE CHANGE + BN BAR  (GIF)
# =============================================================================

def plot_yearly_relative_shoreline_and_bn(
    cascades_by_label,
    hist_nourish_volume_by_year,
    cs_series,
    output_dir,
    filename_prefix,
    flip_sign=True,
    make_gif=True,
    gif_duration_seconds=4,
):
    """
    One PNG per year: relative shoreline change from t=0 (upper panel)
    + historical BN volume bar (lower panel).
    X-axis: GIS domain IDs 1–90.
    """
    yearly_dir = os.path.join(output_dir, filename_prefix)
    os.makedirs(yearly_dir, exist_ok=True)

    gis_ids = np.arange(FIRST_FILE_NUMBER, LAST_FILE_NUMBER + 1)   # 1..90

    # Build relative change matrices (real domains only)
    sc_real_by_label = {}
    max_nt = None
    for label, cascade in cascades_by_label.items():
        sc_full = build_relative_shoreline_change_matrix(cascade, to_meters=True,
                                                         flip_sign=flip_sign)
        sc_real = sc_full[:, START_REAL_INDEX:END_REAL_INDEX]   # [nt, 90]
        sc_real_by_label[label] = sc_real
        max_nt = sc_real.shape[0] if max_nt is None else min(max_nt, sc_real.shape[0])

    if max_nt is None:
        print("No cascades for yearly relative-shoreline plotting.")
        return []

    # Fixed global y-limits
    all_vals = [v for sc in sc_real_by_label.values()
                for v in sc[:max_nt, :][np.isfinite(sc[:max_nt, :])]]
    y_min, y_max = (np.nanmin(all_vals), np.nanmax(all_vals)) if all_vals else (-10.0, 10.0)
    y_pad = max(0.15 * (y_max - y_min), 5.0)
    y_min -= y_pad
    y_max += y_pad

    # BN bar scale
    all_bn_m3 = [
        v for yr in range(START_YEAR, END_YEAR + 1)
        if yr in hist_nourish_volume_by_year
        for v in (np.asarray(hist_nourish_volume_by_year[yr], dtype=float) * DOMAIN_LENGTH_M)
        if v > 0
    ]
    max_bn_volume = max(all_bn_m3) if all_bn_m3 else 1.0

    png_files = []

    for t in range(max_nt):
        model_year = START_YEAR + t

        # BN this year — real domains only
        bn_full = (
            np.asarray(hist_nourish_volume_by_year[model_year], dtype=float) * DOMAIN_LENGTH_M
            if model_year in hist_nourish_volume_by_year
            else np.zeros(TOTAL_DOMAINS)
        )
        bn_real     = bn_full[START_REAL_INDEX:END_REAL_INDEX]   # len=90
        active_gis  = gis_ids[bn_real > 0].tolist()

        fig, (ax, ax_bn) = plt.subplots(
            2, 1, figsize=(14, 7), sharex=True,
            gridspec_kw={"height_ratios": [3.2, 1.15]},
            constrained_layout=True,
        )
        fig.patch.set_facecolor("white")
        ax.set_facecolor("white")

        # Geographic annotations (drawn before data so data renders on top)
        ax.set_xlim(FIRST_FILE_NUMBER - 0.5, LAST_FILE_NUMBER + 0.5)
        ax.set_ylim(y_min, y_max)
        add_geographic_annotations(ax)

        # Shoreline change curves
        for label, sc_real in sc_real_by_label.items():
            ax.plot(gis_ids, sc_real[t, :], linewidth=2.0, label=label, zorder=6)

        ax.axhline(0.0, linestyle="--", linewidth=1.2, color="#2c2c2c", alpha=0.65, zorder=3)

        # BN domain markers on the upper panel
        for gis_d in active_gis:
            padded_i = _gis_to_pad(gis_d)
            ax.axvline(gis_d, linestyle=":", linewidth=1.0, color="red", alpha=0.45, zorder=4)
        if active_gis:
            for label, sc_real in sc_real_by_label.items():
                bn_y = [sc_real[t, gis_d - FIRST_FILE_NUMBER] for gis_d in active_gis]
                ax.scatter(active_gis, bn_y, s=35, marker="o", color="red",
                           edgecolor="black", linewidth=0.6, zorder=10,
                           label="Historical BN domain")
                break  # only mark on first scenario

        # Accretion / Erosion labels (axes-fraction y, right-hand side)
        ybot, ytop = ax.get_ylim()
        zero_frac  = (0 - ybot) / (ytop - ybot)
        acc_y = zero_frac + (1 - zero_frac) / 2
        ero_y = zero_frac / 2
        ax.text(1.01, acc_y, "Accretion ▲", transform=ax.transAxes,
                fontsize=9, color="#555555", ha="left", va="center", style="italic")
        ax.text(1.01, ero_y, "Erosion ▼",   transform=ax.transAxes,
                fontsize=9, color="#555555", ha="left", va="center", style="italic")

        ax.set_ylabel(
            f"Relative shoreline change from {START_YEAR} (m)",
            fontsize=11, fontweight="bold",
        )
        ax.set_title(
            f"Hatteras Island  |  Relative shoreline change  |  {model_year}",
            fontsize=13, fontweight="bold", pad=10,
        )
        ax.grid(True, which="major", linestyle=":", linewidth=0.6, alpha=0.4, color="gray")
        ax.spines[["top", "right"]].set_visible(False)
        ax.legend(loc="upper left", fontsize=8.5, frameon=True, framealpha=0.92)

        # BN bar panel
        ax_bn.bar(gis_ids, bn_real, width=0.8, color="tab:blue", alpha=0.85)
        total_bn = np.sum(bn_real)
        ax_bn.set_title(
            (f"Historical BN in {model_year}: YES | Total = {total_bn:,.0f} m³")
            if active_gis
            else f"Historical BN in {model_year}: NO",
            fontsize=11, fontweight="bold",
        )
        ax_bn.set_ylabel("BN volume\n(m³/domain)", fontsize=9, fontweight="bold")
        ax_bn.set_xlabel("GIS Domain ID (1–90)  |  1 = Cape Point  ·  90 = near Rodanthe",
                         fontsize=10, fontweight="bold")
        ax_bn.set_ylim(0, max_bn_volume * 1.25)
        ax_bn.grid(True, axis="y", linestyle=":", linewidth=0.6, alpha=0.4)

        xticks = np.arange(FIRST_FILE_NUMBER, LAST_FILE_NUMBER + 1, DOMAIN_TICK_STEP)
        ax_bn.set_xticks(xticks)
        ax_bn.set_xticklabels([str(i) for i in xticks], rotation=45, ha="right", fontsize=9)

        # One centered label per consecutive group of nourished domains
        for grp in _bn_group_labels(active_gis, bn_real):
            ax_bn.text(grp["x_center"],
                       grp["top_vol"] + max_bn_volume * 0.030,
                       grp["label"],
                       ha="center", va="bottom", fontsize=8.5,
                       color="0.20", fontweight="bold")

        fig_out = os.path.join(yearly_dir, f"{filename_prefix}_{model_year}.png")
        fig.savefig(fig_out, dpi=200, bbox_inches="tight", facecolor="white")
        plt.close(fig)
        png_files.append(fig_out)

    print(f"\nSaved {len(png_files)} yearly relative-shoreline PNGs → {yearly_dir}")
    if make_gif and png_files:
        _save_gif(png_files, output_dir, filename_prefix, gif_duration_seconds)
    return png_files


# =============================================================================
# YEARLY ABSOLUTE SHORELINE POSITION + BN BAR  (GIF)
# =============================================================================

def plot_yearly_absolute_shoreline_and_bn(
    cascades_by_label,
    hist_nourish_volume_by_year,
    output_dir,
    filename_prefix,
    make_gif=True,
    gif_duration_seconds=4,
):
    """
    One PNG per year: raw x_s position with ocean fill (upper panel)
    + historical BN volume bar (lower panel).
    X-axis: GIS domain IDs 1–90.

    Sign convention:
        lower x_s  = seaward / ocean side
        higher x_s = landward / back-barrier side
    Ocean fill drawn from y_min up to shoreline curve.
    """
    yearly_dir = os.path.join(output_dir, filename_prefix)
    os.makedirs(yearly_dir, exist_ok=True)

    gis_ids = np.arange(FIRST_FILE_NUMBER, LAST_FILE_NUMBER + 1)

    sl_real_by_label = {}
    max_nt = None
    for label, cascade in cascades_by_label.items():
        sl_full = build_shoreline_matrix(cascade, to_meters=True)
        sl_real = sl_full[:, START_REAL_INDEX:END_REAL_INDEX]   # [nt, 90]
        sl_real_by_label[label] = sl_real
        max_nt = sl_real.shape[0] if max_nt is None else min(max_nt, sl_real.shape[0])

    if max_nt is None:
        return []

    all_vals = [v for sl in sl_real_by_label.values()
                for v in sl[:max_nt, :][np.isfinite(sl[:max_nt, :])]]
    y_min, y_max = (np.nanmin(all_vals), np.nanmax(all_vals)) if all_vals else (-6000.0, -5000.0)
    y_pad = max(0.10 * (y_max - y_min), 50.0)
    y_min -= y_pad
    y_max += y_pad

    all_bn_m3 = [
        v for yr in range(START_YEAR, END_YEAR + 1)
        if yr in hist_nourish_volume_by_year
        for v in (np.asarray(hist_nourish_volume_by_year[yr], dtype=float) * DOMAIN_LENGTH_M)
        if v > 0
    ]
    max_bn_volume = max(all_bn_m3) if all_bn_m3 else 1.0

    # Reference curve for fill (last key = most management-rich scenario)
    fill_label = list(sl_real_by_label.keys())[-1]

    png_files = []

    for t in range(max_nt):
        model_year = START_YEAR + t

        bn_full = (
            np.asarray(hist_nourish_volume_by_year[model_year], dtype=float) * DOMAIN_LENGTH_M
            if model_year in hist_nourish_volume_by_year
            else np.zeros(TOTAL_DOMAINS)
        )
        bn_real    = bn_full[START_REAL_INDEX:END_REAL_INDEX]
        active_gis = gis_ids[bn_real > 0].tolist()

        fig, (ax, ax_bn) = plt.subplots(
            2, 1, figsize=(14, 7), sharex=True,
            gridspec_kw={"height_ratios": [3.2, 1.15]},
            constrained_layout=True,
        )
        fig.patch.set_facecolor("white")
        ax.set_facecolor("white")

        ax.set_xlim(FIRST_FILE_NUMBER - 0.5, LAST_FILE_NUMBER + 0.5)
        ax.set_ylim(y_min, y_max)
        add_geographic_annotations(ax)

        # Ocean fill
        ocean_curve = sl_real_by_label[fill_label][t, :]
        ax.fill_between(gis_ids, y_min, ocean_curve,
                        color="lightskyblue", alpha=0.22, linewidth=0, zorder=0,
                        label="Ocean / seaward  (lower $x_s$)")
        ax.fill_between(gis_ids, ocean_curve, y_max,
                        color="tan", alpha=0.06, linewidth=0, zorder=0,
                        label="Landward / back-barrier  (higher $x_s$)")

        # Shoreline position curves
        for label, sl_real in sl_real_by_label.items():
            ax.plot(gis_ids, sl_real[t, :], linewidth=2.2, label=label, zorder=6)

        # BN markers
        for gis_d in active_gis:
            ax.axvline(gis_d, linestyle=":", linewidth=1.0, color="red", alpha=0.42, zorder=4)
        if active_gis:
            for label, sl_real in sl_real_by_label.items():
                bn_y = [sl_real[t, gis_d - FIRST_FILE_NUMBER] for gis_d in active_gis]
                ax.scatter(active_gis, bn_y, s=36, marker="o", color="red",
                           edgecolor="black", linewidth=0.6, zorder=12,
                           label="Historical BN domain")
                break

        # Orientation labels
        ax.text(0.012, 0.06, "Ocean / seaward\n(smaller $x_s$)",
                transform=ax.transAxes, fontsize=9.5, fontweight="bold",
                color="navy", ha="left", va="bottom",
                bbox=dict(boxstyle="round,pad=0.22", fc="white", ec="lightskyblue", alpha=0.88),
                zorder=30)
        ax.text(0.012, 0.94, "Landward / back-barrier\n(larger $x_s$)",
                transform=ax.transAxes, fontsize=9.5, fontweight="bold",
                color="saddlebrown", ha="left", va="top",
                bbox=dict(boxstyle="round,pad=0.22", fc="white", ec="tan", alpha=0.88),
                zorder=30)

        ax.set_ylabel(
            "Absolute shoreline position  $x_s$ (m)\n"
            "Higher $x_s$ = landward retreat  |  Lower $x_s$ = seaward advance",
            fontsize=11, fontweight="bold",
        )
        ax.set_title(
            f"Hatteras Island  |  Absolute shoreline position  |  {model_year}",
            fontsize=13, fontweight="bold", pad=10,
        )
        ax.grid(True, which="major", linestyle=":", linewidth=0.6, alpha=0.4, color="gray")
        ax.spines[["top", "right"]].set_visible(False)
        ax.legend(loc="upper right", fontsize=8.5, frameon=True, framealpha=0.92)

        # BN bar panel
        ax_bn.bar(gis_ids, bn_real, width=0.8, color="tab:blue", alpha=0.85)
        total_bn = np.sum(bn_real)
        ax_bn.set_title(
            (f"Historical BN in {model_year}: YES | Total = {total_bn:,.0f} m³")
            if active_gis else f"Historical BN in {model_year}: NO",
            fontsize=11, fontweight="bold",
        )
        ax_bn.set_ylabel("BN volume\n(m³/domain)", fontsize=9, fontweight="bold")
        ax_bn.set_xlabel("GIS Domain ID (1–90)  |  1 = Cape Point  ·  90 = near Rodanthe",
                         fontsize=10, fontweight="bold")
        ax_bn.set_ylim(0, max_bn_volume * 1.25)
        ax_bn.grid(True, axis="y", linestyle=":", linewidth=0.6, alpha=0.4)

        xticks = np.arange(FIRST_FILE_NUMBER, LAST_FILE_NUMBER + 1, DOMAIN_TICK_STEP)
        ax_bn.set_xticks(xticks)
        ax_bn.set_xticklabels([str(i) for i in xticks], rotation=45, ha="right", fontsize=9)

        # One centered label per consecutive group of nourished domains
        for grp in _bn_group_labels(active_gis, bn_real):
            ax_bn.text(grp["x_center"],
                       grp["top_vol"] + max_bn_volume * 0.030,
                       grp["label"],
                       ha="center", va="bottom", fontsize=8.5,
                       color="0.20", fontweight="bold")

        fig_out = os.path.join(yearly_dir, f"{filename_prefix}_{model_year}.png")
        fig.savefig(fig_out, dpi=200, bbox_inches="tight", facecolor="white")
        plt.close(fig)
        png_files.append(fig_out)

    print(f"\nSaved {len(png_files)} yearly absolute-shoreline PNGs → {yearly_dir}")
    if make_gif and png_files:
        _save_gif(png_files, output_dir, filename_prefix, gif_duration_seconds)
    return png_files


def _save_gif(png_files, output_dir, filename_prefix, gif_duration_seconds):
    """Assemble list of PNGs into an animated GIF."""
    try:
        import imageio.v2 as imageio
        images  = [imageio.imread(f) for f in png_files]
        gif_out = os.path.join(output_dir, f"{filename_prefix}.gif")
        imageio.mimsave(gif_out, images, duration=gif_duration_seconds, loop=0)
        print(f"Saved GIF: {gif_out}")
    except ImportError:
        print("imageio not installed — GIF skipped.  (pip install imageio)")


# =============================================================================
# PUBLICATION FIGURE  (matches run-script annotated figure exactly)
# =============================================================================

def plot_publication_rate_figure(
    rate_profiles_by_label,
    cs_series,
    run_name,
    output_dir,
):
    """
    Publication-quality rate comparison figure.

    Matches the annotated figure produced at the end of main() in the run
    script: model line (warm orange) + CoastSat scatter + multi-window LOESS
    + full geographic annotation layer.
    X-axis: GIS domain IDs 1–90.
    """
    gis_ids     = np.arange(FIRST_FILE_NUMBER, LAST_FILE_NUMBER + 1)
    widest_win  = max(LOESS_WINDOW_DOMAINS)
    data_handles = []

    fig, ax = plt.subplots(figsize=(13, 6.0), constrained_layout=True)
    fig.patch.set_facecolor("white")
    ax.set_facecolor("white")

    # Geographic annotations first (data renders on top)
    add_geographic_annotations(ax)

    # CoastSat scatter + multi-window LOESS lines
    for cs in cs_series:
        is_active = cs["active"]
        if not is_active and not PLOT_REFERENCE_PERIOD:
            continue

        scatter_x = cs["transect_along_coast"] / DOMAIN_SPACING_M + FIRST_FILE_NUMBER
        if PLOT_RAW_LRR:
            raw_alpha = RAW_LRR_SCATTER_ALPHA if is_active else RAW_LRR_SCATTER_ALPHA * 0.35
            raw_lbl   = f"{cs['label']} — transect LRR" if is_active else None
            ax.scatter(scatter_x, cs["transect_rates"],
                       color=CS_RAW_COLOR, s=RAW_LRR_SCATTER_SIZE,
                       alpha=raw_alpha, zorder=1, linewidths=0, label=raw_lbl)
            if is_active:
                data_handles.append(
                    Line2D([0], [0], color=CS_RAW_COLOR, marker=".", ms=5,
                           ls="none", alpha=RAW_LRR_SCATTER_ALPHA, label=raw_lbl)
                )

        for idx, win in enumerate(cs["windows"]):
            cs_color = CS_WINDOW_COLORS.get(win["window"], CS_WINDOW_COLOR_DEFAULT)
            lw_base, ls, alpha_factor = (
                LOESS_WINDOW_STYLES[idx] if idx < len(LOESS_WINDOW_STYLES)
                else (1.5, "-", 0.80)
            )
            is_widest = (win["window"] == widest_win)
            gis_x, rate_y = splice_loess_with_raw_south(
                win["gis_x"], win["smoothed"],
                cs["transect_domains"], cs["transect_rates"],
                is_widest_window=is_widest,
            )
            lbl = f"{cs['label']} — LOESS {win['window']}-dom"
            if is_active:
                if is_widest:
                    ax.fill_between(gis_x, rate_y, 0,
                                    where=(rate_y < 0), alpha=0.14, color=cs_color,
                                    interpolate=True)
                    ax.fill_between(gis_x, rate_y, 0,
                                    where=(rate_y >= 0), alpha=0.10, color=cs_color,
                                    interpolate=True)
                ax.plot(gis_x, rate_y, color=cs_color, lw=lw_base,
                        ls=ls, alpha=alpha_factor, zorder=5, label=lbl)
                data_handles.append(
                    Line2D([0], [0], color=cs_color, lw=lw_base, ls=ls,
                           alpha=alpha_factor, label=lbl)
                )
            else:
                ax.plot(gis_x, rate_y, color=cs_color, lw=lw_base * 0.85,
                        ls=ls, alpha=0.40 * alpha_factor, zorder=4)
                data_handles.append(
                    Line2D([0], [0], color=cs_color, lw=lw_base * 0.85, ls=ls,
                           alpha=0.40 * alpha_factor, label=lbl + " (ref)")
                )

    # Model line(s) — one per loaded NPZ
    for label, rate_full in rate_profiles_by_label.items():
        real_rate = rate_full[START_REAL_INDEX:END_REAL_INDEX]
        ax.plot(gis_ids, real_rate,
                color=ANN_MODEL_COLOR, linewidth=2.4, zorder=6, label=label)
        data_handles.append(
            Line2D([0], [0], color=ANN_MODEL_COLOR, lw=2.4, label=label)
        )

    # Zero reference
    ax.axhline(0, color="#2c2c2c", linewidth=1.2, linestyle="--", alpha=0.65, zorder=3)

    # Axes formatting
    ax.set_xlim(FIRST_FILE_NUMBER - 0.5, LAST_FILE_NUMBER + 0.5)
    ax.xaxis.set_major_locator(ticker.MultipleLocator(10))
    ax.xaxis.set_minor_locator(ticker.MultipleLocator(5))
    ax.yaxis.set_major_locator(ticker.MultipleLocator(1))
    ax.tick_params(axis="both", which="major", labelsize=11, direction="in", length=5)
    ax.tick_params(axis="both", which="minor", direction="in", length=3)
    ax.grid(True, which="major", linestyle=":", linewidth=0.6, alpha=0.4, color="gray")
    ax.spines[["top", "right"]].set_visible(False)
    ax.spines[["left", "bottom"]].set_linewidth(1.2)

    # Y-limits from real domain values
    all_vals = []
    for rate_full in rate_profiles_by_label.values():
        r = rate_full[START_REAL_INDEX:END_REAL_INDEX]
        all_vals.extend(r[np.isfinite(r)])
    for cs in cs_series:
        for win in cs["windows"]:
            all_vals.extend(win["smoothed"][np.isfinite(win["smoothed"])])
    if all_vals:
        y_lo, y_hi = np.nanmin(all_vals), np.nanmax(all_vals)
        y_pad = (y_hi - y_lo) * 0.06
        ax.set_ylim(y_lo - y_pad, y_hi + y_pad)

    # Accretion / Erosion labels at axes-fraction y (matches run script)
    ybot, ytop = ax.get_ylim()
    zero_frac  = (0 - ybot) / (ytop - ybot)
    acc_y = (LABEL_ACCRETION_Y if LABEL_ACCRETION_Y is not None
             else zero_frac + (1 - zero_frac) / 2)
    ero_y = (LABEL_EROSION_Y   if LABEL_EROSION_Y   is not None
             else zero_frac / 2)
    ax.text(1.0, acc_y, "Accretion ▲", transform=ax.transAxes,
            fontsize=9.5, color="#555555", ha="right", va="center", style="italic")
    ax.text(1.0, ero_y, "Erosion ▼",   transform=ax.transAxes,
            fontsize=9.5, color="#555555", ha="right", va="center", style="italic")

    # Labels, title, orientation text
    ax.set_xlabel("CASCADE Model Domain (500 m alongshore)",
                  fontsize=12, fontweight="bold", labelpad=8)
    ax.set_ylabel("Shoreline Change Rate (m/yr)",
                  fontsize=12, fontweight="bold", labelpad=8)

    ax.text(0.0, 1.01, "← S  |  Cape Hatteras",
            transform=ax.transAxes, fontsize=9, color="#444444",
            ha="left", va="bottom", style="italic", clip_on=False)
    ax.text(1.0, 1.01, "Pea Island  |  N →",
            transform=ax.transAxes, fontsize=9, color="#444444",
            ha="right", va="bottom", style="italic", clip_on=False)

    be_label = SOURCE_SINK_PRESET
    ax.set_title(
        f"Modeled vs Observed Shoreline Change — Hatteras Island, NC  |  "
        f"{START_YEAR}–{END_YEAR}  |  Hs={Hs_LABEL} m  |  BE={be_label}",
        fontsize=12, fontweight="bold", pad=12, color="#1a2a3a",
    )

    # Legend: data + annotation proxies
    ax.legend(
        handles=data_handles + annotation_legend_handles(),
        loc="lower center", bbox_to_anchor=(0.5, 0.02),
        fontsize=9.5, framealpha=0.95, edgecolor="#cccccc",
        frameon=True, ncol=4,
    )

    # Caption
    fig.text(
        0.012, 0.005,
        f"Model: CASCADE  |  Observed: CoastSat LRR averaged per 500-m domain  |  "
        f"SLR={SEA_LEVEL_RISE_RATE * 1000:.1f} mm/yr  |  Run: {run_name}",
        fontsize=8, color="#666666", ha="left", va="bottom", style="italic",
    )

    fig_out = os.path.join(output_dir, f"{run_name}_annotated.png")
    fig.savefig(fig_out, dpi=300, bbox_inches="tight", facecolor="white")
    print(f"\nSaved publication figure: {fig_out}")
    plt.show()
    return fig_out


# =============================================================================
# MAIN
# =============================================================================

def main():
    print("\n" + "=" * 80)
    print("LOADING CASCADE NPZ OUTPUT(S)")
    print("=" * 80)

    cascades_by_label = {}
    for label, npz_path in NPZ_PATHS_BY_LABEL.items():
        try:
            cascade = load_cascade_from_npz(npz_path)
            cascades_by_label[label] = cascade
            print(f"  ✓ {label}")
        except FileNotFoundError as e:
            print(f"  SKIPPED (not found): {label}\n    {e}")
        except Exception as e:
            print(f"  ERROR loading '{label}': {e}")

    if not cascades_by_label:
        print("\nNo CASCADE runs loaded. Check NPZ_PATHS_BY_LABEL paths.")
        sys.exit(1)

    # ── CoastSat LOESS ────────────────────────────────────────────────────────
    cs_series = load_all_coastsat(active_start_year=START_YEAR)

    # ── Nourishment arrays for BN bar panels ──────────────────────────────────
    print("\nBuilding nourishment display arrays...")
    HIST_NOURISH_ON, HIST_NOURISH_VOLUME = build_nourishment_arrays()

    # ── Shoreline change rates ─────────────────────────────────────────────────
    time_span_years = END_YEAR - START_YEAR
    rate_profiles   = {}

    for label, cascade in cascades_by_label.items():
        sl_m  = build_shoreline_matrix(cascade, to_meters=True)
        denom = time_span_years if time_span_years > 0 else max(sl_m.shape[0] - 1, 1)
        rate  = (sl_m[-1, :] - sl_m[0, :]) / float(denom)
        if FLIP_SIGN_MODEL:
            rate *= -1.0
        rate_profiles[label] = rate

    # Run name for file naming
    run_name = f"HAT_{START_YEAR}_{END_YEAR}_npz_analysis"

    # ── Yearly relative shoreline + BN GIF ────────────────────────────────────
    rel_prefix = f"HAT_{START_YEAR}_{END_YEAR}_yearly_relative_shoreline_BN"
    plot_yearly_relative_shoreline_and_bn(
        cascades_by_label      = cascades_by_label,
        hist_nourish_volume_by_year = HIST_NOURISH_VOLUME,
        cs_series              = cs_series,
        output_dir             = OUTPUT_BASE_DIR,
        filename_prefix        = rel_prefix,
        flip_sign              = True,
        make_gif               = MAKE_YEARLY_BN_GIF,
        gif_duration_seconds   = GIF_DURATION_SECONDS,
    )

    # ── Yearly absolute shoreline + BN GIF ────────────────────────────────────
    abs_prefix = f"HAT_{START_YEAR}_{END_YEAR}_yearly_absolute_shoreline_BN"
    plot_yearly_absolute_shoreline_and_bn(
        cascades_by_label      = cascades_by_label,
        hist_nourish_volume_by_year = HIST_NOURISH_VOLUME,
        output_dir             = OUTPUT_BASE_DIR,
        filename_prefix        = abs_prefix,
        make_gif               = MAKE_YEARLY_BN_GIF,
        gif_duration_seconds   = GIF_DURATION_SECONDS,
    )

    # ── Publication rate figure ────────────────────────────────────────────────
    plot_publication_rate_figure(
        rate_profiles_by_label = rate_profiles,
        cs_series              = cs_series,
        run_name               = run_name,
        output_dir             = OUTPUT_BASE_DIR,
    )


if __name__ == "__main__":
    main()
