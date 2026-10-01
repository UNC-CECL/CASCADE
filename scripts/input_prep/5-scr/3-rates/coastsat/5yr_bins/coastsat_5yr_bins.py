"""
The CoastSat LRR in successive 5-year bins, per GIS domain, for each window of the canonical chain.

    python scripts/input_prep/5-scr/3-rates/coastsat/5yr_bins/coastsat_5yr_bins.py

One folder per window under 3-rates/coastsat/5yr_bins/, tables only;
coastsat_5yr_bins_figure.py draws them. Details: scripts/input_prep/5-scr/3-rates/README.md.

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-09-30
"""

# --- CONFIG ------------------------------------------------------------------
# Resolved through hat_observed_rates.py (2026-09-18)
import sys as _sys
from pathlib import Path as _RP
_sys.path.insert(0, str(next(_q for _q in _RP(__file__).resolve().parents
                             if (_q / "pyproject.toml").exists()) / "scripts"))
from site_layer import hat_observed_rates as _obs  # noqa: E402
LOOKUP_CSV    = str(_obs.TRANSECT_DOMAINS / "transect_domain_lookup.csv")
# Anchored on this file 2026-09-12
from pathlib import Path as _Path

# Repo root, found by searching upward
_PATH_REPO = next(_p for _p in _Path(__file__).resolve().parents
                  if (_p / "pyproject.toml").exists())
ROOT_DATA_DIR = str(_obs.COASTSAT_TIMESERIES)
SITE_FILTER   = "usa_NC"
OUTPUT_DIR    = str(_obs.TIMESERIES_LRR)

# Time periods: (file_tag, figure_title, file_stem, start_date, end_date)

# The canonical chain since 2026-09-17; one folder per window, calendar years inclusive
PERIODS = [
    ("1996_2010", "1996-2010", "1996_2010", "1996-01-01", "2010-12-31"),
    ("2010_2024", "2010-2024", "2010_2024", "2010-01-01", "2024-12-31"),
    ("1996_2024", "1996-2024", "1996_2024", "1996-01-01", "2024-12-31"),
]

# Interval sizes (years): one subfolder per entry

# One subfolder per entry; 5 yr is the minimum that averages through storm recovery
INTERVAL_SIZES_YR = [5]

# Minimum observations per transect per bin

# The only gate a transect must pass
MIN_OBS_PER_BIN = 3

# Statistical quality filters (applied before domain aggregation)

# Off by default (1.0 and 0.0): short bins rarely reach p <= 0.10
MAX_PVALUE = 1.0    # 1.0 = disabled (include all fits regardless of p-value)
MIN_R2     = 0.0    # 0.0 = disabled (include all fits regardless of R2)

# Physical plausibility cap (m/yr)

# Per-transect |LRR| above this is dropped before aggregation
MAX_RATE_M_YR = 50.0

# Minimum bin width relative to the nominal interval

# Bins shorter than this fraction of the interval are dropped
MIN_BIN_FRACTION = 0.75

# Domain aggregation method: "mean" or "median"

# 'mean' matches the rest of the pipeline
DOMAIN_AGG = "mean"

# Color palette

# 'custom' runs cool (oldest) to warm (newest)
PALETTE = "custom"
CUSTOM_COLORS = [
    "#d9f0d3",
    "#a6dba0",
    "#5aae61",
    "#1b7837",
    "#00441b",
]

# Line appearance
LINE_WIDTH = 2.0
LINE_ALPHA = 0.88

# Y-axis range for LRR panel

# None = auto (98th percentile clip). Or fix e.g. (-6, 6) for cross-period comparison.
YLIM_LRR = None

# Y-axis range for the overall LRR bar panel
YLIM_OVERALL = None

# CASCADE buffer domains to exclude
BUFFER_DOMAINS = []

# Community span annotations (kept for backward compatibility)

# Kept for backward compatibility; the annotation system below supersedes it
COMMUNITY_SPANS = [
    (7,  8,  "Buxton"),
    (21, 31, "Avon"),
    (68, 83, "Tri-Village"),
]

# Section 4b: geographic annotation styling, as in every CoastSat figure

ANN_TOWN_SPANS = {
    "Buxton":      (7,  8),
    "Avon":        (21, 31),
    "Tri-Village": (68, 83),   # Salvo / Waves / Rodanthe
}
ANN_VILLAGE_LINES = {"Salvo": 69, "Waves": 74, "Rodanthe": 80}

ANN_PIER_LABEL_Y  = 0.76   # rotated label y for piers (0=bottom, 1=top axes fraction)
ANN_GROIN_LABEL_Y = 0.68   # rotated label y for groin lines
ANN_PIERS = {
    "Avon Pier":     (26, ANN_PIER_LABEL_Y),
    "Rodanthe Pier": (79, ANN_PIER_LABEL_Y),
}
ANN_GROINS        = {"Buxton Groin": 5.5}
ANN_WIMBLE_SHOALS = (60, 74)

ANN_C_TOWN_SPAN    = "#90AFC5"   # steel blue for town span shading
ANN_C_WIMBLE       = "#E0A800"   # gold for Wimble Shoals shading
ANN_C_VILLAGE_LINE = "0.40"      # dark gray for village dividers
ANN_C_PIER         = "#1565C0"   # blue for pier markers
ANN_C_GROIN        = "#B71C1C"   # red for groin line

# LOWESS spatial smoothing overlay

# So the smoothed signal is the dominant visual
LOWESS_OVERLAY      = True
LOWESS_ONLY         = True    # if True, raw line drawn at reduced alpha
RAW_ALPHA_SCALE    = 0.25    # multiplier applied to LINE_ALPHA for raw line
                             # when LOWESS_ONLY is True (0 = hide raw entirely)
LOWESS_FRAC         = 7 / 90  # 7-domain window, the model target's
LOWESS_LINE_WIDTH   = 3.0
LOWESS_LINE_ALPHA   = 0.95

# Physical plausibility cap (m/yr)

# 20 m/yr keeps the real extremes at Oregon Inlet and Rodanthe, drops detection errors
MAX_RATE_M_YR = 50.0

FIG_WIDTH  = 26
FIG_HEIGHT = 10
DPI        = 150


import os
import sys
import glob
import warnings
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.cm as cm
import matplotlib.gridspec as gridspec
import matplotlib.lines as mlines
from datetime import datetime

warnings.filterwarnings("ignore")

sys.path.insert(0, str(next(_q for _q in _RP(__file__).resolve().parents
                            if (_q / "pyproject.toml").exists())
                       / "scripts" / "input_prep" / "5-scr" / "lib"))
import scr_paths  # noqa: E402,F401  (5-scr sibling modules onto sys.path)

# Use the identical compute_lrr as every other CoastSat script in the project
from coastsat_lrr import load_timeseries, compute_lrr

# Optional LOWESS for spatial overlay
_LOWESS_OK = False
if LOWESS_OVERLAY or LOWESS_ONLY:
    try:
        from statsmodels.nonparametric.smoothers_lowess import lowess as sm_lowess
        _LOWESS_OK = True
    except ImportError:
        print("WARNING: statsmodels not found — LOWESS overlay disabled.")


# Color palette builder

_PALETTES = {
    "ocean_to_ember": [
        "#0d3b5e", "#1a6b8a", "#2196b0", "#48cae4",
        "#f4a261", "#e76f51", "#c1440e", "#7b1c00",
    ],
    "indigo_to_coral": [
        "#2c2f6b", "#3a5fc1", "#5e8de8", "#96b8f7",
        "#f7a399", "#f06060", "#c73232", "#8b0000",
    ],
    "navy_to_gold": [
        "#03045e", "#0077b6", "#00b4d8", "#90e0ef",
        "#f4d35e", "#ee964b", "#c77b2a", "#8b5e15",
    ],
}
# -----------------------------------------------------------------------------


# Return exactly n colors interpolated through the chosen palette
def build_color_list(palette_name: str, n: int) -> list:
    if palette_name == "custom":
        return [CUSTOM_COLORS[i % len(CUSTOM_COLORS)] for i in range(n)]

    if palette_name in _PALETTES:
        anchors    = _PALETTES[palette_name]
        rgb_anch   = np.array([mcolors.to_rgb(c) for c in anchors])
        positions  = np.linspace(0, 1, len(anchors))
        targets    = np.linspace(0, 1, n)
        colors = []
        for t in targets:
            idx  = int(np.searchsorted(positions, t, side="right"))
            idx  = np.clip(idx, 1, len(anchors) - 1)
            lo, hi = positions[idx - 1], positions[idx]
            frac = (t - lo) / (hi - lo) if hi > lo else 0.0
            rgb  = (1 - frac) * rgb_anch[idx - 1] + frac * rgb_anch[idx]
            colors.append(tuple(rgb))
        return colors

    try:
        cmap = cm.get_cmap(palette_name)
        return [cmap(i / max(n - 1, 1)) for i in range(n)]
    except ValueError:
        print(f"  WARNING: unknown palette '{palette_name}', falling back to viridis.")
        cmap = cm.get_cmap("viridis")
        return [cmap(i / max(n - 1, 1)) for i in range(n)]


# Data loading

# {transect: CSV path} under the CoastSat root
def collect_csv_map(root_dir: str, site_filter: str = "") -> dict:
    csv_map = {}
    if not os.path.isdir(root_dir):
        print(f"  WARNING: ROOT_DATA_DIR not found: {root_dir}")
        return csv_map
    subfolders = [
        os.path.join(root_dir, d)
        for d in sorted(os.listdir(root_dir))
        if os.path.isdir(os.path.join(root_dir, d))
        and (site_filter == "" or site_filter in d)
    ]
    print(f"Scanning {len(subfolders)} site subfolder(s)...")
    for sf in subfolders:
        for fpath in glob.glob(os.path.join(sf, "*.csv")):
            stem = os.path.splitext(os.path.basename(fpath))[0]
            csv_map[stem] = fpath
    print(f"  -> {len(csv_map)} total time-series CSVs found.\n")
    return csv_map


# Every transect's time series in the lookup
def load_all_transect_data(lookup: pd.DataFrame, csv_map: dict) -> dict:
    all_data  = {}
    n_missing = 0
    for tid in lookup["transect_id"].astype(str).unique():
        if tid not in csv_map:
            n_missing += 1
            continue
        try:
            all_data[tid] = load_timeseries(csv_map[tid])
        except Exception as e:
            print(f"  ERROR loading {tid}: {e}")
            n_missing += 1
    print(f"Loaded {len(all_data)} transect time series "
          f"({n_missing} missing or unreadable).\n")
    return all_data


# Bin generation — whole-year snapping, no partial tail bins

# Divide a period into non-overlapping sequential bins of interval_yr years
def make_bins(period_start: str, period_end: str,
              interval_yr: int, min_bin_fraction: float) -> list:
    start_yr = pd.Timestamp(period_start).year
    end_yr   = pd.Timestamp(period_end).year + 1   # exclusive end

    bins   = []
    cursor = start_yr
    while cursor < end_yr:
        bin_end = cursor + interval_yr
        if bin_end > end_yr:
            bin_end = end_yr   # clamp last bin to period end
        width = bin_end - cursor
        if width >= interval_yr * min_bin_fraction:
            # Label shows the INCLUSIVE year range
            label = f"{cursor}\u2013{bin_end - 1}"
            bins.append((cursor, bin_end, label))
        cursor = cursor + interval_yr   # always advance by full interval

    return bins


# LRR computation per bin — identical to reference script logic

# Build a (bins x domains) LRR matrix
def compute_domain_bin_lrr(
        all_data:        dict,
        lookup:          pd.DataFrame,
        bins:            list,
        period_start:    str,
        period_end:      str,
        min_obs:         int,
        max_pvalue:      float,
        min_r2:          float,
        max_rate:        float,
        agg:             str  = "median",
        buffer_domains:  list = None,
) -> pd.DataFrame:
    if buffer_domains is None:
        buffer_domains = []

    ps_ts = pd.Timestamp(period_start, tz="UTC")
    pe_ts = pd.Timestamp(period_end,   tz="UTC")

    lkp = lookup.copy()
    if buffer_domains:
        lkp = lkp[~lkp["domain_number"].isin(buffer_domains)]
    domains = sorted(lkp["domain_number"].unique())

    # Pre-filter each transect to the period (avoids re-loading per bin)
    period_data = {}
    for tid, df in all_data.items():
        sub = df[(df["date"] >= ps_ts) & (df["date"] <= pe_ts)]
        if len(sub) > 0:
            period_data[tid] = sub

    print(f"    {len(bins)} bins  x  {len(domains)} domains  "
          f"({len(period_data)} transects with data in period)")
    for bs, be, bl in bins:
        print(f"      {bl}  ({bs} -> {be})")

    lrr_matrix = np.full((len(bins), len(domains)), np.nan)

    for i, (bs, be, label) in enumerate(bins):
        bin_start_ts = pd.Timestamp(f"{bs}-01-01", tz="UTC")
        bin_end_ts   = pd.Timestamp(f"{be}-01-01", tz="UTC")

        for j, domain in enumerate(domains):
            tids = lkp[lkp["domain_number"] == domain]["transect_id"].astype(str).values

            valid_rates = []
            for tid in tids:
                if tid not in period_data:
                    continue

                # Slice to this bin's date range
                df_bin = period_data[tid][
                    (period_data[tid]["date"] >= bin_start_ts) &
                    (period_data[tid]["date"] <  bin_end_ts)
                ]

                if len(df_bin) < min_obs:
                    continue

                # Use compute_lrr() — identical to all other scripts
                result = compute_lrr(df_bin)

                if result["lrr_m_yr"] is None or np.isnan(result["lrr_m_yr"]):
                    continue

                # Quality filter 1: statistical significance
                if result["p_value"] > max_pvalue:
                    continue

                # Quality filter 2: goodness of fit
                if result["r_squared"] < min_r2:
                    continue

                # Quality filter 3: physical plausibility
                if abs(result["lrr_m_yr"]) > max_rate:
                    continue

                valid_rates.append(result["lrr_m_yr"])

            if not valid_rates:
                continue

            lrr_matrix[i, j] = (np.median(valid_rates) if agg == "median"
                                  else np.mean(valid_rates))

        n_valid_domains = (~np.isnan(lrr_matrix[i])).sum()
        print(f"      {label}: {n_valid_domains}/{len(domains)} domains with valid LRR")

    bin_labels = [b[2] for b in bins]
    lrr_df     = pd.DataFrame(lrr_matrix, index=bin_labels, columns=domains)
    lrr_df.index.name = "bin"
    return lrr_df


# Single full-period LRR per domain for the bottom summary bar panel
def compute_overall_lrr(
        all_data:       dict,
        lookup:         pd.DataFrame,
        period_start:   str,
        period_end:     str,
        min_obs:        int,
        max_pvalue:     float,
        min_r2:         float,
        max_rate:       float,
        agg:            str  = "median",
        buffer_domains: list = None,
) -> tuple:
    if buffer_domains is None:
        buffer_domains = []

    ps_ts = pd.Timestamp(period_start, tz="UTC")
    pe_ts = pd.Timestamp(period_end,   tz="UTC")

    lkp = lookup.copy()
    if buffer_domains:
        lkp = lkp[~lkp["domain_number"].isin(buffer_domains)]
    domains = sorted(lkp["domain_number"].unique())

    overall = np.full(len(domains), np.nan)
    for j, domain in enumerate(domains):
        tids = lkp[lkp["domain_number"] == domain]["transect_id"].astype(str).values
        rates = []
        for tid in tids:
            if tid not in all_data:
                continue
            df_p = all_data[tid][
                (all_data[tid]["date"] >= ps_ts) &
                (all_data[tid]["date"] <= pe_ts)
            ]
            if len(df_p) < min_obs:
                continue
            result = compute_lrr(df_p)
            if np.isnan(result["lrr_m_yr"]):
                continue
            if result["p_value"]   > max_pvalue:
                continue
            if result["r_squared"] < min_r2:
                continue
            if abs(result["lrr_m_yr"]) > max_rate:
                continue
            rates.append(result["lrr_m_yr"])
        if rates:
            overall[j] = np.median(rates) if agg == "median" else np.mean(rates)

    return np.array(domains), overall


# LOWESS spatial smoothing helper

# Apply LOWESS smoothing along the domain (spatial) axis
def lowess_smooth(x: np.ndarray, y: np.ndarray, frac: float) -> np.ndarray:
    valid = ~np.isnan(y)
    if valid.sum() < 4:
        return y
    smoothed_valid = sm_lowess(y[valid], x[valid],
                                frac=frac, it=0, return_sorted=False)
    out = np.full_like(y, np.nan)
    out[valid] = smoothed_valid
    return out


# Geographic annotation (Section 4b style)

# The full geographic annotation suite on one axes, in the style every CoastSat figure uses
def draw_annotations(ax, domains, show_labels: bool = True):
    d_min, d_max = int(domains.min()), int(domains.max())

    # 1. Wimble Shoals shading (draw first so town spans layer on top)
    ws_lo, ws_hi = ANN_WIMBLE_SHOALS
    lo_c = max(ws_lo, d_min)
    hi_c = min(ws_hi, d_max)
    if lo_c <= hi_c:
        ax.axvspan(lo_c - 0.5, hi_c + 0.5,
                   color=ANN_C_WIMBLE, alpha=0.10, zorder=0)
        if show_labels:
            ax.text(
                (lo_c + hi_c) / 2.0, 1.005, "Wimble Shoals",
                transform=ax.get_xaxis_transform(),
                ha="center", va="bottom",
                fontsize=7, color=ANN_C_WIMBLE, fontstyle="italic",
                clip_on=False,
            )

    # 2. Town span shading + labels (top panel only)
    for name, (d_lo, d_hi) in ANN_TOWN_SPANS.items():
        lo_c = max(d_lo, d_min)
        hi_c = min(d_hi, d_max)
        if lo_c > hi_c:
            continue
        ax.axvspan(lo_c - 0.5, hi_c + 0.5,
                   color=ANN_C_TOWN_SPAN, alpha=0.22, zorder=0)
        if show_labels:
            ax.text(
                (lo_c + hi_c) / 2.0, 1.005, name,
                transform=ax.get_xaxis_transform(),
                ha="center", va="bottom",
                fontsize=8.5, color=ANN_C_TOWN_SPAN,
                fontweight="bold", clip_on=False,
            )

    # 3. Village divider lines inside Tri-Village
    for name, domain in ANN_VILLAGE_LINES.items():
        if d_min <= domain <= d_max:
            ax.axvline(domain, color=ANN_C_VILLAGE_LINE,
                       linewidth=0.8, linestyle=":", zorder=1, alpha=0.75)
            if show_labels:
                ax.text(
                    domain, 0.97, name,
                    transform=ax.get_xaxis_transform(),
                    ha="center", va="top", rotation=90,
                    fontsize=6.5, color=ANN_C_VILLAGE_LINE,
                    clip_on=True,
                )

    # 4. Pier markers
    for name, (domain, label_y) in ANN_PIERS.items():
        if d_min <= domain <= d_max:
            ax.axvline(domain, color=ANN_C_PIER,
                       linewidth=1.3, linestyle="--", zorder=2)
            if show_labels:
                ax.text(
                    domain, label_y, name,
                    transform=ax.get_xaxis_transform(),
                    ha="center", va="center", rotation=90,
                    fontsize=6.5, color=ANN_C_PIER,
                    fontweight="bold", clip_on=True,
                )

    # 5. Groin line
    for name, domain in ANN_GROINS.items():
        if d_min <= float(domain) <= d_max:
            ax.axvline(domain, color=ANN_C_GROIN,
                       linewidth=1.5, linestyle="-.", zorder=2)
            if show_labels:
                ax.text(
                    domain, ANN_GROIN_LABEL_Y, name,
                    transform=ax.get_xaxis_transform(),
                    ha="center", va="center", rotation=90,
                    fontsize=6.5, color=ANN_C_GROIN,
                    fontweight="bold", clip_on=True,
                )


# Single-period figure

# Two-panel figure for one period and interval size
def plot_interval_lines(
        lrr_df:          pd.DataFrame,
        bins:            list,
        overall_domains: np.ndarray,
        overall_lrr:     np.ndarray,
        period_title:    str,
        interval_yr:     int,
        palette:         str,
        line_width:      float,
        line_alpha:      float,
        ylim_lrr,
        ylim_overall,
        out_path:        str,
        dpi:             int = 150,
):
    if lrr_df.empty or np.all(np.isnan(lrr_df.values)):
        print(f"    -> SKIP (no data): {os.path.basename(out_path)}")
        return

    domains = lrr_df.columns.values.astype(int)
    n_bins  = len(bins)
    colors  = build_color_list(palette, n_bins)

    # Y-limit from the true data range, lightly padded
    all_vals = lrr_df.values.flatten()
    valid_vals = all_vals[~np.isnan(all_vals)]
    if ylim_lrr is None:
        if len(valid_vals):
            y_min = float(np.nanmin(valid_vals))
            y_max = float(np.nanmax(valid_vals))
            pad   = max((y_max - y_min) * 0.08, 0.5)
            computed_ylim = (y_min - pad, y_max + pad)
        else:
            computed_ylim = (-5, 5)
    else:
        computed_ylim = ylim_lrr

    fig = plt.figure(figsize=(FIG_WIDTH, FIG_HEIGHT))
    gs  = gridspec.GridSpec(
        2, 1,
        height_ratios=[3, 1],
        hspace=0.18,          # slightly more space so bottom title clears top labels
        left=0.08, right=0.80, top=0.88, bottom=0.09,
    )
    ax_lrr     = fig.add_subplot(gs[0])
    ax_overall = fig.add_subplot(gs[1], sharex=ax_lrr)

    # Top panel: full annotations with labels
    draw_annotations(ax_lrr, domains, show_labels=True)
    # Bottom panel: shading/lines only — no text, so it can't clash with the panel title
    draw_annotations(ax_overall, domains, show_labels=False)

    # Top panel: one line per bin
    legend_handles = []
    n_plotted = 0
    for i, (bs, be, label) in enumerate(bins):
        lrr_vals = lrr_df.iloc[i].values.astype(float)
        n_valid  = (~np.isnan(lrr_vals)).sum()
        if n_valid == 0:
            continue

        color = colors[i]

        # Raw line

        # With LOWESS_ONLY the raw line fades (RAW_ALPHA_SCALE 0 hides it)
        raw_alpha = (line_alpha * RAW_ALPHA_SCALE
                     if (LOWESS_ONLY and _LOWESS_OK) else line_alpha)
        raw_lw    = (line_width * 0.7
                     if (LOWESS_ONLY and _LOWESS_OK) else line_width)
        if raw_alpha > 0.0:
            ax_lrr.plot(
                domains, lrr_vals,
                color=color, linewidth=raw_lw, alpha=raw_alpha,
                zorder=3, solid_capstyle="round",
            )

        # LOWESS smoothed overlay
        if LOWESS_OVERLAY and _LOWESS_OK:
            smoothed = lowess_smooth(domains.astype(float), lrr_vals, LOWESS_FRAC)
            ax_lrr.plot(
                domains, smoothed,
                color=color, linewidth=LOWESS_LINE_WIDTH, alpha=LOWESS_LINE_ALPHA,
                zorder=4, solid_capstyle="round",
            )

        # Legend line weight: use smoothed weight if overlay is active
        leg_lw = (LOWESS_LINE_WIDTH if (LOWESS_OVERLAY and _LOWESS_OK)
                  else line_width + 0.4)
        legend_handles.append(
            mlines.Line2D([], [], color=color, linewidth=leg_lw,
                          label=f"{label}  (n_dom={n_valid})")
        )
        n_plotted += 1

    ax_lrr.axhline(0, color="#333333", linewidth=1.0,
                   linestyle="--", alpha=0.45, zorder=1)

    ax_lrr.set_xlim(domains.min() - 1, domains.max() + 1)
    ax_lrr.set_ylim(computed_ylim)
    ax_lrr.set_ylabel("Shoreline Change Rate — LRR (m/yr)\n"
                       "[+ accretion,  - erosion]", fontsize=10)
    ax_lrr.tick_params(labelsize=9)
    ax_lrr.tick_params(axis="x", labelbottom=False)
    ax_lrr.grid(True, axis="y", alpha=0.22, linewidth=0.7)
    ax_lrr.grid(True, axis="x", alpha=0.10, linewidth=0.5)
    ax_lrr.set_facecolor("#F7F7F7")
    ax_lrr.set_xticks([d for d in domains if d % 10 == 0])

    ax_lrr.legend(
        handles=legend_handles,
        title=f"{interval_yr}-yr bins",
        title_fontsize=8,
        fontsize=8.5,
        loc="upper left",
        bbox_to_anchor=(1.01, 1.0),
        borderaxespad=0,
        framealpha=0.93,
        edgecolor="#CCCCCC",
    )

    # Bottom panel: overall period LRR
    bar_colors = [
        "#CCCCCC" if np.isnan(v) else ("#1a6faf" if v >= 0 else "#c0392b")
        for v in overall_lrr
    ]
    ax_overall.bar(overall_domains, overall_lrr,
                   color=bar_colors, edgecolor="none", width=0.85, zorder=2)
    ax_overall.axhline(0, color="#333333", linewidth=0.9, zorder=3)
    ax_overall.set_title(
        f"Overall LRR — full period (single regression, all obs)",
        fontsize=8, loc="left", pad=3, color="#444444",
        # y slightly below the top of the axes so annotation text above it is clear
    )
    ax_overall.set_ylabel("Overall LRR\n(m/yr)", fontsize=9)
    if ylim_overall is not None:
        ax_overall.set_ylim(ylim_overall)
    ax_overall.tick_params(labelsize=8)
    ax_overall.grid(True, axis="y", alpha=0.22, linewidth=0.7)
    ax_overall.set_facecolor("#F7F7F7")
    ax_overall.set_xlabel("CASCADE Domain  (S -> N)", fontsize=10)
    ax_overall.set_xticks([d for d in domains if d % 10 == 0])

    overlay_note = ""
    smoothing_note = ""
    if LOWESS_OVERLAY and _LOWESS_OK:
        overlay_note = f"  |  LOWESS frac={LOWESS_FRAC:.3f} (7 domains)"
        if LOWESS_ONLY and RAW_ALPHA_SCALE > 0:
            smoothing_note = "  |  thick=smoothed, faint=raw"
        elif LOWESS_ONLY and RAW_ALPHA_SCALE == 0:
            smoothing_note = "  |  LOWESS only (raw hidden)"
        else:
            smoothing_note = "  |  thin=raw, thick=smoothed"
    fig.suptitle(
        f"{period_title}   |   {interval_yr}-Year Intervals  ({n_plotted} bins)"
        f"{overlay_note}{smoothing_note}\n",
        fontsize=11, fontweight="bold", y=0.97,
    )

    plt.savefig(out_path, dpi=dpi, bbox_inches="tight")
    plt.close()
    print(f"    -> Saved: {os.path.relpath(out_path, OUTPUT_DIR)}")


# Run: every window and bin size, the per-domain tables
def main():
    # One folder per window, overwritten in place (deterministic)
    run_dir = OUTPUT_DIR
    os.makedirs(run_dir, exist_ok=True)
    print(f"\nOutput directory: {run_dir}\n{'='*65}\n")

    print(f"Quality filters active:")
    print(f"  Max p-value      : {MAX_PVALUE}")
    print(f"  Min R²           : {MIN_R2}")
    print(f"  Max |LRR|        : {MAX_RATE_M_YR} m/yr")
    print(f"  Min obs per bin  : {MIN_OBS_PER_BIN}")
    print(f"  Min bin fraction : {MIN_BIN_FRACTION}")
    print(f"  Domain agg       : {DOMAIN_AGG}\n")

    print(f"Loading lookup table: {LOOKUP_CSV}")
    lookup = pd.read_csv(LOOKUP_CSV)
    lookup["transect_id"]   = lookup["transect_id"].astype(str)
    lookup["domain_number"] = lookup["domain_number"].astype(int)
    print(f"  {len(lookup)} transect-domain pairs, "
          f"{lookup['domain_number'].nunique()} unique domains.\n")

    csv_map  = collect_csv_map(ROOT_DATA_DIR, SITE_FILTER)
    all_data = load_all_transect_data(lookup, csv_map)

    for interval_yr in INTERVAL_SIZES_YR:

        print(f"\n{'='*65}")
        print(f"Interval: {interval_yr} yr  (non-overlapping, whole-year bins)")
        print(f"{'='*65}")

        for (file_tag, period_title, fig_stem, p_start, p_end) in PERIODS:
            print(f"\n  Period: {period_title}")

            # Generate clean whole-year bins
            bins = make_bins(p_start, p_end, interval_yr, MIN_BIN_FRACTION)
            print(f"    -> {len(bins)} bins after partial-bin filtering")

            # Compute binned LRR matrix
            lrr_df = compute_domain_bin_lrr(
                all_data       = all_data,
                lookup         = lookup,
                bins           = bins,
                period_start   = p_start,
                period_end     = p_end,
                min_obs        = MIN_OBS_PER_BIN,
                max_pvalue     = MAX_PVALUE,
                min_r2         = MIN_R2,
                max_rate       = MAX_RATE_M_YR,
                agg            = DOMAIN_AGG,
                buffer_domains = BUFFER_DOMAINS,
            )

            # Compute overall period LRR for the bottom panel
            ov_domains, ov_lrr = compute_overall_lrr(
                all_data       = all_data,
                lookup         = lookup,
                period_start   = p_start,
                period_end     = p_end,
                min_obs        = MIN_OBS_PER_BIN * 3,   # stricter for full-period fit
                max_pvalue     = MAX_PVALUE,
                min_r2         = MIN_R2,
                max_rate       = MAX_RATE_M_YR,
                agg            = DOMAIN_AGG,
                buffer_domains = BUFFER_DOMAINS,
            )

            # Save CSV: <window>/lrr_bins_<n>yr.csv, rows = bins, cols = GIS
            win_dir = os.path.join(run_dir, file_tag)
            os.makedirs(win_dir, exist_ok=True)
            csv_name = f"lrr_bins_{interval_yr}yr.csv"
            lrr_df.to_csv(os.path.join(win_dir, csv_name))
            print(f"    CSV saved: {file_tag}/{csv_name}")

    print(f"\n{'='*65}")
    print(f"Done. All outputs written to:\n  {run_dir}")
    print(f"{'='*65}")


if __name__ == "__main__":
    main()
