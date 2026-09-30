#!/usr/bin/env python3
"""
Overlay up to four finished CASCADE runs' shoreline change rates against the smoothed CoastSat LRR.

    python scripts/analyze_output/compare_runs/compare_runs.py

Fill RUNS_TO_COMPARE and COMPARISON_NAME first (the list ships empty). Reads
each run's tables/shoreline_change_rate.csv; writes diagnostic, annotated,
two-period and residual figures to output/comparisons/<COMPARISON_NAME>/.
Details: scripts/analyze_output/README.md.

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-09-30
"""
import os
import pathlib
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# House style (site_layer/hat_figure_style.py), applied at import
import sys as _sys
from pathlib import Path as _HP
_sys.path.insert(0, str(next(_q for _q in _HP(__file__).resolve().parents
                             if (_q / "pyproject.toml").exists()) / "scripts"))
from site_layer.hat_figure_style import (apply_style, figsize,  # noqa: E402
                              DOMAIN_AXIS_LABEL, FIG_W_DOUBLE)
apply_style()
import matplotlib.ticker as ticker
import matplotlib.colors as mcolors
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from matplotlib.transforms import blended_transform_factory
from statsmodels.nonparametric.smoothers_lowess import lowess
# Resolved through hat_observed_rates.py (2026-09-18), not typed.
import sys as _sys
from pathlib import Path as _RP
_sys.path.insert(0, str(next(_q for _q in _RP(__file__).resolve().parents
                             if (_q / "pyproject.toml").exists()) / "scripts"))
from site_layer import hat_observed_rates as _obs  # noqa: E402
from site_layer.hat_figure_style import COMPARISONS_ROOT  # noqa: E402


# --- CONFIG ------------------------------------------------------------------
# Domain layout; must match the hindcast runner
NUM_REAL_DOMAINS   = 90
NUM_BUFFER_DOMAINS = 15
FIRST_FILE_NUMBER = 1     # GIS domain IDs: 1–90
LAST_FILE_NUMBER  = FIRST_FILE_NUMBER + NUM_REAL_DOMAINS - 1   # = 90
TOTAL_DOMAINS    = NUM_BUFFER_DOMAINS + NUM_REAL_DOMAINS + NUM_BUFFER_DOMAINS  # 120
START_REAL_INDEX = NUM_BUFFER_DOMAINS             # = 15
END_REAL_INDEX   = START_REAL_INDEX + NUM_REAL_DOMAINS  # = 105
DOMAIN_TICK_STEP    = 5
DOMAIN_SPACING_M    = 500   # metres per CASCADE domain (used to convert window_domains → km)
# Paths resolved from the repo root, never typed
PROJECT_BASE_DIR = next(
    p for p in pathlib.Path(__file__).resolve().parents
    if (p / "pyproject.toml").exists()
)
RAW_RUNS = PROJECT_BASE_DIR / "output" / "raw_runs"
COASTSAT_BASE_DIR = str(_obs.COASTSAT_LRR_ROOT)
# Figures go under output/comparisons/<COMPARISON_NAME>/
COMPARISON_ROOT_DIR = str(COMPARISONS_ROOT)
# Rate CSV columns; lrr_m_yr because the CoastSat target is an OLS rate too
RUN_DOMAIN_COL = "gis_domain"
RUN_RATE_COL   = "lrr_m_yr"
sys.path.insert(0, str(PROJECT_BASE_DIR / "scripts"))
from cascade_pipeline.run_registry import find_run_dir   # noqa: E402
from cascade_pipeline.run_layout import resolve as resolve_run_file  # noqa: E402

# Section 3: runs to compare

# Runs to overlay (fields in README); EMPTY: name live runs before running
RUNS_TO_COMPARE = [
    # dict(
    #     run_name   = "HAT_1984_2004_calibBE_road_bdm_groin",
    #     period     = "1984_2004",
    #     preset     = "calibBE",
    #     label      = "Calibrated (1984)",
    #     start_year = 1984,
    #     sort_key   = 2.5,
    # ),
    # dict(
    #     run_name   = "HAT_2004_2024_calibBE_road_bdm_nourish_groin",
    #     period     = "2004_2024",
    #     preset     = "calibBE",
    #     label      = "Calibrated (2004)",
    #     start_year = 2004,
    #     sort_key   = 3.0,
    # ),
    # A run outside the raw_runs tree entirely -- see run_dir in Section 3.
    # dict(
    #     run_name   = "HAT_2004_2024_L7_Hs2p5",
    #     label      = "Hs2.5 (2004)",
    #     start_year = 2004,
    #     sort_key   = 2.5,
    #     run_dir    = r"D:\shared\HAT_2004_2024_L7_Hs2p5",
    # ),
]

# Name this analysis (controls the comparison folder + filenames)
COMPARISON_NAME = "source_sink_zones"   # <-- EDIT THIS to name your comparison folder

# Section 4: CoastSat datasets

# CoastSat LRR per period; LOWESS at transect resolution, then domain means
COASTSAT_DATASETS = [
    dict(
        label           = "CoastSat (1984–2004)",
        period_start    = 1984,
        csv             = os.path.join(COASTSAT_BASE_DIR, "1984_2004", "transect_lrr_full.csv"),
        domain_col      = "domain_number",
        rate_col        = "lrr_m_yr",
        transect_id_col = "transect_id",
    ),
    dict(
        label           = "CoastSat (2004–2024)",
        period_start    = 2004,
        csv             = os.path.join(
            COASTSAT_BASE_DIR, "2004_2024", "transect_lrr_full.csv"
        ),
        domain_col      = "domain_number",
        rate_col        = "lrr_m_yr",
        transect_id_col = "transect_id",
    ),
]
# Not used anywhere (README)
TRANSECT_DATASETS = [
    dict(
        period_start = 1984,
        csv          = os.path.join(COASTSAT_BASE_DIR, "1984_2004", "transect_lrr_full.csv"),
        domain_col   = "domain_number",   # CASCADE domain ID column (1–90)
        lrr_col      = "lrr_m_yr",        # individual transect LRR column
    ),
    dict(
        period_start = 2004,
        csv          = os.path.join(
            COASTSAT_BASE_DIR, "2004_2024", "transect_lrr_full.csv"
        ),
        domain_col   = "domain_number",
        lrr_col      = "lrr_m_yr",
    ),
]
# CoastSat period drawn solid; None infers it from the runs
ACTIVE_PERIOD_START = None   # 1984 or 2004, or None for auto
# LOWESS widths in domains (1 domain = 500 m); [7, 10] until 2026-09-28
LOWESS_WINDOW_DOMAINS = [7]   # list of 1 or 2 window sizes (domains); [7, 10] until 2026-09-28
# (linewidth, linestyle, active alpha) per width; fill only for '-'
LOWESS_WINDOW_STYLES = [
    (2.0, "-",  1.00),   # the 7-domain window: solid, full opacity, primary reference
]
# The width used as the residuals reference
RESIDUALS_LOWESS_WINDOW = 7

# Section 5: plot options

# Draw the residuals figure
PLOT_RESIDUALS = True
# Draw the two-panel figure (1984-start left, 2004-start right)
PLOT_TWO_PERIOD = True
# Pier and groin label heights, axes fraction 0-1
ANN_PIER_LABEL_Y  = 0.80   # default rotated label y for any pier not given its own override
ANN_PIERS = {
    "Avon Pier":     (26, 0.85),   # (domain, label_y) - label_y is axes-fraction [0,1]
    "Rodanthe Pier": (79, 0.70),
}
ANN_GROIN_LABEL_Y = 0.65
# Accretion / erosion label heights (axes fraction); None = automatic
LABEL_ACCRETION_Y = None   # e.g. 0.80 to pin near the top
LABEL_EROSION_Y   = None   # e.g. 0.15 to pin near the bottom

# Colour palette reference

# Colormap for run colours, sampled light -> dark
RUN_COLORMAP = "YlOrRd"
# Part of the colormap used: skips near-white and near-black
RUN_COLORMAP_RANGE = (0.35, 0.95)
# CoastSat LOWESS colour per width
CS_WINDOW_COLORS = {
     7: "#6BAED6",   # medium sky blue  — 7-domain LOWESS
    10: "#08519C",   # deep ocean blue  — 10-domain LOWESS
}
CS_WINDOW_COLOR_DEFAULT = "#4A7C8E"   # fallback for any unlisted window size
# Transect scatter, styled as in the hindcast script
CS_RAW_COLOR            = "#5BA3C9"    # medium blue
PLOT_RAW_LRR            = True         # set False to hide transect scatter from all figures
RAW_LRR_SOUTHERN_ONLY   = True         # True: scatter only where LOWESS is suppressed
RAW_LRR_SCATTER_SIZE    = 6            # marker area in points²
RAW_LRR_SCATTER_ALPHA   = 0.60         # opacity for active period; ×0.35 for reference period
# Geographic annotations, shared with the hindcast script
ANN_TOWN_SPANS = {
    "Buxton":      (7,   8),
    "Avon":        (21, 31),
    "Tri-Village": (68, 83),
}
ANN_VILLAGE_LINES = {"Salvo": 69, "Waves": 74, "Rodanthe": 80}
ANN_GROINS        = {"Buxton Groin": 5.5}   # boundary between domains 5 and 6
ANN_WIMBLE_SHOALS = (60, 74)
ANN_AVON_SHOALS   = (24, 39)   # Avon Shoals influence zone (same feature type as Wimble Shoals)
# Annotation colours
ANN_C_TOWN_SPAN    = "#90AFC5"
ANN_C_WIMBLE       = "#E0A800"   # amber - both shoal zones share this color
ANN_C_AVON_SHOALS  = "#E0A800"   # same amber as Wimble Shoals (same feature type)
ANN_C_VILLAGE_LINE = "0.40"
ANN_C_PIER         = "#1565C0"
ANN_C_GROIN        = "#B71C1C"
# Domains 1..N show raw scatter instead of LOWESS (Oregon Inlet)
LOWESS_SKIP_SOUTHERN_DOMAINS = 10
# -----------------------------------------------------------------------------


# GIS domain ID (1-based) -> CASCADE padded array index
def _gis_to_pad(gis_id):
    return START_REAL_INDEX + (gis_id - FIRST_FILE_NUMBER)


# Light -> dark colours by sort_key (or list order); an explicit 'color' is kept
def assign_run_colors(runs):
    cmap = plt.get_cmap(RUN_COLORMAP)
    lo, hi = RUN_COLORMAP_RANGE

    # Split into runs that need an auto color vs. runs with an explicit override
    auto_runs   = [r for r in runs if r.get("color") is None]
    fixed_runs  = [r for r in runs if r.get("color") is not None]

    if auto_runs:
        have_all_keys = all(r.get("sort_key") is not None for r in auto_runs)
        if have_all_keys:
            ranked = sorted(auto_runs, key=lambda r: r["sort_key"])
        else:
            ranked = auto_runs   # fall back to list order

        n = len(ranked)
        for i, r in enumerate(ranked):
            # n==1 -> sample at the dark end so a single run isn't washed out
            frac = (lo + hi) / 2.0 if n == 1 else lo + (hi - lo) * (i / (n - 1))
            r["color"] = mcolors.to_hex(cmap(frac))

    return fixed_runs + auto_runs if fixed_runs else auto_runs


# A run's (gis_ids, LRR m/yr, run_dir), resolved by run_registry unless run_dir is given
def load_run_rates(run_name, period=None, preset=None, arm=None, run_dir=None):
    # Resolve the run folder through run_registry, never by joining a path
    if run_dir is None:
        if period is None or preset is None:
            raise ValueError(
                f"run '{run_name}' names neither (period, preset) nor an "
                f"explicit run_dir; one of the two is required."
            )
        kwargs = {"arm": arm} if arm else {}
        run_dir = str(find_run_dir(RAW_RUNS, run_name, period, preset, **kwargs))
    # run_layout finds the rate CSV in either the current or the pre-2026-09-10 layout
    csv_path = str(resolve_run_file(run_dir, "rate_csv", run_name))

    if not os.path.exists(csv_path):
        raise FileNotFoundError(
            f"Rate CSV not found for run '{run_name}'.\n"
            f"Expected: {csv_path}\n"
            f"Run HAT_hindcast_1984_2024_old version.py first to generate this file."
        )

    df = pd.read_csv(csv_path)
    required = {RUN_DOMAIN_COL, RUN_RATE_COL}
    if not required.issubset(df.columns):
        raise ValueError(
            f"Rate CSV for '{run_name}' is missing required columns.\n"
            f"Need: {sorted(required)}  |  Found: {sorted(df.columns)}\n"
            f"  {csv_path}"
        )

    # Keyed on the file's own gis_domain column, never on row order
    real = df[df[RUN_DOMAIN_COL].between(FIRST_FILE_NUMBER, LAST_FILE_NUMBER)]
    real = real.sort_values(RUN_DOMAIN_COL)

    if len(real) != NUM_REAL_DOMAINS:
        print(f"  ⚠️  '{run_name}': expected {NUM_REAL_DOMAINS} real-domain rows, "
              f"got {len(real)} — results may be incomplete.")

    gis_ids   = real[RUN_DOMAIN_COL].values.astype(int)
    rates_myr = real[RUN_RATE_COL].values.astype(float)

    return gis_ids, rates_myr, run_dir


# Median spacing between consecutive transects (m)
def estimate_transect_spacing(along_coast_m):
    arr   = np.sort(along_coast_m)
    diffs = np.diff(arr)
    pos   = diffs[diffs > 0]
    return float(np.median(pos)) if len(pos) else 50.0


# Transect LRR with along-coast distance, each domain's transects spread over its 500 m
def load_transect_data(ds):
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

    # Spread each domain's transects over its 500 m band, for the LOWESS frac
    def _spread(grp):
        n         = len(grp)
        base      = (grp[domain_col].iloc[0] - 1) * DOMAIN_SPACING_M
        offsets   = (np.arange(n) + 0.5) * (DOMAIN_SPACING_M / n)
        grp       = grp.copy()
        grp["along_coast_m"] = base + offsets
        return grp

    df = df.groupby(domain_col, group_keys=False).apply(_spread)

    domain_ids    = df[domain_col].values.astype(int)
    lrr_values    = df[rate_col].values.astype(float)
    along_coast_m = df["along_coast_m"].values.astype(float)

    spacing = estimate_transect_spacing(along_coast_m)
    print(f"  ✓ {ds['label']}: {len(df)} transects  "
          f"est. spacing {spacing:.0f} m  "
          f"LRR range {np.nanmin(lrr_values):+.2f}–{np.nanmax(lrr_values):+.2f} m/yr")
    return domain_ids, lrr_values, along_coast_m


# LOWESS at transect resolution, then averaged to domains; returns (gis_x, smoothed, frac)
def lowess_smooth_transect_to_domains(along_coast_m, lrr, domain_ids, window_domains):
    window_km = window_domains * DOMAIN_SPACING_M / 1000.0
    spacing_m = estimate_transect_spacing(along_coast_m)
    n         = len(along_coast_m)
    frac      = float(np.clip((window_km * 1000.0 / spacing_m) / n, 0.02, 1.0))

    valid = np.isfinite(lrr)
    if valid.sum() < 5:
        print(f"  ⚠️  Too few valid transects ({valid.sum()}) for LOWESS — skipping")
        return None, None, frac

    result            = lowess(lrr[valid], along_coast_m[valid], frac=frac, return_sorted=True)
    smoothed_t        = np.full(n, np.nan)
    smoothed_t[valid] = np.interp(along_coast_m[valid], result[:, 0], result[:, 1])

    # Average smoothed transect values within each CASCADE domain
    dom_agg = (pd.DataFrame({"domain": domain_ids, "smoothed": smoothed_t})
                 .groupby("domain")["smoothed"].mean()
                 .dropna())

    return dom_agg.index.values.astype(int), dom_agg.values, frac


# Drop the LOWESS curve over domains 1..skip_n, where raw scatter is shown instead
def splice_lowess_with_raw_south(win_gis_x, win_smoothed, skip_n=None):
    if skip_n is None:
        skip_n = LOWESS_SKIP_SOUTHERN_DOMAINS
    if skip_n <= 0:
        return win_gis_x, win_smoothed
    mask = win_gis_x > skip_n
    return win_gis_x[mask], win_smoothed[mask]


# Shoals, villages, piers and groin on an axis in GIS domain units
def add_geographic_annotations(ax):
    trans = blended_transform_factory(ax.transData, ax.transAxes)

    # 1a. Avon Shoals influence zone (same feature type as Wimble Shoals)
    alo, ahi = ANN_AVON_SHOALS
    ax.axvspan(alo - 0.5, ahi + 0.5,
               color=ANN_C_AVON_SHOALS, alpha=0.10, zorder=0,
               hatch="///", edgecolor=ANN_C_AVON_SHOALS, linewidth=0)
    ax.text((alo + ahi) / 2.0, 0.04,
            "Avon Shoals\nPosition", transform=trans,
            ha="center", va="bottom", fontsize=7, color="#7A5800", style="italic",
            bbox=dict(boxstyle="round,pad=0.2", fc="white", ec="none", alpha=0.80))

    # 1b. Wimble Shoals influence zone
    wlo, whi = ANN_WIMBLE_SHOALS
    ax.axvspan(wlo - 0.5, whi + 0.5,
               color=ANN_C_WIMBLE, alpha=0.10, zorder=0,
               hatch="///", edgecolor=ANN_C_WIMBLE, linewidth=0)
    ax.text((wlo + whi) / 2.0, 0.04,
            "Wimble Shoals\nPosition", transform=trans,
            ha="center", va="bottom", fontsize=7, color="#7A5800", style="italic",
            bbox=dict(boxstyle="round,pad=0.2", fc="white", ec="none", alpha=0.80))

    # 2. Community spans
    for span_label, (d_lo, d_hi) in ANN_TOWN_SPANS.items():
        ax.axvspan(d_lo - 0.5, d_hi + 0.5,
                   color=ANN_C_TOWN_SPAN, alpha=0.14, zorder=0)
        ax.text((d_lo + d_hi) / 2.0, 0.90,
                span_label, transform=trans,
                ha="center", va="top", fontsize=8, color="0.25", fontweight="bold",
                bbox=dict(boxstyle="round,pad=0.2", fc="white", ec="none", alpha=0.85))

    # 3. Village center lines
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


# Legend handles for the annotation layer
def annotation_legend_handles():
    return [
        Patch(fc=ANN_C_TOWN_SPAN, alpha=0.30, label="Community"),
        Patch(fc=ANN_C_WIMBLE, alpha=0.25, hatch="///",
              edgecolor=ANN_C_WIMBLE, linewidth=0, label="Shoals position (Avon / Wimble)"),
        Line2D([0], [0], color=ANN_C_VILLAGE_LINE, lw=0.9, ls="--", label="Village center"),
        Line2D([0], [0], color=ANN_C_PIER,         lw=1.0, ls="-.", label="Pier"),
        Line2D([0], [0], color=ANN_C_GROIN,        lw=1.1, ls=":",  label="Groin"),
    ]


# Shared axis styling
def _style_ax(ax, ylabel="Shoreline change rate (m/yr)"):
    ax.set_xlim(FIRST_FILE_NUMBER - 0.5, LAST_FILE_NUMBER + 0.5)
    ax.axhline(0.0, color="#2c2c2c", linewidth=1.0, linestyle="--", alpha=0.55, zorder=3)
    ax.xaxis.set_major_locator(ticker.MultipleLocator(10))
    ax.xaxis.set_minor_locator(ticker.MultipleLocator(5))
    ax.tick_params(axis="both", which="major", labelsize=10, direction="in", length=5)
    ax.tick_params(axis="both", which="minor", direction="in", length=3)
    ax.grid(True, which="major", linestyle=":", linewidth=0.6, alpha=0.35, color="gray")
    ax.spines[["top", "right"]].set_visible(False)
    ax.spines[["left", "bottom"]].set_linewidth(1.1)
    ax.set_ylabel(ylabel, fontsize=11, fontweight="bold", labelpad=8)


# Plotting

# Quick multi-run diagnostic figure
def plot_diagnostic(run_data, cs_series, active_period, out_path, comparison_name):
    fig, ax = plt.subplots(figsize=figsize("double", height=3.49))
    fig.patch.set_facecolor("white")
    ax.set_facecolor("white")

    model_handles, cs_handles = _draw_comparison_panel(ax, run_data, cs_series, active_period)

    ax.set_xlabel(DOMAIN_AXIS_LABEL,
                  fontsize=11, fontweight="bold", labelpad=4)
    ax.set_title(
        f"CASCADE Run Comparison — Hatteras Island, NC  |  {comparison_name}",
        fontsize=12, fontweight="bold", pad=12, color="#1a2a3a"
    )

    # Reserve the bottom margin before the legend, so it stays on the canvas (README)
    fig.subplots_adjust(bottom=0.26)

    all_handles = model_handles + cs_handles + annotation_legend_handles()
    fig.legend(handles=all_handles,
               loc="lower center",
               bbox_to_anchor=(0.5, 0.04),
               fontsize=9, framealpha=0.95, edgecolor="#cccccc",
               frameon=True, ncol=4)

    fig.savefig(out_path, dpi=300, facecolor="white")
    plt.close(fig)
    print(f"✓ Saved diagnostic:  {os.path.basename(out_path)}")


# One panel: CoastSat scatter and LOWESS, then every run
def _draw_comparison_panel(ax, run_data, cs_series, active_period,
                            panel_title=None, show_reference_period=True):
    add_geographic_annotations(ax)

    cs_handles = []
    widest_window = max(LOWESS_WINDOW_DOMAINS)   # fill drawn only for the widest window
    for cs in cs_series:
        is_active = cs["period_start"] == active_period
        if not is_active and not show_reference_period:
            continue   # skip the inactive period entirely for this panel
        # --- individual transect scatter (lowest zorder — contextual background) ---
        scatter_x = cs["transect_along_coast"] / DOMAIN_SPACING_M + FIRST_FILE_NUMBER
        if PLOT_RAW_LRR:
            if RAW_LRR_SOUTHERN_ONLY:
                south_mask     = cs["transect_domains"] <= LOWESS_SKIP_SOUTHERN_DOMAINS
                scatter_x_plot = scatter_x[south_mask]
                scatter_y_plot = cs["transect_rates"][south_mask]
                raw_lbl = (f"{cs['label']} — transect LRR (D1-{LOWESS_SKIP_SOUTHERN_DOMAINS})"
                           if is_active else None)
            else:
                scatter_x_plot = scatter_x
                scatter_y_plot = cs["transect_rates"]
                raw_lbl = f"{cs['label']} — transect LRR" if is_active else None
            raw_alpha = RAW_LRR_SCATTER_ALPHA if is_active else RAW_LRR_SCATTER_ALPHA * 0.35
            ax.scatter(scatter_x_plot, scatter_y_plot,
                       color=CS_RAW_COLOR, s=RAW_LRR_SCATTER_SIZE,
                       alpha=raw_alpha, zorder=1, linewidths=0, label=raw_lbl)
            if is_active:
                cs_handles.append(
                    Line2D([0], [0], color=CS_RAW_COLOR, marker=".", ms=5,
                           ls="none", alpha=RAW_LRR_SCATTER_ALPHA, label=raw_lbl)
                )
        # --- LOWESS smoothed curves (transect-based, aggregated to domain resolution) ---
        for idx, win in enumerate(cs["windows"]):
            cs_color  = CS_WINDOW_COLORS.get(win["window"], CS_WINDOW_COLOR_DEFAULT)
            lw_base, ls, alpha_factor = (
                LOWESS_WINDOW_STYLES[idx] if idx < len(LOWESS_WINDOW_STYLES)
                else (1.5, "--", 0.80)
            )
            # With RAW_LRR_SOUTHERN_ONLY, LOWESS starts north of the raw-scatter zone
            if RAW_LRR_SOUTHERN_ONLY:
                w_gis_x, rate = splice_lowess_with_raw_south(win["gis_x"], win["smoothed"])
            else:
                w_gis_x = win["gis_x"]
                rate     = win["smoothed"]
            w_lbl    = f"LOWESS {win['window']}-dom"
            lbl      = f"{cs['label']} — {w_lbl}"
            if is_active:
                if win["window"] == widest_window:
                    ax.fill_between(w_gis_x, rate, 0,
                                    where=(rate < 0),  alpha=0.14, color=cs_color,
                                    interpolate=True)
                    ax.fill_between(w_gis_x, rate, 0,
                                    where=(rate >= 0), alpha=0.10, color=cs_color,
                                    interpolate=True)
                ax.plot(w_gis_x, rate, color=cs_color, linewidth=lw_base,
                        linestyle=ls, alpha=alpha_factor, zorder=4, label=lbl)
                cs_handles.append(
                    Line2D([0], [0], color=cs_color, lw=lw_base, ls=ls,
                           alpha=alpha_factor, label=lbl)
                )
            else:
                ax.plot(w_gis_x, rate, color=cs_color, linewidth=lw_base * 0.85,
                        linestyle=ls, alpha=0.40 * alpha_factor, zorder=3)
                cs_handles.append(
                    Line2D([0], [0], color=cs_color, lw=lw_base * 0.85, ls=ls,
                           alpha=0.40 * alpha_factor, label=lbl + " (ref)")
                )

    model_handles = []
    for run in run_data:
        ax.plot(run["gis_ids"], run["rates"],
                color=run["color"], linewidth=2.4, zorder=5, label=run["label"])
        model_handles.append(
            Line2D([0], [0], color=run["color"], lw=2.4, label=run["label"])
        )

    _style_ax(ax)

    # Compass / orientation labels
    ax.text(0.0, 1.01, "← S  |  Cape Point",
            transform=ax.transAxes, fontsize=9, color="#444444",
            ha="left", va="bottom", style="italic", clip_on=False)
    ax.text(1.0, 1.01, "Pea Island  |  N →",
            transform=ax.transAxes, fontsize=9, color="#444444",
            ha="right", va="bottom", style="italic", clip_on=False)

    if panel_title:
        ax.set_title(panel_title, fontsize=11, fontweight="bold", pad=10, color="#1a2a3a")

    return model_handles, cs_handles


# Publication figure with the geographic annotations
def plot_annotated(run_data, cs_series, active_period, out_path, comparison_name):
    fig, ax = plt.subplots(figsize=figsize("double", height=4.01))
    fig.patch.set_facecolor("white")
    ax.set_facecolor("white")

    model_handles, cs_handles = _draw_comparison_panel(ax, run_data, cs_series, active_period)

    ax.set_xlabel(DOMAIN_AXIS_LABEL,
                  fontsize=12, fontweight="bold", labelpad=4)

    # Accretion / erosion side labels
    all_vals = np.concatenate(
        [r["rates"] for r in run_data] +
        [w["smoothed"][np.isfinite(w["smoothed"])]
         for cs in cs_series for w in cs["windows"]]
    )
    ymin = all_vals.min(); ymax = all_vals.max()
    ypad = (ymax - ymin) * 0.07
    ax.set_ylim(ymin - ypad, ymax + ypad)
    ybot, ytop = ax.get_ylim()
    zero_frac  = (0 - ybot) / (ytop - ybot)
    acc_y = LABEL_ACCRETION_Y if LABEL_ACCRETION_Y is not None else zero_frac + (1 - zero_frac) / 2
    ero_y = LABEL_EROSION_Y   if LABEL_EROSION_Y   is not None else zero_frac / 2
    ax.text(1.0, acc_y, "Accretion ▲",
            transform=ax.transAxes, fontsize=9, color="#555555",
            ha="right", va="center", style="italic")
    ax.text(1.0, ero_y, "Erosion ▼",
            transform=ax.transAxes, fontsize=9, color="#555555",
            ha="right", va="center", style="italic")

    ax.set_title(
        f"CASCADE Run Comparison — Hatteras Island, NC  |  {comparison_name}",
        fontsize=12, fontweight="bold", pad=12, color="#1a2a3a"
    )

    # Reserve the bottom margin before legend and caption: off-canvas placement crashed saves (README)
    fig.subplots_adjust(bottom=0.26)

    # Legend: model runs | CoastSat | geographic annotations
    all_handles = model_handles + cs_handles + annotation_legend_handles()
    fig.legend(handles=all_handles,
               loc="lower center",
               bbox_to_anchor=(0.5, 0.07),
               fontsize=9, framealpha=0.95, edgecolor="#cccccc",
               frameon=True, ncol=4)

    # Caption — figure-fraction, sits below the legend, fully on-canvas.
    fig.text(
        0.5, 0.01,
        f"Model: CASCADE  |  Observed: CoastSat LRR "
        f"(LOWESS {'/'.join(str(w) for w in LOWESS_WINDOW_DOMAINS)}-domain windows)  |  "
        f"Comparison: {comparison_name}",
        fontsize=7.5, color="#666666", ha="center", va="bottom", style="italic",
    )

    fig.savefig(out_path, dpi=300, facecolor="white")
    plt.close(fig)
    print(f"✓ Saved annotated:   {os.path.basename(out_path)}")


# 1984-start runs left, 2004-start right, shared y axis
def plot_two_period(run_data, cs_series, out_path, comparison_name):
    runs_1984 = [r for r in run_data if r["start_year"] == 1984]
    runs_2004 = [r for r in run_data if r["start_year"] == 2004]

    panels = []
    if runs_1984:
        panels.append((1984, runs_1984, "1984-2004"))
    if runs_2004:
        panels.append((2004, runs_2004, "2004-2024"))

    if not panels:
        print("  ⚠️  No runs with a recognized start_year — skipping two-period figure.")
        return

    n_panels = len(panels)
    fig, axes = plt.subplots(
        1, n_panels, figsize=figsize("double", height=8 * FIG_W_DOUBLE / (11 * n_panels)), sharey=True,
    )
    if n_panels == 1:
        axes = [axes]
    fig.patch.set_facecolor("white")

    all_model_handles = []
    all_cs_handles_by_period = {}
    all_vals_chunks = []

    for ax, (period_start, panel_runs, period_label) in zip(axes, panels):
        ax.set_facecolor("white")
        model_handles, cs_handles = _draw_comparison_panel(
            ax, panel_runs, cs_series, period_start,
            panel_title=f"{period_label}",
            show_reference_period=False,   # each panel shows ONLY its own period's CoastSat
        )
        all_model_handles.extend(model_handles)
        # CoastSat handles kept per period, so panel 2's aren't dropped as duplicates
        all_cs_handles_by_period[period_start] = cs_handles
        all_vals_chunks.append(np.concatenate(
            [r["rates"] for r in panel_runs] +
            [w["smoothed"][np.isfinite(w["smoothed"])]
             for cs in cs_series if cs["period_start"] == period_start
             for w in cs["windows"]]
        ) if panel_runs else np.array([]))
        ax.set_xlabel(DOMAIN_AXIS_LABEL,
                      fontsize=11, fontweight="bold", labelpad=4)

    # Shared y-limits across both panels, computed from ALL data in either panel
    all_vals = np.concatenate([c for c in all_vals_chunks if c.size])
    ymin = all_vals.min(); ymax = all_vals.max()
    ypad = (ymax - ymin) * 0.07
    for ax in axes:
        ax.set_ylim(ymin - ypad, ymax + ypad)
        ybot, ytop = ax.get_ylim()
        zero_frac  = (0 - ybot) / (ytop - ybot)
        acc_y = LABEL_ACCRETION_Y if LABEL_ACCRETION_Y is not None else zero_frac + (1 - zero_frac) / 2
        ero_y = LABEL_EROSION_Y   if LABEL_EROSION_Y   is not None else zero_frac / 2
        ax.text(1.0, acc_y, "Accretion ▲",
                transform=ax.transAxes, fontsize=8.5, color="#555555",
                ha="right", va="center", style="italic")
        ax.text(1.0, ero_y, "Erosion ▼",
                transform=ax.transAxes, fontsize=8.5, color="#555555",
                ha="right", va="center", style="italic")

    # Only the leftmost panel gets the y-axis label (sharey hides the rest)
    axes[0].set_ylabel("Shoreline change rate (m/yr)", fontsize=11,
                        fontweight="bold", labelpad=8)

    # Reserve top and bottom margins first so everything stays on the canvas (README)
    fig.subplots_adjust(top=0.88, bottom=0.20, wspace=0.06)

    fig.suptitle(
        f"CASCADE Run Comparison by Period — Hatteras Island, NC  |  {comparison_name}",
        fontsize=13, fontweight="bold", color="#1a2a3a", y=0.97,
    )

    # Combined legend: runs (deduplicated by label), CoastSat entries, annotations
    seen_labels = set()
    dedup_model_handles = []
    for h in all_model_handles:
        if h.get_label() not in seen_labels:
            dedup_model_handles.append(h)
            seen_labels.add(h.get_label())

    combined_cs_handles = []
    for period_start, _, period_label in panels:
        combined_cs_handles.extend(all_cs_handles_by_period.get(period_start, []))

    # One CoastSat legend entry per visual style, period dropped from the label (README)
    def _style_key(h):
        return (h.get_color(), h.get_linestyle(), round(h.get_linewidth(), 2),
                h.get_marker(), round(h.get_alpha() or 1.0, 2))

    seen_styles = set()
    dedup_cs_handles = []
    for h in combined_cs_handles:
        key = _style_key(h)
        if key in seen_styles:
            continue
        seen_styles.add(key)
        # Strip the "CoastSat (1984-2004) - " prefix
        label = h.get_label()
        if " — " in label:
            label = label.split(" — ", 1)[1]
        h.set_label(label)
        dedup_cs_handles.append(h)

    all_handles = dedup_model_handles + dedup_cs_handles + annotation_legend_handles()

    fig.legend(handles=all_handles,
               loc="lower center",
               bbox_to_anchor=(0.5, 0.055),
               fontsize=8.5, framealpha=0.95, edgecolor="#cccccc",
               frameon=True, ncol=4)

    fig.text(
        0.5, 0.005,
        f"Model: CASCADE  |  Observed: CoastSat LRR "
        f"(LOWESS {'/'.join(str(w) for w in LOWESS_WINDOW_DOMAINS)}-domain windows)  |  "
        f"Comparison: {comparison_name}  |  Left/right panels each show their own active period",
        fontsize=7.5, color="#666666", ha="center", va="bottom", style="italic",
    )

    # The figure's own bounds, not "tight": margins are already reserved
    fig.savefig(out_path, dpi=300, facecolor="white")
    plt.close(fig)
    print(f"✓ Saved two-period:  {os.path.basename(out_path)}")


# Each run minus the active CoastSat LOWESS
def plot_residuals(run_data, cs_series, active_period, out_path, comparison_name):
    # Find the active CoastSat series and the designated residuals window
    active_cs = next(
        (cs for cs in cs_series if cs["period_start"] == active_period), None
    )
    if active_cs is None:
        print("  ⚠️  No active CoastSat series found — skipping residuals plot.")
        return

    # Select the configured residuals window; fall back to the last (widest) if missing
    active_win = next(
        (w for w in active_cs["windows"] if w["window"] == RESIDUALS_LOWESS_WINDOW),
        active_cs["windows"][-1],
    )
    cs_gis_x    = active_win["gis_x"]
    cs_smoothed = active_win["smoothed"]
    all_gis_ids = np.arange(FIRST_FILE_NUMBER, LAST_FILE_NUMBER + 1, dtype=int)

    # Interpolate CoastSat to full 1–90 grid (fills any gaps)
    cs_interp = np.interp(
        all_gis_ids.astype(float),
        cs_gis_x.astype(float),
        cs_smoothed,
    )

    fig, ax = plt.subplots(figsize=figsize("double", height=2.99))
    fig.patch.set_facecolor("white")
    ax.set_facecolor("white")
    add_geographic_annotations(ax)

    # Per-run fit statistics, shown in the legend
    run_stats = []
    for run in run_data:
        residual = run["rates"] - cs_interp
        mean_abs = np.nanmean(np.abs(residual))
        rmse     = np.sqrt(np.nanmean(residual ** 2))
        run_stats.append((run, residual, mean_abs, rmse))
        ax.plot(run["gis_ids"], residual,
                color=run["color"], linewidth=2.0, zorder=4,
                label=f"{run['label']}  (MAE={mean_abs:.2f}, RMSE={rmse:.2f})")
        ax.fill_between(run["gis_ids"], residual, 0,
                        where=(residual > 0), alpha=0.08, color=run["color"])
        ax.fill_between(run["gis_ids"], residual, 0,
                        where=(residual < 0), alpha=0.08, color=run["color"])

    _style_ax(ax, ylabel="Model − CoastSat (m/yr)")
    ax.set_xlabel(DOMAIN_AXIS_LABEL,
                  fontsize=11, fontweight="bold", labelpad=4)
    ax.set_title(
        f"Residuals: Model − CoastSat (active period)  |  {comparison_name}",
        fontsize=12, fontweight="bold", pad=12, color="#1a2a3a"
    )

    run_handles = [Line2D([0], [0], color=r["color"], lw=2.0,
                          label=f"{r['label']}  (MAE={mae:.2f}, RMSE={rmse:.2f})")
                   for r, _, mae, rmse in run_stats]
    cs_label_handle = Line2D([0], [0], color="gray", lw=1.0, ls="--",
                              label=f"Reference: {active_cs['label']} LOWESS {active_win['window']}-dom")

    # Reserve the bottom margin before the legend (loc="best" landed on the data)
    fig.subplots_adjust(bottom=0.24)
    fig.legend(handles=run_handles + [cs_label_handle],
               loc="lower center",
               bbox_to_anchor=(0.5, 0.03),
               fontsize=9, framealpha=0.95, edgecolor="#cccccc",
               frameon=True, ncol=2)

    fig.savefig(out_path, dpi=300, facecolor="white")
    plt.close(fig)
    print(f"✓ Saved residuals:   {os.path.basename(out_path)}")
    for run, _, mae, rmse in run_stats:
        print(f"    {run['label']:<20s} MAE={mae:.3f} m/yr   RMSE={rmse:.3f} m/yr")


# Run: load runs and CoastSat, smooth, draw every figure
def main():
    if not RUNS_TO_COMPARE:
        raise SystemExit("RUNS_TO_COMPARE is empty: name the runs to compare in the CONFIG "
                         "block first (fields in scripts/analyze_output/README.md)")
    # Resolve comparison name
    global COMPARISON_NAME
    if COMPARISON_NAME is None:
        COMPARISON_NAME = "_vs_".join(r["label"].replace(" ", "") for r in RUNS_TO_COMPARE)
    comparison_name = COMPARISON_NAME

    # Resolve run colors (auto-gradient unless a run overrides 'color')
    RUNS_TO_COMPARE[:] = assign_run_colors(RUNS_TO_COMPARE)

    # Resolve active CoastSat period
    global ACTIVE_PERIOD_START
    if ACTIVE_PERIOD_START is None:
        years = [r["start_year"] for r in RUNS_TO_COMPARE]
        ACTIVE_PERIOD_START = max(set(years), key=years.count)
        print(f"ACTIVE_PERIOD_START inferred from runs: {ACTIVE_PERIOD_START}")
    active_period = ACTIVE_PERIOD_START

    # Create comparison folder
    out_dir = os.path.join(COMPARISON_ROOT_DIR, comparison_name)
    os.makedirs(out_dir, exist_ok=True)
    print(f"\nComparison:   {comparison_name}")
    print(f"Output dir:   {out_dir}")
    print("=" * 70)

    # Load model runs
    print("\nLoading CASCADE run rate CSVs...")
    run_data = []
    for run_cfg in RUNS_TO_COMPARE:
        try:
            gis_ids, rates, run_dir = load_run_rates(
                run_cfg["run_name"],
                period  = run_cfg.get("period"),
                preset  = run_cfg.get("preset"),
                arm     = run_cfg.get("arm"),
                run_dir = run_cfg.get("run_dir"),
            )
            run_data.append(dict(
                run_name   = run_cfg["run_name"],
                label      = run_cfg["label"],
                color      = run_cfg["color"],
                start_year = run_cfg["start_year"],
                gis_ids    = gis_ids,
                rates      = rates,
                run_dir    = run_dir,
            ))
            print(f"  ✓ {run_cfg['label']:<20s}  "
                  f"rate range {rates.min():.2f}–{rates.max():.2f} m/yr  "
                  f"color={run_cfg['color']}  "
                  f"({run_dir})")
        except FileNotFoundError as e:
            print(f"  ❌ SKIPPED '{run_cfg['run_name']}': {e}")

    if not run_data:
        print("\n❌ No valid runs loaded — RUNS_TO_COMPARE is empty "
              "or every entry failed to resolve. See Section 3.")
        sys.exit(1)

    print(f"\n  {len(run_data)} run(s) loaded successfully.")

    # Load CoastSat transects + apply LOWESS at transect resolution
    print("\nLoading CoastSat transect data...")
    cs_series = []
    for ds in COASTSAT_DATASETS:
        domain_ids, lrr_values, along_coast_m = load_transect_data(ds)
        if domain_ids is None:
            continue
        windows = []
        for w in LOWESS_WINDOW_DOMAINS:
            gis_x, smoothed, frac = lowess_smooth_transect_to_domains(
                along_coast_m, lrr_values, domain_ids, w
            )
            if gis_x is None:
                continue
            print(f"  ✓ LOWESS applied: window={w} domains "
                  f"({w * DOMAIN_SPACING_M / 1000.0:.1f} km)  "
                  f"frac={frac:.3f}  ({ds['label']})")
            windows.append(dict(window=w, gis_x=gis_x, smoothed=smoothed, frac=frac))
        cs_series.append(dict(
            label                 = ds["label"],
            period_start          = ds["period_start"],
            transect_domains      = domain_ids,       # one entry per transect (domain ID)
            transect_rates        = lrr_values,        # one entry per transect (LRR y)
            transect_along_coast  = along_coast_m,    # one entry per transect (physical x)
            windows               = windows,           # LOWESS curves at domain resolution
        ))

    if not cs_series:
        print("  ⚠️  No CoastSat data loaded — plots will show model lines only.")
    else:
        # Stop if a run's CoastSat period failed to load, rather than compare against the wrong one
        loaded_periods = {cs["period_start"] for cs in cs_series}
        needed_periods = {r["start_year"] for r in RUNS_TO_COMPARE}
        missing_periods = needed_periods - loaded_periods
        if missing_periods:
            print("=" * 70)
            print(f"WARNING: CoastSat data for period(s) {sorted(missing_periods)} "
                  f"did NOT load (see 'Transect CSV not found' warning above for "
                  f"the exact path checked).")
            print(f"         Runs with start_year in {sorted(missing_periods)} will "
                  f"be compared against the WRONG period's CoastSat data, or show "
                  f"no CoastSat overlay at all, in figures that need it.")
            print(f"         Fix the path in COASTSAT_DATASETS (Section 4) before "
                  f"trusting any figure that includes these runs.")
            print("=" * 70)

    # Produce figures
    print("\nGenerating figures...")

    plot_diagnostic(
        run_data, cs_series, active_period,
        os.path.join(out_dir, f"{comparison_name}_diagnostic.png"),
        comparison_name,
    )

    plot_annotated(
        run_data, cs_series, active_period,
        os.path.join(out_dir, f"{comparison_name}_annotated.png"),
        comparison_name,
    )

    if PLOT_TWO_PERIOD and cs_series:
        plot_two_period(
            run_data, cs_series,
            os.path.join(out_dir, f"{comparison_name}_two_period.png"),
            comparison_name,
        )

    if PLOT_RESIDUALS and cs_series:
        plot_residuals(
            run_data, cs_series, active_period,
            os.path.join(out_dir, f"{comparison_name}_residuals.png"),
            comparison_name,
        )

    print(f"\n✓ All figures saved to:\n  {out_dir}")
    print("=" * 70)


if __name__ == "__main__":
    main()
