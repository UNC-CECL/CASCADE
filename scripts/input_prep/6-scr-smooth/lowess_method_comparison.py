"""
CoastSat LRR smoothing along the island, transect-based and domain-averaged LOWESS side by side.

    python scripts/input_prep/6-scr-smooth/lowess_method_comparison.py

Both methods for every window size, both periods; overview, smoothed-only,
combined and method-comparison figures. Details: scripts/input_prep/6-scr-smooth/README.md.

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-09-30
"""

# pathlib first: every path in CONFIG is built from PROJECT_BASE_DIR
import pathlib

# Repo root, found by searching upward
_PATH_REPO = next(_p for _p in pathlib.Path(__file__).resolve().parents
                  if (_p / "pyproject.toml").exists())

# Paths anchored on the repo root, never typed (the old literals had broken)
PROJECT_BASE_DIR = next(
    q for q in pathlib.Path(__file__).resolve().parents
    if (q / "pyproject.toml").exists()
)


# --- CONFIG ------------------------------------------------------------------
# Smoothing x-axis

# 'along_coast_m' keeps the window physical; 'transect_id' spaces unevenly
TRANSECT_X_AXIS = "along_coast_m"   # "along_coast_m" | "transect_id"

# Domain-mode inputs

# Resolved through hat_observed_rates.py (2026-09-18), not typed.
import sys as _sys
from pathlib import Path as _RP
_sys.path.insert(0, str(next(_q for _q in _RP(__file__).resolve().parents
                             if (_q / "pyproject.toml").exists()) / "scripts"))
from site_layer import hat_observed_rates as _obs  # noqa: E402
DOMAIN_CSV_1984_2004 = str(_obs.domain_csv(1984, 2004))
DOMAIN_CSV_2004_2024 = str(_obs.domain_csv(2004, 2024))
CS_DOMAIN_COL = "domain_number"
CS_LRR_COL    = "mean_lrr"
CS_STD_COL    = "std_lrr"

# Transect-mode inputs

# Point to your transect_lrr_full.csv files for each period
TRANSECT_CSV_1984_2004 = str(_obs.lrr_csv(1984, 2004))
TRANSECT_CSV_2004_2024 = str(_obs.lrr_csv(2004, 2024))

# Column names in your transect CSV (transect_lrr_full.csv)
T_TRANSECT_ID_COL = "transect_id"    # string IDs — converted to sequential int internally
T_ALONG_COAST_COL = None             # not in CSV — derived from domain position automatically
T_DOMAIN_COL      = "domain_number"
T_LRR_COL         = "lrr_m_yr"
T_STD_COL         = "unc_m_yr"       # uncertainty column; set to None to skip

# Optional: filter to only point_in_polygon domain matches
FILTER_POINT_IN_POLYGON = False
T_MATCH_METHOD_COL      = "match_method"

# Domain geometry
DOMAIN_MIN       = 1
DOMAIN_MAX       = 90
DOMAIN_SPACING_M = 500   # metres per CASCADE domain

# LOWESS window

# Physical window width in km — applies to both modes
LOWESS_WINDOW_KM = 3.5   # primary smoothing window (7 domains)

# Window sizes (km) tested in sensitivity / comparison figures
COMPARE_WINDOWS_KM = [2.5, 3.5, 5.0]   # 5, 7, 10 domains

# Southern boundary guard

# Domains 1..N are dropped from the SMOOTHED series
SKIP_SOUTHERN_DOMAINS = 10

# Geographic annotations

# This block used to hold a copy of the town spans, the village centres, the piers

# Output

# Products live under data/hatteras_init/<stage>/, beside every other input_prep stage's output
OUTPUT_DIR = str(_obs.SMOOTH_METHOD_COMPARISON)
# -----------------------------------------------------------------------------

import os
import sys

# Windows consoles default to cp1252, which cannot encode the arrows and en-dashes in the progress output
if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from matplotlib.lines import Line2D
from matplotlib.transforms import blended_transform_factory
from statsmodels.nonparametric.smoothers_lowess import lowess
import warnings
warnings.filterwarnings("ignore")

# The house figure style and the site's annotation config are siblings in scripts/
sys.path.insert(0, str(PROJECT_BASE_DIR / "scripts"))
from site_layer.hat_figure_style import (            # noqa: E402
    C, C_1984, C_1997, DOMAIN_AXIS_LABEL, INK, INK_MUTED, _title, apply_style,
    caption, figsize, open_frame, save, town_bands)
from site_layer.hatteras_site_config import HATTERAS_ANNOTATIONS as ANN   # noqa: E402

# Subfolders are created automatically in main()

# Style

# One typographic and colour standard for every Hatteras figure, in scripts/site_layer/hat_figure_style.py
apply_style()

# Two periods drawn together are the house vintage pair
C_PERIOD_1984 = C_1984
C_PERIOD_2004 = C_1997

# The withheld-LOWESS band, and town_bands()' village shade (repeated for the legend patch)
C_SKIP_ZONE = C["WATER"]
TOWN_SHADE = "0.94"

# Three LOWESS windows are compared on the sweep figures
C_WINDOWS = [C["BASE"], C["ACCENT"], C["REF"]]

# Frac / window helpers

# Convert a physical window width (km) to a LOWESS frac for n_points
def km_to_frac(window_km, n_points, spacing_m):
    k = (window_km * 1000.0) / spacing_m
    return float(np.clip(k / n_points, 0.02, 1.0))


# Median spacing between consecutive sorted x values (positive diffs only)
def estimate_spacing(x_values):
    arr   = np.sort(np.asarray(x_values, dtype=float))
    diffs = np.diff(arr)
    pos   = diffs[diffs > 0]
    return float(np.median(pos)) if len(pos) else 1.0


# LOWESS frac for domain-space smoothing at a given physical window width
def domain_frac(window_km=LOWESS_WINDOW_KM):
    n = DOMAIN_MAX - DOMAIN_MIN + 1
    return km_to_frac(window_km, n, DOMAIN_SPACING_M)


# LOWESS frac for transect-space smoothing at a given physical window width
def transect_frac(n_transects, spacing_m, window_km=LOWESS_WINDOW_KM):
    return km_to_frac(window_km, n_transects, spacing_m)

# Data loading

# Load domain-averaged LRR summary CSV
def load_domain_csv(path, period_label):
    if path is None or not os.path.exists(path):
        print(f"  Domain CSV ({period_label}): SKIPPED — not found: {path}")
        return None
    df = pd.read_csv(path)

    # Strip any filename prefix accidentally prepended to column names e.g
    df.columns = [c.split(".csv")[-1] if ".csv" in c else c for c in df.columns]

    # Resolve domain column — try configured name then common alternatives
    domain_col = CS_DOMAIN_COL
    if domain_col not in df.columns:
        # Also check for columns that contain the target name as a substring
        matches = [c for c in df.columns if CS_DOMAIN_COL in c or c in CS_DOMAIN_COL]
        fallbacks = matches + ["domain", "Domain", "domain_id", "DOMAIN"]
        for fb in fallbacks:
            if fb in df.columns:
                domain_col = fb
                print(f"  Note: '{CS_DOMAIN_COL}' not found, using '{domain_col}' instead")
                break
        else:
            print(f"  ERROR: cannot find domain column. Available: {list(df.columns)}")
            return None

    for col in [domain_col, CS_LRR_COL, CS_STD_COL]:
        df[col] = pd.to_numeric(df[col], errors="coerce")
    df = df.dropna(subset=[domain_col, CS_LRR_COL])
    df[domain_col] = df[domain_col].astype(int)
    df = df[(df[domain_col] >= DOMAIN_MIN) & (df[domain_col] <= DOMAIN_MAX)]
    df = df[[domain_col, CS_LRR_COL, CS_STD_COL]].rename(columns={
        domain_col:  "domain",
        CS_LRR_COL:  "cs_lrr",
        CS_STD_COL:  "cs_std",
    }).sort_values("domain").reset_index(drop=True)
    print(f"  Domain CSV ({period_label}): {len(df)} domains  "
          f"LRR range {df['cs_lrr'].min():+.2f}–{df['cs_lrr'].max():+.2f} m/yr")
    return df

# Load transect-level LRR CSV
def load_transect_csv(path, period_label):
    if path is None or not os.path.exists(path):
        print(f"  Transect CSV ({period_label}): SKIPPED — not found: {path}")
        return None

    df = pd.read_csv(path)

    # Coerce required numeric columns
    for col in [T_DOMAIN_COL, T_LRR_COL]:
        df[col] = pd.to_numeric(df[col], errors="coerce")
    if T_STD_COL and T_STD_COL in df.columns:
        df[T_STD_COL] = pd.to_numeric(df[T_STD_COL], errors="coerce")

    # Drop rows missing domain or LRR
    df = df.dropna(subset=[T_TRANSECT_ID_COL, T_LRR_COL, T_DOMAIN_COL])
    df[T_DOMAIN_COL] = df[T_DOMAIN_COL].astype(int)
    df = df[(df[T_DOMAIN_COL] >= DOMAIN_MIN) & (df[T_DOMAIN_COL] <= DOMAIN_MAX)]

    # Optional match-method filter
    if FILTER_POINT_IN_POLYGON and T_MATCH_METHOD_COL in df.columns:
        before = len(df)
        df = df[df[T_MATCH_METHOD_COL] == "point_in_polygon"]
        print(f"  point_in_polygon filter: {before} → {len(df)} transects")

    # Sort by domain then transect ID string (zero-padded suffix sorts correctly)
    df = df.sort_values([T_DOMAIN_COL, T_TRANSECT_ID_COL]).reset_index(drop=True)

    # Replace string transect IDs with sequential integers based on sort order
    df["transect_id"] = np.arange(1, len(df) + 1)

    # Derive along_coast_m from domain position if not present in CSV
    if T_ALONG_COAST_COL is None or T_ALONG_COAST_COL not in df.columns:
        def _spread_within_domain(grp):
            n         = len(grp)
            dom_start = (grp[T_DOMAIN_COL].iloc[0] - 1) * DOMAIN_SPACING_M
            offsets   = (np.arange(n) + 0.5) * (DOMAIN_SPACING_M / n)
            grp       = grp.copy()
            grp["along_coast_m"] = dom_start + offsets
            return grp
        df = df.groupby(T_DOMAIN_COL, group_keys=False).apply(_spread_within_domain)
        print(f"  along_coast_m: derived from domain position (not in CSV)")
    else:
        df["along_coast_m"] = pd.to_numeric(df[T_ALONG_COAST_COL], errors="coerce")

    # Standardise remaining column names
    rename = {T_DOMAIN_COL: "domain", T_LRR_COL: "lrr"}
    if T_STD_COL and T_STD_COL in df.columns:
        rename[T_STD_COL] = "lrr_std"
    df = df.rename(columns=rename).reset_index(drop=True)

    spacing = estimate_spacing(df["along_coast_m"].values)
    print(f"  Transect CSV ({period_label}): {len(df)} transects  "
          f"est. spacing {spacing:.1f} m  "
          f"LRR range {df['lrr'].min():+.2f}–{df['lrr'].max():+.2f} m/yr")
    return df

# Smoothing

# LOWESS smoother
def apply_lowess(x, values, frac):
    valid = ~np.isnan(values)
    if valid.sum() < 5:
        return values.copy()
    smoothed = np.full_like(values, np.nan, dtype=float)
    result   = lowess(values[valid], x[valid], frac=frac, return_sorted=True)
    smoothed[valid] = np.interp(x[valid], result[:, 0], result[:, 1])
    return smoothed


# Apply LOWESS in domain space
def smooth_domain_df(df, window_km=LOWESS_WINDOW_KM):
    df   = df.copy()
    frac = domain_frac(window_km)
    df["cs_lrr_smooth"] = apply_lowess(
        df["domain"].values.astype(float),
        df["cs_lrr"].values,
        frac,
    )
    if SKIP_SOUTHERN_DOMAINS > 0:
        df.loc[df["domain"] <= SKIP_SOUTHERN_DOMAINS, "cs_lrr_smooth"] = np.nan
    return df


# Apply LOWESS in transect space
def smooth_transect_df(df, window_km=LOWESS_WINDOW_KM):
    df      = df.copy()
    spacing = estimate_spacing(df["along_coast_m"].values)
    frac    = transect_frac(len(df), spacing, window_km)

    x = (df["along_coast_m"].values.astype(float)
         if TRANSECT_X_AXIS == "along_coast_m"
         else df["transect_id"].values.astype(float))

    df["lrr_smooth"] = apply_lowess(x, df["lrr"].values, frac)
    # Masked on domain, not on along_coast_m
    if SKIP_SOUTHERN_DOMAINS > 0:
        df.loc[df["domain"] <= SKIP_SOUTHERN_DOMAINS, "lrr_smooth"] = np.nan
    df["_x_smooth"]  = x   # stored so plot functions don't recompute
    return df


# Average smoothed (and raw) transect values within each CASCADE domain
def aggregate_to_domains(t_df):
    grp = t_df.groupby("domain")
    domain_df = pd.DataFrame({
        "domain":        grp["lrr"].mean().index,
        "cs_lrr":        grp["lrr"].mean().values,
        "cs_std":        grp["lrr"].std(ddof=1).values,
        "cs_lrr_smooth": grp["lrr_smooth"].mean().values,
    }).reset_index(drop=True)
    return domain_df.sort_values("domain").reset_index(drop=True)

# Annotation helpers

# Villages via town_bands(); shoals, piers and groins from HATTERAS_ANNOTATIONS

# Sentences the figures used to carry on the canvas. They belong in a caption.
CAP_ENDPOINTS = ("Domain 1 is at Cape Point in the south and domain 90 at Pea "
                 "Island in the north; each domain is 500 m of shoreline.")
CAP_MARKS = ("Grey bands name the village spans, amber bands the Avon and "
             "Wimble shoal zones, dash-dot lines the Avon and Rodanthe piers "
             "and the dotted line the Buxton groin.")
CAP_GUARD = ("The LOWESS curve is withheld over the southernmost "
             f"{SKIP_SOUTHERN_DOMAINS} domains (shaded), where a local linear "
             "fit extrapolates rather than smooths and Oregon Inlet dominates; "
             "the hindcast applies the same guard."
             if SKIP_SOUTHERN_DOMAINS > 0 else "")


# The value unchanged (an x transform that does nothing)
def _identity(d):
    return d


# A GIS domain number as along-coast metres
def _domain_to_m(d):
    return (d - 0.5) * DOMAIN_SPACING_M


# Shoal zones, piers and groins from the site config
def _reference_marks(ax, to_x=_identity, label_shoals=True):
    trans = blended_transform_factory(ax.transData, ax.transAxes)
    for name, (lo, hi) in ANN.shoal_zones.items():
        ax.axvspan(to_x(lo - 0.5), to_x(hi + 0.5), color=ANN.color_shoal,
                   alpha=0.10, lw=0, zorder=0)
        if label_shoals:
            ax.text(to_x((lo + hi) / 2.0), 0.02, name, transform=trans,
                    ha="center", va="bottom", fontsize=7, style="italic",
                    color=INK_MUTED, zorder=1, clip_on=True)
    for dom, _lbl_y in ANN.piers.values():
        ax.axvline(to_x(dom), color=ANN.color_pier, lw=0.9, ls="-.",
                   alpha=0.85, zorder=2)
    for dom in ANN.groins.values():
        ax.axvline(to_x(dom), color=ANN.color_groin, lw=0.9, ls=":",
                   alpha=0.85, zorder=2)


# Band from the start of the reach to `hi`, in whatever x-units the axis uses
def _shade_boundary_zone(ax, hi):
    if SKIP_SOUTHERN_DOMAINS > 0:
        ax.axvspan(ax.get_xlim()[0], hi, facecolor=C_SKIP_ZONE, alpha=0.30,
                   lw=0.0, zorder=0)


# Villages, shoals, piers, groins and the guard zone on a domain axis
def add_domain_annotations(ax, label_shoals=True):
    _reference_marks(ax, label_shoals=label_shoals)
    town_bands(ax, strip=0.085)
    _shade_boundary_zone(ax, SKIP_SOUTHERN_DOMAINS + 0.5)


# The same marks on an alongshore-metres axis
def add_transect_annotations(ax, label_shoals=True):
    _reference_marks(ax, to_x=_domain_to_m, label_shoals=label_shoals)
    town_bands(ax, strip=0.085,
               spans={name: (_domain_to_m(lo - 0.5), _domain_to_m(hi + 0.5))
                      for name, (lo, hi) in ANN.town_spans.items()})
    # Domains 1..N occupy [0, N*500) m, so the cut is at N * DOMAIN_SPACING_M.
    _shade_boundary_zone(ax, SKIP_SOUTHERN_DOMAINS * DOMAIN_SPACING_M)


# Legend handles for the marks
def annotation_legend_handles():
    return [
        Patch(facecolor=TOWN_SHADE, edgecolor="none", label="village span"),
        Patch(facecolor=ANN.color_shoal, alpha=0.30, edgecolor="none",
              label="shoal zone"),
        Line2D([0], [0], color=ANN.color_pier, lw=1.0, ls="-.", label="pier"),
        Line2D([0], [0], color=ANN.color_groin, lw=1.0, ls=":", label="groin"),
    ] + ([Patch(facecolor=C_SKIP_ZONE, alpha=0.40, edgecolor="none",
                label="LOWESS withheld")]
         if SKIP_SOUTHERN_DOMAINS > 0 else [])


# A frameless legend below the panels
def _outside_legend(fig, handles, ncol=4):
    fig.legend(handles=handles, loc="outside lower center", ncol=ncol,
               frameon=False, fontsize=7.5)


# Raw domain means across the withheld zone, so it is not simply blank
def draw_raw_in_guard_zone(ax, df, color, label=None, col="cs_lrr"):
    if SKIP_SOUTHERN_DOMAINS <= 0:
        return
    z = df[df["domain"] <= SKIP_SOUTHERN_DOMAINS]
    if z.empty:
        return
    ax.plot(z["domain"], z[col], color=color, lw=1.0, ls=":",
            marker="o", ms=2.5, alpha=0.75, zorder=2, label=label)


# A 'no data' note in an empty panel
def _no_data(ax, period):
    ax.text(0.5, 0.5, f"no data for {period}", transform=ax.transAxes,
            ha="center", va="center", color=INK_MUTED)


# The GIS-domain x axis
def style_domain_axis(ax, is_bottom=True):
    ax.set_xlim(DOMAIN_MIN - 0.5, DOMAIN_MAX + 0.5)
    ax.axhline(0, color=INK_MUTED, lw=0.6, ls="--", zorder=1)
    ax.grid(axis="y")
    ax.set_axisbelow(True)
    open_frame(ax)
    if is_bottom:
        ax.set_xlabel(DOMAIN_AXIS_LABEL)


# The alongshore x axis
def style_transect_axis(ax, x_values, is_bottom=True):
    ax.set_xlim(x_values.min() - 1, x_values.max() + 1)
    ax.axhline(0, color=INK_MUTED, lw=0.6, ls="--", zorder=1)
    ax.grid(axis="y")
    ax.set_axisbelow(True)
    open_frame(ax)
    if is_bottom:
        ax.set_xlabel("along-coast distance (m), south to north"
                      if TRANSECT_X_AXIS == "along_coast_m"
                      else "transect, numbered south to north")

# Domain-space figures Works identically for both modes

# Both periods as two panels, raw optional
def _domain_two_panel(d1984, d2004, show_raw, out_path, cap):
    configs   = [(d1984, "1984–2004", C_1984),
                 (d2004, "2004–2024", C_1997)]
    fig, axes = plt.subplots(2, 1, figsize=figsize("double", height=5.8),
                             sharex=True, constrained_layout=True)
    caption(fig, cap)

    for i, (ax, (df, period, color)) in enumerate(zip(axes, configs)):
        is_bottom = (i == len(configs) - 1)
        _title(ax, i, period)
        if df is None:
            _no_data(ax, period)
            continue
        if show_raw:
            ax.fill_between(df["domain"],
                            df["cs_lrr"] - df["cs_std"],
                            df["cs_lrr"] + df["cs_std"],
                            color=color, alpha=0.12, lw=0, zorder=0,
                            label="±1 s.d. within the domain")
            ax.plot(df["domain"], df["cs_lrr"],
                    color=color, lw=0.7, alpha=0.45, marker="o", ms=1.8,
                    zorder=1,
                    label="domain-averaged linear regression rate, unsmoothed")
        else:
            # Nothing else covers the guard zone on this figure.
            draw_raw_in_guard_zone(
                ax, df, color,
                label="domain-averaged rate where the curve is withheld")
        ax.plot(df["domain"], df["cs_lrr_smooth"],
                color=color, lw=1.8, zorder=3,
                label=f"LOWESS, {LOWESS_WINDOW_KM:g} km window")
        ax.set_ylabel("shoreline change rate (m/yr)")
        style_domain_axis(ax, is_bottom)
        add_domain_annotations(ax, label_shoals=is_bottom)
        if i == 0:
            ax.legend(loc="lower left", ncol=1)

    _outside_legend(fig, annotation_legend_handles(), ncol=5)
    save(fig, out_path, close=True)
    print(f"  Saved: {os.path.basename(out_path)}")


# Raw and smoothed rates by domain, both periods
def plot_domain_overview(d1984, d2004, out_path, method=""):
    _domain_two_panel(
        d1984, d2004, show_raw=True, out_path=out_path,
        cap=("CoastSat shoreline change rate by CASCADE domain, "
             f"{method}. The pale line and band are the unsmoothed "
             "domain mean and its standard deviation over the transects in "
             "the domain; the heavy line is the LOWESS curve at a "
             f"{LOWESS_WINDOW_KM:g} km window. " + CAP_GUARD + " " +
             CAP_ENDPOINTS + " " + CAP_MARKS),
    )


# The smoothed rates alone, both periods
def plot_domain_smoothed_only(d1984, d2004, out_path, method=""):
    _domain_two_panel(
        d1984, d2004, show_raw=False, out_path=out_path,
        cap=("The smoothed CoastSat shoreline change rate alone, by CASCADE "
             f"domain, {method}: the LOWESS curve at a "
             f"{LOWESS_WINDOW_KM:g} km window, with the unsmoothed domain means "
             "shown only where the curve is withheld. " + CAP_GUARD + " " +
             CAP_ENDPOINTS + " " + CAP_MARKS),
    )


# Both periods on one panel, smoothed
def plot_domain_combined(d1984, d2004, out_path, method=""):
    fig, ax = plt.subplots(figsize=figsize("double", height=3.6),
                           constrained_layout=True)
    caption(fig, "Both hindcast periods on one panel, " + method + ": the "
                 f"LOWESS curve at a {LOWESS_WINDOW_KM:g} km window, shaded by "
                 "one standard deviation of the transects in each domain. " +
                 CAP_GUARD + " " + CAP_ENDPOINTS + " " + CAP_MARKS)
    labelled_zone = False
    for df, period, color in [(d1984, "1984–2004", C_1984),
                               (d2004, "2004–2024", C_1997)]:
        if df is None:
            continue
        draw_raw_in_guard_zone(
            ax, df, color,
            label=None if labelled_zone
            else "domain-averaged rate where the curve is withheld")
        labelled_zone = True
        ax.fill_between(df["domain"],
                        df["cs_lrr_smooth"] - df["cs_std"],
                        df["cs_lrr_smooth"] + df["cs_std"],
                        color=color, alpha=0.12, lw=0)
        ax.plot(df["domain"], df["cs_lrr_smooth"],
                color=color, lw=1.8, label=period)
    ax.set_ylabel("shoreline change rate (m/yr)")
    style_domain_axis(ax, is_bottom=True)
    add_domain_annotations(ax)
    handles, _ = ax.get_legend_handles_labels()
    _outside_legend(fig, handles + annotation_legend_handles(), ncol=4)
    save(fig, out_path, close=True)
    print(f"  Saved: {os.path.basename(out_path)}")


# 3-panel bandwidth sensitivity figure for a single period (domain space)
def plot_domain_sensitivity(df, period_label, color, out_path, method=""):
    fracs = [domain_frac(w) for w in COMPARE_WINDOWS_KM]
    fig, axes = plt.subplots(len(COMPARE_WINDOWS_KM), 1,
                             figsize=figsize("double", height=7.2),
                             sharex=True, constrained_layout=True)
    caption(fig, f"LOWESS bandwidth sensitivity, {period_label}, {method}. One "
                 "panel per window: the pale line is the unsmoothed "
                 "domain-averaged linear regression rate and the heavy line "
                 "the LOWESS curve at that window. " + CAP_GUARD + " " +
                 CAP_ENDPOINTS + " " + CAP_MARKS)
    for j, (ax, frac, km) in enumerate(zip(axes, fracs, COMPARE_WINDOWS_KM)):
        is_bottom = (j == len(fracs) - 1)
        ndom = int(round(km * 1000 / DOMAIN_SPACING_M))
        m = smooth_domain_df(df, window_km=km)
        ax.plot(df["domain"], df["cs_lrr"],
                color=color, lw=0.7, alpha=0.40, marker="o", ms=1.8,
                label="domain-averaged rate, unsmoothed")
        ax.plot(m["domain"], m["cs_lrr_smooth"],
                color=color, lw=1.8, label="LOWESS")
        _title(ax, j, f"{km:g} km window · {ndom} domains · frac {frac:.3f}")
        ax.set_ylabel("rate (m/yr)")
        style_domain_axis(ax, is_bottom)
        add_domain_annotations(ax, label_shoals=is_bottom)
        if j == 0:
            ax.legend(loc="lower left")
    _outside_legend(fig, annotation_legend_handles(), ncol=5)
    save(fig, out_path, close=True)
    print(f"  Saved: {os.path.basename(out_path)}")


# All window sizes overlaid in domain space — both periods
def plot_domain_window_comparison(d1984, d2004, out_path, method=""):
    configs   = [(d1984, "1984–2004", C_1984),
                 (d2004, "2004–2024", C_1997)]
    fig, axes = plt.subplots(2, 1, figsize=figsize("double", height=6.2),
                             sharex=True, constrained_layout=True)
    caption(fig, "The three LOWESS windows overlaid, " + method + ", one panel "
                 "per hindcast period. The pale line is the unsmoothed "
                 "domain-averaged linear regression rate. " + CAP_GUARD + " " +
                 CAP_ENDPOINTS + " " + CAP_MARKS)
    for i, (ax, (df, period, pcol)) in enumerate(zip(axes, configs)):
        is_bottom = (i == len(configs) - 1)
        _title(ax, i, period)
        if df is None:
            _no_data(ax, period)
            continue
        ax.plot(df["domain"], df["cs_lrr"],
                color=pcol, lw=0.7, alpha=0.40, marker="o", ms=1.8,
                label="domain-averaged rate, unsmoothed")
        for km, wc in zip(COMPARE_WINDOWS_KM, C_WINDOWS):
            m    = smooth_domain_df(df, window_km=km)
            ndom = int(round(km * 1000 / DOMAIN_SPACING_M))
            ax.plot(m["domain"], m["cs_lrr_smooth"],
                    color=wc, lw=1.8, label=f"{km:g} km ({ndom} domains)")
        ax.set_ylabel("rate (m/yr)")
        style_domain_axis(ax, is_bottom)
        add_domain_annotations(ax, label_shoals=is_bottom)
        if i == 0:
            ax.legend(loc="lower left", ncol=2)
    _outside_legend(fig, annotation_legend_handles(), ncol=5)
    save(fig, out_path, close=True)
    print(f"  Saved: {os.path.basename(out_path)}")

# Transect-space figures (produced only in "transect" mode)

# Raw transect scatter + LOWESS smoothed curve in along-coast space, with domain-averaged LRR ...
def plot_transect_overview(t1984, t2004, d1984, d2004, out_path):
    configs   = [(t1984, d1984, "1984–2004", C_1984),
                 (t2004, d2004, "2004–2024", C_1997)]
    fig, axes = plt.subplots(2, 1, figsize=figsize("double", height=5.8),
                             sharex=True, constrained_layout=True)
    caption(fig, "CoastSat shoreline change rate at the resolution it is "
                 "measured. Dots are the per-transect linear regression rate, "
                 "unsmoothed; the heavy line is the LOWESS curve fitted to "
                 f"those transects at a {LOWESS_WINDOW_KM:g} km window; the "
                 "open markers are the mean of the transects in each 500 m "
                 "CASCADE domain, which is what domain averaging keeps. "
                 "Faint verticals are the domain boundaries. " + CAP_GUARD +
                 " " + CAP_ENDPOINTS + " " + CAP_MARKS)
    for i, (ax, (t_df, d_df, period, color)) in enumerate(zip(axes, configs)):
        is_bottom = (i == len(configs) - 1)
        _title(ax, i, period)
        if t_df is None:
            _no_data(ax, period)
            continue
        x = t_df["_x_smooth"].values

        # Domain boundary lines at every 500 m — drawn first so they sit behind data
        for domain_n in range(DOMAIN_MIN, DOMAIN_MAX + 2):
            boundary_m = (domain_n - 1) * DOMAIN_SPACING_M
            ax.axvline(boundary_m, color="0.88", lw=0.3, ls="-", zorder=0)

        ax.scatter(x, t_df["lrr"], color=color, s=3, alpha=0.20, lw=0,
                   zorder=1,
                   label="per-transect linear regression rate, unsmoothed")
        ax.plot(x, t_df["lrr_smooth"], color=color, lw=1.8, zorder=3,
                label=f"LOWESS on the transects, {LOWESS_WINDOW_KM:g} km window")
        if d_df is not None:
            x_dom = (d_df["domain"].values - 0.5) * DOMAIN_SPACING_M
            ax.plot(x_dom, d_df["cs_lrr"].values,
                    color=INK, lw=0.9, ls="--",
                    marker="o", ms=3, markerfacecolor="white",
                    markeredgewidth=0.9, zorder=4, alpha=0.8,
                    label="mean of the transects in each domain")
        ax.set_ylabel("shoreline change rate (m/yr)")
        style_transect_axis(ax, x, is_bottom)
        add_transect_annotations(ax, label_shoals=is_bottom)
        if i == 0:
            ax.legend(loc="lower left")
    _outside_legend(fig, annotation_legend_handles(), ncol=5)
    save(fig, out_path, close=True)
    print(f"  Saved: {os.path.basename(out_path)}")


# Window sensitivity in transect space — all km windows overlaid, both periods
def plot_transect_window_comparison(t1984, t2004, out_path):
    configs = [(t1984, "1984–2004", C_1984),
               (t2004, "2004–2024", C_1997)]
    fig, axes = plt.subplots(2, 1, figsize=figsize("double", height=6.2),
                             sharex=True, constrained_layout=True)
    caption(fig, "The three LOWESS windows overlaid at transect resolution, "
                 "one panel per hindcast period. Dots are the per-transect "
                 "linear regression rate, unsmoothed. " + CAP_GUARD + " " +
                 CAP_ENDPOINTS + " " + CAP_MARKS)
    for i, (ax, (df, period, pcol)) in enumerate(zip(axes, configs)):
        is_bottom = (i == len(configs) - 1)
        _title(ax, i, period)
        if df is None:
            _no_data(ax, period)
            continue
        x       = df["_x_smooth"].values
        spacing = estimate_spacing(df["along_coast_m"].values)
        ax.scatter(x, df["lrr"], color=pcol, s=3, alpha=0.18, lw=0,
                   label="per-transect rate, unsmoothed")
        for km, wc in zip(COMPARE_WINDOWS_KM, C_WINDOWS):
            frac     = transect_frac(len(df), spacing, km)
            smoothed = apply_lowess(x, df["lrr"].values, frac)
            ndom     = int(round(km * 1000 / DOMAIN_SPACING_M))
            ax.plot(x, smoothed, color=wc, lw=1.8,
                    label=f"{km:g} km ({ndom} domains)")
        ax.set_ylabel("rate (m/yr)")
        style_transect_axis(ax, x, is_bottom)
        add_transect_annotations(ax, label_shoals=is_bottom)
        if i == 0:
            ax.legend(loc="lower left", ncol=2)
    _outside_legend(fig, annotation_legend_handles(), ncol=5)
    save(fig, out_path, close=True)
    print(f"  Saved: {os.path.basename(out_path)}")


# Transect-smoothed windows in domain space Smooths at transect level for each window

# For each window in COMPARE_WINDOWS_KM
def plot_transect_windows_domain_space(t1984, t2004, out_path):
    configs   = [(t1984, "1984–2004", C_1984),
                 (t2004, "2004–2024", C_1997)]
    fig, axes = plt.subplots(2, 1, figsize=figsize("double", height=6.2),
                             sharex=True, constrained_layout=True)
    caption(fig, "The three LOWESS windows overlaid, each fitted to the "
                 "individual transects and then averaged to CASCADE domains, "
                 "one panel per hindcast period. The pale line is the "
                 "unsmoothed mean of the transects in each domain. " +
                 CAP_GUARD + " " + CAP_ENDPOINTS + " " + CAP_MARKS)

    for i, (ax, (t_df, period, pcol)) in enumerate(zip(axes, configs)):
        is_bottom = (i == len(configs) - 1)
        _title(ax, i, period)
        if t_df is None:
            _no_data(ax, period)
            continue

        # Raw domain means — mean of raw transect LRRs, same for all windows
        d_raw = aggregate_to_domains(t_df)
        ax.plot(d_raw["domain"], d_raw["cs_lrr"],
                color=pcol, lw=0.7, alpha=0.40, marker="o", ms=1.8,
                label="mean of the transects in each domain, unsmoothed")

        spacing = estimate_spacing(t_df["along_coast_m"].values)

        # One smoothed curve per window: smooth transects → aggregate to domains
        for km, wc in zip(COMPARE_WINDOWS_KM, C_WINDOWS):
            frac     = transect_frac(len(t_df), spacing, km)
            t_smooth = smooth_transect_df(t_df, window_km=km)
            d_agg    = aggregate_to_domains(t_smooth)
            ndom     = int(round(km * 1000 / DOMAIN_SPACING_M))
            ax.plot(d_agg["domain"], d_agg["cs_lrr_smooth"],
                    color=wc, lw=1.8,
                    label=f"{km:g} km ({ndom} domains, frac {frac:.3f})")

        ax.set_ylabel("rate (m/yr)")
        style_domain_axis(ax, is_bottom)
        add_domain_annotations(ax, label_shoals=is_bottom)
        if i == 0:
            ax.legend(loc="lower left", ncol=2)

    _outside_legend(fig, annotation_legend_handles(), ncol=5)
    save(fig, out_path, close=True)
    print(f"  Saved: {os.path.basename(out_path)}")


# Transect sensitivity (matches plot_transect_windows_domain_space)

# 3-panel bandwidth sensitivity — transect mode
def plot_transect_sensitivity(t_df, period_label, color, out_path, method=""):
    fig, axes = plt.subplots(len(COMPARE_WINDOWS_KM), 1,
                             figsize=figsize("double", height=7.2),
                             sharex=True, constrained_layout=True)
    caption(fig, f"LOWESS bandwidth sensitivity, {period_label}, {method}. One "
                 "panel per window: the LOWESS is fitted to the individual "
                 "transects and then averaged to CASCADE domains, against the "
                 "unsmoothed mean of the transects in each domain. " +
                 CAP_GUARD + " " + CAP_ENDPOINTS + " " + CAP_MARKS)
    spacing = estimate_spacing(t_df["along_coast_m"].values)
    d_raw   = aggregate_to_domains(t_df)

    for j, (ax, km) in enumerate(zip(axes, COMPARE_WINDOWS_KM)):
        is_bottom = (j == len(COMPARE_WINDOWS_KM) - 1)
        frac  = transect_frac(len(t_df), spacing, km)
        ndom  = int(round(km * 1000 / DOMAIN_SPACING_M))
        t_smooth = smooth_transect_df(t_df, window_km=km)
        d_agg    = aggregate_to_domains(t_smooth)
        ax.plot(d_raw["domain"], d_raw["cs_lrr"],
                color=color, lw=0.7, alpha=0.40, marker="o", ms=1.8,
                label="mean of the transects in each domain, unsmoothed")
        ax.plot(d_agg["domain"], d_agg["cs_lrr_smooth"],
                color=color, lw=1.8, label="LOWESS")
        _title(ax, j, f"{km:g} km window · {ndom} domains · frac {frac:.3f}")
        ax.set_ylabel("rate (m/yr)")
        style_domain_axis(ax, is_bottom)
        add_domain_annotations(ax, label_shoals=is_bottom)
        if j == 0:
            ax.legend(loc="lower left")
    _outside_legend(fig, annotation_legend_handles(), ncol=5)
    save(fig, out_path, close=True)
    print(f"  Saved: {os.path.basename(out_path)}")


# Method comparison: transect-based against domain-averaged LOWESS, per window size

# The two methods are the point of these four figures, so they carry the colour
M_LS = ["-", "--", (0, (1, 1.4))]
LBL_TRANSECT = "LOWESS on the individual transects, averaged to domains"
LBL_DOMAIN   = "LOWESS on the domain averages"
LBL_RAW      = "per-transect rates averaged to domains, unsmoothed"

CAP_METHODS = ("Purple is the LOWESS fitted to the individual CoastSat "
               "transects and then averaged to CASCADE domains; grey is the "
               "LOWESS fitted to the domain averages directly. Dots are those "
               "domain averages before smoothing.")


# The raw domain means under a method panel
def _method_raw(ax, t_df):
    d_raw = aggregate_to_domains(t_df)
    ax.plot(d_raw["domain"], d_raw["cs_lrr"], color=INK_MUTED, lw=0,
            marker="o", ms=2.2, alpha=0.45, zorder=1, label=LBL_RAW)


# For each window in COMPARE_WINDOWS_KM, plots both
def plot_method_comparison(t1984, t2004, da1984, da2004, out_path):
    configs   = [(t1984, da1984, "1984–2004"),
                 (t2004, da2004, "2004–2024")]
    fig, axes = plt.subplots(2, 1, figsize=figsize("double", height=6.4),
                             sharex=True, constrained_layout=True)
    caption(fig, "Do the two ways of smoothing the CoastSat rates differ? "
                 "Each panel is one hindcast period and carries both methods "
                 "at all three LOWESS windows, the window as the line style. " +
                 CAP_METHODS + " " + CAP_GUARD + " " + CAP_ENDPOINTS + " " +
                 CAP_MARKS)

    for i, (ax, (t_df, da_df, period)) in enumerate(zip(axes, configs)):
        is_bottom = (i == len(configs) - 1)
        _title(ax, i, period)
        if t_df is None and da_df is None:
            _no_data(ax, period)
            continue

        if t_df is not None:
            _method_raw(ax, t_df)

        for km, ls in zip(COMPARE_WINDOWS_KM, M_LS):
            if t_df is not None:
                t_smooth = smooth_transect_df(t_df, window_km=km)
                d_agg    = aggregate_to_domains(t_smooth)
                ax.plot(d_agg["domain"], d_agg["cs_lrr_smooth"],
                        color=C["ACCENT"], lw=1.7, ls=ls, zorder=4)
            if da_df is not None:
                m = smooth_domain_df(da_df, window_km=km)
                ax.plot(m["domain"], m["cs_lrr_smooth"],
                        color=C["BASE"], lw=1.7, ls=ls, zorder=3)

        ax.set_ylabel("shoreline change rate (m/yr)")
        style_domain_axis(ax, is_bottom)
        add_domain_annotations(ax, label_shoals=is_bottom)

    handles = [
        Line2D([0], [0], color=INK_MUTED, lw=0, marker="o", ms=2.5,
               alpha=0.6, label=LBL_RAW),
        Line2D([0], [0], color=C["ACCENT"], lw=1.7, label=LBL_TRANSECT),
        Line2D([0], [0], color=C["BASE"], lw=1.7, label=LBL_DOMAIN),
    ] + [
        Line2D([0], [0], color=INK, lw=1.2, ls=ls,
               label=f"{km:g} km window ({int(round(km * 1000 / DOMAIN_SPACING_M))} domains)")
        for km, ls in zip(COMPARE_WINDOWS_KM, M_LS)
    ] + annotation_legend_handles()
    _outside_legend(fig, handles, ncol=3)
    save(fig, out_path, close=True)
    print(f"  Saved: {os.path.basename(out_path)}")


# Single-window method comparison
def plot_method_comparison_single(t1984, t2004, da1984, da2004,
                                   window_km, out_path):
    ndom  = int(round(window_km * 1000 / DOMAIN_SPACING_M))
    configs   = [(t1984, da1984, "1984–2004"),
                 (t2004, da2004, "2004–2024")]
    fig, axes = plt.subplots(2, 1, figsize=figsize("double", height=5.8),
                             sharex=True, constrained_layout=True)
    caption(fig, "Do the two ways of smoothing the CoastSat rates differ? "
                 f"Both are shown at a single {window_km:g} km LOWESS window "
                 f"({ndom} CASCADE domains), one panel per hindcast period. " +
                 CAP_METHODS + " " + CAP_GUARD + " " + CAP_ENDPOINTS + " " +
                 CAP_MARKS)

    for i, (ax, (t_df, da_df, period)) in enumerate(zip(axes, configs)):
        is_bottom = (i == len(configs) - 1)
        _title(ax, i, period)
        if t_df is None and da_df is None:
            _no_data(ax, period)
            continue

        if t_df is not None:
            _method_raw(ax, t_df)
            t_smooth = smooth_transect_df(t_df, window_km=window_km)
            d_agg    = aggregate_to_domains(t_smooth)
            ax.plot(d_agg["domain"], d_agg["cs_lrr_smooth"],
                    color=C["ACCENT"], lw=2.0, ls="-", zorder=4,
                    label=LBL_TRANSECT)

        if da_df is not None:
            m = smooth_domain_df(da_df, window_km=window_km)
            ax.plot(m["domain"], m["cs_lrr_smooth"],
                    color=C["BASE"], lw=2.0, ls="--", zorder=3,
                    label=LBL_DOMAIN)

        ax.set_ylabel("shoreline change rate (m/yr)")
        style_domain_axis(ax, is_bottom)
        add_domain_annotations(ax, label_shoals=is_bottom)

    handles, _ = axes[0].get_legend_handles_labels()
    _outside_legend(fig, handles + annotation_legend_handles(), ncol=4)
    save(fig, out_path, close=True)
    print(f"  Saved: {os.path.basename(out_path)}")

# Main Always runs both transect-based and domain-averaged smoothing

# Run: load both periods, every figure
def main():
    print("=" * 65)
    print("CoastSat LRR Smoothing — Hatteras Island")
    print(f"  X-axis      : {TRANSECT_X_AXIS}")
    print(f"  Window      : {LOWESS_WINDOW_KM} km  "
          f"({int(round(LOWESS_WINDOW_KM * 1000 / DOMAIN_SPACING_M))} domains)")
    print(f"  Compare     : {COMPARE_WINDOWS_KM} km  "
          f"({[int(round(w * 1000 / DOMAIN_SPACING_M)) for w in COMPARE_WINDOWS_KM]} domains)")
    print("  Both methods always run — outputs saved to subfolders.")
    print("=" * 65)

    # Create subfolders
    DIR_T       = os.path.join(OUTPUT_DIR, "01_transect_based")
    DIR_D       = os.path.join(OUTPUT_DIR, "02_domain_averaged")
    DIR_C       = os.path.join(OUTPUT_DIR, "03_cascade_inputs")
    DIR_COMPARE = os.path.join(OUTPUT_DIR, "04_method_comparison")
    for d in [DIR_T, DIR_D, DIR_C, DIR_COMPARE]:
        os.makedirs(d, exist_ok=True)

    # Load transect data
    print("\nLoading transect data...")
    t1984_raw = load_transect_csv(TRANSECT_CSV_1984_2004, "1984–2004")
    t2004_raw = load_transect_csv(TRANSECT_CSV_2004_2024, "2004–2024")
    t1984 = smooth_transect_df(t1984_raw) if t1984_raw is not None else None
    t2004 = smooth_transect_df(t2004_raw) if t2004_raw is not None else None
    td1984 = aggregate_to_domains(t1984) if t1984 is not None else None
    td2004 = aggregate_to_domains(t2004) if t2004 is not None else None

    # Load domain-averaged data
    print("\nLoading domain-averaged data...")
    da1984 = load_domain_csv(DOMAIN_CSV_1984_2004, "1984–2004")
    da2004 = load_domain_csv(DOMAIN_CSV_2004_2024, "2004–2024")
    if da1984 is not None: da1984 = smooth_domain_df(da1984)
    if da2004 is not None: da2004 = smooth_domain_df(da2004)

    # 01: Transect-based figures
    print("\n[01] Transect-based figures → 01_transect_based/")

    # Transect overview: raw scatter + LOWESS + domain averages overlaid
    plot_transect_overview(t1984, t2004, td1984, td2004,
        os.path.join(DIR_T, "transect_overview.png"))

    # Domain-space overview using transect-smoothed values
    plot_domain_overview(td1984, td2004,
        os.path.join(DIR_T, "overview_smoothed.png"),
        "smoothed on the individual transects")
    plot_domain_smoothed_only(td1984, td2004,
        os.path.join(DIR_T, "smoothed_only.png"),
        "smoothed on the individual transects")
    plot_domain_combined(td1984, td2004,
        os.path.join(DIR_T, "combined_periods.png"),
        "smoothed on the individual transects")

    # Window comparison in domain space (transect-smoothed)
    plot_transect_windows_domain_space(t1984, t2004,
        os.path.join(DIR_T, "window_comparison.png"))

    # Sensitivity — transect-based
    if t1984 is not None:
        plot_transect_sensitivity(t1984, "1984–2004", C_PERIOD_1984,
            os.path.join(DIR_T, "sensitivity_1984_2004.png"),
            "smoothed on the individual transects")
    if t2004 is not None:
        plot_transect_sensitivity(t2004, "2004–2024", C_PERIOD_2004,
            os.path.join(DIR_T, "sensitivity_2004_2024.png"),
            "smoothed on the individual transects")

    # 02: Domain-averaged figures
    print("\n[02] Domain-averaged figures → 02_domain_averaged/")

    plot_domain_window_comparison(da1984, da2004,
        os.path.join(DIR_D, "window_comparison.png"),
        "smoothed on the domain averages")

    if da1984 is not None:
        plot_domain_sensitivity(da1984, "1984–2004", C_PERIOD_1984,
            os.path.join(DIR_D, "sensitivity_1984_2004.png"),
            "smoothed on the domain averages")
    if da2004 is not None:
        plot_domain_sensitivity(da2004, "2004–2024", C_PERIOD_2004,
            os.path.join(DIR_D, "sensitivity_2004_2024.png"),
            "smoothed on the domain averages")

    # 04: Method comparison (Laura's request)
    print("\n[04] Method comparison figures → 04_method_comparison/")

    plot_method_comparison(t1984, t2004, da1984, da2004,
        os.path.join(DIR_COMPARE, "transect_vs_domain_smoothing.png"))

    # Individual window figures — one per window size
    for km in COMPARE_WINDOWS_KM:
        ndom  = int(round(km * 1000 / DOMAIN_SPACING_M))
        fname = f"transect_vs_domain_{ndom}domains_{km:.1f}km.png"
        plot_method_comparison_single(t1984, t2004, da1984, da2004,
            km, os.path.join(DIR_COMPARE, fname))

    # 03: CASCADE inputs (transect-based)
    print("\n[03] Exporting CASCADE inputs → 03_cascade_inputs/")

    parts = []
    for df, lbl in [(td1984, "1984_2004"), (td2004, "2004_2024")]:
        if df is None:
            continue
        t = df[["domain", "cs_lrr", "cs_std", "cs_lrr_smooth"]].copy()
        t.columns = ["domain"] + [f"{c}_{lbl}" for c in ["cs_lrr", "cs_std", "cs_lrr_smooth"]]
        parts.append(t)
    if parts:
        from functools import reduce
        table = reduce(lambda a, b: a.merge(b, on="domain", how="outer"), parts)
        fname = "cascade_lrr_inputs_transect_based.csv"
        table.sort_values("domain").to_csv(os.path.join(DIR_C, fname), index=False)
        print(f"  Saved: {fname}")

    # Also save full transect-level tables for reference
    for df, lbl in [(t1984, "1984_2004"), (t2004, "2004_2024")]:
        if df is None:
            continue
        cols = ["transect_id", "along_coast_m", "domain", "lrr", "lrr_smooth"]
        if "lrr_std" in df.columns:
            cols.insert(4, "lrr_std")
        fname = f"transect_lrr_{lbl}.csv"
        df[cols].to_csv(os.path.join(DIR_C, fname), index=False)
        print(f"  Saved: {fname}")

    print("\n" + "=" * 65)
    print("Done!  Output structure:")
    print(f"  01_transect_based/   — transect-smoothed figures")
    print(f"  02_domain_averaged/  — domain-averaged figures")
    print(f"  03_cascade_inputs/   — CASCADE-ready CSV + transect tables")
    print("=" * 65)


if __name__ == "__main__":
    main()
