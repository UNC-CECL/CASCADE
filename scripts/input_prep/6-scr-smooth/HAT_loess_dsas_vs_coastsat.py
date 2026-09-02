"""
DSAS vs CoastSat LRR Comparison — Hatteras Island  (SMOOTHED VERSION)
======================================================================
Extends the original comparison script with LOESS smoothing and
multiple visualization options for collaborator review.

New outputs (saved to OUTPUT_DIR)
----------------------------------
  overview_smoothed.png          – 2-panel both periods, raw + LOESS overlay
  smoothed_only_comparison.png   – 2-panel both periods, smoothed lines only
  smoothing_sensitivity.png      – 3-panel showing frac=0.10, 0.15, 0.20 side by side
  combined_sources.png           – single panel combining both periods + both sources
  scatter_smoothed_1978_1997.png – scatter of smoothed values, period 1
  scatter_smoothed_1997_2019.png – scatter of smoothed values, period 2

All original outputs are also regenerated.

comparison_table_smoothed.csv carries the raw values for every domain, but its
LOESS columns are blank across the southern boundary zone (domains
1..SKIP_SOUTHERN_DOMAINS), matching what the figures draw and why.

Smoothing method: LOESS (locally weighted scatterplot smoothing)
  - Applied independently to each series (DSAS and CoastSat)
  - Preserves large-scale spatial patterns while removing per-domain noise
"""

# pathlib must be imported before the CONFIG block because every path below is
# built from PROJECT_BASE_DIR at module level.
import pathlib

# ANCHORED, NOT TYPED. Every path below used to be an absolute literal: the
# output one had lost its drive ("/scripts/input_prep/...") and so wrote its
# figures to C:\scripts\ instead of into the repository, and the input ones
# still spelled the folder "input_preperation" and pointed at a CoastSat tree
# that has since moved under 5-scr. Anchoring on the pyproject.toml at the repo
# root makes all of them follow the checkout and survive this file changing
# depth.
PROJECT_BASE_DIR = next(
    q for q in pathlib.Path(__file__).resolve().parents
    if (q / "pyproject.toml").exists()
)


# ============================================================
# CONFIG
# ============================================================

# --- DSAS CSVs ---
DSAS_CSV_1978_1997 = str(PROJECT_BASE_DIR / "data" / "hatteras_init" / "5-scr" / "scr-dsas-1978-2019"
                      / "dsas_1978_1997_rates.csv")
DSAS_CSV_1997_2019 = str(PROJECT_BASE_DIR / "data" / "hatteras_init" / "5-scr" / "scr-dsas-1978-2019"
                      / "dsas_1997_2019_rates.csv")

DSAS_DOMAIN_COL = "domain_id"
DSAS_LRR_COL    = "MEAN_LRR"
DSAS_STD_COL    = "STD_LRR"

# --- CoastSat CSVs ---
COASTSAT_CSV_1978_1997 = str(PROJECT_BASE_DIR / "scripts" / "input_prep" / "5-scr" / "CoastSat" / "old_time_periods"
                          / "1978_1997" / "domain_lrr_summary.csv")
COASTSAT_CSV_1997_2019 = str(PROJECT_BASE_DIR / "scripts" / "input_prep" / "5-scr" / "CoastSat" / "old_time_periods"
                          / "1997_2019" / "domain_lrr_summary.csv")

CS_DOMAIN_COL = "domain_number"
CS_LRR_COL    = "mean_lrr"
CS_STD_COL    = "std_lrr"

# --- Domain range ---
DOMAIN_MIN = 1
DOMAIN_MAX = 90

# --- LOESS bandwidth (fraction of data used per local fit) ---
# 0.10 = ~9 domains  → more local, preserves more variation
# 0.167 = ~15 domains → recommended default (1.5km smoothing window)
# 0.20 = ~18 domains → smoother, loses finer spatial patterns
LOESS_FRAC = 0.15

# --- Southern boundary guard ---
# Domains 1..N are dropped from the SMOOTHED curves. Oregon Inlet dominates
# that zone, and LOESS is a local linear fit, so at the very edge it
# extrapolates: CoastSat 1978-1997 smooths to -6.21 m/yr at domain 1 where
# the raw value is -0.59 and the local raw spread is -7.01..-0.59. The
# smoothed value lands outside the data it claims to summarise.
#
# Same guard, same width as the hindcast's
# cascade_pipeline/coastsat_loess.py: LoessConfig(skip_southern_domains=10).
# DISPLAY ONLY, as it is there -- the LOESS still fits over all 90 domains,
# so the southern data still pulls the values just north of the cut; only the
# result is withheld. Set to 0 to show LOESS everywhere.
SKIP_SOUTHERN_DOMAINS = 10

# --- Town locations for reference lines ---
TOWNS = {
    "Buxton": 8,
    "Avon": 26,
    "Salvo": 69,
    "Waves": 74,
    "Rodanthe": 80,
}

# --- Output ---
# Products live under data/hatteras_init/<stage>/, beside every other
# input_prep stage's output; only the scripts live under scripts/.
OUTPUT_DIR = str(PROJECT_BASE_DIR / "data" / "hatteras_init" / "6-scr-smooth"
                 / "HAT_loess_dsas_vs_coastsat_output")

# ============================================================
# IMPORTS
# ============================================================
import os
import sys

# Windows consoles default to cp1252, which cannot encode the arrows and
# en-dashes in the closing summary -- the run died there after writing every
# figure. UTF-8 here, matching HAT_loess_method_comparison.py.
if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.transforms import blended_transform_factory
from scipy import stats
from statsmodels.nonparametric.smoothers_lowess import lowess
import warnings
warnings.filterwarnings("ignore")

os.makedirs(OUTPUT_DIR, exist_ok=True)

# ============================================================
# STYLE
# ============================================================
plt.rcParams.update({
    "font.family": "Arial",
    "axes.spines.top": False,
    "axes.spines.right": False,
    "axes.grid": True,
    "grid.alpha": 0.25,
    "grid.linestyle": ":",
})

# Color palette — matches your existing poster style
C_DSAS_1978   = "#1F4E79"   # dark blue
C_CS_1978     = "#5B9BD5"   # light blue (dashed)
C_DSAS_1997   = "#833C00"   # dark red
C_CS_1997     = "#F4A460"   # light red/tan (dashed)
C_NO_DSAS     = "#8C8C8C"   # hatch over domains where DSAS has no data
C_SKIP_ZONE   = "#6A8CAF"   # band over the domains where LOESS is suppressed

# ============================================================
# DATA LOADING
# ============================================================

def load_dsas(path, domain_col, lrr_col, std_col, period_label):
    df = pd.read_csv(path)
    df[domain_col] = pd.to_numeric(df[domain_col], errors="coerce")
    df[lrr_col]    = pd.to_numeric(df[lrr_col],    errors="coerce")
    df[std_col]    = pd.to_numeric(df[std_col],     errors="coerce")
    df = df.dropna(subset=[domain_col, lrr_col])
    df[domain_col] = df[domain_col].astype(int)
    df = df[(df[domain_col] >= DOMAIN_MIN) & (df[domain_col] <= DOMAIN_MAX)]
    df = df[[domain_col, lrr_col, std_col]].rename(columns={
        domain_col: "domain",
        lrr_col:    "dsas_lrr",
        std_col:    "dsas_std",
    }).sort_values("domain").reset_index(drop=True)
    print(f"  DSAS ({period_label}): {len(df)} domains  "
          f"range {df['dsas_lrr'].min():+.2f} to {df['dsas_lrr'].max():+.2f} m/yr")
    return df


def load_coastsat(path, domain_col, lrr_col, std_col, period_label):
    if path is None or not os.path.exists(path):
        print(f"  CoastSat ({period_label}): SKIPPED — not found: {path}")
        return None
    df = pd.read_csv(path)
    df[domain_col] = pd.to_numeric(df[domain_col], errors="coerce")
    df[lrr_col]    = pd.to_numeric(df[lrr_col],    errors="coerce")
    df[std_col]    = pd.to_numeric(df[std_col],     errors="coerce")
    df = df.dropna(subset=[domain_col, lrr_col])
    df[domain_col] = df[domain_col].astype(int)
    df = df[(df[domain_col] >= DOMAIN_MIN) & (df[domain_col] <= DOMAIN_MAX)]
    df = df[[domain_col, lrr_col, std_col]].rename(columns={
        domain_col: "domain",
        lrr_col:    "cs_lrr",
        std_col:    "cs_std",
    }).sort_values("domain").reset_index(drop=True)
    print(f"  CoastSat ({period_label}): {len(df)} domains  "
          f"range {df['cs_lrr'].min():+.2f} to {df['cs_lrr'].max():+.2f} m/yr")
    return df


def merge_datasets(dsas, coastsat):
    merged = dsas.merge(coastsat, on="domain", how="outer")
    merged["difference"] = merged["cs_lrr"] - merged["dsas_lrr"]
    return merged.sort_values("domain").reset_index(drop=True)


# ============================================================
# SMOOTHING
# ============================================================

def apply_loess(domains, values, frac=LOESS_FRAC):
    """Apply LOESS smoothing. Returns smoothed values at same domain positions."""
    valid = ~np.isnan(values)
    if valid.sum() < 5:
        return values.copy()
    smoothed = np.full_like(values, np.nan)
    result = lowess(values[valid], domains[valid], frac=frac, return_sorted=True)
    # Interpolate back to original domain positions
    smoothed[valid] = np.interp(domains[valid], result[:, 0], result[:, 1])
    return smoothed


def add_smoothed_columns(merged, frac=LOESS_FRAC):
    """Add LOESS-smoothed LRR columns to a merged dataframe."""
    d = merged["domain"].values.astype(float)
    merged = merged.copy()
    merged["dsas_lrr_smooth"] = apply_loess(d, merged["dsas_lrr"].values, frac)
    merged["cs_lrr_smooth"]   = apply_loess(d, merged["cs_lrr"].values,   frac)
    merged["diff_smooth"]     = merged["cs_lrr_smooth"] - merged["dsas_lrr_smooth"]
    return merged


# ============================================================
# STATS
# ============================================================

def regression_stats(x, y):
    mask = ~(np.isnan(x) | np.isnan(y))
    x, y = x[mask], y[mask]
    if len(x) < 3:
        return {}
    slope, intercept, r, p, _ = stats.linregress(x, y)
    rmse = np.sqrt(np.mean((y - x) ** 2))
    bias = np.mean(y - x)
    return dict(slope=slope, intercept=intercept, r2=r**2,
                rmse=rmse, bias=bias, p=p, n=len(x))


def stats_annotation(s, prefix=""):
    return (f"{prefix}R²={s['r2']:.2f}  "
            f"RMSE={s['rmse']:.2f} m/yr  "
            f"Bias={s['bias']:+.2f} m/yr  "
            f"n={s['n']}")


# ============================================================
# SHARED AXIS HELPERS
# ============================================================

# Town labels sit in this band, measured down from the top of the axes; the
# stats boxes start below it. Both are in axes fractions so they cannot drift
# into each other when ylim changes.
TOWN_LABEL_Y = 0.985
STATS_BOX_Y  = 0.88


def add_town_lines(ax):
    """Vertical reference line and label per town, pinned to the axes top."""
    trans = blended_transform_factory(ax.transData, ax.transAxes)
    for town, dom in TOWNS.items():
        ax.axvline(dom, color="0.55", lw=0.9, ls="--", alpha=0.7, zorder=0)
        ax.text(dom, TOWN_LABEL_Y, town, transform=trans,
                ha="center", va="top", fontsize=8, color="0.35",
                bbox=dict(boxstyle="round,pad=0.2", fc="white", ec="none", alpha=0.7))


def mask_southern_smoothed(m):
    """Blank the smoothed columns across the southern boundary zone.

    Applied after the fit, never before: the LOESS still sees every domain,
    matching how the hindcast splices this zone out. Raw columns are
    untouched, so the raw series still plots across the whole island.
    """
    if SKIP_SOUTHERN_DOMAINS <= 0:
        return m
    m = m.copy()
    zone = m["domain"] <= SKIP_SOUTHERN_DOMAINS
    for col in ("dsas_lrr_smooth", "cs_lrr_smooth", "diff_smooth"):
        if col in m.columns:
            m.loc[zone, col] = np.nan
    return m


def shade_boundary_zone(ax, label=True):
    """Mark the domains whose smoothed values were withheld."""
    if SKIP_SOUTHERN_DOMAINS <= 0:
        return
    ax.axvspan(DOMAIN_MIN - 0.5, SKIP_SOUTHERN_DOMAINS + 0.5,
               facecolor=C_SKIP_ZONE, alpha=0.10, lw=0.0, zorder=0,
               label="LOESS withheld (Oregon Inlet)" if label else None)


def shade_missing_dsas(ax, m, label=True):
    """Hatch the domains where DSAS has no data, so absence is not read as agreement.

    The 1997-2019 DSAS export covers domains 8-90; without this the DSAS
    series simply starts late while CoastSat runs the full span, which reads
    as the two sources agreeing rather than as one of them being absent.

    Runs are found from the data, so this stays correct if the DSAS coverage
    changes. Only the first run is labelled, to keep one legend entry.
    """
    missing = m["dsas_lrr"].isna().values
    if not missing.any():
        return
    domains = m["domain"].values

    runs, start = [], None
    for i, gap in enumerate(missing):
        if gap and start is None:
            start = i
        elif not gap and start is not None:
            runs.append((start, i - 1)); start = None
    if start is not None:
        runs.append((start, len(missing) - 1))

    trans = blended_transform_factory(ax.transData, ax.transAxes)
    for k, (i0, i1) in enumerate(runs):
        lo, hi = domains[i0] - 0.5, domains[i1] + 0.5
        ax.axvspan(lo, hi, facecolor=C_NO_DSAS, alpha=0.13, hatch="///",
                   edgecolor=C_NO_DSAS, lw=0.0, zorder=0,
                   label="No DSAS data" if (label and k == 0) else None)
        # Backed in white like the town labels: on the difference panels the
        # mean-difference line runs straight through this text otherwise.
        ax.text((lo + hi) / 2, 0.5,
                f"no DSAS\ndata\n(domains {domains[i0]}–{domains[i1]})",
                transform=trans, ha="center", va="center", fontsize=7.5,
                color="0.35", style="italic", zorder=5,
                bbox=dict(boxstyle="round,pad=0.25", fc="white", ec="none",
                          alpha=0.75))


def style_domain_axis(ax):
    ax.set_xlim(DOMAIN_MIN - 0.5, DOMAIN_MAX + 0.5)
    ax.set_xticks(range(DOMAIN_MIN, DOMAIN_MAX + 1, 10))
    ax.set_xticks(range(DOMAIN_MIN, DOMAIN_MAX + 1, 5), minor=True)
    ax.set_xlabel("CASCADE Model Domain (500 m alongshore)",
                  fontsize=11, fontweight="bold")
    ax.axhline(0, color="black", lw=1.1, ls="--", alpha=0.55)
    ax.annotate("Accretion ▲", xy=(DOMAIN_MAX + 0.5, 0.15),
                xycoords=("data", "axes fraction"),
                ha="right", va="bottom", fontsize=8, color="0.4")
    ax.annotate("Erosion ▼", xy=(DOMAIN_MAX + 0.5, 0.10),
                xycoords=("data", "axes fraction"),
                ha="right", va="top", fontsize=8, color="0.4")


# ============================================================
# FIGURE 1 — PRIMARY: Raw + LOESS overlay, both periods
# ============================================================

def plot_overview_smoothed(merged_1978, merged_1997, out_path):
    """
    2-panel figure: raw lines (faded) + LOESS overlay (bold).
    This is the recommended figure for collaborator review.
    """
    fig, axes = plt.subplots(2, 1, figsize=(16, 11), sharex=True)
    fig.suptitle("DSAS vs CoastSat Shoreline Change Rates — Hatteras Island, NC\n"
                 "Raw (faded) + LOESS smoothed (bold)",
                 fontsize=14, fontweight="bold", y=1.01)
    fig.text(0.5, -0.01, "These are the legacy DSAS period pair, not the hindcast's 1984-2004 / 2004-2024 split.",
             ha="center", fontsize=9, color="0.35", style="italic")

    configs = [
        (axes[0], merged_1978, "1978–1997",
         C_DSAS_1978, C_CS_1978),
        (axes[1], merged_1997, "1997–2019",
         C_DSAS_1997, C_CS_1997),
    ]

    for ax, merged, label, c_dsas, c_cs in configs:
        if merged is None:
            ax.text(0.5, 0.5, f"No data — {label}",
                    transform=ax.transAxes, ha="center", fontsize=12)
            continue

        m = add_smoothed_columns(merged)
        m = mask_southern_smoothed(m)
        d = m["domain"]

        # --- Raw (faded, thin) ---
        ax.plot(d, m["dsas_lrr"], color=c_dsas, lw=1.0, alpha=0.25,
                marker="o", ms=2.5, zorder=2)
        ax.plot(d, m["cs_lrr"],   color=c_cs,   lw=1.0, alpha=0.25,
                marker="s", ms=2.5, ls="--", zorder=2)

        # --- ±1 std shading (very faint) ---
        ax.fill_between(d, m["dsas_lrr"] - m["dsas_std"],
                           m["dsas_lrr"] + m["dsas_std"],
                        color=c_dsas, alpha=0.06, zorder=1)
        ax.fill_between(d, m["cs_lrr"] - m["cs_std"],
                           m["cs_lrr"] + m["cs_std"],
                        color=c_cs, alpha=0.06, zorder=1)

        # --- LOESS smoothed (bold) ---
        ax.plot(d, m["dsas_lrr_smooth"], color=c_dsas, lw=2.8,
                label=f"DSAS {label} (smoothed)", zorder=4)
        ax.plot(d, m["cs_lrr_smooth"],   color=c_cs,   lw=2.8, ls="--",
                label=f"CoastSat {label} (smoothed)", zorder=4)

        # Stats on smoothed values
        s_raw    = regression_stats(m["dsas_lrr"].values, m["cs_lrr"].values)
        s_smooth = regression_stats(m["dsas_lrr_smooth"].values,
                                    m["cs_lrr_smooth"].values)

        #txt = (f"Raw:     R²={s_raw['r2']:.2f}  RMSE={s_raw['rmse']:.2f} m/yr  "
               #f"Bias={s_raw['bias']:+.2f} m/yr\n"
               #f"Smoothed: R²={s_smooth['r2']:.2f}  RMSE={s_smooth['rmse']:.2f} m/yr  "
               #f"Bias={s_smooth['bias']:+.2f} m/yr")
        #ax.text(0.01, 0.97, txt, transform=ax.transAxes, fontsize=8.5,
                #va="top", family="monospace",
                #bbox=dict(boxstyle="round", fc="white", alpha=0.88, ec="0.7"))

        ax.set_ylabel("Shoreline Change Rate (m/yr)", fontsize=11, fontweight="bold")
        shade_boundary_zone(ax)
        shade_missing_dsas(ax, m)
        ax.legend(fontsize=9.5, framealpha=0.95, loc="lower right")
        style_domain_axis(ax)
        add_town_lines(ax)

    plt.tight_layout()
    plt.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close()
    print(f"  Saved: {os.path.basename(out_path)}")


# ============================================================
# FIGURE 2 — SMOOTHED ONLY (clean, for presentations)
# ============================================================

def plot_smoothed_only(merged_1978, merged_1997, out_path):
    """
    2-panel: LOESS smoothed lines only, no raw data.
    Cleanest version for presentations or dissertation figures.
    """
    fig, axes = plt.subplots(2, 1, figsize=(16, 10), sharex=True)
    fig.suptitle("DSAS vs CoastSat Shoreline Change Rates — Hatteras Island, NC\n"
                 f"LOESS smoothed (frac={LOESS_FRAC})",
                 fontsize=14, fontweight="bold", y=1.01)
    fig.text(0.5, -0.01, "These are the legacy DSAS period pair, not the hindcast's 1984-2004 / 2004-2024 split.",
             ha="center", fontsize=9, color="0.35", style="italic")

    configs = [
        (axes[0], merged_1978, "1978–1997", C_DSAS_1978, C_CS_1978),
        (axes[1], merged_1997, "1997–2019", C_DSAS_1997, C_CS_1997),
    ]

    for ax, merged, label, c_dsas, c_cs in configs:
        if merged is None:
            continue
        m = add_smoothed_columns(merged)
        m = mask_southern_smoothed(m)
        d = m["domain"]

        ax.plot(d, m["dsas_lrr_smooth"], color=c_dsas, lw=3.0,
                label=f"DSAS {label}")
        ax.plot(d, m["cs_lrr_smooth"],   color=c_cs,   lw=3.0, ls="--",
                label=f"CoastSat {label}")

        # Shaded difference between smoothed series
        ax.fill_between(d, m["dsas_lrr_smooth"], m["cs_lrr_smooth"],
                        where=(m["cs_lrr_smooth"] >= m["dsas_lrr_smooth"]),
                        color=c_cs, alpha=0.15, label="CoastSat > DSAS")
        ax.fill_between(d, m["dsas_lrr_smooth"], m["cs_lrr_smooth"],
                        where=(m["cs_lrr_smooth"] < m["dsas_lrr_smooth"]),
                        color=c_dsas, alpha=0.15, label="DSAS > CoastSat")

        s = regression_stats(m["dsas_lrr_smooth"].values, m["cs_lrr_smooth"].values)
        ax.text(0.01, STATS_BOX_Y, stats_annotation(s), transform=ax.transAxes,
                fontsize=9, va="top",
                bbox=dict(boxstyle="round", fc="white", alpha=0.88, ec="0.7"))

        ax.set_ylabel("Shoreline Change Rate (m/yr)", fontsize=11, fontweight="bold")
        shade_boundary_zone(ax)
        shade_missing_dsas(ax, m)
        ax.legend(fontsize=9.5, framealpha=0.95, loc="lower right", ncol=2)
        style_domain_axis(ax)
        add_town_lines(ax)

    plt.tight_layout()
    plt.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close()
    print(f"  Saved: {os.path.basename(out_path)}")


# ============================================================
# FIGURE 3 — SMOOTHING SENSITIVITY (frac comparison)
# ============================================================

def plot_smoothing_sensitivity(merged, period_label, out_path):
    """
    3-panel showing effect of different LOESS bandwidths.
    Helps collaborators understand smoothing choice.
    """
    fracs = [0.10, 0.15, 0.20]
    labels_frac = ["frac=0.10  (~9 domains)", 
                   "frac=0.15  (~13 domains) ← default",
                   "frac=0.20  (~18 domains)"]

    fig, axes = plt.subplots(3, 1, figsize=(16, 13), sharex=True, sharey=True)
    fig.suptitle(f"LOESS Smoothing Sensitivity — {period_label}\n"
                 f"Effect of bandwidth on DSAS vs CoastSat comparison",
                 fontsize=13, fontweight="bold", y=1.01)

    c_dsas = C_DSAS_1978 if "1978" in period_label else C_DSAS_1997
    c_cs   = C_CS_1978   if "1978" in period_label else C_CS_1997

    for ax, frac, flabel in zip(axes, fracs, labels_frac):
        m = add_smoothed_columns(merged, frac=frac)
        m = mask_southern_smoothed(m)
        d = m["domain"]

        # Raw faded
        ax.plot(d, m["dsas_lrr"], color=c_dsas, lw=0.8, alpha=0.2, zorder=1)
        ax.plot(d, m["cs_lrr"],   color=c_cs,   lw=0.8, alpha=0.2, ls="--", zorder=1)

        # Smoothed
        ax.plot(d, m["dsas_lrr_smooth"], color=c_dsas, lw=2.5,
                label=f"DSAS smoothed", zorder=3)
        ax.plot(d, m["cs_lrr_smooth"],   color=c_cs,   lw=2.5, ls="--",
                label=f"CoastSat smoothed", zorder=3)

        s = regression_stats(m["dsas_lrr_smooth"].values, m["cs_lrr_smooth"].values)
        ax.text(0.01, STATS_BOX_Y,
                f"{flabel}   |   " + stats_annotation(s),
                transform=ax.transAxes, fontsize=8.5, va="top",
                bbox=dict(boxstyle="round", fc="white", alpha=0.88, ec="0.7"))

        ax.set_ylabel("Rate (m/yr)", fontsize=10, fontweight="bold")
        shade_boundary_zone(ax)
        shade_missing_dsas(ax, m)
        ax.legend(fontsize=9, framealpha=0.95, loc="lower right")
        style_domain_axis(ax)
        add_town_lines(ax)

    plt.tight_layout()
    plt.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close()
    print(f"  Saved: {os.path.basename(out_path)}")


# ============================================================
# FIGURE 4 — COMBINED: Both periods on one panel (smoothed)
# ============================================================

def plot_combined_sources(merged_1978, merged_1997, out_path):
    """
    Single panel: all 4 smoothed series together.
    Good for seeing overall pattern and period differences simultaneously.
    """
    fig, ax = plt.subplots(figsize=(16, 6))

    for merged, label, c_dsas, c_cs, ls_cs in [
        (merged_1978, "1978–1997", C_DSAS_1978, C_CS_1978,   "--"),
        (merged_1997, "1997–2019", C_DSAS_1997, C_CS_1997, "--"),
    ]:
        if merged is None:
            continue
        m = add_smoothed_columns(merged)
        m = mask_southern_smoothed(m)
        d = m["domain"]
        ax.plot(d, m["dsas_lrr_smooth"], color=c_dsas, lw=2.5,
                label=f"DSAS {label}")
        ax.fill_between(d, m["dsas_lrr_smooth"] - m["dsas_std"] * 0.5,
                           m["dsas_lrr_smooth"] + m["dsas_std"] * 0.5,
                        color=c_dsas, alpha=0.08)
        ax.plot(d, m["cs_lrr_smooth"], color=c_cs, lw=2.5, ls=ls_cs,
                label=f"CoastSat {label}")

    ax.set_ylabel("Shoreline Change Rate (m/yr)", fontsize=12, fontweight="bold")
    ax.set_title("DSAS vs CoastSat — All Periods Combined (LOESS smoothed)",
                 fontsize=13, fontweight="bold", pad=12)
    ax.legend(fontsize=10, framealpha=0.95, ncol=2)
    style_domain_axis(ax)
    add_town_lines(ax)

    # Caption
    fig.text(0.5, -0.04,
             f"Rates: domain-averaged LRR smoothed with LOESS (frac={LOESS_FRAC}). "
             f"Shading = ±0.5 std of DSAS transects per domain.",
             ha="center", fontsize=8, color="0.4", style="italic")

    plt.tight_layout()
    plt.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close()
    print(f"  Saved: {os.path.basename(out_path)}")


# ============================================================
# FIGURE 5 — SCATTER of smoothed values
# ============================================================

def pearson_r(x, y):
    """Pearson r over the pairs where both series are present."""
    x, y = np.asarray(x, dtype=float), np.asarray(y, dtype=float)
    both = ~np.isnan(x) & ~np.isnan(y)
    if both.sum() < 3:
        return np.nan, int(both.sum())
    return float(np.corrcoef(x[both], y[both])[0, 1]), int(both.sum())


def plot_scatter_smoothed(merged, period_label, out_path):
    m = add_smoothed_columns(merged)
    m = mask_southern_smoothed(m)
    s_raw    = regression_stats(m["dsas_lrr"].values,        m["cs_lrr"].values)
    s_smooth = regression_stats(m["dsas_lrr_smooth"].values, m["cs_lrr_smooth"].values)

    c_dsas = C_DSAS_1978 if "1978" in period_label else C_DSAS_1997
    c_cs   = C_CS_1978   if "1978" in period_label else C_CS_1997

    r_raw, _ = pearson_r(m["dsas_lrr"].values, m["cs_lrr"].values)
    N_EFF_SPANS = 1.0 / LOESS_FRAC

    fig, axes = plt.subplots(1, 2, figsize=(14, 6))

    for ax, x_col, y_col, s, subtitle in [
        (axes[0], "dsas_lrr",        "cs_lrr",        s_raw,    "Raw values"),
        (axes[1], "dsas_lrr_smooth", "cs_lrr_smooth", s_smooth,
         f"LOESS smoothed (frac={LOESS_FRAC}) — spatial structure, not agreement"),
    ]:
        x = m[x_col].values
        y = m[y_col].values
        ax.scatter(x, y, c=m["domain"], cmap="viridis", s=40, alpha=0.75, zorder=3)

        lim_val = max(abs(np.nanmin([x, y])), abs(np.nanmax([x, y]))) * 1.1
        lim = (-lim_val, lim_val)
        ax.plot(lim, lim, "k--", lw=1.2, alpha=0.5, label="1:1 line", zorder=2)

        #if s and s.get('n', 0) > 2:
            #xfit = np.linspace(*lim, 100)
            #yfit = s["slope"] * xfit + s["intercept"]
            #ax.plot(xfit, yfit, color=c_dsas, lw=2.0, label="Best fit", zorder=4)
            #txt = (f"R²={s['r2']:.2f}\nRMSE={s['rmse']:.2f} m/yr\n"
                   #f"Bias={s['bias']:+.2f} m/yr\nSlope={s['slope']:.2f}")
            #ax.text(0.04, 0.96, txt, transform=ax.transAxes, fontsize=9,
                    #va="top", family="monospace",
                    #bbox=dict(boxstyle="round", fc="white", alpha=0.88, ec="0.7"))

        ax.set_xlim(lim); ax.set_ylim(lim)
        ax.set_aspect("equal")
        ax.set_xlabel("DSAS LRR (m/yr)", fontsize=11, fontweight="bold")
        ax.set_ylabel("CoastSat LRR (m/yr)", fontsize=11, fontweight="bold")
        ax.set_title(subtitle, fontsize=11)
        ax.legend(fontsize=9, framealpha=0.95, loc="lower right")

        # Both panels carry their r. The smoothed one also carries the raw r
        # and a count of independent spans, because its own r is not a
        # like-for-like improvement: LOESS runs along domain order, so
        # neighbouring smoothed values are built from overlapping windows and
        # are not independent samples. A window holds frac*n domains, so the
        # number of effectively independent spans is about 1/frac -- roughly
        # 6 here, not the 90 points plotted.
        r_here, n_here = pearson_r(x, y)
        if x_col.endswith("_smooth"):
            note = (f"r = {r_here:+.2f}  over ~{N_EFF_SPANS:.0f} independent spans\n"
                    f"NOT comparable to the raw panel: LOESS correlates\n"
                    f"neighbouring domains. Raw r = {r_raw:+.2f}.")
            fc = "#FFF4E5"
        else:
            note = f"r = {r_here:+.2f}   n = {n_here}"
            fc = "white"
        ax.text(0.04, 0.96, note, transform=ax.transAxes, fontsize=8.5,
                va="top", family="monospace", zorder=6,
                bbox=dict(boxstyle="round", fc=fc, alpha=0.92, ec="0.7"))

        sm = plt.cm.ScalarMappable(cmap="viridis",
             norm=plt.Normalize(vmin=DOMAIN_MIN, vmax=DOMAIN_MAX))
        sm.set_array([])
        plt.colorbar(sm, ax=ax, label="Domain number", shrink=0.75)

    fig.suptitle(f"DSAS vs CoastSat Scatter — {period_label}",
                 fontsize=13, fontweight="bold")
    plt.tight_layout()
    plt.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close()
    print(f"  Saved: {os.path.basename(out_path)}")


# ============================================================
# ORIGINAL PLOTS (kept from original script)
# ============================================================

def plot_line_comparison(merged, period_label, out_path):
    fig, ax = plt.subplots(figsize=(16, 6))
    c_dsas = C_DSAS_1978 if "1978" in period_label else C_DSAS_1997
    c_cs   = C_CS_1978   if "1978" in period_label else C_CS_1997

    ax.plot(merged["domain"], merged["dsas_lrr"], color=c_dsas, lw=2.5,
            marker="o", ms=4, label=f"DSAS LRR {period_label}")
    ax.fill_between(merged["domain"],
                    merged["dsas_lrr"] - merged["dsas_std"],
                    merged["dsas_lrr"] + merged["dsas_std"],
                    color=c_dsas, alpha=0.12)
    ax.plot(merged["domain"], merged["cs_lrr"], color=c_cs, lw=2.5,
            marker="s", ms=4, ls="--", label=f"CoastSat LRR {period_label}")
    ax.fill_between(merged["domain"],
                    merged["cs_lrr"] - merged["cs_std"],
                    merged["cs_lrr"] + merged["cs_std"],
                    color=c_cs, alpha=0.12)
    ax.set_ylabel("LRR (m/yr)", fontsize=12, fontweight="bold")
    ax.set_title(f"DSAS vs CoastSat LRR — {period_label}  (raw, ±1 std)",
                 fontsize=13, fontweight="bold")
    shade_boundary_zone(ax)
    shade_missing_dsas(ax, merged)
    ax.legend(fontsize=10, framealpha=0.95)
    style_domain_axis(ax)
    add_town_lines(ax)
    plt.tight_layout()
    plt.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close()
    print(f"  Saved: {os.path.basename(out_path)}")


def plot_difference(merged, period_label, out_path):
    m    = add_smoothed_columns(merged)
    m = mask_southern_smoothed(m)
    diff = m["difference"]
    diff_smooth = m["diff_smooth"]
    colors = ["#2166ac" if v >= 0 else "#b2182b" for v in diff]

    fig, ax = plt.subplots(figsize=(16, 5))
    ax.bar(m["domain"], diff, color=colors, edgecolor="none",
           width=0.85, alpha=0.55, label="Raw difference")
    ax.plot(m["domain"], diff_smooth, color="black", lw=2.5,
            label=f"LOESS smoothed (frac={LOESS_FRAC})")
    ax.axhline(0, color="black", lw=1.2)
    ax.axhline(diff.mean(), color="grey", lw=1.5, ls="--",
               label=f"Mean difference ({diff.mean():+.2f} m/yr)")

    ax.set_ylabel("CoastSat − DSAS (m/yr)", fontsize=12, fontweight="bold")
    ax.set_title(f"Dataset Difference: CoastSat minus DSAS — {period_label}\n"
                 f"Blue = CoastSat more accretional  |  Red = CoastSat more erosional",
                 fontsize=12, fontweight="bold")
    shade_boundary_zone(ax)
    shade_missing_dsas(ax, m)
    ax.legend(fontsize=10, framealpha=0.95)
    style_domain_axis(ax)
    add_town_lines(ax)
    plt.tight_layout()
    plt.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close()
    print(f"  Saved: {os.path.basename(out_path)}")


# ============================================================
# MAIN
# ============================================================

def main():
    print("=" * 65)
    print("DSAS vs CoastSat LRR Comparison — Hatteras Island (Smoothed)")
    print(f"LOESS bandwidth: frac={LOESS_FRAC}")
    print("=" * 65)

    # Load
    print("\nLoading data...")
    dsas_1978 = load_dsas(DSAS_CSV_1978_1997, DSAS_DOMAIN_COL,
                          DSAS_LRR_COL, DSAS_STD_COL, "1978–1997")
    dsas_1997 = load_dsas(DSAS_CSV_1997_2019, DSAS_DOMAIN_COL,
                          DSAS_LRR_COL, DSAS_STD_COL, "1997–2019")
    cs_1978   = load_coastsat(COASTSAT_CSV_1978_1997, CS_DOMAIN_COL,
                              CS_LRR_COL, CS_STD_COL, "1978–1997")
    cs_1997   = load_coastsat(COASTSAT_CSV_1997_2019, CS_DOMAIN_COL,
                              CS_LRR_COL, CS_STD_COL, "1997–2019")

    # Merge
    merged_1978 = merge_datasets(dsas_1978, cs_1978) if cs_1978 is not None else None
    merged_1997 = merge_datasets(dsas_1997, cs_1997) if cs_1997 is not None else None

    # Stats summary
    #print("\n--- Agreement Statistics ---")
    #for merged, label in [(merged_1978, "1978–1997"), (merged_1997, "1997–2019")]:
        #if merged is None:
            #continue
        #m  = add_smoothed_columns(merged)
        #sr = regression_stats(m["dsas_lrr"].values,        m["cs_lrr"].values)
        #ss = regression_stats(m["dsas_lrr_smooth"].values, m["cs_lrr_smooth"].values)
        #print(f"\n  {label}")
        #print(f"    Raw     — R²={sr['r2']:.3f}  RMSE={sr['rmse']:.3f}  "
              #f"Bias={sr['bias']:+.3f}  Slope={sr['slope']:.3f}")
        #print(f"    Smoothed — R²={ss['r2']:.3f}  RMSE={ss['rmse']:.3f}  "
              #f"Bias={ss['bias']:+.3f}  Slope={ss['slope']:.3f}")

    # Generate all figures
    print("\nGenerating figures...")

    # NEW smoothed figures
    plot_overview_smoothed(merged_1978, merged_1997,
        os.path.join(OUTPUT_DIR, "overview_smoothed.png"))

    plot_smoothed_only(merged_1978, merged_1997,
        os.path.join(OUTPUT_DIR, "smoothed_only_comparison.png"))

    plot_combined_sources(merged_1978, merged_1997,
        os.path.join(OUTPUT_DIR, "combined_sources.png"))

    if merged_1978 is not None:
        plot_smoothing_sensitivity(merged_1978, "1978–1997",
            os.path.join(OUTPUT_DIR, "smoothing_sensitivity_1978_1997.png"))
        plot_scatter_smoothed(merged_1978, "1978–1997",
            os.path.join(OUTPUT_DIR, "scatter_smoothed_1978_1997.png"))

    if merged_1997 is not None:
        plot_smoothing_sensitivity(merged_1997, "1997–2019",
            os.path.join(OUTPUT_DIR, "smoothing_sensitivity_1997_2019.png"))
        plot_scatter_smoothed(merged_1997, "1997–2019",
            os.path.join(OUTPUT_DIR, "scatter_smoothed_1997_2019.png"))

    # Original figures (updated with new styling)
    if merged_1978 is not None:
        plot_line_comparison(merged_1978, "1978–1997",
            os.path.join(OUTPUT_DIR, "line_comparison_1978_1997.png"))
        plot_difference(merged_1978, "1978–1997",
            os.path.join(OUTPUT_DIR, "difference_1978_1997.png"))

    if merged_1997 is not None:
        plot_line_comparison(merged_1997, "1997–2019",
            os.path.join(OUTPUT_DIR, "line_comparison_1997_2019.png"))
        plot_difference(merged_1997, "1997–2019",
            os.path.join(OUTPUT_DIR, "difference_1997_2019.png"))

    # Save comparison table
    table_parts = []
    for merged, label in [(merged_1978, "1978_1997"), (merged_1997, "1997_2019")]:
        if merged is None:
            continue
        m = add_smoothed_columns(merged)
        # Same guard as the figures, so the table cannot be read as evidence
        # the plots withheld. The three *_smooth columns are therefore EMPTY
        # across domains 1..SKIP_SOUTHERN_DOMAINS; the raw columns are not,
        # and still cover the whole island.
        m = mask_southern_smoothed(m)
        cols = ["domain", "dsas_lrr", "cs_lrr", "dsas_std", "cs_std",
                "difference", "dsas_lrr_smooth", "cs_lrr_smooth", "diff_smooth"]
        t = m[cols].copy()
        t.columns = ["domain"] + [f"{c}_{label}" for c in cols if c != "domain"]
        table_parts.append(t)

    if table_parts:
        from functools import reduce
        table = reduce(lambda a, b: a.merge(b, on="domain", how="outer"), table_parts)
        table = table.sort_values("domain")
        out_csv = os.path.join(OUTPUT_DIR, "comparison_table_smoothed.csv")
        table.to_csv(out_csv, index=False)
        print(f"  Saved: comparison_table_smoothed.csv")

    print("\n" + "=" * 65)
    print("Done! Key figures for collaborator review:")
    print("  overview_smoothed.png        ← PRIMARY: raw + smoothed overlay")
    print("  smoothed_only_comparison.png ← clean version for presentations")
    print("  smoothing_sensitivity_*.png  ← shows effect of bandwidth choice")
    print("  combined_sources.png         ← all 4 series on one panel")
    print("  scatter_smoothed_*.png       ← raw vs smoothed scatter comparison")
    print("=" * 65)


if __name__ == "__main__":
    main()
