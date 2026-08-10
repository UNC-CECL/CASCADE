"""
HAT_shoreline_trajectory_classification.py
===========================================
Classifies shoreline trajectory stability across Hatteras Island domains
using CoastSat transect time-series over three periods:
    - Period 1 (calibration): 1984–2004
    - Period 2 (validation):  2004–2024
    - Full record:            1984–2024

For each domain and period, computes:
    - Linear rate of change (LRR) via OLS on annual median chainage
    - Sign consistency (fraction of year-on-year steps in dominant direction)
    - Interannual variability (std dev of annual positions around trend)

Classification scheme (applied per domain, per period pair):
    Persistent Erosion    — both periods erosional  (LRR < -THRESHOLD)
    Persistent Accretion  — both periods accretional (LRR > +THRESHOLD)
    Persistently Stable   — both periods within ±THRESHOLD
    Switching: Acc→Ero    — Period 1 accretional, Period 2 erosional
    Switching: Ero→Acc    — Period 1 erosional,   Period 2 accretional
    Transitioning         — one period stable, one not

Outputs
-------
1.  Along-island classification bar chart (domain × period)
2.  Period 1 vs Period 2 LRR scatter plot (per domain, coloured by class)
3.  Hovmöller heatmap (domain × year annual deviation) with classification overlay
4.  CSV summary table of all metrics

Usage
-----
    python HAT_shoreline_trajectory_classification.py

Dependencies
------------
    pip install pandas numpy matplotlib scipy tqdm
"""

import os
import glob

import numpy as np
import pandas as pd
from scipy import stats
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.patches as mpatches
import matplotlib.gridspec as gridspec
from tqdm import tqdm

# ============================================================
# CONFIG
# ============================================================

ROOT_DATA_DIR = r"/scripts/input_prep/CoastSat/coastsat_timeseries"
LOOKUP_CSV    = r"C:\Users\hanna\PycharmProjects\CASCADE\scripts\input_preperation\CoastSat\transect_domain_lookup.csv"
OUTPUT_DIR    = r"/scripts/input_prep/shoreline_change_patterns/old_drafts/comparison"

SITE_FILTER = "usa_NC"

# Analysis periods
PERIOD_1_START, PERIOD_1_END = 1984, 2004
PERIOD_2_START, PERIOD_2_END = 2004, 2024
FULL_START,     FULL_END     = 1984, 2024

# Stability threshold (m/yr) — LRR within ±THRESHOLD = stable
THRESHOLD = 1.0

# High-variability flag: std dev of annual chainage deviations above this = episodic
VARIABILITY_THRESHOLD = 15.0   # metres

# Minimum observations per transect per year to count
MIN_OBS_PER_YEAR = 2

# Minimum fraction of years with valid data for a domain LRR to be trusted
MIN_YEAR_FRACTION = 0.50

# Gaussian smoothing for Hovmöller (sigma in domain units, 0 = off)
HOVM_SMOOTH_SIGMA = 1.5

NUM_REAL_DOMAINS = 90

# ── Publication annotation colours ───────────────────────────────────────────
ANN_TOWN_SPANS = {
    "Buxton":      (7,  8),
    "Avon":        (21, 31),
    "Tri-Village": (68, 83),
}
ANN_VILLAGE_LINES = {"Salvo": 69, "Waves": 74, "Rodanthe": 80}
ANN_PIERS         = {"Avon Pier": 26, "Rodanthe Pier": 79}
ANN_GROINS        = {"Buxton Groin": 5.5}
ANN_WIMBLE_SHOALS = (60, 74)

ANN_C_TOWN_SPAN    = "#90AFC5"
ANN_C_WIMBLE       = "#E0A800"
ANN_C_VILLAGE_LINE = "0.40"
ANN_C_PIER         = "#1565C0"
ANN_C_GROIN        = "#B71C1C"

# ── Trajectory classification colours ────────────────────────────────────────
CLASS_COLORS = {
    "Persistent Erosion":    "#b2182b",   # strong red
    "Persistent Accretion":  "#2166ac",   # strong blue
    "Persistently Stable":   "#4dac26",   # green
    "Switching: Acc→Ero":    "#f4a582",   # light red-orange
    "Switching: Ero→Acc":    "#92c5de",   # light blue
    "Transitioning":         "#d9d9d9",   # light grey
    "Insufficient Data":     "#ffffff",   # white
}

# ============================================================
# HELPERS
# ============================================================

def find_csvs(root_dir, site_filter):
    csv_map = {}
    for folder in os.listdir(root_dir):
        if not folder.startswith(site_filter):
            continue
        folder_path = os.path.join(root_dir, folder)
        if not os.path.isdir(folder_path):
            continue
        for csv_file in glob.glob(os.path.join(folder_path, "*.csv")):
            tid = os.path.splitext(os.path.basename(csv_file))[0]
            csv_map[tid] = csv_file
    return csv_map


def load_transect(csv_path):
    try:
        df = pd.read_csv(csv_path)
        df.columns = [c.strip() for c in df.columns]
        date_col  = next((c for c in df.columns if c.lower().startswith("date")), None)
        chain_col = next((c for c in df.columns
                          if "chainage" in c.lower() or "cross_dist" in c.lower()), None)
        if date_col is None or chain_col is None:
            return None
        df[date_col] = pd.to_datetime(df[date_col], utc=True).dt.tz_localize(None)
        df = df.dropna(subset=[date_col, chain_col]).sort_values(date_col)
        return df.set_index(date_col)[chain_col]
    except Exception:
        return None


def annual_median(ts, yr_start, yr_end, min_obs=2):
    """Return Series of annual median chainage, index = year, for yr_start..yr_end."""
    years, vals = [], []
    for yr in range(yr_start, yr_end + 1):
        d = ts[ts.index.year == yr]
        if len(d) >= min_obs:
            years.append(yr)
            vals.append(float(d.median()))
    if not years:
        return pd.Series(dtype=float)
    return pd.Series(vals, index=years)


def compute_lrr(ann_series):
    """OLS linear rate of change (m/yr) from annual median series.
    Returns (lrr, r2, n_years) or (nan, nan, 0) if insufficient data."""
    s = ann_series.dropna()
    if len(s) < 4:
        return np.nan, np.nan, len(s)
    slope, intercept, r, p, se = stats.linregress(s.index.values, s.values)
    return slope, r**2, len(s)


def sign_consistency(ann_series):
    """Fraction of year-on-year steps in the dominant direction (0.5–1.0).
    1.0 = monotonic, 0.5 = random walk."""
    s = ann_series.dropna()
    if len(s) < 3:
        return np.nan
    diffs = np.diff(s.values)
    diffs = diffs[diffs != 0]
    if len(diffs) == 0:
        return np.nan
    dominant = np.sign(diffs).mean()   # +1 all accreting, -1 all eroding
    return abs(dominant)               # consistency regardless of direction


def classify(lrr1, lrr2, thr=THRESHOLD):
    """Classify trajectory from Period 1 and Period 2 LRRs."""
    if np.isnan(lrr1) or np.isnan(lrr2):
        return "Insufficient Data"
    e1 = lrr1 < -thr
    a1 = lrr1 >  thr
    s1 = not e1 and not a1
    e2 = lrr2 < -thr
    a2 = lrr2 >  thr
    s2 = not e2 and not a2
    if   e1 and e2:  return "Persistent Erosion"
    elif a1 and a2:  return "Persistent Accretion"
    elif s1 and s2:  return "Persistently Stable"
    elif a1 and e2:  return "Switching: Acc→Ero"
    elif e1 and a2:  return "Switching: Ero→Acc"
    else:            return "Transitioning"


def gaussian_smooth_1d(arr, sigma):
    if sigma <= 0:
        return arr
    from scipy.ndimage import gaussian_filter1d
    mask   = np.isnan(arr)
    filled = arr.copy(); filled[mask] = 0.0
    sm     = gaussian_filter1d(filled, sigma=sigma)
    wt     = gaussian_filter1d((~mask).astype(float), sigma=sigma)
    with np.errstate(invalid="ignore"):
        return np.where(wt > 0.05, sm / wt, np.nan)


def add_geo_annotations(ax, orientation="horizontal", label_side="top",
                         secondary_axis=True):
    """
    Add community spans, Wimble Shoals, piers, groin, village lines.
    orientation: 'horizontal' (domain on x-axis) or 'vertical' (domain on y-axis).
    """
    if orientation == "horizontal":
        span_fn   = ax.axvspan
        line_fn   = ax.axvline
        lim_fn    = ax.get_xlim
        text_x    = lambda mid: mid
        text_y_top    = lambda ymax, frac: ymax - 0.01 * abs(ymax - ax.get_ylim()[0])
        text_va   = "top"
    else:
        span_fn   = ax.axhspan
        line_fn   = ax.axhline
        lim_fn    = ax.get_ylim

    ymin, ymax = ax.get_ylim()
    xmin, xmax = ax.get_xlim()
    yspan = abs(ymax - ymin)

    # Wimble Shoals
    span_fn(ANN_WIMBLE_SHOALS[0] - 0.5, ANN_WIMBLE_SHOALS[1] + 0.5,
            color=ANN_C_WIMBLE, alpha=0.15, zorder=0)

    # Community spans
    for name, (d0, d1) in ANN_TOWN_SPANS.items():
        span_fn(d0 - 0.5, d1 + 0.5, color=ANN_C_TOWN_SPAN, alpha=0.20, zorder=0)
        if orientation == "horizontal":
            mid = (d0 + d1) / 2
            ax.text(mid, ymax - 0.01 * yspan, name,
                    ha="center", va="top", fontsize=6.5,
                    color=ANN_C_TOWN_SPAN, fontweight="bold", zorder=3)

    # Village lines
    for vname, vdom in ANN_VILLAGE_LINES.items():
        line_fn(vdom, color=ANN_C_VILLAGE_LINE, lw=0.6, ls=(0, (3, 4)), zorder=2)

    # Piers
    for pname, pdom in ANN_PIERS.items():
        line_fn(pdom, color=ANN_C_PIER, lw=0.9, ls="--", zorder=2)
        if orientation == "horizontal":
            ax.text(pdom + 0.3, ymin + 0.72 * yspan, pname,
                    rotation=90, va="bottom", ha="left",
                    fontsize=6, color=ANN_C_PIER, zorder=3)

    # Groins
    for gname, gdom in ANN_GROINS.items():
        line_fn(gdom, color=ANN_C_GROIN, lw=0.9, ls="--", zorder=2)
        if orientation == "horizontal":
            ax.text(gdom + 0.3, ymin + 0.62 * yspan, gname,
                    rotation=90, va="bottom", ha="left",
                    fontsize=6, color=ANN_C_GROIN, zorder=3)


# ============================================================
# ANALYSIS
# ============================================================

def compute_domain_metrics(ts_dict, transect_order, domain_per_transect):
    """
    For every domain, compute LRR, sign consistency, and variability
    for Period 1, Period 2, and the full record.
    Returns a DataFrame indexed by domain number.
    """
    domains = sorted(set(domain_per_transect))
    rows = []

    for dom in tqdm(domains, desc="  Computing domain metrics"):
        tids = [t for t, d in zip(transect_order, domain_per_transect) if d == dom]

        # Collect all transect annual series and aggregate to domain median per year
        transect_rows_p1, transect_rows_p2, transect_rows_full = [], [], []
        transect_lrr_p1, transect_lrr_p2 = [], []

        for tid in tids:
            ts = ts_dict.get(tid)
            if ts is None:
                continue
            ann_p1   = annual_median(ts, PERIOD_1_START, PERIOD_1_END)
            ann_p2   = annual_median(ts, PERIOD_2_START, PERIOD_2_END)
            ann_full = annual_median(ts, FULL_START, FULL_END)

            lrr1, r2_1, n1 = compute_lrr(ann_p1)
            lrr2, r2_2, n2 = compute_lrr(ann_p2)

            transect_lrr_p1.append(lrr1)
            transect_lrr_p2.append(lrr2)
            transect_rows_p1.append(ann_p1)
            transect_rows_p2.append(ann_p2)
            transect_rows_full.append(ann_full)

        if not transect_rows_full:
            rows.append({"domain": dom})
            continue

        # Domain-level: median across transects for each year
        def domain_annual(transect_series_list, yr_start, yr_end):
            years = list(range(yr_start, yr_end + 1))
            out = {}
            for yr in years:
                vals = [s[yr] for s in transect_series_list if yr in s.index]
                out[yr] = float(np.nanmedian(vals)) if vals else np.nan
            return pd.Series(out).dropna()

        dom_ann_p1   = domain_annual(transect_rows_p1,   PERIOD_1_START, PERIOD_1_END)
        dom_ann_p2   = domain_annual(transect_rows_p2,   PERIOD_2_START, PERIOD_2_END)
        dom_ann_full = domain_annual(transect_rows_full,  FULL_START,     FULL_END)

        # Domain-level LRR
        dom_lrr1, dom_r2_1, dom_n1 = compute_lrr(dom_ann_p1)
        dom_lrr2, dom_r2_2, dom_n2 = compute_lrr(dom_ann_p2)
        dom_lrr_full, dom_r2_full, dom_n_full = compute_lrr(dom_ann_full)

        # Sign consistency
        dom_sc1 = sign_consistency(dom_ann_p1)
        dom_sc2 = sign_consistency(dom_ann_p2)

        # Variability (std dev of annual values around linear trend)
        def detrended_std(ann_s):
            s = ann_s.dropna()
            if len(s) < 4:
                return np.nan
            slope, intercept, *_ = stats.linregress(s.index.values, s.values)
            resid = s.values - (slope * s.index.values + intercept)
            return float(np.std(resid))

        dom_var_p1   = detrended_std(dom_ann_p1)
        dom_var_p2   = detrended_std(dom_ann_p2)
        dom_var_full = detrended_std(dom_ann_full)

        # Trajectory classification
        traj_class = classify(dom_lrr1, dom_lrr2)
        high_var   = (dom_var_full > VARIABILITY_THRESHOLD
                      if not np.isnan(dom_var_full) else False)

        # Transect-level median LRR (for comparison)
        txn_lrr1_med = float(np.nanmedian(transect_lrr_p1)) if transect_lrr_p1 else np.nan
        txn_lrr2_med = float(np.nanmedian(transect_lrr_p2)) if transect_lrr_p2 else np.nan

        rows.append({
            "domain":         dom,
            # Domain-averaged metrics
            "dom_lrr_p1":     dom_lrr1,
            "dom_lrr_p2":     dom_lrr2,
            "dom_lrr_full":   dom_lrr_full,
            "dom_r2_p1":      dom_r2_1,
            "dom_r2_p2":      dom_r2_2,
            "dom_r2_full":    dom_r2_full,
            "dom_n_p1":       dom_n1,
            "dom_n_p2":       dom_n2,
            "dom_n_full":     dom_n_full,
            "dom_sc_p1":      dom_sc1,
            "dom_sc_p2":      dom_sc2,
            "dom_var_p1":     dom_var_p1,
            "dom_var_p2":     dom_var_p2,
            "dom_var_full":   dom_var_full,
            # Transect-median metrics
            "txn_lrr_p1":     txn_lrr1_med,
            "txn_lrr_p2":     txn_lrr2_med,
            "txn_class":      classify(txn_lrr1_med, txn_lrr2_med),
            # Classification
            "trajectory":     traj_class,
            "high_var":       high_var,
            # Full annual series for Hovmöller (stored as dict)
            "_ann_full":      dom_ann_full.to_dict(),
        })

    df = pd.DataFrame(rows).set_index("domain")
    return df


# ============================================================
# FIGURE 1 — Along-island classification bar chart
# ============================================================

def plot_classification_bar(metrics, out_path):
    fig, axes = plt.subplots(4, 1, figsize=(18, 12), sharex=True,
                             gridspec_kw={"height_ratios": [3, 3, 1, 1]})

    domains = metrics.index.values
    x       = domains  # plot directly in domain space

    # ── Panel 1: Period 1 LRR ────────────────────────────────────────────────
    ax = axes[0]
    lrr1 = metrics["dom_lrr_p1"].values
    txn1 = metrics["txn_lrr_p1"].values
    ax.bar(x, lrr1, width=0.7, color=[
        "#b2182b" if v < -THRESHOLD else "#2166ac" if v > THRESHOLD else "#4dac26"
        if not np.isnan(v) else "#dddddd" for v in lrr1
    ], alpha=0.85, zorder=2, label="Domain median LRR")
    ax.plot(x, txn1, "o", ms=2.5, color="#333333", alpha=0.6, zorder=3,
            label="Transect median LRR")
    ax.axhline(0, color="black", lw=0.8, zorder=1)
    ax.axhline( THRESHOLD, color="#888888", lw=0.7, ls="--", zorder=1)
    ax.axhline(-THRESHOLD, color="#888888", lw=0.7, ls="--", zorder=1)
    ax.set_ylabel("LRR  (m/yr)\n▲ accretion   ▼ erosion", fontsize=8)
    ax.set_title("Period 1  (1984–2004)", fontsize=9, loc="left", pad=3)
    ax.set_ylim(np.nanmin([lrr1, txn1]) * 1.3, np.nanmax([lrr1, txn1]) * 1.3)
    ax.tick_params(labelsize=7)
    ax.legend(fontsize=6, loc="upper right")
    add_geo_annotations(ax, "horizontal")

    # ── Panel 2: Period 2 LRR ────────────────────────────────────────────────
    ax = axes[1]
    lrr2 = metrics["dom_lrr_p2"].values
    txn2 = metrics["txn_lrr_p2"].values
    ax.bar(x, lrr2, width=0.7, color=[
        "#b2182b" if v < -THRESHOLD else "#2166ac" if v > THRESHOLD else "#4dac26"
        if not np.isnan(v) else "#dddddd" for v in lrr2
    ], alpha=0.85, zorder=2)
    ax.plot(x, txn2, "o", ms=2.5, color="#333333", alpha=0.6, zorder=3)
    ax.axhline(0, color="black", lw=0.8, zorder=1)
    ax.axhline( THRESHOLD, color="#888888", lw=0.7, ls="--", zorder=1)
    ax.axhline(-THRESHOLD, color="#888888", lw=0.7, ls="--", zorder=1)
    ax.set_ylabel("LRR  (m/yr)\n▲ accretion   ▼ erosion", fontsize=8)
    ax.set_title("Period 2  (2004–2024)", fontsize=9, loc="left", pad=3)
    ax.set_ylim(np.nanmin([lrr2, txn2]) * 1.3, np.nanmax([lrr2, txn2]) * 1.3)
    ax.tick_params(labelsize=7)
    add_geo_annotations(ax, "horizontal")

    # ── Panel 3: Trajectory classification strip ──────────────────────────────
    ax = axes[2]
    for dom in domains:
        cl = metrics.loc[dom, "trajectory"]
        col = CLASS_COLORS.get(cl, "#ffffff")
        ax.bar(dom, 1, width=0.95, bottom=0, color=col,
               edgecolor="white", linewidth=0.3, zorder=2)
    ax.set_yticks([])
    ax.set_ylabel("Trajectory\nclass", fontsize=7)
    ax.set_ylim(0, 1)
    ax.tick_params(labelsize=7)
    # Wimble + community spans (no text labels, panel too thin)
    ax.axvspan(ANN_WIMBLE_SHOALS[0]-0.5, ANN_WIMBLE_SHOALS[1]+0.5,
               color=ANN_C_WIMBLE, alpha=0.12, zorder=0)
    for name, (d0, d1) in ANN_TOWN_SPANS.items():
        ax.axvspan(d0-0.5, d1+0.5, color=ANN_C_TOWN_SPAN, alpha=0.18, zorder=0)

    # ── Panel 4: Variability strip ────────────────────────────────────────────
    ax = axes[3]
    var_vals = metrics["dom_var_full"].values
    bar_colors = ["#d73027" if v > VARIABILITY_THRESHOLD else "#fee090"
                  if v > VARIABILITY_THRESHOLD * 0.5 else "#74add1"
                  for v in var_vals]
    ax.bar(x, np.nan_to_num(var_vals), width=0.7, color=bar_colors,
           alpha=0.85, zorder=2)
    ax.axhline(VARIABILITY_THRESHOLD, color="#d73027", lw=0.8, ls="--", zorder=1)
    ax.set_ylabel("Variability\n(m std dev)", fontsize=7)
    ax.set_xlabel("CASCADE domain  (1 = Buxton  →  90 = Rodanthe)", fontsize=8)
    ax.tick_params(labelsize=7)
    ax.set_xlim(0.5, NUM_REAL_DOMAINS + 0.5)
    add_geo_annotations(ax, "horizontal")

    # Legend for classification colours
    patches = [mpatches.Patch(color=v, label=k)
               for k, v in CLASS_COLORS.items() if k != "Insufficient Data"]
    axes[2].legend(handles=patches, fontsize=6, loc="upper right",
                   ncol=3, framealpha=0.9, title="Trajectory class", title_fontsize=6)

    fig.suptitle(
        "Hatteras Island  |  Shoreline trajectory classification by domain\n"
        f"CoastSat 1984–2024  ·  stability threshold ±{THRESHOLD} m/yr",
        fontsize=11, fontweight="bold", y=0.98)
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    fig.savefig(out_path, dpi=180, bbox_inches="tight")
    plt.close(fig)
    print(f"  Figure 1 saved → {out_path}")


# ============================================================
# FIGURE 2 — P1 vs P2 LRR scatter
# ============================================================

def plot_lrr_scatter(metrics, out_path):
    fig, axes = plt.subplots(1, 2, figsize=(13, 6), sharey=False)

    for ax, lrr_col_x, lrr_col_y, label_x, label_y, title_suffix in [
        (axes[0], "dom_lrr_p1", "dom_lrr_p2",
         "Period 1 LRR  (1984–2004)  m/yr",
         "Period 2 LRR  (2004–2024)  m/yr",
         "Domain-averaged"),
        (axes[1], "txn_lrr_p1", "txn_lrr_p2",
         "Period 1 LRR  (1984–2004)  m/yr",
         "Period 2 LRR  (2004–2024)  m/yr",
         "Transect-median"),
    ]:
        x_vals = metrics[lrr_col_x].values
        y_vals = metrics[lrr_col_y].values
        classes = metrics["trajectory"].values if "dom" in lrr_col_x \
                  else metrics["txn_class"].values

        colors = [CLASS_COLORS.get(c, "#aaaaaa") for c in classes]

        sc = ax.scatter(x_vals, y_vals, c=colors, s=40, zorder=3,
                        edgecolors="white", linewidths=0.5, alpha=0.9)

        # Annotate a few notable domains
        for dom in metrics.index:
            xi = metrics.loc[dom, lrr_col_x]
            yi = metrics.loc[dom, lrr_col_y]
            if np.isnan(xi) or np.isnan(yi):
                continue
            if abs(xi) > 2.5 or abs(yi) > 2.5:
                ax.text(xi + 0.05, yi + 0.05, str(dom),
                        fontsize=5.5, color="#333333", zorder=4)

        # Quadrant lines
        ax.axhline(0, color="black", lw=0.7, zorder=1)
        ax.axvline(0, color="black", lw=0.7, zorder=1)
        ax.axhline( THRESHOLD, color="#888888", lw=0.6, ls="--", zorder=1)
        ax.axhline(-THRESHOLD, color="#888888", lw=0.6, ls="--", zorder=1)
        ax.axvline( THRESHOLD, color="#888888", lw=0.6, ls="--", zorder=1)
        ax.axvline(-THRESHOLD, color="#888888", lw=0.6, ls="--", zorder=1)

        # 1:1 line — where P1=P2 (perfectly persistent)
        lims = [min(np.nanmin(x_vals), np.nanmin(y_vals)) * 1.1,
                max(np.nanmax(x_vals), np.nanmax(y_vals)) * 1.1]
        ax.plot(lims, lims, color="#555555", lw=0.8, ls=":", zorder=1,
                label="1:1 line (P1 = P2)")

        # Quadrant labels
        for tx, ty, label in [
            (lims[0]*0.7,  lims[1]*0.7,  "Persistent\nErosion"),
            (lims[1]*0.7,  lims[1]*0.7,  "Switching\nEro→Acc"),
            (lims[0]*0.7,  lims[0]*0.7,  "Switching\nAcc→Ero"),
            (lims[1]*0.7,  lims[0]*0.7,  "Persistent\nAccretion"),
        ]:
            ax.text(tx, ty, label, ha="center", va="center",
                    fontsize=6, color="#aaaaaa", style="italic")

        ax.set_xlim(lims); ax.set_ylim(lims)
        ax.set_xlabel(label_x, fontsize=8)
        ax.set_ylabel(label_y, fontsize=8)
        ax.set_title(title_suffix, fontsize=9, fontweight="bold")
        ax.set_aspect("equal")
        ax.tick_params(labelsize=7)
        ax.spines[["top", "right"]].set_visible(False)

    # Shared legend
    patches = [mpatches.Patch(color=v, label=k)
               for k, v in CLASS_COLORS.items() if k != "Insufficient Data"]
    fig.legend(handles=patches, fontsize=7, loc="lower center",
               ncol=3, framealpha=0.9, title="Trajectory class",
               title_fontsize=7, bbox_to_anchor=(0.5, -0.02))

    fig.suptitle(
        "Hatteras Island  |  Period 1 vs Period 2 linear rate of change\n"
        "Quadrant position shows trajectory persistence or switching",
        fontsize=11, fontweight="bold")
    fig.tight_layout(rect=[0, 0.08, 1, 0.95])
    fig.savefig(out_path, dpi=180, bbox_inches="tight")
    plt.close(fig)
    print(f"  Figure 2 saved → {out_path}")


# ============================================================
# FIGURE 3 — Hovmöller heatmap with classification overlay
# ============================================================

def plot_hovmoller(metrics, out_path):
    # Build domain × year matrix from stored annual series
    years   = list(range(FULL_START, FULL_END + 1))
    domains = sorted(metrics.index)

    matrix = np.full((len(domains), len(years)), np.nan)
    for di, dom in enumerate(domains):
        ann_dict = metrics.loc[dom, "_ann_full"]
        if not isinstance(ann_dict, dict):
            continue
        # Subtract 1984 value as within-domain baseline for display
        base = ann_dict.get(1984, np.nan)
        for yi, yr in enumerate(years):
            if yr in ann_dict:
                matrix[di, yi] = ann_dict[yr] - base if not np.isnan(base) else ann_dict[yr]

    # Smooth lightly along domain axis
    if HOVM_SMOOTH_SIGMA > 0:
        from scipy.ndimage import gaussian_filter1d
        for yi in range(len(years)):
            col = matrix[:, yi]
            mask = np.isnan(col)
            if mask.all():
                continue
            filled = col.copy(); filled[mask] = 0.0
            sm = gaussian_filter1d(filled, HOVM_SMOOTH_SIGMA)
            wt = gaussian_filter1d((~mask).astype(float), HOVM_SMOOTH_SIGMA)
            with np.errstate(invalid="ignore"):
                matrix[:, yi] = np.where(wt > 0.05, sm / wt, np.nan)

    # Colour scale
    vabs = np.nanpercentile(np.abs(matrix), 98)
    vabs = np.ceil(vabs / 5) * 5

    from matplotlib.colors import LinearSegmentedColormap
    hovm_cmap = LinearSegmentedColormap.from_list(
        "hovm_ea",
        [(0.00, "#4a0010"), (0.15, "#b2182b"), (0.35, "#ef8a62"),
         (0.50, "#f7f7f7"), (0.65, "#67a9cf"), (0.85, "#2166ac"),
         (1.00, "#08306b")]
    )

    fig = plt.figure(figsize=(16, 9))
    gs  = gridspec.GridSpec(2, 1, height_ratios=[5, 1], hspace=0.08)

    # ── Main heatmap ──────────────────────────────────────────────────────────
    ax_hov = fig.add_subplot(gs[0])
    im = ax_hov.imshow(
        matrix, aspect="auto", origin="lower",
        extent=[years[0] - 0.5, years[-1] + 0.5,
                domains[0]  - 0.5, domains[-1] + 0.5],
        cmap=hovm_cmap,
        vmin=-vabs, vmax=vabs,
        interpolation="nearest"
    )

    # Period divider
    ax_hov.axvline(2004, color="white", lw=1.5, ls="--", zorder=3)
    ax_hov.text(2004.3, domains[-1] - 1, "P1 | P2", color="white",
                fontsize=7, va="top", zorder=4)

    # Geographic annotations (horizontal = domain on y-axis here, so use axhspan)
    ax_hov.axhspan(ANN_WIMBLE_SHOALS[0]-0.5, ANN_WIMBLE_SHOALS[1]+0.5,
                   color=ANN_C_WIMBLE, alpha=0.15, zorder=1)
    for name, (d0, d1) in ANN_TOWN_SPANS.items():
        ax_hov.axhspan(d0-0.5, d1+0.5, color=ANN_C_TOWN_SPAN, alpha=0.15, zorder=1)
        ax_hov.text(years[0] + 0.5, (d0+d1)/2, name,
                    va="center", ha="left", fontsize=6,
                    color=ANN_C_TOWN_SPAN, fontweight="bold", zorder=4)
    for pname, pdom in ANN_PIERS.items():
        ax_hov.axhline(pdom, color=ANN_C_PIER, lw=0.7, ls="--", zorder=2)
    for gname, gdom in ANN_GROINS.items():
        ax_hov.axhline(gdom, color=ANN_C_GROIN, lw=0.7, ls="--", zorder=2)
    ax_hov.axhline(ANN_WIMBLE_SHOALS[0], color=ANN_C_WIMBLE, lw=0.5, ls=":", zorder=2)
    ax_hov.axhline(ANN_WIMBLE_SHOALS[1], color=ANN_C_WIMBLE, lw=0.5, ls=":", zorder=2)

    cb = fig.colorbar(im, ax=ax_hov, orientation="vertical",
                      pad=0.01, shrink=0.85)
    cb.set_label("Deviation from 1984\ndomain position (m)\nblue=accretion · red=erosion",
                 fontsize=7)
    cb.ax.tick_params(labelsize=7)

    ax_hov.set_ylabel("CASCADE domain  (1=Buxton  →  90=Rodanthe)", fontsize=8)
    ax_hov.set_yticks(np.arange(1, NUM_REAL_DOMAINS + 1, 5))
    ax_hov.set_yticklabels([str(d) for d in np.arange(1, NUM_REAL_DOMAINS + 1, 5)],
                            fontsize=6)
    ax_hov.tick_params(axis="x", labelsize=7)
    ax_hov.set_title(
        "Hatteras Island  |  Hovmöller: domain × year shoreline deviation\n"
        "with trajectory classification overlay",
        fontsize=11, fontweight="bold")

    # ── Classification colour strip below heatmap ─────────────────────────────
    ax_cl = fig.add_subplot(gs[1])
    for dom in domains:
        cl  = metrics.loc[dom, "trajectory"]
        col = CLASS_COLORS.get(cl, "#ffffff")
        ax_cl.bar(dom, 1, width=0.95, color=col,
                  edgecolor="white", linewidth=0.3, zorder=2)
        if metrics.loc[dom, "high_var"]:
            ax_cl.bar(dom, 0.25, width=0.95, bottom=0.75,
                      color="#ff7f00", zorder=3, alpha=0.7)

    ax_cl.set_xlim(domains[0] - 0.5, domains[-1] + 0.5)
    ax_cl.set_ylim(0, 1)
    ax_cl.set_yticks([])
    ax_cl.set_xlabel("CASCADE domain", fontsize=8)
    ax_cl.set_xticks(np.arange(1, NUM_REAL_DOMAINS + 1, 5))
    ax_cl.set_xticklabels([str(d) for d in np.arange(1, NUM_REAL_DOMAINS + 1, 5)],
                           fontsize=6)
    ax_cl.tick_params(labelsize=7)

    # Wimble + community on strip
    ax_cl.axvspan(ANN_WIMBLE_SHOALS[0]-0.5, ANN_WIMBLE_SHOALS[1]+0.5,
                  color=ANN_C_WIMBLE, alpha=0.15, zorder=0)
    for name, (d0, d1) in ANN_TOWN_SPANS.items():
        ax_cl.axvspan(d0-0.5, d1+0.5, color=ANN_C_TOWN_SPAN, alpha=0.18, zorder=0)

    patches = [mpatches.Patch(color=v, label=k)
               for k, v in CLASS_COLORS.items() if k != "Insufficient Data"]
    patches.append(mpatches.Patch(color="#ff7f00", alpha=0.7, label="High variability flag"))
    ax_cl.legend(handles=patches, fontsize=6, loc="upper right",
                 ncol=4, framealpha=0.9, title="Trajectory class", title_fontsize=6)

    fig.savefig(out_path, dpi=180, bbox_inches="tight")
    plt.close(fig)
    print(f"  Figure 3 saved → {out_path}")


# ============================================================
# MAIN
# ============================================================

def main():
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    # ── 1. Load lookup ────────────────────────────────────────────────────────
    print("Loading lookup table …")
    lookup = pd.read_csv(LOOKUP_CSV)
    lookup["transect_id"]   = lookup["transect_id"].astype(str)
    lookup["domain_number"] = lookup["domain_number"].astype(int)
    lookup = lookup[
        (lookup["domain_number"] >= 1) &
        (lookup["domain_number"] <= NUM_REAL_DOMAINS)
    ].sort_values(["domain_number", "transect_id"]).reset_index(drop=True)
    transect_order      = lookup["transect_id"].tolist()
    domain_per_transect = lookup["domain_number"].values
    n_transects         = len(transect_order)
    print(f"  {n_transects} transects across {lookup['domain_number'].nunique()} domains.\n")

    # ── 2. Load transects ─────────────────────────────────────────────────────
    print("Discovering CSV files …")
    csv_map = find_csvs(ROOT_DATA_DIR, SITE_FILTER)
    print(f"  Found {len(csv_map)} CSVs on disk.\n")

    print("Loading transect time-series …")
    ts_dict = {}
    for tid in tqdm(transect_order, unit="transect"):
        if tid in csv_map:
            ts = load_transect(csv_map[tid])
            if ts is not None and len(ts) > 0:
                ts_dict[tid] = ts
    print(f"  Loaded {len(ts_dict)} / {n_transects} transects.\n")

    # ── 3. Compute metrics ────────────────────────────────────────────────────
    print("Computing domain metrics …")
    metrics = compute_domain_metrics(ts_dict, transect_order, domain_per_transect)
    print()

    # ── 4. Save CSV ───────────────────────────────────────────────────────────
    csv_out = os.path.join(OUTPUT_DIR, "domain_trajectory_metrics.csv")
    metrics.drop(columns=["_ann_full"]).to_csv(csv_out)
    print(f"  Metrics CSV saved → {csv_out}")

    # Print summary
    print("\n  Trajectory class counts:")
    counts = metrics["trajectory"].value_counts()
    for cls, n in counts.items():
        print(f"    {cls:30s}: {n} domains")

    # ── 5. Figures ────────────────────────────────────────────────────────────
    print("\n" + "=" * 60)
    print("Generating figures …")
    print("=" * 60)

    plot_classification_bar(
        metrics,
        os.path.join(OUTPUT_DIR, "fig1_classification_bar.png"))

    plot_lrr_scatter(
        metrics,
        os.path.join(OUTPUT_DIR, "fig2_lrr_scatter.png"))

    plot_hovmoller(
        metrics,
        os.path.join(OUTPUT_DIR, "fig3_hovmoller_classified.png"))

    print("\nDone.")


if __name__ == "__main__":
    main()
