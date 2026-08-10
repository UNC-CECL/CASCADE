"""
HAT_be_zone_analysis.py
========================
Derives defensible background erosion (BE) source/sink corrections for CASCADE
from the residual between a LOESS-smoothed CoastSat observed shoreline change
rate and the CASCADE base-run LRR.

Philosophy
----------
Corrections are applied only where THREE conditions are simultaneously met:
  1. The residual (smoothed observed rate minus CASCADE LRR) exceeds a
     significance threshold (not just noise)
  2. The signal is spatially coherent across multiple adjacent domains
  3. A physical mechanism can be named for that zone

This minimises the number of free parameters and maximises scientific
defensibility — every correction has a name and a reason.

Workflow
--------
  1. Load CoastSat domain-averaged LRR (P1 and P2) — the observed shoreline
     change rate
  2. Load CASCADE base-run LRR from NPZ (P1 and P2, management included, BE=0)
  3. LOESS-smooth the observed shoreline change rate (10-domain window),
     excluding domains 1-GROIN_EXCLUDE_THROUGH_DOMAIN (Buxton groin influence
     zone) from the fit entirely — those domains pass through with their raw,
     unsmoothed rate
  4. Compute residual = smoothed CoastSat rate - CASCADE LRR, per domain per
     period (the fully raw, unsmoothed residual is also retained, for
     diagnostic comparison only — it no longer drives any decision)
  5. Identify significant zones (|residual| > SIGNIFICANCE_THRESHOLD,
     spatially coherent over >= MIN_ZONE_WIDTH adjacent domains)
  6. Classify each domain: stable correction vs shifting between periods
  7. Output:
       - Diagnostic figures showing raw/smoothed rate, residual, and zone
         identification
       - DOMAIN_BE_RATES dicts for P1, P2, and three forecast scenarios
       - Summary CSV with all metrics and physical zone assignments

Forecast scenarios (for domains where P1 ≠ P2 correction)
----------------------------------------------------------
  "continue" : use P2 correction (current trajectory continues into future)
  "revert"   : use P1 correction (system returns to pre-2004 state)
  "neutral"  : use mean(P1, P2) (no prior on future state)

Usage
-----
  python HAT_be_zone_analysis.py

Dependencies
------------
  pip install pandas numpy matplotlib scipy statsmodels tqdm
"""

import os
import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.nonparametric.smoothers_lowess import lowess
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.gridspec as gridspec
from tqdm import tqdm

# ============================================================
# CONFIG — edit paths and thresholds here
# ============================================================

# CoastSat domain-averaged LRR CSVs
P1_COASTSAT_CSV = r"C:\Users\hanna\PycharmProjects\CASCADE\scripts\input_preperation\CoastSat\1984_2004\domain_lrr_1984_2004_summary.csv"
P2_COASTSAT_CSV = r"C:\Users\hanna\PycharmProjects\CASCADE\scripts\input_preperation\CoastSat\2004_2024\domain_lrr_2004_2024_summary.csv"

# CASCADE base run NPZs (management included, BE=0 throughout)
P1_CASCADE_NPZ  = r"C:\Users\hanna\PycharmProjects\CASCADE\output\raw_runs\HAT_1984_2004_edge_calibrated\HAT_1984_2004_edge_calibrated.npz"
P2_CASCADE_NPZ  = r"C:\Users\hanna\PycharmProjects\CASCADE\output\raw_runs\HAT_2004_2024_edge_calibrated\HAT_2004_2024_edge_calibrated.npz"

OUTPUT_DIR = r"/scripts/input_prep/loess_smooth\output"

# ── Column names in CoastSat CSVs ─────────────────────────────────────────────
LRR_COL    = "median_lrr"   # use median — more robust to outlier transects
DOMAIN_COL = "domain_number"

# ── CASCADE structure ─────────────────────────────────────────────────────────
NUM_REAL_DOMAINS = 90
START_REAL_INDEX = 15        # buffer domains before domain 1
CASCADE_SIGN     = -1        # x_s_TS increases landward = erosion → flip to standard
P1_START, P1_END = 1984, 2004
P2_START, P2_END = 2004, 2024

# ── Correction thresholds ─────────────────────────────────────────────────────
# Minimum smoothed residual magnitude to warrant any correction at all.
# Below this = within model noise/uncertainty, leave at zero.
SIGNIFICANCE_THRESHOLD = 0.5   # m/yr

# Minimum number of adjacent domains with |residual| > threshold to form a zone.
# Prevents correcting isolated noisy domains.
MIN_ZONE_WIDTH = 3   # domains

# If |P1_correction - P2_correction| exceeds this, the zone is "shifting"
# and needs period-specific values + forecast scenarios.
SHIFT_THRESHOLD = 0.75   # m/yr

# ── LOESS smoothing ───────────────────────────────────────────────────────────
# Fraction of data used for each local regression (larger = smoother).
# 10-domain window matches the CoastSat LOESS calibration window used
# throughout the dissertation. As of this version, smoothing is applied to
# the OBSERVED SHORELINE CHANGE RATE itself (before differencing against
# CASCADE) — not to the residual. LOESS_FRAC is calibrated for the full
# 90-domain array; see smooth_shoreline_rate() for how it's re-derived once
# the groin zone is excluded from the fit.
LOESS_FRAC = 0.111  # exactly 10 domains at 90 total
LOESS_WINDOW_DOMAINS = 10  # the actual window width LOESS_FRAC is calibrated to hit

# Domains 1 through this value are excluded ENTIRELY from the shoreline-rate
# LOESS fit (Buxton groin influence zone) — the groin's localized signal
# would otherwise bleed into the smoothed estimate at neighbouring domains.
# These domains always keep their raw, unsmoothed CoastSat rate.
GROIN_EXCLUDE_THROUGH_DOMAIN = 10

# ── Manual overrides ─────────────────────────────────────────────────────────
# Domain-level corrections that override the LOESS-derived value.
# Use sparingly — only where the smoothing window demonstrably under/over-corrects
# and you have a clear physical justification for the different value.
# Format: domain → (p1_override, p2_override, reason)
# Set either value to None to keep the LOESS-derived value for that period.
# Forecast scenarios for overridden domains use the same logic as normal:
#   continue = p2, revert = p1, neutral = mean(p1, p2)
MANUAL_OVERRIDES = {
    # No overrides — pure LOESS-derived values against the current baseline.
    # This script is for DISCOVERY: see what the data says before any
    # manual intervention. Once you have chosen final values (informed by
    # this comparison plus your own judgement), enter them in
    # HAT_be_apply_final_rates.py to generate the final figures and
    # forecast scenarios.
}

# ── Locked domains ────────────────────────────────────────────────────────────
# Domains where the BE rate has already been solved independently (e.g. D1 and
# D90, found by reproducing the buffer-cell shoreline change rate directly).
# These are forced to 0.0 in ALL DOMAIN_BE_RATES outputs and excluded from the
# significance/strategy calculation entirely — the script will not suggest a
# correction for these domains, since you will always supply your own value.
# Add a short note here for your own records of what each locked value is.
LOCKED_DOMAINS = {
    1:  "Solved directly via buffer-cell reproduction (see GIS 1 value)",
    90: "Solved directly via buffer-cell reproduction (see GIS 90 value)",
}

# ── Physical zone definitions ─────────────────────────────────────────────────
# These are your prior hypotheses about where physical mechanisms operate.
# The script will test whether the residual data supports them.
# Format: zone_name → (domain_start, domain_end, mechanism_description)
PHYSICAL_ZONES = {
    "Cape Point / Shoal Dynamics":  (1,  10,  "Cape/shoal attachment-detachment cycle + post-Isabel recovery"),
    "Buxton–Avon Transition":       (9,  20,  "Post-Isabel geomorphic recovery, background SLR erosion"),
    "Avon":                         (21, 31,  "Pier-induced sediment shadow, nourishment interactions"),
    "Mid-island":                   (32, 59,  "Background SLR-driven erosion, minimal local forcing"),
    "Wimble Shoals Influence":      (60, 74,  "Offshore shoal migration — sediment delivery to nearshore"),
    "Tri-Village / Rodanthe":       (75, 83,  "Persistent erosion hotspot, event-driven, infrastructure effects"),
    "Pea Island NWR":               (84, 90,  "Oregon Inlet dynamics, northern Wimble Shoals influence"),
}

# Display-only shortenings for the in-place zone strip labels (fig_be_rates).
# The canonical PHYSICAL_ZONES names above are kept everywhere else — CSV
# comparison, DOMAIN_BE_RATES comments, mechanism lookups — since the longer
# names carry more information there. Add more entries here if you want
# other zones shortened on the chart too.
ZONE_DISPLAY_NAMES = {
    "Cape Point / Shoal Dynamics":  "Cape Hatteras",
    "Buxton–Avon Transition":       "Buxton-Avon",
    "Wimble Shoals Influence":      "Wimble Shoals",
    "Tri-Village / Rodanthe":       "Tri-Village",
}

# ── Annotation colours (publication style) ───────────────────────────────────
ANN_TOWN_SPANS = {
    "Buxton":      (7,  8),
    "Avon":        (21, 31),
    "Tri-Village": (68, 83),
}
ANN_WIMBLE_SHOALS = (60, 74)
ANN_PIERS   = {"Avon Pier": 26, "Rodanthe Pier": 79}
ANN_GROINS  = {"Buxton Groin": 5.5}
ANN_C_TOWN  = "#90AFC5"
ANN_C_WIMBLE = "#E0A800"
ANN_C_PIER  = "#1565C0"
ANN_C_GROIN = "#B71C1C"

# ── Figure font sizes ─────────────────────────────────────────────────────────
# Centralised here so the whole figure can be resized for readability (e.g.
# for sharing with Laura / printing) by editing these six values instead of
# hunting through every plot function.
FONT_SUPTITLE = 18   # main figure title (fig.suptitle)
FONT_TITLE    = 13   # per-panel titles
FONT_LABEL    = 12   # axis labels (x/y)
FONT_TICK     = 11   # tick labels
FONT_LEGEND   = 11   # legend text
FONT_ANNOT    = 10   # small annotation text (town/pier/groin labels in annotate_ax)

# ============================================================
# CASCADE LOADER
# ============================================================

def load_cascade_lrr(npz_path, start_year, end_year):
    """Extract per-domain LRR (m/yr) from CASCADE base-run NPZ."""
    print(f"  Loading: {os.path.basename(npz_path)}")
    data    = np.load(npz_path, allow_pickle=True)
    cascade = data["cascade"][0]
    b3d     = cascade.barrier3d
    n_years = end_year - start_year
    years   = np.arange(start_year, end_year + 1)

    lrr = {}
    for dom in range(1, NUM_REAL_DOMAINS + 1):
        idx   = START_REAL_INDEX + (dom - 1)
        b3d_i = b3d[idx]
        if hasattr(b3d_i, "x_s_TS"):
            xs = np.array(b3d_i.x_s_TS, dtype=float)
        elif hasattr(b3d_i, "_x_s_TS"):
            xs = np.array(b3d_i._x_s_TS, dtype=float)
        else:
            lrr[dom] = np.nan; continue

        xs_m = xs * 10.0 * CASCADE_SIGN   # dam → m, flip sign
        nt   = min(len(xs_m), n_years + 1)
        if nt < 4:
            lrr[dom] = np.nan; continue

        slope, *_ = stats.linregress(years[:nt], xs_m[:nt])
        lrr[dom]  = slope

    return pd.Series(lrr, name="cascade_lrr")


# ============================================================
# SMOOTHING
# ============================================================

def loess_smooth(values, frac=LOESS_FRAC):
    """Apply LOESS smoothing to a domain-indexed array, handling NaNs."""
    domains = np.arange(1, len(values) + 1, dtype=float)
    mask    = ~np.isnan(values)
    if mask.sum() < 5:
        return values.copy()
    smoothed_valid = lowess(values[mask], domains[mask],
                            frac=frac, return_sorted=False)
    out = np.full_like(values, np.nan)
    out[mask] = smoothed_valid
    return out


def smooth_shoreline_rate(raw_rate, exclude_through=GROIN_EXCLUDE_THROUGH_DOMAIN,
                          window_domains=LOESS_WINDOW_DOMAINS):
    """
    LOESS-smooth a domain-indexed OBSERVED shoreline change rate, excluding
    domains 1..exclude_through entirely from the fit (Buxton groin influence
    zone) and passing those domains through unchanged with their raw rate.

    Domains > exclude_through are smoothed using ONLY data from domains
    > exclude_through — the groin-zone values never enter the regression at
    all, so they cannot bleed into the smoothed estimate near the zone
    boundary (this is stronger than merely overwriting the comparison for
    domains 1..exclude_through after smoothing over the full array).

    frac is re-derived here rather than reusing LOESS_FRAC directly: LOESS_FRAC
    (0.111) is calibrated to give a 10-domain window when fit over all 90
    domains. Once the groin zone is excluded, only 80 domains remain in the
    fit, so reusing 0.111 unchanged would narrow the window to ~8.9 domains.
    Recomputing frac = window_domains / n_valid preserves the true ~10-domain
    (~5 km) window this dissertation uses everywhere else.
    """
    n_valid = len(raw_rate) - exclude_through
    frac    = window_domains / n_valid

    fit_input = raw_rate.copy()
    fit_input[:exclude_through] = np.nan                       # groin zone never enters the fit
    smoothed  = loess_smooth(fit_input, frac=frac)              # NaN-aware LOESS
    smoothed[:exclude_through] = raw_rate[:exclude_through]     # restore raw, unsmoothed values
    return smoothed


# ============================================================
# ZONE IDENTIFICATION
# ============================================================

def identify_correction_zones(smoothed_residual, min_width=MIN_ZONE_WIDTH,
                               threshold=SIGNIFICANCE_THRESHOLD):
    """
    Find contiguous runs of domains where |smoothed_residual| > threshold
    and the run is at least min_width domains wide.

    Returns array of booleans: True = correction warranted.
    """
    domains  = np.arange(1, NUM_REAL_DOMAINS + 1)
    sig      = np.abs(smoothed_residual) > threshold
    warranted = np.zeros(NUM_REAL_DOMAINS, dtype=bool)

    # Find contiguous runs
    i = 0
    while i < NUM_REAL_DOMAINS:
        if sig[i]:
            j = i
            while j < NUM_REAL_DOMAINS and sig[j]:
                j += 1
            run_length = j - i
            if run_length >= min_width:
                warranted[i:j] = True
            i = j
        else:
            i += 1
    return warranted


def assign_physical_zone(domain):
    """Return the physical zone name for a given domain number."""
    for zone_name, (d0, d1, _) in PHYSICAL_ZONES.items():
        if d0 <= domain <= d1:
            return zone_name
    return "Unassigned"


# ============================================================
# BE RATE COMPUTATION
# ============================================================

def compute_be_rates(raw_p1, raw_p2, smooth_p1, smooth_p2):
    """
    For each domain, determine the appropriate BE correction strategy.

    Rules:
      - If smoothed residual not significant in either period → BE = 0
      - If significant in one or both periods:
          - If |correction_P1 - correction_P2| < SHIFT_THRESHOLD → stable → use mean
          - Otherwise → shifting → flag for scenario treatment, provide P1/P2/mean
    """
    domains = np.arange(1, NUM_REAL_DOMAINS + 1)
    rows    = []

    sig_p1 = identify_correction_zones(smooth_p1)
    sig_p2 = identify_correction_zones(smooth_p2)

    for i, dom in enumerate(domains):
        r1_raw  = raw_p1[i]
        r2_raw  = raw_p2[i]
        r1_sm   = smooth_p1[i]
        r2_sm   = smooth_p2[i]
        s1      = sig_p1[i]
        s2      = sig_p2[i]

        # Locked domains: force to 0.0, skip all significance/strategy logic.
        # These domains already have an independently-solved BE rate that you
        # will always supply yourself — the script should not suggest a value.
        if dom in LOCKED_DOMAINS:
            strategy = "locked"
            be_hindcast_p1   = 0.0
            be_hindcast_p2   = 0.0
            be_forecast_continue = 0.0
            be_forecast_revert   = 0.0
            be_forecast_neutral  = 0.0
            zone = assign_physical_zone(dom)
            mechanism = f"LOCKED — {LOCKED_DOMAINS[dom]}"
            rows.append({
                "domain":               dom,
                "raw_residual_p1":      r1_raw,
                "raw_residual_p2":      r2_raw,
                "smooth_residual_p1":   r1_sm,
                "smooth_residual_p2":   r2_sm,
                "correction_warranted_p1": False,
                "correction_warranted_p2": False,
                "strategy":             strategy,
                "be_hindcast_p1":       be_hindcast_p1,
                "be_hindcast_p2":       be_hindcast_p2,
                "be_forecast_continue": be_forecast_continue,
                "be_forecast_revert":   be_forecast_revert,
                "be_forecast_neutral":  be_forecast_neutral,
                "physical_zone":        zone,
                "mechanism":            mechanism,
            })
            continue

        # Correction value to apply: use smoothed residual (spatially coherent)
        # Only apply where warranted
        corr_p1 = r1_sm if s1 else 0.0
        corr_p2 = r2_sm if s2 else 0.0

        # Is either period significant?
        any_sig = s1 or s2

        if not any_sig:
            strategy = "zero"
            be_hindcast_p1   = 0.0
            be_hindcast_p2   = 0.0
            be_forecast_continue = 0.0
            be_forecast_revert   = 0.0
            be_forecast_neutral  = 0.0
        else:
            delta = abs(corr_p1 - corr_p2)
            if delta < SHIFT_THRESHOLD:
                strategy = "stable"
                be_mean = float(np.nanmean([corr_p1, corr_p2]))
                be_hindcast_p1   = be_mean
                be_hindcast_p2   = be_mean
                be_forecast_continue = be_mean
                be_forecast_revert   = be_mean
                be_forecast_neutral  = be_mean
            else:
                strategy = "shifting"
                be_hindcast_p1   = corr_p1
                be_hindcast_p2   = corr_p2
                # Forecast scenarios
                be_forecast_continue = corr_p2          # current state continues
                be_forecast_revert   = corr_p1          # reverts to pre-2004 state
                be_forecast_neutral  = float(np.nanmean([corr_p1, corr_p2]))

        zone = assign_physical_zone(dom)
        _, _, mechanism = PHYSICAL_ZONES.get(zone, (None, None, "Unknown"))

        # Apply manual overrides if defined for this domain
        if dom in MANUAL_OVERRIDES:
            ov_p1, ov_p2, ov_reason = MANUAL_OVERRIDES[dom]
            if ov_p1 is not None:
                be_hindcast_p1       = ov_p1
                be_forecast_revert   = ov_p1
            if ov_p2 is not None:
                be_hindcast_p2       = ov_p2
                be_forecast_continue = ov_p2
            # Recalculate neutral as mean of (possibly overridden) P1 and P2
            be_forecast_neutral = float(np.nanmean([be_hindcast_p1, be_hindcast_p2]))
            strategy = strategy + "*"  # flag as manually overridden in comparison
            mechanism = ov_reason

        rows.append({
            "domain":               dom,
            "raw_residual_p1":      r1_raw,
            "raw_residual_p2":      r2_raw,
            "smooth_residual_p1":   r1_sm,
            "smooth_residual_p2":   r2_sm,
            "correction_warranted_p1": s1,
            "correction_warranted_p2": s2,
            "strategy":             strategy,
            "be_hindcast_p1":       be_hindcast_p1,
            "be_hindcast_p2":       be_hindcast_p2,
            "be_forecast_continue": be_forecast_continue,
            "be_forecast_revert":   be_forecast_revert,
            "be_forecast_neutral":  be_forecast_neutral,
            "physical_zone":        zone,
            "mechanism":            mechanism,
        })

    return pd.DataFrame(rows).set_index("domain")


# ============================================================
# ANNOTATION HELPER
# ============================================================

def annotate_ax(ax, ylim):
    ymin, ymax = ylim; yspan = ymax - ymin
    ax.axvspan(ANN_WIMBLE_SHOALS[0]-0.5, ANN_WIMBLE_SHOALS[1]+0.5,
               color=ANN_C_WIMBLE, alpha=0.12, zorder=0)
    for name, (d0, d1) in ANN_TOWN_SPANS.items():
        ax.axvspan(d0-0.5, d1+0.5, color=ANN_C_TOWN, alpha=0.18, zorder=0)
        ax.text((d0+d1)/2, ymax - 0.02*yspan, name, ha="center", va="top",
                fontsize=FONT_ANNOT, color=ANN_C_TOWN, fontweight="bold")
    for pname, pdom in ANN_PIERS.items():
        ax.axvline(pdom, color=ANN_C_PIER, lw=0.9, ls="--", zorder=2)
        ax.text(pdom+0.3, ymin+0.72*yspan, pname, rotation=90,
                va="bottom", fontsize=FONT_ANNOT, color=ANN_C_PIER)
    for gname, gdom in ANN_GROINS.items():
        ax.axvline(gdom, color=ANN_C_GROIN, lw=0.9, ls="--", zorder=2)
        ax.text(gdom+0.3, ymin+0.62*yspan, gname, rotation=90,
                va="bottom", fontsize=FONT_ANNOT, color=ANN_C_GROIN)
    ax.axhline(0, color="black", lw=0.8, ls="--", alpha=0.5)
    ax.axhline( SIGNIFICANCE_THRESHOLD, color="#888", lw=0.6, ls=":", alpha=0.7)
    ax.axhline(-SIGNIFICANCE_THRESHOLD, color="#888", lw=0.6, ls=":", alpha=0.7)
    ax.set_xlim(0.5, NUM_REAL_DOMAINS + 0.5)


def find_zone_runs(zone_series, domains):
    """Collapse a per-domain physical_zone Series into contiguous
    (start_domain, end_domain, zone_name) runs, in domain order."""
    runs = []
    current_zone = None
    run_start = None
    for dom in domains:
        z = zone_series.loc[dom]
        if z != current_zone:
            if current_zone is not None:
                runs.append((run_start, prev_dom, current_zone))
            current_zone, run_start = z, dom
        prev_dom = dom
    runs.append((run_start, prev_dom, current_zone))
    return runs


def label_zone_runs(ax, fig, runs, y=0.5, fontsize_start=FONT_LABEL,
                    fontsize_min=6.5, pad_frac=0.90):
    """
    Place one centered text label per zone run directly on the strip,
    shrinking that label's own font (and only that one) until it fits
    within its own zone's width — so a narrow zone with a long name
    (e.g. "Cape Point / Shoal Dynamics") never bleeds into its neighbour.
    """
    renderer = fig.canvas.get_renderer()
    for d0, d1, name in runs:
        if name is None or name == "Unassigned":
            continue
        display_name = ZONE_DISPLAY_NAMES.get(name, name)
        xc = (d0 + d1) / 2.0
        x0_disp = ax.transData.transform((d0 - 0.5, y))[0]
        x1_disp = ax.transData.transform((d1 + 0.5, y))[0]
        avail_px = abs(x1_disp - x0_disp) * pad_frac

        fs = fontsize_start
        txt = ax.text(xc, y, display_name, ha="center", va="center",
                      fontsize=fs, color="black", zorder=5, clip_on=False,
                      bbox=dict(facecolor="white", alpha=0.78,
                                edgecolor="none", boxstyle="round,pad=0.15"))
        fig.canvas.draw()
        bbox = txt.get_window_extent(renderer=renderer)
        while bbox.width > avail_px and fs > fontsize_min:
            fs -= 0.5
            txt.set_fontsize(fs)
            fig.canvas.draw()
            bbox = txt.get_window_extent(renderer=renderer)


# ============================================================
# FIGURE 1 — Diagnostic: raw residual, smoothed, zone identification
# ============================================================

def plot_diagnostic(cs_p1, cs_p2, casc_p1, casc_p2,
                    cs_p1_smooth, cs_p2_smooth,
                    raw_p1, raw_p2, smooth_p1, smooth_p2,
                    results, out_path):

    domains = np.arange(1, NUM_REAL_DOMAINS + 1)
    fig, axes = plt.subplots(5, 1, figsize=(18, 18), sharex=True,
                             gridspec_kw={"height_ratios": [3, 3, 3, 3, 1.5]})

    # ── Panel 1: CoastSat vs CASCADE LRR — Period 1 ───────────────────────────
    ax = axes[0]
    ax.plot(domains, cs_p1,   "o-", ms=3, lw=1.0, color="#b2182b",
            label="CoastSat LRR (raw)", zorder=3, alpha=0.55)
    ax.plot(domains, cs_p1_smooth, "-", lw=1.8, color="#762a83",
            label=f"CoastSat LOESS-smoothed (D{GROIN_EXCLUDE_THROUGH_DOMAIN+1}+)",
            zorder=4)
    ax.plot(domains, casc_p1, "s-", ms=3, lw=1.0, color="#2166ac",
            label="CASCADE base LRR", zorder=3, alpha=0.8)
    ax.axvline(GROIN_EXCLUDE_THROUGH_DOMAIN + 0.5, color="#762a83",
               lw=0.8, ls=":", alpha=0.7, zorder=2)
    ax.set_ylabel("LRR  (m/yr)", fontsize=FONT_LABEL)
    ax.set_title(f"Period 1  (1984–2004)  |  CoastSat vs CASCADE base run  |  "
                 f"rate smoothing excludes D1-{GROIN_EXCLUDE_THROUGH_DOMAIN} (groin)",
                 fontsize=FONT_TITLE, loc="left", pad=3)
    ylim = (min(np.nanmin([cs_p1, cs_p1_smooth, casc_p1])*1.2, -3),
            max(np.nanmax([cs_p1, cs_p1_smooth, casc_p1])*1.2,  3))
    ax.set_ylim(ylim); ax.legend(fontsize=FONT_LEGEND, loc="upper center")
    ax.tick_params(labelsize=FONT_TICK)
    annotate_ax(ax, ylim)

    # ── Panel 2: CoastSat vs CASCADE LRR — Period 2 ───────────────────────────
    ax = axes[1]
    ax.plot(domains, cs_p2,   "o-", ms=3, lw=1.0, color="#b2182b",
            label="CoastSat LRR (raw)", zorder=3, alpha=0.55)
    ax.plot(domains, cs_p2_smooth, "-", lw=1.8, color="#762a83",
            label=f"CoastSat LOESS-smoothed (D{GROIN_EXCLUDE_THROUGH_DOMAIN+1}+)",
            zorder=4)
    ax.plot(domains, casc_p2, "s-", ms=3, lw=1.0, color="#2166ac",
            label="CASCADE base LRR", zorder=3, alpha=0.8)
    ax.axvline(GROIN_EXCLUDE_THROUGH_DOMAIN + 0.5, color="#762a83",
               lw=0.8, ls=":", alpha=0.7, zorder=2)
    ax.set_ylabel("LRR  (m/yr)", fontsize=FONT_LABEL)
    ax.set_title(f"Period 2  (2004–2024)  |  CoastSat vs CASCADE base run  |  "
                 f"rate smoothing excludes D1-{GROIN_EXCLUDE_THROUGH_DOMAIN} (groin)",
                 fontsize=FONT_TITLE, loc="left", pad=3)
    ylim = (min(np.nanmin([cs_p2, cs_p2_smooth, casc_p2])*1.2, -3),
            max(np.nanmax([cs_p2, cs_p2_smooth, casc_p2])*1.2,  3))
    ax.set_ylim(ylim); ax.legend(fontsize=FONT_LEGEND, loc="upper center")
    ax.tick_params(labelsize=FONT_TICK)
    annotate_ax(ax, ylim)

    # ── Panel 3: Raw vs smoothed-rate residual — Period 1 ────────────────────
    ax = axes[2]
    ax.bar(domains, raw_p1, width=0.7, color="#d9d9d9", alpha=0.6,
           label="Raw residual (unsmoothed both sides)", zorder=2)
    ax.plot(domains, smooth_p1, "-", lw=2.0, color="#b2182b",
            label="Residual (from smoothed rate)", zorder=3)
    # Shade warranted correction zones
    sig = results["correction_warranted_p1"].values
    for i, dom in enumerate(domains):
        if sig[i]:
            ax.axvspan(dom-0.5, dom+0.5, color="#b2182b", alpha=0.12, zorder=0)
    ax.set_ylabel("Residual  (m/yr)\nCoastSat − CASCADE", fontsize=FONT_LABEL)
    ax.set_title(f"Period 1 residual  |  shaded = correction warranted  "
                 f"(|residual| > {SIGNIFICANCE_THRESHOLD} m/yr, ≥{MIN_ZONE_WIDTH} domains wide)",
                 fontsize=FONT_TITLE, loc="left", pad=3)
    all_vals = np.concatenate([raw_p1[~np.isnan(raw_p1)],
                               smooth_p1[~np.isnan(smooth_p1)]])
    ylim = (min(np.nanmin(all_vals)*1.2, -3), max(np.nanmax(all_vals)*1.2, 3))
    ax.set_ylim(ylim); ax.legend(fontsize=FONT_LEGEND, loc="upper center")
    ax.tick_params(labelsize=FONT_TICK)
    annotate_ax(ax, ylim)

    # ── Panel 4: Raw vs smoothed-rate residual — Period 2 ────────────────────
    ax = axes[3]
    ax.bar(domains, raw_p2, width=0.7, color="#d9d9d9", alpha=0.6,
           label="Raw residual (unsmoothed both sides)", zorder=2)
    ax.plot(domains, smooth_p2, "-", lw=2.0, color="#2166ac",
            label="Residual (from smoothed rate)", zorder=3)
    sig = results["correction_warranted_p2"].values
    for i, dom in enumerate(domains):
        if sig[i]:
            ax.axvspan(dom-0.5, dom+0.5, color="#2166ac", alpha=0.12, zorder=0)
    ax.set_ylabel("Residual  (m/yr)\nCoastSat − CASCADE", fontsize=FONT_LABEL)
    ax.set_title(f"Period 2 residual  |  shaded = correction warranted",
                 fontsize=FONT_TITLE, loc="left", pad=3)
    all_vals = np.concatenate([raw_p2[~np.isnan(raw_p2)],
                               smooth_p2[~np.isnan(smooth_p2)]])
    ylim = (min(np.nanmin(all_vals)*1.2, -3), max(np.nanmax(all_vals)*1.2, 3))
    ax.set_ylim(ylim); ax.legend(fontsize=FONT_LEGEND, loc="upper center")
    ax.tick_params(labelsize=FONT_TICK)
    annotate_ax(ax, ylim)

    # ── Panel 5: Strategy strip ───────────────────────────────────────────────
    ax = axes[4]
    # Use .get() with strip of * so any starred variant falls back gracefully
    strat_colors = {"zero": "#ffffff", "stable": "#4dac26",
                    "shifting": "#d73027", "locked": "#999999"}
    for i, dom in enumerate(domains):
        strat = results.loc[dom, "strategy"]
        ax.bar(dom, 1, width=0.9, color=strat_colors.get(strat.rstrip("*"), "#cccccc"),
               edgecolor="#cccccc", linewidth=0.3, zorder=2)
    ax.set_yticks([])
    ax.set_ylabel("BE strategy", fontsize=FONT_LABEL)
    ax.set_xlabel("CASCADE domain  (1 = Buxton  →  90 = Rodanthe)", fontsize=FONT_LABEL)
    ax.set_ylim(0, 1); ax.tick_params(labelsize=FONT_TICK)
    ax.axvspan(ANN_WIMBLE_SHOALS[0]-0.5, ANN_WIMBLE_SHOALS[1]+0.5,
               color=ANN_C_WIMBLE, alpha=0.12, zorder=0)
    for name, (d0, d1) in ANN_TOWN_SPANS.items():
        ax.axvspan(d0-0.5, d1+0.5, color=ANN_C_TOWN, alpha=0.18, zorder=0)

    patches = [
        mpatches.Patch(color="#ffffff", ec="#888", label="No correction (within noise)"),
        mpatches.Patch(color="#4dac26", label="Stable correction (single value P1=P2)"),
        mpatches.Patch(color="#d73027", label="Shifting (period-specific + forecast scenarios)"),
    ]
    ax.legend(handles=patches, fontsize=FONT_LEGEND, loc="upper right", framealpha=0.9)

    fig.suptitle(
        "Hatteras Island  |  BE source/sink zone identification\n"
        f"Significance threshold: ±{SIGNIFICANCE_THRESHOLD} m/yr  ·  "
        f"Min zone width: {MIN_ZONE_WIDTH} domains  ·  "
        f"Shift threshold: {SHIFT_THRESHOLD} m/yr",
        fontsize=FONT_SUPTITLE, fontweight="bold")
    fig.tight_layout(rect=[0, 0, 1, 0.97])
    fig.savefig(out_path, dpi=180, bbox_inches="tight")
    plt.close(fig)
    print(f"  Diagnostic figure saved → {out_path}")


# ============================================================
# FIGURE 2 — Final BE rates: hindcast + forecast scenarios
# ============================================================

def plot_be_rates(results, out_path):
    domains = np.arange(1, NUM_REAL_DOMAINS + 1)
    fig, axes = plt.subplots(3, 1, figsize=(18, 12), sharex=True,
                             gridspec_kw={"height_ratios": [4, 4, 1.5]})

    # ── Panel 1: Hindcast BE rates ────────────────────────────────────────────
    ax = axes[0]
    be_p1 = results["be_hindcast_p1"].values
    be_p2 = results["be_hindcast_p2"].values
    w = 0.38
    ax.bar(domains - w/2, be_p1, width=w, color="#4575b4", alpha=0.85,
           label="P1 hindcast BE  (1984–2004)", zorder=2)
    ax.bar(domains + w/2, be_p2, width=w, color="#d73027", alpha=0.85,
           label="P2 hindcast BE  (2004–2024)", zorder=2)
    ylim = (min(np.nanmin([be_p1, be_p2])*1.3, -2),
            max(np.nanmax([be_p1, be_p2])*1.3,  2))
    ax.set_ylim(ylim)
    ax.legend(fontsize=FONT_LEGEND, loc="lower center", ncol=2,
              frameon=True, framealpha=0.9)
    ax.set_ylabel("BE rate  (m/yr)\n+ = accretion source\n− = erosion sink", fontsize=FONT_LABEL)
    ax.set_title("Hindcast BE rates  |  only applied where physically warranted",
                 fontsize=FONT_TITLE, loc="left", pad=3)
    ax.tick_params(labelsize=FONT_TICK)
    annotate_ax(ax, ylim)

    # Mark shifting zones
    for i, dom in enumerate(domains):
        if results.loc[dom, "strategy"] == "shifting":
            ax.axvspan(dom-0.5, dom+0.5, color="#ff7f00", alpha=0.10, zorder=0)

    # ── Panel 2: Forecast scenario BE rates ───────────────────────────────────
    ax = axes[1]
    be_cont = results["be_forecast_continue"].values
    be_rev  = results["be_forecast_revert"].values
    be_neut = results["be_forecast_neutral"].values

    ax.fill_between(domains,
                    np.minimum(be_rev, be_cont),
                    np.maximum(be_rev, be_cont),
                    color="#aaaaaa", alpha=0.30, label="Scenario range (revert↔continue)")
    ax.plot(domains, be_cont, "o-", ms=3, lw=1.2, color="#d73027",
            label='Scenario "continue" (P2 state)', zorder=3)
    ax.plot(domains, be_rev,  "s-", ms=3, lw=1.2, color="#4575b4",
            label='Scenario "revert" (P1 state)', zorder=3, alpha=0.8)
    ax.plot(domains, be_neut, "^-", ms=3, lw=1.0, color="#333333",
            ls="--", label='Scenario "neutral" (mean)', zorder=3, alpha=0.7)

    ylim = (min(np.nanmin([be_cont, be_rev])*1.3, -2),
            max(np.nanmax([be_cont, be_rev])*1.3,  2))
    ax.set_ylim(ylim)
    ax.legend(fontsize=FONT_LEGEND, loc="lower center", ncol=4,
              frameon=True, framealpha=0.9)
    ax.set_ylabel("BE rate  (m/yr)", fontsize=FONT_LABEL)
    ax.set_title("Forecast scenario BE rates  |  grey band = physical uncertainty range",
                 fontsize=FONT_TITLE, loc="left", pad=3)
    ax.tick_params(labelsize=FONT_TICK)
    annotate_ax(ax, ylim)

    # ── Panel 3: Zone identification strip ────────────────────────────────────
    ax = axes[2]
    zone_palette = {
        "Cape Point / Shoal Dynamics":  "#b2182b",
        "Buxton–Avon Transition":       "#f4a582",
        "Avon":                         "#4393c3",
        "Mid-island":                   "#d9d9d9",
        "Wimble Shoals Influence":      "#E0A800",
        "Tri-Village / Rodanthe":       "#762a83",
        "Pea Island NWR":               "#4dac26",
    }
    for i, dom in enumerate(domains):
        zone  = results.loc[dom, "physical_zone"]
        strat = results.loc[dom, "strategy"]
        color = zone_palette.get(zone, "#eeeeee")
        alpha = 0.85 if strat != "zero" else 0.30
        ax.bar(dom, 1, width=0.9, color=color, alpha=alpha,
               edgecolor="white", linewidth=0.3, zorder=2)

    ax.set_yticks([])
    ax.set_ylabel("Physical zone", fontsize=FONT_LABEL)
    ax.set_xlabel("CASCADE domain  (1 = Buxton  →  90 = Rodanthe)", fontsize=FONT_LABEL)
    ax.set_ylim(0, 1); ax.tick_params(labelsize=FONT_TICK)

    ax.legend(handles=[mpatches.Patch(color="#aaaaaa", alpha=0.3,
                                      label="Faded = no correction applied")],
              fontsize=FONT_ANNOT, loc="lower right", framealpha=0.85)

    fig.suptitle(
        "Hatteras Island  |  Final BE source/sink rates\n"
        "Corrections applied only where residual is significant and physically motivated",
        fontsize=FONT_SUPTITLE, fontweight="bold")
    fig.tight_layout(rect=[0, 0, 1, 0.97])

    # Direct in-place zone labels, fitted to each zone's own width — done
    # AFTER tight_layout so the pixel-width measurements used to size each
    # label reflect the figure's final axes geometry, not a stale pre-layout
    # one.
    zone_runs = find_zone_runs(results["physical_zone"], domains)
    label_zone_runs(ax, fig, zone_runs)

    fig.savefig(out_path, dpi=180, bbox_inches="tight")
    plt.close(fig)
    print(f"  BE rates figure saved → {out_path}")


# ============================================================
# PRINT DOMAIN_BE_RATES DICTS
# ============================================================

def print_be_dicts(results, txt_path=None):
    """Print ready-to-paste DOMAIN_BE_RATES dicts for all scenarios
    and optionally write to a txt file."""
    scenarios = {
        "P1 hindcast":         "be_hindcast_p1",
        "P2 hindcast":         "be_hindcast_p2",
        "Forecast: continue":  "be_forecast_continue",
        "Forecast: revert":    "be_forecast_revert",
        "Forecast: neutral":   "be_forecast_neutral",
    }

    lines = []
    lines.append("=" * 70)
    lines.append("DOMAIN_BE_RATES — ready to paste into CASCADE run script")
    lines.append("=" * 70)

    for label, col in scenarios.items():
        nonzero = (results[col] != 0).sum()
        lines.append(f"")
        lines.append(f"# {label}  ({nonzero} domains with non-zero correction)")
        lines.append("DOMAIN_BE_RATES = {")
        for dom in results.index:
            val = results.loc[dom, col]
            strat = results.loc[dom, "strategy"]
            zone  = results.loc[dom, "physical_zone"]
            if strat == "locked":
                lines.append(f"    {dom:3d}: 0.0,  # LOCKED — use your solved value, not 0.0")
            elif val != 0.0:
                lines.append(f"    {dom:3d}: {val:+.1f},  # {zone}")
            else:
                lines.append(f"    {dom:3d}: 0.0,")
        lines.append("}")

    # Print to console
    print("\n" + "\n".join(lines))

    # Write to txt file
    if txt_path is not None:
        with open(txt_path, "w") as f:
            f.write("\n".join(lines) + "\n")
        print(f"\n  BE rates txt saved → {txt_path}")


# ============================================================
# MAIN
# ============================================================

def main():
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    # ── Load CoastSat LRRs ────────────────────────────────────────────────────
    print("Loading CoastSat LRRs …")
    cs_p1_ser = pd.read_csv(P1_COASTSAT_CSV).set_index(DOMAIN_COL)[LRR_COL]
    cs_p2_ser = pd.read_csv(P2_COASTSAT_CSV).set_index(DOMAIN_COL)[LRR_COL]
    cs_p1 = np.array([cs_p1_ser.get(d, np.nan) for d in range(1, NUM_REAL_DOMAINS+1)])
    cs_p2 = np.array([cs_p2_ser.get(d, np.nan) for d in range(1, NUM_REAL_DOMAINS+1)])
    print(f"  P1: {np.sum(~np.isnan(cs_p1))} valid domains")
    print(f"  P2: {np.sum(~np.isnan(cs_p2))} valid domains")

    # ── Load CASCADE LRRs ─────────────────────────────────────────────────────
    print("\nLoading CASCADE base-run LRRs …")
    casc_p1_ser = load_cascade_lrr(P1_CASCADE_NPZ, P1_START, P1_END)
    casc_p2_ser = load_cascade_lrr(P2_CASCADE_NPZ, P2_START, P2_END)
    casc_p1 = np.array([casc_p1_ser.get(d, np.nan) for d in range(1, NUM_REAL_DOMAINS+1)])
    casc_p2 = np.array([casc_p2_ser.get(d, np.nan) for d in range(1, NUM_REAL_DOMAINS+1)])

    domains = np.arange(1, NUM_REAL_DOMAINS + 1)
    pd.DataFrame({"domain": domains, "casc_p1": casc_p1, "casc_p2": casc_p2}
                 ).to_csv(os.path.join(OUTPUT_DIR, "cascade_base_lrr.csv"), index=False)

    # ── Smooth the OBSERVED shoreline change rate (not the residual) ─────────
    # Updated methodology: LOESS smoothing is applied to the CoastSat rate
    # itself, before differencing against CASCADE, rather than to the
    # residual afterward. Domains 1-GROIN_EXCLUDE_THROUGH_DOMAIN (Buxton
    # groin influence zone) are excluded from the fit and pass through with
    # their raw, unsmoothed rate — see smooth_shoreline_rate().
    print(f"\nSmoothing observed shoreline change rate "
          f"(excluding groin zone D1-{GROIN_EXCLUDE_THROUGH_DOMAIN}) …")
    cs_p1_smooth = smooth_shoreline_rate(cs_p1)
    cs_p2_smooth = smooth_shoreline_rate(cs_p2)

    # ── Raw residual — fully unsmoothed on both sides, kept for diagnostic
    # comparison only. It no longer drives any decision below. ────────────────
    raw_p1 = cs_p1 - casc_p1
    raw_p2 = cs_p2 - casc_p2
    print(f"  Mean |raw residual| P1: {np.nanmean(np.abs(raw_p1)):.2f} m/yr")
    print(f"  Mean |raw residual| P2: {np.nanmean(np.abs(raw_p2)):.2f} m/yr")

    # ── Residual from the smoothed observed rate — THIS is what drives zone
    # identification, significance testing, and BE corrections from here on ──
    print("Computing residual from smoothed shoreline rate …")
    smooth_p1 = cs_p1_smooth - casc_p1
    smooth_p2 = cs_p2_smooth - casc_p2

    # ── Compute BE corrections ────────────────────────────────────────────────
    print("Computing BE corrections …")
    results = compute_be_rates(raw_p1, raw_p2, smooth_p1, smooth_p2)

    # Summary
    n_zero     = (results["strategy"] == "zero").sum()
    n_stable   = (results["strategy"] == "stable").sum()
    n_shifting = (results["strategy"] == "shifting").sum()
    print(f"\n  No correction needed:          {n_zero:3d} domains")
    print(f"  Stable correction (single BE): {n_stable:3d} domains")
    print(f"  Shifting (period-specific):    {n_shifting:3d} domains")

    # ── Save CSV ──────────────────────────────────────────────────────────────
    csv_out = os.path.join(OUTPUT_DIR, "be_zone_metrics.csv")
    results.to_csv(csv_out)
    print(f"\n  Metrics CSV saved → {csv_out}")

    # ── Figures ───────────────────────────────────────────────────────────────
    txt_out = os.path.join(OUTPUT_DIR, "DOMAIN_BE_RATES.txt")
    print_be_dicts(results, txt_path=txt_out)

    print("\nGenerating figures …")
    plot_diagnostic(
        cs_p1, cs_p2, casc_p1, casc_p2,
        cs_p1_smooth, cs_p2_smooth,
        raw_p1, raw_p2, smooth_p1, smooth_p2,
        results,
        os.path.join(OUTPUT_DIR, "fig_be_diagnostic.png"))

    plot_be_rates(
        results,
        os.path.join(OUTPUT_DIR, "fig_be_rates.png"))

    # ── Print dicts ───────────────────────────────────────────────────────────
    print("\nDone.")


if __name__ == "__main__":
    main()
