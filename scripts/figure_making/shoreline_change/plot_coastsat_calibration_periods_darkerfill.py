"""
CoastSat Observed Shoreline Change Rates - Calibration & Test Periods
Hatteras Island, NC

Produces a single publication-quality figure for poster use showing
domain-averaged LRR shoreline change rates for:
  - Calibration Period: 1984–2004  (CoastSat satellite-derived)
  - Test Period:        2004–2024  (CoastSat satellite-derived)

Rates are pre-computed LRR values averaged across all CoastSat transects
within each 500-m CASCADE model domain.
Data source: CoastSat satellite-derived shoreline analysis

Usage:
    Place this script in the same directory as both domain summary CSVs and run:
    python plot_coastsat_calibration_periods_darkerfill.py
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from matplotlib.transforms import blended_transform_factory

# ============================================================================
# CONFIGURATION — edit these if your file names or column names differ
# ============================================================================

# Input CSVs — pre-aggregated domain-level LRR summaries
CSV_P1     = r'C:\Users\hanna\PycharmProjects\CASCADE\scripts\input_preperation\CoastSat\1984_2004\domain_lrr_1984_2004_summary.csv'
CSV_P2     = r'C:\Users\hanna\PycharmProjects\CASCADE\scripts\input_preperation\CoastSat\2004_2024\domain_lrr_2004_2024_summary.csv'

DOMAIN_COL = 'domain_number'
LRR_COL    = 'mean_lrr'   # Column containing domain-averaged LRR (m/yr)

# Period labels
PERIOD_1_LABEL = '1984_2004'
PERIOD_2_LABEL = '2004-2024'

DOMAIN_MIN = 1
DOMAIN_MAX = 90

# --- Geographic annotations ---

# Community shaded spans: label → (domain_lo, domain_hi), inclusive
TOWN_SPANS = {
    "Buxton":      (7,  8),
    "Avon":        (21, 31),
    "Tri-Village": (68, 83),   # Salvo / Waves / Rodanthe combined area
}

# Individual village center lines within Tri-Village
VILLAGE_LINES = {
    "Salvo":    69,
    "Waves":    74,
    "Rodanthe": 80,
}

# Pier reference lines
PIERS = {
    "Avon Pier":     26,
    "Rodanthe Pier": 79,
}

# Groin reference lines
GROINS = {
    "Buxton Groin": 6,
}

# Wimble Shoals offshore shoal influence zone: (domain_lo, domain_hi)
# Spans domains 60–74; overlaps with the southern portion of Tri-Village.
WIMBLE_SHOALS = (60, 74)

# --- Annotation colors (matched to reference script) ---
C_TOWN_SPAN    = "#90AFC5"   # steel blue  — community spans
C_WIMBLE       = "#E0A800"   # amber       — Wimble Shoals fill
C_VILLAGE_LINE = "0.40"      # dark gray   — Salvo / Waves / Rodanthe lines
C_PIER         = "#1565C0"   # medium blue — pier lines
C_GROIN        = "#B71C1C"   # dark red    — groin lines

# Output
OUTPUT_FILE = 'coastsat_calibration_periods_poster.png'

# Color scheme — matched to reference script (smoothing_vs_cascade_final_fixed.py)
# Period 1 (1984–2004): dark blue  |  Period 2 (2004–2024): dark brown-red
COLOR_P1   = '#1F4E79'   # dark blue  (C_CS_1984 from reference)
COLOR_P2   = '#C0392B'   # red
COLOR_ZERO = '#2c2c2c'

# ============================================================================
# HELPER FUNCTIONS
# ============================================================================

def load_csv_safe(path):
    """Load CSV, stripping any filename prefix accidentally prepended to the header.
    Handles a known export artifact where the filename is concatenated onto
    the first column name (e.g. 'domain_lrr_summary.csvdomain_number')."""
    df = pd.read_csv(path)
    df.columns = [c.split('csv')[-1] if 'csv' in c else c for c in df.columns]
    return df


def prepare_domain_rates(df, label):
    """Clip to real island domains, drop NaNs, and sort by domain number."""
    df = df[[DOMAIN_COL, LRR_COL]].dropna().copy()
    df = df[(df[DOMAIN_COL] >= DOMAIN_MIN) & (df[DOMAIN_COL] <= DOMAIN_MAX)]
    df = df.sort_values(DOMAIN_COL).reset_index(drop=True)
    print(f"  {label}: {len(df)} domains  |  "
          f"range {df[DOMAIN_COL].min():.0f}-{df[DOMAIN_COL].max():.0f}  |  "
          f"island mean {df[LRR_COL].mean():+.2f} m/yr")
    return df


def add_annotations(ax):
    """
    Add all geographic reference annotations to an axis.

    Layer order (bottom → top):
      1. Wimble Shoals influence zone  (hatched amber fill, bottom label)
      2. Community shaded spans        (steel-blue fill, top labels)
      3. Village center lines          (dashed gray,     y=0.88)
      4. Pier lines                    (dash-dot blue,   y=0.76, rotated)
      5. Groin lines                   (dotted red,      y=0.76, rotated)

    All label y-positions use blended axes-fraction coordinates so they
    stay fixed relative to the panel height regardless of data range.
    """
    # Blended transform: data x, axes-fraction y
    trans = blended_transform_factory(ax.transData, ax.transAxes)

    # ------------------------------------------------------------------
    # 1. Wimble Shoals influence zone
    #    Placed first so community spans render on top of it.
    #    Label at the bottom of the panel to avoid crowding the top tier.
    # ------------------------------------------------------------------
    wlo, whi = WIMBLE_SHOALS
    ax.axvspan(wlo - 0.5, whi + 0.5,
               color=C_WIMBLE, alpha=0.10, zorder=0,
               hatch="///", edgecolor=C_WIMBLE, linewidth=0)
    ax.text((wlo + whi) / 2.0, 0.04,
            "Wimble Shoals\ninfluence", transform=trans,
            ha="center", va="bottom", fontsize=7, color="#7A5800",
            style="italic",
            bbox=dict(boxstyle="round,pad=0.2", fc="white", ec="none", alpha=0.80))

    # ------------------------------------------------------------------
    # 2. Community / town spans
    #    Labels at y=0.90 (near panel top, clear of the data lines).
    # ------------------------------------------------------------------
    for span_label, (d_lo, d_hi) in TOWN_SPANS.items():
        ax.axvspan(d_lo - 0.5, d_hi + 0.5,
                   color=C_TOWN_SPAN, alpha=0.14, zorder=0)
        ax.text((d_lo + d_hi) / 2.0, 0.90,
                span_label, transform=trans,
                ha="center", va="top", fontsize=8, color="0.25",
                fontweight="bold",
                bbox=dict(boxstyle="round,pad=0.2", fc="white", ec="none", alpha=0.85))

    # ------------------------------------------------------------------
    # 3. Village center lines (within Tri-Village span)
    #    Dashed gray; labels at y=0.79 so they sit clearly below span labels.
    # ------------------------------------------------------------------
    for vname, dom in VILLAGE_LINES.items():
        ax.axvline(dom, color=C_VILLAGE_LINE, lw=0.9, ls="--",
                   alpha=0.65, zorder=1)
        ax.text(dom, 0.84, vname, transform=trans,
                ha="center", va="top", fontsize=7.5, color="0.30",
                bbox=dict(boxstyle="round,pad=0.15", fc="white",
                          ec="none", alpha=0.80))

    # ------------------------------------------------------------------
    # 4. Pier lines
    #    Dash-dot blue; rotated labels at y=0.74.
    # ------------------------------------------------------------------
    for pname, dom in PIERS.items():
        ax.axvline(dom, color=C_PIER, lw=1.0, ls="-.", alpha=0.80, zorder=2)
        ax.text(dom, 0.76, pname, transform=trans,
                ha="center", va="top", fontsize=7, color=C_PIER,
                rotation=90,
                bbox=dict(boxstyle="round,pad=0.15", fc="white",
                          ec="none", alpha=0.80))

    # ------------------------------------------------------------------
    # 5. Groin lines
    #    Dotted dark-red; rotated labels at y=0.74.
    # ------------------------------------------------------------------
    for gname, dom in GROINS.items():
        ax.axvline(dom, color=C_GROIN, lw=1.1, ls=":", alpha=0.85, zorder=2)
        ax.text(dom, 0.76, gname, transform=trans,
                ha="center", va="top", fontsize=7, color=C_GROIN,
                rotation=90,
                bbox=dict(boxstyle="round,pad=0.15", fc="white",
                          ec="none", alpha=0.80))


def annotation_legend_handles():
    """
    Return proxy artists explaining the annotation layer types.
    Append these to the legend handle list so readers can decode
    all reference marks without scanning every label individually.
    """
    return [
        Patch(fc=C_TOWN_SPAN,  alpha=0.30, label="Community"),
        Patch(fc=C_WIMBLE,     alpha=0.25, hatch="///",
              edgecolor=C_WIMBLE, linewidth=0, label="Wimble Shoals influence"),
        Line2D([0], [0], color=C_VILLAGE_LINE, lw=0.9, ls="--",
               label="Village center"),
        Line2D([0], [0], color=C_PIER,  lw=1.0, ls="-.", label="Pier"),
        Line2D([0], [0], color=C_GROIN, lw=1.1, ls=":",  label="Groin"),
    ]


# ============================================================================
# LOAD & PROCESS DATA
# ============================================================================

print("Loading CoastSat domain LRR data...")

df_p1 = load_csv_safe(CSV_P1)
df_p2 = load_csv_safe(CSV_P2)

df_p1 = prepare_domain_rates(df_p1, f'Period 1 ({PERIOD_1_LABEL})')
df_p2 = prepare_domain_rates(df_p2, f'Period 2 ({PERIOD_2_LABEL})')

domains_p1 = df_p1[DOMAIN_COL].values
p1         = df_p1[LRR_COL].values

domains_p2 = df_p2[DOMAIN_COL].values
p2         = df_p2[LRR_COL].values

# ============================================================================
# PLOT
# ============================================================================

fig, ax = plt.subplots(figsize=(13, 6.0))
fig.patch.set_facecolor('white')
ax.set_facecolor('white')

# ------------------------------------------------------------------
# Annotations (drawn first so data lines render on top)
# ------------------------------------------------------------------
add_annotations(ax)

# ------------------------------------------------------------------
# Shaded fill under curves toward zero line
# ------------------------------------------------------------------
ax.fill_between(domains_p1, p1, 0,
                where=(p1 < 0),  alpha=0.20, color=COLOR_P1, interpolate=True)
ax.fill_between(domains_p1, p1, 0,
                where=(p1 >= 0), alpha=0.18, color=COLOR_P1, interpolate=True)
ax.fill_between(domains_p2, p2, 0,
                where=(p2 < 0),  alpha=0.20, color=COLOR_P2, interpolate=True)
ax.fill_between(domains_p2, p2, 0,
                where=(p2 >= 0), alpha=0.18, color=COLOR_P2, interpolate=True)

# ------------------------------------------------------------------
# Main lines
# ------------------------------------------------------------------
ax.plot(domains_p1, p1, color=COLOR_P1, linewidth=2.8, zorder=5,
        label=f'Period 1: {PERIOD_1_LABEL}  (calibration)')
ax.plot(domains_p2, p2, color=COLOR_P2, linewidth=2.8, zorder=5,
        label=f'Period 2: {PERIOD_2_LABEL}  (test)')

# Zero reference line
ax.axhline(0, color=COLOR_ZERO, linewidth=1.2, linestyle='--', alpha=0.65, zorder=3)

# ------------------------------------------------------------------
# Axes formatting
# ------------------------------------------------------------------
ax.set_xlim(DOMAIN_MIN - 0.5, DOMAIN_MAX + 0.5)
ax.xaxis.set_major_locator(ticker.MultipleLocator(10))
ax.xaxis.set_minor_locator(ticker.MultipleLocator(5))
ax.yaxis.set_major_locator(ticker.MultipleLocator(1))
ax.tick_params(axis='both', which='major', labelsize=11, direction='in', length=5)
ax.tick_params(axis='both', which='minor', direction='in', length=3)
ax.grid(True, which='major', linestyle=':', linewidth=0.6, alpha=0.4, color='gray')
ax.spines[['top', 'right']].set_visible(False)
ax.spines[['left', 'bottom']].set_linewidth(1.2)

# Lock ylim to data range before placing accretion/erosion labels
all_vals  = np.concatenate([p1, p2])
ymin_data = all_vals.min()
ymax_data = all_vals.max()
pad       = (ymax_data - ymin_data) * 0.06
ax.set_ylim(ymin_data - pad, ymax_data + pad)

# ------------------------------------------------------------------
# Accretion / Erosion side labels anchored relative to zero line
# ------------------------------------------------------------------
ybot, ytop     = ax.get_ylim()
zero_frac      = (0 - ybot) / (ytop - ybot)
accretion_frac = zero_frac + (1 - zero_frac) / 2
erosion_frac   = zero_frac / 2

ax.text(1.0, accretion_frac, 'Accretion ▲', transform=ax.transAxes,
        fontsize=9.5, color='#555555', ha='right', va='center', style='italic')
ax.text(1.0, erosion_frac,   'Erosion ▼',   transform=ax.transAxes,
        fontsize=9.5, color='#555555', ha='right', va='center', style='italic')

# ------------------------------------------------------------------
# Axis labels
# ------------------------------------------------------------------
ax.set_xlabel('CASCADE Model Domain (500 m alongshore)', fontsize=12,
              fontweight='bold', labelpad=8)
ax.set_ylabel('Shoreline Change Rate (m/yr)', fontsize=12,
              fontweight='bold', labelpad=8)

# ------------------------------------------------------------------
# North / South geographic orientation at top of plot
# ------------------------------------------------------------------
ax.text(0.0, 1.01, '← S  |  Cape Hatteras',
        transform=ax.transAxes,
        fontsize=9, color='#444444', ha='left', va='bottom',
        style='italic', clip_on=False)
ax.text(1.0, 1.01, 'Pea Island  |  N →',
        transform=ax.transAxes,
        fontsize=9, color='#444444', ha='right', va='bottom',
        style='italic', clip_on=False)

# ------------------------------------------------------------------
# Title
# ------------------------------------------------------------------
ax.set_title(
    'CoastSat Shoreline Change Rates — Hatteras Island, North Carolina, USA',
    fontsize=13.5, fontweight='bold', pad=12, color='#1a2a3a'
)

# ------------------------------------------------------------------
# Legend — data series + annotation proxies combined
# ------------------------------------------------------------------
data_handles = [
    Line2D([0], [0], color=COLOR_P1, linewidth=2.8,
           label=f'Period 1: {PERIOD_1_LABEL}  (calibration)'),
    Line2D([0], [0], color=COLOR_P2, linewidth=2.8,
           label=f'Period 2: {PERIOD_2_LABEL}  (test)'),
]
all_handles = data_handles + annotation_legend_handles()
ax.legend(handles=all_handles,
          loc='lower center', bbox_to_anchor=(0.5, 0.02),
          fontsize=9.5, framealpha=0.95, edgecolor='#cccccc',
          frameon=True, ncol=4)

# ------------------------------------------------------------------
# Data source caption
# ------------------------------------------------------------------
fig.text(
    0.012, 0.005,
    'Rates calculated as Linear Regression Rate (LRR) averaged across all CoastSat transects within each 500-m CASCADE model domain. '
    'Shoreline data: CoastSat satellite-derived shorelines (Vos et al., 2019).',
    fontsize=8, color='#666666', ha='left', va='bottom', style='italic',
)

plt.tight_layout(rect=[0, 0.04, 1, 1])
plt.savefig(OUTPUT_FILE, dpi=300, bbox_inches='tight', facecolor='white')
print(f"\nSaved: {OUTPUT_FILE}")
plt.close()
