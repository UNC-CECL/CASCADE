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
    python plot_coastsat_calibration_periods.py
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from matplotlib.lines import Line2D

# ============================================================================
# CONFIGURATION — edit these if your file names or column names differ
# ============================================================================

# Input CSVs — pre-aggregated domain-level LRR summaries
CSV_P1           = r'C:\Users\hanna\PycharmProjects\CASCADE\scripts\input_preperation\CoastSat\1984_2004\domain_lrr_1984_2004_summary.csv'   # Calibration period
CSV_P2           = r'C:\Users\hanna\PycharmProjects\CASCADE\scripts\input_preperation\CoastSat\2004_2024\domain_lrr_2004_2024_summary.csv'   # Test period

DOMAIN_COL       = 'domain_number'
LRR_COL          = 'mean_lrr'        # Column containing domain-averaged LRR (m/yr)

# Period labels
PERIOD_1_LABEL   = '1984-2004'
PERIOD_2_LABEL   = '2004-2024'

DOMAIN_MIN       = 1
DOMAIN_MAX       = 90

# Geographic place labels: (domain_number, label_string)
GEO_LABELS       = [
    (9,  'Buxton'),
    (27, 'Avon'),
    (69, 'Salvo'),
    (75, 'Waves'),
    (79, 'Rodanthe'),
]

OUTPUT_FILE      = '../coastsat_calibration_periods_poster.png'

# Color scheme -- matched to UNC teal poster palette
# Period 1 (calibration): dark teal (consistent with original DSAS plot)
# Period 2 (test):        amber/gold -- distinct and readable on white
COLOR_P1         = '#1a6b8a'   # Dark teal
COLOR_P2         = '#d47c0f'   # Amber/gold
COLOR_ZERO       = '#2c2c2c'

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

# Shaded fill under curves toward zero line
ax.fill_between(domains_p1, p1, 0,
                where=(p1 < 0), alpha=0.10, color=COLOR_P1, interpolate=True)
ax.fill_between(domains_p1, p1, 0,
                where=(p1 >= 0), alpha=0.08, color=COLOR_P1, interpolate=True)
ax.fill_between(domains_p2, p2, 0,
                where=(p2 < 0), alpha=0.09, color=COLOR_P2, interpolate=True)
ax.fill_between(domains_p2, p2, 0,
                where=(p2 >= 0), alpha=0.07, color=COLOR_P2, interpolate=True)

# Main lines
ax.plot(domains_p1, p1, color=COLOR_P1, linewidth=2.8, zorder=5,
        label=f'Period 1: {PERIOD_1_LABEL}  (calibration)')
ax.plot(domains_p2, p2, color=COLOR_P2, linewidth=2.8, zorder=5,
        label=f'Period 2: {PERIOD_2_LABEL}  (test)')

# Zero line
ax.axhline(0, color=COLOR_ZERO, linewidth=1.2, linestyle='--', alpha=0.65, zorder=3)

# Axes formatting
ax.set_xlim(DOMAIN_MIN - 0.5, DOMAIN_MAX + 0.5)
ax.xaxis.set_major_locator(ticker.MultipleLocator(10))
ax.xaxis.set_minor_locator(ticker.MultipleLocator(5))
ax.yaxis.set_major_locator(ticker.MultipleLocator(1))
ax.tick_params(axis='both', which='major', labelsize=11, direction='in', length=5)
ax.tick_params(axis='both', which='minor', direction='in', length=3)
ax.grid(True, which='major', linestyle=':', linewidth=0.6, alpha=0.4, color='gray')
ax.spines[['top', 'right']].set_visible(False)
ax.spines[['left', 'bottom']].set_linewidth(1.2)

# Lock ylim to data range before placing labels
all_vals  = np.concatenate([p1, p2])
ymin_data = all_vals.min()
ymax_data = all_vals.max()
pad       = (ymax_data - ymin_data) * 0.06
ax.set_ylim(ymin_data - pad, ymax_data + pad)

# Geographic place labels -- positioned in axes coords so ylim is unaffected
# Waves staggered lower because it is only 2 domains from Salvo
label_y_axes = {
    'Buxton':   0.97,
    'Avon':     0.97,
    'Salvo':    0.97,
    'Waves':    0.85,
    'Rodanthe': 0.97,
}

xmin_d, xmax_d = ax.get_xlim()
xrange_d = xmax_d - xmin_d

for domain, place in GEO_LABELS:
    x_frac = (domain - xmin_d) / xrange_d
    ax.axvline(domain, color='#444444', linewidth=1.0, linestyle='--',
               alpha=0.5, zorder=2, clip_on=True)
    y_frac = label_y_axes.get(place, 0.95)
    ax.text(x_frac, y_frac, place,
            transform=ax.transAxes,
            fontsize=10, fontweight='bold', color='#222222',
            ha='center', va='top', clip_on=False,
            bbox=dict(boxstyle='round,pad=0.25', facecolor='white',
                      edgecolor='#999999', linewidth=0.8, alpha=0.92))

# Erosion / Accretion side labels anchored to zero line
ybot, ytop     = ax.get_ylim()
zero_frac      = (0 - ybot) / (ytop - ybot)
accretion_frac = zero_frac + (1 - zero_frac) / 2
erosion_frac   = zero_frac / 2

ax.text(1.0, accretion_frac, 'Accretion ▲', transform=ax.transAxes,
        fontsize=9.5, color='#555555', ha='right', va='center', style='italic')
ax.text(1.0, erosion_frac,   'Erosion ▼',   transform=ax.transAxes,
        fontsize=9.5, color='#555555', ha='right', va='center', style='italic')

# Axis labels
ax.set_xlabel('CASCADE Model Domain (500 m alongshore)', fontsize=12,
              fontweight='bold', labelpad=8)
ax.set_ylabel('Shoreline Change Rate (m/yr)', fontsize=12,
              fontweight='bold', labelpad=8)

# North/South geographic orientation at top of plot
ax.text(0.0, 1.01, '← S  |  Cape Hatteras',
        transform=ax.transAxes,
        fontsize=9, color='#444444', ha='left', va='bottom',
        style='italic', clip_on=False)
ax.text(1.0, 1.01, 'Pea Island  |  N →',
        transform=ax.transAxes,
        fontsize=9, color='#444444', ha='right', va='bottom',
        style='italic', clip_on=False)

# Title
ax.set_title(
    'CoastSat Shoreline Change Rates — Hatteras Island, North Carolina, USA',
    fontsize=13.5, fontweight='bold', pad=12, color='#1a2a3a'
)

# Legend
legend_elements = [
    Line2D([0], [0], color=COLOR_P1, linewidth=2.8,
           label=f'Period 1: {PERIOD_1_LABEL}  (calibration)'),
    Line2D([0], [0], color=COLOR_P2, linewidth=2.8,
           label=f'Period 2: {PERIOD_2_LABEL}  (test)'),
]
ax.legend(handles=legend_elements, loc='lower center', bbox_to_anchor=(0.5, 0.02),
          fontsize=10, framealpha=0.95, edgecolor='#cccccc', frameon=True)

# Data source caption
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
