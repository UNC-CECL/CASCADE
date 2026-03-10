"""
DSAS Observed Shoreline Change Rates - Calibration & Test Periods
Hatteras Island, NC

Produces a single publication-quality figure for poster use showing
domain-averaged EPR shoreline change rates for:
  - Calibration Period 1: 1978–1997  (shorelines: 1978, 1997)
  - Test Period:          1997–2019  (shorelines: 1997, 2019)

Rates are averaged across all DSAS transects within each 500-m CASCADE model domain.
Data source: USGS Digital Shoreline Analysis System (DSAS)

Usage:
    Place this script in the same directory as your DSAS intersection CSV and run:
    python plot_dsas_calibration_periods.py
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from matplotlib.lines import Line2D

# ============================================================================
# CONFIGURATION — edit these if your column names or file differ
# ============================================================================

INPUT_CSV = 'All_Shoreline_Transect_Intersections.csv'

TRANSECT_ID_COL  = 'Transects_100m_LineID'
DOMAIN_ID_COL    = 'Transects_100m_AddSpatialJoin_domain_id'
YEAR_COL         = 'Year'
DISTANCE_COL     = 'NEAR_DIST'

# Shoreline years used in each period
PERIOD_1_YEARS   = (1978, 1997)   # Calibration period 1 (uses: 1978, 1987, 1997)
PERIOD_2_YEARS   = (1997, 2019)   # Test period          (uses: 1997, 2009, 2019)

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

OUTPUT_FILE      = 'dsas_calibration_periods_poster.png'

# Poster color scheme — matched to your UNC teal poster
COLOR_P1         = '#1a6b8a'   # Dark teal  (calibration)
COLOR_P2         = '#c0392b'   # Deep red   (test period)
COLOR_ZERO       = '#2c2c2c'
FILL_P1          = '#aed6e8'
FILL_P2          = '#f5b7b1'

# ============================================================================
# LOAD & PROCESS DATA
# ============================================================================

print("Loading DSAS intersection data...")
df = pd.read_csv(INPUT_CSV)
df_clean = df[[TRANSECT_ID_COL, DOMAIN_ID_COL, YEAR_COL, DISTANCE_COL]].dropna()

# Save domain map BEFORE groupby drops the domain column
domain_map = df_clean[[TRANSECT_ID_COL, DOMAIN_ID_COL]].drop_duplicates()

# Average duplicate intersections (shoreline crossing transect twice)
df_clean = df_clean.groupby(
    [TRANSECT_ID_COL, YEAR_COL], as_index=False
)[DISTANCE_COL].mean()

# Pivot to wide format: rows = transect, columns = year
transect_wide = df_clean.pivot(
    index=TRANSECT_ID_COL, columns=YEAR_COL, values=DISTANCE_COL
)

# Calculate end-point rates (EPR) for each period
# Sign convention: negative = erosion, positive = accretion
def calc_epr(wide, y1, y2):
    return -1 * (wide[y2] - wide[y1]) / (y2 - y1)

y1a, y1b = PERIOD_1_YEARS
y2a, y2b = PERIOD_2_YEARS

transect_rates = pd.DataFrame(index=transect_wide.index)
transect_rates['EPR_P1'] = calc_epr(transect_wide, y1a, y1b)
transect_rates['EPR_P2'] = calc_epr(transect_wide, y2a, y2b)

# Map transects to domains and compute domain means
rates_with_domain = transect_rates.merge(
    domain_map, left_index=True, right_on=TRANSECT_ID_COL, how='left'
)
domain_rates = rates_with_domain.groupby(DOMAIN_ID_COL)[['EPR_P1', 'EPR_P2']].mean().reset_index()
domain_rates.columns = ['Domain', 'EPR_P1', 'EPR_P2']
domain_rates = domain_rates.sort_values('Domain')

# Clip to real domains (removes buffer domains at island ends)
domain_rates = domain_rates[
    (domain_rates['Domain'] >= DOMAIN_MIN) & (domain_rates['Domain'] <= DOMAIN_MAX)
].copy()

print(f"  Domains with data: {len(domain_rates)}")
print(f"  Domain range: {domain_rates['Domain'].min()} – {domain_rates['Domain'].max()}")
print(f"\n  Period 1 ({y1a}–{y1b})  island mean: {domain_rates['EPR_P1'].mean():+.2f} m/yr")
print(f"  Period 2 ({y2a}–{y2b})  island mean: {domain_rates['EPR_P2'].mean():+.2f} m/yr")

# ============================================================================
# PLOT
# ============================================================================

fig, ax = plt.subplots(figsize=(13, 6.0))
fig.patch.set_facecolor('white')
ax.set_facecolor('white')

domains = domain_rates['Domain'].values
p1      = domain_rates['EPR_P1'].values
p2      = domain_rates['EPR_P2'].values

# Shaded fill under curves toward zero line
ax.fill_between(domains, p1, 0,
                where=(p1 < 0), alpha=0.10, color=COLOR_P1, interpolate=True)
ax.fill_between(domains, p1, 0,
                where=(p1 >= 0), alpha=0.08, color=COLOR_P1, interpolate=True)
ax.fill_between(domains, p2, 0,
                where=(p2 < 0), alpha=0.09, color=COLOR_P2, interpolate=True)
ax.fill_between(domains, p2, 0,
                where=(p2 >= 0), alpha=0.07, color=COLOR_P2, interpolate=True)

# Main lines
ax.plot(domains, p1, color=COLOR_P1, linewidth=2.8, zorder=5,
        label=f'Period 1: {y1a}–{y1b}  (calibration)')
ax.plot(domains, p2, color=COLOR_P2, linewidth=2.8, zorder=5,
        label=f'Period 2: {y2a}–{y2b}  (test)')

# Zero line
ax.axhline(0, color=COLOR_ZERO, linewidth=1.2, linestyle='--', alpha=0.65, zorder=3)

# ---- Axes formatting ----
ax.set_xlim(DOMAIN_MIN - 0.5, DOMAIN_MAX + 0.5)
ax.xaxis.set_major_locator(ticker.MultipleLocator(10))
ax.xaxis.set_minor_locator(ticker.MultipleLocator(5))
ax.yaxis.set_major_locator(ticker.MultipleLocator(1))
ax.tick_params(axis='both', which='major', labelsize=11, direction='in', length=5)
ax.tick_params(axis='both', which='minor', direction='in', length=3)
ax.grid(True, which='major', linestyle=':', linewidth=0.6, alpha=0.4, color='gray')
ax.spines[['top', 'right']].set_visible(False)
ax.spines[['left', 'bottom']].set_linewidth(1.2)

# ---- Geographic place labels — positioned in axes coords so ylim is unaffected ----
# Stagger Salvo/Rodanthe since they're only 2 domains apart
label_y_axes = {'Buxton': 0.97, 'Avon': 0.97, 'Salvo': 0.97, 'Waves': 0.85, 'Rodanthe': 0.97}

# Lock ylim to actual data range before placing labels
ymin_data = min(p1.min(), p2.min())
ymax_data = max(p1.max(), p2.max())
pad = (ymax_data - ymin_data) * 0.06
ax.set_ylim(ymin_data - pad, ymax_data + pad)

xmin_d, xmax_d = ax.get_xlim()
xrange_d = xmax_d - xmin_d

for domain, place in GEO_LABELS:
    # Convert domain number to axes-fraction x coordinate
    x_frac = (domain - xmin_d) / xrange_d

    # Dashed vertical line clipped within plot
    ax.axvline(domain, color='#444444', linewidth=1.0, linestyle='--',
               alpha=0.5, zorder=2, clip_on=True)

    # Label in axes coordinates — won't affect ylim
    y_frac = label_y_axes.get(place, 0.95)
    ax.text(x_frac, y_frac, place,
            transform=ax.transAxes,
            fontsize=10, fontweight='bold', color='#222222',
            ha='center', va='top', clip_on=False,
            bbox=dict(boxstyle='round,pad=0.25', facecolor='white',
                      edgecolor='#999999', linewidth=0.8, alpha=0.92))

# ---- Erosion / Accretion side labels anchored to zero line ----
ybot, ytop = ax.get_ylim()
zero_frac = (0 - ybot) / (ytop - ybot)
# Place labels equidistant from zero, at the midpoint of each half
accretion_frac = zero_frac + (1 - zero_frac) / 2
erosion_frac   = zero_frac / 2

ax.text(1.0, accretion_frac, 'Accretion ▲', transform=ax.transAxes,
        fontsize=9.5, color='#555555', ha='right', va='center', style='italic')
ax.text(1.0, erosion_frac,   'Erosion ▼',   transform=ax.transAxes,
        fontsize=9.5, color='#555555', ha='right', va='center', style='italic')

ax.set_xlabel('CASCADE Model Domain (500 m alongshore)', fontsize=12, fontweight='bold', labelpad=8)
ax.set_ylabel('Shoreline Change Rate (m/yr)', fontsize=12, fontweight='bold', labelpad=8)

# ---- North/South geographic orientation at top of plot ----
ax.text(0.0, 1.01, '← S  |  Cape Hatteras',
        transform=ax.transAxes,
        fontsize=9, color='#444444', ha='left', va='bottom',
        style='italic', clip_on=False)
ax.text(1.0, 1.01, 'Pea Island  |  N →',
        transform=ax.transAxes,
        fontsize=9, color='#444444', ha='right', va='bottom',
        style='italic', clip_on=False)

# ---- Title ----
ax.set_title(
    'Observed Shoreline Change Rates — Hatteras Island, North Carolina, USA',
    fontsize=13.5, fontweight='bold', pad=12, color='#1a2a3a'
)

# ---- Legend ----
legend_elements = [
    Line2D([0], [0], color=COLOR_P1, linewidth=2.8,
           label='Period 1: 1978–1997  [shorelines: 1978, 1987, 1997]  (calibration)'),
    Line2D([0], [0], color=COLOR_P2, linewidth=2.8,
           label='Period 2: 1997–2019  [shorelines: 1997, 2009, 2019]  (test)'),
]
ax.legend(handles=legend_elements, loc='lower center', bbox_to_anchor=(0.5, 0.02),
          fontsize=10, framealpha=0.95, edgecolor='#cccccc', frameon=True)

# ---- Data source caption ----
fig.text(
    0.012, 0.005,
    'Rates calculated as Linear Regression Rate (LRR) averaged across all 100-m DSAS transects within each 500-m model domain. '
    'Shoreline data: Hapke & Henderson (2015, USGS OFR 2015-1002); Kratzmann et al. (2017, USGS data release, doi:10.5066/F74X55X7).',
    fontsize=12, color='#666666', ha='left', va='bottom', style='italic',
    wrap=True
)

plt.tight_layout(rect=[0, 0.04, 1, 1])
plt.savefig(OUTPUT_FILE, dpi=300, bbox_inches='tight', facecolor='white')
print(f"\nSaved: {OUTPUT_FILE}")
plt.close()
