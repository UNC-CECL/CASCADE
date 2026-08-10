"""
DSAS vs. CoastSat Shoreline Change Rate Comparison
Hatteras Island, NC

Reads comparison_table.csv (one row per domain) and produces a
publication-quality poster figure with four lines:
  - DSAS  1978–1997  (calibration)
  - DSAS  1997–2019  (test)
  - CoastSat  1978–1997  (calibration)
  - CoastSat  1997–2019  (test)

Optional shaded ±1 std bands can be toggled via SHOW_STD_BANDS.

Input CSV columns expected:
  domain
  dsas_lrr_1978_1997, dsas_std_1978_1997
  cs_lrr_1978_1997,   cs_std_1978_1997
  dsas_lrr_1997_2019, dsas_std_1997_2019
  cs_lrr_1997_2019,   cs_std_1997_2019

Usage:
    python plot_shoreline_comparison.py
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from matplotlib.lines import Line2D

# ============================================================================
# CONFIGURATION
# ============================================================================

INPUT_CSV       = r'C:\Users\hanna\PycharmProjects\CASCADE\scripts\input_preperation\CoastSat\dsas_coastsat_comparison_outputs\comparison_table.csv'
OUTPUT_FILE     = r'C:\Users\hanna\PycharmProjects\CASCADE\scripts\figure_making\plots\1978_1997_shoreline_coastsatonly.png'

DOMAIN_MIN      = 1
DOMAIN_MAX      = 90

# Toggle shaded ±1 std bands around each line
SHOW_STD_BANDS  = False

# ---- Dataset-level toggles ----
# Set SHOW_DSAS = False to hide all DSAS lines at once
# Set SHOW_COASTSAT = False to hide all CoastSat lines at once
SHOW_DSAS       = False
SHOW_COASTSAT   = True

# ---- Per-line toggles (only active when dataset toggle above is True) ----
SHOW_DSAS_P1    = True   # DSAS 1978–1997
SHOW_DSAS_P2    = True   # DSAS 1997–2019
SHOW_CS_P1      = True   # CoastSat 1978–1997
SHOW_CS_P2      = True   # CoastSat 1997–2019

# Geographic place labels: (domain_number, label_string)
GEO_LABELS = [
    (9,  'Buxton'),
    (27, 'Avon'),
    (69, 'Salvo'),
    (75, 'Waves'),
    (79, 'Rodanthe'),
]

# Color palette — matched by period; CoastSat is lighter/muted version of DSAS color
# Calibration period (1978–1997): teal family
COLOR_DSAS_P1 = '#1a6b8a'   # deep teal  (DSAS)
COLOR_CS_P1   = '#74bcd4'   # light teal (CoastSat)

# Test period (1997–2019): red family
COLOR_DSAS_P2 = '#c0392b'   # deep red   (DSAS)
COLOR_CS_P2   = '#e8938d'   # light rose (CoastSat)

# Vertical label positions (axes fraction) for geographic labels
# Stagger Salvo/Waves/Rodanthe to avoid overlap
LABEL_Y = {'Buxton': 0.97, 'Avon': 0.97, 'Salvo': 0.97,
           'Waves': 0.08, 'Rodanthe': 0.97}   # Waves pinned to bottom to clear CoastSat spike at d75–76

# ============================================================================
# LOAD DATA
# ============================================================================

print(f"Loading: {INPUT_CSV}")
df = pd.read_csv(INPUT_CSV)
df = df[(df['domain'] >= DOMAIN_MIN) & (df['domain'] <= DOMAIN_MAX)].sort_values('domain')
print(f"  {len(df)} domains loaded (range {df['domain'].min()}–{df['domain'].max()})")

# Helper: safe series → numpy (NaN where missing)
def col(name):
    return df[name].values.astype(float) if name in df.columns else np.full(len(df), np.nan)

domains       = df['domain'].values
dsas_p1       = col('dsas_lrr_1978_1997')
dsas_p1_std   = col('dsas_std_1978_1997')
dsas_p2       = col('dsas_lrr_1997_2019')
dsas_p2_std   = col('dsas_std_1997_2019')
cs_p1         = col('cs_lrr_1978_1997')
cs_p1_std     = col('cs_std_1978_1997')
cs_p2         = col('cs_lrr_1997_2019')
cs_p2_std     = col('cs_std_1997_2019')

# Print summary stats
for label, arr in [('DSAS  1978–1997', dsas_p1), ('DSAS  1997–2019', dsas_p2),
                   ('CoastSat 1978–1997', cs_p1), ('CoastSat 1997–2019', cs_p2)]:
    valid = arr[~np.isnan(arr)]
    if len(valid):
        print(f"  {label:22s}  mean={valid.mean():+.2f} m/yr  n={len(valid)}")

# ============================================================================
# PLOT
# ============================================================================

fig, ax = plt.subplots(figsize=(14, 6.5))
fig.patch.set_facecolor('white')
ax.set_facecolor('white')

def draw_line(x, y, std, color, lw, ls, label, show_band):
    """Plot a line with optional ±1 std shaded band, masking NaNs."""
    mask = ~np.isnan(y)
    if not mask.any():
        return
    ax.plot(x[mask], y[mask], color=color, linewidth=lw,
            linestyle=ls, zorder=5, label=label)
    if show_band and std is not None:
        lo = y[mask] - std[mask]
        hi = y[mask] + std[mask]
        ax.fill_between(x[mask], lo, hi, color=color, alpha=0.10, zorder=2)

show_dsas_p1 = SHOW_DSAS and SHOW_DSAS_P1
show_dsas_p2 = SHOW_DSAS and SHOW_DSAS_P2
show_cs_p1   = SHOW_COASTSAT and SHOW_CS_P1
show_cs_p2   = SHOW_COASTSAT and SHOW_CS_P2

if show_dsas_p1:
    draw_line(domains, dsas_p1, dsas_p1_std, COLOR_DSAS_P1, 2.8, '-',
              'DSAS  1978–1997  (calibration)', SHOW_STD_BANDS)

if show_dsas_p2:
    draw_line(domains, dsas_p2, dsas_p2_std, COLOR_DSAS_P2, 2.8, '-',
              'DSAS  1997–2019  (test)', SHOW_STD_BANDS)

if show_cs_p1:
    draw_line(domains, cs_p1, cs_p1_std, COLOR_CS_P1, 2.0, '--',
              'CoastSat  1978–1997  (calibration)', SHOW_STD_BANDS)

if show_cs_p2:
    draw_line(domains, cs_p2, cs_p2_std, COLOR_CS_P2, 2.0, '--',
              'CoastSat  1997–2019  (test)', SHOW_STD_BANDS)

# Zero line
ax.axhline(0, color='#2c2c2c', linewidth=1.2, linestyle='--', alpha=0.65, zorder=3)

# ---- Axis formatting ----
ax.set_xlim(DOMAIN_MIN - 0.5, DOMAIN_MAX + 0.5)
ax.xaxis.set_major_locator(ticker.MultipleLocator(10))
ax.xaxis.set_minor_locator(ticker.MultipleLocator(5))
ax.yaxis.set_major_locator(ticker.MultipleLocator(1))
ax.tick_params(axis='both', which='major', labelsize=11, direction='in', length=5)
ax.tick_params(axis='both', which='minor', direction='in', length=3)
ax.grid(True, which='major', linestyle=':', linewidth=0.6, alpha=0.4, color='gray')
ax.spines[['top', 'right']].set_visible(False)
ax.spines[['left', 'bottom']].set_linewidth(1.2)

# ---- Set y limits from data ----
all_vals = np.concatenate([v for v in [dsas_p1, dsas_p2, cs_p1, cs_p2]
                           if not np.all(np.isnan(v))])
ymin_data = np.nanmin(all_vals)
ymax_data = np.nanmax(all_vals)
pad = (ymax_data - ymin_data) * 0.06
ax.set_ylim(ymin_data - pad, ymax_data + pad)
ybot, ytop = ax.get_ylim()

# ---- Geographic labels ----
xmin_d, xmax_d = ax.get_xlim()
xrange_d = xmax_d - xmin_d

for domain, place in GEO_LABELS:
    x_frac = (domain - xmin_d) / xrange_d
    ax.axvline(domain, color='#444444', linewidth=1.0, linestyle='--',
               alpha=0.45, zorder=2)
    y_frac = LABEL_Y.get(place, 0.95)
    ax.text(x_frac, y_frac, place,
            transform=ax.transAxes,
            fontsize=10, fontweight='bold', color='#222222',
            ha='center', va='top', clip_on=False,
            bbox=dict(boxstyle='round,pad=0.25', facecolor='white',
                      edgecolor='#999999', linewidth=0.8, alpha=0.92))

# ---- Erosion / Accretion labels ----
zero_frac = (0 - ybot) / (ytop - ybot)
ax.text(1.005, zero_frac + (1 - zero_frac) / 2, 'Accretion ▲',
        transform=ax.transAxes, fontsize=9.5, color='#555555',
        ha='left', va='center', style='italic')
ax.text(1.005, zero_frac / 2, 'Erosion ▼',
        transform=ax.transAxes, fontsize=9.5, color='#555555',
        ha='left', va='center', style='italic')

# ---- Axis labels ----
ax.set_xlabel('CASCADE Model Domain (500 m alongshore)',
              fontsize=12, fontweight='bold', labelpad=8)
ax.set_ylabel('Shoreline Change Rate (m/yr)',
              fontsize=12, fontweight='bold', labelpad=8)

# ---- Geographic orientation ----
ax.text(0.0, 1.015, '← S  |  Cape Hatteras',
        transform=ax.transAxes, fontsize=9, color='#444444',
        ha='left', va='bottom', style='italic', clip_on=False)
ax.text(1.0, 1.015, 'Pea Island  |  N →',
        transform=ax.transAxes, fontsize=9, color='#444444',
        ha='right', va='bottom', style='italic', clip_on=False)

# ---- Title ----
ax.set_title('DSAS vs. CoastSat Shoreline Change Rates — Hatteras Island, NC',
             fontsize=13.5, fontweight='bold', pad=12, color='#1a2a3a')

# ---- Legend — two columns (DSAS | CoastSat) ----
legend_elements = []
if show_dsas_p1:
    legend_elements.append(
        Line2D([0], [0], color=COLOR_DSAS_P1, lw=2.8, ls='-',
               label='DSAS  1978–1997  (calibration)'))
if show_dsas_p2:
    legend_elements.append(
        Line2D([0], [0], color=COLOR_DSAS_P2, lw=2.8, ls='-',
               label='DSAS  1997–2019  (test)'))
if show_cs_p1:
    legend_elements.append(
        Line2D([0], [0], color=COLOR_CS_P1, lw=2.0, ls='--',
               label='CoastSat  1978–1997  (calibration)'))
if show_cs_p2:
    legend_elements.append(
        Line2D([0], [0], color=COLOR_CS_P2, lw=2.0, ls='--',
               label='CoastSat  1997–2019  (test)'))

ax.legend(handles=legend_elements,
          loc='lower center', bbox_to_anchor=(0.5, 0.02),
          ncol=2, fontsize=10,
          framealpha=0.95, edgecolor='#cccccc', frameon=True)

# ---- Caption ----
std_note = ' Shaded bands: ±1 std across transects per domain.' if SHOW_STD_BANDS else ''
caption = (
    'Rates: LRR averaged across transects per 500-m domain.' + std_note + ' '
    'DSAS data: Hapke & Henderson (2015, USGS OFR 2015-1002); Kratzmann et al. (2017, doi:10.5066/F74X55X7). '
    'CoastSat: satellite-derived shorelines (Vos et al. 2019). '
)
fig.text(0.018, 0.005, caption,
         fontsize=7.5, color='#666666', ha='left', va='bottom', style='italic')

plt.tight_layout(rect=[0, 0.05, 1, 1])
plt.savefig(OUTPUT_FILE, dpi=300, bbox_inches='tight', facecolor='white')
print(f"\nSaved: {OUTPUT_FILE}")
plt.close()
