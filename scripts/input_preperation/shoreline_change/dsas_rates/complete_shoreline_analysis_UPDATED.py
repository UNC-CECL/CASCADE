"""
Complete Shoreline Change Rate Analysis - All-in-One
Analyzes DSAS data and creates publication-quality visualizations
Run this once and get everything!
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# ============================================================================
# PART 1: LOAD AND CALCULATE RATES
# ============================================================================

print("="*70)
print("SHORELINE CHANGE RATE ANALYSIS - HATTERAS ISLAND")
print("="*70)
print("\nPART 1: Loading data and calculating rates...")
print("-"*70)

# Load the intersection data
df = pd.read_csv('All_Shoreline_Transect_Intersections.csv')

# Define column names
TRANSECT_ID_COL = 'Transects_100m_LineID'
DOMAIN_ID_COL = 'Transects_100m_AddSpatialJoin_domain_id'
YEAR_COL = 'Year'
DISTANCE_COL = 'NEAR_DIST'

print(f"Total intersection points: {len(df)}")
print(f"Available years: {sorted(df[YEAR_COL].unique())}")
print(f"Unique transects: {df[TRANSECT_ID_COL].nunique()}")
print(f"Unique domains: {df[DOMAIN_ID_COL].nunique()}")

# Clean data
df_clean = df[[TRANSECT_ID_COL, DOMAIN_ID_COL, YEAR_COL, DISTANCE_COL]].dropna()

# Check for duplicate intersections
duplicates = df_clean.groupby([TRANSECT_ID_COL, YEAR_COL]).size()
duplicates = duplicates[duplicates > 1]
if len(duplicates) > 0:
    print(f"⚠ Found {len(duplicates)} transect-year combos with multiple intersections")
    print(f"  (Shorelines crossing transects twice - will average)")

# Handle duplicates by averaging
df_clean_avg = df_clean.groupby([TRANSECT_ID_COL, YEAR_COL], as_index=False)[DISTANCE_COL].mean()

# Pivot to get distances by year for each transect
transect_distances = df_clean_avg.pivot(
    index=TRANSECT_ID_COL, 
    columns=YEAR_COL, 
    values=DISTANCE_COL
)

years = sorted(transect_distances.columns.tolist())
print(f"\nCalculating rates for: {years}")

# Calculate rates for each consecutive time period
transect_rates = pd.DataFrame(index=transect_distances.index)

for i in range(len(years) - 1):
    year1, year2 = years[i], years[i + 1]
    period_name = f'EPR_{year1}_{year2}'
    dt = year2 - year1
    
    # NOTE: -1 multiplier for correct sign convention
    # Negative = erosion, Positive = accretion
    transect_rates[period_name] = -1 * (transect_distances[year2] - transect_distances[year1]) / dt

# Calculate full calibration periods
if 1978 in years and 1997 in years:
    transect_rates['EPR_1978_1997'] = -1 * (transect_distances[1997] - transect_distances[1978]) / 19

if 1997 in years and 2019 in years:
    transect_rates['EPR_1997_2019'] = -1 * (transect_distances[2019] - transect_distances[1997]) / 22

if len(years) >= 2:
    full_period_name = f'EPR_{years[0]}_{years[-1]}'
    transect_rates[full_period_name] = -1 * (transect_distances[years[-1]] - transect_distances[years[0]]) / (years[-1] - years[0])

# Aggregate by domain
transect_domain_map = df_clean[[TRANSECT_ID_COL, DOMAIN_ID_COL]].drop_duplicates()
transect_rates_with_domain = transect_rates.merge(
    transect_domain_map,
    left_index=True,
    right_on=TRANSECT_ID_COL,
    how='left'
)

rate_cols = [col for col in transect_rates.columns if 'EPR_' in col]
domain_rates = transect_rates_with_domain.groupby(DOMAIN_ID_COL)[rate_cols].mean()
domain_rates['Domain'] = domain_rates.index

# Save results
domain_rates.to_csv('domain_shoreline_change_rates.csv', index=False)
print(f"\n✓ Saved: domain_shoreline_change_rates.csv")

# ============================================================================
# PART 2: CREATE PUBLICATION-QUALITY VISUALIZATIONS
# ============================================================================

print("\n" + "="*70)
print("PART 2: Creating publication-quality visualizations...")
print("-"*70)

period_names = {
    'EPR_1978_1987': '1978-1987',
    'EPR_1987_1997': '1987-1997', 
    'EPR_1997_2009': '1997-2009',
    'EPR_2009_2019': '2009-2019'
}

# ============================================================================
# OPTION 1: Sequential color scheme
# ============================================================================

fig, ax = plt.subplots(figsize=(18, 8))

colors_seq = ['#08519c', '#3182bd', '#6baed6', '#c6dbef']
periods = ['EPR_1978_1987', 'EPR_1987_1997', 'EPR_1997_2009', 'EPR_2009_2019']

for i, (col, color) in enumerate(zip(periods, colors_seq)):
    if col in domain_rates.columns:
        label = period_names.get(col, col)
        ax.plot(domain_rates['Domain'], domain_rates[col], 
                label=label, linewidth=2.5, color=color, alpha=0.9, zorder=10-i)

ax.axvspan(70, 90, alpha=0.12, color='#d62728', label='Rodanthe', zorder=0)
ax.axhline(y=0, color='black', linestyle='--', linewidth=1.5, alpha=0.7, zorder=1)

ax.set_xlabel('Domain Number', fontsize=14, fontweight='bold')
ax.set_ylabel('Shoreline Change Rate (m/yr)\n← Erosion | Accretion →', fontsize=13, fontweight='bold')
ax.set_title('Shoreline Change Rates: Temporal Evolution', fontsize=16, fontweight='bold', pad=20)
ax.legend(loc='upper right', fontsize=11, framealpha=0.95, edgecolor='gray')
ax.grid(True, alpha=0.25, linestyle=':', linewidth=0.8)
ax.set_xlim(domain_rates['Domain'].min() - 1, domain_rates['Domain'].max() + 1)

plt.tight_layout()
plt.savefig('v1_sequential_colors.png', dpi=300, bbox_inches='tight')
print("✓ Saved: v1_sequential_colors.png")
plt.close()

# ============================================================================
# OPTION 2: Faceted panels
# ============================================================================

fig, axes = plt.subplots(4, 1, figsize=(16, 12), sharex=True)

colors_distinct = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728']
period_labels = ['1978-1987', '1987-1997', '1997-2009', '2009-2019']

for ax, col, color, label in zip(axes, periods, colors_distinct, period_labels):
    if col in domain_rates.columns:
        ax.plot(domain_rates['Domain'], domain_rates[col], linewidth=2.5, color=color)
        ax.axvspan(77, 83, alpha=0.12, color='red', zorder=0)
        ax.axhline(y=0, color='black', linestyle='--', linewidth=1, alpha=0.6)
        ax.text(0.02, 0.95, label, transform=ax.transAxes, 
                fontsize=12, fontweight='bold', va='top',
                bbox=dict(boxstyle='round', facecolor='white', alpha=0.8, edgecolor=color, linewidth=2))
        ax.set_ylabel('Rate (m/yr)', fontsize=11, fontweight='bold')
        ax.grid(True, alpha=0.25, linestyle=':', linewidth=0.8)
        ax.set_ylim(domain_rates[col].min() - 1, domain_rates[col].max() + 1)

axes[-1].set_xlabel('Domain Number', fontsize=14, fontweight='bold')
axes[-1].set_xticks(range(0, int(domain_rates['Domain'].max()) + 10, 10))
axes[-1].set_xticks(range(0, int(domain_rates['Domain'].max()) + 10, 5), minor=True)
fig.suptitle('Shoreline Change Rates by Time Period', fontsize=16, fontweight='bold', y=0.995)
plt.tight_layout()
plt.savefig('v2_faceted.png', dpi=300, bbox_inches='tight')
print("✓ Saved: v2_faceted.png")
plt.close()

# ============================================================================
# OPTION 3: Early vs Recent comparison
# ============================================================================

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(18, 6), sharey=True)

if 'EPR_1978_1987' in domain_rates.columns:
    ax1.plot(domain_rates['Domain'], domain_rates['EPR_1978_1987'], 
            label='1978-1987', linewidth=3, color='#2166ac', marker='o', markersize=3)
if 'EPR_1987_1997' in domain_rates.columns:
    ax1.plot(domain_rates['Domain'], domain_rates['EPR_1987_1997'], 
            label='1987-1997', linewidth=3, color='#67a9cf', marker='s', markersize=3)

ax1.axvspan(77, 83, alpha=0.12, color='red', zorder=0)
ax1.axhline(y=0, color='black', linestyle='--', linewidth=1.5, alpha=0.7)
ax1.set_xlabel('Domain Number', fontsize=13, fontweight='bold')
ax1.set_ylabel('Shoreline Change Rate (m/yr)\n← Erosion | Accretion →', fontsize=12, fontweight='bold')
ax1.set_title('Earlier Period (1978-1997)', fontsize=14, fontweight='bold', pad=15)
ax1.legend(loc='best', fontsize=11, framealpha=0.95)
ax1.grid(True, alpha=0.25, linestyle=':', linewidth=0.8)
ax1.set_xticks(range(0, int(domain_rates['Domain'].max()) + 10, 10))
ax1.set_xticks(range(0, int(domain_rates['Domain'].max()) + 10, 5), minor=True)

if 'EPR_1997_2009' in domain_rates.columns:
    ax2.plot(domain_rates['Domain'], domain_rates['EPR_1997_2009'], 
            label='1997-2009', linewidth=3, color='#ef8a62', marker='^', markersize=3)
if 'EPR_2009_2019' in domain_rates.columns:
    ax2.plot(domain_rates['Domain'], domain_rates['EPR_2009_2019'], 
            label='2009-2019', linewidth=3, color='#b2182b', marker='d', markersize=3)

ax2.axvspan(77, 83, alpha=0.12, color='red', label='Rodanthe', zorder=0)
ax2.axhline(y=0, color='black', linestyle='--', linewidth=1.5, alpha=0.7)
ax2.set_xlabel('Domain Number', fontsize=13, fontweight='bold')
ax2.set_title('Recent Period (1997-2019)', fontsize=14, fontweight='bold', pad=15)
ax2.legend(loc='best', fontsize=11, framealpha=0.95)
ax2.grid(True, alpha=0.25, linestyle=':', linewidth=0.8)
ax2.set_xticks(range(0, int(domain_rates['Domain'].max()) + 10, 10))
ax2.set_xticks(range(0, int(domain_rates['Domain'].max()) + 10, 5), minor=True)

plt.tight_layout()
plt.savefig('v3_early_vs_recent.png', dpi=300, bbox_inches='tight')
print("✓ Saved: v3_early_vs_recent.png")
plt.close()

# ============================================================================
# OPTION 4: Calibration periods only (cleanest for CASCADE comparison)
# ============================================================================

if 'EPR_1978_1997' in domain_rates.columns and 'EPR_1997_2019' in domain_rates.columns:
    fig, ax = plt.subplots(figsize=(18, 8))
    
    ax.plot(domain_rates['Domain'], domain_rates['EPR_1978_1997'], 
            label='1978-1997 (Calibration Period 1)', linewidth=4, 
            color='#2166ac', marker='o', markersize=4, alpha=0.9)
    
    ax.plot(domain_rates['Domain'], domain_rates['EPR_1997_2019'], 
            label='1997-2019 (Calibration Period 2)', linewidth=4, 
            color='#b2182b', marker='s', markersize=4, alpha=0.9)
    
    ax.axvspan(77, 83, alpha=0.15, color='orange', label='Rodanthe Area', zorder=0)
    ax.axhline(y=0, color='black', linestyle='--', linewidth=1.5, alpha=0.7)
    
    ax.set_xlabel('Domain Number', fontsize=14, fontweight='bold')
    ax.set_ylabel('Shoreline Change Rate (m/yr)\n← Erosion | Accretion →', fontsize=13, fontweight='bold')
    ax.set_title('CASCADE Calibration Periods: Observed Shoreline Change Rates', 
                 fontsize=16, fontweight='bold', pad=20)
    ax.legend(loc='upper right', fontsize=13, framealpha=0.95, edgecolor='gray')
    ax.grid(True, alpha=0.25, linestyle=':', linewidth=0.8)
    ax.set_xlim(domain_rates['Domain'].min() - 1, domain_rates['Domain'].max() + 1)
    ax.set_xticks(range(0, int(domain_rates['Domain'].max()) + 10, 10))
    ax.set_xticks(range(0, int(domain_rates['Domain'].max()) + 10, 5), minor=True)
    
    plt.tight_layout()
    plt.savefig('v4_calibration_periods.png', dpi=300, bbox_inches='tight')
    print("✓ Saved: v4_calibration_periods.png")
    plt.close()

# ============================================================================
# SUMMARY STATISTICS
# ============================================================================

print("\n" + "="*70)
print("SUMMARY STATISTICS")
print("="*70)

print("\nRODANTHE AREA (Domains 77-83):")
rodanthe_summary = domain_rates[(domain_rates['Domain'] >= 77) & (domain_rates['Domain'] <= 83)]
for col in rate_cols:
    if col in rodanthe_summary.columns:
        mean_rate = rodanthe_summary[col].mean()
        min_rate = rodanthe_summary[col].min()
        max_rate = rodanthe_summary[col].max()
        trend = "EROSIONAL ⚠" if mean_rate < -0.5 else ("Erosional" if mean_rate < 0 else "Accretional")
        print(f"  {col:20s}: {mean_rate:+6.2f} m/yr ({trend})")
        print(f"  {'':20s}  Range: {min_rate:+6.2f} to {max_rate:+6.2f} m/yr")

print("\nFULL ISLAND AVERAGES:")
for col in rate_cols:
    if col in domain_rates.columns:
        mean_rate = domain_rates[col].mean()
        print(f"  {col:20s}: {mean_rate:+6.2f} m/yr")

print("\n" + "="*70)
print("ANALYSIS COMPLETE!")
print("="*70)
print("\nGenerated files:")
print("  • domain_shoreline_change_rates.csv (data table)")
print("  • v1_sequential_colors.png (all periods, color gradient)")
print("  • v2_faceted.png (separate panels)")
print("  • v3_early_vs_recent.png (two-panel comparison)")
print("  • v4_calibration_periods.png (CASCADE calibration periods)")
print("\nAll plots use correct sign convention:")
print("  Negative = Erosion (landward) | Positive = Accretion (seaward)")
print("="*70)
