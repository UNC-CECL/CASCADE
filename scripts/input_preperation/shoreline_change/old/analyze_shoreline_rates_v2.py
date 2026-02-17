"""
Shoreline Change Rate Analysis by Domain and Time Period
Analyzes DSAS intersection data to calculate erosion/accretion rates
across different time periods for Hatteras Island domains
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# ============================================================================
# STEP 1: Load and examine the data
# ============================================================================

print("="*70)
print("LOADING DATA")
print("="*70)

# Load the intersection data
df = pd.read_csv('../All_Shoreline_Transect_Intersections.csv')

# Define column names (these are from your actual file)
TRANSECT_ID_COL = 'Transects_100m_LineID'
DOMAIN_ID_COL = 'Transects_100m_AddSpatialJoin_domain_id'
YEAR_COL = 'Year'
DISTANCE_COL = 'NEAR_DIST'

print(f"\nTotal intersection points: {len(df)}")
print(f"Available years: {sorted(df[YEAR_COL].unique())}")
print(f"Number of unique transects: {df[TRANSECT_ID_COL].nunique()}")
print(f"Number of unique domains: {df[DOMAIN_ID_COL].nunique()}")

# Remove any rows with missing critical data
df_clean = df[[TRANSECT_ID_COL, DOMAIN_ID_COL, YEAR_COL, DISTANCE_COL]].dropna()
print(f"\nRows after removing missing data: {len(df_clean)}")

# Check transects per domain
transects_per_domain = df_clean.groupby(DOMAIN_ID_COL)[TRANSECT_ID_COL].nunique()
print(f"Average transects per domain: {transects_per_domain.mean():.1f}")
print(f"Min transects per domain: {transects_per_domain.min()}")
print(f"Max transects per domain: {transects_per_domain.max()}")

# Check for duplicate intersections (same transect-year combo)
duplicates = df_clean.groupby([TRANSECT_ID_COL, YEAR_COL]).size()
duplicates = duplicates[duplicates > 1]
if len(duplicates) > 0:
    print(f"\n⚠ Found {len(duplicates)} transect-year combinations with multiple intersections")
    print(f"  This can happen when shorelines curve back and cross a transect twice")
    print(f"  Will average the distances for these cases")

# ============================================================================
# STEP 2: Calculate rates at transect level
# ============================================================================

print("\n" + "="*70)
print("CALCULATING TRANSECT-LEVEL RATES")
print("="*70)

# Handle duplicates by averaging distances when a transect-year has multiple intersections
df_clean_avg = df_clean.groupby([TRANSECT_ID_COL, YEAR_COL], as_index=False)[DISTANCE_COL].mean()

# Pivot to get distances by year for each transect
transect_distances = df_clean_avg.pivot(
    index=TRANSECT_ID_COL, 
    columns=YEAR_COL, 
    values=DISTANCE_COL
)

print(f"\nYears with data: {transect_distances.columns.tolist()}")

# Get available years and sort them
years = sorted(transect_distances.columns.tolist())
print(f"Time periods to analyze: {years}")

# Calculate rates for each consecutive time period
transect_rates = pd.DataFrame(index=transect_distances.index)

for i in range(len(years) - 1):
    year1, year2 = years[i], years[i + 1]
    period_name = f'EPR_{year1}_{year2}'
    dt = year2 - year1
    
    # Calculate rate in m/yr
    transect_rates[period_name] = (transect_distances[year2] - transect_distances[year1]) / dt
    print(f"  Calculated {period_name}: {dt} years")

# Calculate full calibration periods if available
if 1978 in years and 1997 in years:
    transect_rates['EPR_1978_1997'] = (transect_distances[1997] - transect_distances[1978]) / 19
    print(f"  Calculated EPR_1978_1997: 19 years")

if 1997 in years and 2019 in years:
    transect_rates['EPR_1997_2019'] = (transect_distances[2019] - transect_distances[1997]) / 22
    print(f"  Calculated EPR_1997_2019: 22 years")

# Full time span
if len(years) >= 2:
    full_period_name = f'EPR_{years[0]}_{years[-1]}'
    transect_rates[full_period_name] = (transect_distances[years[-1]] - transect_distances[years[0]]) / (years[-1] - years[0])
    print(f"  Calculated {full_period_name}: {years[-1] - years[0]} years")

# ============================================================================
# STEP 3: Add domain information and aggregate by domain
# ============================================================================

print("\n" + "="*70)
print("AGGREGATING BY DOMAIN")
print("="*70)

# Create a transect-to-domain lookup from the original data (take unique combinations)
transect_domain_map = df_clean[[TRANSECT_ID_COL, DOMAIN_ID_COL]].drop_duplicates()

# Merge with transect rates
transect_rates_with_domain = transect_rates.merge(
    transect_domain_map,
    left_index=True,
    right_on=TRANSECT_ID_COL,
    how='left'
)

# Get only the rate columns
rate_cols = [col for col in transect_rates.columns if 'EPR_' in col]

# Calculate mean rate for each domain (averaging the 5 transects per domain)
domain_rates = transect_rates_with_domain.groupby(DOMAIN_ID_COL)[rate_cols].mean()
domain_rates['Domain'] = domain_rates.index

# Also calculate standard deviation for each domain (measure of variability)
domain_std = transect_rates_with_domain.groupby(DOMAIN_ID_COL)[rate_cols].std()

print(f"\nDomain-level statistics calculated")
print(f"Number of domains with data: {len(domain_rates)}")
print("\nFirst few domains:")
print(domain_rates.head(10))

# Save results
domain_rates.to_csv('domain_shoreline_change_rates.csv', index=False)
print("\n✓ Saved: domain_shoreline_change_rates.csv")

# ============================================================================
# STEP 4: Create plots
# ============================================================================

print("\n" + "="*70)
print("CREATING PLOTS")
print("="*70)

# Define consecutive periods for plotting
consecutive_periods = [(f'EPR_{years[i]}_{years[i+1]}', f'{years[i]}-{years[i+1]}') 
                       for i in range(len(years) - 1)]

# Color scheme
colors = plt.cm.tab10(np.arange(len(consecutive_periods)))

# ============================================================================
# PLOT 1: All domains comparison
# ============================================================================

fig, ax = plt.subplots(figsize=(16, 8))

for (col, label), color in zip(consecutive_periods, colors):
    if col in domain_rates.columns:
        ax.plot(domain_rates['Domain'], domain_rates[col], 
                label=label, linewidth=2.5, 
                marker='o', markersize=4, alpha=0.85, color=color)

# Highlight Rodanthe area (domains 70-90)
ax.axvspan(70, 90, alpha=0.15, color='red', label='Rodanthe', zorder=0)

# Add zero reference line
ax.axhline(y=0, color='black', linestyle='--', linewidth=1.5, alpha=0.7)

# Labels and formatting
ax.set_xlabel('Domain Number', fontsize=14, fontweight='bold')
ax.set_ylabel('Shoreline Change Rate (m/yr)', fontsize=14, fontweight='bold')
ax.set_title('Shoreline Change Rates Across Time Periods - Hatteras Island', 
             fontsize=16, fontweight='bold', pad=20)
ax.legend(loc='best', fontsize=12, framealpha=0.9)
ax.grid(True, alpha=0.3, linestyle=':', linewidth=0.8)
ax.set_xlim(domain_rates['Domain'].min() - 1, domain_rates['Domain'].max() + 1)

# Add horizontal line at y=0 for reference
ax.axhline(y=0, color='gray', linestyle='-', linewidth=0.8, alpha=0.5)

plt.tight_layout()
plt.savefig('shoreline_rates_all_domains.png', dpi=300, bbox_inches='tight')
print("✓ Saved: shoreline_rates_all_domains.png")
plt.close()

# ============================================================================
# PLOT 2: Rodanthe-focused analysis
# ============================================================================

fig, ax = plt.subplots(figsize=(12, 8))

# Filter to Rodanthe area (domains 65-90 to show context)
rodanthe_data = domain_rates[(domain_rates['Domain'] >= 65) & (domain_rates['Domain'] <= 90)]

for (col, label), color in zip(consecutive_periods, colors):
    if col in domain_rates.columns:
        ax.plot(rodanthe_data['Domain'], rodanthe_data[col], 
                label=label, linewidth=3, 
                marker='o', markersize=6, alpha=0.85, color=color)

# Highlight main Rodanthe zone (70-90)
ax.axvspan(70, 90, alpha=0.15, color='red', label='Rodanthe Core', zorder=0)

ax.axhline(y=0, color='black', linestyle='--', linewidth=1.5, alpha=0.7)

ax.set_xlabel('Domain Number', fontsize=14, fontweight='bold')
ax.set_ylabel('Shoreline Change Rate (m/yr)', fontsize=14, fontweight='bold')
ax.set_title('Rodanthe Area: Temporal Variation in Shoreline Change Rates', 
             fontsize=16, fontweight='bold', pad=20)
ax.legend(loc='best', fontsize=12, framealpha=0.9)
ax.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig('rodanthe_temporal_analysis.png', dpi=300, bbox_inches='tight')
print("✓ Saved: rodanthe_temporal_analysis.png")
plt.close()

# ============================================================================
# PLOT 3: Full calibration periods comparison (if available)
# ============================================================================

if 'EPR_1978_1997' in domain_rates.columns and 'EPR_1997_2019' in domain_rates.columns:
    fig, ax = plt.subplots(figsize=(16, 8))
    
    ax.plot(domain_rates['Domain'], domain_rates['EPR_1978_1997'], 
            label='1978-1997 (Calibration 1)', linewidth=3, 
            marker='o', markersize=5, alpha=0.85, color='#2E86AB')
    
    ax.plot(domain_rates['Domain'], domain_rates['EPR_1997_2019'], 
            label='1997-2019 (Calibration 2)', linewidth=3, 
            marker='s', markersize=5, alpha=0.85, color='#A23B72')
    
    # Highlight Rodanthe
    ax.axvspan(70, 90, alpha=0.15, color='red', label='Rodanthe', zorder=0)
    ax.axhline(y=0, color='black', linestyle='--', linewidth=1.5, alpha=0.7)
    
    ax.set_xlabel('Domain Number', fontsize=14, fontweight='bold')
    ax.set_ylabel('Shoreline Change Rate (m/yr)', fontsize=14, fontweight='bold')
    ax.set_title('CASCADE Calibration Periods: Observed Shoreline Change Rates', 
                 fontsize=16, fontweight='bold', pad=20)
    ax.legend(loc='best', fontsize=12, framealpha=0.9)
    ax.grid(True, alpha=0.3, linestyle=':', linewidth=0.8)
    ax.set_xlim(domain_rates['Domain'].min() - 1, domain_rates['Domain'].max() + 1)
    
    plt.tight_layout()
    plt.savefig('calibration_periods_comparison.png', dpi=300, bbox_inches='tight')
    print("✓ Saved: calibration_periods_comparison.png")
    plt.close()

# ============================================================================
# SUMMARY STATISTICS
# ============================================================================

print("\n" + "="*70)
print("SUMMARY STATISTICS")
print("="*70)

print("\nRodanthe Area (Domains 70-90) - Mean Rates:")
rodanthe_summary = domain_rates[(domain_rates['Domain'] >= 70) & (domain_rates['Domain'] <= 90)]
for col in rate_cols:
    if col in rodanthe_summary.columns:
        mean_rate = rodanthe_summary[col].mean()
        min_rate = rodanthe_summary[col].min()
        max_rate = rodanthe_summary[col].max()
        print(f"  {col}: {mean_rate:+.3f} m/yr (range: {min_rate:+.3f} to {max_rate:+.3f})")

print("\nFull Island - Mean Rates:")
for col in rate_cols:
    if col in domain_rates.columns:
        mean_rate = domain_rates[col].mean()
        print(f"  {col}: {mean_rate:+.3f} m/yr")

print("\n" + "="*70)
print("ANALYSIS COMPLETE!")
print("="*70)
print("\nGenerated files:")
print("  1. domain_shoreline_change_rates.csv")
print("  2. shoreline_rates_all_domains.png")
print("  3. rodanthe_temporal_analysis.png")
if 'EPR_1978_1997' in domain_rates.columns:
    print("  4. calibration_periods_comparison.png")

print("\nThese plots show how shoreline change rates vary across domains")
print("and time periods, with special focus on the Rodanthe area.")
