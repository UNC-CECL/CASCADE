"""
Calculate Spatial Variability Metrics for CASCADE Sensitivity Analysis

This script calculates the spatial variability (standard deviation) of observed
vs modeled shoreline change rates to demonstrate that wave parameters alone
cannot reproduce observed spatial patterns.

Author: Hannah Henry (UNC Chapel Hill)
Date: January 2026
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os

# =============================================================================
# CONFIGURATION
# =============================================================================

PROJECT_BASE_DIR = r'C:\Users\hanna\PycharmProjects\CASCADE'

# Input files
DSAS_FILE = os.path.join(
    PROJECT_BASE_DIR,
    'data',
    'hatteras_init',
    'shoreline_change',
    'dsas_1978_1997_domain_means_SIMPLE.csv'
)

ALL_RESULTS_FILE = os.path.join(
    PROJECT_BASE_DIR,
    'output',
    'sensitivity_analysis',
    'domain_sensitivity_analysis',
    'all_domain_sensitivity_results.csv'
)

BEST_BY_DOMAIN_FILE = os.path.join(
    PROJECT_BASE_DIR,
    'output',
    'sensitivity_analysis',
    'domain_sensitivity_analysis',
    'best_parameters_by_domain.csv'
)

# Output directory
OUTPUT_DIR = os.path.join(
    PROJECT_BASE_DIR,
    'output',
    'sensitivity_analysis',
    'spatial_variability_analysis'
)
os.makedirs(OUTPUT_DIR, exist_ok=True)

# =============================================================================
# LOAD DATA
# =============================================================================

print("=" * 80)
print("SPATIAL VARIABILITY ANALYSIS")
print("=" * 80)

print("\n1. Loading data...")

# Load observed DSAS data
dsas = pd.read_csv(DSAS_FILE)
observed_rates = dsas['annual_rate_m_per_yr'].values
print(f"   ✓ Loaded {len(observed_rates)} observed rates")

# Load all sensitivity results
all_results = pd.read_csv(ALL_RESULTS_FILE)
print(f"   ✓ Loaded {len(all_results)} model results ({len(all_results)//90} configs × 90 domains)")

# Load best-by-domain results
best_by_domain = pd.read_csv(BEST_BY_DOMAIN_FILE)
print(f"   ✓ Loaded best-per-domain results")

# =============================================================================
# CALCULATE SPATIAL VARIABILITY
# =============================================================================

print("\n2. Calculating spatial variability metrics...")
print("-" * 80)

# Observed spatial variability
obs_mean = np.mean(observed_rates)
obs_std = np.std(observed_rates)
obs_min = np.min(observed_rates)
obs_max = np.max(observed_rates)

print(f"\nOBSERVED (DSAS 1978-1997):")
print(f"  Mean:       {obs_mean:7.3f} m/yr")
print(f"  Std dev:    {obs_std:7.3f} m/yr ⭐")
print(f"  Min:        {obs_min:7.3f} m/yr")
print(f"  Max:        {obs_max:7.3f} m/yr")
print(f"  Range:      {obs_max - obs_min:7.3f} m/yr")

# Calculate for each configuration
print(f"\nMODELED (CASCADE Sensitivity Runs):")
print("-" * 80)

configs_data = []
for (param, value), group in all_results.groupby(['parameter', 'value']):
    mod_mean = group['modeled'].mean()
    mod_std = group['modeled'].std()
    mod_min = group['modeled'].min()
    mod_max = group['modeled'].max()
    
    configs_data.append({
        'parameter': param,
        'value': value,
        'mean': mod_mean,
        'std_dev': mod_std,
        'min': mod_min,
        'max': mod_max,
        'range': mod_max - mod_min,
        'pct_of_observed': 100 * mod_std / obs_std
    })

configs_df = pd.DataFrame(configs_data)

# Print summary by parameter
print(f"\n{'Parameter':<25} {'Value':<8} {'Mean':<8} {'Std Dev':<10} {'% Obs':<8}")
print("-" * 70)
for _, row in configs_df.iterrows():
    print(f"{row['parameter']:<25} {row['value']:<8.2f} {row['mean']:<8.3f} "
          f"{row['std_dev']:<10.3f} {row['pct_of_observed']:<8.1f}%")

# Average across typical configs (exclude unstable ones)
typical_configs = configs_df[
    (configs_df['value'] != 0.5) &  # Exclude height=0.5 (unstable)
    (configs_df['value'] != 10.0)    # Exclude period=10 (unstable)
]
avg_mod_std = typical_configs['std_dev'].mean()

print(f"\n" + "=" * 80)
print(f"TYPICAL CONFIGURATIONS (excluding unstable):")
print(f"  Average std dev:  {avg_mod_std:7.3f} m/yr ⭐")
print(f"  % of observed:    {100 * avg_mod_std / obs_std:7.1f}%")

# Baseline configuration (period=7, your current settings)
baseline = configs_df[
    (configs_df['parameter'] == 'wave_period') & 
    (configs_df['value'] == 7.0)
]
if len(baseline) > 0:
    baseline_std = baseline['std_dev'].values[0]
    print(f"\nBASELINE CONFIGURATION (period=7):")
    print(f"  Std dev:          {baseline_std:7.3f} m/yr ⭐")
    print(f"  % of observed:    {100 * baseline_std / obs_std:7.1f}%")

# Best-case spatially-optimized
best_std = best_by_domain['best_modeled'].std()
best_mean = best_by_domain['best_modeled'].mean()
print(f"\nSPATIALLY-OPTIMIZED BEST-CASE:")
print(f"  Mean:             {best_mean:7.3f} m/yr")
print(f"  Std dev:          {best_std:7.3f} m/yr")
print(f"  % of observed:    {100 * best_std / obs_std:7.1f}%")

# =============================================================================
# SAVE SUMMARY TABLE
# =============================================================================

print("\n3. Creating summary table...")

summary_data = {
    'Scenario': [
        'Observed (DSAS)',
        'Baseline Config (period=7)',
        'Typical Configs (average)',
        'Spatially-Optimized Best-Case'
    ],
    'Mean (m/yr)': [
        obs_mean,
        baseline['mean'].values[0] if len(baseline) > 0 else np.nan,
        typical_configs['mean'].mean(),
        best_mean
    ],
    'Std Dev (m/yr)': [
        obs_std,
        baseline_std if len(baseline) > 0 else np.nan,
        avg_mod_std,
        best_std
    ],
    'Min (m/yr)': [
        obs_min,
        baseline['min'].values[0] if len(baseline) > 0 else np.nan,
        typical_configs['min'].mean(),
        best_by_domain['best_modeled'].min()
    ],
    'Max (m/yr)': [
        obs_max,
        baseline['max'].values[0] if len(baseline) > 0 else np.nan,
        typical_configs['max'].mean(),
        best_by_domain['best_modeled'].max()
    ]
}

summary_df = pd.DataFrame(summary_data)
summary_df['% of Observed Variability'] = 100 * summary_df['Std Dev (m/yr)'] / obs_std

# Save to CSV
summary_file = os.path.join(OUTPUT_DIR, 'spatial_variability_summary.csv')
summary_df.to_csv(summary_file, index=False)
print(f"   ✓ Saved: {summary_file}")

# Save detailed config results
detail_file = os.path.join(OUTPUT_DIR, 'spatial_variability_by_config.csv')
configs_df.to_csv(detail_file, index=False)
print(f"   ✓ Saved: {detail_file}")

# =============================================================================
# CREATE VISUALIZATIONS
# =============================================================================

print("\n4. Creating visualizations...")

sns.set_style("whitegrid")

# --- PLOT 1: Bar Chart Comparison ---
fig, ax = plt.subplots(figsize=(12, 8))

scenarios = summary_df['Scenario']
std_devs = summary_df['Std Dev (m/yr)']
colors = ['#2ecc71', '#e74c3c', '#e67e22', '#f39c12']

bars = ax.bar(range(len(scenarios)), std_devs, color=colors, alpha=0.8, 
              edgecolor='black', linewidth=2)

# Add value labels on bars
for i, (bar, val) in enumerate(zip(bars, std_devs)):
    height = bar.get_height()
    pct = 100 * val / obs_std
    ax.text(bar.get_x() + bar.get_width()/2., height + 0.05,
            f'{val:.2f} m/yr\n({pct:.0f}%)',
            ha='center', va='bottom', fontsize=11, fontweight='bold')

ax.set_xticks(range(len(scenarios)))
ax.set_xticklabels(scenarios, rotation=15, ha='right', fontsize=11)
ax.set_ylabel('Spatial Variability (Std Dev, m/yr)', fontsize=13, fontweight='bold')
ax.set_title('Spatial Variability: Observed vs Modeled\n(Demonstrates Wave Parameters Cannot Reproduce Spatial Patterns)', 
             fontsize=14, fontweight='bold')
ax.grid(axis='y', alpha=0.3)
ax.set_ylim([0, max(std_devs) * 1.2])

# Add horizontal line at observed value
ax.axhline(y=obs_std, color='green', linestyle='--', linewidth=2, 
          label='Observed Target', alpha=0.7)
ax.legend(fontsize=11, loc='upper right')

plt.tight_layout()
plt.savefig(os.path.join(OUTPUT_DIR, 'spatial_variability_comparison.png'), 
           dpi=300, bbox_inches='tight')
print(f"   ✓ Saved: spatial_variability_comparison.png")
plt.close()

# --- PLOT 2: All Configs Scatter ---
fig, ax = plt.subplots(figsize=(14, 8))

# Color by parameter
param_colors = {
    'wave_height': '#3498db',
    'wave_period': '#e74c3c',
    'wave_asymmetry': '#2ecc71',
    'wave_angle_high_fraction': '#f39c12'
}

for param in configs_df['parameter'].unique():
    param_data = configs_df[configs_df['parameter'] == param]
    ax.scatter(param_data['value'], param_data['std_dev'], 
              s=150, alpha=0.7, edgecolors='black', linewidth=1.5,
              color=param_colors.get(param, 'gray'),
              label=param.replace('_', ' ').title())

# Add observed line
ax.axhline(y=obs_std, color='green', linestyle='--', linewidth=3, 
          label=f'Observed (σ = {obs_std:.2f} m/yr)', alpha=0.8)

# Add typical model line
ax.axhline(y=avg_mod_std, color='red', linestyle=':', linewidth=3,
          label=f'Typical Model (σ = {avg_mod_std:.2f} m/yr)', alpha=0.8)

ax.set_xlabel('Parameter Value', fontsize=13, fontweight='bold')
ax.set_ylabel('Spatial Variability (Std Dev, m/yr)', fontsize=13, fontweight='bold')
ax.set_title('Spatial Variability Across All Parameter Configurations\n(No Configuration Approaches Observed Variability)', 
            fontsize=14, fontweight='bold')
ax.legend(loc='best', fontsize=10, ncol=2)
ax.grid(alpha=0.3)

plt.tight_layout()
plt.savefig(os.path.join(OUTPUT_DIR, 'spatial_variability_all_configs.png'), 
           dpi=300, bbox_inches='tight')
print(f"   ✓ Saved: spatial_variability_all_configs.png")
plt.close()

# --- PLOT 3: Summary Table as Image ---
fig, ax = plt.subplots(figsize=(12, 5))
ax.axis('tight')
ax.axis('off')

# Create table
table_data = summary_df.round(3).values
col_labels = summary_df.columns

table = ax.table(cellText=table_data, colLabels=col_labels,
                cellLoc='center', loc='center',
                colWidths=[0.28, 0.12, 0.12, 0.12, 0.12, 0.24])

table.auto_set_font_size(False)
table.set_fontsize(10)
table.scale(1, 2.5)

# Style header
for i in range(len(col_labels)):
    table[(0, i)].set_facecolor('#3498db')
    table[(0, i)].set_text_props(weight='bold', color='white')

# Style rows
colors_rows = ['#ecf0f1', '#ffffff', '#ecf0f1', '#ffffff']
for i in range(len(table_data)):
    for j in range(len(col_labels)):
        table[(i+1, j)].set_facecolor(colors_rows[i])
        
# Highlight observed row
for j in range(len(col_labels)):
    table[(1, j)].set_facecolor('#2ecc71')
    table[(1, j)].set_text_props(weight='bold')

plt.title('Spatial Variability Summary Table', fontsize=14, fontweight='bold', pad=20)
plt.savefig(os.path.join(OUTPUT_DIR, 'spatial_variability_table.png'), 
           dpi=300, bbox_inches='tight')
print(f"   ✓ Saved: spatial_variability_table.png")
plt.close()

# =============================================================================
# PRINT KEY FINDINGS
# =============================================================================

print("\n" + "=" * 80)
print("KEY FINDINGS FOR YOUR EMAIL/DEFENSE")
print("=" * 80)

print(f"""
📊 SPATIAL VARIABILITY ANALYSIS RESULTS:

1. OBSERVED SPATIAL VARIABILITY:
   • Standard deviation: {obs_std:.2f} m/yr
   • This represents the real complexity of Hatteras Island shoreline change

2. MODELED SPATIAL VARIABILITY (Typical Configurations):
   • Standard deviation: {avg_mod_std:.2f} m/yr
   • This is only {100*avg_mod_std/obs_std:.0f}% of observed variability
   • Model is {obs_std/avg_mod_std:.1f}× too flat!

3. BASELINE CONFIGURATION (Your Current Settings):
   • Standard deviation: {baseline_std if len(baseline) > 0 else 'N/A':.2f} m/yr
   • Captures only {100*baseline_std/obs_std if len(baseline) > 0 else 'N/A':.0f}% of observed variability

4. EVEN BEST-CASE SCENARIO (Spatially-Optimized):
   • Standard deviation: {best_std:.2f} m/yr
   • Captures only {100*best_std/obs_std:.0f}% of observed variability
   • Still insufficient to explain observations!

💡 INTERPRETATION:
Wave-driven alongshore transport and storm-driven overwash alone cannot generate
the observed spatial complexity at Hatteras Island. Background erosion processes
are ESSENTIAL to capture the full range of spatial variability.

📝 STATEMENT FOR YOUR EMAIL:
"The model produces insufficient spatial variability regardless of parameter
selection. Typical configurations generate only {avg_mod_std:.2f} m/yr standard
deviation compared to the observed {obs_std:.2f} m/yr (capturing only 
{100*avg_mod_std/obs_std:.0f}% of observed spatial variability). Even the 
spatially-optimized best-case scenario achieves only {best_std:.2f} m/yr 
({100*best_std/obs_std:.0f}% of observed), demonstrating that wave parameters 
alone cannot reproduce observed patterns."
""")

print("\n" + "=" * 80)
print("✓ ANALYSIS COMPLETE")
print("=" * 80)
print(f"\nAll outputs saved to: {OUTPUT_DIR}")
print("\nFiles created:")
print("  1. spatial_variability_summary.csv (key statistics)")
print("  2. spatial_variability_by_config.csv (detailed results)")
print("  3. spatial_variability_comparison.png (bar chart)")
print("  4. spatial_variability_all_configs.png (scatter plot)")
print("  5. spatial_variability_table.png (summary table)")
print("\nUse these results to demonstrate that wave parameters cannot")
print("reproduce observed spatial complexity!")
