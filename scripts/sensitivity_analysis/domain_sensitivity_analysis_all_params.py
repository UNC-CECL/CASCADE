"""
Domain-by-Domain Sensitivity Analysis - All Parameters

This script analyzes CASCADE performance at each domain for ALL sensitivity runs,
revealing which parameter combinations work best in which locations.

Author: Hannah Henry (UNC Chapel Hill)
Date: January 2026
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os
import glob

# =============================================================================
# CONFIGURATION
# =============================================================================

PROJECT_BASE_DIR = r'C:\Users\hanna\PycharmProjects\CASCADE'
SENSITIVITY_RESULTS_DIR = os.path.join(PROJECT_BASE_DIR, 'output', 'sensitivity_analysis')
OUTPUT_DIR = os.path.join(SENSITIVITY_RESULTS_DIR, 'domain_sensitivity_analysis')

# DSAS observed data
DSAS_FILE = os.path.join(
    PROJECT_BASE_DIR,
    'data',
    'hatteras_init',
    'shoreline_change',
    'dsas_1978_1997_domain_means_SIMPLE.csv'
)

# Path to sensitivity results CSV (has all run names)
SENSITIVITY_CSV = os.path.join(SENSITIVITY_RESULTS_DIR, 'sensitivity_results_20260127_134003.csv')

# Domain setup
NUM_REAL_DOMAINS = 90
NUM_BUFFER_DOMAINS = 15
START_REAL_INDEX = NUM_BUFFER_DOMAINS
END_REAL_INDEX = START_REAL_INDEX + NUM_REAL_DOMAINS
RUN_YEARS = 19

# =============================================================================
# HELPER FUNCTIONS
# =============================================================================

def get_x_s_TS(b3d):
    """Get shoreline time series."""
    if hasattr(b3d, "x_s_TS"):
        return np.array(b3d.x_s_TS)
    elif hasattr(b3d, "_x_s_TS"):
        return np.array(b3d._x_s_TS)
    else:
        raise AttributeError("No x_s_TS found")


def load_cascade_domain_rates(run_dir):
    """Extract CASCADE rates for each domain."""
    npz_files = glob.glob(os.path.join(run_dir, '*.npz'))
    if len(npz_files) == 0:
        return None
    
    try:
        data = np.load(npz_files[0], allow_pickle=True)
        cascade = data["cascade"][0]
        
        b3d_list = cascade.barrier3d
        ndom = len(b3d_list)
        nt = len(get_x_s_TS(b3d_list[0]))
        
        shoreline = np.zeros((nt, ndom))
        for j in range(ndom):
            xs = get_x_s_TS(b3d_list[j])
            shoreline[:, j] = xs
        
        shoreline = shoreline * 10.0  # dam → m
        
        initial = shoreline[0, START_REAL_INDEX:END_REAL_INDEX]
        final = shoreline[-1, START_REAL_INDEX:END_REAL_INDEX]
        rates = (final - initial) / RUN_YEARS
        
        return rates
    except Exception as e:
        print(f"  ⚠️  Error loading {run_dir}: {e}")
        return None


# =============================================================================
# MAIN ANALYSIS
# =============================================================================

print("=" * 80)
print("DOMAIN-BY-DOMAIN SENSITIVITY ANALYSIS - ALL PARAMETERS")
print("=" * 80)

os.makedirs(OUTPUT_DIR, exist_ok=True)

# Load DSAS observed data
print("\n1. Loading DSAS observations...")
dsas = pd.read_csv(DSAS_FILE)
observed_rates = dsas['annual_rate_m_per_yr'].values
print(f"   ✓ Loaded {len(observed_rates)} observations")

# Load sensitivity results summary
print("\n2. Loading sensitivity results...")
sens_results = pd.read_csv(SENSITIVITY_CSV)
print(f"   ✓ Found {len(sens_results)} sensitivity runs")

# Extract domain-level results for each run
print("\n3. Extracting domain-level results for all runs...")
all_domain_data = []

for idx, row in sens_results.iterrows():
    run_name = row['run_name']
    run_dir = row['output_dir']
    param = row['parameter']
    value = row['value']
    
    print(f"   Processing {run_name}...", end='')
    
    domain_rates = load_cascade_domain_rates(run_dir)
    
    if domain_rates is not None:
        # Calculate domain-level errors
        for domain_idx in range(NUM_REAL_DOMAINS):
            domain_num = domain_idx + 1
            obs = observed_rates[domain_idx]
            mod = domain_rates[domain_idx]
            error = mod - obs
            abs_error = abs(error)
            
            all_domain_data.append({
                'run_name': run_name,
                'parameter': param,
                'value': value,
                'domain': domain_num,
                'observed': obs,
                'modeled': mod,
                'error': error,
                'abs_error': abs_error,
                'background_needed': obs - mod
            })
        print(" ✓")
    else:
        print(" ✗ Failed")

# Create comprehensive dataframe
df_all = pd.DataFrame(all_domain_data)
print(f"\n✓ Extracted {len(df_all)} domain-level results ({len(df_all)//NUM_REAL_DOMAINS} runs × {NUM_REAL_DOMAINS} domains)")

# Save detailed results
detail_file = os.path.join(OUTPUT_DIR, 'all_domain_sensitivity_results.csv')
df_all.to_csv(detail_file, index=False)
print(f"✓ Saved detailed results: {detail_file}")

# =============================================================================
# ANALYSIS 1: WHICH PARAMETER VALUES WORK BEST AT EACH DOMAIN?
# =============================================================================

print("\n4. Analyzing which parameters work best at each domain...")

# For each domain, find the best-performing parameter configuration
best_by_domain = []

for domain in range(1, NUM_REAL_DOMAINS + 1):
    domain_data = df_all[df_all['domain'] == domain]
    
    # Find configuration with minimum absolute error
    best_run = domain_data.loc[domain_data['abs_error'].idxmin()]
    
    best_by_domain.append({
        'domain': domain,
        'observed': best_run['observed'],
        'best_parameter': best_run['parameter'],
        'best_value': best_run['value'],
        'best_modeled': best_run['modeled'],
        'best_error': best_run['error'],
        'best_abs_error': best_run['abs_error'],
    })

df_best = pd.DataFrame(best_by_domain)
best_file = os.path.join(OUTPUT_DIR, 'best_parameters_by_domain.csv')
df_best.to_csv(best_file, index=False)
print(f"✓ Saved best parameters by domain: {best_file}")

# =============================================================================
# ANALYSIS 2: PARAMETER PERFORMANCE BY DOMAIN
# =============================================================================

print("\n5. Creating parameter performance heatmaps...")

# Create heatmaps for each parameter showing absolute error across domains
parameters = df_all['parameter'].unique()

for param in parameters:
    param_data = df_all[df_all['parameter'] == param]
    
    # Pivot to create heatmap: rows=parameter values, columns=domains
    pivot = param_data.pivot_table(
        values='abs_error',
        index='value',
        columns='domain',
        aggfunc='mean'
    )
    
    # Create heatmap
    fig, ax = plt.subplots(figsize=(18, 6))
    
    sns.heatmap(pivot, cmap='RdYlGn_r', vmin=0, vmax=3, 
                cbar_kws={'label': 'Absolute Error (m/yr)'},
                ax=ax, linewidths=0.1, linecolor='gray')
    
    ax.set_xlabel('Domain Number', fontsize=12, fontweight='bold')
    ax.set_ylabel(f'{param.replace("_", " ").title()} Value', fontsize=12, fontweight='bold')
    ax.set_title(f'Domain-Level Performance: {param.replace("_", " ").title()}\n(Green = Good, Red = Bad)', 
                 fontsize=14, fontweight='bold')
    
    plt.tight_layout()
    plt.savefig(os.path.join(OUTPUT_DIR, f'heatmap_{param}.png'), dpi=300, bbox_inches='tight')
    print(f"   ✓ Created heatmap for {param}")
    plt.close()

# =============================================================================
# ANALYSIS 3: SPATIAL PATTERNS IN BEST PARAMETERS
# =============================================================================

print("\n6. Analyzing spatial patterns in optimal parameters...")

fig, axes = plt.subplots(4, 1, figsize=(16, 12))

# Plot 1: Observed rates
ax1 = axes[0]
domains = np.arange(1, NUM_REAL_DOMAINS + 1)
ax1.plot(domains, observed_rates, 'o-', linewidth=2, markersize=5, color='darkblue', label='Observed')
ax1.axhline(0, color='gray', linestyle='--', alpha=0.5)
ax1.set_ylabel('Observed Rate\n(m/yr)', fontsize=11, fontweight='bold')
ax1.set_title('Observed Shoreline Change Rates (DSAS)', fontsize=13, fontweight='bold')
ax1.legend()
ax1.grid(alpha=0.3)

# Plot 2: Best absolute error by domain
ax2 = axes[1]
colors_err = ['green' if e < 0.5 else 'orange' if e < 1.0 else 'red' for e in df_best['best_abs_error']]
ax2.bar(domains, df_best['best_abs_error'], color=colors_err, alpha=0.7, edgecolor='black', linewidth=0.5)
ax2.set_ylabel('Best Achievable\nAbs Error (m/yr)', fontsize=11, fontweight='bold')
ax2.set_title('Minimum Absolute Error Across All Parameter Combinations', fontsize=13, fontweight='bold')
ax2.grid(alpha=0.3, axis='y')

# Plot 3: Which parameter matters most by domain
ax3 = axes[2]
param_colors = {
    'wave_height': 'blue',
    'wave_period': 'green',
    'wave_asymmetry': 'orange',
    'wave_angle_high_fraction': 'red'
}

for domain in range(1, NUM_REAL_DOMAINS + 1):
    best_param = df_best.loc[df_best['domain'] == domain, 'best_parameter'].values[0]
    color = param_colors.get(best_param, 'gray')
    ax3.bar(domain, 1, color=color, alpha=0.7, edgecolor='black', linewidth=0.5)

# Add legend
from matplotlib.patches import Patch
legend_elements = [Patch(facecolor=color, label=param.replace('_', ' ').title()) 
                  for param, color in param_colors.items()]
ax3.legend(handles=legend_elements, loc='upper right', fontsize=9)
ax3.set_ylabel('Best Parameter', fontsize=11, fontweight='bold')
ax3.set_title('Which Parameter Configuration Worked Best at Each Domain', fontsize=13, fontweight='bold')
ax3.set_yticks([])
ax3.grid(alpha=0.3, axis='x')

# Plot 4: Background erosion needed (using best config)
ax4 = axes[3]
colors_bg = ['darkred' if x < -2 else 'red' if x < 0 else 'lightblue' if x < 2 else 'darkblue' 
             for x in df_best['observed'] - df_best['best_modeled']]
ax4.bar(domains, df_best['observed'] - df_best['best_modeled'], 
        color=colors_bg, alpha=0.7, edgecolor='black', linewidth=0.5)
ax4.axhline(0, color='black', linestyle='-', linewidth=1.5)
ax4.set_ylabel('Background Erosion\nNeeded (m/yr)', fontsize=11, fontweight='bold')
ax4.set_xlabel('Domain Number', fontsize=12, fontweight='bold')
ax4.set_title('Required Background Erosion (Even with Optimal Wave Parameters)', fontsize=13, fontweight='bold')
ax4.grid(alpha=0.3, axis='y')

plt.tight_layout()
plt.savefig(os.path.join(OUTPUT_DIR, 'spatial_patterns_best_configs.png'), dpi=300, bbox_inches='tight')
print(f"✓ Created spatial patterns plot")
plt.close()

# =============================================================================
# ANALYSIS 4: PARAMETER IMPORTANCE BY DOMAIN
# =============================================================================

print("\n7. Analyzing parameter importance by domain...")

# For each domain, calculate sensitivity to each parameter (range of errors)
param_sensitivity = []

for domain in range(1, NUM_REAL_DOMAINS + 1):
    domain_data = df_all[df_all['domain'] == domain]
    obs = observed_rates[domain - 1]
    
    sensitivities = {}
    for param in parameters:
        param_data = domain_data[domain_data['parameter'] == param]
        if len(param_data) > 0:
            error_range = param_data['abs_error'].max() - param_data['abs_error'].min()
            sensitivities[param] = error_range
    
    param_sensitivity.append({
        'domain': domain,
        'observed': obs,
        **sensitivities
    })

df_sens = pd.DataFrame(param_sensitivity)

# Create stacked area plot showing which parameters matter where
fig, ax = plt.subplots(figsize=(16, 8))

domains = df_sens['domain'].values
bottom = np.zeros(len(domains))

for param in parameters:
    if param in df_sens.columns:
        values = df_sens[param].values
        ax.fill_between(domains, bottom, bottom + values, 
                        label=param.replace('_', ' ').title(),
                        alpha=0.7)
        bottom += values

ax.set_xlabel('Domain Number', fontsize=13, fontweight='bold')
ax.set_ylabel('Parameter Sensitivity (Error Range, m/yr)', fontsize=13, fontweight='bold')
ax.set_title('Parameter Importance by Domain\n(Larger area = More sensitive to that parameter)', 
             fontsize=15, fontweight='bold')
ax.legend(loc='upper left', fontsize=11)
ax.grid(alpha=0.3)

plt.tight_layout()
plt.savefig(os.path.join(OUTPUT_DIR, 'parameter_importance_by_domain.png'), dpi=300, bbox_inches='tight')
print(f"✓ Created parameter importance plot")
plt.close()

# =============================================================================
# ANALYSIS 5: COMPARISON MATRIX
# =============================================================================

print("\n8. Creating parameter comparison matrix...")

# For each parameter, show performance across all domains
fig, axes = plt.subplots(2, 2, figsize=(16, 12))
axes = axes.flatten()

for idx, param in enumerate(parameters):
    if idx >= 4:
        break
    
    ax = axes[idx]
    param_data = df_all[df_all['parameter'] == param]
    
    # Get unique values for this parameter
    param_values = sorted(param_data['value'].unique())
    
    # Plot each configuration
    for pval in param_values:
        val_data = param_data[param_data['value'] == pval].sort_values('domain')
        ax.plot(val_data['domain'], val_data['abs_error'], 
               'o-', label=f'{pval}', alpha=0.7, markersize=4, linewidth=1.5)
    
    ax.set_xlabel('Domain Number', fontsize=11, fontweight='bold')
    ax.set_ylabel('Absolute Error (m/yr)', fontsize=11, fontweight='bold')
    ax.set_title(f'{param.replace("_", " ").title()}', fontsize=12, fontweight='bold')
    ax.legend(loc='best', fontsize=8, ncol=2)
    ax.grid(alpha=0.3)
    ax.set_ylim([0, 6])

plt.suptitle('Domain-by-Domain Performance for Each Parameter', 
             fontsize=15, fontweight='bold')
plt.tight_layout()
plt.savefig(os.path.join(OUTPUT_DIR, 'parameter_comparison_matrix.png'), dpi=300, bbox_inches='tight')
print(f"✓ Created comparison matrix")
plt.close()

# =============================================================================
# STATISTICS
# =============================================================================

print("\n" + "=" * 80)
print("SUMMARY STATISTICS")
print("=" * 80)

print(f"\nOverall Performance:")
print(f"  Total domain-parameter combinations tested: {len(df_all)}")
print(f"  Average absolute error (all configs): {df_all['abs_error'].mean():.3f} m/yr")
print(f"  Best achievable MAE (optimal config per domain): {df_best['best_abs_error'].mean():.3f} m/yr")

print(f"\nBest Parameter Distribution:")
for param in parameters:
    count = (df_best['best_parameter'] == param).sum()
    pct = 100 * count / len(df_best)
    print(f"  {param:30s}: Best in {count:2d} domains ({pct:4.1f}%)")

print(f"\nRemaining Error (Even with Optimal Parameters):")
still_large_error = df_best[df_best['best_abs_error'] > 1.0]
print(f"  Domains with error > 1.0 m/yr: {len(still_large_error)} ({100*len(still_large_error)/len(df_best):.1f}%)")
print(f"  Domains with error > 2.0 m/yr: {(df_best['best_abs_error'] > 2.0).sum()}")
print(f"  Domains with error > 3.0 m/yr: {(df_best['best_abs_error'] > 3.0).sum()}")

print(f"\nDomains with Worst Remaining Error (Need Background Erosion Most):")
worst_domains = df_best.nlargest(10, 'best_abs_error')[['domain', 'observed', 'best_parameter', 
                                                          'best_value', 'best_abs_error']]
print(worst_domains.to_string(index=False))

print("\n" + "=" * 80)
print("✓ ANALYSIS COMPLETE")
print("=" * 80)
print(f"\nAll outputs saved to: {OUTPUT_DIR}")
print("\nKey Files:")
print("  - all_domain_sensitivity_results.csv (complete dataset)")
print("  - best_parameters_by_domain.csv (optimal config per domain)")
print("  - heatmap_[param].png (performance maps for each parameter)")
print("  - spatial_patterns_best_configs.png (4-panel spatial analysis)")
print("  - parameter_importance_by_domain.png (which params matter where)")
print("  - parameter_comparison_matrix.png (detailed comparison)")

print("\n🎯 KEY FINDING:")
print(f"   Even using OPTIMAL wave parameters for each domain,")
print(f"   {(df_best['best_abs_error'] > 1.0).sum()} domains ({100*(df_best['best_abs_error'] > 1.0).sum()/len(df_best):.1f}%) still have error > 1 m/yr")
print(f"   → Background erosion is ESSENTIAL, not optional!")
