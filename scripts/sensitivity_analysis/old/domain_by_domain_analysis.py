"""
Domain-by-Domain CASCADE Performance Analysis

This script compares CASCADE model output to DSAS observations at each domain
to identify spatial patterns in model performance and inform background erosion
calculation.

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

PROJECT_BASE_DIR = r'/'

# Path to CASCADE output from your baseline run
CASCADE_OUTPUT_DIR = os.path.join(
    PROJECT_BASE_DIR,
    'output', 
    'sensitivity_analysis',
    'sens_wave_period_7_20260126_230543'  # Your baseline configuration
)

# Path to DSAS observed data
DSAS_FILE = os.path.join(
    PROJECT_BASE_DIR,
    'data',
    'hatteras_init',
    'shoreline_change',
    'dsas_1978_1997_domain_means_SIMPLE.csv'
)

# Output directory
OUTPUT_DIR = os.path.join(PROJECT_BASE_DIR, 'output', 'spatial_analysis')
os.makedirs(OUTPUT_DIR, exist_ok=True)

# Domain configuration
NUM_BUFFER_DOMAINS = 15
NUM_REAL_DOMAINS = 90
START_REAL_INDEX = NUM_BUFFER_DOMAINS
END_REAL_INDEX = START_REAL_INDEX + NUM_REAL_DOMAINS
RUN_YEARS = 19

# =============================================================================
# LOAD DATA
# =============================================================================

print("=" * 80)
print("DOMAIN-BY-DOMAIN MODEL PERFORMANCE ANALYSIS")
print("=" * 80)

# Load DSAS observed data
print("\n1. Loading DSAS observations...")
dsas = pd.read_csv(DSAS_FILE)
gis_domains = dsas['domain_id'].to_numpy()
observed_rates = dsas['annual_rate_m_per_yr'].to_numpy()
print(f"   ✓ Loaded {len(observed_rates)} observations")

# Load CASCADE model output
print("\n2. Loading CASCADE model output...")

def get_x_s_TS(b3d):
    """Get shoreline time series."""
    if hasattr(b3d, "x_s_TS"):
        return np.array(b3d.x_s_TS)
    elif hasattr(b3d, "_x_s_TS"):
        return np.array(b3d._x_s_TS)
    else:
        raise AttributeError("No x_s_TS found")

# Find NPZ file
npz_files = glob.glob(os.path.join(CASCADE_OUTPUT_DIR, '*.npz'))
if len(npz_files) == 0:
    raise FileNotFoundError(f"No NPZ file found in {CASCADE_OUTPUT_DIR}")

data = np.load(npz_files[0], allow_pickle=True)
cascade = data["cascade"][0]

# Extract model rates
b3d_list = cascade.barrier3d
ndom = len(b3d_list)
nt = len(get_x_s_TS(b3d_list[0]))

shoreline = np.zeros((nt, ndom))
for j in range(ndom):
    xs = get_x_s_TS(b3d_list[j])
    shoreline[:, j] = xs

shoreline = shoreline * 10.0  # dam → m

# Calculate rates for real domains
initial = shoreline[0, START_REAL_INDEX:END_REAL_INDEX]
final = shoreline[-1, START_REAL_INDEX:END_REAL_INDEX]
modeled_rates = (final - initial) / RUN_YEARS

print(f"   ✓ Extracted CASCADE rates for {NUM_REAL_DOMAINS} domains")

# =============================================================================
# DOMAIN-BY-DOMAIN ANALYSIS
# =============================================================================

print("\n3. Calculating domain-by-domain metrics...")

# Create comprehensive dataframe
results = pd.DataFrame({
    'domain': np.arange(1, NUM_REAL_DOMAINS + 1),
    'observed': observed_rates,
    'modeled': modeled_rates,
})

# Calculate metrics
results['error'] = results['modeled'] - results['observed']
results['abs_error'] = np.abs(results['error'])
results['percent_error'] = 100 * results['error'] / (np.abs(results['observed']) + 0.01)  # avoid div by zero
results['background_needed'] = results['observed'] - results['modeled']

# Categorize domains
def categorize_domain(row):
    """Categorize based on observed and modeled behavior."""
    obs = row['observed']
    mod = row['modeled']
    err = row['error']
    
    # Check if signs match
    if np.sign(obs) == np.sign(mod):
        if abs(err) < 0.5:
            return 'Good Match'
        elif abs(err) < 1.0:
            return 'Moderate Match'
        else:
            return 'Wrong Magnitude'
    else:
        return 'Wrong Sign'

results['category'] = results.apply(categorize_domain, axis=1)

# Add spatial zones (approximate)
def get_zone(domain):
    """Approximate geographic zones of Hatteras Island."""
    if domain <= 30:
        return 'North (Rodanthe-Waves-Salvo)'
    elif domain <= 60:
        return 'Central (Avon-Buxton)'
    else:
        return 'South (Frisco-Hatteras)'

results['zone'] = results['domain'].apply(get_zone)

# Save detailed results
results.to_csv(os.path.join(OUTPUT_DIR, 'domain_by_domain_results.csv'), index=False)
print(f"   ✓ Saved detailed results to CSV")

# =============================================================================
# STATISTICS
# =============================================================================

print("\n4. Domain-by-Domain Statistics:")
print("=" * 80)

print(f"\nOverall Metrics:")
print(f"  Mean error (bias): {results['error'].mean():.3f} m/yr")
print(f"  Mean absolute error: {results['abs_error'].mean():.3f} m/yr")
print(f"  RMSE: {np.sqrt((results['error']**2).mean()):.3f} m/yr")
print(f"  Correlation: {np.corrcoef(results['observed'], results['modeled'])[0,1]:.3f}")

print(f"\nDomain Categories:")
for cat in ['Good Match', 'Moderate Match', 'Wrong Magnitude', 'Wrong Sign']:
    count = (results['category'] == cat).sum()
    pct = 100 * count / len(results)
    print(f"  {cat:20s}: {count:2d} domains ({pct:4.1f}%)")

print(f"\nBy Zone:")
for zone in results['zone'].unique():
    zone_data = results[results['zone'] == zone]
    print(f"\n  {zone}:")
    print(f"    Mean error: {zone_data['error'].mean():.3f} m/yr")
    print(f"    MAE: {zone_data['abs_error'].mean():.3f} m/yr")
    print(f"    Correlation: {np.corrcoef(zone_data['observed'], zone_data['modeled'])[0,1]:.3f}")

# =============================================================================
# VISUALIZATIONS
# =============================================================================

print("\n5. Creating gifs...")

sns.set_style("whitegrid")

# ============= PLOT 1: Observed vs Modeled Scatter =============
fig, ax = plt.subplots(figsize=(10, 10))

# Color by category
colors = {
    'Good Match': 'green',
    'Moderate Match': 'orange', 
    'Wrong Magnitude': 'red',
    'Wrong Sign': 'darkred'
}

for cat, color in colors.items():
    mask = results['category'] == cat
    if mask.sum() > 0:
        ax.scatter(results[mask]['observed'], results[mask]['modeled'],
                  c=color, label=cat, s=80, alpha=0.7, edgecolors='black', linewidth=0.5)

# 1:1 line
lim_min = min(results['observed'].min(), results['modeled'].min()) - 1
lim_max = max(results['observed'].max(), results['modeled'].max()) + 1
ax.plot([lim_min, lim_max], [lim_min, lim_max], 'k--', linewidth=2, label='1:1 Line', alpha=0.5)

# Quadrant lines
ax.axhline(0, color='gray', linestyle='-', linewidth=1, alpha=0.3)
ax.axvline(0, color='gray', linestyle='-', linewidth=1, alpha=0.3)

ax.set_xlabel('Observed Rate (m/yr)', fontsize=14, fontweight='bold')
ax.set_ylabel('Modeled Rate (m/yr)', fontsize=14, fontweight='bold')
ax.set_title('Domain-by-Domain: Observed vs Modeled Shoreline Change', 
             fontsize=16, fontweight='bold')
ax.legend(loc='upper left', fontsize=11)
ax.grid(True, alpha=0.3)
ax.set_aspect('equal')

plt.tight_layout()
plt.savefig(os.path.join(OUTPUT_DIR, 'observed_vs_modeled_scatter.png'), dpi=300, bbox_inches='tight')
print(f"   ✓ Saved: observed_vs_modeled_scatter.png")
plt.close()

# ============= PLOT 2: Spatial Patterns Along Island =============
fig, axes = plt.subplots(4, 1, figsize=(16, 12))

domains = results['domain']

# Panel 1: Observed vs Modeled
ax1 = axes[0]
ax1.plot(domains, results['observed'], 'o-', label='Observed (DSAS)', 
         linewidth=2.5, markersize=6, color='darkblue', alpha=0.8)
ax1.plot(domains, results['modeled'], 's-', label='Modeled (CASCADE)', 
         linewidth=2.5, markersize=5, color='red', alpha=0.7)
ax1.axhline(0, color='gray', linestyle='--', linewidth=1, alpha=0.5)
ax1.set_ylabel('Shoreline Change\nRate (m/yr)', fontsize=12, fontweight='bold')
ax1.set_title('Observed vs Modeled Shoreline Change Rates', fontsize=14, fontweight='bold')
ax1.legend(loc='best', fontsize=11)
ax1.grid(True, alpha=0.3)

# Panel 2: Error (Bias)
ax2 = axes[1]
colors_err = ['green' if abs(e) < 0.5 else 'orange' if abs(e) < 1.0 else 'red' 
              for e in results['error']]
ax2.bar(domains, results['error'], color=colors_err, alpha=0.7, edgecolor='black', linewidth=0.5)
ax2.axhline(0, color='black', linestyle='-', linewidth=2)
ax2.set_ylabel('Error\n(Model - Obs, m/yr)', fontsize=12, fontweight='bold')
ax2.set_title('Model Error by Domain (Positive = Overprediction, Negative = Underprediction)', 
              fontsize=14, fontweight='bold')
ax2.grid(True, alpha=0.3, axis='y')

# Panel 3: Absolute Error
ax3 = axes[2]
ax3.bar(domains, results['abs_error'], color='darkorange', alpha=0.7, edgecolor='black', linewidth=0.5)
ax3.axhline(results['abs_error'].mean(), color='red', linestyle='--', linewidth=2, 
           label=f'Mean MAE = {results["abs_error"].mean():.2f} m/yr')
ax3.set_ylabel('Absolute Error\n(m/yr)', fontsize=12, fontweight='bold')
ax3.set_title('Magnitude of Model Error (Regardless of Sign)', fontsize=14, fontweight='bold')
ax3.legend(loc='best', fontsize=11)
ax3.grid(True, alpha=0.3, axis='y')

# Panel 4: Background Erosion Needed
ax4 = axes[3]
colors_bg = ['darkred' if x < -2 else 'red' if x < 0 else 'lightblue' if x < 2 else 'darkblue' 
             for x in results['background_needed']]
ax4.bar(domains, results['background_needed'], color=colors_bg, alpha=0.7, 
        edgecolor='black', linewidth=0.5)
ax4.axhline(0, color='black', linestyle='-', linewidth=2)
ax4.set_ylabel('Background Erosion\nNeeded (m/yr)', fontsize=12, fontweight='bold')
ax4.set_xlabel('Domain Number', fontsize=13, fontweight='bold')
ax4.set_title('Required Background Erosion (Observed - Modeled)', fontsize=14, fontweight='bold')
ax4.grid(True, alpha=0.3, axis='y')

# Add zone labels
for ax in axes:
    ax.axvspan(0, 30.5, alpha=0.05, color='blue', label='North')
    ax.axvspan(30.5, 60.5, alpha=0.05, color='green', label='Central')
    ax.axvspan(60.5, 90.5, alpha=0.05, color='orange', label='South')

plt.tight_layout()
plt.savefig(os.path.join(OUTPUT_DIR, 'spatial_patterns_along_island.png'), dpi=300, bbox_inches='tight')
print(f"   ✓ Saved: spatial_patterns_along_island.png")
plt.close()

# ============= PLOT 3: Error Distribution =============
fig, axes = plt.subplots(2, 2, figsize=(14, 10))

# Histogram of errors
ax1 = axes[0, 0]
ax1.hist(results['error'], bins=30, color='steelblue', edgecolor='black', alpha=0.7)
ax1.axvline(0, color='red', linestyle='--', linewidth=2, label='Zero Error')
ax1.axvline(results['error'].mean(), color='darkred', linestyle='-', linewidth=2, 
           label=f'Mean = {results["error"].mean():.2f} m/yr')
ax1.set_xlabel('Error (Model - Observed, m/yr)', fontsize=12, fontweight='bold')
ax1.set_ylabel('Number of Domains', fontsize=12, fontweight='bold')
ax1.set_title('Distribution of Model Errors', fontsize=13, fontweight='bold')
ax1.legend()
ax1.grid(True, alpha=0.3)

# Histogram of absolute errors
ax2 = axes[0, 1]
ax2.hist(results['abs_error'], bins=30, color='darkorange', edgecolor='black', alpha=0.7)
ax2.axvline(results['abs_error'].mean(), color='darkred', linestyle='-', linewidth=2,
           label=f'Mean = {results["abs_error"].mean():.2f} m/yr')
ax2.set_xlabel('Absolute Error (m/yr)', fontsize=12, fontweight='bold')
ax2.set_ylabel('Number of Domains', fontsize=12, fontweight='bold')
ax2.set_title('Distribution of Absolute Errors', fontsize=13, fontweight='bold')
ax2.legend()
ax2.grid(True, alpha=0.3)

# Error by zone
ax3 = axes[1, 0]
zone_means = results.groupby('zone')['abs_error'].mean().sort_values()
colors = ['lightcoral', 'lightyellow', 'lightgreen']
ax3.barh(zone_means.index, zone_means.values, color=colors, edgecolor='black', linewidth=1.5)
ax3.set_xlabel('Mean Absolute Error (m/yr)', fontsize=12, fontweight='bold')
ax3.set_title('Model Performance by Geographic Zone', fontsize=13, fontweight='bold')
ax3.grid(True, alpha=0.3, axis='x')

# Category pie chart
ax4 = axes[1, 1]
category_counts = results['category'].value_counts()
colors_pie = [colors[cat] for cat in category_counts.index]
ax4.pie(category_counts.values, labels=category_counts.index, autopct='%1.1f%%',
       colors=colors_pie, startangle=90, textprops={'fontsize': 11, 'fontweight': 'bold'})
ax4.set_title('Domain Performance Categories', fontsize=13, fontweight='bold')

plt.tight_layout()
plt.savefig(os.path.join(OUTPUT_DIR, 'error_distribution_analysis.png'), dpi=300, bbox_inches='tight')
print(f"   ✓ Saved: error_distribution_analysis.png")
plt.close()

# ============= PLOT 4: Background Erosion Heatmap =============
fig, ax = plt.subplots(figsize=(16, 6))

# Create heatmap data
bg_needed = results['background_needed'].values.reshape(1, -1)
domains_arr = results['domain'].values.reshape(1, -1)

im = ax.imshow(bg_needed, cmap='RdBu_r', aspect='auto', vmin=-5, vmax=5,
              extent=[0.5, 90.5, 0, 1])

# Add colorbar
cbar = plt.colorbar(im, ax=ax, orientation='horizontal', pad=0.15, aspect=40)
cbar.set_label('Background Erosion Rate (m/yr)', fontsize=13, fontweight='bold')

# Customize
ax.set_xlim([0.5, 90.5])
ax.set_ylim([0, 1])
ax.set_xlabel('Domain Number', fontsize=13, fontweight='bold')
ax.set_yticks([])
ax.set_title('Spatial Pattern of Required Background Erosion\n(Red = Erosion Needed, Blue = Deposition Needed)', 
             fontsize=14, fontweight='bold')

# Add zone boundaries
ax.axvline(30.5, color='black', linestyle='--', linewidth=2, alpha=0.5)
ax.axvline(60.5, color='black', linestyle='--', linewidth=2, alpha=0.5)

# Add zone labels
ax.text(15, 0.5, 'North', ha='center', va='center', fontsize=12, fontweight='bold', 
       bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
ax.text(45, 0.5, 'Central', ha='center', va='center', fontsize=12, fontweight='bold',
       bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
ax.text(75, 0.5, 'South', ha='center', va='center', fontsize=12, fontweight='bold',
       bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))

plt.tight_layout()
plt.savefig(os.path.join(OUTPUT_DIR, 'background_erosion_heatmap.png'), dpi=300, bbox_inches='tight')
print(f"   ✓ Saved: background_erosion_heatmap.png")
plt.close()

# =============================================================================
# IDENTIFY KEY DOMAINS
# =============================================================================

print("\n6. Identifying Key Domains:")
print("=" * 80)

print("\n✓ Best Performing Domains (Smallest Absolute Error):")
best = results.nsmallest(5, 'abs_error')
for _, row in best.iterrows():
    print(f"   Domain {row['domain']:2.0f}: Obs={row['observed']:6.2f}, Mod={row['modeled']:6.2f}, "
          f"Error={row['error']:6.2f} m/yr")

print("\n✗ Worst Performing Domains (Largest Absolute Error):")
worst = results.nlargest(5, 'abs_error')
for _, row in worst.iterrows():
    print(f"   Domain {row['domain']:2.0f}: Obs={row['observed']:6.2f}, Mod={row['modeled']:6.2f}, "
          f"Error={row['error']:6.2f} m/yr")

print("\n⚠️  Domains with Wrong Sign (Model predicts opposite of reality):")
wrong_sign = results[results['category'] == 'Wrong Sign']
if len(wrong_sign) > 0:
    for _, row in wrong_sign.iterrows():
        print(f"   Domain {row['domain']:2.0f}: Obs={row['observed']:6.2f}, Mod={row['modeled']:6.2f}")
else:
    print("   None! Model at least gets the direction right everywhere.")

print("\n" + "=" * 80)
print("ANALYSIS COMPLETE")
print("=" * 80)
print(f"\nAll outputs saved to: {OUTPUT_DIR}")
print(f"  - domain_by_domain_results.csv (detailed metrics)")
print(f"  - observed_vs_modeled_scatter.png (1:1 comparison)")
print(f"  - spatial_patterns_along_island.png (4-panel spatial view)")
print(f"  - error_distribution_analysis.png (statistical summaries)")
print(f"  - background_erosion_heatmap.png (visual of needed corrections)")
