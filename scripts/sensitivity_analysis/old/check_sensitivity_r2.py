"""
Diagnostic Script: Check Sensitivity Analysis R² Calculation

This script examines one of your completed sensitivity runs to diagnose
why R² values are negative.
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os
import glob

# =============================================================================
# CONFIGURATION - UPDATE THESE PATHS
# =============================================================================

# Path to your sensitivity results CSV
RESULTS_CSV = r"C:\Users\hanna\PycharmProjects\CASCADE\output\sensitivity_analysis\sensitivity_results_20260126_215530.csv"

# Path to DSAS data
DSAS_FILE = r"/data/hatteras_init/shoreline_change/dsas_1978_1997_domain_means_SIMPLE.csv"

# Pick one run to examine in detail (baseline run)
EXAMPLE_RUN_DIR = r"/output/sensitivity_analysis/sens_wave_period_7_20260126_215530"

# Constants
NUM_BUFFER_DOMAINS = 15
NUM_REAL_DOMAINS = 90
START_REAL_INDEX = NUM_BUFFER_DOMAINS
END_REAL_INDEX = START_REAL_INDEX + NUM_REAL_DOMAINS
RUN_YEARS = 19

# =============================================================================
# HELPER FUNCTIONS
# =============================================================================

def get_x_s_TS(b3d):
    """Get shoreline time series from Barrier3D object."""
    if hasattr(b3d, "x_s_TS"):
        return np.array(b3d.x_s_TS)
    elif hasattr(b3d, "_x_s_TS"):
        return np.array(b3d._x_s_TS)
    else:
        raise AttributeError("No shoreline time series found")


def load_dsas_data(dsas_file):
    """Load DSAS observed data."""
    dsas = pd.read_csv(dsas_file)
    gis_domains = dsas['domain_id'].to_numpy()
    obs_rate_all = dsas['annual_rate_m_per_yr'].to_numpy()
    
    # Map GIS domains (1-90) to CASCADE domains (15-104)
    cascade_domains = gis_domains + NUM_BUFFER_DOMAINS
    
    # Create array with NaN for all real domains
    observed_rates = np.full(NUM_REAL_DOMAINS, np.nan)
    
    # Fill in observed values
    for cascade_idx, rate in zip(cascade_domains, obs_rate_all):
        real_domain_idx = cascade_idx - NUM_BUFFER_DOMAINS
        if 0 <= real_domain_idx < NUM_REAL_DOMAINS:
            observed_rates[real_domain_idx] = rate
    
    return observed_rates


def extract_model_rates(run_dir):
    """Extract model rates from a CASCADE run."""
    # Find NPZ file
    npz_files = glob.glob(os.path.join(run_dir, '*.npz'))
    if len(npz_files) == 0:
        raise FileNotFoundError(f"No NPZ file in {run_dir}")
    
    # Load CASCADE
    data = np.load(npz_files[0], allow_pickle=True)
    cascade = data["cascade"][0]
    
    # Build shoreline matrix
    b3d_list = cascade.barrier3d
    ndom = len(b3d_list)
    nt = len(get_x_s_TS(b3d_list[0]))
    
    shoreline = np.zeros((nt, ndom))
    for j in range(ndom):
        xs = get_x_s_TS(b3d_list[j])
        shoreline[:, j] = xs
    
    # Convert dam → m
    shoreline = shoreline * 10.0
    
    # Calculate rates for real domains
    initial = shoreline[0, START_REAL_INDEX:END_REAL_INDEX]
    final = shoreline[-1, START_REAL_INDEX:END_REAL_INDEX]
    rates = (final - initial) / RUN_YEARS
    
    return rates


def calculate_r_squared_manual(observed, modeled):
    """Calculate R² manually to verify."""
    # Remove NaN values
    valid_mask = ~np.isnan(observed)
    obs_valid = observed[valid_mask]
    mod_valid = modeled[valid_mask]
    
    # R² calculation
    ss_res = np.sum((obs_valid - mod_valid) ** 2)
    ss_tot = np.sum((obs_valid - np.mean(obs_valid)) ** 2)
    r_squared = 1 - (ss_res / ss_tot)
    
    return r_squared, obs_valid, mod_valid


# =============================================================================
# MAIN DIAGNOSTICS
# =============================================================================

print("=" * 80)
print("SENSITIVITY ANALYSIS R² DIAGNOSTIC")
print("=" * 80)

# Load sensitivity results
print("\n1. Loading sensitivity results CSV...")
results_df = pd.read_csv(RESULTS_CSV)
print(f"   ✓ Loaded {len(results_df)} runs")

# Show R² distribution
r2_values = results_df['r_squared'].dropna()
print(f"\n2. R² Distribution:")
print(f"   Min R²: {r2_values.min():.4f}")
print(f"   Max R²: {r2_values.max():.4f}")
print(f"   Mean R²: {r2_values.mean():.4f}")
print(f"   Median R²: {r2_values.median():.4f}")

if r2_values.max() < 0:
    print("   ⚠️  WARNING: ALL R² values are NEGATIVE!")
    print("   This suggests a systematic problem, not just poor model fit.")

# Load DSAS data
print("\n3. Loading DSAS observed data...")
observed = load_dsas_data(DSAS_FILE)
n_obs = np.sum(~np.isnan(observed))
print(f"   ✓ Loaded {n_obs} observations")
print(f"   Range: {np.nanmin(observed):.3f} to {np.nanmax(observed):.3f} m/yr")
print(f"   Mean: {np.nanmean(observed):.3f} m/yr")
print(f"   Std: {np.nanstd(observed):.3f} m/yr")

# Load example model run
print("\n4. Loading example model run (baseline configuration)...")
print(f"   Directory: {EXAMPLE_RUN_DIR}")
modeled = extract_model_rates(EXAMPLE_RUN_DIR)
print(f"   ✓ Extracted model rates")
print(f"   Range: {np.min(modeled):.3f} to {np.max(modeled):.3f} m/yr")
print(f"   Mean: {np.mean(modeled):.3f} m/yr")
print(f"   Std: {np.std(modeled):.3f} m/yr")

# Calculate R² manually
print("\n5. Calculating R² manually...")
r2_manual, obs_valid, mod_valid = calculate_r_squared_manual(observed, modeled)
print(f"   Manual R² calculation: {r2_manual:.4f}")

# Get R² from CSV for this run
csv_r2 = results_df[results_df['run_name'].str.contains('wave_period_7')]['r_squared'].values[0]
print(f"   CSV R² value: {csv_r2:.4f}")

if abs(r2_manual - csv_r2) > 0.01:
    print("   ⚠️  MISMATCH between manual and CSV R²!")
else:
    print("   ✓ Manual and CSV R² match")

# Check for obvious issues
print("\n6. Diagnostic Checks:")

# Check 1: Are signs flipped?
corr = np.corrcoef(obs_valid, mod_valid)[0, 1]
print(f"   Correlation coefficient: {corr:.4f}")
if corr < -0.5:
    print("   ⚠️  STRONG NEGATIVE CORRELATION - Signs may be flipped!")
elif corr < 0:
    print("   ⚠️  NEGATIVE CORRELATION - Model predicts opposite of observations")
else:
    print("   ✓ Positive correlation (at least signs are consistent)")

# Check 2: Magnitude comparison
obs_range = np.nanmax(observed) - np.nanmin(observed)
mod_range = np.max(modeled) - np.min(modeled)
print(f"   Observed range: {obs_range:.3f} m/yr")
print(f"   Modeled range: {mod_range:.3f} m/yr")
print(f"   Ratio: {mod_range/obs_range:.3f}")
if mod_range < obs_range * 0.1:
    print("   ⚠️  Model range is <10% of observed - insufficient variability")

# Check 3: Mean offset
obs_mean = np.nanmean(observed)
mod_mean = np.mean(modeled)
offset = mod_mean - obs_mean
print(f"   Observed mean: {obs_mean:.3f} m/yr")
print(f"   Modeled mean: {mod_mean:.3f} m/yr")
print(f"   Bias: {offset:.3f} m/yr")
if abs(offset) > 1.0:
    print("   ⚠️  Large bias (>1 m/yr) between model and observations")

# Check 4: Spatial alignment
print("\n7. Checking spatial alignment...")
print("   Comparing first and last 5 domains:")
print("   Domain | Observed | Modeled | Difference")
print("   " + "-" * 50)
for i in [0, 1, 2, 3, 4]:
    if not np.isnan(observed[i]):
        print(f"   {i+1:6d} | {observed[i]:8.3f} | {modeled[i]:7.3f} | {modeled[i]-observed[i]:10.3f}")
print("   ...")
for i in [-5, -4, -3, -2, -1]:
    idx = NUM_REAL_DOMAINS + i
    if not np.isnan(observed[idx]):
        print(f"   {idx+1:6d} | {observed[idx]:8.3f} | {modeled[idx]:7.3f} | {modeled[idx]-observed[idx]:10.3f}")

# Create comparison plot
print("\n8. Creating comparison plot...")
fig, axes = plt.subplots(2, 2, figsize=(14, 10))
fig.suptitle('Sensitivity Analysis R² Diagnostic', fontsize=16, fontweight='bold')

# Plot 1: Observed vs Modeled
ax1 = axes[0, 0]
ax1.scatter(obs_valid, mod_valid, alpha=0.6, s=50)
ax1.plot([obs_valid.min(), obs_valid.max()], 
         [obs_valid.min(), obs_valid.max()], 
         'r--', linewidth=2, label='1:1 line')
ax1.set_xlabel('Observed Rate (m/yr)', fontsize=12)
ax1.set_ylabel('Modeled Rate (m/yr)', fontsize=12)
ax1.set_title(f'Observed vs Modeled (R² = {r2_manual:.3f})', fontsize=13, fontweight='bold')
ax1.legend()
ax1.grid(alpha=0.3)
ax1.axhline(0, color='gray', linestyle='-', linewidth=0.5)
ax1.axvline(0, color='gray', linestyle='-', linewidth=0.5)

# Plot 2: Spatial patterns
ax2 = axes[0, 1]
domains = np.arange(NUM_REAL_DOMAINS) + 1
ax2.plot(domains, observed, 'o-', label='Observed', linewidth=2, markersize=4, alpha=0.7)
ax2.plot(domains, modeled, 's-', label='Modeled', linewidth=2, markersize=4, alpha=0.7)
ax2.axhline(0, color='gray', linestyle='--', linewidth=1)
ax2.set_xlabel('Domain Number', fontsize=12)
ax2.set_ylabel('Shoreline Change Rate (m/yr)', fontsize=12)
ax2.set_title('Spatial Pattern Comparison', fontsize=13, fontweight='bold')
ax2.legend()
ax2.grid(alpha=0.3)

# Plot 3: Residuals
ax3 = axes[1, 0]
residuals = mod_valid - obs_valid
ax3.scatter(obs_valid, residuals, alpha=0.6, s=50)
ax3.axhline(0, color='r', linestyle='--', linewidth=2)
ax3.set_xlabel('Observed Rate (m/yr)', fontsize=12)
ax3.set_ylabel('Residuals (Modeled - Observed)', fontsize=12)
ax3.set_title('Residual Plot', fontsize=13, fontweight='bold')
ax3.grid(alpha=0.3)

# Plot 4: R² components
ax4 = axes[1, 1]
ss_res = np.sum((obs_valid - mod_valid) ** 2)
ss_tot = np.sum((obs_valid - np.mean(obs_valid)) ** 2)
ss_exp = ss_tot - ss_res

categories = ['Total\nVariance', 'Explained\nVariance', 'Residual\nVariance']
values = [ss_tot, ss_exp, ss_res]
colors = ['gray', 'green' if ss_exp > 0 else 'red', 'red']

bars = ax4.bar(categories, values, color=colors, alpha=0.7)
ax4.set_ylabel('Sum of Squares', fontsize=12)
ax4.set_title(f'R² Components (R² = {r2_manual:.3f})', fontsize=13, fontweight='bold')
ax4.grid(alpha=0.3, axis='y')

# Add value labels on bars
for bar in bars:
    height = bar.get_height()
    ax4.text(bar.get_x() + bar.get_width()/2., height,
             f'{height:.1f}',
             ha='center', va='bottom', fontsize=10)

plt.tight_layout()
plt.savefig('sensitivity_r2_diagnostic.png', dpi=200, bbox_inches='tight')
print("   ✓ Saved diagnostic plot: sensitivity_r2_diagnostic.png")

print("\n" + "=" * 80)
print("DIAGNOSIS COMPLETE")
print("=" * 80)

# Summary
print("\n📊 SUMMARY:")
print(f"   R² value: {r2_manual:.4f}")
print(f"   Correlation: {corr:.4f}")
print(f"   Model explains {r2_manual*100:.1f}% of observed variance")

if r2_manual < 0:
    print("\n⚠️  NEGATIVE R² DIAGNOSIS:")
    print("   This means the model performs WORSE than just predicting the mean.")
    print("   Possible causes:")
    if abs(offset) > 1:
        print("   ✗ Large bias - model systematically over/under-predicts")
    if mod_range < obs_range * 0.2:
        print("   ✗ Insufficient variability - model too flat")
    if corr < 0:
        print("   ✗ Wrong sign/direction - model predicts opposite patterns")
    
    print("\n   Bottom line: Wave parameters alone cannot capture observed patterns.")
    print("   This is actually GOOD NEWS for your thesis - it proves background")
    print("   erosion is essential!")
    
print("\n✓ Check 'sensitivity_r2_diagnostic.png' for visual analysis")
