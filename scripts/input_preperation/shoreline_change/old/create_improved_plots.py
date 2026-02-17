"""
Improved Shoreline Change Rate Visualizations
Creates cleaner, more readable plots with better color schemes
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colors as mcolors

# Load the calculated rates
domain_rates = pd.read_csv('../domain_shoreline_change_rates.csv')

# Get available time periods (exclude the 'Domain' column)
rate_cols = [col for col in domain_rates.columns if 'EPR_' in col]

# ============================================================================
# OPTION 1: Better colors with sequential palette
# ============================================================================

fig, ax = plt.subplots(figsize=(18, 8))

# Use a sequential colormap that shows progression through time
# Darker = earlier, lighter = more recent
colors_seq = ['#08519c', '#3182bd', '#6baed6', '#c6dbef']  # Blues, dark to light

# Or use a diverging palette to show temporal change
# colors_seq = ['#d73027', '#fc8d59', '#91bfdb', '#4575b4']  # Red-blue diverging

period_names = {
    'EPR_1978_1987': '1978-1987',
    'EPR_1987_1997': '1987-1997', 
    'EPR_1997_2009': '1997-2009',
    'EPR_2009_2019': '2009-2019'
}

for i, (col, color) in enumerate(zip(['EPR_1978_1987', 'EPR_1987_1997', 'EPR_1997_2009', 'EPR_2009_2019'], colors_seq)):
    if col in domain_rates.columns:
        label = period_names.get(col, col)
        ax.plot(domain_rates['Domain'], domain_rates[col], 
                label=label, linewidth=2.5, 
                color=color, alpha=0.9, zorder=10-i)

# Rodanthe highlight
ax.axvspan(70, 90, alpha=0.12, color='#d62728', label='Rodanthe', zorder=0)
ax.axhline(y=0, color='black', linestyle='--', linewidth=1.5, alpha=0.7, zorder=1)

ax.set_xlabel('Domain Number', fontsize=14, fontweight='bold')
ax.set_ylabel('Shoreline Change Rate (m/yr)', fontsize=14, fontweight='bold')
ax.set_title('Shoreline Change Rates: Temporal Evolution', 
             fontsize=16, fontweight='bold', pad=20)
ax.legend(loc='upper right', fontsize=11, framealpha=0.95, edgecolor='gray')
ax.grid(True, alpha=0.25, linestyle=':', linewidth=0.8)
ax.set_xlim(domain_rates['Domain'].min() - 1, domain_rates['Domain'].max() + 1)

plt.tight_layout()
plt.savefig('improved_v1_sequential_colors.png', dpi=300, bbox_inches='tight')
print("✓ Saved: improved_v1_sequential_colors.png")
plt.close()

# ============================================================================
# OPTION 2: Faceted/subplot approach - one panel per period
# ============================================================================

fig, axes = plt.subplots(4, 1, figsize=(16, 12), sharex=True)

colors_distinct = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728']
periods = ['EPR_1978_1987', 'EPR_1987_1997', 'EPR_1997_2009', 'EPR_2009_2019']
period_labels = ['1978-1987', '1987-1997', '1997-2009', '2009-2019']

for ax, col, color, label in zip(axes, periods, colors_distinct, period_labels):
    if col in domain_rates.columns:
        ax.plot(domain_rates['Domain'], domain_rates[col], 
                linewidth=2.5, color=color)
        
        # Rodanthe highlight
        ax.axvspan(70, 90, alpha=0.12, color='red', zorder=0)
        ax.axhline(y=0, color='black', linestyle='--', linewidth=1, alpha=0.6)
        
        # Period label
        ax.text(0.02, 0.95, label, transform=ax.transAxes, 
                fontsize=12, fontweight='bold', va='top',
                bbox=dict(boxstyle='round', facecolor='white', alpha=0.8, edgecolor=color, linewidth=2))
        
        ax.set_ylabel('Rate (m/yr)', fontsize=11, fontweight='bold')
        ax.grid(True, alpha=0.25, linestyle=':', linewidth=0.8)
        ax.set_ylim(domain_rates[col].min() - 1, domain_rates[col].max() + 1)

axes[-1].set_xlabel('Domain Number', fontsize=14, fontweight='bold')
fig.suptitle('Shoreline Change Rates by Time Period', fontsize=16, fontweight='bold', y=0.995)
plt.tight_layout()
plt.savefig('improved_v2_faceted.png', dpi=300, bbox_inches='tight')
print("✓ Saved: improved_v2_faceted.png")
plt.close()

# ============================================================================
# OPTION 3: Two-panel comparison - Early vs Recent
# ============================================================================

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(18, 6), sharey=True)

# Left panel: Early periods
if 'EPR_1978_1987' in domain_rates.columns:
    ax1.plot(domain_rates['Domain'], domain_rates['EPR_1978_1987'], 
            label='1978-1987', linewidth=3, color='#2166ac', marker='o', markersize=3)
if 'EPR_1987_1997' in domain_rates.columns:
    ax1.plot(domain_rates['Domain'], domain_rates['EPR_1987_1997'], 
            label='1987-1997', linewidth=3, color='#67a9cf', marker='s', markersize=3)

ax1.axvspan(70, 90, alpha=0.12, color='red', zorder=0)
ax1.axhline(y=0, color='black', linestyle='--', linewidth=1.5, alpha=0.7)
ax1.set_xlabel('Domain Number', fontsize=13, fontweight='bold')
ax1.set_ylabel('Shoreline Change Rate (m/yr)', fontsize=13, fontweight='bold')
ax1.set_title('Earlier Period (1978-1997)', fontsize=14, fontweight='bold', pad=15)
ax1.legend(loc='best', fontsize=11, framealpha=0.95)
ax1.grid(True, alpha=0.25, linestyle=':', linewidth=0.8)

# Right panel: Recent periods
if 'EPR_1997_2009' in domain_rates.columns:
    ax2.plot(domain_rates['Domain'], domain_rates['EPR_1997_2009'], 
            label='1997-2009', linewidth=3, color='#ef8a62', marker='^', markersize=3)
if 'EPR_2009_2019' in domain_rates.columns:
    ax2.plot(domain_rates['Domain'], domain_rates['EPR_2009_2019'], 
            label='2009-2019', linewidth=3, color='#b2182b', marker='d', markersize=3)

ax2.axvspan(70, 90, alpha=0.12, color='red', label='Rodanthe', zorder=0)
ax2.axhline(y=0, color='black', linestyle='--', linewidth=1.5, alpha=0.7)
ax2.set_xlabel('Domain Number', fontsize=13, fontweight='bold')
ax2.set_title('Recent Period (1997-2019)', fontsize=14, fontweight='bold', pad=15)
ax2.legend(loc='best', fontsize=11, framealpha=0.95)
ax2.grid(True, alpha=0.25, linestyle=':', linewidth=0.8)

plt.tight_layout()
plt.savefig('improved_v3_early_vs_recent.png', dpi=300, bbox_inches='tight')
print("✓ Saved: improved_v3_early_vs_recent.png")
plt.close()

# ============================================================================
# OPTION 4: Clean overlay with different line styles
# ============================================================================

fig, ax = plt.subplots(figsize=(18, 8))

# Use both color AND line style for differentiation
styles = ['-', '--', '-.', ':']
colors_clean = ['#1b9e77', '#d95f02', '#7570b3', '#e7298a']  # ColorBrewer Dark2

for i, (col, style, color) in enumerate(zip(['EPR_1978_1987', 'EPR_1987_1997', 'EPR_1997_2009', 'EPR_2009_2019'], 
                                             styles, colors_clean)):
    if col in domain_rates.columns:
        label = period_names.get(col, col)
        ax.plot(domain_rates['Domain'], domain_rates[col], 
                label=label, linewidth=2.8, linestyle=style,
                color=color, alpha=0.85)

ax.axvspan(70, 90, alpha=0.1, color='#e31a1c', label='Rodanthe', zorder=0)
ax.axhline(y=0, color='black', linestyle='-', linewidth=1.2, alpha=0.5)

ax.set_xlabel('Domain Number', fontsize=14, fontweight='bold')
ax.set_ylabel('Shoreline Change Rate (m/yr)', fontsize=14, fontweight='bold')
ax.set_title('Shoreline Change Rates Across Time Periods - Hatteras Island', 
             fontsize=16, fontweight='bold', pad=20)
ax.legend(loc='upper right', fontsize=12, framealpha=0.95, edgecolor='gray', fancybox=True)
ax.grid(True, alpha=0.25, linestyle=':', linewidth=0.8)
ax.set_xlim(domain_rates['Domain'].min() - 1, domain_rates['Domain'].max() + 1)

plt.tight_layout()
plt.savefig('improved_v4_line_styles.png', dpi=300, bbox_inches='tight')
print("✓ Saved: improved_v4_line_styles.png")
plt.close()

# ============================================================================
# OPTION 5: Focus on calibration periods only (cleaner)
# ============================================================================

if 'EPR_1978_1997' in domain_rates.columns and 'EPR_1997_2019' in domain_rates.columns:
    fig, ax = plt.subplots(figsize=(18, 8))
    
    ax.plot(domain_rates['Domain'], domain_rates['EPR_1978_1997'], 
            label='1978-1997 (Calibration Period 1)', linewidth=4, 
            color='#2166ac', marker='o', markersize=4, alpha=0.9)
    
    ax.plot(domain_rates['Domain'], domain_rates['EPR_1997_2019'], 
            label='1997-2019 (Calibration Period 2)', linewidth=4, 
            color='#b2182b', marker='s', markersize=4, alpha=0.9)
    
    # Rodanthe highlight
    ax.axvspan(70, 90, alpha=0.15, color='orange', label='Rodanthe Area', zorder=0)
    ax.axhline(y=0, color='black', linestyle='--', linewidth=1.5, alpha=0.7)
    
    ax.set_xlabel('Domain Number', fontsize=14, fontweight='bold')
    ax.set_ylabel('Shoreline Change Rate (m/yr)', fontsize=14, fontweight='bold')
    ax.set_title('CASCADE Calibration Periods: Observed Shoreline Change Rates', 
                 fontsize=16, fontweight='bold', pad=20)
    ax.legend(loc='upper right', fontsize=13, framealpha=0.95, edgecolor='gray')
    ax.grid(True, alpha=0.25, linestyle=':', linewidth=0.8)
    ax.set_xlim(domain_rates['Domain'].min() - 1, domain_rates['Domain'].max() + 1)
    
    plt.tight_layout()
    plt.savefig('improved_v5_calibration_only.png', dpi=300, bbox_inches='tight')
    print("✓ Saved: improved_v5_calibration_only.png")
    plt.close()

print("\n" + "="*70)
print("CREATED 5 IMPROVED VISUALIZATION OPTIONS:")
print("="*70)
print("1. improved_v1_sequential_colors.png - Sequential color scheme")
print("2. improved_v2_faceted.png - Separate panel for each period")
print("3. improved_v3_early_vs_recent.png - Two-panel early vs recent comparison")
print("4. improved_v4_line_styles.png - Different line styles + colors")
print("5. improved_v5_calibration_only.png - Just your two calibration periods")
print("\nReview these and pick the one(s) that work best for your presentation!")
