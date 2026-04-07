"""
Verification: Compare Calculated Rates to Known Values
Checks if your calculated shoreline change rates make sense
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# Load your calculated rates
domain_rates = pd.read_csv('domain_shoreline_change_rates.csv')

print("="*70)
print("VERIFICATION ANALYSIS")
print("="*70)

# ============================================================================
# CHECK 1: Compare to your CASCADE background erosion rate
# ============================================================================

print("\nCHECK 1: Comparison to CASCADE Background Erosion Rate")
print("-"*70)

# Your CASCADE background erosion rate (from your 1978-1997 calibration)
CASCADE_BG_EROSION = -1.091  # dam/yr
CASCADE_BG_EROSION_M = CASCADE_BG_EROSION * 10  # Convert to m/yr

print(f"CASCADE background erosion rate: {CASCADE_BG_EROSION_M:.2f} m/yr")

# Compare to your calculated 1978-1997 average
if 'EPR_1978_1997' in domain_rates.columns:
    mean_1978_1997 = domain_rates['EPR_1978_1997'].mean()
    print(f"Your calculated 1978-1997 average: {mean_1978_1997:.2f} m/yr")
    print(f"Difference: {abs(mean_1978_1997 - CASCADE_BG_EROSION_M):.2f} m/yr")
    
    if abs(mean_1978_1997 - CASCADE_BG_EROSION_M) < 3:
        print("✓ CLOSE MATCH - Good agreement!")
    else:
        print("⚠ DISCREPANCY - May need investigation")

# ============================================================================
# CHECK 2: Sign convention verification
# ============================================================================

print("\n" + "="*70)
print("CHECK 2: Sign Convention Verification")
print("-"*70)

# Rodanthe should be erosional (negative)
rodanthe = domain_rates[(domain_rates['Domain'] >= 77) & (domain_rates['Domain'] <= 83)]

print("\nRodanthe Area (Domains 77-83):")
for col in ['EPR_1978_1987', 'EPR_1987_1997', 'EPR_1997_2009', 'EPR_2009_2019']:
    if col in rodanthe.columns:
        mean_rate = rodanthe[col].mean()
        pct_erosional = (rodanthe[col] < 0).sum() / len(rodanthe) * 100
        
        if mean_rate < 0:
            status = "✓ EROSIONAL (correct)"
        else:
            status = "✗ ACCRETIONAL (WRONG!)"
        
        print(f"{col}: {mean_rate:+.2f} m/yr  |  {pct_erosional:.0f}% domains eroding  |  {status}")

# ============================================================================
# CHECK 3: Magnitude sanity check
# ============================================================================

print("\n" + "="*70)
print("CHECK 3: Magnitude Sanity Check")
print("-"*70)

print("\nTypical barrier island change rates:")
print("  Extreme erosion:  < -5 m/yr")
print("  Moderate erosion: -5 to -1 m/yr")
print("  Stable:           -1 to +1 m/yr")
print("  Moderate accretion: +1 to +5 m/yr")
print("  Extreme accretion:  > +5 m/yr")

print("\nYour data ranges:")
for col in ['EPR_1978_1987', 'EPR_1987_1997', 'EPR_1997_2009', 'EPR_2009_2019']:
    if col in domain_rates.columns:
        min_val = domain_rates[col].min()
        max_val = domain_rates[col].max()
        mean_val = domain_rates[col].mean()
        
        # Check for unrealistic values
        if abs(min_val) > 15 or abs(max_val) > 15:
            flag = "⚠ EXTREME VALUES - Check data"
        else:
            flag = "✓ Reasonable range"
        
        print(f"{col}: {min_val:+.2f} to {max_val:+.2f} m/yr (mean: {mean_val:+.2f})  |  {flag}")

# ============================================================================
# CHECK 4: Spatial pattern verification
# ============================================================================

print("\n" + "="*70)
print("CHECK 4: Known Spatial Patterns")
print("-"*70)

# Oregon Inlet area (domains 1-15) - should be erosional
oregon_inlet = domain_rates[domain_rates['Domain'] <= 15]
if len(oregon_inlet) > 0 and 'EPR_1978_1997' in oregon_inlet.columns:
    oregon_mean = oregon_inlet['EPR_1978_1997'].mean()
    print(f"\nOregon Inlet area (domains 1-15):")
    print(f"  Mean rate 1978-1997: {oregon_mean:+.2f} m/yr")
    if oregon_mean < -2:
        print("  ✓ Strong erosion as expected")
    elif oregon_mean < 0:
        print("  ✓ Erosional as expected")
    else:
        print("  ⚠ Unexpected accretion - check data")

# Buxton area (domains 50-65) - historically accretional
buxton = domain_rates[(domain_rates['Domain'] >= 50) & (domain_rates['Domain'] <= 65)]
if len(buxton) > 0 and 'EPR_1978_1997' in buxton.columns:
    buxton_mean = buxton['EPR_1978_1997'].mean()
    print(f"\nBuxton area (domains 50-65):")
    print(f"  Mean rate 1978-1997: {buxton_mean:+.2f} m/yr")
    if buxton_mean > 0:
        print("  ✓ Accretional as expected")
    else:
        print("  ⚠ Erosional - unexpected but possible")

# Rodanthe verification
if len(rodanthe) > 0 and 'EPR_1978_1997' in rodanthe.columns:
    rodanthe_mean = rodanthe['EPR_1978_1997'].mean()
    print(f"\nRodanthe area (domains 77-83):")
    print(f"  Mean rate 1978-1997: {rodanthe_mean:+.2f} m/yr")
    if rodanthe_mean < -2:
        print("  ✓ Strong erosion as expected")
    elif rodanthe_mean < 0:
        print("  ✓ Erosional as expected")
    else:
        print("  ✗ WRONG - Should be erosional!")

# ============================================================================
# CHECK 5: Temporal trends
# ============================================================================

print("\n" + "="*70)
print("CHECK 5: Temporal Trends")
print("-"*70)

print("\nIs erosion getting worse over time at Rodanthe?")
rodanthe_trends = []
for col in ['EPR_1978_1987', 'EPR_1987_1997', 'EPR_1997_2009', 'EPR_2009_2019']:
    if col in rodanthe.columns:
        mean_rate = rodanthe[col].mean()
        rodanthe_trends.append((col, mean_rate))
        print(f"  {col}: {mean_rate:+.2f} m/yr")

if len(rodanthe_trends) >= 2:
    first_period = rodanthe_trends[0][1]
    last_period = rodanthe_trends[-1][1]
    if last_period < first_period:
        print("  → Erosion accelerating (more negative over time)")
    else:
        print("  → Erosion slowing (less negative over time)")

# ============================================================================
# CREATE VERIFICATION PLOT
# ============================================================================

print("\n" + "="*70)
print("Creating verification plot...")
print("-"*70)

fig, axes = plt.subplots(2, 2, figsize=(16, 10))

# Plot 1: Histogram of rates (1978-1997)
if 'EPR_1978_1997' in domain_rates.columns:
    axes[0, 0].hist(domain_rates['EPR_1978_1997'], bins=30, edgecolor='black', alpha=0.7)
    axes[0, 0].axvline(x=0, color='red', linestyle='--', linewidth=2, label='Zero (stable)')
    axes[0, 0].axvline(x=CASCADE_BG_EROSION_M, color='green', linestyle='--', 
                       linewidth=2, label=f'CASCADE BG ({CASCADE_BG_EROSION_M:.1f} m/yr)')
    axes[0, 0].set_xlabel('Rate (m/yr)', fontweight='bold')
    axes[0, 0].set_ylabel('Number of Domains', fontweight='bold')
    axes[0, 0].set_title('Distribution of Rates (1978-1997)', fontweight='bold')
    axes[0, 0].legend()
    axes[0, 0].grid(True, alpha=0.3)

# Plot 2: Rodanthe temporal evolution
rodanthe_means = []
rodanthe_labels = []
for col in ['EPR_1978_1987', 'EPR_1987_1997', 'EPR_1997_2009', 'EPR_2009_2019']:
    if col in rodanthe.columns:
        rodanthe_means.append(rodanthe[col].mean())
        rodanthe_labels.append(col.replace('EPR_', ''))

axes[0, 1].plot(range(len(rodanthe_means)), rodanthe_means, 'o-', linewidth=3, markersize=10)
axes[0, 1].axhline(y=0, color='black', linestyle='--', linewidth=1, alpha=0.5)
axes[0, 1].set_xticks(range(len(rodanthe_labels)))
axes[0, 1].set_xticklabels(rodanthe_labels, rotation=45, ha='right')
axes[0, 1].set_ylabel('Mean Rate (m/yr)', fontweight='bold')
axes[0, 1].set_title('Rodanthe Erosion Trend Over Time', fontweight='bold')
axes[0, 1].grid(True, alpha=0.3)

# Plot 3: Comparison to CASCADE background rate
if 'EPR_1978_1997' in domain_rates.columns:
    axes[1, 0].scatter(domain_rates['Domain'], domain_rates['EPR_1978_1997'], 
                       alpha=0.6, s=50, label='Calculated rates')
    axes[1, 0].axhline(y=CASCADE_BG_EROSION_M, color='red', linestyle='--', 
                       linewidth=2, label=f'CASCADE BG ({CASCADE_BG_EROSION_M:.1f} m/yr)')
    axes[1, 0].axhline(y=0, color='black', linestyle='-', linewidth=1, alpha=0.3)
    axes[1, 0].axvspan(77, 83, alpha=0.1, color='orange', label='Rodanthe')
    axes[1, 0].set_xlabel('Domain Number', fontweight='bold')
    axes[1, 0].set_ylabel('Rate (m/yr)', fontweight='bold')
    axes[1, 0].set_title('Rates vs CASCADE Background (1978-1997)', fontweight='bold')
    axes[1, 0].legend()
    axes[1, 0].grid(True, alpha=0.3)

# Plot 4: Erosion vs Accretion breakdown
periods = ['1978-1987', '1987-1997', '1997-2009', '2009-2019']
erosional_pct = []
for col in ['EPR_1978_1987', 'EPR_1987_1997', 'EPR_1997_2009', 'EPR_2009_2019']:
    if col in domain_rates.columns:
        pct = (domain_rates[col] < 0).sum() / len(domain_rates) * 100
        erosional_pct.append(pct)

axes[1, 1].bar(range(len(erosional_pct)), erosional_pct, edgecolor='black', alpha=0.7)
axes[1, 1].axhline(y=50, color='red', linestyle='--', linewidth=2, label='50% threshold')
axes[1, 1].set_xticks(range(len(periods)))
axes[1, 1].set_xticklabels(periods, rotation=45, ha='right')
axes[1, 1].set_ylabel('% Domains Eroding', fontweight='bold')
axes[1, 1].set_title('Percentage of Island Experiencing Erosion', fontweight='bold')
axes[1, 1].set_ylim(0, 100)
axes[1, 1].legend()
axes[1, 1].grid(True, alpha=0.3, axis='y')

plt.tight_layout()
plt.savefig('verification_diagnostics.png', dpi=300, bbox_inches='tight')
print("✓ Saved: verification_diagnostics.png")
plt.close()

print("\n" + "="*70)
print("VERIFICATION COMPLETE")
print("="*70)
print("\nReview the diagnostics above and the verification_diagnostics.png plot")
print("to confirm your results make sense!")
