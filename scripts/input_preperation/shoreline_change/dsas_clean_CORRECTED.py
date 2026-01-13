"""
Corrected DSAS cleaning script for Hatteras shoreline change data

Key fixes:
1. Uses LRR (Linear Regression Rate, m/yr) instead of NSM (Net Shoreline Movement, total m)
2. Preserves all domains (no quality filtering that removes entire domains)
3. LRR is already an annual rate, so no need to divide by time period
"""

import pandas as pd
import numpy as np

# Paths – update these for your system
RAW_DSAS_CSV = r"C:\Users\hanna\PycharmProjects\CASCADE\data\hatteras_init\shoreline_change\dsas_1978_1997_rates.csv"
OUTPUT_DOMAIN_CSV = r"C:\Users\hanna\PycharmProjects\CASCADE\data\hatteras_init\shoreline_change\dsas_1978_1997_CLEAN.csv"

# Load
df = pd.read_csv(RAW_DSAS_CSV)

print("="*60)
print("RAW DATA SUMMARY")
print("="*60)
print(f"Total transects: {len(df)}")
print(f"Transects with domain_id: {df['domain_id'].notna().sum()}")
print(f"Unique domains: {df['domain_id'].nunique()}")

# Keep only relevant columns
cols_keep = ["domain_id", "TransectId", "TCD", "LRR", "LR2", "NSM"]
df = df[cols_keep].copy()

# Remove rows without domain_id
df = df[df['domain_id'].notna()].copy()

# IMPORTANT: Keep only domains 1-90 (your study area)
# Domains 91-120 are collaborator's area and not part of your CASCADE model
df = df[df['domain_id'] <= 90].copy()
print(f"\nFiltered to domains 1-90 only (your study area): {len(df)} transects")

print(f"\nLRR (annual rate) statistics:")
print(f"  Mean: {df['LRR'].mean():.3f} m/yr")
print(f"  Std: {df['LRR'].std():.3f} m/yr")
print(f"  Min: {df['LRR'].min():.3f} m/yr")
print(f"  Max: {df['LRR'].max():.3f} m/yr")

print(f"\nLR2 (R²) statistics:")
print(f"  Mean: {df['LR2'].mean():.3f}")
print(f"  Median: {df['LR2'].median():.3f}")
print(f"  Transects with LR2 < 0.3: {(df['LR2'] < 0.3).sum()} ({(df['LR2'] < 0.3).sum()/len(df)*100:.1f}%)")

# Optional quality filter: Remove only very poor regression fits
# This is more lenient than before to preserve more domains
df_filtered = df[df["LR2"] >= 0.3].copy()

print(f"\n" + "="*60)
print("AFTER FILTERING (LR2 >= 0.3)")
print("="*60)
print(f"Transects remaining: {len(df_filtered)} ({len(df_filtered)/len(df)*100:.1f}%)")
print(f"Domains remaining: {df_filtered['domain_id'].nunique()}")

# Check for missing domains after filtering
all_domains_original = set(df['domain_id'].unique())
domains_after_filter = set(df_filtered['domain_id'].unique())
missing_domains = all_domains_original - domains_after_filter

if missing_domains:
    print(f"\n⚠ WARNING: {len(missing_domains)} domains lost due to filtering:")
    print(f"Missing domains: {sorted([int(d) for d in missing_domains])}")
    
    # Add back missing domains using unfiltered data
    print("\nRecovering missing domains with unfiltered data...")
    for domain in sorted(missing_domains):
        domain_data = df[df['domain_id'] == domain]
        print(f"  Domain {int(domain)}: {len(domain_data)} transects (LR2: {domain_data['LR2'].min():.2f}-{domain_data['LR2'].max():.2f})")
        df_filtered = pd.concat([df_filtered, domain_data])
    
    print(f"\n✓ All {len(all_domains_original)} domains preserved!")
else:
    print("\n✓ All domains preserved after filtering!")

# Sort (not strictly required, but nice)
df_filtered = df_filtered.sort_values(["domain_id", "TCD"])

# Aggregate to domain: mean LRR (m/yr) over all transects in that domain
# LRR is already an annual rate, so we just average across transects
domain_means = (
    df_filtered.groupby("domain_id")["LRR"]
    .agg(['mean', 'std', 'count'])
    .reset_index()
    .rename(columns={
        'mean': 'annual_rate_m_per_yr',
        'std': 'std_dev_m_per_yr',
        'count': 'n_transects'
    })
)

print("\n" + "="*60)
print("DOMAIN-LEVEL AGGREGATION")
print("="*60)
print(f"Total domains: {len(domain_means)}")
print(f"\nAnnual rate (LRR) per domain:")
print(f"  Mean: {domain_means['annual_rate_m_per_yr'].mean():.3f} m/yr")
print(f"  Std: {domain_means['annual_rate_m_per_yr'].std():.3f} m/yr")
print(f"  Min: {domain_means['annual_rate_m_per_yr'].min():.3f} m/yr (domain {domain_means.loc[domain_means['annual_rate_m_per_yr'].idxmin(), 'domain_id']:.0f})")
print(f"  Max: {domain_means['annual_rate_m_per_yr'].max():.3f} m/yr (domain {domain_means.loc[domain_means['annual_rate_m_per_yr'].idxmax(), 'domain_id']:.0f})")

# Save for overlay with model
# Simple version (just domain_id and rate)
domain_means[['domain_id', 'annual_rate_m_per_yr']].to_csv(OUTPUT_DOMAIN_CSV, index=False)

# Also save detailed version with stats
OUTPUT_DETAILED = OUTPUT_DOMAIN_CSV.replace('.csv', '_DETAILED.csv')
domain_means.to_csv(OUTPUT_DETAILED, index=False)

print(f"\n✓ Saved simple output: {OUTPUT_DOMAIN_CSV}")
print(f"✓ Saved detailed output: {OUTPUT_DETAILED}")

print("\n" + "="*60)
print("Sample of output:")
print("="*60)
print(domain_means.head(10))

print("\n" + "="*60)
print("DONE!")
print("="*60)
