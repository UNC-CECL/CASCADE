"""
Format the DSAS domain statistics already calculated in GIS for the Python scripts.

    python scripts/figure_making/shoreline/dsas/dsas_from_gis.py

Reads the DSAS 1997-2019 rate table; writes the domain CSVs. Details: scripts/figure_making/shoreline/README.md.

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-09-30
"""

import pandas as pd
import numpy as np

# Anchored 2026-09-22
import sys as _sys
from pathlib import Path as _P
_sys.path.insert(0, str(next(_q for _q in _P(__file__).resolve().parents
                             if (_q / "pyproject.toml").exists()) / "scripts"))
from site_layer.hat_observed_rates import DSAS_ROOT  # noqa: E402

# --- CONFIG ------------------------------------------------------------------
RAW_DSAS_CSV = str(DSAS_ROOT / "dsas_1997_2019_rates.csv")
OUTPUT_DOMAIN_CSV = str(DSAS_ROOT / "dsas_1997_2019_domain_means.csv")

# Domain range for your study area
MIN_DOMAIN = 1
MAX_DOMAIN = 90
# -----------------------------------------------------------------------------

# Load and format data
print("="*70)
print("DSAS DATA FORMATTING - GIS-CALCULATED DOMAIN MEANS")
print("="*70)

# Load data
df = pd.read_csv(RAW_DSAS_CSV)

print(f"\nRaw data loaded: {len(df)} rows")
print(f"Columns: {', '.join(df.columns)}")

# Keep only rows with valid domain_id
df = df[df['domain_id'].notna()].copy()
print(f"After removing missing domain_id: {len(df)} rows")

# Filter to study area domains (1-90)
df = df[(df['domain_id'] >= MIN_DOMAIN) & (df['domain_id'] <= MAX_DOMAIN)].copy()
print(f"After filtering to domains {MIN_DOMAIN}-{MAX_DOMAIN}: {len(df)} rows")

# Check that we have all expected domains
expected_domains = set(range(MIN_DOMAIN, MAX_DOMAIN + 1))
actual_domains = set(df['domain_id'].unique())
missing_domains = expected_domains - actual_domains

if missing_domains:
    print(f"\n⚠️  WARNING: Missing {len(missing_domains)} domains!")
    print(f"Missing: {sorted(missing_domains)}")
else:
    print(f"\n✓ All {len(expected_domains)} domains present")

# Extract relevant columns

# Select columns of interest (you already have domain-level stats from GIS)
output_df = df[[
    'domain_id',
    'MEAN_LRR',      # Mean linear regression rate (m/yr)
    'STD_LRR',       # Standard deviation of LRR
    'N_TRANSECTS',   # Number of transects per domain
    'MIN_LRR',       # Minimum LRR in domain
    'MAX_LRR'        # Maximum LRR in domain
]].copy()

# Rename for clarity and consistency with CASCADE scripts
output_df = output_df.rename(columns={
    'MEAN_LRR': 'annual_rate_m_per_yr',
    'STD_LRR': 'std_dev_m_per_yr',
    'N_TRANSECTS': 'n_transects',
    'MIN_LRR': 'min_rate_m_per_yr',
    'MAX_LRR': 'max_rate_m_per_yr'
})

# Sort by domain
output_df = output_df.sort_values('domain_id').reset_index(drop=True)

# Summary statistics
print("\n" + "="*70)
print("DOMAIN-LEVEL STATISTICS")
print("="*70)

print(f"\nTotal domains: {len(output_df)}")
print(f"\nAnnual rate (LRR) statistics across all domains:")
print(f"  Mean:   {output_df['annual_rate_m_per_yr'].mean():>7.3f} m/yr")
print(f"  Median: {output_df['annual_rate_m_per_yr'].median():>7.3f} m/yr")
print(f"  Std:    {output_df['annual_rate_m_per_yr'].std():>7.3f} m/yr")
print(f"  Min:    {output_df['annual_rate_m_per_yr'].min():>7.3f} m/yr (domain {output_df.loc[output_df['annual_rate_m_per_yr'].idxmin(), 'domain_id']:.0f})")
print(f"  Max:    {output_df['annual_rate_m_per_yr'].max():>7.3f} m/yr (domain {output_df.loc[output_df['annual_rate_m_per_yr'].idxmax(), 'domain_id']:.0f})")

# Erosion vs accretion breakdown
eroding = (output_df['annual_rate_m_per_yr'] < 0).sum()
accreting = (output_df['annual_rate_m_per_yr'] > 0).sum()
stable = (output_df['annual_rate_m_per_yr'] == 0).sum()

print(f"\nDomain classification:")
print(f"  Eroding (< 0):    {eroding:>3} ({eroding/len(output_df)*100:>5.1f}%)")
print(f"  Accreting (> 0):  {accreting:>3} ({accreting/len(output_df)*100:>5.1f}%)")
print(f"  Stable (= 0):     {stable:>3} ({stable/len(output_df)*100:>5.1f}%)")

print(f"\nTransects per domain:")
print(f"  Mean:   {output_df['n_transects'].mean():>5.1f}")
print(f"  Median: {output_df['n_transects'].median():>5.0f}")
print(f"  Min:    {output_df['n_transects'].min():>5.0f}")
print(f"  Max:    {output_df['n_transects'].max():>5.0f}")

# Save output

# Save full version with all statistics
output_df.to_csv(OUTPUT_DOMAIN_CSV, index=False)
print(f"\n✓ Saved full comparison: {OUTPUT_DOMAIN_CSV}")

# Also save simplified version (just domain_id and rate) for quick use
simple_output = OUTPUT_DOMAIN_CSV.replace('.csv', '_SIMPLE.csv')
output_df[['domain_id', 'annual_rate_m_per_yr']].to_csv(simple_output, index=False)
print(f"✓ Saved simple comparison: {simple_output}")

# Preview output
print("\n" + "="*70)
print("SAMPLE OUTPUT (first 10 and last 10 domains):")
print("="*70)
print("\nFirst 10 domains:")
print(output_df.head(10).to_string(index=False))
print("\nLast 10 domains:")
print(output_df.tail(10).to_string(index=False))

# Spatial pattern analysis
print("\n" + "="*70)
print("SPATIAL PATTERNS")
print("="*70)

# Identify areas of significant erosion/accretion
very_erosive = output_df[output_df['annual_rate_m_per_yr'] < -2.0]
very_accretive = output_df[output_df['annual_rate_m_per_yr'] > 2.0]

if len(very_erosive) > 0:
    print(f"\nHighly erosive domains (< -2.0 m/yr): {len(very_erosive)}")
    print(f"  Domain range: {very_erosive['domain_id'].min():.0f} to {very_erosive['domain_id'].max():.0f}")
    print(f"  Mean rate: {very_erosive['annual_rate_m_per_yr'].mean():.3f} m/yr")

if len(very_accretive) > 0:
    print(f"\nHighly accretive domains (> +2.0 m/yr): {len(very_accretive)}")
    print(f"  Domain range: {very_accretive['domain_id'].min():.0f} to {very_accretive['domain_id'].max():.0f}")
    print(f"  Mean rate: {very_accretive['annual_rate_m_per_yr'].mean():.3f} m/yr")

print("\n" + "="*70)
print("DONE!")
print("="*70)
print("\nFormatted data ready for use in CASCADE scripts.")
print("Load with: df = pd.read_csv('dsas_1978_1997_domain_means.csv')")
