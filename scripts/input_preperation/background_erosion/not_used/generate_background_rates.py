"""
Generate Background Erosion Rates from DSAS Data for CASCADE Model

This script takes DSAS shoreline change data (LRR format) and produces
background erosion rates suitable for use in CASCADE.

Usage:
    python generate_background_rates.py

Inputs:
    - DSAS CSV file with domain_id and annual_rate_m_per_yr columns
    
Outputs:
    - Python file with BACKGROUND_EROSION_RATES list
    - CSV file with rates per domain
    - Visualization plot
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.ndimage import uniform_filter1d
from scipy.interpolate import interp1d
import sys
import os

# ============================================================
# CONFIGURATION - EDIT THESE VALUES
# ============================================================

# Input file
DSAS_CSV = r"C:\Users\hanna\PycharmProjects\CASCADE\data\hatteras_init\shoreline_change\dsas_1978_1997_CLEAN_DETAILED.csv"

# Output directory
OUTPUT_DIR = r"/output/background_rates"

# CASCADE domain structure
TOTAL_DOMAINS = 120  # Total including buffers
NUM_REAL_DOMAINS = 90  # Real island domains
NUM_LEFT_BUFFER = 15  # Number of left buffer domains
NUM_RIGHT_BUFFER = 15  # Number of right buffer domains

# IMPORTANT: Domain numbering
# Your DSAS file should contain only GIS domains 1-90 (your study area)
# These map to CASCADE domains 15-104 (with buffers on each side)
# GIS Domain 1 → CASCADE Domain 15 (south end)
# GIS Domain 90 → CASCADE Domain 104 (north end)

# Processing parameters
SMOOTHING_WINDOW = 1  # Larger = more smoothing (removes more storm signal)
SCALING_FACTOR = 0.9  # Reduce rates to partition background vs storm effects (0-1)
FORCE_MASS_BALANCE = True  # Force net rate to zero (closed sediment budget)

# Column names in input CSV
DOMAIN_ID_COL = "domain_id"
ANNUAL_RATE_COL = "annual_rate_m_per_yr"

# ============================================================
# PROCESSING FUNCTIONS
# ============================================================

def load_dsas_data(csv_path, domain_col, rate_col):
    """Load and validate DSAS data."""
    if not os.path.exists(csv_path):
        raise FileNotFoundError(f"Input file not found: {csv_path}")
    
    df = pd.read_csv(csv_path)
    
    # Check required columns
    if domain_col not in df.columns:
        raise ValueError(f"Column '{domain_col}' not found. Available: {df.columns.tolist()}")
    if rate_col not in df.columns:
        raise ValueError(f"Column '{rate_col}' not found. Available: {df.columns.tolist()}")
    
    print("="*70)
    print("LOADED DSAS DATA")
    print("="*70)
    print(f"File: {csv_path}")
    print(f"Records: {len(df)}")
    print(f"Domain range: {df[domain_col].min():.0f} to {df[domain_col].max():.0f}")
    print(f"Rate range: {df[rate_col].min():.3f} to {df[rate_col].max():.3f} m/yr")
    print(f"Mean rate: {df[rate_col].mean():.3f} m/yr")
    
    # Validate domain range
    max_domain = df[domain_col].max()
    if max_domain > 90:
        print(f"\n⚠️  WARNING: Found domains up to {max_domain:.0f}")
        print(f"   Expected only domains 1-90 for your study area.")
        print(f"   Domains 91+ may be from collaborator's area.")
        print(f"   Filtering to domains 1-90 only...")
        df = df[df[domain_col] <= 90].copy()
        print(f"   Filtered to {len(df)} records in domains 1-90")
    
    return df


def process_rates(df, domain_col, rate_col, num_real_domains, 
                  smoothing_window, scaling_factor, force_balance):
    """Process DSAS rates into background erosion rates."""
    
    # Extract data
    domain_ids = df[domain_col].values
    annual_rates = df[rate_col].values
    
    # Create interpolation function to fill any gaps
    interp_func = interp1d(
        domain_ids,
        annual_rates,
        kind='linear',
        bounds_error=False,
        fill_value='extrapolate'
    )
    
    # Create continuous domain array (1 to NUM_REAL_DOMAINS)
    real_domains = np.arange(1, num_real_domains + 1)
    continuous_rates = interp_func(real_domains)
    
    print("\n" + "="*70)
    print("PROCESSING STEPS")
    print("="*70)
    print(f"1. Interpolated to {len(continuous_rates)} continuous domains")
    
    # Smooth to extract background trend
    smoothed_rates = uniform_filter1d(continuous_rates, size=smoothing_window, mode='nearest')
    print(f"2. Applied smoothing (window size: {smoothing_window})")
    
    # Scale to partition background vs storm effects
    scaled_rates = smoothed_rates * scaling_factor
    print(f"3. Scaled by {scaling_factor} to partition background vs storms")
    
    # Check mass balance
    net_rate = np.sum(scaled_rates)
    print(f"4. Net island-wide rate: {net_rate:.3f} m/yr")
    
    # Force mass balance if requested
    if force_balance:
        adjustment = net_rate / num_real_domains
        balanced_rates = scaled_rates - adjustment
        print(f"5. Forced mass balance (adjusted by {-adjustment:.4f} m/yr per domain)")
        final_rates = balanced_rates
    else:
        print(f"5. Mass balance NOT enforced")
        final_rates = scaled_rates
    
    return final_rates, continuous_rates, smoothed_rates


def create_full_domain_array(processed_rates, total_domains, num_left_buffer, num_right_buffer):
    """Create full domain array with buffer zones."""
    
    background_erosion_rates = np.zeros(total_domains)
    
    # Set buffer zones to 0
    background_erosion_rates[:num_left_buffer] = 0.0
    background_erosion_rates[-num_right_buffer:] = 0.0
    
    # Fill in real island domains
    start_idx = num_left_buffer
    end_idx = num_left_buffer + len(processed_rates)
    background_erosion_rates[start_idx:end_idx] = processed_rates
    
    print("\n" + "="*70)
    print("FINAL ARRAY STRUCTURE")
    print("="*70)
    print(f"Total domains: {len(background_erosion_rates)}")
    print(f"  Left buffer (0.0): domains 0-{num_left_buffer-1}")
    print(f"  Real island: domains {num_left_buffer}-{end_idx-1}")
    print(f"  Right buffer (0.0): domains {end_idx}-{total_domains-1}")
    
    return background_erosion_rates, start_idx, end_idx


def save_outputs(background_rates, output_dir, smoothing_window, scaling_factor, force_balance):
    """Save background rates to files."""
    
    os.makedirs(output_dir, exist_ok=True)
    
    # Save as Python file
    py_file = os.path.join(output_dir, "background_erosion_rates.py")
    with open(py_file, 'w') as f:
        f.write("# Background erosion rates calculated from DSAS LRR data\n")
        f.write("# Units: m/yr\n")
        f.write("# Source: Linear Regression Rate (LRR) from DSAS\n")
        f.write(f"# Smoothing window: {smoothing_window} domains\n")
        f.write(f"# Scaling factor: {scaling_factor}\n")
        f.write(f"# Mass balanced: {force_balance}\n")
        f.write(f"# Net rate: {np.sum(background_rates):.4f} m/yr\n\n")
        
        f.write("BACKGROUND_EROSION_RATES = [\n")
        for i, rate in enumerate(background_rates):
            f.write(f"    {rate:7.4f},  # Domain {i}\n")
        f.write("]\n")
    
    print(f"\n✓ Saved Python file: {py_file}")
    
    # Save as CSV
    csv_file = os.path.join(output_dir, "background_erosion_rates.csv")
    df_out = pd.DataFrame({
        'domain_index': np.arange(len(background_rates)),
        'background_rate_m_per_yr': background_rates,
    })
    df_out.to_csv(csv_file, index=False)
    print(f"✓ Saved CSV file: {csv_file}")
    
    return py_file, csv_file


def create_visualization(df, domain_col, rate_col, continuous_rates, smoothed_rates, 
                        final_rates, background_rates, num_real_domains, 
                        num_left_buffer, output_dir):
    """Create diagnostic visualization."""
    
    fig, axes = plt.subplots(3, 1, figsize=(14, 10))
    
    # Plot 1: Original observed data
    ax1 = axes[0]
    real_domains = np.arange(1, num_real_domains + 1)
    ax1.plot(df[domain_col], df[rate_col], 'o-', 
             color='orange', linewidth=2, markersize=4, label='Observed LRR')
    ax1.axhline(y=0, color='gray', linestyle='--', alpha=0.5)
    ax1.set_xlabel('GIS Domain ID (1-90)')
    ax1.set_ylabel('Shoreline Change Rate (m/yr)')
    ax1.set_title('Original Observed Data (LRR from DSAS)')
    ax1.legend()
    ax1.grid(True, alpha=0.3)
    
    # Plot 2: Processing steps
    ax2 = axes[1]
    ax2.plot(real_domains, continuous_rates, 'o-', color='orange', 
             linewidth=1, markersize=3, alpha=0.5, label='Interpolated')
    ax2.plot(real_domains, smoothed_rates, '-', color='blue', 
             linewidth=2, label='Smoothed')
    ax2.plot(real_domains, final_rates, '-', color='green', 
             linewidth=2, label='Final (scaled & balanced)')
    ax2.axhline(y=0, color='gray', linestyle='--', alpha=0.5)
    ax2.set_xlabel('GIS Domain (1-90)')
    ax2.set_ylabel('Shoreline Change Rate (m/yr)')
    ax2.set_title('Processing Steps')
    ax2.legend()
    ax2.grid(True, alpha=0.3)
    
    # Plot 3: Final CASCADE array
    ax3 = axes[2]
    domain_indices = np.arange(len(background_rates))
    end_idx = num_left_buffer + num_real_domains
    colors = ['lightcoral' if i < num_left_buffer or i >= end_idx else 'steelblue' 
              for i in range(len(background_rates))]
    ax3.bar(domain_indices, background_rates, color=colors, edgecolor='black', linewidth=0.5)
    ax3.axhline(y=0, color='gray', linestyle='--', alpha=0.5)
    ax3.axvline(x=num_left_buffer-0.5, color='red', linestyle='--', alpha=0.7, label='Buffer boundaries')
    ax3.axvline(x=end_idx-0.5, color='red', linestyle='--', alpha=0.7)
    ax3.set_xlabel('CASCADE Domain Index (0-based)')
    ax3.set_ylabel('Background Erosion Rate (m/yr)')
    ax3.set_title('Final Background Rates for CASCADE (Buffers=red, Island=blue)')
    ax3.legend()
    ax3.grid(True, alpha=0.3, axis='y')
    
    plt.tight_layout()
    
    plot_file = os.path.join(output_dir, "background_rates_diagnostic.png")
    plt.savefig(plot_file, dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"✓ Saved visualization: {plot_file}")
    return plot_file


def print_statistics(background_rates, num_left_buffer, num_real_domains):
    """Print summary statistics."""
    
    end_idx = num_left_buffer + num_real_domains
    real_island_rates = background_rates[num_left_buffer:end_idx]
    
    print("\n" + "="*70)
    print("SUMMARY STATISTICS")
    print("="*70)
    print(f"Real island domains ({num_real_domains} total):")
    print(f"  Mean rate: {np.mean(real_island_rates):.4f} m/yr")
    print(f"  Std dev: {np.std(real_island_rates):.3f} m/yr")
    print(f"  Min rate: {np.min(real_island_rates):.3f} m/yr (domain {num_left_buffer + np.argmin(real_island_rates)})")
    print(f"  Max rate: {np.max(real_island_rates):.3f} m/yr (domain {num_left_buffer + np.argmax(real_island_rates)})")
    print(f"  Net change: {np.sum(real_island_rates):.4f} m/yr")
    
    # Identify erosion vs accretion zones
    erosion_domains = np.where(real_island_rates < -0.5)[0] + num_left_buffer
    accretion_domains = np.where(real_island_rates > 0.5)[0] + num_left_buffer
    
    print(f"\nErosion hotspots (< -0.5 m/yr): {len(erosion_domains)} domains")
    if len(erosion_domains) > 0:
        print(f"  Domains: {erosion_domains.tolist()}")
    
    print(f"\nAccretion hotspots (> 0.5 m/yr): {len(accretion_domains)} domains")
    if len(accretion_domains) > 0:
        print(f"  Domains: {accretion_domains.tolist()}")


# ============================================================
# MAIN EXECUTION
# ============================================================

def main():
    """Main execution function."""
    
    print("\n" + "="*70)
    print("BACKGROUND EROSION RATE GENERATOR")
    print("="*70)
    
    try:
        # Load data
        df = load_dsas_data(DSAS_CSV, DOMAIN_ID_COL, ANNUAL_RATE_COL)
        
        # Process rates
        final_rates, continuous_rates, smoothed_rates = process_rates(
            df, DOMAIN_ID_COL, ANNUAL_RATE_COL, NUM_REAL_DOMAINS,
            SMOOTHING_WINDOW, SCALING_FACTOR, FORCE_MASS_BALANCE
        )
        
        # Create full domain array
        background_rates, start_idx, end_idx = create_full_domain_array(
            final_rates, TOTAL_DOMAINS, NUM_LEFT_BUFFER, NUM_RIGHT_BUFFER
        )
        
        # Save outputs
        py_file, csv_file = save_outputs(
            background_rates, OUTPUT_DIR, SMOOTHING_WINDOW, 
            SCALING_FACTOR, FORCE_MASS_BALANCE
        )
        
        # Create visualization
        plot_file = create_visualization(
            df, DOMAIN_ID_COL, ANNUAL_RATE_COL, continuous_rates, 
            smoothed_rates, final_rates, background_rates, 
            NUM_REAL_DOMAINS, NUM_LEFT_BUFFER, OUTPUT_DIR
        )
        
        # Print statistics
        print_statistics(background_rates, NUM_LEFT_BUFFER, NUM_REAL_DOMAINS)
        
        print("\n" + "="*70)
        print("SUCCESS!")
        print("="*70)
        print(f"Generated files:")
        print(f"  1. {py_file}")
        print(f"  2. {csv_file}")
        print(f"  3. {plot_file}")
        print("\nNext steps:")
        print(f"  1. Open {os.path.basename(py_file)}")
        print(f"  2. Copy the BACKGROUND_EROSION_RATES list")
        print(f"  3. Paste into your CASCADE configuration file")
        print("="*70 + "\n")
        
    except Exception as e:
        print("\n" + "="*70)
        print("ERROR")
        print("="*70)
        print(f"{type(e).__name__}: {e}")
        print("\nPlease check:")
        print("  1. Input file path is correct")
        print("  2. Column names match your CSV file")
        print("  3. Output directory is writable")
        print("="*70 + "\n")
        sys.exit(1)


if __name__ == "__main__":
    main()
