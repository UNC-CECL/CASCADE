# =============================================================================
# CASCADE Standard Storm File Viewer
# -----------------------------------------------------------------------------
# Description:
#   Clean, standardized viewer for CASCADE storm files with known structure:
#   [Year_Index, Rhigh, Rlow, Wave Period, Duration]
#
# Author: Hannah Henry
# Created: 2025-01-12
# =============================================================================

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path
from typing import Optional

# =============================================================================
# USER CONFIGURATION
# =============================================================================

# Path to your storm file
STORM_PATH = Path(r"/data/hatteras_init/storms/hindcast_storms/storms_1978_1997.npy")

# Start year for converting Year_Index to Calendar_Year
START_YEAR = 1978

# Column names (standard for CASCADE storm files)
COLUMNS = ["Year_Index", "Rhigh", "Rlow", "Wave Period", "Duration"]

# Plotting options
SHOW_PLOTS = True
SAVE_PLOTS = False
OUTPUT_DIR = Path(r"/output/figures")

# Year verification (set to None to skip, or a year to check specific year)
# Useful for verifying against historical storm records
CHECK_YEAR = None  # e.g., 1985 to see all storms in 1985
# CHECK_YEAR = 1985  # Uncomment and set year to check specific year

# =============================================================================
# LOAD AND PROCESS DATA
# =============================================================================

def load_storm_file(path: Path) -> np.ndarray:
    """Load .npy storm file."""
    if not path.exists():
        raise FileNotFoundError(f"Storm file not found: {path}")
    
    arr = np.load(path)
    print(f"Loaded: {path}")
    print(f"Shape: {arr.shape} (rows={arr.shape[0]} storms, cols={arr.shape[1]})")
    return arr


def process_storm_data(arr: np.ndarray, columns: list, start_year: int) -> pd.DataFrame:
    """Convert array to DataFrame and add Calendar_Year."""
    df = pd.DataFrame(arr, columns=columns)
    
    # Convert Year_Index to Calendar_Year
    df['Calendar_Year'] = df['Year_Index'].astype(int) + start_year
    df = df.drop(columns=['Year_Index'])
    
    # Reorder columns to put Calendar_Year first
    cols = ['Calendar_Year'] + [c for c in df.columns if c != 'Calendar_Year']
    df = df[cols]
    
    print(f"Converted Year_Index to Calendar_Year (starting {start_year})")
    
    return df


# =============================================================================
# SUMMARY STATISTICS
# =============================================================================

def check_year(df: pd.DataFrame, year: int):
    """
    Print detailed storm information for a specific year.
    Useful for verification against historical storm records.
    
    Args:
        df: Storm DataFrame
        year: Calendar year to check
    """
    storms = df[df['Calendar_Year'] == year].copy()
    
    print("\n" + "="*80)
    print(f"STORMS IN {year}")
    print("="*80)
    
    if len(storms) == 0:
        print(f"\nNo storms recorded in {year}")
        print("="*80 + "\n")
        return
    
    print(f"\nTotal storms: {len(storms)}")
    print(f"\nRhigh range: {storms['Rhigh'].min():.3f} to {storms['Rhigh'].max():.3f}")
    print(f"Average Rhigh: {storms['Rhigh'].mean():.3f}")
    print(f"Total storm duration: {storms['Duration'].sum():.1f} hours")
    
    # All storms that year
    print("\n" + "-"*80)
    print(f"ALL {len(storms)} STORMS IN {year}:")
    print("-"*80)
    
    # Sort by Rhigh descending
    storms_sorted = storms.sort_values('Rhigh', ascending=False).reset_index(drop=True)
    storms_sorted.index = storms_sorted.index + 1  # Start index at 1 instead of 0
    
    print(storms_sorted[['Calendar_Year', 'Rhigh', 'Rlow', 'Wave Period', 'Duration']].to_string())
    
    # Highlight most intense
    max_storm = storms_sorted.iloc[0]
    print(f"\nMost intense storm: Rhigh={max_storm['Rhigh']:.3f}, Duration={max_storm['Duration']:.1f} hrs")
    
    print("="*80 + "\n")


def print_summary(df: pd.DataFrame):
    """Print comprehensive summary statistics."""
    print("\n" + "="*80)
    print("STORM FILE SUMMARY")
    print("="*80)
    
    # Basic info
    year_min = int(df['Calendar_Year'].min())
    year_max = int(df['Calendar_Year'].max())
    n_years = year_max - year_min + 1
    
    print(f"\nTotal Events: {len(df)}")
    print(f"Time Period: {year_min}–{year_max} ({n_years} years)")
    
    # Storms per year
    storms_per_year = df.groupby('Calendar_Year').size()
    print(f"Average storms per year: {storms_per_year.mean():.1f}")
    print(f"Min/Max storms per year: {storms_per_year.min():.0f} / {storms_per_year.max():.0f}")
    
    # First 10 events
    print("\n" + "-"*80)
    print("FIRST 10 EVENTS:")
    print("-"*80)
    print(df.head(10).to_string(index=False))
    
    # Statistical summary
    print("\n" + "-"*80)
    print("STATISTICAL SUMMARY:")
    print("-"*80)
    summary = df[["Rhigh", "Rlow", "Wave Period", "Duration"]].describe(
        percentiles=[0.01, 0.25, 0.5, 0.75, 0.95, 0.99]
    )
    print(summary.to_string())
    
    # Top 5 most intense storms
    print("\n" + "-"*80)
    print("TOP 5 MOST INTENSE STORMS (by Rhigh):")
    print("-"*80)
    top_storms = df.nlargest(5, 'Rhigh')
    print(top_storms.to_string(index=False))
    
    # Storms per year (first 10 years)
    print("\n" + "-"*80)
    print("STORMS PER YEAR (first 10 years):")
    print("-"*80)
    yearly_counts = storms_per_year.head(10)
    for year, count in yearly_counts.items():
        print(f"  {int(year)}: {count} storms")
    
    print("\n" + "="*80 + "\n")


# =============================================================================
# VISUALIZATION
# =============================================================================

def create_distribution_plots(df: pd.DataFrame, save_dir: Optional[Path] = None):
    """Create four-panel distribution plot."""
    
    fig, axs = plt.subplots(2, 2, figsize=(14, 10))
    
    year_min = int(df['Calendar_Year'].min())
    year_max = int(df['Calendar_Year'].max())
    fig.suptitle(f"Storm Parameter Distributions (N={len(df)} Events, {year_min}–{year_max})", 
                 fontsize=16)
    
    # Panel 1: Rhigh Distribution
    df["Rhigh"].hist(ax=axs[0,0], bins=20, color="#1f77b4", edgecolor='black', alpha=0.8)
    axs[0,0].set_title("1. Storm Intensity (Rhigh) Distribution")
    axs[0,0].set_xlabel("Rhigh Value")
    axs[0,0].set_ylabel("Count")
    axs[0,0].grid(alpha=0.3)
    
    # Panel 2: Wave Period Distribution
    df["Wave Period"].hist(ax=axs[0,1], bins=20, color="#ff7f0e", edgecolor='black', alpha=0.8)
    axs[0,1].set_title("2. Wave Period Distribution")
    axs[0,1].set_xlabel("Wave Period (s)")
    axs[0,1].set_ylabel("Count")
    axs[0,1].grid(alpha=0.3)
    
    # Panel 3: Duration Distribution
    df["Duration"].hist(ax=axs[1,0], bins=20, color="#2ca02c", edgecolor='black', alpha=0.8)
    axs[1,0].set_title("3. Storm Duration Distribution")
    axs[1,0].set_xlabel("Duration (hours)")
    axs[1,0].set_ylabel("Count")
    axs[1,0].grid(alpha=0.3)
    
    # Panel 4: Rlow vs Rhigh Scatter
    axs[1,1].scatter(df["Rlow"], df["Rhigh"], color="#d62728", alpha=0.7, s=30)
    axs[1,1].set_title("4. Rlow vs Rhigh Intensity Comparison")
    axs[1,1].set_xlabel("Rlow Value")
    axs[1,1].set_ylabel("Rhigh Value")
    axs[1,1].grid(True, linestyle='--', alpha=0.6)
    
    plt.tight_layout(rect=[0, 0, 1, 0.96])
    
    if save_dir and SAVE_PLOTS:
        save_dir.mkdir(parents=True, exist_ok=True)
        filepath = save_dir / "storm_distributions.png"
        fig.savefig(filepath, dpi=300, bbox_inches='tight')
        print(f"Saved: {filepath}")
    
    if SHOW_PLOTS:
        plt.show()
    else:
        plt.close()


def create_timeseries_plot(df: pd.DataFrame, save_dir: Optional[Path] = None):
    """Create time series plot of storm intensity."""
    
    fig, ax = plt.subplots(figsize=(15, 6))
    
    # Individual storms (smaller markers)
    ax.plot(df['Calendar_Year'], df['Rhigh'], 'o', 
           label='Individual Storm Rhigh', alpha=0.5, color='gray', markersize=4)
    
    # Yearly average (prominent line)
    yearly_avg = df.groupby('Calendar_Year')['Rhigh'].mean().reset_index()
    ax.plot(yearly_avg['Calendar_Year'], yearly_avg['Rhigh'], 
           '-o', linewidth=2.5, label='Average Rhigh per Year', 
           color='#9467bd', markersize=6)
    
    year_min = int(df['Calendar_Year'].min())
    year_max = int(df['Calendar_Year'].max())
    ax.set_title(f"Storm Intensity (Rhigh) Over Time ({year_min}–{year_max})", fontsize=14)
    ax.set_xlabel("Calendar Year", fontsize=12)
    ax.set_ylabel("Rhigh Value (Storm Intensity)", fontsize=12)
    ax.legend(loc='best', fontsize=10)
    ax.grid(True, linestyle='--', alpha=0.7)
    
    plt.tight_layout()
    
    if save_dir and SAVE_PLOTS:
        save_dir.mkdir(parents=True, exist_ok=True)
        filepath = save_dir / "storm_timeseries.png"
        fig.savefig(filepath, dpi=300, bbox_inches='tight')
        print(f"Saved: {filepath}")
    
    if SHOW_PLOTS:
        plt.show()
    else:
        plt.close()


# =============================================================================
# MAIN EXECUTION
# =============================================================================

def main():
    """Main execution function."""
    print("\n" + "="*80)
    print("CASCADE STANDARD STORM FILE VIEWER")
    print("="*80 + "\n")
    
    # Load data
    arr = load_storm_file(STORM_PATH)
    
    # Process into DataFrame
    df = process_storm_data(arr, COLUMNS, START_YEAR)
    
    # Print summary
    print_summary(df)
    
    # Check specific year if requested
    if CHECK_YEAR is not None:
        check_year(df, CHECK_YEAR)
    
    # Create visualizations
    if SHOW_PLOTS or SAVE_PLOTS:
        print("Creating visualizations...")
        create_distribution_plots(df, OUTPUT_DIR if SAVE_PLOTS else None)
        create_timeseries_plot(df, OUTPUT_DIR if SAVE_PLOTS else None)
        print("✓ Plots complete")
    
    print("\n✓ Analysis complete!\n")
    
    return df


if __name__ == "__main__":
    df = main()
    
    # Interactive year checking
    # Uncomment the lines below to interactively check specific years after running
    # check_year(df, 1985)  # Check 1985
    # check_year(df, 1996)  # Check 1996
    # check_year(df, 1999)  # Check 1999
