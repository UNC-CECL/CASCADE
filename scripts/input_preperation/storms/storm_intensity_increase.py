"""
Storm Intensity Modifier for CASCADE

Modifies existing storm series to test if storm intensity is limiting overwash.
Creates enhanced storm scenarios while preserving storm timing and structure.

Storm Series Format:
  Column 0: Year_Index (integer, 0-based)
  Column 1: Rhigh (dam) - Total water level (surge + runup)
  Column 2: Rlow (dam) - Minimum water level during storm
  Column 3: Wave_Period (seconds)
  Column 4: Duration (hours)
"""

import numpy as np
import os
from pathlib import Path
import matplotlib.pyplot as plt

# ==========================
# CONFIG
# ==========================

# Input storm file
ORIGINAL_STORM_FILE = Path(
    r"C:\Users\hanna\PycharmProjects\CASCADE\data\hatteras_init\storms"
    r"\hindcast_storms\storms_1978_1997.npy"
)

# Output directory for modified storms
OUTPUT_DIR = Path(
    r"C:\Users\hanna\PycharmProjects\CASCADE\data\hatteras_init\storms"
    r"\hindcast_storms\modified"
)

OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# Modification scenarios to test
SCENARIOS = {
    "1.2x_intensity": {
        "rhigh_multiplier": 1.2,
        "description": "20% increase in storm intensity (realistic enhancement)"
    },
    "1.5x_intensity": {
        "rhigh_multiplier": 1.5,
        "description": "50% increase in storm intensity (moderate test)"
    },
    "2.0x_intensity": {
        "rhigh_multiplier": 2.0,
        "description": "100% increase in storm intensity (extreme test)"
    },
    "add_0.15dam": {
        "rhigh_additive": 0.15,  # Add 1.5m to all storms
        "description": "Add 1.5m to all storm surge heights"
    },
}


# ==========================
# FUNCTIONS
# ==========================

def load_storm_series(filepath):
    """Load and verify storm series."""
    storms = np.load(filepath)

    print(f"\nLoaded storm series from: {filepath.name}")
    print(f"  Shape: {storms.shape}")
    print(f"  Number of storms: {storms.shape[0]}")

    if storms.shape[1] != 5:
        print(f"  ⚠️  WARNING: Expected 5 columns, got {storms.shape[1]}")

    return storms


def analyze_storms(storms, label="Original"):
    """Print statistics about storm series."""
    print(f"\n{'=' * 70}")
    print(f"{label} Storm Statistics")
    print(f"{'=' * 70}")

    rhigh_dam = storms[:, 1]
    rhigh_m = rhigh_dam * 10

    print(f"\nRhigh (Total Water Level):")
    print(f"  Min:    {rhigh_dam.min():.3f} dam = {rhigh_m.min():.2f}m")
    print(f"  Max:    {rhigh_dam.max():.3f} dam = {rhigh_m.max():.2f}m")
    print(f"  Mean:   {rhigh_dam.mean():.3f} dam = {rhigh_m.mean():.2f}m")
    print(f"  Median: {np.median(rhigh_dam):.3f} dam = {np.median(rhigh_m):.2f}m")

    # Categorize storms
    weak = np.sum(rhigh_dam < 0.25)
    moderate = np.sum((rhigh_dam >= 0.25) & (rhigh_dam < 0.35))
    strong = np.sum((rhigh_dam >= 0.35) & (rhigh_dam < 0.45))
    extreme = np.sum(rhigh_dam >= 0.45)

    print(f"\nStorm Categories:")
    print(f"  Weak     (Rhigh < 2.5m):  {weak} storms")
    print(f"  Moderate (2.5-3.5m):      {moderate} storms")
    print(f"  Strong   (3.5-4.5m):      {strong} storms")
    print(f"  Extreme  (> 4.5m):        {extreme} storms")

    # Check for potential overwash (need Rhigh > ~4m for 4-5m dunes)
    overwash_potential = np.sum(rhigh_m > 4.0)
    print(f"\nStorms likely to cause overwash (Rhigh > 4m): {overwash_potential}")

    if overwash_potential == 0:
        print("  ⚠️  No storms strong enough to overtop 4-5m engineered dunes!")

    return {
        'rhigh_dam': rhigh_dam,
        'rhigh_m': rhigh_m,
        'weak': weak,
        'moderate': moderate,
        'strong': strong,
        'extreme': extreme,
        'overwash_potential': overwash_potential
    }


def modify_storm_intensity(storms, rhigh_multiplier=None, rhigh_additive=None):
    """
    Modify storm intensity while preserving structure.

    Parameters:
    -----------
    storms : ndarray
        Original storm series
    rhigh_multiplier : float, optional
        Multiply Rhigh by this factor (e.g., 1.2 = 20% increase)
    rhigh_additive : float, optional
        Add this value (in dam) to all Rhigh values

    Returns:
    --------
    modified_storms : ndarray
        Modified storm series
    """
    modified = storms.copy()

    if rhigh_multiplier is not None:
        modified[:, 1] = storms[:, 1] * rhigh_multiplier
        print(f"  Applied {rhigh_multiplier}x multiplier to Rhigh")

    if rhigh_additive is not None:
        modified[:, 1] = storms[:, 1] + rhigh_additive
        print(f"  Added {rhigh_additive} dam ({rhigh_additive * 10}m) to Rhigh")

    # Also adjust Rlow proportionally to maintain storm structure
    if rhigh_multiplier is not None:
        modified[:, 2] = storms[:, 2] * rhigh_multiplier
    elif rhigh_additive is not None:
        modified[:, 2] = storms[:, 2] + rhigh_additive * 0.5  # Less adjustment to Rlow

    # Ensure Rlow doesn't exceed Rhigh (physical constraint)
    modified[:, 2] = np.minimum(modified[:, 2], modified[:, 1] * 0.9)

    return modified


def plot_storm_comparison(original, modified, scenario_name, save_dir):
    """Plot comparison of original vs modified storms."""
    fig, axes = plt.subplots(2, 1, figsize=(12, 8))

    years_orig = original[:, 0]
    rhigh_orig = original[:, 1] * 10  # Convert to meters

    years_mod = modified[:, 0]
    rhigh_mod = modified[:, 1] * 10

    # Plot 1: Time series comparison
    ax1 = axes[0]
    ax1.scatter(years_orig, rhigh_orig, alpha=0.6, s=30, label='Original', color='blue')
    ax1.scatter(years_mod, rhigh_mod, alpha=0.6, s=30, label='Modified', color='red', marker='^')
    ax1.axhline(4.0, color='green', linestyle='--', linewidth=2,
                label='~4m (needed for overwash)', alpha=0.7)
    ax1.set_xlabel('Year Index')
    ax1.set_ylabel('Storm Rhigh (m)')
    ax1.set_title(f'Storm Intensity Comparison: {scenario_name}')
    ax1.legend()
    ax1.grid(alpha=0.3)

    # Plot 2: Histogram comparison
    ax2 = axes[1]
    bins = np.linspace(0, max(rhigh_orig.max(), rhigh_mod.max()) + 0.5, 30)
    ax2.hist(rhigh_orig, bins=bins, alpha=0.5, label='Original', color='blue')
    ax2.hist(rhigh_mod, bins=bins, alpha=0.5, label='Modified', color='red')
    ax2.axvline(4.0, color='green', linestyle='--', linewidth=2,
                label='~4m threshold', alpha=0.7)
    ax2.set_xlabel('Storm Rhigh (m)')
    ax2.set_ylabel('Count')
    ax2.set_title('Distribution of Storm Intensities')
    ax2.legend()
    ax2.grid(alpha=0.3)

    plt.tight_layout()

    plot_path = save_dir / f'comparison_{scenario_name}.png'
    plt.savefig(plot_path, dpi=150)
    plt.close()

    print(f"  Saved comparison plot: {plot_path.name}")


def main():
    """Main execution function."""
    print("=" * 70)
    print("CASCADE STORM INTENSITY MODIFIER")
    print("=" * 70)

    # Load original storms
    original_storms = load_storm_series(ORIGINAL_STORM_FILE)
    original_stats = analyze_storms(original_storms, "Original")

    # Process each scenario
    for scenario_name, scenario_config in SCENARIOS.items():
        print(f"\n{'=' * 70}")
        print(f"Processing Scenario: {scenario_name}")
        print(f"Description: {scenario_config['description']}")
        print(f"{'=' * 70}")

        # Modify storms
        modified_storms = modify_storm_intensity(
            original_storms,
            rhigh_multiplier=scenario_config.get('rhigh_multiplier'),
            rhigh_additive=scenario_config.get('rhigh_additive')
        )

        # Analyze modified storms
        modified_stats = analyze_storms(modified_storms, f"Modified ({scenario_name})")

        # Save modified storm series
        output_file = OUTPUT_DIR / f"storms_1978_1997_{scenario_name}.npy"
        np.save(output_file, modified_storms)
        print(f"\n✓ Saved: {output_file}")

        # Create comparison plot
        plot_storm_comparison(original_storms, modified_storms, scenario_name, OUTPUT_DIR)

    # Summary
    print("\n" + "=" * 70)
    print("MODIFICATION COMPLETE")
    print("=" * 70)
    print(f"\nOutput directory: {OUTPUT_DIR}")
    print(f"Generated {len(SCENARIOS)} modified storm scenarios")

    print("\n" + "=" * 70)
    print("NEXT STEPS:")
    print("=" * 70)
    print("1. Review the comparison plots to verify modifications")
    print("2. Choose a scenario to test (recommend starting with 1.2x_intensity)")
    print("3. Update your CASCADE run script to use the modified storm file:")
    print(f"\n   STORM_FILE = r'{OUTPUT_DIR / 'storms_1978_1997_1.2x_intensity.npy'}'")
    print("\n4. Re-run CASCADE and check for overwash in the output")
    print("5. If 1.2x doesn't produce overwash, try 1.5x or 2.0x")
    print("=" * 70)


if __name__ == "__main__":
    main()