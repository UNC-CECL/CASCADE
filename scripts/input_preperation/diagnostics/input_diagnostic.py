"""
Diagnostic script to check CASCADE input data for Hatteras
Run this BEFORE running your main simulation to identify issues
"""

import numpy as np
import os
import matplotlib.pyplot as plt

# =============================================================================
# CONFIGURATION - Match your main script
# =============================================================================

HATTERAS_BASE = r'C:\Users\hanna\PycharmProjects\CASCADE\data\hatteras_init'
DUNE_LOAD_ALL = r'C:\Users\hanna\PycharmProjects\CASCADE\data\hatteras_init\island_offset\Dune_Offsets_1978_1997_PADDED_135.csv'
ROAD_LOAD_1978 = r'C:\Users\hanna\PycharmProjects\CASCADE\data\hatteras_init\roads\offset\1978_road_setback_2row_FINAL.csv'

START_FILE_NUMBER = 30
START_REAL_DOMAIN_INDEX = 15
END_REAL_DOMAIN_INDEX = 120

start_year = 1978
year_column_index = 0  # 1978 is first column

# =============================================================================
# CHECK 1: Dune Offsets
# =============================================================================

print("=" * 80)
print("CHECK 1: DUNE OFFSETS")
print("=" * 80)

dune_offset_all = np.loadtxt(DUNE_LOAD_ALL, skiprows=1, delimiter=',')
dune_offset = dune_offset_all[:, year_column_index]

print(f"Dune offset shape: {dune_offset.shape}")
print(f"Dune offset range: [{np.min(dune_offset):.2f}, {np.max(dune_offset):.2f}]")
print(f"Dune offset mean: {np.mean(dune_offset):.2f}")
print(f"Dune offset units: (Assuming METERS - verify this!)")
print(f"\nFirst 10 values: {dune_offset[:10]}")
print(f"Last 10 values: {dune_offset[-10:]}")

# Check for anomalies
if np.any(dune_offset > 1000):
    print("\n⚠️ WARNING: Dune offsets > 1000m detected - possible unit error!")
if np.any(dune_offset < -1000):
    print("\n⚠️ WARNING: Dune offsets < -1000m detected - possible unit error!")

# =============================================================================
# CHECK 2: Road Setbacks
# =============================================================================

print("\n" + "=" * 80)
print("CHECK 2: ROAD SETBACKS")
print("=" * 80)

road_setbacks = np.loadtxt(ROAD_LOAD_1978, skiprows=1, delimiter=',')

print(f"Road setbacks shape: {road_setbacks.shape}")
print(f"Road setbacks range (before *10): [{np.min(road_setbacks):.2f}, {np.max(road_setbacks):.2f}]")
print(f"Road setbacks range (after *10): [{np.min(road_setbacks * 10):.2f}, {np.max(road_setbacks * 10):.2f}]")
print(f"\n⚠️ CHECK: Are your road setbacks in METERS or DECAMETERS?")
print(f"   If in meters, remove the '*10' multiplication!")

# =============================================================================
# CHECK 3: Sample Elevation and Dune Files
# =============================================================================

print("\n" + "=" * 80)
print("CHECK 3: SAMPLE ELEVATION FILES")
print("=" * 80)

# Check a few sample domains
sample_domains = [30, 50, 100, 134]  # First, middle, last

for domain_num in sample_domains:
    try:
        # Calculate list index
        if domain_num < START_FILE_NUMBER:
            print(f"\nDomain {domain_num}: BUFFER (using sample_1)")
            elev_file = os.path.join(HATTERAS_BASE, 'buffer', 'sample_1_topography.npy')
            dune_file = os.path.join(HATTERAS_BASE, 'buffer', 'sample_1_dune.npy')
        elif domain_num >= START_FILE_NUMBER and domain_num < START_FILE_NUMBER + 105:
            print(f"\nDomain {domain_num}: REAL DOMAIN")
            elev_file = os.path.join(HATTERAS_BASE, 'topography_dunes', '2009',
                                     f'domain_{domain_num}_resampled_topography_2009.npy')
            dune_file = os.path.join(HATTERAS_BASE, 'dunes', '2009',
                                     f'domain_{domain_num}_resampled_dune_2009.npy')

        # Load elevation data
        if os.path.exists(elev_file):
            elev_data = np.load(elev_file)
            print(f"   Elevation file: EXISTS")
            print(f"   Shape: {elev_data.shape}")
            print(f"   Range: [{np.min(elev_data):.3f}, {np.max(elev_data):.3f}]")
            print(f"   Mean: {np.mean(elev_data):.3f}")
            print(f"   Units: (Should be in DECAMETERS for CASCADE)")

            # Critical check for drowning
            if np.max(elev_data) < 0.3:  # Less than 3m in decameters
                print(f"   ⚠️ WARNING: Maximum elevation is very low! May cause drowning.")
        else:
            print(f"   ❌ ERROR: Elevation file NOT FOUND: {elev_file}")

        # Load dune data
        if os.path.exists(dune_file):
            dune_data = np.load(dune_file)
            print(f"   Dune file: EXISTS")
            print(f"   Shape: {dune_data.shape}")
            print(f"   Range: [{np.min(dune_data):.3f}, {np.max(dune_data):.3f}]")
            print(f"   Mean: {np.mean(dune_data):.3f}")

            if np.max(dune_data) < 0.3:
                print(f"   ⚠️ WARNING: Maximum dune height is very low! May cause drowning.")
        else:
            print(f"   ❌ ERROR: Dune file NOT FOUND: {dune_file}")

    except Exception as e:
        print(f"   ❌ ERROR loading domain {domain_num}: {e}")

# =============================================================================
# CHECK 4: CASCADE Parameter File
# =============================================================================

print("\n" + "=" * 80)
print("CHECK 4: CASCADE PARAMETERS")
print("=" * 80)

param_file = os.path.join(HATTERAS_BASE, "Hatteras-CASCADE-parameters.yaml")
if os.path.exists(param_file):
    print(f"Parameter file EXISTS: {param_file}")
    print("Check these parameters in the YAML file:")
    print("  - BermEl (Berm Elevation)")
    print("  - MHW (Mean High Water)")
    print("  - DShoreface (Shoreface Depth)")
    print("  - Initial sea level settings")
else:
    print(f"❌ ERROR: Parameter file NOT FOUND: {param_file}")

# =============================================================================
# CHECK 5: Unit Consistency Check
# =============================================================================

print("\n" + "=" * 80)
print("CHECK 5: UNIT CONSISTENCY ANALYSIS")
print("=" * 80)

print("\nCASCADE expects:")
print("  - Elevations/PEA_2011: DECAMETERS (1 dam = 10 meters)")
print("  - Distances: METERS")
print("  - Dune offsets: METERS")
print("  - Road setbacks: METERS (then converted to dam internally?)")
print("\nCommon mistakes:")
print("  1. Elevation data in METERS instead of DECAMETERS")
print("  2. Dune offsets in wrong units")
print("  3. Multiplying already-correct values by 10")

# =============================================================================
# RECOMMENDATIONS
# =============================================================================

print("\n" + "=" * 80)
print("RECOMMENDATIONS TO FIX DROWNING")
print("=" * 80)

print("\n1. CHECK ELEVATION UNITS:")
print("   - If your .npy files are in METERS, divide by 10 when loading")
print("   - If already in DECAMETERS, load as-is")

print("\n2. CHECK DUNE OFFSET UNITS:")
print("   - Dune offsets should typically be in METERS")
print("   - Check if values make physical sense (usually < 500m)")

print("\n3. CHECK ROAD SETBACK MULTIPLICATION:")
print("   - Remove '*10' if road_setbacks are already in correct units")
print("   - Or verify CASCADE expects decameters for this parameter")

print("\n4. CHECK SEA LEVEL RISE SETTINGS:")
print("   - sea_level_rise_rate = 0.0056 m/yr")
print("   - Initial sea level should match your datum")

print("\n5. VERIFY BERM ELEVATION:")
print("   - Check BermEl in Hatteras-CASCADE-parameters.yaml")
print("   - Should typically be 0.5-1.5 dam (5-15m) above MHW")

print("\n6. CHECK INITIAL TOPOGRAPHY:")
print("   - Run this on ONE domain first to debug")
print("   - Visualize the cross-shore profile")

print("\n" + "=" * 80)
print("DIAGNOSTIC COMPLETE")
print("=" * 80)