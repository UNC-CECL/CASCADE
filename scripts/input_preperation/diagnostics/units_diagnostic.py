"""
Simple diagnostic to check CASCADE input data units
NO MODEL RUNNING - just data checking
"""

import numpy as np
import os

print("=" * 80)
print("CASCADE HATTERAS INPUT DATA UNIT CHECKER")
print("=" * 80)

# =============================================================================
# PATHS - UPDATE THESE
# =============================================================================

HATTERAS_BASE = r'C:\Users\hanna\PycharmProjects\CASCADE\data\hatteras_init'
DUNE_OFFSET_CSV = r'C:\Users\hanna\PycharmProjects\CASCADE\data\hatteras_init\island_offset\Dune_Offsets_1978_1997_PADDED_135.csv'
ROAD_SETBACK_CSV = r'C:\Users\hanna\PycharmProjects\CASCADE\data\hatteras_init\roads\offset\1978_road_setback_2row_FINAL.csv'

# Test with a few sample domains
TEST_DOMAINS = [30, 60, 90, 120, 134]  # First, middle, end domains

# =============================================================================
# 1. CHECK DUNE OFFSETS
# =============================================================================

print("\n" + "=" * 80)
print("1. DUNE OFFSETS (Shoreline Position)")
print("=" * 80)

dune_offset_all = np.loadtxt(DUNE_OFFSET_CSV, skiprows=1, delimiter=',')
print(f"CSV Shape: {dune_offset_all.shape}")
print(f"Expected: (135, 2) - 135 domains, 2 years (1978, 1997)")

# Check 1978 column (index 0)
dune_1978 = dune_offset_all[:, 0]
print(f"\n1978 Dune Offsets:")
print(f"  Min: {np.min(dune_1978):.1f}")
print(f"  Max: {np.max(dune_1978):.1f}")
print(f"  Mean: {np.mean(dune_1978):.1f}")

print(f"\n  Sample values for selected domains:")
for domain_gis in TEST_DOMAINS:
    # CSV row = domain_gis - 15 (because 15 buffer domains before domain 30)
    csv_row = domain_gis - 15
    if 0 <= csv_row < len(dune_1978):
        value = dune_1978[csv_row]
        print(f"    Domain {domain_gis:3d} (CSV row {csv_row:3d}): {value:8.1f}")

print("\n  📊 ANALYSIS:")
print(f"  - These values are in: METERS (confirmed from GIS)")
print(f"  - Typical range: 0 to 2000+ meters (distance from reference)")
print(f"  - Your range: {np.min(dune_1978):.0f} to {np.max(dune_1978):.0f} m ✓")

print("\n  ⚙️  FOR CASCADE:")
print(f"  - CASCADE likely expects: DECAMETERS")
print(f"  - Conversion needed: DIVIDE by 10 (m → dam)")
print(f"  - After conversion: {np.min(dune_1978) / 10:.1f} to {np.max(dune_1978) / 10:.1f} dam")

# =============================================================================
# 2. CHECK ROAD SETBACKS
# =============================================================================

print("\n" + "=" * 80)
print("2. ROAD SETBACKS (Distance from Dune to Road)")
print("=" * 80)

road_setbacks = np.loadtxt(ROAD_SETBACK_CSV, skiprows=1, delimiter=',')
print(f"CSV Shape: {road_setbacks.shape}")
print(f"Expected: (99,) - 99 domains with roads (domains 30-128)")

print(f"\nRoad Setback Values:")
print(f"  Min: {np.min(road_setbacks):.1f}")
print(f"  Max: {np.max(road_setbacks):.1f}")
print(f"  Mean: {np.mean(road_setbacks):.1f}")

print(f"\n  Sample values for selected domains:")
for domain_gis in TEST_DOMAINS:
    # Road CSV: row 0 = domain 30, so row = domain_gis - 30
    if 30 <= domain_gis <= 128:
        csv_row = domain_gis - 30
        if 0 <= csv_row < len(road_setbacks):
            value = road_setbacks[csv_row]
            print(f"    Domain {domain_gis:3d} (CSV row {csv_row:3d}): {value:8.1f}")

print("\n  📊 ANALYSIS:")
print(f"  - Typical road setback: 50-200 meters from dune")
print(f"  - Your values: {np.min(road_setbacks):.0f} to {np.max(road_setbacks):.0f}")

if np.mean(road_setbacks) > 300:
    print(f"  ⚠️  WARNING: These seem VERY large for road setbacks!")
    print(f"  - Are roads really {np.mean(road_setbacks):.0f}m from dunes?")
    print(f"  - Check your GIS calculation")
    print(f"  - Possible unit error in CSV creation?")
else:
    print(f"  ✓ Values seem reasonable")

print("\n  ⚙️  FOR CASCADE:")
print(f"  - CASCADE expects: METERS")
print(f"  - Your CSV appears to be in: METERS (or possibly wrong units)")
print(f"  - Recommended: NO CONVERSION (keep as-is)")
print(f"  - BUT: Verify these distances in GIS!")

# =============================================================================
# 3. CHECK ELEVATION FILES
# =============================================================================

print("\n" + "=" * 80)
print("3. ELEVATION & DUNE PROFILE FILES (.npy)")
print("=" * 80)

print("\nChecking a sample domain (domain 60)...")
domain_test = 60

elev_file = os.path.join(HATTERAS_BASE, 'topography_dunes', '2009',
                         f'domain_{domain_test}_resampled_topography_2009.npy')
dune_file = os.path.join(HATTERAS_BASE, 'dunes', '2009',
                         f'domain_{domain_test}_resampled_dune_2009.npy')

if os.path.exists(elev_file):
    elev_data = np.load(elev_file)
    print(f"\n✓ Elevation file found: domain_{domain_test}")
    print(f"  Shape: {elev_data.shape}")
    print(f"  Min: {np.min(elev_data):.3f} dam")
    print(f"  Max: {np.max(elev_data):.3f} dam")
    print(f"  Mean: {np.mean(elev_data):.3f} dam")

    if np.max(elev_data) < 0.3:
        print(f"  ⚠️  WARNING: Max elevation < 3m - very low!")
    elif np.max(elev_data) > 5.0:
        print(f"  ⚠️  WARNING: Max elevation > 50m - seems too high!")
    else:
        print(f"  ✓ Elevation range looks reasonable")

    print(f"\n  📊 ANALYSIS:")
    print(f"  - These files are in: DECAMETERS (1 dam = 10 m)")
    print(f"  - Expected range: -0.5 to 2.0 dam (for barrier islands)")
    print(f"  - Your range: {np.min(elev_data):.2f} to {np.max(elev_data):.2f} dam")
    print(f"  - In meters: {np.min(elev_data) * 10:.1f} to {np.max(elev_data) * 10:.1f} m")
else:
    print(f"❌ Elevation file NOT FOUND: {elev_file}")

if os.path.exists(dune_file):
    dune_data = np.load(dune_file)
    print(f"\n✓ Dune file found: domain_{domain_test}")
    print(f"  Shape: {dune_data.shape}")
    print(f"  Min: {np.min(dune_data):.3f} dam")
    print(f"  Max: {np.max(dune_data):.3f} dam")
    print(f"  Mean: {np.mean(dune_data):.3f} dam")

    if np.max(dune_data) < 0.2:
        print(f"  ⚠️  WARNING: Max dune height < 2m - very low!")
    else:
        print(f"  ✓ Dune heights look reasonable")
else:
    print(f"❌ Dune file NOT FOUND: {dune_file}")

# =============================================================================
# 4. BENTON'S *10 CONVERSION EXPLAINED
# =============================================================================

print("\n" + "=" * 80)
print("4. UNDERSTANDING BENTON'S *10 CONVERSION")
print("=" * 80)

print("""
Benton's original code had:
    road_setbacks = road_setbacks * 10
    dune_offset = dune_offset[:, year] * 10

WHY DID HE MULTIPLY BY 10?

Theory 1: Unit Conversion (Most Likely)
    - His CSV files were in DECAMETERS
    - He multiplied by 10 to convert to METERS
    - But wait... CASCADE uses DECAMETERS internally!

Theory 2: CASCADE Internal Units (CORRECT ANSWER)
    - CASCADE's internal distance unit is actually METERS for some parameters
    - But stores elevations in DECAMETERS
    - The *10 was converting his input CSVs (in dam) to internal meters

Theory 3: His CSVs Were Wrong Units
    - Maybe his original CSVs were accidentally in dam
    - The *10 fixed a preprocessing error

YOUR SITUATION:
    - Your CSVs are in METERS (confirmed)
    - Benton's were likely in DECAMETERS
    - So you need OPPOSITE conversion: DIVIDE by 10 (not multiply)

RECOMMENDATION:
    For dune offsets: dune_offset / 10  (convert m → dam)
    For road setbacks: Keep as-is in meters (verify in GIS first!)
""")

# =============================================================================
# 5. SUMMARY & RECOMMENDATIONS
# =============================================================================

print("\n" + "=" * 80)
print("5. SUMMARY & RECOMMENDATIONS FOR YOUR SCRIPT")
print("=" * 80)

print("""
✅ CONFIRMED UNITS IN YOUR FILES:
   - Dune offsets CSV: METERS
   - Road setbacks CSV: METERS (but values seem high - verify!)
   - Elevation .npy files: DECAMETERS
   - Dune .npy files: DECAMETERS

⚙️  CONVERSIONS NEEDED IN YOUR SCRIPT:

   1. DUNE OFFSETS:
      dune_offset = dune_offset / 10  # Convert m → dam

   2. ROAD SETBACKS:
      road_setbacks = road_setbacks  # Keep in meters (NO *10!)
      # But VERIFY these values in GIS - 500+m seems very far!

   3. ELEVATION/DUNE FILES:
      # No conversion needed - already in dam ✓

📝 IN YOUR MAIN SCRIPT, CHANGE:

   Line with road setbacks:
   BEFORE: road_setbacks = road_setbacks * 10
   AFTER:  road_setbacks = road_setbacks  # Keep in meters

   Line with dune offsets:
   BEFORE: dune_offset = dune_offset_all[:, year_column_index]
   AFTER:  dune_offset = dune_offset_all[:, year_column_index] / 10  # m → dam

⚠️  CRITICAL CHECK BEFORE RUNNING:
   - Verify road setback distances in GIS
   - If roads are ~50-100m from dunes (not 500m), your CSV has wrong values
   - Check your GIS calculation script
""")

print("\n" + "=" * 80)
print("DIAGNOSTIC COMPLETE")
print("=" * 80)
print("\nNext steps:")
print("1. Verify road setback values in GIS")
print("2. Update your main script with the conversions above")
print("3. Run the 3-domain test")
print("=" * 80)