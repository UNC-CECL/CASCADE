"""
Fix Corrupted CASCADE YAML Parameter File

Problem: The YAML file has numpy objects that newer PyYAML can't load
Solution: Clean the YAML file by converting all numpy types to native Python types

Run this script once to fix your parameter file.
"""

import yaml
import numpy as np
import os

# Path to your parameter file
PARAM_FILE = r'C:\Users\hanna\PycharmProjects\CASCADE\data\hatteras_init\Hatteras-CASCADE-parameters.yaml'
BACKUP_FILE = PARAM_FILE.replace('.yaml', '_BACKUP.yaml')


def numpy_to_python(obj):
    """
    Recursively convert numpy types to native Python types.
    """
    if isinstance(obj, dict):
        return {key: numpy_to_python(value) for key, value in obj.items()}
    elif isinstance(obj, list):
        return [numpy_to_python(item) for item in obj]
    elif isinstance(obj, np.integer):
        return int(obj)
    elif isinstance(obj, np.floating):
        return float(obj)
    elif isinstance(obj, np.ndarray):
        return obj.tolist()
    elif isinstance(obj, (np.bool_, bool)):
        return bool(obj)
    else:
        return obj


def fix_yaml_file(filepath):
    """
    Fix a corrupted CASCADE parameter YAML file.
    """
    print("=" * 80)
    print("FIXING CASCADE YAML PARAMETER FILE")
    print("=" * 80)
    print(f"File: {filepath}\n")

    # Check if file exists
    if not os.path.exists(filepath):
        print(f"❌ ERROR: File not found: {filepath}")
        return False

    # Create backup
    print("Creating backup...")
    try:
        with open(filepath, 'r') as f:
            original_content = f.read()
        with open(BACKUP_FILE, 'w') as f:
            f.write(original_content)
        print(f"✓ Backup saved: {BACKUP_FILE}")
    except Exception as e:
        print(f"❌ ERROR creating backup: {e}")
        return False

    # Try to load with unsafe loader first
    print("\nReading original file...")
    try:
        with open(filepath, 'r') as f:
            # Use unsafe loader to read numpy objects
            data = yaml.load(f, Loader=yaml.UnsafeLoader)
        print("✓ File loaded successfully")
    except Exception as e:
        print(f"❌ ERROR reading file: {e}")
        print("\nThe file may be too corrupted. You might need to recreate it from scratch.")
        return False

    # Convert all numpy types to Python types
    print("\nConverting numpy types to Python types...")
    cleaned_data = numpy_to_python(data)
    print("✓ Conversion complete")

    # Print what we're about to save
    print("\nPreview of cleaned parameters:")
    print("-" * 80)
    for key, value in list(cleaned_data.items())[:5]:
        print(f"  {key}: {value}")
    print("  ...")
    print("-" * 80)

    # Save cleaned version
    print("\nSaving cleaned file...")
    try:
        with open(filepath, 'w') as f:
            yaml.dump(cleaned_data, f, default_flow_style=False, sort_keys=False)
        print(f"✓ Cleaned file saved: {filepath}")
    except Exception as e:
        print(f"❌ ERROR saving file: {e}")
        # Restore backup
        print("Restoring backup...")
        with open(BACKUP_FILE, 'r') as f:
            backup_content = f.read()
        with open(filepath, 'w') as f:
            f.write(backup_content)
        print("✓ Backup restored")
        return False

    print("\n" + "=" * 80)
    print("SUCCESS! Parameter file fixed.")
    print("=" * 80)
    print("\nYou can now run your CASCADE simulation.")
    print(f"Original file backed up to: {BACKUP_FILE}")

    return True


# Alternative: Create a fresh parameter file from scratch
def create_default_parameters(filepath):
    """
    Create a fresh CASCADE parameter file with default Hatteras values.
    """

    default_params = {
        # Domain configuration
        'BarrierLength': 100,  # Alongshore length (dam)
        'DuneWidth': 20,  # Width of dune domain (dam)
        'LShoreface': 500,  # Shoreface length (dam)

        # Elevation parameters (dam = 10 meters)
        'BermEl': 0.144,  # Berm elevation (1.44 m MHW)
        'DShoreface': 1.0,  # Shoreface depth (10 m)
        'MHW': 0.0,  # Mean high water reference
        'SL': -0.244,  # Sea level relative to MHW

        # Dune parameters
        'rmin': 0.55,  # Minimum dune growth rate
        'rmax': 0.95,  # Maximum dune growth rate
        'Dmaxel': 0.5,  # Maximum dune elevation (5 m)
        'Dmin': 0.01,  # Minimum dune height

        # Overwash parameters
        'Cx': 50,  # Overwash decay constant
        'Qs_min': 1.0,  # Minimum overwash flux
        'MaxUpSlope': 0.25,  # Maximum slope for overwash deposition

        # Shoreface parameters
        'k_sf': 5000.0,  # Shoreface diffusivity
        's_sf_eq': 0.02,  # Equilibrium shoreface slope

        # Bay parameters
        'BayDepth': 0.3,  # Bay depth (3 m)

        # Subsidence
        'SubRate': 0.0,  # Subsidence rate (m/yr)

        # BRIE parameters (alongshore transport)
        'Rat': 1.0,  # Relative AST effect
        'Rin': 0.5,  # Inlet effect
    }

    print("=" * 80)
    print("CREATING DEFAULT CASCADE PARAMETER FILE")
    print("=" * 80)
    print(f"File: {filepath}\n")

    # Save
    try:
        with open(filepath, 'w') as f:
            yaml.dump(default_params, f, default_flow_style=False, sort_keys=False)
        print("✓ Default parameter file created successfully")
        print("\nPreview:")
        print("-" * 80)
        for key, value in list(default_params.items())[:10]:
            print(f"  {key}: {value}")
        print("-" * 80)
        return True
    except Exception as e:
        print(f"❌ ERROR creating file: {e}")
        return False


# ==============================================================================
# MAIN EXECUTION
# ==============================================================================

if __name__ == '__main__':
    print("\n" + "=" * 80)
    print("CASCADE YAML PARAMETER FILE FIXER")
    print("=" * 80)
    print("\nThis script will fix your corrupted YAML parameter file.")
    print("Choose an option:\n")
    print("1. Try to fix the existing file (recommended)")
    print("2. Create a fresh default parameter file")
    print("3. Exit")

    choice = input("\nEnter choice (1/2/3): ").strip()

    if choice == '1':
        success = fix_yaml_file(PARAM_FILE)
        if not success:
            print("\n⚠️  Fixing failed. You may want to try option 2 (create fresh file).")

    elif choice == '2':
        # Create backup first
        if os.path.exists(PARAM_FILE):
            print(f"\nBacking up existing file to {BACKUP_FILE}...")
            with open(PARAM_FILE, 'r') as f:
                with open(BACKUP_FILE, 'w') as bf:
                    bf.write(f.read())
            print("✓ Backup created")

        success = create_default_parameters(PARAM_FILE)
        if success:
            print("\n⚠️  NOTE: Review the default values and adjust as needed!")
            print("         Compare with your backup to ensure correct values.")

    elif choice == '3':
        print("\nExiting...")
    else:
        print("\n❌ Invalid choice")

    print()