import numpy as np
import os


# --- Helper function to check the first bytes of the file ---
def inspect_file_header(filepath):
    """Checks the first few bytes of a file to identify its format."""
    try:
        with open(filepath, 'rb') as f:
            header = f.read(10)

        print(f"\n--- Inspection for: {os.path.basename(filepath)} ---")
        print(f"File path: {filepath}")

        # NumPy magic string starts with \x93NUMPY
        if header.startswith(b'\x93NUMPY'):
            print("Status: ✅ Looks like a valid NumPy (.npy) file based on its header.")
            print("Action: You should be able to load this file normally.")
            return True

        # Check for UTF-8 BOM, which causes the '\xef' error
        elif header.startswith(b'\xef\xbb\xbf'):
            print("Status: ❌ File starts with a UTF-8 Byte Order Mark (BOM, the '\\xef' key).")
            print("Action: This file is likely a text file (CSV/TXT) saved with the wrong extension (.npy).")
            print("Recommended fix: Load it using a text-based function like np.loadtxt or pd.read_csv.")
            return False

        else:
            # Fallback for general corruption or unknown format
            print(f"Status: ❓ File format is unknown or corrupted. First 10 bytes: {header.hex()}")
            print("Action: Check the program that created this file. It might not have used np.save().")
            return False

    except FileNotFoundError:
        print(f"ERROR: File not found at {filepath}")
        return False
    except Exception as e:
        print(f"An unexpected error occurred during inspection: {e}")
        return False


# --- Main execution block to test the problematic file ---

# Replace this list with the actual 'd_file' list from your script
# We only have the path for the problematic file from the error message
problematic_file_path = r'/data/hatteras_init/buffer/sample_1_dune.npy'
d_file_test = [problematic_file_path]

# 1. Run the diagnostic check
is_npy = inspect_file_header(problematic_file_path)

# 2. Attempt to load the files based on the finding
if is_npy:
    # If the inspection passed, try to load it as a NumPy array again
    print("\nAttempting standard np.load() again...")
    try:
        initial_dune_profiles = np.concatenate([np.load(f, allow_pickle=True) for f in d_file_test], axis=0)
        print(f"SUCCESS: Loaded profiles with shape {initial_dune_profiles.shape}")
    except Exception as e:
        print(f"ERROR: Failed to load even after inspection check. Underlying error: {type(e).__name__}: {e}")

else:
    # If the inspection suggested it's a text file, try loading it as text data
    print("\nAttempting to load as a TEXT file (e.g., CSV/TXT) using np.loadtxt...")
    try:
        # Assuming the data is simple, space/comma-separated floats/integers
        initial_dune_profiles = np.concatenate([
            np.loadtxt(f, delimiter=',', ndmin=2) for f in d_file_test
        ], axis=0)
        print(f"SUCCESS: Loaded profiles with shape {initial_dune_profiles.shape} using np.loadtxt.")
        print("Note: You may need to adjust 'delimiter' if your file uses tabs or spaces.")
    except Exception as e:
        print(f"ERROR: Failed to load as text file. Underlying error: {type(e).__name__}: {e}")
        print(
            "Please manually inspect the content of 'sample_1_dune.npy' to determine the correct loading function (e.g., pandas.read_csv).")