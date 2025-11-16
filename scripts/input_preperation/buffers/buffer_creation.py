import os
import shutil
import sys


def copy_and_rename_files_with_absolute_paths():
    """
    Copies three specified NumPy files from their absolute source locations
    to the target buffer directory and renames them.

    *** IMPORTANT: Update the SOURCE_PATHS dictionary below with the actual
        full paths to your three source files. ***
    """

    # ----------------------------------------------------------------------
    # 1. Define Paths and Mappings (You MUST edit the file paths here)
    # ----------------------------------------------------------------------

    # Define the absolute path for the destination buffer folder
    DESTINATION_DIR = r'C:\Users\hanna\PycharmProjects\CASCADE\data\hatteras_init\buffer'

    # Define the source file paths and their corresponding target names
    # !! REPLACE 'PATH_TO_YOUR_FILE' with the actual full path on your system !!
    SOURCE_PATHS = {
        r'C:\Users\hanna\PycharmProjects\CASCADE\data\hatteras_init\dunes\2009\domain_111_resampled_dune_2009.npy': 'sample_1_dune.npy',
        r'C:\Users\hanna\PycharmProjects\CASCADE\data\hatteras_init\topography\2009\domain_111_resampled_topography_2009.npy': 'sample_1_topography.npy',
        r'C:\Users\hanna\PycharmProjects\CASCADE\data\hatteras_init\elevations\2009_pea_hatteras\domain_111_resampled.npy': 'domain_buffer.npy',
    }

    print(f"Target destination directory: {DESTINATION_DIR}")
    print("-" * 40)

    # ----------------------------------------------------------------------
    # 2. Ensure Destination Directory Exists
    # ----------------------------------------------------------------------
    try:
        os.makedirs(DESTINATION_DIR, exist_ok=True)
        print("Success: Destination buffer directory is ready.")
    except Exception as e:
        print(f"ERROR: Could not create or access directory {DESTINATION_DIR}")
        print(f"Reason: {e}")
        sys.exit(1)

    # ----------------------------------------------------------------------
    # 3. Perform Copy and Rename Operation
    # ----------------------------------------------------------------------
    successful_copies = 0

    for source_path, target_name in SOURCE_PATHS.items():
        destination_path = os.path.join(DESTINATION_DIR, target_name)

        # We use os.path.basename to get the original filename for logging
        source_filename = os.path.basename(source_path)

        try:
            # Check if source file actually exists before trying to copy
            if not os.path.exists(source_path):
                print(f"FAIL: Source file '{source_path}' not found. Skipping.")
                continue

            # Perform the binary file copy
            shutil.copyfile(source_path, destination_path)

            print(f"Copied: '{source_filename}' -> '{target_name}'")
            successful_copies += 1

        except Exception as e:
            print(f"ERROR copying '{source_filename}': {e}")

    print("-" * 40)

    if successful_copies == len(SOURCE_PATHS):
        print("✅ All files were copied and renamed successfully!")
    else:
        print(f"⚠️ Operation finished with {successful_copies}/{len(SOURCE_PATHS)} files copied.")


if __name__ == "__main__":
    copy_and_rename_files_with_absolute_paths()