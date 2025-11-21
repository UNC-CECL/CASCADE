# ==============================================================================
# Hatteras CASCADE Dune and Topography Extractor
# ==============================================================================
# Author: Hannah Henry
# Last modified: 2025-08-1
# Description:
# This script calculates the average distance from the dune toe to a baseline
# transect for each barrier island domain and outputs both absolute and relative
# offsets for use in CASCADE model hindcasts.
# ==============================================================================

"""
Purpose:
---------
This script processes .npy elevation arrays (exported from ArcGIS Pro) for CASCADE hindcasting.
It extracts dune crest elevations and island interior elevation profiles, relative to Mean High Water (MHW),
and saves them in decameters as required by CASCADE.

Site: Hatteras Island, NC
Hindcast Year: 1978
Elevation Data Source: 2009 LiDAR-derived DEMs (used as best-available approximation for 1978)

MHW Notes:
----------
All elevation arrays are referenced to NAVD88 in **meters**, so MHW must also be in meters relative to NAVD88.

Berm elevation is set to 1.9 meters, a standard CASCADE value for NC barrier islands.

Output:
--------
- Topography array (landward of dune), shape: (100, n_columns), values in decameters
- Dune height vector, shape: (n_columns,), values in decameters
"""

import os
import numpy as np

# === INPUT PATHS AND PARAMETERS ============================================

LOAD_PATH = r'/data/hatteras_init/elevations/old_study_area_2009'
TOPO_SAVE_PATH = r'/data/hatteras_init/topography_dunes'
DUNE_SAVE_PATH = r'/data/hatteras_init/dunes'

os.makedirs(TOPO_SAVE_PATH, exist_ok=True)
os.makedirs(DUNE_SAVE_PATH, exist_ok=True)

# === KEY PARAMETERS ========================================================

MHW = 0.26       # meters, relative to NAVD88
# Previously used value in legacy scripts: MHW = 0.26
BERM_ELEV = 1.9    # meters, CASCADE standard berm elevation for NC

SEAWARD_EDGE = 'Right'  # 'Top', 'Bottom', 'Left', or 'Right' depending on array orientation. Confirmed by visual inspection: lowest elevations are on the right

# === MAIN FUNCTION =========================================================

def process_topo_arrays(input_path, topo_out_path, dune_out_path,
                        mhw=0.26, berm_elev=1.9, seaward_edge='Right'):
    """
    Processes elevation arrays to extract dune and island interior elevations for CASCADE.

    Parameters:
        input_path (str): Folder with input .npy elevation arrays
        topo_out_path (str): Folder to save processed landward elevation arrays
        dune_out_path (str): Folder to save dune crest height vectors
        mhw (float): Mean High Water in meters (relative to NAVD88)
        berm_elev (float): Elevation of the berm to offset dune heights
        seaward_edge (str): Array orientation (which edge faces the ocean)
    """
    for file in os.listdir(input_path):
        if not file.endswith('.npy'):
            continue

        base_name = os.path.splitext(file)[0]
        raw_array = np.load(os.path.join(input_path, file))

        # Normalize to MHW
        raw_array = raw_array - mhw
        raw_array[raw_array < -1] = -3  # Remove invalid low values

        # Re-orient array based on ocean-facing edge
        if seaward_edge == 'Right':
            raw_array = raw_array.T
        elif seaward_edge == 'Left':
            raw_array = np.flipud(raw_array.T)
        elif seaward_edge == 'Top':
            raw_array = np.flipud(raw_array)
        # 'Bottom' = no change

        # Output arrays
        processed_matrix = np.full((100, raw_array.shape[1]), fill_value=-3.0)
        dune_vector = np.full((raw_array.shape[1],), fill_value=-3.0)

        for i in range(raw_array.shape[1]):
            profile = np.flip(raw_array[:, i])  # Flip to go landward

            try:
                start_beach = next(idx for idx, val in enumerate(profile) if val > 0.5)
            except StopIteration:
                continue  # Skip column with no beach

            end_beach = start_beach + 8
            dune_elev = np.max(profile[start_beach:end_beach])
            dune_index = np.argmax(profile == dune_elev)
            start_island = dune_index + 1
            landward_profile = profile[start_island:]

            processed_matrix[0:len(landward_profile), i] = landward_profile
            dune_adj = dune_elev - berm_elev
            dune_vector[i] = dune_adj if dune_adj > 0 else 0.1

        # Convert meters to decameters for CASCADE
        processed_matrix *= 0.1
        dune_vector *= 0.1

        # Save arrays
        np.save(os.path.join(topo_out_path, f"{base_name}_topography.npy"), processed_matrix)
        np.save(os.path.join(dune_out_path, f"{base_name}_dune.npy"), dune_vector)

# === EXECUTE SCRIPT ========================================================

if __name__ == "__main__":
    process_topo_arrays(
        input_path=LOAD_PATH,
        topo_out_path=TOPO_SAVE_PATH,
        dune_out_path=DUNE_SAVE_PATH,
        mhw=MHW,
        berm_elev=BERM_ELEV,
        seaward_edge=SEAWARD_EDGE
    )
