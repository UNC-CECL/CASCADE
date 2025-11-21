# ==============================================================================
# Hatteras CASCADE Dune Offset Calculator
# ==============================================================================
# Author: Hannah Henry
# Last modified: 2025-08-1
# Description:
# This script calculates the average distance from the dune toe to a baseline
# transect for each barrier island domain and outputs both absolute and relative
# offsets for use in CASCADE model hindcasts.
# ==============================================================================

import pandas as pd
import numpy as np
import copy

# ------------------------------------------------------------------------------
# Load Data
# ------------------------------------------------------------------------------

# Path to input CSV exported from ArcGIS with dune-baseline intersections
Raw_Data_Path = r'/data/hatteras_init/island_offset/1978/1978_duneline_offset_raw.csv'

# Load CSV into DataFrame
Raw_Data = pd.read_csv(Raw_Data_Path)

# Extract only the relevant columns
Data_Dict = {
    'B3D_Grid': Raw_Data['domain_id'],     # Domain ID
    'Distance': Raw_Data['ORIG_LEN'],     # Distance from baseline to dune
    'Transect': Raw_Data['LineID']        # Transect ID within each domain
}
Data_DF = pd.DataFrame(Data_Dict)

# ------------------------------------------------------------------------------
# Organize Data by Domain
# ------------------------------------------------------------------------------

# Define domain IDs from 1 to 92 (inclusive)
B3D_Grids = range(30, 134)

# Split data into a list of domain-specific DataFrames
Grid_Values = []
for i in B3D_Grids:
    Subset = Data_DF[Data_DF['B3D_Grid'] == int(i)]
    Grid_Values.append(copy.deepcopy(Subset))

# ------------------------------------------------------------------------------
# Filter Duplicate Transects and Keep First Entry per Transect
# ------------------------------------------------------------------------------

Filtered_Data = []

for h in range(len(Grid_Values)):
    if Grid_Values[h].empty:
        print(f"Skipping B3D_Grid ID {B3D_Grids[h]} — no data")
        Filtered_Data.append([])  # maintain list alignment
        continue

    Min_Transect_Num = int(min(Grid_Values[h]['Transect']))
    Max_Transect_Num = int(max(Grid_Values[h]['Transect']))
    Transects = range(Min_Transect_Num, Max_Transect_Num + 1)

    Temp_Array = []
    for ll in Transects:
        transect_vals = Grid_Values[h][Grid_Values[h]['Transect'] == ll]
        if not transect_vals.empty:
            Temp_Array.append(copy.deepcopy(transect_vals.iloc[0]))  # Take first instance

    Filtered_Data.append(copy.deepcopy(Temp_Array))

# ------------------------------------------------------------------------------
# Calculate Average Dune Distance per Domain
# ------------------------------------------------------------------------------

Avg_Dune_Toe_Distance = []        # List of average distances
Avg_Dune_Toe_Distance_Dict = {}   # Dictionary for export with domain ID keys

for j in range(len(Grid_Values)):
    focus_data = Filtered_Data[j]
    if not focus_data:  # skip empty entries
        continue

    Temp_Distance_Array = []
    for k in range(len(focus_data)):
        Distance = focus_data[k][1]  # column index 1 = 'Distance'
        Temp_Distance_Array.append(copy.deepcopy(Distance))

    Mean_Distance = int(np.mean(Temp_Distance_Array))
    Avg_Dune_Toe_Distance.append(copy.deepcopy(Mean_Distance))
    Avg_Dune_Toe_Distance_Dict[str(B3D_Grids[j])] = Mean_Distance

# ------------------------------------------------------------------------------
# Calculate Relative Offset (baseline = 0)
# ------------------------------------------------------------------------------

baseline_distance = min(Avg_Dune_Toe_Distance)
Relative_dune_offset = np.subtract(Avg_Dune_Toe_Distance, baseline_distance)

# ------------------------------------------------------------------------------
# Export Results to CSV
# ------------------------------------------------------------------------------

# Absolute offset for each domain
Datum_dune_offset_df = pd.DataFrame(Avg_Dune_Toe_Distance_Dict, index=[0])
Datum_dune_offset_df.to_csv(
    r'C:\Users\hanna\PycharmProjects\CASCADE\data\hatteras_init\island_offset\1978\Dune_1978_Offset.csv',
    index=False
)

# Create a dataframe with domain IDs and relative dune offsets
relative_dune_df = pd.DataFrame({
    'ID': list(range(30, 134)),
    '1978': Relative_dune_offset  # Use the year this file represents
})

# Save with headers for easier loading in the model
relative_dune_df.to_csv(
    r'C:\Users\hanna\PycharmProjects\CASCADE\data\hatteras_init\island_offset\1978\Relative_Dune_1978_Offset.csv',
    index=False
)

# Relative offset in decimeters (for CASCADE input)
Relative_dune_offset_DM = np.round(Relative_dune_offset / 10).astype(int)
pd.DataFrame(Relative_dune_offset_DM).to_csv(
    r'C:\Users\hanna\PycharmProjects\CASCADE\data\hatteras_init\island_offset\1978\Relative_Dune_1978_Offset_DM.csv',
    index=False,
    header=False
)

# ------------------------------------------------------------------------------
# End of script
# ------------------------------------------------------------------------------