# Benton Franklin
# 9/11/2024

import copy
import pandas as pd
import numpy as np

# Script to load in points from ARCPro that are measured reletive to some fixed datum
# and calculate the average road location for each cell reletive to the datum

# Set path to data you want to load
Road_Data_Base_Path = '/Users/rsahrae/PycharmProjects/FAll2024/PeaIsland/Hindcast/'
Road_Data_Name = '1978_duneline_offset_raw.csv'

# Load data
Road_Data = pd.read_csv(str(Road_Data_Base_Path+Road_Data_Name))

# Seperate out data to only include relevant attributes
Data_Dict = {'B3D_Grid':Road_Data['ID'],
             'Distance':Road_Data['ORIG_LEN'],
              'Transect':Road_Data['LineID']}

# Create DF from dictionary
Data_DF = pd.DataFrame(Data_Dict)

B3D_Grids = range(17,58)

Grid_Values = []

for i in (B3D_Grids):
    Subset = Data_DF[Data_DF['B3D_Grid']==int(i)]
    Grid_Values.append(copy.deepcopy(Subset))

# Filter data to remove repeated transects
Filtered_Data = []

for h in range(len(Grid_Values)):

    Min_Transect_Num =int(min(Grid_Values[h]['Transect']))
    Max_Transect_Num = int(max(Grid_Values[h]['Transect']))
    Transects = range(Min_Transect_Num,(Max_Transect_Num+1))
    Temp_Array = []
    for ll in Transects:
        transect_vals = Grid_Values[h][Grid_Values[h]['Transect'] == ll]
        Temp_Array.append(copy.deepcopy(transect_vals.iloc[0]))
    Filtered_Data.append(copy.deepcopy(Temp_Array))

# Calculate the average distance for each B3D grid cell
Avg_Road_Distance = [] # Save array
Avg_Road_Distance_Dict = {} # Save in a dict

for j in range(len(Grid_Values)):
    focus_data = Filtered_Data[j]
    Temp_Distance_Array = []
    for k in range(len(focus_data)):
        Distance = focus_data[k][1]
        Temp_Distance_Array.append(copy.deepcopy(Distance))
    Mean_Distance = int(np.mean(Temp_Distance_Array))
    Avg_Road_Distance.append(copy.deepcopy(Mean_Distance))
    Avg_Road_Distance_Dict[str(B3D_Grids[j])] = copy.deepcopy(Mean_Distance)

Export_DF = pd.DataFrame(Avg_Road_Distance_Dict, index=[0])
Export_DF.to_csv(path_or_buf=(Road_Data_Base_Path+'Avg_NC_12_B3D_Distances_2011.csv'))

# Calculate the relative distance from dune line to roadway

# Load the dune offset location data
Dune_Offset = '/Users/rsahrae/PycharmProjects/FAll2024/PeaIsland/Hindcast/PeaIslandHindcastAttributes/Relative_Dune_1978_Offset_PeaIsland.csv'
Dune_Location_Data = pd.read_csv(Dune_Offset)
Dune_Location_Data = Dune_Location_Data.drop(Dune_Location_Data.columns[0], axis=1)

# Convert data to np arrays for subtraction
Dune_Values = np.array(Dune_Location_Data)
Road_Values = np.array(Export_DF)

Dune_Road_Dif = np.subtract(Road_Values,Dune_Values)
Dune_Road_Dif_DM = np.divide(Dune_Road_Dif,10)
Corrected_Dune_Dif_DM = np.subtract(Dune_Road_Dif_DM,2)
Export_Road_Distance = abs(np.round(Corrected_Dune_Dif_DM))
Export_Road_Distance_DF = pd.DataFrame(Export_Road_Distance)

# Save DF as .csv
Export_Road_Distance_DF.to_csv(path_or_buf=(Road_Data_Base_Path+'Road_Dune_Distance_2011.csv'))