# Benton Franklin
# Last modified 1/22/2025
# Used to generate shoreline locatiom for Ocracoke Cascade model domains
# Script designed to import .csv file generated from ArcGis that shows the intersection
# between the shoreline and some baseline offshore transect

import pandas as pd
import numpy as np
import copy

# Specify location of the wanted data

#Raw_Data_Path = 'C:\\Users\\frank\\OneDrive - University of North Carolina at Chapel Hill\\Chapter 3\\Dune Offsets\\Raw_1976_Offsets.csv'
Raw_Data_Path = 'C:\\Users\\frank\\OneDrive - University of North Carolina at Chapel Hill\\Chapter 3\\Revised_Offshore_Datum\\Shoreline_2020_Final.csv'


Raw_Data = pd.read_csv(Raw_Data_Path)

Data_Dict = {'B3D_Grid':Raw_Data['B3D_Grid'],
             'Distance':Raw_Data['ORIG_LEN'],
             'Transect':Raw_Data['LineID']}
Data_DF = pd.DataFrame(Data_Dict)

Mean_Distance = []

# Seperate out by grid
B3D_Grids = range(11,50)
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

# Calculate distance between the dune and the datum
Avg_Dune_Toe_Distance = [] # Save array
Avg_Dune_Toe_Distance_Dict = {} # Save in a dict

for j in range(len(Grid_Values)):
    focus_data = Filtered_Data[j]
    Temp_Distance_Array = []
    for k in range(len(focus_data)):
        Distance = focus_data[k][1]
        Temp_Distance_Array.append(copy.deepcopy(Distance))
    Mean_Distance = int(np.mean(Temp_Distance_Array))
    Avg_Dune_Toe_Distance.append(copy.deepcopy(Mean_Distance))
    Avg_Dune_Toe_Distance_Dict[str(B3D_Grids[j])] = copy.deepcopy(Mean_Distance)

# Calculate the relative distance

baseline_distance = min(Avg_Dune_Toe_Distance)
Relative_dune_offset = np.subtract(Avg_Dune_Toe_Distance,baseline_distance)

Datum_dune_offset_df = pd.DataFrame(Avg_Dune_Toe_Distance_Dict, index=[0])
x = 10
Datum_dune_offset_df.to_csv(path_or_buf='C:\\Users\\frank\\OneDrive - University of North Carolina at Chapel Hill\\Chapter 3\\Revised_Offshore_Datum\\Shoreline_2020_Datum_Offset.csv')
#Relative_dune_offset_df.to_csv(path_or_buf='C:\\Users\\frank\\OneDrive - University of North Carolina at Chapel Hill\\Chapter 3\\Dune Offsets\\Relative_Dune_2019_Offset.csv')

# Convert and save relative distance in DM for use in CASCADE model runs
#Relative_dune_offset_DM = np.round(np.divide(Relative_dune_offset_df,10))
#Relative_dune_offset_DM.to_csv(path_or_buf='C:\\Users\\frank\\OneDrive - University of North Carolina at Chapel Hill\\Chapter 3\\Dune Offsets\\Relative_Dune_2019_Offset_DM.csv')

'''
# Calculate the mean distance between datum and dune line for each B3D grid
for k in range(11,50):
    B3D_Grid_Set = Raw_Data[Raw_Data['B3D_Grid'] == k]
    Data_Dict[str(k)] = copy.deepcopy(B3D_Grid_Set)
    h = B3D_Grid_Set['ORIG_LEN']
    hm = np.mean(B3D_Grid_Set['ORIG_LEN'])
    Mean_Distance.append(copy.deepcopy(hm))

# Subtract the smallest distance to have distance 0 represent the baseline value
Relative_Distance = np.subtract(Mean_Distance,min(Mean_Distance))

Relative_Distance_dm = np.round((Relative_Distance/10))

Final_Distance = Relative_Distance_dm*10
'''
