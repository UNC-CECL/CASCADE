import copy

import numpy as np
import os
import copy

import pandas as pd

os.chdir("C:\\Users\\frank\\PycharmProjects\\CASCADE\\Run_output")
save_path = 'C:\\Users\\frank\\OneDrive - University of North Carolina at Chapel Hill\\Chapter 2\\Cascade_CSV_Outputs\\Hindcasts\\'


berm_elev_list = [1.27,
                  1.17,
                  1.14,
                  0.87,
                  0.87]

Base_Name_List = ['Geom_1',
                  'Geom_2',
                  'Geom_3',
                  'Geom_4',
                  'Geom_5']

run_name = []
show_name = []

background_erosion_value = [0]

#rmin = [0.1,0.55,0.9]
#rmax = [0.55,0.95,0.95]

rmin = [0.55]
rmax = [0.95]

for geos in range(len(Base_Name_List)):
    for i in range(0,len(rmin)):
        run_name.append(copy.deepcopy(Base_Name_List[geos]+'_'+str(background_erosion_value[0])+'_'+str(rmin[i])+'_OB.npz'))
        run_name.append(copy.deepcopy(Base_Name_List[geos]+'_'+str(background_erosion_value[0])+'_'+str(rmin[i])+'_NB.npz'))


Combined_Dict = {}
for islands in range(len(run_name)):
    Geom_Dict = {}
    output = np.load(run_name[islands], allow_pickle=True)
    cascade = output["cascade"]
    cascade = cascade[0]
    b3d = cascade._barrier3d[0]

    Shoreline_Change_Rate = cascade._brie_coupler.brie.x_s_dt
    all_shoreline_change = cascade._brie_coupler.brie.x_s_save
    all_shoreline_change_0 = cascade._brie_coupler.brie.x_s_save[0][0]
    all_shoreline_change_end = cascade._brie_coupler.brie.x_s_save[0][-1]
    Change_Rate = -(all_shoreline_change_end - all_shoreline_change_0)/len(cascade._brie_coupler.brie.x_s_save[0])
    Combined_Dict[str(run_name[islands])] = copy.deepcopy(Change_Rate)
    b3d_upper_shoreface = b3d.x_s_TS
    b3d_lower_shoreface = b3d.x_t_TS
    Export_Dict_Upper_Shoreface = {'Upper_Shoreline':b3d_upper_shoreface}
    Export_Dict_Lower_Shoreface = {'Lower_Shoreline':b3d_lower_shoreface}

    Export_DF_Upper = pd.DataFrame(Export_Dict_Upper_Shoreface)
    Export_DF_Lower = pd.DataFrame(Export_Dict_Lower_Shoreface)
    #Export_DF_Upper.to_csv(save_path+Base_Name_List[islands]+'_Hindcast_US.csv')
    #Export_DF_Lower.to_csv(save_path+Base_Name_List[islands]+'_Hindcast_LS.csv')

z = 20