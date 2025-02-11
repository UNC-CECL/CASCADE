import copy

import numpy as np
import os
import copy

import pandas as pd

os.chdir("C:\\Users\\frank\\PycharmProjects\\CASCADE\\Run_output")
save_path = 'C:\\Users\\frank\\OneDrive - University of North Carolina at Chapel Hill\\Chapter 2\\Cascade_CSV_Outputs\\Hindcasts\\'


Base_Name_List = ['Geom_1',
                  'Geom_2',
                  'Geom_3',
                  'Geom_4',
                  'Geom_5']

run_name = ['Geom_1_Calibrated_Hindcast.npz',
            'Geom_2_Calibrated_Hindcast.npz',
            'Geom_3_Calibrated_Hindcast.npz',
            'Geom_4_Calibrated_Hindcast_2.npz',
            'Geom_5_Calibrated_Hindcast_2.npz']
#for geos in range(len(Base_Name_List)):
#    temp_run_name = copy.deepcopy(Base_Name_List[geos]+'_Calibrated_Hindcast_2.npz')
#    run_name.append(copy.deepcopy(temp_run_name))


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
    Combined_Dict[str(Base_Name_List[islands])] = copy.deepcopy(Change_Rate)
    b3d_upper_shoreface = b3d.x_s_TS
    b3d_lower_shoreface = b3d.x_t_TS
    Export_Dict_Upper_Shoreface = {'Upper_Shoreline':b3d_upper_shoreface}
    Export_Dict_Lower_Shoreface = {'Lower_Shoreline':b3d_lower_shoreface}

    Export_DF_Upper = pd.DataFrame(Export_Dict_Upper_Shoreface)
    Export_DF_Lower = pd.DataFrame(Export_Dict_Lower_Shoreface)
    Export_DF_Upper.to_csv(save_path+Base_Name_List[islands]+'_Hindcast_US.csv')
    Export_DF_Lower.to_csv(save_path+Base_Name_List[islands]+'_Hindcast_LS.csv')

z = 20