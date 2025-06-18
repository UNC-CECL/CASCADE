import copy

import numpy as np
import os
import copy

import pandas as pd

os.chdir('/Users/ceclmac/PycharmProjects/CASCADE/Run_output')



Base_Name_List = ['Geom_1',
                  'Geom_2',
                  'Geom_3',
                  'Geom_4',
                  'Geom_5']


background_erosion_value = [x / 100.0 for x in range(0,-1500,-25)]


rmin = [0.2]
rmax = [0.75]
run_name = []
for geos in range(len(Base_Name_List)):
    temp_run_name = []
    for i in range(0,len(background_erosion_value)):
        temp_run_name.append(copy.deepcopy(Base_Name_List[geos]+'_'+str(background_erosion_value[i])+'_DGR.30_1994_2017_NS.npz'))
        #temp_run_name.append(copy.deepcopy(Base_Name_List[geos]+'_'+str(background_erosion_value[i])+'_DGR.45_1980_2006_NS.npz'))

    run_name.append(copy.deepcopy(temp_run_name))

Combined_Dict = {}
for islands in range(len(Base_Name_List)):
    Geom_Dict = {}
    for sinks in range(len(background_erosion_value)):
        output = np.load(run_name[islands][sinks], allow_pickle=True)
        cascade = output["cascade"]
        cascade = cascade[0]
        b3d = cascade._barrier3d[0]

        Shoreline_Change_Rate = cascade._brie_coupler.brie.x_s_dt
        all_shoreline_change = cascade._brie_coupler.brie.x_s_save
        all_shoreline_change_0 = cascade._brie_coupler.brie.x_s_save[0][0]
        all_shoreline_change_end = cascade._brie_coupler.brie.x_s_save[0][-1]
        Change_Rate = -(all_shoreline_change_end - all_shoreline_change_0)/len(cascade._brie_coupler.brie.x_s_save[0])
        Geom_Dict[str(background_erosion_value[sinks])] = copy.deepcopy(Change_Rate)
    Combined_Dict[str(Base_Name_List[islands])] = copy.deepcopy(Geom_Dict)

save_path = '/Users/ceclmac/OneDrive - University of North Carolina at Chapel Hill/Chapter 2/Hindcasts/'
Export_DF = pd.DataFrame(Combined_Dict)
#Export_DF.to_csv(save_path+'DGR_45_1980_2006.csv')

z = 20