import copy

import numpy as np
import os
import copy

os.chdir("C:\\Users\\frank\\PycharmProjects\\CASCADE\\Run_output")

#run_name = 'Metompkin_Marsh_Five_Percent_Storms_High_RSLR_0_Hindcast_Test.npz'
run_name = 'Metompkin_Marsh_Five_Percent_Storms_High_RSLR_1_Hindcast_Test_0.npz'
run_name = 'Hog_10_Percent_High_RSLR_0.npz'

marsh_edge_TS = []
forest_edge_TS = []
bmft_elevation_TS = []
marsh_offset_TS = []
shoreline_TS = []
shoreline_toe_TS = []

output = np.load(run_name, allow_pickle=True)
cascade = output["cascade"]
cascade = cascade[0]
b3d = cascade._barrier3d[0]

Shoreline_Change_Rate = cascade._brie_coupler.brie.x_s_dt
all_shoreline_change = cascade._brie_coupler.brie.x_s_save
all_shoreline_change_0 = cascade._brie_coupler.brie.x_s_save[0][0]
all_shoreline_change_end = cascade._brie_coupler.brie.x_s_save[0][-1]
Change_Rate = (all_shoreline_change_end - all_shoreline_change_0)/len(cascade._brie_coupler.brie.x_s_save[0])

shoreline_TS.append(b3d.x_s_TS)
shoreline_toe_TS.append(b3d.x_t_TS)

Final_Shoreline_Distance = (shoreline_TS[0][-1]-shoreline_TS[0][0])

Shoreline_change_Rate = Final_Shoreline_Distance/len(shoreline_TS[0])

z = 20