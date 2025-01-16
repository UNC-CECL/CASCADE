#Test changes and graphs for 1 Cascade instance

import copy

import numpy as np
import pandas as pd
from scipy import stats as st
import os

os.chdir('E:\\Model_Runs')

run_name ='OCR_IH_Status_Quo_S84_Accretional_Sink'

output = np.load(run_name + ".npz", allow_pickle=True)
cascade = output["cascade"]
cascade = cascade[0]
b3d = cascade.barrier3d
ny = np.size(b3d)

Buffer_Domains = 15

Domains_of_Interest = range(Buffer_Domains,(len(b3d)-Buffer_Domains))


for k in Domains_of_Interest:
    Temp_B3D = b3d[k]
    Shoreline_Change_TS = Temp_B3D._ShorelineChangeTS
    #for k in range(0,Temp_B3D.time_index):
    initial_elev = Temp_B3D.DomainTS[0]
    final_elev = Temp_B3D.DomainTS[Temp_B3D.time_index-1]
    distance_traveled = int(abs(np.sum(Shoreline_Change_TS[0:int(Temp_B3D.time_index)])))

    blank_cells = np.full((distance_traveled,int(Temp_B3D.BarrierLength)),-0.3)

    updated_final_array = np.concatenate((blank_cells,final_elev),axis=0)

    if len(initial_elev) > len(updated_final_array):
        print('Add Extra Cells to final array')
        len_difference = len(initial_elev) - len(updated_final_array)
        extra_cells = np.full((len_difference,int(Temp_B3D.BarrierLength)),-0.3)
        updated_final_array = np.concatenate((updated_final_array,extra_cells),axis=0)
    elif len(initial_elev) < len(updated_final_array):
        print('Add extra cells to initial array')
        len_difference = len(updated_final_array) - len(initial_elev)
        extra_cells = np.full((len_difference,int(Temp_B3D.BarrierLength)),-0.3)
        initial_elev = np.concatenate((initial_elev,extra_cells),axis=0)
    else:
        print('They have equal lengths')

    elev_difference = np.subtract(updated_final_array,initial_elev)

    c = 20
x = 10