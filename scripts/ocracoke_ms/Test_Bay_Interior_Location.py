#Test changes and graphs for 1 Cascade instance

import copy

import numpy as np
import pandas as pd
from scipy import stats as st
import os

os.chdir('E:\\Model_Runs')

run_name ='OCR_IL_Status_Quo_S2_Erosional_Sink'

output = np.load(run_name + ".npz", allow_pickle=True)
cascade = output["cascade"]
cascade = cascade[0]
b3d = cascade.barrier3d
shoreface_ts = cascade._brie_coupler.brie.x_s_save
ny = np.size(b3d)

Buffer_Domains = 15

run_name_m ='OCR_IL_Nourishment_S2_Erosional_Sink'

output_n = np.load(run_name_m + ".npz", allow_pickle=True)
cascade_n = output_n["cascade"]
cascade_n = cascade_n[0]
b3d_n = cascade_n.barrier3d
ny = np.size(b3d)
shoreface_ts_n = cascade_n._brie_coupler.brie.x_s_save

focus_xs = shoreface_ts[50]
focus_xn = shoreface_ts_n[50]

Domains_of_Interest = range(Buffer_Domains,(len(b3d)-Buffer_Domains))
x =20

def Calculate_Bay_Distance(cascade, buffer_length, number_barrier3d_models):
    Domains_of_Interest = range(buffer_length, (number_barrier3d_models - buffer_length))
    b3d = cascade.barrier3d

    Domain_TS = []

    for k in Domains_of_Interest:
        Temp_B3D = b3d[k]
        Shoreline_Change_TS = Temp_B3D._ShorelineChangeTS
        Mean_Distance = []
        for years in range((Temp_B3D.time_index-1)):
            Yearly_Vals = []
            raw_elev = Temp_B3D.DomainTS[years]
            flipped_array = copy.deepcopy(np.flip(raw_elev))
            for cols in range((raw_elev.shape[1])):
                Zero_Row_Data = (np.where(raw_elev[:,cols]>0)[0])
                if np.any(Zero_Row_Data):
                    Zero_Data_Index = Zero_Row_Data.max()
                else:
                    Zero_Data_Index = 0
                if years == 0:
                    Shoreline_Change_Dif = 0
                else:
                    Shoreline_Change_Dif = np.sum(Shoreline_Change_TS[:years])

                Zero_Data_Index_Final = Zero_Data_Index - Shoreline_Change_Dif
                Yearly_Vals.append(copy.deepcopy(Zero_Data_Index_Final))
            Mean_Distance.append(copy.deepcopy(np.mean(Yearly_Vals)))
        print(k)
        Domain_TS.append(copy.deepcopy(Mean_Distance))

    return (Domain_TS)


Domain_Bay_Shoreline = Calculate_Bay_Distance(cascade=cascade,
                                              buffer_length=Buffer_Domains,
                                              number_barrier3d_models=ny)


c = 20
x = 10