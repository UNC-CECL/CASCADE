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
ny = np.size(b3d)

Buffer_Domains = 15




Domains_of_Interest = range(Buffer_Domains,(len(b3d)-Buffer_Domains))


def Calculate_Island_Elevation_Metrics(cascade, buffer_length, number_barrier3d_models):
    Domains_of_Interest = range(buffer_length, (number_barrier3d_models - buffer_length))
    b3d = cascade.barrier3d

    Elevation_Change_TS = []
    Initial_Elevation_TS = []
    Final_Elevation_TS = []

    for k in Domains_of_Interest:
        Temp_B3D = b3d[k]
        Shoreline_Change_TS = Temp_B3D._ShorelineChangeTS

        # Find the initial and final island interior elevations
        initial_elev = Temp_B3D.DomainTS[0]
        final_elev = Temp_B3D.DomainTS[Temp_B3D.time_index-1]

        # Find the initial dune crest elevation
        initial_dune = Temp_B3D.DuneDomain[0]
        final_dune_elev = Temp_B3D.DuneDomain[Temp_B3D.time_index-1]

        # Rotate np arrays to concat with
        initial_dune_elev_r = np.rot90(initial_dune,k=3)
        final_dune_elev_r = np.rot90(final_dune_elev,k=3)

        # Add dunes to front of initial elevation
        combined_initial_elev = np.concatenate((initial_dune_elev_r,initial_elev),axis=0)

        distance_traveled = int(abs(np.sum(Shoreline_Change_TS[0:int(Temp_B3D.time_index)])))

        blank_cells = np.full((distance_traveled,int(Temp_B3D.BarrierLength)),-0.3)

        updated_final_array = np.concatenate((np.concatenate((blank_cells,final_dune_elev_r),axis=0),final_elev),axis=0)

        if len(combined_initial_elev) > len(updated_final_array):
            print('Add Extra Cells to final array')
            len_difference = len(combined_initial_elev) - len(updated_final_array)
            extra_cells = np.full((len_difference,int(Temp_B3D.BarrierLength)),-0.3)
            updated_final_array = np.concatenate((updated_final_array,extra_cells),axis=0)
        elif len(combined_initial_elev) < len(updated_final_array):
            print('Add extra cells to initial array')
            len_difference = len(updated_final_array) - len(combined_initial_elev)
            extra_cells = np.full((len_difference,int(Temp_B3D.BarrierLength)),-0.3)
            combined_initial_elev = np.concatenate((combined_initial_elev,extra_cells),axis=0)
        else:
            print('They have equal lengths')

        elev_difference = np.subtract(updated_final_array,combined_initial_elev)

        Elevation_Change_TS.append(copy.deepcopy(elev_difference))
        Initial_Elevation_TS.append(copy.deepcopy(combined_initial_elev))
        Final_Elevation_TS.append(copy.deepcopy(updated_final_array))

    return (Elevation_Change_TS,Initial_Elevation_TS,Final_Elevation_TS)


Elevation_Change_Output, Initial_Elevation_Output, Final_Elevation_Output = Calculate_Island_Elevation_Metrics(cascade=cascade,
                                                                                                        buffer_length=Buffer_Domains,
                                                                                                        number_barrier3d_models=ny)


c = 20
x = 10