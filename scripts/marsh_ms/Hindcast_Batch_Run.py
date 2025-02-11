# Run a batch of different Cascade Runs
import copy

import numpy as np
import os
from cascade.cascade import Cascade

# ###############################################################################
# Input Variables
# ###############################################################################
os.chdir('C:\\Users\\frank\\PycharmProjects\\CASCADE')
data_base_path = 'C:\\Users\\frank\\PycharmProjects\\CASCADE\\data\\marsh_init_data\\'
# Specify variables to use in calling function
# Dune height path name
d_file = []
e_file = []

e_names = ['Metompkin_Marsh_Topo.npy',
           'Metompkin_Bay_Topo.npy',
           'Hog_Topo.npy',
           'Wreck_Topo.npy',
           'Smith_Topo.npy']

d_names = ['Metompkin_Marsh_Dune.npy',
           'Metompkin_Bay_Dune.npy',
           'Hog_Dune.npy',
           'Wreck_Dune.npy',
           'Smith_Dune.npy']

d_names = [
    'GEOM_1_Dunes.npy',
    'GEOM_2_Dunes.npy',
    'GEOM_3_Dunes.npy',
    'GEOM_4_Dunes.npy',
    'GEOM_5_Dunes.npy'
           ]

Base_Name_List = ['Geom_1',
                  'Geom_2',
                  'Geom_3',
                  'Geom_4',
                  'Geom_5']

for runs in range(len(e_names)):
    d_file.append(copy.deepcopy(data_base_path+d_names[runs]))
    e_file.append(copy.deepcopy(data_base_path+e_names[runs]))


marsh_path_list = ["C:\\Users\\frank\\PycharmProjects\\CASCADE\\data\\marsh_init_data\\BMFT_Marsh_Width_500.mat",
                   "C:\\Users\\frank\\PycharmProjects\\CASCADE\\data\\marsh_init_data\\BMFT_Marsh_Width_500.mat",
                   "C:\\Users\\frank\\PycharmProjects\\CASCADE\\data\\marsh_init_data\\BMFT_Marsh_Width_500.mat",
                   "C:\\Users\\frank\\PycharmProjects\\CASCADE\\data\\marsh_init_data\\Marsh_250.mat",
                   "C:\\Users\\frank\\PycharmProjects\\CASCADE\\data\\marsh_init_data\\BMFT_Marsh_Width_500.mat"]
marsh_width_list = [500,500,500,250,500]

berm_elev_list = [1.27,
                  1.17,
                  1.14,
                  0.87,
                  0.87]

'''berm_elev_list = [1.9,
                  1.9,
                  1.9,
                  1.9,
                  1.9]'''


# Storm file path name
c_wd = os.getcwd()
nt_run = 23  # Number of years model will run
run_name = []
#storm_name = 'C:\\Users\\frank\\PycharmProjects\\CASCADE\\data\\marsh_init_data\\VCR_Storms_1994_2020.npy'

storm_name = ['C:\\Users\\frank\\PycharmProjects\\CASCADE\\data\\marsh_init_data\\VCR_Storms_Geom_1_1994_2020.npy',
    'C:\\Users\\frank\\PycharmProjects\\CASCADE\\data\\marsh_init_data\\VCR_Storms_Geom_2_1994_2020.npy',
    'C:\\Users\\frank\\PycharmProjects\\CASCADE\\data\\marsh_init_data\\VCR_Storms_Geom_3_1994_2020.npy',
    'C:\\Users\\frank\\PycharmProjects\\CASCADE\\data\\marsh_init_data\\VCR_Storms_Geom_4_1994_2020.npy',
    'C:\\Users\\frank\\PycharmProjects\\CASCADE\\data\\marsh_init_data\\VCR_Storms_Geom_5_1994_2020.npy'
              ]

background_erosion_value = [x / 10.0 for x in range(0,-100,-5)]


rmin = [0.2]
rmax = [0.75]

for geos in range(len(Base_Name_List)):
    temp_run_name = []
    for i in range(0,len(background_erosion_value)):
        temp_run_name.append(copy.deepcopy(Base_Name_List[geos]+'_'+str(background_erosion_value[i])+'_DGR.45'))
    run_name.append(copy.deepcopy(temp_run_name))

num_of_batches = 1
marsh_dynamics_on = True
run_set_RSLR = False

# Name storms
set_RSLR = 0.00563


number_barrier3d_models = 1
#rmin = 0.55
#rmax = 0.95

# Define function
def Batch_Runs(
    nt,
    name,
    storm_file,
    alongshore_section_count,
    num_cores,
    rmin,  # the remaining variables are arrays
    rmax,
    elevation_file,
    dune_file,
    background_erosion,
    marsh_width,
    marsh_path,
    berm_elevation,
    sea_level_rise_rate=0.008,  # not an array
    sea_level_constant=True,  # not an array
    enable_shoreline_offset=False,
    shoreline_offset=[0],
    marsh_dynamics=False,
    user_inputed_RSLR=False,
    user_inputed_RSLR_rate=[],
):

    # ###############################################################################
    # 9 - connect cascade domains (human management) with AST
    # ###############################################################################

    # --------- INITIALIZE ---------
    datadir = "data/"
    cascade = Cascade(
        datadir,
        name,
        storm_file=storm_file,
        elevation_file=elevation_file,
        dune_file=dune_file,
        parameter_file="Alongshore_Test-parameters.yaml",
        wave_height=1,
        wave_period=7,
        wave_asymmetry=0.8,
        wave_angle_high_fraction=0.2,
        sea_level_rise_rate=sea_level_rise_rate,
        sea_level_rise_constant=sea_level_constant,
        background_erosion=background_erosion,
        alongshore_section_count=alongshore_section_count,
        time_step_count=nt,
        min_dune_growth_rate=rmin,
        max_dune_growth_rate=rmax,
        num_cores=num_cores,
        roadway_management_module=False,  # no roadway management
        alongshore_transport_module=True,  # Is there brie coupling?
        beach_nourishment_module=False,  # no beach nourishment
        community_economics_module=False,  # no community dynamics
        enable_shoreline_offset=enable_shoreline_offset,  # Bool
        shoreline_offset=shoreline_offset,
        marsh_dynamics=marsh_dynamics,
        user_inputed_RSLR = user_inputed_RSLR,
        user_inputed_RSLR_rate = user_inputed_RSLR_rate,
        marsh_path=marsh_path,
        marsh_width=marsh_width,
        berm_elevation=berm_elevation,

    )
    # --------- LOOP ---------
    for time_step in range(nt - 1):

        # Print time step to screen (NOTE: time_index in each model is time_step+1)
        print("\r", "Time Step: ", time_step, end="")
        # print('First Loop')
        cascade.update()
        if cascade.b3d_break:
            break

    # --------- SAVE ---------
    save_directory = "Run_Output/"
    cascade.save(save_directory)

    return cascade

# Call a batch of functions
for j in range(len(run_name)):
    for k in range(len(run_name[0])):
        Batch_Runs(
            nt=nt_run,
            name=run_name[j][k],
            storm_file=storm_name[j],
            alongshore_section_count=number_barrier3d_models,
            num_cores=3,
            rmin=rmin[k],
            rmax=rmax[k],
            elevation_file=e_file[j],
            dune_file=d_file[j],
            background_erosion=background_erosion_value[0],
            sea_level_constant=True,  # not an array
            enable_shoreline_offset=False,
            marsh_dynamics=marsh_dynamics_on,
            sea_level_rise_rate=set_RSLR,
            user_inputed_RSLR=run_set_RSLR,
            user_inputed_RSLR_rate=set_RSLR,
            marsh_path=marsh_path_list[j],
            marsh_width=marsh_width_list[j],
            berm_elevation=berm_elev_list[j]
        )
        os.chdir('C:\\Users\\frank\\PycharmProjects\\CASCADE')
        print('Finished '+str(run_name[j][k]))
