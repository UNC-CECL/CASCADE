# Run a batch of different Cascade Runs that have been calibrated from Hindcast batch run
import copy

import numpy as np
import os
from cascade.cascade import Cascade

# ###############################################################################
# Input Variables
# ###############################################################################
#os.chdir('C:\\Users\\frank\\PycharmProjects\\CASCADE')
os.chdir('/Users/ceclmac/PycharmProjects/CASCADE')

data_base_path = '/Users/ceclmac/PycharmProjects/CASCADE/data/marsh_init_data/'
# Specify variables to use in calling function
# Dune height path name
d_file = []
e_file = []

Start_Year = 1994

e_names = ['Metompkin_Marsh_Topo.npy',
           'Metompkin_Bay_Topo.npy',
           'Hog_Topo.npy',
           'Wreck_Topo.npy',
           'Smith_Topo.npy']

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


marsh_path_list = [data_base_path+"BMFT_Marsh_Width_500.mat",
                   data_base_path+"BMFT_Marsh_Width_500.mat",
                   data_base_path+"BMFT_Marsh_Width_500.mat",
                   data_base_path+"Marsh_250.mat",
                   data_base_path+"BMFT_Marsh_Width_500.mat"]
marsh_width_list = [500,500,500,250,500]

berm_elev_list = 1.27

# Storm file path name
c_wd = os.getcwd()

if Start_Year == 1980:
    nt_run = 26  # Number of years model will run
    storm_name = data_base_path + 'VCR_Storms_Combined_Time_1980_2006.npy'
elif Start_Year == 1994:
    nt_run = 23
    storm_name = data_base_path + 'VCR_Storms_Combined_Time_1994_2020.npy'
else:
    nt_run = 11  # Number of years model will run
    storm_name = data_base_path+'VCR_Storms_Combined_Time_2006_2020.npy'

run_name = []

for geos in range(len(Base_Name_List)):
    if Start_Year == 1980:
        temp_run_name =  copy.deepcopy(Base_Name_List[geos]+'_Calibrated_Hindcast_1980_2006_T2')
        background_erosion_value = [-9.1, -0.75, -5.2, -1.3, -3.25]
    elif Start_Year == 1994:
        temp_run_name =  copy.deepcopy(Base_Name_List[geos]+'_Calibrated_Hindcast_1994_2017_T6')
        background_erosion_value = [-11, -3.75, -5.5, -6.75, -6]
    else:
        temp_run_name =  copy.deepcopy(Base_Name_List[geos]+'_Calibrated_Hindcast_2006_2020_T2')
        background_erosion_value = [-9.1, -0.75, -5.2, -1.3, -3.25]
    run_name.append(copy.deepcopy(temp_run_name))

num_of_batches = 1
marsh_dynamics_on = True
run_set_RSLR = False

# Name storms
set_RSLR = 0.00563


number_barrier3d_models = 1
rmin = 0.2
rmax = 0.75

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
    Batch_Runs(
        nt=nt_run,
        name=run_name[j],
        storm_file=storm_name,
        alongshore_section_count=number_barrier3d_models,
        num_cores=3,
        rmin=rmin,
        rmax=rmax,
        elevation_file=e_file[j],
        dune_file=d_file[j],
        background_erosion=background_erosion_value[j],
        sea_level_constant=True,  # not an array
        enable_shoreline_offset=False,
        marsh_dynamics=marsh_dynamics_on,
        sea_level_rise_rate=set_RSLR,
        user_inputed_RSLR=run_set_RSLR,
        user_inputed_RSLR_rate=set_RSLR,
        marsh_path=marsh_path_list[j],
        marsh_width=marsh_width_list[j],
        berm_elevation = berm_elev_list

    )
    #os.chdir('C:\\Users\\frank\\PycharmProjects\\CASCADE')
    os.chdir('/Users/ceclmac/PycharmProjects/CASCADE')
    print(' Finished '+str(run_name[j]))
