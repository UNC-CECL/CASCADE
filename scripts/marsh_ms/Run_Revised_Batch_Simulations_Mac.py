# Run a batch of different Cascade Runs
import copy

import numpy as np
import os
from cascade.cascade import Cascade

# ###############################################################################
# Input Variables
# ###############################################################################
os.chdir('/Users/ceclmac/PycharmProjects/CASCADE')

storm_intensity = '10%'
RSLR_Rate = 'IL'

data_base_path = '/Users/ceclmac/PycharmProjects/CASCADE/data/marsh_init_data/'
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


marsh_path_list = [data_base_path+"BMFT_Marsh_Width_1500.mat",
                   data_base_path+"BMFT_Marsh_Width_500.mat",
                   data_base_path+"BMFT_Marsh_Width_500.mat",
                   data_base_path+"Marsh_250.mat",
                   data_base_path+"BMFT_Marsh_Width_500.mat"]
marsh_width_list = [1500,500,500,250,500]


berm_elev_list = 1.03

# RSLR Data
if RSLR_Rate == 'IL':
    RSLR_Data = np.load(data_base_path+'Updated_2025_IL_SLR.npy')
elif RSLR_Rate == 'I':
    RSLR_Data = np.load(data_base_path+'Updated_2025_Int_SLR.npy')
elif RSLR_Rate == 'IH':
    RSLR_Data = np.load(data_base_path+'Updated_2025_High_SLR.npy')


# Storm file path name
c_wd = os.getcwd()
nt_run = 125  # Number of years model will run
storm_name = []

if storm_intensity == 'Baseline':
    for k in range(50,100):
        storm_name.append(copy.deepcopy(data_base_path+'Baseline/StormList_'+str(k)+'_baseline_RSD.npy'))
    #storm_name[32] = copy.deepcopy(data_base_path + 'Baseline/StormList_32_baseline_N.npy')
    storm_add = '_Baseline'
elif storm_intensity == '5%':
    for k in range(0,50):
        storm_name.append(copy.deepcopy(data_base_path+'Five_Percent_Increase/StormList_'+str(k)+'_5_percent_increase_RSD.npy'))
    storm_add = '_5'

elif storm_intensity == '10%':
    for k in range(0,50):
        storm_name.append(copy.deepcopy(data_base_path+'Ten_Percent_Increase/StormList_'+str(k)+'_10_percent_increase_RSD.npy'))
    storm_add = '_10'

background_erosion_value = [-9.1, -0.75, -5.2, -2.6, -3.25]

rmin = [0.05]
rmax = [0.55]

run_name = []

for geos in range(len(Base_Name_List)):
    temp_run_name = []
    for i in range(0,len(storm_name)):
        temp_run_name.append(copy.deepcopy(Base_Name_List[geos]+'_'+str(RSLR_Rate)+str(storm_add)+'_S'+str(i)+'_New_Sink'))
    run_name.append(copy.deepcopy(temp_run_name))

num_of_batches = 1
marsh_dynamics_on = [True,False,True,True,True]
run_set_RSLR = True

# Name storms
set_RSLR = RSLR_Data
number_barrier3d_models = 1


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
    save_directory = '/Volumes/BF_Backup/Chapter_2_Runs'
    cascade.save(save_directory)

    return cascade

# Call a batch of functions
for j in range(3,4):#len(run_name)):
    for k in range(0,50):#len(storm_name)):
        Batch_Runs(
            nt=nt_run,
            name=run_name[j][k],
            storm_file=storm_name[k],
            alongshore_section_count=number_barrier3d_models,
            num_cores=3,
            rmin=rmin[0],
            rmax=rmax[0],
            elevation_file=e_file[j],
            dune_file=d_file[j],
            background_erosion=background_erosion_value[j],
            sea_level_constant=False,  # not an array
            enable_shoreline_offset=False,
            marsh_dynamics=marsh_dynamics_on[j],
            sea_level_rise_rate=0.01,
            user_inputed_RSLR=run_set_RSLR,
            user_inputed_RSLR_rate=set_RSLR,
            marsh_path=marsh_path_list[j],
            marsh_width=marsh_width_list[j],
            berm_elevation=berm_elev_list
        )
        #os.chdir('C:\\Users\\frank\\PycharmProjects\\CASCADE')
        os.chdir('/Users/ceclmac/PycharmProjects/CASCADE')

        print(' Finished '+str(run_name[j][k]))
