# Run a batch of different Cascade Runs

import numpy as np
import time

import matplotlib.pyplot as plt

import os
import imageio


from cascade.cascade import Cascade


# ###############################################################################
# Example of Running Function 1
# ###############################################################################
os.chdir("/Users/ceclmac/PycharmProjects/CASCADE")
# Specify variables to use in calling function
# Dune height path name
d_file = "/Users/ceclmac/PycharmProjects/CASCADE/data/marsh_init_data/Smith_Dune.npy"
e_file = "/Users/ceclmac/PycharmProjects/CASCADE/data/marsh_init_data/Smith_Topo.npy"
# Storm file path name
#RSLR_File = '/data/marsh_init_data/Low_SLR.npy'
Num_Storms = 50
c_wd = os.getcwd()
nt_run = 150  # Number of years model will run
run_name = []
storm_name = []
for i in range(0,Num_Storms):
    run_name.append('Hog_10_Percent_High_RSLR_'+str(i))
    storm_name.append('/Users/ceclmac/PycharmProjects/CASCADE/data/marsh_init_data/Ten_Percent_Increase/StormList_'+
                       str(i)+'_10_percent_increase.npy')
num_of_batches = 1
marsh_dynamics_on = True
run_set_RSLR = True

# Name storms


if run_set_RSLR == True:
    set_RSLR = np.load('/Users/ceclmac/PycharmProjects/CASCADE/data/marsh_init_data/Low_SLR.npy')
else:
    set_RSLR = []


number_barrier3d_models = 1
rmin = 0.55
rmax = 0.95

# Call function


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

for j in range(0,num_storms):

    Batch_Runs(
        nt=nt_run,
        name=run_name[j],
        storm_file=storm_name[j],
        alongshore_section_count=number_barrier3d_models,
        num_cores=3,
        rmin=rmin,
        rmax=rmax,
        elevation_file=e_file,
        dune_file=d_file,
        background_erosion=-1.00,
        sea_level_constant=False,  # not an array
        enable_shoreline_offset=False,
        marsh_dynamics=True,
        sea_level_rise_rate=0.004,
        user_inputed_RSLR=run_set_RSLR,
        user_inputed_RSLR_rate=set_RSLR,
    )
    os.chdir("/Users/ceclmac/PycharmProjects/CASCADE")
    print('Finished '+str(run_name[j]))
