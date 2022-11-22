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
os.chdir('/Users/ceclmac/PycharmProjects/CASCADE')
# Specify variables to use in calling function
# Elevation file path name
e_file = "/B3D_Inputs/Marsh_Test_Inputs/InitElevHog.npy"
# Dune height path name
d_file = "/B3D_Inputs/Marsh_Test_Inputs/barrier3d-dunes.npy"
# Storm file path name
#s_file = "/B3D_Inputs/StormTimeSeries_1000_10.npy"
#s_file = "/B3D_Inputs/Default_StormTimeSeries_1000yr.npy"
s_file = "/B3D_Inputs/Altered_Twenty_Five_StormTimeSeries_1000.npy"
c_wd = os.getcwd()
nt_run = 500 # Number of years model will run
run_name = ['Wreck_Marsh_250_S10_1','Wreck_Marsh_250_S10_2','Wreck_Marsh_250_S10_3']
num_of_batches = len(run_name)
rslr_index = [0.0053,.0147,0.026]
#rslr_index = [.0147,0.026]

number_barrier3d_models = 5
rmin = [0.55] * number_barrier3d_models
rmax = [0.95] * number_barrier3d_models
elevation_file = [
    c_wd + "/B3D_Inputs/Wreck_Elevation_6.npy",
    c_wd + "/B3D_Inputs/Wreck_Elevation_5.npy",
    c_wd + "/B3D_Inputs/Wreck_Elevation_4.npy",
    c_wd + "/B3D_Inputs/Wreck_Elevation_3.npy",
    c_wd + "/B3D_Inputs/Wreck_Elevation_2.npy",
]
dune_file = [
    c_wd + "/B3D_Inputs/barrier3d-dunes.npy",
]*5
storm_file = c_wd+s_file

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
    sea_level_rise_rate=0.004,  # not an array
    sea_level_constant=True,  # not an array
    enable_shoreline_offset=False,
    shoreline_offset=[0],
    marsh_dynamics = True,
):

    # ###############################################################################
    # 9 - connect cascade domains (human management) with AST
    # ###############################################################################

    # --------- INITIALIZE ---------
    datadir = "B3D_Inputs/"
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
        community_dynamics_module=False,  # no community dynamics
        enable_shoreline_offset=enable_shoreline_offset,  # Bool
        shoreline_offset=shoreline_offset,
        marsh_dynamics=marsh_dynamics,
    )
    # --------- LOOP ---------
    for time_step in range(nt - 1):

        # Print time step to screen (NOTE: time_index in each model is time_step+1)
        print("\r", "Time Step: ", time_step, end="")
        cascade.update(Time_step = time_step)
        if cascade.b3d_break:
            break

    # --------- SAVE ---------
    save_directory = "Run_Output/"
    cascade.save(save_directory)

    return cascade


for i in range(3):
    Batch_Runs(
        nt=nt_run,
        name=run_name[i],
        storm_file=c_wd + s_file,
        alongshore_section_count=number_barrier3d_models,
        num_cores=3,
        rmin=rmin,
        rmax=rmax,
        elevation_file=elevation_file,
        dune_file=dune_file,
        background_erosion=-1.00,
        sea_level_constant=True,  # not an array
        enable_shoreline_offset=False,
        marsh_dynamics=True,
        sea_level_rise_rate = rslr_index[i],
    )
    os.chdir('/Users/ceclmac/PycharmProjects/CASCADE')


