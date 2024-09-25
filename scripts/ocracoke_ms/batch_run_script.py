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
os.chdir('C:\\Users\\frank\\PycharmProjects\\CASCADE')
# Specify variables to use in calling function

# Dune height path name
e_file = []
d_file = []

for i in range(11,13):
    dune_name = 'C:\\Users\\frank\\PycharmProjects\\CASCADE\\data\\Ocracoke_init_data\\dunes\\Sample_'+str(i)+'_dune.npy'
    elev_name = 'C:\\Users\\frank\\PycharmProjects\\CASCADE\\data\\Ocracoke_init_data\\elevations\\Ocracoke_Revised_Topography_1974_Grid_'+str(i)+'.npy'
    d_file.append(dune_name)
    e_file.append(elev_name)
Num_Storms = 100
c_wd = os.getcwd()
nt_run = 46  # Number of years model will run
run_name = []
s_file = []
for i in range(10,12):
    run_name.append('O_Hindcast_1974_S_'+str(i))
    s_file.append('C:\\Users\\frank\\PycharmProjects\\CASCADE\\data\\Ocracoke_init_data\\Ocracoke_StormList_'+
                       str(i)+'_baseline.npy')
num_of_batches = 1

road_load_name = 'C:\\Users\\frank\\OneDrive - University of North Carolina at Chapel Hill\\Chapter 3\\Revised_Offshore_Datum\\Corrected_Road_Offsets.csv'
dune_load_name = 'C:\\Users\\frank\\OneDrive - University of North Carolina at Chapel Hill\\Chapter 3\\Revised_Offshore_Datum\\Dune_Offsets.csv'


#road_setbacks = np.loadtxt(road_load_name,skiprows=1,delimiter=',')
#dune_offset = np.loadtxt(dune_load_name,skiprows=1,delimiter=',')

road_setbacks =[0,0] #road_setbacks[0,0]*10
dune_offset = [0,0]#dune_offset[0,0]*10

# Start with 1974 shoreline shape
background_threhold_list = [0,0]

# Name storms


# Call function


def alongshore_connected(
    nt,
    name,
    storm_file,
    alongshore_section_count,
    num_cores,
    beach_width_threshold,  # not a parameter in cascade, for triggering: must be list
    rmin,  # the remaining variables are arrays
    rmax,
    elevation_file,
    dune_file,
    dune_design_elevation,
    dune_minimum_elevation,
    road_ele,
    road_width,
    road_setback,
    overwash_filter,
    overwash_to_dune,
    nourishment_volume,
    background_erosion,
    rebuild_dune_threshold,
    roadway_management_on,
    beach_dune_manager_on,
    sea_level_rise_rate=0.0037,  # not an array
    sea_level_constant=True,  # not an array
    trigger_dune_knockdown=False,
    group_roadway_abandonment=None,
    sandbag_management_on = False,
    sandbag_elevation = 5,
    enable_shoreline_offset = False,
    shoreline_offset = [0] ,
):
    # ###############################################################################
    # 9 - connect cascade domains (human management) with AST
    # ###############################################################################

    # --------- INITIALIZE ---------
    datadir = "./data/Ocracoke_init_data/"
    cascade = Cascade(
        datadir,
        name,
        storm_file=storm_file,
        elevation_file=elevation_file,
        dune_file=dune_file,
        parameter_file="Ocracoke-CASCADE-parameters.yaml",
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
        roadway_management_module=roadway_management_on,
        alongshore_transport_module=True,  # couple brie
        beach_nourishment_module=beach_dune_manager_on,
        community_economics_module=False,  # no community dynamics
        road_ele=road_ele,
        road_width=road_width,
        road_setback=road_setback,
        dune_design_elevation=dune_design_elevation,
        dune_minimum_elevation=dune_minimum_elevation,
        trigger_dune_knockdown=trigger_dune_knockdown,
        group_roadway_abandonment=group_roadway_abandonment,
        nourishment_interval=None,  # yrs
        nourishment_volume=nourishment_volume,  # m^3/m
        overwash_filter=overwash_filter,  # % overwash removed
        overwash_to_dune=overwash_to_dune,
        sandbag_management_on = sandbag_management_on,
        sandbag_elevation = sandbag_elevation,
        enable_shoreline_offset=enable_shoreline_offset,
        shoreline_offset=shoreline_offset,
    )

    # --------- LOOP ---------

    # after each year, check the beach width and dune elevation and decide if you want to nourish or rebuild the dune
    # next year with nourish_now parameter; just use first B3D domain, since all berm elevations are equivalent
    dune_rebuild_threshold = rebuild_dune_threshold + (
        cascade.barrier3d[0].BermEl * 10
    )  # if rebuild_dune_threshold=0.3, this is the same threshold for abs. min elevation as in RoadwayManager (m MHW)

    for time_step in range(nt - 1):
        # Print time step to screen (NOTE: time_index in each model is time_step+1)
        print("\r", "Time Step: ", time_step, end="")
        cascade.update()
        if cascade.b3d_break:
            break

        t = cascade.barrier3d[0].time_index
        tmp_rebuild_dune = np.zeros(alongshore_section_count)
        tmp_nourish_now = np.zeros(alongshore_section_count)

        for iB3D in range(alongshore_section_count):
            # don't do any beach/dune management activities if the barrier has become too narrow to sustain a community
            if cascade.community_break[iB3D]:
                pass
            # and only manage beach/dune if it is turned on
            elif beach_dune_manager_on[iB3D]:
                if (
                    cascade.nourishments[iB3D].beach_width[t - 1]
                    < beach_width_threshold[iB3D]
                ):
                    # cascade.nourish_now[iB3D] = 1
                    tmp_nourish_now[iB3D] = 1

                DuneDomainCrest = (
                    cascade.barrier3d[iB3D].DuneDomain[t - 1, :, :].max(axis=1)
                )  # Maximum height of each row in dune domain [dam]
                DuneCrestMin = (
                    np.min(DuneDomainCrest) + cascade.barrier3d[iB3D].BermEl
                ) * 10  # m MHW

                if DuneCrestMin < dune_rebuild_threshold:
                    # cascade.rebuild_dune_now[iB3D] = 1
                    tmp_rebuild_dune[iB3D] = 1

        # only nourish or rebuild dune if all segments fall below threshold (more realistic)
        if np.all(tmp_nourish_now[beach_dune_manager_on]) == 1:
            cascade.nourish_now = tmp_nourish_now
        if np.all(tmp_rebuild_dune[beach_dune_manager_on]) == 1:
            cascade.rebuild_dune_now = tmp_rebuild_dune

    # --------- SAVE ---------
    save_directory = "Run_Output/"
    cascade.save(save_directory)
    os.chdir("..")

    return cascade

for j in range(0,Num_Storms):
    def alongshore_uniform():
        # variables that DO NOT change among runs
        number_barrier3d_models = 2
        beach_width_threshold = [30] * number_barrier3d_models
        rmin = [0.55] * number_barrier3d_models
        rmax = [0.95] * number_barrier3d_models
        elevation_file = e_file  # [e_file] * number_barrier3d_models
        dune_file = d_file  # [d_file] * number_barrier3d_models
        storm_file = s_file[j]
        dune_design_elevation = [2.6] * number_barrier3d_models  # 2 m scenario
        num_cores = 4  # for my laptop, max is 15
        dune_minimum_elevation = 1.1  # m MHW, allow dune to erode down to 0.5 m above the roadway, for roadways only
        road_ele = 1.45  # m MHW (Based on 1997 LIDAR)
        road_width = 20  # m (2 lane road on Ocracoke)
        road_setback = road_setbacks  # m
        overwash_filter = 10  # residental
        overwash_to_dune = 9
        nourishment_volume = 100  # m^3/m
        background_erosion = background_threhold_list  # m/yr, background shoreline erosion
        rebuild_dune_threshold = 1  # m
        sandbag_management_on = [False] * number_barrier3d_models
        sandbag_elevation = 4  # m

        # baseline models for comparison -- all roadways ----------------------------------------
        roads_on = [False] * number_barrier3d_models
        nourishments_on = [False] * number_barrier3d_models
        sea_level_rise_rate = 0.0056
        sea_level_constant = True  # linear

        # Island offsets
        shoreline_offset_enabled = True
        shoreline_offset = dune_offset

        # Island is too narrow for roadway to be relocated. Roadway eaten up by dunes at 73 years
        alongshore_connected(
            nt=nt_run,
            name=run_name[j],
            storm_file=storm_file,
            alongshore_section_count=number_barrier3d_models,
            num_cores=num_cores,
            beach_width_threshold=beach_width_threshold,
            rmin=rmin,
            rmax=rmax,  # rave = 0.75
            elevation_file=elevation_file,
            dune_file=dune_file,
            dune_design_elevation=dune_design_elevation,
            dune_minimum_elevation=dune_minimum_elevation,
            road_ele=road_ele,
            road_width=road_width,
            road_setback=road_setback,
            overwash_filter=overwash_filter,
            overwash_to_dune=overwash_to_dune,
            nourishment_volume=nourishment_volume,
            background_erosion=background_erosion,
            rebuild_dune_threshold=rebuild_dune_threshold,
            roadway_management_on=roads_on,
            beach_dune_manager_on=nourishments_on,
            sea_level_rise_rate=sea_level_rise_rate,
            sea_level_constant=sea_level_constant,
            sandbag_management_on=sandbag_management_on,
            sandbag_elevation=sandbag_elevation,
            shoreline_offset=shoreline_offset,
            enable_shoreline_offset=shoreline_offset_enabled,
        )


    alongshore_uniform()
    os.chdir('C:\\Users\\frank\\PycharmProjects\\CASCADE')
    print('Finished ' + str(run_name[j]))