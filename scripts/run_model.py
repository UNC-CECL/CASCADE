# run file for

# ~******* CASCADE ********~

"""----------------------------------------------------
Copyright (C) 2020 Katherine Anarde
----------------------------------------------------"""

# remember if I move to a different computer to $ pip install -e . in the brie, B3D, and CHOM directories

import numpy as np
import os

from scripts import CASCADE_plotters as CASCADEplt

from cascade.cascade import Cascade

from barrier3d.tools.input_files import (
    yearly_storms,
    gen_dune_height_start,
    gen_alongshore_variable_rmin_rmax,
    shift_storm_intensity,
)
from itertools import compress

# for laptop and desktop, use all but one core; on supercomputer, use all cores; KA Macbook has 15
# num_cores = multiprocessing.cpu_count() - 1

# # ###############################################################################
# # runs
# # ###############################################################################


def RUN_4_CASCADE_noAST_Rave_SLR_pt004_NoHumans(
    nt,
    rmin,
    rmax,
    name,
    storm_file,
    elevation_file,
    dune_file,
):

    # ###############################################################################
    # 4 - CASCADE with only one B3D model and no human dynamics
    # ###############################################################################
    # Use the starting interior domain from the 10,000 yr runs for each dune growth rate and run for 1000 years
    # or until the barrier drowns. All other modules (brie and human dymnamics modules) turned off. Can also use this
    # run script for the 10,000 year runs.

    # --------- INITIALIZE ---------
    datadir = "B3D_Inputs/"
    # datadir = "/Users/KatherineAnardeWheels/PycharmProjects/CASCADE/B3D_Inputs/"
    cascade = Cascade(
        datadir,
        name,
        storm_file=storm_file,
        elevation_file=elevation_file,
        dune_file=dune_file,
        parameter_file="RUN4-CASCADE-parameters.yaml",
        wave_height=1,
        wave_period=7,  # s (lowered from 10 s to reduce k_sf)
        wave_asymmetry=0.8,  # fraction approaching from left
        wave_angle_high_fraction=0.2,  # fraction of waves approaching from higher than 45 degrees
        sea_level_rise_rate=0.004,  # m/yr
        sea_level_rise_constant=True,  # linear SLR
        background_erosion=0.0,
        alongshore_section_count=1,  # only one B3D domain
        time_step_count=nt,
        min_dune_growth_rate=rmin,
        max_dune_growth_rate=rmax,
        num_cores=1,
        roadway_management_module=False,  # no roadway management
        alongshore_transport_module=False,  # no brie coupling
        beach_nourishment_module=False,  # no beach nourishment
        community_dynamics_module=False,  # no community dynamics
    )

    # --------- LOOP ---------

    for time_step in range(nt - 1):
        # Print time step to screen (NOTE: time_index in each model is time_step+1)
        print("\r", "Time Step: ", time_step, end="")
        cascade.update()
        if cascade.b3d_break:
            break

    # --------- SAVE ---------
    save_directory = "Run_Output/"
    cascade.save(save_directory)
    os.chdir("..")

    return cascade


def RUN_6_CASCADE_noAST_Rave_SLR_pt004_Roadways(
    nt,
    rmin,
    rmax,
    name,
    storm_file,
    elevation_file,
    dune_file,
    road_ele=1.7,
    road_width=30,
    road_setback=30,
    dune_design_elevation=3.7,
    dune_minimum_elevation=2.2,
    percent_water_cells_sensitivity=None,
    background_erosion=0.0,
):

    # ###############################################################################
    # 6 - CASCADE with only one B3D model and roadway management
    # ###############################################################################
    # Use the starting interior domain from the 10,000 yr runs for each dune growth rate and run for 1000 years
    # or until the barrier drowns. All other modules (beach nourishment, community dyanmics) turned off.

    # --------- INITIALIZE ---------
    datadir = "B3D_Inputs/"
    cascade = Cascade(
        datadir,
        name,
        storm_file=storm_file,
        elevation_file=elevation_file,
        dune_file=dune_file,
        parameter_file="RUN6-CASCADE-parameters.yaml",
        wave_height=1,
        wave_period=7,
        wave_asymmetry=0.8,
        wave_angle_high_fraction=0.2,
        sea_level_rise_rate=0.004,
        sea_level_rise_constant=True,
        background_erosion=background_erosion,
        alongshore_section_count=1,
        time_step_count=nt,
        min_dune_growth_rate=rmin,
        max_dune_growth_rate=rmax,
        num_cores=1,
        roadway_management_module=True,
        alongshore_transport_module=False,  # no brie coupling
        beach_nourishment_module=False,  # no beach nourishment
        community_dynamics_module=False,  # no community dynamics
        road_ele=road_ele,
        road_width=road_width,
        road_setback=road_setback,
        dune_design_elevation=dune_design_elevation,
        dune_minimum_elevation=dune_minimum_elevation,
    )

    # for sensitivity testing
    if percent_water_cells_sensitivity is not None:
        cascade.roadways[
            0
        ].percent_water_cells_touching_road = percent_water_cells_sensitivity

    # --------- LOOP ---------

    for time_step in range(nt - 1):
        # Print time step to screen (NOTE: time_index in each model is time_step+1)
        print("\r", "Time Step: ", time_step, end="")
        cascade.update()
        if cascade.b3d_break:
            break

    # --------- SAVE ---------
    save_directory = "Run_Output/"
    cascade.save(save_directory)
    os.chdir("..")

    return cascade


def RUN_7_CASCADE_noAST_Rave_variableSLR_NoHumans(
    nt,
    rmin,
    rmax,
    name,
    storm_file,
    elevation_file,
    dune_file,
    sea_level_rise_rate,
    sea_level_constant,
):

    # ###############################################################################
    # 7 - same as RUN 4 but with variable rates of SLR (i.e., no AST, no human dynamics)
    # ###############################################################################

    # --------- INITIALIZE ---------
    datadir = "B3D_Inputs/"
    cascade = Cascade(
        datadir,
        name,
        storm_file=storm_file,
        elevation_file=elevation_file,
        dune_file=dune_file,
        parameter_file="RUN7-CASCADE-parameters.yaml",
        wave_height=1,
        wave_period=7,
        wave_asymmetry=0.8,
        wave_angle_high_fraction=0.2,
        sea_level_rise_rate=sea_level_rise_rate,
        sea_level_rise_constant=sea_level_constant,
        background_erosion=0.0,
        alongshore_section_count=1,
        time_step_count=nt,
        min_dune_growth_rate=rmin,
        max_dune_growth_rate=rmax,
        num_cores=1,
        roadway_management_module=False,  # no roadway dynamics
        alongshore_transport_module=False,  # no brie coupling
        beach_nourishment_module=False,  # no beach nourishment
        community_dynamics_module=False,  # no community dynamics
    )

    # --------- LOOP ---------

    for time_step in range(nt - 1):
        # Print time step to screen (NOTE: time_index in each model is time_step+1)
        print("\r", "Time Step: ", time_step, end="")
        cascade.update()
        if cascade.b3d_break:
            break

    # --------- SAVE ---------
    save_directory = "Run_Output/"
    cascade.save(save_directory)
    os.chdir("..")

    return cascade


def RUN_8_CASCADE_noAST_Rave_SLR_pt004_Nourishment(
    nt,
    rmin,
    rmax,
    name,
    dune_design_elevation,
    storm_file,
    elevation_file,
    dune_file,
    overwash_filter,
    overwash_to_dune,
    nourishment_volume,
    beach_width_threshold,
    background_erosion,
    rebuild_dune_threshold,
):

    # ###############################################################################
    # 8 - nourish beach, rebuild dunes, and remove overwash from barrier interior
    # ###############################################################################

    # --------- INITIALIZE ---------
    datadir = "B3D_Inputs/"
    cascade = Cascade(
        datadir,
        name,
        storm_file=storm_file,
        elevation_file=elevation_file,
        dune_file=dune_file,
        parameter_file="RUN8-CASCADE-parameters.yaml",
        wave_height=1,
        wave_period=7,
        wave_asymmetry=0.8,
        wave_angle_high_fraction=0.2,
        sea_level_rise_rate=0.004,
        sea_level_rise_constant=True,
        background_erosion=background_erosion,
        alongshore_section_count=1,
        time_step_count=nt,
        min_dune_growth_rate=rmin,
        max_dune_growth_rate=rmax,
        num_cores=1,
        roadway_management_module=False,  # no roadway management
        alongshore_transport_module=False,  # no brie coupling
        beach_nourishment_module=True,
        community_dynamics_module=False,  # no community dynamics
        dune_design_elevation=dune_design_elevation,
        nourishment_interval=None,  # yrs
        nourishment_volume=nourishment_volume,  # m^3/m
        overwash_filter=overwash_filter,  # % overwash filtered by development
        overwash_to_dune=overwash_to_dune,  # % overwash bulldozed back to dune
    )

    # --------- LOOP ---------

    iB3D = 0  # we only have one Barrier3D domain here

    # after each year, check the beach width and dune elevation and decide if you want to nourish or rebuild the dune
    # next year with nourish_now parameter
    dune_rebuild_threshold = rebuild_dune_threshold + (
        cascade.barrier3d[iB3D].BermEl * 10
    )  # if rebuild_dune_threshold=0.3, this is the same threshold for abs. min elevation as in RoadwayManager (m MHW)

    for time_step in range(nt - 1):

        # Print time step to screen (NOTE: time_index in each model is time_step+1)
        print("\r", "Time Step: ", time_step, end="")
        cascade.update()
        if cascade.b3d_break:
            break

        # stop managing if the barrier becomes too narrow to sustain a community
        if cascade.community_break[iB3D]:
            pass
        else:
            t = cascade.barrier3d[iB3D].time_index

            if cascade.nourishments[iB3D].beach_width[t - 1] < beach_width_threshold:
                cascade.nourish_now[iB3D] = 1

            DuneDomainCrest = (
                cascade.barrier3d[iB3D].DuneDomain[t - 1, :, :].max(axis=1)
            )  # Maximum height of each row in dune domain [dam]
            # DuneRestart = cascade.barrier3d[iB3D].DuneRestart
            # DuneDomainCrest[DuneDomainCrest < DuneRestart] = DuneRestart
            DuneCrestMin = (
                np.min(DuneDomainCrest) + cascade.barrier3d[iB3D].BermEl
            ) * 10  # m MHW

            if DuneCrestMin < dune_rebuild_threshold:
                cascade.rebuild_dune_now[iB3D] = 1

    # --------- SAVE ---------
    save_directory = "Run_Output/"
    cascade.save(save_directory)
    os.chdir("..")

    return cascade


def RUN_9_CASCADE_Rave_SLR_pt004_AlongshoreVariableManagement(
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
    sea_level_rise_rate=0.004,  # not an array
    sea_level_constant=True,  # not an array
    trigger_dune_knockdown=False,
    group_roadway_abandonment=None,
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
        parameter_file="RUN9-CASCADE-parameters.yaml",
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
        community_dynamics_module=False,  # no community dynamics
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


# # ###############################################################################
# # plotters
# # ###############################################################################


def PLOT_5_Nonlinear_Dynamics_CASCADE_B3Donly_RoadwayManager(
    name_prefix,
    tmax_roadways,
    tmax_sim,
    plot_name,
    run_road_mgmt,
    cross_sections=None,
    text_out=False,
    gif_on=False,
):

    os.chdir(
        "/Users/KatherineAnardeWheels/Research/BARis/UNC/CNH/CASCADE_save_dir/Run_Output"
    )

    # --------- plot ---------
    output = np.load(name_prefix + ".npz", allow_pickle=True)
    cascade = output["cascade"]
    cascade = cascade[0]
    b3d = cascade.barrier3d

    # individual dune growth rates
    ib3d = 0  # this plotter is just for one B3D domain, no alongshore grid cells

    if run_road_mgmt:
        post_storm_dunes = cascade.roadways[ib3d]._post_storm_dunes
        post_storm_ave_interior_height = cascade.roadways[
            ib3d
        ]._post_storm_ave_interior_height
        design_height = cascade.roadways[ib3d]._dune_design_elevation_TS
        rebuild_threshold = cascade.roadways[ib3d]._dune_minimum_elevation_TS
        road_elevation = cascade.roadways[ib3d]._road_ele_TS
        dunes_rebuilt = cascade.roadways[ib3d]._dunes_rebuilt_TS
        road_relocated = cascade.roadways[ib3d]._road_relocated_TS
    else:
        post_storm_dunes = None
        post_storm_ave_interior_height = None
        design_height = None
        rebuild_threshold = None
        road_elevation = None
        dunes_rebuilt = None
        road_relocated = None

    (
        BarrierWidth,
        DuneCrestMean,
        BarrierHeight,
        bh_rate,
        bw_rate,
        sc_rate,
        DuneCrestMin,
        DuneCrestMax,
        shoreline_position,
        shoreface_slope,
        overwash,
    ) = CASCADEplt.plot_nonlinear_stats_RoadwayManager(
        b3d,
        ib3d,
        tmax_roadways=tmax_roadways,
        tmax_sim=tmax_sim,
        post_storm_dunes=post_storm_dunes,
        post_storm_ave_interior_height=post_storm_ave_interior_height,
        design_height=design_height,
        rebuild_threshold=rebuild_threshold,
        road_elevation=road_elevation,
        dunes_rebuilt=dunes_rebuilt,
        road_relocated=road_relocated,
    )

    if text_out:
        # save to text file for use in Matlab
        np.savetxt(
            name_prefix + ".txt",
            (
                np.array(BarrierWidth),
                np.array(BarrierHeight),
                np.array(DuneCrestMean),
            ),
        )

    if gif_on:
        # also make the gif
        directory = "/Users/KatherineAnardeWheels/PycharmProjects/CASCADE/"
        if cascade.roadways is not None:
            CASCADEplt.plot_ElevAnimation_CASCADE(
                cascade,
                ny=1,
                directory=directory,
                TMAX_MGMT=[tmax_roadways],
                name=plot_name,
                TMAX_SIM=tmax_sim,
                beach_management_ny=[0],
                roadway_management_ny=[1],
                z_lim=4,
                # y_lim=[150, 220],
                fig_size=(6, 6),
                fig_eps=False,
            )
        else:
            CASCADEplt.plot_ElevAnimation_CASCADE(
                cascade,
                ny=1,
                directory=directory,
                TMAX_MGMT=[0],
                name=plot_name,
                TMAX_SIM=tmax_sim,
                beach_management_ny=[0],
                roadway_management_ny=[0],
                z_lim=4,
                # y_lim=[150, 220],
                fig_size=(6, 6),
                fig_eps=False,
            )

    if cross_sections is not None:
        time_step = cross_sections
        fig = CASCADEplt.plot_ModelTransects(cascade, time_step, iB3D=0)
        # fig.set_title("RoadwayManager")

    return (
        BarrierWidth,
        DuneCrestMean,
        BarrierHeight,
        bh_rate,
        bw_rate,
        sc_rate,
        DuneCrestMin,
        DuneCrestMax,
        shoreline_position,
        shoreface_slope,
        overwash,
        cascade,
    )


def PLOT_6_Nonlinear_Dynamics_CASCADE_B3Donly_Nourishments(
    name_prefix,
    tmax_management,
    tmax_sim,
    plot_name,
    rebuild_dune_threshold,
    cross_sections=None,
    text_out=False,
    gif_on=False,
):

    os.chdir(
        "/Users/KatherineAnardeWheels/Research/BARis/UNC/CNH/CASCADE_save_dir/Run_Output"
    )

    # --------- plot ---------
    output = np.load(name_prefix + ".npz", allow_pickle=True)
    cascade = output["cascade"]
    cascade = cascade[0]
    b3d = cascade.barrier3d

    ib3d = 0  # this is just B3D, so no alongshore grid cells
    rebuild_threshold = rebuild_dune_threshold + (
        b3d[ib3d].BermEl * 10
    )  # min dune height above the berm [m MHW]

    (
        BarrierWidth,
        DuneCrestMean,
        BarrierHeight,
        bh_rate,
        bw_rate,
        sc_rate,
        DuneCrestMin,
        DuneCrestMax,
        shoreline_position,
        shoreface_slope,
        beach_width,
        overwash,
        dune_toe,
    ) = CASCADEplt.plot_nonlinear_stats_BeachDuneManager(
        b3d,
        ib3d,
        tmax_management=tmax_management,
        tmax_sim=tmax_sim,
        nourishments=cascade.nourishments,
        post_storm_dunes=cascade.nourishments[ib3d]._post_storm_dunes,
        post_storm_x_s=cascade.nourishments[ib3d]._post_storm_x_s,
        post_storm_s_sf=cascade.nourishments[ib3d]._post_storm_s_sf,
        post_storm_ave_interior_width=cascade.nourishments[
            ib3d
        ]._post_storm_ave_interior_width,
        post_storm_ave_interior_height=cascade.nourishments[
            ib3d
        ]._post_storm_ave_interior_height,
        post_storm_beach_width=cascade.nourishments[ib3d]._post_storm_beach_width,
        post_storm_Qow=cascade.nourishments[ib3d]._post_storm_Qow,
        design_elevation=cascade.nourishments[ib3d]._dune_design_elevation,  # m MHW,
        rebuild_threshold=rebuild_threshold,
        dunes_rebuilt=cascade.nourishments[ib3d]._dunes_rebuilt_TS,
    )

    if gif_on:
        directory = "/Users/KatherineAnardeWheels/PycharmProjects/CASCADE/"
        CASCADEplt.plot_ElevAnimation_CASCADE(
            cascade,
            ny=1,
            directory=directory,
            TMAX_MGMT=[tmax_management],
            name=plot_name,
            TMAX_SIM=tmax_sim,
            beach_management_ny=[1],
            roadway_management_ny=[0],
            z_lim=4,
            # y_lim=[150, 220],
            fig_size=(6, 6),
            fig_eps=False,
        )

    if cross_sections is not None:
        time_step = cross_sections
        fig = CASCADEplt.plot_ModelTransects(cascade, time_step, iB3D=0)
        # fig.set_title("BeachDuneManager")

    return (
        BarrierWidth,
        DuneCrestMean,
        BarrierHeight,
        bh_rate,
        bw_rate,
        sc_rate,
        DuneCrestMin,
        DuneCrestMax,
        shoreline_position,
        shoreface_slope,
        beach_width,
        overwash,
        dune_toe,
        cascade,
    )


def PLOT_7_Initial_CNH_Topographies(name_prefix_list):

    os.chdir(
        "/Users/KatherineAnardeWheels/Research/BARis/UNC/CNH/CASCADE_save_dir/Run_Output"
    )

    cascade = []

    for i in range(0, len(name_prefix_list)):
        output = np.load(name_prefix_list[i] + ".npz", allow_pickle=True)
        csc8d = output["cascade"]
        cascade.append(csc8d[0])

    CASCADEplt.fig2_initialCNH_topo(cascade)

    return


def PLOT_9_Nonlinear_Dynamics_CASCADE_AST(
    name_prefix,
    tmax_management,  # an array
    tmax_sim,  # not an array
    plot_name,
    rebuild_dune_threshold,
    beach_management_ny=None,
    roadway_management_ny=None,
    gif_on=False,
    y_lim=None,
    z_lim=3.5,
    fig_size=None,
    time_series_on=True,
    fig_eps=False,
):

    os.chdir(
        "/Users/KatherineAnardeWheels/Research/BARis/UNC/CNH/CASCADE_save_dir/Run_Output"
    )

    # --------- plot ---------
    output = np.load(name_prefix + ".npz", allow_pickle=True)
    cascade = output["cascade"]
    cascade = cascade[0]
    b3d = cascade.barrier3d
    ny = np.size(b3d)

    (
        barrier_width,
        dune_crest_mean,
        barrier_height,
        bh_rate,
        bw_rate,
        sc_rate,
        dune_crest_min,
        dune_crest_max,
        shoreline_position,
        shoreface_slope,
        beach_width,
        overwash,
        dune_toe,
    ) = ([] for i in range(13))

    if time_series_on:
        for ib3d in range(ny):
            if beach_management_ny[ib3d]:

                rebuild_threshold = rebuild_dune_threshold + (
                    b3d[ib3d].BermEl * 10
                )  # min dune height above the berm [m MHW]

                (
                    BarrierWidth,
                    DuneCrestMean,
                    BarrierHeight,
                    BarrierHeightRate,
                    BarrierWidthRate,
                    ShorelineChangeRate,
                    DuneCrestMin,
                    DuneCrestMax,
                    ShorelinePosition,
                    ShorefaceSlope,
                    BeachWidth,
                    Qow,
                    DuneToePosition,
                ) = CASCADEplt.plot_nonlinear_stats_BeachDuneManager(
                    b3d,
                    ib3d,
                    tmax_management=tmax_management[ib3d],
                    tmax_sim=tmax_sim,
                    nourishments=cascade.nourishments,
                    post_storm_dunes=cascade.nourishments[ib3d]._post_storm_dunes,
                    post_storm_x_s=cascade.nourishments[ib3d]._post_storm_x_s,
                    post_storm_s_sf=cascade.nourishments[ib3d]._post_storm_s_sf,
                    post_storm_ave_interior_width=cascade.nourishments[
                        ib3d
                    ]._post_storm_ave_interior_width,
                    post_storm_ave_interior_height=cascade.nourishments[
                        ib3d
                    ]._post_storm_ave_interior_height,
                    post_storm_beach_width=cascade.nourishments[
                        ib3d
                    ]._post_storm_beach_width,
                    post_storm_Qow=cascade.nourishments[ib3d]._post_storm_Qow,
                    design_elevation=cascade.nourishments[
                        ib3d
                    ]._dune_design_elevation,  # m MHW,
                    rebuild_threshold=rebuild_threshold,
                    dunes_rebuilt=cascade.nourishments[ib3d]._dunes_rebuilt_TS,
                )

            elif roadway_management_ny[ib3d]:
                post_storm_dunes = cascade.roadways[ib3d]._post_storm_dunes
                post_storm_ave_interior_height = cascade.roadways[
                    ib3d
                ]._post_storm_ave_interior_height
                design_height = cascade.roadways[ib3d]._dune_design_elevation_TS
                rebuild_threshold = cascade.roadways[ib3d]._dune_minimum_elevation_TS
                road_elevation = cascade.roadways[ib3d]._road_ele_TS
                dunes_rebuilt = cascade.roadways[ib3d]._dunes_rebuilt_TS
                road_relocated = cascade.roadways[ib3d]._road_relocated_TS

                (
                    BarrierWidth,
                    DuneCrestMean,
                    BarrierHeight,
                    BarrierHeightRate,
                    BarrierWidthRate,
                    ShorelineChangeRate,
                    DuneCrestMin,
                    DuneCrestMax,
                    ShorelinePosition,
                    ShorefaceSlope,
                    Qow,
                ) = CASCADEplt.plot_nonlinear_stats_RoadwayManager(
                    b3d,
                    ib3d,
                    tmax_roadways=tmax_management[ib3d],
                    tmax_sim=tmax_sim,
                    post_storm_dunes=post_storm_dunes,
                    post_storm_ave_interior_height=post_storm_ave_interior_height,
                    design_height=design_height,
                    rebuild_threshold=rebuild_threshold,
                    road_elevation=road_elevation,
                    dunes_rebuilt=dunes_rebuilt,
                    road_relocated=road_relocated,
                )

                BeachWidth = [0]  # dummy
                DuneToePosition = [0]  # dummy

            else:

                post_storm_dunes = None
                post_storm_ave_interior_height = None
                design_height = None
                rebuild_threshold = None
                road_elevation = None
                dunes_rebuilt = None
                road_relocated = None

                (
                    BarrierWidth,
                    DuneCrestMean,
                    BarrierHeight,
                    BarrierHeightRate,
                    BarrierWidthRate,
                    ShorelineChangeRate,
                    DuneCrestMin,
                    DuneCrestMax,
                    ShorelinePosition,
                    ShorefaceSlope,
                    Qow,
                ) = CASCADEplt.plot_nonlinear_stats_RoadwayManager(
                    b3d,
                    ib3d,
                    tmax_roadways=tmax_management[ib3d],
                    tmax_sim=tmax_sim,
                    post_storm_dunes=post_storm_dunes,
                    post_storm_ave_interior_height=post_storm_ave_interior_height,
                    design_height=design_height,
                    rebuild_threshold=rebuild_threshold,
                    road_elevation=road_elevation,
                    dunes_rebuilt=dunes_rebuilt,
                    road_relocated=road_relocated,
                )

                BeachWidth = [0]  # dummy
                DuneToePosition = [0]  # dummy

            barrier_width.append(BarrierWidth),
            dune_crest_mean.append(DuneCrestMean),
            barrier_height.append(BarrierHeight),
            bh_rate.append(BarrierHeightRate),
            bw_rate.append(BarrierWidthRate),
            sc_rate.append(ShorelineChangeRate),
            dune_crest_min.append(DuneCrestMin),
            dune_crest_max.append(DuneCrestMax),
            shoreline_position.append(ShorelinePosition),
            shoreface_slope.append(ShorefaceSlope),
            beach_width.append(BeachWidth),
            overwash.append(Qow),
            dune_toe.append(DuneToePosition),

    if gif_on:
        directory = "/Users/KatherineAnardeWheels/PycharmProjects/CASCADE/"
        CASCADEplt.plot_ElevAnimation_CASCADE(
            cascade,
            ny=ny,
            directory=directory,
            TMAX_MGMT=tmax_management,  # an array
            name=plot_name,
            TMAX_SIM=tmax_sim,  # not an array
            beach_management_ny=beach_management_ny,
            roadway_management_ny=roadway_management_ny,
            y_lim=y_lim,
            z_lim=z_lim,
            fig_size=fig_size,
            fig_eps=fig_eps,
        )

    return (
        barrier_width,
        dune_crest_mean,
        barrier_height,
        bh_rate,
        bw_rate,
        sc_rate,
        dune_crest_min,
        dune_crest_max,
        shoreline_position,
        shoreface_slope,
        beach_width,
        overwash,
        dune_toe,
        cascade,
    )


# # ###############################################################################
# # statistics for paper
# # ###############################################################################


def get_roadway_statistics(
    folder_prefix,
    natural_barrier_elev=None,
    natural_barrier_width=None,
    individual_fid=None,
    tmax=None,
    iB3D=0,
):

    # folder_prefix = "Roadway_100sims_1m_lowGR_lowEle"
    folder_path = (
        "/Users/KatherineAnardeWheels/Research/BARis/UNC/CNH/CASCADE_save_dir/Run_Output/"
        + folder_prefix
    )
    os.chdir(folder_path)

    year_abandoned = []
    sim_max = []  # total length of simulation
    road_bulldozed = []  # overwash removal -- number of times?
    overwash_removed = []  # total m^3 of overwash
    dune_rebuilt = []  # rebuild dune -- total number
    road_relocated = []  # road relocation -- total number
    diff_barrier_width = []  # diff between natural barrier height and width at end
    diff_barrier_elev = []

    if individual_fid is not None:
        output = np.load(individual_fid + ".npz", allow_pickle=True)
        cascade = output["cascade"]
        cascade = cascade[0]
        b3d = cascade.barrier3d[iB3D]

        tmax_sim = b3d.time_index - 1
        sim_max.append(tmax_sim)

        # if user specified a window for statistics, then just set tmax to that window
        if tmax is not None:
            tmax_sim = tmax

        final_barrier_width = (
            np.array(b3d.x_b_TS[tmax_sim - 1]) - np.array(b3d.x_s_TS[tmax_sim - 1])
        ) * 10  # m
        final_domain = np.array(b3d.DomainTS[tmax_sim - 1]) * 10
        final_barrier_elev = final_domain[final_domain > 0].mean()  # m MHW
        if natural_barrier_width is not None:
            diff_barrier_width.append(natural_barrier_width - final_barrier_width)
        if natural_barrier_elev is not None:
            diff_barrier_elev.append(natural_barrier_elev - final_barrier_elev)

        year_abandoned.append(cascade.roadways[iB3D]._time_index - 1)
        road_bulldozed.append(
            sum(cascade.roadways[iB3D]._road_overwash_volume[0:tmax_sim] > 0)
        )  # only counts true elements
        overwash_removed.append(
            sum(cascade.roadways[iB3D]._road_overwash_volume[0:tmax_sim])
        )  # m^3
        dune_rebuilt.append(
            int(sum(cascade.roadways[iB3D]._dunes_rebuilt_TS[0:tmax_sim]))
        )
        road_relocated.append(
            int(sum(cascade.roadways[iB3D]._road_relocated_TS[0:tmax_sim]))
        )

    else:
        for filenum in range(100):
            output = np.load(folder_prefix + str(filenum) + ".npz", allow_pickle=True)
            cascade = output["cascade"]
            cascade = cascade[0]
            b3d = cascade.barrier3d[iB3D]

            tmax_sim = b3d.time_index - 1
            sim_max.append(tmax_sim)

            # if user specified a window for statistics, then just set tmax to that window
            if tmax is not None:
                tmax_sim = tmax

            final_barrier_width = (
                np.array(b3d.x_b_TS[tmax_sim - 1]) - np.array(b3d.x_s_TS[tmax_sim - 1])
            ) * 10  # m
            final_domain = np.array(b3d.DomainTS[tmax_sim - 1]) * 10
            final_barrier_elev = final_domain[final_domain > 0].mean()  # m MHW
            if natural_barrier_width is not None:
                diff_barrier_width.append(natural_barrier_width - final_barrier_width)
            if natural_barrier_elev is not None:
                diff_barrier_elev.append(natural_barrier_elev - final_barrier_elev)
            year_abandoned.append(cascade.roadways[iB3D]._time_index - 1)
            road_bulldozed.append(
                sum(cascade.roadways[iB3D]._road_overwash_volume[0:tmax_sim] > 0)
            )  # only counts true elements
            overwash_removed.append(
                sum(cascade.roadways[iB3D]._road_overwash_volume[0:tmax_sim])
            )  # m^3
            dune_rebuilt.append(
                int(sum(cascade.roadways[iB3D]._dunes_rebuilt_TS[0:tmax_sim]))
            )
            road_relocated.append(
                int(sum(cascade.roadways[iB3D]._road_relocated_TS[0:tmax_sim]))
            )

    return (
        year_abandoned,
        sim_max,
        road_bulldozed,
        overwash_removed,
        dune_rebuilt,
        road_relocated,
        diff_barrier_width,
        diff_barrier_elev,
    )


def get_nourishment_statistics(
    folder_prefix,
    natural_barrier_elev=None,
    natural_barrier_width=None,
    individual_fid=None,
    tmax=None,
    iB3D=0,
):

    folder_path = (
        "/Users/KatherineAnardeWheels/Research/BARis/UNC/CNH/CASCADE_save_dir/Run_Output/"
        + folder_prefix
    )
    os.chdir(folder_path)

    year_abandoned = []
    sim_max = []  # total length of simulation
    dune_rebuilt = []  # rebuild dune -- total number
    overwash_filtered_removed = []  # m^3
    beach_nourished = []  # beach nourishments -- total number
    diff_barrier_width = []  # diff between natural barrier height and width at end
    diff_barrier_elev = []

    if individual_fid is not None:
        output = np.load(individual_fid + ".npz", allow_pickle=True)
        cascade = output["cascade"]
        cascade = cascade[0]
        b3d = cascade.barrier3d[iB3D]

        tmax_sim = b3d.time_index - 1
        sim_max.append(tmax_sim)

        # if user specified a window for statistics, then just set tmax to that window
        if tmax is not None:
            tmax_sim = tmax

        final_barrier_width = (
            np.array(b3d.x_b_TS[tmax_sim - 1]) - np.array(b3d.x_s_TS[tmax_sim - 1])
        ) * 10  # m
        final_domain = np.array(b3d.DomainTS[tmax_sim - 1]) * 10
        final_barrier_elev = final_domain[final_domain > 0].mean()  # m MHW
        if natural_barrier_width is not None:
            diff_barrier_width.append(natural_barrier_width - final_barrier_width)
        if natural_barrier_elev is not None:
            diff_barrier_elev.append(natural_barrier_elev - final_barrier_elev)
        year_abandoned.append(cascade.nourishments[iB3D]._time_index - 1)
        overwash_filtered_removed.append(
            sum(cascade.nourishments[iB3D]._overwash_volume_removed[0:tmax_sim])
        )  # m^3
        dune_rebuilt.append(
            int(sum(cascade.nourishments[iB3D]._dunes_rebuilt_TS[0:tmax_sim]))
        )
        beach_nourished.append(
            int(sum(cascade.nourishments[iB3D]._nourishment_TS[0:tmax_sim]))
        )

    else:
        for filenum in range(100):
            output = np.load(folder_prefix + str(filenum) + ".npz", allow_pickle=True)
            cascade = output["cascade"]
            cascade = cascade[0]
            b3d = cascade.barrier3d[iB3D]

            tmax_sim = b3d.time_index - 1
            sim_max.append(tmax_sim)

            # if user specified a window for statistics, then just set tmax to that window
            if tmax is not None:
                tmax_sim = tmax

            final_barrier_width = (
                np.array(b3d.x_b_TS[tmax_sim - 1]) - np.array(b3d.x_s_TS[tmax_sim - 1])
            ) * 10  # m
            final_domain = np.array(b3d.DomainTS[tmax_sim - 1]) * 10
            final_barrier_elev = final_domain[final_domain > 0].mean()  # m MHW
            if natural_barrier_width is not None:
                diff_barrier_width.append(natural_barrier_width - final_barrier_width)
            if natural_barrier_elev is not None:
                diff_barrier_elev.append(natural_barrier_elev - final_barrier_elev)
            year_abandoned.append(cascade.nourishments[iB3D]._time_index - 1)
            overwash_filtered_removed.append(
                sum(cascade.nourishments[iB3D]._overwash_volume_removed[0:tmax_sim])
            )  # m^3
            dune_rebuilt.append(
                int(sum(cascade.nourishments[iB3D]._dunes_rebuilt_TS[0:tmax_sim]))
            )
            beach_nourished.append(
                int(sum(cascade.nourishments[iB3D]._nourishment_TS[0:tmax_sim]))
            )

    return (
        year_abandoned,
        sim_max,
        overwash_filtered_removed,
        dune_rebuilt,
        beach_nourished,
        diff_barrier_width,
        diff_barrier_elev,
    )


# # ###############################################################################
# # record of runs and plots
# # ###############################################################################

# 10,000 year simulations -- only one B3D model (AST off) -------------------------------------------------------


def cascade_10kyr_sensitivity():

    cascade_10kyr_pt45_01 = RUN_4_CASCADE_noAST_Rave_SLR_pt004_NoHumans(
        nt=10000,
        rmin=0.25,  # rave = 0.45 (but not 0.5 spaced like in Reeves et al., 2021 -- arbitrary)
        rmax=0.65,
        name="4-B3D_noAST_Rave_pt45_SLR_pt004_10k-yrs_01",
        storm_file="StormSeries_10kyrs_VCR_Berm1pt9m_Slope0pt04_01.npy",
        elevation_file="InitElevHog.npy",
        dune_file="DuneStart_1000dam.npy",
    )

    cascade_10kyr_pt45_02 = RUN_4_CASCADE_noAST_Rave_SLR_pt004_NoHumans(
        nt=10000,
        rmin=0.25,  # rave = 0.45
        rmax=0.65,
        name="4-B3D_noAST_Rave_pt45_SLR_pt004_10k-yrs_02",
        storm_file="StormSeries_10kyrs_VCR_Berm1pt9m_Slope0pt04_02.npy",
        elevation_file="InitElevHog.npy",
        dune_file="DuneStart_1000dam.npy",
    )

    cascade_10kyr_pt45_03 = RUN_4_CASCADE_noAST_Rave_SLR_pt004_NoHumans(
        nt=10000,
        rmin=0.25,  # rave = 0.45
        rmax=0.65,
        name="4-B3D_noAST_Rave_pt45_SLR_pt004_10k-yrs_03",
        storm_file="StormSeries_10kyrs_VCR_Berm1pt9m_Slope0pt04_03.npy",
        elevation_file="InitElevHog.npy",
        dune_file="DuneStart_1000dam.npy",
    )

    cascade_10kyr_pt45_04 = RUN_4_CASCADE_noAST_Rave_SLR_pt004_NoHumans(
        nt=10000,
        rmin=0.25,  # rave = 0.45
        rmax=0.65,
        name="4-B3D_noAST_Rave_pt45_SLR_pt004_10k-yrs_04",
        storm_file="StormSeries_10kyrs_VCR_Berm1pt9m_Slope0pt04_04.npy",
        elevation_file="InitElevHog.npy",
        dune_file="DuneStart_1000dam.npy",
    )

    cascade_10kyr_pt45_05 = RUN_4_CASCADE_noAST_Rave_SLR_pt004_NoHumans(
        nt=10000,
        rmin=0.25,  # rave = 0.45
        rmax=0.65,
        name="4-B3D_noAST_Rave_pt45_SLR_pt004_10k-yrs_05",
        storm_file="StormSeries_10kyrs_VCR_Berm1pt9m_Slope0pt04_05.npy",
        elevation_file="InitElevHog.npy",
        dune_file="DuneStart_1000dam.npy",
    )

    cascade_10kyr_pt75_01 = RUN_4_CASCADE_noAST_Rave_SLR_pt004_NoHumans(
        nt=10000,
        rmin=0.55,  # rave = 0.75
        rmax=0.95,
        name="4-B3D_noAST_Rave_pt75_SLR_pt004_10k-yrs_01",
        storm_file="StormSeries_10kyrs_VCR_Berm1pt9m_Slope0pt04_01.npy",
        elevation_file="InitElevHog.npy",
        dune_file="DuneStart_1000dam.npy",
    )

    cascade_10kyr_pt75_02 = RUN_4_CASCADE_noAST_Rave_SLR_pt004_NoHumans(
        nt=10000,
        rmin=0.55,  # rave = 0.75
        rmax=0.95,
        name="4-B3D_noAST_Rave_pt75_SLR_pt004_10k-yrs_02",
        storm_file="StormSeries_10kyrs_VCR_Berm1pt9m_Slope0pt04_02.npy",
        elevation_file="InitElevHog.npy",
        dune_file="DuneStart_1000dam.npy",
    )

    cascade_10kyr_pt75_03 = RUN_4_CASCADE_noAST_Rave_SLR_pt004_NoHumans(
        nt=10000,
        rmin=0.55,  # rave = 0.75
        rmax=0.95,
        name="4-B3D_noAST_Rave_pt75_SLR_pt004_10k-yrs_03",
        storm_file="StormSeries_10kyrs_VCR_Berm1pt9m_Slope0pt04_03.npy",
        elevation_file="InitElevHog.npy",
        dune_file="DuneStart_1000dam.npy",
    )

    cascade_10kyr_pt75_04 = RUN_4_CASCADE_noAST_Rave_SLR_pt004_NoHumans(
        nt=10000,
        rmin=0.55,  # rave = 0.75
        rmax=0.95,
        name="4-B3D_noAST_Rave_pt75_SLR_pt004_10k-yrs_04",
        storm_file="StormSeries_10kyrs_VCR_Berm1pt9m_Slope0pt04_04.npy",
        elevation_file="InitElevHog.npy",
        dune_file="DuneStart_1000dam.npy",
    )

    cascade_10kyr_pt75_05 = RUN_4_CASCADE_noAST_Rave_SLR_pt004_NoHumans(
        nt=10000,
        rmin=0.55,  # rave = 0.75
        rmax=0.95,
        name="4-B3D_noAST_Rave_pt75_SLR_pt004_10k-yrs_05",
        storm_file="StormSeries_10kyrs_VCR_Berm1pt9m_Slope0pt04_05.npy",
        elevation_file="InitElevHog.npy",
        dune_file="DuneStart_1000dam.npy",
    )

    def sensitivity_tests_Ian_model():
        cascade_10kyr_pt75_Cbbr0pt5 = RUN_4_CASCADE_noAST_Rave_SLR_pt004_NoHumans(
            nt=10000,
            rmin=0.55,  # rave = 0.75
            rmax=0.95,
            name="4-CASCADE_noAST_Rave_pt75_SLR_pt004_10k-yrs_Cbb0pt5",
            storm_file="StormSeries_10kyrs_VCR_Berm1pt9m_Slope0pt04.npy",
            elevation_file="InitElevHog.npy",
            dune_file="DuneStart_1000dam.npy",
        )

        cascade_10kyr_pt75_old_storms = RUN_4_CASCADE_noAST_Rave_SLR_pt004_NoHumans(
            nt=10000,
            rmin=0.55,  # rave = 0.75
            rmax=0.95,
            name="4-CASCADE_noAST_Rave_pt75_SLR_pt004_10k-yrs_OLD_STORMS",
            storm_file="Default_StormTimeSeries_10k-yr.npy",
            elevation_file="InitElevHog.npy",
            dune_file="DuneStart_1000dam.npy",
        )

        cascade_10kyr_pt75_old_storms_Cbbr0pt5 = (
            RUN_4_CASCADE_noAST_Rave_SLR_pt004_NoHumans(
                nt=10000,
                rmin=0.55,  # rave = 0.75
                rmax=0.95,
                name="4-CASCADE_noAST_Rave_pt75_SLR_pt004_10k-yrs_OLD_STORMS_Cbb0pt5",
                storm_file="Default_StormTimeSeries_10k-yr.npy",
                elevation_file="InitElevHog.npy",
                dune_file="DuneStart_1000dam.npy",
            )
        )

        # manually changed the berm elevation to 2.0 in the yaml
        cascade_10kyr_pt75_old_storms_BermEl2 = (
            RUN_4_CASCADE_noAST_Rave_SLR_pt004_NoHumans(
                nt=10000,
                rmin=0.55,  # rave = 0.75
                rmax=0.95,
                name="4-CASCADE_noAST_Rave_pt75_SLR_pt004_10k-yrs_OLD_STORMS_BermEl2",
                storm_file="Default_StormTimeSeries_10k-yr.npy",
                elevation_file="InitElevHog.npy",
                dune_file="DuneStart_1000dam.npy",
            )
        )


def cascade_10kyr_plots():

    datadir = "/Users/KatherineAnardeWheels/Research/BARis/UNC/CNH/CASCADE_save_dir/Run_Output/"
    tmax_pt45 = [10000, 10000, 10000, 10000, 10000]
    name_prefix_45 = "4-B3D_noAST_Rave_pt45_SLR_pt004_10k-yrs_0"
    tmax_pt75 = [5710, 10000, 2878, 10000, 10000]
    # tmax_pt75_old = [5725, 992, 4870, 10000, 6669]
    name_prefix_75 = "4-B3D_noAST_Rave_pt75_SLR_pt004_10k-yrs_0"

    CASCADEplt.supp_10kyr_timeseries(
        datadir, tmax_pt45, name_prefix_45, vertical_line_1=8757, vertical_line_2=802
    )
    CASCADEplt.supp_10kyr_timeseries(
        datadir, tmax_pt75, name_prefix_75, vertical_line_1=4261, vertical_line_2=829
    )


# record of B3D time series initial conditions (storms, dune growth rates, growth parameters) -------------------
def time_series():

    datadir = "/Users/KatherineAnardeWheels/PycharmProjects/CASCADE/B3D_Inputs/"

    StormSeries_NormDist_10kyrs_01 = yearly_storms(
        datadir=datadir,
        storm_list_name="StormList_20k_VCR_Berm1pt9m_Slope0pt04.csv",  # can by .py or .csv
        mean_yearly_storms=8.3,
        SD_yearly_storms=5.9,
        MHW=0.46,  # m NAVD88
        StormStart=2,
        BermEl=1.9,  # m NAVD88, just used for plotting
        model_years=10000,
        bPlot=True,
        bSave=True,
        output_filename="StormSeries_10kyrs_VCR_Berm1pt9m_Slope0pt04_01",
    )

    StormSeries_NormDist_10kyrs_02 = yearly_storms(
        datadir=datadir,
        storm_list_name="StormList_20k_VCR_Berm1pt9m_Slope0pt04.csv",  # can by .py or .csv
        mean_yearly_storms=8.3,
        SD_yearly_storms=5.9,
        MHW=0.46,  # m NAVD88
        StormStart=2,
        BermEl=1.9,  # m NAVD88, just used for plotting
        model_years=10000,
        bPlot=True,
        bSave=True,
        output_filename="StormSeries_10kyrs_VCR_Berm1pt9m_Slope0pt04_02",
    )

    StormSeries_NormDist_10kyrs_03 = yearly_storms(
        datadir=datadir,
        storm_list_name="StormList_20k_VCR_Berm1pt9m_Slope0pt04.csv",  # can by .py or .csv
        mean_yearly_storms=8.3,
        SD_yearly_storms=5.9,
        MHW=0.46,  # m NAVD88
        StormStart=2,
        BermEl=1.9,  # m NAVD88, just used for plotting
        model_years=10000,
        bPlot=True,
        bSave=True,
        output_filename="StormSeries_10kyrs_VCR_Berm1pt9m_Slope0pt04_03",
    )

    StormSeries_NormDist_10kyrs_04 = yearly_storms(
        datadir=datadir,
        storm_list_name="StormList_20k_VCR_Berm1pt9m_Slope0pt04.csv",  # can by .py or .csv
        mean_yearly_storms=8.3,
        SD_yearly_storms=5.9,
        MHW=0.46,  # m NAVD88
        StormStart=2,
        BermEl=1.9,  # m NAVD88, just used for plotting
        model_years=10000,
        bPlot=True,
        bSave=True,
        output_filename="StormSeries_10kyrs_VCR_Berm1pt9m_Slope0pt04_04",
    )

    StormSeries_NormDist_10kyrs_05 = yearly_storms(
        datadir=datadir,
        storm_list_name="StormList_20k_VCR_Berm1pt9m_Slope0pt04.csv",  # can by .py or .csv
        mean_yearly_storms=8.3,
        SD_yearly_storms=5.9,
        MHW=0.46,  # m NAVD88
        StormStart=2,
        BermEl=1.9,  # m NAVD88, just used for plotting
        model_years=10000,
        bPlot=True,
        bSave=True,
        output_filename="StormSeries_10kyrs_VCR_Berm1pt9m_Slope0pt04_05",
    )

    StormSeries_NormDist_1kyrs_01 = yearly_storms(
        datadir=datadir,
        storm_list_name="StormList_20k_VCR_Berm1pt9m_Slope0pt04.csv",  # can by .py or .csv
        mean_yearly_storms=8.3,
        SD_yearly_storms=5.9,
        MHW=0.46,  # m NAVD88
        StormStart=2,
        BermEl=1.9,  # m NAVD88, just used for plotting
        model_years=1000,
        bPlot=True,
        bSave=True,
        output_filename="StormSeries_1kyrs_VCR_Berm1pt9m_Slope0pt04_01",
    )

    StormSeries_NormDist_1kyrs_02 = yearly_storms(
        datadir=datadir,
        storm_list_name="StormList_20k_VCR_Berm1pt9m_Slope0pt04.csv",  # can by .py or .csv
        mean_yearly_storms=8.3,
        SD_yearly_storms=5.9,
        MHW=0.46,  # m NAVD88
        StormStart=2,
        BermEl=1.9,  # m NAVD88, just used for plotting
        model_years=1000,
        bPlot=True,
        bSave=True,
        output_filename="StormSeries_1kyrs_VCR_Berm1pt9m_Slope0pt04_02",
    )

    def one_hundred_increase_storm_intensity_and_frequency():
        number_storms = 105
        datadir = "/Users/KatherineAnardeWheels/PycharmProjects/CASCADE/B3D_Inputs"

        # for iStorm in range(number_storms):
        for iStorm in range(100, 105):

            output_filename = (
                "StormSeries_1kyrs_VCR_Berm1pt9m_Slope0pt04_FutureScenario"
                + str(iStorm)
            )
            shift_storm_intensity(
                datadir=datadir,
                storm_list_name="StormList_20k_VCR_Berm1pt9m_Slope0pt04.csv",  # can by .py or .csv
                mean_yearly_storms=12,
                SD_yearly_storms=5.9,
                shift=0.15,  # shift the TWL distribution to change intensity, m NAVD88; [-0.15, 0.15] for Reeves et al., 2021
                MHW=0.46,  # m NAVD88
                StormStart=2,
                BermEl=1.9,  # m NAVD88, just used for plotting
                model_years=1000,
                bPlot=False,
                bSave=True,
                output_filename=output_filename,
            )

    def one_hundered_ish_1kyr_storms():
        number_storms = 100
        datadir = "/Users/KatherineAnardeWheels/PycharmProjects/CASCADE/B3D_Inputs"

        for iStorm in range(number_storms):

            output_filename = "StormSeries_1kyrs_VCR_Berm1pt9m_Slope0pt04_" + str(
                iStorm
            )
            yearly_storms(
                datadir=datadir,
                storm_list_name="StormList_20k_VCR_Berm1pt9m_Slope0pt04.csv",  # can by .py or .csv
                mean_yearly_storms=8.3,
                SD_yearly_storms=5.9,
                MHW=0.46,  # m NAVD88
                StormStart=2,
                BermEl=1.9,  # m NAVD88, just used for plotting
                model_years=1000,
                bPlot=False,
                bSave=True,
                output_filename=output_filename,
            )

    def BermEl_2m_sensitivity_test():
        name = "StormTimeSeries_10k-yr.npy"
        yearly_storms(
            datadir=datadir,
            name=name,
            storm_list_name="VCRStormList.npy",
            mean_storm=8.3,
            SD_storm=5.9,
            MHW=0.46,
            StormStart=2,
            BermEl=1.9,
            model_years=10000,  # note, this is the number of storms contained in the MSSM model. probably should make more.
        )

        name = "StormTimeSeries_3000yr.npy"
        yearly_storms(
            datadir=datadir,
            name=name,
            storm_list_name="VCRStormList.npy",
            mean_storm=8.3,
            SD_storm=5.9,
            MHW=0.46,
            StormStart=2,
            BermEl=1.9,
            model_years=3000,
        )

        name = "StormTimeSeries_1000yr.npy"
        yearly_storms(
            datadir=datadir,
            name=name,
            storm_list_name="VCRStormList.npy",
            mean_storm=8.3,
            SD_storm=5.9,
            MHW=0.46,
            StormStart=2,
            BermEl=1.9,
            model_years=1000,
        )

        name = "StormTimeSeries_200yr.npy"
        yearly_storms(
            datadir=datadir,
            name=name,
            storm_list_name="VCRStormList.npy",
            mean_storm=8.3,
            SD_storm=5.9,
            MHW=0.46,
            StormStart=2,
            BermEl=1.9,
            model_years=200,
        )

    name = "DuneStart_1000dam.npy"
    gen_dune_height_start(datadir, name, Dstart=0.5, ny=1000)

    name = "growthparam_1000dam.npy"
    gen_alongshore_variable_rmin_rmax(datadir, name, rmin=0.35, rmax=0.85, ny=1000)


# record of 1,000 yr runs and plots for paper -------------------------------------------------------------------
def cascade_1kyr_runs():
    def SLR_sensitivity():
        """
        Completed runs using new elevations on 6/30/222 and 7/24-25/22
        """
        cascade_pt75_low_SLR0pt008 = RUN_7_CASCADE_noAST_Rave_variableSLR_NoHumans(
            nt=200,
            rmin=0.55,
            rmax=0.95,  # rave = 0.75
            name="7-B3D_Rave_pt75_Natural_low_0pt008SLR",
            storm_file="StormSeries_1kyrs_VCR_Berm1pt9m_Slope0pt04_01.npy",
            elevation_file="b3d_pt75_4261yrs_low-elevations.csv",
            dune_file="barrier3d-default-dunes.npy",
            sea_level_rise_rate=0.008,  # m/yr
            sea_level_constant=True,
        )

        cascade_pt75_low_SLR0pt012 = RUN_7_CASCADE_noAST_Rave_variableSLR_NoHumans(
            nt=200,
            rmin=0.55,
            rmax=0.95,  # rave = 0.75
            name="7-B3D_Rave_pt75_Natural_low_0pt012SLR",
            storm_file="StormSeries_1kyrs_VCR_Berm1pt9m_Slope0pt04_01.npy",
            elevation_file="b3d_pt75_4261yrs_low-elevations.csv",
            dune_file="barrier3d-default-dunes.npy",
            sea_level_rise_rate=0.012,  # m/yr
            sea_level_constant=True,
        )

        cascade_pt75_high_SLR0pt008 = RUN_7_CASCADE_noAST_Rave_variableSLR_NoHumans(
            nt=200,
            rmin=0.55,
            rmax=0.95,  # rave = 0.75
            name="7-B3D_Rave_pt75_Natural_high_0pt008SLR",
            storm_file="StormSeries_1kyrs_VCR_Berm1pt9m_Slope0pt04_01.npy",
            elevation_file="b3d_pt75_829yrs_high-elevations.csv",
            dune_file="barrier3d-default-dunes.npy",
            sea_level_rise_rate=0.008,  # m/yr
            sea_level_constant=True,
        )

        cascade_pt75_high_SLR0pt012 = RUN_7_CASCADE_noAST_Rave_variableSLR_NoHumans(
            nt=200,
            rmin=0.55,
            rmax=0.95,  # rave = 0.75
            name="7-B3D_Rave_pt75_Natural_high_0pt012SLR",
            storm_file="StormSeries_1kyrs_VCR_Berm1pt9m_Slope0pt04_01.npy",
            elevation_file="b3d_pt75_829yrs_high-elevations.csv",
            dune_file="barrier3d-default-dunes.npy",
            sea_level_rise_rate=0.012,  # m/yr
            sea_level_constant=True,
        )

        cascade_pt45_low_SLR0pt008 = RUN_7_CASCADE_noAST_Rave_variableSLR_NoHumans(
            nt=200,
            rmin=0.25,
            rmax=0.65,  # rave = 0.45
            name="7-B3D_Rave_pt45_Natural_low_0pt008SLR",
            storm_file="StormSeries_1kyrs_VCR_Berm1pt9m_Slope0pt04_01.npy",
            elevation_file="b3d_pt45_8757yrs_low-elevations.csv",
            dune_file="barrier3d-default-dunes.npy",
            sea_level_rise_rate=0.008,  # m/yr
            sea_level_constant=True,
        )

        cascade_pt45_low_SLR0pt012 = RUN_7_CASCADE_noAST_Rave_variableSLR_NoHumans(
            nt=200,
            rmin=0.25,
            rmax=0.65,  # rave = 0.45
            name="7-B3D_Rave_pt45_Natural_low_0pt012SLR",
            storm_file="StormSeries_1kyrs_VCR_Berm1pt9m_Slope0pt04_01.npy",
            elevation_file="b3d_pt45_8757yrs_low-elevations.csv",
            dune_file="barrier3d-default-dunes.npy",
            sea_level_rise_rate=0.012,  # m/yr
            sea_level_constant=True,
        )

        cascade_pt45_high_SLR0pt008 = RUN_7_CASCADE_noAST_Rave_variableSLR_NoHumans(
            nt=200,
            rmin=0.25,
            rmax=0.65,  # rave = 0.45
            name="7-B3D_Rave_pt45_Natural_high_0pt008SLR",
            storm_file="StormSeries_1kyrs_VCR_Berm1pt9m_Slope0pt04_01.npy",
            elevation_file="b3d_pt45_802yrs_high-elevations.csv",
            dune_file="barrier3d-default-dunes.npy",
            sea_level_rise_rate=0.008,  # m/yr
            sea_level_constant=True,
        )

        cascade_pt45_high_SLR0pt012 = RUN_7_CASCADE_noAST_Rave_variableSLR_NoHumans(
            nt=200,
            rmin=0.25,
            rmax=0.65,  # rave = 0.45
            name="7-B3D_Rave_pt45_Natural_high_0pt012SLR",
            storm_file="StormSeries_1kyrs_VCR_Berm1pt9m_Slope0pt04_01.npy",
            elevation_file="b3d_pt45_802yrs_high-elevations.csv",
            dune_file="barrier3d-default-dunes.npy",
            sea_level_rise_rate=0.012,  # m/yr
            sea_level_constant=True,
        )

        # start of accelerated SLR scenarios
        # for the accelerated SLR scenario, I had to hard code the parameters that correspond to the
        # Rohling et al. (2013) 68% upper bound for AD2000-2200. SLRR starts at 0.003 m/yr and ends at 0.022 m/yr;
        # matches with the bounds of RCP8.5 SLR by 2100 and 2200
        cascade_pt75_low_SLRacc = RUN_7_CASCADE_noAST_Rave_variableSLR_NoHumans(
            nt=200,
            rmin=0.55,
            rmax=0.95,  # rave = 0.75
            name="7-B3D_Rave_pt75_Natural_low_AccSLR",
            storm_file="StormSeries_1kyrs_VCR_Berm1pt9m_Slope0pt04_01.npy",
            elevation_file="b3d_pt75_4261yrs_low-elevations.csv",
            dune_file="barrier3d-default-dunes.npy",
            sea_level_rise_rate=0.008,  # dummy
            sea_level_constant=False,  # accelerated
        )

        cascade_pt75_high_SLRacc = RUN_7_CASCADE_noAST_Rave_variableSLR_NoHumans(
            nt=200,
            rmin=0.55,
            rmax=0.95,  # rave = 0.75
            name="7-B3D_Rave_pt75_Natural_high_AccSLR",
            storm_file="StormSeries_1kyrs_VCR_Berm1pt9m_Slope0pt04_01.npy",
            elevation_file="b3d_pt75_829yrs_high-elevations.csv",
            dune_file="barrier3d-default-dunes.npy",
            sea_level_rise_rate=0.008,  # m/yr
            sea_level_constant=False,
        )

        cascade_pt45_low_SLRacc = RUN_7_CASCADE_noAST_Rave_variableSLR_NoHumans(
            nt=200,
            rmin=0.25,
            rmax=0.65,  # rave = 0.45
            name="7-B3D_Rave_pt45_Natural_low_AccSLR",
            storm_file="StormSeries_1kyrs_VCR_Berm1pt9m_Slope0pt04_01.npy",
            elevation_file="b3d_pt45_8757yrs_low-elevations.csv",
            dune_file="barrier3d-default-dunes.npy",
            sea_level_rise_rate=0.012,  # m/yr
            sea_level_constant=False,
        )

        cascade_pt45_high_SLRacc = RUN_7_CASCADE_noAST_Rave_variableSLR_NoHumans(
            nt=200,
            rmin=0.25,
            rmax=0.65,  # rave = 0.45
            name="7-B3D_Rave_pt45_Natural_high_AccSLR",
            storm_file="StormSeries_1kyrs_VCR_Berm1pt9m_Slope0pt04_01.npy",
            elevation_file="b3d_pt45_802yrs_high-elevations.csv",
            dune_file="barrier3d-default-dunes.npy",
            sea_level_rise_rate=0.008,  # m/yr
            sea_level_constant=False,
        )

    def natural():
        """
        Completed runs using new elevations on 6/30/222
        """
        cascade_pt75_low = RUN_4_CASCADE_noAST_Rave_SLR_pt004_NoHumans(
            nt=1000,
            rmin=0.55,
            rmax=0.95,  # rave = 0.75
            name="4-B3D_Rave_pt75_Natural_low",
            storm_file="StormSeries_1kyrs_VCR_Berm1pt9m_Slope0pt04_01.npy",
            elevation_file="b3d_pt75_4261yrs_low-elevations.csv",
            dune_file="barrier3d-default-dunes.npy",
        )

        cascade_pt75_high = RUN_4_CASCADE_noAST_Rave_SLR_pt004_NoHumans(
            nt=1000,
            rmin=0.55,
            rmax=0.95,  # rave = 0.75
            name="4-B3D_Rave_pt75_Natural_high",
            storm_file="StormSeries_1kyrs_VCR_Berm1pt9m_Slope0pt04_01.npy",
            elevation_file="b3d_pt75_829yrs_high-elevations.csv",
            dune_file="barrier3d-default-dunes.npy",
        )

        cascade_pt45_low = RUN_4_CASCADE_noAST_Rave_SLR_pt004_NoHumans(
            nt=1000,
            rmin=0.25,
            rmax=0.65,  # rave = 0.45
            name="4-B3D_Rave_pt45_Natural_low",
            storm_file="StormSeries_1kyrs_VCR_Berm1pt9m_Slope0pt04_01.npy",
            elevation_file="b3d_pt45_8757yrs_low-elevations.csv",
            dune_file="barrier3d-default-dunes.npy",
        )

        cascade_pt45_high = RUN_4_CASCADE_noAST_Rave_SLR_pt004_NoHumans(
            nt=1000,
            rmin=0.25,
            rmax=0.65,  # rave = 0.45
            name="4-B3D_Rave_pt45_Natural_high",
            storm_file="StormSeries_1kyrs_VCR_Berm1pt9m_Slope0pt04_01.npy",
            elevation_file="b3d_pt45_802yrs_high-elevations.csv",
            dune_file="barrier3d-default-dunes.npy",
        )

        def averages():
            def one_hundred_natural_runs(
                name_prefix, rmin, rmax, elevation_file, year_start, year_end
            ):

                for iStorm in range(year_start, year_end):
                    name = name_prefix + str(iStorm)
                    storm_file = (
                        "StormSeries_1kyrs_VCR_Berm1pt9m_Slope0pt04_"
                        + str(iStorm)
                        + ".npy"
                    )

                    RUN_4_CASCADE_noAST_Rave_SLR_pt004_NoHumans(
                        nt=1000,
                        rmin=rmin,
                        rmax=rmax,
                        name=name,
                        storm_file=storm_file,
                        elevation_file=elevation_file,
                        dune_file="barrier3d-default-dunes.npy",
                    )

            one_hundred_natural_runs(
                name_prefix="4-B3D_Rave_pt45_Natural_low",
                rmin=0.25,
                rmax=0.65,
                elevation_file="b3d_pt45_8757yrs_low-elevations.csv",
                year_start=22,
                year_end=100,
            )

            one_hundred_natural_runs(
                name_prefix="4-B3D_Rave_pt45_Natural_high",
                rmin=0.25,
                rmax=0.65,
                elevation_file="b3d_pt45_802yrs_high-elevations.csv",
                year_start=22,
                year_end=100,
            )

            one_hundred_natural_runs(
                name_prefix="4-B3D_Rave_pt75_Natural_low",
                rmin=0.55,
                rmax=0.95,  # rave = 0.75
                elevation_file="b3d_pt75_4261yrs_low-elevations.csv",
                year_start=0,
                year_end=50,
            )

            one_hundred_natural_runs(
                name_prefix="4-B3D_Rave_pt75_Natural_high",
                rmin=0.55,
                rmax=0.95,  # rave = 0.75
                elevation_file="b3d_pt75_829yrs_high-elevations.csv",
                year_start=92,
                year_end=100,
            )

    def roadways():
        def pt75():
            def low():
                # Barrier has HEIGHT DROWNED at t = 136 years
                cascade_pt75_h1m_low = RUN_6_CASCADE_noAST_Rave_SLR_pt004_Roadways(
                    nt=1000,
                    rmin=0.55,
                    rmax=0.95,  # rave = 0.75
                    name="6-B3D_Rave_pt75_Roadways_1mDune_20mSetback_20mWidth_low",
                    road_ele=0.6,  # average initial elevation, 0.575 m MHW
                    road_width=20,  # m
                    road_setback=20,  # m
                    dune_design_elevation=1.6,  # m MHW, rebuild to 1 m dune above the roadway
                    dune_minimum_elevation=1.1,  # m MHW, allow dune to erode down to 0.5 m above the roadway
                    storm_file="StormSeries_1kyrs_VCR_Berm1pt9m_Slope0pt04_01.npy",
                    elevation_file="b3d_pt75_4261yrs_low-elevations.csv",
                    dune_file="barrier3d-default-dunes.npy",
                    background_erosion=0.0,
                )

                # Barrier has HEIGHT DROWNED at t = 136 years
                cascade_pt75_h2m_low = RUN_6_CASCADE_noAST_Rave_SLR_pt004_Roadways(
                    nt=1000,
                    rmin=0.55,
                    rmax=0.95,  # rave = 0.75
                    name="6-B3D_Rave_pt75_Roadways_2mDune_20mSetback_20mWidth_low",
                    road_ele=0.6,
                    road_width=20,  # m
                    road_setback=20,  # m
                    dune_design_elevation=2.6,  # m MHW, rebuild to 2 m dune above the roadway
                    dune_minimum_elevation=1.1,  # m MHW, allow dune to erode down to 0.5 m above the roadway
                    storm_file="StormSeries_1kyrs_VCR_Berm1pt9m_Slope0pt04_01.npy",
                    elevation_file="b3d_pt75_4261yrs_low-elevations.csv",
                    dune_file="barrier3d-default-dunes.npy",
                    background_erosion=0.0,
                )

                # Island is too narrow for roadway to be relocated. Roadway eaten up by dunes at 73 years
                cascade_pt75_h2m_low_BE1m = RUN_6_CASCADE_noAST_Rave_SLR_pt004_Roadways(
                    nt=1000,
                    rmin=0.55,
                    rmax=0.95,  # rave = 0.75
                    name="6-B3D_Rave_pt75_Roadways_2mDune_20mSetback_20mWidth_low_BE1m",
                    road_ele=0.6,
                    road_width=20,  # m
                    road_setback=20,  # m
                    dune_design_elevation=2.6,  # m MHW, rebuild to 2 m dune above the roadway
                    dune_minimum_elevation=1.1,  # m MHW, allow dune to erode down to 0.5 m above the roadway
                    storm_file="StormSeries_1kyrs_VCR_Berm1pt9m_Slope0pt04_01.npy",
                    elevation_file="b3d_pt75_4261yrs_low-elevations.csv",
                    dune_file="barrier3d-default-dunes.npy",
                    background_erosion=-1,
                )

                # Barrier has HEIGHT DROWNED at t = 132 years
                cascade_pt75_h3m_low = RUN_6_CASCADE_noAST_Rave_SLR_pt004_Roadways(
                    nt=1000,
                    rmin=0.55,
                    rmax=0.95,  # rave = 0.75
                    name="6-B3D_Rave_pt75_Roadways_3mDune_20mSetback_20mWidth_low",
                    road_ele=0.6,
                    road_width=20,  # m
                    road_setback=20,  # m
                    dune_design_elevation=3.6,  # m MHW, rebuild to 3 m dune above the roadway
                    dune_minimum_elevation=1.1,  # m MHW, allow dune to erode down to 0.5 m above the roadway
                    storm_file="StormSeries_1kyrs_VCR_Berm1pt9m_Slope0pt04_01.npy",
                    elevation_file="b3d_pt75_4261yrs_low-elevations.csv",
                    dune_file="barrier3d-default-dunes.npy",
                    background_erosion=0.0,
                )

            def high():
                # Roadway width drowned at 535 years, 20.0% of road borders water
                cascade_pt75_h1m_high = RUN_6_CASCADE_noAST_Rave_SLR_pt004_Roadways(
                    nt=1000,
                    rmin=0.55,
                    rmax=0.95,  # rave = 0.75
                    name="6-B3D_Rave_pt75_Roadways_1mDune_20mSetback_20mWidth_high",
                    road_ele=2.1,  # average initial elevation, 2.14 m MHW
                    road_width=20,  # m
                    road_setback=20,  # m
                    dune_design_elevation=3.1,  # m MHW, rebuild to 1 m dune above the roadway
                    dune_minimum_elevation=2.6,  # m MHW, allow dune to erode down to 0.5 m above the roadway
                    storm_file="StormSeries_1kyrs_VCR_Berm1pt9m_Slope0pt04_01.npy",
                    elevation_file="b3d_pt75_829yrs_high-elevations.csv",
                    dune_file="barrier3d-default-dunes.npy",
                    background_erosion=0.0,
                )

                # Roadway width drowned at 520 years, 20.0% of road borders water
                # Barrier has HEIGHT DROWNED at t = 571 years
                cascade_pt75_h2m_high = RUN_6_CASCADE_noAST_Rave_SLR_pt004_Roadways(
                    nt=1000,
                    rmin=0.55,
                    rmax=0.95,  # rave = 0.75
                    name="6-B3D_Rave_pt75_Roadways_2mDune_20mSetback_20mWidth_high",
                    road_ele=2.1,
                    road_width=20,  # m
                    road_setback=20,  # m
                    dune_design_elevation=4.1,  # m MHW, rebuild to 2 m dune above the roadway
                    dune_minimum_elevation=2.6,  # m MHW, allow dune to erode down to 0.5 m above the roadway
                    storm_file="StormSeries_1kyrs_VCR_Berm1pt9m_Slope0pt04_01.npy",
                    elevation_file="b3d_pt75_829yrs_high-elevations.csv",
                    dune_file="barrier3d-default-dunes.npy",
                    background_erosion=0.0,
                )

                # Island is to narrow for roadway to be relocated. Roadway eaten up by dunes at 395 years
                cascade_pt75_h3m_high = RUN_6_CASCADE_noAST_Rave_SLR_pt004_Roadways(
                    nt=1000,
                    rmin=0.55,
                    rmax=0.95,  # rave = 0.75
                    name="6-B3D_Rave_pt75_Roadways_3mDune_20mSetback_20mWidth_high",
                    road_ele=2.1,
                    road_width=20,  # m
                    road_setback=20,  # m
                    dune_design_elevation=5.1,  # m MHW, rebuild to 3 m dune above the roadway
                    dune_minimum_elevation=2.6,  # m MHW, allow dune to erode down to 0.5 m above the roadway
                    storm_file="StormSeries_1kyrs_VCR_Berm1pt9m_Slope0pt04_01.npy",
                    elevation_file="b3d_pt75_829yrs_high-elevations.csv",
                    dune_file="barrier3d-default-dunes.npy",
                    background_erosion=0.0,
                )

        def pt45():
            def low():
                def averages():
                    def one_hundred_roadway_runs(
                        name_prefix, year_start, year_end, dune_design_elevation
                    ):
                        for iStorm in range(year_start, year_end):
                            name = name_prefix + str(iStorm)
                            storm_file = (
                                "StormSeries_1kyrs_VCR_Berm1pt9m_Slope0pt04_"
                                + str(iStorm)
                                + ".npy"
                            )

                            RUN_6_CASCADE_noAST_Rave_SLR_pt004_Roadways(
                                nt=1000,
                                rmin=0.25,
                                rmax=0.65,  # rave = 0.45
                                name=name,
                                road_ele=1.6,  # average initial elevation 1.64 m MHW
                                road_width=20,  # m
                                road_setback=20,  # m
                                dune_design_elevation=dune_design_elevation,  # m MHW, rebuild to 1 m dune above the roadway
                                dune_minimum_elevation=2.1,
                                # m MHW, allow dune to erode down to 0.5 m above the roadway
                                storm_file=storm_file,
                                elevation_file="b3d_pt45_8757yrs_low-elevations.csv",
                                dune_file="barrier3d-default-dunes.npy",
                                background_erosion=0.0,
                            )

                    one_hundred_roadway_runs(
                        name_prefix="6-B3D_Rave_pt45_Roadways_1mDune_20mSetback_20mWidth_low",
                        year_start=0,
                        year_end=100,
                        dune_design_elevation=2.6,
                    )

                    one_hundred_roadway_runs(
                        name_prefix="6-B3D_Rave_pt45_Roadways_2mDune_20mSetback_20mWidth_low",
                        year_start=0,
                        year_end=100,
                        dune_design_elevation=3.6,
                    )

                    one_hundred_roadway_runs(
                        name_prefix="6-B3D_Rave_pt45_Roadways_3mDune_20mSetback_20mWidth_low",
                        year_start=0,
                        year_end=100,
                        dune_design_elevation=4.6,
                    )

                # Roadway width drowned at 544 years, 20.0% of road borders water
                cascade_pt45_h1m_low = RUN_6_CASCADE_noAST_Rave_SLR_pt004_Roadways(
                    nt=1000,
                    rmin=0.25,
                    rmax=0.65,  # rave = 0.45
                    name="6-B3D_Rave_pt45_Roadways_1mDune_20mSetback_20mWidth_low",
                    road_ele=1.6,  # average initial elevation 1.64 m MHW
                    road_width=20,  # m
                    road_setback=20,  # m
                    dune_design_elevation=2.6,  # m MHW, rebuild to 1 m dune above the roadway
                    dune_minimum_elevation=2.1,  # m MHW, allow dune to erode down to 0.5 m above the roadway
                    storm_file="StormSeries_1kyrs_VCR_Berm1pt9m_Slope0pt04_01.npy",
                    elevation_file="b3d_pt45_8757yrs_low-elevations.csv",
                    dune_file="barrier3d-default-dunes.npy",
                    background_erosion=0.0,
                )

                # Roadway width drowned at 533 years, 20.0% of road borders water
                cascade_pt45_h2m_low = RUN_6_CASCADE_noAST_Rave_SLR_pt004_Roadways(
                    nt=1000,
                    rmin=0.25,
                    rmax=0.65,  # rave = 0.45
                    name="6-B3D_Rave_pt45_Roadways_2mDune_20mSetback_20mWidth_low",
                    road_ele=1.6,
                    road_width=20,  # m
                    road_setback=20,  # m
                    dune_design_elevation=3.6,  # m MHW, rebuild to 2 m dune above the roadway
                    dune_minimum_elevation=2.1,  # m MHW, allow dune to erode down to 0.5 m above the roadway
                    storm_file="StormSeries_1kyrs_VCR_Berm1pt9m_Slope0pt04_01.npy",
                    elevation_file="b3d_pt45_8757yrs_low-elevations.csv",
                    dune_file="barrier3d-default-dunes.npy",
                    background_erosion=0.0,
                )

                # # Island is too narrow for roadway to be relocated. Roadway eaten up by dunes at 187 years
                # cascade_pt45_h2m_low_BE1m = RUN_6_CASCADE_noAST_Rave_SLR_pt004_Roadways(
                #     nt=1000,
                #     rmin=0.25,
                #     rmax=0.65,  # rave = 0.45
                #     name="6-B3D_Rave_pt45_Roadways_2mDune_20mSetback_20mWidth_low_BE1m",
                #     road_ele=1.6,
                #     road_width=20,  # m
                #     road_setback=20,  # m
                #     dune_design_elevation=3.6,  # m MHW, rebuild to 2 m dune above the roadway
                #     dune_minimum_elevation=2.1,  # m MHW, allow dune to erode down to 0.5 m above the roadway
                #     storm_file="StormSeries_1kyrs_VCR_Berm1pt9m_Slope0pt04_01.npy",
                #     elevation_file="b3d_pt45_8757yrs_low-elevations.csv",
                #     dune_file="barrier3d-default-dunes.npy",
                #     background_erosion=-1,
                # )

                # Roadway width drowned at 322 years, 20.0% of road borders water
                cascade_pt45_h3m_low = RUN_6_CASCADE_noAST_Rave_SLR_pt004_Roadways(
                    nt=1000,
                    rmin=0.25,
                    rmax=0.65,  # rave = 0.45
                    name="6-B3D_Rave_pt45_Roadways_3mDune_20mSetback_20mWidth_low",
                    road_ele=1.6,
                    road_width=20,  # m
                    road_setback=20,  # m
                    dune_design_elevation=4.6,  # m MHW, rebuild to 3 m dune above the roadway
                    dune_minimum_elevation=2.1,  # m MHW, allow dune to erode down to 0.5 m above the roadway
                    storm_file="StormSeries_1kyrs_VCR_Berm1pt9m_Slope0pt04_01.npy",
                    elevation_file="b3d_pt45_8757yrs_low-elevations.csv",
                    dune_file="barrier3d-default-dunes.npy",
                    background_erosion=0.0,
                )

            def high():
                # Island is too narrow for roadway to be relocated. Roadway eaten up by dunes at 650 years
                cascade_pt45_h1m_high = RUN_6_CASCADE_noAST_Rave_SLR_pt004_Roadways(
                    nt=1000,
                    rmin=0.25,
                    rmax=0.65,  # rave = 0.45
                    name="6-B3D_Rave_pt45_Roadways_1mDune_20mSetback_20mWidth_high",
                    road_ele=1.8,  # initial average, 1.83 m MHW
                    road_width=20,  # m
                    road_setback=20,  # m
                    dune_design_elevation=2.8,  # m MHW, rebuild to 1 m dune above the roadway
                    dune_minimum_elevation=2.3,  # m MHW, allow dune to erode down to 0.5 m above the roadway
                    storm_file="StormSeries_1kyrs_VCR_Berm1pt9m_Slope0pt04_01.npy",
                    elevation_file="b3d_pt45_802yrs_high-elevations.csv",
                    dune_file="barrier3d-default-dunes.npy",
                    background_erosion=0.0,
                )

                # Roadway width drowned at 628 years, 20.0% of road borders water
                cascade_pt45_h2m_high = RUN_6_CASCADE_noAST_Rave_SLR_pt004_Roadways(
                    nt=1000,
                    rmin=0.25,
                    rmax=0.65,  # rave = 0.45
                    name="6-B3D_Rave_pt45_Roadways_2mDune_20mSetback_20mWidth_high",
                    road_ele=1.8,
                    road_width=20,  # m
                    road_setback=20,  # m
                    dune_design_elevation=3.8,  # m MHW, rebuild to 2 m dune above the roadway
                    dune_minimum_elevation=2.3,  # m MHW, allow dune to erode down to 0.5 m above the roadway
                    storm_file="StormSeries_1kyrs_VCR_Berm1pt9m_Slope0pt04_01.npy",
                    elevation_file="b3d_pt45_802yrs_high-elevations.csv",
                    dune_file="barrier3d-default-dunes.npy",
                    background_erosion=0.0,
                )

                # Roadway width drowned at 522 years, 20.0% of road borders water
                cascade_pt45_h3m_high = RUN_6_CASCADE_noAST_Rave_SLR_pt004_Roadways(
                    nt=1000,
                    rmin=0.25,
                    rmax=0.65,  # rave = 0.45
                    name="6-B3D_Rave_pt45_Roadways_3mDune_20mSetback_20mWidth_high",
                    road_ele=1.8,
                    road_width=20,  # m
                    road_setback=20,  # m
                    dune_design_elevation=4.8,  # m MHW, rebuild to 3 m dune above the roadway
                    dune_minimum_elevation=2.3,  # m MHW, allow dune to erode down to 0.5 m above the roadway
                    storm_file="StormSeries_1kyrs_VCR_Berm1pt9m_Slope0pt04_01.npy",
                    elevation_file="b3d_pt45_802yrs_high-elevations.csv",
                    dune_file="barrier3d-default-dunes.npy",
                    background_erosion=0.0,
                )

        def roadway_sensitivity_abandonment_criteria():

            # test the sensitivity of varying the number of water cells that border the roadway as a metric to stop
            # managing the road for the most extreme barrier trajectory (high dune growth rate, low barrier)

            # Roadway width drowned at 185 years, 10.0% of road borders water
            cascade_pt45_h2m_low_10percent = RUN_6_CASCADE_noAST_Rave_SLR_pt004_Roadways(
                nt=700,
                rmin=0.25,
                rmax=0.65,  # rave = 0.45
                name="6-B3D_Rave_pt45_Roadways_2mDune_20mSetback_20mWidth_low_10percent",
                road_ele=1.6,  # average initial elevation, 0.575 m MHW
                road_width=20,  # m
                road_setback=20,  # m
                dune_design_elevation=3.6,  # m MHW, rebuild to 2 m dune above the roadway
                dune_minimum_elevation=2.1,  # m MHW, allow dune to erode down to 0.5 m above the roadway
                storm_file="StormSeries_1kyrs_VCR_Berm1pt9m_Slope0pt04_01.npy",
                elevation_file="b3d_pt45_8757yrs_low-elevations.csv",
                dune_file="barrier3d-default-dunes.npy",
                background_erosion=0.0,
                percent_water_cells_sensitivity=0.1,
            )

            # Roadway width drowned at 187 years, 20.0% of road borders water
            cascade_pt45_h2m_low_20percent = RUN_6_CASCADE_noAST_Rave_SLR_pt004_Roadways(
                nt=700,
                rmin=0.25,
                rmax=0.65,  # rave = 0.45
                name="6-B3D_Rave_pt45_Roadways_2mDune_20mSetback_20mWidth_low_20percent",
                road_ele=1.6,
                road_width=20,  # m
                road_setback=20,  # m
                dune_design_elevation=3.6,  # m MHW, rebuild to 2 m dune above the roadway
                dune_minimum_elevation=2.1,  # m MHW, allow dune to erode down to 0.5 m above the roadway
                storm_file="StormSeries_1kyrs_VCR_Berm1pt9m_Slope0pt04_01.npy",
                elevation_file="b3d_pt45_8757yrs_low-elevations.csv",
                dune_file="barrier3d-default-dunes.npy",
                background_erosion=0.0,
                percent_water_cells_sensitivity=0.2,
            )

            # Roadway width drowned at 320 years, 30.0% of road borders water
            cascade_pt45_h2m_low_30percent = RUN_6_CASCADE_noAST_Rave_SLR_pt004_Roadways(
                nt=700,
                rmin=0.25,
                rmax=0.65,  # rave = 0.45
                name="6-B3D_Rave_pt45_Roadways_2mDune_20mSetback_20mWidth_low_30percent",
                road_ele=1.6,
                road_width=20,  # m
                road_setback=20,  # m
                dune_design_elevation=3.6,  # m MHW, rebuild to 2 m dune above the roadway
                dune_minimum_elevation=2.1,  # m MHW, allow dune to erode down to 0.5 m above the roadway
                storm_file="StormSeries_1kyrs_VCR_Berm1pt9m_Slope0pt04_01.npy",
                elevation_file="b3d_pt45_8757yrs_low-elevations.csv",
                dune_file="barrier3d-default-dunes.npy",
                background_erosion=0.0,
                percent_water_cells_sensitivity=0.3,
            )

            # Roadway width drowned at 322 years, 40.0% of road borders water
            cascade_pt45_h2m_low_40percent = RUN_6_CASCADE_noAST_Rave_SLR_pt004_Roadways(
                nt=700,
                rmin=0.25,
                rmax=0.65,  # rave = 0.45
                name="6-B3D_Rave_pt45_Roadways_2mDune_20mSetback_20mWidth_low_40percent",
                road_ele=1.6,
                road_width=20,  # m
                road_setback=20,  # m
                dune_design_elevation=3.6,  # m MHW, rebuild to 2 m dune above the roadway
                dune_minimum_elevation=2.1,  # m MHW, allow dune to erode down to 0.5 m above the roadway
                storm_file="StormSeries_1kyrs_VCR_Berm1pt9m_Slope0pt04_01.npy",
                elevation_file="b3d_pt45_8757yrs_low-elevations.csv",
                dune_file="barrier3d-default-dunes.npy",
                background_erosion=0.0,
                percent_water_cells_sensitivity=0.4,
            )

            # Roadway width drowned at 325 years, 50.0% of road borders water
            cascade_pt45_h2m_low_50percent = RUN_6_CASCADE_noAST_Rave_SLR_pt004_Roadways(
                nt=700,
                rmin=0.25,
                rmax=0.65,  # rave = 0.45
                name="6-B3D_Rave_pt45_Roadways_2mDune_20mSetback_20mWidth_low_50percent",
                road_ele=1.6,
                road_width=20,  # m
                road_setback=20,  # m
                dune_design_elevation=3.6,  # m MHW, rebuild to 2 m dune above the roadway
                dune_minimum_elevation=2.1,  # m MHW, allow dune to erode down to 0.5 m above the roadway
                storm_file="StormSeries_1kyrs_VCR_Berm1pt9m_Slope0pt04_01.npy",
                elevation_file="b3d_pt45_8757yrs_low-elevations.csv",
                dune_file="barrier3d-default-dunes.npy",
                background_erosion=0.0,
                percent_water_cells_sensitivity=0.5,
            )

        def old_overwash_model():
            def pt45():
                # used 3000 year storm time series for the _low and _high runs
                cascade_pt45_low = RUN_6_CASCADE_noAST_Rave_SLR_pt004_Roadways(
                    nt=550,
                    rmin=0.25,
                    rmax=0.65,  # rave = 0.45
                    name="6-B3D_Rave_pt45_550yrs_Natural_low",
                    road_ele=None,
                    road_width=None,
                    road_setback=None,
                    dune_design_elevation=None,
                    dune_minimum_elevation=None,
                    run_road_mgmt=False,
                )

                cascade_pt45_high = RUN_6_CASCADE_noAST_Rave_SLR_pt004_Roadways(
                    nt=750,
                    rmin=0.25,
                    rmax=0.65,  # rave = 0.45
                    name="6-B3D_Rave_pt45_750yrs_Natural_high",
                    road_ele=None,
                    road_width=None,
                    road_setback=None,
                    dune_design_elevation=None,
                    dune_minimum_elevation=None,
                    run_road_mgmt=False,
                )

                # note that this roadway is higher than the pt75 "low" scenario, and equal to the "high" scenario (just higher
                # starting topoagraphy); drowned at 388 years
                cascade_pt45_h2m_low = RUN_6_CASCADE_noAST_Rave_SLR_pt004_Roadways(
                    nt=387,
                    rmin=0.25,
                    rmax=0.65,  # rave = 0.45
                    name="6-B3D_Rave_pt45_387yrs_Roadways_2mDune_20mSetback_20mWidth_low",
                    road_ele=1.7,  # average of NC-12 is 1.3 m NAVD88, berm ele is 1.4 m MHW (I don't want to set it lower than that)
                    road_width=20,  # m
                    road_setback=20,  # m
                    dune_design_elevation=3.7,  # m MHW, rebuild to 2 m dune above the roadway
                    dune_minimum_elevation=2.2,  # m MHW, allow dune to erode down to 0.5 m above the roadway
                    run_road_mgmt=True,
                )

                # note that this roadway is higher than the pt75 "high" scenario (just higher starting topoagraphy); drowned at 515 years
                cascade_pt45_h2m_high = RUN_6_CASCADE_noAST_Rave_SLR_pt004_Roadways(
                    nt=515,
                    rmin=0.25,
                    rmax=0.65,  # rave = 0.45
                    name="6-B3D_Rave_pt45_515yrs_Roadways_2mDune_20mSetback_20mWidth_high",
                    road_ele=1.7,
                    road_width=20,  # m
                    road_setback=20,  # m
                    dune_design_elevation=3.7,  # m MHW, rebuild to 2 m dune above the roadway
                    dune_minimum_elevation=2.2,  # m MHW, allow dune to erode down to 0.5 m above the roadway
                    run_road_mgmt=True,
                )

                # drowned at 370 years
                cascade_pt45_h3m_low = RUN_6_CASCADE_noAST_Rave_SLR_pt004_Roadways(
                    nt=369,
                    rmin=0.25,
                    rmax=0.65,  # rave = 0.45
                    name="6-B3D_Rave_pt45_369yrs_Roadways_3mDune_20mSetback_20mWidth_low",
                    road_ele=1.7,
                    road_width=20,  # m
                    road_setback=20,  # m
                    dune_design_elevation=4.7,  # m MHW, rebuild to 3 m dune above the roadway
                    dune_minimum_elevation=2.2,  # m MHW, allow dune to erode down to 0.5 m above the roadway
                    run_road_mgmt=True,
                )

                cascade_pt45_h3m_high = RUN_6_CASCADE_noAST_Rave_SLR_pt004_Roadways(
                    nt=473,
                    rmin=0.25,
                    rmax=0.65,  # rave = 0.45
                    name="6-B3D_Rave_pt45_473yrs_Roadways_3mDune_20mSetback_20mWidth_high",
                    road_ele=1.7,
                    road_width=20,  # m
                    road_setback=20,  # m
                    dune_design_elevation=4.7,  # m MHW, rebuild to 3 m dune above the roadway
                    dune_minimum_elevation=2.2,  # m MHW, allow dune to erode down to 0.5 m above the roadway
                    run_road_mgmt=True,
                )

                cascade_pt45_h1m_low = RUN_6_CASCADE_noAST_Rave_SLR_pt004_Roadways(
                    nt=503,
                    rmin=0.25,
                    rmax=0.65,  # rave = 0.45
                    name="6-B3D_Rave_pt45_503yrs_Roadways_1mDune_20mSetback_20mWidth_low",
                    road_ele=1.7,
                    road_width=20,  # m
                    road_setback=20,  # m
                    dune_design_elevation=2.7,  # m MHW, rebuild to 1 m dune above the roadway
                    dune_minimum_elevation=2.2,  # m MHW, allow dune to erode down to 0.5 m above the roadway
                    run_road_mgmt=True,
                )

                cascade_pt45_h1m_high = RUN_6_CASCADE_noAST_Rave_SLR_pt004_Roadways(
                    nt=700,
                    rmin=0.25,
                    rmax=0.65,  # rave = 0.45
                    name="6-B3D_Rave_pt45_700yrs_Roadways_1mDune_20mSetback_20mWidth_high",
                    road_ele=1.7,
                    road_width=20,  # m
                    road_setback=20,  # m
                    dune_design_elevation=2.7,  # m MHW, rebuild to 1 m dune above the roadway
                    dune_minimum_elevation=2.2,  # m MHW, allow dune to erode down to 0.5 m above the roadway
                    run_road_mgmt=True,
                )

            def pt75():
                # REMEMBER TO SWITCH TOPOGRAHPHY FILES ------------------------
                # used 1000 year storm time series for the _low and _high runs

                # ACCIDENTALLY reran and saved over this one when checking the new roadways class...
                # will want to check with new plots that they are the same
                cascade_pt75_low = RUN_6_CASCADE_noAST_Rave_SLR_pt004_Roadways(
                    nt=200,
                    rmin=0.55,
                    rmax=0.95,  # rave = 0.75
                    name="6-B3D_Rave_pt75_200yrs_Natural_low",
                    road_ele=None,
                    road_width=None,
                    road_setback=None,
                    dune_design_elevation=None,
                    dune_minimum_elevation=None,
                    run_road_mgmt=False,
                    storm_file="Default_StormTimeSeries_1000yr.npy",
                    elevation_file="b3d_pt75_4929yrs_low-elevations.csv",
                    dune_file="barrier3d-default-dunes.npy",
                )

                cascade_pt75_high = RUN_6_CASCADE_noAST_Rave_SLR_pt004_Roadways(
                    nt=450,
                    rmin=0.55,
                    rmax=0.95,  # rave = 0.75
                    # name="6-B3D_Rave_pt75_450yrs_Natural_high",
                    name="test2",
                    road_ele=None,
                    road_width=None,
                    road_setback=None,
                    dune_design_elevation=None,
                    dune_minimum_elevation=None,
                    run_road_mgmt=False,
                    storm_file="Default_StormTimeSeries_1000yr.npy",
                    elevation_file="b3d_pt75_8793yrs_high-elevations.csv",
                    dune_file="barrier3d-default-dunes.npy",
                )

                # lowered roadway and decreased setback to accommodate low-lying barrier, roadway drowned at 98 from back bay
                cascade_pt75_h2m_low = RUN_6_CASCADE_noAST_Rave_SLR_pt004_Roadways(
                    nt=97,
                    rmin=0.55,
                    rmax=0.95,  # rave = 0.75
                    # name="6-B3D_Rave_pt75_97yrs_Roadways_2mDune_20mSetback_20mWidth_low",
                    name="test",
                    road_ele=1.4,
                    road_width=20,  # m
                    road_setback=20,  # m
                    dune_design_elevation=3.4,  # m MHW, rebuild to 2 m dune above the roadway
                    dune_minimum_elevation=1.9,  # m MHW, allow dune to erode down to 0.5 m above the roadway, v1 = 2.7 m
                    run_road_mgmt=True,
                    storm_file="Default_StormTimeSeries_1000yr.npy",
                    elevation_file="b3d_pt75_4929yrs_low-elevations.csv",
                    dune_file="barrier3d-default-dunes.npy",
                )

                # able to set the highway slightly higher (0.3 m), kept setback the same; drowned at 428 years
                cascade_pt75_h2m_high = RUN_6_CASCADE_noAST_Rave_SLR_pt004_Roadways(
                    nt=427,
                    rmin=0.55,
                    rmax=0.95,  # rave = 0.75
                    name="6-B3D_Rave_pt75_425yrs_Roadways_2mDune_20mSetback_20mWidth_high",
                    road_ele=1.7,
                    road_width=20,  # m
                    road_setback=20,  # m
                    dune_design_elevation=3.7,  # m MHW, rebuild to 2 m dune above the roadway
                    dune_minimum_elevation=2.2,  # m MHW, allow dune to erode down to 0.5 m above the roadway, v1 = 2.7 m
                    run_road_mgmt=True,
                )

                # lowered roadway and decreased setback to accommodate low-lying barrier, drowned at 71 years
                cascade_pt75_h3m_low = RUN_6_CASCADE_noAST_Rave_SLR_pt004_Roadways(
                    nt=70,
                    rmin=0.55,
                    rmax=0.95,  # rave = 0.75
                    # name="6-B3D_Rave_pt75_70yrs_Roadways_3mDune_20mSetback_20mWidth_low",
                    name="test",
                    road_ele=1.4,  # average of NC-12 is 1.3 m NAVD88, berm ele is 1.4 m MHW (I don't want to set it lower than that)
                    road_width=20,  # m
                    road_setback=20,  # m
                    dune_design_elevation=4.4,  # m MHW, rebuild to 3 m dune above the roadway
                    dune_minimum_elevation=1.9,  # m MHW, allow dune to erode down to 0.5 m above the roadway
                    run_road_mgmt=True,
                    storm_file="Default_StormTimeSeries_3000yr.npy",
                    elevation_file="b3d_pt75_4929yrs_low-elevations.csv",
                    dune_file="barrier3d-default-dunes.npy",
                )

                # able to set the highway slightly higher (0.3 m), kept setback the same; drowned at 358 years
                cascade_pt75_h3m_high = RUN_6_CASCADE_noAST_Rave_SLR_pt004_Roadways(
                    nt=357,
                    rmin=0.55,
                    rmax=0.95,  # rave = 0.75
                    name="6-B3D_Rave_pt75_357yrs_Roadways_3mDune_20mSetback_20mWidth_high",
                    road_ele=1.7,
                    road_width=20,  # m
                    road_setback=20,  # m
                    dune_design_elevation=4.7,  # m MHW, rebuild to 3 m dune above the roadway
                    dune_minimum_elevation=2.2,  # m MHW, allow dune to erode down to 0.5 m above the roadway, v1 = 2.7 m
                    run_road_mgmt=True,
                )

                # lowered roadway and decreased setback to accommodate low-lying barrier, roadway drowned at 117 years
                cascade_pt75_h1m_low = RUN_6_CASCADE_noAST_Rave_SLR_pt004_Roadways(
                    nt=116,
                    rmin=0.55,
                    rmax=0.95,  # rave = 0.75
                    name="6-B3D_Rave_pt75_116yrs_Roadways_1mDune_20mSetback_20mWidth_low",
                    road_ele=1.4,
                    road_width=20,  # m
                    road_setback=20,  # m
                    dune_design_elevation=2.4,  # m MHW, rebuild to 1 m dune above the roadway
                    dune_minimum_elevation=1.9,  # m MHW, allow dune to erode down to 0.5 m above the roadway
                    run_road_mgmt=True,
                )

                # able to set the highway slightly higher (0.3 m), kept setback the same; drowned at 439 years
                cascade_pt75_h1m_high = RUN_6_CASCADE_noAST_Rave_SLR_pt004_Roadways(
                    nt=438,
                    rmin=0.55,
                    rmax=0.95,  # rave = 0.75
                    name="6-B3D_Rave_pt75_438yrs_Roadways_1mDune_20mSetback_20mWidth_high",
                    road_ele=1.7,
                    road_width=20,  # m
                    road_setback=20,  # m
                    dune_design_elevation=2.7,  # m MHW, rebuild to 1 m dune above the roadway
                    dune_minimum_elevation=2.2,  # m MHW, allow dune to erode down to 0.5 m above the roadway, v1 = 2.7 m
                    run_road_mgmt=True,
                )

            # OLD versions using the final output from 10,000 year run -------------------------------
            def old_10k_initial_elevation():
                b3d_pt45 = RUN_6_CASCADE_noAST_Rave_SLR_pt004_Roadways(
                    nt=200,
                    rmin=0.25,
                    rmax=0.65,  # rave = 0.45
                    name="6-B3D_Rave_pt45_200yrs_Natural_v2",
                    road_ele=None,
                    road_width=None,
                    road_setback=None,
                    dune_design_elevation=None,
                    dune_minimum_elevation=None,
                )

                # v1 - didn't drown roadway, v2 - roadway drowned at 160, v3 - new class
                # b3d_pt45_h2m, dunes_rebuilt, road_overwash_volume = RUN_6_CASCADE_noAST_Rave_SLR_pt004_Roadways(
                cascade_pt45_h2m = RUN_6_CASCADE_noAST_Rave_SLR_pt004_Roadways(
                    nt=159,  # 200
                    rmin=0.25,
                    rmax=0.65,  # rave = 0.45
                    name="6-B3D_Rave_pt45_200yrs_Roadways_2mDune_40mSetback_20mWidth_v3",
                    road_ele=1.7,
                    road_width=20,  # m
                    road_setback=40,  # m
                    dune_design_elevation=3.7,  # m MHW, rebuild to 2 m dune above the roadway
                    dune_minimum_elevation=2.2,  # m MHW, allow dune to erode down to 0.5 m above the roadway, v1 = 2.7 m
                    run_road_mgmt=True,
                )

                b3d_pt45_h3m = RUN_6_CASCADE_noAST_Rave_SLR_pt004_Roadways(
                    nt=200,
                    rmin=0.25,
                    rmax=0.65,  # rave = 0.45
                    name="6-B3D_Rave_pt45_200yrs_Roadways_3mDune_40mSetback_20mWidth_v2",
                    road_ele=1.7,
                    road_width=20,  # m
                    road_setback=40,  # m
                    dune_design_elevation=4.7,  # m MHW, rebuild to 3 m dune above the roadway
                    dune_minimum_elevation=2.2,  # m MHW, allow dune to erode down to 0.5 m above the roadway, v1 = 3.7 m
                )

                # v1 drowned at 91 years, v2 - roadway drowned at 101 years
                # b3d_pt45_h1m = RUN_6_CASCADE_noAST_Rave_SLR_pt004_Roadways(
                cascade_pt45_h1m = RUN_6_CASCADE_noAST_Rave_SLR_pt004_Roadways(
                    nt=100,  # 90
                    rmin=0.25,
                    rmax=0.65,  # rave = 0.45
                    name="6-B3D_Rave_pt45_200yrs_Roadways_1mDune_40mSetback_20mWidth_v2_classtest2",
                    road_ele=1.7,
                    road_width=20,  # m
                    road_setback=40,  # m
                    dune_design_elevation=2.7,  # m MHW, rebuild to 1 m dune above the roadway
                    dune_minimum_elevation=2.2,  # m MHW, allow dune to erode down to 0.5 m above the roadway, v1 = 2.4 m
                    run_road_mgmt=True,
                )

                # increase road width
                b3d_pt45_h4 = RUN_6_CASCADE_noAST_Rave_SLR_pt004_Roadways(
                    nt=200,
                    rmin=0.25,
                    rmax=0.65,  # rave = 0.45
                    name="6-B3D_Rave_pt45_200yrs_Roadways_2mDune_40mSetback_30mWidth_v2",
                    road_ele=1.7,
                    road_width=30,  # m
                    road_setback=40,  # m
                    dune_design_elevation=3.7,  # m MHW, rebuild to 2 m dune above the roadway
                    dune_minimum_elevation=2.2,  # m MHW, allow dune to erode down to 0.5 m above the roadway, v1 = 2.7 m
                )

                # change setback distance
                # v1 drowned at 183 years
                b3d_pt45_h5 = RUN_6_CASCADE_noAST_Rave_SLR_pt004_Roadways(
                    nt=200,  # 182
                    rmin=0.25,
                    rmax=0.65,  # rave = 0.45
                    name="6-B3D_Rave_pt45_200yrs_Roadways_2mDune_30mSetback_20mWidth",
                    road_ele=1.7,
                    road_width=20,  # m
                    road_setback=30,  # m
                    dune_design_elevation=3.7,  # m MHW, 2 m dune above the roadway
                    dune_minimum_elevation=2.7,  # m MHW, 1 m dune above the roadway
                )

                # REMEMBER TO SWITCH TOPOGRAHPHY FILES ------------------------
                b3d_pt75 = RUN_6_CASCADE_noAST_Rave_SLR_pt004_Roadways(
                    nt=200,
                    rmin=0.55,
                    rmax=0.95,  # rave = 0.75
                    name="6-B3D_Rave_pt75_200yrs_Natural",
                    road_ele=None,
                    road_width=None,
                    road_setback=None,
                    dune_design_elevation=None,
                    dune_minimum_elevation=None,
                    run_road_mgmt=False,
                )

                # v2, roadway drowned at 157 years
                b3d_pt75_h2m = RUN_6_CASCADE_noAST_Rave_SLR_pt004_Roadways(
                    nt=156,  # 200
                    rmin=0.55,
                    rmax=0.95,  # rave = 0.75
                    name="6-B3D_Rave_pt75_200yrs_Roadways_2mDune_40mSetback_20mWidth_v2",
                    road_ele=1.7,
                    road_width=20,  # m
                    road_setback=40,  # m
                    dune_design_elevation=3.7,  # m MHW, rebuild to 2 m dune above the roadway
                    dune_minimum_elevation=2.2,  # m MHW, allow dune to erode down to 0.5 m above the roadway, v1 = 2.7 m
                )

                b3d_pt75_h3m = RUN_6_CASCADE_noAST_Rave_SLR_pt004_Roadways(
                    nt=200,
                    rmin=0.55,
                    rmax=0.95,  # rave = 0.75
                    name="6-B3D_Rave_pt75_200yrs_Roadways_3mDune_40mSetback_20mWidth_v2",
                    road_ele=1.7,
                    road_width=20,  # m
                    road_setback=40,  # m
                    dune_design_elevation=4.7,  # m MHW, rebuild to 3 m dune above the roadway
                    dune_minimum_elevation=2.2,  # m MHW, allow dune to erode down to 0.5 m above the roadway, v1 = 3.7 m
                )

                # v1 drowned at 88 years, v2 drowned at 105 years, v3 roadway couldn't be relocated at 408 years
                b3d_pt75_h1m = RUN_6_CASCADE_noAST_Rave_SLR_pt004_Roadways(
                    nt=407,  # 87, 104
                    rmin=0.55,
                    rmax=0.95,  # rave = 0.75
                    name="6-B3D_Rave_pt75_407yrs_Roadways_1mDune_40mSetback_20mWidth",
                    road_ele=1.7,
                    road_width=20,  # m
                    road_setback=40,  # m
                    dune_design_elevation=2.7,  # m MHW, rebuild to 1 m dune above the roadway
                    dune_minimum_elevation=2.2,  # m MHW, allow dune to erode down to 0.5 m above the roadway, v1 = 2.4 m
                    run_road_mgmt=True,
                )

    def nourishments():

        # note, we keep all other variables the same for comparison to the roadways scenarios except we rebuild if the
        # dune is eroded to 1-m above the berm

        def rebuild_threshold_1m():
            def pt75_low():
                # Community reached minimum width, drowned at 160 years
                cascade_pt75_h2m_low_nourishment_residential_RT1m = RUN_8_CASCADE_noAST_Rave_SLR_pt004_Nourishment(
                    nt=1000,
                    rmin=0.55,
                    rmax=0.95,  # rave = 0.75
                    name="8-B3D_Rave_pt75_Nourishment_2mDune_lowEle_residential_rebuild1m",
                    dune_design_elevation=2.6,
                    # m MHW, keep dune design height the same as 2m dune above the initial "roadway"
                    storm_file="StormSeries_1kyrs_VCR_Berm1pt9m_Slope0pt04_01.npy",
                    elevation_file="b3d_pt75_4261yrs_low-elevations.csv",
                    dune_file="barrier3d-default-dunes.npy",
                    overwash_filter=40,  # corresponds with residential
                    overwash_to_dune=9,
                    nourishment_volume=100,  # m^3/m
                    beach_width_threshold=30,  # m
                    background_erosion=0.0,
                    rebuild_dune_threshold=1,  # m above the berm elevation
                )

                # Community reached minimum width, drowned at 83 years
                cascade_pt75_h2m_low_nourishment_commercial_RT1m = RUN_8_CASCADE_noAST_Rave_SLR_pt004_Nourishment(
                    nt=1000,
                    rmin=0.55,
                    rmax=0.95,  # rave = 0.75
                    name="8-B3D_Rave_pt75_Nourishment_2mDune_lowEle_commercial_rebuild1m",
                    dune_design_elevation=2.6,  # m MHW, keep dune design height the same as 2 m dune above the "roadway"
                    storm_file="StormSeries_1kyrs_VCR_Berm1pt9m_Slope0pt04_01.npy",
                    elevation_file="b3d_pt75_4261yrs_low-elevations.csv",
                    dune_file="barrier3d-default-dunes.npy",
                    overwash_filter=90,  # corresponds with commercial
                    overwash_to_dune=9,
                    nourishment_volume=100,  # m^3/m
                    beach_width_threshold=30,  # m
                    background_erosion=0.0,
                    rebuild_dune_threshold=1,  # m above the berm elevation
                )

                # Community reached minimum width, drowned at 83 years
                cascade_pt75_h2m_low_nourishment_commercial_BE1m_RT1m = RUN_8_CASCADE_noAST_Rave_SLR_pt004_Nourishment(
                    nt=1000,
                    rmin=0.55,
                    rmax=0.95,  # rave = 0.75
                    name="8-B3D_Rave_pt75_Nourishment_2mDune_lowEle_commercial_backerosion1m_rebuild1m",
                    dune_design_elevation=2.6,  # m MHW, keep dune design height the same as 2 m dune above the "roadway"
                    storm_file="StormSeries_1kyrs_VCR_Berm1pt9m_Slope0pt04_01.npy",
                    elevation_file="b3d_pt75_4261yrs_low-elevations.csv",
                    dune_file="barrier3d-default-dunes.npy",
                    overwash_filter=90,  # corresponds with commercial
                    overwash_to_dune=9,
                    nourishment_volume=100,  # m^3/m
                    beach_width_threshold=30,  # m
                    background_erosion=-1.0,  # m/yr, background shoreline erosion
                    rebuild_dune_threshold=1,  # m above the berm elevation
                )

            def pt75_high():
                # Community reached minimum width, drowned at 550
                cascade_pt75_h2m_high_nourishment_residential_RT1m = RUN_8_CASCADE_noAST_Rave_SLR_pt004_Nourishment(
                    nt=1000,
                    rmin=0.55,
                    rmax=0.95,  # rave = 0.75
                    name="8-B3D_Rave_pt75_Nourishment_2mDune_highEle_residential_rebuild1m",
                    dune_design_elevation=4.1,
                    # m MHW, keep dune design height the same as 2m dune above the initial "roadway" for comparison
                    storm_file="StormSeries_1kyrs_VCR_Berm1pt9m_Slope0pt04_01.npy",
                    elevation_file="b3d_pt75_829yrs_high-elevations.csv",
                    dune_file="barrier3d-default-dunes.npy",
                    overwash_filter=40,  # corresponds with residential
                    overwash_to_dune=9,
                    nourishment_volume=100,  # m^3/m
                    beach_width_threshold=30,  # m
                    background_erosion=0.0,
                    rebuild_dune_threshold=1,  # m above the berm elevation
                )

                # Community reached minimum width, drowned at 518; Barrier has HEIGHT DROWNED at t = 580 years
                cascade_pt75_h2m_high_nourishment_commercial_RT1m = RUN_8_CASCADE_noAST_Rave_SLR_pt004_Nourishment(
                    nt=1000,
                    rmin=0.55,
                    rmax=0.95,  # rave = 0.75
                    name="8-B3D_Rave_pt75_Nourishment_2mDune_highEle_commercial_rebuild1m",
                    dune_design_elevation=4.1,
                    # m MHW, keep dune design height the same as 2 m dune above the "roadway"
                    storm_file="StormSeries_1kyrs_VCR_Berm1pt9m_Slope0pt04_01.npy",
                    elevation_file="b3d_pt75_829yrs_high-elevations.csv",
                    dune_file="barrier3d-default-dunes.npy",
                    overwash_filter=90,  # corresponds with commercial
                    overwash_to_dune=9,
                    nourishment_volume=100,  # m^3/m
                    beach_width_threshold=30,  # m
                    background_erosion=0.0,
                    rebuild_dune_threshold=1,  # m above the berm elevation
                )

                # Community reached minimum width, drowned at 518 years
                cascade_pt75_h2m_high_nourishment_commercial_BE1m_RT1m = RUN_8_CASCADE_noAST_Rave_SLR_pt004_Nourishment(
                    nt=1000,
                    rmin=0.55,
                    rmax=0.95,  # rave = 0.75
                    name="8-B3D_Rave_pt75_Nourishment_2mDune_highEle_commercial_backerosion1m_rebuild1m",
                    dune_design_elevation=4.1,
                    # m MHW, keep dune design height the same as 2 m dune above the "roadway"
                    storm_file="StormSeries_1kyrs_VCR_Berm1pt9m_Slope0pt04_01.npy",
                    elevation_file="b3d_pt75_829yrs_high-elevations.csv",
                    dune_file="barrier3d-default-dunes.npy",
                    overwash_filter=90,  # corresponds with commercial
                    overwash_to_dune=9,
                    nourishment_volume=100,  # m^3/m
                    beach_width_threshold=30,  # m
                    background_erosion=-1.0,  # m/yr, background shoreline erosion
                    rebuild_dune_threshold=1,  # m above the berm elevation
                )

            def pt45_low():
                def averages():
                    def one_hundred_nourishment_runs(
                        name_prefix,
                        year_start,
                        year_end,
                        overwash_filter,
                        background_erosion,
                    ):
                        for iStorm in range(year_start, year_end):
                            name = name_prefix + str(iStorm)
                            storm_file = (
                                "StormSeries_1kyrs_VCR_Berm1pt9m_Slope0pt04_"
                                + str(iStorm)
                                + ".npy"
                            )

                            RUN_8_CASCADE_noAST_Rave_SLR_pt004_Nourishment(
                                nt=1000,
                                rmin=0.25,
                                rmax=0.65,  # rave = 0.45
                                name=name,
                                dune_design_elevation=3.6,
                                # m MHW, keep dune design height the same as 2 m dune above the "roadway"
                                storm_file=storm_file,
                                elevation_file="b3d_pt45_8757yrs_low-elevations.csv",
                                dune_file="barrier3d-default-dunes.npy",
                                overwash_filter=overwash_filter,
                                overwash_to_dune=9,
                                nourishment_volume=100,  # m^3/m
                                beach_width_threshold=30,  # m
                                background_erosion=background_erosion,
                                rebuild_dune_threshold=1,  # m above the berm elevation
                            )

                    one_hundred_nourishment_runs(
                        name_prefix="8-B3D_Rave_pt45_Nourishment_2mDune_lowEle_residential_rebuild1m",
                        year_start=0,
                        year_end=100,
                        overwash_filter=40,
                        background_erosion=0.0,
                    )

                    one_hundred_nourishment_runs(
                        name_prefix="8-B3D_Rave_pt45_Nourishment_2mDune_lowEle_commercial_rebuild1m",
                        year_start=0,
                        year_end=100,
                        overwash_filter=90,
                        background_erosion=0.0,
                    )

                    one_hundred_nourishment_runs(
                        name_prefix="8-B3D_Rave_pt45_Nourishment_2mDune_lowEle_commercial_backerosion1m_rebuild1m",
                        year_start=0,
                        year_end=100,
                        overwash_filter=90,
                        background_erosion=-1.0,
                    )

                # Community reached minimum width, drowned at 407 years
                cascade_pt45_h2m_low_nourishment_residential_RT1m = RUN_8_CASCADE_noAST_Rave_SLR_pt004_Nourishment(
                    nt=1000,
                    rmin=0.25,
                    rmax=0.65,  # rave = 0.45
                    name="8-B3D_Rave_pt45_Nourishment_2mDune_lowEle_residential_rebuild1m",
                    dune_design_elevation=3.6,  # m MHW, keep dune design height the same as 2 m dune above the "roadway"
                    storm_file="StormSeries_1kyrs_VCR_Berm1pt9m_Slope0pt04_01.npy",
                    elevation_file="b3d_pt45_8757yrs_low-elevations.csv",
                    dune_file="barrier3d-default-dunes.npy",
                    overwash_filter=40,
                    overwash_to_dune=9,
                    nourishment_volume=100,  # m^3/m
                    beach_width_threshold=30,  # m
                    background_erosion=0.0,
                    rebuild_dune_threshold=1,  # m above the berm elevation
                )

                # Community reached minimum width, drowned at 302 years
                cascade_pt45_h2m_low_nourishment_commercial_RT1m = RUN_8_CASCADE_noAST_Rave_SLR_pt004_Nourishment(
                    nt=1000,
                    rmin=0.25,
                    rmax=0.65,  # rave = 0.45
                    name="8-B3D_Rave_pt45_Nourishment_2mDune_lowEle_commercial_rebuild1m",
                    dune_design_elevation=3.6,  # m MHW, keep dune design height the same as 2 m dune above the "roadway"
                    storm_file="StormSeries_1kyrs_VCR_Berm1pt9m_Slope0pt04_01.npy",
                    elevation_file="b3d_pt45_8757yrs_low-elevations.csv",
                    dune_file="barrier3d-default-dunes.npy",
                    overwash_filter=90,  # corresponds with commercial
                    overwash_to_dune=9,
                    nourishment_volume=100,  # m^3/m
                    beach_width_threshold=30,  # m
                    background_erosion=0.0,
                    rebuild_dune_threshold=1,  # m above the berm elevation
                )

                # Community reached minimum width, drowned at 302 years
                cascade_pt45_h2m_low_nourishment_commercial_BE1m_RT1m = RUN_8_CASCADE_noAST_Rave_SLR_pt004_Nourishment(
                    nt=1000,  # will need to run for longer later, after AGU
                    rmin=0.25,
                    rmax=0.65,  # rave = 0.45
                    name="8-B3D_Rave_pt45_Nourishment_2mDune_lowEle_commercial_backerosion1m_rebuild1m",
                    dune_design_elevation=3.6,  # m MHW, keep dune design height the same as 2 m dune above the "roadway"
                    storm_file="StormSeries_1kyrs_VCR_Berm1pt9m_Slope0pt04_01.npy",
                    elevation_file="b3d_pt45_8757yrs_low-elevations.csv",
                    dune_file="barrier3d-default-dunes.npy",
                    overwash_filter=90,  # corresponds with commercial
                    overwash_to_dune=9,
                    nourishment_volume=100,  # m^3/m
                    beach_width_threshold=30,  # m
                    background_erosion=-1.0,  # m/yr, background shoreline erosion
                    rebuild_dune_threshold=1,  # m above the berm elevation
                )

            def pt45_high():
                # Community reached minimum width, drowned at 544 years; Barrier has HEIGHT DROWNED at t = 574 years
                cascade_pt45_h2m_high_nourishment_residential_RT1m = RUN_8_CASCADE_noAST_Rave_SLR_pt004_Nourishment(
                    nt=1000,
                    rmin=0.25,
                    rmax=0.65,  # rave = 0.45
                    name="8-B3D_Rave_pt45_Nourishment_2mDune_highEle_residential_rebuild1m",
                    dune_design_elevation=3.8,  # m MHW, keep dune design height the same as 2 m dune above the "roadway"
                    storm_file="StormSeries_1kyrs_VCR_Berm1pt9m_Slope0pt04_01.npy",
                    elevation_file="b3d_pt45_802yrs_high-elevations.csv",
                    dune_file="barrier3d-default-dunes.npy",
                    overwash_filter=40,
                    overwash_to_dune=9,
                    nourishment_volume=100,  # m^3/m
                    beach_width_threshold=30,  # m
                    background_erosion=0.0,
                    rebuild_dune_threshold=1,  # m above the berm elevation
                )

                # Community reached minimum width, drowned at 429 years
                cascade_pt45_h2m_high_nourishment_commercial_RT1m = RUN_8_CASCADE_noAST_Rave_SLR_pt004_Nourishment(
                    nt=1000,
                    rmin=0.25,
                    rmax=0.65,  # rave = 0.45
                    name="8-B3D_Rave_pt45_Nourishment_2mDune_highEle_commercial_rebuild1m",
                    dune_design_elevation=3.8,  # m MHW, keep dune design height the same as 2 m dune above the "roadway"
                    storm_file="StormSeries_1kyrs_VCR_Berm1pt9m_Slope0pt04_01.npy",
                    elevation_file="b3d_pt45_802yrs_high-elevations.csv",
                    dune_file="barrier3d-default-dunes.npy",
                    overwash_filter=90,  # corresponds with commercial
                    overwash_to_dune=9,
                    nourishment_volume=100,  # m^3/m
                    beach_width_threshold=30,  # m
                    background_erosion=0.0,
                    rebuild_dune_threshold=1,  # m above the berm elevation
                )

                # Community reached minimum width, drowned at 429 years
                cascade_pt45_h2m_high_nourishment_commercial_BE1m_RT1m = RUN_8_CASCADE_noAST_Rave_SLR_pt004_Nourishment(
                    nt=1000,  # will need to run for longer later, after AGU
                    rmin=0.25,
                    rmax=0.65,  # rave = 0.45
                    name="8-B3D_Rave_pt45_Nourishment_2mDune_highEle_commercial_backerosion1m_rebuild1m",
                    dune_design_elevation=3.8,  # m MHW, keep dune design height the same as 2 m dune above the "roadway"
                    storm_file="StormSeries_1kyrs_VCR_Berm1pt9m_Slope0pt04_01.npy",
                    elevation_file="b3d_pt45_802yrs_high-elevations.csv",
                    dune_file="barrier3d-default-dunes.npy",
                    overwash_filter=90,  # corresponds with commercial
                    overwash_to_dune=9,
                    nourishment_volume=100,  # m^3/m
                    beach_width_threshold=30,  # m
                    background_erosion=-1.0,  # m/yr, background shoreline erosion
                    rebuild_dune_threshold=1,  # m above the berm elevation
                )

        def topo_only():
            # we only run 10 years of the following runs because we use them for plotting the initial topo figure for
            # the CNH simulations
            cascade_pt45_h2m_high_nourishment_commercial = RUN_8_CASCADE_noAST_Rave_SLR_pt004_Nourishment(
                nt=10,
                rmin=0.25,
                rmax=0.65,  # rave = 0.45
                name="8-B3D_Rave_pt45_Nourishment_2mDune_highEle_commercial",
                dune_design_elevation=3.8,  # m MHW, keep dune design height the same as 2 m dune above the "roadway"
                storm_file="StormSeries_1kyrs_VCR_Berm1pt9m_Slope0pt04_01.npy",
                elevation_file="b3d_pt45_802yrs_high-elevations.csv",
                dune_file="barrier3d-default-dunes.npy",
                overwash_filter=90,  # corresponds with commercial
                overwash_to_dune=9,
                nourishment_volume=100,  # m^3/m
                beach_width_threshold=30,  # m
                background_erosion=0.0,
                rebuild_dune_threshold=1,
            )

            cascade_pt75_h2m_high_nourishment_commercial = RUN_8_CASCADE_noAST_Rave_SLR_pt004_Nourishment(
                nt=10,
                rmin=0.55,
                rmax=0.95,  # rave = 0.75
                name="8-B3D_Rave_pt75_Nourishment_2mDune_highEle_commercial",
                dune_design_elevation=4.1,  # m MHW, keep dune design height the same as 2 m dune above the "roadway"
                storm_file="StormSeries_1kyrs_VCR_Berm1pt9m_Slope0pt04_01.npy",
                elevation_file="b3d_pt75_829yrs_high-elevations.csv",
                dune_file="barrier3d-default-dunes.npy",
                overwash_filter=90,  # corresponds with commercial
                overwash_to_dune=9,
                nourishment_volume=100,  # m^3/m
                beach_width_threshold=30,  # m
                background_erosion=0.0,
                rebuild_dune_threshold=1,
            )

        def old_versions():

            # note, we keep all other variables the same for comparison to the roadways scenarios, except here, we test
            # the sensitivity of the dune rebuilding threshold: 1) only rebuild if it is totally wiped out (we specify 0.3 m
            # above the berm) or 2) rebuild if the dune is eroded to 1-m above the berm
            def rebuild_threshold_pt3m():
                # Roadway scenario drowned at 162 years
                # Community reached minimum width, drowned at 178 years
                cascade_pt75_h2m_low_nourishment_residential_RTpt3m = RUN_8_CASCADE_noAST_Rave_SLR_pt004_Nourishment(
                    nt=1000,
                    rmin=0.55,
                    rmax=0.95,  # rave = 0.75
                    name="8-B3D_Rave_pt75_Nourishment_2mDune_lowEle_residential_rebuildpt3m",
                    dune_design_elevation=3.2,
                    # m MHW, keep dune design height the same as 2m dune above the initial "roadway"
                    storm_file="StormSeries_1kyrs_VCR_Berm1pt9m_Slope0pt04_01.npy",
                    elevation_file="b3d_pt75_3284yrs_low-elevations.csv",
                    dune_file="barrier3d-default-dunes.npy",
                    overwash_filter=40,  # corresponds with residential
                    overwash_to_dune=9,
                    nourishment_volume=100,  # m^3/m
                    beach_width_threshold=30,  # m
                    background_erosion=0.0,
                    rebuild_dune_threshold=0.3,  # m above the berm elevation
                )

                # Community reached minimum width, drowned at 80 years
                cascade_pt75_h2m_low_nourishment_commercial_RTpt3m = RUN_8_CASCADE_noAST_Rave_SLR_pt004_Nourishment(
                    nt=1000,
                    rmin=0.55,
                    rmax=0.95,  # rave = 0.75
                    name="8-B3D_Rave_pt75_Nourishment_2mDune_lowEle_commercial_rebuildpt3m",
                    dune_design_elevation=3.2,
                    # m MHW, keep dune design height the same as 2 m dune above the "roadway"
                    storm_file="StormSeries_1kyrs_VCR_Berm1pt9m_Slope0pt04_01.npy",
                    elevation_file="b3d_pt75_3284yrs_low-elevations.csv",
                    dune_file="barrier3d-default-dunes.npy",
                    overwash_filter=90,  # corresponds with commercial
                    overwash_to_dune=9,
                    nourishment_volume=100,  # m^3/m
                    beach_width_threshold=30,  # m
                    background_erosion=0.0,
                    rebuild_dune_threshold=0.3,  # m above the berm elevation
                )

                # Community reached minimum width, drowned at 80 years
                cascade_pt75_h2m_low_nourishment_commercial_BE1m_RTpt3m = RUN_8_CASCADE_noAST_Rave_SLR_pt004_Nourishment(
                    nt=1000,
                    rmin=0.55,
                    rmax=0.95,  # rave = 0.75
                    name="8-B3D_Rave_pt75_Nourishment_2mDune_lowEle_commercial_backerosion1m_rebuildpt3m",
                    dune_design_elevation=3.2,
                    # m MHW, keep dune design height the same as 2 m dune above the "roadway"
                    storm_file="StormSeries_1kyrs_VCR_Berm1pt9m_Slope0pt04_01.npy",
                    elevation_file="b3d_pt75_3284yrs_low-elevations.csv",
                    dune_file="barrier3d-default-dunes.npy",
                    overwash_filter=90,  # corresponds with commercial
                    overwash_to_dune=9,
                    nourishment_volume=100,  # m^3/m
                    beach_width_threshold=30,  # m
                    background_erosion=-1.0,  # m/yr, background shoreline erosion
                    rebuild_dune_threshold=0.3,  # m above the berm elevation
                )

                # Roadway scenario drowned at 404 years
                # Community reached minimum width, drowned at 648 years, barrier HEIGHT DROWNED at t = 710 years
                cascade_pt45_h2m_high_nourishment_residential_RTpt3m = RUN_8_CASCADE_noAST_Rave_SLR_pt004_Nourishment(
                    nt=1000,
                    rmin=0.25,
                    rmax=0.65,  # rave = 0.45
                    name="8-B3D_Rave_pt45_Nourishment_2mDune_highEle_residential_rebuildpt3m",
                    dune_design_elevation=3.7,
                    # m MHW, keep dune design height the same as 2 m dune above the "roadway"
                    storm_file="StormSeries_1kyrs_VCR_Berm1pt9m_Slope0pt04_01.npy",
                    elevation_file="b3d_pt45_802yrs_high-elevations.csv",
                    dune_file="barrier3d-default-dunes.npy",
                    overwash_filter=40,
                    overwash_to_dune=9,
                    nourishment_volume=100,  # m^3/m
                    beach_width_threshold=30,  # m
                    background_erosion=0.0,
                    rebuild_dune_threshold=0.3,  # m above the berm elevation
                )

                # Community reached minimum width, drowned at 426 years; Barrier has HEIGHT DROWNED at t = 452 years
                cascade_pt45_h2m_high_nourishment_commercial_RTpt3m = RUN_8_CASCADE_noAST_Rave_SLR_pt004_Nourishment(
                    nt=1000,
                    rmin=0.25,
                    rmax=0.65,  # rave = 0.45
                    name="8-B3D_Rave_pt45_Nourishment_2mDune_highEle_commercial_rebuildpt3m",
                    dune_design_elevation=3.7,
                    # m MHW, keep dune design height the same as 2 m dune above the "roadway"
                    storm_file="StormSeries_1kyrs_VCR_Berm1pt9m_Slope0pt04_01.npy",
                    elevation_file="b3d_pt45_802yrs_high-elevations.csv",
                    dune_file="barrier3d-default-dunes.npy",
                    overwash_filter=90,  # corresponds with commercial
                    overwash_to_dune=9,
                    nourishment_volume=100,  # m^3/m
                    beach_width_threshold=30,  # m
                    background_erosion=0.0,
                    rebuild_dune_threshold=0.3,  # m above the berm elevation
                )

                # # Community reached minimum width, drowned at 426 years
                # cascade_pt45_h2m_high_nourishment_commercial_BEpt25m_RTpt3m = RUN_8_CASCADE_noAST_Rave_SLR_pt004_Nourishment(
                #     nt=1000,
                #     rmin=0.25,
                #     rmax=0.65,  # rave = 0.45
                #     name="8-B3D_Rave_pt45_Nourishment_2mDune_highEle_commercial_backerosionpt25m_rebuildpt3m",
                #     dune_design_elevation=3.7,  # m MHW, keep dune design height the same as 2 m dune above the "roadway"
                #     storm_file="StormSeries_1kyrs_VCR_Berm1pt9m_Slope0pt04_01.npy",
                #     elevation_file="b3d_pt45_802yrs_high-elevations.csv",
                #     dune_file="barrier3d-default-dunes.npy",
                #     overwash_filter=90,  # corresponds with commercial
                #     overwash_to_dune=9,
                #     nourishment_volume=100,  # m^3/m
                #     beach_width_threshold=30,  # m
                #     background_erosion=-0.25,  # m/yr, background shoreline erosion
                #     rebuild_dune_threshold=0.3,  # m above the berm elevation
                # )

                # Community reached minimum width, drowned at 426 years
                cascade_pt45_h2m_high_nourishment_commercial_BE1m_RTpt3m = RUN_8_CASCADE_noAST_Rave_SLR_pt004_Nourishment(
                    nt=1000,
                    rmin=0.25,
                    rmax=0.65,  # rave = 0.45
                    name="8-B3D_Rave_pt45_Nourishment_2mDune_highEle_commercial_backerosion1m_rebuildpt3m",
                    dune_design_elevation=3.7,
                    # m MHW, keep dune design height the same as 2 m dune above the "roadway"
                    storm_file="StormSeries_1kyrs_VCR_Berm1pt9m_Slope0pt04_01.npy",
                    elevation_file="b3d_pt45_802yrs_high-elevations.csv",
                    dune_file="barrier3d-default-dunes.npy",
                    overwash_filter=90,  # corresponds with commercial
                    overwash_to_dune=9,
                    nourishment_volume=100,  # m^3/m
                    beach_width_threshold=30,  # m
                    background_erosion=-1.0,  # m/yr, background shoreline erosion
                    rebuild_dune_threshold=0.3,  # m above the berm elevation
                )

                def old():
                    # roadway scenario drowned at 404 years
                    # Community reached minimum width, drowned at 496 years
                    cascade_pt45_h2m_low_nourishment_residential = RUN_8_CASCADE_noAST_Rave_SLR_pt004_Nourishment(
                        nt=500,
                        rmin=0.25,
                        rmax=0.65,  # rave = 0.45
                        name="8-B3D_Rave_pt45_Nourishment_2mDune_lowEle_residential",
                        dune_design_elevation=3.7,
                        # m MHW, keep dune design height the same as 2 m dune above the "roadway"
                        storm_file="StormSeries_1kyrs_VCR_Berm1pt9m_Slope0pt04_01.npy",
                        elevation_file="b3d_pt45_8750yrs_low-elevations.csv",
                        dune_file="barrier3d-default-dunes.npy",
                        overwash_filter=40,  # corresponds with commercial
                        nourishment_volume=100,  # m^3/m
                        beach_width_threshold=30,  # m
                        background_erosion=0.0,
                    )

                    # Community reached minimum width, drowned at 421 years
                    # Barrier has HEIGHT DROWNED at t = 458 years
                    cascade_pt45_h2m_low_nourishment_commercial = RUN_8_CASCADE_noAST_Rave_SLR_pt004_Nourishment(
                        nt=500,
                        rmin=0.25,
                        rmax=0.65,  # rave = 0.45
                        name="8-B3D_Rave_pt45_Nourishment_2mDune_lowEle_commercial",
                        dune_design_elevation=3.7,
                        # m MHW, keep dune design height the same as 2 m dune above the "roadway"
                        storm_file="StormSeries_1kyrs_VCR_Berm1pt9m_Slope0pt04_01.npy",
                        elevation_file="b3d_pt45_8750yrs_low-elevations.csv",
                        dune_file="barrier3d-default-dunes.npy",
                        overwash_filter=90,  # corresponds with commercial
                        nourishment_volume=100,  # m^3/m
                        beach_width_threshold=30,  # m
                        background_erosion=0.0,
                    )

                    # Community reached minimum width, drowned at 421 years
                    # Barrier has HEIGHT DROWNED at t = 454 years
                    cascade_pt45_h2m_low_nourishment_commercial_background_erosion_pt25m = RUN_8_CASCADE_noAST_Rave_SLR_pt004_Nourishment(
                        nt=500,
                        rmin=0.25,
                        rmax=0.65,  # rave = 0.45
                        name="8-B3D_Rave_pt45_Nourishment_2mDune_lowEle_commercial_backerosion_pt25m",
                        dune_design_elevation=3.7,
                        # m MHW, keep dune design height the same as 2 m dune above the "roadway"
                        storm_file="StormSeries_1kyrs_VCR_Berm1pt9m_Slope0pt04_01.npy",
                        elevation_file="b3d_pt45_8750yrs_low-elevations.csv",
                        dune_file="barrier3d-default-dunes.npy",
                        overwash_filter=90,  # corresponds with commercial
                        nourishment_volume=100,  # m^3/m
                        beach_width_threshold=30,  # m
                        background_erosion=-0.25,  # m/yr, background shoreline erosion
                    )

                    # Community reached minimum width, drowned at 421 years
                    cascade_pt45_h2m_low_nourishment_commercial_background_erosion_1m = RUN_8_CASCADE_noAST_Rave_SLR_pt004_Nourishment(
                        nt=500,
                        rmin=0.25,
                        rmax=0.65,  # rave = 0.45
                        name="8-B3D_Rave_pt45_Nourishment_2mDune_lowEle_commercial_backerosion_1m",
                        dune_design_elevation=3.7,
                        # m MHW, keep dune design height the same as 2 m dune above the "roadway"
                        storm_file="StormSeries_1kyrs_VCR_Berm1pt9m_Slope0pt04_01.npy",
                        elevation_file="b3d_pt45_8750yrs_low-elevations.csv",
                        dune_file="barrier3d-default-dunes.npy",
                        overwash_filter=90,  # corresponds with commercial
                        nourishment_volume=100,  # m^3/m
                        beach_width_threshold=30,  # m
                        background_erosion=-1.0,  # m/yr, background shoreline erosion
                    )

    def alongshore_variable_management():
        def old_runs():
            def nourishment_pt75_low():
                # these initial conditions drowned at 80 years in nourishments
                number_barrier3d_models = 6
                beach_width_threshold = [30] * number_barrier3d_models
                rmin = [0.55] * number_barrier3d_models
                rmax = [0.95] * number_barrier3d_models
                elevation_file = [
                    "b3d_pt75_3284yrs_low-elevations.csv"
                ] * number_barrier3d_models
                dune_file = ["barrier3d-default-dunes.npy"] * number_barrier3d_models
                dune_design_elevation = [3.2] * number_barrier3d_models
                roads_on = [False] * number_barrier3d_models
                nourishments_on = [True] * number_barrier3d_models

                # all B3D segments drown at 80
                nourishment_only_6AST_low_pt75_comm_BE1m_RT1m = RUN_9_CASCADE_Rave_SLR_pt004_AlongshoreVariableManagement(
                    nt=100,
                    name="9-CASCADE_Rave_pt75_Nourishment_2mDune_lowEle_comm_BE1m_RT1m_6AST",
                    storm_file="StormSeries_1kyrs_VCR_Berm1pt9m_Slope0pt04_01.npy",
                    alongshore_section_count=number_barrier3d_models,  # NOTE: will want to go back to sensitivity modeling
                    num_cores=6,  # for my laptop, max is ?
                    beach_width_threshold=beach_width_threshold,  # m
                    rmin=rmin,
                    rmax=rmax,  # rave = 0.75
                    elevation_file=elevation_file,
                    dune_file=dune_file,
                    dune_design_elevation=dune_design_elevation,
                    dune_minimum_elevation=None,
                    road_ele=None,
                    road_width=None,
                    road_setback=None,
                    overwash_filter=90,  # corresponds with commercial
                    overwash_to_dune=9,
                    nourishment_volume=100,  # m^3/m
                    background_erosion=-1.00,  # m/yr, background shoreline erosion
                    rebuild_dune_threshold=1,  # m above the berm elevation
                    roadway_management_on=roads_on,
                    beach_dune_manager_on=nourishments_on,
                )

            def nourishment_pt45_high_RT1m():
                number_barrier3d_models = 6
                beach_width_threshold = [30] * number_barrier3d_models
                rmin = [0.25] * number_barrier3d_models
                rmax = [0.65] * number_barrier3d_models
                elevation_file = [
                    "b3d_pt45_802yrs_high-elevations.csv"
                ] * number_barrier3d_models
                dune_file = ["barrier3d-default-dunes.npy"] * number_barrier3d_models
                dune_design_elevation = [3.2] * number_barrier3d_models
                roads_on = [False] * number_barrier3d_models
                nourishments_on = [True] * number_barrier3d_models

                nourishment_only_6AST_low_pt75_comm_BE1m_RT1m = RUN_9_CASCADE_Rave_SLR_pt004_AlongshoreVariableManagement(
                    nt=100,
                    name="9-CASCADE_Rave_pt45_Nourishment_2mDune_highEle_res_BE1m_RT1m_6AST",
                    storm_file="StormSeries_1kyrs_VCR_Berm1pt9m_Slope0pt04_01.npy",
                    alongshore_section_count=number_barrier3d_models,  # NOTE: will want to go back to sensitivity modeling
                    num_cores=6,  # for my laptop, max is ?
                    beach_width_threshold=beach_width_threshold,  # m
                    rmin=rmin,
                    rmax=rmax,  # rave = 0.45
                    elevation_file=elevation_file,
                    dune_file=dune_file,
                    dune_design_elevation=dune_design_elevation,
                    dune_minimum_elevation=None,
                    road_ele=None,
                    road_width=None,
                    road_setback=None,
                    overwash_filter=90,  # corresponds with commercial
                    overwash_to_dune=9,
                    nourishment_volume=100,  # m^3/m
                    background_erosion=-1.00,  # m/yr, background shoreline erosion
                    rebuild_dune_threshold=1,  # m above the berm elevation, more realistic is 1 m
                    roadway_management_on=roads_on,
                    beach_dune_manager_on=nourishments_on,
                )

            def nourishment_pt45_high_RTpt3m():
                # nourishments only community drowned at 426 years; Barrier has HEIGHT DROWNED at t = 452 years
                number_barrier3d_models = 6
                beach_width_threshold = [30] * number_barrier3d_models
                rmin = [0.25] * number_barrier3d_models
                rmax = [0.65] * number_barrier3d_models
                elevation_file = [
                    "b3d_pt45_802yrs_high-elevations.csv"
                ] * number_barrier3d_models
                dune_file = ["barrier3d-default-dunes.npy"] * number_barrier3d_models
                dune_design_elevation = [3.2] * number_barrier3d_models
                roads_on = [False] * number_barrier3d_models
                nourishments_on = [True] * number_barrier3d_models

                # Community reached minimum width, drowned at 485 years
                nourishment_only_6AST_low_pt75_comm_BE1m_RT1m = RUN_9_CASCADE_Rave_SLR_pt004_AlongshoreVariableManagement(
                    nt=500,
                    name="9-CASCADE_Rave_pt45_Nourishment_2mDune_highEle_res_BE1m_RT1m_6AST",
                    storm_file="StormSeries_1kyrs_VCR_Berm1pt9m_Slope0pt04_01.npy",
                    alongshore_section_count=number_barrier3d_models,  # NOTE: will want to go back to sensitivity modeling
                    num_cores=6,  # for my laptop, max is ?
                    beach_width_threshold=beach_width_threshold,  # m
                    rmin=rmin,
                    rmax=rmax,  # rave = 0.45
                    elevation_file=elevation_file,
                    dune_file=dune_file,
                    dune_design_elevation=dune_design_elevation,
                    dune_minimum_elevation=None,
                    road_ele=None,
                    road_width=None,
                    road_setback=None,
                    overwash_filter=90,  # corresponds with commercial
                    overwash_to_dune=9,
                    nourishment_volume=100,  # m^3/m
                    background_erosion=0,  # m/yr, background shoreline erosion
                    rebuild_dune_threshold=0.3,  # m above the berm elevation
                    roadway_management_on=roads_on,
                    beach_dune_manager_on=nourishments_on,
                )

            def nourishment_pt75_low_split_natural():
                # these initial conditions drowned at 80 years in nourishments
                number_barrier3d_models = 6
                beach_width_threshold = [30] * number_barrier3d_models
                rmin = [0.55] * number_barrier3d_models
                rmax = [0.95] * number_barrier3d_models
                elevation_file = [
                    "b3d_pt75_3284yrs_low-elevations.csv"
                ] * number_barrier3d_models
                dune_file = ["barrier3d-default-dunes.npy"] * number_barrier3d_models
                dune_design_elevation = [3.2] * number_barrier3d_models
                roads_on = [False] * number_barrier3d_models
                nourishments_on = [True, True, True, False, False, False]

                # very quickly these barriers separate -- MAKES SENSE, but not interesting to show
                nourish_natural_6AST_low_pt75_comm_BE1m_RT1m = RUN_9_CASCADE_Rave_SLR_pt004_AlongshoreVariableManagement(
                    nt=100,
                    name="9-CASCADE_Rave_pt75_2mDune_lowEle_comm_BE1m_RT1m_6AST_nourish_nat",
                    storm_file="StormSeries_1kyrs_VCR_Berm1pt9m_Slope0pt04_01.npy",
                    alongshore_section_count=number_barrier3d_models,  # NOTE: will want to go back to sensitivity modeling
                    num_cores=6,  # for my laptop, max is ?
                    beach_width_threshold=beach_width_threshold,  # m
                    rmin=rmin,
                    rmax=rmax,  # rave = 0.75
                    elevation_file=elevation_file,
                    dune_file=dune_file,
                    dune_design_elevation=dune_design_elevation,
                    dune_minimum_elevation=None,
                    road_ele=None,
                    road_width=None,
                    road_setback=None,
                    overwash_filter=90,  # corresponds with commercial
                    overwash_to_dune=9,
                    nourishment_volume=100,  # m^3/m
                    background_erosion=-1.00,  # m/yr, background shoreline erosion
                    rebuild_dune_threshold=1,  # m above the berm elevation
                    roadway_management_on=roads_on,
                    beach_dune_manager_on=nourishments_on,
                )

        def nourishment_pt75_low_split_roadways_BE1m_3km():

            # variables that DO NOT change among runs
            number_barrier3d_models = 6
            beach_width_threshold = [30] * number_barrier3d_models
            rmin = [0.55] * number_barrier3d_models
            rmax = [0.95] * number_barrier3d_models
            elevation_file = [
                "b3d_pt75_4261yrs_low-elevations.csv"
            ] * number_barrier3d_models
            dune_file = ["barrier3d-default-dunes.npy"] * number_barrier3d_models
            storm_file = "StormSeries_1kyrs_VCR_Berm1pt9m_Slope0pt04_01.npy"
            dune_design_elevation = [2.6] * number_barrier3d_models  # 2 m scenario
            num_cores = 6  # for my laptop, max is 15
            dune_minimum_elevation = 1.1  # m MHW, allow dune to erode down to 0.5 m above the roadway, for roadways only
            road_ele = 0.6  # m MHW
            road_width = 20  # m
            road_setback = 20  # m
            overwash_filter = 90  # commercial
            overwash_to_dune = 9
            nourishment_volume = 100  # m^3/m
            background_erosion = -1.0  # m/yr, background shoreline erosion
            rebuild_dune_threshold = 1  # m

            # baseline models for comparison -- all roadways ----------------------------------------
            roads_on = [True, True, True, True, True, True]
            nourishments_on = [False, False, False, False, False, False]
            sea_level_rise_rate = 0.004
            sea_level_constant = True  # linear

            # Island is too narrow for roadway to be relocated. Roadway eaten up by dunes at 73 years
            roadways_6AST_low_pt75_BE1m = (
                RUN_9_CASCADE_Rave_SLR_pt004_AlongshoreVariableManagement(
                    nt=200,
                    name="9-CASCADE_Rave_pt75_2mDune_lowEle_BE1m_6AST_6roads",
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
                )
            )

            # baseline models for comparison -- all nourishments ----------------------------------------
            roads_on = [False, False, False, False, False, False]
            nourishments_on = [True, True, True, True, True, True]
            sea_level_rise_rate = 0.004
            sea_level_constant = True  # linear

            # Community reached minimum width, drowned at 83 years
            community_6AST_low_pt75_BE1m = RUN_9_CASCADE_Rave_SLR_pt004_AlongshoreVariableManagement(
                nt=200,
                name="9-CASCADE_Rave_pt75_2mDune_lowEle_comm_BE1m_RT1m_6AST_6nourish",
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
            )

            # split management - 1 m background erosion and 0.004 m/yr SLR ----------------------------------------
            roads_on = [False, False, False, True, True, True]
            nourishments_on = [True, True, True, False, False, False]
            sea_level_rise_rate = 0.004
            sea_level_constant = True  # linear

            # All 3: Community reached minimum width, drowned at 83 years
            # All 3: Island is too narrow for roadway to be relocated. Roadway eaten up by dunes at 143 years
            nourishment_roadways_6AST_low_pt75_comm_BE1m_RT1m = RUN_9_CASCADE_Rave_SLR_pt004_AlongshoreVariableManagement(
                nt=200,
                name="9-CASCADE_Rave_pt75_Nourish_2mDune_lowEle_comm_BE1m_RT1m_6AST_3roads",
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
            )

            # split management - 1 m background erosion and accelerated SLR ----------------------------------------
            roads_on = [False, False, False, True, True, True]
            nourishments_on = [True, True, True, False, False, False]
            sea_level_rise_rate = 0.004  # dummy here
            sea_level_constant = False  # accelerated

            # All 3: Community reached minimum width, drowned at 64 years
            # All 3: Barrier has HEIGHT DROWNED at t = 85 years
            nourishment_roadways_6AST_low_pt75_comm_BE1m_RT1m = RUN_9_CASCADE_Rave_SLR_pt004_AlongshoreVariableManagement(
                nt=200,
                name="9-CASCADE_Rave_pt75_Nourish_2mDune_lowEle_comm_BE1m_RT1m_6AST_3roads_AccSLR",
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
            )

        def nourishment_pt45_low_left_roadways_pt75_low_right_BE1m_3km():

            # variables that DO NOT change among runs
            number_barrier3d_models = 6
            beach_width_threshold = [30] * number_barrier3d_models
            rmin = [0.25] * 3 + [0.55] * 3
            rmax = [0.65] * 3 + [0.95] * 3
            elevation_file = ["b3d_pt45_8757yrs_low-elevations.csv"] * 3 + [
                "b3d_pt75_4261yrs_low-elevations.csv"
            ] * 3
            dune_file = ["barrier3d-default-dunes.npy"] * number_barrier3d_models
            storm_file = "StormSeries_1kyrs_VCR_Berm1pt9m_Slope0pt04_01.npy"
            dune_design_elevation = [3.6] * 3 + [2.6] * 3  # 2 m above the roadway
            dune_minimum_elevation = [2.1] * 3 + [
                1.1
            ] * 3  # m MHW, allow dune to erode down to 0.5 m above the roadway, for roadways only (others dummy)
            road_ele = [1.6] * 3 + [
                0.6
            ] * 3  # first 3 are dummys since we are only doing nourishment there
            num_cores = 6  # for my laptop, max is 15
            road_width = 20  # m
            road_setback = 20  # m
            overwash_filter = 90  # commercial
            overwash_to_dune = 9
            nourishment_volume = 100  # m^3/m
            background_erosion = -1.0  # m/yr, background shoreline erosion
            rebuild_dune_threshold = 1  # m, for nourishments only

            # split management - 1 m background erosion and 0.004 m/yr SLR ----------------------------------------
            roads_on = [False, False, False, True, True, True]
            nourishments_on = [True, True, True, False, False, False]
            sea_level_rise_rate = 0.004
            sea_level_constant = True  # linear

            # All 3: Roadway drowned in place at 132 years due to SLR - road cannot be below 0 m MHW
            nourishment_roadways_6AST_low_pt45_pt75_comm_BE1m_RT1m = RUN_9_CASCADE_Rave_SLR_pt004_AlongshoreVariableManagement(
                nt=200,
                name="9-CASCADE_Rave_pt45_pt75_low_split_2mDune_comm_BE1m_RT1m_6AST_3roads",
                storm_file=storm_file,
                alongshore_section_count=number_barrier3d_models,
                num_cores=num_cores,
                beach_width_threshold=beach_width_threshold,
                rmin=rmin,
                rmax=rmax,
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
            )

            # split management - add accelerated SLR ----------------------------------------
            roads_on = [False, False, False, True, True, True]
            nourishments_on = [True, True, True, False, False, False]
            sea_level_rise_rate = 0.004  # dummy here
            sea_level_constant = False  # accelerated

            # All three: Roadway drowned in place at 85 years due to SLR - road cannot be below 0 m MHW
            # All three: Community reached minimum width, drowned at 137 years
            # Barrier has HEIGHT DROWNED at t = 141 years
            nourishment_roadways_6AST_low_pt75_comm_BE1m_RT1m = RUN_9_CASCADE_Rave_SLR_pt004_AlongshoreVariableManagement(
                nt=200,
                name="9-CASCADE_Rave_pt45_pt75_low_split_2mDune_comm_BE1m_RT1m_6AST_3roads_AccSLR",
                storm_file=storm_file,
                alongshore_section_count=number_barrier3d_models,
                num_cores=num_cores,
                beach_width_threshold=beach_width_threshold,
                rmin=rmin,
                rmax=rmax,
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
            )

            # baseline pt45 low -- all nourishment ----------------------------------------
            number_barrier3d_models = 6
            beach_width_threshold = [30] * number_barrier3d_models
            rmin = [0.25] * number_barrier3d_models
            rmax = [0.65] * number_barrier3d_models
            elevation_file = [
                "b3d_pt45_8757yrs_low-elevations.csv"
            ] * number_barrier3d_models
            dune_file = ["barrier3d-default-dunes.npy"] * number_barrier3d_models
            storm_file = "StormSeries_1kyrs_VCR_Berm1pt9m_Slope0pt04_01.npy"
            dune_design_elevation = [
                3.6
            ] * number_barrier3d_models  # 2 m above the roadway
            dune_minimum_elevation = [2.1] * number_barrier3d_models  # dummy
            road_ele = [1.6] * number_barrier3d_models  # dummy
            num_cores = 6  # for my laptop, max is 15
            road_width = 20  # m, dummy
            road_setback = 20  # m, dummy
            overwash_filter = 90  # commercial
            overwash_to_dune = 9
            nourishment_volume = 100  # m^3/m
            background_erosion = -1.0  # m/yr, background shoreline erosion
            rebuild_dune_threshold = 1  # m, for nourishments only

            roads_on = [False, False, False, False, False, False]
            nourishments_on = [True, True, True, True, True, True]
            sea_level_rise_rate = 0.004
            sea_level_constant = True  # linear

            nourishment_6AST_low_pt45_comm_BE1m_RT1m = RUN_9_CASCADE_Rave_SLR_pt004_AlongshoreVariableManagement(
                nt=200,
                name="9-CASCADE_Rave_pt45_low_nourishment_2mDune_comm_BE1m_RT1m_6AST",
                storm_file=storm_file,
                alongshore_section_count=number_barrier3d_models,
                num_cores=num_cores,
                beach_width_threshold=beach_width_threshold,
                rmin=rmin,
                rmax=rmax,
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
            )

        def nourishment_pt45_low_left_middle_pt75_low_roadways_right_pt45_low_roadways_4pt5km():
            # this new run: better illustrates what abandoning in the middle of a domain looks like (Rodanthe)
            # --- natural dune dynamics don't really matter here b/s we are managing all of them ---
            # left: low but wide barrier -- nourishments (pt45 low)
            # middle: low but narrow barrier -- roadways (pt75 low)
            # right: low but wide barrier -- roadways (pt45 low)

            # variables that DO NOT change among runs
            # (NOTE: these variables the same as above -- we maintain a 2 m dune)
            number_barrier3d_models = 9
            beach_width_threshold = [30] * number_barrier3d_models
            rmin = [0.25] * 3 + [0.55] * 3 + [0.25] * 3
            rmax = [0.65] * 3 + [0.95] * 3 + [0.65] * 3
            elevation_file = (
                ["b3d_pt45_8757yrs_low-elevations.csv"] * 3
                + ["b3d_pt75_4261yrs_low-elevations.csv"] * 3
                + ["b3d_pt45_8757yrs_low-elevations.csv"] * 3
            )
            dune_file = ["barrier3d-default-dunes.npy"] * number_barrier3d_models
            storm_file = "StormSeries_1kyrs_VCR_Berm1pt9m_Slope0pt04_01.npy"
            dune_design_elevation = (
                [3.6] * 3 + [2.6] * 3 + [3.6] * 3
            )  # 2 m above the original roadway
            dune_minimum_elevation = (
                [2.1] * 3 + [1.1] * 3 + [2.1] * 3
            )  # m MHW, allow dune to erode down to 0.5 m above the roadway, roadways only (others dummy)
            road_ele = (
                [1.6] * 3 + [0.6] * 3 + [1.6] * 3
            )  # first 3 are dummys since we are only doing nourishment there
            num_cores = 9  # for my laptop, max is 15
            road_width = 20  # m
            road_setback = 20  # m
            overwash_filter = 90  # commercial
            overwash_to_dune = 9
            nourishment_volume = 100  # m^3/m
            background_erosion = -1.0  # m/yr, background shoreline erosion
            rebuild_dune_threshold = 1  # m
            roads_on = [False, False, False, True, True, True, True, True, True]
            nourishments_on = [
                True,
                True,
                True,
                False,
                False,
                False,
                False,
                False,
                False,
            ]

            # 1 m background erosion and 0.004 m/yr SLR ----------------------------------------
            sea_level_rise_rate = 0.004
            sea_level_constant = True  # linear

            # 1353, Barrier has HEIGHT DROWNED at t = 132
            # 6490, Barrier has HEIGHT DROWNED at t = 147 years
            # 3284 is wider and higher (190 m wide vs. 144 m for 4261): no roadway drowned
            # *** 4261, Roadway drowned at 99, 115, 131 due to SLR, road cannot be below 0 m MHW -- created curved shape
            # After I grouped roadway abandonment, all three in the middle are abandoned at 99 years
            AST_3domains_BE1m = RUN_9_CASCADE_Rave_SLR_pt004_AlongshoreVariableManagement(
                nt=200,
                name="9-CASCADE_AST_3domains_BE1m",
                storm_file=storm_file,
                alongshore_section_count=number_barrier3d_models,
                num_cores=num_cores,
                beach_width_threshold=beach_width_threshold,
                rmin=rmin,
                rmax=rmax,
                elevation_file=elevation_file,
                dune_file=dune_file,
                dune_design_elevation=dune_design_elevation,
                dune_minimum_elevation=dune_minimum_elevation,
                road_ele=road_ele,
                road_width=road_width,
                road_setback=road_setback,
                trigger_dune_knockdown=False,  # this didn't really do anything due to the timing of storms
                group_roadway_abandonment=[0, 0, 0, 1, 1, 1, 2, 2, 2],
                overwash_filter=overwash_filter,
                overwash_to_dune=overwash_to_dune,
                nourishment_volume=nourishment_volume,
                background_erosion=background_erosion,
                rebuild_dune_threshold=rebuild_dune_threshold,
                roadway_management_on=roads_on,
                beach_dune_manager_on=nourishments_on,
                sea_level_rise_rate=sea_level_rise_rate,
                sea_level_constant=sea_level_constant,
            )

            # 1 m background erosion and accelerated SLR ----------------------------------------
            sea_level_rise_rate = 0.004  # dummy
            sea_level_constant = False  # accelerated

            # Barrier has HEIGHT DROWNED at t = 71 years (#5 B3D) - 4261
            # even if I change the dune design height to <1 m, this barrier wants to drown due to Acc SLR
            AST_3domains_AccSLR_BE1m = (
                RUN_9_CASCADE_Rave_SLR_pt004_AlongshoreVariableManagement(
                    nt=200,
                    name="9-CASCADE_AST_3domains_BE1m_AccSLR",
                    storm_file=storm_file,
                    alongshore_section_count=number_barrier3d_models,
                    num_cores=num_cores,
                    beach_width_threshold=beach_width_threshold,
                    rmin=rmin,
                    rmax=rmax,
                    elevation_file=elevation_file,
                    dune_file=dune_file,
                    dune_design_elevation=dune_design_elevation,
                    dune_minimum_elevation=dune_minimum_elevation,
                    road_ele=road_ele,
                    road_width=road_width,
                    road_setback=road_setback,
                    group_roadway_abandonment=[0, 0, 0, 1, 1, 1, 2, 2, 2],
                    overwash_filter=overwash_filter,
                    overwash_to_dune=overwash_to_dune,
                    nourishment_volume=nourishment_volume,
                    background_erosion=background_erosion,
                    rebuild_dune_threshold=rebuild_dune_threshold,
                    roadway_management_on=roads_on,
                    beach_dune_manager_on=nourishments_on,
                    sea_level_rise_rate=sea_level_rise_rate,
                    sea_level_constant=sea_level_constant,
                )
            )

            # NO background erosion and accelerated SLR (just talk about this one) ----------------------------
            sea_level_rise_rate = 0.004  # dummy
            sea_level_constant = False  # accelerated
            background_erosion = 0.0

            # Barrier has HEIGHT DROWNED at t = 92 years
            AST_3domains_AccSLR = (
                RUN_9_CASCADE_Rave_SLR_pt004_AlongshoreVariableManagement(
                    nt=200,
                    name="9-CASCADE_AST_3domains_AccSLR",
                    storm_file=storm_file,
                    alongshore_section_count=number_barrier3d_models,
                    num_cores=num_cores,
                    beach_width_threshold=beach_width_threshold,
                    rmin=rmin,
                    rmax=rmax,
                    elevation_file=elevation_file,
                    dune_file=dune_file,
                    dune_design_elevation=dune_design_elevation,
                    dune_minimum_elevation=dune_minimum_elevation,
                    road_ele=road_ele,
                    road_width=road_width,
                    road_setback=road_setback,
                    group_roadway_abandonment=[0, 0, 0, 1, 1, 1, 2, 2, 2],
                    overwash_filter=overwash_filter,
                    overwash_to_dune=overwash_to_dune,
                    nourishment_volume=nourishment_volume,
                    background_erosion=background_erosion,
                    rebuild_dune_threshold=rebuild_dune_threshold,
                    roadway_management_on=roads_on,
                    beach_dune_manager_on=nourishments_on,
                    sea_level_rise_rate=sea_level_rise_rate,
                    sea_level_constant=sea_level_constant,
                )
            )

            # 1 m background erosion and accelerated SLR, middle natural scenario --------------------------------------
            sea_level_rise_rate = 0.004  # dummy
            sea_level_constant = False  # accelerated

            # set middle to no management and lets see what happens
            roads_on = [False, False, False, False, False, False, True, True, True]

            # Roadway width drowned at 137 years, 20.0% of road borders water
            # All: Community reached minimum width, drowned at 137 years
            # Roadway width drowned at 142 years, 20.0% of road borders water
            # Island is too narrow for roadway to be relocated. Roadway eaten up by dunes at 149 years
            # ---- With new roadway abandonment grouping, all drown at 137 (roadway and community) ----
            AST_3domains_AccSLR_BE1m = (
                RUN_9_CASCADE_Rave_SLR_pt004_AlongshoreVariableManagement(
                    nt=200,
                    name="9-CASCADE_AST_3domains_BE1m_AccSLR_nat_middle",
                    storm_file=storm_file,
                    alongshore_section_count=number_barrier3d_models,
                    num_cores=num_cores,
                    beach_width_threshold=beach_width_threshold,
                    rmin=rmin,
                    rmax=rmax,
                    elevation_file=elevation_file,
                    dune_file=dune_file,
                    dune_design_elevation=dune_design_elevation,
                    dune_minimum_elevation=dune_minimum_elevation,
                    road_ele=road_ele,
                    road_width=road_width,
                    road_setback=road_setback,
                    group_roadway_abandonment=[0, 0, 0, 0, 0, 0, 1, 1, 1],
                    overwash_filter=overwash_filter,
                    overwash_to_dune=overwash_to_dune,
                    nourishment_volume=nourishment_volume,
                    background_erosion=background_erosion,
                    rebuild_dune_threshold=rebuild_dune_threshold,
                    roadway_management_on=roads_on,
                    beach_dune_manager_on=nourishments_on,
                    sea_level_rise_rate=sea_level_rise_rate,
                    sea_level_constant=sea_level_constant,
                )
            )

        def averages():
            def one_hundred_thirds_acc_BE1m_ast_runs(
                storm_start=0,
                storm_end=100,
                name_prefix="9-CASCADE_AST_3domains_BE1m_AccSLR",
                storm_prefix="StormSeries_1kyrs_VCR_Berm1pt9m_Slope0pt04_",
            ):

                # variables that DO NOT change among runs
                # (NOTE: these variables the same as above -- we maintain a 2 m dune)
                number_barrier3d_models = 9
                beach_width_threshold = [30] * number_barrier3d_models
                rmin = [0.25] * 3 + [0.55] * 3 + [0.25] * 3
                rmax = [0.65] * 3 + [0.95] * 3 + [0.65] * 3
                elevation_file = (
                    ["b3d_pt45_8757yrs_low-elevations.csv"] * 3
                    + ["b3d_pt75_4261yrs_low-elevations.csv"] * 3
                    + ["b3d_pt45_8757yrs_low-elevations.csv"] * 3
                )
                dune_file = ["barrier3d-default-dunes.npy"] * number_barrier3d_models
                dune_design_elevation = (
                    [3.6] * 3 + [2.6] * 3 + [3.6] * 3
                )  # 2 m above the original roadway
                dune_minimum_elevation = (
                    [2.1] * 3 + [1.1] * 3 + [2.1] * 3
                )  # m MHW, allow dune to erode down to 0.5 m above the roadway, roadways only (others dummy)
                road_ele = (
                    [1.6] * 3 + [0.6] * 3 + [1.6] * 3
                )  # first 3 are dummys since we are only doing nourishment there
                num_cores = 9
                road_width = 20  # m
                road_setback = 20  # m
                overwash_filter = 90  # commercial
                overwash_to_dune = 9
                nourishment_volume = 100  # m^3/m
                background_erosion = -1.0  # m/yr, background shoreline erosion
                rebuild_dune_threshold = 1  # m
                nourishments_on = [
                    True,
                    True,
                    True,
                    False,
                    False,
                    False,
                    False,
                    False,
                    False,
                ]

                sea_level_rise_rate = 0.004  # dummy
                sea_level_constant = False  # accelerated

                for iStorm in range(storm_start, storm_end):
                    name = name_prefix + str(iStorm)
                    storm_file = storm_prefix + str(iStorm) + ".npy"

                    roads_on = [False, False, False, True, True, True, True, True, True]

                    RUN_9_CASCADE_Rave_SLR_pt004_AlongshoreVariableManagement(
                        nt=200,
                        name=name,
                        storm_file=storm_file,
                        alongshore_section_count=number_barrier3d_models,
                        num_cores=num_cores,
                        beach_width_threshold=beach_width_threshold,
                        rmin=rmin,
                        rmax=rmax,
                        elevation_file=elevation_file,
                        dune_file=dune_file,
                        dune_design_elevation=dune_design_elevation,
                        dune_minimum_elevation=dune_minimum_elevation,
                        road_ele=road_ele,
                        road_width=road_width,
                        road_setback=road_setback,
                        group_roadway_abandonment=[0, 0, 0, 1, 1, 1, 2, 2, 2],
                        overwash_filter=overwash_filter,
                        overwash_to_dune=overwash_to_dune,
                        nourishment_volume=nourishment_volume,
                        background_erosion=background_erosion,
                        rebuild_dune_threshold=rebuild_dune_threshold,
                        roadway_management_on=roads_on,
                        beach_dune_manager_on=nourishments_on,
                        sea_level_rise_rate=sea_level_rise_rate,
                        sea_level_constant=sea_level_constant,
                    )

                    # set middle to no management and lets see what happens
                    roads_on = [
                        False,
                        False,
                        False,
                        False,
                        False,
                        False,
                        True,
                        True,
                        True,
                    ]

                    name_prefix = "9-CASCADE_AST_3domains_BE1m_AccSLR_nat_middle"
                    name = name_prefix + str(iStorm)

                    RUN_9_CASCADE_Rave_SLR_pt004_AlongshoreVariableManagement(
                        nt=200,
                        name=name,
                        storm_file=storm_file,
                        alongshore_section_count=number_barrier3d_models,
                        num_cores=num_cores,
                        beach_width_threshold=beach_width_threshold,
                        rmin=rmin,
                        rmax=rmax,
                        elevation_file=elevation_file,
                        dune_file=dune_file,
                        dune_design_elevation=dune_design_elevation,
                        dune_minimum_elevation=dune_minimum_elevation,
                        road_ele=road_ele,
                        road_width=road_width,
                        road_setback=road_setback,
                        group_roadway_abandonment=[0, 0, 0, 0, 0, 0, 1, 1, 1],
                        overwash_filter=overwash_filter,
                        overwash_to_dune=overwash_to_dune,
                        nourishment_volume=nourishment_volume,
                        background_erosion=background_erosion,
                        rebuild_dune_threshold=rebuild_dune_threshold,
                        roadway_management_on=roads_on,
                        beach_dune_manager_on=nourishments_on,
                        sea_level_rise_rate=sea_level_rise_rate,
                        sea_level_constant=sea_level_constant,
                    )

            one_hundred_thirds_acc_BE1m_ast_runs(
                storm_start=0,
                storm_end=100,
                name_prefix="9-CASCADE_AST_3domains_BE1m_AccSLR",
                storm_prefix="StormSeries_1kyrs_VCR_Berm1pt9m_Slope0pt04_",
            )

            one_hundred_thirds_acc_BE1m_ast_runs(
                storm_start=0,
                storm_end=100,
                name_prefix="9-CASCADE_AST_3domains_BE1m_AccSLR_AdaptationScenario_",
                storm_prefix="StormSeries_1kyrs_VCR_Berm1pt9m_Slope0pt04_FutureScenario",
            )


def cascade_1kyr_plots():
    def roadways():

        # rave = 0.75 runs, low
        def pt75_low():
            (
                BarrierWidth_nat,
                DuneCrestMean_nat,
                BarrierHeight_nat,
                bh_rate_nat,
                bw_rate_nat,
                sc_rate_nat,
                DuneCrestMin_nat,
                DuneCrestMax_nat,
                shoreline_position_nat,
                shoreface_slope_nat,
                overwash_nat,
                cascade_nat,
            ) = PLOT_5_Nonlinear_Dynamics_CASCADE_B3Donly_RoadwayManager(
                name_prefix="4-B3D_Rave_pt75_Natural_low",
                tmax_roadways=1000,  # dummy
                tmax_sim=1000,
                plot_name="b3d_pt75_plots_low",
                run_road_mgmt=False,
                gif_on=False,
                cross_sections=[28, 29, 30],
            )

            (
                BarrierWidth_h1m,
                DuneCrestMean_h1m,
                BarrierHeight_h1m,
                bh_rate_h1m,
                bw_rate_h1m,
                sc_rate_h1m,
                DuneCrestMin_h1m,
                DuneCrestMax_h1m,
                shoreline_position_h1m,
                shoreface_slope_h1m,
                overwash_h1m,
                cascade_h1m,
            ) = PLOT_5_Nonlinear_Dynamics_CASCADE_B3Donly_RoadwayManager(
                name_prefix="6-B3D_Rave_pt75_Roadways_1mDune_20mSetback_20mWidth_low",
                tmax_roadways=135,  # Barrier has HEIGHT DROWNED at t = 136 years
                tmax_sim=136,
                plot_name="b3d_pt75_h1m_plots_low",
                run_road_mgmt=True,
            )

            (
                BarrierWidth_h2m,
                DuneCrestMean_h2m,
                BarrierHeight_h2m,
                bh_rate_h2m,
                bw_rate_h2m,
                sc_rate_h2m,
                DuneCrestMin_h2m,
                DuneCrestMax_h2m,
                shoreline_position_h2m,
                shoreface_slope_h2m,
                overwash_h2m,
                cascade_h2m,
            ) = PLOT_5_Nonlinear_Dynamics_CASCADE_B3Donly_RoadwayManager(
                name_prefix="6-B3D_Rave_pt75_Roadways_2mDune_20mSetback_20mWidth_low",
                tmax_roadways=135,  # Barrier has HEIGHT DROWNED at t = 136 years
                tmax_sim=136,
                plot_name="b3d_pt75_h2m_plots_low",
                run_road_mgmt=True,
                gif_on=False,
                cross_sections=[0, 1, 45, 46],
            )

            (
                BarrierWidth_h3m,
                DuneCrestMean_h3m,
                BarrierHeight_h3m,
                bh_rate_h3m,
                bw_rate_h3m,
                sc_rate_h3m,
                DuneCrestMin_h3m,
                DuneCrestMax_h3m,
                shoreline_position_h3m,
                shoreface_slope_h3m,
                overwash_h3m,
                cascade_h3m,
            ) = PLOT_5_Nonlinear_Dynamics_CASCADE_B3Donly_RoadwayManager(
                name_prefix="6-B3D_Rave_pt75_Roadways_3mDune_20mSetback_20mWidth_low",
                tmax_roadways=131,  # Barrier has HEIGHT DROWNED at t = 132 years
                tmax_sim=132,
                plot_name="b3d_pt75_h3m_plots_low",
                run_road_mgmt=True,
            )

            CASCADEplt.plot_nonlinear_stats_mgmt_array4(
                cascade=[cascade_nat, cascade_h1m, cascade_h2m, cascade_h3m],
                DuneCrestMin=[
                    DuneCrestMin_nat,
                    DuneCrestMin_h1m,
                    DuneCrestMin_h2m,
                    DuneCrestMin_h3m,
                ],
                DuneCrestMax=[
                    DuneCrestMax_nat,
                    DuneCrestMax_h1m,
                    DuneCrestMax_h2m,
                    DuneCrestMax_h3m,
                ],
                BarrierHeight=[
                    BarrierHeight_nat,
                    BarrierHeight_h1m,
                    BarrierHeight_h2m,
                    BarrierHeight_h3m,
                ],
                BarrierWidth=[
                    BarrierWidth_nat,
                    BarrierWidth_h1m,
                    BarrierWidth_h2m,
                    BarrierWidth_h3m,
                ],
                DuneCrestMean=[
                    DuneCrestMean_nat,
                    DuneCrestMean_h1m,
                    DuneCrestMean_h2m,
                    DuneCrestMean_h3m,
                ],
                TMAX=[1000, 136, 136, 132],
                tmax_management=[
                    0,
                    135,
                    135,
                    131,
                ],
                shoreline_position=[
                    shoreline_position_nat,
                    shoreline_position_h1m,
                    shoreline_position_h2m,
                    shoreline_position_h3m,
                ],
                overwash=[
                    overwash_nat,
                    overwash_h1m,
                    overwash_h2m,
                    overwash_h3m,
                ],
                dune_toe=None,
                roadways_on=True,
                nourishment_on=False,
                rebuild_threshold=None,  # this comes from the roadways module
                scenarios=[
                    "natural",
                    "1 m",
                    "2 m",
                    "3 m",
                ],
            )

            # roadway statistics
            (
                year_abandoned,
                sim_max,
                road_bulldozed,
                overwash_removed,
                dune_rebuilt,
                road_relocated,
                diff_barrier_width,
                diff_barrier_elev,
            ) = get_roadway_statistics(
                folder_prefix="",
                natural_barrier_elev=BarrierHeight_nat[-1],
                natural_barrier_width=BarrierWidth_nat[-1],
                # individual_fid="6-B3D_Rave_pt75_Roadways_2mDune_20mSetback_20mWidth_low",
                individual_fid="6-B3D_Rave_pt75_Roadways_3mDune_20mSetback_20mWidth_low",
            )

        # rave = 0.75 runs, high
        def pt75_high():
            (
                BarrierWidth_nat,
                DuneCrestMean_nat,
                BarrierHeight_nat,
                bh_rate_nat,
                bw_rate_nat,
                sc_rate_nat,
                DuneCrestMin_nat,
                DuneCrestMax_nat,
                shoreline_position_nat,
                shoreface_slope_nat,
                overwash_nat,
                cascade_nat,
            ) = PLOT_5_Nonlinear_Dynamics_CASCADE_B3Donly_RoadwayManager(
                name_prefix="4-B3D_Rave_pt75_Natural_high",
                tmax_roadways=1000,
                tmax_sim=1000,
                plot_name="b3d_pt75_plots_high",
                run_road_mgmt=False,
            )

            (
                BarrierWidth_h1m,
                DuneCrestMean_h1m,
                BarrierHeight_h1m,
                bh_rate_h1m,
                bw_rate_h1m,
                sc_rate_h1m,
                DuneCrestMin_h1m,
                DuneCrestMax_h1m,
                shoreline_position_h1m,
                shoreface_slope_h1m,
                overwash_h1m,
                cascade_h1m,
            ) = PLOT_5_Nonlinear_Dynamics_CASCADE_B3Donly_RoadwayManager(
                name_prefix="6-B3D_Rave_pt75_Roadways_1mDune_20mSetback_20mWidth_high",
                tmax_roadways=535,  # Roadway width drowned at 535 years, 20.0% of road borders water
                tmax_sim=1000,
                plot_name="b3d_pt75_h1m_plots_high",
                run_road_mgmt=True,
            )

            (
                BarrierWidth_h2m,
                DuneCrestMean_h2m,
                BarrierHeight_h2m,
                bh_rate_h2m,
                bw_rate_h2m,
                sc_rate_h2m,
                DuneCrestMin_h2m,
                DuneCrestMax_h2m,
                shoreline_position_h2m,
                shoreface_slope_h2m,
                overwash_h2m,
                cascade_h2m,
            ) = PLOT_5_Nonlinear_Dynamics_CASCADE_B3Donly_RoadwayManager(
                name_prefix="6-B3D_Rave_pt75_Roadways_2mDune_20mSetback_20mWidth_high",
                tmax_roadways=520,  # Roadway width drowned at 520 years, 20.0% of road borders water
                tmax_sim=571,  # Barrier has HEIGHT DROWNED at t = 571 years
                plot_name="b3d_pt75_h2m_plots_high",
                run_road_mgmt=True,
            )

            (
                BarrierWidth_h3m,
                DuneCrestMean_h3m,
                BarrierHeight_h3m,
                bh_rate_h3m,
                bw_rate_h3m,
                sc_rate_h3m,
                DuneCrestMin_h3m,
                DuneCrestMax_h3m,
                shoreline_position_h3m,
                shoreface_slope_h3m,
                overwash_h3m,
                cascade_h3m,
            ) = PLOT_5_Nonlinear_Dynamics_CASCADE_B3Donly_RoadwayManager(
                name_prefix="6-B3D_Rave_pt75_Roadways_3mDune_20mSetback_20mWidth_high",
                tmax_roadways=395,  # Island is to narrow for roadway to be relocated. Roadway eaten up by dunes at 395 years
                tmax_sim=1000,
                plot_name="b3d_pt75_h3m_plots_high",
                run_road_mgmt=True,
            )

            CASCADEplt.plot_nonlinear_stats_mgmt_array4(
                cascade=[cascade_nat, cascade_h1m, cascade_h2m, cascade_h3m],
                DuneCrestMin=[
                    DuneCrestMin_nat,
                    DuneCrestMin_h1m,
                    DuneCrestMin_h2m,
                    DuneCrestMin_h3m,
                ],
                DuneCrestMax=[
                    DuneCrestMax_nat,
                    DuneCrestMax_h1m,
                    DuneCrestMax_h2m,
                    DuneCrestMax_h3m,
                ],
                BarrierHeight=[
                    BarrierHeight_nat,
                    BarrierHeight_h1m,
                    BarrierHeight_h2m,
                    BarrierHeight_h3m,
                ],
                BarrierWidth=[
                    BarrierWidth_nat,
                    BarrierWidth_h1m,
                    BarrierWidth_h2m,
                    BarrierWidth_h3m,
                ],
                DuneCrestMean=[
                    DuneCrestMean_nat,
                    DuneCrestMean_h1m,
                    DuneCrestMean_h2m,
                    DuneCrestMean_h3m,
                ],
                # TMAX=[
                #     750,
                #     750,
                #     582,
                #     750,
                # ],  # # h3, h1, h2 - 536, 530, 416 roadways drowned
                TMAX=[1000, 1000, 571, 1000],
                tmax_management=[0, 535, 520, 395],
                shoreline_position=[
                    shoreline_position_nat,
                    shoreline_position_h1m,
                    shoreline_position_h2m,
                    shoreline_position_h3m,
                ],
                overwash=[
                    overwash_nat,
                    overwash_h1m,
                    overwash_h2m,
                    overwash_h3m,
                ],
                dune_toe=None,
                roadways_on=True,
                nourishment_on=False,
                rebuild_threshold=None,  # this comes from the roadways module
                scenarios=[
                    "natural",
                    "1 m",
                    "2 m",
                    "3 m",
                ],
            )

        # rave = 0.45 runs, low
        def pt45_low():
            (
                BarrierWidth_nat,
                DuneCrestMean_nat,
                BarrierHeight_nat,
                bh_rate_nat,
                bw_rate_nat,
                sc_rate_nat,
                DuneCrestMin_nat,
                DuneCrestMax_nat,
                shoreline_position_nat,
                shoreface_slope_nat,
                overwash_nat,
                cascade_nat,
            ) = PLOT_5_Nonlinear_Dynamics_CASCADE_B3Donly_RoadwayManager(
                name_prefix="4-B3D_Rave_pt45_Natural_low",
                tmax_roadways=1000,
                tmax_sim=1000,
                plot_name="b3d_pt45_plots_low",
                run_road_mgmt=False,
                gif_on=False,
                cross_sections=[0, 1, 2, 3, 4],
            )

            (
                BarrierWidth_h1m,
                DuneCrestMean_h1m,
                BarrierHeight_h1m,
                bh_rate_h1m,
                bw_rate_h1m,
                sc_rate_h1m,
                DuneCrestMin_h1m,
                DuneCrestMax_h1m,
                shoreline_position_h1m,
                shoreface_slope_h1m,
                overwash_h1m,
                cascade_h1m,
            ) = PLOT_5_Nonlinear_Dynamics_CASCADE_B3Donly_RoadwayManager(
                name_prefix="6-B3D_Rave_pt45_Roadways_1mDune_20mSetback_20mWidth_low",
                tmax_roadways=544,  # Roadway width drowned at 544 years, 20.0% of road borders water
                tmax_sim=1000,
                plot_name="b3d_pt45_h1m_plots_low",
                run_road_mgmt=True,
            )

            (
                BarrierWidth_h2m,
                DuneCrestMean_h2m,
                BarrierHeight_h2m,
                bh_rate_h2m,
                bw_rate_h2m,
                sc_rate_h2m,
                DuneCrestMin_h2m,
                DuneCrestMax_h2m,
                shoreline_position_h2m,
                shoreface_slope_h2m,
                overwash_h2m,
                cascade_h2m,
            ) = PLOT_5_Nonlinear_Dynamics_CASCADE_B3Donly_RoadwayManager(
                name_prefix="6-B3D_Rave_pt45_Roadways_2mDune_20mSetback_20mWidth_low",
                tmax_roadways=533,  # Roadway width drowned at 533 years, 20.0% of road borders water
                tmax_sim=1000,
                plot_name="b3d_pt45_h2m_plots_low",
                run_road_mgmt=True,
            )

            (
                BarrierWidth_h3m,
                DuneCrestMean_h3m,
                BarrierHeight_h3m,
                bh_rate_h3m,
                bw_rate_h3m,
                sc_rate_h3m,
                DuneCrestMin_h3m,
                DuneCrestMax_h3m,
                shoreline_position_h3m,
                shoreface_slope_h3m,
                overwash_h3m,
                cascade_h3m,
            ) = PLOT_5_Nonlinear_Dynamics_CASCADE_B3Donly_RoadwayManager(
                name_prefix="6-B3D_Rave_pt45_Roadways_3mDune_20mSetback_20mWidth_low",
                tmax_roadways=322,  # Roadway width drowned at 322 years, 20.0% of road borders water
                tmax_sim=1000,
                plot_name="b3d_pt45_h3m_plots_low",
                run_road_mgmt=True,
            )

            CASCADEplt.plot_nonlinear_stats_mgmt_array4(
                cascade=[cascade_nat, cascade_h1m, cascade_h2m, cascade_h3m],
                DuneCrestMin=[
                    DuneCrestMin_nat,
                    DuneCrestMin_h1m,
                    DuneCrestMin_h2m,
                    DuneCrestMin_h3m,
                ],
                DuneCrestMax=[
                    DuneCrestMax_nat,
                    DuneCrestMax_h1m,
                    DuneCrestMax_h2m,
                    DuneCrestMax_h3m,
                ],
                BarrierHeight=[
                    BarrierHeight_nat,
                    BarrierHeight_h1m,
                    BarrierHeight_h2m,
                    BarrierHeight_h3m,
                ],
                BarrierWidth=[
                    BarrierWidth_nat,
                    BarrierWidth_h1m,
                    BarrierWidth_h2m,
                    BarrierWidth_h3m,
                ],
                DuneCrestMean=[
                    DuneCrestMean_nat,
                    DuneCrestMean_h1m,
                    DuneCrestMean_h2m,
                    DuneCrestMean_h3m,
                ],
                TMAX=[1000, 1000, 1000, 1000],
                tmax_management=[0, 544, 533, 322],
                shoreline_position=[
                    shoreline_position_nat,
                    shoreline_position_h1m,
                    shoreline_position_h2m,
                    shoreline_position_h3m,
                ],
                overwash=[
                    overwash_nat,
                    overwash_h1m,
                    overwash_h2m,
                    overwash_h3m,
                ],
                dune_toe=None,
                roadways_on=True,
                nourishment_on=False,
                rebuild_threshold=None,  # this comes from the roadways module
                scenarios=[
                    "natural",
                    "1 m",
                    "2 m",
                    "3 m",
                ],
            )

        # rave = 0.45 runs, high
        def pt45_high():
            (
                BarrierWidth_nat,
                DuneCrestMean_nat,
                BarrierHeight_nat,
                bh_rate_nat,
                bw_rate_nat,
                sc_rate_nat,
                DuneCrestMin_nat,
                DuneCrestMax_nat,
                shoreline_position_nat,
                shoreface_slope_nat,
                overwash_nat,
                cascade_nat,
            ) = PLOT_5_Nonlinear_Dynamics_CASCADE_B3Donly_RoadwayManager(
                name_prefix="4-B3D_Rave_pt45_Natural_high",
                tmax_roadways=1000,
                tmax_sim=1000,
                plot_name="b3d_pt45_plots_high",
                run_road_mgmt=False,
                gif_on=False,
                cross_sections=[0, 1, 2, 3, 4],
            )

            (
                BarrierWidth_h1m,
                DuneCrestMean_h1m,
                BarrierHeight_h1m,
                bh_rate_h1m,
                bw_rate_h1m,
                sc_rate_h1m,
                DuneCrestMin_h1m,
                DuneCrestMax_h1m,
                shoreline_position_h1m,
                shoreface_slope_h1m,
                overwash_h1m,
                cascade_h1m,
            ) = PLOT_5_Nonlinear_Dynamics_CASCADE_B3Donly_RoadwayManager(
                name_prefix="6-B3D_Rave_pt45_Roadways_1mDune_20mSetback_20mWidth_high",
                tmax_roadways=650,  # Island is too narrow for roadway to be relocated. Roadway eaten up by dunes at 650 years
                tmax_sim=1000,
                plot_name="b3d_pt45_h1m_plots_high",
                run_road_mgmt=True,
            )

            (
                BarrierWidth_h2m,
                DuneCrestMean_h2m,
                BarrierHeight_h2m,
                bh_rate_h2m,
                bw_rate_h2m,
                sc_rate_h2m,
                DuneCrestMin_h2m,
                DuneCrestMax_h2m,
                shoreline_position_h2m,
                shoreface_slope_h2m,
                overwash_h2m,
                cascade_h2m,
            ) = PLOT_5_Nonlinear_Dynamics_CASCADE_B3Donly_RoadwayManager(
                name_prefix="6-B3D_Rave_pt45_Roadways_2mDune_20mSetback_20mWidth_high",
                tmax_roadways=628,  # Roadway width drowned at 628 years, 20.0% of road borders water
                tmax_sim=1000,
                plot_name="b3d_pt45_h2m_plots_high",
                run_road_mgmt=True,
            )

            # (
            #     BarrierWidth_h3m,
            #     DuneCrestMean_h3m,
            #     BarrierHeight_h3m,
            #     bh_rate_h3m,
            #     bw_rate_h3m,
            #     sc_rate_h3m,
            #     DuneCrestMin_h3m,
            #     DuneCrestMax_h3m,
            #     shoreline_position_h3m,
            #     shoreface_slope_h3m,
            #     overwash_h3m,
            #     cascade_h3m,
            # ) = PLOT_5_Nonlinear_Dynamics_CASCADE_B3Donly_RoadwayManager(
            #     name_prefix="6-B3D_Rave_pt45_Roadways_3mDune_20mSetback_20mWidth_high",
            #     tmin=0,
            #     tmax_roadways=521,
            #     tmax_sim=532,
            #     plot_name="b3d_pt45_h3m_plots_high",
            #     run_road_mgmt=True,
            # )

            (
                BarrierWidth_h3m,
                DuneCrestMean_h3m,
                BarrierHeight_h3m,
                bh_rate_h3m,
                bw_rate_h3m,
                sc_rate_h3m,
                DuneCrestMin_h3m,
                DuneCrestMax_h3m,
                shoreline_position_h3m,
                shoreface_slope_h3m,
                overwash_h3m,
                cascade_h3m,
            ) = PLOT_5_Nonlinear_Dynamics_CASCADE_B3Donly_RoadwayManager(
                name_prefix="6-B3D_Rave_pt45_Roadways_3mDune_20mSetback_20mWidth_high",
                tmax_roadways=522,  # Roadway width drowned at 522 years, 20.0% of road borders water
                tmax_sim=1000,
                plot_name="b3d_pt45_h3m_plots_high",
                run_road_mgmt=True,
            )

            CASCADEplt.plot_nonlinear_stats_mgmt_array4(
                cascade=[cascade_nat, cascade_h1m, cascade_h2m, cascade_h3m],
                DuneCrestMin=[
                    DuneCrestMin_nat,
                    DuneCrestMin_h1m,
                    DuneCrestMin_h2m,
                    DuneCrestMin_h3m,
                ],
                DuneCrestMax=[
                    DuneCrestMax_nat,
                    DuneCrestMax_h1m,
                    DuneCrestMax_h2m,
                    DuneCrestMax_h3m,
                ],
                BarrierHeight=[
                    BarrierHeight_nat,
                    BarrierHeight_h1m,
                    BarrierHeight_h2m,
                    BarrierHeight_h3m,
                ],
                BarrierWidth=[
                    BarrierWidth_nat,
                    BarrierWidth_h1m,
                    BarrierWidth_h2m,
                    BarrierWidth_h3m,
                ],
                DuneCrestMean=[
                    DuneCrestMean_nat,
                    DuneCrestMean_h1m,
                    DuneCrestMean_h2m,
                    DuneCrestMean_h3m,
                ],
                TMAX=[1000, 1000, 1000, 1000],
                tmax_management=[0, 650, 628, 522],
                shoreline_position=[
                    shoreline_position_nat,
                    shoreline_position_h1m,
                    shoreline_position_h2m,
                    shoreline_position_h3m,
                ],
                overwash=[
                    overwash_nat,
                    overwash_h1m,
                    overwash_h2m,
                    overwash_h3m,
                ],
                dune_toe=None,
                roadways_on=True,
                nourishment_on=False,
                rebuild_threshold=None,  # this comes from the roadways module
                scenarios=[
                    "natural",
                    "1 m",
                    "2 m",
                    "3 m",
                ],
            )

            # roadway statistics
            (
                year_abandoned,
                sim_max,
                road_bulldozed,
                overwash_removed,
                dune_rebuilt,
                road_relocated,
                diff_barrier_width,
                diff_barrier_elev,
            ) = get_roadway_statistics(
                folder_prefix="",
                # natural_barrier_elev=BarrierHeight_nat[-1],
                # natural_barrier_width=BarrierWidth_nat[-1],
                # individual_fid="6-B3D_Rave_pt45_Roadways_3mDune_20mSetback_20mWidth_high",
                # individual_fid="6-B3D_Rave_pt45_Roadways_2mDune_20mSetback_20mWidth_high",
                individual_fid="6-B3D_Rave_pt45_Roadways_1mDune_20mSetback_20mWidth_high",
                tmax=100,
            )

        # supplementary material
        def sensitivity_abandonment_criteria():
            (
                BarrierWidth_nat,
                DuneCrestMean_nat,
                BarrierHeight_nat,
                bh_rate_nat,
                bw_rate_nat,
                sc_rate_nat,
                DuneCrestMin_nat,
                DuneCrestMax_nat,
                shoreline_position_nat,
                shoreface_slope_nat,
                overwash_nat,
                cascade_nat,
            ) = PLOT_5_Nonlinear_Dynamics_CASCADE_B3Donly_RoadwayManager(
                name_prefix="6-B3D_Rave_pt75_Natural_low",
                tmax_roadways=500,  # dummy
                tmax_sim=500,
                plot_name="b3d_pt75_plots_low",
                run_road_mgmt=False,
            )

            (
                BarrierWidth_10,
                DuneCrestMean_10,
                BarrierHeight_10,
                bh_rate_10,
                bw_rate_10,
                sc_rate_10,
                DuneCrestMin_10,
                DuneCrestMax_10,
                shoreline_position_10,
                shoreface_slope_10,
                overwash_10,
                cascade_10,
            ) = PLOT_5_Nonlinear_Dynamics_CASCADE_B3Donly_RoadwayManager(
                name_prefix="6-B3D_Rave_pt75_Roadways_2mDune_20mSetback_20mWidth_low_10percent",
                tmax_roadways=158,
                tmax_sim=500,
                plot_name="b3d_pt75_h2m_plots_low_10percent",
                run_road_mgmt=True,
            )

            (
                BarrierWidth_20,
                DuneCrestMean_20,
                BarrierHeight_20,
                bh_rate_20,
                bw_rate_20,
                sc_rate_20,
                DuneCrestMin_20,
                DuneCrestMax_20,
                shoreline_position_20,
                shoreface_slope_20,
                overwash_20,
                cascade_20,
            ) = PLOT_5_Nonlinear_Dynamics_CASCADE_B3Donly_RoadwayManager(
                name_prefix="6-B3D_Rave_pt75_Roadways_2mDune_20mSetback_20mWidth_low_20percent",
                tmax_roadways=161,  # drowned at 162
                tmax_sim=500,
                plot_name="b3d_pt75_h2m_plots_low_20percent",
                run_road_mgmt=True,
            )

            (
                BarrierWidth_30,
                DuneCrestMean_30,
                BarrierHeight_30,
                bh_rate_30,
                bw_rate_30,
                sc_rate_30,
                DuneCrestMin_30,
                DuneCrestMax_30,
                shoreline_position_30,
                shoreface_slope_30,
                overwash_30,
                cascade_30,
            ) = PLOT_5_Nonlinear_Dynamics_CASCADE_B3Donly_RoadwayManager(
                name_prefix="6-B3D_Rave_pt75_Roadways_2mDune_20mSetback_20mWidth_low_30percent",
                tmax_roadways=164,  # drowned at 165
                tmax_sim=500,
                plot_name="b3d_pt75_h2m_plots_low_30percent",
                run_road_mgmt=True,
            )

            (
                BarrierWidth_40,
                DuneCrestMean_40,
                BarrierHeight_40,
                bh_rate_40,
                bw_rate_40,
                sc_rate_40,
                DuneCrestMin_40,
                DuneCrestMax_40,
                shoreline_position_40,
                shoreface_slope_40,
                overwash_40,
                cascade_40,
            ) = PLOT_5_Nonlinear_Dynamics_CASCADE_B3Donly_RoadwayManager(
                name_prefix="6-B3D_Rave_pt75_Roadways_2mDune_20mSetback_20mWidth_low_40percent",
                tmax_roadways=172,  # drowned at 173
                tmax_sim=500,
                plot_name="b3d_pt75_h2m_plots_low_40percent",
                run_road_mgmt=True,
            )

            (
                BarrierWidth_50,
                DuneCrestMean_50,
                BarrierHeight_50,
                bh_rate_50,
                bw_rate_50,
                sc_rate_50,
                DuneCrestMin_50,
                DuneCrestMax_50,
                shoreline_position_50,
                shoreface_slope_50,
                overwash_50,
                cascade_50,
            ) = PLOT_5_Nonlinear_Dynamics_CASCADE_B3Donly_RoadwayManager(
                name_prefix="6-B3D_Rave_pt75_Roadways_2mDune_20mSetback_20mWidth_low_50percent",
                tmax_roadways=181,  # drowned at 182
                tmax_sim=500,
                plot_name="b3d_pt75_h2m_plots_low_50percent",
                run_road_mgmt=True,
            )

            CASCADEplt.plot_nonlinear_stats_low_high_sensitivity(
                cascade=[
                    cascade_10,
                    cascade_20,
                    cascade_30,
                    cascade_40,
                    cascade_50,
                ],
                BarrierHeight=[
                    BarrierHeight_10,
                    BarrierHeight_20,
                    BarrierHeight_30,
                    BarrierHeight_40,
                    BarrierHeight_50,
                ],
                BarrierWidth=[
                    BarrierWidth_10,
                    BarrierWidth_20,
                    BarrierWidth_30,
                    BarrierWidth_40,
                    BarrierWidth_50,
                ],
                TMAX=[500, 500, 500, 500, 500],
                tmax_roadways=[158, 161, 164, 172, 181],
            )

        # roadway statistics
        (
            year_abandoned,
            sim_max,
            drown,
            road_bulldozed,
            overwash_removed,
            dune_rebuilt,
            road_relocated,
            diff_barrier_height,
            diff_barrier_elev,
        ) = get_roadway_statistics(
            folder_prefix="",
            natural_barrier_elev=0.72,
            natural_barrier_width=229,
            individual_fid="6-B3D_Rave_pt75_Roadways_2mDune_20mSetback_20mWidth_low",
        )

    def slr_sensitivity():
        def pt45_high():
            (
                BarrierWidth_pt45_high_SLRacc,  # we only use this one
                DuneCrestMean,
                BarrierHeight,
                bh_rate,
                bw_rate,
                sc_rate,
                DuneCrestMin,
                DuneCrestMax,
                shoreline_position,
                shoreface_slope,
                overwash,
                cascade_pt45_high_SLRacc,  # and this one
            ) = PLOT_5_Nonlinear_Dynamics_CASCADE_B3Donly_RoadwayManager(
                name_prefix="7-B3D_Rave_pt45_Natural_high_AccSLR",
                tmax_sim=200,
                tmax_roadways=200,  # dummy
                plot_name="b3d_pt45_Natural_high_AccSLR",
                run_road_mgmt=False,
            )

            (
                BarrierWidth_pt45_high_0pt012SLR,
                DuneCrestMean,
                BarrierHeight,
                bh_rate,
                bw_rate,
                sc_rate,
                DuneCrestMin,
                DuneCrestMax,
                shoreline_position,
                shoreface_slope,
                overwash,
                cascade_pt45_high_0pt012SLR,
            ) = PLOT_5_Nonlinear_Dynamics_CASCADE_B3Donly_RoadwayManager(
                name_prefix="7-B3D_Rave_pt45_Natural_high_0pt012SLR",
                tmax_sim=200,
                tmax_roadways=200,  # dummy
                plot_name="b3d_pt45_Natural_high_0pt012SLR",
                run_road_mgmt=False,
            )

            (
                BarrierWidth_pt45_high_0pt008SLR,
                DuneCrestMean,
                BarrierHeight,
                bh_rate,
                bw_rate,
                sc_rate,
                DuneCrestMin,
                DuneCrestMax,
                shoreline_position,
                shoreface_slope,
                overwash,
                cascade_pt45_high_0pt008SLR,
            ) = PLOT_5_Nonlinear_Dynamics_CASCADE_B3Donly_RoadwayManager(
                name_prefix="7-B3D_Rave_pt45_Natural_high_0pt008SLR",
                tmax_sim=200,
                tmax_roadways=200,  # dummy
                plot_name="b3d_pt45_Natural_high_0pt008SLR",
                run_road_mgmt=False,
            )

            (
                BarrierWidth_pt45_high_0pt004SLR,
                DuneCrestMean,
                BarrierHeight,
                bh_rate,
                bw_rate,
                sc_rate,
                DuneCrestMin,
                DuneCrestMax,
                shoreline_position,
                shoreface_slope,
                overwash,
                cascade_pt45_high_0pt004SLR,
            ) = PLOT_5_Nonlinear_Dynamics_CASCADE_B3Donly_RoadwayManager(
                name_prefix="4-B3D_Rave_pt45_Natural_high",  # this is the 0.004 case
                tmax_sim=200,
                tmax_roadways=200,  # dummy
                plot_name="b3d_pt45_Natural_high_0pt004SLR",
                run_road_mgmt=False,
            )

            cascade = [
                cascade_pt45_high_0pt004SLR,
                cascade_pt45_high_0pt008SLR,
                cascade_pt45_high_0pt012SLR,
                cascade_pt45_high_SLRacc,
            ]
            TMAX = [200, 200, 200, 200]
            CASCADEplt.fig3_slr_sensitivity(
                cascade,  # lists
                TMAX,
            )

        def pt45_low():
            (
                BarrierWidth_pt45_low_SLRacc,
                DuneCrestMean,
                BarrierHeight,
                bh_rate,
                bw_rate,
                sc_rate,
                DuneCrestMin,
                DuneCrestMax,
                shoreline_position,
                shoreface_slope,
                overwash,
                cascade_pt45_low_SLRacc,
            ) = PLOT_5_Nonlinear_Dynamics_CASCADE_B3Donly_RoadwayManager(
                name_prefix="7-B3D_Rave_pt45_Natural_low_AccSLR",
                tmax_sim=200,
                tmax_roadways=200,  # dummy
                plot_name="b3d_pt45_Natural_low_AccSLR",
                run_road_mgmt=False,
            )

            (
                BarrierWidth_pt45_low_0pt012SLR,
                DuneCrestMean,
                BarrierHeight,
                bh_rate,
                bw_rate,
                sc_rate,
                DuneCrestMin,
                DuneCrestMax,
                shoreline_position,
                shoreface_slope,
                overwash,
                cascade_pt45_low_0pt012SLR,
            ) = PLOT_5_Nonlinear_Dynamics_CASCADE_B3Donly_RoadwayManager(
                name_prefix="7-B3D_Rave_pt45_Natural_low_0pt012SLR",
                tmax_sim=200,
                tmax_roadways=200,  # dummy
                plot_name="b3d_pt45_Natural_low_0pt012SLR",
                run_road_mgmt=False,
            )

            (
                BarrierWidth_pt45_low_0pt008SLR,
                DuneCrestMean,
                BarrierHeight,
                bh_rate,
                bw_rate,
                sc_rate,
                DuneCrestMin,
                DuneCrestMax,
                shoreline_position,
                shoreface_slope,
                overwash,
                cascade_pt45_low_0pt008SLR,
            ) = PLOT_5_Nonlinear_Dynamics_CASCADE_B3Donly_RoadwayManager(
                name_prefix="7-B3D_Rave_pt45_Natural_low_0pt008SLR",
                tmax_sim=200,
                tmax_roadways=200,  # dummy
                plot_name="b3d_pt45_Natural_low_0pt008SLR",
                run_road_mgmt=False,
            )

            (
                BarrierWidth_pt45_low_0pt004SLR,
                DuneCrestMean,
                BarrierHeight,
                bh_rate,
                bw_rate,
                sc_rate,
                DuneCrestMin,
                DuneCrestMax,
                shoreline_position,
                shoreface_slope,
                overwash,
                cascade_pt45_low_0pt004SLR,
            ) = PLOT_5_Nonlinear_Dynamics_CASCADE_B3Donly_RoadwayManager(
                name_prefix="4-B3D_Rave_pt45_Natural_low",  # this is the 0.004 case
                tmax_sim=200,
                tmax_roadways=200,  # dummy
                plot_name="b3d_pt45_Natural_low_0pt004SLR",
                run_road_mgmt=False,
            )

            cascade = [
                cascade_pt45_low_0pt004SLR,
                cascade_pt45_low_0pt008SLR,
                cascade_pt45_low_0pt012SLR,
                cascade_pt45_low_SLRacc,
            ]
            TMAX = [200, 200, 200, 200]
            CASCADEplt.fig3_slr_sensitivity(
                cascade,  # lists
                TMAX,
            )

        def pt75_low():
            (
                BarrierWidth_pt75_low_SLRacc,
                DuneCrestMean,
                BarrierHeight,
                bh_rate,
                bw_rate,
                sc_rate,
                DuneCrestMin,
                DuneCrestMax,
                shoreline_position,
                shoreface_slope,
                overwash,
                cascade_pt75_low_SLRacc,
            ) = PLOT_5_Nonlinear_Dynamics_CASCADE_B3Donly_RoadwayManager(
                name_prefix="7-B3D_Rave_pt75_Natural_low_AccSLR",
                tmax_sim=200,
                tmax_roadways=200,  # dummy
                plot_name="b3d_pt75_Natural_low_AccSLR",
                run_road_mgmt=False,
            )

            (
                BarrierWidth_pt75_low_0pt012SLR,
                DuneCrestMean,
                BarrierHeight,
                bh_rate,
                bw_rate,
                sc_rate,
                DuneCrestMin,
                DuneCrestMax,
                shoreline_position,
                shoreface_slope,
                overwash,
                cascade_pt75_low_0pt012SLR,
            ) = PLOT_5_Nonlinear_Dynamics_CASCADE_B3Donly_RoadwayManager(
                name_prefix="7-B3D_Rave_pt75_Natural_low_0pt012SLR",
                tmax_sim=200,
                tmax_roadways=200,  # dummy
                plot_name="b3d_pt75_Natural_low_0pt012SLR",
                run_road_mgmt=False,
            )

            (
                BarrierWidth_pt75_low_0pt008SLR,
                DuneCrestMean,
                BarrierHeight,
                bh_rate,
                bw_rate,
                sc_rate,
                DuneCrestMin,
                DuneCrestMax,
                shoreline_position,
                shoreface_slope,
                overwash,
                cascade_pt75_low_0pt008SLR,
            ) = PLOT_5_Nonlinear_Dynamics_CASCADE_B3Donly_RoadwayManager(
                name_prefix="7-B3D_Rave_pt75_Natural_low_0pt008SLR",
                tmax_sim=200,
                tmax_roadways=200,  # dummy
                plot_name="b3d_pt75_Natural_low_0pt008SLR",
                run_road_mgmt=False,
            )

            (
                BarrierWidth_pt75_low_0pt004SLR,
                DuneCrestMean,
                BarrierHeight,
                bh_rate,
                bw_rate,
                sc_rate,
                DuneCrestMin,
                DuneCrestMax,
                shoreline_position,
                shoreface_slope,
                overwash,
                cascade_pt75_low_0pt004SLR,
            ) = PLOT_5_Nonlinear_Dynamics_CASCADE_B3Donly_RoadwayManager(
                name_prefix="4-B3D_Rave_pt75_Natural_low",  # this is the 0.004 case
                tmax_sim=200,
                tmax_roadways=200,  # dummy
                plot_name="b3d_pt75_Natural_low_0pt004SLR",
                run_road_mgmt=False,
            )

            cascade = [
                cascade_pt75_low_0pt004SLR,
                cascade_pt75_low_0pt008SLR,
                cascade_pt75_low_0pt012SLR,
                cascade_pt75_low_SLRacc,
            ]
            TMAX = [200, 200, 200, 200]
            CASCADEplt.fig3_slr_sensitivity(
                cascade,  # lists
                TMAX,
            )

        def pt75_high():
            (
                BarrierWidth_pt75_high_SLRacc,
                DuneCrestMean,
                BarrierHeight,
                bh_rate,
                bw_rate,
                sc_rate,
                DuneCrestMin,
                DuneCrestMax,
                shoreline_position,
                shoreface_slope,
                overwash,
                cascade_pt75_high_SLRacc,
            ) = PLOT_5_Nonlinear_Dynamics_CASCADE_B3Donly_RoadwayManager(
                name_prefix="7-B3D_Rave_pt75_Natural_high_AccSLR",
                tmax_sim=200,
                tmax_roadways=200,  # dummy
                plot_name="b3d_pt75_Natural_high_AccSLR",
                run_road_mgmt=False,
            )

            (
                BarrierWidth_pt75_high_0pt012SLR,
                DuneCrestMean,
                BarrierHeight,
                bh_rate,
                bw_rate,
                sc_rate,
                DuneCrestMin,
                DuneCrestMax,
                shoreline_position,
                shoreface_slope,
                overwash,
                cascade_pt75_high_0pt012SLR,
            ) = PLOT_5_Nonlinear_Dynamics_CASCADE_B3Donly_RoadwayManager(
                name_prefix="7-B3D_Rave_pt75_Natural_high_0pt012SLR",
                tmax_sim=200,
                tmax_roadways=200,  # dummy
                plot_name="b3d_pt75_Natural_high_0pt012SLR",
                run_road_mgmt=False,
            )

            (
                BarrierWidth_pt75_high_0pt008SLR,
                DuneCrestMean,
                BarrierHeight,
                bh_rate,
                bw_rate,
                sc_rate,
                DuneCrestMin,
                DuneCrestMax,
                shoreline_position,
                shoreface_slope,
                overwash,
                cascade_pt75_high_0pt008SLR,
            ) = PLOT_5_Nonlinear_Dynamics_CASCADE_B3Donly_RoadwayManager(
                name_prefix="7-B3D_Rave_pt75_Natural_high_0pt008SLR",
                tmax_sim=200,
                tmax_roadways=200,  # dummy
                plot_name="b3d_pt75_Natural_high_0pt008SLR",
                run_road_mgmt=False,
            )

            (
                BarrierWidth_pt75_high_0pt004SLR,
                DuneCrestMean,
                BarrierHeight,
                bh_rate,
                bw_rate,
                sc_rate,
                DuneCrestMin,
                DuneCrestMax,
                shoreline_position,
                shoreface_slope,
                overwash,
                cascade_pt75_high_0pt004SLR,
            ) = PLOT_5_Nonlinear_Dynamics_CASCADE_B3Donly_RoadwayManager(
                name_prefix="4-B3D_Rave_pt75_Natural_high",  # this is the 0.004 case
                tmax_sim=200,
                tmax_roadways=200,  # dummy
                plot_name="b3d_pt75_Natural_high_0pt004SLR",
                run_road_mgmt=False,
            )

            cascade = [
                cascade_pt75_high_0pt004SLR,
                cascade_pt75_high_0pt008SLR,
                cascade_pt75_high_0pt012SLR,
                cascade_pt75_high_SLRacc,
            ]
            TMAX = [200, 200, 200, 200]
            CASCADEplt.fig3_slr_sensitivity(
                cascade,  # lists
                TMAX,
            )

    def nourishments():
        def rebuild_threshold_1m():
            def pt45_low():
                (
                    BarrierWidth_pt45_nat,
                    DuneCrestMean_pt45_nat,
                    BarrierHeight_pt45_nat,
                    bh_rate_pt45_nat,
                    bw_rate_pt45_nat,
                    sc_rate_pt45_nat,
                    DuneCrestMin_pt45_nat,
                    DuneCrestMax_pt45_nat,
                    shoreline_position_pt45_nat,
                    shoreface_slope_pt45_nat,
                    overwash_pt45_nat,
                    cascade_pt45_nat,
                ) = PLOT_5_Nonlinear_Dynamics_CASCADE_B3Donly_RoadwayManager(
                    name_prefix="4-B3D_Rave_pt45_Natural_low",
                    tmax_roadways=1000,
                    tmax_sim=1000,
                    plot_name="b3d_pt45_plots_low",
                    run_road_mgmt=False,
                    gif_on=False,
                )

                (
                    BarrierWidth_pt45_low_40pc,
                    DuneCrestMean_pt45_low_40pc,
                    BarrierHeight_pt45_low_40pc,
                    bh_rate_pt45_low_40pc,
                    bw_rate_pt45_low_40pc,
                    sc_rate_pt45_low_40pc,
                    DuneCrestMin_pt45_low_40pc,
                    DuneCrestMax_pt45_low_40pc,
                    shoreline_position_pt45_low_40pc,
                    shoreface_slope_pt45_low_40pc,
                    beach_width_pt45_low_40pc,
                    overwash_pt45_low_40pc,
                    dune_toe_pt45_low_40pc,
                    cascade_pt45_low_40pc,
                ) = PLOT_6_Nonlinear_Dynamics_CASCADE_B3Donly_Nourishments(
                    name_prefix="8-B3D_Rave_pt45_Nourishment_2mDune_lowEle_residential_rebuild1m",
                    tmax_management=407,  # Community reached minimum width, drowned at 407 years
                    tmax_sim=1000,
                    plot_name="b3d_pt45_Nourishment_2mDune_lowEle_residential_rebuild1m",
                    rebuild_dune_threshold=1,
                )

                (
                    BarrierWidth_pt45_low_90pc,
                    DuneCrestMean_pt45_low_90pc,
                    BarrierHeight_pt45_low_90pc,
                    bh_rate_pt45_low_90pc,
                    bw_rate_pt45_low_90pc,
                    sc_rate_pt45_low_90pc,
                    DuneCrestMin_pt45_low_90pc,
                    DuneCrestMax_pt45_low_90pc,
                    shoreline_position_pt45_low_90pc,
                    shoreface_slope_pt45_low_90pc,
                    beach_width_pt45_low_90pc,
                    overwash_pt45_low_90pc,
                    dune_toe_pt45_low_90pc,
                    cascade_pt45_low_90pc,
                ) = PLOT_6_Nonlinear_Dynamics_CASCADE_B3Donly_Nourishments(
                    name_prefix="8-B3D_Rave_pt45_Nourishment_2mDune_lowEle_commercial_rebuild1m",
                    tmax_management=302,  # Community reached minimum width, drowned at 302 years
                    tmax_sim=1000,
                    plot_name="b3d_pt45_Nourishment_2mDune_lowEle_commercial_rebuild1m",
                    rebuild_dune_threshold=1,
                )

                (
                    BarrierWidth_pt45_low_90pc_backerosion_1m,
                    DuneCrestMean_pt45_low_90pc_backerosion_1m,
                    BarrierHeight_pt45_low_90pc_backerosion_1m,
                    bh_rate_pt45_low_90pc_backerosion_1m,
                    bw_rate_pt45_low_90pc_backerosion_1m,
                    sc_rate_pt45_low_90pc_backerosion_1m,
                    DuneCrestMin_pt45_low_90pc_backerosion_1m,
                    DuneCrestMax_pt45_low_90pc_backerosion_1m,
                    shoreline_position_pt45_low_90pc_backerosion_1m,
                    shoreface_slope_pt45_low_90pc_backerosion_1m,
                    beach_width_pt45_low_90pc_backerosion_1m,
                    overwash_pt45_low_90pc_backerosion_1m,
                    dune_toe_pt45_low_90pc_backerosion_1m,
                    cascade_pt45_low_90pc_backerosion_1m,
                ) = PLOT_6_Nonlinear_Dynamics_CASCADE_B3Donly_Nourishments(
                    name_prefix="8-B3D_Rave_pt45_Nourishment_2mDune_lowEle_commercial_backerosion1m_rebuild1m",
                    tmax_management=302,  # Community reached minimum width, drowned at 302 years
                    tmax_sim=1000,
                    plot_name="b3d_pt45_Nourishment_2mDune_lowEle_commercial_backerosion_1m_rebuild1m",
                    rebuild_dune_threshold=1,
                )

                rebuild_threshold = 1 + (cascade_pt45_low_40pc.barrier3d[0].BermEl * 10)

                def version1_residential():
                    CASCADEplt.plot_nonlinear_stats_mgmt_array4(
                        cascade=[
                            cascade_pt45_nat,
                            cascade_pt45_low_40pc,  # only residential here
                            cascade_pt45_low_90pc,
                            cascade_pt45_low_90pc_backerosion_1m,
                        ],
                        DuneCrestMin=[
                            DuneCrestMin_pt45_nat,
                            DuneCrestMin_pt45_low_40pc,
                            DuneCrestMin_pt45_low_90pc,
                            DuneCrestMin_pt45_low_90pc_backerosion_1m,
                        ],
                        DuneCrestMax=[
                            DuneCrestMax_pt45_nat,
                            DuneCrestMax_pt45_low_40pc,
                            DuneCrestMax_pt45_low_90pc,
                            DuneCrestMax_pt45_low_90pc_backerosion_1m,
                        ],
                        BarrierHeight=[
                            BarrierHeight_pt45_nat,
                            BarrierHeight_pt45_low_40pc,
                            BarrierHeight_pt45_low_90pc,
                            BarrierHeight_pt45_low_90pc_backerosion_1m,
                        ],
                        BarrierWidth=[
                            BarrierWidth_pt45_nat,
                            BarrierWidth_pt45_low_40pc,
                            BarrierWidth_pt45_low_90pc,
                            BarrierWidth_pt45_low_90pc_backerosion_1m,
                        ],
                        DuneCrestMean=[
                            DuneCrestMean_pt45_nat,
                            DuneCrestMean_pt45_low_40pc,
                            DuneCrestMean_pt45_low_90pc,
                            DuneCrestMean_pt45_low_90pc_backerosion_1m,
                        ],
                        shoreline_position=[
                            shoreline_position_pt45_nat,
                            shoreline_position_pt45_low_40pc,
                            shoreline_position_pt45_low_90pc,
                            shoreline_position_pt45_low_90pc_backerosion_1m,
                        ],
                        overwash=[
                            overwash_pt45_nat,
                            overwash_pt45_low_40pc,
                            overwash_pt45_low_90pc,
                            overwash_pt45_low_90pc_backerosion_1m,
                        ],
                        dune_toe=[
                            [0],  # dummy
                            dune_toe_pt45_low_40pc,
                            dune_toe_pt45_low_90pc,
                            dune_toe_pt45_low_90pc_backerosion_1m,
                        ],
                        TMAX=[
                            1000,
                            1000,
                            1000,
                            1000,
                        ],
                        tmax_management=[
                            1000,  # dummy
                            407,
                            302,
                            302,
                        ],
                        roadways_on=False,
                        nourishment_on=True,
                        rebuild_threshold=rebuild_threshold,
                        # min dune height above the berm [m MHW], same as in RoadwayManager
                        scenarios=[
                            "natural",
                            "residential",
                            "commercial",
                            "comm, 1 m/yr",
                        ],
                    )

                def statistics():
                    (
                        year_abandoned,
                        sim_max,
                        overwash_filtered_removed,
                        dune_rebuilt,
                        beach_nourished,
                        _,
                        _,
                    ) = get_nourishment_statistics(
                        folder_prefix="",
                        natural_barrier_elev=None,
                        natural_barrier_width=None,
                        # individual_fid="8-B3D_Rave_pt45_Nourishment_2mDune_lowEle_residential_rebuild1m",
                        # individual_fid="8-B3D_Rave_pt45_Nourishment_2mDune_lowEle_commercial_rebuild1m",
                        individual_fid="8-B3D_Rave_pt45_Nourishment_2mDune_lowEle_commercial_backerosion1m_rebuild1m",
                        tmax=200,
                        iB3D=0,
                    )

                    nourishment_frequency = 200 / np.array(beach_nourished)

            def pt45_high():
                (
                    BarrierWidth_pt45_nat,
                    DuneCrestMean_pt45_nat,
                    BarrierHeight_pt45_nat,
                    bh_rate_pt45_nat,
                    bw_rate_pt45_nat,
                    sc_rate_pt45_nat,
                    DuneCrestMin_pt45_nat,
                    DuneCrestMax_pt45_nat,
                    shoreline_position_pt45_nat,
                    shoreface_slope_pt45_nat,
                    overwash_pt45_nat,
                    cascade_pt45_nat,
                ) = PLOT_5_Nonlinear_Dynamics_CASCADE_B3Donly_RoadwayManager(
                    name_prefix="4-B3D_Rave_pt45_Natural_high",
                    tmax_roadways=1000,
                    tmax_sim=1000,
                    plot_name="b3d_pt45_plots_high",
                    run_road_mgmt=False,
                    gif_on=False,
                )

                (
                    BarrierWidth_pt45_high_40pc,
                    DuneCrestMean_pt45_high_40pc,
                    BarrierHeight_pt45_high_40pc,
                    bh_rate_pt45_high_40pc,
                    bw_rate_pt45_high_40pc,
                    sc_rate_pt45_high_40pc,
                    DuneCrestMin_pt45_high_40pc,
                    DuneCrestMax_pt45_high_40pc,
                    shoreline_position_pt45_high_40pc,
                    shoreface_slope_pt45_high_40pc,
                    beach_width_pt45_high_40pc,
                    overwash_pt45_high_40pc,
                    dune_toe_pt45_high_40pc,
                    cascade_pt45_high_40pc,
                ) = PLOT_6_Nonlinear_Dynamics_CASCADE_B3Donly_Nourishments(
                    name_prefix="8-B3D_Rave_pt45_Nourishment_2mDune_highEle_residential_rebuild1m",
                    tmax_management=544,  # Community reached minimum width, drowned at 544 years
                    tmax_sim=574,  # Barrier has HEIGHT DROWNED at t = 574 years
                    plot_name="b3d_pt45_Nourishment_2mDune_highEle_residential_rebuild1m",
                    gif_on=False,
                    rebuild_dune_threshold=1,
                )

                (
                    BarrierWidth_pt45_high_90pc,
                    DuneCrestMean_pt45_high_90pc,
                    BarrierHeight_pt45_high_90pc,
                    bh_rate_pt45_high_90pc,
                    bw_rate_pt45_high_90pc,
                    sc_rate_pt45_high_90pc,
                    DuneCrestMin_pt45_high_90pc,
                    DuneCrestMax_pt45_high_90pc,
                    shoreline_position_pt45_high_90pc,
                    shoreface_slope_pt45_high_90pc,
                    beach_width_pt45_high_90pc,
                    overwash_pt45_high_90pc,
                    dune_toe_pt45_high_90pc,
                    cascade_pt45_high_90pc,
                ) = PLOT_6_Nonlinear_Dynamics_CASCADE_B3Donly_Nourishments(
                    name_prefix="8-B3D_Rave_pt45_Nourishment_2mDune_highEle_commercial_rebuild1m",
                    tmax_management=429,  # Community reached minimum width, drowned at 429 years
                    tmax_sim=1000,
                    plot_name="b3d_pt45_Nourishment_2mDune_highEle_commercial_rebuild1m",
                    gif_on=False,
                    rebuild_dune_threshold=1,
                )

                (
                    BarrierWidth_pt45_high_90pc_backerosion_1m,
                    DuneCrestMean_pt45_high_90pc_backerosion_1m,
                    BarrierHeight_pt45_high_90pc_backerosion_1m,
                    bh_rate_pt45_high_90pc_backerosion_1m,
                    bw_rate_pt45_high_90pc_backerosion_1m,
                    sc_rate_pt45_high_90pc_backerosion_1m,
                    DuneCrestMin_pt45_high_90pc_backerosion_1m,
                    DuneCrestMax_pt45_high_90pc_backerosion_1m,
                    shoreline_position_pt45_high_90pc_backerosion_1m,
                    shoreface_slope_pt45_high_90pc_backerosion_1m,
                    beach_width_pt45_high_90pc_backerosion_1m,
                    overwash_pt45_high_90pc_backerosion_1m,
                    dune_toe_pt45_high_90pc_backerosion_1m,
                    cascade_pt45_high_90pc_backerosion_1m,
                ) = PLOT_6_Nonlinear_Dynamics_CASCADE_B3Donly_Nourishments(
                    name_prefix="8-B3D_Rave_pt45_Nourishment_2mDune_highEle_commercial_backerosion1m_rebuild1m",
                    tmax_management=429,
                    tmax_sim=1000,
                    plot_name="b3d_pt45_Nourishment_2mDune_highEle_commercial_backerosion_1m_rebuild1m",
                    gif_on=False,
                    rebuild_dune_threshold=1,
                )

                rebuild_threshold = 1 + (
                    cascade_pt45_high_40pc.barrier3d[0].BermEl * 10
                )

                def version1_residential():
                    CASCADEplt.plot_nonlinear_stats_mgmt_array4(
                        cascade=[
                            cascade_pt45_nat,
                            cascade_pt45_high_40pc,  # only residential here
                            cascade_pt45_high_90pc,
                            cascade_pt45_high_90pc_backerosion_1m,
                        ],
                        DuneCrestMin=[
                            DuneCrestMin_pt45_nat,
                            DuneCrestMin_pt45_high_40pc,
                            DuneCrestMin_pt45_high_90pc,
                            DuneCrestMin_pt45_high_90pc_backerosion_1m,
                        ],
                        DuneCrestMax=[
                            DuneCrestMax_pt45_nat,
                            DuneCrestMax_pt45_high_40pc,
                            DuneCrestMax_pt45_high_90pc,
                            DuneCrestMax_pt45_high_90pc_backerosion_1m,
                        ],
                        BarrierHeight=[
                            BarrierHeight_pt45_nat,
                            BarrierHeight_pt45_high_40pc,
                            BarrierHeight_pt45_high_90pc,
                            BarrierHeight_pt45_high_90pc_backerosion_1m,
                        ],
                        BarrierWidth=[
                            BarrierWidth_pt45_nat,
                            BarrierWidth_pt45_high_40pc,
                            BarrierWidth_pt45_high_90pc,
                            BarrierWidth_pt45_high_90pc_backerosion_1m,
                        ],
                        DuneCrestMean=[
                            DuneCrestMean_pt45_nat,
                            DuneCrestMean_pt45_high_40pc,
                            DuneCrestMean_pt45_high_90pc,
                            DuneCrestMean_pt45_high_90pc_backerosion_1m,
                        ],
                        shoreline_position=[
                            shoreline_position_pt45_nat,
                            shoreline_position_pt45_high_40pc,
                            shoreline_position_pt45_high_90pc,
                            shoreline_position_pt45_high_90pc_backerosion_1m,
                        ],
                        overwash=[
                            overwash_pt45_nat,
                            overwash_pt45_high_40pc,
                            overwash_pt45_high_90pc,
                            overwash_pt45_high_90pc_backerosion_1m,
                        ],
                        dune_toe=[
                            [0],  # dummy
                            dune_toe_pt45_high_40pc,
                            dune_toe_pt45_high_90pc,
                            dune_toe_pt45_high_90pc_backerosion_1m,
                        ],
                        TMAX=[
                            1000,  # was 800, switched to 750 to match roadways
                            574,
                            1000,
                            1000,
                        ],
                        tmax_management=[
                            1000,  # dummy
                            544,
                            429,
                            429,
                        ],
                        roadways_on=False,
                        nourishment_on=True,
                        rebuild_threshold=rebuild_threshold,
                        # min dune height above the berm [m MHW], same as in RoadwayManager
                        scenarios=[
                            "natural",
                            "residential",
                            "commercial",
                            "comm, 1 m/yr",
                        ],
                    )

            def pt75_low():
                (
                    BarrierWidth_pt75_nat,
                    DuneCrestMean_pt75_nat,
                    BarrierHeight_pt75_nat,
                    bh_rate_pt75_nat,
                    bw_rate_pt75_nat,
                    sc_rate_pt75_nat,
                    DuneCrestMin_pt75_nat,
                    DuneCrestMax_pt75_nat,
                    shoreline_position_pt75_nat,
                    shoreface_slope_pt75_nat,
                    overwash_pt75_nat,
                    cascade_pt75_nat,
                ) = PLOT_5_Nonlinear_Dynamics_CASCADE_B3Donly_RoadwayManager(
                    name_prefix="4-B3D_Rave_pt75_Natural_low",
                    tmax_roadways=1000,  # dummy
                    tmax_sim=1000,
                    plot_name="b3d_pt75_plots_low",
                    run_road_mgmt=False,
                    gif_on=False,
                )

                (
                    BarrierWidth_pt75_low_40pc,
                    DuneCrestMean_pt75_low_40pc,
                    BarrierHeight_pt75_low_40pc,
                    bh_rate_pt75_low_40pc,
                    bw_rate_pt75_low_40pc,
                    sc_rate_pt75_low_40pc,
                    DuneCrestMin_pt75_low_40pc,
                    DuneCrestMax_pt75_low_40pc,
                    shoreline_position_pt75_low_40pc,
                    shoreface_slope_pt75_low_40pc,
                    beach_width_pt75_low_40pc,
                    overwash_pt75_low_40pc,
                    dune_toe_pt75_low_40pc,
                    cascade_pt75_low_40pc,
                ) = PLOT_6_Nonlinear_Dynamics_CASCADE_B3Donly_Nourishments(
                    name_prefix="8-B3D_Rave_pt75_Nourishment_2mDune_lowEle_residential_rebuild1m",
                    tmax_management=160,  # Community reached minimum width, drowned at 160 years
                    tmax_sim=1000,
                    plot_name="b3d_pt75_Nourishment_2mDune_lowEle_residential_rebuild1m",
                    gif_on=False,
                    rebuild_dune_threshold=1,
                )

                (
                    BarrierWidth_pt75_low_90pc,
                    DuneCrestMean_pt75_low_90pc,
                    BarrierHeight_pt75_low_90pc,
                    bh_rate_pt75_low_90pc,
                    bw_rate_pt75_low_90pc,
                    sc_rate_pt75_low_90pc,
                    DuneCrestMin_pt75_low_90pc,
                    DuneCrestMax_pt75_low_90pc,
                    shoreline_position_pt75_low_90pc,
                    shoreface_slope_pt75_low_90pc,
                    beach_width_pt75_low_90pc,
                    overwash_pt75_low_90pc,
                    dune_toe_pt75_low_90pc,
                    cascade_pt75_low_90pc,
                ) = PLOT_6_Nonlinear_Dynamics_CASCADE_B3Donly_Nourishments(
                    name_prefix="8-B3D_Rave_pt75_Nourishment_2mDune_lowEle_commercial_rebuild1m",
                    tmax_management=83,  # Community reached minimum width, drowned at 83 years
                    tmax_sim=1000,
                    plot_name="b3d_pt75_Nourishment_2mDune_lowEle_commercial_rebuild1m",
                    gif_on=False,
                    rebuild_dune_threshold=1,
                )

                (
                    BarrierWidth_pt75_low_90pc_backerosion_1m,
                    DuneCrestMean_pt75_low_90pc_backerosion_1m,
                    BarrierHeight_pt75_low_90pc_backerosion_1m,
                    bh_rate_pt75_low_90pc_backerosion_1m,
                    bw_rate_pt75_low_90pc_backerosion_1m,
                    sc_rate_pt75_low_90pc_backerosion_1m,
                    DuneCrestMin_pt75_low_90pc_backerosion_1m,
                    DuneCrestMax_pt75_low_90pc_backerosion_1m,
                    shoreline_position_pt75_low_90pc_backerosion_1m,
                    shoreface_slope_pt75_low_90pc_backerosion_1m,
                    beach_width_pt75_low_90pc_backerosion_1m,
                    overwash_pt75_low_90pc_backerosion_1m,
                    dune_toe_pt75_low_90pc_backerosion_1m,
                    cascade_pt75_low_90pc_backerosion_1m,
                ) = PLOT_6_Nonlinear_Dynamics_CASCADE_B3Donly_Nourishments(
                    name_prefix="8-B3D_Rave_pt75_Nourishment_2mDune_lowEle_commercial_backerosion1m_rebuild1m",
                    tmax_management=83,  # Community reached minimum width, drowned at 83 years
                    tmax_sim=1000,
                    plot_name="b3d_pt75_Nourishment_2mDune_lowEle_commercial_backerosion_1m_rebuild1m",
                    gif_on=False,
                    rebuild_dune_threshold=1,
                )

                rebuild_threshold = 1 + (cascade_pt75_low_40pc.barrier3d[0].BermEl * 10)

                def version1_res_vs_commercial():
                    CASCADEplt.plot_nonlinear_stats_mgmt_array4(
                        cascade=[
                            cascade_pt75_nat,
                            cascade_pt75_low_40pc,  # only residential here
                            cascade_pt75_low_90pc,
                            cascade_pt75_low_90pc_backerosion_1m,
                        ],
                        DuneCrestMin=[
                            DuneCrestMin_pt75_nat,
                            DuneCrestMin_pt75_low_40pc,
                            DuneCrestMin_pt75_low_90pc,
                            DuneCrestMin_pt75_low_90pc_backerosion_1m,
                        ],
                        DuneCrestMax=[
                            DuneCrestMax_pt75_nat,
                            DuneCrestMax_pt75_low_40pc,
                            DuneCrestMax_pt75_low_90pc,
                            DuneCrestMax_pt75_low_90pc_backerosion_1m,
                        ],
                        BarrierHeight=[
                            BarrierHeight_pt75_nat,
                            BarrierHeight_pt75_low_40pc,
                            BarrierHeight_pt75_low_90pc,
                            BarrierHeight_pt75_low_90pc_backerosion_1m,
                        ],
                        BarrierWidth=[
                            BarrierWidth_pt75_nat,
                            BarrierWidth_pt75_low_40pc,
                            BarrierWidth_pt75_low_90pc,
                            BarrierWidth_pt75_low_90pc_backerosion_1m,
                        ],
                        DuneCrestMean=[
                            DuneCrestMean_pt75_nat,
                            DuneCrestMean_pt75_low_40pc,
                            DuneCrestMean_pt75_low_90pc,
                            DuneCrestMean_pt75_low_90pc_backerosion_1m,
                        ],
                        shoreline_position=[
                            shoreline_position_pt75_nat,
                            shoreline_position_pt75_low_40pc,
                            shoreline_position_pt75_low_90pc,
                            shoreline_position_pt75_low_90pc_backerosion_1m,
                        ],
                        overwash=[
                            overwash_pt75_nat,
                            overwash_pt75_low_40pc,
                            overwash_pt75_low_90pc,
                            overwash_pt75_low_90pc_backerosion_1m,
                        ],
                        dune_toe=[
                            [0],  # dummy
                            dune_toe_pt75_low_40pc,
                            dune_toe_pt75_low_90pc,
                            dune_toe_pt75_low_90pc_backerosion_1m,
                        ],
                        TMAX=[
                            1000,
                            1000,
                            1000,
                            1000,
                        ],
                        tmax_management=[
                            1000,  # dummy
                            160,
                            83,
                            83,
                        ],
                        roadways_on=False,
                        nourishment_on=True,
                        rebuild_threshold=rebuild_threshold,  # min dune height above the berm [m MHW], same as in RoadwayManager
                        scenarios=[
                            "natural",
                            "residential",
                            "commercial",
                            "comm, 1 m/yr",
                        ],
                    )

            def pt75_high():
                (
                    BarrierWidth_pt75_nat,
                    DuneCrestMean_pt75_nat,
                    BarrierHeight_pt75_nat,
                    bh_rate_pt75_nat,
                    bw_rate_pt75_nat,
                    sc_rate_pt75_nat,
                    DuneCrestMin_pt75_nat,
                    DuneCrestMax_pt75_nat,
                    shoreline_position_pt75_nat,
                    shoreface_slope_pt75_nat,
                    overwash_pt75_nat,
                    cascade_pt75_nat,
                ) = PLOT_5_Nonlinear_Dynamics_CASCADE_B3Donly_RoadwayManager(
                    name_prefix="4-B3D_Rave_pt75_Natural_high",
                    tmax_roadways=1000,  # dummy
                    tmax_sim=1000,
                    plot_name="b3d_pt75_plots_high",
                    run_road_mgmt=False,
                    gif_on=False,
                )

                (
                    BarrierWidth_pt75_high_40pc,
                    DuneCrestMean_pt75_high_40pc,
                    BarrierHeight_pt75_high_40pc,
                    bh_rate_pt75_high_40pc,
                    bw_rate_pt75_high_40pc,
                    sc_rate_pt75_high_40pc,
                    DuneCrestMin_pt75_high_40pc,
                    DuneCrestMax_pt75_high_40pc,
                    shoreline_position_pt75_high_40pc,
                    shoreface_slope_pt75_high_40pc,
                    beach_width_pt75_high_40pc,
                    overwash_pt75_high_40pc,
                    dune_toe_pt75_high_40pc,
                    cascade_pt75_high_40pc,
                ) = PLOT_6_Nonlinear_Dynamics_CASCADE_B3Donly_Nourishments(
                    name_prefix="8-B3D_Rave_pt75_Nourishment_2mDune_highEle_residential_rebuild1m",
                    tmax_management=550,  # Community reached minimum width, drowned at 550
                    tmax_sim=1000,
                    plot_name="b3d_pt75_Nourishment_2mDune_highEle_residential_rebuild1m",
                    gif_on=False,
                    rebuild_dune_threshold=1,
                )

                (
                    BarrierWidth_pt75_high_90pc,
                    DuneCrestMean_pt75_high_90pc,
                    BarrierHeight_pt75_high_90pc,
                    bh_rate_pt75_high_90pc,
                    bw_rate_pt75_high_90pc,
                    sc_rate_pt75_high_90pc,
                    DuneCrestMin_pt75_high_90pc,
                    DuneCrestMax_pt75_high_90pc,
                    shoreline_position_pt75_high_90pc,
                    shoreface_slope_pt75_high_90pc,
                    beach_width_pt75_high_90pc,
                    overwash_pt75_high_90pc,
                    dune_toe_pt75_high_90pc,
                    cascade_pt75_high_90pc,
                ) = PLOT_6_Nonlinear_Dynamics_CASCADE_B3Donly_Nourishments(
                    name_prefix="8-B3D_Rave_pt75_Nourishment_2mDune_highEle_commercial_rebuild1m",
                    tmax_management=518,  # Community reached minimum width, drowned at 518
                    tmax_sim=580,  # Barrier has HEIGHT DROWNED at t = 580 years
                    plot_name="b3d_pt75_Nourishment_2mDune_highEle_commercial_rebuild1m",
                    gif_on=False,
                    rebuild_dune_threshold=1,
                )

                (
                    BarrierWidth_pt75_high_90pc_backerosion_1m,
                    DuneCrestMean_pt75_high_90pc_backerosion_1m,
                    BarrierHeight_pt75_high_90pc_backerosion_1m,
                    bh_rate_pt75_high_90pc_backerosion_1m,
                    bw_rate_pt75_high_90pc_backerosion_1m,
                    sc_rate_pt75_high_90pc_backerosion_1m,
                    DuneCrestMin_pt75_high_90pc_backerosion_1m,
                    DuneCrestMax_pt75_high_90pc_backerosion_1m,
                    shoreline_position_pt75_high_90pc_backerosion_1m,
                    shoreface_slope_pt75_high_90pc_backerosion_1m,
                    beach_width_pt75_high_90pc_backerosion_1m,
                    overwash_pt75_high_90pc_backerosion_1m,
                    dune_toe_pt75_high_90pc_backerosion_1m,
                    cascade_pt75_high_90pc_backerosion_1m,
                ) = PLOT_6_Nonlinear_Dynamics_CASCADE_B3Donly_Nourishments(
                    name_prefix="8-B3D_Rave_pt75_Nourishment_2mDune_highEle_commercial_backerosion1m_rebuild1m",
                    tmax_management=518,  # Community reached minimum width, drowned at 518 years
                    tmax_sim=1000,
                    plot_name="b3d_pt75_Nourishment_2mDune_highEle_commercial_backerosion_1m_rebuild1m",
                    gif_on=False,
                    rebuild_dune_threshold=1,
                )

                rebuild_threshold = 1 + (
                    cascade_pt75_high_40pc.barrier3d[0].BermEl * 10
                )

                def version1_res_vs_commercial():
                    CASCADEplt.plot_nonlinear_stats_mgmt_array4(
                        cascade=[
                            cascade_pt75_nat,
                            cascade_pt75_high_40pc,  # only residential here
                            cascade_pt75_high_90pc,
                            cascade_pt75_high_90pc_backerosion_1m,
                        ],
                        DuneCrestMin=[
                            DuneCrestMin_pt75_nat,
                            DuneCrestMin_pt75_high_40pc,
                            DuneCrestMin_pt75_high_90pc,
                            DuneCrestMin_pt75_high_90pc_backerosion_1m,
                        ],
                        DuneCrestMax=[
                            DuneCrestMax_pt75_nat,
                            DuneCrestMax_pt75_high_40pc,
                            DuneCrestMax_pt75_high_90pc,
                            DuneCrestMax_pt75_high_90pc_backerosion_1m,
                        ],
                        BarrierHeight=[
                            BarrierHeight_pt75_nat,
                            BarrierHeight_pt75_high_40pc,
                            BarrierHeight_pt75_high_90pc,
                            BarrierHeight_pt75_high_90pc_backerosion_1m,
                        ],
                        BarrierWidth=[
                            BarrierWidth_pt75_nat,
                            BarrierWidth_pt75_high_40pc,
                            BarrierWidth_pt75_high_90pc,
                            BarrierWidth_pt75_high_90pc_backerosion_1m,
                        ],
                        DuneCrestMean=[
                            DuneCrestMean_pt75_nat,
                            DuneCrestMean_pt75_high_40pc,
                            DuneCrestMean_pt75_high_90pc,
                            DuneCrestMean_pt75_high_90pc_backerosion_1m,
                        ],
                        shoreline_position=[
                            shoreline_position_pt75_nat,
                            shoreline_position_pt75_high_40pc,
                            shoreline_position_pt75_high_90pc,
                            shoreline_position_pt75_high_90pc_backerosion_1m,
                        ],
                        overwash=[
                            overwash_pt75_nat,
                            overwash_pt75_high_40pc,
                            overwash_pt75_high_90pc,
                            overwash_pt75_high_90pc_backerosion_1m,
                        ],
                        dune_toe=[
                            [0],  # dummy
                            dune_toe_pt75_high_40pc,
                            dune_toe_pt75_high_90pc,
                            dune_toe_pt75_high_90pc_backerosion_1m,
                        ],
                        TMAX=[
                            1000,
                            1000,
                            580,
                            1000,
                        ],
                        tmax_management=[
                            1000,  # dummy
                            550,
                            518,
                            1000,
                        ],
                        roadways_on=False,
                        nourishment_on=True,
                        rebuild_threshold=rebuild_threshold,  # min dune height above the berm [m MHW], same as in RoadwayManager
                        scenarios=[
                            "natural",
                            "residential",
                            "commercial",
                            "comm, 1 m/yr",
                        ],
                    )

    def initial_topos():

        PLOT_7_Initial_CNH_Topographies(
            [
                "8-B3D_Rave_pt45_Nourishment_2mDune_lowEle_commercial_rebuild1m",
                "8-B3D_Rave_pt45_Nourishment_2mDune_highEle_commercial_rebuild1m",
                "8-B3D_Rave_pt75_Nourishment_2mDune_lowEle_commercial_rebuild1m",
                "8-B3D_Rave_pt75_Nourishment_2mDune_highEle_commercial_rebuild1m",
            ]
        )

    def ast_connections():
        def half_nourishment_half_roadways():
            # half nourishment, half roadways, background erosion of 1 m
            # All 3: Community reached minimum width, drowned at 83 years
            # All 3: Island is to narrow for roadway to be relocated. Roadway eaten up by dunes at 143
            (
                barrier_width_linSLR_pt75low,
                dune_crest_mean_linSLR_pt75low,
                barrier_height_linSLR_pt75low,
                bh_rate_linSLR_pt75low,
                bw_rate_linSLR_pt75low,
                sc_rate_linSLR_pt75low,
                dune_crest_min_linSLR_pt75low,
                dune_crest_max_linSLR_pt75low,
                shoreline_position_linSLR_pt75low,
                shoreface_slope_linSLR_pt75low,
                beach_width_linSLR_pt75low,
                overwash_linSLR_pt75low,
                dune_toe_linSLR_pt75low,
                cascade_linSLR_pt75low,
            ) = PLOT_9_Nonlinear_Dynamics_CASCADE_AST(
                name_prefix="9-CASCADE_Rave_pt75_Nourish_2mDune_lowEle_comm_BE1m_RT1m_6AST_3roads",
                tmax_management=[
                    83,
                    83,
                    83,
                    143,
                    143,
                    143,
                ],  # an array, CAN'T BE 200 OR THE ALGORITHM BREAKS
                tmax_sim=200,  # not an array
                plot_name="9-CASCADE_AST_nourishment_roadways_3km",
                rebuild_dune_threshold=1,
                beach_management_ny=[1, 1, 1, 0, 0, 0],
                roadway_management_ny=[0, 0, 0, 1, 1, 1],
                gif_on=False,
                z_lim=4,
                y_lim=[150, 220],
                fig_size=(6, 2.5),
                time_series_on=True,
                fig_eps=True,
            )

            # half nourishment, half roadways, background erosion of 1 m, accelerated SLR
            # All 3: Community reached minimum width, drowned at 64 years;
            # Barrier has HEIGHT DROWNED at t = 85 years
            (
                barrier_width_accSLR_pt75low,
                dune_crest_mean_accSLR_pt75low,
                barrier_height_accSLR_pt75low,
                bh_rate_accSLR_pt75low,
                bw_rate_accSLR_pt75low,
                sc_rate_accSLR_pt75low,
                dune_crest_min_accSLR_pt75low,
                dune_crest_max_accSLR_pt75low,
                shoreline_position_accSLR_pt75low,
                shoreface_slope_accSLR_pt75low,
                beach_width_accSLR_pt75low,
                overwash_accSLR_pt75low,
                dune_toe_accSLR_pt75low,
                cascade_accSLR_pt75low,
            ) = PLOT_9_Nonlinear_Dynamics_CASCADE_AST(
                name_prefix="9-CASCADE_Rave_pt75_Nourish_2mDune_lowEle_comm_BE1m_RT1m_6AST_3roads_AccSLR",
                tmax_management=[
                    64,
                    64,
                    64,
                    84,
                    84,
                    84,
                ],  # an array, CAN'T BE 100 OR THE ALGORITHM BREAKS
                tmax_sim=85,  # not an array
                plot_name="9-CASCADE_AST_nourishment_roadways_3km_AccSLR",
                rebuild_dune_threshold=1,
                beach_management_ny=[1, 1, 1, 0, 0, 0],
                roadway_management_ny=[0, 0, 0, 1, 1, 1],
                gif_on=False,
                z_lim=4,
                y_lim=[150, 220],
                fig_size=(6, 2.5),
                time_series_on=True,
                fig_eps=True,
            )

            # half nourishment, half roadways, background erosion of 1 m -- different elevations (pt45 left, pt75 right)
            # All 3: Roadway drowned in place at 132 years due to SLR - road cannot be below 0 m MHW
            (
                barrier_width_pt45nourish_linSLR,
                dune_crest_mean_pt45nourish_linSLR,
                barrier_height_pt45nourish_linSLR,
                bh_rate_pt45nourish_linSLR,
                bw_rate_pt45nourish_linSLR,
                sc_rate_pt45nourish_linSLR,
                dune_crest_min_pt45nourish_linSLR,
                dune_crest_max_pt45nourish_linSLR,
                shoreline_position_pt45nourish_linSLR,
                shoreface_slope_pt45nourish_linSLR,
                beach_width_pt45nourish_linSLR,
                overwash_pt45nourish_linSLR,
                dune_toe_pt45nourish_linSLR,
                cascade_pt45nourish_linSLR,
            ) = PLOT_9_Nonlinear_Dynamics_CASCADE_AST(
                name_prefix="9-CASCADE_Rave_pt45_pt75_low_split_2mDune_comm_BE1m_RT1m_6AST_3roads",
                tmax_management=[
                    199,
                    199,
                    199,
                    132,
                    132,
                    132,
                ],  # an array, CAN'T BE 200 OR THE ALGORITHM BREAKS
                tmax_sim=200,  # not an array
                plot_name="9-CASCADE_AST_pt45nourishment_pt75roadways_3km",
                rebuild_dune_threshold=1,
                beach_management_ny=[1, 1, 1, 0, 0, 0],
                roadway_management_ny=[0, 0, 0, 1, 1, 1],
                gif_on=False,
                z_lim=4,
                y_lim=[150, 220],
                fig_size=(6, 2.5),
                time_series_on=True,
                fig_eps=True,
            )

            # half nourishment, half roadways, BE of 1 m -- different elevations (pt45 left, pt75 right) -- ACC SLR
            # All three: Roadway drowned in place at 85 years due to SLR - road cannot be below 0 m MHW
            # All three: Community reached minimum width, drowned at 137 years
            # Barrier has HEIGHT DROWNED at t = 141 years
            (
                barrier_width_pt45nourish_accSLR,
                dune_crest_mean_pt45nourish_accSLR,
                barrier_height_pt45nourish_accSLR,
                bh_rate_pt45nourish_accSLR,
                bw_rate_pt45nourish_accSLR,
                sc_rate_pt45nourish_accSLR,
                dune_crest_min_pt45nourish_accSLR,
                dune_crest_max_pt45nourish_accSLR,
                shoreline_position_pt45nourish_accSLR,
                shoreface_slope_pt45nourish_accSLR,
                beach_width_pt45nourish_accSLR,
                overwash_pt45nourish_accSLR,
                dune_toe_pt45nourish_accSLR,
                cascade_pt45nourish_accSLR,
            ) = PLOT_9_Nonlinear_Dynamics_CASCADE_AST(
                name_prefix="9-CASCADE_Rave_pt45_pt75_low_split_2mDune_comm_BE1m_RT1m_6AST_3roads_AccSLR",
                tmax_management=[
                    137,
                    137,
                    137,
                    85,
                    85,
                    85,
                ],  # an array, CAN'T BE 200 OR THE ALGORITHM BREAKS
                tmax_sim=141,  # not an array
                plot_name="9-CASCADE_AST_pt45nourishment_pt75roadways_3km_AccSLR",
                rebuild_dune_threshold=1,
                beach_management_ny=[1, 1, 1, 0, 0, 0],
                roadway_management_ny=[0, 0, 0, 1, 1, 1],
                gif_on=False,
                z_lim=4,
                y_lim=[150, 220],
                fig_size=(6, 2.5),
                time_series_on=True,
                fig_eps=True,
            )

            def ast_baselines():
                # all roadways, linear SLR, backround erosion 1 m (pt75)
                # Island is to narrow for roadway to be relocated. Roadway eaten up by dunes at 73
                (
                    barrier_width_allroads_pt75low,
                    dune_crest_mean_allroads_pt75low,
                    barrier_height_allroads_pt75low,
                    bh_rate_allroads_pt75low,
                    bw_rate_allroads_pt75low,
                    sc_rate_allroads_pt75low,
                    dune_crest_min_allroads_pt75low,
                    dune_crest_max_allroads_pt75low,
                    shoreline_position_allroads_pt75low,
                    shoreface_slope_allroads_pt75low,
                    beach_width_allroads_pt75low,
                    overwash_allroads_pt75low,
                    dune_toe_allroads_pt75low,
                    cascade_allroads_pt75low,
                ) = PLOT_9_Nonlinear_Dynamics_CASCADE_AST(
                    name_prefix="9-CASCADE_Rave_pt75_2mDune_lowEle_BE1m_6AST_6roads",
                    tmax_management=[
                        73,
                        73,
                        73,
                        73,
                        73,
                        73,
                    ],  # an array, CAN'T BE 200 OR THE ALGORITHM BREAKS
                    tmax_sim=200,  # not an array
                    plot_name="9-CASCADE_AST_all_roadways_3km",
                    rebuild_dune_threshold=1,
                    beach_management_ny=[0, 0, 0, 0, 0, 0],
                    roadway_management_ny=[1, 1, 1, 1, 1, 1],
                    gif_on=False,
                    z_lim=4,
                    y_lim=[150, 220],
                    fig_size=(6, 2.5),
                    time_series_on=True,
                )

                # all nourishments, linear SLR, backround erosion 1 m (pt75)
                # Community reached minimum width, drowned at 83 years
                (
                    barrier_width_allnourish_pt75low,
                    dune_crest_mean_allnourish_pt75low,
                    barrier_height_allnourish_pt75low,
                    bh_rate_allnourish_pt75low,
                    bw_rate_allnourish_pt75low,
                    sc_rate_allnourish_pt75low,
                    dune_crest_min_allnourish_pt75low,
                    dune_crest_max_allnourish_pt75low,
                    shoreline_position_allnourish_pt75low,
                    shoreface_slope_allnourish_pt75low,
                    beach_width_allnourish_pt75low,
                    overwash_allnourish_pt75low,
                    dune_toe_allnourish_pt75low,
                    cascade_allnourish_pt75low,
                ) = PLOT_9_Nonlinear_Dynamics_CASCADE_AST(
                    name_prefix="9-CASCADE_Rave_pt75_2mDune_lowEle_comm_BE1m_RT1m_6AST_6nourish",
                    tmax_management=[
                        83,
                        83,
                        83,
                        83,
                        83,
                        83,
                    ],  # an array, CAN'T BE 200 OR THE ALGORITHM BREAKS
                    tmax_sim=200,  # not an array
                    plot_name="9-CASCADE_AST_all_nourishments_3km_pt75low",
                    rebuild_dune_threshold=1,
                    beach_management_ny=[1, 1, 1, 1, 1, 1],
                    roadway_management_ny=[0, 0, 0, 0, 0, 0],
                    gif_on=False,
                    z_lim=4,
                    y_lim=[150, 220],
                    fig_size=(6, 2.5),
                    time_series_on=True,
                )

                # all nourishments pt45, linear SLR, backround erosion 1 m
                (
                    barrier_width_allnourish_pt45low,
                    dune_crest_mean_allnourish_pt45low,
                    barrier_height_allnourish_pt45low,
                    bh_rate_allnourish_pt45low,
                    bw_rate_allnourish_pt45low,
                    sc_rate_allnourish_pt45low,
                    dune_crest_min_allnourish_pt45low,
                    dune_crest_max_allnourish_pt45low,
                    shoreline_position_allnourish_pt45low,
                    shoreface_slope_allnourish_pt45low,
                    beach_width_allnourish_pt45low,
                    overwash_allnourish_pt45low,
                    dune_toe_allnourish_pt45low,
                    cascade_allnourish_pt45low,
                ) = PLOT_9_Nonlinear_Dynamics_CASCADE_AST(
                    name_prefix="9-CASCADE_Rave_pt45_low_nourishment_2mDune_comm_BE1m_RT1m_6AST",
                    tmax_management=[
                        199,
                        199,
                        199,
                        199,
                        199,
                        199,
                    ],  # an array, CAN'T BE 200 OR THE ALGORITHM BREAKS
                    tmax_sim=200,  # not an array
                    plot_name="9-CASCADE_AST_all_nourishments_3km_pt45low",
                    rebuild_dune_threshold=1,
                    beach_management_ny=[1, 1, 1, 1, 1, 1],
                    roadway_management_ny=[0, 0, 0, 0, 0, 0],
                    gif_on=False,
                    z_lim=4,
                    y_lim=[150, 220],
                    fig_size=(6, 2.5),
                    time_series_on=True,
                )

            def ast_time_series():
                iB3D_roadways = 4
                iB3D_community = 1

                # pt75
                CASCADEplt.plot_nonlinear_stats_ast_array3(
                    shoreline_position=[
                        shoreline_position_allroads_pt75low[iB3D_roadways],
                        shoreline_position_linSLR_pt75low[iB3D_roadways],
                        shoreline_position_accSLR_pt75low[iB3D_roadways],
                    ],
                    beach_width=[
                        beach_width_allnourish_pt75low[iB3D_community],
                        beach_width_linSLR_pt75low[iB3D_community],
                        beach_width_accSLR_pt75low[iB3D_community],
                    ],
                    TMAX=[
                        200,
                        200,
                        85,
                    ],
                    tmax_management_nourishments=[
                        83,
                        83,
                        64,
                    ],
                    tmax_management_roadways=[
                        73,
                        143,
                        85,
                    ],
                    scenarios_beach_width=[
                        "comm +1 m/yr",
                        "split mgmt +1 m/yr",
                        "split mgmt +1 m/yr, acc SLR",
                    ],
                    scenarios_shoreline_position=[
                        "roads +1 m/yr",
                        "split mgmt +1 m/yr",
                        "split mgmt +1 m/yr, acc SLR",
                    ],
                )

                # pt45 and pt75 low
                CASCADEplt.plot_nonlinear_stats_ast_array3(
                    shoreline_position=[
                        shoreline_position_allroads_pt75low[iB3D_roadways],
                        shoreline_position_pt45nourish_linSLR[
                            iB3D_roadways
                        ],  # actually pt75 on right side of domain
                        shoreline_position_pt45nourish_accSLR[
                            iB3D_roadways
                        ],  # actually pt75 on right side of domain
                    ],
                    beach_width=[
                        beach_width_allnourish_pt45low[iB3D_community],
                        beach_width_pt45nourish_linSLR[iB3D_community],
                        beach_width_pt45nourish_accSLR[iB3D_community],
                    ],
                    TMAX=[
                        200,
                        200,
                        141,
                    ],
                    tmax_management_nourishments=[
                        199,
                        199,
                        137,
                    ],
                    tmax_management_roadways=[
                        73,
                        131,
                        85,
                    ],
                    scenarios_beach_width=[
                        "comm +1 m/yr",
                        "split mgmt +1 m/yr",
                        "split mgmt +1 m/yr, acc SLR",
                    ],
                    scenarios_shoreline_position=[
                        "roads +1 m/yr",
                        "split mgmt +1 m/yr",
                        "split mgmt +1 m/yr, acc SLR",
                    ],
                )

        def thirds():
            # one third nourishment, low road, high road + 1 m/yr background erosion
            # middle roadway drowned at 99, 115, 131 due to SLR -- created a curved shape
            # after grouping abandonment, all radways drown at 99
            (
                barrier_width_thirds_linSLR,
                dune_crest_mean_thirds_linSLR,
                barrier_height_thirds_linSLR,
                bh_rate_thirds_linSLR,
                bw_rate_thirds_linSLR,
                sc_rate_thirds_linSLR,
                dune_crest_min_thirds_linSLR,
                dune_crest_max_thirds_linSLR,
                shoreline_position_thirds_linSLR,
                shoreface_slope_thirds_linSLR,
                beach_width_thirds_linSLR,
                overwash_thirds_linSLR,
                dune_toe_thirds_linSLR,
                cascade_thirds_linSLR,
            ) = PLOT_9_Nonlinear_Dynamics_CASCADE_AST(
                name_prefix="9-CASCADE_AST_3domains_BE1m",
                tmax_management=[
                    199,
                    199,
                    199,
                    99,
                    99,
                    99,
                    199,
                    199,
                    199,
                ],  # an array, CAN'T BE 200 OR THE ALGORITHM BREAKS
                tmax_sim=200,  # not an array
                plot_name="9-CASCADE_AST_3domains_BE1m",
                rebuild_dune_threshold=1,
                beach_management_ny=[1, 1, 1, 0, 0, 0, 0, 0, 0],
                roadway_management_ny=[0, 0, 0, 1, 1, 1, 1, 1, 1],
                gif_on=False,
                z_lim=4,
                y_lim=[150, 220],
                fig_size=(6, 2.5),
                time_series_on=True,
                fig_eps=True,
            )

            # same as above, but with accelerated SLR and 1m background erosion
            # Barrier has HEIGHT DROWNED at t = 71 years (#5 B3D) - 4261
            (
                barrier_width_thirds_acc,
                dune_crest_mean_thirds_acc,
                barrier_height_thirds_acc,
                bh_rate_thirds_acc,
                bw_rate_thirds_acc,
                sc_rate_thirds_acc,
                dune_crest_min_thirds_acc,
                dune_crest_max_thirds_acc,
                shoreline_position_thirds_acc,
                shoreface_slope_thirds_acc,
                beach_width_thirds_acc,
                overwash_thirds_acc,
                dune_toe_thirds_acc,
                cascade_thirds_acc,
            ) = PLOT_9_Nonlinear_Dynamics_CASCADE_AST(
                name_prefix="9-CASCADE_AST_3domains_BE1m_AccSLR",
                tmax_management=[
                    70,
                    70,
                    70,
                    70,
                    70,
                    70,  # height drowned
                    70,
                    70,
                    70,
                ],  # an array, CAN'T BE 100 OR THE ALGORITHM BREAKS
                tmax_sim=71,  # not an array
                plot_name="9-CASCADE_AST_3domains_BE1m_AccSLR",
                rebuild_dune_threshold=1,
                beach_management_ny=[1, 1, 1, 0, 0, 0, 0, 0, 0],
                roadway_management_ny=[0, 0, 0, 1, 1, 1, 1, 1, 1],
                gif_on=False,
                z_lim=4,
                y_lim=[150, 220],
                fig_size=(6, 2.5),
                time_series_on=True,
                fig_eps=True,
            )

            # same as above, but with accelerated SLR and no background erosion
            # Barrier has HEIGHT DROWNED at t = 92 years
            # NOT USED
            PLOT_9_Nonlinear_Dynamics_CASCADE_AST(
                name_prefix="9-CASCADE_AST_3domains_AccSLR",
                tmax_management=[
                    91,
                    91,
                    91,
                    91,
                    91,
                    91,  # drowned
                    91,
                    91,
                    91,
                ],  # an array, CAN'T BE 100 OR THE ALGORITHM BREAKS
                tmax_sim=92,  # not an array
                plot_name="9-CASCADE_AST_3domains_AccSLR",
                rebuild_dune_threshold=1,
                beach_management_ny=[1, 1, 1, 0, 0, 0, 0, 0, 0],
                roadway_management_ny=[0, 0, 0, 1, 1, 1, 1, 1, 1],
                gif_on=True,
                z_lim=4,
                y_lim=[150, 220],
                fig_size=(6, 2.5),
                time_series_on=False,
            )

            # natural case in the middle
            # Roadway width drowned at 137 years, 20.0% of road borders water
            # All: Community reached minimum width, drowned at 137 years
            # Roadway width drowned at 142 years, 20.0% of road borders water
            # Island is too narrow for roadway to be relocated. Roadway eaten up by dunes at 149 years
            # after grouping abandonment, all roadways drown at 137
            (
                barrier_width_acc_nat,
                dune_crest_mean_acc_nat,
                barrier_height_acc_nat,
                bh_rate_acc_nat,
                bw_rate_acc_nat,
                sc_rate_acc_nat,
                dune_crest_min_acc_nat,
                dune_crest_max_acc_nat,
                shoreline_position_acc_nat,
                shoreface_slope_acc_nat,
                beach_width_acc_nat,
                overwash_acc_nat,
                dune_toe_acc_nat,
                cascade_acc_nat,
            ) = PLOT_9_Nonlinear_Dynamics_CASCADE_AST(
                name_prefix="9-CASCADE_AST_3domains_BE1m_AccSLR_nat_middle",
                tmax_management=[
                    137,
                    137,
                    137,
                    0,
                    0,
                    0,
                    137,
                    137,
                    137,
                ],  # an array, CAN'T BE 100 OR THE ALGORITHM BREAKS
                tmax_sim=200,  # not an array
                plot_name="9-CASCADE_AST_3domains_AccSLR_nat_middle",
                rebuild_dune_threshold=1,
                beach_management_ny=[1, 1, 1, 0, 0, 0, 0, 0, 0],
                roadway_management_ny=[0, 0, 0, 0, 0, 0, 1, 1, 1],
                gif_on=False,
                z_lim=4,
                y_lim=[150, 220],
                fig_size=(6, 2.5),
                time_series_on=True,
                fig_eps=True,
            )

            def ast_time_series():
                iB3D_roadways = 4
                iB3D_community = 1

                # pt45 and pt75 low
                CASCADEplt.plot_nonlinear_stats_ast_array3(
                    shoreline_position=[
                        shoreline_position_allroads_pt75low[
                            iB3D_roadways
                        ],  # this is a dummy
                        shoreline_position_thirds_acc[iB3D_roadways],
                        shoreline_position_acc_nat[iB3D_roadways],
                    ],
                    beach_width=[
                        beach_width_allnourish_pt45low[
                            iB3D_community
                        ],  # this is a dummy
                        beach_width_thirds_acc[iB3D_community],
                        beach_width_acc_nat[iB3D_community],
                    ],
                    TMAX=[
                        200,
                        71,
                        200,
                    ],
                    tmax_management_nourishments=[
                        199,
                        70,
                        137,
                    ],
                    tmax_management_roadways=[
                        73,
                        70,
                        137,
                    ],
                    scenarios_beach_width=[
                        "comm +1 m/yr",
                        "variable mgmt +1 m/yr, acc SLR",
                        "retreat +1 m/yr, acc SLR",
                    ],
                    scenarios_shoreline_position=[
                        "roads +1 m/yr",
                        "variable mgmt +1 m/yr, acc SLR",
                        "retreat +1 m/yr, acc SLR",
                    ],
                )

            def misc_stats():
                def nourishment_stats(cascade, cutoff, iB3D, nourishment_volume):
                    # nourishment statistics
                    nourishments = (
                        cascade.nourishments[iB3D].nourishment_volume_TS[:cutoff]
                        == nourishment_volume
                    )
                    nourishment_frequency = [i for i, x in enumerate(nourishments) if x]
                    nourishment_frequency_pre = [
                        y - x
                        for x, y in zip(
                            nourishment_frequency, nourishment_frequency[1:]
                        )
                    ]
                    mean_pre = np.mean(nourishment_frequency_pre)

                    nourishments = (
                        cascade.nourishments[iB3D].nourishment_volume_TS[cutoff:]
                        == nourishment_volume
                    )
                    nourishment_frequency = [i for i, x in enumerate(nourishments) if x]
                    nourishment_frequency_post = [
                        y - x
                        for x, y in zip(
                            nourishment_frequency, nourishment_frequency[1:]
                        )
                    ]
                    mean_post = np.mean(nourishment_frequency_post)

                    return mean_pre, mean_post

                mean_pre_pt45_linSLR, mean_post_pt45_linSLR = nourishment_stats(
                    cascade_pt45nourish_linSLR,
                    cutoff=132,
                    iB3D=0,
                    nourishment_volume=100,
                )

                mean_pre_thirds_linSLR, mean_post_thirds_linSLR = nourishment_stats(
                    cascade_thirds_linSLR, cutoff=99, iB3D=0, nourishment_volume=100
                )

                mean_pre_thirds_accSLR, mean_post_thirds_accSLR = nourishment_stats(
                    cascade_thirds_acc, cutoff=71, iB3D=0, nourishment_volume=100
                )

                (
                    mean_pre_thirds_nat_accSLR,
                    mean_post_thirds_nat_accSLR,
                ) = nourishment_stats(
                    cascade_acc_nat, cutoff=71, iB3D=0, nourishment_volume=100
                )

            def primary_stats():
                def thirds_acc_slr():
                    # roadway and nourishment statistics for 100 simulations -- thirds acc SLR
                    (
                        year_abandoned_middle_roadway,
                        sim_max_middle_roadway,
                        road_bulldozed_middle_roadway,
                        overwash_removed_middle_roadway,
                        dune_rebuilt_middle_roadway,
                        road_relocated_middle_roadway,
                        _,
                        _,
                    ) = get_roadway_statistics(
                        folder_prefix="9-CASCADE_AST_3domains_BE1m_AccSLR",
                        iB3D=3,  # the middle roadway section
                    )

                    drowned = np.array(sim_max_middle_roadway) < 199
                    percent_drown = sum(drowned)  # 81
                    min_drown = np.min(list(compress(sim_max_middle_roadway, drowned)))
                    max_drown = np.max(list(compress(sim_max_middle_roadway, drowned)))
                    mean_drown = np.mean(
                        list(compress(sim_max_middle_roadway, drowned))
                    )

                    abandoned = ~(
                        np.array(year_abandoned_middle_roadway)
                        == np.array(sim_max_middle_roadway) - 1
                    )  # if b3d drowned before the roadway was abandoned, the last roadway TS would be one behind; filter out
                    total_abandoned = sum(abandoned)
                    min_abandonment = np.min(
                        list(compress(year_abandoned_middle_roadway, abandoned))
                    )  # 46
                    max_abandonment = np.max(
                        list(compress(year_abandoned_middle_roadway, abandoned))
                    )  # 82
                    mean_abandonment = np.mean(
                        list(compress(year_abandoned_middle_roadway, abandoned))
                    )  # 57

                    (
                        year_abandoned_far_roadway,
                        sim_max_far_roadway,
                        road_bulldozed_far_roadway,
                        overwash_removed_far_roadway,
                        dune_rebuilt_far_roadway,
                        road_relocated_far_roadway,
                        _,
                        _,
                    ) = get_roadway_statistics(
                        folder_prefix="9-CASCADE_AST_3domains_BE1m_AccSLR",
                        iB3D=6,  # the far roadway section
                    )

                    drowned = np.array(sim_max_far_roadway) < 199
                    percent_drown = sum(drowned)  # 81
                    min_drown = np.min(list(compress(sim_max_far_roadway, drowned)))
                    max_drown = np.max(list(compress(sim_max_far_roadway, drowned)))
                    mean_drown = np.mean(list(compress(sim_max_far_roadway, drowned)))

                    abandoned = ~(
                        np.array(year_abandoned_far_roadway)
                        == np.array(sim_max_far_roadway) - 1
                    )  # if b3d drowned before the roadway was abandoned, the last roadway TS would be one behind; filter out
                    total_abandoned = sum(abandoned)
                    min_abandonment = np.min(
                        list(compress(year_abandoned_far_roadway, abandoned))
                    )
                    max_abandonment = np.max(
                        list(compress(year_abandoned_far_roadway, abandoned))
                    )
                    mean_abandonment = np.mean(
                        list(compress(year_abandoned_far_roadway, abandoned))
                    )

                    (
                        year_abandoned_community,
                        sim_max_community,
                        overwash_filtered_removed_community,
                        dune_rebuilt_community,
                        beach_nourished_community,
                        _,
                        _,
                    ) = get_nourishment_statistics(
                        folder_prefix="9-CASCADE_AST_3domains_BE1m_AccSLR",
                        iB3D=0,  # community
                    )

                    # NOTE: in some cases the last simulation time step will be different because of drowning of diff segments
                    drowned = np.array(sim_max_community) < 199
                    percent_drown = sum(drowned)  # 81
                    min_drown = np.min(list(compress(sim_max_community, drowned)))
                    max_drown = np.max(list(compress(sim_max_community, drowned)))
                    mean_drown = np.mean(list(compress(sim_max_community, drowned)))

                    abandoned = ~(
                        np.array(year_abandoned_community)
                        == np.array(sim_max_community) - 1
                    )  # if b3d drowned before the community was abandoned, the last community TS would be one behind; filter out
                    total_abandoned = sum(abandoned)  # 26
                    min_abandonment = np.min(
                        list(compress(year_abandoned_community, abandoned))
                    )  # 132
                    max_abandonment = np.max(
                        list(compress(year_abandoned_community, abandoned))
                    )  # 138
                    mean_abandonment = np.mean(
                        list(compress(year_abandoned_community, abandoned))
                    )  # 137

                    # roadway and nourishment statistics for 100 simulations -- nat acc SLR
                    (
                        year_abandoned_community,
                        sim_max_community,
                        overwash_filtered_removed_community,
                        dune_rebuilt_community,
                        beach_nourished_community,
                        _,
                        _,
                    ) = get_nourishment_statistics(
                        folder_prefix="9-CASCADE_AST_3domains_BE1m_AccSLR_nat_middle",
                        iB3D=0,  # community
                    )

                    # NOTE: in some cases the last simulation time step will be different because of drowning of diff segments
                    drowned = np.array(sim_max_community) < 199
                    percent_drown = sum(drowned)  # 54
                    min_drown = np.min(list(compress(sim_max_community, drowned)))  # 85
                    max_drown = np.max(
                        list(compress(sim_max_community, drowned))
                    )  # 153
                    mean_drown = np.mean(
                        list(compress(sim_max_community, drowned))
                    )  # 127

                    abandoned = ~(
                        np.array(year_abandoned_community)
                        == np.array(sim_max_community) - 1
                    )  # if b3d drowned before the community was abandoned, the last community TS would be one behind; filter out
                    total_abandoned = sum(abandoned)  #
                    min_abandonment = np.min(
                        list(compress(year_abandoned_community, abandoned))
                    )  # 132
                    max_abandonment = np.max(
                        list(compress(year_abandoned_community, abandoned))
                    )  # 138
                    mean_abandonment = np.mean(
                        list(compress(year_abandoned_community, abandoned))
                    )  # 137

                    (
                        year_abandoned_far_roadway,
                        sim_max_far_roadway,
                        road_bulldozed_far_roadway,
                        overwash_removed_far_roadway,
                        dune_rebuilt_far_roadway,
                        road_relocated_far_roadway,
                        _,
                        _,
                    ) = get_roadway_statistics(
                        folder_prefix="9-CASCADE_AST_3domains_BE1m_AccSLR_nat_middle",
                        iB3D=6,  # the far roadway section
                    )

                    drowned = np.array(sim_max_far_roadway) < 199
                    percent_drown = sum(drowned)  # 54
                    min_drown = np.min(list(compress(sim_max_far_roadway, drowned)))
                    max_drown = np.max(list(compress(sim_max_far_roadway, drowned)))
                    mean_drown = np.mean(list(compress(sim_max_far_roadway, drowned)))

                    abandoned = ~(
                        np.array(year_abandoned_far_roadway)
                        == np.array(sim_max_far_roadway) - 1
                    )  # if b3d drowned before the roadway was abandoned, the last roadway TS would be one behind; filter out
                    total_abandoned = sum(abandoned)
                    min_abandonment = np.min(
                        list(compress(year_abandoned_far_roadway, abandoned))
                    )  # 109
                    max_abandonment = np.max(
                        list(compress(year_abandoned_far_roadway, abandoned))
                    )  # 161
                    mean_abandonment = np.mean(
                        list(compress(year_abandoned_far_roadway, abandoned))
                    )

                def thirds_acc_slr_storminess():
                    # same as above, but increase in storminess
                    (
                        year_abandoned_middle_roadway,
                        sim_max_middle_roadway,
                        road_bulldozed_middle_roadway,
                        overwash_removed_middle_roadway,
                        dune_rebuilt_middle_roadway,
                        road_relocated_middle_roadway,
                        _,
                        _,
                    ) = get_roadway_statistics(
                        folder_prefix="9-CASCADE_AST_3domains_BE1m_AccSLR_AdaptationScenario",
                        iB3D=3,  # the middle roadway section
                    )

                    drowned = np.array(sim_max_middle_roadway) < 199
                    percent_drown = sum(drowned)  # 70
                    min_drown = np.min(
                        list(compress(sim_max_middle_roadway, drowned))
                    )  # 64
                    max_drown = np.max(
                        list(compress(sim_max_middle_roadway, drowned))
                    )  # 159
                    mean_drown = np.mean(
                        list(compress(sim_max_middle_roadway, drowned))
                    )  # 97

                    abandoned = ~(
                        np.array(year_abandoned_middle_roadway)
                        == np.array(sim_max_middle_roadway) - 1
                    )  # if b3d drowned before the roadway was abandoned, the last roadway TS would be one behind; filter out
                    total_abandoned = sum(abandoned)  # 59
                    min_abandonment = np.min(
                        list(compress(year_abandoned_middle_roadway, abandoned))
                    )  # 46
                    max_abandonment = np.max(
                        list(compress(year_abandoned_middle_roadway, abandoned))
                    )  # 141
                    mean_abandonment = np.mean(
                        list(compress(year_abandoned_middle_roadway, abandoned))
                    )  # 81

                    (
                        year_abandoned_far_roadway,
                        sim_max_far_roadway,
                        road_bulldozed_far_roadway,
                        overwash_removed_far_roadway,
                        dune_rebuilt_far_roadway,
                        road_relocated_far_roadway,
                        _,
                        _,
                    ) = get_roadway_statistics(
                        folder_prefix="9-CASCADE_AST_3domains_BE1m_AccSLR_AdaptationScenario",
                        iB3D=6,  # the far roadway section
                    )

                    drowned = np.array(sim_max_far_roadway) < 199
                    percent_drown = sum(drowned)  # 70
                    min_drown = np.min(
                        list(compress(sim_max_far_roadway, drowned))
                    )  # 64
                    max_drown = np.max(
                        list(compress(sim_max_far_roadway, drowned))
                    )  # 158
                    mean_drown = np.mean(
                        list(compress(sim_max_far_roadway, drowned))
                    )  # 97

                    abandoned = ~(
                        np.array(year_abandoned_far_roadway)
                        == np.array(sim_max_far_roadway) - 1
                    )  # if b3d drowned before the roadway was abandoned, the last roadway TS would be one behind; filter out
                    total_abandoned = sum(abandoned)  # 32
                    min_abandonment = np.min(
                        list(compress(year_abandoned_far_roadway, abandoned))
                    )  # 103
                    max_abandonment = np.max(
                        list(compress(year_abandoned_far_roadway, abandoned))
                    )  # 189
                    mean_abandonment = np.mean(
                        list(compress(year_abandoned_far_roadway, abandoned))
                    )  # 155

                    (
                        year_abandoned_community,
                        sim_max_community,
                        overwash_filtered_removed_community,
                        dune_rebuilt_community,
                        beach_nourished_community,
                        _,
                        _,
                    ) = get_nourishment_statistics(
                        folder_prefix="9-CASCADE_AST_3domains_BE1m_AccSLR_AdaptationScenario",
                        iB3D=0,  # community
                    )

                    # NOTE: in some cases the last simulation time step will be different because of drowning of diff segments
                    drowned = np.array(sim_max_community) < 199
                    percent_drown = sum(drowned)  # 70
                    min_drown = np.min(list(compress(sim_max_community, drowned)))  # 64
                    max_drown = np.max(
                        list(compress(sim_max_community, drowned))
                    )  # 159
                    mean_drown = np.mean(
                        list(compress(sim_max_community, drowned))
                    )  # 97

                    abandoned = ~(
                        np.array(year_abandoned_community)
                        == np.array(sim_max_community) - 1
                    )  # if b3d drowned before the community was abandoned, the last community TS would be one behind; filter out
                    total_abandoned = sum(abandoned)  # 41
                    min_abandonment = np.min(
                        list(compress(year_abandoned_community, abandoned))
                    )  # 127
                    max_abandonment = np.max(
                        list(compress(year_abandoned_community, abandoned))
                    )  # 138
                    mean_abandonment = np.mean(
                        list(compress(year_abandoned_community, abandoned))
                    )  # 136

                    # -------------------------
                    # roadway and nourishment statistics for 100 simulations -- nat acc SLR, increased storminess
                    (
                        year_abandoned_community,
                        sim_max_community,
                        overwash_filtered_removed_community,
                        dune_rebuilt_community,
                        beach_nourished_community,
                        _,
                        _,
                    ) = get_nourishment_statistics(
                        folder_prefix="9-CASCADE_AST_3domains_BE1m_AccSLR_nat_middle_AdaptationScenario",
                        iB3D=0,  # community
                    )

                    # NOTE: in some cases the last simulation time step will be different because of drowning of diff segments
                    drowned = np.array(sim_max_community) < 199
                    percent_drown = sum(drowned)  # 20
                    min_drown = np.min(
                        list(compress(sim_max_community, drowned))
                    )  # 136
                    max_drown = np.max(
                        list(compress(sim_max_community, drowned))
                    )  # 156
                    mean_drown = np.mean(
                        list(compress(sim_max_community, drowned))
                    )  # 144

                    abandoned = ~(
                        np.array(year_abandoned_community)
                        == np.array(sim_max_community) - 1
                    )  # if b3d drowned before the community was abandoned, the last community TS would be one behind; filter out
                    total_abandoned = sum(abandoned)  # 94
                    min_abandonment = np.min(
                        list(compress(year_abandoned_community, abandoned))
                    )  # 123
                    max_abandonment = np.max(
                        list(compress(year_abandoned_community, abandoned))
                    )  # 138
                    mean_abandonment = np.mean(
                        list(compress(year_abandoned_community, abandoned))
                    )  # 135

                    (
                        year_abandoned_far_roadway,
                        sim_max_far_roadway,
                        road_bulldozed_far_roadway,
                        overwash_removed_far_roadway,
                        dune_rebuilt_far_roadway,
                        road_relocated_far_roadway,
                        _,
                        _,
                    ) = get_roadway_statistics(
                        folder_prefix="9-CASCADE_AST_3domains_BE1m_AccSLR_nat_middle_AdaptationScenario",
                        iB3D=6,  # the far roadway section
                    )

                    drowned = np.array(sim_max_far_roadway) < 199
                    percent_drown = sum(drowned)  # 20
                    min_drown = np.min(
                        list(compress(sim_max_far_roadway, drowned))
                    )  # 135
                    max_drown = np.max(
                        list(compress(sim_max_far_roadway, drowned))
                    )  # 155
                    mean_drown = np.mean(
                        list(compress(sim_max_far_roadway, drowned))
                    )  # 144

                    abandoned = ~(
                        np.array(year_abandoned_far_roadway)
                        == np.array(sim_max_far_roadway) - 1
                    )  # if b3d drowned before the roadway was abandoned, the last roadway TS would be one behind; filter out
                    total_abandoned = sum(abandoned)  # 93
                    min_abandonment = np.min(
                        list(compress(year_abandoned_far_roadway, abandoned))
                    )  # 110
                    max_abandonment = np.max(
                        list(compress(year_abandoned_far_roadway, abandoned))
                    )  # 199
                    mean_abandonment = np.mean(
                        list(compress(year_abandoned_far_roadway, abandoned))
                    )  # 159
