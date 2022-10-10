# process model results -- including plotting and statistics -- for
#
# ~******* CASCADE ********~
#
# for the manuscript titled:
#
# "The Future of Developed Barrier Systems: Pathways Toward Uninhabitability, Drowning, and Rebound"
#
"""----------------------------------------------------
Copyright (C) 2022 Katherine Anarde
----------------------------------------------------"""

# please following the instructions at https://github.com/UNC-CECL/CASCADE for installing CASCADE

import numpy as np
import os

from cascade.tools import plotters as cascade_plt
from scripts.pathways_ms import plotters_pathways as pathways_plt

from itertools import compress

# # ###############################################################################
# # plotting functions
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
    ) = cascade_plt.plot_nonlinear_stats_RoadwayManager(
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
            cascade_plt.plot_ElevAnimation_CASCADE(
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
            cascade_plt.plot_ElevAnimation_CASCADE(
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
        fig = cascade_plt.plot_ModelTransects(cascade, time_step, iB3D=0)
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
    ) = cascade_plt.plot_nonlinear_stats_BeachDuneManager(
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
        cascade_plt.plot_ElevAnimation_CASCADE(
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
        fig = cascade_plt.plot_ModelTransects(cascade, time_step, iB3D=0)
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

    pathways_plt.fig3_initialCNH_topo(cascade)

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
                ) = cascade_plt.plot_nonlinear_stats_BeachDuneManager(
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
                ) = cascade_plt.plot_nonlinear_stats_RoadwayManager(
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
                ) = cascade_plt.plot_nonlinear_stats_RoadwayManager(
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
        cascade_plt.plot_ElevAnimation_CASCADE(
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
# # human dynamics statistics
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
# # record of plots and statistics
# # ###############################################################################

# 10,000 year plots -------------------------------------------------------
def cascade_10kyr_plots():

    datadir = "/Users/KatherineAnardeWheels/Research/BARis/UNC/CNH/CASCADE_save_dir/Run_Output/"
    tmax_pt45 = [10000, 10000, 10000, 10000, 10000]
    name_prefix_45 = "4-B3D_noAST_Rave_pt45_SLR_pt004_10k-yrs_0"
    tmax_pt75 = [5710, 10000, 2878, 10000, 10000]
    # tmax_pt75_old = [5725, 992, 4870, 10000, 6669]
    name_prefix_75 = "4-B3D_noAST_Rave_pt75_SLR_pt004_10k-yrs_0"

    pathways_plt.fig2_10kyr_timeseries(
        datadir, tmax_pt45, name_prefix_45, vertical_line_1=8757, vertical_line_2=802
    )
    pathways_plt.fig2_10kyr_timeseries(
        datadir, tmax_pt75, name_prefix_75, vertical_line_1=4261, vertical_line_2=829
    )

# 1,000 year plots -------------------------------------------------------------------
# NOTE: these are organized into functions, but don't actually act as functions -- I just wanted a way to collapse
# and organize my analysis
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

            pathways_plt.fig5_8_plot_human_dynamics_stats_array4(
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

            pathways_plt.fig5_8_plot_human_dynamics_stats_array4(
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

            pathways_plt.fig5_8_plot_human_dynamics_stats_array4(
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

            pathways_plt.fig5_8_plot_human_dynamics_stats_array4(
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
                name_prefix="4-B3D_Rave_pt45_Natural_low",
                tmax_roadways=1000,  # dummy
                tmax_sim=1000,
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
                name_prefix="6-B3D_Rave_pt45_Roadways_2mDune_20mSetback_20mWidth_low_10percent",
                tmax_roadways=462,
                tmax_sim=700,
                plot_name="b3d_pt45_h2m_plots_low_10percent",
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
                name_prefix="6-B3D_Rave_pt45_Roadways_2mDune_20mSetback_20mWidth_low_20percent",
                tmax_roadways=533,
                tmax_sim=700,
                plot_name="b3d_pt45_h2m_plots_low_20percent",
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
                name_prefix="6-B3D_Rave_pt45_Roadways_2mDune_20mSetback_20mWidth_low_30percent",
                tmax_roadways=545,
                tmax_sim=700,
                plot_name="b3d_pt45_h2m_plots_low_30percent",
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
                name_prefix="6-B3D_Rave_pt45_Roadways_2mDune_20mSetback_20mWidth_low_40percent",
                tmax_roadways=548,
                tmax_sim=700,
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
                name_prefix="6-B3D_Rave_pt45_Roadways_2mDune_20mSetback_20mWidth_low_50percent",
                tmax_roadways=553,
                tmax_sim=700,
                plot_name="b3d_pt45_h2m_plots_low_50percent",
                run_road_mgmt=True,
            )

            pathways_plt.supp_sensitivity_road_abandonment(
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
                TMAX=[600, 600, 600, 600, 600],
                tmax_roadways=[462, 533, 545, 548, 553],
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
            pathways_plt.fig4_slr_sensitivity(
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
            pathways_plt.fig4_slr_sensitivity(
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
            pathways_plt.fig4_slr_sensitivity(
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
            pathways_plt.fig4_slr_sensitivity(
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
                    pathways_plt.fig5_8_plot_human_dynamics_stats_array4(
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
                    pathways_plt.fig5_8_plot_human_dynamics_stats_array4(
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
                    pathways_plt.fig5_8_plot_human_dynamics_stats_array4(
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
                    pathways_plt.fig5_8_plot_human_dynamics_stats_array4(
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

        # thresholds_supplement
        pathways_plt.supp_nourishment_thresholds(
            directory="/Users/KatherineAnardeWheels/PycharmProjects/CASCADE/scripts/pathways_ms/data/Nags_Head"
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
        # All 3: Island is too narrow for roadway to be relocated. Roadway eaten up by dunes at 143
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
            # Island is too narrow for roadway to be relocated. Roadway eaten up by dunes at 73
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

        def ast_time_series(
                shoreline_position_allroads_pt75low,
                beach_width_allnourish_pt75low,
                beach_width_allnourish_pt45low
        ):
            iB3D_roadways = 4
            iB3D_community = 1

            # pt75
            pathways_plt.fig11_14_stats_ast_array3(
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
            pathways_plt.fig11_14_stats_ast_array3(
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

        def ast_time_series(
                shoreline_position_allroads_pt75low,
                shoreline_position_thirds_acc,
                shoreline_position_acc_nat,
                beach_width_allnourish_pt45low,
                beach_width_thirds_acc,
                beach_width_acc_nat
        ):
            iB3D_roadways = 4
            iB3D_community = 1

            # pt45 and pt75 low
            pathways_plt.fig11_14_stats_ast_array3(
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

        def misc_stats(
                cascade_pt45nourish_linSLR,
                cascade_thirds_linSLR,
                cascade_thirds_acc,
                cascade_acc_nat
        ):
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
