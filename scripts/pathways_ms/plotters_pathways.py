# plotting functions for
#
# ~******* CASCADE ********~
#
# and for the manuscript titled:
#
# "The Future of Developed Barrier Systems: Pathways Toward Uninhabitability, Drowning, and Rebound"
#
"""----------------------------------------------------
Copyright (C) 2022 Katherine Anarde
----------------------------------------------------"""

import glob
import math
import os

import imageio
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn
from matplotlib import cm
from matplotlib.ticker import AutoMinorLocator

# # ###############################################################################
# # plotters for ms
# # ###############################################################################


def fig5_8_plot_human_dynamics_stats_array4(
    cascade,  # these are lists
    DuneCrestMin,
    DuneCrestMax,
    BarrierHeight,
    BarrierWidth,
    DuneCrestMean,
    shoreline_position,
    overwash,
    TMAX,
    tmax_management,
    dune_toe=None,
    roadways_on=True,
    nourishment_on=False,
    rebuild_threshold=None,
    scenarios=None,
):
    iB3D = 0

    # dunes ---------------------------------------------------------------------------------------------------------- #
    fig, axs = plt.subplots(1, 4, figsize=(10, 3), sharey=True, sharex=True)
    for i in range(len(cascade)):
        if i > 0:
            time = np.arange(0, TMAX[i] - 0.5, 0.5)
            yearly_time = np.arange(0, TMAX[i], 1)
            dune_min = DuneCrestMin[i][0 : len(time)]
            dune_max = DuneCrestMax[i][0 : len(time)]
            mask1 = np.isfinite(dune_min)
            mask2 = np.isfinite(dune_max)
            dunes_rebuilt_cellular_shoreline_change = (
                np.array(
                    cascade[i]
                    .barrier3d[iB3D]
                    ._ShorelineChangeTS[0 : tmax_management[i]]
                )
                < 0
            )  # so really this is just shoreline change; I have to go in and manually delete occurences where it does
            # not coincide with the dune height min (and double check with dune_rebuild_TS)
            indices = [
                i for i, x in enumerate(dunes_rebuilt_cellular_shoreline_change) if x
            ]  # for nourishment this is typically zero
            axs[i].plot(time[mask1], dune_min[mask1])
            axs[i].plot(time[mask2], dune_max[mask2])

            if roadways_on:
                axs[i].plot(
                    yearly_time[0 : tmax_management[i]],
                    cascade[i]
                    .roadways[iB3D]
                    ._dune_minimum_elevation_TS[0 : tmax_management[i]],
                    "-kD",
                    markevery=indices,
                )
                axs[i].plot(
                    yearly_time[0 : tmax_management[i]],
                    cascade[i]
                    .roadways[iB3D]
                    ._dune_design_elevation_TS[0 : tmax_management[i]],
                )
                axs[i].plot(
                    yearly_time[0 : tmax_management[i]],
                    cascade[i].roadways[iB3D]._road_ele_TS[0 : tmax_management[i]],
                )
            if nourishment_on:
                axs[i].hlines(
                    rebuild_threshold,  # in m
                    time[0],
                    # time[-1],
                    tmax_management[i],
                    "red",
                )
                axs[i].hlines(
                    cascade[i].nourishments[iB3D]._dune_design_elevation,  # in m
                    time[0],
                    # time[-1],
                    tmax_management[i],
                    "yellow",
                )

            axs[i].hlines(
                cascade[i].barrier3d[0]._Dmaxel * 10, time[0], time[-1], colors="green"
            )
        else:
            # the natural scenario
            time = np.arange(0, TMAX[i], 1)
            axs[i].plot(time, DuneCrestMin[i][0 : TMAX[i]])
            axs[i].plot(time, DuneCrestMax[i][0 : TMAX[i]])
            axs[i].hlines(
                cascade[i].barrier3d[0]._Dmaxel * 10, time[0], time[-1], colors="green"
            )
        axs[i].set(xlabel="time (yr)")
        axs[i].xaxis.set_minor_locator(AutoMinorLocator())
        # axs[i].tick_params(which='minor', length=4, color='r')
        if roadways_on:
            # axs[i].set_xlim([-15, 765])
            axs[i].set_xlim([-15, 715])
            # axs[i].set_xlim([-15, 1015])
            axs[i].set_ylim([0, 5.25])
        if nourishment_on:
            # axs[i].set_xlim([-15, 1015])
            # axs[i].set_xlim([-15, 765])
            axs[i].set_xlim([-15, 715])
            axs[i].set_ylim([0, 5.25])
    axs[0].set(ylabel="elevation (m MHW)")
    if roadways_on:
        axs[3].legend(
            [
                "dune alongshore min",
                "dune alongshore max",
                "dune rebuild threshold",
                "dune design height",
                "road",
                "dune max-equilibrium",
            ]
        )
    if nourishment_on:
        axs[3].legend(
            [
                "dune along. min",
                "dune along. max",
                "dune rebuild threshold",
                "dune design height",
                "dune max-equil",
            ]
        )
    # plt.tight_layout()

    # interior height, width, and overwash --------------------------------------------------------------------------- #
    fig2, axs2 = plt.subplots(1, 4, figsize=(10, 3), sharex=True)
    color = ["b", "r", "g", "m"]
    for i in range(len(cascade)):
        if i > 0:
            # time = np.arange(0, TMAX[i] - 0.5, 0.5)
            yearly_time = np.arange(0, TMAX[i], 1)

            # NOTE: I decided I didn't want to plot the 0.5 year time step (post storm), just the final after mgmt
            # barrier_height = BarrierHeight[i][0 : len(time)]
            # mask_bh = np.isfinite(barrier_height)
            barrier_height = np.array(cascade[i].barrier3d[0].h_b_TS[0 : TMAX[i]]) * 10

            scts = [
                (x - shoreline_position[i][0])
                for x in shoreline_position[i][0 : TMAX[i]]
            ]
            axs2[2].plot(yearly_time, scts, color[i])

            if roadways_on:
                axs2[0].plot(yearly_time, barrier_height, color[i])
                axs2[1].plot(yearly_time, BarrierWidth[i][0 : TMAX[i]], color[i])
                axs2[3].plot(
                    yearly_time, overwash[i][0 : TMAX[i]], color[i]
                )  # this is just net overwash

            if nourishment_on:
                # plot interannual time steps for barrier width and height
                # axs2[0].plot(
                #     time[mask_bh], barrier_height[mask_bh], color[i], alpha=0.3
                # )
                axs2[0].plot(
                    yearly_time,
                    (np.array(cascade[i].barrier3d[0].h_b_TS[0 : TMAX[i]]) * 10),
                    color[i],
                )

                # barrier_width = BarrierWidth[i][0 : len(time)]
                # mask_bw = np.isfinite(barrier_width)
                # axs2[1].plot(time[mask_bw], barrier_width[mask_bw], color[i], alpha=0.3)
                axs2[1].plot(
                    yearly_time,
                    (
                        np.array(
                            cascade[i].barrier3d[0].InteriorWidth_AvgTS[0 : TMAX[i]]
                        )
                        * 10
                    ),  # m
                    color[i],
                )

                # dtts = np.array(
                #     [(x - shoreline_position[i][0]) for x in dune_toe[i][0 : TMAX[i]]]
                # )
                # mask_dtts = np.isnan(dtts)
                # axs2[2].plot(
                #     yearly_time[~mask_dtts], dtts[~mask_dtts], color[i], alpha=0.3
                # )

                # don't plot interannual overwash, too messy
                # ow = overwash[i][0 : len(time)]
                # mask_ow = np.isfinite(ow)
                # axs2[3].plot(time[mask_ow], ow[mask_ow], color[i], alpha=0.3)
                axs2[3].plot(
                    yearly_time, cascade[i].barrier3d[0].QowTS[0 : TMAX[i]], color[i]
                )

            axs2[0].axvline(x=tmax_management[i], color=color[i])
            axs2[1].axvline(x=tmax_management[i], color=color[i])
            axs2[2].axvline(x=tmax_management[i], color=color[i])
            axs2[3].axvline(x=tmax_management[i], color=color[i])

        else:
            # the natural scenario
            yearly_time = np.arange(0, TMAX[i], 1)
            axs2[0].plot(yearly_time, BarrierHeight[i][0 : TMAX[i]], color[i])
            axs2[1].plot(yearly_time, BarrierWidth[i][0 : TMAX[i]], color[i])
            scts = [
                (x - cascade[i].barrier3d[0]._x_s_TS[0]) * 10
                for x in cascade[i].barrier3d[0]._x_s_TS[0 : TMAX[i]]
            ]
            axs2[2].plot(yearly_time, scts, color[i])
            axs2[3].plot(
                yearly_time, cascade[i].barrier3d[0].QowTS[0 : TMAX[i]], color[i]
            )

        axs2[i].set(xlabel="time (yr)")
        axs2[i].xaxis.set_minor_locator(AutoMinorLocator())
        # axs2[i].tick_params(which='minor', length=4, color='k')

        if roadways_on:
            # axs2[i].set_xlim([-15, 765])
            axs2[i].set_xlim([-15, 1015])
        if nourishment_on:
            axs2[i].set_xlim([-15, 1015])
            # axs2[i].set_xlim([-15, 765])

    axs2[0].set(ylabel="barrier elevation (m MHW)")
    axs2[0].set_ylim([-0.03, 1.75])
    axs2[1].set(ylabel="barrier width (m)")
    axs2[1].set_ylim([-6, 450])
    if roadways_on:
        axs2[2].set_ylim([-10, 900])
        axs2[2].set(ylabel="shoreline position (m)")
    if nourishment_on:
        axs2[2].set_ylim([-30, 900])
        axs2[2].set(ylabel="shoreline position (m)")
    axs2[3].set(ylabel="overwash flux (m$^3$/m)")
    axs2[3].set_ylim([-3, 225])
    if roadways_on:
        scenarios = [
            "natural",
            "1-m dune",
            "2-m dune",
            "3-m dune",
            "end mgmt",
            "drowned",
        ]
        axs2[2].legend(scenarios)
    if nourishment_on:
        axs2[2].legend(scenarios)
    plt.tight_layout()

    # plot a zoom in of nourishments -- this is the only place we plot the dune toe
    if nourishment_on:
        fig3, axs3 = plt.subplots(1, 2, figsize=(10, 3), sharex=True)

        color = ["b", "r", "g", "m"]
        for i in range(len(cascade)):
            if i > 0:
                time = np.arange(0, TMAX[i] - 0.5, 0.5)
                yearly_time = np.arange(0, TMAX[i], 1)

                scts = [
                    (x - shoreline_position[i][0])
                    for x in shoreline_position[i][0 : TMAX[i]]
                ]
                axs3[1].plot(yearly_time, scts, color[i])

                dtts = np.array(
                    [(x - shoreline_position[i][0]) for x in dune_toe[i][0 : TMAX[i]]]
                )
                mask_dtts = np.isnan(dtts)
                axs3[1].plot(
                    yearly_time[~mask_dtts], dtts[~mask_dtts], color[i], alpha=0.3
                )

            else:
                # the natural scenario
                yearly_time = np.arange(0, TMAX[i], 1)
                scts = [
                    (x - cascade[i].barrier3d[0]._x_s_TS[0]) * 10
                    for x in cascade[i].barrier3d[0]._x_s_TS[0 : TMAX[i]]
                ]
                axs3[1].plot(yearly_time, scts, color[i])

            axs3[1].set(ylabel="cross-shore position (m)")
            axs3[1].set_ylim([-30, 50])
            axs3[1].set(xlabel="time (yr)")
            axs3[1].set_xlim([110, 210])
            axs3[1].legend(["shoreline", "dune toe"])

    # # characteristic trajectory plots
    # fig3, axs3 = plt.subplots(1, 4, figsize=(10, 3), sharex=True)
    # cmap = ["Blues", "Oranges", "Greens", "Purples"]
    # for i in range(len(cascade)):
    #
    #     if i > 0:
    #         axs3[0].scatter(
    #             DuneCrestMean[i],
    #             BarrierHeight[i],
    #             c=np.arange(0, np.size(DuneCrestMean[i]), 1),
    #             cmap=cmap[i],
    #         )
    #         axs3[1].scatter(
    #             DuneCrestMean[i],
    #             BarrierWidth[i],
    #             c=np.arange(0, np.size(DuneCrestMean[i]), 1),
    #             cmap=cmap[i],
    #         )
    #
    #     else:
    #         axs3[0].scatter(
    #             DuneCrestMean[i][0 : TMAX[i]],
    #             BarrierHeight[i][0 : TMAX[i]],
    #             c=np.arange(0, np.size(DuneCrestMean[i][0 : TMAX[i]]), 1),
    #             cmap=cmap[i],
    #         )
    #         axs3[1].scatter(
    #             DuneCrestMean[i][0 : TMAX[i]],
    #             BarrierWidth[i][0 : TMAX[i]],
    #             c=np.arange(0, np.size(DuneCrestMean[i][0 : TMAX[i]]), 1),
    #             cmap=cmap[i],
    #         )
    #     axs3[i].set(xlabel="mean dune ele (m MHW)")
    # axs3[0].set(ylabel="barrier height (m MHW)")
    # axs3[1].set(ylabel="barrier width (m)")
    # # axs3[1].legend(["natural", "1-m dune", "2-m dune", "3-m dune"])
    # plt.tight_layout()


def fig11_14_stats_ast_array3(
    shoreline_position,
    beach_width,
    TMAX,
    tmax_management_nourishments,
    tmax_management_roadways,
    scenarios_beach_width,
    scenarios_shoreline_position,
):
    fig, axs = plt.subplots(1, 2, figsize=(10, 3))
    color = ["b", "r", "m"]

    for i in range(len(TMAX)):
        time = np.arange(0, TMAX[i] - 0.5, 0.5)
        yearly_time = np.arange(0, TMAX[i], 1)

        # includes post-storm beach width, before human modifications
        mask = np.isfinite(beach_width[i])
        axs[0].plot(time[mask], beach_width[i][mask], color[i])

        # shoreline position relative to 0 starting position
        scts = [
            (x - shoreline_position[i][0]) for x in shoreline_position[i][0 : TMAX[i]]
        ]
        axs[1].plot(yearly_time, scts, color[i])

    axs[0].set(ylabel="beach width (m)")
    axs[0].set_ylim([8, 52])
    axs[1].set(ylabel="shoreline position (m)")
    # axs2[1].set_ylim([-6, 425])
    axs[0].legend(scenarios_beach_width)
    axs[1].legend(scenarios_shoreline_position)

    # plot when both roadway and nourishment management ceased
    for i in range(len(TMAX)):
        axs[0].axvline(x=tmax_management_roadways[i], color=color[i])
        axs[0].axvline(x=tmax_management_nourishments[i], color=color[i])
        axs[1].axvline(x=tmax_management_roadways[i], color=color[i])
        axs[1].axvline(x=tmax_management_nourishments[i], color=color[i])

    axs[0].set(xlabel="time (yr)")
    # axs[0].set_xlim([0, np.max(tmax_management_nourishments) + 2])
    axs[0].set_xlim([0, np.max(TMAX[:])])
    axs[1].set_xlim([0, np.max(TMAX[:])])
    axs[1].set(xlabel="time (yr)")

    plt.tight_layout()


def supp_sensitivity_road_abandonment(
    cascade,  # these are lists
    BarrierHeight,
    BarrierWidth,
    TMAX,
    tmax_roadways,
):
    # interior height and width
    fig2, axs2 = plt.subplots(1, 3, figsize=(10, 3), sharex=True)
    color = mpl.cm.inferno_r(np.linspace(0, 1, 5))
    for i in range(len(cascade)):
        time = np.arange(0, TMAX[i] - 0.5, 0.5)
        yearly_time = np.arange(0, TMAX[i], 1)
        barrier_height = BarrierHeight[i][0 : len(time)]
        mask = np.isfinite(barrier_height)
        axs2[0].plot(time[mask], barrier_height[mask], color=color[i])
        axs2[1].plot(yearly_time, BarrierWidth[i][0 : TMAX[i]], color=color[i])
        scts = [
            (x - cascade[i].barrier3d[0]._x_s_TS[0]) * 10
            for x in cascade[i].barrier3d[0]._x_s_TS[0 : TMAX[i]]
        ]
        axs2[2].plot(yearly_time, scts, color=color[i])
        axs2[0].axvline(x=tmax_roadways[i], color=color[i])
        axs2[1].axvline(x=tmax_roadways[i], color=color[i])
        axs2[2].axvline(x=tmax_roadways[i], color=color[i])

    axs2[0].set(ylabel="barrier elevation (m MHW)")
    axs2[0].set(xlabel="time (yr)")
    axs2[1].set(ylabel="barrier width (m)")
    axs2[1].set(xlabel="time (yr)")
    axs2[2].set(ylabel="shoreline position (m)")
    axs2[2].set(xlabel="time (yr)")
    axs2[2].legend(["10%", "20%", "30%", "40%", "50%"])
    plt.tight_layout()


def supp_nourishment_thresholds(
    directory,  # directory containing csv files
):
    # directory = "/Users/KatherineAnardeWheels/PycharmProjects/CASCADE/scripts/pathways_ms/data/Nags_Head"

    # use glob to get all the csv files
    csv_files = glob.glob(os.path.join(directory, "*.csv"))
    beach_width = []
    dune_height = []
    labels = [
        1997,
        1998,
        1999,
        2000,
        2001,
        2004,
        2005,
        2008,
        2009,
        2012,
        2013,
        2014,
        2015,
        2016,
        2017,
        2018,
        2019,
    ]

    # loop over the list of csv files
    for f in csv_files:
        # read the csv file
        df = pd.read_csv(f)
        bw = df["Beach Width"]
        filter = bw < 9000  # there seem to be some erroneus values
        beach_width.append(bw[filter])
        dune_height.append(df["Dune Height"])

    # plot beach width
    fig1 = plt.subplots(1, 1, figsize=(10, 3))
    ax = seaborn.boxplot(data=beach_width, flierprops={"marker": "o"})
    ax.set_xticklabels(labels)
    ax.set(ylabel="beach width (m)")
    ax.set(xlabel="time (yr)")

    # plot dune height
    fig2 = plt.subplots(1, 1, figsize=(10, 3))
    ax = seaborn.boxplot(data=dune_height, flierprops={"marker": "o"})
    ax.set_xticklabels(labels)
    ax.set(ylabel="dune height (m)")
    ax.set(xlabel="time (yr)")

    return fig1, fig2


def fig3_initialCNH_topo(
    cascade_model_list,  # must be from a nourishment simulation (i.e., have a beach width)
    km_on=True,
):
    fig, axs = plt.subplots(1, 4, figsize=(10, 3), sharey=True, sharex=True)

    # make the mat image of the beach to the back-barrier; stole this code from the animation plots above
    for iCascade in range(len(cascade_model_list)):
        cascade = cascade_model_list[iCascade]
        barrier3d = cascade.barrier3d

        BarrierLength = barrier3d[0]._BarrierLength
        OriginY = 52
        AniDomainWidth = 120  # end Y
        ny = 1
        iB3D = 0
        AnimateDomain = np.ones([AniDomainWidth + 1, BarrierLength * ny]) * -1
        scts = 0  # make shoreline position relative to zero starting

        # Build beach elevation domain, we only show beach width decreasing in increments of 10 m and we don't
        # illustrate a berm, just a sloping beach up to the elevation of the berm
        BeachWidth = math.floor(
            cascade.nourishments[iB3D].beach_width[0] / 10
        )  # switched from ceil
        BeachDomain = np.zeros(
            [
                BeachWidth,
                BarrierLength,
            ]
        )
        # beach width must be greater than zero
        add = (barrier3d[iB3D].BermEl - barrier3d[iB3D]._SL) / (BeachWidth + 1)
        for i in range(0, BeachWidth):
            BeachDomain[i, :] = (barrier3d[iB3D]._SL + add) * (i + 1)

        # Make animation frame domain
        Domain = barrier3d[iB3D]._DomainTS[0] * 10  # m MHW
        Dunes = (
            barrier3d[iB3D]._DuneDomain[0, :, :] + barrier3d[iB3D]._BermEl
        ) * 10  # m MHW
        Dunes = np.rot90(Dunes)
        Dunes = np.flipud(Dunes)
        Beach = BeachDomain * 10
        Domain = np.vstack([Beach, Dunes, Domain])
        Domain[Domain < 0] = -1  # anything underwater
        widthTS = len(Domain)
        OriginTstart = OriginY + math.floor(scts)  # ceil
        OriginTstop = OriginTstart + widthTS
        xOrigin = iB3D * BarrierLength
        AnimateDomain[OriginTstart:OriginTstop, xOrigin : xOrigin + BarrierLength] = (
            Domain
        )

        # plot
        print(np.max(AnimateDomain))
        cax = axs[iCascade].matshow(
            AnimateDomain, origin="lower", cmap="terrain", vmin=-1, vmax=3.0
        )  # , interpolation='gaussian') # analysis:ignore
        axs[iCascade].xaxis.set_ticks_position("bottom")
        axs[iCascade].set(xlabel="alongshore distance (dam)")
        if km_on:
            axs[iCascade].set(xlabel="alongshore distance (km)")
        axs[iCascade].set_ylim([40, 110])
        axs[iCascade].set_xlim([-1, 50])

    axs[0].set(ylabel="cross-shore distance (dam)")
    # cbar = fig.colorbar(cax)
    # cbar.set_label("elevation (m MHW)", rotation=270)
    # plt.tight_layout()
    if km_on:
        locs, _ = plt.yticks()
        plt.yticks(locs, locs / 100)
        locs, _ = plt.xticks()
        plt.xticks(locs[1:], locs[1:] / 100)
        axs[0].set(ylabel="cross-shore distance (km)")

    # now make the cross-section; stole this from the cross-section code above and modified
    v = 10  # just use the 10th transect
    fig, axs = plt.subplots(1, 2, figsize=(10, 3), sharey=True, sharex=True)

    for iCascade in range(len(cascade_model_list)):
        cascade = cascade_model_list[iCascade]

        sea_level = cascade.barrier3d[iB3D]._SL

        # Create data points
        shoreface_toe_x = (
            cascade.barrier3d[iB3D].x_t_TS[0] - cascade.barrier3d[iB3D].x_t_TS[0]
        )
        shoreface_toe_y = (sea_level - cascade.barrier3d[iB3D].DShoreface) * 10  # m
        shoreline_x = (
            cascade.barrier3d[iB3D].x_s_TS[0] - cascade.barrier3d[iB3D].x_t_TS[0]
        )
        shoreline_y = sea_level * 10  # m
        bay_y = (sea_level - cascade.barrier3d[iB3D]._BayDepth) * 10  # m
        end_of_bay_y = bay_y

        berm_x = shoreline_x + (
            cascade.nourishments[iB3D].beach_width[0] / 10
        )  # beach width (in dam)
        berm_y = (
            cascade.barrier3d[iB3D]._BermEl * 10
        ) + shoreline_y  # convert to meters
        dune_toe_x = berm_x
        dune_toe_y = berm_y

        interior_y = cascade.barrier3d[iB3D]._DomainTS[0]
        interior_y = interior_y[:, v]
        print(
            np.max(cascade.barrier3d[iB3D]._DuneDomain[0, v, :] * 10)
        )  # max dune height
        dunes_y = (
            cascade.barrier3d[iB3D]._DuneDomain[0, v, :]
            + cascade.barrier3d[iB3D]._BermEl
        )
        cross_barrier_y = np.insert(interior_y, 0, dunes_y)
        cross_barrier_y = (cross_barrier_y * 10) + shoreline_y  # Convert to meters
        cross_barrier_x = np.arange(0, len(cross_barrier_y), 1) + dune_toe_x

        end_of_bay_x = (
            cross_barrier_x[-1] + 20
        )  # just add a buffer to the end of the plt

        x = np.hstack(
            [
                shoreface_toe_x,
                shoreline_x,
                berm_x,
                dune_toe_x,
                cross_barrier_x,
                end_of_bay_x,
            ]
        )
        y = np.hstack(
            [
                shoreface_toe_y,
                shoreline_y,
                berm_y,
                dune_toe_y,
                cross_barrier_y,
                end_of_bay_y,
            ]
        )

        # Plot
        if iCascade < 2:  # just the 0.45 cases
            axs[0].plot(x, y)
            axs[0].hlines(
                sea_level * 10, shoreface_toe_x, end_of_bay_x, colors="dodgerblue"
            )
            # axs[0].set_xlim([0, 110])
            axs[0].set(xlabel="cross-shore distance (dam)")
            axs[0].set(ylabel="elevation (m MHW)")
        else:  # 0.75 cases
            axs[1].plot(x, y)
            axs[1].hlines(
                sea_level * 10, shoreface_toe_x, end_of_bay_x, colors="dodgerblue"
            )
            # axs[1].set_xlim([0, 110])
            axs[1].set(xlabel="cross-shore distance (dam)")
    # plt.tight_layout()

    if km_on:
        locs, _ = plt.xticks()
        plt.xticks(locs[1:], locs[1:] / 100)
        axs[0].set(xlabel="cross-shore distance (km)")
        axs[1].set(xlabel="cross-shore distance (km)")
        axs[0].set_xlim([-1, 141])
        axs[1].set_xlim([-1, 141])
        axs[0].legend(["profile A", "profile B"])
        axs[1].legend(["profile C", "profile D"])


def fig4_slr_sensitivity(
    cascade,  # lists
    TMAX,
):
    # col 1 - barrier height (4 different SLR)
    # col 2 - barrier width (4 different SLR)
    # col 3 - shoreline position (4 different SLR)
    # col 4 - overwash flux (4 diff SLR)

    fig1, axs1 = plt.subplots(1, 4, figsize=(10, 3), sharex=True)
    for i in range(len(cascade)):
        time = np.arange(0, TMAX[i], 1)

        RSLR = cascade[i].barrier3d[0].RSLR[0 : TMAX[i]]  # in dam
        SLR = np.cumsum(np.array(RSLR) * 10)

        axs1[0].plot(time, SLR)

    axs1[0].set(ylabel="sea level (m)")
    axs1[0].set(xlabel="time (yr)")
    axs1[0].legend(["4 mm/yr", "8 mm/yr", "12 mm/yr", "Acc"])
    plt.tight_layout()

    fig2, axs2 = plt.subplots(1, 4, figsize=(10, 3), sharex=True)
    for i in range(len(cascade)):
        time = np.arange(0, TMAX[i], 1)

        barrier_width = (
            np.array(cascade[i].barrier3d[0].x_b_TS[0 : TMAX[i]])
            - np.array(cascade[i].barrier3d[0].x_s_TS[0 : TMAX[i]])
        ) * 10

        barrier_height = np.array(cascade[i].barrier3d[0].h_b_TS[0 : TMAX[i]]) * 10

        scts = [
            (x - cascade[i].barrier3d[0].x_s_TS[0]) * 10
            for x in cascade[i].barrier3d[0].x_s_TS[0 : TMAX[i]]
        ]

        overwash_flux = cascade[i].barrier3d[0].QowTS[0 : TMAX[i]]

        axs2[0].plot(time, barrier_height)
        axs2[1].plot(time, barrier_width)
        axs2[2].plot(time, scts)
        axs2[3].plot(time, overwash_flux)

        axs2[i].set(xlabel="time (yr)")

    axs2[0].set(ylabel="barrier elevation (m MHW)")
    axs2[0].set_ylim([-0.03, 1.75])
    axs2[1].set(ylabel="barrier width (m)")
    axs2[1].set_ylim([-6, 450])
    axs2[2].set(ylabel="shoreline position (m)")
    axs2[2].set_ylim([-10, 900])
    axs2[2].legend(["4 mm/yr", "8 mm/yr", "12 mm/yr", "Acc"])
    axs2[3].set(ylabel="overwash flux (m$^3$/m)")
    axs2[3].set_ylim([-3, 225])
    plt.tight_layout()


def fig2_10kyr_timeseries(datadir, tmax, name_prefix, vertical_line_1, vertical_line_2):
    b3d = []
    bw = []
    bh = []
    dune_height = []
    time = []

    ib3d = 0

    for i in range(1, 6):
        full_name_prefix = name_prefix + str(i)
        output = np.load(datadir + full_name_prefix + ".npz", allow_pickle=True)
        cascade = output["cascade"][0]
        b3d.append(cascade.barrier3d)

        BarrierWidth = (
            np.array(cascade.barrier3d[ib3d].InteriorWidth_AvgTS[0 : tmax[i - 1]])
        ) * 10  # m
        bw.append(BarrierWidth)

        BarrierHeight = (
            np.array(cascade.barrier3d[ib3d].h_b_TS[0 : tmax[i - 1]])
        ) * 10  # m
        bh.append(BarrierHeight)

        time.append(np.arange(0, tmax[i - 1], 1))

        dune_height.append(
            np.array(cascade.barrier3d[ib3d]._Hd_AverageTS[0 : tmax[i - 1]]) * 10
        )

    fig1, axs1 = plt.subplots(3, 1, figsize=(5, 7), sharex="col")

    tlow = 0
    thigh = 10000

    axs1[0].plot(
        time[0], bw[0], time[1], bw[1], time[2], bw[2], time[3], bw[3], time[4], bw[4]
    )
    axs1[0].set(ylabel="barrier width (m)")
    axs1[0].set_xlim([tlow, thigh])
    axs1[0].set_ylim([-10, 480])
    axs1[0].axvline(x=vertical_line_1)
    axs1[0].axvline(x=vertical_line_2)

    axs1[1].plot(
        time[0], bh[0], time[1], bh[1], time[2], bh[2], time[3], bh[3], time[4], bh[4]
    )
    axs1[1].set(ylabel="barrier elevation (m MHW)")
    axs1[1].set_xlim([tlow, thigh])
    axs1[1].set_ylim([-0.10, 1.8])
    axs1[1].axvline(x=vertical_line_1)
    axs1[1].axvline(x=vertical_line_2)
    axs1[1].legend(
        [
            "storm sequence 1",
            "storm sequence 2",
            "storm sequence 3",
            "storm sequence 4",
            "storm sequence 5",
        ]
    )

    axs1[2].plot(
        time[0],
        dune_height[0],
        time[1],
        dune_height[1],
        time[2],
        dune_height[2],
        time[3],
        dune_height[3],
        time[4],
        dune_height[4],
    )
    axs1[2].set(ylabel="dune height (m)")
    axs1[2].set(xlabel="years")
    axs1[2].set_xlim([tlow, thigh])


# # ###############################################################################
# # old plotters when exploring nonlinear dynamics
# # ###############################################################################


def nonlinear_animated(CASCADE_b3d, ib3d, tmin, tmax, name):
    # mean dune height ---------
    # sub_domain = CASCADE_b3d[ib3d]._DuneDomain[:, :, :]
    # DuneCrest = sub_domain.max(axis=2)
    # DuneCrestMean = np.mean(DuneCrest, axis=1) * 10
    DuneCrestMean = [
        a * 10 for a in CASCADE_b3d[ib3d]._Hd_AverageTS
    ]  # don't know why this is slightly different, but it is

    # barrier width
    BarrierWidth = (
        np.array(CASCADE_b3d[ib3d].x_b_TS) - np.array(CASCADE_b3d[ib3d].x_s_TS)
    ) * 10

    # change in barrier width (to accommodate lag?)
    scts = [(x - BarrierWidth[0]) for x in BarrierWidth]
    rate = [0]
    for k in range(1, len(scts)):
        rate.append(scts[k] - scts[k - 1])
    bw_rate = (
        rate  # note, np.diff doesn't start with zero rate of change, so we do this calc
    )

    # shoreline change rate
    scts = [(x - CASCADE_b3d[ib3d]._x_s_TS[0]) * 10 for x in CASCADE_b3d[ib3d]._x_s_TS]
    rate = [0]
    for k in range(1, len(scts)):
        rate.append(scts[k] - scts[k - 1])
    sc_rate = rate

    # overwash flux
    Qoverwash = CASCADE_b3d[ib3d].QowTS

    # try an animated scatter and line plot ---------------------
    os.chdir("/Users/KatherineAnardeWheels/PycharmProjects/CASCADE/Output")
    newpath = name + "/test" + "/SimFrames/"
    if not os.path.exists(newpath):
        os.makedirs(newpath)
    os.chdir(newpath)

    nt = tmax - tmin
    colors = mpl.cm.viridis(np.linspace(0, 1, nt))

    for t in range(nt - 1):
        fig = plt.figure()
        ax = plt.axes(projection="3d")
        ax.set_xlabel("Dune height (m)")
        ax.set_zlabel("Barrier width change rate (m/yr)")
        # ax.set_ylabel("Overwash flux ($m^3/m$)")
        ax.set_zlabel("Shoreline change rate (m/yr)")

        ax.set_title(str(tmin + t) + " yrs")

        ax.plot3D(
            DuneCrestMean[tmin : tmin + t],
            bw_rate[tmin : tmin + t],
            sc_rate[tmin : tmin + t],
            "gray",
        )
        # three-dimensional scattered points
        ax.scatter3D(
            DuneCrestMean[tmin : tmin + t],
            bw_rate[tmin : tmin + t],
            sc_rate[tmin : tmin + t],
            c=colors[0:t],
        )
        ax.set_xlim([np.min(DuneCrestMean), np.max(DuneCrestMean)])
        ax.set_ylim([np.min(bw_rate), np.max(bw_rate)])
        ax.set_zlim([np.min(sc_rate), np.max(sc_rate)])

        fig.savefig(str(t))  # dpi=200
        plt.close(fig)

    frames = []

    for filenum in range(nt - 1):
        filename = str(filenum) + ".png"
        frames.append(imageio.imread(filename))
    imageio.mimsave("triplot.gif", frames, "GIF-FI")
    print()
    print("[ * GIF successfully generated * ]")

    # three-dimensional surface
    # ax.plot_surface(
    #     DuneCrestMean[tmin:tmax],
    #     BarrierWidth[tmin:tmax],
    #     sc_rate[tmin:tmax], # this needs to be two dimensional, a function of x and y
    #     rstride=1,
    #     cstride=1,
    #     cmap="viridis",
    #     edgecolor="none",
    # )
    # ax.set_title("surface")

    # export text files for mtm psd plots and coherence


def nonlinear_comparison(
    DuneCrestMean_45,
    BarrierWidth_45,
    DuneCrestMean_55,
    BarrierWidth_55,
    DuneCrestMean_65,
    BarrierWidth_65,
    DuneCrestMean_75,
    BarrierWidth_75,
):
    # fig, axs = plt.subplots(1, 4, figsize=(10, 5), sharey=True)
    plt.subplot(1, 4, 1)
    plt.scatter(
        DuneCrestMean_45,
        BarrierWidth_45,
        c=np.arange(0, np.size(DuneCrestMean_45), 1),
        cmap=cm.viridis,
    )
    # plt.plot(DuneCrestMean[tmin:tmax], BarrierWidth[tmin:tmax], "m")
    plt.xlabel("Ave Dune Height (m)")
    plt.ylabel("Barrier width (m)")
    plt.ylim([100, 500])

    plt.subplot(1, 4, 2)
    plt.scatter(
        DuneCrestMean_55,
        BarrierWidth_55,
        c=np.arange(0, np.size(DuneCrestMean_55), 1),
        cmap=cm.viridis,
    )
    # plt.plot(DuneCrestMean[tmin:tmax], BarrierWidth[tmin:tmax], "m")
    plt.xlabel("Ave Dune Height (m)")
    plt.ylim([100, 500])

    plt.subplot(1, 4, 3)
    plt.scatter(
        DuneCrestMean_65,
        BarrierWidth_65,
        c=np.arange(0, np.size(DuneCrestMean_65), 1),
        cmap=cm.viridis,
    )
    # plt.plot(DuneCrestMean[tmin:tmax], BarrierWidth[tmin:tmax], "m")
    plt.xlabel("Ave Dune Height (m)")
    plt.ylim([100, 500])

    plt.subplot(1, 4, 4)
    plt.scatter(
        DuneCrestMean_75,
        BarrierWidth_75,
        c=np.arange(0, np.size(DuneCrestMean_75), 1),
        cmap=cm.viridis,
    )
    # plt.plot(DuneCrestMean[tmin:tmax], BarrierWidth[tmin:tmax], "m")
    plt.xlabel("Ave Dune Height (m)")
    plt.ylim([100, 500])

    # plt.subplot(1, 2, 1)
    # plt.hist(BarrierWidth_75, density=True, bins=30)  # density=False would make counts
    # plt.ylabel("Probability")
    # plt.xlabel("Barrier Width - high dune growth rate")
    # plt.xlim([100, 350])
    #
    # plt.subplot(1, 2, 2)
    # plt.hist(BarrierWidth_75[9000:9999], density=True, bins=30)
    # plt.xlabel("Barrier Width - high dune growth rate")
    # plt.xlim([100, 350])
    #
    # plt.subplot(1, 2, 1)
    # plt.hist(BarrierWidth_45, density=True, bins=30)  # density=False would make counts
    # plt.ylabel("Probability")
    # plt.xlabel("Barrier Width - low dune growth rate")
    # plt.xlim([320, 460])
    #
    # plt.subplot(1, 2, 2)
    # plt.hist(BarrierWidth_45[9000:9999], density=True, bins=30)
    # plt.xlabel("Barrier Width - low dune growth rate")
    # plt.xlim([320, 460])
