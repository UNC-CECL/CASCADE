# Plotting functions for CASCADE

import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib import cm
import numpy as np
import pandas as pd
import os
import imageio
import math
from scipy import signal

# import plotly.graph_objects as go

# import seaborn as sns


# ===================================================
# 1: Animation Frames of Barrier and Dune Elevation (#4 in Barrier3D_Functions, modified here for sequential
# versions of CASCADE)
#
#       inputs:         - barrier3d (a list containing BMI objects)
#                       - ny (the number of barrier3d subgrids that you want to plot)
#                       - directory (for saving)
#                       - TMAX (the last time index for plotting)
#                       - name (for saving new directory)
#       outputs:        - gif

# barrier3d only
def plot_ElevAnimation(barrier3d, ny, directory, TMAX, name):
    BarrierLength = barrier3d[0]._BarrierLength

    BeachWidth = 6
    OriginY = int(barrier3d[0]._x_s_TS[0] - barrier3d[0]._x_t_TS[0])
    AniDomainWidth = int(
        np.amax(barrier3d[0]._InteriorWidth_AvgTS)
        + BeachWidth
        + np.abs(barrier3d[0]._ShorelineChange)
        + OriginY
        + 35
    )

    os.chdir(directory)
    newpath = "Output/" + name + "/SimFrames/"
    if not os.path.exists(newpath):
        os.makedirs(newpath)
    os.chdir(newpath)

    for t in range(TMAX - 1):

        AnimateDomain = np.ones([AniDomainWidth + 1, BarrierLength * ny]) * -1

        for iB3D in range(ny):
            # Build beach elevation domain
            BeachDomain = np.zeros([BeachWidth, BarrierLength])
            berm = math.ceil(BeachWidth * 0.65)
            BeachDomain[berm : BeachWidth + 1, :] = barrier3d[iB3D]._BermEl
            add = (barrier3d[iB3D]._BermEl - barrier3d[iB3D]._SL) / berm
            for i in range(berm):
                BeachDomain[i, :] = barrier3d[iB3D]._SL + add * i

            # Make animation frame domain
            Domain = barrier3d[iB3D]._DomainTS[t] * 10
            Dunes = (
                barrier3d[iB3D]._DuneDomain[t, :, :] + barrier3d[iB3D]._BermEl
            ) * 10
            Dunes = np.rot90(Dunes)
            Dunes = np.flipud(Dunes)
            Beach = BeachDomain * 10
            Domain = np.vstack([Beach, Dunes, Domain])
            Domain[Domain < 0] = -1
            widthTS = len(Domain)
            scts = [(x - barrier3d[iB3D]._x_s_TS[0]) for x in barrier3d[iB3D]._x_s_TS]
            if scts[t] >= 0:
                OriginTstart = OriginY + math.floor(scts[t])
            else:
                OriginTstart = OriginY + math.ceil(scts[t])
            OriginTstop = OriginTstart + widthTS
            xOrigin = iB3D * BarrierLength
            AnimateDomain[
                OriginTstart:OriginTstop, xOrigin : xOrigin + BarrierLength
            ] = Domain

        # Plot and save
        elevFig1 = plt.figure(figsize=(15, 7))
        ax = elevFig1.add_subplot(111)
        cax = ax.matshow(
            AnimateDomain, origin="lower", cmap="terrain", vmin=-1.1, vmax=5.0
        )  # , interpolation='gaussian') # analysis:ignore
        ax.xaxis.set_ticks_position("bottom")
        elevFig1.colorbar(cax)
        plt.xlabel("Alongshore Distance (dam)")
        plt.ylabel("Cross-Shore Distance (dam)")
        plt.title("Interior Elevation")
        plt.tight_layout()
        timestr = "Time = " + str(t) + " yrs"
        plt.text(1, 1, timestr)
        plt.rcParams.update({"font.size": 20})
        name = "elev_" + str(t)
        elevFig1.savefig(name)  # dpi=200
        plt.close(elevFig1)

    frames = []

    for filenum in range(TMAX - 1):
        filename = "elev_" + str(filenum) + ".png"
        frames.append(imageio.imread(filename))
    imageio.mimsave("elev.gif", frames, "GIF-FI")
    print()
    print("[ * GIF successfully generated * ]")


# changed limits and labels for AGU
def plot_ElevAnimationAGU(barrier3d, ny, directory, TMAX, name):
    BarrierLength = barrier3d[0]._BarrierLength

    BeachWidth = 6
    OriginY = int(barrier3d[0]._x_s_TS[0] - barrier3d[0]._x_t_TS[0])
    OriginY = int(barrier3d[0]._x_s_TS[0] - barrier3d[0]._x_t_TS[0])
    AniDomainWidth = int(
        np.amax(barrier3d[0]._InteriorWidth_AvgTS)
        + BeachWidth
        + np.abs(barrier3d[0]._ShorelineChange)
        + OriginY
        + 35
    )

    os.chdir(directory)
    newpath = "Output/" + name + "/SimFrames/"
    if not os.path.exists(newpath):
        os.makedirs(newpath)
    os.chdir(newpath)

    for t in range(TMAX - 1):  # range(TMAX-2, TMAX-1):# range(TMAX-1):

        AnimateDomain = np.ones([AniDomainWidth + 1, BarrierLength * ny]) * -1

        for iB3D in range(ny):
            # Build beach elevation domain
            BeachDomain = np.zeros([BeachWidth, BarrierLength])
            berm = math.ceil(BeachWidth * 0.65)
            BeachDomain[berm : BeachWidth + 1, :] = barrier3d[iB3D]._BermEl
            add = (barrier3d[iB3D]._BermEl - barrier3d[iB3D]._SL) / berm
            for i in range(berm):
                BeachDomain[i, :] = barrier3d[iB3D]._SL + add * i

            # Make animation frame domain
            Domain = barrier3d[iB3D]._DomainTS[t] * 10
            Dunes = (
                barrier3d[iB3D]._DuneDomain[t, :, :] + barrier3d[iB3D]._BermEl
            ) * 10
            Dunes = np.rot90(Dunes)
            Dunes = np.flipud(Dunes)
            Beach = BeachDomain * 10
            Domain = np.vstack([Beach, Dunes, Domain])
            Domain[Domain < 0] = -1
            widthTS = len(Domain)
            scts = [(x - barrier3d[iB3D]._x_s_TS[0]) for x in barrier3d[iB3D]._x_s_TS]
            if scts[t] >= 0:
                OriginTstart = OriginY + math.floor(scts[t])
            else:
                OriginTstart = OriginY + math.ceil(scts[t])
            OriginTstop = OriginTstart + widthTS
            xOrigin = iB3D * BarrierLength
            AnimateDomain[
                OriginTstart:OriginTstop, xOrigin : xOrigin + BarrierLength
            ] = Domain

        # Plot and save
        elevFig1 = plt.figure(figsize=(15, 7))
        ax = elevFig1.add_subplot(111)
        cax = ax.matshow(
            AnimateDomain, origin="lower", cmap="terrain", vmin=-1.1, vmax=4.0
        )  # , interpolation='gaussian') # analysis:ignore
        ax.xaxis.set_ticks_position("bottom")
        # cbar = elevFig1.colorbar(cax)
        plt.xlabel("Alongshore Distance (dam)")
        plt.ylabel("Cross-Shore Distance (dam)")
        plt.title("Interior Elevation")
        plt.tight_layout()
        timestr = "Time = " + str(t) + " yrs"
        plt.text(1, 1, timestr)
        plt.rcParams.update({"font.size": 20})
        plt.xlim([0, 599])  # hard coded
        plt.ylim([0, 217])
        locs, labels = plt.xticks()
        plt.xticks(locs, ["0", "1", "2", "3", "4", "5", "6"])
        locs, labels = plt.yticks()
        plt.yticks(locs, ["0", "0.5", "1", "1.5", "2", "2.5"])
        # plt.xlim([0,299])
        # plt.ylim([0,277])
        plt.xlim([0, 599])
        plt.ylim([0, 217])
        plt.xlabel("Alongshore Distance (km)")
        plt.ylabel("Cross-Shore Distance (km)")
        name = "elev_" + str(t)
        elevFig1.savefig(name)  # dpi=200
        plt.close(elevFig1)

    frames = []

    for filenum in range(TMAX - 1):
        filename = "elev_" + str(filenum) + ".png"
        frames.append(imageio.imread(filename))
    imageio.mimsave("elev.gif", frames, "GIF-FI")
    print()
    print("[ * GIF successfully generated * ]")


# cascade
def plot_ElevAnimation_CASCADE(
    cascade,
    ny,
    directory,
    TMAX_MGMT,  # an array
    name,
    TMAX_SIM,  # float (same for all)
    beach_management_ny=None,  # list of booleans
    roadway_management_ny=None,
    y_lim=None,
    fig_size=None,
):
    """
    NOTE THAT THE BEACH REPRESENTATION IS BASED ON A MODEL SPECIFIED BEACH WIDTH. We set the beach width for the
    remaining time steps after the community has been abandoned to the last managed beach width in order to not have a
    huge jump in the back-barrier position in Barrier3D. OTHERWISE, it is meaningless to the dynamics Barrier3D.
    """
    barrier3d = cascade.barrier3d

    # set up the domain; here we just use the first grid, but that could break in future runs
    BarrierLength = barrier3d[0].BarrierLength
    if beach_management_ny is not None:
        indices = [
            i for i in range(len(beach_management_ny)) if beach_management_ny[i] == 1
        ]
        iB3D = indices[0]
        MaxBeachWidth = (
            np.max(cascade.nourishments[iB3D].beach_width[0 : TMAX_MGMT[iB3D]]) / 10
        )  # dam
    else:
        MaxBeachWidth = cascade._initial_beach_width[0]
    OriginY = int(barrier3d[0].x_s_TS[0])
    AniDomainWidth = int(
        np.amax(barrier3d[0].InteriorWidth_AvgTS)
        + MaxBeachWidth
        + np.abs(barrier3d[0]._ShorelineChange)
        + OriginY
        + 35
    )

    os.chdir(directory)
    newpath = "Output/" + name + "/SimFrames/"
    if not os.path.exists(newpath):
        os.makedirs(newpath)
    os.chdir(newpath)

    maxMGMT = np.max(TMAX_MGMT)

    for t in range(maxMGMT + 1):

        # start with plotting t=0, then plot post-storm dune, interior, shoreface, etc. before management, treating this
        # as t=0.5 (i.e., this is really the final configuration from storms at t=1,2,3,4,5,... but before human
        # modifications); for the default dune inputs, we should be able to see natural dune growth at
        # t=0.5 before rebuild and storms (i.e., t=0.5 and t=1 should not be the same)
        if 0 < t <= TMAX_SIM:

            # post-storm variables in the BeachDuneManager are: interior, dune height, x_s, s_sf, beach width
            AnimateDomain = np.ones([AniDomainWidth + 1, BarrierLength * ny]) * -1

            for iB3D in range(ny):

                # the nourishment scenario
                if beach_management_ny[iB3D] and (t < TMAX_MGMT[iB3D] + 1):
                    actual_shoreline_post_storm = np.hstack(
                        [
                            0,
                            cascade.nourishments[iB3D]._post_storm_x_s[
                                1 : TMAX_MGMT[iB3D] + 1
                            ],
                        ]
                    )
                    beach_width = (
                        cascade.nourishments[iB3D]._post_storm_beach_width[t] / 10
                    )
                    Domain = cascade.nourishments[iB3D]._post_storm_interior[t] * 10
                    Dunes = (
                        cascade.nourishments[iB3D]._post_storm_dunes[t]
                        + barrier3d[iB3D].BermEl
                    ) * 10

                elif roadway_management_ny[iB3D] and (t < TMAX_MGMT[iB3D] + 1):
                    # the roadways scenario
                    actual_shoreline_post_storm = barrier3d[iB3D].x_s_TS[
                        0 : TMAX_MGMT[iB3D] + 1
                    ]
                    beach_width = cascade._initial_beach_width[iB3D] / 10
                    Domain = cascade.roadways[iB3D]._post_storm_interior[t] * 10
                    Dunes = (
                        cascade.roadways[iB3D]._post_storm_dunes[t]
                        + barrier3d[iB3D].BermEl
                    ) * 10

                else:
                    # the natural scenario
                    actual_shoreline_post_storm = barrier3d[iB3D].x_s_TS[0:TMAX_SIM]
                    beach_width = cascade._initial_beach_width[iB3D] / 10
                    Domain = barrier3d[iB3D].DomainTS[t] * 10
                    Dunes = (
                        barrier3d[iB3D].DuneDomain[t, :, :] + barrier3d[iB3D].BermEl
                    ) * 10

                # Build beach elevation domain, we only show beach width decreasing in increments of 10 m and we don't
                # illustrate a berm, just a sloping beach up to the elevation of the berm
                cellular_dune_toe_post_storm = np.floor(
                    actual_shoreline_post_storm[t] + beach_width
                )
                cellular_shoreline_post_storm = np.floor(actual_shoreline_post_storm[t])
                cellular_beach_width = int(
                    cellular_dune_toe_post_storm - cellular_shoreline_post_storm
                )  # not actual bw

                BeachDomain = np.zeros(
                    [
                        cellular_beach_width,
                        BarrierLength,
                    ]
                )
                if cellular_beach_width == 0:
                    pass
                else:
                    add = (barrier3d[iB3D].BermEl - barrier3d[iB3D].SL) / (
                        cellular_beach_width + 1
                    )
                    for i in range(0, cellular_beach_width):
                        BeachDomain[i, :] = (barrier3d[iB3D].SL + add) * (i + 1)

                # Make animation frame domain
                Dunes = np.rot90(Dunes)
                Dunes = np.flipud(Dunes)
                Beach = BeachDomain * 10
                Domain = np.vstack([Beach, Dunes, Domain])
                Domain[Domain < 0] = -1
                widthTS = len(Domain)
                OriginTstart = int(cellular_shoreline_post_storm)
                OriginTstop = OriginTstart + widthTS
                xOrigin = iB3D * BarrierLength
                AnimateDomain[
                    OriginTstart:OriginTstop, xOrigin : xOrigin + BarrierLength
                ] = Domain

            # Plot and save
            if fig_size is not None:
                elevFig1 = plt.figure(figsize=fig_size)
            else:
                elevFig1 = plt.figure(figsize=(7, 7))
            ax = elevFig1.add_subplot(111)
            cax = ax.pcolormesh(
                AnimateDomain,
                cmap="terrain",
                vmin=-1.1,
                vmax=3.5,
                # edgecolors="w",  # for debugging
                # linewidth=0.01,
            )
            cbar = elevFig1.colorbar(cax)
            cbar.set_label("Elevation (m MHW)", rotation=270)
            plt.xlabel("Alongshore Distance (dam)")
            plt.ylabel("Cross-Shore Distance (dam)")
            timestr = (
                "Time = " + str(t - 0.5) + " yrs"
            )  # we are letting the post-storm output represent 0.5 years
            if y_lim is not None:
                plt.ylim(y_lim)
                plt.text(3, y_lim[0] + 3, timestr, color="w")
            else:
                plt.ylim(bottom=OriginY - 35)
                plt.text(1, OriginY - 33, timestr)
            plt.tight_layout()
            plt.rcParams.update({"font.size": 16})
            name = "elev_" + str(t - 1) + "pt5"
            # elevFig1.tight_layout()
            elevFig1.savefig(name)  # dpi=200
            plt.close(elevFig1)

    for t in range(TMAX_SIM):

        # ok, now the annual time step, which incorporates human modifications to the shoreface, beach, dune, & interior
        AnimateDomain = np.ones([AniDomainWidth + 1, BarrierLength * ny]) * -1

        for iB3D in range(ny):

            actual_shoreline_post_humans = barrier3d[iB3D].x_s_TS[0 : TMAX_SIM + 1]

            if beach_management_ny[iB3D]:
                beach_width = cascade.nourishments[iB3D].beach_width[t] / 10
            else:
                # both roadways and natural scenario
                beach_width = cascade._initial_beach_width[iB3D] / 10

            cellular_dune_toe_post_humans = np.floor(
                actual_shoreline_post_humans[t] + beach_width
            )

            cellular_shoreline_post_humans = np.floor(actual_shoreline_post_humans[t])
            cellular_beach_width = int(
                cellular_dune_toe_post_humans - cellular_shoreline_post_humans
            )  # not actual bw

            if t == maxMGMT + 1:
                cellular_beach_width_final = (
                    cellular_beach_width  # save this for plotting moving forward
                )
            if t > maxMGMT + 1:
                cellular_beach_width = cellular_beach_width_final

            BeachDomain = np.zeros(
                [
                    cellular_beach_width,
                    BarrierLength,
                ]
            )

            if cellular_beach_width == 0:
                pass
            else:
                add = (barrier3d[iB3D].BermEl - barrier3d[iB3D].SL) / (
                    cellular_beach_width + 1
                )
                for i in range(0, cellular_beach_width):
                    BeachDomain[i, :] = (barrier3d[iB3D].SL + add) * (i + 1)

            # Make animation frame domain
            Domain = barrier3d[iB3D].DomainTS[t] * 10
            Dunes = (barrier3d[iB3D].DuneDomain[t, :, :] + barrier3d[iB3D].BermEl) * 10
            Dunes = np.rot90(Dunes)
            Dunes = np.flipud(Dunes)
            Beach = BeachDomain * 10
            Domain = np.vstack([Beach, Dunes, Domain])
            Domain[Domain < 0] = -1
            widthTS = len(Domain)
            OriginTstart = int(cellular_shoreline_post_humans)
            OriginTstop = OriginTstart + widthTS
            xOrigin = iB3D * BarrierLength
            AnimateDomain[
                OriginTstart:OriginTstop, xOrigin : xOrigin + BarrierLength
            ] = Domain

        # Plot and save
        if fig_size is not None:
            elevFig2 = plt.figure(figsize=fig_size)
        else:
            elevFig2 = plt.figure(figsize=(7, 7))
        ax = elevFig2.add_subplot(111)
        cax = ax.pcolormesh(
            AnimateDomain,
            cmap="terrain",
            vmin=-1.1,
            vmax=3.5,
            # edgecolors="w",  # for debugging
            # linewidth=0.01,
        )
        cbar = elevFig2.colorbar(cax)
        cbar.set_label("Elevation (m MHW)", rotation=270)
        plt.xlabel("Alongshore Distance (dam)")
        plt.ylabel("Cross-Shore Distance (dam)")
        timestr = "Time = " + str(t) + " yrs"
        if y_lim is not None:
            plt.ylim(y_lim)
            plt.text(3, y_lim[0] + 3, timestr, color="w")
        else:
            plt.ylim(bottom=OriginY - 35)
            plt.text(1, OriginY - 33, timestr)
        plt.tight_layout()
        plt.rcParams.update({"font.size": 16})
        name = "elev_" + str(t)
        # elevFig2.tight_layout()
        elevFig2.savefig(name)  # dpi=200
        plt.close(elevFig2)

    frames = []

    for filenum in range(TMAX_SIM):
        filename = "elev_" + str(filenum) + ".png"
        frames.append(imageio.imread(filename))

        if (
            filenum
            < TMAX_MGMT[iB3D]
            # and cascade.nourishments[iB3D]._post_storm_interior[TMAX_MGMT] is not None
        ):
            filename = "elev_" + str(filenum) + "pt5" ".png"
            frames.append(imageio.imread(filename))

    imageio.mimsave("elev.gif", frames, "GIF-FI")
    print()
    print("[ * GIF successfully generated * ]")


# modified for use with the RoadwayManager Module
# -- NOTE THAT THE BEACH REPRESENTATION IS NOT REALISTIC (i.e., there is no beach width in the RoadwayManager module) --
def plot_ElevAnimation_Humans_Roadways(
    cascade,
    ny,
    directory,
    TMAX,
    TMAX_roadways,
    name,
):
    barrier3d = cascade.barrier3d

    BarrierLength = barrier3d[0]._BarrierLength

    BeachWidth = 3
    OriginY = int(barrier3d[0]._x_s_TS[0] - barrier3d[0]._x_t_TS[0])
    AniDomainWidth = int(
        np.amax(barrier3d[0]._InteriorWidth_AvgTS)
        + BeachWidth
        + np.abs(barrier3d[0]._ShorelineChange)
        + OriginY
        + 35
    )

    os.chdir(directory)
    newpath = "Output/" + name + "/SimFrames/"
    if not os.path.exists(newpath):
        os.makedirs(newpath)
    os.chdir(newpath)

    for t in range(TMAX - 1):

        # start with plotting t=0, then plot post-storm dune and interior before rebuild, treating this as t=0.5 (i.e.,
        # this is really the final dune from storms at t=1,2,3,4,5,... but before human modifications to the dune
        # (rebuild, move overwash); for the default dune inputs, we should be able to see natural dune growth at
        # t=0.5 before rebuild and storms (i.e., t=0.5 and t=1 should not be the same)
        if t > 0 and cascade.roadways[0]._post_storm_interior[t] is not None:

            AnimateDomain = np.ones([AniDomainWidth + 1, BarrierLength * ny]) * -1

            for iB3D in range(ny):
                # Build beach elevation domain
                BeachDomain = np.zeros([BeachWidth, BarrierLength])
                berm = math.ceil(BeachWidth * 0.65)
                BeachDomain[berm : BeachWidth + 1, :] = barrier3d[iB3D]._BermEl
                add = (barrier3d[iB3D]._BermEl - barrier3d[iB3D]._SL) / berm
                for i in range(berm):
                    BeachDomain[i, :] = barrier3d[iB3D]._SL + add * i

                # Make animation frame domain
                Domain = cascade.roadways[iB3D]._post_storm_interior[t] * 10
                Dunes = (
                    cascade.roadways[iB3D]._post_storm_dunes[t]
                    + barrier3d[iB3D]._BermEl
                ) * 10
                Dunes = np.rot90(Dunes)
                Dunes = np.flipud(Dunes)
                Beach = BeachDomain * 10
                Domain = np.vstack([Beach, Dunes, Domain])
                Domain[Domain < 0] = -1
                widthTS = len(Domain)
                scts = [
                    (x - barrier3d[iB3D]._x_s_TS[0]) for x in barrier3d[iB3D]._x_s_TS
                ]
                # if scts[t] >= 0:
                OriginTstart = OriginY + math.floor(scts[t])
                # else:
                #     OriginTstart = OriginY + math.ceil(scts[t])
                OriginTstop = OriginTstart + widthTS
                xOrigin = iB3D * BarrierLength
                AnimateDomain[
                    OriginTstart:OriginTstop, xOrigin : xOrigin + BarrierLength
                ] = Domain

            # Plot and save
            elevFig1 = plt.figure(figsize=(7, 7))
            ax = elevFig1.add_subplot(111)
            cax = ax.matshow(
                AnimateDomain, origin="lower", cmap="terrain", vmin=-1.1, vmax=5.0
            )  # , interpolation='gaussian') # analysis:ignore
            ax.xaxis.set_ticks_position("bottom")
            cbar = elevFig1.colorbar(cax)
            cbar.set_label("Elevation (m MHW)", rotation=270)
            plt.xlabel("Alongshore Distance (dam)")
            plt.ylabel("Cross-Shore Distance (dam)")
            plt.title("Interior Elevation")
            plt.tight_layout()
            timestr = (
                "Time = " + str(t - 0.5) + " yrs"
            )  # we are letting the post-storm output represent 0.5 years
            plt.text(1, 1, timestr)
            plt.rcParams.update({"font.size": 14})
            name = "elev_" + str(t - 1) + "pt5"
            elevFig1.savefig(name)  # dpi=200
            plt.close(elevFig1)

        AnimateDomain = np.ones([AniDomainWidth + 1, BarrierLength * ny]) * -1

        # annual time steps
        for iB3D in range(ny):
            # Build beach elevation domain
            BeachDomain = np.zeros([BeachWidth, BarrierLength])
            berm = math.ceil(BeachWidth * 0.65)
            BeachDomain[berm : BeachWidth + 1, :] = barrier3d[iB3D]._BermEl
            add = (barrier3d[iB3D]._BermEl - barrier3d[iB3D]._SL) / berm
            for i in range(berm):
                BeachDomain[i, :] = barrier3d[iB3D]._SL + add * i

            # Make animation frame domain
            Domain = barrier3d[iB3D]._DomainTS[t] * 10
            Dunes = (
                barrier3d[iB3D]._DuneDomain[t, :, :] + barrier3d[iB3D]._BermEl
            ) * 10
            Dunes = np.rot90(Dunes)
            Dunes = np.flipud(Dunes)
            Beach = BeachDomain * 10
            Domain = np.vstack([Beach, Dunes, Domain])
            Domain[Domain < 0] = -1
            widthTS = len(Domain)
            scts = [(x - barrier3d[iB3D]._x_s_TS[0]) for x in barrier3d[iB3D]._x_s_TS]
            # if scts[t] >= 0:
            OriginTstart = OriginY + math.floor(scts[t])
            # else:
            #     OriginTstart = OriginY + math.ceil(scts[t])
            OriginTstop = OriginTstart + widthTS
            xOrigin = iB3D * BarrierLength
            AnimateDomain[
                OriginTstart:OriginTstop, xOrigin : xOrigin + BarrierLength
            ] = Domain

        # Plot and save
        elevFig1 = plt.figure(figsize=(7, 7))
        ax = elevFig1.add_subplot(111)
        cax = ax.matshow(
            AnimateDomain, origin="lower", cmap="terrain", vmin=-1.1, vmax=5.0
        )  # , interpolation='gaussian') # analysis:ignore
        ax.xaxis.set_ticks_position("bottom")
        cbar = elevFig1.colorbar(cax)
        cbar.set_label("Elevation (m MHW)", rotation=270)
        plt.xlabel("Alongshore Distance (dam)")
        plt.ylabel("Cross-Shore Distance (dam)")
        plt.tight_layout()
        timestr = "Time = " + str(t) + " yrs"
        plt.text(1, 1, timestr)
        plt.rcParams.update({"font.size": 14})
        name = "elev_" + str(t)
        elevFig1.savefig(name)  # dpi=200
        plt.close(elevFig1)

    frames = []

    for filenum in range(TMAX - 1):
        filename = "elev_" + str(filenum) + ".png"
        frames.append(imageio.imread(filename))

        if (
            filenum < TMAX - 2
            # and cascade.roadways[iB3D]._post_storm_interior[TMAX - 2] is not None
            # and cascade.roadways[iB3D]._post_storm_interior[filenum] is not None
            and filenum < TMAX_roadways
        ):
            filename = "elev_" + str(filenum) + "pt5" ".png"
            frames.append(imageio.imread(filename))

    imageio.mimsave("elev.gif", frames, "GIF-FI")
    print()
    print("[ * GIF successfully generated * ]")


# for use with BeachDuneManager Module
def plot_ElevAnimation_Humans_BeachDuneManager(
    cascade, ny, directory, TMAX, name, TMAX_SIM  # management max  # simulation max
):
    """
    NOTE THAT THE BEACH REPRESENTATION IS BASED ON A MODEL SPECIFIED BEACH WIDTH. We set the beach width for the
    remaining time steps after the community has been abandoned to the last managed beach width in order to not have a
    huge jump in the back-barrier position in Barrier3D. OTHERWISE, it is meaningless to the dynamics Barrier3D.
    """
    barrier3d = cascade.barrier3d

    # set up the domain; here we just use the first grid, but that could break in future runs
    BarrierLength = barrier3d[0].BarrierLength
    MaxBeachWidth = np.max(cascade.nourishments[0].beach_width[0:TMAX]) / 10  # dam
    OriginY = int(barrier3d[0].x_s_TS[0])
    AniDomainWidth = int(
        np.amax(barrier3d[0].InteriorWidth_AvgTS)
        + MaxBeachWidth
        + np.abs(barrier3d[0]._ShorelineChange)
        + OriginY
        + 35
    )

    os.chdir(directory)
    newpath = "Output/" + name + "/SimFrames/"
    if not os.path.exists(newpath):
        os.makedirs(newpath)
    os.chdir(newpath)

    for t in range(TMAX + 1):

        # start with plotting t=0, then plot post-storm dune, interior, shoreface, etc. before management, treating this
        # as t=0.5 (i.e., this is really the final configuration from storms at t=1,2,3,4,5,... but before human
        # modifications); for the default dune inputs, we should be able to see natural dune growth at
        # t=0.5 before rebuild and storms (i.e., t=0.5 and t=1 should not be the same)

        if t > 0 and cascade.nourishments[0]._post_storm_interior[t] is not None:

            # post-storm variables in the BeachDuneManager are: interior, dune height, x_s, s_sf, beach width
            AnimateDomain = np.ones([AniDomainWidth + 1, BarrierLength * ny]) * -1

            for iB3D in range(ny):

                actual_shoreline_post_storm = np.hstack(
                    [0, cascade.nourishments[iB3D]._post_storm_x_s[1 : TMAX + 1]]
                )

                # Build beach elevation domain, we only show beach width decreasing in increments of 10 m and we don't
                # illustrate a berm, just a sloping beach up to the elevation of the berm
                cellular_dune_toe_post_storm = np.floor(
                    actual_shoreline_post_storm[t]
                    + (cascade.nourishments[iB3D]._post_storm_beach_width[t] / 10)
                )
                cellular_shoreline_post_storm = np.floor(actual_shoreline_post_storm[t])
                cellular_beach_width = int(
                    cellular_dune_toe_post_storm - cellular_shoreline_post_storm
                )  # not actual bw

                BeachDomain = np.zeros(
                    [
                        cellular_beach_width,
                        BarrierLength,
                    ]
                )
                if cellular_beach_width == 0:
                    pass
                else:
                    add = (barrier3d[iB3D].BermEl - barrier3d[iB3D].SL) / (
                        cellular_beach_width + 1
                    )
                    for i in range(0, cellular_beach_width):
                        BeachDomain[i, :] = (barrier3d[iB3D].SL + add) * (i + 1)

                # Make animation frame domain
                Domain = cascade.nourishments[iB3D]._post_storm_interior[t] * 10
                Dunes = (
                    cascade.nourishments[iB3D]._post_storm_dunes[t]
                    + barrier3d[iB3D].BermEl
                ) * 10
                Dunes = np.rot90(Dunes)
                Dunes = np.flipud(Dunes)
                Beach = BeachDomain * 10
                Domain = np.vstack([Beach, Dunes, Domain])
                Domain[Domain < 0] = -1
                widthTS = len(Domain)
                OriginTstart = int(cellular_shoreline_post_storm)
                OriginTstop = OriginTstart + widthTS
                xOrigin = iB3D * BarrierLength
                AnimateDomain[
                    OriginTstart:OriginTstop, xOrigin : xOrigin + BarrierLength
                ] = Domain

            # Plot and save
            elevFig1 = plt.figure(figsize=(7, 7))
            ax = elevFig1.add_subplot(111)
            cax = ax.pcolormesh(
                AnimateDomain,
                cmap="terrain",
                vmin=-1.1,
                vmax=5.0,
                # edgecolors="w",  # for debugging
                # linewidth=0.01,
            )
            cbar = elevFig1.colorbar(cax)
            cbar.set_label("Elevation (m MHW)", rotation=270)
            plt.xlabel("Alongshore Distance (dam)")
            plt.ylabel("Cross-Shore Distance (dam)")
            plt.ylim(bottom=OriginY - 35)
            plt.tight_layout()
            timestr = (
                "Time = " + str(t - 0.5) + " yrs"
            )  # we are letting the post-storm output represent 0.5 years
            plt.text(1, OriginY - 33, timestr)
            plt.rcParams.update({"font.size": 16})
            name = "elev_" + str(t - 1) + "pt5"
            # elevFig1.tight_layout()
            elevFig1.savefig(name)  # dpi=200
            plt.close(elevFig1)

    for t in range(TMAX_SIM):

        # ok, now the annual time step, which incorporates human modifications to the shoreface, beach, dune, & interior
        AnimateDomain = np.ones([AniDomainWidth + 1, BarrierLength * ny]) * -1

        for iB3D in range(ny):

            actual_shoreline_post_humans = barrier3d[iB3D].x_s_TS[0:TMAX_SIM]

            # Build beach elevation domain, we only show beach width decreasing in increments of 10 m and we don't
            # illustrate a berm, just a sloping beach up to the elevation of the berm

            # if math.isnan(cascade.nourishments[iB3D].beach_width[t]):
            #     cellular_dune_toe_post_humans = np.floor(
            #         actual_shoreline_post_humans[t]
            #         + (cascade.nourishments[iB3D].beach_width[TMAX] / 10)
            #     )  # just set to the last TMAX so not weird jolt in shoreline position
            #     # cellular_dune_toe_post_humans = np.floor(
            #     #     actual_shoreline_post_humans[t]
            #     # )  # make bw zero
            # else:
            cellular_dune_toe_post_humans = np.floor(
                actual_shoreline_post_humans[t]
                + (cascade.nourishments[iB3D].beach_width[t] / 10)
            )

            cellular_shoreline_post_humans = np.floor(actual_shoreline_post_humans[t])
            cellular_beach_width = int(
                cellular_dune_toe_post_humans - cellular_shoreline_post_humans
            )  # not actual bw
            BeachDomain = np.zeros(
                [
                    cellular_beach_width,
                    BarrierLength,
                ]
            )
            if cellular_beach_width == 0:
                pass
            else:
                add = (barrier3d[iB3D].BermEl - barrier3d[iB3D].SL) / (
                    cellular_beach_width + 1
                )
                for i in range(0, cellular_beach_width):
                    BeachDomain[i, :] = (barrier3d[iB3D].SL + add) * (i + 1)

            # Make animation frame domain
            Domain = barrier3d[iB3D].DomainTS[t] * 10
            Dunes = (barrier3d[iB3D].DuneDomain[t, :, :] + barrier3d[iB3D].BermEl) * 10
            Dunes = np.rot90(Dunes)
            Dunes = np.flipud(Dunes)
            Beach = BeachDomain * 10
            Domain = np.vstack([Beach, Dunes, Domain])
            Domain[Domain < 0] = -1
            widthTS = len(Domain)
            OriginTstart = int(cellular_shoreline_post_humans)
            OriginTstop = OriginTstart + widthTS
            xOrigin = iB3D * BarrierLength
            AnimateDomain[
                OriginTstart:OriginTstop, xOrigin : xOrigin + BarrierLength
            ] = Domain

        # Plot and save
        elevFig2 = plt.figure(figsize=(7, 7))
        ax = elevFig2.add_subplot(111)
        cax = ax.pcolormesh(
            AnimateDomain,
            cmap="terrain",
            vmin=-1.1,
            vmax=5.0,
            # edgecolors="w",  # for debugging
            # linewidth=0.01,
        )
        cbar = elevFig2.colorbar(cax)
        cbar.set_label("Elevation (m MHW)", rotation=270)
        plt.xlabel("Alongshore Distance (dam)")
        plt.ylabel("Cross-Shore Distance (dam)")
        plt.ylim(bottom=OriginY - 35)
        plt.tight_layout()
        timestr = "Time = " + str(t) + " yrs"
        plt.text(1, OriginY - 33, timestr)
        plt.rcParams.update({"font.size": 16})
        name = "elev_" + str(t)
        # elevFig2.tight_layout()
        elevFig2.savefig(name)  # dpi=200
        plt.close(elevFig2)

    frames = []

    for filenum in range(TMAX_SIM):
        filename = "elev_" + str(filenum) + ".png"
        frames.append(imageio.imread(filename))

        if (
            filenum < TMAX
            and cascade.nourishments[iB3D]._post_storm_interior[TMAX] is not None
        ):
            filename = "elev_" + str(filenum) + "pt5" ".png"
            frames.append(imageio.imread(filename))

    imageio.mimsave("elev.gif", frames, "GIF-FI")
    print()
    print("[ * GIF successfully generated * ]")


# ===================================================
# 2: Shoreline positions over time (#6 in Barrier3D_Functions)
#
#       inputs:         - x_s_TS, x_b_TS (shoreline and back-barrier time series for one B3D subgrid)
#       outputs:        - fig


def plot_ShorelinePositions(x_s_TS, x_b_TS):
    scts = [(x - x_s_TS[0]) * -10 for x in x_s_TS]
    bscts = [(x - x_s_TS[0]) * -10 for x in x_b_TS]
    plt.figure()
    plt.plot(scts, "b")
    plt.plot(bscts, "g")
    fig = plt.gcf()
    fig.set_size_inches(14, 8)
    plt.ylabel("Shoreline Position (m)")
    plt.xlabel("Year")
    plt.show()


# ===================================================
# 3 B3D cross-shore transect for one subgrid every 100 m for last time step (#5 in Barrier3D_Functions)
#
#       inputs:         - barrier3d (a list containing BMI objects)
#                       - TMAX (the last time index for plotting)
#       outputs:        - fig


def plot_XShoreTransects(barrier3d, TMAX):
    # Build beach elevation
    BW = 6  # beach width (dam) for illustration purposes
    BeachX = np.zeros(BW)
    berm = math.ceil(BW * 0.5)
    BeachX[berm : BW + 1] = barrier3d._BermEl
    add = (barrier3d._BermEl - barrier3d._SL) / berm
    for i in range(berm):
        BeachX[i] = barrier3d._SL + add * i
        # Plot full transects
    plt.figure()

    for v in range(0, barrier3d._BarrierLength, 20):
        CrossElev = barrier3d._InteriorDomain[:, v]
        Dunes = barrier3d._DuneDomain[TMAX, v, :] + barrier3d._BermEl
        CrossElev1 = np.insert(CrossElev, 0, Dunes)
        CrossElev2 = np.insert(CrossElev1, 0, BeachX)
        CrossElev = CrossElev2 * 10  # Convert to meters
        plt.plot(CrossElev)
    fig = plt.gcf()
    fig.set_size_inches(14, 6)
    plt.hlines(barrier3d._SL, -1, len(CrossElev + 1), colors="dodgerblue")
    plt.xlabel("Cross-Shore Distance (dam)")
    plt.ylabel("Elevation (m)")
    plt.title("Cross-shore Topo Transects")
    plt.show()
    name = "Output/Profiles"
    # fig.savefig(name)


# ===================================================
# 4: Cross-shore transects for both brie and B3d; for second function, just cascade (b3d) which includes beach width
#
#       inputs:         - barrier3d (a list of B3D objects)
#                       - brieLTA (a brie object with LTA model on)
#                       - time_step (a list of time steps to plot)
#                       - iB3D (an integer corresponding to the B3D subdomain / brie alongshore grid)
#       outputs:        - fig


def plot_ModelTransects_b3d_brie(b3d, brieLTA, time_step, iB3D):
    plt.figure(figsize=(10, 5))
    colors = mpl.cm.viridis(np.linspace(0, 1, b3d[0]._TMAX))

    plt.subplot(2, 1, 1)

    for t in time_step:

        # Sea level
        SL = b3d[iB3D]._SL + (t * b3d[iB3D]._RSLR[t])

        # Create data points
        Tx = b3d[iB3D]._x_t_TS[t] - b3d[iB3D]._x_t_TS[0]
        Ty = (SL - b3d[iB3D]._DShoreface) * 10
        Sx = b3d[iB3D]._x_s_TS[t] - b3d[iB3D]._x_t_TS[0]
        Sy = SL * 10
        Bx = b3d[iB3D]._x_b_TS[t] - b3d[iB3D]._x_t_TS[0]
        By = (SL - b3d[iB3D]._BayDepth) * 10
        My = By

        BW = 6  # beach width (dam) for illustration purposes
        BeachX = np.zeros(BW)
        berm = math.ceil(BW * 0.5)
        BeachX[berm : BW + 1] = b3d[iB3D]._BermEl
        add = (b3d[iB3D]._BermEl - b3d[iB3D]._SL) / berm
        for i in range(berm):
            BeachX[i] = b3d[iB3D]._SL + add * i

        v = 0  # just use the first transect
        CrossElev = b3d[iB3D]._DomainTS[t]
        CrossElev = CrossElev[:, v]
        Dunes = b3d[iB3D]._DuneDomain[t, v, :] + b3d[iB3D]._BermEl
        CrossElev1 = np.insert(CrossElev, 0, Dunes)
        CrossElev2 = np.insert(CrossElev1, 0, BeachX)
        CrossElev = (CrossElev2 * 10) + Sy  # Convert to meters
        xCrossElev = np.arange(0, len(CrossElev), 1) + Sx

        Mx = xCrossElev[-1] + 20  # just add a buffer to the end of the plt

        x = np.hstack([Tx, Sx, xCrossElev, Mx])
        y = np.hstack([Ty, Sy, CrossElev, My])

        # Plot
        plt.plot(x, y, color=colors[t])
        plt.hlines(SL * 10, Tx, Mx, colors="dodgerblue")
        # plt.hlines((b3d[iB3D]._SL + (t * b3d[iB3D]._RSLR[t])) * 10, Tx, Sx, colors='dodgerblue')
        # plt.hlines((b3d[iB3D]._SL + (t * b3d[iB3D]._RSLR[t])) * 10, Bx + 2, xCrossElev[-1]+20, colors='dodgerblue')  # KA: scrappy fix
        # plt.xlim([0, 500])
        plt.rcParams.update({"font.size": 20})

        plt.ylabel("Elevation (m)")
        plt.title("Profile Evolution")
        plt.show()

    plt.subplot(2, 1, 2)

    for t in time_step:
        # BRIE

        # Sea level
        SL = b3d[iB3D]._SL + (t * b3d[iB3D]._RSLR[t])

        # Create data points
        Tx = (
            brieLTA._x_t_save[iB3D, t] - brieLTA._x_t_save[iB3D, 0]
        ) / 10  # convert to dam and subtract start position
        Ty = (SL * 10) - brieLTA._d_sf  # units of m
        Sx = (brieLTA._x_s_save[iB3D, t] - brieLTA._x_t_save[iB3D, 0]) / 10
        Sy = SL * 10
        Bx = (brieLTA._x_b_save[iB3D, t] - brieLTA._x_t_save[iB3D, 0]) / 10
        By = (SL * 10) - brieLTA._bb_depth
        Hx1 = Sx
        Hy1 = brieLTA._h_b_save[iB3D, t] + Sy
        Hx2 = Bx
        Hy2 = brieLTA._h_b_save[iB3D, t] + Sy
        Mx = Bx + 20
        My = By

        x = [Tx, Sx, Hx1, Hx2, Bx, Mx]
        y = [Ty, Sy, Hy1, Hy2, By, My]

        # Plot
        plt.plot(x, y, color=colors[t])
        # plt.xlim([0, 500])
        plt.rcParams.update({"font.size": 20})
        plt.hlines(SL * 10, Tx, Mx, colors="dodgerblue")
        # plt.hlines((b3d[iB3D]._SL + (t * b3d[iB3D]._RSLR[t])) * 10, Tx, Sx, colors='dodgerblue')
        # plt.hlines((b3d[iB3D]._SL + (t * b3d[iB3D]._RSLR[t])) * 10, Bx, Mx, colors='dodgerblue')

        plt.xlabel("Alongshore Distance (dam)")
        plt.ylabel("Elevation (m)")
        plt.show()


def plot_ModelTransects(cascade, time_step, iB3D):
    plt.figure(figsize=(10, 5))
    fig = plt.subplot(1, 1, 1)
    legend_t = []

    for t in time_step:

        # Sea level
        sea_level = cascade.barrier3d[iB3D]._SL + (t * cascade.barrier3d[iB3D]._RSLR[t])

        # Create data points
        shoreface_toe_x = (
            cascade.barrier3d[iB3D].x_t_TS[t] - cascade.barrier3d[iB3D].x_t_TS[0]
        )
        shoreface_toe_y = (sea_level - cascade.barrier3d[iB3D].DShoreface) * 10  # m
        shoreline_x = (
            cascade.barrier3d[iB3D].x_s_TS[t] - cascade.barrier3d[iB3D].x_t_TS[0]
        )
        shoreline_y = sea_level * 10  # m
        bay_y = (sea_level - cascade.barrier3d[iB3D]._BayDepth) * 10  # m
        end_of_bay_y = bay_y

        if cascade.nourishments[iB3D].beach_width[t] is not None:
            berm_x = shoreline_x + (
                cascade.nourishments[iB3D].beach_width[t] / 10
            )  # beach width (in dam)
        else:
            berm_x = shoreline_x + (
                int(cascade.barrier3d[iB3D].BermEl / cascade.barrier3d[iB3D]._beta)
            )  # initial beach width (in dam)
        berm_y = (
            cascade.barrier3d[iB3D]._BermEl * 10
        ) + shoreline_y  # convert to meters
        dune_toe_x = berm_x
        dune_toe_y = berm_y

        v = 10  # just use 10th first transect
        interior_y = cascade.barrier3d[iB3D]._DomainTS[t]
        interior_y = interior_y[:, v]
        dunes_y = (
            cascade.barrier3d[iB3D]._DuneDomain[t, v, :]
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
        plt.plot(x, y)
        plt.hlines(sea_level * 10, shoreface_toe_x, end_of_bay_x, colors="black")
        # NOTE: the berm elevation is relative to the MHW, so everything that relies on it is m MHW; confirmed with Ian
        # that the InteriorDomain is m MHW (i.e., m NAVD88 - MHW [in NAVD88])
        plt.rcParams.update({"font.size": 20})
        legend_t.append(str(t))

    plt.xlabel("Cross-shore position (dam)")
    plt.ylabel("Elevation (m MHW)")
    plt.title("Profile Evolution")
    plt.legend(legend_t)

    return fig


# ===================================================
# 5: Statistics from B3D
#
#       inputs:         - b3d (a list of B3D objects)
#                       - iB3D (an integer corresponding to the B3D subdomain / brie alongshore grid)
#                       - TMAX (the last time index for plotting)
#                       - iB3D
#       outputs:        - fig


def plot_statistics(b3d, iB3D, TMAX):
    colors = mpl.cm.viridis(np.linspace(0, 1, 10))
    plt.figure(figsize=(10, 10))

    # A Shoreface Slope
    plt.subplot(3, 2, 1)
    ssfTS = b3d[iB3D]._s_sf_TS
    plt.plot(ssfTS, color=colors[1])
    plt.hlines(b3d[iB3D]._s_sf_eq, 0, TMAX, colors="black", linestyles="dashed")
    # plt.ylabel(r'$\alpha$')
    plt.ylabel("Shoreface Slope")
    plt.legend(["B3D sub-grid #" + str(iB3D)])
    plt.rcParams["legend.loc"] = "lower right"

    # A Interior Width
    plt.subplot(3, 2, 2)
    aiw = [a * 10 for a in b3d[iB3D]._InteriorWidth_AvgTS]
    plt.plot(aiw, color=colors[2])
    plt.ylabel("Avg. Width (m)")  # Average Interior Width
    plt.legend(["B3D sub-grid #" + str(iB3D)])
    plt.rcParams["legend.loc"] = "lower right"

    # A Shoreline Change
    scts = [(x - b3d[iB3D]._x_s_TS[0]) * 10 for x in b3d[iB3D]._x_s_TS]
    plt.subplot(3, 2, 3)
    plt.plot(scts, color=colors[3])
    plt.ylabel("Shoreline Position (m)")
    plt.legend(["B3D sub-grid #" + str(iB3D)])
    plt.rcParams["legend.loc"] = "lower right"

    # A Dune Height
    aHd = [a * 10 for a in b3d[iB3D]._Hd_AverageTS]
    plt.subplot(3, 2, 4)
    plt.plot(aHd, color=colors[4])
    plt.ylabel("Avg. Dune Height (m)")  # Average Dune Height
    plt.legend(["B3D sub-grid #" + str(iB3D)])
    plt.rcParams["legend.loc"] = "lower right"

    # Qoverwash for entire B3D grid
    plt.subplot(3, 2, 5)
    Qoverwash = np.zeros(np.size(b3d[0]._QowTS[0:TMAX]))
    for iGrid in range(len(b3d)):
        Qoverwash = Qoverwash + (
            np.array(b3d[iGrid]._QowTS[0:TMAX]) * (b3d[iGrid]._BarrierLength * 10)
        )  # m^3/yr
    Qoverwash = Qoverwash / (len(b3d) * b3d[0]._BarrierLength * 10)  # m^3/m/yr
    plt.plot(Qoverwash, color=colors[5])
    plt.ylabel(r"Qow ($m^3/m/yr$)")
    plt.xlabel("Year")

    # Shoreface flux for entire B3D grid
    plt.subplot(3, 2, 6)
    Qsf = np.zeros(np.size(b3d[0]._QsfTS[0:TMAX]))
    for iGrid in range(len(b3d)):
        Qsf = Qsf + (
            np.array(b3d[iGrid]._QsfTS[0:TMAX]) * (b3d[iGrid]._BarrierLength * 10)
        )  # m^3/yr
    Qsf = Qsf / (len(b3d) * b3d[0]._BarrierLength * 10)  # m^3/m/yr
    plt.plot(Qsf, color=colors[6])
    plt.ylabel(r"Qsf ($m^3/m/yr$)")
    plt.xlabel("Year")

    plt.show()


def plot_statisticsAGU(b3d, brieLTA, iB3D, iBRIE, TMAX, TMAX_BRIE):
    # make a pretty colormap
    # my_cmap = sns.color_palette("husl", 9, as_cmap=True)
    # my_cmap = sns.color_palette("flare", as_cmap=True)
    # my_cmap = sns.color_palette
    # my_cmap = sns.color_palette("Set2", 9, as_cmap=True)
    iC1 = 0
    iC2 = 5
    iC3 = 9
    # colors = my_cmap(np.linspace(0, 1, 10))
    colors = mpl.cm.vidiris(np.linspace(0, 1, 10))
    fig = plt.figure(figsize=(16, 8))

    # A Shoreface Slope
    plt.subplot(2, 3, 1)
    ssfTS = b3d[iB3D]._s_sf_TS
    plt.plot(ssfTS, color=colors[iC1])
    plt.hlines(b3d[iB3D]._s_sf_eq, 0, TMAX, colors="black", linestyles="dashed")

    ssfTS = brieLTA._s_sf_save[iBRIE, :]
    plt.plot(ssfTS, color=colors[iC2], linewidth=2)
    plt.hlines(
        brieLTA._s_sf_eq,
        0,
        len(brieLTA._s_sf_save[iBRIE, :]),
        colors="black",
        linestyles="dashed",
    )

    plt.ylabel("Shoreface Slope")
    plt.xlabel("Year")
    plt.xlim([0, 1000])
    plt.legend(["CASCADE", "BRIE (no coupling)"])
    plt.rcParams["legend.loc"] = "lower right"
    plt.rcParams.update({"font.size": 15})

    # A Interior Width
    plt.subplot(2, 3, 2)
    aiw = [a * 10 for a in b3d[iB3D]._InteriorWidth_AvgTS]
    plt.plot(aiw, color=colors[iC1])

    aiw = brieLTA._x_b_save[iBRIE, :] - brieLTA._x_s_save[iBRIE, :]
    plt.plot(aiw, color=colors[iC2], linewidth=3)
    plt.rcParams["legend.loc"] = "lower right"
    plt.xlabel("Year")
    plt.xlim([0, 1000])
    plt.ylabel("Barrier Width (m)")
    plt.rcParams.update({"font.size": 15})

    # A Shoreline Change
    plt.subplot(2, 3, 3)
    scts = [(x - b3d[iB3D]._x_s_TS[0]) * 10 for x in b3d[iB3D]._x_s_TS]
    plt.plot(scts, color=colors[iC1], linewidth=2)

    scts = [(x - brieLTA._x_s_save[iBRIE, 0]) for x in brieLTA._x_s_save[iBRIE, :]]
    plt.plot(scts, color=colors[iC2], linewidth=2)
    plt.xlabel("Year")
    plt.xlim([0, 1000])
    plt.ylim([0, 500])
    plt.ylabel("Shoreline Position (m)")
    plt.rcParams["legend.loc"] = "lower right"
    plt.rcParams.update({"font.size": 15})

    # Qoverwash for entire B3D grid
    plt.subplot(2, 3, 4)
    Qoverwash = np.zeros(np.size(b3d[0]._QowTS[0:TMAX]))
    for iGrid in range(len(b3d)):
        Qoverwash = Qoverwash + (
            np.array(b3d[iGrid]._QowTS[0:TMAX]) * (b3d[iGrid]._BarrierLength * 10)
        )  # m^3/yr
    Qoverwash = Qoverwash / (len(b3d) * b3d[0]._BarrierLength * 10)  # m^3/m/yr
    plt.plot(Qoverwash, color=colors[iC1])
    # movingavg = np.convolve(Qoverwash, np.ones((50,))/50, mode='valid')
    # movingavg = [i * 10 for i in movingavg]
    # plt.plot(movingavg, 'r--')
    df = pd.DataFrame(data=Qoverwash)
    plt.plot(df.rolling(window=50).mean(), color=colors[iC3])

    QoverwashLTA = brieLTA._Qoverwash[0:TMAX_BRIE] / (
        brieLTA._ny * brieLTA._dy
    )  # from brie in m^3/yr --> m^3/m/yr
    plt.plot(brieLTA._t[0:TMAX_BRIE], QoverwashLTA, color=colors[iC2], linewidth=2)
    plt.ylabel(r"Qow ($m^3/m$)")
    plt.xlabel("Year")
    plt.xlim([0, 1000])
    plt.rcParams.update({"font.size": 15})

    # shoreface flux for entire grid
    plt.subplot(2, 3, 5)
    Qsf = np.zeros(np.size(b3d[0]._QsfTS[0:TMAX]))
    for iGrid in range(len(b3d)):
        Qsf = Qsf + (
            np.array(b3d[iGrid]._QsfTS[0:TMAX]) * (b3d[iGrid]._BarrierLength * 10)
        )  # m^3/yr
    Qsf = Qsf / (len(b3d) * b3d[0]._BarrierLength * 10)  # m^3/m/yr
    plt.plot(Qsf, color=colors[iC1])
    df = pd.DataFrame(data=Qsf)
    plt.plot(df.rolling(window=50).mean(), color=colors[iC3])

    QsfLTA = brieLTA._Qshoreface[0:TMAX_BRIE] / (
        brieLTA._ny * brieLTA._dy
    )  # from brie in m^3/yr --> m^3/m/yr
    plt.plot(brieLTA._t[0:TMAX_BRIE], QsfLTA, color=colors[iC2], linewidth=2)
    plt.ylabel(r"Qsf ($m^3/m$)")
    plt.xlabel("Year")
    plt.rcParams.update({"font.size": 15})
    plt.xlim([0, 1000])

    # A Dune Height
    aHd = [a * 10 for a in b3d[iB3D]._Hd_AverageTS]
    plt.subplot(2, 3, 6)
    plt.plot(aHd, color=colors[iC1])
    plt.ylabel("Avg. Dune Height (m)")  # Average Dune Height
    plt.xlim([0, 1000])
    plt.xlabel("Year")
    plt.rcParams.update({"font.size": 15})

    plt.show()
    fig.tight_layout()


# ===================================================
# 6: Statistics from BRIE with LTA
#
#       inputs:         - brieLTA (brie model with LTA on)
#                       - iB3D (an integer corresponding to the B3D subdomain / brie alongshore grid)
#                       - TMAX (the last time index for plotting and calculations)
#       outputs:        - fig


def plot_statistics_BRIE(brieLTA, iB3D, TMAX):
    colors = mpl.cm.viridis(np.linspace(0, 1, 10))
    plt.figure(figsize=(10, 10))

    # A Shoreface Slope
    plt.subplot(3, 2, 1)
    ssfTS = brieLTA._s_sf_save[iB3D, :]
    plt.plot(ssfTS, color=colors[1])
    plt.hlines(
        brieLTA._s_sf_eq,
        0,
        len(brieLTA._s_sf_save[iB3D, :]),
        colors="black",
        linestyles="dashed",
    )
    # plt.ylabel(r'$\alpha$')
    plt.ylabel("Shoreface Slope")
    plt.legend(["BRIE sub-grid #" + str(iB3D)])
    plt.rcParams["legend.loc"] = "lower right"

    # A Interior Width
    plt.subplot(3, 2, 2)
    aiw = brieLTA._x_b_save[iB3D, :] - brieLTA._x_s_save[iB3D, :]
    plt.plot(aiw, color=colors[2])
    plt.ylabel("Barrier Width (m)")
    plt.legend(["BRIE sub-grid #" + str(iB3D)])
    plt.rcParams["legend.loc"] = "lower right"

    # A Shoreline Change
    scts = [(x - brieLTA._x_s_save[iB3D, 0]) for x in brieLTA._x_s_save[iB3D, :]]
    plt.subplot(3, 2, 3)
    plt.plot(scts, color=colors[3])
    plt.ylabel("Shoreline Position (m)")
    plt.legend(["BRIE sub-grid #" + str(iB3D)])
    plt.rcParams["legend.loc"] = "lower right"

    # Barrier Height
    aHd = brieLTA._h_b_save[iB3D, :]
    plt.subplot(3, 2, 4)
    plt.plot(aHd, color=colors[4])
    plt.ylabel("Barrier Height (m)")  # Average Dune Height
    plt.legend(["BRIE sub-grid #" + str(iB3D)])
    plt.rcParams["legend.loc"] = "upper right"
    plt.xlabel("Year")

    # Qoverwash for entire BRIE grid
    plt.subplot(3, 2, 5)
    QoverwashLTA = brieLTA._Qoverwash[0:TMAX] / (
        brieLTA._ny * brieLTA._dy
    )  # from brie in m^3/yr --> m^3/m/yr
    plt.plot(brieLTA._t[0:TMAX], QoverwashLTA, color=colors[5])
    plt.ylabel(r"Qow ($m^3/m/yr$)")
    plt.xlabel("Year")

    # Shoreface flux for entire B3D grid
    plt.subplot(3, 2, 6)
    QsfLTA = brieLTA._Qshoreface[0:TMAX] / (
        brieLTA._ny * brieLTA._dy
    )  # from brie in m^3/yr --> m^3/m/yr
    plt.plot(QsfLTA, color=colors[6])
    plt.ylabel(r"Qsf ($m^3/m/yr$)")
    plt.xlabel("Year")

    plt.show()


# ===================================================
# 7: Calculate shoreline change periodicity
#
#       inputs:         - x_s_TS (an single shoreline time series for one B3D subdomain)
#       outputs:        - Periodicity (period of punctuated retreat)
#                       - AvgFastDur (average duration of fast periods)
#                       - AvgSlowDur (average duration of slow periods)
#                       - Punc (boolean for punctuated retreat)


def calc_ShorelinePeriodicity(x_s_TS):
    # Shoreline Change & Change Rate Over Time
    scts = [(x - x_s_TS[0]) * 10 for x in x_s_TS]

    # Filter
    win = 31
    poly = 3
    der1 = signal.savgol_filter(scts, win, poly, deriv=1)

    HitDown = []  # Slow-downs
    HitUp = []  # Speed-ups

    window1 = 3  # Maximum allowed length for gaps in slow periods
    window2 = 30  # Minimum length required for slow periods, including gaps
    buffer = 3
    thresh1 = 0.5  # Max slope for slow periods (if shoreline change rate is below this, immobile, otherwise transgressive)
    thresh2 = 1

    # Find slow periods
    der_under = np.where(der1 < thresh1)[0]

    if len(der_under) > 0:

        gaps = np.diff(der_under) > window1
        peak_start = np.insert(der_under[1:][gaps], 0, der_under[0])
        peak_stop = np.append(der_under[:-1][gaps], der_under[-1])

        for n in range(len(peak_stop)):
            if peak_stop[n] - peak_start[n] > window2:
                if len(HitDown) == 0:
                    if peak_start[n] > buffer:
                        HitDown.append(peak_start[n])
                    if peak_stop[n] < len(scts) - buffer:
                        HitUp.append(peak_stop[n])
                else:
                    gap_length = peak_start[n] - HitUp[-1]
                    gap_slope = (scts[peak_start[n]] - scts[HitUp[-1]]) / gap_length
                    if gap_length < window2 and gap_slope < thresh2:
                        if peak_stop[n] < len(scts) - buffer:
                            HitUp[-1] = peak_stop[n]
                        else:
                            del HitUp[-1]
                    else:
                        if peak_start[n] > buffer:
                            HitDown.append(peak_start[n])
                        if peak_stop[n] < len(scts) - buffer:
                            HitUp.append(peak_stop[n])

    ### CALCULATE STATS

    Jumps = len(HitDown)
    Slows = len(HitUp)
    SlowDur = []
    FastDur = []

    if Jumps > 0 and Slows > 0:
        DownFirst = HitDown[0] < HitUp[0]
    elif Jumps == 0 and Slows > 1:
        DownFirst = True
    else:
        DownFirst = False

    if Jumps >= 2 or Slows >= 2:
        if Jumps >= 2 and Slows >= 2:
            Periodicity = (np.mean(np.diff(HitDown)) + np.mean(np.diff(HitUp))) / 2
        elif Jumps >= Slows:
            Periodicity = np.mean(np.diff(HitDown))
        else:
            Periodicity = np.mean(np.diff(HitUp))
        if DownFirst:
            for n in range(Slows):
                SlowDur.append(HitUp[n] - HitDown[n])
            for n in range(Jumps - 1):
                FastDur.append(HitDown[n + 1] - HitUp[n])
        else:
            for n in range(Slows - 1):
                SlowDur.append(HitUp[n + 1] - HitDown[n])
            for n in range(Jumps):
                FastDur.append(HitDown[n] - HitUp[n])
    else:
        Periodicity = 0
        if Jumps == 1 and Slows == 1:
            if DownFirst:
                SlowDur.append(HitUp[0] - HitDown[0])
            else:
                FastDur.append(HitDown[0] - HitUp[0])

    AvgFastDur = np.mean(FastDur)
    if np.isnan(AvgFastDur):
        AvgFastDur = 0
    AvgSlowDur = np.mean(SlowDur)
    if np.isnan(AvgSlowDur):
        AvgSlowDur = 0

    if len(SlowDur) >= 2 and len(FastDur) >= 2:
        Punc = 1
    else:
        Punc = 0

    return Periodicity, AvgFastDur, AvgSlowDur, Punc


# ===================================================
# 8: Dune Height Over Time for CASCADE
#
#       inputs:         - barrier3d (a list containing BMI objects)
#                       - TMAX (the last time index for plotting)
#       outputs:        - fig (dune domain for all sub-grids)


def plot_dune_domain(b3d, TMAX):
    DuneCrest = []
    Dmax = []

    for iB3D in range(len(b3d)):
        sub_domain = b3d[iB3D]._DuneDomain[0:TMAX, :, :]
        DuneCrest.append(sub_domain.max(axis=2))
        Dmax.append(b3d[iB3D]._Dmax)

    DuneCrest = np.hstack(DuneCrest).astype(float)
    Dmax = np.max(Dmax)

    duneFig = plt.figure(figsize=(14, 8))
    plt.rcParams.update({"font.size": 13})
    ax = duneFig.add_subplot(111)
    ax.matshow(
        (DuneCrest) * 10,
        origin="lower",
        cmap="bwr",
        aspect="auto",
        vmin=0,
        vmax=Dmax * 10,
    )
    cax = ax.xaxis.set_ticks_position("bottom")  # analysis:ignore
    # cbar = duneFig.colorbar(cax)
    # cbar.set_label('Dune Height Above Berm Elevation (m)', rotation=270)
    plt.xlabel("Alongshore Distance (dam)")
    plt.ylabel("Year")
    plt.title("Dune Height (m)")


# ===================================================
# 9: Average shoreline change rate over time for CASCADE
#
#       inputs:         - b3d (a list containing BMI objects)
#                       - TMAX (the last time index for plotting)
#       outputs:        - fig (dune domain for all sub-grids)


def plot_ShorelineChangeRate(b3d):
    ave_rate = []

    for iB3D in range(len(b3d)):
        scts = [(x - b3d[iB3D]._x_s_TS[0]) * 10 for x in b3d[iB3D]._x_s_TS]
        rate = [0]
        for k in range(1, len(scts)):
            rate.append(scts[k] - scts[k - 1])
        ave_rate.append(rate)

    df = pd.DataFrame(data=ave_rate)
    ave = df.mean(axis=0)

    ratefig = plt.figure()
    plt.plot(ave)
    fig = plt.gcf()
    fig.set_size_inches(14, 5)
    plt.xlabel("Year")
    plt.ylabel("Shoreline Erosion Rate (m/yr)")
    plt.show()


def plot_ShorelineChangeRate_AGU(b3d1, b3d2):
    plt.figure(figsize=(14, 5))

    ave_rate1 = []
    ave_scts1 = []

    for iB3D in range(len(b3d1)):
        scts = [(x - b3d1[iB3D]._x_s_TS[0]) * 10 for x in b3d1[iB3D]._x_s_TS]
        rate = [0]
        for k in range(1, len(scts)):
            rate.append(scts[k] - scts[k - 1])
        ave_rate1.append(rate)
        ave_scts1.append(scts)

    df = pd.DataFrame(data=ave_rate1)
    ave1 = df.mean(axis=0)
    df = pd.DataFrame(data=ave_scts1)
    aveTS1 = df.mean(axis=0)

    ave_rate2 = []
    ave_scts2 = []

    for iB3D in range(len(b3d2)):
        scts = [(x - b3d2[iB3D]._x_s_TS[0]) * 10 for x in b3d2[iB3D]._x_s_TS]
        rate = [0]
        for k in range(1, len(scts)):
            rate.append(scts[k] - scts[k - 1])
        ave_rate2.append(rate)
        ave_scts2.append(scts)

    df = pd.DataFrame(data=ave_rate2)
    ave2 = df.mean(axis=0)
    df = pd.DataFrame(data=ave_scts2)
    aveTS2 = df.mean(axis=0)

    plt.subplot(1, 3, 1)
    plt.plot(aveTS1)
    plt.plot(aveTS2)
    plt.xlabel("Year")
    plt.ylabel("Average Shoreline Position (m)")
    plt.legend(["homogenous dune growth", "alongshore variable dune growth"])
    plt.rcParams["legend.loc"] = "best"

    plt.subplot(1, 3, 2)
    plt.plot(ave1)
    plt.plot(ave2)
    plt.xlabel("Year")
    plt.ylabel("Average Shoreline Erosion Rate (m/yr)")
    plt.rcParams.update({"font.size": 15})

    # # Qoverwash for entire B3D grid
    # plt.subplot(1, 3, 3)
    #
    # Qoverwash = np.zeros(np.size(b3d[0]._QowTS[0:TMAX]))
    # for iGrid in range(len(b3d)):
    #     Qoverwash = Qoverwash + (np.array(b3d[iGrid]._QowTS[0:TMAX]) * (b3d[iGrid]._BarrierLength * 10)) # m^3/yr
    # Qoverwash = Qoverwash / (len(b3d) * b3d[0]._BarrierLength * 10) # m^3/m/yr
    # plt.plot(Qoverwash)
    # plt.ylabel(r'Qow ($m^3/m/yr$)')
    # plt.xlabel('Year')

    plt.show()


# ===================================================
# 10: Differences in punctuated retreat for CASCADE vs B3D (AST model vs no AST)
#
#       inputs:         - CASCADE_b3d (a list containing BMI objects from a CASCADE run)
#                       - b3d_only (a list containing BMI objects from individual b3d runs, corresponding to the growth
#                       rates used in the CASCADE model above)
#                       - ny (number of alongshore cells)
#       outputs:        - fig (statistical comparison of punctuated retreat)


def plot_punctuated_difference(CASCADE_b3d, b3d_only, ny):
    # Calculate shoreline change periodicity from CASCADE model
    Punc = []
    Period = []
    AvgFastDur = []
    AvgSlowDur = []
    ShorelinePosition = []
    sc_rate = []
    var_sc_rate = []
    mean_sc_rate = []

    for iB3D in range(ny):
        tmpPeriod, tmpAvgFastDur, tmpAvgSlowDur, tmpPunc = calc_ShorelinePeriodicity(
            CASCADE_b3d[iB3D]._x_s_TS
        )
        Punc.append(tmpPunc)
        Period.append(tmpPeriod)
        AvgFastDur.append(tmpAvgFastDur)
        AvgSlowDur.append(tmpAvgSlowDur)
        ShorelinePosition.append(CASCADE_b3d[iB3D]._x_s_TS)

        # shoreline change rate
        scts = [
            (x - CASCADE_b3d[iB3D]._x_s_TS[0]) * 10 for x in CASCADE_b3d[iB3D]._x_s_TS
        ]  # m/yr
        rate = [0]
        for k in range(1, len(scts)):
            rate.append(scts[k] - scts[k - 1])
        sc_rate.append(rate)
        var_sc_rate.append(np.var(rate))  # variance of shoreline change rate
        mean_sc_rate.append(np.mean(rate))  # mean of shoreline change rate

    # Calculate shoreline change periodicity from b3d only model
    Punc_b3D = []
    Period_b3D = []
    AvgFastDur_b3D = []
    AvgSlowDur_b3D = []
    ShorelinePosition_b3D = []
    sc_rate_b3d = []
    var_sc_rate_b3d = []
    mean_sc_rate_b3d = []

    for iB3D in range(ny):
        tmpPeriod, tmpAvgFastDur, tmpAvgSlowDur, tmpPunc = calc_ShorelinePeriodicity(
            b3d_only[iB3D]._x_s_TS
        )
        Punc_b3D.append(tmpPunc)
        Period_b3D.append(tmpPeriod)
        AvgFastDur_b3D.append(tmpAvgFastDur)
        AvgSlowDur_b3D.append(tmpAvgSlowDur)
        ShorelinePosition_b3D.append(b3d_only[iB3D]._x_s_TS)

        # shoreline change rate
        scts = [
            (x - b3d_only[iB3D]._x_s_TS[0]) * 10 for x in b3d_only[iB3D]._x_s_TS
        ]  # m/yr
        rate = [0]
        for k in range(1, len(scts)):
            rate.append(scts[k] - scts[k - 1])
        sc_rate_b3d.append(rate)
        var_sc_rate_b3d.append(np.var(rate))  # variance of shoreline change rate
        mean_sc_rate_b3d.append(np.mean(rate))  # mean of shoreline change rate

    colors = mpl.cm.viridis(np.linspace(0, 1, 10))
    # plt.figure(figsize=(10, 8))
    plt.figure(figsize=(17, 8))

    # punctuated retreat
    plt.subplot(2, 4, 1)
    plt.scatter(
        np.arange(500, 500 * (ny + 1), 500), Punc, color=colors[1], marker="p", s=500
    )
    plt.scatter(
        np.arange(500, 500 * (ny + 1), 500),
        Punc_b3D,
        color=colors[6],
        marker="o",
        s=200,
    )
    plt.ylabel("Punctuated Retreat")
    plt.xlabel("alongshore (m)")
    plt.legend(["CASCADE (AST)", "B3D (no AST)"])
    plt.rcParams["legend.loc"] = "center"

    # period
    plt.subplot(2, 4, 2)
    plt.scatter(
        np.arange(500, 500 * (ny + 1), 500), Period, color=colors[1], marker="p", s=500
    )
    plt.scatter(
        np.arange(500, 500 * (ny + 1), 500),
        Period_b3D,
        color=colors[6],
        marker="o",
        s=200,
    )
    plt.ylabel("Period (years)")
    plt.xlabel("alongshore (m)")

    # fast duration
    plt.subplot(2, 4, 3)
    plt.scatter(
        np.arange(500, 500 * (ny + 1), 500),
        AvgFastDur,
        color=colors[1],
        marker="p",
        s=500,
    )
    plt.scatter(
        np.arange(500, 500 * (ny + 1), 500),
        AvgFastDur_b3D,
        color=colors[6],
        marker="o",
        s=200,
    )
    plt.ylabel("Ave Fast Duration (years)")
    plt.xlabel("alongshore (m)")

    # slow duration
    plt.subplot(2, 4, 4)
    plt.scatter(
        np.arange(500, 500 * (ny + 1), 500),
        AvgSlowDur,
        color=colors[1],
        marker="p",
        s=500,
    )
    plt.scatter(
        np.arange(500, 500 * (ny + 1), 500),
        AvgSlowDur_b3D,
        color=colors[6],
        marker="o",
        s=200,
    )
    plt.ylabel("Ave Slow Duration (years)")
    plt.xlabel("alongshore (m)")

    # shoreline change rate variance
    plt.subplot(2, 4, 5)
    plt.scatter(
        np.arange(500, 500 * (ny + 1), 500),
        var_sc_rate,
        color=colors[1],
        marker="p",
        s=500,
    )
    plt.scatter(
        np.arange(500, 500 * (ny + 1), 500),
        var_sc_rate_b3d,
        color=colors[6],
        marker="o",
        s=200,
    )
    plt.ylabel("Variance SCR (m/yr)")
    plt.xlabel("alongshore (m)")

    # shoreline change rate mean
    plt.subplot(2, 4, 6)
    plt.scatter(
        np.arange(500, 500 * (ny + 1), 500),
        mean_sc_rate,
        color=colors[1],
        marker="p",
        s=500,
    )
    plt.scatter(
        np.arange(500, 500 * (ny + 1), 500),
        mean_sc_rate_b3d,
        color=colors[6],
        marker="o",
        s=200,
    )
    plt.ylabel("Mean SCR (m/yr)")
    plt.xlabel("alongshore (m)")

    plt.tight_layout()

    # figure of shoreline change plots
    plt.figure(figsize=(17, 8))

    for iB3D in range(ny):
        plt.subplot(2, 6, iB3D + 1)
        plt.plot(ShorelinePosition[iB3D], color=colors[1])
        plt.plot(ShorelinePosition_b3D[iB3D], color=colors[6])
        plt.xlabel("Year")
        plt.ylabel("Average Shoreline Position (m)")
    plt.legend(["CASCADE", "B3D only"])
    plt.rcParams["legend.loc"] = "upper left"
    plt.tight_layout()

    # figure of shoreline change rate
    plt.figure(figsize=(17, 8))

    for iB3D in range(ny):
        plt.subplot(2, 6, iB3D + 1)
        plt.plot(sc_rate[iB3D], color=colors[1])
        plt.plot(sc_rate_b3d[iB3D], color=colors[6])
        plt.xlabel("Year")
        plt.ylabel("Shoreline change rate (m/yr)")
    plt.legend(["CASCADE", "B3D only"])
    plt.rcParams["legend.loc"] = "upper left"
    plt.tight_layout()

    # figure of histograms
    plt.figure(figsize=(17, 8))

    for iB3D in range(ny):
        plt.subplot(2, 6, iB3D + 1)
        plt.hist(sc_rate[iB3D], color=colors[1], density=True, bins=30)
        plt.hist(sc_rate_b3d[iB3D], color=colors[6], density=True, bins=30)
        plt.xlabel("shoreline change rate [m/yr]")
        plt.ylabel("density")
    plt.legend(["CASCADE", "B3D only"])
    plt.rcParams["legend.loc"] = "lower right"
    plt.tight_layout()


# ===================================================
# 11: Statistics from the human dynamics modules; combines time series of human modified variables such that
# the 0.5 year time step represents post-storm morphology and pre-human modifications to 1) dune and barrier height
# variables for the Roadway Module, and 2) dune height, barrier height, barrier width, x_s, s_sf, beach width for the
# BeachDuneManager Module
#
#       inputs:         - CASCADE_b3d (a list of B3D objects)
#                       - ib3d(an integer corresponding to the B3D subdomain / brie alongshore grid)
#                       - tmax_roadways (time step where roadway is abandoned)
#                       - tmax_sim (time step where simulation ends)
#                       - nourishments (a list of nourishments objects, BeachDuneManager only)
#       outputs:        - fig


def combine_post_storm_human_time_series(
    tmax_sim, post_storm_statistic, human_modified_statistic
):

    time = np.arange(0, tmax_sim - 0.5, 0.5)
    combined_statistic = [None] * (len(time))
    combined_statistic[::2] = human_modified_statistic
    combined_statistic[1::2] = post_storm_statistic[1:]
    combined_statistic = np.array(combined_statistic)

    return combined_statistic


def plot_nonlinear_stats_RoadwayManager(
    CASCADE_b3d,
    ib3d,
    tmax_roadways,
    tmax_sim,
    post_storm_dunes=None,
    post_storm_interior=None,
    design_height=None,
    rebuild_threshold=None,
    road_elevation=None,
    dunes_rebuilt=None,
    road_relocated=None,
):

    # variables that need to be combined and plotted: dune height (min, max) and barrier height

    # mean dune height (using both dune domain columns) ---------

    # Maximum height of each row in DuneDomain
    DuneDomainCrest = CASCADE_b3d[ib3d].DuneDomain[0:tmax_sim, :, :].max(axis=2)
    DuneRestart = CASCADE_b3d[ib3d].DuneRestart
    DuneDomainCrest[DuneDomainCrest < DuneRestart] = DuneRestart
    DuneCrestMean = (
        np.mean(DuneDomainCrest, axis=1) + CASCADE_b3d[ib3d].BermEl
    ) * 10  # m MHW
    DuneCrestMin = (
        np.min(DuneDomainCrest, axis=1) + CASCADE_b3d[ib3d].BermEl
    ) * 10  # m MHW
    DuneCrestMax = (
        np.max(DuneDomainCrest, axis=1) + CASCADE_b3d[ib3d].BermEl
    ) * 10  # m MHW

    # these are the dune dynamics saved prior to human modifications (essentially a 0.5 year time step); create
    # a combined time series of pre- and post-human modifications to the dune line
    if post_storm_dunes is not None:
        post_DuneCrestMean = [None]
        post_DuneCrestMin = [None]
        post_DuneCrestMax = [None]

        # same calculation as Hd_AverageTS, but here average post-storm dune-height for each time step
        for t in range(1, tmax_roadways):
            if (
                np.size(post_storm_dunes[t], 1) == 1
            ):  # if the dune domain is only one row, don't take the max
                DuneDomainCrest = post_storm_dunes[t]
            else:
                DuneDomainCrest = post_storm_dunes[t].max(
                    axis=1
                )  # Maximum height of each row in DuneDomain
            DuneDomainCrest[
                DuneDomainCrest < CASCADE_b3d[ib3d].DuneRestart
            ] = CASCADE_b3d[ib3d].DuneRestart
            post_DuneCrestMean.append(
                (np.mean(DuneDomainCrest) + CASCADE_b3d[ib3d].BermEl) * 10
            )  # relative to MHW
            post_DuneCrestMax.append(
                (np.max(DuneDomainCrest) + CASCADE_b3d[ib3d].BermEl) * 10
            )  # relative to MHW
            post_DuneCrestMin.append(
                (np.min(DuneDomainCrest) + CASCADE_b3d[ib3d].BermEl) * 10
            )  # relative to MHW

        # now, append NaNs to the rest of the 0.5 year time step arrays
        empty_nans = np.empty(tmax_sim - tmax_roadways)
        empty_nans[:] = np.nan
        # post_DuneCrestMean = np.concatenate((post_DuneCrestMean, empty_nans))
        post_DuneCrestMax = np.concatenate((post_DuneCrestMax, empty_nans))
        post_DuneCrestMin = np.concatenate((post_DuneCrestMin, empty_nans))

        time = np.arange(0, tmax_sim - 0.5, 0.5)
        combined_DuneCrestMin = [None] * (len(time))
        combined_DuneCrestMin[::2] = DuneCrestMin
        combined_DuneCrestMin[1::2] = post_DuneCrestMin[1:]  # post_DuneCrestMean[1:]
        combined_DuneCrestMin = np.array(combined_DuneCrestMin)

        combined_DuneCrestMax = [None] * (len(time))
        combined_DuneCrestMax[::2] = DuneCrestMax
        combined_DuneCrestMax[1::2] = post_DuneCrestMax[1:]
        combined_DuneCrestMax = np.array(combined_DuneCrestMax)

    # barrier width
    BarrierWidth = (
        np.array(CASCADE_b3d[ib3d].x_b_TS[0:tmax_sim])
        - np.array(CASCADE_b3d[ib3d].x_s_TS[0:tmax_sim])
    ) * 10

    # change in barrier width
    bwts = [(x - BarrierWidth[0]) for x in BarrierWidth[0:tmax_sim]]
    rate = [0]
    for k in range(1, len(bwts)):
        rate.append(bwts[k] - bwts[k - 1])
    bw_rate = (
        rate  # note, np.diff doesn't start with zero rate of change, so we do this calc
    )

    # average interior height
    BarrierHeight = []
    for t in range(0, tmax_sim):
        bh_array = np.array(CASCADE_b3d[ib3d].DomainTS[t]) * 10
        BarrierHeight.append(bh_array[bh_array > 0].mean())

    # average interior height post-storm; create a combined time series of pre- and post-human modifications to the
    # barrier interior (changes are going to be small because just affecting the road)
    if post_storm_interior is not None:
        post_BarrierHeight = [None]

        for t in range(1, tmax_roadways):
            bh_array = np.array(post_storm_interior[t]) * 10
            post_BarrierHeight.append(bh_array[bh_array > 0].mean())

        time = np.arange(0, tmax_sim - 0.5, 0.5)

        empty_nans = np.empty(tmax_sim - tmax_roadways)
        empty_nans[:] = np.nan
        post_BarrierHeight = np.concatenate((post_BarrierHeight, empty_nans))

        combined_BarrierHeight = [None] * (len(time))
        combined_BarrierHeight[::2] = BarrierHeight
        combined_BarrierHeight[1::2] = post_BarrierHeight[1:]
        combined_BarrierHeight = np.array(combined_BarrierHeight)

    # change in barrier height
    bhts = [(x - BarrierHeight[0]) for x in BarrierHeight[0:tmax_sim]]
    rate = [0]
    for k in range(1, len(bhts)):
        rate.append(bhts[k] - bhts[k - 1])
    bh_rate = (
        rate  # note, np.diff doesn't start with zero rate of change, so we do this calc
    )
    if post_storm_interior is not None:

        time_mgmt = np.arange(0, tmax_roadways - 0.5, 0.5)
        bhts = [
            (x - combined_BarrierHeight[0])
            for x in combined_BarrierHeight[0 : len(time_mgmt)]
        ]
        rate = [0]
        for k in range(1, len(bhts)):
            rate.append(bhts[k] - bhts[k - 1])
        management_bh_rate = rate

        time_sim = np.arange(0, tmax_sim - 0.5, 0.5)

        empty_nans = np.empty(len(time_sim) - len(time_mgmt))
        empty_nans[:] = np.nan
        management_bh_rate = np.concatenate((management_bh_rate, empty_nans))

        combined_bh_rate = management_bh_rate
        combined_bh_rate[
            ::2
        ] = bh_rate  # populate the nans later with the yearly bh rate
        combined_bh_rate[0 : len(time_mgmt)] = management_bh_rate[
            0 : len(time_mgmt)
        ]  # just replace the bh_rate during management with the mgmt_bh_rate
        combined_bh_rate = np.array(combined_bh_rate)

    # # shoreline change rate
    scts = [
        (x - CASCADE_b3d[ib3d].x_s_TS[0]) * 10
        for x in CASCADE_b3d[ib3d].x_s_TS[0:tmax_sim]
    ]
    rate = [0]
    for k in range(1, len(scts)):
        rate.append(scts[k] - scts[k - 1])
    sc_rate = rate

    # overwash flux
    Qoverwash = CASCADE_b3d[ib3d].QowTS[0:tmax_sim]

    # individual time series
    plt.figure(figsize=(10, 5))
    plt.subplot(3, 3, 1)

    if post_storm_dunes is not None:
        time = np.arange(0, tmax_sim - 0.5, 0.5)
        yearly_time = np.arange(0, tmax_sim, 1)
        mask1 = np.isfinite(combined_DuneCrestMin)
        mask2 = np.isfinite(combined_DuneCrestMax)
        plt.plot(time[mask1], combined_DuneCrestMin[mask1])
        plt.plot(time[mask2], combined_DuneCrestMax[mask2])
        plt.hlines(CASCADE_b3d[ib3d]._Dmaxel * 10, time[0], time[-1], colors="black")
        plt.plot(
            yearly_time[0:tmax_roadways], rebuild_threshold[0:tmax_roadways]
        )  # , color="red")
        plt.plot(
            yearly_time[0:tmax_roadways], design_height[0:tmax_roadways]
        )  # , color="black")
        plt.plot(
            yearly_time[0:tmax_roadways], road_elevation[0:tmax_roadways]
        )  # , color="gray")
        plt.legend(["min", "max", "rebuild", "design", "road", "max-equil"])

    else:
        plt.plot(DuneCrestMin)
        plt.plot(DuneCrestMax)
        plt.legend(["min", "max"])
    plt.ylabel("dune elevation (m MHW)")
    plt.xlabel("Time (yr)")

    plt.subplot(3, 3, 2)
    if post_storm_interior is not None:
        time = np.arange(0, tmax_sim - 0.5, 0.5)
        mask = np.isfinite(combined_BarrierHeight)
        plt.plot(time[mask], combined_BarrierHeight[mask], "m")
    else:
        plt.plot(BarrierHeight, "m")
    plt.ylabel("Ave interior elevation (m MHW)")
    plt.xlabel("Time (yr)")

    plt.subplot(3, 3, 3)
    plt.plot(BarrierWidth, "m")
    plt.ylabel("Barrier width (m)")
    plt.xlabel("Time (yr)")

    # when are dunes rebuilt
    plt.subplot(3, 3, 4)
    if dunes_rebuilt is not None:
        plt.plot(dunes_rebuilt[0:tmax_sim], "m")
        plt.ylabel("Dunes Rebuilt")
        plt.xlabel("Time (yr)")

    # road relocated
    plt.subplot(3, 3, 5)
    if road_relocated is not None:
        plt.plot(road_relocated[0:tmax_sim], "m")
        plt.ylabel("Road Relocated")
        plt.xlabel("Time (yr)")

    # # dune height vs barrier height
    # plt.subplot(3, 3, 4)
    # plt.scatter(
    #     DuneCrestMean,
    #     BarrierHeight,
    #     c=np.arange(0, np.size(DuneCrestMean), 1),
    #     cmap=cm.viridis,
    # )
    # plt.xlabel("Ave dune ele (m MHW)")
    # plt.ylabel("Ave interior ele (m MHW)")

    # # dune height vs barrier width
    # plt.subplot(3, 3, 5)
    # plt.scatter(
    #     DuneCrestMean,
    #     BarrierWidth,
    #     c=np.arange(0, np.size(DuneCrestMean), 1),
    #     cmap=cm.viridis,
    # )
    # plt.xlabel("Ave dune ele (m MHW)")
    # plt.ylabel("Barrier width (m)")

    plt.subplot(3, 3, 6)
    plt.plot(sc_rate, "m")
    plt.ylabel("Shoreline change rate (m/yr)")
    plt.xlabel("Time (yr)")

    plt.subplot(3, 3, 7)
    plt.plot(Qoverwash, "m")
    plt.ylabel("Overwash flux ($m^3/m$)")
    plt.xlabel("Time (yr)")

    # plt.plot(bh_rate)
    # plt.ylabel("Interior elevation change rate (m MHW/yr)")
    # plt.xlabel("Time (yr)")

    # plt.plot(bw_rate)
    # plt.ylabel("Barrier width change rate (m/yr)")
    # plt.xlabel("Time (yr)")

    plt.subplot(3, 3, 8)
    plt.plot(np.array(CASCADE_b3d[ib3d].x_s_TS[0:tmax_sim]) * 10, "m")
    plt.ylabel("Shoreline position (m)")
    plt.xlabel("Time (yr)")

    plt.subplot(3, 3, 9)
    time = np.arange(0, tmax_sim - 0.5, 0.5)
    plt.plot(CASCADE_b3d[ib3d].s_sf_TS[0:tmax_sim], "m")
    plt.hlines(CASCADE_b3d[ib3d]._s_sf_eq, time[0], time[-1], colors="black")
    plt.ylabel("shoreface slope")
    plt.xlabel("Time (yr)")

    # dune height vs shoreline position
    # plt.plot(DuneCrestMean, np.array(CASCADE_b3d[ib3d]._x_s_TS) * 10)
    # plt.plot(
    #     DuneCrestMean[0:tmax],
    #     np.array(CASCADE_b3d[ib3d]._x_s_TS[0:tmax]) * 10,
    #     "m",
    # )
    # plt.scatter(
    #     DuneCrestMean,
    #     np.array(CASCADE_b3d[ib3d]._x_s_TS) * 10,
    #     c=np.arange(0, np.size(DuneCrestMean), 1),
    #     cmap="Greens",
    # )
    # plt.xlabel("Ave Dune Height (m)")
    # plt.ylabel("Shoreline position (m)")

    # dune height vs shoreline change rate
    # plt.scatter(DuneCrestMean, sc_rate)
    # plt.plot(DuneCrestMean[0:tmax], sc_rate[0:tmax], "m")
    # plt.scatter(
    #     DuneCrestMean,
    #     sc_rate,
    #     c=np.arange(0, np.size(DuneCrestMean), 1),
    #     cmap=cm.viridis,  # "Greens"
    # )
    # plt.xlabel("Ave Dune Height (m)")
    # plt.ylabel("Shoreline change rate (m/yr)")

    # dune height vs overwash flux
    # plt.scatter(DuneCrestMean, Qoverwash)
    # plt.plot(DuneCrestMean[0:tmax], Qoverwash[0:tmax], "m")
    # plt.scatter(
    #     DuneCrestMean,
    #     Qoverwash,
    #     c=np.arange(0, np.size(DuneCrestMean), 1),
    #     cmap=cm.viridis,
    # )
    # plt.xlabel("Ave Dune Height (m)")
    # plt.ylabel("Qoverwash ($m^3/m$)")

    # # dune height vs barrier width change rate
    # plt.subplot(2, 5, 6)
    # # plt.subplot(5, 2, 2)
    # plt.scatter(
    #     DuneCrestMean,
    #     np.array(bw_rate),
    #     c=np.arange(0, np.size(DuneCrestMean), 1),
    #     cmap=cm.viridis,
    # )
    # plt.xlabel("Ave Dune Height (m)")
    # plt.ylabel("Change in barrier width (m/yr)")

    # # now try my hand at a 3d plot
    # from mpl_toolkits import mplot3d
    #
    # plt.figure()
    # ax = plt.axes(projection="3d")
    #
    # # three-dimensional line
    # # ax.plot3D(DuneCrestMean, BarrierWidth, sc_rate, "gray")
    # ax.plot3D(
    #     DuneCrestMean[0:tmax], BarrierWidth[0:tmax], sc_rate[0:tmax], "gray"
    # )
    #
    # # three-dimensional scattered points
    # ax.scatter3D(
    #     DuneCrestMean[0:tmax],
    #     BarrierWidth[0:tmax],
    #     sc_rate[0:tmax],
    #     c=sc_rate[0:tmax],
    #     cmap="Greens",
    # )
    # ax.set_title(str(0) + "-" + str(tmax) + " yrs")
    # ax.set_xlabel("Dune height (m)")
    # ax.set_ylabel("Barrier width (m)")
    # ax.set_zlabel("Shoreline change rate (m/yr)")
    #
    # # try Qoverwash as the x axis instead of dune height
    # fig = plt.figure()
    # ax = plt.axes(projection="3d")
    #
    # # three-dimensional line
    # # ax.plot3D(DuneCrestMean, BarrierWidth, sc_rate, "gray")
    # ax.plot3D(Qoverwash[0:tmax], BarrierWidth[0:tmax], sc_rate[0:tmax], "gray")
    #
    # ax.scatter3D(
    #     Qoverwash[0:tmax],
    #     BarrierWidth[0:tmax],
    #     sc_rate[0:tmax],
    #     c=sc_rate[0:tmax],
    #     cmap="Greens",
    # )
    # ax.set_title(str(0) + "-" + str(tmax) + " yrs")
    # ax.set_xlabel("Overwash flux ($m^3/m$)")
    # ax.set_ylabel("Barrier width (m)")
    # ax.set_zlabel("Shoreline change rate (m/yr)")
    #
    # # try Qoverwash, dune height, and barrier width change
    # fig = plt.figure()
    # ax = plt.axes(projection="3d")
    #
    # ax.plot3D(
    #     DuneCrestMean[0:tmax], Qoverwash[0:tmax], bw_rate[0:tmax], "gray"
    # )
    #
    # # three-dimensional scattered points
    # ax.scatter3D(
    #     DuneCrestMean[0:tmax],
    #     Qoverwash[0:tmax],
    #     bw_rate[0:tmax],
    #     c=bw_rate[0:tmax],
    #     cmap="Greens",
    # )
    # ax.set_title(str(0) + "-" + str(tmax) + " yrs")
    # ax.set_xlabel("Dune height (m)")
    # ax.set_ylabel("Overwash flux ($m^3/m$)")
    # ax.set_zlabel("Barrier width change rate (m/yr)")

    plt.tight_layout()

    if post_storm_dunes is not None:
        DuneCrestMin = combined_DuneCrestMin
        DuneCrestMax = combined_DuneCrestMax

    if post_storm_interior is not None:
        BarrierHeight = combined_BarrierHeight
        bh_rate = combined_bh_rate

    return (
        BarrierWidth,
        DuneCrestMean,
        BarrierHeight,
        bw_rate,
        bh_rate,
        sc_rate,
        DuneCrestMin,
        DuneCrestMax,
    )


def plot_nonlinear_stats_BeachDuneManager(
    CASCADE_b3d,
    ib3d,
    tmax_management,
    tmax_sim,
    nourishments,
    post_storm_dunes=None,
    post_storm_x_s=None,
    post_storm_s_sf=None,
    post_storm_ave_interior_width=None,
    post_storm_ave_interior_height=None,
    post_storm_beach_width=None,
    post_storm_Qow=None,
    design_elevation=None,
    rebuild_threshold=None,
    dunes_rebuilt=None,
):
    # if the post-storm variables are not supplied (essentially a 0.5 year time step), then only the human-modified
    # statistics are plotted (the 1 year time step)

    # variables that need to be combined and plotted: dune height, barrier width, barrier height, x_s, s_sf,
    # beach width, Qow
    empty_nans = np.empty(tmax_sim - tmax_management)
    empty_nans[:] = np.nan

    # dune height (using both dune domain columns) --------------------------------------------------

    # Maximum height of each row in DuneDomain
    DuneDomainCrest = CASCADE_b3d[ib3d].DuneDomain[0:tmax_sim, :, :].max(axis=2)
    DuneRestart = CASCADE_b3d[ib3d].DuneRestart
    DuneDomainCrest[DuneDomainCrest < DuneRestart] = DuneRestart
    DuneCrestMean = (
        np.mean(DuneDomainCrest, axis=1) + CASCADE_b3d[ib3d].BermEl
    ) * 10  # m MHW
    DuneCrestMin = (
        np.min(DuneDomainCrest, axis=1) + CASCADE_b3d[ib3d].BermEl
    ) * 10  # m MHW
    DuneCrestMax = (
        np.max(DuneDomainCrest, axis=1) + CASCADE_b3d[ib3d].BermEl
    ) * 10  # m MHW

    # these are the dune dynamics saved prior to human modifications (essentially a 0.5 year time step)
    if post_storm_dunes is not None:

        post_storm_DuneCrestMin = [None]
        post_storm_DuneCrestMax = [None]

        # same calculation as Hd_AverageTS
        for t in range(1, tmax_management):
            if (
                np.size(post_storm_dunes[t], 1) == 1
            ):  # if the dune domain is only one row, don't take the max
                DuneDomainCrest = post_storm_dunes[t]
            else:
                DuneDomainCrest = post_storm_dunes[t].max(
                    axis=1
                )  # Maximum height of each row in DuneDomain
            DuneDomainCrest[
                DuneDomainCrest < CASCADE_b3d[ib3d].DuneRestart
            ] = CASCADE_b3d[ib3d].DuneRestart
            post_storm_DuneCrestMax.append(
                (np.max(DuneDomainCrest) + CASCADE_b3d[ib3d].BermEl) * 10
            )  # relative to MHW
            post_storm_DuneCrestMin.append(
                (np.min(DuneDomainCrest) + CASCADE_b3d[ib3d].BermEl) * 10
            )  # relative to MHW

        # now, append NaNs to the rest of the 0.5 year time step arrays
        post_storm_DuneCrestMax = np.concatenate((post_storm_DuneCrestMax, empty_nans))
        post_storm_DuneCrestMin = np.concatenate((post_storm_DuneCrestMin, empty_nans))

        combined_DuneCrestMin = combine_post_storm_human_time_series(
            tmax_sim, post_storm_DuneCrestMin, DuneCrestMin
        )
        combined_DuneCrestMax = combine_post_storm_human_time_series(
            tmax_sim, post_storm_DuneCrestMax, DuneCrestMax
        )

    # barrier width --------------------------------------------------
    # note that here, and everywhere in the drowning paper, barrier width refers to the average interior width, and
    # not x_b - x_s, which incorporates changes in beach width and the dune line
    BarrierWidth = (
        np.array(CASCADE_b3d[ib3d].InteriorWidth_AvgTS[0:tmax_sim])
    ) * 10  # m

    # barrier width prior to human modifications (essentially a 0.5 year time step)
    if post_storm_ave_interior_width is not None:

        post_storm_ave_interior_width[0] = np.nan
        post_storm_BarrierWidth = (
            np.array(post_storm_ave_interior_width[0:tmax_management])
        ) * 10

        post_storm_BarrierWidth = np.concatenate((post_storm_BarrierWidth, empty_nans))
        combined_BarrierWidth = combine_post_storm_human_time_series(
            tmax_sim, post_storm_BarrierWidth, BarrierWidth
        )  # m

    # change in barrier width --------------------------------------------------
    bwts = [(x - BarrierWidth[0]) for x in BarrierWidth[0:tmax_sim]]
    rate = [0]
    for k in range(1, len(bwts)):
        rate.append(bwts[k] - bwts[k - 1])
    bw_rate = rate  # note, np.diff doesn't start with zero rate of change, so we do this calc  # m/yr

    # barrier height --------------------------------------------------
    BarrierHeight = (np.array(CASCADE_b3d[ib3d].h_b_TS[0:tmax_sim])) * 10  # m

    # before human modifications
    if post_storm_ave_interior_height is not None:

        post_storm_ave_interior_height[0] = np.nan
        post_storm_BarrierHeight = (
            np.array(post_storm_ave_interior_height[0:tmax_management])
        ) * 10  # m

        post_storm_BarrierHeight = np.concatenate(
            (post_storm_BarrierHeight, empty_nans)
        )
        combined_BarrierHeight = combine_post_storm_human_time_series(
            tmax_sim, post_storm_BarrierHeight, BarrierHeight
        )  # m

    # change in barrier height --------------------------------------------------
    bhts = [(x - BarrierHeight[0]) for x in BarrierHeight[0:tmax_sim]]
    rate = [0]
    for k in range(1, len(bhts)):
        rate.append(bhts[k] - bhts[k - 1])
    bh_rate = rate  # note, np.diff doesn't start with zero rate of change, so we do this calc  # m/yr

    # shoreline position and change rate --------------------------------------------------
    shoreline_position = np.array(CASCADE_b3d[ib3d].x_s_TS[0:tmax_sim]) * 10  # m
    scts = [(x - shoreline_position[0]) for x in shoreline_position]
    rate = [0]
    for k in range(1, len(scts)):
        rate.append(scts[k] - scts[k - 1])
    sc_rate = rate  # m/yr

    # before human modifications
    if post_storm_x_s is not None:

        post_storm_shoreline_position = np.array(post_storm_x_s[0:tmax_management])
        post_storm_shoreline_position[0] = np.nan
        post_storm_shoreline_position = post_storm_shoreline_position * 10  # m
        post_storm_shoreline_position = np.concatenate(
            (post_storm_shoreline_position, empty_nans)
        )
        combined_shoreline_position = combine_post_storm_human_time_series(
            tmax_sim, post_storm_shoreline_position, shoreline_position
        )  # m

    # overwash flux --------------------------------------------------
    net_overwash = np.array(CASCADE_b3d[ib3d].QowTS[0:tmax_sim])  # m^3/m

    # before human modifications
    if post_storm_Qow is not None:

        post_storm_Qow[0] = np.nan
        post_storm_overwash = np.array(
            post_storm_Qow[0:tmax_management]
        )  # this is the total overwash before removal
        post_storm_overwash = np.concatenate((post_storm_overwash, empty_nans))
        combined_overwash = combine_post_storm_human_time_series(
            tmax_sim, post_storm_overwash, net_overwash
        )

    # shoreface slope -----------------------------------------
    shoreface_slope = CASCADE_b3d[ib3d].s_sf_TS[0:tmax_sim]
    equilibrium_slope = CASCADE_b3d[ib3d]._s_sf_eq

    # before human modifications
    if post_storm_s_sf is not None:

        post_storm_s_sf[0] = np.nan
        post_storm_shoreface_slope = post_storm_s_sf[0:tmax_management]
        post_storm_shoreface_slope = np.concatenate(
            (post_storm_shoreface_slope, empty_nans)
        )
        combined_shoreface_slope = combine_post_storm_human_time_series(
            tmax_sim, post_storm_shoreface_slope, shoreface_slope
        )

    # beach width -----------------------------------------
    beach_width = nourishments[ib3d].beach_width[0:tmax_sim]  # m

    # before human modifications
    if post_storm_beach_width is not None:

        post_storm_beach_width[0] = np.nan
        post_storm_beach_width = post_storm_beach_width[0:tmax_management]
        post_storm_beach_width = np.concatenate((post_storm_beach_width, empty_nans))
        combined_beach_width = combine_post_storm_human_time_series(
            tmax_sim, post_storm_beach_width, beach_width
        )

    # individual time series
    plt.figure(figsize=(10, 8))
    # plt.subplot(4, 3, 1)

    full_time = np.arange(0, tmax_sim - 0.5, 0.5)

    plt.subplot(2, 2, 1)
    if post_storm_dunes is not None:
        full_time = np.arange(0, tmax_sim - 0.5, 0.5)
        time_mgmt = np.arange(0, tmax_management, 1)
        mask1 = np.isfinite(combined_DuneCrestMin)
        mask2 = np.isfinite(combined_DuneCrestMax)
        plt.plot(full_time[mask1], combined_DuneCrestMin[mask1])
        plt.plot(full_time[mask2], combined_DuneCrestMax[mask2])
        if rebuild_threshold is not None:
            plt.hlines(rebuild_threshold, time_mgmt[0], time_mgmt[-1], colors="green")
        plt.hlines(design_elevation, time_mgmt[0], time_mgmt[-1], colors="red")
        plt.hlines(
            CASCADE_b3d[ib3d]._Dmaxel * 10, full_time[0], full_time[-1], colors="black"
        )
        if rebuild_threshold is not None:
            plt.legend(["min", "max", "rebuild", "design", "max-equil"])
        else:
            plt.legend(["min", "max", "design", "max-equil"])
    else:
        plt.plot(DuneCrestMin)
        plt.plot(DuneCrestMax)
        plt.legend(["min", "max"])
    plt.ylabel("dune elevation (m MHW)")
    plt.xlabel("Time (yr)")

    # when are dunes rebuilt
    # plt.subplot(4, 3, 4)
    plt.subplot(2, 2, 2)
    if dunes_rebuilt is not None:
        plt.plot(dunes_rebuilt, "k")
        plt.ylabel("Dunes Rebuilt")
        plt.xlabel("Time (yr)")

    # nourishment volumes
    # plt.subplot(4, 3, 5)
    plt.subplot(2, 2, 3)
    barrier_length = CASCADE_b3d[ib3d].BarrierLength * 10
    plt.plot(nourishments[ib3d].nourishment_volume_TS, "k")
    plt.plot(
        nourishments[ib3d].rebuild_dune_volume_TS / barrier_length, "--"
    )  # m^3 to m^3/m
    plt.ylabel("Nourishment volume (m$^3$/m)")
    plt.legend(["shoreface", "dunes"])
    plt.xlabel("Time (yr)")

    # shoreface slope
    plt.subplot(2, 2, 4)
    # plt.subplot(4, 3, 10)
    if post_storm_s_sf is not None:
        mask = np.isfinite(combined_shoreface_slope)
        plt.plot(full_time[mask], combined_shoreface_slope[mask], "m")
        plt.plot(shoreface_slope, "k")
        plt.legend(["includes post-storm", "mgmt only"])
    else:
        plt.plot(shoreface_slope, "k")
    plt.hlines(equilibrium_slope, full_time[0], full_time[-1], colors="green")
    plt.ylabel("shoreface slope")
    plt.xlabel("Time (yr)")

    plt.tight_layout()

    plt.figure(figsize=(10, 8))

    # shoreline change rate
    plt.subplot(2, 2, 3)
    # plt.subplot(4, 3, 6)
    plt.plot(sc_rate, "k")
    plt.ylabel("Shoreline change rate (m/yr)")
    plt.xlabel("Time (yr)")

    # beach width
    plt.subplot(2, 2, 1)
    # plt.subplot(4, 3, 7)
    if post_storm_beach_width is not None:
        plt.plot(full_time, combined_beach_width, "m")
        plt.plot(beach_width, "k")
        # plt.legend(["includes post-storm", "mgmt only"])
    else:
        plt.plot(beach_width, "k")
    plt.ylabel("Beach width (m)")
    plt.xlabel("Time (yr)")
    plt.xlim([0, tmax_sim])

    # shoreline position
    plt.subplot(2, 2, 2)
    # plt.subplot(4, 3, 9)
    if post_storm_x_s is not None:
        mask = np.isfinite(combined_shoreline_position)
        plt.plot(full_time[mask], combined_shoreline_position[mask], "m")
        combined_dune_toe = np.array(combined_shoreline_position) + np.array(
            combined_beach_width
        )
        mask = np.isfinite(combined_dune_toe)
        plt.plot(full_time[mask], combined_dune_toe[mask], "--m")
        plt.legend(["shoreline", "dune toe"])
        # plt.legend(["includes post-storm", "mgmt only"])
    dune_toe = np.array(shoreline_position) + np.array(beach_width)
    plt.plot(shoreline_position, "k")
    plt.plot(dune_toe, "--k")
    plt.ylabel("Cross-shore position (m)")
    plt.xlabel("Time (yr)")
    # plt.legend(["shoreline", "dune toe"])

    plt.tight_layout()

    plt.figure(figsize=(10, 8))

    # barrier height
    # plt.subplot(4, 3, 2)
    plt.subplot(2, 2, 1)
    if post_storm_BarrierHeight is not None:
        mask = np.isfinite(combined_BarrierHeight)
        plt.plot(full_time[mask], combined_BarrierHeight[mask], "m")
        # plt.legend(["includes post-storm", "mgmt only"])
    plt.plot(BarrierHeight, "k")
    plt.ylabel("barrier height (m MHW)")
    plt.xlabel("Time (yr)")

    # barrier width
    plt.subplot(2, 2, 2)
    # plt.subplot(4, 3, 3)
    if post_storm_BarrierWidth is not None:
        mask = np.isfinite(combined_BarrierWidth)
        plt.plot(full_time[mask], combined_BarrierWidth[mask], "m")
        # plt.legend(["includes post-storm", "mgmt only"])
    plt.plot(BarrierWidth, "k")
    plt.ylabel("barrier width (m)")
    plt.xlabel("Time (yr)")

    # overwash flux
    # plt.subplot(4, 3, 8)
    plt.subplot(2, 2, 3)
    if post_storm_Qow is not None:
        plt.plot(full_time, combined_overwash, "m")
        # plt.legend(["includes post-storm", "mgmt only"])
    plt.plot(net_overwash, "k")
    plt.ylabel("Overwash flux ($m^3/m$)")
    plt.xlabel("Time (yr)")

    plt.tight_layout()

    if post_storm_dunes is not None:
        DuneCrestMin = combined_DuneCrestMin
        DuneCrestMax = combined_DuneCrestMax

    if post_storm_BarrierHeight is not None:
        BarrierHeight = combined_BarrierHeight

    if post_storm_BarrierWidth is not None:
        BarrierWidth = combined_BarrierWidth

    # if post_storm_x_s is not None:
    # shoreline_position = combined_shoreline_position

    if post_storm_beach_width is not None:
        beach_width = combined_beach_width
        # dune_toe = combined_dune_toe

    if post_storm_s_sf is not None:
        shoreface_slope = combined_shoreface_slope

    if post_storm_Qow is not None:
        overwash = combined_overwash

    return (
        BarrierWidth,
        DuneCrestMean,
        BarrierHeight,
        bh_rate,
        bw_rate,
        sc_rate,
        DuneCrestMin,
        DuneCrestMax,
        shoreline_position,  # not combined
        shoreface_slope,
        beach_width,
        overwash,
        dune_toe,  # not combined
    )


# ===================================================
# The remaining figures are all used to create paper figures for "Pathways to barrier drowning from coastal management"


def plot_nonlinear_stats_mgmt_array4(
    cascade,  # these are lists
    DuneCrestMin,
    DuneCrestMax,
    BarrierHeight,
    BarrierWidth,
    DuneCrestMean,
    shoreline_position,
    overwash,
    dune_toe,
    TMAX,
    tmax_management,
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
            )
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
        if roadways_on:
            axs[i].set_xlim([-15, 765])
            axs[i].set_ylim([0, 5.2])
        if nourishment_on:
            # axs[i].set_xlim([-15, 515])
            axs[i].set_xlim([-15, 765])
            axs[i].set_ylim([0, 5.2])
    axs[0].set(ylabel="elevation (m MHW)")
    if roadways_on:
        axs[3].legend(
            [
                "dune along. min",
                "dune along. max",
                "dune rebuild",
                "dune design",
                "road",
                "dune max-equil",
            ]
        )
    if nourishment_on:
        axs[3].legend(
            [
                "dune along. min",
                "dune along. max",
                "dune rebuild",
                "dune design",
                "dune max-equil",
            ]
        )
    plt.tight_layout()

    # interior height and width -------------------------------------------------------------------------------------- #
    fig2, axs2 = plt.subplots(1, 4, figsize=(10, 3), sharex=True)
    color = ["b", "r", "g", "m"]
    for i in range(len(cascade)):

        if i > 0:
            time = np.arange(0, TMAX[i] - 0.5, 0.5)
            yearly_time = np.arange(0, TMAX[i], 1)

            barrier_height = BarrierHeight[i][0 : len(time)]
            mask_bh = np.isfinite(barrier_height)

            scts = [
                (x - shoreline_position[i][0])
                for x in shoreline_position[i][0 : TMAX[i]]
            ]
            axs2[2].plot(yearly_time, scts, color[i])

            if roadways_on:
                axs2[0].plot(time[mask_bh], barrier_height[mask_bh], color[i])
                axs2[1].plot(yearly_time, BarrierWidth[i][0 : TMAX[i]], color[i])
                axs2[3].plot(yearly_time, overwash[0 : TMAX[i]], color[i])

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

                barrier_width = BarrierWidth[i][0 : len(time)]
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
        if roadways_on:
            axs2[i].set_xlim([-15, 765])
        if nourishment_on:
            # axs2[i].set_xlim([-15, 515])
            axs2[i].set_xlim([-15, 765])

    axs2[0].set(ylabel="barrier elevation (m MHW)")
    axs2[0].set_ylim([-0.03, 1.75])
    axs2[1].set(ylabel="barrier width (m)")
    axs2[1].set_ylim([-6, 425])
    if roadways_on:
        axs2[2].set_ylim([-10, 615])
        axs2[2].set(ylabel="shoreline position (m)")
    if nourishment_on:
        axs2[2].set_ylim([-30, 615])
        axs2[2].set(ylabel="cross-shore position (m)")
        axs2[2].legend(["shoreline", "dune toe"])
    axs2[3].set(ylabel="overwash flux (m$^3$/m)")
    axs2[3].set_ylim([-3, 225])
    if roadways_on:
        scenarios = ["natural", "1-m dune", "2-m dune", "3-m dune"]
        axs2[2].legend(scenarios)
    if nourishment_on:
        axs2[2].legend(scenarios)
    plt.tight_layout()

    # plot a zoom in of nourishments
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


def plot_nonlinear_stats_low_high_sensitivity(
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


# def summary_table_stats_low_high(cascade):
#     fig = go.Figure(
#         data=[
#             go.Table(
#                 header=dict(values=["A Scores", "B Scores"]),
#                 cells=dict(values=[[100, 90, 80, 90], [95, 85, 75, 95]]),
#             )
#         ]
#     )
#     fig.show()


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


def fig2_initialCNH_topo(
    cascade_model_list,  # must be from a nourishments simulation (i.e., have a beach width)
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
            cascade.nourishments[iB3D]._post_storm_beach_width[0] / 10
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
        AnimateDomain[
            OriginTstart:OriginTstop, xOrigin : xOrigin + BarrierLength
        ] = Domain

        # plot
        print(np.max(AnimateDomain))
        cax = axs[iCascade].matshow(
            AnimateDomain, origin="lower", cmap="terrain", vmin=-1.1, vmax=3.0
        )  # , interpolation='gaussian') # analysis:ignore
        axs[iCascade].xaxis.set_ticks_position("bottom")
        axs[iCascade].set(xlabel="alongshore distance (dam)")
        axs[iCascade].set_ylim([40, 110])

    axs[0].set(ylabel="cross-shore distance (dam)")
    cbar = fig.colorbar(cax)
    cbar.set_label("elevation (m MHW)", rotation=270)
    # plt.tight_layout()

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
            axs[0].set_xlim([0, 110])
            axs[0].set(xlabel="cross-shore distance (dam)")
            axs[0].set(ylabel="elevation (m MHW)")
        else:  # 0.75 cases
            axs[1].plot(x, y)
            axs[1].hlines(
                sea_level * 10, shoreface_toe_x, end_of_bay_x, colors="dodgerblue"
            )
            axs[1].set_xlim([0, 110])
            axs[1].set(xlabel="cross-shore distance (dam)")
    # plt.tight_layout()


def fig3_slr_sensitivity(
    cascade,  # lists
    TMAX,
):
    # col 1 - barrier height (4 different SLR)
    # col 2 - barrier width (4 different SLR)
    # col 3 - shoreline position (4 different SLR)

    fig1, axs1 = plt.subplots(1, 4, figsize=(10, 3), sharex=True)
    for i in range(len(cascade)):
        time = np.arange(0, TMAX[i], 1)

        RSLR = cascade[i].barrier3d[0].RSLR[0 : TMAX[i]]  # in dam
        SLR = np.cumsum(np.array(RSLR) * 10)

        axs1[0].plot(time, SLR)

    axs1[0].set(ylabel="sea level (m)")
    axs1[0].set(xlabel="time (yr)")
    axs1[0].legend(["0.004 m/yr", "0.008", "0.012", "Acc"])
    plt.tight_layout()

    fig2, axs2 = plt.subplots(1, 4, figsize=(10, 3), sharex=True)
    for i in range(len(cascade)):
        time = np.arange(0, TMAX[i], 1)

        BarrierWidth = (
            np.array(cascade[i].barrier3d[0].x_b_TS[0 : TMAX[i]])
            - np.array(cascade[i].barrier3d[0].x_s_TS[0 : TMAX[i]])
        ) * 10

        # average interior height
        BarrierHeight = []
        for t in range(0, TMAX[i]):
            bh_array = np.array(cascade[i].barrier3d[0]._DomainTS[t]) * 10  # dam to m
            BarrierHeight.append(bh_array[bh_array > 0].mean())

        scts = [
            (x - cascade[i].barrier3d[0].x_s_TS[0]) * 10
            for x in cascade[i].barrier3d[0].x_s_TS[0 : TMAX[i]]
        ]

        axs2[0].plot(time, BarrierHeight)
        axs2[1].plot(time, BarrierWidth)
        axs2[2].plot(time, scts)

        axs2[i].set(xlabel="time (yr)")

    axs2[0].set(ylabel="barrier height (m MHW)")
    axs2[0].set_ylim([0.15, 1.75])
    axs2[1].set(ylabel="barrier width (m)")
    axs2[1].set_ylim([50, 425])
    axs2[2].set(ylabel="shoreline position (m)")
    axs2[2].set_ylim([-10, 480])
    axs2[2].legend(["0.004 m/yr", "0.008", "0.012", "Acc"])
    plt.tight_layout()


def supp_10kyr_timeseries():

    datadir = "/Users/KatherineAnardeWheels/Research/BARis/UNC/CNH/CASCADE_save_dir/Run_Output/"

    b3d = []
    bw = []
    bh = []
    # xs = []
    # storm_count = []
    dune_height = []
    time = []

    ib3d = 0
    tmax = [5725, 992, 4870, 10000, 6669]

    for i in range(1, 6):
        name_prefix = "4-CASCADE_noAST_Rave_pt75_SLR_pt004_10k-yrs_0" + str(i)
        output = np.load(datadir + name_prefix + ".npz", allow_pickle=True)
        cascade = output["cascade"][0]
        b3d.append(cascade.barrier3d)

        BarrierWidth = (
            np.array(cascade.barrier3d[ib3d].x_b_TS[0 : tmax[i - 1]])
            - np.array(cascade.barrier3d[ib3d].x_s_TS[0 : tmax[i - 1]])
        ) * 10
        bw.append(BarrierWidth)

        BarrierHeight = []
        for t in range(0, tmax[i - 1]):
            bh_array = np.array(cascade.barrier3d[ib3d].DomainTS[t]) * 10  # dam to m
            BarrierHeight.append(bh_array[bh_array > 0].mean())
        bh.append(BarrierHeight)

        # xs.append(np.array(cascade.barrier3d[ib3d].x_s_TS[0 : tmax[i - 1]]) * 10)
        # storm_count.append(cascade.barrier3d[ib3d]._StormCount[0 : tmax[i - 1]])

        time.append(np.arange(0, tmax[i - 1], 1))

        dune_height.append(
            np.array(cascade.barrier3d[ib3d]._Hd_AverageTS[0 : tmax[i - 1]]) * 10
        )

    fig1, axs1 = plt.subplots(3, 1, figsize=(5, 7), sharex=True)

    tlow = 0
    thigh = 10000

    axs1[0].plot(
        time[0], bw[0], time[1], bw[1], time[2], bw[2], time[3], bw[3], time[4], bw[4]
    )
    axs1[0].set(ylabel="barrier width (m)")
    axs1[0].set_xlim([tlow, thigh])
    axs1[0].set_ylim([-10, 480])
    axs1[0].axvline(x=829)
    axs1[0].axvline(x=3284)

    axs1[1].plot(
        time[0], bh[0], time[1], bh[1], time[2], bh[2], time[3], bh[3], time[4], bh[4]
    )
    axs1[1].set(ylabel="barrier height (m MHW)")
    axs1[1].set_xlim([tlow, thigh])
    axs1[1].set_ylim([-0.10, 1.8])
    axs1[1].axvline(x=829)
    axs1[1].axvline(x=3284)
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

    # 0.45 case

    b3d = []
    bw = []
    bh = []
    # xs = []
    # storm_count = []
    dune_height = []
    time = []

    ib3d = 0
    tmax = [10000, 10000, 10000, 10000, 10000]

    for i in range(1, 6):
        name_prefix = "4-CASCADE_noAST_Rave_pt45_SLR_pt004_10k-yrs_0" + str(i)
        output = np.load(datadir + name_prefix + ".npz", allow_pickle=True)
        cascade = output["cascade"][0]
        b3d.append(cascade.barrier3d)

        BarrierWidth = (
            np.array(cascade.barrier3d[ib3d].x_b_TS[0 : tmax[i - 1]])
            - np.array(cascade.barrier3d[ib3d].x_s_TS[0 : tmax[i - 1]])
        ) * 10
        bw.append(BarrierWidth)

        BarrierHeight = []
        for t in range(0, tmax[i - 1]):
            bh_array = np.array(cascade.barrier3d[ib3d].DomainTS[t]) * 10  # dam to m
            BarrierHeight.append(bh_array[bh_array > 0].mean())
        bh.append(BarrierHeight)

        # xs.append(np.array(cascade.barrier3d[ib3d].x_s_TS[0 : tmax[i - 1]]) * 10)
        # storm_count.append(cascade.barrier3d[ib3d]._StormCount[0 : tmax[i - 1]])

        time.append(np.arange(0, tmax[i - 1], 1))

        dune_height.append(
            np.array(cascade.barrier3d[ib3d]._Hd_AverageTS[0 : tmax[i - 1]]) * 10
        )

    fig1, axs1 = plt.subplots(3, 1, figsize=(5, 7), sharex=True)

    tlow = 0
    thigh = 10000

    axs1[0].plot(
        time[0], bw[0], time[1], bw[1], time[2], bw[2], time[3], bw[3], time[4], bw[4]
    )
    axs1[0].set(ylabel="barrier width (m)")
    axs1[0].set_xlim([tlow, thigh])
    axs1[0].set_ylim([-10, 480])
    axs1[0].axvline(x=802)
    axs1[0].axvline(x=8750)

    axs1[1].plot(
        time[0], bh[0], time[1], bh[1], time[2], bh[2], time[3], bh[3], time[4], bh[4]
    )
    axs1[1].set(ylabel="barrier height (m MHW)")
    axs1[1].set_xlim([tlow, thigh])
    axs1[1].set_ylim([-0.10, 1.8])
    axs1[1].axvline(x=802)
    axs1[1].axvline(x=8750)
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


## OLD CODE THAT I NEED TO FIX EVENTUALLY:

# plot_BRIE_frames
# CASCADEplt.plot_BRIE_frames(brieLTA._y,
#                             brieLTA._x_t_save,
#                             brieLTA._x_s_save,
#                             brieLTA._x_b_save,
#                             brieLTA._nt,
#                             brieLTA._ny,
#                             brieLTA._dy,
#                             brieLTA._dtsave,
#                             brieLTA._inlet_idx,
#                             'True',
#                             'test')

# y = brieLTA._y
# x_t = brieLTA._x_t_save
# x_s = brieLTA._x_s_save
# x_b = brieLTA._x_b_save
# nt = brieLTA._nt
# ny = brieLTA._ny
# dy = brieLTA._dy
# dtsave = brieLTA._dtsave
# inlet_idx = brieLTA._inlet_idx
# make_gif = 'True'
# file_name = 'test'
# directory = "/Users/KatherineAnardeWheels/PycharmProjects/CASCADE/"
#
# # create time array
# t = np.arange(0, nt, dtsave)
#
# if make_gif:
#     os.chdir(directory)
#     os.mkdir(directory+'/GIF')
#
# fig, axes = plt.subplots()
# frames = []
# ymax = int(math.ceil(np.max(x_b[:, -1]) / 100.0)) * 100
# ymin = int(math.floor(np.min(x_t[:, 0]) / 100.0)) * 100
#
# for i in np.arange(0, np.size(t)):
#
#     # plot initial conditions
#     axes.plot(y / 1000, x_t[:, i], color="cornflowerblue")
#     axes.plot(y / 1000, x_s[:, i], color="gold")
#     axes.plot(y / 1000, x_b[:, i], color="teal")
#     axes.fill_between(y / 1000, ymin, x_t[:, i], color="royalblue")
#     axes.fill_between(y / 1000, x_b[:, i], ymax, color="teal", alpha=0.6)
#     axes.fill_between(y / 1000, x_t[:, i], x_s[:, i], color="cornflowerblue", alpha=0.6)
#     axes.fill_between(y / 1000, x_s[:, i], x_b[:, i], color="gold", alpha=0.6)
#     axes.legend(['x_t', 'x_s', 'x_b'])
#     plt.xlabel('Alongshore (km)')
#     plt.ylabel('Cross-shore (m)')
#     plt.title('Time = ' + str(int(t[i])) + ' yr')
#     axes.set_xlim(0, (ny-1) * dy / 1000)
#     # axes.set_ylim(-2000, 6000)  # KA, will need to update later - placeholder
#     # ratio = 0.4   # aspect ratio
#     # axes.set_aspect(0.05)
#     # axes.margins(0.5)
#     # axes.set_aspect(1/axes.get_data_ratio())
#     axes.set_ylim(ymin, ymax)
#
#     # Here I chose to only mark the inlets (all indices), and not where the barrier volume is less than 0...
#     # need to go back and fix
#     # if np.size(inlet_idx) != 0 :
#     #    axes.plot(y[np.hstack(inlet_idx)]/1000, x_s[np.hstack(inlet_idx)], 'kD')
#
#     if make_gif:
#         filename = 'GIF/year_' + str(int(t[i])) + '.png'
#         fig.savefig(filename, dpi=200)
#
#         # use imageio to create gif
#         frames.append(imageio.imread(filename))
#         os.remove(filename)
#
#     axes.clear()  # alternatively, plt.pause(0.0001)
#
# # make gif from saved png files
# if make_gif:
#     imageio.mimsave(name + '.gif', frames, 'GIF-FI')
#     os.rmdir(directory + '/GIF')
