# plotting functions for
#
# ~******* CASCADE ********~
#
# and for the manuscript titled:
#
# "Lexi's very cool outwash module"
#
"""----------------------------------------------------
Copyright (C) 2022 Lexi Van Blunk
----------------------------------------------------"""

import matplotlib.pyplot as plt
import numpy as np
import os
import imageio
from matplotlib.ticker import AutoMinorLocator

# -------------------------------------------elevation gif--------------------------------------------------------------
def plot_ElevAnimation(dunes, interior, directory, start, stop, freq, berm_el):
    """

    :param dunes: dune domain array from cascade for example cascade_outwash50.barrier3d[0]._DuneDomain
    :param interior: interior array from cascade for example cascade_outwash50.barrier3d[0]._DomainTS
    :param directory: folder where we will save outputs
    :param start: start time index
    :param stop: end time index
    :param freq: interval for plotting the domains
    :param berm_el: berm elevation in dam
    :return:
    """
    if not os.path.exists(directory):
        os.makedirs(directory)
    os.chdir(directory)

    # os.chdir(directory)
    newpath = "Elevations/"
    if not os.path.exists(newpath):
        os.makedirs(newpath)
    os.chdir(newpath)

    for t in range(start, stop, freq):
        dunes_TS = np.transpose(dunes[t]) + berm_el
        interior_TS = interior[t]
        animate_domain = np.vstack([dunes_TS, interior_TS])

        # Plot and save
        plt.rcParams.update({"font.size": 15})
        elevFig1 = plt.figure(figsize=(15, 7))
        ax = elevFig1.add_subplot(111)
        cax = ax.matshow(
            animate_domain * 10,
            cmap="terrain",
            vmin=-3, vmax=6
        )  # , interpolation='gaussian') # analysis:ignore
        ax.xaxis.set_ticks_position("bottom")
        xtick_max = np.shape(animate_domain)[1]  # n_cols = x
        x_ticks = np.array(range(0, xtick_max, 5))
        x_tick_labels = x_ticks * 10
        ytick_max = np.shape(animate_domain)[0]  # n_rows = y
        y_ticks = np.array(range(0, ytick_max, 5))
        y_tick_labels = y_ticks * 10
        plt.xticks(x_ticks, x_tick_labels)
        plt.yticks(y_ticks, y_tick_labels)

        cbar = elevFig1.colorbar(cax)
        cbar.set_label('m MHW', rotation=270, labelpad=15)
        plt.xlabel("Alongshore Distance (m)")
        plt.ylabel("Cross-Shore Distance (m)")
        plt.title("Time = " + str(t) + " years")
        plt.tight_layout()
        # timestr = "Time = " + str(t) + " hrs"
        # plt.text(1, 1, timestr)
        name = "elev_" + str(t)
        elevFig1.savefig(name, facecolor='w')  # dpi=200
        plt.close(elevFig1)

    frames = []

    for filenum in range(start, stop, freq):
        filename = "elev_" + str(filenum) + ".png"
        frames.append(imageio.imread(filename))
    imageio.mimsave("elev.gif", frames, fps=2)
    # imageio.mimsave("elev.gif", frames, "GIF-FI")
    print("[ * elevation GIF successfully generated * ]")


def plot_dune_domain(b3d, TMAX):
    """
    :param b3d: a list containing BMI objects from cascade
    :param TMAX: the last time index for plotting
    :return: dune figure

    Dune Height Over Time for CASCADE

    """
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
        # vmin=0,
        # vmax=Dmax * 10,
    )
    ax.xaxis.set_ticks_position("bottom")  # analysis:ignore
    cax = ax.xaxis.set_ticks_position("bottom")  # analysis:ignore
    cbar = duneFig.colorbar()
    # cbar.set_label('Dune Height Above Berm Elevation (m)', rotation=270)
    plt.xlabel("Alongshore Distance (dam)")
    plt.ylabel("Year")
    plt.title("Dune Height (m)")

    return duneFig


def plot_ElevAnimation_CASCADE(
    cascade,
    directory,
    TMAX_MGMT,
    name,
    TMAX_SIM,
    ny=1,
    beach_management_ny=None,  # list of bool the length of ny, or None for all False
    roadway_management_ny=None,
    y_lim=(150, 250),
    z_lim=6,
    z_bot=0,
    fig_size=None,
    fig_eps=False,
    km_on=True,
):
    """
    :param cascade: a cascade model object
    :param directory: for saving (string)
    :param TMAX_MGMT: last time index that the b3d subgrid was managed (list)
    :param name: for saving (string)
    :param TMAX_SIM: last time index that the b3d subgrid was managed (list)
    :param ny: the number of barrier3d subgrids that you want to plot
    :param beach_management_ny: booleans for b3d subgrids that used
        BeachDuneManager (list)
    :param roadway_management_ny: booleans for b3d subgrids that used
        RoadwayManager (list)
    :param y_lim: y limits [low, high] of plot [dam]
    :param z_lim: z limit of plot [m MHW]
    :param fig_size: size of plot, e.g., (6, 2.5)
    :param fig_eps: output files in eps format -- default is png [best for gif]
    :return: gif

    Animation Frames of Barrier and Dune Elevation (#4 in Barrier3D_Functions,
    modified here for CASCADE)

    NOTE THAT THE BEACH REPRESENTATION IS BASED ON A MODEL SPECIFIED BEACH WIDTH.
    We set the beach width for the remaining time steps after the community has
    been abandoned to the last managed beach width in order to not have a huge
    jump in the back-barrier position in Barrier3D. OTHERWISE, it is meaningless
    to the dynamics Barrier3D.
    """
    if beach_management_ny is None:
        beach_management_ny = [False] * ny
    if roadway_management_ny is None:
        roadway_management_ny = [False] * ny

    barrier3d = cascade.barrier3d

    # set up the domain; here we just use the first grid, but that could break in
    # future runs
    BarrierLength = barrier3d[0].BarrierLength
    if np.any(beach_management_ny):
        indices = [i for i in range(ny) if beach_management_ny[i] == 1]
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
    newpath = name
    if not os.path.exists(newpath):
        os.makedirs(newpath)
    os.chdir(newpath)

    maxMGMT = np.max(TMAX_MGMT)
    plt.rcParams.update({"font.size": 15})

    if np.any(beach_management_ny) or np.any(roadway_management_ny):
        for t in range(maxMGMT + 1):
            # start with plotting t=0, then plot post-storm dune, interior,
            # shoreface, etc. before management, treating this as t=0.5 (i.e.,
            # this is really the final configuration from storms at t=1,2,3,4,5,...
            # but before human modifications); for the default dune inputs, we
            # should be able to see natural dune growth at t=0.5 before rebuild
            # and storms (i.e., t=0.5 and t=1 should not be the same)
            if 0 < t < TMAX_SIM:
                # post-storm variables in the BeachDuneManager are: interior,
                # dune height, x_s, s_sf, beach width
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

                    # Build beach elevation domain, we only show beach width decreasing
                    # in increments of 10 m and we don't illustrate a berm, just a
                    # sloping beach up to the elevation of the berm
                    cellular_dune_toe_post_storm = np.floor(
                        actual_shoreline_post_storm[t] + beach_width
                    )
                    cellular_shoreline_post_storm = np.floor(
                        actual_shoreline_post_storm[t]
                    )
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
                    vmin=z_bot,
                    vmax=z_lim,
                    # edgecolors="w",  # for debugging
                    # linewidth=0.01,
                )
                cbar = elevFig1.colorbar(cax)
                cbar.set_label("elevation (m MHW)", rotation=270)
                plt.xlabel("alongshore distance (dam)")
                plt.ylabel("cross-shore distance (dam)")
                timestr = (
                    "Time = " + str(t - 0.5) + " yrs"
                )  # we are letting the post-storm output represent 0.5 years
                if km_on:
                    locs, _ = plt.yticks()
                    plt.yticks(locs, locs / 100)
                    locs, _ = plt.xticks()
                    plt.xticks(locs[0:-1], locs[0:-1] / 100)
                    plt.xlabel("alongshore distance (km)")
                    plt.ylabel("cross-shore distance (km)")
                    ax.yaxis.set_minor_locator(AutoMinorLocator())
                if y_lim is not None:
                    plt.ylim(y_lim)
                    plt.text(3, y_lim[0] + 3, timestr, color="w")
                else:
                    plt.ylim(bottom=OriginY - 35)
                    plt.text(1, OriginY - 33, timestr)
                plt.tight_layout()
                plt.rcParams.update({"font.size": 15})
                # elevFig1.tight_layout()
                if fig_eps:
                    name = "elev_" + str(t - 1) + "pt5.eps"
                    elevFig1.savefig(name, format="eps")
                else:
                    name = "elev_" + str(t - 1) + "pt5"
                    elevFig1.savefig(name)  # dpi=200
                plt.close(elevFig1)

    for t in range(TMAX_SIM):
        # ok, now the annual time step, which incorporates human modifications to
        # the shoreface, beach, dune, & interior
        AnimateDomain = np.ones([AniDomainWidth + 1, BarrierLength * ny]) * -1

        for iB3D in range(ny):
            actual_shoreline_post_humans = barrier3d[iB3D].x_s_TS[0 : TMAX_SIM + 1]

            if beach_management_ny[iB3D]:
                beach_width = cascade.nourishments[iB3D].beach_width[t] / 10
            else:
                # both roadways and natural scenario
                beach_width = cascade._initial_beach_width[iB3D] / 10
                # beach_width = cascade.outwash[iB3D]._beach_width[t] / 10

            if beach_management_ny[iB3D] and np.isnan(beach_width):
                beach_width = (
                    cascade.nourishments[iB3D].beach_width[t - 1] / 10
                )  # this can happen if the barrier height drowns

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
            Domain = np.flip(np.vstack([Beach, Dunes, Domain]),1)
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
            vmin=z_bot,
            vmax=z_lim,
            # edgecolors="w",  # for debugging
            # linewidth=0.01,
        )
        cbar = elevFig2.colorbar(cax)
        cbar.set_label("elevation (m MHW)", rotation=270)
        plt.xlabel("alongshore distance (dam)")
        plt.ylabel("cross-shore distance (dam)")
        timestr = "Time = " + str(t) + " yrs"
        if km_on:
            locs, _ = plt.yticks()
            plt.yticks(locs, locs / 100)
            locs, _ = plt.xticks()
            plt.xticks(locs[0:-1], locs[0:-1] / 100)
            plt.xlabel("alongshore distance (km)")
            plt.ylabel("cross-shore distance (km)")
            ax.yaxis.set_minor_locator(AutoMinorLocator())
        if y_lim is not None:
            plt.ylim(y_lim)
            plt.text(3, y_lim[0] + 3, timestr, color="w")
        else:
            plt.ylim(bottom=OriginY - 35)
            plt.text(1, OriginY - 33, timestr)
        plt.tight_layout()
        plt.rcParams.update({"font.size": 15})
        # elevFig2.tight_layout()
        if fig_eps:
            name = "elev_" + str(t) + ".eps"
            elevFig2.savefig(name, format="eps")
        else:
            name = "elev_" + str(t)
            elevFig2.savefig(name)  # dpi=200
        plt.close(elevFig2)

    frames = []

    for filenum in range(TMAX_SIM):
        if fig_eps:
            filename = "elev_" + str(filenum) + ".eps"
        else:
            filename = "elev_" + str(filenum) + ".png"
        frames.append(imageio.imread(filename))

        # management simulations
        if np.any(beach_management_ny) or np.any(roadway_management_ny):
            tmax_management = np.array(TMAX_MGMT)
            max_pt5_year = np.max(tmax_management[tmax_management < TMAX_SIM - 1])
            if filenum < max_pt5_year:
                if fig_eps:
                    filename = "elev_" + str(filenum) + "pt5" ".eps"
                else:
                    filename = "elev_" + str(filenum) + "pt5" ".png"
                frames.append(imageio.imread(filename))

    imageio.mimsave("elev.gif", frames, fps=2)
    # imageio.mimsave("elev.gif", frames, "GIF-FI")
    print()
    print("[ * GIF successfully generated * ]")