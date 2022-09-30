# plotting functions for
#
# ~******* CASCADE ********~
#
"""----------------------------------------------------
Copyright (C) 2022 Katherine Anarde
----------------------------------------------------"""

import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib import cm
import numpy as np
import pandas as pd
import os
import imageio
import math
from scipy import signal

# # ###############################################################################
# # modified from barrier3d
# # ###############################################################################

def plot_ShorelinePositions(x_s_TS, x_b_TS):
    """
    :param x_s_TS: shoreline time series for one b3d subgrid
    :param x_b_TS: back-barrier time series for on b3d subgrid
    :return: fig

    Shoreline positions over time (#6 in Barrier3D_Functions)

    """
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

    return fig

def plot_XShoreTransects(barrier3d, TMAX):
    """
    :param barrier3d: a list containing BMI objects from cascade
    :param TMAX: the last time index for plotting
    :return: fig

    B3D cross-shore transect for one subgrid every 100 m for last time step (#5 in Barrier3D_Functions)

    """

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

    return fig

def plot_ModelTransects(cascade, time_step, iB3D):
    """
    :param cascade: cascade object
    :param time_step: a list of time steps to plot
    :param iB3D: an integer corresponding to the B3D subdomain
    :return: fig

    4: Cross-shore transects from one B3D model (from CASCADE), which includes beach width

    Function plots model transects over time. Note that the barrier interior narrows due to SLR, but also when the
    shoreline erodes one or more full dam cells and the dunes are forced to migrate into the interior (i.e., one or more
    barrier interior cells become the dune line). Hence, because this plotter shows the migration of the shoreline for
    each time step (non-integer multiples) and always two rows of dunes, it appears that the barrier island is moving
    landward when really it is eroding the dune line.

    """
    fig, axs = plt.subplots(1, 1, figsize=(5, 3), sharey=True, sharex=True)
    legend_t = []

    for t in time_step:

        # Sea level
        sea_level = cascade.barrier3d[iB3D]._SL + (
            t * cascade.barrier3d[iB3D]._RSLR[t]
        )  # dam

        # Create data points
        shoreface_toe_x = (
            cascade.barrier3d[iB3D].x_t_TS[t] - cascade.barrier3d[iB3D].x_t_TS[0]
        )
        # shoreface_toe_x = cascade.barrier3d[iB3D].x_t_TS[t]  # dam
        shoreface_toe_y = (sea_level - cascade.barrier3d[iB3D].DShoreface) * 10  # m
        shoreline_x = (
            cascade.barrier3d[iB3D].x_s_TS[t] - cascade.barrier3d[iB3D].x_t_TS[0]
        )  # dam
        # shoreline_x = cascade.barrier3d[iB3D].x_s_TS[t]  # dam
        shoreline_y = sea_level * 10  # m
        bay_y = (sea_level - cascade.barrier3d[iB3D]._BayDepth) * 10  # m
        end_of_bay_y = bay_y

        # if cascade.nourishments[iB3D].beach_width[t] is None:
        if np.isnan(cascade.nourishments[iB3D].beach_width[t]):
            berm_x = shoreline_x + (
                int(cascade.barrier3d[iB3D].BermEl / cascade.barrier3d[iB3D]._beta)
            )  # initial beach width (in dam)
        else:
            berm_x = shoreline_x + (
                cascade.nourishments[iB3D].beach_width[t] / 10
            )  # beach width (in dam)
        berm_y = (
            cascade.barrier3d[iB3D]._BermEl * 10
        ) + shoreline_y  # convert to meters
        dune_toe_x = berm_x
        dune_toe_y = berm_y

        v = 10  # just use 10th transect
        interior_y = cascade.barrier3d[iB3D]._DomainTS[t]  # dam MHW
        interior_y = interior_y[:, v]
        dunes_y = (
            cascade.barrier3d[iB3D]._DuneDomain[t, v, :]
            + cascade.barrier3d[iB3D]._BermEl
        )  # dam MHW
        cross_barrier_y = np.insert(interior_y, 0, dunes_y)
        cross_barrier_y = (
            cross_barrier_y * 10
        ) + shoreline_y  # Convert to meters, with SLR included
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
        legend_t.append("year " + str(t))

    plt.xlabel("cross-shore distance (dam)")
    plt.ylabel("elevation (m MHW)")
    plt.legend(legend_t)
    plt.tight_layout()

    return fig

def calc_ShorelinePeriodicity(x_s_TS):
    """
    :param x_s_TS: a single shoreline time series for one B3D model
    :return: [
    Periodicity (period of punctuated retreat),
    AvgFastDur (average duration of fast periods),
    AvgSlowDur (average duration of slow periods),
    Punc (boolean for punctuated retreat)
    ]

    Calculate shoreline change periodicity

    """
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

    # CALCULATE STATS
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

    return [Periodicity, AvgFastDur, AvgSlowDur, Punc]

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
        vmin=0,
        vmax=Dmax * 10,
    )
    ax.xaxis.set_ticks_position("bottom")  # analysis:ignore
    # cax = ax.xaxis.set_ticks_position("bottom")  # analysis:ignore
    # cbar = duneFig.colorbar(cax)
    # cbar.set_label('Dune Height Above Berm Elevation (m)', rotation=270)
    plt.xlabel("Alongshore Distance (dam)")
    plt.ylabel("Year")
    plt.title("Dune Height (m)")

    return duneFig

def plot_ShorelineChangeRate(b3d):
    """
    :param b3d: list containing BMI objects from cascade
    :return: figure

    Average shoreline change rate over time for CASCADE

    """
    ave_rate = []

    for iB3D in range(len(b3d)):
        scts = [(x - b3d[iB3D]._x_s_TS[0]) * 10 for x in b3d[iB3D]._x_s_TS]
        rate = [0]
        for k in range(1, len(scts)):
            rate.append(scts[k] - scts[k - 1])
        ave_rate.append(rate)

    df = pd.DataFrame(data=ave_rate)
    ave = df.mean(axis=0)

    plt.figure()
    plt.plot(ave)
    fig = plt.gcf()
    fig.set_size_inches(14, 5)
    plt.xlabel("Year")
    plt.ylabel("Shoreline Erosion Rate (m/yr)")
    plt.show()

    return fig

def plot_punctuated_difference(CASCADE_b3d, b3d_only, ny):
    """
    :param CASCADE_b3d: a list containing BMI objects from a CASCADE run
    :param b3d_only: a list of BMI objects from ind b3d runs, corr. to the growth rates used in the CASCADE model above
    :param ny: number of alongshore cells
    :return: figure

    Differences in punctuated retreat for CASCADE vs B3D (AST model vs no AST)

    """

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

# # ###############################################################################
# # new to cascade
# # ###############################################################################

def plot_ElevAnimation_CASCADE(
        cascade,
        directory,
        TMAX_MGMT,
        name,
        TMAX_SIM,
        ny=1,
        beach_management_ny=[False],  # list of booleans the length of ny
        roadway_management_ny=[False],
        y_lim=[150, 250],
        z_lim=3.5,
        fig_size=None,
        fig_eps=False,
):
    """
    :param cascade: a cascade model object
    :param directory: for saving (string)
    :param TMAX_MGMT: last time index that the b3d subgrid was managed (list)
    :param name: for saving (string)
    :param TMAX_SIM: last time index that the b3d subgrid was managed (list)
    :param ny: the number of barrier3d subgrids that you want to plot
    :param beach_management_ny: booleans for b3d subgrids that used BeachDuneManager (list)
    :param roadway_management_ny: booleans for b3d subgrids that used RoadwayManager (list)
    :param y_lim: y limits [low, high] of plot [dam]
    :param z_lim: z limit of plot [m MHW]
    :param fig_size: size of plot, e.g., (6, 2.5)
    :param fig_eps: output files in eps format -- default is png [best for gif]
    :return: gif

    Animation Frames of Barrier and Dune Elevation (#4 in Barrier3D_Functions, modified here for CASCADE)

    NOTE THAT THE BEACH REPRESENTATION IS BASED ON A MODEL SPECIFIED BEACH WIDTH. We set the beach width for the
    remaining time steps after the community has been abandoned to the last managed beach width in order to not have a
    huge jump in the back-barrier position in Barrier3D. OTHERWISE, it is meaningless to the dynamics Barrier3D.

    """
    barrier3d = cascade.barrier3d

    # set up the domain; here we just use the first grid, but that could break in future runs
    BarrierLength = barrier3d[0].BarrierLength
    if np.any(beach_management_ny):
        indices = [i for i in range(ny) if beach_management_ny[i] == 1]
        iB3D = indices[0]
        MaxBeachWidth = (
                np.max(cascade.nourishments[iB3D].beach_width[0: TMAX_MGMT[iB3D]]) / 10
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

    if np.any(beach_management_ny) or np.any(roadway_management_ny):
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
                                1: TMAX_MGMT[iB3D] + 1
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
                                                      0: TMAX_MGMT[iB3D] + 1
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
                    OriginTstart:OriginTstop, xOrigin: xOrigin + BarrierLength
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
                if y_lim is not None:
                    plt.ylim(y_lim)
                    plt.text(3, y_lim[0] + 3, timestr, color="w")
                else:
                    plt.ylim(bottom=OriginY - 35)
                    plt.text(1, OriginY - 33, timestr)
                plt.tight_layout()
                plt.rcParams.update({"font.size": 11})
                # elevFig1.tight_layout()
                if fig_eps:
                    name = "elev_" + str(t - 1) + "pt5.eps"
                    elevFig1.savefig(name, format="eps")
                else:
                    name = "elev_" + str(t - 1) + "pt5"
                    elevFig1.savefig(name)  # dpi=200
                plt.close(elevFig1)

    for t in range(TMAX_SIM):

        # ok, now the annual time step, which incorporates human modifications to the shoreface, beach, dune, & interior
        AnimateDomain = np.ones([AniDomainWidth + 1, BarrierLength * ny]) * -1

        for iB3D in range(ny):

            actual_shoreline_post_humans = barrier3d[iB3D].x_s_TS[0: TMAX_SIM + 1]

            if beach_management_ny[iB3D]:
                beach_width = cascade.nourishments[iB3D].beach_width[t] / 10
            else:
                # both roadways and natural scenario
                beach_width = cascade._initial_beach_width[iB3D] / 10

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
            Domain = np.vstack([Beach, Dunes, Domain])
            Domain[Domain < 0] = -1
            widthTS = len(Domain)
            OriginTstart = int(cellular_shoreline_post_humans)
            OriginTstop = OriginTstart + widthTS
            xOrigin = iB3D * BarrierLength
            AnimateDomain[
            OriginTstart:OriginTstop, xOrigin: xOrigin + BarrierLength
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
            vmax=z_lim,
            # edgecolors="w",  # for debugging
            # linewidth=0.01,
        )
        cbar = elevFig2.colorbar(cax)
        cbar.set_label("elevation (m MHW)", rotation=270)
        plt.xlabel("alongshore distance (dam)")
        plt.ylabel("cross-shore distance (dam)")
        timestr = "Time = " + str(t) + " yrs"
        if y_lim is not None:
            plt.ylim(y_lim)
            plt.text(3, y_lim[0] + 3, timestr, color="w")
        else:
            plt.ylim(bottom=OriginY - 35)
            plt.text(1, OriginY - 33, timestr)
        plt.tight_layout()
        plt.rcParams.update({"font.size": 11})
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

        # find the maximum 0.5 year
        tmax_management = np.array(TMAX_MGMT)
        max_pt5_year = np.max(tmax_management[tmax_management < TMAX_SIM - 1])
        if (np.any(beach_management_ny) or np.any(roadway_management_ny)) and (
                filenum < max_pt5_year
        ):
            if fig_eps:
                filename = "elev_" + str(filenum) + "pt5" ".eps"
            else:
                filename = "elev_" + str(filenum) + "pt5" ".png"
            frames.append(imageio.imread(filename))

    imageio.mimsave("elev.gif", frames, "GIF-FI")
    print()
    print("[ * GIF successfully generated * ]")

def combine_post_storm_human_time_series(
    tmax_sim, post_storm_statistic, human_modified_statistic
):
    """
    :param tmax_sim: max simulation time
    :param post_storm_statistic: whatever post-storm (0.5 yr time step) statistic you want to combine
    :param human_modified_statistic: whatever post-human response (1 yr time step) statistic you want to combine
    :return: the combined [0, 0.5, 1, 1.5 ...] yr statistics

    Combines the immediate post storm, and human recovery effort, time series into a single array

    """
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
    post_storm_ave_interior_height=None,
    design_height=None,
    rebuild_threshold=None,
    road_elevation=None,
    dunes_rebuilt=None,
    road_relocated=None,
):
    """
    :param CASCADE_b3d: a list of B3D objects
    :param ib3d: an integer corresponding to the B3D subdomain / brie alongshore grid that you want statistics for
    :param tmax_roadways: time step where roadway is abandoned
    :param tmax_sim: time step where simulation ends
    :param post_storm_dunes: post-storm dune heights from RoadwayManager -- before human mods
    :param post_storm_ave_interior_height: post-storm ave interior height from RoadwayManager
    :param design_height: dune design height
    :param rebuild_threshold: threshold for rebuilding the dune
    :param road_elevation: road elevation time series (in m MHW)
    :param dunes_rebuilt: when the dunes were rebuilt (boolean time series)
    :param road_relocated: when the road was relocated (boolean time series)
    :return: (BarrierWidth, DuneCrestMean, BarrierHeight, bh_rate, bw_rate, sc_rate, DuneCrestMin,DuneCrestMax,
        shoreline_position, shoreface_slope, overwash)

    Statistics from the human dynamics modules; combines time series of human modified variables such that
    the 0.5 year time step represents post-storm morphology and pre-human modifications to dune and barrier height
    variables for the RoadwayManager Module

    """

    # if the post-storm variables are not supplied (essentially a 0.5 year time step), then only the human-modified
    # statistics are plotted (the 1 year time step)

    # variables that need to be combined and plotted: dune height, barrier height/interior
    empty_nans = np.empty(tmax_sim - tmax_roadways)
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
            post_storm_DuneCrestMax.append(
                (np.max(DuneDomainCrest) + CASCADE_b3d[ib3d].BermEl) * 10
            )  # relative to MHW
            post_storm_DuneCrestMin.append(
                (np.min(DuneDomainCrest) + CASCADE_b3d[ib3d].BermEl) * 10
            )  # relative to MHW

        # post_DuneCrestMean = np.concatenate((post_DuneCrestMean, empty_nans))
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
    # not x_b - x_s, which incorporates changes in beach width and the dune line in the nourishment module
    BarrierWidth = (
        np.array(CASCADE_b3d[ib3d].InteriorWidth_AvgTS[0:tmax_sim])
    ) * 10  # m

    # change in barrier width
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
            np.array(post_storm_ave_interior_height[0:tmax_roadways])
        ) * 10  # m

        post_storm_BarrierHeight = np.concatenate(
            (post_storm_BarrierHeight, empty_nans)
        )
        combined_BarrierHeight = combine_post_storm_human_time_series(
            tmax_sim, post_storm_BarrierHeight, BarrierHeight
        )  # m

    # change in barrier height
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

    # overwash flux --------------------------------------------------
    overwash = np.array(CASCADE_b3d[ib3d].QowTS[0:tmax_sim])  # m^3/m

    # shoreface slope -----------------------------------------
    shoreface_slope = CASCADE_b3d[ib3d].s_sf_TS[0:tmax_sim]
    equilibrium_slope = CASCADE_b3d[ib3d]._s_sf_eq

    # individual time series ------------------------------------
    plt.figure(figsize=(10, 8))
    # plt.subplot(4, 3, 1)

    full_time = np.arange(0, tmax_sim - 0.5, 0.5)

    # dunes
    plt.subplot(2, 2, 1)
    if post_storm_dunes is not None:
        time_mgmt = np.arange(0, tmax_roadways, 1)
        mask1 = np.isfinite(combined_DuneCrestMin)
        mask2 = np.isfinite(combined_DuneCrestMax)
        plt.plot(full_time[mask1], combined_DuneCrestMin[mask1])
        plt.plot(full_time[mask2], combined_DuneCrestMax[mask2])
        plt.plot(time_mgmt, rebuild_threshold[0:tmax_roadways], color="green")
        plt.plot(time_mgmt, design_height[0:tmax_roadways], color="red")
        plt.hlines(
            CASCADE_b3d[ib3d]._Dmaxel * 10, full_time[0], full_time[-1], colors="black"
        )
        plt.plot(time_mgmt, road_elevation[0:tmax_roadways], color="purple")
        plt.legend(["min", "max", "rebuild", "design", "road", "max-equil"])

    else:
        plt.plot(DuneCrestMin)
        plt.plot(DuneCrestMax)
        plt.legend(["min", "max"])
    plt.ylabel("dune elevation (m MHW)")
    plt.xlabel("Time (yr)")

    # when are dunes rebuilt
    plt.subplot(2, 2, 2)
    if dunes_rebuilt is not None:
        plt.plot(dunes_rebuilt[0:tmax_sim], "k")
        plt.ylabel("Dunes Rebuilt")
        plt.xlabel("Time (yr)")

    # road relocated
    plt.subplot(2, 2, 3)
    if road_relocated is not None:
        plt.plot(road_relocated[0:tmax_sim], "k")
        plt.ylabel("Road Relocated")
        plt.xlabel("Time (yr)")

    plt.tight_layout()

    plt.figure(figsize=(10, 8))

    # shoreline position
    plt.subplot(2, 2, 1)
    plt.plot(shoreline_position, "k")
    plt.ylabel("Cross-shore position (m)")
    plt.xlabel("Time (yr)")

    # shoreline change rate
    plt.subplot(2, 2, 2)
    plt.plot(sc_rate, "k")
    plt.ylabel("Shoreline change rate (m/yr)")
    plt.xlabel("Time (yr)")

    plt.tight_layout()

    plt.figure(figsize=(10, 8))

    # barrier height
    plt.subplot(2, 2, 1)
    if post_storm_ave_interior_height is not None:
        mask = np.isfinite(combined_BarrierHeight)
        plt.plot(full_time[mask], combined_BarrierHeight[mask], "m")
        # plt.legend(["includes post-storm", "mgmt only"])
    plt.plot(BarrierHeight, "k")
    plt.ylabel("barrier height (m MHW)")
    plt.xlabel("Time (yr)")

    # barrier width
    plt.subplot(2, 2, 2)
    plt.plot(BarrierWidth, "k")
    plt.ylabel("barrier width (m)")
    plt.xlabel("Time (yr)")

    # overwash flux
    plt.subplot(2, 2, 3)
    plt.plot(overwash, "k")
    plt.ylabel("Overwash flux ($m^3/m$)")
    plt.xlabel("Time (yr)")

    # shoreface slope
    plt.subplot(2, 2, 4)
    plt.plot(shoreface_slope, "k")
    plt.hlines(equilibrium_slope, full_time[0], full_time[-1], colors="green")
    plt.ylabel("shoreface slope")
    plt.xlabel("Time (yr)")

    plt.tight_layout()

    if post_storm_dunes is not None:
        DuneCrestMin = combined_DuneCrestMin
        DuneCrestMax = combined_DuneCrestMax

    if post_storm_ave_interior_height is not None:
        BarrierHeight = combined_BarrierHeight

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
    """

    :param CASCADE_b3d: a list of B3D objects
    :param ib3d: an integer corresponding to the B3D subdomain / brie alongshore grid that you want statistics for
    :param tmax_management: time step when community management ends
    :param tmax_sim: time step where simulation ends
    :param nourishments: a list of nourishment objects (i.e., from the BeachDuneManager)
    :param post_storm_dunes: post-storm dune heights from BeachDuneManager -- before human mods
    :param post_storm_x_s: post-storm shoreline positions from BeachDuneManager
    :param post_storm_s_sf: post-storm shoreface slope ...
    :param post_storm_ave_interior_width: post-storm average interior width ...
    :param post_storm_ave_interior_height: post-storm average interior height ...
    :param post_storm_beach_width: post-storm beach width ...
    :param post_storm_Qow: post-storm overwash flux ...
    :param design_elevation: dune design height
    :param rebuild_threshold: dune rebuild threshold
    :param dunes_rebuilt: when the dunes were rebuilt
    :return: (BarrierWidth, DuneCrestMean, BarrierHeight, bh_rate, bw_rate, sc_rate, DuneCrestMin, DuneCrestMax,
        shoreline_position, shoreface_slope, beach_width, overwash, dune_toe)

    Statistics from the human dynamics modules; combines time series of human modified variables such that
    the 0.5 year time step represents post-storm morphology and pre-human modifications to dune height,
    barrier height, barrier width, x_s, s_sf, beach width for the BeachDuneManager Module

    """

    # if the post-storm variables are not supplied (essentially a 0.5 year time step), then only the human-modified
    # statistics are plotted (the 1-year time step)

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

    # individual time series ------------------------------------
    plt.figure(figsize=(10, 8))
    # plt.subplot(4, 3, 1)

    full_time = np.arange(0, tmax_sim - 0.5, 0.5)

    # dunes
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
    plt.subplot(2, 2, 2)
    if dunes_rebuilt is not None:
        plt.plot(dunes_rebuilt, "k")
        plt.ylabel("Dunes Rebuilt")
        plt.xlabel("Time (yr)")

    # nourishment volumes
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

    # beach width
    plt.subplot(2, 2, 1)
    if post_storm_beach_width is not None:
        mask = np.isfinite(combined_beach_width)
        plt.plot(full_time[mask], combined_beach_width[mask], "m")
        plt.plot(beach_width, "k")
        # plt.legend(["includes post-storm", "mgmt only"])
    else:
        plt.plot(beach_width, "k")
    plt.ylabel("Beach width (m)")
    plt.xlabel("Time (yr)")
    plt.xlim([0, tmax_sim])

    # shoreline position
    plt.subplot(2, 2, 2)
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

    # shoreline change rate
    plt.subplot(2, 2, 3)
    plt.plot(sc_rate, "k")
    plt.ylabel("Shoreline change rate (m/yr)")
    plt.xlabel("Time (yr)")

    plt.tight_layout()

    plt.figure(figsize=(10, 8))

    # barrier height
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
    if post_storm_BarrierWidth is not None:
        mask = np.isfinite(combined_BarrierWidth)
        plt.plot(full_time[mask], combined_BarrierWidth[mask], "m")
        # plt.legend(["includes post-storm", "mgmt only"])
    plt.plot(BarrierWidth, "k")
    plt.ylabel("barrier width (m)")
    plt.xlabel("Time (yr)")

    # overwash flux
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
