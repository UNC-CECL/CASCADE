# Plotting functions for

# ~******* CASCADE ********~

"""----------------------------------------------------
Copyright (C) 2020 Katherine Anarde
Full copyright notice located in main CASCADE.py file
----------------------------------------------------"""

# Simulation Functions Included:


import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib import cm
import numpy as np
import pandas as pd
import os
import imageio
import math

# import seaborn as sns
from scipy import signal

# ===================================================
# 1: Animation Frames of Barrier and Dune Elevation (#4 in Barrier3D_Functions, modified here for CASCADE)
#
#       inputs:         - barrier3d (a list containing BMI objects)
#                       - ny (the number of barrier3d subgrids that you want to plot)
#                       - directory (for saving)
#                       - TMAX (the last time index for plotting)
#                       - name (for saving new directory)
#       outputs:        - gif


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


# added a 0.5 time step for post-storm impacts before human modifications
def plot_ElevAnimation_Humans(
    cascade,
    ny,
    directory,
    TMAX,
    name,
):

    barrier3d = cascade.barrier3d

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

        # start with t=0, then post-storm dune and interior before rebuild, treating this as t=0.5 (i.e., this is really the
        # final dune from storms at t=1,2,3,4,5,... but before human modifications to the dune (rebuild, move overwash);
        # for the high dune scenarios, we should be able to see natural dune growth at t=0.5 before rebuild (i.e., t=0.5
        # and t=1 should not be the same)
        if t > 0:
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
                Domain = cascade.post_storm_interior[t][iB3D] * 10
                Dunes = (
                    cascade.post_storm_dunes[t][iB3D] + barrier3d[iB3D]._BermEl
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
            timestr = (
                "Time = " + str(t - 0.5) + " yrs"
            )  # we are letting the post-storm output represent 0.5 years
            plt.text(1, 1, timestr)
            plt.rcParams.update({"font.size": 20})
            name = "elev_" + str(t - 1) + "pt5"
            elevFig1.savefig(name)  # dpi=200
            plt.close(elevFig1)

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
# #3 B3D cross-Shore Transect for one subgrid every 100 m for last time step (#5 in Barrier3D_Functions)
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
# #4: Cross-shore transects for both brie and B3d
#
#       inputs:         - barrier3d (a list of B3D objects)
#                       - brieLTA (a brie object with LTA model on)
#                       - time_step (a list of time steps to plot)
#                       - iB3D (an integer corresponding to the B3D subdomain / brie alongshore grid)
#       outputs:        - fig


def plot_ModelTransects(b3d, brieLTA, time_step, iB3D):

    plt.figure(figsize=(10, 5))
    colors = mpl.cm.viridis(np.linspace(0, 1, b3d[0]._TMAX))

    plt.subplot(2, 1, 1)

    for t in time_step:

        # Sea level
        SL = b3d[iB3D]._SL + (t * b3d[iB3D]._RSLR[t])

        # Create data points
        Tx = (
            b3d[iB3D]._x_t_TS[t] - b3d[iB3D]._x_t_TS[0]
        )  # subtract start position so starts at zero
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


# ===================================================
# #5: Statistics from B3D
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
# #6: Statistics from BRIE with LTA
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
# #7: Calculate shoreline change periodicity
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
# #8: Dune Height Over Time for CASCADE
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
# #9: Average shoreline change rate over time for CASCADE
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
# #10: Differences in punctuated retreat for CASCADE vs B3D (AST model vs no AST)
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


def plot_nonlinear_stats(
    CASCADE_b3d,
    ib3d,
    tmin,
    tmax,
    post_storm_dunes=None,
    post_storm_interior=None,
    design_height=None,
    rebuild_threshold=None,
):

    # mean dune height (using both dune domain columns) ---------
    # DuneCrestMean = [
    #     (a + CASCADE_b3d[ib3d].BermEl) * 10 for a in CASCADE_b3d[ib3d]._Hd_AverageTS
    # ]  # relative to NAVD88
    DuneDomainCrest = (
        CASCADE_b3d[ib3d]._DuneDomain[:, :, :].max(axis=2)
    )  # Maximum height of each row in DuneDomain
    # DuneDomainCrest = CASCADE_b3d[ib3d]._DuneDomain[
    #     :, :, 0
    # ]  # Just front row of dune domain
    DuneDomainCrest[DuneDomainCrest < CASCADE_b3d[ib3d]._DuneRestart] = CASCADE_b3d[
        ib3d
    ]._DuneRestart
    DuneCrestMean = (np.mean(DuneDomainCrest, axis=1) + CASCADE_b3d[ib3d].BermEl) * 10
    DuneCrestMin = (np.min(DuneDomainCrest, axis=1) + CASCADE_b3d[ib3d].BermEl) * 10
    DuneCrestMax = (np.max(DuneDomainCrest, axis=1) + CASCADE_b3d[ib3d].BermEl) * 10

    if post_storm_dunes is not None:

        post_DuneCrestMean = [None]
        post_DuneCrestMin = [None]
        post_DuneCrestMax = [None]

        # same calculation as Hd_AverageTS, but here average post-storm dune-height for each time step
        for t in range(1, CASCADE_b3d[ib3d].time_index):
            DuneDomainCrest = post_storm_dunes[t][ib3d].max(
                axis=1
            )  # Maximum height of each row in DuneDomain (typically front row: could also just force this to be front row, more robust)
            DuneDomainCrest[
                DuneDomainCrest < CASCADE_b3d[ib3d]._DuneRestart
            ] = CASCADE_b3d[ib3d]._DuneRestart
            post_DuneCrestMean.append(
                (np.mean(DuneDomainCrest) + CASCADE_b3d[ib3d].BermEl) * 10
            )  # relative to NAVD88
            post_DuneCrestMax.append(
                (np.max(DuneDomainCrest) + CASCADE_b3d[ib3d].BermEl) * 10
            )  # relative to NAVD88
            post_DuneCrestMin.append(
                (np.min(DuneDomainCrest) + CASCADE_b3d[ib3d].BermEl) * 10
            )  # relative to NAVD88

    # barrier width
    BarrierWidth = (
        np.array(CASCADE_b3d[ib3d].x_b_TS) - np.array(CASCADE_b3d[ib3d].x_s_TS)
    ) * 10

    # change in barrier width
    bwts = [(x - BarrierWidth[0]) for x in BarrierWidth]
    rate = [0]
    for k in range(1, len(bwts)):
        rate.append(bwts[k] - bwts[k - 1])
    bw_rate = (
        rate  # note, np.diff doesn't start with zero rate of change, so we do this calc
    )

    # average interior height
    BarrierHeight = []
    for t in range(len(CASCADE_b3d[ib3d]._DomainTS)):
        bh_array = np.array(CASCADE_b3d[ib3d]._DomainTS[t]) * 10
        BarrierHeight.append(bh_array[bh_array > 0].mean())

    # average interior height post-storm (before human modifications)
    if post_storm_interior is not None:
        post_BarrierHeight = [None]

        for t in range(1, CASCADE_b3d[ib3d].time_index):
            bh_array = np.array(post_storm_interior[t][ib3d]) * 10
            post_BarrierHeight.append(bh_array[bh_array > 0].mean())

    # change in barrier height
    bhts = [(x - BarrierHeight[0]) for x in BarrierHeight]
    rate = [0]
    for k in range(1, len(bhts)):
        rate.append(bhts[k] - bhts[k - 1])
    bh_rate = (
        rate  # note, np.diff doesn't start with zero rate of change, so we do this calc
    )

    # # shoreline change rate
    # scts = [(x - CASCADE_b3d[ib3d]._x_s_TS[0]) * 10 for x in CASCADE_b3d[ib3d]._x_s_TS]
    # rate = [0]
    # for k in range(1, len(scts)):
    #     rate.append(scts[k] - scts[k - 1])
    # sc_rate = rate

    # overwash flux
    Qoverwash = CASCADE_b3d[ib3d].QowTS

    plt.figure(figsize=(10, 5))

    # individual time series
    plt.subplot(2, 5, 1)

    if post_storm_dunes is not None:
        # time = np.arange(0.5, CASCADE_b3d[ib3d].time_index - 1, 1)
        # plt.scatter(time, post_DuneCrestMean[1:])
        time = np.arange(0, CASCADE_b3d[ib3d].time_index - 0.5, 0.5)
        combined = [None] * (len(time))
        combined[::2] = DuneCrestMin
        combined[1::2] = post_DuneCrestMin[1:]  # post_DuneCrestMean[1:]
        plt.hlines(rebuild_threshold, time[0], time[-1], colors="red")
        plt.hlines(design_height, time[0], time[-1], colors="black")
        plt.hlines(CASCADE_b3d[ib3d]._Dmaxel * 10, time[0], time[-1], colors="green")
        plt.plot(time, combined)
        plt.legend(["Min dune height", "rebuild", "design", "max-equil"])

    else:
        plt.plot(DuneCrestMean)
    # plt.plot(DuneCrestMean[tmin:tmax], "m")
    plt.ylabel("Ave Dune Height (m NAVD88)")
    plt.xlabel("Time (yr)")

    # plt.subplot(2, 5, 2)
    # # plt.subplot(5, 2, 3)
    # plt.plot(sc_rate)
    # # plt.plot(sc_rate[tmin:tmax], "m")
    # plt.ylabel("Shoreline change rate (m/yr)")
    # plt.xlabel("Time (yr)")

    plt.subplot(2, 5, 2)
    plt.plot(BarrierHeight)
    plt.ylabel("Ave barrier interior height (m NAVD88)")
    plt.xlabel("Time (yr)")

    plt.subplot(2, 5, 3)
    # plt.subplot(5, 2, 5)
    plt.plot(BarrierWidth)
    # plt.plot(BarrierWidth[tmin:tmax], "m")
    plt.ylabel("Barrier width (m)")
    plt.xlabel("Time (yr)")

    # plt.subplot(2, 5, 4)
    # # plt.subplot(5, 2, 7)
    # plt.plot(Qoverwash)
    # # plt.plot(Qoverwash[tmin:tmax], "m")
    # plt.ylabel("Overwash flux ($m^3/m$)")
    # plt.xlabel("Time (yr)")

    plt.subplot(2, 5, 4)
    # plt.subplot(5, 2, 7)
    plt.plot(bh_rate)
    # plt.plot(Qoverwash[tmin:tmax], "m")
    plt.ylabel("Barrier height change rate (m/yr)")
    plt.xlabel("Time (yr)")

    plt.subplot(2, 5, 5)
    # plt.subplot(5, 2, 9)
    plt.plot(bw_rate)
    plt.ylabel("Barrier width change rate (m/yr)")
    plt.xlabel("Time (yr)")

    # plt.subplot(2, 5, 5)
    # plt.plot(np.array(CASCADE_b3d[ib3d]._x_s_TS) * 10)
    # # plt.plot(np.array(CASCADE_b3d[ib3d]._x_s_TS[tmin:tmax]) * 10, "m")
    # plt.ylabel("Shoreline position (m)")
    # plt.xlabel("Time (yr)")

    # # dune height vs shoreline position
    # plt.subplot(2, 5, 6)
    # # plt.plot(DuneCrestMean, np.array(CASCADE_b3d[ib3d]._x_s_TS) * 10)
    # # plt.plot(
    # #     DuneCrestMean[tmin:tmax],
    # #     np.array(CASCADE_b3d[ib3d]._x_s_TS[tmin:tmax]) * 10,
    # #     "m",
    # # )
    # plt.scatter(
    #     DuneCrestMean,
    #     np.array(CASCADE_b3d[ib3d]._x_s_TS) * 10,
    #     c=np.arange(0, np.size(DuneCrestMean), 1),
    #     cmap="Greens",
    # )
    # plt.xlabel("Ave Dune Height (m)")
    # plt.ylabel("Shoreline position (m)")

    # # dune height vs shoreline change rate
    # plt.subplot(2, 5, 7)
    # # plt.subplot(5, 2, 4)
    # # plt.scatter(DuneCrestMean, sc_rate)
    # # plt.plot(DuneCrestMean[tmin:tmax], sc_rate[tmin:tmax], "m")
    # plt.scatter(
    #     DuneCrestMean,
    #     sc_rate,
    #     c=np.arange(0, np.size(DuneCrestMean), 1),
    #     cmap=cm.viridis,  # "Greens"
    # )
    # plt.xlabel("Ave Dune Height (m)")
    # plt.ylabel("Shoreline change rate (m/yr)")

    # # dune height vs overwash flux
    # plt.subplot(2, 5, 8)
    # # plt.subplot(5, 2, 6)
    # # plt.scatter(DuneCrestMean, Qoverwash)
    # # plt.plot(DuneCrestMean[tmin:tmax], Qoverwash[tmin:tmax], "m")
    # plt.scatter(
    #     DuneCrestMean,
    #     Qoverwash,
    #     c=np.arange(0, np.size(DuneCrestMean), 1),
    #     cmap=cm.viridis,
    # )
    # plt.xlabel("Ave Dune Height (m)")
    # plt.ylabel("Qoverwash ($m^3/m$)")

    # dune height vs barrier height
    plt.subplot(2, 5, 6)
    plt.scatter(
        DuneCrestMean,
        BarrierHeight,
        c=np.arange(0, np.size(DuneCrestMean), 1),
        cmap=cm.viridis,
    )
    # plt.plot(DuneCrestMean[tmin:tmax], BarrierWidth[tmin:tmax], "m")
    plt.xlabel("Ave Dune Height (m)")
    plt.ylabel("Ave barrier interior height (m)")

    # dune height vs barrier width
    plt.subplot(2, 5, 7)
    # plt.subplot(5, 2, 8)
    # plt.plot(DuneCrestMean, BarrierWidth)
    plt.scatter(
        DuneCrestMean,
        BarrierWidth,
        c=np.arange(0, np.size(DuneCrestMean), 1),
        cmap=cm.viridis,
    )
    # plt.plot(DuneCrestMean[tmin:tmax], BarrierWidth[tmin:tmax], "m")
    plt.xlabel("Ave Dune Height (m)")
    plt.ylabel("Barrier width (m)")

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

    plt.tight_layout()

    # # now try my hand at a 3d plot
    # from mpl_toolkits import mplot3d
    #
    # plt.figure()
    # ax = plt.axes(projection="3d")
    #
    # # three-dimensional line
    # # ax.plot3D(DuneCrestMean, BarrierWidth, sc_rate, "gray")
    # ax.plot3D(
    #     DuneCrestMean[tmin:tmax], BarrierWidth[tmin:tmax], sc_rate[tmin:tmax], "gray"
    # )
    #
    # # three-dimensional scattered points
    # ax.scatter3D(
    #     DuneCrestMean[tmin:tmax],
    #     BarrierWidth[tmin:tmax],
    #     sc_rate[tmin:tmax],
    #     c=sc_rate[tmin:tmax],
    #     cmap="Greens",
    # )
    # ax.set_title(str(tmin) + "-" + str(tmax) + " yrs")
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
    # ax.plot3D(Qoverwash[tmin:tmax], BarrierWidth[tmin:tmax], sc_rate[tmin:tmax], "gray")
    #
    # ax.scatter3D(
    #     Qoverwash[tmin:tmax],
    #     BarrierWidth[tmin:tmax],
    #     sc_rate[tmin:tmax],
    #     c=sc_rate[tmin:tmax],
    #     cmap="Greens",
    # )
    # ax.set_title(str(tmin) + "-" + str(tmax) + " yrs")
    # ax.set_xlabel("Overwash flux ($m^3/m$)")
    # ax.set_ylabel("Barrier width (m)")
    # ax.set_zlabel("Shoreline change rate (m/yr)")
    #
    # # try Qoverwash, dune height, and barrier width change
    # fig = plt.figure()
    # ax = plt.axes(projection="3d")
    #
    # ax.plot3D(
    #     DuneCrestMean[tmin:tmax], Qoverwash[tmin:tmax], bw_rate[tmin:tmax], "gray"
    # )
    #
    # # three-dimensional scattered points
    # ax.scatter3D(
    #     DuneCrestMean[tmin:tmax],
    #     Qoverwash[tmin:tmax],
    #     bw_rate[tmin:tmax],
    #     c=bw_rate[tmin:tmax],
    #     cmap="Greens",
    # )
    # ax.set_title(str(tmin) + "-" + str(tmax) + " yrs")
    # ax.set_xlabel("Dune height (m)")
    # ax.set_ylabel("Overwash flux ($m^3/m$)")
    # ax.set_zlabel("Barrier width change rate (m/yr)")

    return BarrierWidth, DuneCrestMean, BarrierHeight


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
