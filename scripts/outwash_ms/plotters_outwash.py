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



# -------------------------------------------elevation gif--------------------------------------------------------------
def plot_ElevAnimation(elev, directory, start, stop):
    if not os.path.exists(directory):
        os.makedirs(directory)
    os.chdir(directory)

    # os.chdir(directory)
    newpath = "Elevations/"
    if not os.path.exists(newpath):
        os.makedirs(newpath)
    os.chdir(newpath)

    for t in range(start, stop):
        AnimateDomain = elev[t]

        # Plot and save
        elevFig1 = plt.figure(figsize=(15, 7))
        ax = elevFig1.add_subplot(111)
        cax = ax.matshow(
            AnimateDomain,
            # origin="upper",
            cmap="seismic",
            # vmin=-0.000002, vmax=0.000002
        )  # , interpolation='gaussian') # analysis:ignore
        ax.xaxis.set_ticks_position("bottom")
        elevFig1.colorbar(cax)
        plt.xlabel("Alongshore Distance (dam)")
        plt.ylabel("Cross-Shore Distance (dam)")
        plt.title("Elevation change (dam)")
        plt.tight_layout()
        timestr = "Time = " + str(t) + " hrs"
        plt.text(1, 1, timestr)
        plt.rcParams.update({"font.size": 15})
        name = "elev_" + str(t)
        elevFig1.savefig(name, facecolor='w')  # dpi=200
        plt.close(elevFig1)

    frames = []

    for filenum in range(start, stop):
        filename = "elev_" + str(filenum) + ".png"
        frames.append(imageio.imread(filename))
    imageio.mimsave("elev.gif", frames, fps=2)
    # imageio.mimsave("elev.gif", frames, "GIF-FI")
    print("[ * elevation GIF successfully generated * ]")


# -------------------------------------------discharge gif--------------------------------------------------------------
def plot_DischargeAnimation(dis, directory, start, stop):
    if not os.path.exists(directory):
        os.makedirs(directory)
    os.chdir(directory)

    newpath = "Discharges/"
    if not os.path.exists(newpath):
        os.makedirs(newpath)
    os.chdir(newpath)

    # for t in range(TMAX - 1):
    for t in range(start, stop):
        AnimateDomain = dis[t]

        # Plot and save
        elevFig1 = plt.figure(figsize=(15, 7))
        ax = elevFig1.add_subplot(111)
        cax = ax.matshow(
            AnimateDomain,
            # origin="upper",
            cmap="jet_r",
            # vmin=0, vmax=20,
        )  # , interpolation='gaussian') # analysis:ignore
        ax.xaxis.set_ticks_position("bottom")
        elevFig1.colorbar(cax)
        plt.xlabel("Alongshore Distance (dam)")
        plt.ylabel("Cross-Shore Distance (dam)")
        plt.title("Discharge (dam^3/hr)")
        plt.tight_layout()
        timestr = "Time = " + str(t) + " hrs"
        plt.text(1, 1, timestr)
        plt.rcParams.update({"font.size": 15})
        name = "dis_" + str(t)
        elevFig1.savefig(name, facecolor='w')  # dpi=200
        plt.close(elevFig1)

    frames = []

    for filenum in range(start, stop):
        filename = "dis_" + str(filenum) + ".png"
        frames.append(imageio.imread(filename))
    imageio.mimsave("dis.gif", frames, fps=2)
    # imageio.mimsave("dis.gif", frames, "GIF-FI")
    print()
    print("[ * discharge GIF successfully generated * ]")


# ---------------------------------------------------slope gif----------------------------------------------------------
def plot_SlopeAnimation(slope, directory, start, stop):
    os.chdir(directory)
    newpath = "Slopes/"
    if not os.path.exists(newpath):
        os.makedirs(newpath)
    os.chdir(newpath)

    for t in range(start, stop):
        AnimateDomain = slope[t]

        # Plot and save
        elevFig1 = plt.figure(figsize=(15, 7))
        ax = elevFig1.add_subplot(111)
        cax = ax.matshow(
            AnimateDomain,
            # origin="upper",
            # cmap="jet_r",
            # vmin=-0.1, vmax=0.1,
        )  # , interpolation='gaussian') # analysis:ignore
        ax.xaxis.set_ticks_position("bottom")
        elevFig1.colorbar(cax)
        plt.xlabel("Alongshore Distance (dam)")
        plt.ylabel("Cross-Shore Distance (dam)")
        plt.title("S2 Slopes")
        plt.tight_layout()
        timestr = "Time = " + str(t) + " hrs"
        plt.text(1, 1, timestr)
        plt.rcParams.update({"font.size": 15})
        name = "slope_" + str(t)
        elevFig1.savefig(name, facecolor='w')  # dpi=200
        plt.close(elevFig1)

    frames = []

    for filenum in range(start, stop):
        filename = "slope_" + str(filenum) + ".png"
        frames.append(imageio.imread(filename))
    imageio.mimsave("dis.gif", frames, fps=2)
    # imageio.mimsave("slopes.gif", frames, "GIF-FI")
    print()
    print("[ * slope GIF successfully generated * ]")


# -------------------------------------------qs2 gif--------------------------------------------------------------------
def plot_Qs2Animation(qs2, directory, start, stop):
    os.chdir(directory)
    newpath = "Qs2/"
    if not os.path.exists(newpath):
        os.makedirs(newpath)
    os.chdir(newpath)

    # for t in range(TMAX - 1):
    for t in range(start, stop):
        AnimateDomain = qs2[t]

        # Plot and save
        elevFig1 = plt.figure(figsize=(15, 7))
        ax = elevFig1.add_subplot(111)
        cax = ax.matshow(
            AnimateDomain,
            # origin="upper",
            cmap="jet_r",
            # vmin=-0.005, vmax=0.05,
        )  # , interpolation='gaussian') # analysis:ignore
        ax.xaxis.set_ticks_position("bottom")
        elevFig1.colorbar(cax)
        plt.xlabel("Alongshore Distance (dam)")
        plt.ylabel("Cross-Shore Distance (dam)")
        plt.title("Qs2 (dam^3)")
        plt.tight_layout()
        timestr = "Time = " + str(t) + " hrs"
        plt.text(1, 1, timestr)
        plt.rcParams.update({"font.size": 15})
        name = "qs2_" + str(t)
        elevFig1.savefig(name, facecolor='w')  # dpi=200
        plt.close(elevFig1)

    frames = []

    for filenum in range(start, stop):
        filename = "qs2_" + str(filenum) + ".png"
        frames.append(imageio.imread(filename))
    imageio.mimsave("dis.gif", frames, fps=2)
    # imageio.mimsave("qs2.gif", frames, "GIF-FI")
    print()
    print("[ * Qs2 GIF successfully generated * ]")


# ---------------------- Sed out array ---------------------------------------------------------------------------------
def plot_SedOutAnimation(sedout, directory, start, stop):
    os.chdir(directory)
    newpath = "SedOut/"
    if not os.path.exists(newpath):
        os.makedirs(newpath)
    os.chdir(newpath)

    # for t in range(TMAX - 1):
    for t in range(start, stop):
        AnimateDomain = sedout[t]

        # Plot and save
        elevFig1 = plt.figure(figsize=(15, 7))
        ax = elevFig1.add_subplot(111)
        cax = ax.matshow(
            AnimateDomain,
            # origin="upper",
            cmap="jet_r",
            # vmin=min_v, vmax=max_v,
        )  # , interpolation='gaussian') # analysis:ignore
        ax.xaxis.set_ticks_position("bottom")
        elevFig1.colorbar(cax)
        plt.xlabel("Alongshore Distance (dam)")
        plt.ylabel("Cross-Shore Distance (dam)")
        plt.title("Sed Flux Out (dam^3)")
        plt.tight_layout()
        timestr = "Time = " + str(t) + " hrs"
        plt.text(1, 1, timestr)
        plt.rcParams.update({"font.size": 15})
        name = "sedout_" + str(t)
        elevFig1.savefig(name, facecolor='w')  # dpi=200
        plt.close(elevFig1)

    frames = []

    for filenum in range(start, stop):
        filename = "sedout_" + str(filenum) + ".png"
        frames.append(imageio.imread(filename))
    imageio.mimsave("dis.gif", frames, fps=2)
    # imageio.mimsave("sedout.gif", frames, "GIF-FI")
    print()
    print("[ * SedOut GIF successfully generated * ]")


# ---------------------- Sed out array ---------------------------------------------------------------------------------
def plot_SedInAnimation(sedin, directory, start, stop):
    os.chdir(directory)
    newpath = "SedIn/"
    if not os.path.exists(newpath):
        os.makedirs(newpath)
    os.chdir(newpath)

    # for t in range(TMAX - 1):
    for t in range(start, stop):
        AnimateDomain = sedin[t]

        # Plot and save
        elevFig1 = plt.figure(figsize=(15, 7))
        ax = elevFig1.add_subplot(111)
        cax = ax.matshow(
            AnimateDomain,
            # origin="upper",
            cmap="jet_r",
            # vmin=min_v, vmax=max_v,
        )  # , interpolation='gaussian') # analysis:ignore
        ax.xaxis.set_ticks_position("bottom")
        elevFig1.colorbar(cax)
        plt.xlabel("Alongshore Distance (dam)")
        plt.ylabel("Cross-Shore Distance (dam)")
        plt.title("Sed Flux In (dam^3)")
        plt.tight_layout()
        timestr = "Time = " + str(t) + " hrs"
        plt.text(1, 1, timestr)
        plt.rcParams.update({"font.size": 15})
        name = "sedin_" + str(t)
        elevFig1.savefig(name, facecolor='w')  # dpi=200
        plt.close(elevFig1)

    frames = []

    for filenum in range(start, stop):
        filename = "sedin_" + str(filenum) + ".png"
        frames.append(imageio.imread(filename))
    imageio.mimsave("dis.gif", frames, fps=2)
    # imageio.mimsave("sedin.gif", frames, "GIF-FI")
    print()
    print("[ * SedIn GIF successfully generated * ]")


def plot_dischargeComp(discharge_array, directory, start, stop, bay_level):
    os.chdir(directory)
    newpath = "dis_comparison/"
    if not os.path.exists(newpath):
        os.makedirs(newpath)
    os.chdir(newpath)

    # for t in range(TMAX - 1):
    for t in range(start, stop):
        # Plot and save
        elevFig1 = plt.figure(figsize=(15, 7))
        ax = elevFig1.add_subplot(111)
        x = range(len(discharge_array[0, 0, :]))
        y = np.ones(len(x))*np.mean(discharge_array[t, 0, :])
        y2 = np.ones(len(x))*np.mean(discharge_array[t, 1, :])
        ax.plot(discharge_array[t, 0, :], label="expected discharge")
        ax.plot(discharge_array[t, 1, :], label="actual discharge")
        ax.plot(x, y, label="average expected discharge", color="k", linestyle="dashed")
        ax.plot(x, y2, label="average actual discharge", linestyle="dashed")
        ax.xaxis.set_ticks_position("bottom")
        plt.xlabel("Alongshore Distance (dam)")
        plt.ylabel("Discharge (dam^3/hr)")
        plt.title("Discharge Comparison at the First Dune Line (dam^3/hr)")
        ax.legend(loc="upper left")
        plt.tight_layout()
        full_text = "Time = " + str(t) + "; Bay level = " + str(round(bay_level[t], 3)) + " dam"
        plt.text(0.5, 0.99, full_text, horizontalalignment='center',
        verticalalignment='top', transform=ax.transAxes)
        plt.rcParams.update({"font.size": 15})
        name = "dis_comparison_" + str(t)
        elevFig1.savefig(name, facecolor='w')  # dpi=200
        plt.close(elevFig1)

    frames = []

    for filenum in range(start, stop):
        filename = "dis_comparison_" + str(filenum) + ".png"
        frames.append(imageio.imread(filename))
    imageio.mimsave("discomp.gif", frames, fps=2)
    # imageio.mimsave("sedin.gif", frames, "GIF-FI")
    print()
    print("[ * Discharge comparison GIF successfully generated * ]")


def plot_FRarray(FR_array, directory, start, stop):
    if not os.path.exists(directory):
        os.makedirs(directory)
    os.chdir(directory)

    newpath = "new_FRA_cells/"
    if not os.path.exists(newpath):
        os.makedirs(newpath)
    os.chdir(newpath)

    # for t in range(TMAX - 1):
    for t in range(start, stop):
        AnimateDomain = FR_array[t]

        # Plot and save
        elevFig1 = plt.figure(figsize=(15, 7))
        ax = elevFig1.add_subplot(111)
        cax = ax.matshow(
            AnimateDomain,
            # origin="upper",
            cmap="binary",
            vmin=0, vmax=1,
        )  # , interpolation='gaussian') # analysis:ignore
        ax.xaxis.set_ticks_position("bottom")
        elevFig1.colorbar(cax)
        plt.xlabel("Alongshore Distance (dam)")
        plt.ylabel("Cross-Shore Distance (dam)")
        plt.title("Cells where flow routing occurs")
        plt.tight_layout()
        timestr = "Time = " + str(t)
        plt.text(1, 1, timestr)
        plt.rcParams.update({"font.size": 15})
        name = "FR_array_" + str(t)
        elevFig1.savefig(name, facecolor='w')  # dpi=200
        plt.close(elevFig1)

    frames = []

    for filenum in range(start, stop):
        filename = "FR_array_" + str(filenum) + ".png"
        frames.append(imageio.imread(filename))
    imageio.mimsave("FR_arrays.gif", frames, fps=2)
    # imageio.mimsave("sedin.gif", frames, "GIF-FI")
    print()
    print("[ * FR array GIF successfully generated * ]")


def plot_FR_slopesarray(FR_slopes_array, directory, start, stop):
    if not os.path.exists(directory):
        os.makedirs(directory)
    os.chdir(directory)

    newpath = "FRA_slope_cells/"
    if not os.path.exists(newpath):
        os.makedirs(newpath)
    os.chdir(newpath)

    # for t in range(TMAX - 1):
    for t in range(start, stop):
        AnimateDomain = FR_slopes_array[t]

        # Plot and save
        elevFig1 = plt.figure(figsize=(15, 7))
        ax = elevFig1.add_subplot(111)
        cax = ax.matshow(
            AnimateDomain,
            # origin="upper",
            cmap="binary",
            # vmin=min_v, vmax=max_v,
        )  # , interpolation='gaussian') # analysis:ignore
        ax.xaxis.set_ticks_position("bottom")
        elevFig1.colorbar(cax)
        plt.xlabel("Alongshore Distance (dam)")
        plt.ylabel("Cross-Shore Distance (dam)")
        plt.title("Cells where flow routing occurs")
        plt.tight_layout()
        timestr = "Time = " + str(t)
        plt.text(1, 1, timestr)
        plt.rcParams.update({"font.size": 15})
        name = "FR_slopes_array_" + str(t)
        elevFig1.savefig(name, facecolor='w')  # dpi=200
        plt.close(elevFig1)

    frames = []

    for filenum in range(start, stop):
        filename = "FR_slopes_array_" + str(filenum) + ".png"
        frames.append(imageio.imread(filename))
    imageio.mimsave("FR_slopes_arrays.gif", frames, fps=2)
    # imageio.mimsave("sedin.gif", frames, "GIF-FI")
    print()
    print("[ * FR slopes array GIF successfully generated * ]")


# -------------------------------------------b3d domain plot------------------------------------------------------------
def plot_ModelTransects(b3d, time_step):
    plt.figure(figsize=(10, 5))
    fig = plt.subplot(1, 1, 1)
    legend_t = []

    for t in time_step:
        # Sea level
        sea_level = b3d._SL + (t * b3d._RSLR[t])

        # Create data points
        shoreface_toe_x = (b3d.x_t_TS[t] - b3d.x_t_TS[0])
        shoreface_toe_y = (sea_level - b3d.DShoreface) * 10  # m
        shoreline_x = (b3d.x_s_TS[t] - b3d.x_t_TS[0])
        shoreline_y = sea_level * 10  # m
        bay_y = (sea_level - b3d._BayDepth) * 10  # m
        end_of_bay_y = bay_y

        # if cascade.nourishments[iB3D].beach_width[t] is not None:
        #     berm_x = shoreline_x + (
        #         cascade.nourishments[iB3D].beach_width[t] / 10
        #     )  # beach width (in dam)
        # else:
        berm_x = shoreline_x + (
            int(b3d.BermEl / b3d._beta)
        )  # initial beach width (in dam)
        # end of un-shifted
        berm_y = (
                         b3d._BermEl * 10
                 ) + shoreline_y  # convert to meters
        dune_toe_x = berm_x
        dune_toe_y = berm_y

        v = 10  # just use 10th transect
        interior_y = b3d._DomainTS[t]  # dam MHW
        interior_y = interior_y[:, v]
        dunes_y = (b3d._DuneDomain[t, v, :] + b3d._BermEl)  # dam MHW
        cross_barrier_y = np.insert(interior_y, 0, dunes_y)
        cross_barrier_y = (cross_barrier_y * 10) + shoreline_y  # Convert to meters with sea level rise included
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
        # plt.rcParams.update({"font.size": 20})
        # legend_t.append(str(t))
        plt.plot(shoreface_toe_x, shoreface_toe_y, 'r-o', label='shoreface toe')
        plt.plot(shoreline_x, shoreline_y, 'g-o', label='shoreline')
        plt.plot(berm_x, berm_y, 'b-o', label='berm')
        plt.plot(dune_toe_x, dune_toe_y, 'c-o', label='dune toe')
        plt.plot(cross_barrier_x, cross_barrier_y, 'm', label='cross barrier')
        plt.plot(end_of_bay_x, end_of_bay_y, 'k-o', label='end of bay')
        plt.legend(loc='lower right')

    plt.xlabel("Cross-shore position (dam)")
    plt.ylabel("Elevation (m MHW)")
    plt.title("Profile Evolution")
    # plt.legend(legend_t)

    return fig

