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
            vmin=-3, vmax=4
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
    imageio.mimsave("elev.gif", frames, fps=1)
    # imageio.mimsave("elev.gif", frames, "GIF-FI")
    print("[ * elevation GIF successfully generated * ]")