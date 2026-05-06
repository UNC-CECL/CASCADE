"""
creates individual plots and a GIF with a constant frame of reference
"""

import numpy as np
import os
import matplotlib.pyplot as plt
import imageio
import math


def plot_domains_dam(casc_model, save_dir, newfolder, buffer):
    """

    :param casc_model: string
        location of the cascade instance (.npz)
    :param save_dir: string
        path to save the files
    :param newfolder: string
        new folder name (figures will go in this folder which will be located in the save_dir)
    :param buffer: int [dam]
        extra space at top (marsh-side) of the cascade plots so the barrier is surrounded by water
    :return:
    """
    for t in range(int(casc_model.TMAX)-1):
    # for t in range(3):

        # determine the necessary width of the figure (maximum of the shoreline position + barrier width at each TS)
        total_width_list = []
        for ts in range(len(casc_model.DomainTS)):
            ts_width = np.shape(casc_model.DomainTS[ts])[0]  # dam
            xs_pos = casc_model.x_s_TS[t] / 10  # dam
            end_width = ts_width + xs_pos
            total_width_list.append(end_width)
        max_width_dam = max(total_width_list)
        plot_width = math.ceil(max_width_dam + buffer)

        # initialize figure using bay depth of -0.3 dam
        domain = np.ones([plot_width, casc_model.BarrierLength]) * -0.3  # [dam]
        start_width = math.floor(casc_model.x_s_TS[t] / 10)  # [dam] rounded down to nearest integer
        domain_width = np.shape(casc_model.DomainTS[t])[0]  # [dam]
        end_width = start_width + domain_width  # [dam]
        casc_model.DomainTS[t][casc_model.DomainTS[t] <=0 ] = -0.3  # set all values < 0 to bay depth
        domain[start_width:end_width, :] = casc_model.DomainTS[t]

        # plotting the domains for comparison
        plt.rcParams.update({"font.size": 11})

        # model 1
        elevFig1 = plt.figure(figsize=(10, 8))
        ax1 = elevFig1.add_subplot(111)
        cax1 = ax1.pcolormesh(
            domain * 10,
            cmap="terrain",
            vmin=-3,
            vmax=6,
            # edgecolors="w",  # for debugging
            # linewidth=0.01,
        )
        cbar1 = elevFig1.colorbar(cax1)
        cbar1.set_label("elevation (m MHW)", rotation=270)
        ax1.set_xlabel("alongshore distance (dam)")
        ax1.set_ylabel("cross-shore distance (dam)")
        # timestr = (
        #         "Time = " + str(t) + " yrs"
        # )
        # ax1.text(0, 0, timestr)
        ax1.set_title("Year {0}".format(t))
        ax1.set_ylim([0, plot_width])

        # Save
        newpath = os.path.join(save_dir, newfolder)
        if not os.path.exists(newpath):
            os.makedirs(newpath)

        name = os.path.join(newpath, "DomainTS_{0}.png".format(t))
        elevFig1.savefig(name)  # dpi=200
        plt.close(elevFig1)


def create_gif(casc_model, folder_path):
    frames = []
    for filenum in range(int(casc_model.TMAX)-1):
        filename = os.path.join(folder_path, "DomainTS_{0}.png".format(filenum))
        frames.append(imageio.imread(filename))
    imageio.mimsave(os.path.join(folder_path, "DomainTS.gif"), frames, "GIF-FI", fps=5)
    print()
    print("[ * GIF successfully generated * ]")