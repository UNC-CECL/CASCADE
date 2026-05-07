"""
creates individual plots and a GIF with a constant frame of reference
"""

import numpy as np
import os
import matplotlib.pyplot as plt
import imageio
import math


def plot_domains(casc_model, iB3D, save_dir, newfolder, buffer, units):
    """

    :param casc_model: string
        location of the cascade instance (.npz)
    :param save_dir: string
        path to save the files
    :param newfolder: string
        new folder name (figures will go in this folder which will be located in the save_dir)
    :param buffer: int [dam]
        extra space on the oceanside of the cascade plots so the barrier is surrounded by water
    :param units: string
        "m" or "km" to change the axes units
    :return:
    """
    for t in range(int(casc_model.barrier3d[iB3D].TMAX)-1):
    # for t in range(1):

        # determine the necessary width of the figure (maximum of the shoreline position + barrier width at each TS)
        total_width_list = []
        for ts in range(len(casc_model.barrier3d[iB3D].DomainTS)):
            ts_width = np.shape(casc_model.barrier3d[iB3D].DomainTS[ts])[0]  # dam
            xs_pos = casc_model.barrier3d[iB3D].x_s_TS[t] / 10  # dam
            end_width = ts_width + xs_pos
            total_width_list.append(end_width)
        max_width_dam = max(total_width_list)
        plot_width = math.ceil(max_width_dam + buffer)

        # initialize figure using bay depth of -0.3 dam
        domain = np.ones([plot_width, casc_model.barrier3d[iB3D].BarrierLength]) * -0.3  # [dam]
        start_width = math.floor(casc_model.barrier3d[iB3D].x_s_TS[t] / 10)  # [dam] rounded down to nearest integer

        # add the barrier dune domain
        dunes = (casc_model.barrier3d[iB3D].DuneDomain[t] + casc_model.barrier3d[iB3D].BermEl)  # [dam]
        dunes = np.rot90(dunes)   # left column (oceanside) is now the bottom row
        dunes = np.flipud(dunes)  # switch rows so oceanside is now on top to match orientation of the interior

        # add the barrier interior domain - oriented high elevation (ocean) on top, marsh on bottom
        int_domain = casc_model.barrier3d[iB3D].DomainTS[t]
        # int_domain[int_domain <= 0] = -0.3  # set all values <= 0 to water
        full_domain = np.vstack([dunes, int_domain])

        # set the barrier domain in the water domain
        domain_width = np.shape(full_domain)[0]  # [dam]
        end_width = start_width + domain_width   # [dam]
        domain[start_width:end_width, :] = full_domain

        # plotting the domains for comparison
        plt.rcParams.update({"font.size": 11})

        # model 1
        elevFig1 = plt.figure(figsize=(10, 8))
        ax1 = elevFig1.add_subplot(111)
        cax1 = ax1.pcolormesh(
            domain * 10,  # [m]
            cmap="terrain",
            vmin=-1,
            vmax=6,
        )
        cbar1 = elevFig1.colorbar(cax1)
        cbar1.set_label("elevation (m MHW)", rotation=270)

        # add text and set title
        ax1.text(0.5, 0.01, 'Ocean', c="white", horizontalalignment='center',
         verticalalignment='bottom', transform=ax1.transAxes)
        ax1.text(0.5, 0.98, 'Marsh/Bay', c="white", horizontalalignment='center',
             verticalalignment='top', transform=ax1.transAxes)
        ax1.set_title("Year {0}".format(t))

        # set axis labels
        if units == "km":
            xticks = ax1.get_xticks()
            xlabels = xticks / 100
            ax1.set_xticklabels(xlabels)
            yticks = ax1.get_yticks()
            ylabels = yticks / 100
            ax1.set_yticklabels(ylabels)
            ax1.set_xlabel("alongshore distance (km)")
            ax1.set_ylabel("cross-shore distance (km)")
        elif units == "m":
            xticks = ax1.get_xticks()
            xlabels = xticks * 10
            ax1.set_xticklabels(xlabels.astype(int))
            yticks = ax1.get_yticks()
            ylabels = yticks * 10
            ax1.set_yticklabels(ylabels.astype(int))
            ax1.set_xlabel("alongshore distance (m)")
            ax1.set_ylabel("cross-shore distance (m)")
        else:
            ax1.set_xlabel("alongshore distance (dam)")
            ax1.set_ylabel("cross-shore distance (dam)")


        # Save
        newpath = os.path.join(save_dir, newfolder)
        if not os.path.exists(newpath):
            os.makedirs(newpath)

        name = os.path.join(newpath, "DomainTS_{0}.png".format(t))
        elevFig1.savefig(name)  # dpi=200
        plt.close(elevFig1)


def create_gif(casc_model, iB3D, folder_path):
    frames = []
    for filenum in range(int(casc_model.barrier3d[iB3D].TMAX)-1):
        filename = os.path.join(folder_path, "DomainTS_{0}.png".format(filenum))
        frames.append(imageio.imread(filename))
    imageio.mimsave(os.path.join(folder_path, "DomainTS.gif"), frames, "GIF-FI", fps=5)
    print()
    print("[ * GIF successfully generated * ]")

def plot_domains_with_marsh_transect(casc_model, iB3D, save_dir, newfolder, buffer, units):
    """

    :param casc_model: string
        location of the cascade instance (.npz)
    :param save_dir: string
        path to save the files
    :param newfolder: string
        new folder name (figures will go in this folder which will be located in the save_dir)
    :param buffer: int [dam]
        extra space on the oceanside of the cascade plots so the barrier is surrounded by water
    :param units: string
        "m" or "km" to change the axes units
    :return:
    """
    for t in range(int(casc_model.barrier3d[iB3D].TMAX)-1):
    # for t in range(20):

        # determine the necessary width of the figure (maximum of the shoreline position + barrier width at each TS)
        total_width_list = []
        for ts in range(len(casc_model.barrier3d[iB3D].DomainTS)):
            ts_width = np.shape(casc_model.barrier3d[iB3D].DomainTS[ts])[0]  # dam
            xs_pos = casc_model.barrier3d[iB3D].x_s_TS[t] / 10  # dam
            end_width = ts_width + xs_pos
            total_width_list.append(end_width)
        max_width_dam = max(total_width_list)
        plot_width = math.ceil(max_width_dam + buffer)

        # initialize figure using bay depth of -0.3 dam
        domain = np.ones([plot_width, casc_model.barrier3d[iB3D].BarrierLength]) * -0.3  # [dam]
        start_width = math.floor(casc_model.barrier3d[iB3D].x_s_TS[t] / 10)  # [dam] rounded down to nearest integer

        # add the barrier dune domain
        dunes = (casc_model.barrier3d[iB3D].DuneDomain[t] + casc_model.barrier3d[iB3D].BermEl)  # [dam]
        dunes = np.rot90(dunes)   # left column (oceanside) is now the bottom row
        dunes = np.flipud(dunes)  # switch rows so oceanside is now on top to match orientation of the interior

        # add the barrier interior domain - oriented high elevation (ocean) on top, marsh on bottom
        int_domain = casc_model.barrier3d[iB3D].DomainTS[t]
        # int_domain[int_domain <= 0] = -0.3  # set all values <= 0 to water
        full_domain = np.vstack([dunes, int_domain])

        # set the barrier domain in the water domain
        domain_width = np.shape(full_domain)[0]  # [dam]
        end_width = start_width + domain_width   # [dam]
        domain[start_width:end_width, :] = full_domain

        # add the transect
        BB_transect_m1 = np.flip(casc_model.barrierbmft[iB3D].bmftc_BB.elevation[casc_model.barrierbmft[iB3D].bmftc_BB.startyear + t - 1, int(casc_model.barrierbmft[iB3D].bmftc_BB.Marsh_edge[casc_model.barrierbmft[iB3D].bmftc_BB.startyear + t]):])
        if casc_model.barrierbmft[iB3D].x_b_TS_ML[t] < 0:
            ML_transect_m1 = np.append(np.ones([abs(int(casc_model.barrierbmft[iB3D].x_b_TS_ML[t]))]) * casc_model.barrierbmft[iB3D].bmftc_ML.elevation[casc_model.barrierbmft[iB3D].bmftc_ML.startyear + t - 1, 1], casc_model.barrierbmft[iB3D].bmftc_ML.elevation[casc_model.barrierbmft[iB3D].bmftc_ML.startyear + t - 1, :])
        elif casc_model.barrierbmft[iB3D].x_b_TS_ML[t] > 0:
            ML_transect_m1 = casc_model.barrierbmft[iB3D].bmftc_ML.elevation[casc_model.barrierbmft[iB3D].bmftc_ML.startyear + t - 1, int(casc_model.barrierbmft[iB3D].x_b_TS_ML[t]):]

        whole_transect_m1 = np.append(BB_transect_m1, ML_transect_m1)

        # Extract marsh & forest locations to plot
        x_forest_m1 = [casc_model.barrierbmft[iB3D].bmftc_BB.B - int(casc_model.barrierbmft[iB3D].bmftc_BB.Forest_edge[casc_model.barrierbmft[iB3D].bmftc_BB.startyear + t]), len(BB_transect_m1) - int(casc_model.barrierbmft[iB3D].x_b_TS_ML[t]) + int(casc_model.barrierbmft[iB3D].bmftc_ML.Forest_edge[casc_model.barrierbmft[iB3D].bmftc_ML.startyear + t])]
        x_marsh_m1 = [len(BB_transect_m1) - 1, len(BB_transect_m1) - int(casc_model.barrierbmft[iB3D].x_b_TS_ML[t]) + int(casc_model.barrierbmft[iB3D].bmftc_ML.Marsh_edge[casc_model.barrierbmft[iB3D].bmftc_ML.startyear + t])]
        y_forest_m1 = whole_transect_m1[x_forest_m1]
        y_marsh_m1 = whole_transect_m1[x_marsh_m1]

        # plot
        plt.rcParams.update({"font.size": 11})
        elevFig1 = plt.figure(figsize=(8, 8))

        # domains on subplot 1 -----------------------------------
        nrows = 3
        ncols = 9  # at 5, marsh figure right boundary aligns with domain right boundary (before colorbar)
        ax1 = plt.subplot2grid(shape=(nrows, ncols), loc=(0, 0), rowspan=2, colspan=ncols)
        # shape = (#rows, #cols), loc=(row index, col index), rowspan = # rows to cover
        # we want the domain to take up 2/3 of the figure space
        cax1 = ax1.pcolormesh(
            domain * 10,  # [m]
            cmap="terrain",
            vmin=-1,
            vmax=6,
        )
        cbar1 = elevFig1.colorbar(cax1)
        cbar1.set_label("elevation (m MHW)", rotation=270)

        # add text and set title
        ax1.text(0.5, 0.01, 'Ocean', c="white", horizontalalignment='center',
         verticalalignment='bottom', transform=ax1.transAxes)
        ax1.text(0.5, 0.98, 'Marsh/Bay', c="white", horizontalalignment='center',
             verticalalignment='top', transform=ax1.transAxes)
        ax1.set_title("Year {0}".format(t))

        # set axis labels
        if units == "km":
            xticks = ax1.get_xticks()
            xlabels = xticks / 100
            ax1.set_xticklabels(xlabels)
            yticks = ax1.get_yticks()
            ylabels = yticks / 100
            ax1.set_yticklabels(ylabels)
            ax1.set_xlabel("alongshore distance (km)")
            ax1.set_ylabel("cross-shore distance (km)")
        elif units == "m":
            xticks = ax1.get_xticks()
            xlabels = xticks * 10
            ax1.set_xticklabels(xlabels.astype(int))
            yticks = ax1.get_yticks()
            ylabels = yticks * 10
            ax1.set_yticklabels(ylabels.astype(int))
            ax1.set_xlabel("alongshore distance (m)")
            ax1.set_ylabel("cross-shore distance (m)")
        else:
            ax1.set_xlabel("alongshore distance (dam)")
            ax1.set_ylabel("cross-shore distance (dam)")

        # transects on subplot 2 -----------------------------------
        ax2 = plt.subplot2grid(shape=(nrows, ncols), loc=(nrows-1, 0), rowspan=1, colspan=ncols-1)
        # shape = (#rows, #cols), loc=(row index, col index), rowspan = # rows to cover
        # we want the transect to take up 1/3 of the figure space
        ax2.plot(whole_transect_m1, c="black")
        ax2.scatter(x_forest_m1, y_forest_m1, c="green")
        ax2.scatter(x_marsh_m1, y_marsh_m1, c="brown")
        ax2.set_ylim(-2.1, 4)
        ax2.set_xlabel("Cross-shore Distance (m)")
        ax2.set_ylabel("Elevation (m)")

        # add text and set title
        ax2.text(0.01, 0.05, 'Ocean', c="grey", horizontalalignment='left',
         verticalalignment='center', transform=ax2.transAxes)
        ax2.text(0.99, 0.05, 'Mainland', c="grey", horizontalalignment='right',
             verticalalignment='center', transform=ax2.transAxes)

        plt.tight_layout()

        # Save
        newpath = os.path.join(save_dir, newfolder)
        if not os.path.exists(newpath):
            os.makedirs(newpath)

        name = os.path.join(newpath, "DomainTS_{0}.png".format(t))
        elevFig1.savefig(name)  # dpi=200
        plt.close(elevFig1)
