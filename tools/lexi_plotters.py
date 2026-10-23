# plotting functions for
#
# ~******* CASCADE ********~
#
"""----------------------------------------------------
Copyright (C) 2022 Katherine Anarde
----------------------------------------------------"""

import math
import os
from matplotlib.lines import Line2D

import matplotlib.pyplot as plt
import numpy as np
from matplotlib import colormaps


# new from Lexi
def plot_start_end_domains(
        cascade_b3d,
        save_dir,
        min_z,
        max_z,
        figsize=[10,8]
):
    """
    plots the start cascade domain (model year 1) with the end domain (dunes and interior)
    and saves them to a specified directory
    :param max_z: min elevation (m MHW) for plot
    :param min_z: max elevation (m MHW) for plot
    :param cascade_b3d: barrier3d object of a cascade class
    :param save_dir: folder path to save the plots. makes it if it does not exist

    :return: nothing, saves the plots to specified location
    """
    # plot parameters
    plt.rcParams["font.size"] = 12
    plt.ioff()

    # create the save folder if it does not already exist
    if not os.path.exists(save_dir):
        os.makedirs(save_dir)

    # loop through each domain and plot model years 1 and the last
    ib3ds = len(cascade_b3d)

    for i in range(ib3ds):
        # make 1 plot per domain and save it after
        fig1 = plt.figure(figsize=figsize)
        # create and format the plot
        xlabel = "alongshore distance (dam)"
        ylabel = "cross-shore distance (dam)"
        minz = min_z
        maxz = max_z

        # interiors in decameters
        first_model_year = cascade_b3d[i].DomainTS[1]  # model year 1
        last_model_year = cascade_b3d[i].DomainTS[-1]  # last model year

        # dunes in decameters above berm elevation
        first_model_year_dunes = cascade_b3d[i].DuneDomain[1]  # model year 1
        last_model_year_dunes = cascade_b3d[i].DuneDomain[-1]  # last model year

        # berm elev
        berm_elev = cascade_b3d[i].BermEl  # dam MHW

        # TMAX
        tmax = cascade_b3d[i].TMAX

        # convert domains into m MHW
        first_model_year = first_model_year * 10
        last_model_year = last_model_year * 10
        first_model_year_dunes = (first_model_year_dunes + berm_elev) * 10
        last_model_year_dunes = (last_model_year_dunes + berm_elev) * 10

        # flip dunes
        first_model_year_dunes = np.transpose(first_model_year_dunes)
        last_model_year_dunes = np.transpose(last_model_year_dunes)

        # ---------------- plotting the domains --------------------------------
        n_plots = 2

        # model year 1 domain
        ax1 = fig1.add_subplot(1, n_plots, 1)  # rows, cols, plot number
        mat1 = ax1.matshow(
            np.vstack((first_model_year_dunes, first_model_year)),
            cmap="terrain",
            vmin=minz,
            vmax=maxz,
        )
        cbar = fig1.colorbar(mat1)
        cbar.set_label('m MHW', rotation=270, labelpad=10)
        ax1.set_title("model year 1")
        ax1.set_ylabel(ylabel)
        ax1.set_xlabel(xlabel)
        plt.gca().xaxis.tick_bottom()

        # last model year domain
        ax2 = fig1.add_subplot(1, n_plots, 2)  # rows, cols, plot number
        mat1 = ax2.matshow(
            np.vstack((last_model_year_dunes, last_model_year)),
            cmap="terrain",
            vmin=minz,
            vmax=maxz,
        )
        cbar = fig1.colorbar(mat1)
        cbar.set_label('m MHW', rotation=270, labelpad=10)
        ax2.set_title("model year {0}".format(tmax))
        ax2.set_ylabel(ylabel)
        ax2.set_xlabel(xlabel)
        plt.gca().xaxis.tick_bottom()

        # once all years are plotted, save the figure
        fig1.suptitle("iB3D {0}".format(i))
        fig1.tight_layout()
        fig1.savefig(os.path.join(save_dir, "iB3D_{0}.png".format(i)))
        plt.close(fig1)

    print("finished. saved {0} figures to {1}".format(ib3ds, save_dir))


def plot_overwash_flux(
        cascade_b3d,
        save_dir,
        figsize=[10,8]
):
    """
    plots overwash along all domains on one figure
    :param cascade_b3d:
    :param save_dir:
    :param figsize:
    :return:
    """
    # plot parameters
    plt.rcParams["font.size"] = 12
    plt.ioff()

    # create the save folder if it does not already exist
    if not os.path.exists(save_dir):
        os.makedirs(save_dir)

    # tmax is same for all
    tmax = cascade_b3d[0].TMAX
    blues_cmap = colormaps['Blues']
    blues = blues_cmap(np.linspace(0.2, 1, tmax))  # n evenly spaced colors
    markers = list(Line2D.markers.keys())

    # loop through each ib3d and add the overwash flux to an array
    ib3ds = len(cascade_b3d)
    all_overwash = np.zeros([ib3ds,tmax])  # nrows = ib3ds, ncols = nyears
    for i in range(ib3ds):
        overwash = cascade_b3d[i].QowTS  # m^3/m, average for the domain for all model years
        all_overwash[i] = overwash

    # plot each column, which corresponds to a year and has one value per domain
    fig1 = plt.figure(figsize=figsize)
    ax1 = fig1.add_subplot(111)  # rows, cols, plot number
    xlabel = "iB3D domain"
    ylabel = "Overwash Flux (m^3/m)"
    for y in range(tmax):
        ax1.plot(all_overwash[:,y], c=blues[y], marker=markers[y], label="year {0}".format(y))
    # other plot features
    ax1.legend(loc="upper right")
    ax1.set_ylabel(ylabel)
    ax1.set_xlabel(xlabel)
    ax1.set_xticks(np.arange(0,ib3ds,1))
    ax1.set_ylim(0,100)
    # save the figure
    fig1.tight_layout()
    fig1.savefig(os.path.join(save_dir, "overwash.png"))
    plt.close(fig1)

    return all_overwash


def compare_modeled_observed_domains(
        cascade_b3d,
        save_dir,
        observed_domains,
        observed_dunes_list,
        min_z,
        max_z,
        figsize=[15, 8]
):
    # plot parameters
    plt.rcParams["font.size"] = 12
    plt.ioff()

    # create the save folder if it does not already exist
    if not os.path.exists(save_dir):
        os.makedirs(save_dir)

    # loop through each domain and plot the first and last model year, plus observed year
    ib3ds = len(cascade_b3d)

    for i in range(ib3ds):
        # make 1 plot per domain and save it after
        fig1 = plt.figure(figsize=figsize)
        # create and format the plot
        xlabel = "alongshore distance (dam)"
        ylabel = "cross-shore distance (dam)"
        minz = min_z
        maxz = max_z

        # interiors in decameters
        first_model_year = cascade_b3d[i].DomainTS[0]  # initial domain
        last_model_year = cascade_b3d[i].DomainTS[-1]  # last model year
        observed = np.load(observed_domains[i]) # observed

        # dunes in decameters above berm elevation
        first_model_year_dunes = cascade_b3d[i].DuneDomain[0]  # initial domain
        last_model_year_dunes = cascade_b3d[i].DuneDomain[-1]  # last model year
        observed_dunes = np.load(observed_dunes_list[i])  # observed

        # berm elev
        berm_elev = cascade_b3d[i].BermEl  # dam MHW

        # TMAX
        tmax = cascade_b3d[i].TMAX

        # convert domains into m MHW
        first_model_year = first_model_year * 10
        last_model_year = last_model_year * 10
        observed = observed * 10
        first_model_year_dunes = (first_model_year_dunes + berm_elev) * 10
        last_model_year_dunes = (last_model_year_dunes + berm_elev) * 10
        observed_dunes = (observed_dunes + berm_elev) * 10

        # flip dunes
        first_model_year_dunes = np.transpose(first_model_year_dunes)
        last_model_year_dunes = np.transpose(last_model_year_dunes)
        observed_dunes = np.transpose(observed_dunes)

        # ---------------- plotting the domains --------------------------------
        n_plots = 3

        # model year 1 domain
        ax1 = fig1.add_subplot(1, n_plots, 1)  # rows, cols, plot number
        mat1 = ax1.matshow(
            np.vstack((first_model_year_dunes, first_model_year)),
            cmap="terrain",
            vmin=minz,
            vmax=maxz,
        )
        cbar = fig1.colorbar(mat1)
        cbar.set_label('m MHW', rotation=270, labelpad=10)
        ax1.set_title("model year 0")
        ax1.set_ylabel(ylabel)
        ax1.set_xlabel(xlabel)
        plt.gca().xaxis.tick_bottom()

        # last model year domain
        ax2 = fig1.add_subplot(1, n_plots, 2)  # rows, cols, plot number
        mat1 = ax2.matshow(
            np.vstack((last_model_year_dunes, last_model_year)),
            cmap="terrain",
            vmin=minz,
            vmax=maxz,
        )
        cbar = fig1.colorbar(mat1)
        cbar.set_label('m MHW', rotation=270, labelpad=10)
        ax2.set_title("model year {0}".format(tmax-1))
        ax2.set_ylabel(ylabel)
        ax2.set_xlabel(xlabel)
        plt.gca().xaxis.tick_bottom()

        # observed results
        ax3 = fig1.add_subplot(1, n_plots, 3)  # rows, cols, plot number
        mat1 = ax3.matshow(
            np.vstack((observed_dunes, observed)),
            cmap="terrain",
            vmin=minz,
            vmax=maxz,
        )
        cbar = fig1.colorbar(mat1)
        cbar.set_label('m MHW', rotation=270, labelpad=10)
        ax3.set_title("observed".format(tmax))
        ax3.set_ylabel(ylabel)
        ax3.set_xlabel(xlabel)
        plt.gca().xaxis.tick_bottom()

        # once all years are plotted, save the figure
        fig1.suptitle("iB3D {0}".format(i))
        fig1.tight_layout()
        fig1.savefig(os.path.join(save_dir, "iB3D_{0}.png".format(i)))
        plt.close(fig1)

    print("finished. saved {0} figures to {1}".format(ib3ds, save_dir))

def plot_linear_dune_dif(
        cascade_b3d,
        save_dir,
        observed_dunes_list,
        min_z,
        max_z,
        figsize=[15, 8]
):
    # plot parameters
    plt.rcParams["font.size"] = 12
    plt.ioff()

    # create the save folder if it does not already exist
    if not os.path.exists(save_dir):
        os.makedirs(save_dir)

    # loop through each domain and plot the first and last model year, plus observed year
    ib3ds = len(cascade_b3d)


    # intialize the arrays for tracking over the entire island
    modeled_dif_island = []
    observed_dif_island = []

    # make individual domain plots and one showing all domains
    for i in range(ib3ds):
        # make 1 plot per domain and save it after
        fig1 = plt.figure(figsize=figsize)
        # create and format the plot
        xlabel = "alongshore distance (dam)"
        ylabel = "change in dune height (m)"
        minz = min_z
        maxz = max_z

        # dunes in decameters above berm elevation
        first_model_year_dunes = cascade_b3d[i].DuneDomain[0]  # initial conditions
        last_model_year_dunes = cascade_b3d[i].DuneDomain[-1]  # last model year
        observed_dunes = np.load(observed_dunes_list[i])  # observed

        # convert into m
        first_model_year_dunes = first_model_year_dunes * 10
        last_model_year_dunes = last_model_year_dunes * 10
        observed_dunes = observed_dunes * 10

        # only use the crests
        first_model_year_dunes = np.max(first_model_year_dunes, axis=1)
        last_model_year_dunes = np.max(last_model_year_dunes, axis=1)
        if len(np.shape(observed_dunes)) == 2:
            observed_dunes = np.max(observed_dunes, axis=1)

        # differences
        modeled_dif = last_model_year_dunes - first_model_year_dunes
        observed_dif = observed_dunes - first_model_year_dunes
        modeled_dif_island.append(modeled_dif)
        observed_dif_island.append(observed_dif)


        # ---------------- plotting the dunes --------------------------------

        # model year 1 domain
        ax1 = fig1.add_subplot(1, 1, 1)  # rows, cols, plot number
        ax1.plot(modeled_dif, label="modeled")
        ax1.plot(observed_dif, label="observed")
        ax1.axhline(0, c="black", ls="dashed")
        ax1.set_ylabel(ylabel)
        ax1.set_xlabel(xlabel)
        ax1.set_ylim(minz, maxz)
        ax1.legend(loc="upper right")


        # once all years are plotted, save the figure
        fig1.suptitle("iB3D {0}".format(i))
        fig1.tight_layout()
        fig1.savefig(os.path.join(save_dir, "iB3D_{0}.png".format(i)))
        plt.close(fig1)


    # plot the full changes
    modeled_dif_island = np.concatenate(modeled_dif_island)
    observed_dif_island = np.concatenate(observed_dif_island)
    fig2 = plt.figure(figsize=figsize)
    ax1 = fig2.add_subplot(1, 1, 1)
    ax1.plot(modeled_dif_island, label="modeled")
    ax1.plot(observed_dif_island, label="observed")
    ax1.axhline(0, c="black", ls="dashed")
    ax1.set_ylabel(ylabel)
    ax1.set_xlabel("alongshore distance, north to south (dam)")
    ax1.set_ylim(minz, maxz)
    ax1.legend(loc="upper right")
    # plot domain edges as vertical lines
    domain_len_dam = len(modeled_dif)
    domain_edges = np.arange(0,len(modeled_dif_island)+domain_len_dam,domain_len_dam)
    plt.vlines(domain_edges, ymin=minz, ymax=maxz, colors="black", ls="dotted", linewidth=0.5)
    # plot domain labels as text at the top of the figure, halfway between the domain limits
    start_text = domain_len_dam / 2  # first label at 25 dam
    end_text = (ib3ds+1)*domain_len_dam - start_text
    text_pos = np.arange(start_text, end_text, domain_len_dam)
    d_text = 0  # ib3d domain
    for t in range(len(text_pos)):
        plt.text(x=text_pos[t], y=minz+0.5, s="{0}".format(d_text), va="center", ha="center")
        d_text += 1

    # once all years are plotted, save the figure
    fig2.suptitle("difference in dune heights")
    fig2.tight_layout()
    fig2.savefig(os.path.join(save_dir, "dune_height_dif_island.png"))
    plt.close(fig2)

    print("finished. saved {0} figures to {1}".format(ib3ds, save_dir))

    # # save numpy arrays for comparison
    # np.save(os.path.join(save_dir, "dune_height_dif_island"), observed_dif_island)


def plot_shoreline_all_domains(cascade, years, save_dir, figsize=[10,4], seaward_pos=True):
    """
    :param cascade: cascade class loaded from the output npz file
    :param years: list of the year(s) you want to plot
    :param figsize: size of the figure
    :param seaward_pos: if True, an increase in shoreline position
        corresponds to seaward movement of the barrier
    :return: fig

    Shoreline positions for a model year over each domain

    """
    # create the save folder if it does not already exist
    if not os.path.exists(save_dir):
        os.makedirs(save_dir)

    plt.rcParams["font.size"] = 16
    fig1 = plt.figure(figsize=figsize)
    for year in years:
        shorelines = []
        # group all domains from the same year in one list
        for iB3D in range(cascade._ny):
            b3d = cascade.barrier3d[iB3D]
            if seaward_pos:
                x_s = b3d.x_s_TS[year] * -10  # in b3d, shorelines are in dam
                ylab = "Shoreline Position (m) \n <--- landward || seaward --->"
            else:
                x_s = b3d.x_s_TS[year] * 10  # in b3d, shorelines are in dam
                ylab = "Shoreline Position (m) \n <--- seaward || landward --->"
            shorelines.append(x_s)  # all shorelines per domain
        plt.plot(shorelines, marker='o', markersize=4, label="year {0}".format(year))
    plt.ylabel(ylab)
    plt.xlabel("iB3D Domain")
    plt.legend()
    fig1.tight_layout()
    fig1.savefig(os.path.join(save_dir, "shorelines_all_domains.png"))
    plt.close()

def plot_shoreline_change_all_domains(cascade, years, coastsat, save_dir, title, figsize=[10,4], seaward_pos=True):
    """
    :param cascade: cascade class loaded from the output npz file
    :param years: list of two model years to use for calculating the difference
    :param coastsat: list of observed change (in same orientation as iB3D domains)
    :return: fig

    Shoreline positions for a model year over each domain

    """
    # create the save folder if it does not already exist
    if not os.path.exists(save_dir):
        os.makedirs(save_dir)

    plt.rcParams["font.size"] = 16
    fig1 = plt.figure(figsize=figsize)
    shoreline_change = []
    year1 = years[0]
    year2 = years[1]
    # group all domains in one list
    for iB3D in range(cascade._ny):
        b3d = cascade.barrier3d[iB3D]
        xs1 = b3d.x_s_TS[year1] * 10  # in b3d, shorelines are in dam
        xs2 = b3d.x_s_TS[year2] * 10
        # in b3d, shorelines are positive as they move landward, but coastsat uses the opposite sign
        if seaward_pos:
            dif = xs1 - xs2
        else:
            dif = xs2 - xs1
        shoreline_change.append(dif)  # all shorelines per domain
    # plot the modeled and observed data
    plt.plot(shoreline_change, marker='o', markersize=4, label="modeled")
    plt.plot(coastsat, marker='o', markersize=4, label="observed")
    plt.ylabel("Change in Shoreline Position (m) \n <--- landward || seaward --->")
    plt.xlabel("iB3D Domain")
    plt.axhline(y=0,color="grey", linewidth=0.5)
    plt.legend()
    plt.title(title)
    fig1.tight_layout()
    fig1.savefig(os.path.join(save_dir, "shoreline_change_comparison.png"))
    plt.close()


def plot_avg_increase_marsh_elev(
        cascade,
        comparison_year=None,
        figsize=[15,6]
):
    """
    calculates the average change in elevation between "marsh" cells
    for comparison with an expected change over the time frame
    plots the average for each ib3d domain
    :param cascade: cascade class loaded from the output npz file
    :param comparison_year: specific model year for comparison to initial
        if None, use the last time step
    :param figsize: size of the output figure
    :return: fig
    """
    plt.rcParams["font.size"] = 16

    # set barrier3d and marsh classes for easier calling
    ib3ds = cascade._ny
    b3d = cascade.barrier3d
    marsh = cascade.marsh
    avg_elev_change_ib3d = np.ones(ib3ds)*np.nan

    # loop through each domain and find the average difference in the marsh elevation
    for i in range(ib3ds):
        initial_elev = b3d[i].DomainTS[0]
        if comparison_year is None:
            final_elev = b3d[i].DomainTS[-1]
        else:
            final_elev = b3d[i].DomainTS[comparison_year]

        # initial cells that were identified as marsh cells
        min_marsh_elev = marsh[i]._m_min  # dam MHW
        max_marsh_elev = marsh[i]._m_max  # dam MHW
        marsh_cells = np.where((initial_elev > min_marsh_elev) & (initial_elev < max_marsh_elev))
        # if there were no marsh cells to begin with, just keep the array value
        if len(marsh_cells[0])>0:
            # compare the elevations of the initial marsh cells to their final and take average
            initial_marsh_elev = initial_elev[marsh_cells]
            final_marsh_elev = final_elev[marsh_cells]
            avg_elev_change = np.mean(final_marsh_elev - initial_marsh_elev)  # positive means they grew
            # convert to m and add to array
            avg_elev_change_m = avg_elev_change * 10
            avg_elev_change_ib3d[i] = avg_elev_change_m

    # plot the avg per domain and the avg of all domains
    overall_avg = round(np.nanmean(avg_elev_change_ib3d),3)
    domains = np.arange(0,ib3ds,1)
    fig1 = plt.figure(figsize=figsize)
    plt.plot(avg_elev_change_ib3d, marker="o")
    plt.axhline(y=overall_avg, label="overall avg = {0}".format(overall_avg), color="grey", linewidth=0.5)
    plt.xticks(domains, domains)
    plt.ylim(-0.2,0.2)
    plt.xlabel("iB3D Domain")
    plt.ylabel("avg change in elevation (m)")
    plt.title("avg change in elevation of cells initially defined as 'marsh'")
    plt.legend(loc="lower right")
    plt.tight_layout()

