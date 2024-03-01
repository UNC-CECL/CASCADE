# plotting the time series of dune crest elevations
# adapted from Cascade Plotters

# Lexi Van Blunk
# 2/23/2024

from matplotlib import pyplot as plt
import numpy as np


rname_array = ["r025", "r035"]

for rname in rname_array:
    storm_interval = 20        # 20 or 10 years
    config = 4                 # 1, 2, 3, or 4
    storm_num = 1

    # location of the npz files
    datadir_b3d = "C:/Users/Lexi/PycharmProjects/CASCADE/cascade/data/outwash_data/storms/slope0pt03/run_output/{0}/overwash_only/".format(rname)
    datadir_100 = "C:/Users/Lexi/PycharmProjects/CASCADE/cascade/data/outwash_data/storms/slope0pt03/run_output/{0}/outwash100/".format(rname)
    datadir_50 = "C:/Users/Lexi/PycharmProjects/CASCADE/cascade/data/outwash_data/storms/slope0pt03/run_output/{0}/outwash50/".format(rname)
    datadir_0 = "C:/Users/Lexi/PycharmProjects/CASCADE/cascade/data/outwash_data/storms/slope0pt03/run_output/{0}/outwash0/".format(rname)


    # Barrier3d only
    filename_b3d = "config{0}_b3d_startyr1_interval{1}yrs_Slope0pt03_{2}.npz".format(config, storm_interval, storm_num)
    file_b3d = datadir_b3d + filename_b3d
    b3d_numpy_obj = np.load(file_b3d, allow_pickle=True)
    b3d_obj = b3d_numpy_obj["cascade"][0]
    b3d = b3d_obj.barrier3d

    # 100% variables
    filename_100 = "config{0}_outwash100_startyr1_interval{1}yrs_Slope0pt03_{2}.npz".format(config, storm_interval, storm_num)
    file_100 = datadir_100 + filename_100
    outwash100_numyp_obj = np.load(file_100, allow_pickle=True)
    outwash100_obj = outwash100_numyp_obj["cascade"][0]
    outwash100 = outwash100_obj.barrier3d

    # 50% variables
    filename_50 = "config{0}_outwash50_startyr1_interval{1}yrs_Slope0pt03_{2}.npz".format(config, storm_interval, storm_num)
    file_50 = datadir_50 + filename_50
    outwash50_numyp_obj = np.load(file_50, allow_pickle=True)
    outwash50_obj = outwash50_numyp_obj["cascade"][0]
    outwash50 = outwash50_obj.barrier3d

    # 0% variables
    filename_0 = "config{0}_outwash0_startyr1_interval{1}yrs_Slope0pt03_{2}.npz".format(config, storm_interval, storm_num)
    file_0 = datadir_0 + filename_0
    outwash0_numpy_obj = np.load(file_0, allow_pickle=True)
    outwash0_obj = outwash0_numpy_obj["cascade"][0]
    outwash0 = outwash0_obj.barrier3d

    ### ------------------------------------------- Dunes ------------------------------------------------------------------
    TMAX = 101
    vmin = 0
    vmax = 6
    fontsize = 12

    # Barrier3d only
    DuneCrest = []

    for iB3D in range(len(b3d)):
        sub_domain = b3d[iB3D]._DuneDomain[0:TMAX, :, :]
        DuneCrest.append(sub_domain.max(axis=2))

    DuneCrest = np.hstack(DuneCrest).astype(float)

    duneFig = plt.figure(figsize=(20, 14))
    plt.rcParams.update({"font.size": fontsize})
    ax = duneFig.add_subplot(221)
    cax = ax.matshow(
        np.flip(DuneCrest,1) * 10,
        origin="lower",
        cmap="bwr",
        aspect="auto",
        vmin=vmin,
        vmax=vmax,
    )
    ax.xaxis.set_ticks_position("bottom")  # analysis:ignore
    # cax = ax.xaxis.set_ticks_position("bottom")  # analysis:ignore
    cbar = duneFig.colorbar(cax)
    cbar.set_label('Dune Height Above Berm Elevation (m)', rotation=270, labelpad=25)
    plt.xlabel("Alongshore Distance (m)")
    plt.ylabel("Year")
    plt.title("overwash only", weight="bold")
    # plt.hlines(20, -0.5, 49.5, color="k", linestyles='dashed', linewidth=1)

    xtick_max = np.shape(DuneCrest)[1]  # n_cols = x
    x_ticks = np.array(range(0, xtick_max, 10))
    x_tick_labels = x_ticks * 10
    plt.xticks(x_ticks, x_tick_labels)

    # cascade 100%
    DuneCrest = []

    for iB3D in range(len(b3d)):
        sub_domain = outwash100[iB3D]._DuneDomain[0:TMAX, :, :]
        DuneCrest.append(sub_domain.max(axis=2))

    DuneCrest = np.hstack(DuneCrest).astype(float)

    plt.rcParams.update({"font.size": fontsize})
    ax = duneFig.add_subplot(224)
    cax = ax.matshow(
        np.flip(DuneCrest,1) * 10,
        origin="lower",
        cmap="bwr",
        aspect="auto",
        vmin=vmin,
        vmax=vmax,
    )
    ax.xaxis.set_ticks_position("bottom")  # analysis:ignore
    # cax = ax.xaxis.set_ticks_position("bottom")  # analysis:ignore
    cbar = duneFig.colorbar(cax)
    cbar.set_label('Dune Height Above Berm Elevation (m)', rotation=270, labelpad=25)
    plt.xlabel("Alongshore Distance (m)")
    plt.ylabel("Year")
    plt.title("100% washout to shoreface", weight="bold")
    plt.xticks(x_ticks, x_tick_labels)

    # cascade 50%
    DuneCrest = []

    for iB3D in range(len(b3d)):
        sub_domain = outwash50[iB3D]._DuneDomain[0:TMAX, :, :]
        DuneCrest.append(sub_domain.max(axis=2))

    DuneCrest = np.hstack(DuneCrest).astype(float)

    plt.rcParams.update({"font.size": fontsize})
    ax = duneFig.add_subplot(223)
    cax = ax.matshow(
        np.flip((DuneCrest),1) * 10,
        origin="lower",
        cmap="bwr",
        aspect="auto",
        vmin=vmin,
        vmax=vmax,
    )
    ax.xaxis.set_ticks_position("bottom")  # analysis:ignore
    # cax = ax.xaxis.set_ticks_position("bottom")  # analysis:ignore
    cbar = duneFig.colorbar(cax)
    cbar.set_label('Dune Height Above Berm Elevation (m)', rotation=270, labelpad=25)
    plt.xlabel("Alongshore Distance (m)")
    plt.ylabel("Year")
    plt.title("50% washout to shoreface", weight="bold")
    plt.xticks(x_ticks, x_tick_labels)

    # cascade 0%
    DuneCrest = []

    for iB3D in range(len(b3d)):
        sub_domain = outwash0[iB3D]._DuneDomain[0:TMAX, :, :]
        DuneCrest.append(sub_domain.max(axis=2))

    DuneCrest = np.hstack(DuneCrest).astype(float)

    plt.rcParams.update({"font.size": fontsize})
    ax = duneFig.add_subplot(222)
    cax = ax.matshow(
        np.flip(DuneCrest,1) * 10,
        origin="lower",
        cmap="bwr",
        aspect="auto",
        vmin=vmin,
        vmax=vmax,
    )
    ax.xaxis.set_ticks_position("bottom")  # analysis:ignore
    # cax = ax.xaxis.set_ticks_position("bottom")  # analysis:ignore
    cbar = duneFig.colorbar(cax)
    cbar.set_label('Dune Height Above Berm Elevation (m)', rotation=270, labelpad=25)
    plt.xlabel("Alongshore Distance (m)")
    plt.ylabel("Year")
    plt.title("0% washout to shoreface", weight="bold")
    plt.xticks(x_ticks, x_tick_labels)
    # plt.hlines(20, -0.5, 49.5, color="k", linestyles='dashed', linewidth=1)


    plt.subplots_adjust(hspace=0.5, wspace=0.2)