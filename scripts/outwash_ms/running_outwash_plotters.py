# running plotters_outwash to make dune and elevation figures

# Lexi Van Blunk
# 2/23/2024
# 3/16/2025 updated to the new directory (re-runs)

from matplotlib import pyplot as plt
import numpy as np
import math
from cascade.tools.outwash_plotters import plot_ElevAnimation_CASCADE as plt_CASCADE
from cascade.tools.outwash_plotters import plot_Elev_CASCADE_subplots as subplts_CASCADE

elevation_plots_on = False
dune_plots_on = False
animated_plotters_on = False
CASCADE_elev_subplots_on = True

rname_array = ["r025", "r035"]

for rname in rname_array:
    storm_interval = 20        # 20 or 10 years (we did not do 10 years)
    config = 4                 # 1, 2, 3, or 4
    storm_num = 1              # we ran 100 storms, so this can be any integer 1-100 (inclusive)

    # location of the npz files
    datadir_b3d = "C:/Users/Lexi/PycharmProjects/CASCADE/data/outwash_data/storms/slope0pt03/rerun_output/{0}/overwash_only/".format(rname)
    datadir_100 = "C:/Users/Lexi/PycharmProjects/CASCADE/data/outwash_data/storms/slope0pt03/rerun_output/{0}/outwash100/".format(rname)
    datadir_50 = "C:/Users/Lexi/PycharmProjects/CASCADE/data/outwash_data/storms/slope0pt03/rerun_output/{0}/outwash50/".format(rname)
    datadir_0 = "C:/Users/Lexi/PycharmProjects/CASCADE/data/outwash_data/storms/slope0pt03/rerun_output/{0}/outwash0/".format(rname)


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

    if elevation_plots_on:
        TMAX = 101
        vmin = -3
        vmax = 6
        fontsize = 12
        shrink = 0.5

        fig1 = plt.figure()
        fig1.suptitle('Baseline, {0}'.format(rname), weight="bold")
        # text is distance along the horizontal axis (left to right), then distance on the vertical axis (bottom to top)
        fig1.text(0.5, 0.42, 'barrier length (m)', ha='center', va='center', fontsize=12)  # x label
        fig1.text(0.08, 0.72, 'barrier width (m)', ha='center', va='center', rotation='vertical',
                  fontsize=12)  # y label

        fig2 = plt.figure()
        fig2.suptitle('100% Outwash, {0}'.format(rname), weight="bold")
        fig2.text(0.5, 0.42, 'barrier length (m)', ha='center', va='center', fontsize=12)  # x label
        fig2.text(0.08, 0.72, 'barrier width (m)', ha='center', va='center', rotation='vertical',
                  fontsize=12)  # y label

        fig3 = plt.figure()
        fig3.suptitle('50% Outwash, {0}'.format(rname), weight="bold")
        fig3.text(0.5, 0.42, 'barrier length (m)', ha='center', va='center', fontsize=12)  # x label
        fig3.text(0.08, 0.72, 'barrier width (m)', ha='center', va='center', rotation='vertical',
                  fontsize=12)  # y label

        fig4 = plt.figure()
        plt.axis('off')
        fig4.suptitle('0% Outwash, {0}'.format(rname), weight="bold")
        fig4.text(0.5, 0.42, 'barrier length (m)', ha='center', va='center', fontsize=12)  # x label
        fig4.text(0.08, 0.72, 'barrier width (m)', ha='center', va='center', rotation='vertical',
                  fontsize=12)  # y label

        years = [0, 1, 20, 21, 40, 41, 60, 61, 80, 81, 100]
        n_plots = len(years)
        plot_num = 1
        n_rows = 2
        n_cols = int(math.ceil(n_plots/2))

        for year in years:
            # overwash only
            if year <= b3d[0].TMAX:
                dunes = np.transpose(b3d[0]._DuneDomain[year]) + b3d[0].BermEl
                interior = b3d[0]._DomainTS[year]
                domain = np.vstack([dunes, interior])

                ax1 = fig1.add_subplot(n_rows, n_cols, plot_num)
                mat = ax1.matshow(
                    np.flip(domain) * 10,
                    cmap="terrain",
                    vmin=vmin,
                    vmax=vmax,
                )
                ax1.set_title(year)
                plt.gca().xaxis.tick_bottom()
                xtick_max = np.shape(domain)[1]  # n_cols = x
                x_ticks = np.array(range(0, xtick_max, 10))
                x_tick_labels = x_ticks * 10
                ytick_max = np.shape(domain)[0]  # n_rows = y
                y_ticks = np.array(range(0, ytick_max, 10))
                y_tick_labels = y_ticks * 10
                plt.xticks(x_ticks, x_tick_labels)
                plt.yticks(y_ticks, y_tick_labels)

            # 100% Outwash
            if year <= outwash100[0].TMAX:
                dunes = np.transpose(outwash100[0]._DuneDomain[year]) + outwash100[0].BermEl
                interior = outwash100[0]._DomainTS[year]
                domain = np.vstack([dunes, interior])

                ax1 = fig2.add_subplot(n_rows, n_cols, plot_num)
                mat = ax1.matshow(
                    np.flip(domain) * 10,
                    cmap="terrain",
                    vmin=vmin,
                    vmax=vmax,
                )
                ax1.set_title(year)
                plt.gca().xaxis.tick_bottom()
                xtick_max = np.shape(domain)[1]  # n_cols = x
                x_ticks = np.array(range(0, xtick_max, 10))
                x_tick_labels = x_ticks * 10
                ytick_max = np.shape(domain)[0]  # n_rows = y
                y_ticks = np.array(range(0, ytick_max, 10))
                y_tick_labels = y_ticks * 10
                plt.xticks(x_ticks, x_tick_labels)
                plt.yticks(y_ticks, y_tick_labels)

            # 50% Outwash
            if year <= outwash50[0].TMAX:
                dunes = np.transpose(outwash50[0]._DuneDomain[year]) + outwash50[0].BermEl
                interior = outwash50[0]._DomainTS[year]
                domain = np.vstack([dunes, interior])

                ax1 = fig3.add_subplot(n_rows, n_cols, plot_num)
                mat = ax1.matshow(
                    np.flip(domain) * 10,
                    cmap="terrain",
                    vmin=vmin,
                    vmax=vmax,
                )
                ax1.set_title(year)
                plt.gca().xaxis.tick_bottom()
                xtick_max = np.shape(domain)[1]  # n_cols = x
                x_ticks = np.array(range(0, xtick_max, 10))
                x_tick_labels = x_ticks * 10
                ytick_max = np.shape(domain)[0]  # n_rows = y
                y_ticks = np.array(range(0, ytick_max, 10))
                y_tick_labels = y_ticks * 10
                plt.xticks(x_ticks, x_tick_labels)
                plt.yticks(y_ticks, y_tick_labels)

            # 0% Outwash
            if year <= outwash0[0].TMAX:
                dunes = np.transpose(outwash0[0]._DuneDomain[year]) + outwash0[0].BermEl
                interior = outwash0[0]._DomainTS[year]
                domain = np.vstack([dunes, interior])

                ax1 = fig4.add_subplot(n_rows, n_cols, plot_num)
                mat = ax1.matshow(
                    np.flip(domain) * 10,
                    cmap="terrain",
                    vmin=vmin,
                    vmax=vmax,
                )
                ax1.set_title(year)
                plt.gca().xaxis.tick_bottom()
                xtick_max = np.shape(domain)[1]  # n_cols = x
                x_ticks = np.array(range(0, xtick_max, 10))
                x_tick_labels = x_ticks * 10
                ytick_max = np.shape(domain)[0]  # n_rows = y
                y_ticks = np.array(range(0, ytick_max, 10))
                y_tick_labels = y_ticks * 10
                plt.xticks(x_ticks, x_tick_labels)
                plt.yticks(y_ticks, y_tick_labels)

            plot_num += 1

        horiz_loc = 0.93
        vert_loc = 0.45
        width = 0.01
        height = 0.45

        cb_ax = fig1.add_axes([horiz_loc, vert_loc, width, height])
        cbar = fig1.colorbar(mat, cax=cb_ax)
        cbar.set_label('m MHW', rotation=270, labelpad=5)
        fig1.subplots_adjust(top=1.25, hspace=-0.75, wspace=0.5)

        cb_ax = fig2.add_axes([horiz_loc, vert_loc, width, height])
        cbar = fig2.colorbar(mat, cax=cb_ax)
        cbar.set_label('m MHW', rotation=270, labelpad=5)
        fig2.subplots_adjust(top=1.25, hspace=-0.75, wspace=0.5)

        cb_ax = fig3.add_axes([horiz_loc, vert_loc, width, height])
        cbar = fig3.colorbar(mat, cax=cb_ax)
        cbar.set_label('m MHW', rotation=270, labelpad=5)
        fig3.subplots_adjust(top=1.25, hspace=-0.75, wspace=0.5)

        cb_ax = fig4.add_axes([horiz_loc, vert_loc, width, height])
        cbar = fig4.colorbar(mat, cax=cb_ax)
        cbar.set_label('m MHW', rotation=270, labelpad=5)
        fig4.subplots_adjust(top=1.25, hspace=-0.75, wspace=0.5)


    ### ----------------------------------- Dune Crest Figures ---------------------------------------

    if dune_plots_on:

        TMAX = 101
        vmin = 0
        vmax = 6
        line_color = "white"  # line to mark outwash years

        # Barrier3d only
        DuneCrest = []

        for iB3D in range(len(b3d)):
            sub_domain = b3d[iB3D]._DuneDomain[0:TMAX, :, :]
            DuneCrest.append(sub_domain.max(axis=2))

        DuneCrest = np.hstack(DuneCrest).astype(float)

        duneFig = plt.figure(figsize=(15, 10))
        duneFig.suptitle('{0}'.format(rname), weight="bold", fontsize=12)
        plt.rcParams.update({"font.size": 12})
        ax = duneFig.add_subplot(221)
        cax = ax.matshow(
            np.flip(DuneCrest, 1) * 10,
            origin="lower",
            cmap="bwr",
            aspect="auto",
            vmin=vmin,
            vmax=vmax,
        )
        ax.xaxis.set_ticks_position("bottom")  # analysis:ignore
        # cax = ax.xaxis.set_ticks_position("bottom")  # analysis:ignore
        cbar = duneFig.colorbar(cax)
        cbar.set_label('Dune Height Above Berm Elevation (m)', rotation=270, labelpad=25, fontsize=10)
        plt.xlabel("Alongshore Distance (m)")
        plt.ylabel("Year")
        plt.title("baseline", weight="bold")
        # plt.hlines(20, -0.5, 49.5, color="k", linestyles='dashed', linewidth=1)

        xtick_max = np.shape(DuneCrest)[1]  # n_cols = x
        x_ticks = np.array(range(0, xtick_max, 10))
        x_tick_labels = x_ticks * 10
        plt.xticks(x_ticks, x_tick_labels)

        # outwash 100%
        DuneCrest = []

        for iB3D in range(len(b3d)):
            sub_domain = outwash100[iB3D]._DuneDomain[0:TMAX, :, :]
            DuneCrest.append(sub_domain.max(axis=2))

        DuneCrest = np.hstack(DuneCrest).astype(float)

        plt.rcParams.update({"font.size": 12})
        ax = duneFig.add_subplot(222)
        cax = ax.matshow(
            np.flip(DuneCrest, 1) * 10,
            origin="lower",
            cmap="bwr",
            aspect="auto",
            vmin=vmin,
            vmax=vmax,
        )
        ax.xaxis.set_ticks_position("bottom")  # analysis:ignore
        # cax = ax.xaxis.set_ticks_position("bottom")  # analysis:ignore
        cbar = duneFig.colorbar(cax)
        cbar.set_label('Dune Height Above Berm Elevation (m)', rotation=270, labelpad=25, fontsize=10)
        plt.xlabel("Alongshore Distance (m)")
        plt.ylabel("Year")
        plt.title("100% outwash to shoreface", weight="bold")
        plt.xticks(x_ticks, x_tick_labels)
        plt.hlines([1, 21, 41, 61, 81], -0.5, 49.5, colors=line_color, linestyles='solid')

        # outwash 50%
        DuneCrest = []

        for iB3D in range(len(b3d)):
            sub_domain = outwash50[iB3D]._DuneDomain[0:TMAX, :, :]
            DuneCrest.append(sub_domain.max(axis=2))

        DuneCrest = np.hstack(DuneCrest).astype(float)

        plt.rcParams.update({"font.size": 12})
        ax = duneFig.add_subplot(223)
        cax = ax.matshow(
            np.flip((DuneCrest), 1) * 10,
            origin="lower",
            cmap="bwr",
            aspect="auto",
            vmin=vmin,
            vmax=vmax,
        )
        ax.xaxis.set_ticks_position("bottom")  # analysis:ignore
        # cax = ax.xaxis.set_ticks_position("bottom")  # analysis:ignore
        cbar = duneFig.colorbar(cax)
        cbar.set_label('Dune Height Above Berm Elevation (m)', rotation=270, labelpad=25, fontsize=10)
        plt.xlabel("Alongshore Distance (m)")
        plt.ylabel("Year")
        plt.title("50% outwash to shoreface", weight="bold")
        plt.xticks(x_ticks, x_tick_labels)
        plt.hlines([1, 21, 41, 61, 81], -0.5, 49.5, colors=line_color, linestyles='solid')

        # outwash 0%
        DuneCrest = []

        for iB3D in range(len(b3d)):
            sub_domain = outwash0[iB3D]._DuneDomain[0:TMAX, :, :]
            DuneCrest.append(sub_domain.max(axis=2))

        DuneCrest = np.hstack(DuneCrest).astype(float)

        plt.rcParams.update({"font.size": 12})
        ax = duneFig.add_subplot(224)
        cax = ax.matshow(
            np.flip(DuneCrest, 1) * 10,
            origin="lower",
            cmap="bwr",
            aspect="auto",
            vmin=vmin,
            vmax=vmax,
        )
        ax.xaxis.set_ticks_position("bottom")  # analysis:ignore
        # cax = ax.xaxis.set_ticks_position("bottom")  # analysis:ignore
        cbar = duneFig.colorbar(cax)
        cbar.set_label('Dune Height Above Berm Elevation (m)', rotation=270, labelpad=25, fontsize=10)
        plt.xlabel("Alongshore Distance (m)")
        plt.ylabel("Year")
        plt.title("0% outwash to shoreface", weight="bold")
        plt.xticks(x_ticks, x_tick_labels)
        if rname == "r025":
            plt.hlines([1, 21], -0.5, 49.5, colors=line_color, linestyles='solid')
        else:
            plt.hlines([1, 21, 41], -0.5, 49.5, colors=line_color, linestyles='solid')
        # plt.hlines(20, -0.5, 49.5, color="k", linestyles='dashed', linewidth=1)

        plt.subplots_adjust(hspace=0.5, wspace=0.3)

    ### ----------------------------------- Animated CASCADE figures ---------------------------------------

    if animated_plotters_on:
        ymin = 125
        ymax = 225
        # Baseline
        plt_CASCADE(
            cascade=b3d_obj,
            directory=r"C:\Users\Lexi\Documents\Research\Outwasher Paper\figures\CASCADE_elevations\{0}".format(rname),
            TMAX_MGMT=b3d[0].TMAX,
            name="beach_baseline/",
            TMAX_SIM=b3d[0].TMAX,
            ny=1,
            beach_management_ny=None,  # list of bool the length of ny, or None for all False
            roadway_management_ny=None,
            y_lim=(ymin,ymax),
            z_lim=6,
            z_bot=-3,
            fig_size=None,
            fig_eps=False,
            km_on=False,
            meters_on=True,
        )

        # 100% outwash to shoreface
        plt_CASCADE(
            cascade=outwash100_obj,
            directory=r"C:\Users\Lexi\Documents\Research\Outwasher Paper\figures\CASCADE_elevations\{0}".format(rname),
            TMAX_MGMT=outwash100[0].TMAX,
            name="beach_outwash100/",
            TMAX_SIM=outwash100[0].TMAX,
            ny=1,
            beach_management_ny=None,  # list of bool the length of ny, or None for all False
            roadway_management_ny=None,
            y_lim=(ymin,ymax),
            z_lim=6,
            z_bot=-3,
            fig_size=None,
            fig_eps=False,
            km_on=False,
            meters_on=True,
        )

        # 50% outwash to shoreface
        plt_CASCADE(
            cascade=outwash50_obj,
            directory=r"C:\Users\Lexi\Documents\Research\Outwasher Paper\figures\CASCADE_elevations\{0}".format(rname),
            TMAX_MGMT=outwash50[0].TMAX,
            name="beach_outwash50/",
            TMAX_SIM=outwash50[0].TMAX,
            ny=1,
            beach_management_ny=None,  # list of bool the length of ny, or None for all False
            roadway_management_ny=None,
            y_lim=(ymin,ymax),
            z_lim=6,
            z_bot=-3,
            fig_size=None,
            fig_eps=False,
            km_on=False,
            meters_on=True,
        )

        # 0% outwash to shoreface
        plt_CASCADE(
            cascade=outwash0_obj,
            directory=r"C:\Users\Lexi\Documents\Research\Outwasher Paper\figures\CASCADE_elevations\{0}".format(rname),
            TMAX_MGMT=outwash0[0].TMAX,
            name="beach_outwash0/",
            TMAX_SIM=outwash0[0].TMAX,
            ny=1,
            beach_management_ny=None,  # list of bool the length of ny, or None for all False
            roadway_management_ny=None,
            y_lim=(ymin,ymax),
            z_lim=6,
            z_bot=-3,
            fig_size=None,
            fig_eps=False,
            km_on=False,
            meters_on=True,
        )

    ### ----------------------------------- CASCADE elevation subplots ---------------------------------------
    # CASCADE elevation plots
    if CASCADE_elev_subplots_on:
        ymin = 150
        ymax = 210
        # plot_years = [0,1,20,21,40,41,60,61,80,81,100]
        # plot_years = [0,1,20,21,40,41,100]
        # plot_years = [0,1,21,25,41,50,61,81,100]
        plot_years = [0, 1, 21, 41, 61, 81, 100]

        # Baseline
        subplts_CASCADE(
            cascade=b3d_obj,
            directory=r"C:\Users\Lexi\Documents\Research\Outwasher Paper\figures\elev_subplots\{0}".format(rname),
            TMAX_MGMT=b3d[0].TMAX,
            name="baseline_outwashtext_outyrs",
            TMAX_SIM=b3d[0].TMAX,
            plot_timesteps=plot_years,  # array with the timesteps you want to plot
            main_plot_title="baseline",
            ny=1,
            beach_management_ny=None,  # list of bool the length of ny, or None for all False
            roadway_management_ny=None,
            y_lim=(ymin, ymax),
            z_lim=6,
            z_bot=-3,
            fig_size=None,
            fig_eps=False,
            km_on=False,
            meters_on=True,
        )

        # 100% outwash
        subplts_CASCADE(
            cascade=outwash100_obj,
            directory=r"C:\Users\Lexi\Documents\Research\Outwasher Paper\figures\elev_subplots\{0}".format(rname),
            TMAX_MGMT=outwash100[0].TMAX,
            name="outwash100_outwashtext_outyrs",
            TMAX_SIM=outwash100[0].TMAX,
            plot_timesteps=plot_years,  # array with the timesteps you want to plot
            main_plot_title="100% outwash to shoreface",
            ny=1,
            beach_management_ny=None,  # list of bool the length of ny, or None for all False
            roadway_management_ny=None,
            y_lim=(ymin, ymax),
            z_lim=6,
            z_bot=-3,
            fig_size=None,
            fig_eps=False,
            km_on=False,
            meters_on=True,
        )

        # 50% outwash
        subplts_CASCADE(
            cascade=outwash50_obj,
            directory=r"C:\Users\Lexi\Documents\Research\Outwasher Paper\figures\elev_subplots\{0}".format(rname),
            TMAX_MGMT=outwash50[0].TMAX,
            name="outwash50_outwashtext_outyrs",
            TMAX_SIM=outwash50[0].TMAX,
            plot_timesteps=plot_years,  # array with the timesteps you want to plot
            main_plot_title="50% outwash to shoreface",
            ny=1,
            beach_management_ny=None,  # list of bool the length of ny, or None for all False
            roadway_management_ny=None,
            y_lim=(ymin, ymax),
            z_lim=6,
            z_bot=-3,
            fig_size=None,
            fig_eps=False,
            km_on=False,
            meters_on=True,
        )

        # 0%
        subplts_CASCADE(
            cascade=outwash0_obj,
            directory=r"C:\Users\Lexi\Documents\Research\Outwasher Paper\figures\elev_subplots\{0}".format(rname),
            TMAX_MGMT=outwash0[0].TMAX,
            name="outwash0_outwashtext_outyrs",
            TMAX_SIM=outwash0[0].TMAX,
            plot_timesteps=plot_years,  # array with the timesteps you want to plot
            main_plot_title="0% outwash to shoreface",
            ny=1,
            beach_management_ny=None,  # list of bool the length of ny, or None for all False
            roadway_management_ny=None,
            y_lim=(ymin, ymax),
            z_lim=6,
            z_bot=-3,
            fig_size=None,
            fig_eps=False,
            km_on=False,
            meters_on=True,
        )