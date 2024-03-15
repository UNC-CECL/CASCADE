# loading in the npz files and saving variables that we want to plot together

# Lexi Van Blunk
# 7/26/2023
# updated 2/23/2024

import numpy as np
from matplotlib import pyplot as plt

# ---------------------------------- set model parameters that change per run ------------------------------------------
rname_array = ["r025", "r035"]
for rname in rname_array:
    storm_interval = 20        # 20 or 10 years
    # rname = "r025"             # "r025" or "r035"
    config = 4                 # 1, 2, 3, or 4

    # location of the npz files
    datadir_b3d = "C:/Users/Lexi/PycharmProjects/CASCADE/cascade/data/outwash_data/storms/slope0pt03/run_output/{0}/overwash_only/".format(rname)
    datadir_100 = "C:/Users/Lexi/PycharmProjects/CASCADE/cascade/data/outwash_data/storms/slope0pt03/run_output/{0}/outwash100/".format(rname)
    datadir_50 = "C:/Users/Lexi/PycharmProjects/CASCADE/cascade/data/outwash_data/storms/slope0pt03/run_output/{0}/outwash50/".format(rname)
    datadir_0 = "C:/Users/Lexi/PycharmProjects/CASCADE/cascade/data/outwash_data/storms/slope0pt03/run_output/{0}/outwash0/".format(rname)

    # initialize empty arrays
    drowning_array_b3d = np.zeros(100)
    drowning_array_100 = np.zeros(100)
    drowning_array_50 = np.zeros(100)
    drowning_array_0 = np.zeros(100)

    drown_year_array_b3d = np.zeros(100)
    drown_year_array_100 = np.zeros(100)
    drown_year_array_50 = np.zeros(100)
    drown_year_array_0 = np.zeros(100)

    shoreline_pos_array_b3d = np.zeros(100, dtype=object)
    shoreline_pos_array_100 = np.zeros(100, dtype=object)
    shoreline_pos_array_50 = np.zeros(100, dtype=object)
    shoreline_pos_array_0 = np.zeros(100, dtype=object)

    shoreface_slope_array_b3d = np.zeros(100, dtype=object)
    shoreface_slope_array_100 = np.zeros(100, dtype=object)
    shoreface_slope_array_50 = np.zeros(100, dtype=object)
    shoreface_slope_array_0 = np.zeros(100, dtype=object)

    avg_int_height_array_b3d = np.zeros(100, dtype=object)
    avg_int_height_array_100 = np.zeros(100, dtype=object)
    avg_int_height_array_50 = np.zeros(100, dtype=object)
    avg_int_height_array_0 = np.zeros(100, dtype=object)

    avg_int_width_array_b3d = np.zeros(100, dtype=object)
    avg_int_width_array_100 = np.zeros(100, dtype=object)
    avg_int_width_array_50 = np.zeros(100, dtype=object)
    avg_int_width_array_0 = np.zeros(100, dtype=object)

    overwash_array_b3d = np.zeros(100, dtype=object)
    overwash_array_100 = np.zeros(100, dtype=object)
    overwash_array_50 = np.zeros(100, dtype=object)
    overwash_array_0 = np.zeros(100, dtype=object)

    outwash_array_b3d = np.zeros(100, dtype=object)
    outwash_array_100 = np.zeros(100, dtype=object)
    outwash_array_50 = np.zeros(100, dtype=object)
    outwash_array_0 = np.zeros(100, dtype=object)

    for storm_num in range(1, 101):
        # b3d variables
        filename_b3d = "config{0}_b3d_startyr1_interval{1}yrs_Slope0pt03_{2}.npz".format(config, storm_interval, storm_num)
        file_b3d = datadir_b3d + filename_b3d
        b3d = np.load(file_b3d, allow_pickle=True)
        b3d_obj = b3d["cascade"][0]
        # drowning variable
        drowning_array_b3d[storm_num-1] = b3d_obj.barrier3d[0].drown_break
        if b3d_obj.barrier3d[0].TMAX < 100:
            drown_year_array_b3d[storm_num-1] = b3d_obj.barrier3d[0].TMAX
        # shoreline position
        m_xsTS = np.subtract(b3d_obj.barrier3d[0].x_s_TS, b3d_obj.barrier3d[0].x_s_TS[0])
        m_xsTS = np.multiply(m_xsTS, 10)
        shoreline_pos_array_b3d[storm_num-1] = m_xsTS
        # shoreface slope
        sfTS = b3d_obj.barrier3d[0].s_sf_TS
        shoreface_slope_array_b3d[storm_num - 1] = sfTS
        # average elevation
        hbTS = np.array(b3d_obj.barrier3d[0].h_b_TS) * 10  # m MHW
        avg_int_height_array_b3d[storm_num - 1] = hbTS
        # width
        avg_int_width_TS = np.array(b3d_obj.barrier3d[0].InteriorWidth_AvgTS) * 10  # meters
        avg_int_width_array_b3d[storm_num - 1] = avg_int_width_TS
        # overwash
        QowTS = b3d_obj.barrier3d[0].QowTS
        overwash_array_b3d[storm_num - 1] = QowTS
        # outwash
        QoutTS = b3d_obj.outwash[0]._outwash_flux_TS
        outwash_array_b3d[storm_num - 1] = QoutTS

        # 100% variables
        filename_100 = "config{0}_outwash100_startyr1_interval{1}yrs_Slope0pt03_{2}.npz".format(config, storm_interval, storm_num)
        file_100 = datadir_100 + filename_100
        outwash100 = np.load(file_100, allow_pickle=True)
        outwash100_obj = outwash100["cascade"][0]
        drowning_array_100[storm_num-1] = outwash100_obj.barrier3d[0].drown_break
        if outwash100_obj.barrier3d[0].TMAX < 100:
            drown_year_array_100[storm_num-1] = outwash100_obj.barrier3d[0].TMAX
        # shoreline position
        m_xsTS100 = np.subtract(outwash100_obj.barrier3d[0].x_s_TS, outwash100_obj.barrier3d[0].x_s_TS[0])
        m_xsTS100 = np.multiply(m_xsTS100, 10)
        shoreline_pos_array_100[storm_num-1] = m_xsTS100
        # shoreface slope
        sfTS100 = outwash100_obj.barrier3d[0].s_sf_TS
        shoreface_slope_array_100[storm_num - 1] = sfTS100
        # average elevation
        hbTS100 = np.array(outwash100_obj.barrier3d[0].h_b_TS) * 10  # m MHW
        avg_int_height_array_100[storm_num - 1] = hbTS100
        # width
        avg_int_width_TS100 = np.array(outwash100_obj.barrier3d[0].InteriorWidth_AvgTS) * 10  # meters
        avg_int_width_array_100[storm_num - 1] = avg_int_width_TS100
        # overwash
        QowTS100 = outwash100_obj.barrier3d[0].QowTS
        overwash_array_100[storm_num - 1] = QowTS100
        # outwash
        QoutTS100 = outwash100_obj.outwash[0]._outwash_flux_TS
        outwash_array_100[storm_num - 1] = QoutTS100

        # 50% variables
        filename_50 = "config{0}_outwash50_startyr1_interval{1}yrs_Slope0pt03_{2}.npz".format(config, storm_interval, storm_num)
        file_50 = datadir_50 + filename_50
        outwash50 = np.load(file_50, allow_pickle=True)
        outwash50_obj = outwash50["cascade"][0]
        drowning_array_50[storm_num-1] = outwash50_obj.barrier3d[0].drown_break
        if outwash50_obj.barrier3d[0].TMAX < 100:
            drown_year_array_50[storm_num - 1] = outwash50_obj.barrier3d[0].TMAX
        # shoreline position
        m_xsTS50 = np.subtract(outwash50_obj.barrier3d[0].x_s_TS, outwash50_obj.barrier3d[0].x_s_TS[0])
        m_xsTS50 = np.multiply(m_xsTS50, 10)
        shoreline_pos_array_50[storm_num-1] = m_xsTS50
        # shoreface slope
        sfTS50 = outwash50_obj.barrier3d[0].s_sf_TS
        shoreface_slope_array_50[storm_num - 1] = sfTS50
        # average elevation
        hbTS50 = np.array(outwash50_obj.barrier3d[0].h_b_TS) * 10  # m MHW
        avg_int_height_array_50[storm_num - 1] = hbTS50
        # width
        avg_int_width_TS50 = np.array(outwash50_obj.barrier3d[0].InteriorWidth_AvgTS) * 10  # meters
        avg_int_width_array_50[storm_num - 1] = avg_int_width_TS50
        # overwash
        QowTS50 = outwash50_obj.barrier3d[0].QowTS
        overwash_array_50[storm_num - 1] = QowTS50
        # outwash
        QoutTS50 = outwash50_obj.outwash[0]._outwash_flux_TS
        outwash_array_50[storm_num - 1] = QoutTS50

        # 0% variables
        filename_0 = "config{0}_outwash0_startyr1_interval{1}yrs_Slope0pt03_{2}.npz".format(config, storm_interval, storm_num)
        file_0 = datadir_0 + filename_0
        outwash0 = np.load(file_0, allow_pickle=True)
        outwash0_obj = outwash0["cascade"][0]
        drowning_array_0[storm_num-1] = outwash0_obj.barrier3d[0].drown_break
        if outwash0_obj.barrier3d[0].TMAX < 100:
            drown_year_array_0[storm_num - 1] = outwash0_obj.barrier3d[0].TMAX
        # shoreline position
        m_xsTS0 = np.subtract(outwash0_obj.barrier3d[0].x_s_TS, outwash0_obj.barrier3d[0].x_s_TS[0])
        m_xsTS0 = np.multiply(m_xsTS0, 10)
        shoreline_pos_array_0[storm_num-1] = m_xsTS0
        # shoreface slope
        sfTS0 = outwash0_obj.barrier3d[0].s_sf_TS
        shoreface_slope_array_0[storm_num - 1] = sfTS0
        # average elevation
        hbTS0 = np.array(outwash0_obj.barrier3d[0].h_b_TS) * 10  # meters MHW
        avg_int_height_array_0[storm_num - 1] = hbTS0
        # width
        avg_int_width_TS0 = np.array(outwash0_obj.barrier3d[0].InteriorWidth_AvgTS) * 10  # meters
        avg_int_width_array_0[storm_num - 1] = avg_int_width_TS0
        # overwash
        QowTS0 = outwash0_obj.barrier3d[0].QowTS
        overwash_array_0[storm_num - 1] = QowTS0
        # outwash
        QoutTS0 = outwash0_obj.outwash[0]._outwash_flux_TS
        outwash_array_0[storm_num - 1] = QoutTS0

    # because we have 100 storms and drowns = 1, the sum of the array is the percent that drown
    percent_drown_b3d = np.sum(drowning_array_b3d)
    print("{0}% of barriers drown for the overwash only scenario".format(int(percent_drown_b3d)))

    percent_drown_100 = np.sum(drowning_array_100)
    print("{0}% of barriers drown for the 100% outwash to shoreface scenario".format(int(percent_drown_100)))

    percent_drown_50 = np.sum(drowning_array_50)
    print("{0}% of barriers drown for the 50% outwash to shoreface scenario".format(int(percent_drown_50)))

    percent_drown_0 = np.sum(drowning_array_0)
    print("{0}% of barriers drown for the 0% outwash to shoreface scenario".format(int(percent_drown_0)))

    # histogram of years that the barriers drown (three outwash scenarios only, each on separate plot)
    plt.rcParams.update({"font.size": 12})
    bins = 100
    fig1 = plt.figure()
    fig1.suptitle('{0}'.format(rname), weight="bold")
    ax1 = fig1.add_subplot(131)
    ax1.hist(drown_year_array_100, bins=bins)
    ax1.set_title("100% washout")
    ax1.set_ylabel("number of barriers that drown")
    ax1.set_xlabel("drown year")
    plt.gca().xaxis.tick_bottom()
    ax1.set_xlim(left=1, right=100)
    ax1.set_ylim(bottom=0, top=40)

    ax1 = fig1.add_subplot(132)
    ax1.hist(drown_year_array_50, bins=bins)
    ax1.set_title("50% washout")
    ax1.set_ylabel("frequency")
    ax1.set_xlabel("drown year")
    plt.gca().xaxis.tick_bottom()
    ax1.set_xlim(left=1, right=100)
    ax1.set_ylim(bottom=0, top=40)
    ax1.set_ylabel(None)

    ax1 = fig1.add_subplot(133)
    ax1.hist(drown_year_array_0, bins=bins)
    ax1.set_title("0% washout")
    ax1.set_ylabel("frequency")
    ax1.set_xlabel("drown year")
    plt.gca().xaxis.tick_bottom()
    ax1.set_xlim(left=1, right=100)
    ax1.set_ylim(bottom=0, top=40)
    ax1.set_ylabel(None)

    # plotting shoreline position for all 100 storms, each scenario on separate plot
    bot = -200
    top = 400
    fig3 = plt.figure()
    fig3.suptitle('{0}'.format(rname), weight="bold")
    ax1 = fig3.add_subplot(221)
    ax1.set_title("no washout")
    ax1.set_ylim(bottom=bot, top=top)
    ax2 = fig3.add_subplot(222)
    ax2.set_title("100% washout")
    ax2.set_ylim(bottom=bot, top=top)
    ax3 = fig3.add_subplot(223)
    ax3.set_title("50% washout")
    ax3.set_ylim(bottom=bot, top=top)
    ax4 = fig3.add_subplot(224)
    ax4.set_title("0% washout")
    ax4.set_ylim(bottom=bot, top=top)

    fig3.text(0.5, 0.04, 'Time (years)', ha='center')
    fig3.text(0.06, 0.5, 'Shoreline Position (m)', va='center', rotation='vertical')
    fig3.subplots_adjust(hspace=0.4)

    for x in range(100):
        ax1.plot(shoreline_pos_array_b3d[x])
        ax2.plot(shoreline_pos_array_100[x])
        ax3.plot(shoreline_pos_array_50[x])
        ax4.plot(shoreline_pos_array_0[x])

    # plotting first storm shoreline position, shoreface slope, avg int. height and width,
    # overwash and outwash flux, all four scenarios on same plot

    fig8 = plt.figure()
    fig8.suptitle('{0}'.format(rname), weight="bold")
    ax1 = fig8.add_subplot(231)
    ls = "dashed"
    ax1.plot(shoreline_pos_array_b3d[0])
    ax1.plot(shoreline_pos_array_100[0], linestyle=ls)
    ax1.plot(shoreline_pos_array_50[0], linestyle=ls)
    ax1.plot(shoreline_pos_array_0[0], linestyle=ls)
    ax1.legend(["no outwash", "100%", "50%", "0%"], prop={'size': 9})
    ax1.set_ylabel("Shoreline Position (m)")
    ax1.set_xlabel("Simulation Years")
    ax1.set_ylim(bottom=-60, top=160)

    ax2 = fig8.add_subplot(232)
    ax2.plot(shoreface_slope_array_b3d[0])
    ax2.plot(shoreface_slope_array_100[0], linestyle=ls)
    ax2.plot(shoreface_slope_array_50[0], linestyle=ls)
    ax2.plot(shoreface_slope_array_0[0], linestyle=ls)
    ax2.set_ylabel("Shoreface Slope")
    ax2.set_xlabel("Simulation Years")
    ax2.set_ylim(top=0.02)

    ax3 = fig8.add_subplot(235)
    ax3.plot(avg_int_height_array_b3d[0])
    ax3.plot(avg_int_height_array_100[0], linestyle=ls)
    ax3.plot(avg_int_height_array_50[0], linestyle=ls)
    ax3.plot(avg_int_height_array_0[0], linestyle=ls)
    ax3.set_ylabel("Average Interior Elevation (m MHW)")
    ax3.set_xlabel("Simulation Years")
    ax3.set_ylim(top=1.6)

    ax4 = fig8.add_subplot(234)
    ax4.plot(avg_int_width_array_b3d[0])
    ax4.plot(avg_int_width_array_100[0], linestyle=ls)
    ax4.plot(avg_int_width_array_50[0], linestyle=ls)
    ax4.plot(avg_int_width_array_0[0], linestyle=ls)
    ax4.set_ylabel("Average Interior Width (m)")
    ax4.set_xlabel("Simulation Years")
    ax4.set_ylim(top=310)

    ax5 = fig8.add_subplot(233)
    ax5.plot(overwash_array_b3d[0])
    ax5.plot(overwash_array_100[0], linestyle=ls)
    ax5.plot(overwash_array_50[0], linestyle=ls)
    ax5.plot(overwash_array_0[0], linestyle=ls)
    ax5.set_ylabel("Overwash Flux (m3/m)")
    ax5.set_xlabel("Simulation Years")
    ax5.set_ylim(top=150)

    ax6 = fig8.add_subplot(236)
    ax6.plot(outwash_array_b3d[0])
    ax6.plot(outwash_array_100[0], linestyle=ls)
    ax6.plot(outwash_array_50[0], linestyle=ls)
    ax6.plot(outwash_array_0[0], linestyle=ls)
    ax6.set_ylabel("Outwash Flux (m3/m)")
    ax6.set_xlabel("Simulation Years")
    ax6.set_ylim(top=1200)

    # fig8.text(0.5, 0.04, 'Time (years)', ha='center')
    # fig8.text(0.06, 0.5, 'Shoreline Position (m)', va='center', rotation='vertical')
    fig8.subplots_adjust(hspace=0.3, wspace=0.3)
    # fig8.tight_layout()
