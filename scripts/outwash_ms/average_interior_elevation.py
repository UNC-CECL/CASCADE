# loading in the npz files and saving variables that we want to plot together

# Lexi Van Blunk
# 5/15/2025
# calculating the average interior elevation using DomainTS instead of the h_b_TS variable because I think it is
# messed up

import numpy as np
from matplotlib import pyplot as plt

# ---------------------------------- set model parameters that change per run ------------------------------------------
rname_array = ["r025", "r035"]
# rname_array = ["r035"]
for rname in rname_array:
    storm_interval = 20        # 20 or 10 years
    config = 4                 # 1, 2, 3, or 4

    # Display stats on console/show plots
    migration_stats = False
    plotters = False
    geomoetry_stats = True
    dune_crest_stats = False
    flux_stats = False

    # location of the npz files
    datadir_b3d = "C:/Users/Lexi/PycharmProjects/CASCADE/data/outwash_data/storms/slope0pt03/rerun_output/{0}/overwash_only/".format(rname)
    datadir_100 = "C:/Users/Lexi/PycharmProjects/CASCADE/data/outwash_data/storms/slope0pt03/rerun_output/{0}/outwash100/".format(rname)
    datadir_50 = "C:/Users/Lexi/PycharmProjects/CASCADE/data/outwash_data/storms/slope0pt03/rerun_output/{0}/outwash50/".format(rname)
    datadir_0 = "C:/Users/Lexi/PycharmProjects/CASCADE/data/outwash_data/storms/slope0pt03/rerun_output/{0}/outwash0/".format(rname)

    # initialize empty arrays
    tmax_array_b3d = []
    tmax_array_100 = []
    tmax_array_50 = []
    tmax_array_0 = []

    tmax_first_drown_100 = []
    tmax_first_drown_50 = []
    tmax_first_drown_0 = []

    avg_elev_TS_array_b3d = []
    avg_elev_TS_array_100 = []
    avg_elev_TS_array_50 = []
    avg_elev_TS_array_0 = []

    avg_elev_b3d_array = np.zeros(100, dtype=object)
    avg_elev_100_array = np.zeros(100, dtype=object)
    avg_elev_50_array = np.zeros(100, dtype=object)
    avg_elev_0_array = np.zeros(100, dtype=object)

    avg_last_elev_b3d_array = np.zeros(100, dtype=object)
    avg_last_elev_100_array = np.zeros(100, dtype=object)
    avg_last_elev_50_array = np.zeros(100, dtype=object)
    avg_last_elev_0_array = np.zeros(100, dtype=object)

    for storm_num in range(1, 101):
        # print(storm_num)
        # b3d variables
        filename_b3d = "config{0}_b3d_startyr1_interval{1}yrs_Slope0pt03_{2}.npz".format(config, storm_interval, storm_num)
        file_b3d = datadir_b3d + filename_b3d
        b3d = np.load(file_b3d, allow_pickle=True)
        b3d_obj = b3d["cascade"][0]
        tmax_b3d = b3d_obj.barrier3d[0].TMAX
        tmax_array_b3d.append(tmax_b3d)

        for t in range(0, tmax_b3d):
            # average elevation
            domain_TS = b3d_obj.barrier3d[0].DomainTS[t] * 10
            avg_elev_TS = np.average(domain_TS)
            avg_elev_TS_array_b3d.append(avg_elev_TS)

        last_domain_TS_avg = np.average(domain_TS)
        avg_elev_b3d = np.average(avg_elev_TS_array_b3d)
        avg_elev_b3d_array[storm_num-1] = avg_elev_b3d
        avg_last_elev_b3d_array[storm_num-1] = last_domain_TS_avg

        # 100% variables
        filename_100 = "config{0}_outwash100_startyr1_interval{1}yrs_Slope0pt03_{2}.npz".format(config, storm_interval, storm_num)
        file_100 = datadir_100 + filename_100
        outwash100 = np.load(file_100, allow_pickle=True)
        outwash100_obj = outwash100["cascade"][0]
        tmax_100 = outwash100_obj.barrier3d[0].TMAX
        tmax_array_100.append(tmax_100)

        # average elevation
        if tmax_100 < 101:
            max_t = tmax_100 + 1
        else:
            max_t = tmax_100
        for t in range(0, max_t):
            # average elevation
            domain_TS = outwash100_obj.barrier3d[0].DomainTS[t] * 10
            avg_elev_TS = np.average(domain_TS)
            avg_elev_TS_array_100.append(avg_elev_TS)

        last_domain_TS_avg = np.average(domain_TS)  #domain TS will be the domain before drowning
        avg_elev_100 = np.average(avg_elev_TS_array_100)
        avg_elev_100_array[storm_num-1] = avg_elev_100
        avg_last_elev_100_array[storm_num-1] = last_domain_TS_avg

        # 50% variables
        filename_50 = "config{0}_outwash50_startyr1_interval{1}yrs_Slope0pt03_{2}.npz".format(config, storm_interval, storm_num)
        file_50 = datadir_50 + filename_50
        outwash50 = np.load(file_50, allow_pickle=True)
        outwash50_obj = outwash50["cascade"][0]
        tmax_50 = outwash50_obj.barrier3d[0].TMAX
        tmax_array_50.append(tmax_50)

        # average elevation
        if tmax_50 < 101:
            max_t = tmax_50 + 1
        else:
            max_t = tmax_50
        for t in range(0, max_t):
            # average elevation
            domain_TS = outwash50_obj.barrier3d[0].DomainTS[t] * 10
            avg_elev_TS = np.average(domain_TS)
            avg_elev_TS_array_50.append(avg_elev_TS)

        last_domain_TS_avg = np.average(domain_TS)
        avg_elev_50 = np.average(avg_elev_TS_array_50)
        avg_elev_50_array[storm_num-1] = avg_elev_50
        avg_last_elev_50_array[storm_num-1] = last_domain_TS_avg

        # 0% variables
        filename_0 = "config{0}_outwash0_startyr1_interval{1}yrs_Slope0pt03_{2}.npz".format(config, storm_interval, storm_num)
        file_0 = datadir_0 + filename_0
        outwash0 = np.load(file_0, allow_pickle=True)
        outwash0_obj = outwash0["cascade"][0]
        tmax_0 = outwash0_obj.barrier3d[0].TMAX
        tmax_array_0.append(tmax_0)

        # average elevation
        if tmax_0 < 101:
            max_t = tmax_0 + 1
        else:
            max_t = tmax_0
        for t in range(0, max_t):
            # average elevation
            domain_TS = outwash0_obj.barrier3d[0].DomainTS[t] * 10
            avg_elev_TS = np.average(domain_TS)
            avg_elev_TS_array_0.append(avg_elev_TS)

        # last_domain_TS_avg = np.average(outwash0_obj.barrier3d[0].DomainTS[-1]*10)
        last_domain_TS_avg = np.average(domain_TS)
        avg_elev_0 = np.average(avg_elev_TS_array_0)
        avg_elev_0_array[storm_num-1] = avg_elev_0
        avg_last_elev_0_array[storm_num-1] = last_domain_TS_avg


    # printing geometry stats
    if geomoetry_stats:
        # print avg height and width stats
        print("avg geometry stats, {0}".format(rname))
        print("baseline \n avg interior height: {0} \n avg last interior height: {1}"
              .format(np.round(np.average(avg_elev_b3d_array), 2),
                      np.round(np.average(avg_last_elev_b3d_array), 2)))
        print("100% outwash \n avg interior height: {0} \n avg last interior height: {1}"
              .format(np.round(np.average(avg_elev_100_array), 2),
                      np.round(np.average(avg_last_elev_100_array), 2)))
        print("50% outwash \n avg interior height: {0} \n avg last interior height: {1}"
              .format(np.round(np.average(avg_elev_50_array), 2),
                      np.round(np.average(avg_last_elev_50_array), 2)))
        print("0% outwash \n avg interior height: {0} \n avg last interior height: {1}"
              .format(np.round(np.average(avg_elev_0_array), 2),
                      np.round(np.average(avg_last_elev_0_array), 2)))


    # if plotters:
    #     storm_num = 1
    #     fig9 = plt.figure()
    #     # fig9.suptitle('overwash storm {0} - {1}'.format(storm_num, rname), weight="bold")
    #     # fig9.suptitle('{0}'.format(rname), weight="bold")
    #     ax1 = fig9.add_subplot(111)
    #     ls = "solid"
    #     ax1.plot(shoreline_pos_array_b3d[storm_num-1])
    #     ax1.plot(shoreline_pos_array_100[storm_num-1], linestyle=ls)
    #     ax1.plot(shoreline_pos_array_50[storm_num-1], linestyle=ls)
    #     ax1.plot(shoreline_pos_array_0[storm_num-1], linestyle=ls)
    #     ax1.legend(["baseline", "100% outwash to shoreface", "50% outwash to shoreface", "0% outwash to shoreface"],
    #                prop={'size': 9}, loc="upper left")
    #     ax1.set_ylabel("Shoreline Position (m)")
    #     ax1.set_xlabel("Simulation Years")
    #     ax1.set_ylim(bottom=-60, top=160)
    #     ax1.set_title('{0}'.format(rname), weight="bold")
    #     fig9.subplots_adjust(hspace=0.3)
    #
    #     # plotting one year of overwash and outwash flux
    #     fig10 = plt.figure()
    #     # fig9.suptitle('overwash storm {0} - {1}'.format(storm_num, rname), weight="bold")
    #     # fig9.suptitle('{0}'.format(rname), weight="bold")
    #     ax1 = fig10.add_subplot(211)
    #     ls = "solid"
    #     ax1.plot(overwash_array_b3d[storm_num-1])
    #     ax1.plot(overwash_array_100[storm_num-1], linestyle=ls)
    #     ax1.plot(overwash_array_50[storm_num-1], linestyle=ls)
    #     ax1.plot(overwash_array_0[storm_num-1], linestyle=ls)
    #     ax1.set_ylabel("Overwash Flux (m3/m)")
    #     ax1.set_xlabel("Simulation Years")
    #     ax1.set_ylim(top=150)
    #     ax1.set_title(rname, weight="bold")
    #     ax1.legend(["baseline", "100% outwash to shoreface", "50% outwash to shoreface", "0% outwash to shoreface"],
    #                prop={'size': 9}, loc="upper left")
    #
    #     max_t_b3d = tmax_array_b3d[storm_num-1]
    #     max_t_100 = tmax_array_100[storm_num - 1]
    #     max_t_50 = tmax_array_50[storm_num - 1]
    #     max_t_0 = tmax_array_0[storm_num - 1]
    #
    #     outwash_array_b3d_storm = outwash_array_b3d[storm_num-1]
    #     outwash_array_b3d_storm[max_t_b3d+1:] = np.nan
    #     outwash_array_100_storm = outwash_array_100[storm_num - 1]
    #     outwash_array_100_storm[max_t_100+1:] = np.nan
    #     outwash_array_50_storm = outwash_array_50[storm_num - 1]
    #     outwash_array_50_storm[max_t_50+1:] = np.nan
    #     outwash_array_0_storm = outwash_array_0[storm_num - 1]
    #     outwash_array_0_storm[max_t_0+1:] = np.nan




