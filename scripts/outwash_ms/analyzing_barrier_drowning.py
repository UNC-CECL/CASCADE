# loading in the npz files and saving variables that we want to plot together

# Lexi Van Blunk
# 7/26/2023
# updated 2/23/2024
# NOTE: I realized there are some instances where the code only behaves correctly because we have both 100 storms and
# 100 years in each scenario. There are a few loops inside loops that would not work if these values differed, because
# I am referencing two different things (model time years and number of storm runs) all using one loop (num storm runs)

import numpy as np
from matplotlib import pyplot as plt

# ---------------------------------- set model parameters that change per run ------------------------------------------
rname_array = ["r025", "r035"]
for rname in rname_array:
    storm_interval = 20        # 20 or 10 years
    config = 4                 # 1, 2, 3, or 4

    # Display stats on console/show plots
    migration_stats = False
    plotters = True
    geomoetry_stats = True

    # location of the npz files
    datadir_b3d = "C:/Users/Lexi/PycharmProjects/CASCADE/cascade/data/outwash_data/storms/slope0pt03/run_output/{0}/overwash_only/".format(rname)
    datadir_100 = "C:/Users/Lexi/PycharmProjects/CASCADE/cascade/data/outwash_data/storms/slope0pt03/run_output/{0}/outwash100/".format(rname)
    datadir_50 = "C:/Users/Lexi/PycharmProjects/CASCADE/cascade/data/outwash_data/storms/slope0pt03/run_output/{0}/outwash50/".format(rname)
    datadir_0 = "C:/Users/Lexi/PycharmProjects/CASCADE/cascade/data/outwash_data/storms/slope0pt03/run_output/{0}/outwash0/".format(rname)

    # initialize empty arrays
    tmax_array_b3d = []
    tmax_array_100 = []
    tmax_array_50 = []
    tmax_array_0 = []

    tmax_first_drown_100 = []
    tmax_first_drown_50 = []
    tmax_first_drown_0 = []

    drowning_array_b3d = np.zeros(100)
    drowning_array_100 = np.zeros(100)
    drowning_array_50 = np.zeros(100)
    drowning_array_0 = np.zeros(100)

    drown_year_array_b3d = []
    drown_year_array_100 = []
    drown_year_array_50 = []
    drown_year_array_0 = []

    shoreline_pos_array_b3d = np.zeros(100, dtype=object)
    shoreline_pos_array_100 = np.zeros(100, dtype=object)
    shoreline_pos_array_50 = np.zeros(100, dtype=object)
    shoreline_pos_array_0 = np.zeros(100, dtype=object)

    end_net_migration_b3d = []
    end_net_migration_100 = []
    end_net_migration_50 = []
    end_net_migration_0 = []

    end_net_migration_no_drowns_b3d = []
    end_net_migration_no_drowns_100 = []
    end_net_migration_no_drowns_50 = []
    end_net_migration_no_drowns_0 = []

    end_net_migration_only_drowns_b3d = []
    end_net_migration_only_drowns_100 = []
    end_net_migration_only_drowns_50 = []
    end_net_migration_only_drowns_0 = []

    shoreface_slope_array_b3d = np.zeros(100, dtype=object)
    shoreface_slope_array_100 = np.zeros(100, dtype=object)
    shoreface_slope_array_50 = np.zeros(100, dtype=object)
    shoreface_slope_array_0 = np.zeros(100, dtype=object)

    int_height_array_b3d = np.zeros(100, dtype=object)
    int_height_array_100 = np.zeros(100, dtype=object)
    int_height_array_50 = np.zeros(100, dtype=object)
    int_height_array_0 = np.zeros(100, dtype=object)

    avg_int_height_array_b3d = np.zeros(100, dtype=object)
    avg_int_height_array_100 = np.zeros(100, dtype=object)
    avg_int_height_array_50 = np.zeros(100, dtype=object)
    avg_int_height_array_0 = np.zeros(100, dtype=object)

    int_width_array_b3d = np.zeros(100, dtype=object)
    int_width_array_100 = np.zeros(100, dtype=object)
    int_width_array_50 = np.zeros(100, dtype=object)
    int_width_array_0 = np.zeros(100, dtype=object)

    avg_int_width_array_b3d = np.zeros(100, dtype=object)
    avg_int_width_array_100 = np.zeros(100, dtype=object)
    avg_int_width_array_50 = np.zeros(100, dtype=object)
    avg_int_width_array_0 = np.zeros(100, dtype=object)

    dune_crest_array_b3d = np.zeros(100, dtype=object)
    dune_crest_array_100 = np.zeros(100, dtype=object)
    dune_crest_array_50 = np.zeros(100, dtype=object)
    dune_crest_array_0 = np.zeros(100, dtype=object)

    avg_dune_crest_array_b3d = []
    avg_dune_crest_array_100 = []
    avg_dune_crest_array_50 = []
    avg_dune_crest_array_0 = []

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
        tmax_array_b3d.append(b3d_obj.barrier3d[0].TMAX)
        drowning_array_b3d[storm_num-1] = b3d_obj.barrier3d[0].drown_break
        if drowning_array_b3d[storm_num-1] == 1:
            drown_year_array_b3d.append(b3d_obj.barrier3d[0].TMAX)
        # shoreline position
        m_xsTS = np.subtract(b3d_obj.barrier3d[0].x_s_TS, b3d_obj.barrier3d[0].x_s_TS[0])
        m_xsTS = np.multiply(m_xsTS, 10)
        m_xsTS = m_xsTS[~np.isnan(m_xsTS)]
        shoreline_pos_array_b3d[storm_num-1] = m_xsTS
        if len(shoreline_pos_array_b3d[storm_num-1]) < 100:
            drowning_array_b3d[storm_num - 1] = 1
            drown_year_array_b3d.append(b3d_obj.barrier3d[0].TMAX)
        # average net migration including drownings
        end_net_migration_b3d.append(shoreline_pos_array_b3d[storm_num-1][-1])
        # average net migration excluding drownings
        if drowning_array_b3d[storm_num - 1] == 0:
            end_net_migration_no_drowns_b3d.append(shoreline_pos_array_b3d[storm_num-1][-1])
        else:  # average net migration only including drownings
            end_net_migration_only_drowns_b3d.append(shoreline_pos_array_b3d[storm_num - 1][-1])
        # shoreface slope
        sfTS = b3d_obj.barrier3d[0].s_sf_TS
        shoreface_slope_array_b3d[storm_num - 1] = sfTS
        # average elevation
        hbTS = np.array(b3d_obj.barrier3d[0].h_b_TS) * 10  # m MHW
        int_height_array_b3d[storm_num - 1] = hbTS
        avg_int_height_array_b3d[storm_num - 1] = np.mean(hbTS)
        # width
        int_width_TS = np.array(b3d_obj.barrier3d[0].InteriorWidth_AvgTS) * 10  # meters
        int_width_array_b3d[storm_num - 1] = int_width_TS
        avg_int_width_array_b3d[storm_num - 1] = np.mean(int_width_TS)
        # dune crest
        sub_domain = b3d_obj.barrier3d[0]._DuneDomain[0:tmax_array_b3d[storm_num-1], :, :]
        dune_crest_array_b3d[storm_num - 1] = sub_domain.max(axis=2)
        avg_dune_crest_array_b3d.append(np.average(dune_crest_array_b3d[storm_num - 1], axis=0))
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
        # shoreline position - used to also inform us of barrier drowning stats
        m_xsTS100 = np.subtract(outwash100_obj.barrier3d[0].x_s_TS, outwash100_obj.barrier3d[0].x_s_TS[0])
        m_xsTS100 = np.multiply(m_xsTS100, 10)
        m_xsTS100 = m_xsTS100[~np.isnan(m_xsTS100)]
        shoreline_pos_array_100[storm_num-1] = m_xsTS100
        tmax_first_drown_100 = len(m_xsTS100)  # based on the first drowning event if there is one

        if tmax_first_drown_100 < 101:
            drowning_array_100[storm_num - 1] = 1
            if tmax_first_drown_100 < outwash100_obj.barrier3d[0].TMAX:
                drown_year_array_100.append(tmax_first_drown_100)
                tmax_array_100.append(tmax_first_drown_100)
            else:
                drown_year_array_100.append(outwash100_obj.barrier3d[0].TMAX)
                tmax_array_100.append(outwash100_obj.barrier3d[0].TMAX)
        else:
            tmax_array_100.append(outwash100_obj.barrier3d[0].TMAX)

        # average net migration
        end_net_migration_100.append(shoreline_pos_array_100[storm_num-1][-1])
        # average net migration excluding drownings
        if drowning_array_100[storm_num - 1] == 0:
            end_net_migration_no_drowns_100.append(shoreline_pos_array_100[storm_num-1][-1])
        else:  # average net migration only including drownings
            end_net_migration_only_drowns_100.append(shoreline_pos_array_100[storm_num - 1][-1])
        # shoreface slope
        sfTS100 = outwash100_obj.barrier3d[0].s_sf_TS
        shoreface_slope_array_100[storm_num - 1] = sfTS100
        # average elevation
        hbTS100 = np.array(outwash100_obj.barrier3d[0].h_b_TS) * 10  # m MHW
        hbTS100 = hbTS100[~np.isnan(hbTS100)]  # hbTS is the average elevation at each model time step
        int_height_array_100[storm_num - 1] = hbTS100  # this is an array with array components
        avg_int_height_array_100[storm_num - 1] = np.mean(hbTS100)
        # width
        avg_int_width_TS100 = np.array(outwash100_obj.barrier3d[0].InteriorWidth_AvgTS) * 10  # meters
        int_width_array_100[storm_num - 1] = avg_int_width_TS100
        avg_int_width_array_100[storm_num - 1] = np.mean(avg_int_width_TS100)
        # dune crest
        sub_domain100 = outwash100_obj.barrier3d[0]._DuneDomain[0:tmax_array_100[storm_num-1], :, :]
        dune_crest_array_100[storm_num - 1] = sub_domain100.max(axis=2)
        avg_dune_crest_array_100.append(np.average(dune_crest_array_100[storm_num - 1], axis=0))
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
        # shoreline position
        m_xsTS50 = np.subtract(outwash50_obj.barrier3d[0].x_s_TS, outwash50_obj.barrier3d[0].x_s_TS[0])
        m_xsTS50 = np.multiply(m_xsTS50, 10)
        m_xsTS50 = m_xsTS50[~np.isnan(m_xsTS50)]
        shoreline_pos_array_50[storm_num-1] = m_xsTS50
        tmax_first_drown_50 = len(m_xsTS50)  # based on the first drowning event if there is one
        if tmax_first_drown_50 < 101:
            drowning_array_50[storm_num - 1] = 1
            if tmax_first_drown_50 < outwash50_obj.barrier3d[0].TMAX:
                drown_year_array_50.append(tmax_first_drown_50)
                tmax_array_50.append(tmax_first_drown_50)
            else:
                drown_year_array_50.append(outwash50_obj.barrier3d[0].TMAX)
                tmax_array_50.append(outwash50_obj.barrier3d[0].TMAX)
        else:
            tmax_array_50.append(outwash50_obj.barrier3d[0].TMAX)
        # average net migration
        end_net_migration_50.append(shoreline_pos_array_50[storm_num-1][-1])
        # average net migration excluding drownings
        if drowning_array_50[storm_num - 1] == 0:
            end_net_migration_no_drowns_50.append(shoreline_pos_array_50[storm_num-1][-1])
        else:  # average net migration only including drownings
            end_net_migration_only_drowns_50.append(shoreline_pos_array_50[storm_num - 1][-1])
        # shoreface slope
        sfTS50 = outwash50_obj.barrier3d[0].s_sf_TS
        shoreface_slope_array_50[storm_num - 1] = sfTS50
        # average elevation
        hbTS50 = np.array(outwash50_obj.barrier3d[0].h_b_TS) * 10  # m MHW
        hbTS50 = hbTS50[~np.isnan(hbTS50)]
        int_height_array_50[storm_num - 1] = hbTS50
        avg_int_height_array_50[storm_num - 1] = np.mean(hbTS50)
        # dune crest
        sub_domain50 = outwash50_obj.barrier3d[0]._DuneDomain[0:tmax_array_50[storm_num - 1], :, :]
        dune_crest_array_50[storm_num - 1] = sub_domain50.max(axis=2)
        avg_dune_crest_array_50.append(np.average(dune_crest_array_50[storm_num - 1], axis=0))
        # width
        avg_int_width_TS50 = np.array(outwash50_obj.barrier3d[0].InteriorWidth_AvgTS) * 10  # meters
        int_width_array_50[storm_num - 1] = avg_int_width_TS50
        avg_int_width_array_50[storm_num - 1] = np.mean(avg_int_width_TS50)
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
        # shoreline position
        m_xsTS0 = np.subtract(outwash0_obj.barrier3d[0].x_s_TS, outwash0_obj.barrier3d[0].x_s_TS[0])
        m_xsTS0 = np.multiply(m_xsTS0, 10)
        m_xsTS0 = m_xsTS0[~np.isnan(m_xsTS0)]
        shoreline_pos_array_0[storm_num-1] = m_xsTS0
        tmax_first_drown_0 = len(m_xsTS0)  # based on the first drowning event if there is one
        if tmax_first_drown_0 < 101:
            drowning_array_0[storm_num - 1] = 1
            if tmax_first_drown_0 < outwash0_obj.barrier3d[0].TMAX:
                drown_year_array_0.append(tmax_first_drown_0)
                tmax_array_0.append(tmax_first_drown_0)
            else:
                drown_year_array_0.append(outwash0_obj.barrier3d[0].TMAX)
                tmax_array_0.append(outwash0_obj.barrier3d[0].TMAX)
        else:
            tmax_array_0.append(outwash0_obj.barrier3d[0].TMAX)
        # average net migration
        end_net_migration_0.append(shoreline_pos_array_0[storm_num-1][-1])
        # average net migration excluding drownings
        if drowning_array_0[storm_num - 1] == 0:
            end_net_migration_no_drowns_0.append(shoreline_pos_array_0[storm_num - 1][-1])
        else:  # average net migration only including drownings
            end_net_migration_only_drowns_0.append(shoreline_pos_array_0[storm_num - 1][-1])
        # shoreface slope
        sfTS0 = outwash0_obj.barrier3d[0].s_sf_TS
        shoreface_slope_array_0[storm_num - 1] = sfTS0
        # average elevation
        hbTS0 = np.array(outwash0_obj.barrier3d[0].h_b_TS) * 10  # meters MHW
        hbTS0 = hbTS0[~np.isnan(hbTS0)]
        int_height_array_0[storm_num - 1] = hbTS0
        avg_int_height_array_0[storm_num - 1] = np.mean(hbTS0)
        # dune crest
        sub_domain0 = outwash0_obj.barrier3d[0]._DuneDomain[0:tmax_array_0[storm_num - 1], :, :]
        dune_crest_array_0[storm_num - 1] = sub_domain0.max(axis=2)
        avg_dune_crest_array_0.append(np.average(dune_crest_array_0[storm_num - 1], axis=0))
        # width
        avg_int_width_TS0 = np.array(outwash0_obj.barrier3d[0].InteriorWidth_AvgTS) * 10  # meters
        int_width_array_0[storm_num - 1] = avg_int_width_TS0
        avg_int_width_array_0[storm_num - 1] = np.mean(avg_int_width_TS0)
        # overwash
        QowTS0 = outwash0_obj.barrier3d[0].QowTS
        overwash_array_0[storm_num - 1] = QowTS0
        # outwash
        QoutTS0 = outwash0_obj.outwash[0]._outwash_flux_TS
        outwash_array_0[storm_num - 1] = QoutTS0

    DuneCrest_b3d = np.average(np.vstack(avg_dune_crest_array_b3d).astype(float), axis=0) * 10
    print("average dune crest height overwash only: {0}".format(np.average(DuneCrest_b3d)))
    DuneCrest_100 = np.average(np.vstack(avg_dune_crest_array_100).astype(float), axis=0) * 10
    print("average dune crest height 100% outwash: {0}".format(np.average(DuneCrest_100)))
    DuneCrest_50 = np.average(np.vstack(avg_dune_crest_array_50).astype(float), axis=0) * 10
    print("average dune crest height 50% outwash: {0}".format(np.average(DuneCrest_50)))
    DuneCrest_0 = np.average(np.vstack(avg_dune_crest_array_0).astype(float), axis=0) * 10
    print("average dune crest height 0% outwash: {0}".format(np.average(DuneCrest_0)))

    fig1 = plt.figure()
    plt.title("{0}".format(rname))
    plt.plot(DuneCrest_b3d, label="overwash only")
    plt.plot(DuneCrest_100, label="100% outwash")
    plt.plot(DuneCrest_50, label="50% outwash")
    plt.plot(DuneCrest_0, label="0% outwash")
    plt.ylim(bottom=0, top=5)
    plt.xlabel("Alongshore Barrier Length (m)")
    plt.ylabel("Average Dune Crest Elevation (m)")
    xtick_max = np.shape(DuneCrest_0)[0]
    x_ticks = np.array(range(0, xtick_max, 10))
    x_tick_labels = x_ticks * 10
    plt.xticks(x_ticks, x_tick_labels)
    plt.legend()

    # because we have 100 storms and drowns = 1, the sum of the array is the percent that drown
    percent_drown_b3d = np.sum(drowning_array_b3d)
    print("{0}% of barriers drown for the overwash only scenario".format(int(percent_drown_b3d)))
    percent_drown_100 = np.sum(drowning_array_100)
    print("{0}% of barriers drown for the 100% outwash to shoreface scenario".format(int(percent_drown_100)))
    percent_drown_50 = np.sum(drowning_array_50)
    print("{0}% of barriers drown for the 50% outwash to shoreface scenario".format(int(percent_drown_50)))
    percent_drown_0 = np.sum(drowning_array_0)
    print("{0}% of barriers drown for the 0% outwash to shoreface scenario".format(int(percent_drown_0)))

    # average drown years
    if len(drown_year_array_b3d) > 0:
        avg_drown_year_b3d = np.average(drown_year_array_b3d)
        min_drown_year_b3d = np.min(drown_year_array_b3d)
        max_drown_year_b3d = np.max(drown_year_array_b3d)
        print("overwash only, n = {3} \n average drown year: {0} \n min drown year is: {1} \n max drown year is: {2}"
              "".format(np.round(avg_drown_year_b3d, 0), np.round(min_drown_year_b3d, 0),
                        np.round(max_drown_year_b3d, 0), len(drown_year_array_b3d)))
    if len(drown_year_array_100) > 0:
        avg_drown_year_100 = np.average(drown_year_array_100)
        min_drown_year_100 = np.min(drown_year_array_100)
        max_drown_year_100 = np.max(drown_year_array_100)
        print("100% outwash, n = {3} \n average drown year: {0} \n min drown year is: {1} \n max drown year is: {2}"
              "".format(np.round(avg_drown_year_100, 0), np.round(min_drown_year_100, 0),
                        np.round(max_drown_year_100, 0), len(drown_year_array_100)))
    if len(drown_year_array_50) > 0:
        avg_drown_year_50 = np.average(drown_year_array_50)
        min_drown_year_50 = np.min(drown_year_array_50)
        max_drown_year_50 = np.max(drown_year_array_50)
        print(" 50% outwash, n = {3} \n average drown year: {0} \n min drown year is: {1} \n max drown year is: {2}"
              "".format(np.round(avg_drown_year_50, 0), np.round(min_drown_year_50, 0), np.round(max_drown_year_50, 0),
                        len(drown_year_array_50)))
    if len(drown_year_array_0) > 0:
        avg_drown_year_0 = np.average(drown_year_array_0)
        min_drown_year_0 = np.min(drown_year_array_0)
        max_drown_year_0 = np.max(drown_year_array_0)
        print("0% outwash, n = {3} \n average drown year: {0} \n min drown year is: {1} \n max drown year is: {2}"
              "".format(np.round(avg_drown_year_0, 0), np.round(min_drown_year_0, 0), np.round(max_drown_year_0, 0),
                        len(drown_year_array_0)))

    total_drown_years = []
    total_drown_years = drown_year_array_b3d + drown_year_array_100 + drown_year_array_50 + drown_year_array_0
    total_drown_years = np.array(total_drown_years)
    unique_drown_years = np.unique(total_drown_years)
    bar_b3d = []
    bar_100 = []
    bar_50 = []
    bar_0 = []

    for year in unique_drown_years:
        count_b3d = 0
        count_100 = 0
        count_50 = 0
        count_0 = 0
        for a in range(len(drown_year_array_b3d)):
            if drown_year_array_b3d[a] == year:
                count_b3d += 1
        bar_b3d.append(count_b3d)
        for a in range(len(drown_year_array_100)):
            if drown_year_array_100[a] == year:
                count_100 += 1
        bar_100.append(count_100)
        for a in range(len(drown_year_array_50)):
            if drown_year_array_50[a] == year:
                count_50 += 1
        bar_50.append(count_50)
        for a in range(len(drown_year_array_0)):
            if drown_year_array_0[a] == year:
                count_0 += 1
        bar_0.append(count_0)

    bar_b3d = np.array(bar_b3d)
    bar_100 = np.array(bar_100)
    bar_50 = np.array(bar_50)
    bar_0 = np.array(bar_0)

    # change to a string so that the x axis only plots these years rather than all years between the min and max values
    unique_drown_years_strings = []
    for year in unique_drown_years:
        unique_drown_years_strings.append(str(year))

    # printing migration stats
    if migration_stats:
        # print avg net migration stats
        print("avg outwash stats")
        print("overwash only \n avg net migration: {0} \n avg net migration excluding drown years: {1}"
              .format(np.round(np.average(end_net_migration_b3d)), np.round(np.average(end_net_migration_no_drowns_b3d))))
        if len(end_net_migration_only_drowns_b3d) > 0:
            print("avg net migration of barriers that drown: {0}".format(np.average(end_net_migration_only_drowns_b3d)))
        print("100% outwash \n avg net migration: {0} \n avg net migration excluding drown years: {1} "
              "\n avg net migration of barriers that drown: {2}".format(np.round(np.average(end_net_migration_100)),
                                                                    np.round(np.average(end_net_migration_no_drowns_100)),
                                                                    np.round(np.average(end_net_migration_only_drowns_100))))
        print("50% outwash \n avg net migration: {0} \n avg net migration excluding drown years: {1} "
              "\n avg net migration of barriers that drown: {2}".format(np.round(np.average(end_net_migration_50)),
                                                                    np.round(np.average(end_net_migration_no_drowns_50)),
                                                                    np.round(np.average(end_net_migration_only_drowns_50))))
        print("0% outwash \n avg net migration: {0} \n avg net migration excluding drown years: {1} "
              "\n avg net migration of barriers that drown: {2}".format(np.round(np.average(end_net_migration_0)),
                                                                    np.round(np.average(end_net_migration_no_drowns_0)),
                                                                    np.round(np.average(end_net_migration_only_drowns_0))))

        # print max net migration stats
        print("max outwash stats")
        print("overwash only \n max net migration: {0} \n max net migration excluding drown years: {1}"
              .format(np.round(np.max(end_net_migration_b3d)), np.round(np.max(end_net_migration_no_drowns_b3d))))
        if len(end_net_migration_only_drowns_b3d) > 0:
            print("max net migration of barriers that drown: {0}".format(np.max(end_net_migration_only_drowns_b3d)))
        print("100% outwash \n max net migration: {0} \n max net migration excluding drown years: {1} "
              "\n max net migration of barriers that drown: {2}".format(np.round(np.max(end_net_migration_100)),
                                                                    np.round(np.max(end_net_migration_no_drowns_100)),
                                                                    np.round(np.max(end_net_migration_only_drowns_100))))
        print("50% outwash \n max net migration: {0} \n max net migration excluding drown years: {1} "
              "\n max net migration of barriers that drown: {2}".format(np.round(np.max(end_net_migration_50)),
                                                                    np.round(np.max(end_net_migration_no_drowns_50)),
                                                                    np.round(np.max(end_net_migration_only_drowns_50))))
        print("0% outwash \n max net migration: {0} \n max net migration excluding drown years: {1} "
              "\n max net migration of barriers that drown: {2}".format(np.round(np.max(end_net_migration_0)),
                                                                    np.round(np.max(end_net_migration_no_drowns_0)),
                                                                    np.round(np.max(end_net_migration_only_drowns_0))))
        # print min net migration stats
        print("min outwash stats")
        print("overwash only \n min net migration: {0} \n min net migration excluding drown years: {1}"
              .format(np.round(np.min(end_net_migration_b3d)), np.round(np.min(end_net_migration_no_drowns_b3d))))
        if len(end_net_migration_only_drowns_b3d) > 0:
            print("min net migration of barriers that drown: {0}".format(np.min(end_net_migration_only_drowns_b3d)))
        print("100% outwash \n min net migration: {0} \n min net migration excluding drown years: {1} "
              "\n min net migration of barriers that drown: {2}".format(np.round(np.min(end_net_migration_100)),
                                                                    np.round(np.min(end_net_migration_no_drowns_100)),
                                                                    np.round(np.min(end_net_migration_only_drowns_100))))
        print("50% outwash \n min net migration: {0} \n min net migration excluding drown years: {1} "
              "\n min net migration of barriers that drown: {2}".format(np.round(np.min(end_net_migration_50)),
                                                                    np.round(np.min(end_net_migration_no_drowns_50)),
                                                                    np.round(np.min(end_net_migration_only_drowns_50))))
        print("0% outwash \n min net migration: {0} \n excluding drown years: {1} "
              "\n min net migration of barriers that drown: {2}".format(np.round(np.min(end_net_migration_0)),
                                                                    np.round(np.min(end_net_migration_no_drowns_0)),
                                                                    np.round(np.min(end_net_migration_only_drowns_0))))

    # printing geometry stats
    if geomoetry_stats:
        # print avg height and width stats
        print("avg geometry stats")
        print("overwash only \n avg interior height: {0} \n avg interior width: {1}"
              .format(np.round(np.average(avg_int_height_array_b3d), 3),
                      np.round(np.average(avg_int_width_array_b3d))))
        print("100% outwash \n avg interior height: {0} \n avg interior width: {1}"
              .format(np.round(np.average(avg_int_height_array_100), 3),
                      np.round(np.average(avg_int_width_array_100))))
        print("50% outwash \n avg interior height: {0} \n avg interior width: {1}"
              .format(np.round(np.average(avg_int_height_array_50), 3),
                      np.round(np.average(avg_int_width_array_50))))
        print("0% outwash \n avg interior height: {0} \n avg interior width: {1}"
              .format(np.round(np.average(avg_int_height_array_0), 3),
                      np.round(np.average(avg_int_width_array_0))))

    if plotters:
        # histogram of years that the barriers drown (three outwash scenarios only, each on separate plot)
        plt.rcParams.update({"font.size": 12})
        # bins = 100
        # fig1 = plt.figure()
        # fig1.suptitle('{0}'.format(rname), weight="bold")
        # ax1 = fig1.add_subplot(131)
        # ax1.hist(drown_year_array_100, bins=bins)
        # ax1.set_title("100% washout")
        # ax1.set_ylabel("number of barriers that drown")
        # ax1.set_xlabel("drown year")
        # plt.gca().xaxis.tick_bottom()
        # ax1.set_xlim(left=1, right=100)
        # ax1.set_ylim(bottom=0, top=40)
        #
        # ax1 = fig1.add_subplot(132)
        # ax1.hist(drown_year_array_50, bins=bins)
        # ax1.set_title("50% washout")
        # ax1.set_ylabel("frequency")
        # ax1.set_xlabel("drown year")
        # plt.gca().xaxis.tick_bottom()
        # ax1.set_xlim(left=1, right=100)
        # ax1.set_ylim(bottom=0, top=40)
        # ax1.set_ylabel(None)
        #
        # ax1 = fig1.add_subplot(133)
        # ax1.hist(drown_year_array_0, bins=bins)
        # ax1.set_title("0% washout")
        # ax1.set_ylabel("frequency")
        # ax1.set_xlabel("drown year")
        # plt.gca().xaxis.tick_bottom()
        # ax1.set_xlim(left=1, right=100)
        # ax1.set_ylim(bottom=0, top=40)
        # ax1.set_ylabel(None)

        # stacked bar charts showing same data as historams
        fig2, ax = plt.subplots()
        ax.bar(unique_drown_years_strings, bar_b3d, label="overwash only")
        ax.bar(unique_drown_years_strings, bar_100, label="100% outwash", bottom=bar_b3d)
        ax.bar(unique_drown_years_strings, bar_50, label="50% outwash", bottom=bar_b3d+bar_100)
        ax.bar(unique_drown_years_strings, bar_0, label="0% outwash", bottom=bar_b3d+bar_100+bar_50)
        plt.xlabel("Drown Year")
        plt.ylabel("Number of Barriers that Drown")
        plt.legend()
        plt.ylim(top=65)
        plt.title('{0}'.format(rname), weight="bold")

        # plotting shoreline position for all 100 storms, each scenario on separate plot
        # shoreline position: finding the extremes
        min_b3d_migration = min(end_net_migration_b3d)
        min_b3d_migration_index = end_net_migration_b3d.index(min_b3d_migration)
        min_b3d_shoreline_position_array = shoreline_pos_array_b3d[min_b3d_migration_index]
        max_b3d_migration = max(end_net_migration_b3d)
        max_b3d_migration_index = end_net_migration_b3d.index(max_b3d_migration)
        max_b3d_shoreline_position_array = shoreline_pos_array_b3d[max_b3d_migration_index]

        # for the outwash scenarios, we want to plot the min and maxes for those storms that do not drown
        # as well as some storms that drown
        min_100_migration = min(end_net_migration_no_drowns_100)
        min_100_migration_index = end_net_migration_100.index(min_100_migration)
        min_100_shoreline_position_array = shoreline_pos_array_100[min_100_migration_index]
        max_100_migration = max(end_net_migration_no_drowns_100)
        max_100_migration_index = end_net_migration_100.index(max_100_migration)
        max_100_shoreline_position_array = shoreline_pos_array_100[max_100_migration_index]

        min_50_migration = min(end_net_migration_no_drowns_50)
        min_50_migration_index = end_net_migration_50.index(min_50_migration)
        min_50_shoreline_position_array = shoreline_pos_array_50[min_50_migration_index]
        max_50_migration = max(end_net_migration_no_drowns_50)
        max_50_migration_index = end_net_migration_50.index(max_50_migration)
        max_50_shoreline_position_array = shoreline_pos_array_50[max_50_migration_index]

        min_0_migration = min(end_net_migration_no_drowns_0)
        min_0_migration_index = end_net_migration_0.index(min_0_migration)
        min_0_shoreline_position_array = shoreline_pos_array_0[min_0_migration_index]
        max_0_migration = max(end_net_migration_no_drowns_0)
        max_0_migration_index = end_net_migration_0.index(max_0_migration)
        max_0_shoreline_position_array = shoreline_pos_array_0[max_0_migration_index]

        bot = -225
        top = 425
        fig3 = plt.figure()
        fig3.suptitle('{0}'.format(rname), weight="bold")
        ax1 = fig3.add_subplot(221)
        ax1.set_title("overwash only")
        ax1.set_ylim(bottom=bot, top=top)
        ax2 = fig3.add_subplot(222)
        ax2.set_title("100% outwash")
        ax2.set_ylim(bottom=bot, top=top)
        ax3 = fig3.add_subplot(223)
        ax3.set_title("50% outwash")
        ax3.set_ylim(bottom=bot, top=top)
        ax4 = fig3.add_subplot(224)
        ax4.set_title("0% outwash")
        ax4.set_ylim(bottom=bot, top=top)

        fig3.text(0.5, 0.04, 'Time (years)', ha='center')
        fig3.text(0.06, 0.5, 'Shoreline Position (m)', va='center', rotation='vertical')
        fig3.subplots_adjust(hspace=0.4)

        alpha = 0.1

        for x in range(100):
            ax1.plot(shoreline_pos_array_b3d[x], alpha=alpha)
            ax2.plot(shoreline_pos_array_100[x], alpha=alpha)
            ax3.plot(shoreline_pos_array_50[x], alpha=alpha)
            ax4.plot(shoreline_pos_array_0[x], alpha=alpha)

        min_color = "black"
        max_color = "red"

        ax1.plot(min_b3d_shoreline_position_array, label="min migration", color=min_color, zorder=0)
        ax1.plot(max_b3d_shoreline_position_array, label="max migration", color=max_color, zorder=0)

        ax2.plot(min_100_shoreline_position_array, label="min migration", color=min_color, zorder=0)
        ax2.plot(max_100_shoreline_position_array, label="max migration", color=max_color, zorder=0)

        ax3.plot(min_50_shoreline_position_array, label="min migration", color=min_color, zorder=0)
        ax3.plot(max_50_shoreline_position_array, label="max migration", color=max_color, zorder=0)

        ax4.plot(min_0_shoreline_position_array, label="min migration", color=min_color, zorder=0)
        ax4.plot(max_0_shoreline_position_array, label="max migration", color=max_color, zorder=0)

        plot_drown_year_array = [21, 41, 61, 81]
        colors = ["lime", "fuchsia", "blue", "blueviolet"]
        zorder_list = [3, 2, 1, 0]
        for plot_drown_year in plot_drown_year_array:
            color_index = plot_drown_year_array.index(plot_drown_year)
            color = colors[color_index]
            zorders = zorder_list[plot_drown_year_array.index(plot_drown_year)]
            if plot_drown_year in tmax_array_100:
                index100 = tmax_array_100.index(plot_drown_year)
                ax2.plot(shoreline_pos_array_100[index100], zorder=zorders, label="drown year: {0}".format(plot_drown_year), color=color)
            if plot_drown_year in tmax_array_50:
                index50 = tmax_array_50.index(plot_drown_year)
                ax3.plot(shoreline_pos_array_50[index50], zorder=zorders, label="drown year: {0}".format(plot_drown_year), color=color)
            if plot_drown_year in tmax_array_0:
                index0 = tmax_array_0.index(plot_drown_year)
                ax4.plot(shoreline_pos_array_0[index0], zorder=zorders, label="drown year: {0}".format(plot_drown_year), color=color)

        ax1.legend(fontsize=8, ncol=2, loc="upper center")
        ax2.legend(fontsize=8, ncol=3, loc="upper center")
        ax3.legend(fontsize=8, ncol=3, loc="upper center")
        ax4.legend(fontsize=8, ncol=3, loc="lower center")

        # plotting first storm shoreline position, shoreface slope, avg int. height and width,
        # overwash and outwash flux, all four scenarios on same plot

        storm_num = 1
        fig8 = plt.figure()
        fig8.suptitle('{0}'.format(rname), weight="bold")
        ax1 = fig8.add_subplot(231)
        ls = "dashed"
        ax1.plot(shoreline_pos_array_b3d[storm_num-1])
        ax1.plot(shoreline_pos_array_100[storm_num-1], linestyle=ls)
        ax1.plot(shoreline_pos_array_50[storm_num-1], linestyle=ls)
        ax1.plot(shoreline_pos_array_0[storm_num-1], linestyle=ls)
        ax1.legend(["no outwash", "100%", "50%", "0%"], prop={'size': 9})
        ax1.set_ylabel("Shoreline Position (m)")
        ax1.set_xlabel("Simulation Years")
        ax1.set_ylim(bottom=-60, top=160)

        ax2 = fig8.add_subplot(232)
        ax2.plot(shoreface_slope_array_b3d[storm_num-1])
        ax2.plot(shoreface_slope_array_100[storm_num-1], linestyle=ls)
        ax2.plot(shoreface_slope_array_50[storm_num-1], linestyle=ls)
        ax2.plot(shoreface_slope_array_0[storm_num-1], linestyle=ls)
        ax2.set_ylabel("Shoreface Slope")
        ax2.set_xlabel("Simulation Years")
        ax2.set_ylim(top=0.02)

        ax3 = fig8.add_subplot(235)
        ax3.plot(int_height_array_b3d[storm_num-1])
        ax3.plot(int_height_array_100[storm_num-1], linestyle=ls)
        ax3.plot(int_height_array_50[storm_num-1], linestyle=ls)
        ax3.plot(int_height_array_0[storm_num-1], linestyle=ls)
        ax3.set_ylabel("Average Interior Elevation (m MHW)")
        ax3.set_xlabel("Simulation Years")
        ax3.set_ylim(top=1.6)

        ax4 = fig8.add_subplot(234)
        ax4.plot(int_width_array_b3d[storm_num-1])
        ax4.plot(int_width_array_100[storm_num-1], linestyle=ls)
        ax4.plot(int_width_array_50[storm_num-1], linestyle=ls)
        ax4.plot(int_width_array_0[storm_num-1], linestyle=ls)
        ax4.set_ylabel("Average Interior Width (m)")
        ax4.set_xlabel("Simulation Years")
        ax4.set_ylim(top=310)

        ax5 = fig8.add_subplot(233)
        ax5.plot(overwash_array_b3d[storm_num-1])
        ax5.plot(overwash_array_100[storm_num-1], linestyle=ls)
        ax5.plot(overwash_array_50[storm_num-1], linestyle=ls)
        ax5.plot(overwash_array_0[storm_num-1], linestyle=ls)
        ax5.set_ylabel("Overwash Flux (m3/m)")
        ax5.set_xlabel("Simulation Years")
        ax5.set_ylim(top=150)

        ax6 = fig8.add_subplot(236)
        ax6.plot(outwash_array_b3d[storm_num-1])
        ax6.plot(outwash_array_100[storm_num-1], linestyle=ls)
        ax6.plot(outwash_array_50[storm_num-1], linestyle=ls)
        ax6.plot(outwash_array_0[storm_num-1], linestyle=ls)
        ax6.set_ylabel("Outwash Flux (m3/m)")
        ax6.set_xlabel("Simulation Years")
        ax6.set_ylim(top=1200)

        fig8.subplots_adjust(hspace=0.3, wspace=0.3)
