# loading in the npz files and saving variables that we want to plot together

# Lexi Van Blunk
# 5/15/2025
# calculating the average interior elevation using DomainTS instead of the h_b_TS variable because I think it is
# more clear

import numpy as np
from matplotlib import pyplot as plt

# ---------------------------------- set model parameters that change per run ------------------------------------------
rname_array = ["r025", "r035"]
for rname in rname_array:
    storm_interval = 20  # 20 or 10 years
    config = 4  # 1, 2, 3, or 4

    # Display stats on console/show plots
    migration_stats = False
    plotters = False
    geomoetry_stats = True
    dune_crest_stats = False
    flux_stats = False

    # location of the npz files
    datadir_b3d = "C:/Users/Lexi/PycharmProjects/CASCADE/data/outwash_data/storms/slope0pt03/rerun_output/{}/overwash_only/".format(
        rname
    )
    datadir_100 = "C:/Users/Lexi/PycharmProjects/CASCADE/data/outwash_data/storms/slope0pt03/rerun_output/{}/outwash100/".format(
        rname
    )
    datadir_50 = "C:/Users/Lexi/PycharmProjects/CASCADE/data/outwash_data/storms/slope0pt03/rerun_output/{}/outwash50/".format(
        rname
    )
    datadir_0 = "C:/Users/Lexi/PycharmProjects/CASCADE/data/outwash_data/storms/slope0pt03/rerun_output/{}/outwash0/".format(
        rname
    )

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

    no_drowns_last_elev_b3d_array = np.zeros(100, dtype=object)
    no_drowns_last_elev_100_array = np.zeros(100, dtype=object)
    no_drowns_last_elev_50_array = np.zeros(100, dtype=object)
    no_drowns_last_elev_0_array = np.zeros(100, dtype=object)

    no_drowns_last_width_b3d = []
    no_drowns_last_width_100 = []
    no_drowns_last_width_50 = []
    no_drowns_last_width_0 = []

    for storm_num in range(1, 101):
        # print(storm_num)
        # b3d variables
        filename_b3d = "config{}_b3d_startyr1_interval{}yrs_Slope0pt03_{}.npz".format(
            config, storm_interval, storm_num
        )
        file_b3d = datadir_b3d + filename_b3d
        b3d = np.load(file_b3d, allow_pickle=True)
        b3d_obj = b3d["cascade"][0]
        tmax_b3d = b3d_obj.barrier3d[0].TMAX
        tmax_array_b3d.append(tmax_b3d)

        for t in range(0, tmax_b3d):
            # average elevation
            domain_TS = (
                b3d_obj.barrier3d[0].DomainTS[t] * 10
            )  # domain for the current time step, t
            avg_elev_TS = np.average(
                domain_TS
            )  # average value of the domain for this time step
            avg_elev_TS_array_b3d.append(
                avg_elev_TS
            )  # array of the average elevations at each time step

        last_domain_TS_avg = np.average(
            domain_TS
        )  # the current value of domain_TS is for the tmax
        avg_elev_b3d = np.average(avg_elev_TS_array_b3d)
        avg_elev_b3d_array[storm_num - 1] = avg_elev_b3d
        avg_last_elev_b3d_array[storm_num - 1] = last_domain_TS_avg

        if tmax_b3d == 101:
            # interior domain: remove a full row if all values are less than 0
            check = 1
            while check == 1:
                if all(x <= 0 for x in domain_TS[-1, :]):  # bay_depth = -3 m
                    domain_TS = np.delete(domain_TS, -1, axis=0)
                else:
                    check = 0

            domain_TS[domain_TS < 0] = 0  # set all values less than 0 to 0
            new_ave_interior_height = np.average(domain_TS)  # all in m MHW
            no_drowns_last_elev_b3d_array[storm_num - 1] = new_ave_interior_height

            # width
            int_width_TS = (
                np.array(b3d_obj.barrier3d[0].InteriorWidth_AvgTS) * 10
            )  # meters
            last_int_width = int_width_TS[tmax_b3d - 1]
            no_drowns_last_width_b3d.append(last_int_width)

        # 100% variables
        filename_100 = (
            "config{}_outwash100_startyr1_interval{}yrs_Slope0pt03_{}.npz".format(
                config, storm_interval, storm_num
            )
        )
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

        last_domain_TS_avg = np.average(
            domain_TS
        )  # domain TS will be the domain before drowning
        avg_elev_100 = np.average(avg_elev_TS_array_100)
        avg_elev_100_array[storm_num - 1] = avg_elev_100
        avg_last_elev_100_array[storm_num - 1] = last_domain_TS_avg
        if tmax_100 == 101:
            # interior domain: remove a full row if all values are less than 0
            check = 1
            while check == 1:
                if all(x <= 0 for x in domain_TS[-1, :]):  # bay_depth = -3 m
                    domain_TS = np.delete(domain_TS, -1, axis=0)
                else:
                    check = 0

            domain_TS[domain_TS < 0] = 0  # set all values less than 0 to 0
            new_ave_interior_height = np.average(domain_TS)  # all in m MHW
            no_drowns_last_elev_100_array[storm_num - 1] = new_ave_interior_height

            # width
            int_width_TS = (
                np.array(outwash100_obj.barrier3d[0].InteriorWidth_AvgTS) * 10
            )  # meters
            last_int_width = int_width_TS[tmax_100 - 1]
            no_drowns_last_width_100.append(last_int_width)

        # 50% variables
        filename_50 = (
            "config{}_outwash50_startyr1_interval{}yrs_Slope0pt03_{}.npz".format(
                config, storm_interval, storm_num
            )
        )
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
        avg_elev_50_array[storm_num - 1] = avg_elev_50
        avg_last_elev_50_array[storm_num - 1] = last_domain_TS_avg
        if tmax_50 == 101:
            # interior domain: remove a full row if all values are less than 0
            check = 1
            while check == 1:
                if all(x <= 0 for x in domain_TS[-1, :]):  # bay_depth = -3 m
                    domain_TS = np.delete(domain_TS, -1, axis=0)
                else:
                    check = 0

            domain_TS[domain_TS < 0] = 0  # set all values less than 0 to 0
            new_ave_interior_height = np.average(domain_TS)  # all in m MHW
            no_drowns_last_elev_50_array[storm_num - 1] = new_ave_interior_height

            # width
            int_width_TS = (
                np.array(outwash50_obj.barrier3d[0].InteriorWidth_AvgTS) * 10
            )  # meters
            last_int_width = int_width_TS[tmax_50 - 1]
            no_drowns_last_width_50.append(last_int_width)

        # 0% variables
        filename_0 = (
            "config{}_outwash0_startyr1_interval{}yrs_Slope0pt03_{}.npz".format(
                config, storm_interval, storm_num
            )
        )
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
        avg_elev_0_array[storm_num - 1] = avg_elev_0
        avg_last_elev_0_array[storm_num - 1] = last_domain_TS_avg

        if tmax_0 == 101:
            # interior domain: remove a full row if all values are less than 0
            check = 1
            while check == 1:
                if all(x <= 0 for x in domain_TS[-1, :]):  # bay_depth = -3 m
                    domain_TS = np.delete(domain_TS, -1, axis=0)
                else:
                    check = 0

            new_ave_interior_height = np.average(domain_TS)  # all in m MHW
            no_drowns_last_elev_0_array[storm_num - 1] = new_ave_interior_height

            # width
            int_width_TS = (
                np.array(outwash0_obj.barrier3d[0].InteriorWidth_AvgTS) * 10
            )  # meters
            last_int_width = int_width_TS[-1]
            no_drowns_last_width_0.append(last_int_width)

    avg_no_drowns_b3d = np.average(
        no_drowns_last_elev_b3d_array[np.nonzero(no_drowns_last_elev_b3d_array)]
    )
    avg_no_drowns_100 = np.average(
        no_drowns_last_elev_100_array[np.nonzero(no_drowns_last_elev_100_array)]
    )
    avg_no_drowns_50 = np.average(
        no_drowns_last_elev_50_array[np.nonzero(no_drowns_last_elev_50_array)]
    )
    avg_no_drowns_0 = np.average(
        no_drowns_last_elev_0_array[np.nonzero(no_drowns_last_elev_0_array)]
    )

    avg_no_drowns_last_width_b3d = np.average(no_drowns_last_width_b3d)
    avg_no_drowns_last_width_100 = np.average(no_drowns_last_width_100)
    avg_no_drowns_last_width_50 = np.average(no_drowns_last_width_50)
    avg_no_drowns_last_width_0 = np.average(no_drowns_last_width_0)

    # printing geometry stats
    if geomoetry_stats:
        # print non-drowning values
        print(f"avg geometry stats, {rname}")
        print(
            f"baseline \n avg interior height: {np.round(avg_no_drowns_b3d, 2)}"
        )
        print(
            "100% outwash \n avg interior height: {}".format(
                np.round(avg_no_drowns_100, 2)
            )
        )
        print(
            "50% outwash \n avg interior height: {}".format(
                np.round(avg_no_drowns_50, 2)
            )
        )
        print(
            f"0% outwash \n avg interior height: {np.round(avg_no_drowns_0, 2)}"
        )
        # width
        print(f"avg geometry stats, {rname}")
        print(
            "baseline \n avg interior width: {}".format(
                np.round(avg_no_drowns_last_width_b3d, 2)
            )
        )
        print(
            "100% outwash \n avg interior width: {}".format(
                np.round(avg_no_drowns_last_width_100, 2)
            )
        )
        print(
            "50% outwash \n avg interior width: {}".format(
                np.round(avg_no_drowns_last_width_50, 2)
            )
        )
        print(
            "0% outwash \n avg interior width: {}".format(
                np.round(avg_no_drowns_last_width_0, 2)
            )
        )
