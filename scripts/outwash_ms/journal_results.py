# Lexi Van Blunk
# 11/22/2025
# compiling analysis specifically used for the outwasher journal article
# Barrier vulnerability following outwash: A balance of overwash and dune gap recovery

# --------------------------------- description of variables ----------------------------------------------------------
# there are eight main scenarios: overwash only, 100, 50, or 0% washout to shoreface and two dune growth rates each
# each scenario was run with 100 different storm sequences for a model duration of 100 years or until drowning
# a single storm sequence consists of an overwash storm series and potentially outwash storm series
# summary of each scenario:
# overwash only, r=0.25, 100 storm sequences
# overwash only, r=0.35, 100 storm sequences
# 100% WTS, r=0.25, 100 storm sequences
# 100% WTS, r=0.35, 100 storm sequences
# 50% WTS, r=0.25, 100 storm sequences
# 50% WTS, r=0.35, 100 storm sequences
# 0% WTS, r=0.25, 100 storm sequences
# 0% WTS, r=0.35, 100 storm sequences

# the variables we analyzed include: percent drowned, overwash flux, outwash flux, total length of dune gaps,
# final interior elevation, final interior width, final shoreline position
# the final interior elevation, final interior width, final shoreline position only use the storm sequences that
# do not result in drowning

# percent drowned: 1 value per storm sequence, (drowned = 1, no drown = 0). average over the storm sequences (1 average)

# overwash and outwash flux: there is a flux at each model year (can be 0), so there are 100 values per storm sequence.
# first, we average these values to get 1 for a storm sequence. then average over all storm sequences to get one value
# for a given scenario (2 averages)

# dune crest elevation: the DuneDomain variable returns all model years, and then we use the maximum from either row
# to get an array that has n rows = model years and n cols = alongshore length with the dune crest elevations for all
# years. we average this array to get a single dune crest elevation for each storm sequence (1st average).
# then, we average over all storm sequences to get one value for each scenario (2nd average)

# total length of dune gaps: similar to dune crest elevation. for each model year, we start with the Dune Crest value
# and count the number of cells beneath a threshold (basically the berm elevation). then, we multilpy by 10 to get the
# total length of the dune gaps (1 cell length = 1 dam or 10 m). this gives us one value per model year, so we
# average over all years to get one per storm sequence. the average over all storm sequences (2 averages)

# final interior elevation: only for barriers that do not drown. we use the domain from the last time step and remove
# any rows that are completely below 0. then, we set all values less than 0 to 0. then average all cells in the domain
# to get one value per storm sequence. lastly, average over all storm sequences (2 averages)

# final interior width: only for barriers that do not drown. Average width is calculated during the model runs using
# Barrier3D’s FindWidths function. This takes the average of each column (across-shore) width. The cross-shore width is
# calculated as the number of consecutive cells from the interior to the ocean that are above sea level (here -3 m).
# therefore we already have the average barrier width for all model years and just use the final timestep value to get
# a single value for each storm sequence (already averaged). then average over all storm sequences (2 averages)

# final shoreline position: only for barriers that do not drown. each model year has a shoreline position, so the last
# value is the position at the end of the model simulation for a single storm sequence. we average over all storm
# sequences to get the final value per scenario (1 average)

import numpy as np
from matplotlib import pyplot as plt

# ---------------------------------- set model parameters that change per run ------------------------------------------
rname_array = ["r025", "r035"]
for rname in rname_array:
    storm_interval = 20  # 20 or 10 years
    config = 4  # 1, 2, 3, or 4

    # Display stats on console/show plots
    drowning_stats = False
    migration_stats = False
    plotters = False
    geomoetry_stats = False
    dune_stats = True
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

    drowning_array_b3d = np.zeros(100)
    drowning_array_100 = np.zeros(100)
    drowning_array_50 = np.zeros(100)
    drowning_array_0 = np.zeros(100)

    end_net_migration_no_drowns_b3d = []
    end_net_migration_no_drowns_100 = []
    end_net_migration_no_drowns_50 = []
    end_net_migration_no_drowns_0 = []

    avg_int_height_array_b3d = np.zeros(100, dtype=object)
    avg_int_height_array_100 = np.zeros(100, dtype=object)
    avg_int_height_array_50 = np.zeros(100, dtype=object)
    avg_int_height_array_0 = np.zeros(100, dtype=object)

    avg_int_width_array_b3d = np.zeros(100, dtype=object)
    avg_int_width_array_100 = np.zeros(100, dtype=object)
    avg_int_width_array_50 = np.zeros(100, dtype=object)
    avg_int_width_array_0 = np.zeros(100, dtype=object)

    dune_crest_array_b3d = np.zeros(100, dtype=object)
    dune_crest_array_100 = np.zeros(100, dtype=object)
    dune_crest_array_50 = np.zeros(100, dtype=object)
    dune_crest_array_0 = np.zeros(100, dtype=object)

    avg_dune_cells_array_b3d = np.zeros(100)
    avg_dune_cells_array_100 = np.zeros(100)
    avg_dune_cells_array_50 = np.zeros(100)
    avg_dune_cells_array_0 = np.zeros(100)

    overwash_array_b3d = np.zeros(100, dtype=object)
    overwash_array_100 = np.zeros(100, dtype=object)
    overwash_array_50 = np.zeros(100, dtype=object)
    overwash_array_0 = np.zeros(100, dtype=object)

    avg_overwash_array_b3d = np.zeros(100, dtype=object)
    avg_overwash_array_100 = np.zeros(100, dtype=object)
    avg_overwash_array_50 = np.zeros(100, dtype=object)
    avg_overwash_array_0 = np.zeros(100, dtype=object)

    avg_outwash_array_b3d = np.zeros(100, dtype=object)
    avg_outwash_array_100 = np.zeros(100, dtype=object)
    avg_outwash_array_50 = np.zeros(100, dtype=object)
    avg_outwash_array_0 = np.empty(100, dtype=object)

    for storm_num in range(1, 101):

        # b3d variables ------------------------------------------------------------------------------------------------
        filename_b3d = "config{}_b3d_startyr1_interval{}yrs_Slope0pt03_{}.npz".format(
            config, storm_interval, storm_num
        )
        file_b3d = datadir_b3d + filename_b3d
        b3d = np.load(file_b3d, allow_pickle=True)
        b3d_obj = b3d["cascade"][0]
        tmax_b3d = b3d_obj.barrier3d[0].TMAX + 1
        tmax_array_b3d.append(tmax_b3d)

        # drowning variable
        drowning_array_b3d[storm_num - 1] = b3d_obj.barrier3d[0].drown_break

        # shoreline position - only for barriers that do not drown
        m_xsTS = np.subtract(
            b3d_obj.barrier3d[0].x_s_TS, b3d_obj.barrier3d[0].x_s_TS[0]
        )
        m_xsTS = np.multiply(m_xsTS, 10)
        if b3d_obj.barrier3d[0].drown_break == 0:
            end_net_migration_no_drowns_b3d.append(m_xsTS[-1])

        # average elevation - only for barriers that do not drown
        if b3d_obj.barrier3d[0].drown_break == 0:
            domain_TS = b3d_obj.barrier3d[0].DomainTS[-1] * 10  # m MHW
            check = 1
            # remove all alongshore rows that are completely submerged less than 0
            while check == 1:
                if all(x <= 0 for x in domain_TS[-1, :]):
                    domain_TS = np.delete(domain_TS, -1, axis=0)
                else:
                    check = 0
            domain_TS[domain_TS < 0] = 0  # set all values less than 0 to 0
            new_ave_interior_height = np.average(
                domain_TS
            )  # all in m MHW, average value of the domain for this storm sequence
            avg_int_height_array_b3d[storm_num - 1] = (
                new_ave_interior_height  # array of the average elevations for all storm sequences
            )

        # width - only for barriers that do not drown
        if b3d_obj.barrier3d[0].drown_break == 0:
            int_width_TS = (
                np.array(b3d_obj.barrier3d[0].InteriorWidth_AvgTS) * 10
            )  # meters, already averaged values
            last_int_width = int_width_TS[-1]  # the width at the final time step
            # avg_width_last_TS_b3d.append(last_int_width)
            avg_int_width_array_b3d[storm_num - 1] = int_width_TS[
                -1
            ]  # array of the average widths for all storm sequences

        # dune crest
        sub_domain = b3d_obj.barrier3d[0]._DuneDomain[
            0:tmax_b3d, :, :
        ]  # dune domains for current storm sequence, all model years
        dune_crest = sub_domain.max(
            axis=2
        )  # n rows are simulation years, n cols are alongshore length (dune crests)
        dune_crest_avg = (
            np.average(dune_crest) * 10
        )  # average of all years, converted to meters
        dune_crest_array_b3d[storm_num - 1] = (
            dune_crest_avg  # array of the average dune crests for all storm sequences
        )
        # avg_dune_crest_array_b3d.append(np.average(dune_crest_array_b3d[storm_num-1], axis=0))  # for dune crest plots

        # dune gaps
        initial_gap_height = np.min(
            sub_domain[0]
        )  # get the lowest value of the dunes at year 0 without the berm elev
        berm_el = b3d_obj.barrier3d[0].BermEl
        dune_crest_elev = dune_crest + berm_el  # dam
        # here, a dune gap is defined as any cell less than or equal to the initial gap height + berm elevation
        dune_gap_limit = (
            initial_gap_height + berm_el
        )  # the minimum elevation of a dune cell
        # if a barrier drowns, all dune cells are at the berm elevation, so here, if all the dune crest cells for a
        # single year are at the berm elevation, we change the elevation to -0.3 (water) for plotting purposes
        for row in range(np.shape(dune_crest_elev)[0]):
            if np.all(dune_crest_elev[row] == berm_el) == True:
                dune_crest_elev[row] = -0.3  # dam

        n_dune_gap_cells_array = np.zeros([np.shape(dune_crest_elev)[0], 1])
        # count the number of dune gap cells in each row (each row = a model year)
        for year in range(np.shape(dune_crest_elev)[0]):
            dune_gap_row = dune_crest_elev[year]
            n_dune_gap_cells = len(dune_gap_row[dune_gap_row <= dune_gap_limit])
            n_dune_gap_cells_array[year] = (
                n_dune_gap_cells  # save the number of dune gap cells per model year
            )
        avg_dune_cells = np.average(
            n_dune_gap_cells_array
        )  # average the number dune gaps over all model years (for a single storm series)
        avg_dune_cells_array_b3d[storm_num - 1] = (
            avg_dune_cells * 10
        )  # each cell is 1 dam, so the total number of cells
        # gives the total length of gaps, then multiply by 10 to convert to meters

        # overwash
        QowTS = b3d_obj.barrier3d[
            0
        ].QowTS  # m^3/m, all model years for this storm sequence
        Qow_avg = np.mean(
            QowTS
        )  # m^3/m/yr average overwash flux for this storm sequence
        # overwash_array_b3d[storm_num-1] = np.array(QowTS)
        avg_overwash_array_b3d[storm_num - 1] = (
            Qow_avg  # array of the average overwash fluxes
        )

        # 100% variables  ---------------------------------------------------------------------------------------------
        filename_100 = (
            "config{}_outwash100_startyr1_interval{}yrs_Slope0pt03_{}.npz".format(
                config, storm_interval, storm_num
            )
        )
        file_100 = datadir_100 + filename_100
        outwash100 = np.load(file_100, allow_pickle=True)
        outwash100_obj = outwash100["cascade"][0]
        tmax_100 = outwash100_obj.barrier3d[0].TMAX + 1

        # drowning variable
        drowning_array_100[storm_num - 1] = outwash100_obj.barrier3d[0].drown_break

        # shoreline position - only for barriers that do not drown
        m_xsTS = np.subtract(
            outwash100_obj.barrier3d[0].x_s_TS, outwash100_obj.barrier3d[0].x_s_TS[0]
        )
        m_xsTS = np.multiply(m_xsTS, 10)
        if outwash100_obj.barrier3d[0].drown_break == 0:
            end_net_migration_no_drowns_100.append(m_xsTS[-1])

        # average elevation - only for barriers that do not drown
        if outwash100_obj.barrier3d[0].drown_break == 0:
            domain_TS = outwash100_obj.barrier3d[0].DomainTS[-1] * 10  # m MHW
            check = 1
            # remove all alongshore rows that are completely submerged less than 0
            while check == 1:
                if all(x <= 0 for x in domain_TS[-1, :]):
                    domain_TS = np.delete(domain_TS, -1, axis=0)
                else:
                    check = 0
            domain_TS[domain_TS < 0] = 0  # set all values less than 0 to 0
            new_ave_interior_height = np.average(
                domain_TS
            )  # m MHW, avg value of the domain for this storm sequence
            avg_int_height_array_100[storm_num - 1] = (
                new_ave_interior_height  # array of the avg elev for all storm sequences
            )

        # width - only for barriers that do not drown
        if outwash100_obj.barrier3d[0].drown_break == 0:
            int_width_TS = (
                np.array(outwash100_obj.barrier3d[0].InteriorWidth_AvgTS) * 10
            )  # meters, already averaged
            last_int_width = int_width_TS[-1]  # the width at the final time step
            avg_int_width_array_100[storm_num - 1] = int_width_TS[
                -1
            ]  # array of the avg widths for all storm sequences

        # dune crest
        sub_domain = outwash100_obj.barrier3d[0]._DuneDomain[
            0:tmax_100, :, :
        ]  # dune domains for current storm sequence, all model years
        dune_crest = sub_domain.max(
            axis=2
        )  # n rows are simulation years, n cols are alongshore length (dune crests)
        dune_crest_avg = (
            np.average(dune_crest) * 10
        )  # average of all years, converted to meters
        dune_crest_array_100[storm_num - 1] = (
            dune_crest_avg  # array of the average dune crests for all storm sequences
        )

        # dune gaps
        initial_gap_height = np.min(
            sub_domain[0]
        )  # get the lowest value of the dunes at year 0 without the berm elev
        berm_el = outwash100_obj.barrier3d[0].BermEl
        dune_crest_elev = dune_crest + berm_el  # dam
        # here, a dune gap is defined as any cell less than or equal to the initial gap height + berm elevation
        dune_gap_limit = (
            initial_gap_height + berm_el
        )  # the minimum elevation of a dune cell
        # if a barrier drowns, all dune cells are at the berm elevation, so here, if all the dune crest cells for a
        # single year are at the berm elevation, we change the elevation to -0.3 (water) for plotting purposes
        for row in range(np.shape(dune_crest_elev)[0]):
            if np.all(dune_crest_elev[row] == berm_el) == True:
                dune_crest_elev[row] = -0.3  # dam

        n_dune_gap_cells_array = np.zeros([np.shape(dune_crest_elev)[0], 1])
        # count the number of dune gap cells in each row (each row = a model year)
        # here, a dune gap is defined as any cell less than or equal to the initial gap height + berm elevation
        for year in range(np.shape(dune_crest_elev)[0]):
            dune_gap_row = dune_crest_elev[year]
            n_dune_gap_cells = len(dune_gap_row[dune_gap_row <= dune_gap_limit])
            n_dune_gap_cells_array[year] = (
                n_dune_gap_cells  # save the number of dune gap cells per model year
            )
        avg_dune_cells = np.average(
            n_dune_gap_cells_array
        )  # average the number dune gaps over all model years (for a single storm series)
        avg_dune_cells_array_100[storm_num - 1] = (
            avg_dune_cells * 10
        )  # each cell is 1 dam, so the total number of cells
        # gives the total length of gaps, then multiply by 10 to convert to meters

        # overwash
        QowTS = outwash100_obj.barrier3d[
            0
        ].QowTS  # m^3/m, all model years for this storm sequence
        Qow_avg = np.mean(
            QowTS
        )  # m^3/m/yr average overwash flux for this storm sequence
        avg_overwash_array_100[storm_num - 1] = (
            Qow_avg  # array of the average overwash fluxes
        )

        # outwash
        QoutTS = outwash100_obj.outwash[
            0
        ]._outwash_flux_TS  # m^3/m, all model years for this storm sequence
        QoutTS = QoutTS[
            QoutTS != 0
        ]  # remove all nonzero values (they should have been initialized as nans)
        Qout_avg = np.mean(
            QoutTS
        )  # m^3/m/yr average outwash flux for this storm sequence
        avg_outwash_array_100[storm_num - 1] = (
            Qout_avg  # array of the average outwash fluxes
        )

        # 50% variables ------------------------------------------------------------------------------------------------
        filename_50 = (
            "config{}_outwash50_startyr1_interval{}yrs_Slope0pt03_{}.npz".format(
                config, storm_interval, storm_num
            )
        )
        file_50 = datadir_50 + filename_50
        outwash50 = np.load(file_50, allow_pickle=True)
        outwash50_obj = outwash50["cascade"][0]
        tmax_50 = outwash50_obj.barrier3d[0].TMAX + 1

        # drowning variable
        drowning_array_50[storm_num - 1] = outwash50_obj.barrier3d[0].drown_break

        # shoreline position - only for barriers that do not drown
        m_xsTS = np.subtract(
            outwash50_obj.barrier3d[0].x_s_TS, outwash50_obj.barrier3d[0].x_s_TS[0]
        )
        m_xsTS = np.multiply(m_xsTS, 10)
        if outwash50_obj.barrier3d[0].drown_break == 0:
            end_net_migration_no_drowns_50.append(m_xsTS[-1])

        # average elevation - only for barriers that do not drown
        if outwash50_obj.barrier3d[0].drown_break == 0:
            domain_TS = outwash50_obj.barrier3d[0].DomainTS[-1] * 10  # m MHW
            check = 1
            # remove all alongshore rows that are completely submerged less than 0
            while check == 1:
                if all(x <= 0 for x in domain_TS[-1, :]):
                    domain_TS = np.delete(domain_TS, -1, axis=0)
                else:
                    check = 0
            domain_TS[domain_TS < 0] = 0  # set all values less than 0 to 0
            new_ave_interior_height = np.average(
                domain_TS
            )  # all in m MHW, average value of the domain for this storm sequence
            avg_int_height_array_50[storm_num - 1] = (
                new_ave_interior_height  # array of the average elevations for all storm sequences
            )

        # width - only for barriers that do not drown
        if outwash50_obj.barrier3d[0].drown_break == 0:
            int_width_TS = (
                np.array(outwash50_obj.barrier3d[0].InteriorWidth_AvgTS) * 10
            )  # meters, already averaged values
            last_int_width = int_width_TS[-1]  # the width at the final time step
            avg_int_width_array_50[storm_num - 1] = int_width_TS[
                -1
            ]  # array of the average widths for all storm sequences

        # dune crest
        sub_domain = outwash50_obj.barrier3d[0]._DuneDomain[
            0:tmax_50, :, :
        ]  # dune domains for current storm sequence, all model years
        dune_crest = sub_domain.max(
            axis=2
        )  # n rows are simulation years, n cols are alongshore length (dune crests)
        dune_crest_avg = (
            np.average(dune_crest) * 10
        )  # average of all years, converted to meters
        dune_crest_array_50[storm_num - 1] = (
            dune_crest_avg  # array of the average dune crests for all storm sequences
        )

        # dune gaps
        initial_gap_height = np.min(
            sub_domain[0]
        )  # get the lowest value of the dunes at year 0 without the berm elev
        berm_el = outwash50_obj.barrier3d[0].BermEl
        dune_crest_elev = dune_crest + berm_el  # dam
        # here, a dune gap is defined as any cell less than or equal to the initial gap height + berm elevation
        dune_gap_limit = (
            initial_gap_height + berm_el
        )  # the minimum elevation of a dune cell
        # if a barrier drowns, all dune cells are at the berm elevation, so here, if all the dune crest cells for a
        # single year are at the berm elevation, we change the elevation to -0.3 (water) for plotting purposes
        for row in range(np.shape(dune_crest_elev)[0]):
            if np.all(dune_crest_elev[row] == berm_el) == True:
                dune_crest_elev[row] = -0.3  # dam

        n_dune_gap_cells_array = np.zeros([np.shape(dune_crest_elev)[0], 1])
        # count the number of dune gap cells in each row (each row = a model year)
        # here, a dune gap is defined as any cell less than or equal to the initial gap height + berm elevation
        for year in range(np.shape(dune_crest_elev)[0]):
            dune_gap_row = dune_crest_elev[year]
            n_dune_gap_cells = len(dune_gap_row[dune_gap_row <= dune_gap_limit])
            n_dune_gap_cells_array[year] = (
                n_dune_gap_cells  # save the number of dune gap cells per model year
            )
        avg_dune_cells = np.average(
            n_dune_gap_cells_array
        )  # average the number dune gaps over all model years (for a single storm series)
        avg_dune_cells_array_50[storm_num - 1] = (
            avg_dune_cells * 10
        )  # each cell is 1 dam, so the total number of cells
        # gives the total length of gaps, then multiply by 10 to convert to meters

        # overwash
        QowTS = outwash50_obj.barrier3d[
            0
        ].QowTS  # m^3/m, all model years for this storm sequence
        Qow_avg = np.mean(
            QowTS
        )  # m^3/m/yr average overwash flux for this storm sequence
        avg_overwash_array_50[storm_num - 1] = (
            Qow_avg  # array of the average overwash fluxes
        )

        # outwash
        QoutTS = outwash50_obj.outwash[
            0
        ]._outwash_flux_TS  # m^3/m, all model years for this storm sequence
        QoutTS = QoutTS[QoutTS != 0]
        Qout_avg = np.mean(
            QoutTS
        )  # m^3/m/yr average outwash flux for this storm sequence
        avg_outwash_array_50[storm_num - 1] = (
            Qout_avg  # array of the average outwash fluxes
        )

        # 0% variables -------------------------------------------------------------------------------------------------
        filename_0 = (
            "config{}_outwash0_startyr1_interval{}yrs_Slope0pt03_{}.npz".format(
                config, storm_interval, storm_num
            )
        )
        file_0 = datadir_0 + filename_0
        outwash0 = np.load(file_0, allow_pickle=True)
        outwash0_obj = outwash0["cascade"][0]
        tmax_0 = outwash0_obj.barrier3d[0].TMAX + 1

        # drowning variable
        drowning_array_0[storm_num - 1] = outwash0_obj.barrier3d[0].drown_break

        # shoreline position - only for barriers that do not drown
        m_xsTS = np.subtract(
            outwash0_obj.barrier3d[0].x_s_TS, outwash0_obj.barrier3d[0].x_s_TS[0]
        )
        m_xsTS = np.multiply(m_xsTS, 10)
        if outwash0_obj.barrier3d[0].drown_break == 0:
            end_net_migration_no_drowns_0.append(m_xsTS[-1])

        # average elevation - only for barriers that do not drown
        if outwash0_obj.barrier3d[0].drown_break == 0:
            domain_TS = outwash0_obj.barrier3d[0].DomainTS[-1] * 10  # m MHW
            check = 1
            # remove all alongshore rows that are completely submerged less than 0
            while check == 1:
                if all(x <= 0 for x in domain_TS[-1, :]):
                    domain_TS = np.delete(domain_TS, -1, axis=0)
                else:
                    check = 0
            domain_TS[domain_TS < 0] = 0  # set all values less than 0 to 0
            new_ave_interior_height = np.average(
                domain_TS
            )  # all in m MHW, average value of the domain for this storm sequence
            avg_int_height_array_0[storm_num - 1] = (
                new_ave_interior_height  # array of the average elevations for all storm sequences
            )

        # width - only for barriers that do not drown
        if outwash0_obj.barrier3d[0].drown_break == 0:
            int_width_TS = (
                np.array(outwash0_obj.barrier3d[0].InteriorWidth_AvgTS) * 10
            )  # meters, already averaged values
            last_int_width = int_width_TS[-1]  # the width at the final time step
            avg_int_width_array_0[storm_num - 1] = int_width_TS[
                -1
            ]  # array of the average widths for all storm sequences

        # dune crest
        sub_domain = outwash0_obj.barrier3d[0]._DuneDomain[
            0:tmax_0, :, :
        ]  # dune domains for current storm sequence, all model years
        dune_crest = sub_domain.max(
            axis=2
        )  # n rows are simulation years, n cols are alongshore length (dune crests)
        dune_crest_avg = (
            np.average(dune_crest) * 10
        )  # average of all years, converted to meters
        dune_crest_array_0[storm_num - 1] = (
            dune_crest_avg  # array of the average dune crests for all storm sequences
        )

        # dune gaps
        initial_gap_height = np.min(
            sub_domain[0]
        )  # get the lowest value of the dunes at year 0 without the berm elev
        berm_el = outwash0_obj.barrier3d[0].BermEl
        dune_crest_elev = dune_crest + berm_el  # dam
        # here, a dune gap is defined as any cell less than or equal to the initial gap height + berm elevation
        dune_gap_limit = (
            initial_gap_height + berm_el
        )  # the minimum elevation of a dune cell
        # if a barrier drowns, all dune cells are at the berm elevation, so here, if all the dune crest cells for a
        # single year are at the berm elevation, we change the elevation to -0.3 (water) for plotting purposes
        for row in range(np.shape(dune_crest_elev)[0]):
            if np.all(dune_crest_elev[row] == berm_el) == True:
                dune_crest_elev[row] = -0.3  # dam

        n_dune_gap_cells_array = np.zeros([np.shape(dune_crest_elev)[0], 1])
        # count the number of dune gap cells in each row (each row = a model year)
        # here, a dune gap is defined as any cell less than or equal to the initial gap height + berm elevation
        for year in range(np.shape(dune_crest_elev)[0]):
            dune_gap_row = dune_crest_elev[year]
            n_dune_gap_cells = len(dune_gap_row[dune_gap_row <= dune_gap_limit])
            n_dune_gap_cells_array[year] = (
                n_dune_gap_cells  # save the number of dune gap cells per model year
            )
        avg_dune_cells = np.average(
            n_dune_gap_cells_array
        )  # average the number dune gaps over all model years (for a single storm series)
        avg_dune_cells_array_0[storm_num - 1] = (
            avg_dune_cells * 10
        )  # each cell is 1 dam, so the total number of cells
        # gives the total length of gaps, then multiply by 10 to convert to meters

        # overwash
        QowTS = outwash0_obj.barrier3d[
            0
        ].QowTS  # m^3/m, all model years for this storm sequence
        Qow_avg = np.mean(
            QowTS
        )  # m^3/m/yr average overwash flux for this storm sequence
        avg_overwash_array_0[storm_num - 1] = (
            Qow_avg  # array of the average overwash fluxes
        )

        # outwash
        QoutTS = outwash0_obj.outwash[
            0
        ]._outwash_flux_TS  # m^3/m, all model years for this storm sequence
        QoutTS = QoutTS[QoutTS != 0]
        Qout_avg = np.mean(
            QoutTS
        )  # m^3/m/yr average outwash flux for this storm sequence
        avg_outwash_array_0[storm_num - 1] = (
            Qout_avg  # array of the average outwash fluxes
        )

    # ------------------------------------------------------------------------------------------------------------------
    # now that we have a value for all storm sequences, find average for each scenario

    # drowning stats - because we have 100 storms and drowns = 1, the sum of the array is the percent that drown
    percent_drown_b3d = np.sum(drowning_array_b3d)
    percent_drown_100 = np.sum(drowning_array_100)
    percent_drown_50 = np.sum(drowning_array_50)
    percent_drown_0 = np.sum(drowning_array_0)

    # shoreline position - only for barriers that do not drown
    avg_shoreline_b3d = round(np.average(end_net_migration_no_drowns_b3d), 0)
    avg_shoreline_100 = round(np.average(end_net_migration_no_drowns_100), 0)
    avg_shoreline_50 = round(np.average(end_net_migration_no_drowns_50), 0)
    avg_shoreline_0 = round(np.average(end_net_migration_no_drowns_0), 0)

    # interior elevation - only for barriers that do not drown
    # remove nonzero values first
    avg_int_height_array_b3d = avg_int_height_array_b3d[avg_int_height_array_b3d != 0]
    avg_int_height_array_100 = avg_int_height_array_100[avg_int_height_array_100 != 0]
    avg_int_height_array_50 = avg_int_height_array_50[avg_int_height_array_50 != 0]
    avg_int_height_array_0 = avg_int_height_array_0[avg_int_height_array_0 != 0]

    avg_height_last_TS_b3d = round(np.average(avg_int_height_array_b3d), 2)
    avg_height_last_TS_100 = round(np.average(avg_int_height_array_100), 2)
    avg_height_last_TS_50 = round(np.average(avg_int_height_array_50), 2)
    avg_height_last_TS_0 = round(np.average(avg_int_height_array_0), 2)

    # interior width - only for barriers that do not drown
    avg_int_width_array_b3d = avg_int_width_array_b3d[avg_int_width_array_b3d != 0]
    avg_int_width_array_100 = avg_int_width_array_100[avg_int_width_array_100 != 0]
    avg_int_width_array_50 = avg_int_width_array_50[avg_int_width_array_50 != 0]
    avg_int_width_array_0 = avg_int_width_array_0[avg_int_width_array_0 != 0]

    avg_width_last_TS_b3d = round(np.average(avg_int_width_array_b3d), 0)
    avg_width_last_TS_100 = round(np.average(avg_int_width_array_100), 0)
    avg_width_last_TS_50 = round(np.average(avg_int_width_array_50), 0)
    avg_width_last_TS_0 = round(np.average(avg_int_width_array_0), 0)

    # dune crest
    avg_dune_crest_array_b3d = round(np.average(dune_crest_array_b3d), 2)
    avg_dune_crest_array_100 = round(np.average(dune_crest_array_100), 2)
    avg_dune_crest_array_50 = round(np.average(dune_crest_array_50), 2)
    avg_dune_crest_array_0 = round(np.average(dune_crest_array_0), 2)

    # dune gap lengths
    avg_dune_cells_array_b3d = round(np.average(avg_dune_cells_array_b3d), 0)
    avg_dune_cells_array_100 = round(np.average(avg_dune_cells_array_100), 0)
    avg_dune_cells_array_50 = round(np.average(avg_dune_cells_array_50), 0)
    avg_dune_cells_array_0 = round(np.average(avg_dune_cells_array_0), 0)

    # overwash
    avg_overwash_array_b3d = round(np.average(avg_overwash_array_b3d), 0)
    avg_overwash_array_100 = round(np.average(avg_overwash_array_100), 0)
    avg_overwash_array_50 = round(np.average(avg_overwash_array_50), 0)
    avg_overwash_array_0 = round(np.average(avg_overwash_array_0), 0)

    # outwash - note there is no b3d value since there is no outwash for that scenario
    avg_outwash_array_100 = round(np.average(avg_outwash_array_100), 0)
    avg_outwash_array_50 = round(np.average(avg_outwash_array_50), 0)
    avg_outwash_array_0 = round(np.average(avg_outwash_array_0), 0)

    # printing results to console
    print(f"{rname} \n")

    if drowning_stats:
        print(
            "{}% of barriers drown for the baseline scenario".format(
                int(percent_drown_b3d)
            )
        )
        print(
            "{}% of barriers drown for the 100% outwash to shoreface scenario".format(
                int(percent_drown_100)
            )
        )
        print(
            "{}% of barriers drown for the 50% outwash to shoreface scenario".format(
                int(percent_drown_50)
            )
        )
        print(
            "{}% of barriers drown for the 0% outwash to shoreface scenario".format(
                int(percent_drown_0)
            )
        )
        print()

    if migration_stats:
        print(
            "baseline avg net migration excluding drown years: {}".format(
                avg_shoreline_b3d
            )
        )
        print(
            "100% outwash avg net migration excluding drown years: {}".format(
                avg_shoreline_100
            )
        )
        print(
            "50% outwash avg net migration excluding drown years: {}".format(
                avg_shoreline_50
            )
        )
        print(
            "0% outwash avg net migration excluding drown years: {}".format(
                avg_shoreline_0
            )
        )
        print()

    if geomoetry_stats:
        # elevation
        print(
            "baseline avg interior elevation excluding drown years: {}".format(
                avg_height_last_TS_b3d
            )
        )
        print(
            "100% outwash avg interior elevation excluding drown years: {}".format(
                avg_height_last_TS_100
            )
        )
        print(
            "50% outwash avg interior elevation excluding drown years: {}".format(
                avg_height_last_TS_50
            )
        )
        print(
            "0% outwash avg interior elevation excluding drown years: {}".format(
                avg_height_last_TS_0
            )
        )
        # width
        print(
            "baseline avg interior width excluding drown years: {}".format(
                avg_width_last_TS_b3d
            )
        )
        print(
            "100% outwash avg interior width excluding drown years: {}".format(
                avg_width_last_TS_100
            )
        )
        print(
            "50% outwash avg interior width excluding drown years: {}".format(
                avg_width_last_TS_50
            )
        )
        print(
            "0% outwash avg interior width excluding drown years: {}".format(
                avg_width_last_TS_0
            )
        )
        print()

    # dune crest avg elevations and plots and dune gaps
    if dune_stats:
        print(f"average dune crest height baseline: {avg_dune_crest_array_b3d}")
        print(
            "average dune crest height 100% outwash: {}".format(
                avg_dune_crest_array_100
            )
        )
        print(f"average dune crest height 50% outwash: {avg_dune_crest_array_50}")
        print(f"average dune crest height 0% outwash: {avg_dune_crest_array_0}")

        # dune gaps
        print(
            "The average total length of dune gaps for the b3d scenario is: {} m".format(
                avg_dune_cells_array_b3d
            )
        )
        print(
            "The average total length of dune gaps for the 100% scenario is: {} m".format(
                avg_dune_cells_array_100
            )
        )
        print(
            "The average total length of dune gaps for the 50% scenario is: {} m".format(
                avg_dune_cells_array_50
            )
        )
        print(
            "The average total length of dune gaps for the 0% scenario is: {} m".format(
                avg_dune_cells_array_0
            )
        )
        print()

    # average overwash and outwash
    if flux_stats:
        print(f"average overwash per event baseline: {avg_overwash_array_b3d}")
        print(f"average overwash per event 100%: {avg_overwash_array_100}")
        print(f"average overwash per event 50%: {avg_overwash_array_50}")
        print(f"average overwash per event 0%: {avg_overwash_array_0}")

        print(f"average outwash per event 100%: {avg_outwash_array_100}")
        print(f"average outwash per event 50%: {avg_outwash_array_50}")
        print(f"average outwash per event 0%: {avg_outwash_array_0}")
        print()
