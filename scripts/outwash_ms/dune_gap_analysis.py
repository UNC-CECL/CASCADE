# running plotters_outwash to make dune and elevation figures

# Lexi Van Blunk
# 2/23/2024
# 3/16/2025 updated to the new directory (re-runs)

from matplotlib import pyplot as plt
import numpy as np

rname_array = ["r025", "r035"]
# rname_array = ["r025"]

for rname in rname_array:
    storm_interval = 20        # 20 or 10 years (we did not do 10 years)
    config = 4                 # 1, 2, 3, or 4
    avg_dune_cells_array_b3d = np.zeros(100)
    avg_dune_cells_array_100 = np.zeros(100)
    avg_dune_cells_array_50 = np.zeros(100)
    avg_dune_cells_array_0 = np.zeros(100)

    for storm_num in range(1,101):

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

        ### ----------------------------------- Dune Gap Figures ---------------------------------------

        TMAX = 101

        # Barrier3d only
        DuneCrest = []

        for iB3D in range(len(b3d)):
            sub_domain = b3d[iB3D]._DuneDomain[0:TMAX, :, :]
            initial_gap_elev = np.min(sub_domain[0])
            DuneCrest.append(sub_domain.max(axis=2))

        DuneCrest = np.hstack(DuneCrest).astype(float)
        berm_el = b3d[0].BermEl  # dam
        DuneCrest = DuneCrest + berm_el

        for row in range(np.shape(DuneCrest)[0]):
            if np.all(DuneCrest[row] == berm_el) == True:
                DuneCrest[row] = -0.3

        n_dune_gap_cells_array = np.zeros([np.shape(DuneCrest)[0],1])
        # count the number of dune gap cells in each row (each row = a model year)
        for year in range(np.shape(DuneCrest)[0]):
            dune_gap_row = DuneCrest[year]
            n_dune_gap_cells = len(dune_gap_row[dune_gap_row<=initial_gap_elev+berm_el])
            n_dune_gap_cells_array[year] = n_dune_gap_cells
        avg_dune_cells = np.average(n_dune_gap_cells_array)
        avg_dune_cells_array_b3d[storm_num-1] = avg_dune_cells * 10  # each cell is 1 dam, so the total number of cells
        # gives the total length of gaps, then multiply by 10 to convert to meters

        # outwash 100%
        DuneCrest = []

        for iB3D in range(len(outwash100)):
            sub_domain = outwash100[iB3D]._DuneDomain[0:TMAX, :, :]
            initial_gap_elev = np.min(sub_domain[0])
            DuneCrest.append(sub_domain.max(axis=2))

        DuneCrest = np.hstack(DuneCrest).astype(float)
        berm_el = outwash100[0].BermEl  # dam
        DuneCrest = DuneCrest + berm_el

        for row in range(np.shape(DuneCrest)[0]):
            if np.all(DuneCrest[row] == berm_el) == True:
                DuneCrest[row] = -0.3

        n_dune_gap_cells_array = np.zeros([np.shape(DuneCrest)[0],1])
        # count the number of dune gap cells in each row (each row = a model year)
        for year in range(np.shape(DuneCrest)[0]):
            dune_gap_row = DuneCrest[year]
            n_dune_gap_cells = len(dune_gap_row[dune_gap_row<=initial_gap_elev+berm_el])
            n_dune_gap_cells_array[year] = n_dune_gap_cells
        avg_dune_cells = np.average(n_dune_gap_cells_array)
        avg_dune_cells_array_100[storm_num-1] = avg_dune_cells * 10  # each cell is 1 dam, so the total number of cells
        # gives the total length of gaps, then multiply by 10 to convert to meters

        # outwash 50%
        DuneCrest = []

        for iB3D in range(len(outwash50)):
            sub_domain = outwash50[iB3D]._DuneDomain[0:TMAX, :, :]
            initial_gap_elev = np.min(sub_domain[0])
            DuneCrest.append(sub_domain.max(axis=2))

        DuneCrest = np.hstack(DuneCrest).astype(float)
        berm_el = outwash50[0].BermEl  # dam
        DuneCrest = DuneCrest + berm_el

        for row in range(np.shape(DuneCrest)[0]):
            if np.all(DuneCrest[row] == berm_el) == True:
                DuneCrest[row] = -0.3

        n_dune_gap_cells_array = np.zeros([np.shape(DuneCrest)[0],1])
        # count the number of dune gap cells in each row (each row = a model year)
        for year in range(np.shape(DuneCrest)[0]):
            dune_gap_row = DuneCrest[year]
            n_dune_gap_cells = len(dune_gap_row[dune_gap_row<=initial_gap_elev+berm_el])
            n_dune_gap_cells_array[year] = n_dune_gap_cells
        avg_dune_cells = np.average(n_dune_gap_cells_array)
        avg_dune_cells_array_50[storm_num-1] = avg_dune_cells * 10  # each cell is 1 dam, so the total number of cells
        # gives the total length of gaps, then multiply by 10 to convert to meters

        # outwash 0%
        DuneCrest = []

        for iB3D in range(len(outwash0)):
            sub_domain = outwash0[iB3D]._DuneDomain[0:TMAX, :, :]
            initial_gap_elev = np.min(sub_domain[0])
            DuneCrest.append(sub_domain.max(axis=2))

        DuneCrest = np.hstack(DuneCrest).astype(float)
        berm_el = outwash0[0].BermEl  # dam
        DuneCrest = DuneCrest + berm_el

        for row in range(np.shape(DuneCrest)[0]):
            if np.all(DuneCrest[row] == berm_el) == True:
                DuneCrest[row] = -0.3

        n_dune_gap_cells_array = np.zeros([np.shape(DuneCrest)[0],1])
        # count the number of dune gap cells in each row (each row = a model year)
        for year in range(np.shape(DuneCrest)[0]):
            dune_gap_row = DuneCrest[year]
            n_dune_gap_cells = len(dune_gap_row[dune_gap_row<=initial_gap_elev+berm_el])
            n_dune_gap_cells_array[year] = n_dune_gap_cells
        avg_dune_cells = np.average(n_dune_gap_cells_array)
        avg_dune_cells_array_0[storm_num-1] = avg_dune_cells * 10  # each cell is 1 dam, so the total number of cells
        # gives the total length of gaps, then multiply by 10 to convert to meters

    avg_dune_cells_b3d = np.round(np.average(avg_dune_cells_array_b3d),0)
    avg_dune_cells_100 = np.round(np.average(avg_dune_cells_array_100),0)
    avg_dune_cells_50 = np.round(np.average(avg_dune_cells_array_50),0)
    avg_dune_cells_0 = np.round(np.average(avg_dune_cells_array_0),0)

    print("{0}".format(rname))
    print("The average number of dune gap cells for the b3d% scenario is: {0} m".format(avg_dune_cells_b3d))
    print("The average number of dune gap cells for the 100% scenario is: {0} m".format(avg_dune_cells_100))
    print("The average number of dune gap cells for the 50% scenario is: {0} m".format(avg_dune_cells_50))
    print("The average number of dune gap cells for the 0% scenario is: {0} m".format(avg_dune_cells_0))


