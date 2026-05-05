"""
tests if the Cascade version of BarrierBMFT matches the output from the original version of BarrierBMFT

Lexi (Van Blunk) Fiegelist
Last updated: 5 May 2026
"""

import numpy as np
import os
import pandas as pd
import matplotlib.pyplot as plt
import imageio


def check_init_values(bmft_instance, cascade_instance):
    """
    compares the initilization of two classes
    here, we use it to compare the mainland Bmftc, back-barrier Bmftc, and BarrierBMFT classes
    :param bmft_instance: the class instance using the original BarrierBMFT code
    :param cascade_instance: the class instance using Cascade with marsh dynamics enabled
    :return: dataframe with all attributes and whether the two instances are equal or not
    """
    marsh_cascade_attributes = dir(cascade_instance)
    marsh_attributes = dir(bmft_instance)
    # check if number of attributes are the same
    if len(marsh_cascade_attributes) == len(marsh_attributes):
        len_same = True
    else:
        len_same = False
    print("same number of attributes: {0}".format(len_same))

    # check if values are the same and return a dataframe with the attribute and whether they match
    # need different code for arrays
    # turn the cascade isntance into a dictionary so we can access keys and values
    marsh_cascade_dict = vars(cascade_instance)
    marsh_dict = vars(bmft_instance)
    keys = marsh_cascade_dict.keys()
    # keys should be the same for all runs
    key_list = []
    for k in keys:
        # print(k)
        marsh_cascade_val = marsh_cascade_dict[k]
        marsh_val = marsh_dict[k]
        if isinstance(marsh_cascade_val, np.ndarray) or isinstance(marsh_cascade_val, list):
            try:
                error = np.testing.assert_array_almost_equal(marsh_val, marsh_cascade_val)
                if error is None:
                    value_same = True
                else:
                    value_same = False
            except AssertionError:
                value_same = False
            except ValueError:  # try looping through each individual item to compare
                value_same_at_index = 0  # keep track of how many items are equal
                for i in range(len(marsh_cascade_val)):
                    # check if the individual item is an array or list
                    if isinstance(marsh_cascade_val[i], np.ndarray) or isinstance(marsh_cascade_val[i], list):
                        error = np.testing.assert_array_almost_equal(marsh_val[i], marsh_cascade_val[i])
                        if error is None:
                            value_same_at_index += 1  # if they are the same, add to list
                    else:  # if not an array or list, compare values directly
                        if marsh_cascade_val[i] == marsh_val[i]:
                            value_same_at_index += 1
                if value_same_at_index == len(marsh_cascade_val):  # if all the values are the same,
                    # then the attributes match
                    value_same = True
                else:
                    value_same = False
            key_list.append([k, value_same])


        else:
            if marsh_cascade_val == marsh_val:
                value_same = True
            else:
                value_same = False
            key_list.append([k, value_same])
    key_df = pd.DataFrame(key_list, columns=["attribute", "value"])

    return key_df


def plot_2_cascade_transects(model1, model2, save_dir, newfolder):
    for t in range(int(model1.bmftc.dur)-1):

        # CASCADE results
        # Combine transects
        BB_transect_m1 = np.flip(model1.bmftc_BB.elevation[model1.bmftc_BB.startyear + t - 1, int(model1.bmftc_BB.Marsh_edge[model1.bmftc_BB.startyear + t]):])
        if model1.x_b_TS_ML[t] < 0:
            ML_transect_m1 = np.append(np.ones([abs(int(model1.x_b_TS_ML[t]))]) * model1.bmftc_ML.elevation[model1.bmftc_ML.startyear + t - 1, 1], model1.bmftc_ML.elevation[model1.bmftc_ML.startyear + t - 1, :])
        elif model1.x_b_TS_ML[t] > 0:
            ML_transect_m1 = model1.bmftc_ML.elevation[model1.bmftc_ML.startyear + t - 1, int(model1.x_b_TS_ML[t]):]

        whole_transect_m1 = np.append(BB_transect_m1, ML_transect_m1)

        # Extract marsh & forest locations to plot
        x_forest_m1 = [model1.bmftc_BB.B - int(model1.bmftc_BB.Forest_edge[model1.bmftc_BB.startyear + t]), len(BB_transect_m1) - int(model1.x_b_TS_ML[t]) + int(model1.bmftc_ML.Forest_edge[model1.bmftc_ML.startyear + t])]
        x_marsh_m1 = [len(BB_transect_m1) - 1, len(BB_transect_m1) - int(model1.x_b_TS_ML[t]) + int(model1.bmftc_ML.Marsh_edge[model1.bmftc_ML.startyear + t])]
        y_forest_m1 = whole_transect_m1[x_forest_m1]
        y_marsh_m1 = whole_transect_m1[x_marsh_m1]

        # Combine transects
        BB_transect_m2 = np.flip(model2.bmftc_BB.elevation[model2.bmftc_BB.startyear + t - 1, int(model2.bmftc_BB.Marsh_edge[model2.bmftc_BB.startyear + t]):])
        if model2.x_b_TS_ML[t] < 0:
            ML_transect_m2 = np.append(np.ones([abs(int(model2.x_b_TS_ML[t]))]) * model2.bmftc_ML.elevation[model2.bmftc_ML.startyear + t - 1, 1], model2.bmftc_ML.elevation[model2.bmftc_ML.startyear + t - 1, :])
        elif model2.x_b_TS_ML[t] > 0:
            ML_transect_m2 = model2.bmftc_ML.elevation[model2.bmftc_ML.startyear + t - 1, int(model2.x_b_TS_ML[t]):]

        whole_transect_m2 = np.append(BB_transect_m2, ML_transect_m2)

        # Extract marsh & forest locations to plot
        x_forest_m2 = [model2.bmftc_BB.B - int(model2.bmftc_BB.Forest_edge[model2.bmftc_BB.startyear + t]), len(BB_transect_m2) - int(model2.x_b_TS_ML[t]) + int(model2.bmftc_ML.Forest_edge[model2.bmftc_ML.startyear + t])]
        x_marsh_m2 = [len(BB_transect_m2) - 1, len(BB_transect_m2) - int(model2.x_b_TS_ML[t]) + int(model2.bmftc_ML.Marsh_edge[model2.bmftc_ML.startyear + t])]
        y_forest_m2 = whole_transect_m2[x_forest_m2]
        y_marsh_m2 = whole_transect_m2[x_marsh_m2]
        sizes = np.ones(2) * 100

        # Plot and save
        transectFig = plt.figure(figsize=(15, 7))
        plt.rcParams.update({"font.size": 11})
        plt.plot(whole_transect_m2, c="black")
        plt.plot(whole_transect_m1, c="magenta", ls="dashed")
        plt.scatter(x_forest_m1, y_forest_m1, c="green", marker="x", s=sizes)
        plt.scatter(x_marsh_m1, y_marsh_m1, c="brown", marker="x", s=sizes)
        plt.scatter(x_forest_m2, y_forest_m2, c="green")
        plt.scatter(x_marsh_m2, y_marsh_m2, c="brown")
        plt.ylim(-2.1, 4)
        plt.xlabel("Cross-shore Distance (m)")
        plt.ylabel("Elevation (m)")
        plt.tight_layout()
        timestr = "Time = " + str(t) + " yrs"
        plt.text(0, 0, timestr)
        newpath = os.path.join(save_dir, newfolder)
        if not os.path.exists(newpath):
            os.makedirs(newpath)
        plt.text(0, 0, timestr)
        name = os.path.join(newpath, "transect_{0}.png".format(t))
        transectFig.savefig(name)  # dpi=200
        plt.close(transectFig)
    #
    frames = []
    for filenum in range(int(model1.bmftc.dur)-1):
        filename = os.path.join(newpath, "transect_{0}.png".format(filenum))
        frames.append(imageio.imread(filename))
    imageio.mimsave(os.path.join(newpath, "transect.gif"), frames, "GIF-FI", fps=5)
    print()
    print("[ * GIF successfully generated * ]")

    # ===========
    plt.show()


def plot_1_casc_1_barbmftc_transects(model1, bmftc_ml, bmftc_bb, xbts_ml, save_dir, newfolder):
    for t in range(int(model1.bmftc.dur)-1):

        # Combine transects - cascade results
        BB_transect_m1 = np.flip(model1.bmftc_BB.elevation[model1.bmftc_BB.startyear + t - 1, int(model1.bmftc_BB.Marsh_edge[model1.bmftc_BB.startyear + t]):])
        if model1.x_b_TS_ML[t] < 0:
            ML_transect_m1 = np.append(np.ones([abs(int(model1.x_b_TS_ML[t]))]) * model1.bmftc_ML.elevation[model1.bmftc_ML.startyear + t - 1, 1], model1.bmftc_ML.elevation[model1.bmftc_ML.startyear + t - 1, :])
        elif model1.x_b_TS_ML[t] > 0:
            ML_transect_m1 = model1.bmftc_ML.elevation[model1.bmftc_ML.startyear + t - 1, int(model1.x_b_TS_ML[t]):]

        whole_transect_m1 = np.append(BB_transect_m1, ML_transect_m1)

        # Extract marsh & forest locations to plot
        x_forest_m1 = [model1.bmftc_BB.B - int(model1.bmftc_BB.Forest_edge[model1.bmftc_BB.startyear + t]), len(BB_transect_m1) - int(model1.x_b_TS_ML[t]) + int(model1.bmftc_ML.Forest_edge[model1.bmftc_ML.startyear + t])]
        x_marsh_m1 = [len(BB_transect_m1) - 1, len(BB_transect_m1) - int(model1.x_b_TS_ML[t]) + int(model1.bmftc_ML.Marsh_edge[model1.bmftc_ML.startyear + t])]
        y_forest_m1 = whole_transect_m1[x_forest_m1]
        y_marsh_m1 = whole_transect_m1[x_marsh_m1]

        # Combine transects - barrierbmftc results
        BB_transect_m2 = np.flip(bmftc_bb.elevation[bmftc_bb.startyear + t - 1, int(bmftc_bb.Marsh_edge[bmftc_bb.startyear + t]):])
        if xbts_ml[t] < 0:
            ML_transect_m2 = np.append(np.ones([abs(int(xbts_ml[t]))]) * bmftc_ml.elevation[bmftc_ml.startyear + t - 1, 1], bmftc_ml.elevation[bmftc_ml.startyear + t - 1, :])
        elif xbts_ml[t] > 0:
            ML_transect_m2 = bmftc_ml.elevation[bmftc_ml.startyear + t - 1, int(xbts_ml[t]):]

        whole_transect_m2 = np.append(BB_transect_m2, ML_transect_m2)

        # Extract marsh & forest locations to plot
        x_forest_m2 = [bmftc_bb.B - int(bmftc_bb.Forest_edge[bmftc_bb.startyear + t]), len(BB_transect_m2) - int(xbts_ml[t]) + int(bmftc_ml.Forest_edge[bmftc_ml.startyear + t])]
        x_marsh_m2 = [len(BB_transect_m2) - 1, len(BB_transect_m2) - int(xbts_ml[t]) + int(bmftc_ml.Marsh_edge[bmftc_ml.startyear + t])]
        y_forest_m2 = whole_transect_m2[x_forest_m2]
        y_marsh_m2 = whole_transect_m2[x_marsh_m2]

        # Plot and save
        sizes = np.ones(2) * 100
        transectFig = plt.figure(figsize=(15, 7))
        plt.rcParams.update({"font.size": 11})
        plt.plot(whole_transect_m2, c="black")
        plt.plot(whole_transect_m1, c="magenta", ls="dashed")
        plt.scatter(x_forest_m1, y_forest_m1, c="green", marker="x", s=sizes)
        plt.scatter(x_marsh_m1, y_marsh_m1, c="brown", marker="x", s=sizes)
        plt.scatter(x_forest_m2, y_forest_m2, c="green")
        plt.scatter(x_marsh_m2, y_marsh_m2, c="brown")
        plt.ylim(-2.1, 4)
        plt.xlabel("Cross-shore Distance (m)")
        plt.ylabel("Elevation (m)")
        plt.tight_layout()
        timestr = "Time = " + str(t) + " yrs"
        plt.text(0, 0, timestr)
        newpath = os.path.join(save_dir, newfolder)
        if not os.path.exists(newpath):
            os.makedirs(newpath)
        plt.text(0, 0, timestr)
        name = os.path.join(newpath, "transect_{0}.png".format(t))
        transectFig.savefig(name)  # dpi=200
        plt.close(transectFig)

    frames = []
    for filenum in range(int(model1.bmftc.dur)-1):
        filename = os.path.join(newpath, "transect_{0}.png".format(filenum))
        frames.append(imageio.imread(filename))
    imageio.mimsave(os.path.join(newpath, "transect.gif"), frames, "GIF-FI", fps=5)
    print()
    print("[ * GIF successfully generated * ]")

    # ===========
    plt.show()


# ----------------------------------------------------------------------------------------------------------------------
plots_on = False

# location of the npz files
cascade_datadir = r"C:/Users/Lexi\Documents\UNC\BarrierBMFT\tests\cascade_runs"
barrierbmft_datadir = r"C:/Users/Lexi\Documents\UNC\BarrierBMFT\tests\orig_barrierbmft"

# ======================================================================================================================
# enter classes to compare
# ======================================================================================================================
cascade_compl = os.path.join(cascade_datadir, "marsh_cascade_TS20_edit2.npz")  # this includes ALL classes (bmft and b3d)
cascade_compl2 = os.path.join(cascade_datadir, "marsh_cascade_TS20_edit2_moved_marsh_update.npz")  # this includes ALL classes (bmft and b3d)
# note: got an error when trying to save full Bmftc class, so had to save mainland/backbarrier bmftc and b3d classes
# seperately instead of having them all in one Bmftc class
bb_complete = os.path.join(barrierbmft_datadir, "backbarrier_orig_TS20_edit2.npz")
ml_complete = os.path.join(barrierbmft_datadir, "mainland_orig_TS20_edit2.npz")
b3d_bmft_complete = os.path.join(barrierbmft_datadir, "b3d_orig_TS20_edit2.npz")  # b3d in bmft

# ======================================================================================================================
# load the classes
# ======================================================================================================================
cascade_comp_inst = np.load(cascade_compl, allow_pickle=True)['cascade'][0]
cascade_comp_inst2 = np.load(cascade_compl2, allow_pickle=True)['cascade'][0]
bb_comp_inst = np.load(bb_complete, allow_pickle=True)['cascade'][0]
ml_comp_inst = np.load(ml_complete, allow_pickle=True)['cascade'][0]

b3d_cascade_inst = cascade_comp_inst.barrier3d[0]
b3d_cascade_inst2 = cascade_comp_inst2.barrier3d[0]
b3d_bmft_inst = np.load(b3d_bmft_complete, allow_pickle=True)['cascade'][0]

# ======================================================================================================================
# use the check init values function
# ======================================================================================================================
# comparing the cascade and barrierbmft classes
check_backbarrier_df = check_init_values(bmft_instance=bb_comp_inst,
                                             cascade_instance=cascade_comp_inst._barrierbmft[0].bmftc_BB)
check_mainland_df = check_init_values(bmft_instance=ml_comp_inst,
                                          cascade_instance=cascade_comp_inst._barrierbmft[0].bmftc_ML)
check_b3d_df = check_init_values(bmft_instance=b3d_bmft_inst, cascade_instance=b3d_cascade_inst)
# comparing the two cascade classes (before and after mocing the marsh_update() function to end of cascade.update()
check_bb_cascade_df = check_init_values(bmft_instance=cascade_comp_inst._barrierbmft[0].bmftc_BB,
                                             cascade_instance=cascade_comp_inst2._barrierbmft[0].bmftc_BB)
check_ml_cascade_df = check_init_values(bmft_instance=cascade_comp_inst._barrierbmft[0].bmftc_ML,
                                          cascade_instance=cascade_comp_inst2._barrierbmft[0].bmftc_ML)
check_b3d_cascade_df = check_init_values(bmft_instance=b3d_cascade_inst, cascade_instance=b3d_cascade_inst2)

# ======================================================================================================================
# save the dataframes as csv files
# ======================================================================================================================
check_backbarrier_df.to_csv(r"C:\Users\Lexi\Documents\UNC\BarrierBMFT\tests\comparison_results\backbarrier_check.csv", index=False)
check_mainland_df.to_csv(r"C:\Users\Lexi\Documents\UNC\BarrierBMFT\tests\comparison_results\mainland_check.csv", index=False)
check_b3d_df.to_csv(r"C:\Users\Lexi\Documents\UNC\BarrierBMFT\tests\comparison_results\b3d_check.csv", index=False)

check_bb_cascade_df.to_csv(r"C:\Users\Lexi\Documents\UNC\BarrierBMFT\tests\comparison_results\backbarrier_check_cascades.csv", index=False)
check_ml_cascade_df.to_csv(r"C:\Users\Lexi\Documents\UNC\BarrierBMFT\tests\comparison_results\mainland_check_cascades.csv", index=False)
check_b3d_cascade_df.to_csv(r"C:\Users\Lexi\Documents\UNC\BarrierBMFT\tests\comparison_results\b3d_check_cascades.csv", index=False)

# ======================================================================================================================
# plotting the data
# ======================================================================================================================

# since I could not save the barrierbmft class, we needed to save the self._x_b_TS_ML variable during the run so that
# we could use it here
xbts_ml = np.load(r"C:\Users\Lexi\Documents\UNC\BarrierBMFT\tests\xbTS_ML.npy", allow_pickle=True)

if plots_on:
    plot_2_cascade_transects(
        model1=cascade_comp_inst._barrierbmft[0],
        model2=cascade_comp_inst2._barrierbmft[0],
        save_dir=r"C:/Users/Lexi/Documents/UNC/BarrierBMFT/",
        newfolder="test_cascades"
    )

    plot_1_casc_1_barbmftc_transects(
        model1=cascade_comp_inst._barrierbmft[0],
        bmftc_ml=ml_comp_inst,
        bmftc_bb=bb_comp_inst,
        xbts_ml=xbts_ml,
        save_dir=r"C:/Users/Lexi/Documents/UNC/BarrierBMFT/tests",
        newfolder="edit2_comparison"
    )

    # plotting the domains for comparison
    elevFig1 = plt.figure(figsize=(7, 5))
    ax = elevFig1.add_subplot(111)
    cax = ax.pcolormesh(
        b3d_cascade_inst.DomainTS[-1]*10,
        cmap="terrain",
        vmin=-3,
        vmax=6,
        # edgecolors="w",  # for debugging
        # linewidth=0.01,
    )
    cbar = elevFig1.colorbar(cax)
    cbar.set_label("elevation (m MHW)", rotation=270)
    plt.xlabel("alongshore distance (dam)")
    plt.ylabel("cross-shore distance (dam)")
    timestr = (
            "Time = " + str(6) + " yrs"
    )  # we are letting the post-storm output represent 0.5 years
    elevFig1 = plt.figure(figsize=(7, 5))
    ax = elevFig1.add_subplot(111)
    cax = ax.pcolormesh(
        b3d_cascade_inst2.DomainTS[-1]*10,
        cmap="terrain",
        vmin=-3,
        vmax=6,
        # edgecolors="w",  # for debugging
        # linewidth=0.01,
    )
    cbar = elevFig1.colorbar(cax)
    cbar.set_label("elevation (m MHW)", rotation=270)
    plt.xlabel("alongshore distance (dam)")
    plt.ylabel("cross-shore distance (dam)")
    timestr = (
            "Time = " + str(6) + " yrs"
    )  # we are letting the post-storm output represent 0.5 years