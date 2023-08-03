# Lexi Van Blunk
# last updated 8/3/2023

# this code breaks up the "observed" configurations into interior, dune, and beach domains for input into CASCADE

# NOTE: configurations 1 and 2 vary slightly here because I was still using the domains from the thesis rather than
# the ones from the new larger DEM. In the future, I would suggest using these newer domains since they are consistent
# with configurations 3 and 4

import numpy as np

berm_el = 0.11
save_dir = "C:/Users/Lexi/PycharmProjects/CASCADE/scripts/outwash_ms/configurations/"
save_now = False

# ------------------------------------------------ load data -----------------------------------------------------------

section1 = np.load(r"C:\Users\Lexi\PycharmProjects\CASCADE\scripts\outwash_ms\configurations\config1\config1_observed_pre.npy")
section2 = np.load(r"C:\Users\Lexi\PycharmProjects\CASCADE\scripts\outwash_ms\configurations\config2\config2_observed_pre.npy")
section3 = np.load(r"C:\Users\Lexi\PycharmProjects\CASCADE\scripts\outwash_ms\configurations\config3\config3_observed_pre.npy")
section4 = np.load(r"C:\Users\Lexi\PycharmProjects\CASCADE\scripts\outwash_ms\configurations\config4\config4_observed_pre.npy")


# --------------------------------------------- section 1 --------------------------------------------------------------
section1_int = section1[0:30]
section1_dunes = section1[30:35]
section1_beach = section1[35:]

# take only two rows for the dunes because B3D only allows two dune rows
section1_dunes = section1_dunes[1:3]

full_1 = np.append(section1_int, section1_dunes, 0)
full_1 = np.append(full_1, section1_beach, 0)

if save_now:
    # save the domains as needed for B3D
    save_dir_1 = save_dir + "config1/"

    # flip the interior for correct orientation in B3D
    interior_b3d_input_1 = np.flip(section1_int)
    np.save(save_dir_1 + "NCB-default-elevation-config1-damMHW.npy", interior_b3d_input_1)

    # flip the dunes and put them into one long array
    dunes_b3d_1 = np.flip(section1_dunes) - berm_el
    dunes_input_1 = np.append(dunes_b3d_1[0], dunes_b3d_1[1], 0)
    np.save(save_dir_1 + "NCB-default-dunes-config1-dam.npy", dunes_input_1)

    # since the beach is not used in B3D, it can be saved normally (outwash orientation)
    np.save(save_dir_1 + "NCB-default-beach-config1-damMHW.npy", section1_beach)

    # save the full input just because
    np.save(save_dir_1 + "NCB-default_full_config1-damMHW.npy", full_1)


# --------------------------------------------- section 2 --------------------------------------------------------------
section2_int = section2[0:27]
section2_dunes = section2[27:31]
section2_beach = section2[31:]

# manipulate the dunes into 2 rows
section2_dunes[0:2, 0:4] += 0.1
section2_dunes[0:2, 0:2] = section2_dunes[3, 0:2]
section2_dunes[0, 2:5] = section2_dunes[2, 2:5]
section2_dunes = section2_dunes[0:2]
section2_dunes[0, 0:20] += 0.05
section2_dunes[1, 0:16] += 0.1

full_2 = np.append(section2_int, section2_dunes, 0)
full_2 = np.append(full_2, section2_beach, 0)

if save_now:
    # save the domains as needed for B3D
    save_dir_2 = save_dir + "config2/"

    # flip the interior for correct orientation in B3D
    interior_b3d_input_2 = np.flip(section2_int)
    np.save(save_dir_2 + "NCB-default-elevation-config2-damMHW.npy", interior_b3d_input_2)

    # flip the dunes and put them into one long array
    dunes_b3d_2 = np.flip(section2_dunes) - berm_el
    dunes_input_2 = np.append(dunes_b3d_2[0], dunes_b3d_2[1], 0)
    np.save(save_dir_2 + "NCB-default-dunes-config2-dam.npy", dunes_input_2)

    # since the beach is not used in B3D, it can be saved normally (outwash orientation)
    np.save(save_dir_2 + "NCB-default-beach-config2-damMHW.npy", section2_beach)

    # save the full input just because
    np.save(save_dir_2 + "NCB-default_full_config2-damMHW.npy", full_2)

# --------------------------------------------- section 3 --------------------------------------------------------------

section3_int = section3[:23, :]
section3_dunes = section3[23:28, :]
section3_beach = section3[28:, :]

# manipulate the dunes into 2 rows that represent the true dune line
section3_stitched_dunes = np.zeros([2, np.shape(section3_dunes)[1]])
section3_stitched_dunes[:, 0:4] = section3[26:28, 0:4]
section3_stitched_dunes[:, 4:9] = section3[25:27, 4:9]
section3_stitched_dunes[:, 9:22] = section3[24:26, 9:22]
section3_stitched_dunes[:, 22:] = section3[23:25, 22:]

# # filling in the beach under the dunes
extra_beach = np.zeros([3, np.shape(section3_dunes)[1]])
extra_beach[:, 22:] = section3[25:28, 22:]
extra_beach[:, 0:4] = section3[28:31, 0:4]
extra_beach[:, 4:9] = section3[27:30, 4:9]
extra_beach[:, 9:22] = section3[26:29, 9:22]
section3_beach[0:2, 0:9] = section3[30:32, 0:9]
section3_beach = np.append(extra_beach, section3_beach, 0)

full_3 = np.append(section3_int, section3_stitched_dunes, 0)
full_3 = np.append(full_3, section3_beach, 0)

if save_now:
    # save the domains as needed for B3D
    save_dir_3 = save_dir + "config3/"

    # flip the interior for correct orientation in B3D
    interior_b3d_input_3 = np.flip(section3_int)
    np.save(save_dir_3 + "NCB-default-elevation-config3-damMHW.npy", interior_b3d_input_3)

    # flip the dunes and put them into one long array
    dunes_b3d_3 = np.flip(section3_dunes) - berm_el
    dunes_input_3 = np.append(dunes_b3d_3[0], dunes_b3d_3[1], 0)
    np.save(save_dir_3 + "NCB-default-dunes-config3-dam.npy", dunes_input_3)

    # since the beach is not used in B3D, it can be saved normally (outwash orientation)
    np.save(save_dir_3 + "NCB-default-beach-config3-damMHW.npy", section3_beach)

    # save the full input just because
    np.save(save_dir_3 + "NCB-default_full_config3-damMHW.npy", full_3)

# --------------------------------------------- section 4 --------------------------------------------------------------

section4_int = section4[:26]
section4_dunes = section4[28:30]
section4_beach = section4[31:]
full_4 = np.append(section4_int, section4_dunes, 0)
full_4 = np.append(full_4, section4_beach, 0)

beach4_slope = (np.mean(section4_beach[0, :]) - np.mean(section4_beach[8, :])) / 9

if save_now:
    # save the domains as needed for B3D
    save_dir_4 = save_dir + "config4/"

    # flip the interior for correct orientation in B3D
    interior_b3d_input_4 = np.flip(section4_int)
    np.save(save_dir_4 + "NCB-default-elevation-config4-damMHW.npy", interior_b3d_input_4)

    # flip the dunes and put them into one long array
    dunes_b3d_4 = np.flip(section4_dunes) - berm_el
    dunes_input_4 = np.append(dunes_b3d_4[0], dunes_b3d_4[1], 0)
    np.save(save_dir_4 + "NCB-default-dunes-config4-dam.npy", dunes_input_4)

    # since the beach is not used in B3D, it can be saved normally (outwash orientation)
    np.save(save_dir_4 + "NCB-default-beach-config4-damMHW.npy", section4_beach)

    # save the full input just because
    np.save(save_dir_4 + "NCB-default_full_config4-damMHW.npy", full_4)