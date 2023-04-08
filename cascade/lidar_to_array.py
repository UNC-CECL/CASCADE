# conversion of Benton's Lidar to Barrier3D R code to Python
# some changes have been made to the original code to better suit North Core Banks
# Lexi Van Blunk
# 2/16/23

import numpy as np
import pandas as pd
import os
from matplotlib import pyplot as plt

plt.rcParams['figure.figsize'] = (8,6)
plt.rcParams.update({"font.size": 15})


# load the data as a dataframe
os.chdir('C:/Users/Lexi/Documents/Research/Outwasher/GIS')
prestorm_points_df = pd.read_csv("prestorm_points.csv")

# get the column that contains the elevations
prestorm_elevations = prestorm_points_df['grid_code']  # same as RawMatrix[,3]

num_cols_domain = 150  # cross-shore side of rectangle from GIS (barrier width)
num_rows_domain = int(len(prestorm_elevations)/num_cols_domain)  # this will be the max length of the island
MHW = 0.36  # meters NAVD88? (subtract 0.36 m from NAVD88 elevations to get MHW elevations)

prestorm_elevations = prestorm_elevations - MHW  # convert elevations from NAVD88 to MHW

# organizing the list of elevations into an array
elev_array = np.zeros([num_rows_domain, num_cols_domain])
start = 0
end = start + num_cols_domain
for row in range(num_rows_domain):
    elev_array[row, :] = prestorm_elevations[start:end]
    start += num_cols_domain
    end += num_cols_domain

# set all values less than 0 to -3 m
# tried range(135) also
for r in range(135):
    for c in range(num_cols_domain):
        if elev_array[r, c] <= 0:
            elev_array[r, c] = -3

# plot the adjusted domain
# fig1 = plt.figure()
# ax1 = fig1.add_subplot(111)
# mat = ax1.matshow(
#     elev_array,
#     cmap="terrain",
#     vmin=-3.0,
#     vmax=6.0,
# )
# cbar = fig1.colorbar(mat)
# cbar.set_label('m MHW', rotation=270, labelpad=15)
# ax1.set_title("Pre-Storm Elevation")
# ax1.set_ylabel("barrier width (dam)")
# ax1.set_xlabel("barrier length (dam)")
# plt.gca().xaxis.tick_bottom()

# convert to dam and make our two segments
elev_array = elev_array / 10
elev_array = elev_array[100:-1]
section1_pre = elev_array[:, 20:70]
section2_pre = elev_array[:, 70:120]

# # plot section 1
# fig1 = plt.figure()
# ax1 = fig1.add_subplot(111)
# mat = ax1.matshow(
#     section1_pre * 10,
#     cmap="terrain",
#     vmin=-3.0,
#     vmax=6.0,
# )
# cbar = fig1.colorbar(mat)
# cbar.set_label('m MHW', rotation=270, labelpad=15)
# ax1.set_title("Left Pre-Storm Elevation")
# ax1.set_ylabel("barrier width (dam)")
# ax1.set_xlabel("barrier length (dam)")
# plt.gca().xaxis.tick_bottom()

# # plot section 2
# fig1 = plt.figure()
# ax1 = fig1.add_subplot(111)
# mat = ax1.matshow(
#     section2_pre * 10,
#     cmap="terrain",
#     vmin=-3.0,
#     vmax=6.0,
# )
# cbar = fig1.colorbar(mat)
# cbar.set_label('m MHW', rotation=270, labelpad=15)
# ax1.set_title("Right Pre-Storm Elevation")
# ax1.set_ylabel("barrier width (dam)")
# ax1.set_xlabel("barrier length (dam)")
# plt.gca().xaxis.tick_bottom()


### plotting the post storm
post_points_df = pd.read_csv("post_storm_points.csv")

# get the column that contains the elevations
post_elevations = post_points_df['grid_code']  # same as RawMatrix[,3]

num_cols_domain = 150  # perpendicular side of rectangle
num_rows_domain = int(len(post_elevations)/num_cols_domain)  # this will be the max length of the island
post_elevations = post_elevations - MHW  # convert elevations from NAVD88 to MHW

# plot the entire domain out
post_elev_array = np.zeros([num_rows_domain, num_cols_domain])
start = 0
end = start + num_cols_domain
for row in range(num_rows_domain):
    post_elev_array[row, :] = post_elevations[start:end]
    start += num_cols_domain
    end += num_cols_domain

# set all values less than 0 to -3 m
for r in range(135):
    for c in range(num_cols_domain):
        if post_elev_array[r, c] <= 0:
            post_elev_array[r, c] = -3

# fig2 = plt.figure()
# ax2 = fig2.add_subplot(111)
# mat = ax2.matshow(
#     post_elev_array,
#     cmap="terrain",
#     vmin=-3.0,
#     vmax=6.0,
# )
# cbar = fig2.colorbar(mat)
# cbar.set_label('m MHW', rotation=270, labelpad=15)
# ax2.set_title("Post-Storm Elevation")
# ax2.set_ylabel("barrier width (dam)")
# ax2.set_xlabel("barrier length (dam)")
# plt.gca().xaxis.tick_bottom()


# convert to dam
post_elev_array = post_elev_array / 10
post_elev_array = post_elev_array[100:-1]
section1_post = post_elev_array[:, 19:69]
section2_post = post_elev_array[:, 69:119]

# # plot left post section
# fig2 = plt.figure()
# ax2 = fig2.add_subplot(111)
# mat = ax2.matshow(
#     # post_elev_array[120:-1]*10,
#     section1_post*10,
#     cmap="terrain",
#     vmin=-3.0,
#     vmax=6.0,
# )
# cbar = fig2.colorbar(mat)
# cbar.set_label('m MHW', rotation=270, labelpad=15)
# ax2.set_title("Left Post-Storm Elevation")
# ax2.set_ylabel("barrier width (dam)")
# ax2.set_xlabel("barrier length (dam)")
# plt.gca().xaxis.tick_bottom()
#
#
# # plot right post-section
# fig2 = plt.figure()
# ax2 = fig2.add_subplot(111)
# mat = ax2.matshow(
#     section2_post*10,
#     cmap="terrain",
#     vmin=-3.0,
#     vmax=6.0,
# )
# cbar = fig2.colorbar(mat)
# cbar.set_label('m MHW', rotation=270, labelpad=15)
# ax2.set_title("Right Post-Storm Elevation")
# ax2.set_ylabel("barrier width (dam)")
# ax2.set_xlabel("barrier length (dam)")
# plt.gca().xaxis.tick_bottom()

# re-formatting the left section
pre_left = section1_pre[37:78]
post_left = section1_post[36:77]

difference_left = post_left-pre_left
total_erosion_left = np.sum(difference_left)*10
back_erosion_left = np.sum(difference_left[3:30, 14:38])*10

# fig2 = plt.figure()
# ax2 = fig2.add_subplot(111)
# mat = ax2.matshow(
#     pre_left*10,
#     cmap="terrain",
#     vmin=-3.0,
#     vmax=6.0,
# )
# cbar = fig2.colorbar(mat)
# cbar.set_label('m MHW', rotation=270, labelpad=15)
# ax2.set_title("Left Zoom Pre-Storm Elevation")
# ax2.set_ylabel("barrier width (dam)")
# ax2.set_xlabel("barrier length (dam)")
# plt.gca().xaxis.tick_bottom()
#
# fig2 = plt.figure()
# ax2 = fig2.add_subplot(111)
# mat = ax2.matshow(
#     post_left*10,
#     cmap="terrain",
#     vmin=-3.0,
#     vmax=6.0,
# )
# cbar = fig2.colorbar(mat)
# cbar.set_label('m MHW', rotation=270, labelpad=15)
# ax2.set_title("Left Zoom Post-Storm Elevation")
# ax2.set_ylabel("barrier width (dam)")
# ax2.set_xlabel("barrier length (dam)")
# plt.gca().xaxis.tick_bottom()
#
# fig3 = plt.figure()
# ax3 = fig3.add_subplot(111)
# mat = ax3.matshow(
#     difference_left*10,
#     cmap="seismic",
#     vmin=-5,
#     vmax=5,
# )
# cbar = fig3.colorbar(mat)
# cbar.set_label('meters', rotation=270, labelpad=15)
# ax3.set_title("Difference Map")
# ax3.set_ylabel("barrier width (dam)")
# ax3.set_xlabel("barrier length (dam)")
# plt.gca().xaxis.tick_bottom()


# break up the domain
section1_pre_int = pre_left[0:30]
# section1_pre_int = pre_left[0:32]
# section1_pre_int[24:, 24:27] = 0.085
# section1_pre_int[27, 37] = 0.085
# section1_pre_int[28, 37] = 0.085
# section1_pre_int[29, 39] = 0.085
# section1_pre_int[28, 38:40] = 0.085
# section1_pre_dunes = pre_left[32:35]
section1_pre_dunes = pre_left[30:35]
section1_pre_dunes = section1_pre_dunes[1:3]
# section1_pre_dunes[0:2, 24] = 0.085
section1_pre_beach = pre_left[35:]
# section1_pre_beach2 = section1_pre[73:77]
# section1_pre_beach3 = section1_pre[73:78]

full = np.append(section1_pre_int, section1_pre_dunes, 0)
full = np.append(full, section1_pre_beach, 0)

fig2 = plt.figure()
ax2 = fig2.add_subplot(111)
mat = ax2.matshow(
    full*10,
    cmap="terrain",
    vmin=-3.0,
    vmax=6.0,
)
cbar = fig2.colorbar(mat)
cbar.set_label('m MHW', rotation=270, labelpad=15)
ax2.set_title("Pre-Storm Elevation")
ax2.set_ylabel("barrier width (dam)")
ax2.set_xlabel("barrier length (dam)")
plt.gca().xaxis.tick_bottom()


interior_b3d_input = np.flip(section1_pre_int)
# np.save(r"C:\Users\Lexi\PycharmProjects\CASCADE\cascade\data\outwash_data/NCB-default-elevation-left-damMHW-test", interior_b3d_input)

berm_el = 0.11
dunes_b3d = np.flip(section1_pre_dunes) - berm_el
dunes_input = np.append(dunes_b3d[0], dunes_b3d[1], 0)
# np.save(r"C:\Users\Lexi\PycharmProjects\CASCADE\cascade\data\outwash_data/NCB-default-dunes-left-dam-test", dunes_input)

# the beach should be loaded in with the ocean on the bottom
# np.save(r"C:\Users\Lexi\PycharmProjects\CASCADE\cascade\data\outwash_data/NCB-default-beach-left-dam-test", section1_pre_beach)

# re-formatting the right section
pre_right = section2_pre[36:76]
post_right = section2_post[35:75]

difference = post_right-pre_right
total_erosion_right = np.sum(difference)*10
back_erosion_right = np.sum(difference[4:28, 16:42])*10

fig2 = plt.figure()
ax2 = fig2.add_subplot(111)
mat = ax2.matshow(
    pre_right*10,
    cmap="terrain",
    vmin=-3.0,
    vmax=6.0,
)
cbar = fig2.colorbar(mat)
cbar.set_label('m MHW', rotation=270, labelpad=15)
ax2.set_title("Right Edited Pre-Storm Elevation")
ax2.set_ylabel("barrier width (dam)")
ax2.set_xlabel("barrier length (dam)")
plt.gca().xaxis.tick_bottom()

fig2 = plt.figure()
ax2 = fig2.add_subplot(111)
mat = ax2.matshow(
    post_right*10,
    cmap="terrain",
    vmin=-3.0,
    vmax=6.0,
)
cbar = fig2.colorbar(mat)
cbar.set_label('m MHW', rotation=270, labelpad=15)
ax2.set_title("Right Edited Post-Storm Elevation")
ax2.set_ylabel("barrier width (dam)")
ax2.set_xlabel("barrier length (dam)")
plt.gca().xaxis.tick_bottom()

fig3 = plt.figure()
ax3 = fig3.add_subplot(111)
mat = ax3.matshow(
    difference*10,
    cmap="seismic",
    vmin=-5,
    vmax=5,
)
cbar = fig3.colorbar(mat)
cbar.set_label('meters', rotation=270, labelpad=15)
ax3.set_title("Difference Map")
ax3.set_ylabel("barrier width (dam)")
ax3.set_xlabel("barrier length (dam)")
plt.gca().xaxis.tick_bottom()

section2_pre_int = pre_right[0:28]
section2_pre_dunes = pre_right[28:32]
section2_pre_dunes[0:2, 0:4] += 0.1
# section2_pre_dunes[0:2, 0:2] = section2_pre_dunes[3, 0:2]
# section2_pre_dunes[0, 2:5] = section2_pre_dunes[2, 2:5]
section2_pre_dunes = section2_pre_dunes[0:2]
section2_pre_dunes[0, 0:20] += 0.05
section2_pre_dunes[1, 0:16] += 0.1
# section2_pre_beach = pre[32:]
section2_pre_beach = np.load(r"C:\Users\Lexi\PycharmProjects\CASCADE\cascade\data\outwash_data\NCB-default-beach.npy")

full = np.append(section2_pre_int, section2_pre_dunes, 0)
full = np.append(full, section2_pre_beach, 0)

fig2 = plt.figure()
ax2 = fig2.add_subplot(111)
mat = ax2.matshow(
    full*10,
    cmap="terrain",
    vmin=-3.0,
    vmax=6.0,
)
cbar = fig2.colorbar(mat)
cbar.set_label('m MHW', rotation=270, labelpad=15)
ax2.set_title("Pre-Storm Elevation")
ax2.set_ylabel("barrier width (dam)")
ax2.set_xlabel("barrier length (dam)")
plt.gca().xaxis.tick_bottom()

interior_b3d_input = np.flip(section2_pre_int)
# np.save(r"C:\Users\Lexi\PycharmProjects\CASCADE\cascade\data\outwash_data/NCB-default-elevation-right-damMHW-test", interior_b3d_input)

berm_el = 0.11
dunes_b3d = np.flip(section2_pre_dunes) - berm_el
dunes_input = np.append(dunes_b3d[0], dunes_b3d[1], 0)
# np.save(r"C:\Users\Lexi\PycharmProjects\CASCADE\cascade\data\outwash_data/NCB-default-dunes-right-dam-test", dunes_input)

# the beach should be loaded in with the ocean on the bottom
# np.save(r"C:\Users\Lexi\PycharmProjects\CASCADE\cascade\data\outwash_data/NCB-default-beach-right-mod", section2_pre_beach)
