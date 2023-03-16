# conversion of Benton's Lidar to Barrier3D R code to Python
# some changes have been made to the original code to better suit North Core Banks
# Lexi Van Blunk
# 2/16/23

import numpy as np
import pandas as pd
import os
from matplotlib import pyplot as plt
from collections import Counter

def marsh_averages(
        elev_array,
        int_width,
        length,
        block_size):

    # adjusting for a block size that does not divide evenly into the domain
    b_size = block_size
    extra_vert_cells = int_width % b_size  # gives the extra row cells
    extra_lat_cells = length % b_size  # gives the extra column cells

    # calculating how many times we will shift the block to calculate blocks of average slopes
    if extra_vert_cells == 0:
        n_shifts_vert = int(int_width / b_size)
    else:
        n_shifts_vert = int((int_width - extra_vert_cells) / b_size) + 1
    if extra_lat_cells == 0:
        n_shifts_lat = int(length / b_size)
    else:
        n_shifts_lat = int((length - extra_lat_cells) / b_size) + 1

    # Shift through the entire column and then move block up to the next row
    # start right before the dune line

    for v in range(n_shifts_vert):
        if v == n_shifts_vert-1 and extra_vert_cells != 0:  # this is the last shift
            start_row = stop_row
            stop_row = start_row + extra_vert_cells
        else:
            start_row = v * b_size
            stop_row = v * b_size + b_size
        for l in range(n_shifts_lat):
            if l == n_shifts_lat-1 and extra_lat_cells != 0:
                start_col = end_col
                end_col = start_col + extra_lat_cells + 1
            else:
                start_col = l*b_size
                end_col = l*b_size + b_size
            S = np.mean(elev_array[start_row:(stop_row+1), start_col:end_col])
            if S < -0.1:
                elev_array[start_row:(stop_row+1), start_col:end_col] = -0.3

    return elev_array

# load the data as a dataframe
os.chdir('C:/Users/Lexi/Documents/Research/Outwasher/GIS')
prestorm_points_df = pd.read_csv("prestorm_points.csv")

# get the column of elevations
prestorm_elevations = prestorm_points_df['grid_code']  # same as RawMatrix[,3]

num_cols_domain = 150  # perpendicular side of rectangle
num_rows_domain = int(len(prestorm_elevations)/num_cols_domain)  # this will be the max length of the island
MHW = 0.36  # meters

prestorm_elevations = prestorm_elevations - MHW  # convert elevations from NAVD88 to MHW

# plot the entire domain out
elev_array = np.zeros([num_rows_domain, num_cols_domain])
start = 0
end = start + num_cols_domain
for row in range(num_rows_domain):
    elev_array[row, :] = prestorm_elevations[start:end]
    start += num_cols_domain
    end += num_cols_domain

# set all values less than 0 to -3 m
for r in range(135):
    for c in range(num_cols_domain):
        if elev_array[r, c] <= 0:
            elev_array[r, c] = -3

# convert to dam
elev_array = elev_array / 10
elev_array = elev_array[100:-1]
section1_pre = elev_array[:, 20:70]
section2_pre = elev_array[:, 70:120]

fig1 = plt.figure()
ax1 = fig1.add_subplot(111)
mat = ax1.matshow(
    section1_pre * 10,
    cmap="terrain",
    vmin=-3.0,
    vmax=6.0,
)
cbar = fig1.colorbar(mat)
cbar.set_label('m MHW', rotation=270, labelpad=15)
ax1.set_title("Pre-Storm Elevation")
ax1.set_ylabel("barrier width (dam)")
ax1.set_xlabel("barrier length (dam)")
plt.gca().xaxis.tick_bottom()

### plotting the post storm
post_points_df = pd.read_csv("post_storm_points.csv")

# get the column of elevations
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

# for r in range(175, num_rows_domain):
#     for c in range(num_cols_domain):
#         if post_elev_array[r, c] <= 0:
#             post_elev_array[r, c] = -3

# convert to dam
post_elev_array = post_elev_array / 10
post_elev_array = post_elev_array[100:-1]
section1_post = post_elev_array[:, 19:69]
section2_post = post_elev_array[:, 69:119]


fig2 = plt.figure()
ax2 = fig2.add_subplot(111)
mat = ax2.matshow(
    # post_elev_array[120:-1]*10,
    section1_post*10,
    cmap="terrain",
    vmin=-3.0,
    vmax=6.0,
)
cbar = fig2.colorbar(mat)
cbar.set_label('m MHW', rotation=270, labelpad=15)
ax2.set_title("Post-Storm Elevation")
ax2.set_ylabel("barrier width (dam)")
ax2.set_xlabel("barrier length (dam)")
plt.gca().xaxis.tick_bottom()

pre = section1_pre[37:76]
post = section1_post[36:75]
difference = post-pre
total_erosion = np.sum(difference)

fig2 = plt.figure()
ax2 = fig2.add_subplot(111)
mat = ax2.matshow(
    pre*10,
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

fig2 = plt.figure()
ax2 = fig2.add_subplot(111)
mat = ax2.matshow(
    post*10,
    cmap="terrain",
    vmin=-3.0,
    vmax=6.0,
)
cbar = fig2.colorbar(mat)
cbar.set_label('m MHW', rotation=270, labelpad=15)
ax2.set_title("Post-Storm Elevation")
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

## break up the domain
section1_pre_int = pre[0:30]
section1_pre_int[24:, 24:27] = 0.085
section1_pre_int[27, 37] = 0.085
section1_pre_int[28, 37] = 0.085
section1_pre_int[29, 39] = 0.085
section1_pre_int[28, 38:40] = 0.085
section1_pre_dunes = pre[30:35]
section1_pre_dunes = section1_pre_dunes[1:3]
section1_pre_dunes[0:2, 24] = 0.085
section1_pre_beach = pre[35:]
section1_pre_beach2 = section1_pre[73:77]
section1_pre_beach3 = section1_pre[73:78]

full = np.append(section1_pre_int, section1_pre_dunes, 0)
full = np.append(full, section1_pre_beach3, 0)

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
np.save(r"C:\Users\Lexi\PycharmProjects\CASCADE\cascade\data\outwash_data/NCB-default-elevation_mod", interior_b3d_input)

# berm_el = np.mean(section1_pre_beach[0:2])
berm_el = 0.11
# max_dunes = np.append(section1_pre_dunes[3, 0:24], section1_pre_dunes[1, 24:])
# dunes_b3d = np.flip(max_dunes) - berm_el
dunes_b3d = np.flip(section1_pre_dunes) - berm_el
dunes_input = np.append(dunes_b3d[0], dunes_b3d[1], 0)
# np.save(r"C:\Users\Lexi\PycharmProjects\CASCADE\cascade\data\outwash_data/NCB-default-dunes2", dunes_input)
#
# # the beach should be loaded in with the ocean on the bottom
# np.save(r"C:\Users\Lexi\PycharmProjects\CASCADE\cascade\data\outwash_data/NCB-default-beach4", section1_pre_beach3)
# beach_slope = np.mean(section1_pre_beach[0, :] - section1_pre_beach[-1, :])/len(section1_pre_beach)