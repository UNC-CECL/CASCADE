# conversion of Benton's Lidar to Barrier3D R code to Python
# some changes have been made to the original code to better suit North Core Banks
# Lexi Van Blunk
# 2/16/23

import numpy as np
import pandas as pd
import os
from matplotlib import pyplot as plt

# load the data as a dataframe
os.chdir('C:/Users/Lexi/Documents/Research/ArcPro/')
points_df = pd.read_csv("prestorm_points_lidar.csv")

# get the column of elevations
elevations = points_df['grid_code']  # same as RawMatrix[,3]

num_cols_domain = 150  # perpendicular side of rectangle
num_rows_domain = int(len(elevations)/num_cols_domain)  # this will be the max length of the island
MHW = 0.46
elevations = elevations - MHW  # convert elevations from NAVD88 to MHW

# set all values less than 0 to -3 m
for index, val in enumerate(elevations):
    if val < -0.5:
        elevations[index] = -3

# plot the entire domain out
elev_array = np.zeros([num_rows_domain, num_cols_domain])
start = 0
end = start + num_cols_domain
for row in range(num_rows_domain):
    elev_array[row, :] = elevations[start:end] - 0.46
    start += num_cols_domain
    end += num_cols_domain

fig1 = plt.figure()
ax1 = fig1.add_subplot(111)
mat = ax1.matshow(
    elev_array,
    cmap="terrain",
    vmin=-3.0
    # vmax=3.0,
)
cbar = fig1.colorbar(mat)
cbar.set_label('m MHW', rotation=270, labelpad=15)
ax1.set_title("Initial Elevation")
ax1.set_ylabel("barrier width (dam)")
ax1.set_xlabel("barrier length (dam)")
plt.gca().xaxis.tick_bottom()

# set the desired island width and length
island_width = 60
island_length = 50

# initialize interior and dune array
interior_array = np.zeros([island_width, island_length])
dune_array = np.zeros(island_length)  # currently, this is only one row, but it could be 2-3

start = 0
end = start + num_cols_domain
for l in range(island_length):
    z = elevations[start:end]  # (one cross-shore row)
    revz = np.flip(z)  # flipped so that we start at the ocean side
    start_beach = np.min(np.where(revz > 0.46))  # finding the lowest index that is greater than 0.46
    end_beach = start_beach + 15  # at NCB, the beach + dunes are about 150 m long in some areas
    dunes = revz[start_beach:end_beach]
    dune_max_location = np.argmax(dunes)  # finding the index of the tallest dune
    dune_elev = np.max(dunes)  # finding the value of the tallest dune

    # start interior immediatly after duneline
    # I might change this: will talk to katherine about it
    start_island = dune_max_location + start_beach + 1
    end_island = start_island + island_width
    one_elev_col = revz[start_island:end_island]

    # Add values to larger matrix
    interior_array[:, island_length-l-1] = one_elev_col
    dune_array[island_length-l-1] = dune_elev

    # increment to the next row in the domain
    start += num_cols_domain
    end += num_cols_domain

# convert to dam
interior_array = interior_array / 10
dune_array= dune_array / 10

# plot the interior domain
fig2 = plt.figure()
ax2 = fig2.add_subplot(111)
mat = ax2.matshow(
    interior_array*10,
    cmap="terrain",
    vmin=-3.0
    # vmax=3.0,
)
cbar = fig2.colorbar(mat)
cbar.set_label('m MHW', rotation=270, labelpad=15)
ax2.set_title("Initial Elevation")
ax2.set_ylabel("barrier width (dam)")
ax2.set_xlabel("barrier length (dam)")
plt.gca().xaxis.tick_bottom()