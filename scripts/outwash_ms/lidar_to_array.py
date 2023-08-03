# Lexi Van Blunk
# last updated 8/3/2023

# this code loads the csv file of elevation points from ArcGIS and organizes them into python arrays
# we also break up the domains into the 4 "observed" configurations and save them to the
# scripts/outwash_ms/configurations folder

# SEE NOTEBOOK: plotting observed configurations, which does the same thing but visualizes each configuration into
# subplots displaying the prestorm, poststorm, and elevation change plots (and shows boxes where we calculate erosion)
# we also calculate erosion in the notebook

import numpy as np
import pandas as pd
import os
from matplotlib import pyplot as plt

plt.rcParams.update({"font.size": 15})


# --------------------------------------- load the data as a dataframe ------------------------------------------------

os.chdir('C:/Users/Lexi/Documents/Research/Outwasher Paper/GIS data/')

# PRESTORM DATA
prestorm_points_df = pd.read_csv("raster_points.csv")

# get the column from the csv that contains the elevations
prestorm_elevations = prestorm_points_df['grid_code']

alongshore_length_feature_class = 4000  # in meters, along-shore side of rectangle from GIS (barrier length)
num_cols_domain = int(alongshore_length_feature_class / 10)  # we specified 10x10 m cells
num_rows_domain = int(len(prestorm_elevations)/num_cols_domain)  # this is the width of the island
MHW = 0.36  # meters NAVD88? (subtract 0.36 m from NAVD88 elevations to get MHW elevations)

prestorm_elevations = prestorm_elevations - MHW  # convert elevations from NAVD88 to MHW

# POSTSTORM DATA
post_points_df = pd.read_csv("post_storm_points.csv")
post_elevations = post_points_df['grid_code']
post_elevations = post_elevations - MHW  # convert elevations from NAVD88 to MHW


# ---------------------------- organize the list of elevations into an array--------------------------------------------

elev_array = np.zeros([num_rows_domain, num_cols_domain])
start = 0
end = start + num_cols_domain
for row in range(num_rows_domain):
    elev_array[row, :] = prestorm_elevations[start:end]
    start += num_cols_domain
    end += num_cols_domain

post_elev_array = np.zeros([num_rows_domain, num_cols_domain])
start = 0
end = start + num_cols_domain
for row in range(num_rows_domain):
    post_elev_array[row, :] = post_elevations[start:end]
    start += num_cols_domain
    end += num_cols_domain

# set all values less than 0 to -3 m
for r in range(150):
    for c in range(num_cols_domain):
        if post_elev_array[r, c] <= -0.5:
            post_elev_array[r, c] = -3


# ------------------------------------------ plotting the domains ------------------------------------------------------
xlabel = "alongshore distance (m)"
ylabel = "cross-shore distance (m)"

# PRESTORM
fig1 = plt.figure()
ax1 = fig1.add_subplot(111)
mat = ax1.matshow(
    elev_array,
    cmap="terrain",
    vmin=-3.0,
    vmax=6.0,
)
cbar = fig1.colorbar(mat)
cbar.set_label('m MHW', rotation=270, labelpad=10)
ax1.set_title("Pre-Storm Elevation")
ax1.set_ylabel(ylabel)
ax1.set_xlabel(xlabel)
plt.gca().xaxis.tick_bottom()

# add the rectangles around the 4 configurations
# (x,y), length, width (up = neg)
ax1.add_patch(plt.Rectangle((98,187), 50, -40, lw=2, ec="k", fc="none"))            # black
ax1.add_patch(plt.Rectangle((148,187), 50, -40, lw=2, ec="m", fc="none"))           # purple
ax1.add_patch(plt.Rectangle((198,187), 50, -40, lw=2, ec="tab:orange", fc="none"))  # blue
ax1.add_patch(plt.Rectangle((248,187), 50, -40, lw=2, ec="tab:blue", fc="none"))    # orange

# edit the tick marks to be in meters
xtick_max = np.shape(elev_array)[1]  # n_cols = x
x_ticks = np.array(range(0, xtick_max, 50))
x_tick_labels = x_ticks * 10
ytick_max = np.shape(elev_array)[0]  # n_rows = y
y_ticks = np.array(range(0, ytick_max, 25))
y_tick_labels = y_ticks * 10
plt.xticks(x_ticks, x_tick_labels)
plt.yticks(y_ticks, y_tick_labels)
plt.tight_layout()

# POSTSTORM
fig2 = plt.figure()
ax2 = fig2.add_subplot(111)
mat = ax2.matshow(
    post_elev_array,
    cmap="terrain",
    vmin=-3.0,
    vmax=6.0,
)
cbar = fig2.colorbar(mat)
cbar.set_label('m MHW', rotation=270, labelpad=10)
ax2.set_title("Post-Storm Elevation")
ax2.set_ylabel(ylabel)
ax2.set_xlabel(xlabel)
plt.gca().xaxis.tick_bottom()

# add the rectangles around the 4 configurations
ax2.add_patch(plt.Rectangle((98,187), 50, -40, lw=2, ec="k", fc="none"))            # black
ax2.add_patch(plt.Rectangle((148,187), 50, -40, lw=2, ec="m", fc="none"))           # purple
ax2.add_patch(plt.Rectangle((198,187), 50, -40, lw=2, ec="tab:orange", fc="none"))  # blue
ax2.add_patch(plt.Rectangle((248,187), 50, -40, lw=2, ec="tab:blue", fc="none"))    # orange

# edit the tick marks to be in meters
xtick_max = np.shape(post_elev_array)[1]  # n_cols = x
x_ticks = np.array(range(0, xtick_max, 50))
x_tick_labels = x_ticks * 10
ytick_max = np.shape(post_elev_array)[0]  # n_rows = y
y_ticks = np.array(range(0, ytick_max, 25))
y_tick_labels = y_ticks * 10
plt.xticks(x_ticks, x_tick_labels)
plt.yticks(y_ticks, y_tick_labels)
plt.tight_layout()

# --------------------------- convert to dam and make the 4 configurations ---------------------------------------------
elev_array = elev_array / 10
section1_pre = elev_array[147:188, 98:148]
section2_pre = elev_array[147:188, 148:198]
section3_pre = elev_array[147:188, 248:298]
section4_pre = elev_array[147:188, 198:248]  # validation configuration

np.save(r"C:\Users\Lexi\PycharmProjects\CASCADE\scripts\outwash_ms\configurations\config1_observed_pre.npy", section1_pre)
np.save(r"C:\Users\Lexi\PycharmProjects\CASCADE\scripts\outwash_ms\configurations\config2_observed_pre.npy", section2_pre)
np.save(r"C:\Users\Lexi\PycharmProjects\CASCADE\scripts\outwash_ms\configurations\config3_observed_pre.npy", section3_pre)
np.save(r"C:\Users\Lexi\PycharmProjects\CASCADE\scripts\outwash_ms\configurations\config4_observed_pre.npy", section4_pre)

post_elev_array = post_elev_array / 10
section1_post = post_elev_array[147:188, 98:148]
section2_post = post_elev_array[147:188, 148:198]
section3_post = post_elev_array[147:188, 248:298]
section4_post = post_elev_array[147:188, 198:248]

np.save(r"C:\Users\Lexi\PycharmProjects\CASCADE\scripts\outwash_ms\configurations\config1_observed_post.npy", section1_post)
np.save(r"C:\Users\Lexi\PycharmProjects\CASCADE\scripts\outwash_ms\configurations\config2_observed_post.npy", section2_post)
np.save(r"C:\Users\Lexi\PycharmProjects\CASCADE\scripts\outwash_ms\configurations\config3_observed_post.npy", section3_post)
np.save(r"C:\Users\Lexi\PycharmProjects\CASCADE\scripts\outwash_ms\configurations\config4_observed_post.npy", section4_post)
