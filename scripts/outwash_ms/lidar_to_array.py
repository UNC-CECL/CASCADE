# conversion of Benton's Lidar to Barrier3D R code to Python
# some changes have been made to the original code to better suit North Core Banks
# Lexi Van Blunk
# 2/16/23, updated 4/8/23
# import matplotlib.patches
import numpy as np
import pandas as pd
import os
from matplotlib import pyplot as plt
# import heapq

plt.rcParams.update({"font.size": 15})
# plt.rcParams['figure.constrained_layout.use'] = True


# ### load the data as a dataframe--------------------------------------------------------------------------------------

# os.chdir('C:/Users/Lexi/Documents/Research/Outwasher/GIS')
# prestorm_points_df = pd.read_csv("prestorm_points.csv")  # this data is in NAVD88

os.chdir('C:/Users/Lexi/Documents/Research/Outwasher Paper/GIS data/')
prestorm_points_df = pd.read_csv("raster_points.csv")

# get the column that contains the elevations
prestorm_elevations = prestorm_points_df['grid_code']  # same as RawMatrix[,3]

alongshore_length_feature_class = 4000
# alongshore_length_feature_class = 1500  # in meters, along-shore side of rectangle from GIS (barrier length)
num_cols_domain = int(alongshore_length_feature_class / 10)  # we specified 10x10 m cells
num_rows_domain = int(len(prestorm_elevations)/num_cols_domain)  # this will be the width of the island
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

# # set all values less than 0 to -3 m
# # tried range(135) also
# for r in range(150):
#     for c in range(num_cols_domain):
#         if elev_array[r, c] <= 0:
#             elev_array[r, c] = -3

# plot the adjusted domain
fig1 = plt.figure()
ax1 = fig1.add_subplot(111)
# ax1 = fig1.add_subplot(121)
mat = ax1.matshow(
    elev_array,
    cmap="terrain",
    vmin=-3.0,
    vmax=6.0,
)
cbar = fig1.colorbar(mat)
cbar.set_label('m MHW', rotation=270, labelpad=10)
ax1.set_title("Pre-Storm Elevation")
ax1.set_ylabel("barrier width (m)")
ax1.set_xlabel("barrier length (m)")
plt.gca().xaxis.tick_bottom()

# ax1.add_patch(plt.Rectangle((20,175), 50, -40, lw=2, ec="k", fc="none"))  # (x,y), length, width (up = neg)
# ax1.add_patch(plt.Rectangle((70,175), 50, -40, lw=2, ec="m", fc="none"))
ax1.add_patch(plt.Rectangle((98,187), 50, -40, lw=2, ec="k", fc="none"))  # (x,y), length, width (up = neg)
ax1.add_patch(plt.Rectangle((148,187), 50, -40, lw=2, ec="m", fc="none"))
ax1.add_patch(plt.Rectangle((198,187), 50, -40, lw=2, ec="tab:orange", fc="none"))
ax1.add_patch(plt.Rectangle((248,187), 50, -40, lw=2, ec="tab:blue", fc="none"))
# ls = linestyle, lw = line width, ec = edge color, fc = fill color

xtick_max = np.shape(elev_array)[1]  # n_cols = x
x_ticks = np.array(range(0, xtick_max, 25))
x_tick_labels = x_ticks * 10
ytick_max = np.shape(elev_array)[0]  # n_rows = y
y_ticks = np.array(range(0, ytick_max, 25))
y_tick_labels = y_ticks * 10
plt.xticks(x_ticks, x_tick_labels)
plt.yticks(y_ticks, y_tick_labels)
fig1.subplots_adjust(wspace=0.25)

# print(np.mean(elev_array[168:171, 0:25]))
#
# convert to dam and make our two segments
elev_array = elev_array / 10
# elev_array = elev_array[100:-1]
# section1_pre = elev_array[:, 20:70]
# section2_pre = elev_array[:, 70:120]
# section3_pre = elev_array[147:188, 248:298]
section4_pre = elev_array[147:188, 198:248]

# print(np.max(elev_array))

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

# plot section 3
# fig1 = plt.figure()
# ax1 = fig1.add_subplot(111)
# mat = ax1.matshow(
#     section3_pre * 10,
#     cmap="terrain",
#     vmin=-3.0,
#     vmax=6.0,
# )
# cbar = fig1.colorbar(mat)
# cbar.set_label('m MHW', rotation=270, labelpad=15)
# ax1.set_title("Pre-Storm Elevation")
# ax1.set_ylabel("barrier width (m)")
# ax1.set_xlabel("barrier length (m)")
# plt.gca().xaxis.tick_bottom()
#
# plt.hlines(21.5, -0.5, 49.5, color="k", linestyles='solid', linewidth=3)
#
# xtick_max = np.shape(section3_pre)[1]  # n_cols = x
# x_ticks = np.array(range(0, xtick_max, 5))
# x_tick_labels = x_ticks * 10
# ytick_max = np.shape(section3_pre)[0]  # n_rows = y
# y_ticks = np.array(range(0, ytick_max, 5))
# y_tick_labels = y_ticks * 10
# plt.xticks(x_ticks, x_tick_labels)
# plt.yticks(y_ticks, y_tick_labels)

# plot section 4
fig1 = plt.figure()
ax1 = fig1.add_subplot(111)
mat = ax1.matshow(
    section4_pre * 10,
    cmap="terrain",
    vmin=-3.0,
    vmax=6.0,
)
cbar = fig1.colorbar(mat)
cbar.set_label('m MHW', rotation=270, labelpad=15)
ax1.set_title("Pre-Storm Elevation")
ax1.set_ylabel("barrier width (m)")
ax1.set_xlabel("barrier length (m)")
plt.gca().xaxis.tick_bottom()

# plt.hlines(21.5, -0.5, 49.5, color="k", linestyles='solid', linewidth=3)

xtick_max = np.shape(section4_pre)[1]  # n_cols = x
x_ticks = np.array(range(0, xtick_max, 5))
x_tick_labels = x_ticks * 10
ytick_max = np.shape(section4_pre)[0]  # n_rows = y
y_ticks = np.array(range(0, ytick_max, 5))
y_tick_labels = y_ticks * 10
plt.xticks(x_ticks, x_tick_labels)
plt.yticks(y_ticks, y_tick_labels)


### plotting the post storm
post_points_df = pd.read_csv("post_storm_points.csv")
#
# get the column that contains the elevations
post_elevations = post_points_df['grid_code']  # same as RawMatrix[,3]

# num_cols_domain = 150  # perpendicular side of rectangle
# num_rows_domain = int(len(post_elevations)/num_cols_domain)  # this will be the max length of the island
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
for r in range(150):
    for c in range(num_cols_domain):
        if post_elev_array[r, c] <= -0.5:
            post_elev_array[r, c] = -3

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
ax2.set_ylabel("barrier width (m)")
ax2.set_xlabel("barrier length (m)")
plt.gca().xaxis.tick_bottom()
# ax2.add_patch(plt.Rectangle((20,175), 50, -40, lw=2, ec="k", fc="none"))
# ax2.add_patch(plt.Rectangle((70,175), 50, -40, lw=2, ec="m", fc="none"))
ax2.add_patch(plt.Rectangle((98,187), 50, -40, lw=2, ec="k", fc="none"))  # (x,y), length, width (up = neg)
ax2.add_patch(plt.Rectangle((148,187), 50, -40, lw=2, ec="m", fc="none"))
ax2.add_patch(plt.Rectangle((198,187), 50, -40, lw=2, ec="tab:orange", fc="none"))
ax2.add_patch(plt.Rectangle((248,187), 50, -40, lw=2, ec="tab:blue", fc="none"))

xtick_max = np.shape(post_elev_array)[1]  # n_cols = x
x_ticks = np.array(range(0, xtick_max, 25))
x_tick_labels = x_ticks * 10
ytick_max = np.shape(post_elev_array)[0]  # n_rows = y
y_ticks = np.array(range(0, ytick_max, 25))
y_tick_labels = y_ticks * 10
plt.xticks(x_ticks, x_tick_labels)
plt.yticks(y_ticks, y_tick_labels)

# ax2.tick_params(left=False)
# ax2.set(yticklabels=[])
# plt.subplots_adjust(hspace=2)

# convert to dam
post_elev_array = post_elev_array / 10
# post_elev_array = post_elev_array[100:-1]
# section1_post = post_elev_array[:, 19:69]
# section2_post = post_elev_array[:, 69:119]
# section3_post = post_elev_array[147:188, 248:298]
section4_post = post_elev_array[147:188, 198:248]
#
# plot post section
fig2 = plt.figure()
ax2 = fig2.add_subplot(111)
mat = ax2.matshow(
    # post_elev_array[120:-1]*10,
    section4_post*10,
    cmap="terrain",
    vmin=-3.0,
    vmax=6.0,
)
cbar = fig2.colorbar(mat)
cbar.set_label('m MHW', rotation=270, labelpad=15)
ax2.set_title("Post-Storm Elevation")
ax2.set_ylabel("barrier width (m)")
ax2.set_xlabel("barrier length (m)")
plt.gca().xaxis.tick_bottom()

# plt.hlines(21.5, -0.5, 49.5, color="k", linestyles='solid', linewidth=3)

xtick_max = np.shape(section4_post)[1]  # n_cols = x
x_ticks = np.array(range(0, xtick_max, 5))
x_tick_labels = x_ticks * 10
ytick_max = np.shape(section4_post)[0]  # n_rows = y
y_ticks = np.array(range(0, ytick_max, 5))
y_tick_labels = y_ticks * 10
plt.xticks(x_ticks, x_tick_labels)
plt.yticks(y_ticks, y_tick_labels)
# fig2.subplots_adjust(wspace=0.25)

dif4 = section4_post - section4_pre
total_erosion_4 = np.sum(dif4)*10*100  # cubic meters
left_channel_erosion4 = np.sum(dif4[5:29, 2:21])*10*100  # cubic meters
right_channel_erosion4 = np.sum(dif4[10:29, 22:39])*10*100  # cubic meters
one_box_erosion4 = np.sum(dif4[5:29, 2:39])*10*100  # cubic meters
# back_barrier_erosion_3 = np.sum(dif3[0:23, :])*10*100  # cubic meters
# back_barrier_erosion_3_withdunes = np.sum(dif3[0:25, :])*10*100  # cubic meters
# better representation of back-barrier and dune line erosion
# back_barrier_erosion_3 = np.sum(dif3[0:22, :])*10*100  # cubic meters
# back_barrier_erosion_3_withdunes = np.sum(dif3[0:29, :])*10*100  # cubic meters

# plot difference
fig2 = plt.figure()
ax2 = fig2.add_subplot(111)
mat = ax2.matshow(
    # post_elev_array[120:-1]*10,
    dif4*10,
    cmap="seismic",
    vmin=-5,
    vmax=5,
)
cbar = fig2.colorbar(mat)
cbar.set_label('m MHW', rotation=270, labelpad=15)
ax2.set_title("Elevation Change")
ax2.set_ylabel("barrier width (m)")
ax2.set_xlabel("barrier length (m)")
plt.gca().xaxis.tick_bottom()
# plt.hlines(21.5, -0.5, 49.5, color="k", linestyles='solid', linewidth=3)
ax2.add_patch(plt.Rectangle((1, 29), 38, -25, lw=3, ec="k", fc="none"))
# ax2.add_patch(plt.Rectangle((1,29), 19, -25, lw=3, ec="k", fc="none"))
# ax2.add_patch(plt.Rectangle((21,29), 18, -20, lw=3, ec="k", fc="none"))

xtick_max = np.shape(dif4)[1]  # n_cols = x
x_ticks = np.array(range(0, xtick_max, 5))
x_tick_labels = x_ticks * 10
ytick_max = np.shape(dif4)[0]  # n_rows = y
y_ticks = np.array(range(0, ytick_max, 5))
y_tick_labels = y_ticks * 10
plt.xticks(x_ticks, x_tick_labels)
plt.yticks(y_ticks, y_tick_labels)


# ### breaking up the domain -------------------------------------------------------------------------------------------

# break up the domain into interior, dunes, and beach
# will likely take some manipulating with the dunes
section4_int = section4_pre[:26]
section4_beach = section4_pre[31:]
# section4_dunes = section4_pre[26:31]
section4_dunes = section4_pre[28:30]
full4 = np.append(section4_int, section4_dunes, 0)
full4 = np.append(full4, section4_beach, 0)

# section3_int = section3_pre[:23, :]
# section3_dunes = section3_pre[23:28, :]
# # section3_new_dunes = np.zeros([2, np.shape(section3_dunes)[1]])
# section3_beach = section3_pre[28:, :]
#
# extra_beach = np.zeros([3, np.shape(section3_dunes)[1]])
# extra_beach[:, 22:end] = section3_pre[25:28, 22:end]
# extra_beach[:, 0:4] = section3_pre[28:31, 0:4]
# extra_beach[:, 4:9] = section3_pre[27:30, 4:9]
# extra_beach[:, 9:22] = section3_pre[26:29, 9:22]
#
# section3_beach[0:2, 0:9] = section3_pre[30:32, 0:9]
#
# section3_beach = np.append(extra_beach, section3_beach, 0)
#
# # # dune manipulation
# # for c in range(np.shape(section3_dunes)[1]):
# #     max_dune1, max_dune2 = heapq.nlargest(2, section3_dunes[:, c])
# #     section3_new_dunes[0, c] = max_dune2
# #     section3_new_dunes[1, c] = max_dune1
#
# # second option
# section3_stitched_dunes = np.zeros([2, np.shape(section3_dunes)[1]])
# section3_stitched_dunes[:, 0:4] = section3_pre[26:28, 0:4]
# section3_stitched_dunes[:, 4:9] = section3_pre[25:27, 4:9]
# section3_stitched_dunes[:, 9:22] = section3_pre[24:26, 9:22]
# section3_stitched_dunes[:, 22:] = section3_pre[23:25, 22:]
#
# # full3_option1 = np.append(section3_int, section3_new_dunes, 0)
# # full3_option1 = np.append(full3_option1, section3_beach, 0)
#
# full3_option2 = np.append(section3_int, section3_stitched_dunes, 0)
# full3_option2 = np.append(full3_option2, section3_beach, 0)

fig1 = plt.figure()
ax1 = fig1.add_subplot(111)
mat = ax1.matshow(
    full4 * 10,
    cmap="terrain",
    vmin=-3.0,
    vmax=6.0,
)
cbar = fig1.colorbar(mat)
cbar.set_label('m MHW', rotation=270, labelpad=15)
ax1.set_title("Pre-Storm Elevation: Model Input")
ax1.set_ylabel("barrier width (m)")
ax1.set_xlabel("barrier length (m)")
plt.gca().xaxis.tick_bottom()

# plt.hlines(25.5, -0.5, 49.5, color="k", linestyles='solid', linewidth=3)
# plt.hlines(27.5, -0.5, 49.5, color="k", linestyles='solid', linewidth=3)

xtick_max = np.shape(full4)[1]  # n_cols = x
x_ticks = np.array(range(0, xtick_max, 5))
x_tick_labels = x_ticks * 10
ytick_max = np.shape(full4)[0]  # n_rows = y
y_ticks = np.array(range(0, ytick_max, 5))
y_tick_labels = y_ticks * 10
plt.xticks(x_ticks, x_tick_labels)
plt.yticks(y_ticks, y_tick_labels)

beach4_slope = (np.mean(section4_beach[0, :]) - np.mean(section4_beach[8, :])) / 9


# ### save the domains -------------------------------------------------------------------------------------------------
# interior_b3d_input = np.flip(section4_int)
# np.save(r"C:\Users\Lexi\PycharmProjects\CASCADE\cascade\data\outwash_data/NCB-default-elevation-config4-damMHW", interior_b3d_input)
#
# berm_el = 0.11
# dunes_b3d = np.flip(section4_dunes) - berm_el
# dunes_input = np.append(dunes_b3d[0], dunes_b3d[1], 0)
# np.save(r"C:\Users\Lexi\PycharmProjects\CASCADE\cascade\data\outwash_data/NCB-default-dunes-config4-dam", dunes_input)
#
# # the beach should be loaded in with the ocean on the bottom
# np.save(r"C:\Users\Lexi\PycharmProjects\CASCADE\cascade\data\outwash_data/NCB-default-beach-config4-damMHW", section4_beach)


# ### save the domains -------------------------------------------------------------------------------------------------
# interior_b3d_input = np.flip(section3_int)
# np.save(r"C:\Users\Lexi\PycharmProjects\CASCADE\cascade\data\outwash_data/NCB-default-elevation-config3-damMHW", interior_b3d_input)
#
# berm_el = 0.11
# dunes_b3d = np.flip(section3_stitched_dunes) - berm_el
# dunes_input = np.append(dunes_b3d[0], dunes_b3d[1], 0)
# np.save(r"C:\Users\Lexi\PycharmProjects\CASCADE\cascade\data\outwash_data/NCB-default-dunes-config3-dam", dunes_input)
#
# # the beach should be loaded in with the ocean on the bottom
# np.save(r"C:\Users\Lexi\PycharmProjects\CASCADE\cascade\data\outwash_data/NCB-default-beach-config3-damMHW", section3_beach)



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

# # re-formatting the left section
# pre_left = section1_pre[37:78]
# post_left = section1_post[36:77]
#
# difference_left = post_left-pre_left
# total_erosion_left = np.sum(difference_left)*10
# back_erosion_left = np.sum(difference_left[3:30, 14:38])*10
#
# # shrink=0.4
# #
# # fig1 = plt.figure()
# # fig1.tight_layout()
# # fig1.suptitle('Configuration 1: Observed Erosion', weight="bold")
# # ax1 = fig1.add_subplot(131)
# # mat = ax1.matshow(
# #     pre_left*10,
# #     cmap="terrain",
# #     vmin=-3.0,
# #     vmax=6.0,
# # )
# # cbar = fig1.colorbar(mat, shrink=shrink)
# # cbar.set_label('m MHW', rotation=270, labelpad=5)
# # ax1.set_title("Pre-Storm Elevation")
# # ax1.set_ylabel("barrier width (m)")
# # ax1.set_xlabel("barrier length (m)")
# # plt.gca().xaxis.tick_bottom()
# # xtick_max = np.shape(pre_left)[1]  # n_cols = x
# # x_ticks = np.array(range(0, xtick_max, 10))
# # x_tick_labels = x_ticks * 10
# # ytick_max = np.shape(pre_left)[0]  # n_rows = y
# # y_ticks = np.array(range(0, ytick_max, 10))
# # y_tick_labels = y_ticks * 10
# # plt.xticks(x_ticks, x_tick_labels)
# # plt.yticks(y_ticks, y_tick_labels)
#
# #
# # ax2 = fig1.add_subplot(132)
# # mat = ax2.matshow(
# #     post_left*10,
# #     cmap="terrain",
# #     vmin=-3.0,
# #     vmax=6.0,
# # )
# # cbar = fig1.colorbar(mat, shrink=shrink)
# # cbar.set_label('m MHW', rotation=270, labelpad=5)
# # ax2.set_title("Post-Storm Elevation")
# # ax2.set_ylabel("barrier width (m)")
# # ax2.set_xlabel("barrier length (m)")
# # plt.gca().xaxis.tick_bottom()
# # xtick_max = np.shape(post_left)[1]  # n_cols = x
# # x_ticks = np.array(range(0, xtick_max, 10))
# # x_tick_labels = x_ticks * 10
# # ytick_max = np.shape(post_left)[0]  # n_rows = y
# # y_ticks = np.array(range(0, ytick_max, 10))
# # y_tick_labels = y_ticks * 10
# # plt.xticks(x_ticks, x_tick_labels)
# # plt.yticks(y_ticks, y_tick_labels)
# # # ax2.tick_params(left=False)
# # # ax2.set(yticklabels=[])
# # # fig2.subplots_adjust(wspace=0.3, hspace=0)
# #
# # # fig3 = plt.figure()
# # ax3 = fig1.add_subplot(133)
# # mat = ax3.matshow(
# #     difference_left*10,
# #     cmap="seismic",
# #     vmin=-5,
# #     vmax=5,
# # )
# # # Configuration 1 \n
# # cbar = fig1.colorbar(mat, shrink=shrink)
# # cbar.set_label('meters', rotation=270, labelpad=5)
# # ax3.set_title("Elevation Change")
# # ax3.set_ylabel("barrier width (m)")
# # ax3.set_xlabel("barrier length (m)")
# # plt.gca().xaxis.tick_bottom()
# #
# # ax3.add_patch(plt.Rectangle((13,30), 24, -27, lw=3, ec="k", fc="none"))
# #
# # xtick_max = np.shape(difference_left)[1]  # n_cols = x
# # x_ticks = np.array(range(0, xtick_max, 10))
# # x_tick_labels = x_ticks * 10
# # ytick_max = np.shape(difference_left)[0]  # n_rows = y
# # y_ticks = np.array(range(0, ytick_max, 10))
# # y_tick_labels = y_ticks * 10
# # plt.xticks(x_ticks, x_tick_labels)
# # plt.yticks(y_ticks, y_tick_labels)
# # fig1.subplots_adjust(top=1.25, wspace=0.4, hspace=0)


# # break up the domain --------------------------------------------------------------------------------------------------
# section1_pre_int = pre_left[0:30]
# # section1_pre_int = pre_left[0:32]
# # section1_pre_int[24:, 24:27] = 0.085
# # section1_pre_int[27, 37] = 0.085
# # section1_pre_int[28, 37] = 0.085
# # section1_pre_int[29, 39] = 0.085
# # section1_pre_int[28, 38:40] = 0.085
# # section1_pre_dunes = pre_left[32:35]
# section1_pre_dunes = pre_left[30:35]
#
# section1_pre_dunes = section1_pre_dunes[1:3]
# # section1_pre_dunes[0:2, 24] = 0.085
# section1_pre_beach = pre_left[35:]
# # section1_pre_beach2 = section1_pre[73:77]
# # section1_pre_beach3 = section1_pre[73:78]
#
# full = np.append(section1_pre_int, section1_pre_dunes, 0)
# full = np.append(full, section1_pre_beach, 0)
#
# # fig2 = plt.figure()
# # fig2.set_tight_layout
# # ax1 = fig2.add_subplot(121)
# # mat = ax1.matshow(
# #     full*10,
# #     cmap="terrain",
# #     vmin=-3.0,
# #     vmax=6.0,
# # )
# # # cbar = fig1.colorbar(mat)
# # # cbar.set_label('m MHW', rotation=270, labelpad=15)
# # ax1.set_title("Configuration 1")
# # ax1.set_ylabel("barrier width (m)")
# # ax1.set_xlabel("barrier length (m)")
# # plt.gca().xaxis.tick_bottom()
# # # ax1.add_patch(plt.Rectangle((20,175), 50, -40, lw=2, ec="k", fc="none"))
# # # ax1.add_patch(plt.Rectangle((70,175), 50, -40, lw=2, ec="m", fc="none"))
# # # ls = linestyle, lw = line width, ec = edge color, fc = fill color
# #
# # xtick_max = np.shape(full)[1]  # n_cols = x
# # x_ticks = np.array(range(0, xtick_max, 10))
# # x_tick_labels = x_ticks * 10
# # ytick_max = np.shape(full)[0]  # n_rows = y
# # y_ticks = np.array(range(0, ytick_max, 5))
# # y_tick_labels = y_ticks * 10
# # plt.xticks(x_ticks, x_tick_labels)
# # plt.yticks(y_ticks, y_tick_labels)
# # # fig2.subplots_adjust(wspace=0, hspace=0)
#
#
# interior_b3d_input = np.flip(section1_pre_int)
# # np.save(r"C:\Users\Lexi\PycharmProjects\CASCADE\cascade\data\outwash_data/NCB-default-elevation-left-damMHW-test", interior_b3d_input)
#
# berm_el = 0.11
# dunes_b3d = np.flip(section1_pre_dunes) - berm_el
# dunes_input = np.append(dunes_b3d[0], dunes_b3d[1], 0)
# # np.save(r"C:\Users\Lexi\PycharmProjects\CASCADE\cascade\data\outwash_data/NCB-default-dunes-left-dam-test", dunes_input)
#
# # the beach should be loaded in with the ocean on the bottom
# # np.save(r"C:\Users\Lexi\PycharmProjects\CASCADE\cascade\data\outwash_data/NCB-default-beach-left-dam-test", section1_pre_beach)
#
#
# # re-formatting the right section --------------------------------------------------------------------------------------
# pre_right = section2_pre[34:76]
# post_right = section2_post[33:75]
#
# difference_right = post_right-pre_right
# total_erosion_right = np.sum(difference_right)*10
# back_erosion_right = np.sum(difference_right[4:30, 16:42])*10
#
# # shrink=0.4
# #
# # fig2 = plt.figure()
# # fig2.tight_layout()
# # fig2.suptitle('Configuration 2: Observed Erosion', weight="bold")
# # ax1 = fig2.add_subplot(131)
# # mat = ax1.matshow(
# #     pre_right*10,
# #     cmap="terrain",
# #     vmin=-3.0,
# #     vmax=6.0,
# # )
# # cbar = fig2.colorbar(mat, shrink=shrink)
# # cbar.set_label('m MHW', rotation=270, labelpad=5)
# # ax1.set_title("Pre-Storm Elevation")
# # ax1.set_ylabel("barrier width (m)")
# # ax1.set_xlabel("barrier length (m)")
# # plt.gca().xaxis.tick_bottom()
# # xtick_max = np.shape(pre_right)[1]  # n_cols = x
# # x_ticks = np.array(range(0, xtick_max, 10))
# # x_tick_labels = x_ticks * 10
# # ytick_max = np.shape(pre_right)[0]  # n_rows = y
# # y_ticks = np.array(range(0, ytick_max, 10))
# # y_tick_labels = y_ticks * 10
# # plt.xticks(x_ticks, x_tick_labels)
# # plt.yticks(y_ticks, y_tick_labels)
# #
# #
# # ax2 = fig2.add_subplot(132)
# # mat = ax2.matshow(
# #     post_right*10,
# #     cmap="terrain",
# #     vmin=-3.0,
# #     vmax=6.0,
# # )
# # cbar = fig2.colorbar(mat, shrink=shrink)
# # cbar.set_label('m MHW', rotation=270, labelpad=5)
# # ax2.set_title("Post-Storm Elevation")
# # ax2.set_ylabel("barrier width (m)")
# # ax2.set_xlabel("barrier length (m)")
# # plt.gca().xaxis.tick_bottom()
# # xtick_max = np.shape(post_right)[1]  # n_cols = x
# # x_ticks = np.array(range(0, xtick_max, 10))
# # x_tick_labels = x_ticks * 10
# # ytick_max = np.shape(post_right)[0]  # n_rows = y
# # y_ticks = np.array(range(0, ytick_max, 10))
# # y_tick_labels = y_ticks * 10
# # plt.xticks(x_ticks, x_tick_labels)
# # plt.yticks(y_ticks, y_tick_labels)
# # # ax2.tick_params(left=False)
# # # ax2.set(yticklabels=[])
# # # fig2.subplots_adjust(wspace=0.3, hspace=0)
# #
# # # fig3 = plt.figure()
# # ax3 = fig2.add_subplot(133)
# # mat = ax3.matshow(
# #     difference_right*10,
# #     cmap="seismic",
# #     vmin=-5,
# #     vmax=5,
# # )
# # # Configuration 1 \n
# # cbar = fig2.colorbar(mat, shrink=shrink)
# # cbar.set_label('meters', rotation=270, labelpad=5)
# # ax3.set_title("Elevation Change")
# # ax3.set_ylabel("barrier width (m)")
# # ax3.set_xlabel("barrier length (m)")
# # plt.gca().xaxis.tick_bottom()
# #
# # ax3.add_patch(plt.Rectangle((15,30), 25, -24, lw=3, ec="k", fc="none"))
# #
# # xtick_max = np.shape(difference_right)[1]  # n_cols = x
# # x_ticks = np.array(range(0, xtick_max, 10))
# # x_tick_labels = x_ticks * 10
# # ytick_max = np.shape(difference_right)[0]  # n_rows = y
# # y_ticks = np.array(range(0, ytick_max, 10))
# # y_tick_labels = y_ticks * 10
# # plt.xticks(x_ticks, x_tick_labels)
# # plt.yticks(y_ticks, y_tick_labels)
# # fig2.subplots_adjust(top=1.25, wspace=0.4, hspace=0)
# pre_right = section2_pre[34:78]
# section2_pre_int = pre_right[0:30]
# section2_pre_dunes = pre_right[30:32]
# section2_pre_dunes[0:2, 0:4] += 0.1
# # section2_pre_dunes[0:2, 0:2] = section2_pre_dunes[3, 0:2]
# # section2_pre_dunes[0, 2:5] = section2_pre_dunes[2, 2:5]
# section2_pre_dunes = section2_pre_dunes[0:2]
# section2_pre_dunes[0, 0:20] += 0.05
# section2_pre_dunes[1, 0:16] += 0.1
# # section2_pre_beach = pre[32:]
# # section2_pre_beach = np.load(r"C:\Users\Lexi\PycharmProjects\CASCADE\cascade\data\outwash_data\NCB-default-beach.npy")
# # section2_pre_beach = pre_right[32:]
# section2_pre_beach = pre_right[34:]
#
# full_right = np.append(section2_pre_int, section2_pre_dunes, 0)
# full_right = np.append(full_right, section2_pre_beach, 0)
#
# fig2 = plt.figure()
# ax1 = fig2.add_subplot(121)
# mat = ax1.matshow(
#     full_right*10,
#     cmap="terrain",
#     vmin=-3.0,
#     vmax=6.0,
# )
# cbar = fig2.colorbar(mat, shrink=0.7)
# cbar.set_label('m MHW', rotation=270, labelpad=15)
# #
# # plt.hlines(29.5, -0.5, 49.5, color="k", linestyles='solid', linewidth=3)
# # plt.hlines(31.5, -0.5, 49.5, color="k", linestyles='solid', linewidth=3)
# #
# # ax1.set_title("Configuration 1: Model Input")
# # # ax2.set_ylabel("barrier width (m)")
# # ax1.set_xlabel("barrier length (m)")
# # plt.gca().xaxis.tick_bottom()
# # fig2.set_figheight(15)
# # # ax2.add_patch(plt.Rectangle((20,175), 50, -40, lw=2, ec="k", fc="none"))
# # # ax2.add_patch(plt.Rectangle((70,175), 50, -40, lw=2, ec="m", fc="none"))
# #
# # xtick_max = np.shape(full_right)[1]  # n_cols = x
# # x_ticks = np.array(range(0, xtick_max, 10))
# # x_tick_labels = x_ticks * 10
# # ytick_max = np.shape(full_right)[0]  # n_rows = y
# # y_ticks = np.array(range(0, ytick_max, 5))
# # y_tick_labels = y_ticks * 10
# # plt.xticks(x_ticks, x_tick_labels)
# # plt.yticks(y_ticks, y_tick_labels)
# # # ax2.tick_params(left=False)
# #
# #
# #
# # ax2 = fig2.add_subplot(122)
# # mat = ax2.matshow(
# #     full_right*10,
# #     cmap="terrain",
# #     vmin=-3.0,
# #     vmax=6.0,
# # )
# # cbar = fig2.colorbar(mat, shrink=0.7)
# # cbar.set_label('m MHW', rotation=270, labelpad=15)
# #
# # ax2.set_title("Configuration 2: Model Input")
# # ax2.set_ylabel("barrier width (m)")
# # ax2.set_xlabel("barrier length (m)")
# # plt.gca().xaxis.tick_bottom()
# # fig2.set_figheight(15)
# # plt.hlines(29.5, -0.5, 49.5, color="k", linestyles='solid', linewidth=3)
# # plt.hlines(31.5, -0.5, 49.5, color="k", linestyles='solid', linewidth=3)
# # # ax2.add_patch(plt.Rectangle((20,175), 50, -40, lw=2, ec="k", fc="none"))
# # # ax2.add_patch(plt.Rectangle((70,175), 50, -40, lw=2, ec="m", fc="none"))
# #
# # xtick_max = np.shape(full_right)[1]  # n_cols = x
# # x_ticks = np.array(range(0, xtick_max, 10))
# # x_tick_labels = x_ticks * 10
# # ytick_max = np.shape(full_right)[0]  # n_rows = y
# # y_ticks = np.array(range(0, ytick_max, 5))
# # y_tick_labels = y_ticks * 10
# # plt.xticks(x_ticks, x_tick_labels)
# # plt.yticks(y_ticks, y_tick_labels)
# # # ax2.tick_params(left=False)
# # # ax2.set(yticklabels=[])
# #
# # plt.subplots_adjust(hspace=0.5)
#
#
#
# interior_b3d_input = np.flip(section2_pre_int)
# # np.save(r"C:\Users\Lexi\PycharmProjects\CASCADE\cascade\data\outwash_data/NCB-default-elevation-right-damMHW-test", interior_b3d_input)
#
# berm_el = 0.11
# dunes_b3d = np.flip(section2_pre_dunes) - berm_el
# dunes_input = np.append(dunes_b3d[0], dunes_b3d[1], 0)
# # np.save(r"C:\Users\Lexi\PycharmProjects\CASCADE\cascade\data\outwash_data/NCB-default-dunes-right-dam-test", dunes_input)
#
# # the beach should be loaded in with the ocean on the bottom
# # np.save(r"C:\Users\Lexi\PycharmProjects\CASCADE\cascade\data\outwash_data/NCB-default-beach-right-damMHW", section2_pre_beach)
