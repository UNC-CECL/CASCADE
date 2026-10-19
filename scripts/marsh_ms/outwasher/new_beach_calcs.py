# this is derived from the "plotting best model configurations and avg configuration 4" jupyter notebook

import numpy as np
import os
from matplotlib import pyplot as plt

os.chdir(r"C:\Users\agfig\OneDrive - University of North Carolina at Chapel Hill\post-thesis\tuning")

# load configurations
config1 = np.load(
    r"config1\cx_0pt45_ki_9pt5E-3.npz",
    allow_pickle=True,
)
config1_obj = config1["cascade"][0]
config2 = np.load(
    r"config2\cx_1pt05_ki_8pt5E-3.npz",
    allow_pickle=True,
)
config2_obj = config2["cascade"][0]

config3 = np.load(
    r"config3\cx_6pt5_ki_8pt25.npz",
    allow_pickle=True,
)
config3_obj = config3["cascade"][0]

config4 = np.load(
    r"config4\c_0pt0134_ki_8pt75E-3.npz",
    allow_pickle=True,
)
config4_obj = config4["cascade"][0]

# =================================================================
# ------------------------ configuration 1 ------------------------
# =================================================================
plt.rcParams.update({"font.size": 15})
xlabel = "alongshore distance (m)"
ylabel = "cross-shore distance (m)"
shrink = 0.15

prestorm1 = config1_obj.outwash[0]._initial_full_domain
poststorm1 = config1_obj.outwash[0]._full_domain

elev_change1 = poststorm1 - prestorm1

# plot pre-storm section -------------------------------------------------------------------------------------------------
fig1 = plt.figure(figsize=[18, 13])
fig1.suptitle("Configuration 1: Modeled Erosion", weight="bold")
ax1 = fig1.add_subplot(131)
mat1_pre = ax1.matshow(
    prestorm1 * 10,
    cmap="terrain",
    vmin=-3.0,
    vmax=6.0,
)
cbar = fig1.colorbar(mat1_pre, shrink=shrink)
cbar.set_label("m MHW", rotation=270, labelpad=10)
ax1.set_title("Pre-Storm Elevation")
ax1.set_ylabel(ylabel)
ax1.set_xlabel(xlabel)
plt.gca().xaxis.tick_bottom()

xtick_max = np.shape(prestorm1)[1]  # n_cols = x
x_ticks = np.array(range(0, xtick_max, 10))
x_tick_labels = x_ticks * 10
plt.xticks(x_ticks, x_tick_labels)
ytick_max = np.shape(prestorm1)[0]  # n_rows = y
y_ticks = np.array(range(0, ytick_max, 10))
y_tick_labels = y_ticks * 10
plt.yticks(y_ticks, y_tick_labels)


# plot post-storm section -------------------------------------------------------------------------------------------------
ax1 = fig1.add_subplot(132)
mat1_post = ax1.matshow(
    poststorm1 * 10,
    cmap="terrain",
    vmin=-3.0,
    vmax=6.0,
)
cbar = fig1.colorbar(mat1_post, shrink=shrink)
cbar.set_label("m MHW", rotation=270, labelpad=10)
ax1.set_title("Post-Storm Elevation")
ax1.set_ylabel(ylabel)
ax1.set_xlabel(xlabel)
plt.gca().xaxis.tick_bottom()

xtick_max = np.shape(prestorm1)[1]  # n_cols = x
x_ticks = np.array(range(0, xtick_max, 10))
x_tick_labels = x_ticks * 10
plt.xticks(x_ticks, x_tick_labels)
ytick_max = np.shape(prestorm1)[0]  # n_rows = y
y_ticks = np.array(range(0, ytick_max, 10))
y_tick_labels = y_ticks * 10
plt.yticks(y_ticks, y_tick_labels)

# plot elevation change --------------------------------------------------------------------------------------------------
ax1 = fig1.add_subplot(133)
mat1_dif = ax1.matshow(
    elev_change1 * 10,
    cmap="seismic_r",
    vmin=-5,
    vmax=5,
)
cbar = fig1.colorbar(mat1_dif, shrink=shrink)
cbar.set_label("m MHW", rotation=270, labelpad=10)
ax1.set_title("Elevation Change")
ax1.set_ylabel(ylabel)
ax1.set_xlabel(xlabel)
plt.gca().xaxis.tick_bottom()

xtick_max = np.shape(prestorm1)[1]  # n_cols = x
x_ticks = np.array(range(0, xtick_max, 10))
x_tick_labels = x_ticks * 10
plt.xticks(x_ticks, x_tick_labels)
ytick_max = np.shape(prestorm1)[0]  # n_rows = y
y_ticks = np.array(range(0, ytick_max, 10))
y_tick_labels = y_ticks * 10
plt.yticks(y_ticks, y_tick_labels)

fig1.subplots_adjust(wspace=0.35, top=1.5)

plt.hlines(31.5, -0.5, 49.5, color="k", linestyles="solid", linewidth=2)

channel_erosion1_m = np.sum(elev_change1[0:32] * 10)  # m
channel_erosion1_m3 = channel_erosion1_m * 10 * 10
print(round(channel_erosion1_m3, -2))

print(
    "outwash volume is: ", int(round(config1_obj.outwash[0]._outwash_TS[1], -2)), " m3"
)

# beach erosion calcs
beach_erosion1_m = np.sum(elev_change1[32:] * 10)  # m
beach_erosion1_m3 = beach_erosion1_m * 10 * 10  # m3

print(
    "beach volume is: ", int(round(beach_erosion1_m3, -2)), " m3"
)

# =================================================================
# ------------------------ configuration 2 ------------------------
# =================================================================
prestorm2 = config2_obj.outwash[0]._initial_full_domain
poststorm2 = config2_obj.outwash[0]._full_domain

elev_change2 = poststorm2 - prestorm2

# plot pre-storm section -------------------------------------------------------------------------------------------------
fig2 = plt.figure(figsize=[18, 13])
fig2.suptitle("Configuration 2: Modeled Erosion", weight="bold")
ax2 = fig2.add_subplot(131)
mat2_pre = ax2.matshow(
    prestorm2 * 10,
    cmap="terrain",
    vmin=-3.0,
    vmax=6.0,
)
cbar = fig2.colorbar(mat2_pre, shrink=shrink)
cbar.set_label("m MHW", rotation=270, labelpad=10)
ax2.set_title("Pre-Storm Elevation")
ax2.set_ylabel(ylabel)
ax2.set_xlabel(xlabel)
plt.gca().xaxis.tick_bottom()

xtick_max = np.shape(prestorm2)[1]  # n_cols = x
x_ticks = np.array(range(0, xtick_max, 10))
x_tick_labels = x_ticks * 10
plt.xticks(x_ticks, x_tick_labels)
ytick_max = np.shape(prestorm2)[0]  # n_rows = y
y_ticks = np.array(range(0, ytick_max, 10))
y_tick_labels = y_ticks * 10
plt.yticks(y_ticks, y_tick_labels)


# plot post-storm section -------------------------------------------------------------------------------------------------
ax2 = fig2.add_subplot(132)
mat2_pre = ax2.matshow(
    poststorm2 * 10,
    cmap="terrain",
    vmin=-3.0,
    vmax=6.0,
)
cbar = fig2.colorbar(mat2_pre, shrink=shrink)
cbar.set_label("m MHW", rotation=270, labelpad=10)
ax2.set_title("Post-Storm Elevation")
ax2.set_ylabel(ylabel)
ax2.set_xlabel(xlabel)
plt.gca().xaxis.tick_bottom()

xtick_max = np.shape(prestorm2)[1]  # n_cols = x
x_ticks = np.array(range(0, xtick_max, 10))
x_tick_labels = x_ticks * 10
plt.xticks(x_ticks, x_tick_labels)
ytick_max = np.shape(prestorm2)[0]  # n_rows = y
y_ticks = np.array(range(0, ytick_max, 10))
y_tick_labels = y_ticks * 10
plt.yticks(y_ticks, y_tick_labels)

# plot elevation change --------------------------------------------------------------------------------------------------
ax2 = fig2.add_subplot(133)
mat2_dif = ax2.matshow(
    elev_change2 * 10,
    cmap="seismic_r",
    vmin=-5,
    vmax=5,
)
cbar = fig2.colorbar(mat2_dif, shrink=shrink)
cbar.set_label("m MHW", rotation=270, labelpad=10)
ax2.set_title("Elevation Change")
ax2.set_ylabel(ylabel)
ax2.set_xlabel(xlabel)
plt.gca().xaxis.tick_bottom()

xtick_max = np.shape(prestorm2)[1]  # n_cols = x
x_ticks = np.array(range(0, xtick_max, 10))
x_tick_labels = x_ticks * 10
plt.xticks(x_ticks, x_tick_labels)
ytick_max = np.shape(prestorm2)[0]  # n_rows = y
y_ticks = np.array(range(0, ytick_max, 10))
y_tick_labels = y_ticks * 10
plt.yticks(y_ticks, y_tick_labels)

fig2.subplots_adjust(wspace=0.35, top=1.5)

plt.hlines(31.5, -0.5, 49.5, color="k", linestyles="solid", linewidth=2)

channel_erosion2_m = np.sum(elev_change2[0:32] * 10)  # m
channel_erosion2_m3 = channel_erosion2_m * 10 * 10
print(round(channel_erosion2_m3, -2))

print(
    "outwash volume is: ", int(round(config2_obj.outwash[0]._outwash_TS[1], -2)), " m3"
)

# beach erosion calcs
beach_erosion2_m = np.sum(elev_change2[32:] * 10)  # m
beach_erosion2_m3 = beach_erosion2_m * 10 * 10  # m3

print(
    "beach volume is: ", int(round(beach_erosion2_m3, -2)), " m3"
)

# =================================================================
# ------------------------ configuration 3 ------------------------
# =================================================================
prestorm3 = config3_obj.outwash[0]._initial_full_domain
poststorm3 = config3_obj.outwash[0]._full_domain

elev_change3 = poststorm3 - prestorm3

# plot pre-storm section -------------------------------------------------------------------------------------------------
fig3 = plt.figure(figsize=[18, 13])
fig3.suptitle("Configuration 3: Modeled Erosion", weight="bold")
ax3 = fig3.add_subplot(131)
mat3_pre = ax3.matshow(
    prestorm3 * 10,
    cmap="terrain",
    vmin=-3.0,
    vmax=6.0,
)
cbar = fig3.colorbar(mat3_pre, shrink=shrink)
cbar.set_label("m MHW", rotation=270, labelpad=10)
ax3.set_title("Pre-Storm Elevation")
ax3.set_ylabel(ylabel)
ax3.set_xlabel(xlabel)
plt.gca().xaxis.tick_bottom()

xtick_max = np.shape(prestorm3)[1]  # n_cols = x
x_ticks = np.array(range(0, xtick_max, 10))
x_tick_labels = x_ticks * 10
plt.xticks(x_ticks, x_tick_labels)
ytick_max = np.shape(prestorm3)[0]  # n_rows = y
y_ticks = np.array(range(0, ytick_max, 10))
y_tick_labels = y_ticks * 10
plt.yticks(y_ticks, y_tick_labels)

# plt.hlines(24.5, -0.5, 49.5, color="k", linestyles='solid', linewidth=1)

# plot post-storm section -------------------------------------------------------------------------------------------------
ax3 = fig3.add_subplot(132)
mat3_pre = ax3.matshow(
    poststorm3 * 10,
    cmap="terrain",
    vmin=-3.0,
    vmax=6.0,
)
cbar = fig3.colorbar(mat3_pre, shrink=shrink)
cbar.set_label("m MHW", rotation=270, labelpad=10)
ax3.set_title("Post-Storm Elevation")
ax3.set_ylabel(ylabel)
ax3.set_xlabel(xlabel)
plt.gca().xaxis.tick_bottom()

xtick_max = np.shape(prestorm3)[1]  # n_cols = x
x_ticks = np.array(range(0, xtick_max, 10))
x_tick_labels = x_ticks * 10
plt.xticks(x_ticks, x_tick_labels)
ytick_max = np.shape(prestorm3)[0]  # n_rows = y
y_ticks = np.array(range(0, ytick_max, 10))
y_tick_labels = y_ticks * 10
plt.yticks(y_ticks, y_tick_labels)
#
# plt.hlines(24.5, -0.5, 49.5, color="k", linestyles='solid', linewidth=1)

# plot elevation change --------------------------------------------------------------------------------------------------
ax3 = fig3.add_subplot(133)
mat3_dif = ax3.matshow(
    elev_change3 * 10,
    cmap="seismic_r",
    vmin=-5,
    vmax=5,
)
cbar = fig3.colorbar(mat3_dif, shrink=shrink)
cbar.set_label("m MHW", rotation=270, labelpad=10)
ax3.set_title("Elevation Change")
ax3.set_ylabel(ylabel)
ax3.set_xlabel(xlabel)
plt.gca().xaxis.tick_bottom()

xtick_max = np.shape(prestorm3)[1]  # n_cols = x
x_ticks = np.array(range(0, xtick_max, 10))
x_tick_labels = x_ticks * 10
plt.xticks(x_ticks, x_tick_labels)
ytick_max = np.shape(prestorm3)[0]  # n_rows = y
y_ticks = np.array(range(0, ytick_max, 10))
y_tick_labels = y_ticks * 10
plt.yticks(y_ticks, y_tick_labels)

fig3.subplots_adjust(wspace=0.35, top=1.5)

plt.hlines(24.5, -0.5, 49.5, color="k", linestyles="solid", linewidth=2)

back_barrier_dune_erosion3_m = np.sum(elev_change3[0:25] * 10)  # m
back_barrier_dune_erosion3_m3 = back_barrier_dune_erosion3_m * 10 * 10
print(round(back_barrier_dune_erosion3_m3, -2))

print(
    "outwash volume is: ", int(round(config3_obj.outwash[0]._outwash_TS[1], -2)), " m3"
)

# beach erosion calcs
beach_erosion3_m = np.sum(elev_change3[25:] * 10)  # m
beach_erosion3_m3 = beach_erosion3_m * 10 * 10  # m3

print(
    "beach volume is: ", int(round(beach_erosion3_m3, -2)), " m3"
)

# =================================================================
# ------------------------ configuration 4 ------------------------
# =================================================================
prestorm4 = config4_obj.outwash[0]._initial_full_domain
poststorm4 = config4_obj.outwash[0]._full_domain

elev_change4 = poststorm4 - prestorm4

# plot pre-storm section -------------------------------------------------------------------------------------------------
fig4 = plt.figure(figsize=[18, 13])
fig4.suptitle("Configuration 4: Modeled Erosion", weight="bold")
ax4 = fig4.add_subplot(131)
mat4_pre = ax4.matshow(
    prestorm4 * 10,
    cmap="terrain",
    vmin=-3.0,
    vmax=6.0,
)
cbar = fig4.colorbar(mat4_pre, shrink=shrink)
cbar.set_label("m MHW", rotation=270, labelpad=10)
ax4.set_title("Pre-Storm Elevation")
ax4.set_ylabel(ylabel)
ax4.set_xlabel(xlabel)
plt.gca().xaxis.tick_bottom()

xtick_max = np.shape(prestorm4)[1]  # n_cols = x
x_ticks = np.array(range(0, xtick_max, 10))
x_tick_labels = x_ticks * 10
plt.xticks(x_ticks, x_tick_labels)
ytick_max = np.shape(prestorm4)[0]  # n_rows = y
y_ticks = np.array(range(0, ytick_max, 10))
y_tick_labels = y_ticks * 10
plt.yticks(y_ticks, y_tick_labels)

# plot post-storm section -------------------------------------------------------------------------------------------------
ax4 = fig4.add_subplot(132)
mat4_post = ax4.matshow(
    poststorm4 * 10,
    cmap="terrain",
    vmin=-3.0,
    vmax=6.0,
)
cbar = fig4.colorbar(mat4_post, shrink=shrink)
cbar.set_label("m MHW", rotation=270, labelpad=10)
ax4.set_title("Post-Storm Elevation")
ax4.set_ylabel(ylabel)
ax4.set_xlabel(xlabel)
plt.gca().xaxis.tick_bottom()

xtick_max = np.shape(prestorm4)[1]  # n_cols = x
x_ticks = np.array(range(0, xtick_max, 10))
x_tick_labels = x_ticks * 10
plt.xticks(x_ticks, x_tick_labels)
ytick_max = np.shape(prestorm4)[0]  # n_rows = y
y_ticks = np.array(range(0, ytick_max, 10))
y_tick_labels = y_ticks * 10
plt.yticks(y_ticks, y_tick_labels)
#
# plt.hlines(24.5, -0.5, 49.5, color="k", linestyles='solid', linewidth=1)

# plot elevation change --------------------------------------------------------------------------------------------------
ax4 = fig4.add_subplot(133)
mat4_dif = ax4.matshow(
    elev_change4 * 10,
    cmap="seismic_r",
    vmin=-5,
    vmax=5,
)
cbar = fig4.colorbar(mat4_dif, shrink=shrink)
cbar.set_label("m MHW", rotation=270, labelpad=10)
ax4.set_title("Elevation Change")
ax4.set_ylabel(ylabel)
ax4.set_xlabel(xlabel)
plt.gca().xaxis.tick_bottom()

xtick_max = np.shape(prestorm4)[1]  # n_cols = x
x_ticks = np.array(range(0, xtick_max, 10))
x_tick_labels = x_ticks * 10
plt.xticks(x_ticks, x_tick_labels)
ytick_max = np.shape(prestorm4)[0]  # n_rows = y
y_ticks = np.array(range(0, ytick_max, 10))
y_tick_labels = y_ticks * 10
plt.yticks(y_ticks, y_tick_labels)

fig4.subplots_adjust(wspace=0.35, top=1.5)

plt.hlines(27.5, -0.5, 49.5, color="k", linestyles="solid", linewidth=2)

back_barrier_dune_erosion4_m = np.sum(elev_change4[0:28] * 10)  # m
back_barrier_dune_erosion4_m3 = back_barrier_dune_erosion4_m * 10 * 10
print(round(back_barrier_dune_erosion4_m3, -2))

print(
    "outwash volume is: ", int(round(config4_obj.outwash[0]._outwash_TS[1], -2)), " m3"
)

# beach erosion calcs
beach_erosion4_m = np.sum(elev_change4[28:] * 10)  # m
beach_erosion4_m3 = beach_erosion4_m * 10 * 10  # m3

print(
    "beach volume is: ", int(round(beach_erosion4_m3, -2)), " m3"
)