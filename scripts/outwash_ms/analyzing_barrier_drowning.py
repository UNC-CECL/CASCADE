# loading in the npz files and saving variables that we want to plot together

# Lexi Van Blunk
# 7/26/2023

import numpy as np
from matplotlib import pyplot as plt

# ---------------------------------- set model parameters that change per run ------------------------------------------
storm_interval = 20        # 20 or 10 years
rname = "r025"             # "r025" or "r035"
config = 4                 # 1, 2, 3, or 4

# location of the npz files
datadir_b3d = "C:/Users/Lexi/PycharmProjects/CASCADE/cascade/data/outwash_data/storms/slope0pt03/run_output/{0}/overwash_only/".format(rname)
datadir_100 = "C:/Users/Lexi/PycharmProjects/CASCADE/cascade/data/outwash_data/storms/slope0pt03/run_output/{0}/outwash100/".format(rname)
datadir_50 = "C:/Users/Lexi/PycharmProjects/CASCADE/cascade/data/outwash_data/storms/slope0pt03/run_output/{0}/outwash50/".format(rname)
datadir_0 = "C:/Users/Lexi/PycharmProjects/CASCADE/cascade/data/outwash_data/storms/slope0pt03/run_output/{0}/outwash0/".format(rname)

# initialize empty arrays
drowning_array_b3d = np.zeros(100)
drowning_array_100 = np.zeros(100)
drowning_array_50 = np.zeros(100)
drowning_array_0 = np.zeros(100)

drown_year_array_b3d = np.zeros(100)
drown_year_array_100 = np.zeros(100)
drown_year_array_50 = np.zeros(100)
drown_year_array_0 = np.zeros(100)

for storm_num in range(1, 101):
    # b3d variables
    filename_b3d = "config{0}_b3d_startyr1_interval{1}yrs_Slope0pt03_{2}.npz".format(config, storm_interval, storm_num)
    file_b3d = datadir_b3d + filename_b3d
    b3d = np.load(file_b3d, allow_pickle=True)
    b3d_obj = b3d["cascade"][0]
    drowning_array_b3d[storm_num-1] = b3d_obj.barrier3d[0].drown_break
    if b3d_obj.barrier3d[0].TMAX < 100:
        drown_year_array_b3d[storm_num-1] = b3d_obj.barrier3d[0].TMAX

    # 100% variables
    filename_100 = "config{0}_outwash100_startyr1_interval{1}yrs_Slope0pt03_{2}.npz".format(config, storm_interval, storm_num)
    file_100 = datadir_100 + filename_100
    outwash100 = np.load(file_100, allow_pickle=True)
    outwash100_obj = outwash100["cascade"][0]
    drowning_array_100[storm_num-1] = outwash100_obj.barrier3d[0].drown_break
    if outwash100_obj.barrier3d[0].TMAX < 100:
        drown_year_array_100[storm_num-1] = outwash100_obj.barrier3d[0].TMAX

    # 50% variables
    filename_50 = "config{0}_outwash50_startyr1_interval{1}yrs_Slope0pt03_{2}.npz".format(config, storm_interval, storm_num)
    file_50 = datadir_50 + filename_50
    outwash50 = np.load(file_50, allow_pickle=True)
    outwash50_obj = outwash50["cascade"][0]
    drowning_array_50[storm_num-1] = outwash50_obj.barrier3d[0].drown_break
    if outwash50_obj.barrier3d[0].TMAX < 100:
        drown_year_array_50[storm_num - 1] = outwash50_obj.barrier3d[0].TMAX

    # 0% variables
    filename_0 = "config{0}_outwash0_startyr1_interval{1}yrs_Slope0pt03_{2}.npz".format(config, storm_interval, storm_num)
    file_0 = datadir_0 + filename_0
    outwash0 = np.load(file_0, allow_pickle=True)
    outwash0_obj = outwash0["cascade"][0]
    drowning_array_0[storm_num-1] = outwash0_obj.barrier3d[0].drown_break
    if outwash0_obj.barrier3d[0].TMAX < 100:
        drown_year_array_0[storm_num - 1] = outwash0_obj.barrier3d[0].TMAX

# because we have 100 storms and drowns = 1, the sum of the array is the percent that drown
percent_drown_b3d = np.sum(drowning_array_b3d)
print("{0}% of barriers drown for the overwash only scenario".format(int(percent_drown_b3d)))

percent_drown_100 = np.sum(drowning_array_100)
print("{0}% of barriers drown for the 100% outwash to shoreface scenario".format(int(percent_drown_100)))

percent_drown_50 = np.sum(drowning_array_50)
print("{0}% of barriers drown for the 50% outwash to shoreface scenario".format(int(percent_drown_50)))

percent_drown_0 = np.sum(drowning_array_0)
print("{0}% of barriers drown for the 0% outwash to shoreface scenario".format(int(percent_drown_0)))

# histogram of years that the barriers drown
plt.rcParams.update({"font.size": 12})
bins = 100
fig1 = plt.figure()
ax1 = fig1.add_subplot(131)
ax1.hist(drown_year_array_100, bins=bins)
ax1.set_title("100% washout")
ax1.set_ylabel("number of barriers that drown")
ax1.set_xlabel("drown year")
plt.gca().xaxis.tick_bottom()
ax1.set_xlim(left=1, right=100)
ax1.set_ylim(bottom=0, top=40)

ax1 = fig1.add_subplot(132)
ax1.hist(drown_year_array_50, bins=bins)
ax1.set_title("50% washout")
ax1.set_ylabel("frequency")
ax1.set_xlabel("drown year")
plt.gca().xaxis.tick_bottom()
ax1.set_xlim(left=1, right=100)
ax1.set_ylim(bottom=0, top=40)
ax1.set_ylabel(None)

ax1 = fig1.add_subplot(133)
ax1.hist(drown_year_array_0, bins=bins)
ax1.set_title("0% washout")
ax1.set_ylabel("frequency")
ax1.set_xlabel("drown year")
plt.gca().xaxis.tick_bottom()
ax1.set_xlim(left=1, right=100)
ax1.set_ylim(bottom=0, top=40)
ax1.set_ylabel(None)

# --------------------------------------------- r = 0.35 ---------------------------------------------------------------
rname = "r035"

# location of the npz files
datadir_b3d = "C:/Users/Lexi/PycharmProjects/CASCADE/cascade/data/outwash_data/storms/slope0pt03/run_output/{0}/overwash_only/".format(rname)
datadir_100 = "C:/Users/Lexi/PycharmProjects/CASCADE/cascade/data/outwash_data/storms/slope0pt03/run_output/{0}/outwash100/".format(rname)
# datadir_50 = "C:/Users/Lexi/PycharmProjects/CASCADE/cascade/data/outwash_data/storms/slope0pt03/run_output/{0}/outwash50/".format(rname)
datadir_0 = "C:/Users/Lexi/PycharmProjects/CASCADE/cascade/data/outwash_data/storms/slope0pt03/run_output/{0}/outwash0/".format(rname)

# initialize empty arrays
drowning_array_b3d = np.zeros(100)
drowning_array_100 = np.zeros(100)
drowning_array_50 = np.zeros(100)
drowning_array_0 = np.zeros(100)

drown_year_array_b3d = np.zeros(100)
drown_year_array_100 = np.zeros(100)
drown_year_array_50 = np.zeros(100)
drown_year_array_0 = np.zeros(100)

for storm_num in range(1, 101):
    # b3d variables
    filename_b3d = "config{0}_b3d_startyr1_interval{1}yrs_Slope0pt03_{2}.npz".format(config, storm_interval, storm_num)
    file_b3d = datadir_b3d + filename_b3d
    b3d = np.load(file_b3d, allow_pickle=True)
    b3d_obj = b3d["cascade"][0]
    drowning_array_b3d[storm_num-1] = b3d_obj.barrier3d[0].drown_break
    if b3d_obj.barrier3d[0].TMAX < 100:
        drown_year_array_b3d[storm_num-1] = b3d_obj.barrier3d[0].TMAX

    # 100% variables
    filename_100 = "config{0}_outwash100_startyr1_interval{1}yrs_Slope0pt03_{2}.npz".format(config, storm_interval, storm_num)
    file_100 = datadir_100 + filename_100
    outwash100 = np.load(file_100, allow_pickle=True)
    outwash100_obj = outwash100["cascade"][0]
    drowning_array_100[storm_num-1] = outwash100_obj.barrier3d[0].drown_break
    if outwash100_obj.barrier3d[0].TMAX < 100:
        drown_year_array_100[storm_num-1] = outwash100_obj.barrier3d[0].TMAX

    # # 50% variables
    # filename_50 = "config{0}_outwash50_startyr1_interval{1}yrs_Slope0pt03_{2}.npz".format(config, storm_interval, storm_num)
    # file_50 = datadir_50 + filename_50
    # outwash50 = np.load(file_50, allow_pickle=True)
    # outwash50_obj = outwash50["cascade"][0]
    # drowning_array_50[storm_num-1] = outwash50_obj.barrier3d[0].drown_break
    # if outwash50_obj.barrier3d[0].TMAX < 100:
    #     drown_year_array_50[storm_num - 1] = outwash50_obj.barrier3d[0].TMAX

    # 0% variables
    filename_0 = "config{0}_outwash0_startyr1_interval{1}yrs_Slope0pt03_{2}.npz".format(config, storm_interval, storm_num)
    file_0 = datadir_0 + filename_0
    outwash0 = np.load(file_0, allow_pickle=True)
    outwash0_obj = outwash0["cascade"][0]
    drowning_array_0[storm_num-1] = outwash0_obj.barrier3d[0].drown_break
    if outwash0_obj.barrier3d[0].TMAX < 100:
        drown_year_array_0[storm_num - 1] = outwash0_obj.barrier3d[0].TMAX

# because we have 100 storms and drowns = 1, the sum of the array is the percent that drown
percent_drown_b3d = np.sum(drowning_array_b3d)
print("{0}% of barriers drown for the overwash only scenario".format(int(percent_drown_b3d)))

percent_drown_100 = np.sum(drowning_array_100)
print("{0}% of barriers drown for the 100% outwash to shoreface scenario".format(int(percent_drown_100)))

# percent_drown_50 = np.sum(drowning_array_50)
# print("{0}% of barriers drown for the 50% outwash to shoreface scenario".format(int(percent_drown_50)))

percent_drown_0 = np.sum(drowning_array_0)
print("{0}% of barriers drown for the 0% outwash to shoreface scenario".format(int(percent_drown_0)))

# histogram of years that the barriers drown
plt.rcParams.update({"font.size": 12})
bins = 100
fig2 = plt.figure()
ax1 = fig2.add_subplot(131)
ax1.hist(drown_year_array_100, bins=bins)
ax1.set_title("100% washout")
ax1.set_ylabel("number of barriers that drown")
ax1.set_xlabel("drown year")
plt.gca().xaxis.tick_bottom()
ax1.set_xlim(left=1, right=100)
ax1.set_ylim(bottom=0, top=40)

# ax1 = fig2.add_subplot(132)
# ax1.hist(drown_year_array_50, bins=bins)
# ax1.set_title("50% washout")
# ax1.set_ylabel("frequency")
# ax1.set_xlabel("drown year")
# plt.gca().xaxis.tick_bottom()
# ax1.set_xlim(left=1, right=100)
# ax1.set_ylim(bottom=0, top=40)
# ax1.set_ylabel(None)

ax1 = fig2.add_subplot(133)
ax1.hist(drown_year_array_0, bins=bins)
ax1.set_title("0% washout")
ax1.set_ylabel("frequency")
ax1.set_xlabel("drown year")
plt.gca().xaxis.tick_bottom()
ax1.set_xlim(left=1, right=100)
ax1.set_ylim(bottom=0, top=40)
ax1.set_ylabel(None)
