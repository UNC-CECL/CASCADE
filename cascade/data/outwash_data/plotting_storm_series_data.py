import numpy as np
from matplotlib import pyplot as plt
import pandas as pd

datadir = r"C:\Users\Lexi\PycharmProjects\CASCADE\cascade\data\outwash_data/"

cascade_storms_file = "cascade-default-storms.npy"
slope_03_file = "StormSeries_100yrs_NCB_Berm1pt46m_Slope0pt03_01.npy"
slope_002_file = "StormSeries_100yrs_NCB_Berm1pt46m_Slope0pt002_01.npy"


# reminder that the arrays are storm year, rhigh (dam), rlow (dam), alpha, duration (hrs)
cascade_storm_series = np.load(datadir+cascade_storms_file)[0:864, :]
cascade_years = cascade_storm_series[:, 0]
cascade_rhigh = cascade_storm_series[:, 1] * 10
cascade_rlow = cascade_storm_series[:, 2] * 10
cascade_dur = cascade_storm_series[:, 4]

slope_03_storm_series = np.load(datadir+slope_03_file)
slope_03_years = slope_03_storm_series[:, 0]
slope_03_rhigh = slope_03_storm_series[:, 1] * 10
slope_03_rlow = slope_03_storm_series[:, 2] * 10
slope_03_dur = slope_03_storm_series[:, 4]

slope_002_storm_series = np.load(datadir+slope_002_file)
slope_002_years = slope_002_storm_series[:, 0]
slope_002_rhigh = slope_002_storm_series[:, 1] * 10
slope_002_rlow = slope_002_storm_series[:, 2] * 10
slope_002_dur = slope_002_storm_series[:, 4]

cascade_unique_years = np.unique(cascade_years)
cascade_storms_per_year = []
for year in cascade_unique_years:
    cascade_storms_per_year = np.append(cascade_storms_per_year, sum(cascade_years == year))

slope_03_unique_years = np.unique(slope_03_years)
slope_03_storms_per_year = []
for year in slope_03_unique_years:
    slope_03_storms_per_year = np.append(slope_03_storms_per_year, sum(slope_03_years == year))

slope_002_unique_years = np.unique(slope_002_years)
slope_002_storms_per_year = []
for year in slope_002_unique_years:
    slope_002_storms_per_year = np.append(slope_002_storms_per_year, sum(slope_002_years == year))

# make an array that stores average data values
# [n storms per year, rhigh, rlow, duration]
# row 0 = cascade
# row 1 = slope 0.03
# row 2 = slope 0.002

cols = ["n storms per year", "rhigh", "rlow", "duration"]
rows = ["cascade", "slope 0.03", "slope 0.002"]

avg_n_storms_cascade = np.mean(cascade_storms_per_year)
avg_n_storms_slope_03 = np.mean(slope_03_storms_per_year)
avg_n_storms_slope_002 = np.mean(slope_002_storms_per_year)

data_array = np.zeros([3, 4])
# n storms per year
data_array[0, 0] = avg_n_storms_cascade
data_array[1, 0] = avg_n_storms_slope_03
data_array[2, 0] = avg_n_storms_slope_002
# rhigh
data_array[0, 1] = np.mean(cascade_rhigh)
data_array[1, 1] = np.mean(slope_03_rhigh)
data_array[2, 1] = np.mean(slope_002_rhigh)
# rlow
data_array[0, 2] = np.mean(cascade_rlow)
data_array[1, 2] = np.mean(slope_03_rlow)
data_array[2, 2] = np.mean(slope_002_rlow)
# duration
data_array[0, 3] = np.mean(cascade_dur)
data_array[1, 3] = np.mean(slope_03_dur)
data_array[2, 3] = np.mean(slope_002_dur)

df = pd.DataFrame(data=data_array, index=rows, columns=cols)

# plotting the data
plotters = False
if plotters:
    # figure 1: histograms of n storms per year
    fig1 = plt.figure()
    ax1 = fig1.add_subplot(311)
    ax1.bar(cascade_unique_years, cascade_storms_per_year, label="cascade")
    plt.xticks(np.arange(0,100,5))
    ax1.set_ylim(0,42)
    ax1.set_ylabel("number of storms", weight='bold')
    # ax1.set_xlabel("storm year", weight='bold')
    plt.legend()

    ax2 = fig1.add_subplot(312)
    ax2.bar(slope_03_unique_years, slope_03_storms_per_year, label="slope 0.03", color="orange")
    plt.xticks(np.arange(0,100,5))
    ax2.set_ylim(0,42)
    ax2.set_ylabel("number of storms", weight='bold')
    # ax2.set_xlabel("storm year", weight='bold')
    plt.legend()

    ax3 = fig1.add_subplot(313)
    ax3.bar(slope_002_unique_years, slope_002_storms_per_year, label="slope 0.002", color="green")
    plt.xticks(np.arange(0,100,5))
    ax3.set_ylim(0,42)
    ax3.set_ylabel("number of storms", weight='bold')
    ax3.set_xlabel("storm year", weight='bold')
    plt.legend()

    # figure 2: n storms per year line plots
    plt.figure(2)
    plt.plot(cascade_unique_years, cascade_storms_per_year, label="cascade")
    plt.plot(slope_03_unique_years, slope_03_storms_per_year, label="slope 0.03", color="orange")
    plt.plot(slope_002_unique_years, slope_002_storms_per_year, label="slope 0.002", color="green")
    plt.xticks(np.arange(0,100,5))
    plt.legend()
    plt.scatter(cascade_unique_years, cascade_storms_per_year, label="cascade")
    plt.scatter(slope_03_unique_years, slope_03_storms_per_year, label="slope 0.03", color="orange")
    plt.scatter(slope_002_unique_years, slope_002_storms_per_year, label="slope 0.002", color="green")
    plt.xticks(np.arange(0,100,5))
    plt.ylabel("number of storms", weight='bold')
    plt.xlabel("storm year", weight='bold')

    # figure 3: R high plots
    plt.figure(3)
    plt.scatter(cascade_years, cascade_rhigh, label="cascade")
    plt.scatter(slope_03_years, slope_03_rhigh, label="slope 0.03", color="orange")
    plt.scatter(slope_002_years, slope_002_rhigh, label="slope 0.002", color="green")
    plt.ylabel("R-high (m)", weight='bold')
    plt.xlabel("storm year", weight='bold')
    plt.xticks(np.arange(0,100,5))
    plt.legend()

    # figure 4: R low plots
    plt.figure(4)
    plt.scatter(cascade_years, cascade_rlow, label="cascade")
    plt.scatter(slope_03_years, slope_03_rlow, label="slope 0.03", color="orange")
    plt.scatter(slope_002_years, slope_002_rlow, label="slope 0.002", color="green")
    plt.ylabel("R-low (m)", weight='bold')
    plt.xlabel("storm year", weight='bold')
    plt.xticks(np.arange(0,100,5))
    plt.legend()

    # figure 5: duration plots
    plt.figure(5)
    plt.scatter(cascade_years, cascade_dur, label="cascade")
    plt.scatter(slope_03_years, slope_03_dur, label="slope 0.03", color="orange")
    plt.scatter(slope_002_years, slope_002_dur, label="slope 0.002", color="green")
    plt.ylabel("Duration (hrs)", weight='bold')
    plt.xlabel("storm year", weight='bold')
    plt.xticks(np.arange(0,100,5))
    plt.legend()