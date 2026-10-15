# Lexi (Van Blunk) Fiegelist
# 10/14/2026

# code for comparing the offsets created by using the dune line (see CASCADE Input Prep
# document for more details) and the method of just using the coastsat curvature for the
# start year (here, that is 2004)
import os
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from scipy import stats
import math

plt.rcParams["font.size"]=14
save_dir = r"C:\Users\agfig\OneDrive - University of North Carolina at Chapel Hill\UNC\data\CoastSat\offsets"
if not os.path.exists(save_dir):
    os.makedirs(save_dir)

# set the year
year = 2004
# set the domains used in cascade
d_start = 3
d_end = 25

# --------------------------------- dune offsets (CSV is in meters) ---------------------------------
# processing script already includes only the domains I wanted to include in cascade (ordered 3-25)
dune_offsets_file = r"C:\Users\agfig\model\calibration\offsets\Island_Dune_Offsets_{0}_CASCADE_Input.csv".format(year)
dune_offsets = np.loadtxt(dune_offsets_file, skiprows=1, delimiter=',')
# the dune offsets also assume that a more positive value is more landward (positive movement is landward), so the largest value in this array is the most landward
dune_offsets = (dune_offsets - np.max(dune_offsets))*-1

# --------------------------------- coastsat offsets (also in meters) ---------------------------------
# Part 1 - load data that was already processed using shoreline_endpoint_template.py
# these are the average position from a baseline first averaged over the start year (all positions recorded in 2004,
# giving 1 value per transect), then averaged over the transects that fell into a common domain (giving 1 value per domain)
coastsat_file = r"C:\Users\agfig\OneDrive - University of North Carolina at Chapel Hill\UNC\data\CoastSat\step2_output\endpoint_per_zone_2004_2014.csv"
coastsat_df = pd.read_csv(coastsat_file)
# only include domains that I want
coastsat_df = coastsat_df[(coastsat_df["zone_id"]>=d_start)&(coastsat_df["zone_id"]<=d_end)]
# mean start position column is for 2004, note they are also ordered with domain 3 first
coastsat_offsets = coastsat_df["mean_position_start_m"].values  # distance from the baseline to the 2004 shoreline

# Part 2 - shift these based on the offset of each transect baseline from one another
# I exported the XY points of each coastsat baseline into state plane coords (meters)
xy_points = pd.read_excel(os.path.join(save_dir, "baseline_points_stateplane.xls"))
xy_points = xy_points.sort_values(by="id")
xy_points = xy_points[(xy_points["zone_id"]>=d_start)&(xy_points["zone_id"]<=d_end)]
# group by zone and average the values
xy_per_domain = xy_points.groupby("zone_id").agg(
    domain=("zone_id", "mean"),
    mean_x=("POINT_X", "mean"),
    mean_y=("POINT_Y", "mean"))
xy_per_domain["domain"] = xy_per_domain["domain"].astype(int)
x = xy_per_domain["mean_x"].values
y = xy_per_domain["mean_y"].values
# only use a subset of values for the fit bc these are the most shoreline-parallel
xfit = xy_per_domain["mean_x"].values[2:8]
yfit = xy_per_domain["mean_y"].values[2:8]
fit = stats.linregress(xfit, yfit)
slope = fit.slope
intercept = fit.intercept
r2 = fit.rvalue ** 2
fit_y = slope*x + intercept
# plot to see what everything looks like
fig2 = plt.figure(figsize = (8,10))
plt.scatter(x, y, label="baseline points per domain")
plt.plot(x, fit_y, label="shoreline parallel fit")
plt.xlabel("state plane (m)")
plt.ylabel("state plane (m)")
plt.legend()
fig2.tight_layout()
fig2.savefig(os.path.join(save_dir, "coastsat_baseline_point_offsets.png"))
# find perpendicular distance of each point to the shoreline-parallel line to get approximate offsets from each other
a = slope
b = -1
c = intercept
distances = (a*x + b*y + c) / math.sqrt(a*a + b*b)  # meters

# Part 3 - combine the offset from each baseline with their offsets from each other
rel_coastsat_offsets = coastsat_offsets + distances
# reference to the minimum offset position
final_coastsat_offsets = rel_coastsat_offsets - np.min(rel_coastsat_offsets)
# convert to cascade orientation where positive is seaward and save
coastsat_offsets_cascade = (final_coastsat_offsets*-1) - np.min((final_coastsat_offsets*-1))
rows = len(coastsat_offsets_cascade)
cols = 2
# to make it clear for anyone opening the file, save the domains with the values
coastsat_offsets_cascade_save = np.zeros([rows,cols])
coastsat_offsets_cascade_save[:,0] = np.arange(d_start, d_end+1)
coastsat_offsets_cascade_save[:,1] = coastsat_offsets_cascade
np.save(
    os.path.join(save_dir, "coastsat_offsets_casc_{0}.npy".format(year)),
    coastsat_offsets_cascade_save
)
np.savetxt(
    os.path.join(save_dir, "coastsat_offsets_casc_{0}.csv".format(year)),
    coastsat_offsets_cascade_save,
    delimiter = ",",
    fmt='%.2f'
)


# make a plot
labels = np.arange(d_start, d_end+1)
y_vals = np.arange(0,len(dune_offsets))
fig1 = plt.figure(figsize = (5,10))
plt.plot(dune_offsets, y_vals, label="dune offsets")
plt.plot(final_coastsat_offsets, y_vals, label="coastsat offsets")
plt.xlabel(r"$\leftarrow$ landward || seaward $\rightarrow$" +  "\n initial position (m)")
plt.ylabel("domain")
plt.yticks(y_vals, labels)
plt.legend()
fig1.tight_layout()

fig1.savefig(os.path.join(save_dir, "dune_vs_coastsat_offsets_{0}.png".format(year)))