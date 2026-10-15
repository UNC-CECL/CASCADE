# Lexi (Van Blunk) Fiegelist
# 10/14/2026

# this code makes figures from the outputs from folders 1 and 2 in the same directory

import os
import pandas as pd
from matplotlib import pyplot as plt

plt.rcParams["font.size"]=14
start_year = 1984
end_year = 2025
datadir = r"C:\Users\agfig\OneDrive - University of North Carolina at Chapel Hill\UNC\data\CoastSat\step2_output"
os.chdir(datadir)
# set the domains used in cascade
d_start = 3
d_end = 25

# endpoint per zone
endpoint_file = "endpoint_per_zone_{0}_{1}.csv".format(start_year, end_year)
endpoints = pd.read_csv(endpoint_file)
endpoints = endpoints[(endpoints["zone_id"]>=d_start)&(endpoints["zone_id"]<=d_end)]
# rates per zone
lrr_rates_file = "rates_per_zone_{0}_{1}.csv".format(start_year, end_year)
lrr_rates = pd.read_csv(lrr_rates_file)
lrr_rates = lrr_rates[(lrr_rates["zone_id"]>=d_start)&(lrr_rates["zone_id"]<=d_end)]

# plots
domains = lrr_rates["zone_id"].values
y_endpoints = endpoints["mean_change_m"].values
y_lrr = lrr_rates["mean_rate_m_yr"].values
# set ylims
y_end_lim = 200
y_lrr_lim = 5.5
# endpoints
fig1 = plt.figure(figsize=(10, 6))
ax = fig1.add_subplot(111)
ax.plot(domains, y_endpoints, marker="o", markersize=4)
plt.axhline(0, color="grey", linewidth=0.5)
ax.set(
    xlabel="domain",
    ylabel="change in shoreline position (m)",
    title="{0} to {1}".format(start_year, end_year),
    xticks=domains,
    ylim=[-y_end_lim, y_end_lim],
)
fig1.tight_layout()
# lrrs
fig2 = plt.figure(figsize=(10, 6))
ax = fig2.add_subplot(111)
ax.plot(domains, y_lrr, marker="o", markersize=4)
plt.axhline(0, color="grey", linewidth=0.5)
ax.set(
    xlabel="domain",
    ylabel="LRR rates (m/yr)",
    title="{0} to {1}".format(start_year, end_year),
    xticks=domains,
    ylim=[-y_lrr_lim, y_lrr_lim],
)
fig2.tight_layout()