import pandas as pd
import numpy as np
import os
import matplotlib.pyplot as plt


def calculate_slope(x1, x2, elev):
    y2 = elev[x2]
    y1 = elev[x1]
    slope = (y2-y1)/(x2-x1)
    return round(slope,3), y1, y2

# set directory to the text files folder
directory = r'C:\Users\agfig\OneDrive - University of North Carolina at Chapel Hill\UNC\data\DEM\berm_slope_analysis'
file_2004 = os.path.join(directory, "2004_transects.xls")
file_2014 = os.path.join(directory, "2014_transects.xls")
file_2020 = os.path.join(directory, "2020_transects.xls")
df_2004 = pd.read_excel(file_2004)
df_2014 = pd.read_excel(file_2014)
df_2020 = pd.read_excel(file_2020)

domains = np.arange(26,2,-1)
transects = np.unique(df_2004["transectID"].values)
berm_elevs = [1.5,1.4,1.5,1.5,1.5,2.1,2,2,2.2,2.2,2.1,2.2,2.2,2.3,2.1,2.3,2.3,2.2,2,2.3,2,2,2.3,1.9]

plt.rcParams["font.size"] = 14
save_on = False
save_slopes_on = False
plot_grid = False
save_berms_on = False
d = 0
beach_range = 15  # points which are spaced every 1 meter
beach_slopes_2014 = []
beach_slopes_2020 = []
version = "final_berms"

savedir = os.path.join(directory, "results/odd_only")
if not os.path.exists(savedir):
    os.makedirs(savedir)

# plot and save transect figures
for t in transects:
    if t % 2 != 0:
        # get transect elevations
        elev_2004 = df_2004[df_2004["transectID"] == t].sort_values(by='pointID').elevation.values
        elev_2014 = df_2014[df_2014["transectID"] == t].sort_values(by='pointID').elevation.values
        elev_2020 = df_2020[df_2020["transectID"] == t].sort_values(by='pointID').elevation.values
        # assume beach elevation for streamlined beach slope calcs
        beach_elev = 1
        # find the first position of the value closest to 1 in the elev arrays
        dif_2014 = elev_2014 - beach_elev
        dif_2020 = elev_2020 - beach_elev
        # first value closest to 1
        first_beach_2014 = dif_2014[dif_2014>0][0]
        first_beach_2020 = dif_2020[dif_2020 > 0][0]
        # index position of that value
        beach_2014_pos = np.where(dif_2014==first_beach_2014)[0][0]
        beach_2020_pos = np.where(dif_2020 == first_beach_2020)[0][0]
        # use the surrounding 50 cells to calc beach slope
        # 2014
        # positive side
        if t == 7 or t == 9 or t == 47:
            pos1_2014 = beach_2014_pos + 10
        elif t == 23:
            pos1_2014 = beach_2014_pos + 20
        elif t == 31:
            pos1_2014 = beach_2014_pos + 5
        else:
            pos1_2014 = beach_2014_pos + beach_range
        # negative side
        if t == 15 or t==17 or t==19 or t==25 or t==47:
            pos2_2014 = beach_2014_pos - 5
        elif t == 27 or t==29:
            pos2_2014 = beach_2014_pos - 1
        # elif t == 29:
        #     pos2_2014 = beach_2014_pos - 25
        else:
            pos2_2014 = beach_2014_pos - beach_range
        # 2020
        # positive side
        if t == 5 or t == 13 or t == 21 or t==29 or t==35 or t==37 or t==39 or t==41:
            pos1_2020 = beach_2020_pos + 10
        elif t==27 or t==33 or t == 43:
            pos1_2020 = beach_2020_pos + 5
        else:
            pos1_2020 = beach_2020_pos + beach_range
        # negative side
        if t == 45 or t == 47:
            pos2_2020 = beach_2020_pos - 5
        else:
            pos2_2020 = beach_2020_pos - beach_range
        # calculate slope of specific year
        beach_slope_2014, y1_2014, y2_2014 = calculate_slope(x1=pos1_2014, x2=pos2_2014, elev=elev_2014)
        beach_slopes_2014.append(beach_slope_2014)
        beach_slope_2020, y1_2020, y2_2020 = calculate_slope(x1=pos1_2020, x2=pos2_2020, elev=elev_2020)
        beach_slopes_2020.append(beach_slope_2020)
        # create the plot
        fig = plt.subplots(1, 1, figsize=(12, 6))
        # plot elevations
        plt.plot(elev_2004, label="2004")
        plt.plot(elev_2014, ls="dotted", label="2014")
        plt.plot(elev_2020, ls="dotted", label="2020")
        # plot approximate beach elevation and values used to calculate the slope
        plt.hlines(y=beach_elev, xmin=0, xmax=500, color='darkgrey', linestyle='solid', linewidth=0.25, label="water line")
        plt.plot([pos1_2014, pos2_2014], [y1_2014, y2_2014], color="coral", linewidth=0.5)
        plt.scatter([pos1_2014, pos2_2014], [y1_2014, y2_2014], color="coral", s=15)
        plt.plot([pos1_2020, pos2_2020], [y1_2020, y2_2020], color="forestgreen", linewidth=0.5)
        plt.scatter([pos1_2020, pos2_2020], [y1_2020, y2_2020], color="forestgreen", s=15)
        # plot apprxomate berm elevation
        plt.hlines(y=berm_elevs[d], xmin=0, xmax=500, color='goldenrod', linestyle='solid', linewidth=0.5, label="berm")
        # plot labels
        plt.text(10, 4.5, "2014 beach slope = {0}".format(beach_slope_2014))
        plt.text(10, 4.0, "2020 beach slope = {0}".format(beach_slope_2020))
        # plt.title("transect: {0}, domain: {1}".format(t, domains[d]))
        plt.title("domain: {0}".format(domains[d]))
        plt.ylabel("elevation [m NAVD88]")
        plt.xlabel("alongshore position [m]")
        plt.ylim(-5, 5)
        plt.xlim(0, 500)
        plt.legend(loc="upper right")
        # minor ticks just to try to identify berms
        if plot_grid:
            plt.minorticks_on()
            plt.grid(visible=True, which='minor', axis="y", color="lightgrey", linewidth=1)
            plt.grid(visible=True, which='major', axis="y", color="darkgrey", linewidth=1)
        # save and close figure
        if save_on:
            plt.savefig(os.path.join(savedir, "transect_{0}_{1}.png".format(t, version)))
        if t!=17:
            plt.close()
        # move to next domain
        d += 1

# plot the beach slopes per domain
fig2 = plt.figure(figsize=(15, 6))
# 2014
plt.scatter(domains, beach_slopes_2014, label=2014, color="tab:orange")
avg_slope_2014 = round(np.mean(beach_slopes_2014),3)
plt.hlines(y=avg_slope_2014, xmin=min(domains), xmax=max(domains), linestyle="dashed", linewidth=0.75, color="tab:orange")
# 2020
plt.scatter(domains, beach_slopes_2020, label=2020, color="tab:green")
avg_slope_2020 = round(np.mean(beach_slopes_2020),3)
plt.hlines(y=avg_slope_2020, xmin=min(domains), xmax=max(domains), linestyle="dashed", linewidth=0.75, color="tab:green")
# plot labels
plt.ylabel("beach slopes")
plt.xlabel("domain")
plt.ylim(0, 0.11)
plt.legend()
plt.text(min(domains)-0.75, 0.011, "2014 average slope = {0}".format(avg_slope_2014), color="orangered")
plt.text(min(domains)-0.75, 0.005, "2020 average slope = {0}".format(avg_slope_2020), color="tab:green")
plt.xticks(ticks=domains, labels=domains)
plt.grid(visible=True, axis="x", color="lightgrey", linewidth=0.5, linestyle="dotted")
if save_slopes_on:
    plt.savefig(os.path.join(savedir, "domain_beach_slopes_{0}.png".format(version)))
plt.show()

# plot the berm elevation per domain
fig3 = plt.figure(figsize=(15, 6))
# berm points
plt.scatter(domains, berm_elevs, color="black")
avg_elev = round(np.mean(berm_elevs),1)
plt.hlines(y=avg_elev, xmin=min(domains), xmax=max(domains), linestyle="dashed", linewidth=0.75, color="darkgrey")
# plot labels
plt.ylabel("berm elevation (m NAVD88)")
plt.xlabel("domain")
plt.ylim(0, 4)
# plt.legend()
plt.text(min(domains)-0.75, 0.1, "average berm elev = {0} m NAVD88".format(avg_elev), color="black")
plt.xticks(ticks=domains, labels=domains)
plt.grid(visible=True, axis="x", color="lightgrey", linewidth=0.5, linestyle="dotted")
if save_berms_on:
    plt.savefig(os.path.join(savedir, "domain_berm_elevs_{0}.png".format(version)))
plt.show()