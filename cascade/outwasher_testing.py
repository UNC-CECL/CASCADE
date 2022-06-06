# this is a simplified domain where we are just trying to get realistic results
# uses a simpler interior domain
# testing for right dune and beach domains
# trying to correct flow routing algorithm
# also uses sound levels from chris

import numpy as np
import math
import os
import imageio
from barrier3d import Barrier3d
import matplotlib.pyplot as plt
import imageio
import csv

# -------------------------------------------------------------------------------------------------------------------
def dunes(b3d, n_rows=1, n_gaps=2, dune_height=0.3):
    length = b3d._BarrierLength
    width = n_rows
    gap_height = np.mean(b3d.InteriorDomain[0, :])
    dune_domain = np.zeros([width, length])
    dune_domain[0, :] = dune_height
    # gap = 1
    # gap1 = 1 / n_gaps * 0.5 * length
    gap_locations = np.zeros(n_gaps)
    # gap_locations[0] = gap1
    # while gap <= n_gaps-1:
    #     gap_locations[gap] = gap_locations[gap-1]+1/n_gaps*length
    #     gap += 1
    # gap_locations = np.round(gap_locations, decimals=1)
    # gap_locations = gap_locations.astype(int)
    # random gaps
    for g in range(n_gaps):
        gap_locations[g] = np.random.randint(0, length-1)
    gap_locations = gap_locations.astype(int)
    for loc in gap_locations:
        dune_domain[0, loc] = gap_height
        dune_domain[0, loc-1] = gap_height
        dune_domain[0, loc+1] = gap_height
    dune_domain[:, :] = dune_domain[0, :]
    return dune_domain, dune_height

def outwasher(b3d, storm_series, runID, dune_domain, dune_height):
    # set Barrier3D variables
    cx = b3d._Cx                                        # multiplier with the average slope of the interior for constant "C" in inundation transport rule
    beach_elev = b3d._BermEl                            # [dam MHW]
    length = b3d._BarrierLength                         # [dam] length of barrier
    # substep = b3d._OWss_i                               # 2 for inundation
    substep = 1
    nn = b3d._nn                                        # flow routing constant
    max_slope = b3d._MaxUpSlope                         # max slope that sediment can go uphill, previously Slim
    ki = b3d._Ki                                        # sediment transport coefficient
    mm = b3d._mm                                        # inundation overwash coefficient
    sea_level = b3d._SL                                 # equal to 0 dam
    q_min = b3d._Qs_min                                 # [m^3 / hr]? Minimum discharge needed for sediment transport (0.001)
    # bay_depth = -b3d._BayDepth                          # [dam MHW] Depth of bay behind island segment, currently set to 0.3

    # Set other variables
    qs_lost_total = 0  # previously OWloss
    # numstorm = int(len(storm_series))
    numstorm = 1
    bay_depth = storm_series[1][0]  # should be the first bay elevation
    interior_domain = np.empty([30, length])
    for row in range(len(interior_domain)):
        if row == 0:
            interior_domain[row, :] = bay_depth
        else:
            interior_domain[row, :] = interior_domain[row-1, :] + 0.01  # giving the back barrier a slope of 0.01 dam

    if numstorm > 0:
        # ### Individual Storm Impacts
        for n in range(numstorm):
            dur = storm_series[2]  # [hr] duration of the storm

            # Determine Sediment And Water Routing Rules
            # ### Set Domain
            add = 10

            # output the initial full domain before sediment movement
            if n == 0:
                beach_domain = np.ones([add, length]) * beach_elev  # [dam MHW] 10 rows
                # we actually want the beach to have a slope, but keep the first few rows the berm elevation
                # we want the beach slope to be 0.04 dam
                for b in range(len(beach_domain)):
                    if b >= 3:
                        beach_domain[b, :] = beach_domain[b-1, :] - 0.04
                full_domain = np.append(interior_domain, dune_domain, 0)  # [dam MHW]
                full_domain = np.append(full_domain, beach_domain, 0)  # [dam MHW]
                plt.figure(1)
                # plt.subplot(2,2,1)
                plt.matshow(
                    full_domain,
                    origin="upper",
                    cmap="seismic",
                    vmin=-0.5, vmax=0.5,
                )
                plt.xlabel('barrier length (dam)')
                plt.ylabel('barrier width (dam)')
                plt.title("Initial Elevation $(dam)$")
                plt.colorbar()
                plt.show()
                plt.savefig("C:/Users/Lexi/Documents/Research/Outwasher/Output/Test_years/0_domain_{0}".format(runID))
                # plt.close(1)

            # testing scenarios
            # full_domain[15, :] = 0.7  # testing wall scenario
            # full_domain[15, 25:30] = 0  # testing wall scenario
            # full_domain = full_domain[5:, :]  # testing removing the bay

            width = np.shape(full_domain)[0]  # width is the number of rows in the full domain
            duration = dur * substep  # from previous code
            Elevation = np.zeros([duration, width, length])
            # elevation at the first time step is set to the full domain
            Elevation[0, :, :] = full_domain

            # Initialize Memory Storage Arrays
            Discharge = np.zeros([duration, width, length])
            SedFluxIn = np.zeros([duration, width, length])
            SedFluxOut = np.zeros([duration, width, length])
            qs_lost = 0  # the array for storing the sediment leaving the last row
            # currently reset every storm, but could do full cumulative

            # get the average slope of the inerior using first and last rows of just the interior domain
            # should be negative because the first row is close to bay and last row close to "dunes"
            Si = np.mean((interior_domain[-1, :] - interior_domain[0, :]) / 20)

            # plot the bay elevation throughout each storm with sea level and beach elevation references
            plt.figure(2)
            x = range(0, duration)
            sea_level_line = sea_level * np.ones(len(x))
            beach_elev_line = beach_elev * np.ones(len(x))
            dune_elev_line = dune_height * np.ones(len(x))
            plt.figure(2)
            # plt.subplot(2,2,3)
            plt.plot(x, storm_series[1], label='storm {0}'.format(n + 1))
            # if we have multiple storms, will only need to plot these once
            plt.plot(x, sea_level_line, 'k', linestyle='dashed', label='sea level')
            plt.plot(x, beach_elev_line, 'k', linestyle='dotted', label='beach elevation')
            plt.plot(x, dune_elev_line, 'k', linestyle='dashdot', label='dune elevation')
            plt.legend()
            plt.xlabel("storm duration")
            plt.ylabel("dam MHW")
            plt.title("bay elevation over the course of each storm")
            plt.show()
            plt.savefig("C:/Users/Lexi/Documents/Research/Outwasher/Output/Test_years/{0}_hydrograph".format(runID))
            # plt.close(2)


            # ### Run Flow Routing Algorithm
            for TS in range(duration):
                # ### Set Water at Dune Crest
                # the velocity here assumes dune overtopping (Larson 2004), probably need to change
                # also, we think there is an error in units
                # overtop_flow can be dam^3 because we are multiplying by 1 to convert
                # our initial discharge amount starts at the first row and is later distributed down the rows/cols
                Rexcess = storm_series[1][TS]  # [dam]
                if Rexcess < 0:
                    overtop_vel = 0
                else:
                    overtop_vel = np.sqrt(2 * 9.8 * (Rexcess * 10)) / 10  # (dam/s)
                overtop_flow = overtop_vel * Rexcess * 3600  # (dam^2/hr), do I need to multiply by 1 dam?
                Discharge[TS, 0, :] = overtop_flow  # (dam^3/hr)

                # Begin with elevation from previous timestep
                if TS > 0:
                    Elevation[TS, 1:, :] = Elevation[TS - 1, 1:, :]

                # Loop through the rows, excluding the last one (because there is nowhere for the water/sed to go)
                # need to include the last row to keep track of how much is leaving the system
                for d in range(width):
                    # Discharge for each TS, row 1, and all cols set above
                    Discharge[TS, d, :][Discharge[TS, d, :] < 0] = 0

                    # Loop through each col of the specified row
                    for i in range(length):
                        # if we have discharge, set Qo equal to that value
                        if Discharge[TS, d, i] > 0:
                            Q0 = Discharge[TS, d, i]  # (dam^3/hr)

                            # ### Calculate Slopes
                            # if we are not at the last row, do normal calculations
                            if d != width-1:
                                if i > 0:  # i = 0 means there are no cols to the left
                                    S1 = (Elevation[TS, d, i] - Elevation[TS, d + 1, i - 1]) / (math.sqrt(2))
                                    S1 = np.nan_to_num(S1)
                                else:
                                    S1 = 0

                                S2 = Elevation[TS, d, i] - Elevation[TS, d + 1, i]
                                S2 = np.nan_to_num(S2)

                                if i < (length - 1):  # i at the end length means there are no cols to the right
                                    S3 = (Elevation[TS, d, i] - Elevation[TS, d + 1, i + 1]) / (math.sqrt(2))
                                    S3 = np.nan_to_num(S3)
                                else:
                                    S3 = 0
                            # if at the last row, defining the start of the water as the same elevation as the bay
                            # we just need a slope for sed flux out calculations
                            else:
                                if i > 0:
                                    S1 = (Elevation[TS, d, i] - bay_depth) / (math.sqrt(2))
                                    S1 = np.nan_to_num(S1)
                                else:
                                    S1 = 0

                                S2 = Elevation[TS, d, i] - bay_depth
                                S2 = np.nan_to_num(S2)

                                if i < (length - 1):  # i at the end length means there are no cols to the right
                                    S3 = (Elevation[TS, d, i] - bay_depth) / (math.sqrt(2))
                                    S3 = np.nan_to_num(S3)
                                else:
                                    S3 = 0

                            # ### Calculate Discharge To Downflow Neighbors
                            # One or more slopes positive (we have downhill flow)
                            if S1 > 0 or S2 > 0 or S3 > 0:

                                # flow does not go uphill (when theres a downhill option)
                                if S1 < 0:
                                    S1 = 0
                                if S2 < 0:
                                    S2 = 0
                                if S3 < 0:
                                    S3 = 0

                                Q1 = (
                                        Q0
                                        * S1 ** nn
                                        / (
                                                S1 ** nn
                                                + S2 ** nn
                                                + S3 ** nn
                                        )
                                )
                                Q2 = (
                                        Q0
                                        * S2 ** nn
                                        / (
                                                S1 ** nn
                                                + S2 ** nn
                                                + S3 ** nn
                                        )
                                )
                                Q3 = (
                                        Q0
                                        * S3 ** nn
                                        / (
                                                S1 ** nn
                                                + S2 ** nn
                                                + S3 ** nn
                                        )
                                )

                                Q1 = np.nan_to_num(Q1)
                                Q2 = np.nan_to_num(Q2)
                                Q3 = np.nan_to_num(Q3)

                            # No slopes positive, one or more equal to zero
                            # np downhill slopes, but some that are 0
                            elif S1 == 0 or S2 == 0 or S3 == 0:

                                # start by counting the number (1, 2, or 3) of slopes that are 0
                                s_zero = 0
                                if S1 == 0:
                                    s_zero += 1
                                if S2 == 0:
                                    s_zero += 1
                                if S3 == 0:
                                    s_zero += 1

                                # dividing the initial discharge equally among the 0 slopes
                                Qx = Q0 / s_zero
                                Qx = np.nan_to_num(Qx)

                                if S1 == 0 and i > 0:
                                    Q1 = Qx
                                else:
                                    Q1 = 0
                                if S2 == 0:
                                    Q2 = Qx
                                else:
                                    Q2 = 0
                                if S3 == 0 and i < (length - 1):
                                    Q3 = Qx
                                else:
                                    Q3 = 0

                            # All slopes negative
                            # all uphill options (likely our case for outwasher)
                            else:

                                Q1 = (
                                        Q0
                                        * abs(S1) ** (-nn)
                                        / (
                                                abs(S1) ** (-nn)
                                                + abs(S2) ** (-nn)
                                                + abs(S3) ** (-nn)
                                        )
                                )
                                Q2 = (
                                        Q0
                                        * abs(S2) ** (-nn)
                                        / (
                                                abs(S1) ** (-nn)
                                                + abs(S2) ** (-nn)
                                                + abs(S3) ** (-nn)
                                        )
                                )
                                Q3 = (
                                        Q0
                                        * abs(S3) ** (-nn)
                                        / (
                                                abs(S1) ** (-nn)
                                                + abs(S2) ** (-nn)
                                                + abs(S3) ** (-nn)
                                        )
                                )

                                Q1 = np.nan_to_num(Q1)
                                Q2 = np.nan_to_num(Q2)
                                Q3 = np.nan_to_num(Q3)

                                # we set a maximum uphill slope that sediment can still be transported to
                                if S1 > max_slope:
                                    Q1 = 0
                                else:
                                    Q1 = Q1 * (1 - (abs(S1) / max_slope))

                                if S2 > max_slope:
                                    Q2 = 0
                                else:
                                    Q2 = Q2 * (1 - (abs(S2) / max_slope))

                                if S3 > max_slope:
                                    Q3 = 0
                                else:
                                    Q3 = Q3 * (1 - (abs(S3) / max_slope))

                            ### Update discharge
                            # discharge is defined for the next row, so we do not need to include the last row
                            # the first row of discharge was already defined
                            if d != width-1:
                                # Cell 1
                                if i > 0:
                                    Discharge[TS, d + 1, i - 1] = Discharge[TS, d + 1, i - 1] + Q1

                                # Cell 2
                                Discharge[TS, d + 1, i] = Discharge[TS, d + 1, i] + Q2

                                # Cell 3
                                if i < (length - 1):
                                    Discharge[TS, d + 1, i + 1] = Discharge[TS, d + 1, i + 1] + Q3

                            # ### Calculate Sed Movement
                            # fluxLimit = b3d._Dmax
                            # typically max sediment comes from dune, but bay sediment probably negligible
                            # fluxLimit = 0.001  # [dam^3/hr]
                            fluxLimit = 1  # [dam^3/hr]
                            # the Q values must be between limits set by barrier 3d
                            # if the flow is greater than the min value required for sed movement, then sed movement
                            # will be calculated. Then, if sed flow < 0, it's set to 0, and if it's greater than the
                            # flux limit, it is set to the flux limit

                            # all Qs in [dam^3/hr]
                            C = cx * Si  # 10 x the avg slope (from Murray)
                            if Q1 > q_min:
                                Qs1 = ki * (Q1 * (S1 + C)) ** mm
                                if Qs1 < 0:
                                    Qs1 = 0
                                elif Qs1 > fluxLimit:
                                    Qs1 = fluxLimit
                            else:
                                Qs1 = 0

                            if Q2 > q_min:
                                Qs2 = ki * (Q2 * (S2 + C)) ** mm
                                if Qs2 < 0:
                                    Qs2 = 0
                                elif Qs2 > fluxLimit:
                                    Qs2 = fluxLimit
                            else:
                                Qs2 = 0

                            if Q3 > q_min:
                                Qs3 = ki * (Q3 * (S3 + C)) ** mm
                                if Qs3 < 0:
                                    Qs3 = 0
                                elif Qs3 > fluxLimit:
                                    Qs3 = fluxLimit
                            else:
                                Qs3 = 0

                            Qs1 = np.nan_to_num(Qs1)
                            Qs2 = np.nan_to_num(Qs2)
                            Qs3 = np.nan_to_num(Qs3)

                            # ### Calculate Net Erosion/Accretion
                            # if we are at the bay, or any of the next 10 are at the bay, we should not be moving sediment
                            if Elevation[TS, d, i] <= sea_level:
                                    #or any(z < sea_level for z in Elevation[TS, d + 1: d + 10, i]):

                                Elevation[TS, d, i] = Elevation[TS, d, i]

                            else:
                                # If cell is interior, elevation change is determined by difference between
                                # flux in vs. flux out
                                # sed flux in goes to the next row, and is used for determing flux out at current row
                                # so we need a flux in for the last row, which will be its own variable
                                if d != width - 1:
                                    if i > 0:
                                        SedFluxIn[TS, d + 1, i - 1] += Qs1

                                    SedFluxIn[TS, d + 1, i] += Qs2

                                    if i < (length - 1):
                                        SedFluxIn[TS, d + 1, i + 1] += Qs3
                                # Qs1,2,3 calculated for current row
                                Qs_out = Qs1 + Qs2 + Qs3
                                SedFluxOut[TS, d, i] = Qs_out

                                # END OF DOMAIN LOOPS

                # ### Update Elevation After Every Storm Hour
                ElevationChange = (SedFluxIn[TS, :, :] - SedFluxOut[TS, :, :]) / substep
                Elevation[TS, :, :] = Elevation[TS, :, :] + ElevationChange

                # Calculate and save volume of sediment leaving the island for every hour
                qs_lost = qs_lost + sum(SedFluxOut[TS, width-1, :]) / substep  # [dam^3] previously OWloss
                qs_lost_total = qs_lost_total + sum(SedFluxOut[TS, width - 1, :]) / substep  # [dam^3] previously OWloss

            with open('C:/Users/Lexi/Documents/Research/Outwasher/Output/Test_years/sediment_tracking.txt', 'a') as f:
                if n == numstorm - 1:
                    f.write(runID)
                    f.write('\n')
                    f.write("The total sediment volume lost was {0} dam^3".format(round(qs_lost_total, 3)))
                    f.write('\n')

            print("The sediment volume lost in storm {0} was {1} dam^3".format(n+1, round(qs_lost, 3)))
            if n == numstorm-1:
                print("The total sediment volume lost was {0} dam^3".format(round(qs_lost_total, 3)))


            # ### Update Interior Domain After Every Storm
            # use the last TS, not including the first (0th) row
            InteriorUpdate = Elevation[-1, 1:, :]

            # Remove all rows of bay without any deposition from the domain
            # check = 1
            # while check == 1:
            #     if all(x <= bay_depth for x in InteriorUpdate[-1, :]):
            #         InteriorUpdate = np.delete(InteriorUpdate, (-1), axis=0)
            #     else:
            #         check = 0

            # Update interior domain
            # b3d._InteriorDomain = np.flip(InteriorUpdate)
            full_domain[1:, :] = InteriorUpdate
            # Update Domain widths
            # DomainWidth = np.shape(b3d._InteriorDomain)[0]

            # plot post storm elevation
            plt.figure(3)
            # plt.subplot(2,2,1)
            plt.matshow(
                full_domain[:, :],
                origin="upper",
                cmap="seismic",
                vmin=-0.5, vmax=0.5,
            )
            plt.xlabel('barrier length (dam)')
            plt.ylabel('barrier width (dam)')
            plt.title("Elevation after storm {0} $(dam)$".format(n + 1))
            plt.colorbar()
            plt.show()
            plt.savefig("C:/Users/Lexi/Documents/Research/Outwasher/Output/Test_years/{0}_domain_{1}".format(n + 1, runID))

    # Record storm data
    b3d._StormCount.append(numstorm)
    return Discharge, ElevationChange, full_domain, qs_lost_total


# ----------------------------------------------------------------------------------------------------------------------
# running outwasher

# importing Chris' bay data
with open(r"C:\Users\Lexi\Documents\Research\Outwasher\chris stuff\sound_data.txt", newline='') as csvfile:
    sound_data = list(csv.reader(csvfile))[0]
sound_data = [float(s)/10 for s in sound_data]  # [dam MSL] Chris' sound elevations were in m MSL, so converted to dam

# storm series is year the storm occured, the bay elevation for every time step, and the duration of the storm
storm_series = [1, sound_data, len(sound_data)]
b3d = Barrier3d.from_yaml("C:/Users/Lexi/PycharmProjects/Barrier3d/tests/test_params/")
runID = "simplified_domain"
# print(storm_series[1])
# print(len(storm_series[1]))
# plt.plot(range(len(storm_series[1])), storm_series[1])

dune_domain, dune_height = dunes(b3d, n_gaps=10)
discharge_hydro, elev_change, domain, qs_out = outwasher(b3d, storm_series, runID, dune_domain, dune_height)
cross_section = np.mean(domain, 1)
# plotting cross section
# plt.subplot(2,2,4)
plt.figure(4)
plt.plot(range(len(domain)), cross_section)
plt.xlabel("barrier width from bay to ocean (dam)")
plt.ylabel("average alongshore elevation (dam)")
plt.title("Cross shore elevation profile")
plt.show()
plt.savefig("C:/Users/Lexi/Documents/Research/Outwasher/Output/Test_years/cross_shore_{0}".format(runID))


# making the elevation gif
frames = []
for i in range(2):
    filename = "C:/Users/Lexi/Documents/Research/Outwasher/Output/Test_years/" + str(i) +"_domain_{0}.png".format(runID)
    frames.append(imageio.imread(filename))
imageio.mimwrite("C:/Users/Lexi/Documents/Research/Outwasher/Output/Test_years/test_{0}.gif".format(runID), frames, format= '.gif', fps = 1)

# ----------------------------------------------------------------------------------------------------------------------
# # discharge gif
# def plot_ElevAnimation(discharge, directory, TMAX, name):
#     os.chdir(directory)
#     newpath = "Output/" + name + "/SimFrames/"
#     if not os.path.exists(newpath):
#         os.makedirs(newpath)
#     os.chdir(newpath)
#
#     # for t in range(TMAX - 1):
#     for t in range(TMAX):
#         AnimateDomain = discharge[t]
#
#         # Plot and save
#         elevFig1 = plt.figure(figsize=(15, 7))
#         ax = elevFig1.add_subplot(111)
#         cax = ax.matshow(
#             AnimateDomain, origin="upper", cmap="jet_r", vmin=0, vmax=5000,
#         )  # , interpolation='gaussian') # analysis:ignore
#         ax.xaxis.set_ticks_position("bottom")
#         elevFig1.colorbar(cax)
#         plt.xlabel("Alongshore Distance (dam)")
#         plt.ylabel("Cross-Shore Distance (dam)")
#         plt.title("Discharge (dam^3/hr)")
#         plt.tight_layout()
#         timestr = "Time = " + str(t) + " yrs"
#         plt.text(1, 1, timestr)
#         plt.rcParams.update({"font.size": 20})
#         name = "dis_" + str(t)
#         elevFig1.savefig(name)  # dpi=200
#         plt.close(elevFig1)
#
#     frames = []
#
#     for filenum in range(TMAX - 1):
#         filename = "dis_" + str(filenum) + ".png"
#         frames.append(imageio.imread(filename))
#     imageio.mimsave("dis.gif", frames, "GIF-FI")
#     print()
#     print("[ * GIF successfully generated * ]")
#
# TMAX = storm_series[0][2]
# name = runID
# plot_ElevAnimation(discharge_hydro, r"C:\Users\Lexi\Documents\Research\Outwasher", TMAX, name)

