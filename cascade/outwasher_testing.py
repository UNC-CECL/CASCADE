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
def outwasher(b3d, storm_series, runID):
    # set Barrier3D variables
    cx = b3d._Cx                                        # multiplier with the average slope of the interior for constant "C" in inundation transport rule
    berm_el = b3d._BermEl                               # [dam MHW]
    beach_elev = b3d._BermEl                            # [dam MHW]
    length = b3d._BarrierLength                         # [dam] length of barrier
    substep = b3d._OWss_i                               # 2 for inundation
    # substep = 1
    nn = b3d._nn                                        # flow routing constant
    max_slope = -b3d._MaxUpSlope                         # max slope that sediment can go uphill, previously Slim (0.25)
    ki = b3d._Ki                                        # sediment transport coefficient
    mm = b3d._mm                                        # inundation overwash coefficient
    sea_level = b3d._SL                                 # equal to 0 dam
    q_min = b3d._Qs_min                                 # [m^3 / hr]? Minimum discharge needed for sediment transport (0.001)
    # bay_depth = -b3d._BayDepth                          # [dam MHW] Depth of bay behind island segment, currently set to 0.3
    Si = (np.mean(b3d.InteriorDomain[-10, :]) - np.mean(b3d.InteriorDomain[0, :])) / len(b3d.InteriorDomain)
    avg_slope = b3d._BermEl / 20                        # how it is defined in barrier 3D which is much smaller than when
                                                        # you calculate the slope using the avg of the first and last rows
    # setting up dune domain using b3d
    dune_domain = b3d.DuneDomain
    dune_width = b3d._DuneWidth
    dune_crest = b3d._DuneDomain[b3d._time_index, :, :].max(axis=1)  # not sure about this it is 0.0075
    # dune_crest used to be DuneDomainCrest
    dune_restart = b3d._DuneRestart                     # currently set to 0.0075 also
    max_dune = b3d._Dmaxel - b3d._BermEl                # [dam MHW]
    # max_dune_elev = b3d._Dmaxel
    Hd_avgTS = b3d._Hd_AverageTS
    dune_crest[dune_crest < dune_restart] = dune_restart
    Hd_avgTS.append(
        np.mean(dune_crest)
    )  # Store average pre-storm dune-height for time step
    time_b3d = b3d._time_index
    c1 = b3d._C1
    c2 = b3d._C2
    Hd_loss_TS = b3d._Hd_Loss_TS
    Rhigh = max(storm_series[1])

    # Set other variables
    qs_lost_total = 0  # previously OWloss
    # numstorm = int(len(storm_series))
    numstorm = 1
    bay_depth = storm_series[1][0]  # should be the first bay elevation which is 0

    # initializing our barrier interior
    interior_domain = np.zeros([30, length])
    for row in range(len(interior_domain)):
        if row == 0:
            interior_domain[row, :] = 0
        else:
            # interior_domain[row, :] = interior_domain[row-1, :] + avg_slope  # giving the back barrier an equal slope
            interior_domain[row, :] = interior_domain[row-1, :] + -Si  # giving the back barrier an equal slope

    if numstorm > 0:
        # ### Individual Storm Impacts
        for n in range(numstorm):
            # setting up the storm series for substep of 2
            if substep == 2:
                dur = 2*storm_series[2]
                sound_level = np.zeros(dur)  # initializing the new bay level array
                avgs = []
                # we need to double the duration of the storm with the same magnitudes as before
                # making the avg water level array
                count = 0
                for time in range(len(storm_series[1]) - 1):
                    avg_level = (storm_series[1][time] + storm_series[1][time + 1]) / 2
                    avgs.append(avg_level)
                # for the new time series, if the time step is even, we keep the value from the original time series
                # if the time step is odd, we take the value from the avgs array
                for time2 in range(dur):
                    if time2 % 2 == 0:
                        sound_level[time2] = storm_series[1][int(time2 / 2)]
                    elif time2 == dur - 1:
                        sound_level[time2] = avgs[-1]
                    else:
                        count += 1
                        sound_level[time2] = avgs[time2 - count]
                        # if we are at time 1, need to access avgs time 0 (time - 1)
                        # if we are at time 3, need to access avgs time 1 (time - 2)
                        # if we are at time 5, need to access avgs time 2 (time - 3)
                storm_series[1] = sound_level
            else:
                dur = storm_series[2]  # [hr] duration of the storm

            # ### Dune Erosion
            DuneChange = np.zeros(
                [length, dune_width]
            )  # Vector storing dune height change for this storm

            # Find overwashed dunes and gaps
            # currently changed Rhigh[n] to max of the bay depth in the storm series, but this might have to go into the
            # TS loop so that Rhigh is set to whatever the hydrograph is at for that time step
            Dow = [index for index, value in enumerate((dune_crest + berm_el))
                if value < Rhigh]  # bay depth used to be Rhigh[n]
            gaps = b3d.DuneGaps(
                dune_crest, Dow, berm_el, Rhigh
            )  # Finds location and Rexcess of continuous gaps in dune ridge, baydepth used to be Rhigh[n]

            for d in range(len(Dow)):  # Loop through each overwashed dune cell
                for w in range(dune_width):

                    # Calculate dune elevation loss
                    Rnorm = Rhigh / (
                            dune_domain[time_b3d, Dow[d], w] + berm_el
                    )  # Rhigh relative to pre-storm dune elevation
                    Dloss = Rnorm / (
                            c1 + (Rnorm * (Rnorm - c2))
                    )  # Amount of dune crest elevation change normalized by pre-storm dune elevation
                    # (i.e. a percent change), from Goldstein and Moore (2016)

                    # Set new dune height
                    InitDElev = (
                            dune_domain[time_b3d, Dow[d], w] + berm_el
                    )
                    NewDElev = InitDElev * (1 - Dloss)  # Calculate new dune elevation from storm lowering
                    if NewDElev < berm_el:
                        NewDElev = berm_el
                    dune_domain[time_b3d, Dow[d], w] = (
                            NewDElev - berm_el
                    )  # Convert elevation to height above berm

                    DuneChange[Dow[d], w] = InitDElev - NewDElev

                    # If dune is lowered to ~ zero, allow for chance of regrowth by raising dune height to 5 cm
                    if dune_domain[time_b3d, Dow[d], w] < dune_restart:
                        if dune_restart < max_dune:
                            dune_domain[time_b3d, Dow[d], w] = dune_restart
                        else:
                            dune_domain[time_b3d, Dow[d], w] = (max_dune)  # Restart height can't be greater than Dmax

            # Dune Height Diffusion
            dune_domain = b3d.DiffuseDunes(dune_domain, time_b3d)
            dune_domain[dune_domain < sea_level] = dune_restart

            DuneLoss = np.sum(DuneChange) / length
            Hd_TSloss = (DuneChange.max(axis=1) / dur
            )  # Average height of dune loss for each substep during storm
            Hd_loss_TS[time_b3d, :] = Hd_loss_TS[time_b3d, :] + DuneChange.max(axis=1)

            # ### Overwash
            Iow = 0  # Count of dune gaps in inundation regime
            dunes_prestorm = dune_crest
            for q in range(len(gaps)):
                start = gaps[q][0]
                stop = gaps[q][1]
                gapwidth = stop - start + 1
                meandune = (sum(dunes_prestorm[start: stop + 1]) / gapwidth) + berm_el  # Average elevation of dune gap
                # Determine number of gaps in inundation regime
                # if Rlow[n] > meandune:
                #     Iow += 1
                Iow += 1  # all inundation regime fo rus

            # Determine Sediment And Water Routing Rules
            # ### Set Domain
            if n == 0:
                beach_domain = np.ones([7, length]) * beach_elev  # [dam MHW] 7 rows
                shoreface_domain = np.zeros([6, length])
                # we actually want the beach to have a slope, but keep the first few rows the berm elevation
                # we want the beach slope to be 0.004 m = 0.0004 dam
                for b in range(len(beach_domain)):
                    if b >= 3:
                        beach_domain[b, :] = beach_domain[b-1, :] - 0.0004
                m_shoreface = beach_domain[-1, 0]/len(shoreface_domain)
                for s in range(len(shoreface_domain)):
                    if s == 0:
                        shoreface_domain[s, :] = beach_domain[-1, 0] - m_shoreface  # slope of shoreface in b3d
                    else:
                        shoreface_domain[s, :] = shoreface_domain[s-1, :] - m_shoreface
                dune_domain_full = np.transpose(dune_domain[0]) + berm_el

                full_domain = np.append(interior_domain, dune_domain_full, 0)  # [dam MHW]
                full_domain = np.append(full_domain, beach_domain, 0)  # [dam MHW]
                full_domain = np.append(full_domain, shoreface_domain, 0)
                # full_domain[0, :, :] = np.vstack([interior_domain, dune_domain_full, beach_domain, shoreface_domain])
                np.save("C:/Users/Lexi/Documents/Research/Outwasher/full_domain", full_domain)

                # testing scenarios
                # full_domain[5:10, 20:25] = 0.5
                # full_domain[15:25, 30:40] = 0.5
                # full_domain[25:30, 5:10] = 0.5
                # full_domain[:, 20:22] = 0
                # full_domain[:, 0:2] = 0
                # full_domain[:, 38:40] = 0
                # full_domain[15, :] = 0.7  # testing wall scenario
                # full_domain[15, 25:30] = 0  # testing wall scenario
                # full_domain = full_domain[5:, :]  # testing removing the bay

                # plt.plot(range(length), dune_domain_full[0])
                # plt.title("Dune Domain")
                # plt.ylim(0, 0.25)
                # plt.ylabel("dam MHW")
                # plt.xlabel("alongshore distance (dam)")
                # plt.show()

                # plot the initial full domain before sediment movement
                fig1 = plt.figure()
                ax1 = fig1.add_subplot(111)
                mat = ax1.matshow(
                    full_domain,
                    origin="upper",
                    cmap="Greens",
                    vmin=0, vmax=0.25,
                )
                fig1.colorbar(mat)
                ax1.set_title("Initial Elevation $(dam)$")
                ax1.set_ylabel("barrier width (dam)")
                ax1.set_xlabel("barrier length (dam)")
                plt.savefig("C:/Users/Lexi/Documents/Research/Outwasher/Output/Test_years/0_domain_{0}".format(runID))
                plt.show()

                # plotting cross section
                cross_section = np.mean(full_domain, 1)
                cross_section = np.flip(cross_section)
                b3d_cross = np.mean(b3d.InteriorDomain, 1)
                fig4 = plt.figure()
                ax4 = fig4.add_subplot(111)
                ax4.plot(range(len(full_domain)), cross_section, label="pre-storm")
                # ax4.plot(range(len(b3d.InteriorDomain)), b3d_cross, label="B3D")
                # ax4.set_xlabel("barrier width from ocean to bay (dam)")
                # ax4.set_ylabel("average alongshore elevation (dam)")
                # ax4.set_title("Cross shore elevation profile\n "
                #               "shoreface slope = {0}, back barrier slope = {1}".format(round(m_shoreface,3), round(Si,3)))
                # ax4.legend()
                # plt.savefig(
                #     "C:/Users/Lexi/Documents/Research/Outwasher/Output/Test_years/cross_shore_{0}".format(runID))


            width = np.shape(full_domain)[0]  # width is the number of rows in the full domain
            # duration = dur * substep  # from previous code
            duration = dur  # we already multiplied dur by the substep above
            Elevation = np.zeros([duration, width, length])
            # elevation at the first time step is set to the full domain
            Elevation[0, :, :] = full_domain

            dunes = dunes_prestorm + berm_el

            # Initialize Memory Storage Arrays
            Discharge = np.zeros([duration, width, length])
            SedFluxIn = np.zeros([duration, width, length])
            SedFluxOut = np.zeros([duration, width, length])
            elev_change_array = np.zeros([duration, width, length])
            slopes_array = np.zeros([duration, width, length])
            rexcess_dict = {}
            qs2_array = np.zeros([duration, width, length])
            qs_lost = 0  # the array for storing the sediment leaving the last row
            # currently reset every storm, but could do full cumulative

            # get the average slope of the inerior using first and last rows of just the interior domain
            # should be negative because the first row is close to bay and last row close to "dunes"
            # Si = np.mean((interior_domain[-1, :] - interior_domain[0, :]) / 20)

            # plot the bay elevation throughout each storm with sea level and beach elevation references
            x = range(0, duration)
            sea_level_line = sea_level * np.ones(len(x))
            beach_elev_line = beach_elev * np.ones(len(x))
            dune_elev_line = max(dune_domain_full[0]) * np.ones(len(x))

            fig2 = plt.figure()
            ax2 = fig2.add_subplot(111)
            ax2.plot(x, storm_series[1], label='storm {0}'.format(n + 1))
            # if we have multiple storms, will only need to plot these once
            ax2.plot(x, sea_level_line, 'k', linestyle='dashed', label='sea level')
            ax2.plot(x, beach_elev_line, 'k', linestyle='dotted', label='beach elevation')
            ax2.plot(x, dune_elev_line, 'k', linestyle='dashdot', label='dune elevation')
            ax2.legend()
            ax2.set_xlabel("storm duration")
            ax2.set_ylabel("dam MHW")
            ax2.set_title("bay elevation over the course of each storm")
            plt.savefig("C:/Users/Lexi/Documents/Research/Outwasher/Output/Test_years/{0}_hydrograph".format(runID))


            # ### Run Flow Routing Algorithm
            for TS in range(duration):
                print(TS)
                # ### Set Water at Dune Crest
                # the velocity here assumes dune overtopping (Larson 2004), probably need to change
                # also, we think there is an error in units
                # overtop_flow can be dam^3 because we are multiplying by 1 to convert
                # our initial discharge amount starts at the first row and is later distributed down the rows/cols
                Rhigh1 = storm_series[1][TS]  # [dam]
                # Rexcess = Rhigh1 - berm_el
                Rexcess = Rhigh1 - np.mean(dune_domain_full)  # this was changed from: Rexcess = Rhigh - (Hmean + bermel)
                # Rexcess = gaps[q][2]  # (m) i think this is actually dam for us
                if Rexcess < 0:
                    overtop_vel = 0
                else:
                    overtop_vel = np.sqrt(2 * (9.8/10) * Rexcess)  # (dam/s) previously Vdune
                overtop_flow = overtop_vel * Rexcess * 3600  # (dam^2/hr), do I need to multiply by 1 dam? prev Qdune
                Discharge[TS, 0, :] = overtop_flow  # (dam^3/hr)
                rexcess_dict[round(Rexcess, 5)] = round(overtop_flow, 5)
                # Begin with elevation from previous timestep
                if TS > 0:
                    # Elevation[TS, :, :] = Elevation[TS - 1, :, :]
                    Elevation[TS, 1:, :] = Elevation[TS - 1, 1:, :]
                    # Elevation[TS, 0, :] = dunes - \
                    #                       (Hd_TSloss / substep * TS
                    #     )  # Reduce dune in height linearly over course of storm

                # Loop through the rows, excluding the last one (because there is nowhere for the water/sed to go)
                # need to include the last row to keep track of how much is leaving the system
                for d in range(width-1):
                # for d in range(width): uncomment if we are letting sedimnet out the last row
                    # Discharge for each TS, row 1, and all cols set above
                    Discharge[TS, d, :][Discharge[TS, d, :] < 0] = 0

                    # Loop through each col of the specified row
                    for i in range(length):
                        # if we have discharge, set Qo equal to that value
                        if Discharge[TS, d, i] > 0:
                            Q0 = Discharge[TS, d, i]  # (dam^3/hr)

                            # ### Calculate Slopes
                            # if we are not at the last row, do normal calculations
                            # if d != width-1: uncomment if we are letting sed out the last row
                            if i > 0:  # i = 0 means there are no cols to the left
                                S1 = (Elevation[TS, d, i] - Elevation[TS, d + 1, i - 1]) / (math.sqrt(2))
                                S1 = np.nan_to_num(S1)
                            else:
                                S1 = 0

                            S2 = Elevation[TS, d, i] - Elevation[TS, d + 1, i]
                            S2 = np.nan_to_num(S2)
                            slopes_array[TS, d, i] = S2

                            if i < (length - 1):  # i at the end length means there are no cols to the right
                                S3 = (Elevation[TS, d, i] - Elevation[TS, d + 1, i + 1]) / (math.sqrt(2))
                                S3 = np.nan_to_num(S3)
                            else:
                                S3 = 0
                            # if at the last row, apply the same slope that the beach slope has
                            # else:
                            #     if i > 0:
                            #         S1 = 0.0004 / (math.sqrt(2))
                            #         S1 = np.nan_to_num(S1)
                            #     else:
                            #         S1 = 0
                            #
                            #     S2 = 0.0004
                            #     S2 = np.nan_to_num(S2)
                            #
                            #     if i < (length - 1):  # i at the end length means there are no cols to the right
                            #         S3 = 0.0004 / (math.sqrt(2))
                            #         S3 = np.nan_to_num(S3)
                            #     else:
                            #         S3 = 0

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
                                if S1 < max_slope:
                                    Q1 = 0
                                else:
                                    Q1 = Q1 * (1 - (abs(S1) / abs(max_slope)))

                                if S2 < max_slope:
                                    Q2 = 0
                                else:
                                    Q2 = Q2 * (1 - (abs(S2) / abs(max_slope)))

                                if S3 < max_slope:
                                    Q3 = 0
                                else:
                                    Q3 = Q3 * (1 - (abs(S3) / abs(max_slope)))

                            ### Update discharge
                            # discharge is defined for the next row, so we do not need to include the last row
                            # the first row of discharge was already defined
                            # if d != width-1: # uncomment
                            # Cell 1
                            if i > 0:
                                Discharge[TS, d + 1, i - 1] = Discharge[TS, d + 1, i - 1] + Q1

                            # Cell 2
                            Discharge[TS, d + 1, i] = Discharge[TS, d + 1, i] + Q2

                            # Cell 3
                            if i < (length - 1):
                                Discharge[TS, d + 1, i + 1] = Discharge[TS, d + 1, i + 1] + Q3

                            # ### Calculate Sed Movement
                            fluxLimit = max_dune  # [dam MHW] dmaxel - bermel
                            # typically max sediment comes from dune, but bay sediment probably negligible
                            # there is a source of sediment at the dune tha can be distributed to the beach
                            # fluxLimit = 0.001  # [dam^3/hr]
                            # fluxLimit = 1  # [dam^3/hr]
                            # the Q values must be between limits set by barrier 3d
                            # all Qs in [dam^3/hr]
                            # C = cx * Si  # 10 x the avg slope (from Murray)
                            C = 0.72  # directly from barrier3d
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

                            qs2_array[TS, d, i] = Qs2

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
                                # if d != width - 1:
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
                elev_change_array[TS] = ElevationChange

                # Calculate and save volume of sediment leaving the island for every hour
                # the last row is the accumulation cell for keeping track of sed out
                qs_lost = qs_lost + sum(SedFluxIn[TS, width - 1, :]) / substep  # [dam^3] previously OWloss
                qs_lost_total = qs_lost_total + sum(SedFluxIn[TS, width - 1, :]) / substep  # [dam^3] previously OWloss
                # OWloss = OWloss + np.sum(SedFluxOut[TS, 0, :]) / substep
                # if we are just letting it go out of the cell uncomment
                # qs_lost = qs_lost + sum(SedFluxOut[TS, width-1, :]) / substep  # [dam^3] previously OWloss
                # qs_lost_total = qs_lost_total + sum(SedFluxOut[TS, width - 1, :]) / substep  # [dam^3] previously OWloss

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
            check = 1
            while check == 1:
                if all(x <= bay_depth for x in InteriorUpdate[-1, :]):
                    InteriorUpdate = np.delete(InteriorUpdate, (-1), axis=0)
                else:
                    check = 0

            # Update interior domain
            # b3d._InteriorDomain = np.flip(InteriorUpdate)
            full_domain[1:, :] = InteriorUpdate
            # Update Domain widths
            # DomainWidth = np.shape(b3d._InteriorDomain)[0]

            # not entirely sure how to update the dune domain at this moment
            # dune_domain = b3d.update_dune_domain()

            # plotting cross section
            cross_section2 = np.mean(full_domain, 1)
            cross_section2 = np.flip(cross_section2)
            b3d_cross = np.mean(b3d.InteriorDomain, 1)
            # fig4 = plt.figure()
            # ax4 = fig4.add_subplot(111)
            ax4.plot(range(len(full_domain)), cross_section2, label="post storm")
            # ax4.plot(range(len(b3d.InteriorDomain)), b3d_cross, label="B3D")
            ax4.set_xlabel("barrier width from ocean to bay (dam)")
            ax4.set_ylabel("average alongshore elevation (dam)")
            ax4.set_title("Cross shore elevation profile\n "
                          "shoreface slope = {0}, back barrier slope = {1}".format(round(m_shoreface, 3), round(Si, 3)))
            ax4.legend()
            plt.savefig(
                "C:/Users/Lexi/Documents/Research/Outwasher/Output/Test_years/cross_shore_{0}".format(runID))

            # plot post storm elevation
            fig3 = plt.figure()
            ax3 = fig3.add_subplot(111)
            mat2 = ax3.matshow(
                full_domain[:, :],
                origin="upper",
                cmap="Greens",
                vmin=0, vmax=0.25,
            )
            ax3.set_xlabel('barrier length (dam)')
            ax3.set_ylabel('barrier width (dam)')
            ax3.set_title("Elevation after storm {0} $(dam)$".format(n + 1))
            fig3.colorbar(mat2)
            plt.savefig("C:/Users/Lexi/Documents/Research/Outwasher/Output/Test_years/{0}_domain_{1}".format(n + 1, runID))

    # Record storm data
    b3d._StormCount.append(numstorm)
    return Discharge, elev_change_array, full_domain, qs_lost_total, slopes_array, rexcess_dict, qs2_array


# --------------------------------------------running outwasher---------------------------------------------------------
# importing Chris' bay data
with open(r"C:\Users\Lexi\Documents\Research\Outwasher\chris stuff\sound_data.txt", newline='') as csvfile:
    sound_data = list(csv.reader(csvfile))[0]
sound_data = [float(s)/10-0.054 for s in sound_data]  # [dam MHW] Chris' sound elevations were in m MSL,
# so converted to NAVD88 then MHW and dam
sound_data = [s+0.05 for s in sound_data]  # [dam MHW] just increasing the values
# setting all negative values to 0
sound_data = sound_data[20:]
sound_data[0] = 0
# np.save("C:/Users/Lexi/Documents/Research/Outwasher/sound_data", sound_data)

# storm series is year the storm occured, the bay elevation for every time step, and the duration of the storm
storm_series = [1, sound_data, len(sound_data)]
b3d = Barrier3d.from_yaml("C:/Users/Lexi/PycharmProjects/Barrier3d/tests/test_params/")
runID = "lowered_rexcess_b3d_C"

discharge, elev_change, domain, qs_out, slopes2, dictionary, qs2 = outwasher(b3d, storm_series, runID)

# np.save("C:/Users/Lexi/Documents/Research/Outwasher/discharge", discharge)
# plt.matshow(slopes2[1], cmap="jet_r")
# plt.title('S2 Slopes')
# plt.colorbar()
# plt.matshow(discharge[1], cmap="jet_r")
# plt.title('Discharge (dam^3/hr)')
# plt.colorbar()
# ----------------------------------making the elevation gif------------------------------------------------------------
frames = []
for i in range(2):
    filename = "C:/Users/Lexi/Documents/Research/Outwasher/Output/Test_years/" + str(i) +"_domain_{0}.png".format(runID)
    frames.append(imageio.imread(filename))
imageio.mimwrite("C:/Users/Lexi/Documents/Research/Outwasher/Output/Test_years/test_{0}.gif".format(runID), frames, format= '.gif', fps = 1)


# -------------------------------------------elevation gif--------------------------------------------------------------
def plot_ElevAnimation(elev, directory, TMAX, name):
    os.chdir(directory)
    newpath = "Output/Elevations/" + name + "//"
    if not os.path.exists(newpath):
        os.makedirs(newpath)
    os.chdir(newpath)

    for t in range(TMAX):
        AnimateDomain = elev[t]

        # Plot and save
        elevFig1 = plt.figure(figsize=(15, 7))
        ax = elevFig1.add_subplot(111)
        cax = ax.matshow(
            # AnimateDomain, origin="upper", cmap="jet_r", vmin=0, vmax=0.5,
            AnimateDomain, origin="upper", cmap="seismic",
            # vmin=-0.00005, vmax=0.00005
        )  # , interpolation='gaussian') # analysis:ignore
        ax.xaxis.set_ticks_position("bottom")
        elevFig1.colorbar(cax)
        plt.xlabel("Alongshore Distance (dam)")
        plt.ylabel("Cross-Shore Distance (dam)")
        plt.title("Elevation change (dam^3)")
        plt.tight_layout()
        timestr = "Time = " + str(t) + " hrs"
        plt.text(1, 1, timestr)
        plt.rcParams.update({"font.size": 15})
        name = "elev_" + str(t)
        elevFig1.savefig(name)  # dpi=200
        plt.close(elevFig1)

    frames = []

    for filenum in range(TMAX):
        filename = "elev_" + str(filenum) + ".png"
        frames.append(imageio.imread(filename))
    # imageio.mimsave("elev.gif", frames, fps=2)
    imageio.mimsave("elev.gif", frames, "GIF-FI")
    print("[ * elevation GIF successfully generated * ]")

# -------------------------------------------discharge gif--------------------------------------------------------------
def plot_DischargeAnimation(dis, directory, TMAX, name):
    os.chdir(directory)
    newpath = "Output/Discharges/" + name + "/"
    if not os.path.exists(newpath):
        os.makedirs(newpath)
    os.chdir(newpath)

    # for t in range(TMAX - 1):
    for t in range(TMAX):
        AnimateDomain = dis[t]

        # Plot and save
        elevFig1 = plt.figure(figsize=(15, 7))
        ax = elevFig1.add_subplot(111)
        cax = ax.matshow(
            AnimateDomain, origin="upper", cmap="jet_r",
            vmin=0, vmax=200,
        )  # , interpolation='gaussian') # analysis:ignore
        ax.xaxis.set_ticks_position("bottom")
        elevFig1.colorbar(cax)
        plt.xlabel("Alongshore Distance (dam)")
        plt.ylabel("Cross-Shore Distance (dam)")
        plt.title("Discharge (dam^3/hr)")
        plt.tight_layout()
        timestr = "Time = " + str(t) + " hrs"
        plt.text(1, 1, timestr)
        plt.rcParams.update({"font.size": 15})
        name = "dis_" + str(t)
        elevFig1.savefig(name)  # dpi=200
        plt.close(elevFig1)

    frames = []

    for filenum in range(TMAX):
        filename = "dis_" + str(filenum) + ".png"
        frames.append(imageio.imread(filename))
    # imageio.mimsave("dis.gif", frames, fps=2)
    imageio.mimsave("dis.gif", frames, "GIF-FI")
    print()
    print("[ * dishcarge GIF successfully generated * ]")

# ---------------------------------------------------slope gif----------------------------------------------------------
def plot_SlopeAnimation(slope, directory, TMAX, name):
    os.chdir(directory)
    newpath = "Output/Slopes/" + name + "/"
    if not os.path.exists(newpath):
        os.makedirs(newpath)
    os.chdir(newpath)

    for t in range(TMAX):
        AnimateDomain = slope[t]

        # Plot and save
        elevFig1 = plt.figure(figsize=(15, 7))
        ax = elevFig1.add_subplot(111)
        cax = ax.matshow(
            AnimateDomain, origin="upper", cmap="jet_r",
            # vmin=-0.005, vmax=0.05,
        )  # , interpolation='gaussian') # analysis:ignore
        ax.xaxis.set_ticks_position("bottom")
        elevFig1.colorbar(cax)
        plt.xlabel("Alongshore Distance (dam)")
        plt.ylabel("Cross-Shore Distance (dam)")
        plt.title("S2 Slopes")
        plt.tight_layout()
        timestr = "Time = " + str(t) + " hrs"
        plt.text(1, 1, timestr)
        plt.rcParams.update({"font.size": 15})
        name = "slope_" + str(t)
        elevFig1.savefig(name)  # dpi=200
        plt.close(elevFig1)

    frames = []

    for filenum in range(TMAX):
        filename = "slope_" + str(filenum) + ".png"
        frames.append(imageio.imread(filename))
    # imageio.mimsave("dis.gif", frames, fps=2)
    imageio.mimsave("slopes.gif", frames, "GIF-FI")
    print()
    print("[ * slope GIF successfully generated * ]")

# -------------------------------------------qs2 gif--------------------------------------------------------------------
def plot_SlopeAnimation(qs2, directory, TMAX, name):
    os.chdir(directory)
    newpath = "Output/Qs2/" + name + "/"
    if not os.path.exists(newpath):
        os.makedirs(newpath)
    os.chdir(newpath)

    # for t in range(TMAX - 1):
    for t in range(TMAX):
        AnimateDomain = qs2[t]

        # Plot and save
        elevFig1 = plt.figure(figsize=(15, 7))
        ax = elevFig1.add_subplot(111)
        cax = ax.matshow(
            AnimateDomain, origin="upper", cmap="jet_r",
            # vmin=-0.005, vmax=0.05,
        )  # , interpolation='gaussian') # analysis:ignore
        ax.xaxis.set_ticks_position("bottom")
        elevFig1.colorbar(cax)
        plt.xlabel("Alongshore Distance (dam)")
        plt.ylabel("Cross-Shore Distance (dam)")
        plt.title("Qs2 (dam^3)")
        plt.tight_layout()
        timestr = "Time = " + str(t) + " hrs"
        plt.text(1, 1, timestr)
        plt.rcParams.update({"font.size": 15})
        name = "qs2_" + str(t)
        elevFig1.savefig(name)  # dpi=200
        plt.close(elevFig1)

    frames = []

    for filenum in range(TMAX):
        filename = "qs2_" + str(filenum) + ".png"
        frames.append(imageio.imread(filename))
    # imageio.mimsave("dis.gif", frames, fps=2)
    imageio.mimsave("qs2.gif", frames, "GIF-FI")
    print()
    print("[ * Qs2 GIF successfully generated * ]")

# -------------------------------------------b3d domain plot------------------------------------------------------------
def plot_ModelTransects(b3d, time_step):
    plt.figure(figsize=(10, 5))
    fig = plt.subplot(1, 1, 1)
    legend_t = []

    for t in time_step:

        # Sea level
        sea_level = b3d._SL + (t * b3d._RSLR[t])

        # Create data points
        shoreface_toe_x = (b3d.x_t_TS[t] - b3d.x_t_TS[0])
        shoreface_toe_y = (sea_level - b3d.DShoreface) * 10  # m
        shoreline_x = (b3d.x_s_TS[t] - b3d.x_t_TS[0])
        shoreline_y = sea_level * 10  # m
        bay_y = (sea_level - b3d._BayDepth) * 10  # m
        end_of_bay_y = bay_y

        # if cascade.nourishments[iB3D].beach_width[t] is not None:
        #     berm_x = shoreline_x + (
        #         cascade.nourishments[iB3D].beach_width[t] / 10
        #     )  # beach width (in dam)
        # else:
        berm_x = shoreline_x + (
            int(b3d.BermEl / b3d._beta)
        )  # initial beach width (in dam)
        # end of un-shifted
        berm_y = (
            b3d._BermEl * 10
        ) + shoreline_y  # convert to meters
        dune_toe_x = berm_x
        dune_toe_y = berm_y

        v = 10  # just use 10th transect
        interior_y = b3d._DomainTS[t]  # dam MHW
        interior_y = interior_y[:, v]
        dunes_y = ( b3d._DuneDomain[t, v, :] + b3d._BermEl)  # dam MHW
        cross_barrier_y = np.insert(interior_y, 0, dunes_y)
        cross_barrier_y = (cross_barrier_y * 10) + shoreline_y  # Convert to meters with sea level rise included
        cross_barrier_x = np.arange(0, len(cross_barrier_y), 1) + dune_toe_x

        end_of_bay_x = (
            cross_barrier_x[-1] + 20
        )  # just add a buffer to the end of the plt


        x = np.hstack(
            [
                shoreface_toe_x,
                shoreline_x,
                berm_x,
                dune_toe_x,
                cross_barrier_x,
                end_of_bay_x,
            ]
        )
        y = np.hstack(
            [
                shoreface_toe_y,
                shoreline_y,
                berm_y,
                dune_toe_y,
                cross_barrier_y,
                end_of_bay_y,
            ]
        )

        # Plot
        plt.plot(x, y)
        plt.hlines(sea_level * 10, shoreface_toe_x, end_of_bay_x, colors="black")
        # NOTE: the berm elevation is relative to the MHW, so everything that relies on it is m MHW; confirmed with Ian
        # that the InteriorDomain is m MHW (i.e., m NAVD88 - MHW [in NAVD88])
        # plt.rcParams.update({"font.size": 20})
        # legend_t.append(str(t))
        plt.plot(shoreface_toe_x, shoreface_toe_y, 'r-o', label='shoreface toe')
        plt.plot(shoreline_x, shoreline_y, 'g-o', label='shoreline')
        plt.plot(berm_x, berm_y, 'b-o', label='berm')
        plt.plot(dune_toe_x, dune_toe_y, 'c-o', label='dune toe')
        plt.plot(cross_barrier_x, cross_barrier_y, 'm', label='cross barrier')
        plt.plot(end_of_bay_x, end_of_bay_y, 'k-o', label='end of bay')
        plt.legend(loc='lower right')


    plt.xlabel("Cross-shore position (dam)")
    plt.ylabel("Elevation (m MHW)")
    plt.title("Profile Evolution")
    # plt.legend(legend_t)

    return fig

# change if substep is 2
# TMAX = storm_series[2]
TMAX = 2*storm_series[2]
name = runID
dir = r"C:\Users\Lexi\Documents\Research\Outwasher"
plot_ElevAnimation(elev_change, dir, TMAX, name)
plot_DischargeAnimation(discharge, dir, TMAX, name)
plot_SlopeAnimation(slopes2, dir, TMAX, name)
plot_SlopeAnimation(qs2, dir, TMAX, name)
# time_step = [0]
# plot_ModelTransects(b3d, time_step)