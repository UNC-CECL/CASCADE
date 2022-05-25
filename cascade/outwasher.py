import numpy as np
import math
import os
import imageio
from barrier3d import Barrier3d
import matplotlib.pyplot as plt
import imageio

b3d = Barrier3d.from_yaml("C:/Users/Lexi/PycharmProjects/Barrier3d/tests/test_params/")
# b3d.update()

# -------------------------------------------------------------------------------------------------------------------
# Next steps
    # currently have triangle bayhigh
    # add dunes back in?
    # find better velocity: derive something from the pressure/surface water gradient (Bernoulli)?
        # each bayhigh level has an associated gradient which drives velocity
    # check flux limit ~line 300
    # backwater effects?

# have a storm that starts at 4 m (bay is -3)

#  create synthetic storm series
#  first position is year, second is max bayhigh level (dam), third is duration (hours)
storm_series = [[1, 0.6, 65], [1, 0.8, 57], [1, 0.5, 80]]
# storm_series = [[1, 0.7, 65]]

def outwasher(b3d, storm_series):
    # Set other variables
    OWloss = 0
    numstorm = int(len(storm_series))

    # If we have storms, then we will alter the full domain
    if numstorm > 0:

        # ### Individual Storm Impacts
        for n in range(numstorm):
            bayhigh = storm_series[n][1]  # dam
            dur = storm_series[n][2]  # hr

            # Determine Sediment And Water Routing Rules
            # we have decided to always do inundation regime which is 1 in barrier3d code
            inundation = 1
            substep = b3d._OWss_i  # equals 2 for inundation
            # ### Set Domain
            add = 10
            interior_domain = np.flip(b3d._InteriorDomain, 0)  # dam MHW
            if n == 0:
                beach_domain = np.ones([add, b3d._BarrierLength]) * b3d._BermEl  # dam MHW
                full_domain = np.append(interior_domain, beach_domain, 0)  # dam MHW
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
                plt.savefig("C:/Users/Lexi/Documents/Research/Outwasher/Output/Test_years/0_domain_hydro")
            # testing scenarios
            # full_domain[15, :] = 0.7  # testing wall scenario
            # full_domain[15, 25:30] = 0  # testing wall scenario
            # full_domain = full_domain[5:, :]  # testing removing the bay

            width = np.shape(full_domain)[0] # width is the number of rows in the full domain
            duration = dur * substep  # from previous code
            Elevation = np.zeros([duration, width, b3d._BarrierLength])
            # elevation at the first time step is set to the full domain
            Elevation[0, :, :] = full_domain

            # Initialize Memory Storage Arrays
            Discharge = np.zeros([duration, width, b3d._BarrierLength])
            SedFluxIn = np.zeros([duration, width, b3d._BarrierLength])
            SedFluxOut = np.zeros([duration, width, b3d._BarrierLength])
            Q1_list = np.zeros([duration, width - 1, b3d._BarrierLength])
            Q2_list = np.zeros([duration, width - 1, b3d._BarrierLength])
            Q3_list = np.zeros([duration, width - 1, b3d._BarrierLength])

            # get the average slope of the inerior using first and last rows of just the interior domain
            # should be negative because the first row is close to bay and last row close to "dunes"
            Si = np.mean((interior_domain[-1, :] - interior_domain[0, :]) / 20)

            # ### Set Water at Dune Crest
            # the velocity here assumes dune overtopping (Larson 2004), probably need to change
            # also, we think there is an error in units
            bayhigh_TS = []
            for q in range(duration):
                if q <= duration/2:
                    m = bayhigh/(duration/2)
                    bayhigh_TS.append(m*q + 0.3)
                else:
                    m = -bayhigh/(duration/2)
                    b = -m*duration
                    bayhigh_TS.append(m*q + b + 0.3)
            plt.figure(2)
            plt.plot(range(duration), bayhigh_TS)
            plt.title("bayhigh over the course of each storm")
            if n == 0:
                leg = []
            leg.append("storm {0}".format(n+1))
            plt.legend(leg)
            C = b3d._Cx * Si  # 10 x the avg slope (from Murray)
            # overtop_flow can be dam^3 because we are multiplying by 1 to convert
            # our initial discharge amount starts at the first row and is later distributed down the rows/cols
            # Discharge[:, 0, :] = overtop_flow  # (dam^3/hr)
            # C = b3d._Cx * Si  # 10 x the avg slope (from Murray)

            # ### Run Flow Routing Algorithm
            for TS in range(duration):
                Rexcess = bayhigh_TS[TS]  # dam
                overtop_vel = np.sqrt(2 * 9.8 * (Rexcess * 10)) / 10  # (dam/s)
                overtop_flow = overtop_vel * Rexcess * 3600  # (dam^2/hr), do I need to multiply by 1 dam?
                Discharge[TS, 0, :] = overtop_flow  # (dam^3/hr)


                # Begin with elevation from previous timestep
                if TS > 0:
                    Elevation[TS, 1:, :] = Elevation[TS - 1, 1:, :]

                # Loop through the rows, excluding the last one (because there is nowhere for the water/sed to go)
                for d in range(width - 1):
                    # Discharge for each TS, row 1, and all cols set above
                    Discharge[TS, d, :][Discharge[TS, d, :] < 0] = 0

                    # Loop through each col of the specified row
                    for i in range(b3d._BarrierLength):
                        # if we have discharge, set Qo equal to that value
                        if Discharge[TS, d, i] > 0:
                            Q0 = Discharge[TS, d, i]  # (dam^3/hr)

                            # ### Calculate Slopes
                            # see drawing for more conceptual understanding of this
                            #
                            if i > 0:  # i = 0 means there are no cols to the left
                                S1 = (
                                             Elevation[TS, d, i]
                                             - Elevation[TS, d + 1, i - 1]
                                     ) / (math.sqrt(2))
                                S1 = np.nan_to_num(S1)
                            else:
                                S1 = 0

                            S2 = Elevation[TS, d, i] - Elevation[TS, d + 1, i]
                            S2 = np.nan_to_num(S2)

                            if i < (b3d._BarrierLength - 1):  # i at the end length means there are no cols to the right
                                S3 = (
                                             Elevation[TS, d, i]
                                             - Elevation[TS, d + 1, i + 1]
                                     ) / (math.sqrt(2))
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
                                        * S1 ** b3d._nn
                                        / (
                                                S1 ** b3d._nn
                                                + S2 ** b3d._nn
                                                + S3 ** b3d._nn
                                        )
                                )
                                Q2 = (
                                        Q0
                                        * S2 ** b3d._nn
                                        / (
                                                S1 ** b3d._nn
                                                + S2 ** b3d._nn
                                                + S3 ** b3d._nn
                                        )
                                )
                                Q3 = (
                                        Q0
                                        * S3 ** b3d._nn
                                        / (
                                                S1 ** b3d._nn
                                                + S2 ** b3d._nn
                                                + S3 ** b3d._nn
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
                                if S3 == 0 and i < (b3d._BarrierLength - 1):
                                    Q3 = Qx
                                else:
                                    Q3 = 0

                            # All slopes negative
                            # all uphill options (likely our case for outwasher)
                            else:

                                Q1 = (
                                        Q0
                                        * abs(S1) ** (-b3d._nn)
                                        / (
                                                abs(S1) ** (-b3d._nn)
                                                + abs(S2) ** (-b3d._nn)
                                                + abs(S3) ** (-b3d._nn)
                                        )
                                )
                                Q2 = (
                                        Q0
                                        * abs(S2) ** (-b3d._nn)
                                        / (
                                                abs(S1) ** (-b3d._nn)
                                                + abs(S2) ** (-b3d._nn)
                                                + abs(S3) ** (-b3d._nn)
                                        )
                                )
                                Q3 = (
                                        Q0
                                        * abs(S3) ** (-b3d._nn)
                                        / (
                                                abs(S1) ** (-b3d._nn)
                                                + abs(S2) ** (-b3d._nn)
                                                + abs(S3) ** (-b3d._nn)
                                        )
                                )

                                Q1 = np.nan_to_num(Q1)
                                Q2 = np.nan_to_num(Q2)
                                Q3 = np.nan_to_num(Q3)

                                # we set a maximum uphill slope that sediment can still be transported to
                                Slim = b3d._MaxUpSlope

                                if S1 > Slim:
                                    Q1 = 0
                                else:
                                    Q1 = Q1 * (1 - (abs(S1) / Slim))

                                if S2 > Slim:
                                    Q2 = 0
                                else:
                                    Q2 = Q2 * (1 - (abs(S2) / Slim))

                                if S3 > Slim:
                                    Q3 = 0
                                else:
                                    Q3 = Q3 * (1 - (abs(S3) / Slim))

                            Q1_list[TS, d, i] = Q1
                            Q2_list[TS, d, i] = Q2
                            Q3_list[TS, d, i] = Q3

                            ### Update discharge
                            # Cell 1
                            if i > 0:
                                Discharge[TS, d + 1, i - 1] = Discharge[TS, d + 1, i - 1] + Q1

                            # Cell 2
                            Discharge[TS, d + 1, i] = Discharge[TS, d + 1, i] + Q2

                            # Cell 3
                            if i < (b3d._BarrierLength - 1):
                                Discharge[TS, d + 1, i + 1] = Discharge[TS, d + 1, i + 1] + Q3

# --------------------------------------------------Clear until here----------------------------------------------------
                            # ### Calculate Sed Movement
                            # fluxLimit = b3d._Dmax
                            # typically max sediment comes from dune, but bay sediment probably negligible
                            fluxLimit = 1
                            # the Q values must be between limits set by barrier 3d
                            # if the flow is greater than the min value required for sed movement, then sed movement
                            # will be calculated. Then, if sed flow < 0, it's set to 0, and if it's greater than the
                            # flux limit, it is set to the flux limit

                            # all Qs in dam^3/hr

                            if Q1 > b3d._Qs_min:
                                Qs1 = b3d._Ki * (Q1 * (S1 + C)) ** b3d._mm
                                if Qs1 < 0:
                                    Qs1 = 0
                                elif Qs1 > fluxLimit:
                                    Qs1 = fluxLimit
                            else:
                                Qs1 = 0

                            if Q2 > b3d._Qs_min:
                                Qs2 = b3d._Ki * (Q2 * (S2 + C)) ** b3d._mm
                                if Qs2 < 0:
                                    Qs2 = 0
                                elif Qs2 > fluxLimit:
                                    Qs2 = fluxLimit
                            else:
                                Qs2 = 0

                            if Q3 > b3d._Qs_min:
                                Qs3 = b3d._Ki * (Q3 * (S3 + C)) ** b3d._mm
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
                            if Elevation[TS, d, i] <= b3d._SL or any(
                                    z < b3d._SL for z in Elevation[TS, d + 1: d + 10, i]
                            ):
                                Elevation[TS, d, i] = Elevation[TS, d, i]
                            # currently the beach elevation is set to b3d._BermEl so using that here
                            # determine interior sediment distribution by flux in and out
                            # otherwise, at the beach, sediment distributed by exponential decay
                            # elif Elevation[TS, d, i] != b3d._BermEl or any(
                            #         z != b3d._BermEl for z in Elevation[TS, d + 1: d + 10, i]):
                            elif d < 45:
                                # If cell is interior, elevation change is determined by difference between
                                # flux in vs. flux out
                                if i > 0:
                                    SedFluxIn[TS, d + 1, i - 1] += Qs1

                                SedFluxIn[TS, d + 1, i] += Qs2

                                if i < (b3d._BarrierLength - 1):
                                    SedFluxIn[TS, d + 1, i + 1] += Qs3

                                Qs_out = Qs1 + Qs2 + Qs3
                                SedFluxOut[TS, d, i] = Qs_out

                            else:  # If beach cell, exponentially decay dep. of remaining sed across beach
                                #
                                # if inundation == 0:
                                #     Cbb = b3d._Cbb_r
                                # else:
                                Cbb = b3d._Cbb_i

                                # this should be nonzero after first row i think
                                Qs0 = SedFluxIn[TS, d, i] * Cbb

                                Qs1 = Qs0 * Q1 / (Q1 + Q2 + Q3)
                                Qs2 = Qs0 * Q2 / (Q1 + Q2 + Q3)
                                Qs3 = Qs0 * Q3 / (Q1 + Q2 + Q3)

                                Qs1 = np.nan_to_num(Qs1)
                                Qs2 = np.nan_to_num(Qs2)
                                Qs3 = np.nan_to_num(Qs3)

                                if Qs1 < b3d._Qs_bb_min:
                                    Qs1 = 0
                                elif Qs1 > fluxLimit:
                                    Qs1 = fluxLimit
                                if Qs2 < b3d._Qs_bb_min:
                                    Qs2 = 0
                                elif Qs2 > fluxLimit:
                                    Qs2 = fluxLimit
                                if Qs3 < b3d._Qs_bb_min:
                                    Qs3 = 0
                                elif Qs3 > fluxLimit:
                                    Qs3 = fluxLimit

                                if i > 0:
                                    SedFluxIn[TS, d + 1, i - 1] += Qs1

                                SedFluxIn[TS, d + 1, i] += Qs2

                                if i < (b3d._BarrierLength - 1):
                                    SedFluxIn[TS, d + 1, i + 1] += Qs3

                                Qs_out = Qs1 + Qs2 + Qs3
                                SedFluxOut[TS, d, i] = Qs_out
                                # END OF DOMAIN LOOPS

                # ### Update Elevation After Every Storm Hour
                ElevationChange = (SedFluxIn[TS, :, :] - SedFluxOut[TS, :, :]) / substep
                Elevation[TS, :, :] = Elevation[TS, :, :] + ElevationChange

                # Calculate and save volume of sediment deposited on/behind the island for every hour
                # (all four methods below should equal the same!)
                # OWloss = OWloss + np.sum(ElevationChange[1:,:])
                # OWloss = OWloss + (np.sum(SedFluxIn[TS,1:,:]) - np.sum(SedFluxOut[TS,1:,:])) / substep
                # OWloss = OWloss + np.sum(SedFluxIn[TS,1,:]) / substep
                OWloss = OWloss + np.sum(SedFluxOut[TS, 0, :]) / substep

            # ### Update Interior Domain After Every Storm
            # use the last TS, not including the first (0th) row
            InteriorUpdate = Elevation[-1, 1:, :]

            # Remove all rows of bay without any deposition from the domain
            check = 1
            while check == 1:
                if all(x <= -b3d._BayDepth for x in InteriorUpdate[-1, :]):
                    InteriorUpdate = np.delete(InteriorUpdate, (-1), axis=0)
                else:
                    check = 0

            # Update interior domain
            b3d._InteriorDomain = np.flip(InteriorUpdate)
            full_domain[1:, :] = InteriorUpdate
            # Update Domain widths
            # DomainWidth = np.shape(b3d._InteriorDomain)[0]

            # plots
            plt.matshow(
                full_domain[:, :],
                origin="upper",
                cmap="seismic",
                vmin=-0.5, vmax=0.5,
            )
            plt.xlabel('barrier length (dam)')
            plt.ylabel('barrier width (dam)')
            plt.title("Elevation after storm {0} $(dam)$".format(n+1))
            plt.colorbar()
            plt.savefig("C:/Users/Lexi/Documents/Research/Outwasher/Output/Test_years/{0}_domain_hydro".format(n+1))

    # plt.matshow(Discharge[0, :, :],
    #             origin="upper",
    #             cmap='jet_r')
    # # plt.xlabel('barrier length (dam)')
    # # plt.ylabel('barrier width (dam)')
    # avg = np.mean(full_domain[0, :])
    # h_water = round(avg + bayhigh, 2)
    # plt.title("total discharge in dam$^3$/hr \n TS=0 \n height of water = {0} dam".format(h_water))
    # plt.colorbar()

    # Record storm data
    b3d._StormCount.append(numstorm)
    return Discharge, ElevationChange, full_domain

# ----------------------------------------------------------------------------------------------------------------------
# running outwasher
discharge_hydro, elev_change, domain = outwasher(b3d, storm_series)

# making the elevation gif
frames = []
for i in range(4):
    filename = "C:/Users/Lexi/Documents/Research/Outwasher/Output/Test_years/" + str(i) +"_domain_hydro.png"
    frames.append(imageio.imread(filename))
imageio.mimwrite("C:/Users/Lexi/Documents/Research/Outwasher/Output/Test_years/test_hydro.gif", frames, format= '.gif', fps = 1)

# ----------------------------------------------------------------------------------------------------------------------
# # discharge gif
# def plot_ElevAnimation(discharge, directory, TMAX, name):
#     # length = b3d[0]._BarrierLength
#
#     # BeachWidth = 6
#     # OriginY = int(b3d[0]._x_s_TS[0] - b3d[0]._x_t_TS[0])
#     # AniDomainWidth = int(
#     #     np.amax(b3d[0]._InteriorWidth_AvgTS)
#     #     + BeachWidth
#     #     + np.abs(b3d[0]._ShorelineChange)
#     #     + OriginY
#     #     + 35
#     # )
#
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
#             AnimateDomain, origin="upper", cmap="jet_r", vmin=0, vmax=12000,
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
# name = "discharge_hydro"
# plot_ElevAnimation(discharge_hydro, r"C:\Users\Lexi\Documents\Research\Outwasher", TMAX, name)

