import numpy as np
import math
from barrier3d import Barrier3d
from scripts import CASCADE_plotters as cplots
import matplotlib.pyplot as plt

b3d = Barrier3d.from_yaml("C:/Users/Lexi/PycharmProjects/Barrier3d/tests/test_params/")
# b3d.update()

# -------------------------------------------------------------------------------------------------------------------
# nothing happens when there are no storms
# when there are storms, go into loop
# we will have a different storm series for outwash flow

#  create synthetic storm series
#  first position is year, second is bayhigh level (dam), third is duration (hours?)
storm_series = [[1, 0.1, 40], [1, 0.2, 57], [1, 0.05, 65]]

def outwasher(b3d, storm_series):
    ### Set Domain

    # check b3d interior domain
    plt.matshow(
        b3d._InteriorDomain * 10,
        origin="lower",
        cmap="terrain",
        # vmin=-1.1,
        # vmax=4.0,
    )
    plt.colorbar()

    # check our full domain
    # plt.matshow(
    #     full_domain * 10,
    #     # np.flip(full_domain, 0) * 10,
    #     origin="upper",
    #     cmap="terrain",
    #     # vmin=-1.1,
    #     # vmax=4.0,
    # )
    # plt.colorbar()

    # Set other variables
    OWloss = 0
    numstorm = int(len(storm_series))

    # If we have storms, then we will alter the full domain
    if numstorm > 0:

        # ### Individual Storm Impacts
        for n in range(numstorm):
            bayhigh = storm_series[n][1]
            dur = storm_series[n][2]

            # Determine Sediment And Water Routing Rules
            # we have decided to always do inundation regime which is 1 in barrier3d code
            inundation = 1
            substep = b3d._OWss_i  # OWss for inundation

            ### Set Domain
            add = 10
            interior_domain = np.flip(b3d._InteriorDomain, 0)
            beach_domain = np.ones([add, b3d._BarrierLength]) * b3d._BermEl
            full_domain = np.append(interior_domain, beach_domain, 0)
            # # interior domain in b3d is ocean to bay, so flip this for outwash flow:
            # interior_domain = np.flip(b3d._InteriorDomain, 0)
            # # bay domain is 10 rows set to a specified depth from barrier3d:
            # bay_domain = np.ones([add, b3d._BarrierLength]) * - b3d._BayDepth
            # # new domain that we are adding after the interior domain set to the berm elevation in b3d:
            # beach_domain = np.ones([add, b3d._BarrierLength]) * b3d._BermEl
            # # combining the bay and interior domain:
            # bayint_domain = np.append(bay_domain, interior_domain, 0)
            # # getting the full domain (bay, interior, beach):
            # full_domain = np.append(bayint_domain, beach_domain, 0)
            # # width is the number of rows in the full domain
            width = np.shape(full_domain)[0]
            duration = dur * substep  # from previous code
            Elevation = np.zeros([duration, width, b3d._BarrierLength])
            # elevation at the first time step is set to the full domain
            Elevation[0, :, :] = full_domain

            # Initialize Memory Storage Arrays
            Discharge = np.zeros([duration, width, b3d._BarrierLength])
            SedFluxIn = np.zeros([duration, width, b3d._BarrierLength])
            SedFluxOut = np.zeros([duration, width, b3d._BarrierLength])

            # get the average slope of the inerior using first and last rows of just the interior domain
            # should be negative because the first row is close to bay and last row close to "dunes"
            Si = np.mean((interior_domain[-1, :] - interior_domain[0, :]) / 20)

            # ### Set Water at Dune Crest
            # these equations will probably need to be changed
            Rexcess = bayhigh
            overtop_vel = np.sqrt(2 * 9.8 * (Rexcess * 10)) / 10  # (dam/s)
            overtop_flow = overtop_vel * Rexcess * 3600  # (dam^3/hr)
            # our initial discharge amount starts at the first row and is later distributed down the rows/cols
            Discharge[:, 0, :] = overtop_flow
            C = b3d._Cx * Si  # 10 x the avg slope (from Murray)

            # ### Run Flow Routing Algorithm
            for TS in range(duration):
                # Begin with elevation from previous timestep
                if TS > 0:
                    Elevation[TS, 1:, :] = Elevation[
                        TS - 1, 1:, :
                    ]

                # Loop through the rows, excluding the last one (because there is nowhere for the water/sed to go)
                for d in range(width - 1):
                    # if any of the discharge values are less than 0, set them equal to 0
                    Discharge[TS, d, :][Discharge[TS, d, :] < 0] = 0

                    # Loop through each col of the specified row
                    for i in range(b3d._BarrierLength):
                        # if we have discharge, set Qo equal to that value
                        if Discharge[TS, d, i] > 0:
                            Q0 = Discharge[TS, d, i]

                            # ### Calculate Slopes
                            # see drawing for more conceptual understanding of this
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


                            # ### Calculate Sed Movement
                            # fluxLimit = b3d._Dmax
                            # typically max sediment comes from dune, but bay sediment probably negligible
                            fluxLimit = b3d._Dmaxel/10  # play with this number
                            # the Q values must be between limits set by barrier 3d
                            # if the flow is greater than the min value required for sed movement, then sed movement
                            # will be calculated. Then, if sed flow < 0, it's set to 0, and if it's greater than the
                            # flux limit, it is set to the flux limit

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

                            # check if bay
                            # once we get to land check the 10
                            # if it is less than 10 exponential decay all of it
                            # use same exp decay for beach

                            if Elevation[TS, d, i] > b3d._SL or any(
                                z > b3d._SL for z in Elevation[TS, d + 1: d + 10, i]
                            ):  # If cell is subaerial, elevation change is determined by difference between
                                # flux in vs. flux out
                                if i > 0:
                                    SedFluxIn[TS, d + 1, i - 1] += Qs1

                                SedFluxIn[TS, d + 1, i] += Qs2

                                if i < (b3d._BarrierLength - 1):
                                    SedFluxIn[TS, d + 1, i + 1] += Qs3

                                Qs_out = Qs1 + Qs2 + Qs3
                                SedFluxOut[TS, d, i] = Qs_out

                            else:  # If cell is subaqeous, exponentially decay dep. of remaining sed across bay

                                if inundation == 0:
                                    Cbb = b3d._Cbb_r
                                else:
                                    Cbb = b3d._Cbb_i

                                # when is this nonzero
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
                if inundation == 1:
                    ElevationChange = (
                        SedFluxIn[TS, :, :] - SedFluxOut[TS, :, :]
                    ) / substep
                else:
                    ElevationChange = (
                        SedFluxIn[TS, :, :] - SedFluxOut[TS, :, :]
                    ) / substep
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
            full_domain[1:, :] = InteriorUpdate

    # plots
    plt.matshow(
        #b3d.InteriorDomain * 10,
        full_domain * 10,
        origin="upper",
        cmap="terrain",
        # vmin=-1.1,
        # vmax=4.0,
    )
    plt.colorbar()


outwasher(b3d, storm_series)
