import numpy as np
import math
from barrier3d import Barrier3d

b3d = Barrier3d.from_yaml("tests/test_params/")

# --------------------------------------------------------------------------------------------------------------------
#  set discharge- currently set to the Qdune value?
Qo = 0.5  # dam^3/hr function of time and full domain (3D array fun t, x, y)

# from the Barrier3D flow routing. Instead of having DuneCrestDomain, using interior
interior_crest = Barrier3d._InteriorDomain[Barrier3d.time_index, :, :].max(
    axis=1
)  # Maximum height of each row in InteriorDomain
interior_crest[interior_crest < Barrier3d.DuneRestart] = Barrier3d.DuneRestart

# Barrier3d.Hd_AverageTS.append(
#     np.mean(interior_crest)
# )  # Store average pre-storm dune-height for time step

OWloss = 0
DuneLoss = 0
numstorm = 0

if Barrier3d.time_index >= Barrier3d.StormStart:
    # Select number of storms for this time step from normal distribution
    TSloc = np.argwhere(Barrier3d.StormSeries[:, 0] == Barrier3d.time_index)
    numstorm = int(len(TSloc))  # analysis:ignore

    if numstorm > 0:
        start = TSloc[0, 0]
        stop = TSloc[-1, 0] + 1
        Rhigh = b3d._StormSeries[start:stop, 1]
        Rlow = b3d._StormSeries[start:stop, 2]
        dur = np.array(b3d._StormSeries[start:stop, 4], dtype="int")

        # ### Individual Storm Impacts
        for n in range(numstorm):  # Loop through each individual storm

            # ###########################################
            # ### Dune Erosion

            # Find overwashed dunes and gaps
            Dow = [
                index
                for index, value in enumerate((interior_crest + b3d._BermEl))
                if value < Rhigh[n]
            ]
            gaps = Barrier3d.DuneGaps(
                interior_crest, Dow, b3d._BermEl, Rhigh[n]
            )  # Finds location and Rexcess of continuous gaps in dune ridge

            # ###########################################
            # ### Overwash

            Iow = 0  # Count of dune gaps in inundation regime
            interior_prestorm = interior_crest
            for q in range(len(gaps)):
                start = gaps[q][0]
                stop = gaps[q][1]
                gapwidth = stop - start + 1
                meandune = (
                    sum(interior_prestorm[start : stop + 1]) / gapwidth
                ) + b3d._BermEl  # Average elevation of dune gap

                # Determine number of gaps in inundation regime
                Iow += 1
                # if Rlow[n] > meandune:
                #     Iow += 1

            # Determine Sediment And Water Routing Rules
            inundation_count = 0
            if (
                len(gaps) > 0 #and Iow / len(gaps) >= Barrier3d._threshold_in
            ):  # If greater than threshold % of dune gaps are inunundation regime, use inun. regime routing
                # we changed this to only use inundation regime rules
                inundation = 1
                substep = b3d.OWss_i
                inundation_count += 1

            # Set Domain
            duration = dur[n] * substep
            add = 10
            interior_domain = Barrier3d.InteriorDomain  # elevation cell
            bay_domain = np.ones([add, b3d._BarrierLength]) * - b3d._BayDepth
            beach_domain = np.ones([add, b3d._BarrierLength]) * b3d._BermEl
            full_domain = beach_domain.append(interior_domain, 0)
            full_domain = full_domain.append(bay_domain, 0)
            flipped_domain = np.flip(full_domain, 0)

            width = np.shape(flipped_domain)[0]
            Elevation = np.zeros([duration, width, b3d._BarrierLength])
            Elevation[0, :, :] = flipped_domain


            # Initialize Memory Storage Arrays
            # I believe the barrier length should not have changed just the distribution of rows
            Discharge = np.zeros([duration, width, b3d._BarrierLength])
            SedFluxIn = np.zeros([duration, width, b3d._BarrierLength])
            SedFluxOut = np.zeros([duration, width, b3d._BarrierLength])

            Rin = 0  # (dam^3/t) Infiltration Rate, volume of overwash flow lost per m cross-shore per time

            # Set Water at Dune Crest
            for q in range(len(gaps)):
                start = gaps[q][0]
                stop = gaps[q][1]
                Rexcess = gaps[q][2]  # (m)

                # Calculate discharge through each dune cell
                Vdune = np.sqrt(2 * 9.8 * (Rexcess * 10)) / 10  # (dam/s)
                Qdune = Vdune * Rexcess * 3600  # (dam^3/hr)

                # Set discharge at dune gap
                Discharge[:, 0, start:stop] = Qdune

                if inundation == 1:  # Inundation regime
                    Rin = b3d.Rin_i

                    # # Find average slope of interior
                    # AvgSlope = self._BermEl / InteriorWidth_Avg
                    #
                    # # Enforce max average interior slope
                    # AvgSlope = min(self._MaxAvgSlope, AvgSlope)

                    # Representative average slope of interior (made static - represent. of 200-m wide barrier)
                    #AvgSlope = Barrier3d._BermEl / 20
                    Si = 0.007  # directional slope positive bc uphill
                    Slim = b3d.MaxUpSlope
                    C = b3d.Cx * Si  # 10 x the avg slope


            # ### Run Flow Routing Algorithm
            for TS in range(duration):
                if TS > 0:
                    Elevation[TS, 1:, :] = Elevation[
                        TS - 1, 1:, :
                    ]  # Begin timestep with elevation from end of last
                    # Elevation[TS, 0, :] = Dunes - (
                    #     Hd_TSloss / substep * TS
                    # )  # Reduce dune in height linearly over course of storm

                for d in range(width - 1):
                    # Reduce discharge across row via infiltration
                    if d > 0:
                        Discharge[TS, d, :][Discharge[TS, d, :] > 0] -= Rin
                    Discharge[TS, d, :][Discharge[TS, d, :] < 0] = 0

                    for i in range(b3d._BarrierLength):
                        if Discharge[TS, d, i] > 0:

                            Q0 = Discharge[TS, d, i]

                            # ### Calculate Slopes
                            if i > 0:
                                S1 = (
                                    Elevation[TS, d, i]
                                    - Elevation[TS, d + 1, i - 1]
                                ) / (math.sqrt(2))
                                S1 = np.nan_to_num(S1)
                            else:
                                S1 = 0

                            S2 = Elevation[TS, d, i] - Elevation[TS, d + 1, i]
                            S2 = np.nan_to_num(S2)

                            if i < (b3d._BarrierLength - 1):
                                S3 = (
                                    Elevation[TS, d, i]
                                    - Elevation[TS, d + 1, i + 1]
                                ) / (math.sqrt(2))
                                S3 = np.nan_to_num(S3)
                            else:
                                S3 = 0

                            # ### Calculate Discharge To Downflow Neighbors
                            # One or more slopes positive
                            if S1 > 0 or S2 > 0 or S3 > 0:

                                if S1 < 0:
                                    S1 = 0
                                if S2 < 0:
                                    S2 = 0
                                if S3 < 0:
                                    S3 = 0

                                Q1 = (
                                    Q0
                                    * S1 ** b3d.nn
                                    / (
                                        S1 ** b3d.nn
                                        + S2 ** b3d.nn
                                        + S3 ** b3d.nn
                                    )
                                )
                                Q2 = (
                                    Q0
                                    * S2 ** b3d.nn
                                    / (
                                        S1 ** b3d.nn
                                        + S2 ** b3d.nn
                                        + S3 ** b3d.nn
                                    )
                                )
                                Q3 = (
                                    Q0
                                    * S3 ** b3d.nn
                                    / (
                                        S1 ** b3d.nn
                                        + S2 ** b3d.nn
                                        + S3 ** b3d.nn
                                    )
                                )

                                Q1 = np.nan_to_num(Q1)
                                Q2 = np.nan_to_num(Q2)
                                Q3 = np.nan_to_num(Q3)

                            # No slopes positive, one or more equal to zero
                            elif S1 == 0 or S2 == 0 or S3 == 0:

                                pos = 0
                                if S1 == 0:
                                    pos += 1
                                if S2 == 0:
                                    pos += 1
                                if S3 == 0:
                                    pos += 1

                                Qx = Q0 / pos
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
                            else:

                                Q1 = (
                                    Q0
                                    * abs(S1) ** (-b3d.nn)
                                    / (
                                        abs(S1) ** (-b3d.nn)
                                        + abs(S2) ** (-b3d.nn)
                                        + abs(S3) ** (-b3d.nn)
                                    )
                                )
                                Q2 = (
                                    Q0
                                    * abs(S2) ** (-b3d.nn)
                                    / (
                                        abs(S1) ** (-b3d.nn)
                                        + abs(S2) ** (-b3d.nn)
                                        + abs(S3) ** (-b3d.nn)
                                    )
                                )
                                Q3 = (
                                    Q0
                                    * abs(S3) ** (-b3d.nn)
                                    / (
                                        abs(S1) ** (-b3d.nn)
                                        + abs(S2) ** (-b3d.nn)
                                        + abs(S3) ** (-b3d.nn)
                                    )
                                )

                                Q1 = np.nan_to_num(Q1)
                                Q2 = np.nan_to_num(Q2)
                                Q3 = np.nan_to_num(Q3)

                                # MaxUpSlope = 0.25  # dam

                                if S1 > b3d.MaxUpSlope:
                                    Q1 = 0
                                else:
                                    Q1 = Q1 * (1 - (abs(S1) / b3d.MaxUpSlope))

                                if S2 > b3d.MaxUpSlope:
                                    Q2 = 0
                                else:
                                    Q2 = Q2 * (1 - (abs(S2) / b3d.MaxUpSlope))

                                if S3 > b3d.MaxUpSlope:
                                    Q3 = 0
                                else:
                                    Q3 = Q3 * (1 - (abs(S3) / b3d.MaxUpSlope))


                            # ### Calculate Sed Movement
                            fluxLimit = Barrier3d.Dmax
                            if Q1 > b3d.Qs_min:
                                Qs1 = b3d.Ki * (Q1 * (S1 + C)) ** b3d.mm
                                if Qs1 < 0:
                                    Qs1 = 0
                                elif Qs1 > fluxLimit:
                                    Qs1 = fluxLimit
                            else:
                                Qs1 = 0

                            if Q2 > b3d.Qs_min:
                                Qs2 = b3d.Ki * (Q2 * (S2 + C)) ** b3d.mm
                                if Qs2 < 0:
                                    Qs2 = 0
                                elif Qs2 > fluxLimit:
                                    Qs2 = fluxLimit
                            else:
                                Qs2 = 0

                            if Q3 > b3d.Qs_min:
                                Qs3 = b3d.Ki * (Q3 * (S3 + C)) ** b3d.mm
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
                            if Elevation[TS, d, i] > Barrier3d._SL or any(
                                z > Barrier3d._SL
                                for z in Elevation[TS, d + 1: d + 10, i]
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

                                Qs0 = SedFluxIn[TS, d, i] * Cbb

                                Qs1 = Qs0 * Q1 / (Q1 + Q2 + Q3)
                                Qs2 = Qs0 * Q2 / (Q1 + Q2 + Q3)
                                Qs3 = Qs0 * Q3 / (Q1 + Q2 + Q3)

                                Qs1 = np.nan_to_num(Qs1)
                                Qs2 = np.nan_to_num(Qs2)
                                Qs3 = np.nan_to_num(Qs3)

                                if Qs1 < b3d.Qs_bb_min:
                                    Qs1 = 0
                                elif Qs1 > fluxLimit:
                                    Qs1 = fluxLimit
                                if Qs2 < b3d.Qs_bb_min:
                                    Qs2 = 0
                                elif Qs2 > fluxLimit:
                                    Qs2 = fluxLimit
                                if Qs3 < b3d.Qs_bb_min:
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
            InteriorUpdate = Elevation[-1, 1:, :]

            # Remove all rows of bay without any deposition from the domain
            check = 1
            while check == 1:
                if all(x <= -b3d._BayDepth for x in InteriorUpdate[-1, :]):
                    InteriorUpdate = np.delete(InteriorUpdate, (-1), axis=0)
                else:
                    check = 0

            # Update interior domain
            Barrier3d._InteriorDomain = InteriorUpdate