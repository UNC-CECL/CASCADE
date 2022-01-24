# For this class, we want to create bay to ocean flow and sediment transport
# Use Murray/Paola to carve channels through dune gaps and transport sediment?
import numpy as np

# --------------------------------------------------------------------------------------------------------------------
# 1. increase bay level by standard amount which will also change the interior domain in B3d
# this will only need to happen once, not every time step
# combine the dune domain with the interior domain and flip them so that the dune is at the back
def combine_domain(interior_domain, dune_domain):
    whole_domain = np.append(interior_domain, dune_domain, axis=0)
    flipped_domain = np.flip(whole_domain, 0)
    return flipped_domain

# this and everything after will need to update every time step
def bay_surge(flipped_domain, bay_level):
    flipped_domain = flipped_domain - bay_level
    return flipped_domain, bay_level

# --------------------------------------------------------------------------------------------------------------------
# 2. route water through interior: create channels in the interior domain to/through the dune domain using Murray/Paola
# seems like in the paper they track bed elevation and sediment transport (outputs of this function?)
# will we need to add some kind of force to get water to flow uphill, or will we need to submerge the entire interior?

def flow_routing(flipped_domain, barrier3d, Iow, DuneDomainCrest, Dow):
    # Determine Sediment And Water Routing Rules
    start = TSloc[0, 0]
    stop = TSloc[-1, 0] + 1
    Rlow = barrier3d._StormSeries[start:stop, 2]
    dur = np.array(barrier3d._StormSeries[start:stop, 4], dtype="int")
    gaps = barrier3d.DuneGaps(DuneDomain=flipped_domain, Dow=Dow, bermel=barrier3d.BermEl, Rhigh=barrier3d._StormSeries[start:stop, 1])
    if (
            len(gaps) > 0 and Iow / len(gaps) >= barrier3d._threshold_in
    ):  # If greater than threshold % of dune gaps are inunundation regime, use inun. regime routing
        inundation = 1
        substep = barrier3d._OWss_i
        barrier3d._InundationCount += 1
    else:
        inundation = 0
        substep = barrier3d._OWss_r
        barrier3d._RunUpCount += 1

    # Set Domain
    add = 10
    duration = dur[n] * substep
    width = (
            np.shape(barrier3d._InteriorDomain)[0] + 1 + add
    )  # (dam) Add one for Dunes (really a row for setting water elevation and 25 for bay)
    Elevation = np.zeros([duration, width, barrier3d._BarrierLength])
    Dunes_prestorm = DuneDomainCrest
    Dunes = Dunes_prestorm + barrier3d._BermEl
    Bay = np.ones([add, barrier3d._BarrierLength]) * -barrier3d._BayDepth
    Elevation[0, :, :] = np.vstack([Dunes, barrier3d._InteriorDomain, Bay])

    # Initialize Memory Storage Arrays
    Discharge = np.zeros([duration, width, barrier3d._BarrierLength])
    SedFluxIn = np.zeros([duration, width, barrier3d._BarrierLength])
    SedFluxOut = np.zeros([duration, width, barrier3d._BarrierLength])

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
            Rin = barrier3d._Rin_i

            # # Find average slope of interior
            # AvgSlope = self._BermEl / InteriorWidth_Avg
            #
            # # Enforce max average interior slope
            # AvgSlope = min(self._MaxAvgSlope, AvgSlope)

            # Representative average slope of interior (made static - represent. of 200-m wide barrier)
            AvgSlope = barrier3d._BermEl / 20

            C = barrier3d._Cx * AvgSlope  # Momentum constant

        else:  # Run-up regime
            Rin = barrier3d._Rin_r

    # ### Run Flow Routing Algorithm
    Hd_TSloss = (DuneChange.max(axis=1) / dur[n])

    for TS in range(duration):

        ShrubDomainWidth = np.shape(barrier3d._ShrubDomainFemale)[0]
        DeadDomainWidth = np.shape(barrier3d._ShrubDomainDead)[0]

        # this will need to be changed because we are not reducing dune height first
        # we will be finding gaps and routing water that way
        if TS > 0:
            Elevation[TS, 1:, :] = Elevation[
                                   TS - 1, 1:, :
                                   ]  # Begin timestep with elevation from end of last
            Elevation[TS, 0, :] = Dunes - (
                    Hd_TSloss / substep * TS
            )  # Reduce dune in height linearly over course of storm

        for d in range(width - 1):
            # Reduce discharge across row via infiltration
            if d > 0:
                Discharge[TS, d, :][Discharge[TS, d, :] > 0] -= Rin
            Discharge[TS, d, :][Discharge[TS, d, :] < 0] = 0

            for i in range(barrier3d._BarrierLength):
                if Discharge[TS, d, i] > 0:

                    Q0 = Discharge[TS, d, i]

                    # ### Calculate Slopes
                    if i > 0:
                        S1 = (
                                     Elevation[TS, d, i]
                                     - Elevation[TS, d + 1, i - 1]
                             ) / (np.sqrt(2))
                        S1 = np.nan_to_num(S1)
                    else:
                        S1 = 0

                    S2 = Elevation[TS, d, i] - Elevation[TS, d + 1, i]
                    S2 = np.nan_to_num(S2)

                    if i < (barrier3d._BarrierLength - 1):
                        S3 = (
                                     Elevation[TS, d, i]
                                     - Elevation[TS, d + 1, i + 1]
                             ) / (np.sqrt(2))
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
                                * S1 ** barrier3d._nn
                                / (
                                        S1 ** barrier3d._nn
                                        + S2 ** barrier3d._nn
                                        + S3 ** barrier3d._nn
                                )
                        )
                        Q2 = (
                                Q0
                                * S2 ** barrier3d._nn
                                / (
                                        S1 ** barrier3d._nn
                                        + S2 ** barrier3d._nn
                                        + S3 ** barrier3d._nn
                                )
                        )
                        Q3 = (
                                Q0
                                * S3 ** barrier3d._nn
                                / (
                                        S1 ** barrier3d._nn
                                        + S2 ** barrier3d._nn
                                        + S3 ** barrier3d._nn
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
                        if S3 == 0 and i < (barrier3d._BarrierLength - 1):
                            Q3 = Qx
                        else:
                            Q3 = 0

                    # All slopes negative
                    else:

                        Q1 = (
                                Q0
                                * abs(S1) ** (-barrier3d._nn)
                                / (
                                        abs(S1) ** (-barrier3d._nn)
                                        + abs(S2) ** (-barrier3d._nn)
                                        + abs(S3) ** (-barrier3d._nn)
                                )
                        )
                        Q2 = (
                                Q0
                                * abs(S2) ** (-barrier3d._nn)
                                / (
                                        abs(S1) ** (-barrier3d._nn)
                                        + abs(S2) ** (-barrier3d._nn)
                                        + abs(S3) ** (-barrier3d._nn)
                                )
                        )
                        Q3 = (
                                Q0
                                * abs(S3) ** (-barrier3d._nn)
                                / (
                                        abs(S1) ** (-barrier3d._nn)
                                        + abs(S2) ** (-barrier3d._nn)
                                        + abs(S3) ** (-barrier3d._nn)
                                )
                        )

                        Q1 = np.nan_to_num(Q1)
                        Q2 = np.nan_to_num(Q2)
                        Q3 = np.nan_to_num(Q3)

                        # MaxUpSlope = 0.25  # dam

                        if S1 > barrier3d._MaxUpSlope:
                            Q1 = 0
                        else:
                            Q1 = Q1 * (1 - (abs(S1) / barrier3d._MaxUpSlope))

                        if S2 > barrier3d._MaxUpSlope:
                            Q2 = 0
                        else:
                            Q2 = Q2 * (1 - (abs(S2) / barrier3d._MaxUpSlope))

                        if S3 > barrier3d._MaxUpSlope:
                            Q3 = 0
                        else:
                            Q3 = Q3 * (1 - (abs(S3) / barrier3d._MaxUpSlope))

                    # ### Reduce Overwash Through Shrub Cells and Save Discharge
                    if barrier3d._Shrub_ON:
                        # Cell 1
                        if i > 0:
                            if (
                                    d < ShrubDomainWidth
                                    and barrier3d._ShrubPercentCover[d, i - 1]
                                    > 0
                            ):
                                Q1 = Q1 * (
                                        (1 - barrier3d._Qshrub_max)
                                        * barrier3d._ShrubPercentCover[d, i - 1]
                                )
                            elif (
                                    d < DeadDomainWidth
                                    and barrier3d._DeadPercentCover[d, i - 1] > 0
                            ):
                                Q1 = Q1 * (
                                        (1 - barrier3d._Qshrub_max * 0.66)
                                        * barrier3d._DeadPercentCover[d, i - 1]
                                )  # Dead shrubs block 2/3 of what living shrubs block
                            Discharge[TS, d + 1, i - 1] = (
                                    Discharge[TS, d + 1, i - 1] + Q1
                            )

                        # Cell 2
                        # (ShrubDomainWidth - 1)?
                        if (
                                d < ShrubDomainWidth
                                and barrier3d._ShrubPercentCover[d, i] > 0
                        ):
                            Q2 = Q2 * (
                                    (1 - barrier3d._Qshrub_max)
                                    * barrier3d._ShrubPercentCover[d, i]
                            )
                        elif (
                                d < DeadDomainWidth
                                and barrier3d._DeadPercentCover[d, i] > 0
                        ):
                            Q2 = Q2 * (
                                    (1 - barrier3d._Qshrub_max * 0.66)
                                    * barrier3d._DeadPercentCover[d, i]
                            )
                        Discharge[TS, d + 1, i] = (
                                Discharge[TS, d + 1, i] + Q2
                        )

                        # Cell 3
                        if i < (barrier3d._BarrierLength - 1):
                            if (
                                    d < ShrubDomainWidth
                                    and barrier3d._ShrubPercentCover[d, i + 1]
                                    > 0
                            ):
                                Q3 = Q3 * (
                                        (1 - barrier3d._Qshrub_max)
                                        * barrier3d._ShrubPercentCover[d, i + 1]
                                )
                            elif (
                                    d < DeadDomainWidth
                                    and barrier3d._DeadPercentCover[d, i + 1] > 0
                            ):
                                Q3 = Q3 * (
                                        (1 - barrier3d._Qshrub_max * 0.66)
                                        * barrier3d._DeadPercentCover[d, i + 1]
                                )
                            Discharge[TS, d + 1, i + 1] = (
                                    Discharge[TS, d + 1, i + 1] + Q3
                            )
                    else:
                        # Cell 1
                        if i > 0:
                            Discharge[TS, d + 1, i - 1] = (
                                    Discharge[TS, d + 1, i - 1] + Q1
                            )

                        # Cell 2
                        Discharge[TS, d + 1, i] = (
                                Discharge[TS, d + 1, i] + Q2
                        )

                        # Cell 3
                        if i < (barrier3d._BarrierLength - 1):
                            Discharge[TS, d + 1, i + 1] = (
                                    Discharge[TS, d + 1, i + 1] + Q3
                            )

                    # ### Calculate Sed Movement
                    fluxLimit = barrier3d._Dmax

                    # Run-up Regime
                    if inundation == 0:
                        if Q1 > barrier3d._Qs_min and S1 >= 0:
                            Qs1 = barrier3d._Kr * Q1
                            if Qs1 > fluxLimit:
                                Qs1 = fluxLimit
                        else:
                            Qs1 = 0

                        if Q2 > barrier3d._Qs_min and S2 >= 0:
                            Qs2 = barrier3d._Kr * Q2
                            if Qs2 > fluxLimit:
                                Qs2 = fluxLimit
                        else:
                            Qs2 = 0

                        if Q3 > barrier3d._Qs_min and S3 >= 0:
                            Qs3 = barrier3d._Kr * Q3
                            if Qs3 > fluxLimit:
                                Qs3 = fluxLimit
                        else:
                            Qs3 = 0

                    # Inundation Regime - Murray and Paola (1994, 1997) Rule 3 with flux limiter
                    else:
                        if Q1 > barrier3d._Qs_min:
                            Qs1 = barrier3d._Ki * (Q1 * (S1 + C)) ** barrier3d._mm
                            if Qs1 < 0:
                                Qs1 = 0
                            elif Qs1 > fluxLimit:
                                Qs1 = fluxLimit
                        else:
                            Qs1 = 0

                        if Q2 > barrier3d._Qs_min:
                            Qs2 = barrier3d._Ki * (Q2 * (S2 + C)) ** barrier3d._mm
                            if Qs2 < 0:
                                Qs2 = 0
                            elif Qs2 > fluxLimit:
                                Qs2 = fluxLimit
                        else:
                            Qs2 = 0

                        if Q3 > barrier3d._Qs_min:
                            Qs3 = barrier3d._Ki * (Q3 * (S3 + C)) ** barrier3d._mm
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
                    if Elevation[TS, d, i] > barrier3d._SL or any(
                            z > barrier3d._SL
                            for z in Elevation[TS, d + 1: d + 10, i]
                    ):  # If cell is subaerial, elevation change is determined by difference between
                        # flux in vs. flux out
                        if i > 0:
                            SedFluxIn[TS, d + 1, i - 1] += Qs1

                        SedFluxIn[TS, d + 1, i] += Qs2

                        if i < (barrier3d._BarrierLength - 1):
                            SedFluxIn[TS, d + 1, i + 1] += Qs3

                        Qs_out = Qs1 + Qs2 + Qs3
                        SedFluxOut[TS, d, i] = Qs_out

                    else:  # If cell is subaqeous, exponentially decay dep. of remaining sed across bay

                        if inundation == 0:
                            Cbb = barrier3d._Cbb_r
                        else:
                            Cbb = barrier3d._Cbb_i

                        Qs0 = SedFluxIn[TS, d, i] * Cbb

                        Qs1 = Qs0 * Q1 / (Q1 + Q2 + Q3)
                        Qs2 = Qs0 * Q2 / (Q1 + Q2 + Q3)
                        Qs3 = Qs0 * Q3 / (Q1 + Q2 + Q3)

                        Qs1 = np.nan_to_num(Qs1)
                        Qs2 = np.nan_to_num(Qs2)
                        Qs3 = np.nan_to_num(Qs3)

                        if Qs1 < barrier3d._Qs_bb_min:
                            Qs1 = 0
                        elif Qs1 > fluxLimit:
                            Qs1 = fluxLimit
                        if Qs2 < barrier3d._Qs_bb_min:
                            Qs2 = 0
                        elif Qs2 > fluxLimit:
                            Qs2 = fluxLimit
                        if Qs3 < barrier3d._Qs_bb_min:
                            Qs3 = 0
                        elif Qs3 > fluxLimit:
                            Qs3 = fluxLimit

                        if i > 0:
                            SedFluxIn[TS, d + 1, i - 1] += Qs1

                        SedFluxIn[TS, d + 1, i] += Qs2

                        if i < (barrier3d._BarrierLength - 1):
                            SedFluxIn[TS, d + 1, i + 1] += Qs3

                        Qs_out = Qs1 + Qs2 + Qs3
                        SedFluxOut[TS, d, i] = Qs_out

                    # ### Saline Flooding
                    (
                        barrier3d._ShrubDomainFemale,
                        barrier3d._ShrubDomainMale,
                        barrier3d._ShrubDomainDead,
                    ) = barrier3d.SalineFlooding(
                        ShrubDomainWidth,
                        barrier3d._ShrubDomainAll,
                        barrier3d._ShrubDomainFemale,
                        barrier3d._ShrubDomainMale,
                        barrier3d._ShrubDomainDead,
                        d,
                        i,
                        Q0,
                    )

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

        # Update amount of burial/erosion for each shrub
        barrier3d._BurialDomain = barrier3d.UpdateBurial(
            barrier3d._BurialDomain,
            ElevationChange,
            ShrubDomainWidth,
            barrier3d._ShrubDomainAll,
        )
# --------------------------------------------------------------------------------------------------------------------
# 3. track sediment transport at the dune gaps using Nienhuis?
# only at the dune gap itself and does not alter the gap width

# --------------------------------------------------------------------------------------------------------------------
# 4. add a beach and deposit sediment (in different ways to see what happens?)

# other questions:
# will we allow any migration?
# will this be added into cascade so that it only occurs at a specific time interval?

# then we will make the class

# class Outwasher:
#     def __init__(self):  # add anything that will be an attribute to the class
#         self.____
#
#     def update(self, b3d):  # use the above functions to get outputs, will likely be using some b3d variables
#        self.___
