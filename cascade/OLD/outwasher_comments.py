# For this class, we want to create bay to ocean flow and sediment transport
# Use Murray/Paola to carve channels through dune gaps and transport sediment?
import numpy as np
from barrier3d import Barrier3d
b3d = Barrier3d.from_yaml("tests/test_params/")

# --------------------------------------------------------------------------------------------------------------------
# 1. increase bay level by standard amount which will also change the interior domain in B3d
# this will only need to happen once, not every time step
# combine the dune domain with the interior domain and flip them so that the dune is at the back
# def combine_domain(interior_domain, dune_domain):
#     whole_domain = np.append(interior_domain, dune_domain, axis=0)
#     flipped_domain = np.flip(whole_domain, 0)
#     return flipped_domain
#
# # this and everything after will need to update every time step
# def bay_surge(flipped_domain, bay_level):
#     flipped_domain = flipped_domain - bay_level
#     return flipped_domain, bay_level

# --------------------------------------------------------------------------------------------------------------------
# 2. route water through interior: create channels in the interior domain to/through the dune domain using Murray/Paola
# seems like in the paper they track bed elevation and sediment transport (outputs of this function?)
# will we need to add some kind of force to get water to flow uphill, or will we need to submerge the entire interior?

# Currently happening:
# for each storm: erode dunes, route water from ocean to bay, use sed transport to route sand

#within our flow routing we are going to assume inundation overwash
# 1. set inundation rules
# 2. create water domain (970 in other code)- one in front of interior domain and one behind
# ------ sensitivity model for how many cells we need in front and behind
# ------ also assume no dune rn to see how flow routing works
# 3. create a wall of discharge water from the bay (constant mag)
# ------ later we will make this into a hydrograph

# set inundation regime rules
Rin = b3d.Rin      # from paper
Si = 0.007      # directional slope positive bc uphill
# how to find/calc this? what is the elevation of interior domain?
# what does the interior domain variable hold?
Slim = 0.25     # currently set in the yaml
n = b3d.nn       # currently set in the yaml

# creating the water domains
interior_domain = Barrier3d.InteriorDomain # elevation cell

add = 10
bay_domain = np.ones([add, Barrier3d._BarrierLength]) * -Barrier3d._BayDepth
beach_domain = np.ones([add, Barrier3d._BarrierLength]) * Barrier3d._BermEl

full_domain = beach_domain.append(interior_domain, 0)
full_domain = full_domain.append(bay_domain, 0)

#  set discharge- currently set to the Qdune value?

Qo = 0.5  # dam^3/hr function of time and full domain (3D array fun t, x, y)

# water flow rules Murray
# might want to calc slopes and calc discharge to neighbors based on slopes 1043
Qi = ((Qo - Rin)*abs(Si)**-n/sum(abs(Si)**-n))*(1-(abs(Si)/Slim))

# sediment transport Murray
Ki = b3d.Ki     # from paper
C = 0.72        # order of the average slope (Murray), supposed to be 10x average slope of barrier
m = b3d.mm      # >1 usually 2.5 (Murray), from yaml
Qsi = Ki*(Qi*(Si+C))**m

# ###########################################
# ### Storms

# OWloss = 0
# DuneLoss = 0
# numstorm = 0
#
# if self._time_index >= self._StormStart:
#     # Select number of storms for this time step from normal distribution
#     TSloc = np.argwhere(self._StormSeries[:, 0] == self._time_index)
#     numstorm = int(len(TSloc))  # analysis:ignore
#
#     if numstorm > 0:
#         start = TSloc[0, 0]
#         stop = TSloc[-1, 0] + 1
#         Rhigh = self._StormSeries[start:stop, 1]
#         Rlow = self._StormSeries[start:stop, 2]
#         dur = np.array(self._StormSeries[start:stop, 4], dtype="int")
#
#         # ### Individual Storm Impacts
#         for n in range(numstorm):  # Loop through each individual storm
#
#             # ###########################################
#             # ### Dune Erosion
#
#             DuneChange = np.zeros(
#                 [self._BarrierLength, self._DuneWidth]
#             )  # Vector storing dune height change for this storm
#
#             # Find overwashed dunes and gaps
#             Dow = [
#                 index
#                 for index, value in enumerate((DuneDomainCrest + self._BermEl))
#                 if value < Rhigh[n]
#             ]
#             gaps = self.DuneGaps(
#                 DuneDomainCrest, Dow, self._BermEl, Rhigh[n]
#             )  # Finds location and Rexcess of continuous gaps in dune ridge
#
#             for d in range(len(Dow)):  # Loop through each overwashed dune cell
#                 for w in range(self._DuneWidth):
#
#                     # Calculate dune elevation loss
#                     Rnorm = Rhigh[n] / (
#                             self._DuneDomain[self._time_index, Dow[d], w]
#                             + self._BermEl
#                     )  # Rhigh relative to pre-storm dune elevation
#                     Dloss = Rnorm / (
#                             self._C1 + (Rnorm * (Rnorm - self._C2))
#                     )  # Amount of dune crest elevation change normalized by pre-storm dune elevation
#                     # (i.e. a percent change), from Goldstein and Moore (2016)
#
#                     # Set new dune height
#                     InitDElev = (
#                             self._DuneDomain[self._time_index, Dow[d], w]
#                             + self._BermEl
#                     )
#                     NewDElev = InitDElev * (
#                             1 - Dloss
#                     )  # Calculate new dune elevation from storm lowering
#                     if NewDElev < self._BermEl:
#                         NewDElev = self._BermEl
#                     self._DuneDomain[self._time_index, Dow[d], w] = (
#                             NewDElev - self._BermEl
#                     )  # Convert elevation to height above berm
#
#                     DuneChange[Dow[d], w] = InitDElev - NewDElev
#
#                     # If dune is lowered to ~ zero, allow for chance of regrowth by raising dune height to 5 cm
#                     if (
#                             self._DuneDomain[self._time_index, Dow[d], w]
#                             < self._DuneRestart
#                     ):
#                         if self._DuneRestart < self._Dmax:
#                             self._DuneDomain[
#                                 self._time_index, Dow[d], w
#                             ] = self._DuneRestart
#                         else:
#                             self._DuneDomain[
#                                 self._time_index, Dow[d], w
#                             ] = (
#                                 self._Dmax
#                             )  # Restart height can't be greater than Dmax
#
#             # Dune Height Diffusion
#             self._DuneDomain = self.DiffuseDunes(
#                 self._DuneDomain, self._time_index
#             )
#             self._DuneDomain[self._DuneDomain < self._SL] = self._DuneRestart
#
#             DuneLoss = np.sum(DuneChange) / self._BarrierLength
#             Hd_TSloss = (
#                     DuneChange.max(axis=1) / dur[n]
#             )  # Average height of dune loss for each substep during storm
#
#             self._Hd_Loss_TS[self._time_index, :] = self._Hd_Loss_TS[
#                                                     self._time_index, :
#                                                     ] + DuneChange.max(axis=1)
#
#             # ###########################################
#             # ### Overwash
#
#             Iow = 0  # Count of dune gaps in inundation regime
#             Dunes_prestorm = DuneDomainCrest
#             for q in range(len(gaps)):
#                 start = gaps[q][0]
#                 stop = gaps[q][1]
#                 gapwidth = stop - start + 1
#                 meandune = (
#                                    sum(Dunes_prestorm[start: stop + 1]) / gapwidth
#                            ) + self._BermEl  # Average elevation of dune gap
#
#                 # Determine number of gaps in inundation regime
#                 if Rlow[n] > meandune:
#                     Iow += 1
#
#             # Determine Sediment And Water Routing Rules
#
#             if (
#                     len(gaps) > 0 and Iow / len(gaps) >= self._threshold_in
#             ):  # If greater than threshold % of dune gaps are inunundation regime, use inun. regime routing
#                 inundation = 1
#                 substep = self._OWss_i
#                 self._InundationCount += 1
#             else:
#                 inundation = 0
#                 substep = self._OWss_r
#                 self._RunUpCount += 1
#
#             # Set Domain
#             add = 10
#             duration = dur[n] * substep
#             width = (
#                     np.shape(self._InteriorDomain)[0] + 1 + add
#             )  # (dam) Add one for Dunes (really a row for setting water elevation and 25 for bay)
#             Elevation = np.zeros([duration, width, self._BarrierLength])
#             Dunes = Dunes_prestorm + self._BermEl
#             Bay = np.ones([add, self._BarrierLength]) * -self._BayDepth
#             Elevation[0, :, :] = np.vstack([Dunes, self._InteriorDomain, Bay])
#
#             # Initialize Memory Storage Arrays
#             Discharge = np.zeros([duration, width, self._BarrierLength])
#             SedFluxIn = np.zeros([duration, width, self._BarrierLength])
#             SedFluxOut = np.zeros([duration, width, self._BarrierLength])
#
#             Rin = 0  # (dam^3/t) Infiltration Rate, volume of overwash flow lost per m cross-shore per time
#
#             # Set Water at Dune Crest
#             for q in range(len(gaps)):
#                 start = gaps[q][0]
#                 stop = gaps[q][1]
#                 Rexcess = gaps[q][2]  # (m)
#
#                 # Calculate discharge through each dune cell
#                 Vdune = math.sqrt(2 * 9.8 * (Rexcess * 10)) / 10  # (dam/s)
#                 Qdune = Vdune * Rexcess * 3600  # (dam^3/hr)
#
#                 # Set discharge at dune gap
#                 Discharge[:, 0, start:stop] = Qdune
#
#                 if inundation == 1:  # Inundation regime
#                     Rin = self._Rin_i
#
#                     # # Find average slope of interior
#                     # AvgSlope = self._BermEl / InteriorWidth_Avg
#                     #
#                     # # Enforce max average interior slope
#                     # AvgSlope = min(self._MaxAvgSlope, AvgSlope)
#
#                     # Representative average slope of interior (made static - represent. of 200-m wide barrier)
#                     AvgSlope = self._BermEl / 20
#
#                     C = self._Cx * AvgSlope  # Momentum constant
#
#                 else:  # Run-up regime
#                     Rin = self._Rin_r
#
#             # ### Run Flow Routing Algorithm
#             for TS in range(duration):
#
#                 ShrubDomainWidth = np.shape(self._ShrubDomainFemale)[0]
#                 DeadDomainWidth = np.shape(self._ShrubDomainDead)[0]
#
#                 if TS > 0:
#                     Elevation[TS, 1:, :] = Elevation[
#                                            TS - 1, 1:, :
#                                            ]  # Begin timestep with elevation from end of last
#                     Elevation[TS, 0, :] = Dunes - (
#                             Hd_TSloss / substep * TS
#                     )  # Reduce dune in height linearly over course of storm
#
#                 for d in range(width - 1):
#                     # Reduce discharge across row via infiltration
#                     if d > 0:
#                         Discharge[TS, d, :][Discharge[TS, d, :] > 0] -= Rin
#                     Discharge[TS, d, :][Discharge[TS, d, :] < 0] = 0
#
#                     for i in range(self._BarrierLength):
#                         if Discharge[TS, d, i] > 0:
#
#                             Q0 = Discharge[TS, d, i]
#
#                             # ### Calculate Slopes
#                             if i > 0:
#                                 S1 = (
#                                              Elevation[TS, d, i]
#                                              - Elevation[TS, d + 1, i - 1]
#                                      ) / (math.sqrt(2))
#                                 S1 = np.nan_to_num(S1)
#                             else:
#                                 S1 = 0
#
#                             S2 = Elevation[TS, d, i] - Elevation[TS, d + 1, i]
#                             S2 = np.nan_to_num(S2)
#
#                             if i < (self._BarrierLength - 1):
#                                 S3 = (
#                                              Elevation[TS, d, i]
#                                              - Elevation[TS, d + 1, i + 1]
#                                      ) / (math.sqrt(2))
#                                 S3 = np.nan_to_num(S3)
#                             else:
#                                 S3 = 0
#
#                             # ### Calculate Discharge To Downflow Neighbors
#                             # One or more slopes positive
#                             if S1 > 0 or S2 > 0 or S3 > 0:
#
#                                 if S1 < 0:
#                                     S1 = 0
#                                 if S2 < 0:
#                                     S2 = 0
#                                 if S3 < 0:
#                                     S3 = 0
#
#                                 Q1 = (
#                                         Q0
#                                         * S1 ** self._nn
#                                         / (
#                                                 S1 ** self._nn
#                                                 + S2 ** self._nn
#                                                 + S3 ** self._nn
#                                         )
#                                 )
#                                 Q2 = (
#                                         Q0
#                                         * S2 ** self._nn
#                                         / (
#                                                 S1 ** self._nn
#                                                 + S2 ** self._nn
#                                                 + S3 ** self._nn
#                                         )
#                                 )
#                                 Q3 = (
#                                         Q0
#                                         * S3 ** self._nn
#                                         / (
#                                                 S1 ** self._nn
#                                                 + S2 ** self._nn
#                                                 + S3 ** self._nn
#                                         )
#                                 )
#
#                                 Q1 = np.nan_to_num(Q1)
#                                 Q2 = np.nan_to_num(Q2)
#                                 Q3 = np.nan_to_num(Q3)
#
#                             # No slopes positive, one or more equal to zero
#                             elif S1 == 0 or S2 == 0 or S3 == 0:
#
#                                 pos = 0
#                                 if S1 == 0:
#                                     pos += 1
#                                 if S2 == 0:
#                                     pos += 1
#                                 if S3 == 0:
#                                     pos += 1
#
#                                 Qx = Q0 / pos
#                                 Qx = np.nan_to_num(Qx)
#
#                                 if S1 == 0 and i > 0:
#                                     Q1 = Qx
#                                 else:
#                                     Q1 = 0
#                                 if S2 == 0:
#                                     Q2 = Qx
#                                 else:
#                                     Q2 = 0
#                                 if S3 == 0 and i < (self._BarrierLength - 1):
#                                     Q3 = Qx
#                                 else:
#                                     Q3 = 0
#
#                             # All slopes negative
#                             else:
#
#                                 Q1 = (
#                                         Q0
#                                         * abs(S1) ** (-self._nn)
#                                         / (
#                                                 abs(S1) ** (-self._nn)
#                                                 + abs(S2) ** (-self._nn)
#                                                 + abs(S3) ** (-self._nn)
#                                         )
#                                 )
#                                 Q2 = (
#                                         Q0
#                                         * abs(S2) ** (-self._nn)
#                                         / (
#                                                 abs(S1) ** (-self._nn)
#                                                 + abs(S2) ** (-self._nn)
#                                                 + abs(S3) ** (-self._nn)
#                                         )
#                                 )
#                                 Q3 = (
#                                         Q0
#                                         * abs(S3) ** (-self._nn)
#                                         / (
#                                                 abs(S1) ** (-self._nn)
#                                                 + abs(S2) ** (-self._nn)
#                                                 + abs(S3) ** (-self._nn)
#                                         )
#                                 )
#
#                                 Q1 = np.nan_to_num(Q1)
#                                 Q2 = np.nan_to_num(Q2)
#                                 Q3 = np.nan_to_num(Q3)
#
#                                 # MaxUpSlope = 0.25  # dam
#
#                                 if S1 > self._MaxUpSlope:
#                                     Q1 = 0
#                                 else:
#                                     Q1 = Q1 * (1 - (abs(S1) / self._MaxUpSlope))
#
#                                 if S2 > self._MaxUpSlope:
#                                     Q2 = 0
#                                 else:
#                                     Q2 = Q2 * (1 - (abs(S2) / self._MaxUpSlope))
#
#                                 if S3 > self._MaxUpSlope:
#                                     Q3 = 0
#                                 else:
#                                     Q3 = Q3 * (1 - (abs(S3) / self._MaxUpSlope))
#
#                             # ### Reduce Overwash Through Shrub Cells and Save Discharge
#                             if self._Shrub_ON:
#                                 # Cell 1
#                                 if i > 0:
#                                     if (
#                                             d < ShrubDomainWidth
#                                             and self._ShrubPercentCover[d, i - 1]
#                                             > 0
#                                     ):
#                                         Q1 = Q1 * (
#                                                 (1 - self._Qshrub_max)
#                                                 * self._ShrubPercentCover[d, i - 1]
#                                         )
#                                     elif (
#                                             d < DeadDomainWidth
#                                             and self._DeadPercentCover[d, i - 1] > 0
#                                     ):
#                                         Q1 = Q1 * (
#                                                 (1 - self._Qshrub_max * 0.66)
#                                                 * self._DeadPercentCover[d, i - 1]
#                                         )  # Dead shrubs block 2/3 of what living shrubs block
#                                     Discharge[TS, d + 1, i - 1] = (
#                                             Discharge[TS, d + 1, i - 1] + Q1
#                                     )
#
#                                 # Cell 2
#                                 # (ShrubDomainWidth - 1)?
#                                 if (
#                                         d < ShrubDomainWidth
#                                         and self._ShrubPercentCover[d, i] > 0
#                                 ):
#                                     Q2 = Q2 * (
#                                             (1 - self._Qshrub_max)
#                                             * self._ShrubPercentCover[d, i]
#                                     )
#                                 elif (
#                                         d < DeadDomainWidth
#                                         and self._DeadPercentCover[d, i] > 0
#                                 ):
#                                     Q2 = Q2 * (
#                                             (1 - self._Qshrub_max * 0.66)
#                                             * self._DeadPercentCover[d, i]
#                                     )
#                                 Discharge[TS, d + 1, i] = (
#                                         Discharge[TS, d + 1, i] + Q2
#                                 )
#
#                                 # Cell 3
#                                 if i < (self._BarrierLength - 1):
#                                     if (
#                                             d < ShrubDomainWidth
#                                             and self._ShrubPercentCover[d, i + 1]
#                                             > 0
#                                     ):
#                                         Q3 = Q3 * (
#                                                 (1 - self._Qshrub_max)
#                                                 * self._ShrubPercentCover[d, i + 1]
#                                         )
#                                     elif (
#                                             d < DeadDomainWidth
#                                             and self._DeadPercentCover[d, i + 1] > 0
#                                     ):
#                                         Q3 = Q3 * (
#                                                 (1 - self._Qshrub_max * 0.66)
#                                                 * self._DeadPercentCover[d, i + 1]
#                                         )
#                                     Discharge[TS, d + 1, i + 1] = (
#                                             Discharge[TS, d + 1, i + 1] + Q3
#                                     )
#                             else:
#                                 # Cell 1
#                                 if i > 0:
#                                     Discharge[TS, d + 1, i - 1] = (
#                                             Discharge[TS, d + 1, i - 1] + Q1
#                                     )
#
#                                 # Cell 2
#                                 Discharge[TS, d + 1, i] = (
#                                         Discharge[TS, d + 1, i] + Q2
#                                 )
#
#                                 # Cell 3
#                                 if i < (self._BarrierLength - 1):
#                                     Discharge[TS, d + 1, i + 1] = (
#                                             Discharge[TS, d + 1, i + 1] + Q3
#                                     )
#
#                             # ### Calculate Sed Movement
#                             fluxLimit = self._Dmax
#
#                             # Run-up Regime
#                             if inundation == 0:
#                                 if Q1 > self._Qs_min and S1 >= 0:
#                                     Qs1 = self._Kr * Q1
#                                     if Qs1 > fluxLimit:
#                                         Qs1 = fluxLimit
#                                 else:
#                                     Qs1 = 0
#
#                                 if Q2 > self._Qs_min and S2 >= 0:
#                                     Qs2 = self._Kr * Q2
#                                     if Qs2 > fluxLimit:
#                                         Qs2 = fluxLimit
#                                 else:
#                                     Qs2 = 0
#
#                                 if Q3 > self._Qs_min and S3 >= 0:
#                                     Qs3 = self._Kr * Q3
#                                     if Qs3 > fluxLimit:
#                                         Qs3 = fluxLimit
#                                 else:
#                                     Qs3 = 0
#
#                             # Inundation Regime - Murray and Paola (1994, 1997) Rule 3 with flux limiter
#                             else:
#                                 if Q1 > self._Qs_min:
#                                     Qs1 = self._Ki * (Q1 * (S1 + C)) ** self._mm
#                                     if Qs1 < 0:
#                                         Qs1 = 0
#                                     elif Qs1 > fluxLimit:
#                                         Qs1 = fluxLimit
#                                 else:
#                                     Qs1 = 0
#
#                                 if Q2 > self._Qs_min:
#                                     Qs2 = self._Ki * (Q2 * (S2 + C)) ** self._mm
#                                     if Qs2 < 0:
#                                         Qs2 = 0
#                                     elif Qs2 > fluxLimit:
#                                         Qs2 = fluxLimit
#                                 else:
#                                     Qs2 = 0
#
#                                 if Q3 > self._Qs_min:
#                                     Qs3 = self._Ki * (Q3 * (S3 + C)) ** self._mm
#                                     if Qs3 < 0:
#                                         Qs3 = 0
#                                     elif Qs3 > fluxLimit:
#                                         Qs3 = fluxLimit
#                                 else:
#                                     Qs3 = 0
#
#                             Qs1 = np.nan_to_num(Qs1)
#                             Qs2 = np.nan_to_num(Qs2)
#                             Qs3 = np.nan_to_num(Qs3)
#
#                             # ### Calculate Net Erosion/Accretion
#                             if Elevation[TS, d, i] > self._SL or any(
#                                     z > self._SL
#                                     for z in Elevation[TS, d + 1: d + 10, i]
#                             ):  # If cell is subaerial, elevation change is determined by difference between
#                                 # flux in vs. flux out
#                                 if i > 0:
#                                     SedFluxIn[TS, d + 1, i - 1] += Qs1
#
#                                 SedFluxIn[TS, d + 1, i] += Qs2
#
#                                 if i < (self._BarrierLength - 1):
#                                     SedFluxIn[TS, d + 1, i + 1] += Qs3
#
#                                 Qs_out = Qs1 + Qs2 + Qs3
#                                 SedFluxOut[TS, d, i] = Qs_out
#
#                             else:  # If cell is subaqeous, exponentially decay dep. of remaining sed across bay
#
#                                 if inundation == 0:
#                                     Cbb = self._Cbb_r
#                                 else:
#                                     Cbb = self._Cbb_i
#
#                                 Qs0 = SedFluxIn[TS, d, i] * Cbb
#
#                                 Qs1 = Qs0 * Q1 / (Q1 + Q2 + Q3)
#                                 Qs2 = Qs0 * Q2 / (Q1 + Q2 + Q3)
#                                 Qs3 = Qs0 * Q3 / (Q1 + Q2 + Q3)
#
#                                 Qs1 = np.nan_to_num(Qs1)
#                                 Qs2 = np.nan_to_num(Qs2)
#                                 Qs3 = np.nan_to_num(Qs3)
#
#                                 if Qs1 < self._Qs_bb_min:
#                                     Qs1 = 0
#                                 elif Qs1 > fluxLimit:
#                                     Qs1 = fluxLimit
#                                 if Qs2 < self._Qs_bb_min:
#                                     Qs2 = 0
#                                 elif Qs2 > fluxLimit:
#                                     Qs2 = fluxLimit
#                                 if Qs3 < self._Qs_bb_min:
#                                     Qs3 = 0
#                                 elif Qs3 > fluxLimit:
#                                     Qs3 = fluxLimit
#
#                                 if i > 0:
#                                     SedFluxIn[TS, d + 1, i - 1] += Qs1
#
#                                 SedFluxIn[TS, d + 1, i] += Qs2
#
#                                 if i < (self._BarrierLength - 1):
#                                     SedFluxIn[TS, d + 1, i + 1] += Qs3
#
#                                 Qs_out = Qs1 + Qs2 + Qs3
#                                 SedFluxOut[TS, d, i] = Qs_out
#
#                             # ### Saline Flooding
#                             (
#                                 self._ShrubDomainFemale,
#                                 self._ShrubDomainMale,
#                                 self._ShrubDomainDead,
#                             ) = self.SalineFlooding(
#                                 ShrubDomainWidth,
#                                 self._ShrubDomainAll,
#                                 self._ShrubDomainFemale,
#                                 self._ShrubDomainMale,
#                                 self._ShrubDomainDead,
#                                 d,
#                                 i,
#                                 Q0,
#                             )
#
#                 # ### Update Elevation After Every Storm Hour
#                 if inundation == 1:
#                     ElevationChange = (
#                                               SedFluxIn[TS, :, :] - SedFluxOut[TS, :, :]
#                                       ) / substep
#                 else:
#                     ElevationChange = (
#                                               SedFluxIn[TS, :, :] - SedFluxOut[TS, :, :]
#                                       ) / substep
#                 Elevation[TS, :, :] = Elevation[TS, :, :] + ElevationChange
#
#                 # Calculate and save volume of sediment deposited on/behind the island for every hour
#                 # (all four methods below should equal the same!)
#                 # OWloss = OWloss + np.sum(ElevationChange[1:,:])
#                 # OWloss = OWloss + (np.sum(SedFluxIn[TS,1:,:]) - np.sum(SedFluxOut[TS,1:,:])) / substep
#                 # OWloss = OWloss + np.sum(SedFluxIn[TS,1,:]) / substep
#                 OWloss = OWloss + np.sum(SedFluxOut[TS, 0, :]) / substep
#
#                 # Update amount of burial/erosion for each shrub
#                 self._BurialDomain = self.UpdateBurial(
#                     self._BurialDomain,
#                     ElevationChange,
#                     ShrubDomainWidth,
#                     self._ShrubDomainAll,
#                 )
#
#             # ### Update Interior Domain After Every Storm
#             InteriorUpdate = Elevation[-1, 1:, :]
#
#             # Remove all rows of bay without any deposition from the domain
#             check = 1
#             while check == 1:
#                 if all(x <= -self._BayDepth for x in InteriorUpdate[-1, :]):
#                     InteriorUpdate = np.delete(InteriorUpdate, (-1), axis=0)
#                 else:
#                     check = 0
#
#             # Update interior domain
#             self._InteriorDomain = InteriorUpdate
#
#             # Update Domain widths
#             DomainWidth = np.shape(self._InteriorDomain)[0]
#             (
#                 self._ShrubDomainFemale,
#                 self._ShrubDomainMale,
#                 self._ShrubDomainAll,
#                 self._ShrubPercentCover,
#                 self._BurialDomain,
#                 self._ShrubDomainDead,
#                 self._DeadPercentCover,
#             ) = self.UpdateShrubDomains(
#                 DomainWidth,
#                 ShrubDomainWidth,
#                 self._ShrubDomainFemale,
#                 self._ShrubDomainMale,
#                 self._ShrubDomainAll,
#                 self._ShrubPercentCover,
#                 self._BurialDomain,
#                 self._ShrubDomainDead,
#                 self._DeadPercentCover,
#             )
#
# # Record storm data
# self._StormCount.append(numstorm)
#
# if self._interior_noise_on:  # add noise to flat parts of interior domain
#     if self._time_index > self._StormStart:
#         self._InteriorDomain = self.InteriorNoise(
#             self._InteriorDomain, InteriorWidth
#         )
# plotters
# plt.subplot(1,2,2)
# plt.imshow(
#     cascade.barrier3d[0]._DomainTS[i_min_pt75] * 10,
#     origin=“lower”,
#     cmap=“terrain”,
#     vmin=-1.1,
#     vmax=4.0,
# )
# plt.colorbar()


# --------------------------------------------------------------------------------------------------------------------
# 3. track sediment transport at the dune gaps using Nienhuis?
# only at the dune gap itself and does not alter the gap width

# --------------------------------------------------------------------------------------------------------------------
# 4. add a beach and deposit sediment (in different ways to see what happens?)



# then we will make the class

# class Outwasher:
#     def __init__(self):  # add anything that will be an attribute to the class
#         self.____
#
#     def update(self, b3d):  # use the above functions to get outputs, will likely be using some b3d variables
#        self.___
