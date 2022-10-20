# uses a simpler interior domain
# uses chris' sound levels
# trying to match dune discharges from outwasher new flow routing starting from the first row

import numpy as np
import math
import os
import imageio
from barrier3d import Barrier3d
import matplotlib.pyplot as plt
import imageio
import csv

# -------------------------------------------------------------------------------------------------------------------
# ### Storm Series Converter
# def storm_converter(storm_series):
#     dur = 2 * storm_series[2]
#     sound_level = np.zeros(dur)  # initializing the new bay level array
#     avgs = []
#     # we need to double the duration of the storm with the same magnitudes as before
#     # making the avg water level array
#     count = 0
#     for time in range(len(storm_series[1]) - 1):
#         avg_level = (storm_series[1][time] + storm_series[1][time + 1]) / 2
#         avgs.append(avg_level)
#     # for the new time series, if the time step is even, we keep the value from the original time series
#     # if the time step is odd, we take the value from the avgs array
#     for time2 in range(dur):
#         if time2 % 2 == 0:
#             sound_level[time2] = storm_series[1][int(time2 / 2)]
#         elif time2 == dur - 1:
#             sound_level[time2] = avgs[-1]
#         else:
#             count += 1
#             sound_level[time2] = avgs[time2 - count]
#             # if we are at time 1, need to access avgs time 0 (time - 1)
#             # if we are at time 3, need to access avgs time 1 (time - 2)
#             # if we are at time 5, need to access avgs time 2 (time - 3)
#     storm_series[1] = sound_level
#     return storm_series, dur

# ### Storm Series Converter
def bay_converter(storms, substep):
    # storms will be the sound data, so storm_series[1]
    new_ss = []
    num_new_vals = substep - 1  # number of values we are adding between existing values
    for s in range(len(storms)-1):
        new_ss.append(storms[s])  # add the starting value back in
        inc = (storms[s+1]-storms[s])/substep  # increment for the substeps
        for a in range(num_new_vals):
            new_ss.append(storms[s]+(a+1)*inc)
            # example, we want to do a substep of 3, so we are adding 2 new values
            # if our original storm series is [0, 0.15], then the increment will be 0.05
            # we will add 0 into the new storm series, then 0+1*0.05 = 0.05
            # then 0+2*0.05 = 0.1
            # the next iteration will start with 0.15, and that will be appended as well
    new_ss.append(storms[-1])  # make sure to include our last value
    return new_ss


# ### Calculate Slopes
def calculate_slopes(row, col, domain_width, elev_array, domain_length, time_step, slopes_array, beachface):
    """
    :param d: incremental width (row)
    :param domain_width: width of the input domain
    :param elev_array: array storing the elevations
    :param domain_length: length of the input domain
    :param time_step: time step of the storm
    :param slopes_array: array storing the S2 slopes
    :param beachface: slope of the beachface, used to calc the last row slope, so this must be the beachface slope
    :return: S1, S2, S3, slopes_array
    """
    # ### Calculate Slopes
    # if we are not at the last row, do normal calculations
    if row != domain_width - 1:  # uncomment if we are letting sed out the last row
        # tab everything until else statement
        if col > 0:  # i = 0 means there are no cols to the left
            S1 = (elev_array[time_step, row, col] - elev_array[time_step, row + 1, col - 1]) / (math.sqrt(2))
            S1 = np.nan_to_num(S1)
        else:
            S1 = 0

        S2 = elev_array[time_step, row, col] - elev_array[time_step, row + 1, col]
        S2 = np.nan_to_num(S2)
        # slopes_array[TS, d, i] = S2

        if col < (domain_length - 1):  # i at the end length means there are no cols to the right
            S3 = (elev_array[time_step, row, col] - elev_array[time_step, row + 1, col + 1]) / (math.sqrt(2))
            S3 = np.nan_to_num(S3)
        else:
            S3 = 0
    # if at the last row, apply the same slope as the beachface slope
    else:
        if col > 0:
            S1 = beachface / (math.sqrt(2))
            S1 = np.nan_to_num(S1)
        else:
            S1 = 0

        S2 = beachface
        S2 = np.nan_to_num(S2)

        if col < (domain_length - 1):  # i at the end length means there are no cols to the right
            S3 = beachface / (math.sqrt(2))
            S3 = np.nan_to_num(S3)
        else:
            S3 = 0

    if col == 0 and S2 < 0 and S3 < 0:
        # this is greater than the max slope, so no sediment will go to outside
        S1 = -999
    if col == domain_length - 1 and S2 < 0 and S1 < 0:
        S3 = -999

    slopes_array[time_step, row, col] = S2
    return S1, S2, S3, slopes_array


# ### Calculate Discharges
def calculate_discharges(col, S1, S2, S3, Q0, nn, domain_length, max_slope):
    """
    :param S1: slope to bottom left cell
    :param S2: slope to bottom cell
    :param S3: slope to bottom right cell
    :param Q0: initial discharge
    :param nn: parameter
    :param domain_length: length in the alongshore
    :param max_slope: maximum slope that sediment can be transported up
    :return: Q1, Q2, Q3
    """


    # One or more slopes positive (we have downhill flow)
    if S1 > 0 or S2 > 0 or S3 > 0:

        # flow does not go uphill (when theres a downhill option)
        # and a zero slope yields a zero discharge with this equation

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
    # no downhill slopes, but some that are 0
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

        if S1 == 0 and col > 0:
            Q1 = Qx
        else:
            Q1 = 0
        if S2 == 0:
            Q2 = Qx
        else:
            Q2 = 0
        if S3 == 0 and col < (domain_length - 1):
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
    return Q1, Q2, Q3

class Outwasher:
    """Starts an outwash event

        Examples
        --------
        # >>> from cascade.outwasher_reorganized import Outwasher
        # >>> outwash = Outwasher()
        # >>> outwash.update(barrier3d, storm_series, runID)
        """

    def __init__(
            self,
            b3d,
            runID,
            path,
            substep=2,
            Cx=10,
            Ki=7.5E-3,
    ):
        # make a folder where all graphs will be saved for that run
        self._runID = runID
        self._newpath = path + self._runID + "/"
        if not os.path.exists(self._newpath):
            os.makedirs(self._newpath)

        ### Set Barrier3D variables
        self._berm_el = b3d._BermEl,  # [dam MHW]
        self._beach_elev = b3d._BermEl  # [dam MHW]
        self._length = b3d._BarrierLength  # [dam] length of barrier
        self._substep = substep
        self._nn = b3d._nn  # flow routing constant
        self._max_slope = -b3d._MaxUpSlope  # max slope that sediment can go uphill, previously Slim (0.25)
        # ki = b3d._Ki  # sediment transport coefficient
        self._ki = Ki
        self._cx = Cx
        self._mm = b3d._mm  # inundation overwash coefficient
        # mm = 6
        self._sea_level = b3d._SL  # equal to 0 dam
        self._bay_depth = -b3d._BayDepth  # [dam MHW] Depth of bay behind island segment, currently set to 0.3
        self._Si = (np.mean(b3d.InteriorDomain[-10, :]) - np.mean(b3d.InteriorDomain[0, :])) / len(b3d.InteriorDomain)
        # avg_slope = b3d._BermEl / 20                        # how it is defined in barrier 3D which is much smaller than when
        # you calculate the slope using the avg of the first and last rows
        # setting up dune domain using b3d
        self._dune_domain = b3d.DuneDomain[b3d._time_index - 1, :, :]
        self._dune_crest = self._dune_domain.max(axis=1)  # dune_crest used to be DuneDomainCrest

        # initializing our barrier interior
        # give it an arbitrary width of 30 dam
        self._interior_domain = np.zeros([30, self._length])
        # setting the elevations of each row in the domain
        for row in range(len(self._interior_domain)):
            if row == 0:
                self._interior_domain[row, :] = 0  # the domain (bay) starts at an elevation of 0 as opposed to -0.3 dam in B3D
            else:
                self._interior_domain[row, :] = self._interior_domain[row - 1, :] - self._Si  # giving the back barrier an equal slope
                # (Si is negative, so when we subtract the negative value, we are increasing the elevation)

    def update(
            self,
            storm_series,
            b3d,
            fudge_fac=1
    ):
        ### Set other variables
        q_min = b3d._Qs_min  # [m^3 / hr]? Minimum discharge needed for sediment transport (0.001)
        qs_lost_total = 0  # previously OWloss
        # numstorm = int(len(storm_series))
        numstorm = 1

        if numstorm > 0:
            # ### Individual Storm Impacts
            for n in range(numstorm):
                # setting up the storm series for substep of 2
                if self._substep != 1:
                    updated_bay_levels = bay_converter(storm_series[1], substep=self._substep)
                    dur = len(updated_bay_levels)
                    storm_series[1] = updated_bay_levels
                    storm_series[2] = dur
                    # this has been updated to hopefully be used on any storm series for any substep

                else:
                    dur = storm_series[2]  # [hr] duration of the storm

                # ### Overwash
                # Iow = 0  # Count of dune gaps in inundation regime
                # for q in range(len(gaps)):
                #     Iow += 1  # all inundation regime for us

                # Determine Sediment And Water Routing Rules
                # ### Set Domain
                if n == 0:  # we only need to initialize the domain for the first storm because we want the effects of
                    # storm 1 to stay through the remaining storms
                    # we added a beach (and "beachface") to the domain with a width of 7 dam based on general beach widths
                    beach_domain = np.ones([7, self._length]) * self._beach_elev  # [dam MHW] 7 rows
                    beachface_domain = np.zeros([6, self._length])
                    # we actually want the beach to have a slope, but keep the first few rows the berm elevation
                    # we give the beach slope to be 0.004 m = 0.0004 dam
                    m_beach = 0.0004
                    for b in range(len(beach_domain)):
                        if b >= 3:
                            beach_domain[b, :] = beach_domain[b - 1, :] - m_beach  # m_beach is positive (downhill)
                    self._m_beachface = beach_domain[-1, 0] / len(beachface_domain)  # positive (downhill)
                    for s in range(len(beachface_domain)):
                        if s == 0:
                            beachface_domain[s, :] = beach_domain[-1, 0] - self._m_beachface  # slope of beachface
                        else:
                            beachface_domain[s, :] = beachface_domain[s - 1, :] - self._m_beachface

                    # the dune domain is being taken from B3D, but has 2 rows with length rows, so it needs to be transposed
                    # I believe each row starts off exactly the same, but now I am not sure
                    dune_domain_full = np.transpose(self._dune_domain) + self._berm_el
                    # the full domain of outwasher starts with the interior domain, then the dune, beach, and beachface
                    full_domain = np.append(self._interior_domain, dune_domain_full, 0)  # [dam MHW]
                    full_domain = np.append(full_domain, beach_domain, 0)  # [dam MHW]
                    full_domain = np.append(full_domain, beachface_domain, 0)
                    full_domain = full_domain[20:-1]
                    np.save(self._newpath + "full_domain", full_domain)

                    # plot the initial full domain before sediment movement
                    # fig1 = plt.figure()
                    # ax1 = fig1.add_subplot(111)
                    # mat = ax1.matshow(
                    #     full_domain,
                    #     extent=[0, self._length, 20, np.shape(full_domain)[0]+20],
                    #     # origin="upper",
                    #     cmap="Greens",
                    #     vmin=0, vmax=0.25,
                    # )
                    # fig1.colorbar(mat)
                    # ax1.set_title("Initial Elevation $(dam)$")
                    # ax1.set_ylabel("barrier width (dam)")
                    # ax1.set_xlabel("barrier length (dam)")
                    # plt.gca().xaxis.tick_bottom()
                    # # plt.savefig(self._newpath + "0_domain")
                    # plt.show()

                    # # plotting pre-storm cross section for row 21 (a gap where overwash occurs)
                    # cross_section = np.mean(full_domain, 1)
                    # cross_section = np.flip(cross_section)
                    # fig4 = plt.figure()
                    # ax4 = fig4.add_subplot(111)
                    # ax4.plot(range(len(full_domain)), cross_section, label="pre-storm")
                    # dune_gap_el = np.flip(full_domain[:, 21])
                    # ax4.plot(range(len(full_domain)), dune_gap_el, label="dune gap (21) pre", linestyle="dashed")

                width = np.shape(full_domain)[0]  # width is the number of rows in the full domain
                domain_array = np.zeros([2, width, self._length])
                domain_array[0, :, :] = full_domain
                # duration = dur * substep  # from previous code
                duration = dur  # we already multiplied dur by the substep above
                Elevation = np.zeros([duration, width, self._length])
                # elevation at the first time step is set to the full domain
                Elevation[0, :, :] = full_domain
                OW_TS = []

                # Initialize Memory Storage Arrays
                Discharge = np.zeros([duration, width, self._length])
                SedFluxIn = np.zeros([duration, width, self._length])
                SedFluxOut = np.zeros([duration, width, self._length])
                elev_change_array = np.zeros([duration, width, self._length])
                slopes_array = np.zeros([duration, width, self._length])
                rexcess_dict = {}
                qs2_array = np.zeros([duration, width, self._length])
                qs_lost = 0  # the array for storing the sediment leaving the last row
                # qs_lost is reset every storm, but qs_lost_total is cumulative through all the storms

                # get the average slope of the interior using first and last rows of just the interior domain
                # should be negative because the first row is close to bay and last row close to "dunes"
                # Si = np.mean((interior_domain[-1, :] - interior_domain[0, :]) / 20)

                # plot the bay elevation throughout each storm with sea level and beach elevation references
                # x = range(0, duration)
                # sea_level_line = self._sea_level * np.ones(len(x))
                # beach_elev_line = self._beach_elev * np.ones(len(x))
                # dune_elev_line = max(self._dune_crest + self._berm_el) * np.ones(len(x))

                # fig2 = plt.figure()
                # ax2 = fig2.add_subplot(111)
                # ax2.plot(x, storm_series[1], label='storm {0}'.format(n + 1))
                # # if we have multiple storms, will only need to plot these once
                # ax2.plot(x, sea_level_line, 'k', linestyle='dashed', label='sea level')
                # ax2.plot(x, beach_elev_line, 'k', linestyle='dotted', label='beach elevation')
                # ax2.plot(x, dune_elev_line, 'k', linestyle='dashdot', label='max dune elevation')
                # ax2.set_xlabel("storm duration")
                # ax2.set_ylabel("dam MHW")
                # ax2.set_title("bay elevation over the course of each storm")
                # ax2.legend()
                # plt.show()
                # plt.savefig(self._newpath + "hydrograph")
                # plt.close()

                # ### Run Flow Routing Algorithm
                for TS in range(duration):
                    # Begin with elevation from previous timestep
                    if TS > 0:
                        Elevation[TS, :, :] = Elevation[TS - 1, :, :]
                        # Elevation[TS, 1:, :] = Elevation[TS - 1, 1:, :]  # in B3D, the first row is dunes, so that is
                        # updated from teh dune domain (I think, line 983 in lexi_b3d)
                    print(TS)
                    # get dune crest out here first (dont remember what this means 9/16/2022)
                    bayhigh = storm_series[1][TS]  # [dam]
                    dune_gap = np.min(self._dune_crest + self._berm_el)  # elevation of the lowest dune gap
                    # we only want discharge and sediment transport if the bay level is high enough to move through the dune
                    # gaps. If it is not, nothing happens.
                    if bayhigh <= dune_gap:
                        Discharge[TS, :, :] = 0
                        # nothing happens
                        # sediment suspension stuff?
                    else:
                        # ### DUNES
                        OW_TS.append(TS)
                        int_width = np.shape(self._interior_domain)[0]
                        # max_dune, self._dune_domain, self._dune_crest, gaps, D_not_ow = dune_erosion(b3d,
                        #     self._length, self._berm_el, self._dune_domain, self._dune_crest,
                        #     bayhigh)  # fines the overwashed dune segments, dune

                        # ### finding the OW cells --------------------------------------------------------------------
                        dune_elev = self._dune_crest + self._berm_el
                        Dow = [index for index, value in enumerate(dune_elev) if
                               value < bayhigh]  # bayhigh used to be Rhigh
                        # D_not_ow = [index for index, value in enumerate(dune_elev) if value > bayhigh]
                        gaps = b3d.DuneGaps(self._dune_crest, Dow, self._berm_el, bayhigh)
                        max_dune = b3d._Dmaxel - b3d._BermEl  # [dam MHW]

                        rexcess_array = np.zeros(self._length)
                        for col in range(self._length):
                            rexcess_array[col] = bayhigh - Elevation[TS, int_width-1-20, col]
                        rexcess = np.mean(rexcess_array)
                        overtop_vel = math.sqrt(2 * 9.8 * (rexcess * 10)) / 10  # (dam/s)
                        overtop_flow = overtop_vel * rexcess * 3600  # (dam^3/hr)
                        Discharge[TS, 0, :] = overtop_flow  # (dam^3/hr)

                        # Back barrier flow starts at a specified value to see another value at the dune gaps based on B3D (CHECK)
                        # Loop through the rows
                        # need to include the last row to keep track of how much is leaving the system
                        # for d in range(width-1):
                        # uncomment if we are letting sediment out the last row
                        for d in range(width):  # for letting sediment out, discharge scenarios 2 and 3
                            # for d in range(int_width + 1, width):  # if we are using scenario 1 discharge
                            # Discharge for each TS, row 1, and all cols set above
                            Discharge[TS, d, :][Discharge[TS, d, :] < 0] = 0
                            # Loop through each col of the specified row
                            for i in range(self._length):
                                # if we have discharge, set Qo equal to that value
                                if Discharge[TS, d, i] > 0:
                                    Q0 = Discharge[TS, d, i]  # (dam^3/hr)
                                    # ### Calculate Slopes
                                    S1, S2, S3, slopes_array = calculate_slopes(d, i, width,
                                                                                Elevation, self._length, TS,
                                                                                slopes_array, self._m_beachface)

                                    # ---------for all scenarios calc Qs------------------------------------------------
                                    Q1, Q2, Q3 = calculate_discharges(i, S1, S2, S3, Q0,
                                                                      self._nn, self._length, self._max_slope)

                                    ### Update Discharge
                                    # discharge is defined for the next row, so we do not need to include the last row
                                    # the first row of discharge was already defined
                                    if d != width - 1:  # uncomment and tab until calculate sed movement
                                        # Cell 1
                                        if i > 0:
                                            Discharge[TS, d + 1, i - 1] = Discharge[TS, d + 1, i - 1] + Q1
                                        # Cell 2
                                        Discharge[TS, d + 1, i] = Discharge[TS, d + 1, i] + Q2
                                        # Cell 3
                                        if i < (self._length - 1):
                                            Discharge[TS, d + 1, i + 1] = Discharge[TS, d + 1, i + 1] + Q3
                                    # --------------------------------------------------------------------------------------

                                    # ### Calculate Sed Movement
                                    fluxLimit = max_dune  # [dam MHW] dmaxel - bermel
                                    # all Qs in [dam^3/hr]
                                    C = self._cx * abs(self._Si)  # 10 x the avg slope (from Murray) normal way
                                    # -------------setting different regimes for sed flux momentum (C)----------------------
                                    # if d < int_width:  # in the back barrier with uphill slopes (include dunes)
                                    #     # we would have slower velocities with deeper water = less sed flux
                                    #     depth = bayhigh - Elevation[TS, d, i]  # dimensionless although its tech. dam
                                    #     C = 1/depth*3E-3  # dimensionless although tech. 1/dam
                                    #     # C = Elevation[TS, d, i]
                                    # else:
                                    #     # avg_slope = np.mean(Elevation[TS, int_width+2, :] - Elevation[TS, -1, :]) \
                                    #     #             / (width-int_width)
                                    #     avg_slope = berm_el/(width-int_width)
                                    #     # we could change this to just be the berm elevation - final elevation
                                    #     C = cx * avg_slope
                                    # --------------------------------------------------------------------------------------
                                    if Q1 > q_min:
                                        Qs1 = self._ki * (Q1 * (S1 + C)) ** self._mm
                                        if Qs1 < 0:
                                            Qs1 = 0
                                        elif Qs1 > fluxLimit:
                                            Qs1 = fluxLimit
                                    else:
                                        Qs1 = 0

                                    if Q2 > q_min:
                                        Qs2 = self._ki * (Q2 * (S2 + C)) ** self._mm
                                        if Qs2 < 0:
                                            Qs2 = 0
                                        elif Qs2 > fluxLimit:
                                            Qs2 = fluxLimit
                                    else:
                                        Qs2 = 0

                                    if Q3 > q_min:
                                        Qs3 = self._ki * (Q3 * (S3 + C)) ** self._mm
                                        if Qs3 < 0:
                                            Qs3 = 0
                                        elif Qs3 > fluxLimit:
                                            Qs3 = fluxLimit
                                    else:
                                        Qs3 = 0

                                    multi = 1
                                    Qs1 = np.nan_to_num(Qs1) * multi
                                    Qs2 = np.nan_to_num(Qs2) * multi
                                    Qs3 = np.nan_to_num(Qs3) * multi

                                    qs2_array[TS, d, i] = Qs2

                                    # ### Calculate Net Erosion/Accretion
                                    # flux in vs. flux out
                                    # sed flux in goes to the next row, and is used for determining flux out at current row
                                    # so we need a flux in for the last row, which will be its own variable
                                    if d != width - 1:  # uncomment, tab next two ifs
                                        if i > 0:
                                            SedFluxIn[TS, d + 1, i - 1] += Qs1

                                        SedFluxIn[TS, d + 1, i] += Qs2

                                        if i < (self._length - 1):
                                            SedFluxIn[TS, d + 1, i + 1] += Qs3
                                    # Qs1,2,3 calculated for current row
                                    Qs_out = Qs1 + Qs2 + Qs3
                                    SedFluxOut[TS, d, i] = Qs_out

                                    # END OF DOMAIN LOOPS

                        # print("discharge at dune gaps after flow routing:", Discharge[TS, int_width, :])
                        # ### Update Elevation After Every Storm Hour
                        ElevationChange = (SedFluxIn[TS, :, :] - SedFluxOut[TS, :, :]) / self._substep
                        Elevation[TS, :, :] = Elevation[TS, :, :] + ElevationChange
                        elev_change_array[TS] = ElevationChange

                        # Calculate and save volume of sediment leaving the island for every hour
                        # OWloss = OWloss + np.sum(SedFluxOut[TS, 0, :]) / substep

                        # if we are just letting it go out of the cell uncomment
                        qs_lost = qs_lost + sum(SedFluxOut[TS, width - 1, :]) / self._substep  # [dam^3] previously OWloss
                        # qs_lost is reset to zero within the storm loop
                        qs_lost_total = qs_lost_total + sum(
                            SedFluxOut[TS, width - 1, :]) / self._substep  # [dam^3] previously OWloss
                        # qs_lost_total is initialized to zero outside all loops

                with open('C:/Users/Lexi/Documents/Research/Outwasher/Output/sediment_tracking.txt', 'a') as f:
                    if n == numstorm - 1:
                        f.write(self._runID)
                        f.write('\n')
                        f.write("The total sediment volume lost was {0} dam^3".format(round(qs_lost_total, 3)))
                        f.write('\n')

                print("The sediment volume lost in storm {0} was {1} dam^3".format(n + 1, round(qs_lost, 3)))
                if n == numstorm - 1:
                    print("The total sediment volume lost was {0} dam^3".format(round(qs_lost_total, 3)))

                # ### Update Interior Domain After Every Storm
                # use the last TS
                InteriorUpdate = Elevation[-1, :, :]

                # Remove all rows of bay without any deposition from the domain
                check = 1
                while check == 1:
                    if all(x <= self._bay_depth for x in InteriorUpdate[-1, :]):
                        InteriorUpdate = np.delete(InteriorUpdate, (-1), axis=0)
                    else:
                        check = 0

                # Update interior domain
                # int_update_b3d = InteriorUpdate[0:int_width, :]
                # b3d._InteriorDomain = np.flip(int_update_b3d)
                # full_domain[:, :] = InteriorUpdate
                domain_array[1, :, :] = InteriorUpdate
                # Update Domain widths
                # DomainWidth = np.shape(b3d._InteriorDomain)[0]

                # # plotting post-storm cross section
                # cross_section2 = np.mean(full_domain, 1)
                # cross_section2 = np.flip(cross_section2)
                # ax4.plot(range(len(full_domain)), cross_section2, label="post storm", linestyle="dashed")
                # # dune_gap_el = np.flip(full_domain[:, 21])
                # # ax4.plot(range(len(full_domain)), dune_gap_el, label="dune gap (21) post", linestyle="dashed")
                # ax4.set_xlabel("barrier width from ocean to bay (dam)")
                # ax4.set_ylabel("average alongshore elevation (dam)")
                # # ax4.set_title("Cross shore elevation profile for col 21\n "
                # #               "beachface slope = {0}, back barrier slope = {1}".format(round(m_beachface, 3), round(Si, 3)))
                # ax4.set_title("Average cross shore elevation profile\n "
                #               "beachface slope = {0}, back barrier slope = {1}".format(round(m_beachface, 3),
                #                                                                        round(self._Si, 3)))
                # ax4.legend()
                # plt.show()
                # # plt.savefig(newpath + "cross_shore_21")
                # plt.savefig(self._newpath + "avg_cross_shore")

                # # plot post storm elevation
                # fig3 = plt.figure()
                # ax3 = fig3.add_subplot(111)
                # mat2 = ax3.matshow(
                #     full_domain[:, :],
                #     # origin="upper",
                #     cmap="Greens",
                #     # vmin=-0.2, vmax=0.4,
                # )
                # ax3.set_xlabel('barrier length (dam)')
                # ax3.set_ylabel('barrier width (dam)')
                # ax3.set_title("Elevation after storm {0} $(dam)$".format(n + 1))
                # plt.gca().xaxis.tick_bottom()
                # fig3.colorbar(mat2)
                # plt.savefig(self._newpath + "{0}_domain".format(n + 1))

        # Record storm data
        b3d._StormCount.append(numstorm)
        return Discharge, elev_change_array, full_domain, qs_lost_total, slopes_array, rexcess_dict, qs2_array, \
               storm_series, SedFluxOut, SedFluxIn, domain_array, OW_TS


# --------------------------------------------running outwasher---------------------------------------------------------
# importing Chris' bay data

# ### start of the actual code

# with open(r"C:\Users\Lexi\Documents\Research\Outwasher\chris stuff\sound_data.txt", newline='') as csvfile:
#     sound_data = list(csv.reader(csvfile))[0]
# sound_data = [float(s) / 10 - 0.054 for s in sound_data]  # [dam MHW] Chris' sound elevations were in m MSL,
# # so converted to NAVD88 then MHW and dam
# sound_data = [s + 0.05 for s in sound_data]  # [dam MHW] just increasing the values
# # setting all negative values to 0
# sound_data = sound_data[20:]
# for index, value in enumerate(sound_data):
#     #     # smaller used 0.05
#     if value > 0.220:
#         sound_data[index] = 0.220
# sound_data[0] = 0
# # np.save("C:/Users/Lexi/Documents/Research/Outwasher/sound_data", sound_data)
#
# # storm series is year the storm occured, the bay elevation for every time step, and the duration of the storm
# storm_series = [1, sound_data, len(sound_data)]
# path = "C:/Users/Lexi/Documents/Research/Outwasher/Output/edgesedited_bay220limited/"
# # runID = "test"
# runID = "dynamic_discharge_AVG_FACTOR15_Kie-3_substep4"
# # the number in runID is 0.__
# # ss in runID stands for storm series
# # syndunes = synthetic dunes
# # sedout = sediment fully leaves the system offshore
# # see flux notes for changes to sediment fluxes (around lines 555 in code)
# ## NOW USING REGULAR SOUND DATA, SED OUT AND EDITED EDGES
# b3d = Barrier3d.from_yaml("C:/Users/Lexi/PycharmProjects/Barrier3d/tests/test_params/")
# b3d.update()
# b3d.update_dune_domain()
# # from cascade.outwasher_reorganized import Outwasher
# outwash = Outwasher(b3d, runID, path, substep=4, Cx=10, Ki=7.5E-3,)
# discharge, elev_change, domain, qs_lost, slopes2, dictionary, qs2, avg_initial_cross, storm_elev, sedout, sedin, domain_array \
#     = outwash.update(storm_series)

# domain_change = domain_array[1] - domain_array[0]
# fig5 = plt.figure()
# ax5 = fig5.add_subplot(111)
# mat5 = ax5.matshow(
#     domain_change,
#     # origin="upper",
#     cmap="seismic",
#     vmin=-0.2, vmax=0.2,
# )
# ax5.set_xlabel('barrier length (dam)')
# ax5.set_ylabel('barrier width (dam)')
# ax5.set_title("Elevation Change")
# plt.gca().xaxis.tick_bottom()
# fig5.colorbar(mat5)
# plt.savefig("C:/Users/Lexi/Documents/Research/Outwasher/Output/edgesedited_bay220limited/" + runID + "/elev_change_domain")


# fig5 = plt.figure()
# ax5 = fig5.add_subplot(111)
# cols = range(np.size(qs2, 1))
# for col in cols:
#     line = qs2[7, :, col]
#     ax5.plot(cols, line, label="column {0}".format(col))
# ax5.legend()
# ax5.set_ylabel("Qs2 (dam3/hr)")
# ax5.set_xlabel("Cross-shore Distance from Bay to Ocean (dam)")
# ax5.set_title("{0} \n Qs2 at time 7".format(runID))
# plt.savefig("C:/Users/Lexi/Documents/Research/Outwasher/Output/" + runID + "/cross_shore_qs2".format(runID))

# fig6 = plt.figure()
# ax6 = fig6.add_subplot(111)
# ax6.plot(avg_initial_cross)
# x = len(avg_initial_cross)
# for index, value in enumerate(storm_elev):
#     if index < 12:
#         ax6.plot(range(x), np.ones(x) * value, linestyle="dashed", label="TS = {0}".format(index))
# ax6.set_title("sound level for first 12 TS")
# ax6.legend()
# ax6.set_ylabel("Elevation (dam)")
# ax6.set_xlabel("Cross-shore Distance from Bay to Ocean (dam)")
# plt.savefig("C:/Users/Lexi/Documents/Research/Outwasher/Output/edgesedited_bay220limited/" + runID + "/baylevels")
# np.save("C:/Users/Lexi/Documents/Research/Outwasher/discharge", discharge)
# np.save(newpath + "discharge", discharge)
# plt.matshow(slopes2[1], cmap="jet_r")
# plt.title('S2 Slopes')
# plt.colorbar()
# plt.matshow(discharge[1], cmap="jet_r")
# plt.title('Discharge (dam^3/hr)')
# plt.colorbar()
# ----------------------------------making the elevation gif------------------------------------------------------------
# frames = []
# for i in range(2):
#     filename = "C:/Users/Lexi/Documents/Research/Outwasher/Output/edgesedited_bay220limited/" + runID + "/" + str(i) + "_domain.png"
#     frames.append(imageio.imread(filename))
# imageio.mimwrite("C:/Users/Lexi/Documents/Research/Outwasher/Output/edgesedited_bay220limited/{0}/test.gif".format(runID), frames, format='.gif',
#                  fps=1)


# -------------------------------------------elevation gif--------------------------------------------------------------
def plot_ElevAnimation(elev, directory, start, stop):
    os.chdir(directory)
    newpath = "Elevations/"
    if not os.path.exists(newpath):
        os.makedirs(newpath)
    os.chdir(newpath)

    for t in range(start, stop):
        AnimateDomain = elev[t]

        # Plot and save
        elevFig1 = plt.figure(figsize=(15, 7))
        ax = elevFig1.add_subplot(111)
        cax = ax.matshow(
            AnimateDomain,
            # origin="upper",
            cmap="seismic",
            # vmin=-0.000002, vmax=0.000002
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
        elevFig1.savefig(name, facecolor='w')  # dpi=200
        plt.close(elevFig1)

    frames = []

    for filenum in range(start, stop):
        filename = "elev_" + str(filenum) + ".png"
        frames.append(imageio.imread(filename))
    imageio.mimsave("elev.gif", frames, fps=2)
    # imageio.mimsave("elev.gif", frames, "GIF-FI")
    print("[ * elevation GIF successfully generated * ]")


# -------------------------------------------discharge gif--------------------------------------------------------------
def plot_DischargeAnimation(dis, directory, start, stop):
    os.chdir(directory)
    newpath = "Discharges/"
    if not os.path.exists(newpath):
        os.makedirs(newpath)
    os.chdir(newpath)

    # for t in range(TMAX - 1):
    for t in range(start, stop):
        AnimateDomain = dis[t]

        # Plot and save
        elevFig1 = plt.figure(figsize=(15, 7))
        ax = elevFig1.add_subplot(111)
        cax = ax.matshow(
            AnimateDomain,
            # origin="upper",
            cmap="jet_r",
            # vmin=0, vmax=20,
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
        elevFig1.savefig(name, facecolor='w')  # dpi=200
        plt.close(elevFig1)

    frames = []

    for filenum in range(start, stop):
        filename = "dis_" + str(filenum) + ".png"
        frames.append(imageio.imread(filename))
    imageio.mimsave("dis.gif", frames, fps=2)
    # imageio.mimsave("dis.gif", frames, "GIF-FI")
    print()
    print("[ * discharge GIF successfully generated * ]")


# ---------------------------------------------------slope gif----------------------------------------------------------
def plot_SlopeAnimation(slope, directory, start, stop):
    os.chdir(directory)
    newpath = "Slopes/"
    if not os.path.exists(newpath):
        os.makedirs(newpath)
    os.chdir(newpath)

    for t in range(start, stop):
        AnimateDomain = slope[t]

        # Plot and save
        elevFig1 = plt.figure(figsize=(15, 7))
        ax = elevFig1.add_subplot(111)
        cax = ax.matshow(
            AnimateDomain,
            # origin="upper",
            # cmap="jet_r",
            # vmin=-0.1, vmax=0.1,
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
        elevFig1.savefig(name, facecolor='w')  # dpi=200
        plt.close(elevFig1)

    frames = []

    for filenum in range(start, stop):
        filename = "slope_" + str(filenum) + ".png"
        frames.append(imageio.imread(filename))
    imageio.mimsave("dis.gif", frames, fps=2)
    # imageio.mimsave("slopes.gif", frames, "GIF-FI")
    print()
    print("[ * slope GIF successfully generated * ]")


# -------------------------------------------qs2 gif--------------------------------------------------------------------
def plot_Qs2Animation(qs2, directory, start, stop):
    os.chdir(directory)
    newpath = "Qs2/"
    if not os.path.exists(newpath):
        os.makedirs(newpath)
    os.chdir(newpath)

    # for t in range(TMAX - 1):
    for t in range(start, stop):
        AnimateDomain = qs2[t]

        # Plot and save
        elevFig1 = plt.figure(figsize=(15, 7))
        ax = elevFig1.add_subplot(111)
        cax = ax.matshow(
            AnimateDomain,
            # origin="upper",
            cmap="jet_r",
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
        elevFig1.savefig(name, facecolor='w')  # dpi=200
        plt.close(elevFig1)

    frames = []

    for filenum in range(start, stop):
        filename = "qs2_" + str(filenum) + ".png"
        frames.append(imageio.imread(filename))
    imageio.mimsave("dis.gif", frames, fps=2)
    # imageio.mimsave("qs2.gif", frames, "GIF-FI")
    print()
    print("[ * Qs2 GIF successfully generated * ]")


# ---------------------- Sed out array ----------------------------------------------------------------------------------
def plot_SedOutAnimation(sedout, directory, start, stop):
    os.chdir(directory)
    newpath = "SedOut/"
    if not os.path.exists(newpath):
        os.makedirs(newpath)
    os.chdir(newpath)

    # for t in range(TMAX - 1):
    for t in range(start, stop):
        AnimateDomain = sedout[t]

        # Plot and save
        elevFig1 = plt.figure(figsize=(15, 7))
        ax = elevFig1.add_subplot(111)
        cax = ax.matshow(
            AnimateDomain,
            # origin="upper",
            cmap="jet_r",
            # vmin=min_v, vmax=max_v,
        )  # , interpolation='gaussian') # analysis:ignore
        ax.xaxis.set_ticks_position("bottom")
        elevFig1.colorbar(cax)
        plt.xlabel("Alongshore Distance (dam)")
        plt.ylabel("Cross-Shore Distance (dam)")
        plt.title("Sed Flux Out (dam^3)")
        plt.tight_layout()
        timestr = "Time = " + str(t) + " hrs"
        plt.text(1, 1, timestr)
        plt.rcParams.update({"font.size": 15})
        name = "sedout_" + str(t)
        elevFig1.savefig(name, facecolor='w')  # dpi=200
        plt.close(elevFig1)

    frames = []

    for filenum in range(start, stop):
        filename = "sedout_" + str(filenum) + ".png"
        frames.append(imageio.imread(filename))
    imageio.mimsave("dis.gif", frames, fps=2)
    # imageio.mimsave("sedout.gif", frames, "GIF-FI")
    print()
    print("[ * SedOut GIF successfully generated * ]")


# ---------------------- Sed out array ----------------------------------------------------------------------------------
def plot_SedInAnimation(sedin, directory, start, stop):
    os.chdir(directory)
    newpath = "SedIn/"
    if not os.path.exists(newpath):
        os.makedirs(newpath)
    os.chdir(newpath)

    # for t in range(TMAX - 1):
    for t in range(start, stop):
        AnimateDomain = sedin[t]

        # Plot and save
        elevFig1 = plt.figure(figsize=(15, 7))
        ax = elevFig1.add_subplot(111)
        cax = ax.matshow(
            AnimateDomain,
            # origin="upper",
            cmap="jet_r",
            # vmin=min_v, vmax=max_v,
        )  # , interpolation='gaussian') # analysis:ignore
        ax.xaxis.set_ticks_position("bottom")
        elevFig1.colorbar(cax)
        plt.xlabel("Alongshore Distance (dam)")
        plt.ylabel("Cross-Shore Distance (dam)")
        plt.title("Sed Flux In (dam^3)")
        plt.tight_layout()
        timestr = "Time = " + str(t) + " hrs"
        plt.text(1, 1, timestr)
        plt.rcParams.update({"font.size": 15})
        name = "sedin_" + str(t)
        elevFig1.savefig(name, facecolor='w')  # dpi=200
        plt.close(elevFig1)

    frames = []

    for filenum in range(start, stop):
        filename = "sedin_" + str(filenum) + ".png"
        frames.append(imageio.imread(filename))
    imageio.mimsave("dis.gif", frames, fps=2)
    # imageio.mimsave("sedin.gif", frames, "GIF-FI")
    print()
    print("[ * SedIn GIF successfully generated * ]")


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
        dunes_y = (b3d._DuneDomain[t, v, :] + b3d._BermEl)  # dam MHW
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


# # change if substep is 2
# # TMAX = storm_series[2]
# # TMAX = 2*storm_series[2]
# TMAX = 20
# dir = "C:/Users/Lexi/Documents/Research/Outwasher/Output/edgesedited_bay220limited/" + runID + "/"
# plot_ElevAnimation(elev_change, dir, TMAX)
# plot_DischargeAnimation(discharge, dir, TMAX)
# plot_SlopeAnimation(slopes2, dir, TMAX)
# # plot_Qs2Animation(qs2, dir, TMAX)
# plot_SedOutAnimation(sedout, dir, TMAX)
# plot_SedInAnimation(sedin, dir, TMAX)
# # time_step = [0]
# # plot_ModelTransects(b3d, time_step)
