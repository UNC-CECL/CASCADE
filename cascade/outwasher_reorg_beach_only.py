# uses a simpler interior domain
# uses chris' sound levels
# trying to match dune discharges from outwasher new flow routing starting from the first row

import numpy as np
import math
import os
import matplotlib.pyplot as plt
import imageio

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

# ### Dune Erosion
def dune_erosion(b3d, dune_length, berm_el, dune_domain, dune_crest, bayhigh):
    dune_width = b3d._DuneWidth
    dune_restart = b3d._DuneRestart  # currently set to 0.0075
    max_dune = b3d._Dmaxel - b3d._BermEl  # [dam MHW]
    Hd_avgTS = b3d._Hd_AverageTS
    dune_crest[dune_crest < dune_restart] = dune_restart
    Hd_avgTS.append(np.mean(dune_crest))  # Store average pre-storm dune-height for time step
    c1 = b3d._C1
    c2 = b3d._C2
    # Hd_loss_TS = b3d._Hd_Loss_TS
    # Rhigh = max(storm_series[1])  # Highest elevation of the landward margin of runup. Just using max of storm series
    DuneChange = np.zeros([dune_length, dune_width])  # Vector storing dune height change for this storm
    # dune = dune_crest+berm_el

    # Find overwashed dunes and gaps
    # currently changed Rhigh to current bay depth in the storm series
    dune_elev = dune_crest + berm_el
    Dow = [index for index, value in enumerate(dune_elev) if value < bayhigh]  # bayhigh used to be Rhigh
    D_not_ow = [index for index, value in enumerate(dune_elev) if value > bayhigh]
    gaps = b3d.DuneGaps(dune_crest, Dow, berm_el, bayhigh)
    # Finds location and Rexcess of continuous gaps in dune ridge
    for ow_cell in range(len(Dow)):  # Loop through each overwashed dune cell
        for w in range(dune_width):
            # Calculate dune elevation loss
            ## dune domain at time index 1 is all zeros
            Rnorm = bayhigh / (dune_domain[Dow[ow_cell], w] + berm_el)  # bayhigh relative to pre-storm dune elevation
            Dloss = Rnorm / (c1 + (Rnorm * (Rnorm - c2)))  # Amount of dune crest elevation change normalized by
            # pre-storm dune elevation (increased by a factor by LVB) not sure this is still true
            # (i.e. a percent change), from Goldstein and Moore (2016)
            # Set new dune height
            InitDElev = (dune_domain[Dow[ow_cell], w] + berm_el)
            NewDElev = InitDElev * (1 - Dloss)  # Calculate new dune elevation from storm lowering
            if NewDElev < berm_el:
                NewDElev = berm_el
            dune_domain[Dow[ow_cell], w] = (NewDElev[0] - berm_el[0])  # Convert elevation to height above berm
            # LVB got rid of parenthesis around NewDElev - berm_el
            row = Dow[ow_cell]
            DuneChange[Dow[ow_cell], w] = InitDElev - NewDElev

            # If dune is lowered to ~ zero, allow for chance of regrowth by raising dune height to 5 cm
            if dune_domain[Dow[ow_cell], w] < dune_restart:
                if dune_restart < max_dune:
                    dune_domain[Dow[ow_cell], w] = dune_restart
                else:
                    dune_domain[Dow[ow_cell], w] = (max_dune)  # Restart height can't be greater than Dmax

    # Dune Height Diffusion
    # dune_domain = b3d.DiffuseDunes(dune_domain, 0)
    # dune_domain[dune_domain < sea_level] = dune_restart
    # # DuneLoss = np.sum(DuneChange) / length  # I think duneloss is only used for shoreline change
    # # Hd_TSloss = (DuneChange.max(axis=1) / dur)  # Average height of dune loss for each substep during storm
    # Hd_TSloss = (DuneChange.max(axis=1) / 1)  # Average height of dune loss for each substep during storm
    # # I think we should change dur to 1 because this will change with each TS now
    # # it was previously calculated one time and then each TS would have the same value
    # # this is usually used to reduce dune height lineraly over the course of the storm
    # # but I think it was super high, so I stopped using it
    # Hd_loss_TS[0, :] = Hd_loss_TS[0, :] + DuneChange.max(axis=1)
    return max_dune, dune_domain, dune_crest, gaps, D_not_ow


# ### Calculate Slopes
def calculate_slopes(row, col, domain_width, elev_array, domain_length, time_step, slopes_array, beachface):
    """
    :param d: incremental width (row)
    :param domain_width: width of the input domain
    :param elev_array: array storing the elevations
    :param domain_length: length of the input domain
    :param time_step: time step of the storm
    :param slopes_array: array storing the S2 slopes
    :param beachface: slope of the beachface
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
        self._ki = Ki
        self._cx = Cx
        self._mm = b3d._mm  # inundation overwash coefficient
        # mm = 6
        self._sea_level = b3d._SL  # equal to 0 dam
        self._bay_depth = -b3d._BayDepth  # [dam MHW] Depth of bay behind island segment, currently set to 0.3
        self._Si = (np.mean(b3d.InteriorDomain[-10, :]) - np.mean(b3d.InteriorDomain[0, :])) / len(b3d.InteriorDomain)
        self._dune_domain = b3d.DuneDomain[b3d._time_index - 1, :, :]
        self._dune_crest = self._dune_domain.max(axis=1)  # dune_crest used to be DuneDomainCrest

        # # initializing our barrier interior
        # # give it an arbitrary width of 30 dam
        # self._interior_domain = np.zeros([30, self._length])
        # # setting the elevations of each row in the domain
        # for row in range(len(self._interior_domain)):
        #     if row == 0:
        #         self._interior_domain[row, :] = 0  # the domain (bay) starts at an elevation of 0 as opposed to -0.3 dam in B3D
        #     else:
        #         self._interior_domain[row, :] = self._interior_domain[row - 1, :] - self._Si  # giving the back barrier an equal slope
        #         # (Si is negative, so when we subtract the negative value, we are increasing the elevation)

    def update(
            self,
            storm_series,
            b3d,
            m_beach
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
                    # beach_domain = np.ones([7, self._length]) * self._beach_elev  # [dam MHW] 7 rows
                    # beachface_domain = np.zeros([6, self._length])
                    # # we actually want the beach to have a slope, but keep the first few rows the berm elevation
                    # # we give the beach slope to be 0.004 m = 0.0004 dam
                    # m_beach = 0.0004
                    # for b in range(len(beach_domain)):
                    #     if b >= 3:
                    #         beach_domain[b, :] = beach_domain[b - 1, :] - m_beach  # m_beach is positive (downhill)
                    #     beach_domain[b, :] = beach_domain[b - 1, :] - m_beach  # m_beach is positive (downhill)
                    # self._m_beachface = beach_domain[-1, 0] / len(beachface_domain)  # positive (downhill)
                    # for s in range(len(beachface_domain)):
                    #     if s == 0:
                    #         beachface_domain[s, :] = beach_domain[-1, 0] - self._m_beachface  # slope of beachface
                    #     else:
                    #         beachface_domain[s, :] = beachface_domain[s - 1, :] - self._m_beachface
                    # ### ------------edited for just a beach with same slope as beachface------------------------------
                    beach_domain = np.ones([12, self._length]) * self._beach_elev  # [dam MHW] 7 rows
                    for b in range(len(beach_domain)):
                        beach_domain[b, :] = beach_domain[b - 1, :] - m_beach  # m_beach is positive (downhill)

                    # the dune domain is being taken from B3D, but has 2 rows with length rows, so it needs to be transposed
                    # I believe each row starts off exactly the same, but now I am not sure
                    dune_domain_full = np.transpose(self._dune_domain) + self._berm_el
                    # the full domain of outwasher starts with the interior domain, then the dune, beach, and beachface
                    full_domain = beach_domain
                    # full_domain = np.append(self._interior_domain, dune_domain_full, 0)  # [dam MHW]
                    # full_domain = np.append(full_domain, beach_domain, 0)  # [dam MHW]
                    # full_domain = np.append(full_domain, beachface_domain, 0)
                    # full_domain[0, :, :] = np.vstack([interior_domain, dune_domain_full, beach_domain, beachface_domain])
                    np.save(self._newpath + "full_domain", full_domain)

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
                    dune_gap = np.min(self._dune_crest + self._berm_el)  # elevation of the dune gap
                    # we only want discharge and sediment transport if the bay level is high enough to move through the dune
                    # gaps. If it is not, nothing happens.
                    if bayhigh <= dune_gap:
                        Discharge[TS, :, :] = 0
                        # nothing happens
                        # sediment suspension stuff?
                    else:
                        # ### DUNES
                        OW_TS.append(TS)
                        # int_width = np.shape(self._interior_domain)[0]
                        max_dune, self._dune_domain, self._dune_crest, gaps, D_not_ow = dune_erosion(b3d,
                            self._length, self._berm_el, self._dune_domain, self._dune_crest,
                            bayhigh)  # fines the overwashed dune segments, dune
                        # gaps, and lowers dunes (gaps only?) based on the excess water height
                        # returns the new dune domain values, so we need to update that in the full domain
                        # dune_domain_full = np.transpose(self._dune_domain) + self._berm_el
                        # full_domain[int_width, :] = dune_domain_full[0, :]
                        # full_domain[int_width + 1] = dune_domain_full[1, :]
                        # Elevation[TS, int_width, :] = dune_domain_full[0, :]
                        # Elevation[TS, int_width + 1, :] = dune_domain_full[1, :]

                        dunes_prestorm = self._dune_crest
                        dunes = dunes_prestorm + self._berm_el  # can be used later for reducing dune height with storm
                        # Elevation[TS, int_width:(int_width+2), :] = dunes - \
                        #                               (Hd_TSloss / substep * TS)
                        # Reduce dune in height linearly over course of storm

                        # -------------------------using Rexcess to set discharge levels--------------------------------
                        # loops through each of the gaps and get their starts/stop index, as well as the Rexcess
                        rexcess_tot = 0
                        for q in range(len(gaps)):
                            start = gaps[q][0]
                            stop = gaps[q][1]
                            Rexcess = gaps[q][2]  # (dam)
                            rexcess_tot += Rexcess
                            # Calculate discharge through each dune cell
                            overtop_vel = math.sqrt(2 * 9.8 * (Rexcess * 10)) / 10  # (dam/s)
                            factor = 15
                            overtop_flow = overtop_vel * Rexcess * 3600 * factor  # (dam^3/hr)
                            # Set discharge at dune gap
                            Discharge[:, 0, start:stop] = overtop_flow  # (dam^3/hr) in B3D, no dunes so start at row 0
                            # in Outwasher, dunes located in the 2 rows after interior domain
                            # (cells 0-29, but int_width = 30)
                            # Discharge[TS, int_width, start:stop] = overtop_flow  # (dam^3/hr)
                            # Discharge[TS, int_width + 1, start:stop] = overtop_flow  # (dam^3/hr)
                        # ----------------------------------------------------------------------------------------------
                        # # 1. method of setting discharge through the interior: conservation
                        # # we evenly distribute over each cell the total discharge through the dune gaps
                        # avg_discharge = sum(Discharge[TS, int_width, :]) / self._length
                        # Discharge[TS, 0:int_width, :] = avg_discharge  # (dam^3/hr)
                        # # ----------------------------------------------------------------------------------------------
                        # # 2. method of setting discharge through the interior: fudge factor
                        # # we take the known discharge through the dune gaps, and multiply it by a factor for initial
                        # # discharge at the first row
                        # # we will then ultimately overwrite the dune gap values
                        # # factor = 10
                        # # avg_rexcess = rexcess_tot/len(gaps)
                        # # overtop_vel = math.sqrt(2 * 9.8 * (avg_rexcess * 10)) / 10  # (dam/s)
                        # # overtop_flow = overtop_vel * avg_rexcess * 3600  # (dam^3/hr)
                        # # Discharge[TS, 0, :] = overtop_flow*factor  # (dam^3/hr)
                        # # --------------------------------------------------------------------------------------------------
                        # # 3. just setting the discharge arbitrarily at the first row
                        # # Discharge[TS, 0, :] = 150

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
                                                                                slopes_array, m_beach)
                                    # ### Calculate Discharge To Downflow Neighbors
                                    # do we want water only routed through dune gaps?
                                    # Discharge[TS, int_width, D_not_ow] = 0  # (dam^3/hr)
                                    # Discharge[TS, int_width + 1, D_not_ow] = 0  # (dam^3/hr)
                                    # if Discharge[TS, d, i] == 0:
                                    #     Q0 = 0

                                    # --------------------------------------------------------------------------------------
                                    Q1, Q2, Q3 = calculate_discharges(i, S1, S2, S3, Q0,
                                                                      self._nn, self._length, self._max_slope)
                                    # if we are using scenario 1, we do not want to re-assign discharges for the back barrier
                                    # if d > int_width:
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
                                    # For scenarios 2 and 3, we want to calc Q throughout entire domain
                                    # ### Update Discharge
                                    # # discharge is defined for the next row, so we do not need to include the last row
                                    # # the first row of discharge was already defined
                                    # if d != width - 1:  # uncomment and tab until calculate sed movement
                                    #     # Cell 1
                                    #     if i > 0:
                                    #         Discharge[TS, d + 1, i - 1] = Discharge[TS, d + 1, i - 1] + Q1
                                    #     # Cell 2
                                    #     Discharge[TS, d + 1, i] = Discharge[TS, d + 1, i] + Q2
                                    #     # Cell 3
                                    #     if i < (self._length - 1):
                                    #         Discharge[TS, d + 1, i + 1] = Discharge[TS, d + 1, i + 1] + Q3
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
                full_domain[:, :] = InteriorUpdate
                domain_array[1, :, :] = full_domain
                # Update Domain widths
                # DomainWidth = np.shape(b3d._InteriorDomain)[0]

        # Record storm data
        b3d._StormCount.append(numstorm)
        return Discharge, elev_change_array, full_domain, qs_lost_total, slopes_array, rexcess_dict, qs2_array, \
               storm_series, SedFluxOut, SedFluxIn, domain_array, OW_TS


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
        elevFig1.savefig(name)  # dpi=200
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
        elevFig1.savefig(name)  # dpi=200
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
        elevFig1.savefig(name)  # dpi=200
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
        elevFig1.savefig(name)  # dpi=200
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
        elevFig1.savefig(name)  # dpi=200
        plt.close(elevFig1)

    frames = []

    for filenum in range(start, stop):
        filename = "sedout_" + str(filenum) + ".png"
        frames.append(imageio.imread(filename))
    imageio.mimsave("dis.gif", frames, fps=2)
    # imageio.mimsave("sedout.gif", frames, "GIF-FI")
    print()
    print("[ * SedOut GIF successfully generated * ]")


# ---------------------- Sed in array ----------------------------------------------------------------------------------
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
        elevFig1.savefig(name)  # dpi=200
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

