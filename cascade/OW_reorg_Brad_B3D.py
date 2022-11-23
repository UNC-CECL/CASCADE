# moving boundary conditions for flow routing
# uphill = no flow routing
# downhill = flow routing

import numpy as np
import math
import os
import imageio
import matplotlib.pyplot as plt
import imageio
import copy


def bay_converter(storms, substep):
    """ takes a hydrograph and linearly interpolates a specified number of additional points based on the substep value
    :param storms: list or array of bay level
    :param substep: number of extra points you want
    :returns new list of bay levels
    (ex, substep 2 = flow route every 30 minutes by interpolating one extra value)
    """
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


def dunes(length, berm_el, n_rows, n_gaps, dune_height=0.25):
    """
    creates a dune line with gaps
    :param length: alongshore length of the barrier, int
    :param berm_el: height of the berm [dam MHW], float
    :param n_rows: number of dune rows (width of dunes), int
    :param n_gaps: number of dune gaps, int
    :param dune_height: desired height of dune gaps, float
    :return: dune domain, numpy array
    """
    width = n_rows
    gap_height = berm_el[0]
    dune_domain = np.zeros([width, length])
    dune_domain[0, :] = dune_height
    # centered evenly-spaced gaps
    gap = 1
    gap1 = 1 / n_gaps * 0.5 * length
    gap_locations = np.zeros(n_gaps)
    gap_locations[0] = gap1
    while gap <= n_gaps-1:
        gap_locations[gap] = gap_locations[gap-1]+1/n_gaps*length
        gap += 1
    gap_locations = np.round(gap_locations, decimals=1)
    gap_locations = gap_locations.astype(int)
    # random gaps
    # for g in range(n_gaps):
    #     gap_locations[g] = np.random.randint(0, length-1)
    # gap_locations = gap_locations.astype(int)
    # for both methods:
    for loc in gap_locations:
        dune_domain[0, loc] = gap_height
        dune_domain[0, loc-1] = gap_height
        dune_domain[0, loc+1] = gap_height
    dune_domain[:, :] = dune_domain[0, :]
    return dune_domain


def DuneGaps(DuneDomain, Dow, Rhigh):
    """Returns tuple of [gap start index, stop index, avg Rexcess of each gap,
    alpha: ratio of TWL / dune height]
    :param DuneDomain: numpy array with dune elevations [dam MHW]
    :param Dow: list of cell indices where outwash occurs
    :param Rhigh: current bay level [dam MHW]
    """
    gaps = []
    start = 0
    #        stop = 0
    i = start
    while i < (len(Dow) - 1):
        adjacent = Dow[i + 1] - Dow[i]
        if adjacent == 1:
            i = i + 1
        else:
            stop = i
            x = DuneDomain[Dow[start]: (Dow[stop] + 1)]
            Hmean = sum(x) / float(len(x))
            Rexcess = Rhigh - Hmean
            # alpha = Rhigh / (Hmean + bermel)
            gaps.append([Dow[start], Dow[stop], Rexcess])

            start = stop + 1
            i = start
    # if i > 0:
        # stop = i - 1
    stop = i
    x = DuneDomain[Dow[start]: (Dow[stop] + 1)]
    if len(x) > 0:
        Hmean = sum(x) / float(len(x))
        Rexcess = Rhigh - Hmean
        # alpha = Rhigh / (Hmean + bermel)
        gaps.append([Dow[start], Dow[stop], Rexcess])
    return gaps


def FR_slopes(truth_array, avg_slope_array, domain, width, length, time_step):
    """
    takes the elevations and differentiates uphill and downhill regions based on an average slope of
    a block of cells
    :param truth_array: empty array that is filled with 1s and 0s based on uphill versus downhill slopes
    :param avg_slope_array: currently empty, for each cell, stores the average value of the S1, S2, and S3 values
    :param domain: elevation array
    :param width: cross-shore barrier width
    :param length: alongshore barrier length
    :param time_step: current time step of the storm
    :return truth_array: array of 1s and 0s indicating downhill and uphill (respectively) slopes for future flow
            routing determination
    """
    # ### Calculate Slopes
    for row in range(width):  # for letting sediment out, discharge scenarios 2 and 3
        for col in range(length):
            # if we are not at the last row, do normal calculations
            if row != width - 1:  # uncomment if we are letting sed out the last row
                # tab everything until else statement
                if col > 0:  # i = 0 means there are no cols to the left
                    S1 = (domain[row, col] - domain[row + 1, col - 1]) / (math.sqrt(2))
                    S1 = np.nan_to_num(S1)
                else:
                    S1 = 0

                S2 = domain[row, col] - domain[row + 1, col]
                S2 = np.nan_to_num(S2)
                # slopes_array[TS, d, i] = S2

                if col < (length - 1):  # i at the end length means there are no cols to the right
                    S3 = (domain[row, col] - domain[row + 1, col + 1]) / (math.sqrt(2))
                    S3 = np.nan_to_num(S3)
                else:
                    S3 = 0
            # if at the last row, apply the same slope as the beachface slope
            else:
                if col > 0:  # i = 0 means there are no cols to the left
                    S1 = (domain[row-1, col] - domain[row, col - 1]) / (math.sqrt(2))
                    S1 = np.nan_to_num(S1)
                else:
                    S1 = 0

                S2 = domain[row-1, col] - domain[row, col]
                S2 = np.nan_to_num(S2)
                # slopes_array[TS, d, i] = S2

                if col < (length - 1):  # i at the end length means there are no cols to the right
                    S3 = (domain[row-1, col] - domain[row, col + 1]) / (math.sqrt(2))
                    S3 = np.nan_to_num(S3)
                else:
                    S3 = 0

            if col == 0 and S2 < 0 and S3 < 0:
                # this is greater than the max slope, so no sediment will go to outside
                S1 = -999
            if col == length - 1 and S2 < 0 and S1 < 0:
                S3 = -999
# ----------------------------for singular cell slope analysis----------------------------------------------------------
#             if S1 < 0 and S2 < 0 and S3 < 0:
#                 truth_array[time_step, row, col] = 0
#             else:
#                 truth_array[time_step, row, col] = 1
# ----------------------- group cell slope analysis --------------------------------------------------------------------
            if col == 0:
                avg_slope_array[time_step, row, col] = (S2 + S3) / 2
            elif col == length-1:
                avg_slope_array[time_step, row, col] = (S1 + S2) / 2
            else:
                avg_slope_array[time_step, row, col] = (S1 + S2 + S3) / 3

            # do something with modulus operator
    b_size = 10
    extra_vert_cells = width % b_size
    extra_lat_cells = length % b_size
    if extra_vert_cells == 0:
        n_shifts_vert = int((width - extra_vert_cells) / b_size)
    else:
        n_shifts_vert = int((width - extra_vert_cells) / b_size) + 1
    if extra_lat_cells == 0:
        n_shifts_lat = int((length - extra_lat_cells) / b_size)
    else:
        n_shifts_lat = int((length - extra_lat_cells) / b_size) + 1

    for v in range(n_shifts_vert):
        if v == n_shifts_vert-1 and extra_vert_cells != 0:
            start_row = end_row
            end_row = start_row + extra_vert_cells + 1
        else:
            start_row = v*b_size
            end_row = v*b_size + b_size
        for l in range(n_shifts_lat):
            if l == n_shifts_lat-1 and extra_lat_cells != 0:
                start_col = end_col
                end_col = start_col + extra_lat_cells + 1
            else:
                start_col = l*b_size
                end_col = l*b_size + b_size
            S = np.mean(avg_slope_array[time_step, start_row:end_row, start_col:end_col])
            if S < 0:
                truth_array[time_step, start_row:end_row, start_col:end_col] = 0
            else:
                truth_array[time_step, start_row:end_row, start_col:end_col] = 1

    return truth_array, avg_slope_array


def calculate_slopes(row, col, domain_width, elev_array, domain_length, time_step):
    """
    calculates the slopes at each cell
    :param row: current row of the domain, int
    :param col: current column of the domain, int
    :param domain_width: cross-shore barrier width, int
    :param elev_array: array storing the elevations
    :param domain_length: alongshore barrier length, int
    :param time_step: current time step of the storm
    :return: S1, S2, S3 (floats)
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
        if col > 0:  # i = 0 means there are no cols to the left
            S1 = (elev_array[time_step, row-1, col] - elev_array[time_step, row, col - 1]) / (math.sqrt(2))
            S1 = np.nan_to_num(S1)
        else:
            S1 = 0

        S2 = elev_array[time_step, row-1, col] - elev_array[time_step, row, col]
        S2 = np.nan_to_num(S2)
        # slopes_array[TS, d, i] = S2

        if col < (domain_length - 1):  # i at the end length means there are no cols to the right
            S3 = (elev_array[time_step, row-1, col] - elev_array[time_step, row, col + 1]) / (math.sqrt(2))
            S3 = np.nan_to_num(S3)
        else:
            S3 = 0

    if col == 0 and S2 < 0 and S3 < 0:
        # this is greater than the max slope, so no sediment will go to outside
        S1 = -999
    if col == domain_length - 1 and S2 < 0 and S1 < 0:
        S3 = -999

    return S1, S2, S3


# ### Calculate Discharges
def calculate_discharges(col, S1, S2, S3, Q0, nn, domain_length, max_slope):
    """
    calculates the discharge at each cell
    :param col: current column of the domain
    :param S1: slope to bottom left cell
    :param S2: slope to bottom cell
    :param S3: slope to bottom right cell
    :param Q0: initial discharge
    :param nn: parameter for flow routing
    :param domain_length: alongshore length of the barrier
    :param max_slope: maximum slope that sediment can be transported uphill
    :return: Q1, Q2, Q3 (floats)
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
            datadir,
            time_step_count,
            berm_elev,
            barrier_length,
            sea_level,
            bay_depth,
            interior_domain,
            dune_domain,
            outwash_storm_series,
            # runID,  # eventually this will go away
            # path,  # eventually this will go away
            substep=20,
            sediment_flux_coefficient_Cx=10,
            sediment_flux_coefficient_Ki=2E-3,  # b3d = 7.5E-6 for inundation
            max_slope=-0.25
    ):

        # initial variables
        self._berm_el = berm_elev,  # [dam MHW]
        self._beach_elev = self._berm_el  # [dam MHW]
        self._length = barrier_length  # [dam] length of barrier
        self._substep = substep
        self._max_slope = max_slope
        self._ki = sediment_flux_coefficient_Ki
        self._cx = sediment_flux_coefficient_Cx
        self._sea_level = sea_level  # equal to 0 dam
        self._bay_depth = -bay_depth  # [dam MHW] Depth of bay behind island segment, currently set to 0.3
        self._Si = (np.mean(interior_domain[-10, :]) - np.mean(interior_domain[0, :])) / len(interior_domain)
        # avg_slope = b3d._BermEl / 20    # how it is defined in barrier 3D which is much smaller than when
        # you calculate the slope using the avg of the first and last rows
        # setting up dune domain using b3d
        self._dune_domain = dune_domain
        self._dune_crest = self._dune_domain.max(axis=1)  # dune_crest used to be DuneDomainCrest
        # self._dune_domain = dunes(self._length, self._berm_el, n_rows=2, n_gaps=1, dune_height=0.25)
        # initializing our barrier interior
        # give it an arbitrary width of 30 dam
        # self._interior_domain = np.zeros([30, self._length])
        int_domain = np.flip(interior_domain)
        # initializing our barrier interior
        self._interior_domain = int_domain[5:-1, :]
        # self._interior_domain = int_domain[10:-1, :]
        self._initial_storm_series = np.load(datadir + outwash_storm_series)
        self._final_storm_series = np.load(datadir + outwash_storm_series)  # will be modified 4 substeps, if specified above

        # post-storm (overwash) variables, before outwash modifications
        self._post_b3d_interior_domain = [None] * time_step_count
        self._post_b3d_ave_interior_height = [None] * time_step_count
        self._post_b3d_dunes = [None] * time_step_count
        self._post_storm_x_s = [None] * time_step_count
        self._post_storm_s_sf = [None] * time_step_count

        # output variables
        self._m_beachface = []  # slope of the beach face
        self._Qs_shoreface = np.zeros(time_step_count)

    def update(
            self,
            b3d,
            fudge_fac=1
    ):

        # initialize tracking and other b3d variables
        self._time_index = b3d._time_index
        q_min = b3d._Qs_min  # [m^3 / hr]? Minimum discharge needed for sediment transport (0.001)
        qs_lost_total = 0  # previously OWloss

        # save post-storm dune, interior, shoreline, shoreface params before outwash modifications (a 0.5 yr time step)
        self._post_b3d_interior[self._time_index - 1] = copy.deepcopy(
            b3d.InteriorDomain
        )
        self._post_b3d_dunes[self._time_index - 1] = copy.deepcopy(
            b3d.DuneDomain[self._time_index - 1, :, :]
        )
        self._post_b3d_ave_interior_height[self._time_index - 1] = copy.deepcopy(
            b3d.h_b_TS[-1]
        )
        self._post_storm_x_s[self._time_index - 1] = copy.deepcopy(b3d.x_s)
        self._post_storm_s_sf[self._time_index - 1] = copy.deepcopy(
            b3d.s_sf_TS[-1]
        )

        # only simulate on outwash years (max of one outwash event per model year)
        storm_index = np.argwhere(self._initial_storm_series[:, 0, 0] == self._time_index - 1)
        numstorm = int(len(storm_index))  # we will always only have one storm per year (if that)

        if numstorm > 0:

            # get the hydrograph and duration for this time step
            storm_series = self._initial_storm_series[storm_index]

            # allow for substeps to better simulate morphodynamics (subset=20 equates to 3 min intervals)
            # NOTE: this has been updated to hopefully be used on any storm series for any substep
            # NOTE: LEXI YOU NEED TO MAKE A NEW VARIABLE FOR THE SUBSTEP STORM
            if self._substep != 1:
                updated_bay_levels = bay_converter(storm_series[1], substep=self._substep)
                dur = len(updated_bay_levels)
                storm_series[1] = updated_bay_levels
                storm_series[2] = dur
                self._final_storm_series[storm_index] = storm_series  # save substep version

            else:
                dur = storm_series[2]  # [hr] duration of the storm

            # merge the interior domain, dunes, and create a beach --------------------------------------------------
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

            # flow routing --------------------------------------------------
            # initialize domain variables
            int_width = np.shape(self._interior_domain)[0]
            front_Si = (self._berm_el - np.mean(full_domain[-1, :])) / len(full_domain[int_width+2:-1])
            width = np.shape(full_domain)[0]  # width is the number of rows in the full domain
            domain_array = np.zeros([2, width, self._length])
            domain_array[0, :, :] = full_domain
            # duration = dur * substep  # from previous code
            duration = dur  # we already multiplied dur by the substep above
            Elevation = np.zeros([duration, width, self._length])
            # elevation at the first time step is set to the full domain
            Elevation[0, :, :] = full_domain
            OW_TS = []

            # initialize arrays for flow routing
            Discharge = np.zeros([duration, width, self._length])
            SedFluxIn = np.zeros([duration, width, self._length])
            SedFluxOut = np.zeros([duration, width, self._length])
            elev_change_array = np.zeros([duration, width, self._length])
            truth_array = np.zeros([duration, width, self._length])
            avg_slope_array = np.zeros([duration, width, self._length])
            qs2_array = np.zeros([duration, width, self._length])

            # route the flow
            for TS in range(duration):
                # Begin with elevation from previous timestep
                if TS > 0:
                    Elevation[TS, :, :] = Elevation[TS - 1, :, :]
                    # Elevation[TS, 1:, :] = Elevation[TS - 1, 1:, :]  # in B3D, the first row is dunes, so that is
                    # updated from teh dune domain (I think, line 983 in lexi_b3d)
                print(TS)

                # Need to calculate slopes immediately
                FR_array, avg_slope_array = FR_slopes(truth_array, avg_slope_array, Elevation[TS], width, self._length, TS)

                # get dune crest out here first (dont remember what this means 9/16/2022)
                bayhigh = storm_series[1][TS]  # [dam]

                row = int_width
                # for row in np.arange(int_width, -1, -1):
                #     if max(FR_array[TS, row, :]) > 0 and max(FR_array[TS, row + 1, :]) > 0:
                #         start_row = row
                #         counter = 0
                counter = 0
                while counter == 0 and row >= 0:
                    if max(FR_array[TS, row, :]) > 0 and max(FR_array[TS, row + 1, :]) > 0:
                        row = row - 1
                    else:
                        counter = 1
                start_row = row + 1
                gap_row = Elevation[TS, start_row, :]
                if bayhigh <= min(gap_row):
                    Discharge[TS, :, :] = 0
                else:
                    # ### DUNES
                    OW_TS.append(TS)
                    # ### finding the OW cells --------------------------------------------------------------------
                    Dow = [index for index, value in enumerate(gap_row) if
                           value < bayhigh]  # bayhigh used to be Rhigh
                    gaps = DuneGaps(gap_row, Dow, bayhigh)
                    max_dune = b3d._Dmaxel - b3d._BermEl  # [dam MHW]

                    for q in range(len(gaps)):
                        start = gaps[q][0]
                        stop = gaps[q][1]
                        Rexcess = gaps[q][2]  # (dam)
                        # Calculate discharge through each dune cell
                        overtop_vel = math.sqrt(2 * 9.8 * (Rexcess * 10)) / 10  # (dam/s)
                        overtop_flow = overtop_vel * Rexcess * 3600  # (dam^3/hr)
                        # Set discharge at dune gap
                        Discharge[TS, start_row, start:stop + 1] = overtop_flow  # (dam^3/hr)

                    for d in range(start_row, width):  # for letting sediment out, discharge scenarios 2 and 3
                        # for d in range(int_width + 1, width):  # if we are using scenario 1 discharge
                        # Discharge for each TS, row 1, and all cols set above
                        Discharge[TS, d, :][Discharge[TS, d, :] < 0] = 0
                        # Loop through each col of the specified row
                        for i in range(self._length):
                            # ### Calculate Slopes
                            S1, S2, S3 = calculate_slopes(d, i, width, Elevation, self._length, TS)

                            # if we have discharge, set Qo equal to that value
                            if Discharge[TS, d, i] > 0:
                                Q0 = Discharge[TS, d, i]  # (dam^3/hr)

                                # ---------for all scenarios calc Q vals--------------------------------------------
                                Q1, Q2, Q3 = calculate_discharges(i, S1, S2, S3, Q0,
                                                                  b3d._nn, self._length, self._max_slope)

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
                                if d < int_width:
                                    # C = self._cx * abs(self._Si)  # 10 x the avg slope (from Murray) normal way
                                    C = 0
                                else:
                                    C = self._cx * abs(front_Si[0])

                                if Q1 > q_min:
                                    Qs1 = self._ki * (Q1 * (S1 + C)) ** b3d._mm
                                    if Qs1 < 0:
                                        Qs1 = 0
                                    elif Qs1 > fluxLimit:
                                        Qs1 = fluxLimit
                                else:
                                    Qs1 = 0
                                if Q2 > q_min:
                                    Qs2 = self._ki * (Q2 * (S2 + C)) ** b3d._mm
                                    if Qs2 < 0:
                                        Qs2 = 0
                                    elif Qs2 > fluxLimit:
                                        Qs2 = fluxLimit
                                else:
                                    Qs2 = 0

                                if Q3 > q_min:
                                    Qs3 = self._ki * (Q3 * (S3 + C)) ** b3d._mm
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
                    qs_lost_total = qs_lost_total + sum(
                        SedFluxOut[TS, width - 1, :]) / self._substep  # [dam^3]

            # output variables --------------------------------------------------
            self._Qs_shoreface = qs_lost_total
            print("The sediment volume lost in storm was {0} dam^3".format(round(qs_lost_total, 3)))

            # WE LEFT OFF HERE: updating barrier3d variables based on the flow routing algorithm
            # # update Barrier3D class variables
            # new_ave_interior_height = np.average(
            #     new_xyz_interior_domain[
            #         new_xyz_interior_domain >= barrier3d.SL
            #         ]  # all in dam MHW
            # )
            # barrier3d.h_b_TS[
            #     -1
            # ] = new_ave_interior_height  # slightly altered due to roadway
            # barrier3d.InteriorDomain = new_xyz_interior_domain
            # barrier3d.DomainTS[self._time_index - 1] = new_xyz_interior_domain

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

        # Record storm data
        b3d._StormCount.append(numstorm)
        # return Discharge, elev_change_array, Elevation, qs_lost_total, slopes_array, rexcess_dict, qs2_array, \
        #        self._storm_series, SedFluxOut, SedFluxIn, domain_array, OW_TS, dis_comp_array, FR_array


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
        plt.title("Elevation change (dam)")
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


# ---------------------- Sed out array ---------------------------------------------------------------------------------
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


# ---------------------- Sed out array ---------------------------------------------------------------------------------
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


def plot_dischargeComp(discharge_array, directory, start, stop, bay_level):
    os.chdir(directory)
    newpath = "dis_comparison/"
    if not os.path.exists(newpath):
        os.makedirs(newpath)
    os.chdir(newpath)

    # for t in range(TMAX - 1):
    for t in range(start, stop):
        # Plot and save
        elevFig1 = plt.figure(figsize=(15, 7))
        ax = elevFig1.add_subplot(111)
        x = range(len(discharge_array[0, 0, :]))
        y = np.ones(len(x))*np.mean(discharge_array[t, 0, :])
        y2 = np.ones(len(x))*np.mean(discharge_array[t, 1, :])
        ax.plot(discharge_array[t, 0, :], label="expected discharge")
        ax.plot(discharge_array[t, 1, :], label="actual discharge")
        ax.plot(x, y, label="average expected discharge", color="k", linestyle="dashed")
        ax.plot(x, y2, label="average actual discharge", linestyle="dashed")
        ax.xaxis.set_ticks_position("bottom")
        plt.xlabel("Alongshore Distance (dam)")
        plt.ylabel("Discharge (dam^3/hr)")
        plt.title("Discharge Comparison at the First Dune Line (dam^3/hr)")
        ax.legend(loc="upper left")
        plt.tight_layout()
        full_text = "Time = " + str(t) + "; Bay level = " + str(round(bay_level[t], 3)) + " dam"
        plt.text(0.5, 0.99, full_text, horizontalalignment='center',
        verticalalignment='top', transform=ax.transAxes)
        plt.rcParams.update({"font.size": 15})
        name = "dis_comparison_" + str(t)
        elevFig1.savefig(name, facecolor='w')  # dpi=200
        plt.close(elevFig1)

    frames = []

    for filenum in range(start, stop):
        filename = "dis_comparison_" + str(filenum) + ".png"
        frames.append(imageio.imread(filename))
    imageio.mimsave("discomp.gif", frames, fps=2)
    # imageio.mimsave("sedin.gif", frames, "GIF-FI")
    print()
    print("[ * Discharge comparison GIF successfully generated * ]")


def plot_FRarray(FR_array, directory, start, stop):
    os.chdir(directory)
    newpath = "FRA_cells/"
    if not os.path.exists(newpath):
        os.makedirs(newpath)
    os.chdir(newpath)

    # for t in range(TMAX - 1):
    for t in range(start, stop):
        AnimateDomain = FR_array[t]

        # Plot and save
        elevFig1 = plt.figure(figsize=(15, 7))
        ax = elevFig1.add_subplot(111)
        cax = ax.matshow(
            AnimateDomain,
            # origin="upper",
            cmap="binary",
            # vmin=min_v, vmax=max_v,
        )  # , interpolation='gaussian') # analysis:ignore
        ax.xaxis.set_ticks_position("bottom")
        elevFig1.colorbar(cax)
        plt.xlabel("Alongshore Distance (dam)")
        plt.ylabel("Cross-Shore Distance (dam)")
        plt.title("Cells where flow routing occurs")
        plt.tight_layout()
        timestr = "Time = " + str(t)
        plt.text(1, 1, timestr)
        plt.rcParams.update({"font.size": 15})
        name = "FR_array_" + str(t)
        elevFig1.savefig(name, facecolor='w')  # dpi=200
        plt.close(elevFig1)

    frames = []

    for filenum in range(start, stop):
        filename = "FR_array_" + str(filenum) + ".png"
        frames.append(imageio.imread(filename))
    imageio.mimsave("FR_arrays.gif", frames, fps=2)
    # imageio.mimsave("sedin.gif", frames, "GIF-FI")
    print()
    print("[ * FRA array GIF successfully generated * ]")



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
