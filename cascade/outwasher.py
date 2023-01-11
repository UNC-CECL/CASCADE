import numpy as np
import math
import copy

from .beach_dune_manager import shoreface_nourishment


def bay_converter(storms, substep):
    """ repeats the hydrograph value for a specific number of specified times (substep)
    :param storms: list or array of bay levels for this year
    :param substep: number of extra points you want
    :returns new list of bay levels
    """
    new_ss = []
    num_new_vals = substep
    for s in range(len(storms)):
        for a in range(num_new_vals):
            new_ss.append(storms[s])
    return new_ss


def DuneGaps(DuneDomain, Dow, Rhigh):
    """Returns tuple of [gap start index, stop index, avg Rexcess of each gap,
    alpha: ratio of TWL / dune height]
    :param DuneDomain: numpy array with dune elevations [dam MHW]
    :param Dow: list of cell indices where outwash occurs
    :param Rhigh: current bay level [dam MHW]
    """
    gaps = []
    start = 0
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

    stop = i
    x = DuneDomain[Dow[start]: (Dow[stop] + 1)]
    if len(x) > 0:
        Hmean = sum(x) / float(len(x))
        Rexcess = Rhigh - Hmean
        # alpha = Rhigh / (Hmean + bermel)
        gaps.append([Dow[start], Dow[stop], Rexcess])
    return gaps


def FR_slopes(truth_array, avg_slope_array, domain, width, length, time_step, block_size):
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

            # averaging
            if col == 0:
                avg_slope_array[time_step, row, col] = (S2 + S3) / 2
            elif col == length-1:
                avg_slope_array[time_step, row, col] = (S1 + S2) / 2
            else:
                avg_slope_array[time_step, row, col] = (S1 + S2 + S3) / 3

            # do something with modulus operator
    b_size = block_size
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


def flow_routing_corrections(truth_array, width, length, time_step):
    """
    Removes the disconnected downhill cells from the flow routing slopes array.
    :param truth_array: array of 1s and 0s differentiating uphill and downhill slopes
    :param width: cross-shore barrier width
    :param length: alongshore barrier length
    :param time_step: current time step of the storm
    :return: new_truth_array
    """
    row = width
    for w in range(width, -1, -1):
        if np.array_equal(truth_array[time_step, w-1, :], np.ones(length)) == True:
            row = row - 1

    for w in range(row, -1, -1):
        for l in range(length):
            # possibly change this to look at all 3 cells
            if truth_array[time_step, w+1, l] == 0:
                truth_array[time_step, w, l] = 0

    counter = 0
    while counter == 0 and row >= 0:
        if max(truth_array[time_step, row, :]) > 0 and max(truth_array[time_step, row + 1, :]) > 0:
            row = row - 1
        else:
            counter = 1
    start_row = row + 1

    return truth_array, start_row


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
        # >>> from cascade.outwasher import Outwasher
        # >>> outwash = Outwasher()
        # >>> outwash.update(barrier3d)
        """

    def __init__(
            self,
            datadir,
            outwash_years,
            outwash_bay_levels,
            time_step_count,
            berm_elev,
            barrier_length,
            sea_level,
            bay_depth,
            interior_domain,
            dune_domain,
            block_size=5,
            substep=20,
            sediment_flux_coefficient_Cx=10,
            sediment_flux_coefficient_Ki=2E-3,  # b3d = 7.5E-6 for inundation
            max_slope=-0.25,
            shoreface_on=True
    ):

        # initial variables
        self._shoreface_on = shoreface_on
        self._block_size = block_size
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
        # initializing our barrier interior
        int_domain = np.flip(interior_domain)
        # initializing our barrier interior
        self._interior_domain = int_domain
        # loading the bay levels for each time step
        self._outwash_years = np.load(datadir + outwash_years)
        self._initial_bay_levels = np.load(datadir + outwash_bay_levels)
        self._final_bay_levels = []
        self._time_index = 0

        # post-storm (overwash) variables, before outwash modifications
        self._post_b3d_interior_domain = [None] * time_step_count
        self._post_b3d_ave_interior_height = [None] * time_step_count
        self._post_b3d_dunes = [None] * time_step_count
        self._post_storm_x_s = [None] * time_step_count
        self._post_storm_s_sf = [None] * time_step_count

        # output variables
        self._m_beachface = []  # slope of the beach face
        self._OW_TS = []  # array for storing overwashed time steps
        self._initial_full_domain = []
        self._full_dunes = []
        self._full_domain = []
        self._Qs_shoreface = np.zeros(time_step_count)  # dam^3
        self._Qs_shoreface_per_length = np.zeros(time_step_count)  # dam^3/dam
        self._discharge = np.zeros(time_step_count, dtype=object)  # dam^3/substep
        self._elevation_change = np.zeros(time_step_count, dtype=object)
        self._flow_routing_cellular_array = np.zeros(time_step_count, dtype=object)
        self._post_outwash_beach_domain = np.zeros(time_step_count, dtype=object)

    def update(
            self,
            b3d
    ):

        # initialize tracking and other b3d variables
        self._time_index = b3d.time_index
        q_min = b3d._Qs_min  # [m^3 / hr]? Minimum discharge needed for sediment transport (0.001)
        qs_lost_total = 0  # previously OWloss

        # save post-storm dune, interior, shoreline, shoreface params before outwash modifications (a 0.5 yr time step)
        self._post_b3d_interior_domain[self._time_index - 1] = copy.deepcopy(
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
        storm_index = np.argwhere(self._outwash_years[:] == self._time_index - 1)
        numstorm = int(len(storm_index))

        if numstorm > 0:
            bay_index = storm_index[0, 0]
            storm_series = self._initial_bay_levels[bay_index]  # just the bay levels for this outwash year
            # allow for substeps (run the same bay level X times with proportional sed transport)
            # to better simulate morphodynamics
            if self._substep != 1:
                updated_bay_levels = bay_converter(storm_series, substep=self._substep)
                storm_series = np.asarray(updated_bay_levels)
            dur = len(storm_series)  # [hr] duration of the storm

            # merge the interior domain, dunes, and create a beach --------------------------------------------------
            # we added a beach (and "beachface") to the domain with a width of 7 dam based on general beach widths
            beach_domain = np.ones([7, self._length]) * self._beach_elev  # [dam MHW] 7 rows
            beachface_domain = np.zeros([6, self._length])
            # we give the beach slope to be 0.004 m = 0.0004 dam
            m_beach = 0.0004
            # we want the beach to have a slope, but keep the first few rows the berm elevation
            for b in range(len(beach_domain)):
                if b >= 3:
                    beach_domain[b, :] = beach_domain[b - 1, :] - m_beach  # m_beach is positive (downhill)
            self._m_beachface = beach_domain[-1, 0] / len(beachface_domain)  # positive (downhill)
            for s in range(len(beachface_domain)):
                if s == 0:
                    beachface_domain[s, :] = beach_domain[-1, 0] - self._m_beachface  # slope of beachface
                else:
                    beachface_domain[s, :] = beachface_domain[s - 1, :] - self._m_beachface

            # the dune domain is being taken from B3D, but is a set of tuples, so it needs to be transposed
            dune_domain_full = np.flip(np.transpose(self._dune_domain) + self._berm_el)
            self._full_dunes = copy.deepcopy(dune_domain_full)
            # the full domain of outwasher starts with the interior domain, then the dune, beach, and beachface
            full_domain = np.append(self._interior_domain, dune_domain_full, 0)  # [dam MHW]
            full_domain = np.append(full_domain, beach_domain, 0)  # [dam MHW]
            full_domain = np.append(full_domain, beachface_domain, 0)
            self._initial_full_domain = copy.deepcopy(full_domain)

            # flow routing --------------------------------------------------
            # initialize domain variables
            int_width = np.shape(self._interior_domain)[0]
            front_Si = (self._berm_el - np.mean(full_domain[-1, :])) / len(full_domain[int_width+2:-1])  # slope of the back barrier
            width = np.shape(full_domain)[0]  # width is the number of rows in the full domain
            duration = dur  # we already multiplied dur by the substep above
            Elevation = np.zeros([duration, width, self._length])
            # elevation at the first time step is set to the full domain
            Elevation[0, :, :] = full_domain
            FR_array = []
            ElevationChange = 0

            # initialize arrays for flow routing
            Discharge = np.zeros([duration, width, self._length])
            SedFluxIn = np.zeros([duration, width, self._length])
            SedFluxOut = np.zeros([duration, width, self._length])
            truth_array = np.zeros([duration, width, self._length])
            avg_slope_array = np.zeros([duration, width, self._length])

            # route the flow
            for TS in range(duration):
                # Begin with elevation from previous timestep
                if TS > 0:
                    Elevation[TS, :, :] = Elevation[TS - 1, :, :]  # initial elevation is same as previous TS domain
                print(TS)

                # need to calculate grouped averaged slopes over the domain
                FR_array, avg_slope_array = FR_slopes(
                    truth_array,
                    avg_slope_array,
                    Elevation[TS],
                    width,
                    self._length,
                    TS,
                    block_size=self._block_size
                )

                # remove any isolated pockets of downhill slopes
                FR_array, start_row = flow_routing_corrections(FR_array, width, self._length, TS)

                # determine the initial discharge (for each TS, a start row, and all cols) location and flow value
                bayhigh = storm_series[TS]  # [dam]
                gap_row = Elevation[TS, start_row, :]
                if bayhigh <= min(gap_row):
                    Discharge[TS, :, :] = 0
                else:
                    self._OW_TS.append(TS)
                    # find overwashed cells
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
                    for l in range(self._length):
                        if FR_array[TS, start_row, l] == 0:
                            Discharge[TS, start_row, l] = 0

                    # begin flow routing algorithm
                    for d in range(start_row, width):
                        Discharge[TS, d, :][Discharge[TS, d, :] < 0] = 0
                        # Loop through each col of the specified row
                        for i in range(self._length):
                            # ### Calculate Slopes
                            S1, S2, S3 = calculate_slopes(d, i, width, Elevation, self._length, TS)

                            # if we have discharge, set Qo equal to that value
                            if Discharge[TS, d, i] > 0:
                                Q0 = Discharge[TS, d, i]  # (dam^3/hr)
                                Q1, Q2, Q3 = calculate_discharges(i, S1, S2, S3, Q0,
                                                                  b3d._nn, self._length, self._max_slope)

                                ### Update Discharge
                                # discharge is defined for the next row, so we do not need to include the last row
                                # the first row of discharge was already defined
                                if d != width - 1:
                                    # Cell 1
                                    if i > 0:
                                        Discharge[TS, d + 1, i - 1] = Discharge[TS, d + 1, i - 1] + Q1
                                    # Cell 2
                                    Discharge[TS, d + 1, i] = Discharge[TS, d + 1, i] + Q2
                                    # Cell 3
                                    if i < (self._length - 1):
                                        Discharge[TS, d + 1, i + 1] = Discharge[TS, d + 1, i + 1] + Q3

                                # ### Calculate Sed Movement
                                fluxLimit = max_dune  # [dam MHW] dmaxel - bermel
                                # all Qs in [dam^3/hr]
                                if d < int_width:
                                    C = 0
                                else:
                                    C = self._cx * abs(front_Si[0])

                                # bottom left cell
                                if Q1 > q_min:
                                    Qs1 = self._ki * (Q1 * (S1 + C)) ** b3d._mm
                                    if Qs1 < 0:
                                        Qs1 = 0
                                    elif Qs1 > fluxLimit:
                                        Qs1 = fluxLimit
                                else:
                                    Qs1 = 0
                                # bottom center cell
                                if Q2 > q_min:
                                    Qs2 = self._ki * (Q2 * (S2 + C)) ** b3d._mm
                                    if Qs2 < 0:
                                        Qs2 = 0
                                    elif Qs2 > fluxLimit:
                                        Qs2 = fluxLimit
                                else:
                                    Qs2 = 0
                                # bottom right cell
                                if Q3 > q_min:
                                    Qs3 = self._ki * (Q3 * (S3 + C)) ** b3d._mm
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
                                # flux in vs. flux out
                                # sed flux in goes to the next row, and is used for determining flux out at current row
                                # so we need a flux in for the last row, which will be its own variable
                                if d != width - 1:
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

                    # Calculate and save volume of sediment leaving the island for every hour
                    qs_lost_total = qs_lost_total + sum(
                        SedFluxOut[TS, width - 1, :]) / self._substep  # [dam^3]

            # update barrier3d interior and dune domain class variables ---------------------------------------
            # interior domain: remove all rows of bay without any deposition from the domain
            post_outwash_full_domain = Elevation[-1, :, :]  # use the last TS
            check = 1
            while check == 1:
                if all(x <= -self._bay_depth for x in post_outwash_full_domain[0, :]):
                    post_outwash_full_domain = np.delete(post_outwash_full_domain, 0, axis=0)
                else:
                    check = 0

            # domain variables we want to save
            self._full_domain = Elevation[-1]
            post_outwash_interior_domain = Elevation[-1, 0:int_width, :]
            post_outwash_dune_domain = Elevation[-1, int_width:int_width+2, :] - self._berm_el
            post_outwash_beach_domain = Elevation[-1, int_width+2:-1, :]
            self._post_outwash_beach_domain[self._time_index - 1] = post_outwash_beach_domain

            new_ave_interior_height = np.average(
                post_outwash_interior_domain[
                    post_outwash_interior_domain >= b3d._SL
                    ]  # all in dam MHW
            )
            b3d.h_b_TS[-1] = new_ave_interior_height
            b3d.InteriorDomain = np.flip(post_outwash_interior_domain)
            b3d.DomainTS[self._time_index - 1] = np.flip(post_outwash_interior_domain)
            b3d.DuneDomain[self._time_index - 1, :, :] = copy.deepcopy(
               np.flip(np.transpose(post_outwash_dune_domain))
            )

            # "nourish" the shoreface with the washout ----------------------------------------------------
            if self._shoreface_on:
                self._Qs_shoreface[self._time_index - 1] = qs_lost_total * 1000  # m^3
                self._Qs_shoreface_per_length[self._time_index - 1] = (qs_lost_total / self._length) * 100  # m^3/m
            else:
                self._Qs_shoreface[self._time_index - 1] = 0
                self._Qs_shoreface_per_length[self._time_index - 1] = 0
            dummy_beach_width = 0

            (
                b3d.x_s,  # save over class variables
                b3d.s_sf_TS[-1],
                _,  # this is just the change in shoreline position
            ) = shoreface_nourishment(
                b3d.x_s,  # in dam
                b3d.x_t,  # in dam
                self._Qs_shoreface_per_length[self._time_index - 1] / 100,  # convert m^3/m to dam^3/dam
                b3d.h_b_TS[-1],  # in dam
                b3d.DShoreface,  # in dam
                dummy_beach_width,  # in dam
            )
            b3d.x_s_TS[-1] = b3d.x_s

            # other class variables that we want to save
            self._final_bay_levels = storm_series
            self._discharge[self._time_index - 1] = Discharge
            # self._flow_routing_slopes_array[self._time_index - 1] = FR_slopes_array
            self._flow_routing_cellular_array[self._time_index - 1] = FR_array
            self._elevation_change[self._time_index - 1] = ElevationChange
            print("The outwash storm has ended.")

        # return
