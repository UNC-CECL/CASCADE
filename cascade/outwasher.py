import numpy as np
import math
import copy
from matplotlib import pyplot as plt
from collections import Counter

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


def calculate_slopes(
        truth_array,
        avg_slope_array,
        domain,
        width,
        int_width,
        length,
        time_step,
        s1_vals,
        s2_vals,
        s3_vals,
        block_size):
    """
    takes the elevations and differentiates uphill and downhill regions based on an average slope of
    a block of cells
    :param truth_array: empty array that is filled with 1s and 0s based on uphill versus downhill slopes
    :param avg_slope_array: currently empty, for each cell, stores the average value of the S1, S2, and S3 values
    :param domain: elevation array
    :param width: cross-shore barrier width
    :param length: alongshore barrier length
    :param time_step: current time step of the storm
    :param s1_vals: array used to store the S1 slope values for each time step, row, and column
    :param s2_vals: array used to store the S2 slope values for each time step, row, and column
    :param s3_vals: array used to store the S3 slope values for each time step, row, and column
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

            # if col == 0 and S2 < 0 and S3 < 0:
            if col == 0:
                # this is greater than the max slope, so no sediment will go to outside
                S1 = -999
            # if col == length - 1 and S2 < 0 and S1 < 0:
            if col == length - 1:
                S3 = -999

            s1_vals[time_step, row, col] = S1
            s2_vals[time_step, row, col] = S2
            s3_vals[time_step, row, col] = S3

            # averaging S1, S2, and S3 for each cell
            if col == 0:
                avg_slope_array[time_step, row, col] = (S2 + S3) / 2
            elif col == length-1:
                avg_slope_array[time_step, row, col] = (S1 + S2) / 2
            else:
                avg_slope_array[time_step, row, col] = (S1 + S2 + S3) / 3

    # adjusting for a block size that does not divide evenly into the domain
    b_size = block_size
    extra_vert_cells = int_width % b_size  # gives the extra row cells
    extra_lat_cells = length % b_size  # gives the extra column cells

    # calculating how many times we will shift the block to calculate blocks of average slopes
    if extra_vert_cells == 0:
        n_shifts_vert = int(int_width / b_size)
    else:
        n_shifts_vert = int((int_width - extra_vert_cells) / b_size) + 1
    if extra_lat_cells == 0:
        n_shifts_lat = int(length / b_size)
    else:
        n_shifts_lat = int((length - extra_lat_cells) / b_size) + 1

    # Shift through the entire column and then move block up to the next row
    # start right before the dune line
    for v in range(n_shifts_vert):
        if v == 0:
            bot_row = int_width - 1  # start at the bay side dune line and move toward the bay
            top_row = bot_row - b_size
        elif v == n_shifts_vert-1 and extra_vert_cells != 0:  # this is the last shift
            bot_row = top_row
            top_row = bot_row - extra_vert_cells
        else:
            bot_row = top_row
            top_row = bot_row - b_size
        for l in range(n_shifts_lat):
            if l == n_shifts_lat-1 and extra_lat_cells != 0:
                start_col = end_col
                end_col = start_col + extra_lat_cells + 1
            else:
                start_col = l*b_size
                end_col = l*b_size + b_size
            S = np.mean(avg_slope_array[time_step, (top_row+1):(bot_row+1), start_col:end_col])
            if S < 0:
                truth_array[time_step, (top_row+1):(bot_row+1), start_col:end_col] = 0
            else:
                truth_array[time_step, (top_row+1):(bot_row+1), start_col:end_col] = 1

    return truth_array, avg_slope_array, s1_vals, s2_vals, s3_vals


def dune_bay_comparison(dune_flow_type, dune_domain, elev_array, time_step, int_width):
    # determine how we will initiate flow routing with the dunes
    # we will either compare the bay level to the highest dune line or the first dune line
    if dune_flow_type.lower() == "full":
        # find the highest elevation dune row to determine if the bay level is high enough for outwash
        max_dunes_indeces = np.argmax(dune_domain, axis=0)  # returns the max row index for each column
        max_dunes_index = Counter(max_dunes_indeces).most_common(1)[0][0]  # returns most frequent max row
        max_dunes_row = int_width + max_dunes_index  # tells you the row with respect to the whole domain
        max_dune_elevs = elev_array[time_step, max_dunes_row, :]  # gives you the elevation of that row
    else:
        # use the first row of the dune gaps
        max_dune_elevs = elev_array[time_step, int_width, :]

    return max_dune_elevs


def dune_flow_routing_gaps(dune_flow_type, n_dunes, elev_array, time_step, int_width, FR_array, bayhigh):
    if dune_flow_type.lower() == "full":
        for dune in range(n_dunes):
            dune_gap_row = elev_array[time_step, int_width + dune, :]
            Dow = [index for index, value in enumerate(dune_gap_row) if
                   value < bayhigh]  # bayhigh used to be Rhigh
            # assign the dune gaps to the flow routing array
            for val in Dow:
                FR_array[time_step, int_width + dune, val] = 1

        # identify which row we start comparing dune gaps to downhill slopes within the domain
        start_row_connectivity = int_width + n_dunes - 2
        # if the interior is 150 m it extends from row 0 to 149, and row 150 is the first dune row
        # if we add the number of dunes, for example 3 rows of dunes, we are saying the start row is 153, which is
        # none row past the end of the dune line so we need to subtract 1 to get to oceanside dunes (152), BUT
        # we base the rest of the flow routing connectivity to that row, so we start flow routing connectivity one row
        # before that at row 151
    else:
        dune_gap_row = elev_array[time_step, int_width, :]
        Dow = [index for index, value in enumerate(dune_gap_row) if
               value < bayhigh]  # bayhigh used to be Rhigh
        # assign the dune gaps to the flow routing array
        for val in Dow:
            FR_array[time_step, int_width, val] = 1

        # when looking at the first dune line, we will always start one row in front of the dune line
        start_row_connectivity = int_width - 1

    return FR_array, start_row_connectivity


def flow_routing_corrections(truth_array, width, length, time_step, elevation, bayhigh):
    """
    Removes the disconnected downhill cells from the flow routing slopes array.
    :param truth_array: array of 1s and 0s differentiating uphill and downhill slopes
    :param width: the first row that we compare to a sample row
    :param length: alongshore barrier length
    :param time_step: current time step of the storm
    :return: new_truth_array
    """
    start_row = width
    fr_col_array = []
    fr_row_array = []

    for w in range(start_row, -1, -1):
        for l in range(length):
            if truth_array[time_step, w+1, l] == 0 or elevation[time_step, w, l] > bayhigh:
                truth_array[time_step, w, l] = 0

    # now that we have our domain, we are going to find the first flow routing row in each column
    for l in range(length):
        if np.max(truth_array[time_step, :, l]) == 1:  # determine if there is flow routing in this column
            fr_row = np.min(np.where(truth_array[time_step, :, l] == 1))  # find the row flow routing starts in the col
            fr_col = l
            # save the row and columns where fow routing begins
            fr_row_array.append(fr_row)
            fr_col_array.append(fr_col)

    return truth_array, fr_row_array, fr_col_array


def initialize_flow_routing(fr_rows, fr_cols, timestep, elevation, bayhigh, discharge_array):
    start = 0
    i = start

    velocities = []
    flows = []


    while i < (len(fr_cols) - 1):
        adjacent = fr_cols[i + 1] - fr_cols[i]
        if adjacent == 1 and fr_rows[i] == fr_rows[i+1]:
            i = i + 1
        else:
            stop = i
            x = elevation[timestep, fr_rows[i], fr_cols[start]: (fr_cols[stop] + 1)]
            Hmean = sum(x) / float(len(x))
            Rexcess = bayhigh - Hmean
            overtop_vel = math.sqrt(2 * 9.8 * (Rexcess * 10)) / 10  # (dam/s)
            overtop_flow = overtop_vel * Rexcess * 3600  # (dam^3/hr)
            discharge_array[timestep, fr_rows[i], fr_cols[start]:(fr_cols[stop] + 1)] = overtop_flow  # (dam^3/hr)
            overtop_vel_mps = overtop_vel * 10  # m/s
            overtop_flow_cms = overtop_flow / 3600 * 1000  # (m^3/s)
            velocities.append(overtop_vel_mps)
            flows.append(overtop_flow_cms)

            start = stop + 1
            i = start

    # for the last item in the domain
    stop = i
    x = elevation[timestep, fr_rows[i], fr_cols[start]: (fr_cols[stop] + 1)]
    if len(x) > 0:
        Hmean = sum(x) / float(len(x))
        Rexcess = bayhigh - Hmean
        overtop_vel = math.sqrt(2 * 9.8 * (Rexcess * 10)) / 10  # (dam/s)
        overtop_flow = overtop_vel * Rexcess * 3600  # (dam^3/hr)
        discharge_array[timestep, fr_rows[i], fr_cols[start]:(fr_cols[stop] + 1)] = overtop_flow  # (dam^3/hr)
        overtop_vel_mps = overtop_vel * 10  # m/s
        overtop_flow_cms = overtop_flow / 3600 * 1000  # (m^3/s)
        velocities.append(overtop_vel_mps)
        flows.append(overtop_flow_cms)

    return discharge_array, velocities, flows


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
    # all uphill options
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
            outwash_storms_file,
            time_step_count,
            berm_elev,
            barrier_length,
            sea_level,
            bay_depth,
            interior_domain,
            dune_domain,
            substep=20,
            sediment_flux_coefficient_Ki=7.5E-3,  # b3d = 7.5E-6 for inundation
            percent_washout_to_shoreface=100,
            outwash_beach_file=None,
            dune_flow_dynamics="full",
    ):

        # initial variables
        self._percent_washout_to_shoreface = percent_washout_to_shoreface
        self._berm_el = berm_elev,  # [dam MHW]
        self._beach_elev = self._berm_el  # [dam MHW]
        self._length = barrier_length  # [dam] length of barrier
        self._substep = substep
        self._max_slope = -0.25
        self._ki = sediment_flux_coefficient_Ki
        self._cx = 10
        self._sea_level = sea_level  # equal to 0 dam
        self._bay_depth = -bay_depth  # [dam MHW] Depth of bay behind island segment, currently set to 0.3
        self._dune_flow_dynamics = dune_flow_dynamics

        # setting up dune domain using b3d
        self._dune_domain = dune_domain
        self._dune_crest = self._dune_domain.max(axis=1)  # dune_crest used to be DuneDomainCrest
        # initializing our barrier interior
        self._interior_domain = interior_domain
        # loading the bay levels for each time step
        self._outwash_storms = np.load(datadir + outwash_storms_file, allow_pickle=True)
        self._final_bay_levels = []
        self._time_index = 0
        self._outwash_beach = np.load(datadir + outwash_beach_file, allow_pickle=True)

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
        self.velocities = []
        self.flows = []

    def update(
            self,
            b3d
    ):

        ## ADD CHECK FOR OUTWASHER YEAR
        if b3d._time_index - 1 in self._outwash_storms[:, 0]:

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

            storm_index = np.argwhere(self._outwash_storms[:, 0] == self._time_index - 1)
            numstorm = int(len(storm_index))

            if numstorm > 0:
                bay_index = storm_index[0, 0]
                storm_series = self._outwash_storms[bay_index, 1]  # just the bay levels for this outwash year
                # allow for substeps (run the same bay level X times with proportional sed transport)
                # to better simulate morphodynamics
                if self._substep != 1:
                    updated_bay_levels = bay_converter(storm_series, substep=self._substep)
                    storm_series = np.asarray(updated_bay_levels)
                dur = len(storm_series)  # [hr] duration of the storm

                # merge the interior domain, dunes, and create a beach -------------------------------------------------
                if self._outwash_beach is None:
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
                        # the first row of the beach face depends on the last row of the beach
                        if s == 0:
                            beachface_domain[s, :] = beach_domain[-1, 0] - self._m_beachface  # slope of beachface
                        else:
                            beachface_domain[s, :] = beachface_domain[s - 1, :] - self._m_beachface
                else:
                    beach_domain = self._outwash_beach
                    m_beach = np.mean(beach_domain[0, 0] - beach_domain[-1, 0]) / len(beach_domain)

                # the dune domain is being taken from B3D, but is a set of tuples, so it needs to be transposed
                dune_domain_full = np.flip(np.transpose(self._dune_domain) + self._berm_el)
                n_dune_rows = np.shape(dune_domain_full)[0]
                self._full_dunes = copy.deepcopy(dune_domain_full)

                # the full domain of outwasher starts with the interior domain, then the dune, beach, and beachface
                self._interior_domain = np.flip(self._interior_domain)
                full_domain = np.append(self._interior_domain, dune_domain_full, 0)  # [dam MHW]
                full_domain = np.append(full_domain, beach_domain, 0)  # [dam MHW]
                if self._outwash_beach is None:
                    full_domain = np.append(full_domain, beachface_domain, 0)
                self._initial_full_domain = copy.deepcopy(full_domain)

                # flow routing --------------------------------------------------
                # initialize domain variables
                int_width = np.shape(self._interior_domain)[0]

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
                s1_array = np.zeros([duration, width, self._length])
                s2_array = np.zeros([duration, width, self._length])
                s3_array = np.zeros([duration, width, self._length])

                # route the flow
                for TS in range(duration):
                    # Begin with elevation from previous timestep
                    if TS > 0:
                        Elevation[TS, :, :] = Elevation[TS - 1, :, :]  # initial elevation is same as previous TS domain
                    print("Outwasher Time Step: ", TS)

                    if TS == 0:
                        plt.rcParams['figure.figsize'] = (8, 6)
                        plt.rcParams.update({"font.size": 15})
                        fig1 = plt.figure()
                        ax1 = fig1.add_subplot(111)
                        mat = ax1.matshow(
                            Elevation[TS] * 10,
                            cmap="terrain",
                            vmin=-3.0, vmax=3.0,
                        )
                        cbar = fig1.colorbar(mat)
                        cbar.set_label('m MHW', rotation=270, labelpad=15)
                        ax1.set_title("Initial Elevation")
                        ax1.set_ylabel("barrier width (dam)")
                        ax1.set_xlabel("barrier length (dam)")
                        plt.gca().xaxis.tick_bottom()


                    # need to calculate grouped averaged slopes over the domain
                    FR_array, avg_slope_array, s1_array, s2_array, s3_array = calculate_slopes(
                        truth_array,
                        avg_slope_array,
                        Elevation[TS],
                        width,
                        int_width,
                        self._length,
                        TS,
                        s1_array,
                        s2_array,
                        s3_array,
                        block_size=5
                    )

                    # get the hydrograph for this time step
                    bayhigh = storm_series[TS]  # [dam]

                    # determine how we will initiate flow routing with the dunes
                    # we will either compare the bay level to the highest dune line or the first dune line
                    max_dune_elevs = dune_bay_comparison(
                        dune_flow_type=self._dune_flow_dynamics,
                        dune_domain=dune_domain_full,
                        elev_array=Elevation,
                        time_step=TS,
                        int_width=int_width)


                    # determine if the bay level is high enough to overtop the dune gaps
                    if bayhigh <= min(max_dune_elevs):
                        Discharge[TS, :, :] = 0
                    else:
                        self._OW_TS.append(TS)

                        # initialize the flow routing array based on all the dunes or just the first row
                        FR_array, start_row = dune_flow_routing_gaps(
                            dune_flow_type=self._dune_flow_dynamics,
                            n_dunes=n_dune_rows,
                            elev_array=Elevation,
                            time_step=TS,
                            int_width=int_width,
                            FR_array=FR_array,
                            bayhigh=bayhigh)


                        # if TS == 0 or TS == 100:
                        #     fig2 = plt.figure()
                        #     ax2 = fig2.add_subplot(111)
                        #     mat = ax2.matshow(
                        #         FR_array[TS, 0:int_width+n_dune_rows-1, :],
                        #         cmap="binary",
                        #         # vmin=-3.0, vmax=3.0,
                        #     )
                        #     # cbar = fig2.colorbar(mat)
                        #     # cbar.set_label('m MHW', rotation=270, labelpad=15)
                        #     ax2.set_title("Averaged blocks of slopes")
                        #     ax2.set_ylabel("barrier width (dam)")
                        #     ax2.set_xlabel("barrier length (dam)")
                        #     ax2.text(2, 7, 'black = downhill \n white = uphill',
                        #              bbox={'facecolor': 'white', 'pad': 1, 'edgecolor': 'none'})
                        #     plt.gca().xaxis.tick_bottom()

                        # connect the dune gaps to the areas of downhill flow
                        # start at the dune line closest to the ocean
                        FR_array, flow_rows, flow_cols = flow_routing_corrections(
                            FR_array, start_row, self._length, TS, Elevation, bayhigh)

                        max_dune = b3d._Dmaxel - b3d._BermEl  # [dam MHW]

                        # Set discharge at dune gap
                        Discharge, self.velocities, self.flows = initialize_flow_routing(
                            flow_rows,
                            flow_cols,
                            TS,
                            Elevation,
                            bayhigh,
                            Discharge,
                        )

                        # begin flow routing algorithm
                        # this should be min of the flow_rows
                        start_flow_route = min(flow_rows)
                        for d in range(start_flow_route, width):
                            Discharge[TS, d, :][Discharge[TS, d, :] < 0] = 0
                            # Loop through each col of the specified row
                            for i in range(self._length):
                                # ### Calculate Slopes
                                S1 = s1_array[TS, d, i]
                                S2 = s2_array[TS, d, i]
                                S3 = s3_array[TS, d, i]

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
                                        C = self._cx * m_beach

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
                    if all(x <= -0.3 for x in post_outwash_full_domain[0, :]):
                        post_outwash_full_domain = np.delete(post_outwash_full_domain, 0, axis=0)
                    else:
                        check = 0
                # while check == 1:
                #     if all(x <= self._bay_depth for x in post_outwash_full_domain[0, :]):
                #         post_outwash_full_domain = np.delete(post_outwash_full_domain, 0, axis=0)
                #     else:
                #         check = 0

                # domain variables we want to save
                self._full_domain = Elevation[-1]
                post_outwash_interior_domain = Elevation[-1, 0:int_width, :]
                post_outwash_dune_domain = Elevation[-1, int_width:int_width+n_dune_rows, :] - self._berm_el
                post_outwash_beach_domain = Elevation[-1, int_width+n_dune_rows:-1, :]
                self._post_outwash_beach_domain[self._time_index - 1] = post_outwash_beach_domain

                new_ave_interior_height = np.average(
                    post_outwash_interior_domain[
                        post_outwash_interior_domain >= b3d.SL
                        ]  # all in dam MHW
                )
                b3d.h_b_TS[-1] = new_ave_interior_height
                b3d.InteriorDomain = np.flip(post_outwash_interior_domain)
                b3d.DomainTS[self._time_index - 1] = np.flip(post_outwash_interior_domain)
                b3d.DuneDomain[self._time_index - 1, :, :] = copy.deepcopy(
                   np.flip(np.transpose(post_outwash_dune_domain))
                )

                # "nourish" the shoreface with the washout ----------------------------------------------------
                if self._percent_washout_to_shoreface > 0:
                    self._Qs_shoreface[self._time_index - 1] = qs_lost_total * 1000 \
                                                * self._percent_washout_to_shoreface/100  # m^3
                    self._Qs_shoreface_per_length[self._time_index - 1] = (qs_lost_total / self._length) * 100 \
                                                * self._percent_washout_to_shoreface/100  # m^3/m
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
                # self._flow_routing_slopes_array[self._time_index - 1] = FR_s
                self._flow_routing_cellular_array[self._time_index - 1] = FR_array
                self._elevation_change[self._time_index - 1] = ElevationChange



