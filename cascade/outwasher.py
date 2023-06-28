import numpy as np
import math
import copy
from matplotlib import pyplot as plt

from .beach_dune_manager import shoreface_nourishment, beach_width_dune_dynamics


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
        domain,
        width,
        length,
        time_step,
        s1_vals,
        s2_vals,
        s3_vals,
        sL_left_vals,
        sL_right_vals
):
    """
    takes the elevations and differentiates uphill and downhill regions based on an average slope of
    a block of cells
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
    for row in range(width):
        for col in range(length):
            # if we are not at the last row, do normal calculations
            if row != width - 1:
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
            # if at the last row, apply the same slope as the previous row slopes
            else:
                if col > 0:  # i = 0 means there are no cols to the left
                    S1 = (domain[row - 1, col] - domain[row, col - 1]) / (math.sqrt(2))
                    S1 = np.nan_to_num(S1)
                else:
                    S1 = 0

                S2 = domain[row - 1, col] - domain[row, col]
                S2 = np.nan_to_num(S2)

                if col < (length - 1):  # i at the end length means there are no cols to the right
                    S3 = (domain[row - 1, col] - domain[row, col + 1]) / (math.sqrt(2))
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

            # lateral transport going INTO cells
            # where sL is positive (downhill) into the center cell from the left or right
            # sL_left is coming from the left cell into the center cell
            # sL_right is coming from the right cell into the center cell
            if col == 0:
                sL_left = 0
                sL_right = domain[row, col + 1] - domain[row, col]
            elif col == length - 1:
                sL_left = domain[row, col - 1] - domain[row, col]
                sL_right = 0
            else:
                sL_left = domain[row, col - 1] - domain[row, col]
                sL_right = domain[row, col + 1] - domain[row, col]

            s1_vals[time_step, row, col] = S1
            s2_vals[time_step, row, col] = S2
            s3_vals[time_step, row, col] = S3
            sL_left_vals[time_step, row, col] = sL_left
            sL_right_vals[time_step, row, col] = sL_right

    return s1_vals, s2_vals, s3_vals, sL_left_vals, sL_right_vals


def dune_bay_comparison(dune_flow_type, n_dunes, elev_array, time_step, int_width, bayhigh):
    # determine how we will initiate flow routing with the dunes
    # we will either compare the bay level to the highest dune line or the first dune line
    if dune_flow_type.lower() == "full":
        counter = 0
        for d in range(n_dunes):
            dune_elevs = elev_array[time_step, int_width + d, :]
            if bayhigh > min(dune_elevs):
                counter += 1

        if counter == n_dunes:
            start_bay = True
        else:
            start_bay = False

    else:
        # use the first row of the dune gaps
        dune_elevs = elev_array[time_step, int_width, :]
        if bayhigh > min(dune_elevs):
            start_bay = True
        else:
            start_bay = False

    return start_bay


def dune_flow_routing_gaps(
        dune_flow_type,
        n_dunes,
        elev_array,
        time_step,
        int_width,
        bayhigh,
        discharge_array,
        Q0_array,
        Q1_array,
        Q2_array,
        Q3_array,
        s1_array,
        s2_array,
        s3_array,
        nn,
        max_slope,
        length,
        bay_array,
        downhill_array,
        endcell_array,
):
    if dune_flow_type.lower() == "full":
        for dune in range(n_dunes):
            dune_gap_row = int_width + dune
            dune_gap_elevs = elev_array[time_step, dune_gap_row, :]
            Dow = [index for index, value in enumerate(dune_gap_elevs) if
                   value < bayhigh]  # bayhigh used to be Rhigh
            if len(Dow) > 0:
                # assign the dune gaps to the flow routing array
                start = 0
                i = start
                velocities = []
                flows = []

                while i < (len(Dow) - 1):
                    adjacent = Dow[i + 1] - Dow[i]
                    if adjacent == 1:
                        i = i + 1
                    else:
                        stop = i
                        x = elev_array[time_step, dune_gap_row, Dow[start]: (Dow[stop] + 1)]
                        Hmean = sum(x) / float(len(x))
                        Rexcess = bayhigh - Hmean
                        if Rexcess < 0:
                            Rexcess = 0
                        overtop_vel = math.sqrt(2 * 9.8 * (Rexcess * 10)) / 10  # (dam/s)
                        overtop_flow = overtop_vel * Rexcess * 3600  # (dam^3/hr)
                        discharge_array[time_step, dune_gap_row, Dow[start]:(Dow[stop] + 1)] = overtop_flow
                        Q0_array[time_step, dune_gap_row, Dow[start]:(Dow[stop] + 1)] = discharge_array[time_step,
                                                                                        dune_gap_row,
                                                                                        Dow[start]:(Dow[
                                                                                                        stop] + 1)]  # (dam^3/hr)
                        overtop_vel_mps = overtop_vel * 10  # m/s
                        overtop_flow_cms = overtop_flow / 3600 * 1000  # (m^3/s)
                        velocities.append(overtop_vel_mps)
                        flows.append(overtop_flow_cms)

                        start = stop + 1
                        i = start

                # for the last item in the domain
                stop = i
                x = elev_array[time_step, dune_gap_row, Dow[start]: (Dow[stop] + 1)]
                if len(x) > 0:
                    Hmean = sum(x) / float(len(x))
                    Rexcess = bayhigh - Hmean
                    if Rexcess < 0:
                        Rexcess = 0
                    overtop_vel = math.sqrt(2 * 9.8 * (Rexcess * 10)) / 10  # (dam/s)
                    overtop_flow = overtop_vel * Rexcess * 3600  # (dam^3/hr)
                    discharge_array[time_step, dune_gap_row, Dow[start]:(Dow[stop] + 1)] = overtop_flow  # (dam^3/hr)
                    Q0_array[time_step, dune_gap_row, Dow[start]:(Dow[stop] + 1)] = discharge_array[time_step,
                                                                                    dune_gap_row,
                                                                                    Dow[start]:(Dow[stop] + 1)]
                    overtop_vel_mps = overtop_vel * 10  # m/s
                    overtop_flow_cms = overtop_flow / 3600 * 1000  # (m^3/s)
                    velocities.append(overtop_vel_mps)
                    flows.append(overtop_flow_cms)

                # assign Q1, Q2, Q3 values at the Dow locations without adding to the overall discharge array
                for ow in Dow:
                    bay_array[time_step, dune_gap_row, ow] = 1
                    downhill_array[time_step, dune_gap_row, ow] = 1
                    endcell_array[time_step, dune_gap_row, ow] = 1
                    Q0 = Q0_array[time_step, dune_gap_row, ow]
                    S1 = s1_array[time_step, dune_gap_row, ow]
                    S2 = s2_array[time_step, dune_gap_row, ow]
                    S3 = s3_array[time_step, dune_gap_row, ow]
                    Q1_array[time_step, dune_gap_row, ow], Q2_array[time_step, dune_gap_row, ow], \
                    Q3_array[time_step, dune_gap_row, ow] = calculate_discharges(
                        ow, S1, S2, S3, Q0, nn, length, max_slope)
            else:
                discharge_array[time_step] = 0

    else:
        dune_gap_row = int_width
        dune_gap_elevs = elev_array[time_step, int_width, :]
        Dow = [index for index, value in enumerate(dune_gap_elevs) if
               value < bayhigh]  # bayhigh used to be Rhigh
        # assign the dune gaps to the flow routing array
        # assign the dune gaps to the flow routing array
        start = 0
        i = start
        velocities = []
        flows = []

        while i < (len(Dow) - 1):
            adjacent = Dow[i + 1] - Dow[i]
            if adjacent == 1:
                i = i + 1
            else:
                stop = i
                x = elev_array[time_step, dune_gap_row, Dow[start]: (Dow[stop] + 1)]
                Hmean = sum(x) / float(len(x))
                Rexcess = bayhigh - Hmean
                if Rexcess < 0:
                    Rexcess = 0
                overtop_vel = math.sqrt(2 * 9.8 * (Rexcess * 10)) / 10  # (dam/s)
                overtop_flow = overtop_vel * Rexcess * 3600  # (dam^3/hr)
                discharge_array[time_step, dune_gap_row, Dow[start]:(Dow[stop] + 1)] = overtop_flow  # (dam^3/hr)
                Q0_array[time_step, dune_gap_row, Dow[start]:(Dow[stop] + 1)] = discharge_array[time_step, dune_gap_row,
                                                                                Dow[start]:(Dow[stop] + 1)]
                overtop_vel_mps = overtop_vel * 10  # m/s
                overtop_flow_cms = overtop_flow / 3600 * 1000  # (m^3/s)
                velocities.append(overtop_vel_mps)
                flows.append(overtop_flow_cms)

                start = stop + 1
                i = start

        # for the last item in the domain
        stop = i
        x = elev_array[time_step, dune_gap_row, Dow[start]: (Dow[stop] + 1)]
        if len(x) > 0:
            Hmean = sum(x) / float(len(x))
            Rexcess = bayhigh - Hmean
            if Rexcess < 0:
                Rexcess = 0
            overtop_vel = math.sqrt(2 * 9.8 * (Rexcess * 10)) / 10  # (dam/s)
            overtop_flow = overtop_vel * Rexcess * 3600  # (dam^3/hr)
            discharge_array[time_step, dune_gap_row, Dow[start]:(Dow[stop] + 1)] = overtop_flow  # (dam^3/hr)
            Q0_array[time_step, dune_gap_row, Dow[start]:(Dow[stop] + 1)] = discharge_array[time_step, dune_gap_row,
                                                                            Dow[start]:(Dow[stop] + 1)]
            overtop_vel_mps = overtop_vel * 10  # m/s
            overtop_flow_cms = overtop_flow / 3600 * 1000  # (m^3/s)
            velocities.append(overtop_vel_mps)
            flows.append(overtop_flow_cms)

        for ow in Dow:
            bay_array[time_step, dune_gap_row, ow] = 1
            downhill_array[time_step, dune_gap_row, ow] = 1
            endcell_array[time_step, dune_gap_row, ow] = 1
            Q0 = Q0_array[time_step, dune_gap_row, ow]
            S1 = s1_array[time_step, dune_gap_row, ow]
            S2 = s2_array[time_step, dune_gap_row, ow]
            S3 = s3_array[time_step, dune_gap_row, ow]
            Q1_array[time_step, dune_gap_row, ow], Q2_array[time_step, dune_gap_row, ow], \
            Q3_array[time_step, dune_gap_row, ow] = calculate_discharges(
                ow, S1, S2, S3, Q0, nn, length, max_slope)

    return discharge_array, velocities, flows, Q0_array, Q1_array, Q2_array, Q3_array, bay_array, downhill_array, \
           endcell_array


def check_upstream_discharge(
        start_row,
        timestep,
        elev_array,
        discharge,
        s1_array,
        s2_array,
        s3_array,
        bayhigh,
        length,
        bay_array,
        endcells,
        int_width,
        init_dis
):
    cont = 1
    while start_row >= 1 and cont == 1:
        gaps = np.argwhere(bay_array[timestep, start_row] == 2)
        gap_index = []
        gap_index_left = []
        gap_index_right = []

        for i in range(len(gaps)):
            gap_index.append(gaps[i][0])
            gap_index_left.append(gaps[i][0] - 1)
            gap_index_right.append(gaps[i][0] + 1)

        if -1 in gap_index_left:
            gap_index_left.remove(-1)

        if 50 in gap_index_right:
            gap_index_right.remove(50)

        s1_vals = []
        s2_vals = []
        s3_vals = []
        min_slope = 0

        s1 = s1_array[timestep, start_row - 1, gap_index_right]
        for s in range(len(s1)):
            if s1[s] >= min_slope:
                s1_vals.append(gap_index_right[s])

        s2 = s2_array[timestep, start_row - 1, gap_index]
        for s in range(len(s2)):
            if s2[s] >= min_slope:
                s2_vals.append(gap_index[s])

        s3 = s3_array[timestep, start_row - 1, gap_index_left]
        for s in range(len(s3)):
            if s3[s] >= min_slope:
                s3_vals.append(gap_index_left[s])

        downhill_cells = s1_vals + s2_vals + s3_vals
        downhill_cells = np.asarray(downhill_cells)
        downhill_cells = np.unique(downhill_cells)
        downhill_cells_list = []

        for d in downhill_cells:
            if elev_array[timestep, start_row - 1, d] < bayhigh:
                downhill_cells_list.append(d)
                bay_array[timestep, start_row - 1, d] = 2
                endcells[timestep, start_row - 1, d] = 2

        new_total_cells = len(downhill_cells_list)

        if new_total_cells > 0:
            start_row = start_row - 1
            # cont = 0
        else:
            start_row = start_row
            cont = 0

    # find the starting cells
    fr_row_array = []
    fr_col_array = []
    for l in range(length):
        if np.max(bay_array[timestep, :, l]) == 2:  # determine if there is flow routing in this column
            fr_row = np.min(np.where(bay_array[timestep, :, l] == 2))
            fr_col = l
            # if we are at the first column, just check top and diag right cols
            if fr_col == 0 and bay_array[timestep, fr_row - 1, l] != 2 and bay_array[timestep, fr_row - 1, l + 1] != 2:
                fr_row_array.append(fr_row)
                fr_col_array.append(fr_col)
            # if we are at the last column, just check top and diag left
            elif fr_col == length - 1 and bay_array[timestep, fr_row - 1, l] != 2 and \
                    bay_array[timestep, fr_row - 1, l - 1] != 2:
                fr_row_array.append(fr_row)
                fr_col_array.append(fr_col)
            # otherwise, check all three (top, diag left, diag right)
            elif bay_array[timestep, fr_row - 1, l - 1] != 2 and \
                    bay_array[timestep, fr_row - 1, l] != 2 and \
                    bay_array[timestep, fr_row - 1, l + 1] != 2:
                # save the row and columns where fow routing begins
                fr_row_array.append(fr_row)
                fr_col_array.append(fr_col)

    end_cells = len(fr_row_array)
    if end_cells > 0:
        new_dis = sum(discharge[timestep, int_width, :]) / end_cells
        discharge[timestep] = 0
        discharge[timestep, fr_row_array, fr_col_array] = new_dis
        init_dis[timestep, fr_row_array, fr_col_array] = new_dis
        endcells[timestep, fr_row_array, fr_col_array] = 3

    return discharge, bay_array, endcells, init_dis


def check_underwater(
        start_row,
        timestep,
        elev_array,
        bayhigh,
        length,
        bay_array,
        downhill_array,
        endcell_array,
):
    cont = 1
    while start_row >= 1 and cont == 1:
        gaps = np.argwhere(bay_array[timestep, start_row])
        gap_index = []
        gap_index_left = []
        gap_index_right = []

        for i in range(len(gaps)):
            gap_index.append(gaps[i][0])
            gap_index_left.append(gaps[i][0] - 1)
            gap_index_right.append(gaps[i][0] + 1)

        if -1 in gap_index_left:
            gap_index_left.remove(-1)

        if length in gap_index_right:
            gap_index_right.remove(length)

        connected_cells = gap_index + gap_index_left + gap_index_right
        connected_cells = np.asarray(connected_cells)
        connected_cells = np.unique(connected_cells)
        underwater_cells_list = []

        # add the cells that are above the bay level
        for c in connected_cells:
            if elev_array[timestep, start_row - 1, c] < bayhigh:
                underwater_cells_list.append(c)

        # if there are underwater cells, continue to backtrack through the rows
        new_total_cells = len(underwater_cells_list)

        if new_total_cells > 0:
            bay_array[timestep, start_row - 1, underwater_cells_list] = 1
            downhill_array[timestep, start_row - 1, underwater_cells_list] = 1
            endcell_array[timestep, start_row - 1, underwater_cells_list] = 1
            start_row = start_row - 1
        else:
            cont = 0

    if start_row == 0:
        bay_sufficient = True
    else:
        bay_sufficient = False

    return bay_sufficient, bay_array, downhill_array, endcell_array


def underwater_corrections(
        underwater_array,
        downhill_array,
        endcell_array,
        width,
        TS,
        length
):
    # starting at the second row, check to make sure all bay cells touch a previous bay cell, if not, set cell to 0
    for row in range(1, width):
        gaps = np.argwhere(underwater_array[TS, row])
        gap_index = []
        for i in range(len(gaps)):
            gap_index.append(gaps[i][0])
        if len(gap_index) > 0:
            for g in gap_index:
                if g == 0:
                    if underwater_array[TS, row - 1, g] == 0 and underwater_array[TS, row - 1, g + 1] == 0:
                        underwater_array[TS, row, g] = 0
                        downhill_array[TS, row, g] = 0
                        endcell_array[TS, row, g] = 0
                elif g == length - 1:
                    if underwater_array[TS, row - 1, g] == 0 and underwater_array[TS, row - 1, g - 1] == 0:
                        underwater_array[TS, row, g] = 0
                        downhill_array[TS, row, g] = 0
                        endcell_array[TS, row, g] = 0
                else:
                    if underwater_array[TS, row - 1, g + 1] == 0 and underwater_array[TS, row - 1, g] == 0 and \
                            underwater_array[TS, row - 1, g - 1] == 0:
                        underwater_array[TS, row, g] = 0
                        downhill_array[TS, row, g] = 0
                        endcell_array[TS, row, g] = 0
    return underwater_array, downhill_array, endcell_array


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


def flow_routing_with_lateral_transport(
        duration,
        width,
        int_width,
        length,
        elevation,
        storm_series,
        n_dune_rows,
        dune_flow_dynamics,
        nn,
        max_slope,
        OW_TS,
        Dmaxel,
        BermEl,
        C,
        q_min,
        ki,
        mm,
        k_lat,
        substep,
        qs_lost_total
):
    # initialize arrays for flow routing
    Discharge = np.zeros([duration, width, length])
    SedFluxIn = np.zeros([duration, width, length])
    SedFluxOut = np.zeros([duration, width, length])
    s1_array = np.zeros([duration, width, length])
    s2_array = np.zeros([duration, width, length])
    s3_array = np.zeros([duration, width, length])
    sL_left_array = np.zeros([duration, width, length])
    sL_right_array = np.zeros([duration, width, length])
    Q0_array = np.zeros([duration, width, length])
    Q1_array = np.zeros([duration, width, length])
    Q2_array = np.zeros([duration, width, length])
    Q3_array = np.zeros([duration, width, length])
    QsL_left_in = np.zeros([duration, width, length])
    QsL_right_in = np.zeros([duration, width, length])
    QsL_total = np.zeros([duration, width, length])
    elev_change_array = np.zeros([duration, width, length])
    underwater_array = np.zeros([duration, width, length])
    downhill_array = np.zeros([duration, width, length])
    endcell_array = np.zeros([duration, width, length])
    init_discharge_array = np.zeros([duration, width, length])
    velocities = np.zeros(duration, dtype=object)
    flows = np.zeros(duration, dtype=object)

    # route the flow
    for TS in range(duration):
        # Begin with elevation from previous timestep
        if TS > 0:
            elevation[TS, :, :] = elevation[TS - 1, :, :]  # initial elevation is same as previous TS domain
        # print("Outwasher Time Step: ", TS)
        # print('\r                                                                      ', end='')
        # print("\r", "Outwasher Time Step: ", TS, end="")
        # if TS == duration - 1:
        #     print("\n")

        # get the hydrograph for this time step
        bayhigh = storm_series[TS]  # [dam]

        # determine how we will initiate flow routing with the dunes
        # we will either compare the bay level to the highest dune line or the first dune line
        start_bay = dune_bay_comparison(
            dune_flow_type=dune_flow_dynamics,
            n_dunes=n_dune_rows,
            elev_array=elevation,
            time_step=TS,
            int_width=int_width,
            bayhigh=bayhigh)

        # first, check to see if the bay level exceeds the dune gaps
        if start_bay is False:
            Discharge[TS, :, :] = 0
        else:

            # need to calculate and store slopes over the domain
            s1_array, s2_array, s3_array, sL_left_array, sL_right_array = calculate_slopes(
                domain=elevation[TS],
                width=width,
                length=length,
                time_step=TS,
                s1_vals=s1_array,
                s2_vals=s2_array,
                s3_vals=s3_array,
                sL_left_vals=sL_left_array,
                sL_right_vals=sL_right_array
            )

            # initialize the flow routing array based on all the dunes or just the first row
            [Discharge, velocities[TS], flows[TS], Q0_array, Q1_array, Q2_array, Q3_array,
             underwater_array, downhill_array, endcell_array] = dune_flow_routing_gaps(
                dune_flow_type=dune_flow_dynamics,
                n_dunes=n_dune_rows,
                elev_array=elevation,
                time_step=TS,
                int_width=int_width,
                bayhigh=bayhigh,
                discharge_array=Discharge,
                Q0_array=Q0_array,
                Q1_array=Q1_array,
                Q2_array=Q2_array,
                Q3_array=Q3_array,
                s1_array=s1_array,
                s2_array=s2_array,
                s3_array=s3_array,
                nn=nn,
                max_slope=max_slope,
                length=length,
                bay_array=underwater_array,
                downhill_array=downhill_array,
                endcell_array=endcell_array,
            )

            # we check to see if the dune gaps line up through the entire dune line, if not, we do not
            # consider it a dune gap (if set to full)
            if dune_flow_dynamics.lower() == "full":  # we will not route flow through the dunes
                # start_flow_route = int_width + n_dune_rows - 1
                for col in range(length):
                    if np.any(Discharge[TS, int_width:(int_width + n_dune_rows), col] == 0):
                        Discharge[TS, int_width:(int_width + n_dune_rows), col] = 0
                        underwater_array[TS, int_width:(int_width + n_dune_rows), col] = 0
                        downhill_array[TS, int_width:(int_width + n_dune_rows), col] = 0
                        endcell_array[TS, int_width:(int_width + n_dune_rows), col] = 0

            # set the first row of dunes to 2 for future flow routing
            for col in range(length):
                if Discharge[TS, int_width, col] > 0:
                    downhill_array[TS, int_width, col] = 2
                    endcell_array[TS, int_width, col] = 2

            # check to see if the bay level is high enough to reach the dune gaps through the entire
            # domain (check if we route any flow)
            bay_sufficient, underwater_array, downhill_array, endcell_array = check_underwater(
                start_row=int_width,
                timestep=TS,
                elev_array=elevation,
                bayhigh=bayhigh,
                length=length,
                bay_array=underwater_array,
                downhill_array=downhill_array,
                endcell_array=endcell_array,
            )

            # look for downhill cells connected to the dune gaps
            if bay_sufficient is True:
                OW_TS.append(TS)
                start_row = int_width

                # first remove any dune gaps not connected to the bay cells
                underwater_array, downhill_array, endcell_array = underwater_corrections(
                    underwater_array=underwater_array,
                    downhill_array=downhill_array,
                    endcell_array=endcell_array,
                    width=width,
                    TS=TS,
                    length=length
                )

                Discharge, downhill_array, endcell_array, init_discharge_array = check_upstream_discharge(
                    start_row=start_row,
                    timestep=TS,
                    elev_array=elevation,
                    discharge=Discharge,
                    s1_array=s1_array,
                    s2_array=s2_array,
                    s3_array=s3_array,
                    bayhigh=bayhigh,
                    length=length,
                    bay_array=downhill_array,
                    endcells=endcell_array,
                    int_width=int_width,
                    init_dis=init_discharge_array,
                )

            else:
                Discharge[TS] = 0

            max_dune = Dmaxel - BermEl  # [dam MHW]

            # find min row that has nonzero discharge
            start_flow_row = []
            for col in range(length):
                if np.max(Discharge[TS, :, col]) > 0:
                    start_flow = np.min(np.argwhere(Discharge[TS, :, col]))
                    start_flow_row.append(start_flow)

            if len(start_flow_row) > 0:
                start_flow_route = min(start_flow_row)
                for d in range(start_flow_route, width):
                    for i in range(length):

                        # ### Calculate Slopes
                        S1 = s1_array[TS, d, i]
                        S2 = s2_array[TS, d, i]
                        S3 = s3_array[TS, d, i]

                        # if we have discharge, set Qo equal to that value
                        if Discharge[TS, d, i] > 0:
                            Q0 = Discharge[TS, d, i]  # (dam^3/hr)
                            Q1, Q2, Q3 = calculate_discharges(i, S1, S2, S3, Q0,
                                                              nn, length, max_slope)

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
                                if i < (length - 1):
                                    Discharge[TS, d + 1, i + 1] = Discharge[TS, d + 1, i + 1] + Q3

                            # ### Calculate Sed Movement
                            # ki = self._ki
                            # C = cx * m_beach
                            fluxLimit = max_dune  # [dam MHW] dmaxel - bermel
                            # all Qs in [dam^3/hr]
                            # bottom left cell
                            if Q1 > q_min:
                                Qs1 = ki * (Q1 * (S1 + C)) ** mm
                                if Qs1 < 0:
                                    Qs1 = 0
                                elif Qs1 > fluxLimit:
                                    Qs1 = fluxLimit
                            else:
                                Qs1 = 0
                            # bottom center cell
                            if Q2 > q_min:
                                Qs2 = ki * (Q2 * (S2 + C)) ** mm
                                if Qs2 < 0:
                                    Qs2 = 0
                                elif Qs2 > fluxLimit:
                                    Qs2 = fluxLimit
                            else:
                                Qs2 = 0
                            # bottom right cell
                            if Q3 > q_min:
                                Qs3 = ki * (Q3 * (S3 + C)) ** mm
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
                            # SED OUT CURRENT ROW CELL
                            Qs_out = Qs1 + Qs2 + Qs3

                            limit = -0.3
                            if elevation[TS, d, i] - Qs_out < limit and Qs_out > 0:  # dam
                                new_loss = elevation[TS, d, i] + abs(limit)  # dam
                                Qs1 = (Qs1 / Qs_out) * new_loss
                                Qs2 = (Qs2 / Qs_out) * new_loss
                                Qs3 = (Qs3 / Qs_out) * new_loss
                                Qs_out = Qs1 + Qs2 + Qs3

                            Qs1 = np.nan_to_num(Qs1)
                            Qs2 = np.nan_to_num(Qs2)
                            Qs3 = np.nan_to_num(Qs3)
                            SedFluxOut[TS, d, i] = Qs_out

                            # LATERAL TRANSPORT: still need to add end column constraints
                            # k_lat = self._k_lat  # from one of the runs in Paola and Murray 1997

                            if sL_left_array[TS, d, i] > 0:
                                QsL_left_in[TS, d, i] = k_lat * Qs_out * sL_left_array[TS, d, i]
                                SedFluxOut[TS, d, i - 1] += QsL_left_in[TS, d, i]
                            else:
                                QsL_left_in[TS, d, i] = 0

                            if sL_right_array[TS, d, i] > 0:
                                QsL_right_in[TS, d, i] = k_lat * Qs_out * sL_right_array[TS, d, i]
                                SedFluxOut[TS, d, i + 1] += QsL_right_in[TS, d, i]
                            else:
                                QsL_right_in[TS, d, i] = 0

                            QsL_total[TS, d, i] = QsL_left_in[TS, d, i] + QsL_right_in[TS, d, i]

                            # SED INTO NEXT ROW CELLS
                            if d != width - 1:
                                if i > 0:
                                    SedFluxIn[TS, d + 1, i - 1] += Qs1

                                SedFluxIn[TS, d + 1, i] += Qs2

                                if i < (length - 1):
                                    SedFluxIn[TS, d + 1, i + 1] += Qs3

                            # END OF DOMAIN LOOPS

                # Update Elevation After Every Storm Hour
                ElevationChange = (SedFluxIn[TS] - SedFluxOut[TS] + QsL_total[TS]) / substep
                elevation[TS] = elevation[TS] + ElevationChange
                elev_change_array[TS] = ElevationChange

                # Calculate and save volume of sediment leaving the island for every hour
                qs_lost_total = qs_lost_total + sum(
                    SedFluxOut[TS, width - 1, :]) / substep  # [dam^3]

    return (elevation, qs_lost_total, Discharge, elev_change_array, underwater_array, downhill_array, endcell_array,
            init_discharge_array, OW_TS, velocities, flows)


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
            beta,
            interior_domain,
            dune_domain,
            substep=20,
            sediment_flux_coefficient_Ki=7.5E-3,  # b3d = 7.5E-6 for inundation
            percent_washout_to_shoreface=100,
            outwash_beach_file=None,
            dune_flow_dynamics="full",
            C=0.0134,
            initial_beach_width=0
    ):
        """
         Parameters
         ----------

        initial_beach_width: int, optional
             Initial beach width [m]
        """

        # initial variables
        self._percent_washout_to_shoreface = percent_washout_to_shoreface
        self._berm_el = berm_elev  # [dam MHW]
        self._beach_slope = beta
        self._beach_elev = self._berm_el  # [dam MHW]
        self._length = barrier_length  # [dam] length of barrier
        self._substep = substep
        self._max_slope = -0.25
        self._ki = sediment_flux_coefficient_Ki
        self._k_lat = 0.5
        # self._cx = cx
        self._C = C
        self._mm = 2
        self._sea_level = sea_level  # equal to 0 dam
        self._bay_depth = bay_depth  # [dam MHW] Depth of bay behind island segment, currently set to 0.3
        self._dune_flow_dynamics = dune_flow_dynamics

        # setting up dune domain using b3d
        self._dune_domain = dune_domain  # [dam above berm elevation]
        self._dune_crest = self._dune_domain.max(axis=1)  # dune_crest used to be DuneDomainCrest
        # initializing our barrier interior
        self._interior_domain = interior_domain  # [dam MHW]
        # loading the bay levels for each time step
        self._outwash_storms = np.load(datadir + outwash_storms_file, allow_pickle=True)
        self._final_bay_levels = []
        self._time_index = 0
        self._outwash_beach = np.load(datadir + outwash_beach_file, allow_pickle=True)

        # beach/shoreface "nourishment" variables
        self._beach_width_threshold = 0  # m, triggers dune migration to turn back on
        self._beach_width = [np.nan] * time_step_count
        self._beach_width[0] = initial_beach_width  # m
        self._dune_migration_on = [np.nan] * time_step_count
        self._dune_migration_on[0] = False

        # post-storm (overwash) variables, before outwash modifications
        self._post_b3d_interior_domain = [None] * time_step_count
        self._post_b3d_ave_interior_height = [None] * time_step_count
        self._post_b3d_dunes = [None] * time_step_count
        self._post_storm_x_s = [None] * time_step_count
        self._post_storm_s_sf = [None] * time_step_count

        # output variables
        self._m_beachface = []  # slope of the beach face
        self._OW_TS = []  # array for storing outwashed time steps
        self._outwash_TS = np.zeros(time_step_count, dtype=object)  # array for storing outwash volume (m^3)
        self._outwash_flux_TS = np.zeros(time_step_count, dtype=object)  # array for storing outwash flux (m^3/m)
        self._initial_full_domain = []
        self._full_dunes = []
        self._full_domain = []
        self._Qs_shoreface = np.zeros(time_step_count)  # dam^3
        self._Qs_shoreface_per_length = np.zeros(time_step_count)  # dam^3/dam
        self._discharge = np.zeros(time_step_count, dtype=object)  # dam^3/substep
        self._elevation_change = np.zeros(time_step_count, dtype=object)
        self._underwater_array = np.zeros(time_step_count, dtype=object)
        self._downhill_array = np.zeros(time_step_count, dtype=object)
        self._endcell_array = np.zeros(time_step_count, dtype=object)
        self._post_outwash_beach_domain = np.zeros(time_step_count, dtype=object)
        self._velocities = np.zeros(time_step_count, dtype=object)
        self._flows = np.zeros(time_step_count, dtype=object)
        self._initial_discharge = np.zeros(time_step_count, dtype=object)

    def update(
            self,
            b3d
    ):

        # reduce beach width by the amount of post-storm shoreline change; if the
        # beach width reaches zero, turn dune migration in B3D back on -- otherwise
        # keep it off (we don't want the dune line to prograde
        # because we have fake houses there!)
        self._time_index = b3d.time_index

        change_in_shoreline = (b3d.x_s_TS[-1] - b3d.x_s_TS[-2]) * 10  # m
        self._beach_width[self._time_index - 1] = (
            self._beach_width[self._time_index - 2] - change_in_shoreline
        )
        self._beach_width[self._time_index - 1] = beach_width_dune_dynamics(
            current_beach_width=self._beach_width[self._time_index - 1],
            beach_width_last_year=self._beach_width[self._time_index - 2],
            beach_width_threshold=self._beach_width_threshold,
            barrier3d=b3d,
            time_index=self._time_index,
        )

        # keep track of dune migration
        self._dune_migration_on[self._time_index - 1] = b3d.dune_migration_on

        # CHECK FOR WHETHER IT IS AN OUTWASHER YEAR
        if self._time_index - 1 in self._outwash_storms[:, 0]:
            print("\n start outwash storm")

            # initialize tracking and other b3d variables
            # self._time_index = b3d.time_index
            q_min = b3d._Qs_min  # [m^3 / hr]? Minimum discharge needed for sediment transport (0.001)
            qs_lost_total = 0  # previously OWloss

            # save post-storm dune, interior, shoreline, shoreface params before outwash mods (a 0.5 yr time step)
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
                # we added a beach (and "beachface") to domain with a width of 7 dam based on general beach widths
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
                m_beach = self._beach_slope

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

            # flow routing ----------------------------------------------------------------------------------------
            # initialize domain variables
            int_width = np.shape(self._interior_domain)[0]

            width = np.shape(full_domain)[0]  # width is the number of rows in the full domain
            duration = dur  # we already multiplied dur by the substep above
            Elevation = np.zeros([duration, width, self._length])
            # elevation at the first time step is set to the full domain
            Elevation[0, :, :] = full_domain
            # FR_array = []
            # ElevationChange = 0

            (Elevation, qs_lost_total, Discharge, elev_change_array, underwater_array, downhill_array,
             endcell_array, init_discharge_array, self._OW_TS, self._velocities[self._time_index - 1],
             self._flows[self._time_index - 1],
             ) = flow_routing_with_lateral_transport(
                duration=duration,
                width=width,
                int_width=int_width,
                length=self._length,
                elevation=Elevation,
                storm_series=storm_series,
                n_dune_rows=n_dune_rows,
                dune_flow_dynamics=self._dune_flow_dynamics,
                nn=b3d._nn,
                max_slope=self._max_slope,
                OW_TS=self._OW_TS,
                Dmaxel=b3d._Dmaxel,
                BermEl=b3d._BermEl,
                C=self._C,
                q_min=q_min,
                ki=self._ki,
                mm=self._mm,
                k_lat=self._k_lat,
                substep=self._substep,
                qs_lost_total=qs_lost_total)

            # update barrier3d interior and dune domain class variables for outwash modifications -----------------

            # in barrier3d the elevation array is just the interior domain, here it is the full domain
            # domain variables we want to save
            self._full_domain = Elevation[-1]
            post_outwash_interior_domain = Elevation[-1, 0:int_width, :]
            post_outwash_dune_domain = Elevation[-1, int_width:int_width + n_dune_rows, :] - self._berm_el
            post_outwash_beach_domain = Elevation[-1, int_width + n_dune_rows:-1, :]
            self._post_outwash_beach_domain[self._time_index - 1] = post_outwash_beach_domain

            # interior domain: remove all rows of bay without any deposition from the domain
            # we check the first row rather than the last row (which is used in B3D) because we have not flipped
            # the domain yet
            check = 1
            while check == 1:
                if all(x <= -self._bay_depth for x in post_outwash_interior_domain[0, :]):  # bay_depth = 0.3
                    post_outwash_interior_domain = np.delete(post_outwash_interior_domain, 0, axis=0)
                else:
                    check = 0

            new_ave_interior_height = np.average(
                post_outwash_interior_domain[
                    post_outwash_interior_domain >= b3d.SL
                    ]  # all in dam MHW
            )
            b3d.h_b_TS[-1] = new_ave_interior_height
            # flip the Interior and Dune Domains to put back into Barrier3D
            b3d.InteriorDomain = np.flip(post_outwash_interior_domain)
            b3d.DomainTS[self._time_index - 1] = np.flip(post_outwash_interior_domain)
            # b3d.DuneDomain[self._time_index - 1, :, :] = copy.deepcopy(
            #     np.flip(np.transpose(post_outwash_dune_domain))
            # )
            b3d.DuneDomain[self._time_index - 1, :, :] = np.flip(np.transpose(post_outwash_dune_domain))

            # "nourish" the shoreface with the washout ------------------------------------------------------------
            if self._percent_washout_to_shoreface > 0:
                self._Qs_shoreface[self._time_index - 1] = qs_lost_total * 1000 \
                                                           * self._percent_washout_to_shoreface / 100  # m^3
                self._Qs_shoreface_per_length[self._time_index - 1] = (qs_lost_total / self._length) * 100 \
                                                            * self._percent_washout_to_shoreface / 100  # m^3/m
                (
                    b3d.x_s,  # save over class variables
                    b3d.s_sf_TS[-1],
                    self._beach_width[self._time_index - 1],  # this is just the change in shoreline position
                ) = shoreface_nourishment(
                    b3d.x_s,  # in dam
                    b3d.x_t,  # in dam
                    self._Qs_shoreface_per_length[self._time_index - 1] / 100,  # convert m^3/m to dam^3/dam
                    b3d.h_b_TS[-1],  # in dam
                    b3d.DShoreface,  # in dam
                    self._beach_width[self._time_index - 1] / 10,  # convert m to dam
                )
                self._beach_width[self._time_index - 1] *= 10  # convert dam back to m
                b3d.x_s_TS[-1] = b3d.x_s

                # set dune migration off after nourishment (we don't want the dune line
                # to prograde if the beach width was previously less than threshold)
                b3d.dune_migration_on = False

            # keep track of dune migration
            self._dune_migration_on[self._time_index - 1] = b3d.dune_migration_on

            # final update to outwash and domain variables --------------------------------------------------------
            # recalculate and save DomainWidth and InteriorWidth
            _, _, new_ave_interior_width = b3d.FindWidths(
                b3d.InteriorDomain, b3d.SL
            )
            b3d.InteriorWidth_AvgTS[-1] = new_ave_interior_width

            # replace value of x_b_TS with new InteriorWidth_Avg
            b3d.x_b_TS[-1] = b3d.x_s_TS[-1] + new_ave_interior_width

            # save erosion volumes and fluxes
            vol_m = qs_lost_total * 1000  # dam3 to m3
            flux = qs_lost_total / self._length * 100  # dam3/dam to m3/m
            self._outwash_TS[self._time_index - 1] = vol_m  # array for storing outwash volume (m^3)
            self._outwash_flux_TS[self._time_index - 1] = flux  # array for storing outwash flux (m^3/m)

            # other class variables that we want to save
            self._final_bay_levels = storm_series
            self._discharge[self._time_index - 1] = Discharge
            self._elevation_change[self._time_index - 1] = elev_change_array
            self._underwater_array[self._time_index - 1] = underwater_array
            self._downhill_array[self._time_index - 1] = downhill_array
            self._endcell_array[self._time_index - 1] = endcell_array
            self._initial_discharge[self._time_index - 1] = init_discharge_array

            print(" end outwash storm \n", end="")
