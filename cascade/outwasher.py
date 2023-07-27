"""Simulate outwash events within CASCADE

This module provides functions for modifying a barrier segment from Barrier3D --
consisting of 1+ rows of dune cells, a separate interior grid, and an idealized
shoreface -- to replicate erosion from outwash events,
including:

1. back-barrier, dune, and beach erosion,
2. potential beach and shoreface nourishment at the end of the outwash event,


References
----------

.. [1] Ashton, A. D., & Lorenzo-Trueba, J. (2018). Morphodynamics of barrier response
       to sea-level rise. In Barrier dynamics and response to changing climate
       (pp. 277-304). Springer, Cham.

.. [2] Murray, A. B., & Paola, C. (1994). A cellular model of braided rivers.
        Nature, 371(6492), 54–57.

.. [3] Murray, A. B., & Paola, C. (1997). Properties of a cellular braided-stream model.
        Earth Surface Properties and Landforms, 22(11), 1001–1025.


Notes
-----
Barrier3D does not resolve the beach. Here, we have two beaches:
1. part of the outwash domain that is eroded with the interior and dunes. This beach
recovers between outwash events, as we re-initialize its domain each model time step
that Outwasher is called.
2. representative of shoreline position and is modified dynamically
via nourishment and shoreface dynamics, following the beach nourishment module.

"""


import numpy as np
import math
import copy

from .beach_dune_manager import shoreface_nourishment, beach_width_dune_dynamics


def repeat_hydrograph_for_substep(storms, substep):
    """ repeat each hydrograph value for the specified substep to ensure
    the hydrodynamic forcing is the same for the morphological substep

    :param storms: list or array of bay levels [dam MHW] for this year
    :param substep: number of times we want each value to appear
    :returns new list of bay levels, new_ss [dam MHW]
    """
    new_ss = []
    num_new_vals = substep
    for s in range(len(storms)):
        for a in range(num_new_vals):
            new_ss.append(storms[s])
    return new_ss


def compare_bay_height_to_dune_gap_height(
        dune_flow_type,
        n_dunes,
        elev_array,
        time_step,
        int_width,
        bayhigh
):
    """
    determine how we will initiate flow routing with the dunes
    we will either compare the bay level to all the dune lines or just the first dune line

    :param dune_flow_type: string, "full" or "first" to specify rules for initializing flow routing through the dunes
    :param n_dunes: int, number of dune rows
    :param elev_array: numpy array of elevations [dam]
    :param time_step: int, current time step of the outwash storm
    :param int_width: int, interior barrier width [dam]
    :param bayhigh: float, current bay elevation [dam MHW]
    :return: continue_flow_routing: bool, continue flow routing rules (True) or end this outwash storm time step (False)
    """

    if dune_flow_type.lower() == "full":
        # check both [all] dune rows
        counter = 0
        for d in range(n_dunes):
            dune_elevs = elev_array[time_step, int_width + d, :]
            if bayhigh > min(dune_elevs):
                counter += 1

        # if the bay level exceeds the min of all the dune lines, we will continue flow routing
        if counter == n_dunes:
            continue_flow_routing = True
        # if the bay level does not exceed the min of all the dune lines, we will set the discharge to 0 and nothing
        # else will happen this time step
        else:
            continue_flow_routing = False

    else:
        # only check the first (most landward) row of the dune gaps
        dune_elevs = elev_array[time_step, int_width, :]
        if bayhigh > min(dune_elevs):
            continue_flow_routing = True
        else:
            continue_flow_routing = False

    return continue_flow_routing


def calculate_slopes(
        elev_array,
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
    calculate the slope from each cell to its downstream cells AND
    the lateral slopes from either side of a cell into that cell

    :param elev_array: numpy array of elevations [dam]
    :param width: int, cross-shore barrier width [dam]
    :param length: int, alongshore barrier length [dam]
    :param time_step: int, current time step of the outwash storm
    :param s1_vals: zeros array, used to store the S1 slope values for each time step, row, and column
    :param s2_vals: zeros array, used to store the S2 slope values for each time step, row, and column
    :param s3_vals: zeros array, used to store the S3 slope values for each time step, row, and column
    :param sL_left_vals: zeros array, used to store the lateral slopes coming from the left cell into the center cell
    :param sL_right_vals: zeros array, used to store the lateral slopes coming from the right cell into the center cell
    :return updated (nonzero) s1_vals, s2_vals, s3_vals, sL_left_vals, and sL_right_vals arrays
    """
    # ### Calculate Slopes
    for row in range(width):
        for col in range(length):
            # if we are not at the last row, do normal calculations
            if row != width - 1:
                if col > 0:  # col = 0 means there are no cols to the left
                    S1 = (elev_array[row, col] - elev_array[row + 1, col - 1]) / (math.sqrt(2))
                    S1 = np.nan_to_num(S1)
                else:
                    # prevent the domain from losing sediment on its sides
                    # -999 is greater than the max slope, so no sediment will move that direction
                    S1 = -999

                S2 = elev_array[row, col] - elev_array[row + 1, col]
                S2 = np.nan_to_num(S2)

                if col < (length - 1):  # col = length-1 means there are no cols to the right
                    S3 = (elev_array[row, col] - elev_array[row + 1, col + 1]) / (math.sqrt(2))
                    S3 = np.nan_to_num(S3)
                else:
                    # prevent the domain from losing sediment on its sides
                    # -999 is greater than the max slope, so no sediment will move that direction
                    S3 = -999

            # if at the last row, apply the same slope as the previous row slopes
            else:
                if col > 0:
                    S1 = (elev_array[row - 1, col] - elev_array[row, col - 1]) / (math.sqrt(2))
                    S1 = np.nan_to_num(S1)
                else:
                    S1 = -999

                S2 = elev_array[row - 1, col] - elev_array[row, col]
                S2 = np.nan_to_num(S2)

                if col < (length - 1):  # i at the end length means there are no cols to the right
                    S3 = (elev_array[row - 1, col] - elev_array[row, col + 1]) / (math.sqrt(2))
                    S3 = np.nan_to_num(S3)
                else:
                    S3 = -999

            # lateral transport going INTO cells
            # where sL is positive (downhill) into the center cell from the left or right
            if col == 0:
                # there are no cells to the left, so sediment cannot come from the left
                sL_left = 0
                sL_right = elev_array[row, col + 1] - elev_array[row, col]
            elif col == length - 1:
                # there are no cells to the right, so sediment cannot come from the right
                sL_left = elev_array[row, col - 1] - elev_array[row, col]
                sL_right = 0
            else:
                # the middle cells have adjacent cells on both sides
                sL_left = elev_array[row, col - 1] - elev_array[row, col]
                sL_right = elev_array[row, col + 1] - elev_array[row, col]

            s1_vals[time_step, row, col] = S1
            s2_vals[time_step, row, col] = S2
            s3_vals[time_step, row, col] = S3
            sL_left_vals[time_step, row, col] = sL_left
            sL_right_vals[time_step, row, col] = sL_right

    return s1_vals, s2_vals, s3_vals, sL_left_vals, sL_right_vals


def calculate_discharge_at_dune_gaps(
        dune_flow_type,
        n_dunes,
        elev_array,
        time_step,
        int_width,
        bayhigh,
        discharge_array,
        bay_array,
        downhill_array,
        start_cell_array,
):
    """
    calculate the expected discharge at the dune gaps from gravity driven flow

    :param dune_flow_type: string, "full" or "first" to specify rules for initializing flow routing through the dunes
    :param n_dunes: int, number of dune rows
    :param elev_array: numpy array of elevations [dam]
    :param time_step: int, current time step of the outwash storm
    :param int_width: int, interior barrier width [dam]
    :param bayhigh: float, current bay elevation [dam MHW]
    :param discharge_array: zeros array, used to store the total discharge (sum of Q1, Q2, Q3) for each cell [dam^3]
    :param bay_array: zeros array, stores the values for each cell (1 for underwater, 2 for downhill, 3 for start cell)
    to indicate flow routing initialization
    :param downhill_array: zeros array, builds on the underwater array,
    indicates the cells that are both underwater and downhill
    :param start_cell_array: zeros array, builds on the underwater and downhill array,
    indicates the cells that are underwater, downhill, and start cells for flow routing
    :return: updated discharge_array, bay_array, downhill_array, and start_cell_array
    """
    if dune_flow_type.lower() == "full":
        # slightly altered Barrier3D's DuneGaps function
        for dune in range(n_dunes):
            dune_gap_row = int_width + dune
            dune_gap_elevs = elev_array[time_step, dune_gap_row, :]
            # find and return (Dow) the index of the dune cells that are lower than the bay elevation
            Dow = [index for index, value in enumerate(dune_gap_elevs) if
                   value < bayhigh]  # bayhigh used to be Rhigh
            if len(Dow) > 0:
                start = 0
                i = start
                velocities = []  # not used in calculations, but here for debugging
                flows = []  # not used in calculations, but here for debugging

                # check for a section of adjacent submerged dune cells acting as one dune gap
                while i < (len(Dow) - 1):
                    adjacent = Dow[i + 1] - Dow[i]
                    if adjacent == 1:
                        i = i + 1
                    else:
                        stop = i
                        # calculate the average elevation of the dune gap across all its cells
                        x = elev_array[time_step, dune_gap_row, Dow[start]: (Dow[stop] + 1)]
                        Hmean = sum(x) / float(len(x))
                        # find the bay level above the average dune gap elevation
                        Rexcess = bayhigh - Hmean
                        # if Rexcess < 0:
                        #     Rexcess = 0
                        # calculate the velocity (based on gravitational flow) and flows (Q=VA)
                        overtop_vel = math.sqrt(2 * 9.8 * (Rexcess * 10)) / 10  # (dam/s)
                        overtop_flow = overtop_vel * Rexcess * 3600  # (dam^3/hr)
                        # assign the flow to each cell in the dune gap
                        discharge_array[time_step, dune_gap_row, Dow[start]:(Dow[stop] + 1)] = overtop_flow

                        # convert velocity to m/s and flow to cubic m/s
                        overtop_vel_mps = overtop_vel * 10  # m/s
                        overtop_flow_cms = overtop_flow / 3600 * 1000  # (m^3/s)
                        velocities.append(overtop_vel_mps)
                        flows.append(overtop_flow_cms)

                        # start searching for dune gap at the next cell
                        start = stop + 1
                        i = start

                # for the last item in the domain
                stop = i
                x = elev_array[time_step, dune_gap_row, Dow[start]: (Dow[stop] + 1)]
                if len(x) > 0:
                    Hmean = sum(x) / float(len(x))
                    Rexcess = bayhigh - Hmean
                    # if Rexcess < 0:
                    #     Rexcess = 0
                    overtop_vel = math.sqrt(2 * 9.8 * (Rexcess * 10)) / 10  # (dam/s)
                    overtop_flow = overtop_vel * Rexcess * 3600  # (dam^3/hr)
                    discharge_array[time_step, dune_gap_row, Dow[start]:(Dow[stop] + 1)] = overtop_flow  # (dam^3/hr)

                    # convert velocity to m/s and flow to cubic m/s
                    overtop_vel_mps = overtop_vel * 10  # m/s
                    overtop_flow_cms = overtop_flow / 3600 * 1000  # (m^3/s)
                    velocities.append(overtop_vel_mps)
                    flows.append(overtop_flow_cms)

                # calculate Q1, Q2, Q3 values at the Dow locations without adding to the overall discharge array
                for ow in Dow:
                    bay_array[time_step, dune_gap_row, ow] = 1
                    # should we set the downhill and start_cell array here??
                    downhill_array[time_step, dune_gap_row, ow] = 1
                    start_cell_array[time_step, dune_gap_row, ow] = 1
            else:
                # if none of the dune cells are overwashed (bay level too low), we set all the discharge to 0
                # for this outwash storm time step
                discharge_array[time_step] = 0

    else:
        # only base on the first (most landward) dune row
        dune_gap_row = int_width
        dune_gap_elevs = elev_array[time_step, int_width, :]
        Dow = [index for index, value in enumerate(dune_gap_elevs) if
               value < bayhigh]  # bayhigh used to be Rhigh
        if len(Dow) > 0:
            start = 0
            i = start
            velocities = []  # not used in calculations, but here for debugging
            flows = []  # not used in calculations, but here for debugging

            # check for a section of adjacent submerged dune cells acting as one dune gap
            while i < (len(Dow) - 1):
                adjacent = Dow[i + 1] - Dow[i]
                if adjacent == 1:
                    i = i + 1
                else:
                    stop = i
                    # calculate the average elevation of the dune gap across all its cells
                    x = elev_array[time_step, dune_gap_row, Dow[start]: (Dow[stop] + 1)]
                    Hmean = sum(x) / float(len(x))
                    # find the bay level above the average dune gap elevation
                    Rexcess = bayhigh - Hmean
                    # if Rexcess < 0:
                    #     Rexcess = 0
                    # calculate the velocity (based on gravitational flow) and flows (Q=VA)
                    overtop_vel = math.sqrt(2 * 9.8 * (Rexcess * 10)) / 10  # (dam/s)
                    overtop_flow = overtop_vel * Rexcess * 3600  # (dam^3/hr)
                    # assign the flow to each cell in the dune gap
                    discharge_array[time_step, dune_gap_row, Dow[start]:(Dow[stop] + 1)] = overtop_flow
                    # convert velocity to m/s and flow to cubic m/s
                    overtop_vel_mps = overtop_vel * 10  # m/s
                    overtop_flow_cms = overtop_flow / 3600 * 1000  # (m^3/s)
                    velocities.append(overtop_vel_mps)
                    flows.append(overtop_flow_cms)
                    # start searching for dune gap at the next cell
                    start = stop + 1
                    i = start

            # for the last item in the domain
            stop = i
            x = elev_array[time_step, dune_gap_row, Dow[start]: (Dow[stop] + 1)]
            if len(x) > 0:
                Hmean = sum(x) / float(len(x))
                Rexcess = bayhigh - Hmean
                # if Rexcess < 0:
                #     Rexcess = 0
                overtop_vel = math.sqrt(2 * 9.8 * (Rexcess * 10)) / 10  # (dam/s)
                overtop_flow = overtop_vel * Rexcess * 3600  # (dam^3/hr)
                discharge_array[time_step, dune_gap_row, Dow[start]:(Dow[stop] + 1)] = overtop_flow  # (dam^3/hr)
                overtop_vel_mps = overtop_vel * 10  # m/s
                overtop_flow_cms = overtop_flow / 3600 * 1000  # (m^3/s)
                velocities.append(overtop_vel_mps)
                flows.append(overtop_flow_cms)

            for ow in Dow:
                bay_array[time_step, dune_gap_row, ow] = 1
                # should we set the downhill and start_cell array here??
                downhill_array[time_step, dune_gap_row, ow] = 1
                start_cell_array[time_step, dune_gap_row, ow] = 1

    return discharge_array, bay_array, downhill_array, start_cell_array


def calculate_discharges(
        col,
        S1,
        S2,
        S3,
        Q0,
        nn,
        domain_length,
        max_slope
):
    """
    calculates the discharge at each cell following Murray and Paola (1994, 1997)

    :param col: current column of the domain
    :param S1: slope to bottom left cell
    :param S2: slope to bottom cell
    :param S3: slope to bottom right cell
    :param Q0: initial discharge [dam^3]
    :param nn: float = 0.5, flow parameter
    :param domain_length: int, alongshore barrier length [dam]
    :param max_slope: float, max uphill slope value that still allows sediment transport
    :return: Q1, Q2, Q3: floats, [dam^3]
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


def check_underwater(
        start_row,
        time_step,
        elev_array,
        bayhigh,
        length,
        bay_array,
        downhill_array,
        start_cell_array,
):
    """
    check for hydraulic connection from the bay to the dune gaps: there must be a path of submerged cells

    :param start_row: int, most landward dune row (equal to the interior width)
    :param time_step: int, current time step of the outwash storm
    :param elev_array: numpy array of elevations [dam]
    :param bayhigh: float, current bay elevation [dam MHW]
    :param length: int, alongshore barrier length [dam]
    :param bay_array: stores the values for each cell (1 for underwater, 2 for downhill, 3 for start cell)
    to indicate flow routing initialization
    :param downhill_array: builds on the underwater array,
    indicates the cells that are both underwater and downhill
    :param start_cell_array: builds on the underwater and downhill array,
    indicates the cells that are underwater, downhill, and start cells for flow routing
    :return: bay_sufficient (boolean indicating if there is a hydraulic connection from the bay to the dunes)
    :return: updated bay_array, downhill_array, and start_cell_array
    """
    cont = 1
    current_row = start_row
    # backtrack from the dunes to the bay (0th row) to see if there are underwater cells connected to the dunes
    # here we start at the dune gaps and look upstream (one cell diagonally left, center, and diagonally right)
    while current_row >= 1 and cont == 1:
        # find the submerged cells in our row (nonzero in bay_array)
        gaps = np.argwhere(bay_array[time_step, current_row])
        gap_index = []
        gap_index_left = []
        gap_index_right = []

        # get the left and right indices of the submerged cells (hydraulically connected to our submerged cell)
        for i in range(len(gaps)):
            gap_index.append(gaps[i][0])
            gap_index_left.append(gaps[i][0] - 1)
            gap_index_right.append(gaps[i][0] + 1)

        # in case our gap exists at the first or last column, remove the left and right cells from the list
        if -1 in gap_index_left:
            gap_index_left.remove(-1)

        if length in gap_index_right:
            gap_index_right.remove(length)

        # create one list of all the UNIQUE indices that are connected to our submerged cells
        connected_cells = gap_index + gap_index_left + gap_index_right
        connected_cells = np.asarray(connected_cells)
        connected_cells = np.unique(connected_cells)
        underwater_cells_list = []

        # find the cells from the upstream row that are submerged
        for c in connected_cells:
            if elev_array[time_step, current_row - 1, c] < bayhigh:
                underwater_cells_list.append(c)

        new_total_cells = len(underwater_cells_list)
        # if there are submerged cells, continue to backtrack through the rows
        if new_total_cells > 0:
            # assign the upstream, submerged cells their value of 1
            bay_array[time_step, current_row - 1, underwater_cells_list] = 1
            downhill_array[time_step, current_row - 1, underwater_cells_list] = 1
            start_cell_array[time_step, current_row - 1, underwater_cells_list] = 1
            current_row = current_row - 1
        else:
            # if there are no submerged cells in the upstream row, stop backtracking
            cont = 0

    # in order to route flow, we must have submerged cells to the 0th row
    if current_row == 0:
        bay_sufficient = True
    else:
        bay_sufficient = False

    return bay_sufficient, bay_array, downhill_array, start_cell_array


def clear_floating_submerged_cells(
        underwater_array,
        downhill_array,
        start_cell_array,
        width,
        time_step,
        length
):
    """
    remove any "floating" submerged cells. these are submerged cells that started at a dune gap, but did not reach the
    bay. the previous function, check_underwater, only checks if there is A dune gap that is hydraulically connected
    to the bay; it does not check every dune gap. this function clears any submerged cells that are not connected to
    the bay

    :param underwater_array: indicates the cells that are underwater (1 in the array)
    :param downhill_array: builds on the underwater array,
    indicates the cells that are both underwater and downhill
    :param start_cell_array: builds on the underwater and downhill array,
    indicates the cells that are underwater, downhill, and start cells for flow routing
    :param width: int, cross-shore barrier width [dam]
    :param time_step: int, current time step of the outwash storm
    :param length: int, alongshore barrier length [dam]
    :return: updated underwater_array, downhill_array, start_cell_array
    """
    # starting at the second row, check to make sure all bay cells touch a previous bay cell, if not, set cell to 0
    for row in range(1, width):
        # find the nonzero (submerged) cells
        submerged_cells = np.argwhere(underwater_array[time_step, row])
        submerged_cells_index = []
        for i in range(len(submerged_cells)):
            submerged_cells_index.append(submerged_cells[i][0])
        # if len(submerged_cells_index) > 0:
        for g in submerged_cells_index:
            if g == 0:
                # at the left-most column, only check center and diagonal right cells
                if underwater_array[time_step, row - 1, g] == 0 and underwater_array[time_step, row - 1, g + 1] == 0:
                    underwater_array[time_step, row, g] = 0
                    downhill_array[time_step, row, g] = 0
                    start_cell_array[time_step, row, g] = 0
            elif g == length - 1:
                # at the right-most column, only check center and diagonal left cells
                if underwater_array[time_step, row - 1, g] == 0 and underwater_array[time_step, row - 1, g - 1] == 0:
                    underwater_array[time_step, row, g] = 0
                    downhill_array[time_step, row, g] = 0
                    start_cell_array[time_step, row, g] = 0
            else:
                # in middle columns, check diagonal left, center and diagonal right cells
                if underwater_array[time_step, row - 1, g + 1] == 0 and underwater_array[time_step, row - 1, g] == 0 and \
                        underwater_array[time_step, row - 1, g - 1] == 0:
                    underwater_array[time_step, row, g] = 0
                    downhill_array[time_step, row, g] = 0
                    start_cell_array[time_step, row, g] = 0
    return underwater_array, downhill_array, start_cell_array


def initialize_discharge_at_start_cells(
        start_row,
        time_step,
        elev_array,
        discharge,
        s1_array,
        s2_array,
        s3_array,
        bayhigh,
        length,
        bay_array,
        start_cell_array,
        int_width,
        init_dis
):
    """
    starting at the most landward dune line, find the cells that have downhill slopes toward the dune gaps, and
    set the initial discharge at the selected start cells using a conservation of flow rule

    :param start_row: int, equal to the interior width
    :param time_step: int, current time step of the outwash storm
    :param elev_array: numpy array of elevations [dam]
    :param discharge: numpy array, stores the total discharge (sum of Q1, Q2, Q3) for each cell [dam^3]
    :param s1_array: numpy array, S1 slope values for each time step, row, and column
    :param s2_array: numpy array, S2 slope values for each time step, row, and column
    :param s3_array: numpy array, S3 slope values for each time step, row, and column
    :param bayhigh: float, current bay elevation [dam MHW]
    :param length: int, alongshore barrier length [dam]
    :param bay_array: stores the values for each cell (1 for underwater, 2 for downhill, 3 for start cell)
    to indicate flow routing initialization (EQUAL TO DOWNHILL ARRAY HERE)
    :param start_cell_array: builds on the underwater and downhill array,
    indicates the cells that are underwater, downhill, and start cells for flow routing
    :param int_width: int, interior barrier width [dam]
    :param init_dis: zeros array, stores the initial discharge values after specifying the start cells
    used for plotting purposes only, can be deleted
    :return: updated discharge, bay_array, start_cell_array, and init_dis
    """
    cont = 1
    current_row = start_row

    # find the downhill slopes, starting at the most landward dune line
    while current_row >= 1 and cont == 1:
        downhill_cells = np.argwhere(bay_array[time_step, current_row] == 2)
        downhill_index = []
        downhill_index_left = []
        downhill_index_right = []

        for i in range(len(downhill_cells)):
            downhill_index.append(downhill_cells[i][0])
            downhill_index_left.append(downhill_cells[i][0] - 1)
            downhill_index_right.append(downhill_cells[i][0] + 1)

        # in case our cell exists at the first or last column, remove the left and right cells from the list
        if -1 in downhill_index_left:
            downhill_index_left.remove(-1)

        if 50 in downhill_index_right:
            downhill_index_right.remove(50)

        # we will store the slope values that are downhill (positive) in these arrays
        s1_vals = []
        s2_vals = []
        s3_vals = []
        # since we only want downhill or zero slope values, the minimum slope allowable is 0
        min_slope = 0

        # get the previously calculated S1, S2, and S3 values for the arrays (left, right, and regular index)
        # since we are looking upstream from the current cell, the S1 value comes from the right index and the S3 value
        # comes from the left index
        s1 = s1_array[time_step, current_row - 1, downhill_index_right]
        for s in range(len(s1)):
            if s1[s] >= min_slope:
                s1_vals.append(downhill_index_right[s])

        s2 = s2_array[time_step, current_row - 1, downhill_index]
        for s in range(len(s2)):
            if s2[s] >= min_slope:
                s2_vals.append(downhill_index[s])

        s3 = s3_array[time_step, current_row - 1, downhill_index_left]
        for s in range(len(s3)):
            if s3[s] >= min_slope:
                s3_vals.append(downhill_index_left[s])

        # create one list of all the UNIQUE cells that are downhill
        downhill_cells = s1_vals + s2_vals + s3_vals
        downhill_cells = np.asarray(downhill_cells)
        downhill_cells = np.unique(downhill_cells)
        submerged_downhill_cells_list = []

        # checking that these cells are also submerged
        for d in downhill_cells:
            if bay_array[time_step, current_row - 1, d] > 0:
                submerged_downhill_cells_list.append(d)
                bay_array[time_step, current_row - 1, d] = 2
                start_cell_array[time_step, current_row - 1, d] = 2

        new_total_cells = len(submerged_downhill_cells_list)

        # if we have submerged and downhill cells, continue backtracking
        if new_total_cells > 0:
            current_row = current_row - 1
        else:
            # if we do not, then stop backtracking and move onto locating the start cells for flow routing
            cont = 0

    # find the starting cells: the most landward, downhill cells in each column not hydraulically connected to a more
    # landward, downhill cell
    start_cell_row_array = []
    start_cell_col_array = []
    for l in range(length):
        if np.max(bay_array[time_step, :, l]) == 2:  # determine if there is a downhill cell in this column
            start_cell_row = np.min(np.where(bay_array[time_step, :, l] == 2))  # find the most landward cell row
            start_cell_col = l
            # if we are at the first column, just check if upstream diag right column is downhill
            if start_cell_col == 0 and bay_array[time_step, start_cell_row - 1, l + 1] != 2:
                # save the row and column of the cell
                start_cell_row_array.append(start_cell_row)
                start_cell_col_array.append(start_cell_col)
            # if we are at the last column, just check if upstream diag left column is downhill
            elif start_cell_col == length - 1 and bay_array[time_step, start_cell_row - 1, l - 1] != 2:
                # save the row and column of the cell
                start_cell_row_array.append(start_cell_row)
                start_cell_col_array.append(start_cell_col)
            # otherwise, check both (diag left, diag right)
            elif bay_array[time_step, start_cell_row - 1, l - 1] != 2 and \
                    bay_array[time_step, start_cell_row - 1, l + 1] != 2:
                # save the row and column of the cell
                start_cell_row_array.append(start_cell_row)
                start_cell_col_array.append(start_cell_col)

    # count the number of start cells for flow distribution
    start_cells = len(start_cell_row_array)
    if start_cells > 0:
        # apply conservation of flow from the discharge expected at the dune gaps to the start cells
        new_dis = sum(discharge[time_step, int_width, :]) / start_cells
        discharge[time_step] = 0  # re-set our discharge array
        discharge[time_step, start_cell_row_array, start_cell_col_array] = new_dis
        init_dis[time_step, start_cell_row_array, start_cell_col_array] = new_dis
        start_cell_array[time_step, start_cell_row_array, start_cell_col_array] = 3

    return discharge, bay_array, start_cell_array, init_dis


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
    """
    perform flow routing and sediment transport

    :param duration:
    :param width:
    :param int_width:
    :param length:
    :param elevation:
    :param storm_series:
    :param n_dune_rows:
    :param dune_flow_dynamics:
    :param nn:
    :param max_slope:
    :param OW_TS:
    :param Dmaxel:
    :param BermEl:
    :param C:
    :param q_min:
    :param ki:
    :param mm:
    :param k_lat:
    :param substep:
    :param qs_lost_total:
    :return:
    """
    # initialize arrays for flow routing
    Discharge = np.zeros([duration, width, length])
    SedFluxIn = np.zeros([duration, width, length])
    SedFluxOut = np.zeros([duration, width, length])
    s1_array = np.zeros([duration, width, length])
    s2_array = np.zeros([duration, width, length])
    s3_array = np.zeros([duration, width, length])
    sL_left_array = np.zeros([duration, width, length])
    sL_right_array = np.zeros([duration, width, length])
    QsL_left_in = np.zeros([duration, width, length])
    QsL_right_in = np.zeros([duration, width, length])
    QsL_total = np.zeros([duration, width, length])
    elev_change_array = np.zeros([duration, width, length])
    underwater_array = np.zeros([duration, width, length])
    downhill_array = np.zeros([duration, width, length])
    start_cell_array = np.zeros([duration, width, length])
    init_discharge_array = np.zeros([duration, width, length])  # used for plotting

    # route the flow
    for TS in range(duration):
        # Begin with elevation from previous timestep
        if TS > 0:
            elevation[TS, :, :] = elevation[TS - 1, :, :]  # initial elevation is same as previous TS domain

        # get the hydrograph for this time step
        bayhigh = storm_series[TS]  # [dam]

        # determine how we will initiate flow routing with the dunes
        # we will either compare the bay level to the highest dune line or the first (most landward) dune line
        continue_flow_routing = compare_bay_height_to_dune_gap_height(
            dune_flow_type=dune_flow_dynamics,
            n_dunes=n_dune_rows,
            elev_array=elevation,
            time_step=TS,
            int_width=int_width,
            bayhigh=bayhigh)

        # first, check to see if the bay level exceeds the dune gaps
        if continue_flow_routing is False:
            Discharge[TS, :, :] = 0
        else:
            # need to calculate and store slopes over the domain
            s1_array, s2_array, s3_array, sL_left_array, sL_right_array = calculate_slopes(
                elev_array=elevation[TS],
                width=width,
                length=length,
                time_step=TS,
                s1_vals=s1_array,
                s2_vals=s2_array,
                s3_vals=s3_array,
                sL_left_vals=sL_left_array,
                sL_right_vals=sL_right_array
            )

            # initialize the flow routing array based on all the dunes or just the first (most landward) row
            [Discharge, underwater_array, downhill_array, start_cell_array] = calculate_discharge_at_dune_gaps(
                dune_flow_type=dune_flow_dynamics,
                n_dunes=n_dune_rows,
                elev_array=elevation,
                time_step=TS,
                int_width=int_width,
                bayhigh=bayhigh,
                discharge_array=Discharge,
                bay_array=underwater_array,
                downhill_array=downhill_array,
                start_cell_array=start_cell_array,
            )

            # we check to see if the dune gaps line up through the entire dune line, if not, we do not
            # consider it a dune gap
            if dune_flow_dynamics.lower() == "full":  # we will not route flow through the dunes
                for col in range(length):
                    # we check all the rows in each column and if there is a discharge value = 0
                    # in any column position, we set the entire column to 0
                    if np.any(Discharge[TS, int_width:(int_width + n_dune_rows), col] == 0):
                        Discharge[TS, int_width:(int_width + n_dune_rows), col] = 0
                        # we also update our arrays as they will no longer be submerged
                        underwater_array[TS, int_width:(int_width + n_dune_rows), col] = 0

                        # should these have been set earlier? TBD
                        downhill_array[TS, int_width:(int_width + n_dune_rows), col] = 0
                        start_cell_array[TS, int_width:(int_width + n_dune_rows), col] = 0

            # set the first row of dunes as downhill
            for col in range(length):
                if Discharge[TS, int_width, col] > 0:
                    downhill_array[TS, int_width, col] = 2
                    start_cell_array[TS, int_width, col] = 2

            # check to see if the bay level is high enough to reach the dune gaps through the entire
            # domain (check for hydrualic connection from the bay (first row) to the dune line)
            bay_sufficient, underwater_array, downhill_array, start_cell_array = check_underwater(
                start_row=int_width,
                time_step=TS,
                elev_array=elevation,
                bayhigh=bayhigh,
                length=length,
                bay_array=underwater_array,
                downhill_array=downhill_array,
                start_cell_array=start_cell_array,
            )

            if bay_sufficient is True:
                OW_TS.append(TS)
                start_row = int_width

                # remove any submerged cells not connected to the bay cells
                underwater_array, downhill_array, start_cell_array = clear_floating_submerged_cells(
                    underwater_array=underwater_array,
                    downhill_array=downhill_array,
                    start_cell_array=start_cell_array,
                    width=width,
                    time_step=TS,
                    length=length
                )

                # since we assume that bed slopes drive water slopes, we want to focus the flow on the downhill cells
                # look for downhill cells connected to the dune gaps. then, we select the most landward cells not
                # influenced by any adjacent cells as the start cells and initalize discharge
                Discharge, downhill_array, start_cell_array, init_discharge_array = initialize_discharge_at_start_cells(
                    start_row=start_row,
                    time_step=TS,
                    elev_array=elevation,
                    discharge=Discharge,
                    s1_array=s1_array,
                    s2_array=s2_array,
                    s3_array=s3_array,
                    bayhigh=bayhigh,
                    length=length,
                    bay_array=downhill_array,
                    start_cell_array=start_cell_array,
                    int_width=int_width,
                    init_dis=init_discharge_array,
                )

            else:
                # if there is not a hydraulic connection from the bay to the dune gaps, set discharge to 0
                Discharge[TS] = 0

            max_dune = Dmaxel - BermEl  # [dam MHW]

            # find min row that has nonzero discharge and start flow routing!
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

    return (elevation, qs_lost_total, Discharge, elev_change_array, underwater_array, downhill_array, start_cell_array,
            init_discharge_array, OW_TS)


class Outwasher:
    """Start an outwash event

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
            percent_washout_to_shoreface=100,
            outwash_beach_file=None,
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
        self._substep = 100
        self._max_slope = -0.25
        self._ki = 8.75E-3
        self._k_lat = 0.5
        self._C = 0.0134
        self._mm = 2
        self._sea_level = sea_level  # equal to 0 dam
        self._bay_depth = bay_depth  # [dam MHW] Depth of bay behind island segment, currently set to 0.3
        self._dune_flow_dynamics = "full"

        # setting up dune domain using b3d
        self._dune_domain = dune_domain  # [dam above berm elevation]
        self._dune_crest = self._dune_domain.max(axis=1)  # dune_crest used to be DuneDomainCrest

        # initializing our barrier interior
        self._interior_domain = interior_domain  # [dam MHW]

        # loading the bay levels for each time step
        self._outwash_storms = np.load(datadir + outwash_storms_file, allow_pickle=True)
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
        self._outwash_TS = np.zeros(time_step_count, dtype=object)  # array for storing outwash volume (m^3)
        self._outwash_flux_TS = np.zeros(time_step_count, dtype=object)  # array for storing outwash flux (m^3/m)
        self._initial_full_domain = []
        self._full_dunes = []
        self._full_domain = []
        self._Qs_shoreface = np.zeros(time_step_count)  # dam^3
        self._Qs_shoreface_per_length = np.zeros(time_step_count)  # dam^3/dam
        self._elevation_change = np.zeros(time_step_count, dtype=object)
        self._post_outwash_beach_domain = np.zeros(time_step_count, dtype=object)

    def update(
            self,
            b3d
    ):

        self._time_index = b3d.time_index
        OW_TS = []

        # reduce beach width by the amount of post-storm shoreline change; if the
        # beach width reaches zero, turn dune migration in B3D back on -- otherwise
        # keep it off (we don't want the dune line to prograde
        # because we have fake houses there!)
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
            # allow for substeps to better simulate morphodynamics
            # (run the same bay level 100 times with proportionally reduced sed transport)
            if self._substep != 1:
                updated_bay_levels = repeat_hydrograph_for_substep(storm_series, substep=self._substep)
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

            (Elevation, qs_lost_total, Discharge, elev_change_array, underwater_array, downhill_array,
             start_cell_array, init_discharge_array, OW_TS) = flow_routing_with_lateral_transport(
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
                OW_TS=OW_TS,
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
            post_outwash_beach_domain = Elevation[-1, int_width + n_dune_rows:, :]
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

            print(" end outwash storm \n", end="")
