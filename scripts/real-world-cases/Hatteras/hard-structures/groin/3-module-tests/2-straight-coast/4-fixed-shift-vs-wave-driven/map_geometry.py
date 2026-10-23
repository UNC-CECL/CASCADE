#!/usr/bin/env python3
"""
The straight-coast maps' geometry: each 500 m cell drawn as a rectangle along the coast, moved shore-normal.

    from map_geometry import shoreline, draw_coast

Shared by straight_coast_fixed_shift_vs_wave_driven_map_figure.py and _map_gif.py; draws nothing on its own.
BRIE keeps one cross-shore position per cell and moves the whole cell, so a cell is a rectangle with its
long side on the coast. Here cell k (offset from the groin, + north) spans alongshore distance (k - 1) to k
cells from the groin face along the start coast, and its change moves it along the coast's seaward normal,
so the step between two cells is square to the coast. Before 2026-10-16 the maps moved each cell east-west,
which drew the steps horizontal and the cells as parallelograms on a turned coast.

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-10-16
"""
from __future__ import annotations

import numpy as np
from matplotlib.patches import Polygon

REACH_M = 6000.0                     # how far land and sea are filled either side of the shoreline


# Unit vectors of a coast at theta: alongshore going north, and the seaward (east-facing) normal
def coast_axes(theta):
    r = np.deg2rad(theta)
    return np.array([-np.sin(r), np.cos(r)]), np.array([np.cos(r), np.sin(r)])


# A shoreline as a polyline: each cell's straight edge, set seaward of the start coast by its change
def shoreline(theta, offs, seaward, dy):
    along, normal = coast_axes(theta)
    pts = [s * along + e * normal for k, e in zip(offs, seaward) for s in ((k - 1) * dy, k * dy)]
    pts = np.array(pts)
    return pts[:, 0], pts[:, 1]


# Land behind the shoreline, sea in front, and each cell's change since the start as a rectangle
def draw_coast(ax, theta, offs, seaward, dy, land, sea, built, cut):
    along, normal = coast_axes(theta)
    line = np.c_[shoreline(theta, offs, seaward, dy)]
    ax.add_patch(Polygon(np.r_[line, (line - REACH_M * normal)[::-1]], closed=True, lw=0, zorder=0, **land))
    ax.add_patch(Polygon(np.r_[line, (line + REACH_M * normal)[::-1]], closed=True, lw=0, zorder=0, **sea))
    for k, e in zip(offs, seaward):
        if e:
            a, b = (k - 1) * dy * along, k * dy * along
            ax.add_patch(Polygon([a, b, b + e * normal, a + e * normal], closed=True, lw=0, zorder=1,
                                 **(built if e > 0 else cut)))
