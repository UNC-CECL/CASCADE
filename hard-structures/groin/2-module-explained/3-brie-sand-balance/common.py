#!/usr/bin/env python3
"""
Shared by the five question folders here: an unmodified BRIE with only its alongshore step on, the same step with one
rate per boundary, and the drawing helpers. Imported, not run.

Every run is a real `brie.Brie` with the barrier and inlet models off, so `update()` sets x_s_dt to zero and the only
thing that moves the shoreline is BRIE's alongshore step (brie.py 1293-1328), unchanged. Option A waves (Hs 2.0 m,
Tp 7.5 s, asymmetry 0.6, high-angle fraction 0.5), 30 sections of 500 m on BRIE's periodic coast, 1-year steps.
"One rate per boundary" re-solves BRIE's step with each boundary's rate shared by the two sections it joins (BRIE's own
rates, its Crank-Nicolson weighting, its periodic coast); it conserves to round-off. BRIE is not edited.

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-10-15
"""
from __future__ import annotations

import sys
import warnings
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402,F401
import numpy as np  # noqa: E402
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch, Rectangle  # noqa: E402,F401
from scipy.sparse import diags  # noqa: E402
from scipy.sparse.linalg import spsolve  # noqa: E402

REPO = next(p for p in Path(__file__).resolve().parents if (p / "pyproject.toml").exists())
sys.path.insert(0, str(REPO / "scripts"))
from brie import Brie  # noqa: E402
from site_layer.hat_figure_style import C, INK, INK_MUTED  # noqa: E402

warnings.filterwarnings("ignore", category=RuntimeWarning)

# --- CONFIG ------------------------------------------------------------------
NY, DY, YEARS = 30, 500.0, 50
WAVES = dict(wave_height=2.0, wave_period=7.5, wave_asymmetry=0.6, wave_angle_high_fraction=0.5)   # option A
I = np.arange(NY)
SHAPES = {                                   # seaward + (m)
    "straight": np.zeros(NY),
    "gentle bump": 100 * np.exp(-0.5 * ((I - 15) / 3.0) ** 2),
    "sharp step": np.where((I >= 15) & (I < 18), 100.0, 0.0),
}
SHAPE_C = {"straight": INK_MUTED, "gentle bump": C["LATE"], "sharp step": C["GROIN"]}
SHARED_C = C["REF"]
SAND, SEA, SAND_EDGE = "#f3dca6", "#cfe3f1", "#b8913f"
GOOD, BAD, GAIN = C["REF"], C["GROIN"], C["ADDED"]
TRUCK_M3 = 10.0                              # one dump truck
LO, HI = 11, 21                              # sections drawn around the headland (sections 15-17)
HEADLAND = SHAPES["sharp step"]
# -----------------------------------------------------------------------------


# A real Brie with only its alongshore step on, started from `shape` (seaward +)
def brie_on(shape, years=YEARS):
    b = Brie(name="sand", barrier_model=False, inlet_model=False, ast_model=True, **WAVES,
             alongshore_section_length=DY, alongshore_section_count=NY, time_step=1.0,
             time_step_count=years + 2, sea_level_rise_rate=0.0)
    b.base = float(np.mean(b._x_s))
    b._x_s = b.base - np.asarray(shape, float)          # BRIE's x_s is landward +
    return b


def seaward(b, x=None):
    return -((b._x_s if x is None else x) - b.base)


# BRIE's per-domain rate (brie.py 820-826 for the angle, 1294-1309 for the rate): from the boundary to the north
def rates(b, x):
    theta = 180 * np.arctan2(x[np.r_[1:NY, 0]] - x, b._dy) / np.pi
    k = np.maximum(1, np.minimum(b._wave_climl, np.round(90 - theta).astype(int)))
    return np.maximum(0, b._coast_diff[k] * b._dt / 2 / b._dy ** 2), theta


# One year of BRIE's step with one rate per boundary: boundary i|i+1 uses domain i's rate on both sides
def shared_step(b, x):
    r, _ = rates(b, x)
    rs = np.roll(r, 1)
    L = (diags([rs[1:], r[:-1]], [-1, 1], shape=(NY, NY)) + diags([-(r + rs)], [0])).tolil()
    L[0, NY - 1], L[NY - 1, 0] = rs[0], r[NY - 1]
    L = L.tocsr()
    eye = diags([np.ones(NY)], [0])
    return spsolve((eye - L).tocsc(), x + L @ x)


# Summed shoreline (seaward +) each year: real Brie.update(), and the same start with one rate per boundary
def run(shape, years=YEARS):
    b = brie_on(shape, years)
    x = b._x_s.copy()
    real, shared, prof = [seaward(b).sum()], [seaward(b).sum()], {0: seaward(b).copy()}
    for y in range(1, years + 1):
        b.update()
        x = shared_step(b, x)
        real.append(seaward(b).sum())
        shared.append(seaward(b, x).sum())
        if y in (years,):
            prof[y] = seaward(b).copy()
            prof["shared"] = seaward(b, x).copy()
    return np.array(real) - real[0], np.array(shared) - shared[0], prof, b


def depth(b):
    return b._h_b_crit + b._d_sf                       # active depth that turns metres of shoreline into volume


def box(ax, x, y, w, h, text, fc, ec=INK, fs=7.2, weight="normal"):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.02,rounding_size=0.06", fc=fc, ec=ec, lw=0.9))
    ax.text(x + w / 2, y + h / 2, text, ha="center", va="center", fontsize=fs, color=INK, linespacing=1.35,
            fontweight=weight)


def blank(ax):
    ax.set_xticks([])
    ax.set_yticks([])
    for s in ax.spines.values():
        s.set_visible(False)


def bubble(ax, x, y, text, ec=INK, fc="white", fs=10, ha="center", weight="normal", color=INK):
    ax.text(x, y, text, ha=ha, va="center", fontsize=fs, color=color, fontweight=weight, linespacing=1.3,
            bbox=dict(boxstyle="round,pad=0.45", fc=fc, ec=ec, lw=1.2))


def arrow(ax, x0, y0, x1, y1, col=INK, lw=2.5, style="-|>", ms=18, ls="-"):
    ax.add_patch(FancyArrowPatch((x0, y0), (x1, y1), arrowstyle=style, mutation_scale=ms, color=col, lw=lw,
                                 linestyle=ls, shrinkA=0, shrinkB=0))


# Plan view, north to the right: land below the shoreline, sea above; sections drawn as columns
def coast(ax, s, scale=1.0, lo=LO, hi=HI, outline=INK, base=-60):
    cells = np.arange(lo, hi)
    for i in cells:
        h = s[i] * scale
        ax.add_patch(Rectangle((i - 0.5, base), 1, h - base, fc=SAND, ec="none"))
        ax.add_patch(Rectangle((i - 0.5, h), 1, 400, fc=SEA, ec="none"))
    xs = np.repeat(np.r_[cells - 0.5, cells[-1] + 0.5], 2)[1:-1]
    ax.plot(xs, np.repeat(s[cells] * scale, 2), color=outline, lw=2.2, solid_joinstyle="miter")
    for i in cells[1:]:
        ax.plot([i - 0.5, i - 0.5], [base, min(s[i - 1], s[i]) * scale], color=SAND_EDGE, lw=0.8, ls=(0, (2, 2)))
    ax.set_xlim(lo - 0.5, hi - 0.5)


# The headland in year 0 with BRIE's own numbers
def headland_numbers():
    b = brie_on(HEADLAND, years=1)
    r, theta = rates(b, b._x_s.copy())
    s = seaward(b)
    edges = {}
    for e in (14, 17):                         # boundary e|e+1
        d = s[e + 1] - s[e]
        edges[e] = (2 * r[e] * d + 0.0, -2 * r[e + 1] * d + 0.0)   # south section's book, north section's book
    return b, r, theta, s, edges
