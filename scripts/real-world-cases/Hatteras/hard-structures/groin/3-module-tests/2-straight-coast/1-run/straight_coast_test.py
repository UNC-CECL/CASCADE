#!/usr/bin/env python3
"""
Four groin modules on a straight coast at a range of orientations, in BRIE's alongshore solve alone.

    python straight_coast_test.py [--years 50]   ->  <module>/tables, <module>/figures, comparison/

An emulator of BRIE's implicit shoreline-diffusion step (its own coast_diff table, sparse
indices, clip at zero) on an infinitely long straight coast at angle theta0: the shoreline is
the tilted line plus a periodic perturbation eta, so the tilt never meets the periodic wrap.
One groin sits at the face between the middle two domains. The modules:

    1 fixed-shift groin    (the source/sink dipole) GroinCallback: -M updrift, +M downdrift, whatever the coast does
    2 trapping, pinned     BlockingGroinCallback: cancels blocking_fraction of each cell's own diffusive coupling
    3 trapping, conserving the same, one face coefficient on both sides
    4 trapping, drift      cancels blocking_fraction of the wave-climate net drift Q_net across the face
    5 trapping, hybrid     module 4 plus module 2's row correction: blocking_fraction of BRIE's r[HI] - r[LO] mismatch at the face
    6 wave-driven groin    (trapping, corrected) module 4 plus a leak correction: near the groin, gives back what BRIE's row-scaled
                           solve creates or destroys, 2 (r[i] - r[i-1]) (x[i] - x[i-1]) per cell, BRIE untouched
    7 corrected, length    module 6 with a groin length L: blocking_fraction (1 - a/L)^(2/3), a the updrift
                           cell's build-out, so the groin passes everything once the fillet reaches its tip

Q_net is the wave pdf convolved with BRIE's own _coast_qs, aligned as BRIE aligns coast_diff;
coast_diff = -(1/depth) dQ_net/dtheta to within one angle bin, so module 4 uses the same physics
BRIE does. No storms, no Barrier3D, no sea level. Details: README.md beside this script.

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-10-08
"""
from __future__ import annotations

# Non-code files beside this script live in the workspace on D: (hatteras_paths.workspace_twin)
import sys as _hat_sys  # noqa: E402
from pathlib import Path as _HatPath  # noqa: E402
_HAT_DIR = next(_p for _p in _HatPath(__file__).resolve().parents if (_p / "hatteras_paths.py").exists())
if str(_HAT_DIR) not in _hat_sys.path:
    _hat_sys.path.insert(0, str(_HAT_DIR))
from hatteras_paths import workspace_twin  # noqa: E402
import argparse
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402
from scipy.sparse import csr_matrix  # noqa: E402
from scipy.sparse.linalg import spsolve  # noqa: E402

from brie import Brie  # noqa: E402

HERE = Path(__file__).resolve().parent
# Hatteras/ holds the scripts; data/ and output/ live in the workspace on D: (hatteras_paths)
import sys as _hat_sys  # noqa: E402
_HAT_DIR = next(p for p in HERE.parents  if (p / "hatteras_paths.py").exists())
if str(_HAT_DIR) not in _hat_sys.path:
    _hat_sys.path.insert(0, str(_HAT_DIR))
from hatteras_paths import WORKSPACE as REPO  # noqa: E402
sys.path[:0] = [str(_HAT_DIR)]
from site_layer.hat_figure_style import (  # noqa: E402
    C, INK, INK_MUTED, _letter_inside, apply_style, figsize, open_frame, record_caption, save)

# --- CONFIG ------------------------------------------------------------------
DY_M, DT_YR = 500.0, 1.0                 # fixed by the CASCADE coupler
CLIMATE = dict(Hs=2.0, Tp=7.5, asym=0.6, ahf=0.5)   # option A
NY = 61
LO, HI = 30, 31                          # the groin is the face between these two
ORIENTATIONS = tuple(range(-40, 45, 5))  # shoreline angle theta0, degrees
SHOWN = (-30, -10, 0, 10, 20, 30)        # orientations drawn as profiles
PROFILE_YEARS = (1, 5, 10, 25, 50)
BLOCKING_FRACTION = 0.6                  # share of the transport the groin holds; the 10-05 pin
                                         # (cascade.groin.BlockingGroinCallback's blocking_fraction)
M = 12.0                                 # fixed-shift groin, the best option-A dipole (2026-09-29)
# Versions 2-5 and 7 were retired on 2026-10-12; their results live in ../../archive/2-straight-coast/
MODULES = {
    "dipole": dict(folder="fixed_shift", label="Fixed-shift", col=C["ADDED"]),
    "pinned": dict(folder="2-trapping-pinned", archived=True, label="Trapping, pinned", col=C["ACCENT"]),
    "conserving": dict(folder="3-trapping-conserving", archived=True, label="Trapping, conserving", col=C["REF"]),
    "drift": dict(folder="4-trapping-drift", archived=True, label="Trapping, drift", col=C["LATE"]),
    "hybrid": dict(folder="5-trapping-hybrid", archived=True, label="Trapping, hybrid", col=C["LOCATOR"]),
    "corrected": dict(folder="wave_driven", label="Wave-driven", col=INK),
    "length": dict(folder="7-corrected-length", archived=True, label="Corrected, with length", col="#e7298a"),
}
CORRECTION_HALF_WIDTH = 5                # module 6: cells either side of the groin face that get the leak correction
GROIN_LENGTH_M = 110.0                   # module 7: Buxton's digitized groins are 90-155 m
ARCHIVE = workspace_twin(HERE.parents[1]) / "archive" / "2-straight-coast"
# -----------------------------------------------------------------------------


# Where a version's tables and figures live: here (tables/, figures/), or its own folder in the archive
def module_dir(spec):
    return ARCHIVE / spec["folder"] if spec.get("archived") else HERE


# coast_diff, Q_net, depth and BRIE's sparse indices from a real Brie instance
def brie_tables():
    b = Brie(ast_model=True, barrier_model=False, inlet_model=False, b3d=True,
             wave_height=CLIMATE["Hs"], wave_period=CLIMATE["Tp"],
             wave_asymmetry=CLIMATE["asym"], wave_angle_high_fraction=CLIMATE["ahf"],
             alongshore_section_length=DY_M, alongshore_section_count=NY,
             time_step=DT_YR, time_step_count=10, save_spacing=1)
    step = 180.0 / (b._wave_climl - 1)
    pdf = b._angles.pdf(b._angle_array) * np.deg2rad(step)
    conv = np.convolve(pdf, b._coast_qs, mode="full")
    npad = len(b._coast_qs) - 1
    first = npad - npad // 2
    g = conv[first:first + len(pdf)]
    return dict(coast_diff=b._coast_diff.copy(), g=g, depth=b._h_b_crit + b._d_sf,
                di=b._di.copy(), dj=b._dj.copy(), climl=b._wave_climl)


def angle_index(theta_deg, climl):
    return np.clip(np.round(90 - theta_deg).astype(int), 1, climl - 1)


# Net drift through a face toward higher domain index, m3/yr (minus BRIE's convolved G)
def drift(t, theta_deg):
    return -t["g"][angle_index(np.atleast_1d(theta_deg), t["climl"])]


# What BRIE's row-scaled solve adds beyond a face-centred one, with the sign reversed, near the groin.
# Row i couples to both neighbours with r[i]; a conserving solve uses r[i-1] on its lower face, so the
# difference is (r[i] - r[i-1]) times the lower step; x2 for the two Crank-Nicolson halves.
def leak_correction(r, step, half_width=None):
    half_width = CORRECTION_HALF_WIDTH if half_width is None else half_width
    near = np.abs(np.arange(r.size) - (LO + 0.5)) <= half_width
    return np.where(near, 2 * (r - np.roll(r, 1)) * np.roll(step, 1), 0.0)


# Module 7's share of blocking_fraction still held: GENESIS bypassing, falling to zero as the fillet reaches the tip
def length_factor(advance_m, length_m=GROIN_LENGTH_M):
    return max(0.0, 1.0 - max(0.0, advance_m) / length_m) ** (2.0 / 3.0)


def run(t, theta0, module, years, blocking_fraction=BLOCKING_FRACTION, correction_fraction=BLOCKING_FRACTION):
    slope = np.tan(np.deg2rad(theta0))
    eta = np.zeros(NY)
    f0 = float(drift(t, theta0)[0])
    up, down = (LO, HI) if f0 > 0 else (HI, LO)       # drift toward higher index: LO is updrift
    rows, profiles = [], {}
    applied_cum = 0.0
    for year in range(1, years + 1):
        step = slope * DY_M + (np.roll(eta, -1) - eta)          # x[i+1] - x[i], landward +
        theta = 180.0 * np.arctan2(step, DY_M) / np.pi
        r = np.maximum(0.0, t["coast_diff"][angle_index(theta, t["climl"])] * DT_YR / 2 / DY_M ** 2)
        off = step[LO]
        dx = np.zeros(NY)
        if module == "dipole":
            dx[up] -= M
            dx[down] += M
        elif module == "pinned":
            dx[LO] -= blocking_fraction * 2 * r[LO] * off
            dx[HI] += blocking_fraction * 2 * r[HI] * off
        elif module == "conserving":
            dx[LO] -= blocking_fraction * 2 * r[LO] * off
            dx[HI] += blocking_fraction * 2 * r[LO] * off
        elif module in ("drift", "hybrid", "corrected", "length"):
            # the drift arriving at the groin: the angle of the updrift approach face, not of the step
            approach = theta[HI] if up == HI else theta[LO - 1]
            q = blocking_fraction * float(drift(t, approach)[0]) * DT_YR / (DY_M * t["depth"])
            if module == "length":
                q *= length_factor(-eta[up])
            dx[LO] -= q
            dx[HI] += q
            if module == "hybrid":
                # the pinned module's row correction: BRIE counts the face with r[HI] in HI's row, r[LO] in LO's
                dx[HI] += correction_fraction * 2 * (r[HI] - r[LO]) * off
            if module in ("corrected", "length"):
                dx += leak_correction(r, step)
        rhs = eta + r * (np.roll(eta, -1) - 2 * eta + np.roll(eta, 1)) + dx
        dv = np.r_[-r[-1], -r[1:], 1 + 2 * r, -r[0:-1], -r[0]]
        eta = spsolve(csr_matrix((dv, (t["di"], t["dj"]))), rhs)
        applied_cum += dx.sum()
        rows.append(dict(
            theta0_deg=theta0, year=year, drift0_m3_yr=f0,
            updrift_seaward_m=-eta[up], downdrift_seaward_m=-eta[down],
            applied_up_m=-dx[up], applied_down_m=-dx[down],
            applied_net_cum_m=-applied_cum,                     # seaward +: sand the module added
            reach_change_m=-eta.sum(),                          # sum over domains, seaward +
            solve_error_m=-(eta.sum() - applied_cum),           # what the solve itself added
            face_angle_deg=float(theta[LO]), r_lo=float(r[LO]), r_hi=float(r[HI])))
        if year in PROFILE_YEARS:
            prof = -eta.copy()                                  # seaward +, updrift at positive offsets
            profiles[year] = prof if up == HI else prof[(LO + HI - np.arange(NY)) % NY]
    return pd.DataFrame(rows), profiles


# Pelnard-Considere fillet at a fully blocking groin: 2 tan(a) sqrt(D t / pi), tan(a) = Q0 / (D depth)
def pelnard_considere(t, theta0, years):
    d = float(t["coast_diff"][angle_index(np.array([theta0]), t["climl"])][0])
    q0 = abs(float(drift(t, theta0)[0]))
    if d <= 0:
        return np.nan, np.nan
    tan_a = q0 / (d * t["depth"])
    return 2 * tan_a * np.sqrt(d * years / np.pi), tan_a


YEAR_RAMP = plt.get_cmap("Blues")(np.linspace(0.35, 0.95, len(PROFILE_YEARS)))   # light early, dark late


# The letter inside the corner, a short title centred above
def label(ax, i, text):
    _letter_inside(ax, i)
    ax.set_title(text, loc="center")


# The drift as an arrow along the top of a panel, from the updrift (right) side toward the groin
def drift_arrow(ax, q):
    ax.annotate("", xy=(0.56, 0.93), xytext=(0.9, 0.93), xycoords="axes fraction",
                arrowprops=dict(arrowstyle="-|>", color=INK_MUTED, lw=1.0))
    ax.text(0.73, 0.86, f"drift {abs(q) / 1e6:.2f} M m$^3$/yr", transform=ax.transAxes, ha="center", va="top",
            fontsize=6.5, color=INK_MUTED)


def module_figures(name, frame, profiles):
    spec = MODULES[name]
    out = module_dir(spec) / "figures"
    out.mkdir(parents=True, exist_ok=True)
    dom = np.arange(NY) - LO - 0.5                              # domains from the groin face
    fig, axes = plt.subplots(2, 3, figsize=figsize("double", height=5.4), sharex=True,
                             constrained_layout=True)
    lim = max(np.abs(profiles[(th, y)]).max() for th in SHOWN for y in PROFILE_YEARS)
    lim = max(5.0, np.ceil(lim / 10) * 10) * 1.15
    for i, (ax, th) in enumerate(zip(axes.flat, SHOWN)):
        for c, y in zip(YEAR_RAMP, PROFILE_YEARS):
            ax.plot(dom, profiles[(th, y)], color=c, lw=1.3)
        ax.axvline(0, color=C["GROIN"], lw=0.9)
        ax.axhline(0, color=INK_MUTED, lw=0.5)
        ax.set_xlim(-12, 12)
        ax.set_ylim(-lim, lim)
        drift_arrow(ax, frame[(frame.theta0_deg == th) & (frame.year == 1)].drift0_m3_yr.iloc[0])
        ax.grid(axis="y")
        open_frame(ax)
        label(ax, i, f"Coast at {th:+d}$^\\circ$")
    for ax in axes[:, 0]:
        ax.set_ylabel("Shoreline change\n(m, seaward +)")
    for ax in axes[1]:
        ax.set_xlabel("Domains from the groin (updrift right)")
    fig.suptitle(f"Does the {spec['label'].lower()} groin hold sand updrift on a straight coast?")
    handles = [Line2D([], [], color=c, lw=1.3, label=f"year {y}") for c, y in zip(YEAR_RAMP, PROFILE_YEARS)]
    handles.append(Line2D([], [], color=C["GROIN"], lw=0.9, label="groin"))
    fig.legend(handles=handles, loc="outside lower center", ncol=6, fontsize=7, frameon=False)
    png = workspace_twin(out) / f"straight_coast_{spec['folder']}_profiles.png"
    save(fig, png, dpi=250, close=True)
    last = frame[frame.year == frame.year.max()].set_index("theta0_deg")
    record_caption(png, (
        f"**Tests:** what the {spec['label'].lower()} groin does on its own, on a straight coast with nothing "
        "else going on, and whether it behaves like a groin: sand held on the updrift side, taken from the "
        "downdrift side, and no sand created or destroyed. **How to read:** BRIE's alongshore solve alone "
        "(option A waves; no storms, no Barrier3D, no sea-level rise), one groin between the two middle "
        "domains. Each panel is one shoreline orientation; lines are the shoreline change from the straight "
        "start at years 1, 5, 10, 25 and 50 (light to dark), seaward positive, with the updrift side on the "
        "right; the arrow is the wave-climate net longshore drift on the undisturbed coast, which runs toward "
        "the groin from the right. A groin that traps drift builds out on the right and cuts back on the "
        "left. Red line: the groin face. **Shows:** after 50 years the cell just updrift has moved "
        f"{last.loc[-30, 'updrift_seaward_m']:+.0f} m at -30$^\\circ$, {last.loc[0, 'updrift_seaward_m']:+.0f} m at "
        f"0$^\\circ$ and {last.loc[30, 'updrift_seaward_m']:+.0f} m at +30$^\\circ$; the module's own sand "
        f"balance over 50 years ranges from {last['applied_net_cum_m'].min():+.0f} to "
        f"{last['applied_net_cum_m'].max():+.0f} m of summed shoreline across orientations, and the total "
        f"change in the reach (module plus solve) from {last['reach_change_m'].min():+.0f} to "
        f"{last['reach_change_m'].max():+.0f} m (comparison figures, and tables/)."))


def comparison_figures(t, summary, frames):
    out = workspace_twin(ARCHIVE) / "comparison" / "figures"   # all seven versions: archived
    out.mkdir(parents=True, exist_ok=True)
    th = np.arange(-45, 46)
    dif = t["coast_diff"][angle_index(th, t["climl"])]
    q = drift(t, th)
    k = np.flatnonzero(np.diff(np.sign(q)) != 0)
    zero_q = float(np.interp(0, q[k[0]:k[0] + 2][::-1], th[k[0]:k[0] + 2][::-1])) if len(k) else np.nan
    fig, axes = plt.subplots(2, 3, figsize=figsize("double", height=5.8), constrained_layout=True)
    ax = axes[0, 0]
    ax.plot(th, -q / 1e6, color=INK, lw=1.5, label="net drift toward the groin")
    ax.axhline(0, color=INK_MUTED, lw=0.5)
    ax.set_ylabel("Net drift (M m$^3$/yr)")
    ax2 = ax.twinx()
    ax2.plot(th, dif / 1e3, color=C["BASE"], lw=1.2, ls="--", label="BRIE diffusivity")
    frozen = th[dif <= 0]
    if len(frozen):
        ax.axvspan(frozen.min() - 0.5, frozen.max() + 0.5, color=C["BASE_FILL"], alpha=0.8, lw=0)
    ax2.set_ylabel("BRIE diffusivity\n(10$^3$ m$^2$/yr, dashed)", color=C["BASE"], fontsize=7)
    ax2.tick_params(labelsize=6.5)
    if np.isfinite(zero_q):
        ax.axvline(zero_q, color=INK_MUTED, lw=0.6, ls=":")
        ax.text(zero_q - 1.5, ax.get_ylim()[0] + 0.04 * np.ptp(ax.get_ylim()), f"zero drift\n{zero_q:+.0f}$^\\circ$",
                ha="right", va="bottom",
                fontsize=6.3, color=INK_MUTED)
    ax.set_xlabel("Shoreline orientation ($^\\circ$)")
    open_frame(ax)
    label(ax, 0, "What drives the groin")
    panels = (("updrift_seaward_m", "Updrift cell after 50 yr\n(m, seaward +)", "Sand held updrift"),
              ("applied_net_cum_m", "Sand the module adds, 50 yr\n(m of summed shoreline)", "The module's own balance"),
              ("solve_error_m", "Sand the solve adds, 50 yr\n(m of summed shoreline)", "BRIE's solve"),
              ("reach_change_m", "Total change in the reach, 50 yr\n(m of summed shoreline)", "Module plus solve"))
    for k, (col, ylab, title) in enumerate(panels):
        ax = axes.flat[k + 1]
        for name, spec in MODULES.items():
            s_ = summary[summary.module == name]
            ax.plot(s_.theta0_deg, s_[col], color=spec["col"], lw=1.6, marker="o", ms=3)
        ax.axhline(0, color=INK_MUTED, lw=0.6)
        ax.axvline(0, color=INK_MUTED, lw=0.5, ls=(0, (2, 2)))
        ax.set_xlabel("Shoreline orientation ($^\\circ$)")
        ax.set_ylabel(ylab, fontsize=7)
        ax.grid(axis="y")
        open_frame(ax)
        label(ax, k + 1, title)
    ax = axes.flat[5]
    ax.axis("off")
    handles = [Line2D([], [], color=s_["col"], lw=1.6, marker="o", ms=3, label=s_["label"]) for s_ in MODULES.values()]
    handles += [Line2D([], [], color=INK, lw=1.5, label="net drift toward the groin (a)"),
                Line2D([], [], color=C["BASE"], lw=1.2, ls="--", label="BRIE diffusivity (a)"),
                plt.Rectangle((0, 0), 1, 1, color=C["BASE_FILL"], alpha=0.8, label="diffusivity at or below zero:\ncoast frozen (a)")]
    ax.legend(handles=handles, loc="center", fontsize=7.5, frameon=False)
    fig.suptitle("Does each version hold sand updrift, and is sand conserved?")
    png = workspace_twin(out) / "straight_coast_modules_vs_orientation.png"
    save(fig, png, dpi=250, close=True)
    sm = summary.set_index(["module", "theta0_deg"])
    tot = summary.groupby("module")["reach_change_m"].agg(lambda v: v.abs().max())
    record_caption(png, (
        "**Tests:** across shoreline orientations from -40 to +40 degrees, whether each version of the groin "
        "module holds sand updrift, and whether the sand balance closes, for the module on its own, for BRIE's "
        "alongshore solve, and for the two together. **How to read:** 50 years of BRIE's alongshore solve "
        "alone (option A waves), one groin, every value read at year 50. (a) What drives the modules: the "
        "wave-climate net longshore drift toward the groin (solid) and BRIE's shoreline diffusivity (dashed); "
        "where the diffusivity is at or below zero BRIE clips it to zero (shaded), which freezes a straight coast. "
        "(b) How far the cell just updrift of the groin has built out. (c) The sand the module itself adds, as "
        "shoreline change summed over the reach; zero means the module's own books balance. (d) The sand "
        "BRIE's row-scaled solve adds on top, non-zero wherever the groin makes the diffusivity change along "
        "the coast. (e) The total, (c) plus (d): zero means the groin, in the model, neither creates nor "
        f"destroys sand. Dipole M {M:g} m/yr; trapping versions blocking fraction {BLOCKING_FRACTION}. **Shows:** the net drift is "
        f"{-drift(t, 0)[0] / 1e6:.2f} M m$^3$/yr toward the groin on a grid-aligned coast and reaches zero "
        f"near {zero_q:+.0f} degrees, while the diffusivity is about zero there, so the pinned and conserving "
        "versions, which block only the diffusive part of the transport, do nothing at 0 degrees. Only the "
        f"pinned version stays near zero in (e) (largest |total| {tot['pinned']:.0f} m against "
        f"{tot['conserving']:.0f} m conserving, {tot['dipole']:.0f} m dipole, {tot['drift']:.0f} m drift, {tot['hybrid']:.0f} m hybrid and {tot['corrected']:.0f} m corrected), "
        "because its own imbalance offsets the solve's; the drift version traps the right quantity but sits on "
        "a solve that does not conserve sand, and the hybrid, which adds the pinned version's row correction "
        "at the groin face to the drift version, does not close the balance either, because the solve loses "
        "sand along the whole bent fillet, not only at the groin face. The corrected version gives that sand back "
        "within five cells of the groin each year, and stays near zero while trapping on the updrift side."))

    # Year-50 profiles, the versions side by side at the shown orientations
    fig, axes = plt.subplots(2, 3, figsize=figsize("double", height=5.4), sharex=True,
                             constrained_layout=True)
    dom = np.arange(NY) - LO - 0.5
    f0 = frames["dipole"][0]
    for i, (ax, th0) in enumerate(zip(axes.flat, SHOWN)):
        for name, spec in MODULES.items():
            ax.plot(dom, frames[name][1][(th0, 50)], color=spec["col"], lw=1.5)
        ax.axvline(0, color=C["GROIN"], lw=0.9)
        ax.axhline(0, color=INK_MUTED, lw=0.5)
        ax.set_xlim(-12, 12)
        lo_, hi_ = ax.get_ylim()
        ax.set_ylim(lo_, hi_ + 0.18 * (hi_ - lo_))
        drift_arrow(ax, f0[(f0.theta0_deg == th0) & (f0.year == 1)].drift0_m3_yr.iloc[0])
        ax.grid(axis="y")
        open_frame(ax)
        label(ax, i, f"Coast at {th0:+d}$^\\circ$")
    for ax in axes[:, 0]:
        ax.set_ylabel("Shoreline change after\n50 yr (m, seaward +)")
    for ax in axes[1]:
        ax.set_xlabel("Domains from the groin (updrift right)")
    fig.suptitle("Where does each version put the sand after 50 years?")
    fig.legend(handles=handles[:len(MODULES)] + [Line2D([], [], color=C["GROIN"], lw=0.9, label="groin")],
               loc="outside lower center", ncol=4, fontsize=7, frameon=False)
    png = workspace_twin(out) / "straight_coast_profiles_year50.png"
    save(fig, png, dpi=250, close=True)
    record_caption(png, (
        "**Tests:** side by side, where each version of the groin module puts the sand after 50 years, at six "
        "shoreline orientations. **How to read:** BRIE's alongshore solve alone (option A waves), one groin; "
        "shoreline change from the straight start, seaward positive, with the updrift side on the right; the "
        "arrow is the net longshore drift toward the groin. A groin that traps drift builds out on the right "
        "and cuts back on the left. Each panel has its own vertical range. **Shows:** on coasts turned "
        "against the drift (negative angles) every trapping version builds out updrift; at 0 degrees only the "
        "drift version and the fixed-shift groin change the coast; on coasts turned with the drift (positive angles) the "
        "pinned version builds out on the wrong (downdrift) side, and the conserving and drift versions cut "
        "back both sides, the sand lost through BRIE's solve "
        f"(summary table: comparison/tables/straight_coast_summary_year50.csv; e.g. at +20 degrees the "
        f"conserving version's updrift cell is at {sm.loc[('conserving', 20), 'updrift_seaward_m']:+.0f} m)."))


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--years", type=int, default=50)
    a = ap.parse_args()
    apply_style()
    t = brie_tables()
    print(f"depth {t['depth']:.2f} m; drift at 0 deg {drift(t, 0)[0]:+.0f} m3/yr")

    frames, summary = {}, []
    for name, spec in MODULES.items():
        rows, profiles = [], {}
        for th in ORIENTATIONS:
            f, p = run(t, th, name, a.years)
            rows.append(f)
            profiles.update({(th, y): v for y, v in p.items()})
        frame = pd.concat(rows, ignore_index=True)
        tables = module_dir(spec) / "tables"
        tables.mkdir(parents=True, exist_ok=True)
        frame.to_csv(tables / f"straight_coast_{spec['folder']}_by_year.csv", index=False,
                     float_format="%.4f")
        pd.DataFrame([dict(theta0_deg=th, year=y, domain_from_groin=int(i - LO), seaward_m=v)
                      for (th, y), prof in profiles.items() for i, v in enumerate(prof)]).to_csv(
            tables / f"straight_coast_{spec['folder']}_profiles.csv", index=False, float_format="%.4f")
        module_figures(name, frame, profiles)
        frames[name] = (frame, profiles)
        last = frame[frame.year == a.years].copy()
        last.insert(0, "module", name)
        summary.append(last)

    summary = pd.concat(summary, ignore_index=True)
    bench = []
    for th in ORIENTATIONS:
        f, _ = run(t, th, "drift", a.years, blocking_fraction=1.0)
        y_pc, tan_a = pelnard_considere(t, th, a.years)
        bench.append(dict(theta0_deg=th, model_updrift_m=float(f.updrift_seaward_m.iloc[-1]),
                          pelnard_considere_m=y_pc, tan_alpha=tan_a,
                          linear_ok=bool(tan_a < 0.2) if np.isfinite(tan_a) else False))
    (workspace_twin(HERE) / "tables").mkdir(parents=True, exist_ok=True)
    summary.to_csv(workspace_twin(HERE) / "tables" / "straight_coast_summary_year50.csv",
                   index=False, float_format="%.3f")
    pd.DataFrame(bench).to_csv(workspace_twin(HERE) / "tables" / "pelnard_considere_check.csv",
                               index=False, float_format="%.3f")
    comparison_figures(t, summary, frames)

    pd.set_option("display.width", 200)
    print(summary.pivot(index="theta0_deg", columns="module",
                        values="updrift_seaward_m").round(1).to_string())
    print(summary.pivot(index="theta0_deg", columns="module",
                        values="applied_net_cum_m").round(1).to_string())
    print(summary.pivot(index="theta0_deg", columns="module",
                        values="solve_error_m").round(1).to_string())
    print(pd.DataFrame(bench).round(2).to_string(index=False))


if __name__ == "__main__":
    main()
