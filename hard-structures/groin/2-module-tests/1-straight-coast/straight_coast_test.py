#!/usr/bin/env python3
"""
Four groin modules on a straight coast at a range of orientations, in BRIE's alongshore solve alone.

    python straight_coast_test.py [--years 50]   ->  <module>/tables, <module>/figures, comparison/

An emulator of BRIE's implicit shoreline-diffusion step (its own coast_diff table, sparse
indices, clip at zero) on an infinitely long straight coast at angle theta0: the shoreline is
the tilted line plus a periodic perturbation eta, so the tilt never meets the periodic wrap.
One groin sits at the face between the middle two domains. The modules:

    1 source/sink dipole   GroinCallback: -M updrift, +M downdrift, whatever the coast does
    2 trapping, pinned     BlockingGroinCallback: cancels b of each cell's own diffusive coupling
    3 trapping, conserving the same, one face coefficient on both sides
    4 trapping, drift      cancels b of the wave-climate net drift Q_net across the face

Q_net is the wave pdf convolved with BRIE's own _coast_qs, aligned as BRIE aligns coast_diff;
coast_diff = -(1/depth) dQ_net/dtheta to within one angle bin, so module 4 uses the same physics
BRIE does. No storms, no Barrier3D, no sea level. Details: README.md beside this script.

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-10-08
"""
from __future__ import annotations

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
REPO = next(p for p in HERE.parents if (p / "pyproject.toml").exists())
sys.path[:0] = [str(REPO / "scripts")]
from site_layer.hat_figure_style import (  # noqa: E402
    C, INK, INK_MUTED, _title, apply_style, figsize, open_frame, record_caption, save)

# --- CONFIG ------------------------------------------------------------------
DY_M, DT_YR = 500.0, 1.0                 # fixed by the CASCADE coupler
CLIMATE = dict(Hs=2.0, Tp=7.5, asym=0.6, ahf=0.5)   # option A
NY = 61
LO, HI = 30, 31                          # the groin is the face between these two
ORIENTATIONS = tuple(range(-40, 45, 5))  # shoreline angle theta0, degrees
SHOWN = (-30, -10, 0, 10, 20, 30)        # orientations drawn as profiles
PROFILE_YEARS = (1, 5, 10, 25, 50)
B = 0.6                                  # trapping fraction, the 10-05 pin
M = 12.0                                 # dipole, the best option-A dipole (2026-09-29)
MODULES = {
    "dipole": dict(folder="1-source-sink-dipole", label="Source/sink (dipole)", col=C["ADDED"]),
    "pinned": dict(folder="2-trapping-pinned", label="Trapping, pinned", col=C["ACCENT"]),
    "conserving": dict(folder="3-trapping-conserving", label="Trapping, conserving", col=C["REF"]),
    "drift": dict(folder="4-trapping-drift", label="Trapping, drift", col=C["LATE"]),
}
# -----------------------------------------------------------------------------


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


def run(t, theta0, module, years, b=B):
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
            dx[LO] -= b * 2 * r[LO] * off
            dx[HI] += b * 2 * r[HI] * off
        elif module == "conserving":
            dx[LO] -= b * 2 * r[LO] * off
            dx[HI] += b * 2 * r[LO] * off
        elif module == "drift":
            # the drift arriving at the groin: the angle of the updrift approach face, not of the step
            approach = theta[HI] if up == HI else theta[LO - 1]
            q = b * float(drift(t, approach)[0]) * DT_YR / (DY_M * t["depth"])
            dx[LO] -= q
            dx[HI] += q
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


def module_figures(name, frame, profiles):
    spec = MODULES[name]
    out = HERE / spec["folder"] / "figures"
    out.mkdir(parents=True, exist_ok=True)
    dom = np.arange(NY) - LO - 0.5                              # domains from the groin face
    fig, axes = plt.subplots(2, 3, figsize=figsize("double", height=5.2), sharex=True,
                             constrained_layout=True)
    shades = plt.cm.Greys(np.linspace(0.35, 1.0, len(PROFILE_YEARS)))
    lim = max(np.abs(profiles[(th, y)]).max() for th in SHOWN for y in PROFILE_YEARS)
    lim = max(5.0, np.ceil(lim / 10) * 10)
    for i, (ax, th) in enumerate(zip(axes.flat, SHOWN)):
        for c, y in zip(shades, PROFILE_YEARS):
            ax.plot(dom, profiles[(th, y)], color=c, lw=1.2)
        ax.axvline(0, color=C["GROIN"], lw=0.8)
        ax.axhline(0, color=INK_MUTED, lw=0.5)
        ax.set_xlim(-12, 12)
        ax.set_ylim(-lim, lim)
        q = frame[(frame.theta0_deg == th) & (frame.year == 1)].drift0_m3_yr.iloc[0]
        ax.text(0.03, 0.04, f"drift {q / 1e3:+.0f}k m$^3$/yr", transform=ax.transAxes,
                fontsize=6.5, color=INK_MUTED)
        ax.grid(axis="y")
        open_frame(ax)
        _title(ax, i, f"coast at {th:+d}$^\\circ$")
    for ax in axes[:, 0]:
        ax.set_ylabel("Shoreline change\n(m, seaward +)")
    for ax in axes[1]:
        ax.set_xlabel("Domains from the groin (updrift right)")
    handles = [Line2D([], [], color=c, lw=1.2, label=f"year {y}") for c, y in zip(shades, PROFILE_YEARS)]
    fig.legend(handles=handles, loc="outside upper center", ncol=5, fontsize=7, frameon=False)
    png = out / f"straight_coast_{spec['folder']}_profiles.png"
    save(fig, png, dpi=250, close=True)
    record_caption(png, (
        f"{spec['label']} groin on a straight coast in BRIE's alongshore solve alone (option A "
        "waves, no storms, no Barrier3D), at six shoreline orientations. Lines: shoreline change "
        "from the straight start, seaward positive, at years 1, 5, 10, 25 and 50 (light to dark). "
        "Red line: the groin face. The x axis runs in domains from the face with the updrift side "
        "on the right at every orientation. 'drift' is the wave-climate net longshore transport on "
        "the undisturbed coast."))


def comparison_figures(t, summary, frames):
    out = HERE / "comparison" / "figures"
    out.mkdir(parents=True, exist_ok=True)
    th = np.arange(-45, 46)
    dif = t["coast_diff"][angle_index(th, t["climl"])]
    fig, axes = plt.subplots(2, 2, figsize=figsize("double", height=5.6), constrained_layout=True)
    ax = axes[0, 0]
    ax.plot(th, drift(t, th) / 1e3, color=INK, lw=1.4)
    ax.axhline(0, color=INK_MUTED, lw=0.5)
    ax.set_ylabel("Net drift (10$^3$ m$^3$/yr,\ntoward higher domains +)")
    ax2 = ax.twinx()
    ax2.plot(th, dif / 1e3, color=C["BASE"], lw=1.2, ls="--")
    ax2.plot(th, np.maximum(dif, 0) / 1e3, color=C["BASE"], lw=0.0)
    ax2.set_ylabel("BRIE diffusivity\n(10$^3$ m$^2$/yr, dashed)", color=C["BASE"])
    ax.set_xlabel("Shoreline orientation ($^\\circ$)")
    open_frame(ax)
    _title(ax, 0, "What drives each module")

    for k, (col, ylab, title) in enumerate((
            ("updrift_seaward_m", "Updrift cell, year 50\n(m, seaward +)", "Sand held updrift"),
            ("applied_net_cum_m", "Module's net, 50 yr\n(m of shoreline, + added)", "Sand the module created"),
            ("solve_error_m", "Solve's net, 50 yr\n(m of shoreline, + added)", "Sand BRIE's solve created"))):
        ax = axes.flat[k + 1]
        for name, spec in MODULES.items():
            s = summary[summary.module == name]
            ax.plot(s.theta0_deg, s[col], color=spec["col"], lw=1.6, marker="o", ms=3)
        ax.axhline(0, color=INK_MUTED, lw=0.5)
        ax.axvline(0, color=INK_MUTED, lw=0.5, ls=(0, (2, 2)))
        ax.set_xlabel("Shoreline orientation ($^\\circ$)")
        ax.set_ylabel(ylab)
        ax.grid(axis="y")
        open_frame(ax)
        _title(ax, k + 1, title)
    handles = [Line2D([], [], color=s["col"], lw=1.6, marker="o", ms=3, label=s["label"])
               for s in MODULES.values()]
    fig.legend(handles=handles, loc="outside upper center", ncol=4, fontsize=7, frameon=False)
    png = out / "straight_coast_modules_vs_orientation.png"
    save(fig, png, dpi=250, close=True)
    record_caption(png, (
        "Four groin modules on a straight coast at orientations -40 to +40 degrees, 50 years of "
        "BRIE's alongshore solve alone (option A waves). (a) The wave-climate net drift (solid) "
        "and BRIE's diffusivity (dashed; BRIE clips negative values to zero, which freezes the "
        "coast). (b) Seaward change of the updrift cell beside the groin. (c) Sand the module "
        "itself added, as summed shoreline change across the reach (zero for a conserving "
        "module). (d) Sand BRIE's row-scaled solve added on top, which is non-zero wherever the "
        f"groin makes the diffusivity vary along the coast. Dipole M {M:g}; trapping b {B}."))

    # Year-50 profiles, modules side by side at the shown orientations
    fig, axes = plt.subplots(2, 3, figsize=figsize("double", height=5.2), sharex=True,
                             constrained_layout=True)
    dom = np.arange(NY) - LO - 0.5
    for i, (ax, th0) in enumerate(zip(axes.flat, SHOWN)):
        for name, spec in MODULES.items():
            ax.plot(dom, frames[name][1][(th0, 50)], color=spec["col"], lw=1.5)
        ax.axvline(0, color=C["GROIN"], lw=0.8)
        ax.axhline(0, color=INK_MUTED, lw=0.5)
        ax.set_xlim(-12, 12)
        ax.grid(axis="y")
        open_frame(ax)
        _title(ax, i, f"coast at {th0:+d}$^\\circ$, year 50")
    for ax in axes[:, 0]:
        ax.set_ylabel("Shoreline change\n(m, seaward +)")
    for ax in axes[1]:
        ax.set_xlabel("Domains from the groin (updrift right)")
    fig.legend(handles=handles, loc="outside upper center", ncol=4, fontsize=7, frameon=False)
    png = out / "straight_coast_profiles_year50.png"
    save(fig, png, dpi=250, close=True)
    record_caption(png, (
        "The four groin modules after 50 years on a straight coast at six orientations, BRIE's "
        "alongshore solve alone (option A waves). Shoreline change from the straight start, "
        "seaward positive, updrift side on the right; red line, the groin face; each panel has its "
        "own y range."))


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
        tables = HERE / spec["folder"] / "tables"
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
        f, _ = run(t, th, "drift", a.years, b=1.0)
        y_pc, tan_a = pelnard_considere(t, th, a.years)
        bench.append(dict(theta0_deg=th, model_updrift_m=float(f.updrift_seaward_m.iloc[-1]),
                          pelnard_considere_m=y_pc, tan_alpha=tan_a,
                          linear_ok=bool(tan_a < 0.2) if np.isfinite(tan_a) else False))
    (HERE / "comparison" / "tables").mkdir(parents=True, exist_ok=True)
    summary.to_csv(HERE / "comparison" / "tables" / "straight_coast_summary_year50.csv",
                   index=False, float_format="%.3f")
    pd.DataFrame(bench).to_csv(HERE / "comparison" / "tables" / "pelnard_considere_check.csv",
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
