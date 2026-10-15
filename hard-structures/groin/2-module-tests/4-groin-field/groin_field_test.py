#!/usr/bin/env python3
"""
Four groins along a straight coast: which groin module behaves like a groin field?

    python groin_field_test.py [--years 50]   ->  tables/, figures/

The straight-coast emulator of BRIE's alongshore step (../1-straight-coast/1-run/straight_coast_test.py,
imported unchanged), with four groins at a spacing of 2, 5, 10 or 20 domains instead of one, and a
single groin at the middle face as the reference. Four modules: the source/sink dipole, the pinned
trapping groin, the drift groin and the corrected drift groin. Each run is scored at the last year
on four things a groin field should do: conserve sand, hold sand updrift of every groin (a
sawtooth), erode the coast just downdrift of the field, and let close groins act as one. No storms,
no Barrier3D, no sea level. Details: README.md beside this script.

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-10-09
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
from scipy.sparse import csr_matrix
from scipy.sparse.linalg import spsolve

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "1-straight-coast" / "1-run"))
import straight_coast_test as sc  # noqa: E402
from site_layer.hat_figure_style import (  # noqa: E402
    C, INK, INK_MUTED, _letter_inside, apply_style, figsize, open_frame, record_caption, save)

# --- CONFIG ------------------------------------------------------------------
NY = 161                                  # 80 km, room for a 60-domain field and its lee
N_GROINS = 4
SPACINGS = (2, 5, 10, 20)                 # domains between neighbouring groins (1 to 10 km)
ORIENTATIONS = (-30, -20, -10, 0, 10, 20, 30)
MODULES = ("dipole", "pinned", "drift", "corrected", "length")
MIN_DRIFT_SHARE = 0.10                    # orientations scored: drift at least this share of the drift at 0 deg
CONSERVE_TOL = 0.10                       # |net change| at most this share of the gross change
HOLD_M = 1.0                              # an offset across a groin, or a lee cut, counts past 1 m
LEE_CELLS = 5                             # cells just downdrift of the last groin searched for the lee cut
MERGE_FAR_MIN = 0.5                       # at the widest spacing, each downdrift groin holds at least half an isolated one
SHOWN = (-20, 0, 20)                      # orientations drawn as profiles
MARKERS = dict(dipole="o", pinned="s", drift="^", corrected="D", length="v")
TESTS = dict(conserves="Conserves sand", sawtooth="Holds sand at\nevery groin",
             lee_cut="Cuts the coast\nbelow the field", merges="Close groins\nact as one")
# -----------------------------------------------------------------------------


# The groin faces: face k sits between domains lo[k] and lo[k] + 1, the field centred on the coast
def field_faces(n, spacing):
    first = NY // 2 - (n - 1) * spacing // 2 - 1
    return [first + k * spacing for k in range(n)]


# Module 6's leak correction near any groin; overlapping zones get it once, since it is per cell
def leak_correction(r, step, faces):
    i = np.arange(r.size)
    near = np.zeros(r.size, bool)
    for lo in faces:
        near |= np.abs(i - (lo + 0.5)) <= sc.CORRECTION_HALF_WIDTH
    return np.where(near, 2 * (r - np.roll(r, 1)) * np.roll(step, 1), 0.0)


# One run: straight_coast_test.run() with a list of groin faces; returns the final shoreline, seaward +
def run_field(t, theta0, module, faces, years, blocking_fraction=sc.BLOCKING_FRACTION, groin_length_m=sc.GROIN_LENGTH_M,
              record_years=()):
    slope = np.tan(np.deg2rad(theta0))
    eta = np.zeros(NY)
    up_is_hi = float(sc.drift(t, theta0)[0]) < 0          # drift toward lower index: the higher cell is updrift
    snaps = {}
    for year in range(1, years + 1):
        step = slope * sc.DY_M + (np.roll(eta, -1) - eta)       # x[i+1] - x[i], landward +
        theta = 180.0 * np.arctan2(step, sc.DY_M) / np.pi
        r = np.maximum(0.0, t["coast_diff"][sc.angle_index(theta, t["climl"])] * sc.DT_YR / 2 / sc.DY_M ** 2)
        dx = np.zeros(NY)
        for lo in faces:
            hi = lo + 1
            up, down = (hi, lo) if up_is_hi else (lo, hi)
            off = step[lo]
            if module == "dipole":
                dx[up] -= sc.M
                dx[down] += sc.M
            elif module == "pinned":
                dx[lo] -= blocking_fraction * 2 * r[lo] * off
                dx[hi] += blocking_fraction * 2 * r[hi] * off
            elif module in ("drift", "corrected", "length"):
                approach = theta[hi] if up_is_hi else theta[lo - 1]
                q = blocking_fraction * float(sc.drift(t, approach)[0]) * sc.DT_YR / (sc.DY_M * t["depth"])
                if module == "length":
                    q *= sc.length_factor(-eta[up], groin_length_m)
                dx[lo] -= q
                dx[hi] += q
        if module in ("corrected", "length"):
            dx += leak_correction(r, step, faces)
        rhs = eta + r * (np.roll(eta, -1) - 2 * eta + np.roll(eta, 1)) + dx
        dv = np.r_[-r[-1], -r[1:], 1 + 2 * r, -r[0:-1], -r[0]]
        eta = spsolve(csr_matrix((dv, (t["di"], t["dj"]))), rhs)
        if year in record_years:
            snaps[year] = -eta.copy()
    return snaps if record_years else -eta   # record_years: {year: shoreline}, for long runs


# The four scores of one field run, against the isolated groin's offset
def score(x, faces, up_is_hi, iso_offset):
    order = sorted(faces, reverse=up_is_hi)                  # updrift-most groin first
    offsets = [x[lo + 1] - x[lo] if up_is_hi else x[lo] - x[lo + 1] for lo in order]
    last = order[-1]
    lee = (np.arange(last, last - LEE_CELLS, -1) if up_is_hi
           else np.arange(last + 1, last + 1 + LEE_CELLS)) % NY
    total, gross = float(x.sum()), float(np.abs(x).sum())
    return dict(
        net_change_m=total, gross_change_m=gross,
        budget_share=total / gross if gross > 0 else 0.0,
        offset_min_m=float(min(offsets)), offset_mean_m=float(np.mean(offsets)),
        groins_holding=int(sum(o > HOLD_M for o in offsets)),
        lee_min_m=float(x[lee].min()),
        downdrift_ratio=(float(np.mean(offsets[1:]) / iso_offset)
                         if len(offsets) > 1 and iso_offset > HOLD_M else np.nan),   # only where one groin traps updrift
        **{f"offset_{k + 1}_m": float(o) for k, o in enumerate(offsets)})


# Year-N profiles: spacing by row, orientation by column, one line per module
def profile_figure(prof, years):
    fig, axes = plt.subplots(len(SPACINGS), len(SHOWN), figsize=figsize("double", height=8.4), sharex="row",
                             constrained_layout=True)
    centre = NY // 2 - 0.5
    for r_, sp in enumerate(SPACINGS):
        faces = field_faces(N_GROINS, sp)
        for c_, th in enumerate(SHOWN):
            ax = axes[r_, c_]
            for m in MODULES:
                p = prof[(prof.module == m) & (prof.theta0_deg == th) & (prof.spacing == sp)]
                ax.plot(p.domain - centre, p.seaward_m, color=sc.MODULES[m]["col"], lw=1.3,
                        marker=MARKERS[m], ms=2.6, markevery=4)
            for lo in faces:
                ax.axvline(lo + 0.5 - centre, color=C["GROIN"], lw=0.8)
            ax.axhline(0, color=INK_MUTED, lw=0.5)
            half = (N_GROINS - 1) * sp / 2 + 12                 # each row zoomed to its field and the fillets beside it
            ax.set_xlim(-half, half)
            ax.grid(axis="y")
            open_frame(ax)
            _letter_inside(ax, r_ * len(SHOWN) + c_)
            if r_ == 0:
                ax.set_title(f"Coast at {th:+d}$^\\circ$", loc="center")
        axes[r_, 0].set_ylabel(f"Groins {sp} domain{'s' if sp > 1 else ''} apart\n"
                               f"({sp * sc.DY_M / 1000:g} km)\nchange (m, seaward +)", fontsize=7.5)
    axes[-1, len(SHOWN) // 2].set_xlabel("Domains from the field centre (updrift right)")
    handles = [Line2D([], [], color=sc.MODULES[m]["col"], lw=1.3, marker=MARKERS[m], ms=3.5,
                      label=sc.MODULES[m]["label"]) for m in MODULES]
    handles.append(Line2D([], [], color=C["GROIN"], lw=0.8, label="groin"))
    fig.legend(handles=handles, loc="outside lower center", ncol=6, fontsize=7, frameon=False)
    fig.suptitle(f"Where does each version put the sand along a field of {N_GROINS} groins?")
    png = HERE / "figures" / "groin_field_profiles.png"
    save(fig, png, dpi=250, close=True)
    return png


# The scorecard: module by row, test by column, each cell orientation by spacing, filled = passes
def scorecard_figure(field, card):
    th_scored = sorted(field[field.scored].theta0_deg.unique())
    fig, axes = plt.subplots(len(MODULES), len(TESTS), figsize=figsize("double", height=8.2),
                             constrained_layout=True)
    for r_, m in enumerate(MODULES):
        for c_, (test, title) in enumerate(TESTS.items()):
            ax = axes[r_, c_]
            sub = field[(field.module == m) & field.scored]
            for i, th in enumerate(th_scored):
                for j, sp in enumerate(SPACINGS):
                    ok = bool(sub[(sub.theta0_deg == th) & (sub.spacing == sp)][test].iloc[0])
                    ax.add_patch(plt.Rectangle((i + 0.08, j + 0.08), 0.84, 0.84, lw=0.8,
                                               edgecolor=sc.MODULES[m]["col"],
                                               facecolor=sc.MODULES[m]["col"] if ok else "white"))
            ax.set_xlim(0, len(th_scored))
            ax.set_ylim(0, len(SPACINGS))
            ax.set_aspect("equal")
            ax.set_xticks(np.arange(len(th_scored)) + 0.5, [f"{t:+d}" for t in th_scored], fontsize=6)
            ax.set_yticks(np.arange(len(SPACINGS)) + 0.5, [str(sp) for sp in SPACINGS], fontsize=6)
            ax.tick_params(length=0)
            for sp_ in ax.spines.values():
                sp_.set_visible(False)
            ax.text(1.0, 1.04, f"{int(card.loc[m, test])}/{int(card.loc[m, 'cases'])}", transform=ax.transAxes,
                    ha="right", va="bottom", fontsize=7.5, color=INK, fontweight="bold")
            if r_ == 0:
                ax.set_title(title, loc="left", fontsize=8)
            if c_ == 0:
                ax.set_ylabel(f"{sc.MODULES[m]['label']}\nspacing (domains)", fontsize=7.5)
            if r_ == len(MODULES) - 1:
                ax.set_xlabel("Coast orientation ($^\\circ$)", fontsize=7)
    fig.suptitle("Which version behaves like a groin field?")
    png = HERE / "figures" / "groin_field_scorecard.png"
    save(fig, png, dpi=250, close=True)
    return png


# Run: every module, orientation and spacing, plus the isolated groin; score and write
def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--years", type=int, default=50)
    a = ap.parse_args()
    sc.NY = NY                                               # BRIE's tables and sparse indices at this rig's length
    t = sc.brie_tables()
    q0 = abs(float(sc.drift(t, 0)[0]))

    rows, profiles = [], []
    for module in MODULES:
        for th in ORIENTATIONS:
            up_is_hi = float(sc.drift(t, th)[0]) < 0
            iso_face = [NY // 2 - 1]
            x_iso = run_field(t, th, module, iso_face, a.years)
            lo = iso_face[0]
            iso_offset = x_iso[lo + 1] - x_iso[lo] if up_is_hi else x_iso[lo] - x_iso[lo + 1]
            layouts = [(0, iso_face, x_iso)] + [
                (sp, field_faces(N_GROINS, sp), None) for sp in SPACINGS]
            for sp, faces, x in layouts:
                if x is None:
                    x = run_field(t, th, module, faces, a.years)
                s = score(x, faces, up_is_hi, iso_offset)
                rows.append(dict(module=module, theta0_deg=th, spacing=sp, n_groins=len(faces),
                                 drift_m3_yr=float(sc.drift(t, th)[0]),
                                 scored=abs(float(sc.drift(t, th)[0])) >= MIN_DRIFT_SHARE * q0,
                                 iso_offset_m=iso_offset, **s))
                profiles.append(pd.DataFrame(dict(
                    module=module, theta0_deg=th, spacing=sp, domain=np.arange(NY),
                    is_groin_lo=np.isin(np.arange(NY), faces), seaward_m=x)))
    frame = pd.DataFrame(rows)

    # Pass or fail per run, then the merging test across spacings
    field = frame[frame.spacing > 0].copy()
    field["conserves"] = field.budget_share.abs() <= CONSERVE_TOL
    field["sawtooth"] = field.groins_holding == N_GROINS
    field["lee_cut"] = field.lee_min_m < -HOLD_M
    piv = field.pivot_table(index=["module", "theta0_deg"], columns="spacing", values="downdrift_ratio")
    merge = ((piv[min(SPACINGS)] < piv[max(SPACINGS)]) & (piv[max(SPACINGS)] >= MERGE_FAR_MIN)).rename("merges")
    field = field.join(merge, on=["module", "theta0_deg"])

    out = HERE / "tables"
    out.mkdir(parents=True, exist_ok=True)
    frame.to_csv(out / "groin_field_runs.csv", index=False, float_format="%.3f")
    field.to_csv(out / "groin_field_scores.csv", index=False, float_format="%.3f")
    pd.concat(profiles).to_csv(out / "groin_field_profiles.csv", index=False, float_format="%.3f")
    tests = ["conserves", "sawtooth", "lee_cut", "merges"]
    card = field[field.scored].groupby("module")[tests].sum().reindex(MODULES)
    card["cases"] = field[field.scored].groupby("module").size().reindex(MODULES)
    card.to_csv(out / "groin_field_scorecard.csv")

    # Figures and their captions
    apply_style()
    prof = pd.concat(profiles)
    png = profile_figure(prof, a.years)
    fs = field.set_index(["module", "theta0_deg", "spacing"])
    record_caption(png, (
        f"**Tests:** where each version of the groin module puts the sand when {N_GROINS} groins share one "
        "straight coast, and whether the result looks like a groin field: a sawtooth, each groin holding sand on "
        "its updrift side, and the coast cut back just below the last groin. **How to read:** BRIE's alongshore "
        f"solve alone (option A waves; no storms, no Barrier3D, no sea-level rise), {NY} domains, {a.years} years. "
        f"Rows: the groins {', '.join(map(str, SPACINGS[:-1]))} and {SPACINGS[-1]} domains apart "
        f"({SPACINGS[0] * sc.DY_M / 1000:g} to {SPACINGS[-1] * sc.DY_M / 1000:g} km), each row zoomed to its field; "
        "columns: three shoreline orientations. Lines "
        "are the shoreline change from the straight start, seaward positive, with the updrift side on the right "
        "(the net drift runs right to left at every orientation shown); red lines are the groins. Dipole "
        f"M {sc.M:g} m/yr; trapping versions blocking fraction {sc.BLOCKING_FRACTION}. **Shows:** at -20 degrees "
        "every version builds a sawtooth; at 0 degrees the pinned version does nothing, since it blocks only the "
        "diffusive part of the transport and the diffusivity is zero there; at +20 degrees the pinned version "
        "builds out below each groin instead of above it, and the drift version's field loses "
        f"{-fs.loc[('drift', 20, SPACINGS[-1]), 'net_change_m']:.0f} m of summed shoreline at the widest spacing, against "
        f"{fs.loc[('corrected', 20, SPACINGS[-1]), 'net_change_m']:+.0f} m for the corrected version. "
        "Each panel has its own vertical range (tables/groin_field_profiles.csv)."))
    png = scorecard_figure(field, card)
    th_scored = sorted(field[field.scored].theta0_deg.unique())
    record_caption(png, (
        f"**Tests:** which version of the groin module behaves like a field of {N_GROINS} groins, on four tests, "
        f"across {len(th_scored)} shoreline orientations ({th_scored[0]:+d} to {th_scored[-1]:+d} degrees, where the "
        "drift is at least a tenth of its grid-aligned value) and four spacings. **How to read:** one row per "
        "version, one column per test; inside each block, columns are orientations and rows are spacings in "
        f"domains; a filled square passes. The count is passes out of cases, after {a.years} years. Conserves sand: "
        f"the net change across the coast is at most {CONSERVE_TOL:.0%} of the gross change. Holds sand at every "
        f"groin: the cell just updrift of each groin is more than {HOLD_M:g} m seaward of the cell just downdrift. "
        f"Cuts the coast below the field: some cell among the {LEE_CELLS} just downdrift of the last groin has "
        f"retreated more than {HOLD_M:g} m. Close groins act as one: each downdrift groin holds less, relative to a "
        f"single isolated groin, when the groins are {SPACINGS[0]} domains apart than when they are {SPACINGS[-1]} apart, "
        f"and at {SPACINGS[-1]} apart "
        f"holds at least {MERGE_FAR_MIN:.0%} of what the isolated groin holds; it is one result per orientation, "
        "repeated down the block. **Shows:** the corrected version passes the first three tests in "
        f"{int(card.loc['corrected', 'conserves'])}, {int(card.loc['corrected', 'sawtooth'])} and "
        f"{int(card.loc['corrected', 'lee_cut'])} of {int(card.loc['corrected', 'cases'])} cases, more than any "
        "other; it fails conservation only at +20 degrees, where so little sand moves that a small imbalance is a "
        "large share, and the fourth test only at +30 degrees, where the drift is a tenth of its grid-aligned value. "
        "Even so, its groins "
        f"{SPACINGS[0]} domains apart each still hold most of what an isolated groin holds: no version has a groin "
        "length to cap the offset. The pinned version "
        f"conserves sand in {int(card.loc['pinned', 'conserves'])} cases, four of them at 0 degrees where it moves "
        f"nothing, and holds sand at every groin in only {int(card.loc['pinned', 'sawtooth'])}. The dipole and "
        f"drift versions conserve sand in {int(card.loc['dipole', 'conserves'])} and "
        f"{int(card.loc['drift', 'conserves'])} of {int(card.loc['drift', 'cases'])} cases "
        "(tables/groin_field_scorecard.csv)."))

    pd.set_option("display.width", 220)
    print(card.to_string())
    print(field[field.scored].pivot_table(index=["theta0_deg", "spacing"], columns="module",
                                           values=["offset_min_m", "budget_share", "lee_min_m",
                                                   "downdrift_ratio"]).round(2).to_string())


if __name__ == "__main__":
    main()
