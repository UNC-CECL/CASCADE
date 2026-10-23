#!/usr/bin/env python3
"""
Does BRIE's shoreline step conserve sand at a sharp bend? One figure, from BRIE's own step and rate table.

    python brie_sand_balance_figure.py [--years 50]
        ->  figures/brie_sand_balance_at_a_step.png

BRIE's alongshore step (brie.py:1290-1328, copied line for line) is run with nothing added or removed on three
periodic coasts: straight, a gentle bump, and a groin-like step. Each cell reads its rate from the face to its
northern neighbour and uses it for both faces, so at a bend the two cells either side of a face use different rates
for it. The same step with one rate per face (both cells use it) is the comparison. Rates come from a real Brie
instance with option A waves. A diagnostic only: BRIE is unchanged.

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-10-14
"""
from __future__ import annotations

# Non-code files beside this script live in the workspace on D: (hatteras_paths.workspace_twin)
import sys as _hat_sys  # noqa: E402
from pathlib import Path as _HatPath  # noqa: E402
_HAT_DIR = next(_p for _p in _HatPath(__file__).resolve().parents if (_p / "hatteras_paths.py").exists())
if str(_HAT_DIR) not in _hat_sys.path:
    _hat_sys.path.insert(0, str(_HAT_DIR))
from hatteras_paths import workspace_twin  # noqa: E402
import sys
from pathlib import Path
import argparse

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from scipy.sparse import csr_matrix  # noqa: E402
from scipy.sparse.linalg import spsolve  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "1-run"))   # the test, in ../1-run/
import straight_coast_test as sc  # noqa: E402
from brie import Brie  # noqa: E402
from site_layer.hat_figure_style import (  # noqa: E402
    C, INK, INK_MUTED, _letter_inside, apply_style, figsize, open_frame, record_caption, save)

# --- CONFIG ------------------------------------------------------------------
OUT = workspace_twin(Path(__file__).resolve().parent) / "figures" / "brie_sand_balance_at_a_step.png"
NY = 30                                     # cells, 500 m each
STEP_M = 100.0                              # the step's build-out, m
BRIE_C, FACE_C = C["GROIN"], C["REF"]
COASTS = {                                  # seaward +, m
    "straight": np.zeros(NY),
    f"gentle {STEP_M:.0f} m bump": STEP_M * np.exp(-((np.arange(NY) - 15) / 4.0) ** 2),
    f"{STEP_M:.0f} m step": np.r_[np.zeros(15), STEP_M * np.ones(3), np.zeros(12)],
}
COAST_LS = {"straight": ":", f"gentle {STEP_M:.0f} m bump": "--", f"{STEP_M:.0f} m step": "-"}
# -----------------------------------------------------------------------------

b = Brie(alongshore_section_length=500, alongshore_section_count=NY, wave_height=2, wave_period=7.5,
         wave_asymmetry=0.6, wave_angle_high_fraction=0.5)


# Each cell's rate, from the angle of its face to its northern neighbour, as brie.py (x landward +)
def rates(x):
    theta = 180.0 * np.arctan2(np.roll(x, -1) - x, b._dy) / np.pi
    k = np.maximum(1, np.minimum(b._wave_climl, np.round(90 - theta).astype(int)))
    return np.maximum(0.0, b._coast_diff[k] * b._dt / 2 / b._dy ** 2)


# One year: "brie" as brie.py writes it, "face" with one rate per face shared by the two cells it joins
def step(x, form):
    r = rates(x)
    if form == "brie":
        dv = np.r_[-r[-1], -r[1:], 1 + 2 * r, -r[0:-1], -r[0]]
        rhs = x + r * (np.roll(x, -1) - 2 * x + np.roll(x, 1))
        return spsolve(csr_matrix((dv, (b._di, b._dj))), rhs)
    n = np.arange(NY)
    rp, rm = r, np.roll(r, 1)
    A = csr_matrix((np.r_[1 + rp + rm, -rp, -rm], (np.r_[n, n, n], np.r_[n, (n + 1) % NY, (n - 1) % NY])),
                   shape=(NY, NY))
    return spsolve(A, x + rp * (np.roll(x, -1) - x) - rm * (x - np.roll(x, 1)))


def run(seaward, form, years):
    x = -seaward.copy()
    out = [x.copy()]
    for _ in range(years):
        x = step(x, form)
        out.append(x.copy())
    return -np.array(out)                   # back to seaward +


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--years", type=int, default=50)
    a = ap.parse_args()
    apply_style()
    yrs = np.arange(a.years + 1)
    runs = {(name, form): run(x0, form, a.years) for name, x0 in COASTS.items() for form in ("brie", "face")}
    total = {key: v.sum(axis=1) - v[0].sum() for key, v in runs.items()}
    step_name = f"{STEP_M:.0f} m step"
    bump_name = f"gentle {STEP_M:.0f} m bump"

    # the rate each cell uses for each face, at the start of the step run
    x0 = -COASTS[step_name]
    r = rates(x0)
    south_uses, north_uses = r, np.roll(r, -1)          # face i|i+1: cell i uses r[i], cell i+1 uses r[i+1]
    faces = np.arange(12, 19)
    bump_max_deg = np.degrees(np.arctan(np.abs(np.diff(COASTS[bump_name])).max() / b._dy))

    fig, axs = plt.subplots(1, 3, figsize=figsize("double", height=2.9), constrained_layout=True)
    km = (np.arange(NY) + 0.5) * b._dy / 1000

    ax = axs[0]
    ax.step(km, COASTS[step_name], where="mid", color=INK_MUTED, lw=1.0, label="start")
    ax.step(km, runs[(step_name, "brie")][-1], where="mid", color=BRIE_C, lw=1.5, label="BRIE as written")
    ax.step(km, runs[(step_name, "face")][-1], where="mid", color=FACE_C, lw=1.2, ls="--",
            label="one rate")
    ax.set_xlim(4, 11)
    ax.set_xlabel("Alongshore (km, south to north)")
    ax.set_ylabel("Shoreline position (m, seaward +)")
    ax.legend(loc="center left", bbox_to_anchor=(0.0, 0.6), fontsize=6.5, frameon=False)
    ax.grid(axis="y")
    open_frame(ax)
    _letter_inside(ax, 0)
    ax.set_title(f"The step, start and year {a.years}", fontsize=8.5)

    ax = axs[1]
    s_r, n_r = south_uses[faces] * 1e3, north_uses[faces] * 1e3
    ax.vlines(faces, np.minimum(s_r, n_r), np.maximum(s_r, n_r), color=INK_MUTED, lw=1.0, zorder=1)
    ax.plot(faces, s_r, "o", color=INK_MUTED, ms=6, label="cell to the south", zorder=2)
    ax.plot(faces, n_r, "s", color=BRIE_C, ms=5, label="cell to the north", zorder=3)
    ax.set_xticks(faces, [f"{i + 1}|{i + 2}" for i in faces], fontsize=6.5)
    ax.set_xlim(faces[0] - 0.6, faces[-1] + 0.6)
    ax.set_ylim(-0.3, max(s_r.max(), n_r.max()) * 1.2)
    ax.set_xlabel("Boundary between cells (south|north)")
    ax.set_ylabel("Rate used for that boundary (per 1,000)")
    ax.legend(loc="upper right", fontsize=6.5, frameon=False)
    ax.grid(axis="y")
    open_frame(ax)
    _letter_inside(ax, 1)
    ax.set_title("Same boundary, two rates (year 0)", fontsize=8.5)

    ax = axs[2]
    for name in COASTS:
        ax.plot(yrs, total[(name, "brie")], color=BRIE_C, lw=1.5, ls=COAST_LS[name], label=f"BRIE, {name}")
    ax.plot(yrs, total[(step_name, "face")], color=FACE_C, lw=1.2, label=f"one rate, {step_name}")
    ax.axhline(0, color=INK, lw=0.6)
    ax.set_xlim(0, a.years)
    ax.set_xlabel("Year")
    ax.set_ylabel("Sand made (+) or lost (−)\n(m of summed shoreline)")
    ax.legend(loc="lower left", fontsize=6.5, frameon=False)
    ax.grid(axis="y")
    open_frame(ax)
    _letter_inside(ax, 2)
    ax.set_title("Total sand, nothing added", fontsize=8.5)

    fig.suptitle("BRIE's shoreline step does not balance sand at a sharp bend", fontsize=10)
    OUT.parent.mkdir(parents=True, exist_ok=True)
    save(fig, OUT, dpi=250, close=True)
    end = {key: v[-1] for key, v in total.items()}
    mismatch = np.abs(south_uses - north_uses)
    record_caption(OUT, (
        "**Tests:** whether BRIE's alongshore shoreline step keeps the island's total sand constant when nothing "
        "is added or removed. BRIE's step (brie.py:1290-1328) is copied line for line and run on a periodic coast "
        f"of {NY} 500 m cells, with the rate table of a real Brie instance (Hs 2 m, Tp 7.5 s, asymmetry 0.6, "
        "high-angle fraction 0.5). It is compared with the same step given one rate per boundary, used by both "
        "cells it joins. BRIE is unchanged. **How to read:** (a) the step coast at the start and after "
        f"{a.years} years. (b) For each boundary near the step, the rate the cell to its south uses for it and "
        "the rate the cell to its north uses; BRIE gives each cell the rate of its northern boundary and uses it "
        "for both its boundaries, so the two differ wherever the angle changes from one boundary to the next. "
        "(c) The change in summed shoreline (seaward +), which is the change in total sand: zero if sand is "
        "conserved. **Shows:** the straight coast conserves sand exactly "
        f"({end[('straight', 'brie')]:+.2f} m after {a.years} years) and the gentle bump nearly "
        f"({end[(bump_name, 'brie')]:+.2f} m; its steepest angle is {bump_max_deg:.1f} degrees). At the "
        f"{STEP_M:.0f} m step the two sides of a boundary use rates up to {mismatch[faces].max() * 1e3:.2f} "
        "per 1,000 apart. With these waves BRIE's rate on a coast facing straight into the wave climate is zero "
        "(its diffusivity there is negative and clipped), so at each boundary beside the step one side moves "
        "and the other does not. The total changes by "
        f"{end[(step_name, 'brie')]:+.2f} m in {a.years} years with nothing added or removed. With one rate per "
        f"boundary it changes by {end[(step_name, 'face')]:+.2f} m. The two shorelines in (a) differ by up to "
        f"{np.abs(runs[(step_name, 'brie')][-1] - runs[(step_name, 'face')][-1]).max():.0f} m at single cells."))
    for key, v in end.items():
        print(f"{key[0]:20s} {key[1]:5s} total change {v:+.2f} m")
    print(f"wrote {OUT.relative_to(sc.REPO)}")


if __name__ == "__main__":
    main()
