#!/usr/bin/env python3
"""
The corrected trapping groin (straight-coast module 6) on the real Buxton planform, against the
pinned groin the hindcast uses and the drift groin without the correction.

    python corrected_groin_real_planform.py [--years 30]   ->  tables/, figures/

BRIE's alongshore solve alone (option A waves; no Barrier3D, storms, sea level or source/sink),
started from the BRIE shoreline after year 1 of the adopted 1996-2025 edgeBE no-groin matrix run.
One groin between GIS 5 and GIS 6, held at one blocking fraction for the whole run (no failure
schedule). Every groin run is compared with a no-groin run from the same start. BRIE is not
changed: each module only adds to the yearly shoreline change, as CASCADE's groin callbacks do.

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

HERE = Path(__file__).resolve().parent
REPO = next(p for p in HERE.parents if (p / "pyproject.toml").exists())
sys.path[:0] = [str(HERE.parents[1] / "1-straight-coast" / "1-run"), str(REPO / "scripts")]
import straight_coast_test as sc  # noqa: E402
from site_layer.hat_figure_style import (  # noqa: E402
    C, DOMAIN_AXIS_LABEL, INK_MUTED, _title, apply_style, figsize, open_frame, record_caption, save)

# --- CONFIG ------------------------------------------------------------------
PLANFORM_RUN = (REPO / "output" / "raw_runs" / "matrix" / "1996_2025" / "edgeBE"
                / "HAT_1996_2025_edgeBE_offsetmetres_road_bdm_nourish_nogroin")
START_YEAR = 1996
NUM_BUFFER = 15                          # BRIE domains before GIS 1
D, U = NUM_BUFFER + 5 - 1, NUM_BUFFER + 6 - 1   # GIS 5 and GIS 6, the cells either side of the groin
SHOWN_GIS = (1, 20)                      # alongshore window drawn in (a)
PROFILE_YEARS = (10, 30)
MODULES = ("pinned", "drift", "corrected")      # straight-coast module names and colours
OUT = HERE
# -----------------------------------------------------------------------------
assert U == D + 1


# BRIE x_s (m, landward +) after the first model year of the no-groin run
def load_planform():
    npz = np.load(PLANFORM_RUN / f"{PLANFORM_RUN.name}.npz", allow_pickle=True)
    return np.asarray(npz["cascade"][0].brie._x_s_save, float)[:, 1].copy()


# BRIE's row-scaled Crank-Nicolson solve on a periodic coast, with a module's yearly addition
def run(x0, t, module, years, blocking_fraction=sc.BLOCKING_FRACTION):
    ny = x0.size
    x = x0.copy()
    states, transfers, corrections = [x.copy()], [], []
    nxt, prv = np.r_[1:ny, 0], np.r_[ny - 1, 0:ny - 1]
    near = np.abs(np.arange(ny) - (D + 0.5)) <= sc.CORRECTION_HALF_WIDTH
    for _ in range(years):
        step = x[nxt] - x                                    # x[i+1] - x[i]
        theta = 180.0 * np.arctan2(step, sc.DY_M) / np.pi
        r = np.maximum(0.0, t["coast_diff"][sc.angle_index(theta, t["climl"])] * sc.DT_YR / 2 / sc.DY_M ** 2)
        dx = np.zeros(ny)
        q = corr = 0.0
        if module == "pinned":                              # BlockingGroinCallback
            dx[D] -= blocking_fraction * 2 * r[D] * (x[U] - x[D])
            dx[U] -= blocking_fraction * 2 * r[U] * (x[D] - x[U])
            q = -dx[U]
        elif module in ("drift", "corrected"):
            toward_lower = float(sc.drift(t, theta[U])[0]) < 0
            approach = theta[U] if toward_lower else theta[D - 1]
            q = blocking_fraction * float(sc.drift(t, approach)[0]) * sc.DT_YR / (sc.DY_M * t["depth"])
            dx[D] -= q
            dx[U] += q
            if module == "corrected":
                c = np.where(near, 2 * (r - r[prv]) * step[prv], 0.0)
                dx += c
                corr = float(np.abs(c).sum())
        rhs = x + r * (x[nxt] - 2 * x + x[prv]) + dx
        rows = np.r_[np.arange(ny), np.arange(ny), np.arange(ny)]
        cols = np.r_[np.arange(ny), nxt, prv]
        a = csr_matrix((np.r_[1 + 2 * r, -r, -r], (rows, cols)), shape=(ny, ny))
        x = spsolve(a, rhs)
        states.append(x.copy())
        transfers.append(abs(q))
        corrections.append(corr)
    return np.array(states), np.array(transfers), np.array(corrections)


def gis(pad):
    return pad - NUM_BUFFER + 1


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--years", type=int, default=30)
    a = ap.parse_args()
    apply_style()
    t = sc.brie_tables()
    x0 = load_planform()
    base, _, _ = run(x0, t, "none", a.years)
    years = np.arange(a.years + 1)                         # model years from the start
    face0 = 180 / np.pi * np.arctan2(x0[U] - x0[D], sc.DY_M)
    up0 = 180 / np.pi * np.arctan2(x0[U + 1] - x0[U], sc.DY_M)
    print(f"start: groin face {face0:+.1f} deg, updrift approach {up0:+.1f} deg, "
          f"drift there {float(sc.drift(t, up0)[0]) / 1e6:+.2f} M m3/yr (toward higher index +)")

    rows, prof, res = [], [], {}
    for m in MODULES:
        st, tr, co = run(x0, t, m, a.years)
        rel = -(st - base)                                    # seaward + change against no groin
        res[m] = (rel, tr, co)
        for k, y in enumerate(years):
            rows.append(dict(module=m, model_year=int(y), updrift_gis6_m=rel[k, U], downdrift_gis5_m=rel[k, D],
                             gap_change_m=rel[k, U] - rel[k, D], island_sand_m=float(rel[k].sum()),
                             module_transfer_m=float(tr[k - 1]) if k else 0.0,
                             leak_correction_m=float(co[k - 1]) if k else 0.0))
        for y in PROFILE_YEARS:
            for p in range(NUM_BUFFER, x0.size - NUM_BUFFER):
                prof.append(dict(module=m, model_year=y, gis=gis(p), seaward_vs_no_groin_m=rel[y, p]))
    (OUT / "tables").mkdir(parents=True, exist_ok=True)
    frame = pd.DataFrame(rows)
    frame.to_csv(OUT / "tables" / "corrected_groin_by_year.csv", index=False, float_format="%.3f")
    pd.DataFrame(prof).to_csv(OUT / "tables" / "corrected_groin_profiles.csv", index=False, float_format="%.3f")
    # The correction with no groin: what it does to the natural bends inside its window
    only, _, _ = run(x0, t, "corrected", a.years, blocking_fraction=0.0)
    only_sand = -(only - base).sum(axis=1)
    pd.DataFrame(dict(model_year=years, island_sand_m=only_sand)).to_csv(
        OUT / "tables" / "correction_without_groin.csv", index=False, float_format="%.3f")
    groin_part = frame[frame.module == "corrected"].island_sand_m.to_numpy() - only_sand
    last = frame[frame.model_year == years[-1]].set_index("module")
    print(last[["updrift_gis6_m", "downdrift_gis5_m", "gap_change_m", "island_sand_m"]].round(1).to_string())

    # Figure: where the sand goes, the gap through time, the island's sand balance, the module's moves
    fig, axes = plt.subplots(2, 2, figsize=figsize("double", height=6.0), constrained_layout=True)
    g0, g1 = SHOWN_GIS
    sl = slice(NUM_BUFFER + g0 - 1, NUM_BUFFER + g1)
    gx = np.arange(g0, g1 + 1)
    ax = axes[0, 0]
    for m in MODULES:
        col = sc.MODULES[m]["col"]
        ax.plot(gx, res[m][0][PROFILE_YEARS[0], sl], color=col, lw=1.0, ls=(0, (3, 2)))
        ax.plot(gx, res[m][0][PROFILE_YEARS[1], sl], color=col, lw=1.7, marker="o", ms=2.5)
    ax.axvline(5.5, color=C["GROIN"], lw=1.0)
    ax.axhline(0, color=INK_MUTED, lw=0.6)
    ax.set_xlim(g0 - 0.5, g1 + 0.5)
    ax.set_xticks(np.arange(g0, g1 + 1, 2) if g0 % 2 else np.arange(g0 + 1, g1 + 1, 2))
    ax.set_xlabel(DOMAIN_AXIS_LABEL)
    ax.set_ylabel("Shoreline change against no groin\n(m, seaward +)")
    ax.grid(axis="y")
    open_frame(ax)
    _title(ax, 0, "Where does the sand go?")

    ax = axes[0, 1]
    for m in MODULES:
        ax.plot(years, frame[frame.module == m].gap_change_m, color=sc.MODULES[m]["col"], lw=1.6)
    ax.axhline(0, color=INK_MUTED, lw=0.6)
    ax.set_xlabel("Model year (0 = 1996 start)")
    ax.set_ylabel("GIS 6 minus GIS 5, against no groin\n(m, updrift seaward +)")
    ax.grid(axis="y")
    open_frame(ax)
    _title(ax, 1, "How big a step does the groin hold?")

    ax = axes[1, 0]
    for m in MODULES:
        ax.plot(years, frame[frame.module == m].island_sand_m, color=sc.MODULES[m]["col"], lw=1.6)
    ax.plot(years, only_sand, color=sc.MODULES["corrected"]["col"], lw=1.0, ls=(0, (1, 1.5)))
    ax.axhline(0, color=INK_MUTED, lw=0.6)
    ax.set_xlabel("Model year (0 = 1996 start)")
    ax.set_ylabel("Sand gained across the island\n(m of summed shoreline)")
    ax.grid(axis="y")
    open_frame(ax)
    _title(ax, 2, "Is sand conserved?")

    ax = axes[1, 1]
    for m in MODULES:
        f = frame[(frame.module == m) & (frame.model_year > 0)]
        ax.plot(f.model_year, f.module_transfer_m, color=sc.MODULES[m]["col"], lw=1.6)
    f = frame[(frame.module == "corrected") & (frame.model_year > 0)]
    ax.plot(f.model_year, f.leak_correction_m, color=sc.MODULES["corrected"]["col"], lw=1.0, ls=(0, (1, 1.5)))
    ax.set_xlabel("Model year (0 = 1996 start)")
    ax.set_ylabel("Shoreline moved per year (m)")
    ax.set_ylim(bottom=0)
    ax.grid(axis="y")
    open_frame(ax)
    _title(ax, 3, "How much does each module move?")

    handles = [Line2D([], [], color=sc.MODULES[m]["col"], lw=1.6, label=sc.MODULES[m]["label"]) for m in MODULES]
    handles += [Line2D([], [], color=INK_MUTED, lw=1.7, marker="o", ms=2.5, label=f"year {PROFILE_YEARS[1]} (a)"),
                Line2D([], [], color=INK_MUTED, lw=1.0, ls=(0, (3, 2)), label=f"year {PROFILE_YEARS[0]} (a)"),
                Line2D([], [], color=sc.MODULES["corrected"]["col"], lw=1.0, ls=(0, (1, 1.5)),
                       label="corrected: leak returned (d),\ncorrection with no groin (c)"),
                Line2D([], [], color=C["GROIN"], lw=1.0, label="groin")]
    fig.legend(handles=handles, loc="outside lower center", ncol=4, fontsize=7, frameon=False)
    fig.suptitle("The corrected groin on the Buxton shoreline")
    (OUT / "figures").mkdir(parents=True, exist_ok=True)
    png = OUT / "figures" / "corrected_groin_real_planform.png"
    save(fig, png, dpi=250, close=True)

    fc = frame[frame.module == "corrected"].set_index("model_year")
    record_caption(png, (
        "**Tests:** whether the corrected trapping groin (straight-coast module 6: blocks a fraction of the "
        "wave-driven drift arriving at the groin, and gives back near the groin the sand BRIE's row-scaled "
        "solve creates or destroys) behaves on the real Buxton shoreline as it did on straight coasts, against "
        "the pinned groin the hindcast uses and the drift groin without the correction. **How to read:** "
        "BRIE's alongshore solve alone, under option A waves, with no Barrier3D, storms, sea level or "
        f"source/sink, started from the BRIE shoreline after the first year of the adopted {START_YEAR}-2025 "
        f"edgeBE no-groin run and run {a.years} years. One groin between GIS 5 and GIS 6 (red), with "
        f"f$_{{block}}$ = {sc.BLOCKING_FRACTION} throughout: no failure schedule and no fit to the observed "
        "fillet. Every value is the groin run minus a no-groin run from the same start. At the start the "
        f"coast at the groin faces {face0:+.0f} degrees and the updrift approach {up0:+.0f} degrees in "
        "BRIE's frame, and the drift runs south, so GIS 6 is updrift. (a) Shoreline change against no groin "
        f"along GIS {g0}-{g1} (GIS 1 at Cape Point) in years {PROFILE_YEARS[0]} (dashed) and "
        f"{PROFILE_YEARS[1]} (solid). (b) The step across the groin: GIS 6 minus GIS 5, seaward positive. "
        "(c) Sand gained across the whole island, summed shoreline metres; zero means the groin neither "
        "creates nor destroys sand. Dotted: the correction run with no groin, which shows what it does on its "
        "own to the shoreline's natural bends inside its window; the corrected curve minus the dotted one is "
        "the groin's own sand balance. (d) How much shoreline each module moves per year at the groin, and "
        f"(dotted) the leak the corrected module returns within {sc.CORRECTION_HALF_WIDTH} cells of the "
        f"groin. **Shows:** after {a.years} years the corrected groin holds GIS 6 "
        f"{last.loc['corrected', 'updrift_gis6_m']:+.0f} m and GIS 5 {last.loc['corrected', 'downdrift_gis5_m']:+.0f} m "
        f"against no groin and changes the island's sand by {last.loc['corrected', 'island_sand_m']:+.0f} m, of "
        f"which {only_sand[-1]:+.0f} m is the correction acting on the natural bend at the groin and "
        f"{groin_part[-1]:+.0f} m is the groin's own; "
        f"the drift groin without the correction holds GIS 6 {last.loc['drift', 'updrift_gis6_m']:+.0f} m but "
        f"changes the island's sand by {last.loc['drift', 'island_sand_m']:+.0f} m; the pinned groin holds "
        f"GIS 6 {last.loc['pinned', 'updrift_gis6_m']:+.0f} m with {last.loc['pinned', 'island_sand_m']:+.0f} m. "
        f"The corrected module moves {fc.module_transfer_m.iloc[1]:.0f} m in the first year and returns "
        f"{fc.leak_correction_m.iloc[1]:.0f} m of leak; by year {a.years}, {fc.module_transfer_m.iloc[-1]:.0f} m "
        f"and {fc.leak_correction_m.iloc[-1]:.0f} m. The leak it returns is larger than the drift it blocks "
        "because Buxton's shoreline already bends at the groin, so the correction also changes BRIE's no-groin "
        f"behaviour there, by up to {np.abs(-(only[-1] - base[-1])).max():.0f} m of shoreline. As volume "
        f"(one {sc.DY_M:.0f} m cell, {t['depth']:.1f} m active depth), the uncorrected drift groin's "
        f"{last.loc['drift', 'island_sand_m']:+.0f} m is "
        f"{last.loc['drift', 'island_sand_m'] * sc.DY_M * t['depth'] / 1e6:+.1f} M m$^3$ over {a.years} years, "
        f"and the pinned groin's {last.loc['pinned', 'island_sand_m'] * sc.DY_M * t['depth'] / 1e6:+.1f} M m$^3$."))
    print(f"wrote {png.relative_to(REPO)}")


if __name__ == "__main__":
    main()
