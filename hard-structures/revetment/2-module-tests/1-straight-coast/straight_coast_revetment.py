"""
Straight-coast test of the simple revetment: full CASCADE, a wall vs no wall.

    python straight_coast_revetment.py run          # both cases, in parallel
    python straight_coast_revetment.py run --case wall
    python straight_coast_revetment.py figures

Inputs are the CASCADE test-suite files (tests/test_coupling), copied into each run
folder under output/raw_runs/experiments/revetment/. Tables and figures are
written beside this script. Needs cascade, barrier3d, pandas, matplotlib.

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-10-12
"""

from __future__ import annotations

import argparse
import shutil
import subprocess
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
REPO = next(p for p in HERE.parents if (p / "cascade").is_dir())
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(REPO / "scripts"))
sys.path.insert(0, str(HERE))

# --- CONFIG ------------------------------------------------------------------
STUDY = "2026-10-12-straight-coast"
RUN_ROOT = REPO / "output" / "raw_runs" / "experiments" / "revetment" / STUDY
INPUTS = REPO / "tests" / "test_coupling"
CASES = ("nowall", "wall")

N_DOMAINS = 15
WALL_DOMAINS = [6, 7, 8]          # the middle three, 0-based
BUILD_YEAR = 1                    # first model year
CREST_HEIGHT_M = 1.2              # above the berm; Dmax here is 1.5 m
START_GAP_M = None                # None = Barrier3D's implicit beach width
YEARS = 60
BACKGROUND_EROSION = -1.0         # m/yr, every domain, so the shoreline reaches the wall

# Waves: the CASCADE test defaults, diffusive at 0 deg (option A is not; README)
WAVES = dict(wave_height=1, wave_period=7, wave_asymmetry=0.8,
             wave_angle_high_fraction=0.2)
# -----------------------------------------------------------------------------


# Build a straight-coast CASCADE in its own run folder
def build_cascade(run_dir: Path):
    from cascade.cascade import Cascade

    run_dir.mkdir(parents=True, exist_ok=True)
    for f in INPUTS.iterdir():
        if f.is_file():
            shutil.copy(f, run_dir)
    return Cascade(
        str(run_dir),
        name="straight_coast",
        storm_file="StormSeries_1kyrs_VCR_Berm1pt9m_Slope0pt04_01.npy",
        elevation_file="b3d_pt75_3284yrs_low-elevations.csv",
        dune_file="pathways-dunes.npy",
        parameter_file="ast-barrier3d-parameters.yaml",
        **WAVES,
        sea_level_rise_rate=0.004,
        sea_level_rise_constant=True,
        background_erosion=BACKGROUND_EROSION,
        alongshore_section_count=N_DOMAINS,
        time_step_count=YEARS + 1,
        min_dune_growth_rate=0.55,
        max_dune_growth_rate=0.95,
        num_cores=1,
        roadway_management_module=False,
        alongshore_transport_module=True,
        beach_nourishment_module=False,
        community_economics_module=False,
    )


# One row per domain for the year just run
def domain_rows(cascade):
    rows = []
    for d, b3d in enumerate(cascade.barrier3d):
        t = b3d.time_index - 1
        crest = b3d.DuneDomain[t, :, :].max(axis=1)
        rows.append(dict(
            year=t, domain=d, x_s_m=b3d.x_s * 10,
            crest_mean_m=float(crest.mean()) * 10,
            qow_m3_per_m=float(b3d.QowTS[-1]),
            dune_cells_migrated=int(b3d.ShorelineChangeTS[t]),
            interior_width_m=float(b3d.InteriorWidth_AvgTS[-1]) * 10,
            dune_migration_on=bool(b3d.dune_migration_on),
        ))
    return rows


# Run one case and write its two tables
def run_case(case: str) -> None:
    from revetment_schedule import RevetmentSchedule

    run_dir = RUN_ROOT / case
    cascade = build_cascade(run_dir)
    wall = None
    if case == "wall":
        wall = RevetmentSchedule(WALL_DOMAINS, BUILD_YEAR, CREST_HEIGHT_M, START_GAP_M)

    # Year loop: wall before, CASCADE, wall after, record
    rows = []
    started = time.time()
    for _ in range(YEARS):
        if wall:
            wall.before_update(cascade)
        cascade.update()
        if cascade.b3d_break:
            print(f"{case}: Barrier3D break in year {cascade.barrier3d[0].time_index - 1}")
            break
        if wall:
            wall.after_update(cascade)
        rows.extend(domain_rows(cascade))

    # Write: every domain every year, and the wall's own record
    tables = HERE / "tables"
    tables.mkdir(exist_ok=True)
    pd.DataFrame(rows).to_csv(tables / f"domains_by_year_{case}.csv", index=False)
    if wall:
        pd.DataFrame(wall.records).to_csv(tables / "wall_record.csv", index=False)
    print(f"{case}: {len(rows) // N_DOMAINS} yr in {time.time() - started:.0f} s")


# Both cases at once, one process each
def run_all() -> None:
    procs = [subprocess.Popen([sys.executable, __file__, "run", "--case", c])
             for c in CASES]
    codes = [p.wait() for p in procs]
    if any(codes):
        sys.exit(f"a case failed: {dict(zip(CASES, codes))}")


# Checks and figures from the tables
def figures() -> None:
    import matplotlib.pyplot as plt
    from site_layer.hat_figure_style import C, apply_style, caption, figsize, save

    tables = HERE / "tables"
    base = pd.read_csv(tables / "domains_by_year_nowall.csv")
    wall = pd.read_csv(tables / "domains_by_year_wall.csv")
    rec = pd.read_csv(tables / "wall_record.csv")
    first, last = base.year.min(), base.year.max()

    # Shoreline change per domain (seaward positive) and its total
    def change(df):
        x = df.pivot(index="year", columns="domain", values="x_s_m")
        return -(x.loc[last] - x.loc[first])
    dx_base, dx_wall = change(base), change(wall)

    # Checks, printed and written
    held = wall[wall.domain.isin(WALL_DOMAINS)]
    pinned = rec[rec.pinned_m > 0]
    checks = pd.DataFrame([
        ("dune cells migrated in wall domains", int(held.dune_cells_migrated.abs().sum())),
        ("lowest crest in wall domains (m)", round(rec.crest_min_m.min(), 3)),
        ("crest floor (m)", CREST_HEIGHT_M),
        ("first year the shoreline is held", int(pinned.year.min()) if len(pinned) else "none"),
        ("shoreline put back by the wall, summed (m)", round(rec.pinned_m.sum(), 2)),
        ("sand added per wall domain, mean (m3/m)",
            round(rec.groupby("domain").sand_added_m3_per_m.sum().mean(), 1)),
        ("sand added per wall domain, last 10 yr mean (m3/m/yr)", round(
            rec[rec.year > last - 10].groupby("domain").sand_added_m3_per_m.sum().mean() / 10, 2)),
        ("overwash in wall domains, no wall (m3/m)",
            round(base[base.domain.isin(WALL_DOMAINS)].qow_m3_per_m.sum(), 1)),
        ("overwash in wall domains, wall (m3/m)", round(held.qow_m3_per_m.sum(), 1)),
        ("summed shoreline change, no wall (m)", round(dx_base.sum(), 2)),
        ("summed shoreline change, wall (m)", round(dx_wall.sum(), 2)),
        ("wall-domain change vs no wall, mean (m)", round(
            (dx_wall - dx_base)[WALL_DOMAINS].mean(), 2)),
        ("neighbour change vs no wall, max (m)", round(
            (dx_wall - dx_base).drop(WALL_DOMAINS).max(), 2)),
        ("neighbour change vs no wall, min (m)", round(
            (dx_wall - dx_base).drop(WALL_DOMAINS).min(), 2)),
    ], columns=["check", "value"])
    checks.to_csv(tables / "checks.csv", index=False)
    print(checks.to_string(index=False))

    apply_style()
    fig, axes = plt.subplots(2, 3, figsize=figsize("double", aspect=0.68))
    shade = dict(color=C["BASE_FILL"], alpha=0.5, lw=0)
    mid = WALL_DOMAINS[len(WALL_DOMAINS) // 2]

    # (a) gap through time
    ax = axes[0, 0]
    for d, g in rec.groupby("domain"):
        ax.plot(g.year, g.gap_m, color=C["ACCENT"], lw=1)
    ax.axhline(0, color=C["INK_MUTED"], lw=0.6)
    ax.set(xlabel="model year", ylabel="shoreline-to-wall gap (m)",
           title="a  The gap closes, then holds at 0")

    # (b) sand the wall needs to hold the shoreline
    ax = axes[0, 1]
    for d, g in rec.groupby("domain"):
        ax.plot(g.year, g.sand_added_m3_per_m.cumsum(), color=C["ADDED"], lw=1)
    ax.set(xlabel="model year", ylabel="cumulative sand added (m³/m)",
           title="b  Sand to hold the shoreline")

    # (c) dune crest in a wall domain vs the same domain with no wall
    ax = axes[0, 2]
    for df, col, lab in ((base, C["BASE"], "no wall"), (wall, C["ACCENT"], "wall")):
        g = df[df.domain == mid]
        ax.plot(g.year, g.crest_mean_m, color=col, lw=1, label=lab)
    ax.axhline(CREST_HEIGHT_M, color=C["INK_MUTED"], lw=0.6, ls="--")
    ax.set(xlabel="model year", ylabel="mean dune crest above berm (m)",
           title=f"c  Dune crest, domain {mid}")
    ax.legend(frameon=False)

    # (d) cumulative overwash in the wall domains
    ax = axes[1, 0]
    for df, col, lab in ((base, C["BASE"], "no wall"), (wall, C["ACCENT"], "wall")):
        q = df[df.domain.isin(WALL_DOMAINS)].groupby("year").qow_m3_per_m.mean().cumsum()
        ax.plot(q.index, q.values, color=col, lw=1, label=lab)
    ax.set(xlabel="model year", ylabel="cumulative overwash (m³/m)",
           title="d  Overwash in the wall domains")
    ax.legend(frameon=False)

    # (e) shoreline relative to no wall through time, wall domain and neighbours
    ax = axes[1, 1]
    xb = base.pivot(index="year", columns="domain", values="x_s_m")
    xw = wall.pivot(index="year", columns="domain", values="x_s_m")
    ramp = {mid: C["ACCENT"], min(WALL_DOMAINS) - 1: C["LATE"],
            min(WALL_DOMAINS) - 3: C["EARLY"], 0: C["BASE"]}
    for d, col in ramp.items():
        ax.plot(xb.index, -(xw[d] - xb[d]), color=col, lw=1, label=f"domain {d}")
    ax.axhline(0, color=C["INK_MUTED"], lw=0.6)
    ax.set(xlabel="model year", ylabel="seaward of no wall (m)",
           title="e  Shoreline vs no wall")
    ax.legend(frameon=False, fontsize=6.5)

    # (f) shoreline change along the coast
    ax = axes[1, 2]
    ax.axvspan(min(WALL_DOMAINS) - 0.5, max(WALL_DOMAINS) + 0.5, **shade)
    ax.plot(dx_base.index, dx_base.values, "o-", color=C["BASE"], ms=3, lw=1, label="no wall")
    ax.plot(dx_wall.index, dx_wall.values, "o-", color=C["ACCENT"], ms=3, lw=1, label="wall")
    ax.set(xlabel="domain (wall shaded)", ylabel=f"shoreline change {first}–{last} (m)",
           title="f  Change along the coast")
    ax.legend(frameon=False, loc="lower center", bbox_to_anchor=(0.5, 0.06))

    fig.tight_layout()
    caption(fig, (
        f"The simple revetment on a straight coast: {N_DOMAINS} domains, a wall on "
        f"domains {WALL_DOMAINS[0]}-{WALL_DOMAINS[-1]} from year {BUILD_YEAR} (crest "
        f"floor {CREST_HEIGHT_M} m above the berm, dashed), {BACKGROUND_EROSION} m/yr "
        "background erosion everywhere, CASCADE test waves (Hs 1 m). (a) The "
        "shoreline-to-wall gap in each wall domain; the wall is repaired, so once the "
        "gap reaches 0 the shoreline is held at the wall. (b) The sand that takes, "
        "added through CASCADE's nourishment path. (c) Mean dune crest in the middle "
        "wall domain. (d) Cumulative overwash, mean of the wall domains. (e) Shoreline "
        "position relative to the no-wall run in the middle wall domain, a neighbour "
        "and two domains further away. (f) Shoreline change over the run, seaward "
        "positive."))
    out = save(fig, HERE / "straight_coast_revetment.png", close=True)
    print("wrote", out[0])


# Run or plot
def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("step", choices=["run", "figures"])
    p.add_argument("--case", choices=CASES)
    a = p.parse_args()
    if a.step == "figures":
        figures()
    elif a.case:
        run_case(a.case)
    else:
        run_all()


if __name__ == "__main__":
    main()
