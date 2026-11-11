"""
Storms on one dune cell as the model holds it: 10 m columns, the dune's two rows, the domain toe, the water.

    python collision_columns_gif.py [--gis 12] [--cell 25]
        ->  figures/collision_9_columns.gif   (and its last frame, collision_9_columns.png)

The cell's own Barrier3D inputs at the 1996 start (dune-topo v3: the dune height above the berm, both dune rows
alike as Barrier3D starts them, and the first 15 interior rows), its domain's one toe (dune-toes,
dune_toe_cells.csv), the berm (1.7 m NAVD88), and two storms of the 1997-2009 storm series: Dennis (the largest of
1999) and Isabel (the largest of 2003). Each storm goes through the module's own functions at the hindcast's
settings (Larson alone, Cs 1.9e-5, a 20 m dune): collision_cells, larson_volume, lower_dunes. The module acts once
per storm, after Barrier3D's overwash step; the frames that lower the dune only animate that one step. The beach
and shoreface are not on Barrier3D's grid: they are drawn as the berm in front of the dune. Dune growth between the
storms is left out, so the second storm meets the dune the first one left.

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-11-04
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
from matplotlib.animation import PillowWriter  # noqa: E402

_HERE = Path(__file__).resolve()
_HAT_DIR = next(_p for _p in _HERE.parents if (_p / "hatteras_paths.py").exists())
_REPO = next(_p for _p in _HERE.parents if (_p / "pyproject.toml").exists())
for _p in (_HAT_DIR, _REPO):
    if str(_p) not in sys.path:
        sys.path.insert(0, str(_p))
from hatteras_paths import WORKSPACE, workspace_twin  # noqa: E402
from cascade.dune_collision import collision_cells, larson_volume, lower_dunes  # noqa: E402
from site_layer import hat_topo_version as tv  # noqa: E402
from site_layer.hat_figure_style import C, INK, INK_MUTED, apply_style, caption, figsize, open_frame, save  # noqa: E402

# --- CONFIG ------------------------------------------------------------------
OUT = workspace_twin(_HERE.parent) / "figures"
START = "1996-start"
STORMS = (WORKSPACE / "data" / "hatteras_init" / "3-env-forcings" / "3-storms" / "hindcast_storms"
          / "1997_2009" / "1997_2009_storms_v3_240_summary.csv")
STORM_YEARS = ((1999, "Dennis"), (2003, "Isabel"))   # the largest storm of each year
MHW = 0.36                      # m NAVD88
BERM = 1.7                      # m NAVD88, Barrier3D's berm on the Hatteras domains
CS = 1.9e-5                     # the hindcast's Cs (Larson alone)
DUNE_W = 20.0                   # m, the two 10 m dune rows
RESTART = 0.075                 # m, Barrier3D's smallest dune height
BEACH_X = 30.0                  # m of berm drawn in front of the dune (not on the grid)
INTERIOR_ROWS = 15              # interior columns drawn behind the dune
BASE = -0.5                     # m NAVD88, the bottom of the columns as drawn
C_SAND, C_DUNE, C_BEACH = "#e8dcc0", "#d9c69a", "#f1ead8"
C_CUT, C_TOE = "#6A51A3", "#B2182B"
RISE, LOWER, FALL, HOLD = 7, 8, 4, 8
FPS, DPI = 4.0, 110
# -----------------------------------------------------------------------------


def dam(z):
    return (np.asarray(z, dtype=float) - MHW) / 10.0


def inputs(gis, cell):
    v = tv.product_dir(START) / "dune-topo" / "v3"
    dune = float(np.load(v / "dunes" / f"domain_{gis}_dune.npy")[cell]) * 10.0       # m above the berm
    topo = np.load(v / "topography" / f"domain_{gis}_topography.npy")[:INTERIOR_ROWS, cell] * 10.0 + MHW
    toe = float(pd.read_csv(tv.product_dir(START) / "dune-toes" / "dune_toe_cells.csv")
                .groupby("gis").toe_m_navd88.first()[gis])
    s = pd.read_csv(STORMS, parse_dates=["StartTime"]).assign(r_m=lambda d: d.Rhigh * 10.0 + MHW)
    storms = []
    for year, name in STORM_YEARS:
        y = s[s.StartTime.dt.year == year]
        storms.append((name, y.loc[y.r_m.idxmax()]))
    return dune, topo, toe, storms


def one_storm(rows_m, toe, st):
    """The module on one cell: rows (m above the berm) after the storm, the volume (m3/m), in collision?"""
    rows = np.asarray(rows_m, float)[None, :] / 10.0
    crest = rows.max(axis=1) + dam(BERM)
    t = np.array([dam(toe)])
    r = float(dam(st.r_m))
    hit = bool(collision_cells(crest, crest, t, r)[0])
    if not hit:
        return np.asarray(rows_m, float), 0.0, False
    v = larson_volume(r, t, st.period, st.duration, CS)
    new = lower_dunes(rows, v, t, dam(BERM), DUNE_W / 10.0, RESTART / 10.0)[0] * 10.0
    return new, float(v[0]), True


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    ap.add_argument("--gis", type=int, default=12)
    ap.add_argument("--cell", type=int, default=25)
    a = ap.parse_args()
    dune, topo, toe, storms = inputs(a.gis, a.cell)
    rows0 = np.array([dune, dune])
    seq, rows = [], rows0
    for name, st in storms:
        new, vol, hit = one_storm(rows, toe, st)
        seq.append((name, st, rows, new, vol, hit))
        rows = new
    crest0 = BERM + rows0.max()

    apply_style()
    fig, ax = plt.subplots(figsize=figsize("double", 0.45), layout="constrained")
    xs_int = 20.0 + 10.0 * np.arange(INTERIOR_ROWS)
    x_max = 20.0 + 10.0 * INTERIOR_ROWS

    def frame(rows_now, rows_before, water, title, cut_label=None):
        ax.clear()
        # the berm in front of the dune: the beach, not on the grid
        ax.add_patch(plt.Rectangle((-BEACH_X, BASE), BEACH_X, BERM - BASE, facecolor=C_BEACH, edgecolor="none",
                                   hatch="///", lw=0, zorder=1))
        ax.text(-BEACH_X / 2, (BASE + BERM) / 2, "berm 1.7 m\n(beach: not on the grid)", ha="center", va="center",
                fontsize=7, color=INK_MUTED, zorder=6)
        # the two dune rows: 10 m columns from the berm to the dune top
        for k, (h_now, h_before) in enumerate(zip(rows_now, rows_before)):
            x0 = 10.0 * k
            ax.add_patch(plt.Rectangle((x0, BASE), 10.0, BERM - BASE, facecolor=C_SAND, edgecolor=INK_MUTED,
                                       lw=0.5, zorder=2))
            ax.add_patch(plt.Rectangle((x0, BERM), 10.0, h_now, facecolor=C_DUNE, edgecolor=INK, lw=0.9, zorder=3))
            if h_before > h_now + 1e-6:
                ax.add_patch(plt.Rectangle((x0, BERM + h_now), 10.0, h_before - h_now, facecolor=C_CUT, alpha=0.55,
                                           edgecolor=C_CUT, lw=0.8, zorder=3,
                                           label=cut_label if k == 0 else None))
        ax.plot([0, 20], [BERM + rows0.max()] * 2, color=INK_MUTED, lw=0.8, ls="--", zorder=4,
                label="dune before the storms")
        # the interior columns
        for x0, z in zip(xs_int, topo):
            ax.add_patch(plt.Rectangle((x0, BASE), 10.0, z - BASE, facecolor=C_SAND, edgecolor=INK_MUTED, lw=0.5,
                                       zorder=2))
        # the toe: one number for the domain, a level across the dune
        ax.plot([-BEACH_X, 20], [toe, toe], color=C_TOE, lw=1.1, ls=(0, (4, 2)), zorder=5)
        ax.text(-BEACH_X + 1, toe + 0.1, f"GIS {a.gis} toe {toe:.2f} m", fontsize=7.5, color=C_TOE, zorder=6)
        if water is not None:
            ax.fill_between([-BEACH_X, 0], BASE, water, color=C["WATER"], alpha=0.7, lw=0, zorder=4)
            ax.text(-BEACH_X + 1, water + 0.1 if water < toe - 0.4 or water > toe + 0.5 else water - 0.35,
                    f"water {water:.2f} m", fontsize=7.5, color=C["LATE"], zorder=6)
        ax.axvline(0, color=INK, lw=0.4, zorder=0)
        ax.text(10, -0.35, "dune\n2 × 10 m", ha="center", va="bottom", fontsize=6.5, color=INK_MUTED, zorder=6)
        ax.text((20 + x_max) / 2, -0.35, "interior, 10 m columns", ha="center", va="bottom", fontsize=6.5,
                color=INK_MUTED, zorder=6)
        ax.set(xlim=(-BEACH_X, x_max), ylim=(BASE, crest0 + 1.6),
               xlabel="Cross-shore distance landward of the dune front (m)", ylabel="Elevation (m NAVD88)")
        ax.legend(loc="upper right", frameon=False, fontsize=7.5)
        ax.set_title(title, fontsize=9, loc="center")
        open_frame(ax)

    out = OUT / "collision_9_columns.gif"
    writer = PillowWriter(fps=FPS)
    total = 0.0
    with writer.saving(fig, str(out), dpi=DPI):
        frame(rows0, rows0, None, f"GIS {a.gis}, cell {a.cell}, as the model holds it: dune top {crest0:.2f} m, "
                                  f"toe {toe:.2f} m")
        for _ in range(HOLD):
            writer.grab_frame()
        for name, st, before, after, vol, hit in seq:
            label = f"removed by the module ({name}: {vol:.1f} m³ per m of shore)"
            for h in np.linspace(MHW + 0.2, st.r_m, RISE):
                note = (f"{name} ({st.StartTime:%b %Y}): the water rises, below the toe: nothing happens"
                        if h < toe else f"{name}: the water passes the toe and meets the dune: collision")
                frame(before, before, h, note)
                writer.grab_frame()
            above = st.r_m - toe
            frame(before, before, st.r_m, f"{name} peak {st.r_m:.2f} m, held {st.duration:.0f} h (T {st.period:.1f} s): "
                                           f"{above:.2f} m above the toe, below the {BERM + before.max():.2f} m top")
            for _ in range(HOLD):
                writer.grab_frame()
            if hit:
                for f in np.linspace(1 / LOWER, 1, LOWER):
                    mid = before + f * (after - before)
                    frame(mid, before, st.r_m, f"Larson: 4 × {CS:g} × ({above:.2f} m)² × {st.duration:.0f} h / "
                                               f"{st.period:.1f} s = {vol:.1f} m³/m, off the 20 m dune", label)
                    writer.grab_frame()
            for h in np.linspace(st.r_m, MHW + 0.2, FALL):
                frame(after, before, h, f"{name} passes", label)
                writer.grab_frame()
            drop = before.max() - after.max()
            total += drop
            frame(after, before, None, f"After {name}: dune top {BERM + before.max():.2f} → {BERM + after.max():.2f} m "
                                       f"(−{drop:.2f} m)", label)
            for _ in range(HOLD):
                writer.grab_frame()
        last = seq[-1][3]
        frame(last, rows0, None, f"Both storms: dune top {crest0:.2f} → {BERM + last.max():.2f} m (−{total:.2f} m); "
                                 f"the sand leaves the model", "removed by the module, both storms")
        for _ in range(HOLD + 4):
            writer.grab_frame()

    d, i = seq[0], seq[1]
    caption(fig, (
        f"**How to read:** one dune cell as Barrier3D holds it (GIS {a.gis}, cell {a.cell}, the 1996 start): 10 m "
        f"columns, the dune as two 10 m rows {rows0[0]:.2f} m above the {BERM} m berm, the interior behind it, the "
        f"beach in front only as the berm (not on the grid), and the domain's one toe ({toe:.2f} m NAVD88, dashed "
        f"red). In the GIF the two largest storms of 1999 and 2003 come in turn; between toe and dune top the module "
        f"takes Larson's volume (Cs {CS:g}) off both rows. **Shows:** {d[0]} ({d[1].r_m:.2f} m, "
        f"{d[1].r_m - toe:.2f} m above the toe, {d[1].duration:.0f} h) takes {d[4]:.1f} m³/m, "
        f"{d[2].max() - d[3].max():.2f} m off the top; {i[0]} ({i[1].r_m:.2f} m, {i[1].r_m - toe:.2f} m above the "
        f"toe, {i[1].duration:.0f} h) takes {i[4]:.1f} m³/m, {i[2].max() - i[3].max():.2f} m. Erosion goes with the "
        f"square of the water above the toe. Dune growth between the storms is left out."))
    save(fig, OUT / "collision_9_columns", close=True)
    print(out)
    for name, st, before, after, vol, hit in seq:
        print(f"{name}: {st.r_m:.2f} m, {st.r_m - toe:.2f} above toe, {st.duration:.0f} h, T {st.period:.1f}: "
              f"{vol:.2f} m3/m, top {BERM + before.max():.2f} -> {BERM + after.max():.2f}")


if __name__ == "__main__":
    main()
