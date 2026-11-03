"""
One storm on one dune, in cross-section: the water rises past the toe, holds, and the module lowers the dune.

    python collision_storm_profile_gif.py [--width 25]
        ->  figures/collision_7_one_storm_profile.gif   (and its last frame, collision_7_one_storm_profile.png)

The example cell of collision_explained_figures.py (crest 7.3 m, toe 2.4 m NAVD88) behind a beach of the
given width, through the largest 1999 storm of the hindcast storm series (Dennis), with the module's own
functions at the hindcast's settings. The module acts once per storm, after Barrier3D's overwash step: the
storm's whole volume comes off at once. The frames that lower the dune in steps only animate that one step;
the eroded sand leaves the domain (the module does not yet place it).

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-11-03
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib.animation import PillowWriter  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parent))
import collision_explained_figures as E  # noqa: E402  (the example cell, settings, profile, storm series)
from cascade.dune_collision import beach_factor, larson_volume  # noqa: E402
from site_layer.hat_figure_style import C, INK, INK_MUTED, apply_style, caption, figsize, open_frame, save  # noqa: E402

# --- CONFIG ------------------------------------------------------------------
WIDTH = E.WIDTH                 # m, the example cell's dry beach (GIS 10-12: 0-20 m)
YEAR = 1999
RISE, LOWER, FALL = 8, 8, 4     # frames: water rising, the dune lowered, water falling
FPS = 4.0
DPI = 110
HOLD = 8                        # extra frames on the key states (merged into longer frames)
# -----------------------------------------------------------------------------


# Run: the storm, the profile before and after, the frames
def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    ap.add_argument("--width", type=float, default=WIDTH)
    args = ap.parse_args()
    w = args.width

    st = E.storms()
    yr = st[st.calendar_year == YEAR]
    s = yr.loc[yr.r_m.idxmax()]
    after = E.rows_after(s.r_m, s.period, s.duration, w)
    crest0, crest1 = E.BERM + max(E.ROWS), E.BERM + max(after)
    beach = float(beach_factor(np.array([w / 10.0]), E.W0 / 10.0)[0])
    vol = float(larson_volume(E.dam(s.r_m), np.array([E.dam(E.TOE)]), s.period, s.duration, E.CS,
                              np.array([beach]))[0])

    x0, z0, toe_x = E.profile(E.ROWS, w)
    shore_x = toe_x + w

    apply_style()
    fig, ax = plt.subplots(figsize=figsize("double", 0.45), layout="constrained")

    def frame(water, frac, note):
        ax.clear()
        rows = np.asarray(E.ROWS) + frac * (np.asarray(after) - np.asarray(E.ROWS))
        x, z, _ = E.profile(rows, w)
        ax.fill_between(x, z, E.MHW - 2, color=E.C_SAND, lw=0, zorder=1)
        if frac > 0:
            ax.fill_between(x0, z0, z, where=z0 > z + 1e-6, color=E.C_MOD_FILL, lw=0, zorder=1,
                            label=f"removed by the module: {vol:.1f} m³ per m of shore")
        ax.plot(x0, z0, color=INK_MUTED, lw=0.7, ls="--", zorder=2, label="dune before the storm")
        ax.plot(x, z, color=INK, lw=1.0, zorder=3)
        if water is not None:
            # the water: from the sea up the beach to where it meets the profile
            xs = np.linspace(x.min(), x.max(), 600)
            zs = np.interp(xs, x, z)
            seaward = xs >= np.interp(water, z[6:][::-1], x[6:][::-1])
            ax.fill_between(xs, np.where(seaward, zs, water), water, where=seaward, color=E.C_WATER,
                            alpha=0.75, lw=0, zorder=2)
            ax.text(x.max() - 1, water + 0.12, f"water {water:.2f} m", ha="right", va="bottom", fontsize=8,
                    color=C["LATE"])
        ax.plot(toe_x, E.TOE, "o", ms=4, color=INK, zorder=4)
        ax.annotate(f"toe {E.TOE} m", (toe_x, E.TOE), xytext=(toe_x + 3, E.TOE - 1.0), fontsize=7.5,
                    arrowprops=dict(arrowstyle="-", lw=0.6))
        ax.annotate("", (toe_x, E.MHW - 0.6), (shore_x, E.MHW - 0.6),
                    arrowprops=dict(arrowstyle="<->", lw=0.7, color=INK_MUTED))
        ax.text((toe_x + shore_x) / 2, E.MHW - 0.55, f"beach {w:.0f} m", ha="center", va="bottom",
                fontsize=7.5, color=INK_MUTED)
        ax.set(xlim=(-10, x0.max()), ylim=(E.MHW - 1.5, 9.4), xlabel="Cross-shore distance (m, seaward →)",
               ylabel="Elevation (m NAVD88)")
        ax.legend(loc="upper left", frameon=False, fontsize=7.5)
        ax.set_title(note, fontsize=9, loc="center")
        open_frame(ax)

    out = E.OUT / "collision_7_one_storm_profile.gif"
    writer = PillowWriter(fps=FPS)
    with writer.saving(fig, str(out), dpi=DPI):
        frame(None, 0, f"Before the storm: crest {crest0:.1f} m, toe {E.TOE} m, {w:.0f} m of beach")
        for _ in range(HOLD // 2):
            writer.grab_frame()
        writer.grab_frame()
        for i, h in enumerate(np.linspace(E.MHW + 0.3, s.r_m, RISE)):
            if h < E.TOE:
                note = f"{s.StartTime:%d %b %Y}: the water rises, still below the toe: nothing happens"
            else:
                note = "The water passes the toe and reaches the dune face: the collision regime"
            frame(h, 0, note)
            writer.grab_frame()
        frame(s.r_m, 0, f"Peak {s.r_m:.2f} m, held {s.duration:.0f} h (T {s.period:.1f} s): "
              f"{s.r_m - E.TOE:.2f} m of water above the toe, below the {crest0:.1f} m crest")
        for _ in range(HOLD):
            writer.grab_frame()
        for frac in np.linspace(1 / LOWER, 1, LOWER):
            frame(s.r_m, frac, f"Larson: 4 Cs (R − toe)² t / T × beach factor {beach:.2f} = {vol:.1f} m³/m, "
                               f"off the 20 m dune")
            writer.grab_frame()
        for h in np.linspace(s.r_m, E.MHW + 0.3, FALL):
            frame(h, 1, "The storm passes")
            writer.grab_frame()
        frame(None, 1, f"After: crest {crest0:.2f} → {crest1:.2f} m (−{crest0 - crest1:.2f} m); "
                       f"the sand leaves the domain")
        for _ in range(HOLD + 4):
            writer.grab_frame()

    caption(fig, f"**Tests:** nothing new; the last frame of `collision_7_one_storm_profile.gif`, one storm on "
            f"one dune in cross-section. **How to read:** the example cell of the other figures (two 10 m dune "
            f"rows {E.ROWS[0]:.1f} and {E.ROWS[1]:.1f} m above the {E.BERM} m berm, crest {crest0:.1f} m, toe "
            f"{E.TOE} m NAVD88) behind a {w:.0f} m beach, schematic in shape, through the largest {YEAR} storm "
            f"of the hindcast storm series ({s.StartTime:%d %b}, {s.r_m:.2f} m held {s.duration:.0f} h, T "
            f"{s.period:.1f} s), Cs {E.CS:g}, W0 {E.W0:.0f} m. Dashed, the dune before; purple, what the module "
            f"removes. In the GIF the water rises past the toe, holds, the module lowers both dune rows in "
            f"proportion (its one step per storm, animated), and the water falls. **Shows:** {vol:.1f} m³ per m "
            f"of shore, {crest0 - crest1:.2f} m off the crest; the beach factor {beach:.2f} scales Larson's "
            f"volume for the {w:.0f} m beach. The sand is removed from the model, not yet put on the beach.")
    save(fig, E.OUT / "collision_7_one_storm_profile", close=True)
    print(f"{out}\n{s.StartTime:%d %b %Y} {s.r_m:.2f} m {s.duration:.0f} h: {vol:.2f} m3/m, crest {crest0:.2f} -> {crest1:.2f} m")


if __name__ == "__main__":
    main()
