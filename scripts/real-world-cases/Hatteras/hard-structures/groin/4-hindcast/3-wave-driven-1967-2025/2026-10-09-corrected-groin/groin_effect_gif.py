#!/usr/bin/env python3
"""
The groin's effect on the shoreline, year by year: the model with and without the corrected groin, on a map.

    python groin_effect_gif.py   ->  figures/corrected_groin_effect.gif

Left: the 1.5 km true-scale window on the Buxton groins (UTM 18N, north up), the island filled to
the observed CoastSat shoreline, and the modelled shoreline without a groin (grey) and with the
two-number corrected groin (0.12, off from 2017; black), placed as plan_view_figure.py places them.
The strip between the two is shaded: sand where the groin holds the shoreline seaward, blue where
it pushes it landward. Right: the gap across the groins for both runs and the observations. The
calibration run (1996-2008) plays first, then the test run (2009-2024), which restarts from the
observed shoreline. plan_view_figure.py and results_figures.py are imported unchanged.

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-10-09
"""
from __future__ import annotations

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from matplotlib.animation import PillowWriter  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402
from matplotlib.patches import Patch  # noqa: E402
from shapely.geometry import box  # noqa: E402

import plan_view_figure as P  # noqa: E402
from site_layer.hat_figure_style import C, INK, INK_MUTED, MAP_TEXT_DARK, apply_style, figsize, open_frame  # noqa: E402

rf, pv, SR = P.rf, P.pv, P.SR

# --- CONFIG ------------------------------------------------------------------
OUT = rf.FIGS / "corrected_groin_effect.gif"
NONE, GROIN = "no groin", f"corrected, {rf.CHOSEN} to 2016"
HELD_C, CUT_C = C["ADDED"], C["LATE"]       # the strip: the groin holds sand seaward / pushes landward
FPS = 2
HOLD = 6
# -----------------------------------------------------------------------------


# Per-transect modelled positions (map x, y, and chainage) in year y, in alongshore order
def positions(data, path, period, y, start, dirs):
    x = np.load(path)
    k = y - period
    seaward = -(0.5 * (x[k] + x[k + 1]) - x[0])
    rows = []
    for tid, s0 in start.items():
        c = s0 + seaward[P.PAD + int(data.dom.loc[tid]) - 1]
        o, dx, dy = dirs.loc[tid]
        rows.append((tid, o[0] + c * dx, o[1] + c * dy, c, data.order.loc[tid]))
    return pd.DataFrame(rows, columns=["transect_id", "px", "py", "chain", "o"]).sort_values("o")


# Shade between the two runs, transect to transect: sand where the groin run is seaward, blue where landward
def strip(ax, a, b):
    for i in range(len(a) - 1):
        if abs(a["o"].iloc[i + 1] - a["o"].iloc[i]) > pv.GAP_M:
            continue
        d = (b["chain"].iloc[i] + b["chain"].iloc[i + 1]) - (a["chain"].iloc[i] + a["chain"].iloc[i + 1])
        quad = np.array([[a["px"].iloc[i], a["py"].iloc[i]], [a["px"].iloc[i + 1], a["py"].iloc[i + 1]],
                         [b["px"].iloc[i + 1], b["py"].iloc[i + 1]], [b["px"].iloc[i], b["py"].iloc[i]]])
        ax.fill(quad[:, 0], quad[:, 1], color=HELD_C if d > 0 else CUT_C, alpha=0.75, lw=0, zorder=3.3)


def main():
    apply_style()
    data = pv.Data()
    cx, cy = np.vstack(data.groins).mean(axis=0)
    x0, x1, y0, y1 = cx - pv.HALF_M, cx + pv.HALF_M, cy - pv.HALF_M, cy + pv.HALF_M
    ref = data.ref[data.ref["px"].between(x0 - 400, x1 + 400) & data.ref["py"].between(y0 - 400, y1 + 400)]
    origins = data.ann.drop_duplicates("transect_id").set_index("transect_id")[["origin_x", "origin_y", "dir_x", "dir_y"]]
    keep = [t for t in ref["transect_id"] if t in origins.index]
    dirs = pd.DataFrame([((origins.loc[t, "origin_x"], origins.loc[t, "origin_y"]), origins.loc[t, "dir_x"],
                          origins.loc[t, "dir_y"]) for t in keep], index=keep, columns=["o", "dx", "dy"])
    starts = {p: P.start_positions(data, p, dirs.index).dropna() for p in (rf.CAL, rf.TEST)}
    files = {p: rf.run_files(p) for p in (rf.CAL, rf.TEST)}
    frames = [(y, rf.CAL) for y in range(1996, 2009)] + [(y, rf.TEST) for y in range(2009, 2025)]
    gaps = {p: {lab: SR.model_change(files[p][lab]) for lab in (NONE, GROIN)} for p in (rf.CAL, rf.TEST)}
    obs = {p: SR.coastsat_change(p) for p in (rf.CAL, rf.TEST)}
    glim = (min(min(g.min() for g in gaps[p].values()) for p in gaps) - 8,
            max(max(g.max() for g in gaps[p].values()) for p in gaps) + 8)

    fig = plt.figure(figsize=figsize("double", height=4.6))
    gs = fig.add_gridspec(1, 2, width_ratios=(1.15, 1), left=0.01, right=0.97, top=0.88, bottom=0.2, wspace=0.18)
    ax, axg = fig.add_subplot(gs[0]), fig.add_subplot(gs[1])
    handles = [Line2D([], [], color=rf.OBS_C, lw=P.OBS_LW, label="observed shoreline (CoastSat)"),
               Line2D([], [], color=C["BASE"], lw=1.6, label="model, no groin"),
               Line2D([], [], color=INK, lw=1.6, label=f"model, corrected groin ({rf.CHOSEN}, off from 2017)"),
               Patch(facecolor=HELD_C, alpha=0.75, label="groin holds the shoreline seaward"),
               Patch(facecolor=CUT_C, alpha=0.75, label="groin pushes it landward"),
               Line2D([], [], color=C["GROIN"], lw=2.0, label="Buxton groins")]
    fig.legend(handles=handles, loc="lower center", ncol=3, fontsize=6.6, frameon=False)

    def draw(y, period):
        ax.clear()
        ax.set_facecolor(pv.OCEAN)
        pv.fill(ax, data.lands[y].intersection(box(x0 - 50, y0 - 50, x1 + 50, y1 + 50)), pv.SAND)
        a = positions(data, files[period][NONE], period, y, starts[period], dirs)
        b = positions(data, files[period][GROIN], period, y, starts[period], dirs)
        strip(ax, a, b)
        for df, col in ((a, C["BASE"]), (b, INK)):
            ax.plot(*pv.shoreline(df, data.order), color=col, lw=1.6, zorder=3.6)
        ax.plot(*data.lines[y], color=rf.OBS_C, lw=P.OBS_LW, zorder=3.5)
        pv.map_furniture(ax, data, x0, x1, y0, y1, 500, "plain", MAP_TEXT_DARK, segments=1)
        off = {lab: P.domain_offsets(data, files[period][lab], period, y, starts[period]) for lab in (NONE, GROIN)}
        held6 = (b.assign(g=b["transect_id"].map(data.dom)).query("g == 6")["chain"].mean()
                 - a.assign(g=a["transect_id"].map(data.dom)).query("g == 6")["chain"].mean())
        held5 = (b.assign(g=b["transect_id"].map(data.dom)).query("g == 5")["chain"].mean()
                 - a.assign(g=a["transect_id"].map(data.dom)).query("g == 5")["chain"].mean())
        ax.text(0.97, 0.36, f"groin's effect\nGIS 6: {held6:+.0f} m\nGIS 5: {held5:+.0f} m", transform=ax.transAxes,
                ha="right", va="bottom", fontsize=7, color=INK, zorder=9,
                bbox=dict(facecolor="white", alpha=0.85, edgecolor="none", boxstyle="round,pad=0.3"))
        del off
        axg.clear()
        for p, tag in ((rf.CAL, "calibration run"), (rf.TEST, "test run")):
            yrs = p + np.arange(len(gaps[p][NONE]))
            upto = yrs <= y if p == period else (yrs <= 9999 if p < period else yrs < 0)
            for lab, col in ((NONE, C["BASE"]), (GROIN, INK)):
                axg.plot(yrs + 0.5, gaps[p][lab], color="0.88", lw=1.0, zorder=1)
                axg.plot(yrs[upto] + 0.5, gaps[p][lab][upto], color=col, lw=1.8, zorder=2)
            o = obs[p][(obs[p].index >= p) & (obs[p].index < p + len(yrs))]
            shown = o[o.index <= y] if p == period else (o if p < period else o.iloc[:0])
            axg.plot(shown.index + 0.5, shown.values, "o-", color=rf.OBS_C, lw=0.9, ms=3, mfc="white", zorder=3)
        axg.axvline(2009, color=INK_MUTED, lw=0.7, ls=(0, (1, 2)))
        axg.text(2009.3, glim[1] - 3, "test run restarts\nfrom observed", fontsize=6, color=INK_MUTED, va="top")
        axg.axvline(2017.5, color=INK_MUTED, lw=0.7, ls=(0, (3, 2)))
        axg.text(2017.8, glim[0] + 3, "2017 fill;\ngroin off", fontsize=6, color=INK_MUTED, va="bottom")
        axg.axvline(y + 0.5, color=C["GROIN"], lw=0.8)
        axg.axhline(0, color=INK_MUTED, lw=0.5)
        axg.set_xlim(1996, 2025)
        axg.set_ylim(*glim)
        axg.set_xlabel("Year")
        axg.set_ylabel("Change in the gap across the groins\nsince each run's start (m)", fontsize=7.5)
        axg.grid(axis="y")
        open_frame(axg)
        axg.set_title("The gap: groin (black) against no groin (grey)", fontsize=8)
        run = "calibration run, from 1996" if period == rf.CAL else "test run, from 2009"
        state = ("groin on" if y < 2017 else "groin off (2017 fill)")
        fig.suptitle(f"{y}: the shoreline with and without the corrected groin ({run}; {state})", fontsize=9.5)

    writer = PillowWriter(fps=FPS)
    with writer.saving(fig, str(OUT), dpi=120):
        for y, period in frames:
            draw(y, period)
            writer.grab_frame()
        for _ in range(HOLD):
            writer.grab_frame()
    plt.close(fig)
    print("wrote", OUT.relative_to(rf.cgm.REPO))


if __name__ == "__main__":
    main()
