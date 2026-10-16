#!/usr/bin/env python3
"""
The corrected groin on and off, year by year, zoomed on Buxton: one GIF per period, calendar-year setup.

    python groin_on_off_gif.py   ->  figures/corrected_groin_on_off_<calibration|test>.gif

Left: the 1.5 km true-scale window on the groins (UTM 18N, north up), the island filled to the observed
CoastSat shoreline, and the modelled shoreline without a groin (grey) and with the corrected groin
(black); the strip between them is shaded where the groin moves the shoreline. Right: shoreline change
since the period's start along GIS 2-12 for both runs and CoastSat, the groin face between GIS 5 and 6.
Calibration: 1996-2009 at strength 0.16. Test: 2009-2025 at 0.10, off from the 2017 fill. Row k of a
run is calendar year start + k; the start line is the start year's calendar-mean CoastSat shoreline.

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-10-13
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
from site_layer.hatteras_site_config import last_model_year  # noqa: E402

rf, pv, SR, cgm = P.rf, P.pv, P.SR, P.rf.cgm

# --- CONFIG ------------------------------------------------------------------
RUNS = {rf.CAL: ("calibration", 0.16, "1996-2009, strength 0.16"),
        rf.TEST: ("test", 0.10, "2009-2025, strength 0.10, off from the 2017 fill")}
PROFILE_GIS = range(2, 13)                  # profile domains, south to north; GIS 1 left out (the solved end moves it 150 m)
HELD_C, CUT_C = C["ADDED"], C["LATE"]       # groin run seaward / landward of no groin
FPS = 2
HOLD = 6
# -----------------------------------------------------------------------------


# Seaward change since the start in every model domain for calendar year y (landward + in the matrix)
def seaward(path, period, y):
    x = np.load(path)
    return -(x[y - period] - x[0])


# Per-transect modelled positions (map x, y, chainage) in year y, in alongshore order
def positions(data, path, period, y, start, dirs):
    sw = seaward(path, period, y)
    rows = []
    for tid, s0 in start.items():
        c = s0 + sw[P.PAD + int(data.dom.loc[tid]) - 1]
        o, dx, dy = dirs.loc[tid]
        rows.append((tid, o[0] + c * dx, o[1] + c * dy, c, data.order.loc[tid]))
    return pd.DataFrame(rows, columns=["transect_id", "px", "py", "chain", "o"]).sort_values("o")


# Shade between the two runs, transect to transect: sand where the groin run is seaward, blue where landward
def strip(ax, a, b):
    for i in range(len(a) - 1):
        if abs(a["o"].iloc[i + 1] - a["o"].iloc[i]) > pv.GAP_M:
            continue
        d = (b["chain"].iloc[i] + b["chain"].iloc[i + 1]) - (a["chain"].iloc[i] + a["chain"].iloc[i + 1])
        if abs(d) < 1.0:
            continue
        quad = np.array([[a["px"].iloc[i], a["py"].iloc[i]], [a["px"].iloc[i + 1], a["py"].iloc[i + 1]],
                         [b["px"].iloc[i + 1], b["py"].iloc[i + 1]], [b["px"].iloc[i], b["py"].iloc[i]]])
        ax.fill(quad[:, 0], quad[:, 1], color=HELD_C if d > 0 else CUT_C, alpha=0.75, lw=0, zorder=3.3)


# Observed change since the start per GIS domain in year y (mean over the domain's transects), seaward +
def observed_profile(data, y, start):
    ann = data.ann[data.ann["year"] == y].set_index("transect_id")["median"]
    t = [t for t in start.index if t in ann.index]
    d = pd.Series(ann.loc[t].to_numpy() - start.loc[t].to_numpy(), index=t)
    return d.groupby(data.dom.loc[t].astype(int).to_numpy()).mean()


def make(period, data, frame):
    name, s, label = RUNS[period]
    cgm.observed_change(period)                               # sets the calendar start window
    last = last_model_year(period)
    files = {"none": SR.baseline_file(period), "groin": rf.two_number_file(period, s)}
    x0, x1, y0, y1, dirs = frame
    start_map = P.start_positions(data, period, dirs.index).dropna()
    prof_t = data.dom[data.dom.isin(PROFILE_GIS)].index
    start_prof = P.start_positions(data, period, prof_t).dropna()
    years = list(range(period, last + 1))
    idx = [P.PAD + g - 1 for g in PROFILE_GIS]
    prof = {k: np.array([seaward(p, period, y)[idx] for y in years]) for k, p in files.items()}
    obs = {y: observed_profile(data, y, start_prof) for y in years}
    lo = min(min(v.min() for v in prof.values()), min(o.min() for o in obs.values())) - 10
    hi = max(max(v.max() for v in prof.values()), max(o.max() for o in obs.values())) + 10

    fig = plt.figure(figsize=figsize("double", height=4.6))
    gs = fig.add_gridspec(1, 2, width_ratios=(1.0, 1.15), left=0.01, right=0.97, top=0.88, bottom=0.22, wspace=0.16)
    ax, axp = fig.add_subplot(gs[0]), fig.add_subplot(gs[1])
    handles = [Line2D([], [], color=rf.OBS_C, lw=P.OBS_LW, label="observed (CoastSat)"),
               Line2D([], [], color=C["BASE"], lw=1.6, label="model, groin off (no groin)"),
               Line2D([], [], color=INK, lw=1.6, label=f"model, groin on ({label.split(', ', 1)[1]})"),
               Patch(facecolor=HELD_C, alpha=0.75, label="groin holds the shoreline seaward"),
               Patch(facecolor=CUT_C, alpha=0.75, label="groin pushes it landward"),
               Line2D([], [], color=C["GROIN"], lw=2.0, label="Buxton groins")]
    fig.legend(handles=handles, loc="lower center", ncol=3, fontsize=6.6, frameon=False)
    gx = np.array(list(PROFILE_GIS))

    def draw(k, y):
        ax.clear()
        ax.set_facecolor(pv.OCEAN)
        pv.fill(ax, data.lands[y].intersection(box(x0 - 50, y0 - 50, x1 + 50, y1 + 50)), pv.SAND)
        a = positions(data, files["none"], period, y, start_map, dirs)
        b = positions(data, files["groin"], period, y, start_map, dirs)
        strip(ax, a, b)
        for df, col in ((a, C["BASE"]), (b, INK)):
            ax.plot(*pv.shoreline(df, data.order), color=col, lw=1.6, zorder=3.6)
        ax.plot(*data.lines[y], color=rf.OBS_C, lw=P.OBS_LW, zorder=3.5)
        pv.map_furniture(ax, data, x0, x1, y0, y1, 500, "plain", MAP_TEXT_DARK, segments=1)
        ax.set_title("1.5 km around the groins, true scale", fontsize=8)

        axp.clear()
        no, gr = prof["none"][k], prof["groin"][k]
        axp.fill_between(gx, no, gr, where=gr >= no, color=HELD_C, alpha=0.6, lw=0, interpolate=True)
        axp.fill_between(gx, no, gr, where=gr < no, color=CUT_C, alpha=0.6, lw=0, interpolate=True)
        axp.plot(gx, no, color=C["BASE"], lw=1.8, marker="o", ms=3)
        axp.plot(gx, gr, color=INK, lw=1.8, marker="o", ms=3)
        o = obs[y].reindex(gx)
        axp.plot(gx, o.to_numpy(), color=rf.OBS_C, lw=1.2, marker="s", ms=3.5, mfc="white")
        axp.axvline(5.5, color=C["GROIN"], lw=2.0)
        axp.text(5.6, hi - 4, "groins\n(GIS 5|6)", fontsize=6.3, color=C["GROIN"], va="top")
        axp.axhline(0, color=INK_MUTED, lw=0.5)
        e6, e5 = gr[gx == 6][0] - no[gx == 6][0], gr[gx == 5][0] - no[gx == 5][0]
        axp.text(0.98, 0.04, f"groin on minus off:\nGIS 6 (north) {e6:+.0f} m\nGIS 5 (south) {e5:+.0f} m",
                 transform=axp.transAxes, ha="right", va="bottom", fontsize=6.8, color=INK,
                 bbox=dict(facecolor="white", alpha=0.85, edgecolor="none", boxstyle="round,pad=0.3"))
        axp.set_xlim(gx[0] - 0.4, gx[-1] + 0.4)
        axp.set_ylim(lo, hi)
        axp.set_xticks(gx)
        axp.set_xlabel("GIS domain (south to north; drift runs north to south)")
        axp.set_ylabel(f"Shoreline change since {period} (m, seaward +)", fontsize=7.5)
        axp.grid(axis="y")
        open_frame(axp)
        axp.set_title("Along the coast", fontsize=8)
        state = "" if period == rf.CAL else ("; groin on" if y < 2017 else "; groin off (2017 fill)")
        fig.suptitle(f"{y}: the {name} run, groin on vs off ({label}){state}", fontsize=9)

    out = rf.FIGS / f"corrected_groin_on_off_{name}.gif"
    writer = PillowWriter(fps=FPS)
    with writer.saving(fig, str(out), dpi=120):
        for k, y in enumerate(years):
            draw(k, y)
            writer.grab_frame()
        for _ in range(HOLD):
            writer.grab_frame()
    plt.close(fig)
    print("wrote", out.relative_to(cgm.REPO))


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
    for period in (rf.CAL, rf.TEST):
        make(period, data, (x0, x1, y0, y1, dirs))


if __name__ == "__main__":
    main()
