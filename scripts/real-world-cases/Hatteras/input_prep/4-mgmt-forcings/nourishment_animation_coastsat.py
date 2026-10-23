"""
The CoastSat shoreline through each fill, month by month, as an animation.

    python scripts/input_prep/4-mgmt-forcings/nourishment_animation_coastsat.py

One GIF per project the hindcast fires (MP4 too when ffmpeg is on the path).
Each frame is a rolling median shoreline against the median of the 12 months
before placement: (top) both lines straightened along a fitted baseline, gain
and loss shaded; (bottom) the change per transect. For true positions see
nourishment_animation_map_coastsat.py. Dates, footprints and the
pre-fill window are those of nourishment_extent_coastsat.py.

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-10-04
"""
from __future__ import annotations

import shutil
import sys
from pathlib import Path

import numpy as np
import pandas as pd

_HERE = Path(__file__).resolve()
sys.path.insert(0, str(_HERE.parent))
import nourishment_extent_coastsat as X  # noqa: E402

import matplotlib  # noqa: E402
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.animation import FFMpegWriter, FuncAnimation, PillowWriter  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402
from matplotlib.patches import Patch  # noqa: E402

import coastsat_lrr as cl  # noqa: E402
from site_layer.hat_figure_style import (C_1984, C_1984_FILL, C_1997, C_1997_FILL, DOMAIN_AXIS_LABEL, INK,  # noqa: E402
                                         INK_MUTED, apply_style, open_frame, record_caption, save, structures)
from site_layer.hatteras_site_config import HATTERAS_NOURISHMENT_PROJECTS  # noqa: E402

# --- CONFIG ------------------------------------------------------------------
OUT_DIR = X.OUT_DIR / "animations"
MONTHS_BEFORE = 12          # first frame, months before placement starts
MONTHS_AFTER = 24           # last frame, months after placement ends
HALF_WINDOW_DAYS = 45       # each frame is the median of passes within this many days (a 3-month window)
MIN_PASSES = 2              # passes a transect needs in a frame window
# Straightened on purpose: with the cross-shore axis stretched ~16x, the coast's real bends read as cliffs
# (tried true shape 2026-10-04 and reverted); the map animations carry true positions
BASELINE_DEG = 2            # polynomial order of the baseline fitted through the pre-fill line
FPS = 4
DPI = 110
STILL_MONTHS_AFTER = 6      # the still saved beside each GIF, months after placement ends
# -----------------------------------------------------------------------------


# Every pass at every transect in the panel window, one long table
def load_passes(t):
    rows = []
    for tid in t.index:
        df = cl.load_timeseries(str(X.timeseries_path(tid)))
        rows.append(df.assign(transect_id=tid)[["transect_id", "date", "chainage_m"]])
    return pd.concat(rows, ignore_index=True)


# Rolling-median position per transect at each frame date (rows: frame, columns: transect)
def frame_positions(passes, order, frames):
    half = pd.Timedelta(days=HALF_WINDOW_DAYS)
    out = np.full((len(frames), len(order)), np.nan)
    col = {tid: j for j, tid in enumerate(order)}
    for i, f in enumerate(frames):
        w = passes[(passes.date >= f - half) & (passes.date < f + half)]
        g = w.groupby("transect_id").chainage_m.agg(["median", "size"])
        g = g[g["size"] >= MIN_PASSES]
        for tid, v in g["median"].items():
            out[i, col[tid]] = v
    return out


def status(f, start, end):
    if f < start:
        return "before placement", C_1984
    if f <= end:
        return "placement underway", INK
    return "after placement", C_1997


def build(p, d, frame_tbl, filled):
    t, (b0, b1, a0, a1) = X.project_change(p, d, frame_tbl)
    t = t[t.before_m.notna()].sort_values("x")
    start, end = b1, a0 - pd.Timedelta(days=1)
    frames = pd.date_range((start - pd.DateOffset(months=MONTHS_BEFORE)).normalize(),
                           (end + pd.DateOffset(months=MONTHS_AFTER)).normalize(), freq="MS", tz="UTC")
    passes = load_passes(t)
    pos = frame_positions(passes, list(t.index), frames)
    x = t.x.to_numpy()
    ref = t.before_m.to_numpy()
    base = np.polyval(np.polyfit(x, ref, BASELINE_DEG), x)
    ref_off = ref - base
    off = pos - base[None, :]
    change = pos - ref[None, :]
    ctrl = ~t.domain_number.isin(filled).to_numpy()
    return dict(t=t, x=x, ref_off=ref_off, off=off, change=change, frames=frames, start=start, end=end,
                ctrl=ctrl, n_passes=len(passes))


def animate(p, d, A):
    first, last = min(p.gis_domains), max(p.gis_domains)
    x, ref_off, off, change, frames = A["x"], A["ref_off"], A["off"], A["change"], A["frames"]
    lo_a = np.nanpercentile(np.r_[ref_off, off.ravel()], 0.5)
    hi_a = np.nanpercentile(np.r_[ref_off, off.ravel()], 99.5)
    pad = 0.08 * (hi_a - lo_a)
    lim_c = np.nanpercentile(np.abs(change), 99.0) * 1.1

    fig = plt.figure(figsize=(7.2, 5.4))
    ax_a = fig.add_axes([0.10, 0.55, 0.87, 0.33])
    ax_b = fig.add_axes([0.10, 0.17, 0.87, 0.30], sharex=ax_a)
    for ax in (ax_a, ax_b):
        ax.axvspan(first - 0.5, last + 0.5, color="0.92", lw=0, zorder=0)
        ax.grid(axis="y")
        open_frame(ax)
    ax_a.plot(x, ref_off, color=C_1984, lw=1.1, zorder=3)
    (now_line,) = ax_a.plot([], [], color=INK, lw=1.1, zorder=4)
    ax_a.set_ylim(lo_a - pad, hi_a + pad)
    ax_a.set_ylabel("Position from a\nfitted baseline (m)")
    ax_a.text(0.005, 0.97, "ocean ↑", transform=ax_a.transAxes, fontsize=7.5, color=INK_MUTED, va="top")
    ax_a.text(0.005, 0.03, "land ↓", transform=ax_a.transAxes, fontsize=7.5, color=INK_MUTED, va="bottom")
    plt.setp(ax_a.get_xticklabels(), visible=False)
    structures(ax_a, label=True)
    ax_b.axhline(0, color=INK, lw=0.6, zorder=3)
    (ch_line,) = ax_b.plot([], [], color=INK, lw=1.0, zorder=4)
    (bg_line,) = ax_b.plot([], [], color=INK_MUTED, lw=0.9, ls=(0, (3, 2)), zorder=3)
    ax_b.set_ylim(-lim_c, lim_c)
    ax_b.set_xlim(x.min() - 0.3, x.max() + 0.3)
    ax_b.set_ylabel("Change from the\npre-fill line (m)")
    ax_b.set_xlabel(DOMAIN_AXIS_LABEL)
    fig.text(0.10, 0.965, f"{p.name} {p.year}: model GIS {first}–{last}, placed {d.start_date} to {d.end_date}",
             fontsize=10, color=INK, va="top")
    date_txt = fig.text(0.97, 0.915, "", fontsize=11, fontweight="bold", color=INK, ha="right", va="top")
    stat_txt = fig.text(0.10, 0.915, "", fontsize=9.5, fontweight="bold", va="top")
    handles = [Line2D([], [], color=C_1984, lw=1.3, label="Pre-fill: median of the 12 months before"),
               Line2D([], [], color=INK, lw=1.3, label="This month: 3-month rolling median"),
               Patch(color=C_1997_FILL, label="Seaward of pre-fill"), Patch(color=C_1984_FILL, label="Landward"),
               Patch(color="0.92", label="Model footprint"),
               Line2D([], [], color=INK_MUTED, lw=0.9, ls=(0, (3, 2)), label="Background change")]
    fig.legend(handles=handles, loc="lower center", ncol=3, frameon=False, bbox_to_anchor=(0.5, 0.0), fontsize=7.5)
    fills = []

    def draw(i):
        for f in fills:
            f.remove()
        fills.clear()
        o = off[i]
        now_line.set_data(x, o)
        fills.append(ax_a.fill_between(x, ref_off, o, where=o >= ref_off, color=C_1997_FILL, lw=0,
                                       interpolate=True, zorder=1))
        fills.append(ax_a.fill_between(x, ref_off, o, where=o < ref_off, color=C_1984_FILL, lw=0,
                                       interpolate=True, zorder=1))
        ch_line.set_data(x, change[i])
        bg = np.nanmedian(change[i][A["ctrl"]]) if np.isfinite(change[i][A["ctrl"]]).any() else np.nan
        bg_line.set_data([x.min(), x.max()], [bg, bg])
        f = frames[i]
        date_txt.set_text(f.strftime("%b %Y"))
        s, col = status(f, A["start"], A["end"])
        stat_txt.set_text(s)
        stat_txt.set_color(col)
        return [now_line, ch_line, bg_line, date_txt, stat_txt, *fills]

    return fig, draw


def main():
    apply_style()
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    frame_tbl = X.transect_frame()
    dates = pd.read_csv(X.DATES_CSV).set_index("project")
    filled = {g for q in HATTERAS_NOURISHMENT_PROJECTS for g in q.gis_domains}
    has_ffmpeg = shutil.which("ffmpeg") is not None
    for p in sorted(HATTERAS_NOURISHMENT_PROJECTS, key=lambda q: (q.year, min(q.gis_domains))):
        d = dates.loc[p.name]
        A = build(p, d, frame_tbl, filled)
        stem = f"nourishment_animation_{p.year}_{p.name.split()[0].lower()}"
        fig, draw = animate(p, d, A)
        anim = FuncAnimation(fig, draw, frames=len(A["frames"]), blit=False)
        gif = OUT_DIR / f"{stem}.gif"
        anim.save(gif, writer=PillowWriter(fps=FPS), dpi=DPI)
        if has_ffmpeg:
            anim.save(OUT_DIR / f"{stem}.mp4", writer=FFMpegWriter(fps=FPS), dpi=DPI)
        # A still a few months after placement, for the caption and for review
        k = int(np.searchsorted(A["frames"], A["end"] + pd.DateOffset(months=STILL_MONTHS_AFTER)))
        draw(min(k, len(A["frames"]) - 1))
        still = save(fig, OUT_DIR / f"{stem}_still", close=True)[0]
        empty = np.isnan(A["off"]).all(axis=1).sum()
        record_caption(still, (
            f"Still from {gif.name} ({len(A['frames'])} monthly frames, {A['frames'][0]:%b %Y} to "
            f"{A['frames'][-1]:%b %Y}, {FPS} frames per second). {p.name}, {p.year}, in CoastSat. Each frame is the "
            f"median shoreline over ±{HALF_WINDOW_DAYS} days of its date (at least {MIN_PASSES} passes per transect; "
            "gaps where a transect has fewer). Top: that line (black) and the median of the 12 months before "
            "placement (red), each measured from a degree-2 baseline fitted through the red line, so the coast is "
            "straightened and the cross-shore axis stretched (about 16x). Straightening removes the coast's own "
            "shape on purpose, such as the ~150 m step at the Buxton groin, so only the change is enlarged; the "
            "map animations (nourishment_animation_map_*.gif) show the lines at their true positions. Blue "
            "shading is seaward of the pre-fill line, red "
            "landward. Bottom: the change from the pre-fill line per transect; dashed is the median change of the "
            "transects outside every fill footprint. The tag at top left says whether the frame falls before, "
            "during or after placement (dates from nourishment_placement_dates.csv). Grey columns are the model "
            "footprint. Monthly medians carry seasonal swings, so compare a frame with the one 12 months later "
            "rather than with its neighbour."))
        print(f"{gif.name}: {len(A['frames'])} frames, {gif.stat().st_size / 1e6:.1f} MB, "
              f"{len(A['t'])} transects, {A['n_passes']} passes, {empty} empty frames"
              + (", mp4 written" if has_ffmpeg else ""))


if __name__ == "__main__":
    main()
