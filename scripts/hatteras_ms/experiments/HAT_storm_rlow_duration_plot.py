"""
How the equivalent duration is set, drawn for Isabel 2003 and the March 2018 nor'easter.

    python scripts/hatteras_ms/experiments/HAT_storm_rlow_duration_plot.py

The storm's hourly water level against the block Barrier3D runs (24 h trim and
equivalent duration), and the flow over the berm whose area the equivalent block keeps. Details: scripts/hatteras_ms/experiments/README.md.

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-10-01
"""
from __future__ import annotations

import contextlib
import io
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

_HERE = Path(__file__).resolve()
PROJECT_ROOT = next(_p for _p in _HERE.parents if (_p / "pyproject.toml").exists())
sys.path.insert(0, str(PROJECT_ROOT / "scripts"))
sys.path.insert(0, str(_HERE.parent))
import HAT_storm_max_duration as MD  # noqa: E402
import HAT_storm_event_splitting as ES  # noqa: E402
import HAT_storm_rlow_duration as R  # noqa: E402
from site_layer.hat_figure_style import (C, INK, INK_MUTED, apply_style, caption,  # noqa: E402
                                         figsize, open_frame, save, title)

# --- CONFIG ------------------------------------------------------------------
STORMS = [("Hurricane Isabel, 2003", "2003-09-18 12:00", (1996, 2010)),
          ("March 2018 nor'easters", "2018-03-04 12:00", (2010, 2024))]
CONTEXT_H = 6            # hours drawn either side of the event
OUT = R.EXP_DIR / "figures" / "equivalent_duration_isabel_mar2018.png"
# -----------------------------------------------------------------------------


# The split12 event holding `when`, its hourly record, and the record around it
def event(w, when):
    fn = MD.builder_functions()
    with contextlib.redirect_stdout(io.StringIO()):
        df = MD.merged_record(fn, w)
    systems = pd.read_csv(R.EXP_DIR / "storms" / R.wtag(w) / f"{R.wtag(w)}_systems_untrimmed_summary.csv",
                          parse_dates=["StartTime", "EndTime"])
    above = df[df["TWL"] > MD.BERM]
    t = pd.Timestamp(when)
    for _, ev in systems.iterrows():
        hrs = above[(above.index >= ev.StartTime) & (above.index <= ev.EndTime)]
        for g in ES._pieces(hrs.index, R.SPLIT_GAP_H, MD.MIN_DUR):
            piece = hrs.iloc[g]
            if piece.index[0] - pd.Timedelta(days=2) <= t <= piece.index[-1] + pd.Timedelta(days=2) \
                    and abs(piece["TWL"].idxmax() - t) < pd.Timedelta(days=2):
                lo, hi = piece.index[0] - pd.Timedelta(hours=CONTEXT_H), piece.index[-1] + pd.Timedelta(hours=CONTEXT_H)
                return piece, df.loc[lo:hi]
    raise SystemExit(f"no event near {when}")


# One storm's two panels: water level with both blocks, flow over the berm with the equivalent block
def draw(ax_w, ax_q, piece, around, i, name):
    mhw = MD.MHW
    berm = MD.BERM - mhw
    peak_t = piece["TWL"].idxmax()
    hours = lambda idx: (idx - peak_t) / pd.Timedelta(hours=1)  # noqa: E731
    peak = piece["TWL"].max() - mhw
    full = len(piece)
    k = int(np.argmax(piece["TWL"].values))
    lo = min(max(0, k - R.TRIM_H // 2), full - R.TRIM_H) if full > R.TRIM_H else 0
    trim = piece.iloc[lo:lo + R.TRIM_H] if full > R.TRIM_H else piece
    t0, t1 = hours(trim.index[0]), hours(trim.index[-1]) + 1
    d_eq = R.equivalent_duration(piece["TWL"].values)

    # The equivalent block filled, the 24 h block as an outline on top so neither hides the other
    def blocks(ax, y0, y1):
        ax.fill_between([-d_eq / 2, d_eq / 2], y0, y1, color=C["ACCENT_FILL"], lw=0, zorder=1)
        ax.add_patch(plt.Rectangle((t0, y0), t1 - t0, y1 - y0, fill=False, edgecolor=C["BASE"],
                                   lw=1.4, ls=(0, (3, 1.5)), zorder=5))

    # (top) the water level and the blocks Barrier3D runs
    blocks(ax_w, berm, peak)
    ax_w.axhline(berm, color=INK_MUTED, lw=0.8, ls="--", zorder=3)
    ax_w.plot(hours(around.index), around["TWL"] - mhw, color=INK, lw=1.2, zorder=4)
    ax_w.text(0.01, berm - 0.08, "berm", transform=ax_w.get_yaxis_transform(), ha="left", va="top",
              color=INK_MUTED, fontsize=8)
    ax_w.set_ylabel("Total water level (m MHW)")
    title(ax_w, i, name)

    # (bottom) flow over the berm, relative to the peak: the equivalent block has the same area
    rel = (np.clip(piece["TWL"].values - MD.BERM, 0, None) / (piece["TWL"].max() - MD.BERM)) ** R.FLOW_EXPONENT
    x = hours(piece.index)
    blocks(ax_q, 0, 1)
    ax_q.fill_between(x + 0.5, 0, rel, color=C["WATER"], alpha=0.85, lw=0, zorder=3)
    ax_q.plot(x + 0.5, rel, color=INK, lw=0.9, zorder=4)
    ax_q.set_ylim(0, 1.08)
    ax_q.set_ylabel("Flow over the berm\n(fraction of peak)")
    ax_q.set_xlabel("Hours from the peak")
    for ax in (ax_w, ax_q):
        open_frame(ax)
        ax.set_xlim(hours(around.index[0]), hours(around.index[-1]))
    return dict(name=name, peak=peak, hours_above=full, d_eq=d_eq, start=piece.index[0], end=piece.index[-1])


# Both storms side by side, and the caption
def main():
    apply_style()
    fig, axes = plt.subplots(2, 2, figsize=figsize("double", aspect=0.62), constrained_layout=True,
                             gridspec_kw=dict(height_ratios=[1.3, 1]))
    facts = []
    for j, (name, when, w) in enumerate(STORMS):
        piece, around = event(w, when)
        facts.append(draw(axes[0, j], axes[1, j], piece, around, j, name))
    handles = [Line2D([], [], color=INK, lw=1.2, label="Storm water level (Duck gauge + Stockdon runup)"),
               Patch(facecolor="none", edgecolor=C["BASE"], lw=1.4, ls=(0, (3, 1.5)),
                     label="Adopted: the peak held for the 24 h around it"),
               Patch(color=C["ACCENT_FILL"], label="Equivalent duration: the peak held for the same flow"),
               Patch(color=C["WATER"], label="Hourly flow over the berm (same area as the purple block)"),
               Line2D([], [], color=INK_MUTED, lw=0.8, ls="--", label="Berm, 1.34 m MHW")]
    fig.legend(handles=handles, loc="outside lower center", ncol=2, frameon=False)
    a, b = facts
    caption(fig, (
        "How the equivalent storm duration is set. Barrier3D runs every storm as a block: the peak water "
        "level held for `duration` hours. Top row: the hourly total water level (Duck gauge water level + "
        "Stockdon R2, beach slope 0.06) through each storm, with the block the adopted series runs (grey dashed outline, "
        "the 24 h around the peak) and the block the equivalent duration runs (purple, centred on the peak "
        "for display). Bottom row: the flow over the berm each hour, which in Barrier3D goes as the height "
        "above the berm to the power 1.5, as a fraction of the peak hour's flow (blue). The purple block has the "
        "same area as the blue, so it passes the same water. The March 2018 event is one split12 event "
        "because the water never stayed below the berm for 12 h between the early-March nor'easters. "
        f"(a) {a['name']}: {a['hours_above']} h above the berm ({a['start']:%d %b %H:%M}–{a['end']:%d %b %H:%M}), "
        f"peak {a['peak']:.2f} m MHW, equivalent {a['d_eq']} h. "
        f"(b) {b['name']}: {b['hours_above']} h above the berm ({b['start']:%d %b %H:%M}–{b['end']:%d %b %H:%M}), "
        f"peak {b['peak']:.2f} m MHW, equivalent {b['d_eq']} h. "
        "Experiment storms-and-overwash/2026-10-01-rlow-and-duration; not adopted."))
    save(fig, OUT, close=True)
    for f in facts:
        print(f"{f['name']}: {f['start']} .. {f['end']}  {f['hours_above']} h above berm, "
              f"peak {f['peak']:.2f} m MHW, equivalent {f['d_eq']} h")
    print(OUT)


if __name__ == "__main__":
    main()
