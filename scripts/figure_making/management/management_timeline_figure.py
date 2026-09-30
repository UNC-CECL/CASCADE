#!/usr/bin/env python3
"""
The management of NC-12 and the beach as a timeline, domain against year: the runs and the whole record.

    python scripts/figure_making/management/management_timeline_figure.py

Events come from hatteras_site_config, period ends from HATTERAS_PERIODS;
writes timeline_1996_2024.png and timeline_1984_2024.png to
output/figures/3-model-inputs/4-management/. Details: scripts/figure_making/management/README.md.

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-09-30
"""
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import Rectangle
import matplotlib.ticker as ticker

# House style (site_layer/hat_figure_style.py), applied at import
import sys as _sys
from pathlib import Path as _P
_sys.path.insert(0, str(next(_q for _q in _P(__file__).resolve().parents
                             if (_q / "pyproject.toml").exists()) / "scripts"))
from site_layer.hat_figure_style import (apply_style, C, C_1984, C_1997, C_1984_FILL,
                              C_1997_FILL, INK, INK_MUTED, GRID_C, figsize,
                              record_caption)
apply_style()

from site_layer.hatteras_site_config import (HATTERAS_NOURISHMENT_PROJECTS,
                                  HATTERAS_ROAD_EVENTS,
                                  HATTERAS_PERIODS,
                                  HATTERAS_ANNOTATIONS)

# Paths resolved from the repo root (ORGANIZATION.md rule 5)
_PATH_REPO = next(_p for _p in _P(__file__).resolve().parents
                  if (_p / "pyproject.toml").exists())
from site_layer import hat_figure_style as _hs  # noqa: E402
# --- CONFIG ------------------------------------------------------------------
OUT = _hs.figure_dir("inputs", "4-management")

# The period starts in use; ends come from HATTERAS_PERIODS
PERIOD_STARTS = (1996, 2010)
PERIODS = [(st, HATTERAS_PERIODS[st]["end_year"]) for st in PERIOD_STARTS]
RECORD_Y0 = 1984                       # the first year the record covers
PW = 0.6                               # point-event width, in years

# Colour is for management and period; zones are greys
C_N = C["ADDED"]                       # sediment added
C_R = C["ROAD"]                        # NC-12
C_VILLAGE = "0.86"
C_INTERV = "0.955"

_T = HATTERAS_ANNOTATIONS.town_spans
COMMUNITIES = [
    ('Cape Point',                              1,  6,  C_INTERV),
    ('Buxton',                                  *_T["Buxton"], C_VILLAGE),
    ('Buxton–Avon\n(inter-village)',       9, 20,  C_INTERV),
    ('Avon',                                    *_T["Avon"], C_VILLAGE),
    ('Avon–Tri-Village\n(inter-village)', 32, 67,  C_INTERV),
    ('Tri-Village: Salvo /\nWaves / Rodanthe',  *_T["Tri-Village"], C_VILLAGE),
    ('Pea Island NWR',                         84, 90,  C_INTERV),
]

# The map's encoding: fill orange, relocation NC-12 ink, the bridge a hatch
STYLE = {
    'nourish': dict(facecolor=C_N, edgecolor='none'),
    'road': dict(facecolor=C_R, edgecolor='none', alpha=0.85),
    'bridge': dict(facecolor='none', edgecolor=INK, hatch='////', lw=0.0),
}
# -----------------------------------------------------------------------------


# (year0, year1, gis_lo, gis_hi, kind, label) for the window, from the config
def events_in(y0, y1):
    out = []
    for p in sorted(HATTERAS_NOURISHMENT_PROJECTS, key=lambda q: q.year):
        if p.enabled and y0 <= p.year <= y1:
            g = sorted(p.gis_domains)
            out.append((p.year, p.year, g[0], g[-1], 'nourish',
                        f"{p.name.split()[0]} fill {p.year}\n"
                        f"GIS {g[0]}–{g[-1]}"))
    for e in HATTERAS_ROAD_EVENTS:
        if not (y0 <= e.year <= y1):
            continue
        if hasattr(e, 'displacement_m'):
            g, d = sorted(e.displacement_m), e.displacement_m
            out.append((e.year, e.year, g[0], g[-1], 'road',
                        f"relocated {e.year}\n"
                        f"GIS {g[0]}–{g[-1]}, "
                        f"{min(d.values()):.0f}–{max(d.values()):.0f} m"))
        else:
            g = sorted(e.gis_domains)
            out.append((e.year, y1, g[0], g[-1], 'bridge',
                        f"bridge {e.year}\n"
                        f"road off GIS {g[0]}–{g[-1]}"))
    return out


# Does a segment miss every rectangle? (sampled)
def _clear(x0, y0_, x1, y1_, rects, steps=64):
    for i in range(steps + 1):
        t = i / steps
        x, y = x0 + (x1 - x0) * t, y0_ + (y1_ - y0_) * t
        for rx0, rx1, ry0, ry1 in rects:
            if rx0 <= x <= rx1 and ry0 <= y <= ry1:
                return False
    return True


# Place each event's label clear of bars, other labels and leaders
def place(fig, ax, items, bars, y0, y1, fontsize=8, pad_x=0.45, pad_y=2.0):
    renderer = fig.canvas.get_renderer()
    inv = ax.transData.inverted()
    taken = list(bars)
    placed = []
    for x_anchor, y_anchor, text in items:
        probe = ax.text(0, 0, text, fontsize=fontsize)
        bb = probe.get_window_extent(renderer=renderer)
        probe.remove()
        (px0, py0), (px1, py1) = inv.transform([[bb.x0, bb.y0], [bb.x1, bb.y1]])
        w, h = px1 - px0, py1 - py0
        # Test against every box but the label's own bar
        others = [r for r in taken
                  if not (r[0] <= x_anchor <= r[1] and r[2] <= y_anchor <= r[3])]

        cands = []
        for dy in (0, 11, -11, 17, -17, 22, -22, 28, -28, 33, -33, 40, -40):
            for side in (-1, +1):        # left of the bar first, then right
                cx = x_anchor + side * pad_x
                bx0 = cx - w if side < 0 else cx
                by0 = y_anchor + dy - h / 2
                box = (bx0 - pad_x, bx0 + w + pad_x,
                       by0 - pad_y / 2, by0 + h + pad_y / 2)
                if box[0] < y0 - 0.6 or box[1] > y1 + 0.6:
                    continue
                if box[2] < 0 or box[3] > 92:
                    continue
                cands.append((dy, side, cx, box))

        chosen = None
        for dy, side, cx, box in cands:
            if any(not (box[1] < t[0] or box[0] > t[1]
                        or box[3] < t[2] or box[2] > t[3]) for t in others):
                continue
            # and the leader this position would need, against the same set
            if dy and not _clear(cx - side * 0.12, y_anchor + dy,
                                 x_anchor, y_anchor, others):
                continue
            chosen = (dy, side, cx, box)
            break
        # Never drop a label: if nothing is clear, take the first that fits
        if chosen is None and cands:
            chosen = cands[0]
        if chosen is None:
            continue

        dy, side, cx, box = chosen
        ax.text(cx, y_anchor + dy, text, fontsize=fontsize, color=INK,
                ha='right' if side < 0 else 'left', va='center', zorder=11)
        taken.append(box)
        placed.append((cx, y_anchor + dy, x_anchor, y_anchor, side))
    for cx, cy, x_anchor, y_anchor, side in placed:
        if abs(cy - y_anchor) > 1e-9:
            ax.plot([cx - side * 0.12, x_anchor], [cy, y_anchor],
                    color=INK_MUTED, lw=0.5, zorder=4)


# One timeline from y0; pre_run draws the years before the first period start
def build(y0, stem, pre_run=False):
    y1 = max(en for _, en in PERIODS)
    events = events_in(y0, y1)
    fig, ax = plt.subplots(figsize=figsize("double", height=3.9))
    # Explicit margins keep the figure inside the 190 mm width
    fig.subplots_adjust(left=0.115, right=0.790, bottom=0.235, top=0.97)

    for _, d0, d1, fc in COMMUNITIES:
        ax.axhspan(d0 - 0.5, d1 + 0.5, color=fc, zorder=1, lw=0)
    for _, d0, _, _ in COMMUNITIES[1:]:
        ax.axhline(d0 - 0.5, color=GRID_C, lw=0.5, zorder=2)

    # Period bar in the vintage pair, matching the map's alignments
    BAR_Y0, BAR_Y1 = -9.0, -4.5
    segs = [(st, en, fill, ink, f"Period {i}  ({st}–{en})")
            for i, ((st, en), fill, ink) in enumerate(
                zip(PERIODS, (C_1984_FILL, C_1997_FILL), (C_1984, C_1997)), 1)]
    if pre_run:
        # No years on this one: the axis and the dashed break carry them
        segs.insert(0, (y0, PERIODS[0][0], "0.93", INK_MUTED, "before the runs"))
    for st, en, fill, ink, label in segs:
        ax.add_patch(Rectangle((st, BAR_Y0), en - st, BAR_Y1 - BAR_Y0, fc=fill,
                               ec=ink, lw=0.7, zorder=5, clip_on=False))
        ax.text((st + en) / 2, (BAR_Y0 + BAR_Y1) / 2, label, fontsize=7.5,
                color=ink, ha='center', va='center', zorder=6, clip_on=False)
    for st, _, _, _, _ in segs[1:]:
        ax.axvline(st, color=INK_MUTED, lw=0.8, ls=(0, (4, 3)), zorder=9)

    bars = []
    for (yr0, yr1, d0, d1, kind, _) in events:
        point = yr0 == yr1
        bx = yr0 - PW / 2 if point else yr0
        w = PW if point else yr1 - yr0
        ax.add_patch(Rectangle((bx, d0 - 0.5), w, (d1 + 0.5) - (d0 - 0.5),
                               zorder=5, **STYLE[kind]))
        bars.append((bx, bx + w, d0 - 0.5, d1 + 0.5))

    ax.set_xlim(y0 - 0.6, y1 + 0.6)
    ax.set_ylim(-10.0, 95.0)
    place(fig, ax, [(yr0, (d0 + d1) / 2, lab)
                    for (yr0, _, d0, d1, _, lab) in events], bars, y0, y1)

    for lbl, d0, d1, fc in COMMUNITIES:
        village = fc == C_VILLAGE
        ax.text(y1 + 1.0, (d0 + d1) / 2, lbl, fontsize=7.5,
                color=INK if village else INK_MUTED, va='center', ha='left',
                style='normal' if village else 'italic', clip_on=False, zorder=12)

    ax.set_xlabel('Year', labelpad=4)
    ax.set_ylabel('CASCADE domain\n(1 = S / Cape Point  →  90 = N / Oregon Inlet)',
                  fontsize=9)
    ax.xaxis.set_major_locator(ticker.MultipleLocator(4))
    ax.xaxis.set_minor_locator(ticker.MultipleLocator(1))
    ax.yaxis.set_major_locator(ticker.FixedLocator([0, 20, 40, 60, 80, 90]))
    ax.grid(axis='x', which='major', color=GRID_C, lw=0.5, zorder=1.5)
    for spine in ('top', 'right'):
        ax.spines[spine].set_visible(False)
    ax.set_yticks([t for t in ax.get_yticks() if 0 <= t <= 90])
    ax.tick_params(labelsize=8)

    ax.legend(handles=[mpatches.Patch(label='beach nourishment', **STYLE['nourish']),
                       mpatches.Patch(label='NC-12 relocation', **STYLE['road']),
                       mpatches.Patch(facecolor='none', edgecolor=INK, hatch='////',
                                      label='road removed after the bridge')],
              loc='upper center', bbox_to_anchor=(0.5, -0.17), ncol=3,
              frameon=False, fontsize=8, handlelength=1.6, handleheight=1.0,
              columnspacing=2.0)

    png = OUT / f"{stem}.png"
    fig.savefig(png, dpi=300, facecolor='white')
    fig.savefig(OUT / "supporting" / f"{stem}.pdf", facecolor='white')
    plt.close(fig)
    print(f"Saved: {png}")
    return png, y1


_SHARED = ("by domain and year: the beach-nourishment projects, the NC-12 "
           "relocations and the Rodanthe bridge, against the community zones and "
           "the run periods. Nourishment and relocation are single years, drawn as "
           "a narrow slice; the bridge persists to the end of the record. The same "
           "events as management_footprint.png, against time rather than against "
           "the reach.")

# Run: build the runs and the record timelines and record their captions
def main():
    runs_png, y1 = build(PERIODS[0][0], f"timeline_{PERIODS[0][0]}_"
                                        f"{max(en for _, en in PERIODS)}")
    record_png, _ = build(RECORD_Y0, f"timeline_{RECORD_Y0}_{y1}", pre_run=True)

    # the caption lives beside the figure, not on it (figure_making/STYLE.md)
    record_caption(runs_png,
        f"The management the model applies, {PERIODS[0][0]}–{y1}, " + _SHARED +
        " The axis starts at the first period start, so every mark is an event "
        "some run applies; the 1989 Pea Island relocation precedes both periods "
        f"and is absent here. timeline_{RECORD_Y0}_{y1}.png is the same figure "
        "over the whole record and does carry it.")
    record_caption(record_png,
        f"The management record, {RECORD_Y0}–{y1}, " + _SHARED +
        f" The axis covers the whole record rather than the modelled window, so "
        f"{RECORD_Y0}–{PERIODS[0][0]} is drawn as the stretch before the runs "
        "and the 1989 Pea Island relocation appears. That relocation is not "
        "applied by any run: it precedes both periods, and it reaches the model "
        "as the starting road position rather than as an event. "
        f"timeline_{PERIODS[0][0]}_{y1}.png is the modelled window alone.")


if __name__ == "__main__":
    main()
