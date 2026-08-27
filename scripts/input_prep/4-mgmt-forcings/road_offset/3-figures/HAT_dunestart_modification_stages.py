r"""
HAT_dunestart_modification_stages.py
===============================================================================
The dune-start setback carries TWO modifications on top of the measurement, one
at each edge of the island. This script draws the island BEFORE each of them, so
the progression can be read stage by stage instead of taken on trust.

    setback_dunestart_m           raw measurement            1984: 1 negative, 1 drown
       |  OCEAN-SIDE MOVE -- negative setbacks floored to interior row 0
    setback_dunestart_floored_m   ocean-side applied         1984: 0 negative, 0 drown
       |  BAY-SIDE MOVE -- roadways drowning on a wet bayside row moved seaward
    setback_model_m          both applied, MODEL-FACING 1984: 0 negative, 0 drown

  Counts are from the run: 1984 on 1984-start/v1, 2004 on 2004-start/v1, and
  2004 is 0 negative / 0 drown at every stage. They have moved twice. The
  gap-filled DEM took the bayside drowns that the BAY-SIDE move existed for to
  zero, leaving every stage-0 drown a negative being tested on its wrapped row;
  then the 1984 vintage was re-measured on its OWN product (2026-08-26) and the
  1984 negatives went 6 -> 1, GIS 85 alone. The bay-side move is therefore a
  no-op under both topographies rather than a step with work to do -- a property
  of the DEMs, so the figure counts it from the data rather than asserting it.

  stage 0  HAT_dunestart_stage0_raw.png          before BOTH moves
  stage 1  HAT_dunestart_stage1_ocean_floor.png  ocean-side only, before the
                                                 bay-side move
  stage 2  ../HAT_dunestart_road_on_domains.png  the existing figure, both
                                                 applied -- not redrawn here

HOW FAR each ocean-side move actually is, cropped to the two stretches where it
fires, is HAT_oceanfloor_offset_check.py -> HAT_dunestart_oceanfloor_check.png.
At 90 domains these stage figures cannot show a 10 m move; that one can.

WHY THIS IMPORTS RATHER THAN COPIES
-----------------------------------
The drown test, the interiors, the canvas and the palette all come from
HAT_road_placement_on_domains.py by import. The whole value of a stage figure is
that stage 2 is the SAME test as stages 0 and 1; a transcribed copy that drifted
by one row would make the progression a fiction. Nothing about the test is
re-implemented here -- only the negative-setback case, which stage 2 cannot
contain by construction.

HOW A NEGATIVE SETBACK IS DRAWN  (stage 0 only)
-----------------------------------------------
A negative setback has two positions, and only ONE of them is drawn on the plan
view:

  TRUE      where the road was measured -- seaward of interior row 0, out in the
            dune/beach that the interior array does not cover. Drawn on a sand
            band below the island, on a y-axis extended past 0. This is the only
            position the plan view shows.
  WRAPPED   where CASCADE would actually put it. `int(-70/10) = -7` and
            `xyz_interior_grid[-7:-5, :]` is valid Python indexing from the
            LANDWARD end, so the road is bulldozed into the bay with no error
            raised. NOT drawn on the plan view -- a second mark for one road,
            in a place no measurement supports, reads as two roads rather than
            as one road and its consequence.

The wrap is still what decides those domains' drown state, so it is not lost:
panel C reports their percentages from the WRAPPED rows -- that is what the
model would test -- and marks them apart so the number is never read as a
measurement of the true position. The header text says so on the figure itself.

REQUIREMENTS
------------
  numpy, matplotlib
===============================================================================
"""

from __future__ import annotations

import csv
import importlib.util
import sys
from pathlib import Path

import numpy as np

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from matplotlib.lines import Line2D
from matplotlib.patheffects import withStroke
from matplotlib.transforms import blended_transform_factory

# =============================================================================
# SHARED CODE -- imported, never transcribed
# =============================================================================

PROJECT_ROOT = Path(__file__).resolve().parents[5]
PLACEMENT = (PROJECT_ROOT / "scripts" / "input_prep" / "4-mgmt-forcings"
             / "road_offset" / "1-produce" / "HAT_road_placement_on_domains.py")


def load_placement():
    spec = importlib.util.spec_from_file_location("hat_placement", PLACEMENT)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


P = load_placement()

INIT_ROOT = PROJECT_ROOT / "data" / "hatteras_init"
DUNESTART = INIT_ROOT / "4-mgmt-forcing" / "road_offset" / "dunestart_offset"
OUT_DIR = DUNESTART / "modifications"
DOMAINS_CSV_FMT = "{year}/RoadOffset_{year}_domains.csv"

YEARS = (1984, 2004)
CELL = P.CELL_SIZE_M

# Sand, for the strip seaward of interior row 0 that the interior array does not
# cover. Deliberately not the water grey -- a road measured out there is on the
# beach, not in the sound, and the two must not read the same.
C_SAND = "#e8dcc0"
C_WRAP = "#0b0b0b"

STAGES = [
    dict(key="stage0", column="setback_dunestart_m",
         png="HAT_dunestart_stage0_raw.png",
         title="Stage 0 — the raw measurement, before both moves",
         blurb="Nothing applied. Negative setbacks are drawn where they were "
               "MEASURED — on the sand band, seaward of interior row 0. "
               "CASCADE would not put them there: int(-70/10) = -7 and "
               "xyz_interior_grid[-7:-5] indexes from the LANDWARD end, so the "
               "road is bulldozed into the bay with no error raised. Their "
               "drown percentages below are computed from those wrapped rows.",
         nextfix="The OCEAN-SIDE move (stage 1) floors every negative to "
                 "interior row 0."),
    dict(key="stage1", column="setback_dunestart_floored_m",
         png="HAT_dunestart_stage1_ocean_floor.png",
         title="Stage 1 — ocean-side move applied, before the bay-side move",
         blurb="Negative setbacks have been floored to interior row 0. The "
               "roadways that width-drown on a wet BAYSIDE row are still where "
               "they were measured, so CASCADE would stop managing them at "
               "t=0.",
         nextfix="The BAY-SIDE move (stage 2) relocates each drowning roadway "
                 "to the nearest viable row seaward. See "
                 "../HAT_dunestart_road_on_domains.png."),
]


# =============================================================================
# DATA
# =============================================================================

def read_stage(year: int, column: str) -> dict:
    """One stage's setback per domain, from the per-domain diagnostics CSV.

    Filtered on `setback_model_m` being finite, which is exactly the in-span
    set the model-facing file carries -- so all three stages cover the same 82
    domains and the progression compares like with like.
    """
    path = DUNESTART / DOMAINS_CSV_FMT.format(year=year)
    if not path.is_file():
        return {}
    out = {}
    with open(path, newline="") as f:
        for r in csv.DictReader(f):
            try:
                model = float(r["setback_model_m"])
                value = float(r[column])
            except (KeyError, TypeError, ValueError):
                continue
            if np.isfinite(model) and np.isfinite(value):
                out[int(r["domain"])] = value
    return out


def place_stage(interiors: dict, setbacks: dict) -> dict:
    """
    Where each road sits at this stage, and what bulldoze's test says about it.

    Non-negative setbacks are handed to the imported `place_road` unchanged, so
    stages 0/1 and the stage-2 figure cannot drift apart. Only the negative case
    is handled here, because stage 2 has none by construction.
    """
    positive = {d: v for d, v in setbacks.items() if v >= 0}
    out = P.place_road(interiors, positive)
    for d in out:
        out[d].update(negative=False, true_m=out[d]["start_m"], wrap_m=None)

    for d, sb in sorted(setbacks.items()):
        if sb >= 0:
            continue
        a = interiors.get(d)
        if a is None:
            continue
        n = a.shape[0]
        start = int(sb / CELL)          # bulldoze's truncation, toward zero
        end = start + P.ROAD_WIDTH_CELLS
        if abs(start) >= n:             # past the array even wrapped
            continue

        # Exactly what numpy does with these indices -- negative subscripts
        # count from the landward end. Not a simulation of it.
        sea = P._wet_fraction(a[start - 1, :])
        bay = P._wet_fraction(a[end + 1, :])
        out[d] = dict(
            setback_m=sb,
            start_m=(n + start) * CELL,        # where it really lands
            true_m=sb,                         # where it was measured
            wrap_m=(n + start) * CELL,
            island_m=n * CELL,
            headroom_m=n * CELL - sb,
            seaside=sea, bayside=bay,
            drowned=bool(sea > P.DROWN_PCT or bay > P.DROWN_PCT),
            crashes=False, governing=max(sea, bay),
            negative=True,
        )
    return out


# =============================================================================
# PANELS
# =============================================================================

def draw_island(ax, fig, shown, crop_rows, year, placed, label_sections,
                floor_m):
    colour = P.C_YEAR[year]
    ax.set_facecolor(P.WATER)
    im = ax.imshow(np.ma.masked_invalid(shown), aspect="auto", origin="lower",
                   extent=[0.5, len(P.DOMAINS) + 0.5, -CELL / 2,
                           crop_rows * CELL - CELL / 2],
                   cmap=P.LAND_CMAP, norm=Normalize(P.LAND_VMIN, P.LAND_VMAX),
                   interpolation="nearest")
    cax = ax.inset_axes([1.012, 0.0, 0.014, 1.0])
    cb = fig.colorbar(im, cax=cax, extend="max")
    cb.set_label("interior elev (m MHW)", fontsize=8.5)
    cb.ax.tick_params(labelsize=8)
    cb.outline.set_edgecolor(P.INK_MUTED)

    # The strip seaward of interior row 0. Only drawn when a road is out there,
    # so stage 1 keeps the same axes as stage 2 and the two stack cleanly.
    # Tested against -CELL/2, not 0: that is where the interior image starts, so
    # a stage with no negatives has floor_m == -CELL/2 and must draw no band.
    if floor_m < -CELL / 2:
        # No in-band caption: any text here sits on top of the drowned roads
        # this band exists to show. The figure legend names the colour instead.
        ax.axhspan(floor_m, -CELL / 2, color=C_SAND, lw=0, zorder=1)
        ax.axhline(-CELL / 2, color=P.INK_SECOND, lw=1.0, ls=(0, (4, 2)),
                   zorder=5)

    halo = [withStroke(linewidth=4.4, foreground=P.SURFACE)]
    seg = {False: ([], []), True: ([], [])}
    for d, p in placed.items():
        xs, ys = seg[bool(p["drowned"])]
        xs += [d - 0.5, d + 0.5, np.nan]
        ys += [p["true_m"], p["true_m"], np.nan]

    ax.plot(seg[False][0], seg[False][1], color=colour, lw=2.6,
            solid_capstyle="butt", path_effects=halo, zorder=6)
    if seg[True][0]:
        ax.plot(seg[True][0], seg[True][1], color=P.C_DROWN, lw=3.4,
                solid_capstyle="butt", path_effects=halo, zorder=7)

    # The wrapped position is NOT drawn on the plan view. It was, and it put a
    # second mark for one road in a place no measurement supports, which read as
    # two roads rather than one road and its consequence. The wrap still governs
    # the drown state of these domains -- that is where it belongs, and panel C
    # marks them.
    n_neg = sum(1 for p in placed.values() if p.get("wrap_m") is not None)
    drowned = [d for d, p in placed.items() if p["drowned"]]
    notes = []
    if n_neg:
        notes.append(f"{n_neg} NEGATIVE — bulldozed into the bay")
    if drowned:
        notes.append(f"{len(drowned)} of {len(placed)} drown at initialisation")
    ax.text(0.5, 0.985,
            "\n".join(notes) if notes
            else f"0 of {len(placed)} drown at initialisation",
            transform=ax.transAxes, ha="center", va="top", fontsize=9,
            color=P.C_DROWN if notes else P.INK_SECOND, zorder=11,
            bbox=dict(fc=P.SURFACE, ec=P.C_DROWN if notes else P.INK_MUTED,
                      alpha=0.95, pad=3.0))

    # Section names sit on a second row so the count above them stays clear.
    # x in data, y in axes fraction -- the panels do not share a y range once a
    # sand band is added, and a data-coordinate y would drift between them.
    tr = blended_transform_factory(ax.transData, ax.transAxes)
    for (lo, hi), name in P.SECTIONS:
        ax.axvline(hi + 0.5, color=P.SURFACE, lw=1.0, alpha=0.55, zorder=4)
        if label_sections:
            ax.text((lo + hi) / 2, 0.855, name, transform=tr, ha="center",
                    va="top", fontsize=8, color=P.INK_SECOND, zorder=9,
                    bbox=dict(fc=P.SURFACE, ec="none", alpha=0.8, pad=1.6))

    ax.set_title(f"{year} road", loc="left", fontsize=11.5, weight="semibold",
                 color=colour, pad=4)
    ax.set_ylabel("m landward of\ninterior row 0", fontsize=9)
    ax.set_ylim(floor_m, crop_rows * CELL - CELL / 2)
    plt.setp(ax.get_xticklabels(), visible=False)


def draw_drown_panel(ax, placed_by_year):
    ax.axhspan(P.DROWN_PCT * 100, 104, color=P.C_DROWN, alpha=0.07, lw=0,
               zorder=1)
    ax.axhline(P.DROWN_PCT * 100, color=P.C_DROWN, lw=1.3, ls=(0, (4, 3)),
               zorder=3, label=f"drown threshold, {P.DROWN_PCT * 100:.0f}%")
    for year, pl in sorted(placed_by_year.items()):
        xs = sorted(pl)
        ax.plot(xs, [pl[d]["governing"] * 100 for d in xs],
                color=P.C_YEAR[year], lw=1.8, marker="o", ms=3.0, zorder=5,
                label=f"{year}")
        bad = [d for d in xs if pl[d]["drowned"]]
        if bad:
            ax.plot(bad, [pl[d]["governing"] * 100 for d in bad], lw=0,
                    marker="v", ms=9, mfc=P.C_DROWN, mec=P.SURFACE, mew=1.2,
                    zorder=7)
        # Negatives are plotted from the WRAPPED rows, so they are marked apart
        # -- the number is real, but it does not describe the measured position.
        neg = [d for d in xs if pl[d].get("negative")]
        if neg:
            ax.plot(neg, [pl[d]["governing"] * 100 for d in neg], lw=0,
                    marker="X", ms=8.5, mfc=C_WRAP, mec=P.SURFACE, mew=1.0,
                    zorder=8,
                    label="negative — tested on WRAPPED rows" if year ==
                          min(placed_by_year) else None)
    ax.set_ylabel("% of bordering cells\nat or below 0 m MHW", fontsize=9)
    ax.set_xlabel("Barrier3D / GIS domain   "
                  "(1 = Cape Point / south  →  90 = Rodanthe / north)")
    ax.set_xlim(0.5, len(P.DOMAINS) + 0.5)
    ax.set_ylim(-4, 104)
    ax.grid(axis="y", color=P.INK_MUTED, alpha=0.22, lw=0.7)
    ax.set_axisbelow(True)
    ax.legend(loc="upper left", fontsize=8.5, ncol=4, framealpha=0.92)
    ax.set_title("bulldoze's own drown test — the worse of the two rows "
                 "BORDERING the road (road_start − 1, road_end + 1). Above the "
                 "line, CASCADE stops managing this roadway",
                 loc="left", fontsize=10)


# =============================================================================
# FIGURE
# =============================================================================

def build_figure(stage: dict, per, crop_rows) -> dict:
    print(f"\n{'=' * 84}")
    print(f"{stage['key'].upper()} -- {stage['column']}")
    print("=" * 84)

    placed = {}
    for year in YEARS:
        sb = read_stage(year, stage["column"])
        if not sb:
            print(f"  [skip] {year}: no {stage['column']} in the domains CSV")
            continue
        placed[year] = place_stage(per[year]["interiors"], sb)

    if not placed:
        print("  [skip] nothing to draw")
        return {}

    floor_m = min([p["true_m"] for pl in placed.values() for p in pl.values()]
                  + [0.0])
    floor_m = min(floor_m - 30.0, -CELL / 2) if floor_m < 0 else -CELL / 2

    fig = plt.figure(figsize=(16.5, 10.4))
    gs = fig.add_gridspec(3, 1, height_ratios=[1.25, 1.25, 0.85], hspace=0.19,
                          left=0.065, right=0.905, top=0.858, bottom=0.062)
    axes = [fig.add_subplot(gs[0])]
    axes.append(fig.add_subplot(gs[1], sharex=axes[0]))
    axes.append(fig.add_subplot(gs[2], sharex=axes[0]))

    for ax, year, lab in ((axes[0], YEARS[0], True), (axes[1], YEARS[1], False)):
        if year in placed:
            draw_island(ax, fig, per[year]["shown"], crop_rows, year,
                        placed[year], lab,
                        floor_m)

    draw_drown_panel(axes[2], placed)

    fig.text(0.065, 0.985, f"NC-12, dune-start method — {stage['title']}",
             fontsize=14, va="top", weight="semibold")
    fig.text(0.065, 0.957,
             "Each panel on its own period's topography — "
             + ", ".join(f"{y} on {P.topo_label(y)}" for y in YEARS)
             + ". They are different islands, not one island twice. "
             f"Road drawn where bulldoze puts it: road_start = "
             f"int(setback / 10 m).\n{stage['blurb']}\nNEXT: "
             f"{stage['nextfix']}",
             fontsize=9, color=P.INK_SECOND, va="top", linespacing=1.5)

    handles = [
        Line2D([], [], color=P.C_1984, lw=2.6, label="1984 road"),
        Line2D([], [], color=P.C_2004, lw=2.6, label="2004 road"),
        Line2D([], [], color=P.C_DROWN, lw=3.4,
               label="drowns at initialisation"),
        Line2D([], [], color=P.WATER, lw=8, label="off-island / sentinel water"),
    ]
    if floor_m < -CELL / 2:
        # The wrap marker carries no legend entry -- the subtitle already says
        # what the hollow marks are, and the entry was the longest item in the
        # row. The band keeps its swatch.
        handles.append(Line2D([], [], color=C_SAND, lw=8,
                              label="seaward of interior row 0"))
    fig.legend(handles=handles, loc="upper left", bbox_to_anchor=(0.065, 0.898),
               ncol=6, fontsize=8.5, framealpha=0.0, borderpad=0.4,
               columnspacing=1.6, handlelength=2.6)

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    out_png = OUT_DIR / stage["png"]
    fig.savefig(out_png, dpi=150, bbox_inches="tight", facecolor=P.SURFACE)
    plt.close(fig)
    print(f"  [out] {out_png}")

    summary = {}
    for year, pl in sorted(placed.items()):
        neg = [d for d, p in pl.items() if p.get("negative")]
        drowned = [d for d, p in pl.items() if p["drowned"]]
        summary[year] = dict(n=len(pl), negative=neg, drowned=drowned)
        print(f"  {year}: {len(pl)} domains | {len(neg)} negative | "
              f"{len(drowned)} drown")
        if neg:
            print(f"      negative : {neg}")
        if drowned:
            print(f"      drowning : {drowned}")
    return summary


def main() -> int:
    # PER VINTAGE (2026-08-26). This called P.load_interiors() with no
    # argument -- one island under both panels, which the caption used to state
    # outright. load_interiors() now requires a year, so this could not survive
    # the change silently.
    per, crop_rows, _max_rows = P.load_years(YEARS)
    for year in YEARS:
        print(f"  {year}: {len(per[year]['interiors'])} interiors from "
              f"{P.topo_label(year)}")
    print(f"  drown test imported from {PLACEMENT.name}")

    results = {s["key"]: build_figure(s, per, crop_rows) for s in STAGES}

    print(f"\n{'=' * 84}")
    print("PROGRESSION -- domains failing at each stage")
    print("=" * 84)
    print(f"  {'stage':<34} {'year':>6} {'negative':>9} {'drowning':>9}")
    rows = [("stage 0  raw measurement", "stage0"),
            ("stage 1  + ocean-side floor", "stage1")]
    for label, key in rows:
        for year, s in sorted(results.get(key, {}).items()):
            print(f"  {label:<34} {year:>6} {len(s['negative']):>9} "
                  f"{len(s['drowned']):>9}")
    print(f"  {'stage 2  + bay-side relocation':<34} {'(see':>6} "
          f"{'../HAT_dunestart_road_on_domains.png)':>9}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
