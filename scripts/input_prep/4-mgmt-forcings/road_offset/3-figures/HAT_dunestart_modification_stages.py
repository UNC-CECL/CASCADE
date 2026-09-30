"""
The dune-start setback before each of its two modifications, drawn stage by stage.

    python scripts/input_prep/4-mgmt-forcings/road_offset/3-figures/HAT_dunestart_modification_stages.py

One figure per stage from the placement module's drawing code. Details: scripts/input_prep/4-mgmt-forcings/road_offset/README.md.

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-09-18
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
from matplotlib.lines import Line2D
from matplotlib.patheffects import withStroke

# Shared code -- imported, never transcribed

PROJECT_ROOT = next(_p for _p in Path(__file__).resolve().parents
                    if (_p / "pyproject.toml").exists())
PLACEMENT = (PROJECT_ROOT / "scripts" / "input_prep" / "4-mgmt-forcings"
             / "road_offset" / "1-produce" / "HAT_road_placement_on_domains.py")


# HAT_road_placement_on_domains.py, loaded as a module
def load_placement():
    spec = importlib.util.spec_from_file_location("hat_placement", PLACEMENT)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


P = load_placement()

# The house style, through the module that already resolves it
sys.path.insert(0, str(PROJECT_ROOT / "scripts"))
from site_layer.hat_figure_style import (  # noqa: E402
    C, DOMAIN_AXIS_LABEL, caption, figsize, open_frame, save,
    spines_for_image, town_bands, _title)

# --- CONFIG ------------------------------------------------------------------
INIT_ROOT = PROJECT_ROOT / "data" / "hatteras_init"
from site_layer import hat_topo_version as _tv  # noqa: E402
DUNESTART = _tv.ROAD_SETBACK_ROOT
OUT_DIR = DUNESTART / "modifications"
DOMAINS_CSV_FMT = "measured/{year}/RoadOffset_{year}_domains.csv"

YEARS = (1984, 2004)
CELL = P.CELL_SIZE_M

# Sand, for the strip seaward of interior row 0 that the interior array does not cover
C_SAND = C["ADDED_FILL"]
C_WRAP = C["INK"]

STAGES = [
    # `png` names are fixed (referenced elsewhere); the labels say what was done
    dict(key="stage0", column="setback_dunestart_m",
         png="HAT_dunestart_stage0_raw.png",
         title="the setback as measured, before either correction",
         blurb="Nothing applied. Negative setbacks are drawn where they were "
               "measured, on the sand band seaward of interior row 0. CASCADE "
               "would not put them there: int(-70/10) = -7 and "
               "xyz_interior_grid[-7:-5] indexes from the landward end, so the "
               "road is bulldozed into the bay with no error raised. Their "
               "drown percentages in the lower panel are computed from those "
               "landward-indexed rows.",
         nextfix="The ocean-side correction floors every negative setback to "
                 "interior row 0."),
    dict(key="stage1", column="setback_dunestart_floored_m",
         png="HAT_dunestart_stage1_ocean_floor.png",
         title="after the ocean-side correction, before the bay-side one",
         blurb="Negative setbacks have been floored to interior row 0. The "
               "roadways that width-drown on a wet bayside row are still where "
               "they were measured, so CASCADE would stop managing them at "
               "t=0.",
         nextfix="The bay-side correction relocates each drowning roadway to "
                 "the nearest viable row seaward, giving the model-facing "
                 "setback drawn in HAT_dunestart_road_on_domains.png."),
]
# -----------------------------------------------------------------------------


# One stage's setback per domain, from the per-domain diagnostics CSV
def read_stage(year: int, column: str) -> dict:
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


# Where each road sits at this stage, and what bulldoze's test says about it
def place_stage(interiors: dict, setbacks: dict) -> dict:
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

        # Exactly what numpy does with these indices -- negative subscripts count from the landward end
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


# Panels

# One vintage's island at one stage
def draw_island(ax, fig, shown, crop_rows, year, placed, panel_index,
                floor_m):
    colour = P.C_YEAR[year]
    ax.set_facecolor(P.NODATA)
    im = ax.imshow(np.ma.masked_invalid(shown), aspect="auto", origin="lower",
                   extent=[0.5, len(P.DOMAINS) + 0.5, -CELL / 2,
                           crop_rows * CELL - CELL / 2],
                   cmap=P.LAND_CLASS_CMAP, norm=P.LAND_CLASS_NORM,
                   interpolation="nearest")
    cax = ax.inset_axes([1.012, 0.0, 0.014, 1.0])
    cb = fig.colorbar(im, cax=cax, spacing="uniform",
                      ticks=P.LAND_CLASS_BOUNDS[1:-1])
    cb.set_label("elevation (m MHW)")
    cb.outline.set_edgecolor(P.INK_MUTED)
    cb.outline.set_linewidth(0.6)

    # The strip seaward of interior row 0
    if floor_m < -CELL / 2:
        # No in-band caption: the legend names the colour
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

    # The wrapped position is NOT drawn on the plan view

    # The three village spans, named once on the upper panel.
    ax.set_xlim(0.5, len(P.DOMAINS) + 0.5)
    if panel_index == 0:
        town_bands(ax, strip=0.075, shade=P.SURFACE)

    spines_for_image(ax)
    _title(ax, panel_index, f"NC-12 in {year}")
    ax.set_ylabel("m landward of\ninterior row 0")
    ax.set_ylim(floor_m, crop_rows * CELL - CELL / 2)
    plt.setp(ax.get_xticklabels(), visible=False)


# The drown test per domain, against the threshold
def draw_drown_panel(ax, placed_by_year, panel_index=2):
    ax.axhspan(P.DROWN_PCT * 100, 104, color=P.C_DROWN, alpha=0.07, lw=0,
               zorder=1)
    ax.axhline(P.DROWN_PCT * 100, color=P.C_DROWN, lw=1.1, ls=(0, (4, 3)),
               zorder=3, label=f"threshold, {P.DROWN_PCT * 100:.0f}%")
    for year, pl in sorted(placed_by_year.items()):
        xs = sorted(pl)
        ax.plot(xs, [pl[d]["governing"] * 100 for d in xs],
                color=P.C_YEAR[year], lw=1.2, marker="o", ms=2.0, zorder=5,
                label=f"{year}")
        bad = [d for d in xs if pl[d]["drowned"]]
        if bad:
            ax.plot(bad, [pl[d]["governing"] * 100 for d in bad], lw=0,
                    marker="v", ms=6, mfc=P.C_DROWN, mec=P.SURFACE, mew=0.8,
                    zorder=7)
        # Negatives are tested where numpy actually lands them, at the landward end of the array
        neg = [d for d in xs if pl[d].get("negative")]
        if neg:
            ax.plot(neg, [pl[d]["governing"] * 100 for d in neg], lw=0,
                    marker="X", ms=6, mfc=C_WRAP, mec=P.SURFACE, mew=0.8,
                    zorder=8,
                    label="negative setback, tested from the landward end"
                          if year == min(placed_by_year) else None)
    ax.set_ylabel("% of bordering cells\nat or below 0 m MHW")
    ax.set_xlabel(DOMAIN_AXIS_LABEL)
    ax.set_xlim(0.5, len(P.DOMAINS) + 0.5)
    ax.set_ylim(-4, 104)
    ax.grid(axis="y")
    ax.set_axisbelow(True)
    town_bands(ax, label=False)
    open_frame(ax)
    ax.legend(loc="upper left", ncol=2, fontsize=7)
    _title(ax, panel_index, "wet cells bordering the road")


# Figure

# One stage's figure, both vintages
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

    fig = plt.figure(figsize=figsize("double", height=6.2))
    gs = fig.add_gridspec(3, 1, height_ratios=[1.25, 1.25, 0.85], hspace=0.30,
                          left=0.105, right=0.870, top=0.955, bottom=0.115)
    axes = [fig.add_subplot(gs[0])]
    axes.append(fig.add_subplot(gs[1], sharex=axes[0]))
    axes.append(fig.add_subplot(gs[2], sharex=axes[0]))

    for i, (ax, year) in enumerate(((axes[0], YEARS[0]), (axes[1], YEARS[1]))):
        if year in placed:
            draw_island(ax, fig, per[year]["shown"], crop_rows, year,
                        placed[year], i, floor_m)

    draw_drown_panel(axes[2], placed, panel_index=2)

    handles = [
        Line2D([], [], color=P.C_1984, lw=2.6, label="NC-12 in 1984"),
        Line2D([], [], color=P.C_2004, lw=2.6, label="NC-12 in 2004"),
        Line2D([], [], color=P.C_DROWN, lw=3.4,
               label="drowns at initialisation"),
        Line2D([], [], color=P.NODATA, lw=8, label="outside the extraction"),
    ]
    if floor_m < -CELL / 2:
        # The landward-index marker carries no legend entry
        handles.append(Line2D([], [], color=C_SAND, lw=8,
                              label="seaward of interior row 0"))
    fig.legend(handles=handles, loc="lower center",
               bbox_to_anchor=(0.5, -0.004), ncol=5, frameon=False,
               columnspacing=1.4, handlelength=2.2)

    counts = "; ".join(
        f"{year}: {sum(1 for p in pl.values() if p.get('negative'))} negative, "
        f"{sum(1 for p in pl.values() if p['drowned'])} drowning of {len(pl)}"
        for year, pl in sorted(placed.items()))
    caption(fig, (
        f"The dune-start setback for NC-12, {stage['title']}. Domain 1 is at "
        f"Cape Point in the south and domain 90 at Pea Island in the north; "
        f"the shaded spans are the villages (Buxton, Avon, Tri-Village). "
        f"(a, b) the road as roadway_manager.bulldoze places it on each "
        f"period's own extraction, "
        + ", ".join(f"{y} on {P.topo_label(y)}" for y in YEARS if y in placed)
        + f". Interior elevation is in classes relative to mean high water; "
          f"cells outside the extraction carry no data and are drawn grey. "
          f"The sand band below interior row 0 is the beach the interior array "
          f"does not cover. {stage['blurb']} (c) bulldoze's own drown test: "
          f"the wetter of the two rows bordering the bulldozed band "
          f"(road_start − 1, road_end + 1); above "
          f"{P.DROWN_PCT * 100:.0f}% CASCADE stops managing the roadway. "
          f"At this stage {counts}. {stage['nextfix']}"))

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    out_png = OUT_DIR / stage["png"]
    save(fig, out_png)
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


# Run: one figure per stage
def main() -> int:
    # Per vintage: each year on its own island (2026-08-26)
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
