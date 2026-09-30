"""
Where each setback method actually puts NC-12 on the Barrier3D interiors CASCADE starts from.

    python scripts/input_prep/4-mgmt-forcings/road_offset/1-produce/HAT_road_placement_on_domains.py

The road placed as roadway_manager.bulldoze would place it from each method's
setback CSV, with the drown test, one figure per method; other scripts import
its loaders and palette. Details: scripts/input_prep/4-mgmt-forcings/road_offset/README.md.

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-09-18
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.lines import Line2D
from matplotlib.patheffects import withStroke


PROJECT_ROOT = next(_p for _p in Path(__file__).resolve().parents
                    if (_p / "pyproject.toml").exists())
INIT_ROOT = PROJECT_ROOT / "data" / "hatteras_init"

# The topography version is NOT hardcoded here any more
sys.path.insert(0, str(Path(__file__).resolve().parents[4]))
from site_layer.hat_topo_version import (topo_dirs, array_name,  # noqa: E402
                             product_for_year)

# House style from site_layer/hat_figure_style.py; the old palette names stay as aliases
from site_layer import hat_figure_style as _HS  # noqa: E402
from site_layer.hat_figure_style import (  # noqa: E402
    C, DOMAIN_AXIS_LABEL, apply_style, caption, elevation_cmap, figsize,
    open_frame, save, spines_for_image, town_bands, _title)

apply_style()

# NOR IS THE PRODUCT (2026-08-26)
_TOPO_CACHE: dict[int, tuple] = {}


# (topography dir, dunes dir, version) for one road vintage, cached
def topo_for_year(year: int):
    if year not in _TOPO_CACHE:
        _TOPO_CACHE[year] = topo_dirs(product_for_year(year))
    return _TOPO_CACHE[year]


# `1984-start/v1`, for captions and stdout
def topo_label(year: int) -> str:
    return f"{product_for_year(year)}/{topo_for_year(year)[2]}"


# Filenames from array_name(), the same definition the extractor writes with
from site_layer import hat_topo_version as _tv  # noqa: E402
# --- CONFIG ------------------------------------------------------------------
ROADS_ROOT = _tv.ROADS_ROOT

# Each entry produces one figure, beside that method's own data
METHODS = {
    "old": dict(
        label=("setback taken as the minimum road elevation minus the "
               "minimum dune elevation, independently per domain "
               "(superseded)"),
        short="independent minima",
        root=_tv.LEGACY_SETBACK_ROOT,   # was old_method_offset/ (2026-09-18)
        setback="{year}/RoadSetback_{year}.csv",
        detail=None,
        png="HAT_old_method_road_on_domains.png",
    ),
    "dunestart": dict(
        label="setback measured landward from the dune start",
        short="dune start",
        # measured/ only: the two measured starts
        root=ROADS_ROOT / "dunestart_offset" / "measured",
        setback="{year}/RoadSetback_{year}_dunestart.csv",
        detail="{year}/RoadOffset_{year}_domains.csv",
        png="HAT_dunestart_road_on_domains.png",
    ),
}

YEARS = (1984, 2004)
DOMAINS = list(range(1, 91))
ALONG_COLS = 50
CELL_SIZE_M = 10.0
SENTINEL_DAM = -0.3          # extractor's water sentinel, decametres MHW
ROAD_WIDTH_CELLS = 2         # bulldoze's 20 m band: int(road_width 20 / dx 10)
DISPLAY_CROSS_SHORE_M = 900.0

# The drown test, transcribed from roadway_manager.bulldoze

# A roadway width-drowns when water cells BORDER it
DROWN_THRESHOLD_M = 0.0
DROWN_PCT = 0.2

# The palette, re-pointed at the house one (2026-09-10)

# These names are the module's public palette
C_1984, C_2004 = _HS.C_1984, _HS.C_1997
C_YEAR = {1984: C_1984, 2004: C_2004}

# The drowned state: ACCENT purple, neither vintage; drowned domains also get a count
C_DROWN = C["ACCENT"]

# INK_SECOND is used for TEXT by the importers, INK_MUTED for rules, outlines and grid
INK_MUTED, INK_SECOND = _HS.INK_MUTED, _HS.INK
SURFACE = "white"                 # halo / label backing; the house page colour
WATER = C["WATER"]                # cells at or below MHW
NODATA = C["BASE_FILL"]           # outside the extraction: no data, not water

# Elevation is drawn in CLASSES, not a ramp (house rule)
LAND_CLASS_CMAP, LAND_CLASS_NORM, LAND_CLASS_BOUNDS = elevation_cmap()
LAND_CMAP = LinearSegmentedColormap.from_list(
    "hat_land_ramp", list(LAND_CLASS_CMAP.colors)[1:])
LAND_VMIN, LAND_VMAX = 0.0, 4.0

# The alongshore reaches, kept for HAT_oceanfloor_offset_check.py
SECTIONS = [((1, 6), "Cape Pt"), ((7, 8), "Bux"), ((9, 20), "Buxton-Avon"),
            ((21, 31), "Avon"), ((32, 67), "Avon-Tri-Village / Wimble Shoals"),
            ((68, 83), "Tri-Village"), ((84, 90), "Pea Is.")]
# -----------------------------------------------------------------------------


# Load

# {domain: value} from a two-row CSV
def read_two_row(path: Path) -> dict:
    if not path.is_file():
        return {}
    raw = np.loadtxt(path, delimiter=",")
    if raw.ndim != 2 or raw.shape[0] != 2:
        return {}
    return {int(k): float(v) for k, v in zip(raw[0], raw[1])}


# The Barrier3D interiors THIS vintage's setbacks were measured against
def load_interiors(year: int) -> dict:
    topo = topo_for_year(year)[0]
    out = {}
    for d in DOMAINS:
        p = topo / array_name("topography", d)
        if p.is_file():
            out[d] = np.load(p)
    if not out:
        raise SystemExit(f"no topography found in {topo}")
    return out


# {year
def load_years(years=YEARS) -> tuple[dict, int, int]:
    per = {}
    for year in years:
        interiors = load_interiors(year)
        canvas, max_rows = build_canvas(interiors)
        per[year] = dict(interiors=interiors, canvas=canvas, max_rows=max_rows)

    max_rows = max(p["max_rows"] for p in per.values())
    crop_rows = min(max_rows, int(DISPLAY_CROSS_SHORE_M / CELL_SIZE_M))
    for p in per.values():
        c = p["canvas"]
        if c.shape[0] < crop_rows:
            pad = np.full((crop_rows - c.shape[0], c.shape[1]), np.nan)
            c = np.vstack([c, pad])
        p["shown"] = c[:crop_rows, :]
    return per, crop_rows, max_rows


# Fraction of a border row at or below the threshold -- EVERY cell counted
def _wet_fraction(row: np.ndarray) -> float:
    return float((row * 10.0 <= DROWN_THRESHOLD_M).mean())


# Every domain's interior side by side on one canvas
def build_canvas(interiors: dict) -> tuple:
    max_rows = max(a.shape[0] for a in interiors.values())
    canvas = np.full((max_rows, len(DOMAINS) * ALONG_COLS), np.nan)
    for d, a in interiors.items():
        c0 = (d - 1) * ALONG_COLS
        canvas[:a.shape[0], c0:c0 + a.shape[1]] = np.where(
            a > SENTINEL_DAM + 1e-6, a * 10.0, np.nan)
    return canvas, max_rows


# bulldoze's own placement, and bulldoze's own drown test, per domain
def place_road(interiors: dict, setbacks: dict) -> dict:
    out = {}
    for d, sb in sorted(setbacks.items()):
        a = interiors.get(d)
        if a is None:
            continue
        n, ncols = a.shape
        start = int(sb / CELL_SIZE_M)            # truncation, as bulldoze does
        end = start + ROAD_WIDTH_CELLS           # exclusive

        # Bulldoze indexes road_end + 1 with no bounds check, so a road this far back is not "drowned"
        crashes = end + 1 >= n
        if crashes:
            sea = bay = np.nan
            drowned = True
        else:
            bay = _wet_fraction(a[end + 1, :])
            sea = _wet_fraction(a[start - 1, :]) if start > 0 else 0.0
            drowned = bool(sea > DROWN_PCT or bay > DROWN_PCT)

        out[d] = dict(setback_m=sb, start_m=start * CELL_SIZE_M,
                      island_m=n * CELL_SIZE_M,
                      headroom_m=n * CELL_SIZE_M - sb,
                      seaside=sea, bayside=bay, drowned=drowned,
                      crashes=crashes,
                      governing=(np.nan if crashes else max(sea, bay)))
    return out


# Panels

# One vintage's island with the placed road
def draw_island(ax, fig, shown, crop_rows, year, placed, panel_index,
                label_towns=False):
    colour = C_YEAR[year]
    ax.set_facecolor(NODATA)
    im = ax.imshow(np.ma.masked_invalid(shown), aspect="auto", origin="lower",
                   extent=[0.5, len(DOMAINS) + 0.5, -CELL_SIZE_M / 2,
                           crop_rows * CELL_SIZE_M - CELL_SIZE_M / 2],
                   cmap=LAND_CLASS_CMAP, norm=LAND_CLASS_NORM,
                   interpolation="nearest")
    cax = ax.inset_axes([1.012, 0.0, 0.014, 1.0])
    cb = fig.colorbar(im, cax=cax, spacing="uniform",
                      ticks=LAND_CLASS_BOUNDS[1:-1])
    cb.set_label("elevation (m MHW)")
    cb.outline.set_edgecolor(INK_MUTED)
    cb.outline.set_linewidth(0.6)

    # Only road_start is drawn, split by drown state so a failure shows in place
    halo = [withStroke(linewidth=4.4, foreground=SURFACE)]
    seg = {False: ([], []), True: ([], [])}
    for d, p in placed.items():
        xs, ys = seg[bool(p["drowned"])]
        xs += [d - 0.5, d + 0.5, np.nan]
        ys += [p["start_m"], p["start_m"], np.nan]

    ax.plot(seg[False][0], seg[False][1], color=colour, lw=2.6,
            solid_capstyle="butt", path_effects=halo, zorder=6)
    if seg[True][0]:
        ax.plot(seg[True][0], seg[True][1], color=C_DROWN, lw=3.4,
                solid_capstyle="butt", path_effects=halo, zorder=7)

    # No marker on the plan views: the accent segment is the signal

    # The three village spans, from hatteras_site_config
    ax.set_xlim(0.5, len(DOMAINS) + 0.5)
    if label_towns:
        town_bands(ax, strip=0.075, shade=SURFACE)

    spines_for_image(ax)
    _title(ax, panel_index, f"NC-12 in {year}")
    ax.set_ylabel("m landward of\ninterior row 0")
    plt.setp(ax.get_xticklabels(), visible=False)


# Domains whose true setback was negative and got floored to 0
def read_floored(spec: dict, year: int) -> set:
    if not spec.get("detail"):
        return set()
    p = spec["root"] / spec["detail"].format(year=year)
    if not p.is_file():
        return set()
    import csv as _csv
    with open(p, newline="") as f:
        return {int(r["domain"]) for r in _csv.DictReader(f)
                if "NEGATIVE" in (r.get("flags") or "")}


# One method's figure, both vintages
def build_figure(name: str, spec: dict, per: dict, crop_rows,
                 max_rows) -> None:
    print(f"\n{'=' * 88}")
    print(f"{name.upper()} -- where the road lands on the Barrier3D interiors")
    print("=" * 88)

    placed, floored = {}, {}
    for year in YEARS:
        sb = read_two_row(spec["root"] / spec["setback"].format(year=year))
        if not sb:
            print(f"  [skip] {year}: no RoadSetback CSV under {spec['root']}")
            continue
        # Each vintage is placed on ITS OWN interiors
        placed[year] = place_road(per[year]["interiors"], sb)
        floored[year] = read_floored(spec, year)
        print(f"  {year}: {topo_label(year)} interiors")

    if not placed:
        print(f"  [skip] {name}: no setback files found")
        return

    # Anything cropped out of the plan view that is actually road would make the picture a lie
    deepest = max(p["start_m"] + ROAD_WIDTH_CELLS * CELL_SIZE_M
                  for pl in placed.values() for p in pl.values())
    if deepest > crop_rows * CELL_SIZE_M:
        print(f"  [warn] road reaches {deepest:.0f} m but the panels crop at "
              f"{crop_rows * CELL_SIZE_M:.0f} m -- raise DISPLAY_CROSS_SHORE_M")
    else:
        print(f"  deepest road band {deepest:.0f} m, panels crop at "
              f"{crop_rows * CELL_SIZE_M:.0f} m -- nothing cropped is road")

    fig = plt.figure(figsize=figsize("double", height=9.4))
    gs = fig.add_gridspec(5, 1, height_ratios=[1.25, 1.25, 0.95, 0.72, 0.72],
                          hspace=0.30, left=0.105, right=0.870,
                          top=0.962, bottom=0.078)
    ax84 = fig.add_subplot(gs[0])
    ax04 = fig.add_subplot(gs[1], sharex=ax84)
    ax_sb = fig.add_subplot(gs[2], sharex=ax84)
    ax_mv = fig.add_subplot(gs[3], sharex=ax84)
    ax_w = fig.add_subplot(gs[4], sharex=ax84)

    for i, (ax, year) in enumerate(((ax84, YEARS[0]), (ax04, YEARS[1]))):
        if year in placed:
            draw_island(ax, fig, per[year]["shown"], crop_rows, year,
                        placed[year], i, label_towns=(i == 0))

    # (c) setback against the island it has to fit inside

    # One band per vintage: the 1984 and 2004 islands differ by up to 80 m
    for year in YEARS:
        if year not in per:
            continue
        interiors = per[year]["interiors"]
        width_x = sorted(interiors)
        width_y = [interiors[d].shape[0] * CELL_SIZE_M for d in width_x]
        first = year == YEARS[0]
        # Lines, not a shaded band under each width
        ax_sb.plot(width_x, width_y, color=C_YEAR[year], lw=0.9,
                   ls="-" if first else (0, (4, 2)), alpha=0.8, zorder=2,
                   label=f"{year} island width")
    for year, pl in placed.items():
        xs = sorted(pl)
        ax_sb.plot(xs, [pl[d]["setback_m"] for d in xs], color=C_YEAR[year],
                   lw=1.8 if year == YEARS[0] else 1.3, zorder=6,
                   label=f"{year} setback")
    ax_sb.set_ylabel("m landward of\ninterior row 0")
    ax_sb.grid(axis="y")
    ax_sb.set_axisbelow(True)
    ax_sb.set_ylim(0, min(DISPLAY_CROSS_SHORE_M, max(width_y) * 1.05))
    town_bands(ax_sb, label=False)
    open_frame(ax_sb)
    ax_sb.legend(loc="upper left", ncol=2, fontsize=7)
    _title(ax_sb, 2, "setback against island width")
    plt.setp(ax_sb.get_xticklabels(), visible=False)

    # (d) where the road moved between the two periods

    # Ported from the retired 3-figures/island_wide/HAT_plot_road_on_b3d_domains .py
    move_median = None
    if len(placed) == 2:
        ya, yb = sorted(placed)
        common = sorted(set(placed[ya]) & set(placed[yb]))
        move = np.array([placed[yb][d]["setback_m"] - placed[ya][d]["setback_m"]
                         for d in common])
        move_median = float(np.median(move))
        ax_mv.axhline(0, color=INK_MUTED, lw=0.8, zorder=3)
        ax_mv.bar(common, np.where(move >= 0, move, 0.0), width=0.86,
                  color=C_YEAR[yb], linewidth=0, zorder=5)
        ax_mv.bar(common, np.where(move < 0, move, 0.0), width=0.86,
                  color=C_YEAR[ya], linewidth=0, zorder=5)
    ax_mv.set_ylabel(f"{max(placed)} \u2212 {min(placed)}\nsetback (m)")
    ax_mv.grid(axis="y")
    ax_mv.set_axisbelow(True)
    town_bands(ax_mv, label=False)
    open_frame(ax_mv)
    _title(ax_mv, 3, "change in setback between the periods")
    plt.setp(ax_mv.get_xticklabels(), visible=False)

    # (e) what the bulldozed band actually lands on

    # The series is a PERCENTAGE, so the threshold has to be scaled too
    ax_w.axhspan(DROWN_PCT * 100, 104, color=C_DROWN, alpha=0.07, lw=0,
                 zorder=1)
    ax_w.axhline(DROWN_PCT * 100, color=C_DROWN, lw=1.1, ls=(0, (4, 3)),
                 zorder=3, label=f"threshold, {DROWN_PCT * 100:.0f}%")
    for year, pl in sorted(placed.items()):
        xs = sorted(pl)
        ax_w.plot(xs, [pl[d]["governing"] * 100 for d in xs],
                  color=C_YEAR[year], lw=1.2, marker="o", ms=2.0, zorder=5,
                  label=f"{year}")
        bad = [d for d in xs if pl[d]["drowned"]]
        if bad:
            ax_w.plot(bad, [pl[d]["governing"] * 100 for d in bad], lw=0,
                      marker="v", ms=6, mfc=C_DROWN, mec=SURFACE, mew=0.8,
                      zorder=7)
    ax_w.set_ylabel("% of bordering cells\nat or below 0 m MHW")
    ax_w.set_xlabel(DOMAIN_AXIS_LABEL)
    ax_w.set_xlim(0.5, len(DOMAINS) + 0.5)
    ax_w.set_ylim(-4, 104)
    ax_w.grid(axis="y")
    ax_w.set_axisbelow(True)
    town_bands(ax_w, label=False)
    open_frame(ax_w)
    ax_w.legend(loc="upper left", ncol=3, fontsize=7)
    _title(ax_w, 4, "wet cells bordering the road")

    # Legend and caption
    fig.legend(handles=[
        Line2D([], [], color=C_1984, lw=2.6, label="NC-12 in 1984"),
        Line2D([], [], color=C_2004, lw=2.6, label="NC-12 in 2004"),
        Line2D([], [], color=C_DROWN, lw=3.4,
               label="drowns at initialisation"),
        Line2D([], [], color=NODATA, lw=8, label="outside the extraction"),
    ], loc="lower center", bbox_to_anchor=(0.5, -0.004), ncol=4, frameon=False,
        columnspacing=1.6, handlelength=2.4)

    frame_note = (
        "The setback was measured in the same frame it is drawn in, interior "
        "row 0, so the drawing and the measurement agree."
        if name == "dunestart" else
        "The setback was measured against the same-year digitised dune line "
        "but CASCADE applies it landward of interior row 0, so the drawing "
        "and the measurement are in different frames.")
    n_floored = sum(len(v) for v in floored.values())
    floor_note = (f" {n_floored} domain-year(s) had a negative true setback and "
                  f"were floored to 0, putting the road on interior row 0."
                  if n_floored else "")
    move_note = ("" if move_median is None else
                 f" Median change between the periods {move_median:+.0f} m; "
                 f"bars above zero are further inland by {max(placed)}, below "
                 f"zero closer to the dune.")
    drown_note = "; ".join(
        f"{sum(1 for p in pl.values() if p['drowned'])} of {len(pl)} in {year}"
        for year, pl in sorted(placed.items()))

    caption(fig, (
        f"NC-12 placed on the Barrier3D interiors CASCADE initialises with, "
        f"from the {spec['label']}. Domain 1 is at Cape Point in the south and "
        f"domain 90 at Pea Island in the north; the shaded spans are the "
        f"villages (Buxton, Avon, Tri-Village). (a, b) the road as "
        f"roadway_manager.bulldoze places it, road_start = int(setback / "
        f"{CELL_SIZE_M:.0f} m), on each period's OWN extraction \u2014 "
        + ", ".join(f"{y} on {topo_label(y)}" for y in YEARS if y in per)
        + ". These are different islands: 65 of 90 domains differ in interior "
          "shape, so the two panels are not one island drawn twice. Interior "
          "elevation is shown in classes relative to mean high water; cells "
          "outside the extraction carry no data and are drawn grey. "
          "(c) the same setback against the island width it has to fit "
          "inside. (d) the change in setback between the two periods."
        + move_note +
        " (e) bulldoze's own drown test: the wetter of the two rows BORDERING "
        f"the bulldozed band (road_start \u2212 1, road_end + 1). Above "
        f"{DROWN_PCT * 100:.0f}% of bordering cells at or below 0 m MHW "
        f"CASCADE stops managing the roadway, and the road is drawn in the "
        f"accent colour wherever that happens \u2014 {drown_note}. "
        + frame_note + floor_note))

    out_png = spec["root"] / spec["png"]
    save(fig, out_png)
    plt.close(fig)
    print(f"\n[out] {out_png}")

    for year, pl in sorted(placed.items()):
        drowned = [d for d, p in pl.items() if p["drowned"]]
        crash = [d for d, p in pl.items() if p["crashes"]]
        print(f"\n  {year}: {len(pl)} domains placed | "
              f"{len(drowned)} DROWN at initialisation")
        if crash:
            print(f"    [warn] road_end+1 beyond the array in {crash} -- "
                  f"bulldoze would raise IndexError, not drown")
        if drowned:
            print(f"    {'GIS':>4} {'seaside':>8} {'bayside':>8}   "
                  f"(fails above {DROWN_PCT:.2f} on either side)")
            for d in drowned:
                p = pl[d]
                sea = " n/a" if np.isnan(p["seaside"]) else f"{p['seaside']:.2f}"
                bay = " n/a" if np.isnan(p["bayside"]) else f"{p['bayside']:.2f}"
                side = ("both" if p["seaside"] > DROWN_PCT
                        and p["bayside"] > DROWN_PCT else
                        "bayside" if p["bayside"] > DROWN_PCT else "seaside")
                print(f"    {d:>4} {sea:>8} {bay:>8}   {side}")
        # ASCII hyphen: the Windows console is cp1252
        print("    least headroom (island width - setback):")
        for d in sorted(pl, key=lambda d: pl[d]["headroom_m"])[:3]:
            print(f"      GIS {d:>2}: setback {pl[d]['setback_m']:>5.0f} m | "
                  f"island {pl[d]['island_m']:>5.0f} m | "
                  f"headroom {pl[d]['headroom_m']:>5.0f} m")

    if len(placed) == 2:
        ya, yb = sorted(placed)
        a = {d for d, p in placed[ya].items() if p["drowned"]}
        b = {d for d, p in placed[yb].items() if p["drowned"]}
        print(f"\n  drowned in both years : {sorted(a & b)}")
        print(f"  {ya} only              : {sorted(a - b)}")
        print(f"  {yb} only              : {sorted(b - a)}")
    return {y: {d for d, p in pl.items() if p["drowned"]}
            for y, pl in placed.items()}


# Run: load both vintages, one figure per method
def main() -> int:
    per, crop_rows, max_rows = load_years()
    for year in YEARS:
        interiors = per[year]["interiors"]
        print(f"{year}: {len(interiors)} interiors from {topo_label(year)} | "
              f"island width "
              f"{min(a.shape[0] for a in interiors.values()) * 10}"
              f"-{per[year]['max_rows'] * 10} m")

    drowned = {}
    for name, spec in METHODS.items():
        got = build_figure(name, spec, per, crop_rows, max_rows)
        if got:
            drowned[name] = got

    # The comparison the two figures exist to support
    if len(drowned) == 2:
        (na, da), (nb, db) = drowned.items()
        print(f"\n{'=' * 88}")
        print("DROWN AT INITIALISATION -- method against method")
        print("=" * 88)
        for year in sorted(set(da) & set(db)):
            a, b = da[year], db[year]
            print(f"  {year}: {na} {len(a):>2} | {nb} {len(b):>2} | "
                  f"both {len(a & b):>2}")
            print(f"    {na} only      : {sorted(a - b)}")
            print(f"    {nb} only      : {sorted(b - a)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
