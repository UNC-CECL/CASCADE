r"""
HAT_road_placement_on_domains.py
===============================================================================
Where each setback method actually puts NC-12 on the Barrier3D interiors CASCADE
initialises with -- the road placed exactly where `roadway_manager.bulldoze`
would put it from that method's model-facing RoadSetback CSV.

  a  1984 road on the 1984-start interiors (product+version resolved at run time)
  b  2004 road on the 2004-start interiors -- a DIFFERENT island, not the same
     one twice: 65 of 90 domains differ in interior shape between the products
  c  setback against island width, one line per period
  d  road movement between the two periods, under this method
  e  bulldoze's own drown test -- does CASCADE keep managing this roadway

ONE SCRIPT, BOTH METHODS, ON PURPOSE
------------------------------------
Every method in METHODS is drawn by this same code and each figure lands beside
its own data. Two scripts would drift, and a drifted comparison is worse than no
comparison -- it looks like a result. Add a method by adding a dict entry, not
by copying this file.

  old        old_method_offset/<year>/RoadSetback_<year>.csv
  dunestart  dunestart_offset/<year>/RoadSetback_<year>_dunestart.csv

Read-only with respect to the forcing: writes a PNG, a PDF beside it and a
CAPTIONS.md entry into each method's folder.

THE STYLE IS THE HOUSE STYLE, AND THIS MODULE RE-EXPORTS IT
-----------------------------------------------------------
Since 2026-09-10 every colour and every type size here comes from
scripts/hat_figure_style.py, through apply_style(). The palette names this file
used to own -- LAND_CMAP, SURFACE, WATER, INK_MUTED, INK_SECOND, C_1984 /
C_2004 / C_YEAR, C_DROWN, SECTIONS -- are KEPT and now resolve to their house
equivalents, because HAT_dunestart_modification_stages.py,
HAT_method_comparison_figures.py, HAT_oceanfloor_offset_check.py and
HAT_road_geojson_map.py all read them off this module as `P.<name>`.

THE FRAME CAVEAT -- WHICH APPLIES TO ONE METHOD AND NOT THE OTHER
------------------------------------------------------------------
This figure always draws the setback landward of INTERIOR ROW 0 of that
period's own extraction, because that is the reference CASCADE applies it
against
(`roadway_manager.py:99`):

    road_start = int(road_setback / dy)      rows landward of interior row 0

For `dunestart` that is also the frame the number was MEASURED in, so drawing
and measuring agree.

For `old` they do not: that method measured against the same-year digitised dune
line, a different feature from a different year. That mismatch is the REFERENCE
component isolated in 2-audit/HAT_road_method_diagnostic.py. It is not a reason
to draw the figure differently -- a mis-referenced setback still lands wherever
CASCADE lands it -- but it is the reason the two figures differ, and it should be
read as a property of the method rather than of the island.

HOW THE ROAD IS PLACED  (transcribed from bulldoze, not approximated)
---------------------------------------------------------------------
  road_start = int(setback / 10 m)                  truncation, not rounding
  band       = rows road_start .. road_start + 1    20 m, ROAD_WIDTH_CELLS = 2

The same int() is applied to every one of the 50 alongshore profiles in a
domain, so the road is a straight line across an island whose back edge is not.

WHAT THE INTERIOR ARRAYS ARE
----------------------------
domain_<d>_topography.npy, shape (rows, 50), values in DECAMETRES relative
to MHW. -0.30 dam (= -3 m) is the extractor's water sentinel, not an elevation.
Multiplied by 10 for display; masked, not ramped, where it is sentinel -- water
is a state, not a small magnitude.

REQUIREMENTS
------------
  numpy, matplotlib
===============================================================================
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

# =============================================================================
# CONFIG
# =============================================================================

PROJECT_ROOT = next(_p for _p in Path(__file__).resolve().parents
                    if (_p / "pyproject.toml").exists())
INIT_ROOT = PROJECT_ROOT / "data" / "hatteras_init"

# The topography version is NOT hardcoded here any more. It was "2009_v3", and
# when the dune windows were re-picked into 2009_v4 this kept drawing v3
# interiors under v4 setbacks without erroring. See hat_topo_version.py.
# parents[4] IS scripts/ -- hat_topo_version.py moved there 2026-08-20.
sys.path.insert(0, str(Path(__file__).resolve().parents[4]))
from hat_topo_version import (topo_dirs, array_name,  # noqa: E402
                             product_for_year)

# The house style. Everything typographic and every colour comes from here now;
# nothing in this file re-decides a font size or a grey. The names this module
# used to define for its own palette survive below as aliases, because three
# other scripts import this file for them.
import hat_figure_style as _HS  # noqa: E402
from hat_figure_style import (  # noqa: E402
    C, DOMAIN_AXIS_LABEL, apply_style, caption, elevation_cmap, figsize,
    open_frame, save, spines_for_image, town_bands, _title)

apply_style()

# NOR IS THE PRODUCT (2026-08-26). Until today this line read
#
#     TOPO_DIR, DUNE_DIR, TOPO_RUN_NAME = topo_dirs()
#
# at module level -- one topography, resolved once, drawn under BOTH panels of
# a two-vintage figure. That was correct while a single extraction served every
# period. It stopped being correct when the tree went period-first, and the
# figure said so in its own caption: "the SAME interiors back both road panels,
# since one topography serves every period."
#
# They no longer do. All 90 domains differ between 1984-start and 2004-start
# and 65 have a different interior SHAPE, so the 1984 panel was drawing a 1984
# setback -- measured from row 0 of the 1984-start extraction -- against a
# 2004-start island. Same class of error as the v3/v4 one above, and just as
# silent: the picture is plausible, it is simply of somewhere else.
#
# Resolved per vintage now, through the one YEAR_PRODUCT mapping.
_TOPO_CACHE: dict[int, tuple] = {}


def topo_for_year(year: int):
    """(topography dir, dunes dir, version) for one road vintage, cached."""
    if year not in _TOPO_CACHE:
        _TOPO_CACHE[year] = topo_dirs(product_for_year(year))
    return _TOPO_CACHE[year]


def topo_label(year: int) -> str:
    """`1984-start/v1`, for captions and stdout."""
    return f"{product_for_year(year)}/{topo_for_year(year)[2]}"


# array_name() is the single definition of these filenames - the same one
# the extractor writes with. Nothing here spells a name.
ROADS_ROOT = INIT_ROOT / "4-mgmt-forcing" / "road_offset"

# Each entry produces one figure, beside that method's own data.
#   setback  the MODEL-FACING file -- the one CASCADE reads, already floored
#            where the method floors it, because the question is where the model
#            puts the road, not what the method measured before clamping.
#   detail   optional per-domain CSV carrying a `flags` column, used only to
#            mark domains whose true value was negative and got floored to 0.
METHODS = {
    "old": dict(
        label=("setback taken as the minimum road elevation minus the "
               "minimum dune elevation, independently per domain "
               "(superseded)"),
        short="independent minima",
        root=ROADS_ROOT / "old_method_offset",
        setback="{year}/RoadSetback_{year}.csv",
        detail=None,
        png="HAT_old_method_road_on_domains.png",
    ),
    "dunestart": dict(
        label="setback measured landward from the dune start",
        short="dune start",
        root=ROADS_ROOT / "dunestart_offset",
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

# --- the drown test, transcribed from roadway_manager.bulldoze --------------
# A roadway width-drowns when water cells BORDER it. Three details matter and
# none of them is what you would guess from "the road is in water":
#
#   * the rows tested are the NEIGHBOURS of the bulldozed band, not the band --
#     road_start - 1 (seaside) and road_end + 1 (bayside);
#   * "water" is elevation <= 0 m MHW (drown_threshold=0), NOT the extractor's
#     -3 m sentinel. Real land that merely sits below MHW counts;
#   * it fires if EITHER side exceeds the fraction, strictly greater than.
#
# roadway_manager.py:50-51 and 125-165. RoadwayManager re-asserts 0.2 at :531.
DROWN_THRESHOLD_M = 0.0
DROWN_PCT = 0.2

# --- the palette, re-pointed at the house one (2026-09-10) -------------------
# These names are the module's public palette: HAT_dunestart_modification_
# stages.py, HAT_method_comparison_figures.py, HAT_oceanfloor_offset_check.py
# and HAT_road_geojson_map.py all read them off this module as `P.<name>`. They
# are kept, and each one now RESOLVES to its house equivalent, so the figures in
# this folder and the figures elsewhere in the project are one palette.
#
# The earlier vintage is the house red and the later the house blue, which is
# the vintage pair everywhere in this project; the pair used here before
# (#2a78d6 blue 1984 / #eb6834 orange 2004) had 1984 blue, i.e. exactly
# inverted against every other two-vintage figure in the repo.
C_1984, C_2004 = _HS.C_1984, _HS.C_1997
C_YEAR = {1984: C_1984, 2004: C_2004}

# The drowned state. It was a dark crimson picked to separate from the old
# blue/orange pair; against the house vintage RED it would now read as "1984".
# ACCENT purple is the one house colour that is neither vintage, neither
# reference nor fabricated ground, and it separates from both poles on hue and
# from BASE on luminance. Drowned domains also carry a marker and a count, so
# the state is never colour-alone.
C_DROWN = C["ACCENT"]

# INK_SECOND is used for TEXT by the importers, INK_MUTED for rules, outlines
# and grid; the house rule puts text in INK and rules in INK_MUTED.
INK_MUTED, INK_SECOND = _HS.INK_MUTED, _HS.INK
SURFACE = "white"                 # halo / label backing; the house page colour
WATER = C["WATER"]                # cells at or below MHW
NODATA = C["BASE_FILL"]           # outside the extraction: no data, not water

# Elevation is drawn in CLASSES, not a ramp (house rule): a linear ramp over
# 0-4 m renders the whole back-barrier as one tone. LAND_CLASS_* is what the
# panels here use. LAND_CMAP survives as a CONTINUOUS ramp in the same house
# colours, because the two out-of-scope importers pair it with a plain
# Normalize(LAND_VMIN, LAND_VMAX); handing them the discrete list under a
# linear norm would paint 0-0.5 m land in the water colour.
LAND_CLASS_CMAP, LAND_CLASS_NORM, LAND_CLASS_BOUNDS = elevation_cmap()
LAND_CMAP = LinearSegmentedColormap.from_list(
    "hat_land_ramp", list(LAND_CLASS_CMAP.colors)[1:])
LAND_VMIN, LAND_VMAX = 0.0, 4.0

# The alongshore reaches. Kept because HAT_oceanfloor_offset_check.py reads
# them off this module to draw its own dividers. The figures BELOW no longer
# use them: the house `town_bands()` shades the three village spans from
# hatteras_site_config, which is the one authority on where the towns are.
SECTIONS = [((1, 6), "Cape Pt"), ((7, 8), "Bux"), ((9, 20), "Buxton-Avon"),
            ((21, 31), "Avon"), ((32, 67), "Avon-Tri-Village / Wimble Shoals"),
            ((68, 83), "Tri-Village"), ((84, 90), "Pea Is.")]


# =============================================================================
# LOAD
# =============================================================================

def read_two_row(path: Path) -> dict:
    if not path.is_file():
        return {}
    raw = np.loadtxt(path, delimiter=",")
    if raw.ndim != 2 or raw.shape[0] != 2:
        return {}
    return {int(k): float(v) for k, v in zip(raw[0], raw[1])}


def load_interiors(year: int) -> dict:
    """The Barrier3D interiors THIS vintage's setbacks were measured against.

    `year` is required. It used to be absent, and a caller that cannot name a
    year is a caller drawing two vintages on one island -- see topo_for_year.
    """
    topo = topo_for_year(year)[0]
    out = {}
    for d in DOMAINS:
        p = topo / array_name("topography", d)
        if p.is_file():
            out[d] = np.load(p)
    if not out:
        raise SystemExit(f"no topography found in {topo}")
    return out


def load_years(years=YEARS) -> tuple[dict, int, int]:
    """{year: {interiors, canvas, shown}}, plus the crop shared by all panels.

    The panels are stacked and share an x axis, so they must share a y extent
    too -- otherwise a taller 1984 island would read as a wider barrier rather
    than a taller canvas. The crop is therefore the max over vintages, and a
    shorter vintage's canvas is NaN-padded up to it rather than drawn shorter.
    """
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


def _wet_fraction(row: np.ndarray) -> float:
    """
    Fraction of a border row at or below the threshold -- EVERY cell counted.

    A no-data mask (domain_<d>_nodata.npy) exists beside the topography and
    was briefly used here to drop never-surveyed cells from both the numerator
    and the denominator. That is not the test bulldoze runs. bulldoze reads the
    literal array, where the extractor has already written no-data back as the
    water sentinel, so an unsurveyed cell IS a water cell as far as the model is
    concerned -- and a road CASCADE would stop managing has to show as one here.
    Excluding no-data flipped GIS 14/77/78/79/80 to pass while the model still
    drowns them. Deliberately reverted: the figure reports what the run does.
    """
    return float((row * 10.0 <= DROWN_THRESHOLD_M).mean())


def build_canvas(interiors: dict) -> tuple:
    max_rows = max(a.shape[0] for a in interiors.values())
    canvas = np.full((max_rows, len(DOMAINS) * ALONG_COLS), np.nan)
    for d, a in interiors.items():
        c0 = (d - 1) * ALONG_COLS
        canvas[:a.shape[0], c0:c0 + a.shape[1]] = np.where(
            a > SENTINEL_DAM + 1e-6, a * 10.0, np.nan)
    return canvas, max_rows


def place_road(interiors: dict, setbacks: dict) -> dict:
    """
    bulldoze's own placement, and bulldoze's own drown test, per domain.

    The drown test is transcribed rather than approximated -- it decides whether
    CASCADE gives up managing this roadway, so an approximation of it would be a
    figure about a model we are not running. Interior values are decametres MHW
    and the model compares `grid * dz` against the threshold in metres, so the
    * 10 here is bulldoze's, not a display convenience.
    """
    out = {}
    for d, sb in sorted(setbacks.items()):
        a = interiors.get(d)
        if a is None:
            continue
        n, ncols = a.shape
        start = int(sb / CELL_SIZE_M)            # truncation, as bulldoze does
        end = start + ROAD_WIDTH_CELLS           # exclusive

        # bulldoze indexes road_end + 1 with no bounds check, so a road this
        # far back is not "drowned", it is an IndexError at t=0.
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


# =============================================================================
# PANELS
# =============================================================================

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

    # Only road_start is drawn. The 20 m band's landward edge sat 2 cells away,
    # which at this vertical scale read as line weight, not as width.
    #
    # Segments are split by drown state so a failing domain is crimson in place,
    # rather than being annotated off to one side -- the question "where does
    # this road fail" is answered on the island, not in a caption.
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

    # No marker on the plan views: the accent segment IS the signal, and a
    # triangle on top of a 1-domain-wide line obscured the thing it pointed at.
    # The state stays recoverable without colour -- panel (e) names every
    # failing domain on the same x-axis, and the caption carries the count.

    # The three village spans, from hatteras_site_config -- a named strip
    # against the landward edge rather than a full-height wash, which on an
    # image panel would hide the island it is meant to locate.
    ax.set_xlim(0.5, len(DOMAINS) + 0.5)
    if label_towns:
        town_bands(ax, strip=0.075, shade=SURFACE)

    spines_for_image(ax)
    _title(ax, panel_index, f"NC-12 in {year}")
    ax.set_ylabel("m landward of\ninterior row 0")
    plt.setp(ax.get_xticklabels(), visible=False)


def read_floored(spec: dict, year: int) -> set:
    """Domains whose true setback was negative and got floored to 0."""
    if not spec.get("detail"):
        return set()
    p = spec["root"] / spec["detail"].format(year=year)
    if not p.is_file():
        return set()
    import csv as _csv
    with open(p, newline="") as f:
        return {int(r["domain"]) for r in _csv.DictReader(f)
                if "NEGATIVE" in (r.get("flags") or "")}


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
        # Each vintage is placed on ITS OWN interiors. Passing one shared dict
        # here is the bug this signature exists to make impossible.
        placed[year] = place_road(per[year]["interiors"], sb)
        floored[year] = read_floored(spec, year)
        print(f"  {year}: {topo_label(year)} interiors")

    if not placed:
        print(f"  [skip] {name}: no setback files found")
        return

    # Anything cropped out of the plan view that is actually road would make the
    # picture a lie, so check rather than assume.
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

    # --- (c) setback against the island it has to fit inside ----------------
    # ONE BAND PER VINTAGE. This was a single shaded band, drawn from the one
    # shared interiors dict, with both setback curves over it -- which read as
    # "here is the island, here are two roads on it". There are two islands.
    # The 1984 and 2004 widths differ by up to 80 m in places, and a setback
    # that fits one can run off the other.
    for year in YEARS:
        if year not in per:
            continue
        interiors = per[year]["interiors"]
        width_x = sorted(interiors)
        width_y = [interiors[d].shape[0] * CELL_SIZE_M for d in width_x]
        first = year == YEARS[0]
        # Lines, not a shaded band under each width. Two 10%-alpha bands, one
        # red and one blue, overprint to a purple wash across the whole panel
        # -- and purple is the drowned colour everywhere else in this figure.
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

    # --- (d) where the road moved between the two periods -------------------
    # Ported from the retired 3-figures/island_wide/HAT_plot_road_on_b3d_domains
    # .py, which drew this for one method only.
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

    # --- (e) what the bulldozed band actually lands on ----------------------
    # The series is a PERCENTAGE, so the threshold has to be scaled too -- at
    # DROWN_PCT it would sit on 0.2% and read as zero.
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

    # --- legend and caption -------------------------------------------------
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
        # ASCII hyphen, not U+2212: the Windows console is cp1252 and a
        # unicode minus raises UnicodeEncodeError here. Figure text is fine.
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

    # The comparison the two figures exist to support, stated as numbers so it
    # does not have to be eyeballed off two PNGs.
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
