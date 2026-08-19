r"""
HAT_road_placement_on_domains.py
===============================================================================
Where each setback method actually puts NC-12 on the Barrier3D interiors CASCADE
initialises with -- the road placed exactly where `roadway_manager.bulldoze`
would put it from that method's model-facing RoadSetback CSV.

  A  1984 road on the extracted interiors (version resolved at run time)
  B  2004 road on the SAME interiors
  C  setback against island width -- how much headroom each domain has
  D  road movement between the two periods, under this method
  E  bulldoze's own drown test -- does CASCADE keep managing this roadway

ONE SCRIPT, BOTH METHODS, ON PURPOSE
------------------------------------
Every method in METHODS is drawn by this same code and each figure lands beside
its own data. Two scripts would drift, and a drifted comparison is worse than no
comparison -- it looks like a result. Add a method by adding a dict entry, not
by copying this file.

  old        old_method_offset/<year>/RoadSetback_<year>.csv
  dunestart  dunestart_offset/<year>/RoadSetback_<year>_dunestart.csv

Read-only with respect to the forcing: writes only a PNG into each method's
folder.

THE FRAME CAVEAT -- WHICH APPLIES TO ONE METHOD AND NOT THE OTHER
------------------------------------------------------------------
This figure always draws the setback landward of INTERIOR ROW 0 of the 2009 DEM,
because that is the reference CASCADE applies it against
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
domain_<d>_topography_2009.npy, shape (rows, 50), values in DECAMETRES relative
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
from matplotlib.colors import LinearSegmentedColormap, Normalize
from matplotlib.lines import Line2D
from matplotlib.patheffects import withStroke

# =============================================================================
# CONFIG
# =============================================================================

PROJECT_ROOT = Path(__file__).resolve().parents[5]
INIT_ROOT = PROJECT_ROOT / "data" / "hatteras_init"

# The topography version is NOT hardcoded here any more. It was "2009_v3", and
# when the dune windows were re-picked into 2009_v4 this kept drawing v3
# interiors under v4 setbacks without erroring. See hat_topo_version.py.
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from hat_topo_version import topo_dirs  # noqa: E402

TOPO_DIR, DUNE_DIR, TOPO_RUN_NAME = topo_dirs()
ROADS_ROOT = INIT_ROOT / "4-mgmt-forcing" / "road_offset"

# Each entry produces one figure, beside that method's own data.
#   setback  the MODEL-FACING file -- the one CASCADE reads, already floored
#            where the method floors it, because the question is where the model
#            puts the road, not what the method measured before clamping.
#   detail   optional per-domain CSV carrying a `flags` column, used only to
#            mark domains whose true value was negative and got floored to 0.
METHODS = {
    "old": dict(
        label="OLD method (min road − min dune, digitised dune line)",
        short="old method",
        root=ROADS_ROOT / "old_method_offset",
        setback="{year}/RoadSetback_{year}.csv",
        detail=None,
        png="HAT_old_method_road_on_domains.png",
    ),
    "dunestart": dict(
        label="DUNE-START method (per-profile, measured from interior row 0)",
        short="dune-start method",
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

C_1984, C_2004 = "#2a78d6", "#eb6834"
C_YEAR = {1984: C_1984, 2004: C_2004}
# Dark crimson, not a bright red. Separation was computed, not eyeballed: the
# obvious #d62728 measures OKLab dE 11.8 against the 2004 orange, under the 15
# floor -- which is why an earlier figure in this repo rejected status red
# outright. This one clears it by going DARK rather than saturated: dE 26.7 vs
# orange, 32.7 vs blue, worst simulated CVD 26.4. Drowned domains also carry a
# marker and a label, so the state is never colour-alone.
C_DROWN = "#8c0d24"
INK_MUTED, INK_SECOND, SURFACE, WATER = "#8a8a85", "#52514e", "#fcfcfb", "#d8dce0"

LAND_CMAP = LinearSegmentedColormap.from_list(
    "land_seq", ["#eaf3e6", "#c3e0b9", "#8fc487", "#55a05e", "#217a45",
                 "#0d4d2b"])
LAND_VMIN, LAND_VMAX = 0.0, 4.0

SECTIONS = [((1, 6), "Cape Pt"), ((7, 8), "Bux"), ((9, 20), "Buxton-Avon"),
            ((21, 31), "Avon"), ((32, 67), "Avon-Tri-Village / Wimble Shoals"),
            ((68, 83), "Tri-Village"), ((84, 90), "Pea Is.")]

plt.rcParams.update({
    "font.size": 10,
    "axes.edgecolor": INK_MUTED, "axes.labelcolor": "#0b0b0b",
    "text.color": "#0b0b0b", "xtick.color": INK_SECOND,
    "ytick.color": INK_SECOND,
    "figure.facecolor": SURFACE, "axes.facecolor": SURFACE,
})


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


def load_interiors() -> dict:
    out = {}
    for d in DOMAINS:
        p = TOPO_DIR / f"domain_{d}_topography_2009.npy"
        if p.is_file():
            out[d] = np.load(p)
    if not out:
        raise SystemExit(f"no topography found in {TOPO_DIR}")
    return out


def _wet_fraction(row: np.ndarray) -> float:
    """
    Fraction of a border row at or below the threshold -- EVERY cell counted.

    A no-data mask (domain_<d>_nodata_2009.npy) exists beside the topography and
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

def draw_island(ax, fig, shown, crop_rows, year, placed, label_sections,
                short):
    colour = C_YEAR[year]
    ax.set_facecolor(WATER)
    im = ax.imshow(np.ma.masked_invalid(shown), aspect="auto", origin="lower",
                   extent=[0.5, len(DOMAINS) + 0.5, -CELL_SIZE_M / 2,
                           crop_rows * CELL_SIZE_M - CELL_SIZE_M / 2],
                   cmap=LAND_CMAP, norm=Normalize(LAND_VMIN, LAND_VMAX),
                   interpolation="nearest")
    cax = ax.inset_axes([1.012, 0.0, 0.014, 1.0])
    cb = fig.colorbar(im, cax=cax, extend="max")
    cb.set_label("interior elev (m MHW)", fontsize=8.5)
    cb.ax.tick_params(labelsize=8)
    cb.outline.set_edgecolor(INK_MUTED)

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

    # No marker on the plan views: the crimson segment IS the signal, and a
    # triangle on top of a 1-domain-wide line obscured the thing it pointed at.
    # The state stays recoverable without colour -- each panel carries a count,
    # and panel D names every failing domain on the same x-axis.
    drowned = [d for d, p in placed.items() if p["drowned"]]
    if drowned:
        ax.text(0.997, 0.06,
                f"{len(drowned)} of {len(placed)} domains drown at "
                f"initialisation",
                transform=ax.transAxes, ha="right", va="bottom", fontsize=9,
                color=C_DROWN, zorder=9,
                bbox=dict(fc=SURFACE, ec=C_DROWN, alpha=0.95, pad=3.0))

    for (lo, hi), name in SECTIONS:
        ax.axvline(hi + 0.5, color=SURFACE, lw=1.0, alpha=0.55, zorder=4)
        if label_sections:
            ax.text((lo + hi) / 2, crop_rows * CELL_SIZE_M * 0.97, name,
                    ha="center", va="top", fontsize=8, color=INK_SECOND,
                    zorder=9,
                    bbox=dict(fc=SURFACE, ec="none", alpha=0.8, pad=1.6))

    ax.set_title(f"{year} road, {short}", loc="left", fontsize=11.5,
                 weight="semibold", color=colour, pad=4)
    ax.set_ylabel("m landward of\ninterior row 0", fontsize=9)
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


def build_figure(name: str, spec: dict, interiors, shown, crop_rows,
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
        placed[year] = place_road(interiors, sb)
        floored[year] = read_floored(spec, year)

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

    fig = plt.figure(figsize=(16.5, 14.6))
    gs = fig.add_gridspec(5, 1, height_ratios=[1.25, 1.25, 0.95, 0.72, 0.72],
                          hspace=0.17, left=0.065, right=0.905,
                          top=0.888, bottom=0.050)
    ax84 = fig.add_subplot(gs[0])
    ax04 = fig.add_subplot(gs[1], sharex=ax84)
    ax_sb = fig.add_subplot(gs[2], sharex=ax84)
    ax_mv = fig.add_subplot(gs[3], sharex=ax84)
    ax_w = fig.add_subplot(gs[4], sharex=ax84)

    for ax, year, lab in ((ax84, YEARS[0], True), (ax04, YEARS[1], False)):
        if year in placed:
            draw_island(ax, fig, shown, crop_rows, year, placed[year], lab,
                        spec['short'])

    # --- (C) setback against the island it has to fit inside ----------------
    width_x = sorted(interiors)
    width_y = [interiors[d].shape[0] * CELL_SIZE_M for d in width_x]
    ax_sb.fill_between(width_x, 0, width_y, color=INK_MUTED, alpha=0.16, lw=0,
                       zorder=1, label="island width (interior rows)")
    ax_sb.plot(width_x, width_y, color=INK_SECOND, lw=1.0, zorder=2)
    for year, pl in placed.items():
        xs = sorted(pl)
        ax_sb.plot(xs, [pl[d]["setback_m"] for d in xs], color=C_YEAR[year],
                   lw=2.2 if year == YEARS[0] else 1.6, zorder=6,
                   label=f"{year} setback")
    ax_sb.set_ylabel("m landward of\ninterior row 0", fontsize=9)
    ax_sb.grid(axis="y", color=INK_MUTED, alpha=0.22, lw=0.7)
    ax_sb.set_axisbelow(True)
    ax_sb.set_ylim(0, min(DISPLAY_CROSS_SHORE_M, max(width_y) * 1.05))
    ax_sb.legend(loc="upper left", fontsize=8.5, ncol=3, framealpha=0.92)
    ax_sb.set_title("(C)  The setback against the island it must fit inside — "
                    "where the lines converge, the road is near the back",
                    loc="left", fontsize=10)
    plt.setp(ax_sb.get_xticklabels(), visible=False)

    # --- (D) where the road moved between the two periods -------------------
    # Ported from the retired 3-figures/island_wide/HAT_plot_road_on_b3d_domains
    # .py, which drew this for one method only.
    if len(placed) == 2:
        ya, yb = sorted(placed)
        common = sorted(set(placed[ya]) & set(placed[yb]))
        move = np.array([placed[yb][d]["setback_m"] - placed[ya][d]["setback_m"]
                         for d in common])
        ax_mv.axhline(0, color=INK_SECOND, lw=1.1, zorder=3)
        ax_mv.bar(common, np.where(move >= 0, move, 0.0), width=0.86,
                  color=C_YEAR[yb], linewidth=0, zorder=5)
        ax_mv.bar(common, np.where(move < 0, move, 0.0), width=0.86,
                  color=C_YEAR[ya], linewidth=0, zorder=5)
        ax_mv.text(0.997, 0.9,
                   f"up = further inland by {yb};  down = closer to the dune."
                   f"   median {np.median(move):+.0f} m",
                   transform=ax_mv.transAxes, ha="right", va="top",
                   fontsize=8.5, color=INK_SECOND,
                   bbox=dict(fc=SURFACE, ec=INK_MUTED, alpha=0.92, pad=2.6))
    ax_mv.set_ylabel(f"{max(placed)} − {min(placed)}\nsetback (m)", fontsize=9)
    ax_mv.grid(axis="y", color=INK_MUTED, alpha=0.22, lw=0.7)
    ax_mv.set_axisbelow(True)
    ax_mv.set_title("(D)  Road movement between the two periods, under this "
                    "method", loc="left", fontsize=10)
    plt.setp(ax_mv.get_xticklabels(), visible=False)

    # --- (D) what the bulldozed band actually lands on ----------------------
    # The series is a PERCENTAGE, so the threshold has to be scaled too -- at
    # DROWN_PCT it would sit on 0.2% and read as zero.
    ax_w.axhspan(DROWN_PCT * 100, 104, color=C_DROWN, alpha=0.07, lw=0,
                 zorder=1)
    ax_w.axhline(DROWN_PCT * 100, color=C_DROWN, lw=1.3, ls=(0, (4, 3)),
                 zorder=3, label=f"drown threshold, {DROWN_PCT * 100:.0f}%")
    for year, pl in sorted(placed.items()):
        xs = sorted(pl)
        ax_w.plot(xs, [pl[d]["governing"] * 100 for d in xs],
                  color=C_YEAR[year], lw=1.8, marker="o", ms=3.0, zorder=5,
                  label=f"{year}")
        bad = [d for d in xs if pl[d]["drowned"]]
        if bad:
            ax_w.plot(bad, [pl[d]["governing"] * 100 for d in bad], lw=0,
                      marker="v", ms=9, mfc=C_DROWN, mec=SURFACE, mew=1.2,
                      zorder=7)
    ax_w.set_ylabel("% of bordering cells\nat or below 0 m MHW", fontsize=9)
    ax_w.set_xlabel("Barrier3D / GIS domain   "
                    "(1 = Cape Point / south  →  90 = Rodanthe / north)")
    ax_w.set_xlim(0.5, len(DOMAINS) + 0.5)
    ax_w.set_ylim(-4, 104)
    ax_w.grid(axis="y", color=INK_MUTED, alpha=0.22, lw=0.7)
    ax_w.set_axisbelow(True)
    ax_w.legend(loc="upper left", fontsize=8.5, ncol=3, framealpha=0.92)
    ax_w.set_title("(E)  bulldoze's own drown test — the worse of the two rows "
                   "BORDERING the road (road_start − 1, road_end + 1). Above "
                   "the line, CASCADE stops managing this roadway",
                   loc="left", fontsize=10)

    # --- header -------------------------------------------------------------
    frame_note = (
        "Measured in the SAME frame it is drawn in — interior row 0 — so this "
        "figure and the measurement agree."
        if name == "dunestart" else
        "Note the frame — this method measured against the digitised dune line, "
        "but CASCADE applies the number landward of interior row 0.")
    n_floored = sum(len(v) for v in floored.values())
    floor_note = (f"  {n_floored} domain-year(s) had a NEGATIVE true setback, "
                  f"floored to 0 (road on interior row 0)." if n_floored else "")

    fig.text(0.065, 0.985,
             f"NC-12 placed by the {spec['label']}, on the Barrier3D interiors "
             f"CASCADE initialises with",
             fontsize=14, va="top", weight="semibold")
    fig.text(0.065, 0.962,
             f"{TOPO_RUN_NAME} topography; the SAME interiors back both road panels, "
             "since one topography serves every period.\n"
             "Road drawn where bulldoze puts it: road_start = int(setback / "
             "10 m). Crimson = the roadway WIDTH-DROWNS at initialisation — "
             f">{DROWN_PCT * 100:.0f}% of the cells bordering it sit at or "
             "below 0 m MHW, so CASCADE stops managing it.\n"
             + frame_note + floor_note,
             fontsize=9, color=INK_SECOND, va="top", linespacing=1.5)
    fig.legend(handles=[
        Line2D([], [], color=C_1984, lw=2.6, label="1984 road"),
        Line2D([], [], color=C_2004, lw=2.6, label="2004 road"),
        Line2D([], [], color=C_DROWN, lw=3.4,
               label="drowns at initialisation"),
        Line2D([], [], color=WATER, lw=8, label="off-island / sentinel water"),
    ], loc="upper left", bbox_to_anchor=(0.065, 0.930), ncol=4, fontsize=8.5,
        framealpha=0.0, borderpad=0.4, columnspacing=1.6, handlelength=2.6)

    out_png = spec["root"] / spec["png"]
    out_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_png, dpi=150, bbox_inches="tight", facecolor=SURFACE)
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
    interiors = load_interiors()
    canvas, max_rows = build_canvas(interiors)
    crop_rows = min(max_rows, int(DISPLAY_CROSS_SHORE_M / CELL_SIZE_M))
    shown = canvas[:crop_rows, :]
    print(f"{len(interiors)} interiors from {TOPO_DIR.parent.name} | "
          f"island width {min(a.shape[0] for a in interiors.values()) * 10}"
          f"-{max_rows * 10} m")

    drowned = {}
    for name, spec in METHODS.items():
        got = build_figure(name, spec, interiors, shown, crop_rows, max_rows)
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
