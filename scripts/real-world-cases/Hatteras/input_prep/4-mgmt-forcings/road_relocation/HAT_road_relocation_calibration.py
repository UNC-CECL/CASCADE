"""
NC-12 relocation in the 1996-2009 calibration period: the figures of HAT_road_relocation_distance.py,
redrawn for that window only, as one three-panel figure and as each panel on its own.

    python scripts/real-world-cases/Hatteras/input_prep/4-mgmt-forcings/road_relocation/HAT_road_relocation_calibration.py

Only the relocations the model fires inside the window are drawn (HATTERAS_ROAD_EVENTS, years 1996-2009): the
1999 Buxton move (GIS 9-14). The 1989 Pea Island / S-curves move (GIS 84-87) is before 1996, so in this window the
road there is already on its 2008 line and does not move.

Reads (nothing is re-measured):
    road_relocation/1978_2008/road_relocation_1978_2008.csv   the measured move per domain
    nc12_1978.geojson (the road to 1999 at Buxton), nc12_2008.geojson (from 2000)
    the 2008 NOAA aerial mosaic (hat_lidar_strips.imagery), under the Buxton map; cached in supporting/

Writes to data/hatteras_init/4-mgmt-forcing/road_relocation/1996_2009_calibration/:
    road_relocation_1996_2009.png              (a) island overview, (b) Buxton map, (c) the move per domain
    road_relocation_1996_2009_overview.png     (a) alone
    road_relocation_1996_2009_buxton.png       (b) alone
    road_relocation_1996_2009_by_domain.png    (c) alone
    supporting/ (PDFs, CAPTIONS.md, the aerial cache)

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-11-20
"""

import sys
from pathlib import Path

import geopandas as gpd
import matplotlib

matplotlib.use("Agg")
import matplotlib.patheffects as pe  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from matplotlib.colors import BoundaryNorm, LightSource, ListedColormap  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402
from matplotlib.patches import Patch, Rectangle  # noqa: E402
from shapely.affinity import rotate as shapely_rotate  # noqa: E402
from shapely.ops import unary_union  # noqa: E402

_HAT_DIR = next(_p for _p in Path(__file__).resolve().parents if (_p / "hatteras_paths.py").exists())
if str(_HAT_DIR) not in sys.path:
    sys.path.insert(0, str(_HAT_DIR))

from site_layer import hat_figure_style as S_  # noqa: E402
from site_layer import hat_topo_version as _tv  # noqa: E402
from site_layer.hatteras_site_config import (  # noqa: E402
    HATTERAS_ANNOTATIONS, HATTERAS_ROAD_EVENTS, RelocationEvent, round_to_cell)

# --- CONFIG ------------------------------------------------------------------
PERIOD = (1996, 2009)                    # the calibration window
LINE_BEFORE, LINE_AFTER = 1978, 2008     # the digitised centrelines: the road before and after the move
SITE_NAME = "North of Buxton"
EVENT_NOTE = ("washed out by Hurricane Dennis, September 1999; relocated road in place by 16 March 2000 "
              "(NCDOT 2015)")
# -----------------------------------------------------------------------------

# The look: one line-weight ladder and one type scale for every panel (pt)
LW_ROAD, LW_ROAD_OVERVIEW = 1.1, 0.8     # the subject: heaviest
LW_FRAME, LW_DOMAIN, LW_GRID = 0.5, 0.35, 0.3
FS, FS_AXIS = 7.0, 8.0                   # labels and ticks; axis titles
ROAD_C = "#AE017E"                       # NC-12, the magenta of the lidar figures (hat_lidar_strips.ROAD_C)
OLD_DASH = (0, (5, 2.5))
NOT_MOVED_C, NO_ROAD_C, EDGE_C = "0.90", "white", "0.55"
WATER_C = "#DCE9F2"
HIGHLIGHT_C = "#4292C6"                  # the moved stretch on the overview, and the measured bars
STEP_M = 10                              # the move in whole Barrier3D cells, as forced
STEPS = np.arange(10, 90, STEP_M)

OUT = _tv.road_relocation_dir(LINE_BEFORE, LINE_AFTER).parent / f"{PERIOD[0]}_{PERIOD[1]}_calibration"
MEASURED_CSV = _tv.road_relocation_file(LINE_BEFORE, LINE_AFTER)
STEM = f"road_relocation_{PERIOD[0]}_{PERIOD[1]}"
AERIAL_YEAR, AERIAL_RES_M = 2008, 1.5
AERIAL_CACHE = OUT / "supporting" / f"buxton_aerial_{AERIAL_YEAR}.npz"



def set_style():
    """The house style plus this figure's type and line scale. Called again before drawing: importing the lidar
    module (for the aerial) resets matplotlib's settings."""
    S_.apply_style()
    plt.rcParams.update({"font.size": FS, "axes.labelsize": FS_AXIS, "xtick.labelsize": FS, "ytick.labelsize": FS,
                         "legend.fontsize": FS, "axes.linewidth": LW_FRAME, "xtick.major.width": LW_FRAME,
                         "ytick.major.width": LW_FRAME, "xtick.major.size": 2.5, "ytick.major.size": 2.5,
                         "savefig.dpi": 300, "figure.dpi": 100})


set_style()
INK, INK_MUTED = S_.INK, S_.INK_MUTED
HALO = [pe.withStroke(linewidth=2.2, foreground="white")]
MAP_HALO = [pe.withStroke(linewidth=1.6, foreground="0.12")]     # white type on the photograph (house map style)


# The data ----------------------------------------------------------------------

events = [e for e in HATTERAS_ROAD_EVENTS
          if isinstance(e, RelocationEvent) and PERIOD[0] <= e.year <= PERIOD[1]]
skipped = [e for e in HATTERAS_ROAD_EVENTS
           if isinstance(e, RelocationEvent) and not PERIOD[0] <= e.year <= PERIOD[1]]
assert len(events) == 1, f"expected the 1999 Buxton event only, got {[e.year for e in events]}"
EVENT = events[0]
MOVED = sorted(EVENT.displacement_m)                               # GIS 9-14
SITE = list(range(MOVED[0] - 1, MOVED[-1] + 2))                    # GIS 8-15: the ends, where the lines cross
print(f"In {PERIOD[0]}-{PERIOD[1]}: {EVENT.note}; left out: " + "; ".join(e.note for e in skipped))

measured = pd.read_csv(MEASURED_CSV).set_index("domain")
mean_m = {g: float(measured.loc[g, "mean_signed_landward_m"]) for g in MOVED}      # measured mean, landward +
model_m = {g: round_to_cell(float(EVENT.displacement_m[g])) for g in MOVED}         # what the model is forced with
road_domains = set(measured.index.astype(int))                                     # where NC-12 runs

cmap = ListedColormap(plt.get_cmap("Blues")(np.linspace(0.20, 0.85, len(STEPS))))
norm = BoundaryNorm(np.append(STEPS - STEP_M / 2, STEPS[-1] + STEP_M / 2), cmap.N)


CRS = "EPSG:6347"
from site_layer import hat_observed_rates as _obs  # noqa: E402   the domain file the measurement used
domains = gpd.read_file(_obs.DOMAIN_BOXES)
domains =domains.dropna(subset=["domain_id"]).astype({"domain_id": int}).sort_values("domain_id").to_crs(CRS)
road_before = gpd.read_file(_tv.road_line_file(LINE_BEFORE)).to_crs(CRS)
road_after = gpd.read_file(_tv.road_line_file(LINE_AFTER)).to_crs(CRS)
site = domains[domains.domain_id.isin(SITE)].copy()
site_area = unary_union(site.geometry).buffer(5)
old_road = gpd.GeoDataFrame(geometry=[unary_union(road_before.geometry).intersection(site_area)], crs=CRS)
SB = site.total_bounds                                              # the Buxton map's frame: GIS 8-15, no more


# Drawing helpers -----------------------------------------------------------------

def draw_line(ax, frame, **kw):
    """A line layer drawn part by part (geopandas will not take a dash tuple on multi-part lines)."""
    for geom in frame.geometry:
        for part in getattr(geom, "geoms", [geom]):
            if part.is_empty or part.geom_type != "LineString":
                continue
            ax.plot(*part.xy, solid_capstyle="round", dash_capstyle="butt", **kw)


def scale_bar(ax, length_m, x_frac, y_frac, segments=2, height_frac=0.010):
    """A segmented black-and-white bar, labelled at its ends and middle."""
    x0, x1 = ax.get_xlim()
    y0, y1 = ax.get_ylim()
    xs, ys, h = x0 + x_frac * (x1 - x0), y0 + y_frac * (y1 - y0), height_frac * (y1 - y0)
    seg = length_m / segments
    for k in range(segments):
        ax.add_patch(Rectangle((xs + k * seg, ys), seg, h, facecolor=INK if k % 2 == 0 else "white",
                               edgecolor=INK, lw=0.5, zorder=20))
    unit_km = length_m >= 1000
    for k in range(segments + 1):
        v = k * seg / (1000 if unit_km else 1)
        txt = f"{v:g}" + ((" km" if unit_km else " m") if k == segments else "")
        ax.text(xs + k * seg, ys + 1.6 * h, txt, ha="center", va="bottom", fontsize=FS - 0.5, color=INK,
                zorder=20, path_effects=HALO)


def north_arrow(ax, x, y, length=0.05, direction="up"):
    dx, dy = (0, length) if direction == "up" else (length, 0)
    ax.annotate("", xy=(x + dx, y + dy), xytext=(x, y), xycoords="axes fraction",
                arrowprops=dict(arrowstyle="-|>", color=INK, lw=0.8, mutation_scale=9, shrinkA=0, shrinkB=0),
                zorder=20)
    tx, ty = (x, y + length + 0.012) if direction == "up" else (x + length + 0.008, y)
    ax.text(tx, ty, "N", transform=ax.transAxes, ha="center" if direction == "up" else "left",
            va="bottom" if direction == "up" else "center", fontsize=FS, fontweight="bold", color=INK, zorder=20,
            path_effects=HALO)


def map_frame(ax):
    ax.set_xticks([])
    ax.set_yticks([])
    for sp in ax.spines.values():
        sp.set_visible(True)
        sp.set_linewidth(LW_FRAME)
        sp.set_color(INK)


def letter(ax, i, x=0.0, y=1.01):
    ax.text(x, y, f"({chr(97 + i)})", transform=ax.transAxes, ha="left", va="bottom", fontsize=FS_AXIS,
            fontweight="bold", color=INK)


def road_handles():
    return [Line2D([], [], color=ROAD_C, lw=LW_ROAD, ls=OLD_DASH, label="NC-12 before the move (to 1999)"),
            Line2D([], [], color=ROAD_C, lw=LW_ROAD, label="NC-12 after the move (from 2000)")]


# (a) The island: the one stretch that moved -----------------------------------------

_origin = tuple(unary_union(domains.geometry.tolist()).centroid.coords[0])


def to_strip(frame):
    """Turned 90 degrees clockwise: south left, north right, the ocean below."""
    return frame.set_geometry(frame.geometry.apply(lambda g: shapely_rotate(g, -90, origin=_origin)), crs=frame.crs)


STRIP = to_strip(domains.assign(cls=np.where(domains.domain_id.isin(MOVED), "moved",
                                             np.where(domains.domain_id.isin(road_domains), "road", "none"))))
STRIP_ROAD = to_strip(gpd.GeoDataFrame(geometry=[unary_union(road_after.geometry).intersection(
    unary_union(domains.geometry).buffer(50))], crs=CRS))
STRIP_B = STRIP.total_bounds


def draw_overview(ax):
    mx0, my0, mx1, my1 = STRIP_B
    h = my1 - my0
    for cls, colour in (("none", NO_ROAD_C), ("road", NOT_MOVED_C), ("moved", HIGHLIGHT_C)):
        STRIP[STRIP.cls == cls].plot(ax=ax, facecolor=colour, edgecolor=EDGE_C, lw=LW_DOMAIN, zorder=1)
    draw_line(ax, STRIP_ROAD, color=ROAD_C, lw=LW_ROAD_OVERVIEW, zorder=5)
    ax.set_xlim(mx0 - 0.075 * (mx1 - mx0), mx1 + 0.075 * (mx1 - mx0))
    ax.set_ylim(my0 - 0.42 * h, my1 + 0.42 * h)
    ax.set_aspect("equal")
    ax.set_axis_off()
    # Domain numbers above the strip, horizontal, every tenth and the moved stretch's ends
    for _, row in STRIP.iterrows():
        g = int(row.domain_id)
        if g in (1, *range(20, 91, 10)) or g in (MOVED[0], MOVED[-1]):
            p = row.geometry.representative_point()
            moved = g in MOVED
            ax.text(p.x, my1 + 0.06 * h, str(g), ha="center", va="bottom", fontsize=FS - 0.5,
                    color=INK if moved else INK_MUTED, fontweight="bold" if moved else "normal")
    ax.text(mx1 + 0.012 * (mx1 - mx0), my1 + 0.06 * h, "GIS", ha="left", va="bottom", fontsize=FS - 0.5,
            color=INK_MUTED, style="italic")
    # The moved stretch, named once
    m = STRIP[STRIP.domain_id.isin(MOVED)].total_bounds
    yb = my1 + 0.30 * h
    ax.annotate("", xy=(m[0], yb), xytext=(m[2], yb), annotation_clip=False,
                arrowprops=dict(arrowstyle="|-|", lw=LW_FRAME + 0.3, color=INK, mutation_scale=2.5))
    ax.text(m[2] + 0.012 * (mx1 - mx0), yb, f"{SITE_NAME}: NC-12 moved landward, 1999–2000",
            ha="left", va="center", fontsize=FS, color=INK)
    # One row of place names under the strip; the ends beside it
    for name, (lo, hi) in HATTERAS_ANNOTATIONS.town_spans.items():
        b = STRIP[STRIP.domain_id.between(lo, hi)].total_bounds
        ax.text((b[0] + b[2]) / 2, my0 - 0.06 * h, name, ha="center", va="top", fontsize=FS - 0.5,
                color=INK_MUTED)
    ax.text(mx0 - 0.01 * (mx1 - mx0), (my0 + my1) / 2 + 0.25 * h, "Cape\nPoint", ha="right", va="center",
            fontsize=FS - 0.5, color=INK_MUTED)
    ax.text(mx1 + 0.01 * (mx1 - mx0), (my0 + my1) / 2 - 0.15 * h, "Pea\nIsland", ha="left", va="center",
            fontsize=FS - 0.5, color=INK_MUTED)
    ax.text(0.46 * mx0 + 0.54 * mx1, my0 - 0.06 * h, "Atlantic Ocean", ha="center", va="top", fontsize=FS,
            color=INK_MUTED, style="italic")
    scale_bar(ax, 10000, 0.66, 0.04, segments=2, height_frac=0.025)
    north_arrow(ax, 0.90, 0.07, length=0.03, direction="right")


def overview_handles():
    return [Patch(facecolor=HIGHLIGHT_C, edgecolor=EDGE_C, lw=LW_DOMAIN, label="NC-12 relocated"),
            Patch(facecolor=NOT_MOVED_C, edgecolor=EDGE_C, lw=LW_DOMAIN,
                  label=f"NC-12 did not move, {PERIOD[0]}–{PERIOD[1]}"),
            Patch(facecolor=NO_ROAD_C, edgecolor=EDGE_C, lw=LW_DOMAIN, label="no NC-12"),
            Line2D([], [], color=ROAD_C, lw=LW_ROAD_OVERVIEW, label="NC-12")]


# (b) The Buxton move at true scale, on the 2009 lidar ----------------------------------

def load_aerial():
    """The 2008 NOAA NGS mosaic over the Buxton map's frame (cached), lightened and part-greyed: a backdrop."""
    if AERIAL_CACHE.is_file():
        rgb = np.load(AERIAL_CACHE)["rgb"]
    else:
        from site_layer import hat_lidar_strips as L
        rgb = L.imagery(AERIAL_YEAR, OUT / "supporting").read(*FRAME, AERIAL_RES_M, CRS)
        AERIAL_CACHE.parent.mkdir(parents=True, exist_ok=True)
        np.savez_compressed(AERIAL_CACHE, rgb=rgb)
    img = np.clip(rgb / 255.0, 0, 1)
    grey = img.mean(axis=-1, keepdims=True)
    img = 0.55 * img + 0.45 * grey                                  # half the colour
    img = 0.62 * img + 0.38                                         # lightened toward white
    img[rgb.sum(axis=-1) == 0] = matplotlib.colors.to_rgb(WATER_C)   # beyond the mosaic: the ocean
    return img


_pad = 0.012 * (SB[3] - SB[1])
FRAME = (SB[0] - _pad, SB[2] + _pad, SB[1] - _pad, SB[3] + _pad)   # x0, x1, y0, y1
AERIAL = load_aerial()


def draw_buxton(ax):
    ax.imshow(AERIAL, extent=FRAME, interpolation="bilinear", zorder=0)
    site.boundary.plot(ax=ax, color="white", lw=1.6, alpha=0.9, zorder=3)      # outlines that read on a photograph
    site.boundary.plot(ax=ax, color="0.35", lw=0.45, zorder=4)
    draw_line(ax, old_road, color=ROAD_C, lw=LW_ROAD, ls=OLD_DASH, zorder=11)
    draw_line(ax, road_after, color=ROAD_C, lw=LW_ROAD, zorder=12)
    ax.set_xlim(FRAME[0], FRAME[1])
    ax.set_ylim(FRAME[2], FRAME[3])
    ax.set_aspect("equal")
    map_frame(ax)
    for _, row in site.iterrows():
        g = int(row.domain_id)
        gx0, gy0, gx1, gy1 = row.geometry.bounds
        cy = gy0 + 0.16 * (gy1 - gy0) if g == SITE[0] else (gy0 + gy1) / 2   # GIS 8: below the road's bend
        moved = g in MOVED
        label = f"GIS {g}\n{model_m[g]:.0f} m" if moved else f"GIS {g}\nnot moved"
        ax.text(gx0 + 0.035 * (gx1 - gx0), cy, label, ha="left", va="center", fontsize=FS, linespacing=1.25,
                color="white", fontweight="bold" if moved else "normal", zorder=15, path_effects=MAP_HALO)
    ax.text(0.965, 0.30, "Atlantic Ocean", transform=ax.transAxes, rotation=90, ha="center", va="center",
            fontsize=FS, color="#4A6A85", style="italic", zorder=15)
    scale_bar(ax, 500, 0.56, 0.03, segments=2, height_frac=0.008)
    north_arrow(ax, 0.945, 0.06, length=0.045)


def buxton_handles():
    return [*road_handles(),
            Line2D([], [], color="0.35", lw=0.6, label="GIS domain (label: modelled move)")]


# (c) The move per domain, rows aligned with the map's domains -------------------------

ROW_Y = {int(r.domain_id): r.geometry.centroid.y for _, r in site.iterrows()}
ROW_H = {int(r.domain_id): r.geometry.bounds[3] - r.geometry.bounds[1] for _, r in site.iterrows()}


def draw_bars(ax, share_map=True):
    y = np.array([ROW_Y[g] for g in MOVED])
    hgt = 0.55 * np.array([ROW_H[g] for g in MOVED])
    ax.barh(y, [mean_m[g] for g in MOVED], height=hgt, color=HIGHLIGHT_C, edgecolor="none", zorder=2)
    for yy, hh, g in zip(y, hgt, MOVED):
        ax.plot([model_m[g]] * 2, [yy - 0.62 * hh, yy + 0.62 * hh], color=INK, lw=1.2, solid_capstyle="butt",
                zorder=3)
    for g in (SITE[0], SITE[-1]):
        ax.text(1.5, ROW_Y[g], "not moved", ha="left", va="center", fontsize=FS, color=INK_MUTED,
                style="italic")
    if share_map:
        pad = 0.012 * (SB[3] - SB[1])
        ax.set_ylim(SB[1] - pad, SB[3] + pad)
    ax.set_yticks([ROW_Y[g] for g in SITE], [f"GIS {g}" for g in SITE])
    ax.set_xlim(0, 90)
    ax.set_xticks(range(0, 91, 20))
    ax.set_xlabel("NC-12 moved landward (m)")
    ax.grid(axis="x", color="0.88", lw=LW_GRID, zorder=0)
    ax.axvline(0, color=INK, lw=LW_FRAME, zorder=4)
    for sp in ("top", "right", "left"):
        ax.spines[sp].set_visible(False)
    ax.tick_params(axis="y", length=0)


def bar_handles():
    return [Patch(facecolor=HIGHLIGHT_C, label="measured (domain mean)"),
            Line2D([], [], color=INK, lw=1.2, label="in the model (nearest 10 m)")]


# Figures --------------------------------------------------------------------------------

set_style()
W = S_.FIG_W_DOUBLE
ov_aspect = 1.84 * (STRIP_B[3] - STRIP_B[1]) / (1.15 * (STRIP_B[2] - STRIP_B[0]))   # draw_overview's limits
map_aspect = (SB[3] - SB[1]) / (SB[2] - SB[0])


def save(fig, out, caption):
    S_.save(fig, out, bbox_inches="tight", pad_inches=0.04, close=True)
    S_.record_caption(out, caption)
    print(f"wrote {out}")


CAP_EVENT = (f"The one NC-12 relocation in the {PERIOD[0]}–{PERIOD[1]} calibration period, {SITE_NAME.lower()} "
             f"(GIS {MOVED[0]}–{MOVED[-1]}): the road {EVENT_NOTE}. The 1989 relocation north of Rodanthe "
             "(GIS 84–87) is before the window and not drawn; the road there is already on its 2008 line in 1996.")
CAP_A = ("the island in place, turned 90° clockwise (south left, north right, the Atlantic below; turning keeps "
         "distances, so the scale bar holds), the 90 GIS domains shaded by whether NC-12 moved in the window; "
         "magenta, NC-12 from 2000")
CAP_B = (f"GIS {SITE[0]}–{SITE[-1]} at true scale, north up, on the NOAA NGS ortho-rectified colour mosaic of "
         "26-27 March 2008 (InPort 48695; lightened); magenta, NC-12 before the move (dashed, the road to 1999, digitised "
         f"{LINE_BEFORE}) and after (solid, from 2000, digitised {LINE_AFTER}); domains outlined and labelled by the "
         f"landward move the model makes in 1999, the measured domain mean rounded to the nearest 10 m cell. GIS "
         f"{SITE[0]} and {SITE[-1]}, at the ends of the move where the two lines cross, are not moved")
CAP_C = ("the move per domain, each row level with its domain in the map: bars, the measured landward move of the "
         "centreline, the mean of samples every 5 m; black ticks, the move the model makes (the mean rounded to the "
         f"nearest 10 m cell). Measurements from {MEASURED_CSV.name}")

# The combined figure: (a) across the top; (b) map and (c) bars below, sharing the domains' rows
map_w = 3.15
map_h = map_w * map_aspect
ov_h = (W - 0.1) * ov_aspect
fig = plt.figure(figsize=(W, ov_h + map_h + 1.45))
H = fig.get_size_inches()[1]
top = 0.10
ax_a = fig.add_axes([0.0, (H - top - ov_h) / H, 1.0, ov_h / H])
row_y0 = 0.95
ax_b = fig.add_axes([0.06, row_y0 / H, map_w / W, map_h / H])
ax_c = fig.add_axes([(0.06 * W + map_w + 0.75) / W, row_y0 / H, (W - 0.06 * W - map_w - 1.0) / W, map_h / H])
draw_overview(ax_a)
draw_buxton(ax_b)
draw_bars(ax_c)
letter(ax_a, 0, x=0.0, y=0.92)
letter(ax_b, 1, x=0.0, y=1.012)
letter(ax_c, 2, x=-0.16, y=1.012)
fig.legend(handles=[*overview_handles()[:3], *road_handles(), *bar_handles()], loc="lower center",
           bbox_to_anchor=(0.5, 0.0), ncol=4, frameon=False, handlelength=2.4, columnspacing=1.4)
save(fig, OUT / f"{STEM}.png",
     f"{CAP_EVENT} (a) {CAP_A}. (b) {CAP_B}. (c) {CAP_C}.")

# Each panel on its own
fig, ax = plt.subplots(figsize=(W, W * ov_aspect + 0.45), layout="constrained")
draw_overview(ax)
fig.legend(handles=overview_handles(), loc="outside lower center", ncol=4, frameon=False)
save(fig, OUT / f"{STEM}_overview.png", f"{CAP_EVENT} The {CAP_A}.")

fig, ax = plt.subplots(figsize=(3.6, 3.6 * map_aspect + 0.6), layout="constrained")
draw_buxton(ax)
fig.legend(handles=buxton_handles(), loc="outside lower center", ncol=1, frameon=False)
save(fig, OUT / f"{STEM}_buxton.png", f"{CAP_EVENT} {CAP_B[0].upper()}{CAP_B[1:]}.")

fig, ax = plt.subplots(figsize=(3.4, 3.0), layout="constrained")
draw_bars(ax, share_map=False)
ax.set_ylim(min(ROW_Y.values()) - 0.6 * max(ROW_H.values()), max(ROW_Y.values()) + 0.6 * max(ROW_H.values()))
fig.legend(handles=bar_handles(), loc="outside lower center", ncol=2, frameon=False)
save(fig, OUT / f"{STEM}_by_domain.png",
     f"{CAP_EVENT} {CAP_C[0].upper()}{CAP_C[1:]}; rows south (bottom) to north (top), as in "
     f"{STEM}_buxton.")

# The first draft's separate figures, replaced by the four above
for old in ("alongshore", "domain_map"):
    for f in (OUT / f"{STEM}_{old}.png", OUT / "supporting" / f"{STEM}_{old}.pdf"):
        if f.is_file():
            f.unlink()
            print(f"removed {f.name}")
