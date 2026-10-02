"""
The generic site figures: where Hatteras is, how the 90 domains tile it, one domain as a grid, the forcing.

    python scripts/figure_making/island/study_area_figures.py [--vector] [--only NAME] [--talk]

Drawn from the model's own inputs; --vector uses the island outline instead of
satellite tiles (no network). Writes into output/figures/1-site/, 2-observations/
and 3-model-inputs/ (or talk/). Details: scripts/figure_making/island/README.md.

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-09-30
"""
from __future__ import annotations

import argparse
import math
import typing
import sys
import tempfile
from pathlib import Path

import geopandas as gpd
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.legend_handler import HandlerBase
from matplotlib.lines import Line2D
from matplotlib.patches import Patch, Polygon, Rectangle
from matplotlib.transforms import Affine2D
from shapely.affinity import rotate as shp_rotate, translate as shp_translate
from shapely.geometry import LineString, MultiLineString

REPO = next(p for p in Path(__file__).resolve().parents if (p / "pyproject.toml").exists())
sys.path.insert(0, str(REPO / "scripts"))

from site_layer.hat_figure_style import (  # noqa: E402
    apply_style, C, C_1984, C_1997, C_1984_FILL, C_1997_FILL, INK, INK_MUTED, GRID_C, CELL_M,
    DOMAIN_AXIS_LABEL, FIG_W_DOUBLE, FIG_W_SINGLE,
    figsize, save, record_caption, _title, _letter_inside, _scalebar, _halo, _north_arrow,
    spines_for_image, open_frame, elevation_cmap, ELEV_WATER, town_bands, structures,
    MAP_HALO, MAP_TEXT, MAP_TEXT_DARK, GROIN_MARKER, MAP_LEGEND, spaced_caps, water_label, place_label,
    map_label, scale_bar_km, north_dart, groin_marker, groin_handle, letter_at, top_row, axes_point,
)
from site_layer import hat_topo_version as tv  # noqa: E402
from site_layer.hatteras_site_config import (  # noqa: E402
    HATTERAS_ANNOTATIONS as ANN, HATTERAS_COMMUNITY_ZONES, HATTERAS_DOMAINS,
    HATTERAS_NOURISHMENT_PROJECTS, HATTERAS_ROAD_EVENTS, SCORE_INTERIOR_GIS,
)

# --- CONFIG ------------------------------------------------------------------
CRS = "EPSG:26918"
INIT = REPO / "data" / "hatteras_init"
from site_layer.hat_observed_rates import DOMAIN_BOXES  # noqa: E402
from site_layer import hat_map_layers as _ml  # noqa: E402
OUTLINE = _ml.ISLAND_OUTLINE
from site_layer.hat_figure_style import FIGURES_ROOT as FIG_ROOT, figure_dir  # noqa: E402

# Which numbered output/figures/ folder each figure belongs to
SUBJECT = {
    "study_area": ("site",),
    "site_overview": ("site",),
    "domain_framework": ("site",),
    "domain_framework_vertical": ("site",),
    "domain_metrics": ("site",),
    "domain_grid": ("site",),
    "domain_schematic": ("site",),
    "reach_elevation": ("site",),
    "forcing_timeline_1984": ("inputs", "3-forcing"),
    "forcing_timeline_1996": ("inputs", "3-forcing"),
    "management_footprint": ("inputs", "4-management"),
    "observed_rates": ("observations", "shoreline"),
    "dune_lines": ("observations", "duneline"),
}
TALK = False        # set by talk_mode(); the projector set mirrors the layout


# `output/figures/<numbered subject>/.../<name>.png`, the folder created
def fig_path(name):
    d = figure_dir(*SUBJECT[name])
    if TALK:
        d = FIG_ROOT / "talk" / d.relative_to(FIG_ROOT)
    d.mkdir(parents=True, exist_ok=True)
    return d / f"{name}.png"
TILE_CACHE = Path(tempfile.gettempdir()) / "hat_tile_cache"

FIRST = HATTERAS_DOMAINS.first_gis_id
LAST = FIRST + HATTERAS_DOMAINS.num_real_domains - 1
N_BUFFER = HATTERAS_DOMAINS.num_buffer_domains
EXAMPLE_GIS = 45
TOPO_PRODUCT = "2004-start"

# The vector base map: pale water, ivory land, grey edges
WATER_MAP = "#e9eff4"
LAND = "#ede9df"
LAND_EDGE = "0.55"
ROAD_ON_IMAGERY = "#ffd23f"   # a road-map yellow: legible on sand and asphalt
PROFILE_GROUND = "#cdb26f"    # a cross-shore profile's ground fill: a profile, not a topography map, so not terrain
ROLE = {                       # fills for domain_framework
    "boundary": "0.40",        # the end domains: dark, few, the boundary condition
    "community": "#e6b39a",    # the community zones: the settlement tint of a topographic map
}
HALO = _halo(2.0)
# White type on imagery needs a dark halo, the opposite of the house one
DARK_HALO = [mpl.patheffects.withStroke(linewidth=1.8, foreground="0.15")]
# -----------------------------------------------------------------------------


# Layers

# The domain boxes, island outline and NC-12 lines, in UTM
def load_layers():
    dom = gpd.read_file(DOMAIN_BOXES).to_crs(CRS)
    dom["ID"] = dom["domain_id"].astype(int)
    dom = dom[(dom.ID >= FIRST) & (dom.ID <= LAST)].sort_values("ID").reset_index(drop=True)
    outline = gpd.read_file(OUTLINE).to_crs(CRS)
    roads = {v: gpd.read_file(tv.road_line_file(v)).to_crs(CRS) for v in tv.ROAD_LINE_VINTAGES}
    return dom, outline, roads


# A polygon's exterior as an (n, 2) array
def xy2d(geom):
    return np.asarray(geom.exterior.coords)[:, :2]


# The rotation that lays the reach out left to right, and which way is seaward in it
class Frame:

    def __init__(self, dom, outline):
        c = dom.geometry.centroid
        xy = np.c_[c.x.values, c.y.values]
        p1, p2 = xy[0], xy[-1]
        # The boxes are axis-aligned in UTM, north along the reach, so the frame is a quarter turn
        self.u = np.array([0.0, 1.0]) if (p2 - p1)[1] > 0 else np.array([0.0, -1.0])
        self.theta = math.atan2(self.u[1], self.u[0])
        self.origin = tuple(xy.mean(0))
        self.step = float(np.median(np.abs(np.diff(xy[:, 1]))))     # 500 m, the box's alongshore side
        self.centroids = xy
        # the island's land centroid lies soundward of the reach axis
        land_c = self.pts([outline.geometry.union_all().centroid.coords[0]])[0]
        cen = self.pts(xy)
        self.seaward = np.array([0.0, -1.0]) if land_c[1] > cen[:, 1].mean() else np.array([0.0, 1.0])

    @property
    def deg(self):
        return -math.degrees(self.theta)

    def geoms(self, gs):
        return gs.rotate(self.deg, origin=self.origin)

    def pts(self, xy):
        xy = np.asarray(xy, float).reshape(-1, 2) - self.origin
        c, s = math.cos(-self.theta), math.sin(-self.theta)
        return np.c_[c * xy[:, 0] - s * xy[:, 1], s * xy[:, 0] + c * xy[:, 1]] + self.origin

    # UTM bounding box that covers a rotated-frame window
    def unrotate_bounds(self, x0, x1, y0, y1):
        corners = np.array([[x0, y0], [x1, y0], [x1, y1], [x0, y1]], float) - self.origin
        c, s = math.cos(self.theta), math.sin(self.theta)
        back = np.c_[c * corners[:, 0] - s * corners[:, 1],
                     s * corners[:, 0] + c * corners[:, 1]] + self.origin
        return back[:, 0].min(), back[:, 1].min(), back[:, 0].max(), back[:, 1].max()

    # Unit vector of true north in the rotated frame
    def north(self):
        return np.array([math.sin(self.theta), math.cos(self.theta)])

    def image_transform(self, ax):
        return Affine2D().rotate_deg_around(*self.origin, self.deg) + ax.transData

    # Rotated-frame point on the reach axis at a (fractional) GIS id
    def along(self, gis):
        cen = self.pts(self.centroids)
        i = np.asarray(gis, float) - FIRST
        idx = np.arange(len(cen))
        return np.c_[np.interp(i, idx, cen[:, 0]), np.interp(i, idx, cen[:, 1])]

    # (x0, x1, y0, y1) of the reach with asymmetric padding, rotated frame
    def window(self, dom, pad_along_km, pad_sea_km, pad_sound_km):
        x0, y0, x1, y1 = self.geoms(dom.geometry).total_bounds
        x0 -= pad_along_km * 1000
        x1 += pad_along_km * 1000
        if self.seaward[1] < 0:
            y0 -= pad_sea_km * 1000
            y1 += pad_sound_km * 1000
        else:
            y0 -= pad_sound_km * 1000
            y1 += pad_sea_km * 1000
        return x0, x1, y0, y1


# A north arrow pointing to true north in a rotated, equal-aspect frame
def north_arrow_rotated(ax, frame, x=0.965, y=0.10, length=0.16):
    n = frame.north()
    x0, x1 = ax.get_xlim()
    y0, y1 = ax.get_ylim()
    L = length * (y1 - y0)
    dx, dy = n[0] * L / (x1 - x0), n[1] * L / (y1 - y0)
    x = x - max(dx, 0) - 0.02
    y = y - min(dy, 0) + 0.0
    tip = (x + dx, y + dy)
    ax.annotate("", xy=tip, xytext=(x, y), xycoords="axes fraction", textcoords="axes fraction",
                arrowprops=dict(arrowstyle="-|>", color=INK, lw=1.0, shrinkA=0, shrinkB=0,
                                mutation_scale=11), zorder=20)
    lab = (tip[0] + n[0] * 0.016 * (y1 - y0) / (x1 - x0), tip[1] + n[1] * 0.016)
    ax.text(*lab, "N", transform=ax.transAxes, ha="center", va="center", fontsize=8.5,
            fontweight="bold", color=INK, zorder=20,
            bbox=dict(facecolor="white", alpha=0.8, edgecolor="none", boxstyle="square,pad=0.1"))


# Imagery

# Tile mosaic for a UTM window
def tiles(bounds_utm, zoom, source, t_crs=CRS):
    import contextily as cx
    from rasterio.warp import transform_bounds
    cx.set_cache_dir(str(TILE_CACHE))
    w, s, e, n = transform_bounds(CRS, "EPSG:3857", *bounds_utm)
    img, ext = cx.bounds2img(w, s, e, n, zoom=zoom, source=source, ll=False)
    if t_crs == "EPSG:3857":
        return img, ext
    return cx.warp_tiles(img, ext, t_crs=t_crs)


# The satellite tile provider
def imagery_source():
    import contextily as cx
    return cx.providers.Esri.WorldImagery


# The plain canvas tile provider
def canvas_source():
    import contextily as cx
    return cx.providers.Esri.WorldGrayCanvas


# The panel letter inside a frame but clear of it
def letter_corner(ax, i, x=0.06, y=0.94, ha="left"):
    ax.text(x, y, f"({chr(ord('a') + i)})", transform=ax.transAxes, ha=ha, va="top",
            fontsize=10, fontweight="bold", color=INK, zorder=20,
            bbox=dict(facecolor="white", alpha=0.85, edgecolor="none", boxstyle="square,pad=0.15"))


# The tile attribution once for the whole figure, in the bottom margin
def credit_figure(fig, text, x=0.995, y=0.004):
    fig.text(x, y, text, ha="right", va="bottom", fontsize=6.5, color=INK_MUTED, zorder=25)


# Tile attribution: the one text a licence puts on the canvas
def credit(ax, text, loc="lower right"):
    x, ha = (0.988, "right") if loc.endswith("right") else (0.012, "left")
    y, va = (0.018, "bottom") if loc.startswith("lower") else (0.982, "top")
    ax.text(x, y, text, transform=ax.transAxes, ha=ha, va=va, fontsize=8, color=INK_MUTED,
            zorder=25, bbox=dict(facecolor="white", alpha=0.7, edgecolor="none", boxstyle="square,pad=0.15"))


# Shared drawing

# Set an axes to a map window: equal aspect, no ticks
def map_axes(ax, window):
    x0, x1, y0, y1 = window
    ax.set_xlim(x0, x1)
    ax.set_ylim(y0, y1)
    ax.set_aspect("equal")
    ax.set_xticks([])
    ax.set_yticks([])
    spines_for_image(ax)


# What draw_reach drew that a caller needs back
class ReachArt(typing.NamedTuple):

    road_c: str
    village_out_m: float


# The rotated reach
def draw_reach(ax, frame, dom, outline, road, vector, window, zoom=12,
               label_every=10, label_villages=True, water_labels=True, box_lw=0.35,
               arrow_xy=(0.965, 0.20), ocean_x=0.55, sound_x=0.20, piers=True, scalebar=True,
               numbers=True, arrow=True, show_road=True, water=None, ends=True, groin=True):
    map_axes(ax, window)
    x0, x1, y0, y1 = window
    if vector:
        ax.set_facecolor(WATER_MAP)
        frame.geoms(outline.geometry).plot(ax=ax, facecolor=LAND, edgecolor=LAND_EDGE, lw=0.3, zorder=1)
        box_c, road_c, water_c = INK, C["ROAD"], INK_MUTED
    else:
        img, (l, r, b, t) = tiles(frame.unrotate_bounds(x0, x1, y0, y1), zoom, imagery_source())
        im = ax.imshow(img, extent=(l, r, b, t), transform=frame.image_transform(ax), zorder=0,
                       interpolation="bilinear")
        if water is not None:
            # one flat water colour; the imagery shows only inside the island outline
            ax.set_facecolor(water)
            im.set_clip_path(outline_patch(ax, frame, outline))
        box_c, road_c, water_c = "white", ROAD_ON_IMAGERY, "white"
        credit(ax, "Imagery: Esri World Imagery")

    frame.geoms(dom.geometry).plot(ax=ax, facecolor="none", edgecolor=box_c, lw=box_lw, zorder=4)
    if show_road:
        frame.geoms(road.geometry).plot(ax=ax, color=road_c, lw=0.9, zorder=5)

    sea = frame.seaward
    cen = frame.pts(frame.centroids)
    for i, g in enumerate(dom.ID.values):
        if numbers and (g == FIRST or g % label_every == 0):
            p = cen[i] + sea * 1500
            ax.text(p[0], p[1], str(g), ha="center", va="center", fontsize=8, color=INK,
                    zorder=8, path_effects=HALO)
    village_out_m = 0.0
    if label_villages:
        village_out_m = village_row(ax, frame, dom, vector)
    for name, pos in (ANN.groins.items() if groin else ()):
        p = frame.along(pos)[0] + sea * 1150          # in the water, just off the box edge
        ax.plot(p[0], p[1], marker="|", ms=7, mew=1.4, color=INK, zorder=9, ls="none")
    for name, (pos, _) in (ANN.piers.items() if piers else ()):
        p = frame.along(pos)[0]
        ax.plot(p[0], p[1], marker="o", ms=3.2, mfc="white", mec=INK, mew=0.8, zorder=9, ls="none")
    if ends:
        p = cen[0] + np.array([-700, 0]) + sea * 1500
        ax.text(p[0], p[1], "Cape\nPoint", ha="right", va="center", fontsize=8, color=INK, zorder=8,
                path_effects=HALO)
        p = cen[-1] + np.array([700, 0]) - sea * 1500
        ax.text(p[0], p[1], "Pea\nIsland", ha="left", va="center", fontsize=8, color=INK, zorder=8,
                path_effects=HALO)
    if water_labels:
        ocean_y, sound_y = (0.05, 0.95) if sea[1] < 0 else (0.95, 0.05)
        ax.text(ocean_x, ocean_y, "ATLANTIC OCEAN", transform=ax.transAxes, ha="center", va="center",
                fontsize=8, color=water_c, alpha=0.9, zorder=8)
        ax.text(sound_x, sound_y, "PAMLICO SOUND", transform=ax.transAxes, ha="center", va="center",
                fontsize=8, color=water_c, alpha=0.9, zorder=8)
    if scalebar:
        _scalebar(ax, 10_000, show_cells=False)
    if arrow:
        north_arrow_rotated(ax, frame, x=arrow_xy[0], y=arrow_xy[1])
    return ReachArt(road_c=road_c, village_out_m=village_out_m)


# The island outline as a clip path in the rotated frame
def outline_patch(ax, frame, outline):
    from matplotlib.patches import PathPatch
    from matplotlib.path import Path as MPath
    paths = []
    for g in frame.geoms(outline.geometry):
        for poly in getattr(g, "geoms", [g]):
            for ring in [poly.exterior, *poly.interiors]:
                paths.append(MPath(np.asarray(ring.coords)[:, :2], closed=True))
    return PathPatch(MPath.make_compound_path(*paths), transform=ax.transData, facecolor="none",
                     edgecolor="none")


# The five village names in one row along the sound side, each with a leader down to its domains
def village_row(ax, frame, dom, vector, clear_m=1300.0, leader_gap_m=700.0,
                pad_m=600.0, text_kw=None):
    sea = frame.seaward
    rdom = frame.geoms(dom.geometry)
    x0, y0, x1, y1 = rdom.total_bounds
    sound_edge_all = y1 if sea[1] < 0 else y0
    y_row = sound_edge_all - sea[1] * clear_m
    names = {"Buxton": ANN.town_spans["Buxton"], "Avon": ANN.town_spans["Avon"],
             "Salvo": (ANN.village_lines["Salvo"],) * 2, "Waves": (ANN.village_lines["Waves"],) * 2,
             "Rodanthe": (ANN.village_lines["Rodanthe"],) * 2}
    lc = "white" if not vector else INK_MUTED
    # Salvo, Waves and Rodanthe are too close to label separately: one row with leaders (README)
    sized = measure_m(ax, list(names), fontsize=8, fontstyle="italic")
    placed = []
    for (name, (lo, hi)), (hw, hh) in zip(names.items(), sized):
        boxes = rdom[(dom.ID >= lo) & (dom.ID <= hi)]
        bx0, by0, bx1, by1 = boxes.total_bounds
        placed.append(dict(x=float(frame.along((lo + hi) / 2)[0][0]), hw=hw,
                           name=name, anchor_x=(bx0 + bx1) / 2,
                           edge=by1 if sea[1] < 0 else by0))
    separate_x(ax, placed, pad_m)
    out = -sea[1]                      # +1/-1 in y, away from the reach
    for d in placed:
        kw = {"color": INK, "path_effects": HALO, **(text_kw or {})}
        ax.text(d["x"], y_row, d["name"], ha="center", va="center", fontsize=8,
                zorder=8, fontstyle="italic", **kw)
        if abs(y_row - d["edge"]) > leader_gap_m:
            ax.plot([d["x"], d["anchor_x"]],
                    [y_row - out * 330, d["edge"] - out * 80],
                    color=lc, lw=0.5, zorder=7)
    # How far out the names actually reach, so whatever is placed beyond them can be told to stay (README)
    hh = max(h for _, h in sized)
    return abs((y_row + out * hh) - sound_edge_all)


# UTM 18N coordinates on a quarter-turned map
def utm_frame(ax, frame, every_x_km=10, every_y_km=5, fontsize=8):
    ox, oy = frame.origin
    x0, x1 = ax.get_xlim()
    y0, y1 = ax.get_ylim()
    # x_rot = ox ± (northing - oy); y_rot = oy ∓ (easting - ox), for theta = ±90 deg
    sgn = 1.0 if frame.theta > 0 else -1.0
    n_lo = min(oy + sgn * (x0 - ox), oy + sgn * (x1 - ox))
    n_hi = max(oy + sgn * (x0 - ox), oy + sgn * (x1 - ox))
    e_lo = min(ox - sgn * (y0 - oy), ox - sgn * (y1 - oy))
    e_hi = max(ox - sgn * (y0 - oy), ox - sgn * (y1 - oy))
    step_x, step_y = every_x_km * 1000, every_y_km * 1000
    norths = np.arange(math.ceil(n_lo / step_x) * step_x, n_hi, step_x)
    easts = np.arange(math.ceil(e_lo / step_y) * step_y, e_hi, step_y)
    ax.set_xticks([ox + sgn * (n - oy) for n in norths])
    ax.set_xticklabels([f"{n / 1000:.0f}" for n in norths], fontsize=fontsize)
    ax.set_yticks([oy - sgn * (e - ox) for e in easts])
    ax.set_yticklabels([f"{e / 1000:.0f}" for e in easts], fontsize=fontsize)
    ax.tick_params(length=2.5, width=0.5, pad=2, colors=INK, direction="out")
    ax.set_xlabel("UTM 18N northing (km)", fontsize=9, labelpad=2)
    ax.set_ylabel("easting (km)", fontsize=9, labelpad=2)


# The framework legend, one column inside the map over open water
def framework_legend(ax, handles, corner="lower right", anchor=None):
    anchor = anchor or {"lower right": (0.985, 0.025), "lower left": (0.015, 0.025),
                        "upper right": (0.985, 0.975), "upper left": (0.015, 0.975)}[corner]
    ax.legend(handles=handles, loc=corner, bbox_to_anchor=anchor, ncol=1, **MAP_LEGEND)


# SHORT labels
def framework_handles():
    return [
        Patch(facecolor="white", edgecolor=INK, lw=0.6,
              label=f"interior, GIS {SCORE_INTERIOR_GIS[0]}–{SCORE_INTERIOR_GIS[1]}"),
        Patch(facecolor=ROLE["boundary"], edgecolor=INK, lw=0.6,
              label=f"ends, GIS {FIRST} and {LAST}"),
        Patch(facecolor=ROLE["community"], edgecolor=INK, lw=0.6, label="community zone"),
        Line2D([], [], color=C["ROAD"], lw=1.0, label="NC-12, 2008"),
        groin_handle(),
    ]


NE_STATES = _ml.NE_STATES
LOCATOR_C = C_1984          # the reach on the locator: the one strong colour on the inset
INSET_LON = (-80.6, -74.6)
INSET_LAT = (33.4, 37.3)


# The south-eastern US coast with the reach marked
def regional_inset(ax, outline, vector=True, lat_side="right", tick_labels=True,
                   locator_c=LOCATOR_C, locator_lw=1.6):
    states = gpd.read_file(NE_STATES)
    lon0, lon1 = INSET_LON
    lat0, lat1 = INSET_LAT
    ax.set_xlim(lon0, lon1)
    ax.set_ylim(lat0, lat1)
    ax.set_aspect(1.0 / math.cos(math.radians((lat0 + lat1) / 2)))
    ax.set_facecolor(WATER_MAP)
    states.plot(ax=ax, facecolor=LAND, edgecolor="white", lw=0.5, zorder=1)      # state borders in white
    states.dissolve().plot(ax=ax, facecolor="none", edgecolor=LAND_EDGE, lw=0.4, zorder=2)   # the coast
    o = outline.to_crs("EPSG:4326")
    o.plot(ax=ax, facecolor=locator_c, edgecolor=locator_c, lw=locator_lw, zorder=4)
    bx0, by0, bx1, by1 = o.total_bounds
    pad = 0.18
    ax.add_patch(Rectangle((bx0 - pad, by0 - pad), bx1 - bx0 + 2 * pad, by1 - by0 + 2 * pad,
                           facecolor="none", edgecolor=INK, lw=0.7, zorder=5))
    # Locator names sit inside one graticule cell each
    ax.text(bx0 - pad - 0.42, (by0 + by1) / 2, "Hatteras\nIsland", ha="right", va="center", fontsize=8,
            color=INK, zorder=6, path_effects=HALO)
    # Only the two names a reader needs
    w_in = ax.get_position().width * ax.figure.get_size_inches()[0]
    if w_in >= 1.05:
        ax.text(-79.0, 35.2, "North\nCarolina", ha="center", va="center", fontsize=8, color=INK_MUTED,
                zorder=6, path_effects=HALO)
    # At 8 pt a small locator holds two names, not three
    if w_in >= 1.45:
        ax.text(-75.45, 34.35, "Atlantic\nOcean", ha="center", va="center", fontsize=8,
                color=INK_MUTED, fontstyle="italic", zorder=6, path_effects=HALO)
    # the graticule: 2 degree lines, labelled on the frame
    lons = np.arange(-80, -74, 2)
    lats = np.arange(34, 38, 2)
    for x in lons:
        ax.axvline(x, color="white", lw=0.5, zorder=3)
    for y in lats:
        ax.axhline(y, color="white", lw=0.5, zorder=3)
    ax.set_xlim(lon0, lon1)
    ax.set_ylim(lat0, lat1)
    if not tick_labels:
        ax.set_xticks([])
        ax.set_yticks([])
        spines_for_image(ax)
        return
    ax.set_xticks(lons)
    ax.set_xticklabels([f"{abs(x)}°W" for x in lons], fontsize=8)
    ax.set_yticks(lats)
    ax.set_yticklabels([f"{y}°N" for y in lats], fontsize=8)
    ax.tick_params(length=2, width=0.5, pad=1.5, colors=INK)
    # The latitude labels go on whichever side has open water beside them
    if lat_side == "left":
        ax.yaxis.tick_left()
    else:
        ax.yaxis.tick_right()
    ax.xaxis.tick_top()
    spines_for_image(ax)


# A double-column figure whose height follows the window's aspect, so an equal-aspect map fills it
def reach_figure(window, panel_frac=(0.004, 0.006, 0.992, 0.988)):
    x0, x1, y0, y1 = window
    aspect = (x1 - x0) / (y1 - y0)
    fx, fy, fw, fh = panel_frac
    height = FIG_W_DOUBLE * fw / aspect / fh
    fig = plt.figure(figsize=figsize("double", height=height))
    return fig, fig.add_axes([fx, fy, fw, fh])


# Figure 1: study area

# Village names just above their own domains on the sound side, shifted off any land, short leaders
def study_area_villages(ax, frame, dom, land, clear_m=1500.0, pad_m=500.0, row=("Salvo", "Waves", "Rodanthe"),
                        shift_m=None, text_kw=None, leader_c="white"):
    from shapely.geometry import box as shp_box
    sea = frame.seaward
    out = -sea[1]                      # +1/-1 in y, toward the sound
    rdom = frame.geoms(dom.geometry)
    names = {"Buxton": ANN.town_spans["Buxton"], "Avon": ANN.town_spans["Avon"],
             "Salvo": (ANN.village_lines["Salvo"],) * 2, "Waves": (ANN.village_lines["Waves"],) * 2,
             "Rodanthe": (ANN.village_lines["Rodanthe"],) * 2}
    sized = measure_m(ax, list(names), fontsize=8, fontstyle="italic")
    placed = []
    for (name, (lo, hi)), (hw, hh) in zip(names.items(), sized):
        bx0, by0, bx1, by1 = rdom[(dom.ID >= lo) & (dom.ID <= hi)].total_bounds
        # the sound-side edge of the boxes around the village, so the name clears its neighbours too
        near = rdom[(dom.ID >= lo - 3) & (dom.ID <= hi + 3)].total_bounds
        edge = near[3] if out > 0 else near[1]
        placed.append(dict(x=(bx0 + bx1) / 2, hw=hw, hh=hh, name=name, anchor_x=(bx0 + bx1) / 2,
                           anchor_y=by1 if out > 0 else by0, y=edge + out * clear_m))
    # the three close villages share one row
    y_row = max(d["y"] * out for d in placed if d["name"] in row) * out
    for d in placed:
        if d["name"] in row:
            d["y"] = y_row
    # the Cape landmass around Buxton is wider than the island outline: a set shift off it
    for d in placed:
        d["x"] += (shift_m or {}).get(d["name"], 0.0)
    separate_x(ax, placed, pad_m)
    # no part of a name on land: slide it toward open water (to the right, along the reach)
    for d in placed:
        for _ in range(60):
            b = shp_box(d["x"] - d["hw"] - 150, d["y"] - d["hh"], d["x"] + d["hw"] + 150, d["y"] + d["hh"])
            if not b.intersects(land):
                break
            d["x"] += 250.0
    for d in placed:
        ax.text(d["x"], d["y"], d["name"], ha="center", va="center", fontstyle="italic",
                **(text_kw or MAP_TEXT))
        lx = min(max(d["anchor_x"], d["x"] - d["hw"]), d["x"] + d["hw"])
        ax.plot([lx, d["anchor_x"]], [d["y"] - out * (d["hh"] + 120), d["anchor_y"] + out * 60],
                color=leader_c, lw=0.5, zorder=7)


# The study-area labels: one rule, white type with a thin dark halo
# (each part can be switched off, and the dark rule passed for a pale map)
def study_area_labels(ax, frame, dom, vector, outline, text_kw=None, leader_c="white",
                      shift_m=None, numbers=True, villages=True, ends=True, water=True,
                      ocean_x=0.55, sound_x=0.20, groin=True, groin_label=True):
    tk = text_kw or MAP_TEXT
    land = frame.geoms(outline.geometry).union_all()
    sea = frame.seaward
    cen = frame.pts(frame.centroids)
    for i, g in enumerate(dom.ID.values):
        if numbers and (g == FIRST or g % 10 == 0):
            p = cen[i] + sea * 1500
            ax.text(p[0], p[1], str(g), ha="center", va="center", **tk)
    if villages:
        study_area_villages(ax, frame, dom, land,
                            shift_m={"Buxton": 4000.0} if shift_m is None else shift_m,
                            text_kw=tk, leader_c=leader_c)
    if ends:
        p = cen[0] + np.array([-700, 0]) + sea * 1500
        ax.text(p[0], p[1], "Cape\nPoint", ha="right", va="center", fontstyle="italic", **tk)
        p = cen[-1] + np.array([700, 0]) - sea * 1500
        ax.text(p[0], p[1], "Pea\nIsland", ha="left", va="center", fontstyle="italic", **tk)
    if water:
        ocean_y, sound_y = (0.05, 0.95) if sea[1] < 0 else (0.95, 0.05)
        ax.text(ocean_x, ocean_y, spaced_caps("Atlantic Ocean"), transform=ax.transAxes, ha="center",
                va="center", **tk)
        ax.text(sound_x, sound_y, spaced_caps("Pamlico Sound"), transform=ax.transAxes, ha="center",
                va="center", **tk)
    # The groin field: a white-edged bar off the beach, its name out in the ocean on a leader
    for name, pos in (ANN.groins.items() if groin else ()):
        p = frame.along(pos)[0] + sea * 1150
        groin_marker(ax, p[0], p[1])
        if not groin_label:
            continue
        end = p + sea * 2700
        ax.text(end[0], end[1] + sea[1] * 250, "Buxton groins", ha="center",
                va="top" if sea[1] < 0 else "bottom", fontstyle="italic", **tk)
        ax.plot([p[0], end[0]], [p[1] + sea[1] * 450, end[1]], color=leader_c, lw=0.5, zorder=7)

# The reach on imagery with domains, NC-12, villages and the regional inset
def fig_study_area(dom, outline, roads, frame, vector):
    window = frame.window(dom, pad_along_km=4.0, pad_sea_km=3.5, pad_sound_km=11.5)
    fig, ax = reach_figure(window)
    fw, fh = fig.get_size_inches()
    road_c = draw_reach(ax, frame, dom, outline, roads[2008], vector, window, piers=False,
                        scalebar=False, arrow=False, numbers=False, label_villages=False,
                        water_labels=False, ends=False, groin=False).road_c
    study_area_labels(ax, frame, dom, vector, outline)
    scale_bar_km(ax)
    # the regional inset in the upper right, over the open sound; square in inches
    ih, edge = 0.58, 0.018
    iw = ih * fh / fw
    ax_in = fig.add_axes([0.996 - edge * fh / fw - iw, 0.994 - edge - ih, iw, ih])
    regional_inset(ax_in, outline, vector, tick_labels=False, locator_c=C["LOCATOR"], locator_lw=2.4)
    letter_corner(ax_in, 1)
    handles = [
        Patch(facecolor="none", edgecolor=INK, lw=0.6, label="model domain, 500 m alongshore"),
        Line2D([], [], color=road_c, lw=1.2, label="NC-12"),
        groin_handle(),
    ]
    # legend left of the inset, (a) and the north arrow on the same top line
    top_row(ax, ax_in, handles, letter=0, north=frame.north())
    out = save(fig, fig_path("study_area"), vector=False, dpi=300)
    record_caption(out[0],
        f"Study area. (a) The modelled reach, rotated so it runs south to north from left to right: "
        f"{LAST - FIRST + 1} Barrier3D domains (GIS 1 at Cape Point, GIS 90 at the southern end of Pea "
        "Island), each a 500 m alongshore by 2000 m cross-shore box in UTM 18N, numbered every tenth; the "
        "NC-12 centreline as digitised on 2008 imagery; the villages of Buxton, Avon, Salvo, Waves and "
        "Rodanthe; the Buxton groin field (red, white-edged bar drawn just offshore, labelled Buxton groins). The map is turned "
        "a quarter turn, north to the right, so the UTM-aligned boxes are level and the reach steps down the "
        "page where the coast bends. (b) Hatteras Island (teal) on the North Carolina coast; the box marks the reach; "
        "the white lines are the 2° graticule (80°, 78° and 76° W; 34° and 36° N). "
        + ("Land is the island outline shapefile. " if vector else "Basemap: Esri World Imagery. ")
        + "Inset: Natural Earth 10 m coastline and state boundaries.")
    plt.close(fig)
    return out[0]


# Figure 2: domain framework

# Per domain, from the extraction arrays (m NAVD88, -10 water)
def domain_metrics(product=TOPO_PRODUCT):
    d, _ = tv.npy_dirs(product)
    gis, width, crest = [], [], []
    for g in range(FIRST, LAST + 1):
        p = d / f"domain_{g}.npy"
        if not p.exists():
            continue
        a = np.load(p)
        land = a > 0.0
        gis.append(g)
        width.append(np.median(land.sum(1)) * CELL_M)
        crest.append(np.nanmedian(np.nanmax(np.where(land, a, np.nan), axis=1)))
    return np.array(gis), np.array(width), np.array(crest)


# The end domains and the community zones filled on the reach map
def draw_role_fills(ax, frame, dom):
    rdom = frame.geoms(dom.geometry)
    zones = set()
    for lo, hi in HATTERAS_COMMUNITY_ZONES:
        zones.update(range(lo, hi + 1))
    for g, geom in zip(dom.ID.values, rdom.values):
        if g in (FIRST, LAST):
            fc = ROLE["boundary"]
        elif g in zones:
            fc = ROLE["community"]
        else:
            continue
        ax.add_patch(mpl.patches.Polygon(xy2d(geom), closed=True, facecolor=fc, edgecolor="none", zorder=3))


# The buffer domains stepped out due north and south of the reach
def draw_buffers(ax, frame, dom, label=True):
    step = frame.u * frame.step
    sea = frame.seaward
    cen = frame.pts(frame.centroids)
    for end, sign in ((0, -1), (len(dom) - 1, +1)):
        base = dom.geometry.iloc[end]
        for k in range(1, N_BUFFER + 1):
            g = shp_translate(base, xoff=sign * k * step[0], yoff=sign * k * step[1])
            g = shp_rotate(g, frame.deg, origin=frame.origin)
            ax.add_patch(mpl.patches.Polygon(xy2d(g), closed=True, facecolor="none", edgecolor=INK_MUTED,
                                             lw=0.35, ls=(0, (2, 1.5)), zorder=4))
        if label:
            p = cen[end] + np.array([sign * (N_BUFFER / 2 + 0.5) * frame.step, 0]) + sea * 3300
            ax.text(p[0], p[1], f"{N_BUFFER} buffer domains", ha="center", va="center", fontsize=8,
                    color=INK_MUTED, zorder=8)


# The domains by role over the island outline, with the inset
def fig_domain_framework(dom, outline, roads, frame, vector):
    window = frame.window(dom, pad_along_km=3.5, pad_sea_km=7.5, pad_sound_km=13.0)
    fig, ax = reach_figure(window, panel_frac=(0.048, 0.10, 0.94, 0.885))
    fw, fh = fig.get_size_inches()
    draw_reach(ax, frame, dom, outline, roads[2008], vector=True, window=window,
               label_villages=False, water_labels=False, numbers=False, ends=False, groin=False,
               piers=False, scalebar=False, arrow=False)
    draw_role_fills(ax, frame, dom)
    utm_frame(ax, frame)
    study_area_labels(ax, frame, dom, True, outline, text_kw=MAP_TEXT_DARK, leader_c=INK_MUTED,
                      shift_m={}, ocean_x=0.60, sound_x=0.25)
    # The regional inset in the upper right corner, the legend and arrow beside it
    ih, iw = 1.3 / fh, 1.3 / fw
    pos = ax.get_position()
    ax_in = fig.add_axes([pos.x1 - 0.06 / fw - iw, pos.y1 - 0.06 / fh - ih, iw, ih])
    regional_inset(ax_in, outline, vector, tick_labels=False, locator_c=C["LOCATOR"], locator_lw=2.4)
    letter_corner(ax_in, 1)
    top_row(ax, ax_in, framework_handles(), letter=0, north=frame.north(), text_kw=MAP_TEXT_DARK)
    out = save(fig, fig_path("domain_framework"), vector=True, dpi=300)
    record_caption(out[0],
        f"The domain framework. (a) The {LAST - FIRST + 1} Barrier3D domain boxes over the island outline, "
        "south to north from left to right (a quarter turn, north to the right, so the UTM-aligned boxes are "
        "level; UTM 18N northing runs along the frame and easting up it), with the roles the hindcast "
        f"assigns them: the interior domains scored against the CoastSat record (GIS "
        f"{SCORE_INTERIOR_GIS[0]}–{SCORE_INTERIOR_GIS[1]}), the end domains that carry the alongshore "
        f"boundary condition (GIS {FIRST}, {LAST}), and the community zones (Buxton 7–8, Avon 21–31, "
        "Salvo–Waves–Rodanthe 68–83) where NC-12 is a maintained street network rather than a relocatable "
        f"road. (The {N_BUFFER} buffer domains beyond each end, which pad the alongshore transport solve, "
        "are not drawn.) The Buxton groin field is the red, white-edged bar drawn offshore. (b) Hatteras "
        "Island (teal) on the North Carolina coast; the box marks the reach; the white lines are the 2° "
        "graticule (80°, 78° and 76° W; 34° and 36° N); Natural Earth 10 m coastline and state boundaries.")
    plt.close(fig)
    return out[0]


# The domain framework NORTH UP, as a portrait figure
def fig_domain_framework_vertical(dom, outline, roads, vector, panel_in=6.5):
    boxes = dom.geometry
    bx0, by0, bx1, by1 = boxes.total_bounds
    pad_n = 2600.0          # room for the end names clear of the frame
    n0, n1 = by0 - pad_n, by1 + pad_n

    fw = FIG_W_SINGLE
    left, right, top, bottom = 0.52, 0.10, 0.30, 0.48      # inches
    pw = fw - left - right
    win_h = n1 - n0
    win_w = win_h * pw / panel_in                          # equal scale in both directions
    # Names on the ocean side, domain numbers on the sound side (README)
    e0 = bx0 - 8000.0      # 8 pt names need a kilometre more of ocean
    e1 = e0 + win_w
    fh = top + panel_in + bottom

    fig = plt.figure(figsize=figsize(fw, height=fh))
    ax = fig.add_axes([left / fw, bottom / fh, pw / fw, panel_in / fh])
    ax.set_xlim(e0, e1)
    ax.set_ylim(n0, n1)
    ax.set_aspect("equal")
    ax.set_facecolor(WATER_MAP)
    spines_for_image(ax)

    # Land, clipped to the window (the outline file also holds Ocracoke and the mainland, 20 km west (README)
    from shapely.geometry import box as _box
    win = _box(e0, n0, e1, n1)
    land = outline.geometry.intersection(win)
    gpd.GeoSeries(land[~land.is_empty], crs=outline.crs).plot(
        ax=ax, facecolor=LAND, edgecolor=LAND_EDGE, lw=0.3, zorder=1)

    # the domains, by role
    zones = set()
    for lo, hi in HATTERAS_COMMUNITY_ZONES:
        zones.update(range(lo, hi + 1))
    boxes.plot(ax=ax, facecolor="white", edgecolor=INK, lw=0.3, zorder=3)
    for g, geom in zip(dom.ID.values, boxes.values):
        fc = ROLE["boundary"] if g in (FIRST, LAST) else (ROLE["community"] if g in zones else None)
        if fc:
            ax.add_patch(mpl.patches.Polygon(xy2d(geom), closed=True, facecolor=fc, edgecolor=INK,
                                             lw=0.3, zorder=4))
    roads[2008].clip(win).plot(ax=ax, color=C["ROAD"], lw=0.9, zorder=5)

    # Every tenth domain numbered on the SOUND side, the villages named on the OCEAN side beside (README)
    cen = np.c_[boxes.centroid.x.values, boxes.centroid.y.values]
    for i, g in enumerate(dom.ID.values):
        if g == FIRST or g % 10 == 0:
            map_label(ax, dom.geometry.iloc[i].bounds[0] - 500, cen[i, 1], str(g), MAP_TEXT_DARK,
                      ha="right")
    names = {"Buxton": ANN.town_spans["Buxton"], "Avon": ANN.town_spans["Avon"],
             "Salvo": (ANN.village_lines["Salvo"],) * 2, "Waves": (ANN.village_lines["Waves"],) * 2,
             "Rodanthe": (ANN.village_lines["Rodanthe"],) * 2}
    for name, (lo, hi) in names.items():
        sel = dom[(dom.ID >= lo) & (dom.ID <= hi)]
        y = sel.geometry.centroid.y.mean()
        east = sel.total_bounds[2]
        place_label(ax, east + 1100, y, name, MAP_TEXT_DARK, ha="left")
        ax.plot([east + 120, east + 900], [y, y], color=INK_MUTED, lw=0.5, zorder=7)
    # North up: the shore runs north-south, so the groin bar lies east-west
    for name, pos in ANN.groins.items():
        i = int(np.clip(round(pos) - FIRST, 0, len(dom) - 1))
        groin_marker(ax, dom.geometry.iloc[i].bounds[2] + 250,
                     np.interp(pos - FIRST, np.arange(len(cen)), cen[:, 1]), marker="_")
    # The ends and the water bodies the end names
    place_label(ax, bx0 - 300, by0 - 800, "Cape Point", MAP_TEXT_DARK, ha="right", va="top")
    place_label(ax, bx1 + 1100, by1 + 400, "Pea Island", MAP_TEXT_DARK, ha="left", va="bottom")
    for x_, y_, words in ((0.84, 0.22, ("Atlantic", "Ocean")), (0.13, 0.50, ("Pamlico", "Sound"))):
        ax.text(x_, y_, "\n".join(spaced_caps(w) for w in words), transform=ax.transAxes,
                ha="center", va="center", **MAP_TEXT_DARK)

    # North is up here
    north_dart(ax, axes_point(ax, 0.92, 0.48), north=(0.0, 1.0), arrow_m=1600.0, text_kw=MAP_TEXT_DARK)
    ticks_e = np.arange(math.ceil(e0 / 5000) * 5000, e1, 5000)
    ticks_n = np.arange(math.ceil(n0 / 5000) * 5000, n1, 5000)
    ax.set_xticks(ticks_e)
    ax.set_xticklabels([f"{t / 1000:.0f}" for t in ticks_e], fontsize=8)
    ax.set_yticks(ticks_n)
    ax.set_yticklabels([f"{t / 1000:.0f}" for t in ticks_n], fontsize=8)
    ax.tick_params(length=2.5, width=0.5, pad=2, colors=INK)
    ax.set_xlabel("UTM 18N easting (km)", fontsize=9, labelpad=2)
    ax.set_ylabel("northing (km)", fontsize=9, labelpad=2)

    # The locator in the upper left, over the open sound (README)
    side = 1.20
    ax_in = fig.add_axes([(left + 0.30) / fw, (bottom + 0.985 * panel_in - side) / fh,
                          side / fw, side / fh])
    regional_inset(ax_in, outline, vector, lat_side="left", tick_labels=False, locator_c=C["LOCATOR"],
                   locator_lw=2.4)
    # The legend under the locator, in the open sound; the groin bar as drawn here
    handles = framework_handles()
    handles[-1].set_marker("_")
    framework_legend(ax, handles, "upper left",
                     anchor=(0.015, 0.985 - side / panel_in - 0.03))
    _title(ax, 0, "")
    letter_corner(ax_in, 1)

    out = save(fig, fig_path("domain_framework_vertical"), vector=True, dpi=300)
    record_caption(out[0],
        f"The domain framework, north up. (a) The {LAST - FIRST + 1} Barrier3D domain boxes over the island "
        "outline in their true orientation, the reach running south to north up the page from Cape Point "
        "(GIS 1) to the southern end of Pea Island (GIS 90), with the Atlantic to the east; every tenth "
        "domain is numbered on the ocean side and the villages are named on the sound side. The roles the "
        f"hindcast assigns the domains: the interior domains scored against the CoastSat record (GIS "
        f"{SCORE_INTERIOR_GIS[0]}–{SCORE_INTERIOR_GIS[1]}), the end domains that carry the alongshore "
        f"boundary condition (GIS {FIRST}, {LAST}), and the community zones (Buxton 7–8, Avon 21–31, "
        "Salvo–Waves–Rodanthe 68–83) where NC-12 is a maintained street network rather than a relocatable "
        f"road. (The {N_BUFFER} buffer domains beyond each end, which pad the alongshore transport solve, "
        "are not drawn.) The Buxton groin field is the red, white-edged bar drawn offshore. Coordinates are "
        "UTM zone 18N in kilometres, at equal scale in both directions. (b) Hatteras Island (teal) on the "
        "North Carolina coast; the white lines are the 2° graticule (80°, 78° and 76° W; 34° and 36° N); "
        "Natural Earth 10 m coastline and state boundaries. This is the content of `domain_framework` "
        "unturned; a manuscript wants one of the two, not both.")
    plt.close(fig)
    return out[0]


# Study-area imagery and the domain framework as one two-panel figure, one map of the reach
def fig_site_overview(dom, outline, roads, frame, vector):
    window = frame.window(dom, pad_along_km=3.5, pad_sea_km=7.5, pad_sound_km=13.0)
    x0, x1, y0, y1 = window
    aspect = (x1 - x0) / (y1 - y0)
    pw = 0.94
    panel_in = FIG_W_DOUBLE * pw / aspect
    gap_in, top_in, bottom_in = 0.12, 0.06, 0.42
    fh = top_in + 2 * panel_in + gap_in + bottom_in
    fig = plt.figure(figsize=figsize("double", height=fh))
    fw = fig.get_size_inches()[0]
    ax_a = fig.add_axes([0.048, (bottom_in + panel_in + gap_in) / fh, pw, panel_in / fh])
    ax_b = fig.add_axes([0.048, bottom_in / fh, pw, panel_in / fh])
    # (a) imagery, numbered boxes, villages, the locator
    road_c = draw_reach(ax_a, frame, dom, outline, roads[2008], vector, window, piers=False,
                        scalebar=False, arrow=False, numbers=False, label_villages=False,
                        water_labels=False, ends=False, groin=False).road_c
    study_area_labels(ax_a, frame, dom, vector, outline,
                      text_kw=MAP_TEXT_DARK if vector else MAP_TEXT,
                      leader_c=INK_MUTED if vector else "white",
                      shift_m={} if vector else None, ocean_x=0.60, sound_x=0.35)
    scale_bar_km(ax_a, text_kw=MAP_TEXT_DARK if vector else MAP_TEXT)
    ih, iw = 1.3 / fh, 1.3 / fw
    pos = ax_a.get_position()
    ax_in = fig.add_axes([pos.x1 - 0.06 / fw - iw, pos.y1 - 0.06 / fh - ih, iw, ih])
    regional_inset(ax_in, outline, vector, tick_labels=False, locator_c=C["LOCATOR"], locator_lw=2.4)
    letter_corner(ax_in, 1)
    # (a) letter on the inset's top line, the arrow beside the inset; no legend in (a)
    top_row(ax_a, ax_in, None, letter=0, north=frame.north(),
            text_kw=MAP_TEXT_DARK if vector else None)
    # (b) the framework on the same window, with the coordinate frame
    draw_reach(ax_b, frame, dom, outline, roads[2008], vector=True, window=window, label_villages=False,
               water_labels=False, piers=False, scalebar=False, numbers=False, arrow=False,
               ends=False, groin=False)
    draw_role_fills(ax_b, frame, dom)
    utm_frame(ax_b, frame)
    study_area_labels(ax_b, frame, dom, True, outline, text_kw=MAP_TEXT_DARK, leader_c=INK_MUTED,
                      numbers=False, villages=False, water=False, groin_label=False)
    north_dart(ax_b, axes_point(ax_b, 0.78, 0.90), north=frame.north(), text_kw=MAP_TEXT_DARK)
    handles = framework_handles()
    handles[3] = Line2D([], [], color=INK, lw=1.0, label="NC-12, 2008 alignment" + ("" if vector else " (yellow in a)"))
    framework_legend(ax_b, handles, "lower right" if frame.seaward[1] < 0 else "upper right")
    letter_corner(ax_b, 2)
    out = save(fig, fig_path("site_overview"), vector=False, dpi=300)
    record_caption(out[0],
        f"The study reach. (a) Hatteras Island from Cape Point (GIS 1, left) to the southern end of Pea "
        f"Island (GIS 90, right) on Esri World Imagery, turned a quarter turn with north to the right: the "
        f"{LAST - FIRST + 1} Barrier3D domain boxes, 500 m alongshore by 2000 m cross-shore, numbered every "
        "tenth; the NC-12 centreline as digitised on 2008 imagery; the villages; the Buxton groin field "
        "(red, white-edged bar drawn offshore, labelled Buxton groins). (b) The reach on the North Carolina "
        "coast (Natural Earth 10 m), Hatteras Island in teal; the white lines are the 2° graticule (80°, 78° "
        "and 76° W; 34° and 36° N). (c) The same "
        "frame with the roles the hindcast assigns the domains: the interior domains scored against "
        f"CoastSat (GIS {SCORE_INTERIOR_GIS[0]}–{SCORE_INTERIOR_GIS[1]}), the end domains carrying the "
        f"alongshore boundary condition (GIS {FIRST}, {LAST}) and the community zones where roadway "
        "management is off; UTM 18N northing along the frame, easting up it. The 15 buffer domains beyond "
        "each end are not drawn.")
    plt.close(fig)
    return out[0]


# Island width and highest cell per domain, from the extraction
def fig_domain_metrics():
    gis, width, crest = domain_metrics()
    fig = plt.figure(figsize=figsize("double", height=3.6))
    ax_a = fig.add_axes([0.085, 0.585, 0.895, 0.355])
    ax_b = fig.add_axes([0.085, 0.125, 0.895, 0.355])
    for ax in (ax_a, ax_b):
        ax.set_xlim(FIRST - 0.5, LAST + 0.5)
        open_frame(ax)
        ax.grid(axis="y", color=GRID_C, lw=0.5)
    ax_a.plot(gis, width, color=INK, lw=1.0)
    ax_a.set_ylabel("island width (m)")
    ax_a.set_ylim(0, None)
    ax_a.tick_params(labelbottom=False)
    ax_b.plot(gis, crest, color=INK, lw=1.0)
    ax_b.set_ylabel("highest cell (m NAVD88)")
    ax_b.set_ylim(0, None)
    ax_b.set_xlabel(DOMAIN_AXIS_LABEL)
    for i, ax in ((0, ax_a), (1, ax_b)):
        town_bands(ax, label=(ax is ax_a))
        structures(ax, label=(ax is ax_a))
        _title(ax, i, "")
    out = save(fig, fig_path("domain_metrics"), vector=True, dpi=300)
    record_caption(out[0],
        f"What the {TOPO_PRODUCT} extraction gives each domain. (a) Island width: the median across the 50 "
        "alongshore rows of the count of cells above 0 m NAVD88, so capped at the 2000 m box where the "
        "island is wider (Cape Point). (b) The highest cell in each row, median across rows, from the same "
        "arrays. Village bands and the groin and piers are marked as in every alongshore figure; GIS 1 is "
        "Cape Point, 90 the southern end of Pea Island.")
    plt.close(fig)
    return out[0]


# One domain: its box on imagery, the 10 m array, road cells, mean profile
def fig_domain_grid(dom, roads, vector, gis=EXAMPLE_GIS):
    d, _ = tv.npy_dirs(TOPO_PRODUCT)
    a = np.load(d / f"domain_{gis}.npy")                    # (alongshore rows, cross-shore cols), m NAVD88
    mask = np.load(tv.road_mask_file(2008, gis)).astype(bool)
    nrow, ncol = a.shape
    water = a <= -9.9
    prof = np.where(water, np.nan, a).mean(0)               # mean profile over land cells
    crest_col = int(np.nanargmax(prof))
    road_cols = np.where(mask.any(0))[0]
    poly = dom[dom.ID == gis].geometry.iloc[0]
    bx0, by0, bx1, by1 = poly.bounds
    assert abs((bx1 - bx0) - ncol * CELL_M) < 1 and abs((by1 - by0) - nrow * CELL_M) < 1, \
        f"domain box {poly.bounds} is not the array's {ncol}x{nrow} cells"

    fig = plt.figure(figsize=figsize("double", height=4.3))
    ax_a = fig.add_axes([0.03, 0.56, 0.44, 0.41])
    ax_b = fig.add_axes([0.53, 0.56, 0.44, 0.41])
    ax_c = fig.add_axes([0.075, 0.095, 0.90, 0.36])

    # (a) the box, north up
    padx, pady = 250, 450
    map_axes(ax_a, (bx0 - padx, bx1 + padx, by0 - pady, by1 + pady))
    if vector:
        ax_a.set_facecolor(C["WATER"])
        road_c, box_c, text_c = C["ROAD"], INK, INK
    else:
        img, (l, r, b, t) = tiles((bx0 - padx, by0 - pady, bx1 + padx, by1 + pady), 15, imagery_source())
        ax_a.imshow(img, extent=(l, r, b, t), zorder=0, interpolation="bilinear")
        credit(ax_a, "Imagery: Esri World Imagery")
        road_c, box_c, text_c = ROAD_ON_IMAGERY, "white", "white"
    ax_a.add_patch(mpl.patches.Polygon(xy2d(poly), closed=True, facecolor="none", edgecolor=box_c,
                                       lw=1.0, zorder=4))
    roads[2008].plot(ax=ax_a, color=road_c, lw=1.0, zorder=5)
    tk = MAP_TEXT_DARK if vector else MAP_TEXT
    map_label(ax_a, (bx0 + bx1) / 2, by1 + 60, f"GIS {gis}: 2000 m × 500 m domain box", tk, va="bottom")
    scale_bar_km(ax_a, 500, segments=1, unit="m", text_kw=tk, y=0.10)
    north_dart(ax_a, axes_point(ax_a, 0.94, 0.78), north=(0.0, 1.0), arrow_m=220.0, text_kw=tk)
    letter_corner(ax_a, 0)

    # (b) the array the model reads
    cmap, norm, bounds = elevation_cmap()
    ext = (0, ncol * CELL_M, nrow * CELL_M, 0)
    im = ax_b.imshow(a, cmap=cmap, norm=norm, extent=ext, interpolation="nearest", zorder=1)
    rr, cc = np.where(mask)
    for r_, c_ in zip(rr, cc):
        ax_b.add_patch(Rectangle((c_ * CELL_M, r_ * CELL_M), CELL_M, CELL_M, facecolor=C["ROAD"],
                                 edgecolor="none", zorder=3))
    ax_b.axvline((crest_col + 0.5) * CELL_M, color=INK, lw=0.7, ls=(0, (2, 1.5)), zorder=4)
    ax_b.set_aspect("equal")
    ax_b.set_xlim(0, ncol * CELL_M)
    ax_b.set_ylim(nrow * CELL_M, 0)
    ax_b.set_xticks([0, 500, 1000, 1500, 2000])
    ax_b.set_yticks([0, 250, 500])
    ax_b.set_xlabel("cross-shore (m)  →  ocean")
    ax_b.set_ylabel("alongshore (m)")
    spines_for_image(ax_b)
    cb = fig.colorbar(im, ax=ax_b, orientation="horizontal", fraction=0.06, pad=0.32, aspect=35,
                      ticks=bounds[1:-1])
    cb.set_label("elevation (m NAVD88)", fontsize=8)
    cb.ax.tick_params(labelsize=8)
    cb.outline.set_linewidth(0.5)
    ax_b.text(0.02, 0.06, "sound", transform=ax_b.transAxes, fontsize=8, color=INK, ha="left", va="bottom")
    ax_b.text(0.98, 0.06, "ocean", transform=ax_b.transAxes, fontsize=8, color=INK, ha="right", va="bottom")
    _title(ax_b, 1, "")

    # (c) the mean cross-shore profile
    x = (np.arange(ncol) + 0.5) * CELL_M
    ax_c.fill_between(x, -1.0, np.nan_to_num(prof, nan=-1.0), color=PROFILE_GROUND, lw=0, zorder=1)
    ax_c.plot(x, prof, color=INK, lw=1.0, zorder=3)
    ax_c.axhline(0, color=INK_MUTED, lw=0.6, zorder=2)
    rx0, rx1 = road_cols.min() * CELL_M, (road_cols.max() + 1) * CELL_M
    ax_c.axvspan(rx0, rx1, color=C["ROAD"], alpha=0.25, lw=0, zorder=2)
    ax_c.axvline((crest_col + 0.5) * CELL_M, color=INK, lw=0.7, ls=(0, (2, 1.5)), zorder=4)
    ax_c.annotate("dune crest", xy=((crest_col + 0.5) * CELL_M, np.nanmax(prof)), xytext=(-6, 4),
                  textcoords="offset points", ha="right", va="top", fontsize=8, color=INK)
    ax_c.text((rx0 + rx1) / 2, np.nanmax(prof) * 0.62, "NC-12", ha="center", va="bottom", fontsize=8,
              color=INK, rotation=90)
    land_x = x[~np.isnan(prof)]
    ax_c.text(land_x[0] + 10, 0.08, "0 m NAVD88", fontsize=8, color=INK_MUTED, va="bottom")
    ax_c.set_xlim(0, ncol * CELL_M)
    ax_c.set_ylim(-1.0, max(np.nanmax(prof) + 1.0, 4.0))
    ax_c.set_xlabel("cross-shore (m)  →  ocean")
    ax_c.set_ylabel("mean elevation (m NAVD88)")
    open_frame(ax_c)
    ax_c.grid(axis="y", color=GRID_C, lw=0.5)
    _title(ax_c, 2, "")

    out = save(fig, fig_path("domain_grid"), vector=False, dpi=300)
    record_caption(out[0],
        f"One domain as the model reads it (GIS {gis}, north of Avon, {TOPO_PRODUCT} product). "
        "(a) The 2000 m by 500 m box the domain's elevation array is cut from, with the 2008 NC-12 "
        f"centreline, north up; the 500 m scale bar is 50 model cells. (b) The array itself: {nrow} alongshore rows by {ncol} cross-shore columns "
        "of 10 m cells, m NAVD88, the sound at the left and the ocean at the right; the cells the 2008 road "
        "alignment rasterises onto are black, and the dashed line is the column of the mean-profile dune "
        "crest. Water cells are the extractor's -10 m sentinel. (c) The mean cross-shore profile over land "
        "cells, with the road corridor shaded and the crest marked. The dune search, the two dune rows and "
        "interior row 0 are placed on this array by the extractor; the road setback the model spends is "
        "measured from that row.")
    plt.close(fig)
    return out[0]


# Figure 4: forcing timeline
from site_layer import hat_env_forcings as _env  # noqa: E402
STORM_ROOT = _env.HINDCAST_STORMS
RSLR_RECORD = _env.RSLR_RECORD_FILE
RSLR_FITS = _env.RSLR_RATES_CSV
WINDOW_COLOUR = {(1984, 2004): C_1984, (2004, 2024): C_1997,
                 (1996, 2010): "#ef8a62", (2010, 2024): "#67a9cf"}

# The two hindcast chains, one figure each (README)
CHAINS = {
    "forcing_timeline_1984": {"windows": [(1984, 2004), (2004, 2024)], "years": (1983, 2025)},
    "forcing_timeline_1996": {"windows": [(1996, 2010), (2010, 2024)], "years": (1995, 2025)},
}


# One row per storm, 1984-2023, with its calendar year
def storm_record():
    import pandas as pd
    # the series the hindcast runs on (hat_env_forcings.DEFAULT_STORM_VARIANT)
    from site_layer import hat_env_forcings as _envf
    a = pd.read_csv(_envf.storm_summary_file(1984, 2004))
    b = pd.read_csv(_envf.storm_summary_file(2004, 2024))
    s = pd.concat([a[a.calendar_year < 2004], b], ignore_index=True)
    s["Rhigh_m"] = s["Rhigh"] * 10.0          # Barrier3D stores runup in decametres
    return s


# The Duck sea-level record for the forcing timeline
def sea_level_record():
    import pandas as pd
    lines = RSLR_RECORD.read_text(encoding="utf-8").splitlines()
    first = next(i for i, ln in enumerate(lines) if ln.strip().startswith("Year"))
    d = pd.read_csv(RSLR_RECORD, skiprows=first + 1, header=None, usecols=[0, 1, 2],
                    names=["year", "month", "msl"])
    d = d.dropna()
    d["t"] = d.year + (d.month - 0.5) / 12
    return d


# The chain's windows as horizontal bars
def windows_bar(ax, windows):
    for i, (a, b) in enumerate(windows):
        ax.barh(i, b - a, left=a, height=0.62, color=WINDOW_COLOUR[(a, b)], lw=0)
        ax.text(a + 0.4, i, f"{a}–{b}", ha="left", va="center", fontsize=8, color="white",
                fontweight="bold")
    ax.set_ylim(len(windows) - 0.4, -0.6)
    ax.set_yticks([])
    ax.set_ylabel("hindcast\nwindows")
    for side in ("top", "right", "left"):
        ax.spines[side].set_visible(False)


# The forcing record over ONE chain
def fig_forcing_timeline(name):
    chain = CHAINS[name]
    windows, years = chain["windows"], chain["years"]
    s = storm_record()
    sl = sea_level_record()
    import pandas as pd
    fits = pd.read_csv(RSLR_FITS)

    # Cut to the years the chain covers, not the axis margins
    span = (windows[0][0], windows[-1][1])

    def in_span(y):
        return span[0] <= y <= span[1]

    s = s[s.calendar_year.between(*span)]
    sl = sl[sl.t.between(*span)]
    fits = fits[[(int(r.start_year), int(r.end_year)) in windows for r in fits.itertuples()]]

    fig = plt.figure(figsize=figsize("double", height=6.4))
    hs = [0.40, 1.0, 0.9, 1.15, 1.35]                      # panel heights, relative
    gap, top, bottom = 0.24, 0.28, 0.42
    fh = fig.get_size_inches()[1]
    axes = []
    y = fh - top
    for h in hs:
        h_in = h * (fh - top - bottom - gap * (len(hs) - 1)) / sum(hs)
        y -= h_in
        axes.append(fig.add_axes([0.085, y / fh, 0.895, h_in / fh]))
        y -= gap
    ax_w, ax_n, ax_r, ax_s, ax_m = axes
    for ax in axes:
        ax.set_xlim(*years)
        open_frame(ax)
    for ax in axes[:-1]:
        ax.tick_params(labelbottom=False)

    windows_bar(ax_w, windows)
    _title(ax_w, 0, "")

    # storms per year, and the largest runup of the year
    counts = s.groupby("calendar_year").size()
    ax_n.bar(counts.index, counts.values, width=0.8, color=C["BASE"], lw=0)
    ax_n.set_ylabel("storms\nper year")
    ax_n.grid(axis="y", color=GRID_C, lw=0.5)
    _title(ax_n, 1, "")
    rmax = s.groupby("calendar_year")["Rhigh_m"].max()
    ax_r.vlines(rmax.index, 0, rmax.values, color=C["BASE"], lw=0.8)
    ax_r.plot(rmax.index, rmax.values, "o", ms=2.6, color=INK, mec="none")
    ax_r.set_ylabel("largest R$_{high}$\n(m above MHW)")
    ax_r.set_ylim(0, None)
    ax_r.grid(axis="y", color=GRID_C, lw=0.5)
    _title(ax_r, 2, "")

    # sea level: the monthly record and the chain's two fitted trends
    ax_s.plot(sl.t, sl.msl, color=C["BASE_FILL"], lw=0.6, zorder=1)
    ann = sl.groupby("year")["msl"].mean()
    ax_s.plot(ann.index + 0.5, ann.values, color=INK, lw=1.0, zorder=2)
    for _, f in fits.iterrows():
        a, b = int(f.start_year), int(f.end_year)
        t = np.array([a, b], float)
        ax_s.plot(t, f.slope_m_yr * t + f.intercept_m, color=WINDOW_COLOUR[(a, b)], lw=1.6, zorder=3,
                  label=f"{a}–{b}: {f.slope_mm_yr:.1f} ± {f.ci95_mm_yr:.1f} mm/yr")
    ax_s.set_ylabel("monthly mean\nsea level (m)")
    ax_s.grid(axis="y", color=GRID_C, lw=0.5)
    ax_s.legend(loc="upper left", ncol=2, frameon=False, fontsize=8, handlelength=1.4,
                columnspacing=1.0)
    _title(ax_s, 3, "")

    # management: what happened where, GIS domain against year
    ax_m.set_ylim(FIRST - 0.5, LAST + 0.5)
    for nm, (lo, hi) in ANN.town_spans.items():
        ax_m.axhspan(lo - 0.5, hi + 0.5, color="0.94", lw=0, zorder=0)
        ax_m.text(years[0] + 0.4, (lo + hi) / 2, nm, ha="left", va="center", fontsize=8,
                  color=INK_MUTED, zorder=1)
    drew_reloc = drew_bridge = drew_fill = False
    for e in HATTERAS_ROAD_EVENTS:
        if not in_span(e.year):
            continue                        # an event outside this chain's span
        if hasattr(e, "displacement_m"):
            gis = sorted(e.displacement_m)
            ax_m.add_patch(Rectangle((e.year - 0.5, gis[0] - 0.5), 1.0, gis[-1] - gis[0] + 1,
                                     facecolor=C["ROAD"], edgecolor="none", zorder=3))
            drew_reloc = True
            # An early event's label would run off the left of the axis, so it goes to the right of its bar (README)
            side = "left" if e.year - years[0] < 8 else "right"
            ax_m.text(e.year + (0.8 if side == "left" else -0.8), (gis[0] + gis[-1]) / 2,
                      f"{e.year}, GIS {gis[0]}–{gis[-1]}", ha=side, va="center", fontsize=8,
                      color=INK, zorder=4)
        else:
            gis = sorted(e.gis_domains)
            ax_m.add_patch(Rectangle((e.year - 0.5, gis[0] - 0.5), years[1] - e.year + 0.5,
                                     gis[-1] - gis[0] + 1, facecolor="none", edgecolor=INK_MUTED,
                                     hatch="////", lw=0.5, zorder=3))
            drew_bridge = True
            ax_m.text(e.year - 0.8, (gis[0] + gis[-1]) / 2, f"bridge, {e.year}", ha="right",
                      va="center", fontsize=8, color=INK, zorder=4)
    for pr in HATTERAS_NOURISHMENT_PROJECTS:
        if not in_span(pr.year):
            continue
        gis = sorted(pr.gis_domains)
        ax_m.add_patch(Rectangle((pr.year - 0.5, gis[0] - 0.5), 1.0, gis[-1] - gis[0] + 1,
                                 facecolor=C["ADDED"], edgecolor="none", zorder=3))
        drew_fill = True
        ax_m.text(pr.year - 0.8, (gis[0] + gis[-1]) / 2, f"{pr.name.split()[0]}, {pr.year}",
                  ha="right", va="center", fontsize=8, color=INK, zorder=4)
    ax_m.set_ylabel(DOMAIN_AXIS_LABEL.replace(" (", "\n("))
    ax_m.set_xlabel("year")
    ax_m.set_yticks([1, 30, 60, 90])
    handles = []
    if drew_reloc:
        handles.append(Patch(facecolor=C["ROAD"], label="NC-12 relocation"))
    if drew_fill:
        handles.append(Patch(facecolor=C["ADDED"], label="beach nourishment"))
    if drew_bridge:
        handles.append(Patch(facecolor="none", edgecolor=INK_MUTED, hatch="////",
                             label="road removed (bridge)"))
    ax_m.legend(handles=handles, loc="center", ncol=3, frameon=False, fontsize=8,
                handlelength=1.4, columnspacing=1.0, bbox_to_anchor=(0.5, 0.52))
    _title(ax_m, 4, "")

    out = save(fig, fig_path(name), vector=True, dpi=300)
    w = " and ".join(f"{a}–{b}" for a, b in windows)
    reloc = "NC-12 relocations (black, the domains moved), " if drew_reloc else ""
    bridge = ("and the Rodanthe bridge after which the road is removed from GIS 82–88 "
              "(hatched) " if drew_bridge else "")
    record_caption(out[0],
        f"The forcing record over the {windows[0][0]} chain, {span[0]}–{span[1]}. "
        f"(a) The two hindcast windows, {w}; they tile at {windows[0][1]}, and an end year is a "
        "boundary, not a simulated year. (b) Storms per calendar year in the hindcast storm series "
        "(3-env-forcings/3-storms/hindcast_storms, v3, 72 h separation), and (c) the largest storm "
        "runup R_high of each year, m above MHW. (d) Monthly mean sea level at the Duck gauge "
        "(NOAA 8651370, seasonal cycle removed; grey monthly, black annual mean) with the linear "
        "trend fitted over each window, in the window's colour; the model rounds these to "
        f"0.001 m/yr. (e) {reloc}beach nourishment (orange, the domains filled) {bridge}as the "
        "model prescribes them, against the village spans. The companion figure draws the other "
        "chain, on its own span: the two share no axis.")
    plt.close(fig)
    return out[0]


# The observed CoastSat rates of each window
def fig_observed_rates():
    import pandas as pd
    from site_layer import hat_observed_rates as obs
    pairs = (((1984, 2004), (2004, 2024)), ((1996, 2010), (2010, 2024)))
    fig = plt.figure(figsize=figsize("double", height=3.9))
    ax_a = fig.add_axes([0.085, 0.55, 0.895, 0.385])
    ax_b = fig.add_axes([0.085, 0.105, 0.895, 0.385])
    for i, (ax, pair) in enumerate(zip((ax_a, ax_b), pairs)):
        ax.set_xlim(FIRST - 0.5, LAST + 0.5)
        ax.axhline(0, color=INK_MUTED, lw=0.6, zorder=1)
        for w, colour, fill in zip(pair, (C_1984, C_1997), (C_1984_FILL, C_1997_FILL)):
            d = pd.read_csv(obs.domain_csv(*w))
            g = d.domain_number.values
            ax.fill_between(g, d.mean_lrr - d.std_lrr, d.mean_lrr + d.std_lrr, color=fill, alpha=0.5,
                            lw=0, zorder=2)
            ax.plot(g, d.mean_lrr, color=colour, lw=1.1, zorder=3, label=f"{w[0]}–{w[1]}")
        ax.set_ylabel("shoreline change (m/yr)")
        ax.grid(axis="y", color=GRID_C, lw=0.5)
        open_frame(ax)
        ax.legend(loc="lower right", ncol=2, frameon=False, fontsize=8, handlelength=1.6)
        town_bands(ax, label=(i == 0))
        structures(ax, label=(i == 0))
        _title(ax, i, "")
    ax_a.tick_params(labelbottom=False)
    ax_b.set_xlabel(DOMAIN_AXIS_LABEL)
    lo = min(ax_a.get_ylim()[0], ax_b.get_ylim()[0])
    hi = max(ax_a.get_ylim()[1], ax_b.get_ylim()[1])
    ax_a.set_ylim(lo, hi)
    ax_b.set_ylim(lo, hi)
    out = save(fig, fig_path("observed_rates"), vector=True, dpi=300)
    record_caption(out[0],
        "Observed shoreline change, the calibration target. Linear regression rate of the CoastSat "
        "shoreline (m/yr, negative landward) per domain, the mean over the domain's transects with a "
        "band of one standard deviation across them, for (a) the two matrix periods 1984–2004 and "
        "2004–2024 and (b) the two overlapping windows 1996–2010 and 2010–2024. The earlier window of "
        "each pair is red, the later blue. Village bands, the Buxton groin (solid) and the piers "
        "(dotted) as in every alongshore figure; GIS 1 is Cape Point, 90 the southern end of Pea Island.")
    plt.close(fig)
    return out[0]


# Figure 6: the dune lines, 1984-2023
DUNE_VINTAGES = (1984, 1997, 2004, 2009, 2023)
VINTAGE_COLOUR = {1984: C_1984, 1997: "#ef8a62", 2004: "#7f7f7f", 2009: "#67a9cf", 2023: C_1997}
ZOOMS = {"Buxton": (4, 7), "Avon": (29, 32), "GIS 50–53": (50, 53), "Rodanthe": (82, 85)}   # detail panels


# Every dune-line vintage, in UTM
def dune_lines():
    return {v: gpd.read_file(tv.duneline_geojson(v)).to_crs(CRS) for v in DUNE_VINTAGES}


# Mean easting of each vintage's line inside each domain box, m
def line_position(lines, dom):
    pos = {}
    for v, gs in lines.items():
        line = gs.geometry.union_all()
        xs = []
        for geom in dom.geometry:
            seg = line.intersection(geom)
            xs.append(np.nan if seg.is_empty else seg.centroid.x)
        pos[v] = np.array(xs)
    return pos


# The digitised dune lines on the reach, with zooms
def fig_dune_lines(dom, roads, vector):
    lines = dune_lines()
    pos = line_position(lines, dom)
    fig = plt.figure(figsize=figsize("double", height=5.4))
    ax_a = fig.add_axes([0.085, 0.63, 0.895, 0.30])
    zoom_axes = [fig.add_axes([0.035 + i * 0.2425, 0.03, 0.2125, 0.50]) for i in range(len(ZOOMS))]

    ax_a.set_xlim(FIRST - 0.5, LAST + 0.5)
    ax_a.axhline(0, color=VINTAGE_COLOUR[1984], lw=1.0, zorder=2, label="1984")
    for v in DUNE_VINTAGES[1:]:
        ax_a.plot(dom.ID.values, pos[v] - pos[1984], color=VINTAGE_COLOUR[v], lw=1.0, zorder=3, label=str(v))
    ax_a.set_ylabel("dune line relative\nto 1984 (m, + seaward)")
    ax_a.set_xlabel(DOMAIN_AXIS_LABEL)
    ax_a.grid(axis="y", color=GRID_C, lw=0.5)
    open_frame(ax_a)
    ax_a.legend(loc="lower center", bbox_to_anchor=(0.5, 1.0), ncol=5, frameon=False, fontsize=8,
                handlelength=1.6, columnspacing=1.2)
    town_bands(ax_a)
    structures(ax_a)
    for name, (lo, hi) in ZOOMS.items():
        ax_a.axvspan(lo - 0.5, hi + 0.5, facecolor="none", edgecolor=INK, lw=0.6, ls=(0, (2, 1.5)), zorder=4)
    _title(ax_a, 0, "")

    for i, (ax, (name, (lo, hi))) in enumerate(zip(zoom_axes, ZOOMS.items()), start=1):
        boxes = dom[(dom.ID >= lo) & (dom.ID <= hi)]
        x0, y0, x1, y1 = boxes.total_bounds
        # the lines sit in the seaward half of the boxes; frame that
        xs = np.concatenate([pos[v][(dom.ID >= lo) & (dom.ID <= hi)] for v in DUNE_VINTAGES])
        cx_ = np.nanmean(xs)
        half_w = (y1 - y0) / 2 * 0.60
        map_axes(ax, (cx_ - half_w, cx_ + half_w, y0, y1))
        if vector:
            ax.set_facecolor("white")
            road_c = C["ROAD"]
        else:
            img, (l, r, b, t) = tiles((cx_ - half_w, y0, cx_ + half_w, y1), 16, imagery_source())
            ax.imshow(img, extent=(l, r, b, t), zorder=0, interpolation="bilinear")
            road_c = ROAD_ON_IMAGERY
        tk = MAP_TEXT_DARK if vector else MAP_TEXT
        for geom in boxes.geometry:
            ax.add_patch(mpl.patches.Polygon(xy2d(geom), closed=True, facecolor="none",
                                             edgecolor="white" if not vector else INK_MUTED, lw=0.4, zorder=2))
        for v in DUNE_VINTAGES:
            lines[v].plot(ax=ax, color=VINTAGE_COLOUR[v], lw=1.1, zorder=4)
        roads[2008].plot(ax=ax, color=road_c, lw=0.8, zorder=3)
        # The numbers on the OCEAN side, clear of the panel letter in the other corner, and only where (README)
        for g, geom in zip(boxes.ID.values, boxes.geometry):
            yc = geom.centroid.y
            if y0 + 90 < yc < y1 - 90:
                map_label(ax, cx_ + half_w * 0.94, yc, str(g), tk, ha="right")
        # a place name in italics; a domain range is not a place, so it stays upright
        (place_label if name[:3] != "GIS" else map_label)(ax, 0.5, 0.965, name, tk, va="top",
                                                          transform=ax.transAxes, zorder=9)
        scale_bar_km(ax, 250, segments=1, unit="m", text_kw=tk, y=0.08)
        # The arrow sits between two domain numbers, not beside one
        north_dart(ax, axes_point(ax, 0.90, 0.46), north=(0.0, 1.0), arrow_m=110.0, text_kw=tk)
        letter_corner(ax, i)

    if not vector:
        credit_figure(fig, "Imagery: Esri World Imagery")
    out = save(fig, fig_path("dune_lines"), vector=False, dpi=300)
    record_caption(out[0],
        "The digitised dune line through the record. (a) Each vintage's position relative to the 1984 "
        "line, per domain: the mean easting of the line inside the domain box, positive seaward; the "
        "coast trends 8° from north so easting differences are cross-shore movement to within one per "
        "cent. 1997, 2004, 2009 and 2023 are the imagery vintages behind the 1996, 2004, 2010 and 2024 "
        "period boundaries (hat_topo_version.DUNE_LINE_FOR_YEAR); the 1967 line covers only Cape Point "
        "and is not drawn. Dashed boxes mark the detail panels: (b) Buxton, GIS 4–7; (c) Avon, GIS 29–32, "
        "where the 2004 line lies seaward of every other; (d) GIS 50–53, the largest 2004 retreat; (e) "
        "Rodanthe, GIS 82–85, the S-curves. Each shows the five lines with the domain boxes numbered and "
        "the 2008 NC-12 centreline, north up. " + ("" if vector else "Basemap: Esri World Imagery."))
    plt.close(fig)
    return out[0]


# The Barrier3D parameters the runs use, from the yaml
def model_params():
    import yaml
    with open(INIT / "Hatteras-CASCADE-parameters.yaml", encoding="utf-8") as fh:
        return yaml.safe_load(fh)


# The extractor's output for one domain
def processed_domain(gis, product=TOPO_PRODUCT):
    root = tv.dune_topo_root(product)
    ver = (root / "CURRENT").read_text(encoding="utf-8").strip()
    topo = np.load(root / ver / "topography" / f"domain_{gis}_topography.npy")
    dune = np.load(root / ver / "dunes" / f"domain_{gis}_dune.npy")
    return topo, dune, ver


# One domain's road setback (m) for a start year
def setback_m(gis, year=2004):
    a = np.loadtxt(tv.road_setback_file(year), delimiter=",")
    return float(a[1][a[0] == gis][0])


# One domain as Barrier3D holds it: profile and plan on one cross-shore axis
def fig_domain_schematic(gis=EXAMPLE_GIS):
    prm = model_params()
    topo, dune, ver = processed_domain(gis)
    sb = setback_m(gis)
    dam = 10.0
    interior = topo.mean(1) * dam                    # m MHW, ocean first
    dune_h = dune.mean() * dam
    berm = prm["BermEl"] - prm["MHW"]                # BermEl is m NAVD88 in the yaml; MHW 0.36
    bay = -prm["BayDepth"]
    n_dune = 2
    beach_w = 60.0                                   # schematic berm width
    sf_show = 260.0                                  # how much shoreface to draw
    sf_slope = prm["DShoreface"] / prm["LShoreface"]
    x_dune0 = beach_w                                # the first dune row's seaward edge
    x_int0 = beach_w + n_dune * CELL_M               # interior row 0's seaward edge
    x_end = x_int0 + len(interior) * CELL_M
    xlim = (-sf_show, x_end + 130)

    fig = plt.figure(figsize=figsize("double", height=4.6))
    ax_a = fig.add_axes([0.075, 0.50, 0.86, 0.46])
    ax_b = fig.add_axes([0.075, 0.10, 0.86, 0.30])
    cax = fig.add_axes([0.925, 0.10, 0.012, 0.30])

    # (a) the cross-section, ocean at the left
    x_int = x_int0 + (np.arange(len(interior)) + 0.5) * CELL_M
    x_all = np.concatenate([[-sf_show, 0, beach_w], [x_dune0, x_dune0, x_int0, x_int0], x_int, [x_end, x_end + 120]])
    z_all = np.concatenate([[-sf_show * sf_slope, 0, berm], [berm, dune_h, dune_h, interior[0]], interior, [interior[-1], bay]])
    ax_a.axhspan(bay - 1.0, 0, color=C["WATER"], alpha=0.45, lw=0, zorder=0)
    ax_a.fill_between(x_all, bay - 1.0, z_all, color="#e0cf96", lw=0, zorder=1)
    ax_a.plot(x_all, z_all, color=INK, lw=1.0, zorder=3)
    ax_a.axhline(0, color=C_1997, lw=0.7, zorder=2)
    x_road = x_int0 + sb
    z_road = np.interp(x_road, x_int, interior)
    ax_a.add_patch(Rectangle((x_road - CELL_M / 2, z_road), CELL_M, 0.25, facecolor=C["ROAD"], lw=0, zorder=5))
    lab = dict(fontsize=8, color=INK, zorder=6)
    leader = dict(arrowstyle="-", lw=0.5, color=INK_MUTED)
    ax_a.text(-sf_show + 8, bay - 0.7, "shoreface (BRIE)\nslope D / L", ha="left", va="bottom", **lab)
    ax_a.text(-sf_show + 8, 0.12, "MHW", ha="left", va="bottom", fontsize=8, color=C_1997, zorder=6)
    ax_a.annotate("shoreline x$_s$", xy=(0, 0), xytext=(-70, 2.4), ha="right", arrowprops=leader, **lab)
    ax_a.annotate(f"berm, {berm:.2f} m", xy=(beach_w / 2, berm / 2), xytext=(beach_w + 40, -1.5), ha="left", arrowprops=leader, **lab)
    ax_a.annotate(f"dune, {n_dune} rows ({dune_h:.1f} m here)", xy=(x_dune0 + CELL_M, dune_h), xytext=(x_int0 + 130, dune_h + 0.9), ha="left", arrowprops=leader, **lab)
    ax_a.annotate("interior row 0", xy=(x_int[0], interior[0]), xytext=(x_int0 + 130, 4.1), ha="left", arrowprops=leader, **lab)
    ax_a.annotate(f"NC-12, setback {sb:.0f} m from row 0", xy=(x_road, z_road + 0.25), xytext=(x_road + 80, 2.5), ha="left", arrowprops=leader, **lab)
    ax_a.text(x_int0 + 1000, dune_h - 0.7, f"interior, {len(interior)} rows × {CELL_M:.0f} m (Barrier3D)",
              ha="center", **lab)
    ax_a.text(x_end + 60, bay + 0.3, f"bay, {prm['BayDepth']:.0f} m", ha="center", va="bottom", **lab)
    ax_a.set_xlim(*xlim)
    ax_a.set_ylim(bay - 1.0, dune_h + 2.6)
    ax_a.set_ylabel("elevation (m MHW)")
    ax_a.tick_params(labelbottom=False)
    open_frame(ax_a)
    _title(ax_a, 0, "")

    # (b) the same domain in plan, on the same cross-shore axis
    cmap, norm, bounds = elevation_cmap()
    grid = np.vstack([np.tile(dune * dam, (n_dune, 1)), topo * dam]).T      # (alongshore, cross-shore)
    ext = (x_dune0, x_end, 0, grid.shape[0] * CELL_M)
    im = ax_b.imshow(grid, cmap=cmap, norm=norm, extent=ext, interpolation="nearest", aspect="auto", zorder=1)
    ax_b.axvline(x_int0, color=INK, lw=0.6, ls=(0, (2, 1.5)), zorder=4)
    ax_b.add_patch(Rectangle((x_road - CELL_M / 2, 0), CELL_M, grid.shape[0] * CELL_M, facecolor=C["ROAD"], lw=0, zorder=3))
    ax_b.set_xlim(*xlim)
    ax_b.set_ylim(0, grid.shape[0] * CELL_M)
    ax_b.set_yticks([0, 250, 500])
    ax_b.set_xlabel("cross-shore (m)  ←  ocean")
    ax_b.set_ylabel("alongshore (m)")
    ax_b.set_facecolor("white")
    open_frame(ax_b)
    ax_b.text(x_dune0 - 8, grid.shape[0] * CELL_M / 2, "dune rows", ha="right", va="center", fontsize=8, color=INK, rotation=90)
    ax_b.text(x_road + 12, grid.shape[0] * CELL_M - 15, "NC-12", ha="left", va="top", fontsize=8, color=INK)
    cb = fig.colorbar(im, cax=cax, ticks=bounds[1:-1])
    cax.set_title("m MHW", fontsize=9, color=INK, pad=3)
    cb.ax.tick_params(labelsize=8)
    cb.outline.set_linewidth(0.5)
    _title(ax_b, 1, "")

    out = save(fig, fig_path("domain_schematic"), vector=True, dpi=300)
    record_caption(out[0],
        f"How one domain is built (GIS {gis}, {TOPO_PRODUCT} dune-topo {ver}). (a) The cross-section the "
        "coupled model carries, ocean at the left, elevations in m above MHW: the BRIE shoreface of slope "
        f"D/L = {prm['DShoreface']:.1f} m / {prm['LShoreface']:.0f} m (only its top is drawn), the berm at "
        f"BermEl = {prm['BermEl']:.1f} m NAVD88 ({berm:.2f} m MHW), two dune rows carrying the extracted "
        f"dune height (mean {dune_h:.1f} m here), the {len(interior)}-row Barrier3D interior from the "
        f"extraction (alongshore mean shown), and the bay at {prm['BayDepth']:.0f} m depth. NC-12 sits "
        f"{sb:.0f} m landward of interior row 0, the setback the roadway manager reads. The beach width is "
        "schematic; the model sets it from the shoreline position. (b) The same domain in plan on the "
        "same cross-shore axis, 50 cells (500 m) alongshore: the dune rows left of the dashed line, the "
        "interior cells in the elevation classes, the road cells black.")
    plt.close(fig)
    return out[0]


# Figure 8: management footprint

# True when `half` ("sea" or "sound") is the LOWER half of a domain box on the page
def half_is_lower(frame, half):
    return (frame.seaward[1] < 0) == (half == "sea")


# A legend entry that IS a domain box with one half marked, drawn the way `tint_domains` marks it
class HalfBox:

    def __init__(self, color, lower, alpha, hatch=None, label="", whole=False):
        self.color, self.lower, self.alpha, self.hatch = color, lower, alpha, hatch
        self.whole = whole
        self._label = label

    def get_label(self):
        return self._label


# Draws a HalfBox
class HalfBoxHandler(HandlerBase):

    def create_artists(self, legend, orig, xdescent, ydescent, width, height,
                       fontsize, trans):
        x, y = -xdescent, -ydescent
        if orig.hatch:
            kw = dict(facecolor="none", edgecolor=INK, hatch=orig.hatch, lw=0.0)
        else:
            kw = dict(facecolor=orig.color, edgecolor="none", alpha=orig.alpha)
        h = height if orig.whole else height / 2
        marked = mpl.patches.Rectangle(
            (x, y if (orig.whole or orig.lower) else y + height / 2), width, h, **kw)
        box = mpl.patches.Rectangle((x, y), width, height, facecolor="none",
                                    edgecolor=INK, lw=0.35)
        for a in (marked, box):
            a.set_transform(trans)
        return [marked, box]


# A rotated road line with a per-domain LANDWARD displacement applied
def relocate_line(frame, rotated, displacement_m):
    if not displacement_m:
        return rotated
    sea = frame.seaward
    cen = frame.pts(frame.centroids)
    xs, ids = cen[:, 0], np.arange(FIRST, FIRST + len(cen), dtype=float)
    if xs[0] > xs[-1]:                       # np.interp needs increasing x
        xs, ids = xs[::-1], ids[::-1]

    def moved(coords):
        a = np.asarray(coords, float)[:, :2]
        gis = np.rint(np.interp(a[:, 0], xs, ids)).astype(int)
        d = np.array([displacement_m.get(int(g), 0.0) for g in gis])
        return np.c_[a[:, 0] - sea[0] * d, a[:, 1] - sea[1] * d]

    def apply(geom):
        if geom.geom_type == "LineString":
            return LineString(moved(geom.coords))
        if geom.geom_type == "MultiLineString":
            return MultiLineString([moved(g.coords) for g in geom.geoms])
        return geom

    return rotated.apply(apply)


# Tint the seaward or soundward HALF of each selected domain box
def tint_domains(ax, frame, rdom, sel, half, **kw):
    whole = half == "all"
    lower = whole or half_is_lower(frame, half)
    for geom in rdom[sel]:
        x0, y0, x1, y1 = geom.bounds
        mid = y1 if whole else (y0 + y1) / 2
        lo, hi = (y0, mid) if lower else (mid, y1)
        ax.add_patch(mpl.patches.Rectangle((x0, lo), x1 - x0, hi - lo, **kw))


# Half-width and half-height of each string, in DATA metres
def measure_m(ax, texts, fontsize, **kw):
    renderer = ax.figure.canvas.get_renderer()
    p0 = ax.transData.transform((0.0, 0.0))
    p1 = ax.transData.transform((1000.0, 0.0))
    m_per_px = 1000.0 / abs(p1[0] - p0[0])
    out = []
    for text in texts:
        probe = ax.text(0, 0, text, fontsize=fontsize, ha="center", va="center", **kw)
        bb = probe.get_window_extent(renderer=renderer)
        probe.remove()
        out.append((bb.width * m_per_px / 2, bb.height * m_per_px / 2))
    return out


# Push labels apart along the reach until none overlaps, in place
def separate_x(ax, placed, pad_m):
    xlo, xhi = ax.get_xlim()
    margin = pad_m / 2
    placed.sort(key=lambda d: d["x"])
    for i in range(1, len(placed)):
        need = placed[i - 1]["x"] + placed[i - 1]["hw"] + placed[i]["hw"] + pad_m
        placed[i]["x"] = max(placed[i]["x"], need)
    for i in range(len(placed) - 1, -1, -1):
        placed[i]["x"] = min(placed[i]["x"], xhi - margin - placed[i]["hw"])
        if i:
            room = placed[i]["x"] - placed[i]["hw"] - pad_m - placed[i - 1]["hw"]
            placed[i - 1]["x"] = min(placed[i - 1]["x"], room)
    for d in placed:
        d["x"] = max(d["x"], xlo + margin + d["hw"])
    return placed


# Annotation labels beside the reach, on as few straight rows as they fit
def label_lanes(fig, ax, frame, items, side, y_edge, first_offset_m,
                row_step_m=None, pad_m=1100.0, fontsize=8, leader_gap_m=600.0,
                clear_m=0.0, clear_pad_m=700.0):
    if not items:
        return y_edge
    sea = frame.seaward
    out = side * sea[1]                     # +1/-1 in y, away from the reach
    sized = measure_m(ax, [t for _, _, t in items], fontsize)
    placed = [dict(x=float(frame.along(gis_mid)[0][0]), anchor=anchor, text=text,
                   hw=hw, hh=hh)
              for (gis_mid, anchor, text), (hw, hh) in zip(items, sized)]

    # Row assignment for the labels (README)
    xlo, xhi = ax.get_xlim()
    placed.sort(key=lambda d: d["x"])
    avail = (xhi - xlo) - 2 * pad_m
    rows, row, used = [], [], 0.0
    for d in placed:
        need = 2 * d["hw"] + (pad_m if row else 0.0)
        if row and used + need > avail:
            rows.append(row)
            row, used, need = [], 0.0, 2 * d["hw"]
        row.append(d)
        used += need
    if row:
        rows.append(row)

    for row in rows:
        separate_x(ax, row, pad_m)

    # One row's worth of clear space is the tallest label plus half a line of air, so two rows can (README)
    if row_step_m is None:
        row_step_m = max(d["hh"] for d in placed) * 2.9

    # The first row clears both the offset asked for and anything already occupying the space (the (README)
    first = max(first_offset_m,
                clear_m + clear_pad_m + max(d["hh"] for d in placed))

    y_out = y_edge
    for ri, row in enumerate(rows):
        y_row = y_edge + out * (first + ri * row_step_m)
        for d in row:
            ax.text(d["x"], y_row, d["text"], ha="center", va="center",
                    fontsize=fontsize, color=INK, zorder=9, path_effects=HALO)
            # attach to the NEAREST part of the block, not its midpoint
            a = d["anchor"]
            ax_lo, ax_hi, ay = (a if len(a) == 3 else (a[0], a[0], a[1]))
            attach = min(max(d["x"], ax_lo), ax_hi)
            gap = math.hypot(d["x"] - attach, y_row - ay)
            if gap > leader_gap_m:
                ax.plot([d["x"], attach],
                        [y_row - out * (d["hh"] + 150), ay + out * 260],
                        color=INK_MUTED, lw=0.5, zorder=7)
            edge = y_row + out * d["hh"]
            y_out = max(y_out, edge) if out > 0 else min(y_out, edge)
    return y_out


# Two panels on one reach
def fig_management_footprint(dom, outline, roads, frame):
    # One window for both panels
    window = frame.window(dom, pad_along_km=3.0, pad_sea_km=5.2, pad_sound_km=4.4)
    x0, x1, y0, y1 = window
    aspect = (x1 - x0) / (y1 - y0)

    LEFT, PW = 0.004, 0.992
    panel_h_in = FIG_W_DOUBLE * PW / aspect
    GAP_IN, LEG_IN = 0.10, 0.58          # between panels; the legend strip
    fig_h = 2 * panel_h_in + GAP_IN + LEG_IN
    fig = plt.figure(figsize=(FIG_W_DOUBLE, fig_h), facecolor="white")
    ax_a = fig.add_axes([LEFT, (LEG_IN + GAP_IN + panel_h_in) / fig_h,
                         PW, panel_h_in / fig_h])
    ax_b = fig.add_axes([LEFT, LEG_IN / fig_h, PW, panel_h_in / fig_h])

    sea = frame.seaward
    rdom = frame.geoms(dom.geometry)
    bx0, by0, bx1, by1 = rdom.total_bounds
    sea_edge = by0 if sea[1] < 0 else by1

    # The period starts this figure speaks to, and the dune line each reads (README)
    PERIOD_STARTS = (1996, 2010)
    early_start, late_start = PERIOD_STARTS
    early_line = tv.road_line_for_year(early_start)
    late_line = tv.road_line_for_year(late_start)

    # Tint opaque enough to read over land and water, light enough to keep the edges
    TINT_FILL = 0.45
    TINT_ROAD = 0.30

    # Where a label's leader stops
    def anchor_at(gis):
        sel = dom.ID.isin(gis).values
        _x0, _y0, _x1, _y1 = rdom[sel].total_bounds
        return (float(_x0), float(_x1), _y0 if sea[1] < 0 else _y1)

    # The shared basemap
    def panel(ax, i, name, scalebar, show_road):
        draw_reach(ax, frame, dom, outline, roads[2008], vector=True, window=window,
                   label_villages=False, water_labels=False, piers=False, numbers=False,
                   ends=False, groin=False, scalebar=False, arrow=False, show_road=show_road)
        study_area_labels(ax, frame, dom, True, outline, text_kw=MAP_TEXT_DARK, leader_c=INK_MUTED,
                          shift_m={}, water=False, groin_label=False)
        if scalebar:
            scale_bar_km(ax, text_kw=MAP_TEXT_DARK)
            north_dart(ax, axes_point(ax, 0.47, 0.075), north=frame.north(), text_kw=MAP_TEXT_DARK)
        # Upper right (2026-09-17)
        RIGHT_X, GAP = 0.985, 0.012
        probe = ax.text(0, 0, name, fontsize=9)
        w = (probe.get_window_extent(renderer=fig.canvas.get_renderer()).width
             / ax.get_window_extent().width)
        probe.remove()
        ax.text(RIGHT_X, 0.94, name, transform=ax.transAxes, ha="right", va="top",
                fontsize=9, color=INK, zorder=20,
                bbox=dict(facecolor="white", alpha=0.85, edgecolor="none",
                          boxstyle="square,pad=0.2"))
        letter_corner(ax, i, x=RIGHT_X - w - GAP, y=0.94, ha="right")

    # (a) Beach nourishment; NC-12 is not drawn here (README)
    panel(ax_a, 0, "beach nourishment", scalebar=False, show_road=False)
    fills = []
    for p_ in sorted(HATTERAS_NOURISHMENT_PROJECTS, key=lambda q: q.year):
        if not p_.enabled:
            continue
        gis = sorted(p_.gis_domains)
        sel = dom.ID.isin(gis).values
        tint_domains(ax_a, frame, rdom, sel, "all", facecolor=C["ADDED"],
                     edgecolor="none", alpha=TINT_FILL, zorder=3)
        # Terse labels: width is what makes them collide
        fills.append(((gis[0] + gis[-1]) / 2, anchor_at(gis),
                      f"{p_.name.split()[0]} fill, {p_.year}\n"
                      f"GIS {gis[0]}\u2013{gis[-1]}, "
                      f"{p_.volume_cubic_yards / 1e6:.1f} M yd\u00b3"))
    label_lanes(fig, ax_a, frame, fills, +1, sea_edge, first_offset_m=2500)

    # (b) NC-12 ROADWAY
    panel(ax_b, 1, "NC-12 roadway", scalebar=True, show_road=False)
    # Both vintages, the earlier in red, so a relocation shows as the gap (README)
    before_start = {}
    for _e in HATTERAS_ROAD_EVENTS:
        if hasattr(_e, "displacement_m") and _e.year <= early_start:
            before_start.update(_e.displacement_m)
    early_geom = relocate_line(frame, frame.geoms(roads[early_line].geometry),
                               before_start)
    early_geom.plot(ax=ax_b, color=C_1984, lw=0.9, zorder=5)
    frame.geoms(roads[late_line].geometry).plot(ax=ax_b, color=C_1997, lw=0.9, zorder=6)
    # The bridge span CONTAINS the 1989 relocation, so it is a hatch laid over the tint rather than (README)
    events = []
    for e in HATTERAS_ROAD_EVENTS:
        if hasattr(e, "displacement_m"):
            gis = sorted(e.displacement_m)
            sel = dom.ID.isin(gis).values
            tint_domains(ax_b, frame, rdom, sel, "all", facecolor=C["ROAD"],
                         edgecolor="none", alpha=TINT_ROAD, zorder=3)
            d = e.displacement_m
            events.append(((gis[0] + gis[-1]) / 2, anchor_at(gis),
                           f"relocated {e.year}\n"
                           f"GIS {gis[0]}\u2013{gis[-1]}, "
                           f"{min(d.values()):.0f}\u2013{max(d.values()):.0f} m"))
        else:
            gis = sorted(e.gis_domains)
            sel = dom.ID.isin(gis).values
            tint_domains(ax_b, frame, rdom, sel, "all", facecolor="none",
                         edgecolor=INK, hatch="////", lw=0.0, zorder=6)
            events.append(((gis[0] + gis[-1]) / 2, anchor_at(gis),
                           f"bridge {e.year}\n"
                           f"road off GIS {gis[0]}\u2013{gis[-1]}"))
    label_lanes(fig, ax_b, frame, events, +1, sea_edge, first_offset_m=2500)

    # One legend for both panels; each swatch is a domain box marked as the panels mark one
    handles = [Line2D([], [], color=C_1984, lw=1.2,
                      label=f"(b) NC-12, {early_start} alignment"),
               Line2D([], [], color=C_1997, lw=1.2,
                      label=f"(b) NC-12, {late_start} alignment"),
               HalfBox(C["ADDED"], True, TINT_FILL, whole=True,
                       label="(a) domains nourished"),
               HalfBox(C["ROAD"], True, TINT_ROAD, whole=True,
                       label="(b) road relocated"),
               HalfBox(None, True, TINT_ROAD, hatch="////", whole=True,
                       label="(b) road removed after the bridge"),
               groin_handle()]
    fig.legend(handles=handles, loc="lower center", bbox_to_anchor=(0.5, 0.004),
               ncol=3, frameon=False, fontsize=8, handlelength=1.3,
               handleheight=1.25, columnspacing=1.6, labelspacing=0.6,
               handler_map={HalfBox: HalfBoxHandler()})
    out = save(fig, fig_path("management_footprint"), vector=True, dpi=300)
    record_caption(out[0],
        "Where the reach has been managed, as the model prescribes it, on one reach in two "
        "panels. (a) The three beach-nourishment projects in the record (Rodanthe 2014, Buxton "
        "and Avon 2022), each tinted across the domains it was spread over, with the reported "
        "project total; NC-12 is not drawn in (a), though it is why the Buxton and Rodanthe "
        "fills run past their villages, along the road corridor. (b) NC-12: the "
        f"{early_start} (red) and {late_start} (blue) alignments, digitised from the 1978 and "
        "2008 imagery, each carrying the relocations that precede its period start \u2014 so the "
        f"red line has the 1989 Pea Island move applied, stepped per domain the way the model "
        "applies it. The two therefore coincide at GIS 84\u201387, where the road did not move "
        "between the starts, and separate "
        "at GIS 9\u201314, the 1999 relocation south of Avon, which is the one that falls between "
        "them. Both events are tinted and labelled with the measured displacement "
        "range; and hatched where the Jug Handle Bridge carries the road off the barrier from "
        "2022 (GIS 82\u201388). Every mark is drawn to WHOLE domains, which is the resolution the "
        "model applies them at, not the placement geometry. The two panels share one window, so "
        "a domain is the same place on both: the Rodanthe fill and the bridge cover the same "
        "ground. South to north from left to right; the Buxton groin field is the red bar at "
        "GIS 5\u20136, and the scale bar and north arrow in (b) serve both panels.")
    plt.close(fig)
    return out[0]


# Figure 9: the reach at 10 m, both products

# Domains lo..hi of one extraction as one mosaic (m NAVD88, water -10), ocean at the top, cropped
def reach_mosaic(product, dom, lo, hi, margin_cells=8):
    d, _ = tv.npy_dirs(product)
    boxes = dom[(dom.ID >= lo) & (dom.ID <= hi)]
    n = hi - lo + 1
    lefts = np.array([geom.bounds[0] for geom in boxes.geometry])
    e0 = lefts.min()
    ncol_all = int(round((lefts.max() - e0) / CELL_M)) + 200
    canvas = np.full((ncol_all, n * 50), -10.0)
    for k, g in enumerate(range(lo, hi + 1)):
        a = np.load(d / f"domain_{g}.npy")[::-1]           # rows north to south -> south to north
        c0 = int(round((lefts[k] - e0) / CELL_M))
        canvas[c0:c0 + a.shape[1], k * 50:(k + 1) * 50] = a.T
    # The shear
    jc = np.arange(n) * 50 + 25
    slope, icpt = np.polyfit(jc, lefts, 1)
    js = np.arange(n * 50)
    shift = np.round((slope * js + icpt - (slope * jc.mean() + icpt)) / CELL_M).astype(int)
    sheared = np.full_like(canvas, -10.0)
    for jj, sh in enumerate(shift):
        col = canvas[:, jj]
        if sh > 0:
            sheared[:-sh, jj] = col[sh:]
        elif sh < 0:
            sheared[-sh:, jj] = col[:sh]
        else:
            sheared[:, jj] = col
    land_rows = np.where((sheared > 0).any(1))[0]
    r0 = max(land_rows.min() - margin_cells, 0)
    r1 = min(land_rows.max() + margin_cells + 1, ncol_all)
    crop = sheared[r0:r1][::-1]                            # ocean (high easting) at the top
    return crop, (e0 + r1 * CELL_M, e0 + r0 * CELL_M)


# The whole reach's elevation, in three stacked sections
def fig_reach_elevation(dom):
    groups = ((1, 30), (31, 60), (61, 90))
    products = ("1984-start", "2004-start")
    cmap, norm, bounds = elevation_cmap()
    mosaics = {(g, p): reach_mosaic(p, dom, *g) for g in groups for p in products}
    # every strip of a group shares one cross-shore window: the union of its two products
    rows = {g: max(mosaics[(g, p)][0].shape[0] for p in products) for g in groups}
    strip_in = {g: FIG_W_DOUBLE * 0.905 / (30 * 50 / rows[g]) for g in groups}
    gap_in, group_gap_in, top_in, bottom_in = 0.04, 0.30, 0.12, 1.05
    fh = top_in + sum(2 * strip_in[g] + gap_in for g in groups) + 2 * group_gap_in + bottom_in
    fig = plt.figure(figsize=figsize("double", height=fh))
    y = fh - top_in
    axes = []
    for gi, (lo, hi) in enumerate(groups):
        for pi, product in enumerate(products):
            y -= strip_in[(lo, hi)]
            ax = fig.add_axes([0.075, y / fh, 0.905, strip_in[(lo, hi)] / fh])
            axes.append(ax)
            strip, (top, bottom) = mosaics[((lo, hi), product)]
            im = ax.imshow(strip, cmap=cmap, norm=norm, extent=(lo - 0.5, hi + 0.5, bottom, top),
                           interpolation="nearest", aspect="auto", zorder=1)
            ax.set_facecolor(ELEV_WATER)
            ax.set_ylim(bottom, bottom + rows[(lo, hi)] * CELL_M)
            ax.set_yticks([])
            ax.set_xlim(lo - 0.5, hi + 0.5)
            spines_for_image(ax)
            ax.text(0.055 if pi == 0 else 0.004, 0.94, product.replace("-start", " start"), transform=ax.transAxes,
                    ha="left", va="top", fontsize=8, color=INK, zorder=5,
                    bbox=dict(facecolor="white", alpha=0.85, edgecolor="none", boxstyle="square,pad=0.15"))
            ax.set_xticks(range(((lo + 4) // 5) * 5, hi + 1, 5))
            if pi == 0:
                ax.tick_params(labelbottom=False, bottom=False)
                letter_corner(ax, gi)
                y -= gap_in
            else:
                ax.set_xlabel(DOMAIN_AXIS_LABEL if gi == 2 else "")
        y -= group_gap_in
    for ax in axes[1::2]:
        for name, (lo_, hi_) in ANN.town_spans.items():
            ax.add_patch(Rectangle((lo_ - 0.5, 0), hi_ - lo_ + 1, 0.07, transform=ax.get_xaxis_transform(),
                                   facecolor="0.94", edgecolor="0.6", lw=0.4, zorder=5))
    cax = fig.add_axes([0.35, 0.42 / fh, 0.35, 0.10 / fh])
    cb = fig.colorbar(im, cax=cax, orientation="horizontal", ticks=bounds[1:-1])
    cb.set_label("elevation (m NAVD88)", fontsize=8)
    cb.ax.tick_params(labelsize=8)
    cb.outline.set_linewidth(0.5)
    out = save(fig, fig_path("reach_elevation"), vector=False, dpi=300)
    record_caption(out[0],
        "The reach at the model's 10 m resolution, before the dune search: the per-domain elevation "
        "arrays of the 1984-start extraction (the 2009–2014 lidar with the 1996 ALACE survey grafted on "
        "the ocean side, upper strip of each pair) and the 2004-start extraction (the 2009–2014 lidar "
        "alone, lower strip), in the elevation classes of the house style, m NAVD88, water below 0 m. "
        "Each domain is 50 cells (500 m) alongshore by 200 cells (2000 m) cross-shore. The domains are "
        "placed at their eastings and the strip is then sheared by the linear alongshore trend of the box "
        "eastings, the straightening the extractor itself applies, so the coast is continuous and runs level "
        "while its residual bend remains; cross-shore distances within a column are unchanged, and the strips "
        "are cropped to the rows that hold land. "
        "Ocean at the top; (a) GIS 1–30, (b) 31–60, (c) 61–90, south to north from left to right, at equal "
        "scale within a strip. The light bands under the lower strips are the village spans.")
    plt.close(fig)
    return out[0]


# Projector settings: heavier lines and fills for a slide seen from afar
def talk_mode():
    global WATER_MAP, TALK
    WATER_MAP = "#f3f6f9"
    TALK = True
    for key in ("font.size", "axes.labelsize", "xtick.labelsize", "ytick.labelsize", "legend.fontsize",
                "axes.titlesize"):
        if key in mpl.rcParams:
            try:
                mpl.rcParams[key] = float(mpl.rcParams[key]) + 1.0
            except (TypeError, ValueError):
                pass
    mpl.rcParams["lines.linewidth"] = float(mpl.rcParams["lines.linewidth"]) * 1.3
    mpl.rcParams["axes.linewidth"] = float(mpl.rcParams["axes.linewidth"]) * 1.3


# name -> a callable taking the shared layers by keyword
FIGURES = {
    "study_area": lambda dom, outline, roads, frame, vector: fig_study_area(dom, outline, roads, frame, vector),
    "domain_framework": lambda dom, outline, roads, frame, vector: fig_domain_framework(dom, outline, roads, frame, vector),
    "domain_metrics": lambda **k: fig_domain_metrics(),
    "domain_framework_vertical": lambda dom, outline, roads, frame, vector: fig_domain_framework_vertical(dom, outline, roads, vector),
    "site_overview": lambda dom, outline, roads, frame, vector: fig_site_overview(dom, outline, roads, frame, vector),
    "domain_grid": lambda dom, outline, roads, frame, vector: fig_domain_grid(dom, roads, vector),
    "forcing_timeline_1984": lambda **k: fig_forcing_timeline("forcing_timeline_1984"),
    "forcing_timeline_1996": lambda **k: fig_forcing_timeline("forcing_timeline_1996"),
    "observed_rates": lambda **k: fig_observed_rates(),
    "dune_lines": lambda dom, outline, roads, frame, vector: fig_dune_lines(dom, roads, vector),
    "domain_schematic": lambda **k: fig_domain_schematic(),
    "management_footprint": lambda dom, outline, roads, frame, vector: fig_management_footprint(dom, outline, roads, frame),
    "reach_elevation": lambda dom, **k: fig_reach_elevation(dom),
}


# Run: draw every figure, or the one named with --only
def main():
    ap = argparse.ArgumentParser(description="The generic Hatteras site figures.")
    ap.add_argument("--vector", action="store_true", help="island outline instead of imagery (no network)")
    ap.add_argument("--only", choices=list(FIGURES), default=None)
    ap.add_argument("--talk", action="store_true",
                    help="projector version under talk/: water barely off-white, type one point "
                         "larger, lines heavier")
    args = ap.parse_args()
    apply_style()
    if args.talk:
        talk_mode()
    FIG_ROOT.mkdir(parents=True, exist_ok=True)
    dom, outline, roads = load_layers()
    frame = Frame(dom, outline)
    print(f"reach axis {math.degrees(frame.theta):.1f} deg from east; {frame.step:.0f} m between centroids; "
          f"ocean at the {'bottom' if frame.seaward[1] < 0 else 'top'}")
    todo = [args.only] if args.only else list(FIGURES)
    for name in todo:
        print("wrote", FIGURES[name](dom=dom, outline=outline, roads=roads, frame=frame, vector=args.vector))


if __name__ == "__main__":
    main()
