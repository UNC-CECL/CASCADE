"""
Maps of the full-record CoastSat LRR, Cape Point to the Virginia line: an overview and four regional zooms.

    python scripts/input_prep/5-scr/3-rates/coastsat/lrr/coastsat_obx_lrr_maps.py

Reads the table coastsat_obx_lrr.py wrote; fits nothing. Details: scripts/input_prep/5-scr/3-rates/README.md.

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-09-27
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import geopandas as gpd
from shapely.geometry import LineString

_REPO = next(p for p in Path(__file__).resolve().parents
             if (p / "pyproject.toml").exists())
sys.path.insert(0, str(_REPO / "scripts"))
from site_layer import hat_observed_rates as obs  # noqa: E402
from site_layer.hat_figure_style import (  # noqa: E402
    apply_style, figsize, open_frame, save, caption, spines_for_image,
    _scalebar, _title, INK, INK_MUTED, GRID_C)

import matplotlib.pyplot as plt  # noqa: E402
import matplotlib.colors as mcolors  # noqa: E402
import contextily as ctx  # noqa: E402

# --- CONFIG ------------------------------------------------------------------
TAG = "1984_2025"
DIR = obs.COASTSAT_LRR_ROOT / "{0}_obx".format(TAG)
TABLE = DIR / "supporting" / "coastsat_lrr_obx_{0}_full.csv".format(TAG)   # needs the seaward ends
UTM = 32618
# Label-free relief basemap: every name on the maps is placed here
BASEMAP = ctx.providers.Esri.WorldShadedRelief

V_HALF = 3.0          # m/yr, colour scale; beyond it saturates (arrows on the bar)
X_HALF = 4.0          # m/yr, profile axis; beyond it a triangle at the edge
OVERLAP_KM = 1.0

# Latitudes of the places named on the maps and used to split the regions
from coastsat_obx_lrr import PLACES  # noqa: E402  (one list for both scripts)
# Stretches rather than towns, named on the maps only (italic)
MAP_AREAS = {"Pea Island": 35.700, "Bodie Island spit": 35.800}
# (letter, stem, title, south place or None = start, north place or None = end)
REGIONS = [
    ("A", "cape_point_to_salvo", "Cape Point to Salvo", None, 35.560),
    # B runs ~15 km past Oregon Inlet, so not "to Oregon Inlet"
    ("B", "rodanthe_to_south_nags_head", "Rodanthe to South Nags Head", 35.560, 35.900),
    ("C", "nags_head_to_duck", "Nags Head to Duck", 35.900, 36.200),
    ("D", "corolla_to_virginia", "Corolla to the Virginia line", 36.200, None),
]
# -----------------------------------------------------------------------------


# The table, with each transect as a line
def load():
    t = pd.read_csv(TABLE, keep_default_na=False, na_values=[""])
    t["flag"] = t["flag"].fillna("")
    lines = [LineString([(a, b), (c, d)]) for a, b, c, d in
             t[["origin_lon", "origin_lat", "seaward_lon", "seaward_lat"]].to_numpy()]
    g = gpd.GeoDataFrame(t, geometry=lines, crs=4326).to_crs(UTM)
    o = np.array([ls.coords[0] for ls in g.geometry])
    g["x0"], g["y0"] = o[:, 0], o[:, 1]
    g["ymid"] = [0.5 * (ls.coords[0][1] + ls.coords[1][1]) for ls in g.geometry]
    return g


# Alongshore km at a latitude
def km_at_lat(g, lat):
    return float(g["alongshore_km"].iloc[(g["origin_lat"] - lat).abs().argmin()])


# A 1 km running median along one side
def running_median(x, y, side):
    med = np.full_like(y, np.nan, dtype=float)
    for i in range(len(x)):
        w = (np.abs(x - x[i]) <= 0.5) & (side == side[i]) & np.isfinite(y)
        if w.sum() >= 5:
            med[i] = np.median(y[w])
    return med


# Drawing

NORM = mcolors.Normalize(-V_HALF, V_HALF)
CMAP = plt.get_cmap("RdBu")      # red erosion (negative), blue accretion


# Transects coloured by rate, the basemap, places, scale and arrow
def draw_map(ax, g, lw, places, state_y, label_pt=7):
    for (x0, y0), (x1, y1), v in zip(
            [ls.coords[0] for ls in g.geometry],
            [ls.coords[1] for ls in g.geometry], g["lrr_m_yr"]):
        ax.plot([x0, x1], [y0, y1], color=CMAP(NORM(np.clip(v, -V_HALF, V_HALF))),
                lw=lw, solid_capstyle="butt", zorder=3)
    # Equal aspect by widening x to fill the panel ("datalim")
    ax.set_aspect("equal", adjustable="datalim")
    ax.figure.canvas.draw()
    ax.set_xlim(*ax.get_xlim())
    ax.set_ylim(*ax.get_ylim())
    ctx.add_basemap(ax, source=BASEMAP, crs=f"EPSG:{UTM}", attribution_size=4,
                    zorder=0)
    greyscale_basemap(ax)
    y_lo, y_hi = ax.get_ylim()
    for name, lat in places.items():
        i = (g["origin_lat"] - lat).abs().argmin()
        if abs(g["origin_lat"].iloc[i] - lat) > 0.01:
            continue                        # the place is outside this panel
        x, y = g["x0"].iloc[i], g["y0"].iloc[i]
        if not (y_lo < y < y_hi):
            continue
        ax.plot([x - 250, x - 900], [y, y], color=INK_MUTED, lw=0.5, zorder=4)
        ax.text(x - 1000, y, name, ha="right", va="center", fontsize=label_pt,
                color=INK, style="italic" if (name == "Oregon Inlet" or name in MAP_AREAS)
                else "normal", zorder=5)
    # The NC/VA state line, where the product ends: drawn across the panel.
    if y_lo < state_y < y_hi:
        ax.axhline(state_y, color=INK, lw=0.6, ls=(0, (4, 2)), zorder=4)
        x_l = ax.get_xlim()[0] + 0.03 * (ax.get_xlim()[1] - ax.get_xlim()[0])
        # Above the line, clear of box D's edge
        ax.text(x_l, state_y, "North Carolina / Virginia state line", ha="left", va="bottom",
                fontsize=label_pt, color=INK, zorder=5,
                bbox=dict(fc="white", ec="none", pad=1.0, alpha=0.8))
    spines_for_image(ax)
    ax.set_xticks([])
    latitude_ticks(ax, g)


# Recolour the basemap tiles to light grey
def greyscale_basemap(ax, lo=0.80, hi=0.96):
    im = ax.images[-1]
    a = np.asarray(im.get_array(), dtype=float)
    rgb = a[..., :3] / (255.0 if a.max() > 1 else 1.0)
    lum = rgb @ np.array([0.299, 0.587, 0.114])
    q1, q9 = np.percentile(lum, [2, 98])
    g = lo + (hi - lo) * np.clip((lum - q1) / max(q9 - q1, 1e-6), 0, 1)
    out = np.dstack([g, g, g] + ([a[..., 3] / (255.0 if a.max() > 1 else 1.0)]
                                 if a.shape[-1] == 4 else []))
    im.set_data(out)


# Latitude on the map's left edge (the reviewer
def latitude_ticks(ax, g):
    from pyproj import Transformer
    to_utm = Transformer.from_crs(4326, UTM, always_xy=True)
    lon = float(g["origin_lon"].median())
    y_lo, y_hi = ax.get_ylim()
    span = g["origin_lat"].max() - g["origin_lat"].min()
    step = 0.2 if span > 0.8 else 0.05
    lats = np.arange(np.floor(g["origin_lat"].min() / step) * step - step,
                     g["origin_lat"].max() + 2 * step, step)
    ys = np.array([to_utm.transform(lon, la)[1] for la in lats])
    keep = (ys > y_lo) & (ys < y_hi)
    fmt = "{0:.1f}°N" if step >= 0.1 else "{0:.2f}°N"
    ax.set_yticks(ys[keep], labels=[fmt.format(la) for la in lats[keep]])
    ax.tick_params(axis="y", left=True, labelleft=True, length=3, width=0.5,
                   labelsize=7, colors=INK)


# The rate profile beside a map
def draw_profile(ax, g, inlet_km):
    x = g["alongshore_km"].to_numpy()
    y = g["lrr_m_yr"].to_numpy()
    n = g["ymid"].to_numpy()
    side = x > inlet_km
    med = running_median(x, y, side)
    ax.axvline(0, color=INK, lw=0.6, zorder=2)
    # Every transect drawn alike; flags stay in the table
    ax.scatter(np.clip(y, -X_HALF, X_HALF), n, s=14, edgecolors="0.55", linewidths=0.25,
               c=CMAP(NORM(np.clip(np.nan_to_num(y), -V_HALF, V_HALF))), zorder=3)
    for s_ in (~side, side):
        ax.plot(np.clip(med[s_], -X_HALF, X_HALF), n[s_], color=INK, lw=0.7, zorder=4)
    out = np.isfinite(y) & (np.abs(y) > X_HALF)
    if out.any():
        # Small and in the scale's end colour
        for sign, marker, colour in ((-1, "<", CMAP(0.0)), (1, ">", CMAP(1.0))):
            s_ = out & (np.sign(y) == sign)
            ax.scatter(np.full(s_.sum(), sign * X_HALF), n[s_], marker=marker,
                       s=14, lw=0, c=[colour], zorder=5, clip_on=False)
    ax.set_xlim(-X_HALF, X_HALF)
    ax.set_xlabel("LRR (m/yr)")
    ax.grid(True, axis="x", color=GRID_C, lw=0.4)
    ax.tick_params(axis="y", left=False, labelleft=False)
    open_frame(ax)
    ax.spines["left"].set_visible(False)
    off = sorted(zip(x[out], y[out]), key=lambda t: -abs(t[1]))
    return off


# A cartographic north arrow
def north_arrow(ax, inset_in=0.12, height_in=0.46, width_in=0.18, n_pt=11):
    from matplotlib.patches import Polygon, Rectangle
    from matplotlib.transforms import ScaledTranslation
    fig = ax.figure
    tr = fig.dpi_scale_trans + ScaledTranslation(1.0, 1.0, ax.transAxes)
    h, w, notch = height_in, width_in / 2, 0.28 * height_in
    n_h = n_pt / 72.0                        # height of the N, inches
    pad, gap = 0.08, 0.04
    box_w = 2 * w + 2 * pad
    box_h = pad + h + gap + n_h + pad
    bx, by = -inset_in - box_w, -inset_in - box_h      # box lower-left
    x0, y0 = bx + box_w / 2, by + pad                   # arrow base centre
    ax.add_patch(Rectangle((bx, by), box_w, box_h, transform=tr, fc="white",
                           ec="black", lw=0.5, zorder=8, clip_on=False))
    ax.add_patch(Polygon([(x0, y0 + h), (x0 - w, y0), (x0, y0 + notch), (x0 + w, y0)],
                         closed=True, transform=tr, fc="black", ec="black",
                         lw=0.6, zorder=9, clip_on=False))
    ax.text(x0, y0 + h + gap, "N", transform=tr, ha="center", va="bottom",
            fontsize=n_pt, fontweight="bold", color="black", zorder=9)


# The shared horizontal colourbar
def colourbar(fig, axes):
    sm = plt.cm.ScalarMappable(norm=NORM, cmap=CMAP)
    cb = fig.colorbar(sm, ax=axes, orientation="horizontal", location="bottom",
                      extend="both", shrink=0.6, aspect=35, pad=0.02)
    cb.set_label("Shoreline change rate, LRR (m/yr)")
    cb.outline.set_linewidth(0.5)
    return cb


PROFILE_W = 2.6       # in, panel (b)
CHROME_W = 0.5        # in, margins between and around the panels
CHROME_H = 1.6        # in, titles, x labels and colour bar


# One region: map and profile, sized from the region's shape
def figure(g, inlet_km, title_a, lw, places, state_y, boxes=None, height=9.0):
    # Size the map panel from the region's own shape
    pad = 0.02 * (g["y0"].max() - g["y0"].min())
    y_lo, y_hi = g["y0"].min() - pad, g["y0"].max() + pad
    ends = np.array([ls.coords[1][0] for ls in g.geometry])
    x_lo, x_hi = g["x0"].min() - 8000, max(g["x0"].max(), ends.max()) + 3000
    ratio = (x_hi - x_lo) / (y_hi - y_lo)                 # map width / height
    fig_w = figsize("double")[0]
    map_w = fig_w - PROFILE_W - CHROME_W
    map_h = min(map_w / ratio, height - CHROME_H)
    map_w = map_h * ratio
    fig, (am, ap) = plt.subplots(
        1, 2, sharey=True, figsize=figsize("double", height=map_h + CHROME_H),
        gridspec_kw={"width_ratios": [map_w, PROFILE_W]}, constrained_layout=True)
    am.set_ylim(y_lo, y_hi)
    am.set_xlim(x_lo, x_hi)
    draw_map(am, g, lw, places, state_y)
    if boxes:
        for letter, (x_lo, x_hi, y_lo, y_hi) in boxes.items():
            am.add_patch(plt.Rectangle((x_lo, y_lo), x_hi - x_lo, y_hi - y_lo,
                                       fill=False, ec=INK, lw=0.7, zorder=6))
            am.text(x_hi + 600, 0.5 * (y_lo + y_hi), letter, fontsize=10,
                    fontweight="bold", color=INK, va="center", zorder=6)
    span = am.get_xlim()[1] - am.get_xlim()[0]
    _scalebar(am, length_m=10_000 if span > 25_000 else 5_000)
    north_arrow(am)                    # top right: open water in every panel
    off = draw_profile(ap, g, inlet_km)
    _title(am, 0, title_a)
    _title(ap, 1, "Rate by transect")
    colourbar(fig, [am, ap])
    return fig, off


# The caption clause for rates off the scale
def off_clause(off):
    if not off:
        return ""
    return (" Beyond ±{0:g} m/yr in (b), marked with a triangle at the edge: "
            .format(X_HALF)
            + "; ".join("{0:+.1f} m/yr at {1:.1f} km".format(v, k) for k, v in off[:6])
            + ("; and {0} more".format(len(off) - 6) if len(off) > 6 else "") + ".")


COMMON = (
    "Full-record CoastSat shoreline change rate, 1 January 1984 to 31 December "
    "2025: the ordinary least-squares slope of every tidally corrected shoreline "
    "position on each CoastSat transect, no outlier filter, no weighting; "
    "positive is seaward (accretion). (a) Each transect drawn at its position "
    "and length, coloured on a fixed scale that saturates at ±{v:g} m/yr. "
    "(b) The same rates against northing, on (a)'s vertical axis; black line, "
    "the median within 0.5 km of coast either side, computed separately north "
    "and south of Oregon Inlet. Latitude on the left edge is exact at the "
    "coast. Basemap Esri World Shaded Relief (no labels), shown in greyscale; every "
    "name is placed from the town's latitude. Table and method: coastsat_lrr_obx_1984_2025.csv and README.md "
    "in this folder.").format(v=V_HALF)


# Run: the overview and every region
def main():
    apply_style()
    g = load()
    gaps = np.diff(g["alongshore_km"].to_numpy())
    i = int(np.argmax(gaps))
    inlet_km = 0.5 * (g["alongshore_km"].iloc[i] + g["alongshore_km"].iloc[i + 1])

    # region extents in km, from the split latitudes
    reg = []
    for letter, stem, title, lat_lo, lat_hi in REGIONS:
        k_lo = 0.0 if lat_lo is None else km_at_lat(g, lat_lo) - OVERLAP_KM / 2
        k_hi = g["alongshore_km"].max() if lat_hi is None else km_at_lat(g, lat_hi) + OVERLAP_KM / 2
        sub = g[(g["alongshore_km"] >= k_lo) & (g["alongshore_km"] <= k_hi)]
        reg.append((letter, stem, title, sub))

    boxes = {}
    for letter, _, _, sub in reg:
        boxes[letter] = (sub["x0"].min() - 1500, sub["x0"].max() + 1500,
                         sub["y0"].min(), sub["y0"].max())

    state_y = float(g["y0"].iloc[-1])      # usa_NC_0049_0230 ends at the line
    over_places = {k: PLACES[k] for k in ("Cape Point", "Buxton", "Rodanthe", "Oregon Inlet",
                                          "Nags Head", "Duck", "Corolla")}
    fig, off = figure(g, inlet_km, "Cape Point to the Virginia line", lw=2.5,
                      places=over_places, state_y=state_y, boxes=boxes, height=9.4)
    caption(fig, COMMON + " Overview: {n} transects over {km:.0f} km of coast; the "
            "boxes A-D are the regional maps.".format(n=len(g), km=g["alongshore_km"].max())
            + off_clause(off))
    out = save(fig, DIR / "lrr_obx_{0}_map_overview.png".format(TAG), close=True, dpi=300)
    print("Saved", out[0].name)

    for letter, stem, title, sub in reg:
        fig, off = figure(sub, inlet_km, "{0}: {1}".format(letter, title), lw=1.6,
                          places={**PLACES, **MAP_AREAS}, state_y=state_y, height=9.0)
        caption(fig, COMMON + " Region {l}, {t}: {n} transects, {a:.1f}-{b:.1f} km "
                "from Cape Point.".format(l=letter, t=title, n=len(sub),
                                          a=sub["alongshore_km"].min(),
                                          b=sub["alongshore_km"].max())
                + off_clause(off))
        out = save(fig, DIR / "lrr_obx_{0}_map_{1}_{2}.png".format(TAG, letter, stem),
                   close=True, dpi=300)
        print("Saved", out[0].name, len(sub), "transects")


if __name__ == "__main__":
    main()
