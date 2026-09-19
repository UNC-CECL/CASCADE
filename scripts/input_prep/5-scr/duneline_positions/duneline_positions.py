"""
duneline_positions.py
==============================================================================
Where did the dune line sit across Hatteras Island in 1997, 2009 and 2023 --
the lines that stand for the model years 1996, 2010 and 2024? Positions, not
change: a set of maps and alongshore profiles, designed by interview with
Hannah on 2026-09-18.

THE FIGURES   data/hatteras_init/5-scr/4-comparisons/duneline_positions/
    overview/duneline_positions_overview.png
        The island in three north-up segments side by side (south: Cape Point
        to Avon, GIS 1-30; central: GIS 31-60; north: the Tri-Village to GIS
        90), the three dune lines over the island outline and the 90 domain
        boxes, villages named, a scale bar and north arrow per segment, and a
        locator map. Orientation: at this scale the lines overlap; the zooms
        show the metres.
    zooms/zoom_<site>.png, zooms/duneline_positions_zooms.png
        Each a window 3 domains (1.5 km) alongshore by 950 m across (650 m
        landward, 300 m seaward of the 2023 line), every panel the same extent
        and scale, north-up over the 2023 NOAA orthomosaic (D:, read through
        its overviews): the three dune lines, edged for contrast, and NC-12
        white with the year by dash pattern. Named sites Buxton (GIS 1-15),
        Avon (21-31), the Tri-Village (68-83) and Mirlo Beach / the S-curves
        (84-90), each window centred on the site's domain with the largest
        |net dune change| 1997-2023; plus the domain outside them with the
        largest change. The rule and the choices are in
        supporting/zoom_sites.csv. The combined page puts every site on one
        sheet.
    context/dune_to_nc12.png
        Distance from the dune line to the NC-12 centreline per domain, one
        line per year: along each 100 m transect, the road's station minus
        the dune's (both measured from the offshore datum by
        duneline_to_raw_offsets.intersect, the function that builds the dune
        stations), positive where the road is landward of the dune. Road
        line per year: 1978 export for 1997, 2008 export for 2009 (the
        model's ROAD_LINE_FOR_YEAR), today's NC-12 for 2023
        (road_offset/raw_offset/current/).
    context/beach_width.png
        Dune line to CoastSat shoreline per domain, one line per year: along
        each CoastSat transect, the shoreline position (the mean within +/-6
        months of the dune-line image date, 3-rates/coastsat/endpoint) minus
        the distance at which the dune line crosses that transect.
    supporting/   PDFs, CAPTIONS.md, and the tables behind every figure.

YEAR COLOURS: one ordered ramp, light grey 1997 -> slate 2009 -> ink 2023, so
the order reads at a glance and red / blue stay free for seaward / landward.
Every figure is also published to output/figures/shoreline/duneline_positions/.

USAGE
    python scripts/input_prep/5-scr/duneline_positions/duneline_positions.py
    python scripts/input_prep/5-scr/duneline_positions/duneline_positions.py --no-imagery
==============================================================================
"""

from __future__ import annotations

import argparse
import importlib.util
import sys
from pathlib import Path

import numpy as np
import pandas as pd

_REPO = next(_p for _p in Path(__file__).resolve().parents
             if (_p / "pyproject.toml").exists())
sys.path.insert(0, str(_REPO / "scripts"))

import geopandas as gpd  # noqa: E402
import matplotlib  # noqa: E402
matplotlib.use("Agg")
import matplotlib.patheffects as pe  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402
from matplotlib.ticker import MultipleLocator  # noqa: E402
from shapely.geometry import box  # noqa: E402
from shapely.ops import unary_union  # noqa: E402

from site_layer.hat_figure_style import (  # noqa: E402
    DOMAIN_AXIS_LABEL, INK, INK_MUTED, _north_arrow, _scalebar, _title,
    apply_style, caption, figsize, figure_dir, open_frame, save, structures,
    support_dir, town_bands,
)
from site_layer.hat_map_layers import ISLAND_OUTLINE, NC_COAST  # noqa: E402
from site_layer.hat_observed_rates import (  # noqa: E402
    DOMAIN_BOXES, DUNELINE_POSITIONS, TRANSECT_LAYER, coastsat_endpoint_csv,
    transect_lookup,
)
from site_layer.hat_topo_version import (  # noqa: E402
    ROAD_LINE_CURRENT, dune_line_for_year, dune_raw_file, duneline_geojson,
    road_line_file, road_line_for_year,
)
from site_layer.hatteras_site_config import HATTERAS_ANNOTATIONS  # noqa: E402

CRS = "EPSG:6347"                         # NAD83(2011) UTM 18N, the 2023 mosaic's
IMAGERY = Path("D:/Hatteras_GIS/Aerial/2023/2023_full_aerial.tif")
IMAGERY_LABEL = "2023 NOAA NGS orthomosaic"
N = 90
PERIOD_YEARS = (1996, 2010, 2024)
YEAR_C = {1997: "#9a9a9a", 2009: "#4b6a88", 2023: "#141414"}   # the ordered ramp
LINE_LW = 1.3
ROAD_LW = 0.9
OUT = DUNELINE_POSITIONS
PUBLISH = figure_dir("shoreline", "duneline_positions")

SEGMENTS = [("South: Cape Point to Avon", 1, 30),
            ("Central: Avon to the Tri-Village", 31, 60),
            ("North: the Tri-Village to GIS 90", 61, 90)]
NAMED_SITES = [("buxton", "Buxton", 1, 15),
               ("avon", "Avon", 21, 31),
               ("trivillage", "Tri-Village", 68, 83),
               ("mirlo", "Mirlo Beach S-curves", 84, 90)]
PICK_WIDTH = 10                           # domains in a data-picked reach
EXT_M = 400.0                             # CoastSat transects extended landward


def _import(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


D2R = _import("duneline_to_raw_offsets", _REPO / "scripts" / "input_prep"
              / "2-brie-offset" / "1-produce" / "duneline_to_raw_offsets.py")


# -----------------------------------------------------------------------------
# inputs
# -----------------------------------------------------------------------------
def vintages():
    """[(period year, dune vintage, road vintage label, road path), ...]"""
    out = []
    for y in PERIOD_YEARS:
        v = dune_line_for_year(y)
        if y == 2024:
            out.append((y, v, "current", ROAD_LINE_CURRENT))
        else:
            rv = road_line_for_year(y)
            out.append((y, v, str(rv), road_line_file(rv)))
    return out


def _union(path, crs=CRS):
    g = gpd.read_file(path)
    g = g.set_geometry(g.geometry.force_2d()) if hasattr(g.geometry, "force_2d") else g
    return unary_union(g.to_crs(crs).geometry)


def load_layers():
    dom = gpd.read_file(DOMAIN_BOXES).to_crs(CRS)
    dom["gis"] = dom["domain_id"].astype(int)
    dom = dom[dom["gis"].between(1, N)].sort_values("gis")
    outline = gpd.read_file(ISLAND_OUTLINE).to_crs(CRS)
    lines, roads = {}, {}
    for y, v, rlabel, rpath in vintages():
        lines[v] = _union(duneline_geojson(v))
        roads[v] = (rlabel, _union(rpath))
    return dom, outline, lines, roads


# -----------------------------------------------------------------------------
# measurements
# -----------------------------------------------------------------------------
def road_stations(tr, road):
    """Per transect: the station of the SEAWARD-most road crossing (the one
    nearest the dune). duneline_to_raw_offsets.intersect takes the landward-
    most, right for a dune line but not for a road: at Buxton (GIS 8-9) the
    transects also cross the leg of NC-12 that turns west toward Frisco, and
    the landward-most rule put the road 1.6-2 km inland there (first draw,
    2026-09-18). At Rodanthe the 2022 Jug Handle bridge is the only crossing,
    so it is kept."""
    rows = []
    for _, t in tr.iterrows():
        x = t.geometry.intersection(road)
        pts = [] if x.is_empty else ([x] if x.geom_type == "Point" else
                                     [g for g in getattr(x, "geoms", []) if g.geom_type == "Point"])
        st = [t.geometry.project(p) for p in pts]
        rows.append((t.domain_id, t.LineID, float(min(st)) if st else np.nan, len(pts)))
    return pd.DataFrame(rows, columns=["domain_id", "LineID", "road_station_m", "n_crossings"])


def dune_to_road():
    """Per 100 m transect: road station minus dune station (m), per year."""
    tr = D2R.load_transects()
    rows = []
    for y, v, rlabel, rpath in vintages():
        road = _union(rpath, tr.crs)
        rst = road_stations(tr, road)
        raw = (pd.read_csv(dune_raw_file(v), encoding="utf-8-sig")
               .drop_duplicates(subset=["domain_id", "LineID"])
               [["domain_id", "LineID", "ORIG_LEN"]]
               .rename(columns={"ORIG_LEN": "dune_station_m"}))
        m = raw.merge(rst[["domain_id", "LineID", "road_station_m", "n_crossings"]],
                      on=["domain_id", "LineID"], how="left")
        m["dune_to_road_m"] = m["road_station_m"] - m["dune_station_m"]
        m["period_year"], m["dune_vintage"], m["road_line"] = y, v, rlabel
        rows.append(m)
    t = pd.concat(rows, ignore_index=True)
    d = (t.groupby(["dune_vintage", "domain_id"])["dune_to_road_m"]
         .agg(["mean", "min", "max", "count"]).reset_index()
         .rename(columns={"domain_id": "gis", "mean": "dune_to_road_m",
                          "min": "min_m", "max": "max_m", "count": "n_transects"}))
    return t, d


def beach_width():
    """Per CoastSat transect: shoreline chainage minus the chainage where the
    dune line crosses the transect (m), per year."""
    lk = pd.read_csv(transect_lookup())
    lk = lk[lk["domain_number"].between(1, N)]
    import pyogrio
    ids = tuple(lk["transect_id"])
    layer = pyogrio.read_dataframe(TRANSECT_LAYER, columns=["id"],
                                   where=f"id IN ({','.join(repr(i) for i in ids)})")
    layer = layer.to_crs(CRS).set_index("id")
    # shoreline positions at the three dates, from the stored endpoint products
    a = pd.read_csv(coastsat_endpoint_csv(1996, 2010, "transect")).set_index("transect_id")
    b = pd.read_csv(coastsat_endpoint_csv(2010, 2024, "transect")).set_index("transect_id")
    shore = {int(a["start_vintage"].iloc[0]): a["position_start_m"],
             int(a["end_vintage"].iloc[0]): a["position_end_m"],
             int(b["end_vintage"].iloc[0]): b["position_end_m"]}
    dom, _, lines, _ = load_layers()
    rows = []
    from shapely.geometry import LineString
    for tid, gis in zip(lk["transect_id"], lk["domain_number"]):
        if tid not in layer.index:
            continue
        geom = layer.loc[tid, "geometry"]
        # CoastSat transects often START seaward of the dune line (256-470 of
        # 906 missed it in the first draw, 2026-09-18), so extend each one
        # EXT_M landward along its own direction and measure from the ORIGINAL
        # origin: a dune landward of it gets a negative chainage.
        c = np.asarray(geom.coords)
        u = (c[-1] - c[0]) / np.linalg.norm(c[-1] - c[0])
        ext = LineString([tuple(c[0] - EXT_M * u)] + [tuple(p) for p in c])
        for v, line in lines.items():
            x = ext.intersection(line)
            if x.is_empty:
                dune_ch, n = np.nan, 0
            else:
                pts = [x] if x.geom_type == "Point" else [g for g in getattr(x, "geoms", [])
                                                          if g.geom_type == "Point"]
                ch = [ext.project(p) - EXT_M for p in pts]
                dune_ch, n = (float(np.min(ch)), len(pts)) if ch else (np.nan, 0)
            s = shore[v].get(tid, np.nan)
            rows.append((tid, int(gis), v, dune_ch, n, s, s - dune_ch))
    t = pd.DataFrame(rows, columns=["transect_id", "gis", "dune_vintage",
                                    "dune_chainage_m", "n_crossings",
                                    "shoreline_chainage_m", "beach_width_m"])
    d = (t.groupby(["dune_vintage", "gis"])["beach_width_m"]
         .agg(["mean", "min", "max", "count"]).reset_index()
         .rename(columns={"mean": "beach_width_m", "min": "min_m", "max": "max_m",
                          "count": "n_transects"}))
    return t, d


def pick_reach(named):
    """The PICK_WIDTH-domain reach outside the named sites with the largest
    mean |net dune change| 1997-2023."""
    from site_layer.hat_observed_rates import dune_endpoint_csv
    ch = (pd.read_csv(dune_endpoint_csv(1996, 2024, "domain"))
          .set_index("domain_number")["mean_change_m"].abs())
    taken = {g for *_, lo, hi in named for g in range(lo, hi + 1)}
    best = None
    for lo in range(1, N - PICK_WIDTH + 2):
        hi = lo + PICK_WIDTH - 1
        if taken & set(range(lo, hi + 1)):
            continue
        score = float(ch.loc[lo:hi].mean())
        if best is None or score > best[0]:
            best = (score, lo, hi)
    return best


# -----------------------------------------------------------------------------
# drawing
# -----------------------------------------------------------------------------
def _halo(width=2.4, color="white"):
    return [pe.withStroke(linewidth=width, foreground=color, alpha=0.9)]


def _draw_geom(ax, geom, **kw):
    parts = getattr(geom, "geoms", [geom])
    for g in parts:
        if g.geom_type == "LineString":
            x, y = g.xy
            ax.plot(x, y, **kw)
            kw.pop("label", None)


def _year_legend(fig, roads=True, loc="outside lower center"):
    h = [Line2D([], [], color=YEAR_C[v], lw=LINE_LW) for v in YEAR_C]
    lab = [f"dune line {v} (for {y})" for v, y in zip(YEAR_C, PERIOD_YEARS)]
    if roads:
        h += [Line2D([], [], color=INK_MUTED, lw=ROAD_LW, ls=(0, (4, 2)))]
        lab += ["NC-12, same year (dashed)"]
    fig.legend(h, lab, loc=loc, ncol=len(h), frameon=False)


def _imagery(ax, bounds, max_px=1800):
    import rasterio
    from rasterio.enums import Resampling
    from rasterio.windows import from_bounds
    x0, y0, x1, y1 = bounds
    with rasterio.open(IMAGERY) as r:
        win = from_bounds(x0, y0, x1, y1, r.transform)
        scale = max(win.width, win.height) / max_px
        shape = (3, max(1, int(win.height / scale)), max(1, int(win.width / scale)))
        img = r.read([1, 2, 3], window=win, out_shape=shape,
                     resampling=Resampling.average, boundless=True, fill_value=255)
    img = np.where((img == 0).all(axis=0, keepdims=True), 255, img)   # no-data -> white
    ax.imshow(np.transpose(img, (1, 2, 0)), extent=(x0, x1, y0, y1),
              origin="upper", zorder=0, interpolation="bilinear")


def _map_axes(ax, bounds):
    x0, y0, x1, y1 = bounds
    ax.set_xlim(x0, x1)
    ax.set_ylim(y0, y1)
    ax.set_aspect("equal")
    ax.set_xticks([])
    ax.set_yticks([])
    for s in ax.spines.values():
        s.set_color(INK)
        s.set_linewidth(0.6)


def _nice_bar(span_m):
    for L in (100, 200, 250, 500, 1000, 2000, 2500, 5000):
        if L >= span_m / 6:
            return L
    return 5000


def _panel_title(ax, i, text):
    """Letter and title left-aligned on one line, so a narrow map panel
    cannot overlap them (the house _title centres the text)."""
    ax.set_title(f"({chr(97 + i)})  {text}", loc="left", fontsize=9,
                 fontweight="normal", pad=4)


def overview_figure(dom, outline, lines):
    island = (outline.geometry.union_all() if hasattr(outline.geometry, "union_all")
              else unary_union(outline.geometry))
    # ONE extent for every segment (centred on each), so the three panels share
    # a scale and line up; the ocean side carries the village brackets
    ext = []
    for _, lo, hi in SEGMENTS:
        x0, y0, x1, y1 = dom[dom["gis"].between(lo, hi)].total_bounds
        ext.append((x0, y0, x1, y1))
    w = max(e[2] - e[0] for e in ext) + 2600      # west pad + east label room
    h = max(e[3] - e[1] for e in ext) + 800
    fig, axes = plt.subplots(1, 3, figsize=figsize("double", height=6.4),
                             constrained_layout=True)
    for i, (ax, (title, lo, hi), e) in enumerate(zip(axes, SEGMENTS, ext)):
        seg = dom[dom["gis"].between(lo, hi)]
        cx, cy = (e[0] + e[2]) / 2 + 700, (e[1] + e[3]) / 2
        b = (cx - w / 2, cy - h / 2, cx + w / 2, cy + h / 2)
        _map_axes(ax, b)
        gpd.GeoSeries([island], crs=CRS).clip(box(*b)).plot(
            ax=ax, color="0.93", edgecolor="0.62", linewidth=0.5, zorder=1)
        seg.boundary.plot(ax=ax, color="0.74", linewidth=0.35, zorder=2)
        for k, (v, line) in enumerate(lines.items()):
            _draw_geom(ax, line.intersection(box(*b)), color=YEAR_C[v], lw=0.9,
                       zorder=4 + k, solid_capstyle="round")
        east = seg.total_bounds[2]
        for _, r in seg.iterrows():
            if r.gis % 5 == 0 or r.gis in (lo, hi):
                c = r.geometry.centroid
                ax.text(r.geometry.bounds[2] + 120, c.y, str(r.gis), fontsize=6,
                        color=INK_MUTED, va="center", ha="left", zorder=6)
        # villages: a bracket on the ocean side, clear of the domain numbers
        for name, (a, z) in HATTERAS_ANNOTATIONS.town_spans.items():
            if hi < a or z < lo:
                continue
            vs = seg[seg["gis"].between(max(a, lo), min(z, hi))]
            yb0, yb1 = vs.total_bounds[1], vs.total_bounds[3]
            xb = east + 900
            ax.plot([xb, xb], [yb0, yb1], color=INK, lw=0.8, zorder=6,
                    solid_capstyle="butt")
            for yy in (yb0, yb1):
                ax.plot([xb - 120, xb], [yy, yy], color=INK, lw=0.8, zorder=6)
            ax.text(xb + 150, (yb0 + yb1) / 2, name, fontsize=7, color=INK,
                    style="italic", rotation=90, ha="left", va="center", zorder=7)
        _scalebar(ax, 2000, show_cells=False)
        _north_arrow(ax, x=0.86, y=0.07, length=0.035)
        _panel_title(ax, i, title)
    # locator, in the empty sound west of the south segment
    loc = axes[0].inset_axes([0.03, 0.66, 0.34, 0.32])
    coast = gpd.read_file(NC_COAST).to_crs(CRS)
    coast.plot(ax=loc, color="0.86", edgecolor="0.6", linewidth=0.3)
    for (title, lo, hi), ls in zip(SEGMENTS, ("-", "--", ":")):
        bb = dom[dom["gis"].between(lo, hi)].total_bounds
        loc.add_patch(plt.Rectangle((bb[0] - 800, bb[1]), bb[2] - bb[0] + 1600,
                                    bb[3] - bb[1], fill=False, ec=INK, lw=0.7, ls=ls))
    tb = dom.total_bounds
    loc.set_xlim(tb[0] - 30000, tb[2] + 12000)
    loc.set_ylim(tb[1] - 15000, tb[3] + 15000)
    loc.set_aspect("equal")
    loc.set_xticks([]); loc.set_yticks([])
    loc.set_facecolor("white")
    for s in loc.spines.values():
        s.set_linewidth(0.5)
    _year_legend(fig, roads=False)
    caption(fig, (
        "Where the digitized dune line sat along Hatteras Island in 1997, 2009 and "
        "2023, the lines that stand for the model years 1996, 2010 and 2024, "
        "north-up in three segments at one common scale: (a) Cape Point to Avon, "
        "GIS 1–30; (b) Avon to the Tri-Village, GIS 31–60; (c) the Tri-Village to "
        "GIS 90. Light grey 1997, slate 2009, black 2023, over the island outline "
        "(pale fill) and the 500 m model domains (thin boxes, every fifth numbered "
        "on the ocean side); brackets mark the villages. At this scale the three "
        "lines largely overlap; the zoom figures show the separation in metres. "
        "The inset in (a) locates the three segments (solid, dashed, dotted) on "
        "the North Carolina coast. NAD83(2011) / UTM 18N."))
    return _save(fig, "overview", "duneline_positions_overview")


# THE ZOOMS (reworked 2026-09-18 after the first draw). At 5-15 domains a panel
# was 5-7 km tall and the three lines sat on top of each other; the point of a
# zoom is the tens of metres between them. Each site now shows a 3-domain
# (1.5 km) WINDOW centred on the domain of that site with the largest
# |net dune change| 1997-2023, every panel at the same extent and scale.
ZOOM_HALF = 1                             # domains either side of the centre
ROAD_DASH = {1997: (0, (1, 1.5)), 2009: (0, (3, 1.8)), 2023: (0, (7, 2.5))}
EDGE = {1997: "#1a1a1a", 2009: "#1a1a1a", 2023: "white"}


def _change():
    from site_layer.hat_observed_rates import dune_endpoint_csv
    return (pd.read_csv(dune_endpoint_csv(1996, 2024, "domain"))
            .set_index("domain_number")["mean_change_m"])


def zoom_sites():
    """[(key, title, centre GIS, site span lo, hi, why), ...]"""
    ch = _change()
    out = []
    for key, title, lo, hi in NAMED_SITES:
        span = ch.loc[max(lo, 1 + ZOOM_HALF):min(hi, N - ZOOM_HALF)]
        c = int(span.abs().idxmax())
        out.append((key, title, c, lo, hi,
                    f"named site GIS {lo}-{hi}; centred on its largest |net dune "
                    f"change| 1997-2023, GIS {c} ({ch.loc[c]:+.1f} m)"))
    taken = {g for *_, lo, hi in NAMED_SITES
             for g in range(lo - ZOOM_HALF, hi + ZOOM_HALF + 1)}
    free = ch[[g for g in ch.index
               if g not in taken and 1 + ZOOM_HALF <= g <= N - ZOOM_HALF]]
    c = int(free.abs().idxmax())
    out.append(("picked", "Largest change elsewhere", c, c, c,
                f"the domain outside every named site with the largest |net dune "
                f"change| 1997-2023, GIS {c} ({ch.loc[c]:+.1f} m)"))
    return out


ZOOM_LAND_M, ZOOM_SEA_M = 650.0, 300.0     # window either side of the dune


def _zoom_extent(dom, lines, roads, centre):
    """A fixed window, 3 domains alongshore by 950 m cross-shore, set on the
    2023 dune line: 650 m landward (NC-12 is 270 m behind the dune on
    average, so it is usually in view) and 300 m seaward. The first rework
    widened each window to take every road in, which at Buxton (the Frisco
    turn) and Rodanthe (the Jug Handle) stretched every panel and squeezed
    the lines back together."""
    seg = dom[dom["gis"].between(centre - ZOOM_HALF, centre + ZOOM_HALF)]
    x0, y0, x1, y1 = seg.total_bounds
    part = lines[2023].intersection(box(x0 - 3000, y0, x1 + 3000, y1))
    cx = part.centroid.x
    return (cx - ZOOM_LAND_M, y0 - 40, cx + ZOOM_SEA_M, y1 + 40)


def zoom_panel(ax, dom, lines, roads, b, imagery=True):
    _map_axes(ax, b)
    if imagery:
        _imagery(ax, b, max_px=1400)
    else:
        ax.set_facecolor("0.96")
    clip = box(*b)
    seg = dom[dom.geometry.intersects(clip)]
    seg.boundary.plot(ax=ax, color="white", linewidth=0.5, alpha=0.6, zorder=2)
    for v in YEAR_C:
        _draw_geom(ax, roads[v][1].intersection(clip), color="white", lw=1.0,
                   ls=ROAD_DASH[v], zorder=3, alpha=0.95)
    for k, v in enumerate(YEAR_C):
        _draw_geom(ax, lines[v].intersection(clip), color=YEAR_C[v], lw=1.6,
                   zorder=5 + k, solid_capstyle="round",
                   path_effects=[pe.withStroke(linewidth=2.8, foreground=EDGE[v])])
    for _, r in seg.iterrows():
        c = r.geometry.centroid
        if b[1] < c.y < b[3]:
            ax.text(b[2] - 0.03 * (b[2] - b[0]), c.y, f"GIS {r.gis}", fontsize=6.5,
                    color="white", ha="right", va="center", zorder=9,
                    path_effects=_halo(1.8, "black"))
    _scalebar(ax, 200, show_cells=False)
    _north_arrow(ax, x=0.08, y=0.80, length=0.06)


def _zoom_legend(fig):
    h = [Line2D([], [], color=YEAR_C[v], lw=1.6,
                path_effects=[pe.withStroke(linewidth=2.8, foreground=EDGE[v])])
         for v in YEAR_C]
    lab = [f"dune line {v} (for {y})" for v, y in zip(YEAR_C, PERIOD_YEARS)]
    h += [Line2D([], [], color="0.35", lw=1.0, ls=ROAD_DASH[v]) for v in YEAR_C]
    lab += ["NC-12, 1978 alignment (1997)", "NC-12, 2008 alignment (2009)",
            "NC-12, today (2023)"]
    fig.legend(h, lab, loc="outside lower center", ncol=3, frameon=False,
               fontsize=7.5)


def zoom_figures(dom, lines, roads, sites, imagery):
    exts = {s[0]: _zoom_extent(dom, lines, roads, s[2]) for s in sites}
    # one size for every panel: the widest window sets the width
    W = max(e[2] - e[0] for e in exts.values())
    H = max(e[3] - e[1] for e in exts.values())
    boxes = {}
    for k, e in exts.items():
        cx, cy = (e[0] + e[2]) / 2, (e[1] + e[3]) / 2
        boxes[k] = (cx - W / 2, cy - H / 2, cx + W / 2, cy + H / 2)
    written = []
    for key, title, c, lo, hi, why in sites:
        fig, ax = plt.subplots(figsize=figsize("single", aspect=min(H / W, 1.6) + 0.15),
                               constrained_layout=True)
        zoom_panel(ax, dom, lines, roads, boxes[key], imagery)
        ax.set_title(f"{title}, GIS {c - ZOOM_HALF}–{c + ZOOM_HALF}", loc="left",
                     fontsize=9, pad=4)
        _zoom_legend(fig)
        caption(fig, _zoom_caption(f"{title}, GIS {c - ZOOM_HALF}–{c + ZOOM_HALF} "
                                   f"({why})", imagery))
        written += _save(fig, "zooms", f"zoom_{key}")
    n, ncol = len(sites), 3
    nrow = int(np.ceil(n / ncol))
    pw = 7.48 / ncol
    fig, axes = plt.subplots(nrow, ncol, figsize=(7.48, nrow * (pw * H / W + 0.35) + 0.8),
                             constrained_layout=True, squeeze=False)
    for i, (ax, (key, title, c, lo, hi, why)) in enumerate(zip(axes.flat, sites)):
        zoom_panel(ax, dom, lines, roads, boxes[key], imagery)
        _panel_title(ax, i, f"{title}, GIS {c - ZOOM_HALF}–{c + ZOOM_HALF}")
    for ax in list(axes.flat)[n:]:
        ax.axis("off")
    _zoom_legend(fig)
    caption(fig, ("The dune line in 1997, 2009 and 2023 at five sites, each a "
                  "window three model domains (1.5 km) alongshore by 950 m across, "
                  "650 m landward and 300 m seaward of the 2023 line, at one scale: "
                  + "; ".join(f"({chr(97 + i)}) {t}, centred on GIS {c}"
                              for i, (_, t, c, *_r) in enumerate(sites))
                  + ". Each named site's window is centred on its domain with the "
                  "largest net dune change 1997–2023; (e) is the domain with the "
                  "largest change outside the named sites. "
                  + _zoom_caption(None, imagery, sheet=True)))
    written += _save(fig, "zooms", "duneline_positions_zooms")
    return written


def _zoom_caption(what, imagery, sheet=False):
    head = "" if sheet else f"The dune line at {what}, in 1997, 2009 and 2023. "
    return (head + "Solid lines: the digitized dune lines, light grey 1997, slate "
            "2009 and black 2023 (standing in for the model years 1996, 2010 and "
            "2024), edged for contrast. White dashed lines: the NC-12 centreline, "
            "dotted for the 1978 alignment (used for 1997), short dashes for 2008 "
            "(used for 2009) and long dashes for today's NCDOT alignment (2023). "
            "Thin white boxes are the 500 m model domains. "
            + (f"Background: the {IMAGERY_LABEL}, so the 2023 line can be checked "
               "against the dune it traces; the 1997 and 2009 lines are not traced "
               "from this image. " if imagery else "")
            + "Scale bar in metres; north up; NAD83(2011) / UTM 18N.")


def context_figure(table, value, ylabel, stem, what, how):
    fig, ax = plt.subplots(figsize=figsize("double", aspect=0.40),
                           constrained_layout=True)
    ax.set_xlim(0.5, N + 0.5)
    vals = table[value].to_numpy(float)
    lo = np.nanmin(vals)
    hi = np.nanmax(vals)
    step = 50.0 if hi - lo > 200 else 20.0
    ax.set_ylim(min(0.0, np.floor(lo / step) * step), np.ceil(hi / step) * step)
    town_bands(ax, label=True)
    ax.axhline(0, color=INK_MUTED, lw=0.6, zorder=2)
    for k, v in enumerate(YEAR_C):
        t = table[table["dune_vintage"] == v].set_index("gis").reindex(range(1, N + 1))
        ax.plot(t.index, t[value], color=YEAR_C[v], lw=1.2, zorder=4 + k,
                label=f"{v} (for {PERIOD_YEARS[k]})")
    structures(ax, label=True)
    ax.xaxis.set_major_locator(MultipleLocator(10))
    ax.xaxis.set_minor_locator(MultipleLocator(5))
    ax.yaxis.set_major_locator(MultipleLocator(step))
    ax.yaxis.grid(True, zorder=0)
    ax.set_axisbelow(True)
    open_frame(ax)
    ax.set_xlabel(DOMAIN_AXIS_LABEL)
    ax.set_ylabel(ylabel)
    fig.legend(*ax.get_legend_handles_labels(), loc="outside lower center",
               ncol=3, frameon=False)
    means = "; ".join(f"{v} {table[table.dune_vintage == v][value].mean():.0f} m"
                      for v in YEAR_C)
    caption(fig, (f"{what} by GIS domain (1 at Cape Point, 90 at Pea Island) in "
                  "1997, 2009 and 2023, the dune lines that stand for the model years "
                  f"1996, 2010 and 2024: {how} Each line is the domain mean. Island "
                  f"means: {means}. Village spans are shaded; the solid hairline is "
                  "the Buxton groin and the dotted hairlines are the Avon and "
                  "Rodanthe piers."))
    return _save(fig, "context", stem)


def _save(fig, sub, stem):
    out = save(fig, OUT / sub / stem)
    out += save(fig, PUBLISH / stem)
    plt.close(fig)
    return out


# -----------------------------------------------------------------------------
def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description="dune-line positions, 1997/2009/2023")
    ap.add_argument("--no-imagery", action="store_true",
                    help="draw the zooms without the 2023 orthomosaic")
    a = ap.parse_args(argv)
    imagery = not a.no_imagery and IMAGERY.is_file()
    if not a.no_imagery and not imagery:
        print(f"NOTE: {IMAGERY} not found; zooms drawn without imagery")

    apply_style()
    sup = support_dir(OUT)
    dom, outline, lines, roads = load_layers()

    t, d = dune_to_road()
    t.to_csv(sup / "dune_to_nc12_transects.csv", index=False, float_format="%.2f")
    d.to_csv(sup / "dune_to_nc12_domains.csv", index=False, float_format="%.2f")
    bt, bd = beach_width()
    bt.to_csv(sup / "beach_width_transects.csv", index=False, float_format="%.2f")
    bd.to_csv(sup / "beach_width_domains.csv", index=False, float_format="%.2f")

    sites = zoom_sites()
    pd.DataFrame([dict(site=k, title=t_, centre_gis=c, first_gis=c - ZOOM_HALF,
                       last_gis=c + ZOOM_HALF, site_span=f"{lo}-{hi}", chosen_by=why)
                  for k, t_, c, lo, hi, why in sites]).to_csv(sup / "zoom_sites.csv", index=False)

    written = overview_figure(dom, outline, lines)
    written += zoom_figures(dom, lines, roads, sites, imagery)
    written += context_figure(
        d, "dune_to_road_m", "Dune line to NC-12 (m)", "dune_to_nc12",
        "Distance from the dune line to the NC-12 centreline",
        "along each 100 m transect, the road's distance from the fixed offshore "
        "datum minus the dune line's (both found by the same line-transect "
        "intersection), positive where the road lies landward of the dune. The "
        "road line for 1997 is the 1978 export and for 2009 the 2008 export (the "
        "alignments the model uses); for 2023 it is today's NCDOT alignment.")
    written += context_figure(
        bd, "beach_width_m", "Dune line to shoreline (m)", "beach_width",
        "Beach width, the distance from the dune line to the CoastSat shoreline,",
        "along each CoastSat transect, the shoreline position (the mean satellite "
        "position within six months of the dune-line image date) minus the "
        "distance at which the dune line crosses that transect, positive where "
        "the shoreline lies seaward of the dune. The 2023 image date is not known "
        "and is assumed to be 1 July.")

    for v in YEAR_C:
        dd = d[d.dune_vintage == v]["dune_to_road_m"]
        bb = bd[bd.dune_vintage == v]["beach_width_m"]
        print(f"{v}  dune->NC-12 mean {dd.mean():6.1f} m (min {dd.min():6.1f})   "
              f"beach width mean {bb.mean():6.1f} m (min {bb.min():6.1f})")
    print("zoom centres: " + ", ".join(f"{k} GIS {c}" for k, _, c, *_ in sites))
    print(f"{len(written) // 2} figures written (+ published copies)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
