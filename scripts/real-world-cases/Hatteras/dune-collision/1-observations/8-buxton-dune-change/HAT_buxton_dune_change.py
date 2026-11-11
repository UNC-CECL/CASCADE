"""
North of Buxton (GIS 10-15): every 1996-2009 lidar survey, to see what happened to the dunes that the collision
module cannot lower in time (they fall 1-2.5 m at GIS 10-12 in 1998-1999, before or with Dennis, while GIS 13-15 hold).

    python scripts/real-world-cases/Hatteras/dune-collision/1-observations/8-buxton-dune-change/HAT_buxton_dune_change.py

Asked 2026-11-11 (dune-collision/4-hindcast/TESTING_2026-11-09_to_11.md, section 4): were the dunes scarped in place
(storms), truncated near the water level (overtopped), moved, or reshaped with the 1999 NC-12 relocation (GIS 9-14)?

Every survey on one datum, the 2009 lidar's: each survey in its own vertical unit and with the model's vertical offset
(survey_sources.csv of the 1996 start, from HAT_survey_dune_crests.SURVEY_VERTICAL_OFFSET_M: 1996 -0.24, 1997 -0.34,
1999 -0.32 m), then minus what is left of its difference from the 2009 lidar on NC-12 (road_offset_vs_2009_m, the
island's road pixels outside GIS 9-13, where the road moved), so every survey sits within 0.02 m of 2009 on the road.
That last step goes beyond the model, which leaves 2000 (+0.10 m) and 2001 (-0.10 m) as measured. All resampled
bilinearly onto the model clips' 1 m grid; the surveys are 9 ft (about 2.7 m; 1996, 1998, 1999, 2004), 6 ft (2005), 3 m (1997, 2000),
4 m (2001) and 1 m (2009) grids, so narrow dunes look smoother in the coarse years. 2001 is unclassified lidar:
buildings and vegetation can stand out as false highs.

Writes, in the workspace's dune-collision/1-observations/8-buxton-dune-change/:
- buxton_dunes_by_survey.png: the coast unrolled along NC-12, one strip per survey (terrain colours, hillshade, the
  cross-shore stretched), the berm (1.7 m) contour, NC-12 in its
  1978 (dashed) and 2008 (solid) alignments, the domain boundaries
- buxton_dune_change_between_surveys.png: each survey minus the one before
- buxton_dune_profiles.png: one cross-shore profile per domain (the 100 m transect nearest its middle), every survey,
  with the storm series' 1998 and 1999 peak water levels (Bonnie 2.64 m, Dennis 3.65 m NAVD88)
- supporting/buxton_datum_check.csv: each survey's unit, offsets and its difference from 2009 on NC-12 in GIS 15
  (where the road did not move), after the correction

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-11-11
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

_HERE = Path(__file__).resolve()
_HAT_DIR = next(p for p in _HERE.parents if (p / "hatteras_paths.py").exists())
sys.path.insert(0, str(_HAT_DIR))
sys.path.insert(0, str(_HAT_DIR / "input_prep" / "1-barrier3d-domains" / "2-survey-crests"))
from hatteras_paths import WORKSPACE  # noqa: E402
import HAT_survey_dune_crests as S  # noqa: E402

# --- CONFIG ------------------------------------------------------------------
OUT = WORKSPACE / "dune-collision" / "1-observations" / "8-buxton-dune-change"
START = "1996-start"
GIS = range(10, 16)
RES_M = 1.0                                     # the common grid: the model clips' 1 m lattice
PAD_M = 150.0                                   # around the domains
MAP_LAND_M, MAP_SEA_M = 120.0, 260.0            # the plan-view maps (maps(), changes(); not drawn by main)
STRIP_LAND_M, STRIP_SEA_M = 60.0, 180.0         # the unrolled strips: this far landward and seaward of NC-12 (2008)
STRIP_STEP_M = 1.0
# NC-12's relocation in GIS 9-14 (hatteras_site_config RelocationEvent, 1999): after Hurricane Dennis (late August to
# early September 1999), about 3,500 ft just north of Buxton, and in place in an NCDOT photo of 16 March 2000 (NC 12
# Buxton to Avon Feasibility Study, NCDOT, December 2015, section 2.1 and Figure 4; the management timeline's 1999
# row). The September 1999 lidar, flown just after Dennis, still has the old alignment (washed out); the February-
# March 2000 survey is the first with the new one. The strips draw the road where it was in each survey
RELOCATED_BY_SURVEY = 2000
ROAD_C, BOUNDARY_C = "#AE017E", "0.45"           # NC-12 (magenta), the domain boundaries (grey)
MONTHS = {1996: "Aug-Oct", 1997: "Sep-Oct", 1998: "Sep", 1999: "Sep", 2000: "Feb-Mar", 2001: "Jan-Mar",
          2004: "Apr-Sep", 2005: "Aug-Nov", 2009: "Aug"}   # the flights (survey_sources.csv)
WASHED_OUT_SURVEY = 1999
TERRAIN_RANGE = (-1.0, 9.0)                     # m NAVD88, the terrain colour scale (blue below about 0)
ELEV_RANGE = (0.0, 9.0)                         # m NAVD88, the map colour scale
DIFF_MAX = 3.0                                  # m, the change colour scale
BERM_M = 1.7                                    # Barrier3D's berm, the module's toe [m NAVD88]
STORMS = (("Bonnie 1998 peak", 2.64, ":"), ("Dennis 1999 peak", 3.65, "--"))   # the storm series' peak water levels
ROAD_CHECK_GIS = (15,)                          # NC-12 did not move here in 1999
SOURCES = S.tv.product_dir(START) / "survey-crests" / "survey_sources.csv"
TRANSECTS = WORKSPACE / "data" / "hatteras_init" / "2-brie-offset" / "transects" / "transects_100m.geojson"
ROADS = {"1978": S.tv.MGMT_ROOT / "road_offset" / "raw_offset" / "1978",
         "2008": S.tv.MGMT_ROOT / "road_offset" / "raw_offset" / "2008"}
# -----------------------------------------------------------------------------


def domains():
    import geopandas as gpd
    d = gpd.read_file(S.GIS_ROOT / "domains.geojson")
    return d[d.domain_id.isin(GIS)]


def grid():
    """One 1 m grid over the model's own domain clips of GIS 10-15 (their CRS and pixel lattice: the grid
    S.raster_on_clip expects), padded by PAD_M."""
    from affine import Affine
    clips = S.product(S.STARTS[START][0]).gapfill_1m
    grids = [S.clip_grid(clips / f"clip_domain_{g}_filled.tif") for g in GIS]
    crs = grids[0][2]
    x0 = min(t.c for t, _, _ in grids) - PAD_M
    y1 = max(t.f for t, _, _ in grids) + PAD_M
    x1 = max(t.c + w * t.a for t, (h, w), _ in grids) + PAD_M
    y0 = min(t.f + h * t.e for t, (h, w), _ in grids) - PAD_M
    return Affine(RES_M, 0, x0, 0, -RES_M, y1), (int(round((y1 - y0) / RES_M)), int(round((x1 - x0) / RES_M))), crs


def road_lines(crs):
    import geopandas as gpd
    out = {}
    for name, folder in ROADS.items():
        files = sorted(folder.glob("*.geojson")) + sorted(folder.glob("*.shp"))
        if files:
            out[name] = gpd.read_file(files[0]).to_crs(crs)
    return out


def surveys():
    """Each survey on the 2009 datum, on the common grid: {year: (array m NAVD88, flown, row of the datum table)}."""
    import rasterio
    src_tab = pd.read_csv(SOURCES).set_index("survey")
    entries = S.STARTS[START][1]
    dom = domains()
    g = grid()
    out = {}
    for year in sorted(entries):
        kind, path, flown, *_ = entries[year]
        files = sorted(path.parent.glob(path.name)) if "*" in path.name else [path]
        srcs = [rasterio.open(f) for f in files]
        try:
            r = src_tab.loc[year] if year in src_tab.index else None
            unit = "m" if kind == "mosaic" else (r.vertical_unit if r is not None else S.stated_unit(srcs[0]))
            scale = S.US_FT if unit == "ft" else 1.0
            shift = float(r.vertical_offset_applied_m) if r is not None else 0.0
            resid = float(r.road_offset_vs_2009_m) if r is not None and pd.notna(pd.to_numeric(
                r.road_offset_vs_2009_m, errors="coerce")) else 0.0
            z = S.raster_on_clip(srcs, g) * scale + shift - resid
            z[(z < -5) | (z > 30)] = np.nan
            out[year] = (z, flown, dict(survey=year, flown=flown, unit=unit, model_offset_m=shift,
                                        road_residual_removed_m=-resid, resolution_m=float(srcs[0].res[0] * (
                                            S.US_FT if srcs[0].crs and srcs[0].crs.linear_units.lower().startswith(
                                                ("foot", "us survey")) else 1.0))))
        finally:
            for s in srcs:
                s.close()
    return out, g, dom


def datum_check(surv, g, roads):
    from rasterio.features import rasterize
    t, shape, crs = g
    dom = domains().to_crs(crs)
    keep = rasterize([(geom, 1) for geom in dom[dom.domain_id.isin(ROAD_CHECK_GIS)].geometry], out_shape=shape,
                     transform=t, fill=0).astype(bool)
    road = rasterize([(geom.buffer(4), 1) for geom in roads["2008"].geometry], out_shape=shape, transform=t,
                     fill=0).astype(bool) & keep
    ref = surv[2009][0]
    rows = []
    for year, (z, _, info) in surv.items():
        d = (z - ref)[road & np.isfinite(z) & np.isfinite(ref)]
        rows.append(info | dict(road_px_gis15=int(d.size),
                                minus_2009_on_road_gis15_m=float(np.median(d)) if d.size else np.nan))
    t_ = pd.DataFrame(rows)
    (OUT / "supporting").mkdir(parents=True, exist_ok=True)
    t_.round(3).to_csv(OUT / "supporting" / "buxton_datum_check.csv", index=False)
    print(t_.round(3).to_string(index=False))
    return t_


def _decorate(ax, g, dom, roads, lw=0.6):
    t, (h, w), crs = g
    for _, r in dom.to_crs(crs).iterrows():
        x, y = r.geometry.exterior.xy
        ax.plot(x, y, color="0.25", lw=0.4)
    for name, ls in (("1978", (0, (4, 2))), ("2008", "-")):
        for geom in (roads[name].geometry if name in roads else []):
            for part in getattr(geom, "geoms", [geom]):
                x, y = part.xy
                ax.plot(x, y, color="#1A1A1A", lw=lw, ls=ls)
    ax.set_xlim(*_xlim(g, roads))
    ax.set_ylim(t.f + h * t.e, t.f)
    ax.set_aspect("equal")
    ax.set_xticks([])
    ax.set_yticks([])


def _xlim(g, roads):
    """The beach-dune-road zone only: from MAP_LAND_M landward of the westernmost NC-12 to MAP_SEA_M seaward of the
    easternmost, within the grid's rows (the ocean is to the east)."""
    from shapely.geometry import box
    t, (h, w), _ = g
    win = box(t.c, t.f + h * t.e, t.c + w * t.a, t.f)
    xs = [c for name in roads for geom in roads[name].geometry
          for c in (lambda b: (b[0], b[2]))(geom.intersection(win).bounds) if not geom.intersection(win).is_empty]
    return min(xs) - MAP_LAND_M, max(xs) + MAP_SEA_M


def _labels(ax, g, dom, roads):
    t, _, crs = g
    x0 = _xlim(g, roads)[0]
    for _, r in dom.to_crs(crs).iterrows():
        c = r.geometry.centroid
        ax.text(x0 + 20, c.y, f"{int(r.domain_id)}", fontsize=7, va="center", color="0.15", fontweight="bold")


def maps(surv, g, dom, roads, S_, plt):
    from matplotlib.colors import LightSource
    years = sorted(surv)
    t, shape, crs = g
    ls = LightSource(azdeg=315, altdeg=45)
    cmap = plt.get_cmap("viridis")
    fig, axs = plt.subplots(1, len(years), figsize=(S_.FIG_W_DOUBLE, 8.6), layout="constrained")
    extent = (t.c, t.c + shape[1] * t.a, t.f + shape[0] * t.e, t.f)
    for ax, y in zip(axs, years):
        z = surv[y][0]
        zz = np.where(np.isfinite(z), z, ELEV_RANGE[0])
        rgb = ls.shade(np.clip(zz, *ELEV_RANGE), cmap=cmap, vert_exag=2, blend_mode="soft",
                       vmin=ELEV_RANGE[0], vmax=ELEV_RANGE[1], dx=RES_M, dy=RES_M)
        rgb[~np.isfinite(z)] = 1.0
        rgb[np.isfinite(z) & (z < 0.36)] = (0.85, 0.92, 0.97, 1.0)          # below MHW: water
        ax.imshow(rgb, extent=extent, interpolation="nearest")
        ax.contour(np.where(np.isfinite(z), z, np.nan), levels=[BERM_M], colors="white", linewidths=0.5,
                   extent=extent, origin="upper")
        _decorate(ax, g, dom, roads)
        ax.set_title(f"{y}{'*' if y == 2001 else ''}\n{surv[y][1]}", fontsize=7)
    _labels(axs[0], g, dom, roads)
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=plt.Normalize(*ELEV_RANGE))
    fig.colorbar(sm, ax=axs, orientation="horizontal", location="bottom", shrink=0.5,
                 label="Elevation (m NAVD88, every survey on the 2009 datum)")
    out = OUT / "buxton_dunes_by_survey.png"
    S_.save(fig, out, bbox_inches="tight", close=True)
    S_.record_caption(out, (
        "GIS 10-15, north of Buxton, in every 1996-2009 lidar survey (flight dates under each year), elevation with "
        "hillshade on one colour scale, every survey on the 2009 lidar's datum (the model's vertical offsets, then each "
        "survey's remaining difference from 2009 on NC-12 removed; supporting/buxton_datum_check.csv). White line: the "
        "berm, 1.7 m NAVD88, the collision module's toe; light blue: below mean high water; black: NC-12 in its 1978 "
        "(dashed) and 2008 (solid) alignments (relocated in GIS 9-14 in 1999); thin grey: the domain boundaries, "
        "numbered. The surveys are on grids of 1 m (2009) to 9 ft (about 2.7 m), all resampled to the model clips' 1 m grid. *2001 is "
        "unclassified lidar: buildings and vegetation can stand out as false highs."))
    print(f"wrote {out}")


class Strip:
    """The coast unrolled along NC-12 (2008 alignment): alongshore distance north from GIS 10's south edge against
    cross-shore distance from the road, seaward positive, sampled every STRIP_STEP_M; the cross-shore axis is drawn
    stretched so the dunes, tens of metres wide over 3 km of coast, can be seen."""

    def __init__(self, g, roads, dom):
        from shapely.geometry import box
        from shapely.ops import linemerge
        t, (h, w), crs = g
        line = roads["2008"].geometry.intersection(box(t.c, t.f + h * t.e, t.c + w * t.a, t.f)).union_all()
        if line.geom_type == "MultiLineString":
            line = linemerge(line)
        if line.geom_type != "LineString":
            line = max(line.geoms, key=lambda x: x.length)
        c = np.asarray(line.coords)
        if c[0, 1] > c[-1, 1]:                                        # run south to north
            line = type(line)(c[::-1])
        d = dom.to_crs(crs)
        s0 = line.project(d[d.domain_id == min(GIS)].geometry.iloc[0].centroid) - 250
        s1 = line.project(d[d.domain_id == max(GIS)].geometry.iloc[0].centroid) + 250
        self.s = np.arange(max(s0, 0), min(s1, line.length), STRIP_STEP_M)
        self.n = np.arange(-STRIP_LAND_M, STRIP_SEA_M + STRIP_STEP_M, STRIP_STEP_M)
        p = np.array([line.interpolate(v).coords[0] for v in self.s])
        tan = np.gradient(p, axis=0)
        tan /= np.hypot(tan[:, 0], tan[:, 1])[:, None]
        nrm = np.column_stack([tan[:, 1], -tan[:, 0]])               # right of north-going: east, the ocean
        nrm *= np.sign(nrm[:, 0])[:, None]
        xy = p[None, :, :] + self.n[:, None, None] * nrm[None, :, :]
        self.col, self.row = (xy[..., 0] - t.c) / t.a - 0.5, (xy[..., 1] - t.f) / t.e - 0.5
        # domain boundaries and middles along the road
        # each domain's stretch of the road: where NC-12 runs through the domain's box
        self.edges, self.mids, self.spans = [], {}, {}
        for _, r in d.iterrows():
            part = line.intersection(r.geometry)
            pts = [c for p in getattr(part, "geoms", [part]) if not p.is_empty for c in p.coords]
            if not pts:
                continue
            ss = [line.project(type(r.geometry.centroid)(*c[:2])) for c in pts]
            self.spans[int(r.domain_id)] = (min(ss), max(ss))
            self.mids[int(r.domain_id)] = (min(ss) + max(ss)) / 2
            self.edges += [min(ss), max(ss)]
        self.road78 = None
        if "1978" in roads:
            from shapely.geometry import Point
            old = roads["1978"].geometry.union_all()
            self.road78 = np.array([(Point(v).distance(old) * (1 if old.distance(Point(*(np.array(v) + 5 * k)))
                                                                 < Point(v).distance(old) else -1))
                                    for v, k in zip(p, nrm)])

    def boundaries(self):
        """One alongshore position per boundary between neighbouring domains (the mean of the two domains' facing
        edges), plus the outer edges."""
        sp = sorted(self.spans.values())
        return [sp[0][0], *[(a[1] + b[0]) / 2 for a, b in zip(sp[:-1], sp[1:])], sp[-1][1]]

    def sample(self, z):
        from scipy.ndimage import map_coordinates
        return map_coordinates(z, [self.row, self.col], order=1, mode="constant", cval=np.nan)

    def draw(self, ax, year):
        """The road where it was in that survey, the domain boundaries, and the axis limits (seaward down)."""
        ext = (self.s[0] / 1000, self.s[-1] / 1000, self.n[0], self.n[-1])
        old = year < RELOCATED_BY_SURVEY and self.road78 is not None
        y = self.road78 if old else np.zeros_like(self.s)
        ax.plot(self.s / 1000, y, color=ROAD_C, lw=1.0, ls=(0, (4, 2)) if year == WASHED_OUT_SURVEY else "-",
                solid_capstyle="butt", zorder=8)
        for e in self.boundaries()[1:-1]:
            ax.axvline(e / 1000, color=BOUNDARY_C, lw=0.6, zorder=6)
        ax.set_xlim(ext[0], ext[1])
        ax.set_ylim(ext[3], ext[2])                            # seaward down: the ocean at the bottom
        ax.set_yticks([0, 100])
        ax.tick_params(axis="both", labelsize=7, length=2.5, width=0.5)
        for sp in ax.spines.values():
            sp.set_linewidth(0.5)
        return ext

    def top_axis(self, ax):
        """The domain numbers once, along the top of the first strip, between the boundaries."""
        b = self.boundaries()
        sec = ax.secondary_xaxis("top")
        sec.set_xticks([e / 1000 for e in b[1:-1]], labels=[])
        sec.tick_params(length=3, width=0.5)
        for gis, m in self.mids.items():
            ax.text(m / 1000, 1.06, f"{gis}", transform=ax.get_xaxis_transform(), ha="center", va="bottom",
                    fontsize=7.5)
        ax.text(1.0, 1.06, " GIS domain", transform=ax.transAxes, ha="left", va="bottom", fontsize=7,
                color="0.3", style="italic")


def _row_label(ax, text, sub=None):
    ax.set_ylabel(text if not sub else f"{text}\n{sub}", rotation=0, ha="right", va="center", fontsize=8,
                  labelpad=8, linespacing=1.1)


def _legend(fig, extra=()):
    from matplotlib.lines import Line2D
    h = [Line2D([], [], color=ROAD_C, lw=1.0, label="NC-12"),
         Line2D([], [], color=ROAD_C, lw=1.0, ls=(0, (4, 2)), label="NC-12, washed out by Dennis (1999)"),
         Line2D([], [], color=BOUNDARY_C, lw=0.6, label="Domain boundary"), *extra]
    fig.legend(handles=h, loc="outside lower center", ncol=len(h), fontsize=7, frameon=False, handlelength=2.2,
               columnspacing=1.6)


def _colorbar(fig, axs, mappable, label, extend="neither"):
    cb = fig.colorbar(mappable, ax=axs, orientation="horizontal", location="bottom", shrink=0.45, aspect=40,
                      pad=0.012, extend=extend)
    cb.set_label(label, fontsize=8)
    cb.ax.tick_params(labelsize=7, width=0.5, length=2.5)
    cb.outline.set_linewidth(0.5)


def strips(surv, g, dom, roads, S_, plt):
    from matplotlib.colors import LightSource
    from matplotlib.lines import Line2D
    st = Strip(g, roads, dom)
    years = sorted(surv)
    ls = LightSource(azdeg=315, altdeg=45)
    cmap = plt.get_cmap("terrain")
    fig, axs = plt.subplots(len(years), 1, figsize=(S_.FIG_W_DOUBLE, 0.82 * len(years) + 1.3), sharex=True,
                            layout="constrained")
    for ax, y in zip(axs, years):
        z = st.sample(surv[y][0])
        zz = np.where(np.isfinite(z), z, TERRAIN_RANGE[0])
        rgb = ls.shade(np.clip(zz, *TERRAIN_RANGE), cmap=cmap, vert_exag=1, blend_mode="soft",
                       vmin=TERRAIN_RANGE[0], vmax=TERRAIN_RANGE[1], dx=STRIP_STEP_M, dy=STRIP_STEP_M)
        rgb[~np.isfinite(z)] = 1.0
        ext = (st.s[0] / 1000, st.s[-1] / 1000, st.n[0], st.n[-1])
        ax.imshow(rgb, extent=ext, origin="lower", aspect="auto", interpolation="nearest")
        ax.contour(st.s / 1000, st.n, np.where(np.isfinite(z), z, np.nan), levels=[BERM_M], colors="white",
                   linewidths=0.4)
        st.draw(ax, y)
        _row_label(ax, f"{y}{'*' if y == 2001 else ''}", MONTHS.get(y))
    st.top_axis(axs[0])
    axs[-1].set_xlabel("Alongshore distance along NC-12, south to north (km)", fontsize=8)
    fig.supylabel("Cross-shore distance from NC-12 (2008 alignment), seaward down (m)", fontsize=8)
    _colorbar(fig, axs, plt.cm.ScalarMappable(cmap=cmap, norm=plt.Normalize(*TERRAIN_RANGE)), "Elevation (m NAVD88)")
    _legend(fig, (Line2D([], [], color="0.6", lw=0.6, label="Berm, 1.7 m (module toe; white on the map)"),))
    out = OUT / "buxton_dunes_by_survey.png"
    S_.save(fig, out, bbox_inches="tight", close=True)
    S_.record_caption(out, (
        "Lidar elevation of GIS 10-15, north of Buxton, in every 1996-2009 survey (rows; flight months beneath each "
        "year), the coast unrolled along NC-12: alongshore distance north along the road's 2008 alignment against "
        "cross-shore distance from it, seaward down (the ocean at the bottom), the cross-shore axis stretched so the "
        "dunes can be seen. Elevation with hillshade, every survey on the 2009 lidar's datum (the model's vertical "
        "offsets, then each survey's remaining difference from 2009 on NC-12 removed; supporting/buxton_datum_check.csv). "
        "White contour: the berm, 1.7 m NAVD88, the collision module's toe. Magenta: NC-12 where it was at each survey: "
        "the old alignment through 1998, washed out by Hurricane Dennis in the September 1999 survey (dashed), and the "
        "relocated alignment from 2000 (about 3,500 ft moved landward just north of Buxton, in place by 16 March 2000; "
        "NCDOT 2015). Grey vertical lines: the domain boundaries, where NC-12 passes from one domain into the next, "
        "numbered along the top. Survey grids of 1 m (2009) to 9 ft (about 2.7 m), resampled to 1 m. *2001 is "
        "unclassified lidar: buildings and vegetation can stand out as false highs."))
    print(f"wrote {out}")
    return st


def strip_changes(surv, st, S_, plt):
    years = sorted(surv)
    pairs = list(zip(years[:-1], years[1:]))
    fig, axs = plt.subplots(len(pairs), 1, figsize=(S_.FIG_W_DOUBLE, 0.82 * len(pairs) + 1.3), sharex=True,
                            layout="constrained")
    for ax, (a, b) in zip(axs, pairs):
        za, zb = st.sample(surv[a][0]), st.sample(surv[b][0])
        land = (np.nan_to_num(za, nan=-9) > 0.36) | (np.nan_to_num(zb, nan=-9) > 0.36)
        d = np.where(land & np.isfinite(za) & np.isfinite(zb), zb - za, np.nan)
        ext = (st.s[0] / 1000, st.s[-1] / 1000, st.n[0], st.n[-1])
        im = ax.imshow(d, extent=ext, origin="lower", aspect="auto", cmap="RdBu", vmin=-DIFF_MAX, vmax=DIFF_MAX,
                       interpolation="nearest")
        st.draw(ax, b)
        _row_label(ax, f"{a}\u2013{str(b)[2:]}{'*' if 2001 in (a, b) else ''}")
    st.top_axis(axs[0])
    axs[-1].set_xlabel("Alongshore distance along NC-12, south to north (km)", fontsize=8)
    fig.supylabel("Cross-shore distance from NC-12 (2008 alignment), seaward down (m)", fontsize=8)
    _colorbar(fig, axs, im, "Elevation change between surveys (m; red lower, blue higher)", extend="both")
    _legend(fig)
    out = OUT / "buxton_dune_change_between_surveys.png"
    S_.save(fig, out, bbox_inches="tight", close=True)
    S_.record_caption(out, (
        "Elevation change of GIS 10-15 between successive lidar surveys (rows), the coast unrolled along NC-12 as in "
        "buxton_dunes_by_survey (cross-shore stretched, seaward down), on the 2009 datum, land only (above mean high "
        "water in either survey); red lower, blue higher, clipped at +-3 m. Magenta: NC-12 where it was at the later "
        "survey (dashed: washed out by Dennis, September 1999; relocated landward by March 2000, NCDOT 2015). Grey "
        "vertical lines: the domain boundaries, numbered along the top. The survey grids differ (1 m to about 2.7 m), "
        "so changes of a few decimetres along steep dune faces are partly resampling. *2001 is unclassified lidar."))
    print(f"wrote {out}")


def changes(surv, g, dom, roads, S_, plt):
    years = sorted(surv)
    pairs = list(zip(years[:-1], years[1:]))
    t, shape, crs = g
    extent = (t.c, t.c + shape[1] * t.a, t.f + shape[0] * t.e, t.f)
    fig, axs = plt.subplots(1, len(pairs), figsize=(S_.FIG_W_DOUBLE, 8.6), layout="constrained")
    for ax, (a, b) in zip(axs, pairs):
        za, zb = surv[a][0], surv[b][0]
        land = (np.nan_to_num(za, nan=-9) > 0.36) | (np.nan_to_num(zb, nan=-9) > 0.36)
        d = np.where(land & np.isfinite(za) & np.isfinite(zb), zb - za, np.nan)
        im = ax.imshow(d, extent=extent, cmap="RdBu", vmin=-DIFF_MAX, vmax=DIFF_MAX, interpolation="nearest")
        _decorate(ax, g, dom, roads)
        ax.set_title(f"{a}→{b}{'*' if 2001 in (a, b) else ''}", fontsize=7)
    _labels(axs[0], g, dom, roads)
    fig.colorbar(im, ax=axs, orientation="horizontal", location="bottom", shrink=0.5, extend="both",
                 label="Elevation change between surveys (m)   ← lower  |  higher →")
    out = OUT / "buxton_dune_change_between_surveys.png"
    S_.save(fig, out, bbox_inches="tight", close=True)
    S_.record_caption(out, (
        "GIS 10-15: each survey minus the one before, on the 2009 datum (as buxton_dunes_by_survey), land only (above "
        "mean high water in either survey); red lower, blue higher, clipped at +-3 m. NC-12 in its 1978 (dashed) and "
        "2008 (solid) alignments. The survey grids differ (1 m to about 2.7 m), so changes of a few decimetres along "
        "steep dune faces are partly resampling. *2001 is unclassified lidar."))
    print(f"wrote {out}")


def profiles(surv, g, dom, roads, S_, plt):
    """One cross-shore profile per domain, in the strips' frame: distance from NC-12's 2008 alignment, seaward to the
    right, so the panels line up with each other and with the road positions."""
    import geopandas as gpd
    from matplotlib.lines import Line2D
    from scipy.ndimage import map_coordinates
    t, shape, crs = g
    tr = gpd.read_file(TRANSECTS).to_crs(crs)
    d = dom.to_crs(crs)
    new_road, old_road = roads["2008"].geometry.union_all(), roads["1978"].geometry.union_all()
    years = sorted(surv)
    colours = dict(zip(years, plt.get_cmap("viridis")(np.linspace(0.0, 0.9, len(years)))))
    fig, axs = plt.subplots(3, 2, figsize=(S_.FIG_W_DOUBLE, 7.4), sharex=True, sharey=True, layout="constrained")
    for k, (ax, gis) in enumerate(zip(axs.ravel(), GIS)):
        poly = d[d.domain_id == gis].geometry.iloc[0]
        inside = tr[tr.intersects(poly)]
        line = inside.iloc[(inside.geometry.centroid.y - poly.centroid.y).abs().argsort().iloc[0]].geometry
        c = np.asarray(line.coords, float)
        u = (c[-1] - c[0]) / np.hypot(*(c[-1] - c[0]))         # ocean end to land end
        s = np.arange(0, line.length, 1.0)
        xy = c[0] + s[:, None] * u
        col, row = (xy[:, 0] - t.c) / t.a - 0.5, (xy[:, 1] - t.f) / t.e - 0.5
        at = lambda geom: line.project(line.intersection(geom)) if line.intersects(geom) else np.nan
        s_new, s_old = at(new_road), at(old_road)
        x = s_new - s                                            # seaward of the 2008 road, positive
        for y in years:
            z = map_coordinates(surv[y][0], [row, col], order=1, mode="constant", cval=np.nan)
            ax.plot(x, z, color=colours[y], lw=1.3 if y in (1998, 1999, 2004) else 0.8, zorder=3)
        for _, zw, lstyle in STORMS:
            ax.axhline(zw, color="#2171B5", lw=0.7, ls=lstyle, zorder=2)
        ax.axhline(BERM_M, color="0.55", lw=0.6, zorder=2)
        ax.axvline(0, color=ROAD_C, lw=1.0, zorder=4)
        if np.isfinite(s_old) and abs(s_new - s_old) > 3:
            ax.axvline(s_new - s_old, color=ROAD_C, lw=1.0, ls=(0, (4, 2)), zorder=4)
        ax.set_xlim(-STRIP_LAND_M, STRIP_SEA_M)
        ax.set_ylim(-1, 8.5)
        ax.text(0.02, 0.96, f"({'abcdef'[k]}) GIS {gis}", transform=ax.transAxes, fontsize=8, fontweight="bold",
                va="top")
        ax.grid(color="0.93", lw=0.4)
        ax.spines[["top", "right"]].set_visible(False)
        for sp in ("left", "bottom"):
            ax.spines[sp].set_linewidth(0.5)
        ax.tick_params(labelsize=7, length=2.5, width=0.5)
    for ax in axs[:, 0]:
        ax.set_ylabel("Elevation (m NAVD88)", fontsize=8)
    for ax in axs[-1]:
        ax.set_xlabel("Distance seaward of NC-12 (2008 alignment) (m)", fontsize=8)
    h = [Line2D([], [], color=colours[y], lw=1.3 if y in (1998, 1999, 2004) else 0.8,
                label=f"{y}{'*' if y == 2001 else ''}") for y in years]
    h += [Line2D([], [], color=ROAD_C, lw=1.0, label="NC-12 from 2000"),
          Line2D([], [], color=ROAD_C, lw=1.0, ls=(0, (4, 2)), label="NC-12 to 1999"),
          Line2D([], [], color="#2171B5", lw=0.7, ls=":", label="Bonnie peak, 1998 (2.64 m)"),
          Line2D([], [], color="#2171B5", lw=0.7, ls="--", label="Dennis peak, 1999 (3.65 m)"),
          Line2D([], [], color="0.55", lw=0.6, label="Berm, 1.7 m (module toe)")]
    fig.legend(handles=h, loc="outside lower center", ncol=7, fontsize=7, frameon=False, handlelength=2.0,
               columnspacing=1.3)
    out = OUT / "buxton_dune_profiles.png"
    S_.save(fig, out, bbox_inches="tight", close=True)
    S_.record_caption(out, (
        "One cross-shore lidar profile per domain, GIS 10-15 (the 100 m transect nearest each domain's middle), in every "
        "1996-2009 survey on the 2009 datum (lighter to darker, 1996 to 2009; 1998, 1999 and 2004 thicker), against "
        "distance seaward of NC-12's 2008 alignment (the ocean to the right; the same frame as buxton_dunes_by_survey). "
        "Magenta: NC-12 from 2000 (solid, after its relocation landward just north of Buxton, in place by 16 March 2000; "
        "NCDOT 2015) and before (dashed, where the old road crosses the transect). Blue: the storm series' peak water "
        "levels in 1998 (Bonnie, 2.64 m NAVD88, dotted) and 1999 (Dennis, 3.65 m, dashed); grey: the berm, 1.7 m, the "
        "collision module's toe. *2001 is unclassified lidar."))
    print(f"wrote {out}")


def main():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from site_layer import hat_figure_style as S_
    S_.apply_style()
    OUT.mkdir(parents=True, exist_ok=True)
    surv, g, dom = surveys()
    roads = road_lines(g[2])
    datum_check(surv, g, roads)
    st = strips(surv, g, dom, roads, S_, plt)
    strip_changes(surv, st, S_, plt)
    profiles(surv, g, dom, roads, S_, plt)


if __name__ == "__main__":
    main()
