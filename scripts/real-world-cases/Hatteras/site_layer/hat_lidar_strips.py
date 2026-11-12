"""
The 1996-2009 lidar surveys on one datum, and the coast unrolled along NC-12, for any stretch of the island.

Shared by the figure scripts that look at what the dunes did survey by survey (first
dune-collision/1-observations/8-dune-change-by-survey/). Moved here 2026-11-12 from 8-buxton-dune-change's script, so
that any stretch (a Region) can be drawn the same way.

    from site_layer import hat_lidar_strips as L
    r = L.REGIONS["rodanthe"]
    surv, g, dom = L.surveys(r)                  # {year: (array m NAVD88, flown, info)}, the grid, the domains
    st = L.Strip(r, g, L.road_lines(g[2]), dom)  # the coast unrolled along NC-12's 2008 alignment
    z = st.sample(surv[1999][0])                 # (cross-shore, alongshore)

Every survey on the 2009 lidar's datum: each survey in its own vertical unit and with the model's vertical offset
(survey_sources.csv of the 1996 start, from HAT_survey_dune_crests.SURVEY_VERTICAL_OFFSET_M: 1996 -0.24, 1997 -0.34,
1999 -0.32 m), then minus what is left of its difference from the 2009 lidar on NC-12 (road_offset_vs_2009_m, the
island's road pixels outside GIS 9-13, where the road moved), so every survey sits within about 0.02-0.05 m of 2009 on
the road (datum_check). That last step goes beyond the model, which leaves 2000 (+0.10 m) and 2001 (-0.10 m) as
measured. All resampled bilinearly onto the model clips' 1 m grid; the surveys are 9 ft (about 2.7 m; 1996, 1998, 1999,
2004), 6 ft (2005), 3 m (1997, 2000), 4 m (2001) and 1 m (2009) grids, so narrow dunes look smoother in the coarse
years. 2001 is unclassified lidar: buildings and vegetation can stand out as false highs; the 1996-2004 surveys also
show the houses in the villages.

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-11-12
"""
from __future__ import annotations

import sys
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd

_HAT_DIR = next(p for p in Path(__file__).resolve().parents if (p / "hatteras_paths.py").exists())
for _p in (_HAT_DIR, _HAT_DIR / "input_prep" / "1-barrier3d-domains" / "2-survey-crests"):
    if str(_p) not in sys.path:
        sys.path.insert(0, str(_p))
from hatteras_paths import WORKSPACE  # noqa: E402
import HAT_survey_dune_crests as S  # noqa: E402   the survey registry, the clips and the reprojection

# --- CONFIG ------------------------------------------------------------------
START = "1996-start"
RES_M = 1.0                                     # the common grid: the model clips' 1 m lattice
PAD_M = 150.0                                   # around the domains
STRIP_STEP_M = 1.0
BASELINE_SIGMA_M = 250.0                        # the strips' baseline: NC-12 (2008) smoothed, Gaussian sigma (~1 km)
BASELINE_NOTE = ("the coast unrolled along a baseline, NC-12's 2008 alignment smoothed over about 1 km (Gaussian, "
                 "sigma 250 m), so a curving road does not bend the coast: alongshore distance north along the baseline "
                 "against cross-shore distance from it, seaward down (the ocean at the bottom), the cross-shore axis "
                 "stretched so the dunes can be seen")
STRIP_YLABEL = "Cross-shore distance from the baseline (NC-12 smoothed), seaward down (m)"
STRIP_XLABEL = "Alongshore distance along the baseline, south to north (km)"
MHW_M = 0.36                                    # mean high water, m NAVD88
BERM_M = 1.7                                    # Barrier3D's berm, the module's toe [m NAVD88]
STORMS = (("Bonnie 1998 peak", 2.64, ":"), ("Dennis 1999 peak", 3.65, "--"))   # the storm series' peak water levels
MONTHS = {1996: "Aug-Oct", 1997: "Sep-Oct", 1998: "Sep", 1999: "Sep", 2000: "Feb-Mar", 2001: "Jan-Mar",
          2004: "Apr-Sep", 2005: "Aug-Nov", 2009: "Aug"}   # the flights (survey_sources.csv)
IMAGERY_YEARS = (1993, 1995, 1996, 1997, 1998, 2002, 2004, 2008, 2009)    # the photographs on D: around 1996-2009
IMAGERY_LABELS = {1993: "2 Sep", 1995: "18 Oct", 1996: "14 Oct", 1997: "12 Oct", 1998: "2 Oct", 2002: "18 Oct",
                  2004: "25 May", 2008: "26-27 Mar", 2009: "30 May"}
ROAD_C, BOUNDARY_C = "#AE017E", "0.45"           # NC-12 (magenta), the domain boundaries (grey)
TERRAIN_RANGE = (-1.0, 9.0)                     # m NAVD88, the terrain colour scale (blue below about 0)
DIFF_MAX = 3.0                                  # m, the change colour scale
SOURCES = S.tv.product_dir(START) / "survey-crests" / "survey_sources.csv"
TRANSECTS = WORKSPACE / "data" / "hatteras_init" / "2-brie-offset" / "transects" / "transects_100m.geojson"
ROADS = {"1978": S.tv.MGMT_ROOT / "road_offset" / "raw_offset" / "1978",
         "2008": S.tv.MGMT_ROOT / "road_offset" / "raw_offset" / "2008"}
# -----------------------------------------------------------------------------


@dataclass(frozen=True)
class Region:
    """A stretch of the island to draw. relocated_by_survey: the first survey with NC-12 on its 2008 (relocated)
    alignment, earlier ones drawn on the 1978 line (None: the road did not move in 1996-2009); washed_out_survey: the
    survey in which the old road is drawn dashed, washed out."""
    name: str
    label: str
    gis: range
    road_check_gis: tuple
    sea_m: float = 180.0                        # the strips reach this far seaward of NC-12 (2008)
    land_m: float = 60.0                        # ... and this far landward
    relocated_by_survey: int | None = None
    washed_out_survey: int | None = None
    note: str = ""                              # a sentence for the captions


# NC-12's relocation in GIS 9-14 (hatteras_site_config RelocationEvent, 1999): after Hurricane Dennis (late August to
# early September 1999), about 3,500 ft just north of Buxton, in place in an NCDOT photo of 16 March 2000 (NC 12 Buxton
# to Avon Feasibility Study, NCDOT, December 2015, section 2.1 and Figure 4; the management timeline's 1999 row). The
# September 1999 lidar, flown just after Dennis, still has the old alignment (washed out); the February-March 2000
# survey is the first with the new one.
REGIONS = {
    "buxton": Region("buxton", "north of Buxton", range(10, 16), (15,), relocated_by_survey=2000,
                     washed_out_survey=1999,
                     note=("NC-12 where it was at each survey: the old alignment through 1998, washed out by Hurricane "
                           "Dennis in the September 1999 survey (dashed), and the relocated alignment from 2000 (about "
                           "3,500 ft moved landward just north of Buxton, in place by 16 March 2000; NCDOT 2015).")),
    # Avon's dune stands 150-250 m out, behind the oceanfront houses
    "avon_south": Region("avon_south", "Avon south", range(21, 27), (21, 22, 23, 24, 25, 26), sea_m=280.0,
                         note=("NC-12 on its 2008 alignment (it did not move here in 1996-2009). The 1996-2004 "
                               "surveys show the houses between the road and the dune.")),
    # GIS 88-89's dune stands near 180 m out
    "rodanthe": Region("rodanthe", "Rodanthe", range(82, 90), (82, 83, 84, 85, 86, 87, 88, 89), sea_m=240.0,
                       note="NC-12 on its 2008 alignment (it did not move here in 1996-2009)."),
}


def domains(region):
    import geopandas as gpd
    d = gpd.read_file(S.GIS_ROOT / "domains.geojson")
    return d[d.domain_id.isin(region.gis)]


def grid(region):
    """One 1 m grid over the model's own domain clips of the region (their CRS and pixel lattice: the grid
    S.raster_on_clip expects), padded by PAD_M."""
    from affine import Affine
    clips = S.product(S.STARTS[START][0]).gapfill_1m
    grids = [S.clip_grid(clips / f"clip_domain_{g}_filled.tif") for g in region.gis]
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


def surveys(region):
    """Each survey on the 2009 datum, on the region's grid: ({year: (array m NAVD88, flown, info)}, grid, domains)."""
    import rasterio
    src_tab = pd.read_csv(SOURCES).set_index("survey")
    entries = S.STARTS[START][1]
    dom = domains(region)
    g = grid(region)
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


def datum_check(region, surv, g, roads):
    """Each survey's median difference from 2009 on NC-12's 2008 line in the region's road-check domains, after the
    corrections: near zero for every survey if the datum is right."""
    from rasterio.features import rasterize
    t, shape, crs = g
    dom = domains(region).to_crs(crs)
    keep = rasterize([(geom, 1) for geom in dom[dom.domain_id.isin(region.road_check_gis)].geometry],
                     out_shape=shape, transform=t, fill=0).astype(bool)
    road = rasterize([(geom.buffer(4), 1) for geom in roads["2008"].geometry], out_shape=shape, transform=t,
                     fill=0).astype(bool) & keep
    ref = surv[2009][0]
    rows = []
    for year, (z, _, info) in surv.items():
        d = (z - ref)[road & np.isfinite(z) & np.isfinite(ref)]
        rows.append(info | dict(region=region.name, road_check_gis=" ".join(map(str, region.road_check_gis)),
                                road_px=int(d.size), minus_2009_on_road_m=float(np.median(d)) if d.size else np.nan))
    return pd.DataFrame(rows)


class Strip:
    """The coast unrolled along a baseline, NC-12's 2008 alignment smoothed over about 1 km (BASELINE_SIGMA_M):
    alongshore distance north from the region's south end against cross-shore distance from the baseline, seaward
    positive, sampled every STRIP_STEP_M; drawn with the cross-shore axis stretched so the dunes, tens of metres wide over
    a few km of coast, can be seen. NC-12 itself is drawn at its true offset from the baseline (road08, road78).

    Why smoothed (2026-11-12): unrolled along the raw road, a road that curves (Rodanthe's S-curves turn about 40
    degrees in 4 km) is straightened and everything around it bent with it, so the shoreline and the house rows came out
    wavy and fanned. The smoothed baseline follows the coast's trend, not the road's bends."""

    def __init__(self, region, g, roads, dom):
        from scipy.ndimage import gaussian_filter1d
        from shapely.geometry import box, LineString
        from shapely.ops import linemerge
        self.region = region
        t, (h, w), crs = g
        road = roads["2008"].geometry.intersection(box(t.c, t.f + h * t.e, t.c + w * t.a, t.f)).union_all()
        if road.geom_type == "MultiLineString":
            road = linemerge(road)
        if road.geom_type != "LineString":
            road = max(road.geoms, key=lambda x: x.length)
        c = np.asarray(road.coords)
        if c[0, 1] > c[-1, 1]:                                        # run south to north
            road = LineString(c[::-1])
        step = 5.0
        q = np.array([road.interpolate(v).coords[0] for v in np.arange(0, road.length, step)])
        sm = np.column_stack([gaussian_filter1d(q[:, k], BASELINE_SIGMA_M / step, mode="nearest") for k in (0, 1)])
        line = self.baseline = LineString(sm)
        d = dom.to_crs(crs)
        s0 = line.project(d[d.domain_id == min(region.gis)].geometry.iloc[0].centroid) - 250
        s1 = line.project(d[d.domain_id == max(region.gis)].geometry.iloc[0].centroid) + 250
        self.s = np.arange(max(s0, 0), min(s1, line.length), STRIP_STEP_M)
        self.n = np.arange(-region.land_m, region.sea_m + STRIP_STEP_M, STRIP_STEP_M)
        p = np.array([line.interpolate(v).coords[0] for v in self.s])
        tan = np.gradient(p, axis=0)
        tan /= np.hypot(tan[:, 0], tan[:, 1])[:, None]
        nrm = np.column_stack([tan[:, 1], -tan[:, 0]])               # right of north-going: east, the ocean
        nrm *= np.sign(nrm[:, 0])[:, None]
        xy = p[None, :, :] + self.n[:, None, None] * nrm[None, :, :]
        self.col, self.row = (xy[..., 0] - t.c) / t.a - 0.5, (xy[..., 1] - t.f) / t.e - 0.5
        # each domain's stretch of the road: where NC-12 runs through the domain's box
        self.edges, self.mids, self.spans = [], {}, {}
        for _, r in d.iterrows():
            part = line.intersection(r.geometry)
            pts = [c for p_ in getattr(part, "geoms", [part]) if not p_.is_empty for c in p_.coords]
            if not pts:
                continue
            ss = [line.project(type(r.geometry.centroid)(*c[:2])) for c in pts]
            self.spans[int(r.domain_id)] = (min(ss), max(ss))
            self.mids[int(r.domain_id)] = (min(ss) + max(ss)) / 2
            self.edges += [min(ss), max(ss)]
        self._p, self._nrm = p, nrm
        self.road08 = self.offset(roads["2008"].geometry.union_all())     # NC-12 (2008), seaward of the baseline (m)
        self.road78 = None
        if region.relocated_by_survey is not None and "1978" in roads:
            self.road78 = self.offset(roads["1978"].geometry.union_all())

    def offset(self, geom, reach=400.0):
        """Where a line (a road) crosses each cross-shore normal of the baseline: its distance seaward of the baseline
        (m, + seaward), NaN where it does not cross within +-reach."""
        from shapely.geometry import LineString
        out = np.full(len(self.s), np.nan)
        for i, (v, k) in enumerate(zip(self._p, self._nrm)):
            x = LineString([v - reach * k, v + reach * k]).intersection(geom)
            if x.is_empty:
                continue
            pts = [np.asarray(g_.coords[0]) for g_ in getattr(x, "geoms", [x])]
            dd = [float(np.dot(pt - v, k)) for pt in pts]
            out[i] = dd[int(np.argmin(np.abs(dd)))]
        return out

    def extent(self):
        return self.s[0] / 1000, self.s[-1] / 1000, self.n[0], self.n[-1]

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
        r = self.region
        ext = self.extent()
        old = self.road78 is not None and year < r.relocated_by_survey
        y = self.road78 if old else self.road08
        ax.plot(self.s / 1000, y, color=ROAD_C, lw=1.0, ls=(0, (4, 2)) if year == r.washed_out_survey else "-",
                solid_capstyle="butt", zorder=8)
        for e in self.boundaries()[1:-1]:
            ax.axvline(e / 1000, color=BOUNDARY_C, lw=0.6, zorder=6)
        ax.set_xlim(ext[0], ext[1])
        ax.set_ylim(ext[3], ext[2])                            # seaward down: the ocean at the bottom
        ax.set_yticks([0, 100, 200] if r.sea_m > 200 else [0, 100])
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


def row_label(ax, text, sub=None):
    ax.set_ylabel(text if not sub else f"{text}\n{sub}", rotation=0, ha="right", va="center", fontsize=8,
                  labelpad=8, linespacing=1.1)


def colorbar(fig, axs, mappable, label, extend="neither"):
    cb = fig.colorbar(mappable, ax=axs, orientation="horizontal", location="bottom", shrink=0.45, aspect=40,
                      pad=0.012, extend=extend)
    cb.set_label(label, fontsize=8)
    cb.ax.tick_params(labelsize=7, width=0.5, length=2.5)
    cb.outline.set_linewidth(0.5)


def imagery(year, cache_dir):
    """One year's photographs as site_layer.hat_aerial_imagery.Imagery, with two additions: 2002, which exists only as
    the USGS delivery's zipped frames (unzipped once into cache_dir, the drive untouched), and 2009, Hannah's Google
    Earth captures, of which only those with a coordinate system are used."""
    import glob
    import zipfile
    import rasterio
    from rasterio.warp import transform_bounds
    from shapely.geometry import box
    from site_layer import hat_aerial_imagery as AI
    files = AI.year_files(year)
    date = AI.year_date(year)
    if not files and year == 2002:
        cache_dir.mkdir(parents=True, exist_ok=True)
        for z in sorted(glob.glob(str(AI.AERIAL_ROOT / "_delivered" / "usgs_henderson_2025_draft" /
                                      "Hatteras_Georeferenced_DRAFT" / "2002" / "2002_georef_TIF" / "*.zip"))):
            with zipfile.ZipFile(z) as zz:
                for n in zz.namelist():
                    out = cache_dir / Path(n).name
                    if not out.exists():
                        out.write_bytes(zz.read(n))
        files, date = sorted(cache_dir.glob("2002_*.tif")), "2002-10-18"
    if not files and year == 2009:
        for f in sorted(glob.glob(str(AI.AERIAL_ROOT / "2009" / "google_earth_captures" / "raw_ge" / "*.jpg"))):
            with rasterio.open(f) as src:
                if src.crs is not None:
                    files.append(Path(f))
        date = "2009-05-30"
    im = AI.Imagery.__new__(AI.Imagery)
    im.year, im.date, im.files, im.foot, im.has_mosaic = year, date, files, [], False
    for f in files:
        with rasterio.open(f) as src:
            b = transform_bounds(src.crs, "EPSG:3725", *src.bounds)
            unit = 0.3048 if "foot" in str(src.crs).lower() or "us_survey_feet" in str(src.crs).lower() else 1.0
            im.foot.append((box(*b), src.res[0] * unit))
    im.has_mosaic = bool(files) and files[0].name.endswith("_full_aerial.tif")
    return im


def imagery_rows(st, g, cache_dir, years=IMAGERY_YEARS, min_cover=0.05):
    """The photographs of each year sampled into the strip: [(year, date, (cross-shore, alongshore, 3) in 0-1)],
    leaving out years with no photograph over at least min_cover of the strip."""
    from scipy.ndimage import map_coordinates
    t, (h, w), crs = g
    x0, x1, y0, y1 = t.c, t.c + w * t.a, t.f + h * t.e, t.f
    rows = []
    for y in years:
        im = imagery(y, cache_dir)
        if not im.files:
            continue
        rgb = im.read(x0, x1, y0, y1, RES_M, crs)                        # (H, W, 3) uint8, 0 where no photograph
        strip = np.stack([map_coordinates(rgb[..., c].astype(float), [st.row, st.col], order=1, mode="constant",
                                          cval=0) for c in range(3)], axis=-1)
        if (strip.sum(axis=-1) > 0).mean() < min_cover:
            continue
        rows.append((y, im.date, np.clip(strip / 255.0, 0, 1)))
    return rows
