"""
How the lidar dune toes (HAT_dune_toes.py) sit against the dune lines Hannah digitised from aerial imagery.

    python scripts/real-world-cases/Hatteras/dune-collision/1-observations/6-toes-vs-dunelines/HAT_dune_toes_vs_dunelines.py

Each lidar survey is paired with the dune line traced from imagery of the same year (2023's line with the 2020
lidar, the last public survey). On every 10 m dune cell the toe, beach start and crest are picked as HAT_dune_toes.py
picks them, put back on the map, and measured against where the dune line crosses that cell's profile: positive
is the line landward of the feature. The lidar profile's elevation at the line is read too. The profiles run
east-west (the domain clips are north-up UTM), so where the shore is oblique (the cape) distances are along the
profile, not shore-normal.

Writes to the workspace's dune-collision/1-observations/6-toes-vs-dunelines/: dune_toes_vs_dunelines.png, and
supporting/dune_toes_vs_dunelines.csv (one row per survey cell) and supporting/dune_toes_vs_dunelines_summary.csv.
The toe pick itself is the input step's (input_prep/1-barrier3d-domains/3-dune-toes/HAT_dune_toes.py).

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-11-03
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

_HERE = Path(__file__).resolve()
_HAT_DIR = next(_p for _p in _HERE.parents if (_p / "hatteras_paths.py").exists())
for _path in (_HAT_DIR, _HAT_DIR / "input_prep" / "1-barrier3d-domains" / "3-dune-toes"):
    if str(_path) not in sys.path:
        sys.path.insert(0, str(_path))
import HAT_dune_toes as T  # noqa: E402

S, X = T.S, T.S.X

# --- CONFIG ------------------------------------------------------------------
DUNELINE_DIR = T.S.REPO / "data" / "hatteras_init" / "2-brie-offset" / "dunelines"
# (start whose frames are used, lidar survey year, dune line vintage, months between imagery and lidar)
PAIRS = (("1996-start", 1997, 1997, "Oct 1997 photo / Sep-Oct 1997 lidar"),
         ("1996-start", 2004, 2004, "May 2004 image / Apr-Sep 2004 lidar"),
         ("1996-start", 2009, 2009, "May 2009 image / Aug 2009 lidar"),
         ("2009-start", 2020, 2023, "2023 image / 2019-20 lidar (3 years apart)"))
OUT_DIR = T.S.REPO / "dune-collision" / "1-observations" / "6-toes-vs-dunelines"
# The imagery figure: surveys with an aerial image from the same months (label, file glob under IMAGERY_ROOT),
# and the domains drawn: GIS 22 and 85 agree, 12 is Buxton, 34 and 66 the two stretches that disagree
IMAGERY_ROOT = T.S.GIS_ROOT / "aerial_imagery"
MAP_IMAGERY = {1997: ("1997 lidar on the 12 Oct 1997 photo\n(1997 line)", "1997/usgs_henderson/frames/*.tif"),
               2020: ("2020 lidar on the Jan-Apr 2020 NOAA image\n(2023 line)",
                      "2020/noaa_2020_NGS_NaturalColorImagery_HurricaneFlorence_J1472575/*.tif")}
MAP_SITES = (12, 22, 34, 66, 85)
MAP_HALF_M = 150.0
LINE_C, TOE_C, CREST_C, START_C = "#00e5ff", "#ff2d55", "#ffd60a", "#ffffff"
ZONES = (("cape", 1, 7), ("buxton_avon_road", 8, 20), ("avon", 21, 31), ("avon_salvo_road", 32, 67),
         ("tri_village", 68, 83), ("rodanthe_north", 84, 90))
# -----------------------------------------------------------------------------


# The 10 m cell-centre map coordinates, through the same orient / flip / shear / trim as S.profiles
def frame_coords(grid, shear, c0):
    t, (h, w), _ = grid
    cols, rows = np.meshgrid(np.arange(w) + 0.5, np.arange(h) + 0.5)
    out = []
    for arr in (t.c + cols * t.a, t.f + rows * t.e):
        a = X.orient_ocean_right(S.downsample(arr), X.OCEAN_LOC)[:, ::-1]
        zs = np.full_like(a, np.nan)
        for i in range(a.shape[0]):
            k = int(shear[i])
            zs[i, :a.shape[1] - k] = a[i, k:]
        out.append(zs[:, c0:])
    return out


def duneline(vintage, crs):
    import geopandas as gpd
    return gpd.read_file(DUNELINE_DIR / f"duneline_{vintage}.geojson").to_crs(crs).geometry.union_all()


# Where the line crosses one east-west profile: the crossing nearest the cell's crest, its easting
def crossing(line, y, x0, x1, x_near):
    from shapely.geometry import LineString
    hit = line.intersection(LineString([(min(x0, x1), y), (max(x0, x1), y)]))
    xs = [p.x for p in getattr(hit, "geoms", [hit]) if not p.is_empty and p.geom_type == "Point"]
    if not xs:
        return np.nan, 0
    return min(xs, key=lambda x: abs(x - x_near)), len(xs)


def compare():
    import json
    import rasterio
    rows = []
    for start, year, vintage, _ in PAIRS:
        prod, surveys = S.STARTS[start]
        windows = json.loads(S.tv.picks_dir(start).joinpath(
            f"HAT_dune_search_windows_{S.tv.topo_dirs(start)[2]}.json").read_text(encoding="utf-8"))
        kind, path, *_ = surveys[year]
        files = sorted(path.parent.glob(path.name)) if "*" in path.name else [path]
        srcs = [rasterio.open(f) for f in files]
        clips = S.product(prod).gapfill_1m
        grids = {g: S.clip_grid(clips / f"clip_domain_{g}_filled.tif") for g in T.GIS}
        line = duneline(vintage, grids[1][2].to_wkt())
        status = pd.read_csv(S.tv.product_dir(start) / "survey-crests" / "survey_dune_crests.csv")
        status = status[status.survey.astype(str) == str(year)].set_index(["gis", "cell"]).status
        try:
            unit = "m" if kind == "mosaic" else S.unit_of(srcs, start, grids)[0]
            scale = S.US_FT if unit == "ft" else 1.0
            off = S.SURVEY_VERTICAL_OFFSET_M.get(start, {}).get(year, 0.0)
            for g in T.GIS:
                shear, c0, i0, i1 = S.frame(start, g, windows)
                raw = S.mosaic_on_clip(srcs[0], grids[g]) if kind == "mosaic" else S.raster_on_clip(srcs, grids[g])
                z = S.profiles(S.downsample(raw * scale + off), shear, c0)
                xx, yy = frame_coords(grids[g], shear, c0)
                r = T.crest_and_toe(z, i0, i1, X)
                for cell in range(T.CELLS):
                    if not np.isfinite(r["crest_idx"][cell]):
                        continue
                    ci, sb = int(r["crest_idx"][cell]), int(r["beach_start_idx"][cell])
                    y = np.nanmean(yy[cell])
                    xl, n = crossing(line, y, np.nanmax(xx[cell]), np.nanmin(xx[cell]), xx[cell, ci])
                    row = dict(survey=year, line=vintage, gis=g, cell=cell, crossings=n,
                               status=status.get((g, cell), "missing"),
                               crest_m_navd88=r["crest"][cell] + T.MHW_NAVD88_M,
                               toe_m_navd88=r["toe"][cell] + T.MHW_NAVD88_M,
                               line_minus_crest_m=xx[cell, ci] - xl, line_minus_beach_start_m=xx[cell, sb] - xl,
                               # map positions (UTM 18N m) for the imagery figure: one northing per profile
                               y=y, crest_x=xx[cell, ci], beach_start_x=xx[cell, sb], line_x=xl)
                    if np.isfinite(r["toe_idx"][cell]):
                        row["toe_x"] = xx[cell, int(r["toe_idx"][cell])]
                        row["line_minus_toe_m"] = row["toe_x"] - xl
                    # The lidar's elevation at the line: along the profile, eastings fall landward
                    ok = np.isfinite(xx[cell]) & (z[cell] > X.SENTINEL_WATER_M + 1e-9)
                    if np.isfinite(xl) and ok.sum() > 2:
                        row["z_at_line_m_navd88"] = np.interp(-xl, -xx[cell][ok], z[cell][ok]) + T.MHW_NAVD88_M
                    rows.append(row)
                print(f"{year} GIS {g}", end="\r", flush=True)
        finally:
            for s in srcs:
                s.close()
        print(f"{year} done        ", flush=True)
    t = pd.DataFrame(rows)
    # The toes HAT_dune_toes.py keeps: crest ok, toe above the floor and below the crest
    t["toe_kept"] = ((t.status == "ok") & (t.toe_m_navd88 > T.TOE_MIN_M_NAVD88)
                     & (t.toe_m_navd88 < t.crest_m_navd88 - T.TOE_BELOW_CREST_M))
    t["zone"] = pd.cut(t.gis, [z[1] - 0.5 for z in ZONES] + [90.5], labels=[z[0] for z in ZONES])
    return t


def summary(t):
    k = t[t.toe_kept & t.line_minus_toe_m.notna()]
    def stats(d):
        e = d.line_minus_toe_m
        return pd.Series(dict(cells=len(d), median_m=e.median(), iqr_lo_m=e.quantile(.25), iqr_hi_m=e.quantile(.75),
                              within_10m=(e.abs() <= 10).mean(), within_20m=(e.abs() <= 20).mean(),
                              line_minus_crest_median_m=d.line_minus_crest_m.median(),
                              line_minus_beach_start_median_m=d.line_minus_beach_start_m.median(),
                              toe_median_m_navd88=d.toe_m_navd88.median(),
                              z_at_line_median_m_navd88=d.z_at_line_m_navd88.median()))
    by_survey = k.groupby("survey").apply(stats)
    by_zone = k.groupby(["survey", "zone"], observed=True).apply(stats)
    return by_survey, by_zone


def figure(t, by_survey, out_dir):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from site_layer import hat_figure_style as fs
    fs.apply_style()
    k = t[t.toe_kept & t.line_minus_toe_m.notna()]
    fig, axes = plt.subplots(2, 1, figsize=(7.2, 6.4), constrained_layout=True)
    colors = dict(zip([p[1] for p in PAIRS], ("#b2182b", "#d6a01d", "#2166ac", "#5a5a5a")))
    for (_, year, vintage, label) in PAIRS:
        d = k[k.survey == year]
        med = d.groupby("gis").line_minus_toe_m.median()
        axes[0].plot(med.index, med.values, "-o", ms=2.5, lw=1, color=colors[year], label=f"{vintage} line vs {year} lidar")
        axes[1].hist(d.line_minus_toe_m.clip(-60, 60), bins=np.arange(-60, 61, 10), histtype="step", lw=1.5,
                     color=colors[year], label=f"{year}: median {d.line_minus_toe_m.median():+.0f} m")
    axes[0].axhline(0, color="k", lw=0.8)
    axes[0].set(xlabel="GIS domain", ylabel="Dune line minus lidar toe (m)\n+ = line landward of the toe",
                title="Where the imagery dune line sits against the lidar toe, per domain (median of cells)")
    axes[0].legend(fontsize=7, ncol=2, frameon=False)
    axes[1].axvline(0, color="k", lw=0.8)
    axes[1].set(xlabel="Dune line minus lidar toe (m), 10 m cells, clipped at ±60 m", ylabel="Dune cells",
                title="All cells with a kept toe")
    axes[1].legend(fontsize=7, frameon=False)
    png = out_dir / "dune_toes_vs_dunelines.png"
    fig.savefig(png, dpi=200)
    plt.close(fig)
    return png


# The aerial image of one survey under a square window, on the clips' CRS, RGB uint8
def image_under(files, crs, bounds, res=0.5):
    import rasterio
    from rasterio.transform import from_origin
    from rasterio.warp import reproject, transform_bounds, Resampling
    x0, y0, x1, y1 = bounds
    w, h = int((x1 - x0) / res), int((y1 - y0) / res)
    dst_t = from_origin(x0, y1, res, res)
    out = np.zeros((3, h, w), np.uint8)
    for f in files:
        with rasterio.open(f) as src:
            sb = transform_bounds(crs, src.crs, *bounds)
            if sb[2] < src.bounds.left or sb[0] > src.bounds.right or sb[3] < src.bounds.bottom or sb[1] > src.bounds.top:
                continue
            a = np.zeros_like(out)
            reproject(rasterio.band(src, [1, 2, 3]), a, src_nodata=0, dst_transform=dst_t, dst_crs=crs, dst_nodata=0,
                      resampling=Resampling.bilinear)
            empty = (out == 0).all(axis=0)
            out[:, empty] = a[:, empty]
        if not (out == 0).all(axis=0).any():
            break
    return out


def maps(t):
    """Each survey's toes, crests and beach starts on its aerial image, with the dune line, at MAP_SITES."""
    import glob
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    from shapely.geometry import box
    from site_layer import hat_figure_style as fs
    fs.apply_style()
    crs = S.clip_grid(S.product(S.STARTS["2009-start"][0]).gapfill_1m / "clip_domain_1_filled.tif")[2]
    fig, axes = plt.subplots(len(MAP_IMAGERY), len(MAP_SITES), figsize=(2.3 * len(MAP_SITES), 2.6 * len(MAP_IMAGERY) + 0.6),
                             constrained_layout=True)
    for r, (year, (label, pattern)) in enumerate(MAP_IMAGERY.items()):
        files = sorted(glob.glob(str(IMAGERY_ROOT / pattern)))
        vintage = next(p[2] for p in PAIRS if p[1] == year)
        line = duneline(vintage, crs.to_wkt())
        for c, g in enumerate(MAP_SITES):
            ax = axes[r, c]
            d = t[(t.survey == year) & (t.gis == g)]
            k = d[d.toe_kept]
            yc = d.y.median()
            xc = np.nanmedian(np.r_[k.toe_x.to_numpy(), d.crest_x.to_numpy()])
            half = MAP_HALF_M
            bounds = (xc - half, yc - half, xc + half, yc + half)
            img = image_under(files, crs, bounds)
            ax.imshow(np.moveaxis(img, 0, -1), extent=(bounds[0], bounds[2], bounds[1], bounds[3]), zorder=0)
            clipped = line.intersection(box(*bounds))
            for part in getattr(clipped, "geoms", [clipped]):
                if part.is_empty or part.geom_type != "LineString":
                    continue
                xs, ys = part.xy
                ax.plot(xs, ys, color=LINE_C, lw=1.6, zorder=3, solid_capstyle="round")
            ax.plot(d.beach_start_x, d.y, "^", ms=3, color=START_C, mec="k", mew=0.3, zorder=4)
            ax.plot(d.crest_x, d.y, "s", ms=3, color=CREST_C, mec="k", mew=0.3, zorder=4)
            ax.plot(k.toe_x, k.y, "o", ms=3.6, color=TOE_C, mec="k", mew=0.3, zorder=5)
            ax.plot(d[~d.toe_kept].toe_x, d[~d.toe_kept].y, "o", ms=3.6, mfc="none", mec=TOE_C, mew=0.8, zorder=5)
            ax.set_xlim(bounds[0], bounds[2])
            ax.set_ylim(bounds[1], bounds[3])
            ax.set_xticks([])
            ax.set_yticks([])
            ax.set_aspect("equal")
            # 50 m scale bar, bottom left
            sx, sy = bounds[0] + 0.06 * 2 * half, bounds[1] + 0.06 * 2 * half
            ax.plot([sx, sx + 50], [sy, sy], color="white", lw=2.5, zorder=6)
            ax.text(sx + 25, sy + 0.025 * 2 * half, "50 m", color="white", ha="center", va="bottom", fontsize=7, zorder=6)
            med = k.line_minus_toe_m.median()
            ax.set_title(f"GIS {g}: line {med:+.0f} m from toe" if np.isfinite(med) else f"GIS {g}", fontsize=8)
            if c == 0:
                ax.set_ylabel(label, fontsize=8)
    handles = [Line2D([], [], color=LINE_C, lw=1.6, label="dune line traced from imagery"),
               Line2D([], [], ls="none", marker="o", ms=4, color=TOE_C, mec="k", mew=0.3, label="lidar toe (kept)"),
               Line2D([], [], ls="none", marker="o", ms=4, mfc="none", mec=TOE_C, label="lidar toe (screened out)"),
               Line2D([], [], ls="none", marker="s", ms=4, color=CREST_C, mec="k", mew=0.3, label="lidar crest"),
               Line2D([], [], ls="none", marker="^", ms=4, color=START_C, mec="k", mew=0.3, label="lidar beach start")]
    fig.legend(handles=handles, loc="lower center", ncol=5, frameon=False, fontsize=7, bbox_to_anchor=(0.5, -0.04))
    png = OUT_DIR / "dune_toes_on_imagery.png"
    fig.savefig(png, dpi=220, bbox_inches="tight")
    plt.close(fig)
    return png


def main():
    import argparse
    ap = argparse.ArgumentParser(description=__doc__.strip().splitlines()[0])
    ap.add_argument("action", nargs="?", default="all", choices=["all", "compare", "maps"])
    a = ap.parse_args()
    if a.action == "maps":
        print(f"wrote {maps(pd.read_csv(OUT_DIR / 'supporting' / 'dune_toes_vs_dunelines.csv'))}")
        return
    t = compare()
    by_survey, by_zone = summary(t)
    out_dir = OUT_DIR
    out_dir.mkdir(parents=True, exist_ok=True)
    sup = out_dir / "supporting"
    sup.mkdir(parents=True, exist_ok=True)
    t.round(3).to_csv(sup / "dune_toes_vs_dunelines.csv", index=False)
    pd.concat({"all": by_survey.assign(zone="all").set_index("zone", append=True), "zone": by_zone}).round(3).to_csv(
        sup / "dune_toes_vs_dunelines_summary.csv")
    pd.set_option("display.width", 220)
    print(by_survey.round(2).to_string())
    print(by_zone[["cells", "median_m", "iqr_lo_m", "iqr_hi_m", "within_10m"]].round(2).to_string())
    print(f"wrote {figure(t, by_survey, out_dir)}")
    if a.action == "all":
        print(f"wrote {maps(t)}")


if __name__ == "__main__":
    main()
