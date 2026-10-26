"""
The CoastSat shoreline through each fill, month by month, at its true map position on imagery.

    python scripts/input_prep/4-mgmt-forcings/nourishment_animation_map_coastsat.py

One GIF per project the hindcast fires. Left: the whole footprint +/- 2
domains at true scale, the band between the pre-fill line and this month's
line shaded by direction. Right: a true-scale zoom on the transect with the
largest gain, where the two lines separate. Frames and windows are those of
nourishment_animation_coastsat.py.

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-10-04
"""
from __future__ import annotations

import sys
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd

_HERE = Path(__file__).resolve()
sys.path.insert(0, str(_HERE.parent))
import nourishment_animation_coastsat as AN  # noqa: E402
import nourishment_extent_coastsat as X  # noqa: E402
from nourishment_shorelines_coastsat import transect_geometry  # noqa: E402

import matplotlib  # noqa: E402
matplotlib.use("Agg")
import matplotlib.patheffects as pe  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.animation import FuncAnimation, PillowWriter  # noqa: E402
from matplotlib.collections import PolyCollection  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402
from matplotlib.patches import Patch  # noqa: E402

from site_layer.hat_figure_style import C_1984, C_1997, INK, apply_style, record_caption, save  # noqa: E402
from site_layer.hat_observed_rates import DOMAIN_BOXES  # noqa: E402
from site_layer.hatteras_site_config import HATTERAS_NOURISHMENT_PROJECTS  # noqa: E402

# --- CONFIG ------------------------------------------------------------------
OUT_DIR = AN.OUT_DIR
MAP_CRS = "EPSG:26918"
PAD_DOMAINS = 2             # domains drawn either side of the footprint in the overview
ZOOM_TALL_M = 600           # height of the true-scale zoom
OVERVIEW_ZOOM = 14          # tile zoom levels
DETAIL_ZOOM = 17
GAIN_FILL = "#4fa3e0"       # band seaward of the pre-fill line
LOSS_FILL = "#ef6f5e"       # band landward of it
TILE_CACHE = Path(tempfile.gettempdir()) / "hat_tile_cache"
DARK = [pe.withStroke(linewidth=2.2, foreground="black")]
# -----------------------------------------------------------------------------


# Imagery for a box, warped to MAP_CRS once
def basemap(bx, zoom):
    import contextily as cx
    from rasterio.warp import transform_bounds
    cx.set_cache_dir(str(TILE_CACHE))
    w, s, e, n = transform_bounds(MAP_CRS, "EPSG:3857", *bx)
    img, ext = cx.bounds2img(w, s, e, n, zoom=zoom, source=cx.providers.Esri.WorldImagery, ll=False)
    return cx.warp_tiles(img, ext, t_crs=MAP_CRS)


# A north-up box of the given aspect (height/width) around a set of points, padded
def box_around(e, n, aspect, pad=1.08):
    cx_, cy_ = (np.nanmin(e) + np.nanmax(e)) / 2, (np.nanmin(n) + np.nanmax(n)) / 2
    half_w = max(np.nanmax(e) - np.nanmin(e), (np.nanmax(n) - np.nanmin(n)) / aspect) / 2 * pad
    return (cx_ - half_w, cy_ - half_w * aspect, cx_ + half_w, cy_ + half_w * aspect)


# Quads between consecutive transects, and whether each is a gain
def band(ref_e, ref_n, now_e, now_n, change):
    polys, gain = [], []
    for j in range(len(ref_e) - 1):
        q = np.array([[ref_e[j], ref_n[j]], [ref_e[j + 1], ref_n[j + 1]],
                      [now_e[j + 1], now_n[j + 1]], [now_e[j], now_n[j]]])
        if np.isfinite(q).all():
            polys.append(q)
            gain.append(np.nanmean(change[j:j + 2]) >= 0)
    return polys, gain


def scalebar(ax, bx, length, label):
    sx, sy = bx[0] + 0.07 * (bx[2] - bx[0]), bx[1] + 0.04 * (bx[3] - bx[1])
    ax.plot([sx, sx + length], [sy, sy], color="white", lw=3, solid_capstyle="butt", zorder=8)
    ax.text(sx + length / 2, sy + 0.012 * (bx[3] - bx[1]), label, color="white", ha="center", va="bottom",
            fontsize=7.5, zorder=8, path_effects=DARK)


def north(ax):
    ax.annotate("N", xy=(0.9, 0.97), xytext=(0.9, 0.90), xycoords="axes fraction", ha="center", va="center",
                color="white", fontsize=8.5, fontweight="bold", zorder=8, path_effects=DARK,
                arrowprops=dict(arrowstyle="-|>", color="white", lw=1.2))


def setup_axis(ax, bx, img, ext):
    ax.imshow(img, extent=ext, zorder=0, interpolation="bilinear")
    ax.set_xlim(bx[0], bx[2])
    ax.set_ylim(bx[1], bx[3])
    ax.set_aspect("equal")
    ax.set_xticks([])
    ax.set_yticks([])


def main():
    import geopandas as gpd
    apply_style()
    frame_tbl = X.transect_frame()
    dates = pd.read_csv(X.DATES_CSV).set_index("project")
    filled = {g for q in HATTERAS_NOURISHMENT_PROJECTS for g in q.gis_domains}
    geom = transect_geometry(set(frame_tbl.index))
    domains = gpd.read_file(DOMAIN_BOXES).to_crs(MAP_CRS).set_index("domain_id").sort_index()
    for p in sorted(HATTERAS_NOURISHMENT_PROJECTS, key=lambda q: (q.year, min(q.gis_domains))):
        d = dates.loc[p.name]
        A = AN.build(p, d, frame_tbl, filled)
        first, last = min(p.gis_domains), max(p.gis_domains)
        d_lo, d_hi = max(domains.index.min(), first - PAD_DOMAINS), min(domains.index.max(), last + PAD_DOMAINS)
        t = A["t"].join(pd.DataFrame(geom, index=["x0", "y0", "ux", "uy"]).T, how="inner")
        keep = t.domain_number.between(d_lo, d_hi).to_numpy()
        cols = [list(A["t"].index).index(tid) for tid in t.index[keep]]
        t = t[keep]
        x0, y0, ux, uy = (t[c].to_numpy() for c in ("x0", "y0", "ux", "uy"))
        ref = t.before_m.to_numpy()
        pos = A["change"][:, cols] + ref[None, :]
        ref_e, ref_n = x0 + ref * ux, y0 + ref * uy
        now_e, now_n = x0[None, :] + pos * ux[None, :], y0[None, :] + pos * uy[None, :]
        change = A["change"][:, cols]
        frames = A["frames"]

        # Overview box over the footprint +/- PAD_DOMAINS; zoom on the peak running-median gain inside the footprint
        dbx = domains.loc[d_lo:d_hi].total_bounds
        bx_o = box_around(np.r_[dbx[0], dbx[2]], np.r_[dbx[1], dbx[3]], aspect=2.2, pad=1.02)
        k_end = int(np.searchsorted(frames, A["end"] + pd.DateOffset(months=AN.STILL_MONTHS_AFTER)))
        k_end = min(k_end, len(frames) - 1)
        sm = pd.Series(np.nanmean(change[k_end:], axis=0)).rolling(X.RUN_TRANSECTS, center=True, min_periods=3).median()
        sm[~t.domain_number.between(first, last).to_numpy()] = np.nan
        j = int(np.nanargmax(sm.to_numpy()))
        cx_, cy_ = ref_e[j] + 0.5 * np.nanmean(change[k_end:, j]) * ux[j], ref_n[j] + 0.5 * np.nanmean(change[k_end:, j]) * uy[j]
        bx_z = (cx_ - ZOOM_TALL_M / 2 / 2.2, cy_ - ZOOM_TALL_M / 2, cx_ + ZOOM_TALL_M / 2 / 2.2, cy_ + ZOOM_TALL_M / 2)

        fig = plt.figure(figsize=(6.4, 7.0))
        ax_o = fig.add_axes([0.02, 0.12, 0.46, 0.78])
        ax_z = fig.add_axes([0.52, 0.12, 0.46, 0.78])
        setup_axis(ax_o, bx_o, *basemap(bx_o, OVERVIEW_ZOOM))
        setup_axis(ax_z, bx_z, *basemap(bx_z, DETAIL_ZOOM))
        for ax in (ax_o, ax_z):
            domains.loc[d_lo:d_hi].plot(ax=ax, facecolor="none", edgecolor="white",
                                                                      lw=0.5, alpha=0.6, zorder=2)
            domains.loc[first:last].plot(ax=ax, facecolor="none", edgecolor="#ffd23f", lw=1.3, zorder=2)
            ax.plot(ref_e, ref_n, color=C_1984, lw=1.4 if ax is ax_z else 0.9, zorder=4)
            north(ax)
        for dd in range(d_lo, d_hi + 1):
            c = domains.loc[dd].geometry.centroid
            ax_o.text(c.x, c.y, str(dd), color="white", fontsize=6.5, ha="center", va="center", zorder=3,
                      fontweight="bold" if first <= dd <= last else "normal", path_effects=DARK)
        ax_o.add_patch(plt.Rectangle((bx_z[0], bx_z[1]), bx_z[2] - bx_z[0], bx_z[3] - bx_z[1], fill=False,
                                     edgecolor="white", lw=1.0, ls="--", zorder=6))
        scalebar(ax_o, bx_o, 1000, "1 km")
        scalebar(ax_z, bx_z, 100, "100 m")
        ax_o.text(0.04, 0.97, f"GIS {d_lo}–{d_hi}", transform=ax_o.transAxes,
                  color="white", fontsize=8, va="top", fontweight="bold", path_effects=DARK, zorder=8)
        ax_z.text(0.04, 0.97, f"Zoom, GIS {int(t.domain_number.iloc[j])}", transform=ax_z.transAxes,
                  color="white", fontsize=8, va="top", fontweight="bold", path_effects=DARK, zorder=8)
        col_o = PolyCollection([], lw=0, alpha=0.75, zorder=3)
        col_z = PolyCollection([], lw=0, alpha=0.45, zorder=3)
        ax_o.add_collection(col_o)
        ax_z.add_collection(col_z)
        (now_o,) = ax_o.plot([], [], color=INK, lw=0.7, zorder=5)
        (now_z,) = ax_z.plot([], [], color=INK, lw=1.5, zorder=5)
        fig.text(0.02, 0.975, f"{p.name} {p.year}: model GIS {first}–{last}", fontsize=10, color=INK, va="top")
        fig.text(0.02, 0.945, f"placed {d.start_date} to {d.end_date}", fontsize=8.5, color=INK, va="top")
        date_txt = fig.text(0.98, 0.975, "", fontsize=12, fontweight="bold", color=INK, ha="right", va="top")
        stat_txt = fig.text(0.98, 0.94, "", fontsize=9, fontweight="bold", ha="right", va="top")
        handles = [Line2D([], [], color=C_1984, lw=1.4, label="Pre-fill: median of the 12 months before"),
                   Line2D([], [], color=INK, lw=1.4, label="This month: 3-month rolling median"),
                   Patch(color=GAIN_FILL, label="Seaward of pre-fill"), Patch(color=LOSS_FILL, label="Landward"),
                   Patch(facecolor="none", edgecolor="#ffd23f", lw=1.3, label="Model footprint domains")]
        fig.legend(handles=handles, loc="lower center", ncol=2, frameon=False, bbox_to_anchor=(0.5, 0.0),
                   fontsize=7.5)

        def draw(i):
            polys, gain = band(ref_e, ref_n, now_e[i], now_n[i], change[i])
            colours = [GAIN_FILL if g else LOSS_FILL for g in gain]
            for col in (col_o, col_z):
                col.set_verts(polys)
                col.set_facecolors(colours)
            now_o.set_data(now_e[i], now_n[i])
            now_z.set_data(now_e[i], now_n[i])
            date_txt.set_text(frames[i].strftime("%b %Y"))
            s, c = AN.status(frames[i], A["start"], A["end"])
            stat_txt.set_text(s)
            stat_txt.set_color(c if c != C_1997 else "#2166ac")
            return [col_o, col_z, now_o, now_z, date_txt, stat_txt]

        stem = f"nourishment_animation_map_{p.year}_{p.name.split()[0].lower()}"
        gif = OUT_DIR / f"{stem}.gif"
        FuncAnimation(fig, draw, frames=len(frames), blit=False).save(gif, writer=PillowWriter(fps=AN.FPS),
                                                                       dpi=AN.DPI)
        draw(k_end)
        still = save(fig, OUT_DIR / f"{stem}_still", close=True)[0]
        record_caption(still, (
            f"Still from {gif.name} ({len(frames)} monthly frames, {frames[0]:%b %Y} to {frames[-1]:%b %Y}, "
            f"{AN.FPS} frames per second). {p.name}, {p.year}, in CoastSat at true map position (UTM 18N, north up) "
            "on current Esri World Imagery. Red is the median shoreline of the 12 months before placement; black "
            f"is the frame's median over ±{AN.HALF_WINDOW_DAYS} days. The band between them is blue where the "
            "frame's line is seaward, red where it is landward. Left: the footprint ± "
            f"{PAD_DOMAINS} domains, model footprint domains outlined in yellow, dashed box the zoom. Right: a "
            f"{ZOOM_TALL_M} m tall true-scale zoom centred on the transect with the largest gain in the months "
            "after placement. A shift of 20-40 m is a few pixels in the overview; the band shows where it is."))
        print(f"{gif.name}: {len(frames)} frames, {gif.stat().st_size / 1e6:.1f} MB, {len(t)} transects")


if __name__ == "__main__":
    main()
