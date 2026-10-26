"""
The CoastSat shoreline the year before each fill and the year after, drawn as shorelines and as a difference.

    python scripts/input_prep/4-mgmt-forcings/nourishment_shorelines_coastsat.py   # after nourishment_extent_coastsat.py

One figure per project: (a) both shorelines straightened along a fitted
baseline with the cross-shore axis stretched, gain and loss shaded; (b) the
per-transect change; (c) both lines on imagery at true scale where the gain
peaks. Windows, dates and the background are those of nourishment_extent_coastsat.py.

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-10-04
"""
from __future__ import annotations

import json
import sys
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd
import pyproj

_HERE = Path(__file__).resolve()
sys.path.insert(0, str(_HERE.parent))
import nourishment_extent_coastsat as X  # noqa: E402

import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402
from matplotlib.patches import Patch  # noqa: E402

from site_layer import hat_observed_rates as obs  # noqa: E402
from site_layer.hat_figure_style import (C, C_1984, C_1984_FILL, C_1997, C_1997_FILL, DOMAIN_AXIS_LABEL,  # noqa: E402
                                         INK, INK_MUTED, _title, apply_style, open_frame, record_caption,
                                         save, structures)
from site_layer.hatteras_site_config import HATTERAS_NOURISHMENT_PROJECTS  # noqa: E402

# --- CONFIG ------------------------------------------------------------------
OUT_DIR = X.OUT_DIR / "shorelines"
MAP_CRS = "EPSG:26918"
BASELINE_DEG = 2            # polynomial order of the baseline fitted through the before line
ZOOM_HALF_M = 450           # half-width of the true-scale imagery zoom
TILE_ZOOM = 17
TILE_CACHE = Path(tempfile.gettempdir()) / "hat_tile_cache"
# -----------------------------------------------------------------------------


# Origin and seaward unit vector per transect, in MAP_CRS
def transect_geometry(wanted):
    to_utm = pyproj.Transformer.from_crs("EPSG:4326", MAP_CRS, always_xy=True)
    out = {}
    with open(obs.TRANSECT_LAYER, encoding="utf-8") as fh:
        layer = json.load(fh)
    for feat in layer["features"]:
        tid = str(feat["properties"]["id"]).replace("-", "_")
        if tid not in wanted:
            continue
        (lon0, lat0), (lon1, lat1) = feat["geometry"]["coordinates"][0], feat["geometry"]["coordinates"][-1]
        x0, y0 = to_utm.transform(lon0, lat0)
        x1, y1 = to_utm.transform(lon1, lat1)
        n = float(np.hypot(x1 - x0, y1 - y0))
        out[tid] = (x0, y0, (x1 - x0) / n, (y1 - y0) / n)
    return out


# Map points of both lines, and each line's offset from a smooth baseline through the before line
def add_map_coords(t, geom):
    g = pd.DataFrame(geom, index=["x0", "y0", "ux", "uy"]).T
    t = t.join(g, how="inner")
    for k in ("before", "after"):
        t[f"{k}_e"] = t.x0 + t[f"{k}_m"] * t.ux
        t[f"{k}_n"] = t.y0 + t[f"{k}_m"] * t.uy
    ok = t.before_m.notna()
    # baseline: offset along the transect as a smooth function of alongshore position
    coef = np.polyfit(t.x[ok], t.before_m[ok], BASELINE_DEG)
    base = np.polyval(coef, t.x)
    t["before_off"] = t.before_m - base
    t["after_off"] = t.after_m - base
    return t


def draw(p, d, t, summary_row, windows):
    first, last = min(p.gis_domains), max(p.gis_domains)
    filled = {g for q in HATTERAS_NOURISHMENT_PROJECTS for g in q.gis_domains}
    ctrl = t[~t.domain_number.isin(filled)].change_m.dropna()
    bg = float(ctrl.median())
    b0, b1, a0, a1 = windows
    win_b = f"{b0.date()} to {(b1 - pd.Timedelta(days=1)).date()}"
    win_a = f"{a0.date()} to {(a1 - pd.Timedelta(days=1)).date()}"

    fig = plt.figure(figsize=(7.2, 5.6))
    ax_a = fig.add_axes([0.08, 0.56, 0.58, 0.36])
    ax_b = fig.add_axes([0.08, 0.17, 0.58, 0.30], sharex=ax_a)
    ax_c = fig.add_axes([0.70, 0.17, 0.29, 0.75])

    # (a) the two shorelines, straightened, ocean up
    x = t.x.to_numpy()
    bo, ao = t.before_off.to_numpy(), t.after_off.to_numpy()
    ax_a.axvspan(first - 0.5, last + 0.5, color="0.92", lw=0, zorder=0)
    ax_a.fill_between(x, bo, ao, where=ao >= bo, color=C_1997_FILL, lw=0, interpolate=True, zorder=1)
    ax_a.fill_between(x, bo, ao, where=ao < bo, color=C_1984_FILL, lw=0, interpolate=True, zorder=1)
    ax_a.plot(x, bo, color=C_1984, lw=1.1, zorder=3)
    ax_a.plot(x, ao, color=C_1997, lw=1.1, zorder=4)
    ax_a.set_ylabel("Shoreline position\nfrom a fitted baseline (m)")
    ax_a.text(0.01, 0.96, "ocean ↑", transform=ax_a.transAxes, fontsize=8, color=INK_MUTED, va="top")
    ax_a.text(0.01, 0.04, "land ↓", transform=ax_a.transAxes, fontsize=8, color=INK_MUTED, va="bottom")
    ax_a.grid(axis="y")
    open_frame(ax_a)
    plt.setp(ax_a.get_xticklabels(), visible=False)
    _title(ax_a, 0, "Shorelines, the year before and the year after")
    structures(ax_a, label=True)

    # (b) the change at each transect
    ch = t.change_m.to_numpy()
    w = np.median(np.diff(x)) * 0.9
    ax_b.axvspan(first - 0.5, last + 0.5, color="0.92", lw=0, zorder=0)
    ax_b.bar(x, ch, width=w, color=np.where(ch >= 0, C_1997, C_1984), lw=0, zorder=2)
    ax_b.axhline(0, color=INK, lw=0.6, zorder=3)
    ax_b.axhline(bg, color=INK_MUTED, lw=0.9, ls=(0, (3, 2)), zorder=3)
    ax_b.set_ylabel("After − before (m)")
    ax_b.set_xlabel(DOMAIN_AXIS_LABEL)
    ax_b.set_xlim(x.min() - 0.3, x.max() + 0.3)
    ax_b.grid(axis="y")
    lo_, hi_ = ax_b.get_ylim()
    ax_b.set_ylim(lo_, hi_ + 0.18 * (hi_ - lo_))
    open_frame(ax_b)
    _title(ax_b, 1, "Change at each transect")
    of, ol = summary_row.observed_first_gis, summary_row.observed_last_gis
    if pd.notna(of):
        ax_b.annotate("", xy=(of - 0.5, 0.94), xytext=(ol + 0.5, 0.94), xycoords=("data", "axes fraction"),
                      arrowprops=dict(arrowstyle="|-|", color=INK, lw=1.0, mutation_scale=3), annotation_clip=False)
        ax_b.text((of + ol) / 2, 0.96, f"CoastSat core, GIS {int(of)}–{int(ol)}", transform=ax_b.get_xaxis_transform(),
                  ha="center", va="bottom", fontsize=7.5, color=INK)

    # (c) true scale on imagery, centred on the largest smoothed gain inside the footprint
    import contextily as cx
    from rasterio.warp import transform_bounds
    cx.set_cache_dir(str(TILE_CACHE))
    sm = t.change_m.rolling(X.RUN_TRANSECTS, center=True, min_periods=3).median()
    sm = sm.where(t.domain_number.between(first, last))
    c = t.loc[sm.idxmax()]
    ce, cn = (c.before_e + c.after_e) / 2, (c.before_n + c.after_n) / 2
    bx = (ce - ZOOM_HALF_M, cn - ZOOM_HALF_M * 2.2, ce + ZOOM_HALF_M, cn + ZOOM_HALF_M * 2.2)
    W, S, E, N = transform_bounds(MAP_CRS, "EPSG:3857", *bx)
    img, ext = cx.bounds2img(W, S, E, N, zoom=TILE_ZOOM, source=cx.providers.Esri.WorldImagery, ll=False)
    img, ext = cx.warp_tiles(img, ext, t_crs=MAP_CRS)
    ax_c.imshow(img, extent=ext, zorder=0)
    near = t[(t.before_e.between(bx[0] - 300, bx[2] + 300)) & (t.before_n.between(bx[1] - 300, bx[3] + 300))]
    ax_c.plot(near.before_e, near.before_n, color=C_1984, lw=1.6, zorder=3)
    ax_c.plot(near.after_e, near.after_n, color=C_1997, lw=1.6, zorder=4)
    ax_c.set_xlim(bx[0], bx[2])
    ax_c.set_ylim(bx[1], bx[3])
    ax_c.set_aspect("equal")
    ax_c.set_xticks([])
    ax_c.set_yticks([])
    sx, sy = bx[0] + 60, bx[1] + 80
    ax_c.plot([sx, sx + 200], [sy, sy], color="white", lw=3, solid_capstyle="butt", zorder=6)
    ax_c.text(sx + 100, sy + 35, "200 m", color="white", ha="center", va="bottom", fontsize=8, zorder=6,
              path_effects=X_DARK)
    ax_c.text(0.04, 0.97, f"(c) GIS {int(c.domain_number)}, true scale", transform=ax_c.transAxes,
              color="white", fontsize=8.5, va="top", fontweight="bold", path_effects=X_DARK)
    ax_c.text(0.04, 0.925, f"{c.change_m:+.0f} m here", transform=ax_c.transAxes,
              color="white", fontsize=8, va="top", path_effects=X_DARK)

    fig.suptitle(f"{p.name} {p.year}: model GIS {first}–{last}, placed {d.start_date} to {d.end_date}",
                 x=0.08, ha="left", y=0.985, fontsize=10.5, color=INK)
    handles = [Line2D([], [], color=C_1984, lw=1.4, label=f"Before: median {win_b}"),
               Line2D([], [], color=C_1997, lw=1.4, label=f"After: median {win_a}"),
               Patch(color=C_1997_FILL, label="Seaward gain"), Patch(color=C_1984_FILL, label="Landward loss"),
               Patch(color="0.92", label="Model footprint"),
               Line2D([], [], color=INK_MUTED, lw=0.9, ls=(0, (3, 2)), label=f"Background change, {bg:+.1f} m")]
    fig.legend(handles=handles, loc="lower center", ncol=3, frameon=False, bbox_to_anchor=(0.5, 0.0),
               fontsize=8)
    stem = f"nourishment_shorelines_{p.year}_{p.name.split()[0].lower()}"
    out = save(fig, OUT_DIR / stem, close=True)[0]
    record_caption(out, (
        f"{p.name}, {p.year}, in CoastSat. (a) The median shoreline over {win_b} (red, before placement) and over "
        f"{win_a} (blue, after), each transect's position measured from a degree-{BASELINE_DEG} baseline fitted "
        "through the before line, so the coast is straightened and the cross-shore axis is stretched; blue shading "
        "is seaward gain, red landward loss. (b) After minus before at each transect; the dashed line is the median "
        "change of the transects outside every fill footprint, and the bracket the half-peak core from "
        "nourishment_extent_coastsat.py. (c) Both lines on Esri World Imagery at true scale, centred on the largest "
        f"5-transect running-median gain inside the footprint (GIS {int(c.domain_number)}); the imagery is current, "
        "not of either window. Grey columns are the model footprint (HATTERAS_NOURISHMENT_PROJECTS)."))
    return out


X_DARK = None


def main():
    global X_DARK
    import matplotlib.patheffects as pe
    X_DARK = [pe.withStroke(linewidth=2.0, foreground="black")]
    apply_style()
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    frame = X.transect_frame()
    dates = pd.read_csv(X.DATES_CSV).set_index("project")
    summary = pd.read_csv(X.OUT_DIR / "nourishment_extent_coastsat_summary.csv").set_index("project")
    geom = transect_geometry(set(frame.index))
    for p in sorted(HATTERAS_NOURISHMENT_PROJECTS, key=lambda q: (q.year, min(q.gis_domains))):
        t, windows = X.project_change(p, dates.loc[p.name], frame)
        t = add_map_coords(t, geom).sort_values("x")
        print(draw(p, dates.loc[p.name], t, summary.loc[p.name], windows))


if __name__ == "__main__":
    main()
