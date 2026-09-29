"""
offset_build_figures.py
==============================================================================
How the BRIE shoreline offset is built from a digitised dune line, step by
step, for the builds the current runs read (1996 and 2010 starts).

    python scripts/figure_making/pipeline/2-brie-offset/offset_build_figures.py

Writes output/figures/pipeline/2-brie-offset/offset_build_<year>.png.

THE STEPS DRAWN (the producers, not re-implemented)
    1. duneline_to_raw_offsets.py intersects the dune line with the 100 m
       transects. Each transect starts on the offshore datum line and runs
       west across the island; the STATION of a crossing is its distance
       along the transect from the datum. The raw file the build kept is read,
       not recomputed.
    2. island_offset_hybrid.py averages the ~5 transects of each 500 m domain
       and subtracts the smallest domain mean, so the offset is 0 at the most
       seaward domain and positive landward.
    3. cascade_pipeline.hindcast.pad_offset_ring closes BRIE's periodic line
       with 15 buffer domains per side (a cubic Hermite from GIS 90 back
       round to GIS 1). The padded file IS what the runner hands Cascade.

Every path resolves through site_layer.hat_topo_version: the dune-line vintage
from DUNE_LINE_FOR_YEAR, the build from offset_version (env > CURRENT > the
only v<n>). The ocean is on the right in the plan panels (easting across).
"""

from __future__ import annotations

import sys
from pathlib import Path

import matplotlib
import matplotlib.ticker

matplotlib.use("Agg")
import geopandas as gpd  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402
from matplotlib.patches import Rectangle  # noqa: E402

REPO = next(p for p in Path(__file__).resolve().parents if (p / "pyproject.toml").exists())
sys.path.insert(0, str(REPO / "scripts"))

from site_layer import hat_topo_version as tv  # noqa: E402
from site_layer.hat_observed_rates import DOMAIN_BOXES  # noqa: E402
from site_layer.hat_figure_style import (  # noqa: E402
    apply_style, C, C_1997, INK, INK_MUTED, DOMAIN_AXIS_LABEL, figsize, figure_dir, save,
    record_caption, _title, open_frame, town_bands,
)
from cascade_pipeline.hindcast import pad_offset_ring  # noqa: E402

OUT = figure_dir("pipeline", "2-brie-offset")
YEARS = (1996, 2010)
ZOOM_GIS = (44, 47)
BUFFERS = 15
DATUM_X = 460198.45        # the offshore datum line the transects start on (EPSG:3725)
UNSTABLE_DEG = 42.0


def load_transects():
    t = gpd.read_file(tv.TRANSECT_FILE_100M)
    t.columns = [c.split(".")[-1] for c in t.columns]
    t = t.loc[:, ~t.columns.duplicated()]
    t["LineID"] = t["LineID"].astype(int)
    return t


def load_build(year):
    vintage = tv.dune_line_for_year(year)
    version = tv.offset_version(year)
    build = tv.offset_build_dir(year, version)
    raw = pd.read_csv(build / f"{vintage}_duneline_offset_raw.csv")
    unpadded = pd.read_csv(tv.offset_file(year, "unpadded", version=version))
    padded = np.loadtxt(tv.offset_file(year, "padded", version=version), skiprows=1)
    line = gpd.read_file(tv.DUNELINE_DIR / f"duneline_{vintage}.geojson")
    return vintage, version, build, raw, unpadded, padded, line


def fig_offset_build(year, transects, boxes):
    vintage, version, build, raw, unpadded, padded, line = load_build(year)
    line = line.to_crs(transects.crs)
    raw = raw[raw.domain_id.between(1, 90)]
    real = unpadded.iloc[:, 1].to_numpy()
    means = raw.groupby("domain_id").ORIG_LEN.mean()
    baseline = means.min()
    closed = pad_offset_ring(real, BUFFERS)
    if np.abs(closed - padded).max() > 1e-6:
        raise RuntimeError(f"{year}: padded file is not pad_offset_ring(unpadded) - build out of step")

    fig = plt.figure(figsize=figsize("double", height=8.0), constrained_layout=True)
    gs = fig.add_gridspec(2, 3, width_ratios=[0.5, 0.55, 1.9], height_ratios=[1.35, 1])

    # (a) the reach in plan
    ax_a = fig.add_subplot(gs[0, 0])
    km = 1000.0
    for g in transects.geometry[::3]:
        x, y = g.xy
        ax_a.plot(np.array(x) / km, np.array(y) / km, color="0.8", lw=0.3)
    lx, ly = line.geometry.iloc[0].xy
    ax_a.plot(np.array(lx) / km, np.array(ly) / km, color=C_1997, lw=1.0)
    ax_a.axvline(DATUM_X / km, color=INK, lw=0.8, ls=(0, (3, 2)))
    zb = boxes[boxes.domain_id.between(*ZOOM_GIS)].total_bounds
    ax_a.add_patch(Rectangle((zb[0] / km, zb[1] / km), (zb[2] - zb[0]) / km, (zb[3] - zb[1]) / km,
                             fill=False, ec=C["ACCENT"], lw=1.0))
    ax_a.set_aspect("equal")
    ax_a.set_xlabel("easting (km)")
    ax_a.set_ylabel("northing (km, UTM 18N)")
    ax_a.tick_params(labelsize=7)
    open_frame(ax_a)
    _title(ax_a, 0, "")

    # (b) the zoom: transects, crossings, domain boxes; the datum is off to the right
    ax_b = fig.add_subplot(gs[0, 1])
    sub = raw[raw.domain_id.between(*ZOOM_GIS)]
    xc = sub.x.median()
    x_lo, x_hi = xc - 320, xc + 320
    ax_b.set_xlim(x_lo, x_hi)
    zb2 = boxes[boxes.domain_id.between(*ZOOM_GIS)]
    ax_b.set_ylim(zb2.total_bounds[1] - 20, zb2.total_bounds[3] + 20)
    cols = [C_1997, C["ACCENT"]]
    for k, (_, bx) in enumerate(zb2.sort_values("domain_id").iterrows()):
        b = bx.geometry.bounds
        ax_b.axhspan(b[1], b[3], color="0.94" if k % 2 else "white", lw=0, zorder=0)
        ax_b.text(x_lo + 12, (b[1] + b[3]) / 2, f"GIS {int(bx.domain_id)}", ha="left", va="center",
                  fontsize=7.5, color=INK_MUTED)
    ids = set(sub.LineID)
    for _, tr in transects[transects.LineID.isin(ids)].iterrows():
        x, y = tr.geometry.xy
        ax_b.plot(x, y, color="0.6", lw=0.6, zorder=1)
    ax_b.plot(lx, ly, color=C_1997, lw=1.3, zorder=2)
    for k, (gid, g) in enumerate(sub.groupby("domain_id")):
        ax_b.plot(g.x, g.y, "o", color=cols[k % 2], ms=4, zorder=3)
    one = sub.iloc[len(sub) // 2]
    ax_b.annotate("", xy=(one.x, one.y + 14), xytext=(x_hi, one.y + 14),
                  arrowprops=dict(arrowstyle="<-", lw=0.7, color=INK))
    ax_b.text(x_hi - 8, one.y + 24, f"{one.ORIG_LEN / 1000:.2f} km\nto datum", ha="right", va="bottom",
              fontsize=7)
    ax_b.set_aspect("equal")
    ax_b.set_xlabel("easting (m)  ·  ocean →")
    ax_b.tick_params(labelleft=False, labelsize=7)
    ax_b.xaxis.set_major_locator(matplotlib.ticker.MaxNLocator(3))
    open_frame(ax_b)
    ax_b.legend(handles=[Line2D([], [], color=C_1997, lw=1.3, label=f"dune line {vintage}"),
                         Line2D([], [], color="0.6", lw=0.6, label="100 m transects"),
                         Line2D([], [], color=C["ACCENT"], marker="o", ls="", ms=4, label="crossing")],
                loc="upper center", bbox_to_anchor=(0.5, -0.1), frameon=False, fontsize=7)
    _title(ax_b, 1, "")

    # (c) every transect's station, and the domain means
    ax_c = fig.add_subplot(gs[0, 2])
    ax_c.plot(raw.domain_id + (raw.groupby("domain_id").cumcount() - 2) * 0.12, raw.ORIG_LEN, ".",
              color="0.65", ms=2.5, label="transect station")
    ax_c.step(means.index, means.values, where="mid", color=INK, lw=1.2, label="domain mean")
    ax_c.axhline(baseline, color=C["ACCENT"], lw=0.8, ls=(0, (3, 2)), label="smallest mean (the zero)")
    ax_c.invert_yaxis()
    ax_c.set_ylabel("distance from the offshore datum (m)")
    ax_c.set_xlabel(DOMAIN_AXIS_LABEL)
    ax_c.set_xlim(0.5, 90.5)
    town_bands(ax_c)
    open_frame(ax_c)
    ax_c.legend(frameon=False, fontsize=7, loc="center right")
    _title(ax_c, 2, "Transect stations, domain means")

    # (d) zeroed and closed: the model input
    ax_d = fig.add_subplot(gs[1, :])
    gis_pad = np.arange(1 - BUFFERS, 91 + BUFFERS)
    ax_d.axvspan(gis_pad[0] - 0.5, 0.5, color="0.93", lw=0)
    ax_d.axvspan(90.5, gis_pad[-1] + 0.5, color="0.93", lw=0)
    ax_d.plot(gis_pad[:BUFFERS + 1], padded[:BUFFERS + 1], color=C["ACCENT"], lw=1.2, ls=(0, (3, 1.5)))
    ax_d.plot(gis_pad[-BUFFERS - 1:], padded[-BUFFERS - 1:], color=C["ACCENT"], lw=1.2, ls=(0, (3, 1.5)),
              label="buffer closure (cubic Hermite, periodic)")
    ax_d.plot(np.arange(1, 91), real, color=INK, lw=1.4, label="GIS 1-90: domain mean − smallest mean")
    ang = np.degrees(np.arctan2(np.diff(np.r_[padded, padded[0]]), 500.0))
    ax_d.set_ylabel("shoreline offset (m, landward +)")
    ax_d.set_xlabel("padded domain, numbered as GIS (buffers shaded)")
    ax_d.set_xlim(gis_pad[0] - 0.5, gis_pad[-1] + 0.5)
    open_frame(ax_d)
    ax_d.legend(frameon=False, fontsize=7.5, loc="upper right")
    _title(ax_d, 3, f"The model input: {len(padded)} values, one per 500 m domain")

    out = save(fig, OUT / f"offset_build_{year}.png")
    plt.close(fig)
    rel = build.relative_to(REPO).as_posix()
    record_caption(out[0],
        f"How the BRIE shoreline offset for the {year} start is built from the {vintage} dune line "
        f"(build {rel}). (a) The reach in plan, north up, ocean at the right: the digitised dune line (blue), "
        "every third 100 m transect (grey) and the offshore datum line they start on (dashed); the purple box "
        f"is panel b. (b) GIS {ZOOM_GIS[0]}-{ZOOM_GIS[1]}: each transect's crossing with the dune line (points, "
        "alternating colour by domain; shaded bands are the 500 m domain boxes). The station of a crossing is "
        "its distance along the transect from the datum, measured exactly in shapely "
        "(duneline_to_raw_offsets.py; the landward-most crossing where a transect meets the line twice); the "
        "arrow gives one transect's station, the datum lying off the panel to the right. "
        "(c) Every transect's station (grey) and the mean of the ~5 transects in each domain (black), drawn "
        "with the datum at the top so landward is down. The smallest domain mean (dashed, "
        f"{baseline:.0f} m) becomes zero. (d) The model input: each domain mean minus that zero, so the offset "
        f"runs from 0 to {real.max():.0f} m landward, and {BUFFERS} buffer domains each side (shaded) that "
        "close BRIE's periodic shoreline with a cubic Hermite matched to the island's end slopes "
        "(cascade_pipeline.hindcast.pad_offset_ring). The padded file equals that closure to 1e-6 m, so this "
        f"curve is exactly what Cascade receives, in metres. The steepest angle between neighbouring padded "
        f"domains is {np.abs(ang).max():.1f}° against BRIE's ~{UNSTABLE_DEG:.0f}° anti-diffusive limit. "
        "GIS 1 is Cape Point, GIS 90 Pea Island.")
    return out


def main():
    apply_style()
    import geopandas as _g
    boxes = _g.read_file(DOMAIN_BOXES).to_crs("EPSG:3725")
    transects = load_transects()
    for year in YEARS:
        out = fig_offset_build(year, transects, boxes)
        print(out[0].relative_to(REPO))


if __name__ == "__main__":
    main()
