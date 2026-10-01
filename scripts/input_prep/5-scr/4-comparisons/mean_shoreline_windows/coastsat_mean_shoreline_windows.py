"""
One period's CoastSat mean shoreline over two windows, differenced transect by transect and domain by domain.

    python scripts/input_prep/5-scr/4-comparisons/mean_shoreline_windows/coastsat_mean_shoreline_windows.py
    python scripts/input_prep/5-scr/4-comparisons/mean_shoreline_windows/coastsat_mean_shoreline_windows.py --periods 1996

The 3-yr calendar window against the DEM-centred +/-1 yr window; figures,
an island map, a datum check and provenance per period. Details: scripts/input_prep/5-scr/4-comparisons/README.md.

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-09-29
"""

from __future__ import annotations

import argparse
import sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd

_REPO = next(_p for _p in Path(__file__).resolve().parents
             if (_p / "pyproject.toml").exists())
sys.path.insert(0, str(_REPO / "scripts"))
sys.path.insert(0, str(_REPO / "scripts" / "input_prep" / "5-scr" / "lib"))
import scr_paths  # noqa: E402,F401

import coastsat_mean_shoreline as cms  # noqa: E402
from site_layer import hat_figure_style as fs  # noqa: E402
from site_layer.hat_observed_rates import (  # noqa: E402
    MEAN_SHORELINE_WINDOWS, mean_shoreline_csv,
)

# --- CONFIG ------------------------------------------------------------------
# period start -> (the 3-yr calendar window, the 2-yr DEM-centred window)
PERIODS = {
    1996: (cms.Window.from_years(1995, 1997), cms.Window.centred_on("alace_1996")),
    2010: (cms.Window.from_years(2009, 2011), cms.Window.centred_on("usace_2009")),
}
# The unmodified input grey, the modification under test purple (style rule).
C_CAL, C_DEM = fs.C["BASE"], fs.C["ACCENT"]
# -----------------------------------------------------------------------------


# The stored per-transect means of one window, included transects only
def load(window):
    df = pd.read_csv(mean_shoreline_csv(*window.key))
    return df[df["included"]].set_index("transect_id")


# Per transect and per domain
def compare(cal, dem):
    a, b = load(cal), load(dem)
    both = a.index.intersection(b.index)
    t = pd.DataFrame({
        "domain_number": a.loc[both, "domain_number"].astype(int),
        "mean_3yr_m": a.loc[both, "mean_chainage_m"],
        "mean_2yr_m": b.loc[both, "mean_chainage_m"],
        "n_obs_3yr": a.loc[both, "n_obs"].astype(int),
        "n_obs_2yr": b.loc[both, "n_obs"].astype(int),
        "se_3yr_m": a.loc[both, "se_chainage_m"],
        "se_2yr_m": b.loc[both, "se_chainage_m"],
    })
    t["diff_m"] = t["mean_2yr_m"] - t["mean_3yr_m"]
    only = {"only_3yr": sorted(set(a.index) - set(both)),
            "only_2yr": sorted(set(b.index) - set(both))}
    g = t.groupby("domain_number")
    d = pd.DataFrame({
        "n_transects": g.size(),
        "diff_m": g["diff_m"].mean(),
        "diff_min_m": g["diff_m"].min(),
        "diff_max_m": g["diff_m"].max(),
        "se_3yr_m": g["se_3yr_m"].median(),
        "se_2yr_m": g["se_2yr_m"].median(),
        "median_n_obs_3yr": g["n_obs_3yr"].median(),
        "median_n_obs_2yr": g["n_obs_2yr"].median(),
    })
    d["diff_minus_island_mean_m"] = d["diff_m"] - d["diff_m"].mean()
    return t.sort_index(), d, only


# The two windows' profiles and their difference
def figure(t, d, cal, dem, period, folder):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D

    fs.apply_style()
    fig, axes = plt.subplots(2, 1, figsize=fs.figsize("double", aspect=0.62),
                             sharex=True, layout="constrained",
                             gridspec_kw={"height_ratios": [1.6, 1]})
    ax = axes[0]
    ax.axhspan(-fs.CELL_M, fs.CELL_M, color="0.95", lw=0, zorder=0)
    ax.axhline(0, color=fs.INK, lw=0.6, zorder=1)
    ax.scatter(t["domain_number"], t["diff_m"], s=5, color=fs.C["ACCENT_FILL"],
               lw=0, zorder=2)
    ax.plot(d.index, d["diff_m"], "-o", color=C_DEM, ms=2.4, lw=1.1, zorder=3)
    ax.set_xlim(0.5, 90.5)
    lim = max(fs.CELL_M * 1.5, float(np.nanmax(np.abs(t["diff_m"]))) * 1.08)
    # headroom above the data for the village strip, below it for the legend
    ax.set_ylim(-lim * 1.3, lim * 1.3)
    ax.set_ylabel("2-yr − 3-yr position (m)\n+ seaward")
    fs.town_bands(ax, strip=0.07)
    fs._title(ax, 0, "Mean shoreline, DEM-centred 2-yr minus calendar 3-yr")
    ax.legend([Line2D([], [], color=C_DEM, marker="o", ms=2.4, lw=1.1),
               Line2D([], [], ls="none", marker="o", ms=2.4, mfc=fs.C["ACCENT_FILL"], mec="none"),
               matplotlib.patches.Patch(color="0.95")],
              ["domain mean", "transect", "±1 model cell (10 m)"],
              loc="lower left", ncol=3, frameon=False, fontsize=7.5)

    ax = axes[1]
    ax.plot(d.index, d["se_3yr_m"], "-", color=C_CAL, lw=1.1,
            label="3-yr, {0}".format(cal.span))
    ax.plot(d.index, d["se_2yr_m"], "-", color=C_DEM, lw=1.1,
            label="2-yr, {0}".format(dem.span))
    # headroom so the legend sits clear of the lines
    ax.set_ylim(0, float(np.nanmax(d[["se_3yr_m", "se_2yr_m"]].to_numpy())) * 1.35)
    ax.set_ylabel("standard error (m)")
    ax.set_xlabel(fs.DOMAIN_AXIS_LABEL)
    ax.legend(loc="upper left", ncol=2, frameon=False, fontsize=7.5)
    fs._title(ax, 1, "Sampling error of each window's mean")

    a = SURVEYS[dem.anchor]
    fs.caption(fig, (
        "The {p} period's CoastSat mean shoreline over two windows: the calendar "
        "{c} (the window the shoreline island offset v1 was built from) and "
        "{lo} to {hi}, ±1 yr of the {s} flights (flown {f0} to {f1}), the start "
        "DEM's lidar. (a) The 2-yr mean minus the 3-yr mean along each CoastSat "
        "transect, positive where the DEM-centred line is seaward: transects as "
        "light dots, the domain mean as the line; the grey band is ±1 Barrier3D "
        "cell (10 m). Island mean {m:+.1f} m; after removing it, the domain means "
        "vary by sd {sd:.1f} m, the part an offset zeroed on its own minimum can "
        "see. (b) The median standard error of the transect means in each domain, "
        "for each window. The windows share satellite passes, so the two means "
        "are not independent and no significance is attached to (a). Village "
        "spans along the top of (a).".format(
            p=period, c=cal.span, lo=dem.lo_iso, hi=dem.hi_iso, s=a["survey"],
            f0=a["flown"][0], f1=a["flown"][1], m=d["diff_m"].mean(),
            sd=d["diff_minus_island_mean_m"].std())))
    out = fs.save(fig, folder / "mean_shoreline_windows_{0}.png".format(period), close=True)
    print("  figure -> {0}".format(out[0]))


SURVEYS = cms.SURVEY_ANCHORS


# The provenance beside a period's comparison
def write_provenance(t, d, only, cal, dem, period, folder, built_on):
    a = SURVEYS[dem.anchor]
    big = d.reindex(d["diff_minus_island_mean_m"].abs().sort_values(ascending=False).index).head(5)
    big_rows = "\n".join("| GIS {0} | {1:+.1f} | {2:+.1f} | {3:.1f} / {4:.1f} |".format(
        int(i), r["diff_m"], r["diff_minus_island_mean_m"], r["se_3yr_m"], r["se_2yr_m"])
        for i, r in big.iterrows())
    over_cell = int((d["diff_m"].abs() > fs.CELL_M).sum())
    only_txt = ("; ".join("{0}: {1}".format(k.replace("_", " "), ", ".join(v))
                          for k, v in only.items() if v) or "none")
    text = """# mean_shoreline_windows/{p} -- the {p} mean shoreline, 3-yr calendar vs 2-yr DEM-centred

Written by `scripts/input_prep/5-scr/4-comparisons/mean_shoreline_windows/coastsat_mean_shoreline_windows.py`
on {built_on}.

## The two windows

| | window | folder | positions per transect (median) |
|---|---|---|---|
| 3-yr | calendar {cal_span} | `1-observations/mean_shoreline/{cal_label}/` | {n3:.0f} |
| 2-yr | {dem_lo} to {dem_hi}, ±1 yr of the {survey} (flown {f0} to {f1}; {src}) | `1-observations/mean_shoreline/{dem_label}/` | {n2:.0f} |

The 3-yr line built the shoreline island offset v1; the 2-yr line builds v2
(Hannah, 2026-09-29). Both are read as stored; nothing is re-averaged here.

## The difference (2-yr minus 3-yr, + seaward)

| | |
|---|---|
| transects in both windows | {n_both} (not in both: {only_txt}) |
| per transect | mean {tm:+.2f} m, sd {tsd:.2f} m, p5 {tp5:+.1f} m, p95 {tp95:+.1f} m, largest {tmax:.1f} m |
| per domain | mean {dm:+.2f} m, range {dmin:+.1f} to {dmax:+.1f} m |
| per domain, island mean removed | sd {dsd:.2f} m |
| domains beyond one 10 m cell | {over_cell} of {n_dom} |
| median standard error, 3-yr / 2-yr | {se3:.2f} / {se2:.2f} m |

The five domains that move most once the island mean is removed:

| domain | difference (m) | minus island mean (m) | SE 3-yr / 2-yr (m) |
|---|---|---|---|
{big_rows}

## How to read it

- **Only the alongshore-varying part reaches the model.** The island offset
  is zeroed on its own minimum, so the island-mean shift ({dm:+.1f} m) cancels;
  the {dsd:.1f} m sd after removing it is what can change BRIE's shoreline shape.
- **No significance is attached.** The windows overlap and share passes, so
  their means are not independent samples; the standard errors are shown so a
  difference can be read against the sampling noise of either mean.
- **Positive is seaward.** CoastSat chainage grows offshore along every
  transect here, so a positive difference puts the DEM-centred line seaward
  of the calendar one.

## Files

| file | what it is |
|---|---|
| `mean_shoreline_windows_{p}.png` | (a) the difference per transect and domain, (b) the standard error of each window's mean |
| `mean_shoreline_windows_{p}_datum.png` | both lines as the offset build sees them: distance from the shared offshore datum along the 100 m model transects, per domain, six sections on a 2 x 3 grid (the layout of 2-brie-offset's duneline_vs_shoreline figure), each with a strip of the difference beside it |
| `mean_shoreline_windows_{p}_island.png` | the whole island in six north-up segments of 15 domains, the 3-yr line coloured by the difference (blue seaward, red landward) |
| `supporting/datum_stations_{p}.csv` | per domain: each window's station from the offshore datum and the seaward shift of the 2-yr line |
| `mean_shoreline_windows_{p}_lines_largest.png` | both lines on the photograph nearest the DEM survey, the six domains that move most once the island mean is removed (middle 250 m of each) |
| `mean_shoreline_windows_{p}_lines_sites.png` | the same at the centre domains of the six mean_shoreline imagery sites |
| `supporting/domain_comparison_{p}.csv` | per domain: n transects, mean / min / max difference, difference minus the island mean, median SE and positions for each window |
| `supporting/transect_comparison_{p}.csv` | per transect: both means, n, SE, and the difference |
| `supporting/CAPTIONS.md`, `supporting/*.pdf` | the caption and the vector copy |
""".format(
        p=period, built_on=built_on, cal_span=cal.span, cal_label=cal.label,
        dem_lo=dem.lo_iso, dem_hi=dem.hi_iso, dem_label=dem.label,
        survey=a["survey"], f0=a["flown"][0], f1=a["flown"][1], src=a["source"],
        n3=t["n_obs_3yr"].median(), n2=t["n_obs_2yr"].median(),
        n_both=len(t), only_txt=only_txt,
        tm=t["diff_m"].mean(), tsd=t["diff_m"].std(), tp5=t["diff_m"].quantile(0.05),
        tp95=t["diff_m"].quantile(0.95), tmax=t["diff_m"].abs().max(),
        dm=d["diff_m"].mean(), dmin=d["diff_m"].min(), dmax=d["diff_m"].max(),
        dsd=d["diff_minus_island_mean_m"].std(), over_cell=over_cell, n_dom=len(d),
        se3=t["se_3yr_m"].median(), se2=t["se_2yr_m"].median(), big_rows=big_rows)
    (folder / "PROVENANCE.md").write_text(text, encoding="utf-8")


# The two lines themselves, on a photograph (Hannah, 2026-09-29

# The photograph nearest the start DEM's survey
PHOTO_YEAR = {1996: 1996, 2010: 2008}
N_LARGEST = 6
LINE_LAND_M, LINE_SEA_M = 45.0, 35.0       # panel reach either side of the lines
# The middle half of each domain
LINE_ALONG_M = 250.0
# On the photographs the two windows are the house RdBu pair, 3-yr red and 2-yr blue (Hannah, 2026-09-29
C_CAL_PHOTO, C_DEM_PHOTO = fs.C_1984, fs.C_1997
HALO_LW = 1.5


# The window's mean line in alongshore order, as coastsat_mean_shoreline strung it
def line_points(window):
    return cms.line_vertices(pd.read_csv(mean_shoreline_csv(*window.key)))


# One panel per domain
def lines_figure(d, cal, dem, period, folder, gis, tag, what):
    import geopandas as gpd
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.patheffects as pe
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    import coastsat_mean_shoreline_on_imagery as moi
    from site_layer.hat_observed_rates import DOMAIN_BOXES

    fs.apply_style()
    boxes = gpd.read_file(DOMAIN_BOXES).to_crs(moi.CRS)
    boxes["gis"] = np.arange(1, len(boxes) + 1)          # south -> north
    R = moi._review_module()
    year = PHOTO_YEAR[period]
    im = R.Imagery(year, moi.CRS)
    if year in moi.PHOTO_SOURCES:
        im.date = moi.PHOTO_SOURCES[year]["date"]
    a3, a2 = line_points(cal), line_points(dem)

    ext = []
    for g in gis:
        b0, b1 = boxes.loc[boxes["gis"] == g].total_bounds[[1, 3]]
        y0, y1 = (b0 + b1 - LINE_ALONG_M) / 2, (b0 + b1 + LINE_ALONG_M) / 2
        xs = pd.concat([a3.loc[a3["y"].between(y0, y1), "x"], a2.loc[a2["y"].between(y0, y1), "x"]])
        ext.append((xs.min() - LINE_LAND_M, y0, xs.max() + LINE_SEA_M, y1))
    widths = [b[2] - b[0] for b in ext]
    # the panel height the page width allows, so the figure has no dead space
    scale = (fs.FIG_W_DOUBLE - 0.5) / sum(widths)           # in per m
    panel_h = min(LINE_ALONG_M * scale, 6.5)
    fig = plt.figure(figsize=(fs.FIG_W_DOUBLE, panel_h + 1.15), layout="constrained")
    axes = fig.subplots(1, len(gis), gridspec_kw={"width_ratios": widths})
    for i, (ax, g, b) in enumerate(zip(axes, gis, ext)):
        x0, y0, x1, y1 = b
        img = im.read(x0, x1, y0, y1, 0.4, moi.CRS).copy()
        img[img.max(axis=2) == 0] = 255
        ax.imshow(img, extent=(x0, x1, y0, y1), origin="upper", zorder=0,
                  interpolation="bilinear")
        for pts, c in ((a3, C_CAL_PHOTO), (a2, C_DEM_PHOTO)):
            s = pts[pts["y"].between(y0 - 60, y1 + 60)]
            ax.plot(s["x"], s["y"], "-o", color=c, lw=0.9, ms=2.6, mec="white",
                    mew=0.4, zorder=3,
                    path_effects=[pe.withStroke(linewidth=HALO_LW, foreground="white")])
        ax.set_xlim(x0, x1)
        ax.set_ylim(y0, y1)
        ax.set_aspect("equal")
        ax.set_xticks([])
        ax.set_yticks([])
        fs.spines_for_image(ax)
        ax.set_title("({0})  GIS {1}\n{2:+.1f} m".format(
            chr(97 + i), g, d.loc[g, "diff_m"]), loc="left", fontsize=8.5, pad=3)
    fs._scalebar(axes[0], 20.0, show_cells=False)
    fs._north_arrow(axes[-1], x=0.75, y=0.06)
    fig.legend([Line2D([], [], color=C_CAL_PHOTO, lw=0.9, marker="o", ms=2.6, mec="white", mew=0.4),
                Line2D([], [], color=C_DEM_PHOTO, lw=0.9, marker="o", ms=2.6, mec="white", mew=0.4)],
               ["3-yr, calendar {0}".format(cal.span), "2-yr, {0}".format(dem.span)],
               loc="outside lower center", ncol=2, frameon=False)
    fs.caption(fig, (
        "The {p} period's CoastSat mean shoreline over two windows, drawn on "
        "{photo}: red, the calendar {c} mean (the shoreline island offset "
        "v1's line); blue, the {lo} to {hi} mean, ±1 yr of the {s} (v2's "
        "line). {what} Each panel is the middle 250 m of one Barrier3D domain, "
        "north up, all at one scale; the panel width follows the lines. Dots "
        "are the per-transect means (~50 m apart) that the lines join. Under "
        "each letter, the whole domain's mean difference, 2-yr minus 3-yr, "
        "positive seaward. Photograph frames of different exposure meet at "
        "straight seams. White is outside the photographs.".format(
            p=period, photo=moi.photo_ref([im], dem), c=cal.span, lo=dem.lo_iso,
            hi=dem.hi_iso, s=SURVEYS[dem.anchor]["survey"], what=what)))
    out = fs.save(fig, folder / "mean_shoreline_windows_{0}_lines_{1}.png".format(period, tag),
                  vector=False, close=True)
    print("  figure -> {0}".format(out[0].name))


# The whole island in six north-up segments of 15 domains (Hannah, 2026-09-29
SIX_SEGMENTS = [(1, 15), (16, 30), (31, 45), (46, 60), (61, 75), (76, 90)]
DIFF_CMAP = "RdBu"            # blue seaward, red landward; the house RdBu poles


# The difference along the island on a map
def island_figure(t, d, cal, dem, period, folder):
    import geopandas as gpd
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.collections import LineCollection
    from shapely.geometry import box
    from site_layer.hat_map_layers import ISLAND_OUTLINE
    from site_layer.hat_observed_rates import DOMAIN_BOXES

    fs.apply_style()
    crs = cms.TARGET_CRS
    boxes = gpd.read_file(DOMAIN_BOXES).to_crs(crs)
    boxes["gis"] = np.arange(1, len(boxes) + 1)          # south -> north
    island = gpd.read_file(ISLAND_OUTLINE).to_crs(crs)
    a3 = line_points(cal).set_index("transect_id")
    pts = a3.join(t[["diff_m"]], how="inner").sort_values(["site", "transect_number"])
    xy = pts[["x", "y"]].to_numpy()
    seg_xy = np.stack([xy[:-1], xy[1:]], axis=1)
    seg_c = (pts["diff_m"].to_numpy()[:-1] + pts["diff_m"].to_numpy()[1:]) / 2
    vmax = float(np.ceil(np.nanmax(np.abs(t["diff_m"])) / 2.5) * 2.5)
    norm = matplotlib.colors.Normalize(-vmax, vmax)

    # ONE extent for every panel, centred on each segment, so all six share a scale
    ext = [boxes[boxes["gis"].between(lo, hi)].total_bounds for lo, hi in SIX_SEGMENTS]
    w = max(e[2] - e[0] for e in ext) + 900
    h = max(e[3] - e[1] for e in ext) + 300
    panel_w = (fs.FIG_W_DOUBLE - 0.2) / 6
    fig, axes = plt.subplots(1, 6, figsize=(fs.FIG_W_DOUBLE, panel_w * h / w + 0.95),
                             layout="constrained")
    for i, (ax, (lo, hi), e) in enumerate(zip(axes, SIX_SEGMENTS, ext)):
        cx, cy = (e[0] + e[2]) / 2 + 250, (e[1] + e[3]) / 2
        b = (cx - w / 2, cy - h / 2, cx + w / 2, cy + h / 2)
        island.clip(box(*b)).plot(ax=ax, color="0.93", edgecolor="0.62", lw=0.4, zorder=1)
        seg = boxes[boxes["gis"].between(lo, hi)]
        seg.boundary.plot(ax=ax, color="0.78", lw=0.3, zorder=2)
        # A grey edge under the colour, so agreeing (white) stretches stay visible
        ax.add_collection(LineCollection(seg_xy, colors="0.45", lw=3.3, zorder=3,
                                         capstyle="round"))
        lc = LineCollection(seg_xy, cmap=DIFF_CMAP, norm=norm, lw=2.4, zorder=4,
                            capstyle="round")
        lc.set_array(seg_c)
        ax.add_collection(lc)
        for _, r in seg.iterrows():
            if r.gis % 5 == 0 or r.gis in (lo, hi):
                ax.text(r.geometry.bounds[2] + 80, r.geometry.centroid.y, str(r.gis),
                        fontsize=6, color=fs.INK_MUTED, va="center", ha="left", zorder=6)
        ax.set_xlim(b[0], b[2])
        ax.set_ylim(b[1], b[3])
        ax.set_aspect("equal")
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_title("({0})  GIS {1}–{2}".format(chr(97 + i), lo, hi), loc="left",
                     fontsize=8.5, pad=3)
    fs._scalebar(axes[0], 1000.0, show_cells=False)
    fs._north_arrow(axes[0], x=0.12, y=0.84, length=0.03)
    sm = plt.cm.ScalarMappable(cmap=DIFF_CMAP, norm=norm)
    cb = fig.colorbar(sm, ax=list(axes), location="bottom", shrink=0.5, aspect=35, pad=0.01)
    cb.set_label("2-yr minus 3-yr mean shoreline (m); blue seaward, red landward")
    cb.outline.set_linewidth(0.5)
    fs.caption(fig, (
        "Where the {p} period's CoastSat mean shoreline moves when the window "
        "changes from the calendar {c} to {lo} – {hi} (±1 yr of the {s}), along "
        "the whole island in six north-up segments of 15 Barrier3D domains, at "
        "one common scale: (a) GIS 1–15 from Cape Point, through (f) GIS 76–90 "
        "at the north end. The line is the 3-yr mean shoreline as it lies on the "
        "ground, coloured between each pair of CoastSat transects (~50 m) by the "
        "2-yr minus 3-yr difference: blue where the DEM-centred line is seaward, "
        "red where it is landward, white where they agree. At this scale 10 m is "
        "a fraction of the line's width (edged in grey so it stays visible where "
        "the difference is near zero), so the two lines are not drawn separately; "
        "mean_shoreline_windows_{p}_lines_*.png show them apart. Island outline "
        "pale grey, the 500 m model domains as thin boxes, every fifth numbered on "
        "the ocean side. Island mean {m:+.1f} m.".format(
            p=period, c=cal.span, lo=dem.lo_iso, hi=dem.hi_iso,
            s=SURVEYS[dem.anchor]["survey"], m=d["diff_m"].mean())))
    out = fs.save(fig, folder / "mean_shoreline_windows_{0}_island.png".format(period),
                  vector=False, close=True)
    print("  figure -> {0}".format(out[0].name))


# The two lines as the offset build sees them (Hannah, 2026-09-29

OFFSET_PRODUCER = (_REPO / "scripts" / "input_prep" / "2-brie-offset" / "1-produce"
                   / "duneline_to_raw_offsets.py")
# SECTIONS and GRID_COLS as compare_offset_sources.py sets them (its 2 x 3 note
SECTIONS = ((1, 15), (16, 30), (31, 45), (46, 60), (61, 75), (76, 90))
GRID_COLS = 3


# duneline_to_raw_offsets.py, loaded as a module
def _offset_producer():
    import importlib.util
    spec = importlib.util.spec_from_file_location("duneline_to_raw_offsets", OFFSET_PRODUCER)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


# Per domain, each window's mean station from the offshore datum
def datum_stations(windows):
    import geopandas as gpd
    from site_layer.hat_observed_rates import mean_shoreline_geojson
    prod = _offset_producer()
    tr = prod.load_transects()
    out = {}
    for w in windows:
        line = gpd.read_file(mean_shoreline_geojson(*w.key)).to_crs(tr.crs).iloc[0].geometry
        raw = prod.intersect(tr, line)
        out[w.label] = raw.groupby("domain_id")["ORIG_LEN"].mean()
    return pd.DataFrame(out)


# Both windows in the offset producer's datum
def datum_figure(cal, dem, period, folder):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import matplotlib.ticker as mticker

    st = datum_stations([cal, dem])
    s3, s2 = st[cal.label], st[dem.label]
    # Stations grow LANDWARD, so 3-yr minus 2-yr is + where the 2-yr line is seaward
    diff = s3 - s2
    st["seaward_shift_2yr_m"] = diff
    st.index.name = "gis_domain"
    st.to_csv(fs.support_dir(folder) / "datum_stations_{0}.csv".format(period),
              float_format="%.2f")
    v1_raw = None
    try:
        from site_layer import hat_topo_version as tv
        # the calendar window built shoreline v1, which keeps its raw file
        v1_raw = (tv.offset_build_dir(period, "v1", "shoreline")
                  / "{0}_shoreline_offset_raw.csv".format(cal.label))
        chk = (pd.read_csv(v1_raw).drop_duplicates(["domain_id", "LineID"])
               .groupby("domain_id")["ORIG_LEN"].mean())
        print("  3-yr stations vs the stored v1 raw file: max |diff| {0:.3f} m".format(
            float((chk - s3).abs().max())))
    except Exception as exc:                       # a check, never a blocker
        print("  (v1 raw file not checked: {0})".format(exc))

    fs.apply_style()
    # Each section a pair: the two profiles and a strip of their difference
    n_rows = int(np.ceil(len(SECTIONS) / GRID_COLS))
    lim = float(np.ceil(max(np.nanmax(np.abs(diff)), fs.CELL_M) / 5.0) * 5.0) + 2.0
    fig = plt.figure(figsize=(fs.FIG_W_DOUBLE, 3.9 * n_rows + 1.2), layout="constrained")
    gs = fig.add_gridspec(n_rows, 2 * GRID_COLS, width_ratios=[2.6, 1.0] * GRID_COLS)
    fill_sea, fill_land = fs.C_1997_FILL, fs.C_1984_FILL     # RdBu: blue seaward
    for i, (lo, hi) in enumerate(SECTIONS):
        r, c = divmod(i, GRID_COLS)
        ax = fig.add_subplot(gs[r, 2 * c])
        sx = fig.add_subplot(gs[r, 2 * c + 1], sharey=ax)
        sl = (st.index >= lo) & (st.index <= hi)
        dom = st.index.to_numpy()[sl]
        va, vb, dd = s3.to_numpy()[sl], s2.to_numpy()[sl], diff.to_numpy()[sl]

        _town_bands_y(ax, lo, hi)
        ax.plot(va, dom, color=C_CAL, lw=1.3, zorder=4,
                label="3-yr, calendar {0} (offset v1)".format(cal.span))
        ax.plot(vb, dom, color=C_DEM, lw=1.3, ls=(0, (4, 2)), zorder=5,
                label="2-yr, {0} (offset v2)".format(dem.span))
        ax.set_ylim(lo - 0.5, hi + 0.5)
        ax.invert_xaxis()          # landward left, ocean right
        ax.xaxis.set_major_locator(mticker.MaxNLocator(3))
        ax.yaxis.set_major_locator(mticker.MaxNLocator(integer=True))
        ax.tick_params(axis="x", labelsize=7)
        fs._title(ax, i, "GIS {0}–{1}".format(lo, hi))

        # the strip: + = 2-yr seaward, drawn to the RIGHT like the ocean
        for s_lo, s_hi in _town_spans(lo, hi):
            sx.axhspan(s_lo - 0.5, s_hi + 0.5, color="0.94", lw=0, zorder=0)
        for xc in (-fs.CELL_M, fs.CELL_M):          # one Barrier3D cell either way
            sx.axvline(xc, color=fs.INK_MUTED, lw=0.5, ls=(0, (1, 2)), zorder=1)
        sx.axvline(0, color=fs.INK, lw=0.6, zorder=1)
        sx.fill_betweenx(dom, 0, dd, where=dd >= 0, interpolate=True, color=fill_sea,
                         lw=0, zorder=2, label="2-yr seaward of 3-yr")
        sx.fill_betweenx(dom, 0, dd, where=dd < 0, interpolate=True, color=fill_land,
                         lw=0, zorder=2, label="2-yr landward of 3-yr")
        sx.plot(dd, dom, color=C_DEM, lw=0.9, zorder=3)
        sx.set_xlim(-lim, lim)
        sx.xaxis.set_major_locator(mticker.FixedLocator([-10, 0, 10]))
        sx.tick_params(axis="x", labelsize=7)
        sx.tick_params(axis="y", labelleft=False)
        sx.set_xlabel("2-yr − 3-yr (m)", fontsize=7)
    fig.supylabel(fs.DOMAIN_AXIS_LABEL, fontsize=9)
    fig.supxlabel("Distance from the offshore datum (m)   —   landward ←   |   → ocean",
                  fontsize=9)
    h1, l1 = fig.axes[0].get_legend_handles_labels()
    h2, l2 = fig.axes[1].get_legend_handles_labels()
    fig.legend(h1 + h2, l1 + l2, loc="outside upper center", ncol=2, fontsize=7,
               frameon=False)
    fs.caption(fig, (
        "The {p} period's CoastSat mean shoreline over two windows, as the "
        "island-offset build sees it: the calendar {c} mean (solid grey, the "
        "line behind the shoreline offset v1) and the {lo} to {hi} mean, ±1 yr "
        "of the {s} (dashed purple, the line behind v2), in six consecutive "
        "sections of Hatteras on a 2 x 3 grid, alongshore up the page, south at "
        "the bottom. In each wide panel both are the distance from the shared "
        "offshore datum along the 100 m model transects, averaged per 500 m "
        "domain, with the ocean on the right; they are computed with the offset "
        "build's own intersection (duneline_to_raw_offsets.intersect), so the two "
        "go through identical code. At this scale the two profiles print as one "
        "line, so the narrow strip beside each panel draws the difference itself "
        "on the same domain axis: 2-yr minus 3-yr, blue to the right where the "
        "2-yr line is seaward, red to the left where it is landward, the dotted "
        "lines ±1 Barrier3D cell (10 m); every strip shares one scale. Island "
        "mean {m:+.1f} m; domains beyond one cell: {n10} of {n}. Village spans "
        "shaded. Per-domain numbers in datum_stations_{p}.csv.".format(
            p=period, c=cal.span, lo=dem.lo_iso, hi=dem.hi_iso,
            s=SURVEYS[dem.anchor]["survey"], m=diff.mean(),
            n10=int((diff.abs() > fs.CELL_M).sum()), n=len(diff))))
    out = fs.save(fig, folder / "mean_shoreline_windows_{0}_datum.png".format(period),
                  close=True)
    print("  figure -> {0}".format(out[0].name))


# The village spans that touch GIS lo..hi, from the one owner
def _town_spans(lo, hi):
    try:
        from site_layer.hatteras_site_config import HATTERAS_ANNOTATIONS
    except ImportError:
        return []
    return [(a, z) for a, z in HATTERAS_ANNOTATIONS.town_spans.values()
            if not (z + 0.5 < lo - 0.5 or a - 0.5 > hi + 0.5)]


# Village spans against a vertical alongshore axis, as compare_offset_sources draws them
def _town_bands_y(ax, lo, hi):
    try:
        from site_layer.hatteras_site_config import HATTERAS_ANNOTATIONS
        spans = HATTERAS_ANNOTATIONS.town_spans
    except ImportError:
        return
    for name, (s_lo, s_hi) in spans.items():
        if s_hi + 0.5 < lo - 0.5 or s_lo - 0.5 > hi + 0.5:
            continue
        ax.axhspan(s_lo - 0.5, s_hi + 0.5, color="0.94", lw=0, zorder=0)
        mid = (max(s_lo - 0.5, lo - 0.5) + min(s_hi + 0.5, hi + 0.5)) / 2
        ax.text(0.015, mid, name, transform=ax.get_yaxis_transform(), ha="left",
                va="center", fontsize=6.5, color=fs.INK_MUTED, zorder=1, clip_on=True)


# The largest differences on the photographs
def draw_lines(d, cal, dem, period, folder):
    largest = (d["diff_minus_island_mean_m"].abs().sort_values(ascending=False)
               .head(N_LARGEST).index)
    import coastsat_mean_shoreline_on_imagery as moi
    lines_figure(d, cal, dem, period, folder, sorted(int(g) for g in largest),
                 "largest",
                 "The {0} domains whose difference is largest once the island "
                 "mean ({1:+.1f} m) is removed, south to north.".format(
                     N_LARGEST, d["diff_m"].mean()))
    lines_figure(d, cal, dem, period, folder, [c for c, _ in moi.DEFAULT_SITES],
                 "sites",
                 "The centre domains of the six sites of the mean_shoreline "
                 "imagery figures, south to north.")


# Run: every period asked for
def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    ap.add_argument("--periods", nargs="+", type=int, default=sorted(PERIODS),
                    choices=sorted(PERIODS))
    ap.add_argument("--no-lines", action="store_true",
                    help="skip the line-on-photograph figures (they need D: and rasterio)")
    a = ap.parse_args(argv)
    built_on = datetime.now(timezone.utc).strftime("%Y-%m-%d")
    for period in a.periods:
        cal, dem = PERIODS[period]
        print("{0}: {1} vs {2}".format(period, cal.label, dem.label))
        t, d, only = compare(cal, dem)
        folder = MEAN_SHORELINE_WINDOWS / str(period)
        sup = fs.support_dir(folder)
        t.to_csv(sup / "transect_comparison_{0}.csv".format(period))
        d.to_csv(sup / "domain_comparison_{0}.csv".format(period))
        figure(t, d, cal, dem, period, folder)
        island_figure(t, d, cal, dem, period, folder)
        datum_figure(cal, dem, period, folder)
        if not a.no_lines:
            draw_lines(d, cal, dem, period, folder)
        write_provenance(t, d, only, cal, dem, period, folder, built_on)
        print("  island mean {0:+.2f} m, alongshore sd {1:.2f} m".format(
            d["diff_m"].mean(), d["diff_minus_island_mean_m"].std()))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
