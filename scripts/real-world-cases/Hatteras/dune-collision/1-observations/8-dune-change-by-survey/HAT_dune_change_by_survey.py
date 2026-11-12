"""
What the dunes did, survey by survey, 1996-2009, in the stretches where the collision module misses the measured crest
drops: north of Buxton (GIS 10-15), Avon south (GIS 21-26) and Rodanthe (GIS 82-89). And are those drops real dune loss,
or a picking artifact of the survey crests the model is scored on?

    python scripts/real-world-cases/Hatteras/dune-collision/1-observations/8-dune-change-by-survey/HAT_dune_change_by_survey.py [region ...]

Regions: buxton, avon_south, rodanthe (site_layer/hat_lidar_strips.py REGIONS; default all three). The crest check is
drawn when all three are run.

Merged 2026-11-12 from 8-buxton-dune-change (2026-11-11: asked because the module cannot lower GIS 10-12 in time;
dune-collision/4-hindcast/TESTING_2026-11-09_to_11.md, section 4) and 9-hotspot-dune-change (2026-11-12: the hotspots
of the hotspot erosion test, sections 17-18). The surveys, the datum and the unrolled frame are in
site_layer/hat_lidar_strips.py.

The crest check: the survey crests (input_prep/1-barrier3d-domains/2-survey-crests/HAT_survey_dune_crests.py) are each
profile's highest point inside the start's hand-picked dune window, shifted landward with the shore. If the dune moves
landward out of that window, or a new one is built elsewhere, the pick lands on the old dune's remnant and reads low.
Here, independently of the windows: for every alongshore line, the highest lidar point between NC-12 (+8 m) and the
beach (the first point below mean high water), and its distance seaward of the road.

Writes, in the workspace's dune-collision/1-observations/8-dune-change-by-survey/, per region in its own subfolder
(1-buxton/, 2-avon-south/, 3-rodanthe/, south to north):
- <region>_dunes_by_survey.png: the coast unrolled along NC-12 smoothed over about 1 km, one strip per survey (terrain colours, hillshade, the
  cross-shore stretched), the berm (1.7 m) contour, NC-12 where it was, the domain boundaries
- <region>_dune_change_between_surveys.png: each survey minus the one before
- <region>_dune_profiles.png: one cross-shore profile per domain (the 100 m transect nearest its middle), every survey,
  with the storm series' 1998 and 1999 peak water levels (Bonnie 2.64 m, Dennis 3.65 m NAVD88)
- <region>_aerial_by_year.png: the aerial photographs of 1993-2009 in the same frame
- <region>_plan_lidar.png, <region>_plan_aerial.png: the surveys and the photographs as north-up maps, true geometry
- <region>_plan_dune_change.png: each survey minus the one before, as north-up maps
and across them:
- hotspot_crest_check.png: per domain, the picked crest against the highest point between NC-12 and the beach
- supporting/hotspot_crest_check.csv, supporting/datum_check.csv (each survey's difference from 2009 on NC-12 per
  region, after the corrections), supporting/imagery_2002_unzipped/ (the 2002 frames, unzipped from the USGS delivery)

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-11-12
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

_HERE = Path(__file__).resolve()
_HAT_DIR = next(p for p in _HERE.parents if (p / "hatteras_paths.py").exists())
sys.path.insert(0, str(_HAT_DIR))
from hatteras_paths import WORKSPACE  # noqa: E402
from site_layer import hat_lidar_strips as L  # noqa: E402

# --- CONFIG ------------------------------------------------------------------
OUT = WORKSPACE / "dune-collision" / "1-observations" / "8-dune-change-by-survey"
IMAGERY_CACHE = OUT / "supporting" / "imagery_2002_unzipped"
# One subfolder per region, south to north (2026-11-12); the island-wide crest check stays at the top
REGION_ORDER = ("buxton", "avon_south", "rodanthe")


def region_dir(r):
    return OUT / f"{REGION_ORDER.index(r.name) + 1}-{r.name.replace('_', '-')}"
SEAWARD_OF_ROAD_M = 8.0          # the independent crest: searched from this far seaward of NC-12 (off the shoulder)
PLAN_LAND_M, PLAN_SEA_EXTRA_M = 120.0, 20.0   # the plan maps: landward of NC-12, beyond the strips' reach
PLAN_AERIAL_RES_M = 1.0
THICK = (1998, 1999, 2004)       # the profiles drawn thicker: before and after Dennis, and after Isabel
CREST_TABLE = (WORKSPACE / "data" / "hatteras_init" / "1-barrier3d-domains" / "1996-start" / "survey-crests"
               / "survey_dune_crests.csv")
DATUM_NOTE = ("every survey on the 2009 lidar's datum (the model's vertical offsets, then each survey's remaining "
              "difference from 2009 on NC-12 removed; supporting/datum_check.csv)")
# -----------------------------------------------------------------------------


def _title(r):
    return f"GIS {r.gis[0]}-{r.gis[-1]}, {r.label}"


def strips(r, surv, st, S_, plt):
    from matplotlib.colors import LightSource
    from matplotlib.lines import Line2D
    years = sorted(surv)
    ls = LightSource(azdeg=315, altdeg=45)
    cmap = plt.get_cmap("terrain")
    fig, axs = plt.subplots(len(years), 1, figsize=(S_.FIG_W_DOUBLE, 0.82 * len(years) + 1.3), sharex=True,
                            layout="constrained")
    for ax, y in zip(axs, years):
        z = st.sample(surv[y][0])
        zz = np.where(np.isfinite(z), z, L.TERRAIN_RANGE[0])
        rgb = ls.shade(np.clip(zz, *L.TERRAIN_RANGE), cmap=cmap, vert_exag=1, blend_mode="soft",
                       vmin=L.TERRAIN_RANGE[0], vmax=L.TERRAIN_RANGE[1], dx=L.STRIP_STEP_M, dy=L.STRIP_STEP_M)
        rgb[~np.isfinite(z)] = 1.0
        ax.imshow(rgb, extent=st.extent(), origin="lower", aspect="auto", interpolation="nearest")
        ax.contour(st.s / 1000, st.n, np.where(np.isfinite(z), z, np.nan), levels=[L.BERM_M], colors="white",
                   linewidths=0.4)
        st.draw(ax, y)
        L.row_label(ax, f"{y}{'*' if y == 2001 else ''}", L.MONTHS.get(y))
    st.top_axis(axs[0])
    axs[-1].set_xlabel(L.STRIP_XLABEL, fontsize=8)
    fig.supylabel(L.STRIP_YLABEL, fontsize=8)
    L.colorbar(fig, axs, plt.cm.ScalarMappable(cmap=cmap, norm=plt.Normalize(*L.TERRAIN_RANGE)),
               "Elevation (m NAVD88)")
    _legend(fig, r, (Line2D([], [], color="0.6", lw=0.6, label="Berm, 1.7 m (module toe; white on the map)"),))
    out = region_dir(r) / f"{r.name}_dunes_by_survey.png"
    S_.save(fig, out, bbox_inches="tight", close=True)
    S_.record_caption(out, (
        f"Lidar elevation of {_title(r)}, in every 1996-2009 survey (rows; flight months beneath each year), "
        f"{L.BASELINE_NOTE}. Elevation with hillshade, {DATUM_NOTE}. White contour: the berm, 1.7 m NAVD88, the collision module's toe. "
        f"Magenta: {r.note} Grey vertical lines: the domain boundaries, where NC-12 passes from one domain into the "
        "next, numbered along the top. Survey grids of 1 m (2009) to 9 ft (about 2.7 m), resampled to 1 m. *2001 is "
        "unclassified lidar: buildings and vegetation can stand out as false highs."))
    print(f"wrote {out}")


def _legend(fig, r, extra=()):
    from matplotlib.lines import Line2D
    h = [Line2D([], [], color=L.ROAD_C, lw=1.0, label="NC-12")]
    if r.washed_out_survey:
        h.append(Line2D([], [], color=L.ROAD_C, lw=1.0, ls=(0, (4, 2)),
                        label=f"NC-12, washed out ({r.washed_out_survey})"))
    h += [Line2D([], [], color=L.BOUNDARY_C, lw=0.6, label="Domain boundary"), *extra]
    fig.legend(handles=h, loc="outside lower center", ncol=len(h), fontsize=7, frameon=False, handlelength=2.2,
               columnspacing=1.6)


def strip_changes(r, surv, st, S_, plt):
    years = sorted(surv)
    pairs = list(zip(years[:-1], years[1:]))
    fig, axs = plt.subplots(len(pairs), 1, figsize=(S_.FIG_W_DOUBLE, 0.82 * len(pairs) + 1.3), sharex=True,
                            layout="constrained")
    for ax, (a, b) in zip(axs, pairs):
        za, zb = st.sample(surv[a][0]), st.sample(surv[b][0])
        land = (np.nan_to_num(za, nan=-9) > L.MHW_M) | (np.nan_to_num(zb, nan=-9) > L.MHW_M)
        d = np.where(land & np.isfinite(za) & np.isfinite(zb), zb - za, np.nan)
        im = ax.imshow(d, extent=st.extent(), origin="lower", aspect="auto", cmap="RdBu", vmin=-L.DIFF_MAX,
                       vmax=L.DIFF_MAX, interpolation="nearest")
        st.draw(ax, b)
        L.row_label(ax, f"{a}–{str(b)[2:]}{'*' if 2001 in (a, b) else ''}")
    st.top_axis(axs[0])
    axs[-1].set_xlabel(L.STRIP_XLABEL, fontsize=8)
    fig.supylabel(L.STRIP_YLABEL, fontsize=8)
    L.colorbar(fig, axs, im, "Elevation change between surveys (m; red lower, blue higher)", extend="both")
    _legend(fig, r)
    out = region_dir(r) / f"{r.name}_dune_change_between_surveys.png"
    S_.save(fig, out, bbox_inches="tight", close=True)
    S_.record_caption(out, (
        f"Elevation change of {_title(r)}, between successive lidar surveys (rows), the coast unrolled as in "
        f"{r.name}_dunes_by_survey (along NC-12 smoothed over about 1 km, cross-shore stretched, seaward down), on the 2009 datum, land only (above mean high "
        f"water in either survey); red lower, blue higher, clipped at +-3 m. Magenta: NC-12 at the later survey; "
        f"{r.note} Grey vertical lines: the domain boundaries, numbered along the top. The survey grids differ (1 m to "
        "about 2.7 m), so changes of a few decimetres along steep dune faces are partly resampling; where the 2005 "
        "survey leaves out buildings, the houses show as red. *2001 is unclassified lidar."))
    print(f"wrote {out}")


def profiles(r, surv, g, dom, roads, S_, plt):
    """One cross-shore profile per domain, in the strips' frame: distance from NC-12's 2008 alignment, seaward to the
    right, so the panels line up with each other and with the road positions."""
    import geopandas as gpd
    from matplotlib.lines import Line2D
    from scipy.ndimage import map_coordinates
    t, shape, crs = g
    tr = gpd.read_file(L.TRANSECTS).to_crs(crs)
    d = dom.to_crs(crs)
    new_road = roads["2008"].geometry.union_all()
    old_road = roads["1978"].geometry.union_all() if r.relocated_by_survey and "1978" in roads else None
    years = sorted(surv)
    colours = dict(zip(years, plt.get_cmap("viridis")(np.linspace(0.0, 0.9, len(years)))))
    gl = list(r.gis)
    nrow = -(-len(gl) // 2)
    fig, axs = plt.subplots(nrow, 2, figsize=(S_.FIG_W_DOUBLE, 2.3 * nrow + 0.5), sharex=True, sharey=True,
                            layout="constrained", squeeze=False)
    for k, (ax, gis) in enumerate(zip(axs.ravel(), gl)):
        poly = d[d.domain_id == gis].geometry.iloc[0]
        inside = tr[tr.intersects(poly)]
        line = inside.iloc[(inside.geometry.centroid.y - poly.centroid.y).abs().argsort().iloc[0]].geometry
        c = np.asarray(line.coords, float)
        u = (c[-1] - c[0]) / np.hypot(*(c[-1] - c[0]))         # ocean end to land end
        s = np.arange(0, line.length, 1.0)
        xy = c[0] + s[:, None] * u
        col, row = (xy[:, 0] - t.c) / t.a - 0.5, (xy[:, 1] - t.f) / t.e - 0.5
        at = lambda geom: line.project(line.intersection(geom)) if geom is not None and line.intersects(geom) \
            else np.nan
        s_new, s_old = at(new_road), at(old_road)
        x = s_new - s                                            # seaward of the 2008 road, positive
        for y in years:
            z = map_coordinates(surv[y][0], [row, col], order=1, mode="constant", cval=np.nan)
            ax.plot(x, z, color=colours[y], lw=1.3 if y in THICK else 0.8, zorder=3)
        for _, zw, lstyle in L.STORMS:
            ax.axhline(zw, color="#2171B5", lw=0.7, ls=lstyle, zorder=2)
        ax.axhline(L.BERM_M, color="0.55", lw=0.6, zorder=2)
        ax.axvline(0, color=L.ROAD_C, lw=1.0, zorder=4)
        if np.isfinite(s_old) and abs(s_new - s_old) > 3:
            ax.axvline(s_new - s_old, color=L.ROAD_C, lw=1.0, ls=(0, (4, 2)), zorder=4)
        ax.set_xlim(-r.land_m, r.sea_m)
        ax.set_ylim(-1, 9)
        ax.text(0.02, 0.96, f"({chr(97 + k)}) GIS {gis}", transform=ax.transAxes, fontsize=8, fontweight="bold",
                va="top")
        ax.grid(color="0.93", lw=0.4)
        ax.spines[["top", "right"]].set_visible(False)
        for sp in ("left", "bottom"):
            ax.spines[sp].set_linewidth(0.5)
        ax.tick_params(labelsize=7, length=2.5, width=0.5)
    for ax in axs.ravel()[len(gl):]:
        ax.axis("off")
    for ax in axs[:, 0]:
        ax.set_ylabel("Elevation (m NAVD88)", fontsize=8)
    for ax in axs[-1]:
        ax.set_xlabel("Distance seaward of NC-12 (2008 alignment) (m)", fontsize=8)
    h = [Line2D([], [], color=colours[y], lw=1.3 if y in THICK else 0.8,
                label=f"{y}{'*' if y == 2001 else ''}") for y in years]
    if old_road is not None:
        h += [Line2D([], [], color=L.ROAD_C, lw=1.0, label=f"NC-12 from {r.relocated_by_survey}"),
              Line2D([], [], color=L.ROAD_C, lw=1.0, ls=(0, (4, 2)), label=f"NC-12 to {r.relocated_by_survey - 1}")]
    else:
        h.append(Line2D([], [], color=L.ROAD_C, lw=1.0, label="NC-12"))
    h += [Line2D([], [], color="#2171B5", lw=0.7, ls=":", label="Bonnie peak, 1998 (2.64 m)"),
          Line2D([], [], color="#2171B5", lw=0.7, ls="--", label="Dennis peak, 1999 (3.65 m)"),
          Line2D([], [], color="0.55", lw=0.6, label="Berm, 1.7 m (module toe)")]
    fig.legend(handles=h, loc="outside lower center", ncol=7, fontsize=7, frameon=False, handlelength=2.0,
               columnspacing=1.3)
    out = region_dir(r) / f"{r.name}_dune_profiles.png"
    S_.save(fig, out, bbox_inches="tight", close=True)
    S_.record_caption(out, (
        f"One cross-shore lidar profile per domain, {_title(r)} (the 100 m transect nearest each domain's middle), in "
        "every 1996-2009 survey on the 2009 datum (lighter to darker, 1996 to 2009; 1998, 1999 and 2004 thicker), "
        f"against distance seaward of NC-12's 2008 alignment (the ocean to the right; the same frame as "
        f"{r.name}_dunes_by_survey). Magenta: {r.note} Blue: the storm series' peak water levels in 1998 (Bonnie, "
        "2.64 m NAVD88, dotted) and 1999 (Dennis, 3.65 m, dashed); grey: the berm, 1.7 m, the collision module's toe. "
        "*2001 is unclassified lidar."))
    print(f"wrote {out}")


def imagery_strips(r, st, g, S_, plt):
    """The aerial photographs in the strips' frame, one row per year with a photograph over the region."""
    from matplotlib.lines import Line2D
    rows = L.imagery_rows(st, g, IMAGERY_CACHE)
    fig, axs = plt.subplots(len(rows), 1, figsize=(S_.FIG_W_DOUBLE, 0.82 * len(rows) + 1.0), sharex=True,
                            layout="constrained", squeeze=False)
    axs = axs[:, 0]
    for ax, (y, date, img) in zip(axs, rows):
        img = img.copy()
        img[img.sum(axis=-1) == 0] = 1.0
        ax.imshow(img, extent=st.extent(), origin="lower", aspect="auto", interpolation="bilinear")
        st.draw(ax, y)
        ax.lines[0].set_alpha(0.45)                    # the road line see-through: the road in the photograph shows
        L.row_label(ax, str(y), L.IMAGERY_LABELS.get(y, date[5:] if len(date) > 4 else None))
    st.top_axis(axs[0])
    axs[-1].set_xlabel(L.STRIP_XLABEL, fontsize=8)
    fig.supylabel(L.STRIP_YLABEL, fontsize=8)
    fig.legend(handles=[Line2D([], [], color=L.ROAD_C, lw=1.0, label="NC-12" + (
        " (relocation records)" if r.relocated_by_survey else " (2008 alignment)")),
                        Line2D([], [], color=L.BOUNDARY_C, lw=0.6, label="Domain boundary")],
               loc="outside lower center", ncol=2, fontsize=7, frameon=False, handlelength=2.2)
    out = region_dir(r) / f"{r.name}_aerial_by_year.png"
    S_.save(fig, out, bbox_inches="tight", close=True)
    years = ", ".join(str(x[0]) for x in rows)
    S_.record_caption(out, (
        f"Aerial photographs of {_title(r)} ({years}; rows, flight dates beside each year), in the frame of "
        f"{r.name}_dunes_by_survey: {L.BASELINE_NOTE} (north-up maps: {r.name}_plan_aerial). Sources: "
        "1993-1998 and 2002, USGS georeferenced scans of USACE Field Research Facility photography (Henderson, in "
        "prep., doi:10.5066/P1CXBCDW; 2002 from the draft delivery); 2004 and 2009, Google Earth Pro historical "
        "imagery captured and georeferenced by Hannah Henry; 2008, NOAA NGS ortho-rectified colour mosaic (InPort "
        f"48695); years with no photograph over this stretch are left out. Magenta, see-through: {r.note} Grey: the "
        "domain boundaries. White: no photograph. Georeferencing errors of a few metres (USGS: about 1.2 m) can offset "
        "a photograph against the road line."))
    print(f"wrote {out}")


def _plan_window(r, g, roads):
    """The plan maps' window: the region's rows of the grid, and from PLAN_LAND_M landward of the westernmost NC-12
    (either alignment) to the strip's reach plus PLAN_SEA_EXTRA_M seaward of the easternmost (the ocean is east)."""
    from shapely.geometry import box
    t, (h, w), _ = g
    win = box(t.c, t.f + h * t.e, t.c + w * t.a, t.f)
    xs = [c for name in roads for geom in roads[name].geometry if not geom.intersection(win).is_empty
          for c in (lambda b: (b[0], b[2]))(geom.intersection(win).bounds)]
    return min(xs) - PLAN_LAND_M, max(xs) + r.sea_m + PLAN_SEA_EXTRA_M, t.f + h * t.e + L.PAD_M, t.f - L.PAD_M


# The plan maps share one fixed layout per region (2026-11-12), so they stack as layers on slides: the same figure
# size, PLAN_COLS panels at the same positions, the colour bar and legend in the same places (an empty slot or no colour
# bar is left blank, not closed up). The change maps sit in their later survey's slot, under that survey's map
PLAN_COLS = 9
PLAN_IN = dict(left=0.12, right=0.12, gap=0.07, top=0.42, bottom=1.05)       # inches


def _plan_axes(S_, plt, window):
    """The fixed layout: figure, the PLAN_COLS panel axes (south-north map panels, left to right) and the colour-bar
    axes; nothing laid out by matplotlib, and saved without trimming, so every plan figure of a region matches."""
    x0, x1, y0, y1 = window
    m = PLAN_IN
    W = S_.FIG_W_DOUBLE
    pw = (W - m["left"] - m["right"] - m["gap"] * (PLAN_COLS - 1)) / PLAN_COLS
    ph = pw * (y1 - y0) / (x1 - x0)
    H = m["top"] + ph + m["bottom"]
    fig = plt.figure(figsize=(W, H))
    axs = [fig.add_axes([(m["left"] + k * (pw + m["gap"])) / W, m["bottom"] / H, pw / W, ph / H])
           for k in range(PLAN_COLS)]
    for ax in axs:
        ax.axis("off")                                   # switched back on by the panels that are drawn
    cax = fig.add_axes([0.275, 0.62 / H, 0.45, 0.09 / H])
    cax.axis("off")
    return fig, axs, cax


def _plan_panel(ax):
    ax.axis("on")


def _plan_colorbar(fig, cax, mappable, label, extend="neither"):
    cax.axis("on")
    cb = fig.colorbar(mappable, cax=cax, orientation="horizontal", extend=extend)
    cb.set_label(label, fontsize=8)
    cb.ax.tick_params(labelsize=7, width=0.5, length=2.5)
    cb.outline.set_linewidth(0.5)


def _plan_legend(fig, handles):
    H = fig.get_size_inches()[1]
    fig.legend(handles=handles, loc="center", bbox_to_anchor=(0.5, 0.17 / H), ncol=len(handles), fontsize=7,
               frameon=False, handlelength=2.2, columnspacing=1.6)


def _plan_save(S_, fig, out):
    S_.save(fig, out, dpi=300, close=True)                 # no trimming: every plan figure keeps the same frame


def _legend_handles(r, extra=()):
    from matplotlib.lines import Line2D
    h = [Line2D([], [], color=L.ROAD_C, lw=1.0, label="NC-12")]
    if r.washed_out_survey:
        h.append(Line2D([], [], color=L.ROAD_C, lw=1.0, ls=(0, (4, 2)), label=f"NC-12, washed out ({r.washed_out_survey})"))
    return h + [Line2D([], [], color=L.BOUNDARY_C, lw=0.6, label="Domain boundary"), *extra]


def _plan_decorate(ax, r, year, dom, roads, window, crs, first):
    """Domain outlines (numbered on the first panel), NC-12 where it was that year, north arrow and scale bar."""
    x0, x1, y0, y1 = window
    for _, row in dom.to_crs(crs).iterrows():
        x, y = row.geometry.exterior.xy
        ax.plot(x, y, color=L.BOUNDARY_C, lw=0.4, zorder=5)
        if first:
            cy = row.geometry.centroid.y
            if y0 < cy < y1:
                ax.text(x0 + 0.04 * (x1 - x0), cy, f"{int(row.domain_id)}", fontsize=6.5, va="center",
                        color="0.1", fontweight="bold", zorder=9,
                        bbox=dict(boxstyle="round,pad=0.15", fc="white", ec="none", alpha=0.7))
    old = r.relocated_by_survey is not None and year < r.relocated_by_survey and "1978" in roads
    for geom in roads["1978" if old else "2008"].geometry:
        for part in getattr(geom, "geoms", [geom]):
            x, y = part.xy
            ax.plot(x, y, color=L.ROAD_C, lw=0.9, ls=(0, (4, 2)) if year == r.washed_out_survey else "-", zorder=8)
    ax.set_xlim(x0, x1)
    ax.set_ylim(y0, y1)
    ax.set_aspect("equal")
    ax.set_xticks([])
    ax.set_yticks([])
    for sp in ax.spines.values():
        sp.set_linewidth(0.5)
    if first:
        bar = 500.0
        xb, yb = x1 - 0.10 * (x1 - x0), y0 + 0.03 * (y1 - y0)
        ax.plot([xb, xb], [yb, yb + bar], color="0.1", lw=1.5, zorder=9, solid_capstyle="butt")
        ax.text(xb - 0.04 * (x1 - x0), yb + bar / 2, "500 m", fontsize=6, rotation=90, va="center", ha="right",
                zorder=9, bbox=dict(boxstyle="round,pad=0.1", fc="white", ec="none", alpha=0.7))
        ax.annotate("N", xy=(x1 - 0.12 * (x1 - x0), y1 - 0.02 * (y1 - y0)),
                    xytext=(x1 - 0.12 * (x1 - x0), y1 - 0.09 * (y1 - y0)), ha="center", va="center", fontsize=7,
                    fontweight="bold", arrowprops=dict(arrowstyle="-|>", color="0.1", lw=0.8), zorder=9)


def plan_lidar(r, surv, g, dom, roads, S_, plt):
    """North-up maps of every survey, side by side: true geometry, no stretch."""
    from matplotlib.colors import LightSource
    from matplotlib.lines import Line2D
    t, shape, crs = g
    window = _plan_window(r, g, roads)
    extent = (t.c, t.c + shape[1] * t.a, t.f + shape[0] * t.e, t.f)
    years = sorted(surv)
    ls = LightSource(azdeg=315, altdeg=45)
    cmap = plt.get_cmap("terrain")
    fig, axs, cax = _plan_axes(S_, plt, window)
    for k, (ax, y) in enumerate(zip(axs, years)):
        _plan_panel(ax)
        z = surv[y][0]
        zz = np.where(np.isfinite(z), z, L.TERRAIN_RANGE[0])
        rgb = ls.shade(np.clip(zz, *L.TERRAIN_RANGE), cmap=cmap, vert_exag=1, blend_mode="soft",
                       vmin=L.TERRAIN_RANGE[0], vmax=L.TERRAIN_RANGE[1], dx=L.RES_M, dy=L.RES_M)
        rgb[~np.isfinite(z)] = 1.0
        ax.imshow(rgb, extent=extent, interpolation="nearest")
        ax.contour(np.where(np.isfinite(z), z, np.nan), levels=[L.BERM_M], colors="white", linewidths=0.3,
                   extent=extent, origin="upper")
        _plan_decorate(ax, r, y, dom, roads, window, crs, k == 0)
        ax.set_title(f"{y}{'*' if y == 2001 else ''}\n{L.MONTHS.get(y, '')}", fontsize=7)
    _plan_colorbar(fig, cax, plt.cm.ScalarMappable(cmap=cmap, norm=plt.Normalize(*L.TERRAIN_RANGE)),
                   "Elevation (m NAVD88)")
    _plan_legend(fig, _legend_handles(r, (Line2D([], [], color="0.6", lw=0.6, label="Berm, 1.7 m (white on the map)"),)))
    out = region_dir(r) / f"{r.name}_plan_lidar.png"
    _plan_save(S_, fig, out)
    S_.record_caption(out, (
        f"Lidar elevation of {_title(r)}, in every 1996-2009 survey (columns; flight months under each year), as maps "
        f"with north up and true geometry (no stretch; the strips of {r.name}_dunes_by_survey show the same surveys "
        f"unrolled along the coast). Elevation with hillshade, {DATUM_NOTE}. White contour: the berm, 1.7 m NAVD88. "
        f"Magenta: {r.note} Grey: the domain outlines, numbered on the first map. *2001 is unclassified lidar."))
    print(f"wrote {out}")


def plan_changes(r, surv, g, dom, roads, S_, plt):
    """North-up maps of each survey minus the one before: the plan view of strip_changes."""
    t, shape, crs = g
    window = _plan_window(r, g, roads)
    extent = (t.c, t.c + shape[1] * t.a, t.f + shape[0] * t.e, t.f)
    years = sorted(surv)
    pairs = list(zip(years[:-1], years[1:]))
    fig, axs, cax = _plan_axes(S_, plt, window)
    for k, (a, b) in enumerate(pairs):
        ax = axs[years.index(b)]                                # under the later survey's map in the lidar figure
        _plan_panel(ax)
        za, zb = surv[a][0], surv[b][0]
        land = (np.nan_to_num(za, nan=-9) > L.MHW_M) | (np.nan_to_num(zb, nan=-9) > L.MHW_M)
        d = np.where(land & np.isfinite(za) & np.isfinite(zb), zb - za, np.nan)
        ax.set_facecolor("0.93")                                # no data / water, so white stays "no change"
        im = ax.imshow(d, extent=extent, cmap="RdBu", vmin=-L.DIFF_MAX, vmax=L.DIFF_MAX, interpolation="nearest")
        _plan_decorate(ax, r, b, dom, roads, window, crs, k == 0)
        ax.set_title(f"{a}–{str(b)[2:]}{'*' if 2001 in (a, b) else ''}\n ", fontsize=7)   # two lines, as the lidar titles
    _plan_colorbar(fig, cax, im, "Elevation change between surveys (m; red lower, blue higher)", extend="both")
    _plan_legend(fig, _legend_handles(r))
    out = region_dir(r) / f"{r.name}_plan_dune_change.png"
    _plan_save(S_, fig, out)
    S_.record_caption(out, (
        f"Elevation change of {_title(r)} between successive lidar surveys (columns), as maps with north up and true "
        f"geometry (the plan view of {r.name}_dune_change_between_surveys), each in the slot of its later survey's map "
        f"in {r.name}_plan_lidar (the first slot left empty), on the 2009 datum, land only (above mean "
        "high water in either survey); red lower, blue higher, clipped at ±3 m; light grey, water or no data. "
        f"Magenta: NC-12 at the later survey; {r.note} Grey: the domain outlines, numbered on the first map. The "
        "survey grids differ (1 m to about 2.7 m), so changes of a few decimetres along steep dune faces are partly "
        "resampling; where the 2005 survey leaves out buildings, the houses show as red. *2001 is unclassified lidar."))
    print(f"wrote {out}")


def plan_aerial(r, g, dom, roads, S_, plt):
    """North-up maps of each year's photographs, side by side."""
    from matplotlib.lines import Line2D
    t, shape, crs = g
    window = _plan_window(r, g, roads)
    x0, x1, y0, y1 = window
    rows = []
    for y in L.IMAGERY_YEARS:
        im = L.imagery(y, IMAGERY_CACHE)
        if not im.files:
            continue
        rgb = im.read(x0, x1, y0, y1, PLAN_AERIAL_RES_M, crs)
        if (rgb.sum(axis=-1) > 0).mean() < 0.05:
            continue
        rows.append((y, im.date, rgb))
    fig, axs, _ = _plan_axes(S_, plt, window)
    for k, (ax, (y, date, rgb)) in enumerate(zip(axs, rows)):
        _plan_panel(ax)
        img = np.clip(rgb / 255.0, 0, 1)
        img[rgb.sum(axis=-1) == 0] = 1.0
        ax.imshow(img, extent=(x0, x1, y0, y1), interpolation="bilinear")
        _plan_decorate(ax, r, y, dom, roads, window, crs, k == 0)
        for ln in ax.lines:
            if ln.get_color() == L.ROAD_C:
                ln.set_alpha(0.5)
        ax.set_title(f"{y}\n{L.IMAGERY_LABELS.get(y, date[5:] if len(date) > 4 else '')}", fontsize=7)
    _plan_legend(fig, [Line2D([], [], color=L.ROAD_C, lw=1.0, label="NC-12"),
                       Line2D([], [], color=L.BOUNDARY_C, lw=0.6, label="Domain outline")])
    out = region_dir(r) / f"{r.name}_plan_aerial.png"
    _plan_save(S_, fig, out)
    S_.record_caption(out, (
        f"Aerial photographs of {_title(r)} ({', '.join(str(x[0]) for x in rows)}; columns, flight dates under each "
        f"year), as maps with north up and true geometry (no stretch; {r.name}_aerial_by_year shows them unrolled "
        "along the coast). Sources as in that figure: USGS georeferenced scans of USACE Field Research Facility "
        "photography 1993-1998 and 2002 (Henderson, in prep., doi:10.5066/P1CXBCDW), Google Earth Pro captures "
        "georeferenced by Hannah Henry (2004, 2009), NOAA NGS ortho-rectified mosaic (2008). Magenta, see-through: "
        f"{r.note} Grey: the domain outlines, numbered on the first map. White: no photograph."))
    print(f"wrote {out}")


def independent_crest(st, surv):
    """Per survey and alongshore line: the highest point from SEAWARD_OF_ROAD_M seaward of NC-12 (2008) to the beach
    (the first point below MHW seaward of it), and its distance seaward of the road."""
    out = {}
    n = st.n
    for y, (z, *_) in surv.items():
        zs = st.sample(z)                                   # (cross-shore, alongshore)
        zc, dist = np.full(zs.shape[1], np.nan), np.full(zs.shape[1], np.nan)
        for j in range(zs.shape[1]):
            if not np.isfinite(st.road08[j]):
                continue
            ok = n >= st.road08[j] + SEAWARD_OF_ROAD_M
            col = zs[:, j]
            below = np.where(ok & (np.nan_to_num(col, nan=99) < L.MHW_M))[0]
            end = below[0] if below.size else len(n)
            idx = np.where(ok[:end] & np.isfinite(col[:end]))[0]
            if idx.size == 0:
                continue
            k = idx[np.nanargmax(col[idx])]
            zc[j], dist[j] = col[k], n[k] - st.road08[j]
        out[y] = (zc, dist)
    return out


def crest_rows(r, st, surv, picked):
    indep = independent_crest(st, surv)
    rows = []
    for gis, (a, b) in sorted(st.spans.items()):
        m = (st.s >= a) & (st.s <= b)
        for y, (zc, dist) in indep.items():
            p = picked[(picked.gis == gis) & (picked.survey == y)]
            rows.append(dict(region=r.label, gis=gis, survey=y,
                             picked_crest_m=float(p.crest_m_navd88.median()) if len(p) else np.nan,
                             share_at_landward_edge=float(p.at_landward_edge.mean()) if len(p) else np.nan,
                             window_shift_m=float(p.window_shift_m.median()) if len(p) else np.nan,
                             independent_crest_m=float(np.nanmedian(zc[m])) if np.isfinite(zc[m]).any() else np.nan,
                             independent_crest_seaward_of_road_m=float(np.nanmedian(dist[m]))
                             if np.isfinite(dist[m]).any() else np.nan))
    return rows


def crest_check_figure(t, S_, plt):
    """Per region (rows) and domain (colours): the picked crest against the independent one, and where it lies."""
    from matplotlib.lines import Line2D
    regions = list(dict.fromkeys(t.region))
    fig, axs = plt.subplots(len(regions), 3, figsize=(S_.FIG_W_DOUBLE, 2.5 * len(regions) + 0.8),
                            layout="constrained", sharex=True, squeeze=False)
    cm = plt.get_cmap("viridis")
    for i, label in enumerate(regions):
        rr = t[t.region == label]
        gl = sorted(rr.gis.unique())
        for k, g in enumerate(gl):
            d = rr[rr.gis == g].sort_values("survey")
            c = cm(k / max(len(gl) - 1, 1) * 0.9)
            axs[i, 0].plot(d.survey, d.picked_crest_m, "-o", color=c, ms=3, lw=1.2)
            axs[i, 0].plot(d.survey, d.independent_crest_m, "--", color=c, lw=1.0)
            axs[i, 1].plot(d.survey, d.independent_crest_m - d.picked_crest_m, "-o", color=c, ms=3, lw=1.2,
                           label=f"GIS {g}")
            axs[i, 2].plot(d.survey, d.independent_crest_seaward_of_road_m, "-o", color=c, ms=3, lw=1.2)
        axs[i, 0].set_ylabel(f"{label}\nCrest (m NAVD88)", fontsize=8)
        axs[i, 1].set_ylabel("Independent minus\npicked crest (m)", fontsize=8)
        axs[i, 2].set_ylabel("Independent crest,\nseaward of NC-12 (m)", fontsize=8)
        axs[i, 1].axhline(0, color="0.55", lw=0.6)
        axs[i, 1].legend(fontsize=6, ncol=2, frameon=False, loc="upper left")
        for ax in axs[i]:
            ax.grid(axis="y", color="0.92", lw=0.5)
            ax.spines[["top", "right"]].set_visible(False)
            for x in (1998.7, 1999.7, 2003.7):
                ax.axvline(x, color="0.75", lw=0.6, ls=":")
    for ax in axs[-1]:
        ax.set_xlabel("Survey")
        ax.set_xticks([1996, 1998, 2000, 2002, 2004, 2006, 2008])
    for k, ax in enumerate(axs.ravel()):
        ax.text(-0.04, 1.02, f"({chr(97 + k)})", transform=ax.transAxes, fontsize=9, fontweight="bold",
                ha="right", va="bottom")
    fig.legend(handles=[Line2D([], [], color="0.3", marker="o", ms=3, lw=1.2,
                               label="Picked crest (the survey crests the model is scored on)"),
                        Line2D([], [], color="0.3", lw=1.0, ls="--", label="Highest point between NC-12 and the beach")],
               loc="outside upper center", ncol=2, fontsize=7.5, frameon=False)
    out = OUT / "hotspot_crest_check.png"
    S_.save(fig, out, bbox_inches="tight", close=True)
    S_.record_caption(out, (
        "Are the lidar crest drops real? Per domain (colours) of Buxton (GIS 10-15), Avon south (GIS 21-26) and "
        "Rodanthe (GIS 82-89), every 1996-2009 survey on the 2009 datum: (left) the domain-median picked crest (each "
        "profile's highest point inside the start's dune window, shifted landward with the shore: the survey crests "
        "the model is scored against; solid) and, independently of any window, the domain-median highest point "
        f"between {SEAWARD_OF_ROAD_M:.0f} m seaward of NC-12 (2008 alignment) and the beach (the first point below "
        "mean high water; dashed); (middle) the independent minus the picked crest; (right) how far seaward of NC-12 "
        "the independent crest lies. Dotted verticals: Bonnie (1998), Dennis (1999), Isabel (2003). At Buxton the "
        "2008 road is the relocated one, so before 2000 the old road and its dune lie seaward of it. At Avon south "
        "(1996-2004) and GIS 82-83 the highest point is often a rooftop, not a dune. *2001 is unclassified lidar."))
    print(f"wrote {out}")


def main():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from site_layer import hat_figure_style as S_
    S_.apply_style()
    names = sys.argv[1:] or list(L.REGIONS)
    (OUT / "supporting").mkdir(parents=True, exist_ok=True)
    picked = pd.read_csv(CREST_TABLE)
    picked = picked[picked.status == "ok"]
    crest, datum = [], []
    for name in names:
        r = L.REGIONS[name]
        surv, g, dom = L.surveys(r)
        roads = L.road_lines(g[2])
        datum.append(L.datum_check(r, surv, g, roads))
        st = L.Strip(r, g, roads, dom)
        strips(r, surv, st, S_, plt)
        strip_changes(r, surv, st, S_, plt)
        profiles(r, surv, g, dom, roads, S_, plt)
        imagery_strips(r, st, g, S_, plt)
        plan_lidar(r, surv, g, dom, roads, S_, plt)
        plan_changes(r, surv, g, dom, roads, S_, plt)
        plan_aerial(r, g, dom, roads, S_, plt)
        crest += crest_rows(r, st, surv, picked)
        print(f"  {r.label} done")
    pd.concat(datum).round(3).to_csv(OUT / "supporting" / "datum_check.csv", index=False)
    if set(names) == set(L.REGIONS):
        t = pd.DataFrame(crest)
        t.round(3).to_csv(OUT / "supporting" / "hotspot_crest_check.csv", index=False)
        crest_check_figure(t, S_, plt)


if __name__ == "__main__":
    main()
