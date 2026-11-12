"""
Where along the island, and when, the model without the collision module misses the lidar dunes.

    python scripts/real-world-cases/Hatteras/dune-collision/1-observations/HAT_collision_motivation_alongshore.py

Five figures beside collision_motivation.png, each with a calibration (1996-2009) and a test (2009-2025) part,
both from the module-off controls:

    collision_offset_heatmap.png     model minus lidar crest per GIS domain and survey, under a strip map of the island
    collision_hotspot_crests.png     crest through time at each hotspot, model band against lidar dots
    collision_alongshore_crests.png  crest along the island at each period's end (2009, 2020)
    collision_alongshore_crests_<period>_all.png   the same in every survey, one panel each (supplement)
    collision_drivers_scatter.png    each domain's crest offset against the storms that put its dunes in
                                     collision and against its dry beach width

Model: cal_control_g16 (2026-10-26-dune-collision; Hs 3.0, StormStart 1, groin 0.16, module off) and the
2009-2025 test's control (2026-10-29-collision-test-2009-2025, module off, the 2014, 2017 and 2022 fills in).
Lidar: survey_dune_crests.csv of each start (HAT_survey_dune_crests.py), read and set against the model state
its model_state names as cascade_pipeline/plotting/dune_crest_comparison.py does (a cell whose dune window is
now beach counts as the berm). The 2009-start surveys are 2011 (post-Irene), 2012 (post-Sandy), 2014,
2016-2019 and 2020 (Feb; the last public lidar of the island, searched 2026-11-03).

Collision hotspots (since 2026-11-03; collision_hotspots), found in each period on its own: stretches of at least
two neighbouring domains (one-domain gaps bridged) whose modelled crest stands more than HOT_OFF_M above the
lidar's in at least HOT_PERSIST of the surveys after the start, HOT_OFF_M on average, AND whose lidar crest fell at
least HOT_DROP_M from the period's first survey to its last while the model's fell at least HOT_DROP_M less: the
crest lowering the model did not reproduce, the loss the module is there to supply. Domains too high without that
missed loss (dune growth or rebuilding too fast, a start mismatch) are not hotspots; the per-domain table marks
them too_high_other. Before 2026-11-03: the 4 runs with the largest mean offset over either period.

Collision exposure: per cell and storm of the run's own storm series, the storm's peak runup (R_high, m NAVD88)
above the cell's dune toe (the module's toe, dune-toes/dune_toe_cells.csv) and below the model's crest at the
start of that storm's year (Sallenger's 2000 collision regime), counted over the run and averaged over the
domain's cells: the storms the model saw reach its dunes and, with the module off, did nothing about.

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-11-03
"""
from __future__ import annotations

import functools
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

_HERE = Path(__file__).resolve()
_HAT_DIR = next(_p for _p in _HERE.parents if (_p / "hatteras_paths.py").exists())
for _p in (_HAT_DIR, _HERE.parent):
    if str(_p) not in sys.path:
        sys.path.insert(0, str(_p))
from hatteras_paths import WORKSPACE as PROJECT_ROOT  # noqa: E402

import matplotlib  # noqa: E402
import matplotlib.patheffects  # noqa: E402
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.colors import BoundaryNorm, ListedColormap  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402
from matplotlib.patches import Patch  # noqa: E402
from matplotlib.ticker import FixedLocator, FormatStrFormatter, MultipleLocator  # noqa: E402

from cascade_pipeline.plotting import dune_crest_comparison as dc  # noqa: E402
from cascade_pipeline.saved_runs import load_saved_cascade  # noqa: E402
from site_layer import hat_topo_version as tv  # noqa: E402
from site_layer.hat_figure_style import (DOMAIN_AXIS_LABEL, INK, INK_MUTED, apply_style, record_caption,  # noqa: E402
                                         save, town_bands)
from site_layer.hatteras_site_config import HATTERAS_ANNOTATIONS, HATTERAS_DOMAINS  # noqa: E402
import HAT_collision_motivation as M  # noqa: E402

# --- CONFIG ------------------------------------------------------------------
# The topic folders under 1-observations/ (3-mechanism is HAT_collision_motivation.py's)
WHERE_DIR = M.ROOT_DIR / "1-mismatch-along-island"     # the heatmap and the crests along the island
HOT_DIR = M.ROOT_DIR / "2-hotspots"
DRIVER_DIR = M.ROOT_DIR / "4-storms-and-beach-width"    # the offset against collision storms and beach width
RUNS = M.RUNS
PERIODS = {
    "calibration": dict(label="Calibration, 1996–2009", start="1996-start",
                        run=RUNS / "2026-10-26-dune-collision" / "runs" / "cal_control_g16" / "1996_2009" / "edgeBE"
                        / "HAT_1996-2009_managed_groin_edgeBE"),
    "test": dict(label="Test, 2009–2025", start="2009-start",
                 run=RUNS / "2026-10-29-collision-test-2009-2025" / "runs" / "control" / "2009_2025" / "edgeBE"
                 / "HAT_2009-2025_managed_groin_edgeBE_nocollision"),
}
MIN_CELLS = 10                   # a domain-survey counts with this many surveyed cells
# Collision hotspots (collision_hotspots, 2026-11-03): persistently too high AND observed crest lowering missed
HOT_OFF_M, HOT_PERSIST, HOT_DROP_M = 1.0, 0.5, 1.0
# The survey per period drawn along the island: each period's last (2026-11-03; was three per period, too busy)
END_SURVEY = {"calibration": 2009, "test": 2020}
UNCLASSIFIED = {2001}
SURVEY_NOTES = {1999: "post-Dennis", 2011: "post-Irene", 2012: "post-Sandy", 2004: "post-Isabel"}
# Offset bands as the run overview's (dune_crest_comparison), RdBu with a neutral middle
OFFSET_BINS_M = dc.OFFSET_BINS_M
NO_SURVEY = dc.NO_SURVEY_COLOR
MODEL_C, LIDAR_C, HIGH_C, LOW_C = M.MODEL_C, M.LIDAR_C, "#B2182B", "#2166AC"
EXPOSURE_C = M.WATER_C
BEACH_C = "#B08D57"              # sand
HOT_SHADE = "#FBE3D6"
# The heatmap's island strip, a line map: white island, dark coast, light-grey domain lines (every tenth darker)
STRIP_LAND, STRIP_EDGE, STRIP_GRID, STRIP_GRID_10 = "white", "0.25", "0.72", "0.45"
STRIP_GAP = 0.035                # figure fraction between the map and panel (a)'s axes: room for its title
STRIP_HALO = [matplotlib.patheffects.withStroke(linewidth=2.2, foreground="white")]   # labels over the coast
SERIES_TOWN_C = "#F7F7F7"                                  # the small multiples' village bands
SHOAL_C = HATTERAS_ANNOTATIONS.color_shoal                 # the house shoal amber
SHOAL_INK = "#8A6A00"                                       # its name, darker to read as text
STRIP_CUT = "dashed"               # how the land cut west of Buxton reads (island_strip): open, fade, break, dashed
STRIP_FADE, STRIP_BEYOND_M = 2.5, 1500.0   # fade length (domains); dashed coast drawn this far past the cut (m)
STRIP_SOUTH_DOMAINS = 6.0        # domains' worth of island drawn south of GIS 1: Cape Hatteras
UNDER_C = "#D9D9D9"             # the alongshore figure's "model too low" fill: neutral, so the red reads
OFFSET_CLIP_M = 4.0              # the driver figures' offset axis; a bar beyond it is cut and labelled
HOT_DEF = (f"the stretches of at least two neighbouring GIS domains (a one-domain gap bridged), found in each "
           f"period on its own, where the modelled crest stands more than {HOT_OFF_M:g} m above the lidar's in at "
           f"least half the surveys after the period's start ({HOT_OFF_M:g} m on average) and the lidar crest fell "
           f"at least {HOT_DROP_M:g} m from the period's first survey to its last while the model's fell at least "
           f"{HOT_DROP_M:g} m less: observed crest lowering the model did not reproduce")
# -----------------------------------------------------------------------------


def load(key):
    p = PERIODS[key]
    run = p["run"]
    meta = json.loads((run / f"{run.name}_run_metadata.json").read_text(encoding="utf-8"))
    c = load_saved_cascade(run / f"{run.name}.npz")
    crests = dc.model_crests(c)                                   # (states, 90, cells), m NAVD88
    dom = HATTERAS_DOMAINS.for_total(len(c.barrier3d))
    b0 = c.barrier3d[dom.gis_to_pad(1)]
    mhw = b0._MHW * 10
    berm = (b0._BermEl + b0._MHW) * 10
    fmy = int(meta["period"]["first_model_year"])
    start_year = int(meta["period"]["start_year"])
    surveys = dc.survey_table(start_year, berm)
    table = dc.compare(crests, surveys, fmy)
    table = table[table.cells >= MIN_CELLS].copy()
    table["when"] = table.survey.map(surveys.groupby("survey").when.first())
    storms = dc.big_storms(PROJECT_ROOT / "data" / "hatteras_init" / meta["identity"]["storm_file"])
    toes = pd.read_csv(tv.product_dir(p["start"]) / "dune-toes" / "dune_toe_cells.csv").set_index(["gis", "cell"])
    ss = np.asarray(b0._StormSeries, float)
    return dict(key=key, label=p["label"], crests=crests, years=fmy + np.arange(crests.shape[0]), fmy=fmy,
                start_year=start_year, surveys=surveys, table=table, storms=storms, toes=toes, mhw=mhw,
                storm_series=ss)


def mean_offset(d):
    """Per GIS, the crest offset averaged over the surveys after the start."""
    t = d["table"]
    return t[t.survey > d["start_year"]].groupby("gis").model_minus_survey_m.mean().reindex(range(1, 91))


def span(h):
    return f"GIS {h['lo']}" if h["lo"] == h["hi"] else f"GIS {h['lo']}–{h['hi']}"


def hotspots(data):
    """{period: [hotspot]} from collision_hotspots, each hotspot a dict (period, lo, hi and its stretch means), and
    the per-domain table of both periods."""
    sets, tables = {}, []
    for key, d in data.items():
        m, runs, _ = collision_hotspots(d)
        tables.append(m.assign(period=key).reset_index().rename(columns={"index": "gis"}))
        sets[key] = [dict(period=key, lo=lo, hi=hi, **m.loc[lo:hi, ["mean_m", "lidar_drop_m", "missed_loss_m"]]
                          .mean().round(3).to_dict()) for lo, hi in runs]
    return sets, pd.concat(tables)


def stretches(flags, min_len=2):
    """Runs of True along GIS 1-90 (one-domain gaps bridged), at least min_len domains: [(lo, hi)]."""
    f = np.asarray(flags, bool).copy()
    for i in range(1, len(f) - 1):
        if not f[i] and f[i - 1] and f[i + 1]:
            f[i] = True
    out, i = [], 0
    while i < len(f):
        if f[i]:
            j = i
            while j + 1 < len(f) and f[j + 1]:
                j += 1
            if j - i + 1 >= min_len:
                out.append((i + 1, j + 1))
            i = j + 1
        else:
            i += 1
    return out


def collision_hotspots(d):
    """Per period (proposed 2026-11-03): a domain is persistently too high when its modelled crest stands more
    than HOT_OFF_M above the lidar's in at least HOT_PERSIST of the surveys after the start, HOT_OFF_M on
    average; it is a collision hotspot when, besides, the lidar crest fell at least HOT_DROP_M from the period's
    first survey to its last and the model's fell at least HOT_DROP_M less (the loss the model missed).
    Returns (per-domain table, collision stretches, other too-high stretches)."""
    t = d["table"]
    post = t[t.survey > d["start_year"]]
    g = post.groupby("gis").model_minus_survey_m
    m = pd.DataFrame({"persist": g.apply(lambda x: (x > HOT_OFF_M).mean()), "mean_m": g.mean()})
    by = t.sort_values("survey").groupby("gis")
    first, last = by.first(), by.last()
    m["lidar_drop_m"] = first.survey_median_m - last.survey_median_m
    m["model_drop_m"] = first.model_median_m - last.model_median_m
    m["missed_loss_m"] = m.lidar_drop_m - m.model_drop_m
    m = m.reindex(range(1, 91))
    high = (m.persist >= HOT_PERSIST) & (m.mean_m >= HOT_OFF_M)
    col = high & (m.lidar_drop_m >= HOT_DROP_M) & (m.missed_loss_m >= HOT_DROP_M)
    m["class"] = np.where(col, "collision", np.where(high, "too_high_other", ""))
    col_runs = stretches(col.to_numpy())
    covered = np.zeros(90, bool)
    for lo, hi in col_runs:
        covered[lo - 1:hi] = True
    other_runs = stretches((high & ~col).to_numpy() & ~covered, min_len=1)
    return m, col_runs, other_runs


# The looser hotspot rule (2026-11-12): the drop to any later survey, persistence counted from it, 0.75 m
LOOSE_OFF_M, LOOSE_PERSIST, LOOSE_DROP_M = 0.75, 0.5, 0.75


def collision_hotspots_loose(d):
    """As collision_hotspots, loosened for drops the strict rule cannot see. Per domain, against the period's first
    survey, the survey after it where the missed loss (lidar crest loss less the model's) is largest among those
    whose lidar loss is at least LOOSE_DROP_M: a drop the dune later recovered from (rebuilt) still counts. The
    overprediction is then judged from that survey on (more than LOOSE_OFF_M in at least LOOSE_PERSIST of them,
    LOOSE_OFF_M on average), so a late-starting mismatch is not diluted by the years before. Returns
    (per-domain table, stretches of two or more domains, single domains meeting it)."""
    t = d["table"].sort_values("survey")
    rows = {}
    for gis, s in t.groupby("gis"):
        s = s.reset_index(drop=True)
        if len(s) < 2:
            continue
        lidar_drop = s.survey_median_m.iloc[0] - s.survey_median_m
        missed = lidar_drop - (s.model_median_m.iloc[0] - s.model_median_m)
        ok = (lidar_drop >= LOOSE_DROP_M) & (s.index > 0)
        k = int(missed.where(ok).idxmax()) if ok.any() else int(missed.iloc[1:].idxmax())
        after = s.model_minus_survey_m.iloc[k:]
        rows[gis] = dict(drop_survey=int(s.survey.iloc[k]), lidar_drop_m=lidar_drop.iloc[k],
                         missed_loss_m=missed.iloc[k], persist=(after > LOOSE_OFF_M).mean(), mean_m=after.mean())
    m = pd.DataFrame.from_dict(rows, orient="index").reindex(range(1, 91))
    hit = ((m.persist >= LOOSE_PERSIST) & (m.mean_m >= LOOSE_OFF_M) & (m.lidar_drop_m >= LOOSE_DROP_M)
           & (m.missed_loss_m >= LOOSE_DROP_M)).fillna(False).to_numpy()
    runs = stretches(hit)
    covered = np.zeros(90, bool)
    for lo, hi in runs:
        covered[lo - 1:hi] = True
    m["hit"] = hit
    return m, runs, [i + 1 for i in np.flatnonzero(hit & ~covered)]


def exposure(d):
    """Per GIS: storms per cell in collision (toe < R_high < model crest at the start of the storm's year)."""
    ss, crests, toes = d["storm_series"], d["crests"], d["toes"]
    rows = []
    for g in range(1, 91):
        toe = np.array([toes.loc[(g, c), "toe_m_navd88"] if (g, c) in toes.index else np.nan
                        for c in range(crests.shape[2])])
        n_col = np.zeros(crests.shape[2])
        n_over = np.zeros(crests.shape[2])
        for t, r_high in zip(ss[:, 0].astype(int), ss[:, 1]):
            if not 1 <= t <= crests.shape[0] - 1:
                continue
            r = r_high * 10 + d["mhw"]
            crest = crests[t - 1, g - 1]
            n_col += (toe < r) & (r < crest)
            n_over += r >= crest
        ok = np.isfinite(toe)
        rows.append(dict(gis=g, collision_storms_per_cell=float(n_col[ok].mean()) if ok.any() else np.nan,
                         overtopping_storms_per_cell=float(n_over.mean())))
    return pd.DataFrame(rows).set_index("gis")


def beach_width(d):
    return d["toes"].groupby("gis").beach_width_m.median().reindex(range(1, 91))


def gis_axis(ax, label=True):
    ax.set_xlim(0.5, 90.5)
    ax.xaxis.set_major_locator(FixedLocator([1, 10, 20, 30, 40, 50, 60, 70, 80, 90]))
    ax.xaxis.set_minor_locator(MultipleLocator(5))
    if label:
        ax.set_xlabel(DOMAIN_AXIS_LABEL)


def shade_hotspots(ax, hots, label=False):
    for h in hots:
        ax.axvspan(h["lo"] - 0.5, h["hi"] + 0.5, color=HOT_SHADE, lw=0, zorder=0)
        if label:
            ax.text((h["lo"] + h["hi"]) / 2, 1.0, span(h).replace("GIS ", ""), transform=ax.get_xaxis_transform(),
                    ha="center", va="bottom", fontsize=6.5, color=INK)


@functools.lru_cache(maxsize=1)
def strip_layers():
    """The island and its domains in the strip's frame: turned a quarter so north points right (the ocean below,
    the sound above), in units of one domain (the boxes stand ~502.5 m apart due north in UTM), so each domain
    sits over its own column of the plot below. Returns (land rings, coast lines,
    {gis: [rings of the box inside the land]}, the cut's y, the coast just beyond the cut)."""
    import geopandas as gpd
    from shapely.geometry import box
    from site_layer import hat_map_layers as ml
    from site_layer.hat_observed_rates import DOMAIN_BOXES
    dom = gpd.read_file(DOMAIN_BOXES).to_crs(26918)
    dom["ID"] = dom.domain_id.astype(int)
    dom = dom[(dom.ID >= 1) & (dom.ID <= 90)].sort_values("ID")
    step, n0 = np.polyfit(dom.ID.to_numpy(), dom.geometry.centroid.y.to_numpy(), 1)   # northing = n0 + step*gis
    e_mid = float(dom.geometry.centroid.x.mean())

    def to_ax(xy):                                                # (easting, northing) -> (gis, -east in domains)
        xy = np.asarray(xy)[:, :2]
        return np.c_[(xy[:, 1] - n0) / step, -(xy[:, 0] - e_mid) / step]

    def rings(geom):
        return [to_ax(g.exterior.coords) for g in getattr(geom, "geoms", [geom])
                if g.geom_type == "Polygon" and not g.is_empty]

    def lines(geom):
        parts = getattr(geom, "geoms", [geom])
        return [to_ax(g.coords) for g in parts if g.geom_type == "LineString" and not g.is_empty]

    x0, _, x1, _ = dom.total_bounds
    # the reach and, south of GIS 1, Cape Hatteras (its point ~2 km south): westward only to a box's width
    # beyond the domains, where the land is cut; the coast is drawn, the cut is not
    full = gpd.read_file(ml.ISLAND_OUTLINE).to_crs(26918).geometry.union_all()
    win = box(x0 - 2000, n0 - STRIP_SOUTH_DOMAINS * step, x1 + 2000, n0 + 90.5 * step)
    land = full.intersection(win)
    coast = full.boundary.intersection(win.buffer(-1.0))
    beyond = full.boundary.intersection(box(x0 - 2000 - STRIP_BEYOND_M, n0 - STRIP_SOUTH_DOMAINS * step,
                                            x0 - 2000, n0 + 90.5 * step))
    y_cut = -(x0 - 2000 - e_mid) / step
    return (rings(land), lines(coast), {int(i): rings(g.intersection(land)) for i, g in zip(dom.ID, dom.geometry)},
            y_cut, lines(beyond))


def strip_ylim():
    land = strip_layers()[0]
    ys = np.concatenate([r[:, 1] for r in land])
    return ys.min() - 4.6, ys.max() + 3.5         # room below for the shoal bands and their names


def island_strip(ax, cut=None):
    """Hatteras Island above a GIS-domain axis (strip_layers), as a line map: the coast, the domains cut to it
    with every tenth numbered and its line darker, and the villages. No water drawn. `cut` (STRIP_CUT by default):
    how the land cut off west of Buxton reads, "open" (the coast stops), "fade" (the coast fades out to the cut),
    "break" (a wavy break line across the land at the cut) or "dashed" (the coast continues dashed past it)."""
    from matplotlib.collections import LineCollection
    cut = cut or STRIP_CUT
    land, coast, doms, y_cut, beyond = strip_layers()
    # the cape lies left of GIS 1, past the axis: drawn into the margin (clip_on=False), the axis unchanged
    for r in land:
        ax.fill(*r.T, facecolor=STRIP_LAND, edgecolor="none", zorder=1, clip_on=False)
    for gis, rs in doms.items():
        tenth = gis == 1 or gis % 10 == 0
        for r in rs:
            ax.plot(*r.T, color=STRIP_GRID_10 if tenth else STRIP_GRID, lw=0.45 if tenth else 0.3,
                    zorder=2.5 if tenth else 2)
    for c in coast:
        if cut == "fade":
            seg = np.stack([c[:-1], c[1:]], axis=1)
            a = np.clip((y_cut - seg[:, :, 1].mean(axis=1)) / STRIP_FADE, 0, 1)
            rgba = np.tile(matplotlib.colors.to_rgba(STRIP_EDGE), (len(seg), 1))
            rgba[:, 3] = a
            ax.add_collection(LineCollection(seg, colors=rgba, lw=0.7, zorder=3, clip_on=False))
        else:
            ax.plot(*c.T, color=STRIP_EDGE, lw=0.7, zorder=3, clip_on=False)
    if cut == "fade":                                             # the land fades too, to the cut
        from matplotlib.colors import LinearSegmentedColormap
        x_lo = min(r[:, 0].min() for r in land)
        ax.imshow(np.linspace(0, 1, 50)[:, None], extent=(x_lo - 1, 12, y_cut - STRIP_FADE, y_cut + 0.05),
                  origin="lower", aspect="auto", zorder=2.8, clip_on=False,
                  cmap=LinearSegmentedColormap.from_list("f", [(1, 1, 1, 0), (1, 1, 1, 1)]))
    elif cut == "break":                                          # pair the coast ends at the cut, wave between
        ends = sorted(pt[0] for c in coast for pt in (c[0], c[-1]) if abs(pt[1] - y_cut) < 0.2)
        for xa, xb in zip(ends[0::2], ends[1::2]):
            t = np.linspace(0, 1, 60)
            ax.plot(xa + (xb - xa) * t, y_cut - 0.35 + 0.28 * np.sin(2 * np.pi * 2.5 * t), color=STRIP_EDGE,
                    lw=0.7, zorder=3, clip_on=False)
    elif cut == "dashed":
        for c in beyond:
            ax.plot(*c.T, color=STRIP_EDGE, lw=0.6, ls=(0, (2.5, 1.5)), zorder=3, clip_on=False)
    tip = min(land, key=lambda r: r[:, 0].min())
    tip = tip[tip[:, 0].argmin()]
    ax.text(tip[0] - 0.5, tip[1], "Cape\nHatteras", ha="right", va="center", fontsize=6.5, color=INK_MUTED,
            clip_on=False, path_effects=STRIP_HALO)
    for gis, rs in doms.items():
        if rs and (gis == 1 or gis % 10 == 0):
            ax.text(gis, min(r[:, 1].min() for r in rs) - 0.4, str(gis), ha="center", va="top", fontsize=6,
                    color=INK_MUTED, zorder=4, path_effects=STRIP_HALO)
    # the offshore shoals (the site's shoal zones, footprints picked by eye), hatched on the ocean side
    for name, (lo, hi) in HATTERAS_ANNOTATIONS.shoal_zones.items():
        y = min(r[:, 1].min() for g in range(lo, hi + 1) for r in doms[g]) - 1.9
        ax.fill_between([lo - 0.5, hi + 0.5], y - 0.55, y, facecolor=(*matplotlib.colors.to_rgb(SHOAL_C), 0.18),
                        edgecolor=SHOAL_C, hatch="////", lw=0, zorder=2)
        ax.text((lo + hi) / 2, y - 0.75, name, ha="center", va="top", fontsize=6.5, color=SHOAL_INK,
                style="italic", zorder=4, path_effects=STRIP_HALO)
    for name, (lo, hi) in HATTERAS_ANNOTATIONS.town_spans.items():
        top = max(r[:, 1].max() for g in range(lo, hi + 1) for r in doms[g])
        ax.plot([lo - 0.5, hi + 0.5], [top + 0.8] * 2, color=INK_MUTED, lw=1.6, solid_capstyle="butt", zorder=4)
        ax.text((lo + hi) / 2, top + 1.3, name, ha="center", va="bottom", fontsize=6.5, color=INK_MUTED, zorder=4,
                path_effects=STRIP_HALO)
    ax.annotate("N", xy=(0.99, 0.85), xytext=(0.945, 0.85), xycoords="axes fraction", ha="right", va="center",
                fontsize=7, color=INK_MUTED, arrowprops=dict(arrowstyle="-|>", color=INK_MUTED, lw=0.8))
    # the water, named where it is open: the ocean below the south end, the sound above the middle
    ax.text(0.10, 0.04, "Atlantic Ocean", transform=ax.transAxes, ha="left", va="bottom", fontsize=7,
            color=INK_MUTED, style="italic")
    ax.text(0.58, 0.80, "Pamlico Sound", transform=ax.transAxes, ha="center", va="center", fontsize=7,
            color=INK_MUTED, style="italic")
    ax.set_ylim(*strip_ylim())
    ax.set_facecolor("none")
    ax.tick_params(axis="both", which="both", left=False, labelleft=False, bottom=False, labelbottom=False)
    for sp in ax.spines.values():
        sp.set_visible(False)


def title(ax, letter, text):
    ax.set_title(f"({letter})", loc="left", fontweight="bold")
    ax.set_title(text, loc="center")


# ---- 1: heatmap ----------------------------------------------------------------------------------------------
def fig_heatmap(data, hots, mark_hotspots=True, out_png=None):
    """The offset heatmap. mark_hotspots outlines the collision hotspots (collision_hotspots) on it, over the
    surveys after the start; out_png set writes only that image (a preview: no table, no caption)."""
    colors = plt.get_cmap("RdBu_r", len(OFFSET_BINS_M) + 1)(np.arange(len(OFFSET_BINS_M) + 1))
    colors[len(colors) // 2] = (0.97, 0.97, 0.97, 1)
    cmap = ListedColormap(colors[1:-1])
    cmap.set_under(colors[0])
    cmap.set_over(colors[-1])
    cmap.set_bad(NO_SURVEY)
    norm = BoundaryNorm(OFFSET_BINS_M, cmap.N)
    # Time runs up the whole figure (Hannah, 2026-11-03): the test period on top, calibration below, and the
    # oldest survey at the bottom of each panel, so 1996 is at the foot and 2020 at the head
    periods = [data["test"], data["calibration"]]
    n_rows = [d["table"].survey.nunique() for d in periods]
    # The island on top, the domains over their columns (2026-11-03)
    # the strip's height in survey rows, so the island is drawn near true shape (1 domain ~ 0.069 in wide; a row ~0.3 in)
    y0, y1 = strip_ylim()
    fig, (ax_map, *axs) = plt.subplots(3, 1, figsize=(7.48, 7.1 + (y1 - y0) * 0.069), sharex=True,
                                       gridspec_kw=dict(height_ratios=[(y1 - y0) * 0.069 / 0.3] + n_rows,
                                                        hspace=0.32))
    island_strip(ax_map)
    out = []
    for ax, (letter, d) in zip(axs, zip("ab", periods)):
        piv = d["table"].pivot(index="survey", columns="gis", values="model_minus_survey_m").reindex(
            columns=range(1, 91))
        piv = piv.sort_index()
        im = ax.pcolormesh(np.arange(0.5, 91), np.arange(len(piv) + 1), np.ma.masked_invalid(piv.to_numpy()),
                           cmap=cmap, norm=norm, edgecolors="white", linewidth=0.3)
        ax.set_ylim(0, len(piv))                                   # time runs up
        ax.set_yticks(np.arange(len(piv)) + 0.5)
        ax.set_yticklabels([f"{s}{'*' if s in UNCLASSIFIED else ''}"
                            + (f" {SURVEY_NOTES[s]}" if s in SURVEY_NOTES else "") for s in piv.index], fontsize=7)
        gis_axis(ax, label=letter == "b")
        ax.set_ylabel("Lidar survey")
        ax.set_title(f"({letter})", loc="left", fontweight="bold")
        ax.set_title(d["label"], loc="center")
        out.append(piv.stack(dropna=False).rename("model_minus_lidar_m").reset_index().assign(period=d["key"]))
        if mark_hotspots:
            from matplotlib.patches import Rectangle
            _, col_runs, _ = collision_hotspots(d)                   # the other too-high domains: not drawn
            y_lo = int((piv.index <= d["start_year"]).sum())         # the scored surveys: after the start
            for lo, hi in col_runs:
                ax.add_patch(Rectangle((lo - 0.5, y_lo), hi - lo + 1, len(piv) - y_lo, facecolor="none",
                                       edgecolor=INK, lw=1.6, zorder=5, clip_on=False))
            for lo, hi in col_runs:
                ax.text((lo + hi) / 2, len(piv) - 0.08, f"{lo}–{hi}", ha="center", va="top", fontsize=6.5,
                        fontweight="bold", color=INK, zorder=6, path_effects=STRIP_HALO)
    # with the hotspots marked, the colour bar takes the left part of the row and their key the right
    cb = fig.colorbar(im, ax=axs, orientation="horizontal", fraction=0.035, pad=0.13, extend="both",
                      ticks=OFFSET_BINS_M, **(dict(shrink=0.52, anchor=(0.04, 1.0)) if mark_hotspots else {}))
    cb.set_label("Modelled minus observed dune crest (m)\n(domain median; positive, overpredicted)" if mark_hotspots
                 else "Model minus lidar crest, domain median (m): red, model dunes too high")
    cb.ax.tick_params(labelsize=7)
    # the map sits just over (a)'s title, not a full panel gap above it (the crop takes the space left on top)
    m, a = ax_map.get_position(), axs[0].get_position()
    ax_map.set_position([m.x0, a.y1 + STRIP_GAP, m.width, m.height])
    if mark_hotspots:                                              # what the outline means, beside the colour bar
        cbp = cb.ax.get_position()
        fig.legend(handles=[
            Patch(facecolor="none", edgecolor=INK, lw=1.6,
                  label="Collision hotspot:\nobserved crest lowering not reproduced")],   # thresholds: caption
            loc="center left", bbox_to_anchor=(cbp.x1 + 0.06, cbp.y0 + cbp.height / 2), ncol=1, frameon=False,
            fontsize=7.5, handlelength=2.2)
    if out_png is not None:                                        # a preview: the image only
        save(fig, out_png, dpi=300, facecolor="white", bbox_inches="tight", close=True)
        print(f"wrote {out_png}")
        return
    out_png = WHERE_DIR / "collision_offset_heatmap.png"
    save(fig, out_png, dpi=300, facecolor="white", bbox_inches="tight", close=True)
    pd.concat(out).round(3).to_csv(WHERE_DIR / "supporting" / "collision_offset_heatmap.csv", index=False)
    record_caption(out_png, (
        "Where and when the modelled dune crests leave the lidar's, without the collision module. Each cell: the "
        "median modelled crest minus the median lidar crest over a GIS domain's surveyed 10 m cells (at least "
        f"{MIN_CELLS}), against the model state of the 1 January the survey is compared with; red, the model's "
        "dunes too high; white, within 0.5 m; grey, fewer surveyed cells. Boxed (over the surveys after each "
        f"period's start), the collision hotspots: {HOT_DEF}. Time runs up the figure, oldest survey "
        "at the bottom. (a) Test, 2009–2025, the test control (the 2014, 2017 and 2022 fills in); (b) calibration, "
        "1996–2009, the calibration control (Hs 3.0, StormStart 1, groin 0.16). Both with the module off. A cell whose dune window is now beach counts as "
        "the berm. *2001, unclassified lidar. 2011 and 2012 were gridded from the ground returns at 1 m; 2020 "
        "(flown Feb 2020 at Buxton) is the last public lidar of the island. Top: Hatteras Island turned a quarter so "
        "north points right (the ocean side below), true to scale, from Cape Hatteras (left of GIS 1), divided into its 90 GIS domains, each over its "
        "column, every tenth numbered; the coast dashed where the island continues west beyond the map; hatched "
        "offshore, the alongshore footprints of Avon and Wimble Shoals (picked by eye); bars, Buxton, "
        "Avon and the Tri-Village (Salvo, Waves, Rodanthe)."))
    print(f"wrote {out_png}")


# ---- 2: hotspot time series -------------------------------------------------------------------------------------
def fig_hotspots(data, hots):
    """Crest through time at each collision hotspot, in the period it belongs to: calibration hotspots in the
    left column, test hotspots in the right."""
    keys = list(data)                                              # calibration, test
    n = max(len(hots[k]) for k in keys)
    fig, axs = plt.subplots(n, 2, figsize=(7.48, 1.95 * n + 0.6), sharey=True,
                            gridspec_kw=dict(width_ratios=[13, 16]), squeeze=False)
    rows, letter = [], iter("abcdefghij")
    for i in range(n):
        for j, key in enumerate(keys):
            ax, d = axs[i, j], data[key]
            if i >= len(hots[key]):
                ax.set_visible(False)
                continue
            h = hots[key][i]
            gs = range(h["lo"], h["hi"] + 1)
            mod = np.concatenate([d["crests"][:, g - 1, :] for g in gs], axis=1)
            q = np.percentile(mod, [25, 50, 75], axis=1)
            ax.fill_between(d["years"], q[0], q[2], color=MODEL_C, alpha=0.22, lw=0)
            ax.plot(d["years"], q[1], color=MODEL_C, lw=1.8)
            s = d["surveys"][d["surveys"].gis.isin(gs)]
            for sv, grp in s.groupby("survey"):
                qq = grp.crest.quantile([0.25, 0.5, 0.75]).to_numpy()
                ax.errorbar(grp.when.iloc[0], qq[1], yerr=[[qq[1] - qq[0]], [qq[2] - qq[1]]], fmt="o", ms=4,
                            color=LIDAR_C, mfc="white" if sv in UNCLASSIFIED else LIDAR_C, elinewidth=1,
                            capsize=0, zorder=5)
                rows.append(dict(hotspot=span(h), gis=f"{h['lo']}-{h['hi']}", period=key, survey=sv,
                                 lidar_p25=qq[0], lidar_median=qq[1], lidar_p75=qq[2]))
            for t, lab in d["storms"]:
                ax.axvline(t, color=INK_MUTED, lw=0.6, ls=":")
                if i == 0:
                    ax.text(t - 0.15, 0.97, lab, rotation=90, ha="right", va="top", fontsize=6, color=INK_MUTED,
                            transform=ax.get_xaxis_transform())
            lo_y = d["years"][0] - 1
            ax.set_xlim(lo_y - 0.3, d["years"][-1] + 0.3)
            ax.xaxis.set_major_locator(FixedLocator(range(int(np.ceil(lo_y / 2) * 2), int(d["years"][-1]) + 1, 2)))
            ax.xaxis.set_minor_locator(MultipleLocator(1))
            ax.xaxis.set_major_formatter(FormatStrFormatter("%d"))
            # heading: the hotspot on the left, the crest lowering the model missed on the right
            ax.set_title(f"({next(letter)}) {span(h)}", loc="left", fontweight="bold", fontsize=8.5)
            ax.set_title(f"missed lowering {h['missed_loss_m']:.1f} m", loc="right", fontsize=7, color=INK_MUTED)
            if j == 0 or not axs[i, 0].get_visible():
                ax.set_ylabel("Crest (m NAVD88)")
    fig.legend(handles=[Line2D([], [], color=MODEL_C, lw=2, label="Model without the module (median; band 25th–75th)"),
                        Line2D([], [], color=LIDAR_C, marker="o", ls="none", label="Lidar (median; bar 25th–75th)")],
               loc="lower center", ncol=2, frameon=False, fontsize=7.5)
    fig.tight_layout(rect=(0, 0.035, 1, 0.965), h_pad=0.8)
    # the period over each column, above the panel headings (a title would share their offset)
    for j, key in enumerate(keys):
        p = axs[0, j].get_position()
        fig.text((p.x0 + p.x1) / 2, 0.985, data[key]["label"], ha="center", va="top", fontsize=10.5, color=INK)
    out = HOT_DIR / "collision_hotspot_crests.png"
    save(fig, out, dpi=300, facecolor="white", close=True)
    pd.DataFrame(rows).round(3).to_csv(HOT_DIR / "supporting" / "collision_hotspot_crests.csv", index=False)
    record_caption(out, (
        f"Dune crests at the collision hotspots, {HOT_DEF}; each in the period it was found in, calibration (left) "
        "and test (right). Line and band: the model without the collision module, median and 25th–75th percentile "
        "over the hotspot's 10 m cells, each 1 January; dots and bars: each lidar survey on the same cells, at the "
        "middle of the months flown (hollow, 2001 unclassified lidar; a cell whose dune window is now beach counts "
        "as the berm). Right of each heading: the crest lowering the model missed, the lidar's crest loss from "
        "the period's first survey to its last less the model's, averaged over the hotspot's domains. Dotted: each "
        "period's three largest storms by R_high in the model's storm series. The test period has no lidar after "
        "Feb 2020."))
    print(f"wrote {out}")


# ---- 3: alongshore crests ---------------------------------------------------------------------------------------
def crest_panel(ax, d, sv, lw=1.4):
    """Model and lidar crests along the island for one survey, the gap filled red (model too high) or grey.
    Returns the per-domain table and (domains > 1 m too high, domains compared)."""
    s = d["surveys"][d["surveys"].survey == sv]
    state = int(s.state_year.iloc[0]) - d["fmy"]
    g_rows = []
    for g, grp in s.groupby("gis"):
        if len(grp) < MIN_CELLS:
            continue
        m = d["crests"][state, g - 1, grp.cell.to_numpy()]
        g_rows.append(dict(gis=g, model_median=np.median(m), lidar_median=grp.crest.median()))
    t = pd.DataFrame(g_rows, columns=["gis", "model_median", "lidar_median"]).set_index("gis").reindex(range(1, 91))
    x, mod, lid = t.index.to_numpy(), t.model_median.to_numpy(), t.lidar_median.to_numpy()
    ok = np.isfinite(mod) & np.isfinite(lid)
    ax.fill_between(x, lid, mod, where=ok & (mod > lid), step="mid", color=HIGH_C, alpha=0.35, lw=0)
    ax.fill_between(x, lid, mod, where=ok & (mod < lid), step="mid", color=UNDER_C, lw=0)
    ax.plot(x, mod, color=MODEL_C, lw=lw, drawstyle="steps-mid")
    ax.plot(x, lid, color=LIDAR_C, lw=lw, drawstyle="steps-mid")
    t.insert(0, "model_state", f"{int(s.state_year.iloc[0])}-01-01")
    return t, (int((mod - lid > 1)[ok].sum()), int(ok.sum()))


def crest_legend(fig, y=0.0):
    fig.legend(handles=[Line2D([], [], color=MODEL_C, lw=2, label="Model without the module"),
                        Line2D([], [], color=LIDAR_C, lw=2, label="Lidar"),
                        Patch(color=HIGH_C, alpha=0.35, lw=0, label="Model dunes too high"),
                        Patch(color=UNDER_C, lw=0, label="Model dunes too low")],
               loc="lower center", bbox_to_anchor=(0.5, y), ncol=4, frameon=False, fontsize=7.5)


def fig_alongshore(data, hots):
    """One panel per period, its last survey: model and lidar crests along the island, the gap filled."""
    fig, axs = plt.subplots(2, 1, figsize=(7.48, 5.0), sharex=True, sharey=True)
    rows, gaps = [], {}
    for ax, letter, d in zip(axs, "ab", data.values()):
        sv = END_SURVEY[d["key"]]
        t, gaps[d["key"]] = crest_panel(ax, d, sv)
        town_bands(ax, where="top", strip=0.07, label=True, fontsize=6.5)
        gis_axis(ax, label=letter == "b")
        ax.set_ylabel("Dune crest (m NAVD88)")
        ax.set_title(f"({letter})", loc="left", fontweight="bold")
        ax.set_title(f"{d['label']}: model on 1 Jan {t.model_state.iloc[0][:4]} vs. the {sv} lidar", loc="center")
        rows += [dict(period=d["key"], survey=sv, gis=g, **r) for g, r in t.dropna(how="all").iterrows()]
    axs[0].set_ylim(0, 11.5)
    crest_legend(fig)
    fig.tight_layout(rect=(0, 0.04, 1, 1), h_pad=0.8)
    out = WHERE_DIR / "collision_alongshore_crests.png"
    save(fig, out, dpi=300, facecolor="white", close=True)
    pd.DataFrame(rows).round(3).to_csv(WHERE_DIR / "supporting" / "collision_alongshore_crests.csv", index=False)
    cal, test = gaps["calibration"], gaps["test"]
    record_caption(out, (
        "Dune crests along the island at the end of each period, the model without the collision module against "
        "the lidar. (a) Calibration: the model on 1 January 2010 against the 2009 lidar; (b) test: the model on "
        "1 January 2020 against the 2020 lidar (the last public survey; flown Feb 2020 at Buxton). Each line is the "
        "median crest over a GIS domain's surveyed 10 m cells (domains with at least "
        f"{MIN_CELLS}); red fill, the model's dunes above the lidar's; grey, below. The model stands more than 1 m "
        f"too high in {cal[0]} of {cal[1]} domains in (a) and {test[0]} of {test[1]} in (b). Bars along the top: "
        "Buxton, Avon and the Tri-Village (Salvo, Waves, Rodanthe). The crests through time and in every survey "
        "are in collision_offset_heatmap.png."))
    print(f"wrote {out}; >1 m too high: calibration {cal}, test {test}")


# ---- 3s: alongshore crests in every survey (supplementary small multiples) ------------------------------------------
def fig_alongshore_series(data):
    """Per period, one thin panel per lidar survey, oldest at the top: collision_alongshore_crests_<period>_all.png."""
    for d in data.values():
        surveys = sorted(int(v) for v in d["table"].survey.unique())
        n = len(surveys)
        fig, axs = plt.subplots(n, 1, figsize=(7.48, 0.85 * n + 0.8), sharex=True, sharey=True)
        rows, counts = [], []
        for i, (ax, sv) in enumerate(zip(axs, surveys)):
            # villages: faint bands through every panel, named once above the top one
            for name, (lo, hi_) in HATTERAS_ANNOTATIONS.town_spans.items():
                ax.axvspan(lo - 0.5, hi_ + 0.5, color=SERIES_TOWN_C, lw=0, zorder=0)
                if i == 0:
                    ax.text((lo + hi_) / 2, 1.04, name, transform=ax.get_xaxis_transform(), ha="center",
                            va="bottom", fontsize=6.5, color=INK_MUTED)
            t, (hi, n_dom) = crest_panel(ax, d, sv, lw=1.0)
            counts.append((sv, hi, n_dom))
            # the survey in bold, the model date it is set against in grey after it
            note = f" {SURVEY_NOTES[sv]}" if sv in SURVEY_NOTES else ""
            star = "*" if sv in UNCLASSIFIED else ""
            lab = ax.text(0.006, 0.95, f"({'abcdefghijkl'[i]}) {sv}{star}{note}", transform=ax.transAxes,
                          ha="left", va="top", fontsize=7.5, fontweight="bold", color=INK, path_effects=STRIP_HALO)
            ax.annotate(f"   model, 1 Jan {t.model_state.iloc[0][:4]}", xy=(1, 0), xycoords=lab, ha="left",
                        va="bottom", fontsize=6.5, color=INK_MUTED, path_effects=STRIP_HALO)
            ax.text(0.994, 0.95, f"{hi} of {n_dom} domains >1 m too high", transform=ax.transAxes, ha="right",
                    va="top", fontsize=6.5, color=INK_MUTED, path_effects=STRIP_HALO)
            gis_axis(ax, label=i == n - 1)
            ax.set_yticks([2, 5, 8])
            # a light frame: no top or right border, tick marks on the bottom panel only
            for side in ("top", "right"):
                ax.spines[side].set_visible(False)
            for side in ("left", "bottom"):
                ax.spines[side].set_linewidth(0.6)
                ax.spines[side].set_color(INK_MUTED)
            ax.tick_params(axis="both", width=0.6, color=INK_MUTED, labelsize=7)
            if i < n - 1:
                ax.tick_params(axis="x", which="both", bottom=False)
            rows += [dict(period=d["key"], survey=sv, gis=g, **r) for g, r in t.dropna(how="all").iterrows()]
        axs[0].set_ylim(0.5, 12.6)                  # the crests reach ~10 m; the top ~2 m keeps the labels clear
        fig.supylabel("Dune crest (m NAVD88)", fontsize=9)
        crest_legend(fig)
        fig.tight_layout(rect=(0.01, 0.035, 1, 0.985), h_pad=0.1)
        name = f"collision_alongshore_crests_{d['key']}_all"
        out = WHERE_DIR / f"{name}.png"
        save(fig, out, dpi=300, facecolor="white", close=True)
        pd.DataFrame(rows).round(3).to_csv(WHERE_DIR / "supporting" / f"{name}.csv", index=False)
        first, last = counts[0], counts[-1]
        record_caption(out, (
            f"Dune crests along the island in every lidar survey of the {d['label'].lower()} period, the model "
            "without the collision module against the lidar, oldest survey at the top. Each panel: the median "
            f"crest over a GIS domain's surveyed 10 m cells (domains with at least {MIN_CELLS}) in the lidar (blue) "
            "and in the model (grey) at the 1 January the survey is compared with; red fill, the model's dunes "
            "above the lidar's; grey, below. Top right of each panel, the domains whose modelled crest stands more "
            f"than 1 m above the lidar's: {first[1]} of {first[2]} in {first[0]}, {last[1]} of {last[2]} in "
            f"{last[0]}. " + ("*2001, unclassified lidar. " if any(c[0] in UNCLASSIFIED for c in counts) else "")
            + "Shaded bands: Buxton, Avon and the Tri-Village (Salvo, Waves, Rodanthe). The supplement to "
            "collision_alongshore_crests.png; the same offsets as one picture in collision_offset_heatmap.png."))
        print(f"wrote {out}: " + ", ".join(f"{sv} {h}/{m}" for sv, h, m in counts))


# ---- 4 and 5: a driver above the offset --------------------------------------------------------------------------
def fig_drivers(data, hots):
    """One scatter per period and driver: each domain's mean crest offset against the storms that put its dunes in
    collision (top) and against its dry beach width (bottom); the collision hotspots in red."""
    keys = list(data)
    fig, axs = plt.subplots(2, len(keys), figsize=(7.48, 5.6), sharey=True)
    rows, stats = [], {}
    drivers = (("exposure", "Storms in collision per cell", EXPOSURE_C),
               ("beach", "Dry beach width (m)", BEACH_C))
    letter = iter("abcdefgh")
    for i, (kind, xlab, _) in enumerate(drivers):
        for j, key in enumerate(keys):
            d = data[key]
            ax = axs[i, j]
            off = mean_offset(d)
            x = exposure(d).collision_storms_per_cell if kind == "exposure" else beach_width(d)
            cls = collision_hotspots(d)[0]["class"].reindex(off.index).fillna("")
            in_run = pd.Series(False, index=off.index)
            for h in hots[key]:
                in_run.loc[h["lo"]:h["hi"]] = True
            cls = cls.where(~(cls.eq("collision") & ~in_run), "single")   # meets the rule alone, not a stretch
            ok = x.notna() & off.notna()
            y = off.clip(-OFFSET_CLIP_M, OFFSET_CLIP_M)
            for c, colour, size, z, lab in (("", "0.72", 14, 2, "Other domains"),
                                            ("too_high_other", "0.35", 16, 3, "Too high, no missed lowering"),
                                            ("single", "none", 22, 4, "Meets the hotspot rule alone"),
                                            ("collision", HIGH_C, 26, 5, "Collision hotspot")):
                sel = ok & cls.eq(c)
                ax.scatter(x[sel], y[sel], s=size, facecolor=colour, lw=0.9 if c == "single" else 0.4,
                           edgecolor=HIGH_C if c == "single" else "white", zorder=z,
                           label=lab if (i, j) == (0, 0) else None)
            for g in off[ok & (off.abs() > OFFSET_CLIP_M)].index:   # off the axis: drawn at the edge, value written
                ax.annotate(f"GIS {g} ({off[g]:.1f})", (x[g], y[g]), xytext=(4, -3 if off[g] > 0 else 3),
                            textcoords="offset points", fontsize=6, color=INK_MUTED,
                            va="top" if off[g] > 0 else "bottom")
            for h in hots[key]:                                        # each hotspot named once
                pts = [g for g in range(h["lo"], h["hi"] + 1) if ok.get(g, False)]
                if pts:
                    g = max(pts, key=lambda k: y[k])
                    ax.annotate(span(h).replace("GIS ", ""), (x[g], y[g]), xytext=(5, 2), textcoords="offset points",
                                fontsize=6.5, fontweight="bold", color=HIGH_C, path_effects=STRIP_HALO, zorder=5)
            ax.axhline(0, color=INK, lw=0.6, zorder=1)
            ax.axhline(HOT_OFF_M, color=INK_MUTED, lw=0.6, ls=(0, (3, 2)), zorder=1)
            if kind == "beach":
                ax.axvline(M.W0_M, color=INK, lw=0.6, ls="--", zorder=1)
                ax.text(M.W0_M, OFFSET_CLIP_M - 0.1, f" {M.W0_M:.0f} m", fontsize=6.5, color=INK, va="top")
            rho = pd.Series(x[ok]).rank().corr(pd.Series(off[ok]).rank())
            stats[(kind, key)] = (rho, int(ok.sum()))
            ax.text(0.98, 0.04, f"Spearman r = {rho:.2f}", transform=ax.transAxes, ha="right", fontsize=7,
                    color=INK)
            ax.set_ylim(-OFFSET_CLIP_M - 0.3, OFFSET_CLIP_M + 0.3)
            ax.set_xlabel(xlab)
            ax.set_title(f"({next(letter)})", loc="left", fontweight="bold")
            if i == 0:
                ax.set_title(d["label"], loc="center", pad=8)
            if j == 0:
                ax.set_ylabel("Model minus lidar crest (m)")
            rows.append(pd.DataFrame({"period": key, "driver": kind, "gis": off.index, "driver_value": x.to_numpy(),
                                      "mean_offset_m": off.to_numpy(),
                                      "class": cls.to_numpy()}))
    fig.legend(loc="lower center", ncol=4, frameon=False, fontsize=7)
    fig.tight_layout(rect=(0, 0.05, 1, 1), h_pad=1.0)
    out = DRIVER_DIR / "collision_drivers_scatter.png"
    save(fig, out, dpi=300, facecolor="white", close=True)
    pd.concat(rows).round(3).to_csv(DRIVER_DIR / "supporting" / "collision_drivers_scatter.csv", index=False)
    txt = "; ".join(f"{k} {p}: Spearman r {r:.2f} ({n} domains)" for (k, p), (r, n) in stats.items())
    print(f"collision_drivers_scatter: {txt}")
    record_caption(out, (
        "What sets where the model's dunes stand too high. Each point is one GIS domain: its model minus lidar "
        "crest, averaged over the surveys after the period's start (the domain median of the surveyed cells, at "
        f"least {MIN_CELLS}; beyond ±{OFFSET_CLIP_M:g} m drawn at the edge, value written), against (a, b) the "
        "storms that put its dunes in collision (Sallenger 2000: the storm's peak runup in the run's own storm series "
        "above the cell's dune toe and below the model's crest at the start of that storm's year; per cell, averaged "
        "over the domain), which Barrier3D lowers none of, and (c, d) its dry beach width (beach start, 0.5 m above "
        f"MHW, to the dune toe; domain median), dashed at the {M.W0_M:.0f} m beyond which the collision module's "
        f"beach factor is zero. Dotted, the hotspot rule's {HOT_OFF_M:g} m. Red, the collision hotspots: {HOT_DEF}; "
        "open red, a domain meeting that rule alone; dark grey, domains too high without that missed lowering. Model: the module-off controls. "
        f"Rank correlations: {txt}."))
    print(f"wrote {out}")


def main():
    apply_style()
    for d in (WHERE_DIR, HOT_DIR, DRIVER_DIR):
        (d / "supporting").mkdir(parents=True, exist_ok=True)
    data = {k: load(k) for k in PERIODS}
    hots, by_domain = hotspots(data)
    print("collision hotspots:", {k: [(h["lo"], h["hi"]) for h in v] for k, v in hots.items()})
    pd.DataFrame([h for v in hots.values() for h in v]).to_csv(HOT_DIR / "supporting" / "collision_hotspots.csv",
                                                               index=False)
    by_domain.round(3).to_csv(HOT_DIR / "supporting" / "collision_hotspots_by_domain.csv", index=False)
    fig_heatmap(data, hots)
    fig_hotspots(data, hots)
    fig_alongshore(data, hots)
    fig_alongshore_series(data)
    fig_drivers(data, hots)


if __name__ == "__main__":
    main()
