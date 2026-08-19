r"""
HAT_road_domain_views.py
===============================================================================
One domain at a time: what CASCADE actually DOES to this grid.

Merged from three scripts that drew the same geometry three ways
(HAT_plot_road_initialization.py, HAT_plot_barrier3d_grid.py,
HAT_browse_road_domains.py). They shared a config block, a geometry() and a
drown test, and each copy was a place for those to drift. The drawing code below
is ported VERBATIM; only the config, the CLI and the setback selection are new.

MODES
-----
  --mode map        land/water plan view. The road as ONE row index applied to
                    all 50 profiles, against an island whose landward edge is
                    not straight, plus the two rows bulldoze tests.
  --mode section    the assembled Barrier3D cross-section -- shoreface, beach
                    wedge, dune rows, interior, bay -- with the REAL bulldoze()
                    run on a copy of the real interior.
  --mode both       both figures per domain (default).
  --browse          interactive walk instead of writing files (n/p/w/q).
                    Needs a GUI backend; everything else runs headless.

WHICH SETBACK IS DRAWN -- the figures are only honest if this matches the run
-----------------------------------------------------------------------------
  --method legacy      old_method_offset/<year>/RoadSetback_<year>.csv
                       what hatteras_site_config.py spends today
  --method dunestart   dunestart_offset/<year>/RoadSetback_<year>_dunestart.csv
                       measured landward of interior row 0, the reference
                       roadway_manager.py:99 actually uses

Both files are 2 rows x 82 cols, GIS IDs then metres, so this is a drop-in.
Change the default the SAME DAY you change hatteras_site_config.py:78,91.

USAGE
-----
    python HAT_road_domain_views.py --domains 52 --year 1984
    python HAT_road_domain_views.py --domains drowning --mode map
    python HAT_road_domain_views.py --browse --year 2004 --start 52
===============================================================================
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import matplotlib
# Backend is chosen in main(), not at import: --browse needs an interactive one
# and forcing Agg here would silently kill the window.
import matplotlib.pyplot as plt
from matplotlib.colors import (LinearSegmentedColormap, ListedColormap,
                               TwoSlopeNorm)
from matplotlib.patches import Rectangle, Patch
from matplotlib.lines import Line2D
from matplotlib.widgets import Button


def _find_project_root(start: Path) -> Path:
    """
    Walk up until a directory holds data\\hatteras_init.

    NOT parents[N]. These files have moved twice, and the old parents[4]
    silently resolved to scripts\\ -- every data path below it was wrong, and so
    was the sys.path entry the `import cascade.roadway_manager` depends on.
    """
    for p in [start, *start.parents]:
        if (p / "data" / "hatteras_init").is_dir():
            return p
    raise SystemExit(f"cannot find data\\hatteras_init above {start}")


PROJECT_ROOT = _find_project_root(Path(__file__).resolve())
sys.path.insert(0, str(PROJECT_ROOT))
import cascade.roadway_manager as rm          # noqa: E402  the real thing

# =============================================================================
# CONFIG -- every value traced to its source
# =============================================================================

DATA = PROJECT_ROOT / "data" / "hatteras_init"
ROADS_ROOT = DATA / "4-mgmt-forcing" / "road_offset"
# Topography version resolved from the extractor, not hardcoded -- it was
# "2009_v3" and kept drawing v3 interiors under v4 setbacks after the re-pick,
# with no error. See hat_topo_version.py.
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from hat_topo_version import topo_dirs  # noqa: E402

TOPO_DIR, DUNE_DIR, TOPO_RUN_NAME = topo_dirs()
RUN = TOPO_DIR.parent
# One file, both periods -- road elevation is a property of the 2009 surface.
ROAD_ELEV_CSV = DATA / "4-mgmt-forcing" / "road_elevation" / "RoadElevation.csv"

_SETBACK_SOURCES = {
    "legacy": ROADS_ROOT / "old_method_offset" / "{year}" / "RoadSetback_{year}.csv",
    "dunestart": (ROADS_ROOT / "dunestart_offset" / "{year}"
                  / "RoadSetback_{year}_dunestart.csv"),
}
# Matches hatteras_site_config.py:78,91. Change both together, or these
# figures stop describing the model you run.
SETBACK_METHOD = "dunestart"
# Rebound in main() from --method; the ported drawing code reads this global.
SETBACK = {y: Path(str(_SETBACK_SOURCES[SETBACK_METHOD]).format(year=y))
           for y in (1984, 2004)}

OUT_DIR = ROADS_ROOT / "figures"

# Hatteras-CASCADE-parameters.yaml
BERM_NAVD, MHW_M = 1.7, 0.36
BETA, BAY_DEPTH_M = 0.04, 3.0
DUNE_WIDTH_M = 20.0                       # barrier3d-default-parameters.yaml
BERM_MHW = BERM_NAVD - MHW_M              # load_input.py:241 -> 1.34 m MHW

# roadway_manager / runner
ROAD_WIDTH_M, DX, DY, DZ = 20.0, 10, 10, 10
DROWN_THRESHOLD, PCT_WATER = 0.0, 0.2     # roadway_manager.py:50-51, :531
SL_DAM = 0.0                              # barrier3d.py:1165, never changes
ROAD_ELEVATION_FALLBACK = 1.45
FIRST_ROAD_DOMAIN, LAST_ROAD_DOMAIN = 9, 90
# Villages: a RoadwayManager is constructed (cascade.py:428) but never updated.
VILLAGE_RANGES = [(21, 31), (68, 83)]

# --- colour -------------------------------------------------------------
C_ROAD, C_BAY, C_SEA = "#B71C1C", "#1565C0", "#FF8C00"   # validated set
C_INK = "#1a1a2e"
C_LAND, C_WATER = "#E8DCC0", "#B7D3E8"                   # recessive field
C_SKIP = "#9aa0a6"                                       # the untested row
C_DUNE_MARK, C_PROF = "#8D6E33", "#9aa0a6"

# Diverging about MHW: two hues, neutral midpoint, no rainbow.
ELEV_CMAP = LinearSegmentedColormap.from_list("mhw_diverging", [
    (0.00, "#0B3B60"), (0.30, "#2E75B6"), (0.46, "#B7D3E8"),
    (0.50, "#F2EFE6"),
    (0.56, "#D9C089"), (0.78, "#A9843F"), (1.00, "#5A4220"),
])


def is_village(gis: int) -> bool:
    return any(a <= gis <= b for a, b in VILLAGE_RANGES)


def load_interior(gis: int):
    fs = sorted(TOPO_DIR.glob(f"domain_{gis}_topography_*.npy"))
    return np.load(fs[0]) * DZ if fs else None      # dam -> m MHW


def load_setbacks(year: int) -> dict:
    p = SETBACK[year]
    if not p.exists():
        raise SystemExit(f"\n[stop] setback file not found:\n    {p}\n"
                         f"       (--method {SETBACK_METHOD})\n")
    a = np.loadtxt(p, delimiter=",")
    return dict(zip(a[0].astype(int), a[1]))


def geometry(interior: np.ndarray, setback_m: float) -> dict:
    """
    Reproduce bulldoze()'s indexing exactly.

    Superset of the two originals: the section view used only road_start /
    road_end / border / sea / bay / drowns, the map view also needs land, edge
    and water_frac. ONE implementation, so the drown verdict cannot differ
    between two figures of the same domain.
    """
    n_rows, n_prof = interior.shape
    road_start = int(setback_m / DY)
    road_end = road_start + int(ROAD_WIDTH_M / DX)
    border = road_end + 1

    land = interior > SL_DAM
    # where land stops, per profile -- barrier3d FindWidths' rule
    edge = np.array([next((i for i, v in enumerate(interior[:, c])
                           if v <= SL_DAM), n_rows) - 1
                     for c in range(n_prof)])
    edge = np.maximum(edge, 0)

    g = dict(n_rows=n_rows, n_prof=n_prof, road_start=road_start,
             road_end=road_end, border=border, land=land, edge=edge,
             water_frac=(~land).mean(axis=1), off_grid=border >= n_rows)
    g["bay"] = (1.0 if g["off_grid"]
                else float((interior[border, :] <= DROWN_THRESHOLD).mean()))
    g["sea"] = (float((interior[road_start - 1, :] <= DROWN_THRESHOLD).mean())
                if road_start > 0 else 0.0)
    g["drowns"] = g["sea"] > PCT_WATER or g["bay"] > PCT_WATER
    return g


# --- mode: map -------------------------------------------------------------

def plot_domain(gis: int, year: int, setback_m: float, out_dir: Path):
    interior = load_interior(gis)
    if interior is None:
        print(f"  [skip] GIS {gis}: no topography array")
        return None
    g = geometry(interior, setback_m)
    n_rows, n_prof = g["n_rows"], g["n_prof"]

    fig = plt.figure(figsize=(14, 8.6), facecolor="white")
    gs = fig.add_gridspec(2, 2, height_ratios=[3.0, 1.0],
                          width_ratios=[2.5, 1.0], hspace=0.28, wspace=0.16)
    ax = fig.add_subplot(gs[0, 0])
    ax_edge = fig.add_subplot(gs[0, 1], sharey=ax)
    ax_w = fig.add_subplot(gs[1, :])

    # ---- the field: land vs water, deliberately recessive ----
    ax.imshow(g["land"].astype(int), cmap=ListedColormap([C_WATER, C_LAND]),
              aspect="auto", interpolation="nearest", origin="upper",
              extent=[-0.5, n_prof - 0.5, n_rows - 0.5, -0.5], zorder=1)

    # ---- the island's landward edge: the jagged line ----
    ax.step(np.arange(n_prof), g["edge"], where="mid", color=C_INK, lw=2.0,
            zorder=4)

    # ---- the road: one row index, applied to every profile ----
    ax.axhspan(g["road_start"] - 0.5, g["road_end"] - 0.5, color=C_ROAD,
               alpha=0.30, zorder=3, lw=0)
    for r in (g["road_start"] - 0.5, g["road_end"] - 0.5):
        ax.axhline(r, color=C_ROAD, lw=2.0, zorder=5)

    # ---- the row bulldoze skips ----
    ax.axhline(g["road_end"], color=C_SKIP, lw=2.0, ls=(0, (1, 2)), zorder=5)

    # ---- the two rows bulldoze actually tests ----
    if not g["off_grid"]:
        ax.axhline(g["border"], color=C_BAY, lw=2.0, zorder=6)
    if g["road_start"] > 0:
        ax.axhline(g["road_start"] - 1, color=C_SEA, lw=2.0, zorder=6)

    # Direct labels, placed INSIDE the axes. The sea-side colour WARNed on
    # contrast against a light surface in the validator, so it carries a label
    # rather than relying on the legend -- and the label sits on an opaque
    # white plate so the pale field underneath cannot erode it further.
    # (Placing these outside the axes clipped them against the right panel.)
    lab = dict(ha="right", va="center", fontsize=8.5, fontweight="bold",
               zorder=8,
               bbox=dict(boxstyle="round,pad=0.28", fc="white", ec="none",
                         alpha=0.96))
    x_lab = n_prof - 1.2
    if g["road_start"] > 0:
        ax.text(x_lab, g["road_start"] - 1,
                f"sea-side test — {g['sea'] * 100:.0f}% water",
                color=C_SEA, **lab)
    if not g["off_grid"]:
        ax.text(x_lab, g["border"],
                f"bay-side test — {g['bay'] * 100:.0f}% water",
                color=C_BAY, **lab)
    ax.text(x_lab, (g["road_start"] + g["road_end"] - 1) / 2,
            f"road — rows {g['road_start']}–{g['road_end'] - 1}",
            color=C_ROAD, **lab)

    ax.set_xlim(-0.5, n_prof - 0.5)
    ax.set_ylim(n_rows - 0.5, -0.5)
    ax.set_xlabel("alongshore profile  (50 cells = 500 m)")
    ax.set_ylabel("cross-shore row  (0 = behind the dune  ->  landward)")
    ax.set_title(f"The road is one row index applied to every profile;\n"
                 f"the island's landward edge is not a straight line",
                 fontsize=10, color=C_INK, pad=8)

    # ---- right: how many profiles END at each row ----
    # y is the SAME cross-shore row as the map (sharey), so a bar here lines up
    # with the row it describes. Counting profiles rather than tracing the edge
    # again keeps this panel from restating the map's black line.
    counts = np.bincount(g["edge"], minlength=n_rows)[:n_rows]
    ax_edge.barh(np.arange(n_rows), counts, height=0.85, color=C_INK,
                 zorder=3)
    ax_edge.axhspan(g["road_start"] - 0.5, g["road_end"] - 0.5, color=C_ROAD,
                    alpha=0.30, lw=0, zorder=1)
    short = int((g["edge"] < g["road_start"]).sum())
    ax_edge.set_xlabel("profiles whose land ends at this row")
    ax_edge.set_title(f"{short} of {n_prof} profiles run out of\nland before "
                      f"the road", fontsize=9.5, color=C_INK, pad=8)
    ax_edge.grid(alpha=0.25, axis="x")
    ax_edge.tick_params(labelleft=False)   # shares the main panel's y axis

    # ---- bottom: water fraction down the cross-shore ----
    rows = np.arange(n_rows)
    ax_w.fill_between(rows, g["water_frac"] * 100, color=C_WATER, zorder=1)
    ax_w.plot(rows, g["water_frac"] * 100, color=C_INK, lw=1.6, zorder=2)
    ax_w.axhline(PCT_WATER * 100, color=C_INK, ls="--", lw=1.4, zorder=3)
    ax_w.text(0.4, PCT_WATER * 100 + 3,
              f"bulldoze() drowns the road above {PCT_WATER * 100:.0f}%",
              fontsize=8.5, color=C_INK)
    ax_w.axvspan(g["road_start"] - 0.5, g["road_end"] - 0.5, color=C_ROAD,
                 alpha=0.30, lw=0, zorder=0)
    if not g["off_grid"]:
        ax_w.axvline(g["border"], color=C_BAY, lw=2.0, zorder=4)
    if g["road_start"] > 0:
        ax_w.axvline(g["road_start"] - 1, color=C_SEA, lw=2.0, zorder=4)
    ax_w.set_xlim(-0.5, n_rows - 0.5)
    ax_w.set_ylim(0, 105)
    ax_w.set_xlabel("cross-shore row")
    ax_w.set_ylabel("% of profiles\nthat are water")
    ax_w.grid(alpha=0.25)

    # ---- legend: identity is never colour-alone ----
    handles = [
        Patch(facecolor=C_LAND, label="land (> 0 m MHW)"),
        Patch(facecolor=C_WATER, label="water (<= 0 m MHW: sound, marsh, "
                                       "or NoData sentinel)"),
        Patch(facecolor=C_ROAD, alpha=0.30, label="road footprint"),
        Line2D([], [], color=C_BAY, lw=2.0, label="bay-side test row"),
        Line2D([], [], color=C_SEA, lw=2.0, label="sea-side test row"),
        Line2D([], [], color=C_SKIP, lw=2.0, ls=(0, (1, 2)),
               label="road_end: never tested"),
        Line2D([], [], color=C_INK, lw=2.0, label="island's landward edge"),
    ]
    ax.legend(handles=handles, loc="upper left", fontsize=8, framealpha=0.94,
              ncol=2)

    verdict = "DROWNS in year 1" if g["drowns"] else "road fits"
    vcol = C_ROAD if g["drowns"] else C_INK
    note = " (village: manager built but never updated)" if is_village(gis) else ""
    fig.suptitle(
        f"GIS {gis}  |  {year} setback {setback_m:.0f} m  |  {verdict}{note}\n"
        f"array is {n_rows} rows tall, but land ends at row "
        f"{g['edge'].min()}-{g['edge'].max()} (median {int(np.median(g['edge']))})",
        fontsize=13, fontweight="bold", color=vcol, y=0.98)

    out_dir.mkdir(parents=True, exist_ok=True)
    p = out_dir / f"road_init_GIS{gis:03d}_{year}.png"
    fig.savefig(p, dpi=150, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    return p, g


def plot_overview(year: int, setbacks: dict, out_dir: Path):
    gis, road, e_lo, e_med, e_hi, drown, managed = [], [], [], [], [], [], []
    for d in range(FIRST_ROAD_DOMAIN, LAST_ROAD_DOMAIN + 1):
        interior = load_interior(d)
        if interior is None or d not in setbacks:
            continue
        g = geometry(interior, float(setbacks[d]))
        gis.append(d)
        road.append(g["road_start"] * DY)
        e_lo.append(g["edge"].min() * DY)
        e_med.append(np.median(g["edge"]) * DY)
        e_hi.append(g["edge"].max() * DY)
        drown.append(g["drowns"])
        managed.append(not is_village(d))
    gis = np.array(gis); road = np.array(road, float)
    e_lo = np.array(e_lo, float); e_med = np.array(e_med, float)
    e_hi = np.array(e_hi, float)
    drown = np.array(drown); managed = np.array(managed)

    fig, ax = plt.subplots(figsize=(15, 6.5), facecolor="white")
    ax.fill_between(gis, e_lo, e_hi, color=C_LAND, zorder=1,
                    label="island: thinnest to widest profile")
    ax.plot(gis, e_med, color=C_INK, lw=2.0, zorder=3,
            label="island: median profile")
    ax.plot(gis, road, color=C_ROAD, lw=2.0, zorder=4, label="road position")

    bad = drown & managed
    if bad.any():
        ax.plot(gis[bad], road[bad], lw=0, marker="v", ms=10,
                color=C_ROAD, mec="white", mew=1.2, zorder=6,
                label="drowns in year 1 (managed)")
        for d, r in zip(gis[bad], road[bad]):
            ax.annotate(str(d), (d, r), xytext=(0, 12),
                        textcoords="offset points", ha="center", fontsize=8,
                        fontweight="bold", color=C_ROAD)
    vil = drown & ~managed
    if vil.any():
        ax.plot(gis[vil], road[vil], lw=0, marker="v", ms=9, mfc="none",
                color=C_ROAD, mew=1.6, zorder=5,
                label="would drown, but village (never updated)")

    ax.set_xlabel("GIS domain  (9 = Buxton  ->  90 = Rodanthe)")
    ax.set_ylabel("distance behind the dune line (m)")
    ax.set_title(f"Where the {year} road sits, against how far the 2009 island "
                 f"actually extends\n"
                 f"Above the shaded band = the road is behind the island",
                 fontsize=12, fontweight="bold", color=C_INK)
    ax.legend(loc="upper left", fontsize=9, framealpha=0.94)
    ax.grid(alpha=0.25)
    ax.set_xlim(gis.min() - 0.5, gis.max() + 0.5)

    out_dir.mkdir(parents=True, exist_ok=True)
    p = out_dir / f"road_init_overview_{year}.png"
    fig.savefig(p, dpi=150, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    return p


# --- mode: section ---------------------------------------------------------

def load_grid(gis: int):
    t = sorted(TOPO_DIR.glob(f"domain_{gis}_topography_*.npy"))
    d = sorted(DUNE_DIR.glob(f"domain_{gis}_dune_*.npy"))
    if not t:
        return None
    interior = np.load(t[0]) * DZ                       # dam -> m MHW
    dune_h = np.load(d[0]) * DZ if d else None          # dam -> m ABOVE BERM
    return interior, dune_h


def road_elevation(gis: int) -> float:
    if not ROAD_ELEV_CSV.exists():
        return ROAD_ELEVATION_FALLBACK
    a = np.loadtxt(ROAD_ELEV_CSV, delimiter=",")
    m = dict(zip(a[0].astype(int), a[1]))
    return float(m.get(gis, ROAD_ELEVATION_FALLBACK))


def run_real_bulldoze(interior, dune_h, road_ele, setback_m):
    """Call cascade's bulldoze() on a COPY -- it mutates the grid in place."""
    grid = (interior / DZ).copy()                       # back to dam
    n_dune_cells = max(1, int(DUNE_WIDTH_M / DX))
    dune = np.tile((dune_h / DZ)[:, None], (1, n_dune_cells)) if dune_h is not None \
        else np.full((interior.shape[1], n_dune_cells), 0.02)
    new_dune, new_grid, removed, drown = rm.bulldoze(
        time_index=1, xyz_interior_grid=grid, yxz_dune_grid=dune.copy(),
        road_ele=road_ele, road_width=ROAD_WIDTH_M, road_setback=setback_m,
        dx=DX, dy=DY, dz=DZ, drown_threshold=DROWN_THRESHOLD,
        percent_water_cells_touching_road=PCT_WATER,
    )
    return new_grid * DZ, new_dune * DZ, float(removed), bool(drown)


def plot_grid(gis, year, setback_m, out_dir):
    got = load_grid(gis)
    if got is None:
        print(f"  [skip] GIS {gis}: no arrays")
        return None
    interior, dune_h = got
    g = geometry(interior, setback_m)
    road_ele = road_elevation(gis)
    n_rows, n_prof = g["n_rows"], g["n_prof"]

    n_beach = max(1, int(BERM_MHW / BETA / DX))
    n_dune = max(1, int(DUNE_WIDTH_M / DX))
    x_int = np.arange(n_rows) * DX                     # interior cross-shore, m
    x_dune = -np.arange(n_dune, 0, -1) * DX
    x_beach = x_dune[0] - np.arange(n_beach, 0, -1) * DX

    after, dune_after, removed, drown = run_real_bulldoze(
        interior, dune_h, road_ele, setback_m)

    fig = plt.figure(figsize=(15, 11.5), facecolor="white")
    gsp = fig.add_gridspec(3, 1, height_ratios=[1.15, 1.0, 1.0], hspace=0.42)
    ax_map, ax_sec, ax_bd = (fig.add_subplot(gsp[i, 0]) for i in range(3))

    # ---------------- 1. the grid as an array, elevation-coloured ----------
    dune_elev = (dune_h + BERM_MHW) if dune_h is not None else None
    block = interior.T.copy()                           # (profiles, rows)
    if dune_elev is not None:
        block = np.hstack([np.tile(dune_elev[:, None], (1, n_dune)), block])
    x0 = (x_dune[0] if dune_elev is not None else x_int[0]) - DX / 2
    vmax = max(2.0, float(np.nanpercentile(interior, 99.5)))
    norm = TwoSlopeNorm(vmin=min(-BAY_DEPTH_M, float(interior.min())),
                        vcenter=0.0, vmax=vmax)
    im = ax_map.imshow(block, cmap=ELEV_CMAP, norm=norm, aspect="auto",
                       origin="lower", interpolation="nearest",
                       extent=[x0, x_int[-1] + DX / 2, -0.5, n_prof - 0.5])
    if dune_elev is not None:
        ax_map.axvline(-DX / 2, color=C_INK, lw=1.6, ls="--")
        # inside the axes: at n_prof*1.02 this collided with the panel title
        ax_map.text(x_dune[0], n_prof * 0.94, "dune\ncells", fontsize=8,
                    color="white", fontweight="bold", ha="center", va="top",
                    zorder=6)
    ax_map.add_patch(Rectangle((g["road_start"] * DX - DX / 2, -0.5),
                               ROAD_WIDTH_M, n_prof, facecolor="none",
                               edgecolor=C_ROAD, lw=2.2, zorder=5))
    ax_map.set_xlabel("cross-shore distance from the dune toe (m)")
    ax_map.set_ylabel("alongshore profile")
    ax_map.set_title(f"1 · The grid CASCADE holds — dune cells + interior, "
                     f"coloured by elevation about MHW", fontsize=11,
                     color=C_INK, loc="left")
    cb = fig.colorbar(im, ax=ax_map, pad=0.01, fraction=0.03)
    cb.set_label("elevation (m, MHW-relative)", fontsize=9)

    # ---------------- 2. the cross-shore section --------------------------
    ax_sec.axhspan(norm.vmin, 0, color="#DCEAF6", zorder=0)
    ax_sec.axhline(0, color=C_BAY, lw=1.6, zorder=2)
    ax_sec.text(x_beach[0], 0.06, "MHW", color=C_BAY, fontsize=8.5,
                fontweight="bold", va="bottom")

    ax_sec.fill_between(x_beach, 0, np.linspace(0, BERM_MHW, n_beach),
                        color="#EFE3C4", zorder=1)
    ax_sec.text(x_beach.mean(), BERM_MHW * 0.35,
                f"beach\nint(BermEl/β)={n_beach} cells", fontsize=7.5,
                ha="center", color=C_INK, zorder=3)

    for p in range(n_prof):                              # every profile, faint
        ax_sec.step(x_int, interior[:, p], where="mid", color=C_PROF,
                    lw=0.5, alpha=0.55, zorder=2)
    med = np.median(interior, axis=1)
    ax_sec.step(x_int, med, where="mid", color=C_INK, lw=2.2, zorder=5,
                label="interior, alongshore median")
    if dune_elev is not None:
        ax_sec.step(np.r_[x_dune, x_int[0]],
                    np.r_[np.full(n_dune, np.median(dune_elev)), med[0]],
                    where="mid", color=C_DUNE_MARK, lw=2.2, zorder=5,
                    label=f"dune (median crest "
                          f"{np.median(dune_elev):.2f} m MHW)")
    ax_sec.axhline(-BAY_DEPTH_M, color=C_BAY, lw=1.4, ls=":", zorder=2)
    ax_sec.text(x_int[-1], -BAY_DEPTH_M + 0.08,
                f"BayDepth {BAY_DEPTH_M:.0f} m", color=C_BAY, fontsize=8,
                ha="right", va="bottom")

    ax_sec.add_patch(Rectangle((g["road_start"] * DX - DX / 2, road_ele - 0.12),
                               ROAD_WIDTH_M, 0.24, facecolor=C_ROAD,
                               edgecolor="white", lw=1.0, zorder=7))
    # The road label sits BELOW the road box and to the left of the test-row
    # labels; placing it above collided with the rotated sea-side label.
    ax_sec.annotate(f"road — {road_ele:.2f} m MHW, setback {setback_m:.0f} m",
                    xy=(g["road_start"] * DX + ROAD_WIDTH_M / 2, road_ele),
                    xytext=(g["road_start"] * DX - 30, road_ele - 1.6),
                    color=C_ROAD, fontsize=9, fontweight="bold", ha="right",
                    arrowprops=dict(arrowstyle="->", color=C_ROAD, lw=1.6,
                                    connectionstyle="arc3,rad=-0.2"),
                    bbox=dict(boxstyle="round,pad=0.28", fc="white",
                              ec="none", alpha=0.95), zorder=9)
    for x, c, txt in ((g["road_start"] - 1, C_SEA, f"sea-side {g['sea']*100:.0f}%"),
                      (g["border"], C_BAY, f"bay-side {g['bay']*100:.0f}%")):
        if 0 <= x < n_rows:
            ax_sec.axvline(x * DX, color=c, lw=2.0, zorder=6)
            ax_sec.text(x * DX, vmax + 0.45, txt, color=c, fontsize=8.5,
                        fontweight="bold", ha="center", va="top",
                        bbox=dict(boxstyle="round,pad=0.22", fc="white",
                                  ec="none", alpha=0.95), zorder=9)
    ax_sec.set_xlim(x_beach[0], x_int[-1] + DX / 2)
    ax_sec.set_ylim(min(-BAY_DEPTH_M - 0.4, interior.min() - 0.2), vmax + 0.6)
    ax_sec.set_xlabel("cross-shore distance from the dune toe (m)")
    ax_sec.set_ylabel("elevation (m MHW)")
    ax_sec.set_title("2 · The same grid as a cross-shore section — all 50 "
                     "profiles, with the road drawn where CASCADE puts it",
                     fontsize=11, color=C_INK, loc="left")
    ax_sec.legend(loc="upper right", fontsize=8.5, framealpha=0.94)
    ax_sec.grid(alpha=0.2)

    # ---------------- 3. what the real bulldoze() did ---------------------
    med_a = np.median(after, axis=1)
    ax_bd.axhspan(norm.vmin, 0, color="#DCEAF6", zorder=0)
    ax_bd.axhline(0, color=C_BAY, lw=1.4, zorder=2)
    ax_bd.step(x_int, med, where="mid", color=C_PROF, lw=2.4, zorder=4,
               label="before bulldoze()")
    ax_bd.step(x_int, med_a, where="mid", color=C_INK, lw=2.2, zorder=5,
               label="after bulldoze()")
    lo, hi = np.minimum(med, med_a), np.maximum(med, med_a)
    ax_bd.fill_between(x_int, lo, hi, step="mid", color=C_ROAD, alpha=0.30,
                       zorder=3, label="cells bulldoze() rewrote")
    ax_bd.axhline(road_ele, color=C_ROAD, lw=1.4, ls="--", zorder=6)
    ax_bd.text(x_int[-1], road_ele + 0.05, f"road_ele {road_ele:.2f}",
               color=C_ROAD, fontsize=8.5, ha="right", va="bottom")
    ax_bd.set_xlim(x_int[0] - DX / 2, x_int[-1] + DX / 2)
    ax_bd.set_xlabel("cross-shore distance from the dune toe (m)")
    ax_bd.set_ylabel("elevation (m MHW)")
    # Count over the ACTUAL rewritten block (2 road rows x 50 profiles), not
    # the median profile -- summing metres down a median is not a quantity.
    blk_b = interior[g["road_start"]:g["road_end"], :]
    blk_a = after[g["road_start"]:g["road_end"], :]
    n_up = int((blk_a > blk_b + 1e-9).sum())
    n_dn = int((blk_a < blk_b - 1e-9).sum())
    n_cells = blk_b.size
    was_water = int((blk_b <= 0).sum())
    ax_bd.set_title(
        f"3 · What the real bulldoze() did in year 1 — "
        f"{'roadway_drown = True' if drown else 'road survives'} · "
        f"sand to dune {removed * 1000:,.0f} m³\n"
        f"of the {n_cells} road cells: {n_up} raised, {n_dn} lowered to "
        f"road_ele — and {was_water} of them were below MHW beforehand",
        fontsize=11, color=C_ROAD if drown else C_INK, loc="left")
    ax_bd.legend(loc="upper right", fontsize=8.5, framealpha=0.94)
    ax_bd.grid(alpha=0.2)

    verdict = "DROWNS in year 1" if g["drowns"] else "road fits"
    note = "  (village: manager built but never updated)" if is_village(gis) else ""
    fig.suptitle(f"GIS {gis}  ·  {year} setback {setback_m:.0f} m  ·  {verdict}{note}",
                 fontsize=14, fontweight="bold",
                 color=C_ROAD if g["drowns"] else C_INK, y=0.995)

    out_dir.mkdir(parents=True, exist_ok=True)
    p = out_dir / f"b3d_grid_GIS{gis:03d}_{year}.png"
    fig.savefig(p, dpi=140, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    return p, g, drown, removed


# --- mode: browse ----------------------------------------------------------

def ensure_interactive_backend():
    """
    Get a backend that opens a REAL window and blocks in plt.show().

    PyCharm is the reason this exists. With "Show plots in tool window" on
    (Settings > Tools > Python Scientific), PyCharm swaps the backend for
    'module://backend_interagg', which draws into the SciView panel instead of a
    window. plt.show() then returns IMMEDIATELY and no key_press_event is ever
    delivered -- so the browser renders the first domain, never receives a
    keystroke, and exits. That looks like "it only plots GIS 9".

    Jupyter's 'module://matplotlib_inline.backend_inline' behaves the same way.

    Testing for the literal string 'agg' does not catch either of them, so match
    on the module:// prefix instead and switch to a windowing backend.
    """
    current = matplotlib.get_backend()
    low = current.lower()
    ok = not (low == "agg" or low.startswith("module://") or "inline" in low)
    if ok:
        return current, None

    for cand in ("TkAgg", "QtAgg", "Qt5Agg", "WxAgg"):
        try:
            matplotlib.use(cand, force=True)
            importlib.reload(plt)
            return matplotlib.get_backend(), current
        except Exception:
            continue
    return current, current       # could not switch; caller warns


def draw(fig, axes, gis, year, setback_m):
    """Redraw both panels in place. Returns the geometry dict, or None."""
    ax_map, ax_prof = axes
    ax_map.clear(); ax_prof.clear()

    got = load_grid(gis)
    if got is None:
        ax_map.text(0.5, 0.5, f"GIS {gis}: no topography array",
                    ha="center", va="center", transform=ax_map.transAxes)
        return None
    interior, dune_h = got
    g = geometry(interior, setback_m)
    road_ele = road_elevation(gis)
    n_rows, n_prof = g["n_rows"], g["n_prof"]

    vmax = max(2.0, float(np.nanpercentile(interior, 99.5)))
    norm = matplotlib.colors.TwoSlopeNorm(
        vmin=min(-BAY_DEPTH_M, float(interior.min())), vcenter=0.0, vmax=vmax)

    # ---- map: alongshore x, cross-shore y, row 0 at the bottom ----
    ax_map.imshow(interior, cmap=ELEV_CMAP, norm=norm, aspect="auto",
                  origin="lower", interpolation="nearest",
                  extent=[-0.5, n_prof - 0.5, -0.5, n_rows - 0.5])

    edge = np.array([next((i for i, v in enumerate(interior[:, c]) if v <= 0),
                          n_rows) - 1 for c in range(n_prof)])
    edge = np.maximum(edge, 0)
    ax_map.step(np.arange(n_prof), edge, where="mid", color=C_INK, lw=2.0,
                zorder=5, label="island's landward edge")

    ax_map.axhspan(g["road_start"] - 0.5, g["road_end"] - 0.5,
                   color=C_ROAD, alpha=0.30, zorder=3, lw=0)
    ax_map.axhline(g["road_end"], color="#9aa0a6", lw=1.8, ls=(0, (1, 2)),
                   zorder=4)
    if g["road_start"] > 0:
        ax_map.axhline(g["road_start"] - 1, color=C_SEA, lw=2.0, zorder=6)
    if not g["off_grid"]:
        ax_map.axhline(g["border"], color=C_BAY, lw=2.0, zorder=6)

    ax_map.set_xlabel("alongshore cell  (50 = 500 m)")
    ax_map.set_ylabel("cross-shore cell  (0 = behind the dune, landward up)")
    ax_map.set_xlim(-0.5, n_prof - 0.5)
    ax_map.set_ylim(-0.5, n_rows - 0.5)
    handles = [
        Patch(facecolor=C_ROAD, alpha=0.30, label="road cells"),
        Line2D([], [], color=C_SEA, lw=2.0,
               label=f"sea-side test — {g['sea'] * 100:.0f}% water"),
        Line2D([], [], color=C_BAY, lw=2.0,
               label=f"bay-side test — {g['bay'] * 100:.0f}% water"),
        Line2D([], [], color="#9aa0a6", lw=1.8, ls=(0, (1, 2)),
               label="road_end: never tested"),
        Line2D([], [], color=C_INK, lw=2.0, label="island's landward edge"),
    ]
    ax_map.legend(handles=handles, loc="upper left", fontsize=8,
                  framealpha=0.94)

    # ---- profile: elevation x, cross-shore y (shared) ----
    y = np.arange(n_rows)
    ax_prof.plot(interior, y, color="0.78", lw=0.6, zorder=2)
    med = np.median(interior, axis=1)
    ax_prof.plot(med, y, color=C_INK, lw=2.2, zorder=5,
                 label="median profile")
    ax_prof.axvline(0, color=C_BAY, lw=1.6, zorder=3, label="MHW")
    ax_prof.axvline(road_ele, color=C_ROAD, lw=1.6, ls="--", zorder=4,
                    label=f"road_ele {road_ele:.2f} m")
    ax_prof.axhspan(g["road_start"] - 0.5, g["road_end"] - 0.5,
                    color=C_ROAD, alpha=0.30, zorder=1, lw=0)
    if g["road_start"] > 0:
        ax_prof.axhline(g["road_start"] - 1, color=C_SEA, lw=2.0, zorder=6)
    if not g["off_grid"]:
        ax_prof.axhline(g["border"], color=C_BAY, lw=2.0, zorder=6)
    ax_prof.set_xlabel("elevation (m MHW)")
    ax_prof.set_xlim(min(-BAY_DEPTH_M - 0.3, interior.min() - 0.2),
                     vmax + 0.4)
    ax_prof.legend(loc="upper right", fontsize=8, framealpha=0.9)
    ax_prof.grid(alpha=0.25)
    plt.setp(ax_prof.get_yticklabels(), visible=False)

    # ---- the numbers, on the figure so a screenshot carries them ----
    short = int((edge < g["road_start"]).sum())
    verdict = "DROWNS in year 1" if g["drowns"] else "road fits"
    vcol = C_ROAD if g["drowns"] else C_INK
    village = "  ·  village: manager built but never updated" \
        if is_village(gis) else ""
    fig.suptitle(
        f"GIS {gis}  ·  {year} setback {setback_m:.0f} m  ·  road_ele "
        f"{road_ele:.2f} m MHW  ·  {verdict}{village}\n"
        f"road rows {g['road_start']}–{g['road_end'] - 1} of {n_rows}  ·  "
        f"land ends at row {edge.min()}–{edge.max()} (median "
        f"{int(np.median(edge))})  ·  {short}/{n_prof} profiles run out of land "
        f"before the road\n"
        f"n/→ next   p/← prev   d next drowning   w write PNG   q quit",
        fontsize=11, color=vcol)
    return g



# =============================================================================
# CLI
# =============================================================================

def resolve_domains(spec: str, setbacks: dict) -> list:
    """`all`, `drowning`, or an explicit comma / range list."""
    span = [d for d in range(FIRST_ROAD_DOMAIN, LAST_ROAD_DOMAIN + 1)
            if d in setbacks]
    if spec == "all":
        return span
    if spec == "drowning":
        out = []
        for d in span:
            interior = load_interior(d)
            if interior is not None and geometry(interior, setbacks[d])["drowns"]:
                out.append(d)
        return out
    got = []
    for part in spec.split(","):
        part = part.strip()
        if "-" in part:
            a, b = part.split("-")
            got += list(range(int(a), int(b) + 1))
        elif part:
            got.append(int(part))
    return [d for d in got if d in setbacks]


def browse(domains, year, setbacks, start, out_dir):
    """Keyboard walk: n/p next/prev, w write png, q quit."""
    idx = domains.index(start) if start in domains else 0
    fig, axes = plt.subplots(2, 1, figsize=(13, 9),
                             gridspec_kw=dict(height_ratios=[3, 1]))
    state = {"i": idx}

    def show():
        gis = domains[state["i"]]
        for ax in np.atleast_1d(axes).ravel():
            ax.clear()
        draw(fig, axes, gis, year, setbacks[gis])
        fig.canvas.draw_idle()

    def on_key(ev):
        if ev.key in ("n", "right"):
            state["i"] = (state["i"] + 1) % len(domains)
            show()
        elif ev.key in ("p", "left"):
            state["i"] = (state["i"] - 1) % len(domains)
            show()
        elif ev.key == "w":
            gis = domains[state["i"]]
            p = Path(out_dir) / f"browse_GIS{gis:03d}_{year}.png"
            fig.savefig(p, dpi=150, bbox_inches="tight", facecolor="white")
            print(f"  [out] {p}")
        elif ev.key == "q":
            plt.close(fig)

    fig.canvas.mpl_connect("key_press_event", on_key)
    show()
    print("  n/p = next/prev   w = write png   q = quit")
    plt.show()
    return 0


def main() -> int:
    # Declared up front: --method's default READS SETBACK_METHOD, and Python
    # rejects a global statement that comes after a use in the same scope.
    global SETBACK_METHOD, SETBACK

    ap = argparse.ArgumentParser(
        description="Per-domain views of where NC-12 lands in CASCADE.")
    ap.add_argument("--mode", default="both",
                    choices=["map", "section", "both"])
    ap.add_argument("--browse", action="store_true",
                    help="interactive walk; needs a GUI backend")
    ap.add_argument("--method", default=SETBACK_METHOD,
                    choices=sorted(_SETBACK_SOURCES))
    ap.add_argument("--year", type=int, default=1984, choices=[1984, 2004])
    ap.add_argument("--domains", default="drowning",
                    help="all | drowning | 52 | 9-15,52,74")
    ap.add_argument("--start", type=int, default=None,
                    help="--browse only: domain to open on")
    ap.add_argument("--out", default=str(OUT_DIR))
    a = ap.parse_args()

    SETBACK_METHOD = a.method
    SETBACK = {y: Path(str(_SETBACK_SOURCES[a.method]).format(year=y))
               for y in (1984, 2004)}

    if a.browse:
        ensure_interactive_backend()
    else:
        matplotlib.use("Agg")

    setbacks = load_setbacks(a.year)
    domains = resolve_domains(a.domains, setbacks)
    out_dir = Path(a.out)
    out_dir.mkdir(parents=True, exist_ok=True)

    print("=" * 84)
    print(f"NC-12 per-domain views  |  {a.year}  |  --method {a.method}")
    print("=" * 84)
    print(f"  setback file : {SETBACK[a.year]}")
    print(f"  topography   : {RUN.name}")
    shown = domains if len(domains) <= 20 else domains[:20] + ["..."]
    print(f"  domains      : {len(domains)} -> {shown}")
    if not domains:
        print("\n  nothing to draw.")
        return 0

    if a.browse:
        return browse(domains, a.year, setbacks, a.start, out_dir)

    for d in domains:
        if a.mode in ("map", "both"):
            plot_domain(d, a.year, setbacks[d], out_dir)
        if a.mode in ("section", "both"):
            plot_grid(d, a.year, setbacks[d], out_dir)
    if a.mode in ("map", "both") and len(domains) > 1:
        plot_overview(a.year, setbacks, out_dir)
    print(f"\n  written to {out_dir}")
    return 0


if __name__ == "__main__":
    sys.exit(main())

