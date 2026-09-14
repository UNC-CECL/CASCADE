"""
HAT_plot_topo_retained.py

What the DEM's unsurveyed ground becomes once it is a CASCADE input domain, and
what Barrier3D does with it at t = 0.

THE ANSWER TO "ARE THERE CELLS THAT STAY NO-DATA AT MODEL START"
-----------------------------------------------------------------
No. There is no such state to stay in. Barrier3D has one float per cell and no
representation for "unknown", so the extractor writes every unsurveyed cell to
SENTINEL_WATER_M and CASCADE reads it as an elevation of exactly -3.0 m MHW.
The `<stem>_nodata.npy` mask that records which cells those were is a sidecar:
hat_topo_version.domain_arrays() hands Cascade() the topography and dune paths
only, so nothing in the run ever opens it.

So the question is not whether no-data survives. It is what the model believes
instead, and the answer is: open water, at the bottom of the clamp, in the
middle of the barrier.

WHY THAT IS NOT COSMETIC
------------------------
barrier3d.FindWidths - transcribed in roadway.interior_widths - measures the
island as the run of land from interior row 0 to THE FIRST WATER CELL. Land
behind a water cell is invisible to it. One unsurveyed cell at row k therefore
truncates that profile's island at row k-1 and discards every real, measured
cell behind it.

Panel (c) is that cost. It is an UPPER BOUND on the damage, not an estimate:
it compares the width Barrier3D actually sees against the width it would see if
every unsurveyed cell turned out to be land. Nobody knows that they are - that
is what unsurveyed means - so the true loss is somewhere between zero and the
bar drawn. The bound is still worth having, because it is the number that would
have to be small for the truncation not to matter.

The same conflation is what drowned three roadways at t = 0 in an earlier
product: roadway_manager.bulldoze drowns a road when more than 20% of the cells
flanking it sit at or below 0 m MHW, and an unsurveyed cell passes that test.
predict_drowning() is run here against this product's own setbacks and the
verdict is printed.

THE PANELS
----------
a  The CASCADE input domain for all 90 domains: the 2 dune rows Barrier3D
   builds from dunes/domain_<N>_dune.npy, then the interior rows from
   topography/. Ocean at the bottom. This is the whole stack the model starts
   from, in the order it starts from it.

b  The seaward 300 m of the same stack, so the dune rows and the first interior
   rows are actually resolvable. At island scale two 10 m rows are one pixel.

c  Island width Barrier3D sees, and the upper bound on what unsurveyed cells
   cost it.

d  Unsurveyed cells per domain: in the DEM, and still there in the input domain.

DUNE ROWS
---------
dunes/domain_<N>_dune.npy is one height above berm per profile, (50,). Barrier3D
runs DuneWidth = 2, and row 1 is a copy of row 0 - see the dune-rows note in
the extractor. Both rows are drawn. Their elevation is BERM_ELEV + height; the
berm is 1.7 m NAVD88, so 1.34 m MHW.

No dune cell in this product is unsurveyed: all 4500 carry a measured height.
That is checked at run time, not assumed, and the count is printed.

INPUT   <product>/npy-arrays/domain_<N>.npy                    m NAVD88
        <product>/dune-topo/<version>/topography/domain_<N>_topography.npy   dam
        <product>/dune-topo/<version>/topography/domain_<N>_nodata.npy      bool
        <product>/dune-topo/<version>/dunes/domain_<N>_dune.npy             dam

        Product and version resolve through scripts/hat_topo_version.py.
        Do not hardcode either.

OUTPUT  <product>/dune-topo/<version>/figures/HAT_topo_retained_<version>.png
"""

import csv
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import BoundaryNorm, ListedColormap
from matplotlib.patches import Patch

REPO = next(
    _p for _p in Path(__file__).resolve().parents
    if (_p / "pyproject.toml").exists())   # 1-extraction/nodata_audit/ since 2026-09-09
sys.path.insert(0, str(REPO / "scripts"))
import hat_topo_version as htv  # noqa: E402
from cascade_pipeline import roadway  # noqa: E402
from hat_elevation_products import product as elevation_product  # noqa: E402

# =============================================================================
# CONFIG
# =============================================================================

TOPO_PRODUCT = "1984-start"
DEM_PRODUCT = "2009-2014-1996"     # alongshore georeferencing only
VERSION_OVERRIDE = None            # None -> hat_topo_version resolves it

MHW_M = 0.36
BERM_NAVD_M = 1.7                  # RUN_MANIFEST BERM_ELEV_NAVD_M
BERM_MHW_M = BERM_NAVD_M - MHW_M
RAW_NODATA_MAX = -9.0
BEACH_START_THR_M = 0.50
WATER_CLAMP_M = -3.0
SENTINEL_M = -3.0
GRID_M = 10.0
DAM_TO_M = 10.0
DUNE_ROWS = 2                      # Barrier3D DuneWidth

ZOOM_ROWS = 30                     # rows drawn in panel (b), incl. the dune rows

# categories
LAND, OPEN, WET, GAP, DUNE = 0, 1, 2, 3, 4
COLORS = ["#e6dcc8", "#c3d9e6", "#1f6fb4", "#d7191c", "#9c7a3c"]
LABELS = ["interior, subaerial ($z$ > MHW)",
          "beyond the island envelope",
          "interior, $z \\leq$ MHW (measured)",
          "unsurveyed $\\rightarrow$ enters the model as $-$3.0 m water",
          "dune rows (DuneWidth = 2)"]
CMAP = ListedColormap(COLORS)
NORM = BoundaryNorm([-0.5, 0.5, 1.5, 2.5, 3.5, 4.5], CMAP.N)

plt.rcParams.update({
    "font.size": 9.5,
    "axes.linewidth": 0.7,
    "axes.titlesize": 10.5,
    "axes.labelsize": 9.5,
    "xtick.direction": "out",
    "ytick.direction": "out",
    "xtick.major.width": 0.7,
    "ytick.major.width": 0.7,
    "xtick.labelsize": 9,
    "ytick.labelsize": 9,
    "legend.frameon": False,
    "legend.fontsize": 9,
})



# Every output of this folder lands under one directory beside the extraction it
# describes, rather than being scattered through the run folder it did not
# produce. audit_dir() is the only place that name is spelled.
AUDIT_SUBDIR = "nodata-audit"


def audit_dir(topo_dir):
    """<product>/dune-topo/<version>/nodata-audit/, created on demand."""
    d = topo_dir.parent / AUDIT_SUBDIR
    d.mkdir(parents=True, exist_ok=True)
    return d


# =============================================================================
# CLASSIFY
# =============================================================================

def classify_input(raw):
    """Raw ocean-first m NAVD88 -> (unsurveyed, sub-MHW) counts in the envelope.

    The envelope rule is the one HAT_plot_dem_holes.py uses, so panel (d)'s
    'in the DEM' bars mean the same thing that figure's red does.
    """
    nod = raw <= RAW_NODATA_MAX
    z = raw - MHW_M
    land = (~nod) & (z > 0.0)

    zc = np.where(nod, WATER_CLAMP_M, z)
    zc[zc < WATER_CLAMP_M] = WATER_CLAMP_M
    above = zc > BEACH_START_THR_M
    beach_start = np.where(above.any(axis=1), above.argmax(axis=1), -1)
    n_cross = raw.shape[1]
    last_land = np.where(land.any(axis=1),
                         n_cross - 1 - land[:, ::-1].argmax(axis=1), -1)

    gap = wet = 0
    for p in range(raw.shape[0]):
        if beach_start[p] < 0 or last_land[p] < beach_start[p]:
            continue
        s = slice(beach_start[p], last_land[p] + 1)
        gap += int(nod[p, s].sum())
        wet += int(((~nod[p, s]) & (z[p, s] <= 0)).sum())
    return gap, wet


def stack_domain(topo_m, nodata):
    """The CASCADE input domain as one categorical array, ocean-first.

    Rows 0..DUNE_ROWS-1 are the dune, then the interior. Returns the codes and
    the per-profile last-land row of the INTERIOR part, in interior numbering.
    """
    land = topo_m > 0.0
    last_land = np.where(land.any(axis=0),
                         topo_m.shape[0] - 1 - land[::-1, :].argmax(axis=0), -1)
    idx = np.arange(topo_m.shape[0])[:, None]
    inside = (idx <= last_land[None, :]) & (last_land[None, :] >= 0)

    interior = np.full(topo_m.shape, OPEN, dtype=np.int8)
    interior[land] = LAND
    interior[inside & ~land & ~nodata] = WET
    interior[inside & nodata] = GAP

    dune = np.full((DUNE_ROWS, topo_m.shape[1]), DUNE, dtype=np.int8)
    return np.vstack([dune, interior]), last_land


def widths_and_bound(topo_m, nodata):
    """(width Barrier3D sees, width if every unsurveyed cell were land), cells.

    The first is roadway.interior_widths verbatim - barrier3d.FindWidths, which
    stops at the first cell at or below sea level. The second re-runs it on a
    copy with the unsurveyed cells lifted above the threshold, which is the
    most land those cells could possibly be hiding.
    """
    seen = roadway.interior_widths(topo_m)
    lifted = np.where(nodata, 1.0, topo_m)
    return seen, roadway.interior_widths(lifted)


def read_origins(dem_product):
    out = {}
    P = elevation_product(dem_product)
    with (P.resampled_10m / "resample_audit.csv").open() as f:
        for r in csv.DictReader(f):
            out[int(r["domain"])] = (float(r["origin_x"]), float(r["origin_y"]))
    return out


def read_setbacks():
    """1984 road setback in metres per GIS domain, or None if unavailable.

    Only used for the t = 0 drowning check, which is a printout. A missing file
    downgrades that check rather than failing the figure.
    """
    p = (REPO / "data" / "hatteras_init" / "4-mgmt-forcing" / "road_offset"
         / "dunestart_offset" / "1984" / "RoadSetback_1984_dunestart.csv")
    if not p.is_file():
        return None
    rows = list(csv.reader(p.open()))
    if len(rows) < 2:
        return None
    ids = [int(float(v)) for v in rows[0]]
    vals = [float(v) for v in rows[1]]
    return dict(zip(ids, vals))


# =============================================================================
# MAIN
# =============================================================================

def main():
    topo_dir, dune_dir, version = htv.topo_dirs(TOPO_PRODUCT, VERSION_OVERRIDE)
    arr_dir, _ = htv.npy_dirs(TOPO_PRODUCT)
    out_png = audit_dir(topo_dir) / f"HAT_topo_retained_{version}.png"
    print(f"product {TOPO_PRODUCT}, version {version}")

    origins = read_origins(DEM_PRODUCT)
    domains = sorted(origins)
    setbacks = read_setbacks()

    stacks, stats = {}, []
    dune_sentinel = 0
    for n in domains:
        topo = np.load(topo_dir / htv.array_name("topography", n)) * DAM_TO_M
        nod = np.load(topo_dir / htv.array_name("nodata", n))
        dune = np.load(dune_dir / htv.array_name("dune", n)) * DAM_TO_M
        dune_sentinel += int((dune <= SENTINEL_M + 1e-9).sum())

        stacks[n], last_land = stack_domain(topo, nod)
        seen, bound = widths_and_bound(topo, nod)

        raw = np.load(arr_dir / f"domain_{n}.npy").astype(float)[::-1, ::-1]
        in_gap, in_wet = classify_input(raw)

        idx = np.arange(topo.shape[0])[:, None]
        inside = (idx <= last_land[None, :]) & (last_land[None, :] >= 0)
        out_gap = int((inside & nod).sum())
        out_wet = int((inside & ~nod & (topo <= 0.0)).sum())

        drown = None
        if setbacks and n in setbacks:
            drown = roadway.predict_drowning(topo, setbacks[n])

        stats.append(dict(domain=n, in_gap=in_gap, in_wet=in_wet,
                          out_gap=out_gap, out_wet=out_wet,
                          seen=seen, bound=bound, rows=topo.shape[0],
                          drown=drown))

    n_along = stacks[domains[0]].shape[1]

    # --- shared alongshore axis ---------------------------------------------
    oy_all = np.array([origins[n][1] for n in domains])
    north_max = oy_all.max()
    nx = int(round((north_max - (oy_all.min() - n_along * GRID_M)) / GRID_M))

    def col0(oy):
        return nx - int(round((north_max - oy) / GRID_M)) - n_along

    XLIM = (0, nx * GRID_M / 1000.0)
    dom_km = np.array([(col0(oy) + n_along / 2.0) * GRID_M / 1000.0
                       for oy in oy_all])
    ticks = list(range(0, len(domains), 10))

    max_rows = max(s.shape[0] for s in stacks.values())
    grid = np.full((max_rows, nx), -1, dtype=np.int8)
    for i, n in enumerate(domains):
        s = stacks[n]
        x0 = col0(oy_all[i])
        grid[:s.shape[0], x0:x0 + n_along] = s
    gm = np.ma.masked_where(grid < 0, grid)
    zoom = np.ma.masked_where(grid[:ZOOM_ROWS] < 0, grid[:ZOOM_ROWS])

    # --- totals --------------------------------------------------------------
    T = {k: sum(s[k] for s in stats)
         for k in ("in_gap", "in_wet", "out_gap", "out_wet")}
    seen_all = np.concatenate([s["seen"] for s in stats])
    bound_all = np.concatenate([s["bound"] for s in stats])
    lost_all = bound_all - seen_all
    lost_m = np.array([float((s["bound"] - s["seen"]).mean()) * GRID_M
                       for s in stats])
    seen_m = np.array([float(s["seen"].mean()) * GRID_M for s in stats])

    # =========================================================================
    # DRAW
    # =========================================================================
    fig = plt.figure(figsize=(17.5, 13.8))
    gs = fig.add_gridspec(4, 1, height_ratios=[1.50, 1.05, 0.85, 0.80],
                          hspace=0.42, left=0.070, right=0.986,
                          top=0.840, bottom=0.050)

    def panel_title(ax, letter, text):
        ax.set_title(f"({letter})  {text}", loc="left", fontsize=10.5)

    # --- (a) full stack ------------------------------------------------------
    axA = fig.add_subplot(gs[0])
    axA.imshow(gm, cmap=CMAP, norm=NORM, origin="lower", aspect="auto",
               interpolation="nearest",
               extent=[XLIM[0], XLIM[1], 0, max_rows * GRID_M])
    axA.set_facecolor("white")
    axA.axhline(DUNE_ROWS * GRID_M, color="#111111", lw=1.0, ls=(0, (5, 3)),
                zorder=5)
    panel_title(axA, "a", "The CASCADE input domain: 2 dune rows, then the "
                          "interior. Ocean at the bottom, as Barrier3D indexes it")
    axA.set_ylabel("distance landward of\nthe dune toe (m)")
    axA.set_xlim(*XLIM)
    axA.tick_params(labelbottom=False)

    axt = axA.secondary_xaxis("top")
    axt.set_xticks(dom_km[ticks])
    axt.set_xticklabels([str(domains[t]) for t in ticks])
    axt.set_xlabel("Barrier3D domain", labelpad=3)

    # --- (b) zoom ------------------------------------------------------------
    axB = fig.add_subplot(gs[1], sharex=axA)
    axB.imshow(zoom, cmap=CMAP, norm=NORM, origin="lower", aspect="auto",
               interpolation="nearest",
               extent=[XLIM[0], XLIM[1], 0, ZOOM_ROWS * GRID_M])
    axB.set_facecolor("white")
    axB.axhline(DUNE_ROWS * GRID_M, color="#111111", lw=1.3, ls=(0, (5, 3)),
                zorder=5)
    axB.annotate("dune rows 0-1", xy=(XLIM[1] * 0.004, DUNE_ROWS * GRID_M * 2.1),
                 fontsize=9.5, color="#111111", zorder=6, va="center")
    panel_title(axB, "b", f"The seaward {ZOOM_ROWS * GRID_M:.0f} m of the same "
                          f"stack: at island scale the dune is one pixel")
    axB.set_ylabel("distance landward of\nthe dune toe (m)")
    axB.tick_params(labelbottom=False)

    # --- (c) widths ----------------------------------------------------------
    axC = fig.add_subplot(gs[2], sharex=axA)
    w = n_along * GRID_M / 1000.0 * 0.92
    axC.bar(dom_km, seen_m, width=w, color="#e6dcc8", edgecolor="0.55", lw=0.4,
            label="island width Barrier3D sees (FindWidths)")
    axC.bar(dom_km, lost_m, width=w, bottom=seen_m, color=COLORS[GAP],
            alpha=0.75, lw=0,
            label="upper bound on width lost to unsurveyed cells")
    axC.set_ylabel("mean island width\nper domain (m)")
    axC.grid(axis="y", lw=0.4, color="0.88")
    axC.set_axisbelow(True)
    axC.legend(loc="upper right", ncol=2, frameon=True, framealpha=0.95,
               edgecolor="0.8", fancybox=False)
    axC.set_ylim(0, (seen_m + lost_m).max() * 1.32)
    panel_title(axC, "c", "What the truncation costs: FindWidths stops at the "
                          f"first water cell, so {int((lost_all > 0).sum()):,} of "
                          f"{lost_all.size:,} profiles lose land behind an "
                          f"unsurveyed cell")
    axC.tick_params(labelbottom=False)

    # --- (d) retention -------------------------------------------------------
    axD = fig.add_subplot(gs[3], sharex=axA)
    a_in = np.array([s["in_gap"] for s in stats], float)
    a_out = np.array([s["out_gap"] for s in stats], float)
    axD.bar(dom_km, a_in, width=w, facecolor="none", edgecolor=COLORS[GAP],
            lw=0.9, label="unsurveyed in the DEM, within the envelope")
    axD.bar(dom_km, a_out, width=w * 0.48, color=COLORS[GAP], alpha=0.65, lw=0,
            label="still there in the input domain")
    axD.set_ylabel("unsurveyed\n(cells)")
    axD.grid(axis="y", lw=0.4, color="0.88")
    axD.set_axisbelow(True)
    axD.legend(loc="upper right", ncol=2, frameon=True, framealpha=0.95,
               edgecolor="0.8", fancybox=False)
    axD.set_ylim(0, a_in.max() * 1.34)
    panel_title(axD, "d", f"Unsurveyed cells reaching the model: "
                          f"{T['out_gap']:,} of {T['in_gap']:,} in the DEM "
                          f"({100 * T['out_gap'] / T['in_gap']:.1f}%). None is "
                          f"flagged - all enter as $-$3.0 m water")
    axD.set_xlabel("alongshore distance (km of UTM northing from the southern "
                   "edge of domain 1; south $\\rightarrow$ north)")

    for ax in (axA, axB, axC, axD):
        for s in ax.spines.values():
            s.set_color("0.35")

    fig.text(0.070, 0.982,
             f"Unsurveyed ground in the CASCADE input domains: "
             f"{TOPO_PRODUCT}, dune-topo {version}",
             fontsize=13, fontweight="bold", va="top", ha="left")
    fig.text(0.070, 0.958,
             "Barrier3D has no representation for \"unknown\", so every "
             "unsurveyed cell is written to the water sentinel and read as an "
             "elevation of exactly $-$3.0 m MHW. The nodata mask is a sidecar "
             "the run never opens.",
             fontsize=9.5, color="0.30", va="top", ha="left")

    handles = [Patch(facecolor=c, edgecolor="0.45", lw=0.5, label=l)
               for c, l in zip(COLORS, LABELS)]
    fig.legend(handles=handles, loc="upper left", ncol=3, frameon=False,
               fontsize=10, bbox_to_anchor=(0.070, 0.933))

    out_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_png, dpi=170, facecolor="white")
    print(f"wrote {out_png}\n")

    # --- printout ------------------------------------------------------------
    print("at model start")
    print("  cells flagged no-data in what CASCADE reads   : 0 "
          "(no such state exists)")
    print(f"  unsurveyed cells entering as -3.0 m water     : {T['out_gap']:,}"
          f"  of {T['in_gap']:,} in the DEM envelope")
    print(f"  dune cells resting on unsurveyed ground       : {dune_sentinel} "
          f"of {len(domains) * n_along}")
    print()
    print("island width Barrier3D sees (FindWidths, stops at first water cell)")
    print(f"  profiles truncated by an unsurveyed cell      : "
          f"{int((lost_all > 0).sum()):,} of {lost_all.size:,}")
    if (lost_all > 0).any():
        print(f"  land hidden behind one, upper bound           : "
              f"median {np.median(lost_all[lost_all > 0]) * GRID_M:.0f} m, "
              f"max {lost_all.max() * GRID_M:.0f} m")
    print("  worst domains, mean width lost per profile:")
    for i in np.argsort(-lost_m)[:6]:
        print(f"    domain {domains[i]:>3}  {lost_m[i]:>6.0f} m  "
              f"(sees {seen_m[i]:.0f} m)")

    if setbacks:
        dr = [s for s in stats if s["drown"] and s["drown"].get("drowns")]
        wall = [s for s in stats if s["drown"] and s["drown"].get("wall")]
        print()
        print("roadway drowning test at t = 0 (roadway_manager.bulldoze)")
        print(f"  domains drowning at t = 0                    : {len(dr)}"
              + (f"  -> {[s['domain'] for s in dr]}" if dr else ""))
        if wall:
            print(f"  setback unusable                             : "
                  f"{[(s['domain'], s['drown']['wall']) for s in wall]}")
    else:
        print("\n  [skip] no 1984 setback CSV found; drowning test not run")


if __name__ == "__main__":
    main()
