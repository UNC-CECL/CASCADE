"""
HAT_plot_dem_holes.py

Where the nodata and sub-MHW holes are in a Barrier3D start DEM.

WHAT THIS ANSWERS
-----------------
"Do I have a lot of cells by my dunes that are under MHW, or data gaps?"
Counting them per domain (the audit CSVs) says how many. This says WHERE, and
separates the two things that both look like "low" on an elevation ramp:

    a cell that is water because it is the ocean or the sound   - expected
    a cell that is water INSIDE the island                      - a hole

Only the second is coloured loudly. Everything else is deliberately pale, so a
figure of a clean DEM is a quiet figure.

THE CATEGORIES, per profile, in the extractor's own frame
---------------------------------------------------------
Profiles are read ocean-first (arr[:, ::-1], OCEAN_LOC="right"), z = raw - MHW,
exactly as HAT_dune_topo_extractor.load_domain does.

    beach_start  first cell with z > 0.50 m       (BEACH_START_THR_M)
    last_land    last cell with z > 0

Cells seaward of beach_start or landward of last_land are open water. Cells
BETWEEN them that are not land are holes:

    land          z > 0                             pale sand
    open water    outside [beach_start, last_land]  pale blue
    wet hole      z <= 0, valid data, inside        strong blue
    gap           nodata (raw <= -9), inside        strong red

Nodata is never a step on an elevation ramp here, for the reason FIGURES.md
gives: "not surveyed" is not a low elevation, and conflating the two is what
drowned three roadways at t=0.

The two affected categories are #1f6fb4 blue against #d7191c red. The red is
saturated on purpose. Unsurveyed cells are the category a reader must not skim
past - they are the ones that become a fictitious elevation downstream - and
they are also the rarer of the two at 0.32% of cells, so at this scale they
have to hold their own against the blue at single-pixel widths. A muted
#d6604d, sampled off the same RdBu ramp as the blue, was tried and is the
better choice for a figure meant to sit quietly in a page of body text; it lost
too much at one-cell width here.

Under a red-green deficiency the red darkens towards olive while the blue
holds, so the pair still parts. Magenta against this blue, the first draft, is
the pair to avoid: it holds neither the hue nor the luminance gap.

ORIENTATION
-----------
Ocean is at the BOTTOM of panels A and B, landward upward. That is the
convention HAT_dune_topo_extractor.pick_window draws for picking a dune search
window, so this figure and the picker read the same way round.

THE THREE PANELS
----------------
A  the island unrolled. Alongshore runs left-right; cross-shore is UTM easting
   with the island trend removed by a cubic fit through the 90 domain origins,
   so a 45 km arc lies flat instead of drifting 5 km across the panel. The
   detrend is a rigid per-domain shift - no cell is resampled, and cross-shore
   distances within a domain are untouched. This is the locator panel: it keeps
   each domain's own shoreline shape, which panel B removes. The cross-shore
   axis is metres landward of the seaward edge of the detrended strip, so its
   zero is a drawing origin and not a landform.

B  the same cells straightened: every profile shifted so its own beach_start
   sits at cross-shore 0. The dune band becomes a horizontal stripe instead of
   following the shoreline curve, so a hole IN THE DUNES is separable from a
   hole 500 m behind them. Panel A cannot show that; the shoreline moves.

C  per-domain percentages, dune band vs interior.

ONE ALONGSHORE AXIS
-------------------
All three panels share x, in kilometres of UTM northing measured from the
southern edge of domain 1. Every domain is painted at its true northing, so a
vertical line means the same place in all three panels. Domains are spaced
~504 m and carry 500 m of data, so there are ~4 m unpainted seams between them;
that is real, not a plotting artefact.

Left to right is south -> north, which is the model alongshore direction
(ALONGSHORE_FLIP = True flips the raster north-at-top rows so profile index
increases northward). Within a domain, profile p is raster row 49 - p.

INPUT   data/hatteras_init/1-barrier3d-domains/<PRODUCT>/npy-arrays/domain_<N>.npy
        the arrays CASCADE reads - m NAVD88, -10 nodata - not the .tif, so what
        is drawn is what the model ingests.
        Georeferencing comes from the elevation product resample_audit.csv.

OUTPUT  <elevation product>/figures/HAT_<slug>_holes.png
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

REPO = next(_p for _p in Path(__file__).resolve().parents
            if (_p / "pyproject.toml").exists())
sys.path.insert(0, str(REPO / "scripts"))
from hat_elevation_products import product  # noqa: E402

# =============================================================================
# CONFIG
# =============================================================================

TOPO_PRODUCT = "1984-start"        # which npy-arrays folder
DEM_PRODUCT = "2009-2014-1996"     # which elevation product it came from
SLUG = "1984dem"

MHW_M = 0.36               # m NAVD88
RAW_NODATA_MAX = -9.0      # raw <= this is nodata; raw nodata is exactly -10
BEACH_START_THR_M = 0.50   # m MHW, strict '>'
WATER_CLAMP_M = -3.0       # m MHW
DUNE_WIN_CELLS = 8         # DEFAULT_WINDOW_PX, the extractor default window
GRID_M = 10.0

STRAIGHT_ROWS = 120        # cross-shore cells drawn in panel B, from beach start
DETREND_DEG = 3            # polynomial in northing removed from easting, panel A
A_PAD_CELLS = 6            # blank cells kept above/below the unrolled strip

# categories
LAND, OPEN, WET, GAP = 0, 1, 2, 3
COLORS = ["#e6dcc8", "#c3d9e6", "#1f6fb4", "#d7191c"]
LABELS = ["subaerial ($z$ > MHW)",
          "open water, outside the island envelope",
          "$z \\leq$ MHW, within the envelope",
          "unsurveyed, within the envelope"]
CMAP = ListedColormap(COLORS)
NORM = BoundaryNorm([-0.5, 0.5, 1.5, 2.5, 3.5], CMAP.N)

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

INIT = REPO / "data" / "hatteras_init"
ARR_DIR = INIT / "1-barrier3d-domains" / TOPO_PRODUCT / "1-extraction" / "npy-arrays"
P = product(DEM_PRODUCT)
OUT_PNG = P.figures / f"HAT_{SLUG}_holes.png"


# =============================================================================
# CLASSIFY
# =============================================================================

def classify(raw):
    """(n_along, n_cross) ocean-first raw m NAVD88 -> codes, beach_start, last_land.

    beach_start / last_land are -1 on a profile with no land at all. None exist
    in this product, but a forecast domain could have one.
    """
    nod = raw <= RAW_NODATA_MAX
    z = raw - MHW_M
    land = (~nod) & (z > 0.0)

    zc = np.where(nod, WATER_CLAMP_M, z)
    zc[zc < WATER_CLAMP_M] = WATER_CLAMP_M
    above = zc > BEACH_START_THR_M
    beach_start = np.where(above.any(axis=1), above.argmax(axis=1), -1)

    n_cross = raw.shape[1]
    rev = land[:, ::-1]
    last_land = np.where(land.any(axis=1), n_cross - 1 - rev.argmax(axis=1), -1)

    idx = np.arange(n_cross)[None, :]
    inside = ((idx >= beach_start[:, None]) & (idx <= last_land[:, None])
              & (beach_start[:, None] >= 0))

    cat = np.full(raw.shape, OPEN, dtype=np.int8)
    cat[land] = LAND
    cat[inside & ~land & ~nod] = WET
    cat[inside & nod] = GAP
    return cat, beach_start, last_land


def load_domain(n):
    """Raw array, ocean-first cross-shore, profile index increasing NORTHWARD.

    The .npy is raster order: row 0 north, column 199 east = ocean. [:, ::-1]
    puts the ocean first, which is what the extractor does. [::-1] on the rows
    is ALONGSHORE_FLIP, so profile 0 is the southern edge of the domain.
    """
    a = np.load(ARR_DIR / f"domain_{n}.npy").astype(float)
    return a[::-1, ::-1]


def read_origins():
    """domain -> (origin_x, origin_y) from the product resample audit."""
    out = {}
    with (P.resampled_10m / "resample_audit.csv").open() as f:
        for r in csv.DictReader(f):
            out[int(r["domain"])] = (float(r["origin_x"]), float(r["origin_y"]))
    return out


# =============================================================================
# BUILD
# =============================================================================

def main():
    origins = read_origins()
    domains = sorted(origins)
    print(f"{len(domains)} domains from {ARR_DIR}")

    cats, bstarts, llands = {}, {}, {}
    for n in domains:
        c, b, l = classify(load_domain(n))
        cats[n], bstarts[n], llands[n] = c, b, l

    n_along, n_cross = cats[domains[0]].shape

    # --- the shared alongshore axis -----------------------------------------
    # One 10 m grid in UTM northing. Domain n occupies columns [x0, x0+50),
    # counting from the SOUTH, so column index increases northward like the
    # model profile index does.
    ox_all = np.array([origins[n][0] for n in domains])
    oy_all = np.array([origins[n][1] for n in domains])
    north_max, north_min = oy_all.max(), oy_all.min() - n_along * GRID_M
    nx = int(round((north_max - north_min) / GRID_M))
    x_km = np.arange(nx + 1) * GRID_M / 1000.0

    def col0(oy):
        """Southernmost column of the domain whose north edge is at oy."""
        return nx - int(round((north_max - oy) / GRID_M)) - n_along

    # --- panel A: unroll the island -----------------------------------------
    # Remove the island trend from easting so the strip lies flat. The fit is
    # evaluated once per domain, so each block is shifted rigidly - no cell is
    # resampled and no cross-shore distance changes.
    coef = np.polyfit(oy_all, ox_all, DETREND_DEG)
    resid = ox_all - np.polyval(coef, oy_all)          # west edge, detrended
    j_all = np.round((resid - resid.min()) / GRID_M).astype(int)
    ny = int(j_all.max() + n_cross + 2 * A_PAD_CELLS)
    grid = np.full((ny, nx), -1, dtype=np.int8)        # -1 = no domain here

    for i, n in enumerate(domains):
        # cats[n] is (profile S->N, cross ocean-first). Flip cross so index
        # increases eastward, then transpose to (easting, alongshore).
        block = cats[n][:, ::-1].T
        j0 = int(j_all[i]) + A_PAD_CELLS
        x0 = col0(oy_all[i])
        grid[j0:j0 + n_cross, x0:x0 + n_along] = block

    # Row 0 is the WEST (sound) edge. Flip so row 0 is the ocean edge, then
    # origin="lower" puts the ocean at the bottom with landward upward, the
    # same way round as the extractor's pick_window.
    grid = grid[::-1]
    gm = np.ma.masked_where(grid < 0, grid)
    a_km = ny * GRID_M / 1000.0

    # --- panel B: straighten every profile ----------------------------------
    strt = np.full((STRAIGHT_ROWS, nx), -1, dtype=np.int8)
    for i, n in enumerate(domains):
        c, b = cats[n], bstarts[n]
        x0 = col0(oy_all[i])
        for p in range(n_along):
            if b[p] < 0:
                continue
            seg = c[p, b[p]:b[p] + STRAIGHT_ROWS]
            strt[:seg.size, x0 + p] = seg
    sm = np.ma.masked_where(strt < 0, strt)

    # --- panel C: per-domain percentages ------------------------------------
    dune_gap, dune_wet, int_gap, int_wet = [], [], [], []
    for n in domains:
        c, b, l = cats[n], bstarts[n], llands[n]
        dg = dw = dt = ig = iw = it = 0
        for p in range(n_along):
            if b[p] < 0:
                continue
            d = c[p, b[p]:b[p] + DUNE_WIN_CELLS]
            dt += d.size
            dg += int((d == GAP).sum())
            dw += int((d == WET).sum())
            i = c[p, b[p] + DUNE_WIN_CELLS:l[p] + 1]
            it += i.size
            ig += int((i == GAP).sum())
            iw += int((i == WET).sum())
        dune_gap.append(100 * dg / max(dt, 1))
        dune_wet.append(100 * dw / max(dt, 1))
        int_gap.append(100 * ig / max(it, 1))
        int_wet.append(100 * iw / max(it, 1))

    dune_tot = np.array(dune_wet) + np.array(dune_gap)
    int_tot = np.array(int_wet) + np.array(int_gap)

    # island-wide dune-band figure, for the panel B annotation
    dband = np.concatenate([cats[n][p, bstarts[n][p]:bstarts[n][p] + DUNE_WIN_CELLS]
                            for n in domains for p in range(n_along)
                            if bstarts[n][p] >= 0])
    dune_pct_island = 100 * np.isin(dband, (WET, GAP)).mean()

    tot = sum(cats[n].size for n in domains)
    ngap = sum(int((cats[n] == GAP).sum()) for n in domains)
    nwet = sum(int((cats[n] == WET).sum()) for n in domains)
    ngap_pct, nwet_pct = 100 * ngap / tot, 100 * nwet / tot

    # =========================================================================
    # DRAW
    # =========================================================================
    fig = plt.figure(figsize=(17.5, 11.6))
    gs = fig.add_gridspec(3, 1, height_ratios=[1.05, 1.85, 1.05], hspace=0.33,
                          left=0.070, right=0.986, top=0.850, bottom=0.070)

    dom_km = np.array([(col0(oy) + n_along / 2.0) * GRID_M / 1000.0
                       for oy in oy_all])
    ticks = list(range(0, len(domains), 10))
    XLIM = (0, nx * GRID_M / 1000.0)

    def panel_title(ax, letter, text):
        ax.set_title(f"({letter})  {text}", loc="left", fontsize=10.5)
        ax.title.set_position((0.0, 1.0))

    # --- a -------------------------------------------------------------------
    axA = fig.add_subplot(gs[0])
    axA.imshow(gm, cmap=CMAP, norm=NORM, origin="lower", aspect="auto",
               interpolation="nearest", extent=[XLIM[0], XLIM[1], 0, a_km])
    axA.set_facecolor("white")
    panel_title(axA, "a", "Domain mosaic, detrended in UTM easting; "
                          "ocean at the base of each profile")
    axA.set_ylabel("cross-shore (km)\ndetrended easting")
    axA.set_xlim(*XLIM)
    axA.tick_params(labelbottom=False)

    axt = axA.secondary_xaxis("top")
    axt.set_xticks(dom_km[ticks])
    axt.set_xticklabels([str(domains[t]) for t in ticks])
    axt.set_xlabel("Barrier3D domain", labelpad=3)

    # --- b -------------------------------------------------------------------
    axB = fig.add_subplot(gs[1], sharex=axA)
    axB.imshow(sm, cmap=CMAP, norm=NORM, origin="lower", aspect="auto",
               interpolation="nearest",
               extent=[XLIM[0], XLIM[1], 0, STRAIGHT_ROWS * GRID_M])
    axB.set_facecolor("white")
    for yv in (0, DUNE_WIN_CELLS * GRID_M):
        axB.axhline(yv, color="#111111", lw=1.1, ls=(0, (5, 3)), zorder=5)
    axB.annotate(f"dune search band, 0-{DUNE_WIN_CELLS * GRID_M:.0f} m: "
                 f"{dune_pct_island:.2f}% of cells $\\leq$ MHW or unsurveyed",
                 xy=(XLIM[1] * 0.004, DUNE_WIN_CELLS * GRID_M * 1.7),
                 fontsize=9.5, color="#111111", zorder=6, va="center")
    panel_title(axB, "b", "Profiles aligned on their own beach start "
                          "($z$ > 0.50 m first crossing), shore-normal")
    axB.set_ylabel("distance landward of\nbeach start (m)")
    axB.tick_params(labelbottom=False)

    # --- c -------------------------------------------------------------------
    axC = fig.add_subplot(gs[2], sharex=axA)
    w = n_along * GRID_M / 1000.0 * 0.88
    axC.bar(dom_km, int_wet, width=w, color=COLORS[WET], alpha=0.55, lw=0,
            label="interior: $z \\leq$ MHW")
    axC.bar(dom_km, int_gap, width=w, bottom=int_wet, color=COLORS[GAP],
            alpha=0.55, lw=0, label="interior: unsurveyed")
    axC.plot(dom_km, dune_tot, color="#111111", lw=1.2, marker="o", ms=3.0,
             label="dune search band, both categories")
    axC.set_ylabel("affected cells (%)")
    axC.set_xlabel("alongshore distance (km of UTM northing from the southern "
                   "edge of domain 1; south $\\rightarrow$ north)")
    axC.set_ylim(0, int_tot.max() * 1.30)
    axC.grid(axis="y", lw=0.4, color="0.88")
    axC.set_axisbelow(True)
    axC.legend(loc="upper right", ncol=3, frameon=True, framealpha=0.95,
               edgecolor="0.8", fancybox=False)
    panel_title(axC, "c", "Affected cell fraction per domain, "
                          "dune search band and interior compared")

    for i in np.argsort(-int_tot)[:4]:
        axC.annotate(f"{domains[i]}", (dom_km[i], int_tot[i]),
                     textcoords="offset points", xytext=(0, 5), ha="center",
                     fontsize=9, color="#a50f15")

    axt2 = axC.secondary_xaxis("top")
    axt2.set_xticks(dom_km[ticks])
    axt2.set_xticklabels([str(domains[t]) for t in ticks])
    axt2.set_xlabel("Barrier3D domain", labelpad=3)

    for ax in (axA, axB, axC):
        for s in ax.spines.values():
            s.set_color("0.35")

    # --- key and caption -----------------------------------------------------
    handles = [Patch(facecolor=c, edgecolor="0.45", lw=0.5, label=l)
               for c, l in zip(COLORS, LABELS)]
    fig.legend(handles=handles, loc="upper center", ncol=4, frameon=False,
               fontsize=10, bbox_to_anchor=(0.5, 0.935))
    fig.text(0.070, 0.975,
             f"Unsurveyed and sub-MHW cells in the {TOPO_PRODUCT} "
             f"initial-condition DEM ({DEM_PRODUCT})",
             fontsize=13, fontweight="bold", va="top", ha="left")
    fig.text(0.070, 0.949,
             f"{len(domains)} Barrier3D domains, {n_along} alongshore profiles "
             f"each ($n$ = {len(domains) * n_along}); 10 m grid; "
             f"MHW = {MHW_M:.2f} m NAVD88",
             fontsize=10, color="0.30", va="top", ha="left")

    OUT_PNG.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT_PNG, dpi=170, facecolor="white")
    print(f"wrote {OUT_PNG}")

    print(f"  interior no-data   {ngap:>9,}  ({ngap_pct:.2f}% of all cells)")
    print(f"  interior <= MHW    {nwet:>9,}  ({nwet_pct:.2f}%)")
    print(f"  dune band, island  {dune_pct_island:.2f}%")
    print(f"  dune band worst    {dune_tot.max():.2f}%  "
          f"(domain {domains[int(np.argmax(dune_tot))]})")
    print(f"  interior worst     {int_tot.max():.2f}%  "
          f"(domain {domains[int(np.argmax(int_tot))]})")
    print(f"  detrend residual   {resid.max() - resid.min():.0f} m "
          f"(raw easting spread {ox_all.max() - ox_all.min():.0f} m)")


if __name__ == "__main__":
    main()
