"""
HAT_road_on_island.py   (v2 - fixed canvas crop)
===============================================================================
Draw NC-12 on the Hatteras Island CASCADE initialization, using the same canvas
construction as plot_initialization_poster.py (dune offsets applied, real
alongshore planform), and check the prescribed road setbacks against the road's
own signature in the DEM.

FIXED IN v2
-----------
v1 cropped the canvas at a hard-coded CROSS_SHORE_TRIM_M = 900 m. The dune
offsets place each domain at origin = offset_cells, and along Hatteras those
offsets span far more than 900 m (the island curves hard from Cape Point to
Rodanthe). So domains 1-66 sat past row 90 and were sliced off before anything
was drawn -- they rendered as uniform set_bad() light blue, and the GIS 9-15
zoom panel came back empty because its rows had already been deleted.

v2 auto-crops to the actual data extent and shifts `offs` into the cropped
frame so every overlay follows. Zoom panel limits are derived from where land
actually is in that window rather than a fixed 45-row guess, and each domain's
own dune line (row 0 of its grid = InteriorDomain[0], the line the setback is
measured from) is drawn in the zoom panels.

TWO QUESTIONS, TWO FIGURES
--------------------------
FIG 1  HAT_road_on_island_<YEAR>.png
       Where does the setback put the road, relative to the actual island?

FIG 2  HAT_road_profile_check_<YEAR>.png
       Is the setback pointing at the road corridor? NC-12 is a graded,
       shore-parallel surface, so alongshore averaging PRESERVES it -- it
       should appear as a flat shelf near ROAD_ELEVATION. If the setback band
       does not sit on that shelf, the setback is measured from the wrong line.

WHAT THIS SCRIPT DOES *NOT* ANSWER
----------------------------------
Whether the road fits in CASCADE's InteriorDomain at runtime. These .npy files
are the full raw extraction (barrier + sound); Barrier3D trims them. Use
HAT_road_placement_check.py against a live Cascade object for that.

REFERENCE FRAME (the thing being tested)
----------------------------------------
Barrier3D measures road_setback in METRES landward from InteriorDomain[0], i.e.
row 0 of the topography array. Not MHW. Not the shoreline. Not the dune toe.

Note the unit split in the run script, which is correct but easy to break:
    shoreline_offset -> DECAMETRES (dune_offset_all / 10.0)
    road_setback     -> METRES     (raw)
===============================================================================
"""

import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.colors import FuncNorm

# =============================================================================
# CONFIGURATION
# =============================================================================

PROJECT_BASE_DIR   = r"/"
HATTERAS_DATA_BASE = os.path.join(PROJECT_BASE_DIR, "data", "hatteras_init")
OUTPUT_DIR         = os.path.join(PROJECT_BASE_DIR, "scripts", "input_prep",
                                  "road_offset", "visualization", "comparison")

# --- Which period ------------------------------------------------------------
CHECK_YEAR = 1984

DUNE_OFFSET_FILE = os.path.join(
    HATTERAS_DATA_BASE, "2-brie-offset", f"hindcast_{CHECK_YEAR}",
    "Island_Dune_Offsets_1984_PADDED_120.csv",
)
DUNE_OFFSET_COLUMN = 0     # 0 if single-column; set if the file has year columns

ROAD_SETBACK_FILE = os.path.join(
    HATTERAS_DATA_BASE, "roads", "processed_offset", str(CHECK_YEAR),
    f"RoadSetback_{CHECK_YEAR}.csv",
)

# --- Topography --------------------------------------------------------------
TOPO_DUNE_INIT_YEAR = "2009"
TOPO_DUNE_SUBFOLDER = os.path.join("2009", "2009")

# --- Domain geometry ---------------------------------------------------------
NUM_REAL_DOMAINS   = 90
NUM_BUFFER_DOMAINS = 15
TOTAL_DOMAINS      = NUM_BUFFER_DOMAINS + NUM_REAL_DOMAINS + NUM_BUFFER_DOMAINS
START_REAL_INDEX   = NUM_BUFFER_DOMAINS
END_REAL_INDEX     = START_REAL_INDEX + NUM_REAL_DOMAINS
FIRST_FILE_NUMBER  = 1
FIRST_ROAD_DOMAIN  = 9
LAST_ROAD_DOMAIN   = 90

DAM_TO_M = 10.0
DX = DY = 10.0

FLIPLR_DOMAIN = True    # matches plot_initialization_poster.py

# --- Road parameters (from main() in the run script) -------------------------
ROAD_ELEVATION = 1.45    # m NAVD88
ROAD_WIDTH     = 20.0    # m

# --- Road-shelf detection ----------------------------------------------------
SHELF_TOL_M       = 0.25   # |elev - ROAD_ELEVATION| within this = candidate
SHELF_FLAT_M      = 0.12   # max |d(elev)/d(row)| for a row to count as "flat"
SHELF_ALONG_STD_M = 0.45   # max alongshore std for a row to count as "graded"

# --- Prescribed relocations (from HISTORICAL_ROAD_EVENTS) --------------------
HISTORICAL_RELOCATIONS = {
    1989: {84: 163.0, 85: 165.0, 86: 205.0, 87: 113.0},
    1999: {9: 73.0, 10: 97.0, 11: 129.0, 12: 126.0, 13: 125.0,
           14: 126.0, 15: 106.0},
}
RELOC_STYLE = {1989: dict(color="#1565C0", marker="v"),
               1999: dict(color="#7F77DD", marker="^")}

PROFILE_DOMAINS = [9, 11, 12, 13, 14, 15, 16, 84]

# --- Display -----------------------------------------------------------------
ELEV_MIN_M    = -1.0
ELEV_MAX_M    =  4.0
SEA_LEVEL     =  0.0
SEA_LEVEL_POS =  0.35

# Cross-shore crop of the composite canvas.
#   None  -> auto-crop to the rows that actually contain data (recommended;
#            the dune offsets span kilometres, so any fixed guess will slice
#            the south end of the island off - that was the v1 bug).
#   float -> hard crop at this many metres from the canvas top.
CROSS_SHORE_TRIM_M = None
CROSS_SHORE_PAD_M  = 60      # breathing room around the auto-crop (m)

# Orientation of the cross-shore axis.
#   True  -> dune line / ocean side at the BOTTOM, sound at the top.
#            The alongshore axis already runs GIS 1 -> 90 = south -> north,
#            i.e. north to the right. That is a north-up map rotated 90 deg
#            clockwise, which puts east -- the ocean side -- at the bottom.
#            This is the orientation that reads like the real island.
#   False -> ocean at the top (row 0 at the top; the v2 default).
OCEAN_AT_BOTTOM = True

ZOOM_WINDOWS = [(9, 16, "1999 relocation + GIS 16 (crash suspect)"),
                (84, 90, "1989 relocation - Pea Island")]

# --- Style (HAT_hindcast_1984_2024.py palette) -------------------------------
C_MODEL  = "#FF8C00"
C_TOWN   = "#90AFC5"
C_WIMBLE = "#E0A800"
C_PIER   = "#1565C0"
C_GROIN  = "#B71C1C"
C_INK    = "#1a1a2e"

ISLAND_SECTIONS = [
    ("Cape Point",    1,  6,  None),
    ("Buxton",        7,  8,  C_TOWN),
    ("Avon",         21, 31,  C_TOWN),
    ("Wimble Shoals",32, 67,  C_WIMBLE),
    ("Tri-Village",  68, 83,  C_TOWN),
    ("Pea Island",   84, 90,  None),
]


# =============================================================================
# LOAD
# =============================================================================

def load_offsets():
    arr = np.loadtxt(DUNE_OFFSET_FILE, skiprows=1, delimiter=",")
    if arr.ndim == 2:
        arr = arr[:, DUNE_OFFSET_COLUMN]
    return np.round(arr / DAM_TO_M).astype(int)     # m -> cells


def load_setbacks():
    raw = np.loadtxt(ROAD_SETBACK_FILE, delimiter=",")
    if raw.ndim == 2 and raw.shape[0] == 2:
        ids = raw[0].astype(int)
        vals = raw[1]
        print(f"  Setback file is 2-row (IDs + values): "
              f"GIS {ids.min()}-{ids.max()}, {len(ids)} domains")
        expected = np.arange(FIRST_ROAD_DOMAIN, LAST_ROAD_DOMAIN + 1)
        if not np.array_equal(ids, expected):
            print(f"  !! IDs are NOT contiguous {FIRST_ROAD_DOMAIN}-"
                  f"{LAST_ROAD_DOMAIN}; missing "
                  f"{sorted(set(expected) - set(ids))}")
        return {int(g): float(v) for g, v in zip(ids, vals)}

    vals = np.atleast_1d(raw.squeeze())
    gis_ids = list(range(FIRST_ROAD_DOMAIN, LAST_ROAD_DOMAIN + 1))
    if vals.size != len(gis_ids):
        print(f"  !! setback file has {vals.size} values, expected "
              f"{len(gis_ids)}; mapping positionally, alignment NOT guaranteed")
    return {g: float(v) for g, v in zip(gis_ids, vals)}


def load_domain(gis_id):
    path = os.path.join(HATTERAS_DATA_BASE, "topography", TOPO_DUNE_SUBFOLDER,
                        f"domain_{gis_id}_topography_{TOPO_DUNE_INIT_YEAR}.npy")
    if not os.path.exists(path):
        return None
    arr = np.load(path).astype(float)
    if arr.ndim == 1:
        arr = arr.reshape(-1, 1)
    return arr * DAM_TO_M       # dam -> m


def road_rows(setback_m, width_m=ROAD_WIDTH, dy=DY):
    r0 = int(setback_m / dy)
    return r0, r0 + int(width_m / dy)


# =============================================================================
# ROAD-SHELF DETECTION
# =============================================================================

def find_road_shelf(grid):
    """
    Row index of the most road-like feature in a domain's own grid frame.

    NC-12 is graded, flat, and shore-parallel, so in an alongshore-mean profile
    it should be (a) near ROAD_ELEVATION, (b) cross-shore flat, and (c) low
    alongshore variance.
    """
    prof = np.nanmean(grid, axis=1)
    along_std = np.nanstd(grid, axis=1)
    grad = np.abs(np.gradient(prof))

    near = np.abs(prof - ROAD_ELEVATION) <= SHELF_TOL_M
    flat = grad <= SHELF_FLAT_M
    graded = along_std <= SHELF_ALONG_STD_M

    cand = np.where(near & flat & graded)[0]
    if len(cand) == 0:
        alt = np.where(flat)[0]
        if len(alt) == 0:
            return None
        return int(alt[np.argmin(np.abs(prof[alt] - ROAD_ELEVATION))])
    return int(cand[0])     # most seaward candidate


# =============================================================================
# MAIN
# =============================================================================

def main():
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    print("=" * 78)
    print(f"ROAD ON ISLAND  |  setbacks {CHECK_YEAR}  |  topo {TOPO_DUNE_INIT_YEAR}")
    print("=" * 78)

    offsets_all = load_offsets()
    setbacks = load_setbacks()

    gis_ids, grids, offs = [], [], []
    for g in range(FIRST_FILE_NUMBER, FIRST_FILE_NUMBER + NUM_REAL_DOMAINS):
        grid = load_domain(g)
        if grid is None:
            print(f"  WARNING: domain {g} topography not found")
            continue
        gis_ids.append(g)
        grids.append(np.fliplr(grid) if FLIPLR_DOMAIN else grid)
        pad = START_REAL_INDEX + (g - FIRST_FILE_NUMBER)
        offs.append(int(offsets_all[pad]) if pad < len(offsets_all) else 0)

    if not grids:
        raise RuntimeError("No topography loaded - check TOPO_DUNE_SUBFOLDER.")

    print(f"\n  Loaded {len(grids)} domains")
    print(f"  Domain grid rows: {min(g.shape[0] for g in grids)}-"
          f"{max(g.shape[0] for g in grids)}  "
          f"({min(g.shape[0] for g in grids) * DX:.0f}-"
          f"{max(g.shape[0] for g in grids) * DX:.0f} m)")
    print(f"  Dune offsets: {min(offs)}-{max(offs)} cells "
          f"({min(offs) * DX:.0f}-{max(offs) * DX:.0f} m)  "
          f"-> spread {(max(offs) - min(offs)) * DX:.0f} m")
    print(f"    seawardmost: GIS {gis_ids[int(np.argmin(offs))]}  |  "
          f"landwardmost: GIS {gis_ids[int(np.argmax(offs))]}")

    # ---- shelf detection (in each domain's OWN frame) --------------------
    print("\n" + "-" * 78)
    print("ROAD SIGNATURE IN THE DEM vs PRESCRIBED SETBACK")
    print("-" * 78)
    print(f"{'GIS':>4} {'setback':>8} {'sb_row':>7} {'shelf_row':>10} "
          f"{'delta_m':>8} {'shelf_elev':>11}")
    shelf_rows, delta_rows = {}, []
    for g, grid in zip(gis_ids, grids):
        if not (FIRST_ROAD_DOMAIN <= g <= LAST_ROAD_DOMAIN):
            continue
        row = find_road_shelf(grid)
        shelf_rows[g] = row
        sb = setbacks.get(g, np.nan)
        sb_row = int(sb / DY) if not np.isnan(sb) else np.nan
        if row is None:
            print(f"{g:>4} {sb:>8.1f} {sb_row:>7} {'none':>10} "
                  f"{'-':>8} {'-':>11}")
            continue
        d = (sb_row - row) * DY
        delta_rows.append(dict(gis=g, setback_m=sb, sb_row=sb_row,
                               shelf_row=row, delta_m=d))
        elev = float(np.nanmean(grid, axis=1)[row])
        print(f"{g:>4} {sb:>8.1f} {sb_row:>7} {row:>10} {d:>+8.0f} "
              f"{elev:>11.2f}")

    if delta_rows:
        dd = pd.DataFrame(delta_rows)
        dd.to_csv(os.path.join(OUTPUT_DIR,
                               f"HAT_road_shelf_delta_{CHECK_YEAR}.csv"),
                  index=False)
        d = dd.delta_m.values
        print("\n  DELTA = (setback row - DEM road shelf row) x 10 m")
        print(f"    mean {d.mean():+.0f} m   median {np.median(d):+.0f} m   "
              f"std {d.std():.0f} m   range {d.min():+.0f} to {d.max():+.0f} m")
        print()
        if abs(np.median(d)) < 15:
            print("    -> Setbacks land on the road shelf. Reference frame OK.")
        elif d.std() < 25:
            print(f"    -> CONSTANT OFFSET of ~{np.median(d):+.0f} m across "
                  f"domains.")
            print("       Signature of a wrong reference line: every setback is")
            print("       measured from something ~{:.0f} m {} of row 0."
                  .format(abs(np.median(d)),
                          "seaward" if np.median(d) > 0 else "landward"))
        else:
            print("    -> Scattered, not a constant raw_offset. Either the shelf")
            print("       detector is catching dunes/marsh in some domains, or")
            print("       the setbacks are inconsistent domain to domain. Read")
            print("       the profile figure before trusting this number.")

    # ---- canvas ----------------------------------------------------------
    max_rows = max(g.shape[0] for g in grids)
    canvas_rows = max(offs) + max_rows + 5
    total_cols = sum(g.shape[1] for g in grids)
    canvas = np.full((canvas_rows, total_cols), np.nan)

    col_starts, col_widths, cursor = [], [], 0
    for grid, off in zip(grids, offs):
        nr, nc = grid.shape
        col_starts.append(cursor)
        col_widths.append(nc)
        end = min(off + nr, canvas_rows)
        canvas[off:end, cursor:cursor + nc] = grid[:end - off, :]
        cursor += nc
    col_starts = np.array(col_starts)

    # ---- CROP (this was the v1 bug) --------------------------------------
    rows_with_data = np.where(~np.isnan(canvas).all(axis=1))[0]
    if len(rows_with_data) == 0:
        raise RuntimeError("Canvas is entirely NaN - check offsets/topo paths.")

    print("\n" + "-" * 78)
    print("CANVAS")
    print("-" * 78)
    print(f"  Full canvas: {canvas.shape[0]} rows "
          f"({canvas.shape[0] * DX:.0f} m) x {canvas.shape[1]} cols")
    print(f"  Rows containing data: {rows_with_data[0]}-{rows_with_data[-1]} "
          f"({rows_with_data[0] * DX:.0f}-{rows_with_data[-1] * DX:.0f} m)")

    if CROSS_SHORE_TRIM_M is not None:
        r_lo, r_hi = 0, min(int(CROSS_SHORE_TRIM_M / DX), canvas_rows)
        print(f"  Hard crop at CROSS_SHORE_TRIM_M = {CROSS_SHORE_TRIM_M} m "
              f"-> rows {r_lo}-{r_hi}")
        if r_hi <= rows_with_data[0]:
            print("  !! The crop is ENTIRELY ABOVE the data. Everything will")
            print("     render as NaN. Set CROSS_SHORE_TRIM_M = None.")
    else:
        pad = int(CROSS_SHORE_PAD_M / DX)
        r_lo = max(int(rows_with_data[0]) - pad, 0)
        r_hi = min(int(rows_with_data[-1]) + pad + 1, canvas_rows)
        print(f"  Auto-crop to data + {CROSS_SHORE_PAD_M:.0f} m pad "
              f"-> rows {r_lo}-{r_hi} ({(r_hi - r_lo) * DX:.0f} m tall)")

    canvas = canvas[r_lo:r_hi, :]
    offs = [o - r_lo for o in offs]     # shift overlays into the cropped frame

    def domain_centre(g):
        i = gis_ids.index(g)
        return col_starts[i] + col_widths[i] / 2.0

    def domain_span(g):
        i = gis_ids.index(g)
        return col_starts[i], col_starts[i] + col_widths[i]

    # ---- colormap (from the poster script) -------------------------------
    cmap = plt.cm.terrain.copy()
    cmap.set_bad(color="#b0cfe8")

    def _fwd(x):
        r = np.where(x < SEA_LEVEL,
                     SEA_LEVEL_POS * (x - ELEV_MIN_M) / (SEA_LEVEL - ELEV_MIN_M),
                     SEA_LEVEL_POS + (1.0 - SEA_LEVEL_POS) * x / ELEV_MAX_M)
        return np.where(np.isnan(x), np.nan, r)

    def _inv(x):
        return np.where(x < SEA_LEVEL_POS,
                        ELEV_MIN_M + (x / SEA_LEVEL_POS) * (SEA_LEVEL - ELEV_MIN_M),
                        (x - SEA_LEVEL_POS) / (1.0 - SEA_LEVEL_POS) * ELEV_MAX_M)

    norm = FuncNorm((_fwd, _inv), vmin=ELEV_MIN_M, vmax=ELEV_MAX_M)

    # =====================================================================
    # FIGURE 1: island planform + road
    # =====================================================================
    plt.rcParams.update({
        "font.family": "DejaVu Sans", "font.size": 10,
        "figure.facecolor": "white", "axes.facecolor": "white",
        "text.color": C_INK, "axes.labelcolor": C_INK,
        "xtick.color": C_INK, "ytick.color": C_INK,
    })

    fig = plt.figure(figsize=(20, 12), facecolor="white")
    gs = fig.add_gridspec(2, 2, height_ratios=[1.3, 1], hspace=0.26,
                          wspace=0.13, left=0.05, right=0.93, top=0.92,
                          bottom=0.06)

    ax = fig.add_subplot(gs[0, :])
    ax.set_facecolor("#b0cfe8")
    im = ax.pcolormesh(np.ma.masked_invalid(canvas), cmap=cmap, norm=norm,
                       shading="auto", rasterized=True)

    def draw_road(a, zoom=False):
        for g, off in zip(gis_ids, offs):
            c0, c1 = domain_span(g)
            if zoom:
                # each domain's own dune line = row 0 of its grid = the line
                # road_setback is measured from
                a.plot([c0, c1], [off, off], color="k", lw=1.3, zorder=7,
                       solid_capstyle="butt")
            if not (FIRST_ROAD_DOMAIN <= g <= LAST_ROAD_DOMAIN):
                continue
            sb = setbacks.get(g)
            if sb is None:
                continue
            r0, r1 = road_rows(sb)
            a.add_patch(mpatches.Rectangle(
                (c0, off + r0), c1 - c0, max(r1 - r0, 1),
                facecolor=C_GROIN, edgecolor="none",
                alpha=0.95 if zoom else 0.85, zorder=6))
            if zoom and shelf_rows.get(g) is not None:
                a.plot([c0, c1], [off + shelf_rows[g]] * 2, color=C_MODEL,
                       lw=2.2, zorder=8, solid_capstyle="butt")

        for year, per in HISTORICAL_RELOCATIONS.items():
            st = RELOC_STYLE.get(year, dict(color=C_PIER, marker="v"))
            xs, ys = [], []
            for g, sb in sorted(per.items()):
                if g not in gis_ids:
                    continue
                off = offs[gis_ids.index(g)]
                r0, r1 = road_rows(sb)
                xs.append(domain_centre(g))
                ys.append(off + (r0 + r1) / 2.0)
            if xs:
                a.plot(xs, ys, marker=st["marker"], ms=10 if zoom else 8, lw=0,
                       color=st["color"], mec="k", mew=0.6, zorder=9,
                       label=f"{year} prescribed relocation")

    draw_road(ax)

    for name, lo, hi, colr in ISLAND_SECTIONS:
        if lo not in gis_ids or hi not in gis_ids:
            continue
        c0, _ = domain_span(lo)
        _, c1 = domain_span(hi)
        if colr:
            ax.axvspan(c0, c1, color=colr, alpha=0.16, zorder=1)
        ax.text((c0 + c1) / 2, 0.965, name, transform=ax.get_xaxis_transform(),
                ha="center", va="top", fontsize=9, color=C_INK,
                bbox=dict(fc="white", ec="none", alpha=0.82, pad=2))

    for lo, hi, _ in ZOOM_WINDOWS:
        if lo in gis_ids and hi in gis_ids:
            c0, _ = domain_span(lo)
            _, c1 = domain_span(hi)
            ax.add_patch(mpatches.Rectangle(
                (c0, 0), c1 - c0, canvas.shape[0], fill=False,
                edgecolor=C_INK, lw=1.4, ls="--", zorder=10))

    ticks = list(range(0, NUM_REAL_DOMAINS, 5))
    ax.set_xticks([domain_centre(gis_ids[i]) for i in ticks if i < len(gis_ids)])
    ax.set_xticklabels([str(gis_ids[i]) for i in ticks if i < len(gis_ids)],
                       fontsize=9)
    ax.set_xlim(0, canvas.shape[1])
    ax.set_ylim(*((0, canvas.shape[0]) if OCEAN_AT_BOTTOM
                  else (canvas.shape[0], 0)))
    step = max(int(canvas.shape[0] / 8), 1)
    yt = np.arange(0, canvas.shape[0] + 1, step)
    ax.set_yticks(yt)
    ax.set_yticklabels([f"{int((v + r_lo) * DX)}" for v in yt], fontsize=9)
    ax.set_ylabel(
        "Cross-shore position (m), shared canvas frame\n"
        + ("ocean below  ->  sound above" if OCEAN_AT_BOTTOM
           else "ocean above  ->  sound below"),
        fontsize=11)
    ax.set_xlabel("GIS domain  (S -> N,  Cape Point to Rodanthe)", fontsize=11)
    ax.set_title(f"NC-12 on the CASCADE initialization  |  setbacks {CHECK_YEAR}, "
                 f"topography {TOPO_DUNE_INIT_YEAR}, dune offsets applied",
                 fontsize=14, fontweight="bold", pad=10)

    h, l = ax.get_legend_handles_labels()
    d = dict(zip(l, h))
    d[f"NC-12 footprint ({CHECK_YEAR} setback)"] = mpatches.Patch(color=C_GROIN)
    d["road shelf detected in DEM"] = mpatches.Patch(color=C_MODEL)
    d["dune line (row 0, setback origin)"] = mpatches.Patch(color="k")
    ax.legend(d.values(), d.keys(), loc="lower left", fontsize=9,
              framealpha=0.93, facecolor="white", edgecolor="#cccccc")

    cax = fig.add_axes([0.94, 0.50, 0.010, 0.42])
    cb = fig.colorbar(im, cax=cax)
    cb.set_label("Elevation (m NAVD88)", fontsize=10, rotation=270, labelpad=14)
    cb.set_ticks([-1, 0, 1, 2, 3, 4])
    cb.outline.set_edgecolor("#cccccc")

    # ---- zoom panels (data-driven limits) ----
    for k, (lo, hi, label) in enumerate(ZOOM_WINDOWS):
        az = fig.add_subplot(gs[1, k])
        az.set_facecolor("#b0cfe8")
        if lo not in gis_ids or hi not in gis_ids:
            az.axis("off")
            continue
        c0, _ = domain_span(lo)
        _, c1 = domain_span(hi)
        az.pcolormesh(np.ma.masked_invalid(canvas), cmap=cmap, norm=norm,
                      shading="auto", rasterized=True)
        draw_road(az, zoom=True)

        c0i, c1i = int(c0), int(c1)
        win = canvas[:, c0i:c1i]
        with np.errstate(invalid="ignore"):
            colmax = np.nanmax(win, axis=1)
        land = np.where(colmax > 0.0)[0]
        if len(land):
            ylo = max(int(land[0]) - 3, 0)
            yhi = min(int(land[-1]) + 3, canvas.shape[0])
        else:
            ylo, yhi = 0, canvas.shape[0]

        az.set_xlim(c0, c1)
        az.set_ylim(*((ylo, yhi) if OCEAN_AT_BOTTOM else (yhi, ylo)))
        az.set_xticks([domain_centre(g) for g in range(lo, hi + 1)
                       if g in gis_ids])
        az.set_xticklabels([str(g) for g in range(lo, hi + 1) if g in gis_ids],
                           fontsize=9)
        step = max((yhi - ylo) // 7, 1)
        yt = np.arange(ylo, yhi, step)
        az.set_yticks(yt)
        az.set_yticklabels([f"{int((v + r_lo) * DX)}" for v in yt], fontsize=8)
        az.set_xlabel("GIS domain", fontsize=10)
        az.set_ylabel("Cross-shore position (m)", fontsize=10)
        az.set_title(f"GIS {lo}-{hi}  |  {label}", fontsize=11, pad=6)

    out1 = os.path.join(OUTPUT_DIR, f"HAT_road_on_island_{CHECK_YEAR}.png")
    fig.savefig(out1, dpi=190, facecolor="white", bbox_inches="tight")
    plt.close(fig)

    # =====================================================================
    # FIGURE 2: cross-shore profiles + road shelf (each domain's own frame)
    # =====================================================================
    doms = [g for g in PROFILE_DOMAINS if g in gis_ids]
    ncol, nrow = 4, int(np.ceil(len(doms) / 4))
    fig, axes = plt.subplots(nrow, ncol, figsize=(4.6 * ncol, 3.6 * nrow),
                             squeeze=False, facecolor="white")

    for k, g in enumerate(doms):
        a = axes[k // ncol][k % ncol]
        grid = grids[gis_ids.index(g)]
        prof = np.nanmean(grid, axis=1)
        std = np.nanstd(grid, axis=1)
        x = np.arange(len(prof)) * DX

        a.fill_between(x, prof - std, prof + std, color="#d9c7a3", alpha=0.5,
                       zorder=1, label="+/- 1 sd alongshore")
        a.fill_between(x, -2, prof, color="#e8dcc0", zorder=1)
        a.plot(x, prof, color=C_INK, lw=1.6, zorder=3)
        a.axhline(0, color=C_PIER, lw=1.0, ls=":", zorder=2)
        a.axhline(ROAD_ELEVATION, color=C_MODEL, lw=1.5, zorder=4,
                  label=f"ROAD_ELEVATION {ROAD_ELEVATION} m")

        sb = setbacks.get(g)
        if sb is not None:
            r0, r1 = road_rows(sb)
            a.axvspan(r0 * DX, r1 * DX, color=C_GROIN, alpha=0.40, zorder=5,
                      label=f"{CHECK_YEAR} setback {sb:.0f} m")

        row = shelf_rows.get(g)
        if row is not None:
            a.axvline(row * DX, color=C_MODEL, lw=2.2, ls="--", zorder=6,
                      label=f"DEM road shelf (row {row})")
            if sb is not None:
                d = (int(sb / DY) - row) * DY
                a.annotate(f"delta {d:+.0f} m", xy=(row * DX, ROAD_ELEVATION),
                           xytext=(6, 14), textcoords="raw_offset points",
                           fontsize=9, fontweight="bold",
                           bbox=dict(fc="white", ec=C_MODEL, alpha=0.9, pad=2))

        for year, per in HISTORICAL_RELOCATIONS.items():
            if g in per:
                st = RELOC_STYLE.get(year, dict(color=C_PIER))
                rr0, rr1 = road_rows(per[g])
                a.axvspan(rr0 * DX, rr1 * DX, color=st["color"], alpha=0.35,
                          zorder=5, label=f"{year} reloc {per[g]:.0f} m")

        a.set_title(f"GIS {g}", fontsize=11, fontweight="bold")
        a.set_xlabel("Distance landward of row 0 (m)", fontsize=9)
        a.set_ylabel("Elevation (m NAVD88)", fontsize=9)
        a.set_xlim(0, 400)
        a.set_ylim(-2, 4)
        a.tick_params(labelsize=8)
        a.grid(alpha=0.2)
        a.legend(fontsize=6.5, loc="upper right", framealpha=0.9)

    for k in range(len(doms), nrow * ncol):
        axes[k // ncol][k % ncol].axis("off")

    fig.suptitle(
        "Does the setback land on the road? Alongshore-mean profiles, each in "
        "its own domain frame.\nNC-12 is graded and shore-parallel, so it "
        "survives alongshore averaging as a flat shelf near 1.45 m. "
        "Red band should sit on the dashed line.",
        fontsize=12, y=0.99)
    fig.tight_layout(rect=[0, 0, 1, 0.93])
    out2 = os.path.join(OUTPUT_DIR, f"HAT_road_profile_check_{CHECK_YEAR}.png")
    fig.savefig(out2, dpi=190, facecolor="white")
    plt.close(fig)

    print("\n" + "=" * 78)
    print(f"Saved:\n  {out1}\n  {out2}")
    print("=" * 78)


if __name__ == "__main__":
    main()
