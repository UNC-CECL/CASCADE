r"""
HAT_setback_from_lines.py
===============================================================================
Road + dune geojsons in, CASCADE-ready setback CSV out. One script, one run per
vintage.

    setback_i = road_cell_i - dune_cell_i          both lines, SAME YEAR
    S         = median over the 50 profiles
    then make S fit the interior CASCADE actually has

This replaces HAT_rasterize_road_to_domains.py + HAT_road_setback_extract.py for
the job of producing a setback file. Those two still earn their keep for mask QC
and the island mosaic; they are not needed for this.

WHY IT COLLAPSES TO SOMETHING THIS SHORT
----------------------------------------
Everything the long chain did turned out to be either wrong or unnecessary:

  start_island   Measuring the 1984 road against the 2009 LiDAR dune returns
                 (S - R), where R is the 1984->2009 retreat. Where R > S the original
                 road lands seaward of the modern dune and the setback goes
                 negative -- which is what GIS 10-13 and 84-86 did, and those are
                 exactly the blocks NCDOT relocated in 1999 and 1989. The
                 negatives were the retreat. Both lines from the same aerials
                 cancels it.

  the shear      Rolls both masks by the same shear[i], so it cancels in the
                 difference. It buys nothing here -- and the shear plus the c0
                 water trim can DROP the 1984 dune line off the seaward end,
                 because that line sits seaward of the 2009 one by the retreat.
                 Measured raw, nothing can be dropped.

  the obliquity  The grids are north-up (rot = 0.00 on all 131) and the island
                 trends NNW, so both lines cross each domain diagonally. They
                 are diagonal TOGETHER, so the diagonal cancels in the
                 difference too. What does not cancel is the 1/cos(theta)
                 distance inflation -- 1% at 8 deg, 4% at 16, 11% at 26 -- but
                 the interior width carries the same factor, so the setback and
                 the grid stay consistent. Only rotated clips fix it. Note it,
                 do not correct it.

  elevation      A ~24 m mask on a 10 m resampled DEM averages the road crown
                 into the shoulder, and a 1984 line sampled in a 2009 DEM sits
                 on whatever that ground is now. It measured the DEM, not the
                 road. Dropped.

WHAT SURVIVES: two masks, one subtraction, and making the answer fit.

MASKS ARE GEOMETRY, NOT ELEVATION
---------------------------------
rasterize() uses the raster's transform and never reads a pixel. A 1984 line
lying over 2009 surf zone or a LiDAR gap still gets a mask. Position is what is
needed; elevation is not, and NoData is invisible here.

FITTING TO CASCADE
------------------
bulldoze() does:
    road_start = int(setback / 10)
    road_end   = road_start + int(road_width / 10)
    reads xyz_interior_grid[road_end + 1, :]

so S must satisfy  int(S/10) + int(w/10) + 2 <= interior_rows, and S >= 0.
A negative S is worse than an error: int(-110/10) = -11 and numpy WRAPS, so the
road gets bulldozed onto the sound side and the run finishes looking normal.

Both clamps are reported per domain. A clamp is a statement that the model
cannot hold the road where the aerials put it -- that is a result, not a
formatting step, and it belongs in the methods.

REQUIREMENTS
------------
  HAT_dune_topo_extractor.py importable (for OCEAN_LOC / paths / the run it
  belongs to), geopandas, rasterio, numpy, matplotlib
===============================================================================
"""

from __future__ import annotations

import csv
import os
import re
import sys
from pathlib import Path

import numpy as np

EXTRACTOR_DIR = None        # None = this script's folder
_ed = Path(EXTRACTOR_DIR) if EXTRACTOR_DIR else Path(__file__).resolve().parent
if str(_ed) not in sys.path:
    sys.path.insert(0, str(_ed))
try:
    import HAT_dune_topo_extractor as ext  # noqa: E402
except ModuleNotFoundError:
    raise SystemExit(
        f"\nCannot import HAT_dune_topo_extractor from:\n    {_ed}\n\n"
        f"This script reads OCEAN_LOC, CELL_SIZE_M and the run's topography "
        f"paths from it, so the setback lands in the same frame and belongs to "
        f"the same extraction. Put this file next to the extractor, or set "
        f"EXTRACTOR_DIR at the top.\n"
    )

import geopandas as gpd          # noqa: E402
import rasterio                  # noqa: E402
from rasterio.features import rasterize  # noqa: E402
import matplotlib                # noqa: E402
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.transforms import blended_transform_factory  # noqa: E402
from matplotlib.colors import FuncNorm  # noqa: E402


# =============================================================================
# CONFIG
# =============================================================================

GIS_DIR = ext.INIT_ROOT / "gis"

# The per-domain rasters gis-export-npy.py converted. Their affine is the
# grid; taking it from the raster means alignment is exact by construction
# rather than reconstructed.
CLIPRESAMPLE_ROOT = Path(
    r"C:\Users\hanna\OneDrive - University of North Carolina at Chapel Hill"
    r"\Ch1_CASCADE_hatteras\Hatteras_CASCADE_Input\domain_elevation"
    r"\2009_domain_clipresample"
)
TIF_GLOB = "resampled_*.tif"

# One entry per vintage. Road and dune MUST be the same year -- that is the
# whole point; a mismatch reintroduces the retreat.
VINTAGES = {
    1984: {"road": GIS_DIR / "nc12_1984.geojson",
           "dune": GIS_DIR / "duneline_1984.geojson"},
    2004: {"road": GIS_DIR / "nc12_2004.geojson",
           "dune": GIS_DIR / "duneline_2004.geojson"},
}
RUN_VINTAGES = [1984]        # which to process this run

OUT_ROOT = ext.RUN_DIR / "roads"     # -> dune_topo\2009_v2\roads\<year>\

FIRST_ROAD_DOMAIN, LAST_ROAD_DOMAIN = 9, 90

# Must match ROAD_WIDTH in HAT_hindcast_1984_2024_newplot.py.
ROAD_WIDTH_M = 20.0

# Both lines get the same buffer, so it cancels in the difference. It is here to
# keep a zero-width centreline from skipping cells diagonally on a 10 m grid.
LINE_BUFFER_M = 6.0
ALL_TOUCHED = True

AGG = "median"               # across the 50 profiles: "median" | "mean"
MIN_PROFILES = 10            # flag a domain measured on fewer
MAX_STD_M = 40.0             # flag profiles disagreeing by more than this

# Per-domain raw-lidar/geojson vs. CASCADE-frame QC panels.
# "flagged" = only domains with a FLOORED/CLAMPED/THIN/SCATTER flag (cheap,
# and exactly where the fit and the ground are most likely to disagree).
# "all" makes one per domain (slow -- 80ish panels). Or pass a list, e.g.
# QC_PANEL_DOMAINS = [10, 11, 12, 84, 85, 86]
QC_PANEL_DOMAINS = "flagged"

# True-geography island mosaic (one figure, every domain at its own real
# easting/northing -- no offset file needed). Companion to the per-domain
# panels above: those check one domain in detail, this checks the whole
# island reads sanely at once, the way HAT_road_island_mosaic.py's "geo" mode
# did for the original pipeline.
MAKE_ISLAND_MOSAIC = True
MOSAIC_ELEV_MIN_M, MOSAIC_ELEV_MAX_M = -1.0, 4.0
MOSAIC_SEA_LEVEL, MOSAIC_SEA_LEVEL_POS = 0.0, 0.35   # colormap sea-level pivot


C_MODEL, C_TOWN, C_PIER, C_GROIN, C_INK = ("#FF8C00", "#90AFC5", "#1565C0",
                                           "#B71C1C", "#1a1a2e")
C_DUNE, C_ROAD = "#8D6E33", C_PIER   # dune = sand-brown, road = reuse pier blue
SECTIONS = [((1, 6), "Cape Point"), ((7, 8), "Buxton"), ((9, 20), "Buxton-Avon"),
            ((21, 31), "Avon"), ((32, 67), "Avon-Tri-Village / Wimble Shoals"),
            ((68, 83), "Tri-Village"), ((84, 90), "Pea Island / N. Rodanthe")]


# =============================================================================
# HELPERS
# =============================================================================

def domain_id_from(text):
    m = re.search(r"(\d+)", str(text))
    return int(m.group(1)) if m else None


def find_domain_tifs(root):
    out = {}
    root = Path(root)
    if not root.exists():
        raise FileNotFoundError(f"CLIPRESAMPLE_ROOT not found:\n    {root}")
    for sub in sorted(root.iterdir()):
        if not sub.is_dir():
            continue
        d = domain_id_from(sub.name)
        tifs = sorted(sub.glob(TIF_GLOB))
        if d is not None and tifs:
            out[d] = tifs[0]
    return out


def load_line(path, label):
    if not Path(path).exists():
        raise FileNotFoundError(f"{label} not found:\n    {path}")
    gdf = gpd.read_file(path)
    if gdf.crs is None:
        raise ValueError(f"{label} has no CRS. Cannot align it to anything.")
    return gdf


def crs_to_metres(src):
    """EPSG:2264 is US survey FEET. A 6.0 buffer there is 1.8 m, not 6 m."""
    try:
        return float(src.crs.linear_units_factor[1])
    except Exception:
        return 1.0


def burn(gdf, src):
    g = gdf.to_crs(src.crs)
    to_m = crs_to_metres(src)
    geoms = [geom.buffer(LINE_BUFFER_M / to_m) for geom in g.geometry]
    return rasterize([(gm, 1) for gm in geoms],
                     out_shape=(src.height, src.width),
                     transform=src.transform, fill=0,
                     all_touched=ALL_TOUCHED, dtype="uint8")


def ocean_first(arr):
    """Index 0 = ocean, matching the frame the extractor works in."""
    a = ext.orient_ocean_right(arr.astype(float), ext.OCEAN_LOC)
    return a[:, ::-1]


def interior_rows_of(d):
    """Rows CASCADE will see, from the run's own saved topography."""
    p = Path(ext.TOPO_SAVE_PATH) / f"domain_{d}_topography_{ext.DEM_YEAR}.npy"
    if not p.exists():
        for cand in Path(ext.TOPO_SAVE_PATH).glob(f"domain_{d}_topography_*.npy"):
            p = cand
            break
        else:
            return -1
    return int(np.load(p, mmap_mode="r").shape[0])


# =============================================================================
# MEASURE
# =============================================================================

def measure(road_mask, dune_mask):
    """
    Per-profile setback: road cell minus dune cell, both raw, ocean-first.

    Median of each mask recovers its line's own position, so the shared buffer
    cancels. Both lines are diagonal together, so the obliquity cancels. No
    shear, no trim, no start_island -- nothing that could drop the seaward line.
    """
    n = min(road_mask.shape[0], dune_mask.shape[0])
    sb = np.full(n, np.nan)
    r_loc = np.full(n, -1, dtype=int)
    d_loc = np.full(n, -1, dtype=int)
    for i in range(n):
        rh = np.where(road_mask[i])[0]
        dh = np.where(dune_mask[i])[0]
        if rh.size == 0 or dh.size == 0:
            continue
        r_loc[i] = int(np.median(rh))
        d_loc[i] = int(np.median(dh))
        sb[i] = (r_loc[i] - d_loc[i]) * ext.CELL_SIZE_M
    return sb, r_loc, d_loc


def fit_to_cascade(S, interior_rows):
    """
    Make S something bulldoze() can index. Returns (S_fit, flag).

    A negative S is not a small error, it is a silent one: int(-110/10) = -11
    and numpy wraps the index to the BACK of the barrier, so the road gets
    bulldozed onto the sound side and the run completes looking normal.
    """
    if not np.isfinite(S):
        return np.nan, "NO_DATA"
    if interior_rows <= 0:
        return np.nan, "NO_INTERIOR"

    w = int(ROAD_WIDTH_M / ext.CELL_SIZE_M)
    max_rows = interior_rows - 2 - w            # leave bulldoze its check row
    S_max = max_rows * ext.CELL_SIZE_M

    if S < 0:
        return 0.0, f"FLOORED({S:.0f})"
    if S_max < 0:
        return np.nan, "TOO_NARROW"
    if S > S_max:
        return float(S_max), f"CLAMPED({S:.0f}->{S_max:.0f})"
    return float(S), ""


# =============================================================================
# MAIN
# =============================================================================

def run_vintage(year, tifs):
    paths = VINTAGES[year]
    print("\n" + "=" * 84)
    print(f"NC-12 {year} SETBACK  |  extraction run {ext.RUN_NAME}")
    print("=" * 84)
    road = load_line(paths["road"], f"road {year}")
    dune = load_line(paths["dune"], f"dune line {year}")
    print(f"  road : {paths['road'].name}  CRS {road.crs.to_string()}")
    print(f"  dune : {paths['dune'].name}  CRS {dune.crs.to_string()}")
    print(f"  Both are {year}. That is the point: measured against the 2009")
    print(f"  LiDAR dune instead, the setback would carry the 1984->2009")
    print(f"  retreat and go negative wherever the shoreline overtook the road.")

    rows = []
    for d in sorted(tifs):
        if not (FIRST_ROAD_DOMAIN <= d <= LAST_ROAD_DOMAIN):
            continue
        with rasterio.open(tifs[d]) as src:
            rm = ocean_first(burn(road, src)) > 0
            dm = ocean_first(burn(dune, src)) > 0

        sb, r_loc, d_loc = measure(rm, dm)
        good = sb[np.isfinite(sb)]
        S_raw = (float(np.median(good) if AGG == "median" else np.mean(good))
                 if good.size else np.nan)

        n_rows = interior_rows_of(d)
        S_fit, flag = fit_to_cascade(S_raw, n_rows)

        flags = [flag] if flag else []
        if good.size and good.size < MIN_PROFILES:
            flags.append(f"THIN({good.size})")
        if good.size and float(np.std(good)) > MAX_STD_M:
            flags.append(f"SCATTER({np.std(good):.0f})")

        rows.append(dict(
            domain=d, year=year,
            setback_raw_m=round(S_raw, 1) if np.isfinite(S_raw) else np.nan,
            setback_m=round(S_fit, 1) if np.isfinite(S_fit) else np.nan,
            interior_rows=n_rows,
            rows_needed=(int(S_fit / ext.CELL_SIZE_M)
                         + int(ROAD_WIDTH_M / ext.CELL_SIZE_M) + 2
                         if np.isfinite(S_fit) else -1),
            n_profiles=int(good.size), n_total=int(len(sb)),
            std_m=round(float(np.std(good)), 1) if good.size else np.nan,
            min_m=round(float(np.min(good)), 1) if good.size else np.nan,
            max_m=round(float(np.max(good)), 1) if good.size else np.nan,
            flags=",".join(flags), _sb=sb, _rloc=r_loc, _dloc=d_loc,
        ))

    if not rows:
        print("  nothing measured")
        return None

    # ---- report ----
    print(f"\n{'GIS':>4} {'raw':>8} {'fitted':>8} {'need':>5} {'rows':>5} "
          f"{'prof':>5} {'sd':>5}  flag")
    for r in rows:
        f_ = lambda v, w=8: ("--".rjust(w) if not np.isfinite(v)
                             else f"{v:>{w}.0f}")
        print(f"{r['domain']:>4} {f_(r['setback_raw_m'])} "
              f"{f_(r['setback_m'])} {r['rows_needed']:>5} "
              f"{r['interior_rows']:>5} {r['n_profiles']:>5} "
              f"{f_(r['std_m'], 5)}  {r['flags']}")

    sb = np.array([r["setback_m"] for r in rows if np.isfinite(r["setback_m"])])
    raw = np.array([r["setback_raw_m"] for r in rows
                    if np.isfinite(r["setback_raw_m"])])
    print()
    if raw.size:
        print(f"  raw    : min {raw.min():+.0f} | median {np.median(raw):+.0f} "
              f"| max {raw.max():+.0f} m")
    if sb.size:
        print(f"  fitted : min {sb.min():+.0f} | median {np.median(sb):+.0f} "
              f"| max {sb.max():+.0f} m")

    fl = [(r["domain"], r["flags"]) for r in rows if r["flags"]]
    if fl:
        print(f"\n  {len(fl)} domain(s) adjusted or flagged:")
        for d, f in fl:
            print(f"    GIS {d:>2}: {f}")
        print()
        print("    FLOORED(x)  raw setback was negative -- the 1984 road sits")
        print("                seaward of the 1984 dune line. Check the two")
        print("                geojsons in that domain; they may cross.")
        print("    CLAMPED(a->b)  the aerials put the road further back than")
        print("                the modelled barrier can hold. NOT a formatting")
        print("                step: it says the 2009 grid is too narrow for")
        print("                the road's real position, which is the width")
        print("                anachronism. Belongs in the methods.")
        print("    SCATTER(x)  profiles disagree -- the road or the dune line")
        print("                is not shore-parallel here, or a spur was")
        print("                digitized in.")
        print("    THIN(n)     measured on few profiles -- a line does not span")
        print("                this domain.")
    else:
        print("  No domain needed adjusting.")

    # ---- write ----
    out = Path(OUT_ROOT) / str(year)
    out.mkdir(parents=True, exist_ok=True)

    dom_rows = [{k: v for k, v in r.items() if not k.startswith("_")}
                for r in rows]
    p = out / f"RoadSetback_{year}_domains.csv"
    with open(p, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(dom_rows[0].keys()))
        w.writeheader()
        w.writerows(dom_rows)
    print(f"\n[out] {p}")

    prof = [dict(domain=r["domain"], profile=i,
                 setback_m=("" if not np.isfinite(v) else round(float(v), 1)))
            for r in rows for i, v in enumerate(r["_sb"])]
    p = out / f"RoadSetback_{year}_profiles.csv"
    with open(p, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(prof[0].keys()))
        w.writeheader()
        w.writerows(prof)
    print(f"[out] {p}")

    expected = list(range(FIRST_ROAD_DOMAIN, LAST_ROAD_DOMAIN + 1))
    by = {r["domain"]: r for r in rows}
    missing = [d for d in expected
               if d not in by or not np.isfinite(by[d]["setback_m"])]
    if missing:
        print(f"\n  NOT WRITING RoadSetback_{year}.csv: no setback for "
              f"domain(s) {missing}.")
        print(f"  The runner fills road_setbacks_full[] BY POSITION and drops")
        print(f"  the IDs with skiprows=1, so one gap shifts every domain north")
        print(f"  of it and nothing reports it.")
    else:
        ids = np.array(expected, dtype=float)
        vals = np.array([by[d]["setback_m"] for d in expected], dtype=float)
        p = out / f"RoadSetback_{year}.csv"
        np.savetxt(p, np.vstack([ids, vals]), delimiter=",", fmt="%.3f")
        print(f"[out] {p}")
        print(f"      2 x {len(ids)}  (row 0 = domain IDs, row 1 = setback m)")

    if QC_PANEL_DOMAINS:
        panel_dir = out / "qc_panels"
        panel_dir.mkdir(parents=True, exist_ok=True)
        by_d = {r["domain"]: r for r in rows}
        if QC_PANEL_DOMAINS == "flagged":
            targets = [d_ for d_, r in by_d.items() if r["flags"]]
        elif QC_PANEL_DOMAINS == "all":
            targets = list(by_d.keys())
        else:
            targets = list(QC_PANEL_DOMAINS)
        print(f"\n  {len(targets)} raw-vs-CASCADE panel(s) -> {panel_dir}")
        for d_ in sorted(targets):
            if d_ not in by_d or d_ not in tifs:
                continue
            p = qc_domain_diagnostic(d_, year, tifs, road, dune, by_d[d_],
                                     panel_dir)
            print(f"    [out] {p.name}")

    qc_figure(rows, year, out)

    if MAKE_ISLAND_MOSAIC:
        island_mosaic(year, tifs, road, dune, rows, out)

    return rows


def _elev_for_display(tif_path):
    """Read band 1 for a background image only -- never used for measurement."""
    with rasterio.open(tif_path) as src:
        elev = src.read(1).astype(float)
        if src.nodata is not None:
            elev[elev == src.nodata] = np.nan
        return elev, src


def qc_domain_diagnostic(d, year, tifs, road, dune, row, out_dir):
    """
    Two-panel check: does the raw LiDAR + geojson support the number that
    setback_m turned into.

    Left  (map view, native CRS): the actual road/dune lines drawn over the
           real elevation surface -- "do these lines trace the real road and
           the real dune, in the right place, this year".
    Right (CASCADE frame, ocean-first): the same elevation and the rasterized
           masks reoriented the way measure() actually saw them, with the
           per-profile road/dune positions, the fitted setback, the raw
           (pre-fit) setback if different, and the interior-row limit --
           "does the fitted placement make sense against the ground CASCADE
           will actually run on".

    Purely diagnostic -- nothing here changes setback_m. It exists to catch
    things a clean number can hide: a mask that's rasterized off the true
    line, a FLOORED(0) domain where 0 m lands the road between the two dune
    rows instead of behind them, or a CLAMPED domain where the "fitted"
    position is nowhere near where the aerials actually show the road.
    """
    tif = tifs[d]
    elev, src = _elev_for_display(tif)

    road_d = road.to_crs(src.crs)
    dune_d = dune.to_crs(src.crs)
    b = rasterio.transform.array_bounds(src.height, src.width, src.transform)
    extent = (b[0], b[2], b[1], b[3])  # left, right, bottom, top

    rm = burn(road, src) > 0
    dm = burn(dune, src) > 0
    elev_o = ocean_first(elev)
    rm_o = ocean_first(rm.astype(float))
    dm_o = ocean_first(dm.astype(float))

    r_loc, d_loc = row.get("_rloc"), row.get("_dloc")

    fig, (axL, axR) = plt.subplots(1, 2, figsize=(14, 6))

    # ---- left: map view, raw lines on raw ground ----
    axL.imshow(elev, cmap="gist_earth", extent=extent, origin="upper")
    road_d.plot(ax=axL, color=C_ROAD, lw=2.2, label=f"road {year} (raw)")
    dune_d.plot(ax=axL, color=C_DUNE, lw=2.2, label=f"dune {year} (raw)")
    axL.set_title(f"GIS {d}  |  {year}  |  map view (native CRS)", fontsize=10)
    axL.legend(loc="upper right", fontsize=8, framealpha=0.9)
    axL.set_xlabel("Easting (m)")
    axL.set_ylabel("Northing (m)")
    axL.set_aspect("equal")

    # ---- right: CASCADE frame, masks + fitted placement ----
    axR.imshow(elev_o, cmap="gist_earth", origin="upper", aspect="auto")
    axR.imshow(np.ma.masked_where(dm_o == 0, dm_o), cmap="Greens",
               alpha=0.45, origin="upper", aspect="auto")
    axR.imshow(np.ma.masked_where(rm_o == 0, rm_o), cmap="Blues",
               alpha=0.45, origin="upper", aspect="auto")

    if r_loc is not None and d_loc is not None:
        prof = np.arange(len(r_loc))
        good = (r_loc >= 0) & (d_loc >= 0)
        axR.plot(d_loc[good], prof[good], "o", ms=3, color=C_DUNE,
                 label="dune line (per profile)")
        axR.plot(r_loc[good], prof[good], "o", ms=3, color=C_ROAD,
                 label="road line (per profile)")

    w = int(ROAD_WIDTH_M / ext.CELL_SIZE_M)
    if np.isfinite(row["setback_m"]):
        s_col = row["setback_m"] / ext.CELL_SIZE_M
        axR.axvline(s_col, color=C_GROIN, lw=2, ls="-",
                    label=f"fitted road start ({row['setback_m']:.0f} m)")
        axR.axvline(s_col + w, color=C_GROIN, lw=1, ls="--",
                    label=f"fitted road end (+{ROAD_WIDTH_M:.0f} m)")
    if (np.isfinite(row["setback_raw_m"])
            and row["setback_raw_m"] != row["setback_m"]):
        axR.axvline(row["setback_raw_m"] / ext.CELL_SIZE_M, color=C_INK,
                    lw=1.4, ls=":",
                    label=f"raw setback ({row['setback_raw_m']:.0f} m)")
    if row["interior_rows"] > 0:
        axR.axvline(row["interior_rows"], color="0.3", lw=1.4, ls="-.",
                    label=f"interior limit ({row['interior_rows']} rows)")

    flag_txt = row["flags"] if row["flags"] else "no flags"
    axR.set_title(f"CASCADE frame  |  {flag_txt}", fontsize=10)
    axR.set_xlabel("cross-shore cell (0 = ocean)")
    axR.set_ylabel("alongshore profile")
    axR.legend(loc="upper right", fontsize=7, framealpha=0.9)

    fig.suptitle(f"Setback QC \u2014 GIS domain {d}, {year}", fontsize=12, y=1.02)
    fig.tight_layout()
    p = Path(out_dir) / f"qc_panel_GIS{d:03d}_{year}.png"
    fig.savefig(p, dpi=140, bbox_inches="tight")
    plt.close(fig)
    return p


def _mosaic_norm():
    """Terrain colormap with the sea-level pivot from HAT_road_on_island.py's
    poster palette, so the coastline reads clearly instead of as a smear."""
    cmap = plt.cm.terrain.copy()
    cmap.set_bad(color="#b0cfe8")

    def fwd(x):
        r = np.where(
            x < MOSAIC_SEA_LEVEL,
            MOSAIC_SEA_LEVEL_POS * (x - MOSAIC_ELEV_MIN_M)
            / (MOSAIC_SEA_LEVEL - MOSAIC_ELEV_MIN_M),
            MOSAIC_SEA_LEVEL_POS + (1 - MOSAIC_SEA_LEVEL_POS)
            * x / MOSAIC_ELEV_MAX_M)
        return np.where(np.isnan(x), np.nan, r)

    def inv(x):
        return np.where(
            x < MOSAIC_SEA_LEVEL_POS,
            MOSAIC_ELEV_MIN_M + (x / MOSAIC_SEA_LEVEL_POS)
            * (MOSAIC_SEA_LEVEL - MOSAIC_ELEV_MIN_M),
            (x - MOSAIC_SEA_LEVEL_POS) / (1 - MOSAIC_SEA_LEVEL_POS)
            * MOSAIC_ELEV_MAX_M)

    return cmap, FuncNorm((fwd, inv), vmin=MOSAIC_ELEV_MIN_M,
                          vmax=MOSAIC_ELEV_MAX_M)


def _stitch_geo(tifs, road, dune, domains):
    """
    True map: every domain at its own easting/northing, read straight from
    its tif's affine -- no offset file needed. Canvas axis0 tracks each
    raster's column dimension (easting), axis1 tracks its row dimension
    (northing), matching HAT_road_island_mosaic.py's stitch_geo().

    The canvas origin is rounded ONCE per domain (ri0/cj0), then stepped by
    whole cells with arange -- not rounded cell-by-cell. That one-time
    rounding is what HAT_road_island_mosaic.py had to add after per-cell
    rounding produced alternating striping on domains whose origin sat a
    half-cell off the shared grid (banker's rounding: round(0.5)=0,
    round(1.5)=2, round(2.5)=2, ...).

    Returns elev, dune_mask, road_mask (all shape (n_e, n_n)) and a dict of
    per-domain canvas index ranges for placing labels/markers consistently
    with what was actually drawn, rather than re-deriving positions.
    """
    cell_m = ext.CELL_SIZE_M
    xs, ys, shapes, transforms = [], [], {}, {}
    for d in domains:
        with rasterio.open(tifs[d]) as src:
            t = src.transform
            h, w = src.height, src.width
        shapes[d], transforms[d] = (h, w), t
        for i, j in ((0, 0), (0, w), (h, 0), (h, w)):
            xs.append(t.c + j * t.a + i * t.b)
            ys.append(t.f + j * t.d + i * t.e)
    xmin, xmax, ymin, ymax = min(xs), max(xs), min(ys), max(ys)

    n_e = int(np.ceil((xmax - xmin) / cell_m)) + 1
    n_n = int(np.ceil((ymax - ymin) / cell_m)) + 1
    elev = np.full((n_e, n_n), np.nan)
    dune_c = np.zeros((n_e, n_n), dtype=bool)
    road_c = np.zeros((n_e, n_n), dtype=bool)
    dom_extent = {}

    for d in domains:
        t = transforms[d]
        h, w = shapes[d]
        with rasterio.open(tifs[d]) as src:
            z = src.read(1).astype(float)
            if src.nodata is not None:
                z[z == src.nodata] = np.nan
            rm = burn(road, src) > 0
            dm = burn(dune, src) > 0

        ri0 = int(round((t.c + 0.5 * t.a - xmin) / cell_m))
        cj0 = int(round((t.f + 0.5 * t.e - ymin) / cell_m))
        ris = ri0 + np.arange(w)
        cjs = cj0 - np.arange(h)

        ok_r = (ris >= 0) & (ris < n_e)
        ok_c = (cjs >= 0) & (cjs < n_n)
        if not ok_r.any() or not ok_c.any():
            continue
        rr, cc = ris[ok_r], cjs[ok_c]
        elev[np.ix_(rr, cc)] = z[np.ix_(ok_c, ok_r)].T
        dune_c[np.ix_(rr, cc)] |= dm[np.ix_(ok_c, ok_r)].T
        road_c[np.ix_(rr, cc)] |= rm[np.ix_(ok_c, ok_r)].T
        dom_extent[d] = (float(rr.mean()), float(cc.mean()))

    return elev, dune_c, road_c, dom_extent


def _crop_to_data(elev, *masks):
    rows = np.where(~np.isnan(elev).all(axis=1))[0]
    cols = np.where(~np.isnan(elev).all(axis=0))[0]
    if not rows.size or not cols.size:
        return (elev,) + masks + (0, 0)
    r0, r1 = rows[0], rows[-1] + 1
    c0, c1 = cols[0], cols[-1] + 1
    cropped = tuple(m[r0:r1, c0:c1] for m in masks)
    return (elev[r0:r1, c0:c1],) + cropped + (r0, c0)


def island_mosaic(year, tifs, road, dune, rows, out_dir):
    """
    One true-geography planview of the whole run: every domain's real
    elevation at its real easting/northing, road and dune lines burned on
    top, and flagged domains (FLOORED/CLAMPED/THIN/SCATTER) marked directly
    on the ground they were flagged on.

    Deliberately the "geo" mode only (true map, no offset file). The
    model-frame planform -- domains stacked at Island_Dune_Offsets the way
    CASCADE actually sees them, like HAT_road_on_island.py's Figure 1 -- is a
    genuinely different question ("does the model's planform look right",
    not "is the road where the aerials put it") and needs the dune-offset
    CSV path plugged in; it isn't part of this script.
    """
    domains = sorted(d for d in tifs
                     if FIRST_ROAD_DOMAIN <= d <= LAST_ROAD_DOMAIN)
    if not domains:
        print("  [mosaic] no domains in range, skipping")
        return None

    elev, dune_c, road_c, dom_extent = _stitch_geo(tifs, road, dune, domains)
    ec, dc, rc, r0, c0 = _crop_to_data(elev, dune_c, road_c)

    cmap, norm = _mosaic_norm()
    fig, ax = plt.subplots(figsize=(20, 7))
    ax.set_facecolor("#b0cfe8")
    ax.imshow(np.ma.masked_invalid(ec), cmap=cmap, norm=norm, origin="upper",
             aspect="auto", interpolation="nearest")

    di, dj = np.where(dc)
    ri, rj = np.where(rc)
    ax.plot(dj, di, lw=0, marker="s", ms=0.7, color=C_DUNE, alpha=0.9,
           label=f"dune line {year}")
    ax.plot(rj, ri, lw=0, marker="s", ms=0.7, color=C_ROAD, alpha=0.9,
           label=f"road {year}")

    by_d = {r["domain"]: r for r in rows}
    for d, (rr, cc) in dom_extent.items():
        flag = by_d.get(d, {}).get("flags", "")
        if not flag:
            continue
        x, y = cc - c0, rr - r0
        ax.plot(x, y, marker="x", ms=9, mew=2, color=C_GROIN, zorder=9)
        ax.annotate(str(d), (x, y), fontsize=6, color=C_INK, ha="center",
                   xytext=(0, 7), textcoords="offset points")

    ax.set_title(
        f"NC-12 {year} on the real ground  |  domains {domains[0]}-"
        f"{domains[-1]}  |  \u00d7 = flagged (FLOORED/CLAMPED/THIN/SCATTER)",
        fontsize=12, fontweight="bold")
    ax.legend(loc="lower right", fontsize=8, framealpha=0.9)
    ax.set_xticks([])
    ax.set_yticks([])
    fig.tight_layout()
    p = Path(out_dir) / f"HAT_island_mosaic_{year}.png"
    fig.savefig(p, dpi=170, bbox_inches="tight")
    plt.close(fig)
    print(f"[out] {p}")
    return p


def qc_figure(rows, year, out_dir):
    rows = sorted(rows, key=lambda r: r["domain"])
    d = np.array([r["domain"] for r in rows])
    fig, (ax0, ax1) = plt.subplots(2, 1, figsize=(13, 8), sharex=True,
                                   gridspec_kw={"hspace": 0.12})

    for k, ((lo, hi), label) in enumerate(SECTIONS):
        for ax in (ax0, ax1):
            ax.axvspan(lo - 0.5, hi + 0.5, color=C_TOWN,
                       alpha=0.18 if k % 2 == 0 else 0.08, lw=0, zorder=0)
        tr = blended_transform_factory(ax0.transData, ax0.transAxes)
        ax0.text((lo + hi) / 2, 0.96, label, transform=tr, ha="center",
                 va="top", fontsize=7,
                 bbox=dict(boxstyle="round,pad=0.2", fc="white", ec="none",
                           alpha=0.85))

    for r in rows:
        v = r["_sb"][np.isfinite(r["_sb"])]
        if v.size:
            ax0.plot(np.full(v.size, r["domain"]), v, lw=0, marker=".", ms=2,
                     color="0.78", zorder=2)
    ax0.plot(d, [r["setback_raw_m"] for r in rows], color="0.45", lw=1.2,
             ls="--", zorder=3, label=f"raw  ({year} road vs {year} dune)")
    ax0.plot(d, [r["setback_m"] for r in rows], color=C_GROIN, lw=1.7,
             marker="o", ms=3, zorder=4, label="fitted to the interior")
    ax0.axhline(0, color=C_PIER, lw=1.0, ls=":", zorder=1)
    ax0.set_ylabel("NC-12 setback from\nthe dune line (m)")
    ax0.legend(loc="upper right", fontsize=8, framealpha=0.9)
    ax0.set_title(f"NC-12 {year} setbacks  —  both lines from the {year} "
                  f"aerials  —  run {ext.RUN_NAME}", fontsize=12)

    have = np.array([r["interior_rows"] for r in rows], dtype=float)
    need = np.array([r["rows_needed"] for r in rows], dtype=float)
    ax1.plot(d, have, color=C_PIER, lw=1.6, marker="o", ms=3,
             label="interior rows available")
    ax1.plot(d, need, color=C_MODEL, lw=1.6, marker="s", ms=3,
             label="rows needed (setback + road + check row)")
    clamped = np.array(["CLAMPED" in r["flags"] for r in rows])
    if clamped.any():
        ax1.plot(d[clamped], need[clamped], lw=0, marker="x", ms=11, mew=2.2,
                 color=C_GROIN, label="clamped to fit")
    ax1.set_ylabel("cross-shore rows")
    ax1.set_xlabel("GIS domain  (1 = Cape Point / south  ->  90 = Rodanthe / north)")
    ax1.legend(loc="upper right", fontsize=8, framealpha=0.9)
    ax1.set_xlim(d.min() - 0.5, d.max() + 0.5)

    fig.tight_layout()
    p = Path(out_dir) / f"HAT_setback_qc_{year}.png"
    fig.savefig(p, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"[out] {p}")


def main():
    print("=" * 84)
    print("NC-12 SETBACKS FROM ROAD + DUNE GEOJSONS")
    print("=" * 84)
    print(f"  extraction run : {ext.RUN_NAME}")
    print(f"  topography     : {ext.TOPO_SAVE_PATH}")
    print(f"  OCEAN_LOC      : {ext.OCEAN_LOC}   CELL_SIZE_M: {ext.CELL_SIZE_M}")
    print(f"  out            : {OUT_ROOT}")

    tifs = find_domain_tifs(CLIPRESAMPLE_ROOT)
    print(f"  {len(tifs)} domain raster(s)")
    if not tifs:
        return
    with rasterio.open(next(iter(tifs.values()))) as s:
        to_m = crs_to_metres(s)
        print(f"  raster CRS     : {s.crs.to_string()[:40]}...  "
              f"unit {to_m:.4f} m")
        print(f"  cell           : {abs(s.transform.a) * to_m:.2f} m  |  "
              f"buffer {LINE_BUFFER_M} m -> {LINE_BUFFER_M / to_m:.2f} CRS units")

    for year in RUN_VINTAGES:
        if year not in VINTAGES:
            print(f"\n[skip] {year}: not in VINTAGES")
            continue
        try:
            run_vintage(year, tifs)
        except FileNotFoundError as e:
            print(f"\n[skip] {year}: {e}")

    print("\n" + "=" * 84)
    print("NEXT: point HAT_hindcast_1984_2024_newplot.py at this run")
    print("=" * 84)
    print(f"  ROAD_SETBACK_FILE = "
          f"r\"{OUT_ROOT / str(RUN_VINTAGES[0]) / f'RoadSetback_{RUN_VINTAGES[0]}.csv'}\"")
    print(f"  topography        = r\"{ext.TOPO_SAVE_PATH}\"")
    print(f"  dunes             = r\"{ext.DUNE_SAVE_PATH}\"")
    print("  All three from the same run.")


if __name__ == "__main__":
    main()
