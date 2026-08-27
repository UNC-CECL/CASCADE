r"""
HAT_rasterize_road_to_domains.py
===============================================================================
Burn a road geojson onto the exact per-domain grids produced by the CASCADE
Clip & Resample Domains tool, so the setback can be measured along the same
profiles that define the interior.

THE ONLY SCRIPT THAT MASKS THE ROAD. Set ROAD_YEAR in CONFIG and run it:

    python HAT_rasterize_road_to_domains.py

The masks it writes are consumed by HAT_road_offset_from_dune_start.py (the
setback) and by the audit/ and figures/ scripts.

The road-side twin of gis-export-npy.py: same folder walk, same grids, same
array orientation -- rasterizing NC-12 instead of converting the DEM.

THE ONLY THING THAT MATTERS
---------------------------
Grid alignment. The road mask must be cell-for-cell identical to the elevation
array: same shape, same affine transform, same rotation (if any).

This script never guesses. It opens each domain's resampled_*.tif and takes the
affine straight from the raster, so alignment holds whether the Clip & Resample
tool produced north-up or shoreline-rotated grids.

OUTPUT LAYOUT
-------------
Everything for one road vintage lands in ONE folder, so comparing vintages
means comparing two directories:

  data\hatteras_init\4-mgmt-forcing\roads\raster\1984\
      RUN_MANIFEST.txt                        every setting that made this folder
      masks\
          domain_9_road_1984.npy              <- the setback scripts read these
          domain_10_road_1984.npy
          ...
      HAT_road_mask_diagnostics_1984.csv      per-domain numbers, incl. elevation
      HAT_road_mask_summary_1984.png          all domains on one page
      figures\
          domain_009_road_mask.png            per-domain map + profile

NOTE ON THE ELEVATION COLUMNS
-----------------------------
road_elev_* is the elevation of the cells the mask landed on, MHW-relative
(the extractor subtracts MHW_M = 0.36 before anything else, so this matches
the frame CASCADE runs in).

Do NOT read a low value as "misregistered" without looking at the figures. Two
reasons the crown may not survive:

  1. ROAD_BUFFER_M = 6 plus all_touched=True gives a ~24 m wide mask on 10 m
     cells. NC-12 is ~8 m. Most of every "road" cell is shoulder and adjacent
     ground, and the median follows them.
  2. The DEM was resampled to 10 m. An 8 m road inside a 10 m cell is averaged
     with whatever else is in that cell.

So road_elev may be measuring the DEM's resolution rather than the road. The
LOW_ELEV flag is INFORMATIONAL. What actually proves registration is the
figures: does the mask trace the island, and does it sit where you digitized it.

The setback does not depend on any of this -- it comes from the geojson
geometry. The elevation only matters if you want a per-domain road_ele.

REQUIREMENTS
------------
  rasterio, geopandas, shapely, numpy, matplotlib
===============================================================================
"""

import csv
import os
import re
from datetime import datetime
from pathlib import Path

import numpy as np
import geopandas as gpd
import rasterio
from rasterio.features import rasterize

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.transforms import blended_transform_factory

# =============================================================================
# CONFIG
# =============================================================================

# Derived from this file's own location, never hardcoded: this constant had been
# reduced to Path("/") -- every path below resolved to the filesystem root and
# the script died on "missing clipresample root". scripts/input_prep/
# 4-mgmt-forcings/road_offset/1-produce/<this file> -> parents[5].
PROJECT_ROOT = Path(__file__).resolve().parents[5]
INIT_ROOT = PROJECT_ROOT / "data" / "hatteras_init"
DOMAIN_ROOT = INIT_ROOT / "1-barrier3d-domains"
ROADS_ROOT = INIT_ROOT / "4-mgmt-forcing" / "road_offset"

# Same root gis-export-npy.py walks: one subfolder per domain, each holding
# a resampled_*.tif. The REPO copy, not the OneDrive original, so a run does
# not need the network drive.
# MOVED 2026-08-25: 1-barrier3d-domains went period-first and the pre-90-domain
# legacy went under superseded/. Path repointed so this script keeps reading
# EXACTLY what it read before - no road number moves because of the reorg.
# MOVED 2026-08-26 out of superseded/. These are LIVE INPUTS, not an
# archive - four scripts read them, including the one that builds the
# road masks the dune/topo extractor requires. Keeping them under a
# directory called 'superseded' invited exactly the deletion this move
# prevents.
CLIPRESAMPLE_ROOT = DOMAIN_ROOT / "domain-clips-1m"
TIF_GLOB = "resampled_*.tif"

# A GEOMETRY REFERENCE, NOT "THE ARRAYS THE EXTRACTOR READS" (corrected
# 2026-08-26). That is what this comment used to claim, and it stopped being
# true when 1-barrier3d-domains went period-first: the extractor now reads
# <product>/npy-arrays, one per period, and nothing reads this directory as an
# input to a forcing.
#
# What the cross-check needs is the GRID, not the elevations - the mask must
# land on the same rows and columns the extractor will later shear and trim.
# Every one of these arrays is (50, 200), identical in npy-arrays_2009_unfilled
# and in both products, because all three are clipped from the same
# resampled_*.tif boxes and the fills change values, never extents. So this
# stays a valid alignment reference and is deliberately NOT repointed at a
# product: doing that would make the road masks - which BOTH periods share -
# depend on one period's DEM.
#
# It held 131 arrays against the live 90; the extras were pre-90-domain legacy
# and were purged on 2026-08-26, leaving domains 1-90. find_elev_npy() looks
# up by domain id, so the count never mattered here either way.
ELEV_NPY_DIR = DOMAIN_ROOT / "npy-arrays_2009_unfilled"

# --- ROAD ---------------------------------------------------------------
# THE ONE LINE TO CHANGE PER VINTAGE. Edit it and re-run; do NOT save a
# per-year copy of this file. Everything below is derived from it, and the
# output lands in its own road_offset\raster\<year>\ folder with a RUN_MANIFEST, so
# the vintages stay separable without a second script.
#
# This file used to be driven by HAT_rasterize_road_run.py, which patched these
# globals at import. That driver is gone; its settings are folded in here, so
# there is ONE rasterization implementation and one place for the year and the
# paths. The rules below are the ones road_offset\raster\1984\RUN_MANIFEST.txt
# records, so a 2004 run is made by the same rules as the 1984 masks on disk.
#
# Overridable from the shell so building BOTH vintages needs no edit and no
# second copy of this file -- which is the rule above, not an exception to it:
#     HAT_ROAD_YEAR=1984 python HAT_rasterize_road_to_domains.py
# With the variable unset the constant below is what runs, so the file still
# reads as "the one line to change".
ROAD_YEAR = int(os.environ.get("HAT_ROAD_YEAR", 2004))
ROAD_GEOJSON = ROADS_ROOT / "raw_offset" / str(ROAD_YEAR) / f"nc12_{ROAD_YEAR}.geojson"

# --- OUTPUT LAYOUT ------------------------------------------------------
OUT_ROOT = ROADS_ROOT / "raster" / str(ROAD_YEAR)
MASK_DIR = OUT_ROOT / "masks"          # <- the setback scripts read from here
FIG_DIR = OUT_ROOT / "figures"
DIAG_CSV = OUT_ROOT / f"HAT_road_mask_diagnostics_{ROAD_YEAR}.csv"
SUMMARY_FIG = OUT_ROOT / f"HAT_road_mask_summary_{ROAD_YEAR}.png"
MANIFEST_PATH = OUT_ROOT / "RUN_MANIFEST.txt"

SAVE_SUMMARY_FIG = True
SAVE_DOMAIN_FIGS = True
# Per-domain map figures. None = every road domain (slow, ~110 figures);
# a list = just those. Defaults to the domains this investigation cares about:
# the 1999 relocation block, the crash suspect, the LOW_ELEV heartland, the
# NODATA cases, and the 1989 Pea Island block.
QC_DOMAINS = [9, 11, 13, 15, 16, 22, 35, 74, 78, 79, 84, 85, 86, 87]

# --- RASTERIZATION ------------------------------------------------------
# Widen the road before burning. The geojson is a zero-width centerline, and a
# zero-width line through a 10 m grid can skip cells diagonally, leaving gaps a
# profile falls straight through. NC-12 is ~2 lanes plus shoulders.
# Set 0.0 to burn the bare centerline: with ALL_TOUCHED the path stays
# continuous, and road_elev then samples only the cells the road crosses. That
# is the cheap test for whether the crown is resolved at 10 m.
ROAD_BUFFER_M = 6.0
ALL_TOUCHED = True

# --- DATUM / DIAGNOSTICS ------------------------------------------------
MHW_M = 0.36              # cascade_export_npy -> extractor subtracts this
BERM_ELEV_NAVD_M = 1.70   # m NAVD88; = 1.34 m MHW
NODATA_VALUE = -10.0      # gis-export-npy.py: nodata_to_value=-10
ROAD_ELEV_MIN_M = 0.6     # MHW; INFORMATIONAL flag only -- see the header
FIRST_ROAD_DOMAIN, LAST_ROAD_DOMAIN = 9, 90

# Display only. Must match OCEAN_LOC in HAT_dune_topo_extractor.py so the
# figures show the same orientation the extractor works in.
OCEAN_LOC = "right"

# --- STYLE (HAT_hindcast_1984_2024.py palette) --------------------------
C_MODEL, C_TOWN, C_WIMBLE, C_PIER, C_GROIN = ("#FF8C00", "#90AFC5", "#E0A800",
                                              "#1565C0", "#B71C1C")
C_INK = "#1a1a2e"

SECTIONS = [
    ((1, 6),   "Cape Point"),
    ((7, 8),   "Buxton"),
    ((9, 20),  "Buxton-Avon"),
    ((21, 31), "Avon"),
    ((32, 67), "Avon-Tri-Village / Wimble Shoals"),
    ((68, 83), "Tri-Village"),
    ((84, 90), "Pea Island / N. Rodanthe"),
]


# =============================================================================
# HELPERS
# =============================================================================

def domain_id_from(text):
    m = re.search(r"(\d+)", str(text))
    return int(m.group(1)) if m else None


def section_for(d):
    for (lo, hi), label in SECTIONS:
        if lo <= d <= hi:
            return label
    return ""


def find_domain_tifs(root):
    """Walk the clipresample folders exactly as gis-export-npy.py does."""
    out = {}
    root = Path(root)
    if not root.exists():
        raise FileNotFoundError(f"CLIPRESAMPLE_ROOT not found: {root}")
    for sub in sorted(root.iterdir()):
        if not sub.is_dir():
            continue
        d = domain_id_from(sub.name)
        if d is None:
            continue
        tifs = sorted(sub.glob(TIF_GLOB))
        if not tifs:
            print(f"  [skip] no {TIF_GLOB} in {sub.name}")
            continue
        if len(tifs) > 1:
            print(f"  [warn] {sub.name}: {len(tifs)} matches; using {tifs[0].name}")
        out[d] = tifs[0]
    return out


def find_elev_npy(npy_dir, domain_id, folder_name):
    if not npy_dir:
        return None
    p = Path(npy_dir)
    for cand in (f"domain_{domain_id}.npy",
                 f"{folder_name}_resampled.npy",
                 f"domain_{domain_id}_resampled.npy",
                 f"{domain_id}_resampled.npy"):
        if (p / cand).exists():
            return p / cand
    return None


def load_road():
    if not Path(ROAD_GEOJSON).exists():
        raise FileNotFoundError(f"ROAD_GEOJSON not found: {ROAD_GEOJSON}")
    gdf = gpd.read_file(ROAD_GEOJSON)
    if gdf.crs is None:
        raise ValueError(f"{ROAD_GEOJSON} has no CRS. Cannot align it to "
                         f"anything. Set one in GIS and re-export.")
    print(f"  road: {len(gdf)} feature(s), CRS {gdf.crs.to_string()}")
    return gdf


def crs_to_metres(src):
    """
    Metres per linear unit of src's CRS. EPSG:2264 (NC State Plane) is US
    SURVEY FEET, so buffering by 6.0 there is 6 ft = 1.8 m, not 6 m.
    """
    try:
        return float(src.crs.linear_units_factor[1])
    except Exception:
        print("  !! could not read the raster's linear unit; assuming metres")
        return 1.0


def burn(road_gdf, src):
    """Rasterize into src's exact grid. The affine carries any rotation."""
    g = road_gdf.to_crs(src.crs)
    geoms = list(g.geometry)
    if ROAD_BUFFER_M > 0:
        geoms = [geom.buffer(ROAD_BUFFER_M / crs_to_metres(src))
                 for geom in geoms]
    return rasterize(
        [(geom, 1) for geom in geoms],
        out_shape=(src.height, src.width),
        transform=src.transform,
        fill=0,
        all_touched=ALL_TOUCHED,
        dtype="uint8",
    )


def ocean_first(arr):
    """
    (n_along, n_cross) with index 0 = ocean, matching what the extractor works
    in after orient_ocean_right() -> [:, ::-1]. Display only.
    """
    if OCEAN_LOC == "right":
        return arr[:, ::-1]
    if OCEAN_LOC == "left":
        return arr
    if OCEAN_LOC == "bottom":
        return np.rot90(arr)[:, ::-1]
    if OCEAN_LOC == "top":
        return np.rot90(arr, -1)[:, ::-1]
    raise ValueError(f"OCEAN_LOC must be right/left/top/bottom, got {OCEAN_LOC!r}")


# =============================================================================
# FIGURES
# =============================================================================

def domain_figure(d, elev_raw, mask_raw, rec, fig_dir: Path):
    """
    Confirm, by eye, that the mask landed on the island where it should.

    Left  : elevation map, ocean at the bottom, road mask overlaid.
    Right : alongshore-mean cross-shore profile, with the road's cross-shore
            span shaded and the berm marked.
    """
    z = ocean_first(elev_raw) - MHW_M
    m = ocean_first(mask_raw) > 0
    z_disp = np.where(z <= NODATA_VALUE - MHW_M + 1e-6, np.nan, z)
    n_along, n_cross = z.shape

    fig, (ax_map, ax_prof) = plt.subplots(
        1, 2, figsize=(13, 8), sharey=True,
        gridspec_kw={"width_ratios": [1.3, 1.0], "wspace": 0.06})

    finite = z_disp[np.isfinite(z_disp)]
    vmax = float(np.percentile(finite, 99)) if finite.size else 3.0

    im = ax_map.imshow(np.ma.masked_invalid(z_disp.T), aspect="auto",
                       origin="lower", extent=[-0.5, n_along - 0.5,
                                               -0.5, n_cross - 0.5],
                       cmap="terrain", vmin=-1.0, vmax=max(vmax, 2.0))
    ai, ci = np.where(m)
    ax_map.plot(ai, ci, lw=0, marker="s", ms=2.5, color=C_GROIN,
                label=f"NC-12 {ROAD_YEAR} mask")
    ax_map.set_xlabel("alongshore cell")
    ax_map.set_ylabel("cross-shore cell  (0 = ocean, landward up)")
    ax_map.legend(loc="upper right", fontsize=9, framealpha=0.9)

    y = np.arange(n_cross)
    with np.errstate(invalid="ignore"):
        import warnings
        with warnings.catch_warnings():
            # columns landward of the clip are all NoData -> all-NaN slice
            warnings.simplefilter("ignore", RuntimeWarning)
            prof = np.nanmean(z_disp, axis=0)
    ax_prof.plot(z_disp.T, y, color="0.8", lw=0.5)
    ax_prof.plot(prof, y, color="k", lw=2.0, label="alongshore-mean")
    ax_prof.axvline(BERM_ELEV_NAVD_M - MHW_M, color=C_GROIN, ls=":", lw=1.4,
                    label=f"berm {BERM_ELEV_NAVD_M} m NAVD88")
    ax_prof.axvline(0, color=C_PIER, ls="--", lw=1.0, label="MHW")
    if ci.size:
        ax_prof.axhspan(ci.min(), ci.max(), color=C_GROIN, alpha=0.18,
                        zorder=0, label="road cross-shore span")
    ax_prof.set_xlabel("elev (m MHW)")
    ax_prof.set_xlim(-1.5, max(vmax, 2.0) + 0.5)
    ax_prof.set_ylim(-0.5, n_cross - 0.5)
    ax_prof.legend(loc="upper right", fontsize=8, framealpha=0.9)
    plt.setp(ax_prof.get_yticklabels(), visible=False)

    try:
        fig.colorbar(im, ax=[ax_map, ax_prof], location="right",
                     fraction=0.035, pad=0.02, label="elev (m MHW)")
    except (TypeError, ValueError):
        fig.colorbar(im, ax=ax_prof, label="elev (m MHW)")

    sec = section_for(d)
    fig.suptitle(
        f"domain_{d}  ({sec})  —  NC-12 {ROAD_YEAR} mask\n"
        f"{rec['n_cells']} cells on {rec['n_profiles']} profiles  |  "
        f"cross-shore {rec['road_cs_min']}-{rec['road_cs_max']} "
        f"(median {rec['road_cs_median']})  |  "
        f"road elev median {rec['road_elev_mhw']:.2f} m MHW"
        f"{'  |  ' + rec['flags'] if rec['flags'] else ''}",
        fontsize=11)

    fig_dir.mkdir(parents=True, exist_ok=True)
    fig.savefig(fig_dir / f"domain_{d:03d}_road_mask.png", dpi=130,
                bbox_inches="tight")
    plt.close(fig)


def summary_figure(recs, path: Path):
    """Every domain on one page: where the road sits, what it sits on."""
    recs = sorted([r for r in recs if r["n_cells"] > 0],
                  key=lambda r: r["domain"])
    if not recs:
        return
    d = np.array([r["domain"] for r in recs])

    fig, (ax0, ax1, ax2) = plt.subplots(3, 1, figsize=(13, 9.5), sharex=True,
                                        gridspec_kw={"hspace": 0.12})

    for k, ((lo, hi), label) in enumerate(SECTIONS):
        if not ((d >= lo) & (d <= hi)).any():
            continue
        for ax in (ax0, ax1, ax2):
            ax.axvspan(lo - 0.5, hi + 0.5, color=C_TOWN,
                       alpha=0.18 if k % 2 == 0 else 0.08, lw=0, zorder=0)
        tr = blended_transform_factory(ax0.transData, ax0.transAxes)
        ax0.text((lo + hi) / 2, 0.96, label, transform=tr, ha="center",
                 va="top", fontsize=8,
                 bbox=dict(boxstyle="round,pad=0.2", fc="white", ec="none",
                           alpha=0.85))

    # 1) where the road sits, cross-shore
    ax0.fill_between(d, [r["road_cs_min"] for r in recs],
                     [r["road_cs_max"] for r in recs],
                     color=C_GROIN, alpha=0.25, step="mid",
                     label="road cross-shore span")
    ax0.step(d, [r["road_cs_median"] for r in recs], where="mid",
             color=C_GROIN, lw=1.6, label="road median cell")
    ax0.step(d, [r["land_cs_max"] for r in recs], where="mid", color=C_PIER,
             lw=1.4, ls="--", label="landwardmost cell above MHW")
    ax0.set_ylabel("cross-shore cell\n(0 = ocean)")
    ax0.legend(loc="upper right", fontsize=8, framealpha=0.9)
    ax0.set_title(f"NC-12 {ROAD_YEAR} rasterized onto the domain grids  —  "
                  f"{len(recs)} domains  —  buffer {ROAD_BUFFER_M:.0f} m",
                  fontsize=12)

    # 2) what it sits on
    ax1.fill_between(d, [r["road_elev_min_mhw"] for r in recs],
                     [r["road_elev_max_mhw"] for r in recs],
                     color=C_MODEL, alpha=0.22, label="min-max")
    ax1.plot(d, [r["road_elev_mhw"] for r in recs], color=C_MODEL, lw=1.6,
             marker="o", ms=3, label="median")
    ax1.axhline(BERM_ELEV_NAVD_M - MHW_M, color=C_GROIN, ls=":", lw=1.4,
                label=f"berm ({BERM_ELEV_NAVD_M} m NAVD88 = "
                      f"{BERM_ELEV_NAVD_M - MHW_M:.2f} m MHW)")
    ax1.axhline(ROAD_ELEV_MIN_M, color="0.5", ls="--", lw=1.0,
                label=f"LOW_ELEV flag ({ROAD_ELEV_MIN_M} m MHW)")
    ax1.axhline(0, color=C_PIER, ls="-", lw=0.8, label="MHW")
    ax1.set_ylabel("elevation under\nthe road mask (m MHW)")
    ax1.legend(loc="upper right", fontsize=7, framealpha=0.9, ncol=2)

    # 3) coverage
    ax2.plot(d, [r["n_profiles"] for r in recs], color=C_PIER, lw=1.6,
             marker="o", ms=3, label="profiles with road (of 50)")
    nd = np.array([r["n_nodata"] for r in recs])
    if (nd > 0).any():
        ax2.plot(d[nd > 0], nd[nd > 0], lw=0, marker="v", ms=7, color=C_GROIN,
                 label="road cells on NoData")
    ax2.set_ylabel("cells")
    ax2.set_xlabel("domain  (1 = Cape Point / south  ->  90 = Rodanthe / north)")
    ax2.legend(loc="upper right", fontsize=8, framealpha=0.9)
    ax2.set_xlim(d.min() - 0.5, d.max() + 0.5)

    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=140, bbox_inches="tight")
    plt.close(fig)
    print(f"[summary] {path}")


# =============================================================================
# OUTPUT
# =============================================================================

def write_diagnostics(recs, path: Path):
    if not recs:
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    cols = [k for k in recs[0] if not k.startswith("_")]
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=cols, extrasaction="ignore")
        w.writeheader()
        w.writerows(recs)
    print(f"[diag] {path}")


def write_manifest(recs, path: Path):
    cfg = {
        "run_time": datetime.now().isoformat(timespec="seconds"),
        "script": Path(__file__).name,
        "ROAD_YEAR": ROAD_YEAR,
        "ROAD_GEOJSON": str(ROAD_GEOJSON),
        "CLIPRESAMPLE_ROOT": str(CLIPRESAMPLE_ROOT),
        "ELEV_NPY_DIR": str(ELEV_NPY_DIR),
        "ROAD_BUFFER_M": ROAD_BUFFER_M,
        "ALL_TOUCHED": ALL_TOUCHED,
        "MHW_M": MHW_M,
        "BERM_ELEV_NAVD_M": BERM_ELEV_NAVD_M,
        "NODATA_VALUE": NODATA_VALUE,
        "ROAD_ELEV_MIN_M": ROAD_ELEV_MIN_M,
        "OCEAN_LOC (display)": OCEAN_LOC,
        "road domains": f"{FIRST_ROAD_DOMAIN}-{LAST_ROAD_DOMAIN}",
    }
    width = max(len(k) for k in cfg)
    lines = [
        "=" * 78,
        f"NC-12 {ROAD_YEAR} rasterization",
        "=" * 78, "",
        "CONTENTS",
        r"  masks\                    road masks -> HAT_road_offset_from_dune_start.py",
        f"  {DIAG_CSV.name}   per-domain numbers incl. elevation",
        f"  {SUMMARY_FIG.name}   all domains on one page",
        r"  figures\                  per-domain map + profile",
        "",
        "READ THIS BEFORE TRUSTING road_elev_*",
        "  ROAD_BUFFER_M plus all_touched gives a mask far wider than NC-12,",
        "  and the DEM is resampled to 10 m. The road crown may not survive",
        "  either. A low road_elev may be the DEM's resolution, not a",
        "  registration error. The figures are what prove registration.",
        "  The SETBACK does not depend on elevation at all -- it comes from",
        "  the geojson geometry.",
        "",
        "SETTINGS",
    ]
    lines += [f"  {k:<{width}} : {v}" for k, v in cfg.items()]

    road = [r for r in recs
            if FIRST_ROAD_DOMAIN <= r["domain"] <= LAST_ROAD_DOMAIN]
    lines += ["", f"DOMAINS: {len(recs)} rasterized, {len(road)} in the road window"]
    for name in ("NO_ROAD", "NODATA", "PATCHY", "LOW_ELEV", "SHAPE_MISMATCH",
                 "NO_ELEV_CHECK"):
        hit = [r["domain"] for r in road if name in (r["flags"] or "")]
        if hit:
            lines.append(f"  {name:<14}: {hit}")
    lines.append("")
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines))
    print(f"[manifest] {path}")


# =============================================================================
# MAIN
# =============================================================================

def main():
    # Check the inputs BEFORE mkdir. Carried over from the retired
    # HAT_rasterize_road_run.py driver, and it is not decoration: the mkdir
    # below is parents=True, so a wrong PROJECT_ROOT silently builds a whole
    # stray tree and the run looks normal. That is exactly how the leftover
    # C:\data\hatteras_init\ came to exist. Fail loudly instead.
    for label, path in [("clipresample root", CLIPRESAMPLE_ROOT),
                        ("elevation npy dir", ELEV_NPY_DIR),
                        ("road geojson", ROAD_GEOJSON)]:
        if not path.exists():
            raise SystemExit(f"\nmissing {label}:\n    {path}\n")

    for p in (OUT_ROOT, MASK_DIR, FIG_DIR):
        p.mkdir(parents=True, exist_ok=True)

    print("=" * 82)
    print(f"RASTERIZE NC-12 {ROAD_YEAR} ONTO DOMAIN GRIDS")
    print("=" * 82)
    road = load_road()

    tifs = find_domain_tifs(CLIPRESAMPLE_ROOT)
    print(f"  found {len(tifs)} domain raster(s)")
    if not tifs:
        return

    with rasterio.open(next(iter(tifs.values()))) as _s:
        _to_m = crs_to_metres(_s)
        try:
            _unit = _s.crs.linear_units
        except Exception:
            _unit = "?"
        print(f"  raster CRS unit: {_unit} ({_to_m:.6f} m per unit)")
        print(f"  cell size : {abs(_s.transform.a):.3f} {_unit} "
              f"= {abs(_s.transform.a) * _to_m:.2f} m")
        print(f"  road buffer: {ROAD_BUFFER_M:.1f} m "
              f"-> {ROAD_BUFFER_M / _to_m:.2f} {_unit} in this CRS")
        if abs(abs(_s.transform.a) * _to_m - 10.0) > 0.5:
            print(f"  !! cell size is {abs(_s.transform.a) * _to_m:.2f} m, not "
                  f"10 m. The extractor and CASCADE assume 10 m (CELL_SIZE_M).")
    print(f"  comparison: {OUT_ROOT}\n")

    folder_names = {domain_id_from(p.parent.name): p.parent.name
                    for p in tifs.values()}

    print(f"{'GIS':>4} {'shape':>11} {'rot':>7} {'cells':>6} {'prof':>5} "
          f"{'cs_med':>7} {'elev':>6} {'nodata':>7}  flag")

    recs, rots = [], []
    for d in sorted(tifs):
        with rasterio.open(tifs[d]) as src:
            mask = burn(road, src)
            t = src.transform
            rot_deg = float(np.degrees(np.arctan2(t.b, t.a)))
            shape_tif = (src.height, src.width)
        rots.append(rot_deg)

        flags, elev = [], None
        npy = find_elev_npy(ELEV_NPY_DIR, d, folder_names.get(d, ""))
        if npy is None:
            flags.append("NO_ELEV_CHECK")
        else:
            elev = np.load(npy).astype(float)
            if elev.shape != shape_tif:
                flags.append(f"SHAPE_MISMATCH{elev.shape}")
                elev = None

        n_cells = int(mask.sum())
        n_prof = int((mask.sum(axis=1) > 0).sum())

        rec = dict(domain=d, section=section_for(d),
                   shape=f"{shape_tif[0]}x{shape_tif[1]}",
                   rot_deg=round(rot_deg, 3),
                   n_cells=n_cells, n_profiles=n_prof, n_nodata=0,
                   road_cs_min=-1, road_cs_median=-1, road_cs_max=-1,
                   road_elev_mhw=np.nan, road_elev_mean_mhw=np.nan,
                   road_elev_std_m=np.nan, road_elev_min_mhw=np.nan,
                   road_elev_max_mhw=np.nan, road_elev_navd=np.nan,
                   land_cs_max=-1, mask_file="", flags="")

        if n_cells == 0:
            flags.append("NO_ROAD")
        else:
            # cross-shore position in the extractor's ocean-first frame
            m_of = ocean_first(mask) > 0
            ci = np.where(m_of)[1]
            rec.update(road_cs_min=int(ci.min()),
                       road_cs_median=int(np.median(ci)),
                       road_cs_max=int(ci.max()))

            if elev is not None:
                z_of = ocean_first(elev) - MHW_M
                raw = z_of[m_of]
                rec["n_nodata"] = int((raw <= NODATA_VALUE - MHW_M + 1e-6).sum())
                good = raw[raw > NODATA_VALUE - MHW_M + 1e-6]
                if good.size:
                    rec.update(
                        road_elev_mhw=round(float(np.median(good)), 3),
                        road_elev_mean_mhw=round(float(np.mean(good)), 3),
                        road_elev_std_m=round(float(np.std(good)), 3),
                        road_elev_min_mhw=round(float(np.min(good)), 3),
                        road_elev_max_mhw=round(float(np.max(good)), 3),
                        road_elev_navd=round(float(np.median(good)) + MHW_M, 3),
                    )
                    if np.median(good) < ROAD_ELEV_MIN_M:
                        flags.append("LOW_ELEV")
                else:
                    flags.append("ALL_NODATA")
                if rec["n_nodata"] > 0.1 * n_cells:
                    flags.append(f"NODATA({rec['n_nodata']})")
                if n_prof < shape_tif[0] * 0.5:
                    flags.append("PATCHY")

                land = np.where((z_of > 0.0).any(axis=0))[0]
                rec["land_cs_max"] = int(land.max()) if land.size else -1

        rec["flags"] = ",".join(flags)
        fn = f"domain_{d}_road_{ROAD_YEAR}.npy"
        np.save(MASK_DIR / fn, mask)
        rec["mask_file"] = fn
        recs.append(rec)

        et = "-" if np.isnan(rec["road_elev_mhw"]) else f"{rec['road_elev_mhw']:.2f}"
        print(f"{d:>4} {rec['shape']:>11} {rot_deg:>7.2f} {n_cells:>6} "
              f"{n_prof:>5} {rec['road_cs_median']:>7} {et:>6} "
              f"{rec['n_nodata']:>7}  {rec['flags']}")

    # --- outputs ---------------------------------------------------------
    print()
    write_diagnostics(recs, DIAG_CSV)
    if SAVE_SUMMARY_FIG:
        summary_figure(recs, SUMMARY_FIG)

    if SAVE_DOMAIN_FIGS:
        want = (QC_DOMAINS if QC_DOMAINS is not None
                else [r["domain"] for r in recs if r["n_cells"] > 0])
        by = {r["domain"]: r for r in recs}
        written, skipped = [], []
        for d in want:
            if d not in tifs:
                skipped.append((d, "no raster"))
                continue
            if by.get(d, {}).get("n_cells", 0) == 0:
                skipped.append((d, "no road in this domain"))
                continue
            npy = find_elev_npy(ELEV_NPY_DIR, d, folder_names.get(d, ""))
            if npy is None:
                skipped.append((d, "elevation array not found"))
                continue
            elev = np.load(npy).astype(float)
            mask = np.load(MASK_DIR / by[d]["mask_file"])
            if elev.shape != mask.shape:
                skipped.append((d, f"shape {elev.shape} vs {mask.shape}"))
                continue
            try:
                domain_figure(d, elev, mask, by[d], FIG_DIR)
                written.append(d)
            except Exception as e:
                skipped.append((d, f"{type(e).__name__}: {e}"))

        print(f"[figures] {FIG_DIR}")
        if written:
            print(f"          wrote {len(written)} of {len(want)}: "
                  + ", ".join(f"domain_{d:03d}_road_mask.png" for d in written))
        if skipped:
            print(f"          SKIPPED {len(skipped)}:")
            for d, why in skipped:
                print(f"            domain {d}: {why}")
        if not written:
            print("          !! nothing was written. Check the reasons above.")

    write_manifest(recs, MANIFEST_PATH)

    # --- verdict ---------------------------------------------------------
    rots = np.array(rots)
    print("\n" + "-" * 82)
    print("GRID ROTATION")
    print("-" * 82)
    if np.allclose(rots, 0.0, atol=0.01):
        print("  All domains rot = 0.00 deg: the grids are NORTH-UP, not")
        print("  rotated to the shoreline. The mask inherits the same")
        print("  transform, so its alignment is still exact.")
        print("  But axis 1 is EAST-WEST, not shore-normal. Where the island")
        print("  trends N-S the error is 1/cos(theta) (~6% at 20 deg); where it")
        print("  bends toward E-W, an E-W profile runs alongshore. This")
        print("  inflates cross-shore distances -- setbacks AND interior")
        print("  widths together -- so it will not show up as OFF_BACK.")
    else:
        print(f"  Rotation varies: {rots.min():+.2f} to {rots.max():+.2f} deg.")
        print("  The grids follow the shoreline; axis 1 is shore-normal.")

    print("\n" + "-" * 82)
    print("REGISTRATION")
    print("-" * 82)
    road_recs = [r for r in recs
                 if FIRST_ROAD_DOMAIN <= r["domain"] <= LAST_ROAD_DOMAIN]
    unver = [r["domain"] for r in road_recs if "NO_ELEV_CHECK" in r["flags"]]
    if unver:
        print(f"  !! NOTHING VERIFIED for {len(unver)} road domain(s): the")
        print(f"     elevation arrays were not found. Fix ELEV_NPY_DIR -- it")
        print(f"     must point at what HAT_dune_topo_extractor.py reads:")
        print(f"       {ELEV_NPY_DIR}")
        return

    hard = [(r["domain"], r["flags"]) for r in road_recs
            if any(k in r["flags"] for k in ("NODATA", "PATCHY", "NO_ROAD",
                                             "SHAPE_MISMATCH", "ALL_NODATA"))]
    low = [r["domain"] for r in road_recs if "LOW_ELEV" in r["flags"]]

    if hard:
        print(f"  {len(hard)} domain(s) with a REAL coverage problem:")
        for d, f in hard:
            print(f"    GIS {d:>2}: {f}")
        print("    NODATA  = road cells outside the domain's own DEM extent.")
        print("              A clip/geojson mismatch, worth fixing.")
        print("    PATCHY  = road on <half the profiles. Raise ROAD_BUFFER_M,")
        print("              or the road genuinely exits the domain.")
    else:
        print("  No coverage problems: every road domain is fully covered.")

    if low:
        print(f"\n  {len(low)} domain(s) flagged LOW_ELEV. THIS IS INFORMATIONAL.")
        print(f"    {low}")
        print("    A ~24 m mask (buffer + all_touched) on a 10 m resampled DEM")
        print("    averages the crown with shoulder and adjacent ground, so a")
        print("    low value may be the DEM's resolution rather than a")
        print("    registration error. Judge registration from the figures:")
        print(f"      {FIG_DIR}")
        print("    Cheap test: set ROAD_BUFFER_M = 0.0 and re-run. If the")
        print("    elevation jumps toward the berm, the buffer was diluting;")
        print("    if it barely moves, the road is not resolved at 10 m.")
        print("    Either way the SETBACK is unaffected -- it comes from the")
        print("    geojson geometry, not the elevation.")

    print(f"\n  Next: HAT_road_offset_from_dune_start.py")
    print(f"    masks in: {MASK_DIR}")


if __name__ == "__main__":
    main()
