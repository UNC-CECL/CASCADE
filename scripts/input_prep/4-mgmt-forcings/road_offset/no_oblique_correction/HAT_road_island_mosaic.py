r"""
HAT_road_island_mosaic.py
===============================================================================
Stitch every domain into one island-wide picture with the NC-12 mask on it, so
the road's placement can be judged by eye rather than by a table of medians.

TWO WAYS TO STITCH, AND THEY ANSWER DIFFERENT QUESTIONS
-------------------------------------------------------
STITCH_MODE = "geo"   (default)
    Places every domain at its TRUE easting/northing, read from each
    resampled_*.tif's own affine transform. This is a real map of Hatteras in
    2009 with the road on it. Nothing is reconstructed -- the transform is
    ground truth, and it needs no offset file.

    Use this to answer: is the road where I digitized it, on the island?

STITCH_MODE = "offset"
    Places every domain at Island_Dune_Offsets_<YEAR>, the way
    plot_initialization_poster.py does and the way CASCADE's shoreline_offset
    positions each Barrier3D profile. This is the MODEL's planform, not the
    geography.

    Use this to answer: does the model's island look like the real one?

STITCH_MODE = "both"
    Draws them one above the other. Where they disagree, the offsets are
    telling CASCADE something the DEM does not say. Worth a look given the
    offsets span 6300 m and drive every domain's initial shoreline position.

ORIENTATION
-----------
North runs LEFT -> RIGHT (matching domain 1 -> 131, south -> north). That is a
north-up map rotated 90 deg clockwise, which puts east -- the ocean side -- at
the BOTTOM. Sound at the top.

The cross-shore axis is exaggerated by default: the island is ~65 km long and
~1 km wide, so a true-aspect plot would be an invisible hairline. ASPECT_EXAG
controls it. Set 1.0 for true aspect if you want to be honest about the shape.

WHAT THIS IS FOR
----------------
Registration. The elevation medians under the road mask came back at 0.3-0.9 m
MHW across most of the island, which flagged 40 domains LOW_ELEV. That may just
be a ~24 m mask (buffer + all_touched) on a 10 m resampled DEM averaging the
crown away, rather than a misregistration. A table cannot separate those. A map
can: if the red line traces the island exactly where NC-12 runs, the mask is
fine and the flag is measuring DEM resolution.

REQUIREMENTS
------------
  rasterio, numpy, matplotlib
===============================================================================
"""

import os
import re
from pathlib import Path

import numpy as np
import rasterio
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.colors import FuncNorm

# =============================================================================
# CONFIG
# =============================================================================

PROJECT_ROOT = Path(r"/")
INIT_ROOT = PROJECT_ROOT / "data" / "hatteras_init"

CLIPRESAMPLE_ROOT = Path(
    r"C:\Users\hanna\OneDrive - University of North Carolina at Chapel Hill"
    r"\Ch1_CASCADE_hatteras\Hatteras_CASCADE_Input\domain_elevation"
    r"\2009_domain_clipresample"
)
TIF_GLOB = "resampled_*.tif"

ELEV_NPY_DIR = INIT_ROOT / "elevations" / "2009_pea_hatteras"

ROAD_YEAR = 1984
MASK_DIR = INIT_ROOT / "roads" / "raster" / str(ROAD_YEAR) / "masks"
MASK_FMT = "domain_{d}_road_{year}.npy"

OUT_DIR = INIT_ROOT / "roads" / "raster" / str(ROAD_YEAR)

# --- what to stitch -----------------------------------------------------
STITCH_MODE = "geo"        # "geo" | "offset" | "both"
DOMAIN_RANGE = (1, 131)    # (1, 90) for just the modelled domains

# --- offset mode only ---------------------------------------------------
DUNE_OFFSET_FILE = (INIT_ROOT / "2-brie-offset" / "hindcast_1984"
                    / "Island_Dune_Offsets_1984_PADDED_120.csv")
DUNE_OFFSET_COLUMN = 0
NUM_BUFFER_DOMAINS = 15    # padding at the front of the offset file
FIRST_FILE_NUMBER = 1

# --- geometry -----------------------------------------------------------
CELL_M = 10.0
OCEAN_LOC = "right"        # must match HAT_dune_topo_extractor.py
ALONG_COLS = 50

# --- display ------------------------------------------------------------
MHW_M = 0.36
NODATA_VALUE = -10.0
ELEV_MIN_M = -1.0          # m MHW
ELEV_MAX_M = 4.0
SEA_LEVEL = 0.0
SEA_LEVEL_POS = 0.35       # colormap position of sea level (poster script)

ASPECT_EXAG = 14.0         # cross-shore stretch; 1.0 = true aspect
FIG_W = 26.0

# Blow-ups. Each is (first_domain, last_domain, label).
ZOOMS = [(9, 16, "1999 relocation block + GIS 16"),
         (32, 45, "LOW_ELEV heartland"),
         (78, 87, "NoData cases + 1989 Pea Island relocation")]

# --- style --------------------------------------------------------------
C_TOWN, C_WIMBLE, C_GROIN, C_INK = "#90AFC5", "#E0A800", "#B71C1C", "#1a1a2e"
SECTIONS = [
    ((1, 6), "Cape Point", None),
    ((7, 8), "Buxton", C_TOWN),
    ((21, 31), "Avon", C_TOWN),
    ((32, 67), "Wimble Shoals", C_WIMBLE),
    ((68, 83), "Tri-Village", C_TOWN),
    ((84, 90), "Pea Island", None),
]


# =============================================================================
# LOADING
# =============================================================================

def domain_id_from(text):
    m = re.search(r"(\d+)", str(text))
    return int(m.group(1)) if m else None


def find_domain_tifs(root):
    out = {}
    for sub in sorted(Path(root).iterdir()):
        if not sub.is_dir():
            continue
        d = domain_id_from(sub.name)
        if d is None:
            continue
        tifs = sorted(sub.glob(TIF_GLOB))
        if tifs:
            out[d] = tifs[0]
    return out


def find_elev_npy(d, folder_name):
    for cand in (f"domain_{d}.npy", f"{folder_name}_resampled.npy",
                 f"domain_{d}_resampled.npy", f"{d}_resampled.npy"):
        p = Path(ELEV_NPY_DIR) / cand
        if p.exists():
            return p
    return None


def load_mask(d):
    p = Path(MASK_DIR) / MASK_FMT.format(d=d, year=ROAD_YEAR)
    return np.load(p) if p.exists() else None


def ocean_first(arr):
    if OCEAN_LOC == "right":
        return arr[:, ::-1]
    if OCEAN_LOC == "left":
        return arr
    if OCEAN_LOC == "bottom":
        return np.rot90(arr)[:, ::-1]
    if OCEAN_LOC == "top":
        return np.rot90(arr, -1)[:, ::-1]
    raise ValueError(OCEAN_LOC)


def collect(lo, hi):
    """Elevation, road mask, and the tif's affine, per domain."""
    tifs = find_domain_tifs(CLIPRESAMPLE_ROOT)
    folder = {domain_id_from(p.parent.name): p.parent.name
              for p in tifs.values()}
    out = {}
    for d in sorted(tifs):
        if not (lo <= d <= hi):
            continue
        npy = find_elev_npy(d, folder.get(d, ""))
        if npy is None:
            print(f"  [skip] domain {d}: no elevation array")
            continue
        elev = np.load(npy).astype(float)
        mask = load_mask(d)
        if mask is None:
            mask = np.zeros_like(elev, dtype=np.uint8)
        elif mask.shape != elev.shape:
            print(f"  [skip] domain {d}: mask {mask.shape} vs elev {elev.shape}")
            continue
        with rasterio.open(tifs[d]) as src:
            t = src.transform
        out[d] = dict(elev=elev, mask=mask > 0, transform=t)
    return out


# =============================================================================
# STITCHING
# =============================================================================

def stitch_geo(doms):
    """
    True map. Each domain sits at its own easting/northing, from the tif's
    affine. Canvas rows = easting (west->east), cols = northing (south->north),
    so that plotting puts north to the right and east (the ocean) at the bottom.
    """
    xs, ys = [], []
    for d, v in doms.items():
        t = v["transform"]
        h, w = v["elev"].shape
        # corners; t.a = +cell (east), t.e = -cell (south) for north-up
        for i, j in ((0, 0), (0, w), (h, 0), (h, w)):
            xs.append(t.c + j * t.a + i * t.b)
            ys.append(t.f + j * t.d + i * t.e)
    xmin, xmax, ymin, ymax = min(xs), max(xs), min(ys), max(ys)

    n_e = int(np.ceil((xmax - xmin) / CELL_M)) + 1
    n_n = int(np.ceil((ymax - ymin) / CELL_M)) + 1
    elev = np.full((n_e, n_n), np.nan)
    road = np.zeros((n_e, n_n), dtype=bool)

    for d, v in doms.items():
        t = v["transform"]
        h, w = v["elev"].shape
        z = v["elev"] - MHW_M
        z = np.where(v["elev"] <= NODATA_VALUE + 1e-6, np.nan, z)

        # Round the ORIGIN once, then step by whole cells. Rounding every cell
        # independently is what striped GIS 78-90: a domain whose easting
        # origin sits a half-cell off xmin hits np.round's banker's rounding
        # (round(0.5)=0, round(1.5)=2, round(2.5)=2 ...), giving duplicated and
        # skipped columns in an alternating pattern. Stepping from a single
        # rounded origin is contiguous by construction.
        ri0 = int(round((t.c + 0.5 * t.a - xmin) / CELL_M))
        cj0 = int(round((t.f + 0.5 * t.e - ymin) / CELL_M))
        ris = ri0 + np.arange(w)          # north-up: t.a = +cell
        cjs = cj0 - np.arange(h)          # north-up: t.e = -cell, row 0 = north

        ok_r = (ris >= 0) & (ris < n_e)
        ok_c = (cjs >= 0) & (cjs < n_n)
        if not ok_r.any() or not ok_c.any():
            continue
        rr = ris[ok_r]
        cc = cjs[ok_c]
        block_z = z[np.ix_(ok_c, ok_r)]          # (rows=along/north, cols=east)
        block_m = v["mask"][np.ix_(ok_c, ok_r)]

        elev[np.ix_(rr, cc)] = block_z.T
        road[np.ix_(rr, cc)] |= block_m.T

    return elev, road, dict(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax,
                            doms=doms, mode="geo")


def stitch_offset(doms):
    """
    Model planform. Each domain is placed at Island_Dune_Offsets, the way
    plot_initialization_poster.py and CASCADE's shoreline_offset do it. Canvas
    rows = cross-shore (ocean-first), cols = alongshore (domain order).
    """
    if not Path(DUNE_OFFSET_FILE).exists():
        raise FileNotFoundError(
            f"offset file not found: {DUNE_OFFSET_FILE}\n"
            f"Set STITCH_MODE = 'geo' to stitch from the tif transforms "
            f"instead -- it needs no offset file."
        )
    raw = np.loadtxt(DUNE_OFFSET_FILE, skiprows=1, delimiter=",")
    if raw.ndim == 2:
        raw = raw[:, DUNE_OFFSET_COLUMN]
    offs_cells = np.round(raw / CELL_M).astype(int)

    ids = sorted(doms)
    keep, offs = [], []
    for d in ids:
        pad = NUM_BUFFER_DOMAINS + (d - FIRST_FILE_NUMBER)
        if 0 <= pad < len(offs_cells):
            keep.append(d)
            offs.append(int(offs_cells[pad]))
    if not keep:
        raise RuntimeError(
            f"No domain in {DOMAIN_RANGE} maps into the offset file "
            f"({len(offs_cells)} values). The PADDED_120 file covers GIS 1-90 "
            f"only; set DOMAIN_RANGE = (1, 90)."
        )
    dropped = [d for d in ids if d not in keep]
    if dropped:
        print(f"  [offset] {len(dropped)} domain(s) have no offset "
              f"(outside GIS 1-90): {dropped[0]}-{dropped[-1]}")

    n_cross = max(doms[d]["elev"].shape[1] for d in keep)
    rows = max(offs) + n_cross + 5
    cols = len(keep) * ALONG_COLS
    elev = np.full((rows, cols), np.nan)
    road = np.zeros((rows, cols), dtype=bool)

    for k, d in enumerate(keep):
        v = doms[d]
        z = ocean_first(v["elev"]) - MHW_M
        z = np.where(ocean_first(v["elev"]) <= NODATA_VALUE + 1e-6, np.nan, z)
        m = ocean_first(v["mask"])
        na, nc = z.shape
        c0 = k * ALONG_COLS
        n = min(na, ALONG_COLS)
        o = offs[k]
        elev[o:o + nc, c0:c0 + n] = z[:n, :].T
        road[o:o + nc, c0:c0 + n] = m[:n, :].T

    return elev, road, dict(keep=keep, offs=offs, mode="offset")


# =============================================================================
# PLOTTING
# =============================================================================

def make_norm():
    cmap = plt.cm.terrain.copy()
    cmap.set_bad(color="#b0cfe8")

    def fwd(x):
        r = np.where(x < SEA_LEVEL,
                     SEA_LEVEL_POS * (x - ELEV_MIN_M) / (SEA_LEVEL - ELEV_MIN_M),
                     SEA_LEVEL_POS + (1 - SEA_LEVEL_POS) * x / ELEV_MAX_M)
        return np.where(np.isnan(x), np.nan, r)

    def inv(x):
        return np.where(x < SEA_LEVEL_POS,
                        ELEV_MIN_M + (x / SEA_LEVEL_POS) * (SEA_LEVEL - ELEV_MIN_M),
                        (x - SEA_LEVEL_POS) / (1 - SEA_LEVEL_POS) * ELEV_MAX_M)

    return cmap, FuncNorm((fwd, inv), vmin=ELEV_MIN_M, vmax=ELEV_MAX_M)


def crop_to_data(elev, road):
    """Trim empty margins so the island fills the frame."""
    rows = np.where(~np.isnan(elev).all(axis=1))[0]
    cols = np.where(~np.isnan(elev).all(axis=0))[0]
    if not rows.size or not cols.size:
        return elev, road, 0, 0
    r0, r1 = rows[0], rows[-1] + 1
    c0, c1 = cols[0], cols[-1] + 1
    return elev[r0:r1, c0:c1], road[r0:r1, c0:c1], r0, c0


def domain_tick_positions(meta, n_cols, c0):
    """(tick_col, label) pairs for the alongshore axis."""
    if meta["mode"] == "offset":
        keep = meta["keep"]
        return [(k * ALONG_COLS + ALONG_COLS / 2 - c0, str(d))
                for k, d in enumerate(keep) if d % 5 == 0]
    ticks = []
    for d, v in meta["doms"].items():
        if d % 5:
            continue
        t = v["transform"]
        h = v["elev"].shape[0]
        north = t.f + (h / 2) * t.e
        ticks.append(((north - meta["ymin"]) / CELL_M - c0, str(d)))
    return [(c, l) for c, l in ticks if 0 <= c <= n_cols]


def plot_panel(ax, elev, road, meta, title):
    cmap, norm = make_norm()
    ec, rc, r0, c0 = crop_to_data(elev, road)
    ax.set_facecolor("#b0cfe8")
    im = ax.imshow(np.ma.masked_invalid(ec), cmap=cmap, norm=norm,
                   aspect="auto", origin="upper", interpolation="nearest")
    ri, ci = np.where(rc)
    ax.plot(ci, ri, lw=0, marker="s", ms=0.6, color=C_GROIN, alpha=0.95)

    ticks = domain_tick_positions(meta, ec.shape[1], c0)
    if ticks:
        ax.set_xticks([t[0] for t in ticks])
        ax.set_xticklabels([t[1] for t in ticks], fontsize=8)
    ax.set_yticks([])
    ax.set_title(title, fontsize=11, pad=6)
    return im, ec, rc, c0


def label_sections(ax, meta, ec_cols, c0):
    if meta["mode"] != "offset":
        return
    keep = meta["keep"]
    idx = {d: k for k, d in enumerate(keep)}
    for (lo, hi), label, colr in [(s[0], s[1], s[2]) for s in SECTIONS]:
        if lo not in idx or hi not in idx:
            continue
        x0 = idx[lo] * ALONG_COLS - c0
        x1 = (idx[hi] + 1) * ALONG_COLS - c0
        if colr:
            ax.axvspan(x0, x1, color=colr, alpha=0.16, zorder=0)
        ax.text((x0 + x1) / 2, 0.02, label, transform=ax.get_xaxis_transform(),
                ha="center", va="bottom", fontsize=8,
                bbox=dict(fc="white", ec="none", alpha=0.8, pad=1.5))


def main():
    Path(OUT_DIR).mkdir(parents=True, exist_ok=True)
    lo, hi = DOMAIN_RANGE

    print("=" * 78)
    print(f"ISLAND MOSAIC  |  NC-12 {ROAD_YEAR}  |  mode = {STITCH_MODE}")
    print("=" * 78)
    doms = collect(lo, hi)
    if not doms:
        print("  nothing to stitch")
        return
    n_road = sum(1 for v in doms.values() if v["mask"].any())
    print(f"  {len(doms)} domains (GIS {min(doms)}-{max(doms)}), "
          f"{n_road} with road")

    panels = []
    if STITCH_MODE in ("geo", "both"):
        e, r, m = stitch_geo(doms)
        span_e = (m["xmax"] - m["xmin"]) / 1000
        span_n = (m["ymax"] - m["ymin"]) / 1000
        print(f"  [geo] canvas {e.shape}  |  "
              f"easting span {span_e:.1f} km, northing span {span_n:.1f} km")
        panels.append((e, r, m,
                       f"True map — domains at their tif easting/northing  "
                       f"(2009 DEM, NC-12 {ROAD_YEAR})"))
    if STITCH_MODE in ("offset", "both"):
        try:
            e, r, m = stitch_offset(doms)
            print(f"  [offset] canvas {e.shape}  |  offsets "
                  f"{min(m['offs'])}-{max(m['offs'])} cells "
                  f"({min(m['offs']) * 10}-{max(m['offs']) * 10} m)")
            panels.append((e, r, m,
                           f"Model planform — domains at "
                           f"Island_Dune_Offsets_{ROAD_YEAR}"))
        except (FileNotFoundError, RuntimeError) as exc:
            print(f"  [offset] skipped: {exc}")

    if not panels:
        return

    h_each = max(3.0, FIG_W * (panels[0][0].shape[0] / panels[0][0].shape[1])
                 * ASPECT_EXAG)
    h_each = min(h_each, 7.0)
    fig, axes = plt.subplots(len(panels), 1,
                             figsize=(FIG_W, h_each * len(panels) + 1.4),
                             squeeze=False, facecolor="white")

    im = None
    for k, (e, r, m, title) in enumerate(panels):
        ax = axes[k][0]
        im, ec, rc, c0 = plot_panel(ax, e, r, m, title)
        label_sections(ax, m, ec.shape[1], c0)
        if k == len(panels) - 1:
            ax.set_xlabel("GIS domain  (south -> north;  ocean at the bottom, "
                          "sound at the top)", fontsize=11)

    handles = [mpatches.Patch(color=C_GROIN, label=f"NC-12 {ROAD_YEAR} mask")]
    axes[0][0].legend(handles=handles, loc="upper right", fontsize=9,
                      framealpha=0.92)

    if im is not None:
        cb = fig.colorbar(im, ax=[a[0] for a in axes], location="right",
                          fraction=0.012, pad=0.008)
        cb.set_label("Elevation (m, MHW-relative)", fontsize=10)
        cb.set_ticks([-1, 0, 1, 2, 3, 4])

    fig.suptitle(f"Hatteras Island — NC-12 {ROAD_YEAR} on the CASCADE domains "
                 f"(cross-shore exaggerated {ASPECT_EXAG:g}x)",
                 fontsize=14, fontweight="bold", y=0.995)

    out = Path(OUT_DIR) / f"HAT_road_island_mosaic_{ROAD_YEAR}_{STITCH_MODE}.png"
    fig.savefig(out, dpi=200, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"\n[mosaic] {out}")

    # --- zooms ----------------------------------------------------------
    e, r, m, _ = panels[0]
    cmap, norm = make_norm()
    for z0, z1, label in ZOOMS:
        sub = {d: v for d, v in doms.items() if z0 <= d <= z1}
        if not sub:
            continue
        ze, zr, zm = (stitch_geo(sub) if m["mode"] == "geo"
                      else stitch_offset(sub))
        ec, rc, _, c0 = crop_to_data(ze, zr)
        if not ec.size:
            continue
        fig, ax = plt.subplots(figsize=(15, 6), facecolor="white")
        ax.set_facecolor("#b0cfe8")
        ax.imshow(np.ma.masked_invalid(ec), cmap=cmap, norm=norm,
                  aspect="auto", origin="upper", interpolation="nearest")
        ri, ci = np.where(rc)
        ax.plot(ci, ri, lw=0, marker="s", ms=2.0, color=C_GROIN,
                label=f"NC-12 {ROAD_YEAR}")
        ticks = domain_tick_positions(zm, ec.shape[1], c0)
        if ticks:
            ax.set_xticks([t[0] for t in ticks])
            ax.set_xticklabels([t[1] for t in ticks], fontsize=9)
        ax.set_yticks([])
        ax.legend(loc="upper right", fontsize=9, framealpha=0.9)
        ax.set_title(f"GIS {z0}-{z1}  |  {label}", fontsize=12, pad=6)
        ax.set_xlabel("GIS domain", fontsize=10)
        p = Path(OUT_DIR) / f"HAT_road_zoom_{z0:03d}_{z1:03d}_{ROAD_YEAR}.png"
        fig.savefig(p, dpi=170, bbox_inches="tight", facecolor="white")
        plt.close(fig)
        print(f"[zoom]   {p}")

    print("\n  If the red line traces the island where NC-12 runs, the mask is")
    print("  registered and the LOW_ELEV flags are measuring DEM resolution,")
    print("  not placement. The setback comes from the geojson geometry, so it")
    print("  does not depend on the elevation either way.")


if __name__ == "__main__":
    main()
