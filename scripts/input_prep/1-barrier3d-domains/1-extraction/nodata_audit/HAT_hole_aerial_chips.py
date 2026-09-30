"""
Reference B: 1996 aerial chips for the holes where references A and C disagree.

    python scripts/input_prep/1-barrier3d-domains/1-extraction/nodata_audit/HAT_hole_aerial_chips.py

Reads hole_verdicts.csv; writes contact sheets of 12 chips each and a blank
aerial_review.csv to fill with HAT_hole_aerial_picker.py. Needs the D: drive. Details: scripts/input_prep/1-barrier3d-domains/1-extraction/nodata_audit/README.md.

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-09-18
"""

import csv
import sys
from collections import defaultdict
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import rasterio
from matplotlib.patches import Polygon
from pyproj import Transformer
from rasterio.windows import from_bounds

REPO = next(
    _p for _p in Path(__file__).resolve().parents
    if (_p / "pyproject.toml").exists())   # 1-extraction/nodata_audit/ since 2026-09-09
sys.path.insert(0, str(REPO / "scripts"))
from site_layer import hat_topo_version as htv  # noqa: E402

# --- CONFIG ------------------------------------------------------------------
TOPO_PRODUCT = "1984-start"
AERIAL_DIR = Path(r"D:\Hatteras_GIS\Aerial\1996_henderson\1996_georef_TIF")
DOMAIN_EPSG = 3725
GRID_M = 10.0

CHIP_HALF_M = 90.0        # half-width of a chip, metres on the ground
CHIP_PX = 420             # rendered chip size, pixels
PER_SHEET = 12            # chips per contact sheet (3 x 4)
POND, DROPOUT, UNKNOWN = "POND", "DROPOUT", "UNKNOWN"

# All outputs in one nodata-audit/ folder beside the extraction; audit_dir() names it
AUDIT_SUBDIR = "nodata-audit"
# -----------------------------------------------------------------------------


# <product>/dune-topo/<version>/nodata-audit/, created on demand
def audit_dir(topo_dir):
    d = topo_dir.parent / AUDIT_SUBDIR
    d.mkdir(parents=True, exist_ok=True)
    return d


# [(path, bounds, crs)] for every 1996 frame, read once
def frame_index():
    out = []
    for p in sorted(AERIAL_DIR.glob("*.tif")):
        with rasterio.open(p) as s:
            out.append((p, s.bounds, s.crs))
    if not out:
        raise SystemExit(f"\nno 1996 frames under {AERIAL_DIR}\n")
    print(f"[aerial] {len(out)} frames, {out[0][2].to_string()[:38]}..., "
          f"units {out[0][2].linear_units}")
    return out


# The chip window from one frame, plus its fraction of black pixels
def read_chip(path, fx, fy, half):
    with rasterio.open(path) as s:
        win = from_bounds(fx - half, fy - half, fx + half, fy + half,
                          s.transform)
        img = s.read(window=win, boundless=True, fill_value=0,
                     out_shape=(s.count, CHIP_PX, CHIP_PX))
    img = np.moveaxis(img, 0, -1)
    black = float((img.max(axis=-1) < 12).mean())
    return img, black


# The covering frame that gives the CLEANEST chip, not the widest margin
def pick_frame(frames, fx, fy, half):
    best = (None, None, 1.1)
    for p, b, _ in frames:
        if not (b.left <= fx <= b.right and b.bottom <= fy <= b.top):
            continue
        img, black = read_chip(p, fx, fy, half)
        if black < best[2]:
            best = (p, img, black)
        if black < 0.001:
            break
    return best


# The four UTM corners of one 10 m domain cell
def cell_corners_utm(npy_row, npy_col, origin):
    ox, oy = origin
    x0, x1 = ox + GRID_M * npy_col, ox + GRID_M * (npy_col + 1)
    y0, y1 = oy - GRID_M * (npy_row + 1), oy - GRID_M * npy_row
    return [(x0, y0), (x1, y0), (x1, y1), (x0, y1)]


# Run: the cleanest frame per conflict hole, the contact sheets, the review template
def main():
    # Headless only when this file is the program
    plt.switch_backend("Agg")
    topo_dir, _, version = htv.topo_dirs(TOPO_PRODUCT)
    vdir = audit_dir(topo_dir)
    outdir = vdir / "aerial_1996_conflicts"
    outdir.mkdir(parents=True, exist_ok=True)

    verdicts = list(csv.DictReader((vdir / "hole_verdicts.csv").open()))
    conflicts = [r for r in verdicts
                 if UNKNOWN not in (r["ncfmp_verdict"], r["shape_verdict"])
                 and r["ncfmp_verdict"] != r["shape_verdict"]]
    print(f"{len(conflicts)} conflicting holes of {len(verdicts)}")

    cells = defaultdict(list)
    with (vdir / "bracketed_hole_cells.csv").open() as f:
        for r in csv.DictReader(f):
            cells[(int(r["domain"]), int(r["profile"]))].append(
                (int(r["npy_row"]), int(r["npy_col"])))

    origins = {}
    from site_layer.hat_elevation_products import product as _elprod
    with _elprod("2009-2014-1996").audit_10m.open() as f:
        for r in csv.DictReader(f):
            origins[int(r["domain"])] = (float(r["origin_x"]),
                                         float(r["origin_y"]))

    frames = frame_index()
    tf = Transformer.from_crs(DOMAIN_EPSG, frames[0][2], always_xy=True)
    # metres -> the frame's linear unit, taken from the CRS rather than assumed
    unit_m = frames[0][2].linear_units_factor[1]
    half = CHIP_HALF_M / unit_m

    chips, missing = [], []
    for r in sorted(conflicts, key=lambda r: (int(r["domain"]), int(r["profile"]))):
        dom, prof = int(r["domain"]), int(r["profile"])
        cs = cells[(dom, prof)]
        polys_utm = [cell_corners_utm(rr, cc, origins[dom]) for rr, cc in cs]
        cx = np.mean([p[0] for poly in polys_utm for p in poly])
        cy = np.mean([p[1] for poly in polys_utm for p in poly])
        fx, fy = tf.transform(cx, cy)

        path, img, black = pick_frame(frames, fx, fy, half)
        if path is None:
            missing.append((dom, prof))
            continue

        # cell outlines in chip pixel coordinates
        polys_px = []
        for poly in polys_utm:
            px = []
            for x, y in poly:
                ax, ay = tf.transform(x, y)
                px.append(((ax - (fx - half)) / (2 * half) * CHIP_PX,
                           (1 - (ay - (fy - half)) / (2 * half)) * CHIP_PX))
            polys_px.append(px)

        chips.append(dict(domain=dom, profile=prof, img=img, polys=polys_px,
                          ncfmp=r["ncfmp_verdict"], shape=r["shape_verdict"],
                          n_cells=len(cs), frame=path.name,
                          black_frac=round(black, 3)))

    if missing:
        print(f"  [warn] no 1996 frame covers {len(missing)} holes: {missing}")
    print(f"  rendered {len(chips)} chips")

    # Contact sheets
    sheets = [chips[i:i + PER_SHEET] for i in range(0, len(chips), PER_SHEET)]
    for k, group in enumerate(sheets, 1):
        ncol = 4
        nrow = int(np.ceil(len(group) / ncol))
        fig, axes = plt.subplots(nrow, ncol, figsize=(4.0 * ncol, 4.75 * nrow))
        axes = np.atleast_1d(axes).ravel()
        for ax in axes[len(group):]:
            ax.axis("off")
        for ax, c in zip(axes, group):
            ax.imshow(c["img"], origin="upper",
                      extent=[0, CHIP_PX, CHIP_PX, 0])
            for px in c["polys"]:
                ax.add_patch(Polygon(px, closed=True, fill=False,
                                     edgecolor="#ffe100", lw=1.8))
            ax.set_xticks([]); ax.set_yticks([])
            ax.set_title(f"D{c['domain']}  p{c['profile']}   "
                         f"{c['n_cells']} cell{'s' if c['n_cells'] > 1 else ''}\n"
                         f"NCFMP {c['ncfmp']}   shape {c['shape']}",
                         fontsize=9.5, color="0.15")
            # 50 m scale bar
            bar = 50.0 / (2 * CHIP_HALF_M) * CHIP_PX
            ax.plot([12, 12 + bar], [CHIP_PX - 14] * 2, color="w", lw=3.2)
            ax.plot([12, 12 + bar], [CHIP_PX - 14] * 2, color="k", lw=1.6)
            ax.text(12, CHIP_PX - 22, "50 m", color="w", fontsize=8,
                    va="bottom", path_effects=None)
        doms = sorted({c["domain"] for c in group})
        fig.suptitle(f"1996 aerial - holes where NCFMP and blob shape disagree "
                     f"  |  sheet {k} of {len(sheets)}  |  domains "
                     f"{doms[0]}-{doms[-1]}\n"
                     f"Is there standing water INSIDE the yellow outline?  "
                     f"water = POND (leave it) - ground = DROPOUT (bridge it) - "
                     f"cannot tell = UNCLEAR",
                     fontsize=11.5, y=0.995, va="top")
        top = 1.0 - (0.36 if nrow == 1 else 0.075 / nrow * 3)
        fig.subplots_adjust(left=0.012, right=0.988, top=top, bottom=0.012,
                            wspace=0.05, hspace=0.26)
        p = outdir / f"sheet_{k:02d}_D{doms[0]}-{doms[-1]}.png"
        fig.savefig(p, dpi=150, facecolor="white")
        plt.close(fig)
        print(f"  wrote {p.name}")

    # Review sheet

    # The picker owns this file: an existing one with verdicts is left alone
    rev = outdir / "aerial_review.csv"
    if rev.is_file():
        done = sum(1 for r in csv.DictReader(rev.open())
                   if any((r[c] or "").strip() for c in r
                          if c.startswith("aerial_verdict")))
        if done:
            print(f"\n[keep] {rev.name} already holds {done} verdicts - not "
                  f"overwritten. Delete it to reset.")
            return
    with rev.open("w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["domain", "profile", "n_cells", "ncfmp_verdict",
                    "shape_verdict", "frame", "black_frac",
                    "aerial_verdict  <- POND / DROPOUT / UNCLEAR", "note"])
        for c in chips:
            w.writerow([c["domain"], c["profile"], c["n_cells"], c["ncfmp"],
                        c["shape"], c["frame"], c["black_frac"], "", ""])
    print(f"\nwrote {rev}")
    print(f"  fill the 'aerial_verdict' column, then re-run "
          f"HAT_test_hole_pond_or_dropout.py to fold it in")

    by_dom = defaultdict(int)
    for c in chips:
        by_dom[c["domain"]] += 1
    print("\n  conflicts per domain: "
          + ", ".join(f"D{d} {n}" for d, n in sorted(by_dom.items())))


if __name__ == "__main__":
    main()
