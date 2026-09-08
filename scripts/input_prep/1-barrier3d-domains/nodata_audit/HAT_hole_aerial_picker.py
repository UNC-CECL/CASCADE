"""
HAT_hole_aerial_picker.py

Interactive review of the 1996 aerial chips: reference B, one keystroke per hole.

WHAT YOU ARE DECIDING
---------------------
For each unsurveyed hole where the NCFMP stamp (A) and the blob shape (C)
disagree, one question:

    is there standing water INSIDE the yellow outline, in 1996?

    w   POND      water. The -3.0 m sentinel is right; leave the DEM alone.
    g   DROPOUT   ground. The lidar failed over land; the cell is bridgeable.
    u   UNCLEAR   no vote. The conservative default keeps it water.

UNCLEAR is a real answer and costs nothing. It is there so you are never forced
to guess, which is the failure mode that would quietly turn a review into a
coin-flip dressed as evidence.

Judge only what is inside the outline. Do not try to reconcile it with the two
votes shown in the title - B is worth having precisely because it is an
independent opinion, and a tiebreaker that has already read the other votes is
not one.

THE TRAP IN THE 1996 IMAGERY, AND WHAT TO DO ABOUT IT
------------------------------------------------------
The 1996 frames have a narrow tonal range. Wet sand, damp marsh and shallow
standing water all land on the same mid-grey, and a dark patch is as often
shadow or dense vegetation as it is water. Texture and a closed edge separate
water better than darkness does.

When 1996 will not resolve, press z to cycle the other sources on the drive:

    2019 NGS      0.30 m colour, covers every hole - the most legible
    2008 IOCM     0.50 m colour
    2014          0.35 m colour
    2018 NAIP     0.60 m colour AND near-infrared, plus an NDWI view
    2004          colour, closest in time to 1996 after 1996 itself

Only 2018 carries a verified NIR band; 2014 and 2016 ship four bands but the
fourth is not infrared, so they get no NDWI. See verify_nir().

NDWI is (green - NIR) / (green + NIR): BLUE is water, RED is land. Water
absorbs near-infrared almost completely, so it separates open water from wet
sand in a way no visible-band image can. It is the view to reach for on exactly
the chips that are hard - but it exists only for 2018.

The catch, and it is a real one: only 1996 is contemporaneous with the survey
whose dropouts are in question. Everything else answers "is there a pond here
NOW", which is strong evidence for a pond that has sat in one place for
decades and weak evidence for a marsh pool that migrates. Let the later
imagery break a tie; do not let it overrule a clear 1996 view. The 2018 NIR
also covers only domains 1-8, which happens to be where 48 of the 57 conflicts
are.

KEYS
----
    w / g / u     verdict, then auto-advance
    left / right  move to another hole without deciding
    n             jump to the next hole with no verdict
    z / x         cycle imagery forward / back for this hole
    r             clear this hole's verdict
    q             quit

Every verdict is written to aerial_review.csv IMMEDIATELY, the same way
HAT_dune_topo_extractor.save_windows writes after every domain. Quit whenever
you like and re-run to resume; 58 holes is more than one sitting.

Any write that would REDUCE the number of verdicts on disk copies the old file
to aerial_review.<timestamp>.bak.csv first. That guard exists because the file
was once blanked by a helper that had checked it was empty earlier in the
session and did not re-check before overwriting. A hand-entered review cannot
be regenerated from anything, so it does not get overwritten silently.

CACHE
-----
Rendering a chip means reading a window from up to 33 scanned frames to find
the one with the least black surround, which is slow enough to feel in an
interactive loop. So chips are rendered once into figures/aerial_1996_conflicts
/chip_cache/ as PNGs plus an index of the outline geometry, and the picker
reads those. Delete the folder to force a rebuild. The PNGs are covered by
.gitignore's *.png rule, like every other figure here.

INPUT   dune-topo/<version>/hole_verdicts.csv
        dune-topo/<version>/bracketed_hole_cells.csv
        the 1996 frames, via HAT_hole_aerial_chips

OUTPUT  dune-topo/<version>/figures/aerial_1996_conflicts/aerial_review.csv
            the aerial_verdict column, filled in
"""

import csv
import json
import sys
from collections import defaultdict
from pathlib import Path

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import rasterio
from matplotlib.patches import Polygon
from pyproj import Transformer
from rasterio.windows import from_bounds

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[3]
sys.path.insert(0, str(REPO / "scripts"))
sys.path.insert(0, str(HERE))
import hat_topo_version as htv  # noqa: E402
import HAT_hole_aerial_chips as chips  # noqa: E402

# Chip half-widths are PER HOLE, not fixed. The truncating holes run from 1 to
# 56 cells, so a fixed 90 m view fits a 3-cell hole in D6 and cuts a 560 m hole
# in D1 clean off the edge - which is what the first version did. Each hole gets
# a close view sized to its own footprint and a wide view for context.
ZOOM_MIN_M = 90.0                 # floor, so a single-cell hole is not a pinhole
ZOOM_CLOSE = 0.85                 # close half-width = this x the hole span
ZOOM_WIDE = 2.6                   # wide view = this x the close one

# --- IMAGERY SOURCES ----------------------------------------------------
# 1996 is the only CONTEMPORANEOUS evidence and stays the primary view: it was
# flown the same year as the ALACE survey whose dropouts are in question. The
# rest are CORROBORATION, and they answer a slightly different question - "is
# there a pond here in 2019" rather than "in 1996". For a pond that has sat in
# the same place for decades that is strong support; for a marsh pool that
# migrates it is weaker. Cycle to them when 1996 is ambiguous, which is often,
# because 1996 is a scanned panchromatic-ish frame where wet sand, damp marsh
# and shallow water all land on the same mid-grey.
#
# Near-infrared is the band that actually settles a hard chip. Water absorbs it
# almost completely, so a pond is near-black in NIR while wet sand stays bright
# - the discrimination the visible bands cannot make. A source with real NIR
# gets an extra NDWI view, (green - NIR) / (green + NIR), positive over water.
#
# WHICH BAND IS NIR IS VERIFIED, NOT CONFIGURED. Three of these datasets ship
# four bands, and only one of them is genuinely RGB+NIR:
#
#     2018 NAIP   band 4 veg/water 3.33 against 2.06 for the best visible  ->  NIR
#     2014        band 4 veg/water 0.86 - BRIGHTER over water              ->  not NIR
#     2016        band 4 veg/water 1.16                                    ->  not NIR
#     2019 NGS    band 4 has 12 distinct values                            ->  alpha mask
#
# Assuming band 4 was NIR produced an NDWI panel for 2014 that called open
# water "land" while its own photo showed a pond. So `nir` below is a CANDIDATE
# index, and verify_nir() has to agree before any NDWI view is written.
#
# NOT Google Earth: its imagery is licensed, bulk tile extraction breaches its
# terms, and it is lower resolution here than the 2019 NGS tiles and has no NIR
# band at all. Nothing it offers is missing from this list.
AERIAL_ROOT = Path(r"D:\Hatteras_GIS\Aerial")
SOURCES = [
    dict(key="1996", label="1996 scanned",  res=1.00, zooms=2, nir=None,
         glob="1996_henderson/1996_georef_TIF/*.tif"),
    dict(key="2019", label="2019 NGS",      res=0.30, zooms=1, nir=None,
         glob="2019/*.tif"),
    dict(key="2008", label="2008 IOCM",     res=0.50, zooms=1, nir=None,
         glob="2008/*.tif"),
    dict(key="2014", label="2014 4-band",   res=0.35, zooms=1, nir=4,
         glob="2014/2014_4BandImagery_J1406387/*.tif"),
    dict(key="2018", label="2018 NAIP",     res=0.60, zooms=1, nir=4,
         glob="2018/2018_4BandImagery_NAIP_NC_J1406390/*.tif"),
    dict(key="2004", label="2004 colour",   res=1.46, zooms=1, nir=None,
         glob="2004_Henry/*.tif"),
]
POND, DROPOUT, UNCLEAR = "POND", "DROPOUT", "UNCLEAR"
KEYMAP = {"w": POND, "g": DROPOUT, "u": UNCLEAR}
COLOR = {POND: "#1f6fb4", DROPOUT: "#d7191c", UNCLEAR: "#7a7a7a", "": "#c9c9c9"}

REVIEW_COL = "aerial_verdict"



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
# THE WINDOW HAS TO BE A REAL WINDOW
# =============================================================================

# Backends that deliver key_press_event to a live window. Anything else - Agg,
# or PyCharm's SciView (module://backend_interagg), which renders a figure as a
# static image in a tool pane - makes plt.show() return immediately and the
# script exit with nothing reviewed and no error. That is not a failure mode a
# user should have to diagnose, so it is checked and fixed here.
INTERACTIVE = {"tkagg", "qtagg", "qt5agg", "qt6agg", "wxagg", "macosx",
               "gtk3agg", "gtk4agg", "nbagg", "webagg"}


def ensure_interactive_backend():
    """Switch to a real windowing backend, or stop with instructions."""
    be = matplotlib.get_backend().lower()
    if be in INTERACTIVE:
        return be
    for cand in ("TkAgg", "QtAgg", "Qt5Agg"):
        try:
            plt.switch_backend(cand)
            print(f"[backend] {be} cannot deliver keystrokes; "
                  f"switched to {cand}")
            return cand.lower()
        except Exception:
            continue
    raise SystemExit(
        f"\nmatplotlib is on '{be}', which draws to a file or a static pane "
        f"and never sees a keypress,\nand no interactive backend could be "
        f"loaded in its place.\n\n"
        f"Run this from a terminal rather than inside an IDE console:\n"
        f"    python {Path(__file__).name}\n\n"
        f"In PyCharm, also turn OFF  Settings > Tools > Python Scientific >\n"
        f"'Show plots in tool window' - SciView renders the figure as an image "
        f"and swallows every key.\n")


def free_the_keys():
    """Drop matplotlib's own bindings for the keys this picker uses.

    Found by collision, not by reading the docs: 'g' toggles a grid over the
    imagery, 'r' resets the view, and left/right walk matplotlib's own view
    history. All three fire alongside the picker's handler, so a verdict could
    also silently rescale the chip you were judging.
    """
    for rc in ("keymap.grid", "keymap.grid_minor", "keymap.home",
               "keymap.back", "keymap.forward", "keymap.save",
               "keymap.yscale", "keymap.xscale", "keymap.zoom",
               "keymap.pan", "keymap.fullscreen"):
        matplotlib.rcParams[rc] = []
    # 'q' is left alone: matplotlib's quit closes the window, which is what the
    # picker wants it to do anyway.
    matplotlib.rcParams["keymap.quit"] = ["q", "ctrl+w", "cmd+w"]


# =============================================================================
# CACHE
# =============================================================================

NIR_MIN_RATIO = 2.5        # veg/water brightness ratio a NIR band must beat
NIR_MIN_MARGIN = 1.4       # and by this factor over the best visible band


def verify_nir(path, cand):
    """Is band `cand` really near-infrared? (ok, ratio, best_visible_ratio).

    NIR is bright over vegetation and near-black over water. So classify from
    the VISIBLE bands only - vegetation as green-dominant and above-median
    brightness, water as the darkest decile - then ask which band separates
    them best. A real NIR band wins by a clear margin; a mislabelled fourth
    band does not, and an alpha mask is not even monotonic.
    """
    with rasterio.open(path) as d:
        a = d.read(out_shape=(d.count, 900, 900)).astype(float)
    if a.shape[0] < cand:
        return False, 0.0, 0.0
    vis = a[:3]
    tot = vis.sum(0)
    ok = tot > 0
    if ok.sum() < 1000:
        return False, 0.0, 0.0
    veg = ok & (vis[1] > vis[0]) & (vis[1] > vis[2]) & (tot > np.percentile(tot[ok], 50))
    wat = ok & (tot < np.percentile(tot[ok], 10))
    if veg.sum() < 200 or wat.sum() < 200:
        return False, 0.0, 0.0
    ratio = [a[i][veg].mean() / max(a[i][wat].mean(), 1.0)
             for i in range(a.shape[0])]
    r_cand = ratio[cand - 1]
    r_vis = max(ratio[:3])
    ok = (r_cand >= NIR_MIN_RATIO and r_cand >= NIR_MIN_MARGIN * r_vis
          and r_cand == max(ratio))
    return ok, r_cand, r_vis


class ImagerySource:
    """One aerial dataset, sampled in the domain CRS regardless of its own.

    Each dataset carries its own projection and linear unit - the 1996 frames
    are State Plane in US survey feet, the rest are UTM in metres - so the
    transform and the metres-to-unit factor are read from each file rather
    than assumed anywhere.
    """

    def __init__(self, cfg):
        import glob as _g
        self.__dict__.update(cfg)
        self.files = sorted(_g.glob(str(AERIAL_ROOT / cfg["glob"])))
        self.boxes = []
        self.crs = None
        for p in self.files:
            with rasterio.open(p) as d:
                self.boxes.append((p, d.bounds))
                self.crs = d.crs
        self.ok = bool(self.files)
        self.nir_note = ""
        if self.ok:
            self.tf = Transformer.from_crs(chips.DOMAIN_EPSG, self.crs,
                                           always_xy=True)
            self.unit_m = self.crs.linear_units_factor[1]
            if self.nir:
                good, rc, rv = verify_nir(self.files[0], self.nir)
                if good:
                    self.nir_note = f"  NIR band {self.nir} (veg/water {rc:.2f})"
                else:
                    self.nir_note = (f"  band {self.nir} FAILS the NIR test "
                                     f"(veg/water {rc:.2f} vs {rv:.2f} visible)"
                                     f" - no NDWI")
                    self.nir = None

    def chip(self, x_utm, y_utm, half_m):
        """(rgb, ndwi|None, fx, fy, half, tile) or None if nothing covers it.

        Picks the covering tile with the least black, for the reason in
        HAT_hole_aerial_chips.pick_frame: scanned frames carry an unexposed
        surround baked into the raster, so bounds margin is not a guide.
        """
        if not self.ok:
            return None
        fx, fy = self.tf.transform(x_utm, y_utm)
        half = half_m / self.unit_m
        best = None
        for p, b in self.boxes:
            if not (b.left <= fx <= b.right and b.bottom <= fy <= b.top):
                continue
            with rasterio.open(p) as d:
                win = from_bounds(fx - half, fy - half, fx + half, fy + half,
                                  d.transform)
                a = d.read(window=win, boundless=True, fill_value=0,
                           out_shape=(d.count, chips.CHIP_PX, chips.CHIP_PX))
            rgb = np.moveaxis(a[:3], 0, -1)
            black = float((rgb.max(axis=-1) < 12).mean())
            if best is None or black < best[0]:
                ndwi = None
                if self.nir and a.shape[0] >= self.nir:
                    g = a[1].astype(float)
                    nir = a[self.nir - 1].astype(float)
                    ndwi = (g - nir) / np.maximum(g + nir, 1e-6)
                best = (black, rgb, ndwi, Path(p).name)
            if black < 0.001:
                break
        if best is None:
            return None
        return best[1], best[2], fx, fy, half, best[3]


def build_cache(vdir, cache):
    """Render every conflict hole in every source that covers it.

    Was 1996 only. Extended because the 1996 frames cannot separate wet sand
    from shallow water, which is most of what makes a hole hard to call, and
    the drive already holds colour at 0.30 m and near-infrared at 0.35 m.
    """
    cache.mkdir(parents=True, exist_ok=True)
    verdicts = list(csv.DictReader((vdir / "hole_verdicts.csv").open()))
    conflicts = [r for r in verdicts
                 if "UNKNOWN" not in (r["ncfmp_verdict"], r["shape_verdict"])
                 and r["ncfmp_verdict"] != r["shape_verdict"]]

    cells = defaultdict(list)
    with (vdir / "bracketed_hole_cells.csv").open() as f:
        for r in csv.DictReader(f):
            cells[(int(r["domain"]), int(r["profile"]))].append(
                (int(r["npy_row"]), int(r["npy_col"])))

    origins = {}
    with (REPO / "data" / "hatteras_init" / "0-elevation" / "2009-2014-1996"
          / "2-resampled-10m" / "resample_audit.csv").open() as f:
        for r in csv.DictReader(f):
            origins[int(r["domain"])] = (float(r["origin_x"]),
                                         float(r["origin_y"]))

    srcs = []
    for cfg in SOURCES:
        src = ImagerySource(cfg)
        miss = "" if src.ok else "   MISSING - skipped"
        print(f"[src] {cfg['key']:6s} {cfg['label']:14s} "
              f"{len(src.files):>3} tiles{src.nir_note}{miss}")
        if src.ok:
            srcs.append(src)

    index = []
    total = len(conflicts)
    for i, r in enumerate(sorted(conflicts,
                                 key=lambda r: (int(r["domain"]),
                                                int(r["profile"]))), 1):
        dom, prof = int(r["domain"]), int(r["profile"])
        cs = cells[(dom, prof)]
        polys_utm = [chips.cell_corners_utm(rr, cc, origins[dom])
                     for rr, cc in cs]
        cx = np.mean([p[0] for poly in polys_utm for p in poly])
        cy = np.mean([p[1] for poly in polys_utm for p in poly])
        span = max(
            max(p[0] for poly in polys_utm for p in poly)
            - min(p[0] for poly in polys_utm for p in poly),
            max(p[1] for poly in polys_utm for p in poly)
            - min(p[1] for poly in polys_utm for p in poly))
        close = max(ZOOM_MIN_M, span * ZOOM_CLOSE)

        entry = dict(domain=dom, profile=prof, n_cells=len(cs),
                     span_m=round(span), ncfmp=r["ncfmp_verdict"],
                     shape=r["shape_verdict"], views=[])
        for src in srcs:
            for zi in range(src.zooms):
                zm = close * (ZOOM_WIDE ** zi)
                got = src.chip(cx, cy, zm)
                if got is None:
                    continue
                rgb, ndwi, fx, fy, half, tile = got
                px = []
                for poly in polys_utm:
                    pts = []
                    for x, y in poly:
                        ax, ay = src.tf.transform(x, y)
                        pts.append([(ax - (fx - half)) / (2 * half)
                                    * chips.CHIP_PX,
                                    (1 - (ay - (fy - half)) / (2 * half))
                                    * chips.CHIP_PX])
                    px.append(pts)
                name = f"D{dom}_p{prof}_{src.key}_z{zi}.png"
                plt.imsave(cache / name, rgb)
                entry["views"].append(dict(
                    file=name, source=src.key, label=src.label, kind="photo",
                    half_m=zm, polys=px, frame=tile, black=0.0))
                if ndwi is not None and zi == 0:
                    nm = f"D{dom}_p{prof}_{src.key}_ndwi.png"
                    plt.imsave(cache / nm, np.clip(ndwi, -0.6, 0.6),
                               cmap="RdBu", vmin=-0.6, vmax=0.6)
                    entry["views"].append(dict(
                        file=nm, source=src.key, label=f"{src.label} NDWI",
                        kind="ndwi", half_m=zm, polys=px, frame=tile,
                        black=0.0))
        if entry["views"]:
            index.append(entry)
        print(f"\r  cached {i}/{total}  ({len(entry['views'])} views)   ",
              end="", flush=True)
    print()
    (cache / "index.json").write_text(json.dumps(index), encoding="utf-8")
    return index


# =============================================================================
# REVIEW FILE
# =============================================================================

def load_review(path, index):
    """Existing verdicts, keyed (domain, profile). Missing file is fine."""
    got = {}
    if path.is_file():
        for r in csv.DictReader(path.open()):
            col = next((c for c in r if c.startswith(REVIEW_COL)), None)
            if col:
                got[(int(r["domain"]), int(r["profile"]))] = (r[col] or "").strip()
    return {(e["domain"], e["profile"]): got.get((e["domain"], e["profile"]), "")
            for e in index}


def count_verdicts(path):
    """How many verdicts the file on disk currently holds. 0 if absent."""
    if not path.is_file():
        return 0
    try:
        return sum(1 for r in csv.DictReader(path.open())
                   if any((r[c] or "").strip() for c in r
                          if c.startswith(REVIEW_COL)))
    except Exception:
        return 0


def save_review(path, index, verdicts):
    """Write the review file, backing it up first if this would LOSE verdicts.

    A review pass is hand-entered judgement that cannot be regenerated from
    anything. This file has already been destroyed once - blanked by a helper
    that had checked it was empty earlier in the same session and did not
    re-check before overwriting - so any write that reduces the verdict count
    now leaves a timestamped copy behind and says so.

    The check is on the count rather than on content because that is the only
    thing that matters here: a write that keeps or adds verdicts is the normal
    path, and a write that drops them is either a deliberate reset or a bug,
    and both deserve a copy on disk.
    """
    have = count_verdicts(path)
    want = sum(1 for v in verdicts.values() if v)
    if have > want:
        import shutil
        from datetime import datetime
        bak = path.with_suffix(
            f".{datetime.now():%Y%m%d_%H%M%S}.bak.csv")
        shutil.copy2(path, bak)
        print(f"[backup] this write drops {have} verdicts to {want}; "
              f"previous file copied to {bak.name}")
    with path.open("w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["domain", "profile", "n_cells", "ncfmp_verdict",
                    "shape_verdict", "frame", "black_frac", REVIEW_COL])
        for e in index:
            v = e["views"][0]
            w.writerow([e["domain"], e["profile"], e["n_cells"], e["ncfmp"],
                        e["shape"], v["frame"], v["black"],
                        verdicts[(e["domain"], e["profile"])]])


# =============================================================================
# PICKER
# =============================================================================

def main():
    backend = ensure_interactive_backend()
    free_the_keys()
    topo_dir, _, version = htv.topo_dirs(chips.TOPO_PRODUCT)
    vdir = audit_dir(topo_dir)
    outdir = vdir / "aerial_1996_conflicts"
    cache = outdir / "chip_cache"
    review = outdir / "aerial_review.csv"

    idx_file = cache / "index.json"
    if idx_file.is_file():
        index = json.loads(idx_file.read_text(encoding="utf-8"))
        print(f"[cache] {len(index)} chips from {cache}")
    else:
        print(f"[cache] building into {cache} - this reads the 1996 frames once")
        index = build_cache(vdir, cache)

    verdicts = load_review(review, index)
    state = {"i": 0, "zoom": 0}
    n = len(index)
    imgs = {}

    fig = plt.figure(figsize=(10.5, 10.8))
    ax = fig.add_axes([0.03, 0.10, 0.94, 0.80])
    ax.set_xticks([]); ax.set_yticks([])
    footer = fig.text(0.5, 0.045, "", ha="center", fontsize=10.5,
                      family="monospace")
    tally = fig.text(0.5, 0.015, "", ha="center", fontsize=10, color="0.35",
                     family="monospace")
    try:
        fig.canvas.manager.set_window_title("1996 aerial - pond or dropout")
    except Exception:
        pass

    def draw():
        e = index[state["i"]]
        z = min(state["zoom"], len(e["views"]) - 1)
        v = e["views"][z]
        state["zoom"] = z
        f = cache / v["file"]
        if f not in imgs:
            imgs[f] = plt.imread(f)
        ax.clear()
        ax.imshow(imgs[f], origin="upper",
                  extent=[0, chips.CHIP_PX, chips.CHIP_PX, 0])
        for px in v["polys"]:
            ax.add_patch(Polygon(px, closed=True, fill=False,
                                 edgecolor="#ffe100", lw=2.0))
        # Clip to the image. Without this a polygon reaching past the chip
        # autoscales the axes and the imagery shrinks into a white field.
        ax.set_xlim(0, chips.CHIP_PX)
        ax.set_ylim(chips.CHIP_PX, 0)
        ax.set_aspect("equal")

        # A scale bar of a round length near a quarter of the view, so it stays
        # useful whether the chip is 180 m or 1.5 km across.
        target = v["half_m"] / 2.0
        step = min([1, 2, 5, 10, 20, 25, 50, 100, 200, 250, 500, 1000],
                   key=lambda s: abs(s - target))
        bar = step / (2 * v["half_m"]) * chips.CHIP_PX
        ax.plot([14, 14 + bar], [chips.CHIP_PX - 16] * 2, color="w", lw=4)
        ax.plot([14, 14 + bar], [chips.CHIP_PX - 16] * 2, color="k", lw=1.8)
        ax.text(14, chips.CHIP_PX - 26, f"{step} m", color="w", fontsize=9,
                va="bottom")
        ax.set_xticks([]); ax.set_yticks([])

        key = (e["domain"], e["profile"])
        cur = verdicts[key]
        era = "contemporaneous" if v["source"] == "1996" else "later - context"
        extra = ("   blue = water, red = land"
                 if v["kind"] == "ndwi" else "")
        ax.set_title(
            f"D{e['domain']}  profile {e['profile']}   "
            f"{e['n_cells']} cell{'s' if e['n_cells'] > 1 else ''}, "
            f"{e.get('span_m', 0)} m   [{state['i'] + 1}/{n}]   "
            f"{2 * v['half_m']:.0f} m across\n"
            f"{v['label']}  ({era})   "
            f"[{z + 1}/{len(e['views'])} views]{extra}\n"
            f"NCFMP says {e['ncfmp']}     blob shape says {e['shape']}",
            fontsize=11.5)
        for s in ax.spines.values():
            s.set_color(COLOR[cur]); s.set_linewidth(4 if cur else 1)
        footer.set_text(f"[ {cur or 'no verdict'} ]      "
                        f"w water/POND   g ground/DROPOUT   u unclear   "
                        f"z next view   x prev view   n next undecided   "
                        f"r clear   q quit")
        footer.set_color(COLOR[cur] if cur else "0.2")
        done = sum(1 for x in verdicts.values() if x)
        c = {k: sum(1 for x in verdicts.values() if x == k)
             for k in (POND, DROPOUT, UNCLEAR)}
        tally.set_text(f"{done}/{n} reviewed    "
                       f"POND {c[POND]}   DROPOUT {c[DROPOUT]}   "
                       f"UNCLEAR {c[UNCLEAR]}")
        fig.canvas.draw_idle()

    def advance():
        nxt = next((j for j in range(state["i"] + 1, n)
                    if not verdicts[(index[j]["domain"], index[j]["profile"])]),
                   None)
        state["i"] = nxt if nxt is not None else min(state["i"] + 1, n - 1)

    def on_key(ev):
        e = index[state["i"]]
        key = (e["domain"], e["profile"])
        if ev.key in KEYMAP:
            verdicts[key] = KEYMAP[ev.key]
            save_review(review, index, verdicts)
            advance()
            state["zoom"] = 0
        elif ev.key == "r":
            verdicts[key] = ""
            save_review(review, index, verdicts)
        elif ev.key == "right":
            state["i"] = min(state["i"] + 1, n - 1)
            state["zoom"] = 0
        elif ev.key == "left":
            state["i"] = max(state["i"] - 1, 0)
            state["zoom"] = 0
        elif ev.key == "n":
            advance()
        elif ev.key == "z":
            state["zoom"] = (state["zoom"] + 1) % len(e["views"])
        elif ev.key == "x":
            state["zoom"] = (state["zoom"] - 1) % len(e["views"])
        elif ev.key == "q":
            plt.close(fig)
            return
        else:
            return
        draw()

    fig.canvas.mpl_connect("key_press_event", on_key)
    draw()
    nv = sum(len(e["views"]) for e in index) / max(len(index), 1)
    print(f"\n{n} holes on {backend}, {nv:.0f} views each. "
          f"w=POND  g=DROPOUT  u=UNCLEAR")
    print("z / x cycle imagery: 1996 close, 1996 wide, 2019 colour, 2008, "
          "2014, 2018 NAIP + NDWI, 2004")
    print(f"writing to {review}")
    print("a window should be open now - click it once so it has keyboard "
          "focus, then press w / g / u")
    plt.show(block=True)

    save_review(review, index, verdicts)
    done = sum(1 for x in verdicts.values() if x)
    if done == 0:
        print("\n[!] the window closed with nothing reviewed.")
        print("    If no window ever appeared, plt.show() returned without "
              "blocking, which means")
        print("    this is running inside an IDE console rather than a "
              "terminal. Run it as:")
        print(f"      python {Path(__file__).relative_to(REPO)}")
        print("    and in PyCharm turn off Settings > Tools > Python "
              "Scientific > 'Show plots in tool window'.")
    c = {k: sum(1 for x in verdicts.values() if x == k)
         for k in (POND, DROPOUT, UNCLEAR)}
    print(f"\n{done}/{n} reviewed   POND {c[POND]}   DROPOUT {c[DROPOUT]}   "
          f"UNCLEAR {c[UNCLEAR]}")
    if done < n:
        print("re-run to resume where you stopped")
    else:
        print("complete - re-run HAT_test_hole_pond_or_dropout.py to fold B in")


if __name__ == "__main__":
    main()
