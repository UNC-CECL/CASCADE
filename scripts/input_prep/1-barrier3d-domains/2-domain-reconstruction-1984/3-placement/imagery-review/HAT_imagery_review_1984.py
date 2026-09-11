#!/usr/bin/env python3
r"""
HAT_imagery_review_1984.py
==============================================================================
The 1984 footprint against the aerial photographs: for every domain where the
footprint adds or removes rows (and a set of unchanged neighbours as controls),
the same window of the island in the 1984 and 1997 photographs, on the 1 m
lidar the rows are cut from, with the two digitized dune lines, NC-12, and BOTH
candidate placements of the rows drawn on each. Under the panels, the
photograph's brightness along the domain's 50 profiles, so the sand-to-
vegetation transitions can be read as a curve without anyone picking them.

WHAT IT DECIDES - NOTHING. It is a review aid (Hannah's colleagues, 2026-09-08:
"look at the aerial imagery for the domains where we add and remove rows ...
to see how the dunes and interior have actually changed"). The judgement is
made by eye and written into `imagery_review_1984.csv`; this script fills the
numbers, draws the figures, and leaves the verdict columns blank. Re-running it
keeps whatever verdicts are already in the sheet.

THE QUESTION THE FIGURES ARE BUILT TO ANSWER (Hannah, 2026-09-08): did the dune
field actually get narrower between 1984 and 1997, and if so from which side -
so that the rows the footprint adds or removes can be placed where the width
was actually lost or gained. v3 books every change behind NC-12 (rows added)
or directly in front of it (rows removed); the seaward alternative books it at
the dune. Both are drawn on every photograph so the reviewer can say which one
the photographs support, per domain:

    extra_width_was     seaward_of_crest | crest_to_road | behind_road | none
                        where the 1984 island was wider (rows added) or
                        narrower (rows removed) than in 1997
    dune_field_change   narrower | wider | same | unclear
    edge_moved          seaward | landward | both | none
    placement_ok        yes | no | unclear      does v3's placement match
    confidence          high | medium | low
    notes               free text

WHAT IS READ OFF THE PHOTOGRAPHS (decided with Hannah, 2026-09-08)
    * the SEAWARD VEGETATION LINE - bright sand to dark vegetation, the
      clearest edge on a greyscale photograph and roughly the dune toe, the
      feature the digitized lines trace (~19 m seaward of interior row 0).
      Preferred over the wet/dry line because it is less sensitive to the
      beach state on the day: 1984-09-19 is days after Hurricane Diana and
      1997-10-12 is a normal autumn beach.
    * the ROAD CENTRELINE as visible in each photograph, so the crest-to-road
      distance can be judged per year. Both NC-12 alignments are drawn (the
      1984 line is the 1978 export, deliberately; see the road-line README),
      and the reviewer should trust the pavement in the photograph over either.

IMAGERY (D:\Hatteras_GIS\Aerial, the USGS Henderson release, doi
10.5066/P1CXBCDW): georeferenced to the Dare County 2007 orthophotos in NC
State Plane feet, ~1 ft pixels, stated horizontal accuracy 1.2 m for both 1984
and 1997 (RSS of the 2007 control, the scan resolution and the fit). 1984 is
one mosaic in UTM 18N at 0.26 m; 1997 is 32 frames, tiled per domain here, no
mosaic built. Other years in the release (1978-2002) can be asked for with
--years; 1996 is there but Hannah judged it poor (2026-09-08), so the default
pair is 1984 and 1997, the year the second dune line was digitized from.
Anything under the 1.2 m accuracy is not evidence; a Barrier3D cell is 10 m.

FRAME. Everything is drawn in map coordinates (EPSG:3725, the 1 m tiles'
frame), so ALONGSHORE_FLIP does not apply. The model rows are placed on the map
exactly as HAT_verify_road_placement_1984.py places them: map x = interior_x -
row * 10 along each profile, where (interior_x, interior_y) is interior row 0
from RoadOffset_1984_profiles.csv for the 82 road domains, and from the same
cell_to_map chain (re-run through the extractor) for GIS 1-8, which that file
does not cover. On the first road domain met, the re-run is checked against
the CSV to the centimetre, so the two sources cannot silently disagree.

BRIGHTNESS STRIP. For each profile and each 10 m cell along it (from 250 m
seaward of row 0 to the landward edge of the window), the mean pixel value of
the 10 x 10 m block; the strip is the median over the 50 profiles with the
25-75 % band, per year, each year scaled to its own 2-98 % range over the
window so that two films of different exposure can share an axis. Sand is
bright, vegetation and water dark; a step down going landward is the
vegetation line. It is a reading aid, not a measurement: nothing is picked,
nothing is written from it.

CONTROLS. Unchanged domains that border a changed one (both sides of every
add/remove run), thinned to --controls evenly along the island. Same
photographs, same reach, so the eye has a local baseline for "no change".

OUTPUTS  2-domain-reconstruction-1984/3-placement/imagery-review/   (the evidence for the placement step, 2026-09-09)
    imagery_review_1984.csv           the review sheet: numbers filled, verdicts blank
    HAT_imagery_review_1984.txt       the report: sources, rules, domain list
    ../../figures/3-placement/imagery-review/rows-{added,removed}/ , unchanged/
        HAT_imagery_review_GIS<N>.png   one per reviewed domain
    figures/CAPTIONS.md               the caption, one section

USAGE
    python HAT_imagery_review_1984.py                      # 52 changed + 10 controls
    python HAT_imagery_review_1984.py --domains 85,63,84   # a pilot
    python HAT_imagery_review_1984.py --years 1984,1997,1996
    python HAT_imagery_review_1984.py --no-controls
    python HAT_imagery_review_1984.py --resume             # after an interrupted run
==============================================================================
"""
from __future__ import annotations

import argparse
import glob
import warnings
import importlib.util
import re
import sys
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd
import geopandas as gpd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from affine import Affine
import rasterio
from rasterio.warp import reproject, transform_bounds, Resampling
from shapely.geometry import box as _box


def _find_root(start: Path) -> Path:
    for p in [start, *start.parents]:
        if (p / "data" / "hatteras_init").is_dir():
            return p
    raise SystemExit("could not locate the project root")


REPO = _find_root(Path(__file__).resolve())
INIT = REPO / "data" / "hatteras_init"
sys.path.insert(0, str(REPO / "scripts"))
sys.path.insert(0, str(REPO / "scripts" / "input_prep" / "0-elevation" / "3-figures"))
from hat_topo_version import (  # noqa: E402
    duneline_shift_dir, insert_figures_dir, insert_figures_dir_for_domain, rows_sign_sub, insert_scope_step)
import HAT_plot_duneline_offset as off  # noqa: E402   house style, lines, tiles, roads

PRODUCT = "1984-start"
CELL_M = 10.0
ROAD_ROWS = 2
HALF_M = 10.0                 # the geojson is a centreline; the model road is 20 m
SCOPE_DIR = INIT / "1-barrier3d-domains" / PRODUCT / "2-domain-reconstruction-1984"
FOOTPRINT_CSV = insert_scope_step(PRODUCT, "2-extent") / "footprint_1984_by_domain.csv"
SHIFT_DIR = duneline_shift_dir(PRODUCT)
ROAD_DIR = INIT / "4-mgmt-forcing" / "road_offset" / "dunestart_offset" / "1984"
ROAD_OFFSET_SCRIPT = (REPO / "scripts" / "input_prep" / "4-mgmt-forcings" / "road_offset"
                      / "1-produce" / "HAT_road_offset_from_dune_start.py")
STEP_DIR = insert_scope_step(PRODUCT, "3-placement", "imagery-review")        # the sheet and the reports (2026-09-09)
SHEET = STEP_DIR / "imagery_review_1984.csv"
REPORT = STEP_DIR / "HAT_imagery_review_1984.txt"
CAPTIONS = SCOPE_DIR / "figures" / "CAPTIONS.md"

AERIAL_ROOT = Path(r"D:\Hatteras_GIS\Aerial")
AERIAL_DOI = "10.5066/P1CXBCDW"
AERIAL_ACCURACY_M = 1.2       # from the USGS metadata, 1984 and 1997 alike
DEFAULT_YEARS = (1984, 1997)
DISPLAY_RES_M = 0.5           # the photographs are resampled to this for drawing and the strip
SEAWARD_STRIP_M = 250.0       # the strip starts this far seaward of row 0
WINDOW_PAD_SEA_M = 120.0
WINDOW_PAD_LAND_M = 150.0
STRETCH_PCT = (2.0, 98.0)

VERDICT_COLS = ["dune_field_change", "edge_moved", "extra_width_was", "placement_ok",
                "confidence", "notes"]
# Written by the window (HAT_imagery_review_gui.py) from the reviewer's clicks on
# each photograph: three features per year, on the profile nearest the click -
#   toe   the seaward vegetation line (beach sand -> dune), what the dune lines trace
#   back  the landward edge of the dune band (dune sand -> flat vegetated interior)
#   road  the seaward edge of the pavement as it is in that year's photograph
# - each as a position and as metres landward of interior row 0; the bands between
# them per year; their change 1997 - 1984; and where the lost / gained width sat,
# given N: the dune band, the strip from the back of the dune to the road, or the
# remainder, behind the road. Measured by the reviewer, not by code. `band_suggests`
# is DERIVED from the picks and labelled so; the verdict columns stay the reviewer's.
PICK_KINDS = ("toe", "back", "road")
PICK_YEARS = (1984, 1997)
MEASURED_COLS = ([f"{k}{str(y)[2:]}_{f}" for y in PICK_YEARS for k in PICK_KINDS
                  for f in ("x", "y", "profile", "from_row0_m")]
                 + [f"{b}{str(y)[2:]}_m" for y in PICK_YEARS for b in ("dune_band", "back_to_road", "toe_to_road")]
                 + ["toe_shift_m", "back_shift_m", "road_shift_m", "toe_shift_cells",
                    "d_dune_band_m", "d_back_to_road_m", "d_toe_to_road_m",
                    "lost_dune_band_m", "lost_back_to_road_m", "lost_behind_road_m",
                    "band_suggests", "reviewed_by", "reviewed_at"])
# The quick review (HAT_imagery_review_quick.py, 2026-09-10): two yes/no/unclear
# answers per domain - is the 1984 road offset right, is N right - kept across re-runs
# like the columns above.
QUICK_COLS = ["offset_ok", "rows_ok"]
VERDICT_VOCAB = {
    "dune_field_change": "narrower | wider | same | unclear",
    "edge_moved": "seaward | landward | both | none",
    "extra_width_was": "seaward_of_crest | crest_to_road | behind_road | none",
    "placement_ok": "yes | no | unclear",
    "confidence": "high | medium | low",
    "notes": "free text",
}

C_ADD, C_REM, INK = off.C_1984, off.C_1997, off.INK
C_ROAD, C_ROAD_OLD = "#1a1a1a", "0.35"
C_SEA_ALT = "#7b3294"         # the seaward alternative placement, outline only


# =============================================================================
# THE IMAGERY
# =============================================================================

def _year_files(year: int) -> list[Path]:
    """The georeferenced tifs for one year, preferring a finished mosaic.

    The release folders are not uniform: 1984 has a mosaic plus the frames it
    was built from (under 1984_georef_TIF/Input; /zip holds archives), 1996/1997
    are frames only, later years are
    mosaics under other names, and the 2025 draft folder holds every year
    again as frames. One rule: a `<year>_full_aerial.tif` wins; otherwise the
    frames named `<year>_MMDD_*.tif` outside any Input/zip/thumbnail folder.
    """
    mosaics = [Path(p) for p in glob.glob(str(AERIAL_ROOT / f"{year}*" / f"{year}_full_aerial.tif"))]
    frames = []
    for pat in (AERIAL_ROOT / f"{year}_henderson" / "**" / f"{year}_*.tif",
                AERIAL_ROOT / "Henderson_Hatteras2025" / "Hatteras_Georeferenced_DRAFT" / str(year)
                / "**" / f"{year}_*.tif"):
        frames += [Path(p) for p in glob.glob(str(pat), recursive=True)]
    # 1984's frames sit under 1984_georef_TIF/Input: georeferenced, 1 ft, the
    # input to the mosaic, and the only thing that fills the mosaic's gaps
    bad = ("zip", "thumbnail")
    frames = [p for p in frames if p.suffix.lower() == ".tif"
              and not any(b in q.lower() for q in p.parts for b in bad)]
    # the 2025 draft duplicates the henderson folders: keep one copy per name
    seen, out = set(), []
    for p in sorted(frames):
        if p.name not in seen:
            seen.add(p.name)
            out.append(p)
    # the mosaic is read first; the frames only fill where it has no pixels
    # (the 1984 mosaic has gaps, e.g. the south of GIS 63)
    return mosaics[:1] + out


def _year_date(year: int) -> str:
    """The flight date from the frame names (YEAR_MMDD_...), or the year."""
    for pat in (AERIAL_ROOT / f"{year}_henderson" / "**" / f"{year}_*",
                AERIAL_ROOT / "Henderson_Hatteras2025" / "**" / f"{year}_*"):
        for p in glob.glob(str(pat), recursive=True):
            m = re.match(rf"{year}_(\d{{2}})(\d{{2}})_", Path(p).name)
            if m:
                return f"{year}-{m.group(1)}-{m.group(2)}"
    return str(year)


class Imagery:
    """One year's files with their footprints in the map frame."""

    def __init__(self, year: int, dst_crs, strict: bool = True):
        self.year = year
        self.date = _year_date(year)
        self.files = _year_files(year)
        self.foot = []
        self.has_mosaic = False
        if not self.files:
            # the drive is external and can disappear mid-session; the batch
            # script stops, the window carries on from its cache
            if strict:
                raise SystemExit(f"\nno georeferenced imagery for {year} under {AERIAL_ROOT}\n")
            print(f"  WARNING: no imagery for {year} under {AERIAL_ROOT} (drive off?); cached tiles only")
            return
        for p in self.files:
            with rasterio.open(p) as src:
                b = transform_bounds(src.crs, dst_crs, *src.bounds)
                self.foot.append((_box(*b), src.res[0] * (0.3048 if "foot" in str(src.crs).lower()
                                                            or "us_survey_feet" in str(src.crs).lower()
                                                            else 1.0)))
        self.has_mosaic = self.files[0].name.endswith("_full_aerial.tif")
        print(f"  {year} ({self.date}): {len(self.files)} file(s)"
              f"{' (mosaic first, frames fill its gaps)' if self.has_mosaic and len(self.files) > 1 else ''}, "
              f"~{np.median([r for _, r in self.foot]):.2f} m pixels")

    def covering(self, win) -> list[Path]:
        return [p for p, (g, _) in zip(self.files, self.foot) if g.intersects(win)]

    def _edge_distance(self, path: Path):
        """Distance to the nearest no-photograph pixel, in metres, on a coarse
        overview of one file (computed once per file, cached).

        The frames carry a dark FRINGE inside their black border (values ~40
        on 1997_1012_040163d, not 0), so "pixel > 0" cannot tell photograph
        from film edge. Instead every pixel of the tile takes the frame it lies
        farthest inside, the ordinary seamline rule, and a fringe a few metres
        wide can never win against a frame that has real photograph there.
        """
        if not hasattr(self, "_dist"):
            self._dist = {}
        if path in self._dist:
            return self._dist[path]
        from scipy.ndimage import distance_transform_edt
        with rasterio.open(path) as src:
            scale = int(max(1, np.ceil(np.sqrt(src.width * src.height / 4e6))))
            h, w = max(1, src.height // scale), max(1, src.width // scale)
            ov = src.read(out_shape=(src.count, h, w), resampling=Resampling.nearest)
            valid = np.pad(ov.max(axis=0) > 0, 1)
            unit_m = src.res[0] * scale * (0.3048 if "foot" in str(src.crs).lower()
                                           or "us_survey_feet" in str(src.crs).lower() else 1.0)
            dist = (distance_transform_edt(valid)[1:-1, 1:-1] * unit_m).astype(np.float32)
            t = src.transform
            ov_t = Affine(t.a * scale, t.b, t.c, t.d, t.e * scale, t.f)
            self._dist[path] = (dist, ov_t, src.crs)
        return self._dist[path]

    def read(self, x0: float, x1: float, y0: float, y1: float, res: float, dst_crs):
        """(H, W, 3) uint8 in the map frame at `res`, 0 where no photograph.

        Where files overlap, each pixel comes from the file it lies farthest
        inside (see _edge_distance). The mosaic is one file among the others,
        so its gaps are filled by the frames and its own pixels win elsewhere.
        """
        W, H = int(round((x1 - x0) / res)), int(round((y1 - y0) / res))
        dst_t = Affine(res, 0.0, x0, 0.0, -res, y1)
        out = np.zeros((3, H, W), dtype=np.uint8)
        best = np.zeros((H, W), dtype=np.float32)
        win = _box(x0, y0, x1, y1)
        for p in self.covering(win):
            with rasterio.open(p) as src:
                sb = transform_bounds(dst_crs, src.crs, x0, y0, x1, y1)
                w = src.window(*sb).round_offsets().round_lengths()
                try:
                    w = w.intersection(rasterio.windows.Window(0, 0, src.width, src.height))
                except rasterio.errors.WindowError:
                    continue            # footprint touches the box, pixels do not
                if w.width <= 0 or w.height <= 0:
                    continue
                arr = src.read(window=w)
                if arr.shape[0] == 1:
                    arr = np.repeat(arr, 3, axis=0)
                elif arr.shape[0] > 3:
                    arr = arr[:3]
                tmp = np.zeros((3, H, W), dtype=np.uint8)
                st, sc = src.window_transform(w), src.crs
                reproject(arr, tmp, src_transform=st, src_crs=sc, dst_transform=dst_t, dst_crs=dst_crs,
                          resampling=Resampling.average if src.res[0] < res * 0.8 else Resampling.bilinear,
                          num_threads=2)
            dist, ov_t, ov_crs = self._edge_distance(p)
            score = np.zeros((H, W), dtype=np.float32)
            reproject(dist, score, src_transform=ov_t, src_crs=ov_crs, dst_transform=dst_t,
                      dst_crs=dst_crs, resampling=Resampling.bilinear, num_threads=2)
            score[tmp.max(axis=0) == 0] = 0.0
            better = score > best
            out[:, better] = tmp[:, better]
            best[better] = score[better]
        out[:, best <= 0] = 0
        return np.moveaxis(out, 0, -1)


# =============================================================================
# THE PROFILES: interior row 0 on the map, for every domain
# =============================================================================

def _load_module(path: Path, name: str):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


class ProfileFrames:
    """(interior_x, interior_y, row0) per profile, per domain."""

    def __init__(self):
        road = pd.read_csv(ROAD_DIR / "RoadOffset_1984_profiles.csv")
        self.csv = road[["domain", "profile", "interior_row0_cell", "interior_x", "interior_y",
                         "road_seaward_cell", "road_landward_cell"]]
        l84 = pd.read_csv(SHIFT_DIR / "duneline_shift_1984_profiles.csv")
        l97 = pd.read_csv(SHIFT_DIR / "duneline_shift_1997_profiles.csv")
        m = l84.merge(l97, on=["domain", "profile"], suffixes=("_84", "_97"))
        if not (m["interior_row0_cell_84"] == m["interior_row0_cell_97"]).all():
            raise SystemExit("row 0 differs between the 1984 and 1997 profile files")
        self.lines = m.rename(columns={"interior_row0_cell_84": "row0",
                                       "duneline_cell_84": "line84_cell",
                                       "duneline_cell_97": "line97_cell"})[
            ["domain", "profile", "row0", "line84_cell", "line97_cell"]]
        self._ext = None
        self._checked = False

    def _extractor(self):
        if self._ext is None:
            print("  loading the extractor for the domains RoadOffset_1984_profiles.csv does not cover ...")
            ro = _load_module(ROAD_OFFSET_SCRIPT, "hat_road_offset")
            ext = ro.load_extractor(PRODUCT)
            self._ext = (ro, ext, ro.load_windows(ext))
        return self._ext

    def _from_extractor(self, d: int) -> pd.DataFrame:
        ro, ext, windows = self._extractor()
        stem = f"domain_{d}"
        dom = ext.load_profiles(ext.LOAD_PATH / f"{stem}.npy")
        z = dom["z"]
        w = windows.get(stem)
        i0, i1 = (ext.default_window(z, dom["start_beach"]) if w is None
                  else (int(w["i0"]), int(w["i1"])))
        _elev, dune_loc = ext.find_dunes(z, dom["start_beach"], i0, i1)
        row0, _lead = ext.interior_row0_line(z, dune_loc)
        geo = ro.domain_georeference(ext, d)
        if geo is None:
            raise SystemExit(f"no resampled tif for domain {d}: cannot place its profiles on the map")
        rows = []
        for i in range(min(ext.ALONG_COLS, z.shape[0])):
            if row0[i] < 0:
                continue
            x, y = ro.cell_to_map(geo, ext, dom, i, int(row0[i]))
            if x == "":
                continue
            rows.append({"domain": d, "profile": i, "interior_row0_cell": int(row0[i]),
                         "interior_x": x, "interior_y": y,
                         "road_seaward_cell": np.nan, "road_landward_cell": np.nan})
        return pd.DataFrame(rows)

    def get(self, d: int) -> pd.DataFrame:
        pr = self.csv[self.csv["domain"] == d]
        if not len(pr):
            pr = self._from_extractor(d)
            if not self._checked:
                # the same chain must reproduce the CSV on a road domain
                d_chk = int(self.csv["domain"].iloc[0])
                a = self._from_extractor(d_chk).set_index("profile")
                b = self.csv[self.csv["domain"] == d_chk].set_index("profile")
                j = a.join(b, lsuffix="_re", rsuffix="_csv", how="inner")
                dx = np.abs(j["interior_x_re"] - j["interior_x_csv"]).max()
                dy = np.abs(j["interior_y_re"] - j["interior_y_csv"]).max()
                if dx > 0.02 or dy > 0.02 or len(j) < 45:
                    raise SystemExit(f"cell_to_map re-run disagrees with RoadOffset_1984_profiles.csv "
                                     f"at GIS {d_chk}: dx {dx:.2f} m, dy {dy:.2f} m, {len(j)} profiles")
                print(f"  extractor chain reproduces the CSV at GIS {d_chk} to {max(dx, dy):.3f} m")
                self._checked = True
        pr = pr.merge(self.lines[self.lines["domain"] == d], on=["domain", "profile"], how="inner")
        if not (pr["interior_row0_cell"] == pr["row0"]).all():
            raise SystemExit(f"GIS {d}: row 0 differs between the map profiles and the dune-line files")
        pr = pr.drop(columns=["interior_row0_cell"])
        pr["r_line84"] = pr["line84_cell"] - pr["row0"]     # cells relative to row 0, negative = seaward
        pr["r_line97"] = pr["line97_cell"] - pr["row0"]
        return pr.sort_values("interior_y").reset_index(drop=True)


# =============================================================================
# THE ROWS ON THE MAP (as HAT_verify_road_placement_1984.py draws them)
# =============================================================================

def _band_poly(pr: pd.DataFrame, r0: float, r1: float) -> np.ndarray:
    y = pr["interior_y"].to_numpy()
    xa = pr["interior_x"].to_numpy() - r0 * CELL_M + CELL_M / 2
    xb = pr["interior_x"].to_numpy() - r1 * CELL_M - CELL_M / 2
    yy = np.concatenate([[y[0] - CELL_M / 2], y, [y[-1] + CELL_M / 2]])
    xa = np.concatenate([[xa[0]], xa, [xa[-1]]])
    xb = np.concatenate([[xb[0]], xb, [xb[-1]]])
    return np.concatenate([np.column_stack([xa, yy]), np.column_stack([xb[::-1], yy[::-1]])])


def _band(ax, pr, r0, r1, **kw):
    ax.add_patch(plt.Polygon(_band_poly(pr, r0, r1), closed=True, **kw))


def placements(t: pd.Series) -> dict:
    """Row ranges (relative to interior row 0, v2 frame) of everything drawn.

    v3 (the live placement): rows added go directly behind the model's road as
    placed under the 1984 setback, rows removed come out directly in front of
    today's pavement; no-road domains use the crest row. The seaward
    alternative: rows added between the 1984 line and row 0 (drawn as the |N|
    cells seaward of row 0), rows removed as rows 0..|N|-1.
    """
    n = int(t["n_cells"])
    out = {"n": n, "anchor": str(t.get("insert_anchor", "") or ""), "v3": None, "sea": None,
           "pav": None, "road": None}
    if n == 0:
        if np.isfinite(t.get("setback_v2_m", np.nan)):
            r = int(t["setback_v2_m"] // CELL_M)
            out["pav"] = (r, r + 1)
        return out
    ins = int(t["insert_row_behind_road"])
    if out["anchor"] == "road":
        r_v2, r_new = int(t["road_row_v2"]), int(t["road_row_new"])
        out["pav"] = (r_v2, r_v2 + 1)
        if n > 0:
            out["road"] = (r_new, r_new + 1)                 # on the v2 cells it covers
            out["v3"] = (r_new + ROAD_ROWS, r_new + ROAD_ROWS + n - 1)
        else:
            # rows removed in front of the road: v3 row r maps to v2 row r
            # before the seam and r + |N| after it
            def v2_row(r):
                return r if r < ins else r - n
            out["road"] = (v2_row(r_new), v2_row(r_new + 1))
            out["v3"] = (ins, ins - n - 1)
    else:                                                   # crest anchor, GIS 1-5, 8
        out["v3"] = (ins, ins + n - 1) if n > 0 else (ins, ins - n - 1)
    out["sea"] = (-n, -1) if n > 0 else (0, -n - 1)
    return out


# =============================================================================
# THE BRIGHTNESS STRIP
# =============================================================================

def luminance(img: np.ndarray) -> np.ndarray:
    L = img.astype(np.float32).mean(axis=-1)
    L[img.max(axis=-1) == 0] = np.nan
    return L


def strip(L: np.ndarray, x0: float, y1: float, res: float, pr: pd.DataFrame,
          r_min: int, r_max: int) -> tuple[np.ndarray, np.ndarray]:
    """Per (profile, cell) mean brightness of the 10 x 10 m block, and the cells."""
    cells = np.arange(r_min, r_max + 1)
    out = np.full((len(pr), len(cells)), np.nan)
    h = CELL_M / 2
    for i, p in pr.iterrows():
        ra = int(round((y1 - (p["interior_y"] + h)) / res))
        rb = int(round((y1 - (p["interior_y"] - h)) / res))
        for k, r in enumerate(cells):
            xc = p["interior_x"] - r * CELL_M
            ca = int(round((xc - h - x0) / res))
            cb = int(round((xc + h - x0) / res))
            if ca < 0 or cb > L.shape[1] or ra < 0 or rb > L.shape[0]:
                continue
            blk = L[ra:rb, ca:cb]
            ok = np.isfinite(blk)
            if ok.sum() >= 0.5 * blk.size:
                out[i, k] = float(np.nanmean(blk))
    return cells, out


def domain_window(pr: pd.DataFrame, pl: dict, bounds) -> tuple[float, float, float, float]:
    """The window drawn for one domain: from a little seaward of the seaward-most
    dune line to WINDOW_PAD_LAND_M behind the landward-most row drawn, clipped to
    the domain box; the box's full alongshore extent plus 10 m."""
    b = bounds
    x_sea = float((pr["interior_x"] - np.minimum(pr[["r_line84", "r_line97"]].min(axis=1), 0) * CELL_M).max())
    r_land = max([rr[1] for rr in (pl["v3"], pl["pav"], pl["road"], pl["sea"]) if rr is not None] + [12])
    x_hi = min(b[2], x_sea + WINDOW_PAD_SEA_M)
    x_lo = max(b[0], float(pr["interior_x"].min()) - r_land * CELL_M - WINDOW_PAD_LAND_M)
    return x_lo, x_hi, b[1] - 10.0, b[3] + 10.0


# =============================================================================
# THE FIGURE, ONE DOMAIN
# =============================================================================

def fig_domain(d: int, t: pd.Series, pr: pd.DataFrame, imagery: list[Imagery], gdf, lines, roads,
               role: str, with_lidar: bool) -> tuple[Path, dict]:
    off.apply_style()
    n = int(t["n_cells"])
    pl = placements(t)
    gd = gdf[gdf["domain_id"].astype(int) == d]
    box = gd.geometry.iloc[0]
    b = box.bounds
    crs = gdf.crs

    x_lo, x_hi, y_lo, y_hi = domain_window(pr, pl, b)
    x_med = float(pr["interior_x"].median())

    n_img = len(imagery)
    n_pan = n_img + (1 if with_lidar else 0)
    ratio = (y_hi - y_lo) / (x_hi - x_lo)
    w_in = 9.4
    h_map = w_in * ratio * 0.92
    fig = plt.figure(figsize=(w_in + 0.4, h_map * n_pan + 2.6 + 1.1), constrained_layout=True)
    gs = fig.add_gridspec(n_pan + 1, 1, height_ratios=[h_map] * n_pan + [2.4])
    axes = [fig.add_subplot(gs[i, 0]) for i in range(n_pan)]
    axs = fig.add_subplot(gs[n_pan, 0])

    drawn = off.clip_for_drawing(lines, box.buffer(30.0))
    curves = {}
    stats = {}
    r_min = -int(SEAWARD_STRIP_M // CELL_M)
    r_max = int(np.ceil((x_med - x_lo) / CELL_M)) + 1

    def overlays(ax, on_photo: bool):
        gd.boundary.plot(ax=ax, color="0.3", linewidth=0.8, zorder=4)
        off.draw_lines(ax, drawn, scale=1.0, style=off.LINE_STYLE)
        off.m.draw_roads(ax, roads, scale=0.9)
        ax.plot(pr["interior_x"], pr["interior_y"], color="white" if on_photo else INK, lw=1.2,
                ls=(0, (4, 2)), zorder=6)
        if pl["pav"] is not None:
            _band(ax, pr, *pl["pav"], facecolor="none", edgecolor="white" if on_photo else C_ROAD_OLD,
                  lw=1.0, ls=(0, (2, 2)), zorder=7)
        if pl["road"] is not None:
            for rr in pl["road"]:
                _band(ax, pr, rr, rr, facecolor=C_ROAD, edgecolor="none", alpha=0.35, zorder=7)
        if pl["v3"] is not None:
            _band(ax, pr, *pl["v3"], facecolor="none", edgecolor=C_ADD if n > 0 else C_REM, lw=2.0,
                  zorder=8, hatch=None if on_photo else "////")
        if pl["sea"] is not None:
            _band(ax, pr, *pl["sea"], facecolor="none", edgecolor=C_SEA_ALT, lw=1.8, ls=(0, (5, 2)),
                  zorder=8)
        ax.set_xlim(x_lo, x_hi)
        ax.set_ylim(y_lo, y_hi)
        ax.set_aspect("equal")
        ax.set_xticks([])
        ax.set_yticks([])
        ax.text(0.99, 0.5, "ocean", transform=ax.transAxes, ha="right", va="center", fontsize=9,
                color="white" if on_photo else off.INK_MUTED, rotation=90, zorder=9)

    # ---- (a, b, ...) the photographs
    for i, im in enumerate(imagery):
        ax = axes[i]
        img = im.read(x_lo, x_hi, y_lo, y_hi, DISPLAY_RES_M, crs)
        L = luminance(img)
        ok = np.isfinite(L)
        cov = float(ok.mean())
        stats[f"cover_{im.year}"] = round(cov, 3)
        if ok.any():
            lo, hi = np.nanpercentile(L, STRETCH_PCT)
            shown = np.clip((img.astype(np.float32) - lo) / max(hi - lo, 1.0), 0, 1)
            shown[~ok] = 0.93
            ax.imshow(shown, extent=(x_lo, x_hi, y_lo, y_hi), origin="upper", zorder=1,
                      interpolation="bilinear")
            cells, S = strip(L, x_lo, y_hi, DISPLAY_RES_M, pr, r_min, r_max)
            S = (S - lo) / max(hi - lo, 1.0)
            curves[im.year] = (cells, S)
        else:
            ax.set_facecolor("0.93")
            ax.text(0.5, 0.5, f"no {im.year} photograph covers this window", transform=ax.transAxes,
                    ha="center", va="center", fontsize=10, color=off.INK_MUTED)
        overlays(ax, on_photo=True)
        off._scalebar(ax, length_m=100.0)
        if i == 0:
            off._north_arrow(ax, x=0.94, y=0.12)
        src = "mosaic" if im.has_mosaic else f"{len(im.covering(_box(x_lo, y_lo, x_hi, y_hi)))} frame(s)"
        off._title(ax, i, f"GIS {d}: {im.date} aerial photograph ({src}, USGS Henderson release)")

    # ---- (c) the lidar the rows are cut from
    if with_lidar:
        ax = axes[n_img]
        arr, extent = off.load_1m(gdf, [d])
        if arr is not None:
            off._hillshade(ax, arr, extent, res=1.0)
        overlays(ax, on_photo=False)
        off._scalebar(ax, length_m=100.0)
        off._title(ax, n_img, f"GIS {d}: the 1 m lidar the model is built on (1996 foredune, 2009 backdune)")

    # ---- (d) the brightness strip
    for im in imagery:
        if im.year not in curves:
            continue
        cells, S = curves[im.year]
        x = x_med - cells * CELL_M
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", RuntimeWarning)      # all-NaN cells outside the photograph
            med = np.nanmedian(S, axis=0)
            p25, p75 = np.nanpercentile(S, 25, axis=0), np.nanpercentile(S, 75, axis=0)
        col = off.LINE_STYLE[1984]["color"] if im.year == 1984 else (
            off.LINE_STYLE[1997]["color"] if im.year == 1997 else "0.3")
        axs.fill_between(x, p25, p75, color=col, alpha=0.15, lw=0)
        axs.plot(x, med, color=col, lw=1.6, label=f"{im.date}")
    y0s, y1s = -0.05, 1.05
    # the lines and rows, as vertical marks in the same frame as the map
    x_l84 = float((pr["interior_x"] - pr["r_line84"] * CELL_M).median())
    x_l97 = float((pr["interior_x"] - pr["r_line97"] * CELL_M).median())
    axs.axvline(x_l84, color=off.LINE_STYLE[1984]["color"], lw=1.4, zorder=3)
    axs.axvline(x_l97, color=off.LINE_STYLE[1997]["color"], lw=1.4, zorder=3)
    axs.axvline(x_med, color=INK, lw=1.0, ls=(0, (4, 2)), zorder=3)

    def span(rr, **kw):
        if rr is None:
            return
        axs.axvspan(x_med - rr[1] * CELL_M - CELL_M / 2, x_med - rr[0] * CELL_M + CELL_M / 2, **kw)
    span(pl["pav"], facecolor="none", edgecolor=C_ROAD_OLD, lw=1.0, ls=(0, (2, 2)), zorder=2)
    if pl["road"] is not None:
        for rr in pl["road"]:
            span((rr, rr), facecolor=C_ROAD, alpha=0.3, lw=0, zorder=1)
    span(pl["v3"], facecolor=C_ADD if n > 0 else C_REM, alpha=0.18, lw=0, zorder=1)
    span(pl["sea"], facecolor="none", edgecolor=C_SEA_ALT, lw=1.6, ls=(0, (5, 2)), zorder=2)
    axs.set_xlim(x_lo, x_hi)
    axs.set_ylim(y0s, y1s)
    axs.set_ylabel("brightness (2–98 % stretch)")
    axs.set_xlabel("easting (m)")
    sec = axs.secondary_xaxis("top", functions=(lambda xx: (x_med - xx), lambda mm: x_med - mm))
    sec.set_xlabel("m landward of interior row 0 (median profile)")
    axs.grid(True, axis="x", alpha=0.4)
    for sp in ("top",):
        axs.spines[sp].set_visible(True)
    off._title(axs, n_pan, f"GIS {d}: median brightness along the 50 profiles, 10 m cells, 25–75 % band")

    # ---- the legend
    what = ("rows added" if n > 0 else "rows removed") if n else "no rows change"
    where_v3 = {"road": "behind NC-12 as placed" if n > 0 else "directly in front of NC-12",
                "crest": "behind the crest row"}.get(pl["anchor"], "")
    handles = [Line2D([0], [0], **dict(off.LINE_STYLE[1984], linewidth=2.0), label="1984 dune line (digitized)"),
               Line2D([0], [0], **dict(off.LINE_STYLE[1997], linewidth=2.0), label="1997 dune line (digitized)")]
    handles += off.m.road_legend_handles(roads)
    handles += [Line2D([0], [0], color=INK, lw=1.2, ls=(0, (4, 2)), label="interior row 0 (crest + 1)")]
    if pl["pav"] is not None:
        handles.append(Patch(facecolor="none", edgecolor=C_ROAD_OLD, lw=1.0, ls=(0, (2, 2)),
                             label="NC-12 rows as measured on the 1996 surface"))
    if pl["road"] is not None:
        handles.append(Patch(facecolor=C_ROAD, alpha=0.35, label="NC-12 rows as placed in v3 (1984 setback)"))
    if n:
        handles.append(Patch(facecolor="none", edgecolor=C_ADD if n > 0 else C_REM, lw=2.0,
                             label=f"v3: {abs(n)} {what} {where_v3}"))
        handles.append(Patch(facecolor="none", edgecolor=C_SEA_ALT, lw=1.8, ls=(0, (5, 2)),
                             label=f"alternative: the same {abs(n)} rows at the dune (seaward placement)"))
    fig.legend(handles=handles, loc="outside lower center", ncol=3, fontsize=8)

    sub = insert_figures_dir_for_domain(PRODUCT, "3-placement", d, under="imagery-review")
    p = sub / f"HAT_imagery_review_GIS{d}.png"
    fig.savefig(p, dpi=200, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    stats.update({"x_line84": round(x_l84, 1), "x_line97": round(x_l97, 1), "x_row0": round(x_med, 1),
                  "window_w_m": round(x_hi - x_lo, 0)})
    return p, stats


# =============================================================================
# THE SHEET AND THE REPORT
# =============================================================================

def choose_controls(tab: pd.DataFrame, n_controls: int) -> list[int]:
    changed = set(int(d) for d in tab.index[tab["n_cells"] != 0])
    none = [int(d) for d in tab.index if int(tab.loc[d, "n_cells"]) == 0]
    neigh = [d for d in none if (d - 1 in changed) or (d + 1 in changed)]
    if n_controls <= 0:
        return []
    if len(neigh) > n_controls:
        idx = np.unique(np.round(np.linspace(0, len(neigh) - 1, n_controls)).astype(int))
        neigh = [neigh[i] for i in idx]
    return neigh


def write_sheet(rows: list[dict], quiet: bool = False) -> Path:
    new = pd.DataFrame(rows).set_index("domain").sort_index()
    for c in VERDICT_COLS + MEASURED_COLS + QUICK_COLS:
        new[c] = ""
    if SHEET.is_file():
        old = pd.read_csv(SHEET, dtype=str).fillna("")
        if "domain" in old:
            old["domain"] = old["domain"].astype(int)
            old = old.set_index("domain")
            kept = 0
            for c in VERDICT_COLS + MEASURED_COLS + QUICK_COLS:      # the reviewer's work survives a re-run
                if c in old:
                    for d in new.index:
                        if d in old.index and str(old.loc[d, c]).strip():
                            new.loc[d, c] = old.loc[d, c]
                            kept += 1
            # domains reviewed earlier but not in this run stay in the sheet
            extra = [d for d in old.index if d not in new.index]
            if extra:
                new = pd.concat([new, old.loc[extra].reindex(columns=new.columns)]).sort_index()
            if kept and not quiet:
                print(f"  kept {kept} verdict entries already in {SHEET.name}")
    new.to_csv(SHEET)
    return SHEET


def write_report(rows: list[dict], imagery: list[Imagery], controls: list[int], figs: list[Path]) -> Path:
    df = pd.DataFrame(rows).set_index("domain").sort_index()
    ch = df[df["role"] == "changed"]
    L = []
    L.append("HAT_imagery_review_1984.txt - the 1984 footprint against the aerial photographs")
    L.append(f"written {datetime.now():%Y-%m-%d %H:%M} by HAT_imagery_review_1984.py")
    L.append("")
    L.append("WHAT THIS IS")
    L.append("  A review aid, not a measurement. One figure per domain: the same window in each")
    L.append("  photograph and on the 1 m lidar, with the digitized dune lines, NC-12, interior row 0,")
    L.append("  and both placements of the rows the footprint adds or removes (v3: behind NC-12 /")
    L.append("  in front of it; the alternative: at the dune). Under them, the brightness along the")
    L.append("  50 profiles. The verdict is written by eye into imagery_review_1984.csv.")
    L.append("")
    L.append("THE QUESTION (Hannah, 2026-09-08)")
    L.append("  Did the dune field actually get narrower between 1984 and 1997, and from which side?")
    L.append("  Where was the 1984 island wider (rows added) or narrower (rows removed) than in 1997:")
    L.append("  seaward of today's crest, between the crest and the road, or behind the road?")
    L.append("  v3 books every change on the road side; the photographs say whether that is where it was.")
    L.append("")
    L.append("IMAGERY")
    for im in imagery:
        L.append(f"  {im.year}: {im.date}, {len(im.files)} file(s), "
                 f"{'one mosaic' if len(im.files) == 1 else 'frames tiled per domain'}; "
                 f"{im.files[0].parent}")
    L.append(f"  USGS Henderson release, doi {AERIAL_DOI}; georeferenced to the Dare County 2007 orthophotos;")
    L.append(f"  stated horizontal accuracy {AERIAL_ACCURACY_M} m (RSS of control, scan resolution and fit).")
    L.append(f"  Drawn at {DISPLAY_RES_M} m in EPSG:3725. A Barrier3D cell is {CELL_M:.0f} m; offsets under")
    L.append(f"  {AERIAL_ACCURACY_M} m are not evidence of anything.")
    L.append("")
    L.append("WHAT TO READ OFF THE PHOTOGRAPHS (decided 2026-09-08)")
    L.append("  * the seaward vegetation line: bright sand to dark vegetation, ~ the dune toe, the feature")
    L.append("    the digitized lines trace (~19 m seaward of interior row 0). Less sensitive to the")
    L.append("    beach state on the day than the wet/dry line (1984-09-19 is days after Hurricane Diana).")
    L.append("  * the road centreline as visible in each photograph. Both NC-12 geojsons are drawn (the")
    L.append("    1984 line is the 1978 export, deliberately); trust the pavement in the photograph.")
    L.append("")
    L.append("THE VERDICT COLUMNS (imagery_review_1984.csv; blank until filled)")
    for c in VERDICT_COLS:
        L.append(f"  {c:20s} {VERDICT_VOCAB[c]}")
    L.append("")
    L.append("ASSUMPTIONS")
    L.append("  * 1997 stands for 1996: the surface is 1996 ALACE, the line and the photograph are a year later.")
    L.append("  * the vegetation line stands for the dune toe; a bulldozed or planted dune breaks that.")
    L.append("  * the model rows are drawn on the map exactly as HAT_verify_road_placement_1984.py draws")
    L.append("    them (map x = interior_x - row x 10 m along each profile); the seaward alternative is")
    L.append("    the |N| cells seaward of row 0 (rows added) or rows 0..|N|-1 (rows removed).")
    L.append("  * profile coordinates for GIS 1-8 come from re-running the extractor's cell_to_map chain,")
    L.append("    checked against RoadOffset_1984_profiles.csv on a road domain to the centimetre.")
    L.append("  * the brightness strip is scaled per year to its own 2-98 % range over the window; the")
    L.append("    curves are comparable in shape, not in level.")
    L.append("")
    L.append("DOMAINS")
    L.append(f"  changed: {len(ch)} ({int((ch['n_cells'] > 0).sum())} add, {int((ch['n_cells'] < 0).sum())} remove): "
             + ", ".join(str(d) for d in ch.index))
    L.append(f"  controls ({len(controls)}, unchanged neighbours of the changed runs, thinned evenly): "
             + ", ".join(str(d) for d in controls))
    L.append("")
    L.append("PER DOMAIN")
    cols = ["role", "action", "n_cells", "shift_m_median", "shift_m_p10", "shift_m_p90", "insert_anchor",
            "rows_behind_road", "setback_v2_m", "setback_new_m"]
    cols += [c for c in df.columns if c.startswith("cover_")]
    with pd.option_context("display.width", 200, "display.max_rows", 500):
        L.append(df[cols].to_string())
    L.append("")
    L.append(f"FIGURES ({len(figs)})")
    for p in figs:
        L.append(f"  {p.relative_to(SCOPE_DIR)}")
    REPORT.write_text("\n".join(L) + "\n", encoding="utf-8")
    return REPORT


def write_caption(imagery: list[Imagery], n_changed: int, n_controls: int) -> None:
    years = ", ".join(im.date for im in imagery)
    head = "## `HAT_imagery_review_GIS<N>.png` (3-placement/imagery-review/rows-added, rows-removed, unchanged)"
    body = (
        f"The 1984 footprint against the aerial photographs, one figure per reviewed domain "
        f"({n_changed} changed, {n_controls} unchanged neighbours as controls). Panels (a) and (b) "
        f"are the same window of the domain in the {years} photographs of the USGS Henderson release "
        f"(georeferenced to the 2007 orthophotos, stated accuracy {AERIAL_ACCURACY_M} m, drawn at "
        f"{DISPLAY_RES_M} m with a 2–98 % stretch), (c) the 1 m lidar the model is built on (grey relief, "
        f"a 1996 foredune on a 2009 backdune). On all three: the 1984 (red) and 1997 (blue) dune lines "
        f"as digitized, NC-12 in 1984 (dashed) and 2004 (solid), interior row 0 on every profile (thin "
        f"dashed), the road rows as measured on the 1996 surface (dotted outline) and as placed in v3 "
        f"under the 1984 setback (dark fill), and both placements of the rows the footprint changes: "
        f"v3's, behind NC-12 where rows are added or directly in front of it where rows are removed "
        f"(solid outline, red added / blue removed), and the seaward alternative, the same rows at the "
        f"dune (purple dashed outline). Domains without a model road (GIS 1–5, 8) anchor v3's rows on "
        f"the crest row. (d) The photograph's brightness along the domain's 50 profiles: the mean of "
        f"each 10 × 10 m cell from 250 m seaward of row 0 to the landward edge of the window, the median "
        f"over profiles (line) and the 25–75 % band, one curve per year, each scaled to its own 2–98 % "
        f"range so that films of different exposure share an axis; bright is sand, dark is vegetation or "
        f"water, and a step down going landward is the vegetation line. The vertical marks and shading "
        f"repeat the lines and rows of the panels above at the median profile. Nothing is picked from "
        f"the curve; it is a reading aid for the verdict columns of `imagery_review_1984.csv`. Window "
        f"from the beach to about 150 m behind the landward-most row drawn; scale bar 100 m = 10 cells; "
        f"ocean to the right."
    )
    upsert_caption(head, body)


def upsert_caption(head: str, body: str) -> None:
    """Replace or append one `## ...` section of CAPTIONS.md."""
    text = CAPTIONS.read_text(encoding="utf-8") if CAPTIONS.is_file() else "# Figure captions\n"
    if head in text:
        i = text.index(head)
        j = text.find("\n## ", i + len(head))
        text = text[:i] + (text[j + 1:] if j >= 0 else "")
    text = text.rstrip("\n") + "\n\n" + head + "\n\n" + body + "\n"
    CAPTIONS.write_text(text, encoding="utf-8")


# =============================================================================

def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[2])
    ap.add_argument("--domains", default="", help="comma-separated GIS ids (default: every changed domain + controls)")
    ap.add_argument("--years", default=",".join(str(y) for y in DEFAULT_YEARS))
    ap.add_argument("--controls", type=int, default=10, help="how many unchanged neighbours to add (default 10)")
    ap.add_argument("--no-controls", action="store_true")
    ap.add_argument("--no-lidar", action="store_true", help="skip the lidar panel")
    ap.add_argument("--resume", action="store_true",
                    help="skip domains whose figure exists and whose row is in the sheet (reuse that row)")
    a = ap.parse_args()

    if not FOOTPRINT_CSV.is_file():
        raise SystemExit(f"{FOOTPRINT_CSV} not found - run HAT_footprint_1984.py first")
    tab = pd.read_csv(FOOTPRINT_CSV).set_index("domain")
    years = [int(y) for y in a.years.split(",") if y.strip()]

    if a.domains:
        ids = [int(x) for x in a.domains.split(",")]
        controls = [d for d in ids if int(tab.loc[d, "n_cells"]) == 0]
    else:
        controls = [] if a.no_controls else choose_controls(tab, a.controls)
        ids = [int(d) for d in tab.index if int(tab.loc[d, "n_cells"]) != 0] + controls
    ids = sorted(set(ids))
    print(f"{len(ids)} domains: {len(ids) - len(controls)} changed, {len(controls)} controls; years {years}")

    print("loading the map layers ...")
    gdf = off.load_domains()
    lines = off.load_lines(gdf.crs)
    roads = off.m.load_roads(gdf.crs, clip_to=gdf.union_all())
    print("indexing the imagery ...")
    imagery = [Imagery(y, gdf.crs) for y in years]
    frames = ProfileFrames()

    old = (pd.read_csv(SHEET).set_index("domain") if (a.resume and SHEET.is_file()) else None)
    rows, figs = [], []
    for k, d in enumerate(ids, 1):
        t = tab.loc[d]
        role = "control" if d in controls else "changed"
        if old is not None and d in old.index and isinstance(old.loc[d, "figure"], str)                 and (SCOPE_DIR / old.loc[d, "figure"]).is_file():
            print(f"[{k}/{len(ids)}] GIS {d} ({role}) - figure on disk, reusing its sheet row", flush=True)
            figs.append(SCOPE_DIR / old.loc[d, "figure"])
            rec = old.loc[d].to_dict()
            rec["domain"] = d
            rec["role"] = role
            rows.append({k_: v for k_, v in rec.items() if k_ not in VERDICT_COLS + MEASURED_COLS + QUICK_COLS})
            continue
        print(f"[{k}/{len(ids)}] GIS {d} ({role}, N={int(t['n_cells']):+d}) ...", flush=True)
        pr = frames.get(d)
        p, st = fig_domain(d, t, pr, imagery, gdf, lines, roads, role, with_lidar=not a.no_lidar)
        figs.append(p)
        rec = {"domain": d, "role": role, "action": t["action"], "n_cells": int(t["n_cells"]),
               "shift_m_median": t["shift_m_median"], "shift_m_p10": t["shift_m_p10"],
               "shift_m_p90": t["shift_m_p90"], "insert_anchor": t.get("insert_anchor", ""),
               "rows_behind_road": t.get("rows_behind_road", ""),
               "setback_v2_m": t.get("setback_v2_m", np.nan), "setback_new_m": t.get("setback_new_m", np.nan),
               "image_dates": "; ".join(im.date for im in imagery),
               "figure": str(p.relative_to(SCOPE_DIR)).replace("\\", "/")}
        rec.update(st)
        rows.append(rec)
        write_sheet(rows, quiet=True)          # after every domain, so an interrupted run can --resume

    write_sheet(rows)
    write_report(rows, imagery, controls, figs)
    write_caption(imagery, len(ids) - len(controls), len(controls))
    print(f"wrote {SHEET}\nwrote {REPORT}\n{len(figs)} figures under {insert_figures_dir(PRODUCT, '3-placement', 'imagery-review')}")


if __name__ == "__main__":
    main()
