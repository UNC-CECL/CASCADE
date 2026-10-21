"""
The georeferenced aerial photographs on D:, read into the map frame.

    from site_layer.hat_aerial_imagery import Imagery, year_files, AERIAL_ROOT

One year's photographs: a finished mosaic first, with the individual frames
filling its gaps. Moved here on 2026-10-19 from the retired 1984 imagery review
(HAT_imagery_review_1984.py), unchanged, so the mean-shoreline figures keep it.

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-10-19
"""
from __future__ import annotations

import glob
import re
from pathlib import Path

import numpy as np
import rasterio
from affine import Affine
from rasterio.warp import reproject, transform_bounds, Resampling
from shapely.geometry import box as _box

from site_layer.hat_gis_root import GIS_ROOT

AERIAL_ROOT = GIS_ROOT / "aerial_imagery"
AERIAL_DOI = "10.5066/P1CXBCDW"
AERIAL_ACCURACY_M = 1.2       # from the USGS metadata, 1984 and 1997 alike


# The georeferenced tifs for one year, preferring a finished mosaic
def year_files(year: int) -> list[Path]:
    mosaics = [Path(p) for p in glob.glob(str(AERIAL_ROOT / f"{year}*" / f"{year}_full_aerial.tif"))]
    frames = []
    for pat in (AERIAL_ROOT / str(year) / "usgs_henderson" / "**" / f"{year}_*.tif",
                AERIAL_ROOT / "_delivered" / "usgs_henderson_2025_draft" / "Hatteras_Georeferenced_DRAFT"
                / str(year) / "**" / f"{year}_*.tif"):
        frames += [Path(p) for p in glob.glob(str(pat), recursive=True)]
    # The frames sit under <year>/usgs_henderson/frames and fill the mosaic's gaps
    bad = ("zip", "thumbnail")
    frames = [p for p in frames if p.suffix.lower() == ".tif"
              and not any(b in q.lower() for q in p.parts for b in bad)]
    # the 2025 draft duplicates the usgs_henderson folders: keep one copy per name
    seen, out = set(), []
    for p in sorted(frames):
        if p.name not in seen:
            seen.add(p.name)
            out.append(p)
    # The mosaic first; frames fill only where it has no pixels
    return mosaics[:1] + out


# The flight date from the frame names (YEAR_MMDD_...), or the year
def year_date(year: int) -> str:
    for pat in (AERIAL_ROOT / str(year) / "usgs_henderson" / "**" / f"{year}_*",
                AERIAL_ROOT / "_delivered" / "usgs_henderson_2025_draft" / "**" / f"{year}_*"):
        for p in glob.glob(str(pat), recursive=True):
            m = re.match(rf"{year}_(\d{{2}})(\d{{2}})_", Path(p).name)
            if m:
                return f"{year}-{m.group(1)}-{m.group(2)}"
    return str(year)


# One year's files with their footprints in the map frame
class Imagery:

    def __init__(self, year: int, dst_crs, strict: bool = True):
        self.year = year
        self.date = year_date(year)
        self.files = year_files(year)
        self.foot = []
        self.has_mosaic = False
        if not self.files:
            # The drive is external and can disappear mid-session
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

    # Distance to the nearest no-photograph pixel, in m, on a coarse overview of one file (cached)
    def _edge_distance(self, path: Path):
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

    # (H, W, 3) uint8 in the map frame at `res`, 0 where no photograph
    def read(self, x0: float, x1: float, y0: float, y1: float, res: float, dst_crs):
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
                st, sc = src.window_transform(w), src.crs
                # A coarse request reads the frame decimated; a fine one at full resolution
                unit = 0.3048 if ("foot" in str(sc).lower() or "us_survey_feet" in str(sc).lower()) else 1.0
                step = int(res / 2 / (src.res[0] * unit))
                if step >= 2:
                    oh, ow = max(1, int(w.height // step)), max(1, int(w.width // step))
                    arr = src.read(window=w, out_shape=(src.count, oh, ow), resampling=Resampling.average)
                    st = st * Affine.scale(w.width / ow, w.height / oh)
                else:
                    arr = src.read(window=w)
                if arr.shape[0] == 1:
                    arr = np.repeat(arr, 3, axis=0)
                elif arr.shape[0] > 3:
                    arr = arr[:3]
                tmp = np.zeros((3, H, W), dtype=np.uint8)
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
