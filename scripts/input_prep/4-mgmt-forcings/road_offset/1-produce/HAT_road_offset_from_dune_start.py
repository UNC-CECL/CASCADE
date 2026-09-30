"""
NC-12 setback and road elevation per domain, measured from interior row 0, the reference CASCADE indexes against.

    python scripts/input_prep/4-mgmt-forcings/road_offset/1-produce/HAT_road_offset_from_dune_start.py

Measures each period's road on its own extraction, floors negatives,
relocates roads that drown at initialisation, and writes the model-facing
RoadSetback CSVs, per-domain and per-profile tables, and an audit. Details: scripts/input_prep/4-mgmt-forcings/road_offset/README.md.

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-09-23
"""
from __future__ import annotations

import csv
import importlib.util
import json
import re
from datetime import datetime
from pathlib import Path

import numpy as np

import sys as _sys
_sys.path.insert(0, str(Path(__file__).resolve().parents[4]))
from site_layer.hat_topo_version import (array_name, topo_dirs,  # noqa: E402
                             YEAR_PRODUCT, road_line_for_year,
                             road_mask_dir, road_mask_file,
                             road_setback_dir)

import matplotlib

# Repo root, found by searching upward
_PATH_REPO = next(_p for _p in Path(__file__).resolve().parents
                  if (_p / "pyproject.toml").exists())

# Only the array helpers are used, never the picker.
matplotlib.use("Agg")


PROJECT_ROOT = Path(str(_PATH_REPO))
INIT_ROOT = PROJECT_ROOT / "data" / "hatteras_init"

# --- CONFIG ------------------------------------------------------------------
# There is now exactly ONE copy of HAT_dune_topo_extractor.py in the repo, and this is it
EXTRACTOR = (PROJECT_ROOT / "scripts" / "input_prep" / "1-barrier3d-domains" / "1-extraction"
             / "HAT_dune_topo_extractor.py")

from site_layer import hat_topo_version as _tv  # noqa: E402
ROADS_ROOT = _tv.ROADS_ROOT
# MASKS ARE KEYED BY LINE VINTAGE, NOT START YEAR (2026-09-15)

# The existing same-year measurements, used as the second reference frame
EXISTING_SETBACK_FMT = _tv.LEGACY_SETBACK_ROOT / "{year}" / "RoadSetback_{year}.csv"

# Offset files, used ONLY to validate delta_m against measured retreat
OFFSET_FMT = None   # superseded by _tv.offset_file(year, "input")

OUT_ROOT = ROADS_ROOT / "dunestart_offset"

YEARS = [1984, 2004]

# EACH YEAR BELONGS TO A TOPOGRAPHY PRODUCT (2026-08-26)
DOMAINS = list(range(1, 91))        # D1 = Cape Point (south) -> D90 = Pea Island (north)

# Road span

# Domains written to the model-facing CSVs
ROAD_SPAN = (9, 90)

# Unstraightened control pass

# A second measurement with STRAIGHTEN = False
CONTROL_UNSTRAIGHTENED = True
# The control picks, repointed after the 2026-08-25 move so nothing moves
CONTROL_WINDOW_JSON = _tv.CONTROL_PICKS_DIR / "HAT_dune_search_windows_2009_pea_hatteras.json"
CONTROL_SUFFIX = "_rawframe"

# Assumed roadway width: one global 20 m; the measured width is a diagnostic
ASSUMED_ROAD_WIDTH_M = 20.0

# Flag thresholds. None of these alter a value; they only label it.
SCATTER_SETBACK_SPREAD_M = 40.0     # p90-p10 of the per-profile setback
SCATTER_ELEV_STD_M = 0.50           # std of road-cell elevation
MIN_ROAD_PROFILES = 25              # of 50; below this the domain is PARTIAL
WATER_FRAC_FLAG = 0.20              # fraction of road cells at sentinel water

# Seaward relocation of drowning roadways

# The second departure from "measure, don't correct"
RELOCATE_DROWNING = True

# bulldoze's own constants, not tuning knobs
DROWN_THRESHOLD_M = 0.0
DROWN_PCT = 0.2
ROAD_WIDTH_CELLS = 2                # int(road_width 20 m / dx 10 m)
# -----------------------------------------------------------------------------


# Extractor import

# Import the corrected extractor, configured for ONE topography product
def load_extractor(product: str | None = None):
    src = EXTRACTOR.read_text(encoding="utf-8")
    name = "hat_extractor"

    if product is not None:
        version = topo_dirs(product)[2]
        # The whole assignment line is replaced, trailing comment and all
        src, n_p = re.subn(r"^TOPO_PRODUCT\s*=.*$",
                           f'TOPO_PRODUCT = "{product}"', src,
                           count=1, flags=re.MULTILINE)
        src, n_v = re.subn(r"^VERSION\s*=.*$",
                           f'VERSION = "{version}"', src,
                           count=1, flags=re.MULTILINE)
        if not (n_p and n_v):
            raise SystemExit(
                f"\ncannot configure {EXTRACTOR.name} for product {product!r}: "
                f"TOPO_PRODUCT matched {n_p} time(s), VERSION {n_v}.\n"
                f"Those two module-level literals are how every path in the "
                f"extractor is derived. If they were renamed or moved into a "
                f"function, this loader has to be updated with them.\n")
        name = f"hat_extractor_{product.replace('-', '_')}"

    spec = importlib.util.spec_from_loader(name, loader=None)
    module = importlib.util.module_from_spec(spec)
    module.__file__ = str(EXTRACTOR)
    _sys.modules[name] = module
    exec(compile(src, str(EXTRACTOR), "exec"), module.__dict__)

    if product is not None and module.TOPO_PRODUCT != product:
        raise SystemExit(
            f"\nextractor configured for {product!r} but reports "
            f"{module.TOPO_PRODUCT!r} after import.\n")

    if not module.ALONGSHORE_FLIP:
        raise SystemExit(
            f"\n{EXTRACTOR} has ALONGSHORE_FLIP = False.\n"
            "The road offsets must be measured in the same alongshore frame as\n"
            "the topography. Either point EXTRACTOR at the corrected copy or\n"
            "set ALONGSHORE_FLIP = True there.\n"
        )
    if not module.STRAIGHTEN:
        print("[warn] extractor has STRAIGHTEN = False; the road mask will not "
              "be sheared, so obliquity stays in the measurement.")
    return module


# Frame alignment now lives in the extractor, called through `ext` (history in readme)


# Geolocation (diagnostics only -- never feeds a setback) Put the measured per-profile cells back on the map

# MOVED TWICE, AND THE SECOND MOVE WAS MISSED UNTIL 2026-08-26
TIF_FMT = _tv.DOMAIN_CLIPS_DIR / "domain_{domain}" / "resampled_domain_{domain}.tif"

WRITE_PROFILE_COORDS = True
_geo_warned: set[str] = set()


# (transform, crs, n_rows, n_cols) for one domain, or None if unavailable
def domain_georeference(ext, domain: int):
    if not WRITE_PROFILE_COORDS:
        return None
    try:
        import rasterio
    except ImportError:
        if "rasterio" not in _geo_warned:
            _geo_warned.add("rasterio")
            print("  [coords] rasterio not installed; profile x/y left blank")
        return None

    # Only OCEAN_LOC = "right" is invertible with the simple chain below
    if ext.OCEAN_LOC != "right":
        if "ocean_loc" not in _geo_warned:
            _geo_warned.add("ocean_loc")
            print(f"  [coords] OCEAN_LOC={ext.OCEAN_LOC!r} is not invertible here; "
                  f"profile x/y left blank")
        return None

    path = Path(str(TIF_FMT).format(domain=domain))
    if not path.is_file():
        if "tif" not in _geo_warned:
            _geo_warned.add("tif")
            print(f"  [coords] no resampled tif under {path.parent.parent}; "
                  f"profile x/y left blank")
        return None
    with rasterio.open(path) as src:
        return src.transform, src.crs, src.shape[0], src.shape[1]


# Aligned-frame (profile, cross-shore cell) -> (x, y) in the tif's CRS
def cell_to_map(geo, ext, dom: dict, profile: int, cross_cell: int):
    if geo is None or cross_cell < 0:
        return "", ""
    transform, _crs, n_rows, n_cols = geo

    j = int(cross_cell) + int(dom["c0"]) + int(dom["shear"][profile])
    col = (n_cols - 1) - j
    row = (n_rows - 1) - profile if ext.ALONGSHORE_FLIP else profile
    if not (0 <= col < n_cols and 0 <= row < n_rows):
        return "", ""

    x, y = transform * (col + 0.5, row + 0.5)
    return round(float(x), 2), round(float(y), 2)


# Per-domain measurement

# Measure one domain's road setback and elevation against its dune start
def measure_domain(ext, domain: int, year: int, windows: dict) -> tuple[dict, list[dict]]:
    stem = f"domain_{domain}"
    dem_path = ext.LOAD_PATH / f"{stem}.npy"
    mask_path = road_mask_file(road_line_for_year(year), domain)

    record = {
        "year": year, "domain": domain,
        "section": ext.section_for(stem),
        "setback_dunestart_m": np.nan,
        "setback_dunestart_floored_m": 0.0,
        # What the model-facing CSV actually carries, and why it differs
        "setback_model_m": np.nan,
        "drowns_at_init": "",
        "relocated_seaward_m": 0.0,
        "n_road_profiles": 0,
        "setback_p10_m": np.nan, "setback_p90_m": np.nan,
        "setback_min_m": np.nan, "setback_max_m": np.nan,
        "setback_center_m": np.nan,
        "measured_road_width_m": np.nan,
        "road_elev_mhw_median": np.nan,
        "road_elev_mhw_mean": np.nan,
        "road_elev_std_m": np.nan,
        "road_elev_navd_median": np.nan,
        "n_road_cells": 0,
        "water_frac": np.nan,
        "obliquity_deg": np.nan,
        "shear_max_cells": 0,
        "lead_trim_rows": 0,
        "window": "",
        "flags": "",
    }
    flags: list[str] = []

    if not dem_path.is_file():
        record["flags"] = "NO_DEM"
        return record, []
    if not mask_path.is_file():
        record["flags"] = "NO_MASK"
        return record, []

    dom = ext.load_profiles(dem_path)
    prof_arr = dom["z"]
    record["obliquity_deg"] = dom["obliquity_deg"]
    record["shear_max_cells"] = int(np.max(dom["shear"]))

    # Same window the topography was extracted with, and the same frame guard the extractor's run pass applies
    w = windows.get(stem)
    if w is None:
        i0, i1 = ext.default_window(prof_arr, dom["start_beach"])
        record["window"] = f"default[{i0},{i1}]"
        flags.append("DEFAULT_WINDOW")
    else:
        if bool(w.get("straightened", False)) != bool(ext.STRAIGHTEN):
            record["flags"] = "WINDOW_FRAME_MISMATCH"
            return record, []
        i0, i1 = int(w["i0"]), int(w["i1"])
        record["window"] = f"[{i0},{i1}]"

    dune_elev, dune_loc = ext.find_dunes(prof_arr, dom["start_beach"], i0, i1)
    if not np.any(dune_loc >= 0):
        record["flags"] = "NO_DUNE"
        return record, []

    row0, lead_trim = ext.interior_row0_line(prof_arr, dune_loc)
    record["lead_trim_rows"] = lead_trim

    mask = ext.align_mask_to_topography(np.load(mask_path), dom)
    if not mask.any():
        record["flags"] = "NO_ROAD"
        return record, []

    n_along = min(ext.ALONG_COLS, prof_arr.shape[0])
    cell = ext.CELL_SIZE_M
    geo = domain_georeference(ext, domain)   # None -> x/y columns stay blank

    profiles: list[dict] = []
    setback_cells: list[float] = []
    center_cells: list[float] = []
    width_cells: list[float] = []
    elev_mhw: list[float] = []
    n_water = 0
    n_cells = 0

    for i in range(n_along):
        road_cells = np.flatnonzero(mask[i])
        if road_cells.size == 0 or row0[i] < 0:
            continue

        # road_setback is the road block's seaward edge, the minimum index ocean-first
        seaward = int(road_cells.min())
        landward = int(road_cells.max())
        center = float(road_cells.mean())

        sb_cells = seaward - int(row0[i])
        setback_cells.append(sb_cells)
        center_cells.append(center - int(row0[i]))
        width_cells.append(landward - seaward + 1)

        profile_elev = prof_arr[i, road_cells]
        wet = ~(profile_elev > ext.SENTINEL_WATER_M + 1e-9)
        n_water += int(wet.sum())
        n_cells += int(road_cells.size)
        elev_mhw.extend(profile_elev[~wet].tolist())

        # Map coordinates of the two cells the setback is measured between
        interior_x, interior_y = cell_to_map(geo, ext, dom, i, int(row0[i]))
        road_x, road_y = cell_to_map(geo, ext, dom, i, seaward)

        profiles.append({
            "year": year, "domain": domain, "profile": i,
            "dune_crest_cell": int(dune_loc[i]),
            "interior_row0_cell": int(row0[i]),
            "road_seaward_cell": seaward,
            "road_landward_cell": landward,
            "road_center_cell": round(center, 2),
            "setback_m": round(sb_cells * cell, 1),
            "road_width_m": round((landward - seaward + 1) * cell, 1),
            "n_road_cells": int(road_cells.size),
            "n_wet_road_cells": int(wet.sum()),
            "interior_x": interior_x, "interior_y": interior_y,
            "road_x": road_x, "road_y": road_y,
        })

    if not setback_cells:
        record["flags"] = "NO_ROAD"
        return record, []

    sb = np.asarray(setback_cells, dtype=float) * cell
    record["n_road_profiles"] = len(sb)
    record["setback_dunestart_m"] = round(float(np.median(sb)), 1)
    record["setback_dunestart_floored_m"] = round(max(float(np.median(sb)), 0.0), 1)
    record["setback_p10_m"] = round(float(np.percentile(sb, 10)), 1)
    record["setback_p90_m"] = round(float(np.percentile(sb, 90)), 1)
    record["setback_min_m"] = round(float(sb.min()), 1)
    record["setback_max_m"] = round(float(sb.max()), 1)
    record["setback_center_m"] = round(
        float(np.median(np.asarray(center_cells) * cell)), 1)
    record["measured_road_width_m"] = round(
        float(np.median(np.asarray(width_cells) * cell)), 1)

    record["n_road_cells"] = n_cells
    record["water_frac"] = round(n_water / n_cells, 3) if n_cells else np.nan
    if elev_mhw:
        arr = np.asarray(elev_mhw, dtype=float)
        record["road_elev_mhw_median"] = round(float(np.median(arr)), 3)
        record["road_elev_mhw_mean"] = round(float(arr.mean()), 3)
        record["road_elev_std_m"] = round(float(arr.std()), 3)
        record["road_elev_navd_median"] = round(
            float(np.median(arr)) + ext.MHW_M, 3)
    else:
        flags.append("ALL_ROAD_CELLS_WET")

    # Flags: labels only, no value is altered
    if record["setback_dunestart_m"] < 0:
        flags.append(f"NEGATIVE({record['setback_dunestart_m']:.0f}->floored 0)")
    if record["n_road_profiles"] < MIN_ROAD_PROFILES:
        flags.append(f"PARTIAL({record['n_road_profiles']}/{n_along})")
    spread = record["setback_p90_m"] - record["setback_p10_m"]
    if spread > SCATTER_SETBACK_SPREAD_M:
        flags.append(f"SCATTER_SETBACK({spread:.0f}m)")
    if (np.isfinite(record["road_elev_std_m"])
            and record["road_elev_std_m"] > SCATTER_ELEV_STD_M):
        flags.append(f"SCATTER_ELEV({record['road_elev_std_m']:.2f})")
    if np.isfinite(record["water_frac"]) and record["water_frac"] > WATER_FRAC_FLAG:
        flags.append(f"IN_WATER({record['water_frac']:.0%})")
    if abs(record["measured_road_width_m"] - ASSUMED_ROAD_WIDTH_M) > 20.0:
        flags.append(f"WIDTH({record['measured_road_width_m']:.0f}m)")

    record["flags"] = ",".join(flags)
    return record, profiles


# Seaward relocation -- see RELOCATE_DROWNING for the decision and its cost

# The interior array CASCADE actually initialises with
def load_saved_interior(ext, domain: int) -> np.ndarray | None:
    # Array names from the resolver: no year tag since 2026-08-26
    p = Path(ext.TOPO_SAVE_PATH) / array_name("topography", domain)
    return np.load(p) if p.is_file() else None


# bulldoze's width-drown test, transcribed
def drown_test(interior: np.ndarray, road_start: int) -> tuple | None:
    n = interior.shape[0]
    end = road_start + ROAD_WIDTH_CELLS
    if road_start < 0 or end + 1 >= n:
        return None
    wet = lambda row: float((row * 10.0 <= DROWN_THRESHOLD_M).mean())
    seaside = wet(interior[road_start - 1, :]) if road_start > 0 else 0.0
    bayside = wet(interior[end + 1, :])
    return seaside, bayside, bool(seaside > DROWN_PCT or bayside > DROWN_PCT)


# Largest road_start strictly seaward of the current one that does not drown
def nearest_viable_seaward(interior: np.ndarray, road_start: int) -> int | None:
    for start in range(road_start - 1, -1, -1):
        t = drown_test(interior, start)
        if t is not None and not t[2]:
            return start
    return None


# Fill `setback_model_m` for every in-span domain, moving the drowned ones
def relocate_drowning(ext, in_span: list[dict], cell: float) -> None:
    for r in in_span:
        floored = float(r["setback_dunestart_floored_m"])
        r["setback_model_m"] = round(floored, 1)

        interior = load_saved_interior(ext, r["domain"])
        if interior is None:
            r["flags"] = ",".join(filter(None, [r["flags"], "NO_SAVED_TOPO"]))
            continue

        start = int(floored / cell)
        t = drown_test(interior, start)
        if t is None:
            # Past the end of the array is an overrun, not a drowning; the audit judges it
            r["flags"] = ",".join(filter(None, [r["flags"], "OVERRUN"]))
            continue

        sea, bay, drowns = t
        r["drowns_at_init"] = "yes" if drowns else "no"
        if not drowns or not RELOCATE_DROWNING:
            continue

        best = nearest_viable_seaward(interior, start)
        if best is None:
            r["flags"] = ",".join(filter(None, [
                r["flags"],
                f"DROWNS_NO_VIABLE_ROW(sea{sea:.0%},bay{bay:.0%})"]))
            continue

        moved_m = (start - best) * cell
        r["setback_model_m"] = round(best * cell, 1)
        r["relocated_seaward_m"] = round(moved_m, 1)
        note = [f"MOVED_SEAWARD({moved_m:.0f}m,sea{sea:.0%},bay{bay:.0%})"]

        # A move inside the domain's own per-profile spread is a re-pick of the alongshore statistic
        if (np.isfinite(r["setback_min_m"])
                and r["setback_model_m"] < r["setback_min_m"]):
            note.append(f"BEYOND_MEASURED(min {r['setback_min_m']:.0f}m)")
        r["flags"] = ",".join(filter(None, [r["flags"]] + note))


# Existing same-year setbacks + offset validation

# Read a 2-row (GIS IDs, values) CASCADE forcing file
def read_two_row_csv(path: Path) -> dict[int, float]:
    if not path.is_file():
        return {}
    raw = np.loadtxt(path, delimiter=",")
    if raw.ndim != 2 or raw.shape[0] != 2:
        raise ValueError(f"{path}: expected 2 rows, got shape {raw.shape}")
    return {int(k): float(v) for k, v in zip(raw[0], raw[1])}


# Write the 2-row format load_padded_series expects
def write_two_row_csv(path: Path, values: dict[int, float]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    ids = sorted(values)
    with open(path, "w", newline="") as f:
        f.write(",".join(f"{i}.000" for i in ids) + "\n")
        f.write(",".join(f"{values[i]:.3f}" for i in ids) + "\n")


# Per-domain dune-line station from the SHARED offshore datum, metres, LANDWARD-positive
def load_stations(year: int) -> np.ndarray | None:
    path = _tv.dune_raw_file_for_year(year, strict=False)
    if path is None or not path.is_file():
        return None
    # First row per (domain, transect), then the mean of the transects in each domain
    seen, sums, counts = set(), {}, {}
    with open(path, newline="", encoding="utf-8-sig") as fh:
        for row in csv.DictReader(fh):
            try:
                dom, line = int(float(row["domain_id"])), int(float(row["LineID"]))
                station = float(row["ORIG_LEN"])
            except (KeyError, TypeError, ValueError):
                continue
            if (dom, line) in seen:
                continue
            seen.add((dom, line))
            sums[dom] = sums.get(dom, 0.0) + station
            counts[dom] = counts.get(dom, 0) + 1
    if not all(d in counts for d in DOMAINS):
        return None
    return np.array([sums[d] / counts[d] for d in DOMAINS], dtype=float)


# The dune-search windows for the product/version an extractor is set to
def load_windows(ext) -> dict:
    if ext.WINDOW_JSON.is_file():
        return json.load(open(ext.WINDOW_JSON))

    picks_dir = ext.WINDOW_JSON.parent
    present = (sorted(p.name for p in picks_dir.glob("HAT_dune_search_windows_*.json"))
               if picks_dir.is_dir() else [])
    raise SystemExit(
        f"\nno dune-search picks for {ext.TOPO_PRODUCT}/{ext.VERSION}.\n"
        f"  looked for : {ext.WINDOW_JSON}\n"
        f"  present    : {', '.join(present) if present else '(none)'}\n\n"
        f"The window file is named from the extractor's VERSION literal, so a\n"
        f"version bump that did not write one lands here. The picks are the only\n"
        f"input to this script that cannot be regenerated, so DO NOT re-pick to\n"
        f"get past this -- seed the new version from the one it derives from:\n\n"
        f"    python -c \"import importlib.util as u; "
        f"s=u.spec_from_file_location('b', r'scripts/input_prep/"
        f"1-barrier3d-domains/1-extraction/nodata_audit/HAT_bridge_dropouts.py'); "
        f"m=u.module_from_spec(s); s.loader.exec_module(m); "
        f"m.carry_picks_forward('<source version>', '{ext.VERSION}')\"\n\n"
        f"That copies the source version's windows, stamps _meta.inherited_from\n"
        f"so the copy is not mistaken for a re-pick, and leaves the source file\n"
        f"untouched. Only correct when the new version did not move any dune --\n"
        f"a bridged or otherwise value-only version. If the new version really\n"
        f"was re-picked, run the extractor's pick pass for it instead.\n")


# Run: every period's extraction, measure, floor, relocate, write the CSVs and audit
def main() -> None:
    # ONE EXTRACTOR PER PRODUCT, both live for the whole run
    exts = {year: load_extractor(YEAR_PRODUCT[year]) for year in YEARS}
    # The picks are loaded up front for the SAME reason, and it is not hypothetical
    windows_by_year = {year: load_windows(exts[year]) for year in YEARS}

    print("=" * 84)
    print("NC-12 road offset measured from the extracted dune start")
    print("=" * 84)
    print(f"  extractor      : {EXTRACTOR.name}")
    for year in YEARS:
        e = exts[year]
        print(f"  {year} -> {e.TOPO_PRODUCT}/{e.VERSION}  "
              f"(ALONGSHORE_FLIP={e.ALONGSHORE_FLIP}, "
              f"STRAIGHTEN={e.STRAIGHTEN})")
        print(f"       DEMs  : {e.LOAD_PATH}")
        print(f"       picks : {e.WINDOW_JSON.name}")
    print(f"  reference      : interior row 0 = dune crest + 1 cell")
    print(f"  output root    : {OUT_ROOT}")

    audit: dict[int, dict] = {}

    for year in YEARS:
        print(f"\n--- {year} " + "-" * 68)
        ext = exts[year]
        windows = windows_by_year[year]
        print(f"  {ext.TOPO_PRODUCT}/{ext.VERSION} | {ext.WINDOW_JSON.name} "
              f"({sum(1 for k in windows if k != '_meta')} domains)")
        mask_dir = road_mask_dir(road_line_for_year(year))
        if not mask_dir.is_dir():
            print(f"  [skip] no mask folder: {mask_dir}")
            continue

        records, all_profiles = [], []
        for domain in DOMAINS:
            record, profiles = measure_domain(ext, domain, year, windows)
            records.append(record)
            all_profiles.extend(profiles)

        with_road = [r for r in records if r["n_road_profiles"] > 0]
        if not with_road:
            print("  [skip] no domain carried road")
            continue

        # Measured everywhere, written only within ROAD_SPAN
        in_span = [r for r in with_road
                   if ROAD_SPAN[0] <= r["domain"] <= ROAD_SPAN[1]]
        excluded = [r for r in with_road if r not in in_span]
        for r in excluded:
            r["flags"] = ",".join(filter(None, [r["flags"], "EXCLUDED_FROM_SPAN"]))
        if not in_span:
            print(f"  [skip] no domain with road inside ROAD_SPAN {ROAD_SPAN}")
            continue

        road_ids = [r["domain"] for r in in_span]
        first_gis, last_gis = min(road_ids), max(road_ids)
        if excluded:
            print(f"  excluded from span   : "
                  f"{[r['domain'] for r in excluded]} "
                  f"(measured and kept in _domains.csv, not written to the "
                  f"model-facing files)")

        # Second reference frame: the existing same-year measurement
        existing = read_two_row_csv(Path(str(EXISTING_SETBACK_FMT).format(year=year)))
        for r in records:
            same = existing.get(r["domain"], np.nan)
            r["setback_legacy_m"] = same
            r["delta_vs_legacy_m"] = (
                round(same - r["setback_dunestart_m"], 1)
                if np.isfinite(same) and np.isfinite(r["setback_dunestart_m"])
                else np.nan)

        out_dir = road_setback_dir(year)

        # Seaward relocation of roadways that drown at initialisation

        # Runs on the floored value, the one CASCADE indexes with
        relocate_drowning(ext, in_span, ext.CELL_SIZE_M)

        # Model-facing files
        write_two_row_csv(
            out_dir / f"RoadSetback_{year}_dunestart.csv",
            {r["domain"]: r["setback_model_m"] for r in in_span},
        )
        elev = {r["domain"]: r["road_elev_mhw_median"] for r in in_span
                if np.isfinite(r["road_elev_mhw_median"])}
        write_two_row_csv(out_dir / f"RoadElevation_{year}_dunestart.csv", elev)

        # Diagnostics: every measured domain, in-span or not
        write_csv(out_dir / f"RoadOffset_{year}_domains.csv", records)
        write_csv(out_dir / f"RoadOffset_{year}_profiles.csv",
                  [p for p in all_profiles
                   if ROAD_SPAN[0] <= p["domain"] <= ROAD_SPAN[1]])

        audit[year] = summarize(year, records, in_span, first_gis, last_gis)
        report(audit[year], records)

    if CONTROL_UNSTRAIGHTENED:
        run_control(exts)

    if audit:
        write_audit(audit, exts)
        print(f"\n[audit] {OUT_ROOT / 'RoadOffset_dunestart_audit.md'}")


# Re-measure with STRAIGHTEN = False so the frame effect is separable
def run_control(exts: dict) -> None:
    if not CONTROL_WINDOW_JSON.is_file():
        print(f"\n[control] skipped: no unstraightened picks at "
              f"{CONTROL_WINDOW_JSON}")
        return

    print(f"\n--- UNSTRAIGHTENED CONTROL " + "-" * 52)
    print(f"  picks: {CONTROL_WINDOW_JSON.name}")

    windows = json.load(open(CONTROL_WINDOW_JSON))
    saved = {year: e.STRAIGHTEN for year, e in exts.items()}
    for e in exts.values():
        e.STRAIGHTEN = False
    try:
        for year in YEARS:
            ext = exts[year]
            mask_dir = road_mask_dir(road_line_for_year(year))
            if not mask_dir.is_dir():
                continue
            records = [measure_domain(ext, d, year, windows)[0] for d in DOMAINS]
            for r in records:
                r.pop("setback_legacy_m", None)
                r.pop("delta_vs_legacy_m", None)
            out = (road_setback_dir(year)
                   / f"RoadOffset_{year}_domains{CONTROL_SUFFIX}.csv")
            write_csv(out, records)
            got = [r for r in records if r["n_road_profiles"] > 0]
            sb = np.array([r["setback_dunestart_m"] for r in got], dtype=float)
            print(f"  {year}: {len(got)} domains, setback median "
                  f"{np.nanmedian(sb):.0f} m, range {np.nanmin(sb):.0f} to "
                  f"{np.nanmax(sb):.0f} m  ->  {out.name}")
    finally:
        for year, was in saved.items():
            exts[year].STRAIGHTEN = was


# A list of dicts as CSV
def write_csv(path: Path, rows: list[dict]) -> None:
    import csv
    if not rows:
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


# Per-period summary statistics of the setbacks
def summarize(year, records, with_road, first_gis, last_gis) -> dict:
    sb = np.array([r["setback_dunestart_m"] for r in with_road], dtype=float)
    delta = np.array([r["delta_vs_legacy_m"] for r in with_road], dtype=float)
    negative = [r for r in with_road if r["setback_dunestart_m"] < 0]

    stations_year = load_stations(year)
    stations_2004 = load_stations(2004)
    retreat = None
    corr = np.nan
    if stations_year is not None and stations_2004 is not None and year != 2004:
        # Stations grow LANDWARD from the offshore datum
        retreat = stations_2004 - stations_year
        # Does delta_vs_legacy behave like retreat? It would correlate positively if so
        by_domain = {d: retreat[d - 1] for d in DOMAINS}
        pairs = [(r["delta_vs_legacy_m"], by_domain[r["domain"]])
                 for r in with_road
                 if np.isfinite(r["delta_vs_legacy_m"])
                 and r["domain"] in by_domain]
        if len(pairs) > 2:
            a, b = np.asarray(pairs, dtype=float).T
            corr = float(np.corrcoef(a, b)[0, 1])

    return {
        "year": year,
        "first_gis": first_gis, "last_gis": last_gis,
        "n_with_road": len(with_road),
        "n_no_road": sum(1 for r in records if r["flags"].startswith("NO_ROAD")),
        "setback_min": float(np.nanmin(sb)), "setback_max": float(np.nanmax(sb)),
        "setback_median": float(np.nanmedian(sb)),
        "n_negative": len(negative),
        "negative_domains": [(r["domain"], r["setback_dunestart_m"]) for r in negative],
        "n_drowned": sum(1 for r in with_road if r["drowns_at_init"] == "yes"),
        "relocated": [(r["domain"], r["setback_dunestart_floored_m"],
                       r["setback_model_m"], r["relocated_seaward_m"],
                       "BEYOND_MEASURED" in r["flags"])
                      for r in with_road if r["relocated_seaward_m"]],
        "unfixable": [r["domain"] for r in with_road
                      if "DROWNS_NO_VIABLE_ROW" in r["flags"]],
        "delta_median": float(np.nanmedian(delta)) if np.isfinite(delta).any() else np.nan,
        "delta_min": float(np.nanmin(delta)) if np.isfinite(delta).any() else np.nan,
        "delta_max": float(np.nanmax(delta)) if np.isfinite(delta).any() else np.nan,
        "retreat_median": float(np.median(retreat)) if retreat is not None else np.nan,
        "delta_retreat_corr": corr,
        "flagged": [(r["domain"], r["flags"]) for r in records if r["flags"]],
    }


# The period summary to the console
def report(s: dict, records: list[dict]) -> None:
    print(f"  road span            : GIS {s['first_gis']}-{s['last_gis']} "
          f"({s['n_with_road']} domains with road, {s['n_no_road']} without)")
    print(f"  setback vs 2009 dune : median {s['setback_median']:.0f} m, "
          f"range {s['setback_min']:.0f} to {s['setback_max']:.0f} m")
    print(f"  delta vs LEGACY file : median {s['delta_median']:.0f} m, "
          f"range {s['delta_min']:.0f} to {s['delta_max']:.0f} m")
    if np.isfinite(s["retreat_median"]):
        print(f"  offset-derived retreat: median {s['retreat_median']:.0f} m")
        print(f"  corr(delta, retreat) : {s['delta_retreat_corr']:+.3f}  "
              f"<- would be near +1 if the two files measured the same feature")
    print(f"  NEGATIVE (floored)   : {s['n_negative']} domain(s)"
          + (f" -> {[d for d, _ in s['negative_domains']]}" if s["n_negative"] else ""))
    print(f"  DROWNS at init       : {s['n_drowned']} domain(s)"
          + ("" if RELOCATE_DROWNING else "   [RELOCATE_DROWNING=False]"))
    if s["relocated"]:
        print(f"  moved seaward to the nearest viable row:")
        for d, was, now, moved, beyond in s["relocated"]:
            print(f"      GIS {d:>2}: {was:>5.0f} -> {now:>5.0f} m "
                  f"({moved:>4.0f} m seaward)"
                  + ("   BEYOND the measured per-profile range" if beyond
                     else "   within the measured range"))
    if s["unfixable"]:
        print(f"  [warn] no viable row seaward, still drowning: {s['unfixable']}")


# The one write-up for the whole tree, covering every year in `audit`
def write_audit(audit: dict, exts: dict) -> None:
    OUT_ROOT.mkdir(parents=True, exist_ok=True)
    path = OUT_ROOT / "RoadOffset_dunestart_audit.md"

    # The no-data explanation prints only when something drowns on this extraction
    total_drowned = sum(s["n_drowned"] for s in audit.values())

    # Any vintage the run did not produce is named, not silently left out
    any_ext = exts[sorted(audit)[0]]
    missing = [y for y in YEARS if y not in audit]
    runs = ", ".join(f"{y} on {exts[y].TOPO_PRODUCT}/{exts[y].VERSION}"
                     for y in sorted(audit))

    lines = [
        "# NC-12 road offset measured from the extracted dune start",
        "",
        f"Generated {datetime.now().isoformat(timespec='seconds')} by "
        "`HAT_road_offset_from_dune_start.py`.",
        "",
        f"Covers {runs}. Each vintage is measured against interior row 0 of "
        "**its own period's extraction** -- the two are different islands, and "
        "65 of 90 domains differ in interior shape between them.",
        "",
        *([f"> **INCOMPLETE RUN.** {', '.join(str(y) for y in missing)} "
           "produced nothing this run, so no section below describes it. The "
           "files on disk for that year are from an earlier run and are NOT "
           "documented here. Re-run to restore the section."]
          if missing else []),
        *([""] if missing else []),
        "| | |",
        "|---|---|",
        f"| Extractor | `{EXTRACTOR.name}` "
        f"(ALONGSHORE_FLIP={any_ext.ALONGSHORE_FLIP}, "
        f"STRAIGHTEN={any_ext.STRAIGHTEN}) |",
        *[f"| {y} topography | `{exts[y].TOPO_PRODUCT}/{exts[y].VERSION}` |"
          for y in sorted(audit)],
        *[f"| {y} DEMs | `{exts[y].LOAD_PATH}` |" for y in sorted(audit)],
        *[f"| {y} picked windows | `{exts[y].WINDOW_JSON.name}` |"
          for y in sorted(audit)],
        f"| Reference | interior row 0 = picked dune crest + 1 cell |",
        f"| Road reference | seaward-most road cell per profile |",
        f"| Alongshore collapse | median over profiles with road |",
        f"| Elevation | median of non-water road cells, m MHW |",
        f"| Assumed road width | {ASSUMED_ROAD_WIDTH_M:.0f} m |",
        "",
        "`road_setback` is metres landward of interior row 0 "
        "(`roadway_manager.py:99`). The road mask is put through the same "
        "orient / alongshore-flip / shear / trim chain as the topography, so "
        "the road and the dune line are compared on identical ground.",
        "",
        "## Negative setbacks are floored, not measured away",
        "",
        "A negative setback means the road falls SEAWARD of interior row 0. It "
        "cannot go into the model-facing CSV: `int(-50/10) = -5` and "
        "`xyz_interior_grid[-5:-3, :]` indexes from the landward end, so the "
        "road would be bulldozed into the bay with no error. `0.0` is not a "
        "\"no road\" sentinel either -- `build_roadway_management_on` uses the "
        "road span, not the value. The model-facing file therefore carries "
        "`max(setback, 0)`; the true signed value is in "
        "`RoadOffset_<year>_domains.csv` under `setback_dunestart_m`.",
        "",
        "## Roadways that drown at initialisation are moved seaward",
        "",
        "A roadway whose flanking rows are >20% water at t=0 is width-drowned "
        "by `roadway_manager.bulldoze` on the first call. `RoadwayManager` sets "
        "`_drown_break`, `cascade.py` never calls `update()` again, and that "
        "domain spends the whole hindcast as an **unmanaged barrier wearing a "
        "road label** -- no overwash removal, no dune rebuilding, no "
        "relocation. To keep the domain managed, the model-facing setback is "
        "moved to the nearest row seaward that passes bulldoze's own test.",
        "",
        "The test transcribed here is bulldoze's: the rows checked are the "
        "NEIGHBOURS of the bulldozed band (`road_start - 1`, `road_end + 1`), "
        "never the band itself, and every cell counts. There is **no cap** on "
        "the distance moved -- the nearest viable row is taken however far it "
        "is. `setback_dunestart_m` is never altered; only "
        "`RoadSetback_<year>_dunestart.csv` carries the moved value, and every "
        "moved domain is flagged `MOVED_SEAWARD`.",
        "",
        "### The assumption this rests on",
        "",
    ]

    if total_drowned == 0:
        lines += [
            f"**Nothing drowns at initialisation** on "
            f"{', '.join(exts[y].TOPO_PRODUCT + '/' + exts[y].VERSION for y in sorted(audit))}"
            " -- 0 "
            "domains in either year -- so no setback was moved and the rest of "
            "this section does not apply to this run.",
            "",
            "That is the point of the gap-filled DEM. On `2009_v4` three "
            "domains per year drowned, and the audit for that version showed "
            "they drowned on **LiDAR coverage gaps, not measured water**: "
            "across the six flanking rows failing at GIS 78/79/80 there were "
            "106 wet cells, of which 105 had never been surveyed and 1 was "
            "genuinely wet. The extractor writes no-data back as the water "
            "sentinel because Barrier3D has no representation for \"unknown\", "
            "so unsurveyed ground read as ocean and drowned the roadway. "
            "Filling those holes from the 2014 NOAA Post-Sandy DEM removes the "
            "cause, and the drown count goes to zero.",
            "",
            "Those figures are quoted from the `2009_v4` audit and were "
            "measured there, not recomputed here -- see "
            "`dunestart_offset_ARCHIVE_2009_v4/`, deleted 2026-08-26 - see "
            "1-barrier3d-domains/archive_purge_20260826.csv.",
            "",
        ]
    else:
        lines += [
            f"On the topography measured here these domains **do not drown on "
            "measured water -- they drown on LiDAR coverage gaps.** Across the "
            "six flanking rows that fail at GIS 78/79/80 there are 106 wet "
            "cells, of which **105 were never surveyed and 1 is genuinely "
            "measured wet.** The extractor writes no-data back as the water "
            "sentinel because Barrier3D has no representation for \"unknown\", "
            "and the drown test counts every cell -- deliberately, so that what "
            "is reported is what the run does.",
            "",
            "NOTE: that 106/105/1 split was measured on `2009_v4` and is not "
            "recomputed per run. Confirm it still holds before quoting it on a "
            "different extraction.",
            "",
            "So this is **not** \"the road was in water, we moved it out\". It "
            "is \"the 2009 DEM has no data there, CASCADE reads no-data as "
            "water, so the road is moved onto surveyed ground to keep the "
            "domain managed\". Any managed-vs-unmanaged result at GIS 78-80 "
            "needs that sentence attached to it.",
            "",
        ]

    lines += [
        "A move that lands inside the domain's own per-profile spread is a "
        "re-pick of the alongshore statistic. A move beyond it is a position no "
        "profile showed -- a different claim, flagged `BEYOND_MEASURED`.",
        "",
        "## `delta_vs_legacy_m` is NOT the dune-line retreat",
        "",
        "It is tempting to read the difference against the legacy "
        "`RoadSetback_<year>.csv` as the retreat between `<year>` and 2009. It "
        "is not, and the reported `corr(delta, retreat)` shows it: a like-for-"
        "like pair of measurements taken in two different years would correlate "
        "near +1 with the dune-line retreat, and the measured correlation is "
        "indistinguishable from zero -- the two files are not tracking the same "
        "feature at all. (It read \"strongly negative\" until 2026-09-22; the "
        "value has been about -0.03, which is no correlation, not a negative "
        "one. A correlation is immune to the zeroing error corrected the same "
        "day, so this conclusion did not depend on it.)",
        "",
        "Three reasons the two files are not commensurable, from "
        "`HAT_setback_from_lines.py` (retired 2026-08-17, git blob "
        "`c37ec03e`):",
        "",
        "1. **Different dune feature.** The legacy file computes "
        "`setback_i = road_cell_i - dune_cell_i` against a digitized same-year "
        "*dune-line geojson*. This script measures against the *DEM dune crest* "
        "found inside the picked search window. Those are different features, "
        "not the same feature at two dates.",
        "2. **Different frame.** The legacy measurement is taken with both lines "
        "\"raw, ocean-first\" (`HAT_setback_from_lines.py:254`) -- i.e. "
        "unstraightened, so it still carries the obliquity smear this script "
        "removes.",
        "3. **The legacy file is already floored.** It prints "
        "`FLOORED(x) raw setback was negative`, so at the domains where the two "
        "disagree most, the legacy value is itself a clamp, not a measurement.",
        "",
        "So `delta_vs_legacy_m` is a *migration diagnostic* -- how much each "
        "domain's forcing changes if you adopt this method -- and nothing more. "
        "Closing the same-year comparison properly needs a same-year dune crest "
        "measured the same way, which needs a same-year DEM.",
        "",
        "## The road reference is a buffered mask edge, not the road",
        "",
        "`road_setback` is measured to the seaward-most cell of the "
        "**rasterized mask**, which `HAT_rasterize_road_to_domains.py` burns with "
        "`ROAD_BUFFER_M = 6.0` and `ALL_TOUCHED = True` -- roughly a 24 m mask "
        "for an ~8 m road. So the reference is not NC-12's centerline, and not "
        "its pavement edge either. `HAT_road_buffer_bias.py` measures the "
        "difference against the source geojson, in the original raster frame "
        "(orient / flip / shear / trim all cancel out of it, since interior row 0 "
        "is common to both):",
        "",
        "| | 1984 | 2004 |",
        "|---|---|---|",
        "| median bias | -6.8 m | -6.8 m |",
        "| p10 / p90 | -10.8 / -2.7 m | -10.8 / -2.8 m |",
        "",
        "Negative means the mask edge sits **seaward** of the centerline, so "
        "every setback here runs about 7 m smaller than a centerline measurement "
        "would.",
        "",
        "### That is close to right, and correcting it would be worse",
        "",
        f"`bulldoze` lays a {ASSUMED_ROAD_WIDTH_M:.0f} m block LANDWARD from "
        "`road_start`, so a block centred on the real road wants "
        f"`road_start = centerline - {ASSUMED_ROAD_WIDTH_M / 2:.0f} m`. The "
        "buffer supplies `centerline - 6.9 m`, leaving a residual misplacement "
        "of **3.1 m -- about a third of a cell**, and an order of magnitude "
        "below the 20-40 m p90 placement error the alongshore collapse to one "
        "scalar already causes (`HAT_method_comparison_figures.py`).",
        "",
        "It is deliberately left alone. Per-profile setbacks are integer cell "
        "differences times the cell size, so domain medians land almost exactly "
        "on cell boundaries, and `bulldoze` truncates with `int()`. Applying the "
        "3.1 m adjustment moves `road_start` a **full cell (10 m) seaward in 83% "
        "of 1984 domains and 90% of 2004 domains** -- trading a 3 m error for a "
        "7 m error in the opposite direction. The buffer is accidentally doing "
        "very nearly the right job; the quantization is what would bite.",
        "",
        "Two independent paths agree on the number: the raster-frame column "
        "comparison above, and the `road_x`/`road_y` columns in "
        "`RoadOffset_<year>_profiles.csv`, which are produced by inverting the "
        "full index chain back to map coordinates and whose shapely distance to "
        "the same geojson is 6.62-6.68 m median (max 12.6 m, no point beyond "
        "60 m). That agreement also validates the index inversion -- an "
        "alongshore-flip error would mirror points within each 500 m domain and "
        "show up as distances in the hundreds of metres.",
        "",
        "D8 is the one extreme outlier (-353 m in 1984). That is the Buxton bend, "
        "where NC-12 runs parallel to the raster rows so a \"seaward-most cell\" "
        "has no physical reading -- the same reason D8 is already "
        "`EXCLUDED_FROM_SPAN`. Its appearance here is a consistency check "
        "passing, not a new problem.",
        "",
    ]
    for year, s in audit.items():
        lines += [
            f"## {year}",
            "",
            f"- Road span: GIS {s['first_gis']}-{s['last_gis']} "
            f"({s['n_with_road']} with road, {s['n_no_road']} without)",
            f"- Setback vs 2009 dune: median {s['setback_median']:.0f} m, "
            f"range {s['setback_min']:.0f} to {s['setback_max']:.0f} m",
            f"- Delta vs the LEGACY `RoadSetback_{year}.csv`: median "
            f"{s['delta_median']:.0f} m, range {s['delta_min']:.0f} to "
            f"{s['delta_max']:.0f} m",
        ]
        if np.isfinite(s["retreat_median"]):
            lines += [
                f"- Offset-derived dune-line retreat {year}->2004: median "
                f"{s['retreat_median']:.0f} m",
                f"- **corr(delta, retreat) = {s['delta_retreat_corr']:+.3f}** "
                f"-- see the caveat below.",
            ]
        lines.append(f"- NEGATIVE, floored to 0: {s['n_negative']} domain(s)")
        if s["negative_domains"]:
            lines += ["", "| GIS | true setback (m) |", "|---:|---:|"]
            lines += [f"| {d} | {v:.0f} |" for d, v in s["negative_domains"]]
        lines.append(f"- DROWNS at initialisation: {s['n_drowned']} domain(s)"
                     + ("" if RELOCATE_DROWNING
                        else "  (RELOCATE_DROWNING = False, none moved)"))
        if s["relocated"]:
            lines += [
                "",
                "| GIS | measured (m) | written (m) | moved seaward (m) | "
                "inside the measured per-profile range? |",
                "|---:|---:|---:|---:|---|",
            ]
            lines += [
                f"| {d} | {was:.0f} | {now:.0f} | {moved:.0f} | "
                + ("**no -- beyond it**" if beyond else "yes") + " |"
                for d, was, now, moved, beyond in s["relocated"]
            ]
        if s["unfixable"]:
            lines.append(f"- No viable row seaward, still drowning: "
                         f"{s['unfixable']}")
        if s["flagged"]:
            lines += ["", "### Flags", "", "| GIS | flags |", "|---:|---|"]
            lines += [f"| {d} | {f} |" for d, f in s["flagged"]]
        lines.append("")

    path.write_text("\n".join(lines), encoding="utf-8")


if __name__ == "__main__":
    main()
