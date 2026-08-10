"""
HAT_road_setback_extract.py
===============================================================================
Compute CASCADE-ready NC-12 setbacks from road masks, using HAT_dune_topo_
extractor.py's own functions so the setback origin cannot drift from the one
CASCADE actually indexes.

WHY IT IMPORTS RATHER THAN REIMPLEMENTS
---------------------------------------
InteriorDomain[0] is defined by a chain inside HAT_dune_topo_extractor.py:

    orient_ocean_right(arr, OCEAN_LOC)      # ocean to the last column
      -> [:, ::-1]                          # ocean-first
      -> [:, c0:c1+1]                       # water-column trim
      -> find_dunes(...)  -> dune_loc
      -> build_interior(...) -> start_island = max(dune_loc) + 1
      -> remove_water_rows(...)             # can drop LEADING rows, moving row 0

Any second script that reconstructs that chain can get it subtly wrong, and
that is precisely the bug this whole exercise is about: setbacks measured
against a 1984 aerial dune line, handed to a grid whose row 0 is a 2009 LiDAR
dune crest, put NC-12 in Pamlico Sound and raised

    IndexError: index 15 is out of bounds for axis 0 with size 15
    roadway_manager.bulldoze -> xyz_interior_grid[road_end + 1, :]

So this script does not reconstruct anything. It imports the module and calls
load_profiles / find_dunes / build_interior with the same picks JSON. The
start_island it gets is the start_island your extraction run got, because it is
the same function on the same input.

Re-running find_dunes/build_interior is cheap and deterministic. Nothing is
re-extracted, nothing is overwritten.

THE TRIM SUBTLETY
-----------------
remove_water_rows() keeps rows[min(keep):max(keep)+1]. min(keep) is not always
0. If the cells immediately landward of the most-landward dune are water in
every profile (a washover throat, a pond), those rows are trimmed and row 0 of
the saved array is NOT start_island -- the origin shifts landward by min(keep),
silently. This script recomputes that shift with the same rule and measures
from start_island + shift. Watch the row0_shift column.

INPUTS
------
  <extractor's LOAD_PATH>/domain_<N>.npy          the DEM domains
  <extractor's WINDOW_JSON>                       your picked dune windows
  ROAD_MASK_DIR/domain_<N>_road_<YEAR>.npy        from HAT_rasterize_road_to_domains.py

OUTPUTS
-------
  OUT_DIR/RoadSetback_<YEAR>.csv            2-row CASCADE file (IDs, setbacks)
  OUT_DIR/RoadSetback_<YEAR>_domains.csv    per-domain stats + QC flags
  OUT_DIR/RoadSetback_<YEAR>_profiles.csv   every profile measurement
  OUT_DIR/HAT_road_setback_qc_<YEAR>.png    alongshore QC figure

REQUIREMENTS
------------
  HAT_dune_topo_extractor.py importable (same folder, or set EXTRACTOR_DIR)
  numpy, matplotlib
===============================================================================
"""

from __future__ import annotations

import csv
import os
import sys
import warnings
from pathlib import Path

import numpy as np

# --- where HAT_dune_topo_extractor.py lives -----------------------------
# This script imports the extractor rather than reimplementing it, so it has to
# find it. Two options:
#
#   None            -> look in THIS script's own folder. Use this if you drop
#                      HAT_road_setback_extract.py next to the extractor.
#   r"C:\...\dir"   -> the folder holding HAT_dune_topo_extractor.py. Use this
#                      to keep the road scripts together in e.g.
#                      scripts\input_prep\road_offset\ while the
#                      extractor stays where it is.
EXTRACTOR_DIR = None

_ed = (Path(EXTRACTOR_DIR) if EXTRACTOR_DIR
       else Path(__file__).resolve().parent)
if str(_ed) not in sys.path:
    sys.path.insert(0, str(_ed))

try:
    import HAT_dune_topo_extractor as ext  # noqa: E402
except ModuleNotFoundError:
    raise SystemExit(
        f"\nCannot import HAT_dune_topo_extractor from:\n    {_ed}\n\n"
        f"This script reuses the extractor's own load_profiles / find_dunes /\n"
        f"build_interior so the setback origin cannot drift from the one\n"
        f"CASCADE indexes. It needs the real module, not a copy.\n\n"
        f"Either put this file next to HAT_dune_topo_extractor.py, or set\n"
        f"EXTRACTOR_DIR at the top of this file to the folder holding it.\n"
    )

import matplotlib                       # noqa: E402
matplotlib.use("Agg")
import matplotlib.pyplot as plt         # noqa: E402
from matplotlib.transforms import blended_transform_factory  # noqa: E402


# =============================================================================
# CONFIG
# =============================================================================

ROAD_YEAR = 1984

# Masks burned onto the RAW domain_<N>.npy grid by HAT_rasterize_road_to_domains.py
ROAD_MASK_DIR = ext.INIT_ROOT / "roads" / "raster" / str(ROAD_YEAR) / "masks"
ROAD_FILE_FMT = "{stem}_road_{year}.npy"

# Road outputs live INSIDE the extraction run that defined their origin, one
# subfolder per road vintage:
#
#   dune_topo\2009_v1\
#       topography\                     <- the grid the setback is measured on
#       dunes\
#       roads\
#           1984\                       <- this run
#               RUN_MANIFEST.txt
#               RoadSetback_1984.csv
#               RoadElevation_1984.csv
#               RoadSetback_1984_domains.csv
#               RoadSetback_1984_profiles.csv
#               HAT_road_setback_qc_1984.png
#               figures\domain_0NN_setback.png
#           1997\                       <- for the relocation displacements
#           2004\
#
# The run folder is not filing preference. The setback is measured from
# start_island, which the dune search windows in THIS run define -- re-pick a
# window in 2009_v2 and every setback in that domain moves. Keeping them inside
# the run means a setback can never be paired with topography it was not
# measured against, which is the frame mismatch this script exists to prevent.
OUT_DIR = ext.RUN_DIR / "roads" / str(ROAD_YEAR)

ROAD_AGG = "median"       # across the 50 profiles: "median" | "mean"

# WHICH ORIGIN EACH PROFILE'S SETBACK IS MEASURED FROM
# ----------------------------------------------------
# "const"        setback_i = road_loc_i - start_island, where
#                start_island = max(dune_loc) + 1 -- ONE constant for the whole
#                domain, exactly what build_interior() uses.
#
# "per_profile"  setback_i = road_loc_i - (dune_loc_i + 1) -- each profile
#                against its OWN dune.
#
# Why this matters: the domain grids are NORTH-UP (rot = 0.00 for all 131 --
# see HAT_rasterize_road_to_domains.py), so cross-shore is due east-west while
# the island trends NNW. The shoreline slides ~80 m in easting across a 500 m
# domain, so dune_loc runs DIAGONALLY -- e.g. GIS 10 goes from cell 4 at
# profile 0 to cell 12 at profile 49.
#
# Measured against a single constant origin, that diagonal produces a linear
# ramp of setbacks (GIS 10: -100 m to -20 m) whose median describes no profile
# and can be negative even though the road is where it should be. Every NEG
# domain is a domain where the shoreline is oblique AND the road hugs the dune.
#
# "per_profile" removes the obliquity from the measurement: it asks how far
# behind ITS dune each bit of asphalt sits, which is the physical quantity
# CASCADE wants, since a Barrier3D domain is an idealized shore-normal profile.
#
# It does NOT fix the obliquity itself. build_interior() still starts every
# profile at the same cell, so a diagonal wedge of real barrier is cut from the
# interior in exactly these domains -- that bias is in your interior widths and
# dune heights too, not just the road. The real fix is rotating the clip boxes
# to shore-normal, which is a re-extraction.
#
# Both are computed and reported. This chooses which goes in RoadSetback_*.csv.
SETBACK_REF = "per_profile"

# FLOORING SMALL NEGATIVES
# ------------------------
# Even measured per-profile, GIS 9-16 come out around -20 m: NC-12 in the
# inter-village stretch runs immediately behind the foredune, and with a ~24 m
# mask (ROAD_BUFFER_M + all_touched) on 10 m resampled cells, the road and the
# dune ARE the same two or three cells. Which is marginally seaward is below
# what the DEM can resolve.
#
# The physically true statement there is "the road is at the dune line", i.e.
# setback = 0, and CASCADE handles that fine: int(0/10) = 0, road_end = 2,
# bulldoze reads row 3.
#
# So a negative shallower than the tolerance is floored to zero and counted.
# Anything DEEPER is not resolution -- it is a real disagreement between the
# dune pick and the road, and the script still refuses to write it rather than
# handing CASCADE a value it would silently wrap to the back of the barrier.
#
# The tolerance is the mask half-width plus a cell: 6 m buffer -> ~12 m of mask
# either side of the centreline, plus 10 m of cell. Not a tuning knob.
SETBACK_FLOOR_M = 0.0
SETBACK_FLOOR_TOL_M = 30.0

FIRST_ROAD_DOMAIN = 9     # GIS domains carrying NC-12; must match the runner
LAST_ROAD_DOMAIN = 90

# MUST match ROAD_WIDTH in HAT_hindcast_1984_2024_newplot.py. Used only to
# predict whether bulldoze() runs off the back of the interior.
ROAD_WIDTH_M = 20.0

# QC thresholds
MIN_ROAD_PROFILES = 10    # flag domains where the road was found on fewer
MAX_ROAD_STD_M = 40.0     # flag domains whose profiles disagree by more
ROAD_ELEV_MIN_M = 0.6     # m MHW; flag road cells sitting on ground below this

# Per-domain diagnostic figures: exactly what measure_road() sees.
# None = every road domain (~82 figures). A list = just those.
SAVE_DOMAIN_FIGS = True
QC_DOMAINS = [9, 10, 11, 12, 13, 16, 22, 35, 74, 84, 85, 87]

# Style (HAT_hindcast_1984_2024.py palette)
C_MODEL, C_TOWN, C_PIER, C_GROIN, C_INK = ("#FF8C00", "#90AFC5", "#1565C0",
                                           "#B71C1C", "#1a1a2e")


# =============================================================================
# HELPERS
# =============================================================================

def leading_water_rows(topo_m: np.ndarray) -> int:
    """
    How many leading rows remove_water_rows() would drop. Same rule as the
    extractor's, so the number matches what was actually saved.
    """
    keep = [r for r in range(topo_m.shape[0])
            if not np.all(topo_m[r, :] == ext.SENTINEL_WATER_M)]
    return min(keep) if keep else 0


def trimmed_rows(topo_m: np.ndarray) -> int:
    """Interior rows CASCADE will see, honouring TRIM_INTERIOR_ROWS."""
    if not ext.TRIM_INTERIOR_ROWS:
        return int(topo_m.shape[0])
    keep = [r for r in range(topo_m.shape[0])
            if not np.all(topo_m[r, :] == ext.SENTINEL_WATER_M)]
    return (max(keep) - min(keep) + 1) if keep else 0


def load_road_mask(elev_path: Path, road_path: Path, c0: int,
                   n_cross_trimmed: int):
    """
    Put the road mask in the SAME frame as dom["z"], using the extractor's own
    orientation function. Any deviation here is the whole bug, so there is no
    clever shortcut.
    """
    if not road_path.exists():
        return None, "NO_MASK"

    raw_shape = np.load(elev_path, mmap_mode="r").shape
    r = np.load(road_path)
    if r.shape != raw_shape:
        raise ValueError(
            f"{elev_path.name}: road mask shape {r.shape} != elevation shape "
            f"{raw_shape}. The mask must be burned onto the RAW domain grid. "
            f"Re-run HAT_rasterize_road_to_domains.py."
        )

    r = ext.orient_ocean_right(r.astype(float), ext.OCEAN_LOC)  # incl. ALONGSHORE_FLIP
    r = r[:, ::-1]                          # ocean-first, matching `raw`
    r = r[:, c0:c0 + n_cross_trimmed]       # same water-column trim as `z`
    mask = r > 0
    if not mask.any():
        return mask, "EMPTY_AFTER_TRIM"
    return mask, ""


def measure_road(road_mask, setback_origin, n_along, dune_loc=None,
                 row0_shift=0):
    """
    Per-profile setback, computed two ways.

    sb_const : road_loc_i - setback_origin, where setback_origin is the single
               start_island (+ trim shift) that build_interior() uses for every
               profile. This is the frame CASCADE's InteriorDomain is in.

    sb_pp    : road_loc_i - (dune_loc_i + 1 + row0_shift), each profile against
               its own dune. On a north-up grid with an oblique shoreline,
               dune_loc runs diagonally and only this one is meaningful.

    road_loc is the first road cell at or landward of the constant origin; if
    the whole road is seaward of it (which the obliquity causes), the
    landwardmost road cell is used instead, so the measurement still refers to
    the road rather than silently vanishing.
    """
    sb_const = np.full(n_along, np.nan)
    sb_pp = np.full(n_along, np.nan)
    loc = np.full(n_along, -1, dtype=int)
    if road_mask is None or setback_origin is None:
        return sb_const, sb_pp, loc

    for i in range(min(n_along, road_mask.shape[0])):
        hits = np.where(road_mask[i])[0]
        if hits.size == 0:
            continue
        landward = hits[hits >= setback_origin]
        idx = int(landward[0]) if landward.size else int(hits[-1])
        loc[i] = idx
        sb_const[i] = (idx - setback_origin) * ext.CELL_SIZE_M
        if dune_loc is not None and i < len(dune_loc) and dune_loc[i] >= 0:
            sb_pp[i] = (idx - (int(dune_loc[i]) + 1 + row0_shift)) \
                * ext.CELL_SIZE_M
    return sb_const, sb_pp, loc


# =============================================================================
# PER-DOMAIN
# =============================================================================

def process(stem: str, elev_path: Path, windows: dict):
    dom = ext.load_profiles(elev_path)
    z, start_beach, c0 = dom["z"], dom["start_beach"], dom["c0"]
    n_along = z.shape[0]

    w = windows.get(stem)
    if w is None:
        i0, i1 = ext.default_window(z, start_beach)
        wsrc = "default"
        print(f"[warn] {stem}: no picked window, using default [{i0}, {i1}]")
    else:
        if w.get("n_cross_trimmed") not in (None, z.shape[1]):
            print(f"[warn] {stem}: trimmed width changed "
                  f"({w.get('n_cross_trimmed')} -> {z.shape[1]}); the saved "
                  f"window may no longer line up. Re-pick.")
        i0, i1 = int(w["i0"]), int(w["i1"])
        wsrc = "picked"

    # --- identical to the extraction run: same functions, same picks ---
    dune_elev, dune_loc = ext.find_dunes(z, start_beach, i0, i1)
    if not np.isfinite(dune_elev).any():
        print(f"[skip] {stem}: no dune found in window [{i0}, {i1}]")
        return None

    topo_m, start_island = ext.build_interior(z, dune_loc)
    if start_island is None:
        raise RuntimeError(
            f"{stem}: build_interior returned start_island=None. This script "
            f"needs USE_CONST_INTERIOR = True -- with False each profile starts "
            f"behind its own dune and there is no single row 0 for a setback."
        )

    shift = leading_water_rows(topo_m) if ext.TRIM_INTERIOR_ROWS else 0
    setback_origin = start_island + shift     # the row bulldoze() indexes from
    n_rows = trimmed_rows(topo_m)

    road_path = Path(ROAD_MASK_DIR) / ROAD_FILE_FMT.format(stem=stem,
                                                           year=ROAD_YEAR)
    mask, mask_note = load_road_mask(elev_path, road_path, c0, z.shape[1])
    sb_const, sb_pp, loc = measure_road(mask, setback_origin, n_along,
                                        dune_loc=dune_loc, row0_shift=shift)
    sb = sb_pp if SETBACK_REF == "per_profile" else sb_const
    good = sb[np.isfinite(sb)]
    good_const = sb_const[np.isfinite(sb_const)]
    good_pp = sb_pp[np.isfinite(sb_pp)]

    # dune_loc's alongshore span is the obliquity, in cells. A north-up grid on
    # an NNW-trending island gives a diagonal dune line; that is what turns a
    # constant origin into a ramp of setbacks.
    dl = np.asarray([d for d in dune_loc if d >= 0], dtype=float)
    dune_span = int(dl.max() - dl.min()) if dl.size else -1

    # --- road elevation, per domain -------------------------------------
    # Two populations, because they answer different questions:
    #   *_all      every cell in the buffered mask. Includes shoulder and
    #              adjacent ground, so it is the roadbed CORRIDOR, not the
    #              crown. Better for registration QC.
    #   *_atroad   only the cell each profile's setback actually lands on
    #              (road_loc). This is the ground CASCADE bulldozes to road_ele,
    #              so it is the one to use if you ever set road_ele per domain.
    ev = {}
    if mask is not None and mask.any():
        cells = z[mask]
        cells = cells[cells > ext.SENTINEL_WATER_M + 1e-9]
        if cells.size:
            ev.update(
                road_elev_mhw=float(np.median(cells)),
                road_elev_mean_mhw=float(np.mean(cells)),
                road_elev_std_m=float(np.std(cells)),
                road_elev_min_mhw=float(np.min(cells)),
                road_elev_max_mhw=float(np.max(cells)),
                n_road_cells=int(cells.size),
            )
    at = [z[i, loc[i]] for i in range(len(loc))
          if loc[i] >= 0 and z[i, loc[i]] > ext.SENTINEL_WATER_M + 1e-9]
    if at:
        at = np.asarray(at, dtype=float)
        ev.update(
            road_elev_atroad_mhw=float(np.median(at)),
            road_elev_atroad_mean_mhw=float(np.mean(at)),
            road_elev_atroad_std_m=float(np.std(at)),
            n_atroad_cells=int(at.size),
        )

    ev.setdefault("road_elev_mhw", np.nan)
    for k in ("road_elev_mean_mhw", "road_elev_std_m", "road_elev_min_mhw",
              "road_elev_max_mhw", "road_elev_atroad_mhw",
              "road_elev_atroad_mean_mhw", "road_elev_atroad_std_m"):
        ev.setdefault(k, np.nan)
    for k in ("n_road_cells", "n_atroad_cells"):
        ev.setdefault(k, 0)

    # NAVD88 for anything that has to leave the extractor's frame. The arrays
    # are MHW-relative (the extractor subtracts MHW_M before anything else), so
    # the MHW columns are the ones CASCADE's grid speaks -- see the note in
    # write_elevation_csv().
    ev["road_elev_navd"] = ev["road_elev_mhw"] + ext.MHW_M
    ev["road_elev_atroad_navd"] = ev["road_elev_atroad_mhw"] + ext.MHW_M

    if good.size:
        val = float(np.median(good) if ROAD_AGG == "median" else np.mean(good))
        # bulldoze(): road_start = int(setback/10)
        #             road_end   = road_start + int(road_width/10)
        #             reads xyz_interior_grid[road_end + 1, :]
        road_end = (int(val / ext.CELL_SIZE_M)
                    + int(ROAD_WIDTH_M / ext.CELL_SIZE_M))
        need = road_end + 2                   # check row + margin
        fits = need <= n_rows
    else:
        val, road_end, need, fits = np.nan, -1, -1, False

    if shift:
        print(f"       [note] {stem}: trim dropped {shift} leading water "
              f"row(s); origin {start_island} -> {setback_origin}")

    return dict(
        domain=ext.domain_number(stem), stem=stem,
        section=ext.section_for(stem),
        window=f"[{i0},{i1}]", window_source=wsrc,
        start_island=start_island, row0_shift=shift,
        setback_origin=setback_origin,
        interior_rows=n_rows,
        setback_m=val,
        setback_ref=SETBACK_REF,
        setback_const_m=(float(np.median(good_const)) if good_const.size
                         else np.nan),
        setback_perprofile_m=(float(np.median(good_pp)) if good_pp.size
                              else np.nan),
        const_std_m=(float(np.std(good_const)) if good_const.size else np.nan),
        perprofile_std_m=(float(np.std(good_pp)) if good_pp.size else np.nan),
        dune_span_cells=dune_span,
        obliquity_deg=(round(float(np.degrees(np.arctan2(
            dune_span * ext.CELL_SIZE_M, n_along * ext.CELL_SIZE_M))), 1)
            if dune_span > 0 else 0.0),
        n_road_hits=int(good.size), n_profiles=n_along,
        road_std_m=float(np.std(good)) if good.size else np.nan,
        road_min_m=float(np.min(good)) if good.size else np.nan,
        road_max_m=float(np.max(good)) if good.size else np.nan,
        n_negative=int((good < 0).sum()),
        # road_elev_mhw and the rest of the elevation stats come in via **ev
        **ev,
        road_end_row=road_end, rows_needed=need, road_fits=bool(fits),
        mask_note=mask_note,
        _sb=sb, _loc=loc, _sb_const=sb_const, _sb_pp=sb_pp,
        _z=z, _mask=mask, _dune_loc=dune_loc, _i0=i0, _i1=i1,
    )


# =============================================================================
# OUTPUT
# =============================================================================

def flags_for(r):
    f = []
    if r["mask_note"]:
        f.append(r["mask_note"])
    if not np.isfinite(r["setback_m"]):
        f.append("NO_ROAD")
        return f
    if r["n_road_hits"] < MIN_ROAD_PROFILES:
        f.append("THIN")
    if r["road_std_m"] > MAX_ROAD_STD_M:
        f.append("SCATTER")
    sb = r["setback_m"]
    if sb < -SETBACK_FLOOR_TOL_M:
        f.append(f"NEG({sb:.0f}m)")
    elif sb < SETBACK_FLOOR_M:
        f.append("FLOORED")
    if np.isfinite(r["road_elev_mhw"]) and r["road_elev_mhw"] < ROAD_ELEV_MIN_M:
        f.append("LOW_ELEV")
    if not r["road_fits"]:
        f.append("OFF_BACK")
    return f


def report(rows):
    by = {r["domain"]: r for r in rows if r["domain"] is not None}
    expected = list(range(FIRST_ROAD_DOMAIN, LAST_ROAD_DOMAIN + 1))

    print("\n" + "=" * 86)
    print(f"ROAD SETBACKS  |  {ROAD_YEAR}  |  run {ext.RUN_NAME}")
    print("=" * 86)
    print(f"  setback_ref = {SETBACK_REF}   "
          f"(const = one start_island for the domain; "
          f"per_profile = each profile's own dune)")
    print(f"  obliq = shoreline angle to the grid's east-west axis. The grids "
          f"are north-up, so")
    print(f"          any obliquity makes dune_loc diagonal and turns a "
          f"constant origin into a ramp.\n")
    print(f"{'GIS':>4} {'setback':>8} {'const':>7} {'per_prof':>9} "
          f"{'obliq':>6} {'origin':>7} {'hits':>5} {'elev':>6} "
          f"{'need':>5} {'rows':>5}  flag")

    missing, flagged, vals = [], [], []
    for d in expected:
        r = by.get(d)
        if r is None:
            missing.append(d)
            print(f"{d:>4} {'--':>8} {'--':>7} {'--':>6} {'--':>5} {'--':>5} "
                  f"{'--':>6} {'--':>5} {'--':>5}  NOT_EXTRACTED")
            continue
        f = flags_for(r)
        if f:
            flagged.append((d, ",".join(f)))
        if not np.isfinite(r["setback_m"]):
            missing.append(d)

        def _f(v, w=8, dp=0):
            return "--".rjust(w) if not np.isfinite(v) else f"{v:>{w}.{dp}f}"
        elt = "--" if not np.isfinite(r["road_elev_mhw"]) else f"{r['road_elev_mhw']:.2f}"
        print(f"{d:>4} {_f(r['setback_m']):>8} {_f(r['setback_const_m'],7)} "
              f"{_f(r['setback_perprofile_m'],9)} {r['obliquity_deg']:>6.1f} "
              f"{r['setback_origin']:>7} {r['n_road_hits']:>5} {elt:>6} "
              f"{r['rows_needed']:>5} {r['interior_rows']:>5}  {','.join(f)}")
        if np.isfinite(r["setback_m"]):
            vals.append(r["setback_m"])

    print()
    if flagged:
        print(f"  {len(flagged)} domain(s) flagged:\n")
        for d, f in flagged:
            print(f"    GIS {d:>2}: {f}")
        print()
        print("    OFF_BACK   setback + road + bulldoze's check row does not fit")
        print("               in the interior. THIS IS THE IndexError, predicted")
        print("               at t=0. It only gets worse as the barrier narrows.")
        print("    LOW_ELEV   road on ground below "
              f"{ROAD_ELEV_MIN_M} m MHW. NC-12 sits near the berm")
        print("               (1.70 m NAVD88 = 1.34 m MHW), so this is marsh,")
        print("               not roadbed -- the mask is misregistered.")
        print("    NEG(n)     road found seaward of row 0 -- misregistered mask,")
        print("               or a dune window picked too far landward.")
        print("    SCATTER    profiles disagree by > "
              f"{MAX_ROAD_STD_M:.0f} m -- road not shore-parallel here, or a")
        print("               spur/side street got digitized in.")
        print("    THIN       road on < "
              f"{MIN_ROAD_PROFILES} profiles -- mask gaps. Raise ROAD_BUFFER_M")
        print("               in HAT_rasterize_road_to_domains.py.")
    else:
        print("  No domains flagged. Every road fits inside its interior.")

    if vals:
        v = np.array(vals)
        print(f"\n  setback ({SETBACK_REF}): min {v.min():.0f} | "
              f"median {np.median(v):.0f} | max {v.max():.0f} m")

    diff = np.array([r["setback_const_m"] - r["setback_perprofile_m"]
                     for r in rows
                     if np.isfinite(r.get("setback_const_m", np.nan))
                     and np.isfinite(r.get("setback_perprofile_m", np.nan))])
    ob = np.array([r["obliquity_deg"] for r in rows])
    if diff.size:
        print(f"\n  const - per_profile: median {np.median(diff):+.0f} m, "
              f"max |{np.abs(diff).max():.0f}| m")
        print(f"  shoreline obliquity: median {np.median(ob):.1f} deg, "
              f"max {ob.max():.1f} deg")
        print(f"  The two agree where the shoreline is square to the grid and")
        print(f"  diverge where it is not. That divergence is the north-up")
        print(f"  clip, and it biases interior_rows and dune heights the same")
        print(f"  way -- not just the road.")

    shifts = [r["row0_shift"] for r in rows if r["row0_shift"]]
    if shifts:
        print(f"\n  {len(shifts)} domain(s) had leading water rows trimmed, "
              f"moving the setback origin.")
        print(f"  Max shift {max(shifts)} rows ({max(shifts) * 10} m). Measured "
              f"from the shifted origin, as CASCADE will.")

    return missing


def domain_setback_figure(r, path: Path):
    """
    Draw exactly what measure_road() sees, in the trimmed ocean-first frame the
    setback is computed in.

    Left  : the domain, ocean at the bottom. Per-profile dune pick (orange),
            the picked search window (orange band), start_island / the setback
            origin (black line), and the road mask (red).
    Right : alongshore-mean profile, same markers, with the berm drawn.
    Bottom: per-profile setback. This is what the median is taken over -- if it
            straddles zero, the origin and the road are interleaved and the
            median is meaningless.

    A NEG domain shows up here immediately: the black origin line sits LANDWARD
    of the red road.
    """
    z, mask = r["_z"], r["_mask"]
    dune_loc, loc, sb = r["_dune_loc"], r["_loc"], r["_sb"]
    n_along, n_cross = z.shape
    origin = r["setback_origin"]

    fig = plt.figure(figsize=(15, 9), facecolor="white")
    gs = fig.add_gridspec(2, 2, height_ratios=[2.4, 1], width_ratios=[1.4, 1],
                          hspace=0.22, wspace=0.06)
    ax_map = fig.add_subplot(gs[0, 0])
    ax_prof = fig.add_subplot(gs[0, 1], sharey=ax_map)
    ax_sb = fig.add_subplot(gs[1, :])

    zz = np.where(z <= ext.SENTINEL_WATER_M + 1e-9, np.nan, z)
    finite = zz[np.isfinite(zz)]
    vmax = float(np.percentile(finite, 99)) if finite.size else 3.0

    im = ax_map.imshow(np.ma.masked_invalid(zz.T), aspect="auto", origin="lower",
                       extent=[-0.5, n_along - 0.5, -0.5, n_cross - 0.5],
                       cmap="terrain", vmin=-1.0, vmax=max(vmax, 2.0))
    ax_map.axhspan(r["_i0"], r["_i1"], color=C_MODEL, alpha=0.18, zorder=1,
                   label=f"dune search window [{r['_i0']},{r['_i1']}]")
    dl = np.where(np.asarray(dune_loc) >= 0, dune_loc, np.nan)
    ax_map.plot(np.arange(len(dl)), dl, lw=0, marker="o", ms=3, color=C_MODEL,
                label="dune pick (per profile)")
    ax_map.axhline(origin, color="k", lw=2.0,
                   label=f"start_island + shift = {origin}  (setback origin)")
    ai, ci = np.where(mask) if mask is not None else (np.array([]), np.array([]))
    ax_map.plot(ai, ci, lw=0, marker="s", ms=2.5, color=C_GROIN,
                label=f"NC-12 {ROAD_YEAR} mask")
    rl = np.where(np.asarray(loc) >= 0, loc, np.nan)
    ax_map.plot(np.arange(len(rl)), rl, lw=1.2, color="w", alpha=0.9,
                label="road cell used")
    ax_map.set_xlabel("alongshore profile")
    ax_map.set_ylabel("cross-shore cell  (0 = ocean, landward up)")
    ax_map.legend(loc="upper left", fontsize=8, framealpha=0.92)

    y = np.arange(n_cross)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        prof = np.nanmean(zz, axis=0)
    ax_prof.plot(zz.T, y, color="0.85", lw=0.4)
    ax_prof.plot(prof, y, color=C_INK, lw=2.0, label="alongshore-mean")
    ax_prof.axhline(origin, color="k", lw=2.0)
    ax_prof.axhspan(r["_i0"], r["_i1"], color=C_MODEL, alpha=0.18)
    if np.isfinite(rl).any():
        ax_prof.axhspan(np.nanmin(rl), np.nanmax(rl), color=C_GROIN, alpha=0.22,
                        label="road span")
    ax_prof.axvline(ext.BERM_ELEV_NAVD_M - ext.MHW_M, color=C_GROIN, ls=":",
                    lw=1.4, label="berm")
    ax_prof.axvline(0, color=C_PIER, ls="--", lw=1.0, label="MHW")
    ax_prof.set_xlabel("elev (m MHW)")
    ax_prof.set_xlim(-1.5, max(vmax, 2.0) + 0.5)
    ax_prof.legend(loc="upper right", fontsize=8, framealpha=0.9)
    plt.setp(ax_prof.get_yticklabels(), visible=False)

    sbc, sbp = r["_sb_const"], r["_sb_pp"]
    gc, gp = np.isfinite(sbc), np.isfinite(sbp)
    ax_sb.plot(np.arange(len(sbc))[gc], np.asarray(sbc)[gc], lw=0, marker="o",
               ms=4, color="0.6", label="vs constant start_island")
    ax_sb.plot(np.arange(len(sbp))[gp], np.asarray(sbp)[gp], lw=0, marker="o",
               ms=4, color=C_GROIN, label="vs each profile's own dune")
    ax_sb.axhline(0, color="k", lw=1.4, label="setback = 0 (the dune)")
    if np.isfinite(r["setback_m"]):
        ax_sb.axhline(r["setback_m"], color=C_MODEL, lw=1.8,
                      label=f"{ROAD_AGG} = {r['setback_m']:.0f} m")
    ax_sb.set_xlabel("alongshore profile")
    ax_sb.set_ylabel("setback (m)")
    ax_sb.legend(loc="upper right", fontsize=8, framealpha=0.9)
    ax_sb.grid(alpha=0.25)

    fig.colorbar(im, ax=[ax_map, ax_prof], location="right", fraction=0.03,
                 pad=0.02, label="elev (m MHW)")

    flags = ",".join(flags_for(r)) or "clean"
    fig.suptitle(
        f"GIS {r['domain']}  ({r['section']})  —  NC-12 {ROAD_YEAR} setback\n"
        f"origin {origin} (start_island {r['start_island']} + shift "
        f"{r['row0_shift']})  |  setback {r['setback_m']:.0f} m "
        f"({SETBACK_REF})  |  needs {r['rows_needed']} of "
        f"{r['interior_rows']} rows  |  {flags}\n"
        f"dune_loc spans {r['dune_span_cells']} cells across the domain "
        f"= {r['obliquity_deg']:.1f} deg of shoreline obliquity  |  "
        f"const {r['setback_const_m']:.0f} m vs per-profile "
        f"{r['setback_perprofile_m']:.0f} m",
        fontsize=11)

    path.mkdir(parents=True, exist_ok=True)
    fig.savefig(path / f"domain_{r['domain']:03d}_setback.png", dpi=130,
                bbox_inches="tight")
    plt.close(fig)


def write_manifest(rows, path: Path):
    """
    Record what made this folder. The outputs no longer sit inside the
    extraction run, so the link back to it has to be written down.
    """
    from datetime import datetime

    cfg = {
        "run_time": datetime.now().isoformat(timespec="seconds"),
        "script": Path(__file__).name,
        "ROAD_YEAR": ROAD_YEAR,
        "extraction run": ext.RUN_NAME,
        "topography": str(ext.TOPO_SAVE_PATH),
        "dunes": str(ext.DUNE_SAVE_PATH),
        "DEM domains": str(ext.LOAD_PATH),
        "dune picks": str(ext.WINDOW_JSON),
        "road masks": str(ROAD_MASK_DIR),
        "OCEAN_LOC": ext.OCEAN_LOC,
        "ALONGSHORE_FLIP": ext.ALONGSHORE_FLIP,
        "USE_CONST_INTERIOR": ext.USE_CONST_INTERIOR,
        "TRIM_INTERIOR_ROWS": ext.TRIM_INTERIOR_ROWS,
        "MHW_M": ext.MHW_M,
        "CELL_SIZE_M": ext.CELL_SIZE_M,
        "ROAD_AGG": ROAD_AGG,
        "ROAD_WIDTH_M": ROAD_WIDTH_M,
        "road domains": f"{FIRST_ROAD_DOMAIN}-{LAST_ROAD_DOMAIN}",
    }
    w = max(len(k) for k in cfg)
    lines = [
        "=" * 78,
        f"NC-12 {ROAD_YEAR} setbacks  |  extraction run {ext.RUN_NAME}",
        "=" * 78, "",
        "CONTENTS",
        f"  RoadSetback_{ROAD_YEAR}.csv            2-row CASCADE file "
        f"(IDs, setback m)",
        f"  RoadElevation_{ROAD_YEAR}.csv          2-row, road elevation "
        f"(m MHW-relative)",
        f"  RoadSetback_{ROAD_YEAR}_domains.csv    per-domain stats + QC flags",
        f"  RoadSetback_{ROAD_YEAR}_profiles.csv   every profile measurement",
        f"  HAT_road_setback_qc_{ROAD_YEAR}.png    alongshore QC",
        r"  figures\                          per-domain diagnostic",
        "",
        "THESE SETBACKS BELONG TO THE EXTRACTION RUN THEY SIT INSIDE.",
        "",
        "  They are measured from start_island, which the dune search windows",
        "  in this run define. Re-pick a window and the setback moves with it.",
        "  Point the hindcast at these three together or none of them:",
        "",
        f"    ROAD_SETBACK_FILE = {OUT_DIR / f'RoadSetback_{ROAD_YEAR}.csv'}",
        f"    topography        = {ext.TOPO_SAVE_PATH}",
        f"    dunes             = {ext.DUNE_SAVE_PATH}",
        "",
        "VERTICAL DATUM",
        f"  RoadElevation_{ROAD_YEAR}.csv is MHW-RELATIVE (the extractor",
        f"  subtracts MHW_M = {ext.MHW_M} before anything else), matching the",
        "  frame CASCADE's interior grid is in. _domains.csv also carries",
        "  road_elev_navd if you need NAVD88.",
        "",
        "SETTINGS",
    ]
    lines += [f"  {k:<{w}} : {v}" for k, v in cfg.items()]

    if rows:
        sb = np.array([r["setback_m"] for r in rows
                       if np.isfinite(r["setback_m"])])
        lines += ["", f"DOMAINS: {len(rows)} measured"]
        if sb.size:
            lines.append(f"  setback (m): min {sb.min():.0f} | "
                         f"median {np.median(sb):.0f} | max {sb.max():.0f}")
        for name in ("OFF_BACK", "NEG", "SCATTER", "LOW_ELEV", "THIN",
                     "NO_ROAD"):
            hit = [r["domain"] for r in rows
                   if any(f.startswith(name) for f in flags_for(r))]
            if hit:
                lines.append(f"  {name:<9}: {hit}")
    lines.append("")

    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines))
    print(f"[road] {path}")


def write_elevation_csv(rows, path: Path):
    """
    Per-domain road elevation, 2-row CASCADE format: row 0 = domain IDs,
    row 1 = elevation. Written in the extractor's MHW-RELATIVE frame.

    DATUM -- READ THIS BEFORE USING IT
    ----------------------------------
    HAT_dune_topo_extractor.py subtracts MHW_M (0.36) before anything else, so
    the interior arrays CASCADE runs on are MHW-relative. bulldoze() compares
    road_ele against those arrays. So road_ele must be MHW-relative too.

    The runner currently has ROAD_ELEVATION = 1.45. If that was meant as NAVD88
    it should be 1.09 MHW-relative; if it was already MHW-relative it is 1.81
    NAVD88. One of those is wrong. This file is MHW-relative; the _domains.csv
    carries a road_elev_navd column if you need the other.

    VINTAGE -- THE BIGGER CAVEAT
    ----------------------------
    These elevations come from the 2009 DEM sampled along the ROAD_YEAR
    alignment. Where NC-12 never moved, that is genuinely the roadbed. Where it
    WAS relocated -- the 11% of the line that shifted >50 m between 1978 and
    1997 -- the original alignment in a 2009 DEM is abandoned ground: buried under
    overwash, bulldozed flat, or bare sand. It is not roadbed, and an
    "average road elevation" there is the elevation of a place a road used to
    be.

    So this is trustworthy for the ~84% of domains where the alignment is
    stable, and misleading in exactly the domains you care most about (9-15,
    84-87). Check road_elev_std_m: a domain whose road cells scatter widely is
    one where the mask is crossing a relocation scar, not a graded surface.

    CASCADE's road_ele may be a scalar in your version of Cascade(). Confirm it
    accepts a per-domain list before wiring this in.
    """
    by = {r["domain"]: r for r in rows}
    expected = [d for d in range(FIRST_ROAD_DOMAIN, LAST_ROAD_DOMAIN + 1)
                if d in by and np.isfinite(by[d]["road_elev_mhw"])]
    if not expected:
        print("\n  no road elevations to write")
        return

    ids = np.array(expected, dtype=float)
    ev = np.array([by[d]["road_elev_mhw"] for d in expected], dtype=float)
    path.parent.mkdir(parents=True, exist_ok=True)
    np.savetxt(path, np.vstack([ids, ev]), delimiter=",", fmt="%.3f")
    print(f"[road] {path}")
    print(f"       2 x {len(ids)}  (row 0 = domain IDs, "
          f"row 1 = road elevation, m MHW-relative)")

    print(f"\n  road elevation (m MHW): min {ev.min():.2f} | "
          f"median {np.median(ev):.2f} | max {ev.max():.2f}")
    print(f"  same in NAVD88        : min {ev.min() + ext.MHW_M:.2f} | "
          f"median {np.median(ev) + ext.MHW_M:.2f} | "
          f"max {ev.max() + ext.MHW_M:.2f}")
    print(f"  the runner's ROAD_ELEVATION = 1.45 sits "
          f"{'inside' if ev.min() <= 1.45 <= ev.max() else 'OUTSIDE'} "
          f"this range if read as MHW-relative")

    scat = [d for d in expected
            if np.isfinite(by[d]["road_elev_std_m"])
            and by[d]["road_elev_std_m"] > 0.5]
    if scat:
        print(f"\n  {len(scat)} domain(s) with road cells scattering > 0.5 m sd:")
        print(f"    {scat}")
        print( "    Not a graded surface. Either the mask is crossing a")
        print( "    relocation scar (the ROAD_YEAR alignment abandoned by 2009),")
        print( "    or the buffer is catching the embankment slope.")


def write_cascade_csv(rows, missing, path: Path):
    if missing:
        print(f"\n  NOT WRITING {path.name}: no usable setback for "
              f"domain(s) {missing}.")
        print("  HAT_hindcast_1984_2024_newplot.py fills road_setbacks_full[]")
        print("  BY POSITION and discards the IDs with skiprows=1, so one gap")
        print("  shifts every domain north of it and nothing reports it.")
        return False

    by = {r["domain"]: r for r in rows}
    expected = list(range(FIRST_ROAD_DOMAIN, LAST_ROAD_DOMAIN + 1))

    # A negative setback is not a small error, it is a silent one. bulldoze()
    # does road_start = int(setback / 10), so -110 gives road_start = -11 and
    # xyz_interior_grid[road_end + 1, :] indexes -8 -- which numpy happily
    # wraps to the BACK of the barrier. No exception. The run completes,
    # having bulldozed a road onto the sound side. That is worse than the
    # IndexError, which at least stopped.
    deep = [d for d in expected
            if by[d]["setback_m"] < -SETBACK_FLOOR_TOL_M]
    if deep:
        print(f"\n  NOT WRITING {path.name}: setback below "
              f"-{SETBACK_FLOOR_TOL_M:.0f} m in {len(deep)} domain(s):")
        for d in deep:
            print(f"    GIS {d:>2}: {by[d]['setback_m']:>7.0f} m  "
                  f"(const {by[d]['setback_const_m']:>6.0f}, "
                  f"obliquity {by[d]['obliquity_deg']:.1f} deg)")
        print()
        print("  CASCADE would NOT raise on these. int(-110/10) = -11, and")
        print("  numpy wraps a negative index to the back of the array, so the")
        print("  road would be bulldozed onto the sound side and the run would")
        print("  finish looking normal. Refusing to write.")
        print()
        print(f"  A negative deeper than the mask width is not resolution -- it")
        print(f"  means the dune pick and the road genuinely disagree about")
        print(f"  which is seaward. Look at roads/figures/domain_0NN_setback.png")
        print(f"  and check where the orange dune picks land relative to the")
        print(f"  red road.")
        return False

    ids = np.array(expected, dtype=float)
    sb = np.array([by[d]["setback_m"] for d in expected], dtype=float)

    floored = [(d, by[d]["setback_m"]) for d in expected
               if by[d]["setback_m"] < SETBACK_FLOOR_M]
    if floored:
        sb = np.maximum(sb, SETBACK_FLOOR_M)
        print(f"\n  Floored {len(floored)} setback(s) to "
              f"{SETBACK_FLOOR_M:.0f} m (within {SETBACK_FLOOR_TOL_M:.0f} m, "
              f"i.e. inside the mask's own width):")
        for d, v in floored:
            print(f"    GIS {d:>2}: {v:>6.0f} -> {SETBACK_FLOOR_M:.0f} m")
        print(f"  These are domains where NC-12 runs immediately behind the")
        print(f"  foredune and a ~24 m mask on 10 m cells cannot separate them.")
        print(f"  'Setback 0' says the road is at the dune line, which is what")
        print(f"  the DEM supports. Worth a line in the methods.")
    path.parent.mkdir(parents=True, exist_ok=True)
    np.savetxt(path, np.vstack([ids, sb]), delimiter=",", fmt="%.3f")
    print(f"\n[road] {path}")
    print(f"       2 x {len(ids)}  (row 0 = domain IDs, row 1 = setback m)")
    return True


def write_csvs(rows, out_dir: Path):
    out_dir.mkdir(parents=True, exist_ok=True)

    dom_rows = [{k: v for k, v in r.items() if not k.startswith("_")}
                for r in rows]
    p = out_dir / f"RoadSetback_{ROAD_YEAR}_domains.csv"
    with open(p, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(dom_rows[0].keys()))
        w.writeheader()
        w.writerows(dom_rows)
    print(f"[road] {p}")

    prof = []
    for r in rows:
        for i, v in enumerate(r["_sb"]):
            prof.append(dict(domain=r["domain"], stem=r["stem"], profile=i,
                             setback_origin=r["setback_origin"],
                             road_loc=int(r["_loc"][i]),
                             setback_m=("" if not np.isfinite(v)
                                        else round(float(v), 1))))
    p = out_dir / f"RoadSetback_{ROAD_YEAR}_profiles.csv"
    with open(p, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(prof[0].keys()))
        w.writeheader()
        w.writerows(prof)
    print(f"[road] {p}")


def qc_figure(rows, path: Path):
    rows = sorted([r for r in rows if r["domain"] is not None],
                  key=lambda r: r["domain"])
    if not rows:
        return
    d = np.array([r["domain"] for r in rows])

    fig, (ax0, ax1) = plt.subplots(2, 1, figsize=(13, 8), sharex=True,
                                   gridspec_kw={"hspace": 0.12})

    for k, ((lo, hi), label) in enumerate(ext.SECTIONS):
        for ax in (ax0, ax1):
            ax.axvspan(lo - 0.5, hi + 0.5, color=C_TOWN,
                       alpha=0.18 if k % 2 == 0 else 0.08, lw=0, zorder=0)
        tr = blended_transform_factory(ax0.transData, ax0.transAxes)
        ax0.text((lo + hi) / 2, 0.97, label, transform=tr, ha="center",
                 va="top", fontsize=8,
                 bbox=dict(boxstyle="round,pad=0.2", fc="white", ec="none",
                           alpha=0.85))

    # per-profile scatter behind the per-domain aggregate
    for r in rows:
        v = r["_sb"][np.isfinite(r["_sb"])]
        if v.size:
            ax0.plot(np.full(v.size, r["domain"]), v, lw=0, marker=".", ms=2,
                     color="0.75", zorder=2)
    sb = np.array([r["setback_m"] for r in rows])
    ax0.plot(d, sb, color=C_GROIN, lw=1.6, marker="o", ms=3, zorder=4,
             label=f"setback ({ROAD_AGG})")
    ax0.axhline(0, color=C_PIER, lw=1.0, ls=":", zorder=1)
    ax0.set_ylabel("NC-12 setback from\nInteriorDomain[0]  (m)")
    ax0.legend(loc="upper right", fontsize=8, framealpha=0.9)
    ax0.set_title(f"NC-12 {ROAD_YEAR} setbacks  —  run {ext.RUN_NAME}  —  "
                  f"{len(rows)} domains", fontsize=12)

    # the crash test, drawn
    need = np.array([r["rows_needed"] for r in rows], dtype=float)
    have = np.array([r["interior_rows"] for r in rows], dtype=float)
    ax1.plot(d, have, color=C_PIER, lw=1.6, marker="o", ms=3,
             label="interior rows available")
    ax1.plot(d, need, color=C_MODEL, lw=1.6, marker="s", ms=3,
             label="rows needed (setback + road + check row)")
    bad = need > have
    if bad.any():
        ax1.plot(d[bad], need[bad], lw=0, marker="x", ms=11, mew=2.2,
                 color=C_GROIN, label="OFF_BACK -> IndexError")
    ax1.fill_between(d, have, need, where=bad, color=C_GROIN, alpha=0.2,
                     step=None)
    ax1.set_ylabel("cross-shore rows")
    ax1.set_xlabel("GIS domain  (1 = Cape Point / south  ->  90 = Rodanthe / north)")
    ax1.legend(loc="upper right", fontsize=8, framealpha=0.9)
    ax1.set_xlim(d.min() - 0.5, d.max() + 0.5)

    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=140, bbox_inches="tight")
    plt.close(fig)
    print(f"[road] {path}")


# =============================================================================
# MAIN
# =============================================================================

def main():
    print("=" * 86)
    print(f"NC-12 {ROAD_YEAR} SETBACK EXTRACTION")
    print("=" * 86)
    print(f"  extractor  : {Path(ext.__file__).name}")
    print(f"  run        : {ext.RUN_NAME}")
    print(f"  DEM domains: {ext.LOAD_PATH}")
    print(f"  picks      : {ext.WINDOW_JSON}")
    print(f"  road masks : {ROAD_MASK_DIR}")
    print(f"  out        : {OUT_DIR}")
    print(f"  settings   : OCEAN_LOC={ext.OCEAN_LOC}, "
          f"ALONGSHORE_FLIP={ext.ALONGSHORE_FLIP}, "
          f"USE_CONST_INTERIOR={ext.USE_CONST_INTERIOR}, "
          f"TRIM_INTERIOR_ROWS={ext.TRIM_INTERIOR_ROWS}")
    print()

    if not ext.USE_CONST_INTERIOR:
        print("[stop] USE_CONST_INTERIOR is False. With per-profile interiors")
        print("       there is no single row 0 for a setback to be measured")
        print("       from. Set it True in the extractor and re-extract.")
        return
    if not Path(ROAD_MASK_DIR).is_dir():
        print(f"[stop] no road masks at {ROAD_MASK_DIR}")
        print("       run HAT_rasterize_road_to_domains.py first")
        return
    if not Path(ext.WINDOW_JSON).exists():
        print(f"[stop] no picks at {ext.WINDOW_JSON}")
        return

    windows = ext.load_windows(Path(ext.WINDOW_JSON))
    names = sorted([n for n in os.listdir(ext.LOAD_PATH)
                    if n.startswith("domain_") and n.endswith(".npy")],
                   key=ext.natural_key)
    print(f"[info] {len(names)} domain file(s)\n")

    rows = []
    n_failed = 0
    for name in names:
        stem = Path(name).stem
        d = ext.domain_number(stem)
        if d is None or not (FIRST_ROAD_DOMAIN <= d <= LAST_ROAD_DOMAIN):
            continue
        try:
            r = process(stem, Path(ext.LOAD_PATH) / name, windows)
        except Exception as e:
            # A per-domain skip is fine for a genuinely bad domain, but the same
            # message on every domain means a code fault, not data. Show it.
            n_failed += 1
            print(f"[skip] {stem}: {type(e).__name__}: {e}")
            if n_failed == 1:
                import traceback
                traceback.print_exc()
                print("  ^ first failure shown in full. If every domain fails "
                      "with this, it is a bug in the script, not your data.")
            continue
        if r is not None:
            rows.append(r)

    if not rows:
        print("[stop] nothing extracted")
        return

    missing = report(rows)
    write_csvs(rows, Path(OUT_DIR))
    qc_figure(rows, Path(OUT_DIR) / f"HAT_road_setback_qc_{ROAD_YEAR}.png")

    if SAVE_DOMAIN_FIGS:
        fig_dir = Path(OUT_DIR) / "figures"
        want = (QC_DOMAINS if QC_DOMAINS is not None
                else [r["domain"] for r in rows])
        by = {r["domain"]: r for r in rows}
        wrote = []
        for d in want:
            if d not in by:
                continue
            try:
                domain_setback_figure(by[d], fig_dir)
                wrote.append(d)
            except Exception as e:
                print(f"  [fig] domain {d}: {type(e).__name__}: {e}")
        print(f"[road] {fig_dir}")
        if wrote:
            print(f"       wrote {len(wrote)} of {len(want)}: "
                  + ", ".join(f"domain_{d:03d}_setback.png" for d in wrote))
    write_elevation_csv(rows,
                        Path(OUT_DIR) / f"RoadElevation_{ROAD_YEAR}.csv")
    write_manifest(rows, Path(OUT_DIR) / "RUN_MANIFEST.txt")
    ok = write_cascade_csv(rows, missing,
                           Path(OUT_DIR) / f"RoadSetback_{ROAD_YEAR}.csv")

    if ok:
        print("\n" + "=" * 86)
        print("NEXT: point HAT_hindcast_1984_2024_newplot.py at this run")
        print("=" * 86)
        print(f"  ROAD_SETBACK_FILE = r\"{OUT_DIR / f'RoadSetback_{ROAD_YEAR}.csv'}\"")
        print(f"  topography        = r\"{ext.TOPO_SAVE_PATH}\"")
        print(f"  dunes             = r\"{ext.DUNE_SAVE_PATH}\"")
        print()
        print("  All three must come from the SAME run folder. A setback from")
        print("  one extraction against topography from another is the frame")
        print("  mismatch this script exists to prevent.")


if __name__ == "__main__":
    main()
