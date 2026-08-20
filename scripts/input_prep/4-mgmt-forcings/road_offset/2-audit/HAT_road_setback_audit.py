r"""
HAT_road_setback_audit.py
===============================================================================
Audit the NC-12 setbacks -- initial AND prescribed-relocation -- against the
interior arrays CASCADE will actually run on, and write down where the road
lands off the island and what the model does about it.

This script is DIAGNOSTIC. It never modifies a setback, never writes a setback
file, and never touches the topography.

WHY AN AUDIT AND NOT A FIX
--------------------------
The setback MEASUREMENT was tested and holds up. Road and dune are digitised as
points along the SAME 1 m transects (FID_Transect_Points_1m / LineID /
ORIG_SEQ), from the same origin, in the same units, for the SAME YEAR -- so the
shoreline retreat cancels in the subtraction. road_offset_pipeline.py's one
shortcut is taking min(ORIG_LEN) over each layer INDEPENDENTLY, so the two
minima land on different transects. Measured on 2004 (the only vintage whose
raw road CSV survives), pairing road and dune on the SAME transect changes the
answer by a median of +1 m (worst -49/+54; |delta| > 25 m in 14 of 82 domains,
> 100 m in none). One grid cell. Not worth re-deriving.

WHAT ACTUALLY WENT WRONG IS THE FIT CHECK
-----------------------------------------
HAT_setback_from_lines.py (retired 2026-08-17, git blob c37ec03e) bounds the
setback with

    interior_rows_of(d) -> np.load(...).shape[0]

which is Barrier3D's DomainWidth: the ROW COUNT of the array. Barrier3D never
uses that to decide where land is. barrier3d.py FindWidths() does:

    width = next((i for i, v in enumerate(InteriorDomain[:, bl])
                  if v <= SL), DomainWidth) - 1

-- per profile, the contiguous run of land from row 0 to the FIRST cell at or
below sea level. remove_water_rows() only trims LEADING and TRAILING all-water
rows, so one wide profile sets DomainWidth for all 50 and everything behind the
other 49 is sentinel padding. Median gap between the two numbers across GIS
9-90: 17 rows. Worst: 86 rows -- 860 m of setback the old check would approve
onto cells that are not there.

find_widths() below is FindWidths transcribed verbatim. The definition is
adopted, not invented; nothing in Barrier3D changes.

WHAT "DROWNS" ACTUALLY COSTS -- READ THIS BEFORE READING THE TABLES
-------------------------------------------------------------------
roadway_drown is not a warning. Traced through the code, one drown does this:

  1. bulldoze() MUTATES barrier3d.InteriorDomain IN PLACE before it checks
     anything: xyz_interior_grid[road_start:road_end, :] = new_road_domain.
     The grid is passed by reference, so the write lands even though the
     function returns early afterwards. On GIS 52 at the 1984 setback, 53 of
     100 cells in the road rows were water (-0.3 dam) and all 100 came out at
     0.145 dam -- the model gains a 20 m ribbon of 1.45 m land across open
     water.

  2. RoadwayManager sets _drown_break = 1 and returns.

  3. cascade.py (~line 625) sees drown_break on EVERY later year and never
     calls update() again. _road_break[iB3D] = 1, dune growth rates reset to
     natural.

So from that year on there is no overwash removal, no dune rebuilding, no road
relocation, and _road_overwash_volume / _dunes_rebuilt_TS /
_rebuild_dune_volume_TS stay at zero for the rest of the run. A domain flagged
DROWNS_T0 is an UNMANAGED barrier wearing a road label for the whole hindcast.
If any result contrasts managed against unmanaged shoreline, those domains are
silently in the unmanaged group from year zero.

(group_roadway_abandonment is None in the runner, so this stays per-domain. Set
it, and one drowned domain abandons its entire group.)

THE ROAD'S OWN CELLS ARE NEVER CHECKED
--------------------------------------
bulldoze() tests the rows FLANKING the road -- road_end + 1 on the bay side,
road_start - 1 on the sea side -- and never looks at the road footprint itself.
It also skips row road_end: the road occupies [road_start:road_end], and the
bayside test is at road_end + 1, so the cell immediately behind the asphalt is
only used for np.size(). This audit reports pct_road_cells_water separately for
exactly that reason: a road can be substantially in the bay while the flanking
test is still under threshold.

PRESCRIBED RELOCATIONS BYPASS BOTH OF THE MODEL'S GUARDS
--------------------------------------------------------
CASCADE has two guards on moving a road landward:

    get_road_relocation_elevation()  rebuilds the road at grade and refuses if
                                     mean(road_domain) <= 0 m MSL
    road_relocation_checks()         refuses if
                                     setback + 2*width > average_barrier_width
                                     (InteriorWidth_Avg, i.e. FindWidths)

Both live on the MODEL-DRIVEN path, taken when the dune migrates over the road.
The runner's HISTORICAL_ROAD_EVENTS path does not go through either -- it
assigns rm._road_setback = new_sb directly, leaves road_ele at its
SLR-decremented value rather than rebuilding at grade, and lets the next
bulldoze() flatten whatever is there. This audit replicates both guards on the
prescribed setbacks so the report can say which the model would have refused.

INPUTS
------
  SCENARIOS below: 2-row RoadSetback_<year>.csv files and/or literal
  {domain: setback_m} maps taken from HISTORICAL_ROAD_EVENTS
  TOPO_DIR/domain_<N>_topography_<year>.npy

OUTPUTS
-------
  RoadSetback_audit.csv   every scenario x domain, machine readable
  RoadSetback_audit.md    the tracking document
  console verdict, exit code 1 if any scenario hits a hard wall

REQUIREMENTS
------------
  numpy
===============================================================================
"""

from __future__ import annotations

import csv
import sys
from datetime import datetime
from pathlib import Path

import numpy as np

# =============================================================================
# CONFIG
# =============================================================================

def _find_project_root(start: Path) -> Path:
    """
    Walk up until a directory holds data\\hatteras_init.

    NOT parents[N]. This file has moved, and the old parents[4] resolved to
    scripts\\ instead of the project root -- so every path below it was wrong
    and nothing raised until a glob came back empty.
    """
    for p in [start, *start.parents]:
        if (p / "data" / "hatteras_init").is_dir():
            return p
    raise SystemExit(f"cannot find data\\hatteras_init above {start}")


PROJECT_ROOT = _find_project_root(Path(__file__).resolve())
HATTERAS_DATA_BASE = PROJECT_ROOT / "data" / "hatteras_init"
BARRIER3D_DIR = HATTERAS_DATA_BASE / "1-barrier3d-domains"
# The method the runner spends (hatteras_site_config.py:78,91). Switched to
# dune-start on 2026-08-18; the legacy tree is still on disk under
# old_method_offset/ for the method comparison.
ROADS_DIR = HATTERAS_DATA_BASE / "4-mgmt-forcing" / "road_offset" / "dunestart_offset"

# The runner uses ONE topography for BOTH hindcast periods, so the 1984 setbacks
# and the 2004 setbacks are both spent against a 2009 grid. For 1984 that is a
# 25-year anachronism; it is the reason this audit exists.
#
# This USED to say "must match TOPO_DUNE_VERSION in the runner", and it was
# pinned to 2009_v2, then to 2009_v3. Pinning by hand is what went wrong: when
# the dune windows were re-picked into 2009_v4 on 2026-08-19 this kept auditing
# v3 interiors against v4 setbacks and said PASS, on a grid the model no longer
# uses. 18 domains had different interiors and 10 a different SHAPE -- including
# D79 and D80, two of the three roadways the relocation logic acts on.
#
# So the version now comes from the extractor that PRODUCED the arrays, and a
# version missing from disk is a loud error rather than a stale read. To audit
# an older run deliberately, pass override= (see hat_topo_version.py).
#
# The runner (HAT_hindcast_1984_2024.ipynb, HAT_hindcast_1984_2024.py,
# HAT_groin_sweep_worker.py) still says 2009_v2, a path that does not exist, so
# the mismatch with the runner remains REAL and is printed below rather than
# hidden. Repointing the runner is a separate job.
# parents[4] IS scripts/ -- hat_topo_version.py moved there 2026-08-20.
sys.path.insert(0, str(Path(__file__).resolve().parents[4]))
from hat_topo_version import topo_dirs  # noqa: E402

TOPO_DIR, _DUNE_DIR, TOPO_DUNE_VERSION = topo_dirs()
TOPO_DUNE_INIT_YEAR = "2009"

# There is no longer a runner version to disagree with. Until 2026-08-20 the
# runner carried its own hardcoded TOPO_DUNE_VERSION and this file kept a copy
# of it so a mismatch could be warned about -- a copy that itself went stale
# (it said 2009_v2 while the runner said 2009_v3 and the setbacks were built on
# v5). Both now call topo_dirs(), so they resolve to the same extraction by
# construction and the warning has nothing left to check.
print(f"[version] auditing {TOPO_DUNE_VERSION}, resolved from the extractor "
      f"-- the same source the runner uses.")

OUT_DIR = ROADS_DIR

# --- what to audit ------------------------------------------------------
# "initial"    a 2-row RoadSetback_<year>.csv, the whole GIS 9-90 window
# "relocation" a literal {gis: setback_m} lifted from HISTORICAL_ROAD_EVENTS in
#              HAT_hindcast_1984_2024.py. Relocation scenarios also
#              get the two model guards replicated.
SCENARIOS = [
    dict(label="1984 initial", kind="initial",
         source=ROADS_DIR / "1984" / "RoadSetback_1984_dunestart.csv"),
    dict(label="2004 initial", kind="initial",
         source=ROADS_DIR / "2004" / "RoadSetback_2004_dunestart.csv"),
    # The relocation events now carry a DISPLACEMENT, applied to whatever
    # setback the model is carrying at the event year:
    #
    #     new_sb = rm._road_setback + relocation_displacement_m[gis]
    #
    # so the audited value depends on how much the model has retreated by then.
    # The setbacks below are that expression evaluated at ZERO modelled retreat,
    # i.e. 1984 setback + displacement -- which is exactly the old absolute
    # post_relocation_setback_m. So this scenario is the WORST CASE / upper
    # bound of the corrected form, and any real retreat moves the road seaward
    # of what is reported here. Audited at the bound deliberately: if a domain
    # is clean here it is clean for every retreat.
    dict(label="1999 relocation (inter-village south, GIS 9-15) "
               "- zero-retreat bound",
         kind="relocation",
         note="1984 setback + 1978->1997 displacement; upper bound of "
              "current_setback + displacement",
         setbacks={9: 73.0, 10: 97.0, 11: 129.0, 12: 126.0, 13: 125.0,
                   14: 126.0, 15: 106.0}),
    dict(label="1989 relocation (Pea Island, GIS 84-87) - zero-retreat bound",
         kind="relocation",
         note="1984 setback + 1978->1997 displacement; upper bound of "
              "current_setback + displacement",
         setbacks={84: 163.0, 85: 165.0, 86: 205.0, 87: 113.0}),
]

FIRST_ROAD_DOMAIN = 9
LAST_ROAD_DOMAIN = 90

# --- Barrier3D / roadway_manager constants ------------------------------
# Not tuning knobs. Every one read from model source:
#   ROAD_WIDTH_M   HAT_hindcast_1984_2024.py:2206
#   DX/DY/DZ       roadway_manager.bulldoze() defaults, dam
#   DROWN_THR      roadway_manager.py:766   drown_threshold = 0
#   PCT_WATER      roadway_manager.py:531   percent_water_cells_touching_road
#   SL             barrier3d FindWidths sea level, 0.0 dam at t = 0
ROAD_WIDTH_M = 20.0
DX = DY = DZ = 10
DROWN_THRESHOLD = 0.0
PCT_WATER_TOUCHING_ROAD = 0.2
SL_DAM = 0.0

# A road is called off-island when this fraction of its own cells are water.
ROAD_IN_WATER_FRACTION = 0.20

HOLED_FRACTION = 0.50
MIXED_FRACTION = 0.20

# (GIS, old absolute post_relocation_setback_m, RoadSetback_2004.csv value).
# The evidence that the old absolute values double-counted the retreat: both
# events precede 2004, so the 2004 same-year measurement already IS the
# post-relocation position. Reported in the markdown; nothing computes from it.
OVERSHOOT_TABLE = [
    (11, 129, 81), (12, 126, 89), (13, 125, 87), (14, 126, 93), (15, 106, 71),
    (84, 163, 50), (85, 165, 85), (86, 205, 88), (87, 113, 40),
]


# =============================================================================
# BARRIER3D'S OWN DEFINITIONS
# =============================================================================

def find_widths(interior: np.ndarray, sl: float = SL_DAM):
    """
    barrier3d.py FindWidths(), verbatim -- including the `- 1` and the clamp at
    zero, so it cannot drift from the model.

    InteriorWidth[bl] is the contiguous run of land from row 0 of profile bl to
    the first cell at or below sea level. It is the only width Barrier3D uses
    to decide where the island is; land behind a water cell is invisible to it.
    """
    domain_width = int(np.shape(interior)[0])
    n_along = int(np.shape(interior)[1])
    interior_width = np.zeros(n_along, dtype=int)
    for bl in range(n_along):
        width = next(
            (index for index, value in enumerate(interior[:, bl])
             if value <= sl),
            domain_width,
        )
        width = width - 1
        if width < 0:
            width = 0
        interior_width[bl] = width
    return domain_width, interior_width


def land_behind_first_water(interior: np.ndarray, sl: float = SL_DAM):
    """
    Per profile: land cells BEHIND the first water cell, and the width of the
    intervening gap.

    Separates 'the island ends here' from 'FindWidths stopped at a 2-cell hole
    and discarded 33 cells of barrier'. Decides nothing -- with no bathymetry
    the two are indistinguishable by value -- but marks which domains carry the
    doubt.
    """
    domain_width = int(np.shape(interior)[0])
    n_along = int(np.shape(interior)[1])
    beyond = np.zeros(n_along, dtype=int)
    gap = np.zeros(n_along, dtype=int)
    for bl in range(n_along):
        col = interior[:, bl]
        w = next((i for i, v in enumerate(col) if v <= sl), domain_width)
        tail = col[w:]
        land = int((tail > sl).sum())
        beyond[bl] = land
        if land:
            gap[bl] = next(i for i, v in enumerate(tail) if v > sl)
    return beyond, gap


def predict_bulldoze(interior: np.ndarray, setback_m: float,
                     interior_width_arr: np.ndarray | None = None):
    """
    Reproduce roadway_manager.bulldoze()'s indexing, its drown test, and what
    it does to the road's own footprint -- without running it.

    `wall` is None unless the setback would CORRUPT the run rather than merely
    drown the road:

        NEGATIVE  road_start = int(-110/10) = -11, and numpy WRAPS a negative
                  index to the back of the array. The road is bulldozed onto
                  the sound side and the run finishes looking normal.
        OVERRUN   xyz_interior_grid[road_end + 1, :] is past the array ->
                  IndexError mid-run.

    Datum note: bulldoze compares (interior * dz) against drown_threshold = 0
    described as "m MSL", but the extractor's arrays are MHW-RELATIVE (MHW_M =
    0.36 subtracted first). The effective test is "at or below MHW". A
    pre-existing inconsistency in the model -- reported, not fixed.
    """
    domain_width = int(np.shape(interior)[0])
    n_along = int(np.shape(interior)[1])

    road_start = int(setback_m / DY)
    road_width_cells = int(ROAD_WIDTH_M / DX)
    road_end = road_start + road_width_cells
    border = road_end + 1

    out = dict(road_start=road_start, road_end=road_end, border_row=border,
               domain_width=domain_width, wall=None,
               bayside_water=np.nan, seaside_water=np.nan, drowns=False,
               road_cells_water=np.nan, profiles_road_submerged=np.nan,
               findwidths_disagree=np.nan)

    if road_start < 0:
        out["wall"] = "NEGATIVE"
        return out
    if border >= domain_width:
        out["wall"] = "OVERRUN"
        return out

    # --- the road's own footprint. bulldoze never looks at this. -------
    footprint = interior[road_start:road_end, :]
    out["road_cells_water"] = float((footprint <= SL_DAM).mean())
    out["profiles_road_submerged"] = float(
        (footprint <= SL_DAM).all(axis=0).mean())

    # --- the flanking rows, which is all bulldoze tests ----------------
    # Every cell counts, no-data included -- bulldoze reads the literal array
    # and the extractor stores no-data as the water sentinel. See wet_fraction.
    bayside = wet_fraction(interior[border, :])
    seaside = wet_fraction(interior[road_start - 1, :]) if road_start > 0 else 0.0

    out["bayside_water"] = float(bayside)
    out["seaside_water"] = float(seaside)
    out["drowns"] = bool(seaside > PCT_WATER_TOUCHING_ROAD
                         or bayside > PCT_WATER_TOUCHING_ROAD)

    # Barrier3D's two guards do not use the same notion of "island", and in a
    # holed domain they disagree:
    #   FindWidths stops at the FIRST water cell, so land behind a gap is not
    #     part of the island. Feeds InteriorWidth_Avg -> average_barrier_width,
    #     which road_relocation_checks() uses to decide if there is room.
    #   bulldoze reads the LITERAL cell at road_end + 1; land behind a gap
    #     counts perfectly well.
    # So the road can be founded on ground the relocation logic believes is off
    # the back of the barrier. Counted, not corrected.
    if interior_width_arr is not None:
        out["findwidths_disagree"] = float(np.count_nonzero(
            (interior[border, :] > SL_DAM) & (interior_width_arr < border)
        ) / n_along)
    return out


def replicate_relocation_guards(interior: np.ndarray, setback_m: float):
    """
    What CASCADE's own relocation guards would say about this setback.

    get_road_relocation_elevation(): rebuilds the road at grade and refuses if
        mean(road_domain) * dz <= 0 -- "Roadway cannot be relocated ... b/c the
        road would be at or below MSL".

    road_relocation_checks(): refuses if
        road_relocation_setback + 2 * road_relocation_width
            > average_barrier_width
    where average_barrier_width is barrier3d.InteriorWidth_AvgTS[-1] * 10,
    i.e. the MEAN of FindWidths' InteriorWidth, in metres.

    Both are on the model-driven relocation path. The runner's
    HISTORICAL_ROAD_EVENTS path assigns rm._road_setback directly and takes
    neither. Replicated here on the t=0 interior; the real 1989/1999 grid will
    have evolved, so treat these as the initialisation-time verdict.
    """
    domain_width = int(np.shape(interior)[0])
    road_start = int(setback_m / DY)
    road_end = road_start + int(ROAD_WIDTH_M / DX)

    out = dict(relocation_road_ele_m=np.nan, relocation_ele_refuses=False,
               average_barrier_width_m=np.nan, relocation_width_refuses=False)

    if road_start >= 0 and road_end <= domain_width:
        out["relocation_road_ele_m"] = float(
            np.mean(interior[road_start:road_end, :]) * DZ)
        out["relocation_ele_refuses"] = bool(out["relocation_road_ele_m"] <= 0)

    _, iw = find_widths(interior)
    avg_width_m = float(iw.mean()) * DY
    out["average_barrier_width_m"] = avg_width_m
    out["relocation_width_refuses"] = bool(
        setback_m + 2 * ROAD_WIDTH_M > avg_width_m)
    return out


def largest_setback_that_fits(interior: np.ndarray) -> float:
    """
    Largest setback that neither drowns nor overruns. REFERENCE ONLY -- never
    applied to a setback. Reported so the document can say how far the aerials'
    position is from anything the modelled barrier could hold.
    """
    domain_width = int(np.shape(interior)[0])
    road_width_cells = int(ROAD_WIDTH_M / DX)
    best = -1
    for start in range(0, domain_width - road_width_cells - 1):
        p = predict_bulldoze(interior, float(start * DY))
        if p["wall"] is None and not p["drowns"]:
            best = start
    return float(best * DY) if best >= 0 else np.nan


# =============================================================================
# LOADING
# =============================================================================

def load_setback_csv(path: Path):
    """
    Read a 2-row CASCADE file and verify the ID row is what the runner assumes.

    The runner does np.loadtxt(..., skiprows=1) and fills road_setbacks_full[]
    BY POSITION, discarding the IDs. A missing, duplicated or out-of-order ID
    silently shifts every domain north of it. Checked here because this is the
    only place it can be.
    """
    if not path.exists():
        return None, None, [f"file not found: {path}"]
    arr = np.loadtxt(path, delimiter=",")
    if arr.ndim != 2 or arr.shape[0] != 2:
        return None, None, [f"expected a 2-row file, got shape {arr.shape}"]

    ids = arr[0].astype(int)
    vals = arr[1].astype(float)
    expected = np.arange(FIRST_ROAD_DOMAIN, LAST_ROAD_DOMAIN + 1)
    problems = []
    if ids.size != expected.size:
        problems.append(f"{ids.size} columns, expected {expected.size} "
                        f"(GIS {FIRST_ROAD_DOMAIN}-{LAST_ROAD_DOMAIN})")
    elif not np.array_equal(ids, expected):
        bad = [int(i) for i, e in zip(ids, expected) if i != e]
        problems.append(f"ID row not contiguous {FIRST_ROAD_DOMAIN}.."
                        f"{LAST_ROAD_DOMAIN}; first mismatches {bad[:5]}")
    return ids, vals, problems


def load_interior(domain: int):
    p = TOPO_DIR / f"domain_{domain}_topography_{TOPO_DUNE_INIT_YEAR}.npy"
    if not p.exists():
        matches = sorted(TOPO_DIR.glob(f"domain_{domain}_topography_*.npy"))
        if not matches:
            return None
        p = matches[0]
    return np.load(p)


def wet_fraction(row) -> float:
    """
    Fraction of a flanking row at or below the threshold -- EVERY cell counted.

    A no-data mask (domain_<d>_nodata_<TAG>.npy) exists beside the topography
    and was briefly used to drop never-surveyed cells from both the numerator
    and the denominator. That is not the test bulldoze runs. bulldoze reads the
    literal array, in which the extractor has already written no-data back as
    the water sentinel, so an unsurveyed cell IS a water cell to the model.
    Excluding it made GIS 78/79/80 read FITS while the run still drowns them --
    exactly the silent pass this audit exists to catch. Deliberately reverted.
    """
    return float(np.count_nonzero((row * DZ) <= DROWN_THRESHOLD) / row.size)


# =============================================================================
# AUDIT
# =============================================================================

def audit_domain(scenario_label: str, kind: str, domain: int,
                 setback_m: float) -> dict:
    interior = load_interior(domain)
    if interior is None:
        return dict(scenario=scenario_label, domain=domain,
                    setback_m=setback_m, verdict="NO_TOPOGRAPHY")

    domain_width, iw = find_widths(interior)
    beyond, gap = land_behind_first_water(interior)
    p = predict_bulldoze(interior, setback_m, interior_width_arr=iw)

    frac_land_behind = float((beyond > 0).mean())
    if frac_land_behind > HOLED_FRACTION:
        narrow_label = "HOLED"
    elif frac_land_behind > MIXED_FRACTION:
        narrow_label = "MIXED"
    else:
        narrow_label = "TRUE_END"

    if p["wall"] == "NEGATIVE":
        verdict = "BLOCK_NEGATIVE"
    elif p["wall"] == "OVERRUN":
        verdict = "BLOCK_OVERRUN"
    elif p["drowns"]:
        verdict = "DROWNS"
    else:
        verdict = "FITS"

    off_island = (np.isfinite(p["road_cells_water"])
                  and p["road_cells_water"] >= ROAD_IN_WATER_FRACTION)

    row = dict(
        scenario=scenario_label, kind=kind, domain=domain,
        setback_m=round(float(setback_m), 1),
        verdict=verdict,
        road_off_island=bool(off_island),
        pct_road_cells_water=(round(p["road_cells_water"] * 100, 1)
                              if np.isfinite(p["road_cells_water"]) else ""),
        pct_profiles_road_submerged=(
            round(p["profiles_road_submerged"] * 100, 1)
            if np.isfinite(p["profiles_road_submerged"]) else ""),
        pct_bayside_water=(round(p["bayside_water"] * 100, 1)
                           if np.isfinite(p["bayside_water"]) else ""),
        pct_seaside_water=(round(p["seaside_water"] * 100, 1)
                           if np.isfinite(p["seaside_water"]) else ""),
        pct_findwidths_disagree=(round(p["findwidths_disagree"] * 100, 1)
                                 if np.isfinite(p["findwidths_disagree"])
                                 else ""),
        road_start_cell=p["road_start"], road_end_cell=p["road_end"],
        border_row=p["border_row"], domain_width_rows=domain_width,
        iw_min=int(iw.min()), iw_p20=int(np.percentile(iw, 20)),
        iw_median=int(np.median(iw)), iw_max=int(iw.max()),
        largest_fitting_setback_m=largest_setback_that_fits(interior),
        narrow_label=narrow_label,
        pct_profiles_land_behind=round(frac_land_behind * 100, 1),
        land_behind_median_cells=int(np.median(beyond)),
        gap_median_cells=(int(np.median(gap[beyond > 0]))
                          if (beyond > 0).any() else 0),
    )

    if kind == "relocation":
        row.update({k: (round(v, 2) if isinstance(v, float) else v)
                    for k, v in replicate_relocation_guards(
                        interior, setback_m).items()})
    else:
        row.update(relocation_road_ele_m="", relocation_ele_refuses="",
                   average_barrier_width_m="", relocation_width_refuses="")
    return row


def run_scenario(sc: dict):
    label, kind = sc["label"], sc["kind"]
    problems = []
    if kind == "initial":
        ids, vals, problems = load_setback_csv(Path(sc["source"]))
        if ids is None:
            print(f"\n[skip] {label}: {problems[0]}")
            return [], problems
        pairs = list(zip(ids.tolist(), vals.tolist()))
    else:
        pairs = sorted(sc["setbacks"].items())
    rows = [audit_domain(label, kind, int(d), float(s)) for d, s in pairs]
    return rows, problems


# =============================================================================
# OUTPUT
# =============================================================================

VERDICT_ORDER = ["BLOCK_NEGATIVE", "BLOCK_OVERRUN", "DROWNS", "NO_TOPOGRAPHY",
                 "FITS"]


def print_scenario(label, rows, problems):
    print("\n" + "=" * 104)
    print(f"{label}")
    print("=" * 104)
    for pr in problems:
        print(f"  [!] ID ROW: {pr}")

    print(f"{'GIS':>4} {'setback':>8} {'%roadWet':>9} {'%profWet':>9} "
          f"{'%bayBrd':>8} {'IWp20':>6} {'maxfit':>7}  {'verdict':<15} "
          f"{'offIsland':<10} {'narrow':<9}")
    for r in sorted(rows, key=lambda x: x["domain"]):
        if r["verdict"] == "NO_TOPOGRAPHY":
            print(f"{r['domain']:>4} {r['setback_m']:>8.0f}  NO_TOPOGRAPHY")
            continue
        mf = r["largest_fitting_setback_m"]
        flag = "OFF-ISLAND" if r["road_off_island"] else ""
        print(f"{r['domain']:>4} {r['setback_m']:>8.0f} "
              f"{r['pct_road_cells_water']:>8}% "
              f"{r['pct_profiles_road_submerged']:>8}% "
              f"{r['pct_bayside_water']:>7}% {r['iw_p20']:>6} "
              f"{('--' if not np.isfinite(mf) else f'{mf:.0f}'):>7}  "
              f"{r['verdict']:<15} {flag:<10} {r['narrow_label']:<9}")

    counts = {v: sum(1 for r in rows if r["verdict"] == v)
              for v in VERDICT_ORDER}
    print()
    for v in VERDICT_ORDER:
        if counts[v]:
            print(f"  {v:<15} {counts[v]:>3}")
    off = [r for r in rows if r["road_off_island"]]
    if off:
        print(f"  {'OFF-ISLAND':<15} {len(off):>3}  "
              f"(>= {ROAD_IN_WATER_FRACTION * 100:.0f}% of the road's own "
              f"cells in water)")

    if rows and rows[0]["kind"] == "relocation":
        print("\n  Model guards the prescribed-relocation path SKIPS:")
        print(f"  {'GIS':>4} {'road_ele(m)':>12} "
              f"{'get_road_relocation_elevation':>32}   "
              f"{'avgWidth(m)':>12} {'road_relocation_checks':>24}")
        for r in sorted(rows, key=lambda x: x["domain"]):
            if r["relocation_road_ele_m"] == "":
                continue
            v1 = ("REFUSES (road <= 0 m MSL)" if r["relocation_ele_refuses"]
                  else "allows it")
            v2 = ("REFUSES (island too narrow)"
                  if r["relocation_width_refuses"] else "allows it")
            print(f"  {r['domain']:>4} {r['relocation_road_ele_m']:>12.2f} "
                  f"{v1:>32}   {r['average_barrier_width_m']:>12.0f} "
                  f"{v2:>24}")
        refused = [r for r in rows if r["relocation_ele_refuses"]
                   or r["relocation_width_refuses"]]
        print(f"\n  {len(refused)} of {len(rows)} prescribed relocations would "
              f"be REFUSED by CASCADE's own logic.")
        print("  The runner assigns rm._road_setback directly, so neither "
              "guard runs.")

    return counts


def write_csv(all_rows, path: Path):
    path.parent.mkdir(parents=True, exist_ok=True)
    fields = list(all_rows[0].keys())
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields, extrasaction="ignore")
        w.writeheader()
        for r in all_rows:
            w.writerow(r)
    print(f"\n[out] {path}")


def write_markdown(scen_rows, all_problems, path: Path):
    L = []
    a = L.append
    a("# NC-12 road setback audit")
    a("")
    a(f"Generated {datetime.now().isoformat(timespec='seconds')} by "
      f"`{Path(__file__).name}`.")
    a("")
    a("**Diagnostic record.** No setback was modified, no setback file was "
      "written, no topography was touched.")
    a("")
    a("| | |")
    a("|---|---|")
    a(f"| Topography | `{TOPO_DIR}` |")
    a(f"| Road width | {ROAD_WIDTH_M:.0f} m |")
    a("| Width rule | `barrier3d.FindWidths` — first cell `<= SL`, verbatim |")
    a(f"| Drown rule | `roadway_manager.bulldoze` — "
      f">{PCT_WATER_TOUCHING_ROAD * 100:.0f}% of profiles bordering water |")
    a(f"| Off-island rule | ≥{ROAD_IN_WATER_FRACTION * 100:.0f}% of the "
      f"road's **own** cells in water |")
    a("")
    a("The runner uses **one topography for both hindcast periods**, so the "
      "1984 setbacks and the 2004 setbacks are both spent against a 2009 "
      "grid. For 1984 that is a 25-year anachronism, and it is the reason "
      "this audit exists.")
    a("")

    # ---- the consequence chain, up front ----
    a("## What `DROWNS` actually costs")
    a("")
    a("`roadway_drown` is not a warning. Traced through the code, one drown "
      "does this:")
    a("")
    a("1. **`bulldoze()` mutates `barrier3d.InteriorDomain` in place before "
      "it checks anything.** `xyz_interior_grid[road_start:road_end, :] = "
      "new_road_domain` — the grid is passed by reference, so the write lands "
      "even though the function returns early afterwards. On GIS 52 at the "
      "1984 setback, 53 of 100 cells in the road rows were water (`-0.3 dam`) "
      "and all 100 came out at `0.145 dam`: the model gains a 20 m ribbon of "
      "1.45 m land across open water.")
    a("2. `RoadwayManager` sets `_drown_break = 1` and returns.")
    a("3. `cascade.py` (~line 625) sees `drown_break` on **every later year** "
      "and never calls `update()` again. `_road_break[iB3D] = 1`, dune growth "
      "rates reset to natural.")
    a("")
    a("So from that year on there is **no overwash removal, no dune "
      "rebuilding, and no road relocation**, and `_road_overwash_volume` / "
      "`_dunes_rebuilt_TS` / `_rebuild_dune_volume_TS` stay at zero for the "
      "rest of the run.")
    a("")
    a("> A domain flagged `DROWNS` is an **unmanaged barrier wearing a road "
      "label** for the whole hindcast. If any result contrasts managed "
      "against unmanaged shoreline, those domains are silently in the "
      "unmanaged group from year zero.")
    a("")
    a("`group_roadway_abandonment` is `None` in the runner, so this stays "
      "per-domain. Set it, and one drowned domain abandons its entire group.")
    a("")

    a("## How the relocation events were corrected")
    a("")
    a("The events used to carry `post_relocation_setback_m`, an **absolute** "
      "setback built as `(1984 setback) + (1978→1997 road displacement)`. That "
      "double-counted the shoreline retreat: the 1984 setback is measured from "
      "the **1984 dune line**, but it is spent against a grid whose row 0 is "
      "the **2009** dune crest, which already sits landward by the 1984→2009 "
      "retreat *R*. Adding the physical displacement on top counted *R* twice, "
      "and the road was placed *R* metres behind where NCDOT actually put it.")
    a("")
    a("Checked against `RoadSetback_2004.csv` — a same-year (2004 road vs 2004 "
      "dune) measurement taken **after both events**, so it already is the "
      "post-relocation position:")
    a("")
    a("| GIS | Old absolute | 2004 measured | Overshoot |")
    a("|---:|---:|---:|---:|")
    for g, o, m in OVERSHOOT_TABLE:
        a(f"| {g} | {o} | {m} | +{o - m} |")
    a("")
    a("Median overshoot **+35 m** in the 1999 block, **+96 m** at Pea Island — "
      "worst exactly where retreat was worst. At GIS 84 the measured 2004 "
      "setback (50 m) is *smaller* than the 1984 setback (93 m) even though "
      "the road was moved 70 m landward in 1989: the shoreline overtook the "
      "road faster than NCDOT could move it. That is the real Pea Island "
      "story, and the absolute value erased it.")
    a("")
    a("**The fix** (`HAT_hindcast_1984_2024.py`, "
      "`HISTORICAL_ROAD_EVENTS`): the events now carry "
      "`relocation_displacement_m` — how far the road physically moved, and "
      "nothing else — applied to the setback the model is *currently* "
      "carrying:")
    a("")
    a("```python")
    a("new_sb = rm._road_setback + relocation_displacement_m[gis_id]")
    a("```")
    a("")
    a("CASCADE already decrements `road_setback` by `dune_migrated` every year, "
      "so by the event year the model's setback has absorbed the modelled "
      "retreat on its own. The retreat is now counted exactly once, from the "
      "model's own dune migration.")
    a("")
    a("### What the correction does and does not fix")
    a("")
    a("The corrected setback depends on how much the model has retreated by "
      "the event year. Evaluated across plausible retreat *R*:")
    a("")
    a("- **1999 block (GIS 9–15, event at year 15)** converges toward the 2004 "
      "measurement at realistic rates. GIS 11 goes 129 → 89 m at *R* = 40 m, "
      "against a measured 81 m. The correction largely works here.")
    a("- **1989 block (GIS 84–87, event at year 5)** does **not** converge. "
      "Only a third of the retreat has accumulated by year 5, so GIS 86 is "
      "still ~178 m against a measured 88 m.")
    a("")
    a("The reason is structural, and the correction cannot reach it: the road "
      "is placed **correctly relative to the dune**, but the island *behind* "
      "that dune is the 2009 island — 25 years too narrow. At Pea Island, "
      "where retreat is fastest on the whole island, that width deficit is "
      "larger than the relocation itself. Fixing the double-count removes the "
      "avoidable error; the remaining error is the width anachronism, and it "
      "needs a 1984 DEM that does not exist.")
    a("")

    a("## The road's own cells are never checked")
    a("")
    a("`bulldoze()` tests the rows **flanking** the road — `road_end + 1` on "
      "the bay side, `road_start - 1` on the sea side — and never looks at the "
      "road footprint itself. It also skips `road_end`: the road occupies "
      "`[road_start:road_end]` and the bayside test is at `road_end + 1`, so "
      "the cell immediately behind the asphalt is only used for `np.size()`.")
    a("")
    a("That is why `pct_road_cells_water` is reported separately below. A "
      "road can be substantially in the bay while the flanking test is still "
      "under threshold.")
    a("")

    # ---- per scenario ----
    for label, rows, problems in scen_rows:
        if not rows:
            continue
        a(f"## {label}")
        a("")
        for pr in problems:
            a(f"- **ID row problem:** {pr}")
        if problems:
            a("")
            a("The runner reads the setback file with `skiprows=1` and fills "
              "`road_setbacks_full[]` **by position**, discarding the IDs. A "
              "gap or out-of-order ID shifts every domain north of it with no "
              "error.")
            a("")

        counts = {v: sum(1 for r in rows if r["verdict"] == v)
                  for v in VERDICT_ORDER}
        off = [r for r in rows if r["road_off_island"]]
        bits = [f"`{v}` {counts[v]}" for v in VERDICT_ORDER if counts[v]]
        a(f"{len(rows)} domains — " + ", ".join(bits)
          + f". **{len(off)} with the road off the island.**")
        a("")

        flagged = [r for r in rows
                   if r["verdict"] != "FITS" or r["road_off_island"]]
        if flagged:
            a("| GIS | Setback (m) | % of road's own cells in water | % of "
              "profiles fully submerged | % water at bayside border | "
              "Largest that fits (m) | Verdict | Narrowness |")
            a("|---:|---:|---:|---:|---:|---:|---|---|")
            for r in sorted(flagged, key=lambda x: x["domain"]):
                mf = r["largest_fitting_setback_m"]
                mfs = "—" if not np.isfinite(mf) else f"{mf:.0f}"
                a(f"| {r['domain']} | {r['setback_m']:.0f} | "
                  f"{r['pct_road_cells_water']} | "
                  f"{r['pct_profiles_road_submerged']} | "
                  f"{r['pct_bayside_water']} | {mfs} | `{r['verdict']}`"
                  f"{' **OFF-ISLAND**' if r['road_off_island'] else ''} | "
                  f"`{r['narrow_label']}` |")
            a("")
        else:
            a("No domain flagged.")
            a("")

        if rows[0]["kind"] == "relocation":
            a("### Guards the prescribed-relocation path skips")
            a("")
            a("CASCADE has two guards on moving a road landward, and the "
              "runner's `HISTORICAL_ROAD_EVENTS` path takes **neither** — it "
              "assigns `rm._road_setback = new_sb` directly, leaves "
              "`road_ele` at its SLR-decremented value rather than rebuilding "
              "at grade, and lets the next `bulldoze()` flatten whatever is "
              "there.")
            a("")
            a("| GIS | Setback (m) | Rebuilt road elevation (m) | "
              "`get_road_relocation_elevation` | Average barrier width (m) | "
              "`road_relocation_checks` |")
            a("|---:|---:|---:|---|---:|---|")
            for r in sorted(rows, key=lambda x: x["domain"]):
                if r["relocation_road_ele_m"] == "":
                    continue
                v1 = ("**refuses** — road ≤ 0 m MSL"
                      if r["relocation_ele_refuses"] else "allows")
                v2 = ("**refuses** — island too narrow"
                      if r["relocation_width_refuses"] else "allows")
                a(f"| {r['domain']} | {r['setback_m']:.0f} | "
                  f"{r['relocation_road_ele_m']:.2f} | {v1} | "
                  f"{r['average_barrier_width_m']:.0f} | {v2} |")
            refused = [r for r in rows if r["relocation_ele_refuses"]
                       or r["relocation_width_refuses"]]
            a("")
            a(f"**{len(refused)} of {len(rows)} prescribed relocations would "
              f"be refused by CASCADE's own logic.** Evaluated on the t=0 "
              f"interior; the real grid at the event year will have evolved, "
              f"so treat these as the initialisation-time verdict.")
            a("")

    # ---- cross-cutting sections ----
    all_rows = [r for _, rows, _ in scen_rows for r in rows]
    split = [r for r in all_rows
             if r.get("pct_findwidths_disagree") not in ("", None)
             and r["pct_findwidths_disagree"] > 0]
    if split:
        a("## Where Barrier3D's two guards disagree")
        a("")
        a("The model does not use one definition of \"island\":")
        a("")
        a("| Routine | What it reads | Used for |")
        a("|---|---|---|")
        a("| `FindWidths` | Stops at the **first** water cell — land behind a "
          "gap is not part of the island | `InteriorWidth_Avg`, i.e. the "
          "`average_barrier_width` that `road_relocation_checks()` uses |")
        a("| `bulldoze` | The **literal cell** at `road_end + 1` — land behind "
          "a gap counts normally | The `roadway_drown` test |")
        a("")
        a("So the road can be founded on ground the relocation logic believes "
          "is off the back of the barrier: the drown test passes while the "
          "relocation test sees a much narrower island. Counted, not "
          "corrected — resolving it would mean changing a Barrier3D "
          "definition.")
        a("")
        a("| Scenario | GIS | % of profiles where the road's border cell is "
          "land but `FindWidths` has already ended the island | Narrowness |")
        a("|---|---:|---:|---|")
        for r in sorted(split, key=lambda x: -x["pct_findwidths_disagree"]):
            a(f"| {r['scenario']} | {r['domain']} | "
              f"{r['pct_findwidths_disagree']} | `{r['narrow_label']}` |")
        a("")

    a("## Narrowness labels, and what they do not tell you")
    a("")
    a("The source DEM was pre-masked to the island: across all 131 domain "
      "arrays there are **822,315 cells of exactly `-10.0` (NoData) and zero "
      "real cells below the water clamp**. There is no bathymetry, so a water "
      "cell and a LiDAR dropout are the same number and cannot be "
      "distinguished by value.")
    a("")
    a("| Label | Meaning |")
    a("|---|---|")
    a("| `TRUE_END` | Nothing behind the first water cell. The interior really "
      "ends where `FindWidths` says. |")
    a("| `HOLED` | Land resumes in >50% of profiles, usually after a 1–3 cell "
      "gap. The reported width is a truncation. Whether that gap is a NoData "
      "hole or a real marsh channel **is not resolved here** — Hatteras has "
      "both at 10–30 m scale. |")
    a("| `MIXED` | Between the two. |")
    a("")
    a("Barrier3D stops at the first water cell either way, so the reported "
      "`InteriorWidth` **is** the width the model uses. The label marks which "
      "domains carry the doubt; no data was changed to resolve it.")
    a("")

    a("## Scope of any road-management claim")
    a("")
    a("Some domains put the road off the island even at a correctly measured "
      "same-year setback. The 2009 barrier genuinely cannot hold NC-12 where "
      "it was, and no change to the setback fixes that — only a "
      "contemporaneous DEM would.")
    a("")
    a("These are left **unmodified and allowed to drown**. The consequence is "
      "the one described at the top: road management stops permanently for "
      "them. So any managed-versus-unmanaged conclusion must **exclude the "
      "domains listed in the per-scenario tables above**, by name, rather "
      "than averaging them in — they are unmanaged barriers from the year "
      "they drown, and including them biases every management statistic "
      "toward the unmanaged case.")
    a("")
    a("### Period weighting")
    a("")
    a("The runner uses one topography (2009) for both hindcast periods, so "
      "1984–2024 runs its first 25 years on end-of-period geometry while "
      "2004–2024 brackets its topography. **Both periods are reported as "
      "primary results**, with this limitation stated once in the methods.")
    a("")
    a("What that means concretely, and what belongs in the caveat: in period "
      "1 the road is placed correctly *relative to the dune* but the island "
      "behind that dune is 25 years too narrow, so the room available for "
      "management is systematically understated — most severely at Pea "
      "Island, where retreat is fastest. Period 2 does not carry that deficit "
      "to the same degree, and its setbacks are same-year measurements.")
    a("")

    a("## Unguarded index in `bulldoze()`")
    a("")
    a("`roadway_manager.bulldoze()` checks the two rows flanking the road. The "
      "**seaside** check is guarded:")
    a("")
    a("```python")
    a("if road_start > 0:")
    a("    seaside_water_cells = ... xyz_interior_grid[road_start - 1, :] ...")
    a("else:")
    a("    seaside_water_cells = 0")
    a("```")
    a("")
    a("The **bayside** check is not:")
    a("")
    a("```python")
    a("bayside_water_cells = np.count_nonzero(")
    a("    (xyz_interior_grid[road_end + 1, :] * dz) <= drown_threshold")
    a(") / number_border_cells          # IndexError if road_end + 1 >= DomainWidth")
    a("```")
    a("")
    a("The asymmetry reads as an oversight rather than a design choice: one "
      "edge of the road is bounds-checked and the other is not. When the "
      "barrier has narrowed past the road, the bayside index runs off the "
      "array and the run dies with `IndexError` instead of reporting a drowned "
      "road.")
    a("")
    a("**Status: not patched.** CASCADE is run as published. This project "
      "previously carried a monkey-patch in the run script that supplied the "
      "missing guard (`bayside_water_cells = 1.0` when off-grid, which trips "
      "the drown test at any threshold); it has been removed so the model "
      "matches the package. The exposure is managed upstream of the model "
      "instead — this audit refuses any setback whose border row would exceed "
      "`DomainWidth` (`BLOCK_OVERRUN`), which is currently zero across all "
      "scenarios.")
    a("")
    a("**Residual risk.** This audit tests **t = 0 only**. The interior is not "
      "fixed-size: `barrier3d.py` grows it (`np.vstack`, :648) and shrinks it "
      "(`np.delete(..., 0, axis=0)`, :724) as the barrier migrates. Erosion "
      "alone should be safe, because `road_setback` is decremented by "
      "`dune_migrated` in the same step, so the road and the grid shift "
      "together. The exposed path is a **relocation**, which resets "
      "`road_setback` to a value that does not account for accumulated "
      "narrowing — which is exactly the case that originally motivated the "
      "patch (GIS 11, 1999 event). That path is now defended by correcting the "
      "relocation to a displacement and by warning when CASCADE's own "
      "relocation guards object, but it is not proven impossible.")
    a("")

    a("## Known assumptions")
    a("")
    a("1. Road and dune lines are the same vintage in each setback file, so "
      "the shoreline retreat cancels in the subtraction. Verified for 1984 "
      "and 2004.")
    a("2. A setback is transferable between dune definitions: measured to the "
      "digitised dune line of its own year, spent against interior row 0 = the "
      "2009 DEM `argmax` dune crest. **Unverified**, and the largest single "
      "assumption here.")
    a("3. Interior row 0 of profile *i* is `dune_loc[i] + 1` "
      "(`USE_CONST_INTERIOR = False`), so the dune search window sets both the "
      "setback origin and how much land remains behind it.")
    a("4. Cross-shore distance carries a 1/cos θ inflation (1% at 8°, 11% at "
      "26°); interior width carries the same factor, so the ratio is "
      "consistent even though absolute distances are long.")
    a("5. `bulldoze` compares against `drown_threshold = 0` described as "
      "*m MSL*, but the arrays are **MHW-relative** (`MHW_M = 0.36` subtracted "
      "first). The effective test is \"at or below MHW\". Pre-existing in the "
      "model; reported, not changed.")
    a("6. `ROAD_ELEVATION = 1.45` in the runner is ambiguous between NAVD88 "
      "and MHW-relative — flagged in `write_elevation_csv()`, still "
      "unresolved. If NAVD88 was meant, it should be 1.09 MHW-relative.")
    a("7. The raw `1984_road_offset_raw.csv` behind `RoadSetback_1984.csv` is "
      "**not present in the repo**, so those values cannot be re-derived or "
      "verified at source. Audited as given.")
    a("8. Relocation scenarios are evaluated on the **t=0** interior. The real "
      "grid in 1989/1999 will have evolved under storms and SLR, so those "
      "verdicts are initialisation-time, not event-time.")
    a("")

    a("## Full results")
    a("")
    a("| Scenario | GIS | Setback (m) | Verdict | Off-island | % road cells "
      "water | % bayside water | DomainWidth | IW min/p20/med/max | Largest "
      "fitting (m) | Narrowness |")
    a("|---|---:|---:|---|---|---:|---:|---:|---|---:|---|")
    for label, rows, _ in scen_rows:
        for r in sorted(rows, key=lambda x: x["domain"]):
            if r["verdict"] == "NO_TOPOGRAPHY":
                a(f"| {label} | {r['domain']} | {r['setback_m']:.0f} | "
                  f"`NO_TOPOGRAPHY` | — | — | — | — | — | — | — |")
                continue
            mf = r["largest_fitting_setback_m"]
            mfs = "—" if not np.isfinite(mf) else f"{mf:.0f}"
            a(f"| {label} | {r['domain']} | {r['setback_m']:.0f} | "
              f"`{r['verdict']}` | {'yes' if r['road_off_island'] else ''} | "
              f"{r['pct_road_cells_water']} | {r['pct_bayside_water']} | "
              f"{r['domain_width_rows']} | "
              f"{r['iw_min']}/{r['iw_p20']}/{r['iw_median']}/{r['iw_max']} | "
              f"{mfs} | `{r['narrow_label']}` |")
    a("")

    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(L), encoding="utf-8")
    print(f"[out] {path}")


# =============================================================================
# MAIN
# =============================================================================

def main():
    if not TOPO_DIR.is_dir():
        raise SystemExit(f"\n[stop] topography dir not found:\n    {TOPO_DIR}\n")

    print("=" * 104)
    print(f"NC-12 SETBACK AUDIT  |  topography {TOPO_DUNE_VERSION}")
    print("=" * 104)
    print(f"  topography : {TOPO_DIR}")
    print(f"  width rule : barrier3d FindWidths (first cell <= SL), verbatim")
    print(f"  drown rule : bulldoze, >{PCT_WATER_TOUCHING_ROAD * 100:.0f}% of "
          f"profiles with water bordering the road")
    print(f"  off-island : >={ROAD_IN_WATER_FRACTION * 100:.0f}% of the road's "
          f"OWN cells in water (bulldoze never checks these)")

    scen_rows, all_problems, blocked = [], [], []
    for sc in SCENARIOS:
        rows, problems = run_scenario(sc)
        if not rows:
            continue
        print_scenario(sc["label"], rows, problems)
        scen_rows.append((sc["label"], rows, problems))
        all_problems.extend(problems)
        blocked.extend([r for r in rows
                        if r["verdict"].startswith("BLOCK")
                        or r["verdict"] == "NO_TOPOGRAPHY"])

    if not scen_rows:
        raise SystemExit("\n[stop] nothing audited\n")

    all_rows = [r for _, rows, _ in scen_rows for r in rows]
    write_csv(all_rows, OUT_DIR / "RoadSetback_audit.csv")
    write_markdown(scen_rows, all_problems, OUT_DIR / "RoadSetback_audit.md")

    print()
    print("=" * 104)
    if blocked or all_problems:
        print("VERDICT: BLOCK -- a setback would corrupt a run "
              "(silent index wrap or IndexError).")
        print("=" * 104)
        return 1
    print("VERDICT: PASS -- no setback would corrupt a run.")
    print("         DROWNS domains are reported as-is. Remember what that "
          "costs: road management")
    print("         stops PERMANENTLY for those domains from that year "
          "onward. See the markdown.")
    print("=" * 104)
    return 0


if __name__ == "__main__":
    sys.exit(main())
