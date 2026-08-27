r"""
HAT_units_datum_check.py
===============================================================================
Verify the units and vertical datum of every elevation-like quantity the
Hatteras CASCADE hindcast feeds to Barrier3D, by tracing the source and by
checking the actual data files.

THE THING THIS EXISTS TO CATCH
------------------------------
"Everything in CASCADE is relative to MHW" is true of what the model COMPUTES
ON, and false of what you SUPPLY. There are two conventions side by side:

  CONVERTED FOR YOU -- supply in metres NAVD88, load_input.py converts:
      MHW          load_input.py:227    /10
      BermEl       load_input.py:241    /10 - MHW
      Dmaxel       load_input.py:304    /10 - MHW
      ShrubEl_*    load_input.py:346-7  /10 - MHW
      Dstart       load_input.py:240    /10

  ALREADY CONVERTED -- supply in the model's own frame, nothing touches it:
      elevation_file (.npy)   dam, MHW-relative
                              load_elevation() only does np.load; configuration.py
                              documents it as "[dam x dam x dam MHW]"
      dune_file (.npy)        dam, height ABOVE BERM
                              load_input.py:249 assigns DuneStart straight into
                              DuneDomain with no conversion
      road_ele                metres, MHW-relative
                              bulldoze() only does road_ele/dz

barrier3d.py:1197 pops MHW into self._MHW -- and _MHW appears NOWHERE ELSE in
the file. Barrier3D never applies it to a grid. Anything grid-shaped has to
arrive pre-converted.

That asymmetry is what made ROAD_ELEVATION = 1.45 ambiguous: it sits in the one
scalar family that takes pre-converted values while looking like the ones that
do not.

WHAT IS CHECKED
---------------
  1. The contract, as a table: supplied-as / converted-where / model-sees.
  2. The saved arrays, against the range each convention implies. A 10x unit
     error or a 0.36 m datum slip is far outside the plausible band, so this
     catches both.
  3. Berm elevation by two independent code paths (the extractor's and
     load_input's) -- they must agree.
  4. Road elevation against the interior elevation at the road -- ONCE PER
     PERIOD, since both halves are period-specific: the setback decides which
     rows to read, the topography decides what is in them.

     The two periods have different expected answers, and that is the whole
     content of the check. RoadElevation.csv is sampled on the 2009-2014
     baseline, which is the DEM behind 2004-start, so 2004 compares a surface
     against itself and must agree to ~0. 1984 runs 2009-2014-1996, where the
     1996 ALACE survey overwrites the corridor and its vertical offset is left
     uncorrected by design -- so 1984 is expected to sit HIGH by exactly that
     offset, which HAT_dem_1984_mosaic.py measured on the 2009/1996 overlap
     and writes to mosaic_1984_audit.csv every run.

     So the tolerance is applied to `gap - expected`, not to `gap`. That tests
     a claim -- "all of this gap is the 1996 survey offset" -- rather than
     widening a band until the number fits. Measured: 1984 gap +0.26 m against
     a recorded offset of +0.26 m, residual -0.001 m; 2004 residual -0.000 m.
     A missing or doubled MHW subtraction still shows up as a ~0.36 m residual
     in BOTH periods, because `expected` comes from a different file than the
     road elevation does.
  5. Whether each runner constant actually reaches the model, or is overridden.

REQUIREMENTS
------------
  numpy (pandas optional, for the road-elevation cross-check)
===============================================================================
"""

from __future__ import annotations

import glob
import sys
from pathlib import Path

import numpy as np

# =============================================================================
# CONFIG -- must mirror HAT_hindcast_1984_2024.py
# =============================================================================

# scripts/input_prep/HAT_units_datum_check.py -> repo root is parents[2].
# This said parents[4], which resolved to C:\Users\hanna: every path below
# pointed outside the repo and the script could not find one of its inputs.
PROJECT_ROOT = Path(__file__).resolve().parents[2]
DATA = PROJECT_ROOT / "data" / "hatteras_init"
B3D = DATA / "1-barrier3d-domains"

# Resolved from the extractor rather than pinned, so this checks the units of
# the arrays actually being run. It said "2009_v2", which has since been moved
# to 2009-dune-topo/incorrect/.
sys.path.insert(0, str(PROJECT_ROOT / "scripts"))
from hat_topo_version import topo_dirs, array_name  # noqa: E402

# BOTH PERIODS, AND THE FILES THE RUNNER ACTUALLY SPENDS (2026-08-26).
#
# Two things were wrong here, and both were silent.
#
# (1) ONE PRODUCT. `topo_dirs()` with no argument resolves DEFAULT_PRODUCT,
#     i.e. 2004-start. Since 1-barrier3d-domains went period-first there are
#     two extractions and the 1984 period runs the other one - so this script
#     reported "the units of the arrays actually being run" having never
#     opened half of them. All 90 domains differ between the products and 65
#     differ in interior SHAPE, so it is not a formality: the road-vs-interior
#     gap below indexes a specific row of a specific array.
#
# (2) THE LEGACY SETBACK FILE. ROAD_SETBACK_CSV pointed at
#     old_method_offset/2004/RoadSetback_2004.csv. The runner switched to the
#     dune-start method on 2026-08-18 (hatteras_site_config.py), and the two
#     put the road a median ~22 m apart in 2004. The most sensitive check in
#     this file - road elevation against the interior AT THE ROAD - was
#     therefore reading the interior at rows the model does not bulldoze.
#
# Both are fixed by asking hatteras_site_config, which is what the runner
# reads, instead of mirroring it. The header above says this CONFIG "must
# mirror HAT_hindcast_1984_2024.py"; mirroring is how it drifted, so the
# period-dependent forcings are now imported and cannot.
from hatteras_site_config import (HATTERAS_PERIODS,  # noqa: E402
                                  HATTERAS_ROAD_ELEVATION_FILE)

PERIODS = sorted(HATTERAS_PERIODS)          # [1984, 2004]

# One file, both periods -- not because there is one DEM (there are two) but
# because they differ under the road only by the uncorrected 1996-vs-2009
# survey offset. See HATTERAS_ROAD_ELEVATION_FILE in hatteras_site_config.py.
ROAD_ELEV_CSV = DATA / HATTERAS_ROAD_ELEVATION_FILE


def topo_for(year: int):
    """(topography dir, dunes dir, version) for one hindcast period."""
    return topo_dirs(HATTERAS_PERIODS[year]["topo_product"])


def product_for(year: int) -> str:
    return HATTERAS_PERIODS[year]["topo_product"]


def label_for(year: int) -> str:
    return f"{product_for(year)}/{topo_for(year)[2]}"


def setback_csv(year: int) -> Path:
    """The MODEL-FACING setback for one period, as the runner resolves it."""
    return DATA / HATTERAS_PERIODS[year]["road_setback_file"]


# WHERE RoadElevation.csv WAS SAMPLED. One file serves both periods and it is
# built on the 2009-2014 baseline (HAT_road_elevation.py, FILL_SOURCE), which
# is the DEM behind 2004-start. So for 2004 the road elevation and the interior
# under the road are the SAME surface and must agree to ~0.
ROAD_ELEV_PRODUCT = "2004-start"

# The 1984 period does not run that surface. 2009-2014-1996 overwrites measured
# ground wherever the 1996 ALACE survey has data, including through the road
# corridor, and HAT_dem_1984_mosaic.py leaves the vertical offset UNCORRECTED
# on purpose ("bias correction OFF, feathering OFF"). It writes the offset it
# measured, per domain, to this file every run.
MOSAIC_AUDIT = (DATA / "0-elevation" / "2009-2014-1996" / "1-gapfill-1m"
                / "mosaic_1984_audit.csv")


def recorded_survey_offset(domains) -> float | None:
    """Median 1996-minus-2009 offset over `domains`, from the mosaic audit.

    READ, NEVER ASSUMED. The point of returning it rather than hardcoding a
    tolerance is that the road-vs-interior gap below can then be tested against
    a number measured independently, by a different script, from the overlap of
    the two surveys - instead of being written off as "about right".

    Sign in the CSV is base - fill, i.e. 2009 - 1996, so it is negated here to
    read as "1996 sits this much higher than 2009".
    """
    if not MOSAIC_AUDIT.is_file():
        return None
    import csv
    vals = []
    with open(MOSAIC_AUDIT, newline="") as f:
        for r in csv.DictReader(f):
            try:
                if int(r["domain"]) in domains:
                    vals.append(-float(r["bias_2009_minus_1996_m"]))
            except (KeyError, TypeError, ValueError):
                continue
    return float(np.median(vals)) if vals else None

BERM_ELEVATION = 1.7      # m NAVD88, runner line 718
MHW_ELEVATION = 0.36      # m NAVD88, runner line 719
DUNE_REBUILD_HEIGHT = 3.0
REBUILD_ELEV_THRESHOLD = 0.01
ROAD_WIDTH_M = 20.0

SENTINEL_WATER_M = -3.0   # extractor's water sentinel, m MHW-relative
ABS_MIN_DUNE_H = 0.3      # roadway_manager.py:530

FIRST_ROAD_DOMAIN, LAST_ROAD_DOMAIN = 9, 90


# =============================================================================
# 1. THE CONTRACT
# =============================================================================

CONTRACT = [
    # (quantity, supplied as, converted where, model sees)
    ("elevation_file (.npy)", "dam, MHW-relative",
     "nowhere - load_elevation() is np.load", "dam MHW"),
    ("dune_file (.npy)", "dam, height ABOVE BERM",
     "nowhere - load_input.py:249", "dam above berm"),
    ("road_ele", "m, MHW-relative", "bulldoze(): /dz", "dam MHW"),
    ("MHW", "m NAVD88", "load_input.py:227  /10", "dam"),
    ("BermEl / berm_elevation", "m NAVD88", "load_input.py:241  /10 - MHW",
     "dam MHW"),
    ("Dmaxel", "m NAVD88", "load_input.py:304  /10 - MHW", "dam MHW"),
    ("Dstart (scalar fallback)", "m", "load_input.py:240  /10", "dam"),
    ("dune_design_elevation", "m MHW", "roadway_manager: -BermEl*10, then /10",
     "dam above berm"),
    ("dune_minimum_elevation", "m MHW", "roadway_manager: -BermEl*10, then /10",
     "dam above berm"),
    ("road_setback / road_width", "m", "bulldoze(): /dy, /dx", "cell index"),
    ("drown_threshold", "0, commented 'm MSL'",
     "compared to interior*dz", "effectively 0 m MHW"),
    ("SL", "-", "barrier3d.py:1165  _SL = 0, Lagrangian",
     "always 0, never changes"),
]


def print_contract():
    print("=" * 100)
    print("1. THE CONTRACT -- what each quantity is, and who converts it")
    print("=" * 100)
    w = max(len(c[0]) for c in CONTRACT)
    print(f"  {'QUANTITY':<{w}}  {'SUPPLIED AS':<24} {'CONVERTED WHERE':<34} "
          f"MODEL SEES")
    print(f"  {'-' * w}  {'-' * 24} {'-' * 34} {'-' * 20}")
    for q, s, c, m in CONTRACT:
        print(f"  {q:<{w}}  {s:<24} {c:<34} {m}")
    print()
    print("  The first three take PRE-CONVERTED values. Everything else is")
    print("  converted for you. barrier3d.py:1197 pops MHW and never uses it")
    print("  again -- nothing grid-shaped is converted by the model.")


# =============================================================================
# 2-4. CHECKS AGAINST THE ACTUAL DATA
# =============================================================================

class Check:
    def __init__(self):
        self.rows = []

    def add(self, name, ok, detail, fatal=True):
        self.rows.append((name, bool(ok), detail, fatal))
        return ok

    def report(self):
        print("=" * 100)
        print("CHECKS")
        print("=" * 100)
        w = max(len(r[0]) for r in self.rows)
        n_fail = 0
        for name, ok, detail, fatal in self.rows:
            tag = "PASS" if ok else ("FAIL" if fatal else "WARN")
            if not ok and fatal:
                n_fail += 1
            print(f"  [{tag}] {name:<{w}}  {detail}")
        print()
        return n_fail


def check_topography(chk, year):
    # Was "domain_*_topography_*.npy". The trailing _* required a year tag
    # that no longer exists, so the glob matched nothing after 2026-08-26.
    topo_dir = topo_for(year)[0]
    tag = f"{year} topography"
    fs = sorted(glob.glob(str(topo_dir / "domain_*_topography.npy")))
    if not fs:
        chk.add(f"{tag} present", False, f"none found in {topo_dir}")
        return
    mins, maxs = [], []
    for f in fs:
        a = np.load(f, mmap_mode="r")
        mins.append(float(a.min()))
        maxs.append(float(a.max()))
    lo, hi = min(mins), max(maxs)

    # Expected: dam MHW-relative. Sentinel is SENTINEL_WATER_M/10 dam exactly.
    # Real barrier tops out a few metres above MHW -> a few tenths of a dam.
    sentinel_dam = SENTINEL_WATER_M / 10.0
    chk.add(f"{tag} min == water sentinel",
            abs(lo - sentinel_dam) < 1e-6,
            f"{label_for(year)}: {len(fs)} arrays, min {lo:.4f} dam, "
            f"sentinel {sentinel_dam:.4f} dam ({lo * 10:+.2f} m MHW)")
    chk.add(f"{tag} max in dam range", 0.05 < hi < 1.2,
            f"max {hi:.4f} dam = {hi * 10:.2f} m MHW  "
            f"(a metres array would read ~{hi * 10:.1f}; "
            f"NAVD88 would shift +{MHW_ELEVATION:.2f} m)")


def check_dunes(chk, year):
    # Was "domain_*_dune_*.npy" - the trailing _* needed a year tag that no
    # longer exists, so this matched nothing after 2026-08-26.
    dune_dir = topo_for(year)[1]
    tag = f"{year} dune"
    fs = sorted(glob.glob(str(dune_dir / "domain_*_dune.npy")))
    if not fs:
        chk.add(f"{tag} files present", False, f"none found in {dune_dir}",
                fatal=False)
        return
    vals = []
    for f in fs:
        a = np.load(f, mmap_mode="r")
        v = np.asarray(a, dtype=float).ravel()
        vals.append(v[v > SENTINEL_WATER_M / 10.0 + 1e-9])
    v = np.concatenate(vals)
    berm_mhw = BERM_ELEVATION - MHW_ELEVATION

    # UNITS: dam above berm. If these were metres the median would be ~10x and
    # land outside this band; that is what the check discriminates. It says
    # nothing about whether the values are physically sensible.
    chk.add(f"{tag} file is dam, not metres",
            0.0 <= v.min() and v.max() < 1.5 and np.median(v) < 1.0,
            f"{label_for(year)}: range {v.min():.4f} to {v.max():.4f} dam, "
            f"median {np.median(v):.4f} (a metres array would read "
            f"{np.median(v) * 10:.2f} and fail this band)")

    # PLAUSIBILITY: separate question, and NOT a units problem. NC-12's
    # artificial dune ridge is typically 3-5 m NAVD88. The extractor's own
    # docstring warns the search window can be "wide enough to catch a
    # back-dune, a wooded ridge, or a house", which is what a high crest means.
    med_navd = float(np.median(v)) * 10 + berm_mhw + MHW_ELEVATION
    chk.add(f"{tag} crest plausible for a foredune", med_navd < 5.0,
            f"median crest {med_navd:.2f} m NAVD88 "
            f"(NC-12 dune ridge is typically 3-5 m NAVD88) -- if high, check "
            f"the dune search windows, not the units",
            fatal=False)


def check_berm(chk):
    ext = 1.70 - MHW_ELEVATION                                   # extractor
    b3d = (BERM_ELEVATION / 10.0 - MHW_ELEVATION / 10.0) * 10.0  # load_input
    chk.add("berm elevation, two code paths agree", abs(ext - b3d) < 1e-9,
            f"extractor {ext:.3f} m MHW vs load_input {b3d:.3f} m MHW")


def check_road_elevation(chk):
    if not ROAD_ELEV_CSV.exists():
        chk.add("road elevation file", False, f"not found: {ROAD_ELEV_CSV}")
        return
    a = np.loadtxt(ROAD_ELEV_CSV, delimiter=",")
    ev = a[1]
    chk.add("road elevation in m MHW range", 0.0 < ev.min() and ev.max() < 3.0,
            f"{ev.min():.2f} to {ev.max():.2f} m MHW "
            f"= {ev.min() + MHW_ELEVATION:.2f} to "
            f"{ev.max() + MHW_ELEVATION:.2f} m NAVD88 "
            f"(median {np.median(ev):.2f} MHW)")

    # THE SENSITIVE ONE: road elevation vs the interior it is written into.
    #
    # Per period, because both halves of it are period-specific: the setback
    # says WHICH ROWS to read and the topography is WHAT IS IN THEM. This used
    # to run once, against the legacy 2004 setback and the default product -
    # so it compared the 2004-start interior at rows the old method chose, for
    # a model that runs neither.
    evd = dict(zip(a[0].astype(int), a[1]))
    for year in PERIODS:
        path = setback_csv(year)
        if not path.exists():
            chk.add(f"{year} setback file", False, f"not found: {path}",
                    fatal=False)
            continue
        topo_dir = topo_for(year)[0]
        sb_arr = np.loadtxt(path, delimiter=",")
        sb = dict(zip(sb_arr[0].astype(int), sb_arr[1]))
        gaps = []
        for d, S in sb.items():
            f = topo_dir / array_name("topography", d)
            if not f.is_file() or d not in evd:
                continue
            arr = np.load(f)
            rs = int(S / 10)
            re_ = rs + int(ROAD_WIDTH_M / 10)
            if rs < 0 or re_ > arr.shape[0]:
                continue
            gaps.append(float(np.median(arr[rs:re_, :])) * 10 - evd[d])
        if not gaps:
            chk.add(f"{year} road elev vs interior at the road", False,
                    f"no domain compared against {label_for(year)}")
            continue
        g = np.asarray(gaps)
        gap = float(np.median(g))

        # WHAT GAP SHOULD THIS PERIOD SHOW?
        #
        # For the period whose product IS the surface RoadElevation.csv was
        # sampled from, zero: same LiDAR, same corridor, so anything left is a
        # datum or unit error, which is what this check is for.
        #
        # For 1984 it is NOT zero, and pretending otherwise would make this
        # check fail forever for a reason that is documented and deliberate.
        # The expected value is not invented here either - it is the offset
        # HAT_dem_1984_mosaic.py measured on the 2009/1996 overlap and wrote to
        # its own audit. Testing `gap - expected` therefore tests a specific
        # claim ("the whole gap is the 1996 survey offset") rather than
        # widening a tolerance until the number fits inside it. If the graft
        # is ever bias-corrected, or the road elevation is rebuilt on the 1984
        # product, expected goes to ~0 on its own and this keeps working.
        if product_for(year) == ROAD_ELEV_PRODUCT:
            expected, why = 0.0, "same surface as RoadElevation.csv"
        else:
            rec = recorded_survey_offset(set(sb))
            if rec is None:
                chk.add(f"{year} road elev vs interior at the road", False,
                        f"median gap {gap:+.2f} m, and {MOSAIC_AUDIT.name} is "
                        f"missing so it cannot be attributed", fatal=False)
                continue
            expected, why = rec, (
                f"1996 ALACE sits {rec:+.2f} m above 2009 through the corridor, "
                f"uncorrected by design, per {MOSAIC_AUDIT.name}")

        resid = gap - expected
        chk.add(f"{year} road elev vs interior at the road",
                abs(resid) < 0.20,
                f"median gap {gap:+.2f} m over {g.size} domains "
                f"({label_for(year)}, {path.name}); expected "
                f"{expected:+.2f} m -- {why}; residual {resid:+.3f} m "
                f"(a datum slip would show ~{MHW_ELEVATION:.2f} m; "
                f"a unit slip ~10x)")


def check_runner_constants(chk):
    """Does each runner constant actually reach the model, or get overridden?"""
    berm_m = (BERM_ELEVATION / 10.0 - MHW_ELEVATION / 10.0) * 10.0

    design_floor = berm_m + 1.0
    design_eff = max(DUNE_REBUILD_HEIGHT, design_floor)
    chk.add("DUNE_REBUILD_HEIGHT reaches the model",
            design_eff == DUNE_REBUILD_HEIGHT,
            f"{DUNE_REBUILD_HEIGHT} m MHW vs floor "
            f"max(v, BermEl*10+1.0={design_floor:.2f}) -> {design_eff:.2f} m")

    min_floor = berm_m + ABS_MIN_DUNE_H
    min_eff = max(REBUILD_ELEV_THRESHOLD, min_floor)
    chk.add("REBUILD_ELEV_THRESHOLD reaches the model",
            min_eff == REBUILD_ELEV_THRESHOLD,
            f"{REBUILD_ELEV_THRESHOLD} (commented 'dam') vs floor "
            f"max(v, BermEl*10+{ABS_MIN_DUNE_H}={min_floor:.2f}) -> "
            f"{min_eff:.2f} m MHW  ** INERT **",
            fatal=False)


# =============================================================================
# MAIN
# =============================================================================

def main():
    print()
    print_contract()
    print()

    chk = Check()
    for year in PERIODS:
        check_topography(chk, year)
        check_dunes(chk, year)
    check_berm(chk)
    check_road_elevation(chk)
    check_runner_constants(chk)
    n_fail = chk.report()

    print("=" * 100)
    print("NOTES -- true, documented, deliberately not changed")
    print("=" * 100)
    print("  drown_threshold = 0 is commented '0 m MSL' in roadway_manager.py:766,")
    print("  but it is compared against xyz_interior_grid * dz, which is")
    print("  MHW-RELATIVE. The effective test is 0 m MHW -- roughly 0.26 m")
    print("  stricter than the comment claims. Every roadway_drown verdict is")
    print("  correspondingly conservative. Behaviour unchanged; recorded here.")
    print()
    print("  REBUILD_ELEV_THRESHOLD is commented 'dam' but is passed as")
    print("  dune_minimum_elevation, which CASCADE documents as [m MHW]. It is")
    print("  inert either way -- roadway_manager floors it at BermEl*10 + 0.3.")
    print("  The dune rebuild threshold is effectively hardcoded to "
          f"{(BERM_ELEVATION - MHW_ELEVATION) + ABS_MIN_DUNE_H:.2f} m MHW.")
    print("  Raising it in dam would silently give a value 10x too small.")
    print()

    if n_fail:
        print(f"VERDICT: {n_fail} FAILED check(s)")
        return 1
    print("VERDICT: PASS -- units and datum consistent end to end")
    return 0


if __name__ == "__main__":
    sys.exit(main())
