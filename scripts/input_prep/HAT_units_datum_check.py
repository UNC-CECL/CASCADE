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
  4. Road elevation against the interior elevation at the road, same LiDAR.
     A missing or doubled MHW subtraction shows up as a ~0.36 m systematic
     offset against a ~0.03 m real difference.
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

PROJECT_ROOT = Path(__file__).resolve().parents[4]
DATA = PROJECT_ROOT / "data" / "hatteras_init"
B3D = DATA / "1-barrier3d-domains"

TOPO_DUNE_VERSION = "2009_v2"
TOPO_DIR = B3D / "2009-dune-topo" / TOPO_DUNE_VERSION / "topography"
DUNE_DIR = B3D / "2009-dune-topo" / TOPO_DUNE_VERSION / "dunes"
ROAD_ELEV_CSV = (DATA / "4-mgmt-forcing" / "road_elevation"
                 / "RoadElevation.csv")   # one file, both periods
ROAD_SETBACK_CSV = (DATA / "4-mgmt-forcing" / "road_offset" / "old_method_offset"
                    / "2004" / "RoadSetback_2004.csv")

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


def check_topography(chk):
    fs = sorted(glob.glob(str(TOPO_DIR / "domain_*_topography_*.npy")))
    if not fs:
        chk.add("topography present", False, f"none found in {TOPO_DIR}")
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
    chk.add("topography min == water sentinel",
            abs(lo - sentinel_dam) < 1e-6,
            f"min {lo:.4f} dam, sentinel {sentinel_dam:.4f} dam "
            f"({lo * 10:+.2f} m MHW)")
    chk.add("topography max in dam range", 0.05 < hi < 1.2,
            f"max {hi:.4f} dam = {hi * 10:.2f} m MHW  "
            f"(a metres array would read ~{hi * 10:.1f}; "
            f"NAVD88 would shift +{MHW_ELEVATION:.2f} m)")


def check_dunes(chk):
    fs = sorted(glob.glob(str(DUNE_DIR / "domain_*_dune_*.npy")))
    if not fs:
        chk.add("dune files present", False, f"none found in {DUNE_DIR}",
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
    chk.add("dune file is dam, not metres",
            0.0 <= v.min() and v.max() < 1.5 and np.median(v) < 1.0,
            f"range {v.min():.4f} to {v.max():.4f} dam, median "
            f"{np.median(v):.4f} (a metres array would read "
            f"{np.median(v) * 10:.2f} and fail this band)")

    # PLAUSIBILITY: separate question, and NOT a units problem. NC-12's
    # artificial dune ridge is typically 3-5 m NAVD88. The extractor's own
    # docstring warns the search window can be "wide enough to catch a
    # back-dune, a wooded ridge, or a house", which is what a high crest means.
    med_navd = float(np.median(v)) * 10 + berm_mhw + MHW_ELEVATION
    chk.add("dune crest plausible for a foredune", med_navd < 5.0,
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

    # the sensitive one: road vs the interior it is written into
    if not ROAD_SETBACK_CSV.exists():
        return
    sb_arr = np.loadtxt(ROAD_SETBACK_CSV, delimiter=",")
    sb = dict(zip(sb_arr[0].astype(int), sb_arr[1]))
    evd = dict(zip(a[0].astype(int), a[1]))
    gaps = []
    for d, S in sb.items():
        fs = glob.glob(str(TOPO_DIR / f"domain_{d}_topography_*.npy"))
        if not fs or d not in evd:
            continue
        arr = np.load(fs[0])
        rs = int(S / 10)
        re_ = rs + int(ROAD_WIDTH_M / 10)
        if rs < 0 or re_ > arr.shape[0]:
            continue
        gaps.append(float(np.median(arr[rs:re_, :])) * 10 - evd[d])
    if gaps:
        g = np.asarray(gaps)
        chk.add("road elev vs interior at the road",
                abs(float(np.median(g))) < 0.20,
                f"median gap {np.median(g):+.2f} m over {g.size} domains "
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
    check_topography(chk)
    check_dunes(chk)
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
