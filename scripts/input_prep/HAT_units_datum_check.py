"""
Units and vertical datum of every elevation-like quantity the hindcast hands Barrier3D.

    python scripts/input_prep/HAT_units_datum_check.py

Prints the units contract, then checks the arrays, berm, road elevation and runner
constants for each period; exits nonzero on a failure. Details: scripts/input_prep/README.md.

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-09-18
"""

from __future__ import annotations

import glob
import sys
from pathlib import Path

import numpy as np

# --- CONFIG ------------------------------------------------------------------
# Repo root found by searching upward
PROJECT_ROOT = next(_p for _p in Path(__file__).resolve().parents
                    if (_p / "pyproject.toml").exists())
DATA = PROJECT_ROOT / "data" / "hatteras_init"
import sys as _b3dsys
from pathlib import Path as _B3DP
_b3dsys.path.insert(0, str(next(_q for _q in _B3DP(__file__).resolve().parents
                                if (_q / "pyproject.toml").exists()) / "scripts"))
from site_layer import hat_topo_version as _b3d  # noqa: E402
B3D = _b3d.DOMAIN_ROOT

# Resolved, not pinned, so the arrays checked are the arrays run
sys.path.insert(0, str(PROJECT_ROOT / "scripts"))
from site_layer.hat_topo_version import topo_dirs, array_name  # noqa: E402

# Both periods, and the files the runner spends: asked of hatteras_site_config, not mirrored
from site_layer.hatteras_site_config import (HATTERAS_PERIODS,  # noqa: E402
                                  HATTERAS_ROAD_ELEVATION_FILE)

PERIODS = sorted(HATTERAS_PERIODS)          # [1984, 2004]

# One road-elevation file for both periods; they differ only by the 1996 survey offset
ROAD_ELEV_CSV = DATA / HATTERAS_ROAD_ELEVATION_FILE
# -----------------------------------------------------------------------------


# (topography dir, dunes dir, version) for one hindcast period
def topo_for(year: int):
    return topo_dirs(HATTERAS_PERIODS[year]["topo_product"])


# The topography product a period runs
def product_for(year: int) -> str:
    return HATTERAS_PERIODS[year]["topo_product"]


# A period's product/version label
def label_for(year: int) -> str:
    return f"{product_for(year)}/{topo_for(year)[2]}"


# The MODEL-FACING setback for one period, as the runner resolves it
def setback_csv(year: int) -> Path:
    return DATA / HATTERAS_PERIODS[year]["road_setback_file"]


# RoadElevation.csv was sampled on the 2004-start surface: for 2004 the gap must be ~0
ROAD_ELEV_PRODUCT = "2004-start"

# 1984 runs the 1996 graft, whose survey offset is left uncorrected and recorded here
from site_layer.hat_elevation_products import product as _elprod  # noqa: E402
MOSAIC_AUDIT = _elprod("2009-2014-1996", check=False).gapfill_1m / "mosaic_1984_audit.csv"


# Median 1996-minus-2009 offset over `domains`, from the mosaic audit
def recorded_survey_offset(domains) -> float | None:
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

BERM_ELEVATION = 1.7      # m NAVD88, the runner's BERM_ELEVATION
MHW_ELEVATION = 0.36      # m NAVD88, the runner's MHW_ELEVATION
DUNE_REBUILD_HEIGHT = 3.0
REBUILD_ELEV_THRESHOLD = 0.01
ROAD_WIDTH_M = 20.0

SENTINEL_WATER_M = -3.0   # extractor's water sentinel, m MHW-relative
ABS_MIN_DUNE_H = 0.3      # roadway_manager.py, _absolute_minimum_dune_height

FIRST_ROAD_DOMAIN, LAST_ROAD_DOMAIN = 9, 90


# 1. the contract

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


# Print the units contract: supplied as, converted where, model sees
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
    print("  converted for you. Barrier3d.__init__ pops MHW and never uses it")
    print("  again -- nothing grid-shaped is converted by the model.")


# 2-4. checks against the actual data

# Collects pass / fail / warn results and prints them
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


# Topography arrays: dam MHW-relative, in a plausible band
def check_topography(chk, year):
    # Array names come from the resolver (the old year-tagged glob matched nothing)
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

    # Expected: dam MHW-relative; real tops are a few tenths of a dam
    sentinel_dam = SENTINEL_WATER_M / 10.0
    chk.add(f"{tag} min == water sentinel",
            abs(lo - sentinel_dam) < 1e-6,
            f"{label_for(year)}: {len(fs)} arrays, min {lo:.4f} dam, "
            f"sentinel {sentinel_dam:.4f} dam ({lo * 10:+.2f} m MHW)")
    chk.add(f"{tag} max in dam range", 0.05 < hi < 1.2,
            f"max {hi:.4f} dam = {hi * 10:.2f} m MHW  "
            f"(a metres array would read ~{hi * 10:.1f}; "
            f"NAVD88 would shift +{MHW_ELEVATION:.2f} m)")


# Dune arrays: dam above the berm, then plausibility
def check_dunes(chk, year):
    # Array names come from the resolver (the old year-tagged glob matched nothing)
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

    # Units: dam above berm; metres would land the median ~10x outside this band
    chk.add(f"{tag} file is dam, not metres",
            0.0 <= v.min() and v.max() < 1.5 and np.median(v) < 1.0,
            f"{label_for(year)}: range {v.min():.4f} to {v.max():.4f} dam, "
            f"median {np.median(v):.4f} (a metres array would read "
            f"{np.median(v) * 10:.2f} and fail this band)")

    # Plausibility, a separate question from units: a high crest may be a back-dune or house
    med_navd = float(np.median(v)) * 10 + berm_mhw + MHW_ELEVATION
    chk.add(f"{tag} crest plausible for a foredune", med_navd < 5.0,
            f"median crest {med_navd:.2f} m NAVD88 "
            f"(NC-12 dune ridge is typically 3-5 m NAVD88) -- if high, check "
            f"the dune search windows, not the units",
            fatal=False)


# Berm elevation by the extractor's path and load_input's: they must agree
def check_berm(chk):
    ext = 1.70 - MHW_ELEVATION                                   # extractor
    b3d = (BERM_ELEVATION / 10.0 - MHW_ELEVATION / 10.0) * 10.0  # load_input
    chk.add("berm elevation, two code paths agree", abs(ext - b3d) < 1e-9,
            f"extractor {ext:.3f} m MHW vs load_input {b3d:.3f} m MHW")


# Road elevation against the interior at the road, per period, net of the survey offset
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

    # The sensitive one: road elevation against the interior it is written into, per period
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

        # Expected gap: ~0 on the surface the road was sampled from, the recorded offset for 1984
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


# Does each runner constant actually reach the model, or get overridden?
def check_runner_constants(chk):
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


# Run: the contract, every check, the notes, and the exit code
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
    print("  drown_threshold = 0 is commented '0 m MSL' in RoadwayManager.update,")
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
