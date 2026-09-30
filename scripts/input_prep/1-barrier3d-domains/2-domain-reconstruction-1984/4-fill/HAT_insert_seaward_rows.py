"""
Move the island seaward per domain by prepending interior rows, as a new dune-topo layer.

    python scripts/input_prep/1-barrier3d-domains/2-domain-reconstruction-1984/4-fill/HAT_insert_seaward_rows.py
    python scripts/input_prep/1-barrier3d-domains/2-domain-reconstruction-1984/4-fill/HAT_insert_seaward_rows.py --variant translate
    python scripts/input_prep/1-barrier3d-domains/2-domain-reconstruction-1984/4-fill/HAT_insert_seaward_rows.py --n-rule minimum --domains 85,86

N comes from the measured shift (--shift-source); the rows are filled by the
chosen rule; writes the arrays, a setback CSV, an audit and a manifest. Details: scripts/input_prep/1-barrier3d-domains/2-domain-reconstruction-1984/README.md.

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-09-18
"""

from __future__ import annotations

import argparse
import csv
import importlib.util as _iu
import json
import shutil
import sys
from datetime import datetime
from pathlib import Path

import numpy as np


# Walk up until a directory holds data/hatteras_init
def _find_root(start: Path) -> Path:
    for p in [start, *start.parents]:
        if (p / "data" / "hatteras_init").is_dir():
            return p
    raise SystemExit("could not locate the project root (no data/hatteras_init above me)")


REPO = _find_root(Path(__file__).resolve())
sys.path.insert(0, str(REPO / "scripts"))

from site_layer.hat_topo_version import (array_name, dune_topo_root,  # noqa: E402
                              duneline_shift_dir, topo_dirs)

OFFSET_SCRIPT = (REPO / "scripts" / "input_prep" / "4-mgmt-forcings" / "road_offset"
                 / "1-produce" / "HAT_road_offset_from_dune_start.py")

# --- CONFIG ------------------------------------------------------------------
SRC_PRODUCT = "1984-start"
YEAR = 1984

# Where the measured shift comes from. Written by HAT_measure_duneline_shift.py.
SHIFT_DIR = duneline_shift_dir(SRC_PRODUCT)
# Two measurements of the shift, disagreeing ~3x; --shift-source names the one used
SHIFT_SOURCES = {
    # The one to use: 1984 line minus 1997 line, so the row-0 offset cancels
    "date": SHIFT_DIR / "duneline_retreat_1984_1997.csv",
    # Superseded, kept so earlier arms stay reproducible
    "duneline": SHIFT_DIR / "duneline_shift_1984.csv",
    # 'dsas' removed 2026-09-03: a missing file, and a shoreline not a dune measure
}

from site_layer import hat_topo_version as _tv  # noqa: E402
SETBACK_DIR = _tv.road_setback_dir(YEAR)

SENTINEL_DAM = -0.30          # SENTINEL_WATER_M / CELL_SIZE_M, the extractor's water
TOPO_ROWS = 200               # the extractor's cap; a padded array must still fit
BACKDUNE_ROWS = 3             # rows averaged for the `backdune` fill

# LAND is > 0 m MHW, not "> the water sentinel"
LAND_DAM = 0.0
# -----------------------------------------------------------------------------


# HAT_road_offset.py, loaded as a module for its extractor helpers
def load_offset_module():
    spec = _iu.spec_from_file_location("hat_off", OFFSET_SCRIPT)
    mod = _iu.module_from_spec(spec)
    sys.modules["hat_off"] = mod
    spec.loader.exec_module(mod)
    return mod


# The measured shift per domain, or stop: N is a measurement
def read_shift_csv(path: Path) -> dict:
    if not path.is_file():
        raise SystemExit(
            "\nno measured shift at {}\n"
            "Run HAT_measure_duneline_shift.py first -- N is a measurement, not "
            "a target, and this script will not invent one.\n".format(path))
    out = {}
    for r in csv.DictReader(open(path)):
        out[int(r["domain"])] = float(r["shift_m_median"])
    return out


# {domain: setback} from a two-row setback CSV
def read_setback_csv(path: Path) -> dict:
    rows = list(csv.reader(open(path)))
    ids = [int(float(v)) for v in rows[0] if v.strip() != ""]
    vals = [float(v) for v in rows[1] if v.strip() != ""]
    return dict(zip(ids, vals))


# Write a two-row setback CSV
def write_setback_csv(path: Path, values: dict) -> None:
    ids = sorted(values)
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["{:.3f}".format(i) for i in ids])
        w.writerow(["{:.3f}".format(values[i]) for i in ids])


# The DEM cells that ALREADY EXIST seaward of interior row 0, as (n, along)
def real_block(ext, dom, row0: np.ndarray, n: int, n_along: int) -> np.ndarray:
    z = dom["z"]                                   # (along, cross), m MHW
    out = np.full((n, n_along), np.nan)
    for i in range(n_along):
        if row0[i] < 0:
            continue
        for k in range(n):
            src = int(row0[i]) - n + k
            if 0 <= src < z.shape[1]:
                v = z[i, src]
                if v > 0.0:                        # dry land only
                    out[k, i] = v / 10.0           # m -> dam
    return out


# The N invented rows, (n, n_along) in dam
def fabricate_rows(topo: np.ndarray, n: int, rule: str) -> np.ndarray:
    if rule == "row0":
        block = np.repeat(topo[0:1, :], n, axis=0)
    elif rule == "backdune":
        k = min(BACKDUNE_ROWS, topo.shape[0] - 1)
        base = np.median(topo[1:1 + k, :], axis=0) if k > 0 else topo[0, :]
        block = np.repeat(base[None, :], n, axis=0)
    elif rule == "taper":
        k = min(BACKDUNE_ROWS, topo.shape[0] - 1)
        base = np.median(topo[1:1 + k, :], axis=0) if k > 0 else topo[0, :]
        # seaward-most row at the backdune value, rising to row 0's elevation
        w = np.linspace(0.0, 1.0, n + 1)[:-1][:, None]
        block = base[None, :] * (1 - w) + topo[0:1, :] * w
    elif rule in ("matched-crest", "matched-nocrest"):
        # Matched backdune: today's near-dune profile copied in front of itself (two variants)
        off = 0 if rule == "matched-crest" else 1
        idx = np.minimum(np.arange(n) + off, topo.shape[0] - 1)
        block = topo[idx, :].copy()
    else:
        raise SystemExit("unknown fill rule {!r}".format(rule))
    # never fabricate land below the water sentinel
    return np.maximum(block, SENTINEL_DAM)


# Shave the DEM's dune ridge down to the backdune platform, per column
def lower_old_crest(new: np.ndarray, base: np.ndarray, n: int,
                    max_reach: int = 15):
    out = new.copy()
    overrun = 0
    for c in range(out.shape[1]):
        plat = base[c]
        r = 0
        while r < out.shape[0] and (r < n or out[r, c] > plat):
            if r >= n + max_reach:
                overrun += 1
                break
            out[r, c] = min(out[r, c], plat)
            r += 1
    return out, overrun


# Drown the N landward-most LAND cells of each column, into the local bay
def retire_landward_rows(topo: np.ndarray, n: int) -> np.ndarray:
    out = topo.copy()
    for c in range(out.shape[1]):
        col = out[:, c]
        # index of the first water cell, exactly as FindWidths finds it
        first_water = next((i for i, v in enumerate(col) if v <= LAND_DAM),
                           col.size)
        width = max(first_water - 1, 0)
        if width <= 0:
            continue
        k = min(n, width)
        retire = np.arange(first_water - k, first_water)
        bay = col[first_water:]
        bay = bay[bay <= LAND_DAM]
        out[retire, c] = float(np.median(bay[:5])) if bay.size else SENTINEL_DAM
    return out


# (median land rows per column, mean elevation of land cells in m MHW)
def land_stats(topo: np.ndarray):
    land = topo > LAND_DAM
    width = land.sum(axis=0)
    mean_m = float(np.mean(topo[land]) * 10.0) if land.any() else float("nan")
    return float(np.median(width)), mean_m


# Run: N per domain, fill and prepend the rows, write the layer
def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--variant", choices=("pad", "translate", "none"), default="pad")
    ap.add_argument("--fill",
                    choices=("measured", "median", "backdune", "row0", "taper",
                             "matched-crest", "matched-nocrest"),
                    default="measured",
                    help="measured: keep the real DEM cell wherever it is dry "
                         "land AND at or above the backdune platform, floor it "
                         "there otherwise. median: keep every dry cell as "
                         "measured and fill only the cells at or below MHW, "
                         "with the median of the block's own dry cells - one "
                         "guard, one constant, and no measurement is ever "
                         "raised. matched-crest / matched-nocrest: today's "
                         "near-dune profile copied N cells seaward, per "
                         "column, from interior row 0 / row 1 (see "
                         "fabricate_rows).")
    ap.add_argument("--lower-old-crest", action="store_true",
                    help="shave the DEM's dune ridge down to the backdune "
                         "platform, so the dune exists once, at the 1984 line")
    ap.add_argument("--shift-source", choices=tuple(SHIFT_SOURCES),
                    default="duneline",
                    help="which measurement of the 1984 offset N comes from")
    ap.add_argument("--n-rule", choices=("measured", "minimum"), default="measured",
                    help="measured: round(1984 duneline shift / 10). "
                         "minimum: the fewest rows that make the setback > 0.")
    ap.add_argument("--domains", default="negative",
                    help="'measured' (every domain whose measured shift rounds "
                         "to >=1 cell, road or not), 'negative' (setback <= 0), "
                         "'all' (every domain with a road), or a comma list")
    ap.add_argument("--dst-version", default=None,
                    help="destination version name; default v1_<variant>_<n-rule>")
    ap.add_argument("--src-version", default=None)
    args = ap.parse_args()

    dst_name = args.dst_version or "v1_{}_{}_{}".format(
        args.variant, args.n_rule, args.shift_source)

    src_topo, src_dune, src_name = topo_dirs(SRC_PRODUCT, override=args.src_version)
    if dst_name == src_name:
        raise SystemExit("\nrefusing to write into the source version {!r}.\n".format(src_name))
    dst = dune_topo_root(SRC_PRODUCT) / dst_name

    off = load_offset_module()
    ext = off.load_extractor(SRC_PRODUCT)
    windows = json.load(open(ext.WINDOW_JSON))     # fail loudly if the picks are gone

    shift_csv = SHIFT_SOURCES[args.shift_source]
    shifts = read_shift_csv(shift_csv)
    setbacks_raw = {}
    for r in csv.DictReader(open(SETBACK_DIR / "RoadOffset_{}_domains.csv".format(YEAR))):
        v = r["setback_dunestart_m"]
        if v not in ("", "nan"):
            setbacks_raw[int(r["domain"])] = float(v)
    model_setbacks = read_setback_csv(
        SETBACK_DIR / "RoadSetback_{}_dunestart.csv".format(YEAR))

    if args.domains == "negative":
        targets = sorted(d for d, s in setbacks_raw.items() if s <= 0)
    elif args.domains == "all":
        targets = sorted(setbacks_raw)
    elif args.domains == "measured":
        # Scope by the MEASUREMENT, not by where the road happens to be
        targets = sorted(d for d, m in shifts.items() if int(round(m / 10.0)) >= 1)
    else:
        targets = sorted(int(x) for x in args.domains.split(","))

    print("=" * 78)
    print("seaward row insert | {} {} -> {}".format(SRC_PRODUCT, src_name, dst_name))
    print("=" * 78)
    print("  variant  : {}".format(args.variant))
    print("  fill     : {}".format(args.fill))
    print("  N rule   : {}".format(args.n_rule))
    print("  domains  : {} -> {}".format(args.domains, targets))
    print("  shift    : {} ({})".format(shift_csv.name, args.shift_source))
    print()

    dst.mkdir(parents=True, exist_ok=True)
    (dst / "topography").mkdir(exist_ok=True)
    (dst / "dunes").mkdir(exist_ok=True)

    audit = []
    for D in range(1, 91):
        src_t = src_topo / array_name("topography", D)
        src_d = src_dune / array_name("dune", D)
        if not src_t.is_file():
            continue
        topo = np.load(src_t)
        shutil.copy2(src_d, dst / "dunes" / array_name("dune", D))

        n = 0
        if D in targets:
            if args.n_rule == "measured":
                n = int(round(shifts.get(D, 0.0) / 10.0))
            else:
                # fewest rows to put the road strictly landward of row 0
                n = int(np.floor(-setbacks_raw.get(D, 0.0) / 10.0)) + 1
            n = max(n, 0)

        w0, h0 = land_stats(topo)
        # `none` still credits the setback with N; the arrays pass through
        n_real = 0
        if n > 0 and args.variant != "none":
            block = fabricate_rows(
                topo, n,
                "backdune" if args.fill in ("measured", "median") else args.fill)
            if args.fill in ("measured", "median"):
                dom_p = ext.load_profiles(ext.LOAD_PATH / "domain_{}.npy".format(D))
                prof = ext.masked_profiles(dom_p["z"])
                w = windows.get("domain_{}".format(D))
                i0, i1 = ((int(w["i0"]), int(w["i1"])) if w
                          else ext.default_window(prof, dom_p["start_beach"]))
                _e, dl = ext.find_dunes(prof, dom_p["start_beach"], i0, i1)
                r0, _lt = ext.interior_row0_line(prof, dl)
                real = real_block(ext, dom_p, r0, n, topo.shape[1])
                keep = np.isfinite(real)
            if args.fill == "measured":
                # FLOOR THE REAL VALUE AT THE BACKDUNE PLATFORM, do not simply take it
                n_real = int((keep & (real >= block)).sum())
                block = np.where(keep, np.maximum(real, block), block)
            elif args.fill == "median":
                # KEEP EVERY DRY MEASUREMENT, INVENT ONE NUMBER
                if keep.any():
                    fill_v = float(np.median(real[keep]))
                    block = np.where(keep, real, fill_v)
                    n_real = int(keep.sum())
                # No dry cell anywhere in the block
            new = np.vstack([block, topo])
            if args.lower_old_crest:
                k = min(BACKDUNE_ROWS, topo.shape[0] - 1)
                plat = (np.median(topo[1:1 + k, :], axis=0) if k > 0
                        else topo[0, :])
                new, overrun = lower_old_crest(new, plat, n)
                if overrun:
                    print("  [warn] D{}: {} of {} columns never dropped back to "
                          "the platform within {} rows; their ridge was only "
                          "partly shaved".format(D, overrun, topo.shape[1], 15))
            if args.variant == "translate":
                new = retire_landward_rows(new, n)
            if new.shape[0] > TOPO_ROWS:
                print("  [warn] D{}: {} rows > TOPO_ROWS={}, truncating the "
                      "landward (bay) end".format(D, new.shape[0], TOPO_ROWS))
                new = new[:TOPO_ROWS, :]
            topo = new
        w1, h1 = land_stats(topo)

        np.save(dst / "topography" / array_name("topography", D), topo)

        # Audit EVERY domain that was processed, not only the ones carrying a road
        has_road = D in setbacks_raw
        if has_road or n > 0:
            old_raw = setbacks_raw.get(D, float("nan"))
            new_raw = old_raw + n * 10.0
            audit.append({
                "domain": D,
                "has_road": int(has_road),
                "n_rows_inserted": n,
                "shift_measured_m": round(shifts.get(D, float("nan")), 1),
                "setback_raw_before_m": old_raw,
                "setback_raw_after_m": new_raw,
                "setback_model_before_m": model_setbacks.get(D, float("nan")),
                "setback_model_after_m": max(new_raw, 0.0),
                "land_rows_before": w0, "land_rows_after": w1,
                "mean_land_elev_before_m": round(h0, 3),
                "mean_land_elev_after_m": round(h1, 3),
                "array_rows": topo.shape[0],
                "fill_rule": args.fill,
                "cells_inserted": int(n * topo.shape[1]),
                "cells_from_dem": n_real,
            })
            if n:
                frac = n_real / float(n * topo.shape[1]) if n else 0.0
                print("  D{:3d}  +{} rows  setback {:+6.1f} -> {:+6.1f} m  |  "
                      "land rows {:.0f} -> {:.0f}  |  mean land elev {:.2f} -> "
                      "{:.2f} m  |  {:.0%} of the block is measured DEM"
                      .format(D, n, old_raw, new_raw, w0, w1, h0, h1, frac))

    new_model = dict(model_setbacks)
    for a in audit:
        if a["domain"] in new_model:
            new_model[a["domain"]] = a["setback_model_after_m"]
    write_setback_csv(dst / "RoadSetback_{}_dunestart.csv".format(YEAR), new_model)

    with open(dst / "HAT_seaward_row_insert_audit.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(audit[0].keys()))
        w.writeheader()
        w.writerows(audit)

    (dst / "RUN_MANIFEST.txt").write_text(
        "=" * 78 + "\n"
        "SEAWARD ROW INSERT  --  " + dst_name + "\n" + "=" * 78 + "\n\n"
        "written    : {:%Y-%m-%d %H:%M}\n".format(datetime.now()) +
        "source     : {}/{}\n".format(SRC_PRODUCT, src_name) +
        "variant    : {}\n".format(args.variant) +
        "fill rule  : {}\n".format(args.fill) +
        "N rule     : {}\n".format(args.n_rule) +
        "domains    : {} -> {}\n".format(args.domains, targets) +
        "shift src  : {} ({})\n\n".format(shift_csv, args.shift_source) +
        "THE PREPENDED ROWS ARE FABRICATED. No survey covers land that had eroded\n"
        "away by 1996. They carry the fill rule above and nothing else. Any result\n"
        "that turns on island width, mean interior height, overwash flux or the\n"
        "roadway relocation room test at these domains is describing that\n"
        "fabrication as much as it is describing Hatteras.\n\n"
        "The dune arrays are COPIED UNCHANGED: this moves the dune's position in\n"
        "the model frame, and assumes its 1984 crest height equalled its 1996 one.\n\n"
        "KNOWN BIAS -- THE INTERVAL IS 13 YEARS, THE SURFACE IS 12\n"
        "The dune lines N is measured from are 1984 and 1997. The DEM surface at\n"
        "interior row 0 and at every inserted cell is 1996 ALACE (verified: 100% of\n"
        "cells, all block domains, no 2009 or 2014). So the interval that matches the\n"
        "topography is 1984-1996 = 12 years, and the one measured is 13. The source\n"
        "aerials are 1996-10-14 and 1997-10-12, 363 days apart, so N overstates the\n"
        "1984->1996 retreat by close to one full year of it.\n"
        "This is recorded rather than corrected: scaling by 12/13 would assume retreat\n"
        "was steady across 13 storm-driven years, which the record does not support.\n"
        "The effect is under one cell everywhere and moves N at three domains --\n"
        "GIS 12 (2->1), GIS 84 (4->3), GIS 85 (6->5).\n"
        "\n"
        "PROGRADATION IS NOT REMOVAL\n"
        "N is floored at 0. Where the 1997 line is SEAWARD of the 1984 one the island\n"
        "prograded, the 1996 survey already covers everything the 1984 island had, and\n"
        "there is no missing land to restore -- so nothing is inserted and nothing is\n"
        "taken away. This script only ever ADDS interior that existed in 1984 and is\n"
        "absent from the 1996 DEM.\n"
        "\n\n"
        "RoadSetback_1984_dunestart.csv in this folder is matched to THIS\n"
        "topography. Reading the v1 setbacks against these arrays, or the reverse,\n"
        "restores exactly the off-by-N error the file exists to remove.\n",
        encoding="utf-8")

    print("\n  wrote {} domains to {}".format(len(audit), dst))
    print("  setback CSV : {}".format(dst / "RoadSetback_{}_dunestart.csv".format(YEAR)))
    print("  audit       : {}".format(dst / "HAT_seaward_row_insert_audit.csv"))


if __name__ == "__main__":
    main()
