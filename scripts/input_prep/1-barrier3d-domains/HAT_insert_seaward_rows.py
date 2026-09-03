r"""
HAT_insert_seaward_rows.py
==============================================================================
Move the island seaward, per domain, by prepending fabricated interior rows.

WHY THIS EXISTS
---------------
The 1984-start topography is a 1996 ALACE beach and foredune grafted onto a
2009 backdune. The NC-12 road lines are 1984-vintage. Between 1984 and 1996 the
island migrated landward, so at the worst domains the 1984 roadbed now sits
SEAWARD of the 1996 dune crest and the measured setback comes out negative:

    GIS 85   road seaward cell 13, interior row 0 = source cell 14
             setback (13 - 14) * 10 = -10 m  ->  floored to 0 in the model CSV

A zero setback is not inert. roadway_manager.road_relocation_checks does
`road_setback += dune_migrated` every year, so the first year with any erosion
drives it negative and relocates the road. In the shipped 1984-2004 run GIS 85
relocates in YEAR 1 and twice more after that; GIS 86 the same; the relocation
year across GIS 10/11/84/85/86 is monotone in the input setback, not in the
physics. That is the artefact this script exists to remove.

WHAT IT DOES
------------
For each selected domain it prepends N rows to the SEAWARD end of the saved
interior array, which moves interior row 0 N cells seaward and turns the road's
setback from (road - row0) into (road - row0) + N*10.

N comes from a MEASUREMENT, not a target: the per-domain cross-shore offset
between the digitized 1984 dune line and interior row 0, measured in the
extractor's own frame by HAT_measure_duneline_shift.py. Island-wide that median
is +18.9 m; at GIS 85 it is +65.9 m with a p10-p90 of +62 to +70.

The control that makes this credible is in the 2004 pair: measured the same way,
the 2004 dune line and 2004-start row 0 agree to +1.4 / -4.1 / -3.6 m on the
stable mid-island (GIS 40/50/60). So "digitized line vs DEM row 0" carries no
material feature offset, and the 1984 number is a date difference rather than a
definitional one.

WHAT IS FABRICATED, AND SAY SO
------------------------------
No survey covers land that was gone by 1996. The N new rows are invented. The
fill rule is explicit and recorded in the manifest; `backdune` (the default)
lays them flat at the median of interior rows 1-3, i.e. a backdune platform, NOT
at row 0's elevation -- at GIS 85 row 0 IS the 4.82 m crest, and copying it
would build a 70 m plateau at crest height.

Note what stays behind: the old row 0 becomes an interior ridge N cells inside
the island. That is a relict foredune, which is a reasonable thing for a
migrating barrier to have, but it is a consequence rather than a choice.

TWO VARIANTS, BECAUSE THE BAY SIDE IS NOT FREE
----------------------------------------------
brie_coupler.offset_shoreline sets x_s per domain and the interior extends
LANDWARD from the dune, so prepending rows does not push the shoreline seaward
-- it pushes the bay edge further into the sound.

    pad         prepend N. Island gets N cells wider. Mean interior height and
                InteriorWidth_AvgTS both rise, which feeds overwash flux and the
                relocation room test.
    translate   prepend N and retire the N landward-most LAND rows per column to
                the water sentinel. Per-column land width is preserved. This is
                the barrier-migration reading: in 1984 the bay edge was also N
                cells seaward.

Neither is free and they are not equivalent; run both and compare before
picking. `none` writes the setback CSV alone and leaves the topography be.

OUTPUT
------
A NEW dune-topo version beside the source. v1 is never written to.

    1984-start/dune-topo/<DST_VERSION>/
        topography/   domain_<N>_topography.npy    (dam)
        dunes/        domain_<N>_dune.npy          copied unchanged
        RUN_MANIFEST.txt
        HAT_seaward_row_insert_audit.csv
        RoadSetback_1984_dunestart.csv             matched to this topography

USAGE
-----
    python HAT_insert_seaward_rows.py                       # measured N, pad
    python HAT_insert_seaward_rows.py --variant translate
    python HAT_insert_seaward_rows.py --n-rule minimum --domains 85,86
==============================================================================
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


def _find_root(start: Path) -> Path:
    for p in [start, *start.parents]:
        if (p / "data" / "hatteras_init").is_dir():
            return p
    raise SystemExit("could not locate the project root (no data/hatteras_init above me)")


REPO = _find_root(Path(__file__).resolve())
sys.path.insert(0, str(REPO / "scripts"))

from hat_topo_version import (array_name, dune_topo_root,  # noqa: E402
                              duneline_shift_dir, topo_dirs)

OFFSET_SCRIPT = (REPO / "scripts" / "input_prep" / "4-mgmt-forcings" / "road_offset"
                 / "1-produce" / "HAT_road_offset_from_dune_start.py")

SRC_PRODUCT = "1984-start"
YEAR = 1984

# Where the measured shift comes from. Written by HAT_measure_duneline_shift.py.
SHIFT_DIR = duneline_shift_dir(SRC_PRODUCT)
# TWO measurements of the same quantity, and they disagree by a factor of ~3
# at the domains this work is about (GIS 85: 65.9 m vs 19.5 m). Neither is
# the truth; --shift-source names which one a build used and the manifest
# records it, so no output can be read without knowing.
SHIFT_SOURCES = {
    # THE ONE TO USE (2026-09-02). 1984 line minus 1997 line, same feature at
    # both ends, so the DEFINITIONAL offset between a digitized line and interior
    # row 0 cancels and what is left is pure date. Measured, not assumed: the
    # 1997 line sits +16.2 m seaward of row 0 island-wide (IQR +12.8 to +21.0),
    # and subtracting that leaves an island-wide date term of +0.8 m -- i.e. the
    # ORIGINAL "duneline" source below was about 85% feature offset.
    "date": SHIFT_DIR / "duneline_retreat_1984_1997.csv",
    # Superseded. Kept so the earlier arms stay reproducible and so the size of
    # the correction stays visible rather than being quietly absorbed.
    "duneline": SHIFT_DIR / "duneline_shift_1984.csv",
    # "dsas" REMOVED 2026-09-03. It pointed at a file that had already been
    # moved into superseded/, so the option failed at runtime; and the
    # estimate measures the SHORELINE, not the dune line, understating dune
    # retreat by ~2.3x (GIS 85: 19.5 m against a measured 58.2 m). A route to
    # a known-wrong number is not worth keeping reachable.
}

SETBACK_DIR = (REPO / "data" / "hatteras_init" / "4-mgmt-forcing" / "road_offset"
               / "dunestart_offset" / str(YEAR))

SENTINEL_DAM = -0.30          # SENTINEL_WATER_M / CELL_SIZE_M, the extractor's water
TOPO_ROWS = 200               # the extractor's cap; a padded array must still fit
BACKDUNE_ROWS = 3             # rows averaged for the `backdune` fill

# LAND is > 0 m MHW, not "> the water sentinel". At these domains the sound side
# is MEASURED marsh and shallow bay sitting between -3 m and 0 m, so a sentinel
# test calls the bay land and reports GIS 85 as 160 rows wide when its island is
# 37. 0 m MHW is the threshold bulldoze() itself uses (drown_threshold = 0) and
# the one HAT_road_domain_views draws.
LAND_DAM = 0.0


def load_offset_module():
    spec = _iu.spec_from_file_location("hat_off", OFFSET_SCRIPT)
    mod = _iu.module_from_spec(spec)
    sys.modules["hat_off"] = mod
    spec.loader.exec_module(mod)
    return mod


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


def read_setback_csv(path: Path) -> dict:
    rows = list(csv.reader(open(path)))
    ids = [int(float(v)) for v in rows[0] if v.strip() != ""]
    vals = [float(v) for v in rows[1] if v.strip() != ""]
    return dict(zip(ids, vals))


def write_setback_csv(path: Path, values: dict) -> None:
    ids = sorted(values)
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["{:.3f}".format(i) for i in ids])
        w.writerow(["{:.3f}".format(values[i]) for i in ids])


def real_block(ext, dom, row0: np.ndarray, n: int, n_along: int) -> np.ndarray:
    """The DEM cells that ALREADY EXIST seaward of interior row 0, as (n, along).

    build_interior fills topo[:, i] from prof_arr[i, row0[i]:], so the cells the
    insert is about to cover are prof_arr[i, row0[i]-n : row0[i]] -- per profile,
    because row 0 is a per-profile cut, not a horizontal one.

    These are real measurements, just excluded from the interior for being
    seaward of the dune pick. At GIS 85 they include cell 13 at 3.44 m, which is
    the road's own seaward cell -- the cell whose elevation decides what
    bulldoze() is scraping. Fabricating over it when the DEM has it measured is
    a loss for nothing.

    Returned in dam, NaN where the source cell is off the array. Cells at or
    below water are left as NaN too: they are 1996 beach, and in 1984 that
    ground was dry island. Those are the ones that genuinely have to be invented.
    """
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


def fabricate_rows(topo: np.ndarray, n: int, rule: str) -> np.ndarray:
    """The N invented rows, (n, n_along) in dam. topo is (rows, along), dam."""
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
    else:
        raise SystemExit("unknown fill rule {!r}".format(rule))
    # never fabricate land below the water sentinel
    return np.maximum(block, SENTINEL_DAM)


def lower_old_crest(new: np.ndarray, base: np.ndarray, n: int,
                    max_reach: int = 15):
    """Shave the DEM's dune ridge down to the backdune platform, per column.

    THE PROBLEM THIS SOLVES. Prepending N rows puts Barrier3D's dune at the 1984
    dune line, but the DEM's own crest is still standing in the interior N cells
    landward -- so the model starts with TWO dunes. That is the same sand counted
    twice: the ridge is at the later position *because* the dune migrated there
    by 1996, so in 1984 it had not formed yet. At GIS 85 it is 4.82 m, the
    tallest thing in the domain, sitting on the road's own cell and shielding it
    from landward.

    THE COST, STATED PLAINLY. This discards a real measurement. The defence is
    that it is a measurement of the wrong YEAR: the whole operation is de-aging
    the surface by 12-25 years, and a 1996 crest is not a 1984 initial condition
    just because it is real. The opposite choice -- keeping it -- is equally
    defensible and is what --no-lower-old-crest gives you. Build both.

    Walks landward from the seaward edge capping at the platform, and stops at
    the first cell already at or below it, so it shaves the ridge and nothing
    else. `max_reach` bounds the walk past the insert: a column whose profile
    never drops back to the platform would otherwise be flattened across the
    whole island, and that is a failure worth hearing about rather than
    absorbing.
    """
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


def retire_landward_rows(topo: np.ndarray, n: int) -> np.ndarray:
    """Drown the N landward-most LAND cells of each column, into the local bay.

    Per column, not per row: the island's landward edge is not a straight line,
    which is the whole reason the road drown test looks at flanking rows rather
    than a single width.

    The retired cells take the elevation of the bay immediately landward of them
    in the SAME column, not the -3 m sentinel. Stamping the sentinel would dig a
    30 m trench along the sound edge of a domain whose real back-barrier is
    measured marsh a few decimetres below MHW, and Barrier3D would read that as
    the island having calved rather than migrated.

    USES BARRIER3D'S OWN WIDTH DEFINITION, not a count of dry cells.

    THIS WAS A BUG AND IT MADE `translate` BEHAVE LIKE `pad` (fixed 2026-09-02).
    The first version counted every cell above 0 m MHW in a column and drowned
    the landward-most n of them. Barrier3D's FindWidths (barrier3d.py:29) does
    something else: it walks from row 0 and STOPS AT THE FIRST cell <= SL, so
    anything beyond an interior water gap is not island at all.

    On GIS 85 the two disagree in all 50 columns -- median 37.5 cells against 44,
    and in the worst column 25 against 52, because the profile dips to -0.02 m at
    row 26 and everything past it is sound-side marsh at 0.03-0.19 m. So the
    cells being drowned were out in that marsh, which the model was never
    counting; the retirement removed nothing while the prepended rows still
    added. Measured effect: t=0 island width rose 263->309 m (D84), 361->402
    (D85), 286->301 (D86) when it should not have moved at all.
    """
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


def land_stats(topo: np.ndarray):
    """(median land rows per column, mean elevation of land cells in m MHW)."""
    land = topo > LAND_DAM
    width = land.sum(axis=0)
    mean_m = float(np.mean(topo[land]) * 10.0) if land.any() else float("nan")
    return float(np.median(width)), mean_m


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--variant", choices=("pad", "translate", "none"), default="pad")
    ap.add_argument("--fill",
                    choices=("measured", "median", "backdune", "row0", "taper"),
                    default="measured",
                    help="measured: keep the real DEM cell wherever it is dry "
                         "land AND at or above the backdune platform, floor it "
                         "there otherwise. median: keep every dry cell as "
                         "measured and fill only the cells at or below MHW, "
                         "with the median of the block's own dry cells - one "
                         "guard, one constant, and no measurement is ever "
                         "raised.")
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
        # Scope by the MEASUREMENT, not by where the road happens to be. Any
        # domain whose 1984 dune line sits >= half a cell seaward of interior
        # row 0 is missing 1984 land from the 1996 survey, and that is true of
        # domains with no NC-12 in them as well. Selecting the two relocation
        # blocks made "unchanged" mean two different things along the island:
        # 30 domains had N >= 1 and were passed over by a scope decision rather
        # than by a measurement. This removes that discontinuity.
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
        # `none` still credits the setback with N: it is the variant that says
        # "the dune-to-road distance was wrong, the island was not", so the
        # correction lands entirely in the CSV and the arrays are passed through.
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
                # FLOOR THE REAL VALUE AT THE BACKDUNE PLATFORM, do not simply
                # take it. The measured cells here are the 1996 surface, and at
                # an eroding domain 1996 is the LATER and LOWER one -- so its
                # elevation is a lower bound on 1984's, not an estimate of it.
                # Taking it raw put a 0.57 m beach cell two rows inside the 1984
                # interior at GIS 85, which is 1996 beach standing where 1984 had
                # dry backdune. The floor keeps the real dune flank (1.98, 3.44)
                # and the road's own cell, and declines to import the beach.
                n_real = int((keep & (real >= block)).sum())
                block = np.where(keep, np.maximum(real, block), block)
            elif args.fill == "median":
                # KEEP EVERY DRY MEASUREMENT, INVENT ONE NUMBER.
                #
                # `measured` above justifies itself with "1996 is a lower bound
                # on 1984" and then RAISES 44% of the block above the DEM to
                # the platform, which is not a lower-bound operation - it
                # asserts the ground was at least backdune height. This rule
                # drops that second step. A dry cell is kept as measured; only
                # the cells with no usable measurement get a value, and that
                # value is the median of the block's OWN dry cells rather than
                # a statistic imported from interior rows 1-3.
                #
                # One guard, one constant, and it never overrides a measurement
                # upward. At GIS 85 it takes the measured share from 47% to
                # 91%: only 28 of 300 cells are at or below MHW, 25 of them in
                # the seaward-most row. NC-12 sits at rows 4-5, which are 100%
                # dry, so the constant does not reach the road at all.
                #
                # WHAT IT ADMITS, deliberately: the 1996 beach ramp. Per-row
                # medians at GIS 85 run -0.00, 0.70, 1.20, 1.84, 3.17, 4.96, so
                # row 1 sits BELOW whatever fills row 0 and the block carries a
                # dip two cells inside the island. That is what the platform
                # floor existed to remove. It is admitted here because it is
                # what the measurement says, and the alternative asserts more
                # than the data supports.
                if keep.any():
                    fill_v = float(np.median(real[keep]))
                    block = np.where(keep, real, fill_v)
                    n_real = int(keep.sum())
                # No dry cell anywhere in the block: nothing to take a median
                # of, so the backdune platform stands as the fallback and
                # n_real stays 0. Does not occur at any of the 38 domains, but
                # a silent nan-median would be worse than a stated fallback.
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

        # Audit EVERY domain that was processed, not only the ones carrying a
        # road. Under --domains measured the scope is the whole island, and a
        # domain with no NC-12 still has rows inserted and still has to be
        # accountable for them.
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
