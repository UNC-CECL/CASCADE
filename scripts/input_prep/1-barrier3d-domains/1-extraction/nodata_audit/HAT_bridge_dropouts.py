"""
Fill the unsurveyed cells that three references agreed are survey dropouts, as a new dune-topo version.

    python scripts/input_prep/1-barrier3d-domains/1-extraction/nodata_audit/HAT_bridge_dropouts.py

Copies SRC_VERSION's arrays and picks to DST_VERSION, linearly bridges each
cleared hole along its profile, and writes the topography, nodata and bridged
masks plus BRIDGE_MANIFEST.txt. Details: scripts/input_prep/1-barrier3d-domains/1-extraction/nodata_audit/README.md.

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-09-18
"""

import csv
import datetime
import json
import shutil
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np

REPO = next(
    _p for _p in Path(__file__).resolve().parents
    if (_p / "pyproject.toml").exists())   # 1-extraction/nodata_audit/ since 2026-09-09
sys.path.insert(0, str(REPO / "scripts"))
from site_layer import hat_topo_version as htv  # noqa: E402
from cascade_pipeline import roadway  # noqa: E402

# --- CONFIG ------------------------------------------------------------------
TOPO_PRODUCT = "1984-start"
SRC_VERSION = "v1"
DST_VERSION = "v2"
DAM_TO_M = 10.0
SENTINEL_DAM = -0.3
EXTRACTOR = (REPO / "scripts" / "input_prep" / "1-barrier3d-domains" / "1-extraction"
             / "HAT_dune_topo_extractor.py")

AUDIT_SUBDIR = "nodata-audit"
# -----------------------------------------------------------------------------


# Give the new version its own pick set, seeded from the source version
def carry_picks_forward(src_version, dst_version):
    picks_dir = htv.picks_dir(TOPO_PRODUCT)
    src = picks_dir / f"HAT_dune_search_windows_{src_version}.json"
    dst = picks_dir / f"HAT_dune_search_windows_{dst_version}.json"

    if not src.is_file():
        print(f"  [picks] no {src.name} to carry forward -- SKIPPED. "
              f"{dst.name} will not exist, and anything resolving picks "
              f"through the extractor's VERSION will fail on it.")
        return None
    if dst.is_file():
        print(f"  [picks] {dst.name} already exists -- left alone")
        return dst, "exists"

    windows = json.loads(src.read_text(encoding="utf-8"))
    meta = windows.setdefault("_meta", {})
    meta["inherited_from"] = src_version
    meta["inherited_by"] = "nodata_audit/HAT_bridge_dropouts.py"
    meta["inherited_at"] = datetime.datetime.now().isoformat(timespec="seconds")
    meta["inherited_note"] = (
        f"NOT re-picked for {dst_version}. {dst_version} is {src_version} with "
        f"unsurveyed cells bridged inside existing arrays; interior shapes and "
        f"interior row 0 are unchanged, so the dune windows are unchanged too. "
        f"This copy exists so {dst_version} owns its pick set and a future "
        f"re-pick cannot overwrite {src_version}'s.")
    dst.write_text(json.dumps(windows, indent=2, sort_keys=True),
                   encoding="utf-8")
    n = sum(1 for k in windows if k != "_meta")
    print(f"  [picks] {src.name} -> {dst.name} ({n} domains, inherited)")
    return dst, "written"


# {(domain, profile)
def cleared_holes(src_run):
    a = src_run / AUDIT_SUBDIR
    verd = {(int(r["domain"]), int(r["profile"])): r["verdict"]
            for r in csv.DictReader((a / "hole_verdicts.csv").open())}
    cells = defaultdict(list)
    for r in csv.DictReader((a / "bracketed_hole_cells.csv").open()):
        k = (int(r["domain"]), int(r["profile"]))
        if verd.get(k) == "DROPOUT":
            cells[k].append(int(r["interior_row"]))
    return {k: sorted(v) for k, v in cells.items()}


# Linear fill of `rows` in one profile
def bridge_profile(col, rows):
    lo, hi = rows[0] - 1, rows[-1] + 1
    if lo < 0 or hi >= col.size:
        return False
    a, b = col[lo], col[hi]
    if a <= SENTINEL_DAM + 1e-9 or b <= SENTINEL_DAM + 1e-9:
        return False
    n = hi - lo
    for k, r in enumerate(rows, start=1):
        col[r] = a + (b - a) * k / n
    return True


# Run: carry the picks forward, bridge every cleared hole, write the version and its manifest
def main():
    src_topo, src_dune, src_ver = htv.topo_dirs(TOPO_PRODUCT, SRC_VERSION)
    src_run = src_topo.parent
    dst_run = src_run.parent / DST_VERSION
    dst_topo, dst_dune = dst_run / "topography", dst_run / "dunes"
    dst_topo.mkdir(parents=True, exist_ok=True)
    dst_dune.mkdir(parents=True, exist_ok=True)
    print(f"{TOPO_PRODUCT}: {src_ver} -> {DST_VERSION}")

    holes = cleared_holes(src_run)
    by_dom = defaultdict(list)
    for (d, p), rows in holes.items():
        by_dom[d].append((p, rows))
    print(f"  {len(holes)} holes cleared as DROPOUT in "
          f"{len(by_dom)} domains: {sorted(by_dom)}")

    n_cells = n_ok = n_skip = 0
    width_gain = {}
    for gis in range(1, 91):
        tname = htv.array_name("topography", gis)
        t = np.load(src_topo / tname)
        nod = np.load(src_topo / htv.array_name("nodata", gis))
        bridged = np.zeros_like(nod)
        before = roadway.interior_widths(t * DAM_TO_M)

        for p, rows in by_dom.get(gis, []):
            col = t[:, p]
            if bridge_profile(col, rows):
                bridged[rows, p] = True
                n_ok += 1
                n_cells += len(rows)
            else:
                n_skip += 1
                print(f"    [skip] D{gis} p{p}: not bracketed by measured "
                      f"land in the saved array")

        after = roadway.interior_widths(t * DAM_TO_M)
        if by_dom.get(gis):
            width_gain[gis] = float((after - before).mean() * DAM_TO_M)

        np.save(dst_topo / tname, t)
        np.save(dst_topo / htv.array_name("nodata", gis), nod)
        np.save(dst_topo / htv.array_name("bridged", gis), bridged)
        dname = htv.array_name("dune", gis)
        shutil.copy2(src_dune / dname, dst_dune / dname)

    print(f"  bridged {n_ok} holes, {n_cells} cells"
          + (f"   ({n_skip} skipped)" if n_skip else ""))
    for d, g in sorted(width_gain.items()):
        print(f"    D{d}: mean island width {g:+.0f} m")

    # Verification
    print("\n  verifying against the source:")
    bad = 0
    for gis in range(1, 91):
        a = np.load(src_topo / htv.array_name("topography", gis))
        b = np.load(dst_topo / htv.array_name("topography", gis))
        m = np.load(dst_topo / htv.array_name("bridged", gis))
        if a.shape != b.shape:
            print(f"    [FAIL] D{gis} shape {a.shape} -> {b.shape}"); bad += 1
        diff = ~np.isclose(a, b)
        if not np.array_equal(diff, m):
            print(f"    [FAIL] D{gis} changed cells != bridged mask "
                  f"({diff.sum()} vs {m.sum()})"); bad += 1
        if m.any() and not (b[m] > 0).all():
            print(f"    [FAIL] D{gis} a bridged cell is not above MHW"); bad += 1
    tot_bridged = sum(int(np.load(dst_topo / htv.array_name("bridged", g)).sum())
                      for g in range(1, 91))
    print(f"    shapes identical, {tot_bridged} cells changed, "
          f"all above MHW, {bad} failures")
    if bad:
        raise SystemExit("\nverification failed - v2 not activated\n")

    # Manifest
    (dst_run / "BRIDGE_MANIFEST.txt").write_text(f"""\
{TOPO_PRODUCT} dune-topo {DST_VERSION}
================================================================
NOT produced by HAT_dune_topo_extractor.py.

{DST_VERSION} is {src_ver} with {n_cells} unsurveyed cells filled by linear
interpolation between the measured cells bracketing them, in {n_ok} holes
across domains {sorted(by_dom)}.

Those holes were cleared as lidar DROPOUTS rather than ponds by three
independent references - the 2014 NCFMP hydro-flattening stamp, the shape of
the nodata blob, and a manual review of 1996 aerial imagery - in
nodata_audit/HAT_test_hole_pond_or_dropout.py. Everything judged POND or left
UNCLEAR still carries the -3.0 m water sentinel.

Interior shapes and interior row 0 are IDENTICAL to {src_ver}, so road setbacks
measured against {src_ver} remain valid.

MASKS
  domain_<N>_nodata.npy    unchanged from {src_ver}: no survey saw this cell
  domain_<N>_bridged.npy   NEW: and its value is now interpolated

  A bridged cell is still unsurveyed. Do not collapse the two masks.

RE-RUNNING THE EXTRACTOR WILL DESTROY THIS
  The extractor's VERSION literal is now "{DST_VERSION}", so a re-run writes a
  fresh UNBRIDGED extraction over these arrays. Recover with:
      python scripts/input_prep/1-barrier3d-domains/1-extraction/nodata_audit/\\
          HAT_bridge_dropouts.py

TO REVERT
  Set dune-topo/CURRENT back to "{src_ver}" AND the extractor's VERSION to
  "{src_ver}". Both, or the extractor literal silently wins.

mean island width gained: """ + ", ".join(
        f"D{d} {g:+.0f} m" for d, g in sorted(width_gain.items())) + "\n",
        encoding="utf-8")

    # Activate

    # Picks first, so the tree never resolves to a version with no pick set
    carry_picks_forward(src_ver, DST_VERSION)

    (dst_run.parent / "CURRENT").write_text(DST_VERSION, encoding="utf-8")
    src = EXTRACTOR.read_text(encoding="utf-8")
    for line in src.splitlines():
        if line.startswith("VERSION "):
            src = src.replace(line, f'VERSION = "{DST_VERSION}"'
                                    f'             # bumped by '
                                    f'nodata_audit/HAT_bridge_dropouts.py', 1)
            break
    EXTRACTOR.write_text(src, encoding="utf-8")

    got = htv.topo_dirs(TOPO_PRODUCT)[2]
    print(f"\n  CURRENT -> {DST_VERSION}, extractor VERSION -> {DST_VERSION}")
    print(f"  topo_dirs('{TOPO_PRODUCT}') now resolves to: {got}")
    if got != DST_VERSION:
        raise SystemExit(f"\nactivation failed: still resolving {got}\n")
    print(f"  wrote {dst_run}")


if __name__ == "__main__":
    main()
