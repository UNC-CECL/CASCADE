"""
HAT_bridge_dropouts.py

Fill the unsurveyed cells that three references agreed are survey dropouts, and
write the result as a NEW dune-topo version.

WHAT IT CHANGES, AND WHAT IT REFUSES TO CHANGE
-----------------------------------------------
Only the cells that HAT_test_hole_pond_or_dropout.py cleared as DROPOUT. Each
such hole is bracketed by MEASURED land on both sides along its own profile, so
the fill is a straight linear interpolation between two measurements. No value
is invented beyond what those two measurements already imply, and nothing that
was judged POND or left UNCLEAR is touched - those keep the -3.0 m sentinel.

    interior row:   r-1      r    r+1   r+2    r+3
    v1 (dam):      0.111   -0.3  -0.3  -0.3   0.116     <- three unsurveyed
    v2 (dam):      0.111   0.112 0.113 0.114  0.116     <- interpolated

THE NODATA MASK IS DELIBERATELY LEFT ALONE
-------------------------------------------
A bridged cell is still a cell no survey ever saw. Its value is now an
interpolation, not a measurement, and any consumer asking "was this measured?"
must still get NO. So `<stem>_nodata.npy` is copied through unchanged and a
SECOND mask, `<stem>_bridged.npy`, records which cells were filled. Two masks,
two different questions:

    nodata   True = no survey saw this cell            (unchanged from v1)
    bridged  True = and its value is now interpolated  (new in v2)

Collapsing them into one would destroy the only record that these values are
inferred, which is the same conflation that put three roadways underwater at
t = 0 in an earlier product.

WHY A NEW VERSION RATHER THAN AN EDIT
--------------------------------------
v1 is what the extractor produced from the DEM. v2 is v1 with a documented,
evidence-backed repair applied on top. Keeping both means the repair can be
audited, reverted, or re-derived, and any figure or run can say which it used.

    v1  extraction, untouched
    v2  v1 + bridged dropouts

Interior SHAPES and interior ROW 0 are identical between them - this only
rewrites values inside existing arrays. That matters because every road setback
is measured from interior row 0, so the 1984 setback CSVs stay valid across the
version bump and do not need re-measuring.

TWO THINGS HAVE TO MOVE FOR v2 TO TAKE EFFECT
----------------------------------------------
hat_topo_version.resolve_version puts the EXTRACTOR'S `VERSION` literal ABOVE
the CURRENT file. Writing CURRENT alone is silently ineffective. So this script
writes both, and says so:

    dune-topo/CURRENT                        -> v2
    HAT_dune_topo_extractor.py  VERSION      -> "v2"

THE CLOBBER RISK, STATED PLAINLY
---------------------------------
Because the extractor now says VERSION = "v2", re-running it would overwrite
these arrays with a fresh, UNBRIDGED extraction. That is recoverable in one
command - this script reads v1 and rebuilds v2 - but it is silent, so
BRIDGE_MANIFEST.txt in the run folder says it too.

INPUT   dune-topo/<src>/topography, dunes
        dune-topo/<src>/nodata-audit/hole_verdicts.csv
        dune-topo/<src>/nodata-audit/bracketed_hole_cells.csv

OUTPUT  dune-topo/<dst>/topography/domain_<N>_{topography,nodata,bridged}.npy
        dune-topo/<dst>/dunes/domain_<N>_dune.npy      (copied unchanged)
        dune-topo/<dst>/BRIDGE_MANIFEST.txt
"""

import csv
import datetime
import json
import shutil
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[4]
sys.path.insert(0, str(REPO / "scripts"))
import hat_topo_version as htv  # noqa: E402
from cascade_pipeline import roadway  # noqa: E402

TOPO_PRODUCT = "1984-start"
SRC_VERSION = "v1"
DST_VERSION = "v2"
DAM_TO_M = 10.0
SENTINEL_DAM = -0.3
EXTRACTOR = (REPO / "scripts" / "input_prep" / "1-barrier3d-domains"
             / "HAT_dune_topo_extractor.py")

AUDIT_SUBDIR = "nodata-audit"


def carry_picks_forward(src_version, dst_version):
    """Give the new version its own pick set, seeded from the source version.

    WHY THIS IS PART OF THE VERSION BUMP AND NOT AN AFTERTHOUGHT

    The extractor derives its window file from the VERSION literal:

        PICK_SET    = RUN_NAME = VERSION
        WINDOW_JSON = PICKS_DIR / f"HAT_dune_search_windows_{PICK_SET}.json"

    So bumping VERSION without writing that file leaves every consumer that
    resolves picks through the extractor pointing at a path that does not
    exist. That is not hypothetical: v1 -> v2 did exactly this on 2026-08-26,
    and `HAT_road_offset_from_dune_start.py` -- which re-derives interior row 0
    from these windows in order to measure the road setbacks -- could not be
    re-run at all afterwards. The setbacks already on disk stayed correct; they
    simply could not be regenerated.

    COPYING IS THE RIGHT ANSWER HERE, NOT SHARING

    The extractor's own guard prescribes this procedure verbatim: "set
    PICK_SET = RUN_NAME and copy the old file to
    HAT_dune_search_windows_<RUN_NAME>.json first". The reason is that
    `save_windows()` writes back to WINDOW_JSON after EVERY domain during a
    pick pass, so pointing a re-pick at a shared file destroys the earlier
    version's picks one domain at a time. A per-version copy makes that
    impossible.

    And the copy is semantically honest for a BRIDGED version: bridging fills
    unsurveyed cells inside existing arrays. It does not move a dune, so it
    cannot move a dune window. v2's windows ARE v1's windows.

    PROVENANCE IS RECORDED RATHER THAN IMPLIED

    A bare copy would leave a file that looks like it was picked for this
    version. `_meta.inherited_from` says otherwise, so a later reader can tell
    an inherited pick set from a re-picked one. It goes INSIDE `_meta` on
    purpose: consumers count domains as "every key that is not `_meta`", so a
    second top-level underscore key would be counted as a domain.

    Returns:
        (destination path, "written" | "exists") or None if the source has no
        pick file to carry.
    """
    picks_dir = htv.product_dir(TOPO_PRODUCT) / "picks"
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


def cleared_holes(src_run):
    """{(domain, profile): [interior rows]} for holes judged DROPOUT."""
    a = src_run / AUDIT_SUBDIR
    verd = {(int(r["domain"]), int(r["profile"])): r["verdict"]
            for r in csv.DictReader((a / "hole_verdicts.csv").open())}
    cells = defaultdict(list)
    for r in csv.DictReader((a / "bracketed_hole_cells.csv").open()):
        k = (int(r["domain"]), int(r["profile"]))
        if verd.get(k) == "DROPOUT":
            cells[k].append(int(r["interior_row"]))
    return {k: sorted(v) for k, v in cells.items()}


def bridge_profile(col, rows):
    """Linear fill of `rows` in one profile. Returns False if not bracketed.

    Refuses rather than guesses. A hole touching either end of the array has no
    measured value on one side, so there is nothing to interpolate between -
    that is an extrapolation, which is the thing this whole exercise exists to
    avoid.
    """
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

    # --- verification --------------------------------------------------------
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

    # --- manifest ------------------------------------------------------------
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
      python scripts/input_prep/1-barrier3d-domains/nodata_audit/\\
          HAT_bridge_dropouts.py

TO REVERT
  Set dune-topo/CURRENT back to "{src_ver}" AND the extractor's VERSION to
  "{src_ver}". Both, or the extractor literal silently wins.

mean island width gained: """ + ", ".join(
        f"D{d} {g:+.0f} m" for d, g in sorted(width_gain.items())) + "\n",
        encoding="utf-8")

    # --- activate ------------------------------------------------------------
    # Picks FIRST. The moment the extractor's VERSION literal moves, every
    # caller that resolves a window file through it looks for this version's
    # pick set. Writing it after the bump would leave a window -- however
    # short -- in which the tree resolves to a path that does not exist.
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
