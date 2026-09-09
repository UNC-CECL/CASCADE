#!/usr/bin/env python3
r"""
HAT_build_footprint_version.py
==============================================================================
Build a 1984-start dune-topo VERSION from the footprint table: the symmetric
1984 footprint (rows added where the 1984 dune line lay seaward of the 1997
line, removed where it lay landward), placed BEHIND THE ROAD (or behind the
crest row where there is no model road) and filled by the COPY rule.

WHAT IS WRITTEN (dune-topo/<dst>/)
    topography/domain_<N>_topography.npy   v2's array with the block inserted
                                           (a copy of the N rows that follow the
                                           insert point) or the rows removed
    topography/domain_<N>_nodata.npy       the same row operation on the mask
    dunes/domain_<N>_dune.npy              copied unchanged: the dune stays put
    RoadSetback_1984_dunestart.csv         the 1984 setbacks: setback_new_m from
                                           the footprint table, the road measured
                                           against the 1984 dune line in the
                                           model's row-0 convention, (road - row 0)
                                           + shift per profile. Never negative, so
                                           NO FLOOR. With the rows behind the road
                                           this moves the model's road N rows
                                           inland; where rows were removed in
                                           front of the road it lands on the old
                                           pavement's first row or the row
                                           seaward of it (2026-09-08). Domains
                                           the footprint has no setback for keep
                                           v2's value.
    HAT_footprint_audit.csv                what was done to every domain
    RUN_MANIFEST.txt, README.md            provenance and the rules

THE RULES, as decided with Hannah 2026-09-07/08 (HAT_footprint_1984.py and
HAT_fill_copy_scope.py carry the argument; this script only applies them)
    N               n_cells in footprint_1984_by_domain.csv
    insert point    insert_row_behind_road: int(setback_new/10) + 2, behind
                    the model's two roadway rows AS PLACED under the 1984
                    setback; crest_row + 1 where there is no model road
                    (GIS 1-5, 8)
    add             block = z[r : r+N] copied cell by cell, inserted at r
    remove          rows r .. r+|N|-1 deleted, r = int(setback_v2/10) - |N|: the
                    |N| rows directly SEAWARD of today's roadway rows (Hannah,
                    2026-09-08: the rows come out of the interior in front of
                    the road); the road's cells and all behind them are kept
    dune array      unchanged
    setback         setback_new_m (the 1984 measurement); v2's value where absent

VERIFIED AFTER WRITING
    unchanged domains are byte-identical to v2; a changed domain has exactly
    rows_before + N rows, its rows before the insert point are identical to
    v2, and its block equals the rows that follow it (add) or its tail is
    v2's tail (remove).

USAGE
    python HAT_build_footprint_version.py --dst-version v3
    python HAT_build_footprint_version.py --dst-version v3 --overwrite
==============================================================================
"""
from __future__ import annotations

import argparse
import csv
import shutil
import sys
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd


def _find_root(start: Path) -> Path:
    for p in [start, *start.parents]:
        if (p / "data" / "hatteras_init").is_dir():
            return p
    raise SystemExit("could not locate the project root")


REPO = _find_root(Path(__file__).resolve())
INIT = REPO / "data" / "hatteras_init"
sys.path.insert(0, str(REPO / "scripts"))
from hat_topo_version import array_name, dune_topo_root, topo_dirs, insert_scope_step# noqa: E402

PRODUCT = "1984-start"
SRC_VERSION = "v2"
SCOPE = INIT / "1-barrier3d-domains" / PRODUCT / "2-domain-reconstruction-1984"
FOOTPRINT_CSV = insert_scope_step(PRODUCT, "2-extent") / "footprint_1984_by_domain.csv"
KINDS = ("topography", "nodata")          # the row operation applies to both


def apply(z: np.ndarray, n: int, r: int) -> tuple[np.ndarray, str]:
    if n > 0:
        block = z[r:r + n].copy()
        if block.shape[0] < n:
            raise SystemExit(f"window {r}..{r + n - 1} runs off the array ({z.shape[0]} rows)")
        return np.concatenate([z[:r], block, z[r:]], axis=0), f"copy of rows {r}..{r + n - 1}"
    if n < 0:
        return np.concatenate([z[:r], z[r - n:]], axis=0), f"rows {r}..{r - n - 1} removed"
    return z, ""


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--dst-version", required=True)
    ap.add_argument("--overwrite", action="store_true")
    args = ap.parse_args()

    src_topo, src_dune, src_name = topo_dirs(PRODUCT, override=SRC_VERSION)
    root = dune_topo_root(PRODUCT)
    dst = root / args.dst_version
    if dst.exists():
        if not args.overwrite:
            raise SystemExit(f"\n{dst} exists. Pass --overwrite to replace it.\n")
        shutil.rmtree(dst)
    (dst / "topography").mkdir(parents=True)
    (dst / "dunes").mkdir()

    tab = pd.read_csv(FOOTPRINT_CSV).set_index("domain")
    audit = []
    changed = 0
    for d, r in tab.iterrows():
        n = int(r["n_cells"])
        ins = int(r["insert_row_behind_road"]) if n != 0 else -1
        rec = {"domain": int(d), "n_cells": n, "insert_anchor": r.get("insert_anchor", ""),
               "insert_row": ins, "rows_before": None, "rows_after": None, "operation": "unchanged"}
        for kind in KINDS:
            src = src_topo / array_name(kind, d)
            a = np.load(src)
            if kind == "topography":
                rec["rows_before"] = int(a.shape[0])
            b, what = apply(a, n, ins)
            np.save(dst / "topography" / array_name(kind, d), b)
            if kind == "topography":
                rec["rows_after"] = int(b.shape[0])
                if what:
                    rec["operation"] = what
                    changed += 1
                # --- verify -------------------------------------------------
                if n == 0:
                    assert np.array_equal(a, b)
                else:
                    assert b.shape[0] == a.shape[0] + n, (d, a.shape, b.shape, n)
                    assert np.array_equal(a[:ins], b[:ins]), (d, "head changed")
                    if n > 0:
                        assert np.array_equal(b[ins:ins + n], b[ins + n:ins + 2 * n]), (d, "block != source")
                    else:
                        assert np.array_equal(b[ins:], a[ins - n:]), (d, "tail changed")
        shutil.copy2(src_dune / array_name("dune", d), dst / "dunes" / array_name("dune", d))
        audit.append(rec)

    with open(dst / "HAT_footprint_audit.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(audit[0].keys()))
        w.writeheader()
        w.writerows(audit)
    # --- the setback CSV: two rows, domain ids then values (the model-facing
    # format hatteras_site_config reads). Start from v2's so the domain list and
    # the format are exactly what the runner expects, then replace every value
    # the footprint has a 1984 setback for.
    src_csv = root / SRC_VERSION / "RoadSetback_1984_dunestart.csv"
    rows = list(csv.reader(open(src_csv, newline="")))
    ids = [int(float(x)) for x in rows[0] if x.strip()]
    vals = [float(x) for x in rows[1] if x.strip()]
    new_vals, changed_sb, kept = [], [], []
    for gid, v in zip(ids, vals):
        sb = tab.loc[gid, "setback_new_m"] if gid in tab.index else np.nan
        if np.isfinite(sb):
            new_vals.append(float(sb))
            if abs(float(sb) - v) > 0.05:
                changed_sb.append((gid, v, float(sb)))
        else:
            new_vals.append(v)
            kept.append(gid)
    assert min(new_vals) >= 0.0, "a 1984 setback came out negative - the model cannot index that"
    with open(dst / "RoadSetback_1984_dunestart.csv", "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow([f"{g:.3f}" for g in ids])
        w.writerow([f"{v:.3f}" for v in new_vals])
    print(f"  setback CSV: {len(changed_sb)} of {len(ids)} road domains take the 1984 setback; "
          f"{len(kept)} keep v2's ({kept}); min {min(new_vals):.1f} m, no floor")
    for gid, a, b in changed_sb:
        if gid in (16, 84, 85, 86):
            print(f"    GIS {gid}: {a:.0f} -> {b:.1f} m")

    add = [a for a in audit if a["n_cells"] > 0]
    rem = [a for a in audit if a["n_cells"] < 0]
    stamp = datetime.now().strftime("%Y-%m-%d %H:%M")
    (dst / "RUN_MANIFEST.txt").write_text(
        "=" * 78 + f"\n1984 FOOTPRINT, BEHIND THE ROAD, COPY FILL  --  {args.dst_version}\n" + "=" * 78 + "\n\n"
        f"written    : {stamp}\n"
        f"source     : {PRODUCT}/{src_name}  (arrays read, never modified)\n"
        f"footprint  : {FOOTPRINT_CSV.relative_to(REPO)}\n"
        f"builder    : scripts/input_prep/1-barrier3d-domains/2-domain-reconstruction-1984/5-build/HAT_build_footprint_version.py\n\n"
        f"domains changed : {changed} of {len(audit)}   "
        f"(+{sum(a['n_cells'] for a in add)} rows at {len(add)} domains, "
        f"-{-sum(a['n_cells'] for a in rem)} rows at {len(rem)} domains)\n\n"
        "RULES\n"
        "  N            trunc(median paired 1997-1984 dune-line shift / 10 m), both signs\n"
        "  insert point int(setback_new/10) + 2, behind the model's two roadway rows AS\n"
        "               PLACED under the 1984 setback; crest_row + 1 with no model road\n"
        "  add          block = a copy of the N rows that follow the insert point\n"
        "  remove       rows insert..insert+|N|-1 deleted, insert = int(setback_v2/10) - |N|:\n"
        "               the |N| rows directly seaward of today's roadway rows; the road's\n"
        "               cells and everything behind them kept\n"
        "  dune array   unchanged; row 0 unmoved\n"
        "  setback      setback_new_m: the 1984 road against the 1984 dune line, row-0\n"
        "               convention, unrounded, NO FLOOR (min " + f"{min(new_vals):.1f}" + " m). The model's\n"
        "               road therefore sits on measured cells at its 1984 distance from\n"
        "               the crest, and the block is directly behind it. v2's value\n"
        "               where the footprint has none.\n"
        "  nodata mask  the same row operation as the topography\n",
        encoding="utf-8")
    (dst / "README.md").write_text(
        f"# `{args.dst_version}` — the 1984 footprint behind the road, copy fill\n\n"
        f"Built {stamp} by `HAT_build_footprint_version.py` from `{src_name}` and "
        f"`2-domain-reconstruction-1984/2-extent/footprint_1984_by_domain.csv`. **Nothing is re-measured here**: the "
        f"footprint (`HAT_footprint_1984.py`) and the fill (`HAT_fill_copy_scope.py`) carry the "
        f"argument; this folder applies them.\n\n"
        f"| | |\n|---|---|\n"
        f"| what | `{src_name}` + the symmetric 1984 footprint: rows added directly behind NC-12 as placed, filled by copying the N rows that follow the insert point; rows removed directly in front of today's NC-12 rows |\n"
        f"| changed | {changed} of {len(audit)} domains: +{sum(a['n_cells'] for a in add)} rows at {len(add)}, "
        f"−{-sum(a['n_cells'] for a in rem)} rows at {len(rem)} |\n"
        f"| dune array | unchanged |\n"
        f"| setback CSV | the 1984 setbacks (`setback_new_m`), no floor: {len(changed_sb)} road domains change, "
        f"GIS 85/86 go from 0 to {tab.loc[85, 'setback_new_m']:.0f}/{tab.loc[86, 'setback_new_m']:.0f} m; "
        f"the road sits on measured cells and the block is directly behind it |\n"
        f"| audit | `HAT_footprint_audit.csv` |\n\n"
        f"Verified on write: unchanged domains byte-identical to `{src_name}`; changed domains have "
        f"rows_before + N rows, identical rows before the insert point, a block equal to the rows "
        f"that follow it (add) or an unchanged tail (remove).\n\n"
        f"To run on it for one run: `HAT_TOPO_VERSION_1984_START={args.dst_version}` AND copy this folder's "
        f"`RoadSetback_1984_dunestart.csv` over the forcing-tree one for the run (restore after) - "
        f"`hatteras_site_config.py` hardcodes that path. `CURRENT` is not changed by building.\n",
        encoding="utf-8")
    (dst / "RoadSetback_README.md").write_text(
        f"`RoadSetback_1984_dunestart.csv`: the 1984 setbacks from "
        f"`2-domain-reconstruction-1984/2-extent/footprint_1984_by_domain.csv` (`setback_new_m` = (road - row 0) + the paired "
        f"1997-1984 dune-line shift, per profile, median per domain), in the two-row model-facing format "
        f"copied from `{src_name}`'s file. No value is negative, so nothing is floored. Domains the "
        f"footprint has no setback for keep `{src_name}`'s value: {kept}.\n", encoding="utf-8")
    print(f"wrote {dst}\n  {changed} of {len(audit)} domains changed; "
          f"+{sum(a['n_cells'] for a in add)} rows at {len(add)}, -{-sum(a['n_cells'] for a in rem)} rows at {len(rem)}")
    print("  verified: unchanged domains identical, blocks equal their source, heads untouched")


if __name__ == "__main__":
    main()
