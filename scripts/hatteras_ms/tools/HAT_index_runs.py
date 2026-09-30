#!/usr/bin/env python3
"""
Rebuild output/raw_runs/run_index.csv from the runs themselves, with one ledger of retired runs.

    python scripts/hatteras_ms/tools/HAT_index_runs.py --check
    python scripts/hatteras_ms/tools/HAT_index_runs.py --dry-run
    python scripts/hatteras_ms/tools/HAT_index_runs.py

Derived, never appended to; vanished rows are retired to the ledger. Details: scripts/hatteras_ms/README.md.

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-09-18
"""
from __future__ import annotations

import argparse
import csv
import sys
from datetime import datetime
from pathlib import Path


# Walk up until a directory holds data/hatteras_init
def _find_root(start: Path) -> Path:
    for p in [start, *start.parents]:
        if (p / "data" / "hatteras_init").is_dir():
            return p
    raise SystemExit("could not locate the project root")


REPO = _find_root(Path(__file__).resolve())
sys.path.insert(0, str(REPO / "scripts"))

from cascade_pipeline.run_registry import (  # noqa: E402
    INDEX_KEY, legacy_arm_to_kind_tag, load_run_index, rebuild_run_index,
    sweep_family)
from site_layer.hat_topo_version import current_topo_versions  # noqa: E402

# --- CONFIG ------------------------------------------------------------------
RAW_RUNS = REPO / "output" / "raw_runs"
INDEX = RAW_RUNS / "run_index.csv"
LEDGER = RAW_RUNS / "retired_runs.csv"
ARCHIVE_GLOBS = ("run_index_archive_*.csv", "run_index_measuredreloc_*.csv")
LEDGER_FIELDS = ["retired_on", "source", "run_name", "kind", "tag", "Hs_m",
                 "timestamp", "start_year", "end_year", "source_sink_preset",
                 "topo_product", "topo_dune_version"]
# -----------------------------------------------------------------------------


# A CSV as a list of dicts, or [] if absent
def read_csv(path: Path) -> list:
    if not path.is_file():
        return []
    with open(path, newline="", encoding="utf-8") as fh:
        return list(csv.DictReader(fh))


# A row from any vintage of the index, carrying kind and tag
def with_kind_tag(row: dict) -> dict:
    row = dict(row)
    if not row.get("kind"):
        kind, tag = legacy_arm_to_kind_tag(row.get("arm", ""))
        family = sweep_family(row.get("run_name", ""))
        if kind == "matrix" and family:
            kind, tag = "sensitivity", family
        row["kind"], row["tag"] = kind, tag
    return row


# A row's index key (run_name, kind, tag)
def key_of(row) -> tuple:
    return tuple(str(row.get(k, "")).strip() for k in INDEX_KEY)


# The index as it stands, every row with kind/tag
def current_rows() -> list:
    return [with_kind_tag(r) for r in read_csv(INDEX)]


# Append retired rows
def append_ledger(entries: list) -> None:
    if not entries:
        return
    existing = [with_kind_tag(r) for r in read_csv(LEDGER)]
    seen = {(r.get("run_name"), r.get("kind"), r.get("tag"), r.get("timestamp"))
            for r in existing}
    fresh = [e for e in entries
             if (e.get("run_name"), e.get("kind"), e.get("tag"),
                 e.get("timestamp")) not in seen]
    if not fresh:
        return
    header_ok = LEDGER.is_file() and read_csv(LEDGER) and \
        "kind" in read_csv(LEDGER)[0]
    if LEDGER.is_file() and not header_ok:
        # Widen once: rewrite with the new header, rows unchanged.
        with open(LEDGER, "w", newline="", encoding="utf-8") as fh:
            w = csv.DictWriter(fh, fieldnames=LEDGER_FIELDS, extrasaction="ignore")
            w.writeheader()
            for r in existing:
                w.writerow(r)
    new = not LEDGER.is_file()
    with open(LEDGER, "a", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=LEDGER_FIELDS, extrasaction="ignore")
        if new:
            w.writeheader()
        for e in fresh:
            w.writerow(e)
    sys.stdout.write(f"  ledger: {len(fresh)} row(s) appended to {LEDGER.name}\n")


# Seed the ledger from the archived copies of the index
def adopt_archives() -> None:
    live = {key_of(r) for r in current_rows()}
    entries = []
    for pattern in ARCHIVE_GLOBS:
        for path in sorted(RAW_RUNS.glob(pattern)):
            for row in read_csv(path):
                row = with_kind_tag(row)
                if key_of(row) in live:
                    continue
                e = dict(row)
                e["retired_on"] = ""
                e["source"] = path.name
                entries.append(e)
    if not entries:
        sys.stdout.write("no archived rows that are not still on disk\n")
        return
    append_ledger(entries)
    sys.stdout.write(f"adopted {len(entries)} archived row(s) into {LEDGER.name}\n")


# Compare the index with disk and, unless asked not to, rewrite it
def rebuild(check: bool = False, dry_run: bool = False) -> int:
    before = current_rows()
    have = {key_of(r): r for r in before}

    # A dry rebuild into a scratch path says what the file WOULD hold.
    target = INDEX if not (check or dry_run) else RAW_RUNS / ".run_index_preview.csv"
    try:
        rows = rebuild_run_index(RAW_RUNS, target,
                                 current_versions=current_topo_versions())
    finally:
        if target != INDEX and target.exists():
            target.unlink()
    disk = {key_of(r): r for r in rows}

    vanished = [r for k, r in have.items() if k not in disk]
    missing = [k for k in disk if k not in have]
    sys.stdout.write(f"{len(disk)} run(s) on disk, {len(before)} row(s) in the index\n")
    if not vanished and not missing:
        sys.stdout.write("  in agreement\n")
    if vanished:
        sys.stdout.write(f"  {len(vanished)} row(s) whose run is gone from disk:\n")
        for r in vanished[:8]:
            sys.stdout.write(f"      {r.get('run_name')}  {r.get('kind')}:{r.get('tag')}  "
                             f"{r.get('timestamp')}\n")
    if missing:
        sys.stdout.write(f"  {len(missing)} run(s) on disk with no row:\n")
        for k in missing[:8]:
            sys.stdout.write(f"      {k[0]}  {k[1]}:{k[2]}\n")
    by_status = {}
    for r in rows:
        by_status[r.get("status", "")] = by_status.get(r.get("status", ""), 0) + 1
    sys.stdout.write("  status: " + ", ".join(f"{k} {v}" for k, v in sorted(by_status.items())) + "\n")

    if check:
        return 1 if (vanished or missing) else 0
    if dry_run:
        sys.stdout.write("\ndry run: nothing written\n")
        return 0
    stamp = datetime.now().strftime("%Y-%m-%d %H:%M")
    append_ledger([dict(r, retired_on=stamp, source="rebuild") for r in vanished])
    sys.stdout.write(f"\nwrote {INDEX.name}: {len(rows)} row(s)\n")
    return 0


# Run: rebuild, compare, or retire, as asked
def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    ap.add_argument("--check", action="store_true",
                    help="compare the index with the runs; exit 1 if they differ")
    ap.add_argument("--dry-run", action="store_true", help="report, write nothing")
    ap.add_argument("--adopt-archives", action="store_true",
                    help="seed the ledger from the archived index copies, then exit")
    args = ap.parse_args()
    if args.adopt_archives:
        adopt_archives()
        return
    raise SystemExit(rebuild(check=args.check, dry_run=args.dry_run))


if __name__ == "__main__":
    main()
