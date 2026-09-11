#!/usr/bin/env python3
r"""
HAT_index_runs.py
==============================================================================
Rebuild output/raw_runs/run_index.csv from the runs themselves, and keep a
single ledger of runs that have been retired.

WHY THE INDEX IS DERIVED (Hannah, 2026-09-11)
    Every run writes its own metadata JSON, and that is the source of truth:
    it is produced by the run, beside the run, from the values the run used.
    The index restates 46 of those 64 facts in one table so that a question
    across runs -- which topography, which preset, what skill -- is one read
    instead of 164.

    A restatement can drift from what it restates. Until now the index was
    maintained incrementally: every run rewrote the whole file through pandas,
    which is how a smoke test on 2026-09-10 truncated floats in the last
    digit of five columns in an UNRELATED row. Rebuilding from the runs makes
    drift impossible rather than unlikely, and `--check` turns "is it right"
    into a question with an answer.

WHAT IS NOT DERIVED
    A run deleted from disk leaves no metadata to rebuild from, so a rebuild
    would drop its row silently. `retired_runs.csv` is the append-only record
    of those: a row that was in the index, whose run is gone. It is written
    here and never rewritten, so the history of what was removed survives a
    rebuild. `--adopt-archives` seeds it from the pre-2026-09-11 archived
    copies of the index.

USAGE
    python HAT_index_runs.py --check          # compare, change nothing
    python HAT_index_runs.py --dry-run
    python HAT_index_runs.py                  # rebuild, retiring vanished rows
    python HAT_index_runs.py --adopt-archives # seed the ledger from the copies
==============================================================================
"""
from __future__ import annotations

import argparse
import csv
import json
import sys
from datetime import datetime
from pathlib import Path


def _find_root(start: Path) -> Path:
    for p in [start, *start.parents]:
        if (p / "data" / "hatteras_init").is_dir():
            return p
    raise SystemExit("could not locate the project root")


REPO = _find_root(Path(__file__).resolve())
sys.path.insert(0, str(REPO / "scripts"))

RAW_RUNS = REPO / "output" / "raw_runs"
INDEX = RAW_RUNS / "run_index.csv"
LEDGER = RAW_RUNS / "retired_runs.csv"
ARCHIVE_GLOBS = ("run_index_archive_*.csv", "run_index_measuredreloc_*.csv")

# What identifies a run. Not the name alone: forcing that is not part of the
# scenario (Hs) scopes the output directory rather than adding a name token,
# so two runs can share a name and differ in what they were forced with.
KEY = ("run_name", "Hs_m", "arm")


def read_csv(path: Path) -> list:
    if not path.is_file():
        return []
    with open(path, newline="", encoding="utf-8") as fh:
        return list(csv.DictReader(fh))


def key_of(row) -> tuple:
    return tuple(str(row.get(k, "")).strip() for k in KEY)


def rows_on_disk() -> dict:
    """{key: index row} rebuilt from each run's metadata, keyed as the index is.

    The row is taken from the index where the run still exists, because the
    index is the shape every consumer already reads and the metadata is its
    source. Where a run is on disk but absent from the index -- a run made
    while the index was missing -- the identity columns are filled from the
    metadata so it is not lost.
    """
    indexed = {key_of(r): r for r in read_csv(INDEX)}
    out = {}
    for meta in sorted(RAW_RUNS.rglob("*_run_metadata.json")):
        try:
            d = json.load(open(meta, encoding="utf-8"))
        except (json.JSONDecodeError, OSError):
            continue
        ident = d.get("identity", {})
        wave = d.get("wave climate", {})
        rel = meta.parent.relative_to(RAW_RUNS).as_posix()
        arm = ""
        parts = rel.split("/")
        if parts and parts[0] == "arms":
            arm = "/".join(parts[1:-3]) if len(parts) > 4 else parts[1]
        probe = (str(ident.get("run_name") or ""),
                 str(wave.get("wave_height_m") or ""),
                 arm or "calibration")
        row = None
        for cand in (probe, (probe[0], probe[1], ""), (probe[0], "", probe[2])):
            if cand in indexed:
                row = dict(indexed[cand])
                break
        if row is None:
            row = {k: "" for k in (indexed and next(iter(indexed.values())) or {})}
            row.update({"run_name": probe[0], "Hs_m": probe[1], "arm": probe[2],
                        "timestamp": ident.get("timestamp") or ""})
        out[key_of(row)] = row
    return out


def write_rows(path: Path, rows: list, fieldnames: list) -> None:
    """Write with the csv module, not pandas: pandas round-trips every float
    through repr and that silently rewrote unrelated rows (2026-09-10)."""
    with open(path, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=fieldnames, extrasaction="ignore")
        w.writeheader()
        for r in rows:
            w.writerow(r)


def append_ledger(entries: list) -> None:
    """Append retired rows. Never rewrites: the ledger is the one record of
    what was removed, and a rebuild must not be able to erase it."""
    if not entries:
        return
    existing = read_csv(LEDGER)
    fields = ["retired_on", "source", "run_name", "Hs_m", "arm", "timestamp",
              "start_year", "end_year", "source_sink_preset", "topo_product",
              "topo_dune_version"]
    seen = {(r.get("run_name"), r.get("Hs_m"), r.get("arm"), r.get("timestamp"))
            for r in existing}
    fresh = [e for e in entries
             if (e.get("run_name"), e.get("Hs_m"), e.get("arm"), e.get("timestamp")) not in seen]
    if not fresh:
        return
    new = not LEDGER.is_file()
    with open(LEDGER, "a", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, extrasaction="ignore")
        if new:
            w.writeheader()
        for e in fresh:
            w.writerow(e)
    sys.stdout.write(f"  ledger: {len(fresh)} row(s) appended to {LEDGER.name}\n")


def adopt_archives() -> None:
    """Seed the ledger from the archived copies of the index.

    Those copies are the only record of runs deleted before the ledger
    existed. Folding them in makes one file answer "what was retired", and
    the copies stay on disk until they are removed deliberately (they are
    tracked, so git holds them either way).
    """
    live = set(rows_on_disk())
    entries = []
    for pattern in ARCHIVE_GLOBS:
        for path in sorted(RAW_RUNS.glob(pattern)):
            for row in read_csv(path):
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

    current = read_csv(INDEX)
    if not current:
        raise SystemExit(f"no index at {INDEX}")
    fieldnames = list(current[0])
    disk = rows_on_disk()
    have = {key_of(r) for r in current}

    vanished = [r for r in current if key_of(r) not in disk]
    missing = [k for k in disk if k not in have]

    sys.stdout.write(f"{len(disk)} run(s) on disk, {len(current)} row(s) in the index\n")
    if not vanished and not missing:
        sys.stdout.write("  in agreement\n")
    if vanished:
        sys.stdout.write(f"  {len(vanished)} row(s) whose run is gone from disk:\n")
        for r in vanished[:8]:
            sys.stdout.write(f"      {r.get('run_name')}  {r.get('timestamp')}\n")
    if missing:
        sys.stdout.write(f"  {len(missing)} run(s) on disk with no row:\n")
        for k in missing[:8]:
            sys.stdout.write(f"      {k[0]}\n")

    if args.check:
        raise SystemExit(1 if (vanished or missing) else 0)
    if args.dry_run:
        sys.stdout.write("\ndry run: nothing written\n")
        return

    # Order: as the index has them, then anything new, so a rebuild of an
    # unchanged tree is a byte-for-byte no-op.
    ordered = [disk[key_of(r)] for r in current if key_of(r) in disk]
    ordered += [disk[k] for k in disk if k not in have]

    stamp = datetime.now().strftime("%Y-%m-%d %H:%M")
    append_ledger([dict(r, retired_on=stamp, source="rebuild") for r in vanished])
    write_rows(INDEX, ordered, fieldnames)
    sys.stdout.write(f"\nwrote {INDEX.name}: {len(ordered)} row(s)\n")


if __name__ == "__main__":
    main()
