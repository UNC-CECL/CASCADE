#!/usr/bin/env python3
r"""
HAT_migrate_run_layout.py
==============================================================================
Move output/raw_runs into the 2026-09-10 layout. Moves only; nothing is
deleted, rewritten or renamed in content.

TWO MOVES, ONE PASS

  1. THE TREE, around each run folder
         <period>/<preset>/<run>_waveHs3   ->  <period>/<preset>/sweeps/waveHs/<run>_waveHs3
         <arm>/<period>/<preset>/<run>     ->  arms/<arm>/<period>/<preset>/<run>
     91 of 164 runs vary one parameter from a scenario run, and the runner
     already encodes which one as a trailing token on the run name, so a
     preset folder can show its scenario runs and one sweeps/ folder instead
     of 41 siblings. Arms move under arms/ so the root shows the two hindcast
     periods and one folder of side experiments.

  2. THE FILES, inside each run folder
         <run>_shoreline_change_rate_REAL_DOMAINS_ONLY.png -> figures/shoreline_change_rate.png
         <run>_shoreline_position_D1-90.gif                -> animations/position_domains_1-90.gif
         <run>_shoreline_change_rate.csv                   -> tables/shoreline_change_rate.csv
     and so on. The archive, the shoreline matrix and both metadata files stay
     at the run root, keeping their run-name prefix: they are globbed across
     runs and travel outside the folder. Files that move lose the prefix, which
     the folder already carries -- worth 56 characters of path, which matters
     because the longest path here was 229 and Windows stops at 260.

WHY IT IS SAFE TO RUN, AND TO INTERRUPT
    Every destination is checked to be free before anything moves, so a name
    collision is reported and nothing happens. Both layouts are readable:
    run_registry.find_run_dir and run_layout.resolve each try the new location
    then the old, so a tree that is half moved -- or interrupted here -- still
    reads. Running it twice is a no-op.

    run_index.csv is NOT touched. Its `arm` column is a logical name, not a
    path, so it is already correct.

USAGE
    python HAT_migrate_run_layout.py --dry-run       # show every move
    python HAT_migrate_run_layout.py --tree-only
    python HAT_migrate_run_layout.py --files-only
    python HAT_migrate_run_layout.py
==============================================================================
"""
from __future__ import annotations

import argparse
import os
import shutil
import sys
from pathlib import Path


def _find_root(start: Path) -> Path:
    for p in [start, *start.parents]:
        if (p / "data" / "hatteras_init").is_dir():
            return p
    raise SystemExit("could not locate the project root")


REPO = _find_root(Path(__file__).resolve())
sys.path.insert(0, str(REPO / "scripts"))

from cascade_pipeline.run_layout import plan_run  # noqa: E402
from cascade_pipeline.run_registry import (  # noqa: E402
    ARMS_DIR, SWEEPS_DIR, sweep_family)

RAW_RUNS = REPO / "output" / "raw_runs"
_PERIOD = ("1984_2004", "2004_2024")


def run_dirs() -> list[Path]:
    """Every run folder under raw_runs, found by its metadata file."""
    return sorted({p.parent for p in RAW_RUNS.rglob("*_run_metadata.json")})


def tree_moves() -> list[tuple[Path, Path]]:
    """(source, destination) for each run FOLDER that changes place."""
    moves = []
    for d in run_dirs():
        rel = d.relative_to(RAW_RUNS).parts
        if ARMS_DIR in rel or SWEEPS_DIR in rel:
            continue                                  # already migrated
        # rel is (<arm parts...>, <period>, <preset>, <run>)
        try:
            i = next(k for k, part in enumerate(rel) if part in _PERIOD)
        except StopIteration:
            continue                                  # not a run path we know
        arm_parts, tail = rel[:i], rel[i:]
        if len(tail) != 3:
            continue
        period, preset, run = tail
        family = sweep_family(run)
        dest = RAW_RUNS
        if arm_parts:
            dest = dest / ARMS_DIR
            for part in arm_parts:
                dest = dest / part
        dest = dest / period / preset
        if family:
            dest = dest / SWEEPS_DIR / family
        dest = dest / run
        if dest != d:
            moves.append((d, dest))
    return moves


def file_moves() -> list[tuple[Path, Path]]:
    """(source, destination) for each FILE that changes place inside a run."""
    moves = []
    for d in run_dirs():
        moves.extend(plan_run(d))
    return moves


def check_free(moves) -> list[str]:
    """Destinations that are already occupied, or collide with each other."""
    problems, seen = [], {}
    for src, dst in moves:
        if dst.exists():
            problems.append(f"destination exists: {dst}")
        if dst in seen:
            problems.append(f"two sources want {dst}: {seen[dst]} and {src}")
        seen[dst] = src
    return problems


def prune_empty(root: Path) -> int:
    """Remove directories left empty by the tree move. Never removes a file."""
    removed = 0
    for d in sorted((p for p in root.rglob("*") if p.is_dir()),
                    key=lambda p: len(p.parts), reverse=True):
        try:
            next(d.iterdir())
        except StopIteration:
            d.rmdir()
            removed += 1
    return removed


def show(moves, limit, label) -> None:
    sys.stdout.write(f"\n{label}: {len(moves)}\n")
    for src, dst in moves[:limit]:
        s = src.relative_to(RAW_RUNS).as_posix()
        t = dst.relative_to(RAW_RUNS).as_posix()
        sys.stdout.write(f"    {s}\n      -> {t}\n")
    if len(moves) > limit:
        sys.stdout.write(f"    ... {len(moves) - limit} more\n")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[1],
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--dry-run", action="store_true", help="show the moves, change nothing")
    ap.add_argument("--tree-only", action="store_true", help="move run folders only")
    ap.add_argument("--files-only", action="store_true", help="move files inside runs only")
    ap.add_argument("--show", type=int, default=8, help="how many moves to list (default 8)")
    a = ap.parse_args()

    # FILES FIRST, then the tree: moving a file into a folder that is about to
    # move is fine, but planning file moves after the tree move would have to
    # re-find every run.
    files = [] if a.tree_only else file_moves()
    tree = [] if a.files_only else tree_moves()

    problems = check_free(files) + check_free(tree)
    if problems:
        sys.stdout.write("REFUSING TO MOVE -- destinations are not free:\n")
        for p in problems[:20]:
            sys.stdout.write(f"    {p}\n")
        raise SystemExit(1)

    show(files, a.show, "files to move inside run folders")
    show(tree, a.show, "run folders to move")

    if a.dry_run:
        sys.stdout.write("\ndry run: nothing was moved\n")
        return
    if not files and not tree:
        sys.stdout.write("\nnothing to do; the tree is already migrated\n")
        return

    for src, dst in files:
        dst.parent.mkdir(parents=True, exist_ok=True)
        shutil.move(os.fspath(src), os.fspath(dst))
    for src, dst in tree:
        dst.parent.mkdir(parents=True, exist_ok=True)
        shutil.move(os.fspath(src), os.fspath(dst))

    emptied = prune_empty(RAW_RUNS)
    sys.stdout.write(f"\nmoved {len(files)} file(s) and {len(tree)} run folder(s); "
                     f"removed {emptied} empty director{'ies' if emptied != 1 else 'y'}\n")
    sys.stdout.write("run_index.csv was not touched: its `arm` column is a "
                     "logical name, not a path\n")


if __name__ == "__main__":
    main()
