#!/usr/bin/env python3
r"""
HAT_migrate_run_layout.py
==============================================================================
Move output/raw_runs into the PURPOSE layout of 2026-09-16. Moves only;
nothing is deleted, rewritten or renamed in content, and run_index.csv is
rebuilt from the runs afterwards rather than edited.

THE MOVES
    <period>/<preset>/<run>                       -> matrix/<period>/<preset>/<run>
    <period>/<preset>/sweeps/<axis>/<run>_<tok>   -> sensitivity/<axis>/<period>/<preset>/<run>_<tok>
    arms/<arm>/<period>/<preset>/<run>            -> versions/<tag>/... or experiments/<tag>/...
                                                     as run_registry.LEGACY_ARMS says
    arms/waveHs<x>/1996_2010/<preset>/<run>       LEFT IN PLACE, and listed. Those are
                                                     the twelve 1996 wave cells filed by
                                                     the 09-01 rule; they are re-run as
                                                     sensitivity cells (the token back in
                                                     the name) and then deleted, because
                                                     renaming a run's files and metadata
                                                     is exactly the in-content edit this
                                                     tool refuses to make.

    Why (Hannah, 2026-09-16): a wave sweep had fanned out into twelve
    top-level arms; nineteen of thirty arms were finished one-off experiments
    nothing marked as finished; and the layout could not tell those from the
    version comparisons that are kept on purpose. A run's folder now says what
    it is FOR. The whole design is in output/raw_runs/README.md.

WHAT ELSE IT WRITES
    experiments/<tag>/NOTE.md   one per experiment set, seeded from the index
                                and the write-ups that already name these
                                runs. Says what was asked, where the answer
                                is, and whether the runs may be deleted. Not
                                overwritten if present.
    run_index.csv               rebuilt (HAT_index_runs), every row with
                                kind, tag and status.

WHY IT IS SAFE TO RUN, AND TO INTERRUPT
    Every destination is checked to be free before anything moves, so a name
    collision is reported and nothing happens. All three layouts are readable:
    run_registry.find_run_dir tries the purpose path, then the 09-10 one, then
    the flat one, so a tree that is half moved -- or interrupted here -- still
    reads. Running it twice is a no-op.

    The 2026-09-10 pass (files inside each run into figures/, animations/,
    tables/) is retained as --files, and is a no-op on a migrated run.

USAGE
    python HAT_migrate_run_layout.py --dry-run       # show every move
    python HAT_migrate_run_layout.py
    python HAT_migrate_run_layout.py --files         # the 09-10 in-run file pass
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
    ARMS_DIR, KIND_DIR, LEGACY_ARMS, SWEEPS_DIR, _PERIOD_DIR, sweep_family)

RAW_RUNS = REPO / "output" / "raw_runs"
PURPOSE_DIRS = set(KIND_DIR.values())

# What each experiment set was for, from the run index and the write-ups
# that cite it. Seeded here so the NOTE.md a folder gets on migration is not
# blank; edit the file, not this table, once the folder exists.
EXPERIMENT_NOTES = {
    "2026-09-02-pea1989": """# 2026-09-02-pea1989

**Question.** The control for the GIS 84-86 seaward-row insert (the Pea Island
1989 relocation): the 1984-2004 calibBE full-management groin-on run on
1984-start **v1** with the setback CSV as shipped, once with the prescribed
relocations (`base`) and once without (`basenoreloc`). Every insert arm it was
the control for was deleted on 2026-09-07 (unmodified topography only).

**Made by.** `scripts/hatteras_ms/experiments/HAT_run_crest_experiment.py`.

**Answer lives in.** `output/experiments/pea1989_crest/` (frozen) and
`data/hatteras_init/1-barrier3d-domains/LINEAGE.md`.

**Runs deletable?** Yes once the crest figures are no longer needed; they are
on v1 and reproducible from the script. Model state was kept because the
comparison reads roadway objects.
""",
    "2026-09-08-behindroad-copy": """# 2026-09-08-behindroad-copy

**Question.** Does the 1984 dune footprint placed directly behind NC-12 (the
dune-topo **v3** layer: v2 + copy fill behind the road) change the calibrated
1984-2004 calibBE full-management groin-on run? One run, on v3 deliberately.

**Answer lives in.** `data/hatteras_init/1-barrier3d-domains/1984-start/
2-domain-reconstruction-1984/6-result/README.md` and `dune-topo/README.md`
(RMSE 0.544 / bias +0.013 against 0.547 / +0.008 on v2).

**Runs deletable?** Keep while v3 is a candidate; it is the only v3 run of
this cell besides `versions/version-pair/v3`.
""",
    "2026-09-14-paramsplit": """# 2026-09-14-paramsplit

**Question.** Did splitting the site configuration on 2026-09-14 (the
`hatteras_site_config_prebe_20260914_*` backups mark the steps) change the
model? Four runs of the same cell, `HAT_1984_2004_calibBE_road_reloc_bdm_nogroin`
on 1984-start v1, under the old code (`oldcode`), the split with the sync on
(`syncon`), a check (`check`) and the final state (`final`).

**Answer.** All four carry interior RMSE 0.523026: bit-identical, the split
changed nothing. Compare with `tools/HAT_compare_rerun.py`.

**Runs deletable?** Yes; the answer is the four identical numbers in
`run_index.csv`, which `retired_runs.csv` keeps.
""",
    "2026-09-14-probe": """# 2026-09-14-probe

**Question.** Probes around the 2026-09-14 change that rounds a prescribed
relocation displacement to whole 10 m cells (see the note in
`hatteras_site_config.py` near `round_to_cell`), all on 1984-start v1:
`natural` and `roadonly` (zeroBE, no relocations), `v1setbacks` (roadway-only
on the v1-era setback CSV), `paired` (calibBE with relocations, rounded) and
`unrounded` (the same without rounding: RMSE 0.523102 against 0.522873).

**Answer lives in.** The recode comparison beside this folder
(`2026-09-14-recode`) and `tools/HAT_compare_rerun.py`.

**Runs deletable?** Yes; each is a single probe whose numbers are in the index.
""",
    "2026-09-14-recode": """# 2026-09-14-recode

**Question.** The whole 1984-2004 relocation arm (12 runs: three presets x
bdm/nobdm x groin/nogroin, relocations on) re-run under the 2026-09-14 code
on the same v1 topography, so the difference against the stored runs is CODE
only -- chiefly the whole-cell rounding of prescribed displacements, which
moved GIS 11 over the drowning line by one cell.

**Made by.** `tools/HAT_rerun_arm.py`; compared with `tools/HAT_compare_rerun.py`.

**Answer lives in.** `hatteras_site_config.py` (the note above `CELL_M`) and
the compare tool's output.

**Runs deletable?** Yes once the stored v1 matrix is itself archived: both
sides of the comparison go together.
""",
    "2026-09-14-currency": """# 2026-09-14-currency

**Question.** Is the calibrated calibBE pair (1984-2004 on v2, 2004-2024 on
2004-start v1) current under today's code? Both re-run and differenced per
domain against the stored matrix runs.

**Answer.** The 1984 run reproduces to exactly zero; the 2004 run moves 56 of
90 domains by at most 0.0084 m/yr (RMSE 0.580353 -> 0.579999). Written up in
`output/comparisons/hindcast_calibrated/README.md`, "Currency, checked
2026-09-14".

**Runs deletable?** Yes; the write-up holds the numbers.
""",
}


def run_dirs() -> list[Path]:
    """Every run folder under raw_runs, found by its metadata file."""
    return sorted({p.parent for p in RAW_RUNS.rglob("*_run_metadata.json")})


def destination(run_dir: Path):
    """(destination, reason) for one run folder; destination None = stay."""
    rel = run_dir.relative_to(RAW_RUNS).parts
    if rel[0] in PURPOSE_DIRS:
        return None, "already in the purpose layout"
    i = next((k for k, part in enumerate(rel) if _PERIOD_DIR.fullmatch(part)),
             None)
    if i is None or len(rel) < i + 3:
        return None, "not a run path this tool knows"
    period, preset, run = rel[i], rel[i + 1], rel[-1]
    head = rel[:i]
    family = sweep_family(run)

    if head and head[0] == ARMS_DIR:
        arm = "/".join(head[1:])
    elif head:
        arm = "/".join(head)              # pre-09-10 loose arm
    else:
        arm = ""

    if arm:
        if arm in LEGACY_ARMS:
            kind, tag = LEGACY_ARMS[arm]
            return RAW_RUNS / KIND_DIR[kind] / tag / period / preset / run, kind
        if any(arm.startswith(f) for f in ("waveHs", "waveTp", "waveahf",
                                           "waveasym")):
            return None, ("a wave cell filed by arm (09-01 rule); re-run as a "
                          "sensitivity cell and delete")
        return (RAW_RUNS / KIND_DIR["experiment"] / arm / period / preset / run,
                "experiment (arm not in LEGACY_ARMS)")
    if family:
        return (RAW_RUNS / KIND_DIR["sensitivity"] / family / period / preset / run,
                "sensitivity")
    return RAW_RUNS / KIND_DIR["matrix"] / period / preset / run, "matrix"


def tree_moves():
    """(source, destination) for each run folder that changes place, and the
    folders that stay with why."""
    moves, stays = [], []
    for d in run_dirs():
        dest, why = destination(d)
        if dest is None:
            stays.append((d, why))
        elif dest != d:
            moves.append((d, dest))
    return moves, stays


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


def write_notes(dry_run: bool) -> int:
    """A NOTE.md per experiment set that has none. Returns how many."""
    root = RAW_RUNS / KIND_DIR["experiment"]
    written = 0
    if not root.is_dir():
        return 0
    for folder in sorted(p for p in root.iterdir() if p.is_dir()):
        note = folder / "NOTE.md"
        if note.exists():
            continue
        text = EXPERIMENT_NOTES.get(folder.name) or (
            f"# {folder.name}\n\n**Question.** (not recorded at migration; "
            f"fill in)\n\n**Answer lives in.**\n\n**Runs deletable?**\n")
        if not dry_run:
            note.write_text(text, encoding="utf-8")
        written += 1
    return written


def show(moves, limit, label) -> None:
    sys.stdout.write(f"\n{label}: {len(moves)}\n")
    for src, dst in moves[:limit]:
        s = src.relative_to(RAW_RUNS).as_posix() if isinstance(src, Path) else src
        t = dst.relative_to(RAW_RUNS).as_posix() if isinstance(dst, Path) else dst
        sys.stdout.write(f"    {s}\n      -> {t}\n")
    if len(moves) > limit:
        sys.stdout.write(f"    ... {len(moves) - limit} more\n")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[1],
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--dry-run", action="store_true", help="show the moves, change nothing")
    ap.add_argument("--files", action="store_true",
                    help="also run the 09-10 in-run file pass (figures/, tables/ ...)")
    ap.add_argument("--show", type=int, default=12, help="how many moves to list (default 12)")
    a = ap.parse_args()

    files = file_moves() if a.files else []
    tree, stays = tree_moves()

    problems = check_free(files) + check_free(tree)
    if problems:
        sys.stdout.write("REFUSING TO MOVE -- destinations are not free:\n")
        for p in problems[:20]:
            sys.stdout.write(f"    {p}\n")
        raise SystemExit(1)

    if files:
        show(files, a.show, "files to move inside run folders")
    show(tree, a.show, "run folders to move")
    left = [(d, why) for d, why in stays if "already" not in why]
    if left:
        sys.stdout.write(f"\nleft in place: {len(left)}\n")
        for d, why in left[:a.show]:
            sys.stdout.write(f"    {d.relative_to(RAW_RUNS).as_posix()}\n      {why}\n")

    if a.dry_run:
        n = write_notes(dry_run=True)
        sys.stdout.write(f"\nwould write {n} NOTE.md file(s)\n")
        sys.stdout.write("dry run: nothing was moved\n")
        return
    if not files and not tree:
        sys.stdout.write("\nnothing to move; the tree is already migrated\n")
    for src, dst in files:
        dst.parent.mkdir(parents=True, exist_ok=True)
        shutil.move(os.fspath(src), os.fspath(dst))
    for src, dst in tree:
        dst.parent.mkdir(parents=True, exist_ok=True)
        shutil.move(os.fspath(src), os.fspath(dst))
    emptied = prune_empty(RAW_RUNS)
    notes = write_notes(dry_run=False)
    sys.stdout.write(f"\nmoved {len(files)} file(s) and {len(tree)} run folder(s); "
                     f"removed {emptied} empty director{'ies' if emptied != 1 else 'y'}; "
                     f"wrote {notes} NOTE.md file(s)\n")

    # The index is DERIVED: rebuild it from the moved runs rather than patch it.
    from HAT_index_runs import rebuild  # beside this file
    sys.stdout.write("\nrebuilding run_index.csv\n")
    rebuild()


if __name__ == "__main__":
    sys.path.insert(0, str(Path(__file__).resolve().parent))
    main()
