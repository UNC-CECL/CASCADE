"""Where a run's files live inside its run folder, and what they are called.

WHY THIS MODULE EXISTS
    A run folder used to be thirteen files in a heap, every one of them
    prefixed with the full run name:

        HAT_1984_2004_calibBE_road_bdm_groin/
            HAT_1984_2004_calibBE_road_bdm_groin.npz
            HAT_1984_2004_calibBE_road_bdm_groin_annotated.png
            HAT_1984_2004_calibBE_road_bdm_groin_shoreline_change_rate.csv
            ... ten more

    Since 2026-09-10 the output is sorted by kind, and the files inside those
    subfolders drop the run-name prefix, which the folder already carries:

        HAT_1984_2004_calibBE_road_bdm_groin/
            HAT_..._run_metadata.json      identity and state stay at the root:
            HAT_..._run_metadata.txt       these are globbed ACROSS runs and
            HAT_....npz                    travel outside the folder, so they
            HAT_..._shoreline_matrix.npy   keep the prefix
            figures/     shoreline_change_rate.png, ..._with_buffers.png
            animations/  displacement_domains_1-90.gif, ...
            tables/      shoreline_change_rate.csv, groin_diagnostics.csv,
                         road_management.csv

    The rename is not cosmetic. The longest path in the runs tree was 229
    characters and Windows gives up at 260; adding a subfolder took it to 240,
    leaving no headroom for a longer run name. The prefix was 56 of those
    characters, and dropping it inside the subfolders brings the worst case
    back to about 184.

THE FALLBACK IS THE POINT
    `resolve()` looks for the new location first and falls back to the old
    flat one. So a half-migrated tree always reads correctly, the migration
    can be interrupted, and a run folder restored from an old backup still
    works. Nothing in the codebase should join a run-folder filename by hand;
    call `resolve()` to read and `write_path()` to write.

USAGE
    from cascade_pipeline.run_layout import resolve, write_path, migrate_run

    p = resolve(run_dir, "rate_csv", run_name)          # read, either layout
    p = write_path(run_dir, "figure_rate", run_name)    # write, new layout
    migrate_run(run_dir)                                # move an old folder
"""
from __future__ import annotations

import os
import re
import shutil
from pathlib import Path

FIGURES = "figures"
ANIMATIONS = "animations"
TABLES = "tables"
SUBFOLDERS = (FIGURES, ANIMATIONS, TABLES)

# kind -> (subfolder or "" for the run root,
#          new filename (no run-name prefix unless it is a root file),
#          legacy filename as a "{run}" template)
KINDS = {
    # --- state and identity: stay at the root, keep the prefix -------------
    "archive":       ("", "{run}.npz", "{run}.npz"),
    "matrix":        ("", "{run}_shoreline_matrix.npy", "{run}_shoreline_matrix.npy"),
    "metadata_json": ("", "{run}_run_metadata.json", "{run}_run_metadata.json"),
    "metadata_txt":  ("", "{run}_run_metadata.txt", "{run}_run_metadata.txt"),

    # --- figures ----------------------------------------------------------
    # "annotated" never said what the figure was: it is the rate comparison
    # drawn across the buffer domains as well as the real ones.
    "figure_rate":          (FIGURES, "shoreline_change_rate.png",
                             "{run}_shoreline_change_rate_REAL_DOMAINS_ONLY.png"),
    "figure_rate_buffers":  (FIGURES, "shoreline_change_rate_with_buffers.png",
                             "{run}_annotated.png"),

    # --- tables -----------------------------------------------------------
    "rate_csv":      (TABLES, "shoreline_change_rate.csv",
                      "{run}_shoreline_change_rate.csv"),
    "groin_csv":     (TABLES, "groin_diagnostics.csv",
                      "{run}_groin_diagnostics.csv"),
    "road_csv":      (TABLES, "road_management.csv",
                      "road_management_summary.csv"),
    "nourishment_csv": (TABLES, "nourishment_log.csv",
                        "{run}_nourishment_log.csv"),
}

# An animation's name is built from its mode and window rather than listed,
# because the windows are open-ended (any GIS range is legal). These turn the
# plotting module's range tag into the new form:
#   D1-90                          -> domains_1-90
#   groinZoom_BuxtonGroin_D1-15    -> buxton_groin_domains_1-15
#   groinSpan_D1-15                -> groin_span_domains_1-15
#   ALL120pad                      -> all_120_padded
_TAG_RULES = (
    (re.compile(r"^groinZoom_(.+?)_D(\d+)-(\d+)$"), lambda m: f"{_snake(m.group(1))}_domains_{m.group(2)}-{m.group(3)}"),
    (re.compile(r"^groinSpan_D(\d+)-(\d+)$"),       lambda m: f"groin_span_domains_{m.group(1)}-{m.group(2)}"),
    (re.compile(r"^ALL(\d+)pad$"),                  lambda m: f"all_{m.group(1)}_padded"),
    (re.compile(r"^D(\d+)-(\d+)$"),                 lambda m: f"domains_{m.group(1)}-{m.group(2)}"),
)


def _snake(text: str) -> str:
    """BuxtonGroin -> buxton_groin; leaves an already-snake tag alone."""
    s = re.sub(r"(?<=[a-z0-9])(?=[A-Z])", "_", str(text))
    return re.sub(r"[^0-9a-zA-Z]+", "_", s).strip("_").lower()


def animation_name(mode: str, range_tag: str) -> str:
    """New-layout filename for one animation, e.g. 'position_domains_1-90.gif'."""
    tag = str(range_tag)
    for pattern, build in _TAG_RULES:
        m = pattern.match(tag)
        if m:
            tag = build(m)
            break
    else:
        tag = _snake(tag)
    return f"{_snake(mode)}_{tag}.gif"


def legacy_animation_name(run_name: str, mode: str, range_tag: str) -> str:
    """The flat name the same animation had before 2026-09-10."""
    return f"{run_name}_shoreline_{mode}_{range_tag}.gif"


def _names(kind: str, run_name: str):
    try:
        sub, new_tmpl, old_tmpl = KINDS[kind]
    except KeyError:
        raise KeyError(
            f"unknown run-folder file kind {kind!r}; known kinds: "
            f"{', '.join(sorted(KINDS))}, plus animations via animation_name()"
        ) from None
    return sub, new_tmpl.format(run=run_name), old_tmpl.format(run=run_name)


def write_path(run_dir, kind: str, run_name: str, make_parent: bool = True) -> Path:
    """Where to WRITE this file: always the new layout. Parent created."""
    sub, new_name, _ = _names(kind, run_name)
    p = Path(run_dir) / sub / new_name if sub else Path(run_dir) / new_name
    if make_parent:
        p.parent.mkdir(parents=True, exist_ok=True)
    return p


def resolve(run_dir, kind: str, run_name: str, must_exist: bool = False):
    """Where to READ this file: new layout if present, else the old flat one.

    Returns the new-layout path when neither exists, so a caller that is about
    to create the file gets the right place. `must_exist=True` raises instead.
    """
    sub, new_name, old_name = _names(kind, run_name)
    run_dir = Path(run_dir)
    candidates = [run_dir / sub / new_name if sub else run_dir / new_name,
                  run_dir / old_name]
    for c in candidates:
        if c.is_file():
            return c
    if must_exist:
        raise FileNotFoundError(
            f"{kind} not found for {run_name} in either layout:\n  "
            + "\n  ".join(str(c) for c in candidates))
    return candidates[0]


def resolve_animation(run_dir, run_name: str, mode: str, range_tag: str,
                      must_exist: bool = False):
    """`resolve` for an animation, whose name is built rather than listed."""
    run_dir = Path(run_dir)
    candidates = [run_dir / ANIMATIONS / animation_name(mode, range_tag),
                  run_dir / legacy_animation_name(run_name, mode, range_tag)]
    for c in candidates:
        if c.is_file():
            return c
    if must_exist:
        raise FileNotFoundError(
            f"animation {mode}/{range_tag} not found for {run_name}")
    return candidates[0]


def animation_write_path(run_dir, mode: str, range_tag: str,
                         make_parent: bool = True) -> Path:
    p = Path(run_dir) / ANIMATIONS / animation_name(mode, range_tag)
    if make_parent:
        p.parent.mkdir(parents=True, exist_ok=True)
    return p


# =============================================================================
# MIGRATION
# =============================================================================

def plan_run(run_dir) -> list[tuple[Path, Path]]:
    """(source, destination) for every file in this run folder that moves.

    Only files whose OLD name is present and whose NEW path is free are
    planned, so running this twice is a no-op and an interrupted migration
    resumes cleanly.
    """
    run_dir = Path(run_dir)
    hits = sorted(run_dir.glob("*_run_metadata.json"))
    if not hits:
        return []
    run_name = hits[0].name[: -len("_run_metadata.json")]
    moves = []
    for kind, (sub, new_tmpl, old_tmpl) in KINDS.items():
        if not sub:
            continue                      # root files do not move
        src = run_dir / old_tmpl.format(run=run_name)
        dst = run_dir / sub / new_tmpl.format(run=run_name)
        if src.is_file() and not dst.exists():
            moves.append((src, dst))
    # animations: anything matching the flat pattern
    flat = re.compile(re.escape(run_name) + r"_shoreline_(\w+?)_(.+)\.gif$")
    for src in sorted(run_dir.glob("*.gif")):
        m = flat.match(src.name)
        if not m:
            continue
        dst = run_dir / ANIMATIONS / animation_name(m.group(1), m.group(2))
        if not dst.exists():
            moves.append((src, dst))
    return moves


def migrate_run(run_dir, dry_run: bool = False) -> list[tuple[Path, Path]]:
    """Move one run folder into the new layout. Returns what was (to be) moved."""
    moves = plan_run(run_dir)
    if dry_run:
        return moves
    for src, dst in moves:
        dst.parent.mkdir(parents=True, exist_ok=True)
        shutil.move(os.fspath(src), os.fspath(dst))
    return moves
