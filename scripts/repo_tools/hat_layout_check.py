# ==============================================================================
# hat_layout_check.py
#
# Does the repository still follow the seven rules in ORGANIZATION.md?
#
# ADVISORY, ALWAYS EXITS ZERO. It produces a worklist, not an obstacle. Every
# tidy-up in this project so far has been someone noticing a case by hand,
# which does not scale and decays between surveys. This is the same survey,
# run on demand.
#
#     python scripts/repo_tools/hat_layout_check.py            every rule
#     python scripts/repo_tools/hat_layout_check.py --rule 5   just one
#     python scripts/repo_tools/hat_layout_check.py --full     every offender, not the
#                                                   first few
#
# WHY IT DOES NOT FAIL THE BUILD
#   Chosen deliberately (Hannah, 2026-09-13): a check that blocks you mid
#   experiment gets disabled, and a disabled check reports nothing. This one
#   is meant to be run when you want a picture.
#
# Author: Hannah A. Henry, UNC CECL
# ==============================================================================

from __future__ import annotations

import argparse
import re
from pathlib import Path

REPO = next(p for p in Path(__file__).resolve().parents
            if (p / "pyproject.toml").exists())

# Trees that hold the project's own work. Vendored code and caches are not
# judged by these rules.
SKIP_PARTS = {".git", "__pycache__", ".nox", ".venv", "node_modules",
              ".pytest_cache", ".ipynb_checkpoints", "coastal_cascade.egg-info"}

# Rule 1: what counts as data rather than code.
DATA_SUFFIXES = {".csv", ".npy", ".npz", ".tif", ".tiff", ".geojson", ".png",
                 ".pdf", ".gif", ".xlsx", ".shp", ".docx"}
# A few data-shaped files belong beside code because they ARE code's input in
# the sense of configuration, or documentation of it.
DATA_ALLOWED = {"reference_yaml_hatteras.yaml"}

# Rule 1, second case: the root of scripts/ holds no files but README.md.
#
# It held eleven on 2026-09-18 -- six site modules, this checker, and four
# hatteras_site_config_prebe_<stamp>.py snapshots of the solved BE field that
# the calibrate step had written beside the config it copied. The snapshots are
# data wearing a .py extension, which is why the suffix test above never saw
# them; they read as stray copies, which is how the equivalent 2026-08-24 one
# came to be discarded, taking the only record of the one-shot solve with it
# and leaving plot_be_zones.py undrawable for three weeks.
#
# The modules then moved into site_layer/ and the checker into repo_tools/, so
# the root is now a list of folders and one README. That is worth holding: a
# reader who opens scripts/ should see the map, not the map plus whatever was
# most recently left lying on it. Anything new at this level belongs in one of
# the folders, or is a snapshot that belongs in the data tree.
ROOT_ALLOWED = {"README.md"}

# Rule 4: the one retirement idiom.
RETIRE_GOOD = re.compile(r"^superseded_\d{8}$")
RETIRE_ANY = re.compile(r"^(old_.*|old|.*_ARCHIVE.*|archived_.*|.*_backup|"
                        r"retired.*|superseded.*)$", re.I)

# Rule 5: a counted depth used to find a root. Only flagged when the name being
# assigned looks like a root, so ordinary parents[] use is not noise.
COUNTED_ANCHOR = re.compile(
    r"^\s*(?:\w+\s*=\s*)?(?:\w*(?:REPO|ROOT|PROJECT|BASE)\w*)\s*=\s*"
    r"[\w.()_]*parents\[\d+\]", re.I | re.M)
BAD_LITERAL = re.compile(r'r"[/\\](?:scripts|data)[/\\]|[A-Z]:\\+Users\\+')

# Rule 7: folders a person opens. Depth 1 and 2 inside these trees.
README_TREES = ("scripts", "data/hatteras_init", "output", "hard-structures")
README_SKIP = re.compile(r"^(superseded_\d{8}|old_.*|__pycache__|figures?|"
                         r"tables|animations|validation|old)$")

# Rule 2: period folders inside the data tree.
VINTAGE = re.compile(r"^\d{4}$")
WINDOW = re.compile(r"^\d{4}_\d{4}$")
PERIOD_PREFIXED = re.compile(r"^(hindcast|period|run)[_-]\d{4}", re.I)


# Retired code is not maintained -- every superseded folder says so -- and a
# rule 5 complaint about a script nobody will run again is noise that makes the
# real ones harder to see. Rules 4 and 7 still apply to these folders; rule 5
# does not (2026-09-14).
RETIRED_PART = re.compile(r"^(superseded.*|old_.*|old|archived_.*|.*_ARCHIVE.*)$",
                          re.I)


def is_retired(path: Path) -> bool:
    return any(RETIRED_PART.match(part) for part in path.parts)


def walk(root: Path, skip_retired: bool = False):
    for path in root.rglob("*"):
        if any(part in SKIP_PARTS for part in path.parts):
            continue
        if skip_retired and is_retired(path):
            continue
        yield path


def rule_1_data_beside_code():
    """Data files under scripts/."""
    out = []
    for path in walk(REPO / "scripts"):
        if path.is_file() and path.suffix.lower() in DATA_SUFFIXES \
                and path.name not in DATA_ALLOWED:
            out.append(path.relative_to(REPO))
    return out


def imported_module_names():
    """Every bare name any code in the repo imports.

    Read once, not once per candidate. Only the code trees are scanned: a .py
    under data/ is itself data -- the BE rate tables are spelled that way --
    and treating one as an importer would let a snapshot vouch for a snapshot.
    """
    names = set()
    for tree in CODE_TREES:
        base = REPO / tree
        if not base.is_dir():
            continue
        for path in walk(base):
            if not (path.is_file() and path.suffix.lower() in IMPORT_SUFFIXES):
                continue
            try:
                text = path.read_text(encoding="utf-8", errors="ignore")
            except OSError:
                continue
            # A notebook is JSON holding the same statements, one per string,
            # so the unanchored pattern reads both without a JSON parse.
            for hit in IMPORT_STMT.findall(text):
                names.add(hit.split(".")[0])
    return names


def rule_1_loose_files_at_scripts_root():
    """Files at the root of scripts/ other than README.md."""
    root = REPO / "scripts"
    out = []
    for path in sorted(root.iterdir()):
        if not path.is_file() or path.name in ROOT_ALLOWED:
            continue
        if path.suffix.lower() == ".py":
            why = "a module belongs in site_layer/, a tool in repo_tools/"
        else:
            why = "not code -- the data tree or a folder README"
        out.append((path.relative_to(REPO), why))
    return out


def rule_2_period_folder_names():
    """Period folders that are neither a vintage nor a window."""
    out = []
    for path in walk(REPO / "data" / "hatteras_init"):
        if path.is_dir() and PERIOD_PREFIXED.match(path.name):
            out.append(path.relative_to(REPO))
    return out


def rule_3_provenance_beside_derived():
    """Period folders holding a file named for a year, with no PROVENANCE.md.

    Only checked where the folder name IS a bare year, which is where the
    period-start naming convention applies.
    """
    out = []
    for path in walk(REPO / "data" / "hatteras_init"):
        if not (path.is_dir() and VINTAGE.match(path.name)):
            continue
        if inside_compliant_retirement(path) or is_retired(path):
            continue          # a retired folder's provenance is in its WHY.md
        if not any(c.suffix for c in path.iterdir() if c.is_file()):
            continue
        if not (path / "PROVENANCE.md").exists():
            out.append(path.relative_to(REPO))
    return out


def inside_compliant_retirement(path: Path) -> bool:
    """Is this already filed under a dated superseded folder?

    What sits INSIDE one keeps its original name on purpose -- that is what it
    was called when it was in use, and renaming it would erase that. So the
    convention applies to the retirement folder, not to its contents.
    """
    return any(RETIRE_GOOD.match(part) for part in path.parts[:-1])


def rule_4_retirement_idioms():
    """Retirement folders that are not superseded_<date>, or lack a note."""
    out = []
    for tree in README_TREES:
        base = REPO / tree
        if not base.is_dir():
            continue
        for path in walk(base):
            if not path.is_dir() or not RETIRE_ANY.match(path.name):
                continue
            if inside_compliant_retirement(path):
                continue
            if not RETIRE_GOOD.match(path.name):
                out.append((path.relative_to(REPO), "not superseded_<date>"))
            elif not any((path / n).exists()
                         for n in ("WHY.md", "README.md")):
                out.append((path.relative_to(REPO), "no WHY.md"))
    return out


def rule_5_paths():
    """Counted-depth anchors, and paths that cannot resolve anywhere."""
    counted, literal = [], []
    for tree in ("scripts", "tests", "hard-structures"):
        base = REPO / tree
        if not base.is_dir():
            continue
        for path in walk(base, skip_retired=True):
            if path.suffix != ".py":
                continue
            try:
                text = path.read_text(encoding="utf-8", errors="ignore")
            except OSError:
                continue
            if COUNTED_ANCHOR.search(text):
                counted.append(path.relative_to(REPO))
            if BAD_LITERAL.search(text):
                literal.append(path.relative_to(REPO))
    return counted, literal


def rule_7_readmes():
    """Folders that need a README of their own.

    ORIENTATION IS INHERITED (2026-09-14). A folder whose PARENT carries a
    README is already explained there, so asking for one in every child turns
    the rule into noise -- one sweep directory alone holds 488 machine-named
    cells, none of which anyone reads a README for. The rule asks at the
    FRONTIER: a folder holding files whose parent explains nothing.

    That makes documenting a tree from the top genuinely finish, rather than
    exposing a new row of demands each time.
    """
    out = []
    for tree in README_TREES:
        base = REPO / tree
        if not base.is_dir():
            continue
        for path in walk(base):
            if not path.is_dir() or README_SKIP.match(path.name):
                continue
            if inside_compliant_retirement(path) or is_retired(path):
                continue
            if len(path.relative_to(base).parts) > 2:
                continue
            if not any(c.is_file() for c in path.iterdir()):
                continue
            if (path / "README.md").exists():
                continue
            if (path.parent / "README.md").exists():
                continue          # the parent explains it
            out.append(path.relative_to(REPO))
    return out


def empty_directories():
    out = []
    for tree in README_TREES + ("tests",):
        base = REPO / tree
        if not base.is_dir():
            continue
        for path in walk(base):
            if path.is_dir() and not any(path.iterdir()):
                out.append(path.relative_to(REPO))
    return out


def show(title, items, rule, limit, note=""):
    print(f"\nRULE {rule}  {title}")
    if not items:
        print("  clean")
        return 0
    print(f"  {len(items)} item(s)" + (f" -- {note}" if note else ""))
    for item in items[:limit]:
        if isinstance(item, tuple):
            print(f"    {item[0]}   ({item[1]})")
        else:
            print(f"    {item}")
    if len(items) > limit:
        print(f"    ... and {len(items) - limit} more (--full to list)")
    return len(items)


def main():
    parser = argparse.ArgumentParser(
        description="report departures from ORGANIZATION.md")
    parser.add_argument("--rule", type=int, action="append",
                        help="check only these rules; repeatable")
    parser.add_argument("--full", action="store_true",
                        help="list every offender, not the first few")
    args = parser.parse_args()
    limit = 10_000 if args.full else 6
    wanted = set(args.rule) if args.rule else None

    print("=" * 74)
    print(f"LAYOUT CHECK   {REPO}")
    print("rules in ORGANIZATION.md; advisory, exits zero")
    print("=" * 74)

    total = 0
    if wanted is None or 1 in wanted:
        total += show("data files under scripts/", rule_1_data_beside_code(),
                      1, limit, "code and data share no tree")
        total += show("loose files at the scripts/ root",
                      rule_1_loose_files_at_scripts_root(), 1, limit,
                      "the root is folders and README.md, nothing else")
    if wanted is None or 2 in wanted:
        total += show("period folders that are neither a vintage nor a window",
                      rule_2_period_folder_names(), 2, limit,
                      "expected <year> or <start>_<end>")
    if wanted is None or 3 in wanted:
        total += show("year-named folders with no PROVENANCE.md",
                      rule_3_provenance_beside_derived(), 3, limit,
                      "which survey is actually behind it")
    if wanted is None or 4 in wanted:
        total += show("retirement folders off the convention",
                      rule_4_retirement_idioms(), 4, limit,
                      "expected superseded_<date>/ with a WHY.md")
    if wanted is None or 5 in wanted:
        counted, literal = rule_5_paths()
        total += show("roots found by counting parent directories",
                      counted, 5, limit, "search upward instead")
        total += show("paths that cannot resolve on any machine",
                      literal, 5, limit, "drive-rooted or a home directory")
    if wanted is None or 7 in wanted:
        total += show("folders with no README.md", rule_7_readmes(), 7, limit)

    empties = empty_directories()
    print(f"\nALSO  empty directories: {len(empties)}")
    for path in empties[:limit]:
        print(f"    {path}")

    print("\n" + "=" * 74)
    print(f"{total} item(s) to look at. Nothing here is an error; "
          f"see ORGANIZATION.md for why each rule exists.")
    print("=" * 74)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
