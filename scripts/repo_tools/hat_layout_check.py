"""
Check the repository against the seven layout rules in ORGANIZATION.md.

    python scripts/repo_tools/hat_layout_check.py            every rule
    python scripts/repo_tools/hat_layout_check.py --rule 5   just one
    python scripts/repo_tools/hat_layout_check.py --full     every offender, not the first few

Advisory: prints a worklist and always exits zero. Details: scripts/repo_tools/README.md.

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-09-22
"""

from __future__ import annotations

import argparse
import re
from pathlib import Path

REPO = next(p for p in Path(__file__).resolve().parents
            if (p / "pyproject.toml").exists())

# --- CONFIG ------------------------------------------------------------------
# Trees the rules skip: vendored code and caches
SKIP_PARTS = {".git", "__pycache__", ".nox", ".venv", "node_modules",
              ".pytest_cache", ".ipynb_checkpoints", "coastal_cascade.egg-info"}

# Rule 1: what counts as data rather than code.
DATA_SUFFIXES = {".csv", ".npy", ".npz", ".tif", ".tiff", ".geojson", ".png",
                 ".pdf", ".gif", ".xlsx", ".shp", ".docx",
                 # .log added 2026-09-22: a log is a product, and products live in output/
                 ".log"}
# Data-shaped files allowed beside code, as its configuration
DATA_ALLOWED = {"reference_yaml_hatteras.yaml"}

# Rule 1, second case: the root of scripts/ holds folders and README.md only
ROOT_ALLOWED = {"README.md"}

# Rule 4: superseded_<date>, optionally with a reason after the date
RETIRE_GOOD = re.compile(r"^superseded_\d{8}(_[A-Za-z0-9][\w-]*)?$")

# output/archive/ files retired material as YYYY-MM-DD_<what>/, its own idiom
ARCHIVE_DATED = re.compile(r"^\d{4}-\d{2}-\d{2}_")
RETIRE_ANY = re.compile(r"^(old_.*|old|.*_ARCHIVE.*|archived_.*|.*_backup|"
                        r"retired.*|superseded.*)$", re.I)

# Rule 5: a counted depth used to find a root, flagged only for root-like names
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
# -----------------------------------------------------------------------------


# Retired folders are exempt from rule 5 (2026-09-14)
RETIRED_PART = re.compile(r"^(superseded.*|old_.*|old|archived_.*|.*_ARCHIVE.*)$",
                          re.I)


# Is any part of the path a retired folder?
def is_retired(path: Path) -> bool:
    return any(RETIRED_PART.match(part) for part in path.parts)


# Every path under root, skipping caches and, if asked, retired folders
def walk(root: Path, skip_retired: bool = False):
    for path in root.rglob("*"):
        if any(part in SKIP_PARTS for part in path.parts):
            continue
        if skip_retired and is_retired(path):
            continue
        yield path


# Extensionless files that are meant to be here
NO_SUFFIX_ALLOWED = {"CURRENT", "LICENSE", "Makefile", "Dockerfile", ".gitignore",
                     ".gitattributes", "MANIFEST.in", "py.typed"}


# Rule 1: data files under scripts/, and files with no extension
def rule_1_data_beside_code():
    out = []
    for path in walk(REPO / "scripts"):
        if not path.is_file():
            continue
        if path.name in DATA_ALLOWED or path.name in NO_SUFFIX_ALLOWED:
            continue
        if path.suffix.lower() in DATA_SUFFIXES:
            out.append(path.relative_to(REPO))
        elif not path.suffix:
            # Untyped: it needs a suffix before any other check can see it
            out.append(path.relative_to(REPO))
    return out


# Every bare name any code in the repo imports
def imported_module_names():
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
            # A notebook's JSON holds the same statements, so one pattern reads both
            for hit in IMPORT_STMT.findall(text):
                names.add(hit.split(".")[0])
    return names


# Rule 1: files at the root of scripts/ other than README.md
def rule_1_loose_files_at_scripts_root():
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


# Rule 2: period folders that are neither a vintage nor a window
def rule_2_period_folder_names():
    out = []
    for path in walk(REPO / "data" / "hatteras_init"):
        if path.is_dir() and PERIOD_PREFIXED.match(path.name):
            out.append(path.relative_to(REPO))
    return out


# Rule 3: year-named period folders with no PROVENANCE.md
def rule_3_provenance_beside_derived():
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


# Is this already filed under a dated superseded folder?
def inside_compliant_retirement(path: Path) -> bool:
    return any(RETIRE_GOOD.match(part) or ARCHIVE_DATED.match(part)
               for part in path.parts[:-1])


# Rule 4: retirement folders that are not superseded_<date>, or lack a note
def rule_4_retirement_idioms():
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


# Rule 5: counted-depth anchors, and paths that cannot resolve anywhere
def rule_5_paths():
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


# Rule 7: folders that need a README of their own
def rule_7_readmes():
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


# Empty folders in the README trees and tests/
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


# Print one rule's findings, up to limit, and return the count
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


# Run: each wanted rule in turn, then the total
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
