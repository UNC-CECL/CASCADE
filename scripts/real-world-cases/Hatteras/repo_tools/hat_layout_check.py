"""
Check the Hatteras/ scripts folder against its layout rules on obx-cascade.

    python scripts/real-world-cases/Hatteras/repo_tools/hat_layout_check.py            every rule
    python scripts/real-world-cases/Hatteras/repo_tools/hat_layout_check.py --rule 1   just one
    python scripts/real-world-cases/Hatteras/repo_tools/hat_layout_check.py --full     every offender

The rules, for a scripts-only folder on a shared branch (data and output live
in the workspace on D:, see ../README.md):

  1. scripts and docs only: no data, figures, logs or other products here
  2. the Hatteras/ root holds folders, README.md, STYLE.md, hatteras_paths.py
     and cascade_patch.py, nothing else
  3. retired code is deleted, not parked in old/, superseded_*/ or archive*/
  4. locations come from hatteras_paths.py: no root found by counting parent
     folders, no pyproject.toml search, no drive-rooted or home-directory paths
  5. nothing beside a script: a script-folder path to a non-code file goes
     through workspace_twin()
  6. a folder (two levels down) has a README.md, or its parent explains it

Advisory: prints a worklist and always exits zero. Rewritten 2026-10-23 for
the move to obx-cascade; the old repository's seven rules (ORGANIZATION.md)
described data/, output/ and scripts/ side by side, which no longer exists.

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-10-23
"""

from __future__ import annotations

import argparse
import re
from pathlib import Path

HATTERAS = next(p for p in Path(__file__).resolve().parents
                if (p / "hatteras_paths.py").exists())

# --- CONFIG ------------------------------------------------------------------
# Trees the rules skip: caches
SKIP_PARTS = {"__pycache__", ".pytest_cache", ".ipynb_checkpoints"}

# Rule 1: what may live here -- code, docs and the configuration scripts read
ALLOWED_SUFFIXES = {".py", ".ipynb", ".md", ".yaml", ".yml", ".patch", ".in"}

# Rule 2: the files allowed at the Hatteras/ root
ROOT_ALLOWED = {"README.md", "STYLE.md", "hatteras_paths.py", "cascade_patch.py"}

# Rule 3: retirement folders
RETIRED = re.compile(r"^(old_.*|old|superseded.*|archived?.*|.*_ARCHIVE.*|.*_backup|retired.*)$", re.I)

# Rule 4: a root found by counting folders (root-like names only), the old
# repository's marker, and paths that only resolve on one machine
COUNTED_ROOT = re.compile(
    r"^\s*\w*(?:REPO|ROOT|PROJECT|BASE)\w*\s*=\s*[\w.()_]*parents\[\d+\]", re.I | re.M)
PYPROJECT_FINDER = re.compile(r'^\s*(\w+)\s*=.*["\']pyproject\.toml["\']\)\.exists\(\)', re.M)
PYPROJECT_ALLOWED = {"REPO_DIR"}          # the runner finds the git repository for its provenance stamp
MACHINE_PATH = re.compile(r'["\'][A-Za-z]:[\\/]|[A-Z]:\\+Users\\+', re.I)
MACHINE_PATH_ALLOWED = {"hatteras_paths.py", "hat_gis_root.py"}   # the two settings, each overridable

# Rule 5: a script-folder name or Path(__file__) joined to a literal
BESIDE = re.compile(
    r'(?<![\w.])(HERE|_HERE|SCRIPT_DIR|THIS_DIR|Path\(__file__\)(?:\.resolve\(\))?)'
    r'(?:\.parent|\.parents\[\d+\])*((?:\s*/\s*f?["\'][^"\']+["\'])+)')
CODE_TARGETS = (".py", ".ipynb")
# Files that are the scripts' own configuration, read from beside them
CONFIG_TARGETS = (".yaml", ".yml")

# Rule 6: folders a person opens, and the ones that never need a README
README_DEPTH = 2
README_SKIP = re.compile(r"^(__pycache__|figures?|tables|logs)$")
# -----------------------------------------------------------------------------


# Every path under Hatteras/, skipping caches
def walk():
    for path in HATTERAS.rglob("*"):
        if any(part in SKIP_PARTS for part in path.parts):
            continue
        yield path


# Every .py file under Hatteras/, with its text
def python_files():
    for path in walk():
        if path.is_file() and path.suffix == ".py":
            yield path, path.read_text(encoding="utf-8", errors="ignore")


# Rule 1: files that are not code or docs
def rule_1_scripts_only():
    return [(p.relative_to(HATTERAS), "belongs in the workspace on D:")
            for p in walk() if p.is_file() and p.suffix.lower() not in ALLOWED_SUFFIXES]


# Rule 2: loose files at the Hatteras/ root
def rule_2_root():
    out = []
    for path in sorted(HATTERAS.iterdir()):
        if path.is_file() and path.name not in ROOT_ALLOWED:
            why = ("a module belongs in site_layer/, a tool in repo_tools/"
                   if path.suffix == ".py" else "a folder README, or the workspace")
            out.append((path.name, why))
    return out


# Rule 3: retirement folders
def rule_3_retired():
    return [(p.relative_to(HATTERAS), "retired code: delete it (git keeps it)")
            for p in walk() if p.is_dir() and RETIRED.match(p.name)]


# Rule 4: locations not taken from hatteras_paths.py
def rule_4_paths():
    counted, pyproject, machine = [], [], []
    for path, text in python_files():
        rel = path.relative_to(HATTERAS)
        if rel.parts[0] == "repo_tools":
            continue          # tools that check or template projects name these on purpose
        if COUNTED_ROOT.search(text):
            counted.append(rel)
        if any(name not in PYPROJECT_ALLOWED for name in PYPROJECT_FINDER.findall(text)):
            pyproject.append(rel)
        if path.name not in MACHINE_PATH_ALLOWED and MACHINE_PATH.search(text):
            machine.append(rel)
    return counted, pyproject, machine


# Rule 5: non-code files reached through the script's own folder
def rule_5_beside():
    out = []
    for path, text in python_files():
        lines = text.splitlines()
        # names put on sys.path are code folders, whatever they are built from
        import_dirs = {n for l in lines if "sys.path" in l
                       for n in re.findall(r'\b([A-Z_][A-Z0-9_]*)\b', l)}
        for line in lines:
            if "workspace_twin(" in line or "sys.path" in line:
                continue
            assigned = re.match(r'\s*([A-Z_][A-Z0-9_]*)\s*=', line)
            if assigned and assigned.group(1) in import_dirs:
                continue
            for m in BESIDE.finditer(line):
                parts = re.findall(r'["\']([^"\']+)["\']', m.group(2))
                if parts[-1].endswith(CODE_TARGETS + CONFIG_TARGETS):
                    continue
                out.append((path.relative_to(HATTERAS), f'{m.group(1)} / "{"/".join(parts)}"'))
    return out


# Rule 6: folders with no README of their own or in the parent
def rule_6_readmes():
    out = []
    for path in walk():
        if not path.is_dir() or README_SKIP.match(path.name):
            continue
        if len(path.relative_to(HATTERAS).parts) > README_DEPTH:
            continue
        if not any(c.is_file() for c in path.iterdir()):
            continue
        if (path / "README.md").exists() or (path.parent / "README.md").exists():
            continue
        out.append(path.relative_to(HATTERAS))
    return out


# Empty folders
def empty_directories():
    return [p.relative_to(HATTERAS) for p in walk() if p.is_dir() and not any(p.iterdir())]


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
    parser = argparse.ArgumentParser(description="report departures from the Hatteras/ layout rules")
    parser.add_argument("--rule", type=int, action="append",
                        help="check only these rules; repeatable")
    parser.add_argument("--full", action="store_true",
                        help="list every offender, not the first few")
    args = parser.parse_args()
    limit = 10_000 if args.full else 6
    wanted = set(args.rule) if args.rule else None

    print("=" * 74)
    print(f"LAYOUT CHECK   {HATTERAS}")
    print("rules in this file's header; advisory, exits zero")
    print("=" * 74)

    total = 0
    if wanted is None or 1 in wanted:
        total += show("files that are not code or docs", rule_1_scripts_only(), 1, limit,
                      "this folder holds scripts only")
    if wanted is None or 2 in wanted:
        total += show("loose files at the Hatteras/ root", rule_2_root(), 2, limit)
    if wanted is None or 3 in wanted:
        total += show("retirement folders", rule_3_retired(), 3, limit)
    if wanted is None or 4 in wanted:
        counted, pyproject, machine = rule_4_paths()
        total += show("roots found by counting parent folders", counted, 4, limit,
                      "use hatteras_paths.py")
        total += show("roots found by searching for pyproject.toml", pyproject, 4, limit,
                      "that marker is the shared repository's; use hatteras_paths.py")
        total += show("paths that only resolve on one machine", machine, 4, limit,
                      "drive-rooted or a home directory")
    if wanted is None or 5 in wanted:
        total += show("non-code files reached beside the script", rule_5_beside(), 5, limit,
                      "wrap the folder in workspace_twin()")
    if wanted is None or 6 in wanted:
        total += show("folders with no README.md", rule_6_readmes(), 6, limit)

    empties = empty_directories()
    print(f"\nALSO  empty directories: {len(empties)}")
    for path in empties[:limit]:
        print(f"    {path}")

    print("\n" + "=" * 74)
    print(f"{total} item(s) to look at. Nothing here is an error.")
    print("=" * 74)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
