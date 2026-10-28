"""
Where the experiment modules live: import this and any experiment imports by name.

    python scripts/real-world-cases/Hatteras/hatteras_ms/experiments/experiment_paths.py   # checks after a move

The experiments are filed in theme folders (waves/, dunes/, storms/, ...) but
import each other across them by bare module name. Importing this puts every
theme folder on sys.path once; run directly it reports any module name that
is in two folders, which would make the import depend on folder order.
Details: scripts/real-world-cases/Hatteras/hatteras_ms/experiments/README.md.

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-10-23
"""
import sys
from pathlib import Path

# --- CONFIG ------------------------------------------------------------------
EXPERIMENTS = Path(__file__).resolve().parent

# Every theme folder: a folder here not starting with "_" or "."
THEME_DIRS = sorted(p for p in EXPERIMENTS.iterdir()
                    if p.is_dir() and not p.name.startswith(("_", ".")))
# -----------------------------------------------------------------------------


# Every module name that is in more than one theme folder
def check():
    seen, clashes = {}, []
    for folder in THEME_DIRS:
        for f in folder.glob("*.py"):
            if f.stem in seen:
                clashes.append(f"{f.stem}: in {seen[f.stem].name}/ and {folder.name}/")
            seen.setdefault(f.stem, folder)
    return clashes


# Importing this module does the work: each folder once
for _folder in THEME_DIRS:
    _s = str(_folder)
    if _s not in sys.path:
        sys.path.insert(0, _s)


if __name__ == "__main__":
    problems = check()
    if problems:
        print("experiment_paths: a module name is in two theme folders --")
        for _p in problems:
            print("   ", _p)
        raise SystemExit(1)
    n = sum(len(list(d.glob("*.py"))) for d in THEME_DIRS)
    print(f"experiment_paths: {n} modules in {len(THEME_DIRS)} theme folders under {EXPERIMENTS}")
