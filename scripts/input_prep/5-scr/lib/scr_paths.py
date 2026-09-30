"""
Where the 5-scr modules live: import this and any 5-scr module imports by name.

    python scripts/input_prep/5-scr/lib/scr_paths.py   # checks the table after a move

Puts every folder in MODULE_DIRS on sys.path once, in order; run directly
it reports any module that is not where the table says. Details: scripts/input_prep/5-scr/lib/README.md.

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-09-30
"""

import sys
from pathlib import Path

# --- CONFIG ------------------------------------------------------------------
SCR = Path(__file__).resolve().parents[1]

# Every module another 5-scr script imports by name, its folder, and who imports it
MODULE_DIRS = {
    # module name             folder                                      imported by
    "coastsat_lrr":           SCR / "lib",                                # lrr, 5yr_bins, extension, dsas
    "coastsat_mean_shoreline": SCR / "1-observations" / "mean_shoreline", # on_imagery, mean_shoreline_windows
    "rates_figures":          SCR / "3-rates",                            # lrr/*, total_change, smoothed_lowess7
    "coastsat_lrr_windows":   SCR / "3-rates" / "coastsat" / "lrr",       # rates_figures, 5yr_bins, net_change
    "duneline_endpoint":      SCR / "3-rates" / "duneline",               # coastsat_endpoint
    "coastsat_vs_duneline":   SCR / "4-comparisons" / "shoreline_vs_duneline",
    "total_change_vs_duneline": SCR / "4-comparisons" / "shoreline_vs_duneline",
    "smoothed_lowess7_vs_duneline": SCR / "4-comparisons" / "shoreline_vs_duneline",
}
# -----------------------------------------------------------------------------


# Every module whose file is not where MODULE_DIRS says
def check():
    missing = []
    for name, folder in MODULE_DIRS.items():
        if not (folder / f"{name}.py").exists():
            missing.append(f"{name}: no {name}.py in {folder}")
    return missing


# Importing this module does the work: each folder once, in declaration order
for _folder in dict.fromkeys(MODULE_DIRS.values()):
    _s = str(_folder)
    if _s not in sys.path:
        sys.path.insert(0, _s)


if __name__ == "__main__":
    problems = check()
    if problems:
        print("scr_paths: the module table is out of date --")
        for _p in problems:
            print("   ", _p)
        raise SystemExit(1)
    print(f"scr_paths: all {len(MODULE_DIRS)} modules found under {SCR}")
