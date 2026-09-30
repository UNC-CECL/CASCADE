"""
scr_paths.py
==============================================================================
Where the 5-scr modules live -- the one place that knows.

WHY THIS EXISTS
    Twelve scripts under 5-scr import a module from a SIBLING folder:
    rates_figures for the house drawing, coastsat_vs_duneline for the chainage
    loader, coastsat_lrr for the OLS fit. Until 2026-09-22 each one built that
    folder's path by hand --

        sys.path.insert(0, str(_REPO / "scripts" / "input_prep" / "5-scr"
                               / "coastsat_vs_duneline"))

    -- so the folder layout was written down in twelve places, in five
    different spellings, and the reorganisation on 2026-09-22 would have had to
    edit all twelve. Worse, a stale one does not raise where it is written: the
    insert succeeds against a directory that no longer exists, and the failure
    surfaces forty lines later as `ModuleNotFoundError: coastsat_vs_duneline`,
    naming the module rather than the path that is wrong.

    Rule 6 of ORGANIZATION.md: a location is decided once, in a resolver, and
    everything else asks. This is that resolver for 5-scr's own modules, as
    site_layer/hat_observed_rates.py is for its data.

USAGE  -- the two lines every 5-scr script that imports a sibling carries:

        sys.path.insert(0, str(_REPO / "scripts" / "input_prep" / "5-scr" / "lib"))
        import scr_paths  # noqa: E402,F401

    The import is the whole point: bringing the module in runs the loop at the
    bottom, which puts every module-bearing folder on sys.path. After it, a
    plain `import rates_figures` resolves no matter which folder the importing
    script sits in.

MOVING A MODULE
    Edit its row in MODULE_DIRS below. Nothing else changes. A name whose
    folder has gone missing is reported by check() with the path that is
    wrong, which is the thing the old hand-built inserts could never say.
==============================================================================

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-09-30
"""

import sys
from pathlib import Path

SCR = Path(__file__).resolve().parents[1]

# Every module another 5-scr script imports BY NAME, and the folder holding it.
# The comment is who imports it, so a move can be checked against real callers.
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


def check():
    """Report any module in MODULE_DIRS whose file is not where it is claimed.

    Returns a list of complaint strings, empty when the table is accurate.
    Run as `python lib/scr_paths.py` after moving anything.
    """
    missing = []
    for name, folder in MODULE_DIRS.items():
        if not (folder / f"{name}.py").exists():
            missing.append(f"{name}: no {name}.py in {folder}")
    return missing


# Importing this module is what does the work. dict.fromkeys keeps the folders
# unique and in declaration order; two modules share a folder.
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
