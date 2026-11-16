"""
The collision designs in one place (2026-11-13): plain names, file tags, folders and runs, so every figure script
names and files them the same way.

A design is two choices: how a storm erodes the dune (No collision / Lower both rows / Lower front row / Face retreat)
and how the two dune rows start in 1996 (Equal rows / Low back row). Figures for a design live in the workspace's
dune-collision/4-hindcast/3-collision-design/<erosion folder>/<start folder>/ and their file names start with the
design's tag (e.g. lowerboth_lowback_cs2e-6_position_1996_2009.png), so a file still reads right outside its folder.

Runs: `run` is the design on the village-migration default (HAT_dune_collision.py DESIGN_GRID, 2026-11-13);
`old_run` the same design before it (village dune lines fixed), which the hotspot-erosion runs were built on.

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-11-13
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import sys

_HERE = Path(__file__).resolve()
sys.path.insert(0, str(_HERE.parent))
from HAT_collision_on_off_position import OUT  # noqa: E402

ROOT = OUT / "3-collision-design"
ALL_DIR = ROOT / "0-all-designs"
TESTS_DIR = ROOT / "_tests"                         # tests before 2026-11-13, village dune lines fixed
TEST_DIRS = {"cs_sweep": TESTS_DIR / "1-faceretreat-vs-lowerboth-cs-sweep",
             "low_back_row": TESTS_DIR / "2-low-back-row-start",
             "search_window": TESTS_DIR / "3-dune-search-window"}

# Erosion: (plain name, tag, folder, colour); ROYGBIV order, blue left for CoastSat
EROSION = {"nocollision": ("No collision", "1-no-collision", "#D7301F"),
           "lowerboth": ("Lower both rows", "2-lower-both-rows", "#F28E2B"),
           "lowerfront": ("Lower front row", "3-lower-front-row", "#1A9850"),
           "faceretreat": ("Face retreat", "4-face-retreat", "#6A3D9A")}
# Start: (plain name, folder, line style)
START = {"equal": ("Equal rows", "1-equal-rows", "-"),
         "lowback": ("Low back row", "2-low-back-row", (0, (5, 2.5)))}


@dataclass(frozen=True)
class Design:
    erosion: str          # EROSION key
    start: str            # START key
    cs: str | None        # "2e-6" etc., None without the module
    run: str              # on the village-migration default
    old_run: str | None   # before it (village dune lines fixed); None if never run

    @property
    def name(self) -> str:
        """Plain name, as in figures and the testing log: 'Lower both rows, low back row, Cs 2e-6'."""
        return ", ".join([EROSION[self.erosion][0], START[self.start][0].lower()]
                         + ([f"Cs {self.cs}"] if self.cs else []))

    @property
    def tag(self) -> str:
        """File-name start: 'lowerboth_lowback_cs2e-6'."""
        return f"{self.erosion}_{self.start}" + (f"_cs{self.cs}" if self.cs else "")

    @property
    def folder(self) -> Path:
        return ROOT / EROSION[self.erosion][1] / START[self.start][1]

    @property
    def colour(self) -> str:
        return EROSION[self.erosion][2]

    @property
    def linestyle(self):
        return START[self.start][2]


# Simple to complex: erosion in order, equal rows before the low back row
DESIGNS = (
    Design("nocollision", "equal", None, "cal_step2_growth_noscrape_vdefault", "cal_step2_growth_noscrape"),
    Design("nocollision", "lowback", None, "cal_step2_growth_noscrape_maxmed_vdefault",
           "cal_step2_growth_noscrape_maxmed"),
    Design("lowerboth", "equal", "2e-6", "cal_toeberm_cs2_vdefault", "cal_toeberm_cs2"),
    Design("lowerboth", "lowback", "2e-6", "cal_toeberm_cs2_maxmed_both_vdefault", None),
    Design("lowerfront", "equal", "2e-6", "cal_toeberm_cs2_front_vdefault", "cal_toeberm_cs2_front"),
    Design("lowerfront", "lowback", "2e-6", "cal_toeberm_cs2_maxmed_vdefault", "cal_toeberm_cs2_maxmed"),
    Design("faceretreat", "equal", "1e-5", "cal_toeberm_lat_cs10_vdefault", "cal_toeberm_lat_cs10"),
    Design("faceretreat", "lowback", "3e-6", "cal_toeberm_lat_cs3_maxmed_vdefault", "cal_toeberm_lat_cs3_maxmed"),
)
BY_OLD_RUN = {d.old_run: d for d in DESIGNS if d.old_run}
BY_RUN = {d.run: d for d in DESIGNS}


def hotspot_dir(old_run: str) -> Path:
    """A design's hotspot-erosion figures (built on its old run)."""
    return BY_OLD_RUN[old_run].folder / "hotspot-erosion"
