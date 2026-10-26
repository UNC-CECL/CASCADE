"""
A slide-sized table of the hindcast's dune and wave settings, for a talk.

    python scripts/figure_making/model/dune_settings_table.py

Writes output/figures/talk/4-model-mechanics/dune_and_wave_settings_table.png
(16:9). Values as adopted 2026-10-13; see MODEL_VARIABLES.md for the full list.

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-10-13
"""
from __future__ import annotations

import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

# Hatteras/ holds the scripts; data/ and output/ live in the workspace on D: (hatteras_paths)
import sys as _hat_sys  # noqa: E402
_HAT_DIR = next(_p for _p in Path(__file__).resolve().parents  if (_p / "hatteras_paths.py").exists())
if str(_HAT_DIR) not in _hat_sys.path:
    _hat_sys.path.insert(0, str(_HAT_DIR))
from hatteras_paths import WORKSPACE as PROJECT_ROOT  # noqa: E402
sys.path.insert(0, str(_HAT_DIR))
from site_layer.hat_figure_style import C, INK, INK_MUTED, apply_style, figure_dir, save  # noqa: E402

# --- CONFIG ------------------------------------------------------------------
OUT = figure_dir("talk", "4-model-mechanics") / "dune_and_wave_settings_table.png"
HEADER = ("Parameter", "Value", "Source")
ROWS = (
    ("Maximum dune elevation\n(Dmaxel)", "5.5 m NAVD88\none value, island-wide",
     "Surveyed dune crest height; top of the\ncalibration overwash plateau (4.5–5.5 m)"),
    ("Dune growth rate\n(rmin – rmax)", "0.55 – 0.95 per year\nrandom per 10 m dune cell",
     "CASCADE high-growth setting;\nresults insensitive to it"),
    ("Dune rebuild height\n(design elevation)", "4.3 m NAVD88\n(3.94 m MHW)",
     "NC-12 safe dune crest\n(Velasquez-Montoya et al., 2021)"),
    ("Dune rebuild trigger", "2.0 m NAVD88\n(1.64 m MHW = berm + 0.3 m)",
     "CASCADE's built-in minimum"),
    ("Wave climate\n(sets alongshore diffusivity)", "Hs 2.0 m · Tp 7.5 s\nasymmetry 0.6 · high-angle 0.5",
     "Calibrated wave climate"),
)
FOOTNOTE = ("Same settings for the calibration (1996–2009) and test (2009–2025) runs. "
            "Rebuilds apply to managed runs only.  Berm 1.7 m NAVD88; MHW 0.36 m NAVD88.")
COL_X = (0.03, 0.36, 0.64)        # left edge of each column, figure fraction
# -----------------------------------------------------------------------------


def main():
    apply_style()
    fig = plt.figure(figsize=(13.33, 7.5))
    ax = fig.add_axes((0, 0, 1, 1))
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.set_axis_off()

    ax.text(COL_X[0], 0.93, "Dune and wave settings", fontsize=26, weight="bold", color=INK, va="center")
    top, bottom = 0.84, 0.13
    row_h = (top - bottom) / (len(ROWS) + 0.6)
    ax.add_patch(plt.Rectangle((0.02, top - 0.6 * row_h), 0.96, 0.6 * row_h, color=C["ACCENT"], lw=0))
    for x, h in zip(COL_X, HEADER):
        ax.text(x, top - 0.3 * row_h, h, fontsize=17, weight="bold", color="white", va="center")
    y = top - 0.6 * row_h
    for i, row in enumerate(ROWS):
        if i % 2 == 0:
            ax.add_patch(plt.Rectangle((0.02, y - row_h), 0.96, row_h, color=C["ACCENT_FILL"], alpha=0.35, lw=0))
        for j, (x, cell) in enumerate(zip(COL_X, row)):
            ax.text(x, y - row_h / 2, cell, fontsize=16 if j < 2 else 15, va="center", linespacing=1.3,
                    color=INK, weight="bold" if j == 0 else "normal")
        y -= row_h
    ax.plot([0.02, 0.98], [y, y], color=INK_MUTED, lw=0.8)
    ax.text(COL_X[0], 0.06, FOOTNOTE, fontsize=12, color=INK_MUTED, va="center")

    OUT.parent.mkdir(parents=True, exist_ok=True)
    save(fig, OUT, vector=False, close=True, dpi=200)
    print(OUT)


if __name__ == "__main__":
    main()
