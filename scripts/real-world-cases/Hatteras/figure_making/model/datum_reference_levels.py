"""
The three vertical references a dune height passes through: NAVD88, MHW and the berm.

    python scripts/figure_making/model/datum_reference_levels.py

Writes output/figures/4-model-mechanics/barrier3d/barrier3d_datums_navd88_mhw_berm.png.
(a) the three zeros close up; (b) a schematic dune profile with the dune cap (Dmaxel)
and a 12 m NAVD88 DEM dune, each given in all three references. Values are the
runner's (HAT_hindcast.py BERM_ELEVATION, MHW_ELEVATION), with
a dune cap 10.5 m above the berm.

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-10-20
"""
from __future__ import annotations

import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

# Hatteras/ holds the scripts; data/ and output/ live in the workspace on D: (hatteras_paths)
import sys as _hat_sys  # noqa: E402
_HAT_DIR = next(_p for _p in Path(__file__).resolve().parents  if (_p / "hatteras_paths.py").exists())
if str(_HAT_DIR) not in _hat_sys.path:
    _hat_sys.path.insert(0, str(_HAT_DIR))
from hatteras_paths import WORKSPACE as PROJECT_ROOT  # noqa: E402
sys.path.insert(0, str(_HAT_DIR))
from site_layer.hat_figure_style import (C, INK, INK_MUTED, _title, apply_style, caption,  # noqa: E402
                                         figsize, figure_dir, open_frame, save)

# --- CONFIG ------------------------------------------------------------------
OUT = figure_dir("mechanics", "barrier3d") / "barrier3d_datums_navd88_mhw_berm.png"
MHW = 0.36          # m NAVD88, Duck NC gauge (NOAA 8651370); runner MHW_ELEVATION
BERM = 1.70         # m NAVD88; runner BERM_ELEVATION
DMAX = 10.5         # m above the berm (the parameter file's Dmaxel 5.5 m NAVD88 is 3.8 m)
DMAXEL = BERM + DMAX  # m NAVD88
DUNE = 12.0         # m NAVD88, the DEM dune in question
C_NAVD = INK_MUTED
C_MHW = C["LATE"]
C_BERM = C["ADDED"]
C_CAP = C["ACCENT"]
SAND = "#ead9b0"
# -----------------------------------------------------------------------------


def three_refs(z: float) -> str:
    """One elevation in NAVD88, MHW and above the berm."""
    return (f"{z:.2f} m NAVD88\n{z - MHW:.2f} m above MHW\n"
            f"{z - BERM:.2f} m above berm")


def level(ax, z, color, ls="-"):
    ax.axhline(z, color=color, lw=1.1, ls=ls, zorder=2)


def span_arrow(ax, x, z0, z1, color, text, side="right"):
    """A double arrow between two levels, its length written beside it."""
    ax.annotate("", xy=(x, z1), xytext=(x, z0),
                arrowprops=dict(arrowstyle="<->", color=color, lw=0.9, shrinkA=0, shrinkB=0))
    dx, ha = (0.012, "left") if side == "right" else (-0.012, "right")
    ax.text(x + dx, (z0 + z1) / 2, text, color=color, ha=ha, va="center", fontsize=8)


def panel_zeros(ax):
    """(a) the three zeros and the distances between them."""
    ax.axhspan(-0.4, MHW, color=C["WATER"], alpha=0.35, lw=0, zorder=0)
    level(ax, 0.0, C_NAVD, "--")
    level(ax, MHW, C_MHW)
    level(ax, BERM, C_BERM)
    for z, name, col in ((0.0, "NAVD88 = 0", C_NAVD), (MHW, "MHW", C_MHW), (BERM, "Berm", C_BERM)):
        ax.text(0.02, z + 0.04, name, color=col, fontsize=8.5, fontweight="bold", va="bottom")
    span_arrow(ax, 0.55, 0.0, MHW, C_MHW, f"{MHW:.2f} m")
    span_arrow(ax, 0.70, MHW, BERM, C_BERM, f"{BERM - MHW:.2f} m\n(berm above MHW)")
    span_arrow(ax, 0.42, 0.0, BERM, C_NAVD, f"{BERM:.2f} m\n(berm in NAVD88)", side="left")
    ax.set_xlim(0, 1)
    ax.set_ylim(-0.3, 2.2)
    ax.set_xticks([])
    ax.set_ylabel("Elevation (m NAVD88)")
    open_frame(ax)
    ax.spines["bottom"].set_visible(False)
    _title(ax, 0, "The three zeros")


def panel_profile(ax):
    """(b) a schematic profile: D and Dmax are both measured up from the berm."""
    x = np.array([0.00, 0.12, 0.30, 0.46, 0.58, 0.70, 1.00])
    z = np.array([-0.3, MHW, BERM, BERM, DUNE, BERM + 0.6, BERM + 0.4])
    ax.axhspan(-0.8, MHW, xmax=0.14, color=C["WATER"], alpha=0.35, lw=0, zorder=0)
    ax.fill_between(x, -0.8, z, color=SAND, lw=0, zorder=1)
    ax.plot(x, z, color=INK, lw=0.9, zorder=3)
    level(ax, 0.0, C_NAVD, "--")
    level(ax, MHW, C_MHW)
    level(ax, BERM, C_BERM)
    level(ax, DMAXEL, C_CAP, "--")
    ax.text(0.99, DMAXEL + 0.15, "Dune cap (Dmaxel)", color=C_CAP, ha="right", va="bottom",
            fontsize=8.5, fontweight="bold")
    ax.text(0.02, BERM + 0.12, "Berm", color=C_BERM, ha="left", va="bottom",
            fontsize=8.5, fontweight="bold")
    ax.text(0.02, MHW + 0.12, "MHW", color=C_MHW, ha="left", va="bottom", fontsize=8.5,
            fontweight="bold")
    ax.text(0.99, -0.12, "NAVD88 = 0", color=C_NAVD, ha="right", va="top", fontsize=8,
            fontweight="bold")
    # the heights Barrier3D's DuneGrowth uses, both from the berm
    span_arrow(ax, 0.22, BERM, DMAXEL, C_CAP, f"Dmax\n{DMAXEL - BERM:.1f} m", side="left")
    span_arrow(ax, 0.76, BERM, DUNE, INK, f"D\n{DUNE - BERM:.1f} m", side="left")
    ax.plot([0.585, 0.76], [DUNE, DUNE], color=INK, lw=0.6, ls=":", zorder=2)
    ax.text(0.555, DUNE - 1.2, "12 m DEM dune", color=INK, ha="right", va="center",
            fontsize=8.5, fontweight="bold")
    ax.text(0.83, 9.6, three_refs(DUNE), color=INK, fontsize=7.5, va="center", linespacing=1.3)
    ax.text(0.99, DMAXEL + 1.3, three_refs(DMAXEL), color=C_CAP, fontsize=7.5, ha="right",
            va="center", linespacing=1.3)
    ax.set_xlim(0, 1)
    ax.set_ylim(-0.8, DMAXEL + 2.4)
    ax.set_xticks([])
    ax.set_xlabel("Cross-shore (schematic)  ocean → sound", color=INK_MUTED)
    ax.set_ylabel("Elevation (m NAVD88)")
    open_frame(ax)
    _title(ax, 1, "Dune heights are measured from the berm")


def main() -> None:
    apply_style()
    fig, (ax_a, ax_b) = plt.subplots(1, 2, figsize=figsize("double", aspect=0.55),
                                     gridspec_kw=dict(width_ratios=(1, 1.6)),
                                     constrained_layout=True)
    panel_zeros(ax_a)
    panel_profile(ax_b)
    caption(fig, (
        f"Vertical references in the Hatteras hindcast. NAVD88 is the DEM's datum; MHW "
        f"({MHW} m NAVD88, Duck gauge NOAA 8651370) is Barrier3D's internal datum, every "
        f"NAVD88 input converted on loading as z/10 - MHW (dam); the berm ({BERM} m NAVD88 = "
        f"{BERM - MHW:.2f} m MHW) is the zero for dune heights. DuneGrowth uses D (dune height "
        f"above the berm) and Dmax = Dmaxel - BermEl, so the MHW offset cancels in both: a "
        f"{DUNE:.0f} m NAVD88 dune is D = {DUNE - BERM:.1f} m and the {DMAXEL} m NAVD88 cap is "
        f"Dmax = {DMAXEL - BERM:.1f} m. Profile schematic, not to cross-shore scale."))
    print(save(fig, OUT, close=True)[0])


if __name__ == "__main__":
    main()
