"""
Modelled against observed overwash for several runs side by side, from each run's own figures/overwash_vs_observed.png
(2026-11-11). Default set: the module off against the median back row with lateral erosion (testing log section 16).

    python scripts/real-world-cases/Hatteras/dune-collision/4-hindcast/3-on-vs-off/HAT_overwash_compare.py

Writes, in the workspace's dune-collision/4-hindcast/3-collision-design/dune-rows/ (D:):
- dune_rows_maxmed_lateral_overwash_vs_off.png: the run figures in a 2 x 2 grid, each headed with its run

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-11-11
"""
from __future__ import annotations

import sys
from pathlib import Path

_HERE = Path(__file__).resolve()
sys.path.insert(0, str(_HERE.parent))
import HAT_collision_on_off_position as P  # noqa: E402
import HAT_collision_row_mode_sweep as R  # noqa: E402

H = R.H

# --- CONFIG ------------------------------------------------------------------
OUT = P.OUT / "3-collision-design" / "dune-rows" / "dune_rows_maxmed_lateral_overwash_vs_off.png"
RUNS = (("cal_step2_growth_noscrape", "Module off"),
        ("cal_toeberm_lat_cs3_maxmed", "Back row at median, lateral, Cs 3e-6"),
        ("cal_toeberm_lat_cs10_maxmed", "Back row at median, lateral, Cs 1e-5"),
        ("cal_toeberm_lat_cs20_maxmed", "Back row at median, lateral, Cs 2e-5"))
NCOL = 2
# -----------------------------------------------------------------------------


def main():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.image as mpimg
    import matplotlib.pyplot as plt
    from site_layer.hat_figure_style import record_caption
    imgs = [(mpimg.imread(H.run_dir(v) / "figures" / "overwash_vs_observed.png"), lab) for v, lab in RUNS]
    h, w = imgs[0][0].shape[:2]
    nrow = -(-len(imgs) // NCOL)
    scale = 9.0 / w                                   # inches per pixel at 9 in per panel
    fig, axs = plt.subplots(nrow, NCOL, figsize=(NCOL * w * scale, nrow * (h * scale + 0.45)))
    for ax, (img, lab) in zip(axs.ravel(), imgs):
        ax.imshow(img)
        ax.set_title(lab, loc="left", fontsize=13, fontweight="bold")
        ax.axis("off")
    for ax in axs.ravel()[len(imgs):]:
        ax.axis("off")
    fig.tight_layout()
    OUT.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT, dpi=200, facecolor="white")
    plt.close(fig)
    record_caption(OUT, (
        "Modelled against observed overwash, 1996-2009, each panel the run's own overwash_vs_observed figure: the "
        "dune collision module off, then on, with the back dune row at the median of each profile's dune window and "
        "lateral erosion (the face retreats through the front row, then the back), toe at the berm, Cs 3e-6, 1e-5 "
        "and 2e-5. All with the road fix, growth 0.02-0.04, groin 0.16, no road rule, Rodanthe 40%."))
    print(f"wrote {OUT}")


if __name__ == "__main__":
    main()
