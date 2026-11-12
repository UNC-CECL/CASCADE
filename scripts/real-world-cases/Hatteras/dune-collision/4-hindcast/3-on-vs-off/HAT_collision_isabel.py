"""
Hurricane Isabel (Sep 2003) along the island: the dune crest change across it in the lidar and in the model with the
collision module off and on, and the module's own share (on minus off).

    python scripts/real-world-cases/Hatteras/dune-collision/4-hindcast/3-on-vs-off/HAT_collision_isabel.py

The interval is the 2001 survey (Jan-Mar 2001, unclassified lidar, so its crests may include buildings and
vegetation) to the post-Isabel 2004 survey (Apr-Sep 2004); Isabel is its largest storm (5.58 m NAVD88 peak water in
the model's storm series). Runs (HAT_dune_collision.py, nothing rerun): cal_step2_growth (module off) and
cal_step3_module_la (Larson, Cs 1.9e-5), identical otherwise. Crests: each run's dune_crest_vs_lidar.csv
(domain medians); GIS 6-90 (the module is off at Cape Point).

Writes, in the workspace's dune-collision/4-hindcast/7-storm-examples/ (D:): collision_isabel.png, the PDF and
supporting/collision_isabel.csv.

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-11-04
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

_HERE = Path(__file__).resolve()
sys.path.insert(0, str(_HERE.parent))
import HAT_collision_on_off_position as P  # noqa: E402
OUT = P.OUT / "7-storm-examples"   # this script's figures (topic subfolder)

sys.path.insert(0, str(P._HAT_DIR / "dune-collision" / "4-hindcast" / "1-calibration-1996-2009"))

# --- CONFIG ------------------------------------------------------------------
BEFORE, AFTER = "2001", "2004"
OFF, ON = "cal_step2_growth", "cal_step3_module_la"
LIDAR_C, OFF_C, ON_C = "#2166AC", "#7F7F7F", "#6A51A3"
PLACES = (("Buxton", 7, 9), ("Avon", 21, 31), ("Salvo", 67, 70), ("Waves", 71, 75), ("Rodanthe", 76, 83))
# -----------------------------------------------------------------------------


def main():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import HAT_dune_collision as C
    from cascade_pipeline.plotting.dune_crest_comparison import TABLE, TABLES
    from site_layer import hat_figure_style as S

    S.apply_style()
    crest = {}
    for v in (OFF, ON):
        t = pd.read_csv(C.run_dir(v) / TABLES / TABLE)
        t = t[t.gis.between(6, 90)].pivot_table(index="gis", columns="survey", values=["model_median_m",
                                                                                         "survey_median_m"])
        crest[v] = t
    d = pd.DataFrame({
        "lidar_change_m": crest[OFF]["survey_median_m"][int(AFTER)] - crest[OFF]["survey_median_m"][int(BEFORE)],
        "off_change_m": crest[OFF]["model_median_m"][int(AFTER)] - crest[OFF]["model_median_m"][int(BEFORE)],
        "on_change_m": crest[ON]["model_median_m"][int(AFTER)] - crest[ON]["model_median_m"][int(BEFORE)],
    })
    d["module_effect_m"] = d.on_change_m - d.off_change_m
    g = d.index.to_numpy()

    fig, (a, b) = plt.subplots(2, 1, figsize=(S.FIG_W_DOUBLE, 4.6), sharex=True,
                               gridspec_kw=dict(height_ratios=(2, 1.1), hspace=0.12))
    for ax in (a, b):
        for name, lo, hi in PLACES:
            ax.axvspan(lo - 0.5, hi + 0.5, color="0.5", alpha=0.07, lw=0, zorder=0)
        ax.axhline(0, color=S.INK, lw=0.6)
        ax.grid(axis="y", color=S.GRID_C, lw=0.5)
        ax.spines[["top", "right"]].set_visible(False)
    a.plot(g, d.lidar_change_m, "o", color=LIDAR_C, ms=3.2, label="Lidar", zorder=4)
    a.plot(g, d.off_change_m, color=OFF_C, lw=1.5, label="Model, module off")
    a.plot(g, d.on_change_m, color=ON_C, lw=1.5, label="Model, collision module on")
    for name, lo, hi in PLACES:
        a.text((lo + hi) / 2, 2.35, name, ha="center", va="top", fontsize=7.5, color=S.INK_MUTED)
    a.set_ylim(-2.6, 2.5)
    a.set_ylabel("Crest change,\n2001 to 2004 (m)")
    a.legend(loc="lower right", ncol=3, frameon=False, fontsize=7.5)
    a.text(-0.075, 1.0, "(a)", transform=a.transAxes, ha="left", va="top", fontweight="bold")
    b.bar(g, d.module_effect_m, width=0.8, color=ON_C)
    b.set_ylabel("Module's effect\n(on − off, m)")
    b.set_ylim(min(-1.2, d.module_effect_m.min() - 0.1), 0.2)
    b.set_xlabel(S.DOMAIN_AXIS_LABEL)
    b.set_xlim(5.5, 90.5)
    b.set_xticks(range(10, 91, 10))
    b.text(-0.075, 1.0, "(b)", transform=b.transAxes, ha="left", va="top", fontweight="bold")
    out = OUT / "collision_isabel.png"
    S.save(fig, out, bbox_inches="tight", close=True)
    d.round(3).to_csv(OUT / "supporting" / "collision_isabel.csv")
    S.record_caption(out, (
        "Hurricane Isabel (September 2003) along Hatteras Island. (a) Change in the domain-median dune crest from "
        "the January-March 2001 lidar survey to the post-Isabel 2004 survey: lidar (dots; the 2001 survey is "
        "unclassified, so some of its crests include buildings and vegetation) and the model with the dune "
        "collision module off (grey) and on (purple; Larson et al. 2004 per storm, Cs 1.9e-5). Isabel is the "
        "interval's largest storm, 5.58 m NAVD88 peak water in the model. (b) The module's own effect, the "
        "crest change with it on minus off. GIS 6-90."))
    print(f"wrote {out}")
    lowered = d.module_effect_m < -0.1
    print(f"module lowered {lowered.sum()} of {len(d)} domains by >0.1 m; mean effect "
          f"{d.module_effect_m.mean():.2f} m, largest {d.module_effect_m.min():.2f} m (GIS {d.module_effect_m.idxmin()})")
    print("mean change: lidar %.2f, off %.2f, on %.2f" % (d.lidar_change_m.mean(), d.off_change_m.mean(),
                                                          d.on_change_m.mean()))
    print("corr with lidar: off %.2f, on %.2f" % (d.lidar_change_m.corr(d.off_change_m),
                                                  d.lidar_change_m.corr(d.on_change_m)))


if __name__ == "__main__":
    main()
