"""
The road-dune rule along Hatteras Island: the model's road (NC 12) behind the dune, and the dunes the model's
road manager rebuilds because the road lies within 150 m of them. Journal style (site_layer.hat_figure_style).

    python scripts/real-world-cases/Hatteras/dune-collision/4-hindcast/3-on-vs-off/HAT_road_rule_strip.py [--trigger 20]

(a) The island outline (the crest-offset map's strip), one column per GIS domain. (b) The island straightened,
each domain along its own cross-shore line, distances measured landward from the dune crest: the beach (shoreline
to dune crest, the 2009 start's toe table, beach_to_crest_m), the 150 m zone, and the model's road at its median
distance behind the dune over the 2009-2025 run final_rd150u (the roadway manager's _road_setback_TS). The dune is
marked where the road lies within the zone. Per-domain thresholds (198 m less each domain's beach) rebuild the same
dunes (HAT_road_rule_threshold.py) and are not drawn.

Writes, in the workspace's dune-collision/4-hindcast/3-on-vs-off/ (D:): road_rule_strip.png, the PDF under
supporting/, and supporting/road_rule_strip.csv.

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
OUT = P.OUT / "4-road-rule"   # this script's figures (topic subfolder)
import HAT_road_rule_threshold as R  # noqa: E402

sys.path.insert(0, str(P._HAT_DIR / "dune-collision" / "1-observations"))

# --- CONFIG ------------------------------------------------------------------
PERIOD = ("HAT_collision_test_2009", "2009-start", "2009–2025")
# The rule drawn: trigger (m behind the dune) -> (the run, the file). 150 m: the 198 m from the shoreline less
# the median beach (the default from 2026-11-04); 20 m: the roadway manager's default before it (hat_run.yaml
# road_dune_rebuild_trigger_m: 20), the same setup without the road rule (step2_module_la, Rodanthe 40%)
RULES = {150.0: ("final_rd150u", "road_rule_strip.png"), 20.0: ("step2_module_la", "road_rule_strip_20m.png")}
Y_LO, Y_HI = -100.0, 650.0           # m landward of the dune crest; roads farther back are drawn at the top edge
GIS_LO, GIS_HI = 6, 90
OCEAN_C, BEACH_C, ZONE_C = "#DCE9F4", "#FFF3C4", "#E3F1E8"     # beach pale yellow, no brown
MAINTAINED_C, NATURAL_C = "#1B7F5A", "#F28C28"                  # maintained green, not maintained orange
PLACES = (("Buxton", 7, 9), ("Avon", 21, 31), ("Salvo", 67, 70), ("Waves", 71, 75), ("Rodanthe", 76, 83))
# -----------------------------------------------------------------------------


def runs_of(mask, x):
    """(start, end) of each run of True along the domains, as domain edges."""
    out, i = [], 0
    while i < len(x):
        if mask[i]:
            j = i
            while j + 1 < len(x) and mask[j + 1]:
                j += 1
            out.append((x[i] - 0.5, x[j] + 0.5))
            i = j + 1
        else:
            i += 1
    return out


def main(TRIGGER_M=150.0):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    from matplotlib.patches import Patch
    from site_layer import hat_topo_version as tv
    from site_layer import hat_figure_style as S
    import HAT_collision_motivation_alongshore as A

    S.apply_style()
    module, start, label = PERIOD
    variant, out_name = RULES[TRIGGER_M]
    run = R.read(module, variant)
    beach = (pd.read_csv(tv.product_dir(start) / "dune-toes" / "dune_toe_cells.csv")
             .groupby("gis").beach_to_crest_m.median())
    g = np.arange(GIS_LO, GIS_HI + 1)
    t = pd.DataFrame(index=g)
    t["beach_m"] = beach.reindex(g).fillna(beach.median())
    t["road_m"] = run.setback_median_m.reindex(g)
    t["maintained"] = t.road_m <= TRIGGER_M
    x = g.astype(float)

    fig = plt.figure(figsize=(S.FIG_W_DOUBLE, 4.6))
    gs = fig.add_gridspec(2, 1, height_ratios=(0.9, 2.2), hspace=0.12)
    ax_m = fig.add_subplot(gs[0])
    ax = fig.add_subplot(gs[1], sharex=ax_m)

    # (a) the island
    A.island_strip(ax_m)
    ax_m.set_xlim(5.5, 90.5)
    # cut the map at panel (b)'s left edge (GIS 6): the strip draws Cape Hatteras and GIS 1-5 into the margin
    from matplotlib.text import Text
    for art in ax_m.get_children():
        if isinstance(art, Text):
            if art.get_transform() == ax_m.transData and art.get_position()[0] < 5.5:
                art.set_visible(False)                         # the cape's name and the GIS 1 label
        elif hasattr(art, "set_clip_on") and art is not ax_m.patch:
            art.set_clip_on(True)
    ax_m.text(0.004, 0.36, "← Cape Hatteras", transform=ax_m.transAxes, ha="left", va="center",
              fontsize=7.5, color=S.INK_MUTED, style="italic")
    ax_m.text(-0.085, 1.0, "(a)", transform=ax_m.transAxes, ha="left", va="top", fontweight="bold", fontsize=10)

    # (b) the island straightened
    xs = np.r_[x[0] - 0.5, x, x[-1] + 0.5]
    shore = -np.r_[t.beach_m.iloc[0], t.beach_m.to_numpy(), t.beach_m.iloc[-1]]
    ax.fill_between(xs, Y_LO, shore, color=OCEAN_C, lw=0, zorder=0)
    ax.fill_between(xs, shore, 0, color=BEACH_C, lw=0, zorder=1)
    ax.plot(xs, shore, color="#8DB3D3", lw=0.6, zorder=2)
    ax.axhspan(0, TRIGGER_M, color=ZONE_C, lw=0, zorder=0.5)
    ax.axhline(TRIGGER_M, color=MAINTAINED_C, lw=0.7, ls=(0, (4, 2.5)), zorder=2)
    ax.axhline(0, color=NATURAL_C, lw=2.4, zorder=3)
    has_road = t.road_m.notna().to_numpy()
    for a, b in runs_of(has_road & t.maintained.to_numpy(), x):
        ax.plot([a, b], [0, 0], color=MAINTAINED_C, lw=4.2, solid_capstyle="butt", zorder=4)
    xr, yr = x[has_road], t.road_m.to_numpy()[has_road]
    ax.plot(xr, np.minimum(yr, Y_HI), color=S.C["ROAD"], lw=1.6, zorder=5)
    far = yr > Y_HI
    if far.any():
        ax.plot(xr[far], np.full(far.sum(), Y_HI), "^", color=S.C["ROAD"], ms=3.2, clip_on=False, zorder=5)
    for name, lo, hi in PLACES:
        ax.axvspan(lo - 0.5, hi + 0.5, color="0.5", alpha=0.06, lw=0, zorder=0.6)
        ax.text((lo + hi) / 2, Y_HI - 15, name, ha="center", va="top", fontsize=7.5, color=S.INK_MUTED)
    ax.text(63.0, TRIGGER_M + 10, f"{TRIGGER_M:g} m", ha="center", va="bottom", fontsize=7.5, color=MAINTAINED_C)
    ax.text(89.8, Y_LO + 12, "Atlantic Ocean", ha="right", va="bottom", fontsize=7.5, style="italic",
            color="#4F7CA6")
    ax.set_xlim(5.5, 90.5)
    ax.set_ylim(Y_LO, Y_HI)
    ax.set_yticks([0, 150, 300, 450, 600])
    ax.set_ylabel("Distance landward of\nthe dune crest (m)")
    ax.set_xticks(range(10, 91, 10))
    ax.set_xlabel(S.DOMAIN_AXIS_LABEL)
    ax.grid(axis="y", color=S.GRID_C, lw=0.5, zorder=0.4)
    ax.spines[["top", "right"]].set_visible(False)
    ax.text(-0.085, 1.02, "(b)", transform=ax.transAxes, ha="left", va="bottom", fontweight="bold", fontsize=10)
    h = [Line2D([], [], color=S.C["ROAD"], lw=1.6), Line2D([], [], color=MAINTAINED_C, lw=4.2),
         Line2D([], [], color=NATURAL_C, lw=2.4), Patch(facecolor=ZONE_C, edgecolor=MAINTAINED_C, lw=0.6, ls="--"),
         Patch(color=BEACH_C)]
    lab = ["NC 12 (model)", "Dune maintained", "Dune not maintained", f"Within {TRIGGER_M:g} m of the dune crest",
           "Beach"]
    ax.legend(h, lab, loc="upper center", bbox_to_anchor=(0.5, -0.17), ncol=5, frameon=False, fontsize=7.5,
              handlelength=1.6, columnspacing=1.2)

    out = OUT / out_name
    S.save(fig, out, bbox_inches="tight", close=True)
    (OUT / "supporting").mkdir(parents=True, exist_ok=True)
    t.round(1).to_csv(OUT / "supporting" / out_name.replace(".png", ".csv"))
    rule = ("The road-dune rule in the model. The threshold is 198 m (650 ft) from the shoreline, the distance "
            "NCDOT/FHWA (2008) list as a monitoring warning sign for NC 12; less the island's median beach width "
            "(about 50 m), the model rebuilds a dune knocked down by storms wherever the road lies within 150 m landward of it."
            if TRIGGER_M == 150 else
            "The road-dune rule before the 198 m criterion: the model's roadway manager rebuilds a dune knocked "
            "down by storms only where the road lies within 20 m landward of it (its default trigger).")
    S.record_caption(out, (
        rule + " (a) Hatteras Island, north to the right, one column per GIS domain (about 0.5 km). (b) The island "
        "straightened, each domain along its own cross-shore line: the beach (shoreline to dune crest, 2009), the "
        f"{TRIGGER_M:g} m zone, and the model's road at its median distance landward of the dune crest over {label} "
        "(triangles: farther than the axis). The dune is drawn dark green where the road lies within the zone and "
        "the dune is maintained."
        + (" Thresholds set per domain (198 m less each domain's beach) maintain the same dunes in the same years."
           if TRIGGER_M == 150 else "")))
    print(f"wrote {out}")
    print(f"maintained: {int((t.maintained & t.road_m.notna()).sum())} of {int(t.road_m.notna().sum())} road domains")


if __name__ == "__main__":
    import argparse
    _p = argparse.ArgumentParser()
    _p.add_argument("--trigger", type=float, default=150.0, choices=sorted(RULES),
                    help="road rule drawn: 150 (198 m from the shoreline, the default) or 20 (the old default)")
    main(_p.parse_args().trigger)
