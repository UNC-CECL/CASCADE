"""
The road-dune rule in plan view: Hatteras Island with, in each GIS domain, the model's dune, the model's road
(NC 12) behind it, and the rebuild thresholds, 150 m behind the dune island-wide and the per-domain gate (198 m
less the domain's beach).

    python scripts/real-world-cases/Hatteras/dune-collision/4-hindcast/3-on-vs-off/HAT_road_rule_plan_view.py

Positions are drawn along each domain's cross-shore line (due east-west, as the model's domains), measured from
the ocean coast of the island outline (map_elements, HAT_island_outline): the dune at the domain's beach width
(shoreline to dune crest, the 2009 start's toe table, beach_to_crest_m), the road at the model's road setback
behind the dune (median over the 2009-2025 run final_rd150u, the roadway manager's _road_setback_TS), the
thresholds at 150 m and at the per-domain gate behind the dune. The map is turned so north points right (ocean
below), as the crest-offset map, and split into three panels.

Writes, in the workspace's dune-collision/4-hindcast/5-road/road-rule/ (D:): road_rule_plan_view.png and
supporting/road_rule_plan_view.csv.

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
OUT = P.OUT / "5-road" / "road-rule"   # this script's figures (topic subfolder)
import HAT_road_rule_threshold as R  # noqa: E402

# --- CONFIG ------------------------------------------------------------------
PERIOD = ("Test 2009–2025", "HAT_collision_test_2009", "final_rd150u", "2009-start")
UNIFORM_M = 150.0
PANELS = ((6, 34), (34, 62), (62, 90))
DUNE_C, ROAD_C, UNIFORM_C, GATE_C = "#F28C28", "#2B2B2B", "#E7298A", "#6A51A3"   # dune orange, no brown
ROAD_IN_C, ROAD_OUT_C, ROAD_BETWEEN_C = "#1B9E77", "#BDBDBD", "#E7298A"
# -----------------------------------------------------------------------------


def geometry():
    """The island and the domain frame (as HAT_collision_motivation_alongshore.strip_layers): north to the right,
    one domain per unit, y = -(easting - mid) / step, so the ocean is below and landward is up."""
    import geopandas as gpd
    from site_layer import hat_map_layers as ml
    from site_layer.hat_observed_rates import DOMAIN_BOXES
    dom = gpd.read_file(DOMAIN_BOXES).to_crs(26918)
    dom["ID"] = dom.domain_id.astype(int)
    dom = dom[(dom.ID >= 1) & (dom.ID <= 90)].sort_values("ID")
    step, n0 = np.polyfit(dom.ID.to_numpy(), dom.geometry.centroid.y.to_numpy(), 1)
    e_mid = float(dom.geometry.centroid.x.mean())
    land = gpd.read_file(ml.ISLAND_OUTLINE).to_crs(26918).geometry.union_all()
    return dom, land, step, n0, e_mid


def ocean_coast(land, northing, e_lo, e_hi):
    """Easting of the ocean coast (the island's east edge) along one northing."""
    from shapely.geometry import LineString
    cut = land.intersection(LineString([(e_lo, northing), (e_hi, northing)]))
    xs = [c[0] for g in getattr(cut, "geoms", [cut]) if not g.is_empty for c in g.coords]
    return max(xs) if xs else np.nan


def main():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    from site_layer import hat_topo_version as tv
    from site_layer.hat_figure_style import record_caption

    title, module, variant, start = PERIOD
    run = R.read(module, variant)
    beach = (pd.read_csv(tv.product_dir(start) / "dune-toes" / "dune_toe_cells.csv")
             .groupby("gis").beach_to_crest_m.median())
    dom, land, step, n0, e_mid = geometry()
    x0, _, x1, _ = dom.total_bounds
    rows = []
    for g in run.index:
        north = n0 + step * g
        coast = ocean_coast(land, north, x0 - 3000, x1 + 3000)
        b = float(beach.get(g, beach.median()))
        gate = max(198.0 - b, 0.0)
        setback = float(run.loc[g, "setback_median_m"])
        dune = coast - b
        rows.append(dict(gis=g, coast_e=coast, beach_m=b, gate_m=gate, road_setback_m=setback,
                         dune_e=dune, road_e=dune - setback, uniform_e=dune - UNIFORM_M, gate_e=dune - gate,
                         rebuilds=int(run.loc[g, "rebuilds"])))
    t = pd.DataFrame(rows).set_index("gis")
    to_y = lambda e: -(np.asarray(e, float) - e_mid) / step          # noqa: E731
    lo_thr = np.minimum(t.gate_m, UNIFORM_M)
    hi_thr = np.maximum(t.gate_m, UNIFORM_M)
    t["zone"] = np.where(t.road_setback_m <= lo_thr, "rebuilt (both rules)",
                         np.where(t.road_setback_m >= hi_thr, "left alone (both rules)", "between (rules differ)"))
    t.loc[t.road_setback_m.isna(), "zone"] = "no road setback in the model"
    zone_c = {"rebuilt (both rules)": ROAD_IN_C, "left alone (both rules)": ROAD_OUT_C,
              "between (rules differ)": ROAD_BETWEEN_C, "no road setback in the model": "white"}

    land_rings = [np.c_[(np.asarray(p.exterior.coords)[:, 1] - n0) / step, to_y(np.asarray(p.exterior.coords)[:, 0])]
                  for p in getattr(land, "geoms", [land]) if p.geom_type == "Polygon"]
    fig, axes = plt.subplots(len(PANELS), 1, figsize=(12, 11.5))
    for ax, (a, b) in zip(axes, PANELS):
        s = t.loc[a:b]
        for r in land_rings:
            ax.fill(*r.T, facecolor="white", edgecolor="0.3", lw=0.7, zorder=1)
        for g in s.index:                                         # each domain's cross-shore line, coast to sound
            ax.plot([g, g], [to_y(s.loc[g, "coast_e"]), to_y(s.loc[g, "coast_e"]) + 12], color="0.88", lw=0.4,
                    zorder=1.5)
        x = s.index.to_numpy()
        ax.step(x, to_y(s.gate_e), where="mid", color=GATE_C, lw=1.6, zorder=3)
        ax.plot(x, to_y(s.uniform_e), color=UNIFORM_C, lw=1.6, ls="--", zorder=3)
        ax.plot(x, to_y(s.dune_e), color=DUNE_C, lw=2.2, zorder=3)
        ax.plot(x, to_y(s.road_e), color=ROAD_C, lw=0.9, zorder=3.5)
        ax.scatter(x, to_y(s.road_e), c=[zone_c[z] for z in s.zone], s=28, edgecolor=ROAD_C, linewidth=0.6,
                   zorder=4)
        for g in s.index[s.rebuilds > 0]:
            ax.annotate(f"{s.loc[g, 'rebuilds']}", (g, to_y(s.loc[g, "road_e"])), xytext=(0, 6),
                        textcoords="offset points", ha="center", fontsize=6.5, color=ROAD_IN_C, zorder=5)
        ys = np.r_[to_y(s.coast_e), to_y(s.road_e)]
        ax.set_xlim(a - 0.6, b + 0.6)
        ax.set_ylim(np.nanmin(ys) - 0.8, np.nanmax(ys) + 1.2)
        ax.set_aspect("equal")
        ax.set_xticks(range(a, b + 1, 2))
        ax.tick_params(labelsize=8, left=False, labelleft=False)
        for sp in ("left", "right", "top"):
            ax.spines[sp].set_visible(False)
        ax.set_xlabel("GIS domain (south → north; north is right, the ocean below)", fontsize=8.5)
        for name, lo, hi in R.VILLAGES:
            if hi >= a and lo <= b:
                ax.axvspan(max(lo, a) - 0.5, min(hi, b) + 0.5, color="#EEF2F6", zorder=0)
                ax.text((max(lo, a) + min(hi, b)) / 2, ax.get_ylim()[1] - 0.2, name, ha="center", va="top",
                        fontsize=8, color="#555555")
    h = [Line2D([], [], color=DUNE_C, lw=2.2), Line2D([], [], color=UNIFORM_C, lw=1.6, ls="--"),
         Line2D([], [], color=GATE_C, lw=1.6),
         Line2D([], [], color=ROAD_C, lw=0.9, marker="o", mfc=ROAD_IN_C, ms=6),
         Line2D([], [], color=ROAD_C, lw=0, marker="o", mfc=ROAD_OUT_C, ms=6),
         Line2D([], [], color=ROAD_C, lw=0, marker="o", mfc=ROAD_BETWEEN_C, ms=6)]
    lab = ["Model dune (the beach's width behind the ocean coast)", "150 m behind the dune (one island-wide trigger)",
           "Per-domain gate (198 m less the domain's beach, behind the dune)",
           "Model road (NC 12): inside both thresholds, dune rebuilt when knocked down (number: rebuilds)",
           "Model road: beyond both, dune left alone", "Model road: between the two, the rules differ"]
    fig.legend(h, lab, loc="lower center", ncol=2, frameon=False, fontsize=8.5, bbox_to_anchor=(0.5, -0.01))
    fig.suptitle("The road-dune rule in plan view: the model's road behind the dune, test 2009–2025",
                 fontsize=13, fontweight="bold", y=0.995)
    fig.tight_layout(rect=(0, 0.07, 1, 0.97))
    out = OUT / "road_rule_plan_view.png"
    fig.savefig(out, dpi=300, facecolor="white", bbox_inches="tight")
    plt.close(fig)
    (OUT / "supporting").mkdir(parents=True, exist_ok=True)
    t.round(1).to_csv(OUT / "supporting" / "road_rule_plan_view.csv")
    record_caption(out, (
        "Hatteras Island in plan view (north to the right, the ocean below; one grid step is one GIS domain, about "
        "500 m), with, along each domain's east-west line: the model's dune (orange), its beach width (shoreline to "
        "dune crest, 2009 start) behind the ocean coast of the island outline; the model's road (dots), its median "
        "setback behind the dune over the 2009-2025 run; and the dune rebuild thresholds, 150 m behind the dune "
        "(pink dashed, one island-wide trigger) and the per-domain gate, 198 m less the domain's beach (purple). "
        "Dots: green, the road inside both thresholds, so its dune is rebuilt when knocked down (number: rebuilds "
        "in 2009-2025); grey, beyond both, left alone; pink, between the two, the only places the rules differ. "
        "The road and dune are the model's, drawn from today's outline, not surveyed positions."))
    print(f"wrote {out}")
    print(t.zone.value_counts().to_string())


if __name__ == "__main__":
    main()
