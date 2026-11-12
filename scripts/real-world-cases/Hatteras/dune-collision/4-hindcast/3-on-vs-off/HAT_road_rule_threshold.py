"""
Where the 198 m road-dune rule switches on along the island: the road's distance behind the dune in each
domain against the rebuild threshold, as per-domain gates (198 m less each domain's beach) and as one island-wide
150 m trigger (198 m less the island's median beach, the roadway manager's own road_dune_rebuild_trigger_m).

    python scripts/real-world-cases/Hatteras/dune-collision/4-hindcast/3-on-vs-off/HAT_road_rule_threshold.py

Runs read (nothing rerun): the final setup with the per-domain gates (HAT_dune_collision.py cal_step6_toeavg_la;
HAT_collision_test_2009.py nocape_rd198_dry_toeavg_la) and with one 150 m trigger (cal_final_rd150u; final_rd150u).
Road setback: each roadway manager's yearly distance from the dune to the road (_road_setback_TS), median and range
over the run; rebuilds: the years it rebuilt the dune (_dunes_rebuilt_TS).

Writes, in the workspace's dune-collision/4-hindcast/5-road/road-rule/ (D:): road_rule_threshold.png and
supporting/road_rule_threshold.csv.

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-11-04
"""
from __future__ import annotations

import importlib
import sys
from pathlib import Path

import numpy as np
import pandas as pd

_HERE = Path(__file__).resolve()
sys.path.insert(0, str(_HERE.parent))
import HAT_collision_on_off_position as P  # noqa: E402
OUT = P.OUT / "5-road" / "road-rule"   # this script's figures (topic subfolder)

for _p in (P._HAT_DIR / "dune-collision" / "4-hindcast" / "1-calibration-1996-2009",
           P._HAT_DIR / "dune-collision" / "4-hindcast" / "2-test-2009-2025"):
    sys.path.insert(0, str(_p))

# --- CONFIG ------------------------------------------------------------------
PERIODS = (("Calibration 1996–2009", "HAT_dune_collision", "cal_step6_toeavg_la", "cal_final_rd150u"),
           ("Test 2009–2025", "HAT_collision_test_2009", "nocape_rd198_dry_toeavg_la", "final_rd150u"))
UNIFORM_M = 150.0
GIS = range(6, 91)                     # Cape Point (GIS 1-5) has no collision module and is left out, as the scores
Y_MAX_M = 400.0                        # roads farther back are drawn at the top edge, with their value
ROAD_C, GATE_C, UNIFORM_C, BAND_C, REBUILD_C = "#2B2B2B", "#6A51A3", "#E7298A", "#FDE0EF", "#1B9E77"
VILLAGES = (("Buxton", 7, 8), ("Avon", 21, 31), ("Tri-Village", 68, 83))
# -----------------------------------------------------------------------------


def read(module, variant):
    from cascade_pipeline.saved_runs import load_saved_cascade
    from site_layer.hatteras_site_config import HATTERAS_DOMAINS
    m = importlib.import_module(module)
    rd = m.run_dir(variant)
    c = load_saved_cascade(rd / f"{rd.name}.npz")
    dom = HATTERAS_DOMAINS.for_total(len(c.barrier3d))
    rows = []
    for g in GIS:
        r = c.roadways[dom.gis_to_pad(g)]
        if r is None:
            continue
        s = np.asarray(r._road_setback_TS, float)
        s = s[s > 0]
        rows.append(dict(gis=g, gate_m=r._dune_mgmt_road_setback_trigger,
                         setback_median_m=np.median(s) if s.size else np.nan,
                         setback_min_m=s.min() if s.size else np.nan, setback_max_m=s.max() if s.size else np.nan,
                         rebuilds=int(np.nansum(r._dunes_rebuilt_TS))))
    return pd.DataFrame(rows).set_index("gis")


def main():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    from matplotlib.patches import Patch
    from site_layer.hat_figure_style import record_caption

    tables = []
    fig, axes = plt.subplots(4, 1, figsize=(11, 9.2), sharex=True,
                             gridspec_kw=dict(height_ratios=(3, 0.8, 3, 0.8), hspace=0.22))
    for k, (title, module, gated, uniform) in enumerate(PERIODS):
        g, u = read(module, gated), read(module, uniform)
        d = g.join(u.rebuilds.rename("rebuilds_150"))
        d["period"] = title
        tables.append(d.reset_index())
        ax, bx = axes[2 * k], axes[2 * k + 1]
        x = d.index.to_numpy()
        lo = np.minimum(d.gate_m, UNIFORM_M)
        hi = np.maximum(d.gate_m, UNIFORM_M)
        between = (d.setback_median_m > lo) & (d.setback_median_m < hi)
        for gx in x[between.to_numpy()]:
            ax.axvspan(gx - 0.5, gx + 0.5, color=BAND_C, zorder=0)
        for name, a, b in VILLAGES:
            ax.axvspan(a - 0.5, b + 0.5, color="#EEF2F6", zorder=-1)
            ax.text((a + b) / 2, 6, name, ha="center", va="bottom", fontsize=8.5, color="#555555")
        ax.step(x, d.gate_m, where="mid", color=GATE_C, lw=2, label="Per-domain gate: 198 m less the domain's beach")
        ax.axhline(UNIFORM_M, color=UNIFORM_C, lw=2, ls="--", label="One island-wide trigger, 150 m")
        med = np.minimum(d.setback_median_m, Y_MAX_M)
        ax.vlines(x, np.minimum(d.setback_min_m, Y_MAX_M), np.minimum(d.setback_max_m, Y_MAX_M), color=ROAD_C,
                  lw=0.9, alpha=0.5)
        ax.plot(x, med, "o", color=ROAD_C, ms=4.5, label="Road behind the dune (median over the run; line, range)",
                zorder=4)
        far = d.setback_median_m[d.setback_median_m > Y_MAX_M]
        if len(far):
            ax.plot(far.index, np.full(len(far), Y_MAX_M), "^", color=ROAD_C, ms=5, clip_on=False, zorder=5)
            ax.text(far.index.min() - 0.8, Y_MAX_M * 0.97, f"roads {far.min():.0f}-{far.max():.0f} m back (off scale)",
                    ha="right", va="top", fontsize=8, color="#555555")
        ax.set_ylim(0, Y_MAX_M)
        ax.set_ylabel("Distance behind\nthe dune (m)")
        ax.set_title(title, loc="left", fontsize=11.5, fontweight="bold")
        ax.grid(axis="y", color="#E6E6E6", lw=0.7)
        ax.spines[["top", "right"]].set_visible(False)
        # rebuilds per domain, the two rules side by side
        bx.bar(x - 0.2, d.rebuilds, width=0.4, color=GATE_C, label="per-domain gates")
        bx.bar(x + 0.2, d.rebuilds_150, width=0.4, color=UNIFORM_C, label="150 m")
        bx.set_ylabel("Dune\nrebuilds", fontsize=9)
        bx.spines[["top", "right"]].set_visible(False)
        bx.set_ylim(0, max(3, d[["rebuilds", "rebuilds_150"]].to_numpy().max() + 1))
        bx.text(0.995, 0.92, f"total: gates {d.rebuilds.sum()}, 150 m {d.rebuilds_150.sum()}",
                transform=bx.transAxes, ha="right", va="top", fontsize=8, color="#555555")
    axes[-1].set_xlabel("GIS domain (south → north)")
    axes[-1].set_xlim(5.5, 90.5)
    h = [Line2D([], [], color=ROAD_C, marker="o", ls="", ms=5), Line2D([], [], color=GATE_C, lw=2),
         Line2D([], [], color=UNIFORM_C, lw=2, ls="--"), Patch(color=BAND_C),
         Patch(color=GATE_C), Patch(color=UNIFORM_C)]
    lab = ["Road behind the dune (median; line, range over the run)",
           "Per-domain gate: 198 m less the domain's beach", "One island-wide trigger, 150 m",
           "Road between the two thresholds (the rules disagree)",
           "Dune rebuilds, per-domain gates", "Dune rebuilds, 150 m trigger"]
    fig.legend(h, lab, loc="lower center", ncol=2, frameon=False, fontsize=9, bbox_to_anchor=(0.5, -0.06))
    fig.suptitle("The road-dune rule along the island: per-domain gates against one 150 m trigger",
                 fontsize=13, fontweight="bold", y=0.955)
    fig.text(0.5, 0.925, "Road below a line: its dune is rebuilt when knocked down. Above: the dune is left alone.",
             ha="center", fontsize=10, color="#444444")
    out = OUT / "road_rule_threshold.png"
    fig.savefig(out, dpi=300, facecolor="white", bbox_inches="tight")
    plt.close(fig)
    t = pd.concat(tables, ignore_index=True)
    (OUT / "supporting").mkdir(parents=True, exist_ok=True)
    t.round(1).to_csv(OUT / "supporting" / "road_rule_threshold.csv", index=False)
    record_caption(out, (
        "The rule's 198 m (650 ft) from the shoreline is the distance NCDOT/FHWA (2008) list as a monitoring "
        "warning sign for NC 12, not an NCDOT dune maintenance rule. In the "
        "model the road manager rebuilds a knocked-down dune only while the road is closer behind the dune than its "
        "threshold. Dots: the road's distance behind the dune in each GIS domain, median over the run, with its range; "
        "roads farther than 400 m are drawn at the top with their value. Purple: the per-domain gate, 198 m less the "
        "domain's median beach width; pink dashed: one island-wide trigger of 150 m (198 m less the island's median "
        "beach of about 50 m). Pink bands: domains whose road lies between the two, the only places the rules can "
        "differ. Bars: the years each domain's dune was rebuilt under each rule. Final setup (Larson alone, Cs 1.9e-5, "
        "one toe per domain, Rodanthe 90%); GIS 6-90."))
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
