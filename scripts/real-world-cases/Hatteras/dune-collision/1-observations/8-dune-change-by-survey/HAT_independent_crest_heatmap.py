"""
Is the crest heatmap comparing the right feature? The model against the picked lidar crest and against an independent
one, island-wide (2026-11-12).

    python scripts/real-world-cases/Hatteras/dune-collision/1-observations/8-dune-change-by-survey/HAT_independent_crest_heatmap.py [plot]

The picked crest (survey_dune_crests.csv, the heatmap's) is the highest point in each 10 m cell's hand-picked dune
window, moved landward with the dune. The independent crest (as HAT_dune_change_by_survey.independent_crest, there for
three regions) needs no window: on each 1 m alongshore line of the coast unrolled along NC-12, the highest point from
SEAWARD_OF_ROAD_M seaward of the 2008 road to the first ground below mean high water. Here it runs over GIS 9-89 (where
NC-12 is) in blocks of BLOCK domains, and each domain takes the median over its lines. The strips also remove each
survey's road residual (hat_lidar_strips.surveys); it is added back, so both crests sit on the picked crests' datum
(the model's vertical offsets alone). In the villages the independent crest can be a rooftop (the houses stand between
the road and the dune at Avon and in the Tri-Village), so they are hatched.

Writes, in the workspace's dune-collision/1-observations/8-dune-change-by-survey/ (D:):
- supporting/independent_crest_by_domain.csv: per domain and survey, the picked and the independent crest
- crest_offset_picked_vs_independent.png: (a) model - picked, (b) model - independent, (c) picked - independent

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-11-12
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

_HERE = Path(__file__).resolve()
sys.path.insert(0, str(_HERE.parent))
import HAT_dune_change_by_survey as D  # noqa: E402

L = D.L
OUT = D.OUT
# --- CONFIG ------------------------------------------------------------------
GIS_LO, GIS_HI = 9, 89                 # where NC-12 runs (no road on GIS 1-8; 90 the locked end)
BLOCK = 6                              # domains per strip
SEA_M = 320.0                          # strip reach seaward of the road: past the widest dune (Avon, ~250 m)
MODEL_RUN = "cal_step2_growth_noscrape"            # no collision, equal rows: the presentation heatmap's
MODEL_LABEL = "No collision, equal rows"
VILLAGES = {"Buxton": (7, 8), "Avon": (21, 31), "Tri-Village": (68, 83)}
# -----------------------------------------------------------------------------


def measure():
    picked = pd.read_csv(D.CREST_TABLE)
    picked = picked[picked.status == "ok"]
    rows = []
    for lo in range(GIS_LO, GIS_HI + 1, BLOCK):
        hi = min(lo + BLOCK - 1, GIS_HI)
        r = L.Region(f"block_{lo}_{hi}", f"GIS {lo}-{hi}", range(lo, hi + 1), (), sea_m=SEA_M)
        surv, g, dom = L.surveys(r)
        roads = L.road_lines(g[2])
        st = L.Strip(r, g, roads, dom)
        indep = D.independent_crest(st, surv)
        for gis, (a, b) in sorted(st.spans.items()):
            if not lo <= gis <= hi:
                continue
            m = (st.s >= a) & (st.s <= b)
            for y, (zc, dist) in indep.items():
                resid = -surv[y][2]["road_residual_removed_m"]      # the strips took it off; back on
                p = picked[(picked.gis == gis) & (picked.survey == y)]
                ok = np.isfinite(zc[m])
                rows.append(dict(gis=gis, survey=y,
                                 picked_crest_m=float(p.crest_m_navd88.median()) if len(p) else np.nan,
                                 independent_crest_m=float(np.nanmedian(zc[m]) + resid) if ok.mean() > 0.2 else np.nan,
                                 independent_seaward_of_road_m=float(np.nanmedian(dist[m])) if ok.any() else np.nan,
                                 lines=int(ok.sum())))
        print(f"GIS {lo}-{hi} done", flush=True)
    t = pd.DataFrame(rows)
    (OUT / "supporting").mkdir(parents=True, exist_ok=True)
    t.round(3).to_csv(OUT / "supporting" / "independent_crest_by_domain.csv", index=False)
    return t


def figure(t):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.patches import Rectangle
    from site_layer import hat_figure_style as S
    sys.path.insert(0, str(D._HAT_DIR / "dune-collision" / "4-hindcast" / "3-on-vs-off"))
    import HAT_collision_on_off_position as P
    import HAT_be_hotspot_zones as Z
    from cascade_pipeline.plotting.dune_crest_comparison import TABLE, TABLES
    S.apply_style()
    model = pd.read_csv(Z.H.run_dir(MODEL_RUN) / TABLES / TABLE)[["gis", "survey", "model_median_m",
                                                                  "survey_median_m"]]
    m = t.merge(model, on=["gis", "survey"], how="left")
    surveys = sorted(m.survey.unique())
    cmap, norm, bins = P.crest_offset_style()
    panels = (("model_median_m", "survey_median_m", f"(a) {MODEL_LABEL} − picked crest (the heatmap as used)"),
              ("model_median_m", "independent_crest_m", f"(b) {MODEL_LABEL} − independent crest (no dune window)"),
              ("survey_median_m", "independent_crest_m", "(c) Picked − independent crest (where the picks differ)"))
    fig = plt.figure(figsize=(S.FIG_W_DOUBLE, 7.0), layout="constrained")
    gs = fig.add_gridspec(5, 1, height_ratios=(0.26, 1, 1, 1, 0.07))
    axs = []
    for k, (a, b, title) in enumerate(panels):
        ax = fig.add_subplot(gs[k + 1], sharex=axs[0] if axs else None)
        axs.append(ax)
        d = (m[a] - m[b]).where(m.gis.between(GIS_LO, GIS_HI))
        grid = m.assign(d=d).pivot_table(index="survey", columns="gis", values="d").reindex(
            index=surveys, columns=range(1, 91))
        im = ax.imshow(np.ma.masked_invalid(grid.values), aspect="auto", cmap=cmap, norm=norm, origin="lower",
                       extent=(0.5, 90.5, -0.5, len(surveys) - 0.5), interpolation="nearest")
        if k > 0:
            for lo, hi in VILLAGES.values():
                ax.add_patch(Rectangle((lo - 0.5, -0.5), hi - lo + 1, len(surveys), facecolor="none",
                                       edgecolor="0.35", hatch="////", lw=0.0, zorder=4))
        ok = d.dropna()
        out_v = m.loc[ok.index, "gis"].map(lambda g: not any(lo <= g <= hi for lo, hi in VILLAGES.values()))
        ax.set_title(title, loc="left", fontsize=8.5, fontweight="bold")
        ax.set_title(f"outside the villages: mean {ok[out_v].mean():+.2f} m, RMSE {np.sqrt((ok[out_v] ** 2).mean()):.2f} m",
                     loc="right", fontsize=7.5, color=S.INK_MUTED)
        ax.set_yticks(range(len(surveys)), [str(s) for s in surveys], fontsize=6.5)
        ax.set_ylabel("Survey", fontsize=8)
        ax.tick_params(axis="x", labelbottom=(k == 2), labelsize=7.5)
    axs[0].set_xlim(0.5, 90.5)
    axs[0].set_xticks([1] + list(range(10, 91, 10)))
    axs[-1].set_xlabel(S.DOMAIN_AXIS_LABEL)
    P._context_strip(fig.add_subplot(gs[0], sharex=axs[0]), S)
    cax = fig.add_subplot(gs[4])
    cb = fig.colorbar(im, cax=cax, orientation="horizontal", extend="both", ticks=bins)
    cb.set_label("Crest difference, domain median (m)   ← first lower   |   first higher →", fontsize=8)
    cb.ax.tick_params(labelsize=7.5)
    out = OUT / "crest_offset_picked_vs_independent.png"
    S.save(fig, out, dpi=300, facecolor="white", close=True)
    S.record_caption(out, (
        f"Is the crest heatmap comparing the right feature? The model ({MODEL_LABEL.lower()}, {MODEL_RUN}) and the "
        "lidar, domain medians per survey, GIS 9-89 (where NC-12 runs). (a) Model minus the picked crest (the highest "
        "point in each 10 m cell's hand-picked dune window, moved landward with the dune; the heatmap as used). (b) "
        "Model minus an independent crest: on each 1 m line of the coast unrolled along NC-12, the highest point from "
        f"{D.SEAWARD_OF_ROAD_M:g} m seaward of the 2008 road to the first ground below mean high water, no window; the "
        "domain median of its lines. (c) Picked minus independent: red where the pick is below the highest dune "
        "between road and beach (a pick on a lower ridge), blue where above (the pick outside that span). Both lidar "
        "crests on the picked crests' datum (the model's vertical offsets). Hatched, the villages (Buxton GIS 7-8, "
        "Avon 21-31, Salvo to Rodanthe 68-83), where the independent crest can be a rooftop. Top right, the mean and "
        "RMSE outside the villages. Above, the villages and the shoals."))
    print(f"wrote {out}")
    out_v = ~m.gis.map(lambda g: any(lo <= g <= hi for lo, hi in VILLAGES.values()))
    dd = (m.survey_median_m - m.independent_crest_m)[out_v]
    print("picked - independent outside villages, by survey:")
    print(m[out_v].assign(d=dd).groupby("survey").d.agg(["mean", "median", "count"]).round(2).to_string())
    for lo, hi, name in ((9, 15, "Buxton 9-15"), (82, 89, "Rodanthe 82-89")):
        s = m[m.gis.between(lo, hi)].assign(d=m.survey_median_m - m.independent_crest_m)
        print(name, s.groupby("survey").d.mean().round(2).to_dict())


def main():
    t = pd.read_csv(OUT / "supporting" / "independent_crest_by_domain.csv") if "plot" in sys.argv[1:] else measure()
    figure(t)


if __name__ == "__main__":
    main()
