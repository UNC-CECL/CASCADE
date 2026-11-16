"""
The dune toe methods side by side on the 1996 transects (2026-11-05): ours (the line method), Automorph (Itzkin et
al. 2020), maximum curvature (Stockdon et al. 2007) and pybeach's four (Beuzen 2019), one profile set, one shoreline.
Each method's script writes its own toes; this one reads them.

    python scripts/real-world-cases/Hatteras/input_prep/9-dune-toes/HAT_toe_methods_compare.py

Reads (D:, 1-transects-1996/): 1-toes/transect_toes_1996.csv (HAT_toe_method_line.py),
3-methods/2-automorph/automorph_transects_1996.csv (HAT_toe_method_automorph.py),
3-methods/3-curvature/curvature_check_transects_1996.csv (HAT_toe_method_curvature.py --check; the 5 m, lower-half
toe) and 3-methods/4-pybeach/pybeach_transects_1996.csv (HAT_toe_method_pybeach.py, its own environment). Per domain
each method's toe is the median of its kept transects, when at least 3.

Writes 3-methods/5-comparison/: toe_methods_transects_1996.csv, toe_methods_domains_1996.csv,
toe_methods_summary_1996.csv (each against ours, and against the model's current toes) and toe_methods_1996.png.

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-11-05
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from HAT_toe_common import (BERM_M, CREST, EXAMPLES, FILL, MHW, MIN_KEPT, M_AUTOMORPH, M_COMPARE,  # noqa: E402
                            M_CURVATURE, M_PYBEACH, PROFILE, TOES, WATER, YEAR)

METHODS = {  # column: (label, colour, marker)
    "ours": ("Ours (line method, crest in the window)", "#1F5A99", "o"),
    "automorph": ("Automorph, code as sent (Itzkin et al. 2020)", "#7B3F99", "s"),
    "automorph_paper": ("Automorph + the paper's checks", "#B48CC8", "s"),
    "automorph_adapted": ("Automorph, adapted for Hatteras", "#3B1566", "P"),
    "curvature": ("Maximum curvature (Stockdon et al. 2007)", "#E07B39", "D"),
    "pb_pd": ("pybeach perpendicular distance", "#6E6E6E", "v"),
    "pb_mc": ("pybeach maximum curvature", "#C9A227", "^"),
    "pb_rr": ("pybeach relative relief", "#2E8B57", "P"),
    "pb_ml": ("pybeach machine learning (barrier island)", "#C2185B", "*"),
}
ALONGSHORE = ("ours", "automorph", "automorph_adapted", "curvature", "pb_ml")   # panel (c); all in the tables


def gather():
    """One row per transect: each method's toe (m NAVD88) and its distance from the shoreline (m)."""
    t = pd.read_csv(TOES / f"transect_toes_{YEAR}.csv")
    c = pd.DataFrame({"transect": t.transect, "gis": t.gis,
                      "ours": t.toe_m_navd88.where(t.kept), "ours_x": t.toe_from_shoreline_m.where(t.kept),
                      "ours_crest": t.crest_m_navd88, "ours_crest_x": t.crest_s_m - t.shoreline_s_m})
    a = pd.read_csv(M_AUTOMORPH / f"automorph_transects_{YEAR}.csv").set_index("transect")
    for col, src in (("automorph", ""), ("automorph_paper", "paper_"), ("automorph_adapted", "adapted_")):
        c[col] = c.transect.map(a[f"{src}toe_m_navd88"])
        c[f"{col}_x"] = c.transect.map(a[f"{src}toe_from_shoreline_m"])
    c["automorph_crest"] = c.transect.map(a.crest_m_navd88)
    c["automorph_crest_x"] = c.transect.map(a.crest_from_shoreline_m)
    k = pd.read_csv(M_CURVATURE / f"curvature_check_transects_{YEAR}.csv").set_index("transect")
    c["curvature"], c["curvature_x"] = c.transect.map(k.curv5_toe), c.transect.map(k.curv5_x)
    pb = M_PYBEACH / f"pybeach_transects_{YEAR}.csv"
    if pb.exists():
        p = pd.read_csv(pb).set_index("transect")
        for m in ("pd", "mc", "rr", "ml"):
            c[f"pb_{m}"], c[f"pb_{m}_x"] = c.transect.map(p[m]), c.transect.map(p[f"{m}_x"])
    else:
        print(f"no {pb.name}: run HAT_toe_method_pybeach.py in its environment; pybeach left out")
    d = pd.read_csv(TOES / f"domain_toes_{YEAR}.csv").set_index("gis")
    return c, d.current_toe_m_navd88


def compare():
    c, current = gather()
    cols = [m for m in METHODS if m in c]
    med = lambda col: c.groupby("gis")[col].median().where(c.groupby("gis")[col].count() >= MIN_KEPT)
    d = pd.DataFrame({m: med(m) for m in cols}).reindex(range(1, 91))
    d.index.name = "gis"
    d["current"] = current
    summ = []
    for m in cols:
        dt = (c[m] - c.ours).dropna()
        dd = (d[m] - d.ours).dropna()
        dc = (d[m] - d.current).dropna()
        summ.append(dict(method=METHODS[m][0], transects_with_toe=int(c[m].notna().sum()),
                         domains_with_median=int(d[m].notna().sum()),
                         transect_within_0p25_pct=100 * (dt.abs() <= 0.25).mean(),
                         transect_within_0p5_pct=100 * (dt.abs() <= 0.5).mean(),
                         domain_mean_diff_m=dd.mean(), domain_rmse_m=float(np.sqrt((dd ** 2).mean())),
                         domains_over_0p5=int((dd.abs() > 0.5).sum()),
                         toes_below_berm=int((c[m] < BERM_M - 1e-6).sum()),
                         vs_current_rmse_m=float(np.sqrt((dc ** 2).mean()))))
    S = pd.DataFrame(summ)
    M_COMPARE.mkdir(parents=True, exist_ok=True)
    c.round(3).to_csv(M_COMPARE / f"toe_methods_transects_{YEAR}.csv", index=False)
    d.round(3).to_csv(M_COMPARE / f"toe_methods_domains_{YEAR}.csv")
    S.round(3).to_csv(M_COMPARE / f"toe_methods_summary_{YEAR}.csv", index=False)
    with pd.option_context("display.width", 220, "display.max_columns", 20):
        print(S.round(2).to_string(index=False))
    dx = (c.automorph_crest_x - c.ours_crest_x).dropna()
    print(f"Automorph crest vs ours: same ridge (within 10 m) {100 * (dx.abs() <= 10).mean():.0f}%, seaward of ours "
          f"> 10 m {int((dx < -10).sum())}, landward > 10 m {int((dx > 10).sum())} (n {len(dx)})")
    figure(c, d, S)


def figure(c, d, S):
    from HAT_toe_method_line import run
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from site_layer import hat_figure_style as fs
    fs.apply_style()
    profiles = run(write=False)[2]
    fig = plt.figure(figsize=(fs.figsize("double")[0], 7.6), layout="constrained")
    gs = fig.add_gridspec(2, 2, height_ratios=(1.0, 1.05))
    ax_ex = [fig.add_subplot(gs[0, 0]), fig.add_subplot(gs[0, 1])]
    ax_al = fig.add_subplot(gs[1, :])
    cols = [m for m in METHODS if m in c]
    ci = c.set_index("transect")

    for k, (ax, tid) in enumerate(zip(ax_ex, EXAMPLES)):
        s, z, p = profiles[tid][:3]
        x = s - s[p["shore"]]
        r = ci.loc[tid]
        xmax = max(140.0, np.nanmax([r.ours_crest_x, r.automorph_crest_x]) + 20)
        v = (x >= -20) & (x <= xmax)
        ax.fill_between(x[v], z[v], -1, color=FILL, lw=0)
        ax.plot(x[v], z[v], color=PROFILE, lw=1.0)
        ax.axhline(MHW, color=WATER, lw=0.6)
        ax.axhline(BERM_M, color=CREST, lw=0.6, ls=(0, (2, 2)))
        ax.plot(r.ours_crest_x, r.ours_crest, "^", ms=8, color=CREST, mec="white", mew=0.5, zorder=5)
        if np.isfinite(r.automorph_crest_x):
            ax.plot(r.automorph_crest_x, r.automorph_crest, "^", ms=8, mfc="white", mec=METHODS["automorph"][1],
                    mew=1.2, zorder=5)
        for j, m in enumerate(cols):
            if np.isfinite(r[m]):
                ax.plot(r[f"{m}_x"], r[m], METHODS[m][2], ms=(8 if m == "pb_ml" else 6.5) - 0.3 * j,
                        color=METHODS[m][1], mec="white", mew=0.5, zorder=6 + j)
        ax.set_xlim(-20, xmax)
        ax.set_ylim(-1, max(7.0, np.nanmax(z[v]) + 0.5))
        ax.set_xlabel("Distance landward of the shoreline (m)")
        if k == 0:
            ax.set_ylabel("Elevation (m NAVD88)")
            ax.text(-18, BERM_M + 0.12, f"Berm ({BERM_M:g} m)", fontsize=7, color=CREST)
        ax.set_title(f"Transect {tid} (GIS {int(r.gis)})", fontsize=8.5, loc="left")
        fs.open_frame(ax)
        ax.text(-0.02, 1.02, "ab"[k], transform=ax.transAxes, fontsize=10, fontweight="bold", ha="right",
                va="bottom")
    h = [plt.Line2D([], [], lw=0, marker="^", ms=7, color=CREST, label="Our crest"),
         plt.Line2D([], [], lw=0, marker="^", ms=7, mfc="white", mec=METHODS["automorph"][1], mew=1.2,
                    label="Automorph crest")]

    g = d.index.to_numpy()
    ax_al.axhline(BERM_M, color=CREST, lw=0.6, ls=(0, (2, 2)))
    for m in [m for m in ALONGSHORE if m in cols]:
        lw, zo = (1.8, 6) if m == "ours" else (0.9, 4)
        ax_al.plot(g, d[m], "-", marker=METHODS[m][2], ms=3, lw=lw, color=METHODS[m][1], zorder=zo,
                   label=METHODS[m][0])
    ax_al.set_xlim(0.5, 90.5)
    ax_al.set_xlabel("GIS domain (south → north)")
    ax_al.set_ylabel("Domain dune toe (m NAVD88)")
    fs.open_frame(ax_al)
    extra = [plt.Line2D([], [], lw=0, marker=METHODS[m][2], ms=5, color=METHODS[m][1], label=METHODS[m][0])
             for m in cols if m not in ALONGSHORE]                   # drawn in (a, b) only
    fig.legend(handles=ax_al.get_legend_handles_labels()[0] + extra + h, fontsize=6.8, ncol=4,
               loc="outside lower center", frameon=False)
    ax_al.text(-0.01, 1.02, "c", transform=ax_al.transAxes, fontsize=10, fontweight="bold", ha="right",
               va="bottom")
    png = M_COMPARE / f"toe_methods_{YEAR}.png"
    fs.save(fig, png, close=True)
    st = S.set_index("method")
    cap = "; ".join(f"{m}: RMSE {st.loc[m].domain_rmse_m:.2f} m, {int(st.loc[m].domains_over_0p5)} domains > 0.5 m"
                    for m in st.index[1:])
    fs.record_caption(png, (
        f"Dune toe methods on the {YEAR} lidar transects, one profile set and one shoreline (mean high water). "
        f"(a, b) Two transects with each method's toe; triangles, our crest (front dune in the model's dune search "
        f"window) and Automorph's (most seaward peak of 2.5-12 m with a 0.6 m drop behind it); dotted, the berm. "
        f"(c) Domain toes, the median of each method's transects (the other methods in (a, b) and the tables). "
        f"Against ours, domain medians: {cap}. Automorph ported from the code of Itzkin et al. (2020); maximum curvature after "
        f"Stockdon et al. (2007), lower half of the rise, 5 m smoothing; pybeach, Beuzen (2019), on our crest and "
        f"shoreline, every toe raised to the berm but Automorph's own."))
    print(f"wrote {png}")


if __name__ == "__main__":
    compare()
