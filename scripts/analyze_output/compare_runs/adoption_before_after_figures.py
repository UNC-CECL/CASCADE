r"""
adoption_before_after_figures.py -- figures of the 2026-09-28 adoption, before vs after
======================================================================================
Hannah, 2026-09-28: "make figures of the before and after comparison".

Reads the tables adoption_before_after.py writes (scores_, cells_, crest_
<side>.csv, crest_lidar_2009.csv) and the runs' own shoreline_change_rate.csv.
edgeBE only: zeroBE is within 0.1 of it on every score (README beside the tables).
The relocation arms are left out: in both windows they score the same as their
non-relocation twins.

    adoption_scorecard.png                     every score, before -> after, per scenario
    adoption_shoreline_alongshore              model LRR vs CoastSat LOESS-7, managed + natural
    adoption_overwash_map_<scenario>           image x domain: hit / miss / false alarm
    adoption_overwash_by_image                 domains overwashed per image, grouped bars
    adoption_dune_crest_2010                   1996-2010 runs' 2010 crest vs the 2009 lidar

WHERE: output/comparisons/adoption_2026-09-28/figures/
======================================================================================
"""
from __future__ import annotations

import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

_HERE = Path(__file__).resolve()
PROJECT_ROOT = next(_p for _p in _HERE.parents if (_p / "pyproject.toml").exists())
sys.path.insert(0, str(PROJECT_ROOT / "scripts"))
sys.path.insert(0, str(_HERE.parent))
from site_layer.hat_figure_style import (C, C_1997, INK_MUTED, DOMAIN_AXIS_LABEL, apply_style,  # noqa: E402
                                         figsize, _title, open_frame, town_bands, caption, save)
import adoption_before_after as AB  # noqa: E402

TABLES = AB.OUT
FIG_DIR = TABLES / "figures"
PRESET = "edgeBE"
SIDES = {"before": ("before (Dmaxel default, 72 h storms)", C["BASE"]),
         "after": ("after (per-cell ceilings, overwash fixes, 24 h storms)", C["ACCENT"])}
OBS_C = C_1997
SCENARIOS = {  # run-name core -> label, in plot order
    "road_bdm": "full management",
    "road_bdm_nourish": "full management",
    "road_bdm_nonourish": "full management, no fill",
    "road_nobdm": "roadway only",
    "noroad_bdm": "beach/dune only",
    "noroad_bdm_nourish": "beach/dune only",
    "noroad_nobdm": "natural",
}
WLABEL = {"1996_2010": "1996–2010", "2010_2024": "2010–2024"}


def core(run):
    return run.split("offsetmetres_", 1)[1].rsplit("_nogroin", 1)[0]


def load(kind):
    frames = []
    for side in SIDES:
        d = pd.read_csv(TABLES / f"{kind}_{side}.csv")
        d["side"] = side
        frames.append(d)
    d = pd.concat(frames)
    d = d[d.preset == PRESET].copy()
    d["core"] = d.run.map(core)
    return d[d.core.isin(SCENARIOS)]


def run_of(d, window, label):
    cores = [k for k, v in SCENARIOS.items() if v == label]
    return d[(d.window == window) & d.core.isin(cores)]


def side_legend(fig, extra=(), points=False, ncol=None):
    h = [Line2D([], [], color=c, lw=0 if points else 1.6, marker="o" if points else None, ms=5, label=l)
         for l, c in SIDES.values()]
    fig.legend(list(extra) + h, [x.get_label() for x in list(extra) + h],
               loc="outside lower center", ncol=ncol or len(h) + len(extra), frameon=False)


# --- 1. scorecard -----------------------------------------------------------------

METRICS = [("rmse_interior_m_yr", "Shoreline RMSE", "m/yr  ·  lower is better", "{:.2f}"),
           ("bias_interior_m_yr", "Shoreline bias", "m/yr  ·  0 is best", "{:+.2f}"),
           ("PSS", "Overwash skill", "POD − POFD  ·  higher is better", "{:.2f}"),
           ("POD", "Overwash hit rate", "POD  ·  higher is better", "{:.2f}"),
           ("POFD", "Overwash false-alarm rate", "POFD  ·  lower is better", "{:.2f}"),
           ("crest_end_median_m", "Dune crest at end of run", "median, m above MHW", "{:.1f}")]
GAP = 1.1   # blank rows between the two window groups


def scorecard():
    d = load("scores")
    rows = []
    for w in WLABEL:
        for lab in dict.fromkeys(SCENARIOS.values()):
            sub = run_of(d, w, lab)
            if len(sub):
                rows.append((w, lab, sub))
    y, yy = [], 0.0
    for i, (w, _, _) in enumerate(rows):
        if i and w != rows[i - 1][0]:
            yy -= GAP
        y.append(yy)
        yy -= 1
    y = np.array(y)
    fig, axes = plt.subplots(3, 2, figsize=figsize("double", height=9.0), sharey=True,
                             constrained_layout=True)
    fig.get_layout_engine().set(hspace=0.07, wspace=0.05)
    for i, (ax, (col, name, hint, fmt)) in enumerate(zip(axes.flat, METRICS)):
        vals = []
        for yv, (w, lab, sub) in zip(y, rows):
            v = {s: sub[sub.side == s][col].iloc[0] for s in SIDES}
            vals += list(v.values())
            ax.annotate("", xy=(v["after"], yv), xytext=(v["before"], yv),
                        arrowprops=dict(arrowstyle="-", color="0.78", lw=2.2), zorder=1)
            for s, (_, c) in SIDES.items():
                ax.plot(v[s], yv, "o", ms=6.5, color=c, mec="white", mew=0.8, zorder=3)
            left = v["after"] < v["before"]
            ax.annotate(fmt.format(v["after"]), (v["after"], yv), xytext=(-7 if left else 7, 0),
                        textcoords="offset points", color=C["ACCENT"], fontsize=7, zorder=4,
                        ha="right" if left else "left", va="center")
        lo, hi = min(vals), max(vals)
        pad = 0.18 * (hi - lo)
        ax.set_xlim(lo - pad, hi + pad)
        if col == "bias_interior_m_yr":
            ax.axvline(0, color=INK_MUTED, lw=0.7, ls=":", zorder=0)
        if col == "crest_end_median_m":
            lid = pd.read_csv(TABLES / "crest_lidar_2009.csv").crest_2009_lidar_m_mhw.median()
            ax.axvline(lid, color=C["REF"], lw=1.0, ls="--", zorder=0)
            ax.text(lid, y.max() + 0.75, "2009 lidar ", color=C["REF"], fontsize=7, ha="right", va="center")
        _title(ax, i, name)
        ax.set_xlabel(hint, color=INK_MUTED, fontsize=7.5)
        open_frame(ax)
        ax.spines["left"].set_visible(False)
        ax.tick_params(axis="y", length=0)
        ax.grid(axis="x", color=C["GRID"], lw=0.4)
        for w in WLABEL:
            ys = [yv for yv, r in zip(y, rows) if r[0] == w]
            ax.axhspan(min(ys) - 0.5, max(ys) + 0.5, color="0.975" if w == "1996_2010" else "white", zorder=-1)
    ax0 = axes[0, 0]
    ax0.set_ylim(y.min() - 0.7, y.max() + 1.1)
    for a in axes[:, 0]:
        a.set_yticks(y, [lab for _, lab, _ in rows])
        for w in WLABEL:
            ys = [yv for yv, r in zip(y, rows) if r[0] == w]
            a.text(-0.52, np.mean(ys), WLABEL[w], transform=a.get_yaxis_transform(), rotation=90,
                   ha="center", va="center", fontweight="bold", fontsize=8.5)
    side_legend(fig, points=True)
    caption(fig, "Every score of the edgeBE matrix before (grey) and after (purple, value printed) the "
                 "2026-09-28 adoption, one row per scenario, grouped by window. The relocation arms score "
                 "as their non-relocation twins and are left out; zeroBE is within 0.1 on every score. "
                 "(a, b) Interior (GIS 2-89) model LRR minus the CoastSat LRR smoothed over 7 domains. "
                 "(c-e) Each assessed image x domain in the window: model overwash counted from the storms "
                 "that ended between the previous image and this one (less 7 days), against presence in "
                 "the imagery. POD = share of observed overwash the model reproduced; POFD = share of "
                 "observed no-overwash where the model overwashed. (f) Median over GIS 1-90 of each "
                 "domain's median dune crest in the final year; dashed, the median crest of the 2009 "
                 "lidar (the 2010 start).")
    save(fig, FIG_DIR / "adoption_scorecard.png", close=True)


# --- 2. shoreline alongshore --------------------------------------------------------

PAIR = ("full management", "natural")


def domain_axis(ax, xlabel=True):
    ax.set_xlim(0.5, 90.5)
    ax.set_xticks([1, 10, 20, 30, 40, 50, 60, 70, 80, 90])
    town_bands(ax, strip=0.07)
    open_frame(ax)
    if xlabel:
        ax.set_xlabel(DOMAIN_AXIS_LABEL)


def shoreline_alongshore():
    d = load("scores")
    fig, axes = plt.subplots(4, 1, figsize=figsize("double", height=9.2), sharex=True,
                             constrained_layout=True)
    fig.get_layout_engine().set(hspace=0.06)
    i = 0
    for w in WLABEL:
        tab = AB.target_table(tuple(int(x) for x in w.split("_"))).set_index("gis_domain")
        tgt = tab.loc[1:90, "target_lrr_m_yr"]
        for lab in PAIR:
            ax = axes[i]
            ax.axhline(0, color=INK_MUTED, lw=0.6, zorder=0)
            ax.plot(tgt.index, tgt.values, color=OBS_C, lw=2.4, zorder=3, solid_capstyle="round")
            sub = run_of(d, w, lab)
            for s, (_, col) in SIDES.items():
                run = sub[sub.side == s].run.iloc[0]
                lrr = pd.read_csv(AB.SIDES[s] / w / PRESET / run / "tables" / "shoreline_change_rate.csv"
                                  ).set_index("gis_domain")["lrr_m_yr"].loc[1:90]
                ax.plot(lrr.index, lrr.values, color=col, lw=1.2, zorder=4 if s == "after" else 2)
            domain_axis(ax, xlabel=(i == 3))
            ax.set_ylabel("m/yr (+ seaward)")
            _title(ax, i, f"{WLABEL[w]}, {lab}")
            ax.grid(axis="y", color=C["GRID"], lw=0.4)
            i += 1
    side_legend(fig, [Line2D([], [], color=OBS_C, lw=2.4, label="CoastSat LRR, 7-domain LOESS")])
    caption(fig, "Model shoreline change rate (OLS slope of the annual shoreline position) alongshore, "
                 "before (grey) and after (purple) the adoption, against the CoastSat target (blue; LRR, "
                 "LOESS over 7 domains, raw for GIS 1-10). edgeBE; the ends at GIS 1 and 90 are solved "
                 "for each side (before 1996 +4.84/+18.25, 2010 +18.87/+24.24; after 1996 +4.35/+19.09, "
                 "2010 +8.00/+21.26 m/yr). Interior RMSE, before -> after: (a) 1.19 -> 1.17, "
                 "(b) 1.14 -> 1.14, (c) 2.31 -> 2.07, (d) 4.18 -> 2.72 m/yr. GIS 1 is Cape Point, "
                 "GIS 90 Pea Island; shaded strips are the villages.")
    save(fig, FIG_DIR / "adoption_shoreline_alongshore.png", close=True)


# --- 3. overwash maps ---------------------------------------------------------------

OUTCOMES = [  # code, label, colour
    (1, "hit: observed and modelled", C["REF"]),
    (2, "miss: observed, not modelled", C_1997),
    (3, "false alarm: modelled, not observed", C["ADDED"]),
    (4, "neither", "0.93"),
]


def outcome_grid(sub):
    sub = sub.copy()
    sub["model"] = sub.model_m3_per_m > 0
    sub["obs"] = sub.observed == 1
    sub["code"] = np.select([sub.obs & sub.model, sub.obs & ~sub.model, ~sub.obs & sub.model], [1, 2, 3], 4)
    g = sub.pivot_table(index="image", columns="gis", values="code", aggfunc="first")
    return g.reindex(columns=range(1, 91)).sort_index()


def overwash_map(label):
    from matplotlib.colors import ListedColormap, BoundaryNorm
    d = load("cells").dropna(subset=["observed"])
    cmap = ListedColormap([c for _, _, c in OUTCOMES])
    cmap.set_bad("white")
    norm = BoundaryNorm([0.5, 1.5, 2.5, 3.5, 4.5], cmap.N)
    heights = [d[d.window == w].image.nunique() for w in WLABEL]
    fig, axes = plt.subplots(2, 2, figsize=figsize("double", height=6.6), sharex=True,
                             gridspec_kw=dict(height_ratios=heights), constrained_layout=True)
    fig.get_layout_engine().set(hspace=0.06, wspace=0.04)
    for r, w in enumerate(WLABEL):
        sub = run_of(d, w, label)
        for c, s in enumerate(SIDES):
            ax = axes[r, c]
            g = outcome_grid(sub[sub.side == s])
            ax.imshow(np.ma.masked_invalid(g.values.astype(float)), cmap=cmap, norm=norm, aspect="auto",
                      interpolation="nearest", extent=(0.5, 90.5, len(g) - 0.5, -0.5))
            ax.hlines(np.arange(len(g) - 1) + 0.5, 0.5, 90.5, color="white", lw=1.2)
            ax.set_yticks(range(len(g)), pd.to_datetime(g.index).strftime("%Y-%m") if c == 0 else [])
            ax.tick_params(axis="y", length=0, labelsize=7)
            ax.set_xticks([1, 10, 20, 30, 40, 50, 60, 70, 80, 90])
            for sp in ax.spines.values():
                sp.set_visible(False)
            _title(ax, 2 * r + c, f"{WLABEL[w]}, {s}")
            if r == 1:
                ax.set_xlabel(DOMAIN_AXIS_LABEL)
    fig.legend([Patch(color=c) for _, _, c in OUTCOMES] + [Patch(facecolor="white", edgecolor="0.7", lw=0.5)],
               [l for _, l, _ in OUTCOMES] + ["not assessed"], loc="outside lower center", ncol=5,
               frameon=False, fontsize=7.5, handlelength=1.2, columnspacing=1.2)
    sc = load("scores")
    stats = []
    for w in WLABEL:
        sub = run_of(sc, w, label)
        for s in SIDES:
            r_ = sub[sub.side == s].iloc[0]
            stats.append(f"{WLABEL[w]} {s}: {r_.hits} hits, {r_.misses} misses, {r_.false_alarms} false "
                         f"alarms (POD {r_.POD:.2f}, POFD {r_.POFD:.2f})")
    caption(fig, f"Overwash in the imagery against the model, {label}, edgeBE: one row per assessed image, "
                 "one column per GIS domain (Cape Point left, Pea Island right). Model overwash is any "
                 "overwash from the storms that ended between the previous image and this one (less 7 "
                 "days). White: domain not assessed in that image. " + "; ".join(stats) + ".")
    save(fig, FIG_DIR / f"adoption_overwash_map_{label.replace(' ', '_')}.png", close=True)


# --- 4. overwash by image -----------------------------------------------------------

def overwash_by_image():
    d = load("cells").dropna(subset=["observed"])
    d["model"] = (d.model_m3_per_m > 0).astype(int)
    fig, axes = plt.subplots(4, 1, figsize=figsize("double", height=9.0), constrained_layout=True)
    fig.get_layout_engine().set(hspace=0.08)
    i = 0
    width = 0.27
    for w in WLABEL:
        for lab in PAIR:
            ax = axes[i]
            sub = run_of(d, w, lab)
            per = sub.groupby(["side", "image"]).agg(obs=("observed", "sum"), model=("model", "sum"),
                                                     n=("gis", "count")).reset_index()
            a = per[per.side == "after"].sort_values("image")
            x = np.arange(len(a))
            ax.bar(x - width, a.obs, width, color=OBS_C, label="observed in imagery")
            for k, (s, (_, col)) in enumerate(SIDES.items()):
                p = per[per.side == s].sort_values("image")
                ax.bar(x + k * width, p.model.values, width, color=col)
            ax.plot(x, a.n, "_", ms=16, color=INK_MUTED, mew=1.0)
            ax.set_xticks(x, pd.to_datetime(a.image).dt.strftime("%Y-%m"), fontsize=7.5)
            ax.set_xlim(-0.6, len(x) - 0.4)
            ax.set_ylim(0, 95)
            ax.set_ylabel("domains")
            ax.grid(axis="y", color=C["GRID"], lw=0.4)
            open_frame(ax)
            _title(ax, i, f"{WLABEL[w]}, {lab}")
            i += 1
    handles = [Patch(color=OBS_C, label="observed in imagery")]
    handles += [Patch(color=c, label=l) for l, c in SIDES.values()]
    handles += [Line2D([], [], color=INK_MUTED, marker="_", lw=0, ms=12, mew=1.0, label="domains assessed")]
    fig.legend(handles, [h.get_label() for h in handles], loc="outside lower center", ncol=2, frameon=False)
    caption(fig, "When the model overwashes, against when the imagery shows it: for each assessed image, "
                 "the number of GIS domains with overwash, observed (blue) and modelled before (grey) and "
                 "after (purple). Model overwash is any overwash from the storms that ended between the "
                 "previous image and this one (less 7 days); a tick marks how many domains the image "
                 "covers. edgeBE. The 2011-08 image follows Irene, whose overwash came from the sound, "
                 "which Barrier3D does not model.")
    save(fig, FIG_DIR / "adoption_overwash_by_image.png", close=True)


# --- 5. dune crest ------------------------------------------------------------------

BDM_CAP_M = 4.0   # cascade/beach_dune_manager.py:612, m above the berm


def dune_crest():
    d = load("crest")
    lid = pd.read_csv(TABLES / "crest_lidar_2009.csv").set_index("gis").crest_2009_lidar_m_mhw
    berm = 1.34
    fig, axes = plt.subplots(2, 1, figsize=figsize("double", height=5.6), sharex=True, sharey=True,
                             constrained_layout=True)
    fig.get_layout_engine().set(hspace=0.06)
    for i, lab in enumerate(PAIR):
        ax = axes[i]
        ax.fill_between(lid.index, 0, lid.values, step=None, color=C["REF"], alpha=0.12, lw=0, zorder=0)
        ax.plot(lid.index, lid.values, color=C["REF"], lw=1.8, zorder=3)
        sub = run_of(d, "1996_2010", lab)
        for s, (_, col) in SIDES.items():
            m = sub[sub.side == s].set_index("gis").crest_end_m_mhw
            ax.plot(m.index, m.values, color=col, lw=1.3, zorder=4 if s == "after" else 2)
        if lab == "full management":
            ax.axhline(berm + BDM_CAP_M, color=C["ADDED"], lw=0.9, ls="--", zorder=1)
        ax.axhline(berm, color=INK_MUTED, lw=0.6, ls=":")
        ax.text(90.3, berm - 0.1, "berm", color=INK_MUTED, fontsize=7, ha="right", va="top")
        ax.set_ylim(0, 9.5)
        domain_axis(ax, xlabel=(i == 1))
        ax.set_ylabel("dune crest (m MHW)")
        ax.grid(axis="y", color=C["GRID"], lw=0.4)
        _title(ax, i, f"1996–2010 {lab} run: dunes in 2010")
    side_legend(fig, [Line2D([], [], color=C["REF"], lw=1.8, label="2009 lidar"),
                      Line2D([], [], color=C["ADDED"], lw=0.9, ls="--",
                             label="beach/dune manager cap on added sand (4 m above berm)")], ncol=2)
    caption(fig, "The dune crest the 1996-2010 runs end on, against the 2009 lidar that starts the "
                 "2010-2024 window: each GIS domain's median crest over its dune cells, m above MHW. "
                 "edgeBE. Before, dunes grew toward Barrier3D's default ceiling (Dmaxel 3.4 m NAVD88) and "
                 "the roadway manager rebuilt them to 3.0 m MHW; elsewhere they collapsed to the restart "
                 "height. After, each cell grows toward its own starting crest (floor 0.5 m above the "
                 "berm). Dashed orange in (a): the beach/dune manager's cap, 4 m above the berm. Since "
                 "2026-09-28 it limits only the overwash sand the manager puts on the dunes; before that "
                 "fix it clipped whole cells and held Avon and the Tri-Village flat at 5.34 m.")
    save(fig, FIG_DIR / "adoption_dune_crest_2010.png", close=True)


C_1997_FILL = C["LATE_FILL"]


def main():
    apply_style()
    for old in ("adoption_overwash_by_domain",):
        for f in (FIG_DIR / f"{old}.png", FIG_DIR / "supporting" / f"{old}.pdf"):
            f.unlink(missing_ok=True)
    scorecard()
    shoreline_alongshore()
    for lab in PAIR:
        overwash_map(lab)
    overwash_by_image()
    dune_crest()
    print("wrote", *sorted(p.name for p in FIG_DIR.glob("*.png")), sep="\n  ")


if __name__ == "__main__":
    main()
