"""
How well does the model reproduce the overwash record, when, where, and for which storms?

    python scripts/input_prep/8-overwash-analysis/4-vs-model/storms_vs_overwash.py

Skill scores per window; every storm the model runs; hits, misses and false
alarms per image and per domain; the model's dunes alongshore against the storms;
and overwash against storm size. Details: scripts/input_prep/8-overwash-analysis/README.md.

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-10-01
"""
from __future__ import annotations

import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402
from matplotlib.patches import Patch  # noqa: E402
from scipy.stats import spearmanr  # noqa: E402

REPO = next(p for p in Path(__file__).resolve().parents if (p / "pyproject.toml").exists())
sys.path.insert(0, str(REPO / "scripts"))
sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(REPO / "scripts" / "input_prep" / "3-env-forcings" / "3-storms"))
import overwash_vs_model as ovm  # noqa: E402
import storm_figures as sf  # noqa: E402
from site_layer import hat_overwash as ow  # noqa: E402
from site_layer import hat_env_forcings as env  # noqa: E402
from site_layer.hat_figure_style import (C, C_1984, C_1997, DOMAIN_AXIS_LABEL, INK, INK_MUTED,  # noqa: E402
                                         apply_style, caption, figsize, open_frame, save, title,
                                         town_bands)
from site_layer.hatteras_site_config import HATTERAS_DOMAINS as DOM  # noqa: E402

# --- CONFIG ------------------------------------------------------------------
ARM = "managed"                       # the matrix run overwash_vs_model draws
HIT_C, MISS_C, FA_C = ovm.CLASS_COLOURS[:3]      # the agreement classes of overwash_vs_model
OBS_C, MOD_C = "0.15", HIT_C
WINDOW_C = {(1996, 2010): C_1984, (2010, 2024): C_1997}   # earlier red, later blue
FEW_ASSESSED = 30                     # an image assessing fewer domains than this is drawn light / hollow
LABEL_IMAGES = {"2004-05-25": "Isabel", "2011-08-27": "Irene", "2018-08-22": "2018", "2013-04-06": "2012-13"}
OUT = ow.VS_MODEL / "storms_vs_overwash_1996_2024.png"
# -----------------------------------------------------------------------------


# Each model year's lowest dune crest per domain (m MHW), GIS 1-90, from the matrix run
def crest_by_year(window):
    name = ovm.WINDOWS[window][ARM]
    c = np.load(ovm.MATRIX / f"{window[0]}_{window[1]}" / "edgeBE" / name / f"{name}.npz",
                allow_pickle=True)["cascade"][0]
    rows = []
    for gis in range(1, 91):
        b = c.barrier3d[DOM.gis_to_pad(gis)]
        dd = np.asarray(b.DuneDomain, dtype=float)                 # (years, cells, width) dam above berm
        crest = dd.max(axis=2) if dd.ndim == 3 else dd
        low = (crest.min(axis=1) + b._BermEl) * 10                  # the lowest point a storm must clear
        for t in range(1, len(low)):
            rows.append(dict(year=window[0] + t - 1, gis=gis, crest_m=low[t]))
    return pd.DataFrame(rows)


# Observed and modelled per assessed (image, domain), with the outcome class
def cells(window, obs):
    d = ovm.compare(window, ARM, obs, ovm.HEADLINE_THR).dropna(subset=["observed"]).copy()
    o, m = d.observed.astype(int), d.model.astype(int)
    d["outcome"] = np.select([(o == 1) & (m == 1), (o == 1) & (m == 0), (o == 0) & (m == 1)],
                             ["hit", "miss", "false_alarm"], "neither")
    d["window_start"] = window[0]
    return d


# Per image: outcome counts, shares, and the largest storm since the last image
def per_image(window, d, storms, obs):
    imgs = ovm.images_in(window, obs).set_index("Obs_ID")
    dated = storms.EndTime - ovm.GRACE              # overwash_vs_model's dating rule
    rows = []
    for oid, g in d.groupby("obs_id"):
        im = imgs.loc[oid]
        s = storms[(dated > im["from"]) & (dated <= im.date)]
        n = g.outcome.value_counts()
        rows.append(dict(window_start=window[0], image=pd.Timestamp(im.date), assessed=len(g),
                         hits=n.get("hit", 0), misses=n.get("miss", 0), false_alarms=n.get("false_alarm", 0),
                         observed_pct=100 * g.observed.mean(), modelled_pct=100 * g.model.mean(),
                         largest_rhigh_m=s.rhigh_m.max() if len(s) else np.nan, storms=len(s)))
    return pd.DataFrame(rows)


# Hit rate, false-alarm rate and skill (hit - false-alarm) for a set of cells
def skill(d):
    n = d.outcome.value_counts()
    h, m, f, c = (n.get(k, 0) for k in ("hit", "miss", "false_alarm", "neither"))
    pod, pofd = h / max(h + m, 1), f / max(f + c, 1)
    return dict(hit_rate=pod, false_alarm_rate=pofd, skill=pod - pofd, hits=h, misses=m, false_alarms=f, neither=c)


def frac_year(t):
    t = pd.Timestamp(t)
    return t.year + (t.dayofyear - 1) / 365.25


# (a) storms against the model's lowest crests, through time
def panel_storms(ax, storms, crests):
    q = crests.groupby("year").crest_m.quantile([0.1, 0.5, 0.9]).unstack()
    for (a, b) in ovm.WINDOWS:
        qq = q.loc[(q.index >= a) & (q.index <= b)]
        x = np.r_[qq.index.values, qq.index.values[-1] + 1]
        ext = lambda s: np.r_[s.values, s.values[-1]]  # noqa: E731
        ax.fill_between(x, ext(qq[0.1]), ext(qq[0.9]), step="post", color=C["BASE_FILL"], lw=0, zorder=1)
        ax.step(x, ext(qq[0.5]), where="post", color=C["BASE"], lw=1.0, zorder=2)
    for typ in sf.TYPES:
        e = storms[storms.type == typ]
        ax.vlines(e.frac_year, sf.BERM_MHW, e.rhigh_m, color=sf.TYPE_COLOURS[typ], lw=0.6, zorder=3)
        ax.scatter(e.frac_year, e.rhigh_m, s=7, color=sf.TYPE_COLOURS[typ], lw=0, zorder=4)
    ax.axhline(sf.BERM_MHW, color=INK_MUTED, lw=0.6, ls="--", zorder=2)
    ax.set_ylabel("Water level (m MHW)")
    ax.set_ylim(sf.BERM_MHW - 0.2, max(storms.rhigh_m.max(), q[0.9].max()) + 0.3)


# (b) the scorecard: hit rate, false-alarm rate and skill per window
def panel_scores(ax, scores):
    names = [("hit_rate", "Hit\nrate"), ("false_alarm_rate", "False-alarm\nrate"), ("skill", "Skill")]
    w = 0.36
    for j, (win, s) in enumerate(scores.items()):
        xs = np.arange(len(names)) + (j - 0.5) * w
        vals = [s[k] for k, _ in names]
        ax.bar(xs, vals, width=w, color=WINDOW_C[win], lw=0, label=f"{win[0]}–{win[1]}")
        for x, v in zip(xs, vals):
            ax.text(x, max(v, 0) + 0.02, f"{v:.2f}", ha="center", va="bottom", fontsize=7, color=INK)
    ax.set_xticks(np.arange(len(names)), [lab for _, lab in names])
    ax.set_ylim(0, 1.05)
    ax.axhline(0, color=INK, lw=0.6)
    ax.set_ylabel("Score (0–1)")
    ax.legend(loc="upper right", frameon=False, fontsize=7, handlelength=1.0)


# (c) per image, the assessed domains split into hits, misses and false alarms
def panel_images(ax, imgs):
    xs = imgs.image.map(frac_year).values
    few = (imgs.assessed < FEW_ASSESSED).values
    for sel, alpha in ((~few, 1.0), (few, 0.35)):
        x = xs[sel]
        h, m, f = (imgs[k].values[sel] for k in ("hits", "misses", "false_alarms"))
        ax.bar(x, h, width=0.45, color=HIT_C, lw=0, alpha=alpha)
        ax.bar(x, m, width=0.45, bottom=h, color=MISS_C, lw=0, alpha=alpha)
        ax.bar(x, f, width=0.45, bottom=h + m, color=FA_C, lw=0, alpha=alpha)
    ax.set_ylabel("Domains")
    ax.set_ylim(0, 90)
    ax.set_xlim(1996, 2025)


# (d) observed and modelled share against the largest storm between images
def panel_size(ax, imgs):
    k = imgs.dropna(subset=["largest_rhigh_m"])
    ax.vlines(k.largest_rhigh_m, k[["observed_pct", "modelled_pct"]].min(axis=1),
              k[["observed_pct", "modelled_pct"]].max(axis=1), color=C["GRID"], lw=0.8, zorder=1)
    kf = k.assessed < FEW_ASSESSED
    for col, colour in (("observed_pct", OBS_C), ("modelled_pct", MOD_C)):
        ax.scatter(k.largest_rhigh_m[~kf], k[col][~kf], s=16, color=colour, zorder=3)
        ax.scatter(k.largest_rhigh_m[kf], k[col][kf], s=16, facecolor="white", edgecolor=colour, lw=1.0, zorder=3)
    for _, r in k.iterrows():
        lab = LABEL_IMAGES.get(str(r.image.date()))
        if lab:
            ax.annotate(lab, (r.largest_rhigh_m, max(r.observed_pct, r.modelled_pct)), xytext=(0, 4),
                        textcoords="offset points", ha="center", fontsize=7, color=INK)
    ax.set_xlabel("Largest storm between\nimages (m MHW)")
    ax.set_ylabel("Domains overwashed (%)")
    ax.set_ylim(0, 104)


# (e) per domain, the images it was a hit, a miss or a false alarm
def panel_alongshore(ax, d):
    t = d.pivot_table(index="gis", columns="outcome", values="obs_id", aggfunc="count", fill_value=0)
    t = t.reindex(range(1, 91), fill_value=0)
    g = t.index.values
    h, m, f = (t.get(k, pd.Series(0, index=t.index)).values for k in ("hit", "miss", "false_alarm"))
    ax.bar(g, h, width=0.85, color=HIT_C, lw=0)
    ax.bar(g, m, width=0.85, bottom=h, color=MISS_C, lw=0)
    ax.bar(g, f, width=0.85, bottom=h + m, color=FA_C, lw=0)
    ax.set_xlim(0.5, 90.5)
    ax.set_ylabel("Images")
    town_bands(ax)


# (f) per domain, the model's lowest crest against the storms of each window
def panel_crest(ax, crests, storms):
    for win, colour in WINDOW_C.items():
        cr = crests[(crests.year >= win[0]) & (crests.year <= win[1])]
        med = cr.groupby("gis").crest_m.median().reindex(range(1, 91))
        ax.plot(med.index, med.values, color=colour, lw=1.2, zorder=3)
        s = storms[storms.pi == list(ovm.WINDOWS).index(win)]
        annual_max = s.groupby("calendar_year").rhigh_m.max().median()
        ax.axhline(annual_max, color=colour, lw=0.9, ls="--", zorder=2)
    ax.axhline(sf.BERM_MHW, color=INK_MUTED, lw=0.6, ls=":", zorder=2)
    ax.set_xlim(0.5, 90.5)
    ax.set_ylabel("m MHW")
    ax.set_xlabel(DOMAIN_AXIS_LABEL)
    town_bands(ax)


# The figure and its caption
def main():
    apply_style()
    obs = pd.read_csv(ow.OBSERVATIONS / "overwash_observations.csv")
    forcing = sf.load_forcing()
    storms = sf.classify(sf.load_record(env.DEFAULT_STORM_VARIANT, forcing), forcing)
    allcells, imgs, crests, scores = [], [], [], {}
    for window, keep_to in zip(ovm.WINDOWS, (2009, 2024)):
        s = storms[storms.pi == list(ovm.WINDOWS).index(window)]
        d = cells(window, obs)
        allcells.append(d)
        imgs.append(per_image(window, d, s, obs))
        scores[window] = skill(d)
        cr = crest_by_year(window)
        crests.append(cr[cr.year <= keep_to])
    d = pd.concat(allcells, ignore_index=True)
    imgs = pd.concat(imgs, ignore_index=True).sort_values("image").reset_index(drop=True)
    crests = pd.concat(crests, ignore_index=True)

    fig = plt.figure(figsize=figsize("double", height=8.6), constrained_layout=True)
    gs = fig.add_gridspec(4, 2, width_ratios=[2.3, 1], height_ratios=[1.15, 1, 0.9, 0.9])
    ax_a, ax_b = fig.add_subplot(gs[0, 0]), fig.add_subplot(gs[0, 1])
    ax_c, ax_d = fig.add_subplot(gs[1, 0], sharex=ax_a), fig.add_subplot(gs[1, 1])
    ax_e = fig.add_subplot(gs[2, :])
    ax_f = fig.add_subplot(gs[3, :], sharex=ax_e)

    panel_storms(ax_a, storms, crests)
    title(ax_a, 0, "Storms the model runs, against its dunes")
    panel_scores(ax_b, scores)
    title(ax_b, 1, "Skill")
    panel_images(ax_c, imgs)
    title(ax_c, 2, "When: each image")
    panel_size(ax_d, imgs)
    title(ax_d, 3, "By storm size")
    panel_alongshore(ax_e, d)
    title(ax_e, 4, "Where: each domain, all images")
    panel_crest(ax_f, crests, storms)
    title(ax_f, 5, "Why: the model's lowest dune crest against a typical year's largest storm")
    for ax in (ax_a, ax_c):
        for x in imgs.image.map(frac_year):
            ax.axvline(x, color=C["GRID"], lw=0.5, zorder=0)
        ax.axvline(2010, color=INK_MUTED, lw=0.6, ls=":", zorder=0)
    plt.setp(ax_a.get_xticklabels(), visible=False)
    plt.setp(ax_e.get_xticklabels(), visible=False)
    for ax in (ax_a, ax_b, ax_c, ax_d, ax_e, ax_f):
        open_frame(ax)

    handles = [Patch(color=HIT_C, label="Hit (observed + modelled)"),
               Patch(color=MISS_C, label="Miss (observed only)"),
               Patch(color=FA_C, label="False alarm (modelled only)"),
               Line2D([], [], color=sf.TYPE_COLOURS["tropical"], marker="o", ms=3, lw=0.6,
                      label="Tropical cyclone ≤ 500 km"),
               Line2D([], [], color=sf.TYPE_COLOURS["other"], marker="o", ms=3, lw=0.6,
                      label="Other storm"),
               Patch(color=C["BASE_FILL"], label="Lowest crest, 10–90% (line: median)"),
               Line2D([], [], color=OBS_C, marker="o", lw=0, ms=4, label="Observed (d)"),
               Line2D([], [], color=MOD_C, marker="o", lw=0, ms=4, label="Modelled (d)"),
               Line2D([], [], color="0.5", marker="o", mfc="white", lw=0, ms=4,
                      label=f"< {FEW_ASSESSED} domains assessed"),
               Line2D([], [], color=C_1984, lw=1.2, label="Crest, 1996–2010 run (f)"),
               Line2D([], [], color=C_1997, lw=1.2, label="Crest, 2010–2024 run (f)"),
               Line2D([], [], color="0.4", lw=0.9, ls="--", label="Typical yearly max storm (f)")]
    fig.legend(handles=handles, loc="outside lower center", ncol=4, frameon=False)

    full = imgs[imgs.assessed >= FEW_ASSESSED]
    k = full.dropna(subset=["largest_rhigh_m"])
    r_o = spearmanr(k.largest_rhigh_m, k.observed_pct).statistic
    r_m = spearmanr(k.largest_rhigh_m, k.modelled_pct).statistic
    r_om = spearmanr(full.observed_pct, full.modelled_pct).statistic
    sc = "; ".join(f"{w[0]}–{w[1]}: {v['hits']} hits, {v['misses']} misses, {v['false_alarms']} false alarms, "
                   f"{v['neither']} correct negatives" for w, v in scores.items())
    caption(fig, (
        f"How well the hindcast reproduces the overwash record, 1996–2024 (managed matrix runs, storm series "
        f"{env.DEFAULT_STORM_VARIANT}). Each assessed image-and-domain pair is a hit (overwash in the image and "
        f"the model), a miss (image only), a false alarm (model only) or a correct negative; modelled = "
        f"overwash > 0 m³/m in a model year whose largest storm falls between the previous image and this one "
        f"(overwash_vs_model's rule). (a) Every storm the model runs (Duck gauge water level + Stockdon runup, "
        f"slope 0.06), drawn from the berm ({sf.BERM_MHW:.2f} m MHW, dashed) to its peak and coloured by type "
        f"(HURDAT2), over the model's own dunes: each domain's lowest dune crest at the end of each model year, "
        f"10th–90th percentile across GIS 1–90 (band) and median (line). The 1996–2010 run supplies 1996–2009, the "
        f"2010–2024 run 2010–2024 (dotted line). (b) Hit rate = hits / (hits + misses); false-alarm rate = false "
        f"alarms / (false alarms + correct negatives); skill = the difference (Peirce skill score; 0 = no better "
        f"than chance, 1 = perfect). {sc}. (c) Each image's assessed domains by outcome; correct negatives are "
        f"the rest up to the domains assessed. (d) The share of assessed domains overwashed, observed and "
        f"modelled, against the largest storm between images; Spearman r with the storm, images assessing "
        f"≥ {FEW_ASSESSED} domains only: observed {r_o:.2f}, modelled {r_m:.2f}; observed against modelled "
        f"{r_om:.2f}. (e) For each domain, the number of images in which it was a hit, a miss or a false alarm, "
        f"both windows together. (f) Each domain's lowest dune crest in the model, median over the run's years, "
        f"per window, with the median of each year's largest storm (dashed, same colour): where the crest sits "
        f"below the dashed line, a typical year's largest storm overwashes that domain. Berm dotted. Images "
        f"assessing fewer than {FEW_ASSESSED} of the 90 domains (Feb 2022, Apr 2023, Oct 2023) are drawn light in "
        f"(c) and hollow in (d)."))
    out = save(fig, OUT, close=True)
    ow.VS_MODEL_TABLES.mkdir(parents=True, exist_ok=True)
    imgs.to_csv(ow.VS_MODEL_TABLES / "storms_vs_overwash_by_image.csv", index=False)
    pd.DataFrame([dict(window=f"{w[0]}_{w[1]}", **v) for w, v in scores.items()]).to_csv(
        ow.VS_MODEL_TABLES / "storms_vs_overwash_skill.csv", index=False)
    for w, v in scores.items():
        print(w, {k: round(x, 2) if isinstance(x, float) else x for k, x in v.items()})
    print(f"spearman: observed {r_o:.2f}  modelled {r_m:.2f}  obs-vs-mod {r_om:.2f}")
    print(out[0])


if __name__ == "__main__":
    main()
