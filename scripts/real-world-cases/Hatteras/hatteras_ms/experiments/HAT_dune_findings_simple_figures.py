"""
Publication versions of the dune ceiling and growth-rate findings, one point per figure.

    python scripts/hatteras_ms/experiments/HAT_dune_findings_simple_figures.py

Same runs and crest cache as HAT_dune_findings_figures.py (run that first), 1996
start, 1996-2009, managed edgeBE, groin off. Three dune settings recur:
current (5.5 m, r 0.55-0.95), best overall fit from the grid study (5.0 m,
r 0.25-0.45; overwash skill first, 2009 crest as tiebreak) and slow growth
(5.0 m, r 0.03-0.07; best dune trajectory with some growth, low-rates study).
Figures dune_fig01-12, numbered in reading order, written to topic folders
(TOPICS): 1-how-dunes-grow (01 growth toward the ceiling), 2-dune-ceiling
(02 ceiling test, 03 Isabel's reach, 04 measured maximum, 05 P99 ceiling with
slow growth), 3-growth-rate (06 low dunes under the current setting, 07 the
ceiling x growth grid, 08 low dunes under three settings, 09 2009 crests against
the lidar, 10 overwash per image and per rate), 4-test-window (11 low dunes,
12 overwash per image, 2009-2025). Model overwash the yearly imagery does not
show is reported as unconfirmed rather than false. Figures 05 and 09-12 also
show the measured ceilings with slow growth (9.03 and 12.1 m, r 0.03-0.07,
HAT_measured_ceiling_slow_growth.py); those runs have the groin on, so they sit
beside the groin-on runs of the three settings, not the groin-off ladder.
Model overwash is scored with overwash_vs_model.compare at its headline
threshold, as the study tables; the overwash series are cached per run.

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
import pandas as pd  # noqa: E402

_HERE = Path(__file__).resolve()
sys.path.insert(0, str(_HERE.parent))
import HAT_dune_findings_figures as F  # noqa: E402  (puts the overwash scripts on the path)
import HAT_dune_ceiling_uniform as U  # noqa: E402
import overwash_vs_model as ovm  # noqa: E402
from site_layer import hat_figure_style as fs  # noqa: E402
from site_layer import hat_overwash as ow  # noqa: E402

# --- CONFIG ------------------------------------------------------------------
OUT = F.EXP_DIR / "figures"
# Topic subfolder of each figure (figures/README.md explains the order)
TOPICS = {"dune_fig01_how_dunes_grow": "1-how-dunes-grow",
          "dune_fig02_ceiling_test": "2-dune-ceiling", "dune_fig03_isabel_reach": "2-dune-ceiling",
          "dune_fig04_measured_ceiling": "2-dune-ceiling", "dune_fig05_measured_ceiling_slow_growth": "2-dune-ceiling",
          "dune_fig06_low_dunes_current": "3-growth-rate", "dune_fig08_low_dunes_settings": "3-growth-rate",
          "dune_fig09_growth_rate_vs_2009_lidar": "3-growth-rate", "dune_fig10_growth_rate_overwash": "3-growth-rate",
          "dune_fig07_ceiling_growth_grid": "3-growth-rate",
          "dune_fig11_test_window_low_dunes": "4-test-window", "dune_fig12_test_window_overwash": "4-test-window"}
OBS_C = fs.INK
SETTINGS = {   # key: (legend label, run, colour)
    "current": ("Current: 5.5 m, r 0.55–0.95", F.CONTROL, "0.45"),
    "best": ("Best overall fit: 5.0 m, r 0.25–0.45", (F.GRID, "c500_r025_045_s0"), fs.C["ACCENT"]),
    "slow": ("Slow growth: 5.0 m, r 0.03–0.07", (F.LOW, "c500_r003_007_s0"), "#e08214"),
}
MEASURED_MAX = ("12.1 m (1996 DEM maximum), r 0.55–0.95", (F.RANGE, "c1210"), fs.C_1984)
# Growth-rate ladder at the 5.0 m ceiling (the best-fit ceiling), then the current setting
LADDER = [(r, F.LADDER[5.0][r]) for r in
          ["0", "0.01-0.03", "0.03-0.07", "0.05-0.15", "0.10-0.20", "0.15-0.30", "0.25-0.45", "0.35-0.55"]]
LADDER_LABELS = ["None", "0.01–0.03", "0.03–0.07", "0.05–0.15", "0.10–0.20", "0.15–0.30", "0.25–0.45",
                 "0.35–0.55"]
ISABEL_IMAGE = "2004-05-25"
# Test window (2009 start): the same three settings, groin on (the runs the test-period check made)
TEST_RUNS = {"current": (F.GRID, "control_s0_groin"), "best": (F.GRID, "c500_r025_045_s0_groin"),
             "slow": (F.LOW, "c500_r003_007_s0_groin")}
# Measured ceilings with slow growth, groin on, both starts: key -> (legend label, run, colour)
MEASURED_SLOW = {"p99": ("Measured P99: 9.03 m, r 0.03–0.07", (F.LOW, "c903_r003_007_s0_groin"), "#d6604d"),
                 "max": ("Measured maximum: 12.1 m, r 0.03–0.07", (F.LOW, "c1210_r003_007_s0_groin"), "#7f0000")}
TEST_IMAGE_NOTES = {"2011-08-27": "(Irene)", "2013-04-06": "(after Sandy)"}
RANGE_VARIANTS = {4.5: "c450", 5.0: "c500", 5.25: "c525", 5.5: "c550", 5.75: "c575", 6.0: "c600", 8.4: "c840",
                  12.1: "c1210"}
# -----------------------------------------------------------------------------


def colour_for(run):
    named = [(r, c) for _, r, c in SETTINGS.values()] + [(r, c) for _, r, c in MEASURED_SLOW.values()] + \
            [(r, SETTINGS[k][2]) for k, r in TEST_RUNS.items()]
    return next((c for r, c in named if r == run), "0.78")


# Groin-on runs of the 2009-2025 figures: the three settings, then the measured ceilings
def groin_on_runs():
    return [(SETTINGS[k][0], TEST_RUNS[k], SETTINGS[k][2]) for k in TEST_RUNS] + list(MEASURED_SLOW.values())


# The groin-on group appended to the 1996 growth-rate panels: slow growth at 5.0 m and the measured ceilings
GROIN_ON_GROUP = [("0.03–0.07\n(5.0 m)", TEST_RUNS["slow"]), ("0.03–0.07\n(9.03 m)", MEASURED_SLOW["p99"][1]),
                  ("0.03–0.07\n(12.1 m)", MEASURED_SLOW["max"][1])]


def save(fig, name, caption):
    path = OUT / TOPICS[name] / f"{name}.png"
    fs.save(fig, path, close=True)
    fs.record_caption(path, caption + " Producer: HAT_dune_findings_simple_figures.py.")
    print(path, flush=True)


def finish(ax, i=None):
    ax.grid(axis="y", color=fs.GRID_C, lw=0.5)
    fs.open_frame(ax)
    if i is not None:
        fs._letter_inside(ax, i)


def survey_markers(ax, key, label="Observed: lidar surveys"):
    pts = F.survey_points(1996)
    ax.plot(pts.year, pts[key], ls="none", marker="o", ms=5, mfc="white", mec=OBS_C, mew=1.3,
            label=label, zorder=6)


# --- overwash per image ------------------------------------------------------------
def qow(study, variant, start=1996):
    cache = F.EXP_DIR / "cache" / f"{study}__{variant}__{start}__qow.npy"
    if not cache.exists():
        c = U.load(F.run_dir(study, variant, start))
        np.save(cache, np.array([np.asarray(b.QowTS) for b in c.barrier3d]).T)
        print(f"cached {cache.name}", flush=True)
    return np.load(cache)


def per_image(run, obs, start=1996):
    q = qow(*run, start)
    orig = ovm.load_qow
    ovm.load_qow = lambda w, arm: (q, "experiment")
    try:
        d = ovm.compare((start, F.END[start]), "managed", obs, ovm.HEADLINE_THR).dropna(subset=["observed"])
    finally:
        ovm.load_qow = orig
    o, m = d.observed.astype(int), d.model.astype(int)
    d["hit"], d["unconfirmed"] = (o == 1) & (m == 1), (o == 0) & (m == 1)
    return d.groupby("image").agg(observed=("observed", "sum"), hits=("hit", "sum"),
                                  unconfirmed=("unconfirmed", "sum")).astype(int) \
        .rename(index=lambda s: str(s)[:10])


# --- figures -----------------------------------------------------------------------
def fig_low_dunes(keys, name, caption):
    fig, ax = plt.subplots(figsize=fs.figsize("single", height=3.0 + 0.25 * len(keys)), layout="constrained")
    ceil = F.crests(*F.CONTROL, 1996)[1]
    ax.axhline(ceil, color=fs.INK_MUTED, lw=0.8, ls="--", label=f"Current ceiling (Dmaxel {ceil:.1f} m)")
    for k in keys:
        label, run, col = SETTINGS[k]
        cr, _ = F.crests(*run, 1996)
        yrs = F.FIRST[1996] + np.arange(cr.shape[0])
        ax.plot(yrs, np.percentile(cr, 10, axis=1), color=col, lw=1.6, label=f"Model: {label.split(': ')[1]}"
                if len(keys) == 1 else f"Model, {label[0].lower() + label[1:]}")
    survey_markers(ax, "p10")
    ax.set(xlabel="Year (1 January)", ylabel="10th-percentile dune crest (m NAVD88)", ylim=(1.5, 6),
           xlim=(1996.5, 2010.5))
    fig.legend(*ax.get_legend_handles_labels(), frameon=False, fontsize=6.5, loc="outside lower center", ncol=1)
    finish(ax)
    save(fig, name, caption)


def fig_measured():
    rng = pd.read_csv(F.EXPS / F.RANGE / "tables" / "scores.csv")
    rng = rng[rng.period == "1996_2009"].set_index("ceiling_m_navd88")
    cur_label, cur_run, cur_c = SETTINGS["current"]
    runs = [("Model, ceiling 5.5 m (current)", cur_run, cur_c, rng.hits.loc[5.5]),
            ("Model, ceiling 12.1 m (1996 DEM maximum)", MEASURED_MAX[1], MEASURED_MAX[2], rng.hits.loc[12.1])]
    fig, (a, b) = plt.subplots(1, 2, figsize=fs.figsize("double", height=3.0), layout="constrained",
                               gridspec_kw=dict(width_ratios=[1.8, 1]))
    for label, run, col, _ in runs:
        cr, _ = F.crests(*run, 1996)
        a.plot(F.FIRST[1996] + np.arange(cr.shape[0]), np.median(cr, axis=1), color=col, lw=1.6, label=label)
    survey_markers(a, "median")
    a.set(xlabel="Year (1 January)", ylabel="Median dune crest (m NAVD88)", ylim=(4, 13), xlim=(1996.5, 2010.5))
    a.legend(frameon=False, fontsize=6.5, loc="center right")
    finish(a, 0)

    labels = ["Observed\n(imagery)", "Model\n5.5 m", "Model\n12.1 m"]
    vals = [int(rng.observed_share.size and 68)] + [int(r[3]) for r in runs]
    cols = [OBS_C] + [r[2] for r in runs]
    bars = b.bar(labels, vals, color=cols, width=0.62, edgecolor=cols)
    bars[0].set_facecolor("white")
    bars[0].set_hatch("////")
    for r, v in zip(bars, vals):
        b.text(r.get_x() + r.get_width() / 2, v + 1.5, str(v), ha="center", fontsize=7.5)
    b.set(ylabel="Overwashed domains, 1996–2009\n(observed; reproduced by the model)", ylim=(0, 78))
    finish(b, 1)
    save(fig, "dune_fig04_measured_ceiling",
         "Dune ceiling set to the measured maximum. (a) Median dune crest of the 4,500 10 m cells each 1 January "
         "from the 1996-start runs with the current ceiling (Dmaxel 5.5 m NAVD88, grey) and with the tallest dune "
         "crest in the 1996 DEM (12.1 m, red), both with growth 0.55-0.95, against the lidar surveys (circles). "
         "(b) Overwash in 1996-2009 as image x domain cases: observed in the imagery (68, hatched) and the number "
         "the model reproduces.")


def box_stats(v, label):
    return dict(label=label, med=np.median(v), q1=np.percentile(v, 25), q3=np.percentile(v, 75),
                whislo=np.percentile(v, 10), whishi=np.percentile(v, 90), fliers=[])


def fig_growth_dunes():
    lid = F.lidar_2009()
    k = lid.year - F.FIRST[1996]
    ok = np.isfinite(lid.obs)
    o = lid.obs[ok]
    runs = LADDER + [("current", SETTINGS["current"][1])] + GROIN_ON_GROUP
    labels = ["2009\nlidar"] + LADDER_LABELS + ["0.55–0.95\n(5.5 m,\ncurrent)"] + [g[0] for g in GROIN_ON_GROUP]
    vals, cols = [o], [OBS_C]
    for _, run in runs:
        cr, _ = F.crests(*run, 1996)
        vals.append(cr[k][ok])
        cols.append(colour_for(run))
    n = len(LADDER)
    pos = np.r_[0, np.arange(1, n + 1) + 0.4, n + 1.8, n + 3.2 + np.arange(len(GROIN_ON_GROUP))]
    fig, (a, b) = plt.subplots(2, 1, figsize=fs.figsize("double", height=4.8), layout="constrained", sharex=True,
                               gridspec_kw=dict(height_ratios=[1.6, 1]))
    a.axhspan(np.percentile(o, 10), np.percentile(o, 90), color="0.93", zorder=0)
    a.axhline(np.median(o), color=OBS_C, lw=0.8, ls="--", zorder=1)
    bp = a.bxp([box_stats(v, l) for v, l in zip(vals, labels)], positions=pos, widths=0.55, patch_artist=True,
               showfliers=False, medianprops=dict(color=fs.INK, lw=1.2))
    for patch, col in zip(bp["boxes"], cols):
        patch.set(facecolor="white" if col == OBS_C else col, edgecolor=fs.INK, lw=0.8)
    bp["boxes"][0].set_hatch("////")
    a.set(ylabel="Dune crest, 1 Jan 2010 (m NAVD88)", ylim=(1.5, 9.5))
    a.text(pos[-1] + 0.45, np.median(o), "lidar median", fontsize=6, va="center", color=OBS_C)
    g0 = pos[-len(GROIN_ON_GROUP)]
    for ax_ in (a, b):
        ax_.axvline(g0 - 0.7, color=fs.INK_MUTED, lw=0.6, ls=":")
    a.text((g0 + pos[-1]) / 2, 1.75, "groin on", ha="center", fontsize=6.5, color=fs.INK_MUTED)
    finish(a, 0)

    design = F.DESIGN_NAVD88
    share = [(v < design).mean() * 100 for v in vals]
    bars = b.bar(pos, share, width=0.55, color=["white"] + cols[1:], edgecolor=fs.INK, lw=0.8)
    bars[0].set_hatch("////")
    b.axhline(share[0], color=OBS_C, lw=0.8, ls="--")
    for x, s in zip(pos, share):
        b.text(x, s + 0.8, f"{s:.0f}", ha="center", fontsize=6.5)
    b.set(ylabel=f"Cells below the NC-12\nrebuild height, {design} m (%)", ylim=(0, max(share) * 1.25))
    b.set_xticks(pos, labels, fontsize=6.5)
    b.set_xlabel("Dune growth rate, rmin–rmax (ceiling 5.0 m unless noted)")
    finish(b, 1)
    handles = [plt.Rectangle((0, 0), 1, 1, facecolor="white", edgecolor=fs.INK, hatch="////")] + \
              [plt.Rectangle((0, 0), 1, 1, facecolor=c, edgecolor=fs.INK) for c in
               (SETTINGS["best"][2], SETTINGS["slow"][2], SETTINGS["current"][2], "0.78",
                MEASURED_SLOW["p99"][2], MEASURED_SLOW["max"][2])]
    fig.legend(handles, ["Observed (2009 lidar)", "Model, best overall fit", "Model, slow growth",
                       "Model, current setting", "Model, other rates", "Model, measured P99 ceiling",
                       "Model, measured maximum ceiling"], frameon=False, fontsize=6.5, loc="outside lower center", ncol=4)
    save(fig, "dune_fig09_growth_rate_vs_2009_lidar",
         "Modelled dune crests on 1 January 2010 from the 1996 start against the 2009 lidar on the same 10 m "
         "cells, for each growth rate (Dmaxel 5.0 m) and the current setting (5.5 m, 0.55-0.95), groin off. Right "
         "of the dotted line, groin on: slow growth (0.03-0.07) with the 5.0 m ceiling and with the ceiling at the "
         "1996 DEM's 99th-percentile (9.03 m) and tallest (12.1 m) crest; compare these three with each other, "
         "as the groin alone moves the slow-growth 10th percentile by about 0.5 m. (a) Crest distribution: boxes "
         "span the 25th-75th percentiles, whiskers the 10th-90th, line the median; grey band and dashed line are "
         "the lidar 10th-90th percentile range and median. (b) Share of cells below the NC-12 rebuild height "
         "(4.3 m NAVD88). Lidar cells where the dune had become beach count as 1.7 m (the berm).")
    return pd.DataFrame(dict(run=labels, median=[np.median(v) for v in vals],
                             p10=[np.percentile(v, 10) for v in vals], share_below_4p3_pct=share))


def fig_growth_overwash():
    obs = pd.read_csv(ow.OBSERVATIONS / "overwash_observations.csv")
    tabs = {k: per_image(run, obs) for k, (_, run, _) in SETTINGS.items()}
    ref = tabs["current"]
    seen = ref[ref.observed > 0].index.tolist()

    fig, (a, b) = plt.subplots(1, 2, figsize=fs.figsize("double", height=3.4), layout="constrained",
                               gridspec_kw=dict(width_ratios=[1.15, 1]))
    x = np.arange(len(seen))
    w = 0.2
    bars = a.bar(x - 1.5 * w, ref.loc[seen].observed, w, color="white", edgecolor=OBS_C, hatch="////",
                 label="Observed (imagery)")
    for j, (k, (label, _, col)) in enumerate(SETTINGS.items()):
        a.bar(x + (j - 0.5) * w, tabs[k].loc[seen].hits, w, color=col, label=f"Model, {label.split(':')[0].lower()}")
    for r in bars:
        a.text(r.get_x() + r.get_width() / 2, r.get_height() + 1, f"{r.get_height():.0f}", ha="center", fontsize=6)
    a.set_xticks(x, [pd.Timestamp(s).strftime("%b %Y") + ("\n(incl. Isabel)" if s == ISABEL_IMAGE else "")
                     for s in seen], fontsize=6.5)
    a.set(xlabel="Image date", ylabel="Overwashed domains (of 90)", ylim=(0, 56))
    a.legend(frameon=False, fontsize=6.5, loc="upper right")
    finish(a, 0)

    rows = []
    groin_on = [((None, run), lab) for lab, run in GROIN_ON_GROUP]
    for (rate, run), lab in list(zip(LADDER, LADDER_LABELS)) + [(("current", SETTINGS["current"][1]), "0.55–0.95\n(5.5 m)")] \
            + groin_on:
        t = per_image(run, obs)
        rows.append(dict(label=lab, run="/".join(run), isabel=t.loc[ISABEL_IMAGE].hits,
                         other=t.drop(ISABEL_IMAGE).hits.sum(), unconfirmed=t.unconfirmed.sum(),
                         colour=colour_for(run)))
    s = pd.DataFrame(rows)
    xs = np.arange(len(s)) + np.r_[np.zeros(len(s) - len(groin_on)), np.full(len(groin_on), 0.8)]
    b.axvline(xs[-len(groin_on)] - 0.9, color=fs.INK_MUTED, lw=0.6, ls=":")
    b.text(xs[-2], 132, "groin on", ha="center", fontsize=6.5, color=fs.INK_MUTED)
    n_isabel, n_other = ref.loc[ISABEL_IMAGE].observed, ref.drop(ISABEL_IMAGE).observed.sum()
    b.bar(xs, s.isabel, 0.6, color=s.colour, edgecolor=fs.INK, lw=0.6, label=f"Reproduced, Isabel image (of {n_isabel})")
    b.bar(xs, s.other, 0.6, bottom=s.isabel, color=s.colour, edgecolor=fs.INK, lw=0.6, hatch="....",
          label=f"Reproduced, other images (of {n_other})")
    b.axhline(n_isabel + n_other, color=OBS_C, lw=0.8, ls="--")
    b.text(xs[-1] + 0.4, n_isabel + n_other + 2, f"observed: {n_isabel + n_other}", ha="right", fontsize=6.5)
    b2 = b.twinx()
    b2.plot(xs, s.unconfirmed, color=fs.INK, lw=1.0, marker="o", ms=3.5, mfc="white",
            label="Model-only overwash, unconfirmed (right axis)")
    b2.set_ylabel("Model-only overwashed domains\n(unconfirmed, ~1 image per year)", fontsize=7)
    b2.set_ylim(0, 140)
    b2.spines["top"].set_visible(False)
    b.set(ylim=(0, 140), ylabel="Observed overwash reproduced (domains)")
    b.set_xticks(xs, s.label.str.replace(chr(10), " "), rotation=55, ha="right", rotation_mode="anchor", fontsize=6)
    b.set_xlabel("Dune growth rate, rmin–rmax (ceiling 5.0 m unless noted)")
    h1, l1 = b.get_legend_handles_labels()
    h2, l2 = b2.get_legend_handles_labels()
    fig.legend(h1 + h2, l1 + l2, frameon=False, fontsize=6.5, loc="outside lower right", ncol=3)
    finish(b, 1)
    save(fig, "dune_fig10_growth_rate_overwash",
         "Overwash in 1996-2009, 1996-start runs. (a) Domains overwashed in each image that shows overwash "
         "(hatched) and the number each dune setting reproduces; the May 2004 image spans 1998-2004 and includes "
         "Hurricane Isabel (2003), 48 of the 68 observed cases. (b) Across growth rates (colours as in a for the "
         "three named settings): observed overwash reproduced, Isabel image (solid) and the other images "
         "(dotted), and model overwash in domains the imagery shows as dry (line, right axis). Right of the dotted "
         "line, groin on: slow growth (0.03-0.07) with the 5.0 m ceiling and with the measured 99th-percentile "
         "(9.03 m, light red) and maximum (12.1 m, dark red) ceilings. With about one "
         "image a year, model-only overwash may include minor events the imagery missed, so it is reported as "
         "unconfirmed rather than as false alarms.")
    return s.drop(columns="colour")


# --- the tests, one figure each ----------------------------------------------------
def fig_growth_explainer():
    berm, ceil, h0 = F.BERM_NAVD88, 5.5, 3.9          # h0: the 1996 10th-percentile crest
    dmax = ceil - berm
    years = np.arange(0, 16)
    fig, ax = plt.subplots(figsize=fs.figsize("single", height=3.0), layout="constrained")
    ax.axhline(ceil, color=fs.INK_MUTED, lw=0.8, ls="--")
    ax.text(15, ceil + 0.03, "Dune ceiling (Dmaxel 5.5 m)", ha="right", va="bottom", fontsize=7, color=fs.INK_MUTED)
    # label position per curve (text x, y), chosen to sit clear of the other curves
    for r, col, lab, xy in ((0.75, SETTINGS["current"][2], "r = 0.75 (current)", (0.4, 5.76)),
                            (0.35, SETTINGS["best"][2], "r = 0.35", (6.3, 5.0)),
                            (0.05, SETTINGS["slow"][2], "r = 0.05", (9.5, 4.62))):
        h = [h0 - berm]
        for _ in years[1:]:
            h.append(h[-1] + r * h[-1] * (1 - h[-1] / dmax))
        z = np.array(h) + berm
        ax.plot(years, z, color=col, lw=2)
        near = np.flatnonzero(z >= ceil - 0.1)
        if near.size:
            k = near[0]
            ax.plot(years[k], z[k], marker="o", ms=5, color=col, zorder=5)
            ax.annotate(f"{lab}: at the ceiling after {years[k]} years", (years[k], z[k]), xytext=xy, fontsize=7,
                        color=col, fontweight="bold", arrowprops=dict(arrowstyle="-", color=col, lw=0.6))
        else:
            ax.text(*xy, f"{lab}: still {ceil - z[-1]:.1f} m below\nthe ceiling after {years[-1]} years",
                    fontsize=7, color=col, fontweight="bold")
    ax.set_xticks(range(0, 16, 3))
    ax.set(xlabel="Years without a storm", ylabel="Height of a low dune (m NAVD88)", xlim=(0, 15), ylim=(3.6, 5.95))
    finish(ax)
    save(fig, "dune_fig01_how_dunes_grow",
         "How fast a low dune grows in Barrier3D. A dune starting at 3.9 m NAVD88 (the lowest 10% of dunes in "
         "1996) grows each year without a storm by r H (1 - H / Dmax), where H is its height above the berm "
         "(1.7 m) and Dmax = 5.5 - 1.7 = 3.8 m, until it reaches the ceiling. Curves for three growth rates r, the "
         "middles of the current range (0.55-0.95), the best-fit range (0.25-0.45) and the slow range "
         "(0.03-0.07), all with the 5.5 m ceiling for comparison. Dots: within 0.1 m of the ceiling. Dunes taller "
         "than the ceiling are cut down to it.")


def fig_ceiling_test():
    rng = pd.read_csv(F.EXPS / F.RANGE / "tables" / "scores.csv")
    rng = rng[(rng.period == "1996_2009") & rng.ceiling_m_navd88.isin(list(RANGE_VARIANTS))] \
        .sort_values("ceiling_m_navd88")
    isabel = float(F.yearly_runup(1996).loc[2003])
    fig, ax = plt.subplots(figsize=fs.figsize("double", height=2.8), layout="constrained")
    x = np.arange(len(rng))
    cols = [SETTINGS["current"][2] if c == 5.5 else ("0.78" if c < isabel else MEASURED_MAX[2])
            for c in rng.ceiling_m_navd88]
    bars = ax.bar(x, rng.hits, 0.62, color=cols, edgecolor=fs.INK, lw=0.6)
    for r, h in zip(bars, rng.hits):
        ax.text(r.get_x() + r.get_width() / 2, h + 1.2, f"{h:.0f}", ha="center", fontsize=7)
    ax.axhline(68, color=OBS_C, lw=0.8, ls="--")
    ax.text(len(rng) - 0.5, 69.5, "Observed: 68", ha="right", fontsize=7)
    split = np.searchsorted(rng.ceiling_m_navd88.to_numpy(), isabel) - 0.5
    ax.axvline(split, color=fs.C_1997, lw=0.9, ls="-.")
    ax.text(split + 0.08, 58, f"Isabel run-up\n{isabel:.2f} m", ha="left", fontsize=6.5, color=fs.C_1997)
    names = {5.5: "\n(current)", 8.4: "\n(DEM P95)", 12.1: "\n(DEM max)"}
    ax.set_xticks(x, [f"{c:g}{names.get(c, '')}" for c in rng.ceiling_m_navd88], fontsize=7)
    ax.set(xlabel="Dune ceiling, Dmaxel (m NAVD88); growth 0.55–0.95",
           ylabel="Observed overwash reproduced\n(domains, 1996–2009)", ylim=(0, 78))
    finish(ax)
    save(fig, "dune_fig02_ceiling_test",
         "Ceiling test (2026-10-13): observed overwash cases in 1996-2009 (image x domain, 68 in all) that the "
         "1996-start run reproduces, for uniform dune ceilings from 4.5 m to the tallest crest in the 1996 DEM "
         "(12.1 m), all with the current growth rate. Any ceiling up to 5.5 m reproduces most of them; above "
         "Isabel's run-up (5.58 m, dash-dot) the dunes grow taller than the only storm able to overtop them and "
         "overwash almost stops. Red: ceilings above Isabel's run-up; dark grey: the current setting.")


def fig_isabel_reach():
    isabel = float(F.yearly_runup(1996).loc[2003])
    k = 2003 - F.FIRST[1996]
    obs = pd.read_csv(ow.OBSERVATIONS / "overwash_observations.csv")
    n_obs = int(per_image(F.CONTROL, obs).loc[ISABEL_IMAGE].observed)
    rows = []
    for c in (4.5, 5.0, 5.5, 5.75, 6.0, 8.4, 12.1):
        cr, _ = F.crests(F.RANGE, RANGE_VARIANTS[c], 1996)
        lowest = cr[k].reshape(90, 50).min(axis=1)
        rows.append(dict(ceiling=c, reachable=int((lowest < isabel).sum())))
    t = pd.DataFrame(rows)
    fig, ax = plt.subplots(figsize=fs.figsize("single", height=3.0), layout="constrained")
    x = np.arange(len(t) + 1)
    vals = [n_obs] + t.reachable.tolist()
    cols = ["white"] + [SETTINGS["current"][2] if c == 5.5 else ("0.78" if c < isabel else MEASURED_MAX[2])
                        for c in t.ceiling]
    bars = ax.bar(x, vals, 0.62, color=cols, edgecolor=fs.INK, lw=0.6)
    bars[0].set_hatch("////")
    for r, v in zip(bars, vals):
        ax.text(r.get_x() + r.get_width() / 2, v + 1.2, str(v), ha="center", fontsize=7)
    ax.set_xticks(x, ["Observed\n(May 2004\nimage)"] + [f"{c:g} m" for c in t.ceiling], fontsize=6.5)
    ax.set(ylabel="Domains (of 90)", ylim=(0, 100),
           xlabel="Model: domains whose lowest dune is below\nIsabel's run-up, by ceiling (growth 0.55–0.95)")
    finish(ax)
    save(fig, "dune_fig03_isabel_reach",
         f"Why the ceiling test has a cliff. Hatched: domains overwashed in the May 2004 image, which spans "
         f"Hurricane Isabel ({n_obs} of 90). Bars: domains whose lowest 10 m dune cell on 1 January 2003, just "
         f"before Isabel, sits below Isabel's run-up ({isabel:.2f} m NAVD88), for each ceiling. With ceilings at "
         "or below 5.5 m almost every domain has a dune Isabel can overtop; once the ceiling passes the run-up, "
         "the dunes have grown above it and only a few domains remain open.")
    return t


def fig_p99_slow():
    obs = pd.read_csv(ow.OBSERVATIONS / "overwash_observations.csv")
    # (bar label, legend label or None for bars only, run, colour, line width)
    runs = [("5.5 m\n0.55–0.95\n(current)", "Model, 5.5 m, r 0.55–0.95 (current)", TEST_RUNS["current"],
             SETTINGS["current"][2], 1.6),
            ("P99\n0.15–0.30", None, (F.DEMG, "p99_r015_030"), "#f4c7a8", 0),
            ("P99\n0.10–0.20", None, (F.DEMG, "p99_r010_020"), "#f4c7a8", 0),
            ("P99\n0.05–0.15", "Model, 9.03 m (P99), r 0.05–0.15 (earlier test)", (F.DEMG, "p99_r005_015"),
             "#f4c7a8", 1.2),
            ("P99\n0.03–0.07", "Model, 9.03 m (P99), r 0.03–0.07", MEASURED_SLOW["p99"][1],
             MEASURED_SLOW["p99"][2], 1.6),
            ("Max\n0.03–0.07", "Model, 12.1 m (maximum), r 0.03–0.07", MEASURED_SLOW["max"][1],
             MEASURED_SLOW["max"][2], 1.6),
            ("5.0 m\n0.03–0.07\n(slow)", "Model, 5.0 m, r 0.03–0.07 (slow growth)", TEST_RUNS["slow"],
             SETTINGS["slow"][2], 1.6)]
    fig, (a, b) = plt.subplots(1, 2, figsize=fs.figsize("double", height=3.3), layout="constrained",
                               gridspec_kw=dict(width_ratios=[1.35, 1.35]))
    for _, label, run, col, lw in runs:
        if label:
            cr, _ = F.crests(*run, 1996)
            a.plot(F.FIRST[1996] + np.arange(cr.shape[0]), np.median(cr, axis=1), color=col, lw=lw, label=label)
    survey_markers(a, "median")
    a.set(xlabel="Year (1 January)", ylabel="Median dune crest (m NAVD88)", ylim=(4, 9.5), xlim=(1996.5, 2010.5))
    a.legend(frameon=False, fontsize=6, loc="upper left", bbox_to_anchor=(0.07, 1.0))
    finish(a, 0)
    vals = [68] + [int(per_image(r[2], obs).hits.sum()) for r in runs]
    cols = ["white"] + [r[3] for r in runs]
    bars = b.bar(np.arange(len(vals)), vals, 0.62, color=cols, edgecolor=fs.INK, lw=0.6)
    bars[0].set_hatch("////")
    for r, v in zip(bars, vals):
        b.text(r.get_x() + r.get_width() / 2, v + 1.2, str(v), ha="center", fontsize=6.5)
    b.set_xticks(np.arange(len(vals)), ["Observed"] + [r[0] for r in runs], fontsize=5.5, rotation=45, ha="right",
                 rotation_mode="anchor")
    b.set(ylabel="Observed overwash reproduced\n(domains, 1996–2009)", ylim=(0, 78))
    finish(b, 1)
    save(fig, "dune_fig05_measured_ceiling_slow_growth",
         "Ceiling at the measured dune height with slower growth, 1996 start. The ceiling set to the 99th "
         "percentile (P99, 9.03 m NAVD88) or the tallest (12.1 m) crest of the 1996 DEM, with growth slowed as far "
         "as the slow-growth setting (0.03-0.07; runs of 2026-10-20), against the current setting and slow growth "
         "with the 5.0 m ceiling. (a) Median crest each 1 January against the lidar surveys (circles). (b) Observed "
         "overwash cases reproduced (of 68, hatched). The 0.03-0.07 and current runs have the groin on; the "
         "earlier P99 runs (0.05-0.30, 2026-10-19) have it off.")


def fig_grid():
    g = pd.read_csv(F.EXPS / F.GRID / "tables" / "stage1_grid_1996.csv")
    g = g[~g.groin]
    ctrl = float(g[g.variant == "control_s0"].skill.iloc[0])
    g = g[g.variant != "control_s0"]
    ceilings = sorted(g.ceiling_m_navd88.unique())
    rates = sorted({(a, b) for a, b in zip(g.rmin, g.rmax)})
    m = np.array([[g[(g.ceiling_m_navd88 == c) & (g.rmin == r[0])].skill.iloc[0] for r in rates]
                  for c in ceilings])
    fig, ax = plt.subplots(figsize=fs.figsize("single", height=2.8), layout="constrained")
    im = ax.imshow(m, cmap="Greys", vmin=0, vmax=0.8, aspect="auto", origin="lower")
    for i, c in enumerate(ceilings):
        for j, r in enumerate(rates):
            v = m[i, j]
            best = c == 5.0 and r == (0.25, 0.45)
            ax.text(j, i, f"{v:.2f}", ha="center", va="center", fontsize=8,
                    color="white" if v > 0.45 else fs.INK, fontweight="bold" if best else "normal")
            if best:
                ax.add_patch(plt.Rectangle((j - 0.48, i - 0.48), 0.96, 0.96, fill=False, ec=SETTINGS["best"][2],
                                           lw=2))
    ax.set_xticks(range(len(rates)), [f"{a:.2f}–{b:.2f}" for a, b in rates])
    ax.set_yticks(range(len(ceilings)), [f"{c:g} m" for c in ceilings])
    ax.set(xlabel="Dune growth rate, rmin–rmax", ylabel="Dune ceiling, Dmaxel")
    cb = fig.colorbar(im, ax=ax, shrink=0.85)
    cb.set_label("Overwash skill, 1996–2009", fontsize=7)
    cb.ax.axhline(ctrl, color=SETTINGS["current"][2], lw=2)
    ax.grid(False)
    save(fig, "dune_fig07_ceiling_growth_grid",
         "Ceiling x growth-rate grid (2026-10-19), 1996-2009 overwash skill (share of observed overwash reproduced "
         "minus share of dry image-domains the model overwashes) for each combination, 1996 start. Ceilings of "
         f"6.0 m and above fail at every growth rate (Isabel). The current setting scores {ctrl:.2f} (line on the "
         "colour bar). The four 5.0-5.5 m cells with growth 0.25-0.55 tie within seed noise (0.008); 5.0 m with "
         "0.25-0.45 (outlined) was chosen because its 2009 dune crests match the lidar best.")


def fig_test_dunes():
    pts = F.survey_points(2009)
    fig, ax = plt.subplots(figsize=fs.figsize("single", height=4.0), layout="constrained")
    ax.axhline(5.5, color=fs.INK_MUTED, lw=0.8, ls="--", label="Current ceiling (Dmaxel 5.5 m)")
    for label, run, col in groin_on_runs():
        cr, _ = F.crests(*run, 2009)
        ax.plot(F.FIRST[2009] + np.arange(cr.shape[0]), np.percentile(cr, 10, axis=1), color=col, lw=1.6,
                label=f"Model, {label[0].lower() + label[1:]}")
    ax.plot(pts.year, pts.p10, ls="none", marker="o", ms=5, mfc="white", mec=OBS_C, mew=1.3,
            label="Observed: lidar surveys", zorder=6)
    ax.set(xlabel="Year (1 January)", ylabel="10th-percentile dune crest (m NAVD88)", ylim=(1.5, 7),
           xlim=(2009.5, 2026.5))
    fig.legend(*ax.get_legend_handles_labels(), frameon=False, fontsize=6.5, loc="outside lower center")
    finish(ax)
    save(fig, "dune_fig11_test_window_low_dunes",
         "Test window, 2009-2025 (not used to choose the settings): low dunes (10th percentile crest) of the "
         "2009-start runs, groin on, against the 2014-2019 lidar surveys on the 2009 dune windows (circles). As "
         "in 1996-2009 the current setting builds the low dunes to the ceiling; slow growth stays nearest the "
         "surveys but ends somewhat low, as the real low dunes partly recovered after 2009. Red: slow "
         "growth with the ceiling at the 1996 DEM's 99th-percentile (9.03 m) and tallest (12.1 m) crest.")


def fig_test_overwash():
    obs = pd.read_csv(ow.OBSERVATIONS / "overwash_observations.csv")
    tabs = {k: per_image(run, obs, 2009) for k, run in TEST_RUNS.items()}
    tabs.update({k: per_image(run, obs, 2009) for k, (_, run, _) in MEASURED_SLOW.items()})
    ref = tabs["current"]
    seen = ref[ref.observed >= 5].index.tolist()
    fig, ax = plt.subplots(figsize=fs.figsize("double", height=3.0), layout="constrained")
    x = np.arange(len(seen))
    w = 0.14
    bars = ax.bar(x - 2.5 * w, ref.loc[seen].observed, w, color="white", edgecolor=OBS_C, hatch="////",
                  label="Observed (imagery)")
    for j, (k, (label, _, col)) in enumerate(list(SETTINGS.items()) + list(MEASURED_SLOW.items())):
        ax.bar(x + (j - 1.5) * w, tabs[k].loc[seen].hits, w, color=col,
               label=f"Model, {label[0].lower() + label[1:]}")
    for r in bars:
        ax.text(r.get_x() + r.get_width() / 2, r.get_height() + 1, f"{r.get_height():.0f}", ha="center", fontsize=6)
    ax.set_xticks(x, [pd.Timestamp(s).strftime("%b %Y") + ("\n" + TEST_IMAGE_NOTES[s] if s in TEST_IMAGE_NOTES
                                                           else "") for s in seen], fontsize=6.5)
    ax.set(xlabel="Image date", ylabel="Overwashed domains (of 90)", ylim=(0, max(ref.observed) * 1.2))
    ax.legend(frameon=False, fontsize=6.5, loc="upper right")
    finish(ax)
    tot = {k: (int(t.hits.sum()), int(t.unconfirmed.sum())) for k, t in tabs.items()}
    save(fig, "dune_fig12_test_window_overwash",
         "Test window, 2009-2025: domains overwashed in each image with at least five overwashed domains "
         "(hatched) and the number each setting reproduces (2009-start runs, groin on). Over all images the "
         f"model reproduces {tot['current'][0]} (current), {tot['best'][0]} (best overall fit) and "
         f"{tot['slow'][0]} (slow growth), {tot['p99'][0]} (9.03 m P99 ceiling, slow growth) and {tot['max'][0]} "
         f"(12.1 m maximum ceiling, slow growth) of {int(ref.observed.sum())} observed cases, with "
         f"{tot['current'][1]}, {tot['best'][1]}, {tot['slow'][1]}, {tot['p99'][1]} and {tot['max'][1]} model-only "
         "(unconfirmed) cases. Most "
         "test-window misses are sound-side overwash, which no dune setting can reach.")
    return pd.concat({k: t for k, t in tabs.items()}, axis=1)


def main():
    fs.apply_style()
    OUT.mkdir(parents=True, exist_ok=True)
    fig_low_dunes(["current"], "dune_fig06_low_dunes_current",
                  "Low dunes through 1996-2009: 10th percentile of the 4,500 10 m dune-crest cells each 1 January, "
                  "1996-start run with the current settings (Dmaxel 5.5 m, growth 0.55-0.95), against the lidar "
                  "surveys on the same cells (circles). The modelled low dunes reach the ceiling within about "
                  "three years; the surveyed low dunes fell from about 4.0 m to 2.4 m. The 2004 dip is Isabel.")
    fig_low_dunes(["current", "best", "slow"], "dune_fig08_low_dunes_settings",
                  "As dune_fig06 with three dune settings: current (5.5 m, 0.55-0.95, grey), best overall fit from "
                  "the ceiling x growth grid (5.0 m, 0.25-0.45, purple; ranked on 1996-2009 overwash skill, 2009 "
                  "crest as tiebreak) and slow growth (5.0 m, 0.03-0.07, orange; best low-dune trajectory among the "
                  "settings with growth).")
    fig_measured()
    t3 = fig_growth_dunes()
    t4 = fig_growth_overwash()
    fig_growth_explainer()
    fig_ceiling_test()
    t7 = fig_isabel_reach()
    fig_p99_slow()
    fig_grid()
    fig_test_dunes()
    t11 = fig_test_overwash()
    t7.to_csv(F.EXP_DIR / "tables" / "dune_fig03_isabel_reach.csv", index=False)
    t11.to_csv(F.EXP_DIR / "tables" / "dune_fig12_test_window_overwash.csv")
    print(t7.to_string(index=False))
    print(t11.to_string())
    (F.EXP_DIR / "tables").mkdir(exist_ok=True)
    t3.round(3).to_csv(F.EXP_DIR / "tables" / "dune_fig09_growth_rate_vs_2009_lidar.csv", index=False)
    t4.to_csv(F.EXP_DIR / "tables" / "dune_fig10_growth_rate_overwash.csv", index=False)
    with pd.option_context("display.width", 200):
        print(t3.round(2).to_string(index=False))
        print(t4.to_string(index=False))


if __name__ == "__main__":
    main()
