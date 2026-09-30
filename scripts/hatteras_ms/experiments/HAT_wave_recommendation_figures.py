"""Figures behind the recommended hindcast wave climate (2026-09-27).

Hannah, 2026-09-27: "given all of the different tests ... what do you suggest
as the best wave parameters to use for the model hindcast. Provide figures and
a clear and detailed explanation". Draws, from the studies already run (no new
runs), into output/raw_runs/experiments/wave-climate/2026-09-27-wave-recommendation/figures/:

  1_score_by_parameter.png      best raw score reachable at each value of each
                                parameter (the complete zeroBE grid, coarse +
                                refine), both windows, both scenarios
  2_agreement_across_tests.png  the best 1996-2010 setting of each search, by
                                parameter
  3_recommended_vs_coastsat.png the recommended setting against CoastSat,
                                rate and position change, both windows
  4_high_angle_roughness.png    domain-to-domain roughness and score against
                                the high-angle fraction
  5_window_2010.png             the recommended setting scored on 2010-2024
                                and on 2010-2020 (before the 2021 CoastSat step)

Scores are the RAW share of the alongshore variation explained (the per-domain
model against the CoastSat LOWESS-10 target, interior GIS 2-89), as the runner
and the matrix runs are scored.

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-09-30
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402

_HERE = Path(__file__).resolve()
sys.path.insert(0, str(_HERE.parent))
import HAT_wave_grid_smoothed_score as G  # noqa: E402
import HAT_metres_2_wave_sensitivity_plot as p2  # noqa: E402
from site_layer.hat_figure_style import (  # noqa: E402
    INK, INK_MUTED, DOMAIN_AXIS_LABEL, _title, apply_style, open_frame, record_caption, save,
    C, STRUCTURE_LABEL_PT, structures, support_dir, town_bands)
from site_layer.hatteras_site_config import (  # noqa: E402
    HATTERAS_ANNOTATIONS, HATTERAS_NOURISHMENT_PROJECTS)

common = G.common
EXP = G.RAW_RUNS / "experiments" / "wave-climate"
OUT = EXP / "2026-09-27-wave-recommendation" / "figures"
K = ["hs", "wave_period_s", "wave_asymmetry", "wave_angle_high_fraction"]
LABEL = {"hs": "Wave height Hs (m)", "wave_period_s": "Wave period Tp (s)",
         "wave_asymmetry": "Asymmetry", "wave_angle_high_fraction": "High-angle fraction"}
REC = {"hs": 2.0, "wave_period_s": 7.5, "wave_asymmetry": 0.6, "wave_angle_high_fraction": 0.5}
ENDS = {1996: (4.84, 17.55), 2010: (18.8, 24.535)}
COLOR = {"natural": "#1b7f6b", "full_management": "#b4501a"}
NAME = {"natural": "Natural", "full_management": "Full management"}
PER = {1996: "1996–2010", 2010: "2010–2024"}


def scored(path):
    t = pd.read_csv(path)
    return t[t.status == "scored"] if "status" in t else t


def zerobe():
    t = scored(EXP / "2026-09-25-wave-grid-smoothed-score" / "tables" / "all_runs.csv")
    return t.drop_duplicates(["scenario", "period_start", *K])


def final():
    t = scored(EXP / "2026-09-27-wave-grid-fixed-ends" / "tables" / "all_runs.csv")
    return t.drop_duplicates(["scenario", "period_start", *K])


def observed_change_smoothed(period):
    """Observed end-minus-start change, smoothed as the target is
    (common.smooth_like_target: LOWESS at common.SMOOTH_DOMAINS, the southern
    10 raw). The 5-scr table carries 0/3/5/10 only, so 7 is built from its raw
    (window 0) column (Hannah, 2026-09-28: the group's range is 7)."""
    f = (p2.INIT_ROOT / "5-scr" / "3-rates" / "coastsat" / "total_change"
         / p2.study.window(period) / "smoothed" / "tables" / "domain_smoothed.csv")
    d = pd.read_csv(f)
    return common.smooth_like_target(d[d.window_domains == 0].set_index("domain_number")["observed_m"])


def draw_shoals(ax, label=True, label_pt=STRUCTURE_LABEL_PT):
    """Shoal zones as faint hatched boxes, as coastsat_lrr_windows.draw_shoals
    (5-scr/3-rates) draws them, named at the bottom when `label`."""
    matplotlib.rcParams["hatch.linewidth"] = 0.5
    for name, (lo, hi) in HATTERAS_ANNOTATIONS.shoal_zones.items():
        kw = dict(transform=ax.get_xaxis_transform(), facecolor="none", zorder=0.8, clip_on=True)
        ax.add_patch(plt.Rectangle((lo - 0.5, 0), hi - lo + 1, 1, hatch="///",
                                   edgecolor=C["ADDED"], lw=0, alpha=0.30, **kw))
        ax.add_patch(plt.Rectangle((lo - 0.5, 0), hi - lo + 1, 1, edgecolor=C["ADDED"],
                                   lw=0.6, alpha=0.55, **kw))
        if label:
            ax.text((lo + hi) / 2, 0.02, name, transform=ax.get_xaxis_transform(),
                    ha="center", va="bottom", fontsize=label_pt, color="#8a620e", zorder=1)


def draw_fills(ax, start, end, label_pt=STRUCTURE_LABEL_PT):
    """A bar over each enabled model-input fill in the window, the year on it,
    just inside the top of the panel (the title sits above the frame)."""
    trans = ax.get_xaxis_transform()
    for p in HATTERAS_NOURISHMENT_PROJECTS:
        if not (p.enabled and start <= p.year <= end):
            continue
        lo, hi = min(p.gis_domains), max(p.gis_domains)
        ax.plot([lo - 0.45, hi + 0.45], [0.90, 0.90], color=INK, lw=2.2,
                solid_capstyle="butt", zorder=7, transform=trans)
        ax.text((lo + hi) / 2, 0.875, f"{p.year} fill", ha="center", va="top",
                fontsize=label_pt, color=INK, zorder=7, transform=trans)


def rec_row(t, sc, p):
    x = t[(t.scenario == sc) & (t.period_start == p)
          & np.logical_and.reduce([np.isclose(t[k], v) for k, v in REC.items()])]
    return None if x.empty else x.iloc[0]


# ---------------------------------------------------------------------------
def fig1():
    # coarse phase only: the one full factorial. Refine values (Hs 1.25, Tp 7.5,
    # ...) were run next to a few settings only, so their "best" is not
    # comparable (it drew a false dip at Tp 7.5).
    z = zerobe()
    z = z[z.phase == "coarse"]
    f, axes = plt.subplots(2, 4, figsize=(17, 8.5), constrained_layout=True)
    rows = []
    for i, p in enumerate((1996, 2010)):
        for j, k in enumerate(K):
            ax = axes[i, j]
            for sc in COLOR:
                x = z[(z.scenario == sc) & (z.period_start == p)]
                best = x.groupby(k).raw_variance_explained.max() * 100
                ax.plot(best.index, best.values, "o-", color=COLOR[sc], lw=2, ms=6, label=NAME[sc])
                for v, s in best.items():
                    rows.append(dict(period=PER[p], scenario=sc, parameter=k, value=v, best_raw_pct=s))
            ax.axvline(REC[k], color=INK, ls=":", lw=1.2)
            ax.axhline(0, color=INK_MUTED, lw=0.6)
            ax.grid(True)
            open_frame(ax)
            if i == 1:
                ax.set_xlabel(LABEL[k])
            if j == 0:
                ax.set_ylabel(f"{PER[p]}\nbest share explained (%)")
            _title(ax, 4 * i + j, "")
    axes[0, 0].legend(frameon=False)
    png = OUT / "1_score_by_parameter.png"
    save(f, png, dpi=250, close=True)
    pd.DataFrame(rows).to_csv(support_dir(png.parent) / f"{png.stem}.csv", index=False)
    record_caption(png, (
        "For each value of each wave parameter, the best raw share of the alongshore "
        "variation explained by any run with that value, the other three free, in the "
        "zeroBE coarse grid (2026-09-25; the only search that ran every combination, "
        "so every value is compared on equal terms; refine values left out). Top 1996-2010, bottom 2010-2024; green natural, orange full "
        "management. Dotted: the recommended value. 0% is a flat line at the observed "
        "mean. A flat curve means the parameter does not constrain the fit; a peak means it does."))
    return png


# ---------------------------------------------------------------------------
def archived_1996():
    rows = []
    tg = common.coastsat_target(1996)
    for md in (G.RAW_RUNS / "archive" / "2026-09-27-fixed-ends-1996-ends-solved-at-hs1" / "runs").glob(
            "*/1996_2010/edgeBE/*/*_run_metadata.json"):
        w = json.loads(md.read_text(encoding="utf-8"))["wave climate"]
        rows.append(dict(scenario=md.parts[-5].split("_", 1)[1], period_start=1996,
                         hs=float(w["wave_height_m"]), wave_period_s=float(w["wave_period_s"]),
                         wave_asymmetry=float(w["wave_asymmetry"]),
                         wave_angle_high_fraction=float(w["wave_angle_high_frac"]),
                         raw_variance_explained=G.score_run(md.parent, tg)["raw_variance_explained"]))
    return pd.DataFrame(rows)


def fig2():
    tests = [
        ("Step 2: one parameter at a time\n(ends zero, 09-24)",
         scored(EXP / "2026-09-25-wave-grid-smoothed-score" / "tables" / "step2_rescored_smoothed.csv")),
        ("Four-parameter grid\n(ends zero, 09-25)", zerobe()),
        ("Shortlist, ends solved\nfor each setting (09-26)",
         scored(EXP / "2026-09-26-wave-shortlist-ends-solved" / "tables" / "all_runs.csv")),
        ("Grid, ends fixed\n(solved at Hs 1, 09-27)", archived_1996()),
        ("Targeted, ends fixed\n(solved at adopted waves)", final()),
    ]
    f, axes = plt.subplots(1, 4, figsize=(17, 6.2), sharey=True, constrained_layout=True)
    rows = []
    for y, (name, t) in enumerate(tests):
        for sc, mk in (("full_management", "o"), ("natural", "^")):
            x = t[(t.scenario == sc) & (t.period_start == 1996)]
            if x.empty:
                continue
            top = x.nlargest(5, "raw_variance_explained")
            b = top.iloc[0]
            rows.append(dict(test=name.replace("\n", " "), scenario=sc, **{k: b[k] for k in K},
                             raw_pct=100 * b.raw_variance_explained))
            for j, k in enumerate(K):
                off = -0.12 if sc == "full_management" else 0.12
                axes[j].scatter(top[k], np.full(len(top), y + off), s=18, color=COLOR[sc], alpha=0.35)
                axes[j].scatter([b[k]], [y + off], s=120, marker=mk, color=COLOR[sc], edgecolor=INK, zorder=5)
    for j, k in enumerate(K):
        ax = axes[j]
        ax.axvline(REC[k], color=INK, ls=":", lw=1.2)
        ax.set_xlabel(LABEL[k])
        ax.grid(axis="x")
        open_frame(ax)
        _title(ax, j, "")
    axes[0].set_yticks(range(len(tests)))
    axes[0].set_yticklabels([n for n, _ in tests])
    axes[0].invert_yaxis()
    handles = [Line2D([], [], marker="o", ls="", ms=10, color=COLOR["full_management"], mec=INK, label="Full management, best"),
               Line2D([], [], marker="^", ls="", ms=10, color=COLOR["natural"], mec=INK, label="Natural, best"),
               Line2D([], [], marker="o", ls="", ms=5, color=INK_MUTED, alpha=0.5, label="next four"),
               Line2D([], [], ls=":", color=INK, label="recommended")]
    f.legend(handles=handles, loc="outside lower center", ncol=4, frameon=False)
    png = OUT / "2_agreement_across_tests.png"
    save(f, png, dpi=250, close=True)
    pd.DataFrame(rows).to_csv(support_dir(png.parent) / f"{png.stem}.csv", index=False)
    record_caption(png, (
        "The best 1996-2010 setting (large markers; small: the next four) of each search, "
        "by parameter, on the raw score. Each search treated the end domains differently "
        "(zero, solved per setting, fixed at values solved at Hs 1, fixed at values solved at "
        "the adopted waves) and covered a different set of settings. Dotted: the recommended "
        "value. The high-angle fraction and the period settle early; Hs rises and asymmetry "
        "falls once sand is supplied at the ends."))
    return png


# ---------------------------------------------------------------------------
def fig3(runs=None, ends=None, png=None, solved_on="the CoastSat LRR"):
    """`runs` {(period, scenario): run folder} and `ends` {period: (GIS 1, GIS 90)}
    draw another end solve (2026-09-28: the position-change solve); by default
    the fixed-ends sweep's option-A runs and ENDS."""
    t = final()
    ends = ends or ENDS
    png = png or OUT / "3_recommended_vs_coastsat.png"
    # Target, observed change and the header scores at common.SMOOTH_DOMAINS
    # (7 since 2026-09-28); the table's scores were made at 10, so re-scored here.
    targets = {p: common.coastsat_target(p) for p in (1996, 2010)}
    obs = {p: observed_change_smoothed(p) for p in (1996, 2010)}
    f, axes = plt.subplots(2, 2, figsize=(17, 10.5), sharex=True, sharey="row", constrained_layout=True)
    rows, head = [], {}
    for j, p in enumerate((1996, 2010)):
        ax_r, ax_p = axes[0, j], axes[1, j]
        ax_r.plot(targets[p].index, targets[p].values, color=INK, lw=2.8, zorder=6)
        ax_p.plot(obs[p].index, obs[p].values, color=INK, lw=2.8, zorder=6)
        parts, parts_p = [], []
        for sc in COLOR:
            if runs is None:
                r = rec_row(t, sc, p)
                if r is None:
                    continue
                run = EXP / "2026-09-27-wave-grid-fixed-ends" / r.run_dir
            else:
                run = runs.get((p, sc))
                if run is None:
                    continue
            rt = pd.read_csv(run / "tables" / "shoreline_change_rate.csv").set_index("gis_domain")
            sc7 = G.score_run(run, targets[p])
            ax_r.plot(rt.index, rt.lrr_m_yr, color=COLOR[sc], lw=1.8, zorder=5)
            ax_p.plot(rt.index, rt.change_rate_m_yr * 14, color=COLOR[sc], lw=1.8, zorder=5)
            name = "managed" if sc == "full_management" else "natural"
            parts.append(f"{name} {sc7['raw_rmse_m_yr']:.2f}")
            d = common.interior(rt.change_rate_m_yr * 14) - common.interior(obs[p])
            rmse_p = float(np.sqrt((d ** 2).mean()))
            parts_p.append(f"{name} {rmse_p:.1f}")
            rows.append(dict(period=PER[p], scenario=sc, lowess_domains=common.SMOOTH_DOMAINS,
                             raw_pct=100 * sc7["raw_variance_explained"],
                             smoothed_pct=100 * sc7["smoothed_variance_explained"],
                             rmse=sc7["raw_rmse_m_yr"], bias=sc7["bias_m_yr"],
                             position_rmse_m=rmse_p, run_dir=str(run)))
        for ax in (ax_r, ax_p):
            ax.axhline(0, color=INK_MUTED, lw=0.6)
            ax.set_xlim(1, 90)
            ax.grid(axis="y")
            open_frame(ax)
            town_bands(ax, label=(ax is ax_r), fontsize=11)
            draw_shoals(ax, label=(ax is ax_r), label_pt=11)
        structures(ax_p, label=True, label_pt=11)
        structures(ax_r, label=False)
        end = p + 14
        draw_fills(ax_r, p, end, label_pt=11)
        draw_fills(ax_p, p, end, label_pt=11)
        _title(ax_r, j, f"{PER[p]} (boundary flux {ends[p][0]:+.1f} / {ends[p][1]:+.1f} m/yr)\n"
                        "RMSE (m/yr): " + ", ".join(parts))
        _title(ax_p, 2 + j, "RMSE (m): " + ", ".join(parts_p))
        ax_p.set_xlabel(DOMAIN_AXIS_LABEL)
    axes[0, 0].set_ylabel("LRR (m/yr)")
    axes[1, 0].set_ylabel("Position change (m)")
    handles = [Line2D([], [], color=INK, lw=2.8, label=f"CoastSat (LOWESS, {common.SMOOTH_DOMAINS} domains)"),
               Line2D([], [], color=COLOR["full_management"], lw=1.8, label="Model, full management"),
               Line2D([], [], color=COLOR["natural"], lw=1.8, label="Model, natural")]
    f.legend(handles=handles, loc="outside lower center", ncol=3, frameon=False,
             title="Hs = 2.0 m, Tp = 7.5 s, asymmetry 0.6, high-angle fraction 0.5")
    save(f, png, dpi=250, close=True)
    pd.DataFrame(rows).to_csv(support_dir(png.parent) / f"{png.stem}.csv", index=False)
    record_caption(png, (
        "The recommended wave climate (Hs 2.0 m, Tp 7.5 s, asymmetry 0.6, high-angle 0.5) "
        f"with the end domains fixed at the values solved for it on {solved_on} (full "
        "management), against CoastSat (black). Top: LRR rate per domain; bottom: position "
        "change, end minus start. Headers: boundary flux at GIS 1 / GIS 90 and the RMSE of "
        "the raw model LRR (top) and position change (bottom), interior GIS 2-89, "
        f"against the target at LOWESS {common.SMOOTH_DOMAINS} domains (southern 10 "
        "raw). Edge source/sink (edgeBE: GIS 1 and 90 only), metres offset (dune line), no "
        "groin, no relocations; the groin and piers are marked, not modelled. Hatched: shoal "
        "zones. Bars: model-input fills in the window (none 1996-2010)."))
    return png


# ---------------------------------------------------------------------------
def fig4():
    # The step-2 high-angle sweep: the one clean series through 0.5 (0.1-0.55,
    # Hs 1.0, Tp 8, asym 0.8, ends zero), so every point differs in the
    # high-angle fraction only.
    s2 = G.step2
    t = pd.read_csv(s2.TABLES_DIR / "all_runs.csv")
    t = t[(t.status == "scored") & (t.period_start == 1996)
          & t.group.isin(["high_angle", "full_management_high_angle", "baseline", "baseline_full_management"])
          & np.isclose(t.hs, 1.0) & np.isclose(t.wave_period_s, 8.0) & np.isclose(t.wave_asymmetry, 0.8)]
    t = t.drop_duplicates(["scenario", "wave_angle_high_fraction"])
    tg = common.coastsat_target(1996)
    rows = []
    for _, r in t.iterrows():
        d = s2.STUDY_DIR / r.run_dir
        v = common.run_rates(d).loc[2:89].values
        sc = G.score_run(d, tg)
        rows.append(dict(scenario=r.scenario, ha=r.wave_angle_high_fraction, jump=np.std(np.diff(v)),
                         spread=np.std(v), raw=100 * sc["raw_variance_explained"],
                         smoothed=100 * sc["smoothed_variance_explained"], r=sc["raw_r"]))
    df = pd.DataFrame(rows).sort_values("ha")
    tgt = tg.loc[2:89].values
    f, axes = plt.subplots(1, 3, figsize=(17, 5.6), constrained_layout=True)
    for sc in COLOR:
        x = df[df.scenario == sc]
        axes[0].plot(x.ha, x.jump, "o-", color=COLOR[sc], lw=2, label=NAME[sc])
        axes[1].plot(x.ha, x.spread, "o-", color=COLOR[sc], lw=2, label=NAME[sc])
        axes[2].plot(x.ha, x.raw, "o-", color=COLOR[sc], lw=2, label=f"{NAME[sc]}, raw")
        axes[2].plot(x.ha, x.smoothed, "o--", color=COLOR[sc], lw=1.4, alpha=0.6, label=f"{NAME[sc]}, smoothed")
    axes[0].axhline(np.std(np.diff(tgt)), color=INK, lw=1.5, ls="--", label="CoastSat target")
    axes[1].axhline(np.std(tgt), color=INK, lw=1.5, ls="--", label="CoastSat target")
    axes[0].set_ylabel("Domain-to-domain jump,\nsd of first difference (m/yr)")
    axes[1].set_ylabel("Alongshore spread of the rate,\nsd over GIS 2-89 (m/yr)")
    axes[2].set_ylabel("Share explained (%)")
    axes[2].set_ylim(-60, 30)
    for k, ax in enumerate(axes):
        ax.axvline(REC["wave_angle_high_fraction"], color=INK, ls=":", lw=1.2)
        ax.axhline(0, color=INK_MUTED, lw=0.6) if k == 2 else None
        ax.set_xlabel(LABEL["wave_angle_high_fraction"])
        ax.grid(True)
        open_frame(ax)
        _title(ax, k, "")
        ax.legend(frameon=False, fontsize=10)
    png = OUT / "4_high_angle_roughness.png"
    save(f, png, dpi=250, close=True)
    df.to_csv(support_dir(png.parent) / f"{png.stem}.csv", index=False)
    record_caption(png, (
        "The high-angle fraction alone, 1996-2010 (step-2 sweep: Hs 1.0 m, Tp 8 s, "
        "asymmetry 0.8, ends zero; every other setting fixed). (a) Domain-to-domain jump "
        "in the modelled rate and (b) its alongshore spread, against the CoastSat target's "
        "(dashed). (c) Raw (solid) and smoothed (dashed) share explained, axis clipped at "
        "-60%. The high-angle fraction sets how strongly the model turns bends in the "
        "measured planform into rate differences: at low fractions the profile is jagged "
        "and swings too widely (spread 2-2.5 m/yr against 1.17 observed); by 0.45-0.5 "
        "the jaggedness is down to CoastSat's but the spread has fallen below half of it. "
        "The score peaks where those two errors balance: 0.45 here (Hs 1, ends zero), "
        "0.5 at Hs 2 with the ends fixed (figure 2)."))
    return png


# ---------------------------------------------------------------------------
def fig5():
    from cascade_pipeline.shoreline import compute_lrr
    from site_layer.hatteras_site_config import HATTERAS_DOMAINS as D
    t = final()
    real = slice(D.start_real_index, D.end_real_index)
    gis = np.arange(D.first_gis_id, D.first_gis_id + D.num_real_domains)
    wins = [(1996, 2010, "1996–2010"), (2010, 2024, "2010–2024"), (2010, 2020, "2010–2020")]
    rows = []
    for sc in COLOR:
        for start, end, lab in wins:
            r = rec_row(t, sc, start)
            if r is None:
                continue
            d = EXP / "2026-09-27-wave-grid-fixed-ends" / r.run_dir
            m = np.load(next(d.glob("*_shoreline_matrix.npy")))
            n = end - start + 1
            lrr, _ = compute_lrr(m[:n], span_years=end - start)
            rates = pd.Series(lrr[real], index=gis)
            tg = common.coastsat_target(start, None if end == start + 14 else end)
            s = common.alongshore_scores(rates, tg)
            mi, oi = common.interior(rates), common.interior(tg)
            rows.append(dict(scenario=sc, window=lab, model_mean=float(mi.mean()), observed_mean=float(oi.mean()),
                             raw_pct=100 * s["variance_explained"], pattern_pct=100 * s["pattern_variance_explained"],
                             r=s["r_alongshore"]))
    df = pd.DataFrame(rows)
    f, axes = plt.subplots(1, 2, figsize=(15, 5.8), constrained_layout=True)
    xs = np.arange(len(wins))
    w = 0.25
    obs = df.drop_duplicates("window").set_index("window").observed_mean.reindex([l for *_, l in wins])
    axes[0].bar(xs - w, obs.values, w, color=INK, label="CoastSat")
    for k, sc in enumerate(COLOR):
        x = df[df.scenario == sc].set_index("window").reindex([l for *_, l in wins])
        axes[0].bar(xs + k * w, x.model_mean.values, w, color=COLOR[sc], label=NAME[sc])
        axes[1].bar(xs + (k - 0.5) * w, x.pattern_pct.values, w, color=COLOR[sc], label=NAME[sc])
    axes[0].set_ylabel("Interior mean rate, GIS 2-89 (m/yr)")
    axes[1].set_ylabel("Pattern-only share explained (%),\nmean bias removed")
    for k, ax in enumerate(axes):
        ax.set_xticks(xs)
        ax.set_xticklabels([l for *_, l in wins])
        ax.axhline(0, color=INK_MUTED, lw=0.8)
        ax.grid(axis="y")
        open_frame(ax)
        _title(ax, k, "")
        ax.legend(frameon=False)
    png = OUT / "5_window_2010.png"
    save(f, png, dpi=250, close=True)
    df.to_csv(support_dir(png.parent) / f"{png.stem}.csv", index=False)
    record_caption(png, (
        "The recommended setting in three windows. (a) Interior mean rate, model against "
        "CoastSat: 2010-2024 CoastSat is lifted to accretion by the island-wide +17 m step "
        "into 2021; on 2010-2020 the observed mean is near zero. (b) The share of the "
        "alongshore pattern explained with each series' mean removed. 2010-2020 is scored on "
        "the model's first 11 annual states against the CoastSat 2010-2020 LRR target."))
    return png


# ---------------------------------------------------------------------------
def fig6():
    """2010-2024: the same setting as 1996-2010 against the one allowed change,
    Hs 2.0 -> 2.5, each on the ends solved for it (2026-09-27, Hannah)."""
    t = scored(EXP / "2026-09-27-wave-grid-fixed-ends" / "tables" / "all_runs.csv")
    t = t[t.period_start == 2010]
    same = t[(t.phase != "final") & np.logical_and.reduce([np.isclose(t[k], v) for k, v in REC.items()])]
    hs25 = t[(t.phase == "final") & np.isclose(t.hs, 2.5) & np.isclose(t.wave_period_s, 7.5)
             & np.isclose(t.wave_asymmetry, 0.6) & np.isclose(t.wave_angle_high_fraction, 0.5)]
    target, obs = common.coastsat_target(2010), p2.observed_change(2010)
    lighter = {"natural": "#8cc7b9", "full_management": "#e0a67f"}
    f, axes = plt.subplots(2, 2, figsize=(17, 10.5), sharex=True, constrained_layout=True)
    rows = []
    for j, sc in enumerate(("full_management", "natural")):
        ax_r, ax_p = axes[0, j], axes[1, j]
        ax_r.plot(target.index, target.values, color=INK, lw=2.8, zorder=6)
        ax_p.plot(obs.index, obs.values, color=INK, lw=2.8, zorder=6)
        head = []
        for lab, x, col, lw, ends in (("same as 1996-2010, Hs 2.0", same, lighter[sc], 2.2, (18.8, 24.5)),
                                      ("Hs 2.5", hs25, COLOR[sc], 2.0, (8.0, 40.4))):
            r = x[x.scenario == sc].iloc[0]
            rt = pd.read_csv(EXP / "2026-09-27-wave-grid-fixed-ends" / r.run_dir / "tables"
                             / "shoreline_change_rate.csv").set_index("gis_domain")
            ax_r.plot(rt.index, rt.lrr_m_yr, color=col, lw=lw, zorder=5 if lab == "Hs 2.5" else 4)
            ax_p.plot(rt.index, rt.change_rate_m_yr * 14, color=col, lw=lw, zorder=5 if lab == "Hs 2.5" else 4)
            head.append(f"{lab}: {100 * r.raw_variance_explained:+.0f}%, RMSE {r.raw_rmse_m_yr:.2f}, "
                        f"bias {r.bias_m_yr:+.2f}")
            rows.append(dict(scenario=sc, setting=lab, raw_pct=100 * r.raw_variance_explained,
                             rmse=r.raw_rmse_m_yr, bias=r.bias_m_yr, gis1_end=ends[0], gis90_end=ends[1],
                             run_dir=r.run_dir))
        for ax in (ax_r, ax_p):
            ax.axhline(0, color=INK_MUTED, lw=0.6)
            ax.set_xlim(1, 90)
            ax.grid(axis="y")
            open_frame(ax)
            town_bands(ax, label=(ax is ax_r), fontsize=11)
        structures(ax_p, label=True, label_pt=11)
        structures(ax_r, label=False)
        _title(ax_r, j, f"{NAME[sc]}, 2010–2024\n" + "\n".join(head))
        _title(ax_p, 2 + j, "")
        ax_p.set_xlabel(DOMAIN_AXIS_LABEL)
    axes[0, 0].set_ylabel("Shoreline change rate, LRR (m/yr)")
    axes[1, 0].set_ylabel("Position change, 2024 minus 2010 (m)")
    handles = [Line2D([], [], color=INK, lw=2.8, label="CoastSat (LOWESS, 10 domains)"),
               Line2D([], [], color="#9a9a9a", lw=2.2, label="Same setting as 1996–2010 (Hs 2.0; ends +18.8 / +24.5)"),
               Line2D([], [], color="#3a3a3a", lw=2.0, label="Hs raised to 2.5 (ends +8.0 / +40.4)")]
    f.legend(handles=handles, loc="outside lower center", ncol=3, frameon=False,
             title="2010–2024, Tp 7.5 s, asymmetry 0.6, high-angle 0.5 (light: same setting; dark: Hs 2.5)")
    png = OUT / "6_hs_change_2010.png"
    save(f, png, dpi=250, close=True)
    pd.DataFrame(rows).to_csv(support_dir(png.parent) / f"{png.stem}.csv", index=False)
    record_caption(png, (
        "2010-2024 under the one allowed change between windows: the same wave climate as "
        "1996-2010 (Hs 2.0 m, light) against Hs raised to 2.5 m (dark), Tp 7.5 s, asymmetry 0.6 "
        "and high-angle 0.5 in both, each with the end rates solved for it (full management, "
        "against CoastSat; GIS 1 / GIS 90 +18.8 / +24.5 and +8.0 / +40.4 m/yr). Left full "
        "management, right natural; top LRR rate, bottom position change, against CoastSat "
        "(black). Header: raw share of the alongshore variation explained, RMSE and mean "
        "interior bias (m/yr); the flat-line RMSE for this window is 1.46 m/yr."))
    return png


def main():
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
    apply_style()
    OUT.mkdir(parents=True, exist_ok=True)
    with plt.rc_context(p2.SCREEN_RC):
        for fn in (fig1, fig2, fig3, fig4, fig5, fig6):
            print(fn().relative_to(EXP))


if __name__ == "__main__":
    main()
