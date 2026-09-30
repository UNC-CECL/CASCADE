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
    structures, support_dir, town_bands)

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
def fig3():
    t = final()
    targets = {p: common.coastsat_target(p) for p in (1996, 2010)}
    obs = {p: p2.observed_change(p) for p in (1996, 2010)}
    f, axes = plt.subplots(2, 2, figsize=(17, 10.5), sharex=True, constrained_layout=True)
    rows, head = [], {}
    for j, p in enumerate((1996, 2010)):
        ax_r, ax_p = axes[0, j], axes[1, j]
        ax_r.plot(targets[p].index, targets[p].values, color=INK, lw=2.8, zorder=6)
        ax_p.plot(obs[p].index, obs[p].values, color=INK, lw=2.8, zorder=6)
        parts = []
        for sc in COLOR:
            r = rec_row(t, sc, p)
            if r is None:
                continue
            rt = pd.read_csv(EXP / "2026-09-27-wave-grid-fixed-ends" / r.run_dir / "tables"
                             / "shoreline_change_rate.csv").set_index("gis_domain")
            ax_r.plot(rt.index, rt.lrr_m_yr, color=COLOR[sc], lw=1.8, zorder=5)
            ax_p.plot(rt.index, rt.change_rate_m_yr * 14, color=COLOR[sc], lw=1.8, zorder=5)
            parts.append(f"{NAME[sc].split()[0].lower()} {100 * r.raw_variance_explained:+.0f}%, "
                         f"bias {r.bias_m_yr:+.2f}")
            rows.append(dict(period=PER[p], scenario=sc, raw_pct=100 * r.raw_variance_explained,
                             smoothed_pct=100 * r.smoothed_variance_explained, bias=r.bias_m_yr,
                             run_dir=r.run_dir))
        for ax in (ax_r, ax_p):
            ax.axhline(0, color=INK_MUTED, lw=0.6)
            ax.set_xlim(1, 90)
            ax.grid(axis="y")
            open_frame(ax)
            town_bands(ax, label=(ax is ax_r), fontsize=11)
        structures(ax_p, label=True, label_pt=11)
        structures(ax_r, label=False)
        _title(ax_r, j, f"{PER[p]}: ends GIS 1 {ENDS[p][0]:+.1f}, GIS 90 {ENDS[p][1]:+.1f} m/yr\n"
                        + ";  ".join(parts))
        _title(ax_p, 2 + j, "")
        ax_p.set_xlabel(DOMAIN_AXIS_LABEL)
    axes[0, 0].set_ylabel("Shoreline change rate, LRR (m/yr)")
    axes[1, 0].set_ylabel("Position change, end minus start (m)")
    handles = [Line2D([], [], color=INK, lw=2.8, label="CoastSat (LOWESS, 10 domains)"),
               Line2D([], [], color=COLOR["full_management"], lw=1.8, label="Model, full management"),
               Line2D([], [], color=COLOR["natural"], lw=1.8, label="Model, natural")]
    f.legend(handles=handles, loc="outside lower center", ncol=3, frameon=False,
             title="Recommended: Hs 2.0 m, Tp 7.5 s, asymmetry 0.6, high-angle 0.5")
    png = OUT / "3_recommended_vs_coastsat.png"
    save(f, png, dpi=250, close=True)
    pd.DataFrame(rows).to_csv(support_dir(png.parent) / f"{png.stem}.csv", index=False)
    record_caption(png, (
        "The recommended wave climate (Hs 2.0 m, Tp 7.5 s, asymmetry 0.6, high-angle 0.5) "
        "with the end domains fixed at the values solved for it (full management), against "
        "CoastSat (black). Top: LRR rate per domain; bottom: position change, end minus "
        "start. Header: raw share of the alongshore variation explained and mean interior "
        "bias (m/yr). Metres offset (dune line), no groin, no relocations."))
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
