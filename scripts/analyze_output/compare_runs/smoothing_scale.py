"""
smoothing_scale.py
==============================================================================
The modelled net change in shoreline position against the change projected
from the CoastSat LRR, with the RATE smoothed at four widths before it is
projected. Built 2026-09-21 (Hannah, by interview, after the same sweep on the
observations alone in 3-rates/coastsat/total_change/<window>/smoothed/).

THE FIGURE IS THE POINT (Hannah, 2026-09-21)
    projected_vs_model_<window>.png: ONE alongshore panel per model period,
    all 90 domains, with the four smoothing widths laid over each other on a
    light-to-dark blue ramp (SMOOTH_RAMP) and the model in black.

    THE BLACK LINE IS THE RUN UNTOUCHED -- nothing about the model responds to
    the smoothing -- so it is the one fixed thing in the panel, and the spread
    of the blue family around it is the whole result.

    The runs are graded against a target that is NOT the raw rate: raw domain
    means over GIS 1-10, a 10-domain LOESS of the transect rates beyond
    (cascade_pipeline.coastsat_loess, via hindcast.build_target_table). The
    darkest curve is that grading window; the palest is no smoothing at all.
    Over GIS 1-10 all four curves coincide, because the splice keeps the raw
    domain means there whatever the window.

    The widths are 0 (raw), 3, 5 and 10 domains (0, 1.5, 2.5, 5.0 km), the
    same four the projected-LRR product uses.

THE TABLE, BEHIND THE FIGURE
    tables/skill_by_window.csv scores every combination in two forms:
    as_graded      the target smoothed, the model raw. What the runner
                   actually does, and what the figures draw. It is skill.csv's
                   coastsat_loess generalised to four widths.
    scale_matched  target and model both smoothed at the same width -- the
                   only form in which the two sides are treated alike. Kept
                   for the record; it is not drawn.

THE NULL, AND WHY IT IS NOT OPTIONAL HERE
    Interior r is 0.05-0.33. That is the regime where a symmetric smoother
    inflates correlation hardest: it strips high-frequency variance that the
    two sides do not share, so r rises and RMSE falls at every window whether
    or not the model is any good. Without a baseline the sweep draws a curve
    that looks like "the model improves at coarser scales" and means nothing.

    So every r is reported beside r_null_p95: the 95th percentile of r over
    N_NULL phase-randomised surrogates of the SAME model series -- same mean,
    same variance, same alongshore autocorrelation, no relation to the target
    -- each put through the identical smoothing and splice. An r above that
    band is skill the smoother cannot manufacture. An r inside it is not.

    The bias is close to smoothing-invariant and needs no null; it is the one
    number in the table that a wider window cannot flatter.

ONE ASYMMETRY, ON THE RECORD
    The target's LOESS is fitted at TRANSECT resolution (~906 points) and then
    averaged to domains. The model exists only at 90 domains, so smoothing it
    means lowess over those 90 values at the same physical width (frac =
    window / n). Same width, coarser resolution. That is why as_graded is
    reported too: it involves no model-side smoothing at all.

SCOPE   the full-period CoastSat target (target_comparison/projected/, the target
        in use, the 1996-2024 LRR x 14 yr in both windows) and the dune line
        beside it, against all three model sets. ends_unsolved is the headline
        -- the zeroBE arm carries no source/sink term in any domain, so
        neither target was fitted anywhere in it and all 90 domains are held
        out. The two solved sets are swept too, which answers whether the edge
        solve still buys anything once the grading is done at 5 km.
        Interior GIS 2-89 throughout, as the run index scores.

OUTPUT  output/comparisons/target_comparison/smoothing_scale/
    projected_vs_model_<window>.png     THE FIGURE: one alongshore panel, the
                                        four smoothing widths on a light-to-
                                        dark ramp over a fixed model line;
                                        PDF and caption under supporting/
    tables/skill_by_window.csv          every width x model set x target x
                                        form: n, bias, RMSE, r, r_null_p95
    tables/domain_values_<window>.csv   the per-domain projection at every
                                        width and each model set's net change
    runs_used.csv, PROVENANCE.md

USAGE
    python scripts/analyze_output/compare_runs/smoothing_scale.py
    python ... --windows 0 3 5 10 --n-null 1000
==============================================================================
"""
from __future__ import annotations

import argparse
import datetime as dt
import sys
from pathlib import Path

import numpy as np
import pandas as pd

_HERE = Path(__file__).resolve()
sys.path.insert(0, str(_HERE.parent))
import rate_windows as rw  # noqa: E402
import target_comparison as tc  # noqa: E402

_REPO = rw._REPO
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402
from matplotlib.ticker import MultipleLocator  # noqa: E402
from statsmodels.nonparametric.smoothers_lowess import lowess  # noqa: E402

from cascade_pipeline.coastsat_loess import spliced_loess_series  # noqa: E402
from site_layer.hat_observed_rates import dune_endpoint_csv, lrr_csv  # noqa: E402
from site_layer.hatteras_site_config import HATTERAS_DOMAINS as DOM  # noqa: E402
from site_layer.hat_figure_style import (  # noqa: E402
    INK, SMOOTH_RAMP, apply_style, caption, figsize, save,
)

OUT_DIR = tc.ROOT_DIR / "smoothing_scale"
SMOOTH_WINDOWS = (0, 3, 5, 10)      # domain units; 10 is the grading window
N_NULL = 1000
SEED = 20260921
FORMS = ("as_graded", "scale_matched")
# Every LRR window on disk, for the structure diagnostic. The target is the
# full-period 1996-2024; the others are there because the same LOESS removes
# very different amounts from them, and 1984-2004 is the field the
# method-comparison figure draws (Hannah, 2026-09-21).
LRR_WINDOWS = ((1984, 2004), (2004, 2024), (1996, 2024), (1996, 2010), (2010, 2024))


# -----------------------------------------------------------------------------
# the targets, at an arbitrary window
# -----------------------------------------------------------------------------
def _along(df, domain_col, order_col):
    """Per-transect along-coast distance in metres, the convention every
    target build here uses: each domain's transects spread evenly across its
    500 m band, ordered within the domain by `order_col`."""
    df = df.sort_values([domain_col, order_col]).reset_index(drop=True)
    rank = df.groupby(domain_col).cumcount()
    n = df.groupby(domain_col)[domain_col].transform("count")
    sp = DOM.domain_spacing_m
    along = ((df[domain_col] - DOM.first_gis_id) * sp
             + (rank + 0.5) * (sp / n)).to_numpy(float)
    return df, along


def coastsat_transects(window):
    """(domain ids, along-coast m, rate m/yr) for the full-period LRR."""
    t = pd.read_csv(lrr_csv(*window))
    t = t[t["domain_number"].between(DOM.first_gis_id, DOM.last_gis_id)].copy()
    t, along = _along(t, "domain_number", "transect_id")
    return t["domain_number"].to_numpy(int), along, t["lrr_m_yr"].to_numpy(float)


def duneline_transects(window):
    """(domain ids, along-coast m, rate m/yr) for the two-survey dune rate,
    read and ordered exactly as rw.load_dune_endpoint_target does."""
    t = pd.read_csv(dune_endpoint_csv(*window, "transect"))
    t, along = _along(t, "domain_number", "line_id")
    return t["domain_number"].to_numpy(int), along, t["rate_m_yr"].to_numpy(float)


def field_structure(window, windows):
    """How much alongshore structure a LOESS of each width takes out of one
    LRR field, and whether there is independent error for it to average.

    Reported as the SD removed IN m/yr, not as a share of variance: the
    domain-mean variance is dominated by the long-wavelength swings, so a
    wiggle that is plainly visible on the figure reads as a few per cent of it
    and the percentage badly undersells the effect (Hannah caught this
    2026-09-21, comparing against the 1984-2004 panels of
    input_prep/6-scr-smooth/loess_method_comparison.py).

    Returns a dict, or None when the window's LRR has not been built.
    """
    path = lrr_csv(*window)
    if not path.exists():
        return None
    t = pd.read_csv(path)
    t = t[t["domain_number"].between(DOM.first_gis_id, DOM.last_gis_id)]
    dom_ids, along, rate = coastsat_transects(window)
    dmean = t.groupby("domain_number")["lrr_m_yr"].mean()
    within = float((t["lrr_m_yr"] - t["domain_number"].map(dmean)).std(ddof=1))

    idx = pd.RangeIndex(DOM.first_gis_id, DOM.last_gis_id + 1)
    y = dmean.reindex(idx).to_numpy(float)
    y = y - np.nanmean(y)
    ac = np.correlate(y, y, "full")[len(y) - 1:]
    ac = ac / ac[0]
    decorr = next((k for k in range(1, len(ac)) if ac[k] < 0.5), len(ac))

    out = dict(window="{}_{}".format(*window), years=window[1] - window[0],
               n_transects=int(len(t)),
               sd_domain_mean_m_yr=round(float(dmean.std(ddof=1)), 4),
               sd_within_domain_m_yr=round(within, 4),
               median_transect_unc_m_yr=round(float(t["unc_m_yr"].median()), 4),
               decorrelation_domains=int(decorr),
               decorrelation_km=round(float(decorr * DOM.domain_spacing_m / 1000.0), 2))
    for w in windows:
        if not w:
            continue
        sm, _ = spliced_loess_series(dom_ids, along, rate, w, skip=rw.SKIP, domains=DOM)
        resid = dmean.reindex(sm.index) - sm
        out[f"sd_removed_w{w:02d}_m_yr"] = round(float(resid.std(ddof=1)), 4)
    for k in (1, 2, 3, 5, 10):
        out[f"autocorr_lag{k:02d}"] = round(float(ac[k]), 3)
    return out


def target_structure(windows, lrr_windows=None):
    """field_structure over every LRR window on disk, so the grading window's
    effect on the TARGET can be read against the other windows -- in
    particular the 1984-2004 field the method-comparison figure draws, where
    the same LOESS removes roughly twice as much."""
    lrr_windows = lrr_windows or LRR_WINDOWS
    rows = [r for r in (field_structure(w, windows) for w in lrr_windows) if r]
    return pd.DataFrame(rows)


def smooth_domain_series(series, window, skip=rw.SKIP):
    """The model's analogue of the target's pass: lowess over the 90 per-domain
    values at the same physical width, with the same GIS 1..skip splice. The
    model has no transect resolution to smooth at -- see the asymmetry note in
    the module docstring."""
    if not window:
        return series.copy()
    y = series.to_numpy(float)
    x = series.index.to_numpy(float)
    ok = np.isfinite(y)
    if ok.sum() < 5:
        return series.copy()
    frac = float(np.clip(window / ok.sum(), 0.02, 1.0))
    res = lowess(y[ok], x[ok], frac=frac, return_sorted=True)
    out = series.copy()
    out.loc[series.index[ok]] = np.interp(x[ok], res[:, 0], res[:, 1])
    keep = series.index[series.index <= skip]
    out.loc[keep] = series.loc[keep]      # splice, as the target has
    return out


def phase_randomise(y, rng):
    """A surrogate with y's mean, variance and alongshore autocorrelation but
    randomised phases, so it carries no relation to the target. The amplitude
    spectrum is kept and only the phases are redrawn."""
    n = len(y)
    mu = y.mean()
    f = np.fft.rfft(y - mu)
    ph = rng.uniform(0.0, 2.0 * np.pi, f.shape)
    ph[0] = 0.0
    if n % 2 == 0:
        ph[-1] = 0.0
    return np.fft.irfft(np.abs(f) * np.exp(1j * ph), n=n) + mu


# -----------------------------------------------------------------------------
def build(observations, models, windows, n_null):
    rng = np.random.default_rng(SEED)
    lo, hi = rw.INTERIOR
    idx = pd.RangeIndex(DOM.first_gis_id, DOM.last_gis_id + 1, name="domain_number")
    interior = (idx >= lo) & (idx <= hi)

    # The CoastSat target is the FULL-PERIOD LRR, the same rate in both model
    # windows (tc.CS_MODE == "projected"), so it is built once.
    cs_parts = coastsat_transects(tc.FULL_WINDOW)
    targets = {}    # (model window, target name, smoothing window) -> Series in m
    for o in observations:
        years = o.window[1] - o.window[0]
        du_parts = duneline_transects(o.window)
        for w in windows:
            cs, _ = spliced_loess_series(*(cs_parts[0], cs_parts[1], cs_parts[2]),
                                         window=w, skip=rw.SKIP, domains=DOM)
            du, _ = spliced_loess_series(*(du_parts[0], du_parts[1], du_parts[2]),
                                         window=w, skip=rw.SKIP, domains=DOM)
            targets[(o.window, "coastsat", w)] = cs.reindex(idx) * years
            targets[(o.window, "duneline", w)] = du.reindex(idx) * years

    rows, values = [], {o.window: pd.DataFrame(index=idx) for o in observations}
    for o in observations:
        years = o.window[1] - o.window[0]
        for w in windows:
            for name in ("coastsat", "duneline"):
                values[o.window][f"target_{name}_w{w:02d}_m"] = targets[(o.window, name, w)]
        for key, folder in tc.MODEL_SETS.items():
            mdf = models[key][0][o.window]
            if mdf is None:
                continue
            raw = pd.Series(mdf["change_rate_m_yr"].to_numpy(float) * years, index=idx)
            for w in windows:
                sm = smooth_domain_series(raw, w)
                values[o.window][f"model_{folder}_w{w:02d}_m"] = sm
                for form in FORMS:
                    m = raw if form == "as_graded" else sm
                    # A raw model is the same series at every window, so the
                    # as_graded null is redrawn per window only because the
                    # TARGET it is scored against changed.
                    for name in ("coastsat", "duneline"):
                        t = targets[(o.window, name, w)]
                        ok = (m.notna() & t.notna() & interior).to_numpy()
                        a, b = m.to_numpy()[ok], t.to_numpy()[ok]
                        d = a - b
                        null = np.empty(n_null)
                        for i in range(n_null):
                            sur = pd.Series(phase_randomise(raw.to_numpy(float), rng),
                                            index=idx)
                            if form == "scale_matched":
                                sur = smooth_domain_series(sur, w)
                            null[i] = np.corrcoef(sur.to_numpy()[ok], b)[0, 1]
                        rows.append(dict(
                            window="{}_{}".format(*o.window), model_ends=folder,
                            target=name, form=form, smoothing_domains=w,
                            smoothing_km=w * DOM.domain_spacing_m / 1000.0,
                            n=int(ok.sum()), bias_m=round(float(d.mean()), 3),
                            rmse_m=round(float(np.sqrt((d ** 2).mean())), 3),
                            r=round(float(np.corrcoef(a, b)[0, 1]), 3),
                            r_null_p95=round(float(np.percentile(null, 95)), 3),
                            r_null_sd=round(float(null.std(ddof=1)), 3)))
    skill = pd.DataFrame(rows)
    skill["r_clears_null"] = skill["r"] > skill["r_null_p95"]
    return skill, values


# -----------------------------------------------------------------------------
def win_label(w):
    """One smoothing width in words, for a legend entry or a column name."""
    return ("unsmoothed rate" if not w else
            f"rate smoothed {w * DOM.domain_spacing_m / 1000.0:g} km"
            + (" (as graded)" if w == rw.TARGET_WINDOW else ""))


def figure(values, skill, window, windows, model_key=tc.UNSOLVED):
    """The alongshore picture (Hannah, 2026-09-21): the model's own net change
    against the projected change from the LRR, ONE PANEL per model period with
    every smoothing width laid over it on a light-to-dark ramp.

    The model line is the run untouched, so it is the one fixed thing in the
    panel: the spread of the blue family around it is the whole result."""
    folder = tc.MODEL_SETS[model_key]
    df = values[window].reset_index()
    x = df["domain_number"].to_numpy(float)
    model = df[f"model_{folder}_w00_m"].to_numpy(float)     # raw: the run as it is
    half = tc.Y_HALF_M
    sk = skill[(skill.window == "{}_{}".format(*window)) & (skill.model_ends == folder)
               & (skill.target == "coastsat") & (skill.form == "as_graded")
               ].set_index("smoothing_domains")

    fig, ax = plt.subplots(constrained_layout=True,
                           figsize=figsize("double", aspect=0.46))
    # Quantity, window, method (Hannah, 2026-09-21). Projected, not total: the
    # target here is the 1996-2024 LRR carried onto a 14-yr window.
    ax.set_title(f"Projected shoreline change vs CASCADE, {window[0]}–{window[1]} "
                 f"(CoastSat LRR 1996–2024 × {window[1] - window[0]} yr, "
                 "every LOESS width)")
    # mean_lrr all-NaN draws the frame, grid and village bands with no sign
    # fill -- four curves share this panel, so the blue/red fill is not
    # available here (the idiom is target_comparison.draw).
    rw.obs.draw_panel(ax, df.assign(mean_lrr=np.nan, std_lrr=0.0), half,
                      label=True, std=False)
    rw.obs.draw_shoals(ax, label=True)
    fills = rw.obs.fills_in(*window)
    if fills:
        rw.obs.draw_fills(ax, fills, half)
    for i, w in enumerate(windows):
        ax.plot(x, df[f"target_coastsat_w{w:02d}_m"].to_numpy(float),
                color=SMOOTH_RAMP[i % len(SMOOTH_RAMP)], lw=1.2, zorder=10 + i,
                solid_capstyle="round")
    ax.plot(x, model, color=INK, lw=tc.LW_MODEL, zorder=20, solid_capstyle="round")
    ax.yaxis.set_major_locator(MultipleLocator(20.0 if half > 60 else 10.0))
    ax.set_xlabel(rw.DOMAIN_AXIS_LABEL)
    ax.set_ylabel(tc.Y_LABEL)
    h = [Line2D([], [], color=SMOOTH_RAMP[i % len(SMOOTH_RAMP)], lw=1.2)
         for i, _ in enumerate(windows)] + [Line2D([], [], color=INK, lw=tc.LW_MODEL)]
    labels = ([f"Projection, {win_label(w)}" for w in windows]
              + [f"CASCADE net change, {tc.MODEL_LABEL[model_key]}"])
    fig.legend(handles=h, labels=labels, loc="outside lower center", ncol=3,
               frameon=False)
    scores = "; ".join(
        "{} {:+.1f} m bias, {:.1f} m RMSE".format(
            "unsmoothed" if not w else f"{w * DOM.domain_spacing_m / 1000.0:g} km",
            sk.loc[w, "bias_m"], sk.loc[w, "rmse_m"]) for w in windows)
    caption(fig, (
        "{}â{}: the CASCADE run's own net change in shoreline position (black) against "
        "the change projected from the CoastSat linear regression rate, by GIS domain "
        "(1 at Cape Point, 90 at Pea Island), in metres over the {}-yr model window, "
        "seaward positive. The RATE is smoothed alongshore BEFORE it is projected, and "
        "the four blue curves are the four widths on a light-to-dark ramp: palest is "
        "the unsmoothed domain means, then {}, darkest is the {}-domain "
        "({:g} km) window the runs are actually graded at. Every curve keeps the raw "
        "domain means over GIS 1â{} â the boundary treatment at Oregon Inlet â so the "
        "four are identical there by construction, and the projection is that rate Ã "
        "{} yr. THE BLACK LINE IS THE SAME RUN THROUGHOUT: nothing about the model "
        "responds to the smoothing, so the spread of the blue family around a fixed "
        "black line is the whole result. {} carries no source/sink term in any domain, "
        "the two ends included, so all 90 domains are the model's own response and the "
        "target was not fitted anywhere in it. Interior GIS {}â{}, model minus the "
        "projection: {}. The bias barely moves because a symmetric smoother preserves "
        "the mean; the RMSE falls because the smoother removes the transect-scale "
        "scatter the model was never going to reproduce, not because the model "
        "improved. The y axis (Â±{:g} m) is the same on every figure in "
        "target_comparison.{} The other two model sets, the dune-line target and the "
        "correlations with their null bands are in tables/skill_by_window.csv."
    ).format(
        window[0], window[1], window[1] - window[0],
        ", ".join(f"{w * DOM.domain_spacing_m / 1000.0:g} km" for w in windows[1:-1]),
        windows[-1], windows[-1] * DOM.domain_spacing_m / 1000.0,
        rw.SKIP, window[1] - window[0],
        tc.MODEL_CLAUSE[model_key][0].upper() + tc.MODEL_CLAUSE[model_key][1:],
        rw.INTERIOR[0], rw.INTERIOR[1], scores, half,
        tc.over_note([(df, [f"target_coastsat_w{w:02d}_m" for w in windows]
                       + [f"model_{folder}_w00_m"], "")], half)))
    out = save(fig, OUT_DIR / "projected_vs_model_{}_{}".format(*window))
    plt.close(fig)
    return out

def structure_section(struct, windows):
    """The part of PROVENANCE.md that says what the grading window is FOR."""
    tgt = struct[struct.window == "{}_{}".format(*tc.FULL_WINDOW)].iloc[0]
    cols = [w for w in windows if w]
    head = ("| LRR window | yr | domain-mean SD | within-domain SD | "
            + " | ".join(f"SD removed, {w} dom" for w in cols) + " | decorrelates |")
    tab = [head, "|" + "---|" * (5 + len(cols))]
    for _, r in struct.iterrows():
        mark = " **(the target)**" if r.window == tgt.window else ""
        tab.append(
            f"| {r.window.replace('_', '–')}{mark} | {int(r.years)} "
            f"| {r.sd_domain_mean_m_yr:.3f} | {r.sd_within_domain_m_yr:.3f} | "
            + " | ".join(f"{r[f'sd_removed_w{w:02d}_m_yr']:.3f}" for w in cols)
            + f" | {r.decorrelation_km:.1f} km |")
    w10 = tgt[f"sd_removed_w{rw.TARGET_WINDOW:02d}_m_yr"]
    return [
        "## What the 10-domain window is actually doing",
        "",
        "Decided 2026-09-21 (Hannah): **the window stays at "
        f"{rw.TARGET_WINDOW} domains, but it is not noise removal and should not be "
        "described as such.** Every run and both edge solves were graded at "
        f"{rw.TARGET_WINDOW}; nothing is re-solved. What changes is the claim.",
        "",
        "All figures in m/yr of alongshore structure removed, NOT as a share of "
        "variance: the domain-mean variance is dominated by the long-wavelength swings, "
        "so a wiggle that is plainly visible on a figure reads as a few per cent of it "
        "and the percentage badly undersells the effect.",
        "",
        *tab,
        "",
        "### Why the target moves so little under smoothing",
        "",
        "**The 1996-2024 LRR is the smoothest of the five fields**, because it is the "
        "longest fit: each transect's OLS slope is the best constrained and there is "
        f"least scatter to take out. At 5 domains the LOESS removes "
        f"{tgt['sd_removed_w05_m_yr']:.3f} m/yr from it against "
        f"{struct[struct.window == '1984_2004'].iloc[0]['sd_removed_w05_m_yr']:.3f} "
        "m/yr from 1984-2004 -- less than half. `projected_vs_model_<window>.png` then "
        "projects over 14 yr rather than 20, which halves the apparent difference "
        "again. So the near-overlapping curves there and the obvious smoothing in "
        "`input_prep/6-scr-smooth/loess_method_comparison.py` are the SAME LOESS on "
        "different fields, not a difference in method.",
        "",
        "### It is still not denoising, for this target",
        "",
        f"- The domain-mean rate decorrelates alongshore at {tgt.decorrelation_km:.1f} km "
        "in every window (r below 0.5 at lag "
        f"{int(tgt.decorrelation_domains)} domains: "
        + ", ".join(f"{k * DOM.domain_spacing_m / 1000:g} km r={tgt[f'autocorr_lag{k:02d}']:+.2f}"
                    for k in (1, 2, 3, 5, 10)) + "). A "
        f"{rw.TARGET_WINDOW}-domain window is "
        f"{rw.TARGET_WINDOW * DOM.domain_spacing_m / 1000.0 / tgt.decorrelation_km:.0f}x "
        "that, so it is well past averaging noise and is cutting into coherent "
        "structure.",
        f"- Scatter between transects inside one domain is only "
        f"{tgt.sd_within_domain_m_yr:.3f} m/yr for this target, and it is LARGER than "
        f"the median per-transect LRR uncertainty ({tgt.median_transect_unc_m_yr:.3f} "
        "m/yr), so even that is not clean estimation noise: it is real sub-500 m "
        "structure, optimistic OLS uncertainties (satellite series violate the "
        "independent-residual assumption), or both.",
        f"- `projected_vs_model_<window>.png` shows where the {w10:.3f} m/yr goes: the "
        "Avon peak (GIS 29-35) and the Wimble peak (GIS 65-72), both shoal-fronted, "
        "both physically anchored features rather than scatter.",
        "",
        "**So the honest statement is: the target is graded at a scale deliberately "
        "coarser than the model resolves, and grading at that scale removes the two "
        "shoal-fronted peaks from the target.** Those are also two of the places the "
        "model misses worst, which is why RMSE improves with the window.",
        "",
        "**This claim is about the FULL-PERIOD target only.** A sub-period LRR carries "
        "roughly twice the estimation noise (within-domain SD "
        f"{struct[struct.window == '1996_2010'].iloc[0]['sd_within_domain_m_yr']:.3f} "
        f"m/yr for 1996-2010 against {tgt.sd_within_domain_m_yr:.3f} here), so for "
        "those targets the LOESS is doing real denoising. Do not carry the statement "
        "across to them. The numbers are in `tables/target_structure.csv`.",
        "",
    ]


def provenance(skill, windows, n_null, runs_used, structure=None):
    def table(form, target):
        s = skill[(skill.form == form) & (skill.target == target)]
        head = ("| window | model set | LOESS | bias (m) | RMSE (m) | r | null r (p95) "
                "| clears null |")
        out = [head, "|" + "---|" * 7]
        for _, x in s.sort_values(["window", "model_ends", "smoothing_domains"]).iterrows():
            out.append(
                f"| {x.window.replace('_', '–')} | {x.model_ends} "
                f"| {'raw' if not x.smoothing_domains else f'{x.smoothing_km:g} km'} "
                f"| {x.bias_m:+.1f} | {x.rmse_m:.1f} | {x.r:.3f} | {x.r_null_p95:.3f} "
                f"| {'yes' if x.r_clears_null else 'no'} |")
        return out

    (OUT_DIR / "PROVENANCE.md").write_text("\n".join([
        "# target_comparison/smoothing_scale - provenance",
        "",
        f"Written {dt.datetime.now():%Y-%m-%d %H:%M} by "
        "scripts/analyze_output/compare_runs/smoothing_scale.py "
        "(2026-09-21, Hannah, by interview).",
        "",
        "## The figure",
        "",
        "`projected_vs_model_<window>.png` is the point: the model's own net change in "
        "shoreline position against the change projected from the CoastSat LRR, ONE "
        "alongshore panel per model period with every smoothing width laid over it on "
        "a light-to-dark blue ramp (" + ", ".join(
            "palest = raw" if not w
            else (f"darkest = {w} domains / {w * DOM.domain_spacing_m / 1000:g} km"
                  if w == windows[-1]
                  else f"{w} domains / {w * DOM.domain_spacing_m / 1000:g} km")
            for w in windows) + "). The rate is smoothed BEFORE it is projected; the "
        f"darkest curve is the {rw.TARGET_WINDOW}-domain window the runs are actually "
        f"graded at. Over GIS 1-{rw.SKIP} all four coincide, because the splice keeps "
        "the raw domain means there whatever the window.",
        "",
        "**The black line is the run untouched, the same in every figure.** Nothing "
        "about the model responds to the smoothing, so the spread of the blue family "
        "around a fixed black line is the whole result. Read it that way: the question "
        "is not whether the model got better, it is how much of the target the model "
        "is being asked to match is real structure and how much is transect-scale "
        "scatter that smoothing removes.",
        "",
        *(structure or []),
        "## The question behind it",
        "",
        "The runs are graded against a 10-domain LOESS of the transect rates (raw "
        f"domain means over GIS 1-{rw.SKIP}). That window has never been examined, and "
        "`../projected/tables/skill.csv` shows it is load-bearing. The "
        "table here sweeps it against the full-period CoastSat target (the 1996-2024 "
        "LRR x 14 yr, the same rate in both windows) and the dune line.",
        "",
        "Model sets and runs are `../projected/runs_used.csv`; the "
        "loaders are `target_comparison.load_model_sets`, so these are the same "
        "runs that folder draws. Interior GIS "
        f"{rw.INTERIOR[0]}-{rw.INTERIOR[1]}, as the run index scores.",
        "",
        "## Read the bias, and r only against its own null",
        "",
        f"Interior r here is small. A symmetric smoother strips high-frequency "
        "variance the two sides do not share, so r rises and RMSE falls at EVERY "
        "window whether or not the model is any good. Every r therefore carries "
        f"`r_null_p95`: the 95th percentile of r over {n_null} phase-randomised "
        "surrogates of that same model series -- same mean, variance and alongshore "
        "autocorrelation, no relation to the target -- each put through the identical "
        "smoothing and splice. Only `r_clears_null` is evidence. The bias is close to "
        "smoothing-invariant and needs no null.",
        "",
        "RMSE has no null column. It falls with the window for the same mechanical "
        "reason and is only comparable BETWEEN model sets at one window, never across "
        "windows as a measure of improvement.",
        "",
        "## Two forms in the table",
        "",
        "- **as_graded** -- target smoothed, model raw: what the runner does, and what "
        "the figures draw.",
        "- **scale_matched** -- both smoothed at the same width, kept for the record "
        "and not drawn. The target's LOESS is "
        "fitted at transect resolution (~906 points) and averaged to domains; the "
        "model exists only at 90 domains, so its pass is lowess over those 90 values "
        "at the same physical width. Same width, coarser grid -- which is why "
        "as_graded, which smooths no model at all, is reported beside it.",
        "",
        "## CoastSat target, as graded",
        "",
        *table("as_graded", "coastsat"),
        "",
        "## CoastSat target, scale matched",
        "",
        *table("scale_matched", "coastsat"),
        "",
        "The dune-line rows are in `tables/skill_by_window.csv`; the figures draw the "
        "CoastSat target only, since the window is a property of that target's build.",
        "",
        "## Runs",
        "",
        "Full provenance, with the topography version and git commit of each, is in "
        "`runs_used.csv`.",
        "",
        "| window | model set | run | arm |",
        "|---|---|---|---|",
        *[f"| {r['window'].replace('_', '–')} | {tc.MODEL_SETS[r['model_ends']]} "
          f"| `{r['run_name']}` | `{r['arm']}` |" for r in runs_used],
        "",
    ]), encoding="utf-8")


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description="Does the grading smoothing window matter?")
    ap.add_argument("--windows", nargs="+", type=int, default=list(SMOOTH_WINDOWS),
                    metavar="N", help="LOESS widths in domain units; 0 = raw.")
    ap.add_argument("--n-null", type=int, default=N_NULL)
    a = ap.parse_args(argv)
    windows = sorted(set(a.windows))

    # The canonical name, not the "full" alias: tc branches on CS_MODE by
    # equality in several places, and an alias would quietly take the wrong arm.
    tc.CS_MODE = "projected"
    tc.OUT_DIR = tc.ROOT_DIR / tc.CS_MODES[tc.CS_MODE]
    apply_style()
    observations = [rw.Observation(w) for w in tc.WINDOWS]
    models = tc.load_model_sets()

    # tc.over_note names every value that runs off the axis; give it words for
    # the columns this script invents.
    tc._COL_NAME.update(
        {f"target_coastsat_w{w:02d}_m": f"projection, {win_label(w)}" for w in windows}
        | {f"model_{f}_w00_m": f"model, {tc.MODEL_LABEL[k]}"
           for k, f in tc.MODEL_SETS.items()})

    skill, values = build(observations, models, windows, a.n_null)
    (OUT_DIR / "tables").mkdir(parents=True, exist_ok=True)
    skill.to_csv(OUT_DIR / "tables" / "skill_by_window.csv", index=False)
    for w, df in values.items():
        df.round(3).to_csv(OUT_DIR / "tables" / "domain_values_{}_{}.csv".format(*w))

    written = []
    for o in observations:
        written += figure(values, skill, o.window, windows)
    runs = [dict(r, model_ends=k) for k, (_, rows) in models.items() for r in rows]
    runs.sort(key=lambda r: (r["window"], r["model_ends"]))
    pd.DataFrame([dict(r, model_ends=tc.MODEL_SETS[r["model_ends"]]) for r in runs]
                 ).to_csv(OUT_DIR / "runs_used.csv", index=False)
    struct = target_structure(windows)
    struct.to_csv(OUT_DIR / "tables" / "target_structure.csv", index=False)
    provenance(skill, windows, a.n_null, runs,
               structure=structure_section(struct, windows))
    cols = ["window", "years", "sd_domain_mean_m_yr", "sd_within_domain_m_yr"] + \
           [f"sd_removed_w{w:02d}_m_yr" for w in windows if w] + ["decorrelation_km"]
    print("\nalongshore structure of each LRR field (m/yr); the target is "
          "{}_{}".format(*tc.FULL_WINDOW))
    print(struct[cols].to_string(index=False))

    show = skill[(skill.target == "coastsat")
                 & (skill.model_ends == tc.MODEL_SETS[tc.UNSOLVED])]
    print(show[["window", "form", "smoothing_km", "bias_m", "rmse_m", "r",
                "r_null_p95", "r_clears_null"]].to_string(index=False))
    for p in written:
        print("wrote   ", Path(p).relative_to(_REPO))
    return 0


if __name__ == "__main__":
    sys.exit(main())
