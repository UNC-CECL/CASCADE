#!/usr/bin/env python3
"""
The corrected-groin hindcast in four figures: the strength over time, the gap across the groins,
how the strength was chosen, and the test period run by run.

    python results_figures.py   ->  figures/corrected_groin_<n>_<what>.png

Reads the runs scored by corrected_groin.py (grid_scores.csv) and their shoreline matrices, the
pinned groin's runs from ../2026-10-08-schedule-refit/, and the observed gap (annual CoastSat and
the wet/dry photo dates). Captions go to figures/supporting/CAPTIONS.md. Details: README.md.

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-10-09
"""
from __future__ import annotations

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402

import corrected_groin as cgm  # noqa: E402
import drift_groin as dg  # noqa: E402
import schedule_refit as sr  # noqa: E402
from score_instant_grid import observed_series  # noqa: E402
from site_layer.hat_figure_style import (  # noqa: E402
    C, INK, INK_MUTED, _letter_inside, apply_style, figsize, open_frame, record_caption, save)

# --- CONFIG ------------------------------------------------------------------
HERE = cgm.HERE
FIGS = HERE / "figures"
CAL, TEST = cgm.CALIBRATION, cgm.TEST
CHOSEN, JOINT, PHOTO = 0.12, 0.19, 0.25               # the strengths each calibration target picks
OBS_C = C["LOCATOR"]                                  # observations: teal, markers only
RUNS = {                                              # label: (colour, line style)
    "no groin": (C["BASE"], "-"),
    "pinned groin (in the runner now)": (C["ACCENT"], "-"),
    f"corrected, {CHOSEN} to 2016": (INK, "-"),
    f"corrected, {JOINT} to 2016": (INK, (0, (3, 2))),
}
EVENTS = ((1969, "built"), (1995, "last repair"), (2003, "Isabel damage"), (2017, "2017 fill"),
          (2022, "2022 fill"), (2026, "south groin\nrebuilt"))
# -----------------------------------------------------------------------------


def two_number(s):
    return f"1969:{s},2017:0"


# The two-number run at strength s; on 1996-2009 schedule A's constant-strength run is the same run
def two_number_file(period, s):
    p = cgm.shoreline_file(cgm.tag("periods", two_number(s)), period, s, 1.0)
    return p if p is not None else cgm.shoreline_file(cgm.tag("instant1996"), period, s, 1.0)


# The shoreline file for each run in a period
def run_files(period):
    return {
        "no groin": sr.baseline_file(period),
        "pinned groin (in the runner now)": sr.shoreline_file("instant2004", period, 0.6, 0.6),
        f"corrected, {CHOSEN} to 2016": two_number_file(period, CHOSEN),
        f"corrected, {JOINT} to 2016": two_number_file(period, JOINT),
    }


# Figure 1: the corrected groin's strength over time, the options, against the structure's history
def strength_over_time():
    fig, ax = plt.subplots(figsize=figsize("double", height=3.3), constrained_layout=True)
    yrs = np.arange(1965, 2028)
    options = (                                       # drawn in this order: the grey line sits under the orange one
        ("2004 failure only (0.2, then 0.16), never off", "1969:0.2,2004:0.16", C["BASE"], "-", 3.2),
        ("three numbers: 0.2, 0.16 from 2004, off from 2017", "1969:0.2,2004:0.16,2017:0", C["ADDED"], "-", 1.4),
        (f"two numbers: {JOINT} (fitted to both records), off from 2017", two_number(JOINT), INK, (0, (3, 2)), 1.4),
        (f"two numbers: {CHOSEN} (fitted to annual CoastSat), off from 2017", two_number(CHOSEN), INK, "-", 2.0),
    )
    for label, spec, col, ls, lw in options:
        ax.step(yrs, [cgm.strength_for(spec, y) if y >= 1969 else 0 for y in yrs], where="post", color=col, ls=ls,
                lw=lw, label=label)
    ax.axvspan(2015, 2021, color=C["BASE_FILL"], alpha=0.6, lw=0)
    ax.text(2019.3, 0.03, "decay\nreported", ha="center", va="center", fontsize=6.5, color=INK_MUTED)
    for y, name in EVENTS:
        ax.axvline(y, color=INK_MUTED, lw=0.6, ls=(0, (1, 2)))
        ax.text(y + 0.3, 0.315, name, rotation=90, fontsize=6.3, color=INK_MUTED, va="top", ha="left")
    for (y0, y1), name in (((1996, 2009), "calibration"), ((2009, 2025), "test")):
        ax.annotate("", xy=(y0, -0.035), xytext=(y1, -0.035), arrowprops=dict(arrowstyle="<->", color=INK, lw=0.8),
                    annotation_clip=False)
        ax.text((y0 + y1) / 2, -0.045, name, ha="center", va="top", fontsize=7, color=INK)
    ax.set_xlim(1965, 2027.5)
    ax.set_ylim(0, 0.32)
    ax.set_ylabel("Strength: share of the drift\nthe groin holds back")
    ax.grid(axis="y")
    open_frame(ax)
    fig.legend(*ax.get_legend_handles_labels(), loc="outside lower center", ncol=2, fontsize=6.6, frameon=False)
    fig.suptitle("How strong is the groin each year in each option?")
    png = FIGS / "corrected_groin_1_strength_over_time.png"
    save(fig, png, dpi=250, close=True)
    record_caption(png, (
        "**Tests:** nothing; this figure shows the options compared in the corrected-groin hindcast. **How to "
        "read:** each line is the corrected groin's strength, the share of the wave-climate net drift arriving at "
        "the groin that it holds back each year (0.12 at Buxton is about 11 m of shoreline a year). Dotted lines: "
        "the structure's recorded events (`../../../1-observations/structure_history.md`); shading: the "
        "structural losses reported for 2015-2021. Arrows: the calibration (1996-2009) and test (2009-2025) "
        "periods. The pinned groin in the runner now is not drawn: its 0.6 is a share of BRIE's diffusive "
        "exchange, not of the drift, so it is not on this axis. **Shows:** the two-number groin holds one strength "
        "until the 2017 fill and nothing after; its strength is fitted on 1996-2009 and its end year is taken "
        "from the record, not fitted."))


# Figure 2: the gap across the groins, observed and modelled, calibration and test
def gap_over_time():
    photos = observed_series()
    fig, axes = plt.subplots(1, 2, figsize=figsize("double", height=3.4), constrained_layout=True,
                             gridspec_kw=dict(width_ratios=(13, 16)))
    for i, period in enumerate((CAL, TEST)):
        ax = axes[i]
        obs = sr.coastsat_change(period)
        files = run_files(period)
        for label, (col, ls) in RUNS.items():
            m = sr.model_change(files[label])
            ax.plot(period + np.arange(len(m)) + 0.5, m, color=col, ls=ls, lw=1.6 if "corrected" in label else 1.3)
        n = len(sr.model_change(files["no groin"]))
        obs = obs[(obs.index >= period) & (obs.index < period + n)]       # this period's years only
        ax.plot(obs.index + 0.5, obs.values, "o", color=OBS_C, ms=3.5, mfc="white", mew=1.0)
        yrs = [y for y in photos.index if period < y <= period + n]
        g0 = np.interp(period, photos.index, photos.values)
        ax.plot(yrs, [photos[y] - g0 for y in yrs], "s", color=OBS_C, ms=5)
        if period == TEST:
            ax.axvline(2017.5, color=INK_MUTED, lw=0.7, ls=(0, (1, 2)))
            ax.text(2017.7, ax.get_ylim()[1] * 0.95, "2017 fill;\ngroin off",
                    fontsize=6.5, color=INK_MUTED, va="top")
        ax.axhline(0, color=INK_MUTED, lw=0.5)
        ax.set_xlim(period - 0.3, period + n + 0.3)
        ax.set_xlabel("Year")
        ax.grid(axis="y")
        open_frame(ax)
        _letter_inside(ax, i)
        ax.set_title(f"{'Calibration' if period == CAL else 'Test'}, {period}-{sr.window(period)[-4:]}",
                     loc="center", fontsize=8.5)
    axes[0].set_ylabel("Change in the gap across the\ngroins since the start (m)")
    handles = [Line2D([], [], color=col, ls=ls, lw=1.5, label=label) for label, (col, ls) in RUNS.items()]
    handles += [Line2D([], [], color=OBS_C, marker="o", ls="none", mfc="white", label="observed, annual CoastSat"),
                Line2D([], [], color=OBS_C, marker="s", ls="none", label="observed, wet/dry photo dates")]
    fig.legend(handles=handles, loc="outside lower center", ncol=3, fontsize=6.8, frameon=False)
    fig.suptitle("Does the corrected groin reproduce the gap across the Buxton groins?")
    png = FIGS / "corrected_groin_2_gap_over_time.png"
    save(fig, png, dpi=250, close=True)
    sc_ = pd.read_csv(HERE / "grid_scores.csv")
    return png, sc_


# Figure 3: how the strength was chosen, and what each choice gives on the test
def choosing_strength(sc_):
    cal = sc_[(sc_.period == CAL) & (sc_.arm == "corrected")].copy()
    cal = cal[((cal.schedule == "instant1996") & (cal.off.fillna("") == "")) | (cal.schedule == "periods")]
    cal["s"] = [r.b if r.schedule == "instant1996" else cgm.parse_periods(r.off)[0][1] for r in cal.itertuples()]
    cal = cal.sort_values("s")
    test = sc_[(sc_.period == TEST) & (sc_.schedule == "periods")].copy()
    test = test[test.off.isin([two_number(s) for s in (CHOSEN, JOINT, PHOTO)])]
    test["s"] = [cgm.parse_periods(o)[0][1] for o in test.off]
    test = test.sort_values("s")
    none = sc_[(sc_.period == TEST) & (sc_.arm == "none")].iloc[0]
    fig, axes = plt.subplots(1, 2, figsize=figsize("double", height=3.0), constrained_layout=True)
    ax = axes[0]
    ax.plot(cal.s, cal.rmse_annual, color=OBS_C, lw=1.6, marker="o", ms=3, label="annual CoastSat")
    ax.plot(cal.s, cal.rmse_photos, color=OBS_C, lw=1.6, ls=(0, (3, 2)), marker="s", ms=3, label="photo dates")
    for s, name in ((CHOSEN, "annual\nbest"), (JOINT, "both"), (PHOTO, "photos\nbest")):
        ax.axvline(s, color=INK_MUTED, lw=0.6, ls=(0, (1, 2)))
        ax.text(s, 63, name, ha="center", va="top", fontsize=6.3, color=INK_MUTED)
    ax.set_ylim(0, 65)
    ax.set_xlabel("Strength (share of the drift held back)")
    ax.set_ylabel("RMS error in the gap,\n1996-2009 (m)")
    ax.legend(loc="lower left", fontsize=6.8, frameon=False)
    ax.grid(axis="y")
    open_frame(ax)
    _letter_inside(ax, 0)
    ax.set_title("Calibration: each record wants a different strength", loc="center", fontsize=8)
    ax = axes[1]
    x = np.arange(len(test))
    for k, (col, mk, name) in enumerate((("rmse_2009_2016", "o", "2009-2016"), ("rmse_2017_2024", "s", "2017-2024"))):
        ax.plot(x + (k - 0.5) * 0.18, test[col], mk, color=INK, mfc="white" if k else INK, ms=6, label=name,
                ls="none")
        ax.axhline(none[col], color=C["BASE"], lw=1.0, ls="-" if k == 0 else (0, (3, 2)))
    ax.text(len(test) - 0.5, none["rmse_2009_2016"] + 1, "no groin, 2009-2016", fontsize=6.3, color=C["BASE"],
            ha="right", va="bottom")
    ax.text(len(test) - 0.5, none["rmse_2017_2024"] - 1, "no groin, 2017-2024", fontsize=6.3, color=C["BASE"],
            ha="right", va="top")
    ax.set_xticks(x, [f"{s:g}\n({n})" for s, n in zip(test.s, ("annual", "both", "photos"))])
    ax.set_xlim(-0.5, len(test) - 0.5)
    ax.set_ylim(0, None)
    ax.set_xlabel("Strength (calibrated on)")
    ax.set_ylabel("Annual RMS error in the gap,\ntest period (m)")
    ax.legend(loc="center right", fontsize=6.8, frameon=False)
    ax.grid(axis="y")
    open_frame(ax)
    _letter_inside(ax, 1)
    ax.set_title("Test: the annual-fitted strength predicts best", loc="center", fontsize=8)
    fig.suptitle("Which strength, and fitted to which observations?")
    png = FIGS / "corrected_groin_3_choosing_strength.png"
    save(fig, png, dpi=250, close=True)
    t = test.set_index("s")
    record_caption(png, (
        "**Tests:** which strength the two-number corrected groin should have, and which observed record to fit "
        "it to. **How to read:** (a) the root-mean-square error in the gap across the groins over the calibration "
        "period, 1996-2009, against the groin's strength, scored on the annual CoastSat gap (solid) and on the "
        "wet/dry photo dates (dashed); dotted lines mark the strength each target picks and the one that balances "
        "them. (b) Each of those three strengths run on the test period, 2009-2025, with the groin off from the "
        "2017 fill: the annual CoastSat error before (filled circles) and after (open squares) the fill; grey "
        "lines are the run without a groin. Only 2009-2016 played no part in any choice. **Shows:** the annual "
        f"record wants {CHOSEN}, the photos {PHOTO}. On the test the annual-fitted strength is best before the fill "
        f"({t.loc[CHOSEN, 'rmse_2009_2016']:.1f} m, against {t.loc[JOINT, 'rmse_2009_2016']:.1f} m for {JOINT}, "
        f"{t.loc[PHOTO, 'rmse_2009_2016']:.1f} m for {PHOTO} and {none['rmse_2009_2016']:.1f} m with no groin) and "
        f"after it ({t.loc[CHOSEN, 'rmse_2017_2024']:.1f} m). The photo-fitted strength does worse than no groin "
        "before the fill. Its pull toward a stronger groin rests on the 2004 survey, where the photos and "
        "CoastSat disagree by 42 m."))


# Figure 4: every version on the test period, before and after the 2017 fill
def test_summary(sc_):
    t = sc_[sc_.period == TEST]
    pin = sr.shoreline_file("instant2004", TEST, 0.6, 0.6)
    pin_split = cgm.split_scores(pin, TEST, sr.coastsat_change(TEST))

    def row(sel):
        r = t[sel]
        return r.iloc[0] if len(r) else None

    rows = [
        ("no groin", row(t.arm == "none")),
        ("pinned groin (in the runner now)", pd.Series(pin_split)),
        ("drift groin, b 0.1 f 0.4", row((t.arm == "drift") & (t.b == 0.1) & (t.f == 0.4))),
        ("corrected, 0.2 then 0.16, never off", row((t.schedule == "ramp1996") & (t.off.fillna("") == "")
                                                    & (t.b == 0.2) & (t.f == 0.8))),
        ("corrected, declining 2015-2021 to 0", row(t.off == "2015-2021r0")),
        ("corrected, 0.2 then 0.16, off from 2017", row(t.off == "1969:0.2,2004:0.16,2017:0")),
        (f"corrected, {JOINT}, off from 2017", row(t.off == two_number(JOINT))),
        (f"corrected, {CHOSEN}, off from 2017", row(t.off == two_number(CHOSEN))),
    ]
    rows = [(n, r) for n, r in rows if r is not None]
    fig, ax = plt.subplots(figsize=figsize("double", height=3.4), constrained_layout=True)
    yv = np.arange(len(rows))[::-1]
    for y, (name, r) in zip(yv, rows):
        a, b = r["rmse_2009_2016"], r["rmse_2017_2024"]
        ax.plot([a, b], [y, y], color=C["BASE_FILL"], lw=2.5, zorder=1)
        ax.plot(a, y, "o", color=INK, ms=6, zorder=2)
        ax.plot(b, y, "s", color=INK, mfc="white", ms=6, zorder=2)
    nb = rows[0][1]
    ax.axvline(nb["rmse_2009_2016"], color=C["BASE"], lw=0.7, ls="-")
    ax.axvline(nb["rmse_2017_2024"], color=C["BASE"], lw=0.7, ls=(0, (3, 2)))
    ax.set_yticks(yv, [n for n, _ in rows], fontsize=7)
    ax.set_xlim(0, None)
    ax.set_xlabel("Annual RMS error in the gap across the groins, test period (m)")
    ax.grid(axis="x")
    open_frame(ax)
    ax.legend(handles=[Line2D([], [], marker="o", color=INK, ls="none", label="2009-2016"),
                       Line2D([], [], marker="s", color=INK, mfc="white", ls="none", label="2017-2024 (after the fill)"),
                       Line2D([], [], color=C["BASE"], lw=0.7, label="no groin, each half")],
              loc="lower right", fontsize=6.8, frameon=False)
    fig.suptitle("How does each groin do on 2009-2025, before and after the 2017 fill?")
    png = FIGS / "corrected_groin_4_test_summary.png"
    save(fig, png, dpi=250, close=True)
    best = rows[-1][1]
    record_caption(png, (
        "**Tests:** how each version of the groin does on the test period, 2009-2025, which was not used to fit "
        "its strength, split at the 2017 Buxton fill. **How to read:** one row per run; the filled circle is the "
        "annual CoastSat error in the gap across the groins for 2009-2016, the open square for 2017-2024; grey "
        "lines mark the run without a groin. The pinned groin is the one in the runner now (b 0.6, f 0.6); the "
        "drift groin is the corrected groin without its leak correction. **Shows:** every groin that stays on "
        "after 2017 misses 2017-2024 by 25-52 m, because it holds the fill's sand at GIS 6 while the observed sand "
        "passed the groin field; switching the groin off at the fill brings that half to 12-17 m. Before the fill, "
        "the pinned groin does best (7.5 m), then the weaker drift and corrected groins (11.5-11.9 m); the "
        "stronger corrected groins do worse. The two-number groin fitted to annual CoastSat "
        f"({CHOSEN}, off from 2017) is best on both halves of the corrected runs ({best['rmse_2009_2016']:.1f} and "
        f"{best['rmse_2017_2024']:.1f} m). The 2017 switch-off was suggested by the 2017-2024 misfit and is "
        "supported by the placement record, so only 2009-2016 is fully independent (README)."))


def main():
    apply_style()
    FIGS.mkdir(parents=True, exist_ok=True)
    strength_over_time()
    png, sc_ = gap_over_time()
    t = sc_[sc_.period == TEST]
    best = t[t.off == two_number(CHOSEN)].iloc[0]
    none = t[t.arm == "none"].iloc[0]
    record_caption(png, (
        "**Tests:** whether the corrected groin reproduces the observed gap across the Buxton groins, in the "
        "calibration period it was fitted to and in the test period it was not. **How to read:** the change in "
        "the gap between GIS 6 (updrift, north) and GIS 5 (downdrift, south) since each period's start, from its "
        "own starting topography; positive means the gap widened. Teal: observed, annual CoastSat (open circles) "
        "and the wet/dry photo dates (squares). Lines: no groin (grey), the pinned groin in the runner now "
        f"(purple), and the two-number corrected groin at {CHOSEN} (solid black, fitted to annual CoastSat) and "
        f"{JOINT} (dashed, fitted to both records), each off from the 2017 fill. **Shows:** without a groin the "
        "gap closes far too fast in both periods. In the test, the pinned groin follows the observed gap until "
        "2017 and then holds it open while the observed gap closes; the corrected groin, off from the fill, "
        f"follows it through ({best['end_model']:.0f} m at the end against {best['end_obs']:.0f} m observed and "
        f"{none['end_model']:.0f} m with no groin)."))
    choosing_strength(sc_)
    test_summary(sc_)
    print("wrote 4 figures to", FIGS.relative_to(cgm.REPO))


if __name__ == "__main__":
    main()
