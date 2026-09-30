"""
target_comparison.py
==============================================================================
Which observation should CASCADE be graded against? The two candidate
targets, the CoastSat shoreline and the digitized dune line, side by side with
the model, as NET CHANGE IN POSITION over each model window, 1996-2010 and
2010-2024. Built 2026-09-19 (Hannah, by interview).

EVERYTHING OVER THE MODEL PERIOD (14 yr per window)
    CoastSat target   the window's LRR (3-rates/coastsat/lrr/<window>, the
                      runner's scoring series) x 14 yr: the projected change.
    Dune-line target  the MEASURED dune-line net change
                      (3-rates/duneline/endpoint/<window>) projected to the
                      model years: its rate over the survey interval (11.6 yr
                      for 1997-10 to 2009-05, 14.1 yr for 2009-05 to 2023-07)
                      x 14 yr.
    Model             the run's own net change over its 14 years, unchanged
                      (endpoint rate x 14 = last annual shoreline minus first).
    Seaward positive, metres. The gap between the two targets is the
    beach-width change they imply (CoastSat minus dune line): solid grey
    where the beach widened, hatched where it narrowed.

RAW AND SMOOTHED
    The lines are the raw domain means. tables/skill.csv scores the model
    against both targets both raw and with the scoring target's LOWESS
    treatment (raw means D1-10, 7-domain LOWESS beyond), the form the runs are
    graded in.

THREE MODEL SETS, ONE SUBFOLDER EACH (Hannah has not chosen the target)
    ends_solved_on_coastsat/   the matrix edgeBE run, end domains solved
                               against the CoastSat target
    ends_solved_on_duneline/   the 09-18 dune edge-solve run (mean3), end
                               domains solved against the dune line
    ends_unsolved/             (2026-09-21, Hannah, for her advisor) the
                               zeroBE arm of the same matrix cell: NO
                               source/sink term in any domain, the two ends
                               included, so all 90 domains are the model's
                               own response and neither target was fitted.
                               Also unsolved_run_and_targets_<window>.png
                               (and _smoothed): that one run against both
                               targets, the paired form with a single line.
    Each holds target_comparison_1996_2010_2024.png: 1996-2010 above
    2010-2024, both targets and the model on one y axis.
    paired/                    (2026-09-19, Hannah) one figure per window,
                               target_and_own_run_<window>.png: each target with ITS OWN
                               run only: CoastSat with the CoastSat-solved
                               run, the dune line with the dune-solved run.
                               The end values each run carries are in the
                               panel titles; the runs have NO other
                               source/sink term (GIS 2-89 are zero, checked
                               from the index: be_nonzero_domains == 2).

OUTPUT   output/comparisons/target_comparison/
    README.md, runs_used.csv, tables/domain_values_<w>.csv, tables/skill.csv,
    <model set>/target_comparison_1996_2010_2024.png (PDF and CAPTIONS.md
    under supporting/)

AS A RATE: --units rate (2026-09-29, Hannah: "do the same net change
    subfolder for target_comparison", choosing to keep the metres figures and
    add m/yr beside them). Every metres column divided by the 14 model years,
    which undoes the x 14 exactly: the CoastSat LRR, the dune line's MEASURED
    rate over its own survey interval (no scaling to the window), and the
    model's endpoint rate. Same figures, same skill table (bias/RMSE in m/yr),
    written to <projected|total_change>/change_rate/ with "_rate" after the
    target mode in every stem. y axis +/-8 m/yr, as model_vs_observed's rates.

The loaders are rate_windows.py's (imported), so the observations and
runs are exactly the ones model_vs_observed draws as rates.

USAGE
    python scripts/analyze_output/compare_runs/target_comparison.py
    python ... --units rate                              # change_rate/ in m/yr
    python ... --coastsat-target total [--units rate]
==============================================================================
"""
from __future__ import annotations

import math
import sys
from pathlib import Path

import numpy as np
import pandas as pd

_HERE = Path(__file__).resolve()
sys.path.insert(0, str(_HERE.parent))
import rate_windows as rw  # noqa: E402  (loaders, runs, style constants)

_REPO = rw._REPO
sys.path.insert(0, str(_REPO / "scripts" / "input_prep" / "5-scr" / "lib"))
import scr_paths  # noqa: E402,F401  (5-scr sibling modules onto sys.path)
from coastsat_vs_duneline import beach_width_handles, shade_beach_width  # noqa: E402

import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.legend_handler import HandlerTuple  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402
from matplotlib.ticker import MultipleLocator  # noqa: E402
from site_layer.hat_figure_style import (  # noqa: E402
    COMPARISONS_ROOT, DOMAIN_AXIS_LABEL, INK, _title, apply_style, caption,
    mark_offaxis,
    figsize, save,
)

obs = rw.obs
ROOT_DIR = COMPARISONS_ROOT / "target_comparison"
# The CoastSat target (2026-09-19, Hannah), in the vocabulary settled by
# interview 2026-09-21: a rate turned into a distance is named by the window
# it was FITTED on, never by the arithmetic.
#
#   projected/    the 1996-2024 LRR x 14 yr in BOTH windows -- the rate is
#                 carried onto windows it was NOT fitted on, so it is a
#                 PROJECTION. THE TARGET IN USE. Paired with runs whose ends
#                 were solved against it
#                 (experiments/end-domain-boundaries/2026-09-19-end-domains-solved-on-lrr-1996-2024).
#   total_change/ each window's OWN LRR x 14 yr, as the runner grades. The
#                 rate is evaluated over the window it was fitted on, so
#                 nothing is extrapolated. Kept for the record.
#
# Was coastsat_full_period_lrr/ and coastsat_subperiod_lrr/ until 2026-09-21;
# the folders described the fit window but not what was done with it.
CS_MODES = {"projected": "projected", "total": "total_change",
            # the pre-2026-09-21 names, kept so old commands still run
            "full": "projected", "subperiod": "total_change"}
CS_MODE_NOUN = {"projected": "Projected shoreline change",
                "total": "Total shoreline change"}
# The same, in rate units: no distance, so nothing is projected or totalled;
# the rate is named by its fit window alone.
CS_MODE_NOUN_RATE = {"projected": "Long-term shoreline change rate",
                     "total": "Shoreline change rate"}
# The fit window is in the method string, so a figure pulled out of its
# folder still says where the rate came from (Hannah, 2026-09-21). "total"
# has no single fit window across the two panels -- each is its own -- so it
# is filled in per window by cs_method().
CS_MODE_METHOD = {"projected": "CoastSat LRR 1996–2024 × 14 yr",
                  "total": "CoastSat LRR {}–{} × {} yr"}
CS_MODE_METHOD_RATE = {"projected": "CoastSat LRR 1996–2024",
                       "total": "CoastSat LRR {}–{}"}


def cs_noun():
    return (CS_MODE_NOUN if _net() else CS_MODE_NOUN_RATE)[CS_MODE]


def _targets_line(o):
    """"CoastSat: ... · dune line: ... measured, scaled to 14 yr" (2026-09-22).

    The header used to name only the CoastSat side, so the red line's dates
    and its real interval -- which is not 14 yr -- were in the caption alone.
    """
    m = o.meta
    return ("CoastSat target: " + cs_method(o.window) + "   ·   "
            f"dune line: {m['start_date']} → {m['end_date']}"
            + (" (assumed)" if bool(m.get("end_date_assumed")) else "")
            + f",  {float(m['interval_yr']):.1f} yr,  "
            + (f"measured, scaled to {o.window[1] - o.window[0]} yr" if _net()
               else "measured rate"))


def cs_method(window):
    """The method string for one window, with its fit window resolved."""
    m = (CS_MODE_METHOD if _net() else CS_MODE_METHOD_RATE)[CS_MODE]
    return m.format(*window, window[1] - window[0]) if CS_MODE == "total" else m
# Everything downstream branches on the canonical pair, never on the alias.
CS_CANON = {"full": "projected", "subperiod": "total",
            "projected": "projected", "total": "total"}
FULL_WINDOW = (1996, 2024)
# Re-solved under option A on 2026-09-27 (Hannah: redo target_comparison for
# the new wave climate and offset); the /10 solve was
# end-domain-boundaries/2026-09-19-end-domains-solved-on-lrr-1996-2024.
# The adopted model since 2026-09-28 (Barrier3D hatteras/adopted, storms v3_trim24);
# the option A pre-adoption solve was end-domain-boundaries/2026-09-27-ends-solved-on-lrr-1996-2024-option-a.
# 2026-09-29: re-solved after the dune-cap fix; before it,
# end-domain-boundaries/2026-09-28-ends-solved-on-lrr-1996-2024-adopted.
FULL_SOLVE_DIR = rw.RAW_RUNS / "experiments" / "end-domain-boundaries/2026-09-29-ends-solved-on-lrr-1996-2024-dunecap"
CS_MODE = "projected"
OUT_DIR = ROOT_DIR / CS_MODES[CS_MODE]
# "net": metres over the 14-yr window (every figure until 2026-09-29).
# "rate": the same lines in m/yr, under <mode>/change_rate/ (see the docstring).
UNITS = "net"
RATE_SUBDIR = "change_rate"
Y_HALF_RATE = 8.0          # m/yr, the model_vs_observed rate figures' range


def _net():
    return UNITS == "net"


def _u():
    """The unit every number on a figure is in."""
    return "m" if _net() else "m/yr"


def _fmt(v, signed=False):
    """A score in the figure's unit: 0.1 m, or 0.01 m/yr."""
    d = 1 if _net() else 2
    return f"{v:+.{d}f}" if signed else f"{v:.{d}f}"


def _x14(prefix=" multiplied by"):
    """' multiplied by 14 yr' on the metres figures, nothing on the rate ones."""
    return f"{prefix} 14 yr" if _net() else ""


def _stem_mode():
    """The target mode as a filename token, '_rate' added in rate units."""
    return CS_MODES[CS_MODE] + ("" if _net() else "_rate")


def dune_label():
    return ("Total dune line change (measured, scaled to 14 yr)" if _net()
            else "Dune line change rate (measured over the survey interval)")


def _quantity():
    """What the y axis is, in caption words."""
    return ("net change in shoreline position over the 14-yr model window" if _net()
            else "shoreline change rate over the 14-yr model window")


def _model_quantity():
    return ("net change over the window" if _net() else
            "endpoint rate (last annual shoreline minus first, over the 14 run years)")


def cs_label():
    """The CoastSat target named by the window its rate was FITTED on
    (the 2026-09-21 vocabulary): PROJECTED when the 1996-2024 rate is carried
    onto a 14-yr half, TOTAL when each window uses its own."""
    return ("CoastSat target — " + cs_noun().lower()
            + (f" ({cs_method(None)})" if CS_MODE == "projected"
               else " (each window's own CoastSat LRR" + _x14(" ×") + ")"))


def cs_clause():
    return ("the full-period 1996–2024 linear regression rate (the same rate in both "
            "windows)" if CS_MODE == "projected" else
            "the window's own linear regression rate (the runner's scoring series)")


def full_period_runs():
    """window -> (run_name, tag): the converged step of the 1996-2024 CoastSat
    edge solve, from its loop_log.csv."""
    log = pd.read_csv(FULL_SOLVE_DIR / "loop_log.csv")
    runs = {}
    for w in rw.WINDOWS:
        hit = log[(log["window"] == w[0]) & log["run_tag"].notna()]
        if w not in WINDOWS or hit.empty:
            runs[w] = None
            continue
        tag = str(hit.iloc[-1]["run_tag"])
        # The tag is <theme>/<study>/<reading>/step<k> since the 2026-09-25
        # regrouping; strip the whole study tag, not just its first part.
        study_tag = FULL_SOLVE_DIR.relative_to(rw.RAW_RUNS / "experiments").as_posix()
        step = tag[len(study_tag) + 1:] if tag.startswith(study_tag + "/") else tag.split("/", 1)[1]
        d = next((FULL_SOLVE_DIR / step / "{}_{}".format(*w) / "edgeBE").glob("HAT_*"))
        runs[w] = (d.name, tag)
    return runs
WINDOWS = [(1996, 2010), (2010, 2024)]
# The third set (2026-09-21, Hannah, for her advisor): the SAME matrix cell
# with NO source/sink term anywhere, not even at the two ends -- the zeroBE
# arm of the run the CoastSat-solved set reads. All 90 domains are then the
# model's own response, so the ends are readable rather than prescribed.
UNSOLVED = "unsolved"
UNSOLVED_PRESET = "zeroBE"
UNSOLVED_RUNS = {   # the option A matrix since 2026-09-27
    (1996, 2010): ("HAT_1996_2010_zeroBE_offsetmetres_road_bdm_nogroin", "calibration"),
    (2010, 2024): ("HAT_2010_2024_zeroBE_offsetmetres_road_bdm_nourish_nogroin", "calibration"),
}
MODEL_SETS = {"coastsat": "ends_solved_on_coastsat",
              rw.MAIN_DUNE: "ends_solved_on_duneline",
              UNSOLVED: "ends_unsolved"}
MODEL_LABEL = {"coastsat": "ends solved on CoastSat",
               rw.MAIN_DUNE: "ends solved on the dune line",
               UNSOLVED: "ends not solved"}
# How each set's run is named in a caption, after "the model's own net change".
MODEL_CLAUSE = {
    "coastsat": "the edgeBE run with its ends solved on CoastSat",
    rw.MAIN_DUNE: "the edgeBE run with its ends solved on the dune line",
    UNSOLVED: ("the zeroBE run of the same matrix cell, which carries no "
               "source/sink term in any domain including the two ends"),
}
LW = 1.1
# Fixed y range, +/- m, every figure here. 80 from 2026-09-19; 100 from
# 2026-09-22 (Hannah), to match the metre figures in 3-rates and
# 4-comparisons so the three trees can be laid side by side. Anything
# beyond it is named by over_note() and marked by mark_offaxis().
Y_HALF_M = 100.0


def over_note(frames_cols, half):
    """' Beyond ±100 m, off the axis: ...' for the (frame, columns, label)
    triples one figure draws, naming each out-of-range domain; '' if none.

    The canvas marker is `mark_offaxis` at each draw site; this is the words
    that go with it.
    """
    hits = []
    for df, cols, label in frames_cols:
        for col in cols:
            v = df.set_index("domain_number")[col]
            for g, x in v[v.abs() > half].items():
                hits.append(f"{label}{_COL_NAME.get(col, col)} {_fmt(x, True)} "
                            f"{_u()} at GIS {g}")
    return (f" Beyond ±{half:g} {_u()}, off the axis and marked with a triangle at "
            "the edge: " + "; ".join(hits) + "."
            if hits else "")


_COL_NAME = {"coastsat_target_m": "CoastSat target", "dune_target_m": "the total dune line change",
             "model_ends_solved_on_coastsat_m": "CoastSat-solved run",
             "model_ends_solved_on_duneline_m": "dune-solved run",
             "model_ends_unsolved_m": "unsolved run"}
LW_MODEL = 1.4
Y_LABEL = "Net change in position (m)"
Y_LABEL_RATE = "Change rate (m/yr)"


def y_label():
    return Y_LABEL if _net() else Y_LABEL_RATE


def y_tick(half):
    return (20.0 if half > 60 else 10.0) if _net() else 2.0


def to_rate(df):
    """Every metres column over the model years: the x 14 undone. Column
    names are kept so the drawing code reads either frame."""
    out = df.copy()
    for c in [c for c in out.columns if c.endswith("_m")]:
        out[c] = out[c] / out["model_years"]
    return out


def _table_names(df):
    """The rate frame's columns renamed for the CSV: *_m -> *_m_yr."""
    return df.rename(columns={c: c[:-2] + "_m_yr" for c in df.columns if c.endswith("_m")})


def window_values(o, mdfs):
    """Per domain, metres over the model period: both targets raw and LOWESS,
    and each model set's net change."""
    s, e = o.window
    years = e - s
    df = pd.DataFrame({"domain_number": np.arange(1, rw.N + 1)})
    df["model_years"] = years
    df["dune_interval_yr"] = o.meta["interval_yr"]
    cs, cs_t = CS_SOURCE.get(o.window, (o.coastsat, o.coastsat_target))
    df["coastsat_lrr_m_yr"] = cs["mean_lrr"].to_numpy(float)
    df["coastsat_target_m"] = df["coastsat_lrr_m_yr"] * years
    df["coastsat_target_lowess_m"] = cs_t["target_lrr_m_yr"].to_numpy(float) * years
    df["dune_rate_m_yr"] = o.endpoint["mean_lrr"].to_numpy(float)
    df["dune_measured_change_m"] = df["dune_rate_m_yr"] * o.meta["interval_yr"]
    df["dune_target_m"] = df["dune_rate_m_yr"] * years
    df["dune_target_lowess_m"] = o.endpoint_target["target_lrr_m_yr"].to_numpy(float) * years
    df["target_difference_m"] = df["coastsat_target_m"] - df["dune_target_m"]
    for key, mdf in mdfs.items():
        col = f"model_{MODEL_SETS[key]}_m"
        df[col] = (mdf["change_rate_m_yr"].to_numpy(float) * years
                   if mdf is not None else np.nan)
    return df


def skill_rows(window, df):
    lo, hi = rw.INTERIOR
    sel = df[df["domain_number"].between(lo, hi)]
    rows = []
    for key, folder in MODEL_SETS.items():
        m = sel[f"model_{folder}_m"]
        for target, col in (("coastsat", "coastsat_target_m"),
                            ("coastsat_lowess", "coastsat_target_lowess_m"),
                            ("duneline", "dune_target_m"),
                            ("duneline_lowess", "dune_target_lowess_m")):
            r = (m - sel[col]).dropna()
            ok = m.notna() & sel[col].notna()
            rows.append({"window": "{}_{}".format(*window), "model_ends": folder,
                         "target": target, "n": int(len(r)),
                         "bias_m": round(float(r.mean()), 3),
                         "rmse_m": round(float(np.sqrt((r ** 2).mean())), 3),
                         "r": round(float(np.corrcoef(m[ok], sel[col][ok])[0, 1]), 3)})
    return rows


def draw(ax, o, df, folder, half, label):
    obs.draw_panel(ax, df.assign(mean_lrr=np.nan, std_lrr=0.0), half,
                   label=label, std=False)
    x = df["domain_number"].to_numpy(float)
    cs, du = df["coastsat_target_m"].to_numpy(float), df["dune_target_m"].to_numpy(float)
    shade_beach_width(ax, x, cs, du)
    ax.plot(x, du, color=rw.C_DUNE_TARGET, lw=LW, zorder=11)
    ax.plot(x, cs, color=rw.C_CS_TARGET, lw=LW, zorder=11)
    ax.plot(x, df[f"model_{folder}_m"], color=INK, lw=LW_MODEL, zorder=12)
    mark_offaxis(ax, x, du, half, color=rw.C_DUNE_TARGET)
    mark_offaxis(ax, x, cs, half, color=rw.C_CS_TARGET)
    mark_offaxis(ax, x, df[f"model_{folder}_m"], half, color=INK)
    obs.draw_shoals(ax, label=label)
    fills = obs.fills_in(*o.window)
    if fills:
        obs.draw_fills(ax, fills, half)


def _pad_title(ax, window):
    """Lift the centred title clear of the fill bars.

    draw_fills puts its bars at 1.025 in axes fractions with the year above
    them, so a title at the default pad lands on "2022 fill". _title() has
    already set the bold letter at the left; re-setting only the centred
    string keeps it and moves both (pad is per-axes in matplotlib).
    """
    if obs.fills_in(*window):
        ax.set_title(ax.get_title(loc="center"), loc="center", pad=20)


def figure(observations, frames, key, folder, half, skill_df):
    fig, axes = plt.subplots(len(WINDOWS), 1, sharex=True, sharey=True,
                             constrained_layout=True, figsize=figsize("double", height=5.6))
    for i, (ax, o) in enumerate(zip(axes, observations)):
        draw(ax, o, frames[o.window], folder, half, label=(i == 0))
        ax.yaxis.set_major_locator(MultipleLocator(y_tick(half)))
        _title(ax, i, "{}, {}–{} ({})".format(
            cs_noun(), *o.window, cs_method(o.window)))
        _pad_title(ax, o.window)
    axes[-1].set_xlabel(DOMAIN_AXIS_LABEL)
    fig.supylabel(y_label(), fontsize=9)
    handles = [Line2D([], [], color=rw.C_CS_TARGET, lw=LW, label=cs_label()),
               Line2D([], [], color=rw.C_DUNE_TARGET, lw=LW, label=dune_label()),
               Line2D([], [], color=INK, lw=LW_MODEL,
                      label=f"CASCADE ({MODEL_LABEL[key]})")] + beach_width_handles()
    # Two columns (2026-09-29): at three, the two long target labels shared a
    # row with the model and the beach-width swatches and ran off the page.
    # Filled column-first, so the lines stack left and the swatches right.
    fig.legend(handles=handles, loc="outside lower center", ncol=2, frameon=False)
    sk = skill_df[skill_df["model_ends"] == folder]
    def _s(w, t):
        r = sk[(sk["window"] == w) & (sk["target"] == t)].iloc[0]
        return f"bias {_fmt(r['bias_m'], True)} {_u()}, RMSE {_fmt(r['rmse_m'])} {_u()}"
    stats = " ".join(
        f"{w.replace('_', '–')}: against CoastSat {_s(w, 'coastsat')}; against the "
        f"dune line {_s(w, 'duneline')}." for w in ("1996_2010", "2010_2024"))
    metas = {o.window: o.meta for o in observations}
    dates = "; ".join(f"{m['start_date']} to {m['end_date']} ({m['interval_yr']:.1f} yr) "
                      f"for {w[0]}–{w[1]}" for w, m in metas.items())
    caption(fig, (
        "The two candidate targets and the CASCADE hindcast as "
        + _quantity().replace("the 14-yr", "each 14-yr") + " by GIS domain (1 at Cape "
        "Point, 90 at Pea Island), seaward positive, domain means. Blue: the "
        f"CoastSat target, {cs_clause()}{_x14(', multiplied by')}. Red: the dune-line target, the measured net "
        f"change between the digitized lines ({dates}; the 2023 date assumed) "
        f"divided by its interval{_x14(' and multiplied by')}. Black: the model's own "
        f"{_model_quantity()}, {MODEL_CLAUSE[key]}, "
        "full management, groin off. The space between the two targets is the "
        "beach-width change they imply: solid grey where the beach widened, hatched "
        "where it narrowed. Interior GIS 2–89, model minus target: " + stats
        + " The LOWESS-smoothed scores, as the runs are graded, are in "
        f"tables/skill.csv. One y axis, ±{half:g} {_u()}, the same on every figure here."
        + over_note([(frames[o.window], ["coastsat_target_m", "dune_target_m",
                                         f"model_{folder}_m"], "{}–{} ".format(*o.window))
                     for o in observations], half)))
    # Stem carries the target mode AND the model set (Hannah, 2026-09-21):
    # six folders wrote this same basename, two targets x three model sets.
    out = save(fig, OUT_DIR / folder
               / f"target_comparison_{_stem_mode()}_{folder}_1996_2010_2024")
    plt.close(fig)
    return out


def _dune_title(window):
    return (f"Total dune line change, {window[0]}–{window[1]} (measured, scaled to 14 yr)"
            if _net() else
            f"Dune line change rate, {window[0]}–{window[1]} (measured over the survey interval)")


def _target_legend():
    return ("Target, net change over 14 yr (seaward / landward)" if _net()
            else "Target, change rate (seaward / landward)")


PAIR_KEYS = (("coastsat", "coastsat_target_m", rw.C_CS_TARGET,
              None,
              "CASCADE, ends solved on CoastSat"),
             (rw.MAIN_DUNE, "dune_target_m", rw.C_DUNE_TARGET,
              "Total dune line change (measured, scaled to 14 yr)",
              "CASCADE, ends solved on the dune line"))
LS_MODEL = "-"
# Observed pale and thick, model dark and thin (Hannah, 2026-09-19: the
# same-hue solid/dashed pair was hard to tell apart). The pale pair are the
# RdBu light poles, the dark pair the house blue and red.
PALE = {"coastsat": "#92c5de", rw.MAIN_DUNE: "#f4a582"}
LW_TARGET_PALE = 2.8
LW_MODEL_DARK = 1.1


def end_values(rows_by_key):
    """{(key, window): (gis1, gis90, n_nonzero)} from the run index."""
    idx = rw.load_run_index(rw.RUN_INDEX)
    out = {}
    for key, rows in rows_by_key.items():
        for r in rows:
            kind, tag = rw.legacy_arm_to_kind_tag(r["arm"])
            h = idx[(idx["run_name"] == r["run_name"]) & (idx["kind"] == kind)
                    & (idx["tag"] == tag)].iloc[0]
            w = tuple(int(x) for x in r["window"].split("_"))
            out[(key, w)] = (float(h["be_rate_gis1_m_yr"]), float(h["be_rate_gis90_m_yr"]),
                             int(h["be_nonzero_domains"]))
    return out


def paired_figure(observations, frames, half, skill_df, ends, smoothed=False):
    """Each target with the run solved on it, ONE FIGURE PER WINDOW, one panel
    per target (Hannah, 2026-09-19, style B of three rendered candidates):
    the target as the house fill (blue seaward, red landward), its run as the
    black line, the misfit the gap between them. Both windows on one y axis.

    smoothed=True (2026-09-19, Hannah): the fill is the target AS GRADED (raw
    domain means over GIS 1-10, the rw.TARGET_WINDOW-domain LOWESS beyond, the form the runs
    and the edge solve are scored against), the raw domain means as dots over
    it; written to paired_smoothed/."""
    sfx = "_lowess" if smoothed else ""
    tgt = {"coastsat": "coastsat" + sfx, rw.MAIN_DUNE: "duneline" + sfx}
    pair_keys = {k for k, *_ in PAIR_KEYS}
    for (key, w), (_, _, n) in ends.items():
        if key in pair_keys and n != 2:
            raise SystemExit(f"{key} {w}: {n} nonzero source/sink domains, expected "
                             "the two ends only; the caption would be wrong")
    sk = skill_df.set_index(["window", "model_ends", "target"])
    names = {"coastsat": ("CoastSat target", "coastsat", "ends_solved_on_coastsat"),
             rw.MAIN_DUNE: ("Total dune line change", "duneline", "ends_solved_on_duneline")}
    out = []
    for o in observations:
        w = "{}_{}".format(*o.window)
        df = frames[o.window]
        x = df["domain_number"].to_numpy(float)
        fig, axes = plt.subplots(2, 1, sharex=True, sharey=True, constrained_layout=True,
                                 figsize=figsize("double", height=5.6))
        for i, (ax, (key, col, _, _, _)) in enumerate(zip(axes, PAIR_KEYS)):
            fill_col = col.replace("_m", "_lowess_m") if smoothed else col
            obs.draw_panel(ax, df.assign(mean_lrr=df[fill_col], std_lrr=0.0), half,
                           label=(i == 0), std=False)
            if smoothed:
                raw = df[col].to_numpy(float)
                ax.scatter(x, raw, s=9, lw=0, alpha=0.8, zorder=11,
                           c=np.where(raw < 0, obs.C_ERODE, obs.C_ACCRETE))
            obs.draw_shoals(ax, label=(i == 0))
            ax.plot(x, df[f"model_{MODEL_SETS[key]}_m"], color=INK, lw=LW_MODEL, zorder=12)
            mark_offaxis(ax, x, df[fill_col], half)
            mark_offaxis(ax, x, df[f"model_{MODEL_SETS[key]}_m"], half, color=INK)
            ax.yaxis.set_major_locator(MultipleLocator(y_tick(half)))
            _title(ax, i, (f"{cs_noun()}, {o.window[0]}–{o.window[1]} "
                           f"({cs_method(o.window)}), and its run"
                           if key == "coastsat" else
                           _dune_title(o.window) + ", and its run"))
        fills = obs.fills_in(*o.window)
        if fills:
            obs.draw_fills(axes[0], fills, half)
            _pad_title(axes[0], o.window)
        axes[-1].set_xlabel(DOMAIN_AXIS_LABEL)
        fig.supylabel(y_label(), fontsize=9)
        handles = [(Line2D([], [], color=obs.C_ACCRETE, lw=1.0),
                    Line2D([], [], color=obs.C_ERODE, lw=1.0))]
        labels = [("Target as graded (raw GIS 1–10, LOWESS beyond)"
                   if smoothed else _target_legend())]
        if smoothed:
            handles.append((Line2D([], [], color=obs.C_ACCRETE, marker="o", ms=3, lw=0),
                            Line2D([], [], color=obs.C_ERODE, marker="o", ms=3, lw=0)))
            labels.append("Raw domain means")
        handles.append(Line2D([], [], color=INK, lw=LW_MODEL))
        labels.append("CASCADE, ends solved on that target")
        fig.legend(handles=handles, labels=labels, loc="outside lower center",
                   ncol=len(handles), frameon=False,
                   handler_map={tuple: HandlerTuple(ndivide=None, pad=0.3)})
        (c1, c90, _), (d1, d90, _) = ends[("coastsat", o.window)], ends[(rw.MAIN_DUNE, o.window)]
        fig.suptitle(_targets_line(o) + chr(10)
                     + "{}–{}: source/sink correction at the end domains (GIS 1 and 90) "
                       "only; none on GIS 2–89".format(*o.window) + chr(10)
                     + f"end terms GIS 1 / 90 (m/yr): (a) {c1:+.1f} / {c90:+.1f}   "
                       f"(b) {d1:+.1f} / {d90:+.1f}", fontsize=8.5)
        m = o.meta
        caption(fig, (
            f"{o.window[0]}–{o.window[1]}: each candidate target with the CASCADE run "
            "calibrated to it, as " + _quantity() + " "
            "by GIS domain (1 at Cape Point, 90 at Pea Island), seaward positive; "
            + ("the targets SMOOTHED as the runs are graded: the raw domain means over "
               f"GIS 1–10 and a {rw.TARGET_WINDOW}-domain LOWESS of the transect values beyond, drawn as the "
               "fill, with the raw domain means as dots. " if smoothed else "domain means. ")
            + "(a) The CoastSat target, " + cs_clause() + _x14(" x") + ", "
            "as the fill (blue seaward, red landward), and the edgeBE run whose "
            "two end domains were solved against it (black). (b) The dune-line target, "
            f"the measured net change between the digitized lines ({m['start_date']} to "
            f"{m['end_date']}, {m['interval_yr']:.1f} yr"
            + (", the end date assumed" if m['end_date_assumed'] else "")
            + ") divided by its interval" + _x14(" and x") + ", as the fill, and the edgeBE run whose "
            "end domains were solved against the dune line (black). The model lines are "
            "each run's own " + ("net change" if _net() else "endpoint rate") + ", unchanged; the gap between line and fill is the "
            "misfit. THE ONLY SOURCE/SINK CORRECTION IN EITHER RUN IS THE BOUNDARY TERM AT "
            "GIS 1 AND GIS 90, whose values are in the figure title; every domain from GIS "
            "2 to 89 carries none, so the interior is the model's own response. Full "
            "management, groin off. Interior GIS 2–89, model minus its own "
            + ("smoothed " if smoothed else "") + "target: (a) "
            "{} {u} bias, {} {u} RMSE; (b) {} {u} bias, {} {u} RMSE. The y axis "
            "(±{:g} {u}) is the same on every {e}figure in target_comparison.{} Scores against the other "
            "target, raw and smoothed, are in tables/skill.csv.".format(
                _fmt(sk.loc[(w, "ends_solved_on_coastsat", tgt["coastsat"]), "bias_m"], True),
                _fmt(sk.loc[(w, "ends_solved_on_coastsat", tgt["coastsat"]), "rmse_m"]),
                _fmt(sk.loc[(w, "ends_solved_on_duneline", tgt[rw.MAIN_DUNE]), "bias_m"], True),
                _fmt(sk.loc[(w, "ends_solved_on_duneline", tgt[rw.MAIN_DUNE]), "rmse_m"]), half,
                over_note([(df, ["coastsat_target_m", "model_ends_solved_on_coastsat_m",
                                 "dune_target_m", "model_ends_solved_on_duneline_m"], "")],
                          half), u=_u(), e="" if _net() else "m/yr ")))
        out += save(fig, OUT_DIR / ("paired_smoothed" if smoothed else "paired")
                    / (f"target_and_own_run_{_stem_mode()}_{w}"
                       f"{'_smoothed' if smoothed else ''}"))
        plt.close(fig)
    return out


def unsolved_figure(observations, frames, half, skill_df, ends, smoothed=False):
    """The UNSOLVED run against both targets, one figure per window (Hannah,
    2026-09-21, for her advisor). The paired form, except that there is only
    one run: the same zeroBE line is drawn in both panels, because no part of
    it was fitted to either target. (a) the CoastSat target as the fill, (b)
    the dune-line target, the run in black over each.

    With no end solve the two end domains are the model's own response too,
    so this is the only figure here whose GIS 1 and GIS 90 mean anything."""
    sfx = "_lowess" if smoothed else ""
    col_model = f"model_{MODEL_SETS[UNSOLVED]}_m"
    for w_, (_, _, n) in ((k[1], v) for k, v in ends.items() if k[0] == UNSOLVED):
        if n != 0:
            raise SystemExit(f"unsolved {w_}: {n} nonzero source/sink domains, "
                             "expected none; the caption would be wrong")
    sk = skill_df.set_index(["window", "model_ends", "target"])
    out = []
    for o in observations:
        w = "{}_{}".format(*o.window)
        df = frames[o.window]
        x = df["domain_number"].to_numpy(float)
        fig, axes = plt.subplots(2, 1, sharex=True, sharey=True, constrained_layout=True,
                                 figsize=figsize("double", height=5.6))
        for i, (ax, (key, col, _, _, _)) in enumerate(zip(axes, PAIR_KEYS)):
            fill_col = col.replace("_m", "_lowess_m") if smoothed else col
            obs.draw_panel(ax, df.assign(mean_lrr=df[fill_col], std_lrr=0.0), half,
                           label=(i == 0), std=False)
            if smoothed:
                raw = df[col].to_numpy(float)
                ax.scatter(x, raw, s=9, lw=0, alpha=0.8, zorder=11,
                           c=np.where(raw < 0, obs.C_ERODE, obs.C_ACCRETE))
            obs.draw_shoals(ax, label=(i == 0))
            ax.plot(x, df[col_model], color=INK, lw=LW_MODEL, zorder=12)
            mark_offaxis(ax, x, df[fill_col], half)
            mark_offaxis(ax, x, df[col_model], half, color=INK)
            ax.yaxis.set_major_locator(MultipleLocator(y_tick(half)))
            _title(ax, i, (f"{cs_noun()}, {o.window[0]}–{o.window[1]} "
                           f"({cs_method(o.window)})"
                           if key == "coastsat" else _dune_title(o.window)))
        fills = obs.fills_in(*o.window)
        if fills:
            obs.draw_fills(axes[0], fills, half)
            _pad_title(axes[0], o.window)
        axes[-1].set_xlabel(DOMAIN_AXIS_LABEL)
        fig.supylabel(y_label(), fontsize=9)
        handles = [(Line2D([], [], color=obs.C_ACCRETE, lw=1.0),
                    Line2D([], [], color=obs.C_ERODE, lw=1.0))]
        labels = [("Target as graded (raw GIS 1–10, LOWESS beyond)"
                   if smoothed else _target_legend())]
        if smoothed:
            handles.append((Line2D([], [], color=obs.C_ACCRETE, marker="o", ms=3, lw=0),
                            Line2D([], [], color=obs.C_ERODE, marker="o", ms=3, lw=0)))
            labels.append("Raw domain means")
        handles.append(Line2D([], [], color=INK, lw=LW_MODEL))
        labels.append("CASCADE, no source/sink anywhere (the same run in both panels)")
        fig.legend(handles=handles, labels=labels, loc="outside lower center",
                   ncol=len(handles), frameon=False,
                   handler_map={tuple: HandlerTuple(ndivide=None, pad=0.3)})
        fig.suptitle(_targets_line(o) + chr(10)
                     + "{}–{}: NO source/sink correction in any domain, the two ends "
                     "included".format(*o.window) + chr(10)
                     + "the same run in both panels; neither target was fitted",
                     fontsize=8.5)
        m = o.meta
        caption(fig, (
            f"{o.window[0]}–{o.window[1]}: the UNCALIBRATED CASCADE run against both "
            "candidate targets, as " + _quantity() + " "
            "by GIS domain (1 at Cape Point, 90 at Pea Island), seaward positive; "
            + ("the targets SMOOTHED as the runs are graded: the raw domain means over "
               f"GIS 1–10 and a {rw.TARGET_WINDOW}-domain LOWESS of the transect values beyond, drawn as the "
               "fill, with the raw domain means as dots. " if smoothed else "domain means. ")
            + "(a) The CoastSat target, " + cs_clause() + _x14(" x") + ", as the fill (blue "
            "seaward, red landward). (b) The dune-line target, the measured net change "
            f"between the digitized lines ({m['start_date']} to {m['end_date']}, "
            f"{m['interval_yr']:.1f} yr"
            + (", the end date assumed" if m['end_date_assumed'] else "")
            + ") divided by its interval" + _x14(" and x") + ", as the fill. THE BLACK LINE IS THE "
            "SAME RUN IN BOTH PANELS: the zeroBE arm of the matrix cell, which carries NO "
            "source/sink term in ANY domain, the two ends included, so every one of the 90 "
            "domains is the model's own response and neither target was fitted. Full "
            "management, groin off. The gap between line and fill is the misfit. Interior "
            "GIS 2–89, model minus " + ("smoothed " if smoothed else "") + "target: (a) "
            "{} {u} bias, {} {u} RMSE; (b) {} {u} bias, {} {u} RMSE. The y axis "
            "(±{:g} {u}) is the same on every {e}figure in target_comparison.{} Scores for every "
            "model set against every target, raw and smoothed, are in tables/skill.csv.".format(
                _fmt(sk.loc[(w, MODEL_SETS[UNSOLVED], "coastsat" + sfx), "bias_m"], True),
                _fmt(sk.loc[(w, MODEL_SETS[UNSOLVED], "coastsat" + sfx), "rmse_m"]),
                _fmt(sk.loc[(w, MODEL_SETS[UNSOLVED], "duneline" + sfx), "bias_m"], True),
                _fmt(sk.loc[(w, MODEL_SETS[UNSOLVED], "duneline" + sfx), "rmse_m"]), half,
                over_note([(df, ["coastsat_target_m", "dune_target_m", col_model], "")],
                          half), u=_u(), e="" if _net() else "m/yr ")))
        out += save(fig, OUT_DIR / MODEL_SETS[UNSOLVED]
                    / (f"unsolved_run_and_targets_{_stem_mode()}_{w}"
                       f"{'_smoothed' if smoothed else ''}"))
        plt.close(fig)
    return out


CS_SOURCE = {}   # window -> (domain frame, LOWESS target frame) when full-period


def load_model_sets():
    """key -> ({window: model frame}, index rows) for the three model sets.
    Factored out of main 2026-09-21 so smoothing_scale.py reads exactly
    the same runs; it depends on CS_MODE, which must be set first."""
    models = {}
    for key in MODEL_SETS:
        if key == UNSOLVED:
            loaded = [rw.load_model(w, UNSOLVED_RUNS.get(w), key, UNSOLVED_PRESET)
                      for w in rw.WINDOWS]
            mdfs, rows = [m for m, _ in loaded], [r for _, r in loaded]
        elif key == "coastsat" and CS_MODE == "projected":
            runs = full_period_runs()
            loaded = [rw.load_model(w, runs[w], key) for w in rw.WINDOWS]
            mdfs, rows = [m for m, _ in loaded], [r for _, r in loaded]
        else:
            mdfs, rows = rw.load_models(key)
        models[key] = ({w: m for w, m in zip(rw.WINDOWS, mdfs)},
                       [r for r in rows if tuple(int(x) for x in r["window"].split("_")) in WINDOWS])
    return models


def main() -> int:
    global CS_MODE, OUT_DIR, UNITS
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--coastsat-target", choices=sorted(CS_MODES), default="projected",
                    help="projected (default): the 1996-2024 LRR x 14 yr, carried onto "
                         "windows it was not fitted on. total: each window's own LRR. "
                         "'full' and 'subperiod' are the pre-2026-09-21 aliases.")
    ap.add_argument("--units", choices=("net", "rate"), default="net",
                    help="net (default): metres over the 14-yr window. rate: the same "
                         "in m/yr, written to <mode>/change_rate/.")
    args = ap.parse_args()
    CS_MODE = CS_CANON[args.coastsat_target]
    UNITS = args.units
    OUT_DIR = ROOT_DIR / CS_MODES[CS_MODE]
    if not _net():
        OUT_DIR = OUT_DIR / RATE_SUBDIR
    apply_style()
    observations = [rw.Observation(w) for w in WINDOWS]
    models = load_model_sets()
    if CS_MODE == "projected":
        full = (obs.load_window(*FULL_WINDOW), rw.load_coastsat_target(FULL_WINDOW))
        CS_SOURCE.update({o.window: full for o in observations})
    frames, skill = {}, []
    for o in observations:
        df = window_values(o, {k: m[0][o.window] for k, m in models.items()})
        if not _net():
            df = to_rate(df)
        frames[o.window] = df
        skill += skill_rows(o.window, df)
    skill_df = pd.DataFrame(skill)

    tables = OUT_DIR / "tables"
    tables.mkdir(parents=True, exist_ok=True)
    for w, df in frames.items():
        (df if _net() else _table_names(df)).round(3).to_csv(
            tables / "domain_values_{}_{}.csv".format(*w), index=False)
    (skill_df if _net() else skill_df.rename(
        columns={"bias_m": "bias_m_yr", "rmse_m": "rmse_m_yr"})).to_csv(
        tables / "skill.csv", index=False)
    pd.DataFrame([dict(r, model_ends=MODEL_SETS[k]) for k, (_, rows) in models.items()
                  for r in rows]).to_csv(OUT_DIR / "runs_used.csv", index=False)

    cols = ["coastsat_target_m", "dune_target_m"] + [f"model_{f}_m" for f in MODEL_SETS.values()]
    # ONE fixed y range for every figure in target_comparison, both versions
    # (Hannah, 2026-09-19); anything beyond it is counted in the caption.
    half = Y_HALF_M if _net() else Y_HALF_RATE
    written = []
    for key, folder in MODEL_SETS.items():
        written += figure(observations, frames, key, folder, half, skill_df)
    ends = end_values({k: rows for k, (_, rows) in models.items()})
    written += paired_figure(observations, frames, half, skill_df, ends)
    written += unsolved_figure(observations, frames, half, skill_df, ends)
    if CS_MODE == "projected":   # the smoothed versions, full-period only (Hannah)
        written += paired_figure(observations, frames, half, skill_df, ends, smoothed=True)
        written += unsolved_figure(observations, frames, half, skill_df, ends, smoothed=True)

    print(skill_df[skill_df["target"].isin(["coastsat", "duneline"])].to_string(index=False))
    print(f"\ny axis +/-{half:g} {_u()}")
    for p in written:
        print("wrote   ", Path(p).relative_to(_REPO))
    return 0


if __name__ == "__main__":
    sys.exit(main())
